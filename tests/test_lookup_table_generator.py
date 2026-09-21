"""Simulator-independent pieces of the lookup-table generator."""

import warnings

import numpy as np
import pytest

from mosplot.lookup_table_generator.simulators.spice_simulators.base_simulator import (
    BaseSimulator,
)
from mosplot.lookup_table_generator.simulators.spice_simulators.hspice_simulator import (
    HspiceSimulator,
)
from mosplot.lookup_table_generator.simulators.spice_simulators.ngspice_simulator import (
    NgspiceSimulator,
)
from mosplot.lookup_table_generator.simulators.spice_simulators.parsers.ngspice import parse_file as parse_ngspice_raw
from mosplot.lookup_table_generator.simulators.spice_simulators.spectre_simulator import (
    SpectreSimulator,
)
from mosplot.lookup_table_generator.simulators.spice_simulators.spice_mosfet_netlist_generator import (
    SpiceMosfetNetlistGenerator,
)
from mosplot.lookup_table_generator.simulators.spice_simulators.utils import (
    list_to_string,
    sweep_axis,
)
from mosplot.lookup_table_generator.table_cleanup import LookupTableCleaner
from mosplot.lookup_table_generator.transistor_sweep import TransistorSweep


class TestUtils:
    def test_list_to_string_skips_none(self):
        assert list_to_string(["a", None, "b"]) == "a\nb"

    def test_sweep_axis_endpoint_included(self):
        axis = sweep_axis((0.0, 1.2, 0.01))
        assert len(axis) == 121
        assert axis[0] == 0.0 and axis[-1] == pytest.approx(1.2)

    def test_sweep_axis_negative_direction(self):
        axis = sweep_axis((0.0, -1.2, -0.01))
        assert len(axis) == 121
        assert axis[-1] == pytest.approx(-1.2)


class TestLookupTableCleaner:
    def _table(self):
        return {
            "nch": {
                "id": np.array([[1.0, 2.0], [3.0, 4.0]]),
                "weff": np.full((2, 2), 5e-6),  # constant nonzero
                "gmbs": np.zeros((2, 2)),  # all zero -> removed
            }
        }

    def test_cleanup_routes_each_kind(self):
        table = self._table()
        params = ["id", "weff", "gmbs"]
        cleaner = LookupTableCleaner(table, params)
        cleaner.clean_lookup_table()

        assert "id" in table["nch"]  # varying data kept
        assert "weff" not in table["nch"]  # constant -> scalar param
        assert cleaner.scalar_params["nch"]["weff"] == pytest.approx(5e-6)
        assert "gmbs" not in table["nch"]  # all-zero removed
        assert params == ["id"]  # parameters_to_save updated in place


class TestSpiceNetlistGenerator:
    def test_generate_netlist(self):
        sweeps = {"nch": TransistorSweep("nmos", (0, 1.2, 0.01), (0, 1.2, 0.06), (0, -0.6, -0.3), [1e-6])}
        gen = SpiceMosfetNetlistGenerator(
            sweeps, {"w": 10e-6}, ("m1", "m1"), ["/inc/models.inc"], None, None
        )
        netlist = gen.generate_netlist("nch", 1e-6, -0.3)
        text = list_to_string(netlist)
        assert ".include '/inc/models.inc'" in text
        assert "VBS NB 0 DC=-0.3" in text
        assert "m1 ND NG 0 NB nch L=1e-06 w=1e-05" in text


NMOS_SWEEP = TransistorSweep("nmos", (0, 1.2, 0.01), (0, 1.2, 0.06), (0, -0.6, -0.3), [1e-6])
PMOS_SWEEP = TransistorSweep("pmos", (0, -1.2, -0.01), (0, -1.2, -0.06), (0, 0.6, 0.3), [1e-6])


class TestParametersToSaveContract:
    """A user must be able to drop any parameter their PDK lacks (e.g. vdssat)
    via parameters_to_save -- without forking the simulator class."""

    def _ngspice(self, params=None):
        # simulator_path="true": validate_paths only checks the binary exists.
        return NgspiceSimulator(simulator_path="true", parameters_to_save=params)

    def test_ngspice_excluded_parameter_leaves_no_netlist_trace(self):
        sim = self._ngspice(["id", "vth", "vdsat", "gm", "gds"])
        text = list_to_string(sim.setup_dc_simulation(NMOS_SWEEP))
        assert "vdssat" not in text
        assert "weff" not in text

    def test_ngspice_included_parameter_gets_save_and_let(self):
        sim = self._ngspice(["id", "vdssat"])
        text = list_to_string(sim.setup_dc_simulation(NMOS_SWEEP))
        assert "save @m1[vdssat]" in text
        assert "let m_vdssat" in text

    def test_ngspice_id_only_still_defines_its_vector(self):
        sim = self._ngspice(["id"])
        text = list_to_string(sim.setup_dc_simulation(NMOS_SWEEP))
        assert "let i_vds" in text
        assert "let m_vth" not in text

    def test_ngspice_pmos_sign_applied(self):
        sim = self._ngspice(["id", "vth"])
        text = list_to_string(sim.setup_dc_simulation(PMOS_SWEEP))
        assert "let m_vth = -abs(@m1[vth])" in text

    def test_unsupported_parameter_warns_not_crashes(self):
        sim = self._ngspice(["id", "not_a_param"])
        with pytest.warns(UserWarning, match="not_a_param"):
            text = list_to_string(sim.setup_dc_simulation(NMOS_SWEEP))
        assert "not_a_param" not in text

    def test_default_parameters_are_self_consistent(self):
        # Every default parameters_to_save entry must exist in the backend's
        # parameter table, so defaults never warn or silently drop columns.
        for sim in (
            self._ngspice(),
            HspiceSimulator(simulator_path="true"),
            SpectreSimulator(simulator_path="true"),
        ):
            with warnings.catch_warnings():
                warnings.simplefilter("error")
                sim.setup_dc_simulation(NMOS_SWEEP)
            assert set(sim.parameter_table) == set(sim.parameters_to_save)

    def test_spectre_falls_back_to_psp_names(self):
        # PSP models (IHP SG13) call them gmb and vdss; the BSIM name wins when both exist.
        sim = SpectreSimulator(simulator_path="true", parameters_to_save=["gm", "gmbs", "vdsat"])
        sim.setup_dc_simulation(NMOS_SWEEP)
        psp = [{"M1:gm": np.ones(6), "M1:gmb": np.arange(6.0), "M1:vdss": np.full(6, 0.1)}]
        results = sim.extract_parameters(psp, n_vgs=2, n_vds=3)
        assert np.array_equal(results["gmbs"], np.arange(6.0).reshape(2, 3))
        assert np.allclose(results["vdsat"], 0.1)
        bsim = [{"M1:gm": np.ones(6), "M1:gmbs": np.full(6, 2.0), "M1:gmb": np.zeros(6)}]
        assert np.all(sim.extract_parameters(bsim, n_vgs=2, n_vds=3)["gmbs"] == 2.0)

    def test_ngspice_falls_back_to_psp_names(self):
        sim = self._ngspice(["gm", "gmbs", "vdsat"])
        text = list_to_string(sim.setup_dc_simulation(NMOS_SWEEP))
        assert "let m_vdss = abs(@m1[vdss])" in text
        cols = [("@m1[gm]", float), ("@m1[gmb]", float), ("v(m_vdss)", float)]
        psp = np.array([(1.0, g, 0.1) for g in np.arange(6.0)], dtype=cols)
        results = sim.extract_parameters([psp], n_vgs=2, n_vds=3)
        assert np.array_equal(results["gmbs"], np.arange(6.0).reshape(2, 3))
        assert np.allclose(results["vdsat"], 0.1)

    def test_ngspice_raw_drops_vectors_the_model_lacks(self, tmp_path):
        # ngspice writes a parameter the model does not have as zeros marked dims=0
        header = ("Title: t\nDate: d\nPlotname: DC\nFlags: real\nNo. Variables: 3\nNo. Points: 2\n"
                  "Variables:\n\t0\tv(v-sweep)\tvoltage\n\t1\t@m1[gmb]\tadmittance dims=0\n"
                  "\t2\t@m1[gmbs]\tadmittance\nBinary:\n")
        raw = tmp_path / "out.raw"
        raw.write_bytes(header.encode() + np.array([[0.1, 0.0, 5e-5], [0.2, 0.0, 6e-5]]).tobytes())
        (arr,), _ = parse_ngspice_raw(str(raw))
        assert arr.dtype.names == ("v(v-sweep)", "@m1[gmbs]")
        assert np.array_equal(arr["@m1[gmbs]"], [5e-5, 6e-5])

    def test_hspice_extract_skips_unsupported_request(self):
        # "weff" is an ngspice-only parameter; requesting it from hspice must
        # warn during setup and not KeyError during extraction.
        sim = HspiceSimulator(simulator_path="true", parameters_to_save=["id", "weff"])
        with pytest.warns(UserWarning, match="weff"):
            sim.setup_dc_simulation(NMOS_SWEEP)
        analysis = {"m_id": np.zeros((3, 2))}
        results = sim.extract_parameters(analysis, n_vgs=2, n_vds=3)
        assert "id" in results and "weff" not in results


class DummySimulator(BaseSimulator):
    """Concrete stand-in: no binary required, temp files under pytest tmp."""

    def __init__(self, tmp_root=None, **overrides):
        cfg = dict(
            device_parameters={"w": 1e-6},
            simulator_path="true",  # always on PATH
            include_paths=None,
            lib_mappings=None,
            raw_spice=None,
            mos_spice_symbols=("m1", "m1"),
            parameters_to_save=["id"],
            temperature=27,
        )
        cfg.update(overrides)
        super().__init__(**cfg)
        self._tmp_root = tmp_root

    def make_temp_files(self):
        import tempfile

        self.tmp_dir = tempfile.mkdtemp(dir=self._tmp_root)

    def build_simulation_command(self, verbose):
        return ["true"]

    def setup_op_simulation(self, sweep):
        return []

    def setup_dc_simulation(self, sweep):
        return []

    def parse_output(self):
        return None

    def extract_parameters(self, analysis, n_vgs, n_vds):
        return {}


class TestBaseSimulator:
    def test_remove_before_make_is_a_noop(self, tmp_path):
        sim = DummySimulator(tmp_root=tmp_path)
        assert sim.tmp_dir is None
        sim.remove_temp_files()  # must not raise or delete anything

    def test_make_then_remove(self, tmp_path):
        import os

        sim = DummySimulator(tmp_root=tmp_path)
        sim.make_temp_files()
        created = sim.tmp_dir
        assert os.path.isdir(created)
        sim.remove_temp_files()
        assert not os.path.isdir(created)
        assert sim.tmp_dir is None

    def test_missing_include_path_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            DummySimulator(tmp_root=tmp_path, include_paths=["/nonexistent/file.inc"])

    def test_missing_binary_raises(self, tmp_path):
        with pytest.raises(ValueError, match="not accessible"):
            DummySimulator(tmp_root=tmp_path, simulator_path="no-such-binary-xyz")

    def test_clone_reproduces_config(self, tmp_path):
        sim = DummySimulator(tmp_root=tmp_path)
        clone = sim.clone()
        assert type(clone) is DummySimulator
        assert clone.parameters_to_save == sim.parameters_to_save
        assert clone is not sim
