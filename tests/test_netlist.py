"""SpectreGenerator output and the netlist-format registry."""

import pytest

from mosplot.optimizer import (
    ISource,
    Instance,
    NetlistGenerator,
    Passive,
    SpectreGenerator,
    VSource,
    available_netlist_formats,
    get_netlist_generator,
    register_netlist_generator,
)


DEVICE_MAP = {"nmos": "nch_lvt", "pmos": "pch_lvt"}


def generate(tmp_path, *, mosfets=(), passives=(), vsources=(), isources=(),
             dimensions=None, passive_params=None, vsource_params=None,
             isource_params=None, context=None, extra_lines=None,
             ports=("VIN", "VOUT", "vdd", "vss")):
    gen = SpectreGenerator("amp", list(ports), ground="vss", context=context)
    path = gen.generate(
        mosfets=list(mosfets),
        passives=list(passives),
        vsources=list(vsources),
        isources=list(isources),
        device_map=DEVICE_MAP,
        dimensions=dimensions or {},
        passive_params=passive_params or {},
        vsource_params=vsource_params or {},
        isource_params=isource_params or {},
        output_path=tmp_path / "design.scs",
        extra_lines=extra_lines,
    )
    return path.read_text().splitlines()


class TestSpectreGenerator:
    def test_subckt_header_and_footer(self, tmp_path):
        lines = generate(tmp_path)
        assert lines[0] == "subckt amp (VIN VOUT vdd vss)"
        assert lines[-1] == "ends amp"

    def test_mosfet_line(self, tmp_path):
        lines = generate(
            tmp_path,
            mosfets=[Instance("M1", "nmos", d="VOUT", g="VIN", s="gnd", b="gnd")],
            dimensions={"M1": {"Width": 4.2e-6, "Length": 200e-9}},
        )
        assert "M1 (VOUT VIN vss vss) nch_lvt l=200n w=4.2u" in lines

    def test_passives_and_external_skipped(self, tmp_path):
        lines = generate(
            tmp_path,
            passives=[
                Passive("CC", "cap", a="x", b="VOUT"),
                Passive("Rz", "res", a="x", b="y"),
                Passive("CL", "cap", a="VOUT", b="gnd", external=True),
            ],
            passive_params={"CC": 1.5e-12, "Rz": 2000.0},
        )
        assert "CC (x VOUT) capacitor c=1.5p" in lines
        assert "Rz (x y) resistor r=2k" in lines
        assert not any("CL" in ln for ln in lines)

    def test_plain_vsource(self, tmp_path):
        lines = generate(
            tmp_path,
            vsources=[
                VSource("VDD", p="vdd", n="vss", supply=True),
                VSource("VB", p="vb", n="vss"),
            ],
            vsource_params={"VB": 0.55},
        )
        assert "VB (vb vss) vsource dc=550m" in lines
        assert not any(ln.startswith("VDD ") for ln in lines)  # supply not emitted

    def test_emit_false_suppresses_source(self, tmp_path):
        lines = generate(
            tmp_path,
            vsources=[VSource("VB", p="vb", n="vss", emit=False)],
            vsource_params={},
        )
        assert not any("VB" in ln for ln in lines)

    def test_mirror_bias_synthesis_nmos_master(self, tmp_path):
        # NMOS master (source at vss): diode replica + reference current
        # pulled from vdd into the bias node.
        lines = generate(
            tmp_path,
            mosfets=[Instance("M3", "nmos", d="x", g="vbn", s="gnd", b="gnd")],
            dimensions={"M3": {"Width": 2e-6, "Length": 400e-9}},
            vsources=[
                VSource("VDD", p="vdd", n="vss", supply=True),
                VSource("VBN", p="vbn", n="vss", mirror="M3"),
            ],
            isource_params={"VBN": 10e-6},
        )
        assert "Mvbn (vbn vbn vss vss) nch_lvt l=400n w=2u" in lines
        assert "Ivbn (vdd vbn) isource dc=10u" in lines
        assert not any("vsource" in ln for ln in lines)

    def test_mirror_bias_synthesis_pmos_master(self, tmp_path):
        # PMOS master (source at vdd): the reference current is pulled from
        # the bias node down to ground instead.
        lines = generate(
            tmp_path,
            mosfets=[Instance("M3", "pmos", d="x", g="vbp", s="vdd", b="vdd")],
            dimensions={"M3": {"Width": 2e-6, "Length": 400e-9}},
            vsources=[
                VSource("VDD", p="vdd", n="vss", supply=True),
                VSource("VBP", p="vbp", n="vss", mirror="M3"),
            ],
            isource_params={"VBP": 10e-6},
        )
        assert "Mvbp (vbp vbp vdd vdd) pch_lvt l=400n w=2u" in lines
        assert "Ivbp (vbp vss) isource dc=10u" in lines

    def test_isource_emitted(self, tmp_path):
        lines = generate(
            tmp_path,
            isources=[ISource("IREF", p="ib", n="vss")],
            isource_params={"IREF": 5e-6},
        )
        assert "IREF (ib vss) isource dc=5u" in lines

    def test_context_parameters_line(self, tmp_path):
        lines = generate(tmp_path, context={"vcm": 0.6, "vdd": 1.2})
        assert lines[0] == "parameters vcm=600m vdd=1.2"

    def test_extra_lines_before_ends(self, tmp_path):
        lines = generate(tmp_path, extra_lines=["bcmfb (vcmfb 0) bsource v=0.6"])
        assert lines[-2] == "bcmfb (vcmfb 0) bsource v=0.6"


class TestRegistry:
    def test_spectre_registered(self):
        assert "spectre" in available_netlist_formats()
        assert get_netlist_generator("spectre") is SpectreGenerator

    def test_name_normalisation(self):
        assert get_netlist_generator("  SPECTRE ") is SpectreGenerator

    def test_unknown_format(self):
        with pytest.raises(ValueError, match="Unsupported netlist format"):
            get_netlist_generator("hspice2000")

    def test_register_requires_subclass(self):
        with pytest.raises(TypeError):
            register_netlist_generator("bad", dict)

    def test_register_custom_generator(self):
        class Custom(SpectreGenerator):
            pass

        register_netlist_generator("custom-test", Custom)
        assert get_netlist_generator("custom-test") is Custom

    def test_empty_name_rejected(self):
        with pytest.raises(ValueError):
            get_netlist_generator("  ")
