"""Real spawn workers: numerical equivalence, shared tables, and cleanup."""

from multiprocessing.shared_memory import SharedMemory
import os

import numpy as np
import pytest

from mosplot.optimizer import Corner, Optimizer
from mosplot.optimizer.parallel import objective, objective_pool
from test_end_to_end import Mirror, PARAMETERS, TARGET_SPECS, COND


class CustomCost(Optimizer):
    def __init__(self, *args, offset, **kwargs):
        super().__init__(*args, **kwargs)
        self.offset = offset

    def compute_cost(self, specs):
        return super().compute_cost(specs) + self.offset


def worker_state():
    from mosplot.optimizer import parallel
    from mosplot.optimizer.corner import _shared_lut, _shared_fast_mosfet

    a, b = parallel._optimizer.corners
    assert a.nmos is b.nmos and a.pmos is b.pmos
    assert _shared_lut.cache_info().currsize == 0
    assert _shared_fast_mosfet.cache_info().currsize == 0
    for device in (a.nmos, a.pmos):
        assert not device._needed_data.flags.writeable
        assert not device._lut._data_stack.flags.writeable
        assert all(np.shares_memory(v, device._lut._data_stack)
                   for v in device._lut._data.values())
        value = device.fm.interpolate(400e-9, 10., .6, 0., device.fm.id_expression)
        assert np.isfinite(value)
    return os.getpid(), [h.name for h in parallel._handles]


@pytest.fixture
def opt(lut_npz_path, tmp_path):
    corners = [Corner(name, str(lut_npz_path), "nch", "pch",
                      conditions=dict(COND, vout_dc=voltage), cache_dir=str(tmp_path))
               for name, voltage in (("tt", .6), ("tt_low", .5))]
    return Optimizer(Mirror(), PARAMETERS, TARGET_SPECS, corners)


def test_spawn_shares_tables_and_cleans_up_on_error(opt):
    # Exercise both polarities, including an auxiliary lookup outside Mirror's topology.
    opt.corners[0].pmos
    candidates = [[8., 200e-9], [12., 400e-9], [18., 700e-9]]
    expected = [opt._objective(x) for x in candidates]
    with pytest.raises(RuntimeError, match="test cleanup"):
        with objective_pool(opt, 2) as pool:
            np.testing.assert_array_equal(list(pool.map(objective, candidates)), expected)
            pid, names = pool.submit(worker_state).result()
            assert pid != os.getpid()
            assert len(names) == 4  # two arrays per device, reused across voltage corners
            raise RuntimeError("test cleanup")
    for name in names:
        with pytest.raises(FileNotFoundError):
            handle = SharedMemory(name=name)
            handle.close()
    # Shared-memory teardown must not invalidate the parent's device tables.
    np.testing.assert_array_equal([opt._objective(x) for x in candidates], expected)


@pytest.mark.parametrize("maxiter", [0, 3])
def test_parallel_optimization_matches_serial_across_restarts(opt, maxiter):
    serial = opt.optimize(maxiter=maxiter, n_restarts=2, seed=7)
    parallel = opt.optimize(maxiter=maxiter, n_restarts=2, seed=7, workers=2)
    np.testing.assert_array_equal(parallel.x, serial.x)
    assert parallel.fun == serial.fun
    assert parallel.nit == serial.nit
    assert all(r.max_residual < 1e-3 for r in opt.corner_results)


@pytest.mark.parametrize("workers", [0, -1, 1.5, True, "4"])
def test_invalid_workers(opt, workers):
    with pytest.raises(ValueError, match="positive integer"):
        opt.optimize(workers=workers)


def test_custom_cost_and_nmos_only_table(lookup_table, tmp_path):
    path = tmp_path / "nmos_only.npz"
    np.savez(path, lookup_table={"nch": lookup_table["nch"]})
    corners = [Corner("tt", str(path), "nch", "unused", conditions=COND, cache_dir=None)]
    opt = CustomCost(Mirror(), PARAMETERS, TARGET_SPECS, corners, offset=123.)
    candidates = [[8., 200e-9], [12., 400e-9], [18., 700e-9]]
    expected = [opt._objective(x) for x in candidates]
    assert all(123. <= cost < 1e6 for cost in expected)
    with objective_pool(opt, 2) as pool:
        np.testing.assert_array_equal(list(pool.map(objective, candidates)), expected)
    assert "pmos" not in vars(corners[0])  # unused table was never loaded in the parent
    serial = opt.optimize(maxiter=2, seed=7)
    parallel = opt.optimize(maxiter=2, seed=7, workers=2)
    np.testing.assert_array_equal(parallel.x, serial.x)
    assert parallel.fun == serial.fun
