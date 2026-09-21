"""Spawn workers once per optimization; share arrays, send only candidate vectors."""

from concurrent.futures import ProcessPoolExecutor
from contextlib import contextmanager, ExitStack
from copy import copy
import multiprocessing as mp
from multiprocessing.shared_memory import SharedMemory
from multiprocessing.util import Finalize
import pickle
import shutil
import sys

import numpy as np

from mosplot.expressions import Expression, build_expressions
from mosplot.interpolation import GmIdTable
from .device import DeviceTable
from .fast_mosfet import FastMosfet


_optimizer = None
_handles = []


def _close_worker():
    global _optimizer
    _optimizer = None
    for handle in _handles:
        handle.close()
    _handles.clear()


def _attach(descriptor):
    name, shape, dtype = descriptor
    handle = SharedMemory(name=name)
    _handles.append(handle)
    array = np.ndarray(shape, dtype=dtype, buffer=handle.buf)
    array.flags.writeable = False
    return array


def _init_worker(optimizer, devices):
    global _optimizer
    Finalize(None, _close_worker, exitpriority=10)
    tables = {}
    for key, (device_state, fm_state, table_state, stack, needed) in devices.items():
        table = GmIdTable.__new__(GmIdTable)
        table.__dict__.update(table_state)
        table._data_stack = _attach(stack)
        table._data = {p: table._data_stack[i] for i, p in enumerate(table.params)}
        fm = FastMosfet.__new__(FastMosfet)
        fm.__dict__.update(fm_state)
        fm._table = table
        vdsat = "vdssat" if "vdssat" in fm._param_names else "vdsat"
        fm.__dict__.update(build_expressions(fm._width, vdsat))
        device = DeviceTable.__new__(DeviceTable)
        device.__dict__.update(device_state)
        device.fm, device._lut = fm, table
        device._needed_data = _attach(needed)
        tables[key] = device
    for corner in optimizer.corners:
        for polarity in ("nmos", "pmos"):
            key = getattr(corner, polarity)
            setattr(corner, polarity, tables[key] if key is not None else None)
    _optimizer = optimizer


def objective(x):
    return _optimizer._objective(x)


@contextmanager
def objective_pool(optimizer, workers):
    if workers == 1:
        yield None
        return

    # Preserve subclass methods and instance state without pickling the loaded tables
    # or a parent-owned executor. Do not rerun the optimizer's constructor in workers.
    worker_optimizer = copy(optimizer)
    worker_optimizer.corners = []
    worker_optimizer.executor = None
    try:
        pickle.dumps(worker_optimizer)
    except (pickle.PicklingError, AttributeError, TypeError) as exc:
        raise ValueError("Parallel optimization requires importable circuit and optimizer "
                         "classes with picklable instance state; otherwise use workers=1.") from exc

    with ExitStack() as cleanup:
        arrays = {}

        def share(array):
            key = id(array)
            if key not in arrays:
                # POSIX allocation can succeed then SIGBUS on the first write if /dev/shm
                # is full (common in containers). Fail before touching an oversized block.
                if sys.platform == "linux" and array.nbytes > shutil.disk_usage("/dev/shm").free:
                    raise MemoryError("Not enough /dev/shm space for lookup tables; "
                                      "increase shared memory capacity or use workers=1.")
                handle = SharedMemory(create=True, size=array.nbytes)
                cleanup.callback(handle.unlink)
                cleanup.callback(handle.close)
                view = np.ndarray(array.shape, dtype=array.dtype, buffer=handle.buf)
                view[:] = array
                arrays[key] = (handle.name, array.shape, array.dtype)
            return arrays[key]

        devices, corners = {}, []
        polarities = {m.kind for m in optimizer.model.MOSFETS}
        # Include auxiliary lookups already used by the model outside its netlist topology.
        polarities.update(p for c in optimizer.corners for p in ("nmos", "pmos")
                          if p in vars(c))
        for original in optimizer.corners:
            corner = copy(original)
            for polarity in ("nmos", "pmos"):
                if polarity not in polarities:
                    setattr(corner, polarity, None)
                    continue
                device = getattr(original, polarity)
                key = (id(device.fm), device.pol)
                if key not in devices:
                    table = device._lut
                    devices[key] = (
                        {k: v for k, v in vars(device).items()
                         if k not in ("fm", "_lut", "_needed_data")},
                        {k: v for k, v in vars(device.fm).items()
                         if k != "_table" and not isinstance(v, Expression)},
                        {k: v for k, v in vars(table).items()
                         if k not in ("_data", "_data_stack")},
                        share(table._data_stack), share(device._needed_data),
                    )
                setattr(corner, polarity, key)
            corners.append(corner)
        worker_optimizer.corners = corners
        # ExitStack shuts down workers before releasing their shared arrays, even on error.
        pool = cleanup.enter_context(ProcessPoolExecutor(
            max_workers=workers, mp_context=mp.get_context("spawn"),
            initializer=_init_worker,
            initargs=(worker_optimizer, devices),
        ))
        yield pool
