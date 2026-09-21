from typing import List

import numpy as np


def list_to_string(lst: List) -> str:
    return "\n".join(s for s in lst if s is not None)


def sweep_axis(sweep_range: tuple) -> np.ndarray:
    """Axis values for a (start, stop, step) sweep tuple, endpoint included."""
    start, stop, step = sweep_range
    n = int(round((stop - start) / step)) + 1
    return np.linspace(start, stop, n)
