import numpy as np
from typing import Any, Tuple, Optional, Union, Dict
import numpy.typing as npt
from .expressions import Expression


def load_lookup_table(path: str) -> dict:
    return np.load(path, allow_pickle=True).tolist()


def tile_arrays(A: npt.NDArray, B: npt.NDArray) -> Tuple[npt.NDArray, npt.NDArray]:
    if A.ndim == 1 and B.ndim == 2:
        if A.shape[0] == B.shape[0]:
            return np.tile(A, (B.shape[1], 1)).T, B
        elif A.shape[0] == B.shape[1]:
            return np.tile(A, (B.shape[0], 1)), B
    elif B.ndim == 1 and A.ndim == 2:
        if B.shape[0] == A.shape[0]:
            return A, np.tile(B, (A.shape[1], 1)).T
        elif B.shape[0] == A.shape[1]:
            return A, np.tile(B, (A.shape[0], 1))
    return A, B


def extract_2d_table(
    lookup_table: dict,
    width: Any,
    lengths: Optional[Union[float, list, npt.NDArray]] = None,
    vsb: Optional[Union[float, Tuple[float, float, float]]] = None,
    vgs: Optional[Union[float, Tuple[float, float, float]]] = None,
    vds: Optional[Union[float, Tuple[float, float, float]]] = None,
    primary: Optional[str] = None,
    parameters: Optional[list] = None,
) -> Tuple[Optional[int], list, dict]:
    # At least two sweep parameters are required.
    params = [lengths is not None, vsb is not None, vgs is not None, vds is not None]
    if sum(params) < 2:
        raise ValueError("Provide at least two parameters.")
    variables = {"lengths": False, "vsb": False, "vgs": False, "vds": False}
    if primary:
        variables[primary] = True

    def get_indices(var: str, target: Any) -> Tuple[Any, Any]:
        data = lookup_table[var]
        if isinstance(target, tuple):
            start, end = target[:2]
            indices = np.where((data >= start) & (data <= end))[0]
            if len(target) == 3:
                step = int(target[2] / (data[1] - data[0]))
                indices = indices[::step]
        elif var == "lengths" and isinstance(target, (list, np.ndarray)):
            mask = np.isin(lookup_table["lengths"], np.array(target))
            indices = np.nonzero(mask)[0]
            indices = np.array(indices, dtype=int)
        else:
            variables[var] = True
            index = (np.abs(data - target)).argmin()
            return np.array([index]), data[index]
        return indices, data[indices]

    indices_and_values = {
        "lengths": get_indices("lengths", lengths)
        if lengths is not None
        else (slice(None), lookup_table["lengths"]),
        "vsb": get_indices("vsb", vsb)
        if vsb is not None
        else (slice(None), lookup_table["vsb"]),
        "vgs": get_indices("vgs", vgs)
        if vgs is not None
        else (slice(None), lookup_table["vgs"]),
        "vds": get_indices("vds", vds)
        if vds is not None
        else (slice(None), lookup_table["vds"]),
    }
    slice_indices = []
    filter_values = []
    for key in ["lengths", "vsb", "vgs", "vds"]:
        slice_indices.append(indices_and_values[key][0])
        filter_values.append(indices_and_values[key][1])
    slice_indices = tuple(slice_indices)

    def slice_me(a: npt.NDArray, slices: Tuple) -> npt.NDArray:
        x = a[slices[0], :, :, :]
        x = x[:, slices[1], :, :]
        x = x[:, :, slices[2], :]
        x = x[:, :, :, slices[3]]
        return x

    extracted_table: Dict[str, Any] = {}
    if parameters is None:
        parameters = lookup_table["parameter_names"]
    if parameters:
        for p in parameters:
            if p in lookup_table:
                x = np.squeeze(slice_me(lookup_table[p], slice_indices))
                extracted_table[p] = x.T if (x.ndim > 1 and x.shape[0] > x.shape[1]) else x
    one_key = next(iter(extracted_table))
    extracted_table["width"] = np.array(width)
    extracted_table["lengths"], _ = tile_arrays(
        filter_values[0], extracted_table[one_key]
    )
    extracted_table["vsb"], _ = tile_arrays(filter_values[1], extracted_table[one_key])
    extracted_table["vgs"], _ = tile_arrays(filter_values[2], extracted_table[one_key])
    extracted_table["vds"], _ = tile_arrays(filter_values[3], extracted_table[one_key])
    secondary_idx = None
    if primary:
        secondary_idx = list(variables.values()).index(False)
    return secondary_idx, filter_values, extracted_table


def evaluate_expression(expression: Expression, table: dict, filter_by_rows: Optional[np.ndarray] = None) -> Tuple[np.ndarray, str]:
    filter_by_rows = np.array([]) if filter_by_rows is None else filter_by_rows
    var_list = []
    for var in expression.variables:
        data = table[var]
        if filter_by_rows.size > 0 and data.ndim > 1:
            var_list.append(np.take(data, filter_by_rows, axis=0))
        else:
            var_list.append(data)
    result = (
        expression.function(*var_list)
        if expression.function is not None
        else var_list[0]
    )
    return result, expression.label
