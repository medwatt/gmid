# imports <<<
from typing import Any, Dict, Optional, Tuple, Union
import numpy as np
import numpy.typing as npt
# >>>


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
    width: float,
    length: Optional[Union[float, list, npt.NDArray]] = None,
    vsb: Optional[Union[float, Tuple[float, float], Tuple[float, float, float]]] = None,
    vgs: Optional[Union[float, Tuple[float, float], Tuple[float, float, float]]] = None,
    vds: Optional[Union[float, Tuple[float, float], Tuple[float, float, float]]] = None,
    primary: Optional[str] = None,
    parameters: Optional[list] = None,
) -> Tuple[Optional[str], Dict, dict]:
    # Ensure at least two sweep parameters are provided.
    if sum(param is not None for param in [length, vsb, vgs, vds]) < 2:
        raise ValueError("Provide at least two parameters.")

    # Process a sweep target and return (indices, values, is_scalar).
    def process_target(
        data: npt.NDArray, target: any, var_name: str
    ) -> Tuple[Union[np.ndarray, slice], Any, bool]:
        if isinstance(target, tuple):
            start, end = target[0], target[1]
            inds = np.where((data >= start) & (data <= end))[0]
            if len(target) == 3:
                step = int(target[2] / (data[1] - data[0]))
                inds = inds[::step]
            return inds, data[inds], False
        if var_name == "length" and isinstance(target, (list, np.ndarray)):
            target_arr = np.array(target)
            mask = np.isin(data, target_arr)
            inds = np.nonzero(mask)[0]
            return inds, data[inds], False
        idx = int(np.abs(data - target).argmin())
        return np.array([idx]), data[idx], True

    # Slice a 4D array using indices for each axis.
    def apply_slice(a: npt.NDArray, slices: Tuple) -> npt.NDArray:
        s = a[slices[0], :, :, :]
        s = s[:, slices[1], :, :]
        s = s[:, :, slices[2], :]
        s = s[:, :, :, slices[3]]
        return s

    sweep_params = {"length": length, "vsb": vsb, "vgs": vgs, "vds": vds}
    indices_dict = {}
    values_dict = {}
    is_scalar = {}

    for var, target in sweep_params.items():
        if target is None:
            indices_dict[var] = slice(None)
            values_dict[var] = lookup_table[var]
            is_scalar[var] = False
        else:
            inds, vals, scalar_flag = process_target(lookup_table[var], target, var)
            indices_dict[var] = inds
            values_dict[var] = vals
            is_scalar[var] = scalar_flag

    # Determine the secondary sweep variable.
    secondary = None
    if primary:
        if sweep_params.get(primary) is None:
            raise ValueError("Primary variable must be provided.")
        non_primary_ranges = [var for var in sweep_params if var != primary and not is_scalar[var]]
        if len(non_primary_ranges) == 1:
            secondary = non_primary_ranges[0]
        elif len(non_primary_ranges) > 1:
            for var in ["length", "vsb", "vgs", "vds"]:
                if var != primary and var in non_primary_ranges:
                    secondary = var
                    break

    slice_indices = (
        indices_dict["length"],
        indices_dict["vsb"],
        indices_dict["vgs"],
        indices_dict["vds"],
    )
    filter_values = {
            "length": values_dict["length"],
            "vsb": values_dict["vsb"],
            "vgs": values_dict["vgs"],
            "vds": values_dict["vds"],
    }

    extracted_table: dict = {}
    if parameters is None:
        parameters = lookup_table.get("parameter_names", [])
    for p in parameters:
        if p in lookup_table:
            data = np.squeeze(apply_slice(lookup_table[p], slice_indices))
            if data.ndim > 1 and data.shape[0] > data.shape[1]:
                data = data.T
            extracted_table[p] = data

    key = next(iter(extracted_table))
    extracted_table["width"] = np.array(width)
    extracted_table["length"], _ = tile_arrays(values_dict["length"], extracted_table[key])
    extracted_table["vsb"], _ = tile_arrays(values_dict["vsb"], extracted_table[key])
    extracted_table["vgs"], _ = tile_arrays(values_dict["vgs"], extracted_table[key])
    extracted_table["vds"], _ = tile_arrays(values_dict["vds"], extracted_table[key])

    return secondary, filter_values, extracted_table
