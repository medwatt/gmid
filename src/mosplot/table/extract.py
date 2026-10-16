from __future__ import annotations

import warnings

import numpy as np
import numpy.typing as npt


def tile_arrays(A: npt.NDArray, B: npt.NDArray) -> tuple[npt.NDArray, npt.NDArray]:
    """
    Broadcast a 1-D axis array to match the shape of a 2-D parameter array.

    This is used after slicing the 4-D lookup table to ensure that the voltage/
    length axes have the same shape as the parameter data, making vectorised
    expression evaluation straightforward.
    """
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
    elif A.ndim == 2 and B.ndim == 2:
        warnings.warn(
            f"tile_arrays: unexpected shapes A={A.shape}, B={B.shape}; returning unchanged"
        )
    return A, B


def extract_2d_table(
    lookup_table: dict,
    length: float | list | npt.NDArray | None = None,
    vbs: float | tuple[float, float] | tuple[float, float, float] | None = None,
    vgs: float | tuple[float, float] | tuple[float, float, float] | None = None,
    vds: float | tuple[float, float] | tuple[float, float, float] | None = None,
    primary: str | None = None,
    parameters: list | None = None,
) -> tuple[str | None, dict, dict]:
    """
    Slice the raw 4-D lookup table (length × vbs × vgs × vds) into a 2-D view.

    Two of the four axes are fixed to scalar values; the remaining two become the
    rows (primary sweep) and columns (secondary sweep) of the returned arrays.

    Returns
    -------
    secondary_var : str | None
        Name of the secondary sweep variable, or None when there is only one
        free axis.
    filter_values : dict
        The actual axis values selected for each of the four sweep variables.
    extracted_table : dict
        Parameter arrays reshaped to 2-D, plus tiled axis arrays for
        length, vbs, vgs, and vds so that expression evaluation works element-wise.
    """
    if sum(param is not None for param in [length, vbs, vgs, vds]) < 2:
        raise ValueError("Provide at least two of: length, vbs, vgs, vds.")

    def process_target(data: np.ndarray, target, var_name: str):
        if isinstance(target, tuple):
            if len(target) < 2:
                raise ValueError("Target tuple must have at least start and stop values.")

            start_val, stop_val = target[0], target[1]
            ascending = start_val <= stop_val

            start_idx = int(np.abs(data - start_val).argmin())
            stop_idx = int(np.abs(data - stop_val).argmin())
            # Index of the requested "stop" endpoint in the table. Must be
            # remembered separately: for a decreasing axis the start/stop
            # indices below get swapped, so inds[-1] no longer points at it.
            stop_end_idx = stop_idx

            if ascending and start_idx > stop_idx:
                start_idx, stop_idx = stop_idx, start_idx
            elif not ascending and start_idx < stop_idx:
                start_idx, stop_idx = stop_idx, start_idx

            inds = (
                np.arange(start_idx, stop_idx + 1)
                if ascending
                else np.arange(start_idx, stop_idx - 1, -1)
            )

            if len(target) == 3:
                step_val = target[2]
                if len(data) < 2:
                    raise ValueError("Data must have at least two elements to determine spacing.")
                dx = data[1] - data[0] if data[0] <= data[-1] else data[-2] - data[-1]
                index_step = int(round(abs(step_val) / dx))
                if index_step == 0:
                    raise ValueError("Step value is too small compared to data spacing.")
                step_idx = index_step if ascending else -index_step
                sliced_inds = inds[::step_idx]
                # Append the true stop endpoint. Using inds[-1] here appended
                # the *start* endpoint on decreasing axes, duplicating the
                # first sample at the very end and making the plotted line
                # jump back and close the curve.
                if sliced_inds[-1] != stop_end_idx:
                    sliced_inds = np.concatenate((sliced_inds, [stop_end_idx]))
                inds = sliced_inds

            return inds, data[inds], False

        if var_name == "length" and isinstance(target, (list, np.ndarray)):
            target_arr = np.array(target)
            mask = np.isin(data, target_arr)
            inds = np.nonzero(mask)[0]
            return inds, data[inds], False

        idx = int(np.abs(data - target).argmin())
        return np.array([idx]), data[idx], True

    def apply_slice(a: npt.NDArray, slices: tuple) -> npt.NDArray:
        s = a[slices[0], :, :, :]
        s = s[:, slices[1], :, :]
        s = s[:, :, slices[2], :]
        s = s[:, :, :, slices[3]]
        return s

    sweep_params = {"length": length, "vbs": vbs, "vgs": vgs, "vds": vds}
    indices_dict: dict = {}
    values_dict: dict = {}
    is_scalar: dict = {}

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

    secondary = None
    if primary:
        if sweep_params.get(primary) is None:
            raise ValueError("Primary variable must be provided.")
        non_primary_ranges = [v for v in sweep_params if v != primary and not is_scalar[v]]
        if len(non_primary_ranges) == 1:
            secondary = non_primary_ranges[0]
        elif len(non_primary_ranges) > 1:
            for v in ["length", "vbs", "vgs", "vds"]:
                if v != primary and v in non_primary_ranges:
                    secondary = v
                    break

    slice_indices = (
        indices_dict["length"],
        indices_dict["vbs"],
        indices_dict["vgs"],
        indices_dict["vds"],
    )
    filter_values = {
        "length": values_dict["length"],
        "vbs": values_dict["vbs"],
        "vgs": values_dict["vgs"],
        "vds": values_dict["vds"],
    }

    params: list = (
        parameters if parameters is not None else list(lookup_table.get("parameter_names", []))
    )
    keys_4d = [p for p in params if p in lookup_table and lookup_table[p].ndim == 4]
    if not keys_4d:
        raise ValueError("No parameter with 4 dimensions found for tiling.")
    reference_key = keys_4d[0]

    extracted_table: dict = {}
    for p in params:
        if p not in lookup_table:
            continue
        data = lookup_table[p]
        if p in keys_4d:
            data = np.squeeze(apply_slice(data, slice_indices))
        if data.ndim > 1 and data.shape[0] > data.shape[1]:
            data = data.T
        extracted_table[p] = data

    for var in ["length", "vbs", "vgs", "vds"]:
        extracted_table[var], _ = tile_arrays(values_dict[var], extracted_table[reference_key])

    return secondary, filter_values, extracted_table
