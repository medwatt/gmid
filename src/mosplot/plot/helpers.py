import numpy as np


def load_lookup_table(path: str) -> dict:
    """Load lookup table from a file."""
    return np.load(path, allow_pickle=True).tolist()


def tile_arrays(A: np.ndarray, B: np.ndarray) -> tuple:
    """Tile A to match the shape of B if needed."""
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


def extract_2d_table(lookup_table: dict, width, parameters=None,
                     lengths=None, vsb=None, vgs=None, vds=None, primary=None):
    """
    Filter the lookup table based on given sweep parameters.

    Args:
        lookup_table (dict): Dictionary of mosfet parameters.
        width: The mosfet width value.
        parameters (list, optional): List of parameter names to extract. Defaults to lookup_table["parameter_names"].
        lengths (float, list, or ndarray, optional): Length(s) of the mosfet.
        vsb (tuple, optional): Source-body voltage range (start, stop, step).
        vgs (tuple, optional): Gate-source voltage range (start, stop, step).
        vds (tuple, optional): Drain-source voltage range (start, stop, step).
        primary (str, optional): Name of the primary sweep variable.

    Returns:
        tuple: (secondary_idx, filter_values, extracted_table)
            - secondary_idx: Index of the secondary sweep variable (if applicable).
            - filter_values: Filtered values for lengths, vsb, vgs, vds.
            - extracted_table: Dictionary with extracted parameter arrays.
    """
    # Require at least two sweep parameters.
    params = [lengths is not None, vsb is not None, vgs is not None, vds is not None]
    if sum(params) < 2:
        raise ValueError("Please provide at least two parameters.")

    variables = {"lengths": False, "vsb": False, "vgs": False, "vds": False}
    if primary:
        variables[primary] = True

    def get_indices(var, target):
        data = lookup_table[var]
        if isinstance(target, tuple):
            # When target is a range: (start, stop, [step])
            start, end = target[:2]
            indices = np.where((data >= start) & (data <= end))[0]
            if len(target) == 3:  # If step is provided.
                step = int(target[2] / (data[1] - data[0]))
                indices = indices[::step]
        elif var == "lengths" and isinstance(target, (list, np.ndarray)):
            mask = np.isin(lookup_table["lengths"], np.array(target))
            indices = np.nonzero(mask)[0]
            indices = np.array(indices, dtype=int)
        else:
            # Mark variable as secondary.
            variables[var] = True
            index = (np.abs(data - target)).argmin()
            return np.array([index]), data[index]
        return indices, data[indices]

    secondary_idx = None
    indices_and_values = {
        "lengths": get_indices("lengths", lengths) if lengths is not None else (slice(None), lookup_table["lengths"]),
        "vsb": get_indices("vsb", vsb) if vsb is not None else (slice(None), lookup_table["vsb"]),
        "vgs": get_indices("vgs", vgs) if vgs is not None else (slice(None), lookup_table["vgs"]),
        "vds": get_indices("vds", vds) if vds is not None else (slice(None), lookup_table["vds"]),
    }

    slice_indices = []
    filter_values = []
    for key in variables.keys():
        slice_indices.append(indices_and_values[key][0])
        filter_values.append(indices_and_values[key][1])
    slice_indices = tuple(slice_indices)

    def slice_me(a, slices):
        x = a[slices[0], :, :, :]
        x = x[:, slices[1], :, :]
        x = x[:, :, slices[2], :]
        x = x[:, :, :, slices[3]]
        return x

    extracted_table = {}
    if not parameters:
        parameters = lookup_table["parameter_names"]

    for p in parameters:
        if p in lookup_table:
            x = np.squeeze(slice_me(lookup_table[p], slice_indices))
            extracted_table[p] = x.T if (x.ndim > 1 and x.shape[0] > x.shape[1]) else x

    one_key = next(iter(extracted_table))
    extracted_table["width"] = np.array(width)
    extracted_table["lengths"], _ = tile_arrays(filter_values[0], extracted_table[one_key])
    extracted_table["vsb"], _ = tile_arrays(filter_values[1], extracted_table[one_key])
    extracted_table["vgs"], _ = tile_arrays(filter_values[2], extracted_table[one_key])
    extracted_table["vds"], _ = tile_arrays(filter_values[3], extracted_table[one_key])

    if primary:
        secondary_idx = list(variables.values()).index(False)

    return secondary_idx, filter_values, extracted_table
