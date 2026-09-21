"""GmIdTable: preprocessed 4-D lookup table for gm/Id-based circuit optimization.

Raw simulator data is swept over (L, vbs, vds, vgs). This module remaps that
data onto a regular grid indexed by gmid = gm/Id instead of vgs, because gmid
is the natural design variable in the gm/Id methodology. Subsequent lookups
are O(1) trilinear interpolation rather than O(n) nearest-neighbour search.

Two forward lookup paths:
  lookup_scalar  -- single operating point; Numba-JIT compiled (~7 µs/call)
  lookup         -- vectorised NumPy; efficient for arrays of query points

One reverse curve-extraction path:
  lookup_curve   -- extract full 1-D parameter curve along one axis for
                    reverse lookup (see FastMosfet.lookup_gmid_for etc.)

Axis validation:
  _require_strictly_increasing is called at construction and cache-load time
  for all four axes.  np.searchsorted-based bracketing requires sorted axes;
  an unsorted axis silently produces wrong lookup results, so we fail loud.

Disk caching is supported: the first run builds and saves the grid as a .npz
file; later runs load it directly, skipping the resampling step.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np

from .interp import _resample_into, _lookup_scalar_nb, _bracket_vec, _interp4d_nan, _bracket_nb, _curve_for_param_nb
from mosplot.interpolation.axes import _AXIS_DIMS, _require_strictly_increasing


# ---------------------------------------------------------------------------
# Cache helpers
# ---------------------------------------------------------------------------


def _as_cache_dir(cache_dir: str | Path) -> Path:
    path = Path(cache_dir).expanduser()
    path.mkdir(parents=True, exist_ok=True)
    return path


def _json_default(obj: Any):
    if isinstance(obj, np.generic):
        return obj.item()
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, Path):
        return str(obj)
    raise TypeError(f"Object of type {type(obj).__name__} is not JSON serializable")


def _canonical_gmid_bounds(gmid_bounds):
    if gmid_bounds is None:
        return None
    return [float(gmid_bounds[0]), float(gmid_bounds[1])]


def _hash_array(hasher, name, value):
    arr = np.asarray(value)
    hasher.update(str(name).encode("utf-8"))
    hasher.update(str(arr.shape).encode("utf-8"))
    hasher.update(str(arr.dtype).encode("utf-8"))
    if arr.dtype == object:
        payload = json.dumps(arr.tolist(), sort_keys=True, default=_json_default)
        hasher.update(payload.encode("utf-8"))
    else:
        hasher.update(np.ascontiguousarray(arr).view(np.uint8))


def _device_fingerprint(lookup_table, device_key):
    dev = lookup_table[device_key]
    hasher = hashlib.sha1()
    hasher.update(str(device_key).encode("utf-8"))

    for axis in ("length", "vbs", "vds", "vgs"):
        if axis in dev:
            _hash_array(hasher, axis, dev[axis])

    parameter_names = list(dev.get("parameter_names", []))
    hasher.update(json.dumps(parameter_names, default=_json_default).encode("utf-8"))
    for param in parameter_names:
        if param in dev:
            _hash_array(hasher, param, dev[param])

    return hasher.hexdigest()[:12]


def _cache_signature(device_key, n_gmid, gmid_bounds, cache_tag=None, source_fingerprint=None):
    """Dict hashed into the cache filename. Change cache_tag to invalidate old caches."""
    return {
        "device_key": str(device_key),
        "n_gmid": int(n_gmid),
        "gmid_bounds": _canonical_gmid_bounds(gmid_bounds),
        "cache_tag": None if cache_tag is None else str(cache_tag),
        "source_fingerprint": source_fingerprint,
    }


def _cache_name(device_key, n_gmid, gmid_bounds, cache_tag=None, source_fingerprint=None):
    signature = _cache_signature(
        device_key, n_gmid, gmid_bounds, cache_tag, source_fingerprint
    )
    payload = json.dumps(signature, sort_keys=True, default=_json_default)
    digest = hashlib.sha1(payload.encode("utf-8")).hexdigest()[:12]
    safe_key = "".join(c if c.isalnum() or c in "._-" else "_" for c in str(device_key))
    return f"{safe_key}_ng{int(n_gmid)}_{digest}.npz"


def _load_lookup_table(lookup_table_or_npz):
    """Accept either an already-loaded lookup_table dict or a path to a .npz file."""
    if isinstance(lookup_table_or_npz, (str, Path)):
        npz_path = Path(lookup_table_or_npz).expanduser()
        return np.load(npz_path, allow_pickle=True)["lookup_table"].item()
    return lookup_table_or_npz


# ---------------------------------------------------------------------------
# GmIdTable
# ---------------------------------------------------------------------------


class GmIdTable:
    """Preprocessed (L, vbs, vds, gmid) -> device-parameter lookup table.

    Construction
    ------------
    The raw simulator sweep is over (L, vbs, vds, vgs). This class remaps it
    to a regular grid indexed by gmid = gm/Id instead of vgs, because gmid is
    the designer's natural control variable in the gm/Id methodology.

    The gmid axis is log-spaced: gmid spans roughly 1–30 S/A from strong to
    weak inversion, and equal spacing in log captures both ends faithfully.

    gm/Id is not monotone over the full vgs sweep -- it peaks near the
    weak-inversion boundary then falls toward strong inversion. For each
    (L, vbs, vds) slice we locate that peak and keep only the monotonically
    decreasing branch. This makes the axis invertible (one gmid -> one device
    state). The post-peak branch (higher |VGS|) is the saturation/design
    region. VGS must be swept from 0 toward VDD (NMOS) or −VDD (PMOS) so
    that this branch is always the post-peak side.

    NaN handling
    ------------
    Not every GMID value is achievable at every (L, VBS, VDS) point -- for
    example, a very small |VDS| limits how deeply into strong inversion the
    device can be driven. Grid cells outside the achievable range stay NaN
    (the output grid is pre-filled with NaN). The interpolation routines
    exclude NaN corners from the weighted average so that lookups near the
    boundary of the valid operating region degrade gracefully.

    Disk caching
    ------------
    Use build_or_load() to transparently cache the prebuilt grid to disk.
    Cache filenames encode device_key, n_gmid, gmid_bounds, and cache_tag;
    if the raw lookup table changes without those parameters changing, pass
    rebuild=True or bump cache_tag.
    """

    def __init__(self, lookup_table, device_key, n_gmid=50, gmid_bounds=None):
        self.device_key = str(device_key)
        self.n_gmid = int(n_gmid)
        self.gmid_bounds = _canonical_gmid_bounds(gmid_bounds)

        dev = lookup_table[device_key]

        lengths = np.asarray(dev["length"], dtype=float)
        vgs_axis = np.asarray(dev["vgs"], dtype=float)
        vds_axis = np.asarray(dev["vds"], dtype=float)
        vbs_axis = np.asarray(dev["vbs"], dtype=float)

        # Sort vbs/vds to ascending order; track permutation for data reordering.
        vbs_order = np.argsort(vbs_axis)
        vbs_axis = vbs_axis[vbs_order]
        vds_order = np.argsort(vds_axis)
        vds_axis = vds_axis[vds_order]

        n_L, n_vbs, n_vds, n_vgs = len(lengths), len(vbs_axis), len(vds_axis), len(vgs_axis)

        # Simulators may store axes as (L, vbs, vds, vgs) or (L, vbs, vgs, vds).
        # Detect the ordering from the actual data shape.
        expected_std = (n_L, n_vbs, n_vds, n_vgs)
        expected_swp = (n_L, n_vbs, n_vgs, n_vds)

        first_param = dev["parameter_names"][0]
        test_shape = np.asarray(dev[first_param]).shape

        if test_shape not in (expected_std, expected_swp):
            raise ValueError(
                f"Data shape {test_shape} does not match "
                f"expected ({expected_std} or {expected_swp})"
            )

        if expected_std == expected_swp:
            raise ValueError(
                f"VGS axis ({n_vgs} points) and VDS axis ({n_vds} points) have "
                f"the same length for {device_key}. Use a unidirectional sweep "
                f"(VGS from 0 toward VDD) with a different number of points on "
                f"each axis to avoid ambiguity."
            )
        swap_inner = test_shape == expected_swp

        # Keep only 4-D parameters (scalars and other shapes are skipped).
        param_names = [
            p
            for p in dev["parameter_names"]
            if np.asarray(dev[p]).ndim == 4
            and np.asarray(dev[p]).shape in (expected_std, expected_swp)
        ]

        # The gmid axis itself is built from gm/id, so these two are the only
        # hard requirements; fail with a clear message instead of a KeyError.
        missing = sorted({"id", "gm"} - set(param_names))
        if missing:
            raise ValueError(
                f"Lookup table for {device_key!r} is missing {missing}. "
                f"'id' and 'gm' must be included in parameters_to_save when "
                f"generating the table; found: {param_names}."
            )
        all_params = param_names + ["vgs"]
        n_params = len(all_params)

        # Apply axis reordering so all data is (L, vbs, vds, vgs).
        raw: dict[str, np.ndarray] = {}
        for p in param_names:
            arr = np.asarray(dev[p], dtype=float)[:, vbs_order, :, :]
            arr = arr.transpose(0, 1, 3, 2) if swap_inner else arr
            raw[p] = arr[:, :, vds_order, :]

        # Broadcast vgs into 4-D so it can be resampled like any other parameter
        # (needed to recover vgs from a gmid query).
        raw["vgs"] = np.broadcast_to(
            vgs_axis[np.newaxis, np.newaxis, np.newaxis, :],
            (n_L, n_vbs, n_vds, n_vgs),
        ).copy()

        with np.errstate(divide="ignore", invalid="ignore"):
            gmid_raw = np.where(
                np.abs(raw["id"]) > 1e-20,
                raw["gm"] / raw["id"],
                np.nan,
            )

        # Determine gmid grid bounds from data if not supplied.
        if gmid_bounds is None:
            finite_pos = gmid_raw[np.isfinite(gmid_raw) & (gmid_raw > 0)]
            if finite_pos.size == 0:
                raise ValueError(f"No positive finite gm/id values for {device_key}")
            g_hi = float(np.percentile(finite_pos, 99))
            g_lo = float(max(np.percentile(finite_pos, 1), 0.5))
        else:
            g_lo, g_hi = float(gmid_bounds[0]), float(gmid_bounds[1])

        if not (np.isfinite(g_lo) and np.isfinite(g_hi) and g_lo > 0 and g_hi > g_lo):
            raise ValueError(f"Invalid gm/id bounds: {(g_lo, g_hi)}")

        # Log-spaced grid: gmid spans orders of magnitude across inversion regions.
        gmid_grid = np.exp(np.linspace(np.log(g_lo), np.log(g_hi), n_gmid))

        raw_stack = np.stack([raw[p] for p in all_params], axis=0)
        out = np.full((n_params, n_L, n_vbs, n_vds, n_gmid), np.nan)

        for li in range(n_L):
            for vi in range(n_vbs):
                for di in range(n_vds):
                    g = gmid_raw[li, vi, di, :]
                    if not np.any(np.isfinite(g) & (g > 0)):
                        continue

                    # For a unidirectional sweep (VGS from 0 toward VDD/−VDD),
                    # the gmid peak separates the sub-threshold turn-on region
                    # (pre-peak, low VGS) from the saturation design region
                    # (post-peak, higher |VGS|). Always take the post-peak branch.
                    peak = int(np.nanargmax(g))
                    g_dec = g[peak:]
                    d_raw = raw_stack[:, li, vi, di, peak:]

                    # Flip to ascending order for searchsorted in _resample_into.
                    g_asc = g_dec[::-1]
                    d_mat = d_raw[:, ::-1]
                    _resample_into(g_asc, d_mat, gmid_grid, out[:, li, vi, di, :])

        if n_vbs == 1:
            # _bracket_nb requires at least two points to form a valid interval.
            eps = 1e-10
            vbs_axis = np.array([vbs_axis[0] - eps, vbs_axis[0] + eps])
            out = np.repeat(out, 2, axis=2)

        # _data_stack is a contiguous C-order array required by the Numba path.
        # _data shares the same underlying memory via views (no copy cost).
        self._data_stack = np.ascontiguousarray(out)
        self._data = {p: self._data_stack[i] for i, p in enumerate(all_params)}
        self._param_idx = {p: i for i, p in enumerate(all_params)}
        self.lengths = lengths
        self.vbs = vbs_axis
        self.vds = vds_axis
        self.gmid = gmid_grid
        self.params = all_params

        # Validate axes at construction time.  np.searchsorted (used by all
        # bracketing, forward and reverse) requires strictly increasing axes.
        # An unsorted axis silently produces wrong interpolation results.
        _require_strictly_increasing(self.lengths, "length")
        _require_strictly_increasing(self.vbs, "vbs")
        _require_strictly_increasing(self.vds, "vds")
        _require_strictly_increasing(self.gmid, "gmid")

    # --- cache API --------------------------------------------------------------

    @staticmethod
    def cache_path(
        cache_dir,
        device_key,
        n_gmid=50,
        gmid_bounds=None,
        cache_tag=None,
        source_fingerprint=None,
    ):
        cache_dir = _as_cache_dir(cache_dir)
        return cache_dir / _cache_name(
            device_key,
            n_gmid,
            gmid_bounds,
            cache_tag,
            source_fingerprint=source_fingerprint,
        )

    @classmethod
    def build_or_load(
        cls,
        lookup_table,
        device_key,
        *,
        cache_dir,
        n_gmid=50,
        gmid_bounds=None,
        cache_tag=None,
        rebuild=False,
        compressed=False,
        verbose=False,
    ):
        """Load from disk if a valid cache exists, otherwise build and save."""
        source_fingerprint = _device_fingerprint(lookup_table, device_key)
        cache_path = cls.cache_path(
            cache_dir,
            device_key,
            n_gmid=n_gmid,
            gmid_bounds=gmid_bounds,
            cache_tag=cache_tag,
            source_fingerprint=source_fingerprint,
        )
        if cache_path.exists() and not rebuild:
            try:
                if verbose:
                    print(f"Loading cached GmIdTable: {cache_path}")
                return cls.load(cache_path)
            except Exception as exc:
                if verbose:
                    print(f"Cache load failed ({exc}). Rebuilding.")

        if verbose:
            print(f"Building GmIdTable for {device_key}")
        table = cls(lookup_table, device_key, n_gmid=n_gmid, gmid_bounds=gmid_bounds)
        table.save(
            cache_path,
            metadata=_cache_signature(
                device_key,
                n_gmid,
                gmid_bounds,
                cache_tag,
                source_fingerprint=source_fingerprint,
            ),
            compressed=compressed,
        )
        if verbose:
            print(f"Saved GmIdTable: {cache_path}")
        return table

    def save(self, path, *, metadata=None, compressed=False):
        """Serialise the prebuilt grid to a .npz file."""
        path = Path(path).expanduser()
        path.parent.mkdir(parents=True, exist_ok=True)
        payload = {
            "__metadata__": np.array(
                json.dumps(metadata or {}, sort_keys=True, default=_json_default),
                dtype=object,
            ),
            "device_key": np.array(getattr(self, "device_key", ""), dtype=object),
            "n_gmid": np.array(getattr(self, "n_gmid", len(self.gmid))),
            "gmid_bounds": np.array(getattr(self, "gmid_bounds", None), dtype=object),
            "lengths": self.lengths,
            "vbs": self.vbs,
            "vds": self.vds,
            "gmid": self.gmid,
            "params": np.array(self.params, dtype=object),
        }
        for p, arr in self._data.items():
            payload[f"data__{p}"] = arr
        saver = np.savez_compressed if compressed else np.savez
        saver(path, **payload)

    @classmethod
    def load(cls, path):
        """Deserialise a prebuilt grid from a .npz file without rebuilding."""
        path = Path(path).expanduser()
        z = np.load(path, allow_pickle=True)
        obj = cls.__new__(cls)
        obj.device_key = str(z["device_key"].item()) if "device_key" in z else ""
        obj.n_gmid = int(z["n_gmid"].item()) if "n_gmid" in z else int(len(z["gmid"]))
        obj.gmid_bounds = z["gmid_bounds"].tolist() if "gmid_bounds" in z else None
        obj.lengths = z["lengths"]
        obj.vbs = z["vbs"]
        obj.vds = z["vds"]
        obj.gmid = z["gmid"]
        obj.params = [str(p) for p in z["params"].tolist()]
        obj._data_stack = np.stack([z[f"data__{p}"] for p in obj.params], axis=0)
        obj._data = {p: obj._data_stack[i] for i, p in enumerate(obj.params)}
        z.close()
        obj._param_idx = {p: i for i, p in enumerate(obj.params)}

        # Same validation as __init__ -- cached axes must still be valid.
        _require_strictly_increasing(obj.lengths, "length")
        _require_strictly_increasing(obj.vbs, "vbs")
        _require_strictly_increasing(obj.vds, "vds")
        _require_strictly_increasing(obj.gmid, "gmid")

        return obj

    # --- lookup API -------------------------------------------------------------

    def lookup_scalar(self, L, gmid, vds, vbs, params=None):
        """Single operating-point lookup (optimizer hot path).

        Returns a dict mapping parameter name -> float. Numba eliminates the
        per-call NumPy dispatch overhead, giving ~5 µs/call vs ~80 µs without.
        """
        if params is None:
            params = self.params
        raw = _lookup_scalar_nb(
            self._data_stack,
            self.lengths,
            self.vbs,
            self.vds,
            self.gmid,
            float(L),
            float(gmid),
            float(vds),
            float(vbs),
        )
        return {p: float(raw[self._param_idx[p]]) for p in params}

    def lookup(self, L, gmid, vds, vbs, params=None):
        """Vectorised lookup: all four inputs may be arrays (broadcast-compatible)."""
        L, gmid, vds, vbs = np.broadcast_arrays(
            np.asarray(L, dtype=float),
            np.asarray(gmid, dtype=float),
            np.asarray(vds, dtype=float),
            np.asarray(vbs, dtype=float),
        )
        shape = L.shape
        il, il1, fl = _bracket_vec(self.lengths, L.ravel())
        iv, iv1, fv = _bracket_vec(self.vbs, vbs.ravel())
        id_, id1_, fd = _bracket_vec(self.vds, vds.ravel())
        ig, ig1, fg = _bracket_vec(self.gmid, gmid.ravel())
        if params is None:
            params = self.params
        return {
            p: _interp4d_nan(
                self._data[p], il, il1, fl, iv, iv1, fv, id_, id1_, fd, ig, ig1, fg
            ).reshape(shape)
            for p in params
        }

    # --- curve extraction (reverse lookup support) ----------------------------
    #
    # lookup_curve extracts the full 1-D parameter curve along one free axis
    # for use by the reverse-lookup API (FastMosfet.lookup_gmid_for etc.).
    #
    # For each raw parameter, the 3 fixed axes are bracketed once (Numba
    # _bracket_nb), then _curve_for_param_nb does an 8-corner trilinear
    # interpolation at every free-axis point -- a single Numba JIT call that
    # produces the full curve in ~20 µs.  The result is used by the 1-D
    # inverse solver to find crossings.

    def _axis_values(self, axis_name: str):
        if axis_name == "length":
            return self.lengths
        if axis_name == "vbs":
            return self.vbs
        if axis_name == "vds":
            return self.vds
        if axis_name == "gmid":
            return self.gmid
        raise KeyError(f"Unknown axis: {axis_name!r}")

    def _fixed_axis_values_for_curve(
        self,
        solve_for: str,
        *,
        length=None,
        gmid=None,
        vds=None,
        vbs=None,
    ):
        if solve_for not in _AXIS_DIMS:
            raise ValueError(
                f"Unknown solve axis: {solve_for!r}. "
                f"Expected one of {sorted(_AXIS_DIMS)}."
            )

        values = {
            "length": length,
            "gmid": gmid,
            "vds": vds,
            "vbs": vbs,
        }

        if values[solve_for] is not None:
            raise ValueError(
                f"Axis {solve_for!r} is being solved for and must not be provided."
            )

        required = [name for name in _AXIS_DIMS if name != solve_for]
        missing = [name for name in required if values[name] is None]

        if missing:
            raise ValueError(
                f"Missing fixed axes for solve_for={solve_for!r}: {missing}"
            )

        return {name: float(values[name]) for name in required}

    def lookup_curve(
        self,
        solve_for: str,
        *,
        params,
        length=None,
        gmid=None,
        vds=None,
        vbs=None,
    ):
        """Return raw parameter curves along one free axis.

        Parameters
        ----------
        solve_for : {"length", "gmid", "vds", "vbs"}
            Axis to leave free.
        params : sequence[str]
            Raw table parameter names to extract as curves. May be empty.
        length, gmid, vds, vbs : float or None
            The three fixed axes must be provided. The solved axis must be None.

        Returns
        -------
        axis : np.ndarray
            Values of the solved axis.
        curves : dict[str, np.ndarray]
            Mapping from raw parameter name to 1-D curve over axis.
        """
        if solve_for not in _AXIS_DIMS:
            raise ValueError(
                f"Unknown solve axis: {solve_for!r}. "
                f"Expected one of {sorted(_AXIS_DIMS)}."
            )

        if isinstance(params, str):
            params = [params]
        else:
            params = list(params)

        for p in params:
            if p not in self._data:
                raise KeyError(f"Unknown table parameter: {p!r}")

        fixed = self._fixed_axis_values_for_curve(
            solve_for,
            length=length,
            gmid=gmid,
            vds=vds,
            vbs=vbs,
        )

        free_axis = self._axis_values(solve_for)
        free_dim = _AXIS_DIMS[solve_for]

        brackets = {}
        for axis_name, value in fixed.items():
            ax = self._axis_values(axis_name)
            i0, frac = _bracket_nb(ax, value)
            brackets[axis_name] = (int(i0), int(i0) + 1, float(frac))

        fixed_axis_names = list(fixed.keys())

        # Pre-compute 8-corner indices and weights for the Numba kernel
        corner_ia = np.empty(8, dtype=np.int64)
        corner_ib = np.empty(8, dtype=np.int64)
        corner_ic = np.empty(8, dtype=np.int64)
        corner_wt = np.empty(8, dtype=np.float64)

        dims = [_AXIS_DIMS[name] for name in fixed_axis_names]
        dfA, dfB, dfC = dims

        i0a, i1a, fa = brackets[fixed_axis_names[0]]
        i0b, i1b, fb = brackets[fixed_axis_names[1]]
        i0c, i1c, fc = brackets[fixed_axis_names[2]]

        w0a = 1.0 - fa; w1a = fa
        w0b = 1.0 - fb; w1b = fb
        w0c = 1.0 - fc; w1c = fc

        corner_ia[0]=i0a; corner_ib[0]=i0b; corner_ic[0]=i0c; corner_wt[0]=w0a*w0b*w0c
        corner_ia[1]=i1a; corner_ib[1]=i0b; corner_ic[1]=i0c; corner_wt[1]=w1a*w0b*w0c
        corner_ia[2]=i0a; corner_ib[2]=i1b; corner_ic[2]=i0c; corner_wt[2]=w0a*w1b*w0c
        corner_ia[3]=i0a; corner_ib[3]=i0b; corner_ic[3]=i1c; corner_wt[3]=w0a*w0b*w1c
        corner_ia[4]=i1a; corner_ib[4]=i1b; corner_ic[4]=i0c; corner_wt[4]=w1a*w1b*w0c
        corner_ia[5]=i1a; corner_ib[5]=i0b; corner_ic[5]=i1c; corner_wt[5]=w1a*w0b*w1c
        corner_ia[6]=i0a; corner_ib[6]=i1b; corner_ic[6]=i1c; corner_wt[6]=w0a*w1b*w1c
        corner_ia[7]=i1a; corner_ib[7]=i1b; corner_ic[7]=i1c; corner_wt[7]=w1a*w1b*w1c

        curves = {
            p: _curve_for_param_nb(
                self._data[p],
                len(free_axis), free_dim,
                corner_ia[0], corner_ia[1], corner_ia[2], corner_ia[3],
                corner_ia[4], corner_ia[5], corner_ia[6], corner_ia[7],
                corner_ib[0], corner_ib[1], corner_ib[2], corner_ib[3],
                corner_ib[4], corner_ib[5], corner_ib[6], corner_ib[7],
                corner_ic[0], corner_ic[1], corner_ic[2], corner_ic[3],
                corner_ic[4], corner_ic[5], corner_ic[6], corner_ic[7],
                corner_wt[0], corner_wt[1], corner_wt[2], corner_wt[3],
                corner_wt[4], corner_wt[5], corner_wt[6], corner_wt[7],
                dfA, dfB, dfC,
            )
            for p in params
        }

        return free_axis, curves


# ---------------------------------------------------------------------------
# Public prebuild helper
# ---------------------------------------------------------------------------


def prebuild_fast_tables(
    lookup_table_or_npz,
    *,
    cache_dir,
    device_keys=None,
    n_gmid=50,
    gmid_bounds=None,
    cache_tag=None,
    rebuild=False,
    compressed=False,
    verbose=True,
):
    """Build and save GmIdTable cache files for all (or selected) devices.

    Parameters
    ----------
    lookup_table_or_npz : dict | str | Path
        Already-loaded lookup_table dict, or a path to a .npz file.
    cache_dir : str | Path
        Directory where .npz cache files are written.
    device_keys : list[str] | None
        Devices to build. If None, all devices in the table are built.
    n_gmid, gmid_bounds
        Passed through to GmIdTable.
    cache_tag : str | None
        Optional tag included in the cache filename. Bump this when the raw
        table is updated but device names and n_gmid settings stay the same.
    rebuild : bool
        If True, overwrite existing cache files.
    compressed : bool
        If True, write compressed .npz files (slower to load).
    verbose : bool
        Print progress messages.

    Returns
    -------
    dict
        Mapping {device_key: cache_path}.
    """
    lookup_table = _load_lookup_table(lookup_table_or_npz)
    if device_keys is None:
        device_keys = list(lookup_table.keys())

    paths = {}
    for device_key in device_keys:
        GmIdTable.build_or_load(
            lookup_table,
            device_key,
            cache_dir=cache_dir,
            n_gmid=n_gmid,
            gmid_bounds=gmid_bounds,
            cache_tag=cache_tag,
            rebuild=rebuild,
            compressed=compressed,
            verbose=verbose,
        )
        # The cache filename includes the data fingerprint, so it must be
        # recomputed with the same fingerprint build_or_load used.
        paths[device_key] = GmIdTable.cache_path(
            cache_dir,
            device_key,
            n_gmid=n_gmid,
            gmid_bounds=gmid_bounds,
            cache_tag=cache_tag,
            source_fingerprint=_device_fingerprint(lookup_table, device_key),
        )

    return paths
