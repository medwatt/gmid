from __future__ import annotations

from typing import overload

from mosplot.expressions import Expression, build_expressions, compute_device_width
from mosplot.interpolation import GmIdTable
from mosplot.optimizer.fast_mosfet_reverse import _lookup_axis_for_impl


class FastMosfet:
    """Fast MOSFET lookup for circuit optimization.

    Wraps GmIdTable to provide expression-based lookups indexed by
    (length, gmid, vds, vbs). Designed for optimization loops where hundreds
    of lookups are performed per iteration.

    Construction
    ------------
        # Build the grid in memory:
        pmos = FastMosfet(lookup_table, "pch_lvt")

        # Build once and cache to disk (subsequent runs load instantly):
        pmos = FastMosfet(lookup_table, "pch_lvt", cache_dir="~/.cache/mosplot")

        # Force a cache rebuild:
        pmos = FastMosfet(lookup_table, "pch_lvt", cache_dir="~/.cache/mosplot", rebuild_cache=True)

    Parameters
    ----------
    cache_dir : str | Path | None
        Directory for cached GmIdTable .npz files. Omit to build in memory.
    rebuild_cache : bool
        Ignore an existing cache and force a rebuild.
    cache_tag : str | None
        Tag included in the cache filename. Bump this when the raw lookup
        table changes but device names and n_gmid settings do not.
    compressed_cache : bool
        Write compressed .npz files (slower to load). Default False.
    verbose : bool
        Print cache load/build messages.
    """

    # Class-level annotations for LSP resolution -- populated by build_expressions()
    length_expression: Expression
    vgs_expression: Expression
    vsg_expression: Expression
    vbs_expression: Expression
    vsb_expression: Expression
    vds_expression: Expression
    vsd_expression: Expression
    vdsat_expression: Expression
    vth_expression: Expression
    vov_expression: Expression
    vstar_expression: Expression
    id_expression: Expression
    gm_expression: Expression
    gmbs_expression: Expression
    gds_expression: Expression
    cgg_expression: Expression
    cgs_expression: Expression
    cgd_expression: Expression
    cbg_expression: Expression
    cdd_expression: Expression
    gmid_expression: Expression
    gain_expression: Expression
    transit_frequency_expression: Expression
    early_voltage_expression: Expression
    inverse_early_voltage_expression: Expression
    current_density_expression: Expression
    rds_expression: Expression

    def __init__(
        self,
        lookup_table,
        mos_name,
        n_gmid=50,
        gmid_bounds=None,
        *,
        cache_dir=None,
        rebuild_cache=False,
        cache_tag=None,
        compressed_cache=False,
        verbose=False,
    ):
        if cache_dir is None:
            self._table = GmIdTable(
                lookup_table,
                mos_name,
                n_gmid=n_gmid,
                gmid_bounds=gmid_bounds,
            )
        else:
            self._table = GmIdTable.build_or_load(
                lookup_table,
                mos_name,
                cache_dir=cache_dir,
                n_gmid=n_gmid,
                gmid_bounds=gmid_bounds,
                cache_tag=cache_tag,
                rebuild=rebuild_cache,
                compressed=compressed_cache,
                verbose=verbose,
            )

        dev = lookup_table[mos_name]
        dev_params = dev.get("device_parameters", {})
        param_names = list(dev.get("parameter_names", []))

        self._width = compute_device_width(dev_params)
        self._param_names = param_names
        self._mos_name = mos_name

        vdsat_var = "vdssat" if "vdssat" in param_names else "vdsat"
        for name, expr in build_expressions(self._width, vdsat_var).items():
            setattr(self, name, expr)

    @overload
    def interpolate(
        self, length: float, gmid: float, vds: float, vbs: float, expression: Expression
    ) -> float: ...
    @overload
    def interpolate(
        self, length: float, gmid: float, vds: float, vbs: float, expression: list[Expression]
    ) -> list[float]: ...

    def interpolate(
        self,
        length: float,
        gmid: float,
        vds: float,
        vbs: float,
        expression: Expression | list[Expression],
    ) -> float | list[float]:
        """Look up one or more expressions at a single operating point.

        Parameters
        ----------
        length : float
            Channel length (metres).
        gmid : float
            Target gm/Id ratio (S/A).
        vds : float
            Drain-source voltage (V).
        vbs : float
            Body-source voltage (V).
        expression : Expression or list[Expression]
            Quantity or quantities to look up.

        Returns
        -------
        float if a single Expression is given, list[float] otherwise.
        """
        if isinstance(expression, list):
            exprs = expression
            single = False
        else:
            exprs = [expression]
            single = True

        needed = list({v for e in exprs for v in e.variables if v != "length"})
        raw = self._table.lookup_scalar(length, gmid, vds, vbs, params=needed)

        results = [expr.function(*[raw[v] for v in expr.variables]) for expr in exprs]
        return results[0] if single else results

    # --- reverse lookup --------------------------------------------------------
    #
    # Public API for finding an axis value that makes one expression reach a
    # target.  Five methods:
    #
    #   lookup_axis_for   -- generic: specify which axis to solve and which
    #                        three axes are fixed.
    #   lookup_gmid_for   -- convenience: solve for gmid
    #   lookup_length_for -- convenience: solve for length
    #   lookup_vds_for    -- convenience: solve for vds
    #   lookup_vbs_for    -- convenience: solve for vbs
    #
    # All delegate to _lookup_axis_for_impl (in fast_mosfet_reverse.py),
    # which extracts the full 1-D parameter curve, evaluates the expression
    # along it, and solves y(axis) == target via the inverse solver.
    #
    # The wrappers contain no logic -- they only fill in the omitted axis and
    # forward the call.
    #
    # Solve method
    # ------------
    # The ``method`` parameter selects the inverse solver:
    #
    #   "scan" (default)
    #       Full 1-D curve scan via the Python _inverse_1d.  Handles
    #       every mode, out_of_range policy, NaN-gap segmentation, and
    #       multiple-crossing detection/warning.  Suitable for any
    #       expression in any operating region.
    #
    #   "fast"
    #       Numba-JIT single-pass scan (_inverse_1d_nb).  ~4× faster
    #       than "scan" (~26 µs vs ~100 µs) but imposes restrictions:
    #       only ``mode="first"`` and ``out_of_range`` in {"nan","clip"}
    #       are supported; multiple-crossing warnings are disabled.
    #       Use when the expression is known to be monotonic (e.g.
    #       current_density vs gmid) and you need maximum throughput
    #       inside an optimization loop.
    #
    #   Passing ``method="fast"`` with an incompatible option raises
    #   ValueError.

    def lookup_axis_for(
        self,
        solve_for: str,
        *,
        expression: Expression,
        target: float,
        length: float | None = None,
        gmid: float | None = None,
        vds: float | None = None,
        vbs: float | None = None,
        mode: str = "first",
        out_of_range: str = "nan",
        warn_multiple: bool = True,
        method: str = "scan",
    ) -> float | list[float]:
        """Find one axis value such that expression(...) == target.

        Parameters
        ----------
        solve_for : {"length", "gmid", "vds", "vbs"}
            Axis to solve for.
        expression : Expression
            Single expression to match.  Lists are not allowed.
        target : float
            Desired expression value.
        length, gmid, vds, vbs : float or None
            Provide the three fixed axes.  Leave the solved axis as None.
        mode : {"first", "last", "nearest", "all"}
            Which crossing to return (only "first" in "fast" mode).
        out_of_range : {"nan", "clip", "raise"}
            What to do if no crossing exists (only "nan"/"clip" in "fast"
            mode).
        warn_multiple : bool
            If True, warn when multiple crossings exist (ignored in "fast"
            mode -- warnings are always suppressed there).
        method : {"scan", "fast"}
            Solve method.  "scan" uses the full Python inverse solver with
            all features.  "fast" uses the Numba-JIT solver for ~4× speed
            but with the restrictions listed above.

        Returns
        -------
        float or list[float]
            Solved axis value.  Returns a list only when mode="all"
            (which requires method="scan").
        """
        return _lookup_axis_for_impl(
            self,
            solve_for,
            expression=expression,
            target=target,
            length=length,
            gmid=gmid,
            vds=vds,
            vbs=vbs,
            mode=mode,
            out_of_range=out_of_range,
            warn_multiple=warn_multiple,
            method=method,
        )

    def lookup_gmid_for(
        self,
        *,
        length: float,
        vds: float,
        vbs: float,
        expression: Expression,
        target: float,
        mode: str = "first",
        out_of_range: str = "nan",
        warn_multiple: bool = True,
        method: str = "scan",
    ) -> float | list[float]:
        """Find the gmid value that gives expression(...) == target.

        Parameters are the same as lookup_axis_for with solve_for="gmid".
        """
        return self.lookup_axis_for(
            "gmid",
            length=length,
            vds=vds,
            vbs=vbs,
            expression=expression,
            target=target,
            mode=mode,
            out_of_range=out_of_range,
            warn_multiple=warn_multiple,
            method=method,
        )

    def lookup_length_for(
        self,
        *,
        gmid: float,
        vds: float,
        vbs: float,
        expression: Expression,
        target: float,
        mode: str = "first",
        out_of_range: str = "nan",
        warn_multiple: bool = True,
        method: str = "scan",
    ) -> float | list[float]:
        """Find the length value that gives expression(...) == target.

        Parameters are the same as lookup_axis_for with solve_for="length".
        """
        return self.lookup_axis_for(
            "length",
            gmid=gmid,
            vds=vds,
            vbs=vbs,
            expression=expression,
            target=target,
            mode=mode,
            out_of_range=out_of_range,
            warn_multiple=warn_multiple,
            method=method,
        )

    def lookup_vds_for(
        self,
        *,
        length: float,
        gmid: float,
        vbs: float,
        expression: Expression,
        target: float,
        mode: str = "first",
        out_of_range: str = "nan",
        warn_multiple: bool = True,
        method: str = "scan",
    ) -> float | list[float]:
        """Find the vds value that gives expression(...) == target.

        Parameters are the same as lookup_axis_for with solve_for="vds".
        """
        return self.lookup_axis_for(
            "vds",
            length=length,
            gmid=gmid,
            vbs=vbs,
            expression=expression,
            target=target,
            mode=mode,
            out_of_range=out_of_range,
            warn_multiple=warn_multiple,
            method=method,
        )

    def lookup_vbs_for(
        self,
        *,
        length: float,
        gmid: float,
        vds: float,
        expression: Expression,
        target: float,
        mode: str = "first",
        out_of_range: str = "nan",
        warn_multiple: bool = True,
        method: str = "scan",
    ) -> float | list[float]:
        """Find the vbs value that gives expression(...) == target.

        Parameters are the same as lookup_axis_for with solve_for="vbs".
        """
        return self.lookup_axis_for(
            "vbs",
            length=length,
            gmid=gmid,
            vds=vds,
            expression=expression,
            target=target,
            mode=mode,
            out_of_range=out_of_range,
            warn_multiple=warn_multiple,
            method=method,
        )
