from __future__ import annotations

from mosplot.util.format import si, si_area
from .optimizer import Optimizer


def _format_spec(key: str, value: float | None) -> str:
    if key == "Area":
        return si_area(value)
    try:
        return si(value, unicode=True)
    except (ValueError, TypeError):
        return str(value)


def _format_target(key: str, spec) -> str:
    mode_sym = {"max": "≥", "min": "≤", "eq": "="}
    sym = mode_sym.get(spec.mode, "")
    val = si_area(spec.target) if key == "Area" else si(spec.target, unicode=True)
    return f"{sym}{val}"


class DesignReport:
    """Generates a human-readable report for an optimized circuit."""

    def __init__(self, optimizer: Optimizer) -> None:
        self.optimizer = optimizer
        self.objective: float | None = (
            optimizer.result.fun if optimizer.result is not None else None
        )

    def report(self) -> str:
        lines: list[str] = []
        opt = self.optimizer

        # --- Transistor dimensions (from the frozen reference design) ---
        dims = opt.frozen.get("dims") if opt.frozen else None
        lines.append("\nTransistor Details:")
        if dims:
            for device, item in dims.items():
                lines.append(f"  {device}:")
                lines.append(f"    Length:  {si(item.get('Length'), unicode=True)}m")
                lines.append(f"    Width:   {si(item.get('Width'), unicode=True)}m")
                lines.append(f"    Area:    {si_area(item.get('Area'))}")
                lines.append(f"    Current: {si(item.get('Current'), unicode=True)}A")
                gmid = item.get("GMID")
                if gmid is not None:
                    lines.append(f"    GmID:    {gmid:.2f}")
        else:
            lines.append("  Not available.")

        # --- Per-corner spec table ---
        lines.append("\nCorner Results:")
        corner_results = [r for r in (opt.corner_results or []) if r is not None and r.specs]
        binding = opt.binding or {}
        target_specs = opt.target_specs

        if not corner_results:
            lines.append("  Not available.")
        else:
            all_keys = list(corner_results[0].specs.keys())
            corner_names = [r.name for r in corner_results]

            # Build formatted cell values
            rows: dict[str, list[str]] = {}
            for key in all_keys:
                rows[key] = [_format_spec(key, r.specs.get(key)) for r in corner_results]

            # Column widths
            key_w = max(len(k) for k in all_keys)
            col_ws = [
                max(len(name), max(len(rows[k][i]) for k in all_keys))
                for i, name in enumerate(corner_names)
            ]
            tgt_w = max(
                (
                    max(len(_format_target(k, target_specs[k])) for k in target_specs)
                    if target_specs
                    else 0
                ),
                len("Target"),
            )
            bind_w = max(
                (max(len(v) for v in binding.values()) if binding else 0),
                len("Binding"),
            )

            def row_line(label: str, cells: list[str], target: str = "", bind: str = "") -> str:
                parts = [f"  {label:<{key_w}}"]
                for cell, w in zip(cells, col_ws):
                    parts.append(f" | {cell:<{w}}")
                parts.append(f" | {target:<{tgt_w}}")
                parts.append(f" | {bind:<{bind_w}}")
                return "".join(parts)

            sep = (
                "  "
                + "-" * key_w
                + "".join("-+-" + "-" * w for w in col_ws)
                + "-+-"
                + "-" * tgt_w
                + "-+-"
                + "-" * bind_w
            )

            lines.append(row_line("Spec", corner_names, "Target", "Binding"))
            lines.append(sep)

            for key in all_keys:
                cells = rows[key]
                if key in target_specs:
                    target_str = _format_target(key, target_specs[key])
                    bind_str = binding.get(key, "--")
                else:
                    target_str = ""
                    bind_str = ""
                lines.append(row_line(key, cells, target_str, bind_str))

        if self.objective is not None:
            lines.append(f"\nFinal cost: {self.objective:.6g}")

        return "\n".join(lines)
