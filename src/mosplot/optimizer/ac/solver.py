"""Public small-signal AC solver API.

The implementation is split across a few small modules:

* :mod:`mosplot.optimizer.ac.model` stores reusable topology compilation.
* :mod:`mosplot.optimizer.ac.analysis` stores concrete analysis objects.
* :mod:`mosplot.optimizer.ac.kernels` stores hot numerical kernels.
* :mod:`mosplot.optimizer.ac._types` stores shared containers and constants.
"""

from __future__ import annotations

import numpy as np

from ._types import GROUND, CapStamp, MosStamp
from .analysis import PortAnalysis, TransferAnalysis
from .model import SmallSignalModel


class SmallSignalSolver:
    """Reduced-MNA small-signal solver for optimizer-grade AC metrics.

    Modeling assumptions
    --------------------
    * ``inputs`` maps ideal independent voltage-source nodes to their AC
      amplitudes. These nodes are removed from the matrix and moved to the RHS.
    * Ground, supply, and ideal bias-source nodes should be named ``"gnd"`` by
      the caller or topology builder.
    * MOS stamps use caller-provided polarity. ``gm`` and optional ``gmb`` are
      drain-to-source current sources controlled by ``vgs`` and ``vbs``.
    * Capacitances are two-terminal lumped capacitances. This is deliberately
      compact for optimizer screening; it is not a full charge-conserving MOS
      capacitance matrix.
    """

    GND = GROUND

    def __init__(
        self,
        model: SmallSignalModel | None = None,
        mos_values: np.ndarray | None = None,
        cap_values: np.ndarray | None = None,
    ) -> None:
        self._model = model
        self._mos_values = mos_values
        self._cap_values = cap_values
        self._mos: list[MosStamp] = []
        self._caps: list[CapStamp] = []

    @classmethod
    def from_model(
        cls,
        model: SmallSignalModel,
        mos_values: np.ndarray,
        cap_values: np.ndarray,
    ) -> SmallSignalSolver:
        """Create a solver from a reusable topology and current numeric values."""

        return cls(model=model, mos_values=mos_values, cap_values=cap_values)

    def reset(self) -> None:
        """Remove all stamps and any cached model reference from this instance."""

        self._model = None
        self._mos_values = None
        self._cap_values = None
        self._mos.clear()
        self._caps.clear()

    def add_mos(
        self,
        *,
        gm: float,
        gds: float,
        d: str,
        g: str,
        s: str,
        b: str,
        gmb: float = 0.0,
    ) -> None:
        """Add a linear MOS small-signal stamp.

        The stamp contributes:

        * ``gds`` between drain and source
        * ``gm * (Vg - Vs)`` flowing from drain to source
        * ``gmb * (Vb - Vs)`` flowing from drain to source

        Use ``"gnd"`` for AC-grounded supply, bias, or substrate nodes.
        """

        self._mos.append(
            MosStamp(
                gm=float(gm),
                gmb=float(gmb),
                gds=float(gds),
                d=d,
                g=g,
                s=s,
                b=b,
            )
        )

    def add_cap(self, *, c: float, a: str, b: str) -> None:
        """Add a two-terminal linear capacitance between nodes ``a`` and ``b``."""

        self._caps.append(CapStamp(c=float(c), a=a, b=b))

    def _as_model_and_values(self) -> tuple[SmallSignalModel, np.ndarray, np.ndarray]:
        """Return a reusable topology model plus numeric value arrays.

        Direct ``add_mos``/``add_cap`` use is converted to the same model path
        used by the topology builder, so there is only one production solve
        implementation to maintain.
        """

        if self._model is not None:
            if self._mos_values is None or self._cap_values is None:
                raise RuntimeError("Compiled small-signal model is missing numeric values.")
            return self._model, self._mos_values, self._cap_values

        if not self._mos and not self._caps:
            raise RuntimeError("No stamps added; call add_mos/add_cap first.")

        model = SmallSignalModel(
            mos_nodes=tuple((m.d, m.g, m.s, m.b) for m in self._mos),
            cap_nodes=tuple((c.a, c.b) for c in self._caps),
        )
        mos_values = np.array([(m.gm, m.gmb, m.gds) for m in self._mos], dtype=float)
        cap_values = np.array([c.c for c in self._caps], dtype=float)
        return model, mos_values, cap_values

    def transfer(
        self,
        *,
        inputs: dict[str, float],
        output: dict[str, float],
        gbw_iters: int = 32,
    ) -> TransferAnalysis:
        """Create a reusable voltage-transfer analysis."""

        model, mos_values, cap_values = self._as_model_and_values()
        return model.transfer(
            mos_values,
            cap_values,
            inputs=inputs,
            output=output,
            gbw_iters=gbw_iters,
        )

    def port(
        self,
        node: str,
        *,
        reference: str = GROUND,
    ) -> PortAnalysis:
        """Create a reusable one-port impedance analysis."""

        model, mos_values, cap_values = self._as_model_and_values()
        return model.port(
            mos_values,
            cap_values,
            node=node,
            reference=reference,
        )
