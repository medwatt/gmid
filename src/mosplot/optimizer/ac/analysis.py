from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from ._types import UNKNOWN_NODE
from .kernels import find_unity_gain, phase_margin, solve_response
from .system import LinearSystem


@dataclass
class TransferAnalysis:
    """Voltage-transfer analysis for a weighted linear output observation."""

    system: LinearSystem
    out_idx: np.ndarray
    weights: np.ndarray
    rhs_g: np.ndarray
    rhs_c: np.ndarray
    gbw_iters: int = 32
    _v: np.ndarray | None = field(default=None, init=False, repr=False)
    _gain: float | None = field(default=None, init=False, repr=False)
    _ugf_hz: float | None = field(default=None, init=False, repr=False)

    def __add__(self, other: "TransferAnalysis") -> "TransferAnalysis":
        """Return a composed transfer that sums this output with another."""

        return _combine_transfers(self, other, 1.0)

    def __sub__(self, other: "TransferAnalysis") -> "TransferAnalysis":
        """Return a composed transfer that subtracts another output."""

        return _combine_transfers(self, other, -1.0)

    def __mul__(self, scale: float) -> "TransferAnalysis":
        """Return this transfer scaled by ``scale``."""

        return TransferAnalysis(
            system=self.system,
            out_idx=self.out_idx.copy(),
            weights=self.weights * float(scale),
            rhs_g=self.rhs_g,
            rhs_c=self.rhs_c,
            gbw_iters=self.gbw_iters,
        )

    def __rmul__(self, scale: float) -> "TransferAnalysis":
        """Return this transfer scaled by ``scale``."""

        return self.__mul__(scale)

    def __neg__(self) -> "TransferAnalysis":
        """Return this transfer with inverted sign."""

        return self * -1.0

    def gain(self) -> float:
        """Return the zero-frequency voltage gain of the output observation."""

        if self._gain is None:
            self._v = self.system.solve_dc(self.rhs_g)
            gain = 0.0
            for out_idx, weight in zip(self.out_idx, self.weights):
                gain += float(weight) * float(self._v[out_idx])
            self._gain = gain
        return self._gain

    def rejection(self, interference: "TransferAnalysis") -> float:
        """Return absolute rejection ratio against another transfer analysis."""

        signal = self.gain()
        interferer = interference.gain()
        return abs(signal / interferer) if abs(interferer) > 1e-30 else 1e12

    def response(self, frequency_hz: float) -> complex:
        """Return the output observation response at one frequency in Hz."""

        return solve_response(
            self.system.G,
            self.system.C,
            self.rhs_g,
            self.rhs_c,
            self.out_idx,
            self.weights,
            2.0 * np.pi * float(frequency_hz),
        )

    def ugf(self) -> float:
        """Return the first unity-gain frequency in Hz."""

        if self._ugf_hz is not None:
            return self._ugf_hz

        dc_gain = self.gain()
        if self.system.has_dynamic_terms(self.rhs_c):
            self._ugf_hz = float(
                find_unity_gain(
                    self.system.G,
                    self.system.C,
                    self.rhs_g,
                    self.rhs_c,
                    self.out_idx,
                    self.weights,
                    self.gbw_iters,
                )
            )
        else:
            self._ugf_hz = float("inf") if abs(dc_gain) >= 1.0 else 0.0
        return self._ugf_hz

    def phase_margin(self) -> float:
        """Return phase margin in degrees at the unity-gain frequency."""

        return float(
            phase_margin(
                self.system.G,
                self.system.C,
                self.rhs_g,
                self.rhs_c,
                self.out_idx,
                self.weights,
                self.ugf(),
                self.gain(),
            )
        )


def _ensure_compatible(left: TransferAnalysis, right: TransferAnalysis) -> None:
    if left.system.node_idx != right.system.node_idx:
        raise ValueError("Cannot compose transfers from different node maps.")
    if not np.array_equal(left.system.G, right.system.G):
        raise ValueError("Cannot compose transfers from different conductance matrices.")
    if not np.array_equal(left.system.C, right.system.C):
        raise ValueError("Cannot compose transfers from different capacitance matrices.")
    if not np.array_equal(left.rhs_g, right.rhs_g):
        raise ValueError("Cannot compose transfers with different input conductance RHS vectors.")
    if not np.array_equal(left.rhs_c, right.rhs_c):
        raise ValueError("Cannot compose transfers with different input capacitance RHS vectors.")


def _combine_transfers(
    left: TransferAnalysis,
    right: TransferAnalysis,
    right_scale: float,
) -> TransferAnalysis:
    _ensure_compatible(left, right)
    return TransferAnalysis(
        system=left.system,
        out_idx=np.concatenate((left.out_idx, right.out_idx)),
        weights=np.concatenate((left.weights, right.weights * right_scale)),
        rhs_g=left.rhs_g,
        rhs_c=left.rhs_c,
        gbw_iters=max(left.gbw_iters, right.gbw_iters),
    )


@dataclass
class PortAnalysis:
    """One-port impedance analysis."""

    system: LinearSystem
    node_idx_probe: int
    reference_idx: int

    def _current_rhs(self, current: float) -> np.ndarray:
        rhs = np.zeros(len(self.system.node_idx), dtype=float)
        if self.node_idx_probe != UNKNOWN_NODE:
            rhs[self.node_idx_probe] -= current
        if self.reference_idx != UNKNOWN_NODE:
            rhs[self.reference_idx] += current
        return rhs

    @staticmethod
    def _voltage_at(solution: np.ndarray, idx: int) -> complex:
        return 0.0 if idx == UNKNOWN_NODE else solution[idx]

    def impedance(self, frequency_hz: float = 0.0, *, test_current: float = 1.0) -> complex:
        """Return impedance seen by a current entering the probe node."""

        if test_current == 0.0:
            raise ValueError("test_current must be non-zero.")

        rhs_g = self._current_rhs(float(test_current))
        if frequency_hz == 0.0:
            v = self.system.solve_dc(rhs_g)
        else:
            rhs_c = np.zeros(len(self.system.node_idx), dtype=float)
            v = self.system.solve_ac(rhs_g, rhs_c, frequency_hz)

        node_v = self._voltage_at(v, self.node_idx_probe)
        reference_v = self._voltage_at(v, self.reference_idx)
        return (node_v - reference_v) / test_current

    def resistance(self, *, test_current: float = 1.0) -> float:
        """Return the DC small-signal resistance."""

        return float(np.real(self.impedance(0.0, test_current=test_current)))
