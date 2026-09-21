"""Reusable topology model for reduced-MNA small-signal analysis."""

from __future__ import annotations

from dataclasses import dataclass
from functools import cached_property, lru_cache

import numpy as np

from ._types import CompiledPort, CompiledTransfer, GROUND, UNKNOWN_NODE
from .analysis import PortAnalysis, TransferAnalysis
from .kernels import assemble_compiled
from .system import LinearSystem


InputItems = tuple[tuple[str, float], ...]
OutputItems = tuple[tuple[str, float], ...]


def canonical_inputs(inputs: dict[str, float]) -> InputItems:
    if not isinstance(inputs, dict):
        raise TypeError("inputs must be a non-empty dict mapping node names to weights.")
    items = tuple(sorted((node, float(value)) for node, value in inputs.items()))
    if not items:
        raise ValueError("At least one AC input must be provided.")
    return items


def canonical_output(output: dict[str, float]) -> OutputItems:
    if not isinstance(output, dict):
        raise TypeError("output must be a non-empty dict mapping node names to weights.")
    items = tuple(sorted((node, float(weight)) for node, weight in output.items()))
    if not items:
        raise ValueError("At least one output node must be provided.")
    return items


def _input_value(node: str, inputs: InputItems) -> float:
    for input_node, value in inputs:
        if node == input_node:
            return value
    return 0.0


def _node_index_and_known(
    node: str,
    node_idx: dict[str, int],
    inputs: InputItems,
) -> tuple[int, float]:
    """Map a node name to either a matrix index or a known-source value."""

    idx = node_idx.get(node)
    if idx is not None:
        return idx, 0.0
    return UNKNOWN_NODE, _input_value(node, inputs)


def _compile_node_maps(
    node_idx: dict[str, int],
    node_groups: tuple[tuple[str, ...], ...],
    inputs: InputItems,
) -> tuple[np.ndarray, np.ndarray]:
    """Compile node-name tuples into index and known-source arrays.

    ``node_groups`` is either a tuple of MOS nodes ``(d, g, s, b)`` or a tuple
    of capacitor nodes ``(a, b)``. The returned arrays have matching shape:

    * ``idx``: matrix row/column indices, or ``UNKNOWN_NODE``
    * ``known``: known-node voltage for the requested AC input
    """

    width = len(node_groups[0]) if node_groups else 0
    idx = np.empty((len(node_groups), width), dtype=np.int64)
    known = np.zeros_like(idx, dtype=float)
    for i, nodes in enumerate(node_groups):
        for j, node in enumerate(nodes):
            idx[i, j], known[i, j] = _node_index_and_known(node, node_idx, inputs)
    return idx, known


@dataclass(frozen=True)
class SmallSignalModel:
    """Reusable small-signal topology with numeric values supplied per solve.

    This object captures only the topology: MOS terminal names and capacitor
    terminal names. It does not store gm/gds/c values. That separation lets
    optimizer code reuse the same compiled node map while changing only numeric
    device values at every candidate point.
    """

    mos_nodes: tuple[tuple[str, str, str, str], ...]
    cap_nodes: tuple[tuple[str, str], ...]

    @cached_property
    def nodes(self) -> tuple[str, ...]:
        """All nodes in first-seen order.

        Stable ordering matters for reproducible matrix layouts and profiles;
        using a set here would make node indices depend on hash iteration order.
        """

        node_order: dict[str, None] = {}
        for nodes in self.mos_nodes:
            for node in nodes:
                node_order.setdefault(node, None)
        for nodes in self.cap_nodes:
            for node in nodes:
                node_order.setdefault(node, None)
        return tuple(node_order)

    def _node_index(self, known_nodes: frozenset[str]) -> dict[str, int]:
        """Assign matrix indices to unknown non-ground nodes."""

        return {
            node: i
            for i, node in enumerate(
                node for node in self.nodes if node != GROUND and node not in known_nodes
            )
        }

    @lru_cache(maxsize=32)
    def _compile_transfer(
        self,
        inputs: InputItems,
        output: OutputItems,
    ) -> CompiledTransfer:
        """Compile topology for a specific known-node set and output observation.

        The same topology can be solved with different input vectors, so the
        compiled map is cached by input nodes and output observation. Typical
        optimizer use has exactly one such tuple, which means this work happens
        once.
        """

        known_nodes = frozenset(node for node, _ in inputs)
        node_idx = self._node_index(known_nodes)
        out_idx = np.empty(len(output), dtype=np.int64)
        weights = np.empty(len(output), dtype=float)
        for i, (node, weight) in enumerate(output):
            if node not in node_idx:
                raise ValueError(
                    f"Output node '{node}' is not an unknown node. Unknown nodes: {list(node_idx)}"
                )
            out_idx[i] = node_idx[node]
            weights[i] = weight

        mos_idx, mos_known = _compile_node_maps(node_idx, self.mos_nodes, inputs)
        cap_idx, cap_known = _compile_node_maps(node_idx, self.cap_nodes, inputs)
        return CompiledTransfer(
            node_idx=node_idx,
            mos_idx=mos_idx,
            mos_known=mos_known,
            cap_idx=cap_idx,
            cap_known=cap_known,
            out_idx=out_idx,
            weights=weights,
        )

    @lru_cache(maxsize=32)
    def _compile_port(self, node: str, reference: str) -> CompiledPort:
        """Compile topology for a one-port small-signal impedance query."""

        if node == reference:
            raise ValueError("Impedance probe node and reference must be different.")

        node_idx = self._node_index(frozenset())
        if node != GROUND and node not in node_idx:
            raise ValueError(f"Impedance probe node '{node}' is not in this topology.")
        if reference != GROUND and reference not in node_idx:
            raise ValueError(f"Impedance reference node '{reference}' is not in this topology.")

        empty_inputs: InputItems = ()
        mos_idx, mos_known = _compile_node_maps(node_idx, self.mos_nodes, empty_inputs)
        cap_idx, cap_known = _compile_node_maps(node_idx, self.cap_nodes, empty_inputs)
        return CompiledPort(
            node_idx=node_idx,
            mos_idx=mos_idx,
            mos_known=mos_known,
            cap_idx=cap_idx,
            cap_known=cap_known,
            node_idx_probe=node_idx.get(node, UNKNOWN_NODE),
            reference_idx=node_idx.get(reference, UNKNOWN_NODE),
        )

    @staticmethod
    def _assemble(
        mos_values: np.ndarray,
        cap_values: np.ndarray,
        compiled: CompiledTransfer | CompiledPort,
    ) -> tuple[LinearSystem, np.ndarray, np.ndarray]:
        G, C, rhs_g, rhs_c = assemble_compiled(
            mos_values,
            cap_values,
            compiled.mos_idx,
            compiled.mos_known,
            compiled.cap_idx,
            compiled.cap_known,
            len(compiled.node_idx),
        )
        system = LinearSystem(G=G, C=C, node_idx=compiled.node_idx)
        return system, rhs_g, rhs_c

    def transfer(
        self,
        mos_values: np.ndarray,
        cap_values: np.ndarray,
        *,
        inputs: dict[str, float],
        output: dict[str, float],
        gbw_iters: int = 32,
    ) -> TransferAnalysis:
        """Create a voltage-transfer analysis for one numeric operating point."""

        input_items = canonical_inputs(inputs)
        compiled = self._compile_transfer(input_items, canonical_output(output))
        system, rhs_g, rhs_c = self._assemble(mos_values, cap_values, compiled)
        return TransferAnalysis(
            system=system,
            out_idx=compiled.out_idx,
            weights=compiled.weights,
            rhs_g=rhs_g,
            rhs_c=rhs_c,
            gbw_iters=gbw_iters,
        )

    def port(
        self,
        mos_values: np.ndarray,
        cap_values: np.ndarray,
        *,
        node: str,
        reference: str = GROUND,
    ) -> PortAnalysis:
        """Create a one-port impedance analysis for one numeric operating point."""

        compiled = self._compile_port(node, reference)
        system, _, _ = self._assemble(mos_values, cap_values, compiled)
        return PortAnalysis(
            system=system,
            node_idx_probe=compiled.node_idx_probe,
            reference_idx=compiled.reference_idx,
        )
