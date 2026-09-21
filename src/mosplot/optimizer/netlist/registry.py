from __future__ import annotations

from .base import NetlistGenerator
from .spectre import SpectreGenerator


_NETLIST_GENERATORS: dict[str, type[NetlistGenerator]] = {}


def _normalize_format_name(format_name: str) -> str:
    if not isinstance(format_name, str):
        raise TypeError("netlist format name must be a string.")
    key = format_name.strip().lower()
    if not key:
        raise ValueError("netlist format name must not be empty.")
    return key


def register_netlist_generator(
    format_name: str,
    generator_cls: type[NetlistGenerator],
) -> None:
    """Register ``generator_cls`` as the backend for ``format_name``."""

    if not isinstance(generator_cls, type) or not issubclass(generator_cls, NetlistGenerator):
        raise TypeError("generator_cls must be a NetlistGenerator subclass.")
    _NETLIST_GENERATORS[_normalize_format_name(format_name)] = generator_cls


def available_netlist_formats() -> tuple[str, ...]:
    """Return registered netlist format names."""

    return tuple(sorted(_NETLIST_GENERATORS))


def get_netlist_generator(format_name: str) -> type[NetlistGenerator]:
    """Return the generator class registered for ``format_name``."""

    key = _normalize_format_name(format_name)
    try:
        return _NETLIST_GENERATORS[key]
    except KeyError as exc:
        available = ", ".join(available_netlist_formats())
        raise ValueError(
            f"Unsupported netlist format '{format_name}'. Available formats: {available}."
        ) from exc


register_netlist_generator("spectre", SpectreGenerator)
