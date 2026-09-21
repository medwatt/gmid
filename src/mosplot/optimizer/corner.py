"""A process / voltage / temperature corner.

A `Corner` owns both the device lookup tables and the operating conditions for one corner.
Process and temperature corners usually carry their own characterised ``.npz`` LUT.
Corners that reuse the same LUT but change only solver/testbench ``conditions`` share the
same prebuilt FastMosfet cache by default.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from functools import cached_property, lru_cache

from mosplot.table.load import load_lookup_table

from .device import DeviceTable


# One in-memory copy per LUT file and device, shared by every Corner that uses it: a config
# with voltage corners (tt, tt_14, tt_10, ...) reuses the same few LUTs many times.
@lru_cache(maxsize=None)
def _shared_lut(path):
    return load_lookup_table(path)


@lru_cache(maxsize=None)
def _shared_fast_mosfet(path, device, cache_dir, cache_tag):
    from .fast_mosfet import FastMosfet

    return FastMosfet(_shared_lut(path), device, cache_dir=cache_dir, cache_tag=cache_tag)


@dataclass
class Corner:
    """One PVT corner: a LUT, the device names in it, and the operating conditions.

    Attributes:
        name: Corner label (e.g. "tt", "ff", "tt_lowvdd"). This is only a report/key label;
            it is not part of the FastMosfet cache identity.
        lut_path: Path to the corner's .npz lookup table.
        nmos_name, pmos_name: Device keys inside the LUT.
        conditions: Operating conditions read by the circuit via cond[...] -- e.g. vdd,
            vout_dc, vin_cm, cout, temp. Voltage corners live entirely here.
        cache_dir: Directory for cached FastMosfet grids.
        cache_tag: Optional manual cache invalidation tag. Leave as None for the normal
            data-derived cache key. The LUT contents already distinguish process/temp
            tables; conditions and corner labels should not create duplicate caches.
    """

    name: str
    lut_path: str
    nmos_name: str
    pmos_name: str
    conditions: dict = field(default_factory=dict)
    cache_dir: str | None = "~/.cache/mosplot/fast_tables"
    cache_tag: str | None = None

    @cached_property
    def nmos(self) -> DeviceTable:
        fm = _shared_fast_mosfet(self.lut_path, self.nmos_name, self.cache_dir, self.cache_tag)
        return DeviceTable(fm, "n")

    @cached_property
    def pmos(self) -> DeviceTable:
        fm = _shared_fast_mosfet(self.lut_path, self.pmos_name, self.cache_dir, self.cache_tag)
        return DeviceTable(fm, "p")

    def cond(self, key, default=None):
        return self.conditions.get(key, default)
