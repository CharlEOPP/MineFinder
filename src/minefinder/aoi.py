"""Area of interest: extent, CRS, and the held-out spatial test block."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Literal

import yaml
from pyproj import CRS
from shapely.geometry import Polygon, box
from shapely.geometry.base import BaseGeometry

PROJECT_CRS = CRS.from_epsg(26913)

Split = Literal["train", "test", "buffer", "outside"]


@dataclass(frozen=True)
class AOI:
    name: str
    crs: CRS
    extent: Polygon
    holdout: Polygon
    holdout_buffer_m: float

    @property
    def train_region(self) -> BaseGeometry:
        """Extent minus the holdout block and its buffer."""
        return self.extent.difference(self.holdout.buffer(self.holdout_buffer_m))

    def split_of(self, geom: BaseGeometry) -> Split:
        """Which split a geometry (point label or chip footprint) belongs to.

        Order matters: anything not wholly in test or wholly in train is
        'buffer' and gets discarded, so no chip can straddle the boundary.
        """
        if not geom.within(self.extent):
            return "outside"
        if geom.within(self.holdout):
            return "test"
        if geom.within(self.train_region):
            return "train"
        return "buffer"


def _box(d: dict, what: str) -> Polygon:
    try:
        minx, miny, maxx, maxy = (float(d[k]) for k in ("minx", "miny", "maxx", "maxy"))
    except (KeyError, TypeError) as e:
        raise ValueError(f"{what}: needs minx, miny, maxx, maxy") from e
    if minx >= maxx or miny >= maxy:
        raise ValueError(f"{what}: min must be less than max, got {d}")
    return box(minx, miny, maxx, maxy)


def load_aoi(path: str | Path) -> AOI:
    """Load and validate an AOI YAML config."""
    cfg = yaml.safe_load(Path(path).read_text())

    try:
        crs = CRS.from_user_input(cfg["crs"])
    except Exception as e:
        raise ValueError(f"crs {cfg.get('crs')!r} is not valid; must be EPSG:26913") from e
    if crs != PROJECT_CRS:
        raise ValueError(f"crs must be EPSG:26913, got {cfg['crs']}")

    extent = _box(cfg["extent"], "extent")
    holdout = _box(cfg["holdout"], "holdout")
    if not holdout.within(extent):
        raise ValueError("holdout block must lie inside the extent")

    buffer_m = float(cfg["holdout_buffer_m"])
    if buffer_m < 0:
        raise ValueError(f"holdout_buffer_m must be >= 0, got {buffer_m}")

    return AOI(cfg["name"], crs, extent, holdout, buffer_m)