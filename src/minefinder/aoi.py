"""Area of interest: extent, CRS, and the held-out spatial test block."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Literal

from pyproj import CRS
from shapely.geometry.base import BaseGeometry
from shapely.geometry import Polygon

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
        raise NotImplementedError

    def split_of(self, geom: BaseGeometry) -> Split:
        """Which split a geometry (point label or chip footprint) belongs to."""
        raise NotImplementedError


def load_aoi(path: str | Path) -> AOI:
    """Load and validate an AOI YAML config."""
    raise NotImplementedError
