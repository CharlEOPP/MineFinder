"""Known-mine labels: normalize USMIN, MRDS, and Colorado inactive mines; merge and dedupe."""
from __future__ import annotations

import geopandas as gpd

from minefinder.aoi import AOI

# Highest location quality first. MRDS locations are often coarse, so it is last.
SOURCE_PRIORITY: tuple[str, ...] = ("usmin", "co_inactive", "mrds")

NORMALIZED_COLUMNS = ["source", "source_id", "name", "geometry"]
MERGED_COLUMNS = ["label_id", "name", "sources", "n_sources", "geometry"]


def normalize(gdf: gpd.GeoDataFrame, source: str, id_col: str, name_col: str) -> gpd.GeoDataFrame:
    """Return one source in the common schema, as points in EPSG:26913."""
    raise NotImplementedError


def merge_labels(
    frames: list[gpd.GeoDataFrame],
    aoi: AOI | None = None,
    dedupe_m: float = 50.0,
) -> gpd.GeoDataFrame:
    """Merge normalized sources, collapsing cross-source duplicates within dedupe_m."""
    raise NotImplementedError
