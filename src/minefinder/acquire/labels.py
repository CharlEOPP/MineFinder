"""Known-mine labels: normalize USMIN, MRDS, and Colorado inactive mines; merge and dedupe."""
from __future__ import annotations

import geopandas as gpd
import pandas as pd
from shapely import STRtree

from minefinder.aoi import AOI, PROJECT_CRS

# Highest location quality first. MRDS locations are often coarse, so it is last.
SOURCE_PRIORITY: tuple[str, ...] = ("usmin", "co_inactive", "mrds")

NORMALIZED_COLUMNS = ["source", "source_id", "name", "feature_type", "geometry"]
MERGED_COLUMNS = ["label_id", "name", "feature_type", "sources", "n_sources", "geometry"]


def normalize(
    gdf: gpd.GeoDataFrame,
    source: str,
    id_col: str,
    name_col: str,
    type_col: str | None = None,
) -> gpd.GeoDataFrame:
    """Return one source in the common schema, as points in EPSG:26913.

    feature_type is lower-cased so "Adit" (USMIN) and "adit" (CGS) compare equal.
    """
    if source not in SOURCE_PRIORITY:
        raise ValueError(f"unknown source {source!r}; expected one of {SOURCE_PRIORITY}")
    if gdf.crs is None:
        raise ValueError(f"{source}: input has no CRS; refusing to guess")

    # Drop empties first, then missing: avoids GeoPandas' notna() ambiguity warning.
    gdf = gdf[~gdf.geometry.is_empty]
    gdf = gdf[gdf.geometry.notna()]
    gdf = gdf.to_crs(PROJECT_CRS)

    # Polygons/lines (e.g. USMIN dump outlines) -> a point guaranteed inside the shape.
    geom = gdf.geometry.where(gdf.geom_type == "Point", gdf.geometry.representative_point())

    return gpd.GeoDataFrame(
        {
            "source": source,
            "source_id": gdf[id_col].astype(str).to_numpy(),
            "name": gdf[name_col].to_numpy(),
            "feature_type": (
                gdf[type_col].astype("string").str.strip().str.lower().astype(object)
                .where(gdf[type_col].notna(), None).to_numpy()
                if type_col
                else None
            ),
        },
        geometry=geom.to_numpy(),
        crs=PROJECT_CRS,
    )


def merge_labels(
    frames: list[gpd.GeoDataFrame],
    aoi: AOI | None = None,
    dedupe_m: float = 50.0,
) -> gpd.GeoDataFrame:
    """Merge normalized sources, collapsing cross-source duplicates within dedupe_m.

    Greedy by priority: every USMIN record starts a cluster; each lower-priority
    record joins the nearest cluster within dedupe_m that doesn't already hold a
    record from its own source, otherwise it starts a new cluster. The cluster
    keeps its anchor's (best source's) location and name.
    """
    for f in frames:
        if f.crs is None or f.crs != PROJECT_CRS:
            raise ValueError(f"all label frames must be EPSG:26913, got {f.crs}; run normalize() first")

    allpts = pd.concat(frames, ignore_index=True)
    if aoi is not None:
        allpts = allpts[allpts.geometry.within(aoi.extent)]

    # Deterministic order regardless of input order: priority, then source_id.
    allpts = allpts.assign(_rank=allpts["source"].map(SOURCE_PRIORITY.index))
    allpts = allpts.sort_values(["_rank", "source_id"], kind="stable")

    clusters: list[dict] = []
    for row in allpts.itertuples(index=False):
        best = None
        if clusters:
            tree = STRtree([c["geometry"] for c in clusters])
            for i in tree.query(row.geometry, predicate="dwithin", distance=dedupe_m):
                c = clusters[i]
                if row.source in c["sources"]:
                    continue  # same-source points never merge
                d = c["geometry"].distance(row.geometry)
                if best is None or d < best[0]:
                    best = (d, c)
        if best is not None:
            best[1]["sources"].add(row.source)
        else:
            clusters.append(
                {
                    "label_id": f"{row.source}:{row.source_id}",
                    "name": row.name,
                    "feature_type": row.feature_type,
                    "sources": {row.source},
                    "geometry": row.geometry,
                }
            )

    out = pd.DataFrame(
        {
            "label_id": [c["label_id"] for c in clusters],
            "name": [c["name"] for c in clusters],
            "feature_type": [c["feature_type"] for c in clusters],
            "sources": [tuple(sorted(c["sources"])) for c in clusters],
            "n_sources": [len(c["sources"]) for c in clusters],
        }
    )
    return gpd.GeoDataFrame(out, geometry=[c["geometry"] for c in clusters], crs=PROJECT_CRS)[MERGED_COLUMNS]
