"""Loader tests. Each source gets a tiny synthetic file in its real on-disk format
(USMIN shapefiles in WGS84, MRDS national CSV, CGS geodatabase-style layers), so
the tests pin down the schema quirks we found in the real downloads."""
from pathlib import Path

import geopandas as gpd
import pandas as pd
import pytest
from pyproj import Transformer
from shapely.geometry import Point, box

from minefinder.acquire.load_labels import (
    CGS_LAYERS,
    build_labels,
    load_co_aml,
    load_mrds,
    load_usmin,
)
from tests.conftest import REPO, X0, Y0

TO_LL = Transformer.from_crs(26913, 4326, always_xy=True)


def ll(x, y):
    """UTM 13N -> lon/lat, so fixtures can be written in metres and stored as WGS84."""
    return TO_LL.transform(x, y)


# ---------------- USMIN ----------------

def write_usmin(d: Path, points, polys=()):
    """points/polys: list of (x, y or geom, FTR_TYPE, FTR_NAME, TOPO_SCALE)."""
    d.mkdir(parents=True, exist_ok=True)
    cols = ["FTR_TYPE", "FTR_NAME", "TOPO_NAME", "TOPO_SCALE"]
    pt = gpd.GeoDataFrame(
        [{"FTR_TYPE": t, "FTR_NAME": n, "TOPO_NAME": "Ward", "TOPO_SCALE": s} for _, _, t, n, s in points],
        columns=cols,
        geometry=[Point(ll(x, y)) for x, y, *_ in points],
        crs="EPSG:4326",
    )
    pt.to_file(d / "CO-point.shp")
    if polys:
        def to_ll(g):
            import numpy as np
            import shapely
            return shapely.transform(g, lambda c: np.c_[TO_LL.transform(c[:, 0], c[:, 1])])
        po = gpd.GeoDataFrame(
            [{"FTR_TYPE": t, "FTR_NAME": n, "TOPO_NAME": "Ward", "TOPO_SCALE": s} for _, t, n, s in polys],
            columns=cols,
            geometry=[to_ll(g) for g, *_ in polys],
            crs="EPSG:4326",
        )
        po.to_file(d / "CO-poly.shp")


def test_usmin_reads_points_and_polygons_into_26913(tmp_path):
    write_usmin(
        tmp_path / "usmin",
        [(X0 + 5_000, Y0 + 5_000, "Adit", "Big Five", "24000")],
        [(box(X0 + 6_000, Y0 + 6_000, X0 + 6_040, Y0 + 6_020), "Tailings - Pond", None, "24000")],
    )
    out = load_usmin(tmp_path / "usmin")
    assert out.crs.to_epsg() == 26913
    assert len(out) == 2
    assert set(out.geom_type) == {"Point"}
    assert set(out["feature_type"]) == {"adit", "tailings - pond"}
    assert (out["source"] == "usmin").all()


def test_usmin_drops_aggregate_pits(tmp_path):
    # Gravel/borrow/sand pits are road-building aggregate, not the hard-rock landforms we target.
    write_usmin(
        tmp_path / "usmin",
        [
            (X0 + 5_000, Y0 + 5_000, "Adit", None, "24000"),
            (X0 + 6_000, Y0 + 5_000, "Gravel Pit", None, "24000"),
            (X0 + 7_000, Y0 + 5_000, "Borrow Pit", None, "24000"),
            (X0 + 8_000, Y0 + 5_000, "Sand Pit", None, "24000"),
        ],
    )
    out = load_usmin(tmp_path / "usmin")
    assert list(out["feature_type"]) == ["adit"]


def test_usmin_drops_62k_duplicates_of_24k_features(tmp_path):
    # ver 10.0 has no reconciliation between map scales: the same shaft drawn on the
    # 1:24k and 1:62.5k quads appears twice, ~10 m apart. Keep the 24k one.
    write_usmin(
        tmp_path / "usmin",
        [
            (X0 + 5_000, Y0 + 5_000, "Mine Shaft", "A", "24000"),
            (X0 + 5_012, Y0 + 5_000, "Mine Shaft", "A", "62500"),  # duplicate -> dropped
            (X0 + 7_000, Y0 + 5_000, "Mine Shaft", "B", "62500"),  # only on 62.5k -> kept
        ],
    )
    out = load_usmin(tmp_path / "usmin", scale_dedupe_m=50)
    assert sorted(out["name"]) == ["A", "B"]
    a = out[out["name"] == "A"].geometry.iloc[0]
    assert a.x == pytest.approx(X0 + 5_000, abs=0.01)


def test_usmin_ids_are_unique_and_stable_across_reads(tmp_path):
    # USMIN has no feature ID column, so we derive one. It must not change between runs.
    write_usmin(
        tmp_path / "usmin",
        [
            (X0 + 5_000, Y0 + 5_000, "Adit", None, "24000"),
            (X0 + 5_100, Y0 + 5_000, "Adit", None, "24000"),
        ],
    )
    a = load_usmin(tmp_path / "usmin")
    b = load_usmin(tmp_path / "usmin")
    assert a["source_id"].is_unique
    assert list(a["source_id"]) == list(b["source_id"])


def test_usmin_ids_do_not_contain_coordinates(tmp_path):
    # IDs flow into review exports; they must not leak locations.
    write_usmin(tmp_path / "usmin", [(X0 + 5_000, Y0 + 5_000, "Adit", None, "24000")])
    sid = load_usmin(tmp_path / "usmin")["source_id"].iloc[0]
    assert str(X0 + 5_000) not in sid and str(Y0 + 5_000) not in sid


# ---------------- MRDS ----------------

def write_mrds(path: Path, rows):
    """rows: (dep_id, name, x, y, state, dev_stat), x/y in UTM metres (None for missing)."""
    path.parent.mkdir(parents=True, exist_ok=True)
    recs = []
    for dep_id, name, x, y, state, dev in rows:
        lon, lat = ll(x, y) if x is not None else (None, None)
        recs.append(dict(dep_id=dep_id, site_name=name, latitude=lat, longitude=lon, state=state, dev_stat=dev))
    pd.DataFrame(recs).to_csv(path, index=False)


def test_mrds_reads_colorado_only_as_wgs84(tmp_path):
    p = tmp_path / "mrds" / "mrds.csv"
    write_mrds(p, [
        (10, "Caribou", X0 + 5_000, Y0 + 5_000, "Colorado", "Past Producer"),
        (11, "Elsewhere", X0 + 5_000, Y0 + 5_000, "Utah", "Producer"),
    ])
    out = load_mrds(p)
    assert list(out["source_id"]) == ["10"]
    pt = out.geometry.iloc[0]
    # Round-trip through WGS84 lands back within centimetres.
    assert pt.x == pytest.approx(X0 + 5_000, abs=0.05)
    assert out.iloc[0]["feature_type"] == "past producer"


def test_mrds_drops_rows_without_coordinates(tmp_path):
    p = tmp_path / "mrds" / "mrds.csv"
    write_mrds(p, [
        (10, "Has coords", X0 + 5_000, Y0 + 5_000, "Colorado", "Occurrence"),
        (11, "No coords", None, None, "Colorado", "Occurrence"),
    ])
    assert list(load_mrds(p)["source_id"]) == ["10"]


def test_mrds_collapses_stacked_records_at_identical_coordinates(tmp_path):
    # MRDS often geocodes several mines to the same section/claim centroid. Stacked
    # copies are one coarse location, not several mines; keep the lowest dep_id.
    p = tmp_path / "mrds" / "mrds.csv"
    write_mrds(p, [
        (30, "C", X0 + 5_000, Y0 + 5_000, "Colorado", "Occurrence"),
        (20, "B", X0 + 5_000, Y0 + 5_000, "Colorado", "Occurrence"),
        (40, "D", X0 + 9_000, Y0 + 5_000, "Colorado", "Occurrence"),
    ])
    out = load_mrds(p)
    assert sorted(out["source_id"]) == ["20", "40"]


# ---------------- CGS USFS AML inventory ----------------

def write_cgs(path: Path, openings, dumps, areas):
    """openings/dumps: (x, y, unique_id, type). areas: (unique_id, sitename)."""
    hole = gpd.GeoDataFrame(
        {"Unique_ID": [u for *_, u, _ in openings], "HTYPE": [t for *_, t in openings]},
        geometry=[Point(x, y) for x, y, *_ in openings], crs="EPSG:26913",
    )
    pile = gpd.GeoDataFrame(
        {"Unique_ID": [u for *_, u, _ in dumps], "PTYPE": [t for *_, t in dumps]},
        geometry=[Point(x, y) for x, y, *_ in dumps], crs="EPSG:26913",
    )
    hdr = gpd.GeoDataFrame(
        {"UniqueID": [u for u, _ in areas], "SITENAME_1": [n for _, n in areas]},
        geometry=[box(X0, Y0, X0 + 10, Y0 + 10)] * len(areas), crs="EPSG:26913",
    )
    hole.to_file(path, layer=CGS_LAYERS["openings"], driver="GPKG")
    pile.to_file(path, layer=CGS_LAYERS["dumps"], driver="GPKG")
    hdr.to_file(path, layer=CGS_LAYERS["areas"], driver="GPKG")


def test_cgs_reads_openings_and_dumps_with_site_names(tmp_path):
    p = tmp_path / "cgs.gpkg"
    write_cgs(
        p,
        openings=[(X0 + 5_000, Y0 + 5_000, "2810146744401100", "adit")],
        dumps=[(X0 + 5_030, Y0 + 4_980, "2810146744401200", "mine dump")],
        areas=[("2810146744401", "WEST MAGNOLIA RIDGE")],
    )
    out = load_co_aml(p)
    assert out.crs.to_epsg() == 26913
    assert (out["source"] == "co_inactive").all()
    assert sorted(out["feature_type"]) == ["adit", "mine dump"]
    # Feature ID = area ID + 3-digit feature number; the area carries the site name.
    assert set(out["name"]) == {"WEST MAGNOLIA RIDGE"}


def test_cgs_drops_dumps_recorded_at_their_openings_location(tmp_path):
    # ~60% of CGS dump points sit exactly on an opening (same GPS fix reused).
    # That is one site recorded twice, not two landforms.
    p = tmp_path / "cgs.gpkg"
    write_cgs(
        p,
        openings=[(X0 + 5_000, Y0 + 5_000, "2810146744401100", "adit")],
        dumps=[
            (X0 + 5_000, Y0 + 5_000, "2810146744401200", "mine dump"),   # coincident -> dropped
            (X0 + 5_040, Y0 + 4_970, "2810146744401201", "mine dump"),   # its own fix -> kept
        ],
        areas=[("2810146744401", "SITE")],
    )
    out = load_co_aml(p)
    assert sorted(out["source_id"]) == ["2810146744401100", "2810146744401201"]


# ---------------- end to end ----------------

def test_build_labels_writes_merged_gpkg_inside_aoi(tmp_path, aoi):
    raw = tmp_path / "raw"
    write_usmin(raw / "usmin", [
        (X0 + 5_000, Y0 + 5_000, "Adit", "Big Five", "24000"),
        (X0 - 3_000, Y0 + 5_000, "Adit", "Outside AOI", "24000"),
    ])
    write_mrds(raw / "mrds" / "mrds.csv", [(10, "Big Five Mine", X0 + 5_020, Y0 + 5_000, "Colorado", "Past Producer")])
    (raw / "co_aml").mkdir()
    write_cgs(raw / "co_aml" / "amli.gpkg",
              openings=[(X0 + 7_000, Y0 + 7_000, "2810146744401100", "vertical shaft")],
              dumps=[], areas=[("2810146744401", "SITE")])

    out_path = tmp_path / "labels.gpkg"
    merged = build_labels(raw, aoi, out_path, dedupe_m=50)

    assert out_path.exists()
    back = gpd.read_file(out_path)
    assert back.crs.to_epsg() == 26913
    assert len(merged) == len(back) == 2            # Big Five (usmin+mrds) and the CGS shaft
    big5 = merged[merged["name"] == "Big Five"].iloc[0]
    assert big5["sources"] == ("mrds", "usmin")
    assert merged.geometry.within(aoi.extent).all()
    # Split is tagged at build time so no later step can skip it.
    assert set(back["split"]) <= {"train", "test", "buffer"}
    assert big5["split"] == aoi.split_of(big5.geometry)


# ---------------- real data smoke test (skipped when data/raw isn't present) ----------------

RAW = REPO / "data" / "raw"


@pytest.mark.skipif(not (RAW / "usmin").exists(), reason="real downloads not present")
def test_real_data_builds(tmp_path):
    from minefinder.aoi import load_aoi

    aoi = load_aoi(REPO / "config" / "aoi_ward_goldhill.yaml")
    merged = build_labels(RAW, aoi, tmp_path / "labels.gpkg")
    assert len(merged) > 500
    assert merged.geometry.within(aoi.extent).all()
    for s in ("usmin", "co_inactive", "mrds"):
        assert merged["sources"].map(lambda t: s in t).any(), f"no labels from {s}"
