import geopandas as gpd
import pandas as pd
import pytest
from pyproj import Transformer
from shapely.geometry import Point, box

from minefinder.acquire.labels import MERGED_COLUMNS, NORMALIZED_COLUMNS, merge_labels, normalize
from tests.conftest import X0, Y0

UTM = "EPSG:26913"


def pts(source, coords, crs=UTM):
    """Normalized-schema frame of points (already in UTM unless crs says otherwise)."""
    return gpd.GeoDataFrame(
        {
            "source": source,
            "source_id": [f"{source}-{i}" for i in range(len(coords))],
            "name": [f"{source} mine {i}" for i in range(len(coords))],
        },
        geometry=[Point(c) for c in coords],
        crs=crs,
    )


def raw(ids, names, geoms, crs):
    return gpd.GeoDataFrame({"DEP_ID": ids, "SITE_NAME": names}, geometry=geoms, crs=crs)


# ---------- normalize ----------

def test_normalize_schema_and_crs():
    g = raw(["10"], ["Big Five"], [Point(X0, Y0)], UTM)
    out = normalize(g, "mrds", id_col="DEP_ID", name_col="SITE_NAME")
    assert list(out.columns) == NORMALIZED_COLUMNS
    assert out.crs.to_epsg() == 26913
    assert out.iloc[0]["source"] == "mrds"
    assert out.iloc[0]["source_id"] == "10"
    assert out.iloc[0]["name"] == "Big Five"


def test_normalize_reprojects_from_geographic():
    lon, lat = -105.5094, 40.0722  # town of Ward (a public place, not a candidate)
    x, y = Transformer.from_crs(4326, 26913, always_xy=True).transform(lon, lat)
    g = raw(["1"], ["x"], [Point(lon, lat)], "EPSG:4326")
    p = normalize(g, "mrds", "DEP_ID", "SITE_NAME").geometry.iloc[0]
    assert p.x == pytest.approx(x, abs=0.01)
    assert p.y == pytest.approx(y, abs=0.01)


def test_normalize_refuses_missing_crs():
    g = raw(["1"], ["x"], [Point(X0, Y0)], None)
    with pytest.raises(ValueError, match="CRS"):
        normalize(g, "mrds", "DEP_ID", "SITE_NAME")


def test_normalize_polygons_become_interior_points():
    # USMIN includes polygon features (e.g. dumps); labels are points.
    poly = box(X0, Y0, X0 + 40, Y0 + 20)
    out = normalize(raw(["7"], ["dump"], [poly], UTM), "usmin", "DEP_ID", "SITE_NAME")
    p = out.geometry.iloc[0]
    assert p.geom_type == "Point"
    assert poly.contains(p)


def test_normalize_drops_null_and_empty_geometries():
    g = raw(["1", "2", "3"], ["a", "b", "c"], [Point(X0, Y0), None, Point()], UTM)
    out = normalize(g, "mrds", "DEP_ID", "SITE_NAME")
    assert list(out["source_id"]) == ["1"]


def test_normalize_rejects_unknown_source():
    with pytest.raises(ValueError, match="source"):
        normalize(raw(["1"], ["a"], [Point(X0, Y0)], UTM), "wikipedia", "DEP_ID", "SITE_NAME")


# ---------- merge ----------

def test_merge_schema():
    out = merge_labels([pts("usmin", [(X0 + 5_000, Y0 + 5_000)])])
    assert list(out.columns) == MERGED_COLUMNS
    assert out.crs.to_epsg() == 26913


def test_cross_source_duplicates_collapse():
    usmin = pts("usmin", [(X0 + 5_000, Y0 + 5_000)])
    mrds = pts("mrds", [(X0 + 5_030, Y0 + 5_000)])  # 30 m away -> same mine
    out = merge_labels([usmin, mrds], dedupe_m=50)
    assert len(out) == 1
    assert out.iloc[0]["sources"] == ("mrds", "usmin")
    assert out.iloc[0]["n_sources"] == 2


def test_far_apart_points_stay_separate():
    usmin = pts("usmin", [(X0 + 5_000, Y0 + 5_000)])
    mrds = pts("mrds", [(X0 + 5_200, Y0 + 5_000)])
    assert len(merge_labels([usmin, mrds], dedupe_m=50)) == 2


def test_same_source_points_never_merge():
    # Two real shafts 20 m apart in one inventory are two features, not a duplicate.
    usmin = pts("usmin", [(X0 + 5_000, Y0 + 5_000), (X0 + 5_020, Y0 + 5_000)])
    assert len(merge_labels([usmin], dedupe_m=50)) == 2


def test_merged_location_comes_from_highest_priority_source():
    # MRDS locations are coarse; USMIN wins regardless of input order.
    mrds = pts("mrds", [(X0 + 5_040, Y0 + 5_000)])
    usmin = pts("usmin", [(X0 + 5_000, Y0 + 5_000)])
    p = merge_labels([mrds, usmin], dedupe_m=50).geometry.iloc[0]
    assert (p.x, p.y) == (X0 + 5_000, Y0 + 5_000)


def test_merged_name_comes_from_highest_priority_source():
    co = pts("co_inactive", [(X0 + 5_010, Y0 + 5_000)])
    mrds = pts("mrds", [(X0 + 5_000, Y0 + 5_000)])
    out = merge_labels([mrds, co], dedupe_m=50)
    assert out.iloc[0]["name"] == "co_inactive mine 0"


def test_three_way_match():
    out = merge_labels(
        [
            pts("usmin", [(X0 + 5_000, Y0 + 5_000)]),
            pts("co_inactive", [(X0 + 5_010, Y0 + 5_010)]),
            pts("mrds", [(X0 + 4_980, Y0 + 5_000)]),
        ],
        dedupe_m=50,
    )
    assert len(out) == 1
    assert out.iloc[0]["sources"] == ("co_inactive", "mrds", "usmin")


def test_merge_refuses_mixed_crs():
    utm = pts("usmin", [(X0 + 5_000, Y0 + 5_000)])
    geo = pts("mrds", [(-105.5, 40.07)], crs="EPSG:4326")
    with pytest.raises(ValueError, match="26913"):
        merge_labels([utm, geo])


def test_merge_clips_to_aoi_extent(aoi):
    inside = (X0 + 5_000, Y0 + 5_000)
    outside = (X0 - 1_000, Y0 + 5_000)
    out = merge_labels([pts("usmin", [inside, outside])], aoi=aoi)
    assert len(out) == 1


def test_label_ids_unique_and_deterministic():
    frames = [
        pts("usmin", [(X0 + 5_000, Y0 + 5_000), (X0 + 7_000, Y0 + 7_000)]),
        pts("mrds", [(X0 + 3_000, Y0 + 3_000)]),
    ]
    a = merge_labels(frames)
    b = merge_labels(list(reversed(frames)))
    assert a["label_id"].is_unique
    pd.testing.assert_frame_equal(
        a.sort_values("label_id").reset_index(drop=True),
        b.sort_values("label_id").reset_index(drop=True),
    )
