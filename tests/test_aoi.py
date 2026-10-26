import copy

import pytest
from shapely.geometry import Point, box

from minefinder.aoi import PROJECT_CRS, load_aoi
from tests.conftest import AOI_DICT, REPO, X0, Y0, write_yaml


# ---------- loading & validation ----------

def test_loads_name_crs_and_geometries(aoi):
    assert aoi.name == "test_aoi"
    assert aoi.crs == PROJECT_CRS
    assert aoi.extent.bounds == (X0, Y0, X0 + 10_000, Y0 + 10_000)
    assert aoi.holdout.bounds == (X0, Y0, X0 + 2_000, Y0 + 2_000)
    assert aoi.holdout_buffer_m == 100


def test_real_config_loads_and_is_valid():
    aoi = load_aoi(REPO / "config" / "aoi_ward_goldhill.yaml")
    assert aoi.crs.to_epsg() == 26913
    assert aoi.holdout.within(aoi.extent)


@pytest.mark.parametrize("bad_crs", ["EPSG:4326", "EPSG:32613", "EPSG:26913x"])
def test_rejects_any_crs_other_than_26913(tmp_path, bad_crs):
    # EPSG:32613 (WGS84 / UTM 13N) is the sneaky one: same zone, different datum.
    d = copy.deepcopy(AOI_DICT)
    d["crs"] = bad_crs
    with pytest.raises(ValueError, match="26913"):
        load_aoi(write_yaml(tmp_path, d))


def test_rejects_holdout_outside_extent(tmp_path):
    d = copy.deepcopy(AOI_DICT)
    d["holdout"] = {"minx": X0 - 500, "miny": Y0, "maxx": X0 + 1_000, "maxy": Y0 + 1_000}
    with pytest.raises(ValueError, match="holdout"):
        load_aoi(write_yaml(tmp_path, d))


def test_rejects_inverted_bounds(tmp_path):
    d = copy.deepcopy(AOI_DICT)
    d["extent"] = {"minx": X0 + 10_000, "miny": Y0, "maxx": X0, "maxy": Y0 + 10_000}
    with pytest.raises(ValueError):
        load_aoi(write_yaml(tmp_path, d))


def test_rejects_negative_buffer(tmp_path):
    d = copy.deepcopy(AOI_DICT)
    d["holdout_buffer_m"] = -1
    with pytest.raises(ValueError, match="buffer"):
        load_aoi(write_yaml(tmp_path, d))


# ---------- train region ----------

def test_train_region_excludes_holdout_and_buffer(aoi):
    tr = aoi.train_region
    assert not tr.intersects(aoi.holdout)
    # A point 50 m outside the holdout edge sits inside the 100 m buffer -> not train.
    assert not tr.contains(Point(X0 + 2_050, Y0 + 1_000))
    # 150 m outside the edge is clear of the buffer -> train.
    assert tr.contains(Point(X0 + 2_150, Y0 + 1_000))


def test_train_region_area(aoi):
    # extent minus (holdout grown by 100 m, clipped to extent, rounded corner)
    full = 10_000 ** 2
    grown = aoi.holdout.buffer(100).intersection(aoi.extent).area
    assert aoi.train_region.area == pytest.approx(full - grown)


# ---------- split assignment (the anti-leakage rule) ----------

@pytest.mark.parametrize(
    "geom, expected",
    [
        (Point(X0 + 1_000, Y0 + 1_000), "test"),        # inside holdout
        (Point(X0 + 2_050, Y0 + 1_000), "buffer"),      # in the gap
        (Point(X0 + 5_000, Y0 + 5_000), "train"),       # well clear
        (Point(X0 - 10, Y0 + 5_000), "outside"),        # off the AOI
        # A chip footprint is "test" only if entirely inside the holdout...
        (box(X0 + 100, Y0 + 100, X0 + 356, Y0 + 356), "test"),
        # ...and a chip straddling the holdout edge is discarded, never train.
        (box(X0 + 1_900, Y0 + 500, X0 + 2_156, Y0 + 756), "buffer"),
        # A chip touching the buffer from outside is also discarded.
        (box(X0 + 2_050, Y0 + 500, X0 + 2_306, Y0 + 756), "buffer"),
        (box(X0 + 4_000, Y0 + 4_000, X0 + 4_256, Y0 + 4_256), "train"),
    ],
)
def test_split_of(aoi, geom, expected):
    assert aoi.split_of(geom) == expected
