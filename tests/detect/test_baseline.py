"""Rule-based baseline. A synthetic landscape with known pits and dumps, plus the two
classic look-alikes: a road embankment (long, linear high) and a ditch (long, linear low)."""
import geopandas as gpd
import numpy as np
import pandas as pd
import pytest
import rasterio
import yaml
from rasterio.transform import from_origin
from shapely.geometry import Point

from minefinder.detect.baseline import (
    DEFAULT_PARAMS,
    extract_from_arrays,
    extract_candidates,
    score_candidates,
    tune,
    evaluate_split,
)
from minefinder.terrain.derivatives import local_relief, LRM_SMALL, LRM_LARGE

X0, Y0 = 450_000, 4_430_000
N = 400

PITS = [(60, 80), (200, 150), (320, 300)]       # (col, row)
DUMPS = [(120, 300), (300, 90)]


def landscape(seed=0):
    rng = np.random.default_rng(seed)
    yy, xx = np.mgrid[0:N, 0:N].astype("float64")
    z = 2800 + 0.2 * xx - 0.1 * yy + rng.normal(0, 0.03, (N, N))
    for c, r in PITS:                            # 1.2 m deep, ~3 m radius bowls
        z -= 1.2 * np.exp(-((xx - c) ** 2 + (yy - r) ** 2) / (2 * 2.0**2))
    for c, r in DUMPS:                           # 3 m high, ~12 m radius mounds
        z += 3.0 * np.exp(-((xx - c) ** 2 + (yy - r) ** 2) / (2 * 8.0**2))
    z += 1.5 * np.exp(-((yy - 220) ** 2) / (2 * 3.0**2)) * ((xx > 20) & (xx < 380))   # road fill, E-W
    z -= 0.8 * np.exp(-((xx - 30) ** 2) / (2 * 1.5**2)) * ((yy > 20) & (yy < 380))    # ditch, N-S
    return z


@pytest.fixture(scope="module")
def arrays():
    z = landscape()
    return local_relief(z, LRM_SMALL), local_relief(z, LRM_LARGE), np.zeros((N, N), "uint8")


T = from_origin(X0, Y0, 1, 1)


def at(c, r):
    return X0 + c + 0.5, Y0 - r - 0.5


def near(cands, c, r, d=3.0):
    x, y = at(c, r)
    return cands[np.hypot(cands.geometry.x - x, cands.geometry.y - y) <= d]


# ---------------- extraction ----------------

def test_finds_every_pit_at_its_deepest_point(arrays):
    cands = extract_from_arrays(*arrays, T)
    pits = cands[cands["kind"] == "pit"]
    for c, r in PITS:
        assert len(near(pits, c, r, d=1.5)) == 1, (c, r)


def test_finds_every_dump(arrays):
    cands = extract_from_arrays(*arrays, T)
    dumps = cands[cands["kind"] == "dump"]
    for c, r in DUMPS:
        assert len(near(dumps, c, r, d=4.0)) == 1, (c, r)


def test_road_fill_and_ditch_are_rejected_as_too_elongated(arrays):
    scored = score_candidates(extract_from_arrays(*arrays, T), DEFAULT_PARAMS)
    x_road = scored.geometry.y.between(Y0 - 226, Y0 - 214) & (scored["kind"] == "dump")
    x_ditch = scored.geometry.x.between(X0 + 25, X0 + 36) & (scored["kind"] == "pit")
    assert not x_road.any() and not x_ditch.any()


def test_candidates_in_seam_mask_are_dropped(arrays):
    s, l, seam = arrays
    seam = seam.copy()
    c, r = PITS[0]
    seam[r - 5 : r + 6, c - 5 : c + 6] = 1
    cands = extract_from_arrays(s, l, seam, T)
    assert len(near(cands, c, r, d=3)) == 0


def test_candidate_table_schema(arrays):
    cands = extract_from_arrays(*arrays, T)
    assert {"kind", "prominence_m", "area_m2", "elongation"} <= set(cands.columns)
    assert cands.crs.to_epsg() == 26913


# ---------------- blockwise == whole ----------------

def write(path, a, dtype, nodata=None, scale=None):
    with rasterio.open(path, "w", driver="GTiff", width=N, height=N, count=1, dtype=dtype,
                       crs="EPSG:26913", transform=T, nodata=nodata) as d:
        if scale:
            d.scales = (scale,)
        d.write(a, 1)
    return path


def test_blockwise_extraction_equals_whole(tmp_path, arrays):
    s, l, seam = arrays
    ps = write(tmp_path / "lrm_small.tif", np.round(s / 0.01).astype("int16"), "int16", -32768, 0.01)
    pl = write(tmp_path / "lrm_large.tif", np.round(l / 0.01).astype("int16"), "int16", -32768, 0.01)
    pm = write(tmp_path / "seam_mask.tif", seam, "uint8")
    # Features longer than the halo can't be seen whole by one block, so both paths drop
    # anything wider than max_extent_m (roads, ridges: not mine features anyway).
    whole = extract_from_arrays(np.round(s / 0.01) * 0.01, np.round(l / 0.01) * 0.01, seam, T, max_extent_m=40)
    blocks = extract_candidates(tmp_path, block=96, halo=48, max_extent_m=40)
    assert len(whole) >= len(PITS) + len(DUMPS)
    key = lambda g: sorted(zip(g["kind"], g.geometry.x.round(1), g.geometry.y.round(1)))
    assert key(blocks) == key(whole)


# ---------------- scoring + tuning ----------------

def test_scoring_applies_thresholds_and_ranks(arrays):
    cands = extract_from_arrays(*arrays, T)
    strict = score_candidates(cands, {**DEFAULT_PARAMS, "pit_min_depth": 5.0})
    assert (strict["kind"] != "pit").all()          # no pit is 5 m deep
    s = score_candidates(cands, DEFAULT_PARAMS)
    assert s["score"].is_monotonic_decreasing


def labels_gdf(points, split):
    return gpd.GeoDataFrame({"sources": ["usmin"] * len(points), "split": [split] * len(points)},
                            geometry=[Point(*at(c, r)) for c, r in points], crs="EPSG:26913")


def test_tuning_never_looks_at_test_labels(arrays):
    # Same candidates, same train labels; wildly different test labels must not change
    # the chosen parameters. This is the anti-leakage rule from decision 0003.
    cands = extract_from_arrays(*arrays, T).assign(split="train")
    train = labels_gdf(PITS + DUMPS, "train")
    test_a = labels_gdf([(5, 5)], "test")
    test_b = labels_gdf(PITS, "test")
    grid = {"pit_min_depth": [0.3, 0.8], "dump_min_height": [0.5, 1.5]}
    pa, _ = tune(cands, pd.concat([train, test_a]), grid=grid, n=5)
    pb, _ = tune(cands, pd.concat([train, test_b]), grid=grid, n=5)
    assert pa == pb


def test_evaluate_test_split_requires_final_flag(arrays, tmp_path):
    cands = extract_from_arrays(*arrays, T).assign(split="test")
    labels = labels_gdf(PITS, "test")
    with pytest.raises(PermissionError, match="final"):
        evaluate_split(cands, labels, DEFAULT_PARAMS, region=None, split="test", final=False,
                       log_path=tmp_path / "log.csv")


def test_final_test_evaluations_are_logged(arrays, tmp_path):
    # Every look at the test block leaves a trace, so repeated peeking is visible.
    from shapely.geometry import box
    cands = extract_from_arrays(*arrays, T).assign(split="test")
    labels = labels_gdf(PITS + DUMPS, "test")
    region = box(X0, Y0 - N, X0 + N, Y0)
    log = tmp_path / "test_evaluations.csv"
    evaluate_split(cands, labels, DEFAULT_PARAMS, region, "test", final=True, log_path=log, ns=(5,), radii=(40,))
    evaluate_split(cands, labels, DEFAULT_PARAMS, region, "test", final=True, log_path=log, ns=(5,), radii=(40,))
    assert len(pd.read_csv(log)) == 2


def test_dump_max_height_drops_natural_knobs(arrays):
    # On real data, "dumps" taller than ~4 m above their 61 m surroundings were natural
    # knobs/outcrops (0 hits on the train cell). The cap is a tunable parameter.
    cands = extract_from_arrays(*arrays, T)
    tall = score_candidates(cands, {**DEFAULT_PARAMS, "dump_max_height": 99.0})
    capped = score_candidates(cands, {**DEFAULT_PARAMS, "dump_max_height": 1.0})
    assert (tall["kind"] == "dump").any()
    assert not (capped["kind"] == "dump").any()        # synthetic dumps are ~2+ m high


def test_tuning_breaks_ties_on_precision_at_2n():
    # With only 50 samples, many settings tie on p@50; the tie goes to the best p@100,
    # a steadier signal, instead of to whichever combination came first in the grid.
    lab = labels_gdf([(10 * i, 0) for i in range(4)], "train")
    xs = [X0 + 10 * i + 0.5 for i in range(4)]
    base = dict(kind="pit", area_m2=20.0, elongation=1.0, lrm_large_here=0.0, max_lrm_large_near=0.0, split="train")
    # Two strong hits (2.0 m), two weaker hits (0.5 m) and a deep miss far away.
    cands = gpd.GeoDataFrame(
        [{**base, "prominence_m": v} for v in (2.0, 2.0, 0.5, 0.5, 3.0)],
        geometry=[Point(xs[0], Y0 - 0.5), Point(xs[1], Y0 - 0.5), Point(xs[2], Y0 - 0.5),
                  Point(xs[3], Y0 - 0.5), Point(X0 + 900, Y0 - 900)], crs="EPSG:26913")
    # min depth 0.3 keeps the weak hits (better p@4); 1.0 drops them. Both tie at p@2.
    best, table = tune(cands, lab, grid={"pit_min_depth": [1.0, 0.3]}, n=2)
    assert best["pit_min_depth"] == 0.3


def test_flat_topped_feature_location_is_deterministic():
    # LRM is stored in whole cm, so a dump top is often a plateau of equal pixels.
    # scipy's maximum_position doesn't promise which tied pixel it returns (it differed
    # between scipy versions on Windows vs Linux), so ties go to the first pixel in
    # row-major order: top-most, then left-most. Translation-invariant, so blocks agree.
    from minefinder.detect.baseline import _components
    v = np.zeros((20, 20))
    v[5:12, 6:13] = 1.0
    v[7:10, 8:11] = 2.0                       # 3x3 plateau at the top
    _, df = _components(v > 0.5, v, "max", max_extent=50)
    assert (df.iloc[0]["row"], df.iloc[0]["col"]) == (7, 8)
    pits = -v
    _, dp = _components(pits < -0.5, pits, "min", max_extent=50)
    assert (dp.iloc[0]["row"], dp.iloc[0]["col"]) == (7, 8)
