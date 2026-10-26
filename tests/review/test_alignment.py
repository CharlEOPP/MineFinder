"""Label-to-landform alignment. Synthetic DEM with pits at known spots, so we
know the right answer: labels on pits, labels offset by a known vector, and
labels placed at random."""
import geopandas as gpd
import numpy as np
import pandas as pd
import pytest
import rasterio
from rasterio.transform import from_origin
from shapely.geometry import Point

from minefinder.review.alignment import alignment_report, label_groups, main

X0, Y0 = 450_000, 4_430_000        # upper-left of the synthetic DEM
N = 800                             # 800 x 800 m at 1 m


@pytest.fixture(scope="module")
def pits():
    rng = np.random.default_rng(1)
    # 60 pit centres away from the edges, at least 40 m apart so they don't overlap
    pts = []
    while len(pts) < 60:
        p = rng.uniform(80, N - 80, 2)
        if all(np.hypot(*(p - q)) > 40 for q in pts):
            pts.append(p)
    return np.array(pts)  # (col, row) in pixels


@pytest.fixture(scope="module")
def dem(tmp_path_factory, pits):
    rng = np.random.default_rng(2)
    yy, xx = np.mgrid[0:N, 0:N].astype("float32")
    z = 2800 + 0.25 * xx - 0.15 * yy + rng.normal(0, 0.05, (N, N)).astype("float32")
    for c, r in pits:  # 2.5 m deep, ~4 m radius bowls
        d = np.hypot(xx - c, yy - r)
        z -= np.where(d < 4, 2.5 * (1 - (d / 4) ** 2), 0).astype("float32")
    path = tmp_path_factory.mktemp("dem") / "dem.tif"
    with rasterio.open(path, "w", driver="GTiff", width=N, height=N, count=1, dtype="float32",
                       crs="EPSG:26913", transform=from_origin(X0, Y0, 1, 1), nodata=-999999.0) as dst:
        dst.write(z, 1)
    return path


def labels_at(pits, dx=0.0, dy=0.0, source="usmin", ftype="mine shaft"):
    """Pixel (col,row) -> map x,y, then offset by (dx, dy) metres (east, north)."""
    geom = [Point(X0 + c + 0.5 + dx, Y0 - r - 0.5 + dy) for c, r in pits]
    return gpd.GeoDataFrame(
        {"label_id": [f"{source}:{i}" for i in range(len(pits))],
         "feature_type": ftype, "sources": source},
        geometry=geom, crs="EPSG:26913",
    )


@pytest.fixture(scope="module")
def mound_dem(tmp_path_factory, pits):
    """Same layout, but 2 m high, ~6 m radius mounds (like waste dumps) instead of pits."""
    rng = np.random.default_rng(3)
    yy, xx = np.mgrid[0:N, 0:N].astype("float32")
    z = 2800 + 0.25 * xx - 0.15 * yy + rng.normal(0, 0.05, (N, N)).astype("float32")
    for c, r in pits:
        d = np.hypot(xx - c, yy - r)
        z += np.where(d < 6, 2.0 * (1 - (d / 6) ** 2), 0).astype("float32")
    path = tmp_path_factory.mktemp("mdem") / "mdem.tif"
    with rasterio.open(path, "w", driver="GTiff", width=N, height=N, count=1, dtype="float32",
                       crs="EPSG:26913", transform=from_origin(X0, Y0, 1, 1), nodata=-999999.0) as dst:
        dst.write(z, 1)
    return path


def row(report, group):
    return report.set_index("group").loc[group]


def test_labels_on_pits_score_far_above_random(dem, pits):
    rep = alignment_report(dem, labels_at(pits), n_random=800, seed=0)
    r = row(rep, "usmin / pit-like")
    assert r["n"] == 60
    assert r["score_at_label_m"] > 2 * r["random_mean_m"]
    assert r["p_at_label"] < 0.01
    assert (r["best_dx_m"], r["best_dy_m"]) == (0, 0)


def test_recovers_a_known_offset(dem, pits):
    # Labels drawn 12 m east and 8 m south of their pits: the fix is (-12, +8).
    rep = alignment_report(dem, labels_at(pits, dx=12, dy=-8), n_random=800, seed=0)
    r = row(rep, "usmin / pit-like")
    assert abs(r["best_dx_m"] - (-12)) <= 1 and abs(r["best_dy_m"] - 8) <= 1
    assert r["p_best_shift"] < 0.01                  # a real offset, not search noise


def test_random_labels_are_not_significant(dem):
    rng = np.random.default_rng(5)
    fake = rng.uniform(80, N - 80, (60, 2))
    rep = alignment_report(dem, labels_at(fake), n_random=800, seed=0)
    r = row(rep, "usmin / pit-like")
    assert r["p_at_label"] > 0.05
    # Searching 81x81 shifts always finds *some* best shift; the null must catch that.
    assert r["p_best_shift"] > 0.05


def test_dumps_are_scored_as_mounds_and_offset_recovered(mound_dem, pits):
    # Dumps stand above their surroundings, so the search looks for high ground, not pits.
    lab = labels_at(pits, dx=-10, dy=6, source="co_inactive", ftype="mine dump")
    r = row(alignment_report(mound_dem, lab, n_random=800, seed=0), "co_inactive / dump")
    assert abs(r["best_dx_m"] - 10) <= 1 and abs(r["best_dy_m"] - (-6)) <= 1
    assert r["p_best_shift"] < 0.01


def test_no_shift_search_for_adits_or_other(dem, pits):
    # An adit is a notch on a slope; "lower ground nearby" is usually just the gully,
    # so a shift search would measure streams, not label error. Score only.
    lab = pd.concat([labels_at(pits, ftype="adit"), labels_at(pits, source="mrds", ftype="past producer")],
                    ignore_index=True)
    rep = alignment_report(dem, lab, n_random=500, seed=0).set_index("group")
    for g in ("usmin / adit", "mrds / other"):
        assert pd.isna(rep.loc[g, "best_dx_m"]) and pd.isna(rep.loc[g, "p_best_shift"])
        assert rep.loc[g, "score_at_label_m"] > 0


def test_skips_labels_without_full_dem_window(dem, pits):
    lab = pd.concat([labels_at(pits), labels_at(np.array([[5.0, 5.0]]))], ignore_index=True)
    rep = alignment_report(dem, lab, n_random=500, seed=0)
    r = row(rep, "usmin / pit-like")
    assert r["n"] == 60 and r["n_skipped"] == 1


def test_groups_by_anchor_source_and_landform_class():
    g = gpd.GeoDataFrame(
        {"label_id": ["usmin:a", "co_inactive:b", "co_inactive:c", "mrds:d"],
         "feature_type": ["mine shaft", "adit", "mine dump", "past producer"],
         "sources": ["usmin", "co_inactive", "co_inactive,usmin", "mrds"]},
        geometry=[Point(0, 0)] * 4, crs="EPSG:26913",
    )
    assert list(label_groups(g)) == ["usmin / pit-like", "co_inactive / adit", "co_inactive / dump", "mrds / other"]


def test_report_contains_no_coordinates(dem, pits, tmp_path):
    # Reports get shared; only aggregate stats are allowed in them.
    rep = alignment_report(dem, labels_at(pits), n_random=300, seed=0)
    text = rep.to_csv()
    assert str(X0)[:4] not in text and str(Y0)[:4] not in text


def test_cli_writes_csv(dem, pits, tmp_path):
    lab = tmp_path / "labels.gpkg"
    labels_at(pits).to_file(lab, driver="GPKG")
    out = tmp_path / "alignment.csv"
    main(["--dem", str(dem), "--labels", str(lab), "--out", str(out), "--n-random", "300"])
    df = pd.read_csv(out)
    assert "usmin / pit-like" in set(df["group"])
