"""precision@N against known mines, with a match radius and one-to-one matching."""
import geopandas as gpd
import numpy as np
import pandas as pd
import pytest
from shapely.geometry import Point, box

from minefinder.detect.evaluate import evaluate, match_candidates, precision_at_n, random_precision

CRS = "EPSG:26913"
X0, Y0 = 450_000, 4_420_000


def pts(xy, **cols):
    return gpd.GeoDataFrame(cols, geometry=[Point(x, y) for x, y in xy], crs=CRS)


def test_candidate_within_radius_is_a_hit():
    labels = pts([(X0, Y0)])
    cands = pts([(X0 + 30, Y0), (X0 + 500, Y0)], score=[2.0, 1.0])
    m = match_candidates(cands, labels, radius_m=40)
    assert list(m["hit"]) == [True, False]


def test_matching_is_one_to_one_by_score():
    # Three candidates on one dump are one find, not three. The best-scored one claims
    # the label; the others are duplicates and count as misses.
    labels = pts([(X0, Y0)])
    cands = pts([(X0 + 5, Y0), (X0 + 10, Y0), (X0 + 15, Y0)], score=[1.0, 3.0, 2.0])
    m = match_candidates(cands, labels, radius_m=40)
    assert m.sort_values("score", ascending=False)["hit"].tolist() == [True, False, False]


def test_a_candidate_takes_its_nearest_free_label():
    labels = pts([(X0, Y0), (X0 + 60, Y0)])
    cands = pts([(X0 + 35, Y0), (X0 + 5, Y0)], score=[2.0, 1.0])
    m = match_candidates(cands, labels, radius_m=40)
    assert m["hit"].all()                      # both labels found, neither double-claimed


def test_precision_at_n_ranks_by_score():
    labels = pts([(X0, Y0)])
    cands = pts([(X0 + 500, Y0), (X0, Y0), (X0 + 900, Y0)], score=[3.0, 2.0, 1.0])
    assert precision_at_n(cands, labels, n=1, radius_m=40) == 0.0
    assert precision_at_n(cands, labels, n=2, radius_m=40) == 0.5


def test_precision_at_n_with_fewer_candidates_than_n_divides_by_n():
    # Returning only 10 candidates shouldn't let a detector look better at @50.
    labels = pts([(X0, Y0)])
    cands = pts([(X0, Y0)], score=[1.0])
    assert precision_at_n(cands, labels, n=50, radius_m=40) == pytest.approx(1 / 50)


def test_random_precision_matches_label_coverage():
    # One label in a 1 km^2 block: a 40 m disk covers pi*40^2/1e6 ~= 0.50% of it.
    region = box(X0 - 500, Y0 - 500, X0 + 500, Y0 + 500)
    labels = pts([(X0, Y0)])
    p = random_precision(region, labels, radius_m=40, n_samples=40_000, seed=0)
    assert p == pytest.approx(np.pi * 40**2 / 1e6, rel=0.15)


def test_evaluate_only_scores_candidates_and_labels_in_the_region():
    region = box(X0, Y0, X0 + 1000, Y0 + 1000)
    labels = pts([(X0 + 100, Y0 + 100), (X0 + 5000, Y0)], sources=["usmin", "usmin"])
    cands = pts([(X0 + 100, Y0 + 100), (X0 + 5000, Y0)], score=[1.0, 2.0])  # 2nd is outside
    rep = evaluate(cands, labels, region, ns=(1,), radii=(40,), n_random=2000)
    r = rep.iloc[0]
    assert r["n_candidates"] == 1 and r["n_labels"] == 1
    assert r["precision"] == 1.0


def test_evaluate_reports_trusted_label_subset():
    # MRDS-only labels are coarse; report precision against USMIN+CGS separately too.
    region = box(X0, Y0, X0 + 1000, Y0 + 1000)
    labels = pts([(X0 + 100, Y0 + 100), (X0 + 600, Y0 + 600)], sources=["mrds", "usmin"])
    cands = pts([(X0 + 100, Y0 + 100), (X0 + 600, Y0 + 600)], score=[2.0, 1.0])
    rep = evaluate(cands, labels, region, ns=(2,), radii=(40,), n_random=2000).set_index("labels")
    assert rep.loc["all", "precision"] == 1.0
    assert rep.loc["usmin+cgs", "precision"] == 0.5


def test_evaluate_report_has_no_coordinates():
    region = box(X0, Y0, X0 + 1000, Y0 + 1000)
    labels = pts([(X0 + 100, Y0 + 100)], sources=["usmin"])
    cands = pts([(X0 + 100, Y0 + 100)], score=[1.0])
    rep = evaluate(cands, labels, region, ns=(1,), radii=(40,), n_random=500)
    assert not {"geometry", "x", "y"} & set(rep.columns)
    num = rep.select_dtypes("number")
    assert (num.abs() < 1e5).all().all()           # nothing coordinate-sized in the report
