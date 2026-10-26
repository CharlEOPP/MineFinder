"""precision@N against known mines: the project's headline metric.

A candidate is a **hit** if it lies within `radius_m` of a known-mine label that no
higher-scored candidate has already claimed (one-to-one matching). So five candidates
on one dump count as one find plus four misses, which is what matters for deciding
whether a field check is worth the hike.

Always read precision next to `random_precision`: the chance that a randomly placed
point lands within `radius_m` of some label in the same region. A detector is only
interesting by how far it beats that (`lift`).

Known-mine precision is a *lower bound* on real precision: a candidate on a real but
unrecorded mine counts as a miss here. Those are exactly the finds we want, and Phase 5
review sorts them out. Why it works this way: docs/decisions/0010-baseline-detector-and-evaluation.md
"""
from __future__ import annotations

import geopandas as gpd
import numpy as np
import pandas as pd
from shapely import STRtree
from shapely.geometry.base import BaseGeometry

TRUSTED = ("usmin", "co_inactive")   # MRDS-only labels are coarse (decision 0006)


def match_candidates(cands: gpd.GeoDataFrame, labels: gpd.GeoDataFrame, radius_m: float = 40.0) -> gpd.GeoDataFrame:
    """Greedy one-to-one match, best score first; each candidate takes its nearest free label."""
    out = cands.sort_values("score", ascending=False, kind="stable").copy()
    hit = np.zeros(len(out), dtype=bool)
    if len(labels):
        lab = labels.geometry.to_numpy()
        tree = STRtree(lab)
        taken = np.zeros(len(lab), dtype=bool)
        for k, g in enumerate(out.geometry.to_numpy()):
            near = tree.query(g, predicate="dwithin", distance=radius_m)
            free = [i for i in near if not taken[i]]
            if free:
                i = min(free, key=lambda j: g.distance(lab[j]))
                taken[i] = True
                hit[k] = True
    out["hit"] = hit
    return out


def precision_at_n(cands, labels, n: int = 50, radius_m: float = 40.0) -> float:
    """Hits among the top-n candidates, divided by n (not by the number returned)."""
    top = cands.sort_values("score", ascending=False, kind="stable").head(n)
    return float(match_candidates(top, labels, radius_m)["hit"].sum()) / n


def random_precision(region: BaseGeometry, labels, radius_m: float = 40.0,
                     n_samples: int = 20_000, seed: int = 0) -> float:
    """Chance that a uniformly random point in `region` is within radius_m of a label."""
    if len(labels) == 0:
        return 0.0
    rng = np.random.default_rng(seed)
    minx, miny, maxx, maxy = region.bounds
    pts: list = []
    while len(pts) < n_samples:
        xs = rng.uniform(minx, maxx, n_samples)
        ys = rng.uniform(miny, maxy, n_samples)
        g = gpd.points_from_xy(xs, ys)
        pts.extend(g[gpd.GeoSeries(g).within(region).to_numpy()])
    pts = np.array(pts[:n_samples], dtype=object)
    tree = STRtree(labels.geometry.to_numpy())
    near = tree.query(pts, predicate="dwithin", distance=radius_m)
    return float(len(np.unique(near[0]))) / len(pts)


def _trusted(labels: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    s = labels["sources"].astype(str)
    return labels[s.str.contains("|".join(TRUSTED))]


def evaluate(cands, labels, region: BaseGeometry, ns=(25, 50, 100), radii=(25, 40, 60),
             n_random: int = 20_000, seed: int = 0) -> pd.DataFrame:
    """precision@n for every n x radius, against all labels and against USMIN+CGS only.
    Candidates and labels outside `region` are ignored. Aggregate numbers only."""
    c = cands[cands.geometry.within(region)]
    lab_all = labels[labels.geometry.within(region)]
    rows = []
    for name, lab in (("all", lab_all), ("usmin+cgs", _trusted(lab_all))):
        for r in radii:
            rand = random_precision(region, lab, r, n_random, seed)
            for n in ns:
                p = precision_at_n(c, lab, n, r)
                hits = int(round(p * n))
                rows.append({
                    "labels": name, "n": n, "radius_m": r,
                    "n_candidates": len(c), "n_labels": len(lab),
                    "hits": hits, "precision": round(p, 4),
                    "recall": round(hits / len(lab), 4) if len(lab) else 0.0,
                    "random_precision": round(rand, 4),
                    "lift": round(p / rand, 2) if rand > 0 else np.nan,
                })
    return pd.DataFrame(rows)
