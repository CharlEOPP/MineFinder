"""Do known-mine labels sit on landforms in the DEM? By how much are they off?

    python -m minefinder.review.alignment          # train-split labels, full mosaic

For every label, read a small DEM window around it, compute local relief, and score
how pit-like the ground is at the label and at every shift up to +/-40 m. Then compare
each label group against random points in the same DEM:

  score_at_label_m   mean depth (pits, adits) or height (dumps) of the strongest local
                     low/high within 4 m of the labels
  p_at_label         how often random point sets of the same size score that well
  best_dx/dy_m       the shift (east, north) where the group's mean local relief is
                     strongest, i.e. the vector to ADD to labels to land on features
  best_relief_m      how far below (pits) / above (dumps) its surroundings the group's
                     mean ground is there. Shift search runs for pits and dumps only.
  p_best_shift       how often random point sets reach that best score with their
                     own best shift. Searching 81x81 shifts always finds *something*;
                     this says whether the best shift is real or just search noise.

The report holds aggregate statistics only, never coordinates (see decision 0005).
Why it works this way: docs/decisions/0008-label-alignment-check.md
"""
from __future__ import annotations

import argparse
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio
from rasterio.windows import Window

from minefinder.terrain.derivatives import local_relief, pit_depth

REPO = Path(__file__).resolve().parents[3]

PIT_LIKE = {"mine shaft", "vertical shaft", "incline shaft", "prospect pit", "prospect hole", "subsidence feature"}
ADIT = {"adit"}
DUMP = {"mine dump", "tailings - pond", "tailings - undifferentiated"}

LRM_WINDOW = 31   # m; wider than any shaft/pit, narrower than a hillside
PIT_RADIUS = 4    # m
MIN_N = 5         # groups smaller than this are not reported

# Which way each landform class departs from its surroundings. +1 = below (pits,
# notches), -1 = above (dumps). The shift search only runs where "nearest strong low
# (or high)" really means "the feature": not for adits (the nearest low is usually
# the gully below) or "other" (mixed MRDS development-status records).
SIGN = {"pit-like": +1, "dump": -1, "adit": +1, "other": +1}
SHIFT_SEARCH = {"pit-like", "dump"}


def landform_class(labels: gpd.GeoDataFrame) -> pd.Series:
    ft = labels["feature_type"].fillna("").str.lower()
    cls = np.select([ft.isin(PIT_LIKE), ft.isin(ADIT), ft.isin(DUMP)], ["pit-like", "adit", "dump"], "other")
    return pd.Series(cls, index=labels.index)


def label_groups(labels: gpd.GeoDataFrame) -> pd.Series:
    """'<anchor source> / <landform class>', e.g. 'usmin / pit-like'.

    The anchor source (the label_id prefix) is the record whose location the label
    uses after the merge, so that's the source whose accuracy we're measuring.
    """
    anchor = labels["label_id"].str.split(":", n=1).str[0]
    return anchor + " / " + landform_class(labels)


def _shift_grids(ds, xs, ys, max_shift: int, sign: int = +1) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """For each usable point, two arrays indexed [dy+S, dx+S] (shift in metres, north/east +):

      depth  pit depth within PIT_RADIUS: tolerant, so a label on a pit's rim still counts.
             Used for "is the label on a pit?". It's flat-topped within 4 m of a pit, so it
             can't locate the centre: ties would resolve to one corner (south-west bias).
      relief negated local relief at the exact pixel: peaks sharply at a pit's centre.
             Used for "which shift lines labels up best?".

    sign=-1 flips relief, so the same machinery finds mounds (dumps) instead of pits.
    """
    S = max_shift
    H = S + PIT_RADIUS + LRM_WINDOW // 2 + 2          # window half-size in px
    depth, relief, ok = [], [], np.zeros(len(xs), dtype=bool)
    for i, (x, y) in enumerate(zip(xs, ys)):
        r, c = ds.index(x, y)
        if r - H < 0 or c - H < 0 or r + H >= ds.height or c + H >= ds.width:
            continue
        z = ds.read(1, window=Window(c - H, r - H, 2 * H + 1, 2 * H + 1)).astype("float64")
        if ds.nodata is not None:
            z[z == ds.nodata] = np.nan
        core = z[H - S - PIT_RADIUS : H + S + PIT_RADIUS + 1, H - S - PIT_RADIUS : H + S + PIT_RADIUS + 1]
        if not np.isfinite(core).all():
            continue
        lrm = sign * local_relief(z, LRM_WINDOW)
        dep = pit_depth(lrm, PIT_RADIUS)
        # Rows run south as the index grows, so flip to make axis 0 = dy (north positive).
        cut = (slice(H - S, H + S + 1), slice(H - S, H + S + 1))
        depth.append(dep[cut][::-1, :].astype("float32"))
        relief.append((-lrm[cut])[::-1, :].astype("float32"))
        ok[i] = True
    empty = np.empty((0, 2 * S + 1, 2 * S + 1), "float32")
    return (np.stack(depth) if depth else empty), (np.stack(relief) if relief else empty), ok


def _random_points(ds, n: int, rng: np.random.Generator, max_shift: int, sign: int = +1):
    """Random points with a full, gap-free window: the null to compare labels against."""
    W, S_, E, N = ds.bounds
    deps, rels, got, tries = [], [], 0, 0
    while got < n and tries < 50:
        tries += 1
        d, r, _ = _shift_grids(ds, rng.uniform(W, E, n), rng.uniform(S_, N, n), max_shift, sign)
        deps.append(d); rels.append(r); got += len(d)
    if got == 0:
        raise RuntimeError("no random point had a full DEM window; is the DEM mostly nodata?")
    return np.concatenate(deps)[:n], np.concatenate(rels)[:n]


def alignment_report(
    dem_path: str | Path,
    labels: gpd.GeoDataFrame,
    n_random: int = 3000,
    max_shift: int = 40,
    seed: int = 0,
    n_null: int = 2000,
    n_null_shift: int = 200,
) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    S = max_shift
    groups = label_groups(labels)
    classes = landform_class(labels)
    rows = []
    with rasterio.open(dem_path) as ds:
        if labels.crs is not None and ds.crs is not None and labels.crs != ds.crs:
            raise ValueError(f"labels CRS {labels.crs} != DEM CRS {ds.crs}; both must be EPSG:26913")
        pools: dict[int, tuple[np.ndarray, np.ndarray]] = {}
        for name in pd.unique(groups):
            mask = groups == name
            cls = classes[mask].iloc[0]
            sign = SIGN[cls]
            if sign not in pools:
                pools[sign] = _random_points(ds, n_random, rng, S, sign)
            rand_dep, rand_rel = pools[sign]
            rand0 = rand_dep[:, S, S]

            sel = labels[mask]
            D, R, ok = _shift_grids(ds, sel.geometry.x.to_numpy(), sel.geometry.y.to_numpy(), S, sign)
            n = len(D)
            if n < MIN_N:
                continue
            at_label = float(D[:, S, S].mean())
            null0 = rand0[rng.integers(0, len(rand0), (n_null, n))].mean(1)

            dx = dy = best = p_best = np.nan
            if cls in SHIFT_SEARCH:
                mean_rel = R.mean(0)
                i, j = np.unravel_index(np.argmax(mean_rel), mean_rel.shape)
                dx, dy, best = int(j - S), int(i - S), round(float(mean_rel[i, j]), 3)
                null_best = np.array([
                    rand_rel[rng.integers(0, len(rand_rel), n)].mean(0).max() for _ in range(n_null_shift)
                ])
                p_best = round(float((null_best >= best).mean()), 4)
            rows.append({
                "group": name,
                "n": n,
                "n_skipped": int((~ok).sum()),
                "score_at_label_m": round(at_label, 3),
                "random_mean_m": round(float(rand0.mean()), 3),
                "p_at_label": round(float((null0 >= at_label).mean()), 4),
                "best_dx_m": dx,
                "best_dy_m": dy,
                "best_relief_m": best,
                "p_best_shift": p_best,
            })
    return pd.DataFrame(rows)


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description="Label-to-landform alignment check")
    ap.add_argument("--dem", default=REPO / "data" / "interim" / "dem_1m.tif", type=Path)
    ap.add_argument("--labels", default=REPO / "data" / "interim" / "labels.gpkg", type=Path)
    ap.add_argument("--out", default=REPO / "data" / "processed" / "alignment.csv", type=Path)
    ap.add_argument("--split", default="train",
                    help="which label split to use (default train: anything we tune from this must not see test)")
    ap.add_argument("--n-random", default=3000, type=int)
    ap.add_argument("--seed", default=0, type=int)
    args = ap.parse_args(argv)

    labels = gpd.read_file(args.labels)
    if "split" in labels and args.split != "all":
        labels = labels[labels["split"] == args.split]
    rep = alignment_report(args.dem, labels, n_random=args.n_random, seed=args.seed)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    rep.to_csv(args.out, index=False)
    with pd.option_context("display.width", 200, "display.max_columns", 20):
        print(rep.to_string(index=False))
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
