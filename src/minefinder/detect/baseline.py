"""Rule-based baseline detector: the number every model has to beat (Phase 3).

    python -m minefinder.detect.baseline extract           # all loose candidates -> outputs/ (gitignored)
    python -m minefinder.detect.baseline tune              # pick thresholds on the TRAIN region only
    python -m minefinder.detect.baseline evaluate          # dev score on the train region
    python -m minefinder.detect.baseline evaluate --final  # the held-out test block (logged, every time)

Two candidate kinds, each read at its own scale (decision 0009):
  pit   connected low in lrm_small (21 m): shafts, prospect pits. Located at the deepest pixel.
  dump  connected high in lrm_large (61 m): waste dumps, tailings. Located at the highest pixel.

Look-alikes are rejected by shape, not by place:
  roads, ditches, ridges, gullies   long: wider than max_extent_m, or too elongated
  dump halos (0009)                 shallow, tiny lows ringing a big high; optionally dropped
                                    if a strong dump is within 12 m and the pit isn't on it
  survey seams, nodata edges        seam_mask

Extraction keeps everything above loose floors and records the measurements. Tuning
and scoring then only filter and rank that table, so trying 100 threshold combinations
costs seconds, not 100 terrain passes.

Why it works this way: docs/decisions/0010-baseline-detector-and-evaluation.md
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from datetime import datetime, timezone
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio
import yaml
from rasterio.windows import Window
from scipy import ndimage as ndi

from minefinder.aoi import PROJECT_CRS, load_aoi
from minefinder.detect.evaluate import evaluate, precision_at_n

REPO = Path(__file__).resolve().parents[3]

# Loose floors for extraction: well below any threshold tuning will try.
LOOSE_PIT_DEPTH = 0.2      # m below the 21 m neighbourhood
LOOSE_DUMP_HEIGHT = 0.3    # m above the 61 m neighbourhood
MIN_PX = 4                 # smaller blobs are noise (and their shape is meaningless)
MAX_EXTENT = 120           # m; anything wider is a road/ridge/gully at this baseline
HALO = 128                 # px of context per block; must be >= MAX_EXTENT
HALO_RADIUS = 12           # m: "is a strong dump right next to this pit?"

DEFAULT_PARAMS = {
    "pit_min_depth": 0.5,
    "pit_min_area": 4,
    "pit_max_area": 400,
    "dump_min_height": 1.0,
    "dump_max_height": 4.0,  # taller "dumps" were natural knobs on the train cell (0010)
    "dump_min_area": 50,
    "dump_max_area": 8000,
    "max_elongation": 3.0,
    "dump_weight": 0.5,      # dump height (m) * weight is compared with pit depth (m)
    "halo_guard": True,
}

DEFAULT_GRID = {
    "pit_min_depth": [0.3, 0.5, 0.8],
    "dump_min_height": [0.5, 1.0, 2.0],
    "dump_max_height": [4.0, 8.0, 99.0],
    "max_elongation": [2.0, 3.0, 4.0],
    "dump_weight": [0.25, 0.5, 1.0],
    "halo_guard": [True, False],
}


# ---------------------------------------------------------------- extraction

def _components(mask: np.ndarray, value: np.ndarray, want: str, max_extent: int):
    """Connected components (8-neighbour) with location of extreme value, area,
    elongation and bbox extent. `want` = 'min' (pits) or 'max' (dumps)."""
    lab, n = ndi.label(mask, structure=np.ones((3, 3), bool))
    if n == 0:
        return lab, pd.DataFrame(columns=["row", "col", "area", "elongation", "ext"])
    idx = np.arange(1, n + 1)
    area = ndi.sum(np.ones_like(lab), lab, idx)
    rr, cc = np.indices(lab.shape)
    m_r, m_c = ndi.mean(rr, lab, idx), ndi.mean(cc, lab, idx)
    v_rr = ndi.mean(rr.astype("float64") ** 2, lab, idx) - m_r**2
    v_cc = ndi.mean(cc.astype("float64") ** 2, lab, idx) - m_c**2
    v_rc = ndi.mean((rr * cc).astype("float64"), lab, idx) - m_r * m_c
    # Eigenvalues of the 2x2 covariance: sqrt(major/minor) = elongation (1 = round).
    tr, det = v_rr + v_cc, v_rr * v_cc - v_rc**2
    disc = np.sqrt(np.maximum(tr**2 / 4 - det, 0))
    lo, hi = np.maximum(tr / 2 - disc, 1e-6), tr / 2 + disc
    # Location = the extreme pixel. LRM is whole centimetres, so extremes are often
    # plateaus of tied pixels, and scipy's *_position doesn't promise which tie it
    # returns (it varies by version). So: find the extreme value, then take the first
    # tied pixel in row-major order. Translation-invariant, so every block agrees.
    ext_fn = ndi.minimum if want == "min" else ndi.maximum
    ext_val = np.asarray(ext_fn(value, lab, idx), dtype="float64")
    full = np.concatenate([[np.nan], ext_val])
    tied = (lab > 0) & (value == full[lab])
    order = rr.astype("int64") * lab.shape[1] + cc
    key = np.where(tied, order, np.iinfo("int64").max)
    first = np.asarray(ndi.minimum(key, lab, idx), dtype="int64")
    pos = np.stack([first // lab.shape[1], first % lab.shape[1]], axis=1)
    objs = ndi.find_objects(lab)
    ext = np.array([max(s[0].stop - s[0].start, s[1].stop - s[1].start) for s in objs])
    df = pd.DataFrame({"row": pos[:, 0], "col": pos[:, 1], "area": area,
                       "elongation": np.sqrt(hi / lo), "ext": ext})
    keep = (df["area"] >= MIN_PX) & (df["ext"] <= max_extent)
    return lab, df[keep]


def _extract(lrm_s, lrm_l, seam, max_extent):
    """Pixel-space candidates for one array (whole raster or one padded block)."""
    s = np.where(np.isfinite(lrm_s), lrm_s, 0.0)
    l = np.where(np.isfinite(lrm_l), lrm_l, 0.0)
    rows = []
    _, pits = _components(s < -LOOSE_PIT_DEPTH, s, "min", max_extent)
    if len(pits):
        near_l = ndi.maximum_filter(np.where(np.isfinite(lrm_l), lrm_l, -np.inf), size=2 * HALO_RADIUS + 1,
                                    mode="constant", cval=-np.inf)
        pits = pits.assign(kind="pit", prominence_m=-s[pits.row, pits.col],
                           lrm_large_here=l[pits.row, pits.col], max_lrm_large_near=near_l[pits.row, pits.col])
        rows.append(pits)
    _, dumps = _components(l > LOOSE_DUMP_HEIGHT, l, "max", max_extent)
    if len(dumps):
        dumps = dumps.assign(kind="dump", prominence_m=l[dumps.row, dumps.col],
                             lrm_large_here=l[dumps.row, dumps.col], max_lrm_large_near=l[dumps.row, dumps.col])
        rows.append(dumps)
    if not rows:
        return pd.DataFrame(columns=["row", "col", "kind"])
    df = pd.concat(rows, ignore_index=True)
    return df[seam[df.row, df.col] == 0]


def _to_gdf(df: pd.DataFrame, transform) -> gpd.GeoDataFrame:
    xs, ys = rasterio.transform.xy(transform, df["row"].to_numpy(), df["col"].to_numpy(), offset="center")
    out = df.rename(columns={"area": "area_m2"}).drop(columns=["row", "col", "ext"])
    return gpd.GeoDataFrame(out.reset_index(drop=True), geometry=gpd.points_from_xy(xs, ys), crs=PROJECT_CRS)


def extract_from_arrays(lrm_s, lrm_l, seam, transform, max_extent_m: int = MAX_EXTENT) -> gpd.GeoDataFrame:
    return _to_gdf(_extract(lrm_s, lrm_l, seam, max_extent_m), transform)


def _read(ds, win, nodata_to_nan=True):
    a = ds.read(1, window=win, boundless=True, fill_value=ds.nodata if ds.nodata is not None else 0)
    if not nodata_to_nan:
        return a
    a = a.astype("float64")
    if ds.nodata is not None:
        a[a == ds.nodata] = np.nan
    return a * (ds.scales[0] or 1.0)


def extract_candidates(terrain_dir: str | Path, block: int = 2048, halo: int = HALO,
                       max_extent_m: int = MAX_EXTENT) -> gpd.GeoDataFrame:
    """Blockwise extraction over the terrain rasters. Each block is read with `halo`
    px of context; a component is kept only by the block whose core holds its location,
    so nothing is counted twice and results equal a whole-raster run."""
    if halo < max_extent_m + 2:
        raise ValueError(f"halo ({halo}) must exceed max_extent_m ({max_extent_m})")
    d = Path(terrain_dir)
    parts = []
    with rasterio.open(d / "lrm_small.tif") as S, rasterio.open(d / "lrm_large.tif") as L, \
         rasterio.open(d / "seam_mask.tif") as M:
        if S.crs.to_epsg() != 26913 or not (S.transform == L.transform == M.transform):
            raise ValueError("terrain rasters must be EPSG:26913 on one grid")
        H, W = S.height, S.width
        for row in range(0, H, block):
            for col in range(0, W, block):
                h, w = min(block, H - row), min(block, W - col)
                win = Window(col - halo, row - halo, w + 2 * halo, h + 2 * halo)
                s, l = _read(S, win), _read(L, win)
                m = M.read(1, window=win, boundless=True, fill_value=1)
                df = _extract(s, l, m, max_extent_m)
                core = df.row.between(halo, halo + h - 1) & df.col.between(halo, halo + w - 1)
                df = df[core].copy()
                df["row"] += row - halo
                df["col"] += col - halo
                parts.append(df)
        transform = S.transform
    allc = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame(columns=["row", "col", "kind"])
    return _to_gdf(allc, transform)


# ---------------------------------------------------------------- scoring + tuning

def score_candidates(cands: gpd.GeoDataFrame, params: dict) -> gpd.GeoDataFrame:
    """Apply thresholds and rank. Score is in 'metres of prominence': pit depth, or
    dump height x dump_weight, so the two kinds share one ranked list."""
    p = {**DEFAULT_PARAMS, **params}
    pit = cands["kind"] == "pit"
    dump = cands["kind"] == "dump"
    ok_pit = pit & (cands["prominence_m"] >= p["pit_min_depth"]) \
        & cands["area_m2"].between(p["pit_min_area"], p["pit_max_area"])
    if p["halo_guard"]:
        # Beside a strong dump but not on it -> likely the 21 m window's halo, not a pit.
        halo = (cands["max_lrm_large_near"] > p["dump_min_height"]) & (cands["lrm_large_here"] < 0.2)
        ok_pit &= ~halo
    ok_dump = dump & cands["prominence_m"].between(p["dump_min_height"], p["dump_max_height"]) \
        & cands["area_m2"].between(p["dump_min_area"], p["dump_max_area"])
    keep = (ok_pit | ok_dump) & (cands["elongation"] <= p["max_elongation"])
    out = cands[keep].copy()
    out["score"] = np.where(out["kind"] == "dump", out["prominence_m"] * p["dump_weight"], out["prominence_m"])
    return out.sort_values("score", ascending=False, kind="stable")


def tune(cands, labels, grid: dict | None = None, n: int = 50, radius_m: float = 40.0):
    """Grid search on the TRAIN split only. Returns (best params, full results table).
    Test candidates and test labels are dropped before anything is computed.
    Ties on p@n (common: only n samples) go to the best p@2n, then to grid order."""
    grid = grid or DEFAULT_GRID
    c = cands[cands["split"] == "train"]
    lab = labels[labels["split"] == "train"]
    keys = list(grid)
    rows = []
    for combo in itertools.product(*(grid[k] for k in keys)):
        params = {**DEFAULT_PARAMS, **dict(zip(keys, combo))}
        scored = score_candidates(c, params)
        rows.append({**dict(zip(keys, combo)),
                     f"train_p@{2 * n}": precision_at_n(scored, lab, 2 * n, radius_m),
                     f"train_p@{n}": precision_at_n(scored, lab, n, radius_m)})
    table = pd.DataFrame(rows)
    order = table.sort_values([f"train_p@{n}", f"train_p@{2 * n}"], ascending=False, kind="stable")
    best = order.iloc[0]
    best_params = {**DEFAULT_PARAMS, **{k: (bool(best[k]) if isinstance(DEFAULT_PARAMS[k], bool) else
                                            type(DEFAULT_PARAMS[k])(best[k])) for k in keys}}
    return best_params, table


def evaluate_split(cands, labels, params, region, split: str, final: bool = False,
                   log_path: str | Path | None = None, ns=(25, 50, 100), radii=(25, 40, 60)) -> pd.DataFrame:
    """Score one split. The test split needs final=True, and every test evaluation is
    appended to a log, so peeking at the test block repeatedly is visible in git."""
    if split == "test" and not final:
        raise PermissionError("scoring the held-out test block needs final=True (--final): "
                              "tune on train first, then look at test once")
    c = score_candidates(cands[cands["split"] == split], params)
    lab = labels[labels["split"] == split]
    rep = evaluate(c, lab, region, ns=ns, radii=radii)
    if split == "test" and log_path is not None:
        log_path = Path(log_path)
        head = rep.query("labels == 'all'")
        summary = {f"p@{r.n}@{r.radius_m}m": r.precision for r in head.itertuples()}
        entry = pd.DataFrame([{
            "utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "detector": "baseline",
            "params_hash": hashlib.sha1(json.dumps(params, sort_keys=True).encode()).hexdigest()[:10],
            "params": json.dumps(params, sort_keys=True),
            **summary,
        }])
        entry.to_csv(log_path, mode="a", header=not log_path.exists(), index=False)
    return rep


# ---------------------------------------------------------------- CLI

def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description="Rule-based baseline detector")
    ap.add_argument("command", choices=["extract", "tune", "evaluate"])
    ap.add_argument("--config", default=REPO / "config" / "aoi_ward_goldhill.yaml", type=Path)
    ap.add_argument("--terrain", default=REPO / "data" / "processed" / "terrain", type=Path)
    ap.add_argument("--labels", default=REPO / "data" / "interim" / "labels.gpkg", type=Path)
    # Candidates carry exact coordinates of possibly-unrecorded mines -> outputs/ (gitignored).
    ap.add_argument("--candidates", default=REPO / "outputs" / "baseline_candidates.gpkg", type=Path)
    ap.add_argument("--params", default=REPO / "config" / "baseline_params.yaml", type=Path)
    ap.add_argument("--final", action="store_true", help="score the held-out test block (logged)")
    args = ap.parse_args(argv)
    aoi = load_aoi(args.config)

    if args.command == "extract":
        c = extract_candidates(args.terrain)
        c["split"] = c.geometry.map(aoi.split_of)
        args.candidates.parent.mkdir(parents=True, exist_ok=True)
        c.to_file(args.candidates, driver="GPKG")
        print(f"wrote {len(c)} loose candidates -> {args.candidates}")
        print(c.groupby(["split", "kind"]).size().to_string())
        return

    cands = gpd.read_file(args.candidates)
    labels = gpd.read_file(args.labels)

    if args.command == "tune":
        best, table = tune(cands, labels)
        out = REPO / "data" / "processed" / "baseline_tuning.csv"
        table.sort_values(table.columns[-1], ascending=False).to_csv(out, index=False)
        args.params.write_text(yaml.safe_dump(best, sort_keys=True))
        print(table.sort_values(table.columns[-1], ascending=False).head(10).to_string(index=False))
        print(f"\nbest (train only) -> {args.params}\n{yaml.safe_dump(best, sort_keys=True)}")
        return

    params = yaml.safe_load(args.params.read_text()) if args.params.exists() else None
    if params is None:
        raise SystemExit(f"{args.params} missing: run `tune` first")
    split, region = ("test", aoi.holdout) if args.final else ("train", aoi.train_region)
    rep = evaluate_split(cands, labels, params, region, split, final=args.final,
                         log_path=REPO / "docs" / "test_evaluations.csv")
    out = REPO / "data" / "processed" / f"baseline_eval_{split}.csv"
    rep.to_csv(out, index=False)
    print(rep.to_string(index=False))
    print(f"\nwrote {out}" + ("  (logged in docs/test_evaluations.csv)" if args.final else ""))


if __name__ == "__main__":
    main()
