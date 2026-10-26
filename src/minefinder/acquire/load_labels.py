"""Read the three raw known-mine inventories, clean source-specific quirks, merge.

Each loader returns the normalized schema from labels.normalize (EPSG:26913 points).
Source quirks found in the real downloads, and why each is handled here, are in
docs/decisions/0006-label-sources-and-loader-cleaning.md.

Run:  python -m minefinder.acquire.load_labels   ->  data/interim/labels.gpkg
"""
from __future__ import annotations

import argparse
import hashlib
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
from shapely import STRtree

from minefinder.acquire.labels import merge_labels, normalize
from minefinder.aoi import AOI, PROJECT_CRS, load_aoi

REPO = Path(__file__).resolve().parents[3]

# USMIN feature types that are aggregate/borrow pits for roads, not hard-rock workings.
AGGREGATE_PIT_WORDS = ("gravel", "borrow", "sand pit", "clay pit", "cinder")

CGS_LAYERS = {
    "openings": "Hole_Mine_Openings_Final",
    "dumps": "PILE_Mine_Tailings_Final",
    "areas": "HDR_Mine_Inventory_Areas_Final",
}

# Two points closer than this are the same GPS fix / geocode (not two features).
COINCIDENT_M = 1.0


def _far_from(a: gpd.GeoDataFrame, b: gpd.GeoDataFrame, radius_m: float) -> np.ndarray:
    """Boolean mask over `a`: True where no point of `b` lies within radius_m."""
    if len(a) == 0 or len(b) == 0:
        return np.ones(len(a), dtype=bool)
    tree = STRtree(b.geometry.to_numpy())
    hits = tree.query(a.geometry.to_numpy(), predicate="dwithin", distance=radius_m)
    near = np.zeros(len(a), dtype=bool)
    near[hits[0]] = True
    return ~near


# ---------------------------------------------------------------- USMIN

def _usmin_id(layer: str, row: pd.Series, x: float, y: float) -> str:
    # USMIN has no feature ID, so derive one from content. Hashed so it is stable
    # across runs but never exposes coordinates if it ends up in an export.
    key = f"{layer}|{row['FTR_TYPE']}|{row['TOPO_NAME']}|{row['TOPO_SCALE']}|{x:.1f}|{y:.1f}"
    return hashlib.sha1(key.encode()).hexdigest()[:12]


def load_usmin(usmin_dir: str | Path, scale_dedupe_m: float = 50.0) -> gpd.GeoDataFrame:
    """USMIN topo-map mine symbols (ver 10.0 state files: CO-point.shp, CO-poly.shp)."""
    usmin_dir = Path(usmin_dir)
    parts = []
    for layer in ("point", "poly"):
        path = usmin_dir / f"CO-{layer}.shp"
        if not path.exists():
            continue
        g = gpd.read_file(path)
        g = g[~g["FTR_TYPE"].str.lower().str.contains("|".join(AGGREGATE_PIT_WORDS), na=False)]
        # Content hash computed in the source CRS so it doesn't depend on our reprojection.
        rp = g.geometry.representative_point()
        ids = [_usmin_id(layer, r, p.x, p.y) for (_, r), p in zip(g.iterrows(), rp)]
        g = g.assign(_id=ids, _scale=g["TOPO_SCALE"].astype(str))
        parts.append(g)
    g = pd.concat(parts, ignore_index=True)

    # Exact-duplicate records (same content) would collide; suffix keeps IDs unique and stable.
    dup = g.groupby("_id").cumcount()
    g["_id"] = g["_id"].where(dup == 0, g["_id"] + "-" + dup.astype(str))

    out = normalize(g, "usmin", "_id", "FTR_NAME", type_col="FTR_TYPE")
    out["_scale"] = g["_scale"].to_numpy()

    # ver 10.0 never reconciled map scales: drop 1:62.5k/1:48k symbols that duplicate a 1:24k one.
    fine = out[out["_scale"] == "24000"]
    coarse = out[out["_scale"] != "24000"]
    coarse = coarse[_far_from(coarse, fine, scale_dedupe_m)]
    return pd.concat([fine, coarse]).drop(columns="_scale").reset_index(drop=True)


# ---------------------------------------------------------------- MRDS

def load_mrds(csv_path: str | Path) -> gpd.GeoDataFrame:
    """USGS MRDS national CSV, filtered to Colorado. Coordinates are WGS84 (per mrds.met)."""
    cols = ["dep_id", "site_name", "latitude", "longitude", "state", "dev_stat"]
    df = pd.read_csv(csv_path, usecols=cols, dtype={"dep_id": str}, low_memory=False)
    df = df[(df["state"] == "Colorado") & df["latitude"].notna() & df["longitude"].notna()]

    # Stacked records at an identical coordinate are one coarse geocode: keep lowest dep_id.
    df = df.assign(_k=df["dep_id"].astype(int)).sort_values("_k")
    df = df.drop_duplicates(subset=["latitude", "longitude"], keep="first")

    g = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df["longitude"], df["latitude"]), crs="EPSG:4326")
    return normalize(g, "mrds", "dep_id", "site_name", type_col="dev_stat")


# ---------------------------------------------------------------- CGS USFS AML inventory

def load_co_aml(gdb_path: str | Path) -> gpd.GeoDataFrame:
    """CGS ON-008-04D Abandoned Mine Land Inventory (USFS lands): openings + dumps."""
    holes = gpd.read_file(gdb_path, layer=CGS_LAYERS["openings"])
    piles = gpd.read_file(gdb_path, layer=CGS_LAYERS["dumps"])
    areas = gpd.read_file(gdb_path, layer=CGS_LAYERS["areas"], ignore_geometry=True)

    # Unique_ID = area UniqueID + 3-digit feature number (HID/PID). The area holds the site name.
    names = dict(zip(areas["UniqueID"].astype(str), areas["SITENAME_1"]))

    def prep(g, type_col):
        uid = g["Unique_ID"].astype(str)
        return g.assign(_id=uid, _name=uid.str[:-3].map(names), _type=g[type_col])

    holes, piles = prep(holes, "HTYPE"), prep(piles, "PTYPE")
    holes_n = normalize(holes, "co_inactive", "_id", "_name", type_col="_type")
    piles_n = normalize(piles, "co_inactive", "_id", "_name", type_col="_type")

    # A dump logged with its opening's GPS fix is the same site recorded twice.
    piles_n = piles_n[_far_from(piles_n, holes_n, COINCIDENT_M)]
    return pd.concat([holes_n, piles_n], ignore_index=True)


# ---------------------------------------------------------------- build

def _find_cgs(co_dir: Path) -> Path:
    found = sorted([*co_dir.rglob("*.gdb"), *co_dir.rglob("*.gpkg")])
    if not found:
        raise FileNotFoundError(f"no .gdb or .gpkg under {co_dir}")
    return found[0]


def build_labels(raw_dir: str | Path, aoi: AOI, out_path: str | Path, dedupe_m: float = 50.0) -> gpd.GeoDataFrame:
    raw_dir = Path(raw_dir)
    frames = [
        load_usmin(raw_dir / "usmin", scale_dedupe_m=dedupe_m),
        load_co_aml(_find_cgs(raw_dir / "co_aml")),
        load_mrds(raw_dir / "mrds" / "mrds.csv"),
    ]
    merged = merge_labels(frames, aoi=aoi, dedupe_m=dedupe_m)
    # Tag every label with its split now, so nothing downstream can forget to.
    merged["split"] = merged.geometry.map(aoi.split_of)

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    # GPKG can't store tuples; write sources as "mrds,usmin" and convert back on read.
    merged.assign(sources=merged["sources"].map(",".join)).to_file(out_path, layer="labels", driver="GPKG")
    return merged


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--config", default=REPO / "config" / "aoi_ward_goldhill.yaml", type=Path)
    ap.add_argument("--raw", default=REPO / "data" / "raw", type=Path)
    ap.add_argument("--out", default=REPO / "data" / "interim" / "labels.gpkg", type=Path)
    args = ap.parse_args(argv)

    import yaml

    aoi = load_aoi(args.config)
    dedupe_m = float(yaml.safe_load(args.config.read_text()).get("labels", {}).get("dedupe_m", 50))
    merged = build_labels(args.raw, aoi, args.out, dedupe_m=dedupe_m)

    print(f"wrote {len(merged)} labels -> {args.out}")
    print(merged["sources"].map(",".join).value_counts().to_string())
    print(merged["split"].value_counts().to_string())


if __name__ == "__main__":
    main()
