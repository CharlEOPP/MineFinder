"""USGS 3DEP 1 m DEM for the AOI: find tiles (TNM Access API), download, mosaic.

    python -m minefinder.acquire.dem --dry-run   # list tiles + total size, download nothing
    python -m minefinder.acquire.dem             # download to data/raw/dem/, write data/interim/dem_1m.tif

Why it works this way: docs/decisions/0007-dem-acquisition-and-mosaic.md
"""
from __future__ import annotations

import argparse
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Iterable

import numpy as np
import rasterio
from pyproj import CRS
from rasterio.enums import Resampling
from rasterio.warp import transform_bounds
from rasterio.windows import Window, from_bounds

from minefinder.aoi import AOI, PROJECT_CRS, load_aoi

REPO = Path(__file__).resolve().parents[3]

TNM_API = "https://tnmaccess.nationalmap.gov/api/v1/products"
DATASET = "Digital Elevation Model (DEM) 1 meter"
UTM_ZONE = 13
NODATA = -999999.0  # 3DEP's own nodata value for 1 m DEMs

# Horizontal CRSs we accept as "the project grid". NAD83(2011) / UTM 13N differs from
# NAD83 / UTM 13N by about a metre or less here, far below our label accuracy (tens
# of metres), so tiles are relabelled, never resampled. See decision 0007.
ACCEPTED_EPSG = {26913, 6342}

# Origin offsets below this are float noise in the file, not a different grid
# (observed: ~3e-5 m on real tiles). 1 mm is 1/1000 of a pixel.
GRID_TOL_M = 1e-3

BLOCK = 1024  # mosaic block size in pixels; keeps memory to a few MB per tile per block

_TITLE = re.compile(r"1 Meter (\d+) (x\d+y\d+) (\S+)")


@dataclass(frozen=True)
class Tile:
    key: str        # 10 km tile cell, e.g. "x45y443" (same cell can come from several surveys)
    zone: int
    project: str    # lidar project, e.g. "CO_DRCOG_2020_B20"
    date: str       # publication date, ISO yyyy-mm-dd
    url: str
    filename: str
    size: int       # bytes, as reported by the API (0 if unknown)


# ---------------------------------------------------------------- find

def query_bbox_lonlat(aoi: AOI, pad_deg: float = 0.001) -> tuple[float, float, float, float]:
    """AOI extent as a lon/lat bbox for the API. Densified so curved edges stay covered."""
    minx, miny, maxx, maxy = aoi.extent.bounds
    w, s, e, n = transform_bounds(PROJECT_CRS, "EPSG:4326", minx, miny, maxx, maxy, densify_pts=21)
    return (w - pad_deg, s - pad_deg, e + pad_deg, n + pad_deg)


def parse_items(items: Iterable[dict]) -> list[Tile]:
    """API items -> Tiles: GeoTIFF only, our UTM zone only, deduped, newest first."""
    seen, tiles = set(), []
    for it in items:
        url = it.get("downloadURL") or ""
        if not url or it.get("format") != "GeoTIFF" or url in seen:
            continue
        m = _TITLE.search(it.get("title", ""))
        if not m or int(m.group(1)) != UTM_ZONE:
            continue
        seen.add(url)
        tiles.append(
            Tile(
                key=m.group(2),
                zone=int(m.group(1)),
                project=m.group(3),
                date=(it.get("publicationDate") or "")[:10],
                url=url,
                filename=url.rsplit("/", 1)[-1],
                size=int(it.get("sizeInBytes") or 0),
            )
        )
    # Newest survey first; ties broken by project name so the order is deterministic.
    return sorted(tiles, key=lambda t: (t.date, t.project), reverse=True)


def _http_get(params: dict) -> dict:
    import requests

    r = requests.get(TNM_API, params=params, timeout=60)
    r.raise_for_status()
    return r.json()


def list_tiles(aoi: AOI, get: Callable[[dict], dict] | None = None, page_size: int = 100) -> list[Tile]:
    get = get or _http_get
    bbox = ",".join(f"{v:.6f}" for v in query_bbox_lonlat(aoi))
    items, offset, total = [], 0, None
    while total is None or offset < total:
        resp = get({"datasets": DATASET, "bbox": bbox, "prodFormats": "GeoTIFF",
                    "max": page_size, "offset": offset})
        total = int(resp.get("total", 0))
        page = resp.get("items", [])
        items.extend(page)
        if not page:
            break
        offset += page_size
    tiles = parse_items(items)
    if not tiles:
        raise RuntimeError(f"no 1 m DEM tiles returned for bbox {bbox}")
    return tiles


# ---------------------------------------------------------------- download

def _http_stream(url: str):
    import requests

    return requests.get(url, stream=True, timeout=120)


def download_tiles(tiles: list[Tile], dest: str | Path, stream: Callable | None = None) -> list[Path]:
    """Download each tile once. Writes to .part and renames on success, so a file
    that exists under its final name is always complete; re-runs skip it."""
    stream = stream or _http_stream
    dest = Path(dest)
    dest.mkdir(parents=True, exist_ok=True)
    paths = []
    for i, t in enumerate(tiles, 1):
        out = dest / t.filename
        if not out.exists():
            print(f"[{i}/{len(tiles)}] downloading {t.filename} ({t.size / 1e6:.0f} MB)")
            part = out.with_suffix(out.suffix + ".part")
            with stream(t.url) as r:
                r.raise_for_status()
                with open(part, "wb") as f:
                    for chunk in r.iter_content(chunk_size=1 << 20):
                        f.write(chunk)
            part.replace(out)
        paths.append(out)
    return paths


# ---------------------------------------------------------------- mosaic

def _check_tile(src: rasterio.DatasetReader, aoi: AOI) -> None:
    crs = CRS.from_wkt(src.crs.to_wkt())
    horiz = crs.sub_crs_list[0] if crs.is_compound else crs
    epsg = horiz.to_epsg(min_confidence=70)
    if epsg not in ACCEPTED_EPSG:
        raise ValueError(f"{src.name}: CRS {horiz.name} (EPSG:{epsg}) is not EPSG:26913 or its NAD83(2011) twin")
    t = src.transform
    minx, _, _, maxy = aoi.extent.bounds
    def whole_metres(v: float) -> bool:
        # Real 3DEP origins carry ~1e-5 m of float noise from the GeoTIFF tags.
        # Anything under GRID_TOL_M is the same grid; reads snap to whole pixels.
        return abs(v - round(v)) < GRID_TOL_M

    on_grid = (
        np.isclose(t.a, 1.0) and np.isclose(t.e, -1.0) and t.b == 0 and t.d == 0
        and whole_metres(t.c - minx) and whole_metres(maxy - t.f)
    )
    if not on_grid:
        raise ValueError(f"{src.name}: not on the 1 m integer-metre grid; would need resampling")


def mosaic(paths: list[str | Path], aoi: AOI, out_path: str | Path) -> Path:
    """Mosaic tiles onto the AOI's 1 m grid. `paths` must be in priority order
    (newest first): earlier tiles win, later ones only fill nodata gaps.

    Also writes, next to the DEM:
      dem_source.tif  uint8: which tile each pixel came from (1-based index; 0 = no data)
      dem_source.csv  index,filename: so seams between surveys are known exactly (0009)
    """
    if len(paths) > 255:
        raise ValueError("more than 255 tiles: source raster is uint8")
    minx, miny, maxx, maxy = aoi.extent.bounds
    width, height = int(round(maxx - minx)), int(round(maxy - miny))
    transform = rasterio.transform.from_origin(minx, maxy, 1.0, 1.0)

    srcs = [rasterio.open(p) for p in paths]
    try:
        for s in srcs:
            _check_tile(s, aoi)

        out_path = Path(out_path)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        profile = dict(
            driver="GTiff", width=width, height=height, count=1, dtype="float32",
            crs=PROJECT_CRS.to_wkt(), transform=transform, nodata=NODATA,
            tiled=True, blockxsize=512, blockysize=512,
            compress="deflate", predictor=3, BIGTIFF="IF_SAFER",
        )
        src_path = out_path.with_name("dem_source.tif")
        src_profile = {**profile, "dtype": "uint8", "nodata": 0, "predictor": 2}
        with rasterio.open(out_path, "w", **profile) as dst, rasterio.open(src_path, "w", **src_profile) as sdst:
            for row in range(0, height, BLOCK):
                for col in range(0, width, BLOCK):
                    win = Window(col, row, min(BLOCK, width - col), min(BLOCK, height - row))
                    wb = rasterio.windows.bounds(win, transform)
                    block = np.full((int(win.height), int(win.width)), NODATA, dtype="float32")
                    who = np.zeros(block.shape, dtype="uint8")
                    for k, s in enumerate(srcs, 1):
                        if s.bounds.left >= wb[2] or s.bounds.right <= wb[0] or \
                           s.bounds.bottom >= wb[3] or s.bounds.top <= wb[1]:
                            continue
                        sw = from_bounds(*wb, transform=s.transform).round_offsets().round_lengths()
                        z = s.read(1, window=sw, boundless=True, fill_value=NODATA).astype("float32")
                        valid = np.isfinite(z) & (z != NODATA)
                        if s.nodata is not None:
                            valid &= z != np.float32(s.nodata)
                        fill = (block == NODATA) & valid
                        block[fill] = z[fill]
                        who[fill] = k
                        if not (block == NODATA).any():
                            break  # block fully covered by higher-priority tiles
                    dst.write(block, 1, window=win)
                    sdst.write(who, 1, window=win)
    finally:
        for s in srcs:
            s.close()

    out_path.with_name("dem_source.csv").write_text(
        "index,filename\n" + "".join(f"{k},{Path(p).name}\n" for k, p in enumerate(paths, 1))
    )
    # Pyramids: QGIS (and anything else) can draw the big mosaic quickly when zoomed out.
    with rasterio.open(out_path, "r+") as d:
        d.build_overviews([2, 4, 8, 16, 32], Resampling.average)
    with rasterio.open(src_path, "r+") as d:
        d.build_overviews([2, 4, 8, 16, 32], Resampling.nearest)
    return out_path


def coverage(path: str | Path) -> float:
    """Fraction of pixels with data."""
    n = valid = 0
    with rasterio.open(path) as ds:
        for _, w in ds.block_windows(1):
            z = ds.read(1, window=w)
            n += z.size
            valid += int((z != ds.nodata).sum())
    return valid / n if n else 0.0


# ---------------------------------------------------------------- CLI

def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description="Download and mosaic 3DEP 1 m DEM tiles for the AOI")
    ap.add_argument("--config", default=REPO / "config" / "aoi_ward_goldhill.yaml", type=Path)
    ap.add_argument("--tiles-dir", default=REPO / "data" / "raw" / "dem", type=Path)
    ap.add_argument("--out", default=REPO / "data" / "interim" / "dem_1m.tif", type=Path)
    ap.add_argument("--dry-run", action="store_true", help="list tiles and total size, download nothing")
    args = ap.parse_args(argv)

    aoi = load_aoi(args.config)
    tiles = list_tiles(aoi)
    print(f"{len(tiles)} tiles over {len({t.key for t in tiles})} cells, "
          f"{sum(t.size for t in tiles) / 1e9:.2f} GB total (newest first):")
    for t in tiles:
        print(f"  {t.key}  {t.date}  {t.project:<28} {t.size / 1e6:6.0f} MB")
    if args.dry_run:
        return

    paths = download_tiles(tiles, args.tiles_dir)
    out = mosaic(paths, aoi, args.out)
    print(f"wrote {out}  coverage {coverage(out):.1%}")


if __name__ == "__main__":
    main()
