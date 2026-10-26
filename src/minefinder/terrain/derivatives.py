"""Terrain derivatives on the 1 m DEM, plus a blockwise builder for the whole AOI.

    python -m minefinder.terrain.derivatives     # data/interim/dem_1m.tif -> data/processed/terrain/

Outputs (all EPSG:26913, same 1 m grid as the DEM, with overviews):
  hillshade.tif   uint8   multidirectional hillshade (0 = nodata)
  slope.tif       uint8   slope in whole degrees (255 = nodata)
  curvature.tif   float32 Laplacian curvature, 1/m (+ = concave-up, e.g. pit floors)
  lrm_small.tif   int16   local relief, 21 m window, centimetres  (shafts, prospect pits)
  lrm_large.tif   int16   local relief, 61 m window, centimetres  (dumps, adit benches)
  seam_mask.tif   uint8   1 where derivatives are unreliable: survey seams, nodata edges

Why: docs/decisions/0009-terrain-derivatives-and-seams.md  (alignment helpers: 0008)
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import rasterio
from rasterio.enums import Resampling
from rasterio.windows import Window

REPO = Path(__file__).resolve().parents[3]

LRM_SMALL = 21      # m: ~2-3x a shaft/prospect-pit diameter (3-10 m)
LRM_LARGE = 61      # m: ~2x a large dump or adit bench (10-30 m); still << a hillside
HALO = LRM_LARGE // 2 + 2
OVERVIEWS = [2, 4, 8, 16, 32]
LRM_SCALE = 0.01    # int16 centimetres: +/-327 m range, 1 cm precision
LRM_NODATA = -32768


# ---------------------------------------------------------------- per-array

def local_relief(z: np.ndarray, window: int = 31) -> np.ndarray:
    """Elevation minus its `window` x `window` moving average (NaN-aware).

    Removes regional slope so small landforms stand out: pits < 0, mounds > 0.
    A plane maps to exactly 0, because a box mean of a plane equals its centre.
    Near edges and nodata, the mean uses only the valid pixels in the window.
    """
    if window % 2 != 1:
        raise ValueError(f"window must be odd, got {window}")
    z = np.asarray(z, dtype="float64")
    valid = np.isfinite(z)
    h = window // 2

    def box_sum(a: np.ndarray) -> np.ndarray:
        # Integral image: any window sum in O(1). Zero padding truncates at edges.
        p = np.pad(a, ((h + 1, h), (h + 1, h)))
        c = p.cumsum(0).cumsum(1)
        return c[window:, window:] - c[:-window, window:] - c[window:, :-window] + c[:-window, :-window]

    s = box_sum(np.where(valid, z, 0.0))
    n = box_sum(valid.astype("float64"))
    with np.errstate(invalid="ignore", divide="ignore"):
        out = z - s / n
    out[~valid] = np.nan
    return out


def pit_depth(lrm: np.ndarray, radius: int = 4) -> np.ndarray:
    """Depth of the deepest local-relief low within `radius` px (0 if none).

    Shafts and prospect pits are a few metres across, so a 4 m disk tolerates a label
    sitting on the pit rim rather than dead centre.
    """
    h, w = lrm.shape
    p = np.pad(np.where(np.isfinite(lrm), lrm, np.inf), radius, constant_values=np.inf)
    lows = np.full(lrm.shape, np.inf)
    # Min filter over a disk: one shifted view per disk offset (49 for r=4).
    for dy in range(-radius, radius + 1):
        for dx in range(-radius, radius + 1):
            if dx * dx + dy * dy <= radius * radius:
                np.minimum(lows, p[radius + dy : radius + dy + h, radius + dx : radius + dx + w], out=lows)
    return np.maximum(0.0, -np.where(np.isfinite(lows), lows, 0.0))


def _gradients(z: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """dz/dx (east) and dz/dy (north) in m/m on a 1 m grid. Rows run south, hence the sign."""
    d_row, d_col = np.gradient(np.asarray(z, dtype="float64"))
    return d_col, -d_row


def slope_deg(z: np.ndarray) -> np.ndarray:
    gx, gy = _gradients(z)
    return np.degrees(np.arctan(np.hypot(gx, gy)))


def hillshade(z: np.ndarray, azimuths=(225, 270, 315, 360), altitude: float = 45.0) -> np.ndarray:
    """Multidirectional hillshade as uint8 (1-255 lit, 0 = nodata).

    Average of several light directions, so features facing any way stay visible
    (a single sun angle hides linear features parallel to it).
    """
    gx, gy = _gradients(z)
    norm = np.sqrt(gx**2 + gy**2 + 1.0)
    nx, ny, nz = -gx / norm, -gy / norm, 1.0 / norm          # unit surface normal (east, north, up)
    alt = np.radians(altitude)
    shade = np.zeros_like(nz)
    for az in azimuths:
        a = np.radians(az)
        lx, ly, lz = np.sin(a) * np.cos(alt), np.cos(a) * np.cos(alt), np.sin(alt)
        shade += np.clip(nx * lx + ny * ly + nz * lz, 0, 1)
    shade /= len(azimuths)
    # A nodata pixel still gets a finite gradient from its neighbours, so mask it too.
    bad = ~np.isfinite(shade) | ~np.isfinite(np.asarray(z, dtype="float64"))
    out = (1 + np.round(np.where(bad, 0.0, shade) * 254)).astype(np.uint8)
    out[bad] = 0
    return out


def curvature(z: np.ndarray) -> np.ndarray:
    """Laplacian (5-point stencil), 1/m. Positive = concave-up (pit floors, gully bottoms)."""
    z = np.asarray(z, dtype="float64")
    p = np.pad(z, 1, constant_values=np.nan)
    return p[:-2, 1:-1] + p[2:, 1:-1] + p[1:-1, :-2] + p[1:-1, 2:] - 4 * z


def _square_filter(a: np.ndarray, radius: int, fn, fill) -> np.ndarray:
    """Separable square min/max filter (2r+1 window); outside the array counts as `fill`."""
    out = a
    for axis in (0, 1):
        pad = [(0, 0), (0, 0)]
        pad[axis] = (radius, radius)
        p = np.pad(out, pad, constant_values=fill)
        n = out.shape[axis]
        acc = np.take(p, range(0, n), axis=axis)
        for k in range(1, 2 * radius + 1):
            acc = fn(acc, np.take(p, range(k, k + n), axis=axis))
        out = acc
    return out


def seam_mask(source: np.ndarray, radius: int) -> np.ndarray:
    """True within `radius` px of any change in source survey, of nodata (source 0), or
    of the raster edge. Local relief there mixes surveys or is truncated, so a detector
    would see straight-line artifacts that aren't terrain."""
    src = np.asarray(source, dtype=np.int16)
    hi = _square_filter(src, radius, np.maximum, 0)
    lo = _square_filter(src, radius, np.minimum, 0)
    return (hi != lo) | (src == 0)


# ---------------------------------------------------------------- blockwise build

def _profile(ref, dtype, nodata, scale=None):
    p = dict(driver="GTiff", width=ref.width, height=ref.height, count=1, dtype=dtype,
             crs=ref.crs, transform=ref.transform, nodata=nodata, tiled=True,
             blockxsize=512, blockysize=512, compress="deflate", BIGTIFF="IF_SAFER")
    if dtype.startswith("float"):
        p["predictor"] = 3
    elif dtype in ("int16", "uint16"):
        p["predictor"] = 2
    return p


def build_terrain(dem_path: str | Path, source_path: str | Path, out_dir: str | Path,
                  block: int = 2048) -> dict[str, Path]:
    """Compute every derivative block by block. Each block is read with a HALO of
    surrounding pixels (boundless; outside the raster = nodata), so results are
    identical to computing on the whole raster at once."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    names = ["hillshade", "slope", "curvature", "lrm_small", "lrm_large", "seam_mask"]
    paths = {n: out_dir / f"{n}.tif" for n in names}

    with rasterio.open(dem_path) as dem, rasterio.open(source_path) as src:
        if dem.crs.to_epsg() != 26913 or src.transform != dem.transform or src.shape != dem.shape:
            raise ValueError("DEM must be EPSG:26913 and the source raster must share its grid")
        specs = {
            "hillshade": _profile(dem, "uint8", 0),
            "slope": _profile(dem, "uint8", 255),
            "curvature": _profile(dem, "float32", -9999.0),
            "lrm_small": _profile(dem, "int16", LRM_NODATA),
            "lrm_large": _profile(dem, "int16", LRM_NODATA),
            "seam_mask": _profile(dem, "uint8", None),
        }
        outs = {n: rasterio.open(paths[n], "w", **specs[n]) for n in names}
        try:
            for n in ("lrm_small", "lrm_large"):
                outs[n].scales = (LRM_SCALE,)
                outs[n].units = ("m",)
            H, W = dem.height, dem.width
            for row in range(0, H, block):
                for col in range(0, W, block):
                    h, w = min(block, H - row), min(block, W - col)
                    big = Window(col - HALO, row - HALO, w + 2 * HALO, h + 2 * HALO)
                    z = dem.read(1, window=big, boundless=True, fill_value=dem.nodata).astype("float64")
                    z[z == dem.nodata] = np.nan
                    s = src.read(1, window=big, boundless=True, fill_value=0)
                    core = (slice(HALO, HALO + h), slice(HALO, HALO + w))
                    valid = np.isfinite(z[core])
                    win = Window(col, row, w, h)

                    outs["hillshade"].write(hillshade(z)[core], 1, window=win)
                    sl = slope_deg(z)[core]
                    ok = valid & np.isfinite(sl)          # slope is NaN beside nodata too
                    outs["slope"].write(np.where(ok, np.clip(np.round(np.nan_to_num(sl)), 0, 90), 255).astype("uint8"), 1, window=win)
                    cv = curvature(z)[core]
                    outs["curvature"].write(np.where(np.isfinite(cv), cv, -9999.0).astype("float32"), 1, window=win)
                    for n, wsz in (("lrm_small", LRM_SMALL), ("lrm_large", LRM_LARGE)):
                        l = local_relief(z, wsz)[core]
                        q = np.where(np.isfinite(l), np.round(l / LRM_SCALE), LRM_NODATA)
                        outs[n].write(np.clip(q, -32767, 32767).astype("int16"), 1, window=win)
                    outs["seam_mask"].write(seam_mask(s, LRM_LARGE // 2)[core].astype("uint8"), 1, window=win)
        finally:
            for o in outs.values():
                o.close()

    for n, p in paths.items():
        with rasterio.open(p, "r+") as d:
            d.build_overviews(OVERVIEWS, Resampling.nearest if n == "seam_mask" else Resampling.average)
    return paths


def mask_fraction(path: str | Path) -> float:
    """Share of pixels set in a 0/1 mask, weighted by pixel count (edge blocks are smaller)."""
    set_px = total = 0
    with rasterio.open(path) as d:
        for _, w in d.block_windows(1):
            a = d.read(1, window=w)
            set_px += int(a.astype(bool).sum())
            total += a.size
    return set_px / total if total else 0.0


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description="Build terrain derivatives for the AOI")
    ap.add_argument("--dem", default=REPO / "data" / "interim" / "dem_1m.tif", type=Path)
    ap.add_argument("--source", default=REPO / "data" / "interim" / "dem_source.tif", type=Path)
    ap.add_argument("--out", default=REPO / "data" / "processed" / "terrain", type=Path)
    args = ap.parse_args(argv)
    if not args.source.exists():
        raise SystemExit(f"{args.source} missing: re-run `python -m minefinder.acquire.dem` to create it")
    paths = build_terrain(args.dem, args.source, args.out)
    for p in paths.values():
        print(f"wrote {p}")
    print(f"seam/edge mask covers {mask_fraction(paths['seam_mask']):.1%} of the AOI")


if __name__ == "__main__":
    main()
