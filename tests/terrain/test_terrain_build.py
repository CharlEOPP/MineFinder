"""Whole-raster derivatives and the blockwise builder."""
import numpy as np
import pytest
import rasterio
from rasterio.transform import from_origin

from minefinder.terrain.derivatives import (
    LRM_LARGE,
    LRM_SMALL,
    build_terrain,
    curvature,
    hillshade,
    local_relief,
    seam_mask,
    slope_deg,
)

X0, Y0 = 450_000, 4_430_000


def test_slope_of_a_45_degree_plane():
    yy, xx = np.mgrid[0:50, 0:50].astype("float64")
    s = slope_deg(xx * 1.0)                      # rises 1 m per 1 m east
    assert np.allclose(s[5:-5, 5:-5], 45.0, atol=1e-6)


def test_slope_of_flat_ground_is_zero():
    assert np.allclose(slope_deg(np.full((20, 20), 2500.0))[2:-2, 2:-2], 0.0)


def test_hillshade_is_uint8_and_flat_ground_is_uniform():
    h = hillshade(np.full((30, 30), 2500.0))
    assert h.dtype == np.uint8
    assert len(np.unique(h[2:-2, 2:-2])) == 1 and h[10, 10] > 0


def test_hillshade_marks_nodata_as_zero():
    z = np.full((30, 30), 2500.0)
    z[10, 10] = np.nan
    assert hillshade(z)[10, 10] == 0


def test_curvature_of_a_bowl_is_its_laplacian():
    # z = a(x^2 + y^2) has Laplacian 4a everywhere: positive = concave-up (pit floor).
    yy, xx = np.mgrid[-20:21, -20:21].astype("float64")
    c = curvature(0.05 * (xx**2 + yy**2))
    assert np.allclose(c[2:-2, 2:-2], 0.2, atol=1e-9)


def test_seam_mask_marks_both_sides_of_a_source_change():
    src = np.ones((40, 40), dtype="uint8")
    src[:, 20:] = 2
    m = seam_mask(src, radius=3)
    assert m[:, 17:23].all()                     # 3 px each side of the 19|20 boundary
    # Away from the seam and from the raster edge (also masked: LRM is truncated there)
    assert not m[4:-4, 4:16].any() and not m[4:-4, 24:-4].any()


def test_seam_mask_also_marks_nodata_edges():
    # LRM windows are truncated next to nodata, so those pixels are unreliable too.
    src = np.ones((40, 40), dtype="uint8")
    src[:, 30:] = 0
    m = seam_mask(src, radius=2)
    assert m[:, 28:30].all() and not m[3:-3, 3:27].any()


# ---------------- blockwise builder ----------------

def write(path, arr, dtype="float32", nodata=-999999.0):
    with rasterio.open(path, "w", driver="GTiff", width=arr.shape[1], height=arr.shape[0], count=1,
                       dtype=dtype, crs="EPSG:26913", transform=from_origin(X0, Y0, 1, 1), nodata=nodata) as d:
        d.write(arr.astype(dtype), 1)
    return path


@pytest.fixture
def terrain_inputs(tmp_path):
    rng = np.random.default_rng(0)
    yy, xx = np.mgrid[0:300, 0:260].astype("float64")
    z = 2600 + 0.3 * xx - 0.2 * yy + 3 * np.sin(xx / 17) * np.cos(yy / 23) + rng.normal(0, 0.05, xx.shape)
    z[200:220, 40:60] = -999999.0                       # a nodata hole
    src = np.ones(z.shape, dtype="uint8")
    src[:, 150:] = 2                                     # a survey seam at column 150
    src[200:220, 40:60] = 0
    return write(tmp_path / "dem.tif", z), write(tmp_path / "dem_source.tif", src, "uint8", 0), z


def test_build_terrain_blockwise_equals_whole_array(tmp_path, terrain_inputs):
    # The key correctness test: small blocks + halo must give identical results to
    # computing on the whole raster at once (no seams at block edges).
    dem, src, _ = terrain_inputs
    out = build_terrain(dem, src, tmp_path / "terrain", block=64)
    with rasterio.open(dem) as d:                 # reference from the same float32 values on disk
        z = d.read(1).astype("float64")
    whole = np.where(z == -999999.0, np.nan, z)
    with rasterio.open(out["lrm_small"]) as d:
        got = d.read(1).astype("float64") * d.scales[0]
        nd = d.nodata
    ref = local_relief(whole, LRM_SMALL)
    valid = np.isfinite(ref)
    assert np.allclose(got[valid], np.round(ref[valid] * 100) / 100, atol=0.0051)
    with rasterio.open(out["lrm_large"]) as d:
        g2 = d.read(1).astype("float64") * d.scales[0]
    ref2 = local_relief(whole, LRM_LARGE)
    assert np.allclose(g2[np.isfinite(ref2)], ref2[np.isfinite(ref2)], atol=0.0051)


def test_build_terrain_outputs_are_26913_with_overviews(tmp_path, terrain_inputs):
    dem, src, _ = terrain_inputs
    out = build_terrain(dem, src, tmp_path / "terrain", block=64)
    assert set(out) == {"hillshade", "slope", "curvature", "lrm_small", "lrm_large", "seam_mask"}
    for p in out.values():
        with rasterio.open(p) as d:
            assert d.crs.to_epsg() == 26913
            assert d.transform == from_origin(X0, Y0, 1, 1)
            assert d.overviews(1)


def test_build_terrain_seam_mask_covers_survey_seam_and_hole(tmp_path, terrain_inputs):
    dem, src, _ = terrain_inputs
    out = build_terrain(dem, src, tmp_path / "terrain", block=64)
    with rasterio.open(out["seam_mask"]) as d:
        m = d.read(1).astype(bool)
    r = LRM_SMALL // 2
    assert m[:, 150 - r : 150 + r].all()              # survey seam
    assert m[200:220, 40:60].all()                    # nodata hole itself
    assert not m[50:150, 60:120].any()                # clean interior


def test_mask_fraction_weights_by_pixels_not_blocks(tmp_path):
    # 600 x 520 raster with 512 px internal tiles: the edge tiles are small. The old CLI
    # averaged per-block means, so a fully masked thin edge strip counted as much as a
    # full 512x512 block and overstated coverage (4.3% printed vs 2.7% real).
    from minefinder.terrain.derivatives import mask_fraction
    m = np.zeros((600, 520), dtype="uint8")
    m[512:, :] = 1                                    # the 88-row strip of edge tiles
    p = tmp_path / "m.tif"
    with rasterio.open(p, "w", driver="GTiff", width=520, height=600, count=1, dtype="uint8",
                       crs="EPSG:26913", transform=from_origin(X0, Y0, 1, 1),
                       tiled=True, blockxsize=512, blockysize=512) as d:
        d.write(m, 1)
    assert mask_fraction(p) == pytest.approx(88 / 600)
