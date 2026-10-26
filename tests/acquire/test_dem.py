"""DEM acquisition tests. No network: the TNM API and downloads are faked, and
tiles are small synthetic GeoTIFFs on the real 1 m grid."""
import json
from pathlib import Path

import numpy as np
import pytest
import rasterio
from pyproj import Transformer
from rasterio.transform import from_origin

from minefinder.acquire.dem import (
    NODATA,
    Tile,
    download_tiles,
    list_tiles,
    mosaic,
    parse_items,
    query_bbox_lonlat,
)
from tests.conftest import X0, Y0


def item(title, url, date, fmt="GeoTIFF", size=100):
    return {
        "title": title,
        "downloadURL": url,
        "publicationDate": date,
        "format": fmt,
        "sizeInBytes": size,
        "boundingBox": {"minX": -105.6, "maxX": -105.5, "minY": 40.0, "maxY": 40.1},
    }


T_NEW = item("USGS 1 Meter 13 x45y443 CO_DRCOG_2020_B20", "https://x/USGS_1M_13_x45y443_CO_DRCOG_2020_B20.tif", "2022-02-10")
T_OLD = item("USGS 1 Meter 13 x45y443 CO_Boulder_2013", "https://x/USGS_1M_13_x45y443_CO_Boulder_2013.tif", "2015-06-01")


# ---------------- query ----------------

def test_query_bbox_covers_every_aoi_corner(aoi):
    minlon, minlat, maxlon, maxlat = query_bbox_lonlat(aoi)
    to_ll = Transformer.from_crs(26913, 4326, always_xy=True)
    for x, y in aoi.extent.exterior.coords:
        lon, lat = to_ll.transform(x, y)
        assert minlon <= lon <= maxlon and minlat <= lat <= maxlat


# ---------------- parse ----------------

def test_parse_items_extracts_tile_key_zone_project_and_date():
    (t,) = parse_items([T_NEW])
    assert t.key == "x45y443"
    assert t.zone == 13
    assert t.project == "CO_DRCOG_2020_B20"
    assert t.date == "2022-02-10"
    assert t.filename == "USGS_1M_13_x45y443_CO_DRCOG_2020_B20.tif"


def test_parse_items_skips_non_geotiff_and_missing_urls():
    bad_fmt = item("USGS 1 Meter 13 x45y444 P", "https://x/a.img", "2020-01-01", fmt="IMG")
    no_url = item("USGS 1 Meter 13 x45y445 P", "", "2020-01-01")
    assert [t.key for t in parse_items([T_NEW, bad_fmt, no_url])] == ["x45y443"]


def test_parse_items_drops_other_utm_zones():
    # A bbox near a zone edge can return zone-12 tiles; they're on a different grid.
    z12 = item("USGS 1 Meter 12 x75y443 P", "https://x/USGS_1M_12_x75y443_P.tif", "2020-01-01")
    assert [t.zone for t in parse_items([T_NEW, z12])] == [13]


def test_parse_items_dedupes_repeated_urls():
    assert len(parse_items([T_NEW, dict(T_NEW)])) == 1


def test_parse_items_orders_newest_first():
    # Order matters: the mosaic lets earlier tiles win, older ones only fill gaps.
    assert [t.project for t in parse_items([T_OLD, T_NEW])] == ["CO_DRCOG_2020_B20", "CO_Boulder_2013"]


# ---------------- API paging ----------------

def test_list_tiles_pages_through_results(aoi):
    pages = {0: {"total": 2, "items": [T_NEW]}, 1: {"total": 2, "items": [T_OLD]}}
    calls = []

    def fake_get(params):
        calls.append(params)
        return pages[int(params["offset"])]

    tiles = list_tiles(aoi, get=fake_get, page_size=1)
    assert len(tiles) == 2 and len(calls) == 2
    assert calls[0]["datasets"] == "Digital Elevation Model (DEM) 1 meter"
    assert calls[0]["prodFormats"] == "GeoTIFF"


def test_list_tiles_raises_when_api_returns_nothing(aoi):
    with pytest.raises(RuntimeError, match="no 1 m DEM tiles"):
        list_tiles(aoi, get=lambda p: {"total": 0, "items": []})


# ---------------- download ----------------

class FakeResp:
    def __init__(self, data):
        self.data = data
        self.headers = {"content-length": str(len(data))}

    def raise_for_status(self):
        pass

    def iter_content(self, chunk_size):
        yield self.data

    def __enter__(self):
        return self

    def __exit__(self, *a):
        pass


def test_download_writes_files_and_skips_complete_ones(tmp_path):
    (t,) = parse_items([T_NEW])
    got = []

    def fake_stream(url):
        got.append(url)
        return FakeResp(b"tifbytes")

    paths = download_tiles([t], tmp_path, stream=fake_stream)
    assert paths[0].read_bytes() == b"tifbytes"
    assert not list(tmp_path.glob("*.part"))       # no half-written leftovers
    download_tiles([t], tmp_path, stream=fake_stream)
    assert len(got) == 1                           # second run skipped the finished file


# ---------------- mosaic ----------------

def write_tile(path: Path, x0, y0, w, h, value, crs="EPSG:26913", hole=None, res=1.0):
    """GeoTIFF with upper-left (x0, y0); `hole` = (row slice, col slice) set to nodata."""
    data = np.full((h, w), value, dtype="float32")
    if hole:
        data[hole] = NODATA
    with rasterio.open(
        path, "w", driver="GTiff", width=w, height=h, count=1, dtype="float32",
        crs=crs, transform=from_origin(x0, y0, res, res), nodata=NODATA,
    ) as dst:
        dst.write(data, 1)
    return path


@pytest.fixture
def small_aoi(tmp_path):
    """A 100 x 100 m AOI so mosaics stay tiny."""
    from minefinder.aoi import load_aoi
    from tests.conftest import write_yaml

    d = {
        "name": "tiny", "crs": "EPSG:26913",
        "extent": {"minx": X0, "miny": Y0, "maxx": X0 + 100, "maxy": Y0 + 100},
        "holdout": {"minx": X0, "miny": Y0, "maxx": X0 + 20, "maxy": Y0 + 20},
        "holdout_buffer_m": 5,
    }
    return load_aoi(write_yaml(tmp_path, d, "tiny.yaml"))


def test_mosaic_is_26913_on_the_aoi_grid(tmp_path, small_aoi):
    a = write_tile(tmp_path / "a.tif", X0 - 10, Y0 + 110, 130, 130, 2000.0)
    out = mosaic([a], small_aoi, tmp_path / "dem.tif")
    with rasterio.open(out) as ds:
        assert ds.crs.to_epsg() == 26913
        assert tuple(ds.bounds) == (X0, Y0, X0 + 100, Y0 + 100)
        assert ds.res == (1.0, 1.0)
        assert (ds.read(1) == 2000.0).all()


def test_mosaic_newest_wins_and_older_fills_gaps(tmp_path, small_aoi):
    # Newer tile has a 10x10 nodata hole; the older survey fills it.
    new = write_tile(tmp_path / "new.tif", X0, Y0 + 100, 100, 100, 2000.0, hole=(slice(0, 10), slice(0, 10)))
    old = write_tile(tmp_path / "old.tif", X0, Y0 + 100, 100, 100, 1999.0)
    with rasterio.open(mosaic([new, old], small_aoi, tmp_path / "dem.tif")) as ds:
        z = ds.read(1)
    assert (z[:10, :10] == 1999.0).all()
    assert (z[10:, :] == 2000.0).all()


def test_mosaic_leaves_uncovered_area_as_nodata(tmp_path, small_aoi):
    half = write_tile(tmp_path / "half.tif", X0, Y0 + 100, 50, 100, 2000.0)  # covers west half only
    with rasterio.open(mosaic([half], small_aoi, tmp_path / "dem.tif")) as ds:
        z = ds.read(1)
        assert ds.nodata == NODATA
    assert (z[:, 50:] == NODATA).all()


def test_mosaic_accepts_nad83_2011_without_resampling(tmp_path, small_aoi):
    # Newer 3DEP projects are tagged NAD83(2011) / UTM 13N (EPSG:6342). It's the same
    # grid to well under a pixel, so we relabel instead of resampling (see decision 0007).
    a = write_tile(tmp_path / "a.tif", X0, Y0 + 100, 100, 100, 2000.0, crs="EPSG:6342")
    with rasterio.open(mosaic([a], small_aoi, tmp_path / "dem.tif")) as ds:
        assert ds.crs.to_epsg() == 26913
        assert (ds.read(1) == 2000.0).all()


@pytest.mark.parametrize("bad_crs", ["EPSG:26912", "EPSG:32613", "EPSG:4326"])
def test_mosaic_refuses_other_crs(tmp_path, small_aoi, bad_crs):
    a = write_tile(tmp_path / "a.tif", X0, Y0 + 100, 100, 100, 2000.0, crs=bad_crs)
    with pytest.raises(ValueError, match="26913"):
        mosaic([a], small_aoi, tmp_path / "dem.tif")


def test_mosaic_snaps_float_noise_in_tile_origin(tmp_path, small_aoi):
    # Real 3DEP tiles store origins like 459993.99996660 (~30 micrometres off a whole
    # metre): float noise in the GeoTIFF tags, not a different grid. Snap, don't refuse.
    a = write_tile(tmp_path / "a.tif", X0 - 3.3e-5, Y0 + 100 + 4.1e-5, 100, 100, 2000.0)
    with rasterio.open(mosaic([a], small_aoi, tmp_path / "dem.tif")) as ds:
        assert (ds.read(1) == 2000.0).all()


def test_mosaic_refuses_centimetre_offset(tmp_path, small_aoi):
    a = write_tile(tmp_path / "a.tif", X0 + 0.01, Y0 + 100, 100, 100, 2000.0)
    with pytest.raises(ValueError, match="grid"):
        mosaic([a], small_aoi, tmp_path / "dem.tif")


def test_mosaic_refuses_misaligned_grid(tmp_path, small_aoi):
    # Half-metre offset would need resampling, which blurs small pits. Fail loudly instead.
    a = write_tile(tmp_path / "a.tif", X0 + 0.5, Y0 + 100.5, 100, 100, 2000.0)
    with pytest.raises(ValueError, match="grid"):
        mosaic([a], small_aoi, tmp_path / "dem.tif")


# ---------------- provenance + overviews (Phase 2) ----------------

def test_mosaic_writes_source_raster_naming_the_winning_tile(tmp_path, small_aoi):
    # Which survey each pixel came from: seams are then known exactly (decision 0009).
    new = write_tile(tmp_path / "new.tif", X0, Y0 + 100, 100, 100, 2000.0, hole=(slice(0, 10), slice(0, 10)))
    old = write_tile(tmp_path / "old.tif", X0, Y0 + 100, 100, 100, 1999.0)
    out = mosaic([new, old], small_aoi, tmp_path / "dem.tif")
    src_path = out.with_name("dem_source.tif")
    with rasterio.open(src_path) as s:
        src = s.read(1)
        assert s.dtypes[0] == "uint8"
    assert (src[:10, :10] == 2).all()          # 1-based index into the tile list: old.tif
    assert (src[10:, :] == 1).all()             # new.tif
    table = (out.with_name("dem_source.csv")).read_text()
    assert "1,new.tif" in table and "2,old.tif" in table


def test_mosaic_source_is_zero_where_no_data(tmp_path, small_aoi):
    half = write_tile(tmp_path / "half.tif", X0, Y0 + 100, 50, 100, 2000.0)
    out = mosaic([half], small_aoi, tmp_path / "dem.tif")
    with rasterio.open(out.with_name("dem_source.tif")) as s:
        assert (s.read(1)[:, 50:] == 0).all()


def test_mosaic_builds_overviews(tmp_path, small_aoi):
    a = write_tile(tmp_path / "a.tif", X0, Y0 + 100, 100, 100, 2000.0)
    with rasterio.open(mosaic([a], small_aoi, tmp_path / "dem.tif")) as ds:
        assert ds.overviews(1)                  # pyramids, so QGIS draws fast zoomed out
