# 0007: Get the DEM from the TNM API, and mosaic newest survey first without ever resampling

**Decision:** `acquire/dem.py` asks USGS's **TNM Access API** for every 3DEP **1 m DEM** GeoTIFF tile touching the AOI. It downloads them to `data/raw/dem/` and mosaics them onto the AOI's exact 1 m grid as `data/interim/dem_1m.tif`. Where surveys overlap, the **newest wins** and older ones only **fill gaps**. Tiles are **never resampled**: any tile that would need it makes the run fail instead.

## In plain words

A **DEM** (digital elevation model) is a grid of ground heights, here one number per square metre. USGS makes them from airborne LiDAR with trees and buildings filtered out ("bare earth"). That's what lets us see a 3 m prospect pit under a forest canopy.

USGS delivers the 1 m DEM as **10 km × 10 km tiles**, named by grid cell (e.g. `x45y443`). Our AOI (21.5 × 19.5 km) touches about 9 cells. Some cells were flown more than once, e.g. a 2013 county survey and a 2020 regional one, so the same cell can come as two files. The **mosaic** stitches them into one seamless raster that covers exactly our AOI.

**Why newest first, with older filling gaps:** newer surveys usually have denser points and better processing, so they win wherever they have data. A survey can stop partway through a tile at the project boundary, though, leaving a hole. Using only the newest tile would leave nodata holes there, and using the oldest would throw away the better data. Layering newest over older gives the best available elevation at every pixel.

**Why never resample:** *resampling* means computing new pixel values on a different grid by interpolating between the old ones. It's unavoidable if you change CRS or shift the grid by a fraction of a pixel, and it smooths the surface slightly every time. A prospect pit is a few pixels across, so smoothing is exactly what makes it disappear. 3DEP tiles already sit on whole-metre UTM 13N coordinates, the same grid as our AOI, so we copy pixels directly. If a tile ever *doesn't* line up, the code stops with an error. It doesn't quietly fix it.

## Formally

- **TNM Access API:** `GET https://tnmaccess.nationalmap.gov/api/v1/products` with `datasets="Digital Elevation Model (DEM) 1 meter"`, `prodFormats=GeoTIFF`, `bbox=minLon,minLat,maxLon,maxLat`, paging with `max` and `offset`. The AOI's bbox is converted to lon/lat with a **densified** boundary transform (21 points per edge), because a UTM rectangle's edges become curves in lon/lat, and corner-only transforms can miss slivers.
- **Tile identity** is parsed from the title `USGS 1 Meter {zone} {xNNyNNN} {project}`. The `key` is the 10 km cell, the `project` is the LiDAR collection, and the `date` is `publicationDate`. Tiles are sorted by `(date, project)` descending, so the order is deterministic.
- **Zone filter:** only UTM zone 13 tiles are kept. A bbox near 108°W can return zone 12 tiles, which are on a different grid.
- **CRS acceptance:** the horizontal part of the tile CRS (3DEP tiles are often compound, *NAD83 / UTM 13N + NAVD88 height*) must be **EPSG:26913** or **EPSG:6342** (NAD83(2011) / UTM 13N). The two are realizations of the same datum about a metre or less apart here, which is below our label accuracy (tens of metres, see 0006). So EPSG:6342 tiles are **relabelled** as 26913 and not resampled. Every other CRS (zone 12, WGS84 UTM, geographic) raises an error. This keeps the "EPSG:26913 everywhere" rule from [0002](0002-crs-epsg-26913.md) without paying for it in blur.
- **Grid alignment check:** pixel size 1 × 1 m, no rotation, and the origin offset from the AOI corner a whole number of metres, **to within 1 mm** (`GRID_TOL_M`). Otherwise raise. The tolerance comes from real data: see *First run* below.
- **Mosaic algorithm (priority fill):** output blocks are 1024 × 1024 px. For each block, read the matching window from each tile in priority order, and copy pixels only where the block is still nodata and the tile has valid data. Stop early once the block is full. Memory stays at a few MB per block regardless of AOI size. A naive `rasterio.merge` of 9+ tiles loads the whole 1.7 GB mosaic into RAM at once.
- **Output:** float32 GeoTIFF, nodata −999999 (3DEP's own value), 512 px internal tiles, DEFLATE with floating-point predictor, BigTIFF if needed. Elevations are metres above **NAVD88**. The vertical datum doesn't matter for terrain shape, since every derivative we compute is relative.
- **Downloads** go to `<name>.part` and are renamed when complete, so a file under its final name is always whole. Re-runs skip finished tiles, which means an interrupted download just picks up where it stopped.

## Why this and not the alternatives

| Option | Why not |
|---|---|
| Click-download from the TNM web viewer | Not reproducible, and easy to miss a tile or grab the wrong project |
| Newest tile only per cell | Leaves holes where a newer survey's boundary cuts through a cell |
| `rasterio.merge` | Loads everything into memory, and its default "first" ordering is easy to get wrong |
| `gdalwarp` to 26913 | Resamples 6342 tiles for a sub-pixel datum difference, which is exactly the blur we're avoiding |
| 3DEP COGs streamed straight from the cloud, or point clouds | Worth doing later. The brief starts with DEMs, and PDAL point-cloud processing comes in a later phase |
| 1/3 arc-second (~10 m) DEM | Far too coarse: a shaft collar disappears at 10 m |

## Where it lives

- `src/minefinder/acquire/dem.py`: `list_tiles`, `parse_items`, `download_tiles`, `mosaic`, `coverage`, the CLI
- Tests: `tests/acquire/test_dem.py`. They fake the API and the downloads, and use synthetic GeoTIFFs to cover newest-wins/older-fills, uncovered area as nodata, 6342 accepted, zone 12 / WGS84 UTM / geographic refused, and a half-metre misalignment refused

**Unverified until your first run:** the API's live response format couldn't be checked from Claude's environment, because the API host isn't reachable from there. The parser follows the documented fields (`title`, `downloadURL`, `publicationDate`, `format`, `sizeInBytes`). Run `--dry-run` first. If it lists tiles with sensible dates and sizes, the format is confirmed.

## First run (2026-10-25)

The live API returned **20 tiles over 9 cells, 3.95 GB**, from four surveys. That confirmed the parser against the real response:

| Project | Published | Notes |
|---|---|---|
| UT_COUTAML_D24 | 2026-06-30 | 1 cell, partial (~21% of its tile) |
| CO_ArapahoRooseveltPikeNF_D23 | 2026-03-20 | 9 cells, but some tiles are nearly empty (~0–1% valid at the west edge) |
| CO_CameronPeakWildfire_2021_D21 | 2022-11-03 | 1 cell, partial (~10%) |
| CO_DRCOG_2020_B20 | 2022-02-10 | all 9 cells, full coverage |

**The "newest wins, older fills gaps" rule matters here.** If we used only the newest survey per cell, much of the AOI would be holes. The full-coverage 2020 DRCOG survey fills them.

**The first mosaic attempt failed the grid check, correctly but for the wrong reason.** Real tile origins look like `459993.99996660`, which is about 30 µm off a whole metre. That's float noise in the GeoTIFF tags, not a different grid. The original tolerance (1 µm) was too strict, so it's now 1 mm. That's 1/1000 of a pixel, and a test proves that a 1 cm offset is still refused. All four inspected tiles are `NAD83 / UTM zone 13N` (EPSG:26913), 10,012 × 10,012 px, with a 6 m overlap on each side.

## The Phase 1 "done when" check (QGIS)

1. Open `data/interim/dem_1m.tif`. Under **Layer Properties → Symbology**, set the render type to **Hillshade** and turn on **Multidirectional**.
2. Add `data/interim/labels.gpkg` and categorize it by `sources`.
3. Zoom in on a few USMIN adits and CGS dumps. The point should sit on or right next to a visible notch, bench or fan. MRDS-only points will often sit on blank slopes, which is expected (see 0006).
4. Make sure the project CRS shows **EPSG:26913**, so QGIS isn't reprojecting on the fly.

If USMIN and CGS points line up with landforms, Phase 1 is done.

## QGIS check findings (2026-10-26)

- **A grid lattice appeared in the hillshade. It's a QGIS display artifact, not data.** With nearest-neighbour resampling slightly past 1:1 zoom, QGIS repeats a row or column every so often, and the hillshade draws each repeat as a straight line. The fix is Symbology → Resampling → **Bilinear** for both zoomed in and zoomed out. Hillshades rendered directly from the raw 2026 and 2020 tiles show no lattice.
- **Surveys agree:** in cell x45y444, the 2026 ArapahoRoosevelt survey minus the 2020 DRCOG survey has a median difference of **+0.10 m** (IQR 0.06–0.13 m).
- **Survey patch edges are visible as straight lines in the hillshade** where a partial newer survey meets an older fill, e.g. the UT_COUTAML_D24 edge in x46y444. The step is about 0.1 m, but a local-relief model will register it as a linear feature. **Phase 2 must handle it**, either by feathering the overlap at patch edges or by masking a few metres either side of seams during detection.
- The 2020 DRCOG data has faint horizontal flight-line streaks. They're normal for LiDAR. Watch for them as false positives in Phase 3.
- **TODO (dem.py):** build overviews (pyramids) at mosaic time, so QGIS draws the 645 MB DEM quickly when zoomed out.

## Revisit when

- **Phase 2+:** switch to the LiDAR point clouds (LAZ/COPC with PDAL) to make our own ground model with finer filtering.
- **The statewide run:** western Colorado is UTM zone 12, which needs either per-zone processing or a CRS decision record (see 0002).
- **Survey age matters:** some mines were reclaimed between surveys. Then keep per-survey DEMs instead of one mosaic, and record which survey each pixel came from.
