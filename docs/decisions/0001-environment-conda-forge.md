# 0001: Use a conda-forge environment, installed with Miniforge

**Decision:** Python dependencies come from the **conda-forge** channel through Miniforge. `pip` installs only our own package (`pip install -e ".[dev]"`).

## In plain words

Most Python libraries are pure Python, and `pip` handles them fine. The geospatial stack is different: GDAL (which reads rasters), PROJ (which does coordinate math) and PDAL (which reads LiDAR point clouds) are large C/C++ programs that every Python geo library sits on top of. If rasterio is built against one copy of GDAL and PDAL against another, things break in confusing ways. On Windows this is the classic "DLL load failed" error.

conda-forge builds all of these together, against one shared copy of each C library, so they agree with each other. Miniforge is the small installer that sets up conda with conda-forge as the default, and includes `mamba`, a faster solver.

## Formally

- **Binary dependency consistency.** `pip` wheels each bundle their own copy of native libraries, whereas conda packages share them. Two pip wheels that each bundle a different GDAL ABI can't safely be loaded in the same process.
- **conda-forge** is a community channel with a single build matrix, so every package in a release is compiled against the same `libgdal`, `proj` and `geos`.
- `environment.yml` is the declarative spec. `mamba env create -f environment.yml` reproduces the environment.

## Why this and not the alternatives

| Option | Why not |
|---|---|
| pip + venv only | PDAL has no reliable pip wheel on Windows, and GDAL wheels conflict with each other |
| Anaconda (full) | Huge, defaults to the `defaults` channel, and mixing channels is the most common way geo envs break |
| Docker | Strong reproducibility, but slow file I/O on Windows and overkill for a prototype; revisit for Phase 5+ |
| OSGeo4W | Fine for QGIS, awkward for a Python project env |

## Choices made along the way

- **We started with a lighter env.** Phase 0/1 only needs `geopandas pyproj shapely pyyaml pytest`. PyTorch and the raster stack get added when Phase 2 needs them, so there's no 3 GB download before it's useful.
- **QGIS gets installed separately** with its standalone Windows installer, not from conda. It's a desktop app, and its conda package drags in a lot of dependencies and tends to conflict with the rest.
- **Activate with `conda activate`, install with `mamba install`.** `mamba activate` needs extra shell setup that Miniforge's `conda init` already covers.

## Update (Phase 1, DEM)

`rasterio` and `requests` are now core dependencies, because `acquire/dem.py` needs them. Install them with **mamba first** (`mamba install -c conda-forge rasterio requests`), then `pip install -e ".[dev]"`. Done in that order, pip sees them as already installed and won't pull a pip rasterio wheel that carries its own GDAL. The `qgis` line was removed from `environment.yml`: QGIS is installed separately, as noted above.

## Where it lives

`environment.yml`, `pyproject.toml`, `.vscode/settings.json` (which tells VS Code to use conda).

## Revisit when

- Adding PDAL/rasterio in Phase 2. If the solve fails on Windows, pin versions or split into `env-core.yml` / `env-ml.yml`.
- Moving to the Scala/Spark scale-up. That has its own build (sbt), separate from this environment.
