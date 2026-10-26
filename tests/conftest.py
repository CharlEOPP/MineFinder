from pathlib import Path

import pytest
import yaml

REPO = Path(__file__).resolve().parents[1]

# Small synthetic AOI (EPSG:26913 metres) so tests don't depend on the real config.
#   extent  : 0..10_000 x 0..10_000 (offset to plausible UTM 13N values)
#   holdout : SW 2 km block, 100 m buffer
X0, Y0 = 450_000, 4_420_000

AOI_DICT = {
    "name": "test_aoi",
    "crs": "EPSG:26913",
    "extent": {"minx": X0, "miny": Y0, "maxx": X0 + 10_000, "maxy": Y0 + 10_000},
    "holdout": {"minx": X0, "miny": Y0, "maxx": X0 + 2_000, "maxy": Y0 + 2_000},
    "holdout_buffer_m": 100,
}


def write_yaml(tmp_path: Path, data: dict, name: str = "aoi.yaml") -> Path:
    p = tmp_path / name
    p.write_text(yaml.safe_dump(data))
    return p


@pytest.fixture
def aoi_path(tmp_path):
    return write_yaml(tmp_path, AOI_DICT)


@pytest.fixture
def aoi(aoi_path):
    from minefinder.aoi import load_aoi

    return load_aoi(aoi_path)
