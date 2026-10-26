"""Terrain derivatives on the 1 m DEM.

Phase 2 grows this module (slope, hillshade, LRM proper). For now it holds the two
pieces the alignment check needs. Decision notes: docs/decisions/0008-label-alignment-check.md
"""
from __future__ import annotations

import numpy as np


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
