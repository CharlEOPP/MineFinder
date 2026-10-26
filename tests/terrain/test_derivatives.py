import numpy as np
import pytest

from minefinder.terrain.derivatives import local_relief, pit_depth


def plane(n=101, sx=0.3, sy=-0.2):
    yy, xx = np.mgrid[0:n, 0:n].astype("float64")
    return 2500 + sx * xx + sy * yy


def test_local_relief_of_a_tilted_plane_is_zero():
    # A box mean of a plane equals the centre value, so relief removes regional slope.
    lrm = local_relief(plane(), window=21)
    assert np.allclose(lrm[10:-10, 10:-10], 0, atol=1e-9)


def test_local_relief_shows_a_pit_as_negative():
    z = plane()
    z[50, 50] -= 3.0
    lrm = local_relief(z, window=21)
    assert lrm[50, 50] < -2.5
    assert abs(lrm[10, 10]) < 1e-9


def test_local_relief_ignores_nan_neighbours():
    # Flat ground: NaN pixels must not leak into neighbours' means.
    # (On a slope, a window with missing pixels is off-centre, so a small bias
    # near nodata edges is expected; alignment skips windows with gaps for that reason.)
    z = plane(sx=0, sy=0)
    z[45:55, 60:70] = np.nan
    lrm = local_relief(z, window=21)
    assert np.isnan(lrm[50, 65])                     # nodata stays nodata
    assert abs(lrm[50, 50]) < 1e-9                    # neighbours not polluted by NaN


def test_local_relief_rejects_even_window():
    with pytest.raises(ValueError, match="odd"):
        local_relief(plane(), window=20)


def test_pit_depth_finds_deepest_point_within_radius():
    lrm = np.zeros((21, 21))
    lrm[10, 13] = -2.0                                # 3 px east of centre
    d = pit_depth(lrm, radius=4)
    assert d[10, 10] == pytest.approx(2.0)            # within 4 px -> seen
    assert d[10, 4] == pytest.approx(0.0)             # 9 px away -> not seen
