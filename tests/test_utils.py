import numpy as np

from kubeviz import utils
from kubeviz.core.extraction import weighted_median_axis0
from kubeviz.core.smooth import smooth_cube
from kubeviz.linesdb import airtovac, vactoair


def test_idl_median_even_odd():
    a = np.array([1.0, 2.0, 3.0, 4.0])
    assert utils.idl_median(a) == 3.0            # upper middle (IDL default)
    assert utils.idl_median(a, even=True) == 2.5
    assert utils.idl_median(np.array([5.0, 1.0, 3.0])) == 3.0
    b = np.array([[1.0, np.nan], [4.0, 2.0], [3.0, 8.0]])
    np.testing.assert_allclose(utils.idl_median(b, axis=0), [3.0, 8.0])


def test_percentile_matches_idl_definition():
    x = np.random.default_rng(0).random(37)
    perc = np.array([16.0, 50.0, 84.0])
    s = np.sort(x)
    pos = 0.01 * perc * (x.size - 1)
    p1 = np.floor(pos).astype(int)
    p2 = np.minimum(p1 + 1, x.size - 1)
    ref = (p2 - pos) * s[p1] + (pos - p1) * s[p2]
    np.testing.assert_allclose(utils.percentile(x, perc), ref)


def test_weighted_median_matches_loop():
    rng = np.random.default_rng(1)
    for _ in range(20):
        a = rng.random(15)
        w = rng.random(15)
        tot = w.sum()
        ind = np.argsort(a)
        s = 0.0
        n = 0
        for i in range(a.size):
            s += w[ind][i]
            n += 1
            if s > tot / 2.0:
                break
        assert utils.weighted_median(a, w) == a[ind][n - 1]
        np.testing.assert_allclose(weighted_median_axis0(a[:, None], w[:, None])[0], a[ind][n - 1])


def test_air_vacuum_roundtrip():
    w = np.array([4861.333, 6562.819, 9068.6])
    v = airtovac(w)
    assert np.all(v > w)
    np.testing.assert_allclose(vactoair(v), w, rtol=1e-8)
    assert abs(airtovac(6562.819) - 6564.632) < 0.005


def _brute_smooth(data, noise, s, ss):
    nw, nrow, ncol = data.shape
    out = np.full(data.shape, np.nan)
    outn = np.full(data.shape, np.nan)
    lo_s, hi_s = (s // 2, s // 2 - 1) if s % 2 == 0 else ((s - 1) // 2, (s - 1) // 2)
    lo_z, hi_z = (ss // 2, ss // 2 - 1) if ss % 2 == 0 else ((ss - 1) // 2, (ss - 1) // 2)
    for k in range(nw):
        for j in range(nrow):
            for i in range(ncol):
                box = data[max(k - lo_z, 0):k + hi_z + 1, max(j - lo_s, 0):j + hi_s + 1, max(i - lo_s, 0):i + hi_s + 1]
                nbox = noise[max(k - lo_z, 0):k + hi_z + 1, max(j - lo_s, 0):j + hi_s + 1, max(i - lo_s, 0):i + hi_s + 1]
                out[k, j, i] = np.nanmedian(box)
                outn[k, j, i] = np.sqrt(np.nansum(nbox ** 2)) / np.sum(np.isfinite(nbox))
    if s % 2 == 0:
        out[:, 0, :] = out[:, :, 0] = np.nan
        outn[:, 0, :] = outn[:, :, 0] = np.nan
    if ss % 2 == 0:
        out[0] = outn[0] = np.nan
    return out, outn


def test_smooth_cube_matches_brute_force():
    rng = np.random.default_rng(2)
    data = rng.random((7, 6, 5))
    noise = rng.random((7, 6, 5)) + 0.1
    data[3, 2, 2] = np.nan
    for s, ss in ((3, 1), (2, 1), (3, 3), (4, 2), (1, 3)):
        got, gotn = smooth_cube(data, noise, s, ss)
        ref, refn = _brute_smooth(data, noise, s, ss)
        np.testing.assert_allclose(got, ref, equal_nan=True)
        np.testing.assert_allclose(gotn, refn, equal_nan=True)
    same, samen = smooth_cube(data, noise, 1, 1)
    np.testing.assert_array_equal(same, data)
