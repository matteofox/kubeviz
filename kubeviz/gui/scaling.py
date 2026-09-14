"""Intensity scaling for the spaxel viewer (the ``zcuts`` cases of ``kubeviz_plotspax``)
and colour tables (``kubeviz_setcolour`` / ``kubeviz_modify_colour``).

Pure numpy, no Qt, so it can be tested headless.
"""
from __future__ import annotations

import numpy as np

from .. import utils

ZCUT_USER_LIN = 1
ZCUT_MINMAX = 2
ZCUT_ZSCALE = 3
ZCUT_HISTEQ = 4
ZCUT_USER_SQRT = 5
ZCUT_USER_LOG = 6
ZCUT_995 = 11
ZCUT_990 = 12
ZCUT_970 = 13
ZCUT_950 = 14

ZCUT_NAMES = {ZCUT_HISTEQ: "HistEq", ZCUT_ZSCALE: "Zscale", ZCUT_MINMAX: "MinMax", ZCUT_995: "99.5%",
              ZCUT_990: "99%", ZCUT_970: "97%", ZCUT_950: "95%", ZCUT_USER_LIN: "User linear",
              ZCUT_USER_SQRT: "User sqrt", ZCUT_USER_LOG: "User log10"}


def zscale_range(image, contrast=0.25):
    """IRAF zscale limits (astropy implementation when available)."""
    img = np.asarray(image, dtype=float)
    finite = img[np.isfinite(img)]
    if finite.size == 0:
        return 0.0, 1.0
    try:
        from astropy.visualization import ZScaleInterval
        lo, hi = ZScaleInterval(contrast=contrast).get_limits(finite)
    except Exception:  # pragma: no cover
        lo, hi = np.percentile(finite, [1, 99])
    if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
        lo, hi = float(finite.min()), float(finite.max())
    return float(lo), float(hi)


def hist_equal(image, nbins=1024):
    """Histogram equalisation to [0,1] (IDL ``hist_equal``). Returns (scaled, (min, max))."""
    img = np.asarray(image, dtype=float)
    finite = np.isfinite(img)
    out = np.full(img.shape, np.nan)
    vals = img[finite]
    if vals.size == 0:
        return out, (0.0, 0.0)
    vmin, vmax = float(vals.min()), float(vals.max())
    if vmin == vmax:
        out[finite] = 0.0
        return out, (vmin, vmax)
    hist, edges = np.histogram(vals, bins=nbins, range=(vmin, vmax))
    cdf = np.cumsum(hist).astype(float)
    cdf /= cdf[-1]
    idx = np.clip(((vals - vmin) / (vmax - vmin) * nbins).astype(int), 0, nbins - 1)
    out[finite] = cdf[idx]
    return out, (vmin, vmax)


def scale_image(image, zcuts: int, zmin_user: float, zmax_user: float):
    """Map ``image`` to [0,1] according to ``zcuts``. Returns ``(scaled, (lo, hi))``;
    NaN stays NaN."""
    img = np.asarray(image, dtype=float)
    finite = np.isfinite(img)
    if not np.any(finite):
        return np.full(img.shape, np.nan), (0.0, 0.0)

    def linear(lo, hi):
        if hi <= lo:
            hi = lo + 1e-30
        with np.errstate(invalid="ignore"):
            return np.clip((img - lo) / (hi - lo), 0, 1), (float(lo), float(hi))

    if zcuts == ZCUT_USER_LIN:
        return linear(zmin_user, zmax_user)
    if zcuts == ZCUT_MINMAX:
        return linear(np.nanmin(img), np.nanmax(img))
    if zcuts == ZCUT_ZSCALE:
        return linear(*zscale_range(img))
    if zcuts == ZCUT_HISTEQ:
        return hist_equal(img)
    if zcuts in (ZCUT_USER_SQRT, ZCUT_USER_LOG):
        lo, hi = zmin_user, zmax_user
        pos = finite & (img > 0)
        out = np.zeros(img.shape)
        out[~finite] = np.nan
        if np.any(pos):
            if lo <= 0:
                lo = float(np.nanmin(img[pos]))
            f = np.sqrt if zcuts == ZCUT_USER_SQRT else np.log10
            flo, fhi = f(lo), f(max(hi, lo * (1 + 1e-9)))
            with np.errstate(all="ignore"):
                out[pos] = np.clip((f(img[pos]) - flo) / (fhi - flo), 0, 1)
        return out, (float(lo), float(hi))
    pct = {ZCUT_995: 99.5, ZCUT_990: 99.0, ZCUT_970: 97.0, ZCUT_950: 95.0}.get(zcuts)
    if pct is not None:
        lo, hi = utils.dataclip(img, pct)
        return linear(lo, hi)
    return linear(np.nanmin(img), np.nanmax(img))


# ----------------------------------------------------------------------------- colour tables
COLOUR_TABLES = [("Grey", 0), ("BlueWhite", 1), ("Heat", 3), ("STD GAMMA-II", 5), ("Rainbow", 13)]
_MPL_NAMES = {0: "gray", 1: "Blues_r", 3: "gist_heat", 5: "gnuplot2", 13: "jet"}


def _bluewhite():
    x = np.linspace(0, 1, 256)
    r = np.clip(2 * x - 1, 0, 1)
    g = np.clip(2 * x - 1, 0, 1) * 0.6 + np.clip(4 * x - 2, 0, 1) * 0.4
    b = np.clip(2 * x, 0, 1)
    return np.stack([r, g, b], axis=1)


def colour_table(ctab: int) -> np.ndarray:
    """256x3 float RGB in [0,1] for an IDL colour table number."""
    try:
        import matplotlib
        cmap = matplotlib.colormaps[_MPL_NAMES.get(ctab, "gray")]
        return np.asarray(cmap(np.linspace(0, 1, 256))[:, :3], dtype=float)
    except Exception:
        x = np.linspace(0, 1, 256)
        if ctab == 1:
            return _bluewhite()
        if ctab == 3:
            return np.stack([np.clip(3 * x, 0, 1), np.clip(3 * x - 1, 0, 1), np.clip(3 * x - 2, 0, 1)], axis=1)
        if ctab == 13:
            r = np.clip(1.5 - abs(4 * x - 3), 0, 1)
            g = np.clip(1.5 - abs(4 * x - 2), 0, 1)
            b = np.clip(1.5 - abs(4 * x - 1), 0, 1)
            return np.stack([r, g, b], axis=1)
        return np.stack([x, x, x], axis=1)


def contrast_lut_indices(x0: float, y0: float) -> np.ndarray:
    """Index remapping used while right-dragging in the spaxel viewer (from ATV via
    IDL kubeviz). ``x0, y0`` are the fractional cursor positions in the window."""
    brightness, contrast = 0.5, 0.5
    x = brightness * 255 * x0
    y = max(contrast * 255 * y0, 2)
    high, low = x + y, x - y
    diff = max(high - low, 1)
    slope = 255.0 / diff
    intercept = -slope * low
    p = (np.arange(256) * slope + intercept).astype(int)
    return np.clip(p, 0, 255)


def render_rgb(scaled, lut_rgb, invert=False, indices=None, bad_color=(1.0, 1.0, 1.0)) -> np.ndarray:
    """Turn a [0,1] scaled image into uint8 RGB using ``lut_rgb`` (256x3, floats);
    NaN pixels get ``bad_color``."""
    lut = lut_rgb[::-1] if invert else lut_rgb
    if indices is not None:
        lut = lut[indices]
    idx = np.zeros(scaled.shape, dtype=int)
    finite = np.isfinite(scaled)
    idx[finite] = np.clip((scaled[finite] * 255).astype(int), 0, 255)
    rgb = lut[idx]
    rgb[~finite] = bad_color
    return (rgb * 255).astype(np.uint8)
