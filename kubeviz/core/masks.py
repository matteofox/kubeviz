"""Spaxel masks (ports of ``kubeviz_clear_select``, ``kubeviz_newmask``,
``kubeviz_deletemask``, ``kubeviz_save_mask``, ``kubeviz_load_mask``,
``kubeviz_optimal_mask`` and the mask-select branches of ``kubeviz_spax_event``).

Masks are stored in ``state.spaxselect`` with shape ``(Nmask, Nrow, Ncol)``; a value
of 1 selects a spaxel, any other positive value acts as a weight for the weighted
average extraction.
"""
from __future__ import annotations

import numpy as np
from astropy.io import fits
from scipy.optimize import least_squares

from .. import utils
from ..constants import MODE_MASK, SPEC_SUM
from .extraction import image_for_imgmode, medianspec


def init_masks(state) -> None:
    if state.spaxselect is None or state.spaxselect.shape[1:] != (state.Nrow, state.Ncol):
        state.spaxselect = np.zeros((max(state.Nmask, 1), state.Nrow, state.Ncol), dtype=float)
        state.Nmask = state.spaxselect.shape[0]
        state.imask = 1


def clear_select(state) -> None:
    state.spaxselect[state.imask - 1] = 0
    state.medspec = np.zeros(state.Nwpix)
    state.nmedspec = np.zeros(state.Nwpix)


def newmask(state) -> None:
    state.spaxselect = np.concatenate(
        [state.spaxselect, np.zeros((1, state.Nrow, state.Ncol), dtype=float)], axis=0)
    state.Nmask += 1
    state.imask = state.Nmask


def deletemask(state) -> None:
    if state.Nmask == 1:
        utils.warn("Last mask: Do not delete (can clear spaxels instead).")
        return
    idx = state.imask - 1
    state.spaxselect = np.delete(state.spaxselect, idx, axis=0)
    state.Nmask -= 1
    # keep mask results aligned with the mask list (IDL left them misaligned)
    for (mode, _err), rs in state.results.items():
        if mode == MODE_MASK:
            rs.delete_mask(idx)
    if state.imask > 1:
        state.imask -= 1


def save_mask(state, fname: str) -> None:
    fits.writeto(fname, np.asarray(state.spaxselect, dtype=np.float32), overwrite=True)
    utils.info(f"Wrote selected spaxel positions to file: {fname}")


def load_mask(state, fname: str) -> None:
    """Load masks from FITS (OVERWRITES the current masks) and switch to mask mode."""
    utils.info(f"Reading spaxels mask file: {fname}")
    m = np.asarray(fits.getdata(fname), dtype=float)
    if m.ndim == 2:
        m = m[None]
    if m.shape[1:] != (state.Nrow, state.Ncol):
        raise ValueError(f"Mask shape {m.shape[1:]} does not match the cube {(state.Nrow, state.Ncol)}")
    state.spaxselect = m
    state.Nmask = m.shape[0]
    state.imask = 1
    state.linefit_mode = MODE_MASK
    state.specmode = SPEC_SUM
    state.cursormode = 2
    medianspec(state)


def select_spaxels(state, col: int, row: int, value: int = 1) -> None:
    """Select (``value=1``) or deselect (``value=0``) spaxels around ``(col,row)`` in
    the current mask according to ``state.maskmode`` / ``state.maskradius``.

    maskmode 0: single spaxel; 1: circle of radius ``maskradius``;
    2: square of side ``maskradius`` (centred; even sides extend one less to the top/right).
    """
    m = state.spaxselect[state.imask - 1]
    nrow, ncol = m.shape
    r = int(state.maskradius)
    if state.maskmode == 0:
        if 0 <= row < nrow and 0 <= col < ncol:
            m[row, col] = value
    elif state.maskmode == 1:
        yy, xx = np.ogrid[:nrow, :ncol]
        dist = np.sqrt((xx - col) ** 2 + (yy - row) ** 2)
        m[dist <= r] = value
    else:
        half = r // 2
        x0, x1 = max(col - half, 0), min(col + half, ncol - 1)
        y0, y1 = max(row - half, 0), min(row + half, nrow - 1)
        if r % 2 == 0:
            x1 = min(x1, col + half - 1)
            y1 = min(y1, row + half - 1)
        if x1 >= x0 and y1 >= y0:
            m[y0:y1 + 1, x0:x1 + 1] = value


def optimal_mask(state, sn_thresh: float = 5.0) -> bool:
    """Select spaxels with S/N > ``sn_thresh`` in the current image
    (port of ``kubeviz_optimal_mask``). Returns True if any spaxel was selected."""
    clear_select(state)
    image, noise = image_for_imgmode(state)
    sn = utils.getsnimage(image, noise)
    sel = sn > sn_thresh
    if not np.any(sel):
        return False
    state.spaxselect[state.imask - 1][sel] = 1
    state.cursormode = 2
    state.linefit_mode = MODE_MASK
    return True


def gauss2d_tilt(params, x, y):
    base, peak, sx, sy, cx, cy, pa = params
    xp = (x - cx) * np.cos(pa) - (y - cy) * np.sin(pa)
    yp = (x - cx) * np.sin(pa) + (y - cy) * np.cos(pa)
    return base + peak * np.exp(-0.5 * ((xp / sx) ** 2 + (yp / sy) ** 2))


def findcentroid(state):
    """Fit a tilted 2D Gaussian to the current image around the current spaxel
    (port of ``kubeviz_findcentroid``, mpfit2dpeak /TILT).

    Returns the parameter vector ``[baseline, peak, sigmaX, sigmaY, centX, centY, PA]``.
    """
    image, noise = image_for_imgmode(state)
    with np.errstate(all="ignore"):
        if state.imgmode == 0:
            weights = (1.0 / image) ** 2
        else:
            weights = 1.0 / noise ** 2
    weights = np.where(np.isfinite(weights), weights, 0.0)
    ok = np.isfinite(image)
    nz = image[ok & (image != 0)]
    baseline = utils.idl_median(nz) if nz.size else 0.0
    peak = float(image[state.row, state.col]) if ok[state.row, state.col] else float(np.nanmax(image))
    p0 = np.array([baseline, peak, 2.0, 2.0, state.col, state.row, 0.0])
    yy, xx = np.mgrid[:state.Nrow, :state.Ncol]
    xs, ys, vals, ws = xx[ok], yy[ok], image[ok], np.sqrt(weights[ok])

    def resid(p):
        return ws * (vals - gauss2d_tilt(p, xs, ys))

    res = least_squares(resid, p0, method="lm" if xs.size >= 7 else "trf")
    A = res.x
    utils.info("Best-fit Gauss parameters: baseline, amplitude, sigmaX, sigmaY, centX, centY:")
    utils.info(" ".join(f"{v:.4g}" for v in A[:6]))
    return A
