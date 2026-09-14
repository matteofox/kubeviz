"""Headless display logic: which image is currently shown and how the fit-result
maps are built (ports of the image-selection part of ``kubeviz_plotspax``,
``kubeviz_linefit_image_update`` and the value lookup of ``kubeviz_plotinfo``).
"""
from __future__ import annotations

import numpy as np

from .. import utils
from ..constants import (NOT_FIT, CUBE_BADPIX, CUBE_DATA, CUBE_LINEFIT, CUBE_LINEFIT_ERR, CUBE_LINEFIT_SN,
                         CUBE_NOISE, CUBE_SN, FIT_GAUSS, IMG_MED1_MINUS_MED2, IMG_MED2_MINUS_MED1,
                         IMG_MEDSUB1, IMG_MEDSUB2, IMG_SLICE, MODE_SPAXEL)
from ..results import moment_reshape
from .extraction import image_for_imgmode

IMGMODE_NAMES = ["Slice", "Sum1", "Median1", "Weighted Avg1", "Weighted Med1", "MedSub1",
                 "Sum2", "Median2", "Weighted Avg2", "Weighted Med2", "MedSub2",
                 "Med2-Med1", "Med1-Med2"]
CUBESEL_NAMES = {CUBE_DATA: "DATA", CUBE_NOISE: "NOISE", CUBE_BADPIX: "BAD", CUBE_SN: "S/N",
                 CUBE_LINEFIT: "FIT", CUBE_LINEFIT_ERR: "FIT ERR", CUBE_LINEFIT_SN: "FIT S/N"}


def bad_image_for_imgmode(state) -> np.ndarray:
    """Bad-pixel image matching the current image mode (True = bad)."""
    im = state.imgmode
    if im == IMG_SLICE:
        return state.badpixelmask[state.wpix]
    if im in (IMG_MEDSUB1,):
        return state.bimg1 & state.badpixelmask[state.wpix]
    if im in (IMG_MEDSUB2,):
        return state.bimg2 & state.badpixelmask[state.wpix]
    if im in (IMG_MED2_MINUS_MED1, IMG_MED1_MINUS_MED2):
        return state.bimg1 & state.bimg2
    if 1 <= im <= 5:
        return state.bimg1
    return state.bimg2


def current_image(state):
    """(image, badimage) to display for ``state.cubesel`` and ``state.imgmode``.

    Bad pixels are NOT yet set to NaN; the caller decides how to render them.
    """
    if state.cubesel in (CUBE_LINEFIT, CUBE_LINEFIT_ERR, CUBE_LINEFIT_SN):
        bad = state.badpixelimg
        if state.cubesel == CUBE_LINEFIT:
            img = state.lineresimg if state.flagmode else state.unm_lineresimg
        elif state.cubesel == CUBE_LINEFIT_ERR:
            img = state.lineerrresimg if state.flagmode else state.unm_lineerrresimg
        else:
            if state.flagmode:
                img = utils.getsnimage(state.lineresimg, state.lineerrresimg, absolute=True)
            else:
                img = utils.getsnimage(state.unm_lineresimg, state.unm_lineerrresimg, absolute=True)
        return np.asarray(img, dtype=float), np.asarray(bad, dtype=bool)

    bad = bad_image_for_imgmode(state)
    data, noise = image_for_imgmode(state)
    if state.cubesel == CUBE_DATA:
        img = data
    elif state.cubesel == CUBE_NOISE:
        img = noise
    elif state.cubesel == CUBE_BADPIX:
        img = bad.astype(float)
    else:
        img = utils.getsnimage(data, noise)
    return np.asarray(img, dtype=float), np.asarray(bad, dtype=bool)


def current_value(state) -> float:
    """Value shown in the info panel for the current spaxel (``kubeviz_plotinfo``)."""
    img, bad = current_image(state)
    r, c = state.row, state.col
    if state.cubesel == CUBE_BADPIX:
        return float(bad[r, c])
    if state.cubesel <= CUBE_SN and bad[r, c]:
        return np.nan
    if state.cubesel == CUBE_LINEFIT_SN:
        err = state.lineerrresimg if state.flagmode else state.unm_lineerrresimg
        val = state.lineresimg if state.flagmode else state.unm_lineresimg
        return float(abs(val[r, c]) / err[r, c]) if err[r, c] > 0 else 0.0
    return float(img[r, c])


def linefit_image_update(state, newpar: str) -> bool:
    """Build the (masked and unmasked) result images for parameter code ``newpar``
    (``'FLAG'``, ``'CHISQ'``, or ``'N3'``/``'B1'``/``'C4'`` = component + parameter
    index). Returns False when there is no valid result for that plane, in which
    case ``state.par_imagebutton`` is left unchanged (IDL reverted the button).
    """
    oldpar = state.par_imagebutton
    if not newpar:
        return False
    rs = state.get_results(MODE_SPAXEL)
    nl = state.Nlines
    if state.linefit_type == FIT_GAUSS:
        n, b, c = rs.n, rs.b, rs.c
        nerr, berr, cerr = rs.nerr, rs.berr, rs.cerr
    else:
        mi = state.mainline_index()
        c, cerr = rs.c, rs.cerr
        n = b = moment_reshape(rs.m, mi, nl)
        nerr = berr = moment_reshape(rs.merr, mi, nl, error=True)
    shape = (state.Nrow, state.Ncol)

    # spaxels never fitted are blank (NaN) in every map; IDL showed them as 0
    unfit = nerr[..., 1, 0] == NOT_FIT
    if newpar == "CHISQ":
        ok = (c[..., 0] == 0) & ~unfit
        img = np.full(shape, np.nan)
        img[ok] = rs.chisq[ok]
        state.lineresimg = img
        state.lineerrresimg = img.copy()
        unm = np.where(unfit, np.nan, rs.chisq)
        state.unm_lineresimg = unm
        state.unm_lineerrresimg = unm.copy()
        state.par_imagebutton = newpar
        return True
    if newpar == "FLAG":
        flag = np.where(unfit, np.nan, n[..., 0].astype(float))
        state.lineresimg = flag
        state.unm_lineresimg = flag.copy()
        ferr = np.where(unfit, np.nan, nerr[..., 0, 0].astype(float))
        state.lineerrresimg = ferr
        state.unm_lineerrresimg = ferr.copy()
        state.par_imagebutton = newpar
        return True

    lt, par = newpar[0], int(newpar[1:])
    if lt == "N":
        ok = n[..., 0] == 0
        plane = n[..., par]
        errplane = 0.5 * (nerr[..., par, 0] - nerr[..., par, 1])
        refplane = nerr[..., par, 0]
    elif lt == "B":
        ok = b[..., 0] == 0
        plane = b[..., par]
        errplane = 0.5 * (berr[..., par, 0] - berr[..., par, 1])
        refplane = berr[..., par, 0]
    else:
        ok = c[..., 0] == 0
        cpar = par - 2
        plane = c[..., cpar]
        errplane = 0.5 * (cerr[..., cpar, 0] - cerr[..., cpar, 1])
        refplane = cerr[..., cpar, 0]
    unfit = refplane == NOT_FIT
    ok = ok & ~unfit
    plane = np.where(unfit, np.nan, plane)
    errplane = np.where(unfit, np.nan, errplane)

    if not np.any(ok):
        state.unm_lineresimg = plane.copy()
        state.unm_lineerrresimg = errplane.copy()
        utils.warn("No spaxels with OK flag.")
        if state.flagmode:
            state.par_imagebutton = oldpar
            return False
        state.par_imagebutton = newpar
        return True
    maxerr = np.nanmax(np.abs(errplane)) if np.any(np.isfinite(errplane)) else 0.0
    if (maxerr == 0.0 or maxerr == 999.0) and np.min(np.abs(refplane)) > 998:
        utils.warn("No valid results in this plane as yet.")
        if state.flagmode:
            state.par_imagebutton = oldpar
            return False
    img = np.full(shape, np.nan)
    eimg = np.full(shape, np.nan)
    img[ok] = plane[ok]
    eimg[ok] = errplane[ok]
    state.lineresimg, state.lineerrresimg = img, eimg
    unm, unme = plane.copy(), errplane.copy()
    unm[state.badpixelimg] = np.nan
    unme[state.badpixelimg] = np.nan
    state.unm_lineresimg, state.unm_lineerrresimg = unm, unme
    state.par_imagebutton = newpar
    return True


def linefit_image_reset(state) -> None:
    state.par_imagebutton = ""
    state.cubesel = CUBE_DATA
    utils.warn("No valid results in this plane for this error method as yet.")
    utils.info("Reverting to data cube")


def current_flag(state) -> int:
    """Flag of the current spaxel/mask for the component selected in the image
    button (``'b'`` key / FLAG button logic)."""
    rs = state.get_results()
    idx = rs.index(col=state.col, row=state.row) if state.linefit_mode == MODE_SPAXEL else rs.index(imask=state.imask)
    if state.linefit_type != FIT_GAUSS:
        mi = state.mainline_index()
        return int(rs.m[idx + (6 * mi + 5,)]) if state.Nlines > 0 else 0
    lt = state.par_imagebutton[:1]
    arr = {"B": rs.b, "C": rs.c}.get(lt, rs.n)
    return int(arr[idx + (0,)])


def toggle_flag(state, all_components: bool = False) -> int:
    """Toggle the manual flag (0 <-> 16) of the current spaxel/mask. With
    ``all_components`` (FLAG button) narrow, broad and continuum are all set; the
    ``'b'`` key sets only the component of the selected result image."""
    from ..constants import FLAG_MANUAL
    rs = state.get_results()
    idx = rs.index(col=state.col, row=state.row) if state.linefit_mode == MODE_SPAXEL else rs.index(imask=state.imask)
    flag = current_flag(state)
    f1 = FLAG_MANUAL if flag == 0 else 0
    if state.linefit_type != FIT_GAUSS:
        mi = state.mainline_index()
        if state.Nlines > 0:
            rs.m[idx + (6 * mi + 5,)] = f1
            rs.merr[idx + (6 * mi + 5,)] = f1
        return f1
    lt = state.par_imagebutton[:1]
    targets = [("N", rs.n, rs.nerr), ("B", rs.b, rs.berr), ("C", rs.c, rs.cerr)]
    if not all_components:
        targets = [t for t in targets if t[0] == (lt if lt in "NBC" and lt else "N")]
    for _, arr, err in targets:
        arr[idx + (0,)] = f1
        err[idx + (0,)] = f1
    return f1
