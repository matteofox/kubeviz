"""Loops over spaxels and masks (ports of ``kubeviz_linefit_fitall``,
``kubeviz_linefit_fitadj``, ``kubeviz_linefit_fitadjall``,
``kubeviz_linefit_Nokadjacentspax`` and ``kubeviz_linefit_guessfromadjacentspax``).

Long loops accept ``should_cancel()`` (polled after every fit) and ``on_progress``
callbacks so that a GUI can interrupt and display progress.
"""
from __future__ import annotations

import numpy as np

from .. import utils
from ..constants import FLAG_MANUAL, NOT_FIT
from .flags import adjflag
from .linefit import dofit, linefit_reset, linefit_resetuser, startminmaxconsistencycheck, testspec

_DX = np.array([-1, 0, 1, 1, 1, 0, -1, -1])
_DY = np.array([1, 1, 1, 0, -1, -1, -1, 0])


def _neighbours(state, x, y, flag):
    xn, yn = x + _DX, y + _DY
    ok = (xn >= 0) & (xn <= state.Ncol - 1) & (yn >= 0) & (yn <= state.Nrow - 1)
    xn, yn = xn[ok], yn[ok]
    good = flag[yn, xn] == 0
    return xn[good], yn[good]


def n_ok_adjacent(state, x, y, flag) -> int:
    xn, _ = _neighbours(state, x, y, flag)
    return int(xn.size)


def guess_from_adjacent(state, x, y, flag, rs) -> bool:
    """Set ``state.gnfit``/``gbfit`` to the inverse-variance weighted mean of the OK
    neighbours' parameters. Returns False when no OK neighbour exists."""
    xn, yn = _neighbours(state, x, y, flag)
    if xn.size == 0:
        return False
    npar = 3 + state.Nlines
    for par in range(1, npar):
        for res, err, gfit, lt in ((rs.n, rs.nerr, state.gnfit, "N"), (rs.b, rs.berr, state.gbfit, "B")):
            e0 = err[yn, xn, par, 0]
            if np.any(e0 == NOT_FIT):
                continue
            values = res[yn, xn, par]
            eok = 0.5 * (np.abs(e0) + np.abs(err[yn, xn, par, 1]))
            with np.errstate(all="ignore"):
                guess = np.nansum(values / eok ** 2) / np.nansum(1.0 / eok ** 2)
            if not np.isfinite(guess):
                guess = 0.0
            gfit[par] = startminmaxconsistencycheck(state, guess, par, lt, silent=True)
    return True


def fitall(state, spaxel: bool = True, should_cancel=None, on_progress=None, on_fit=None) -> bool:
    """Fit every spaxel in ``state.fitallrange`` (or every mask).

    Returns True when the loop completed, False when cancelled.
    """
    state.zoomspax = 1
    completed = True
    if spaxel:
        c0, c1, r0, r1 = [int(v) for v in state.fitallrange]
        total = (c1 - c0 + 1) * (r1 - r0 + 1)
        prog = utils.Progress(total, state.percent_step, label="processed", callback=on_progress)
        cur = (state.col, state.row)
        for col in range(c0, c1 + 1):
            for row in range(r0, r1 + 1):
                state.col, state.row = col, row
                prog.step()
                if testspec(state):
                    dofit(state)
                    if on_fit is not None:
                        on_fit(col, row)
                if should_cancel is not None and should_cancel():
                    completed = False
                    break
            if not completed:
                break
        prog.close()
        state.col, state.row = cur
    else:
        cur = state.imask
        for mask in range(1, state.Nmask + 1):
            state.imask = mask
            utils.info(f"Fitting mask {mask}")
            if testspec(state):
                dofit(state)
                if on_fit is not None:
                    on_fit(None, None)
            if on_progress is not None:
                on_progress(mask / state.Nmask, f"Fitting mask {mask}")
            if should_cancel is not None and should_cancel():
                completed = False
                break
        state.imask = cur
    return completed


def fitadj(state) -> bool:
    """Fit the current spaxel using the neighbours' solution as initial guess."""
    rs = state.get_results()
    nuser, buser = state.gnfit.copy(), state.gbfit.copy()
    ok = guess_from_adjacent(state, state.col, state.row, rs.c[..., 0], rs)
    if ok:
        dofit(state)
        state.gnfit, state.gbfit = nuser, buser
        adjflag(state, state.col, state.row)
    return ok


def fitadjall(state, should_cancel=None, on_progress=None, on_fit=None) -> int:
    """Iteratively refit BAD spaxels adjacent to OK ones, most-OK-neighbours first.
    Returns the number of spaxels whose fit was improved."""
    rs = state.get_results()
    colstart, rowstart = state.col, state.row
    bad = state.badpixelimg

    flag = np.minimum(rs.n[..., 0], rs.b[..., 0])
    for res, lt in ((rs.n, "N"), (rs.b, "B")):
        ys, xs = np.nonzero(res[..., 0] == FLAG_MANUAL)
        for y, x in zip(ys, xs):
            state.col, state.row = int(x), int(y)
            linefit_reset(state, linetypes=(lt,))

    badfit = (~bad) & (flag > 0)
    nbad_init = int(np.sum(badfit))
    hope = np.ones((state.Nrow, state.Ncol), dtype=bool)
    nokneigh = np.zeros((state.Nrow, state.Ncol), dtype=int)
    yy, xx = np.mgrid[:state.Nrow, :state.Ncol]
    c0, c1, r0, r1 = [int(v) for v in state.fitallrange]
    inrange = (xx >= c0) & (xx <= c1) & (yy >= r0) & (yy <= r1)
    hope[~inrange] = False
    badfit[~inrange] = False
    niter = state.Nloopadj if state.Nloopadj > 0 else 100000
    cancelled = False

    while np.sum(badfit & hope) > 0 and niter >= 0:
        niter -= 1
        ys, xs = np.nonzero(badfit & inrange)
        for y, x in zip(ys, xs):
            nokneigh[y, x] = n_ok_adjacent(state, x, y, flag)
        hope[nokneigh == 0] = False
        nfitnow, ii = 0, 0
        for ii in range(8, 0, -1):
            now = hope & (nokneigh == ii) & badfit & inrange
            nfitnow = int(np.sum(now))
            if nfitnow > 0:
                break
        if nfitnow == 0:
            break
        msg = f"[PROGRES] Running step with {ii} neighbours, for {nfitnow} pixels."
        utils.log.info(msg)
        if on_progress is not None:
            on_progress(1.0 - np.sum(badfit & hope) / max(nbad_init, 1), msg)
        ys, xs = np.nonzero(now)
        for y, x in zip(ys, xs):
            state.col, state.row = int(x), int(y)
            if guess_from_adjacent(state, x, y, flag, rs):
                dofit(state)
                if on_fit is not None:
                    on_fit(int(x), int(y))
                newflag = adjflag(state, int(x), int(y))
                flag = np.minimum(rs.n[..., 0], rs.b[..., 0])
                if newflag == 0:
                    hope[max(y - 1, 0):min(y + 2, state.Nrow), max(x - 1, 0):min(x + 2, state.Ncol)] = True
                else:
                    hope[y, x] = False
            if should_cancel is not None and should_cancel():
                cancelled = True
                break
        if cancelled:
            break
        badfit = (~bad) & (flag > 0) & inrange

    linefit_resetuser(state, all_pars=True, startonly=True)
    state.col, state.row = colstart, rowstart
    nbad_final = int(np.sum((~bad) & (flag > 0) & inrange))
    utils.info(f"The fit was improved for {nbad_init - nbad_final} spaxels.")
    return nbad_init - nbad_final
