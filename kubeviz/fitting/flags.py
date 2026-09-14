"""Quality flags (ports of ``kubeviz_linefit_autoflag``, ``kubeviz_linefit_adjflag``
and ``kubeviz_linefit_flagsigclip``).

Flag = sum of binary flags (see :mod:`kubeviz.constants`):
16 manual/bad spaxel, 8 kinematics outside 4 sigma of the map (fitadj only),
4 invalid velocity/dispersion errors, 2 errors above threshold, 1 S/N below threshold.
"""
from __future__ import annotations

import numpy as np

from ..constants import (FIT_GAUSS, FLAG_INVALID_VEL, FLAG_LOW_SN, FLAG_MANUAL, FLAG_SIGCLIP,
                         FLAG_VELERR, MODE_SPAXEL, NO_ERRORS, NOT_FIT)
from .. import utils


def _valid_err(err, i, fitconstr):
    """Velocity/dispersion error plane ``i`` valid (both +/- entries)."""
    e0, e1 = err[..., i, 0], err[..., i, 1]
    if fitconstr:
        return (e0 != NOT_FIT) & (e1 != NOT_FIT)
    return (e0 != 0) & (e0 != NOT_FIT) & (e1 != 0) & (e1 != NOT_FIT)


def _err_below(err, i, maxvelerr):
    return (np.abs(err[..., i, 0]) < maxvelerr) & (np.abs(err[..., i, 1]) < maxvelerr)


def _sn_from_minus_sigma(res, err, nlines):
    """S/N of each line flux relative to the -1 sigma error only (DJW 12-10-2016)."""
    with np.errstate(all="ignore"):
        sn = res[..., 3:3 + nlines] / np.abs(err[..., 3:3 + nlines, 1])
    sn[~np.isfinite(sn)] = 0.0
    return sn


def flagsigclip(state, array, mode):
    """[lower, upper] acceptance range for velocity (mode 1) or width (mode 2) maps."""
    if mode == 1:
        ok = (array != 0) & (array != state.gnlims[1, 0]) & (array != state.gnlims[1, 1])
    else:
        ok = (array > 0) & (array != state.gnlims[2, 0]) & (array != state.gnlims[2, 1])
    arr = array[ok]
    if arr.size > 1:
        p16, p50, p84 = utils.percentile(arr, [16, 50, 84])
        return [p50 - 4 * (p50 - p16), p50 + 4 * (p84 - p50)]
    return [-1e9, 1e9]


def _gauss_spaxel_flags(state, rs, sel=None, sigclip=False):
    """Flag arrays (n, b, combined) for all spaxels or a single (row, col) selection."""
    nl = state.Nlines
    n, b = rs.n, rs.b
    nerr, berr = rs.nerr, rs.berr
    sn_n = _sn_from_minus_sigma(n, nerr, nl)
    sn_b = _sn_from_minus_sigma(b, berr, nl)
    if state.selected_lineset == 1 and not sigclip:
        maxsnn, maxsnb = sn_n.max(axis=-1), sn_b.max(axis=-1)
    else:
        mi = state.mainline_index()
        maxsnn, maxsnb = sn_n[..., mi], sn_b[..., mi]
    mve = state.mask_maxvelerr
    maxerrokn = _err_below(nerr, 1, mve) & _err_below(nerr, 2, mve)
    maxerrokb = _err_below(berr, 1, mve) & _err_below(berr, 2, mve)
    velvalokn = _valid_err(nerr, 1, state.fitconstr) & _valid_err(nerr, 2, state.fitconstr)
    velvalokb = _valid_err(berr, 1, state.fitconstr) & _valid_err(berr, 2, state.fitconstr)
    thr = state.mask_sn_thresh
    setmaskn = (FLAG_LOW_SN * (maxsnn < thr) + FLAG_VELERR * (~maxerrokn) + FLAG_INVALID_VEL * (~velvalokn)).astype(int)
    setmaskb = (FLAG_LOW_SN * (maxsnb < thr) + FLAG_VELERR * (~maxerrokb) + FLAG_INVALID_VEL * (~velvalokb)).astype(int)
    if sigclip:
        nvel = flagsigclip(state, n[..., 1], 1)
        nsig = flagsigclip(state, n[..., 2], 2)
        bvel = flagsigclip(state, b[..., 1], 1)
        bsig = flagsigclip(state, b[..., 2], 2)
        clipokn = (n[..., 1] > nvel[0]) & (n[..., 1] < nvel[1]) & (n[..., 1] != 0) & (n[..., 2] < nsig[1]) & (n[..., 2] > 0)
        clipokb = (b[..., 1] > bvel[0]) & (b[..., 1] < bvel[1]) & (b[..., 1] != 0) & (b[..., 2] < bsig[1]) & (b[..., 2] > 0)
        setmaskn += FLAG_SIGCLIP * (~clipokn)
        setmaskb += FLAG_SIGCLIP * (~clipokb)
    else:
        bad = state.badpixelimg.astype(int)
        setmaskn += FLAG_MANUAL * bad
        setmaskb += FLAG_MANUAL * bad
    return setmaskn, setmaskb


def autoflag(state, doflag=None):
    """Flag all spaxels/masks of the current results (port of ``kubeviz_linefit_autoflag``).
    ``doflag`` (spaxel mode) is a boolean image selecting the spaxels to (re)flag.
    Returns the new narrow-line flag array."""
    rs = state.get_results()
    nl = state.Nlines
    npercs = rs.npercs
    if state.linefit_type == FIT_GAUSS:
        if state.linefit_mode == MODE_SPAXEL:
            if doflag is None:
                doflag = np.ones((state.Nrow, state.Ncol), dtype=int)
            doflag = np.asarray(doflag).astype(int)
            setmaskn, setmaskb = _gauss_spaxel_flags(state, rs)
            setmask = np.minimum(setmaskn, setmaskb)
            newn = doflag * setmaskn + (1 - doflag) * rs.n[..., 0]
            newb = doflag * setmaskb + (1 - doflag) * rs.b[..., 0]
            newc = doflag * setmask + (1 - doflag) * rs.c[..., 0]
            rs.n[..., 0], rs.b[..., 0], rs.c[..., 0] = newn, newb, newc
            for i in range(npercs):
                rs.nerr[..., 0, i] = newn
                rs.berr[..., 0, i] = newb
                rs.cerr[..., 0, i] = newc
            return newn
        # masks
        n, b, nerr, berr = rs.n, rs.b, rs.nerr, rs.berr
        with np.errstate(all="ignore"):
            n1 = 0.5 * (nerr[:, 3:3 + nl, 0] - nerr[:, 3:3 + nl, 1])
            b1 = 0.5 * (berr[:, 3:3 + nl, 0] - berr[:, 3:3 + nl, 1])
            sn_n = n[:, 3:3 + nl] / n1
            sn_b = b[:, 3:3 + nl] / b1
        sn_n[(nerr[:, 3:3 + nl, 0] == 0) & (nerr[:, 3:3 + nl, 1] == 0)] = 0.0
        sn_b[(berr[:, 3:3 + nl, 0] == 0) & (berr[:, 3:3 + nl, 1] == 0)] = 0.0
        sn_n[~np.isfinite(sn_n)] = 0.0
        sn_b[~np.isfinite(sn_b)] = 0.0
        maxsn = np.maximum(sn_n.max(axis=1), sn_b.max(axis=1)) if nl > 0 else np.zeros(n.shape[0])
        mve = state.mask_maxvelerr
        velerrok = ((_valid_err(nerr, 1, False) | _valid_err(berr, 1, False)) &
                    (_valid_err(nerr, 2, False) | _valid_err(berr, 2, False)) &
                    (_err_below(nerr, 1, mve) | _err_below(berr, 1, mve)) &
                    (_err_below(nerr, 2, mve) | _err_below(berr, 2, mve)))
        setok = (maxsn > state.mask_sn_thresh) & velerrok
        flag = (~setok).astype(float)
        rs.n[:, 0], rs.b[:, 0], rs.c[:, 0] = flag, flag, flag
        for i in range(npercs):
            rs.nerr[:, 0, i] = flag
            rs.berr[:, 0, i] = flag
            rs.cerr[:, 0, i] = flag
        return flag

    # moments
    m, merr = rs.m, rs.merr
    if not (np.max(merr[..., 0, 0]) > NO_ERRORS):
        return m[..., 5] if nl > 0 else None
    mve = state.mask_maxvelerr
    for il in range(nl):
        with np.errstate(all="ignore"):
            onesig = 0.5 * (merr[..., 6 * il, 0] - merr[..., 6 * il, 1])
            sn = m[..., 6 * il] / onesig
        sn[(merr[..., 6 * il, 0] == 0) & (merr[..., 6 * il, 1] == 0)] = 0.0
        sn[~np.isfinite(sn)] = 0.0
        velerrok = (_valid_err(merr, 6 * il + 1, False) & _valid_err(merr, 6 * il + 2, False) &
                    _err_below(merr, 6 * il + 1, mve) & _err_below(merr, 6 * il + 2, mve))
        setok = (sn > state.mask_sn_thresh) & velerrok
        if state.linefit_mode == MODE_SPAXEL:
            setok &= ~state.badpixelimg
        m[..., 6 * il + 5] = (~setok).astype(float)
        for i in range(npercs):
            merr[..., 6 * il + 5, i] = (~setok).astype(float)
    return m[..., 5]


def adjflag(state, xfit: int, yfit: int) -> int:
    """Flag a single spaxel after a fit-adjacent attempt, including the 4-sigma clip
    against the velocity/dispersion maps (port of ``kubeviz_linefit_adjflag``).
    Returns the combined (min of narrow/broad) flag."""
    rs = state.get_results()
    if state.linefit_type == FIT_GAUSS:
        setmaskn, setmaskb = _gauss_spaxel_flags(state, rs, sigclip=True)
        fn, fb = int(setmaskn[yfit, xfit]), int(setmaskb[yfit, xfit])
        fc = min(fn, fb)
        rs.n[yfit, xfit, 0], rs.b[yfit, xfit, 0], rs.c[yfit, xfit, 0] = fn, fb, fc
        rs.nerr[yfit, xfit, 0, :] = fn
        rs.berr[yfit, xfit, 0, :] = fb
        rs.cerr[yfit, xfit, 0, :] = fc
        return fc
    autoflag(state)
    mi = state.mainline_index()
    return int(rs.m[yfit, xfit, 6 * mi + 5]) if state.Nlines > 0 else 0
