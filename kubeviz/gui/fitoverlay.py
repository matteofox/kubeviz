"""Curves overplotted in the spectral zoom window (the fit-drawing part of
``kubeviz_plotspeczoom``). Pure numpy so it can be tested without Qt.

Returns a list of ``Overlay`` records: kind is one of ``"narrow"``, ``"broad"``,
``"moment"``, ``"total"``, ``"contline"`` (vertical line at a continuum-window edge).
"""
from __future__ import annotations

import warnings
from dataclasses import dataclass

import numpy as np

from ..constants import CKMS, CUBE_BADPIX, CUBE_NOISE, FIT_GAUSS, MODE_SPAXEL, NOT_FIT, SIGTOFWHM
from ..fitting.linefit import chooselines


@dataclass
class Overlay:
    kind: str
    x: np.ndarray
    y: np.ndarray


def _gauss(state, xx, lam0, dv, sigv, norm):
    if dv == NOT_FIT:
        dv = 0.0
    if sigv == NOT_FIT:
        sigv = 0.0
    pos = lam0 * (1.0 + dv / CKMS)
    width = lam0 * (sigv / CKMS)
    R = float(state.getinstrres(pos))
    instr = pos / (SIGTOFWHM * R) if R > 0 else 0.0
    tot2 = width ** 2 + instr ** 2
    if tot2 <= 0:
        return np.zeros_like(xx)
    return (norm / np.sqrt(2.0 * np.pi * tot2)) * np.exp(-((xx - pos) ** 2 / tot2) / 2.0)


def fit_overlays(state, xx) -> list[Overlay]:
    """Overlays for the wavelength array ``xx`` of the zoom window."""
    out: list[Overlay] = []
    nl = state.Nlines
    if nl == 0 or state.cubesel in (CUBE_NOISE, CUBE_BADPIX):
        return out
    rs = state.get_results()
    idx = rs.index(col=state.col, row=state.row) if state.linefit_mode == MODE_SPAXEL else rs.index(imask=state.imask)
    n, b, c, m = rs.n[idx], rs.b[idx], rs.c[idx], rs.m[idx]
    nfit, bfit, cfit = n[3:3 + nl], b[3:3 + nl], c[1:1 + nl]
    mfit = m[0:6 * nl:6]
    nshow = np.flatnonzero((nfit > 0) & (state.nshow[:nl] == 1))
    bshow = np.flatnonzero((bfit > 0) & (state.bshow[:nl] == 1))
    mshow = np.flatnonzero(state.nshow[:nl] == 1)
    gauss = state.linefit_type == FIT_GAUSS
    if gauss and nshow.size + bshow.size == 0:
        return out
    if not gauss and mshow.size == 0:
        return out

    totfit = np.zeros_like(xx)
    cont = np.full((xx.size, state.lineset_max + 1), np.nan)
    for ls in range(1, state.lineset_max + 1):
        lines_c = chooselines(state, state.linesets, state.lines, lineset=ls)
        if lines_c.size == 0:
            continue
        centline = np.median(state.lines[lines_c])
        inset = (xx >= centline - state.maxwoffb) & (xx <= centline + state.maxwoffr)
        if np.any(inset):
            vals = cfit[lines_c]
            okv = vals[(vals != NOT_FIT) & (vals != 0)]
            cont[inset, ls] = okv.mean() if okv.size else 0.0
        for lo, hi in ((centline - state.continuumfit_maxoff, centline - state.continuumfit_minoff),
                       (centline + state.continuumfit_minoff, centline + state.continuumfit_maxoff)):
            sel = np.flatnonzero((xx >= lo) & (xx <= hi))
            if sel.size:
                for edge in (xx[sel[0]], xx[sel[-1]]):
                    out.append(Overlay("contline", np.array([edge, edge]), np.array([-1e10, 1e10])))
    with np.errstate(all="ignore"), warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        contmed = np.nanmedian(cont, axis=1)
    contmed[~np.isfinite(contmed)] = 0.0
    totfit += contmed

    if gauss:
        for kind, show, res, fluxes in (("narrow", nshow, n, nfit), ("broad", bshow, b, bfit)):
            for line in show:
                ls = int(state.linesets[line])
                g = _gauss(state, xx, state.lines[line], res[1], res[2], fluxes[line])
                totfit += g
                if np.sum(g) > 0:
                    y = g + np.nan_to_num(cont[:, ls]) if state.cshow[line] == 1 else g
                    out.append(Overlay(kind, xx, y))
    else:
        for line in mshow:
            ls = int(state.linesets[line])
            g = _gauss(state, xx, state.lines[line], m[6 * line + 1], m[6 * line + 2], mfit[line])
            totfit += g
            if np.sum(g) > 0:
                y = g + np.nan_to_num(cont[:, ls]) if state.cshow[line] == 1 else g
                out.append(Overlay("moment", xx, y))
    ok = np.abs(totfit) > 1e-30
    if np.any(ok):
        out.append(Overlay("total", xx[ok], totfit[ok]))
    return out
