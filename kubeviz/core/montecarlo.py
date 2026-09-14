"""Bootstrap and Monte Carlo realisations of the data cube (ports of
``kubeviz_readbootstrapcubes``, ``kubeviz_setupmontecarlocubes``,
``kubeviz_createmc1cubes``, ``kubeviz_createmc2cubes``, ``kubeviz_createmc3cubes``).

Realisations are stored as one 4D array ``(N, Nwpix, Nrow, Ncol)`` per method.
Random numbers use ``numpy.random.default_rng`` with the IDL seed values, so runs are
reproducible in Python (but not bit-identical to IDL).
"""
from __future__ import annotations

import os

import numpy as np
from astropy.io import fits

from .. import utils
from ..constants import ERR_BOOTSTRAP, ERR_MC1, ERR_MC2, ERR_MC3, ERR_NOISE, MONTECARLO_PERCS
from .extraction import montecarlonoise
from .smooth import smooth_state_montecarlo


def readbootstrapcubes(state, bootstrap_fname: str, directory: str = "", momwindowfile: str | None = None) -> None:
    path = os.path.join(directory, bootstrap_fname) if directory else bootstrap_fname
    if not os.path.exists(path):
        utils.warn(f"Unable to locate file: {path}")
        state.domontecarlo = ERR_NOISE
        return
    with fits.open(path, memmap=True) as hdul:
        next_ = int(str(hdul[0].header.get("NEXT", 0)).strip() or 0)
        if next_ == 0:
            utils.warn(f"The bootstrap file: {path} appears to be corrupted. Return.")
            state.domontecarlo = ERR_NOISE
            return
        state.Nbootstrap = next_
        utils.info(f"Loading bootstrap data cube: {path}")
        z1, z2 = state.Startwpix, state.Startwpix + state.Nwpix
        y1, y2 = state.Startrow, state.Startrow + state.Nrow
        x1, x2 = state.Startcol, state.Startcol + state.Ncol
        cubes = np.empty((next_, state.Nwpix, state.Nrow, state.Ncol), dtype=np.float32)
        for ext in range(1, next_ + 1):
            cubes[ext - 1] = hdul[ext].section[z1:z2, y1:y2, x1:x2] * state.fluxfac
    state.inbootstrapcubes = cubes

    if momwindowfile is not None:
        maps = []
        with fits.open(momwindowfile) as hd:
            for ext in range(1, next_ + 1):
                if ext >= len(hd):
                    maps.append(None)
                    continue
                m = np.asarray(hd[ext].data, dtype=float)
                if m.ndim != 3 or m.shape[0] != 2:
                    maps.append(None)
                    continue
                if m.shape[1:] != (state.Nrow, state.Ncol):
                    out = np.full((2, state.Nrow, state.Ncol), -999.0)
                    r, c = min(m.shape[1], state.Nrow), min(m.shape[2], state.Ncol)
                    out[:, :r, :c] = m[:, :r, :c]
                    m = out
                maps.append(m)
        state.mom_windowmap_bootstrap = maps


def _finish_setup(state, attr_cubes: str, attr_noise: str, n: int, label: str) -> None:
    utils.info(f"Setting up {label} data cubes...")
    cubes = getattr(state, attr_cubes)
    bad = state.badpixelmask
    for i in range(cubes.shape[0]):
        c = cubes[i]
        c[bad] = 0.0
        c[~np.isfinite(c)] = 0.0
    setattr(state, attr_noise, montecarlonoise(state, "cube"))
    if state.useMonteCarlonoise:
        state.noise = getattr(state, attr_noise)
    state.Nmontecarlo = n
    # results containers for this error method are created lazily by state.get_results


def setupmontecarlocubes(state) -> None:
    """Mask bad pixels in the realisations, compute the MC noise cube and the
    percentile table (port of ``kubeviz_setupmontecarlocubes``)."""
    state.montecarlo_percs = np.array(MONTECARLO_PERCS)
    if state.domontecarlo == ERR_BOOTSTRAP:
        smooth_state_montecarlo(state)
        _finish_setup(state, "bootstrapcubes", "bootstrapnoise", state.Nbootstrap, "bootstrap")
    elif state.domontecarlo == ERR_MC1:
        smooth_state_montecarlo(state)
        _finish_setup(state, "mc1cubes", "mc1noise", state.Nmc1, "MonteCarlo 1")
    elif state.domontecarlo == ERR_MC2:
        _finish_setup(state, "mc2cubes", "mc2noise", state.Nmc2, "MonteCarlo 2")
    elif state.domontecarlo == ERR_MC3:
        _finish_setup(state, "mc3cubes", "mc3noise", state.Nmc3, "MonteCarlo 3")


def createmc1cubes(state) -> None:
    """Unsmoothed data + N(0,1) * unsmoothed noise, smoothed afterwards in setup."""
    utils.info("Creating MonteCarlo 1 cubes...")
    n = state.Nmc1
    out = np.empty((n,) + state.indatacube.shape, dtype=np.float32)
    for i in range(n):
        rng = np.random.default_rng(353 + i * 3)
        out[i] = state.indatacube + rng.standard_normal(state.indatacube.shape) * state.innoisecube
    state.inmc1cubes = out


def createmc2cubes(state) -> None:
    """Smoothed data + N(0,1) * smoothed noise."""
    utils.info("Creating MonteCarlo 2 cubes...")
    n = state.Nmc2
    out = np.empty((n,) + state.datacube.shape, dtype=np.float32)
    for i in range(n):
        rng = np.random.default_rng(353 + i * 3)
        out[i] = state.datacube + rng.standard_normal(state.datacube.shape) * state.noisecube
    state.mc2cubes = out


def createmc3cubes(state) -> None:
    """Smoothed data + N(0,1) * residuals of a first single-component fit, with the
    residual cube randomly shifted along wavelength."""
    from ..fitting.linefit import get_residual_spec

    utils.info("Creating MonteCarlo 3 cubes...")
    saved = (state.domontecarlo, state.pndofit.copy(), state.pbdofit.copy(), state.pcdofit.copy(),
             state.maxwoffb, state.maxwoffr, state.continuumfit_minoff, state.continuumfit_maxoff,
             state.continuumfit_minperc, state.continuumfit_maxperc, state.continuumfit_order)
    state.domontecarlo = ERR_NOISE
    for iline in range(state.Nlines):
        state.pndofit[iline] = 1
        state.pbdofit[iline] = 0
        state.pcdofit[iline] = 1
    state.set_continuumfit_defaults()
    state.set_linefitrange_defaults()

    residual = np.zeros(state.datacube.shape, dtype=float)
    prog = utils.Progress(state.Ncol * state.Nrow, state.percent_step, label="residual spectra")
    for col in range(state.Ncol):
        for row in range(state.Nrow):
            residual[:, row, col] = get_residual_spec(state, col, row)
            prog.step()
    prog.close()

    n = state.Nmc3
    out = np.empty((n,) + state.datacube.shape, dtype=np.float32)
    for i in range(n):
        rng1 = np.random.default_rng(424 + i * 3)
        rng2 = np.random.default_rng(385 + i * 3)
        lshift = int(rng1.standard_normal())
        out[i] = state.datacube + rng2.standard_normal(state.datacube.shape) * np.roll(residual, lshift, axis=0)
    state.mc3cubes = out

    (state.domontecarlo, state.pndofit, state.pbdofit, state.pcdofit,
     state.maxwoffb, state.maxwoffr, state.continuumfit_minoff, state.continuumfit_maxoff,
     state.continuumfit_minperc, state.continuumfit_maxperc, state.continuumfit_order) = saved
