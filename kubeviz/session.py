"""Session start-up and batch mode (ports of the ``kubeviz`` main procedure and
``kubeviz_batchmode``).

:func:`start_session` builds a :class:`~kubeviz.state.State` from the same keywords
the IDL program accepted, loads the cube, initialises the lines and the error
method and optionally loads spaxel masks and a results file. The GUI and the
command-line batch mode both start from here.
"""
from __future__ import annotations

import os

import numpy as np

from . import utils
from .constants import (ERR_BOOTSTRAP, ERR_MC1, ERR_MC2, ERR_MC3, ERR_METHOD_NAMES, ERR_NOISE,
                        FIT_GAUSS, FIT_MOMENTS, MODE_SPAXEL)
from .core.masks import load_mask
from .core.montecarlo import (createmc1cubes, createmc2cubes, createmc3cubes, readbootstrapcubes,
                              setupmontecarlocubes)
from .fitting.fitall import fitall
from .fitting.flags import autoflag
from .fitting.linefit import linefit_init
from .io.cube import load_cube, splitpath
from .io.results import loadres, saveres
from .io.session import SESSION_EXT, load_session, save_session
from .state import State


def start_session(datafile: str | None, noisefile=None, ext=None, noise_ext=None, trim=None,
                  fluxfac=None, smooth=1, specsmooth=1, transpose=False, vacuum=False,
                  logarithmic=False, waveunit=None, redshift=None, do_mc_errors=None, bootstrap=None,
                  nmontecarlo=None, plot_mc_pdf=False, save_mc_pdf=False, use_mc_noise=False,
                  lineset=0, debug=False, mask_sn_thresh=None, mask_maxvelerr=None, mom_thresh=None,
                  instr=None, band=None, fix_ratios=False, fittype="gauss", spmask=None, outdir=None,
                  logfile=None, fitpars=None, momwindowfile=None, gaussinitfile=None,
                  scalenoiseerr=False, resfile=None, scroll=False, nproc=1) -> State:
    """Create and initialise a session (headless). Returns the populated State."""
    utils.setup_logging(level=10 if debug else 20, logfile=logfile)
    utils.log.info("")
    utils.log.info("                              *** KUBEVIZ ***")
    utils.log.info("")
    from . import __version__
    utils.log.info(f"version: {__version__} (Python port of K2.2)")
    utils.log.info("")

    if datafile is None:
        raise ValueError("A datacube file (or a saved session) is required")
    if not os.path.exists(datafile):
        raise FileNotFoundError(f"Unable to locate file: {datafile}")

    if datafile.lower().endswith((".sav", SESSION_EXT)):
        state = load_session(datafile)
        finish_loaded_session(state)
        if spmask:
            load_mask(state, spmask)
        return state

    state = State()
    state.debug = bool(debug)
    state.nproc = int(1 if nproc is None else nproc)
    state.cwdir = os.getcwd() + os.sep
    state.outdir = (outdir if outdir else state.cwdir)
    if not state.outdir.endswith(os.sep):
        state.outdir += os.sep
    os.makedirs(state.outdir, exist_ok=True)

    dataname, datadir = splitpath(datafile)
    state.filename, state.indir = dataname, datadir

    # ------------------------------------------------------------------ error method
    bootstrap_path = None
    if bootstrap is not None:
        do_mc_errors = 1
    if do_mc_errors is not None:
        do_mc_errors = int(do_mc_errors)
        if do_mc_errors == ERR_BOOTSTRAP:
            bname, bdir = splitpath(bootstrap or "")
            if not bname:
                utils.warn("invalid bootstrap filename. Use Noise Cube errors instead.")
                state.domontecarlo = ERR_NOISE
            else:
                state.domontecarlo = ERR_BOOTSTRAP
                bootstrap_path = bdir if bdir else datadir
                state.bootstrap_file = os.path.join(bootstrap_path, bname)
        elif do_mc_errors in (ERR_MC1, ERR_MC2, ERR_MC3):
            state.domontecarlo = do_mc_errors
            n = int(nmontecarlo) if nmontecarlo is not None else 100
            setattr(state, {ERR_MC1: "Nmc1", ERR_MC2: "Nmc2", ERR_MC3: "Nmc3"}[do_mc_errors], n)
        else:
            state.domontecarlo = ERR_NOISE

    res_path = None
    if resfile is not None:
        rname, rdir = splitpath(resfile)
        if not rname:
            utils.warn("Invalid results filename. ")
        else:
            res_path = os.path.join(rdir if rdir else datadir, rname)

    state.plotMonteCarlodistrib = bool(plot_mc_pdf)
    state.saveMonteCarlodistrib = bool(save_mc_pdf)
    state.useMonteCarlonoise = bool(use_mc_noise)
    if instr:
        state.instr = instr
    if band:
        state.band = band
    if fluxfac is not None:
        state.fluxfac = float(fluxfac)
    state.scroll = bool(scroll)
    state.selected_lineset = int(lineset or 0)
    state.smooth = int(smooth) if smooth else 1
    state.specsmooth = int(specsmooth) if specsmooth else 1
    state.transpose = bool(transpose)
    if vacuum:
        state.vacuum = True
    if redshift is not None:
        state.redshift = float(redshift)
    if mask_sn_thresh is not None:
        state.mask_sn_thresh = float(mask_sn_thresh)
    if mask_maxvelerr is not None:
        state.mask_maxvelerr = float(mask_maxvelerr)
    if mom_thresh is not None:
        state.mom_thresh = float(mom_thresh)
    state.fitfixratios = bool(fix_ratios)
    state.scaleNoiseerrors = bool(scalenoiseerr)
    if fittype:
        if fittype == "gauss":
            state.linefit_type = FIT_GAUSS
        elif fittype == "moments":
            state.linefit_type = FIT_MOMENTS
        else:
            utils.warn('fittype keyword value not recognized. Use "gauss" instead. ')
    if fitpars is not None:
        state.set_user_fitpars(fitpars)
    state.zcuts = 4
    state.zoomrange = 32

    # ------------------------------------------------------------------ data
    load_cube(state, datafile, noisefile=noisefile, ext=ext, noise_ext=noise_ext, trim=trim,
              logarithmic=logarithmic, waveunit=waveunit, momwindowfile=momwindowfile,
              gaussinitfile=gaussinitfile)
    linefit_init(state)

    # ------------------------------------------------------------------ Monte Carlo cubes
    if state.domontecarlo == ERR_BOOTSTRAP:
        readbootstrapcubes(state, splitpath(bootstrap)[0], bootstrap_path, momwindowfile=momwindowfile)
    elif state.domontecarlo == ERR_MC1:
        createmc1cubes(state)
    elif state.domontecarlo == ERR_MC2:
        createmc2cubes(state)
    elif state.domontecarlo == ERR_MC3:
        createmc3cubes(state)
    if state.domontecarlo > 0:
        setupmontecarlocubes(state)

    if spmask:
        load_mask(state, spmask)
    if res_path is not None:
        loadres(state, res_path)
    return state


def finish_loaded_session(state: State) -> None:
    """Fix paths of a loaded session and regenerate what a light session dropped."""
    state.cwdir = os.getcwd() + os.sep
    if not os.path.isdir(state.outdir or ""):
        state.outdir = state.cwdir
    if state.noise is None:
        state.noise = state.noisecube
    restore_montecarlo(state)


def restore_montecarlo(state: State) -> None:
    """Recreate the Monte Carlo / bootstrap realisations of a session saved without
    them. MC1/2/3 are deterministic (fixed seeds); bootstrap cubes are re-read from
    ``state.bootstrap_file``. Falls back to noise-cube errors when impossible."""
    if state.domontecarlo == ERR_NOISE or state.montecarlocubes is not None:
        return
    utils.info(f"Regenerating {ERR_METHOD_NAMES[state.domontecarlo]} realisations for the loaded session...")
    if state.domontecarlo == ERR_BOOTSTRAP:
        if state.bootstrap_file and os.path.exists(state.bootstrap_file):
            readbootstrapcubes(state, os.path.basename(state.bootstrap_file), os.path.dirname(state.bootstrap_file))
        else:
            utils.warn("Bootstrap file not available: switching to Noise Cube errors.")
            state.domontecarlo = ERR_NOISE
            state.noise = state.noisecube
            return
    elif state.domontecarlo == ERR_MC1:
        createmc1cubes(state)
    elif state.domontecarlo == ERR_MC2:
        createmc2cubes(state)
    elif state.domontecarlo == ERR_MC3:
        createmc3cubes(state)
    setupmontecarlocubes(state)


def make_demo_cube(outdir: str | None = None, **kw) -> str:
    """Write a synthetic demo cube and return its path (``kubeviz --demo``)."""
    from .synth import make_synthetic_cube
    outdir = outdir or os.getcwd()
    os.makedirs(outdir, exist_ok=True)
    fname = os.path.join(outdir, "kubeviz_demo_cube.fits")
    opts = dict(nx=40, ny=32, redshift=0.02, nan_corner=True)
    opts.update(kw)
    make_synthetic_cube(fname, **opts)
    utils.setup_logging()
    utils.info(f"Synthetic demo cube written to {fname} (redshift 0.02, Halpha + [NII], MUSE-like)")
    return fname


def batchmode(state: State, redshift=None, fit_all_lines: bool = False, lineset=None):
    """Fit all spaxels (or masks), autoflag, save results and session
    (port of ``kubeviz_batchmode``). Returns ``(results_file, session_file)``."""
    if redshift is None and state.redshift == 0:
        raise ValueError("Batch mode does not work if the redshift keyword is not set. ")
    if state.Nlines == 0:
        raise ValueError("No lines identified in the spectrum for the given redshift and lineset. ")
    if fit_all_lines:
        utils.info("Setting fit parameters for all lines in set.")
        for iline in range(state.Nlines):
            state.pndofit[iline] = 1
            state.pcdofit[iline] = 1
    else:
        state.pndofit[0] = 1
        state.pcdofit[0] = 1
    str_err = ERR_METHOD_NAMES[state.domontecarlo]
    str_fit = "Gauss" if state.linefit_type == FIT_GAUSS else "Moments"
    if state.linefit_mode == MODE_SPAXEL:
        utils.info(f"Fitting spaxels with method: {str_fit}")
        utils.info(f"Fitting spaxels with error method: {str_err}")
        str_resmode = "spax"
        state.fitallrange = np.array([0, state.Ncol - 1, 0, state.Nrow - 1])
        fitall(state, spaxel=True)
        state.col, state.row, state.wpix = state.Ncol // 2, state.Nrow // 2, state.Nwpix // 2
        autoflag(state)
    else:
        utils.info(f"Fitting loaded masks with method: {str_fit}")
        utils.info(f"Fitting loaded masks with error method: {str_err}")
        str_resmode = "mask"
        fitall(state, spaxel=False)
    tag = f"{str_fit[:3].lower()}_{str_resmode}"
    resfile = os.path.join(state.outdir, f"{state.basename()}_res_{tag}.fits")
    saveres(state, resfile)
    sessfile = save_session(state, os.path.join(state.outdir, f"{state.basename()}_{tag}{SESSION_EXT}"), light=True)
    utils.info("Quitting...")
    return resfile, sessfile
