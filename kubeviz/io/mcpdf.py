"""Monte Carlo result distributions: FITS tables of the realisations
(``kubeviz_linefit_saveMonteCarlodistrib``) and cumulative-distribution plots
(``kubeviz_linefit_plotMonteCarlodistrib``, PNG instead of EPS)."""
from __future__ import annotations

import os

import numpy as np
from astropy.io import fits

from .. import utils
from ..constants import ERR_BOOTSTRAP, ERR_MC1, ERR_MC2, ERR_MC3, FIT_GAUSS, MODE_SPAXEL

_MONTEDIR = {ERR_BOOTSTRAP: "boot", ERR_MC1: "MC1", ERR_MC2: "MC2", ERR_MC3: "MC3"}
_PLOTDIR = {ERR_BOOTSTRAP: "bootstrap", ERR_MC1: "montecarlo1", ERR_MC2: "montecarlo2", ERR_MC3: "montecarlo3"}
_warned_mpl = False


def save_montecarlo_distrib(state, ctx, bootpars, cbootpars) -> str:
    typestr = "gau" if state.linefit_type == FIT_GAUSS else "mom"
    if state.linefit_mode == MODE_SPAXEL:
        locstr = "spax"
        extstr = f"SPX_{state.col + 1}_{state.row + 1}"
        nexten = state.Ncol * state.Nrow
    else:
        locstr = "mask"
        extstr = f"MASK_{state.imask}"
        nexten = state.Nmask
    fname = os.path.join(state.outdir, f"{state.basename()}_{typestr}_{locstr}_PDF.fits")

    bootpars = np.asarray(bootpars, dtype=float)
    cbootpars = np.asarray(cbootpars, dtype=float) if cbootpars is not None else np.zeros((bootpars.shape[0], 0))
    if state.linefit_type == FIT_GAUSS:
        pdf = np.hstack([bootpars[:, 1:], cbootpars])      # drop the continuum parameter
    else:
        pdf = np.hstack([bootpars, cbootpars])

    hdr = fits.Header()
    hdr["EXTNAME"] = extstr
    hdr["DATACUBE"] = (state.filename, "Name of the original datacube")
    if state.linefit_mode == MODE_SPAXEL:
        hdr["COL"] = (state.col + 1, "Spaxel column number")
        hdr["ROW"] = (state.row + 1, "Spaxel row number")
    else:
        hdr["MASK"] = (state.imask, "Number of the current mask")
    hdr["FITTYPE"] = (int(state.linefit_type), "Type of fit: 0: gauss 1: moments")
    hdr["RESTYPE"] = (int(state.linefit_mode), "Type of result: 0: spaxels 1: masks")
    hdr["ERMETHOD"] = (int(state.domontecarlo), "Err method: 0:Noise, 1:Bootstrap 2,3,4:MC1,2,3")
    hdr["MCNOISE"] = (int(state.useMonteCarlonoise), "MonteCarlo noise: 0: Off 1: On")
    hdr["SPATSMTH"] = (int(state.smooth), "Spatial smoothing applied")
    hdr["SPECSMTH"] = (int(state.specsmooth), "Spectral smoothing applied")
    hdr["LINESET"] = (int(state.selected_lineset), "Lineset used")
    names = state.linenames
    ind = 1
    if state.linefit_type == FIT_GAUSS:
        for comp, fp in (("narrow", ctx.nfitpars), ("broad", ctx.bfitpars)):
            if fp.size and fp[0] == 1:
                hdr[f"COL{ind}"] = f"{comp} line dv"
                ind += 1
            if fp.size and fp[1] == 1:
                hdr[f"COL{ind}"] = f"{comp} line sig(v)"
                ind += 1
            for i in range(2, fp.size):
                if fp[i] == 1:
                    hdr[f"COL{ind}"] = f"{comp} line {names[i - 2]} flux"
                    ind += 1
        for comp, fp in (("narrow", ctx.nfitpars), ("broad", ctx.bfitpars)):
            for i in range(2, fp.size):
                if fp[i] == 1:
                    hdr[f"COL{ind}"] = f"{comp} line {names[i - 2]} continuum flux"
                    ind += 1
    else:
        for i in range(ctx.mfitpars.size):
            if ctx.mfitpars[i] == 1:
                for k, lab in enumerate(["0th (flux)", "1st (dv)", "2nd (sigv)", "3rd (norm Skewness)", "4th (norm Kurtosis)"]):
                    hdr[f"COL{ind + k}"] = f"moment {names[i]} {lab}"
                ind += 5
        for i in range(ctx.mfitpars.size):
            if ctx.mfitpars[i] == 1:
                hdr[f"COL{ind}"] = f"moment {names[i]} continuum"
                ind += 1

    if not os.path.exists(fname):
        prihdr = fits.Header()
        prihdr["NEXT"] = nexten
        prihdr["NMC"] = state.Nmontecarlo
        fits.HDUList([fits.PrimaryHDU(header=prihdr)]).writeto(fname)
    with fits.open(fname, mode="update") as hdul:
        existing = [i for i, h in enumerate(hdul) if h.name == extstr]
        new = fits.ImageHDU(data=pdf, header=hdr, name=extstr)
        if existing:
            hdul[existing[0]] = new
        else:
            hdul.append(new)
        hdul.flush()
    return fname


def plot_montecarlo_distrib(state, dist, percs_boot, val, restype, par) -> str | None:
    """Cumulative distribution plot of one parameter over the realisations."""
    global _warned_mpl
    try:
        import matplotlib
        matplotlib.use("Agg", force=False)
        import matplotlib.pyplot as plt
    except Exception:
        if not _warned_mpl:
            utils.warn("matplotlib not available: Monte Carlo distribution plots disabled.")
            _warned_mpl = True
        return None
    if state.linefit_mode == MODE_SPAXEL:
        locstr = f"spax_{state.col:2d}_{state.row:2d}".replace(" ", "")
    else:
        locstr = f"mask_{state.imask - 1:3d}".replace(" ", "")
    names = state.linenames
    if par == 0:
        strpar = "flag"
    elif restype == "cont":
        strpar = f"flux@{names[par - 1]}"
    elif par == 1:
        strpar = "dv"
    elif par == 2:
        strpar = "sig(v)"
    else:
        strpar = f"{names[par - 3]}_flux"
    montedir = _PLOTDIR.get(state.domontecarlo, "")
    outdir = os.path.join(state.cwdir or ".", montedir)
    os.makedirs(outdir, exist_ok=True)
    fname = os.path.join(outdir, f"{locstr}_{restype}_{strpar}_distr.png")

    dist = np.asarray(dist, dtype=float)
    sx = np.sort(dist)
    y = np.arange(sx.size) / float(sx.size)
    fig, ax = plt.subplots(figsize=(6, 4.5))
    ax.plot(sx, y, "k-", lw=2)
    ax.axvline(val, color="k", lw=2.5)
    styles = [("k", ":"), ("k", ":"), ("b", ":"), ("b", ":"), ("r", "-"), ("g", "--"), ("g", "--"), ("0.5", "--"), ("0.5", "--")]
    for i, pv in enumerate(percs_boot):
        col, ls = styles[i] if i < len(styles) else ("b", "--")
        ax.axvline(pv, color=col, ls=ls, lw=1)
        ax.axhline(0.01 * state.montecarlo_percs[i], color=col, ls=ls, lw=1)
    ax.set_xlabel(f"{restype} {strpar}")
    ax.set_ylabel(f"f({restype} {strpar}>x)")
    ax.set_title(f"Error method: {montedir}")
    ax.set_ylim(0, 1)
    fig.tight_layout()
    fig.savefig(fname, dpi=110)
    plt.close(fig)
    return fname
