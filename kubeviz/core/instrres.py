"""Instrumental spectral resolution (ports of ``kubeviz_instrres_init`` and
``kubeviz_linefit_skylines``).

Three sources, selected by ``state.instrres_mode``:

* :data:`INSTRRES_EXTPOLY` - polynomial ``R(lambda)`` from the data header (KMOS3D
  products), from the KMOS arc templates shipped in ``kubeviz/data/templates/kmos``
  (closest in time to DATE-OBS, per IFU), or hardcoded curves (MUSE manual, VIMOS,
  SAMI, WiFeS);
* :data:`INSTRRES_TEMPLATE` - Gaussian sigma of a template sky-line profile (SINFONI;
  templates not shipped, falls back to the skyline fit);
* :data:`INSTRRES_VARPOLY` - 4th order polynomial fitted to sky lines identified in the
  noise (variance) cube, either at the wavelengths of the shipped OH template or found
  automatically in the variance spectrum.
"""
from __future__ import annotations

import glob
import os
import re

import numpy as np
from astropy.io import fits

from .. import utils
from ..constants import INSTRRES_EXTPOLY, INSTRRES_TEMPLATE, INSTRRES_VARPOLY, SIGTOFWHM
from ..fitting.mpfit import mpfitpeak

DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data", "templates")


def _hget(hdr, *keys, default=None):
    if hdr is None:
        return default
    for k in keys:
        for cand in (k, "HIERARCH " + k, k.replace(" ", "_"), k.replace("_", " ")):
            if cand in hdr:
                return hdr[cand]
    return default


def _read_kmos_template(fname: str):
    """Return dict ifu -> [p0..p4] from a KMOSpoly_arc_*.txt file (microns)."""
    out = {}
    with open(fname) as fh:
        for line in fh:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            parts = s.split()
            try:
                ifu = int(parts[0])
                coeff = [float(v) for v in parts[1:6]]
            except ValueError:
                continue
            out[ifu] = coeff
    return out


def instrres_init(state) -> None:
    """Fill ``state.instrres_*`` for the current instrument (port of ``kubeviz_instrres_init``)."""
    state.instrres_mode = INSTRRES_VARPOLY
    instr = state.instr.lower()
    prihead = state.inprihead if state.inprihead is not None else state.indatahead

    if instr == "kmos3d":
        resorder = _hget(prihead, "ESO K3D RES ORDER", "K3D RES ORDER", default=0)
        try:
            resorder = int(resorder)
        except (TypeError, ValueError):
            resorder = 0
        if resorder > 0:
            if np.sum(state.instrres_extpoly) == 0:
                utils.info("Instrumental resolution coefficients extracted from the primary header.")
            for i in range(resorder + 1):
                c = float(_hget(prihead, f"ESO K3D RES COEFF{i}", f"K3D RES COEFF{i}", default=0.0))
                state.instrres_extpoly[i] = c * (1.0e-4) ** i    # microns -> Angstrom
            state.instrres_mode = INSTRRES_EXTPOLY
        else:
            utils.warn("Polynomial resolution coefficients must be present in the primary header of the official KMOS3D cubes")
        return

    if instr == "kmos":
        resorder = _hget(prihead, "K3D RESORDER", "RESORDER", default=0)
        try:
            resorder = int(resorder)
        except (TypeError, ValueError):
            resorder = 0
        if resorder > 0:
            if np.sum(state.instrres_extpoly) == 0:
                utils.info("Instrumental resolution coefficients extracted from the primary header.")
            for i in range(resorder + 1):
                c = float(_hget(prihead, f"K3D RESP{i}", f"RESP{i}", default=0.0))
                state.instrres_extpoly[i] = c * (1.0e-4) ** i
            state.instrres_mode = INSTRRES_EXTPOLY
            return
        files = sorted(glob.glob(os.path.join(DATA_DIR, "kmos", f"KMOSpoly_arc_{state.band.lower()}_*.txt")))
        if not files:
            utils.warn(f"No KMOS resolution template for band {state.band!r}. Fitting skylines instead.")
            linefit_skylines(state)
            return
        scitime = str(_hget(prihead, "DATE-OBS", default="2014-01-01T00:00:00"))
        best, mindif = None, np.inf
        for f in files:
            m = re.search(r"_(\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2})", os.path.basename(f))
            if not m:
                continue
            try:
                dif = abs(utils.timediff_seconds(scitime, m.group(1)))
            except (ValueError, IndexError):
                dif = np.inf
            if dif < mindif:
                best, mindif = f, dif
        table = _read_kmos_template(best) if best else {}
        if state.ifu in table:
            if np.sum(state.instrres_extpoly) == 0:
                utils.info("Instrumental resolution coefficients extracted from the kubeviz archive")
            p = table[state.ifu]
            for i in range(5):
                state.instrres_extpoly[i] = p[i] * (1.0e-4) ** i
            state.instrres_mode = INSTRRES_EXTPOLY
        else:
            utils.warn(f"IFU {state.ifu} not found in {os.path.basename(best) if best else 'templates'}. Fitting skylines instead.")
            linefit_skylines(state)
        return

    if instr == "sinfoni":
        band = state.band.upper()
        tpl = None
        if band == "J" and state.pixscale == 250:
            tpl = "OH_GaussProf_j250.fits"
        elif band == "H" and state.pixscale == 250:
            tpl = "OH_GaussProf_h250.fits"
        elif band == "K" and state.pixscale in (100, 250):
            tpl = f"OH_GaussProf_k{state.pixscale}.fits"
        else:
            utils.warn(f"Band {state.band} / pixel scale {state.pixscale} not yet supported")
        path = os.path.join(DATA_DIR, "sinfoni", tpl) if tpl else None
        if path and os.path.exists(path):
            data, hdr = fits.getdata(path, header=True)
            crval, cdelt, crpix, naxis = hdr.get("CRVAL1", 0), hdr.get("CDELT1", 0), hdr.get("CRPIX1", 0), hdr.get("NAXIS1", 0)
            if crval == 0 and cdelt == 0 and crpix == 0:
                crval, cdelt, crpix, naxis = hdr.get("CRVAL2", 0), hdr.get("CDELT2", 0), hdr.get("CRPIX2", 0), hdr.get("NAXIS2", 0)
            tl = (np.arange(naxis) + 1 - crpix) * cdelt + crval
            _, params, _ = mpfitpeak(tl, np.asarray(data, dtype=float).ravel(), nterms=3, positive=True)
            state.instrres_tplsig = float(params[2]) * 1.0e4
            state.instrres_mode = INSTRRES_TEMPLATE
        else:
            utils.debug("Instrumental resolution templates missing")
            utils.info("Fitting skylines in the noisecube.")
            linefit_skylines(state)
        return

    if instr == "muse":
        if np.sum(state.instrres_extpoly) == 0:
            utils.info("Instrumental resolution obtained from the MUSE manual (Figure 11)")
        state.instrres_extpoly[:5] = [499.446, -0.201608, 0.000130290, -7.87479e-09, 0.0]
        state.instrres_mode = INSTRRES_EXTPOLY
        return

    if instr in ("vimos", "sami"):
        utils.info("Instrumental resolution is not properly implemented for this instrument.")
        utils.info("Please input the coefficients manually.")
        state.instrres_extpoly[:5] = [2650.0 if instr == "vimos" else 4580.0, 0, 0, 0, 0]
        state.instrres_mode = INSTRRES_EXTPOLY
        return

    if instr == "wifes":
        r0 = {"B": 3000.0, "R": 7000.0}.get(state.band.upper(), 3000.0)
        state.instrres_extpoly[:5] = [r0, 0, 0, 0, 0]
        state.instrres_mode = INSTRRES_EXTPOLY
        return

    linefit_skylines(state)


# ----------------------------------------------------------------------------- skylines
def _find_oh_lines(state, ohspec, lamboh):
    """Wavelengths of clean OH lines in the template overlapping the data."""
    sel = (lamboh > state.wave[0]) & (lamboh < state.wave[-1])
    ohspec = ohspec[sel]
    lamboh = lamboh[sel]
    nbins = 5
    binslim = int(len(lamboh) / nbins)
    zpos2 = []
    for b in range(nbins):
        zspec = ohspec[b * binslim:(b + 1) * binslim]
        fin = zspec[np.isfinite(zspec)]
        if fin.size == 0:
            continue
        maxval = fin.max()
        zline = np.where(zspec < maxval / 50.0, 0.0, zspec)
        zdiff = zline[:-1] - zline[1:]
        zpos = [i for i in range(10, len(zdiff) - 10) if zdiff[i] > 0 and zdiff[i - 1] < 0]
        for i in zpos:
            left = zspec[:i]
            right = zspec[i + 1:]
            leftlim = np.flatnonzero(left[::-1] < zspec[i] / 50.0)
            rightlim = np.flatnonzero(right < zspec[i] / 50.0)
            if leftlim.size and rightlim.size and abs(leftlim[0] - rightlim[0]) < 2:
                zpos2.append(i + b * binslim)
    return lamboh[np.array(zpos2, dtype=int)] if zpos2 else np.zeros(0)


def _find_lines_in_variance(state):
    """Guess skyline wavelengths from peaks in the summed variance spectrum."""
    from .extraction import medianspec
    medianspec(state, sum=True, all=True)
    wave = np.asarray(state.wave, dtype=float)
    var = np.asarray(state.nmedspec, dtype=float) ** 2
    ok = var > 0
    var, wave = var[ok], wave[ok]
    nw = wave.size
    if nw < 50:
        return np.zeros(0)
    from scipy.ndimage import uniform_filter1d
    var_smth = uniform_filter1d(var, 20, mode="nearest")
    var_sub = var - var_smth
    perc_10, perc_80 = utils.percentile(var_sub, 10), utils.percentile(var_sub, 80)
    mode, rms = None, None
    nbins = int(nw / 15.0)
    for _ in range(6):
        if nbins < 5 or not np.isfinite(nbins):
            break
        hh, edges = np.histogram(var_sub, bins=int(nbins), range=(perc_10, perc_80))
        bb = 0.5 * (edges[:-1] + edges[1:])
        peak = int(np.argmax(hh))
        half = 0.5 * hh[peak]
        i = peak
        while i > 0 and hh[i] >= half:
            i -= 1
        lo = i
        i = peak
        while i < hh.size - 1 and hh[i] >= half:
            i += 1
        hi = i
        fwhm_est = bb[min(hi, bb.size - 1)] - bb[max(lo, 0)]
        binsize = fwhm_est / 10.0 if fwhm_est > 0 else (perc_80 - perc_10) / nbins
        nbins = (perc_80 - perc_10) / binsize if binsize > 0 else nbins
        _, params, _ = mpfitpeak(bb, hh.astype(float), nterms=3, positive=True)
        if np.all(np.isfinite(params)):
            mode, rms = params[1], abs(params[2])
    if mode is None or rms is None or rms <= 0:
        return np.zeros(0)
    starts, ends = [], []
    inline = False
    for i in range(10, nw - 10):
        if var_sub[i] - mode >= 10.0 * rms:
            if not inline:
                inline = True
                starts.append(i)
        elif inline:
            inline = False
            ends.append(i)
            if ends[-1] - starts[-1] > 15:      # too broad to be a real line
                starts.pop()
                ends.pop()
    n = min(len(starts), len(ends))
    if n == 0:
        return np.zeros(0)
    return (wave[np.array(starts[:n])] + wave[np.array(ends[:n])]) / 2.0


def linefit_skylines(state) -> None:
    """Fit sky lines in the noise cube to derive ``R(lambda)`` (port of
    ``kubeviz_linefit_skylines``). Sets ``state.instrres_varpoly`` on success."""
    utils.info("Computing the polynomial coefficients for the instrumental resolution.")
    utils.info("This might take a few minutes, but will only happen once..")

    ohfile = os.path.join(DATA_DIR, "ohspec", "kmos_oh_spec.fits")
    ohspec, hdroh = fits.getdata(ohfile, header=True)
    ohspec = np.asarray(ohspec, dtype=float).ravel()
    crval = hdroh["CRVAL1"] * 1.0e4
    cdelt = hdroh["CDELT1"] * 1.0e4
    crpix = hdroh["CRPIX1"]
    lamboh = (np.arange(ohspec.size) + 1 - crpix) * cdelt + crval

    wave_all = np.asarray(state.wave, dtype=float)
    noverlap = np.sum((lamboh >= wave_all[0]) & (lamboh <= wave_all[-1]))
    # IDL compared the two wavelength grids element-wise; require the template to
    # cover at least half of the data range
    if noverlap > 0 and (lamboh.max() >= wave_all[0] and lamboh.min() <= wave_all[-1]) and \
            (min(lamboh.max(), wave_all[-1]) - max(lamboh.min(), wave_all[0])) > 0.5 * (wave_all[-1] - wave_all[0]):
        skywave_all = _find_oh_lines(state, ohspec, lamboh)
        utils.debug(f"The number of identified OH lines is: {skywave_all.size}")
    else:
        if np.sum(wave_all > 6000) < wave_all.size / 2.0:
            utils.info("The spectrum is too blue for a fit of the instr. resolution from skylines.")
            utils.info("Input the instrumental resolution manually.")
            return
        utils.warn("The template OH spec does not overlap enough with the data spectral range.")
        utils.info("The position of the skylines is guessed automatically.")
        skywave_all = _find_lines_in_variance(state)
        if skywave_all.size == 0:
            utils.warn("The number of identified lines in the variance spec is zero.")
            utils.info("Input the instrumental resolution manually.")
            return

    skylinepos, skylineres = [], []
    step = max(1, int(round(np.sqrt(state.Ncol * state.Nrow) / 25.0)))
    for col in range(0, state.Ncol, step):
        for row in range(0, state.Nrow, step):
            var = np.asarray(state.noise[:, row, col], dtype=float) ** 2
            keep = var > 0
            if not np.any(keep):
                continue
            var, wave = var[keep], wave_all[keep]
            if wave.size < 10:
                continue
            lo, hi = wave[int(wave.size / 100.0 * 5.0)], wave[int(wave.size / 100.0 * 95.0)]
            ok = (skywave_all > lo) & (skywave_all < hi)
            if not np.any(ok):
                continue
            skywave = skywave_all[ok]
            skyindx = utils.closest(wave, skywave)
            for i in range(skywave.size):
                a, b = skyindx[i] - 7, skyindx[i] + 7
                if a < 0 or b >= wave.size:
                    continue
                xw, yv = wave[a:b + 1], var[a:b + 1]
                _, A, _ = mpfitpeak(xw, yv, nterms=4, positive=True)
                if not np.all(np.isfinite(A)):
                    continue
                model = A[3] + A[0] * np.exp(-0.5 * ((xw - A[1]) / A[2]) ** 2)
                resid = np.sum(np.abs(yv - model)) / np.sum(yv) * 100.0
                if resid < 10 and abs(skywave[i] - A[1]) < 2 * state.dlambda and A[2] > 0:
                    skylinepos.append(A[1])
                    skylineres.append(A[1] / (A[2] * SIGTOFWHM))
    skylinepos = np.array(skylinepos)
    skylineres = np.array(skylineres)

    pos_bin, res_bin, sig_bin = [], [], []
    for lam in skywave_all:
        this = np.flatnonzero(np.abs(skylinepos - lam) < 2 * state.dlambda)
        if this.size > 5:
            thisres, _ = utils.sigma_clip(skylineres[this], nsig=2)
            if thisres.size == 0:
                continue
            pos_bin.append(lam)
            res_bin.append(utils.idl_median(thisres, even=True))
            sig_bin.append(0.5 * (utils.percentile(thisres, 84) - utils.percentile(thisres, 16)) / np.sqrt(thisres.size))
    if len(pos_bin) > 1:
        pos_bin, res_bin, sig_bin = map(np.array, (pos_bin, res_bin, sig_bin))
        order = min(4, len(pos_bin) - 1)
        w = np.ones_like(sig_bin)
        pos = sig_bin > 0
        w[pos] = 1.0 / sig_bin[pos]
        coeff = np.polynomial.polynomial.polyfit(pos_bin, res_bin, order, w=w)
        state.instrres_varpoly[:] = 0.0
        state.instrres_varpoly[:order + 1] = coeff
        state.instrres_mode = INSTRRES_VARPOLY
        utils.info(f"Instrumental resolution fitted from {len(pos_bin)} sky lines.")
    else:
        utils.warn("Not enough sky lines measured to fit the instrumental resolution.")
