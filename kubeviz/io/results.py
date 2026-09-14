"""Results FITS files, byte-compatible with the IDL ``kubeviz_linefit_saveres`` /
``kubeviz_linefit_loadres`` layout.

Gaussian fits, spaxel mode (``Npl = 3 + Nlines``, ``Ncp = 1 + Nlines``), planes::

    [0, Npl)                     narrow: flag, dv, sig(v), fluxes
    [Npl, 2Npl)                  broad
    [2Npl, 2Npl+Ncp)             continuum: flag, level at each line
    [.., +Npl) [.., +Npl)        narrow +1sig, narrow -1sig
    [.., +Npl) [.., +Npl)        broad +1sig, broad -1sig
    [.., +Ncp) [.., +Ncp)        continuum +1sig, -1sig
    [last]                       reduced chi-square

Mask mode has the same planes without the chi-square, as a 2D image (Nmask x Nplanes).
Moments: per line 6 planes (flag, flux, dv, sigv, skew, kurt), then continuum per line,
then +1sig/-1sig of the 5 moments per line, then continuum +1sig and -1sig.
"""
from __future__ import annotations

import os

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

from .. import utils
from ..constants import ERR_METHOD_NAMES, FIT_GAUSS, MODE_MASK, MODE_SPAXEL, NOT_FIT


def _flag_to_nan(flag):
    out = np.ones(flag.shape, dtype=float)
    out[flag > 0] = np.nan
    return out


def pack_results(state, rs):
    """Return ``(rescube_flagged, rescube_raw)`` as numpy arrays in FITS orientation
    (planes first)."""
    nl = state.Nlines
    npl, ncp = 3 + nl, 1 + nl
    mode = rs.mode
    lead = rs.n.shape[:-1]

    def stack(arrs):
        return np.concatenate([np.moveaxis(a, -1, 0) for a in arrs], axis=0)

    if state.linefit_type == FIT_GAUSS:
        parts = [rs.n, rs.b, rs.c, rs.nerr[..., 0], rs.nerr[..., 1], rs.berr[..., 0], rs.berr[..., 1],
                 rs.cerr[..., 0], rs.cerr[..., 1]]
        raw = stack(parts)
        if mode == MODE_SPAXEL:
            raw = np.concatenate([raw, rs.chisq[None]], axis=0)
        res = raw.copy()
        if mode == MODE_SPAXEL:
            nflag, bflag, cflag = (_flag_to_nan(rs.n[..., 0]), _flag_to_nan(rs.b[..., 0]), _flag_to_nan(rs.c[..., 0]))
            blocks = [(0, npl, nflag), (npl, 2 * npl, bflag), (2 * npl, 2 * npl + ncp, cflag),
                      (2 * npl + ncp, 4 * npl + ncp, nflag), (4 * npl + ncp, 6 * npl + ncp, bflag),
                      (6 * npl + ncp, 6 * npl + 3 * ncp, cflag)]
            for start, end, fl in blocks:
                res[start + 1:end] *= fl[None]
        else:
            flag = _flag_to_nan(rs.n[..., 0])
            res[1:] *= flag[None]
        return res, raw

    # moments
    nplanes = 7 * nl + 12 * nl
    raw = np.zeros((nplanes,) + lead)
    m, merr, c, cerr = np.moveaxis(rs.m, -1, 0), np.moveaxis(rs.merr, -2, 0), np.moveaxis(rs.c, -1, 0), np.moveaxis(rs.cerr, -2, 0)
    for il in range(nl):
        raw[il * 6] = m[il * 6 + 5]
        raw[il * 6 + 1: il * 6 + 6] = m[il * 6: il * 6 + 5]
        raw[6 * nl + il] = c[il + 1]
        raw[7 * nl + 10 * il: 7 * nl + 10 * il + 5] = merr[il * 6: il * 6 + 5, ..., 0]
        raw[7 * nl + 10 * il + 5: 7 * nl + 10 * il + 10] = merr[il * 6: il * 6 + 5, ..., 1]
        raw[17 * nl + il] = cerr[il + 1, ..., 0]
        raw[18 * nl + il] = cerr[il + 1, ..., 1]
    res = raw.copy()
    for il in range(nl):
        bad = rs.m[..., il * 6 + 5] > 0
        for pl in list(range(il * 6 + 1, il * 6 + 6)) + [6 * nl + il] + \
                list(range(7 * nl + 10 * il, 7 * nl + 10 * il + 10)) + [17 * nl + il, 18 * nl + il]:
            res[pl][bad] = np.nan
    return res, raw


def _plane_descriptions(state):
    """P1..PN header cards describing the planes."""
    nl = state.Nlines
    names = state.linenames
    desc = []
    if state.linefit_type == FIT_GAUSS:
        for restype in range(9):
            nres = (3 + nl) if (restype != 2 and restype < 7) else (1 + nl)
            linetype = ["narrow line", "broad line", "continuum", "narrow line", "narrow line",
                        "broad line", "broad line", "continuum", "continuum"][restype]
            datatype = ["", "", "", "plus 1-sig", "minus 1-sig", "plus 1-sig", "minus 1-sig",
                        "plus 1-sig", "minus 1-sig"][restype]
            iscont = restype == 2 or restype >= 7
            for i in range(nres):
                if i == 0:
                    strpar = "flag"
                elif iscont:
                    strpar = f"{names[i - 1]} flux"
                elif i == 1:
                    strpar = "dv"
                elif i == 2:
                    strpar = "sig(v)"
                else:
                    strpar = f"{names[i - 3]} flux"
                desc.append(f"{linetype} {datatype} {strpar}".replace("  ", " "))
        if state.linefit_mode == MODE_SPAXEL:
            desc.append("Reduced chi square, (chisq/dof)")
    else:
        for restype in range(6):
            linetype = ["moment", "continuum", "moment", "moment", "continuum", "continuum"][restype]
            datatype = ["", "", "plus 1-sig", "minus 1-sig", "plus 1-sig", "minus 1-sig"][restype]
            startpar, endpar = [(0, 5), (0, 0), (1, 5), (1, 5), (0, 0), (0, 0)][restype]
            for il in range(nl):
                for i in range(startpar, endpar + 1):
                    if i == 0:
                        strpar = f"{names[il]} flag" if (restype != 1 and restype < 4) else f"{names[il]} flux"
                    else:
                        strpar = f"{names[il]} " + ["0th (flux)", "1st (dv)", "2nd (sigv)",
                                                    "3rd (norm Skewness)", "4th (norm Kurtosis)"][i - 1]
                    desc.append(f"{linetype} {datatype} {strpar}".replace("  ", " "))
    return desc


def results_header(state, rs, redshift_out=-1.0) -> fits.Header:
    hdr = fits.Header()
    hdr["DATACUBE"] = (state.filename, "Name of the original datacube")
    hdr["INSTRUME"] = (state.instr.upper(), "Instrument name")
    hdr["BAND"] = (state.band.upper(), "Band/Filter name")
    hdr["Z_IN"] = (float(state.redshift), "Redshift used to fit the lines")
    hdr["Z_OUT"] = (round(float(redshift_out), 5), "Redshift obtained from the fit")
    hdr["FITTYPE"] = (int(state.linefit_type), "Type of fit: 0: gauss 1: moments")
    hdr["RESTYPE"] = (int(rs.mode), "Type of result: 0: spaxels 1: masks")
    hdr["ERMETHOD"] = (int(state.domontecarlo), "Err method: 0:Noise, 1:Bootstrap 2,3,4:MC1,2,3")
    hdr["SCNOISE"] = (int(state.scaleNoiseerrors), "Scale noise Err. 0:Off 1:On (Only if ERMETH=0)")
    hdr["MCNOISE"] = (int(state.useMonteCarlonoise), "MonteCarlo noise 0:Off 1:On (Only if ERMETH>0)")
    hdr["SPATSMTH"] = (int(state.smooth), "Spatial smoothing applied")
    hdr["SPECSMTH"] = (int(state.specsmooth), "Spectral smoothing applied")
    hdr["LINESET"] = (int(state.selected_lineset), "Lineset used")
    # fit setup (not recorded by the IDL version; needed to reproduce a run)
    hdr["FITRNGB"] = (float(state.maxwoffb), "[A] Fitting range bluewards of the main line")
    hdr["FITRNGR"] = (float(state.maxwoffr), "[A] Fitting range redwards of the main line")
    hdr["CONTMODE"] = (int(state.continuumfit_mode), "Continuum: 0 side bands, 1 fitted with the lines")
    hdr["CONTMINO"] = (float(state.continuumfit_minoff), "[A] Continuum side band min offset")
    hdr["CONTMAXO"] = (float(state.continuumfit_maxoff), "[A] Continuum side band max offset")
    hdr["CONTMINP"] = (float(state.continuumfit_minperc), "Continuum min percentile")
    hdr["CONTMAXP"] = (float(state.continuumfit_maxperc), "Continuum max percentile")
    hdr["CONTORDR"] = (int(state.continuumfit_order), "Continuum polynomial order")
    hdr["FIXRATIO"] = (int(state.fitfixratios), "Fixed line ratios 0:Off 1:On")
    hdr["FITCONST"] = (int(state.fitconstr), "User constraints 0:Off 1:On")
    hdr["SECCOMP"] = (int(state.secondcomp_mode), "2nd component mode 0:free 1:offset 2:width")
    hdr["SNTHRESH"] = (float(state.mask_sn_thresh), "Autoflag S/N threshold")
    hdr["MAXVELER"] = (float(state.mask_maxvelerr), "[km/s] Autoflag max velocity/dispersion error")
    hdr["KVZVERS"] = (state.version, "kubeviz (Python) version")
    if rs.mode == MODE_MASK:
        hdr["NMASK"] = (int(rs.n.shape[0]), "Number of masks stored")
    if rs.mode == MODE_SPAXEL and state.indatahead is not None:
        ih = state.indatahead
        ra, dec = _trim_origin_radec(state)
        hdr["CTYPE1"] = (ih.get("CTYPE1", ""), "TAN projection used")
        hdr["CTYPE2"] = (ih.get("CTYPE2", ""), "TAN projection used")
        hdr["CRPIX1"] = (1, "[pix] Reference pixel in x")
        hdr["CRPIX2"] = (1, "[pix] Reference pixel in y")
        hdr["CRVAL1"] = (ra, "[deg] RA at ref. pixel")
        hdr["CRVAL2"] = (dec, "[deg] DEC at ref. pixel")
        for k, c in (("CD1_1", "[] x-component of East"), ("CD2_1", "[] x-component of North"),
                     ("CD1_2", "[] y-component of East"), ("CD2_2", "[] y-component of North"),
                     ("CDELT1", "[deg] Pixel resolution in x"), ("CDELT2", "[deg] Pixel resolution in y"),
                     ("CUNIT1", "Unit of x-axis"), ("CUNIT2", "Unit of y-axis")):
            if k in ih:
                hdr[k] = (ih[k], c)
        hdr["STARTCOL"] = (int(state.Startcol) + 1, "Index of the first col in the original frame")
        hdr["STARTROW"] = (int(state.Startrow) + 1, "Index of the first row in the original frame")
    for i, d in enumerate(_plane_descriptions(state), start=1):
        hdr[f"P{i}"] = d
    return hdr


def _trim_origin_radec(state):
    """RA/DEC of the (trimmed) origin pixel, as IDL ``xyad`` (0-based pixels)."""
    try:
        w = WCS(state.indatahead).celestial
        ra, dec = w.all_pix2world([[state.Startcol, state.Startrow]], 0)[0]
        return float(ra), float(dec)
    except Exception:
        return 0.0, 0.0


def saveres(state, fname: str, rs=None) -> tuple[str, str]:
    """Write the flagged results to ``fname`` and the unflagged ones to ``*_noflag.fits``."""
    utils.info("Saving FITS results file...")
    if rs is None:
        rs = state.get_results()
    res, raw = pack_results(state, rs)
    hdr = results_header(state, rs)
    rawfname = fname[:-5] + "_noflag.fits" if fname.lower().endswith(".fits") else fname + "_noflag.fits"
    fits.writeto(fname, res, hdr, overwrite=True)
    fits.writeto(rawfname, raw, hdr, overwrite=True)
    utils.info(f"Results written to {fname} and {rawfname}")
    return fname, rawfname


def loadres(state, fname: str) -> bool:
    """Load a kubeviz results file into the matching ``(mode, errmethod)`` container."""
    if not os.path.exists(fname):
        utils.error("Results file does not exist!")
        return False
    rescube, hdr = fits.getdata(fname, header=True)
    rescube = np.asarray(rescube, dtype=float)
    if int(hdr.get("LINESET", -1)) != state.selected_lineset:
        utils.error("The current lineset is not compatible with the saved results!")
        return False
    erm = int(hdr.get("ERMETHOD", -1))
    if erm < 0 or erm > 4:
        utils.error("Unrecognized error method!")
        return False
    if int(hdr.get("SPATSMTH", state.smooth)) != state.smooth or int(hdr.get("SPECSMTH", state.specsmooth)) != state.specsmooth:
        utils.warn("The current spatial or spectral smoothing do not match those for the saved results ")
        utils.warn("Saved results are loaded anyway. ")
    if abs(float(hdr.get("Z_IN", state.redshift)) - state.redshift) > 0.1:
        utils.warn("The current redshift significantly differs from the one stored in the results file. ")
        utils.warn("Saved results are loaded anyway. Use with caution! ")
    restype = int(hdr.get("FITTYPE", 0))
    resmode = int(hdr.get("RESTYPE", 0))
    nl = state.Nlines
    npl, ncp = 3 + nl, 1 + nl

    if resmode == MODE_SPAXEL:
        if rescube.ndim != 3 or rescube.shape[1:] != (state.Nrow, state.Ncol):
            utils.error("Saved results cube does not match current cube dimensions!")
            return False
        nplanes_exp = 6 * npl + 3 * ncp + 1 if restype == FIT_GAUSS else 19 * nl
    else:
        if rescube.ndim != 2:
            utils.error("Saved results cube is corrupted!")
            return False
        nplanes_exp = 6 * npl + 3 * ncp if restype == FIT_GAUSS else 19 * nl
    if rescube.shape[0] != nplanes_exp:
        utils.error("The current number of lines does not match those in the saved results cube!")
        return False

    nmask = rescube.shape[1] if resmode == MODE_MASK else state.Nmask
    if resmode == MODE_MASK and nmask > state.Nmask:
        state.Nmask = nmask
        pad = np.zeros((nmask - state.spaxselect.shape[0], state.Nrow, state.Ncol))
        state.spaxselect = np.concatenate([state.spaxselect, pad], axis=0)
    rs = state.get_results(resmode, erm)
    if resmode == MODE_MASK:
        rs.ensure_nmask(nmask)
        # only overwrite the first nmask entries
        sl = slice(0, nmask)
    else:
        sl = slice(None)

    def planes(a, b):
        return np.moveaxis(rescube[a:b], 0, -1)

    if restype == FIT_GAUSS:
        rs.n[sl] = utils.remove_badvalues(planes(0, npl))
        rs.b[sl] = utils.remove_badvalues(planes(npl, 2 * npl))
        rs.c[sl] = utils.remove_badvalues(planes(2 * npl, 2 * npl + ncp))
        o = 2 * npl + ncp
        rs.nerr[sl][..., 0] = utils.remove_badvalues(planes(o, o + npl), NOT_FIT)
        rs.nerr[sl][..., 1] = utils.remove_badvalues(planes(o + npl, o + 2 * npl), NOT_FIT)
        rs.berr[sl][..., 0] = utils.remove_badvalues(planes(o + 2 * npl, o + 3 * npl), NOT_FIT)
        rs.berr[sl][..., 1] = utils.remove_badvalues(planes(o + 3 * npl, o + 4 * npl), NOT_FIT)
        rs.cerr[sl][..., 0] = utils.remove_badvalues(planes(o + 4 * npl, o + 4 * npl + ncp), NOT_FIT)
        rs.cerr[sl][..., 1] = utils.remove_badvalues(planes(o + 4 * npl + ncp, o + 4 * npl + 2 * ncp), NOT_FIT)
        if resmode == MODE_SPAXEL:
            rs.chisq[...] = utils.remove_badvalues(rescube[-1])
    else:
        for il in range(nl):
            rs.m[sl][..., il * 6 + 5] = utils.remove_badvalues(rescube[il * 6])
            rs.m[sl][..., il * 6: il * 6 + 5] = utils.remove_badvalues(planes(il * 6 + 1, il * 6 + 6))
            rs.c[sl][..., il + 1] = utils.remove_badvalues(rescube[6 * nl + il])
            rs.merr[sl][..., il * 6: il * 6 + 5, 0] = utils.remove_badvalues(planes(7 * nl + 10 * il, 7 * nl + 10 * il + 5), NOT_FIT)
            rs.merr[sl][..., il * 6: il * 6 + 5, 1] = utils.remove_badvalues(planes(7 * nl + 10 * il + 5, 7 * nl + 10 * il + 10), NOT_FIT)
            rs.cerr[sl][..., il + 1, 0] = utils.remove_badvalues(rescube[17 * nl + il], NOT_FIT)
            rs.cerr[sl][..., il + 1, 1] = utils.remove_badvalues(rescube[18 * nl + il], NOT_FIT)

    utils.info("Results file successfully loaded! ")
    if erm != state.domontecarlo:
        utils.warn("Loaded results have been produced with an error method that is not currently set.")
        utils.warn(f"Switch to \"{ERR_METHOD_NAMES[erm]}\" errors to see the results.")
    state.linefit_type = restype
    return True
