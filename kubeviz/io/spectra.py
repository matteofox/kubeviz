"""Saving spectra, cubes and images (ports of ``kubeviz_savespec``,
``kubeviz_savecube`` and ``kubeviz_saveimage``)."""
from __future__ import annotations

import numpy as np
from astropy.io import fits

from .. import utils
from ..constants import (CUBE_SN, SPEC_MEDIAN, SPEC_MEDSUB, SPEC_OPTIMAL, SPEC_SLICE, SPEC_SUM,
                         SPEC_WAVG)
from ..core.extraction import getspec, image_for_imgmode
from .results import _trim_origin_radec


def savespec(state, fname: str) -> None:
    """Spectrum in extension 0 (DATA) and noise in extension 1 (NOISE)."""
    spec, nspec = getspec(state)
    ih = state.indatahead if state.indatahead is not None else fits.Header()
    hdr = fits.Header()
    hdr["EXT_TY"] = ("DATA", "The type of the spectrum in this extension")
    hdr["CRVAL1"] = (float(state.wave[0]), "[A] Wavelength at ref. pixel")
    hdr["CRPIX1"] = (1, "[pix] Reference pixel in x")
    hdr["CDELT1"] = (float(state.dlambda), "[A] Spectral resolution")
    hdr["CTYPE1"] = (ih.get("CTYPE3", "WAVE"), "Coordinate system of x-axis")
    hdr["CUNIT1"] = ("Angstrom", "Unit of x-axis")
    sptype = {SPEC_SLICE: ("SPAXEL", "Single spaxel spectrum from the cube"),
              SPEC_SUM: ("SUM", "Sum spectrum from the spaxel mask"),
              SPEC_MEDIAN: ("MEDIAN", "Median spectrum from the spaxel mask"),
              SPEC_WAVG: ("WAVERAGE", "W.average spectrum from the spaxel mask"),
              SPEC_MEDSUB: ("MEDSUB", "Median subtracted spaxel spectrum"),
              SPEC_OPTIMAL: ("OPTIMAL", "Optimal extraction from the spaxel mask")}[state.specmode]
    hdr["SPTYPE"] = sptype
    if state.specmode in (SPEC_SLICE, SPEC_MEDSUB):
        hdr["X_ORIG"] = (state.col + 1, "X position in the cube (1 based)")
        hdr["Y_ORIG"] = (state.row + 1, "Y position in the cube (1 based)")
    nhdr = hdr.copy()
    nhdr["EXT_TY"] = ("NOISE", "The type of the spectrum in this extension")
    hdul = fits.HDUList([fits.PrimaryHDU(np.asarray(spec, dtype=np.float32), header=hdr),
                         fits.ImageHDU(np.asarray(nspec, dtype=np.float32), header=nhdr)])
    hdul.writeto(fname, overwrite=True)
    utils.info(f"Wrote spectrum to file: {fname}")


def _common_header(state, base: fits.Header, extname: str, twod: bool) -> fits.Header:
    hdr = base.copy() if base is not None else fits.Header()
    ih = state.indatahead if state.indatahead is not None else fits.Header()
    ra, dec = _trim_origin_radec(state)
    for k in ("SIMPLE", "EXTEND", "XTENSION"):
        if k in hdr:
            del hdr[k]
    hdr["EXTNAME"] = (extname, "This extension contains data values")
    hdr["DATACUBE"] = (state.filename, "Name of the original datacube")
    hdr["INSTRUME"] = (state.instr.upper(), "Instrument name")
    hdr["BAND"] = (state.band.upper(), "Band/Filter name")
    hdr["SPATSMTH"] = (int(state.smooth), "Spatial smoothing applied")
    hdr["SPECSMTH"] = (int(state.specsmooth), "Spectral smoothing applied")
    for k, c in (("CTYPE1", "TAN projection used"), ("CTYPE2", "TAN projection used"),
                 ("CD1_1", "[] x-component of East"), ("CD2_1", "[] x-component of North"),
                 ("CD1_2", "[] y-component of East"), ("CD2_2", "[] y-component of North"),
                 ("CDELT1", "[deg] Pixel resolution in x"), ("CDELT2", "[deg] Pixel resolution in y"),
                 ("CUNIT1", "Unit of x-axis"), ("CUNIT2", "Unit of y-axis")):
        if k in ih:
            hdr[k] = (ih[k], c)
    hdr["CRPIX1"] = (1, "[pix] Reference pixel in x")
    hdr["CRPIX2"] = (1, "[pix] Reference pixel in y")
    hdr["CRVAL1"] = (ra, "[deg] RA at ref. pixel")
    hdr["CRVAL2"] = (dec, "[deg] DEC at ref. pixel")
    hdr["STARTCOL"] = (int(state.Startcol) + 1, "Index of the first col in the original frame")
    hdr["STARTROW"] = (int(state.Startrow) + 1, "Index of the first row in the original frame")
    if twod:
        for k in ("NAXIS3", "CRPIX3", "CRVAL3", "CDELT3", "CTYPE3", "CUNIT3", "CD3_3"):
            if k in hdr:
                del hdr[k]
    else:
        hdr["CTYPE3"] = (ih.get("CTYPE3", "WAVE"), "Coordinate system of z-axis")
        hdr["CRPIX3"] = (1, "[pix] Reference pixel in z")
        hdr["CRVAL3"] = (float(state.wave[0]), "[A] Wavelength at ref. pixel")
        hdr["CDELT3"] = (float(state.dlambda), "[A] Spectral resolution")
        hdr["CUNIT3"] = ("Angstrom", "Unit of z-axis")
    return hdr


def _dtype(state):
    return np.float64 if (state.indatahead is not None and state.indatahead.get("BITPIX", -32) == -64) else np.float32


def savecube(state, fname: str) -> None:
    """Current (smoothed) data and noise cubes in extensions DATA and STAT-like."""
    dt = _dtype(state)
    dhdr = _common_header(state, state.indatahead, "DATA", twod=False)
    nhdr = _common_header(state, state.innoisehead, "NOISE", twod=False)
    noise = state.noisecube ** 2 if state.noiseisvar else state.noisecube
    hdul = fits.HDUList([fits.PrimaryHDU(),
                         fits.ImageHDU(np.asarray(state.datacube, dtype=dt), header=dhdr),
                         fits.ImageHDU(np.asarray(noise, dtype=dt), header=nhdr)])
    hdul.writeto(fname, overwrite=True)
    utils.info(f"Wrote cube to file: {fname}")


def saveimage(state, fname: str, data=None, noise=None) -> None:
    """Current image and its noise. When ``data`` is None the image implied by
    ``cubesel``/``imgmode`` (or the line-fit result image) is used."""
    if data is None:
        if state.cubesel <= CUBE_SN:
            data, noise = image_for_imgmode(state)
        else:
            if state.flagmode:
                data, noise = state.lineresimg, state.lineerrresimg
            else:
                data, noise = state.unm_lineresimg, state.unm_lineerrresimg
    dt = _dtype(state)
    dhdr = _common_header(state, state.indatahead, "DATA", twod=True)
    nhdr = _common_header(state, state.innoisehead, "NOISE", twod=True)
    noise_out = np.asarray(noise, dtype=float) ** 2 if state.noiseisvar else np.asarray(noise, dtype=float)
    hdul = fits.HDUList([fits.PrimaryHDU(),
                         fits.ImageHDU(np.asarray(data, dtype=dt), header=dhdr),
                         fits.ImageHDU(np.asarray(noise_out, dtype=dt), header=nhdr)])
    hdul.writeto(fname, overwrite=True)
    utils.info(f"Wrote data and noise images to file: {fname}")
