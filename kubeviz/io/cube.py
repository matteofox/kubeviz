"""Datacube loading (port of ``kubeviz_getdata`` and ``kubeviz_decodetrimstr``).

Instrument detection, extension defaults, trimming, flux scaling, variance -> sigma
conversion, wavelength axis and the optional moment-window / Gaussian initial-guess
maps are handled here. Cubes are returned as ``(Nwpix, Nrow, Ncol)`` float arrays.
"""
from __future__ import annotations

import os
import re

import numpy as np
from astropy.io import fits

from .. import utils
from ..core.smooth import smooth_state_datacube


class CubeLoadError(RuntimeError):
    pass


# ----------------------------------------------------------------------------- helpers
def splitpath(fullname: str) -> tuple[str, str]:
    """IDL ``kubeviz_splitpath``: (basename, dir with trailing separator or '')."""
    fname = os.path.basename(fullname)
    d = os.path.dirname(fullname)
    return fname, (d + os.sep if d else "")


def decode_trimstr(trimstr: str, shape_xyz) -> tuple[list[int] | None, tuple[int, int, int]]:
    """Decode ``'[x1:x2,y1:y2,z1:z2]'`` (1-based, inclusive, ``*`` wildcard).

    Returns ``(P, (startcol, startrow, startwpix))`` where ``P`` is a list of six
    0-based inclusive indices, or ``None`` when the syntax is not understood.
    """
    x_size, y_size, z_size = [int(v) for v in shape_xyz]
    s = trimstr.strip()
    if s.startswith("[") and s.endswith("]"):
        s = s[1:-1]
    parts = [p.strip() for p in s.split(",")]
    if len(parts) != 3:
        utils.warn("Trim keyword syntax not understood. No trim will take place")
        return None, (0, 0, 0)
    out = []
    for part, size in zip(parts, (x_size, y_size, z_size)):
        if part == "*":
            out.extend([0, size - 1])
            continue
        m = re.fullmatch(r"(\d+)\s*:\s*(\d+)", part)
        if not m:
            utils.warn("Trim keyword syntax not understood. No trim will take place")
            return None, (0, 0, 0)
        a, b = int(m.group(1)), int(m.group(2))
        if a < 1 or b > size or b < a:
            utils.warn("Trim keyword syntax not understood. No trim will take place")
            return None, (0, 0, 0)
        out.extend([a - 1, b - 1])
    return out, (out[0], out[2], out[4])


def _hdr_get(hdr, *keys, default=None):
    """First present header keyword among ``keys`` (astropy strips HIERARCH)."""
    if hdr is None:
        return default
    for k in keys:
        if k in hdr:
            return hdr[k]
        hk = "HIERARCH " + k
        if hk in hdr:
            return hdr[hk]
    return default


def _read_section(hdu, P):
    """Read the (possibly trimmed) 3D section of an image HDU as float."""
    if P is None:
        data = hdu.data
    else:
        x1, x2, y1, y2, z1, z2 = P
        data = hdu.section[z1:z2 + 1, y1:y2 + 1, x1:x2 + 1]
    return np.array(data, dtype=np.float64 if hdu.header.get("BITPIX", -32) == -64 else np.float32)


def detect_instrument(prihead) -> tuple[str, bool]:
    """Instrument name (lower case) and whether it works in vacuum."""
    instrume = str(_hdr_get(prihead, "INSTRUME", default="")).strip()
    instrmnt = str(_hdr_get(prihead, "INSTRMNT", default="")).strip()
    if instrume == "KMOS":
        version = _hdr_get(prihead, "VERSION", default=0)
        try:
            v = float(version)
        except (TypeError, ValueError):
            v = 0.0
        if v > 0:
            utils.info("Detected instrument: KMOS")
            utils.info(f"Detected data product from KMOS3D data release {v:4.1f}")
            return "kmos3d", True
        utils.info("Detected instrument: KMOS")
        return "kmos", True
    if instrume == "VIMOS":
        utils.info("Detected instrument: VIMOS")
        return "vimos", False
    if instrume == "SINFONI":
        utils.info("Detected instrument: SINFONI")
        return "sinfoni", True
    if instrume == "MUSE":
        utils.info("Detected instrument: MUSE")
        return "muse", False
    if instrume == "SAMI":
        utils.info("Detected instrument: SAMI")
        return "sami", False
    if instrmnt == "WiFeS":
        utils.info("Detected instrument: WiFeS")
        return "wifes", False
    return instrume.lower(), False


def _default_ext(nextend: int, which: str) -> int:
    if which == "data":
        table = {0: 0, 1: 1, 2: 1}
    else:
        table = {0: 0, 1: 1, 2: 2}
    if nextend not in table:
        raise CubeLoadError(f"FITS data format unrecognized. Please specify the {which} extension.")
    return table[nextend]


def _load_aux_map(path: str, ncol: int, nrow: int, nplanes_expected: int, label: str):
    """Load a (Nplanes, Nrow, Ncol) auxiliary map, adapting the spatial size with [0,0] rigid."""
    if path is None:
        return None
    m = fits.getdata(path, ext=0)
    m = np.asarray(m, dtype=float)
    if m.ndim != 3 or m.shape[0] != nplanes_expected:
        utils.warn(f"File {path} contains wrong number of image planes")
        return None
    mp, mrow, mcol = m.shape
    if (mcol, mrow) != (ncol, nrow):
        utils.warn(f"File {path} image size does not correspond to spatial cube size! Assuming [0,0] is rigid")
        out = np.full((mp, nrow, ncol), -999.0)
        r, c = min(mrow, nrow), min(mcol, ncol)
        out[:, :r, :c] = m[:, :r, :c]
        m = out
    return m


# ----------------------------------------------------------------------------- main entry
def load_cube(state, datafile: str, noisefile: str | None = None, ext: int | None = None,
              noise_ext: int | None = None, trim: str | None = None, logarithmic: bool = False,
              waveunit: str | None = None, momwindowfile: str | None = None,
              gaussinitfile: str | None = None) -> None:
    """Read data and noise cubes into ``state`` (port of ``kubeviz_getdata``).

    ``datafile``/``noisefile`` are full paths. Instrument specific defaults are
    applied unless ``state.instr``, ``state.band``, ``state.fluxfac`` were preset.
    """
    dataname, datadir = splitpath(datafile)
    if noisefile is None:
        noisefile = datafile
    noisename, noisedir = splitpath(noisefile)
    if not noisedir:
        noisedir = datadir
    if not os.path.exists(datafile):
        raise CubeLoadError(f"Unable to locate file: {datafile}")
    if not os.path.exists(noisedir + noisename):
        raise CubeLoadError(f"Unable to locate file: {noisedir + noisename}")

    state.filename = dataname
    state.indir = datadir

    with fits.open(datafile, memmap=True) as fcb, fits.open(noisedir + noisename, memmap=True) as fcbnoi:
        prihead = fcb[0].header
        nextend = len(fcb) - 1
        nextend_noi = len(fcbnoi) - 1

        if not state.instr:
            state.instr, vac = detect_instrument(prihead)
            if vac:
                state.vacuum = True
        else:
            state.instr = state.instr.lower()

        if state.instr == "wifes":
            ext, noise_ext = 0, 1
        if state.instr == "kmos3d":
            ext, noise_ext = 1, 2
        if ext is None:
            ext = _default_ext(nextend, "data")
        if noise_ext is None:
            noise_ext = _default_ext(nextend_noi, "noise")

        hdu = fcb[ext]
        nhdu = fcbnoi[noise_ext]
        ndim_data = hdu.header.get("NAXIS", 0)
        ndim_noi = nhdu.header.get("NAXIS", 0)

        P = None
        state.Startcol = state.Startrow = state.Startwpix = 0
        if trim is not None and ndim_data == 3:
            nx, ny, nz = hdu.header["NAXIS1"], hdu.header["NAXIS2"], hdu.header["NAXIS3"]
            P, (state.Startcol, state.Startrow, state.Startwpix) = decode_trimstr(trim, (nx, ny, nz))

        utils.info(f"Loading data cube: {datafile}")
        if ndim_data == 3 and ndim_noi == 3:
            cube = _read_section(hdu, P)
            noisecube = _read_section(nhdu, P)
        elif ndim_data == 1 and ndim_noi == 1:
            cube = np.asarray(hdu.data, dtype=float).reshape(-1, 1, 1)
            noisecube = np.asarray(nhdu.data, dtype=float).reshape(-1, 1, 1)
        else:
            raise CubeLoadError("Either the input data or the noise is not a 3D datacube. Abort.")
        hdr = hdu.header.copy()
        noisehdr = nhdu.header.copy()

    # ------------------------------------------------------------------ instrument specifics
    instr = state.instr
    if instr == "kmos3d":
        if not state.band:
            state.band = str(_hdr_get(prihead, "OBSBAND", default="")).strip()
            utils.info(f"Detected band: {state.band}")
        state.ifu = -1
        if state.fluxfac == 1:
            state.fluxfac = 0.1
            utils.info("Input unit: 1E-17 W/m^2/um -> Output unit: 1E-17 erg/cm^2/s/A ")
    elif instr == "kmos":
        if not state.band:
            state.band = str(_hdr_get(prihead, "ESO INS FILT1 NAME", default="")).strip()
            utils.info(f"Detected band: {state.band}")
        if state.ifu == 0:
            for j in range(1, 25):
                name = _hdr_get(hdr, f"ESO OCS ARM{j} NAME", default="")
                if str(name).strip() > "0":
                    state.ifu = j
                    utils.info(f"Detected ifu: {j}")
                    break
        else:
            state.ifu = -1
        if state.fluxfac == 1:
            state.fluxfac = 0.1
            utils.info("Input unit: W/m^2/um -> Output unit: erg/cm^2/s/A ")
    elif instr == "sinfoni":
        if not state.band:
            state.band = str(_hdr_get(prihead, "ESO INS FILT1 NAME", default="")).strip()
            utils.info(f"Detected band: {state.band}")
        if state.pixscale == 0:
            cd1 = _hdr_get(prihead, "CDELT1", default=0.0)
            state.pixscale = int(round(2.0 * 3600.0 * 1000.0 * abs(float(cd1))))
            utils.info(f"Detected pixel scale: {state.pixscale} mas")
        if state.fluxfac == 1:
            state.fluxfac = 0.1
            utils.info("Input unit: W/m^2/um -> Output unit: erg/cm^2/s ")
    elif instr == "muse":
        state.noiseisvar = True
        if not state.band:
            state.band = str(_hdr_get(prihead, "ESO INS MODE", default="")).strip()
            utils.info(f"Detected mode: {state.band}")
        if state.fluxfac == 1:
            utils.info("Output unit: 1E-20 erg/cm^2/s ")
    elif instr == "sami":
        state.noiseisvar = True
    elif instr == "wifes":
        state.noiseisvar = True
        if not state.band:
            state.band = str(_hdr_get(prihead, "WIFESARM", default="")).strip()
            utils.info(f"Detected arm: {state.band}")
        if state.fluxfac == 1:
            state.fluxfac = 1e18
            utils.info("Output unit: 1E-18 erg/cm^2/s ")
    else:
        utils.info(f"Instrument: {instr} requires manual input of band and fluxfac parameters, if required.")

    if state.noiseisvar:
        with np.errstate(invalid="ignore"):
            noisecube = np.sqrt(noisecube)

    state.indatacube = cube * state.fluxfac
    state.innoisecube = noisecube * state.fluxfac
    state.indatahead = hdr
    state.innoisehead = noisehdr
    state.inprihead = prihead.copy() if ext != 0 else None

    state.Nwpix, state.Nrow, state.Ncol = state.indatacube.shape

    # ------------------------------------------------------------------ smoothing + bad pixels
    smooth_state_datacube(state)

    # ------------------------------------------------------------------ wavelength axis
    unit = str(hdr.get("CUNIT3", "")).replace(" ", "")
    if unit == "Angstrom":
        waveunit = "ANGSTROMS"
    elif unit in ("um", "MICRON"):
        waveunit = "MICRONS"
    wave_conv = 1.0
    if waveunit is not None:
        wu = waveunit.upper()
        if wu in ("ANGSTROMS", "ANGSTROM", "A"):
            wave_conv = 1.0
        elif wu in ("MICRONS", "MICRON", "UM"):
            wave_conv = 1.0e4
        elif wu == "NM":
            wave_conv = 10.0
        else:
            utils.warn(f"No unit {waveunit} defined")
        utils.debug(f"Units: {waveunit}; conversion {wave_conv}")

    ax = "1" if ndim_data == 1 else "3"
    state.pix0 = float(hdr.get("CRPIX" + ax, 0.0))
    state.lambda0 = float(hdr.get("CRVAL" + ax, 0.0)) * wave_conv
    cdelt = float(hdr.get("CDELT" + ax, 0.0))
    if cdelt > 0:
        state.dlambda = cdelt * wave_conv
    else:
        state.dlambda = float(hdr.get(f"CD{ax}_{ax}", 0.0)) * wave_conv

    if state.lambda0 == 0.0:
        utils.warn("Cannot find wavelength solution in header.")
        utils.warn("Will use pixel number instead of lambda.")
        state.wave = np.arange(state.Nwpix, dtype=float)
    else:
        pix = (np.arange(state.Nwpix, dtype=float) + 1.0) - state.pix0 + state.Startwpix
        if logarithmic:
            state.wave = np.exp(state.lambda0 + state.dlambda * pix)
        else:
            state.wave = state.lambda0 + state.dlambda * pix

    # ------------------------------------------------------------------ transpose
    if state.transpose:
        utils.info("Transposing the datacube...")
        state.datacube = np.ascontiguousarray(np.transpose(state.datacube, (0, 2, 1)))
        state.noisecube = np.ascontiguousarray(np.transpose(state.noisecube, (0, 2, 1)))
        state.indatacube = np.ascontiguousarray(np.transpose(state.indatacube, (0, 2, 1)))
        state.innoisecube = np.ascontiguousarray(np.transpose(state.innoisecube, (0, 2, 1)))
        state.Nwpix, state.Nrow, state.Ncol = state.datacube.shape
        state.badpixelmask = np.ascontiguousarray(np.transpose(state.badpixelmask, (0, 2, 1)))
        state.badpixelimg = np.ascontiguousarray(state.badpixelimg.T)

    state.noise = state.noisecube

    # ------------------------------------------------------------------ auxiliary maps
    state.mom_windowmap = _load_aux_map(momwindowfile, state.Ncol, state.Nrow, 2, "momwindow")
    # WARNING (IDL): hardcoded for 5 planes (dv, sigma, 3 fluxes)
    state.gauss_initmap = _load_aux_map(gaussinitfile, state.Ncol, state.Nrow, 5, "gaussinit")

    # initialise derived images
    state.img1 = np.zeros((state.Nrow, state.Ncol))
    state.img2 = np.zeros((state.Nrow, state.Ncol))
    state.nimg1 = np.zeros((state.Nrow, state.Ncol))
    state.nimg2 = np.zeros((state.Nrow, state.Ncol))
    state.bimg1 = np.zeros((state.Nrow, state.Ncol), dtype=bool)
    state.bimg2 = np.zeros((state.Nrow, state.Ncol), dtype=bool)
    state.col = state.Ncol // 2
    state.row = state.Nrow // 2
    state.wpix = state.Nwpix // 2
    state.fitallrange = np.array([0, state.Ncol - 1, 0, state.Nrow - 1], dtype=int)
