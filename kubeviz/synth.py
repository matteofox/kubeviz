"""Synthetic datacubes with known emission-line parameters, for tests and demos.

:func:`make_synthetic_cube` builds a MUSE-like cube (DATA + STAT extensions, variance
noise) of a rotating disc galaxy emitting Halpha + [NII] (optionally [SII]) on top of a
faint continuum, and returns the truth maps used to generate it.
"""
from __future__ import annotations

import numpy as np
from astropy.io import fits

from .constants import CKMS, SIGTOFWHM
from .linesdb import rest_wavelength

MUSE_RES_POLY = [499.446, -0.201608, 0.000130290, -7.87479e-09, 0.0]


def muse_resolution(lam):
    return np.polynomial.polynomial.polyval(np.asarray(lam, dtype=float), MUSE_RES_POLY)


def make_synthetic_cube(fname: str | None = None, nx: int = 12, ny: int = 10, redshift: float = 0.02,
                        wave0: float = 6000.0, dwave: float = 1.25, nwave: int = 800,
                        noise_level: float = 2.0, cont_level: float = 5.0, flux_peak: float = 800.0,
                        vmax: float = 150.0, sigma0: float = 60.0, nii_ratio: float = 0.35,
                        sii: bool = False, seed: int = 1, nan_corner: bool = False,
                        broad: bool = False):
    """Return ``(hdulist, truth)``; also writes ``fname`` when given.

    truth keys: ``flux_ha, flux_n2b, flux_n2r, vel, sigma, cont`` (2D maps), ``wave``.
    """
    rng = np.random.default_rng(seed)
    wave = wave0 + dwave * np.arange(nwave)
    yy, xx = np.mgrid[:ny, :nx]
    cx, cy = (nx - 1) / 2.0, (ny - 1) / 2.0
    r = np.hypot(xx - cx, yy - cy)
    flux_ha = flux_peak * np.exp(-0.5 * (r / (0.25 * min(nx, ny))) ** 2)
    vel = vmax * np.arctan((xx - cx) / 2.0) / (np.pi / 2) * np.cos(0.0) * 0.8 + 10.0 * (yy - cy) / ny
    sigma = sigma0 + 20.0 * np.exp(-0.5 * (r / 2.0) ** 2)
    cont = cont_level * (1.0 + 0.3 * np.exp(-0.5 * (r / (0.4 * min(nx, ny))) ** 2))

    lines = {"Ha": flux_ha, "n2_r": nii_ratio * flux_ha, "n2_b": nii_ratio * flux_ha / 3.071}
    if sii:
        lines["s2_b"] = 0.2 * flux_ha
        lines["s2_r"] = 0.15 * flux_ha
    cube = np.zeros((nwave, ny, nx))
    truth = {"vel": vel, "sigma": sigma, "cont": cont, "wave": wave, "redshift": redshift}
    for name, fmap in lines.items():
        truth["flux_" + name] = fmap
        lam0 = rest_wavelength(name) * (1.0 + redshift)
        xcen = lam0 * (1.0 + vel / CKMS)
        width = lam0 * sigma / CKMS
        instr = xcen / (SIGTOFWHM * muse_resolution(xcen))
        tot2 = width ** 2 + instr ** 2
        prof = fmap[None] / np.sqrt(2 * np.pi * tot2)[None] * np.exp(-0.5 * (wave[:, None, None] - xcen[None]) ** 2 / tot2[None])
        cube += prof
        if broad:
            bw = lam0 * (3.0 * sigma) / CKMS
            bt2 = bw ** 2 + instr ** 2
            bflux = 0.3 * fmap
            cube += bflux[None] / np.sqrt(2 * np.pi * bt2)[None] * np.exp(-0.5 * (wave[:, None, None] - xcen[None]) ** 2 / bt2[None])
            truth["bflux_" + name] = bflux
    cube += cont[None] * (1.0 + 0.05 * (wave[:, None, None] - wave.mean()) / np.ptp(wave))
    noise = np.full(cube.shape, noise_level) * (1.0 + 0.2 * rng.random((nwave, 1, 1)))
    cube_noisy = cube + rng.standard_normal(cube.shape) * noise
    if nan_corner:
        cube_noisy[:, 0, 0] = np.nan
        noise[:, 0, 0] = np.nan
    truth["cube_noiseless"] = cube

    pri = fits.PrimaryHDU()
    pri.header["INSTRUME"] = "MUSE"
    pri.header["HIERARCH ESO INS MODE"] = "WFM-NOAO-N"
    pri.header["DATE-OBS"] = "2020-01-01T00:00:00.000"
    dhdr = fits.Header()
    dhdr["EXTNAME"] = "DATA"
    dhdr["CTYPE1"], dhdr["CTYPE2"], dhdr["CTYPE3"] = "RA---TAN", "DEC--TAN", "AWAV"
    dhdr["CUNIT1"], dhdr["CUNIT2"], dhdr["CUNIT3"] = "deg", "deg", "Angstrom"
    dhdr["CRPIX1"], dhdr["CRPIX2"], dhdr["CRPIX3"] = 1.0, 1.0, 1.0
    dhdr["CRVAL1"], dhdr["CRVAL2"], dhdr["CRVAL3"] = 150.0, 2.0, wave0
    dhdr["CD1_1"], dhdr["CD1_2"], dhdr["CD2_1"], dhdr["CD2_2"] = -5.5e-5, 0.0, 0.0, 5.5e-5
    dhdr["CD3_3"] = dwave
    dhdr["BUNIT"] = "10**(-20)*erg/s/cm**2/Angstrom"
    data = fits.ImageHDU(cube_noisy.astype(np.float32), header=dhdr, name="DATA")
    shdr = dhdr.copy()
    shdr["EXTNAME"] = "STAT"
    stat = fits.ImageHDU((noise ** 2).astype(np.float32), header=shdr, name="STAT")
    hdul = fits.HDUList([pri, data, stat])
    if fname:
        hdul.writeto(fname, overwrite=True)
    return hdul, truth
