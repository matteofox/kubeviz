"""The multi-Gaussian line model (port of ``kubeviz_ngaussmodel``).

Parameter vector layout (only the parameters being fit are present)::

    a[cpar]                    continuum level (fixed to 0 in SDSS continuum mode)
    a[nparstart]               narrow dv        [km/s]
    a[nparstart+1]             narrow sigma(v)  [km/s]  (intrinsic, instrument added in quadrature)
    a[nparstart+2+i]           narrow flux of line i
    a[bparstart], ...          same for the broad component

Each line ``i`` has rest (observed-frame, redshifted) wavelength ``lambda0``::

    xcen  = lambda0 (1 + dv/c)
    width = lambda0 sigma(v)/c
    sigma_tot^2 = width^2 + (xcen / (2.35482 R(xcen)))^2
    g = flux / sqrt(2 pi sigma_tot^2) exp(-(x-xcen)^2 / (2 sigma_tot^2))
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..constants import CKMS, SIGTOFWHM


@dataclass
class ModelSpec:
    cpar: int
    nparstart: int             # -1 when no narrow component
    bparstart: int             # -1 when no broad component
    nlines: np.ndarray         # observed wavelengths of narrow lines to fit
    blines: np.ndarray         # observed wavelengths of broad lines to fit
    instrres: object           # callable lam -> R (resolution)

    @property
    def npars(self) -> int:
        n = 1
        if self.nparstart >= 0:
            n += 2 + len(self.nlines)
        if self.bparstart >= 0:
            n += 2 + len(self.blines)
        return n


def gaussian_component(x, lambda0, dv, sigv, norm, instrres):
    xcen = lambda0 * (1.0 + dv / CKMS)
    width = lambda0 * (sigv / CKMS)
    R = np.asarray(instrres(xcen), dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        instr_A = np.where(R > 0, xcen / (SIGTOFWHM * R), 0.0)
    totw2 = width ** 2 + instr_A ** 2
    if not np.all(totw2 > 0):
        return np.zeros_like(x)
    return (norm / np.sqrt(2.0 * np.pi * totw2)) * np.exp(-((x - xcen) ** 2 / totw2) / 2.0)


def ngaussmodel(x, a, spec: ModelSpec):
    """Evaluate the model for parameters ``a`` at wavelengths ``x``."""
    x = np.asarray(x, dtype=float)
    fx = np.full(x.shape, a[spec.cpar], dtype=float)
    if spec.nparstart >= 0:
        dv, sigv = a[spec.nparstart], a[spec.nparstart + 1]
        for i, lam0 in enumerate(spec.nlines):
            fx += gaussian_component(x, lam0, dv, sigv, a[spec.nparstart + 2 + i], spec.instrres)
    if spec.bparstart >= 0:
        dv, sigv = a[spec.bparstart], a[spec.bparstart + 1]
        for i, lam0 in enumerate(spec.blines):
            fx += gaussian_component(x, lam0, dv, sigv, a[spec.bparstart + 2 + i], spec.instrres)
    return fx
