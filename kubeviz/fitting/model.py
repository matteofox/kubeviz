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
    """One Gaussian line (kept for callers that evaluate a single component)."""
    return _components(x, np.atleast_1d(np.asarray(lambda0, dtype=float)), dv, sigv,
                       np.atleast_1d(np.asarray(norm, dtype=float)), instrres)


def _components(x, lam0, dv, sigv, norms, instrres):
    """Sum of the Gaussians of one kinematic component, all lines evaluated at once.

    ``lam0`` and ``norms`` are arrays over the lines of the set; the instrumental
    resolution is evaluated for all line centres in a single call.
    """
    xcen = lam0 * (1.0 + dv / CKMS)
    width = lam0 * (sigv / CKMS)
    R = np.asarray(instrres(xcen), dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        instr_A = np.where(R > 0, xcen / (SIGTOFWHM * R), 0.0)
    totw2 = width ** 2 + instr_A ** 2
    ok = totw2 > 0                                # a line with zero total width contributes nothing
    if not np.any(ok):
        return np.zeros(x.shape)
    xcen, totw2, norms = xcen[ok], totw2[ok], norms[ok]
    amp = norms / np.sqrt(2.0 * np.pi * totw2)
    arg = (x[None, :] - xcen[:, None]) ** 2 / (2.0 * totw2[:, None])
    return amp @ np.exp(-arg)


def ngaussmodel(x, a, spec: ModelSpec):
    """Evaluate the model for parameters ``a`` at wavelengths ``x``."""
    x = np.asarray(x, dtype=float)
    a = np.asarray(a, dtype=float)
    fx = np.full(x.shape, a[spec.cpar], dtype=float)
    if spec.nparstart >= 0 and len(spec.nlines):
        s = spec.nparstart
        fx += _components(x, np.asarray(spec.nlines, dtype=float), a[s], a[s + 1],
                          a[s + 2: s + 2 + len(spec.nlines)], spec.instrres)
    if spec.bparstart >= 0 and len(spec.blines):
        s = spec.bparstart
        fx += _components(x, np.asarray(spec.blines, dtype=float), a[s], a[s + 1],
                          a[s + 2: s + 2 + len(spec.blines)], spec.instrres)
    return fx
