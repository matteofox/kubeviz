"""Result containers (port of ``kubeviz_create_rescubes`` and the
``kubeviz_linefit_pointerswitch`` bookkeeping).

One :class:`ResultSet` holds the narrow, broad, continuum and moment results with
their error percentiles for one ``(linefit_mode, error_method)`` pair. The IDL code
kept 5 x 2 x 8 named pointers and swapped them in and out; here a dictionary of
ResultSets keyed by ``(mode, errmethod)`` replaces the pointer switch.

Plane layout (last axis of ``n``/``b``): ``0`` flag, ``1`` dv [km/s], ``2`` sig(v) [km/s],
``3..2+Nlines`` line fluxes.  Continuum ``c``: ``0`` flag, ``1..Nlines`` continuum at
each line. Moments ``m``: 6 planes per line: flux, dv, sigv, skewness, kurtosis, flag.
Errors carry an extra trailing axis of length :data:`N_MONTECARLO_PERCS`; planes 0 and
1 of that axis are the +1 sigma and -1 sigma offsets (for noise-cube errors: +err, -err).

Spatial layout: spaxel results are ``(Nrow, Ncol, ...)`` so that ``res.n[row, col]`` is
the parameter vector of one spaxel; mask results are ``(Nmask, ...)``.
"""
from __future__ import annotations

import numpy as np

from .constants import (FLAG_MANUAL, MODE_MASK, MODE_SPAXEL, N_MONTECARLO_PERCS,
                        NOT_FIT)


class ResultSet:
    """Fit results for one (linefit_mode, errmethod) combination."""

    def __init__(self, mode: int, nlines: int, ncol: int = 0, nrow: int = 0, nmask: int = 1,
                 npercs: int = N_MONTECARLO_PERCS):
        self.mode = mode
        self.nlines = int(nlines)
        self.npercs = int(npercs)
        if mode == MODE_SPAXEL:
            lead = (int(nrow), int(ncol))
        else:
            lead = (max(int(nmask), 1),)
        self._lead = lead
        nl = max(self.nlines, 0)
        nm = max(6 * nl, 1)
        self.n = np.zeros(lead + (3 + nl,), dtype=np.float64)
        self.b = np.zeros(lead + (3 + nl,), dtype=np.float64)
        self.c = np.zeros(lead + (1 + nl,), dtype=np.float64)
        self.m = np.zeros(lead + (nm,), dtype=np.float64)
        self.nerr = np.full(lead + (3 + nl, self.npercs), NOT_FIT, dtype=np.float64)
        self.berr = np.full(lead + (3 + nl, self.npercs), NOT_FIT, dtype=np.float64)
        self.cerr = np.full(lead + (1 + nl, self.npercs), NOT_FIT, dtype=np.float64)
        self.merr = np.full(lead + (nm, self.npercs), NOT_FIT, dtype=np.float64)
        self.chisq = np.zeros(lead, dtype=np.float64)

    # ------------------------------------------------------------------ geometry
    @property
    def nmask(self) -> int:
        return self._lead[0] if self.mode == MODE_MASK else 0

    def ensure_nmask(self, nmask: int) -> None:
        """Grow the mask axis to hold at least ``nmask`` masks (IDL used a fixed 50000)."""
        if self.mode != MODE_MASK or nmask <= self._lead[0]:
            return
        extra = nmask - self._lead[0]

        def grow(arr, fill):
            pad = np.full((extra,) + arr.shape[1:], fill, dtype=arr.dtype)
            return np.concatenate([arr, pad], axis=0)

        self.n, self.b, self.c, self.m = (grow(self.n, 0.0), grow(self.b, 0.0),
                                          grow(self.c, 0.0), grow(self.m, 0.0))
        self.nerr, self.berr = grow(self.nerr, NOT_FIT), grow(self.berr, NOT_FIT)
        self.cerr, self.merr = grow(self.cerr, NOT_FIT), grow(self.merr, NOT_FIT)
        self.chisq = grow(self.chisq, 0.0)
        self._lead = (nmask,)

    def delete_mask(self, imask0: int) -> None:
        """Remove mask index ``imask0`` (0-based) from all containers."""
        if self.mode != MODE_MASK or self._lead[0] <= 1:
            return
        for name in ("n", "b", "c", "m", "nerr", "berr", "cerr", "merr", "chisq"):
            setattr(self, name, np.delete(getattr(self, name), imask0, axis=0))
        self._lead = (self._lead[0] - 1,)

    # ------------------------------------------------------------------ indexing
    def index(self, col=None, row=None, imask=None):
        """Leading index tuple for a spaxel (``col,row``) or a 1-based mask ``imask``."""
        if self.mode == MODE_SPAXEL:
            return (int(row), int(col))
        return (int(imask) - 1,)

    # ------------------------------------------------------------------ reset
    def reset_element(self, idx, linetypes=("N", "B", "C"), pars=None) -> None:
        """Port of ``kubeviz_linefit_reset`` for one spaxel/mask: set the selected
        parameters to NOT_FIT (flag plane to FLAG_MANUAL)."""
        if pars is None:
            pars = range(3 + self.nlines)
        for lt in linetypes:
            for par in pars:
                nullval = NOT_FIT if par > 0 else float(FLAG_MANUAL)
                if lt == "N":
                    self.n[idx + (par,)] = nullval
                    self.nerr[idx + (par,)] = nullval
                elif lt == "B":
                    self.b[idx + (par,)] = nullval
                    self.berr[idx + (par,)] = nullval
                elif lt == "C":
                    cpar = par - 2
                    if cpar > 0:
                        self.c[idx + (cpar,)] = nullval
                        self.cerr[idx + (cpar,)] = nullval

    def reset_region(self, fitallrange=None, linetypes=("N", "B", "C"), pars=None) -> None:
        """Port of ``kubeviz_linefit_resetall``: reset a spaxel range (or the whole
        container for mask mode)."""
        if pars is None:
            pars = range(3 + self.nlines)
        if self.mode == MODE_SPAXEL and fitallrange is not None:
            c0, c1, r0, r1 = [int(v) for v in fitallrange]
            sl = (slice(r0, r1 + 1), slice(c0, c1 + 1))
        else:
            sl = (slice(None),) * len(self._lead)
        for lt in linetypes:
            for par in pars:
                nullval = NOT_FIT if par > 0 else float(FLAG_MANUAL)
                if lt == "N":
                    self.n[sl + (par,)] = nullval
                    self.nerr[sl + (par,)] = nullval
                elif lt == "B":
                    self.b[sl + (par,)] = nullval
                    self.berr[sl + (par,)] = nullval
                elif lt == "C":
                    cpar = par - 2
                    if cpar > 0:
                        self.c[sl + (cpar,)] = nullval
                        self.cerr[sl + (cpar,)] = nullval

    def copy(self) -> "ResultSet":
        new = ResultSet.__new__(ResultSet)
        new.__dict__.update({k: (v.copy() if isinstance(v, np.ndarray) else v)
                             for k, v in self.__dict__.items()})
        return new


def moment_reshape(m: np.ndarray, main_ind: int, nlines: int, error: bool = False) -> np.ndarray:
    """Port of ``kubeviz_moment_reshape``: view the moment cube like a line cube
    (flag, dv, sigv, fluxes) so that the same display code can be used."""
    lead = m.shape[:-2] if error else m.shape[:-1]
    if error:
        out = np.zeros(lead + (3 + nlines, m.shape[-1]), dtype=m.dtype)
        out[..., 0, :] = m[..., main_ind * 6 + 5, :]
        out[..., 1, :] = m[..., main_ind * 6 + 1, :]
        out[..., 2, :] = m[..., main_ind * 6 + 2, :]
        for line in range(nlines):
            out[..., line + 3, :] = m[..., line * 6, :]
    else:
        out = np.zeros(lead + (3 + nlines,), dtype=m.dtype)
        out[..., 0] = m[..., main_ind * 6 + 5]
        out[..., 1] = m[..., main_ind * 6 + 1]
        out[..., 2] = m[..., main_ind * 6 + 2]
        for line in range(nlines):
            out[..., line + 3] = m[..., line * 6]
    return out
