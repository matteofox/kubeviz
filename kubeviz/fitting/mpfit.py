"""A small mpfit-compatible layer over :func:`scipy.optimize.least_squares`.

The IDL code drives ``mpcurvefit`` with a ``parinfo`` array (``fixed``, ``limited``,
``limits``, ``tied``). :class:`ParInfo` mirrors that structure and
:func:`mpcurvefit` returns the same quantities (best-fit parameters, formal 1-sigma
errors from the covariance matrix, chi-square, degrees of freedom, status, error
message, iteration count and the model evaluated at the solution).

Differences from mpfit worth knowing:

* the minimiser is the trust-region-reflective algorithm with bounds instead of the
  bounded Levenberg-Marquardt of mpfit; both converge to the same local minimum for
  the well-posed problems kubeviz sets up;
* a start value outside the user limits is clipped into the allowed range instead of
  raising an error;
* tied parameters are given as ``"P[i]/k"`` strings (the only form used by kubeviz)
  or as callables ``f(p) -> value``.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field

import numpy as np
from scipy.optimize import least_squares


@dataclass
class ParInfo:
    fixed: bool = False
    limited: list = field(default_factory=lambda: [False, False])
    limits: list = field(default_factory=lambda: [-9.0e9, 9.0e9])
    tied: object = None      # None, "P[i]/k" or callable(p_full)->float

    @property
    def is_free(self) -> bool:
        return (not self.fixed) and (self.tied is None)


@dataclass
class MpfitResult:
    params: np.ndarray
    perror: np.ndarray
    chisq: float
    dof: int
    status: int
    errmsg: str
    niter: int
    yfit: np.ndarray | None
    covar: np.ndarray | None = None

    @property
    def ok(self) -> bool:
        return self.status > 0 and self.errmsg == ""


_TIED_RE = re.compile(r"^\s*P\[(\d+)\]\s*(?:([*/])\s*([-+0-9.eE]+))?\s*$")


def _compile_tied(expr):
    if expr is None:
        return None
    if callable(expr):
        return expr
    m = _TIED_RE.match(str(expr))
    if not m:
        raise ValueError(f"Unsupported tied expression: {expr!r}")
    i = int(m.group(1))
    op, k = m.group(2), m.group(3)
    if op is None:
        return lambda p, i=i: p[i]
    k = float(k)
    if op == "/":
        return lambda p, i=i, k=k: p[i] / k
    return lambda p, i=i, k=k: p[i] * k


def mpcurvefit(func, x, y, weights, p0, parinfo=None, itmax: int = 200, quiet: bool = True,
               **kwargs) -> MpfitResult:
    """Weighted least-squares fit of ``func(x, p, **kwargs)`` to ``y``.

    ``weights`` are inverse variances (mpcurvefit convention); pixels with zero weight
    are ignored. Returns an :class:`MpfitResult`.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    p0 = np.array(p0, dtype=float)
    npar = p0.size
    if weights is None:
        weights = np.ones_like(y)
    weights = np.asarray(weights, dtype=float)
    if parinfo is None:
        parinfo = [ParInfo() for _ in range(npar)]
    tied = [_compile_tied(pi.tied) for pi in parinfo]
    free = np.array([pi.is_free for pi in parinfo], dtype=bool)
    ifree = np.flatnonzero(free)
    nfree = ifree.size

    good = np.isfinite(y) & np.isfinite(weights) & (weights > 0) & np.isfinite(x)
    xg, yg, sw = x[good], y[good], np.sqrt(weights[good])
    npts = xg.size

    lower = np.array([pi.limits[0] if pi.limited[0] else -np.inf for pi in parinfo])[ifree]
    upper = np.array([pi.limits[1] if pi.limited[1] else np.inf for pi in parinfo])[ifree]

    def expand(pf):
        p = p0.copy()
        p[ifree] = pf
        for i, t in enumerate(tied):
            if t is not None:
                p[i] = t(p)
        return p

    def model(p):
        return func(xg, p, **kwargs)

    if npts == 0 or npts < nfree:
        return MpfitResult(np.full(npar, -999.0), np.full(npar, -999.0), 0.0, 0, 0,
                           "mpfit: number of parameters exceeds number of data points", 0, None)

    if nfree == 0:
        p = expand(np.zeros(0))
        yfit = model(p)
        chisq = float(np.sum((sw * (yg - yfit)) ** 2))
        full = np.full(x.shape, np.nan)
        full[good] = yfit
        return MpfitResult(p, np.zeros(npar), chisq, npts, 1, "", 0, full)

    start = np.clip(p0[ifree], lower, upper)
    # avoid starting exactly on a bound where the reflective method has no room to move
    span = np.where(np.isfinite(upper - lower), upper - lower, 1.0)
    eps = 1e-8 * np.maximum(np.abs(start), 1.0)
    at_lo = np.isfinite(lower) & (start <= lower)
    at_hi = np.isfinite(upper) & (start >= upper)
    start = np.where(at_lo, lower + np.minimum(eps, 0.5 * span), start)
    start = np.where(at_hi, upper - np.minimum(eps, 0.5 * span), start)

    def resid(pf):
        r = sw * (yg - model(expand(pf)))
        return np.where(np.isfinite(r), r, 1e30)

    errmsg = ""
    try:
        res = least_squares(resid, start, bounds=(lower, upper), method="trf",
                            x_scale="jac", max_nfev=int(itmax) * (nfree + 1) * 2,
                            ftol=1e-10, xtol=1e-10, gtol=1e-10)
    except Exception as exc:  # pragma: no cover - numerical failure path
        return MpfitResult(np.full(npar, -999.0), np.full(npar, -999.0), 0.0, 0, 0,
                           f"mpfit failed: {exc}", 0, None)

    p = expand(res.x)
    yfit_g = model(p)
    chisq = float(np.sum((sw * (yg - yfit_g)) ** 2))
    dof = max(npts - nfree, 1)
    perror = np.zeros(npar)
    covar = None
    # mpfit convention: parameters pegged at a limit are excluded from the covariance
    # and get a zero error (their derivative is meaningless at the bound)
    tol = 1e-9 * np.maximum(np.abs(res.x), 1.0)
    pegged = (np.isfinite(lower) & (res.x <= lower + tol)) | (np.isfinite(upper) & (res.x >= upper - tol))
    unpegged = np.flatnonzero(~pegged)
    try:
        J = res.jac[:, unpegged]
        JTJ = J.T @ J
        cov_free = np.linalg.pinv(JTJ)
        perror[ifree[unpegged]] = np.sqrt(np.clip(np.diag(cov_free), 0, None))
        covar = np.zeros((npar, npar))
        covar[np.ix_(ifree[unpegged], ifree[unpegged])] = cov_free
    except Exception:  # pragma: no cover
        perror[ifree] = 0.0
    if not np.all(np.isfinite(p)):
        errmsg = "mpfit: non-finite parameters"
    status = int(res.status) if res.status > 0 else 0
    if res.status <= 0:
        errmsg = errmsg or f"mpfit: minimisation did not converge ({res.message})"
    yfit = np.full(x.shape, np.nan)
    yfit[good] = yfit_g
    return MpfitResult(p, perror, chisq, dof, status, errmsg, int(res.nfev), yfit, covar)


# ----------------------------------------------------------------------------- peak fitting
def _gauss_terms(x, p, nterms):
    g = p[0] * np.exp(-0.5 * ((x - p[1]) / p[2]) ** 2)
    if nterms >= 4:
        g = g + p[3]
    if nterms >= 5:
        g = g + p[4] * x
    return g


def mpfitpeak(x, y, nterms: int = 3, positive: bool = True, estimates=None):
    """Gaussian peak fit in the spirit of IDL ``mpfitpeak``.

    Returns ``(yfit, params, result)`` with ``params = [height, centre, sigma, (offset,
    (slope))]``.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    ok = np.isfinite(x) & np.isfinite(y)
    xs, ys = x[ok], y[ok]
    if xs.size < nterms:
        return np.full_like(x, np.nan), np.full(nterms, np.nan), None
    if estimates is None:
        order = np.argsort(xs)
        xs_s, ys_s = xs[order], ys[order]
        base = np.min(ys_s) if nterms >= 4 else 0.0
        ymax = np.max(ys_s)
        cen = xs_s[np.argmax(ys_s)]
        dx = np.gradient(xs_s) if xs_s.size > 1 else np.ones(1)
        area = np.sum(dx * np.clip(ys_s - base, 0, None))
        height = ymax - base
        width = area / (height * np.sqrt(2 * np.pi)) if height > 0 else (xs_s[-1] - xs_s[0]) / 4
        if not np.isfinite(width) or width <= 0:
            width = max(np.median(dx), 1e-6)
        estimates = [height, cen, width, base, 0.0][:nterms]
    p0 = np.array(estimates, dtype=float)
    lo = np.full(nterms, -np.inf)
    hi = np.full(nterms, np.inf)
    lo[2] = 1e-12 * max(abs(p0[2]), 1.0)
    if positive:
        lo[0] = 0.0
    p0 = np.clip(p0, lo, hi)

    def resid(p):
        return ys - _gauss_terms(xs, p, nterms)

    try:
        res = least_squares(resid, p0, bounds=(lo, hi), method="trf", x_scale="jac")
    except Exception:
        return np.full_like(x, np.nan), np.full(nterms, np.nan), None
    yfit = _gauss_terms(x, res.x, nterms)
    return yfit, res.x, res
