"""Continuum estimation (ports of ``kubeviz_definecontinuumfitregion``,
``kubeviz_fitcontinuum`` and ``kubeviz_getcontatlambda``)."""
from __future__ import annotations

import numpy as np

from .. import utils


def getcontatlambda(contpars, lam):
    """Polynomial continuum evaluated at ``lam`` (coefficients in increasing order)."""
    contpars = np.asarray(contpars, dtype=float)
    lam = np.asarray(lam, dtype=float)
    return np.polynomial.polynomial.polyval(lam, contpars)


def polyfit_cov(x, y, order: int, weights=None):
    """Polynomial least-squares fit ``y = sum_i a_i x^i`` returning ``(coeffs, sigma)``.

    ``weights`` are inverse variances (IDL ``svdfit`` WEIGHTS). Without weights the
    parameter errors are scaled so that chi-square per degree of freedom is one, as
    ``svdfit`` does when no measurement errors are given.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    A = np.vander(x, order + 1, increasing=True)
    if weights is None:
        sw = np.ones_like(y)
    else:
        sw = np.sqrt(np.asarray(weights, dtype=float))
    Aw = A * sw[:, None]
    yw = y * sw
    coeff, *_ = np.linalg.lstsq(Aw, yw, rcond=None)
    resid = yw - Aw @ coeff
    chisq = float(np.sum(resid ** 2))
    dof = max(x.size - (order + 1), 1)
    try:
        cov = np.linalg.pinv(Aw.T @ Aw)
        sigma = np.sqrt(np.clip(np.diag(cov), 0, None))
    except Exception:  # pragma: no cover
        sigma = np.zeros(order + 1)
    if weights is None:
        sigma = sigma * np.sqrt(chisq / dof)
    return coeff, sigma


def definecontinuumfitregion(state, wave, wave_cen):
    """Indices (into ``wave``) usable for the continuum fit around ``wave_cen``.

    Two windows offset by ``[minoff, maxoff]`` on both sides of the main line are
    used, excluding: 10 A at both ends of the spectrum, the hardcoded bad regions and
    ``+/- maxwoff`` around every line in the spectrum (also lines not being fit).
    """
    wave = np.asarray(wave, dtype=float)
    minoff, maxoff = state.continuumfit_minoff, state.continuumfit_maxoff
    nocuts = ((wave >= wave_cen - maxoff) & (wave <= wave_cen - minoff)) | \
             ((wave >= wave_cen + minoff) & (wave <= wave_cen + maxoff))
    suitable = (wave >= state.wave[0] + 10) & (wave <= state.wave[-1] - 10) & (~utils.excludewave(wave))
    for line in state.lines_all[: state.Nlines_all]:
        suitable &= ~((wave >= line - state.maxwoffb) & (wave <= line + state.maxwoffr))
    region = np.flatnonzero(nocuts & suitable)
    return region, region.size


def fitcontinuum(state, contfitregion, spec, noise, wave, wave_cen):
    """Fit a polynomial of order ``state.continuumfit_order`` to the continuum region.

    Returns ``(contpars, sigcontpars, sigcont)``. Two methods, as in IDL:

    * ``minperc == 0 and maxperc == 100``: inverse-variance weighted fit after
      rejecting noisy channels and spikes;
    * otherwise the "SDSS" method: reject the 20% noisiest channels on each side, keep
      the channels between the (rescaled) percentile limits, unweighted fit.
    """
    order = int(state.continuumfit_order)
    spec_c = np.asarray(spec, dtype=float)[contfitregion]
    wave_c = np.asarray(wave, dtype=float)[contfitregion]
    noise_c = np.asarray(noise, dtype=float)[contfitregion]
    zero = (np.zeros(order + 1), np.zeros(order + 1), 0.0)

    if state.continuumfit_minperc == 0 and state.continuumfit_maxperc == 100:
        with np.errstate(divide="ignore", invalid="ignore"):
            weights = 1.0 / noise_c ** 2
        if noise_c.size == 0:
            return zero
        n90, n10 = utils.percentile(noise_c, 90.0), utils.percentile(noise_c, 10.0)
        s90, s10 = utils.percentile(spec_c, 90.0), utils.percentile(spec_c, 10.0)
        okval = np.isfinite(weights) & np.isfinite(spec_c) & \
            (np.abs(noise_c - utils.idl_median(noise_c)) <= n90 - n10) & \
            (np.abs(spec_c - utils.idl_median(spec_c)) <= s90 - s10)
        if not np.any(okval):
            return zero
        contpars, sigcontpars = polyfit_cov(wave_c[okval], spec_c[okval], order, weights=weights[okval])
        polyspec = getcontatlambda(contpars, wave_c[okval])
        r = spec_c[okval] - polyspec
        sigcont = 0.5 * (utils.percentile(r, 84.0) + abs(utils.percentile(r, 16.0)))
        return contpars, sigcontpars, float(sigcont)

    noisecut_perc = 20.0
    new_range = (state.continuumfit_maxperc - state.continuumfit_minperc) / (0.01 * (100 - noisecut_perc))
    mid_perc = 0.5 * (state.continuumfit_minperc + state.continuumfit_maxperc)
    min_perc = mid_perc - 0.5 * new_range
    max_perc = mid_perc + 0.5 * new_range

    def side(sel):
        if not np.any(sel):
            return np.zeros(0), np.zeros(0)
        s, w, n = spec_c[sel], wave_c[sel], noise_c[sel]
        fin = np.isfinite(s) & np.isfinite(n)
        s, w, n = s[fin], w[fin], n[fin]
        if s.size == 0:
            return np.zeros(0), np.zeros(0)
        ncut = utils.percentile(n, 100 - noisecut_perc)
        ok = n < ncut
        if not np.any(ok):
            return np.zeros(0), np.zeros(0)
        s, w = s[ok], w[ok]
        lo, hi = utils.percentile(s, [min_perc, max_perc])
        ok = (s >= lo) & (s <= hi)
        return s[ok], w[ok]

    spec_b, wave_b = side(wave_c < wave_cen)
    spec_r, wave_r = side(wave_c > wave_cen)
    spec_fit = np.concatenate([spec_b, spec_r])
    wave_fit = np.concatenate([wave_b, wave_r])
    if spec_fit.size == 0:
        return zero
    contpars, sigcontpars = polyfit_cov(wave_fit, spec_fit, order)
    polyspec = getcontatlambda(contpars, wave_fit)
    r = spec_fit - polyspec
    sigcont = 0.5 * (utils.percentile(r, 84.0) + abs(utils.percentile(r, 16.0)))
    return contpars, sigcontpars, float(sigcont)
