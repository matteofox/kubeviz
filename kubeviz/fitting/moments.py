"""Line moments (ports of ``kubeviz_cut_momoutliers`` and ``kubeviz_linefit_moments``).

For each selected line five numbers are produced: integrated flux (0th moment times
the channel width), velocity offset, velocity dispersion (instrumental resolution
subtracted in quadrature), normalised skewness and kurtosis.
"""
from __future__ import annotations

import numpy as np

from ..constants import CKMS, NO_ERRORS, SIGTOFWHM


def cut_momoutliers(good: np.ndarray) -> np.ndarray:
    """Remove isolated pixels (single or pairs) at both ends of an index list so
    that the moment window is a contiguous run. Returns an empty array when fewer
    than three pixels survive."""
    m = np.asarray(good, dtype=int)
    if m.size <= 2:
        return np.zeros(0, dtype=int)
    while m.size > 2:
        if m[1] - 1 != m[0]:
            m = m[1:]
        elif m[2] - 2 != m[0]:
            m = m[2:]
        else:
            break
    while m.size > 2:
        i = m.size - 1
        if m[i] - 1 != m[i - 1]:
            m = m[:i]
        elif m[i] - 2 != m[i - 2]:
            m = m[:i - 1]
        else:
            break
    if m.size < 3:
        return np.zeros(0, dtype=int)
    return m


def linefit_moments(state, ctx, wave, spec, weights, col=None, row=None, boot: int = -1):
    """Compute the moments of every line flagged in ``ctx.mfitpars``.

    Returns ``(pars, sigpars)`` with 5 entries per line; ``sigpars`` is filled with
    :data:`NO_ERRORS` (errors only come from Monte Carlo realisations).
    """
    wave = np.asarray(wave, dtype=float)
    spec = np.asarray(spec, dtype=float)
    weights = np.asarray(weights, dtype=float)

    fixedwin = False
    vcen, vwin = 0.0, 300.0
    if col is not None and row is not None and state.mom_windowmap is not None:
        fixedwin = True
        bmap = None
        if boot >= 0 and state.mom_windowmap_bootstrap is not None and boot < len(state.mom_windowmap_bootstrap):
            bmap = state.mom_windowmap_bootstrap[boot]
        src = bmap if bmap is not None else state.mom_windowmap
        vcen, vwin = float(src[0, row, col]), float(src[1, row, col])

    lines_mom = np.flatnonzero(ctx.mfitpars == 1)
    nlm = lines_mom.size
    pars = np.zeros(5 * nlm)
    sigpars = np.full(5 * nlm, NO_ERRORS)
    sigcont = float(ctx.sigcont) if ctx.sigcont is not None else 0.0

    for j, il in enumerate(lines_mom):
        lambdaz = state.lines[il]
        lambda0 = lambdaz * (1.0 + vcen / CKMS)
        dlambdawin = lambda0 * (vwin / CKMS)
        bis = state.line_bisectors[il] * (1.0 + vcen / CKMS)

        good = np.flatnonzero((wave > bis[0]) & (wave < bis[1]) & (spec > state.mom_thresh * sigcont))
        good = cut_momoutliers(good)
        if good.size > 0:
            sw = spec[good] * weights[good]
            norm = np.sum(sw)
            with np.errstate(all="ignore"):
                mom1 = np.sum(sw * wave[good]) / norm
                mom2 = np.sqrt(np.sum(sw * (wave[good] - mom1) ** 2) / norm)
                mom3 = np.sum(sw * ((wave[good] - mom1) / mom2) ** 3) / norm
                mom4 = np.sum(sw * (((wave[good] - mom1) / mom2) ** 4 - 3)) / norm
                pars[j * 5 + 1] = CKMS * (mom1 / lambdaz - 1.0)
                R = float(state.getinstrres(mom1))
                instr = (mom1 / (SIGTOFWHM * R)) ** 2 if R > 0 else 0.0
                val = mom2 ** 2 - instr
                pars[j * 5 + 2] = (CKMS / lambdaz) * np.sqrt(max(val if np.isfinite(val) else 0.0, 0.0))
                pars[j * 5 + 3] = mom3
                pars[j * 5 + 4] = mom4

        okk_max = np.flatnonzero((wave > bis[0]) & (wave < bis[1]))
        okk_min = np.flatnonzero((wave > lambda0 - dlambdawin) & (wave < lambda0 + dlambdawin))
        if fixedwin:
            okk_final = okk_min
        else:
            okk_diff = np.setdiff1d(okk_max, okk_min)
            below = okk_diff[wave[okk_diff] < lambda0]
            above = okk_diff[wave[okk_diff] > lambda0]
            if below.size > 0:
                lim_b = below[-1]
                while spec[lim_b] > 0 and lim_b >= below[0] and lim_b >= 1:
                    lim_b -= 1
                lim_b += 1
            else:
                lim_b = okk_min.min() if okk_min.size else 0
            if above.size > 0:
                lim_a = above[0]
                while spec[lim_a] > 0 and lim_a <= above[-1] and lim_a <= spec.size - 2:
                    lim_a += 1
                lim_a -= 1
            else:
                lim_a = okk_min.max() if okk_min.size else -1
            okk_final = np.arange(lim_b, lim_a + 1) if lim_a >= lim_b else np.zeros(0, dtype=int)

        mom0 = np.nansum(spec[okk_final]) if okk_final.size > 0 else -99.0
        pars[j * 5 + 0] = mom0 * state.dlambda
    return pars, sigpars
