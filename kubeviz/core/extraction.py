"""Spectrum extraction and collapsed images (ports of ``kubeviz_medianspec``,
``kubeviz_getspec``, ``kubeviz_medsum_image_update`` and ``kubeviz_montecarlonoise``).

Deviation from IDL: sums over a spaxel mask ignore NaN spaxels (IDL ``total`` without
``/nan`` would return NaN for the whole wavelength channel as soon as one masked spaxel
is NaN). A channel is NaN only when every masked spaxel is NaN.
"""
from __future__ import annotations

import numpy as np

from .. import utils
from ..constants import (IMG_MED1, IMG_MED1_MINUS_MED2, IMG_MED2, IMG_MED2_MINUS_MED1,
                         IMG_MEDSUB1, IMG_MEDSUB2, IMG_SLICE, IMG_SUM1, IMG_SUM2, IMG_WAVG1,
                         IMG_WAVG2, IMG_WMED1, IMG_WMED2, SPEC_MEDIAN, SPEC_MEDSUB,
                         SPEC_OPTIMAL, SPEC_SLICE, SPEC_SUM, SPEC_WAVG)

# ----------------------------------------------------------------------------- helpers
def _nansum_or_nan(a, axis):
    """nansum, but NaN where every element along ``axis`` is NaN."""
    with np.errstate(all="ignore"):
        s = np.nansum(a, axis=axis)
        allnan = np.all(~np.isfinite(a), axis=axis)
    s = np.asarray(s, dtype=float)
    s[allnan] = np.nan
    return s


def weighted_median_axis0(values: np.ndarray, weights: np.ndarray) -> np.ndarray:
    """Vectorised ``kubeviz_weighted_median`` along axis 0 of ``values``.

    NaN values are treated as +inf with zero weight (ignored).
    """
    v = np.array(values, dtype=float)
    w = np.array(weights, dtype=float)
    bad = ~np.isfinite(v) | ~np.isfinite(w)
    v[bad] = np.inf
    w[bad] = 0.0
    order = np.argsort(v, axis=0, kind="stable")
    vs = np.take_along_axis(v, order, axis=0)
    ws = np.take_along_axis(w, order, axis=0)
    csum = np.cumsum(ws, axis=0)
    tot = csum[-1]
    # first index where cumulative weight > tot/2
    idx = np.sum(csum <= tot[None] / 2.0, axis=0)
    idx = np.clip(idx, 0, v.shape[0] - 1)
    out = np.take_along_axis(vs, idx[None], axis=0)[0]
    out[~np.isfinite(out)] = np.nan
    out[tot <= 0] = np.nan
    return out


# ----------------------------------------------------------------------------- images
def image_for_imgmode(state, imgmode: int | None = None):
    """(image, noise_image) currently implied by ``state.imgmode`` (data cube only)."""
    im = state.imgmode if imgmode is None else imgmode
    d = state.datacube
    n = state.noise
    k = state.wpix
    with np.errstate(all="ignore"):
        if im == IMG_SLICE:
            return d[k], n[k]
        if im in (IMG_SUM1, IMG_MED1, IMG_WAVG1, IMG_WMED1):
            return state.img1, state.nimg1
        if im == IMG_MEDSUB1:
            return d[k] - state.img1, np.sqrt(n[k] ** 2 + state.nimg1 ** 2)
        if im in (IMG_SUM2, IMG_MED2, IMG_WAVG2, IMG_WMED2):
            return state.img2, state.nimg2
        if im == IMG_MEDSUB2:
            return d[k] - state.img2, np.sqrt(n[k] ** 2 + state.nimg2 ** 2)
        if im == IMG_MED2_MINUS_MED1:
            return state.img2 - state.img1, np.sqrt(state.nimg1 ** 2 + state.nimg2 ** 2)
        if im == IMG_MED1_MINUS_MED2:
            return state.img1 - state.img2, np.sqrt(state.nimg1 ** 2 + state.nimg2 ** 2)
    raise ValueError(f"unknown imgmode {im}")


def _imgmode_props(imgmode: int):
    """(operation, range, weighted) for an image mode; operation None for slice."""
    if imgmode == IMG_SLICE:
        return None, 0, False
    op = "sum" if imgmode in (IMG_SUM1, IMG_WAVG1, IMG_SUM2, IMG_WAVG2) else "med"
    rng = 1 if imgmode in (IMG_SUM1, IMG_MED1, IMG_WAVG1, IMG_WMED1, IMG_MEDSUB1,
                           IMG_MED2_MINUS_MED1) else 2
    weighted = imgmode in (IMG_WAVG1, IMG_WMED1, IMG_WAVG2, IMG_WMED2)
    return op, rng, weighted


def medsum_image_update(state, imgmode: int | None = None, lim=None) -> None:
    """Recompute ``img1/nimg1/bimg1`` or ``img2/...`` for the wavelength range implied
    by ``imgmode`` (port of ``kubeviz_medsum_image_update``)."""
    im = state.imgmode if imgmode is None else imgmode
    op, rng, weighted = _imgmode_props(im)
    if op is None:
        return
    if lim is None:
        lo, hi = (state.wavrange1 if rng == 1 else state.wavrange2)
    else:
        lo, hi = lim
    lo, hi = int(lo), int(hi)
    if hi < lo:
        lo, hi = hi, lo
    shape = (state.Nrow, state.Ncol)
    img = np.zeros(shape)
    nimg = np.zeros(shape)
    bimg = np.zeros(shape, dtype=bool)

    if lo != hi:
        d = state.datacube[lo:hi + 1]
        n = state.noise[lo:hi + 1]
        use_mc = state.domontecarlo > 0 and state.useMonteCarlonoise
        with np.errstate(all="ignore"):
            if op == "sum":
                if not weighted:
                    img = np.sum(d, axis=0)
                    nimg = np.sqrt(np.sum(n ** 2, axis=0)) if not use_mc else \
                        montecarlonoise(state, "image", operation="sum", rng=(lo, hi))
                else:
                    img = np.sum(d / n ** 2, axis=0) / np.sum(1.0 / n ** 2, axis=0)
                    nimg = np.sqrt(1.0 / np.sum(n ** -2, axis=0)) if not use_mc else \
                        montecarlonoise(state, "image", operation="w_avg", rng=(lo, hi))
            else:
                if not weighted:
                    img = utils.idl_median(d, axis=0)
                    nimg = np.sqrt(np.sum(n ** 2, axis=0)) / (1 + hi - lo) if not use_mc else \
                        montecarlonoise(state, "image", operation="med", rng=(lo, hi))
                else:
                    img = weighted_median_axis0(d, 1.0 / n ** 2)
                    nimg = np.sqrt(1.0 / np.sum(n ** -2, axis=0)) if not use_mc else \
                        montecarlonoise(state, "image", operation="w_med", rng=(lo, hi))
        bimg = np.all(state.badpixelmask[lo:hi + 1], axis=0)
        img = utils.remove_badvalues(img)
        nimg = utils.remove_badvalues(nimg)

    if rng == 1:
        state.img1, state.nimg1, state.bimg1 = img, nimg, bimg
    else:
        state.img2, state.nimg2, state.bimg2 = img, nimg, bimg


def update_all_range_images(state) -> None:
    """Refresh both range images (used after a wavelength range changes)."""
    medsum_image_update(state, imgmode=state.imgmode)
    op, rng, _ = _imgmode_props(state.imgmode)
    if state.imgmode in (IMG_MED2_MINUS_MED1, IMG_MED1_MINUS_MED2):
        medsum_image_update(state, imgmode=IMG_MED2 if rng == 1 else IMG_MED1)


def spec_reset_range(state) -> None:
    state.wavrange1[:] = 0
    state.wavrange2[:] = 0
    state.img1 = np.zeros((state.Nrow, state.Ncol))
    state.img2 = np.zeros((state.Nrow, state.Ncol))
    state.nimg1 = np.zeros((state.Nrow, state.Ncol))
    state.nimg2 = np.zeros((state.Nrow, state.Ncol))
    state.wavsel = 1


# ----------------------------------------------------------------------------- Monte Carlo noise
def montecarlonoise(state, mode: str, operation: str | None = None, rng=None) -> np.ndarray:
    """Half the 16-84 percentile spread of the Monte Carlo realisations
    (port of ``kubeviz_montecarlonoise``)."""
    mc = state.montecarlocubes
    if mc is None:
        raise RuntimeError("No Monte Carlo cubes available")
    p_hi, p_lo = state.montecarlo_percs[0], state.montecarlo_percs[1]
    with np.errstate(all="ignore"):
        if mode == "cube":
            q = np.nanpercentile(mc, [p_hi, p_lo], axis=0)
            return 0.5 * (q[0] - q[1])
        if mode == "image":
            lo, hi = int(rng[0]), int(rng[1])
            sub = mc[:, lo:hi + 1]                              # (Nmc, k, Nrow, Ncol)
            n = state.noise[lo:hi + 1]
            if operation == "sum":
                temp = np.sum(sub, axis=1)
            elif operation == "med":
                temp = utils.idl_median(sub, axis=1)
            elif operation == "w_avg":
                temp = np.sum(sub / n[None] ** 2, axis=1) / np.sum(1.0 / n ** 2, axis=0)[None]
            elif operation == "w_med":
                temp = np.stack([weighted_median_axis0(sub[i], 1.0 / n ** 2) for i in range(sub.shape[0])])
            else:
                raise ValueError(operation)
            q = np.nanpercentile(temp, [p_hi, p_lo], axis=0)
            return 0.5 * (q[0] - q[1])
        if mode == "spec":
            q = np.nanpercentile(state.medspec_montecarlo, [p_hi, p_lo], axis=0)
            return 0.5 * (q[0] - q[1])
    raise ValueError(mode)


# ----------------------------------------------------------------------------- spectra
def _extract(mode: int, data: np.ndarray, noise: np.ndarray, mask: np.ndarray,
             weights2d: np.ndarray | None, image: np.ndarray | None, compute_noise: bool):
    """Core extraction of one cube over ``mask`` (bool, Nrow x Ncol).

    Returns (spec, noise_spec or None).
    """
    sel = mask.astype(bool)
    d = data[:, sel]                       # (Nw, Nsel)
    n = noise[:, sel] if compute_noise or mode in (2, 3) else None
    nsel = d.shape[1]
    with np.errstate(all="ignore"):
        if mode == 0:                      # SUM
            spec = _nansum_or_nan(d, axis=1)
            nspec = np.sqrt(_nansum_or_nan(n ** 2, axis=1)) if compute_noise else None
        elif mode == 1:                    # MEDIAN (IDL median: upper middle for even N)
            spec = utils.idl_median(d, axis=1)
            nspec = np.sqrt(_nansum_or_nan(n ** 2, axis=1)) / nsel if compute_noise else None
        elif mode == 2:                    # WEIGHTED AVERAGE
            if weights2d is None:
                w = 1.0 / n ** 2
                spec = _nansum_or_nan(w * d, axis=1) / np.nansum(w, axis=1)
                nspec = np.sqrt(1.0 / np.nansum(n ** -2, axis=1)) if compute_noise else None
            else:
                w = np.broadcast_to(weights2d[sel][None, :], d.shape)
                spec = _nansum_or_nan(w * d, axis=1) / np.nansum(np.where(np.isfinite(d), w, np.nan), axis=1)
                nspec = np.sqrt(np.nansum(w * n ** 2, axis=1) / np.nansum(w, axis=1)) if compute_noise else None
        elif mode == 3:                    # OPTIMAL (Robertson 1986 style)
            frac = image[sel] / np.nansum(image[sel])
            scaled_noi = n / frac[None, :]
            scaled_img = d / frac[None, :]
            spec = np.nansum(scaled_img / scaled_noi ** 2, axis=1) / np.nansum(1.0 / scaled_noi ** 2, axis=1)
            nspec = np.sqrt(1.0 / np.nansum(scaled_noi ** -2, axis=1)) if compute_noise else None
            spec = utils.remove_badvalues(spec)
            if nspec is not None:
                nspec = utils.remove_badvalues(nspec)
        else:
            raise ValueError(mode)
    return spec, nspec


def medianspec(state, mode: int | None = None, all: bool = False, redsh: bool = False,
               sum: bool = False, med: bool = False, wavg: bool = False, optimal: bool = False) -> None:
    """Combine the spectra of the selected spaxels into ``state.medspec``/``nmedspec``
    (and ``medspec_montecarlo``). Port of ``kubeviz_medianspec``.

    ``mode``: 0 sum, 1 median, 2 weighted average, 3 optimal. When None the mode is
    derived from ``state.specmode`` (nothing happens for the single-spaxel mode).
    """
    nw, nrow, ncol = state.Nwpix, state.Nrow, state.Ncol
    weights2d = None
    if all:
        mask = np.ones((nrow, ncol), dtype=bool)
    elif redsh:
        mask = np.zeros((nrow, ncol), dtype=bool)
        mask[nrow // 3: (nrow // 3) * 2 + 1, ncol // 3: (ncol // 3) * 2 + 1] = True
    else:
        raw = np.asarray(state.current_mask(), dtype=float)
        ok = raw > 0
        if np.any(ok) and not np.all(raw[ok] == 1):
            weights2d = raw.copy()
        mask = ok

    if sum:
        mode = 0
    elif med:
        mode = 1
    elif wavg:
        mode = 2
    elif optimal:
        mode = 3
    if mode is None:
        table = {SPEC_SUM: 0, SPEC_MEDIAN: 1, SPEC_WAVG: 2, SPEC_MEDSUB: 1, SPEC_OPTIMAL: 3}
        if state.specmode not in table:
            return
        mode = table[state.specmode]

    if not np.any(mask):
        state.medspec = np.zeros(nw)
        state.nmedspec = np.zeros(nw)
        return

    image = None
    if mode == 3:
        image, _ = image_for_imgmode(state)
        with np.errstate(invalid="ignore"):
            if np.nanmin(np.where(mask, image, np.nan)) < 0:
                utils.warn("The current mask includes spaxels with negative values in the selected image.")
                utils.warn("The optimal method is not reliable. Sum method will be used instead.")
                mode = 0
                state.specmode = SPEC_SUM

    use_mc_noise = state.domontecarlo > 0 and state.useMonteCarlonoise
    spec, nspec = _extract(mode, state.datacube, state.noise, mask, weights2d, image,
                           compute_noise=not use_mc_noise)

    mc = state.montecarlocubes
    need_mc_specs = state.domontecarlo > 0 and mc is not None and (mode != 3 or use_mc_noise)
    if need_mc_specs:
        state.medspec_montecarlo = np.stack([
            _extract(mode, mc[i], state.noise, mask, weights2d, image, compute_noise=False)[0]
            for i in range(mc.shape[0])])
        if use_mc_noise:
            nspec = montecarlonoise(state, "spec")
    elif use_mc_noise:
        # no realisations available for the noise: fall back to propagated noise
        _, nspec = _extract(mode, state.datacube, state.noise, mask, weights2d, image, True)

    state.medspec = np.asarray(spec, dtype=float)
    state.nmedspec = np.asarray(nspec, dtype=float)


def getspec(state):
    """Current 1D spectrum and noise according to ``state.specmode``
    (port of ``kubeviz_getspec``)."""
    if state.specmode == SPEC_SLICE:
        return state.datacube[:, state.row, state.col].astype(float), \
            state.noise[:, state.row, state.col].astype(float)
    if state.specmode == SPEC_MEDSUB:
        spec = state.datacube[:, state.row, state.col] - state.medspec
        nspec = np.sqrt(state.noise[:, state.row, state.col] ** 2 + state.nmedspec ** 2)
        return spec.astype(float), nspec.astype(float)
    return np.asarray(state.medspec, dtype=float), np.asarray(state.nmedspec, dtype=float)


def spectrum_for_fit(state, col: int, row: int):
    """Spectrum, noise and Monte Carlo spectra used by the fitter for spaxel
    ``(col,row)`` or the current mask, following ``kubeviz_linefit_dofit``.

    Returns ``(spec, noise, spec_bootstrap)`` where ``spec_bootstrap`` is
    ``(Nmontecarlo, Nwpix)`` or ``None``.
    """
    nw = state.Nwpix
    mc = state.montecarlocubes if state.domontecarlo > 0 else None
    boot = None
    if state.specmode in (SPEC_SLICE, SPEC_MEDSUB):
        if state.specmode == SPEC_MEDSUB:
            utils.info("Median subtracted spectrum cannot be fit: Fitting whole spectrum")
        spec = state.datacube[:, row, col].astype(float)
        noise = state.noise[:, row, col].astype(float)
        if mc is not None and np.nansum(spec) != 0:
            boot = mc[:, :, row, col].astype(float)
    else:
        medianspec(state)
        spec = np.asarray(state.medspec, dtype=float)
        noise = np.asarray(state.nmedspec, dtype=float)
        if mc is not None and np.nansum(spec) != 0 and state.medspec_montecarlo is not None:
            boot = np.asarray(state.medspec_montecarlo, dtype=float)
    if state.domontecarlo > 0 and boot is None:
        boot = np.zeros((state.Nmontecarlo, nw))
    return spec, noise, boot
