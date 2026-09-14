"""Spatial and spectral median smoothing (port of ``kubeviz_smooth``).

IDL semantics reproduced:

* the data cube is median-filtered over a ``smooth x smooth x specsmooth`` box; the
  noise is combined as ``sqrt(sum(noise^2)) / N_valid`` (error of the mean);
* NaN values are ignored (IDL ``median`` ignores NaN, ``total(/nan)``);
* an even kernel is not centred: the output pixel ``o`` averages the input pixels
  ``[o - s/2, o + s/2 - 1]`` and the first row/column/plane becomes NaN;
* at the edges the box is truncated.

The IDL implementation looped over every pixel; here the filter is vectorised with
sliding windows over wavelength chunks, which is orders of magnitude faster.
"""
from __future__ import annotations

import numpy as np
from numpy.lib.stride_tricks import sliding_window_view

from .. import utils
from ..constants import ERR_BOOTSTRAP, ERR_MC1


def _pad_for_kernel(s: int) -> tuple[int, int]:
    """(pad_before, pad_after) so that window ``o`` covers the IDL box."""
    s = int(s)
    if s <= 1:
        return 0, 0
    if s % 2 == 1:
        h = (s - 1) // 2
        return h, h
    return s // 2, s // 2 - 1


def _chunk_size(nrow: int, ncol: int, s: int, ss: int, budget_bytes: float = 2.5e8) -> int:
    per_plane = max(nrow * ncol * s * s * ss * 8, 1)
    return max(1, min(64, int(budget_bytes // per_plane)))


def smooth_cube(data: np.ndarray, noise: np.ndarray | None, smooth: int, specsmooth: int,
                progress=None, mode: str = "data"):
    """Median-smooth ``data`` (and combine ``noise``) over a box.

    Parameters
    ----------
    data, noise : (Nw, Nrow, Ncol) arrays. ``noise`` may be None.
    smooth, specsmooth : kernel sizes in spatial / spectral pixels (>= 1).
    mode : ``"data"`` returns (smoothed_data, smoothed_noise); ``"median"`` returns
        only the median-smoothed data (used for Monte Carlo / bootstrap cubes).
    """
    smooth = max(int(abs(smooth)), 1)
    specsmooth = max(int(abs(specsmooth)), 1)
    data = np.asarray(data, dtype=np.float32 if data.dtype == np.float32 else float)
    if smooth == 1 and specsmooth == 1:
        if mode == "median":
            return data.copy()
        return data.copy(), (None if noise is None else np.asarray(noise).copy())

    nw, nrow, ncol = data.shape
    plo, phi = _pad_for_kernel(smooth)
    zlo, zhi = _pad_for_kernel(specsmooth)
    even_spat = (smooth % 2 == 0)
    even_spec = (specsmooth % 2 == 0)

    out = np.full(data.shape, np.nan, dtype=data.dtype)
    out_noise = None if (noise is None or mode == "median") else np.full(data.shape, np.nan, dtype=float)

    # pad spatially once (NaN padding == truncated box because NaN are ignored)
    dpad = np.pad(data, ((zlo, zhi), (plo, phi), (plo, phi)), constant_values=np.nan)
    npad = None
    if out_noise is not None:
        npad = np.pad(np.asarray(noise, dtype=float), ((zlo, zhi), (plo, phi), (plo, phi)),
                      constant_values=np.nan)

    chunk = _chunk_size(nrow, ncol, smooth, specsmooth)
    prog = progress
    for k0 in range(0, nw, chunk):
        k1 = min(nw, k0 + chunk)
        # windows for output planes k0..k1-1 need padded planes k0 .. k1-1 + zlo + zhi
        sub = dpad[k0: k1 + zlo + zhi]
        win = sliding_window_view(sub, (specsmooth, smooth, smooth))   # (k, y, x, ss, s, s)
        win = win.reshape(win.shape[:3] + (-1,))
        with np.errstate(all="ignore"):
            out[k0:k1] = np.nanmedian(win, axis=-1)
        if out_noise is not None:
            nsub = npad[k0: k1 + zlo + zhi]
            nwin = sliding_window_view(nsub, (specsmooth, smooth, smooth))
            nwin = nwin.reshape(nwin.shape[:3] + (-1,))
            with np.errstate(all="ignore"):
                nvalid = np.sum(np.isfinite(nwin), axis=-1)
                out_noise[k0:k1] = np.sqrt(np.nansum(nwin ** 2, axis=-1)) / nvalid
        if prog is not None:
            prog.step(k1 - k0)

    # IDL embeds the reduced cube leaving NaN in the first row/col/plane for even kernels
    if even_spat:
        out[:, 0, :] = np.nan
        out[:, :, 0] = np.nan
        if out_noise is not None:
            out_noise[:, 0, :] = np.nan
            out_noise[:, :, 0] = np.nan
    if even_spec:
        out[0] = np.nan
        if out_noise is not None:
            out_noise[0] = np.nan

    if mode == "median":
        return out
    return out, out_noise


def smooth_state_datacube(state) -> None:
    """Port of ``kubeviz_smooth, /datacube``: fill ``state.datacube``/``noisecube``,
    the bad pixel mask and image, and point ``state.noise`` to the smoothed noise."""
    if state.smooth < 1:
        state.smooth = 1
    if state.specsmooth < 1:
        state.specsmooth = 1
    if state.smooth == 1 and state.specsmooth == 1:
        state.datacube = state.indatacube
        state.noisecube = state.innoisecube
    else:
        utils.info(f"Smoothing datacubes by {state.smooth} spatial pixels, and "
                   f"{state.specsmooth} spectral pixels.")
        prog = utils.Progress(state.Nwpix, state.percent_step, label="smoothed", quiet=False)
        state.datacube, state.noisecube = smooth_cube(state.indatacube, state.innoisecube,
                                                      state.smooth, state.specsmooth, progress=prog)
        prog.close()
    update_badpixels(state)
    state.noise = state.noisecube


def update_badpixels(state) -> None:
    """Bad pixel mask (True = bad) from non-finite data or noise; KMOS noise > 1e5 is bad.
    A spaxel is bad (badpixelimg) only if all its wavelength pixels are bad."""
    with np.errstate(invalid="ignore"):
        bad = ~(np.isfinite(state.datacube) & np.isfinite(state.noisecube))
        if state.instr == "kmos":
            bad |= (state.noisecube > 1e5)
    state.badpixelmask = bad
    # IDL: fix(total(mask,3)/Nwpix) -> 1 only when every plane is bad
    state.badpixelimg = np.all(bad, axis=0)


def smooth_state_montecarlo(state) -> None:
    """Port of ``kubeviz_smooth, /montecarlo``: smooth the unsmoothed bootstrap or MC1
    realisations with the current spatial/spectral kernels (median only)."""
    if state.domontecarlo == ERR_BOOTSTRAP:
        src, attr = state.inbootstrapcubes, "bootstrapcubes"
    elif state.domontecarlo == ERR_MC1:
        src, attr = state.inmc1cubes, "mc1cubes"
    else:
        return
    if src is None:
        return
    if state.smooth == 1 and state.specsmooth == 1:
        setattr(state, attr, src.copy())
        return
    out = np.empty_like(src)
    prog = utils.Progress(src.shape[0], state.percent_step, label="cubes smoothed")
    for i in range(src.shape[0]):
        out[i] = smooth_cube(src[i], None, state.smooth, state.specsmooth, mode="median")
        prog.step()
    prog.close()
    setattr(state, attr, out)
