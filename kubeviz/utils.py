"""Small numerical helpers ported from the top of ``kubeviz.pro``.

Where IDL semantics differ from numpy (median of an even number of elements,
NaN handling of ``total``) the IDL behaviour is reproduced explicitly so that results
match the original code.
"""
from __future__ import annotations

import logging
import sys
import time
from datetime import datetime

import numpy as np

log = logging.getLogger("kubeviz")


# ----------------------------------------------------------------------------- logging
def setup_logging(level=logging.INFO, logfile: str | None = None) -> None:
    """Configure the ``kubeviz`` logger to print ``[KUBEVIZ] ...`` style messages.

    ``logfile`` mirrors the IDL ``logunit`` keyword: when given, messages go to that
    file instead of the terminal.
    """
    log.setLevel(level)
    log.handlers.clear()
    handler = logging.FileHandler(logfile) if logfile else logging.StreamHandler(sys.stdout)
    handler.setFormatter(logging.Formatter("%(message)s"))
    log.addHandler(handler)
    log.propagate = False


def info(msg: str) -> None:
    log.info("[KUBEVIZ] " + msg)


def warn(msg: str) -> None:
    log.warning("[WARNING] " + msg)


def error(msg: str) -> None:
    log.error("[ ERROR ] " + msg)


def debug(msg: str) -> None:
    log.debug("[DEBUG] " + msg)


# ----------------------------------------------------------------------------- IDL-like maths
def idl_median(a, axis=None, even: bool = False):
    """IDL ``MEDIAN`` semantics on an array that may contain NaN.

    NaN values are ignored. For an even number of valid elements IDL returns the
    upper of the two central values unless ``/EVEN`` is set, in which case the two
    are averaged (numpy behaviour).
    """
    a = np.asarray(a, dtype=float)
    if even:
        with np.errstate(all="ignore"):
            return np.nanmedian(a, axis=axis)
    if axis is None:
        v = np.sort(a[np.isfinite(a)])
        return v[v.size // 2] if v.size else np.nan
    s = np.sort(a, axis=axis)                    # NaN sorted to the end
    n = np.sum(np.isfinite(a), axis=axis, keepdims=True)
    idx = np.clip(n // 2, 0, a.shape[axis] - 1)
    out = np.take_along_axis(s, idx, axis=axis)
    out = np.squeeze(out, axis=axis)
    out = np.where(np.squeeze(n, axis=axis) > 0, out, np.nan)
    return out


def total_nan(a, axis=None):
    """IDL ``total(a, /nan)``."""
    return np.nansum(a, axis=axis)


def percentile(x, perc):
    """IDL ``kubeviz_percentile``: linear interpolation between order statistics.

    Identical to ``numpy.percentile`` with the default linear method, computed on the
    finite elements only. Returns a scalar for a scalar ``perc``.
    """
    x = np.asarray(x, dtype=float).ravel()
    x = x[np.isfinite(x)]
    if x.size == 0:
        return np.nan if np.ndim(perc) == 0 else np.full(np.shape(perc), np.nan)
    return np.percentile(x, perc)


def weighted_median(array, weight) -> float:
    """IDL ``kubeviz_weighted_median``: first sorted element where the cumulative
    weight exceeds half of the total weight."""
    array = np.asarray(array, dtype=float).ravel()
    weight = np.asarray(weight, dtype=float).ravel()
    tot = np.sum(weight)
    ind = np.argsort(array, kind="stable")
    csum = np.cumsum(weight[ind])
    n = int(np.searchsorted(csum, tot / 2.0, side="right"))  # first index with sum > tot/2
    n = min(n, array.size - 1)
    return float(array[ind][n])


def remove_badvalues(array, repval=0.0):
    """Replace NaN/Inf by ``repval`` (returns a copy)."""
    out = np.array(array, dtype=float, copy=True)
    out[~np.isfinite(out)] = repval
    return out


def sigma_clip(array, nsig=3.0, niter=3):
    """IDL ``kubeviz_sigma_clip``: iterative clipping around the median.

    Returns ``(clipped_values, index)``.
    """
    array = np.asarray(array, dtype=float)
    index = np.arange(array.size)
    for _ in range(niter):
        m = idl_median(array[index])
        s = np.std(array[index], ddof=1) if index.size > 1 else 0.0
        keep = np.abs(array[index] - m) < nsig * s
        if not np.any(keep):
            break
        index = index[keep]
    return array[index], index


def dataclip(image, percentage=95.0):
    """IDL ``kubeviz_dataclip``: [lo, hi] percentile range of the finite pixels."""
    img = np.asarray(image, dtype=float).ravel()
    img = img[np.isfinite(img)]
    if img.size == 0:
        return np.array([0.0, 0.0])
    lo = 50.0 - 0.5 * percentage
    hi = 50.0 + 0.5 * percentage
    return np.percentile(img, [lo, hi])


def closest(array, values):
    """Indices of the elements of ``array`` closest to each of ``values``."""
    array = np.asarray(array, dtype=float)
    values = np.atleast_1d(np.asarray(values, dtype=float))
    return np.array([int(np.argmin(np.abs(array - v))) for v in values], dtype=int)


def set_intersection(a, b):
    return np.intersect1d(np.asarray(a, dtype=int), np.asarray(b, dtype=int))


def set_difference(a, b):
    return np.setdiff1d(np.asarray(a, dtype=int), np.asarray(b, dtype=int))


def excludewave(wave) -> np.ndarray:
    """Boolean mask of wavelengths in hardcoded bad spectral regions
    (kubeviz_excludewave: atmospheric O2 feature 12680-12710 A)."""
    wave = np.asarray(wave, dtype=float)
    bad_start = [12680.0]
    bad_end = [12710.0]
    bad = np.zeros(wave.shape, dtype=bool)
    for s, e in zip(bad_start, bad_end):
        bad |= (wave >= s) & (wave <= e)
    return bad


def getsnimage(im, nim, absolute=False):
    with np.errstate(all="ignore"):
        sn = np.asarray(im, dtype=float) / np.asarray(nim, dtype=float)
    if absolute:
        sn = np.abs(sn)
    return remove_badvalues(sn)


def timediff_seconds(datetime1: str, datetime2: str) -> float:
    """Difference (t2 - t1) in seconds between two ISO ``YYYY-MM-DDThh:mm:ss`` strings."""
    def parse(s):
        return datetime(int(s[0:4]), int(s[5:7]), int(s[8:10]),
                        int(s[11:13]), int(s[14:16]), int(s[17:19]))
    return (parse(datetime2) - parse(datetime1)).total_seconds()


def fmt_num(value, fmt=None) -> str:
    """IDL ``kubeviz_str``: compact string of a number."""
    if fmt is not None:
        return format(value, fmt).strip()
    if isinstance(value, (int, np.integer)):
        return str(int(value))
    return f"{value:g}"


def resultsstring(value, error, symmetric: bool = False) -> str:
    """IDL ``kubeviz_resultsstring``: 'value +/- err' with the -999/-998 sentinels.

    With ``symmetric`` (noise-cube errors, where the two sides are equal) the pair is
    shown as ``value ± err``; Monte Carlo / bootstrap errors keep the ``+up/-down`` form.
    """
    err = np.atleast_1d(np.asarray(error, dtype=float))
    f1, f2 = "9.2f", "8.2f"
    if value > 1e6:
        f1, f2 = "9.1f", "8.1f"
    if err[0] == -999:
        return "Not Fit"
    if err[0] == -998:
        return f"{format(value, f1).strip()} ± no errors"
    if err.size == 1 or (symmetric and err.size == 2):
        return f"{format(value, f1).strip()} ± {format(abs(err[0]), f2).strip()}"
    if err.size == 2:
        return f"{format(value, f1).strip()}+{format(err[0], f2).strip()}/{format(err[1], f2).strip()}"
    return "BUG: Wrong number of elements in error for results string!"


def _fmt_eta(eta: float) -> str:
    eta = int(eta)
    if eta <= 60:
        return f"{eta}s"
    if eta <= 3600:
        return f"{eta // 60}:{eta % 60:02d} mm:ss"
    return f"{eta // 3600}:{(eta % 3600) // 60:02d}:{(eta % 3600) % 60:02d} hh:mm:ss"


class Progress:
    """Terminal progress line with ETA (port of ``kubeviz_statusline``).

    ``callback(fraction, message)`` is called at each reported step so a GUI can
    show the same information. ``should_cancel()`` may be polled by long loops.
    """

    def __init__(self, total: int, percent_step: int = 10, label: str = "processed",
                 callback=None, quiet: bool = False):
        self.total = max(int(total), 1)
        self.percent_step = percent_step
        self.label = label
        self.callback = callback
        self.quiet = quiet
        self.count = 0
        self.t_start = time.time()
        self._last_percent = -1
        self._last_cb_percent = -1

    def _message(self, percent: int) -> str:
        elapsed = time.time() - self.t_start
        if percent < 100:
            eta = elapsed / percent * (100 - percent)
            return f"[PROGRES] {percent}% {self.label} - ETA: {_fmt_eta(eta)}"
        return f"[PROGRES] 100% {self.label} - Time elapsed: {_fmt_eta(elapsed)}"

    def step(self, n: int = 1) -> None:
        self.count += n
        percent_exact = int(100.0 * self.count / self.total)
        percent = percent_exact - percent_exact % self.percent_step
        if percent != self._last_percent and percent > 0:
            self._last_percent = percent
            if not self.quiet:
                sys.stdout.write("\r" + self._message(percent).ljust(70))
                sys.stdout.flush()
        # the GUI callback is finer (every percent) than the terminal line
        if self.callback is not None and percent_exact != self._last_cb_percent and percent_exact > 0:
            self._last_cb_percent = percent_exact
            self.callback(self.count / self.total, self._message(percent_exact))

    def close(self) -> None:
        if not self.quiet and self._last_percent >= 0:
            sys.stdout.write("\n")
            sys.stdout.flush()
