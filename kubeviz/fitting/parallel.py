"""Parallel FIT ALL / FIT ADJ ALL: spaxels are fitted by a pool of forked worker
processes and the per-spaxel results are copied back into the parent's result arrays.

Fork shares the cube copy-on-write, so nothing large is pickled; only the results
travel back (a few hundred numbers per spaxel). The workers never touch Qt. Each
spaxel is fitted by exactly the same code as the sequential loop (``dofit``); the
Monte Carlo realisations use fixed per-realisation seeds and do not depend on the
order of the fits. One pool is forked per FIT ALL / FIT ADJ ALL call; a worker that
dies raises ``BrokenProcessPool`` instead of hanging the loop. Interrupt drops the
queued chunks and lets the running ones finish.
"""
from __future__ import annotations

import multiprocessing as mp
import os
import sys
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait

import numpy as np

from .. import utils
from .linefit import dofit, testspec

_STATE = None          # the State seen by the forked workers (set just before the fork)
_RESULT_KEYS = ("n", "b", "c", "m", "nerr", "berr", "cerr", "merr", "chisq")


def default_nproc() -> int:
    """Workers used when ``state.nproc`` is 0: all cores but one, at most 8."""
    return max(1, min(8, (os.cpu_count() or 2) - 1))


def fork_available() -> bool:
    return sys.platform != "win32" and "fork" in mp.get_all_start_methods()


class FitPool:
    """Forked worker pool bound to ``state``; use as a context manager."""

    def __init__(self, state, nproc: int):
        global _STATE
        self.state = state
        self.nproc = int(nproc)
        _STATE = state
        self.executor = ProcessPoolExecutor(max_workers=self.nproc, mp_context=mp.get_context("fork"))

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()
        return False

    def close(self, cancel: bool = False):
        """Shut the pool down; with ``cancel`` the queued chunks are dropped and only the
        chunks already running finish (well under a second), so the executor is never
        left in a broken state that would stall the interpreter at exit."""
        global _STATE
        if self.executor is not None:
            self.executor.shutdown(wait=True, cancel_futures=cancel)
            self.executor = None
        _STATE = None

    def run(self, func, jobs, chunksize: int, should_cancel=None, on_results=None) -> bool:
        """Run ``func`` on chunks of ``jobs``; ``on_results(list)`` gets each finished
        chunk. Returns False when cancelled."""
        nchunks = max(self.nproc, int(np.ceil(len(jobs) / max(chunksize, 1))))
        chunks = [jobs[i::nchunks] for i in range(nchunks)]
        pending = {self.executor.submit(func, c) for c in chunks if c}
        cancelled = False
        while pending:
            done, pending = wait(pending, timeout=0.2, return_when=FIRST_COMPLETED)
            for fut in done:
                results = fut.result()
                if on_results is not None:
                    on_results(results)
            if should_cancel is not None and should_cancel():
                cancelled = True
                for fut in pending:
                    fut.cancel()
                self.close(cancel=True)
                break
        return not cancelled


def _row_results(rs, col, row):
    return tuple(getattr(rs, k)[row, col].copy() for k in _RESULT_KEYS)


def _store(rs, col, row, values):
    for key, val in zip(_RESULT_KEYS, values):
        getattr(rs, key)[row, col] = val


def _fitall_chunk(spaxels):
    """Worker: fit a list of (col, row) and return their result rows."""
    st = _STATE
    rs = st.get_results()
    out = []
    for col, row in spaxels:
        st.col, st.row = int(col), int(row)
        if testspec(st):
            dofit(st)
        out.append((col, row, _row_results(rs, col, row)))
    return out


def _fitadj_chunk(jobs):
    """Worker for one FIT ADJ ALL pass: ``jobs`` are (col, row, gnfit, gbfit) with the
    start values already derived from the neighbours by the parent."""
    st = _STATE
    rs = st.get_results()
    out = []
    for col, row, gnfit, gbfit in jobs:
        st.col, st.row = int(col), int(row)
        st.gnfit, st.gbfit = np.array(gnfit), np.array(gbfit)
        dofit(st)
        out.append((col, row, _row_results(rs, col, row)))
    return out


def fitall_parallel(state, nproc: int, should_cancel=None, on_progress=None, on_fit=None) -> bool:
    """Fit ``state.fitallrange`` with ``nproc`` forked workers. Returns False when cancelled."""
    c0, c1, r0, r1 = [int(v) for v in state.fitallrange]
    spaxels = [(col, row) for col in range(c0, c1 + 1) for row in range(r0, r1 + 1)]
    if not spaxels:
        return True
    rs = state.get_results()
    cur = (state.col, state.row)
    prog = utils.Progress(len(spaxels), state.percent_step, label="processed", callback=on_progress)
    utils.info(f"Fitting {len(spaxels)} spaxels with {nproc} worker processes")

    def on_results(results):
        for col, row, values in results:
            _store(rs, col, row, values)
        prog.step(len(results))
        if on_fit is not None and results:
            on_fit(int(results[-1][0]), int(results[-1][1]))

    try:
        with FitPool(state, nproc) as pool:
            completed = pool.run(_fitall_chunk, spaxels, 64, should_cancel, on_results)
    finally:
        state.col, state.row = cur
    prog.close()
    return completed


def fitadj_pass_parallel(pool: FitPool, rs, coords, flag, should_cancel=None, on_fit=None):
    """Fit the spaxels ``coords`` of one FIT ADJ ALL pass with the workers of ``pool``.
    The initial guesses come from the OK neighbours (``flag``) as they are now.
    Returns ``(fitted, cancelled)`` with ``fitted`` in the sequential (row-major) order."""
    from .fitall import guess_from_adjacent
    state = pool.state
    jobs = []
    for col, row in coords:
        if guess_from_adjacent(state, col, row, flag, rs):
            jobs.append((col, row, state.gnfit.copy(), state.gbfit.copy()))
    fitted = set()

    def on_results(results):
        for col, row, values in results:
            _store(rs, col, row, values)
            fitted.add((col, row))
        if on_fit is not None and results:
            on_fit(int(results[-1][0]), int(results[-1][1]))

    cur = (state.col, state.row)
    try:
        completed = pool.run(_fitadj_chunk, jobs, 32, should_cancel, on_results) if jobs else True
    finally:
        state.col, state.row = cur
    return [xy for xy in coords if xy in fitted], not completed
