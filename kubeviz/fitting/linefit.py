"""The line-fitting core (ports of ``kubeviz_linefit_init``, ``kubeviz_chooselines``,
``kubeviz_linefit_estimatekinematics``, ``kubeviz_linefit_startminmaxconsistencycheck``,
``kubeviz_linefit_getstartvals``, ``kubeviz_linefit_setparinfo``, ``kubeviz_linefit_fit``,
``kubeviz_linefit_setfixline``, ``kubeviz_linefit_swapnarrowbroad``,
``kubeviz_linefit_keepfit``, ``kubeviz_linefit_keepmom``, ``kubeviz_linefit_fitset``,
``kubeviz_linefit_testspec``, ``kubeviz_linefit_dofit``, ``kubeviz_get_residual_spec``,
``kubeviz_linefit_getredshift``, ``kubeviz_linefit_reset*`` and ``kubeviz_change_redshift``).

The IDL ``fit`` common structure becomes :class:`FitContext`, created per lineset fit.
"""
from __future__ import annotations

import os
from dataclasses import dataclass, field

import numpy as np

from .. import utils
from ..constants import (CKMS, ERR_NOISE, FIT_GAUSS, MODE_MASK, MODE_SPAXEL, NOT_FIT,
                         SIGTOFWHM, SPEC_MEDSUB, SPEC_SLICE)
from ..core.extraction import medianspec, spectrum_for_fit
from ..linesdb import (FIXED_RATIOS, LINEFANCYNAMES, LINENAMES, LINES_REST, LINESET_MAX,
                       LINESETS, airtovac)
from .continuum import definecontinuumfitregion, fitcontinuum, getcontatlambda
from .model import ModelSpec, ngaussmodel
from .moments import linefit_moments
from .mpfit import ParInfo, mpcurvefit


# ============================================================================== context
@dataclass
class FitContext:
    """Port of the IDL ``fit`` structure for one lineset fit."""
    spec: np.ndarray | None = None
    noise: np.ndarray | None = None
    wave: np.ndarray | None = None
    weights: np.ndarray | None = None
    spec_bootstrap: np.ndarray | None = None      # (Nmc, Nokk)
    mfitpars: np.ndarray = field(default_factory=lambda: np.zeros(0, dtype=np.uint8))
    nfitpars: np.ndarray = field(default_factory=lambda: np.zeros(0, dtype=np.uint8))
    bfitpars: np.ndarray = field(default_factory=lambda: np.zeros(0, dtype=np.uint8))
    cfitpars: int = 0
    cpar: int = 0
    nparstart: int = -1
    nparend: int = -1
    bparstart: int = -1
    bparend: int = -1
    contfitregion: np.ndarray | None = None
    contpars: np.ndarray | None = None
    contpars_bootstrap: np.ndarray | None = None
    sigcontpars: np.ndarray | None = None
    sigcont: float | None = None
    pars: np.ndarray | None = None
    sigpars: np.ndarray | None = None
    pars_bootstrap: np.ndarray | None = None
    gpars: np.ndarray | None = None
    parinfo: list | None = None
    chisq: float = 0.0
    maxwoffb: float = 0.0
    maxwoffr: float = 0.0
    okk: np.ndarray | None = None
    specfit: np.ndarray | None = None
    status: int = 0
    errmsg: str = ""
    iter: int = 0
    firstset: bool = False
    linesinset: np.ndarray | None = None

    # ------------------------------------------------------------------ derived
    @property
    def nnfitpars(self) -> int:
        return int(np.sum(self.nfitpars))

    @property
    def nbfitpars(self) -> int:
        return int(np.sum(self.bfitpars))

    @property
    def npars(self) -> int:
        return 1 + self.nnfitpars + self.nbfitpars

    @property
    def nnfitlines(self) -> int:
        return max(self.nnfitpars - 2, 0)

    @property
    def nbfitlines(self) -> int:
        return max(self.nbfitpars - 2, 0)

    def nfit(self) -> np.ndarray:
        """Indices (into state.lines) of the narrow lines being fit."""
        return np.flatnonzero(self.nfitpars[2:] == 1)

    def bfit(self) -> np.ndarray:
        return np.flatnonzero(self.bfitpars[2:] == 1)

    def model_spec(self, state) -> ModelSpec:
        return ModelSpec(cpar=self.cpar, nparstart=self.nparstart, bparstart=self.bparstart,
                         nlines=state.lines[self.nfit()] if self.nnfitlines else np.zeros(0),
                         blines=state.lines[self.bfit()] if self.nbfitlines else np.zeros(0),
                         instrres=state.getinstrres)


# ============================================================================== init
def chooselines(state, linesets, lines, lineset=None):
    """Indices of the lines to fit (port of ``kubeviz_chooselines``)."""
    sel = state.selected_lineset if lineset is None else lineset
    utils.debug(f"selected_lineset: {sel}")
    inspec = np.flatnonzero((lines - state.wave[0] > 10.0) & (state.wave[-1] - lines > 10.0))
    if sel > 0:
        inset = np.flatnonzero(np.asarray(linesets) == sel)
        return np.intersect1d(inset, inspec)
    return inspec


def linefit_init(state) -> None:
    """Identify the lines in the wavelength range, define linesets and bisectors,
    (re)create the result containers and the instrumental resolution."""
    from ..core.instrres import instrres_init
    from ..constants import INSTRRES_EXTPOLY, INSTRRES_TEMPLATE, INSTRRES_VARPOLY

    state.fitallrange = np.array([0, state.Ncol - 1, 0, state.Nrow - 1], dtype=int)

    if state.vacuum:
        if state.Nlines_all == 0:
            utils.info("The instrument works in vacuum. Rest wavelengths converted to vacuum")
        lines_rest = airtovac(LINES_REST)
    else:
        lines_rest = LINES_REST.copy()
    lines = lines_rest * (1.0 + state.redshift)
    ok_all = np.flatnonzero((lines - state.wave[0] >= 0.0) & (state.wave[-1] - lines >= 0.0))

    state.Nlines = 0
    state.lines = np.zeros(0)
    state.linesets = np.zeros(0, dtype=int)
    state.linenames = []
    state.linefancynames = []
    state.lineset_ind = [np.zeros(0, dtype=int) for _ in range(LINESET_MAX)]
    state.lineset_max = LINESET_MAX
    state.line_bisectors = np.zeros((0, 2))

    if ok_all.size > 0:
        state.Nlines_all = int(ok_all.size)
        state.lines_all = lines[ok_all]
        oklines = chooselines(state, LINESETS, lines)
        if oklines.size > 0:
            state.lines = lines[oklines]
            state.linesets = LINESETS[oklines]
            state.linenames = [LINENAMES[i] for i in oklines]
            state.linefancynames = [LINEFANCYNAMES[i] for i in oklines]
            state.Nlines = int(oklines.size)
            state.lineset_ind = [np.flatnonzero(state.linesets == i + 1) for i in range(LINESET_MAX)]
        else:
            utils.warn("No valid selected lines in wavelength range!")
            utils.warn("Check the selected lineset.")
        bis = np.zeros((state.Nlines, 2))
        for i in range(state.Nlines):
            above = state.lines_all[state.lines_all > state.lines[i]]
            below = state.lines_all[(state.lines_all > 0) & (state.lines_all < state.lines[i])]
            bis[i, 0] = state.wave[0] if below.size == 0 else 0.5 * (below.max() + state.lines[i])
            bis[i, 1] = state.wave[-1] if above.size == 0 else 0.5 * (above.min() + state.lines[i])
        state.line_bisectors = bis
    else:
        state.Nlines_all = 0
        state.lines_all = np.zeros(0)
        utils.warn("No valid lines in wavelength range!")
        utils.warn("Check the redshift value.")

    # results containers are recreated (lines may have changed)
    state.reset_all_results()

    if state.spaxselect is None or state.spaxselect.shape[1:] != (state.Nrow, state.Ncol):
        state.spaxselect = np.zeros((max(state.Nmask, 1), state.Nrow, state.Ncol))
        state.Nmask = state.spaxselect.shape[0]
    if state.medspec is None or state.medspec.size != state.Nwpix:
        state.medspec = np.zeros(state.Nwpix)
        state.nmedspec = np.zeros(state.Nwpix)
    shape = (state.Nrow, state.Ncol)
    state.lineresimg = np.zeros(shape)
    state.unm_lineresimg = np.zeros(shape)
    state.lineerrresimg = np.zeros(shape)
    state.unm_lineerrresimg = np.zeros(shape)

    need = {INSTRRES_VARPOLY: np.sum(state.instrres_varpoly) == 0,
            INSTRRES_EXTPOLY: np.sum(state.instrres_extpoly) == 0,
            INSTRRES_TEMPLATE: state.instrres_tplsig == 0}.get(state.instrres_mode, True)
    if need and state.Nlines > 0:
        instrres_init(state)


def change_redshift(state, redshift: float) -> None:
    """Port of ``kubeviz_change_redshift``: reset fit selections and re-initialise."""
    for par in range(1, state.Nlines + 3):
        state.pnfix[par] = 0
        state.pbfix[par] = 0
    state.pndofit[:] = 0
    state.pbdofit[:] = 0
    state.pcdofit[:] = 0
    state.nshow[:] = 0
    state.bshow[:] = 0
    state.cshow[:] = 0
    state.fitconstr = False
    state.fitfixratios = False
    state.par_imagebutton = ""
    if state.cubesel > 3:
        state.cubesel = 0
    state.redshift = float(redshift)
    linefit_init(state)


# ============================================================================== start values
def estimatekinematics(state, ctx: FitContext, gpars, lambda0):
    """Robust [centre, width] guess for the main line (mpfitpeak-style estimate)."""
    glambda = lambda0 * (1.0 + gpars[ctx.nparstart] / CKMS)
    checkrange = 40.0 / (1.0 + state.redshift)
    near = np.flatnonzero(np.abs(ctx.wave - glambda) <= checkrange)
    nx = near.size
    if nx < 5:
        return _guessfail(state)
    x = ctx.wave[near]
    y = ctx.spec[near]
    order = np.argsort(x)
    xs, ys = x[order], y[order]
    fin = np.isfinite(ys)
    if not np.any(fin):
        return _guessfail(state)
    maxy, miny = np.nanmax(ys), np.nanmin(ys)
    maxx, minx = xs[-1], xs[0]
    dx = 0.5 * np.concatenate([[xs[1] - xs[0]], xs[2:] - xs[:-2], [xs[-1] - xs[-2]]])
    totarea = np.nansum(dx * ys)
    av = totarea / (maxx - minx) if maxx > minx else np.nan
    if not np.isfinite(av) or miny == maxy:
        return _guessfail(state)
    wh1 = np.flatnonzero(ys >= av)
    wh2 = np.flatnonzero(ys <= av)
    if wh1.size == 0 or wh2.size == 0:
        return _guessfail(state)
    cent = xs[int(np.nanargmax(ys))]
    peak = maxy - av
    peakarea = totarea - np.nansum(dx * np.minimum(ys, av))
    if peak == 0:
        peak = 0.5 * peakarea
    width = peakarea / (2.0 * abs(peak)) if peak != 0 else 0.0
    if width == 0 or not np.isfinite(width):
        width = float(np.median(dx))
    return np.array([cent, width])


def _guessfail(state):
    utils.debug("GUESS FOR SIGMA FAILED. SIGMA AND CENTROID FIXED")
    return np.array([-1.0, -1.0])


def startminmaxconsistencycheck(state, guess, par, linetype, silent=False):
    """Clip a start value into the user limits when 'fit with constraints' is on."""
    setguess = guess
    update = False
    if state.fitconstr:
        lims = state.gnlims if linetype == "N" else state.gblims
        gfit = state.gnfit if linetype == "N" else state.gbfit
        label = "Narrow" if linetype == "N" else "Broad"
        if lims[par, 0] != NOT_FIT and guess < lims[par, 0]:
            if not silent:
                utils.warn(f"{label} Line Parameter {par:2d} initial guess value = {guess:7.2f} "
                           f"is below user-defined min value of {lims[par, 0]:7.2f}")
                utils.warn(f"Initial guess set to {lims[par, 0]:7.2f}")
            setguess = lims[par, 0]
            gfit[par] = setguess
            update = not silent
        if lims[par, 1] != NOT_FIT and guess > lims[par, 1]:
            if not silent:
                utils.warn(f"{label} Line Parameter {par:2d} initial guess value = {guess:7.2f} "
                           f"is above user-defined max value of {lims[par, 1]:7.2f}")
                utils.warn(f"Initial guess set to {lims[par, 1]:7.2f}")
            setguess = lims[par, 1]
            gfit[par] = setguess
            update = not silent
    if update and state.on_userpars_changed is not None:
        state.on_userpars_changed()
    return setguess


def getstartvals(state, ctx: FitContext, col=None, row=None, silent=False) -> None:
    """Fill ``ctx.gpars`` with the initial guesses (port of ``kubeviz_linefit_getstartvals``)."""
    npars = ctx.npars
    nnl, nbl = ctx.nnfitlines, ctx.nbfitlines
    nfit, bfit = ctx.nfit(), ctx.bfit()
    nmainline = bmainline = None
    nmainlineind = bmainlineind = -1
    if nnl > 0:
        nmainline = state.mainline(linesarr=[state.linenames[i] for i in nfit])
        d = np.abs(state.lines[nfit] - nmainline)
        nmainlineind = int(np.argmin(d)) if d.min() < 2 else 0
    if nbl > 0:
        bmainline = state.mainline(linesarr=[state.linenames[i] for i in bfit])
        d = np.abs(state.lines[bfit] - bmainline)
        bmainlineind = int(np.argmin(d)) if d.min() < 2 else 0

    gpars = np.zeros(npars)
    if state.continuumfit_mode == 1 and ctx.contpars is not None:
        gpars[ctx.cpar] = ctx.contpars[0]

    fluxestimates0 = np.zeros(0)
    if col is not None and row is not None and state.gauss_initmap is not None:
        offsetinkms = float(state.gauss_initmap[0, row, col])
        nwidthinkms = float(state.gauss_initmap[1, row, col])
        fluxestimates0 = np.asarray(state.gauss_initmap[2:, row, col], dtype=float)
        utils.debug(f"Starting values from the external map: {state.gauss_initmap[:, row, col]}")
    else:
        refline = nmainline if nmainline is not None else bmainline
        kin = estimatekinematics(state, ctx, gpars, refline)
        if kin[0] == -1:
            offsetinkms, nwidthinkms = 0.0, 25.0
        else:
            offsetinkms = (kin[0] - refline) / refline * CKMS
            nwidthinkms = (kin[1] / refline) * CKMS
            if offsetinkms < 1e-3:
                offsetinkms = 0.0
            if nwidthinkms < 1e-3:
                nwidthinkms = 0.0
    nvalidflux = int(np.sum(fluxestimates0 != NOT_FIT)) if fluxestimates0.size else 0

    # --- velocity offset (par 1)
    par = 1
    if nbl > 0 and state.secondcomp_mode <= 1:
        noffset, boffset = offsetinkms - 100.0, offsetinkms + 100.0
    else:
        noffset, boffset = offsetinkms, offsetinkms
    if state.gnfit[par] != NOT_FIT:
        noffset = state.gnfit[par]
    if state.gbfit[par] != NOT_FIT:
        boffset = state.gbfit[par]
    if nnl > 0:
        gpars[ctx.nparstart] = startminmaxconsistencycheck(state, noffset, par, "N", silent)
    if nbl > 0:
        gpars[ctx.bparstart] = startminmaxconsistencycheck(state, boffset, par, "B", silent)

    # --- velocity width (par 2)
    par = 2
    bwidthinkms = 3.0 * nwidthinkms if (nbl > 0 and state.secondcomp_mode == 2) else nwidthinkms
    if state.gnfit[par] != NOT_FIT:
        nwidthinkms = state.gnfit[par]
    if state.gbfit[par] != NOT_FIT:
        bwidthinkms = state.gbfit[par]
    if nnl > 0:
        gpars[ctx.nparstart + 1] = startminmaxconsistencycheck(state, nwidthinkms, par, "N", silent)
    if nbl > 0:
        gpars[ctx.bparstart + 1] = startminmaxconsistencycheck(state, bwidthinkms, par, "B", silent)

    # --- fluxes: main line first, other lines scale from it (IDL processed the lines in
    #     order, so lines preceding the main line started from a zero flux guess)
    def flux_guesses(fit_idx, mainind, parstart, mainline, gfit, linetype):
        order = [mainind] + [i for i in range(len(fit_idx)) if i != mainind]
        for i in order:
            line = state.lines[fit_idx[i]]
            par = int(fit_idx[i]) + 3
            p = parstart + 2 + i
            if gfit[par] != NOT_FIT:
                gpars[p] = startminmaxconsistencycheck(state, gfit[par], par, linetype, silent)
                continue
            if nvalidflux > 0 and i < fluxestimates0.size and fluxestimates0[i] != NOT_FIT:
                gpars[p] = fluxestimates0[i]
                continue
            if i == mainind:
                goff_A = (gpars[parstart] / CKMS) * line
                gwidth_A = (gpars[parstart + 1] / CKMS) * mainline
                wr = (ctx.wave >= line + goff_A - 10.0) & (ctx.wave <= line + goff_A + 10.0)
                fluxestimate = float(np.nanmax(ctx.spec[wr])) * np.sqrt(2.0 * np.pi) * gwidth_A if np.any(wr) else 0.0
                if not np.isfinite(fluxestimate):
                    fluxestimate = 0.0
            else:
                fluxestimate = 0.33 * gpars[parstart + 2 + mainind]
            gpars[p] = startminmaxconsistencycheck(state, fluxestimate, par, linetype, silent)

    if nnl > 0:
        flux_guesses(nfit, nmainlineind, ctx.nparstart, nmainline, state.gnfit, "N")
    if nbl > 0:
        flux_guesses(bfit, bmainlineind, ctx.bparstart, bmainline, state.gbfit, "B")

    ctx.gpars = gpars
    utils.debug(f"Starting values: {gpars}")


# ============================================================================== constraints
def setparinfo(state, ctx: FitContext) -> None:
    """Build the ``parinfo`` list (port of ``kubeviz_linefit_setparinfo``)."""
    npars = ctx.npars
    parinfo = [ParInfo() for _ in range(npars)]
    parinfo[ctx.cpar].fixed = (state.continuumfit_mode == 0)

    nfitp = np.flatnonzero(ctx.nfitpars == 1)
    bfitp = np.flatnonzero(ctx.bfitpars == 1)

    if state.fitconstr or state.selected_lineset == 0:
        for lst, fix, lims, start in ((nfitp, state.pnfix, state.gnlims, ctx.nparstart),
                                      (bfitp, state.pbfix, state.gblims, ctx.bparstart)):
            for i, fp in enumerate(lst):
                par = int(fp) + 1
                p = start + i
                if fix[par] == 1:
                    parinfo[p].fixed = True
                if lims[par, 0] != NOT_FIT:
                    parinfo[p].limited[0] = True
                    parinfo[p].limits[0] = float(lims[par, 0])
                elif par > 1:
                    parinfo[p].limited[0] = True
                    parinfo[p].limits[0] = 0.0
                if lims[par, 1] != NOT_FIT:
                    parinfo[p].limited[1] = True
                    parinfo[p].limits[1] = float(lims[par, 1])
    else:
        for n, start in ((ctx.nnfitpars, ctx.nparstart), (ctx.nbfitpars, ctx.bparstart)):
            if n > 2:
                for i in range(1, n):
                    parinfo[start + i].limited[0] = True
                    parinfo[start + i].limits[0] = 0.0

    if state.fitfixratios:
        for lst, start in ((nfitp, ctx.nparstart), (bfitp, ctx.bparstart)):
            if lst.size == 0:
                continue
            pos = {}
            for i in range(2, lst.size):
                name = state.linenames[int(lst[i]) - 2]
                pos[name] = start + i
            for (faint, bright), ratio in FIXED_RATIOS.items():
                if faint in pos and bright in pos:
                    parinfo[pos[faint]].tied = f"P[{pos[bright]}]/{ratio}"
    ctx.parinfo = parinfo


# ============================================================================== fitting
def _errmsg_path(state) -> str:
    return os.path.join(state.cwdir or ".", f"fit_errmsg_{state.basename()}.txt")


def _log_fit_error(state, errmsg: str) -> None:
    try:
        with open(_errmsg_path(state), "a") as fh:
            if state.linefit_mode == MODE_SPAXEL:
                fh.write(f"Spaxel: {state.col} {state.row}\n")
            else:
                fh.write(f"Mask: {state.imask - 1}\n")
            fh.write(f"Fitting Error Message: {errmsg}\n")
    except OSError:
        pass
    utils.debug(f"Fitting Error Message: {errmsg}")


def linefit_fit(state, ctx: FitContext, wave, spec, weights, gpars, savepars=False):
    """Run mpcurvefit on one spectrum. Returns the best-fit parameter vector
    (all NOT_FIT on failure)."""
    spec_model = ctx.model_spec(state)
    res = mpcurvefit(ngaussmodel, wave, spec, weights, gpars, parinfo=ctx.parinfo, itmax=200,
                     quiet=not state.debug, spec=spec_model)
    if res.errmsg:
        _log_fit_error(state, res.errmsg)
        pars = np.full(ctx.npars, NOT_FIT)
        if savepars:
            ctx.pars = pars.copy()
            ctx.sigpars = np.full(ctx.npars, NOT_FIT)
            ctx.errmsg = res.errmsg
            ctx.status = res.status
        return pars
    if savepars:
        ctx.iter = res.niter
        ctx.specfit = res.yfit
        ctx.pars = res.params.copy()
        ctx.sigpars = res.perror.copy()
        ctx.chisq = res.chisq / res.dof if res.dof > 0 else res.chisq
        ctx.status = res.status
        ctx.errmsg = res.errmsg
    if state.debug:
        utils.debug(f"Starting pars: {gpars}")
        utils.debug(f"Pars: {res.params}")
        utils.debug(f"chisq: {res.chisq} dof: {res.dof} status: {res.status} niter: {res.niter}")
    return res.params.copy()


def setfixline(state, ctx: FitContext, linetypes) -> None:
    """Fix line position and width to the current best fit for the following linesets."""
    for lt in linetypes:
        if lt == "N" and ctx.nparstart >= 0:
            state.pnfix[1] = 1
            state.pnfix[2] = 1
            state.gnfit[1] = ctx.pars[ctx.nparstart]
            state.gnfit[2] = ctx.pars[ctx.nparstart + 1]
        elif lt == "B" and ctx.bparstart >= 0:
            state.pbfix[1] = 1
            state.pbfix[2] = 1
            state.gbfit[1] = ctx.pars[ctx.bparstart]
            state.gbfit[2] = ctx.pars[ctx.bparstart + 1]


def swapnarrowbroad(ctx: FitContext) -> None:
    n = ctx.bparstart - ctx.nparstart
    a, b = slice(ctx.nparstart, ctx.bparstart), slice(ctx.bparstart, ctx.bparstart + n)
    for arr in (ctx.pars, ctx.sigpars):
        tmp = arr[a].copy()
        arr[a] = arr[b]
        arr[b] = tmp
    if ctx.pars_bootstrap is not None:
        tmp = ctx.pars_bootstrap[:, a].copy()
        ctx.pars_bootstrap[:, a] = ctx.pars_bootstrap[:, b]
        ctx.pars_bootstrap[:, b] = tmp


# ============================================================================== keep results
def _current_index(state, rs):
    if state.linefit_mode == MODE_SPAXEL:
        return rs.index(col=state.col, row=state.row)
    return rs.index(imask=state.imask)


def keepfit(state, ctx: FitContext) -> None:
    """Copy the Gaussian fit results into the result containers (port of ``kubeviz_linefit_keepfit``)."""
    if state.scaleNoiseerrors and state.domontecarlo == ERR_NOISE and ctx.chisq > 0:
        ctx.sigpars = ctx.sigpars * np.sqrt(ctx.chisq)

    nnfitpars, nbfitpars, ncfitpars = ctx.nnfitpars, ctx.nbfitpars, int(ctx.cfitpars)
    nnl, nbl = ctx.nnfitlines, ctx.nbfitlines

    if nbl > 0 and nnl == nbl:
        pn, pb = ctx.pars, ctx.pars
        if state.secondcomp_mode == 0 and pn[ctx.nparstart + 2] < pb[ctx.bparstart + 2]:
            swapnarrowbroad(ctx)
        elif state.secondcomp_mode == 1 and pn[ctx.nparstart] > pb[ctx.bparstart]:
            swapnarrowbroad(ctx)
        elif state.secondcomp_mode == 2 and pn[ctx.nparstart + 1] > pb[ctx.bparstart + 1]:
            swapnarrowbroad(ctx)

    nfitp = np.flatnonzero(ctx.nfitpars == 1)
    bfitp = np.flatnonzero(ctx.bfitpars == 1)

    rs = state.get_results()
    idx = _current_index(state, rs)
    nres, nerr = rs.n[idx].copy(), rs.nerr[idx].copy()
    bres, berr = rs.b[idx].copy(), rs.berr[idx].copy()
    cres, cerr = rs.c[idx].copy(), rs.cerr[idx].copy()

    mc = state.domontecarlo > 0
    percs = state.montecarlo_percs
    cont_lambda_boot = np.zeros((state.Nmontecarlo, nnl + nbl)) if mc else None
    cont_boot_par = 0
    plot = state.plotMonteCarlodistrib and mc
    if plot:
        from ..io.mcpdf import plot_montecarlo_distrib

    def store_component(fitp, start, res, err, restype, fixarr, dv_ref):
        nonlocal cont_boot_par
        dv = 0.0
        for i, fp in enumerate(fitp):
            par = int(fp) + 1
            p = start + i
            if ctx.pars[p] != NOT_FIT and ctx.sigpars[p] != NOT_FIT:
                res[par] = ctx.pars[p]
                if i == 0:
                    dv = ctx.pars[p]
                if restype == "narrow":
                    skip = (par <= 2 and state.selected_lineset == 0 and not ctx.firstset)
                else:
                    skip = (par <= 2 and state.selected_lineset == 0 and fixarr[par] == 1 and err[par, 0] > 0)
                if not skip:
                    if mc:
                        dist = ctx.pars_bootstrap[:, p]
                        pb = utils.percentile(dist, percs)
                        err[par, :] = pb - ctx.pars[p]
                        if plot:
                            plot_montecarlo_distrib(state, dist, pb, ctx.pars[p], restype, par)
                    else:
                        err[par, 0] = ctx.sigpars[p]
                        err[par, 1] = -ctx.sigpars[p]
            if par > 2 and ncfitpars > 0:
                cpar = par - 2
                line = cpar - 1
                lambda_line = state.lines[line] * (1.0 + dv / CKMS)
                if state.continuumfit_mode == 0:
                    cres[cpar] = float(getcontatlambda(ctx.contpars, lambda_line))
                else:
                    cres[cpar] = ctx.pars[ctx.cpar]
                if mc:
                    for boot in range(state.Nmontecarlo):
                        if state.continuumfit_mode == 0:
                            cont_lambda_boot[boot, cont_boot_par] = getcontatlambda(ctx.contpars_bootstrap[boot], lambda_line)
                        else:
                            cont_lambda_boot[boot, cont_boot_par] = ctx.pars_bootstrap[boot, ctx.cpar]
                    pb = utils.percentile(cont_lambda_boot[:, cont_boot_par], percs)
                    cerr[cpar, :] = pb - cres[cpar]
                    if plot:
                        plot_montecarlo_distrib(state, cont_lambda_boot[:, cont_boot_par], pb, cres[cpar], "cont", cpar)
                    cont_boot_par += 1
                else:
                    cerr[cpar, 0] = ctx.sigcont if state.continuumfit_mode == 0 else ctx.sigpars[ctx.cpar]
                    cerr[cpar, 1] = -cerr[cpar, 0]

    if nnfitpars > 0:
        store_component(nfitp, ctx.nparstart, nres, nerr, "narrow", state.pnfix, None)
    if nbfitpars > 0:
        store_component(bfitp, ctx.bparstart, bres, berr, "broad", state.pbfix, None)

    if nnfitpars + nbfitpars > 0:
        rs.chisq[idx] = ctx.chisq
    rs.n[idx], rs.nerr[idx] = nres, nerr
    rs.b[idx], rs.berr[idx] = bres, berr
    rs.c[idx], rs.cerr[idx] = cres, cerr

    if state.saveMonteCarlodistrib and mc:
        from ..io.mcpdf import save_montecarlo_distrib
        save_montecarlo_distrib(state, ctx, ctx.pars_bootstrap, cont_lambda_boot)


def keepmom(state, ctx: FitContext) -> None:
    """Copy the moment results into the containers (port of ``kubeviz_linefit_keepmom``)."""
    mfitp = np.flatnonzero(ctx.mfitpars == 1)
    ncfitpars = int(ctx.cfitpars)
    rs = state.get_results()
    idx = _current_index(state, rs)
    mres, merr = rs.m[idx].copy(), rs.merr[idx].copy()
    cres, cerr = rs.c[idx].copy(), rs.cerr[idx].copy()
    mc = state.domontecarlo > 0
    percs = state.montecarlo_percs
    cont_lambda_boot = np.zeros((state.Nmontecarlo, mfitp.size)) if mc else None
    plot = state.plotMonteCarlodistrib and mc
    if plot:
        from ..io.mcpdf import plot_montecarlo_distrib

    for i, mf in enumerate(mfitp):
        for j in range(5):
            par = 6 * int(mf) + j
            npar = 5 * i + j
            mres[par] = ctx.pars[npar]
            if mc:
                pb = utils.percentile(ctx.pars_bootstrap[:, npar], percs)
                merr[par, :] = pb - ctx.pars[npar]
            else:
                merr[par, :] = ctx.sigpars[npar]
        if ncfitpars > 0:
            cpar = int(mf) + 1
            lambda_line = state.lines[int(mf)] * (1.0 + mres[6 * int(mf) + 1] / CKMS)
            cres[cpar] = float(getcontatlambda(ctx.contpars, lambda_line))
            if mc:
                for boot in range(state.Nmontecarlo):
                    cont_lambda_boot[boot, i] = getcontatlambda(ctx.contpars_bootstrap[boot], lambda_line)
                pb = utils.percentile(cont_lambda_boot[:, i], percs)
                cerr[cpar, :] = pb - cres[cpar]
                if plot:
                    plot_montecarlo_distrib(state, cont_lambda_boot[:, i], pb, cres[cpar], "cont", cpar)
            else:
                cerr[cpar, 0] = ctx.sigcont if ctx.sigcont is not None else 0.0
                cerr[cpar, 1] = -cerr[cpar, 0]

    rs.m[idx], rs.merr[idx] = mres, merr
    rs.c[idx], rs.cerr[idx] = cres, cerr
    if state.saveMonteCarlodistrib and mc:
        from ..io.mcpdf import save_montecarlo_distrib
        save_montecarlo_distrib(state, ctx, ctx.pars_bootstrap, cont_lambda_boot)


# ============================================================================== fitset / dofit
def fitset(state, wave, spec, noise, linesinset, spec_bootstrap, keep: bool, firstset: bool,
           col=None, row=None, silent=False):
    """Fit one lineset. Returns ``(ctx or None, firstset)``.

    ``firstset`` tracks whether a lineset with valid lines has already been fit; the
    kinematics of subsequent sets are then fixed to the first solution.
    """
    linesinset = np.asarray(linesinset, dtype=int)
    if linesinset.size == 0:
        return None, firstset
    ctx = FitContext()
    ctx.linesinset = linesinset
    wmainline = utils.idl_median(state.lines[linesinset], even=True)
    ctx.maxwoffb, ctx.maxwoffr = state.maxwoffb, state.maxwoffr
    wmin, wmax = wmainline - ctx.maxwoffb, wmainline + ctx.maxwoffr
    useline = linesinset[(state.pndofit[linesinset] + state.pbdofit[linesinset]) >= 1]
    if useline.size > 0:
        wmin = max(state.line_bisectors[useline, 0].min(), wmin)
        wmax = min(state.line_bisectors[useline, 1].max(), wmax)
    okk = np.flatnonzero((wave >= wmin) & (wave <= wmax) & (~utils.excludewave(wave)))
    if okk.size == 0:
        return None, firstset

    nl = state.Nlines
    nfitpars = np.zeros(2 + nl, dtype=np.uint8)
    bfitpars = np.zeros(2 + nl, dtype=np.uint8)
    for line in linesinset:
        if state.pndofit[line] == 1:
            nfitpars[line + 2] = 1
        if state.pbdofit[line] == 1:
            bfitpars[line + 2] = 1
    nnl, nbl = int(nfitpars.sum()), int(bfitpars.sum())
    if nnl + nbl == 0:
        return None, firstset
    if nnl > 0:
        nfitpars[:2] = 1
    if nbl > 0:
        bfitpars[:2] = 1
    ctx.nfitpars, ctx.bfitpars = nfitpars, bfitpars
    ctx.cfitpars = int(state.pcdofit[linesinset[0]] == 1)
    ctx.okk = okk
    ctx.wave = np.asarray(wave, dtype=float)[okk]
    ctx.spec = np.asarray(spec, dtype=float)[okk]
    ctx.noise = np.asarray(noise, dtype=float)[okk]
    ctx.spec_bootstrap = spec_bootstrap[:, okk] if spec_bootstrap is not None else None

    ctx.cpar = 0
    if nnl > 0:
        ctx.nparstart, ctx.nparend = 1, 1 + nnl + 1
    if nbl > 0:
        ctx.bparstart = 1 if ctx.nparend == -1 else ctx.nparend + 1
        ctx.bparend = ctx.bparstart + nbl + 1

    nokk = okk.size
    nmc = state.Nmontecarlo if state.domontecarlo > 0 else 0
    contsub = np.zeros(nokk)
    contsub_boot = np.zeros((nmc, nokk)) if nmc else None
    if ctx.cfitpars == 1:
        region, ncont = definecontinuumfitregion(state, wave, wmainline)
        ctx.contfitregion = region
        if ncont > 0:
            contpars, sigcontpars, sigcont = fitcontinuum(state, region, spec, noise, wave, wmainline)
            contsub = getcontatlambda(contpars, ctx.wave)
            ctx.contpars, ctx.sigcontpars, ctx.sigcont = contpars, sigcontpars, sigcont
            if nmc:
                ctx.contpars_bootstrap = np.full((nmc, state.continuumfit_order + 1), NOT_FIT)
                for boot in range(nmc):
                    cp, _, _ = fitcontinuum(state, region, spec_bootstrap[boot], noise, wave, wmainline)
                    ctx.contpars_bootstrap[boot] = cp
                    contsub_boot[boot] = getcontatlambda(cp, ctx.wave)
        else:
            if state.continuumfit_mode == 0:
                utils.warn("No suitable continuum region found. No fitting takes place.")
                return None, firstset
            ctx.contpars = np.zeros(5)
            ctx.sigcont = 0.0
    with np.errstate(divide="ignore", invalid="ignore"):
        weights = 1.0 / ctx.noise ** 2
    weights[~np.isfinite(weights)] = 0.0
    weights[ctx.noise == 0] = 0.0
    ctx.weights = weights

    if state.linefit_type == FIT_GAUSS:
        target = ctx.spec - contsub if state.continuumfit_mode == 0 else ctx.spec
        getstartvals(state, ctx, col=col, row=row, silent=silent)
        setparinfo(state, ctx)
        linefit_fit(state, ctx, ctx.wave, target, ctx.weights, ctx.gpars, savepars=True)

        # smart second component: revert to a single component when the two are not
        # genuinely separated (Consolandi et al. 2017 criteria)
        if nbl > 0 and not firstset and state.secondcomp_smart:
            vel = np.array([ctx.pars[ctx.nparstart], ctx.pars[ctx.bparstart]])
            sig = np.array([ctx.pars[ctx.nparstart + 1], ctx.pars[ctx.bparstart + 1]])
            with np.errstate(all="ignore"):
                sn = np.array([ctx.pars[ctx.nparstart + 2] / ctx.sigpars[ctx.nparstart + 2],
                               ctx.pars[ctx.bparstart + 2] / ctx.sigpars[ctx.bparstart + 2]])
            if (sig.max() > abs(vel[0] - vel[1]) / 2.0 or abs(vel[0] - vel[1]) < 75
                    or np.nanmin(sn) < 3 or not np.all(np.isfinite(sn))):
                nbl = 0
                ctx.bfitpars[:] = 0
                ctx.bparstart = ctx.bparend = -1
                getstartvals(state, ctx, col=col, row=row, silent=silent)
                setparinfo(state, ctx)
                linefit_fit(state, ctx, ctx.wave, target, ctx.weights, ctx.gpars, savepars=True)

        if nmc:
            pars_boot = np.full((nmc, ctx.npars), NOT_FIT)
            for boot in range(nmc):
                sb = ctx.spec_bootstrap[boot] - contsub_boot[boot] if state.continuumfit_mode == 0 else ctx.spec_bootstrap[boot]
                pars_boot[boot] = linefit_fit(state, ctx, ctx.wave, sb, ctx.weights, ctx.gpars)
            ctx.pars_bootstrap = pars_boot

        if not firstset:
            firstset = True
            ctx.firstset = True
        else:
            ctx.firstset = False
        linetypes = (["N"] if nnl > 0 else []) + (["B"] if nbl > 0 else [])
        if ctx.pars is not None:
            setfixline(state, ctx, linetypes)
        if keep:
            keepfit(state, ctx)
    else:
        ctx.mfitpars = np.roll(np.concatenate([[0, 0], ctx.nfitpars[2:]]), -2).astype(np.uint8)
        if nnl > 0:
            pars, sigpars = linefit_moments(state, ctx, ctx.wave, ctx.spec - contsub, ctx.weights, col=col, row=row)
            ctx.pars, ctx.sigpars = pars, sigpars
            if nmc:
                pars_boot = np.full((nmc, pars.size), NOT_FIT)
                for boot in range(nmc):
                    pb, _ = linefit_moments(state, ctx, ctx.wave, ctx.spec_bootstrap[boot] - contsub_boot[boot],
                                            ctx.weights, col=col, row=row, boot=boot)
                    pars_boot[boot] = pb
                ctx.pars_bootstrap = pars_boot
            if keep:
                keepmom(state, ctx)
    return ctx, firstset


def testspec(state, col=None, row=None) -> bool:
    """Is there any finite data in the spectrum that would be fit? (``kubeviz_linefit_testspec``)"""
    col = state.col if col is None else col
    row = state.row if row is None else row
    if state.specmode in (SPEC_SLICE, SPEC_MEDSUB):
        spec = state.datacube[:, row, col]
    else:
        medianspec(state)
        spec = state.medspec
    return bool(np.any(np.isfinite(spec)))


def dofit(state, col=None, row=None, nosave: bool = False):
    """Fit all linesets for the current spaxel/mask (port of ``kubeviz_linefit_dofit``).

    Returns the list of :class:`FitContext` objects (one per lineset actually fit).
    """
    col = state.col if col is None else col
    row = state.row if row is None else row
    keep = not nosave
    spec, noise, spec_boot = spectrum_for_fit(state, col, row)
    okspec = np.flatnonzero(np.isfinite(spec))
    contexts = []
    if okspec.size == 0:
        return contexts
    saved = (state.pnfix.copy(), state.gnfit.copy(), state.pbfix.copy(), state.gbfit.copy())
    firstset = False
    wave = np.asarray(state.wave, dtype=float)[okspec]
    sb = spec_boot[:, okspec] if spec_boot is not None else None
    for s in range(state.lineset_max):
        ctx, firstset = fitset(state, wave, spec[okspec], noise[okspec], state.lineset_ind[s], sb,
                               keep, firstset, col=col, row=row)
        if ctx is not None:
            contexts.append(ctx)
    state.pnfix, state.gnfit, state.pbfix, state.gbfit = saved
    return contexts


# ============================================================================== residuals / redshift
def get_residual_spec(state, col=None, row=None) -> np.ndarray:
    """Spectrum minus continuum (at the main line) and fitted Gaussians
    (port of ``kubeviz_get_residual_spec``)."""
    col = state.col if col is None else col
    row = state.row if row is None else row
    spec = np.asarray(state.datacube[:, row, col], dtype=float)
    residual = spec.copy()
    if not np.any(np.isfinite(spec)):
        return residual
    contexts = dofit(state, col=col, row=row, nosave=True)
    xx = np.asarray(state.wave, dtype=float)
    for ctx in contexts:
        if ctx.pars is None or ctx.nnfitlines + ctx.nbfitlines == 0:
            continue
        pars = ctx.pars.copy()
        pars[pars == NOT_FIT] = 0.0
        if ctx.nnfitlines > 0 and ctx.contpars is not None and state.continuumfit_mode == 0:
            lambda_line = state.lines[ctx.nfit()[0]] * (1.0 + pars[ctx.nparstart] / CKMS)
            residual -= float(getcontatlambda(ctx.contpars, lambda_line))
        spec_model = ctx.model_spec(state)
        model = ngaussmodel(xx, pars, spec_model)
        if state.continuumfit_mode == 0:
            model = model - pars[ctx.cpar]
        residual -= model
    return residual


def getredshift(state) -> float:
    """Redshift from a fit of the summed central spectrum (port of ``kubeviz_linefit_getredshift``).
    Returns -1 when no significant line is found."""
    saved = dict(domontecarlo=state.domontecarlo, pndofit=state.pndofit.copy(), pbdofit=state.pbdofit.copy(),
                 pcdofit=state.pcdofit.copy(), maxwoffb=state.maxwoffb, maxwoffr=state.maxwoffr,
                 minoff=state.continuumfit_minoff, maxoff=state.continuumfit_maxoff,
                 minperc=state.continuumfit_minperc, maxperc=state.continuumfit_maxperc,
                 order=state.continuumfit_order, pnfix=state.pnfix.copy(), gnfit=state.gnfit.copy(),
                 pbfix=state.pbfix.copy(), gbfit=state.gbfit.copy())
    for iline in range(state.Nlines):
        state.pndofit[iline] = 1
        state.pbdofit[iline] = 0
        state.pcdofit[iline] = 1
    state.domontecarlo = ERR_NOISE
    state.set_continuumfit_defaults()
    state.set_linefitrange_defaults()

    medianspec(state, sum=True, redsh=True)
    spec = np.asarray(state.medspec, dtype=float)
    noise = np.asarray(state.nmedspec, dtype=float)
    wave = np.asarray(state.wave, dtype=float)
    okspec = np.flatnonzero(np.isfinite(spec))
    redshift = -1.0
    main_ind = state.mainline_index()
    if okspec.size > 0:
        firstset = False
        best = None
        for s in range(state.lineset_max):
            ctx, firstset = fitset(state, wave[okspec], spec[okspec], noise[okspec], state.lineset_ind[s],
                                   None, False, firstset)
            if ctx is not None and main_ind in ctx.nfit():
                best = ctx
        if best is not None and best.pars is not None and best.pars[best.nparstart + 2] != NOT_FIT:
            i_main = int(np.flatnonzero(best.nfit() == main_ind)[0])
            flux = best.pars[best.nparstart + 2 + i_main]
            sig = best.sigpars[best.nparstart + 2 + i_main]
            fitres = state.lines[main_ind] * (1.0 + best.pars[best.nparstart] / CKMS)
            if sig > 0 and flux / sig > state.mask_sn_thresh and flux > 0:
                redshift = fitres / state.mainline(redshift=0.0) - 1.0

    state.domontecarlo = saved["domontecarlo"]
    state.pndofit, state.pbdofit, state.pcdofit = saved["pndofit"], saved["pbdofit"], saved["pcdofit"]
    state.maxwoffb, state.maxwoffr = saved["maxwoffb"], saved["maxwoffr"]
    state.continuumfit_minoff, state.continuumfit_maxoff = saved["minoff"], saved["maxoff"]
    state.continuumfit_minperc, state.continuumfit_maxperc = saved["minperc"], saved["maxperc"]
    state.continuumfit_order = saved["order"]
    state.pnfix, state.gnfit, state.pbfix, state.gbfit = saved["pnfix"], saved["gnfit"], saved["pbfix"], saved["gbfit"]
    return float(redshift)


# ============================================================================== resets
def linefit_reset(state, linetypes=("N", "B", "C"), pars=None) -> None:
    """Reset the fit of the current spaxel/mask (``kubeviz_linefit_reset``)."""
    rs = state.get_results()
    rs.reset_element(_current_index(state, rs), linetypes, pars)


def linefit_resetall(state, linetypes=("N", "B", "C"), pars=None) -> None:
    """Reset the fits in the fitall range / current mask (``kubeviz_linefit_resetall``)."""
    rs = state.get_results()
    if state.linefit_mode == MODE_SPAXEL:
        rs.reset_region(state.fitallrange, linetypes, pars)
    else:
        rs.reset_element(_current_index(state, rs), linetypes, pars)


def linefit_resetuser(state, linetypes=("N", "B"), pars=None, all_pars=False, startonly=False) -> None:
    """Reset user start values (and limits) (``kubeviz_linefit_resetuser``)."""
    if all_pars:
        linetypes = ("N", "B")
        pars = range(3 + state.Nlines)
    if pars is None:
        pars = range(3 + state.Nlines)
    for lt in linetypes:
        gfit = state.gnfit if lt == "N" else state.gbfit
        lims = state.gnlims if lt == "N" else state.gblims
        for par in pars:
            if par > 0:
                gfit[par] = NOT_FIT
                if not startonly:
                    lims[par, :] = NOT_FIT
    if state.on_userpars_changed is not None:
        state.on_userpars_changed()
