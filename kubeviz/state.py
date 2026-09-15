"""Session state (port of the ``kubeviz_state`` common block).

Everything the IDL code kept in the ``state`` structure lives here, minus widget ids
and colour tables (GUI concerns handled in :mod:`kubeviz.gui`). Array layout follows
numpy/astropy conventions: cubes are ``(Nwpix, Nrow, Ncol)`` and images ``(Nrow, Ncol)``.

The many ``<method>_<mode>_<x>rescube`` pointers of the IDL code are replaced by the
``results`` dictionary keyed by ``(linefit_mode, error_method)``.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Optional

import numpy as np

from . import __version__
from .constants import (CKMS, ERR_BOOTSTRAP, ERR_MC1, ERR_MC2, ERR_MC3, ERR_NOISE,
                        FIT_GAUSS, INSTRRES_EXTPOLY, INSTRRES_TEMPLATE, INSTRRES_VARPOLY,
                        MAX_POLYCOEFF_INSTRRES, MAXLINES, MODE_MASK, MODE_SPAXEL,
                        MONTECARLO_PERCS, NOT_FIT, SIGTOFWHM)
from .linesdb import airtovac, mainline_rest
from .results import ResultSet


def _f25():
    return np.full(MAXLINES, NOT_FIT, dtype=float)


def _b25():
    return np.zeros(MAXLINES, dtype=np.uint8)


def _f25x2():
    return np.full((MAXLINES, 2), NOT_FIT, dtype=float)


@dataclass
class State:
    # ------------------------------------------------------------------ meta
    version: str = __version__
    debug: bool = False
    percent_step: int = 10
    nproc: int = 1                                  # worker processes for FIT ALL / FIT ADJ ALL (0 = automatic)
    filename: str = ""
    indir: str = ""
    cwdir: str = ""
    outdir: str = ""

    # ------------------------------------------------------------------ data
    indatacube: Optional[np.ndarray] = None     # unsmoothed input cube (Nw, Nrow, Ncol)
    innoisecube: Optional[np.ndarray] = None
    inprihead: Any = None                        # astropy Header or None
    indatahead: Any = None
    innoisehead: Any = None
    datacube: Optional[np.ndarray] = None       # smoothed cube
    noisecube: Optional[np.ndarray] = None      # smoothed noise
    noise: Optional[np.ndarray] = None          # current noise (noisecube or MC noise)
    badpixelmask: Optional[np.ndarray] = None   # bool (Nw, Nrow, Ncol), True = bad
    badpixelimg: Optional[np.ndarray] = None    # bool (Nrow, Ncol), True = whole spaxel bad
    fluxfac: float = 1.0
    noiseisvar: bool = False

    lambda0: float = 0.0
    dlambda: float = 0.0
    pix0: float = 0.0
    wave: Optional[np.ndarray] = None

    medspec: Optional[np.ndarray] = None
    nmedspec: Optional[np.ndarray] = None
    medspec_montecarlo: Optional[np.ndarray] = None   # (Nmontecarlo, Nwpix)

    # ------------------------------------------------------------------ error method
    domontecarlo: int = ERR_NOISE
    Nmontecarlo: int = 0
    montecarlo_percs: np.ndarray = field(default_factory=lambda: np.array(MONTECARLO_PERCS))
    plotMonteCarlodistrib: bool = False
    saveMonteCarlodistrib: bool = False
    useMonteCarlonoise: bool = False
    scaleNoiseerrors: bool = False

    Nbootstrap: int = 0
    bootstrap_file: str = ""                       # path of the bootstrap cubes file
    inbootstrapcubes: Optional[np.ndarray] = None   # (Nboot, Nw, Nrow, Ncol) unsmoothed
    bootstrapcubes: Optional[np.ndarray] = None
    bootstrapnoise: Optional[np.ndarray] = None
    Nmc1: int = 0
    inmc1cubes: Optional[np.ndarray] = None
    mc1cubes: Optional[np.ndarray] = None
    mc1noise: Optional[np.ndarray] = None
    Nmc2: int = 0
    mc2cubes: Optional[np.ndarray] = None
    mc2noise: Optional[np.ndarray] = None
    Nmc3: int = 0
    mc3cubes: Optional[np.ndarray] = None
    mc3noise: Optional[np.ndarray] = None

    # ------------------------------------------------------------------ display scaling
    zmin_ima: float = 0.0
    zmax_ima: float = 0.4
    zmin_spec: float = 0.0
    zmax_spec: float = 0.4

    # ------------------------------------------------------------------ instrument
    instr: str = ""
    band: str = ""
    ifu: int = 0
    vacuum: bool = False
    pixscale: int = 0

    # ------------------------------------------------------------------ navigation
    cubesel: int = 0
    col: int = 0
    row: int = 0
    wpix: int = 0
    Startcol: int = 0
    Startrow: int = 0
    Startwpix: int = 0
    Ncol: int = 0
    Nrow: int = 0
    Nwpix: int = 0

    zoommap: bool = False
    linefitmap: bool = False
    scroll: bool = False
    marker: int = 1
    scale: int = -1
    zoomrange: int = 32
    zcuts: int = 4
    cursormode: int = 1
    specmode: int = 0
    imgmode: int = 0
    flagmode: bool = True
    wavsel: int = 1
    npress: int = 0                                 # 's' key presses while selecting a range
    wavrange1: np.ndarray = field(default_factory=lambda: np.zeros(2, dtype=int))
    wavrange2: np.ndarray = field(default_factory=lambda: np.zeros(2, dtype=int))
    spaxrange: np.ndarray = field(default_factory=lambda: np.zeros(4, dtype=int))
    zoomfac: int = 0
    zoomspax: int = 1
    invert: int = -1
    ctab: int = 0

    img1: Optional[np.ndarray] = None
    img2: Optional[np.ndarray] = None
    nimg1: Optional[np.ndarray] = None
    nimg2: Optional[np.ndarray] = None
    bimg1: Optional[np.ndarray] = None
    bimg2: Optional[np.ndarray] = None
    lineresimg: Optional[np.ndarray] = None
    lineerrresimg: Optional[np.ndarray] = None
    unm_lineresimg: Optional[np.ndarray] = None
    unm_lineerrresimg: Optional[np.ndarray] = None
    par_imagebutton: str = ""

    # ------------------------------------------------------------------ smoothing
    smooth: int = 1
    specsmooth: int = 1
    transpose: bool = False

    # ------------------------------------------------------------------ masks
    Nmask: int = 1
    imask: int = 1                                  # 1-based
    maskradius: int = 1
    maskmode: int = 0                               # 0 single spaxel, 1 circular, 2 square
    spaxselect: Optional[np.ndarray] = None         # (Nmask, Nrow, Ncol) float

    # ------------------------------------------------------------------ fit configuration
    continuumfit_mode: int = 0                      # 0 SDSS sidebands, 1 fit constant with mpfit
    continuumfit_minoff: float = 200.0
    continuumfit_maxoff: float = 500.0
    continuumfit_minperc: float = 40.0
    continuumfit_maxperc: float = 60.0
    continuumfit_order: int = 0
    maxwoffb: float = 80.0
    maxwoffr: float = 80.0
    secondcomp_mode: int = 0                        # 0 free, 1 larger vel, 2 larger sigma
    secondcomp_smart: bool = False
    instrres_mode: int = INSTRRES_VARPOLY
    instrres_tplsig: float = 0.0
    instrres_extpoly: np.ndarray = field(default_factory=lambda: np.zeros(7))
    instrres_varpoly: np.ndarray = field(default_factory=lambda: np.zeros(7))
    max_polycoeff_instrres: int = MAX_POLYCOEFF_INSTRRES
    fitallrange: np.ndarray = field(default_factory=lambda: np.zeros(4, dtype=int))
    Nloopadj: int = 0

    gnfit: np.ndarray = field(default_factory=_f25)      # user start values, narrow
    gbfit: np.ndarray = field(default_factory=_f25)      # user start values, broad
    pnfix: np.ndarray = field(default_factory=_b25)      # fix parameter, narrow
    pbfix: np.ndarray = field(default_factory=_b25)
    gnlims: np.ndarray = field(default_factory=_f25x2)   # user [min,max], narrow
    gblims: np.ndarray = field(default_factory=_f25x2)
    pndofit: np.ndarray = field(default_factory=_b25)    # fit narrow component of line i
    pbdofit: np.ndarray = field(default_factory=_b25)    # fit broad component of line i
    pcdofit: np.ndarray = field(default_factory=_b25)    # fit continuum for line i
    nshow: np.ndarray = field(default_factory=_b25)
    bshow: np.ndarray = field(default_factory=_b25)
    cshow: np.ndarray = field(default_factory=_b25)
    fitconstr: bool = False
    fitfixratios: bool = False

    Nlines_all: int = 0
    Nlines: int = 0
    lines_all: np.ndarray = field(default_factory=lambda: np.zeros(0))    # observed wavelengths, all lines in range
    lines: np.ndarray = field(default_factory=lambda: np.zeros(0))        # observed wavelengths, selected lines
    linesets: np.ndarray = field(default_factory=lambda: np.zeros(0, dtype=int))
    linenames: list = field(default_factory=list)
    linefancynames: list = field(default_factory=list)
    line_bisectors: np.ndarray = field(default_factory=lambda: np.zeros((0, 2)))
    lineset_ind: list = field(default_factory=list)     # per lineset: indices into ``lines`` (may be empty)
    lineset_max: int = 0
    redshift: float = 0.0
    selected_lineset: int = 0
    linefit_mode: int = MODE_SPAXEL
    linefit_type: int = FIT_GAUSS
    mom_windowmap: Optional[np.ndarray] = None          # (2, Nrow, Ncol)
    gauss_initmap: Optional[np.ndarray] = None          # (Nplanes, Nrow, Ncol)
    mom_windowmap_bootstrap: Optional[list] = None      # list of (2, Nrow, Ncol) or None
    mask_sn_thresh: float = 3.0
    mask_maxvelerr: float = 50.0
    mom_thresh: float = 0.5

    # ------------------------------------------------------------------ results
    results: dict = field(default_factory=dict)         # (mode, errmethod) -> ResultSet

    # GUI callbacks (set by the GUI, ignored headless)
    on_userpars_changed: Any = None

    # ================================================================== helpers
    @property
    def ckms(self) -> float:
        return CKMS

    @property
    def shape(self):
        return (self.Nwpix, self.Nrow, self.Ncol)

    def basename(self) -> str:
        """Input file name without the ``.fits`` extension (IDL strmid(...,-5))."""
        fn = self.filename
        return fn[:-5] if fn.lower().endswith(".fits") else fn

    # ------------------------------------------------------------------ results access
    def get_results(self, mode: int | None = None, errmethod: int | None = None) -> ResultSet:
        """Return (creating on demand) the ResultSet for ``(mode, errmethod)``."""
        mode = self.linefit_mode if mode is None else int(mode)
        errmethod = self.domontecarlo if errmethod is None else int(errmethod)
        key = (mode, errmethod)
        rs = self.results.get(key)
        if rs is None or rs.nlines != self.Nlines or (
                mode == MODE_SPAXEL and rs.n.shape[:2] != (self.Nrow, self.Ncol)):
            rs = ResultSet(mode, self.Nlines, ncol=self.Ncol, nrow=self.Nrow, nmask=self.Nmask)
            self.results[key] = rs
        if mode == MODE_MASK:
            rs.ensure_nmask(self.Nmask)
        return rs

    @property
    def res(self) -> ResultSet:
        return self.get_results()

    def reset_all_results(self) -> None:
        self.results = {}

    # ------------------------------------------------------------------ Monte Carlo cubes
    @property
    def montecarlocubes(self) -> Optional[np.ndarray]:
        """4D array (Nmc, Nw, Nrow, Ncol) of the realisations for the active error method."""
        return {ERR_BOOTSTRAP: self.bootstrapcubes, ERR_MC1: self.mc1cubes,
                ERR_MC2: self.mc2cubes, ERR_MC3: self.mc3cubes}.get(self.domontecarlo)

    @property
    def montecarlonoise_cube(self) -> Optional[np.ndarray]:
        return {ERR_BOOTSTRAP: self.bootstrapnoise, ERR_MC1: self.mc1noise,
                ERR_MC2: self.mc2noise, ERR_MC3: self.mc3noise}.get(self.domontecarlo)

    # ------------------------------------------------------------------ lines
    def mainline(self, redshift: float | None = None, vacuum: bool | None = None,
                 linesarr=None) -> float:
        """Observed wavelength of the main line (port of ``kubeviz_getmainline``)."""
        z = self.redshift if redshift is None else redshift
        vac = self.vacuum if vacuum is None else vacuum
        if linesarr is None:
            linesarr = self.linenames if self.Nlines > 0 else None
        rest = mainline_rest(linesarr) if (linesarr is not None and self.Nlines > 0) \
            else mainline_rest(None)
        if vac:
            rest = airtovac(rest)
        return rest * (1.0 + z)

    def mainline_index(self) -> int:
        """Index in ``lines`` of the main line (IDL ``where(abs(lines-mainline) lt 1)``)."""
        if self.Nlines == 0:
            return 0
        d = np.abs(self.lines - self.mainline())
        i = int(np.argmin(d))
        return i

    # ------------------------------------------------------------------ instrumental resolution
    def getinstrres(self, lam=None):
        """Resolution R = lambda/dlambda at ``lam`` (port of ``kubeviz_getinstrres``)."""
        if lam is None:
            lam = self.mainline()
        lam = np.asarray(lam, dtype=float)
        if self.Nlines <= 0:
            return np.zeros_like(lam) if lam.ndim else 0.0
        if self.instrres_mode == INSTRRES_VARPOLY:
            coeff = self.instrres_varpoly[: self.max_polycoeff_instrres + 1]
            return np.polynomial.polynomial.polyval(lam, coeff)
        if self.instrres_mode == INSTRRES_EXTPOLY:
            coeff = self.instrres_extpoly[: self.max_polycoeff_instrres + 1]
            return np.polynomial.polynomial.polyval(lam, coeff)
        if self.instrres_mode == INSTRRES_TEMPLATE:
            return lam / (SIGTOFWHM * self.instrres_tplsig)
        raise ValueError(f"unknown instrres_mode {self.instrres_mode}")

    def instrres_sigma_A(self, lam):
        """Instrumental Gaussian sigma in Angstrom at ``lam``."""
        R = self.getinstrres(lam)
        with np.errstate(divide="ignore", invalid="ignore"):
            return np.where(R > 0, lam / (SIGTOFWHM * R), 0.0)

    # ------------------------------------------------------------------ defaults (IDL kubeviz_set_*_defaults)
    def set_autoflag_defaults(self) -> None:
        self.mask_sn_thresh = 3.0
        self.mask_maxvelerr = 50.0

    def set_linefitrange_defaults(self) -> None:
        self.maxwoffb = 80.0
        self.maxwoffr = 80.0
        self.mom_thresh = 0.5

    def set_continuumfit_defaults(self) -> None:
        self.continuumfit_minoff = 200.0
        self.continuumfit_maxoff = 500.0
        self.continuumfit_minperc = 40.0
        self.continuumfit_maxperc = 60.0
        self.continuumfit_order = 0

    def set_user_fitpars(self, fitpars) -> None:
        """IDL ``kubeviz_set_user_fitpars``: 6 or 8 element vector."""
        fp = [float(v) for v in fitpars]
        if len(fp) == 6:
            (self.continuumfit_minoff, self.continuumfit_maxoff, self.continuumfit_minperc,
             self.continuumfit_maxperc, order, mode) = fp
        elif len(fp) == 8:
            (self.maxwoffb, self.maxwoffr, self.continuumfit_minoff, self.continuumfit_maxoff,
             self.continuumfit_minperc, self.continuumfit_maxperc, order, mode) = fp
        else:
            raise ValueError("fitpars must have 6 or 8 elements")
        self.continuumfit_order = int(order)
        self.continuumfit_mode = int(mode)

    # ------------------------------------------------------------------ misc
    def current_mask(self) -> np.ndarray:
        """2D mask (Nrow, Ncol) of the current spaxel mask."""
        return self.spaxselect[self.imask - 1]
