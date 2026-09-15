"""Main window and controller (ports of ``kubeviz_create``, ``kubeviz_setup_viewers``,
``kubeviz_spax_event``, ``kubeviz_spec_event``, ``kubeviz_speczoom_event``,
``kubeviz_linefit_event``, ``kubeviz_keyboard_handler``, ``kubeviz_plotspax``,
``kubeviz_plotinfo``, ``kubeviz_plotspec``, ``kubeviz_plotspeczoom``).

Layout: spaxel viewer on the left, spectrum (top) and spectral zoom (bottom) on the
right in resizable splitters; the linefit window is a separate top-level window.
"""
from __future__ import annotations

import os

import numpy as np
import logging

from PyQt6.QtCore import Qt, QTimer
from PyQt6.QtGui import QAction, QActionGroup
from PyQt6.QtWidgets import QApplication, QDockWidget, QFileDialog, QLabel, QMainWindow, QMessageBox

from .. import utils
from ..constants import (CUBE_BADPIX, CUBE_DATA, CUBE_LINEFIT, CUBE_LINEFIT_ERR, CUBE_LINEFIT_SN,
                         CUBE_NOISE, CUBE_SN, ERR_BOOTSTRAP, ERR_MC1, ERR_MC2, ERR_MC3, ERR_NOISE,
                         FIT_GAUSS, FIT_MOMENTS, IMG_MED1_MINUS_MED2, IMG_MED2_MINUS_MED1, IMG_MEDSUB1,
                         IMG_MEDSUB2, IMG_SLICE, MODE_MASK, MODE_SPAXEL, NOT_FIT, SPEC_MEDSUB,
                         SPEC_SLICE, SPEC_SUM)
from ..core import masks as masks_mod
from ..core.display import (CUBESEL_NAMES, IMGMODE_NAMES, current_image, current_value,
                            linefit_image_reset, linefit_image_update, toggle_flag)
from ..core.extraction import getspec, medianspec, medsum_image_update, spec_reset_range
from ..core.instrres import linefit_skylines
from ..core.montecarlo import (createmc1cubes, createmc2cubes, createmc3cubes, readbootstrapcubes,
                              setupmontecarlocubes)
from ..core.smooth import smooth_state_datacube, smooth_state_montecarlo
from ..fitting.fitall import fitadj, fitadjall, fitall
from ..fitting.flags import autoflag
from ..fitting.linefit import (change_redshift, chooselines, dofit, linefit_init, linefit_reset,
                               linefit_resetall, linefit_resetuser)
from ..fitting.mpfit import mpfitpeak
from ..io.results import loadres, saveres
from ..io.session import SESSION_EXT, load_session, save_session
from ..io.spectra import savecube, saveimage, savespec
from . import dialogs
from .fitoverlay import fit_overlays
from .linefit_window import LinefitPanels
from .qtenv import apply_app_font
from .scaling import (COLOUR_TABLES, ZCUT_HISTEQ, ZCUT_MINMAX, ZCUT_NAMES, ZCUT_USER_LIN, ZCUT_USER_LOG,
                      ZCUT_USER_SQRT, ZCUT_ZSCALE, ZCUT_950, ZCUT_970, ZCUT_990, ZCUT_995, colour_table,
                      contrast_lut_indices, scale_image)
import pyqtgraph as pg
from .spaxel_view import SpaxelView
from .spectrum_view import SpecZoomView, SpectrumView

UPDATE_NONE, UPDATE_FULL, UPDATE_FAST, UPDATE_RECOMPUTE = 0, 1, 2, 3


def open_cube_interactively(args: dict):
    """File dialog(s) then :func:`start_session` (``kubeviz_selectcube``)."""
    from astropy.io import fits

    from ..session import start_session
    fname, _ = QFileDialog.getOpenFileName(None, "Select a cube for reading", "", "FITS files (*.fits *.fits.gz *.fit);;All files (*)")
    if not fname:
        return None
    noisefile = None
    if not fname.endswith(SESSION_EXT):
        with fits.open(fname) as h:
            nextend = len(h) - 1
        if nextend < 2:
            noisefile, _ = QFileDialog.getOpenFileName(None, "Select the noise cube for reading", os.path.dirname(fname),
                                                       "FITS files (*.fits *.fits.gz *.fit);;All files (*)")
            if not noisefile:
                return None
    kwargs = {k: v for k, v in args.items() if k not in ("datafile", "batch", "fit_all_lines", "noisefile")}
    try:
        return start_session(fname, noisefile=noisefile, **kwargs)
    except Exception as exc:  # pragma: no cover - interactive path
        QMessageBox.critical(None, "kubeviz", f"Could not load {fname}:\n{exc}")
        return None


class KubevizGUI(QMainWindow):
    def __init__(self, state, args: dict | None = None):
        super().__init__()
        self.state = state
        self.args = dict(args or {})
        self.lut = colour_table(state.ctab)
        self.lut_indices = None
        self.dragging = 0
        self.last_range = (0.0, 1.0)
        self._cached_rgb = None
        self._nomove = False
        self._paint_value = 1

        apply_app_font()
        self.setWindowTitle(f"kubeviz: {state.filename}")
        self.setDockNestingEnabled(True)
        self.spax = SpaxelView()
        self.spec = SpectrumView()
        self.zoom = SpecZoomView()
        self.setCentralWidget(self.spax)

        def dock(title, widget, name):
            d = QDockWidget(title, self)
            d.setObjectName(name)
            d.setWidget(widget)
            d.setFeatures(QDockWidget.DockWidgetFeature.DockWidgetMovable | QDockWidget.DockWidgetFeature.DockWidgetFloatable
                          | QDockWidget.DockWidgetFeature.DockWidgetClosable)
            return d
        self.spec_dock = dock("Spectrum", self.spec, "dock_spectrum")
        self.zoom_dock = dock("Spectral zoom", self.zoom, "dock_zoom")
        self.ctrl_dock = dock("Line fitting", None, "dock_controls")
        area = Qt.DockWidgetArea.RightDockWidgetArea
        self.addDockWidget(area, self.spec_dock)
        self.addDockWidget(area, self.zoom_dock)
        self.addDockWidget(area, self.ctrl_dock)
        self.setCorner(Qt.Corner.BottomRightCorner, Qt.DockWidgetArea.RightDockWidgetArea)
        self.setCorner(Qt.Corner.TopRightCorner, Qt.DockWidgetArea.RightDockWidgetArea)
        self.zoom_dock.visibilityChanged.connect(self._zoom_dock_visibility)
        # fit the window to the screen: ~92 % of the available area, at most 1700x1050
        screen = QApplication.primaryScreen()
        avail = screen.availableGeometry() if screen is not None else None
        width = min(1700, int(avail.width() * 0.92)) if avail is not None else 1600
        height = min(1050, int(avail.height() * 0.92)) if avail is not None else 1000
        self.setMinimumSize(800, 500)
        self.resize(width, height)
        self._small_screen = height < 900
        self._docks_sized = False

        # status bar: last log message, progress and interrupt
        sb = self.statusBar()
        sb.setSizeGripEnabled(False)
        sb.setStyleSheet("QStatusBar { min-height: 20px; max-height: 22px; } QStatusBar::item { border: none; }")
        sb.setContentsMargins(4, 0, 4, 0)
        self.status_label = QLabel("")
        self.status_label.setContentsMargins(0, 0, 0, 0)
        sb.addWidget(self.status_label, 1)
        self._cancel_requested = False
        self._log_handler = _StatusLogHandler(self.status_label)
        utils.log.addHandler(self._log_handler)

        self._build_menus()
        self._connect()
        self.linefit = None
        self._install_linefit()
        self.spax.reset_view(state.Ncol, state.Nrow)
        self.zoom_dock.setVisible(bool(state.zoommap))
        self.update_all(UPDATE_FULL)

    def _apply_default_dock_sizes(self):
        """Right column ~52 % of the width; spectrum / zoom / line fitting 28 / 20 / 52 % of the height."""
        width, height = self.width(), self.height()
        self.resizeDocks([self.spec_dock, self.zoom_dock, self.ctrl_dock],
                         [int(height * 0.28), int(height * 0.20), int(height * 0.52)], Qt.Orientation.Vertical)
        self.resizeDocks([self.spec_dock, self.zoom_dock, self.ctrl_dock], [int(width * 0.52)] * 3, Qt.Orientation.Horizontal)

    def showEvent(self, ev):
        super().showEvent(ev)
        # the dock proportions are applied once the window is laid out at its real size
        if not self._docks_sized:
            self._docks_sized = True
            QTimer.singleShot(0, self._apply_default_dock_sizes)

    def _install_linefit(self):
        """(Re)create the line parameter table and the fit controls."""
        self.linefit = LinefitPanels(self)
        self.state.on_userpars_changed = lambda: self.linefit.update_all(update_userpars=True)
        old = self.ctrl_dock.widget()
        if old is not None:
            old.setParent(None)
            old.deleteLater()
        self.ctrl_dock.setWidget(self.linefit.controls)
        self.linefit.interrupt_btn.clicked.connect(self._request_cancel)
        if getattr(self, "_small_screen", False):
            self.linefit.setup_group.btn.setChecked(False)      # fold Fit setup on low displays

    def show_table(self):
        self.ctrl_dock.show()
        self.ctrl_dock.raise_()
        self.linefit.show_results()
        self.state.linefitmap = True

    def _zoom_dock_visibility(self, visible):
        self.state.zoommap = bool(visible)
        if visible:
            self.plotspeczoom()

    def _request_cancel(self):
        self._cancel_requested = True

    # ================================================================== window management
    def show_all(self):
        self.show()
        self.state.linefitmap = self.ctrl_dock.isVisible()

    def closeEvent(self, ev):
        # Qt 6 closes the top-level windows when the application quits, so this runs
        # both for the window close button and for File -> Quit; shut down only once
        self._shutdown()
        ev.accept()

    def _shutdown(self):
        if getattr(self, "_closed", False):
            return
        self._closed = True
        utils.info("Quitting...")
        utils.log.removeHandler(self._log_handler)

    def quit(self):
        self._shutdown()
        QApplication.instance().quit()

    # ================================================================== menus
    def _add_menu(self, menubar, title, items, group=None):
        m = menubar.addMenu(title)
        for it in items:
            if it is None:
                m.addSeparator()
                continue
            label, code = it
            a = QAction(label, self)
            a.triggered.connect(lambda _=False, c=code: self.menu_action(c))
            if group is not None:
                a.setCheckable(True)
                group.addAction(a)
            m.addAction(a)
            self.actions[code] = a
        return m

    def _build_menus(self):
        self.actions = {}
        mb = self.menuBar()
        self._add_menu(mb, "File", [("Open...", "Open"), ("Save Image as PNG...", "Print"), ("Save Image as FITS...", "SaveImage"),
                                    ("Save Cube as FITS...", "SaveCube"), ("Display FITS Header", "DisplayHeader"),
                                    ("Save Session...", "SaveSession"), ("Load Session...", "LoadSession"), None, ("Quit", "Quit")])
        self._add_menu(mb, "Spax.Masks", [("Select", "Select"), ("Deselect", "Deselect"), ("Clear", "Clear"), ("Optimal Mask", "OptimalMask"), None,
                                          ("New Mask", "NewMask"), ("Save...", "SAVE"), ("Load...", "Load"), ("Delete Mask", "DeleteMask"), None,
                                          ("Go to Mask...", "GoToMask"), ("Previous Mask", "PrevMask"), ("Next Mask", "NextMask"), None,
                                          ("Mask Parameters...", "MaskPars")])
        self.err_group = QActionGroup(self)
        self._add_menu(mb, "Errors", [("Use Noise-cube", "NoiseErrors"), ("Use Bootstraps", "BootstrapErrors"), ("Use Monte Carlo 1", "Mc1Errors"),
                                      ("Use Monte Carlo 2", "Mc2Errors"), ("Use Monte Carlo 3", "Mc3Errors")], self.err_group)
        self._add_menu(mb, "Options", [("Smooth parameters...", "SmoothPars"), ("FITALL range...", "FitallRange"),
                                       ("Use Montecarlo Noise (on/off)", "MonteCarloNoise"), ("Save Montecarlo Plots (on/off)", "MonteCarloPlot"),
                                       ("Save Montecarlo PDFs (on/off)", "MonteCarloSave"), ("Scale Noisecube error (on/off)", "NoiseCubeErrScale"),
                                       ("Load Results File...", "LoadResultFile"), None, ("Show linefit window", "ShowLinefit")])
        for code in ("MonteCarloNoise", "MonteCarloPlot", "MonteCarloSave", "NoiseCubeErrScale"):
            self.actions[code].setCheckable(True)
        self._add_menu(mb, "Help", [("What's new", "HelpWhatIsNew"), ("Instructions", "HelpInstructions"),
                                    ("Keyboard shortcuts", "HelpShortcuts"), ("Python port notes", "HelpPython")])
        self._sync_menu_checks()

    def _sync_menu_checks(self):
        st = self.state
        err = {0: "NoiseErrors", 1: "BootstrapErrors", 2: "Mc1Errors", 3: "Mc2Errors", 4: "Mc3Errors"}.get(st.domontecarlo)
        if err in self.actions:
            self.actions[err].setChecked(True)
        self.actions["MonteCarloNoise"].setChecked(bool(st.useMonteCarlonoise))
        self.actions["MonteCarloPlot"].setChecked(bool(st.plotMonteCarlodistrib))
        self.actions["MonteCarloSave"].setChecked(bool(st.saveMonteCarlodistrib))
        self.actions["NoiseCubeErrScale"].setChecked(bool(st.scaleNoiseerrors))
        self.spax.sync_controls(st.cubesel, st.imgmode, st.zcuts, st.ctab, st.invert == 1, st.cursormode)

    def _connect(self):
        self.spax.spaxelPressed.connect(self.on_spaxel_pressed)
        self.spax.spaxelDragged.connect(self.on_spaxel_dragged)
        self.spax.spaxelReleased.connect(self.on_spaxel_released)
        self.spax.contrastDragged.connect(self.on_contrast)
        self.spax.sliceChanged.connect(self.on_slice_changed)
        self.spax.keyPressed.connect(lambda ev: self.keyboard(ev, "spax"))
        self.spax.resetViewRequested.connect(lambda: self.spax.reset_view(self.state.Ncol, self.state.Nrow))
        self.spax.levelsDragged.connect(self._levels_dragged)
        cube_codes = {CUBE_DATA: "Data", CUBE_NOISE: "Noise", CUBE_BADPIX: "BadPixels", CUBE_SN: "SN",
                      CUBE_LINEFIT: "Linefit", CUBE_LINEFIT_ERR: "LineErrors", CUBE_LINEFIT_SN: "LineSN"}
        self.spax.cubeSelected.connect(lambda v: self.menu_action(cube_codes[v]))
        mode_codes = ["Slice", "Sum1", "Median1", "WeightedAvg1", "WeightedMed1", "MedSub1", "Sum2", "Median2",
                      "WeightedAvg2", "WeightedMed2", "MedSub2", "Med2-Med1", "Med1-Med2"]
        self.spax.imgModeSelected.connect(lambda v: self.menu_action(mode_codes[v]))
        zcut_codes = {ZCUT_HISTEQ: "HistEq", ZCUT_ZSCALE: "Zscale", ZCUT_MINMAX: "MinMax", ZCUT_995: "99.5", ZCUT_990: "99.0",
                      ZCUT_970: "97.0", ZCUT_950: "95.0", ZCUT_USER_LIN: "UserLin", ZCUT_USER_SQRT: "UserSqrt", ZCUT_USER_LOG: "UserLog"}
        self.spax.zcutSelected.connect(lambda v: self.menu_action(zcut_codes[v]))
        colour_codes = {num: name for name, num in COLOUR_TABLES}
        self.spax.colourSelected.connect(lambda v: self.menu_action(colour_codes[v]))
        self.spax.invertToggled.connect(lambda on: self.menu_action("Invert") if on != (self.state.invert == 1) else None)
        self.spax.cursorModeSelected.connect(lambda m: self.menu_action({0: "None", 1: "Crosshair", 2: "Select", 3: "Deselect"}[m]))
        self.spax.userCutsRequested.connect(lambda: self.menu_action("UserPars"))
        self.spec.wavelengthClicked.connect(self.on_spec_clicked)
        self.spec.keyPressed.connect(lambda ev: self.keyboard(ev, "spec"))
        self.spec.zminChanged.connect(self._set_zmin_spec)
        self.spec.zmaxChanged.connect(self._set_zmax_spec)
        self.spec.fixScaleToggled.connect(self._toggle_fixscale)
        self.spec.zoomToggled.connect(self._toggle_zoom)
        self.spec.saveRequested.connect(self._save_spec)
        self.spec.rangeSelected.connect(self._select_range)
        self.spec.modeSelected.connect(self._set_specmode)
        self.zoom.wavelengthClicked.connect(self.on_zoom_clicked)
        self.zoom.keyPressed.connect(lambda ev: self.keyboard(ev, "zoom"))
        self.zoom.zoomIn.connect(lambda: self._zoom_range(0.5))
        self.zoom.zoomOut.connect(lambda: self._zoom_range(2.0))

    # ================================================================== display
    def update_all(self, level: int = UPDATE_FULL, nomove: bool = False):
        """The ``update = 1/2/3`` block at the end of the IDL event handlers."""
        st = self.state
        if level == UPDATE_NONE:
            return
        if level == UPDATE_RECOMPUTE:
            if st.imgmode > 0:
                medsum_image_update(st)
                if st.imgmode in (IMG_MED2_MINUS_MED1, IMG_MED1_MINUS_MED2):
                    medsum_image_update(st, imgmode=IMG_MED2_MINUS_MED1 if st.imgmode == IMG_MED1_MINUS_MED2 else IMG_MED1_MINUS_MED2)
            if st.specmode > 0:
                medianspec(st)
            level = UPDATE_FULL
        self.linefit.update_all()
        self.plotspax(fast=(level == UPDATE_FAST))
        self.plotinfo()
        self.plotspec()
        self.plotspeczoom()
        self._sync_menu_checks()

    def plotspax(self, fast: bool = False, fitima_update: str | None = None):
        st = self.state
        if fitima_update is not None:
            st.cubesel = CUBE_LINEFIT
            par_image = fitima_update
        else:
            par_image = st.par_imagebutton
        if not fast or self._cached_rgb is None:
            if st.cubesel in (CUBE_LINEFIT, CUBE_LINEFIT_ERR, CUBE_LINEFIT_SN):
                ok = linefit_image_update(st, par_image)
                if not ok and not st.par_imagebutton:
                    linefit_image_reset(st)
            image, bad = current_image(st)
            image = np.array(image, dtype=float)
            image[bad] = np.nan
            scaled, rng = scale_image(image, st.zcuts, st.zmin_ima, st.zmax_ima)
            self.last_range = rng
            self._cached_scaled = scaled
            self._cached_image = image
            self._cached_rgb = True
            self._render()
        self.spax.set_crosshair(st.col, st.row, visible=st.cursormode >= 1)
        self.spax.set_mask(st.current_mask() if st.cursormode >= 2 else None, visible=st.cursormode >= 2)
        self.spax.set_slice(st.wpix, st.Nwpix, float(st.wave[st.wpix]) if st.wave is not None else None)

    def _lut_display(self):
        lut = self.lut[::-1] if self.state.invert == 1 else self.lut
        if self.lut_indices is not None:
            lut = lut[self.lut_indices]
        return lut

    _LINEAR_ZCUTS = (ZCUT_USER_LIN, ZCUT_MINMAX, ZCUT_ZSCALE, ZCUT_995, ZCUT_990, ZCUT_970, ZCUT_950)

    def _map_label(self) -> str:
        st = self.state
        if st.cubesel in (CUBE_LINEFIT, CUBE_LINEFIT_ERR, CUBE_LINEFIT_SN):
            pb = st.par_imagebutton
            if pb == "FLAG":
                what = "flag"
            elif pb == "CHISQ":
                what = "reduced χ²"
            elif pb:
                lt, par = pb[0], int(pb[1:])
                comp = {"N": "1st comp.", "B": "2nd comp.", "C": "continuum at"}[lt]
                if lt == "C":
                    what = f"{comp} {st.linefancynames[par - 3]}"
                elif par == 1:
                    what = f"{comp} velocity (km/s)"
                elif par == 2:
                    what = f"{comp} dispersion (km/s)"
                else:
                    what = f"{comp} {st.linefancynames[par - 3]} flux"
            else:
                what = "fit result"
            prefix = {CUBE_LINEFIT: "", CUBE_LINEFIT_ERR: "error of ", CUBE_LINEFIT_SN: "S/N of "}[st.cubesel]
            return prefix + what
        mode = IMGMODE_NAMES[st.imgmode] if st.imgmode != IMG_SLICE else f"slice {st.wpix}"
        return f"{CUBESEL_NAMES[st.cubesel]} · {mode}"

    def _render(self):
        """Send the cached image to the viewer with the current colour table and cuts."""
        st = self.state
        if self._cached_scaled is None:
            return
        lut = self._lut_display()
        cmap = pg.ColorMap(np.linspace(0, 1, lut.shape[0]), (lut * 255).astype(np.uint8))
        linear = st.zcuts in self._LINEAR_ZCUTS
        lo, hi = self.last_range
        label = self._map_label()
        if linear:
            self.spax.set_display(self._cached_image, (lo, hi), cmap, label, True)
        else:
            pos = np.array([0.0, 0.25, 0.5, 0.75, 1.0])
            if st.zcuts == ZCUT_HISTEQ:
                vals = np.nanquantile(self._cached_image, pos) if np.any(np.isfinite(self._cached_image)) else pos
            else:
                f, finv = (np.sqrt, np.square) if st.zcuts == ZCUT_USER_SQRT else (np.log10, lambda y: 10.0 ** y)
                lo_p = lo if lo > 0 else 1e-30
                flo, fhi = f(lo_p), f(max(hi, lo_p * (1 + 1e-9)))
                vals = finv(flo + pos * (fhi - flo))
            ticks = [(float(p), _tick_label(v)) for p, v in zip(pos, vals)]
            self.spax.set_display(self._cached_scaled, (0.0, 1.0), cmap, f"{label}  [{ZCUT_NAMES.get(st.zcuts, '')}]", False, ticks)

    def rerender(self):
        """Re-apply the colour table without recomputing the image (contrast drag)."""
        self._render()

    def _levels_dragged(self, lo, hi):
        st = self.state
        st.zmin_ima, st.zmax_ima = float(lo), float(hi)
        st.zcuts = ZCUT_USER_LIN
        self.plotspax()
        self._sync_menu_checks()

    def plotinfo(self):
        st = self.state
        val = current_value(st)
        wcs = "(WCS NOT FOUND)"
        try:
            from astropy.coordinates import Angle
            from astropy.wcs import WCS
            import astropy.units as u
            w = WCS(st.indatahead).celestial
            if w.has_celestial and st.wave is not None and st.wave[0] > 0:
                ra, dec = w.all_pix2world([[st.col + st.Startcol, st.row + st.Startrow]], 0)[0]
                ras = Angle(ra, u.deg).to_string(unit=u.hour, sep=":", precision=2, pad=True)
                decs = Angle(dec, u.deg).to_string(unit=u.deg, sep=":", precision=1, alwayssign=True, pad=True)
                wcs = f"<span style='color:#6b6b6b'>RA</span>&nbsp; {ras} &nbsp;&nbsp; <span style='color:#6b6b6b'>Dec</span>&nbsp; {decs}"
            else:
                wcs = "<span style='color:#6b6b6b'>No spatial WCS</span>"
        except Exception:
            wcs = "<span style='color:#6b6b6b'>No spatial WCS</span>"
        trimmed = bool(st.Startcol or st.Startrow)
        self.spax.set_info(st.col, st.row, st.col + st.Startcol, st.row + st.Startrow, val,
                           f"{st.imask}/{st.Nmask}", wcs, f"{st.smooth}x{st.smooth}x{st.specsmooth}", trimmed)

    def _showspec(self):
        st = self.state
        spec, nspec = getspec(st)
        zspec = st.badpixelmask[:, st.row, st.col].astype(float) if st.specmode in (SPEC_SLICE, SPEC_MEDSUB) else np.zeros(st.Nwpix)
        if st.cubesel == CUBE_NOISE:
            return nspec, None
        if st.cubesel == CUBE_BADPIX:
            return zspec, None
        if st.cubesel == CUBE_SN:
            return spec, nspec
        return spec, None

    def plotspec(self):
        st = self.state
        show, show2 = self._showspec()
        nw = st.Nwpix
        buf = 100 if nw > 400 else 1
        if st.scale == 1:
            yr = (st.zmin_spec, st.zmax_spec)
        else:
            sub = show[buf:nw - buf] if nw > 2 * buf else show
            with np.errstate(all="ignore"):
                yr = (np.nanmin(sub), np.nanmax(sub)) if np.any(np.isfinite(sub)) else (0.0, 1.0)
        titles = {0: f"Spaxel: {st.col},{st.row}", 1: "Sum of selected spaxels", 2: "Median of selected spaxels",
                  3: "Weigh.avg. of selected spaxels", 4: "Median subtracted", 5: "Robertson extracted"}
        title = titles[st.specmode] + f"   Slice: {st.wpix}"
        if st.wave[0] > 0:
            title += f"   λ: {st.wave[st.wpix]:.2f}"
        title += f"   Range: {st.wavsel}"
        w = st.wave
        r1 = (w[st.wavrange1[0]], w[st.wavrange1[1]]) if st.wavrange1[0] != st.wavrange1[1] else (None, None)
        r2 = (w[st.wavrange2[0]], w[st.wavrange2[1]]) if st.wavrange2[0] != st.wavrange2[1] else (None, None)
        start = None
        if st.npress == 1:
            start = w[st.wavrange1[0]] if st.wavsel == 1 else w[st.wavrange2[0]]
        overlays = fit_overlays(st, np.asarray(w, dtype=float)) if st.Nlines > 0 else ()
        self.spec.set_spectrum(w, show, show2, yrange=yr, title=title, marker_x=w[st.wpix] if st.marker == 1 else None,
                               ranges=(r1, r2), startmarker_x=start, overlays=overlays)
        self.spec.set_mode(st.specmode)
        self.spec.set_scale_controls(st.zmin_spec, st.zmax_spec, st.scale == 1)

    def plotspeczoom(self):
        st = self.state
        if not st.zoommap:
            return
        show, show2 = self._showspec()
        x1 = max(st.wpix - st.zoomrange, 0)
        x2 = min(st.wpix + st.zoomrange, st.Nwpix - 1)
        if st.scale == 1:
            yr = (st.zmin_spec, st.zmax_spec)
        else:
            sub = show[x1:x2 + 1]
            with np.errstate(all="ignore"):
                if np.any(np.isfinite(sub)):
                    yr = (min(0.0, 1.1 * np.nanmin(sub)), 1.1 * np.nanmax(sub))
                else:
                    yr = (0.0, 1.0)
        w = st.wave
        r1 = (w[st.wavrange1[0]], w[st.wavrange1[1]]) if st.wavrange1[0] != st.wavrange1[1] else (None, None)
        r2 = (w[st.wavrange2[0]], w[st.wavrange2[1]]) if st.wavrange2[0] != st.wavrange2[1] else (None, None)
        start = None
        if st.npress == 1:
            start = w[st.wavrange1[0]] if st.wavsel == 1 else w[st.wavrange2[0]]
        overlays = fit_overlays(st, np.asarray(w[x1:x2 + 1], dtype=float))
        self.zoom.set_spectrum(w, x1, x2, show, show2, yrange=yr, marker_x=w[st.wpix], ranges=(r1, r2),
                               startmarker_x=start, overlays=overlays, zoomrange=st.zoomrange)

    # ================================================================== spaxel viewer events
    def _clip_spaxel(self, col, row):
        st = self.state
        return int(np.clip(col, 0, st.Ncol - 1)), int(np.clip(row, 0, st.Nrow - 1))

    def on_spaxel_pressed(self, col, row):
        st = self.state
        st.col, st.row = self._clip_spaxel(col, row)
        self.dragging = 1
        if st.cursormode == 2:
            masks_mod.select_spaxels(st, st.col, st.row, 1)
        elif st.cursormode == 3:
            masks_mod.select_spaxels(st, st.col, st.row, 0)
        if st.cursormode < 2:
            self.update_all(UPDATE_FAST)
        else:
            self.plotspax(fast=True)
            self.plotinfo()

    def on_spaxel_dragged(self, col, row):
        if self.dragging != 1:
            return
        st = self.state
        c, r = self._clip_spaxel(col, row)
        if (c, r) == (st.col, st.row):
            return
        self.on_spaxel_pressed(c, r)

    def on_spaxel_released(self):
        st = self.state
        self.dragging = 0
        if st.cursormode in (2, 3):
            medianspec(st)
            self.linefit.update_all()
            self.plotspec()
            self.plotspeczoom()

    def on_contrast(self, x0, y0):
        self.lut_indices = contrast_lut_indices(x0, y0)
        self.rerender()

    def on_slice_changed(self, value):
        st = self.state
        if value == st.wpix:
            return
        st.wpix = int(value)
        if st.imgmode in (IMG_SLICE, IMG_MEDSUB1, IMG_MEDSUB2):
            self.update_all(UPDATE_FULL)
        else:
            self.update_all(UPDATE_FAST)

    # ================================================================== spectrum events
    def _set_wpix_from_wave(self, x):
        st = self.state
        pix = int(np.argmin(np.abs(st.wave - x)))
        st.wpix = int(np.clip(pix, 0, st.Nwpix - 1))
        if st.imgmode in (IMG_SLICE, IMG_MEDSUB1, IMG_MEDSUB2):
            self.update_all(UPDATE_FULL)
        else:
            self.update_all(UPDATE_FAST)

    def on_spec_clicked(self, x):
        self._set_wpix_from_wave(x)

    def on_zoom_clicked(self, x):
        self._set_wpix_from_wave(x)

    def _set_zmin_spec(self, v):
        self.state.zmin_spec = v
        self.update_all(UPDATE_FULL)

    def _set_zmax_spec(self, v):
        self.state.zmax_spec = v
        self.update_all(UPDATE_FULL)

    def _toggle_fixscale(self):
        self.state.scale = -self.state.scale
        self.update_all(UPDATE_FULL)

    def _toggle_zoom(self):
        self.zoom_dock.setVisible(not self.zoom_dock.isVisible())

    def _save_spec(self):
        fname, _ = QFileDialog.getSaveFileName(self, "Save spectrum", self.state.cwdir, "FITS (*.fits)")
        if fname:
            savespec(self.state, fname)

    def _select_range(self, which):
        st = self.state
        if which == 0:
            spec_reset_range(st)
        else:
            st.wavsel = which
        self.update_all(UPDATE_FULL)

    def _set_specmode(self, mode):
        st = self.state
        st.specmode = int(mode)
        st.linefit_mode = MODE_SPAXEL if mode in (SPEC_SLICE, SPEC_MEDSUB) else MODE_MASK
        if mode != SPEC_SLICE:
            medianspec(st)
        self.update_all(UPDATE_FULL)

    def _zoom_range(self, factor):
        st = self.state
        st.zoomrange = int(st.zoomrange * factor)
        st.zoomrange = int(np.clip(st.zoomrange, 1, max(st.Nwpix // 2, 1)))
        self.plotspeczoom()

    # ================================================================== keyboard
    def keyboard(self, ev, source: str):
        """Port of ``kubeviz_keyboard_handler`` plus the window-specific keys."""
        st = self.state
        key = ev.key()
        text = ev.text()
        update = UPDATE_NONE
        if key == Qt.Key.Key_Right:
            st.col = min(st.col + 1, st.Ncol - 1)
            update = UPDATE_FAST
        elif key == Qt.Key.Key_Left:
            st.col = max(st.col - 1, 0)
            update = UPDATE_FAST
        elif key == Qt.Key.Key_Up:
            st.row = min(st.row + 1, st.Nrow - 1)
            update = UPDATE_FAST
        elif key == Qt.Key.Key_Down:
            st.row = max(st.row - 1, 0)
            update = UPDATE_FAST
        elif text == "q":
            self.quit()
            return
        elif text in (",", "<"):
            st.wpix = max(st.wpix - 1, 0)
            update = UPDATE_FULL
        elif text in (".", ">"):
            st.wpix = min(st.wpix + 1, st.Nwpix - 1)
            update = UPDATE_FULL
        elif text == "m":
            st.spaxselect[st.imask - 1, st.row, st.col] = 1
            update = UPDATE_FULL
        elif text == "n":
            st.spaxselect[st.imask - 1, st.row, st.col] = 0
            update = UPDATE_FULL
        elif text == "r":
            masks_mod.clear_select(st)
            update = UPDATE_FULL
        elif text == "k":
            st.imask = max(st.imask - 1, 1)
            update = UPDATE_FULL
        elif text == "l":
            st.imask = min(st.imask + 1, st.Nmask)
            update = UPDATE_FULL
        elif text == "c":
            try:
                A = masks_mod.findcentroid(st)
                self.spax.set_ellipse(A[4], A[5], abs(A[2]), abs(A[3]), True)
            except Exception as exc:
                utils.warn(f"Centroid fit failed: {exc}")
        elif text == "f":
            self._fit_current()
            update = UPDATE_FULL
        elif text == "a":
            self._fitadj_current()
            update = UPDATE_FULL
        elif text == "b":
            toggle_flag(st, all_components=False)
            update = UPDATE_FULL
        elif text == "s" and source in ("spec", "zoom"):
            self._select_wave_range()
            update = UPDATE_FULL
        elif text == "z" and source in ("spec", "zoom"):
            newz = st.wave[st.wpix] / st.mainline(redshift=0.0) - 1.0
            self.do_change_redshift(newz)
            return
        elif text == "g" and source == "zoom":
            self._simple_gauss_fit()
            return
        if update != UPDATE_NONE:
            self.update_all(update)

    def _select_wave_range(self):
        st = self.state
        if st.wavsel == 0:
            st.wavsel = 1
        rng = st.wavrange1 if st.wavsel == 1 else st.wavrange2
        if st.npress == 0:
            rng[0] = rng[1] = st.wpix
            st.npress = 1
        else:
            rng[1] = st.wpix
            if st.imgmode > 0:
                medsum_image_update(st)
                if st.imgmode in (IMG_MED2_MINUS_MED1, IMG_MED1_MINUS_MED2):
                    medsum_image_update(st, imgmode=IMG_MED2_MINUS_MED1 if st.imgmode == IMG_MED1_MINUS_MED2 else IMG_MED1_MINUS_MED2)
            st.npress = 0

    def _simple_gauss_fit(self):
        st = self.state
        show, _ = self._showspec()
        binw = 40
        z0 = st.wpix
        a, b = max(z0 - binw, 0), min(z0 + binw, st.Nwpix - 1)
        xg = st.wave[a:b + 1]
        yg = np.asarray(show[a:b + 1], dtype=float)
        yfit, params, res = mpfitpeak(xg, yg, nterms=4, positive=False)
        if res is None:
            utils.warn("Gauss fit failed")
            return
        utils.info("Best-fit Gauss parameters: amplitude, centroid, width, y-offset:")
        utils.info(" ".join(f"{v:.4g}" for v in params))
        if st.cubesel > CUBE_NOISE:
            utils.warn("Fit applied to data, but alternative cube currently selected!")
        self.plotspeczoom()
        self.zoom.show_gauss_fit(xg, yfit)

    # ================================================================== fitting helpers
    def _autoflag_current(self):
        st = self.state
        if st.linefit_mode == MODE_SPAXEL:
            doflag = np.zeros((st.Nrow, st.Ncol), dtype=int)
            doflag[st.row, st.col] = 1
            autoflag(st, doflag=doflag)
        else:
            autoflag(st)

    def _fit_current(self):
        dofit(self.state)
        self._autoflag_current()

    def _fitadj_current(self):
        st = self.state
        if st.linefit_mode != MODE_SPAXEL:
            utils.warn("Not applicable to masks")
            return
        fitadj(st)
        self._autoflag_current()

    def _run_with_progress(self, title, func):
        """Run a long loop with a progress bar and an Interrupt button in the Status row
        of the line-fitting panel; the crosshair follows the spaxel being fitted."""
        import time
        st = self.state
        lf = self.linefit
        self._cancel_requested = False
        lf.progress_bar.setValue(0)
        lf.progress_bar.setFormat(f"{title} %p%")
        lf.progress_bar.setVisible(True)
        lf.interrupt_btn.setVisible(True)
        self.status_label.setText(f"{title}...")
        QApplication.processEvents()
        last_draw = [0.0]

        def should_cancel():
            QApplication.processEvents()
            return self._cancel_requested

        def progress(fraction, message):
            lf.progress_bar.setValue(int(1000 * min(max(fraction, 0.0), 1.0)))
            self.status_label.setText(message.replace("[PROGRES] ", ""))
            QApplication.processEvents()

        def on_fit(col, row):
            if col is not None:
                self.spax.set_crosshair(col, row, visible=True)
                now = time.time()
                if now - last_draw[0] > 0.25:          # redraw the map a few times per second at most
                    last_draw[0] = now
                    if st.cubesel > CUBE_SN:
                        self.plotspax(fast=True)
                    lf.update_all()
            QApplication.processEvents()

        try:
            result = func(should_cancel, progress, on_fit)
        finally:
            lf.progress_bar.setVisible(False)
            lf.interrupt_btn.setVisible(False)
            if self._cancel_requested:
                self.status_label.setText(f"{title} interrupted")
        return result

    def do_change_redshift(self, z):
        st = self.state
        if abs(z - st.redshift) < 1e-9:
            return
        change_redshift(st, z)
        self._install_linefit()
        self.update_all(UPDATE_FULL)

    # ================================================================== linefit window actions
    def linefit_action(self, code: str, value=None):
        """Port of ``kubeviz_linefit_event`` (button and text events)."""
        st = self.state
        lf = self.linefit
        # ------------------------------------------------------------ text fields
        if code == "REDSHIFT":
            if value is not None and value != NOT_FIT:
                self.do_change_redshift(float(value))
            return
        if code.startswith("SP") and code[2] in "NB":
            par = int(code[3:])
            (st.gnfit if code[2] == "N" else st.gbfit)[par] = value
            lf.update_all()
            return
        if code.startswith("MINP") or code.startswith("MAXP"):
            lt, par = code[4], int(code[5:])
            lims = st.gnlims if lt == "N" else st.gblims
            lims[par, 0 if code.startswith("MINP") else 1] = value
            lf.update_all()
            return
        simple = {"LINEFITMAXOFFB": "maxwoffb", "LINEFITMAXOFFR": "maxwoffr", "MASKSNTHRESH": "mask_sn_thresh",
                  "MASKMAXVELERR": "mask_maxvelerr", "MOMTHRESH": "mom_thresh", "CONTMINOFF": "continuumfit_minoff",
                  "CONTMAXOFF": "continuumfit_maxoff", "CONTMINPERC": "continuumfit_minperc",
                  "CONTMAXPERC": "continuumfit_maxperc"}
        if code in simple:
            if code != "MOMTHRESH" and code != "CONTMAXOFF" and value < 0:
                utils.warn("Value should be > 0!")
            elif code == "CONTMAXPERC" and value > 100:
                utils.warn("continuum max percentile should be <=100!")
            else:
                setattr(st, simple[code], float(value))
            lf.update_all()
            return
        if code == "CONTORDER":
            st.continuumfit_order = int(value)
            lf.update_all()
            return
        if code.startswith("POLYCOEFFPAR"):
            i = int(code[12:])
            if st.instrres_mode == 0:
                st.instrres_varpoly[i] = value
            elif st.instrres_mode == 1:
                st.instrres_extpoly[i] = value
            lf.update_all(update_userpars=True)
            self.plotspeczoom()
            return

        # ------------------------------------------------------------ buttons
        if code == "FLAG":
            toggle_flag(st, all_components=True)
            lf.update_all()
            if st.cubesel > CUBE_SN:
                self.plotspax()
                self.plotinfo()
        elif code == "FIT":
            self._fit_current()
            lf.update_all()
            self.plotspec()
            self.plotspeczoom()
            if st.cubesel > CUBE_SN:
                self.plotspax()
                self.plotinfo()
        elif code == "FITALL":
            def run(cancel, prog, on_fit):
                return fitall(st, spaxel=(st.linefit_mode == MODE_SPAXEL), should_cancel=cancel, on_progress=prog, on_fit=on_fit)
            self._run_with_progress("Fit all", run)
            autoflag(st)
            self.update_all(UPDATE_FULL)
        elif code == "FLAGALL":
            autoflag(st)
            lf.update_all()
            self.plotspax()
            self.plotinfo()
        elif code == "FITADJ":
            self._fitadj_current()
            lf.update_all()
            self.plotspeczoom()
            if st.cubesel > CUBE_SN:
                self.plotspax()
                self.plotinfo()
        elif code == "FITADJALL":
            if st.linefit_mode != MODE_SPAXEL:
                utils.warn("Not applicable to masks")
                return
            st.zoomspax = 1

            def run(cancel, prog, on_fit):
                return fitadjall(st, should_cancel=cancel, on_progress=prog, on_fit=on_fit)
            self._run_with_progress("Fit adjacent all", run)
            self.update_all(UPDATE_FULL)
        elif code == "RESIMAMASK":
            st.flagmode = not st.flagmode
            if st.flagmode:
                st.par_imagebutton = "FLAG"
            self.plotspax()
            self.plotinfo()
            lf.update_all()
        elif code == "SAVE":
            fname, _ = QFileDialog.getSaveFileName(self, "Save results", st.outdir, "FITS (*.fits)")
            if fname:
                saveres(st, fname)
        elif code == "RESETALL":
            linefit_reset(st)
            linefit_resetuser(st, all_pars=True)
            st.set_autoflag_defaults()
            st.set_linefitrange_defaults()
            st.set_continuumfit_defaults()
            lf.update_all(update_userpars=True)
            self.plotspeczoom()
        elif code == "RESETUSER":
            linefit_resetuser(st, all_pars=True)
            lf.update_all(update_userpars=True)
            self.plotspeczoom()
        elif code == "RESETFIT":
            linefit_reset(st)
            lf.update_all()
            self.plotspec()
            self.plotspeczoom()
        elif code == "RESETFITALL":
            linefit_resetall(st)
            if st.cubesel > CUBE_SN and (st.fitallrange[1] - st.fitallrange[0] + 1) == st.Ncol and (st.fitallrange[3] - st.fitallrange[2] + 1) == st.Nrow:
                st.cubesel = CUBE_DATA
            self.update_all(UPDATE_FULL)
        elif code == "MODE":
            if st.linefit_mode == MODE_SPAXEL:
                st.linefit_mode, st.specmode, st.cursormode = MODE_MASK, SPEC_SUM, 2
                medianspec(st)
            else:
                st.linefit_mode, st.specmode, st.cursormode = MODE_SPAXEL, SPEC_SLICE, 1
            self.update_all(UPDATE_FULL)
        elif code in ("PREVMASK", "NEXTMASK"):
            st.imask = max(st.imask - 1, 1) if code == "PREVMASK" else min(st.imask + 1, st.Nmask)
            medianspec(st)
            self.update_all(UPDATE_FULL)
        elif code == "FITSKY":
            st.instrres_mode = 0
            if np.sum(st.instrres_varpoly) == 0:
                linefit_skylines(st)
            lf.update_all(update_userpars=True)
            self.plotspeczoom()
        elif code == "POLYSKY":
            st.instrres_mode = 1
            lf.update_all(update_userpars=True)
            self.plotspeczoom()
        elif code == "TPLSKY":
            st.instrres_mode = 2
            lf.update_all()
            self.plotspeczoom()
        elif code == "TYPE":
            st.linefit_type = FIT_MOMENTS if st.linefit_type == FIT_GAUSS else FIT_GAUSS
            lf.typeswitch()
            self.update_all(UPDATE_FULL)
        elif code == "FITCONSTRAINTS":
            st.fitconstr = not st.fitconstr
        elif code == "FIXRATIOS":
            st.fitfixratios = not st.fitfixratios
        elif code == "CONTMODE":
            st.continuumfit_mode = 1 - st.continuumfit_mode
        elif code in ("FREE_2ND", "HIVEL_2ND", "BROAD_2ND"):
            st.secondcomp_mode = {"FREE_2ND": 0, "HIVEL_2ND": 1, "BROAD_2ND": 2}[code]
            lf.typeswitch()
        elif code == "SMART_2ND":
            st.secondcomp_smart = not st.secondcomp_smart
        elif code.startswith("FIX"):
            lt, par = code[3], int(code[4:])
            arr = st.pnfix if lt == "N" else st.pbfix
            arr[par] = 1 - arr[par]
        elif code.startswith("FIT"):
            lt, par = code[3], int(code[4:])
            line = par - 3
            if lt == "N":
                st.pndofit[line] = 1 - st.pndofit[line]
            elif lt == "B":
                st.pbdofit[line] = 1 - st.pbdofit[line]
            else:
                lines_c = chooselines(st, st.linesets, st.lines, lineset=int(st.linesets[line]))
                if lines_c.size:
                    utils.info("Setting continuum fitting for all lines in set.")
                    for il in lines_c:
                        st.pcdofit[il] = 1 - st.pcdofit[il]
                lf.update_all(update_userpars=True)
        elif code.startswith("RESET"):
            lt, par = code[5], int(code[6:])
            linefit_reset(st, linetypes=(lt,), pars=[par])
            lf.update_all()
            self.plotspeczoom()
        elif code.startswith("SHOW"):
            lt, par = code[4], int(code[5:])
            line = par - 3
            if lt == "N":
                st.nshow[line] = 1 - st.nshow[line]
            elif lt == "B":
                st.bshow[line] = 1 - st.bshow[line]
            else:
                lines_c = chooselines(st, st.linesets, st.lines, lineset=int(st.linesets[line]))
                for il in lines_c:
                    st.cshow[il] = 1 - st.cshow[il]
                lf.update_all(update_userpars=True)
            self.plotspec()
            self.plotspeczoom()
        elif code.startswith("IMAGE"):
            self.plotspax(fitima_update=code[5:])
            self.plotinfo()
            lf.update_all()
            self._sync_menu_checks()

    # ================================================================== menu actions (kubeviz_spax_event BUTT)
    def menu_action(self, code: str):
        st = self.state
        update = UPDATE_FULL
        if code == "Quit":
            self.quit()
            return
        if code == "Open":
            new = open_cube_interactively(self.args)
            if new is None:
                return
            self.replace_state(new)
            return
        if code == "Print":
            fname, _ = QFileDialog.getSaveFileName(self, "Save image", st.cwdir, "PNG (*.png)")
            if fname:
                self.spax.glw.grab().save(fname)
                utils.info(f"Created image file: {fname}")
            return
        if code == "SaveCube":
            fname, _ = QFileDialog.getSaveFileName(self, "Save cube", st.cwdir, "FITS (*.fits)")
            if fname:
                savecube(st, fname)
            return
        if code == "SaveImage":
            fname, _ = QFileDialog.getSaveFileName(self, "Save image", st.cwdir, "FITS (*.fits)")
            if fname:
                saveimage(st, fname)
            return
        if code == "DisplayHeader":
            items = [("Data Header", st.indatahead), ("Noise Header", st.innoisehead)]
            if st.inprihead is not None:
                items.insert(0, ("Primary Header", st.inprihead))
            from PyQt6.QtWidgets import QInputDialog
            choice, ok = QInputDialog.getItem(self, "Select header", "Header:", [n for n, _ in items], 0, False)
            if ok:
                hdr = dict(items)[choice]
                dialogs.header_dialog(hdr, choice, self).show()
            return
        if code == "SaveSession":
            fname, _ = QFileDialog.getSaveFileName(self, "Save session as", st.cwdir, f"kubeviz session (*{SESSION_EXT})")
            if fname:
                save_session(st, fname)
            return
        if code == "LoadSession":
            fname, _ = QFileDialog.getOpenFileName(self, "Load session file", st.cwdir,
                                                   f"kubeviz session (*{SESSION_EXT} *.sav);;All files (*)")
            if fname:
                try:
                    from ..session import finish_loaded_session
                    new = load_session(fname)
                    finish_loaded_session(new)
                    self.replace_state(new)
                except Exception as exc:
                    QMessageBox.critical(self, "kubeviz", f"Could not load session:\n{exc}")
            return
        if code == "ShowLinefit":
            self.show_table()
            return
        if code.startswith("Help"):
            {"HelpWhatIsNew": dialogs.help_whatsnew, "HelpInstructions": dialogs.help_instructions,
             "HelpShortcuts": dialogs.help_shortcuts, "HelpPython": dialogs.help_python}[code](self).show()
            return

        cubes = {"Data": CUBE_DATA, "Noise": CUBE_NOISE, "BadPixels": CUBE_BADPIX, "SN": CUBE_SN, "Linefit": CUBE_LINEFIT,
                 "LineErrors": CUBE_LINEFIT_ERR, "LineSN": CUBE_LINEFIT_SN}
        if code in cubes:
            st.cubesel = cubes[code]
            if st.cubesel > CUBE_SN and not st.par_imagebutton:
                st.par_imagebutton = "FLAG"
        elif code in dict(COLOUR_TABLES):
            st.ctab = dict(COLOUR_TABLES)[code]
            self.lut = colour_table(st.ctab)
            self.lut_indices = None
        elif code == "Invert":
            st.invert = -st.invert
            self.lut_indices = None
        elif code in ("MinMax", "Zscale", "HistEq", "99.5", "99.0", "97.0", "95.0", "UserLin", "UserSqrt", "UserLog"):
            st.zcuts = {"MinMax": ZCUT_MINMAX, "Zscale": ZCUT_ZSCALE, "HistEq": ZCUT_HISTEQ, "99.5": ZCUT_995, "99.0": ZCUT_990,
                        "97.0": ZCUT_970, "95.0": ZCUT_950, "UserLin": ZCUT_USER_LIN, "UserSqrt": ZCUT_USER_SQRT, "UserLog": ZCUT_USER_LOG}[code]
        elif code == "UserPars":
            img = self._cached_image if getattr(self, "_cached_image", None) is not None else current_image(st)[0]
            dlg = dialogs.ZcutParsDialog(img, st.zmin_ima, st.zmax_ima, self)
            if dlg.exec():
                st.zmin_ima, st.zmax_ima = dlg.values()
                st.zcuts = ZCUT_USER_LIN
            else:
                return
        elif code == "None":
            st.cursormode = 0
            update = UPDATE_FAST
        elif code == "Crosshair":
            st.cursormode, st.linefit_mode, st.specmode = 1, MODE_SPAXEL, SPEC_SLICE
            update = UPDATE_FAST
        elif code == "Select":
            st.cursormode, st.linefit_mode = 2, MODE_MASK
            if st.specmode == SPEC_SLICE:
                st.specmode = SPEC_SUM
                medianspec(st)
            update = UPDATE_FAST
        elif code == "Deselect":
            st.cursormode = 3
            update = UPDATE_FAST
        elif code == "Clear":
            masks_mod.clear_select(st)
            update = UPDATE_FAST
        elif code == "OptimalMask":
            if masks_mod.optimal_mask(st):
                if st.specmode == SPEC_SLICE:
                    st.specmode = SPEC_SUM
                medianspec(st)
            update = UPDATE_FAST
        elif code == "SAVE":
            fname, _ = QFileDialog.getSaveFileName(self, "Save mask", st.cwdir, "FITS (*.fits)")
            if fname:
                masks_mod.save_mask(st, fname)
            return
        elif code == "Load":
            fname, _ = QFileDialog.getOpenFileName(self, "Load mask", st.cwdir, "FITS (*.fits)")
            if not fname:
                return
            try:
                masks_mod.load_mask(st, fname)
            except Exception as exc:
                QMessageBox.critical(self, "kubeviz", str(exc))
                return
            update = UPDATE_FAST
        elif code == "NewMask":
            masks_mod.newmask(st)
            update = UPDATE_FAST
        elif code == "DeleteMask":
            masks_mod.deletemask(st)
            update = UPDATE_FAST
        elif code == "GoToMask":
            dlg = dialogs.GotoMaskDialog(st.imask, st.Nmask, self)
            if not dlg.exec():
                return
            st.imask = dlg.value()
            medianspec(st)
            update = UPDATE_FAST
        elif code == "PrevMask":
            st.imask = max(st.imask - 1, 1)
            medianspec(st)
            update = UPDATE_FAST
        elif code == "NextMask":
            st.imask = min(st.imask + 1, st.Nmask)
            medianspec(st)
            update = UPDATE_FAST
        elif code == "MaskPars":
            dlg = dialogs.MaskParsDialog(st.maskmode, st.maskradius, self)
            if dlg.exec():
                st.maskmode, st.maskradius = dlg.values()
            return
        elif code == "Slice":
            st.imgmode = IMG_SLICE
            spec_reset_range(st)
        elif code in IMGMODE_NAMES or code.replace("WeightedAvg", "Weighted Avg").replace("WeightedMed", "Weighted Med") in IMGMODE_NAMES:
            name = code.replace("WeightedAvg", "Weighted Avg").replace("WeightedMed", "Weighted Med")
            st.imgmode = IMGMODE_NAMES.index(name)
            if st.imgmode in (IMG_MED2_MINUS_MED1, IMG_MED1_MINUS_MED2):
                medsum_image_update(st, imgmode=IMG_MED2_MINUS_MED1)
                medsum_image_update(st, imgmode=IMG_MED1_MINUS_MED2)
            else:
                update = UPDATE_RECOMPUTE
        elif code == "SmoothPars":
            dlg = dialogs.SmoothParsDialog(st.smooth, st.specsmooth, self)
            if dlg.exec():
                st.smooth, st.specsmooth = dlg.values()
                smooth_state_datacube(st)
                smooth_state_montecarlo(st)
                linefit_init(st)
                medsum_image_update(st)
                self._install_linefit()
            else:
                return
        elif code == "FitallRange":
            dlg = dialogs.FitallRangeDialog(st.fitallrange, st.Nloopadj, st.Ncol, st.Nrow, self)
            if dlg.exec():
                st.fitallrange = np.array(dlg.result_range, dtype=int)
                st.Nloopadj = dlg.result_nloop
            return
        elif code in ("NoiseErrors", "BootstrapErrors", "Mc1Errors", "Mc2Errors", "Mc3Errors"):
            if not self._switch_error_method(code):
                return
            update = UPDATE_RECOMPUTE
        elif code == "MonteCarloNoise":
            st.useMonteCarlonoise = not st.useMonteCarlonoise
            utils.info(f"Montecarlo Noise switched {'ON' if st.useMonteCarlonoise else 'OFF'}")
            mcn = st.montecarlonoise_cube
            st.noise = mcn if (st.useMonteCarlonoise and mcn is not None) else st.noisecube
            update = UPDATE_RECOMPUTE
        elif code == "MonteCarloPlot":
            st.plotMonteCarlodistrib = not st.plotMonteCarlodistrib
            utils.info(f"Montecarlo plotting switched {'ON' if st.plotMonteCarlodistrib else 'OFF'}")
        elif code == "MonteCarloSave":
            st.saveMonteCarlodistrib = not st.saveMonteCarlodistrib
            utils.info(f"Montecarlo PDF saving switched {'ON' if st.saveMonteCarlodistrib else 'OFF'}")
        elif code == "NoiseCubeErrScale":
            st.scaleNoiseerrors = not st.scaleNoiseerrors
            utils.info(f"Noisecube errors scaling switched {'ON' if st.scaleNoiseerrors else 'OFF'}")
        elif code == "LoadResultFile":
            fname, _ = QFileDialog.getOpenFileName(self, "Load results file", st.cwdir, "FITS (*.fits)")
            if not fname:
                return
            loadres(st, fname)
            self.linefit.typeswitch()
        self.update_all(update)

    def _switch_error_method(self, code) -> bool:
        st = self.state
        target = {"NoiseErrors": ERR_NOISE, "BootstrapErrors": ERR_BOOTSTRAP, "Mc1Errors": ERR_MC1, "Mc2Errors": ERR_MC2, "Mc3Errors": ERR_MC3}[code]
        if st.domontecarlo == target:
            return False
        if target == ERR_NOISE:
            st.domontecarlo = ERR_NOISE
            st.noise = st.noisecube
        elif target == ERR_BOOTSTRAP:
            if st.Nbootstrap == 0:
                fname, _ = QFileDialog.getOpenFileName(self, "Select the bootstrap cubes", st.indir, "FITS (*.fits)")
                if not fname:
                    utils.warn("No bootstrap file selected")
                    return False
                st.domontecarlo = ERR_BOOTSTRAP
                st.bootstrap_file = fname
                readbootstrapcubes(st, os.path.basename(fname), os.path.dirname(fname))
                if st.domontecarlo != ERR_BOOTSTRAP:
                    return False
                setupmontecarlocubes(st)
            else:
                st.domontecarlo = ERR_BOOTSTRAP
                st.Nmontecarlo = st.Nbootstrap
                if st.useMonteCarlonoise:
                    st.noise = st.bootstrapnoise
        else:
            nattr = {ERR_MC1: "Nmc1", ERR_MC2: "Nmc2", ERR_MC3: "Nmc3"}[target]
            if getattr(st, nattr) == 0:
                if target == ERR_MC3 and st.redshift <= 0:
                    utils.warn("Set the correct redshift before calling this method.")
                    return False
                n = st.Nmc1 or st.Nmc2 or st.Nmc3 or 100
                setattr(st, nattr, n)
                st.domontecarlo = target
                QApplication.setOverrideCursor(Qt.CursorShape.WaitCursor)
                try:
                    {ERR_MC1: createmc1cubes, ERR_MC2: createmc2cubes, ERR_MC3: createmc3cubes}[target](st)
                    setupmontecarlocubes(st)
                finally:
                    QApplication.restoreOverrideCursor()
            else:
                st.domontecarlo = target
                st.Nmontecarlo = getattr(st, nattr)
                if st.useMonteCarlonoise:
                    st.noise = st.montecarlonoise_cube
        if st.cubesel > CUBE_SN:
            rs = st.get_results(MODE_SPAXEL, st.domontecarlo)
            if np.max(np.abs(rs.n[..., 1])) == 0.0:
                linefit_image_reset(st)
        return True

    def replace_state(self, new_state):
        """Swap in a new session (Open / Load Session), rebuilding the linefit window."""
        self.state.on_userpars_changed = None
        self.state = new_state
        self.state.on_userpars_changed = lambda: self.linefit.update_all(update_userpars=True)
        self.lut = colour_table(self.state.ctab)
        self.lut_indices = None
        self._cached_rgb = None
        self.setWindowTitle(f"kubeviz: {self.state.filename}")
        self._install_linefit()
        self.spax.reset_view(self.state.Ncol, self.state.Nrow)
        if self.state.specmode > 0:
            medianspec(self.state)
        self.update_all(UPDATE_FULL)


class _StatusLogHandler(logging.Handler):
    """Forward kubeviz log messages to the status bar label."""

    def __init__(self, label):
        super().__init__()
        self.label = label

    def emit(self, record):
        try:
            msg = record.getMessage()
            if msg.strip():
                self.label.setText(msg.replace("[KUBEVIZ] ", "").replace("[PROGRES] ", ""))
        except Exception:  # pragma: no cover
            pass


def _tick_label(v) -> str:
    """Short, readable number for the colour bar ticks."""
    v = float(v)
    if not np.isfinite(v):
        return ""
    if v == 0:
        return "0"
    if 1e-2 <= abs(v) < 1e5:
        return f"{v:.4g}" if abs(v) < 1000 else f"{v:.0f}"
    return f"{v:.2e}"
