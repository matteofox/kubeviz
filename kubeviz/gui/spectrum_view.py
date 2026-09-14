"""Spectrum viewer and spectral zoom widgets (``base2`` and ``base3`` of the IDL code)."""
from __future__ import annotations

import numpy as np
import pyqtgraph as pg
from PyQt6.QtCore import Qt, pyqtSignal
from PyQt6.QtWidgets import (QComboBox, QHBoxLayout, QLabel, QLineEdit, QMenu, QPushButton,
                             QToolButton, QVBoxLayout, QWidget)

COL_DATA = (0, 0, 0)
COL_NOISE = (220, 0, 0)
COL_MARKER = (220, 0, 0)
COL_RANGE1 = (150, 150, 150, 90)
COL_RANGE2 = (215, 215, 215, 110)
COL_NARROW = (220, 0, 0)
COL_BROAD = (0, 0, 230)
COL_TOTAL = (0, 170, 0)
COL_CONT = (0, 0, 230)


class PixelAxis(pg.AxisItem):
    """Top axis of the zoom window labelled in pixel index (linear wavelength assumed)."""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.wave = None

    def set_wave(self, wave):
        self.wave = np.asarray(wave, dtype=float)

    def tickStrings(self, values, scale, spacing):
        if self.wave is None or self.wave.size < 2:
            return [f"{v:g}" for v in values]
        pix = np.interp(values, self.wave, np.arange(self.wave.size))
        return [f"{p:.0f}" for p in pix]


class _ClickablePlot(pg.PlotWidget):
    clickedAtX = pyqtSignal(float)
    keyPressed = pyqtSignal(object)

    def __init__(self, **kw):
        super().__init__(**kw)
        self.setFocusPolicy(Qt.FocusPolicy.StrongFocus)
        self.scene().sigMouseClicked.connect(self._clicked)

    def _clicked(self, ev):
        if ev.button() != Qt.MouseButton.LeftButton:
            return
        vb = self.plotItem.vb
        if vb.sceneBoundingRect().contains(ev.scenePos()):
            self.clickedAtX.emit(vb.mapSceneToView(ev.scenePos()).x())

    def keyPressEvent(self, ev):
        self.keyPressed.emit(ev)


class SpectrumView(QWidget):
    """Full spectrum with controls."""
    wavelengthClicked = pyqtSignal(float)
    keyPressed = pyqtSignal(object)
    zminChanged = pyqtSignal(float)
    zmaxChanged = pyqtSignal(float)
    fixScaleToggled = pyqtSignal()
    zoomToggled = pyqtSignal()
    saveRequested = pyqtSignal()
    rangeSelected = pyqtSignal(int)         # 1, 2 or 0 (reset)
    modeSelected = pyqtSignal(int)          # specmode

    def __init__(self, parent=None):
        super().__init__(parent)
        lay = QVBoxLayout(self)
        lay.setContentsMargins(2, 2, 2, 2)
        self.plot = _ClickablePlot()
        self.plot.setLabel("bottom", "Wavelength (Å)")
        self.plot.plotItem.setMenuEnabled(False)
        self.plot.plotItem.vb.setMouseEnabled(x=True, y=True)
        self.title = self.plot.plotItem.titleLabel
        self.region1 = pg.LinearRegionItem(movable=False, brush=pg.mkBrush(*COL_RANGE1), pen=pg.mkPen(None))
        self.region2 = pg.LinearRegionItem(movable=False, brush=pg.mkBrush(*COL_RANGE2), pen=pg.mkPen(None))
        for r in (self.region1, self.region2):
            r.setZValue(-10)
            r.setVisible(False)
            self.plot.addItem(r)
        self.curve = self.plot.plot([], [], pen=pg.mkPen(COL_DATA, width=1))
        self.curve2 = self.plot.plot([], [], pen=pg.mkPen(COL_NOISE, width=1))
        self.marker = pg.InfiniteLine(angle=90, pen=pg.mkPen(COL_MARKER, width=1))
        self.startmarker = pg.InfiniteLine(angle=90, pen=pg.mkPen(COL_DATA, width=1))
        self.startmarker.setVisible(False)
        self.plot.addItem(self.marker)
        self.plot.addItem(self.startmarker)
        self.plot.clickedAtX.connect(self.wavelengthClicked)
        self.plot.keyPressed.connect(self.keyPressed)
        self.overlay_items = []
        self._user_view = False          # True once the user zoomed/panned: keep that view
        self._wave_id = None
        self.plot.plotItem.vb.sigRangeChangedManually.connect(self._manual_range)
        lay.addWidget(self.plot, stretch=1)

        row = QHBoxLayout()
        self.auto_btn = QPushButton("Auto")
        self.auto_btn.setToolTip("Back to automatic axis ranges (they stay fixed once you zoom or pan)")
        self.auto_btn.clicked.connect(self.reset_view)
        row.addWidget(self.auto_btn)
        row.addWidget(QLabel("Min"))
        self.zmin_edit = QLineEdit("0.0")
        self.zmin_edit.setMaximumWidth(70)
        self.zmin_edit.editingFinished.connect(lambda: self._emit_float(self.zmin_edit, self.zminChanged))
        row.addWidget(self.zmin_edit)
        row.addWidget(QLabel("Max"))
        self.zmax_edit = QLineEdit("0.4")
        self.zmax_edit.setMaximumWidth(70)
        self.zmax_edit.editingFinished.connect(lambda: self._emit_float(self.zmax_edit, self.zmaxChanged))
        row.addWidget(self.zmax_edit)
        self.fix_btn = QPushButton("Fix Scale")
        self.fix_btn.setCheckable(True)
        self.fix_btn.clicked.connect(self.fixScaleToggled)
        row.addWidget(self.fix_btn)
        self.zoom_btn = QPushButton("Toggle Zoom")
        self.zoom_btn.clicked.connect(self.zoomToggled)
        row.addWidget(self.zoom_btn)
        save_btn = QPushButton("Save")
        save_btn.clicked.connect(self.saveRequested)
        row.addWidget(save_btn)

        self.range_btn = QToolButton()
        self.range_btn.setText("Sel Range")
        self.range_btn.setPopupMode(QToolButton.ToolButtonPopupMode.InstantPopup)
        m = QMenu(self.range_btn)
        m.addAction("Range1", lambda: self.rangeSelected.emit(1))
        m.addAction("Range2", lambda: self.rangeSelected.emit(2))
        m.addAction("Reset", lambda: self.rangeSelected.emit(0))
        self.range_btn.setMenu(m)
        row.addWidget(self.range_btn)

        row.addWidget(QLabel("Mode:"))
        self.mode_combo = QComboBox()
        self.mode_combo.addItems(["Spaxel", "Sum", "Median", "W. Avg", "MedSub", "Optimal"])
        self.mode_combo.currentIndexChanged.connect(self.modeSelected)
        row.addWidget(self.mode_combo)
        row.addStretch(1)
        lay.addLayout(row)

    @staticmethod
    def _emit_float(edit, signal):
        try:
            signal.emit(float(edit.text()))
        except ValueError:
            pass

    def set_mode(self, specmode: int) -> None:
        self.mode_combo.blockSignals(True)
        self.mode_combo.setCurrentIndex(int(specmode))
        self.mode_combo.blockSignals(False)

    def set_scale_controls(self, zmin, zmax, fixed: bool) -> None:
        self.zmin_edit.setText(f"{zmin:g}")
        self.zmax_edit.setText(f"{zmax:g}")
        self.fix_btn.setChecked(fixed)

    def _manual_range(self, *args):
        self._user_view = True

    def reset_view(self):
        self._user_view = False
        if self._last is not None:
            self.set_spectrum(**self._last)

    _last = None

    def set_spectrum(self, wave, spec, spec2=None, yrange=None, title="", marker_x=None,
                     ranges=((None, None), (None, None)), startmarker_x=None, overlays=()) -> None:
        self._last = dict(wave=wave, spec=spec, spec2=spec2, yrange=yrange, title=title, marker_x=marker_x,
                          ranges=ranges, startmarker_x=startmarker_x, overlays=overlays)
        self.curve.setData(wave, spec, connect="finite")
        if spec2 is not None:
            self.curve2.setData(wave, spec2, connect="finite")
            self.curve2.setVisible(True)
        else:
            self.curve2.setVisible(False)
        self.plot.setTitle(title, size="10pt")
        vb = self.plot.plotItem.vb
        new_wave = self._wave_id != (id(wave), len(wave))
        if new_wave:
            self._wave_id = (id(wave), len(wave))
            self._user_view = False
        if not self._user_view:
            vb.blockSignals(True)
            vb.setXRange(float(wave[0]), float(wave[-1]), padding=0)
            if yrange is not None and np.all(np.isfinite(yrange)) and yrange[1] > yrange[0]:
                vb.setYRange(float(yrange[0]), float(yrange[1]), padding=0.02)
            vb.blockSignals(False)
        for it in self.overlay_items:
            self.plot.removeItem(it)
        self.overlay_items = []
        for ov in overlays:
            if ov.kind == "contline":
                continue
            pen = {"narrow": pg.mkPen(COL_NARROW, width=1.5), "moment": pg.mkPen(COL_NARROW, width=1.5),
                   "broad": pg.mkPen(COL_BROAD, width=1.5),
                   "total": pg.mkPen(COL_TOTAL, width=1.5, style=Qt.PenStyle.DashLine)}[ov.kind]
            self.overlay_items.append(self.plot.plot(ov.x, ov.y, pen=pen, connect="finite"))
        for region, (lo, hi) in zip((self.region1, self.region2), ranges):
            if lo is not None and hi is not None and lo != hi:
                region.setRegion((min(lo, hi), max(lo, hi)))
                region.setVisible(True)
            else:
                region.setVisible(False)
        if marker_x is not None:
            self.marker.setPos(marker_x)
            self.marker.setVisible(True)
        else:
            self.marker.setVisible(False)
        if startmarker_x is not None:
            self.startmarker.setPos(startmarker_x)
            self.startmarker.setVisible(True)
        else:
            self.startmarker.setVisible(False)


class SpecZoomView(QWidget):
    """Zoom around the current wavelength with the fitted components overplotted."""
    wavelengthClicked = pyqtSignal(float)
    keyPressed = pyqtSignal(object)
    zoomIn = pyqtSignal()
    zoomOut = pyqtSignal()

    def __init__(self, parent=None):
        super().__init__(parent)
        lay = QVBoxLayout(self)
        lay.setContentsMargins(2, 2, 2, 2)
        self.top_axis = PixelAxis(orientation="top")
        self.plot = _ClickablePlot(axisItems={"top": self.top_axis})
        self.plot.plotItem.showAxis("top")
        self.plot.plotItem.setMenuEnabled(False)
        self.region1 = pg.LinearRegionItem(movable=False, brush=pg.mkBrush(*COL_RANGE1), pen=pg.mkPen(None))
        self.region2 = pg.LinearRegionItem(movable=False, brush=pg.mkBrush(*COL_RANGE2), pen=pg.mkPen(None))
        for r in (self.region1, self.region2):
            r.setZValue(-10)
            r.setVisible(False)
            self.plot.addItem(r)
        self.curve = self.plot.plot([], [], pen=pg.mkPen(COL_DATA, width=1), stepMode="center")
        self.curve2 = self.plot.plot([], [], pen=pg.mkPen(COL_NOISE, width=1))
        self.marker = pg.InfiniteLine(angle=90, pen=pg.mkPen(COL_MARKER, width=1, style=Qt.PenStyle.DashLine))
        self.startmarker = pg.InfiniteLine(angle=90, pen=pg.mkPen(COL_DATA, width=2))
        self.startmarker.setVisible(False)
        self.plot.addItem(self.marker)
        self.plot.addItem(self.startmarker)
        self.overlay_items = []
        self.gfit_curve = self.plot.plot([], [], pen=pg.mkPen((220, 0, 0), width=2))
        self.plot.clickedAtX.connect(self.wavelengthClicked)
        self.plot.keyPressed.connect(self.keyPressed)
        lay.addWidget(self.plot, stretch=1)
        row = QHBoxLayout()
        b_in = QPushButton("  +  ")
        b_out = QPushButton("  -  ")
        b_in.clicked.connect(self.zoomIn)
        b_out.clicked.connect(self.zoomOut)
        row.addWidget(b_in)
        row.addWidget(b_out)
        self.zoom_label = QLabel("")
        row.addWidget(self.zoom_label)
        row.addStretch(1)
        lay.addLayout(row)

    def set_spectrum(self, wave_full, x1: int, x2: int, spec, spec2=None, yrange=None, marker_x=None,
                     ranges=((None, None), (None, None)), startmarker_x=None, overlays=(), zoomrange=32) -> None:
        wave = np.asarray(wave_full, dtype=float)
        xx = wave[x1:x2 + 1]
        yy = np.asarray(spec, dtype=float)[x1:x2 + 1]
        self.top_axis.set_wave(wave)
        if xx.size >= 2:
            dx = np.diff(xx)
            edges = np.concatenate([[xx[0] - dx[0] / 2], (xx[:-1] + xx[1:]) / 2, [xx[-1] + dx[-1] / 2]])
            self.curve.setData(edges, np.nan_to_num(yy))
            self.plot.plotItem.vb.setXRange(float(edges[0]), float(edges[-1]), padding=0)
        else:
            self.curve.setData([], [])
        if spec2 is not None:
            self.curve2.setData(xx, np.asarray(spec2, dtype=float)[x1:x2 + 1], connect="finite")
            self.curve2.setVisible(True)
        else:
            self.curve2.setVisible(False)
        if yrange is not None and np.all(np.isfinite(yrange)) and yrange[1] > yrange[0]:
            self.plot.plotItem.vb.setYRange(float(yrange[0]), float(yrange[1]), padding=0.02)
        for region, (lo, hi) in zip((self.region1, self.region2), ranges):
            if lo is not None and hi is not None and lo != hi:
                region.setRegion((min(lo, hi), max(lo, hi)))
                region.setVisible(True)
            else:
                region.setVisible(False)
        if marker_x is not None:
            self.marker.setPos(marker_x)
        if startmarker_x is not None:
            self.startmarker.setPos(startmarker_x)
            self.startmarker.setVisible(True)
        else:
            self.startmarker.setVisible(False)
        for it in self.overlay_items:
            self.plot.removeItem(it)
        self.overlay_items = []
        for ov in overlays:
            if ov.kind == "contline":
                it = pg.InfiniteLine(pos=float(ov.x[0]), angle=90,
                                     pen=pg.mkPen(COL_CONT, width=1, style=Qt.PenStyle.DotLine))
                self.plot.addItem(it)
            else:
                pen = {"narrow": pg.mkPen(COL_NARROW, width=2), "moment": pg.mkPen(COL_NARROW, width=2),
                       "broad": pg.mkPen(COL_BROAD, width=2),
                       "total": pg.mkPen(COL_TOTAL, width=2, style=Qt.PenStyle.DashLine)}[ov.kind]
                it = self.plot.plot(ov.x, ov.y, pen=pen, connect="finite")
            self.overlay_items.append(it)
        self.gfit_curve.setVisible(False)
        self.zoom_label.setText(f"{2 * zoomrange:4d} pixel")

    def show_gauss_fit(self, x, y) -> None:
        self.gfit_curve.setData(x, y)
        self.gfit_curve.setVisible(True)
