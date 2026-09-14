"""Spaxel viewer panel: toolbar, image with an attached interactive colour bar,
wavelength slider and a one-line info strip.

Interaction: left click/drag moves the crosshair or paints the spaxel mask, right-drag
changes brightness/contrast, the wheel zooms around the cursor, middle-drag (or
shift + left drag) pans. Dragging the colour bar handles sets user linear cuts.
"""
from __future__ import annotations

import numpy as np
import pyqtgraph as pg
from PyQt6.QtCore import Qt, pyqtSignal
from PyQt6.QtGui import QColor, QFont, QPainterPath, QPen
from PyQt6.QtWidgets import (QComboBox, QHBoxLayout, QLabel, QPushButton, QSlider, QToolButton,
                             QVBoxLayout, QWidget)

from ..constants import (CUBE_BADPIX, CUBE_DATA, CUBE_LINEFIT, CUBE_LINEFIT_ERR, CUBE_LINEFIT_SN,
                         CUBE_NOISE, CUBE_SN)
from ..core.display import IMGMODE_NAMES
from .scaling import COLOUR_TABLES, ZCUT_NAMES

pg.setConfigOptions(imageAxisOrder="row-major", background="w", foreground="k", antialias=False)

CUBE_ITEMS = [("Data", CUBE_DATA), ("Noise", CUBE_NOISE), ("Bad pixels", CUBE_BADPIX), ("S/N", CUBE_SN),
              ("Fit map", CUBE_LINEFIT), ("Fit error map", CUBE_LINEFIT_ERR), ("Fit S/N map", CUBE_LINEFIT_SN)]
ZCUT_ITEMS = [(name, code) for code, name in ZCUT_NAMES.items()]


class SpaxelViewBox(pg.ViewBox):
    """ViewBox with kubeviz mouse semantics."""
    leftPressed = pyqtSignal(float, float)
    leftDragged = pyqtSignal(float, float)
    leftReleased = pyqtSignal()
    rightDragged = pyqtSignal(float, float)
    keyPressed = pyqtSignal(object)

    def __init__(self):
        super().__init__(lockAspect=True, enableMenu=False)
        self.setMouseMode(pg.ViewBox.PanMode)
        self.invertY(False)

    def mouseClickEvent(self, ev):
        if ev.button() == Qt.MouseButton.LeftButton:
            p = self.mapToView(ev.pos())
            self.leftPressed.emit(p.x(), p.y())
            self.leftReleased.emit()
            ev.accept()
        else:
            ev.ignore()

    def mouseDragEvent(self, ev, axis=None):
        mods = ev.modifiers()
        if ev.button() == Qt.MouseButton.MiddleButton or (
                ev.button() == Qt.MouseButton.LeftButton and mods & Qt.KeyboardModifier.ShiftModifier):
            ev.accept()
            delta = self.mapToView(ev.lastPos()) - self.mapToView(ev.pos())
            self.translateBy(x=delta.x(), y=delta.y())
            return
        if ev.button() == Qt.MouseButton.LeftButton:
            ev.accept()
            p = self.mapToView(ev.pos())
            if ev.isStart():
                self.leftPressed.emit(p.x(), p.y())
            else:
                self.leftDragged.emit(p.x(), p.y())
            if ev.isFinish():
                self.leftReleased.emit()
            return
        if ev.button() == Qt.MouseButton.RightButton:
            ev.accept()
            rect = self.boundingRect()
            if rect.width() > 0 and rect.height() > 0:
                x0 = (ev.pos().x() - rect.left()) / rect.width()
                y0 = 1.0 - (ev.pos().y() - rect.top()) / rect.height()
                self.rightDragged.emit(float(np.clip(x0, 0, 1)), float(np.clip(y0, 0, 1)))
            return
        ev.ignore()

    def keyPressEvent(self, ev):
        self.keyPressed.emit(ev)


def _combo(items, tip):
    c = QComboBox()
    for label, code in items:
        c.addItem(label, code)
    c.setToolTip(tip)
    c.setSizeAdjustPolicy(QComboBox.SizeAdjustPolicy.AdjustToContents)
    return c


class SpaxelView(QWidget):
    spaxelPressed = pyqtSignal(int, int)
    spaxelDragged = pyqtSignal(int, int)
    spaxelReleased = pyqtSignal()
    contrastDragged = pyqtSignal(float, float)
    sliceChanged = pyqtSignal(int)
    keyPressed = pyqtSignal(object)
    resetViewRequested = pyqtSignal()
    levelsDragged = pyqtSignal(float, float)          # user moved the colour bar handles
    cubeSelected = pyqtSignal(int)
    imgModeSelected = pyqtSignal(int)
    zcutSelected = pyqtSignal(int)
    colourSelected = pyqtSignal(int)
    invertToggled = pyqtSignal(bool)
    cursorModeSelected = pyqtSignal(int)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setFocusPolicy(Qt.FocusPolicy.StrongFocus)
        lay = QVBoxLayout(self)
        lay.setContentsMargins(4, 2, 4, 2)
        lay.setSpacing(3)

        # ---------------- toolbar
        tb = QHBoxLayout()
        tb.setSpacing(6)
        self.cube_combo = _combo(CUBE_ITEMS, "Which cube / map to display")
        self.mode_combo = _combo([(n, i) for i, n in enumerate(IMGMODE_NAMES)], "Image mode: single slice or a combination over a wavelength range")
        self.zcut_combo = _combo(ZCUT_ITEMS, "Intensity scaling")
        self.colour_combo = _combo(COLOUR_TABLES, "Colour table")
        self.invert_btn = QToolButton()
        self.invert_btn.setText("Inv")
        self.invert_btn.setCheckable(True)
        self.invert_btn.setToolTip("Invert the colour table")
        self.cursor_btns = []
        for text, tip, mode in (("+", "Crosshair: click / arrows select a spaxel", 1),
                                ("Sel", "Mask select: click / drag to add spaxels to the current mask", 2),
                                ("Desel", "Mask deselect: click / drag to remove spaxels from the mask", 3)):
            b = QToolButton()
            b.setText(text)
            b.setCheckable(True)
            b.setToolTip(tip)
            b.clicked.connect(lambda _=False, m=mode: self.cursorModeSelected.emit(m))
            self.cursor_btns.append(b)
        self.cursor_btns[0].setChecked(True)
        tb.addWidget(QLabel("Cube"))
        tb.addWidget(self.cube_combo)
        tb.addWidget(QLabel("Mode"))
        tb.addWidget(self.mode_combo)
        tb.addWidget(QLabel("Scale"))
        tb.addWidget(self.zcut_combo)
        tb.addWidget(self.colour_combo)
        tb.addWidget(self.invert_btn)
        tb.addSpacing(10)
        tb.addWidget(QLabel("Cursor"))
        for b in self.cursor_btns:
            tb.addWidget(b)
        tb.addStretch(1)
        self.reset_btn = QPushButton("Reset view")
        self.reset_btn.setToolTip("Show the whole field (wheel: zoom, middle-drag / shift-drag: pan)")
        self.reset_btn.clicked.connect(self.resetViewRequested)
        tb.addWidget(self.reset_btn)
        lay.addLayout(tb)
        self.cube_combo.currentIndexChanged.connect(lambda i: self.cubeSelected.emit(self.cube_combo.itemData(i)))
        self.mode_combo.currentIndexChanged.connect(lambda i: self.imgModeSelected.emit(self.mode_combo.itemData(i)))
        self.zcut_combo.currentIndexChanged.connect(lambda i: self.zcutSelected.emit(self.zcut_combo.itemData(i)))
        self.colour_combo.currentIndexChanged.connect(lambda i: self.colourSelected.emit(self.colour_combo.itemData(i)))
        self.invert_btn.toggled.connect(self.invertToggled)

        # ---------------- image + colour bar
        self.glw = pg.GraphicsLayoutWidget()
        self.glw.ci.layout.setContentsMargins(0, 0, 0, 0)
        self.vb = SpaxelViewBox()
        self.glw.addItem(self.vb, 0, 0)
        self.image = pg.ImageItem()
        self.vb.addItem(self.image)
        self.colorbar = pg.ColorBarItem(values=(0, 1), width=18, interactive=True, colorMapMenu=False,
                                        pen=pg.mkPen("k"), hoverPen=pg.mkPen("r"), rounding=1e-9)
        self.colorbar.setImageItem(self.image)
        self.colorbar.sigLevelsChangeFinished.connect(self._levels_finished)
        self.glw.addItem(self.colorbar, 0, 1)
        self.glw.ci.layout.setColumnStretchFactor(0, 1)
        self.vline = pg.InfiniteLine(angle=90, movable=False, pen=pg.mkPen((30, 80, 255), width=1.2))
        self.hline = pg.InfiniteLine(angle=0, movable=False, pen=pg.mkPen((30, 80, 255), width=1.2))
        self.vb.addItem(self.vline, ignoreBounds=True)
        self.vb.addItem(self.hline, ignoreBounds=True)
        self.mask_item = pg.QtWidgets.QGraphicsPathItem()
        self.mask_item.setPen(QPen(QColor(255, 40, 40), 0))
        self.vb.addItem(self.mask_item, ignoreBounds=True)
        self.ellipse_item = pg.QtWidgets.QGraphicsEllipseItem()
        self.ellipse_item.setPen(QPen(QColor(255, 0, 0), 0))
        self.ellipse_item.setVisible(False)
        self.vb.addItem(self.ellipse_item, ignoreBounds=True)
        lay.addWidget(self.glw, stretch=1)
        self._levels_from_code = False

        self.vb.leftPressed.connect(lambda x, y: self.spaxelPressed.emit(int(np.floor(x)), int(np.floor(y))))
        self.vb.leftDragged.connect(lambda x, y: self.spaxelDragged.emit(int(np.floor(x)), int(np.floor(y))))
        self.vb.leftReleased.connect(self.spaxelReleased)
        self.vb.rightDragged.connect(self.contrastDragged)
        self.vb.keyPressed.connect(self.keyPressed)

        # ---------------- wavelength slider
        sl = QHBoxLayout()
        sl.addWidget(QLabel("Wave"))
        self.slider = QSlider(Qt.Orientation.Horizontal)
        self.slider.setTracking(True)
        self.slider.valueChanged.connect(self.sliceChanged)
        sl.addWidget(self.slider, stretch=1)
        self.slice_label = QLabel("")
        self.slice_label.setMinimumWidth(170)
        sl.addWidget(self.slice_label)
        lay.addLayout(sl)

        # ---------------- info strip
        self.info = QLabel(" ")
        f = QFont("Menlo")
        f.setStyleHint(QFont.StyleHint.Monospace)
        f.setPointSize(11)
        self.info.setFont(f)
        self.info.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        lay.addWidget(self.info)

    # ------------------------------------------------------------------ colour bar
    def _levels_finished(self, *args):
        if self._levels_from_code:
            return
        lo, hi = self.colorbar.levels()
        self.levelsDragged.emit(float(lo), float(hi))

    def set_display(self, data, levels, cmap: pg.ColorMap, label: str, linear: bool) -> None:
        """Show ``data`` (NaN = bad, drawn white) with ``levels``; the colour bar is
        interactive only for linear scalings (its ticks are real values then)."""
        self._levels_from_code = True
        try:
            self.image.setImage(data, autoLevels=False)
            self.image.setLookupTable(cmap.getLookupTable(nPts=256))
            self.image.setLevels(levels)
            self.colorbar.setColorMap(cmap)
            self.colorbar.setLevels(levels)
            self.colorbar.interactive = linear
            self.colorbar.axis.setStyle(showValues=linear)
            self.colorbar.axis.setStyle(showValues=linear, tickLength=5 if linear else 0)
            self.colorbar.axis.setLabel(label)
        finally:
            self._levels_from_code = False

    # ------------------------------------------------------------------ overlays
    def reset_view(self, ncol: int, nrow: int) -> None:
        self.vb.setRange(xRange=(0, ncol), yRange=(0, nrow), padding=0.01)

    def set_crosshair(self, col: int, row: int, visible: bool = True) -> None:
        self.vline.setPos(col + 0.5)
        self.hline.setPos(row + 0.5)
        self.vline.setVisible(visible)
        self.hline.setVisible(visible)

    def set_mask(self, mask2d, visible: bool = True) -> None:
        path = QPainterPath()
        if visible and mask2d is not None:
            rows, cols = np.nonzero(np.asarray(mask2d) > 0)
            for r, c in zip(rows, cols):
                path.addRect(float(c), float(r), 1.0, 1.0)
                path.moveTo(float(c), float(r))
                path.lineTo(float(c) + 1.0, float(r) + 1.0)
        self.mask_item.setPath(path)

    def set_ellipse(self, cx, cy, rx, ry, visible=True) -> None:
        self.ellipse_item.setRect(cx - rx + 0.5, cy - ry + 0.5, 2 * rx, 2 * ry)
        self.ellipse_item.setVisible(visible)

    # ------------------------------------------------------------------ controls state
    def set_slice(self, wpix: int, nwpix: int, wave: float | None) -> None:
        self.slider.blockSignals(True)
        self.slider.setRange(0, max(nwpix - 1, 0))
        self.slider.setValue(int(wpix))
        self.slider.blockSignals(False)
        txt = f"slice {wpix}" if wave is None else f"slice {wpix}   λ = {wave:.2f} Å"
        self.slice_label.setText(txt)

    def set_info(self, col, row, pcol, prow, value: float, mask: str, cube: str, wcs: str,
                 smooth: str, imgmode: str) -> None:
        val = "NaN" if not np.isfinite(value) else f"{value:.4g}"
        self.info.setText(f"({col:>4},{row:>4})  phys ({pcol:>4},{prow:>4})  value {val:<12} {wcs}   mask {mask}  smooth {smooth}   {cube} · {imgmode}")

    def sync_controls(self, cubesel, imgmode, zcuts, ctab, invert, cursormode) -> None:
        for combo, value in ((self.cube_combo, cubesel), (self.mode_combo, imgmode),
                             (self.zcut_combo, zcuts), (self.colour_combo, ctab)):
            idx = combo.findData(value)
            if idx >= 0 and idx != combo.currentIndex():
                combo.blockSignals(True)
                combo.setCurrentIndex(idx)
                combo.blockSignals(False)
        self.invert_btn.blockSignals(True)
        self.invert_btn.setChecked(bool(invert))
        self.invert_btn.blockSignals(False)
        for i, b in enumerate(self.cursor_btns):
            b.setChecked(cursormode == i + 1)

    def keyPressEvent(self, ev):
        self.keyPressed.emit(ev)
