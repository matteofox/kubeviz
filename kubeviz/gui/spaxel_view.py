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
from PyQt6.QtGui import QColor, QPainterPath, QPen
from PyQt6.QtWidgets import (QSizePolicy, QComboBox, QHBoxLayout, QLabel, QPushButton, QSlider, QToolButton,
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


SHORT_NAMES = {"Linefit errors": "Fit errors", "Linefit S/N": "Fit S/N", "Weighted Avg1": "W.Avg1",
               "Weighted Med1": "W.Med1", "Weighted Avg2": "W.Avg2", "Weighted Med2": "W.Med2",
               "STD GAMMA-II": "Gamma II", "User linear": "User lin", "User log10": "User log"}


def _combo(items, tip):
    """Toolbar combo that can shrink on small screens and grows with the available width."""
    c = QComboBox()
    for label, code in items:
        c.addItem(SHORT_NAMES.get(label, label), code)
    c.setToolTip(tip)
    c.setSizeAdjustPolicy(QComboBox.SizeAdjustPolicy.AdjustToMinimumContentsLengthWithIcon)
    c.setMinimumContentsLength(4)
    c.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
    c.setMaximumWidth(150)
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
    userCutsRequested = pyqtSignal()

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setFocusPolicy(Qt.FocusPolicy.StrongFocus)
        lay = QVBoxLayout(self)
        lay.setContentsMargins(4, 2, 4, 2)
        lay.setSpacing(3)

        # ---------------- toolbar
        tb = QHBoxLayout()
        tb.setSpacing(4)
        self.cube_combo = _combo(CUBE_ITEMS, "Which cube / map to display")
        self.mode_combo = _combo([(n, i) for i, n in enumerate(IMGMODE_NAMES)], "Image mode: single slice or a combination over a wavelength range")
        self.zcut_combo = _combo(ZCUT_ITEMS, "Intensity scaling")
        self.colour_combo = _combo(COLOUR_TABLES, "Colour table")
        self.invert_btn = QToolButton()
        self.invert_btn.setText("Inv")
        self.invert_btn.setCheckable(True)
        self.invert_btn.setToolTip("Invert the colour table")
        self.cuts_btn = QToolButton()
        self.cuts_btn.setText("Cuts")
        self.cuts_btn.setToolTip("Type the minimum / maximum for the user scalings (or drag the colour bar handles)")
        self.cuts_btn.clicked.connect(self.userCutsRequested)
        self.cursor_btns = []
        for text, tip, mode in (("+", "Crosshair: click / arrows select a spaxel", 1),
                                ("Sel", "Mask select: click / drag to add spaxels to the current mask", 2),
                                ("Desel", "Mask deselect: click / drag to remove spaxels from the mask", 3)):
            b = QToolButton()
            b.setText(text)
            b.setCheckable(True)
            b.setToolTip(tip)
            # clicking the active button switches the cursor off (crosshair) or back to the crosshair (mask modes)
            b.clicked.connect(lambda checked=False, m=mode: self.cursorModeSelected.emit(m if checked else (0 if m == 1 else 1)))
            self.cursor_btns.append(b)
        self.cursor_btns[0].setChecked(True)
        tb.addWidget(self.cube_combo, 1)
        tb.addWidget(self.mode_combo, 1)
        tb.addSpacing(6)
        tb.addWidget(self.zcut_combo, 1)
        tb.addWidget(self.cuts_btn)
        tb.addWidget(self.colour_combo, 1)
        tb.addWidget(self.invert_btn)
        tb.addSpacing(6)
        for b in self.cursor_btns:
            tb.addWidget(b)
        tb.addStretch(2)
        self.reset_btn = QPushButton("Reset")
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
        self.glw.ci.layout.setContentsMargins(0, 10, 0, 8)   # room for the colour bar tick labels
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

        # ---------------- info strip: fields spread along the row
        info_row = QHBoxLayout()
        info_row.setContentsMargins(4, 2, 4, 2)
        self.info_fields = {}
        for key in ("spaxel", "orig", "value", "wcs", "mask", "smooth"):
            lbl = QLabel("")
            lbl.setTextFormat(Qt.TextFormat.RichText)
            lbl.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
            self.info_fields[key] = lbl
            info_row.addWidget(lbl)
            info_row.addStretch(1)
        self.info_fields["orig"].setVisible(False)
        lay.addLayout(info_row)

    # ------------------------------------------------------------------ colour bar
    def _levels_finished(self, *args):
        if self._levels_from_code:
            return
        lo, hi = self.colorbar.levels()
        self.levelsDragged.emit(float(lo), float(hi))

    def set_display(self, data, levels, cmap: pg.ColorMap, label: str, linear: bool, ticks=None) -> None:
        """Show ``data`` (NaN = bad, drawn white) with ``levels``. The colour bar is
        interactive only for linear scalings; for the others ``ticks`` gives
        ``[(position, label), ...]`` with the real data values along the bar."""
        self._levels_from_code = True
        try:
            self.image.setImage(data, autoLevels=False)
            self.image.setLookupTable(cmap.getLookupTable(nPts=256))
            self.image.setLevels(levels)
            self.colorbar.setColorMap(cmap)
            self.colorbar.setLevels(levels)
            self.colorbar.interactive = linear
            self.colorbar.axis.setStyle(showValues=True, tickLength=5)
            self.colorbar.axis.setTicks(None if linear else [list(ticks or [])])
            self.colorbar.axis.setLabel(label, **{"font-size": "12pt"})
            self.colorbar.axis.setWidth(84)      # room for tick text and the label side by side
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
        txt = f"Slice {wpix}" if wave is None else f"Slice {wpix}   λ = {wave:.2f} Å"
        self.slice_label.setText(txt)

    @staticmethod
    def _field(name: str, value: str) -> str:
        return f"<span style='color:#6b6b6b'>{name}</span>&nbsp; {value}"

    def set_info(self, col, row, pcol, prow, value: float, mask: str, wcs: str, smooth: str,
                 trimmed: bool = False) -> None:
        """Fill the info strip; the original-cube position is shown only for a trimmed cube."""
        val = "NaN" if not np.isfinite(value) else f"{value:.4g}"
        f = self.info_fields
        f["spaxel"].setText(self._field("Spaxel", f"({col}, {row})"))
        f["orig"].setText(self._field("Original cube", f"({pcol}, {prow})"))
        f["orig"].setVisible(trimmed)
        f["value"].setText(self._field("Value", val))
        f["wcs"].setText(wcs)
        f["mask"].setText(self._field("Mask", mask))
        f["smooth"].setText(self._field("Smooth", smooth))

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
