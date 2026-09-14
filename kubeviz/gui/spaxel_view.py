"""Spaxel viewer widget: the image with crosshair and mask overlay, colour bar,
wavelength slider and info panel (the ``base1`` window of the IDL code).

Interaction (modernised): mouse wheel zooms around the cursor, middle-drag (or
Shift + left drag) pans, "Reset view" restores the full field. Left click/drag
moves the crosshair or paints the spaxel mask, right-drag changes brightness and
contrast as in the IDL version.
"""
from __future__ import annotations

import numpy as np
import pyqtgraph as pg
from PyQt6.QtCore import Qt, pyqtSignal
from PyQt6.QtGui import QColor, QPainterPath, QPen
from PyQt6.QtWidgets import (QGridLayout, QHBoxLayout, QLabel, QPushButton, QSlider, QVBoxLayout,
                             QWidget)

pg.setConfigOptions(imageAxisOrder="row-major", background="w", foreground="k", antialias=False)


class SpaxelViewBox(pg.ViewBox):
    """ViewBox with kubeviz mouse semantics."""
    leftPressed = pyqtSignal(float, float)      # view coordinates
    leftDragged = pyqtSignal(float, float)
    leftReleased = pyqtSignal()
    rightDragged = pyqtSignal(float, float)     # fractional position in the widget (0..1)
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


class SpaxelView(QWidget):
    """Image + colour bar + slider + info panel."""
    spaxelPressed = pyqtSignal(int, int)
    spaxelDragged = pyqtSignal(int, int)
    spaxelReleased = pyqtSignal()
    contrastDragged = pyqtSignal(float, float)
    sliceChanged = pyqtSignal(int)
    keyPressed = pyqtSignal(object)
    resetViewRequested = pyqtSignal()

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setFocusPolicy(Qt.FocusPolicy.StrongFocus)
        lay = QVBoxLayout(self)
        lay.setContentsMargins(2, 2, 2, 2)

        self.glw = pg.GraphicsLayoutWidget()
        self.vb = SpaxelViewBox()
        self.glw.addItem(self.vb)
        self.image = pg.ImageItem()
        self.vb.addItem(self.image)
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

        self.vb.leftPressed.connect(lambda x, y: self.spaxelPressed.emit(int(np.floor(x)), int(np.floor(y))))
        self.vb.leftDragged.connect(lambda x, y: self.spaxelDragged.emit(int(np.floor(x)), int(np.floor(y))))
        self.vb.leftReleased.connect(self.spaxelReleased)
        self.vb.rightDragged.connect(self.contrastDragged)
        self.vb.keyPressed.connect(self.keyPressed)

        # colour bar
        cb = QHBoxLayout()
        self.colmin = QLabel(" ")
        self.colmin.setMinimumWidth(70)
        self.colmin.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.colorbar = pg.GraphicsLayoutWidget()
        self.colorbar.setFixedHeight(22)
        cvb = self.colorbar.addViewBox(enableMenu=False, enableMouse=False)
        cvb.setMouseEnabled(False, False)
        self.colorbar_item = pg.ImageItem()
        cvb.addItem(self.colorbar_item)
        cvb.setRange(xRange=(0, 256), yRange=(0, 1), padding=0)
        self.colmax = QLabel(" ")
        self.colmax.setMinimumWidth(70)
        cb.addWidget(self.colmin)
        cb.addWidget(self.colorbar, stretch=1)
        cb.addWidget(self.colmax)
        lay.addLayout(cb)

        # wavelength slider
        sl = QHBoxLayout()
        sl.addWidget(QLabel("Wave:"))
        self.slider = QSlider(Qt.Orientation.Horizontal)
        self.slider.setTracking(True)
        self.slider.valueChanged.connect(self.sliceChanged)
        sl.addWidget(self.slider, stretch=1)
        self.slice_label = QLabel("")
        self.slice_label.setMinimumWidth(140)
        sl.addWidget(self.slice_label)
        self.reset_btn = QPushButton("Reset view")
        self.reset_btn.setToolTip("Show the whole field (wheel: zoom, middle-drag / shift-drag: pan)")
        self.reset_btn.clicked.connect(self.resetViewRequested)
        sl.addWidget(self.reset_btn)
        lay.addLayout(sl)

        # info panel
        info = QGridLayout()
        info.setHorizontalSpacing(8)
        self.lbl_image = QLabel(" ")
        self.lbl_phys = QLabel(" ")
        self.lbl_value = QLabel(" ")
        self.lbl_mask = QLabel(" ")
        self.lbl_cube = QLabel(" ")
        self.lbl_cube.setFrameShape(QLabel.Shape.Panel)
        self.lbl_wcs = QLabel(" ")
        self.lbl_smooth = QLabel(" ")
        self.lbl_imgmode = QLabel(" ")
        self.lbl_imgmode.setFrameShape(QLabel.Shape.Panel)
        r = 0
        info.addWidget(QLabel("Image:"), r, 0)
        info.addWidget(self.lbl_image, r, 1)
        info.addWidget(QLabel("Phys:"), r, 2)
        info.addWidget(self.lbl_phys, r, 3)
        info.addWidget(self.lbl_value, r, 4)
        info.addWidget(QLabel("Mask:"), r, 5)
        info.addWidget(self.lbl_mask, r, 6)
        info.addWidget(self.lbl_cube, r, 7)
        r = 1
        info.addWidget(QLabel("WCS:"), r, 0)
        info.addWidget(self.lbl_wcs, r, 1, 1, 4)
        info.addWidget(QLabel("Smooth:"), r, 5)
        info.addWidget(self.lbl_smooth, r, 6)
        info.addWidget(self.lbl_imgmode, r, 7)
        lay.addLayout(info)

    # ------------------------------------------------------------------ API used by the controller
    def set_rgb(self, rgb: np.ndarray) -> None:
        """Display an (Nrow, Ncol, 3) uint8 image; pixel (col,row) covers [col,col+1]x[row,row+1]."""
        self.image.setImage(rgb, autoLevels=False)

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

    def set_colorbar(self, lut_rgb: np.ndarray, lo: float, hi: float) -> None:
        strip = (lut_rgb[None, :, :] * 255).astype(np.uint8)
        self.colorbar_item.setImage(strip, autoLevels=False)
        fmt = lambda v: f"{v:.1f}" if abs(v) < 1e5 else f"{v:.2E}"
        self.colmin.setText(fmt(lo))
        self.colmax.setText(fmt(hi))

    def set_slice(self, wpix: int, nwpix: int, wave: float | None) -> None:
        self.slider.blockSignals(True)
        self.slider.setRange(0, max(nwpix - 1, 0))
        self.slider.setValue(int(wpix))
        self.slider.blockSignals(False)
        txt = f"{wpix}" if wave is None else f"{wpix}   λ = {wave:.2f}"
        self.slice_label.setText(txt)

    def set_info(self, image_coords: str, phys_coords: str, value: float, mask: str, cube: str,
                 wcs: str, smooth: str, imgmode: str) -> None:
        self.lbl_image.setText(image_coords)
        self.lbl_phys.setText(phys_coords)
        self.lbl_value.setText("NaN" if not np.isfinite(value) else f"{value:.4f}")
        self.lbl_mask.setText(mask)
        self.lbl_cube.setText(cube)
        self.lbl_wcs.setText(wcs)
        self.lbl_smooth.setText(smooth)
        self.lbl_imgmode.setText(imgmode)

    def keyPressEvent(self, ev):
        self.keyPressed.emit(ev)
