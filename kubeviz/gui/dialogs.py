"""Small dialogs: smoothing, FITALL range, go-to-mask, mask parameters, zcut
parameters with histogram, FITS header viewer, help texts and the fit progress
dialog with its Interrupt button."""
from __future__ import annotations

import os

import numpy as np
import pyqtgraph as pg
from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import (QApplication, QDialog, QDialogButtonBox, QFormLayout, QHBoxLayout,
                             QLabel, QLineEdit, QMessageBox, QPlainTextEdit, QProgressBar,
                             QPushButton, QRadioButton, QSpinBox, QVBoxLayout)

from .. import __version__

DOC_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "doc")


class SmoothParsDialog(QDialog):
    def __init__(self, spat: int, spec: int, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Smooth Parameters")
        lay = QVBoxLayout(self)
        lay.addWidget(QLabel("WARNING: The current linefit results are reset\nif the smooth parameters are changed!"))
        form = QFormLayout()
        self.spat = QSpinBox()
        self.spat.setRange(1, 99)
        self.spat.setValue(spat)
        self.spec = QSpinBox()
        self.spec.setRange(1, 999)
        self.spec.setValue(spec)
        form.addRow("Spatial smoothing kernel:", self.spat)
        form.addRow("Spectral smoothing kernel:", self.spec)
        lay.addLayout(form)
        bb = QDialogButtonBox(QDialogButtonBox.StandardButton.Apply | QDialogButtonBox.StandardButton.Cancel)
        bb.button(QDialogButtonBox.StandardButton.Apply).clicked.connect(self.accept)
        bb.rejected.connect(self.reject)
        lay.addWidget(bb)

    def values(self):
        return self.spat.value(), self.spec.value()


class FitallRangeDialog(QDialog):
    def __init__(self, fitallrange, nloopadj: int, ncol: int, nrow: int, parent=None):
        super().__init__(parent)
        self.setWindowTitle("FITALL range")
        self.ncol, self.nrow = ncol, nrow
        lay = QVBoxLayout(self)
        lay.addWidget(QLabel("Parameters (pixel coordinates) for the FITALL button"))
        form = QFormLayout()
        self.xmin, self.xmax, self.ymin, self.ymax = (QSpinBox() for _ in range(4))
        for sb, mx in ((self.xmin, ncol - 1), (self.xmax, ncol - 1), (self.ymin, nrow - 1), (self.ymax, nrow - 1)):
            sb.setRange(0, max(mx, 0))
        self.xmin.setValue(int(fitallrange[0]))
        self.xmax.setValue(int(fitallrange[1]))
        self.ymin.setValue(int(fitallrange[2]))
        self.ymax.setValue(int(fitallrange[3]))
        self.nloop = QSpinBox()
        self.nloop.setRange(0, 100000)
        self.nloop.setValue(int(nloopadj))
        xr = QHBoxLayout()
        xr.addWidget(self.xmin)
        xr.addWidget(self.xmax)
        yr = QHBoxLayout()
        yr.addWidget(self.ymin)
        yr.addWidget(self.ymax)
        form.addRow("X range (min/max):", xr)
        form.addRow("Y range (min/max):", yr)
        form.addRow("Niter FITADJ (0 = until no improvement):", self.nloop)
        lay.addLayout(form)
        row = QHBoxLayout()
        apply = QPushButton("Apply")
        default = QPushButton("Defaults")
        cancel = QPushButton("Cancel")
        apply.clicked.connect(self._apply)
        default.clicked.connect(self._defaults)
        cancel.clicked.connect(self.reject)
        row.addWidget(apply)
        row.addWidget(default)
        row.addWidget(cancel)
        lay.addLayout(row)
        self.result_range = None
        self.result_nloop = None

    def _defaults(self):
        self.result_range = [0, self.ncol - 1, 0, self.nrow - 1]
        self.result_nloop = 0
        self.accept()

    def _apply(self):
        r = [self.xmin.value(), self.xmax.value(), self.ymin.value(), self.ymax.value()]
        if r[1] < r[0] or r[3] < r[2]:
            QMessageBox.warning(self, "FITALL range", "At least one value is out of limits")
            return
        self.result_range = r
        self.result_nloop = self.nloop.value()
        self.accept()


class GotoMaskDialog(QDialog):
    def __init__(self, imask: int, nmask: int, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Go To Mask")
        lay = QFormLayout(self)
        self.spin = QSpinBox()
        self.spin.setRange(1, max(nmask, 1))
        self.spin.setValue(imask)
        lay.addRow("Go To Mask:", self.spin)
        bb = QDialogButtonBox(QDialogButtonBox.StandardButton.Apply | QDialogButtonBox.StandardButton.Cancel)
        bb.button(QDialogButtonBox.StandardButton.Apply).clicked.connect(self.accept)
        bb.rejected.connect(self.reject)
        lay.addRow(bb)

    def value(self) -> int:
        return self.spin.value()


class MaskParsDialog(QDialog):
    def __init__(self, maskmode: int, maskradius: int, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Mask Parameters")
        lay = QVBoxLayout(self)
        row = QHBoxLayout()
        self.radios = [QRadioButton("SINGLE"), QRadioButton("CIRCULAR"), QRadioButton("SQUARE")]
        for r in self.radios:
            row.addWidget(r)
        self.radios[int(maskmode)].setChecked(True)
        lay.addLayout(row)
        form = QFormLayout()
        self.radius = QSpinBox()
        self.radius.setRange(1, 999)
        self.radius.setValue(int(maskradius))
        form.addRow("Mask Radius (Size):", self.radius)
        lay.addLayout(form)
        bb = QDialogButtonBox(QDialogButtonBox.StandardButton.Apply | QDialogButtonBox.StandardButton.Cancel)
        bb.button(QDialogButtonBox.StandardButton.Apply).clicked.connect(self.accept)
        bb.rejected.connect(self.reject)
        lay.addWidget(bb)

    def values(self):
        mode = [i for i, r in enumerate(self.radios) if r.isChecked()][0]
        return mode, self.radius.value()


class ZcutParsDialog(QDialog):
    """User min/max cuts with the pixel value histogram (``kubeviz_zcut_parameters``)."""

    def __init__(self, image, zmin: float, zmax: float, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Scale Parameters (zcut)")
        self.resize(460, 380)
        lay = QVBoxLayout(self)
        self.plot = pg.PlotWidget()
        self.plot.setLogMode(y=True)
        self.plot.setTitle("Pixel Distribution")
        lay.addWidget(self.plot, stretch=1)
        vals = np.asarray(image, dtype=float).ravel()
        vals = vals[np.isfinite(vals)]
        if vals.size > 1 and vals.max() > vals.min():
            nbins = max(int(vals.size / 3.0), 10)
            hist, edges = np.histogram(vals, bins=min(nbins, 400))
            self.plot.plot(edges, np.maximum(hist, 0.5), stepMode="center", fillLevel=0.1,
                           brush=(80, 80, 80, 120), pen=pg.mkPen("k"))
        self.lo_line = pg.InfiniteLine(pos=zmin, angle=90, movable=True, pen=pg.mkPen((220, 0, 0), width=2))
        self.hi_line = pg.InfiniteLine(pos=zmax, angle=90, movable=True, pen=pg.mkPen((0, 0, 220), width=2))
        self.plot.addItem(self.lo_line)
        self.plot.addItem(self.hi_line)
        form = QFormLayout()
        self.lo = QLineEdit(f"{zmin:.4g}")
        self.hi = QLineEdit(f"{zmax:.4g}")
        form.addRow("Limits   Low:", self.lo)
        form.addRow("High:", self.hi)
        lay.addLayout(form)
        self.lo.editingFinished.connect(lambda: self._edit_to_line(self.lo, self.lo_line))
        self.hi.editingFinished.connect(lambda: self._edit_to_line(self.hi, self.hi_line))
        self.lo_line.sigPositionChanged.connect(lambda: self.lo.setText(f"{self.lo_line.value():.4g}"))
        self.hi_line.sigPositionChanged.connect(lambda: self.hi.setText(f"{self.hi_line.value():.4g}"))
        bb = QDialogButtonBox(QDialogButtonBox.StandardButton.Apply | QDialogButtonBox.StandardButton.Cancel)
        bb.button(QDialogButtonBox.StandardButton.Apply).clicked.connect(self._apply)
        bb.rejected.connect(self.reject)
        lay.addWidget(bb)

    @staticmethod
    def _edit_to_line(edit, line):
        try:
            line.setValue(float(edit.text()))
        except ValueError:
            pass

    def _apply(self):
        lo, hi = self.values()
        if lo >= hi:
            QMessageBox.warning(self, "zcut", "Min value should be smaller than max value.")
            return
        self.accept()

    def values(self):
        return float(self.lo_line.value()), float(self.hi_line.value())


class TextDialog(QDialog):
    def __init__(self, title: str, text: str, parent=None, width=760, height=600):
        super().__init__(parent)
        self.setWindowTitle(title)
        self.resize(width, height)
        lay = QVBoxLayout(self)
        edit = QPlainTextEdit()
        edit.setReadOnly(True)
        edit.setLineWrapMode(QPlainTextEdit.LineWrapMode.NoWrap)
        font = edit.font()
        font.setFamily("Menlo")
        font.setStyleHint(font.StyleHint.Monospace)
        edit.setFont(font)
        edit.setPlainText(text)
        lay.addWidget(edit)
        bb = QDialogButtonBox(QDialogButtonBox.StandardButton.Close)
        bb.rejected.connect(self.reject)
        bb.clicked.connect(self.accept)
        lay.addWidget(bb)


def header_dialog(header, title: str, parent=None) -> TextDialog:
    text = header.tostring(sep="\n", endcard=False, padding=False) if header is not None else "(no header)"
    return TextDialog(title, text, parent)


def _doc_text(name: str) -> str:
    path = os.path.join(DOC_DIR, name)
    try:
        with open(path, encoding="utf-8", errors="replace") as fh:
            return f"Kubeviz version: {__version__}\n" + fh.read()
    except OSError:
        return f"Documentation file not found: {path}"


def help_whatsnew(parent=None) -> TextDialog:
    return TextDialog("What's new", _doc_text("whatsnew.txt"), parent)


def help_instructions(parent=None) -> TextDialog:
    return TextDialog("Instructions", _doc_text("instructions.txt"), parent)


SHORTCUTS_TEXT = f"""Kubeviz version: {__version__}

Keyboard shortcuts:

cursor keys : step through spatial plane
<  >  (, .) : step through wavelength planes
m           : select spaxel
n           : deselect spaxel
r           : clear all selected spaxels
k           : show previous spaxel mask
l           : show next spaxel mask
c           : fit a 2d gaussian to the current image
b           : toggle the flag for the selected spaxel
f           : Fit using the linefit framework (set parameters first!)
a           : Fit using linefit and the initial guess
              from adjacent spaxels results
q           : exit

In spectrum window only:

s           : Set wavelength (in combination with Sel Range button)
z           : Reset the redshift assuming the current wavelength
              to be that of the mainline

In spectral zoom window only:

g           : Simple Gauss fit at current wavelength

Mouse (spaxel viewer):

left click / drag   : move crosshair, paint or erase the spaxel mask
right drag          : brightness / contrast
wheel               : zoom in / out around the cursor
middle drag, shift+left drag : pan
"""


def help_shortcuts(parent=None) -> TextDialog:
    return TextDialog("Keyboard shortcuts", SHORTCUTS_TEXT, parent, width=640, height=560)


class FitProgressDialog(QDialog):
    """Progress bar with an Interrupt button (``kubeviz_linefit_cancel``)."""

    def __init__(self, title: str = "Fit", parent=None):
        super().__init__(parent)
        self.setWindowTitle(title)
        self.setModal(False)
        self.cancelled = False
        lay = QVBoxLayout(self)
        self.label = QLabel("Fitting...")
        lay.addWidget(self.label)
        self.bar = QProgressBar()
        self.bar.setRange(0, 1000)
        lay.addWidget(self.bar)
        btn = QPushButton(" Interrupt Fit ")
        btn.clicked.connect(self._cancel)
        lay.addWidget(btn)
        self.setMinimumWidth(320)

    def _cancel(self):
        self.cancelled = True

    def should_cancel(self) -> bool:
        QApplication.processEvents()
        return self.cancelled

    def progress(self, fraction: float, message: str) -> None:
        self.bar.setValue(int(1000 * min(max(fraction, 0.0), 1.0)))
        self.label.setText(message.replace("[PROGRES] ", ""))
        QApplication.processEvents()
