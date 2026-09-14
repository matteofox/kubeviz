"""Line-fit panels (``kubeviz_setup_linefit`` / ``kubeviz_linefit_update`` /
``kubeviz_linefit_typeswitch``), redesigned as two dockable widgets:

* ``table``    - the parameter table: one row per line (1st component, with the
  continuum columns appended) plus optional 2nd-component rows and the two
  kinematic rows; it sits in the "Fit results" group at the top of ``controls``;
* ``controls`` - the panel: Fit results, Fit setup, Options and Status groups
  (each can be folded) and the grid of action buttons, always visible.

All user actions go to ``controller.linefit_action(code, value)`` with the IDL uvalue
codes (``'FIT'``, ``'SPN3'``, ``'MINPB1'``, ``'FIXN2'``, ``'IMAGEC4'`` ...).
"""
from __future__ import annotations

import numpy as np
from PyQt6.QtCore import Qt
from PyQt6.QtGui import QFont
from PyQt6.QtWidgets import (QButtonGroup, QCheckBox, QFrame, QGridLayout, QHBoxLayout, QLabel,
                             QLineEdit, QProgressBar, QPushButton, QRadioButton, QScrollArea,
                             QSizePolicy, QToolButton, QVBoxLayout, QWidget)

from .. import utils
from ..constants import (CKMS, ERR_METHOD_NAMES, FIT_GAUSS, INSTRRES_EXTPOLY, INSTRRES_TEMPLATE,
                         INSTRRES_VARPOLY, MODE_SPAXEL, NOT_FIT)

PANEL_STYLE = """
QLabel, QCheckBox, QRadioButton, QPushButton, QToolButton, QLineEdit, QComboBox { font-size: 12px; }
QLineEdit { padding: 1px 3px; }
QPushButton { padding: 3px 6px; }
QPushButton#action { min-height: 30px; font-size: 13px; }
QPushButton#save { min-height: 30px; font-size: 13px; font-weight: bold; }
"""

INSTRRES_MODE_TEXT = {INSTRRES_VARPOLY: "(polynomial fit to cube variance)",
                      INSTRRES_EXTPOLY: "(polynomial fit to external arcs)",
                      INSTRRES_TEMPLATE: "(spectral templates)"}


def _mono_font():
    f = QFont("Menlo")
    f.setStyleHint(QFont.StyleHint.Monospace)
    return f


def _userparstring(v) -> str:
    return "" if v == NOT_FIT else f"{v:.2f}"


def _edit(width=62, text="", tip=""):
    e = QLineEdit(text)
    e.setFixedWidth(width)
    e.setPlaceholderText("auto")
    if tip:
        e.setToolTip(tip)
    return e


def _mono_label(text="", width=None, align_right=True):
    lbl = QLabel(text)
    lbl.setFont(_mono_font())
    if width:
        lbl.setMinimumWidth(width)
    if align_right:
        lbl.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
    return lbl


def _cont_string(value, err):
    err = np.asarray(err, dtype=float)
    if err[0] == NOT_FIT:
        return "Not Fit"
    if err[0] == -998:
        return f"{value:.4f} +/- No Errors"
    return f"{value:.4f}+{err[0]:.4f}/{err[1]:.4f}"


class Collapsible(QWidget):
    """A titled section whose content can be folded away."""

    def __init__(self, title: str, content: QWidget, expanded: bool = True):
        super().__init__()
        lay = QVBoxLayout(self)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.setSpacing(0)
        self.btn = QToolButton()
        self.btn.setText(title)
        self.btn.setCheckable(True)
        self.btn.setChecked(expanded)
        self.btn.setToolButtonStyle(Qt.ToolButtonStyle.ToolButtonTextBesideIcon)
        self.btn.setArrowType(Qt.ArrowType.DownArrow if expanded else Qt.ArrowType.RightArrow)
        self.btn.setStyleSheet("QToolButton { border: none; font-weight: bold; text-align: left; padding: 2px; }")
        self.btn.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        self.content = content
        self.content.setVisible(expanded)
        self.btn.toggled.connect(self._toggle)
        lay.addWidget(self.btn)
        lay.addWidget(content)
        line = QFrame()
        line.setFrameShape(QFrame.Shape.HLine)
        lay.addWidget(line)

    def _toggle(self, on):
        self.content.setVisible(on)
        self.btn.setArrowType(Qt.ArrowType.DownArrow if on else Qt.ArrowType.RightArrow)


class LinefitPanels:
    """Builds ``self.table`` and ``self.controls`` for a controller/state pair."""

    COLS = ["", "comp.", "λ cen", "start", "best", "min", "max", "fix", "reset", "img", "fit", "show",
            "continuum", "img", "fit", "show"]
    _BANDS = ("#f7f7f7", "#e8e8e8")

    def __init__(self, controller):
        self.ctl = controller
        self.state = controller.state
        self.w = {}
        self._building = False
        self.image_group = QButtonGroup()
        self.image_group.setExclusive(True)
        self.table = self._build_table()
        self.controls = self._build_controls()
        self.controls.setStyleSheet(PANEL_STYLE)
        self.typeswitch()
        self.update_all(update_userpars=True)

    # ================================================================== table
    def _build_table(self) -> QWidget:
        self._building = True
        st = self.state
        outer_widget = QWidget()
        outer = QVBoxLayout(outer_widget)
        outer.setContentsMargins(4, 2, 4, 2)
        outer.setSpacing(2)

        hdr = QHBoxLayout()
        self.lbl_type = QLabel("GAUSS")
        self.lbl_type.setStyleSheet("font-weight: bold")
        hdr.addWidget(QLabel("Fit type"))
        hdr.addWidget(self.lbl_type)
        hdr.addSpacing(14)
        hdr.addWidget(QLabel("z ="))
        self.edit_z = _edit(90, f"{st.redshift:.6f}", "Redshift used to place the lines (Enter to apply; resets the fits)")
        self.edit_z.setPlaceholderText("")
        self.edit_z.returnPressed.connect(lambda: self._text("REDSHIFT", self.edit_z))
        hdr.addWidget(self.edit_z)
        hdr.addSpacing(14)
        self.lbl_mode = QLabel("SPAXEL")
        self.lbl_mode.setStyleSheet("font-weight: bold")
        hdr.addWidget(self.lbl_mode)
        self.lbl_sel = _mono_label("", align_right=False)
        hdr.addWidget(self.lbl_sel)
        self.btn_prevmask = QPushButton("-")
        self.btn_nextmask = QPushButton("+")
        for b, code in ((self.btn_prevmask, "PREVMASK"), (self.btn_nextmask, "NEXTMASK")):
            b.setFixedWidth(26)
            b.clicked.connect(lambda _=False, c=code: self.ctl.linefit_action(c))
            hdr.addWidget(b)
        hdr.addStretch(1)
        self.cb_second = QCheckBox("2nd component rows")
        self.cb_second.setToolTip("Show the rows of the second (broad / offset) component")
        self.cb_second.setChecked(bool(np.any(st.pbdofit[:st.Nlines])))
        self.cb_second.toggled.connect(self._toggle_second_rows)
        hdr.addWidget(self.cb_second)
        outer.addLayout(hdr)

        self.grid = QGridLayout()
        self.grid.setContentsMargins(2, 0, 2, 0)
        self.grid.setHorizontalSpacing(5)
        self.grid.setVerticalSpacing(0)
        self.col_labels = []
        for j, name in enumerate(self.COLS):
            lbl = QLabel(name)
            lbl.setAlignment(Qt.AlignmentFlag.AlignCenter)
            f = QFont()
            f.setBold(True)
            lbl.setFont(f)
            self.grid.addWidget(lbl, 0, j)
            self.col_labels.append(lbl)
        self.second_rows = []
        row = 1
        row = self._add_kin_rows(row, 1, "Line offset", "(km/s)", 0)
        row = self._add_kin_rows(row, 2, "Line width", "(km/s)", 1)
        for iline in range(st.Nlines):
            row = self._add_line_rows(row, iline + 3, iline)
        self.grid.setRowStretch(row, 1)
        self.grid.setColumnStretch(len(self.COLS), 1)
        table = QWidget()
        table.setLayout(self.grid)
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setWidget(table)
        scroll.setFrameShape(QFrame.Shape.NoFrame)
        scroll.setMinimumHeight(6 * 26 + 30)          # header + about six rows; grows with the dock
        scroll.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        self.table_scroll = scroll
        outer.addWidget(scroll, stretch=1)
        self._building = False
        self._toggle_second_rows(self.cb_second.isChecked())
        return outer_widget

    def _add_band(self, row, nrows, group):
        band = QFrame()
        band.setStyleSheet(f"background-color: {self._BANDS[group % 2]}; border-radius: 3px;")
        band.setAutoFillBackground(True)
        self.grid.addWidget(band, row, 0, nrows, len(self.COLS))
        band.lower()
        return band

    def _add_param_row(self, row, label, lt, par, comp_text, with_fit_show=True, cont_par=None):
        d = {}
        code = f"{lt}{par}"
        d["label"] = QLabel(label)
        d["label"].setStyleSheet("font-weight: bold")
        self.grid.addWidget(d["label"], row, 0)
        d["comp"] = QLabel(comp_text)
        self.grid.addWidget(d["comp"], row, 1)
        d["lamb"] = _mono_label("", 64)
        self.grid.addWidget(d["lamb"], row, 2)
        d["start"] = _edit(56, tip="Start value (empty = automatic guess)")
        d["start"].editingFinished.connect(lambda c="SP" + code, e=d["start"]: self._text(c, e))
        self.grid.addWidget(d["start"], row, 3)
        d["best"] = _mono_label("Not Fit", 150)
        self.grid.addWidget(d["best"], row, 4)
        d["min"] = _edit(52, tip="Lower limit (used when 'Fit with constraints' is on)")
        d["max"] = _edit(52, tip="Upper limit (used when 'Fit with constraints' is on)")
        d["min"].editingFinished.connect(lambda c="MINP" + code, e=d["min"]: self._text(c, e))
        d["max"].editingFinished.connect(lambda c="MAXP" + code, e=d["max"]: self._text(c, e))
        self.grid.addWidget(d["min"], row, 5)
        self.grid.addWidget(d["max"], row, 6)
        d["fix"] = QCheckBox()
        d["fix"].setToolTip("Keep this parameter fixed at its start value")
        d["fix"].clicked.connect(lambda _=False, c="FIX" + code: self.ctl.linefit_action(c))
        self.grid.addWidget(d["fix"], row, 7, alignment=Qt.AlignmentFlag.AlignCenter)
        d["reset"] = QToolButton()
        d["reset"].setText("×")
        d["reset"].setToolTip("Reset this parameter's result")
        d["reset"].clicked.connect(lambda _=False, c="RESET" + code: self.ctl.linefit_action(c))
        self.grid.addWidget(d["reset"], row, 8, alignment=Qt.AlignmentFlag.AlignCenter)
        d["image"] = QRadioButton()
        d["image"].setToolTip("Show this parameter as a map in the spaxel viewer")
        self.image_group.addButton(d["image"])
        d["image"].clicked.connect(lambda _=False, c="IMAGE" + code: self.ctl.linefit_action(c))
        self.grid.addWidget(d["image"], row, 9, alignment=Qt.AlignmentFlag.AlignCenter)
        if with_fit_show:
            d["fit"] = QCheckBox()
            d["fit"].setToolTip("Fit this component")
            d["fit"].clicked.connect(lambda _=False, c="FIT" + code: self.ctl.linefit_action(c))
            self.grid.addWidget(d["fit"], row, 10, alignment=Qt.AlignmentFlag.AlignCenter)
            d["show"] = QCheckBox()
            d["show"].setToolTip("Overplot this component on the spectrum")
            d["show"].clicked.connect(lambda _=False, c="SHOW" + code: self.ctl.linefit_action(c))
            self.grid.addWidget(d["show"], row, 11, alignment=Qt.AlignmentFlag.AlignCenter)
        self.w[(lt, par)] = d
        if cont_par is not None:
            c = {}
            ccode = f"C{cont_par}"
            c["best"] = _mono_label("Not Fit", 140)
            self.grid.addWidget(c["best"], row, 12)
            c["image"] = QRadioButton()
            c["image"].setToolTip("Show the continuum at this line as a map")
            self.image_group.addButton(c["image"])
            c["image"].clicked.connect(lambda _=False, cc="IMAGE" + ccode: self.ctl.linefit_action(cc))
            self.grid.addWidget(c["image"], row, 13, alignment=Qt.AlignmentFlag.AlignCenter)
            c["fit"] = QCheckBox()
            c["fit"].setToolTip("Fit the continuum for this lineset")
            c["fit"].clicked.connect(lambda _=False, cc="FIT" + ccode: self.ctl.linefit_action(cc))
            self.grid.addWidget(c["fit"], row, 14, alignment=Qt.AlignmentFlag.AlignCenter)
            c["show"] = QCheckBox()
            c["show"].setToolTip("Add the continuum to the overplotted components")
            c["show"].clicked.connect(lambda _=False, cc="SHOW" + ccode: self.ctl.linefit_action(cc))
            self.grid.addWidget(c["show"], row, 15, alignment=Qt.AlignmentFlag.AlignCenter)
            self.w[("C", cont_par)] = c
        return row + 1

    def _add_kin_rows(self, row, par, label, unit, group):
        band = self._add_band(row, 2, group)
        r0 = row
        row = self._add_param_row(row, label, "N", par, "1st", with_fit_show=False)
        row = self._add_param_row(row, unit, "B", par, "2nd", with_fit_show=False)
        self.second_rows.append(([self.w[("B", par)]], band, r0))
        return row

    def _add_line_rows(self, row, par, iline):
        band = self._add_band(row, 2, iline)
        r0 = row
        name = self.state.linefancynames[iline]
        row = self._add_param_row(row, name, "N", par, "1st", cont_par=par)
        row = self._add_param_row(row, "", "B", par, "2nd")
        self.second_rows.append(([self.w[("B", par)]], band, r0))
        return row

    def _toggle_second_rows(self, on: bool):
        for dicts, band, r0 in self.second_rows:
            for d in dicts:
                for wdg in d.values():
                    wdg.setVisible(on)
            self.grid.removeWidget(band)
            self.grid.addWidget(band, r0, 0, 2 if on else 1, len(self.COLS))
            band.lower()

    # ================================================================== controls
    def _build_controls(self) -> QWidget:
        self._building = True
        st = self.state
        panel = QWidget()
        lay = QVBoxLayout(panel)
        lay.setContentsMargins(4, 2, 4, 2)
        lay.setSpacing(2)

        # ---- Fit results (the table)
        self.results_group = Collapsible("Fit results", self.table, expanded=True)
        lay.addWidget(self.results_group, stretch=1)

        # ---- Fit setup
        setup = QWidget()
        g = QGridLayout(setup)
        g.setContentsMargins(6, 2, 6, 2)
        g.setHorizontalSpacing(8)
        g.setVerticalSpacing(2)
        g.addWidget(QLabel("Fit range blue / red (Å)"), 0, 0)
        self.edit_maxwoffb, self.edit_maxwoffr = _edit(62), _edit(62)
        rr = QHBoxLayout()
        rr.addWidget(self.edit_maxwoffb)
        rr.addWidget(self.edit_maxwoffr)
        rr.addStretch(1)
        g.addLayout(rr, 0, 1)
        self.cb_contmode = QCheckBox("Fit continuum with the lines (MPFIT CONT)")
        self.cb_contmode.setToolTip("Off: continuum from side bands, subtracted before the fit. "
                                    "On: a constant fitted together with the lines")
        self.cb_contmode.clicked.connect(lambda: self.ctl.linefit_action("CONTMODE"))
        g.addWidget(self.cb_contmode, 0, 2, 1, 2)
        g.addWidget(QLabel("Side bands offset min / max (Å)"), 1, 0)
        self.edit_cminoff, self.edit_cmaxoff = _edit(62), _edit(62)
        rr = QHBoxLayout()
        rr.addWidget(self.edit_cminoff)
        rr.addWidget(self.edit_cmaxoff)
        rr.addStretch(1)
        g.addLayout(rr, 1, 1)
        g.addWidget(QLabel("Percentiles min / max, order"), 1, 2)
        self.edit_cminperc, self.edit_cmaxperc = _edit(52), _edit(52)
        self.edit_corder = _edit(34)
        rr = QHBoxLayout()
        rr.addWidget(self.edit_cminperc)
        rr.addWidget(self.edit_cmaxperc)
        rr.addWidget(self.edit_corder)
        rr.addStretch(1)
        g.addLayout(rr, 1, 3)
        g.addWidget(QLabel("Autoflag S/N threshold"), 2, 0)
        self.edit_snthresh = _edit(62)
        g.addWidget(self.edit_snthresh, 2, 1)
        g.addWidget(QLabel("Max velocity / dispersion error (km/s)"), 2, 2)
        self.edit_maxvelerr = _edit(62)
        g.addWidget(self.edit_maxvelerr, 2, 3)
        g.setColumnStretch(4, 1)
        for e, code in ((self.edit_maxwoffb, "LINEFITMAXOFFB"), (self.edit_maxwoffr, "LINEFITMAXOFFR"),
                        (self.edit_cminoff, "CONTMINOFF"), (self.edit_cmaxoff, "CONTMAXOFF"),
                        (self.edit_cminperc, "CONTMINPERC"), (self.edit_cmaxperc, "CONTMAXPERC"),
                        (self.edit_corder, "CONTORDER"), (self.edit_snthresh, "MASKSNTHRESH"),
                        (self.edit_maxvelerr, "MASKMAXVELERR")):
            e.setPlaceholderText("")
            e.editingFinished.connect(lambda c=code, ee=e: self._text(c, ee))
        lay.addWidget(Collapsible("Fit setup", setup, expanded=True))

        # ---- Options
        opts = QWidget()
        g = QGridLayout(opts)
        g.setContentsMargins(6, 2, 6, 2)
        g.setVerticalSpacing(2)
        r1 = QHBoxLayout()
        self.cb_constr = QCheckBox("Fit with constraints")
        self.cb_constr.setToolTip("Use the user start values and min/max limits of the table")
        self.cb_constr.clicked.connect(lambda: self.ctl.linefit_action("FITCONSTRAINTS"))
        r1.addWidget(self.cb_constr)
        self.lbl_momthresh = QLabel("Moments threshold")
        self.edit_momthresh = _edit(62)
        self.edit_momthresh.setPlaceholderText("")
        self.edit_momthresh.editingFinished.connect(lambda: self._text("MOMTHRESH", self.edit_momthresh))
        r1.addWidget(self.lbl_momthresh)
        r1.addWidget(self.edit_momthresh)
        self.cb_fixratios = QCheckBox("Fix line ratios ([NII], [OIII], [OI])")
        self.cb_fixratios.clicked.connect(lambda: self.ctl.linefit_action("FIXRATIOS"))
        r1.addWidget(self.cb_fixratios)
        r1.addSpacing(12)
        self.lbl_second = QLabel("2nd component:")
        r1.addWidget(self.lbl_second)
        self.rb_second = [QRadioButton("fainter"), QRadioButton("larger offset"), QRadioButton("larger width")]
        self.second_group = QButtonGroup()
        for i, (rb, code) in enumerate(zip(self.rb_second, ("FREE_2ND", "HIVEL_2ND", "BROAD_2ND"))):
            self.second_group.addButton(rb, i)
            rb.clicked.connect(lambda _=False, c=code: self.ctl.linefit_action(c))
            r1.addWidget(rb)
        self.cb_smart = QCheckBox("smart")
        self.cb_smart.setToolTip("Keep the 2nd component only when two kinematic components are really present")
        self.cb_smart.clicked.connect(lambda: self.ctl.linefit_action("SMART_2ND"))
        r1.addWidget(self.cb_smart)
        r1.addStretch(1)
        g.addLayout(r1, 0, 0)
        r2 = QHBoxLayout()
        r2.setSpacing(4)
        r2.addWidget(QLabel("Instr. resolution:"))
        b = QPushButton("Fit sky lines")
        b.setToolTip("Fit the sky lines in the variance cube")
        b.clicked.connect(lambda: self.ctl.linefit_action("FITSKY"))
        r2.addWidget(b)
        self.btn_poly = QPushButton("Use polynomial")
        self.btn_poly.setToolTip("Polynomial from the header / archive / instrument manual")
        self.btn_poly.clicked.connect(lambda: self.ctl.linefit_action("POLYSKY"))
        r2.addWidget(self.btn_poly)
        self.btn_tpl = QPushButton("Use templates")
        self.btn_tpl.clicked.connect(lambda: self.ctl.linefit_action("TPLSKY"))
        r2.addWidget(self.btn_tpl)
        r2.addSpacing(8)
        r2.addWidget(QLabel("R(main line) ="))
        self.lbl_instrres = _mono_label("", 70, align_right=False)
        r2.addWidget(self.lbl_instrres)
        self.lbl_instrres_mode = QLabel("")
        r2.addWidget(self.lbl_instrres_mode)
        r2.addStretch(1)
        g.addLayout(r2, 1, 0)
        r3 = QHBoxLayout()
        r3.setSpacing(4)
        r3.addWidget(QLabel("Polynomial coefficients:"))
        self.edit_coeff = []
        for i in range(st.max_polycoeff_instrres + 1):
            e = _edit(90)
            e.setPlaceholderText("")
            e.editingFinished.connect(lambda i=i, ee=e: self._text(f"POLYCOEFFPAR{i}", ee))
            self.edit_coeff.append(e)
            r3.addWidget(e)
        r3.addStretch(1)
        g.addLayout(r3, 2, 0)
        lay.addWidget(Collapsible("Options", opts, expanded=False))

        # ---- Status
        status = QWidget()
        g = QHBoxLayout(status)
        g.setContentsMargins(6, 2, 6, 2)
        g.addWidget(QLabel("Errors:"))
        self.lbl_errmethod = QLabel("")
        g.addWidget(self.lbl_errmethod)
        g.addSpacing(10)
        g.addWidget(QLabel("χ²/dof"))
        self.lbl_chisq = _mono_label("0.0000", 70, align_right=False)
        g.addWidget(self.lbl_chisq)
        self.rb_chisq = QRadioButton("map")
        self.image_group.addButton(self.rb_chisq)
        self.rb_chisq.clicked.connect(lambda: self.ctl.linefit_action("IMAGECHISQ"))
        g.addWidget(self.rb_chisq)
        g.addSpacing(10)
        g.addWidget(QLabel("Flag"))
        self.btn_flag = QPushButton("OK")
        self.btn_flag.setFixedWidth(54)
        self.btn_flag.setToolTip("Toggle the fit quality flag of this spaxel / mask")
        self.btn_flag.clicked.connect(lambda: self.ctl.linefit_action("FLAG"))
        g.addWidget(self.btn_flag)
        self.rb_flag = QRadioButton("map")
        self.image_group.addButton(self.rb_flag)
        self.rb_flag.clicked.connect(lambda: self.ctl.linefit_action("IMAGEFLAG"))
        g.addWidget(self.rb_flag)
        g.addSpacing(10)
        self.lbl_mc = QLabel("")
        g.addWidget(self.lbl_mc)
        g.addStretch(1)
        self.progress_bar = QProgressBar()
        self.progress_bar.setRange(0, 1000)
        self.progress_bar.setFixedWidth(220)
        self.progress_bar.setMaximumHeight(16)
        self.progress_bar.setTextVisible(True)
        self.progress_bar.setVisible(False)
        g.addWidget(self.progress_bar)
        self.interrupt_btn = QPushButton("Interrupt")
        self.interrupt_btn.setToolTip("Stop the running loop after the current spaxel")
        self.interrupt_btn.setStyleSheet("background-color: #e06060; color: white; font-weight: bold")
        self.interrupt_btn.setVisible(False)
        g.addWidget(self.interrupt_btn)
        sep = QFrame()
        sep.setFrameShape(QFrame.Shape.HLine)
        lay.addWidget(sep)
        lay.addWidget(status)          # always visible, just above the action buttons

        # ---- actions (always visible): one column per topic, current selection on top, all below
        grid = QGridLayout()
        grid.setContentsMargins(0, 4, 0, 2)
        grid.setHorizontalSpacing(4)
        grid.setVerticalSpacing(3)
        columns = [
            (("FIT", "FIT", "Fit the current spaxel / mask"),
             ("FIT ALL", "FITALL", "Fit all spaxels in the FITALL range (or all masks)")),
            (("FIT ADJ", "FITADJ", "Fit using the neighbours as initial guess"),
             ("FIT ADJ ALL", "FITADJALL", "Refit bad spaxels next to good ones, iteratively")),
            (("RESET FIT", "RESETFIT", "Reset the fit of this spaxel / mask"),
             ("RESET FIT ALL", "RESETFITALL", "Reset all fits")),
            (("AUTO FLAG", "FLAGALL", "Flag all spaxels / masks with the S/N and error thresholds"),
             ("FLAG ON/OFF", "RESIMAMASK", "Show flagged spaxels in the result maps or hide them")),
            (("RESET USER PARS", "RESETUSER", "Reset user start values and limits"),
             ("RESET ALL PARS", "RESETALL", "Reset fit, user values and settings")),
            (("MASK / SPAXEL", "MODE", "Switch mask / spaxel fitting"),
             ("GAUSS / MOMENTS", "TYPE", "Switch Gaussian fits / moments")),
        ]
        for col, items in enumerate(columns):
            for row, (text, code, tip) in enumerate(items):
                btn = QPushButton(text)
                btn.setObjectName("action")
                btn.setToolTip(tip)
                btn.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
                btn.clicked.connect(lambda _=False, c=code: self.ctl.linefit_action(c))
                grid.addWidget(btn, row, col)
            grid.setColumnStretch(col, 1)
        save = QPushButton("SAVE")
        save.setToolTip("Save the results as FITS")
        save.setObjectName("save")
        save.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        save.clicked.connect(lambda: self.ctl.linefit_action("SAVE"))
        grid.addWidget(save, 0, len(columns), 2, 1)
        grid.setColumnStretch(len(columns), 1)
        lay.addLayout(grid)
        self._building = False
        return panel

    def show_results(self):
        """Expand the Fit results group (Options -> Show linefit window)."""
        self.results_group.btn.setChecked(True)

    # ================================================================== events
    def _text(self, code, edit):
        if self._building:
            return
        txt = edit.text().strip()
        if txt.lower() in ("not set", "", "auto"):
            val = NOT_FIT
        else:
            try:
                val = float(txt)
            except ValueError:
                utils.warn(f"Invalid number: {txt!r}")
                self.update_all(update_userpars=True)
                return
        self.ctl.linefit_action(code, val)

    # ================================================================== refresh
    def typeswitch(self):
        st = self.state
        gauss = st.linefit_type == FIT_GAUSS
        self.lbl_type.setText("GAUSS" if gauss else "MOMENTS")
        first, second = {0: ("1st", "2nd"), 1: ("blue", "red"), 2: ("narrow", "broad")}[int(st.secondcomp_mode)]
        for (lt, par), d in self.w.items():
            if lt == "N":
                d["comp"].setText(first if gauss else ("main" if par <= 2 else "flux"))
            elif lt == "B":
                d["comp"].setText(second if gauss else ("" if par <= 2 else "1st/2nd"))
            for key in ("start", "min", "max", "fix"):
                if key in d:
                    d[key].setEnabled(gauss)
        for j in (3, 5, 6, 7, 8):
            self.col_labels[j].setText(self.COLS[j] if gauss else "")
        self.cb_constr.setVisible(gauss)
        self.lbl_momthresh.setVisible(not gauss)
        self.edit_momthresh.setVisible(not gauss)
        self.cb_fixratios.setVisible(gauss)
        self.lbl_second.setVisible(gauss)
        for rb in self.rb_second:
            rb.setVisible(gauss)
        self.cb_smart.setVisible(gauss)

    def update_all(self, update_userpars: bool = False):
        st = self.state
        self._building = True
        try:
            rs = st.get_results()
            if st.linefit_mode == MODE_SPAXEL:
                idx = rs.index(col=st.col, row=st.row)
                self.lbl_mode.setText("SPAXEL")
                self.lbl_sel.setText(f"[{st.col},{st.row}]")
                self.btn_prevmask.setVisible(False)
                self.btn_nextmask.setVisible(False)
            else:
                idx = rs.index(imask=st.imask)
                self.lbl_mode.setText("MASK")
                self.lbl_sel.setText(f"{st.imask}/{st.Nmask}")
                self.btn_prevmask.setVisible(True)
                self.btn_nextmask.setVisible(True)
            n, b, c, m = rs.n[idx], rs.b[idx], rs.c[idx], rs.m[idx]
            nerr, berr, cerr, merr = rs.nerr[idx], rs.berr[idx], rs.cerr[idx], rs.merr[idx]
            gauss = st.linefit_type == FIT_GAUSS
            if gauss:
                lt = st.par_imagebutton[:1]
                flag = {"B": b[0], "C": c[0]}.get(lt, n[0])
                for par in range(1, 3 + st.Nlines):
                    self.w[("N", par)]["best"].setText(utils.resultsstring(n[par], nerr[par, :2]))
                    self.w[("B", par)]["best"].setText(utils.resultsstring(b[par], berr[par, :2]))
                    if par >= 3:
                        self.w[("C", par)]["best"].setText(_cont_string(c[par - 2], cerr[par - 2, :2]))
                self.lbl_chisq.setText(f"{rs.chisq[idx]:.4f}")
            else:
                flag = 0
                if st.Nlines > 0:
                    mi = st.mainline_index()
                    flag = m[6 * mi + 5]
                    for par in (1, 2):
                        self.w[("N", par)]["best"].setText(utils.resultsstring(m[6 * mi + par], merr[6 * mi + par, :2]))
                        self.w[("B", par)]["best"].setText("")
                    for iline in range(st.Nlines):
                        rp, par = 6 * iline, iline + 3
                        self.w[("N", par)]["best"].setText(utils.resultsstring(m[rp], merr[rp, :2]))
                        self.w[("B", par)]["lamb"].setText(utils.resultsstring(m[rp + 1], merr[rp + 1, :2]))
                        self.w[("B", par)]["best"].setText(utils.resultsstring(m[rp + 2], merr[rp + 2, :2]))
                        self.w[("C", par)]["best"].setText(_cont_string(c[iline + 1], cerr[iline + 1, :2]))
                self.lbl_chisq.setText("--")
            for iline in range(st.Nlines):
                par = iline + 3
                if gauss:
                    ndv = 0.0 if nerr[1, 0] == NOT_FIT else n[1]
                    bdv = 0.0 if berr[1, 0] == NOT_FIT else b[1]
                    self.w[("N", par)]["lamb"].setText(f"{st.lines[iline] * (1 + ndv / CKMS):.2f}")
                    self.w[("B", par)]["lamb"].setText(f"{st.lines[iline] * (1 + bdv / CKMS):.2f}")
                else:
                    self.w[("N", par)]["lamb"].setText(f"{st.lines[iline] * (1 + m[6 * iline + 1] / CKMS):.2f}")
            self.set_flag_button(flag)

            self.edit_z.setText(f"{st.redshift:.6f}")
            R = float(st.getinstrres()) if st.Nlines > 0 else 0.0
            self.lbl_instrres.setText(f"{R:.1f}")
            self.lbl_instrres_mode.setText(INSTRRES_MODE_TEXT.get(st.instrres_mode, ""))
            self.lbl_errmethod.setText(ERR_METHOD_NAMES.get(st.domontecarlo, ""))
            onoff = lambda v: "on" if v else "off"  # noqa: E731
            inmap = st.gauss_initmap is not None if gauss else st.mom_windowmap is not None
            self.lbl_mc.setText(f"MC plots {onoff(st.plotMonteCarlodistrib)} · MC PDFs {onoff(st.saveMonteCarlodistrib)} · "
                                f"MC noise {onoff(st.useMonteCarlonoise)} · scale errors {onoff(st.scaleNoiseerrors)} · "
                                f"start-value maps {onoff(inmap)}")
            self.edit_maxwoffb.setText(f"{st.maxwoffb:.2f}")
            self.edit_maxwoffr.setText(f"{st.maxwoffr:.2f}")
            self.edit_momthresh.setText(f"{st.mom_thresh:.2f}")
            self.edit_snthresh.setText(f"{st.mask_sn_thresh:.2f}")
            self.edit_maxvelerr.setText(f"{st.mask_maxvelerr:.2f}")
            self.cb_contmode.setChecked(bool(st.continuumfit_mode))
            self.edit_cminoff.setText(f"{st.continuumfit_minoff:.2f}")
            self.edit_cmaxoff.setText(f"{st.continuumfit_maxoff:.2f}")
            self.edit_cminperc.setText(f"{st.continuumfit_minperc:.2f}")
            self.edit_cmaxperc.setText(f"{st.continuumfit_maxperc:.2f}")
            self.edit_corder.setText(f"{int(st.continuumfit_order)}")
            self.btn_tpl.setEnabled(st.instrres_tplsig > 0)

            self.image_group.setExclusive(False)
            for btn in self.image_group.buttons():
                btn.setChecked(False)
            self.image_group.setExclusive(True)
            pb = st.par_imagebutton
            if pb == "FLAG":
                self.rb_flag.setChecked(True)
            elif pb == "CHISQ":
                self.rb_chisq.setChecked(True)
            elif pb:
                key = (pb[0], int(pb[1:]))
                if key in self.w:
                    self.w[key]["image"].setChecked(True)

            if update_userpars:
                for par in range(1, 3 + st.Nlines):
                    for lt, gfit, lims, fix in (("N", st.gnfit, st.gnlims, st.pnfix), ("B", st.gbfit, st.gblims, st.pbfix)):
                        d = self.w[(lt, par)]
                        d["start"].setText(_userparstring(gfit[par]))
                        d["min"].setText(_userparstring(lims[par, 0]))
                        d["max"].setText(_userparstring(lims[par, 1]))
                        d["fix"].setChecked(bool(fix[par]))
                    if par > 2:
                        line = par - 3
                        self.w[("N", par)]["fit"].setChecked(bool(st.pndofit[line]))
                        self.w[("B", par)]["fit"].setChecked(bool(st.pbdofit[line]))
                        self.w[("C", par)]["fit"].setChecked(bool(st.pcdofit[line]))
                        self.w[("N", par)]["show"].setChecked(bool(st.nshow[line]))
                        self.w[("B", par)]["show"].setChecked(bool(st.bshow[line]))
                        self.w[("C", par)]["show"].setChecked(bool(st.cshow[line]))
                self.cb_constr.setChecked(bool(st.fitconstr))
                self.cb_fixratios.setChecked(bool(st.fitfixratios))
                self.rb_second[int(st.secondcomp_mode)].setChecked(True)
                self.cb_smart.setChecked(bool(st.secondcomp_smart))
                coeffs = st.instrres_varpoly if st.instrres_mode == INSTRRES_VARPOLY else st.instrres_extpoly
                for i, e in enumerate(self.edit_coeff):
                    e.setText(f"{coeffs[i]:g}")
                if np.any(st.pbdofit[:st.Nlines]) and not self.cb_second.isChecked():
                    self.cb_second.setChecked(True)
        finally:
            self._building = False

    def set_flag_button(self, flag):
        bad = flag > 0
        self.btn_flag.setText("BAD" if bad else "OK")
        self.btn_flag.setStyleSheet("background-color: #e06060; color: white; font-weight: bold" if bad
                                    else "background-color: #60c060; color: white; font-weight: bold")
