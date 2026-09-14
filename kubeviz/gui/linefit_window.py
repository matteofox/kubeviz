"""The linefit window (``kubeviz_setup_linefit`` / ``kubeviz_linefit_update`` /
``kubeviz_linefit_typeswitch``).

All user actions are forwarded to ``controller.linefit_action(code, value)`` with the
same codes the IDL event handler used (``'FIT'``, ``'FITALL'``, ``'SPN3'``, ``'MINPB1'``,
``'FIXN2'``, ``'IMAGEC4'``...), so the two implementations stay easy to compare.
"""
from __future__ import annotations

import numpy as np
from PyQt6.QtCore import Qt
from PyQt6.QtGui import QFont
from PyQt6.QtWidgets import (QButtonGroup, QCheckBox, QFrame, QGridLayout, QHBoxLayout, QLabel,
                             QLineEdit, QMainWindow, QPushButton, QRadioButton, QScrollArea,
                             QVBoxLayout, QWidget)

from .. import utils
from ..constants import (CKMS, ERR_METHOD_NAMES, FIT_GAUSS, INSTRRES_EXTPOLY, INSTRRES_TEMPLATE,
                         INSTRRES_VARPOLY, MODE_SPAXEL, NOT_FIT)

INSTRRES_MODE_TEXT = {INSTRRES_VARPOLY: "(Polynomial fit to cube variance)",
                      INSTRRES_EXTPOLY: "(Polynomial fit to external arcs)",
                      INSTRRES_TEMPLATE: "(Spectral Templates)"}


def _userparstring(v) -> str:
    return "Not Set" if v == NOT_FIT else f"{v:.2f}"


def _edit(width=64, text=""):
    e = QLineEdit(text)
    e.setFixedWidth(width)
    return e


def _hline():
    f = QFrame()
    f.setFrameShape(QFrame.Shape.HLine)
    return f


class LinefitWindow(QMainWindow):
    COLS = ["", "component", "lambda_cen", "start", "best", "min", "max", "fix", "reset", "image", "fit", "show"]

    def __init__(self, controller):
        super().__init__()
        self.ctl = controller
        self.state = controller.state
        self.setWindowTitle(f"linefit: {self.state.filename}")
        self.w = {}                  # (linetype, par) -> dict of widgets
        self.image_group = QButtonGroup(self)
        self.image_group.setExclusive(True)
        self._building = False
        self._build()

    # ================================================================== construction
    def _build(self):
        self._building = True
        st = self.state
        central = QWidget()
        outer = QVBoxLayout(central)
        outer.setContentsMargins(4, 4, 4, 4)

        # ---------------- header
        hdr = QHBoxLayout()
        hdr.addWidget(QLabel("Fit Type:"))
        self.lbl_type = QLabel("GAUSS")
        self.lbl_type.setStyleSheet("font-weight: bold")
        hdr.addWidget(self.lbl_type)
        hdr.addSpacing(12)
        hdr.addWidget(QLabel("Error method:"))
        self.lbl_errmethod = QLabel("")
        hdr.addWidget(self.lbl_errmethod)
        hdr.addSpacing(12)
        hdr.addWidget(QLabel("z ="))
        self.edit_z = _edit(90, f"{st.redshift:.6f}")
        self.edit_z.editingFinished.connect(lambda: self._text("REDSHIFT", self.edit_z))
        hdr.addWidget(self.edit_z)
        hdr.addSpacing(12)
        self.lbl_mode = QLabel("SPAXEL")
        self.lbl_mode.setStyleSheet("font-weight: bold")
        hdr.addWidget(self.lbl_mode)
        self.lbl_sel = QLabel("")
        hdr.addWidget(self.lbl_sel)
        self.btn_prevmask = QPushButton("-")
        self.btn_nextmask = QPushButton("+")
        for b, code in ((self.btn_prevmask, "PREVMASK"), (self.btn_nextmask, "NEXTMASK")):
            b.setFixedWidth(28)
            b.clicked.connect(lambda _=False, c=code: self.ctl.linefit_action(c))
            hdr.addWidget(b)
        hdr.addStretch(1)
        outer.addLayout(hdr)

        # ---------------- parameter table
        self.grid = QGridLayout()
        self.grid.setHorizontalSpacing(6)
        self.grid.setVerticalSpacing(2)
        self.col_labels = []
        for j, name in enumerate(self.COLS):
            lbl = QLabel(name)
            lbl.setAlignment(Qt.AlignmentFlag.AlignCenter)
            f = QFont()
            f.setBold(True)
            lbl.setFont(f)
            self.grid.addWidget(lbl, 0, j)
            self.col_labels.append(lbl)
        row = 1
        row = self._add_kin_rows(row, 1, "Line offset", "(km/s)")
        row = self._add_kin_rows(row, 2, "Line width", "(km/s)")
        for iline in range(st.Nlines):
            par = iline + 3
            row = self._add_line_rows(row, par, iline)
        table = QWidget()
        table.setLayout(self.grid)
        outer.addWidget(table)
        outer.addWidget(_hline())

        # ---------------- chisq / flag / autoflag thresholds
        fr = QHBoxLayout()
        fr.addWidget(QLabel("Red-chi-sq:"))
        self.lbl_chisq = QLabel("0.0000")
        self.lbl_chisq.setMinimumWidth(70)
        fr.addWidget(self.lbl_chisq)
        fr.addWidget(QLabel("image:"))
        self.rb_chisq = QRadioButton()
        self.image_group.addButton(self.rb_chisq)
        self.rb_chisq.clicked.connect(lambda: self.ctl.linefit_action("IMAGECHISQ"))
        fr.addWidget(self.rb_chisq)
        fr.addSpacing(16)
        fr.addWidget(QLabel("Flag:"))
        self.btn_flag = QPushButton("OK")
        self.btn_flag.setFixedWidth(60)
        self.btn_flag.setToolTip("Toggle fit status good or bad")
        self.btn_flag.clicked.connect(lambda: self.ctl.linefit_action("FLAG"))
        fr.addWidget(self.btn_flag)
        fr.addWidget(QLabel("image:"))
        self.rb_flag = QRadioButton()
        self.image_group.addButton(self.rb_flag)
        self.rb_flag.clicked.connect(lambda: self.ctl.linefit_action("IMAGEFLAG"))
        fr.addWidget(self.rb_flag)
        fr.addSpacing(16)
        fr.addWidget(QLabel("Autoflag S/N thresh:"))
        self.edit_snthresh = _edit(64, f"{st.mask_sn_thresh:.2f}")
        self.edit_snthresh.editingFinished.connect(lambda: self._text("MASKSNTHRESH", self.edit_snthresh))
        fr.addWidget(self.edit_snthresh)
        fr.addWidget(QLabel("Max Vel/Sig error (km/s):"))
        self.edit_maxvelerr = _edit(64, f"{st.mask_maxvelerr:.2f}")
        self.edit_maxvelerr.editingFinished.connect(lambda: self._text("MASKMAXVELERR", self.edit_maxvelerr))
        fr.addWidget(self.edit_maxvelerr)
        fr.addStretch(1)
        outer.addLayout(fr)
        outer.addWidget(_hline())

        # ---------------- fitting range / continuum
        cr = QGridLayout()
        cr.addWidget(QLabel("Fitting Range (l/r side, Å)"), 0, 0, 1, 2)
        cr.addWidget(QLabel("Continuum Fit:"), 0, 2)
        cr.addWidget(QLabel("Range of Offset (min/max, Å)"), 0, 3, 1, 2)
        cr.addWidget(QLabel("Range of Percentiles included"), 0, 5, 1, 2)
        cr.addWidget(QLabel("Poly. Order"), 0, 7)
        self.edit_maxwoffb = _edit(64)
        self.edit_maxwoffr = _edit(64)
        self.cb_contmode = QCheckBox("MPFIT CONT")
        self.cb_contmode.clicked.connect(lambda: self.ctl.linefit_action("CONTMODE"))
        self.edit_cminoff, self.edit_cmaxoff = _edit(64), _edit(64)
        self.edit_cminperc, self.edit_cmaxperc = _edit(64), _edit(64)
        self.edit_corder = _edit(40)
        for e, code in ((self.edit_maxwoffb, "LINEFITMAXOFFB"), (self.edit_maxwoffr, "LINEFITMAXOFFR"),
                        (self.edit_cminoff, "CONTMINOFF"), (self.edit_cmaxoff, "CONTMAXOFF"),
                        (self.edit_cminperc, "CONTMINPERC"), (self.edit_cmaxperc, "CONTMAXPERC"),
                        (self.edit_corder, "CONTORDER")):
            e.editingFinished.connect(lambda c=code, ee=e: self._text(c, ee))
        cr.addWidget(self.edit_maxwoffb, 1, 0)
        cr.addWidget(self.edit_maxwoffr, 1, 1)
        cr.addWidget(self.cb_contmode, 1, 2)
        cr.addWidget(self.edit_cminoff, 1, 3)
        cr.addWidget(self.edit_cmaxoff, 1, 4)
        cr.addWidget(self.edit_cminperc, 1, 5)
        cr.addWidget(self.edit_cmaxperc, 1, 6)
        cr.addWidget(self.edit_corder, 1, 7)
        cr.setColumnStretch(8, 1)
        outer.addLayout(cr)
        outer.addWidget(_hline())

        # ---------------- fit options
        op = QHBoxLayout()
        self.lbl_constr = QLabel("Fit with constraints?")
        op.addWidget(self.lbl_constr)
        self.cb_constr = QCheckBox()
        self.cb_constr.clicked.connect(lambda: self.ctl.linefit_action("FITCONSTRAINTS"))
        op.addWidget(self.cb_constr)
        self.edit_momthresh = _edit(64, f"{st.mom_thresh:.2f}")
        self.edit_momthresh.editingFinished.connect(lambda: self._text("MOMTHRESH", self.edit_momthresh))
        op.addWidget(self.edit_momthresh)
        op.addSpacing(12)
        self.lbl_fixratios = QLabel("Fix line ratios?")
        op.addWidget(self.lbl_fixratios)
        self.cb_fixratios = QCheckBox()
        self.cb_fixratios.clicked.connect(lambda: self.ctl.linefit_action("FIXRATIOS"))
        op.addWidget(self.cb_fixratios)
        op.addSpacing(12)
        self.lbl_second = QLabel("2nd component:")
        op.addWidget(self.lbl_second)
        self.rb_second = [QRadioButton("FAINTER"), QRadioButton("LARGER OFFSET"), QRadioButton("LARGER WIDTH")]
        self.second_group = QButtonGroup(self)
        for i, (rb, code) in enumerate(zip(self.rb_second, ("FREE_2ND", "HIVEL_2ND", "BROAD_2ND"))):
            self.second_group.addButton(rb, i)
            rb.clicked.connect(lambda _=False, c=code: self.ctl.linefit_action(c))
            op.addWidget(rb)
        self.lbl_smart = QLabel("Smart?")
        op.addWidget(self.lbl_smart)
        self.cb_smart = QCheckBox()
        self.cb_smart.clicked.connect(lambda: self.ctl.linefit_action("SMART_2ND"))
        op.addWidget(self.cb_smart)
        op.addStretch(1)
        outer.addLayout(op)

        st2 = QHBoxLayout()
        st2.addWidget(QLabel("Montecarlo PDFs plot/save:"))
        self.lbl_mcplot = QLabel("OFF")
        st2.addWidget(self.lbl_mcplot)
        st2.addWidget(QLabel("/"))
        self.lbl_mcsave = QLabel("OFF")
        st2.addWidget(self.lbl_mcsave)
        st2.addSpacing(12)
        st2.addWidget(QLabel("Montecarlo noise:"))
        self.lbl_mcnoise = QLabel("OFF")
        st2.addWidget(self.lbl_mcnoise)
        st2.addSpacing(12)
        st2.addWidget(QLabel("Scale Noisecube error:"))
        self.lbl_scnoise = QLabel("OFF")
        st2.addWidget(self.lbl_scnoise)
        st2.addSpacing(12)
        st2.addWidget(QLabel("Input Startvals File:"))
        self.lbl_inputmap = QLabel("OFF")
        st2.addWidget(self.lbl_inputmap)
        st2.addStretch(1)
        outer.addLayout(st2)
        outer.addWidget(_hline())

        # ---------------- instrumental resolution
        ir = QHBoxLayout()
        ir.addWidget(QLabel("Compute Instrumental Resolution:"))
        b = QPushButton("FIT TO VARIANCE")
        b.setToolTip("Fit to lines in variance spectrum")
        b.clicked.connect(lambda: self.ctl.linefit_action("FITSKY"))
        ir.addWidget(b)
        self.btn_poly = QPushButton("USE POLYNOMIAL")
        self.btn_poly.setToolTip("Use the polynomial from the header / archive / instrument manual")
        self.btn_poly.clicked.connect(lambda: self.ctl.linefit_action("POLYSKY"))
        ir.addWidget(self.btn_poly)
        self.btn_tpl = QPushButton("USE TEMPLATES")
        self.btn_tpl.setEnabled(st.instrres_tplsig > 0)
        self.btn_tpl.clicked.connect(lambda: self.ctl.linefit_action("TPLSKY"))
        ir.addWidget(self.btn_tpl)
        ir.addSpacing(12)
        ir.addWidget(QLabel("R (at main line):"))
        self.lbl_instrres = QLabel("")
        self.lbl_instrres.setMinimumWidth(70)
        ir.addWidget(self.lbl_instrres)
        self.lbl_instrres_mode = QLabel("")
        ir.addWidget(self.lbl_instrres_mode)
        ir.addStretch(1)
        outer.addLayout(ir)
        ir2 = QHBoxLayout()
        ir2.addWidget(QLabel("Polynomial coefficients:"))
        self.edit_coeff = []
        for i in range(st.max_polycoeff_instrres + 1):
            e = _edit(110)
            e.editingFinished.connect(lambda i=i, ee=e: self._text(f"POLYCOEFFPAR{i}", ee))
            self.edit_coeff.append(e)
            ir2.addWidget(e)
        ir2.addStretch(1)
        outer.addLayout(ir2)
        outer.addWidget(_hline())

        # ---------------- action buttons
        b1 = QHBoxLayout()
        for text, code, tip in (("FIT", "FIT", "Fit lines"), ("FIT ADJ", "FITADJ", "Fit using adjacent spaxels as initial guess"),
                                ("RESET FIT", "RESETFIT", "Reset fit results"), (None, None, None),
                                ("AUTO FLAG", "FLAGALL", "Auto-Flag all masks / spaxels"), (None, None, None),
                                ("RESET ALL PARS", "RESETALL", "Reset all parameters (fit, user, settings)"),
                                ("RESET USER PARS", "RESETUSER", "Reset user defined constraints")):
            if text is None:
                b1.addSpacing(30)
                continue
            btn = QPushButton(text)
            btn.setToolTip(tip)
            btn.clicked.connect(lambda _=False, c=code: self.ctl.linefit_action(c))
            b1.addWidget(btn)
        b1.addStretch(1)
        outer.addLayout(b1)
        b2 = QHBoxLayout()
        for text, code, tip in (("FIT ALL", "FITALL", "Fit all masks / spaxels"),
                                ("FIT ADJ ALL", "FITADJALL", "Fit bad spaxels using initial guess from adjacent spaxels"),
                                ("RESET FIT ALL", "RESETFITALL", "Reset fit results"), (None, None, None),
                                ("FLAG ON/OFF", "RESIMAMASK", "Switch masking of results on/off"), (None, None, None),
                                ("SAVE", "SAVE", "Save fit results in FITS format"),
                                ("MASK/SPAXEL", "MODE", "Switch MASK/SPAXEL fitting modes"),
                                ("GAUSS/MOMENTS", "TYPE", "Switch GAUSSIAN/MOMENTS fitting types")):
            if text is None:
                b2.addSpacing(30)
                continue
            btn = QPushButton(text)
            btn.setToolTip(tip)
            btn.clicked.connect(lambda _=False, c=code: self.ctl.linefit_action(c))
            b2.addWidget(btn)
        b2.addStretch(1)
        outer.addLayout(b2)
        outer.addStretch(1)

        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setWidget(central)
        self.setCentralWidget(scroll)
        self._building = False
        nrows = 4 + 3 * st.Nlines
        self.resize(1120, min(240 + 34 * nrows + 260, 900))
        self.typeswitch()
        self.update_all(update_userpars=True)

    def _add_param_row(self, row, label, lt, par, comp_text, with_fit_show=True, is_cont=False):
        d = {}
        self.grid.addWidget(QLabel(label), row, 0)
        d["comp"] = QLabel(comp_text)
        self.grid.addWidget(d["comp"], row, 1)
        d["lamb"] = QLabel("")
        self.grid.addWidget(d["lamb"], row, 2)
        code = f"{lt}{par}"
        if not is_cont:
            d["start"] = _edit(70)
            d["start"].editingFinished.connect(lambda c="SP" + code, e=d["start"]: self._text(c, e))
            self.grid.addWidget(d["start"], row, 3)
        d["best"] = QLabel("Not Fit")
        d["best"].setMinimumWidth(150)
        self.grid.addWidget(d["best"], row, 4)
        if not is_cont:
            d["min"] = _edit(64)
            d["max"] = _edit(64)
            d["min"].editingFinished.connect(lambda c="MINP" + code, e=d["min"]: self._text(c, e))
            d["max"].editingFinished.connect(lambda c="MAXP" + code, e=d["max"]: self._text(c, e))
            self.grid.addWidget(d["min"], row, 5)
            self.grid.addWidget(d["max"], row, 6)
            d["fix"] = QCheckBox()
            d["fix"].clicked.connect(lambda _=False, c="FIX" + code: self.ctl.linefit_action(c))
            self.grid.addWidget(d["fix"], row, 7, alignment=Qt.AlignmentFlag.AlignCenter)
            d["reset"] = QPushButton("x")
            d["reset"].setFixedWidth(26)
            d["reset"].setToolTip("Reset this parameter")
            d["reset"].clicked.connect(lambda _=False, c="RESET" + code: self.ctl.linefit_action(c))
            self.grid.addWidget(d["reset"], row, 8, alignment=Qt.AlignmentFlag.AlignCenter)
        d["image"] = QRadioButton()
        d["image"].setToolTip("Show this parameter as a map in the spaxel viewer")
        self.image_group.addButton(d["image"])
        d["image"].clicked.connect(lambda _=False, c="IMAGE" + code: self.ctl.linefit_action(c))
        self.grid.addWidget(d["image"], row, 9, alignment=Qt.AlignmentFlag.AlignCenter)
        if with_fit_show:
            d["fit"] = QCheckBox()
            d["fit"].clicked.connect(lambda _=False, c="FIT" + code: self.ctl.linefit_action(c))
            self.grid.addWidget(d["fit"], row, 10, alignment=Qt.AlignmentFlag.AlignCenter)
            d["show"] = QCheckBox()
            d["show"].clicked.connect(lambda _=False, c="SHOW" + code: self.ctl.linefit_action(c))
            self.grid.addWidget(d["show"], row, 11, alignment=Qt.AlignmentFlag.AlignCenter)
        self.w[(lt, par)] = d
        return row + 1

    def _add_kin_rows(self, row, par, label, unit):
        row = self._add_param_row(row, label, "N", par, "Narrow", with_fit_show=False)
        row = self._add_param_row(row, unit, "B", par, "Broad", with_fit_show=False)
        return row

    def _add_line_rows(self, row, par, iline):
        name = self.state.linefancynames[iline]
        row = self._add_param_row(row, name, "N", par, "Narrow")
        row = self._add_param_row(row, "", "B", par, "Broad")
        row = self._add_param_row(row, "", "C", par, "Continuum", is_cont=True)
        return row

    # ================================================================== events
    def _text(self, code, edit):
        if self._building:
            return
        txt = edit.text().strip()
        if txt.lower() in ("not set", ""):
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
        """Relabel the window for Gaussian or moments fitting (``kubeviz_linefit_typeswitch``)."""
        st = self.state
        gauss = st.linefit_type == FIT_GAUSS
        self.lbl_type.setText("GAUSS" if gauss else "MOMENTS")
        if gauss:
            first, second = {0: ("1st Comp.", "2nd Comp."), 1: ("Blue", "Red"), 2: ("Narrow", "Broad")}[st.secondcomp_mode]
        for (lt, par), d in self.w.items():
            if lt == "N":
                d["comp"].setText((first if gauss else ("Main Line" if par <= 2 else "Flux")))
            elif lt == "B":
                d["comp"].setText((second if gauss else ("" if par <= 2 else "1st/2nd")))
            for key in ("start", "min", "max", "fix"):
                if key in d:
                    d[key].setVisible(gauss)
            if lt == "B" and "image" in d:
                d["image"].setVisible(gauss)
        for j, name in enumerate(self.COLS):
            if j in (3, 5, 6, 7, 8):
                self.col_labels[j].setText(name if gauss else "")
        self.lbl_constr.setText("Fit with constraints?" if gauss else "Moments thresh:")
        self.cb_constr.setVisible(gauss)
        self.edit_momthresh.setVisible(not gauss)
        self.lbl_fixratios.setVisible(gauss)
        self.cb_fixratios.setVisible(gauss)
        self.lbl_second.setVisible(gauss)
        for rb in self.rb_second:
            rb.setVisible(gauss)
        self.lbl_smart.setVisible(gauss)
        self.cb_smart.setVisible(gauss)

    def update_all(self, update_userpars: bool = False):
        """Port of ``kubeviz_linefit_update``."""
        st = self.state
        self._building = True
        try:
            rs = st.get_results()
            if st.linefit_mode == MODE_SPAXEL:
                idx = rs.index(col=st.col, row=st.row)
                self.lbl_mode.setText("SPAXEL")
                self.lbl_sel.setText(f"[{st.col},{st.row}]")
                chisq = f"{rs.chisq[idx]:.4f}"
                self.btn_prevmask.setVisible(False)
                self.btn_nextmask.setVisible(False)
            else:
                idx = rs.index(imask=st.imask)
                self.lbl_mode.setText("MASK")
                self.lbl_sel.setText(f"{st.imask}/{st.Nmask}")
                chisq = f"{rs.chisq[idx]:.4f}"
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
                        cpar = par - 2
                        self.w[("C", par)]["best"].setText(_cont_string(c[cpar], cerr[cpar, :2]))
                self.lbl_chisq.setText(chisq)
            else:
                flag = 0
                if st.Nlines > 0:
                    mi = st.mainline_index()
                    flag = m[6 * mi + 5]
                    for par in (1, 2):
                        self.w[("N", par)]["best"].setText(utils.resultsstring(m[6 * mi + par], merr[6 * mi + par, :2]))
                        self.w[("B", par)]["best"].setText("")
                    for iline in range(st.Nlines):
                        rp, par, cpar = 6 * iline, iline + 3, iline + 1
                        self.w[("N", par)]["best"].setText(utils.resultsstring(m[rp], merr[rp, :2]))
                        self.w[("B", par)]["lamb"].setText(utils.resultsstring(m[rp + 1], merr[rp + 1, :2]))
                        self.w[("B", par)]["best"].setText(utils.resultsstring(m[rp + 2], merr[rp + 2, :2]))
                        self.w[("C", par)]["best"].setText(_cont_string(c[cpar], cerr[cpar, :2]))
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
            self.lbl_instrres.setText(f"{R:.2f}")
            self.lbl_instrres_mode.setText(INSTRRES_MODE_TEXT.get(st.instrres_mode, ""))
            self.lbl_errmethod.setText(ERR_METHOD_NAMES.get(st.domontecarlo, ""))
            onoff = lambda v: "ON" if v else "OFF"
            self.lbl_mcplot.setText(onoff(st.plotMonteCarlodistrib))
            self.lbl_mcsave.setText(onoff(st.saveMonteCarlodistrib))
            self.lbl_mcnoise.setText(onoff(st.useMonteCarlonoise))
            self.lbl_scnoise.setText(onoff(st.scaleNoiseerrors))
            self.lbl_inputmap.setText(onoff(st.gauss_initmap is not None if gauss else st.mom_windowmap is not None))
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

            # image radio
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
        finally:
            self._building = False

    def set_flag_button(self, flag):
        bad = flag > 0
        self.btn_flag.setText("BAD" if bad else "OK")
        self.btn_flag.setStyleSheet("background-color: #e06060; color: white; font-weight: bold" if bad
                                    else "background-color: #60c060; color: white; font-weight: bold")

    def closeEvent(self, ev):
        self.state.linefitmap = False
        ev.accept()


def _cont_string(value, err):
    err = np.asarray(err, dtype=float)
    if err[0] == NOT_FIT:
        return "Not Fit"
    if err[0] == -998:
        return f"{value:.4f} +/- No Errors"
    return f"{value:.4f}+{err[0]:.4f}/{err[1]:.4f}"
