import sys
from PyQt6.QtWidgets import (QDialog, QVBoxLayout, QHBoxLayout, QGridLayout, 
                             QLabel, QCheckBox, QLineEdit, QPushButton, QGroupBox,
                             QScrollArea, QWidget)
from PyQt6.QtCore import Qt
from pykubeviz.core.state import KubevizState
from pykubeviz.core.lines import LINES_DB

class FitLinesDialog(QDialog):
    def __init__(self, state: KubevizState, parent=None):
        super().__init__(parent)
        self.state = state
        self.setWindowTitle("Fit Lines & Limits")
        self.setMinimumWidth(500)
        self.setMinimumHeight(600)
        self.init_ui()
        self.load_state()

    def init_ui(self):
        main_layout = QVBoxLayout(self)

        # Kinematics Limits
        kin_group = QGroupBox("Kinematics Limits")
        kin_layout = QGridLayout(kin_group)
        
        kin_layout.addWidget(QLabel("Velocity (km/s)"), 0, 0)
        self.vel_min = QLineEdit()
        self.vel_max = QLineEdit()
        self.vel_min.setPlaceholderText("Min (-999)")
        self.vel_max.setPlaceholderText("Max (-999)")
        kin_layout.addWidget(self.vel_min, 0, 1)
        kin_layout.addWidget(self.vel_max, 0, 2)
        
        kin_layout.addWidget(QLabel("Sigma (km/s)"), 1, 0)
        self.sig_min = QLineEdit()
        self.sig_max = QLineEdit()
        self.sig_min.setPlaceholderText("Min (-999)")
        self.sig_max.setPlaceholderText("Max (-999)")
        kin_layout.addWidget(self.sig_min, 1, 1)
        kin_layout.addWidget(self.sig_max, 1, 2)
        
        main_layout.addWidget(kin_group)

        # Emission Lines
        lines_group = QGroupBox("Emission Lines (Flux Limits)")
        lines_layout = QVBoxLayout(lines_group)
        
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll_content = QWidget()
        scroll_layout = QGridLayout(scroll_content)
        
        headers = ["Fit?", "Line", "Flux Min", "Flux Max"]
        for col, h in enumerate(headers):
            scroll_layout.addWidget(QLabel(h), 0, col)
            
        self.line_widgets = {}
        
        for i, line in enumerate(LINES_DB, start=1):
            key = line['name']
            label = f"{line['fancy_name']} ({line['rest_wave']})"
            
            cb = QCheckBox()
            cb.setChecked(True)
            scroll_layout.addWidget(cb, i, 0)
            scroll_layout.addWidget(QLabel(label), i, 1)
            
            fmin = QLineEdit()
            fmax = QLineEdit()
            fmin.setPlaceholderText("Min (-999)")
            fmax.setPlaceholderText("Max (-999)")
            scroll_layout.addWidget(fmin, i, 2)
            scroll_layout.addWidget(fmax, i, 3)
            
            self.line_widgets[key] = {'cb': cb, 'min': fmin, 'max': fmax}

        scroll_layout.setRowStretch(len(LINES_DB) + 1, 1)
        scroll.setWidget(scroll_content)
        lines_layout.addWidget(scroll)
        
        # Selection buttons
        sel_btn_layout = QHBoxLayout()
        btn_sel_ha_nii = QPushButton("Select only Ha+[NII]")
        btn_sel_ha_nii.clicked.connect(self.select_ha_nii)
        btn_sel_all = QPushButton("Select All")
        btn_sel_all.clicked.connect(self.select_all)
        btn_desel_all = QPushButton("De-Select All")
        btn_desel_all.clicked.connect(self.deselect_all)
        
        sel_btn_layout.addWidget(btn_sel_ha_nii)
        sel_btn_layout.addWidget(btn_sel_all)
        sel_btn_layout.addWidget(btn_desel_all)
        lines_layout.addLayout(sel_btn_layout)
        
        main_layout.addWidget(lines_group)

        # Buttons
        btn_layout = QHBoxLayout()
        save_btn = QPushButton("Save & Apply")
        save_btn.clicked.connect(self.save_state)
        cancel_btn = QPushButton("Cancel")
        cancel_btn.clicked.connect(self.reject)
        
        btn_layout.addStretch()
        btn_layout.addWidget(save_btn)
        btn_layout.addWidget(cancel_btn)
        
        main_layout.addLayout(btn_layout)

    def load_state(self):
        conf = self.state.fit_lines_config
        
        def set_val(le, val):
            if val is not None:
                le.setText(str(val))
            else:
                le.clear()
                
        set_val(self.vel_min, conf.get('vel', {}).get('min'))
        set_val(self.vel_max, conf.get('vel', {}).get('max'))
        set_val(self.sig_min, conf.get('sigma', {}).get('min'))
        set_val(self.sig_max, conf.get('sigma', {}).get('max'))
        
        for line in LINES_DB:
            key = line['name']
            w = self.line_widgets[key]
            c = conf.get(key, {})
            w['cb'].setChecked(c.get('enabled', True))
            set_val(w['min'], c.get('min'))
            set_val(w['max'], c.get('max'))

    def save_state(self):
        def parse_val(text):
            text = text.strip()
            if not text or text == '-999':
                return None
            try:
                return float(text)
            except ValueError:
                return None
                
        conf = self.state.fit_lines_config
        if 'vel' not in conf: conf['vel'] = {}
        if 'sigma' not in conf: conf['sigma'] = {}
        
        conf['vel']['min'] = parse_val(self.vel_min.text())
        conf['vel']['max'] = parse_val(self.vel_max.text())
        conf['sigma']['min'] = parse_val(self.sig_min.text())
        conf['sigma']['max'] = parse_val(self.sig_max.text())
        
        for line in LINES_DB:
            key = line['name']
            if key not in conf: conf[key] = {}
            w = self.line_widgets[key]
            conf[key]['enabled'] = w['cb'].isChecked()
            conf[key]['min'] = parse_val(w['min'].text())
            conf[key]['max'] = parse_val(w['max'].text())
            
        self.accept()

    def select_ha_nii(self):
        ha_nii_keys = ['Ha', 'n2_b', 'n2_r']
        for key, w in self.line_widgets.items():
            w['cb'].setChecked(key in ha_nii_keys)

    def select_all(self):
        for w in self.line_widgets.values():
            w['cb'].setChecked(True)

    def deselect_all(self):
        for w in self.line_widgets.values():
            w['cb'].setChecked(False)
