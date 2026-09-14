import sys
import os
from PyQt6.QtWidgets import (QMainWindow, QApplication, QWidget, QVBoxLayout, 
                             QHBoxLayout, QGridLayout, QSplitter, QMenuBar, QMenu, QFileDialog, QMessageBox,
                             QPushButton, QLabel, QDoubleSpinBox, QGroupBox, QProgressDialog, QComboBox, QInputDialog)
from PyQt6.QtGui import QAction, QActionGroup
from PyQt6.QtCore import Qt, QTimer
import numpy as np
from pykubeviz.ui.viewers import SpaxelViewer, SpectrumViewer
from pykubeviz.core.state import KubevizState
from pykubeviz.io.data_io import load_datacube, save_fit_results
from pykubeviz.analysis.fitting import (fit_emission_line, fit_all_lines, fit_map, 
                                        reconstruct_fit_spectrum, check_autoflag, 
                                        guess_from_adjacent_spaxels, fit_adj_all)
from pykubeviz.ui.fit_lines_dialog import FitLinesDialog

class KubevizMainWindow(QMainWindow):
    def __init__(self, state: KubevizState):
        super().__init__()
        self.state = state
        self.setWindowTitle("Kubeviz - Python Port")
        self.resize(1440, 800)
        self.rescube = None
        self.viewer_mode = "data"
        
        self.init_ui()
        
    def init_ui(self):
        # Setup menus
        menubar = self.menuBar()
        file_menu = menubar.addMenu('&File')
        file_menu.addAction('Open Datacube...', self.open_datacube)
        file_menu.addAction('Load Fit Results...', self.load_fit_results_action)
        file_menu.addAction('Exit', self.close)
        
        spaxel_menu = menubar.addMenu('Spaxel &Viewer')
        select_data_menu = spaxel_menu.addMenu('Select Data')
        select_data_menu.addAction('Data', lambda: self.set_viewer_mode('data'))
        select_data_menu.addAction('Variance', lambda: self.set_viewer_mode('variance'))
        select_data_menu.addAction('Fit Results', lambda: self.set_viewer_mode('fit_results'))
        
        cmap_menu = spaxel_menu.addMenu('Colormaps')
        cmaps = ['viridis', 'plasma', 'inferno', 'magma', 'cividis', 'turbo', 'grey', 'bipolar']
        for cmap in cmaps:
            cmap_menu.addAction(cmap.capitalize(), lambda checked=False, c=cmap: self.spaxel_viewer.set_colormap(c))
        
        fit_controls_menu = menubar.addMenu('Fit &Controls')
        fit_controls_menu.addAction('Fit Lines', self.open_fit_lines_dialog)
        
        self.error_method_menu = fit_controls_menu.addMenu('Error method')
        self.action_err_none = QAction('Noise Cube', self, checkable=True)
        self.action_err_noise = QAction('Noise Scaling', self, checkable=True)
        
        self.err_action_group = QActionGroup(self)
        self.err_action_group.addAction(self.action_err_none)
        self.err_action_group.addAction(self.action_err_noise)
        self.action_err_none.setChecked(True)
        
        self.error_method_menu.addAction(self.action_err_none)
        self.error_method_menu.addAction(self.action_err_noise)
        
        self.action_err_none.triggered.connect(lambda: self.set_error_method('none'))
        self.action_err_noise.triggered.connect(lambda: self.set_error_method('noise_scaling'))
        
        # Main widget and layout
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        main_layout = QVBoxLayout(central_widget)
        
        # Splitter for adjustable panes
        splitter = QSplitter(Qt.Orientation.Horizontal)
        

        # Left pane: Spaxel Viewer
        self.spaxel_viewer = SpaxelViewer()
        splitter.addWidget(self.spaxel_viewer)
        
        # Right pane: Spectrum Viewer and Controls
        right_pane = QWidget()
        right_layout = QVBoxLayout(right_pane)
        right_layout.setContentsMargins(0, 0, 0, 0)
        
        self.spectrum_viewer = SpectrumViewer()
        self.spectrum_viewer.wavelength_clicked.connect(self.on_wavelength_clicked)
        right_layout.addWidget(self.spectrum_viewer, stretch=1)
        
        # Controls panel below spectrum
        controls_group = QGroupBox("Fitting Controls:")
        controls_grid = QGridLayout(controls_group)
        
        # Redshift input
        controls_grid.addWidget(QLabel("Initial Redshift (z):"), 0, 0)
        self.z_input = QDoubleSpinBox()
        self.z_input.setRange(-0.1, 10.0)
        self.z_input.setDecimals(5)
        self.z_input.setSingleStep(0.001)
        self.z_input.setSpecialValueText("Auto-detect")
        self.z_input.setValue(self.z_input.minimum()) # Auto-detect by default
        controls_grid.addWidget(self.z_input, 0, 1)
        
        # Fit Results Display
        self.fit_res_label = QLabel("Fit Results: None")
        font = self.fit_res_label.font()
        font.setFamily("Monospace")
        font.setPointSize(13)
        self.fit_res_label.setFont(font)
        controls_grid.addWidget(self.fit_res_label, 1, 0, 1, 4)
        
        # Continuum Mode
        controls_grid.addWidget(QLabel("Continuum:"), 2, 0)
        self.cont_mode_cb = QComboBox()
        self.cont_mode_cb.addItem("Side Bands (SDSS)", 0)
        self.cont_mode_cb.addItem("Constant Fit", 1)
        self.cont_mode_cb.currentIndexChanged.connect(self.on_cont_mode_changed)
        controls_grid.addWidget(self.cont_mode_cb, 2, 1)

        # Fit Window
        controls_grid.addWidget(QLabel("Fit Window (\u00c5):"), 2, 2)
        self.fit_window_sp = QDoubleSpinBox()
        self.fit_window_sp.setRange(10.0, 5000.0)
        self.fit_window_sp.setValue(self.state.fit_window)
        self.fit_window_sp.setToolTip("+/- Angstroms around the mainline")
        controls_grid.addWidget(self.fit_window_sp, 2, 3)
        
        # Min/Max Off
        controls_grid.addWidget(QLabel("Min/Max Off:"), 3, 0)
        off_layout = QHBoxLayout()
        self.cont_minoff_sp = QDoubleSpinBox()
        self.cont_minoff_sp.setRange(0, 5000)
        self.cont_minoff_sp.setValue(200.0)
        self.cont_minoff_sp.setToolTip("Min Offset (Å)")
        
        self.cont_maxoff_sp = QDoubleSpinBox()
        self.cont_maxoff_sp.setRange(0, 5000)
        self.cont_maxoff_sp.setValue(500.0)
        self.cont_maxoff_sp.setToolTip("Max Offset (Å)")
        off_layout.addWidget(self.cont_minoff_sp)
        off_layout.addWidget(self.cont_maxoff_sp)
        controls_grid.addLayout(off_layout, 3, 1)
        
        # Min/Max P%
        controls_grid.addWidget(QLabel("Min/Max P%:"), 3, 2)
        perc_layout = QHBoxLayout()
        self.cont_minperc_sp = QDoubleSpinBox()
        self.cont_minperc_sp.setRange(0, 100)
        self.cont_minperc_sp.setValue(40.0)
        self.cont_minperc_sp.setToolTip("Min Percentile")
        
        self.cont_maxperc_sp = QDoubleSpinBox()
        self.cont_maxperc_sp.setRange(0, 100)
        self.cont_maxperc_sp.setValue(60.0)
        self.cont_maxperc_sp.setToolTip("Max Percentile")
        perc_layout.addWidget(self.cont_minperc_sp)
        perc_layout.addWidget(self.cont_maxperc_sp)
        controls_grid.addLayout(perc_layout, 3, 3)
        
        # Buttons
        btn_layout1 = QHBoxLayout()
        self.btn_fit_spaxel = QPushButton("Fit Current Spaxel")
        self.btn_fit_spaxel.clicked.connect(self.fit_current_spaxel)
        
        self.btn_fit_adj = QPushButton("Fit Adj")
        self.btn_fit_adj.clicked.connect(self.fitadj_current_spaxel)
        
        self.btn_save_res = QPushButton("Save Fit Results")
        self.btn_save_res.clicked.connect(self.save_results)
        self.btn_save_res.setEnabled(False)
        
        btn_layout1.addWidget(self.btn_fit_spaxel)
        btn_layout1.addWidget(self.btn_fit_adj)
        btn_layout1.addWidget(self.btn_save_res)
        
        btn_layout2 = QHBoxLayout()
        self.btn_fit_map = QPushButton("Fit Entire Map")
        self.btn_fit_map.clicked.connect(self.fit_entire_map)
        
        self.btn_fit_adj_all = QPushButton("Fit Adj All")
        self.btn_fit_adj_all.clicked.connect(self.fitadj_entire_map)
        
        btn_layout2.addWidget(self.btn_fit_map)
        btn_layout2.addWidget(self.btn_fit_adj_all)
        
        # Autoflag settings
        flag_layout = QHBoxLayout()
        self.btn_autoflag = QPushButton("Autoflag")
        self.btn_autoflag.clicked.connect(self.autoflag_all)
        
        self.mask_sn_thresh_sp = QDoubleSpinBox()
        self.mask_sn_thresh_sp.setRange(0, 100)
        self.mask_sn_thresh_sp.setValue(self.state.mask_sn_thresh)
        self.mask_sn_thresh_sp.valueChanged.connect(self.on_autoflag_thresh_changed)
        
        self.mask_maxvelerr_sp = QDoubleSpinBox()
        self.mask_maxvelerr_sp.setRange(0, 1000)
        self.mask_maxvelerr_sp.setValue(self.state.mask_maxvelerr)
        self.mask_maxvelerr_sp.valueChanged.connect(self.on_autoflag_thresh_changed)
        
        flag_layout.addWidget(self.btn_autoflag)
        flag_layout.addWidget(QLabel("SNR:"))
        flag_layout.addWidget(self.mask_sn_thresh_sp)
        flag_layout.addWidget(QLabel("Max Vel/Sig Err:"))
        flag_layout.addWidget(self.mask_maxvelerr_sp)
        
        controls_grid.addLayout(flag_layout, 4, 0, 1, 4)
        controls_grid.addLayout(btn_layout1, 5, 0, 1, 4)
        controls_grid.addLayout(btn_layout2, 6, 0, 1, 4)
        
        right_layout.addWidget(controls_group, stretch=1)
        splitter.addWidget(right_pane)
        
        # Connect signal
        self.spaxel_viewer.spaxel_clicked.connect(self.on_spaxel_clicked)
        self.spaxel_viewer.slice_changed.connect(self.on_slice_changed)
        
        # Adjust initial sizes
        splitter.setSizes([950, 450])
        splitter.setStretchFactor(0, 2)
        splitter.setStretchFactor(1, 1)
        
        main_layout.addWidget(splitter)
        
    def on_autoflag_thresh_changed(self):
        self.state.mask_sn_thresh = self.mask_sn_thresh_sp.value()
        self.state.mask_maxvelerr = self.mask_maxvelerr_sp.value()
        if hasattr(self.state, 'col'):
            self.update_fit_results_label(self.state.col, self.state.row)
            
    def autoflag_all(self):
        if self.rescube is None:
            QMessageBox.warning(self, "Warning", "No map fit to autoflag.")
            return
            
        # Re-trigger slice drawing which now applies the autoflag mask dynamically
        if self.viewer_mode == 'fit_results':
            z = self.spaxel_viewer.fit_layer_cb.currentIndex()
            self.on_slice_changed(z)
            
        QMessageBox.information(self, "Autoflag", "Autoflag thresholds applied and UI updated.")

    def on_cont_mode_changed(self, index):
        mode = self.cont_mode_cb.itemData(index)
        enabled = (mode == 0)
        self.cont_minoff_sp.setEnabled(enabled)
        self.cont_maxoff_sp.setEnabled(enabled)
        self.cont_minperc_sp.setEnabled(enabled)
        self.cont_maxperc_sp.setEnabled(enabled)
        
    def update_fit_results_label(self, x: int, y: int):
        if self.rescube is None or not np.isfinite(self.rescube[0, y, x]):
            self.fit_res_label.setText("Fit Results: None")
            return
            
        vel = self.rescube[34, y, x]
        vel_err = self.rescube[35, y, x]
        sigma = self.rescube[36, y, x]
        sigma_err = self.rescube[37, y, x]
        
        is_ok = check_autoflag(self.rescube, x, y, self.state.mask_sn_thresh, self.state.mask_maxvelerr)
        status_text = "<span style='color:green;'>OK</span>" if is_ok else "<span style='color:red;'>BAD</span>"
        
        text_lines = [
            f"Velocity: {vel:7.2f} ± {vel_err:5.2f} km/s",
            f"Sigma:    {sigma:6.2f} ± {sigma_err:5.2f} km/s"
        ]
        
        from pykubeviz.core.lines import LINES_DB
        from collections import defaultdict
        
        sets = defaultdict(list)
        for i, line in enumerate(LINES_DB):
            sets[line['set']].append((i, line))
            
        for s_idx in sorted(sets.keys()):
            set_str = []
            for i, line in sets[s_idx]:
                f = self.rescube[i*2, y, x]
                f_err = self.rescube[i*2+1, y, x]
                if np.isfinite(f) and f > 0:
                    set_str.append(f"{line['fancy_name']}: {f:8.2e} ± {f_err:8.2e}")
            if set_str:
                text_lines.append(" | ".join(set_str))
                
        text_lines.append(f"Autoflag Status: <b>{status_text}</b>")
        
        self.fit_res_label.setText("<br>".join(text_lines))

    def update_coord_label(self):
        x, y = self.state.col, self.state.row
        
        if self.viewer_mode in ['data', 'variance']:
            z = self.spaxel_viewer.slice_slider.value()
            wave_val = self.state.wave[z] if (self.state.wave is not None and 0 <= z < len(self.state.wave)) else 0.0
            self.spectrum_viewer.set_current_wavelength(wave_val)
            coord_text = f"X: {x}  Y: {y}  Z: {z}"
        else:
            z = self.spaxel_viewer.fit_layer_cb.currentIndex()
            layer_name = self.spaxel_viewer.fit_layer_cb.itemText(z)
            coord_text = f"X: {x}  Y: {y}  Layer: {layer_name}"
            
        ra_str, dec_str = "--", "--"
        if self.state.wcs is not None:
            try:
                coords = self.state.wcs.pixel_to_world_values(x, y, 0)
                if len(coords) >= 2:
                    ra_str = f"{float(coords[0]):.5f}"
                    dec_str = f"{float(coords[1]):.5f}"
            except Exception:
                pass
                
        if self.viewer_mode in ['data', 'variance']:
            coord_text += f" | RA: {ra_str}  DEC: {dec_str} | Wave: {wave_val:.2f} Å"
        else:
            coord_text += f" | RA: {ra_str}  DEC: {dec_str}"
            
        self.spaxel_viewer.coord_label.setText(coord_text)
        
    def on_wavelength_clicked(self, wave_val: float):
        if self.state.wave is None or self.viewer_mode == 'fit_results':
            return
        
        # Find closest wavelength index
        idx = (np.abs(self.state.wave - wave_val)).argmin()
        self.spaxel_viewer.slice_slider.setValue(int(idx))
        
    def on_spaxel_clicked(self, x: int, y: int):
        self.state.col = x
        self.state.row = y
        self.spaxel_viewer.set_crosshair(x, y)
        self.update_coord_label()
        
        # Extract spectrum
        if self.state.datacube is not None:
            # Check bounds
            if 0 <= x < self.state.ncol and 0 <= y < self.state.nrow:
                flux = self.state.datacube[:, y, x]
                self.spectrum_viewer.set_spectrum(self.state.wave, flux)
                self.spectrum_viewer.plot_widget.setTitle(f"Spectrum at ({x}, {y})")
                
                # Plot fit if exists
                if self.rescube is not None and np.isfinite(self.rescube[0, y, x]):
                    fit_window = self.fit_window_sp.value()
                    z_base = self.z_input.value()
                    best_fit = reconstruct_fit_spectrum(self.state.wave, self.rescube[:, y, x], z_base, fit_window, instrres_extpoly=self.state.instrres_extpoly)
                    self.spectrum_viewer.set_fit(self.state.wave, best_fit)
                else:
                    self.spectrum_viewer.set_fit(self.state.wave, np.array([])) # clear
                    
                self.update_fit_results_label(x, y)
        
    def on_slice_changed(self, z: int):
        if self.state.datacube is None:
            return
            
        if self.viewer_mode == 'data':
            if 0 <= z < self.state.nwpix:
                image = self.state.datacube[z, :, :]
                self.spaxel_viewer.set_image(image)
                wave_val = self.state.wave[z] if self.state.wave is not None else 0
                self.spaxel_viewer.plot_widget.setTitle(f"Spaxel Viewer - Data (Slice {z}, λ={wave_val:.2f}Å)")
        elif self.viewer_mode == 'variance':
            if 0 <= z < self.state.nwpix and self.state.noisecube is not None:
                image = self.state.noisecube[z, :, :] ** 2
                self.spaxel_viewer.set_image(image)
                wave_val = self.state.wave[z] if self.state.wave is not None else 0
                self.spaxel_viewer.plot_widget.setTitle(f"Spaxel Viewer - Variance (Slice {z}, λ={wave_val:.2f}Å)")
        elif self.viewer_mode == 'fit_results':
            if self.rescube is not None and 0 <= z < 56:
                image = self.rescube[z, :, :].copy()
                
                # Dynamically mask out bad spaxels
                for y in range(self.state.nrow):
                    for x in range(self.state.ncol):
                        if not check_autoflag(self.rescube, x, y, self.state.mask_sn_thresh, self.state.mask_maxvelerr):
                            image[y, x] = np.nan
                            
                self.spaxel_viewer.set_image(image)
                layer_name = self.spaxel_viewer.fit_layer_cb.itemText(z)
                self.spaxel_viewer.plot_widget.setTitle(f"Spaxel Viewer - Fit Results ({layer_name})")
                
        self.update_coord_label()

    def set_viewer_mode(self, mode: str):
        if self.state.datacube is None:
            QMessageBox.warning(self, "Warning", "No datacube loaded.")
            return
            
        if mode == 'variance' and self.state.noisecube is None:
            QMessageBox.warning(self, "Warning", "No variance cube loaded.")
            return
            
        if mode == 'fit_results' and self.rescube is None:
            QMessageBox.warning(self, "Warning", "No fit results available yet.")
            return
            
        self.viewer_mode = mode
        if mode == 'data':
            self.spaxel_viewer.slider_stack.setCurrentIndex(0)
            self.on_slice_changed(self.spaxel_viewer.slice_slider.value())
        elif mode == 'variance':
            self.spaxel_viewer.slider_stack.setCurrentIndex(0)
            self.on_slice_changed(self.spaxel_viewer.slice_slider.value())
        elif mode == 'fit_results':
            self.spaxel_viewer.slider_stack.setCurrentIndex(1)
            self.on_slice_changed(self.spaxel_viewer.fit_layer_cb.currentIndex())
        
    def open_datacube(self):
        filepath, _ = QFileDialog.getOpenFileName(self, "Open Datacube", "", "FITS Files (*.fits *.fits.gz);;All Files (*)")
        if filepath:
            success = load_datacube(self.state, filepath)
            if success:
                # Create a median image for the spaxel viewer
                median_image = np.nanmedian(self.state.datacube, axis=0)
                self.spaxel_viewer.set_image(median_image)
                
                # Setup slider
                self.spaxel_viewer.slice_slider.setRange(0, self.state.nwpix - 1)
                self.spaxel_viewer.slice_slider.setEnabled(True)
                self.spaxel_viewer.slice_slider.setValue(self.state.nwpix // 2)
                
                self.setWindowTitle(f"Kubeviz - {self.state.filename}")
                self.rescube = None
                self.btn_save_res.setEnabled(False)
                QMessageBox.information(self, "Success", f"Loaded {self.state.filename}")
            else:
                QMessageBox.critical(self, "Error", "Failed to load the datacube.")

    def get_z_init(self):
        val = self.z_input.value()
        if val <= self.z_input.minimum():
            return None
        return val

    def fit_current_spaxel(self):
        if self.state.datacube is None:
            from PyQt6.QtWidgets import QMessageBox
            QMessageBox.warning(self, "Warning", "No datacube loaded.")
            return
            
        x, y = self.state.col, self.state.row
        flux = self.state.datacube[:, y, x]
        wave = self.state.wave
        noise = self.state.noisecube[:, y, x] if self.state.noisecube is not None else None
        
        try:
            z_init = self.get_z_init()
            is_auto_detect = (z_init is None)
            if is_auto_detect:
                import numpy as np
                max_idx = np.nanargmax(flux)
                guessed_z = wave[max_idx] / 6562.819 - 1.0
                z_init_for_fit = guessed_z
            else:
                z_init_for_fit = z_init
                
            cont_mode = self.cont_mode_cb.itemData(self.cont_mode_cb.currentIndex())
            cont_minoff = self.cont_minoff_sp.value()
            cont_maxoff = self.cont_maxoff_sp.value()
            cont_minperc = self.cont_minperc_sp.value()
            cont_maxperc = self.cont_maxperc_sp.value()
            
            fit_window = self.fit_window_sp.value()
            
            from pykubeviz.analysis.fitting import fit_all_lines
            from pykubeviz.core.lines import LINES_DB
            import numpy as np
            
            best_fit, fit_params = fit_all_lines(
                wave, flux, noise, z_init=z_init_for_fit,
                cont_mode=cont_mode, cont_minoff=cont_minoff, cont_maxoff=cont_maxoff,
                cont_minperc=cont_minperc, cont_maxperc=cont_maxperc,
                lines_config=self.state.fit_lines_config,
                instrres_extpoly=self.state.instrres_extpoly,
                error_method=self.state.error_method,
                fit_window=fit_window
            )
            
            if not fit_params:
                from PyQt6.QtWidgets import QMessageBox
                QMessageBox.warning(self, "Fit Failed", "Could not fit emission lines.")
                return
                
            if self.rescube is None:
                self.rescube = np.full((56, self.state.nrow, self.state.ncol), np.nan)
                
            for i, line in enumerate(LINES_DB):
                self.rescube[i*2, y, x] = fit_params.get(f"flux_{line['name']}", np.nan)
                self.rescube[i*2+1, y, x] = fit_params.get(f"flux_{line['name']}_err", np.nan)
                
            self.rescube[34, y, x] = fit_params.get('vel', np.nan)
            self.rescube[35, y, x] = fit_params.get('vel_err', np.nan)
            self.rescube[36, y, x] = fit_params.get('sigma_kms', np.nan)
            self.rescube[37, y, x] = fit_params.get('sigma_kms_err', np.nan)
            
            for set_id in range(1, 10):
                self.rescube[38 + (set_id-1), y, x] = fit_params.get(f'continuum_set{set_id}', np.nan)
                self.rescube[47 + (set_id-1), y, x] = fit_params.get(f'continuum_set{set_id}_err', np.nan)
            
            self.btn_save_res.setEnabled(True)

            if is_auto_detect:
                vel = fit_params.get('vel', 0.0)
                z_final = (1.0 + vel / self.state.ckms) * (1.0 + z_init_for_fit) - 1.0
                self.z_input.setValue(z_final)
                self.rescube[34, y, x] = 0.0  # Velocity is 0 relative to new systemic z
            
            self.spectrum_viewer.set_fit(wave, best_fit)
            self.update_fit_results_label(x, y)
            
            if self.viewer_mode == 'fit_results':
                z = self.spaxel_viewer.fit_layer_cb.currentIndex()
                self.on_slice_changed(z)
                
        except Exception as e:
            import traceback
            traceback.print_exc()
            from PyQt6.QtWidgets import QMessageBox
            QMessageBox.critical(self, "Fit Error", f"Failed to fit: {e}\n\n" + traceback.format_exc())

    def fitadj_current_spaxel(self):
        if self.state.datacube is None or self.rescube is None:
            from PyQt6.QtWidgets import QMessageBox
            QMessageBox.warning(self, "Warning", "No datacube loaded or map not yet fit.")
            return
            
        x, y = self.state.col, self.state.row
        flux = self.state.datacube[:, y, x]
        wave = self.state.wave
        noise = self.state.noisecube[:, y, x] if self.state.noisecube is not None else None
        
        try:
            from pykubeviz.analysis.fitting import guess_from_adjacent_spaxels, fit_all_lines
            from pykubeviz.core.lines import LINES_DB
            import numpy as np
            
            z_base = self.get_z_init() or 0.0
            z_guess, sig_guess, amp_guess = guess_from_adjacent_spaxels(
                self.rescube, x, y, 
                self.state.mask_sn_thresh, self.state.mask_maxvelerr, z_base
            )
            
            if z_guess is None:
                from PyQt6.QtWidgets import QMessageBox
                QMessageBox.warning(self, "Fit Adj Failed", "No valid adjacent spaxels found.")
                return
                
            cont_mode = self.cont_mode_cb.itemData(self.cont_mode_cb.currentIndex())
            cont_minoff = self.cont_minoff_sp.value()
            cont_maxoff = self.cont_maxoff_sp.value()
            cont_minperc = self.cont_minperc_sp.value()
            cont_maxperc = self.cont_maxperc_sp.value()
            
            fit_window = self.fit_window_sp.value()
            
            best_fit, fit_params = fit_all_lines(
                wave, flux, noise, z_init=z_guess, 
                init_sigma=sig_guess * 6562.819 / 299792.458 if sig_guess else 2.0,
                init_amp_ha=amp_guess,
                cont_mode=cont_mode, cont_minoff=cont_minoff, cont_maxoff=cont_maxoff,
                cont_minperc=cont_minperc, cont_maxperc=cont_maxperc,
                lines_config=self.state.fit_lines_config,
                instrres_extpoly=self.state.instrres_extpoly,
                error_method=self.state.error_method,
                fit_window=fit_window
            )
            
            if not fit_params:
                from PyQt6.QtWidgets import QMessageBox
                QMessageBox.warning(self, "Fit Failed", "Could not fit emission lines.")
                return
                
            for i, line in enumerate(LINES_DB):
                self.rescube[i*2, y, x] = fit_params.get(f"flux_{line['name']}", np.nan)
                self.rescube[i*2+1, y, x] = fit_params.get(f"flux_{line['name']}_err", np.nan)
                
            self.rescube[34, y, x] = fit_params.get('vel', np.nan)
            self.rescube[35, y, x] = fit_params.get('vel_err', np.nan)
            self.rescube[36, y, x] = fit_params.get('sigma_kms', np.nan)
            self.rescube[37, y, x] = fit_params.get('sigma_kms_err', np.nan)
            
            for set_id in range(1, 10):
                self.rescube[38 + (set_id-1), y, x] = fit_params.get(f'continuum_set{set_id}', np.nan)
                self.rescube[47 + (set_id-1), y, x] = fit_params.get(f'continuum_set{set_id}_err', np.nan)
            
            self.btn_save_res.setEnabled(True)
            self.spectrum_viewer.set_fit(wave, best_fit)
            self.update_fit_results_label(x, y)
            
            if self.viewer_mode == 'fit_results':
                z = self.spaxel_viewer.fit_layer_cb.currentIndex()
                self.on_slice_changed(z)
            
        except Exception as e:
            from PyQt6.QtWidgets import QMessageBox
            QMessageBox.critical(self, "Fit Error", f"Failed to fit adj: {e}")

    def fit_entire_map(self):
        if self.state.datacube is None:
            QMessageBox.warning(self, "Warning", "No datacube loaded.")
            return
            
        progress = QProgressDialog("Fitting Entire Map... (0%)", "Cancel", 0, 100, self)
        progress.setWindowModality(Qt.WindowModality.WindowModal)
        progress.setValue(0)
        progress.show()
        
        progress_state = {'last_pct': -1}
        
        def update_progress(completed, total, current_rescube):
            pct = int(completed / total * 100) if total > 0 else 0
            if pct > progress_state['last_pct']:
                progress_state['last_pct'] = pct
                
                # Periodically sync to main rescube
                self.rescube = current_rescube.copy()
                
                text = f"Fitting Entire Map... ({completed}/{total} spaxels - {pct}%)"
                progress.setLabelText(text)
                progress.setValue(pct)
                
                # Process UI events so Cancel button works
                QApplication.processEvents()
                
                if progress.wasCanceled():
                    raise Exception("Cancelled")

        try:
            z_init = self.get_z_init()
            cont_mode = self.cont_mode_cb.itemData(self.cont_mode_cb.currentIndex())
            cont_minoff = self.cont_minoff_sp.value()
            cont_maxoff = self.cont_maxoff_sp.value()
            cont_minperc = self.cont_minperc_sp.value()
            cont_maxperc = self.cont_maxperc_sp.value()
            fit_window = self.fit_window_sp.value()
            
            if z_init is None:
                # Auto-detect on central spaxel to lock in a global redshift
                cx, cy = self.state.ncol // 2, self.state.nrow // 2
                flux_center = self.state.datacube[:, cy, cx]
                if not np.any(np.isfinite(flux_center)):
                    for y_idx in range(self.state.nrow):
                        for x_idx in range(self.state.ncol):
                            if np.any(np.isfinite(self.state.datacube[:, y_idx, x_idx])):
                                flux_center = self.state.datacube[:, y_idx, x_idx]
                                break
                        if np.any(np.isfinite(flux_center)): break
                
                max_idx = np.nanargmax(flux_center)
                guessed_z = self.state.wave[max_idx] / 6562.819 - 1.0
                
                _, fit_params = fit_all_lines(
                    self.state.wave, flux_center, None, z_init=guessed_z,
                    cont_mode=cont_mode, cont_minoff=cont_minoff, cont_maxoff=cont_maxoff,
                    cont_minperc=cont_minperc, cont_maxperc=cont_maxperc,
                    lines_config=self.state.fit_lines_config,
                    instrres_extpoly=self.state.instrres_extpoly,
                    fit_window=fit_window
                )
                if fit_params:
                    vel = fit_params['vel']
                    z_init = (1.0 + vel / self.state.ckms) * (1.0 + guessed_z) - 1.0
                else:
                    z_init = guessed_z
                
                self.z_input.setValue(z_init)
            self.rescube = fit_map(
                self.state.datacube, 
                self.state.wave, 
                self.state.noisecube, 
                z_init=z_init,
                progress_callback=update_progress,
                cont_mode=cont_mode, cont_minoff=cont_minoff, cont_maxoff=cont_maxoff,
                cont_minperc=cont_minperc, cont_maxperc=cont_maxperc,
                lines_config=self.state.fit_lines_config,
                instrres_extpoly=self.state.instrres_extpoly,
                error_method=self.state.error_method,
                fit_window=fit_window
            )
            
            self.btn_save_res.setEnabled(True)
            progress.setValue(100)
            QMessageBox.information(self, "Success", "Map fitting completed successfully.")
        except Exception as e:
            progress.cancel()
            if str(e) == "Cancelled":
                QMessageBox.information(self, "Cancelled", "Map fitting was cancelled by the user. Partial results have been retained.")
                self.btn_save_res.setEnabled(True)
            else:
                QMessageBox.critical(self, "Fit Error", f"Failed to fit map: {e}")

    def fitadj_entire_map(self):
        if self.state.datacube is None or self.rescube is None:
            QMessageBox.warning(self, "Warning", "No datacube loaded or map not yet fit.")
            return
            
        import time
        progress = QProgressDialog("Iterative Fit Adj All...", "Cancel", 0, 0, self)
        progress.setWindowModality(Qt.WindowModality.WindowModal)
        progress.show()
        
        progress_state = {'last_update': 0.0}
        
        def update_progress(it, fixed, total_bad, current_rescube=None):
            QApplication.processEvents()
            if progress.wasCanceled():
                raise Exception("Fit cancelled by user.")
                
            current_time = time.time()
            if current_time - progress_state['last_update'] >= 60.0 or progress_state['last_update'] == 0.0:
                progress_state['last_update'] = current_time
                progress.setLabelText(f"Fit has been improved for {fixed} spaxels so far...")
                if current_rescube is not None:
                    self.rescube = current_rescube.copy()
            
        try:
            cont_mode = self.cont_mode_cb.itemData(self.cont_mode_cb.currentIndex())
            cont_minoff = self.cont_minoff_sp.value()
            cont_maxoff = self.cont_maxoff_sp.value()
            cont_minperc = self.cont_minperc_sp.value()
            cont_maxperc = self.cont_maxperc_sp.value()
            
            kwargs = {
                'cont_mode': cont_mode,
                'cont_minoff': cont_minoff,
                'cont_maxoff': cont_maxoff,
                'cont_minperc': cont_minperc,
                'cont_maxperc': cont_maxperc,
                'lines_config': self.state.fit_lines_config,
                'instrres_extpoly': self.state.instrres_extpoly,
                'error_method': self.state.error_method
            }
            
            z_base = self.get_z_init() or 0.0
            self.rescube, total_fixed = fit_adj_all(
                self.state.datacube, self.state.wave, self.state.noisecube, self.rescube,
                kwargs=kwargs, sn_thresh=self.state.mask_sn_thresh, maxvelerr=self.state.mask_maxvelerr,
                progress_callback=update_progress, max_iter=100000, z_base=z_base
            )
            
            progress.close()
            QMessageBox.information(self, "Success", f"Fit Adj All completed. Fixed {total_fixed} spaxels.")
            
            # update UI if spaxel is clicked
            x, y = self.state.col, self.state.row
            self.update_fit_results_label(x, y)
            if np.isfinite(self.rescube[0, y, x]):
                fit_window = self.fit_window_sp.value()
                z_base = self.z_input.value()
                best_fit = reconstruct_fit_spectrum(self.state.wave, self.rescube[:, y, x], z_base, fit_window, instrres_extpoly=self.state.instrres_extpoly)
                self.spectrum_viewer.set_fit(self.state.wave, best_fit)
                
        except Exception as e:
            progress.cancel()
            QMessageBox.critical(self, "Fit Error", f"Failed to fit adj all: {e}")

    def save_results(self):
        if self.rescube is None:
            return
            
        default_name = f"fit_{self.state.filename}" if self.state.filename else "fit_results.fits"
        filepath, _ = QFileDialog.getSaveFileName(self, "Save Fit Results", default_name, "FITS Files (*.fits *.fits.gz)")
        if filepath:
            success = save_fit_results(filepath, self.rescube, self.z_input.value(), self.state.indatahead)
            if success:
                QMessageBox.information(self, "Success", "Fit results saved successfully.")
            else:
                QMessageBox.critical(self, "Error", "Failed to save fit results.")

    def load_fit_results_action(self):
        if self.state.datacube is None:
            QMessageBox.warning(self, "Warning", "Please load a datacube first.")
            return
            
        filepath, _ = QFileDialog.getOpenFileName(self, "Load Fit Results", "", "FITS Files (*.fits *.fits.gz)")
        if filepath:
            from pykubeviz.io.data_io import load_fit_results
            try:
                rescube, z_base = load_fit_results(filepath)
                if rescube.shape[0] < 56 or rescube.shape[1:] != (self.state.nrow, self.state.ncol):
                    QMessageBox.warning(self, "Warning", "Loaded cube does not match expected dimensions or layers.")
                    return
                self.rescube = rescube
                self.z_input.setValue(z_base)
                self.btn_save_res.setEnabled(True)
                
                # Check autoflag status on the whole map if we want? Or just update UI.
                if self.viewer_mode == 'fit_results':
                    z = self.spaxel_viewer.fit_layer_cb.currentIndex()
                    self.on_slice_changed(z)
                    
                QMessageBox.information(self, "Success", "Fit results loaded successfully.")
            except Exception as e:
                import traceback
                traceback.print_exc()
                QMessageBox.critical(self, "Error", f"Failed to load fit results: {e}")

    def set_error_method(self, method: str):
        self.state.error_method = method

    def open_fit_lines_dialog(self):
        dlg = FitLinesDialog(self.state, self)
        dlg.exec()

    def keyPressEvent(self, event):
        if event.key() == Qt.Key.Key_Z and self.spectrum_viewer.underMouse():
            self.trigger_redshift_identification()
        super().keyPressEvent(event)
        
    def trigger_redshift_identification(self):
        if self.state.wave is None:
            return
            
        # Get the current wavelength from the slider
        idx = self.spaxel_viewer.slice_slider.value()
        if not (0 <= idx < len(self.state.wave)):
            return
            
        current_wave = self.state.wave[idx]
        
        lines = {
            'H-alpha (6562.819)': 6562.819,
            '[NII] red (6583.450)': 6583.450,
            '[NII] blue (6548.050)': 6548.050,
            '[OIII] (5006.843)': 5006.843,
            'H-beta (4861.325)': 4861.325,
            '[SII] (6716.440)': 6716.440,
            '[SII] (6730.810)': 6730.810
        }
        
        item, ok = QInputDialog.getItem(
            self, "Identify Line", 
            f"Current Wavelength: {current_wave:.2f} Å\n\nWhich rest-frame line is this?",
            list(lines.keys()), 0, False
        )
        
        if ok and item:
            rest_wave = lines[item]
            redshift = (current_wave / rest_wave) - 1.0
            
            # Update the spin box
            self.z_input.setValue(redshift)
            QMessageBox.information(self, "Redshift Computed", f"Computed Redshift (z) = {redshift:.5f}\n(Based on {item})")

def main():
    app = QApplication(sys.argv)
    
    # Force US locale so that numerical inputs use '.' as decimal separator
    from PyQt6.QtCore import QLocale
    QLocale.setDefault(QLocale(QLocale.Language.English, QLocale.Country.UnitedStates))
    
    state = KubevizState()
    window = KubevizMainWindow(state)
    window.show()
    sys.exit(app.exec())

if __name__ == '__main__':
    main()
