import pyqtgraph as pg
from PyQt6.QtWidgets import QWidget, QVBoxLayout, QSlider, QLabel, QHBoxLayout, QStackedWidget, QComboBox, QSizePolicy
from PyQt6.QtCore import pyqtSignal, Qt
import numpy as np

class SpaxelViewer(QWidget):
    """
    2D Viewer for the datacube (spatially).
    Displays a selected slice or moment map.
    """
    spaxel_clicked = pyqtSignal(int, int)
    slice_changed = pyqtSignal(int)
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self.layout = QVBoxLayout(self)
        self.layout.setContentsMargins(0, 0, 0, 0)
        
        self.plot_widget = pg.PlotWidget(title="Spaxel Viewer")
        self.plot_widget.setAspectLocked(True, ratio=1.0)
        self.image_item = pg.ImageItem()
        self.plot_widget.addItem(self.image_item)
        
        # Add a crosshair or ROI for selecting spaxels
        blue_pen = pg.mkPen('b', width=1)
        self.v_line = pg.InfiniteLine(angle=90, movable=False, pen=blue_pen)
        self.h_line = pg.InfiniteLine(angle=0, movable=False, pen=blue_pen)
        self.plot_widget.addItem(self.v_line, ignoreBounds=True)
        self.plot_widget.addItem(self.h_line, ignoreBounds=True)
        
        # Connect mouse clicks on the image
        self.image_item.mouseClickEvent = self.on_image_clicked
        
        self.layout.addWidget(self.plot_widget)
        
        # Add slider area
        self.slider_stack = QStackedWidget()
        self.slider_stack.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Maximum)
        
        # Page 0: Wavelength Slider
        page_wave = QWidget()
        page_wave_layout = QHBoxLayout(page_wave)
        page_wave_layout.setContentsMargins(0, 0, 0, 0)
        self.wave_label = QLabel("Wavelength Slice:")
        page_wave_layout.addWidget(self.wave_label)
        self.slice_slider = QSlider(Qt.Orientation.Horizontal)
        self.slice_slider.setEnabled(False)
        self.slice_slider.valueChanged.connect(self._on_slider_value_changed)
        page_wave_layout.addWidget(self.slice_slider)
        
        # Page 1: Fit Results Combo
        page_fit = QWidget()
        page_fit_layout = QHBoxLayout(page_fit)
        page_fit_layout.setContentsMargins(0, 0, 0, 0)
        page_fit_layout.addWidget(QLabel("Fit Result Layer:"))
        self.fit_layer_cb = QComboBox()
        from pykubeviz.core.lines import LINES_DB
        items = []
        for line in LINES_DB:
            items.append(f"Flux {line['fancy_name']}")
            items.append(f"Flux {line['fancy_name']} Error")
            
        items.extend([
            "Velocity", "Velocity Error",
            "Sigma", "Sigma Error"
        ])
        
        for i in range(1, 10):
            items.append(f"Continuum Set {i}")
        for i in range(1, 10):
            items.append(f"Continuum Set {i} Error")
            
        self.fit_layer_cb.addItems(items)
        self.fit_layer_cb.currentIndexChanged.connect(self._on_slider_value_changed)
        page_fit_layout.addWidget(self.fit_layer_cb)
        
        self.slider_stack.addWidget(page_wave)
        self.slider_stack.addWidget(page_fit)
        
        self.layout.addWidget(self.slider_stack)
        
        # Add coordinate label
        self.coord_label = QLabel("X: --  Y: -- | RA: -- DEC: --")
        self.layout.addWidget(self.coord_label)
        
    def _on_slider_value_changed(self, value):
        self.slice_changed.emit(value)
        
    def on_image_clicked(self, event):
        if event.button() == Qt.MouseButton.LeftButton:
            pos = event.pos()
            x, y = int(pos.x()), int(pos.y())
            # Emit signal with coordinates
            self.spaxel_clicked.emit(x, y)
            
    def set_image(self, image: np.ndarray):
        """ Update the displayed 2D image. """
        self.image_item.setImage(image.T, autoLevels=True)
        
    def set_crosshair(self, x: int, y: int):
        """ Update the crosshair position. """
        self.v_line.setPos(x + 0.5)
        self.h_line.setPos(y + 0.5)
        
    def set_colormap(self, cmap_name: str):
        """ Update the colormap of the 2D image. """
        try:
            cmap = pg.colormap.get(cmap_name)
            self.image_item.setColorMap(cmap)
        except Exception as e:
            print(f"Failed to set colormap {cmap_name}: {e}")

class SpectrumViewer(QWidget):
    """
    1D Viewer for the extracted spectrum.
    """
    wavelength_clicked = pyqtSignal(float)
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self.layout = QVBoxLayout(self)
        self.layout.setContentsMargins(0, 0, 0, 0)
        
        self.plot_widget = pg.PlotWidget(title="Spectrum Viewer")
        self.plot_widget.setLabel('bottom', 'Wavelength')
        self.plot_widget.setLabel('left', 'Flux')
        
        self.spectrum_item = self.plot_widget.plot(pen='y')
        self.fit_item = self.plot_widget.plot(pen='r', connect='finite') # for showing fits
        
        self.wave_line = pg.InfiniteLine(angle=90, movable=False, pen=pg.mkPen('b', width=1))
        self.plot_widget.addItem(self.wave_line)
        self.wave_line.setVisible(False)
        
        self.layout.addWidget(self.plot_widget)
        self.plot_widget.scene().sigMouseClicked.connect(self.on_mouse_clicked)
        
    def on_mouse_clicked(self, event):
        if event.button() == Qt.MouseButton.LeftButton:
            pos = event.scenePos()
            if self.plot_widget.sceneBoundingRect().contains(pos):
                mouse_point = self.plot_widget.plotItem.vb.mapSceneToView(pos)
                self.wavelength_clicked.emit(mouse_point.x())
        
    def set_spectrum(self, wave: np.ndarray, flux: np.ndarray):
        """ Update the 1D spectrum plot. """
        self.spectrum_item.setData(wave, flux)
        
    def set_fit(self, wave: np.ndarray, fit_flux: np.ndarray):
        """ Update the fit overlay. """
        self.fit_item.setData(wave, fit_flux)
        
    def set_current_wavelength(self, wave_val: float):
        """ Update the vertical wavelength indicator line. """
        self.wave_line.setPos(wave_val)
        self.wave_line.setVisible(True)
