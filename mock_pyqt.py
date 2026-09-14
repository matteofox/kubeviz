import sys

class MockModule:
    def __getattr__(self, name):
        return MockClass

class MockClass:
    def __init__(self, *args, **kwargs):
        pass
    def __call__(self, *args, **kwargs):
        return MockClass()
    def __getattr__(self, name):
        return MockClass()

sys.modules['PyQt6'] = MockModule()
sys.modules['PyQt6.QtWidgets'] = MockModule()
sys.modules['PyQt6.QtCore'] = MockModule()
sys.modules['PyQt6.QtGui'] = MockModule()
sys.modules['pyqtgraph'] = MockModule()
sys.modules['lmfit'] = MockModule()
sys.modules['lmfit.models'] = MockModule()

import numpy as np

class MockModel:
    def __init__(self, prefix=''):
        self.prefix = prefix
    def __add__(self, other):
        return MockModel(self.prefix)
    def make_params(self):
        return MockParams()
    def fit(self, *args, **kwargs):
        res = MockClass()
        res.best_fit = np.zeros(100)
        res.params = MockParams()
        return res
        
class MockParams:
    def __init__(self):
        pass
    def add(self, *args, **kwargs):
        pass
    def __getitem__(self, key):
        m = MockClass()
        m.set = lambda *a, **k: None
        m.value = 1.0
        m.stderr = 0.1
        return m
    def __setitem__(self, key, value):
        pass

sys.modules['lmfit'].Model = MockModel
sys.modules['lmfit'].Parameters = MockParams
sys.modules['lmfit.models'].GaussianModel = MockModel
sys.modules['lmfit.models'].ConstantModel = MockModel

from pykubeviz.ui.main_window import KubeVizMainWindow

win = KubeVizMainWindow()
# Force state
win.state.datacube = np.zeros((100, 10, 10))
win.state.noisecube = np.zeros((100, 10, 10))
flux = np.random.rand(100)
noise = np.random.rand(100)
win.state.datacube[:, 5, 5] = flux
win.state.noisecube[:, 5, 5] = noise
win.state.wave = np.linspace(4000, 7000, 100)
win.state.col, win.state.row = 5, 5
win.state.nrow, win.state.ncol = 10, 10
win.rescube = np.zeros((10, 10, 10))

# Mock widgets that have value() or currentIndex()
win.get_z_init = lambda: 0.0
win.cont_mode_cb = MockClass()
win.cont_mode_cb.itemData = lambda x: 0
win.cont_minoff_sp = MockClass()
win.cont_minoff_sp.value = lambda: 200.0
win.cont_maxoff_sp = MockClass()
win.cont_maxoff_sp.value = lambda: 500.0
win.cont_minperc_sp = MockClass()
win.cont_minperc_sp.value = lambda: 40.0
win.cont_maxperc_sp = MockClass()
win.cont_maxperc_sp.value = lambda: 60.0
win.spectrum_viewer = MockClass()
win.fit_res_label = MockClass()

win.fit_current_spaxel()
