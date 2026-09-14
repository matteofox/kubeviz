"""GUI smoke tests, run with the offscreen Qt platform."""
import os

import numpy as np
import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
from kubeviz.gui.qtenv import prepare_qt_environment  # noqa: E402

prepare_qt_environment()
pytest.importorskip("PyQt6.QtWidgets")
pytest.importorskip("pyqtgraph")

from PyQt6.QtCore import QEvent, Qt  # noqa: E402
from PyQt6.QtGui import QKeyEvent  # noqa: E402
from PyQt6.QtWidgets import QApplication  # noqa: E402

from kubeviz import constants as C  # noqa: E402
from kubeviz.core.display import current_image, linefit_image_update  # noqa: E402
from kubeviz.gui.controller import UPDATE_FULL, KubevizGUI  # noqa: E402
from kubeviz.gui.fitoverlay import fit_overlays  # noqa: E402
from kubeviz.gui.scaling import (ZCUT_HISTEQ, ZCUT_MINMAX, ZCUT_USER_LOG, ZCUT_ZSCALE, colour_table,  # noqa: E402
                                 render_rgb, scale_image)
from kubeviz.session import start_session  # noqa: E402


@pytest.fixture(scope="module")
def qapp():
    app = QApplication.instance() or QApplication([])
    return app


@pytest.fixture
def gui(qapp, synth_cube):
    fname, truth = synth_cube
    state = start_session(fname, redshift=truth["redshift"], lineset=1)
    g = KubevizGUI(state, {})
    g.show()
    g.linefit.show()
    qapp.processEvents()
    yield g
    g.linefit.close()
    g.close()


def _key(text, key=Qt.Key.Key_unknown):
    return QKeyEvent(QEvent.Type.KeyPress, key, Qt.KeyboardModifier.NoModifier, text)


def test_scaling_functions():
    img = np.random.default_rng(0).normal(size=(20, 30))
    img[0, 0] = np.nan
    for zc in (ZCUT_HISTEQ, ZCUT_MINMAX, ZCUT_ZSCALE):
        s, rng = scale_image(img, zc, 0, 1)
        assert np.isnan(s[0, 0]) and np.nanmin(s) >= 0 and np.nanmax(s) <= 1
        assert rng[1] >= rng[0]
    s, _ = scale_image(np.abs(img) + 0.1, ZCUT_USER_LOG, 0, 3)
    assert np.nanmax(s) <= 1
    rgb = render_rgb(s, colour_table(0))
    assert rgb.shape == (20, 30, 3) and rgb.dtype == np.uint8
    assert tuple(rgb[0, 0]) == (255, 255, 255)          # NaN -> white


def test_gui_builds_and_updates(gui, qapp):
    st = gui.state
    assert gui.spax.image.image.shape == (st.Nrow, st.Ncol, 3)
    assert len(gui.linefit.w) == 2 * 2 + 3 * st.Nlines
    # move crosshair with the keyboard and mouse
    c0 = st.col
    gui.keyboard(_key("", Qt.Key.Key_Right), "spax")
    assert st.col == c0 + 1
    gui.on_spaxel_pressed(3, 4)
    gui.on_spaxel_released()
    assert (st.col, st.row) == (3, 4)
    gui.keyboard(_key("."), "spax")
    assert st.wpix == st.Nwpix // 2 + 1
    # image modes and cube selection through the menu codes
    st.wavrange1[:] = [540, 570]
    gui.menu_action("Sum1")
    assert st.imgmode == C.IMG_SUM1 and st.img1.max() > 0
    gui.menu_action("SN")
    assert st.cubesel == C.CUBE_SN
    gui.menu_action("Zscale")
    gui.menu_action("Heat")
    gui.menu_action("Invert")
    assert st.invert == 1
    gui.on_contrast(0.3, 0.7)
    gui.menu_action("Data")
    gui.menu_action("Slice")
    qapp.processEvents()


def test_gui_fit_and_result_maps(gui, qapp):
    st = gui.state
    lf = gui.linefit
    # tick "fit" for the narrow components and the continuum via the linefit codes
    for par in range(3, 3 + st.Nlines):
        lf.ctl.linefit_action(f"FITN{par}")
    lf.ctl.linefit_action("FITC3")
    assert np.all(st.pndofit[:st.Nlines] == 1) and np.all(st.pcdofit[:st.Nlines] == 1)
    st.col, st.row = 6, 5
    lf.ctl.linefit_action("FIT")
    rs = st.res
    assert rs.n[5, 6, 3] > 0
    assert "Not Fit" not in lf.w[("N", 3)]["best"].text()
    # show the fit in the zoom window
    st.zoommap = True
    gui.zoom.setVisible(True)
    lf.ctl.linefit_action("SHOWN3")
    lf.ctl.linefit_action("SHOWC3")
    st.wpix = int(np.argmin(np.abs(st.wave - st.lines[0])))
    gui.update_all(UPDATE_FULL)
    ov = fit_overlays(st, st.wave[st.wpix - 32: st.wpix + 33])
    kinds = {o.kind for o in ov}
    assert "narrow" in kinds and "total" in kinds
    ov = fit_overlays(st, st.wave[st.wpix - 400: st.wpix + 100])       # continuum windows visible
    assert "contline" in {o.kind for o in ov}
    # result maps
    lf.ctl.linefit_action("IMAGEN3")
    assert st.cubesel == C.CUBE_LINEFIT and st.par_imagebutton == "N3"
    img, bad = current_image(st)
    assert np.isfinite(img[5, 6]) and img[0, 1] == 0            # unfit spaxels are 0 until autoflag
    lf.ctl.linefit_action("FLAG")
    assert rs.n[5, 6, 0] == C.FLAG_MANUAL
    lf.ctl.linefit_action("FLAG")
    assert rs.n[5, 6, 0] == 0
    # fit all with the progress dialog (no user interaction needed offscreen)
    lf.ctl.linefit_action("FITALL")
    assert np.sum(rs.n[..., 3] > 0) > 50
    assert linefit_image_update(st, "N1")
    gui.menu_action("LineErrors")
    gui.menu_action("LineSN")
    gui.update_all(UPDATE_FULL)
    # switch to moments and back
    lf.ctl.linefit_action("TYPE")
    assert st.linefit_type == C.FIT_MOMENTS
    lf.ctl.linefit_action("TYPE")
    # redshift change rebuilds the linefit window
    gui.linefit_action("REDSHIFT", 0.0205)
    assert abs(st.redshift - 0.0205) < 1e-9
    assert gui.linefit.state is st


def test_gui_masks_and_spectrum_modes(gui, qapp):
    st = gui.state
    gui.menu_action("Select")
    assert st.cursormode == 2 and st.linefit_mode == C.MODE_MASK
    st.maskmode, st.maskradius = 1, 2
    gui.on_spaxel_pressed(6, 5)
    gui.on_spaxel_released()
    assert st.spaxselect[0].sum() > 5
    assert st.medspec.max() > 0
    gui._set_specmode(C.SPEC_MEDIAN)
    gui._set_specmode(C.SPEC_WAVG)
    gui.menu_action("NewMask")
    assert st.Nmask == 2
    gui.menu_action("PrevMask")
    assert st.imask == 1
    gui.keyboard(_key("s"), "spec")
    assert st.npress == 1
    st.wpix += 20
    gui.keyboard(_key("s"), "spec")
    assert st.npress == 0 and st.wavrange1[1] - st.wavrange1[0] == 20
    gui.keyboard(_key("g"), "zoom")
    gui.menu_action("Crosshair")
    assert st.cursormode == 1 and st.linefit_mode == C.MODE_SPAXEL
