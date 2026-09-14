import os

import numpy as np
from astropy.io import fits

from kubeviz import constants as C
from kubeviz.core.extraction import medianspec, medsum_image_update
from kubeviz.core.masks import newmask, select_spaxels
from kubeviz.fitting.fitall import fitadjall, fitall
from kubeviz.fitting.flags import autoflag
from kubeviz.fitting.linefit import dofit, getredshift
from kubeviz.io.cube import decode_trimstr
from kubeviz.io.results import loadres, saveres
from kubeviz.io.session import load_session, save_session
from kubeviz.session import batchmode, start_session


def _session(synth_cube, **kw):
    fname, truth = synth_cube
    opts = dict(redshift=truth["redshift"], lineset=1)
    opts.update(kw)
    state = start_session(fname, **opts)
    return state, truth


def test_load_cube_basic(synth_cube):
    state, truth = _session(synth_cube)
    assert state.instr == "muse"
    assert state.noiseisvar
    assert state.shape == (800, 10, 12)
    np.testing.assert_allclose(state.wave, truth["wave"])
    assert state.badpixelimg[0, 0]          # NaN corner spaxel is flagged bad
    assert not state.badpixelimg[5, 5]
    assert state.Nlines == 3
    assert state.linenames == ["Ha", "n2_b", "n2_r"]
    assert state.lineset_ind[0].tolist() == [0, 1, 2]
    assert abs(state.mainline() - 6562.819 * 1.02) < 1e-6
    assert state.instrres_mode == C.INSTRRES_EXTPOLY
    assert 2000 < state.getinstrres(6700.0) < 4000


def test_trim_and_smooth(synth_cube):
    fname, truth = synth_cube
    P, start = decode_trimstr("[3:8,2:7,*]", (12, 10, 800))
    assert P == [2, 7, 1, 6, 0, 799] and start == (2, 1, 0)
    assert decode_trimstr("[0:8,2:7,*]", (12, 10, 800))[0] is None
    state = start_session(fname, redshift=truth["redshift"], lineset=1, trim="[3:8,2:7,101:300]", smooth=3)
    assert state.shape == (200, 6, 6)
    assert state.Startcol == 2 and state.Startrow == 1 and state.Startwpix == 100
    np.testing.assert_allclose(state.wave, truth["wave"][100:300])
    assert np.all(np.isfinite(state.datacube[:, 3, 3]))


def test_single_spaxel_fit_recovers_truth(synth_cube):
    state, truth = _session(synth_cube)
    state.pndofit[:3] = 1
    state.pcdofit[:3] = 1
    col, row = 6, 5
    state.col, state.row = col, row
    ctxs = dofit(state)
    assert len(ctxs) == 1
    ctx = ctxs[0]
    assert ctx.errmsg == ""
    rs = state.res
    n = rs.n[row, col]
    nerr = rs.nerr[row, col]
    assert n[0] == 0
    # dv, sigma, fluxes within ~3 sigma of the truth
    assert abs(n[1] - truth["vel"][row, col]) < 3 * abs(nerr[1, 0]) + 2.0
    assert abs(n[2] - truth["sigma"][row, col]) < 3 * abs(nerr[2, 0]) + 3.0
    assert abs(n[3] - truth["flux_Ha"][row, col]) < 3 * abs(nerr[3, 0]) + 0.02 * truth["flux_Ha"][row, col]
    assert abs(n[5] - truth["flux_n2_r"][row, col]) < 4 * abs(nerr[5, 0]) + 0.05 * truth["flux_n2_r"][row, col]
    # continuum at the lines
    c = rs.c[row, col]
    assert abs(c[1] - truth["cont"][row, col]) < 0.5
    assert rs.chisq[row, col] > 0
    assert nerr[1, 1] == -nerr[1, 0]


def test_fixed_ratios_and_constraints(synth_cube):
    state, truth = _session(synth_cube, fix_ratios=True)
    state.pndofit[:3] = 1
    state.pcdofit[:3] = 1
    state.col, state.row = 5, 5
    dofit(state)
    n = state.res.n[5, 5]
    assert abs(n[4] * 3.071 - n[5]) < 1e-6 * max(abs(n[5]), 1)
    # user constraints: fix sigma to 40 km/s
    state.fitconstr = True
    state.pnfix[2] = 1
    state.gnfit[2] = 40.0
    dofit(state)
    assert abs(state.res.n[5, 5, 2] - 40.0) < 1e-9


def test_fitall_autoflag_save_load(synth_cube, tmp_path):
    state, truth = _session(synth_cube)
    state.pndofit[:3] = 1
    state.pcdofit[:3] = 1
    assert fitall(state, spaxel=True)
    flag = autoflag(state)
    rs = state.res
    good = flag == 0
    assert good[5, 5]
    assert flag[0, 0] >= C.FLAG_MANUAL                       # NaN spaxel is bad
    vel_dev = np.abs(rs.n[..., 1] - truth["vel"])[good]
    assert np.median(vel_dev) < 5.0
    with np.errstate(divide="ignore"):
        flux_dev = (np.abs(rs.n[..., 3] - truth["flux_Ha"]) / np.abs(rs.nerr[..., 3, 0]))[good]
    assert np.median(flux_dev) < 2.0

    fname = str(tmp_path / "res.fits")
    saveres(state, fname)
    assert os.path.exists(fname) and os.path.exists(str(tmp_path / "res_noflag.fits"))
    hdr = fits.getheader(fname)
    assert hdr["NAXIS3"] == 6 * (3 + 3) + 3 * (1 + 3) + 1
    assert hdr["P1"] == "narrow line flag" and hdr["P2"] == "narrow line dv"
    assert hdr["LINESET"] == 1 and hdr["RESTYPE"] == 0 and hdr["FITTYPE"] == 0
    data = fits.getdata(fname)
    assert np.isnan(data[1, 0, 0])                            # flagged spaxel -> NaN

    # load into a fresh session and compare
    state2, _ = _session(synth_cube)
    assert loadres(state2, fname)
    rs2 = state2.res
    np.testing.assert_allclose(rs2.n[good], rs.n[good])
    np.testing.assert_allclose(rs2.nerr[good][..., :2], rs.nerr[good][..., :2])
    np.testing.assert_allclose(rs2.chisq, rs.chisq)

    # session round trip
    sfile = save_session(state, str(tmp_path / "sess.kvz"))
    st3 = load_session(sfile)
    np.testing.assert_array_equal(st3.res.n, rs.n)
    assert st3.redshift == state.redshift


def test_fitadjall_improves_bad_spaxels(synth_cube):
    state, truth = _session(synth_cube)
    state.pndofit[:3] = 1
    state.pcdofit[:3] = 1
    fitall(state)
    autoflag(state)
    rs = state.res
    # spoil a good spaxel and flag it bad
    rs.n[5, 6, 1:] = C.NOT_FIT
    rs.n[5, 6, 0] = C.FLAG_LOW_SN
    rs.b[5, 6, 0] = C.FLAG_LOW_SN
    rs.nerr[5, 6, 1:, :] = C.NOT_FIT
    improved = fitadjall(state)
    assert improved >= 1
    assert rs.n[5, 6, 0] == 0
    assert abs(rs.n[5, 6, 1] - truth["vel"][5, 6]) < 10


def test_mask_mode_and_extraction(synth_cube):
    state, truth = _session(synth_cube)
    state.pndofit[:3] = 1
    state.pcdofit[:3] = 1
    state.maskmode = 1
    state.maskradius = 2
    select_spaxels(state, 6, 5, 1)
    assert state.spaxselect[0].sum() > 5
    state.linefit_mode = C.MODE_MASK
    state.specmode = C.SPEC_SUM
    medianspec(state)
    assert state.medspec.max() > state.datacube[:, 5, 6].max()
    dofit(state)
    rs = state.res
    assert rs.mode == C.MODE_MASK
    assert rs.n[0, 3] > truth["flux_Ha"][5, 6]
    newmask(state)
    assert state.Nmask == 2 and state.imask == 2
    state.get_results().ensure_nmask(2)
    assert state.res.n.shape[0] == 2
    # weighted average and optimal extraction run
    state.imask = 1
    state.specmode = C.SPEC_WAVG
    medianspec(state)
    assert np.all(np.isfinite(state.medspec))
    state.wavrange1[:] = [540, 570]
    state.imgmode = C.IMG_SUM1
    medsum_image_update(state)
    assert state.img1[5, 6] > state.img1[0, 5]
    state.specmode = C.SPEC_OPTIMAL
    medianspec(state)
    assert np.all(np.isfinite(state.medspec))


def test_moments(synth_cube):
    state, truth = _session(synth_cube, fittype="moments")
    state.pndofit[:3] = 1
    state.pcdofit[:3] = 1
    state.col, state.row = 6, 5
    dofit(state)
    m = state.res.m[5, 6]
    assert abs(m[0] - truth["flux_Ha"][5, 6]) / truth["flux_Ha"][5, 6] < 0.15
    assert abs(m[1] - truth["vel"][5, 6]) < 15
    assert state.res.merr[5, 6, 0, 0] == C.NO_ERRORS


def test_montecarlo2_errors(synth_cube):
    state, truth = _session(synth_cube, do_mc_errors=3, nmontecarlo=8)
    assert state.mc2cubes.shape[0] == 8 and state.Nmontecarlo == 8
    state.pndofit[:3] = 1
    state.pcdofit[:3] = 1
    state.col, state.row = 6, 5
    dofit(state)
    rs = state.res
    assert (C.MODE_SPAXEL, C.ERR_MC2) in state.results
    err = rs.nerr[5, 6, 3]
    assert err[0] > 0 and err[1] < 0                      # +1sig / -1sig offsets
    assert err[4] != C.NOT_FIT                            # median offset stored


def test_getredshift(synth_cube):
    state, truth = _session(synth_cube, redshift=0.0205)
    z = getredshift(state)
    assert abs(z - truth["redshift"]) < 5e-4


def test_batchmode(synth_cube, tmp_path):
    fname, truth = synth_cube
    state = start_session(fname, redshift=truth["redshift"], lineset=1, outdir=str(tmp_path))
    resfile, sessfile = batchmode(state, redshift=truth["redshift"], fit_all_lines=True)
    assert os.path.exists(resfile) and os.path.exists(sessfile)
    assert resfile.endswith("_res_gau_spax.fits")
