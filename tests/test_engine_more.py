import os

import numpy as np
from astropy.io import fits

from kubeviz import constants as C
from kubeviz.core.masks import optimal_mask, save_mask, select_spaxels
from kubeviz.fitting.fitall import fitadj, fitall
from kubeviz.fitting.flags import autoflag
from kubeviz.fitting.linefit import dofit, linefit_reset, linefit_resetuser
from kubeviz.io.results import loadres, saveres
from kubeviz.io.spectra import savecube, saveimage, savespec
from kubeviz.session import start_session
from kubeviz.synth import make_synthetic_cube


def _session(synth_cube, **kw):
    fname, truth = synth_cube
    opts = dict(redshift=truth["redshift"], lineset=1)
    opts.update(kw)
    return start_session(fname, **opts), truth


def test_continuum_mode_mpfit(synth_cube):
    state, truth = _session(synth_cube, fitpars=[200, 500, 40, 60, 0, 1])
    assert state.continuumfit_mode == 1
    state.pndofit[:3] = 1
    state.pcdofit[:3] = 1
    state.col, state.row = 6, 5
    ctx = dofit(state)[0]
    assert not ctx.parinfo[0].fixed
    rs = state.res
    assert abs(rs.c[5, 6, 1] - truth["cont"][5, 6]) < 1.0
    assert abs(rs.n[5, 6, 3] - truth["flux_Ha"][5, 6]) < 4 * abs(rs.nerr[5, 6, 3, 0]) + 0.03 * truth["flux_Ha"][5, 6]


def test_broad_component_and_swap(tmp_path):
    fname = str(tmp_path / "broad.fits")
    _, truth = make_synthetic_cube(fname, nx=6, ny=6, redshift=0.02, broad=True, noise_level=0.5)
    state = start_session(fname, redshift=0.02, lineset=1)
    state.pndofit[:3] = 1
    state.pbdofit[:3] = 1
    state.pcdofit[:3] = 1
    state.secondcomp_mode = 2          # broad = larger sigma
    state.col, state.row = 3, 3
    ctx = dofit(state)[0]
    assert ctx.nbfitlines == 3
    rs = state.res
    n, b = rs.n[3, 3], rs.b[3, 3]
    assert b[2] > n[2]                  # broad component has the larger width after ordering
    assert abs(n[2] - truth["sigma"][3, 3]) < 15
    assert abs(b[2] - 3 * truth["sigma"][3, 3]) < 60
    assert b[3] > 0 and rs.berr[3, 3, 3, 0] > 0
    # smart mode with two coincident components should collapse to one
    state.secondcomp_smart = True
    state.secondcomp_mode = 0
    linefit_reset(state)
    ctx = dofit(state)[0]
    assert ctx.nbfitlines == 0
    assert rs.b[3, 3, 0] == C.FLAG_MANUAL or rs.b[3, 3, 3] == C.NOT_FIT


def test_moments_with_montecarlo(synth_cube):
    state, truth = _session(synth_cube, fittype="moments", do_mc_errors=3, nmontecarlo=6)
    state.pndofit[:3] = 1
    state.pcdofit[:3] = 1
    state.col, state.row = 6, 5
    dofit(state)
    rs = state.res
    # percentile offsets: -1sig below +1sig, both defined (with few realisations the
    # data value may lie outside the 16-84 range, so signs are not guaranteed)
    assert rs.merr[5, 6, 0, 1] < rs.merr[5, 6, 0, 0]
    assert np.all(rs.merr[5, 6, :5, :] != C.NOT_FIT)
    autoflag(state)
    assert rs.m[5, 6, 5] == 0


def test_bootstrap_cubes(tmp_path):
    fname = str(tmp_path / "boot_base.fits")
    _, truth = make_synthetic_cube(fname, nx=6, ny=6, redshift=0.02)
    hdul = fits.open(fname)
    data = hdul["DATA"].data
    noise = np.sqrt(hdul["STAT"].data)
    rng = np.random.default_rng(5)
    pri = fits.PrimaryHDU()
    pri.header["NEXT"] = 5
    exts = [fits.ImageHDU((data + rng.standard_normal(data.shape) * noise).astype(np.float32)) for _ in range(5)]
    bfile = str(tmp_path / "boot.fits")
    fits.HDUList([pri] + exts).writeto(bfile)
    state = start_session(fname, redshift=0.02, lineset=1, bootstrap=bfile, use_mc_noise=True)
    assert state.domontecarlo == C.ERR_BOOTSTRAP
    assert state.bootstrapcubes.shape == (5, 800, 6, 6)
    assert state.noise is state.bootstrapnoise
    state.pndofit[:3] = 1
    state.pcdofit[:3] = 1
    state.col, state.row = 3, 3
    dofit(state)
    rs = state.res
    assert (C.MODE_SPAXEL, C.ERR_BOOTSTRAP) in state.results
    assert rs.nerr[3, 3, 3, 0] > 0
    # save/plot MC distributions
    state.saveMonteCarlodistrib = True
    state.outdir = str(tmp_path) + os.sep
    dofit(state)
    pdf = str(tmp_path / "boot_base_gau_spax_PDF.fits")
    assert os.path.exists(pdf)
    with fits.open(pdf) as h:
        assert h[1].name == "SPX_4_4"
        assert h[1].data.shape == (5, 5 + 3)           # dv, sig, 3 fluxes + 3 continua
        assert h[1].header["COL1"] == "narrow line dv"
    state.plotMonteCarlodistrib = True
    state.cwdir = str(tmp_path) + os.sep
    dofit(state)
    assert os.path.isdir(str(tmp_path / "bootstrap"))


def test_mask_results_roundtrip(synth_cube, tmp_path):
    state, truth = _session(synth_cube)
    state.pndofit[:3] = 1
    state.pcdofit[:3] = 1
    state.maskmode = 2
    state.maskradius = 3
    select_spaxels(state, 6, 5, 1)
    assert state.spaxselect[0].sum() == 9
    mfile = str(tmp_path / "mask.fits")
    save_mask(state, mfile)
    state.linefit_mode = C.MODE_MASK
    state.specmode = C.SPEC_MEDIAN
    dofit(state)
    autoflag(state)
    rs = state.res
    assert rs.n[0, 0] == 0
    fname = str(tmp_path / "mres.fits")
    saveres(state, fname)
    hdr = fits.getheader(fname)
    assert hdr["NAXIS"] == 2 and hdr["NAXIS2"] == 6 * 6 + 3 * 4 and hdr["RESTYPE"] == 1
    state2, _ = _session(synth_cube, spmask=mfile)
    assert state2.linefit_mode == C.MODE_MASK and state2.Nmask == 1
    assert loadres(state2, fname)
    np.testing.assert_allclose(state2.get_results(C.MODE_MASK, C.ERR_NOISE).n, rs.n)


def test_optimal_mask_and_fitadj(synth_cube):
    state, truth = _session(synth_cube)
    state.pndofit[:3] = 1
    state.pcdofit[:3] = 1
    state.wavrange1[:] = [540, 570]
    state.imgmode = C.IMG_SUM1
    from kubeviz.core.extraction import medsum_image_update
    medsum_image_update(state)
    assert optimal_mask(state)
    assert state.spaxselect[0][5, 6] == 1 and state.spaxselect[0][0, 0] == 0
    state.linefit_mode = C.MODE_SPAXEL
    state.cursormode = 1
    fitall(state)
    autoflag(state)
    state.col, state.row = 5, 5
    assert fitadj(state)
    linefit_resetuser(state, all_pars=True)
    assert np.all(state.gnfit == C.NOT_FIT)


def test_save_spectrum_cube_image(synth_cube, tmp_path):
    state, truth = _session(synth_cube)
    state.col, state.row = 6, 5
    sp = str(tmp_path / "spec.fits")
    savespec(state, sp)
    with fits.open(sp) as h:
        assert h[0].header["SPTYPE"] == "SPAXEL" and h[1].header["EXT_TY"] == "NOISE"
        assert h[0].data.size == state.Nwpix
    cb = str(tmp_path / "cube.fits")
    savecube(state, cb)
    with fits.open(cb) as h:
        assert h[1].data.shape == (800, 10, 12)
        np.testing.assert_allclose(h[2].data[:, 5, 5], state.noisecube[:, 5, 5] ** 2, rtol=1e-5)
        assert h[1].header["CRVAL3"] == state.wave[0]
    im = str(tmp_path / "img.fits")
    state.imgmode = C.IMG_SLICE
    saveimage(state, im)
    with fits.open(im) as h:
        assert h[1].data.shape == (10, 12)
        assert "NAXIS3" not in h[1].header


def test_skylines_fit_runs_on_generic_instrument(tmp_path):
    """A cube with an unknown instrument triggers the sky-line resolution fit; with a
    smooth synthetic noise cube no lines are found and the code must degrade gracefully."""
    fname = str(tmp_path / "generic.fits")
    hdul, truth = make_synthetic_cube(None, nx=5, ny=5, redshift=0.02)
    hdul[0].header["INSTRUME"] = "FAKEIFU"
    # inject fake sky lines in the variance so the automatic search has something to find
    wave = truth["wave"]
    var = hdul["STAT"].data
    for lam in np.arange(6050, 6950, 60.0):
        var += (40.0 * np.exp(-0.5 * ((wave - lam) / 1.1) ** 2))[:, None, None]
    hdul["STAT"].data = var
    hdul.writeto(fname)
    state = start_session(fname, redshift=0.02, lineset=1)
    assert state.instr == "fakeifu"
    R = state.getinstrres(6700.0)
    assert np.isfinite(R)
    if state.instrres_mode == C.INSTRRES_VARPOLY and np.any(state.instrres_varpoly != 0):
        assert 1000 < R < 6000
