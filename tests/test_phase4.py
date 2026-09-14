"""Phase 4: IDL session import (with a fake readsav structure), Monte Carlo restoration
of light sessions, the results comparison tool and the demo cube."""
import os

import numpy as np
from astropy.io import fits

from kubeviz import constants as C
from kubeviz.cli import build_parser
from kubeviz.compare import compare_results, format_report
from kubeviz.fitting.fitall import fitall
from kubeviz.fitting.flags import autoflag
from kubeviz.fitting.linefit import dofit
from kubeviz.io import idlsession
from kubeviz.io.results import saveres
from kubeviz.io.session import load_session, save_session
from kubeviz.session import make_demo_cube, restore_montecarlo, start_session


def _fit_all(state):
    state.pndofit[:3] = 1
    state.pcdofit[:3] = 1
    fitall(state)
    autoflag(state)


def test_compare_results_tool(synth_cube, tmp_path):
    fname, truth = synth_cube
    st = start_session(fname, redshift=truth["redshift"], lineset=1)
    _fit_all(st)
    fa = str(tmp_path / "a.fits")
    saveres(st, fa)
    # perturb a copy: shift the velocity plane by a fraction of the error
    data, hdr = fits.getdata(fa, header=True)
    data2 = data.copy()
    data2[1] += 0.3 * np.abs(data[2 * 6 + 4 + 1])            # plane 1 = narrow dv, its +1sig plane follows n,b,c blocks
    fb = str(tmp_path / "b.fits")
    fits.writeto(fb, data2, hdr)
    rows, summary = compare_results(fa, fb)
    assert not summary["problems"]
    assert rows[0]["name"] == "narrow line flag" and rows[0]["max_abs"] == 0
    dv = rows[1]
    assert dv["name"] == "narrow line dv" and dv["n_both"] > 50
    assert 0.2 < dv["med_pull"] < 0.4 and dv["frac_1sig"] == 1.0
    report = format_report(rows, summary)
    assert "narrow line dv" in report and "med|d|/err" in report


def test_light_session_restores_montecarlo(synth_cube, tmp_path):
    fname, truth = synth_cube
    st = start_session(fname, redshift=truth["redshift"], lineset=1, do_mc_errors=3, nmontecarlo=4)
    st.pndofit[:3] = 1
    st.pcdofit[:3] = 1
    st.col, st.row = 6, 5
    dofit(st)
    err_before = st.res.nerr[5, 6, 3, 0]
    sfile = save_session(st, str(tmp_path / "light.kvz"), light=True)
    st2 = load_session(sfile)
    assert st2.mc2cubes is None and st2.domontecarlo == C.ERR_MC2
    restore_montecarlo(st2)
    assert st2.mc2cubes is not None and st2.mc2cubes.shape[0] == 4
    np.testing.assert_allclose(st2.mc2cubes, st.mc2cubes)          # deterministic seeds
    assert st2.res.nerr[5, 6, 3, 0] == err_before                 # results survived
    # bootstrap without the file falls back to noise-cube errors
    st2.domontecarlo = C.ERR_BOOTSTRAP
    st2.bootstrapcubes = None
    st2.bootstrap_file = str(tmp_path / "missing.fits")
    restore_montecarlo(st2)
    assert st2.domontecarlo == C.ERR_NOISE


def test_idl_session_import(synth_cube, monkeypatch):
    """Build a fake readsav() structure with IDL axis order and check the mapping."""
    fname, truth = synth_cube
    st = start_session(fname, redshift=truth["redshift"], lineset=1)
    _fit_all(st)
    rs = st.res
    nl = st.Nlines
    nw, nrow, ncol = st.shape

    def idl(arr):
        return np.ascontiguousarray(arr)

    fields = {
        "FILENAME": np.bytes_(st.filename), "INDIR": np.bytes_(st.indir), "OUTDIR": np.bytes_(""),
        "FLUXFAC": np.float32(1.0), "NOISEISVAR": np.uint8(1), "LAMBDA0": st.lambda0, "DLAMBDA": st.dlambda, "PIX0": st.pix0,
        "DOMONTECARLO": np.int32(0), "INSTR": np.bytes_("muse"), "BAND": np.bytes_("WFM"), "VACUUM": np.uint8(0),
        "COL": np.int32(6), "ROW": np.int32(5), "WPIX": np.int32(400), "NCOL": np.int32(ncol), "NROW": np.int32(nrow),
        "NWPIX": np.int32(nw), "REDSHIFT": np.float64(truth["redshift"]), "SELECTED_LINESET": np.int32(1),
        "LINEFIT_MODE": np.uint8(0), "LINEFIT_TYPE": np.uint8(0), "NMASK": np.int32(1), "IMASK": np.int32(1),
        "SMOOTH": np.int32(1), "SPECSMOOTH": np.int32(1), "INSTRRES_MODE": np.int32(1),
        "INSTRRES_EXTPOLY": idl(st.instrres_extpoly), "INSTRRES_VARPOLY": np.zeros(7),
        "INDATACUBE": idl(st.indatacube), "INNOISECUBE": idl(st.innoisecube),
        "DATACUBE": idl(st.datacube), "NOISECUBE": idl(st.noisecube),
        "BADPIXELMASK": idl(st.badpixelmask.astype(np.uint8)), "BADPIXELIMG": idl(st.badpixelimg.astype(np.int16)),
        "WAVE": idl(st.wave), "SPAXSELECT": idl(st.spaxselect), "GNLIMS": idl(st.gnlims.T), "GBLIMS": idl(st.gblims.T),
        "GNFIT": idl(st.gnfit), "GBFIT": idl(st.gbfit), "PNDOFIT": idl(st.pndofit),
        "INDATAHEAD": np.array([np.bytes_(c.image) for c in st.indatahead.cards]),
        "NOI_SP_NRESCUBE": idl(np.moveaxis(rs.n, -1, 0)), "NOI_SP_BRESCUBE": idl(np.moveaxis(rs.b, -1, 0)),
        "NOI_SP_CRESCUBE": idl(np.moveaxis(rs.c, -1, 0)), "NOI_SP_MRESCUBE": idl(np.moveaxis(rs.m, -1, 0)),
        "NOI_SP_NERRRESCUBE": idl(np.transpose(rs.nerr, (3, 2, 0, 1))), "NOI_SP_BERRRESCUBE": idl(np.transpose(rs.berr, (3, 2, 0, 1))),
        "NOI_SP_CERRRESCUBE": idl(np.transpose(rs.cerr, (3, 2, 0, 1))), "NOI_SP_MERRRESCUBE": idl(np.transpose(rs.merr, (3, 2, 0, 1))),
        "SP_CHISQ": idl(rs.chisq), "MAXNMASK": np.int32(50000),
        "NOI_NRESCUBE": idl(np.zeros((3 + nl, 50000))), "NOI_BRESCUBE": idl(np.zeros((3 + nl, 50000))),
        "NOI_CRESCUBE": idl(np.zeros((1 + nl, 50000))), "NOI_MRESCUBE": idl(np.zeros((6 * nl, 50000))),
        "NOI_NERRRESCUBE": idl(np.full((9, 3 + nl, 50000), -999.0)), "NOI_BERRRESCUBE": idl(np.full((9, 3 + nl, 50000), -999.0)),
        "NOI_CERRRESCUBE": idl(np.full((9, 1 + nl, 50000), -999.0)), "NOI_MERRRESCUBE": idl(np.full((9, 6 * nl, 50000), -999.0)),
        "CHISQ": idl(np.zeros(50000)),
    }
    dtype = np.dtype([(k, object) for k in fields])
    rec = np.zeros(1, dtype=dtype)
    for k, v in fields.items():
        rec[k][0] = v
    monkeypatch.setattr(idlsession, "readsav", lambda *a, **k: {"state": rec}, raising=False)
    import scipy.io
    monkeypatch.setattr(scipy.io, "readsav", lambda *a, **k: {"state": rec})

    st2 = idlsession.load_idl_session("fake.sav")
    assert st2.shape == st.shape and st2.instr == "muse" and st2.noiseisvar
    assert st2.Nlines == 3 and st2.linenames == st.linenames
    np.testing.assert_allclose(st2.wave, st.wave)
    np.testing.assert_array_equal(st2.badpixelimg, st.badpixelimg)
    assert st2.indatahead is not None and st2.indatahead["CRVAL3"] == st.indatahead["CRVAL3"]
    rs2 = st2.get_results(C.MODE_SPAXEL, C.ERR_NOISE)
    np.testing.assert_allclose(rs2.n, rs.n)
    np.testing.assert_allclose(rs2.nerr, rs.nerr)
    np.testing.assert_allclose(rs2.chisq, rs.chisq)
    assert (C.MODE_MASK, C.ERR_NOISE) in st2.results and st2.results[(C.MODE_MASK, C.ERR_NOISE)].n.shape[0] == 1


def test_demo_cube_and_cli_flags(tmp_path):
    fname = make_demo_cube(str(tmp_path), nx=8, ny=6)
    assert os.path.exists(fname)
    args = build_parser().parse_args(["--demo", "--batch"])
    assert args.demo and args.batch and args.lineset == 0
    st = start_session(fname, redshift=0.02, lineset=1)
    assert st.Nlines == 3
