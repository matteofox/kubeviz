"""Best-effort loading of IDL kubeviz sessions (``.sav`` files written by ``kubeviz_savesession``).

The IDL ``state`` structure is read with :func:`scipy.io.readsav`, which dereferences
the heap pointers, and its fields are mapped onto :class:`~kubeviz.state.State`.
Array axes are reversed by ``readsav`` (IDL ``[col,row,wave]`` becomes numpy
``(wave,row,col)``), which is exactly the layout used by the Python port; result
containers are transposed to the ``(row, col, par[, perc])`` layout.

Monte Carlo realisation cubes are not restored (they are regenerated on demand, see
:func:`kubeviz.session.restore_montecarlo`). Colour tables and widget ids are ignored.
This code cannot be exercised without IDL-produced files, so every field is read
defensively; anything that fails to convert keeps its default value and is reported.
"""
from __future__ import annotations

import numpy as np
from astropy.io import fits

from .. import utils
from ..constants import ERR_BOOTSTRAP, ERR_MC1, ERR_MC2, ERR_MC3, ERR_NOISE, MODE_MASK, MODE_SPAXEL
from ..results import ResultSet
from ..state import State

_SCALARS = {
    # IDL field -> (State attribute, converter)
    "filename": ("filename", str), "indir": ("indir", str), "outdir": ("outdir", str),
    "fluxfac": ("fluxfac", float), "noiseisvar": ("noiseisvar", bool),
    "lambda0": ("lambda0", float), "dlambda": ("dlambda", float), "pix0": ("pix0", float),
    "domontecarlo": ("domontecarlo", int), "nmontecarlo": ("Nmontecarlo", int),
    "plotmontecarlodistrib": ("plotMonteCarlodistrib", bool), "savemontecarlodistrib": ("saveMonteCarlodistrib", bool),
    "usemontecarlonoise": ("useMonteCarlonoise", bool), "scalenoiseerrors": ("scaleNoiseerrors", bool),
    "zmin_ima": ("zmin_ima", float), "zmax_ima": ("zmax_ima", float),
    "zmin_spec": ("zmin_spec", float), "zmax_spec": ("zmax_spec", float),
    "instr": ("instr", str), "band": ("band", str), "ifu": ("ifu", int), "vacuum": ("vacuum", bool),
    "pixscale": ("pixscale", int), "cubesel": ("cubesel", int), "col": ("col", int), "row": ("row", int),
    "wpix": ("wpix", int), "startcol": ("Startcol", int), "startrow": ("Startrow", int),
    "startwpix": ("Startwpix", int), "ncol": ("Ncol", int), "nrow": ("Nrow", int), "nwpix": ("Nwpix", int),
    "zoommap": ("zoommap", bool), "marker": ("marker", int), "scale": ("scale", int),
    "zoomrange": ("zoomrange", int), "zcuts": ("zcuts", int), "cursormode": ("cursormode", int),
    "specmode": ("specmode", int), "imgmode": ("imgmode", int), "flagmode": ("flagmode", bool),
    "wavsel": ("wavsel", int), "invert": ("invert", int), "ctab": ("ctab", int),
    "smooth": ("smooth", int), "specsmooth": ("specsmooth", int), "transpose": ("transpose", bool),
    "nmask": ("Nmask", int), "imask": ("imask", int), "maskradius": ("maskradius", int), "maskmode": ("maskmode", int),
    "continuumfit_mode": ("continuumfit_mode", int), "continuumfit_minoff": ("continuumfit_minoff", float),
    "continuumfit_maxoff": ("continuumfit_maxoff", float), "continuumfit_minperc": ("continuumfit_minperc", float),
    "continuumfit_maxperc": ("continuumfit_maxperc", float), "continuumfit_order": ("continuumfit_order", int),
    "maxwoffb": ("maxwoffb", float), "maxwoffr": ("maxwoffr", float),
    "secondcomp_mode": ("secondcomp_mode", int), "secondcomp_smart": ("secondcomp_smart", bool),
    "instrres_mode": ("instrres_mode", int), "instrres_tplsig": ("instrres_tplsig", float),
    "nloopadj": ("Nloopadj", int), "fitconstr": ("fitconstr", bool), "fitfixratios": ("fitfixratios", bool),
    "redshift": ("redshift", float), "selected_lineset": ("selected_lineset", int),
    "linefit_mode": ("linefit_mode", int), "linefit_type": ("linefit_type", int),
    "mask_sn_thresh": ("mask_sn_thresh", float), "mask_maxvelerr": ("mask_maxvelerr", float),
    "mom_thresh": ("mom_thresh", float), "nbootstrap": ("Nbootstrap", int),
    "nmc1": ("Nmc1", int), "nmc2": ("Nmc2", int), "nmc3": ("Nmc3", int), "par_imagebutton": ("par_imagebutton", str),
}
_ARRAYS_SAME = ["indatacube", "innoisecube", "datacube", "noisecube", "badpixelimg", "wave", "medspec",
                "nmedspec", "montecarlo_percs", "spaxselect", "img1", "img2", "nimg1", "nimg2",
                "wavrange1", "wavrange2", "fitallrange", "gnfit", "gbfit", "pnfix", "pbfix", "pndofit",
                "pbdofit", "pcdofit", "nshow", "bshow", "cshow", "instrres_extpoly", "instrres_varpoly"]
_METHOD_PREFIX = {"noi": ERR_NOISE, "boot": ERR_BOOTSTRAP, "mc1": ERR_MC1, "mc2": ERR_MC2, "mc3": ERR_MC3}


def _decode(v):
    if isinstance(v, bytes):
        return v.decode("latin-1").strip()
    if isinstance(v, np.ndarray) and v.dtype.kind in "SO" and v.ndim == 0:
        return _decode(v.item())
    return v


class _Rec:
    """Case-insensitive access to the fields of a readsav structure."""

    def __init__(self, rec):
        self.rec = rec
        self.names = {n.lower(): n for n in rec.dtype.names}

    def get(self, name, default=None):
        key = self.names.get(name.lower())
        if key is None:
            return default
        v = self.rec[key]
        if isinstance(v, np.ndarray) and v.shape == (1,):
            v = v[0]
        return v

    def __contains__(self, name):
        return name.lower() in self.names


def _array(v):
    """Dereferenced pointer content as a float array, or None."""
    if v is None:
        return None
    if isinstance(v, np.ndarray) and v.dtype == object:
        if v.size == 1:
            return _array(v.ravel()[0])
        return None
    if hasattr(v, "dtype"):
        return np.asarray(v)
    return None


def _header(v):
    arr = _array(v) if not isinstance(v, np.ndarray) else v
    if arr is None or arr.dtype.kind not in "SU":
        return None
    cards = [(_decode(c) if isinstance(c, bytes) else str(c)) for c in np.asarray(arr).ravel()]
    text = "".join(c.ljust(80)[:80] for c in cards if c.strip())
    try:
        return fits.Header.fromstring(text)
    except Exception:
        return None


def load_idl_session(fname: str) -> State:
    from scipy.io import readsav

    utils.info(f"Loading IDL session: {fname}")
    data = readsav(fname, python_dict=True, verbose=False)
    key = next((k for k in data if k.lower() == "state"), None)
    if key is None:
        raise ValueError(f"{fname} does not contain a kubeviz 'state' structure")
    rec = _Rec(data[key])
    st = State()
    problems = []

    for idl, (attr, conv) in _SCALARS.items():
        if idl not in rec:
            continue
        try:
            v = _decode(rec.get(idl))
            if isinstance(v, np.ndarray):
                v = v.item() if v.size == 1 else v
            setattr(st, attr, conv(v))
        except Exception as exc:
            problems.append(f"{idl}: {exc}")

    for idl in _ARRAYS_SAME:
        if idl not in rec:
            continue
        try:
            arr = _array(rec.get(idl))
            if arr is not None:
                setattr(st, idl if idl not in ("badpixelimg",) else idl, arr)
        except Exception as exc:
            problems.append(f"{idl}: {exc}")

    # 2D limit arrays are stored as [25,2] in IDL -> (2,25) here
    for idl in ("gnlims", "gblims"):
        arr = _array(rec.get(idl))
        if arr is not None and arr.ndim == 2:
            setattr(st, idl, np.ascontiguousarray(arr.T if arr.shape[0] == 2 else arr))
    bm = _array(rec.get("badpixelmask"))
    if bm is not None:
        st.badpixelmask = bm.astype(bool)
    if st.badpixelimg is not None:
        st.badpixelimg = np.asarray(st.badpixelimg).astype(bool)
    if st.spaxselect is not None and st.spaxselect.ndim == 2:
        st.spaxselect = st.spaxselect[None]
    for idl, attr in (("inprihead", "inprihead"), ("indatahead", "indatahead"), ("innoisehead", "innoisehead")):
        st.__dict__[attr] = _header(rec.get(idl))

    # dimensions from the data if the scalars are missing
    if st.datacube is not None:
        st.Nwpix, st.Nrow, st.Ncol = st.datacube.shape
    if st.indatacube is None:
        st.indatacube = st.datacube
    if st.innoisecube is None:
        st.innoisecube = st.noisecube
    st.noise = st.noisecube
    if st.wave is None and st.Nwpix:
        st.wave = np.arange(st.Nwpix, dtype=float)
    st.cwdir = ""
    st.Nmask = max(int(st.Nmask), 1)

    # lines and linesets are recomputed from redshift/lineset (they depend only on the
    # catalogue and the wavelength range) without touching the loaded results
    from ..fitting.linefit import linefit_init
    linefit_init(st, reset_results=False)

    # result containers
    nl = st.Nlines
    for prefix, err in _METHOD_PREFIX.items():
        for mode, tag in ((MODE_SPAXEL, "sp_"), (MODE_MASK, "")):
            n = _array(rec.get(f"{prefix}_{tag}nrescube"))
            if n is None or n.size == 0:
                continue
            try:
                rs = ResultSet(mode, nl, ncol=st.Ncol, nrow=st.Nrow, nmask=st.Nmask)
                if mode == MODE_SPAXEL:
                    conv = lambda a: np.moveaxis(a, 0, -1)                 # (Npar,Nrow,Ncol) -> (Nrow,Ncol,Npar)
                    conve = lambda a: np.transpose(a, (2, 3, 1, 0))        # (Npercs,Npar,Nrow,Ncol)
                else:
                    conv = lambda a: np.ascontiguousarray(a.T)[: st.Nmask]  # (Npar,maxNmask) -> (Nmask,Npar)
                    conve = lambda a: np.transpose(a, (2, 1, 0))[: st.Nmask]
                for name, arr_name, cv in (("n", "nrescube", conv), ("b", "brescube", conv), ("c", "crescube", conv),
                                           ("m", "mrescube", conv), ("nerr", "nerrrescube", conve),
                                           ("berr", "berrrescube", conve), ("cerr", "cerrrescube", conve),
                                           ("merr", "merrrescube", conve)):
                    a = _array(rec.get(f"{prefix}_{tag}{arr_name}"))
                    if a is None:
                        continue
                    a = cv(np.asarray(a, dtype=float))
                    target = getattr(rs, name)
                    if a.shape == target.shape:
                        setattr(rs, name, a.copy())
                    else:
                        problems.append(f"{prefix}_{tag}{arr_name}: shape {a.shape} != {target.shape}")
                chisq = _array(rec.get("sp_chisq" if mode == MODE_SPAXEL else "chisq"))
                if chisq is not None:
                    chisq = np.asarray(chisq, dtype=float)
                    if mode == MODE_MASK:
                        chisq = chisq[: st.Nmask]
                    if chisq.shape == rs.chisq.shape:
                        rs.chisq = chisq.copy()
                st.results[(mode, err)] = rs
            except Exception as exc:
                problems.append(f"results {prefix}/{tag}: {exc}")

    if problems:
        utils.warn("Some fields of the IDL session could not be converted:")
        for p in problems:
            utils.warn("  " + p)
    utils.info("IDL session loaded (Monte Carlo cubes are regenerated on demand).")
    return st
