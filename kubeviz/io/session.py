"""Session save/load (ports of ``kubeviz_savesession`` / ``kubeviz_loadsession``).

Python sessions are gzip-compressed pickles of the :class:`~kubeviz.state.State`
with a format tag (extension ``.kvz``). Loading is tolerant to fields added in later
versions: missing attributes take their default value.
"""
from __future__ import annotations

import gzip
import pickle

from .. import __version__, utils
from ..state import State

FORMAT_TAG = "kubeviz-session-1"
SESSION_EXT = ".kvz"

_TRANSIENT = ("on_userpars_changed",)
_HEAVY = ("inbootstrapcubes", "bootstrapcubes", "inmc1cubes", "mc1cubes", "mc2cubes", "mc3cubes")


def save_session(state: State, fname: str, light: bool = False) -> str:
    """Write the whole state to ``fname``. ``light=True`` drops the Monte Carlo
    realisation cubes (they can be regenerated) to keep the file small."""
    if not fname.endswith(SESSION_EXT):
        fname += SESSION_EXT
    payload = dict(state.__dict__)
    for k in _TRANSIENT:
        payload[k] = None
    if light:
        for k in _HEAVY:
            payload[k] = None
    blob = {"tag": FORMAT_TAG, "version": __version__, "state": payload}
    with gzip.open(fname, "wb", compresslevel=4) as fh:
        pickle.dump(blob, fh, protocol=pickle.HIGHEST_PROTOCOL)
    utils.info("Session saved.")
    return fname


def load_session(fname: str) -> State:
    """Read a ``.kvz`` session (or raise for IDL ``.sav`` files)."""
    if fname.lower().endswith(".sav"):
        return load_idl_session(fname)
    with gzip.open(fname, "rb") as fh:
        blob = pickle.load(fh)
    if not isinstance(blob, dict) or blob.get("tag") != FORMAT_TAG:
        raise ValueError(f"{fname} is not a kubeviz session file")
    state = State()
    known = set(state.__dict__)
    for k, v in blob["state"].items():
        if k in known:
            setattr(state, k, v)
    utils.info(f"Loading session: {fname}")
    return state


def load_idl_session(fname: str) -> State:  # pragma: no cover - phase 4
    raise NotImplementedError(
        "Loading IDL .sav sessions is planned (best effort via scipy.io.readsav) but not implemented yet.")
