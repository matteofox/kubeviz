"""Emission-line database (port of the ``kubeviz_linesdb`` common block).

The order of the lines is significant: it is the order used by ``state.lines`` and
therefore by the planes of the result cubes. Do not reorder.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class Line:
    name: str          # short name used in headers and fixed ratios (IDL linenames)
    fancy: str         # label used in the GUI (IDL linefancynames)
    rest_air: float    # rest-frame AIR wavelength [Angstrom]
    lineset: int       # lineset id (1..9)


# ALL rest wavelengths are AIR wavelengths (table 13.2 Osterbrock for [OI]).
LINES: tuple[Line, ...] = (
    Line("Ha",    "Ha",           6562.819, 1),
    Line("n2_b",  "[NII] blue",   6548.050, 1),
    Line("n2_r",  "[NII] red",    6583.450, 1),
    Line("o3_b",  "[OIII] blue",  4958.911, 2),
    Line("o3_r",  "[OIII] red",   5006.843, 2),
    Line("o2_b",  "[OII] blue",   3726.030, 3),
    Line("o2_r",  "[OII] red",    3728.820, 3),
    Line("Hb",    "Hb",           4861.333, 4),
    Line("s2_b",  "[SII] blue",   6716.440, 5),
    Line("s2_r",  "[SII] red",    6730.810, 5),
    Line("o1_b",  "[OI] blue",    6300.304, 6),
    Line("o1_r",  "[OI] red",     6363.776, 6),
    Line("s3_b",  "[SIII] blue",  9068.600, 7),
    Line("s3_r",  "[SIII] red",   9530.600, 7),
    Line("He1_b", "HeI blue",     5875.621, 8),
    Line("He1_r", "HeI red",      6678.152, 8),
    Line("Lya",   "Lya",          1215.670, 9),
)

LINESETS = np.array([l.lineset for l in LINES], dtype=int)
LINES_REST = np.array([l.rest_air for l in LINES], dtype=float)
LINENAMES = [l.name for l in LINES]
LINEFANCYNAMES = [l.fancy for l in LINES]
LINESET_MAX = int(LINESETS.max())

LINESET_DESCRIPTIONS = {
    0: "All lines",
    1: "Ha + [NII]",
    2: "[OIII]",
    3: "[OII]",
    4: "Hbeta",
    5: "[SII]",
    6: "[OI]",
    7: "[SIII]",
    8: "HeI",
    9: "Lya",
}

# Priority used to define the "main line" of a set of lines (kubeviz_getmainline)
MAINLINE_PRIORITY = ("Ha", "o3_r", "o2_r", "Hb", "s2_r", "o1_b", "He1_b", "s3_b", "Lya")

# Fixed line ratios (Storey & Zeippen 2000): faint = bright / ratio
FIXED_RATIOS = {
    ("n2_b", "n2_r"): 3.071,
    ("o3_b", "o3_r"): 3.013,
    ("o1_r", "o1_b"): 2.997,
}


def rest_wavelength(name: str) -> float:
    """Rest-frame air wavelength of a line by short name."""
    for l in LINES:
        if l.name == name:
            return l.rest_air
    raise KeyError(name)


def airtovac(wave_air):
    """Air -> vacuum wavelength [Angstrom] (IDL astrolib ``airtovac``, Ciddor 1996).

    Wavelengths below 2000 A are returned unchanged, as in astrolib.
    """
    wave_air = np.asarray(wave_air, dtype=float)
    sigma2 = (1.0e4 / wave_air) ** 2
    fact = 1.0 + 5.792105e-2 / (238.0185 - sigma2) + 1.67917e-3 / (57.362 - sigma2)
    out = np.where(wave_air >= 2000.0, wave_air * fact, wave_air)
    return float(out) if out.ndim == 0 else out


def vactoair(wave_vac):
    """Vacuum -> air wavelength [Angstrom] (IDL astrolib ``vactoair``)."""
    wave_vac = np.asarray(wave_vac, dtype=float)
    sigma2 = (1.0e4 / wave_vac) ** 2
    fact = 1.0 + 5.792105e-2 / (238.0185 - sigma2) + 1.67917e-3 / (57.362 - sigma2)
    out = np.where(wave_vac >= 2000.0, wave_vac / fact, wave_vac)
    return float(out) if out.ndim == 0 else out


def mainline_rest(linenames) -> float:
    """Rest air wavelength of the main line among ``linenames`` (kubeviz_getmainline).

    Falls back to Halpha when none of the priority lines is present.
    """
    names = list(linenames) if linenames is not None else []
    for cand in MAINLINE_PRIORITY:
        if cand in names:
            return rest_wavelength(cand)
    return rest_wavelength("Ha")
