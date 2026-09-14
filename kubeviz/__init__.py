"""kubeviz: interactive visualisation and emission-line fitting of IFU datacubes.

Python port of the IDL ``kubeviz`` (Fossati, Wilman & Gerssen). The headless engine
lives in :mod:`kubeviz.core`, :mod:`kubeviz.fitting` and :mod:`kubeviz.io`; the Qt
GUI in :mod:`kubeviz.gui`; the command line entry point in :mod:`kubeviz.cli`.
"""

__version__ = "3.0.0a0"
IDL_VERSION_PORTED = "K2.2"
