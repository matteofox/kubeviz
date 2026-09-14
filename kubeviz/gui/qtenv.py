"""Make PyQt6 find its own Qt plugins.

Some Python distributions (anaconda ships a ``bin/qt.conf`` for its Qt5) redirect the
Qt plugin search path to an incompatible Qt installation, and the application then
aborts with "Could not find the Qt platform plugin". Calling
:func:`prepare_qt_environment` before the first Qt import points Qt6 at the plugins
bundled with the PyQt6 wheel.
"""
from __future__ import annotations

import importlib.util
import os


def prepare_qt_environment() -> None:
    if os.environ.get("KUBEVIZ_NO_QT_ENV_FIX"):
        return
    spec = importlib.util.find_spec("PyQt6")
    if spec is None or not spec.submodule_search_locations:
        return
    for base in spec.submodule_search_locations:
        plugins = os.path.join(base, "Qt6", "plugins")
        if os.path.isdir(os.path.join(plugins, "platforms")):
            os.environ["QT_PLUGIN_PATH"] = plugins
            os.environ["QT_QPA_PLATFORM_PLUGIN_PATH"] = os.path.join(plugins, "platforms")
            return


APP_FONT_POINTS = 12


def apply_app_font(points: int = APP_FONT_POINTS) -> None:
    """One font size for the whole GUI (widgets, plots and dialogs)."""
    from PyQt6.QtWidgets import QApplication
    app = QApplication.instance()
    if app is None:
        return
    f = app.font()
    if f.pointSize() != points:
        f.setPointSize(points)
        app.setFont(f)
