"""GUI entry point."""
from __future__ import annotations

import sys


def run_gui(state, args: dict | None = None) -> int:
    """Start the Qt application on ``state`` (a populated State) or, when ``state`` is
    None, open a file dialog to choose a cube first."""
    from .qtenv import prepare_qt_environment
    prepare_qt_environment()
    from PyQt6.QtWidgets import QApplication

    from .controller import KubevizGUI

    app = QApplication.instance() or QApplication(sys.argv[:1])
    app.setApplicationName("kubeviz")
    if state is None:
        from .controller import open_cube_interactively
        state = open_cube_interactively(args or {})
        if state is None:
            return 1
    gui = KubevizGUI(state, args or {})
    gui.show_all()
    return app.exec()
