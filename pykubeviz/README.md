# PyKubeviz installation & Dependency Guide

This document lists the dependencies and necessary files required to run and install **PyKubeviz** (the Python 3 port of the IDL `kubeviz` software).

---

## 1. System & Python Dependencies

The software is written in Python 3 and requires the following third-party libraries:

*   **numpy**: For numerical array operations and data cube manipulations.
*   **scipy**: For numerical algorithms and interpolation.
*   **astropy**: For FITS file I/O operations and handling astronomical coordinates.
*   **lmfit**: For curve fitting and parameter constraints (tied ratios, lower/upper bounds).
*   **PyQt6**: For the GUI (Graphical User Interface) framework.
*   **pyqtgraph**: For high-performance, interactive 1D (spectrum) and 2D (spaxel/map) visualization.

These dependencies are captured in the root [requirements.txt](file:///Users/matteo/Software/kubeviz/requirements.txt) file:
```text
numpy
scipy
astropy
lmfit
PyQt6
pyqtgraph
```

---

## 2. Necessary Files & Directory Structure

To run the application, the `pykubeviz` directory must contain the following core files and submodules:

*   **[pykubeviz/](file:///Users/matteo/Software/kubeviz/pykubeviz)** - Parent package directory
    *   **[__init__.py](file:///Users/matteo/Software/kubeviz/pykubeviz/__init__.py)** - Package initializer
    *   **[core/](file:///Users/matteo/Software/kubeviz/pykubeviz/core)** - Core state and database modules
        *   **[__init__.py](file:///Users/matteo/Software/kubeviz/pykubeviz/core/__init__.py)**
        *   **[state.py](file:///Users/matteo/Software/kubeviz/pykubeviz/core/state.py)** - Manages the global session state (`KubevizState`)
        *   **[linesdb.py](file:///Users/matteo/Software/kubeviz/pykubeviz/core/linesdb.py)** - Emission line database and metadata
    *   **[io/](file:///Users/matteo/Software/kubeviz/pykubeviz/io)** - Input/Output handling
        *   **[__init__.py](file:///Users/matteo/Software/kubeviz/pykubeviz/io/__init__.py)**
        *   **[data_io.py](file:///Users/matteo/Software/kubeviz/pykubeviz/io/data_io.py)** - Fits file reading & metadata initialization
        *   **[export.py](file:///Users/matteo/Software/kubeviz/pykubeviz/io/export.py)** - Exporter for fits fit results
    *   **[analysis/](file:///Users/matteo/Software/kubeviz/pykubeviz/analysis)** - Fitting and mathematical analysis
        *   **[__init__.py](file:///Users/matteo/Software/kubeviz/pykubeviz/analysis/__init__.py)**
        *   **[fitting.py](file:///Users/matteo/Software/kubeviz/pykubeviz/analysis/fitting.py)** - Main fitting routines using `lmfit` and moment calculation
        *   **[errors.py](file:///Users/matteo/Software/kubeviz/pykubeviz/analysis/errors.py)** - Error estimation and noise simulation
        *   **[calibration.py](file:///Users/matteo/Software/kubeviz/pykubeviz/analysis/calibration.py)** - Wavelength and velocity calibrations
        *   **[utils.py](file:///Users/matteo/Software/kubeviz/pykubeviz/analysis/utils.py)** - Math helpers: clipping, weighted median, data cleaning
    *   **[ui/](file:///Users/matteo/Software/kubeviz/pykubeviz/ui)** - Graphic User Interface component
        *   **[__init__.py](file:///Users/matteo/Software/kubeviz/pykubeviz/ui/__init__.py)**
        *   **[main_window.py](file:///Users/matteo/Software/kubeviz/pykubeviz/ui/main_window.py)** - Main application window and entry point
        *   **[viewers.py](file:///Users/matteo/Software/kubeviz/pykubeviz/ui/viewers.py)** - Plot components for 1D and 2D viewer panes

---

## 3. How to Install & Run

1.  **Navigate to the root directory** of the repository:
    ```bash
    cd /Users/matteo/Software/kubeviz
    ```
2.  **Create and activate a virtual environment** (optional but highly recommended):
    ```bash
    python3 -m venv venv
    source venv/bin/activate
    ```
3.  **Install the dependencies**:
    ```bash
    pip install -r requirements.txt
    ```
4.  **Run the application**:
    ```bash
    python3 -m pykubeviz.ui.main_window
    ```
