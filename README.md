# kubeviz

Interactive visualisation and emission-line fitting of integral-field spectroscopy (IFU)
datacubes, in Python.

kubeviz shows a FITS datacube (x, y, wavelength) in an interactive viewer and fits emission
lines spaxel by spaxel, or over spaxel masks you draw, with Gaussian profiles or moments. It
estimates uncertainties from the noise cube, bootstrap cubes or Monte Carlo realisations,
flags the fits by quality, and writes the results as FITS maps. Everything the GUI does can
also run without it, from the command line or from Python.

The original IDL implementation (V2.2) is kept in [`idl/`](idl/README.md). Results files
have the same layout in both versions and can be exchanged.

## Features

**Data input**
* Instruments recognised from the header: MUSE, KMOS (including KMOS3D data products),
  SINFONI, VIMOS, SAMI and WiFeS. Any other cube with a standard spectral WCS can be read by
  giving the data and noise extensions.
* Noise from an extension of the same file or from a separate file; air or vacuum and
  linear or logarithmic wavelength solutions.
* Load a sub-cube with `--trim '[x1:x2,y1:y2,z1:z2]'` (1-based, IRAF / ds9 convention);
  spatial and spectral smoothing kernels.

**Visualisation**
* Images of a single slice or combined over one or two wavelength ranges: sum, median,
  inverse-variance weighted average or median, median-subtracted slice and the difference
  of the two ranges.
* Data, noise, bad-pixel and S/N cubes, and maps of every fitted quantity with their
  errors and S/N. Spaxels that were never fitted are left blank.
* Continuous zoom and pan. Intensity scalings HistEq, Zscale, MinMax, percentiles and user
  linear, sqrt or log cuts. An interactive colour bar labelled with the displayed quantity
  sets the cuts by dragging its handles.
* RA / Dec of the current spaxel, 2D Gaussian fit to the displayed image, FITS header viewer.
* Spectrum panel with zoom and pan, the fitted components overplotted, and the continuum
  windows marked.

**Spectra and masks**
* Spectrum of a single spaxel, or of a mask as sum, median, weighted average,
  median-subtracted or optimal extraction (Robertson 1986).
* Any number of spaxel masks, painted with the mouse or built automatically from the
  spaxels above S/N 5 in the displayed image; masks are saved and loaded as FITS.

**Line fitting**
* 17 lines in 9 linesets: Hα + [NII], [OIII], [OII], Hβ, [SII], [OI], [SIII], HeI and Lyα.
  The lines of a set share their kinematics; the kinematics of the other sets are fixed to
  the first one fitted.
* Optional second kinematic component, identified as the fainter, the more offset or the
  broader one, with a "smart" mode that keeps it only when two components are really present.
* [NII], [OIII] and [OI] doublet ratios can be tied (Storey & Zeippen 2000). Start values,
  limits and fixed parameters per line; start values can also be read per spaxel from maps.
* Continuum from two side bands (percentile clipping and a polynomial), or a constant fitted
  together with the lines.
* Instrumental resolution from the sky lines in the variance spectrum, from the KMOS
  calibration polynomials of the nearest epoch, from the header or the MUSE manual, or from
  spectral templates.
* Moments (flux, velocity, dispersion) as an alternative to Gaussian fits.

**Uncertainties and quality**
* Errors from the noise cube (optionally rescaled by the reduced χ²), from bootstrap cubes or
  from three Monte Carlo schemes. The parameter distributions can be saved as plots and tables.
* Automatic OK / BAD flags from S/N and velocity-error thresholds, manual flagging, and
  FIT ADJ ALL, which refits bad spaxels from the solutions of their good neighbours.

**Throughput and automation**
* FIT ALL and FIT ADJ ALL on parallel worker processes (`--nproc`), with a progress bar and
  an Interrupt button.
* Batch mode without GUI, and a Python API for scripted analyses.
* Results as FITS cubes with one plane per parameter, error and flag, described in the
  header; sessions saved and reopened as `.kvz` files; spectra, images and cubes exported to
  FITS. IDL `.sav` sessions can be opened read-only.
* `kubeviz-compare` compares two results files plane by plane.

## Installation

kubeviz needs Python 3.10 or newer. The engine depends on numpy, scipy and astropy; the GUI
adds PyQt6 and pyqtgraph. A dedicated conda or virtual environment is recommended.

```bash
git clone https://github.com/matteofox/kubeviz.git
cd kubeviz
pip install ".[gui]"          # or: pip install .  for the engine and batch mode only
```

For development, install in editable mode with the test dependencies:

```bash
pip install -e ".[gui,test]"
```

If the directory that contains this repository is on your `PYTHONPATH` (for example
`~/Software`), the repository folder hides the package and `kubeviz` fails with "cannot
import name `__version__` from 'kubeviz' (unknown location)". Install in compat editable
mode instead:

```bash
pip install -e ".[gui,test]" --config-settings editable_mode=compat
```

## Quick start

```bash
kubeviz --demo                                  # synthetic cube, to try the interface
kubeviz cube.fits --redshift 0.85 --lineset 1   # interactive session
kubeviz session.kvz                             # reopen a saved session
kubeviz --help                                  # all command line options
```

In the viewer, a left click moves the crosshair and shows that spaxel's spectrum, the wheel
zooms and a right drag changes the contrast. Check the redshift and the instrumental
resolution in the Line fitting panel, tick the lines to fit, then press FIT for the current
spaxel or FIT ALL for the whole field. Help → Instructions in the program describes every
control, and Help → Keyboard shortcuts lists the keys.

Batch mode fits every spaxel (or every mask given with `--spmask`) and writes the results
and a session file:

```bash
kubeviz cube.fits --redshift 0.85 --lineset 1 --batch --fit-all-lines --nproc 8 --outdir results/
```

## Using kubeviz from Python

```python
from kubeviz.session import start_session
from kubeviz.fitting.linefit import dofit
from kubeviz.fitting.fitall import fitall
from kubeviz.fitting.flags import autoflag
from kubeviz.io.results import saveres

state = start_session("cube.fits", redshift=0.85, lineset=1)
state.pndofit[:state.Nlines] = 1      # fit the narrow component of every line in the set
state.pcdofit[:state.Nlines] = 1      # and the continuum
fitall(state, nproc=4)                 # or dofit(state) for the current spaxel (state.col/row)
autoflag(state)
saveres(state, "cube_res.fits")
res = state.res                        # ResultSet: res.n[row, col] = [flag, dv, sigma, fluxes...]
```

## Contributing

Contributions are welcome as pull requests. Please run the test suite before submitting;
it uses small synthetic cubes and runs the GUI tests offscreen:

```bash
python -m pytest
```

Add your name to the contributors list below in your first pull request.

## Credits

kubeviz is developed by

* Matteo Fossati
* David J. Wilman

Contributors:

* *(add your name with your first pull request)*

The original IDL kubeviz was written by David J. Wilman, Matteo Fossati and Joris Gerssen,
based on cubeviz by Joris Gerssen.

## Citing kubeviz

If you use kubeviz in your work, please add the following sentence to the acknowledgements
of the resulting publications:

> This work made use of the KUBEVIZ software which is publicly available at
> https://github.com/matteofox/kubeviz/.

and cite [Fossati et al. (2016), MNRAS 455, 2028](http://adsabs.harvard.edu/abs/2016MNRAS.455.2028F).

## License

kubeviz is distributed under the GNU General Public License v3 or later (see `LICENSE`).
