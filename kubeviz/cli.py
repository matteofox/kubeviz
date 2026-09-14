"""Command line interface. Every keyword of the IDL ``kubeviz`` procedure is a flag::

    kubeviz cube.fits --redshift 2.225 --lineset 1 --smooth 2 --bootstrap boot.fits
    kubeviz cube.fits --redshift 0.85 --lineset 1 --batch --fit-all-lines
    kubeviz session.kvz
"""
from __future__ import annotations

import argparse
import sys

from . import __version__


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="kubeviz",
                                description="Interactively examine a FITS datacube (x,y,lambda) and fit emission lines.")
    p.add_argument("datafile", nargs="?", help="datacube FITS file or saved session (.kvz)")
    p.add_argument("--noisefile", help="separate noise cube file")
    p.add_argument("--ext", type=int, help="FITS extension of the data cube (default: auto)")
    p.add_argument("--noise-ext", type=int, dest="noise_ext", help="FITS extension of the noise cube (default: auto)")
    p.add_argument("--trim", help="portion of the cube to load, '[x1:x2,y1:y2,z1:z2]' (1-based, * wildcard)")
    p.add_argument("--fluxfac", type=float, help="factor to multiply the flux scale (default instrument dependent)")
    p.add_argument("--smooth", type=int, default=1, help="spatial smoothing kernel in spaxels (default 1)")
    p.add_argument("--specsmooth", type=int, default=1, help="spectral smoothing kernel in pixels (default 1)")
    p.add_argument("--transpose", action="store_true", help="transpose the spatial axes")
    p.add_argument("--vacuum", action="store_true", help="wavelength solution is in vacuum")
    p.add_argument("--logarithmic", action="store_true", help="wavelength solution in the header is logarithmic")
    p.add_argument("--waveunit", help="wavelength unit if not in CUNIT3 (ANGSTROMS, MICRONS, NM)")
    p.add_argument("--redshift", "-z", type=float, help="redshift of the target")
    p.add_argument("--do-mc-errors", type=int, dest="do_mc_errors", choices=[0, 1, 2, 3, 4],
                   help="error method: 0 noise cube, 1 bootstrap, 2 MC1, 3 MC2, 4 MC3")
    p.add_argument("--bootstrap", help="multi-extension bootstrap cubes file (implies --do-mc-errors 1)")
    p.add_argument("--nmontecarlo", type=int, help="number of Monte Carlo realisations (default 100)")
    p.add_argument("--plot-mc-pdf", action="store_true", dest="plot_mc_pdf", help="save plots of the MC distributions")
    p.add_argument("--save-mc-pdf", action="store_true", dest="save_mc_pdf", help="save the MC distributions in a FITS table")
    p.add_argument("--use-mc-noise", action="store_true", dest="use_mc_noise", help="use the noise cube computed from the MC cubes")
    p.add_argument("--scroll", action="store_true", help="force a scrollbar in the linefit window")
    p.add_argument("--lineset", type=int, default=0,
                   help="1 Ha+[NII], 2 [OIII], 3 [OII], 4 Hb, 5 [SII], 6 [OI], 7 [SIII], 8 HeI, 9 Lya (default 0: all)")
    p.add_argument("--debug", action="store_true", help="increase verbosity")
    p.add_argument("--mask-sn-thresh", type=float, dest="mask_sn_thresh", help="S/N threshold for autoflag (default 3)")
    p.add_argument("--mask-maxvelerr", type=float, dest="mask_maxvelerr", help="max velocity error for autoflag (default 50 km/s)")
    p.add_argument("--mom-thresh", type=float, dest="mom_thresh", help="flux threshold for moments (default 0.5)")
    p.add_argument("--spmask", help="spaxel mask FITS file to load")
    p.add_argument("--instr", help="instrument name (overrides the header)")
    p.add_argument("--band", help="band / filter")
    p.add_argument("--batch", action="store_true", help="run without GUI: fit all spaxels (or masks) and save")
    p.add_argument("--fit-all-lines", action="store_true", dest="fit_all_lines", help="(batch) fit all lines of the lineset")
    p.add_argument("--fix-ratios", action="store_true", dest="fix_ratios", help="fix the [NII], [OIII], [OI] line ratios")
    p.add_argument("--resfile", help="kubeviz results file to load")
    p.add_argument("--fittype", choices=["gauss", "moments"], default="gauss", help="fit type (default gauss)")
    p.add_argument("--outdir", help="output directory (default: current directory)")
    p.add_argument("--logfile", help="redirect the terminal output to this file")
    p.add_argument("--fitpars", type=float, nargs="+",
                   help="6 or 8 numbers: [fitrange blue red] cont minoff maxoff minperc maxperc order mode")
    p.add_argument("--momwindowfile", help="(beta) 2-plane cube with velocity centre/width for the moments window")
    p.add_argument("--gaussinitfile", help="(beta) cube of initial guesses for the Gaussian parameters")
    p.add_argument("--scalenoiseerr", action="store_true", help="rescale the formal errors to chi-square/dof = 1")
    p.add_argument("--version", action="version", version=f"kubeviz {__version__}")
    return p


def main(argv=None) -> int:
    from .session import batchmode, start_session

    args = build_parser().parse_args(argv)
    if args.datafile is None:
        try:
            from .gui.app import run_gui
        except ImportError as exc:
            print(f"[ ERROR ] no datacube file given and the GUI is not available ({exc})", file=sys.stderr)
            return 2
        return run_gui(None, vars(args))

    kwargs = {k: v for k, v in vars(args).items() if k not in ("datafile", "batch", "fit_all_lines")}
    if args.fitpars is not None and len(args.fitpars) not in (6, 8):
        print("[ ERROR ] --fitpars needs 6 or 8 values", file=sys.stderr)
        return 2
    try:
        state = start_session(args.datafile, **kwargs)
    except (FileNotFoundError, ValueError) as exc:
        print(f"[ ERROR ] {exc}", file=sys.stderr)
        return 1

    if args.batch:
        try:
            resfile, sessfile = batchmode(state, redshift=args.redshift, fit_all_lines=args.fit_all_lines)
        except ValueError as exc:
            print(f"[ ERROR ] {exc}", file=sys.stderr)
            return 1
        print(f"[KUBEVIZ] Results: {resfile}")
        print(f"[KUBEVIZ] Session: {sessfile}")
        return 0

    try:
        from .gui.app import run_gui
    except ImportError as exc:
        print(f"[ ERROR ] The GUI requires PyQt6 and pyqtgraph ({exc}). Use --batch for headless operation.",
              file=sys.stderr)
        return 2
    return run_gui(state, vars(args))


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
