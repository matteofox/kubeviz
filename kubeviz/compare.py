"""Compare two kubeviz results files plane by plane (IDL vs Python parity check).

Usage::

    kubeviz-compare idl_res.fits python_res.fits [--tol 0.05]

For every plane the report lists the number of spaxels valid in both files, the
median and maximum absolute difference, and, when a +1 sigma error plane exists for
that parameter, the median of |difference| / error and the fraction of spaxels that
agree within one sigma.
"""
from __future__ import annotations

import argparse
import sys

import numpy as np
from astropy.io import fits

from .constants import FIT_GAUSS


def _plane_names(hdr, nplanes):
    # IDL wrote double spaces ("narrow line  dv"); normalise for matching
    return [" ".join(str(hdr.get(f"P{i + 1}", f"plane {i}")).split()) for i in range(nplanes)]


def _error_plane_index(names, i):
    """Index of the '+1 sig' plane matching value plane ``i`` (Gaussian layout) or None."""
    name = names[i]
    if "1-sig" in name or name.startswith("Reduced"):
        return None
    if name.endswith("flag"):
        return None
    target = name.replace("narrow line ", "narrow line plus 1-sig ").replace("broad line ", "broad line plus 1-sig ") \
        .replace("continuum ", "continuum plus 1-sig ")
    try:
        return names.index(target)
    except ValueError:
        return None


def _align_trim(a, ha, b, hb, problems):
    """Cut the larger frame to the trimmed one using STARTCOL/STARTROW when the two
    results cover different portions of the same cube."""
    if a.ndim != 3 or b.ndim != 3 or a.shape[1:] == b.shape[1:]:
        return a, b
    def cut(big, small, hsmall):
        c0 = int(hsmall.get("STARTCOL", 1)) - 1
        r0 = int(hsmall.get("STARTROW", 1)) - 1
        nr, nc = small.shape[1:]
        if r0 + nr <= big.shape[1] and c0 + nc <= big.shape[2]:
            problems.append(f"trimmed comparison: rows {r0}:{r0 + nr}, cols {c0}:{c0 + nc} of the larger frame")
            return big[:, r0:r0 + nr, c0:c0 + nc]
        return big
    if a.shape[1] * a.shape[2] > b.shape[1] * b.shape[2]:
        return cut(a, b, hb), b
    return a, cut(b, a, ha)


def compare_results(file_a: str, file_b: str, tol: float = 0.05):
    """Return ``(rows, summary)`` where rows is a list of dicts, one per plane."""
    a, ha = fits.getdata(file_a, header=True)
    b, hb = fits.getdata(file_b, header=True)
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    problems = []
    for key in ("LINESET", "FITTYPE", "RESTYPE", "ERMETHOD", "Z_IN", "SPATSMTH", "SPECSMTH"):
        if ha.get(key) != hb.get(key):
            problems.append(f"{key}: {ha.get(key)} vs {hb.get(key)}")
    # a trimmed run (STARTCOL/STARTROW, 1-based) is compared with the matching cut of the other file
    a, b = _align_trim(a, ha, b, hb, problems)
    if a.shape != b.shape:
        nmin = min(a.shape[0], b.shape[0])
        problems.append(f"shape {a.shape} vs {b.shape}; comparing the first {nmin} planes")
        if a.shape[1:] != b.shape[1:]:
            raise ValueError("spatial dimensions differ, cannot compare")
        a, b = a[:nmin], b[:nmin]
    names = _plane_names(ha, a.shape[0])
    gauss = int(ha.get("FITTYPE", FIT_GAUSS)) == FIT_GAUSS
    rows = []
    for i in range(a.shape[0]):
        pa, pb = a[i], b[i]
        both = np.isfinite(pa) & np.isfinite(pb)
        diff = pb[both] - pa[both]
        row = dict(plane=i + 1, name=names[i], n_both=int(both.sum()), n_only_a=int((np.isfinite(pa) & ~np.isfinite(pb)).sum()),
                   n_only_b=int((np.isfinite(pb) & ~np.isfinite(pa)).sum()), med_diff=np.nan, max_abs=np.nan,
                   med_pull=np.nan, frac_1sig=np.nan, frac_within_tol=np.nan)
        if diff.size:
            row["med_diff"] = float(np.median(diff))
            row["max_abs"] = float(np.max(np.abs(diff)))
            scale = np.abs(pa[both])
            with np.errstate(all="ignore"):
                rel = np.abs(diff) / np.where(scale > 0, scale, np.nan)
            row["frac_within_tol"] = float(np.nanmean(rel <= tol)) if np.any(np.isfinite(rel)) else np.nan
            ei = _error_plane_index(names, i) if gauss else None
            if ei is not None:
                err = np.abs(a[ei][both])
                with np.errstate(all="ignore"):
                    pull = np.abs(diff) / np.where(err > 0, err, np.nan)
                if np.any(np.isfinite(pull)):
                    row["med_pull"] = float(np.nanmedian(pull))
                    row["frac_1sig"] = float(np.nanmean(pull <= 1.0))
        rows.append(row)
    summary = dict(problems=problems, file_a=file_a, file_b=file_b)
    return rows, summary


def format_report(rows, summary) -> str:
    out = [f"A: {summary['file_a']}", f"B: {summary['file_b']}"]
    for p in summary["problems"]:
        out.append(f"[WARNING] {p}")
    out.append("")
    out.append(f"{'plane':>5} {'name':<44} {'Nboth':>6} {'onlyA':>5} {'onlyB':>5} {'med(B-A)':>11} {'max|B-A|':>11} {'med|d|/err':>10} {'f(<1sig)':>8} {'f(<tol)':>8}")
    for r in rows:
        out.append(f"{r['plane']:>5} {r['name'][:44]:<44} {r['n_both']:>6} {r['n_only_a']:>5} {r['n_only_b']:>5} "
                   f"{r['med_diff']:>11.4g} {r['max_abs']:>11.4g} {r['med_pull']:>10.3g} {r['frac_1sig']:>8.3f} {r['frac_within_tol']:>8.3f}")
    return "\n".join(out)


def main(argv=None) -> int:
    p = argparse.ArgumentParser(prog="kubeviz-compare", description=__doc__.split("\n\n")[0])
    p.add_argument("file_a", help="reference results file (e.g. IDL)")
    p.add_argument("file_b", help="results file to compare (e.g. Python)")
    p.add_argument("--tol", type=float, default=0.05, help="relative tolerance for f(<tol) (default 0.05)")
    args = p.parse_args(argv)
    try:
        rows, summary = compare_results(args.file_a, args.file_b, tol=args.tol)
    except Exception as exc:
        print(f"[ ERROR ] {exc}", file=sys.stderr)
        return 1
    print(format_report(rows, summary))
    return 0


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
