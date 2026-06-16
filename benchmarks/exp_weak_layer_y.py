"""Experiment: sensitivity of the ACADS weak-layer non-circular FOS to the
starting elevation (y) of the two 'Horiz' sliding-plane points.

The non-circular search moves 'Horiz' points horizontally only (dy=0), so the
y we seed them at is locked in for the whole search. The weak layer spans
y=26.5 (base) to y=27.0 (top). Default seed is y=26.75 (center).

Run from repo root:  PYTHONPATH=. python3 benchmarks/exp_weak_layer_y.py
"""
import io
import contextlib

from xslope.fileio import load_slope_data
from xslope.search import noncircular_search

FILE = "docs/lem/files/xslope_acads_weak_layer.xlsx"
METHOD = "spencer"

# Sweep sliding-plane elevation from base of weak layer up to top.
Y_VALUES = [26.50, 26.55, 26.60, 26.65, 26.70, 26.75, 26.80, 26.90, 27.00]


def run_at_y(y_slide):
    slope_data = load_slope_data(FILE)
    # Override y for the two interior 'Horiz' points; leave Free endpoints alone.
    for p in slope_data['non_circ']:
        if p['Movement'] == 'Horiz':
            p['Y'] = y_slide
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        fs_cache, converged, _ = noncircular_search(
            slope_data, METHOD, num_slices=50, diagnostic=False)
    if fs_cache:
        res = fs_cache[0].get("solver_result")
        fs = res.get("FS") if isinstance(res, dict) else None
        # recover the converged surface x-coords of the horiz points
        return fs, converged
    return None, converged


def main():
    print(f"ACADS Weak Layer — {METHOD} FOS vs. starting sliding-plane y")
    print(f"weak layer: base y=26.5, top y=27.0 ; ref FOS = 1.26\n")
    print(f"  {'y_start':>8} {'FOS':>8} {'converged':>10}")
    for y in Y_VALUES:
        fs, conv = run_at_y(y)
        fs_s = f"{fs:.4f}" if fs is not None else "--"
        print(f"  {y:>8.2f} {fs_s:>8} {str(conv):>10}")


if __name__ == "__main__":
    main()
