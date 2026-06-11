"""Run the LEM benchmark cases and capture per-method critical FOS.

For each case, runs an automated critical-surface search (circular or
non-circular) with each applicable method and prints a results table that maps
straight into docs/benchmark-run-plan.md.

Run from the repo root:  PYTHONPATH=. python3 benchmarks/run_lem.py
"""

from xslope.fileio import load_slope_data
from xslope.search import circular_search, noncircular_search

# Methods that apply to circular vs non-circular surfaces.
CIRC_METHODS = ["oms", "bishop", "janbu", "corps_engineers", "lowe_karafiath", "spencer"]
NONCIRC_METHODS = ["janbu", "corps_engineers", "lowe_karafiath", "spencer"]

CASES = [
    {"id": "LEM-1", "name": "ACADS Simple Slope",
     "file": "docs/lem/files/xslope_acads_simple.xlsx",
     "surface": "circular", "ref": 1.00},
    {"id": "LEM-2b", "name": "Arai & Tagyo Homogeneous",
     "file": "docs/lem/files/xslope_arai_tagyo.xlsx",
     "surface": "circular", "ref": 1.451},
    {"id": "LEM-2", "name": "ACADS Weak Layer",
     "file": "docs/lem/files/xslope_acads_weak_layer.xlsx",
     "surface": "non_circular", "ref": 1.26},
]


def run_case(case, num_slices=20):
    slope_data = load_slope_data(case["file"])
    methods = CIRC_METHODS if case["surface"] == "circular" else NONCIRC_METHODS
    out = {}
    for method in methods:
        try:
            if case["surface"] == "circular":
                fs_cache, converged, _, _ = circular_search(
                    slope_data, method, num_slices=num_slices, diagnostic=False)
            else:
                fs_cache, converged, _ = noncircular_search(
                    slope_data, method, diagnostic=False)
            if fs_cache:
                res = fs_cache[0].get("solver_result")
                fs = res.get("FS") if isinstance(res, dict) else None
            else:
                fs = None
        except Exception as e:
            fs = None
            print(f"    [{method}] ERROR: {e}")
        out[method] = fs
    return case, out


def main():
    for case in CASES:
        print(f"\n=== {case['id']}  {case['name']}  ({case['surface']}) "
              f"ref FOS={case['ref']} ===")
        _, out = run_case(case)
        print(f"  {'Method':<18} {'xslope FOS':>10} {'Ref':>8} {'Diff %':>8}")
        for method, fs in out.items():
            if fs is None:
                print(f"  {method:<18} {'--':>10} {case['ref']:>8.3f} {'--':>8}")
            else:
                diff = 100.0 * (fs - case['ref']) / case['ref']
                print(f"  {method:<18} {fs:>10.3f} {case['ref']:>8.3f} {diff:>8.1f}")


if __name__ == "__main__":
    main()
