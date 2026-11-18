"""The finite element yield coefficient's bisection, on a stand-in solver.

``xslope.seismic.fem_yield_coefficient`` bisects k at full strength, one
``solve_fem`` per trial. A real trial takes minutes, so these checks replace
``xslope.fem.solve_fem`` with a function that stands below a chosen k and fails
above it, and check what the bisection does with the answers:

* it closes on the chosen k to within its tolerance, starting from no seismic
  force with no failing edge given, and from a bracket given;
* every trial carries the seismic force in the direction asked, with its sign,
  and nothing reduced (F = 1);
* an undecided trial is not counted as a failure: where one is the top of the
  final bracket, k_y is the bracket's bottom and is flagged as a lower bound;
* a slope that fails with no seismic force has no yield coefficient, and a model
  with no direction and an unsigned k is refused.

Runs in well under a second.
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

import xslope.fem as fem  # noqa: E402
from xslope.search import AnalysisError  # noqa: E402
from xslope.seismic import fem_yield_coefficient  # noqa: E402


class _StandIn:
    """solve_fem that stands for |k| <= k_true, fails above it, and ends
    'inconclusive' inside ``undecided`` (a (lo, hi) range of |k|)."""

    def __init__(self, k_true, undecided=None):
        self.k_true = k_true
        self.undecided = undecided
        self.calls = []

    def __call__(self, fem_data, F=1.0, **kwargs):
        k = float(fem_data["k_seismic"])
        self.calls.append((k, F, kwargs.get("k0")))
        if self.undecided and self.undecided[0] <= abs(k) <= self.undecided[1]:
            return {"stable": False, "verdict": "INCONCLUSIVE",
                    "exit_reason": "inconclusive", "iterations": 50000}
        stands = abs(k) <= self.k_true
        return {"stable": stands,
                "verdict": "CONVERGED" if stands else "FAILED",
                "exit_reason": "converged" if stands else "disp_limit",
                "iterations": 1000, "max_displacement": 0.01 if stands else 5.0}


def _with(stand_in, **kwargs):
    real = fem.solve_fem
    fem.solve_fem = stand_in
    try:
        return fem_yield_coefficient({"k_seismic": 0.0}, **kwargs)
    finally:
        fem.solve_fem = real


def check_it_closes_on_the_answer():
    fails = []
    s = _StandIn(0.3137)
    r = _with(s, direction="left", tol=0.002)
    if not (r["k_stand"] <= 0.3137 < r["k_fail"]):
        fails.append(f"the final bracket {r['k_stand']}..{r['k_fail']} misses 0.3137")
    if not r["k_fail"] - r["k_stand"] <= 0.002:
        fails.append(f"the final bracket is {r['k_fail'] - r['k_stand']:.4f} wide")
    if abs(r["k_y"] - 0.3137) > 0.001:
        fails.append(f"k_y = {r['k_y']:.4f}, not within 0.001 of 0.3137")
    if r["lower_bound"]:
        fails.append("a bracket closed by decided trials is flagged a lower bound")
    if any(k > 0 for k, _, _ in s.calls):
        fails.append("direction='left' applied a positive k")
    if any(F != 1.0 for _, F, _ in s.calls):
        fails.append("a trial reduced the strength")
    if any(k0 != 1.0 for _, _, k0 in s.calls):
        fails.append("a trial did not start from the K0 = 1 field")
    s = _StandIn(0.62)
    r = _with(s, direction="right", tol=0.005)
    if not (r["k_stand"] <= 0.62 < r["k_fail"]):
        fails.append(f"doubling past the start missed 0.62: {r['k_stand']}..{r['k_fail']}")
    if any(k < 0 for k, _, _ in s.calls):
        fails.append("direction='right' applied a negative k")
    s = _StandIn(0.40)
    r = _with(s, direction="left", k_stand=0.35, k_fail=0.45, tol=0.005)
    if not (r["k_stand"] <= 0.40 < r["k_fail"]):
        fails.append(f"a given bracket closed off 0.40: {r['k_stand']}..{r['k_fail']}")
    if abs(s.calls[0][0]) != 0.35 or abs(s.calls[1][0]) != 0.45:
        fails.append("a given bracket's edges were not solved first")
    return fails


def check_an_undecided_top_is_a_lower_bound():
    fails = []
    s = _StandIn(0.50, undecided=(0.30, 0.50))
    r = _with(s, direction="left", tol=0.005)
    if not r["lower_bound"]:
        fails.append("an undecided top of the bracket is not flagged a lower bound")
    if r["k_y"] != r["k_stand"]:
        fails.append("a lower-bound k_y is not the bracket's bottom")
    if r["k_stand"] >= 0.30:
        fails.append(f"an undecided trial was counted as standing (k_stand {r['k_stand']})")
    return fails


def check_refusals():
    fails = []
    try:
        _with(_StandIn(-1.0), direction="left", tol=0.05)
        fails.append("a slope failing with no seismic force got a yield coefficient")
    except AnalysisError:
        pass
    try:
        _with(_StandIn(0.3))
        fails.append("no direction and an unsigned k was not refused")
    except ValueError:
        pass
    return fails


def run():
    """Failures as a list, for run_tests.py."""
    fails = []
    for check in (check_it_closes_on_the_answer,
                  check_an_undecided_top_is_a_lower_bound, check_refusals):
        try:
            fails += check()
        except Exception as exc:          # a crash is a failure, not a pass
            fails.append(f"{check.__name__} raised {exc!r}")
    return fails


if __name__ == "__main__":
    failures = run()
    if failures:
        print("FAILURES:")
        for f in failures:
            print(f"  - {f}")
        raise SystemExit(1)
    print("Finite element yield coefficient checks passed.")
