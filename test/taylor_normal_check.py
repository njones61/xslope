# Copyright 2025 Norman L. Jones
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""The Taylor series reports the normal reliability index beside the lognormal one.

Duncan, Wright and Brandon (2014) find some justification for either a normal or
a lognormal factor of safety and suggest computing both. Monte Carlo and the
response surface always reported both; the Taylor series, LEM and FEM, reported
only the lognormal beta_ln and its probability of failure. Every number here is
checked against a hand computation from the run's own F_MLV and sigma_F, with
the normal CDF taken from math.erfc rather than the engine's scipy:

  beta_N  = (F_MLV - 1) / sigma_F             P_f,N  = Phi(-beta_N)
  beta_LN = ln(F_MLV / sqrt(1 + V^2)) / sqrt(ln(1 + V^2)),  V = sigma_F / F_MLV
                                               P_f,LN = Phi(-beta_LN)

  A. LEM: reliability_taylor on VP28a's fixed circle (Bishop, search=False,
     3 parameters, 7 fixed-surface solves, well under a second). sigma_F is
     recombined from the per-parameter F+ and F-; the result carries
     beta_normal, pf_normal and pf_lognormal with the lognormal keys unchanged
     in meaning (prob_failure = pf_lognormal = Phi(-beta_ln), reliability =
     1 - prob_failure); the Log text (debug_level=1 prints both the running
     lines and the summary table) has the normal lines beside the lognormal.
  B. FEM: reliability_fem's own assembly, with the SSRM solve replaced by a
     linear stand-in F = F0 + sum(k_i x_i). The perturbation loop, the
     combination and the result dict are the real code; the stand-in makes
     F_MLV and sigma_F closed-form, so no mesh is built and nothing is solved.
  C. Studio: the Taylor summary dialog and status text (LEM and FEM) show
     beta (normal) and Pf (normal) beside the lognormal values, as the Monte
     Carlo summary does.

Runs on its own (``python test/taylor_normal_check.py``) or as a suite row.
"""

import contextlib
import io
import math
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(_HERE)
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from tools import studio_settings                      # noqa: E402

#: VP28a (Slide2 verification #28, case a): a fixed Bishop circle with sigmas on
#: three parameters. beta_N (0.790) and beta_LN (0.772) differ in the third
#: decimal, so a swapped pair shows in the printed text as well as in the dict.
VP28A = os.path.join(_REPO, "docs/verification/files/rocscience/vp028a.xlsx")

#: Exact agreement: the engine and the hand computation do the same arithmetic.
TOL = 1e-9


def _phi(x):
    """The standard normal CDF, from math.erfc (independent of scipy)."""
    return 0.5 * math.erfc(-x / math.sqrt(2.0))


def _hand(f_mlv, sigma_f):
    """Both conventions from F_MLV and sigma_F alone."""
    v = sigma_f / f_mlv
    b_n = (f_mlv - 1.0) / sigma_f
    b_ln = math.log(f_mlv / math.sqrt(1.0 + v * v)) / math.sqrt(math.log(1.0 + v * v))
    return {"beta_normal": b_n, "pf_normal": _phi(-b_n),
            "beta_ln": b_ln, "pf_lognormal": _phi(-b_ln)}


def _close(a, b):
    try:
        return abs(float(a) - float(b)) <= TOL * max(1.0, abs(float(b)))
    except (TypeError, ValueError):
        return False


def _check_dict(res, want, label):
    """The result dict's two conventions against the hand values."""
    fails = []
    for key in ("beta_normal", "pf_normal", "beta_ln", "pf_lognormal"):
        if key not in res:
            fails.append(f"{label}: the result has no '{key}' (keys: "
                         f"{', '.join(sorted(k for k in res if not k.startswith('_')))})")
        elif not _close(res[key], want[key]):
            fails.append(f"{label}: {key} = {res[key]!r}, hand computation gives "
                         f"{want[key]!r}")
    # The lognormal keys keep their meaning.
    if not _close(res.get("prob_failure"), want["pf_lognormal"]):
        fails.append(f"{label}: prob_failure = {res.get('prob_failure')!r} is not the "
                     f"lognormal Phi(-beta_ln) = {want['pf_lognormal']!r}")
    if not _close(res.get("reliability"), 1.0 - want["pf_lognormal"]):
        fails.append(f"{label}: reliability = {res.get('reliability')!r} is not "
                     f"1 - Phi(-beta_ln) = {1.0 - want['pf_lognormal']!r}")
    return fails


def _check_log(text, want, label, running=False):
    """The printed summary carries the normal lines beside the lognormal ones."""
    lines = [f"β_normal: {want['beta_normal']:.3f}",
             f"β_ln: {want['beta_ln']:.3f}",
             f"Probability of failure (normal): {want['pf_normal'] * 100:.2f}%",
             f"Probability of failure (lognormal): {want['pf_lognormal'] * 100:.2f}%"]
    if running:
        lines += [f"β_normal = {want['beta_normal']:.4f}",
                  f"Probability of failure (normal) = {want['pf_normal'] * 100:.2f}%"]
    return [f"{label}: the Log text has no line '{s}'" for s in lines if s not in text]


def leg_lem():
    """A. reliability_taylor on VP28a's fixed circle."""
    from xslope.fileio import load_slope_data
    from xslope.reliability import reliability
    sd = load_slope_data(VP28A)
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        ok, res = reliability(sd, "bishop", circular=True, search=False, debug_level=1)
    if not ok:
        return [f"LEM: VP28a's Taylor run failed: {res}"], None
    fails = []
    # sigma_F recombined from the per-parameter F+ and F-.
    sigma = math.sqrt(sum(((p["F_plus"] - p["F_minus"]) / 2.0) ** 2
                          for p in res["param_info"]))
    if not _close(res["sigma_F"], sigma):
        fails.append(f"LEM: sigma_F = {res['sigma_F']!r}, the per-parameter F+/F- "
                     f"give {sigma!r}")
    want = _hand(float(res["F_MLV"]), sigma)
    fails += _check_dict(res, want, "LEM")
    fails += _check_log(buf.getvalue(), want, "LEM", running=True)
    return fails, res


#: The linear stand-in for the SSRM factor of safety in leg B: F = F0 + sum(k x).
_F0 = 0.2
_K = {"c": 0.02, "phi": 0.03, "gamma": -0.01}
_MAT = {"name": "clay", "option": "mc", "c": 20.0, "phi": 25.0, "gamma": 18.0,
        "sigma_c": 4.0, "sigma_phi": 2.0, "sigma_gamma": 0.9}


def _stand_in_fs(materials):
    m = materials[0]
    return _F0 + sum(k * float(m[p]) for p, k in _K.items())


def leg_fem():
    """B. reliability_fem's assembly with a closed-form stand-in for the SSRM."""
    import xslope.fem as fem
    from xslope.reliability import reliability_fem
    saved = (fem.build_fem_data, fem.solve_ssrm)

    def build(sd, mesh):
        return {"materials": [dict(m) for m in sd["materials"]]}

    def solve(fem_data, **kw):
        return {"converged": True, "FS": _stand_in_fs(fem_data["materials"]),
                "last_solution": None}

    fem.build_fem_data, fem.solve_ssrm = build, solve
    buf = io.StringIO()
    try:
        with contextlib.redirect_stdout(buf):
            ok, res = reliability_fem({"materials": [dict(_MAT)]}, mesh=object(),
                                      check_inputs=False)
    finally:
        fem.build_fem_data, fem.solve_ssrm = saved
    if not ok:
        return [f"FEM: the stand-in run failed: {res}"], None
    f_mlv = _stand_in_fs([_MAT])
    sigma = math.sqrt(sum((k * _MAT["sigma_" + p]) ** 2 for p, k in _K.items()))
    fails = []
    if not _close(res["F_MLV"], f_mlv):
        fails.append(f"FEM: F_MLV = {res['F_MLV']!r}, the stand-in gives {f_mlv!r}")
    if not _close(res["sigma_F"], sigma):
        fails.append(f"FEM: sigma_F = {res['sigma_F']!r}, closed form {sigma!r}")
    want = _hand(f_mlv, sigma)
    fails += _check_dict(res, want, "FEM")
    fails += _check_log(buf.getvalue(), want, "FEM")
    return fails, res


def leg_studio(results):
    """C. Studio's Taylor summary dialog and status text, LEM and FEM."""
    try:
        from studio.main_window import MainWindow
    except ImportError as e:
        print(f"  Studio not importable ({e}); leg C skipped")
        return []
    fails = []
    for app_mode, res in results:
        if res is None:
            continue
        want = _hand(float(res["F_MLV"]), float(res["sigma_F"]))
        label = f"Studio {app_mode.upper()}"
        try:
            summary = MainWindow._reliability_summary(res, "taylor", app_mode)
            status = MainWindow._reliability_status(res, "taylor")
        except KeyError as e:
            fails.append(f"{label}: the summary reads a key the result lacks: {e}")
            continue
        for s in (f"β (normal) = {want['beta_normal']:.3f}",
                  f"β (lognormal) = {want['beta_ln']:.3f}",
                  f"Pf (normal) = {want['pf_normal'] * 100:.2f}%",
                  f"Pf (lognormal) = {want['pf_lognormal'] * 100:.2f}%"):
            if s not in summary:
                fails.append(f"{label}: the summary dialog has no '{s}':\n{summary}")
        for s in (f"β_normal = {want['beta_normal']:.3f}",
                  f"β_ln = {want['beta_ln']:.3f}",
                  f"{want['pf_normal'] * 100:.2f}% (normal)",
                  f"{want['pf_lognormal'] * 100:.2f}% (lognormal)"):
            if s not in status:
                fails.append(f"{label}: the status text has no '{s}': {status}")
    return fails


def _guarded(name, fn, *args):
    try:
        return fn(*args)
    except Exception as e:                            # noqa: BLE001
        import traceback
        traceback.print_exc()
        return [f"{name}: raised {type(e).__name__}: {e}"], None


@studio_settings.isolated
def run():
    fails = []
    lem_fails, lem_res = _guarded("LEM", leg_lem)
    print(f"  A. LEM Taylor series, VP28a            {'ok' if not lem_fails else 'FAIL'}")
    fem_fails, fem_res = _guarded("FEM", leg_fem)
    print(f"  B. FEM Taylor series assembly          {'ok' if not fem_fails else 'FAIL'}")
    try:
        st_fails = leg_studio([("lem", lem_res), ("fem", fem_res)])
    except Exception as e:                            # noqa: BLE001
        st_fails = [f"Studio: raised {type(e).__name__}: {e}"]
    print(f"  C. Studio summary and status text      {'ok' if not st_fails else 'FAIL'}")
    fails += lem_fails + fem_fails + st_fails
    return fails


def main():
    print("Taylor series, normal and lognormal reliability index:")
    failures = run()
    if failures:
        print("\nFAILURES:")
        for f in failures:
            print(f"  - {f}")
        raise SystemExit(1)
    print("\nAll Taylor normal-index checks passed.")


if __name__ == "__main__":
    import matplotlib
    matplotlib.use("Agg")
    main()
