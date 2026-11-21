"""The two-stage seismic analysis's stage 2 strength factor.

The factor multiplies each low-permeability slice's stage 2 undrained strength
(a screening method's strength reduction, Duncan, Wright & Brandon 2014, Table
10.1). Five checks:

* **F scales by exactly the factor.** In stage 2 every low-permeability slice
  carries c = tau_ff and phi = 0. When every slice is one of them, the strength on
  each base is a fixed number and every method's equilibrium involves it only as
  c/F, so F at factor f is exactly f times F at 1. Stage 1 doesn't change. Run on
  Duncan, Wright & Brandon's Fig. 10.7 slope (the Table 10.4 workbook) in every
  method, with both stage 2 strength choices.
* **Free-draining slices keep their drained strength; c/p slices are reduced.**
  On a slice table with a third of its slices stripped of d and psi and a third
  made c/p slices, the factor scales c on the d/psi and c/p slices only, and
  stage 1 doesn't change.
* **The yield coefficient gets it.** k_y with the factor is lower, and the
  two-stage analysis at that k_y with the same factor gives F = 1.
* **The option is checked.** A factor outside 0 < f <= 1, or a dict with a key
  other than 'strength' and 'factor', is refused.
* **The words.** The run label and the Log's stage line name a factor below 1, and
  say nothing about one at 1.

Runs in about ten seconds.
"""
import os
import sys
import warnings

import numpy as np

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, ROOT)

from xslope.advanced import (seismic2_factor, seismic2_strength, seismic2_tag,  # noqa: E402
                             seismic_stages_line, seismic_two_stage)
from xslope.fileio import load_slope_data  # noqa: E402
from xslope.slice import generate_slices  # noqa: E402

MODEL = os.path.join(ROOT, "docs", "verification", "files", "published", "dwb_t10_4_dry.xlsx")
METHODS = ("oms", "bishop", "janbu", "corps", "lowe", "spencer", "mprice")
FACTOR = 0.8


def _slices(soil="Soil 2 two-stage"):
    """The Table 10.4 slope, every zone the named soil, sliced on its circle."""
    sd = load_slope_data(MODEL)
    row = next(m for m in sd["materials"] if m.get("name") == soil)
    sd["materials"] = [dict(row) for _ in sd["materials"]]
    ok, out = generate_slices(sd, circle=sd["circles"][0], num_slices=40)
    if not ok:
        raise RuntimeError(out)
    return sd, out[0]


def _two_stage(df, method, seismic2):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return seismic_two_stage(df.copy(), method, debug_level=0, strength=seismic2)


def check_f_scales_by_the_factor():
    fails = []
    _sd, df = _slices()
    for strength in ("interpolated", "kc1"):
        for method in METHODS:
            ok1, r1 = _two_stage(df, method, strength)
            okf, rf = _two_stage(df, method, {"strength": strength, "factor": FACTOR})
            if not (ok1 and okf):
                fails.append(f"{method}, {strength}: no solution ({r1 if not ok1 else rf})")
                continue
            if not all(r1["seismic_low_k"]):
                fails.append(f"{method}, {strength}: not every slice took an undrained strength")
                continue
            want = FACTOR * r1["FS"]
            if abs(rf["FS"] - want) > 1e-6 * max(1.0, want):
                fails.append(f"{method}, {strength}: F at factor {FACTOR} = {rf['FS']:.6f}, "
                             f"not {FACTOR} x {r1['FS']:.6f} = {want:.6f}")
            if rf["seismic_stage1_FS"] != r1["seismic_stage1_FS"]:
                fails.append(f"{method}, {strength}: the factor changed stage 1 "
                             f"({r1['seismic_stage1_FS']:.6f} -> {rf['seismic_stage1_FS']:.6f})")
            if rf.get("seismic_factor") != FACTOR or r1.get("seismic_factor") != 1.0:
                fails.append(f"{method}, {strength}: the result records factor "
                             f"{rf.get('seismic_factor')!r} and {r1.get('seismic_factor')!r}")
    return fails


def check_free_draining_slices_keep_their_strength():
    fails = []
    _sd, df = _slices()
    if "cp_flag" not in df.columns:
        return ["the slice table has no cp_flag column"]
    k = np.arange(len(df)) % 3
    drained, cp = k == 0, k == 1
    df.loc[drained | cp, "d"] = 0.0
    df.loc[drained | cp, "psi"] = 0.0
    # a c/p slice as the slicer builds one: c the undrained strength, phi, c1, phi1 0
    df.loc[cp, "c"] = 900.0
    df.loc[cp, ["phi", "c1", "phi1"]] = 0.0
    df.loc[cp, "cp_flag"] = True
    one, scaled = df.copy(), df.copy()
    ok1, r1 = seismic_two_stage(one, "spencer", debug_level=0, strength="interpolated")
    okf, rf = seismic_two_stage(scaled, "spencer", debug_level=0,
                                strength={"strength": "interpolated", "factor": FACTOR})
    if not (ok1 and okf):
        return [f"spencer found no solution on the mixed slice table ({r1 if not ok1 else rf})"]
    if not np.allclose(scaled.loc[drained, "c"], df.loc[drained, "c1"]) or \
            not np.allclose(scaled.loc[drained, "phi"], df.loc[drained, "phi1"]):
        fails.append("the factor changed the strength of a slice with no d or psi")
    und = ~drained
    if not np.allclose(scaled.loc[und, "c"], FACTOR * one.loc[und, "c"]):
        fails.append(f"the d/psi and c/p slices' c is not {FACTOR} times its value at factor 1")
    if not np.allclose(scaled.loc[cp, "c"], FACTOR * 900.0):
        fails.append("the c/p slices' strength was not reduced")
    if rf["seismic_stage1_FS"] != r1["seismic_stage1_FS"]:
        fails.append("the factor changed stage 1 on the mixed slice table")
    return fails


def check_the_slicer_flags_cp_slices():
    """Both slice-table builders mark c/p slices: generate_slices, and
    update_slice_materials after a material switches to c/p."""
    from xslope.slice import update_slice_materials
    sd = load_slope_data(MODEL)
    row = next(m for m in sd["materials"] if m.get("name") == "Soil 2 R envelope")
    mc = [dict(row) for _ in sd["materials"]]
    cp = [dict(row, option="cp", c=500.0, cp=20.0, r_elev=30.0) for _ in sd["materials"]]
    fails = []
    for mats, want in ((mc, False), (cp, True)):
        s = dict(sd, materials=mats)
        ok, out = generate_slices(s, circle=s["circles"][0], num_slices=20)
        if not ok:
            return [f"generate_slices failed: {out}"]
        if bool(out[0]["cp_flag"].all()) != want or bool(out[0]["cp_flag"].any()) != want:
            fails.append(f"generate_slices: cp_flag is not {want} on every slice")
    s = dict(sd, materials=mc)
    ok, out = generate_slices(s, circle=s["circles"][0], num_slices=20, material_breakdown=True)
    try:
        df2 = update_slice_materials(out[0].copy(), cp)
    except ValueError as exc:
        return fails + [f"update_slice_materials refused the c/p switch: {exc}"]
    if not bool(df2["cp_flag"].all()):
        fails.append("update_slice_materials left cp_flag False after a switch to c/p")
    return fails


def check_the_yield_coefficient_gets_it():
    from xslope.seismic import yield_coefficient
    sd, _df = _slices()
    fails = []
    ky = {}
    for f in (1.0, FACTOR):
        s2 = "kc1" if f == 1.0 else {"strength": "kc1", "factor": f}
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            bundle = yield_coefficient(sd, "spencer", analysis="single_surface",
                                       surface="circular", num_slices=40, seismic2=s2,
                                       announce=False)
        ky[f] = float(bundle["yield"]["k_y"])
        if bundle["options"]["seismic2"] != s2:
            fails.append(f"the k_y bundle recorded seismic2 = {bundle['options']['seismic2']!r}")
    if not ky[FACTOR] < ky[1.0]:
        fails.append(f"k_y with factor {FACTOR} ({ky[FACTOR]:.4f}) is not below k_y at 1 "
                     f"({ky[1.0]:.4f})")
    sd_k = dict(sd, k_seismic=ky[FACTOR])
    ok, out = generate_slices(sd_k, circle=sd_k["circles"][0], num_slices=40)
    ok, res = _two_stage(out[0], "spencer", {"strength": "kc1", "factor": FACTOR})
    if not ok or abs(res["FS"] - 1.0) > 2e-3:
        fails.append(f"at k_y = {ky[FACTOR]:.4f} the two-stage F with factor {FACTOR} is "
                     f"{res['FS'] if ok else res}, not 1")
    return fails


def check_the_option_is_checked():
    fails = []
    for bad in ({"strength": "kc1", "factor": 0.0}, {"strength": "kc1", "factor": 1.2},
                {"strength": "kc1", "factor": "x"}, {"strength": "kc1", "reduction": 0.8},
                {"strength": False, "factor": 0.8}):
        try:
            seismic2_strength(bad)
            fails.append(f"seismic2 = {bad!r} was accepted")
        except ValueError:
            pass
    if seismic2_strength({"factor": 0.85}) != "interpolated" or \
            seismic2_factor({"factor": 0.85}) != 0.85:
        fails.append("a dict with only a factor is not the interpolated strength at that factor")
    if seismic2_factor("kc1") != 1.0 or seismic2_factor(True) != 1.0:
        fails.append("a plain seismic2 does not carry factor 1")
    return fails


def check_the_words():
    fails = []
    if seismic2_tag({"strength": "kc1", "factor": 0.8}) != \
            " (two-stage seismic, K_c = 1 envelope, times 0.80)":
        fails.append(f"run label: {seismic2_tag({'strength': 'kc1', 'factor': 0.8})!r}")
    if seismic2_tag("interpolated") != " (two-stage seismic)" or seismic2_tag(False) != "":
        fails.append("the run label names a factor that isn't there")
    line = seismic_stages_line({"seismic_stage1_FS": 2.0, "seismic_stage2_FS": 1.2,
                                "seismic_strength": "kc1", "seismic_factor": 0.85},
                               "the surface")
    if "times 0.85" not in line:
        fails.append(f"the Log's stage line doesn't name the factor: {line!r}")
    line1 = seismic_stages_line({"seismic_stage1_FS": 2.0, "seismic_stage2_FS": 1.2,
                                 "seismic_strength": "kc1"}, "the surface")
    if "times" in line1:
        fails.append(f"the Log's stage line names a factor at 1: {line1!r}")
    return fails


def run():
    """Failures as a list, for run_tests.py."""
    fails = []
    for check in (check_f_scales_by_the_factor, check_free_draining_slices_keep_their_strength,
                  check_the_slicer_flags_cp_slices,
                  check_the_yield_coefficient_gets_it, check_the_option_is_checked,
                  check_the_words):
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
    print("Stage 2 strength factor checks passed.")
