"""Checks for continuing a strength reduction run with a higher iteration limit.

What this file locks:

  1. THE BOOKKEEPING, on a stand-in solver (no finite element solve). A run
     keeps the end state of exactly the trials a higher limit could still decide
     -- undecided at the limit, or counted failed while still slowing -- and no
     other; a run whose top trial failed keeps nothing and offers nothing. A
     continuation hands each kept state back to the trial at its F, continues
     the top trial first and walks up through the trials above it that did not
     stand, then bisects; every trial is recorded once, a continued one in
     place with ``resumed_from``; the wall time is both parts; the closing
     summary says where the run was continued from and to what limit, in the
     Run dialog's words; and the bracket is the one a fresh run at the new
     limit finds. A run that cannot be continued is refused in plain words: a
     real failure at the top, a limit not above the one used, and a result that
     has been pickled (the states live in the session only).

  2. THE TRIAL, on Griffiths & Lane Example 1 (coarse, seconds). A trial
     stopped at a small limit and continued to a larger one ends bit-identical
     to the same trial run straight through at the larger limit -- plain loop,
     with the corrector, and accelerated.

  3. THE RUN, on the same model. A run at a small limit continued to a larger
     one reports the bracket and the lower-bound flag of a fresh run at the
     larger limit, and a continuation that ends undecided again can be
     continued again.

  4. STUDIO (offscreen). "Continue with a higher limit…" is on the FEM results
     toolbar only while the displayed run can be continued; pressing it runs
     the continuation in the worker, the results are replaced, and the Run FEM
     dialog takes the new Max iterations per trial, which it now accepts up to
     ten million.

  5. SLOW (``run_slow``; skipped with --quick). FEM-1 and the FEM-3 wall end on
     a failure, keep nothing and offer nothing; the FEM-3 geogrid wall run at
     100,000 iterations and continued to 1,000,000 reports the bracket the
     fresh million-iteration run found, FS 1.996 on [1.9922, 2.0].

Run directly:  PYTHONPATH=. python3 test/ssrm_resume_check.py [--slow]
Exits non-zero on any failure.
"""
import contextlib
import copy
import io
import os
import pickle
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import matplotlib                                           # noqa: E402
matplotlib.use("Agg")

import numpy as np                                          # noqa: E402

import xslope.fem as fem                                    # noqa: E402

FAILURES = []
_HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def check(name, cond, detail=""):
    status = "PASS" if cond else "FAIL"
    print(f"  [{status}] {name}" + (f"  — {detail}" if detail else ""))
    if not cond:
        FAILURES.append(name)


_MODEL = {}


def _griffiths_coarse():
    """Griffiths & Lane Example 1 on the coarse tri6 mesh (the smallest SSRM
    model shipped), built once."""
    if "fd" not in _MODEL:
        from xslope.fileio import load_slope_data
        from xslope.mesh import get_material_polygons, build_mesh_from_polygons
        xlsx = os.path.join(_HERE, "docs", "fem", "files", "xslope_griffiths1.xlsx")
        with contextlib.redirect_stdout(io.StringIO()):
            slope_data = load_slope_data(xlsx)
            mesh = build_mesh_from_polygons(get_material_polygons(slope_data),
                                            target_size=8.0, element_type="tri6")
            _MODEL["fd"] = fem.build_fem_data(slope_data, mesh)
    return _MODEL["fd"]


def _quiet(fn, *a, **kw):
    with contextlib.redirect_stdout(io.StringIO()):
        return fn(*a, **kw)


def _same(a, b):
    return abs(float(a) - float(b)) <= 1e-12 * max(1.0, abs(float(b)))


# ===================== 1. the bookkeeping (stand-in solver) =====================

#: The stand-in slope: stands up to F = 1.62, where a trial needs
#: 200 + 20000 (F - 1) iterations to settle; within 0.025 of that strength a
#: trial the limit stops is still slowing (counted failed), below it undecided.
_FS_TRUE = 1.62


def _stand_in(calls):
    """A solve_fem with the stand-in slope's trials, recording every call."""
    def fake(fem_data, F=1.0, max_iterations=12000, max_iterations_ceiling=50000,
             _keep_resume=False, _resume_state=None, **_kw):
        limit = max(int(max_iterations), int(max_iterations_ceiling or 0))
        start = 0 if _resume_state is None else int(_resume_state["iterations"])
        calls.append({"F": float(F), "limit": limit, "start": start,
                      "state": _resume_state})
        sol = dict(converged=False, stable=False, verdict="FAILED", u_ratio=3.0,
                   u_growth=0.5, max_displacement=0.01 * F, stop_reading=None,
                   creep_reading=None, plateau_iteration=None,
                   budget_extensions=0, iterations=50, exit_reason="diverging",
                   softened_1d_elements=np.zeros(0, dtype=bool),
                   _resume_state=None, F=F)
        if F <= _FS_TRUE + 1e-12:
            need = int(200 + 20000 * (F - 1.0))
            if need <= limit:
                sol.update(converged=True, stable=True, verdict="CONVERGED",
                           u_ratio=None, u_growth=None, iterations=need,
                           exit_reason="converged")
            else:
                slowing = F > _FS_TRUE - 0.025
                sol.update(iterations=limit, verdict="AMBIGUOUS", u_growth=0.05,
                           exit_reason="iteration_cap" if slowing else "inconclusive")
                if _keep_resume:
                    sol["_resume_state"] = {"F": float(F), "iterations": limit,
                                            "iteration": limit - 1}
        return sol
    return fake


def _stand_in_run(calls, **kw):
    real = fem.solve_fem
    fem.solve_fem = _stand_in(calls)
    try:
        return _quiet(fem.solve_ssrm, _griffiths_coarse(), **kw)
    finally:
        fem.solve_fem = real


def _kept_Fs(result):
    return sorted(float(k) for k in ((result.get("resumable") or {})
                                     .get("states") or {}))


def check_bookkeeping():
    print("\n1. the bookkeeping, on a stand-in solver")
    common = dict(F_min=1.0, F_max=2.0, tolerance=0.01, failure_criterion="hybrid",
                  capture_failure_state=False)
    c1 = []
    r1 = _stand_in_run(c1, max_iterations=6000, max_iterations_ceiling=6000, **common)
    check("the run at 6,000 ends on an undecided top: a lower bound",
          r1.get("fs_is_lower_bound") and _same(r1["final_interval"][1], 1.296875),
          f"FS {r1.get('FS')} on {r1.get('final_interval')}")
    want = sorted(float(t["F"]) for t in r1["trials"]
                  if fem._trial_can_continue(t))
    check("it keeps the end state of exactly the trials a higher limit could decide",
          _kept_Fs(r1) == want and want == [1.296875, 1.3125, 1.375, 1.5],
          f"kept {_kept_Fs(r1)}, continuable {want}")
    check("no kept state for a trial that stood or failed",
          not any(t["exit_reason"] in ("converged", "diverging")
                  and any(_same(t["F"], k) for k in _kept_Fs(r1))
                  for t in r1["trials"]))
    check("ssrm_can_continue names the top trial",
          (fem.ssrm_can_continue(r1) or {}).get("F") == 1.296875)
    check("the states are not in the run record",
          "resumable" not in fem.ssrm_run_record(r1)
          and not any("state" in k for t in fem.ssrm_run_record(r1)["trials"]
                      for k in t))

    before = [dict(t) for t in r1["trials"]]
    kept = dict(r1["resumable"]["states"])
    c2 = []
    real = fem.solve_fem
    fem.solve_fem = _stand_in(c2)
    try:
        r2 = _quiet(fem.solve_ssrm, _griffiths_coarse(), resume=r1,
                    max_iterations=20000)
        c3 = []
        fem.solve_fem = _stand_in(c3)
        r3 = _quiet(fem.solve_ssrm, _griffiths_coarse(), max_iterations=20000,
                    max_iterations_ceiling=20000, **common)
    finally:
        fem.solve_fem = real
    check("the continuation finds the bracket a fresh run at 20,000 finds",
          r2.get("final_interval") == r3.get("final_interval")
          and r2.get("FS") == r3.get("FS")
          and r2.get("fs_is_lower_bound") == r3.get("fs_is_lower_bound"),
          f"continued {r2.get('final_interval')}, fresh {r3.get('final_interval')}")
    order = [c["F"] for c in c2]
    check("the top trial is continued first, then the trials above it that "
          "did not stand, lowest first",
          order[:4] == [1.296875, 1.3125, 1.375, 1.5], f"order {order}")
    check("each continued trial is handed its own kept state and starts at "
          "6,000 iterations",
          all(c["state"] is kept[c["F"]] and c["start"] == 6000 for c in c2[:4])
          and all(c["state"] is None and c["start"] == 0 for c in c2[4:]))
    check("the trials that stood or failed are not solved again",
          not any(_same(c["F"], F) for c in c2 for F in (1.0, 2.0, 1.25, 1.28125,
                                                          1.2890625)))
    Fs = [round(float(t["F"]), 12) for t in r2["trials"]]
    check("every trial is recorded once", len(Fs) == len(set(Fs)), f"{Fs}")
    ok_place = True
    for i, t in enumerate(before):
        if float(t["F"]) in (1.296875, 1.3125, 1.375, 1.5):
            n = r2["trials"][i]
            ok_place &= (_same(n["F"], t["F"]) and n.get("resumed_from") == 6000
                         and n["role"] == t["role"] and n["stable"]
                         and n["iterations"] > 6000)
        else:
            ok_place &= (r2["trials"][i] == t)
    check("a continued trial's record replaces the earlier one in place, with "
          "resumed_from = 6,000; the others are as they were", ok_place)
    check("the wall time is both parts",
          r2["elapsed_time"] >= r1["elapsed_time"] > 0.0,
          f"{r1['elapsed_time']:.3f} s then {r2['elapsed_time']:.3f} s")
    sentence = ("Continued from the run that stopped at F = 1.2969, with Max "
                "iterations per trial raised to 20,000.")
    check("the closing summary says where the run was continued from",
          sentence in r2["summary"], r2["summary"][-240:])
    low = r2["summary"].lower()
    check("in the Run dialog's words",
          not any(w in low for w in ("budget", "edge", "sweep", "verdict",
                                     "corrector", "resum", "state")))
    check("a run whose top trial failed keeps nothing and offers nothing",
          "resumable" not in r2 and fem.ssrm_can_continue(r2) is None)
    try:
        _quiet(fem.solve_ssrm, _griffiths_coarse(), resume=r3, max_iterations=90000)
        refused = None
    except ValueError as exc:
        refused = str(exc)
    check("continuing a run that ended on a real failure is refused, in plain words",
          refused is not None and "failed" in refused and "F = 1.6250" in refused,
          repr(refused))

    # The AMBIGUOUS-at-limit top: counted failed while still slowing.
    c4 = []
    r4 = _stand_in_run(c4, max_iterations=12200, max_iterations_ceiling=12200,
                       **common)
    top = fem.ssrm_can_continue(r4) or {}
    check("a top trial the limit stopped while still slowing can be continued",
          not r4.get("fs_is_lower_bound") and top.get("exit_reason") == "iteration_cap"
          and _kept_Fs(r4) == [1.6015625, 1.609375],
          f"FS {r4.get('FS')} top {top.get('F')} kept {_kept_Fs(r4)}")
    c5 = []
    fem.solve_fem = _stand_in(c5)
    try:
        r5 = _quiet(fem.solve_ssrm, _griffiths_coarse(), resume=r4,
                    max_iterations=20000)
    finally:
        fem.solve_fem = real
    check("... and continued, it finds the fresh run's bracket",
          r5.get("final_interval") == r3.get("final_interval")
          and [c["F"] for c in c5][:2] == [1.6015625, 1.609375],
          f"{r5.get('final_interval')} after {[c['F'] for c in c5]}")

    try:
        _quiet(fem.solve_ssrm, _griffiths_coarse(), resume=r1, max_iterations=6000)
        refused = None
    except ValueError as exc:
        refused = str(exc)
    check("a limit not above the one the run used is refused",
          refused is not None and "above the 6,000" in refused, repr(refused))
    thawed = pickle.loads(pickle.dumps(r1))
    why = fem.ssrm_continue_refusal(thawed)
    check("a pickled result carries no states and cannot be continued",
          thawed.get("resumable") == {} and fem.ssrm_can_continue(thawed) is None
          and why is not None and "session" in why, repr(why))
    check("a deep copy keeps the states (the session's own copy)",
          fem.ssrm_can_continue(copy.deepcopy(r1)) is not None)

    reading = {"rule": "slowing_refused"}
    cases = [({"exit_reason": "inconclusive", "stable": False}, True),
             ({"exit_reason": "yield_gate", "stable": False}, True),
             ({"exit_reason": "iteration_cap", "stable": False,
               "verdict": "AMBIGUOUS"}, True),
             ({"exit_reason": "iteration_cap", "stable": False, "verdict": "FAILED",
               "stop_reading": reading}, True),
             ({"exit_reason": "iteration_cap", "stable": False,
               "verdict": "FAILED"}, False),
             ({"exit_reason": "not_slowing", "stable": False}, False),
             ({"exit_reason": "diverging", "stable": False}, False),
             ({"exit_reason": "converged", "stable": True}, False),
             ({"exit_reason": "inconclusive", "stable": True,
               "verdict": "STABLE_STUCK"}, False)]
    check("which trials a higher limit could still decide",
          all(fem._trial_can_continue(t) is want for t, want in cases),
          str([fem._trial_can_continue(t) for t, _ in cases]))


# ===================== 2. the trial (Griffiths, seconds) =====================

def _trial_diff(a, b):
    keys = ("exit_reason", "verdict", "iterations", "converged", "stable",
            "u_ratio", "u_growth", "unbalanced_force_ratio", "budget_extensions",
            "plateau_iteration", "gate_deferrals", "stop_reading", "accelerate")
    bad = [k for k in keys if a.get(k) != b.get(k)]
    du = float(np.max(np.abs(a["displacements"] - b["displacements"])))
    dp = max(float(np.max(np.abs(a["plastic_strains"][i] - b["plastic_strains"][i])))
             for i in a["plastic_strains"])
    return bad, du, dp


def check_trial_continuation():
    print("\n2. a continued trial is the trial run straight through")
    fd = _griffiths_coarse()
    prep = fem._prepare_fem_model(fd)
    for F, n1, n2, solver, accel in ((1.38, 60, 150, "viscoplastic", None),
                                     (1.40, 60, 150, "viscoplastic", None),
                                     (1.38, 100, 400, "viscoplastic", None),
                                     (1.43, 60, 150, None, None),
                                     (1.38, 60, 400, "viscoplastic", True),
                                     (1.43, 100, 1000, None, True)):
        kw = dict(F=F, _prepared=prep, fem_solver=solver, accelerate=accel)
        a = fem.solve_fem(fd, max_iterations=n1, max_iterations_ceiling=n1,
                          _keep_resume=True, **kw)
        st = a.get("_resume_state")
        b = (None if st is None else
             fem.solve_fem(fd, max_iterations=n2, max_iterations_ceiling=n2,
                           _resume_state=st, **kw))
        c = fem.solve_fem(fd, max_iterations=n2, max_iterations_ceiling=n2, **kw)
        name = (f"F = {F}, {n1:,} -> {n2:,}, "
                f"{solver or 'default driver'}{', accelerated' if accel else ''}")
        if b is None:
            check(name, False, f"the first leg kept no state ({a['exit_reason']})")
            continue
        bad, du, dp = _trial_diff(b, c)
        check(name + ": bit-identical to the straight run",
              not bad and du == 0.0 and dp == 0.0
              and b["iterations"] > n1,
              f"{a['exit_reason']} -> {b['exit_reason']} at {b['iterations']:,}; "
              f"differs in {bad}, |du| {du:g}, |d evp| {dp:g}")
    a = fem.solve_fem(fd, F=1.38, max_iterations=60, max_iterations_ceiling=60,
                      _prepared=prep, fem_solver="viscoplastic")
    check("nothing is kept unless asked", a.get("_resume_state") is None)
    st = fem.solve_fem(fd, F=1.38, max_iterations=60, max_iterations_ceiling=60,
                       _prepared=prep, fem_solver="viscoplastic",
                       _keep_resume=True)["_resume_state"]
    try:
        fem.solve_fem(fd, F=1.38, max_iterations=60, max_iterations_ceiling=60,
                      _prepared=prep, fem_solver="viscoplastic", _resume_state=st)
        refused = False
    except ValueError:
        refused = True
    check("a limit the trial has already reached is refused", refused)


# ===================== 3. the run (Griffiths, seconds) =====================

def check_run_continuation():
    print("\n3. a continued run finds the fresh run's bracket")
    fd = _griffiths_coarse()
    common = dict(F_min=1.0, F_max=2.0, tolerance=0.01, failure_criterion="hybrid")
    again = None
    for solver, n1, n2 in (("viscoplastic", 60, 400), ("viscoplastic", 100, 1000),
                           (None, 60, 400)):
        r1 = _quiet(fem.solve_ssrm, fd, max_iterations=n1,
                    max_iterations_ceiling=n1, fem_solver=solver, **common)
        name = f"{solver or 'default driver'}, {n1:,} -> {n2:,}"
        if fem.ssrm_can_continue(r1) is None:
            check(name, False, "the first run cannot be continued: "
                  + str(fem.ssrm_continue_refusal(r1)))
            continue
        r2 = _quiet(fem.solve_ssrm, fd, resume=r1, max_iterations=n2)
        r3 = _quiet(fem.solve_ssrm, fd, max_iterations=n2,
                    max_iterations_ceiling=n2, fem_solver=solver, **common)
        check(name + ": the same bracket and answer as the fresh run",
              r2["final_interval"] == r3["final_interval"]
              and r2["FS"] == r3["FS"]
              and r2["fs_is_lower_bound"] == r3["fs_is_lower_bound"],
              f"continued {r2['FS']} on {r2['final_interval']}, fresh {r3['FS']} "
              f"on {r3['final_interval']}")
        check(name + ": the at-failure field is drawn as on a fresh run",
              ("failure_solution" in r2) == ("failure_solution" in r3))
        if r2["fs_is_lower_bound"] and again is None:
            again = (r2, n2)
    if again is not None:
        r2, n2 = again
        r4 = _quiet(fem.solve_ssrm, fd, resume=r2, max_iterations=4 * n2)
        check("a continuation that ends undecided again can be continued again",
              fem.ssrm_can_continue(r2) is not None
              and r4.get("resumed", {}).get("max_iterations") == 4 * n2,
              f"{r2['FS']} -> {r4['FS']} on {r4['final_interval']}")
    else:
        check("a continuation that ends undecided again can be continued again",
              False, "no continued run ended undecided")


# ===================== 4. Studio (offscreen) =====================

def check_studio():
    print("\n4. Studio: Continue with a higher limit…")
    from PySide6.QtWidgets import QApplication
    app = QApplication.instance() or QApplication([])
    import studio.main_window as mwmod
    from studio.dialogs import MAX_ITERATIONS_LIMIT, RunFemDialog
    fd = _griffiths_coarse()
    common = dict(F_min=1.0, F_max=2.0, tolerance=0.01, failure_criterion="hybrid")
    r1 = _quiet(fem.solve_ssrm, fd, max_iterations=60, max_iterations_ceiling=60,
                **common)
    r_done = _quiet(fem.solve_ssrm, fd, max_iterations=400,
                    max_iterations_ceiling=400, **common)

    def bundle(r):
        return {"fem_data": fd, "solution": r["last_solution"],
                "failure_solution": r.get("failure_solution"), "FS": r["FS"],
                "analysis": "ssrm", "fs_is_lower_bound": r["fs_is_lower_bound"],
                "meta": fem.ssrm_run_record(r, fd, {}), "ssrm_result": r}

    mw = mwmod.MainWindow()
    real_get_int = mwmod.QInputDialog.getInt
    try:
        mw._last_fem_opts = {"analysis": "ssrm", "max_iterations": 60,
                             "max_iterations_ceiling": 60, **common}
        mw.doc.results["fem_solution"] = bundle(r1)
        mw._show_fem_results()
        btn = mw.fem_continue_btn
        check("the button is on the FEM results toolbar for a run that can be "
              "continued", btn is not None and not btn.isHidden()
              and btn.text() == "Continue with a higher limit…"
              and btn.isEnabled(), None if btn is None else btn.text())
        mw.doc.results["fem_solution"] = bundle(r_done)
        mw._update_fem_details_action()
        check("... and hidden for one that cannot", btn.isHidden())
        mw.doc.results["fem_solution"] = bundle(r1)
        mw._update_fem_details_action()
        asked = {}

        def fake_get_int(parent, title, label, value, lo, hi, step):
            asked.update(title=title, label=label, value=value, lo=lo, hi=hi)
            return 400, True
        mwmod.QInputDialog.getInt = staticmethod(fake_get_int)
        mw.continue_fem()
        check("it asks for the new Max iterations per trial, five times the "
              "limit to start", asked.get("value") == 300 and asked.get("lo") == 61
              and "Max iterations per trial" in asked.get("label", ""),
              str({k: v for k, v in asked.items() if k != "label"}))
        t0 = time.time()
        while mw._fem_runner is not None and time.time() - t0 < 120:
            app.processEvents()
            time.sleep(0.02)
        for _ in range(10):
            app.processEvents()
        b2 = mw.doc.results.get("fem_solution") or {}
        r2 = b2.get("ssrm_result") or {}
        check("the continuation ran in the worker and replaced the results",
              r2 is not r1 and (r2.get("resumed") or {}).get("max_iterations") == 400
              and b2.get("FS") == r2.get("FS"),
              f"FS {b2.get('FS')}, resumed {r2.get('resumed')}")
        check("the Run FEM dialog takes the new value for the next run",
              mw._last_fem_opts.get("max_iterations") == 400
              and mw._last_fem_opts.get("max_iterations_ceiling") == 400
              and "continue_from" not in mw._last_fem_opts)
        dlg = RunFemDialog(mw, defaults={"max_iterations": 1_000_000,
                                         "max_iterations_ceiling": 1_000_000})
        check("the dialog accepts a million iterations per trial",
              MAX_ITERATIONS_LIMIT >= 1_000_000
              and dlg.max_iterations.value() == 1_000_000
              and dlg.max_iterations_ceiling.value() == 1_000_000)
        dlg.close()
    finally:
        mwmod.QInputDialog.getInt = real_get_int
        mw.close()


# ===================== 5. slow: FEM-1, the FEM-3 walls =====================

def _grid_fem03():
    import tools.make_block_wall_figures as mb
    from xslope.fileio import load_slope_data
    model = load_slope_data(mb.FEM03_GRID)
    return mb, fem.build_fem_data(model, mb._mesh(model, mb.FEM03_TARGET_SIZE))


def run_slow():
    """FEM-1 and the FEM-3 wall keep nothing; the geogrid wall continued from
    100,000 to 1,000,000 reports the fresh million-iteration bracket."""
    FAILURES.clear()
    import run_tests as RT
    print("\n5. slow: FEM-1 and the FEM-3 walls")
    tag = {"file": "docs/tutorials/files/xslope_ssrm_embankment.xlsx",
           "type": "fem_ssrm", "element_type": "tri6", "target_size": 3.5,
           "tolerance": 0.01, "f_min": 1.0, "f_max": 2.0}
    fd, kw, f_min, f_max, tol = _quiet(RT.build_fem_ssrm_case, tag)
    with RT._force_fast_kernel(fem, False):
        r = _quiet(fem.solve_ssrm, fd, F_min=f_min, F_max=f_max, tolerance=tol, **kw)
    check("FEM-1: FS 1.371 on [1.3672, 1.375], nothing kept, nothing offered",
          r["final_interval"] == (1.3671875, 1.375) and "resumable" not in r
          and fem.ssrm_can_continue(r) is None, f"{r['FS']} {r['final_interval']}")
    import tools.make_block_wall_figures as mb
    from xslope.fileio import load_slope_data
    model = load_slope_data(mb.FEM03_WALL)
    fd = _quiet(fem.build_fem_data, model, _quiet(mb._mesh, model, mb.FEM03_TARGET_SIZE))
    with RT._force_fast_kernel(fem, False):
        r = _quiet(fem.solve_ssrm, fd, F_min=mb.FEM03_F_MIN, F_max=mb.FEM03_F_MAX,
                   tolerance=mb.FEM03_TOLERANCE, failure_criterion=mb.FEM03_CRITERION,
                   max_iterations=mb.FEM03_MAX_ITERATIONS)
    check("FEM-3 wall: FS 1.137 on [1.1328, 1.1406], nothing kept, nothing "
          "offered", r["final_interval"] == (1.1328125, 1.140625)
          and "resumable" not in r and fem.ssrm_can_continue(r) is None,
          f"{r['FS']} {r['final_interval']}")
    mb, fd = _quiet(_grid_fem03)
    with RT._force_fast_kernel(fem, False):
        t0 = time.time()
        r1 = _quiet(fem.solve_ssrm, fd, F_min=mb.FEM03_F_MIN, F_max=mb.FEM03_F_MAX,
                    tolerance=mb.FEM03_TOLERANCE,
                    failure_criterion=mb.FEM03_CRITERION,
                    max_iterations=mb.FEM03_MAX_ITERATIONS)
        s1 = time.time() - t0
        check("FEM-3 geogrid wall at 100,000: FS ≥ 1.56, the top trial can be "
              "continued", r1["fs_is_lower_bound"] and fem.ssrm_can_continue(r1)
              is not None, f"{r1['FS']} {r1['final_interval']} in {s1:.0f} s")
        t1 = time.time()
        r2 = _quiet(fem.solve_ssrm, fd, resume=r1, max_iterations=1_000_000)
        s2 = time.time() - t1
    check("... continued to 1,000,000: FS 1.996 on [1.9922, 2.0], the fresh "
          "million-iteration run's bracket",
          r2["final_interval"] == (1.9921875, 2.0) and r2["FS"] == 1.99609375
          and not r2["fs_is_lower_bound"],
          f"{r2['FS']} {r2['final_interval']}; {s1:.0f} s + {s2:.0f} s")
    return list(FAILURES)


def run():
    FAILURES.clear()
    check_bookkeeping()
    check_trial_continuation()
    check_run_continuation()
    check_studio()
    return list(FAILURES)


if __name__ == "__main__":
    fails = run_slow() if "--slow" in sys.argv else run()
    print(f"\n{'FAILED' if fails else 'PASSED'}: {len(fails)} failure(s)")
    for f in fails:
        print(f"   - {f}")
    sys.exit(1 if fails else 0)
