"""Checks for the FEM solver's guards: the places where a defect in the code must
RAISE instead of being turned into a refusal, a fallback or a default, and the
places where the solver proceeds but must SAY so.

What this file locks:

  1. THE CORRECTOR'S CATCH. A numerical failure inside a Newton corrector attempt
     (LinAlgError, a singular factorization, an overflow) is a refusal, and the
     refusal names the exception's class. A defect in the code (NameError,
     KeyError, AttributeError ...) propagates out of solve_fem. The hold test
     reads the same way.

  2. THE OTHER BROAD CATCHES. A corrupt meta sidecar is ignored with a warning,
     an unrecognized unit system gives no unit, a singular stiffness still raises
     from the factorization, and an element whose stiffness cannot be built stops
     the assembly instead of leaving a hole in K.

  3. NO UNREACHABLE CODE, NO UNBOUND NAME. No function in fem.py has statements
     after a return at the top of its body, and `_nr_joint_slip` handed a state
     with no joint group raises a plain FemInvariantError.

  4. THE PREPARED MODEL. A solve that reuses a prepared model with different
     build options, or after one of its shared arrays changed, raises. A
     continuation of an SSRM run on a model whose content changed (a different
     dict or the same dict edited in place) raises and names what changed.

  5. CARRIED SEEDS. A post-peak bar set, or a Newton seed state, of the wrong size
     raises and names the sizes.

  6. MONOTONICITY. A run whose trials stood above a failure says so and records
     the pair in ``nonmonotone_trials``; its answer is the one it would have been.

  7. THE IN-SITU STATE. A trial that writes into the in-situ state every trial
     starts from is caught before the run returns.

Run directly:  PYTHONPATH=. python3 test/fem_guards_check.py
Each check is seconds on Griffiths & Lane Example 1's coarse mesh; section 1's
integration leg is the longest (a few seconds).
"""
import ast
import contextlib
import io
import json
import os
import sys
import tempfile

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _ROOT)

import numpy as np                                          # noqa: E402

import xslope.fem as fem                                    # noqa: E402

FAILURES = []
_MODEL = {}


def check(label, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {label}" + (f"  — {detail}" if detail else ""))
    if not ok:
        FAILURES.append(label)


def _quiet(fn, *a, **kw):
    with contextlib.redirect_stdout(io.StringIO()):
        return fn(*a, **kw)


def _griffiths_coarse():
    """Griffiths & Lane Example 1 on the coarse tri6 mesh, built once."""
    if "fd" not in _MODEL:
        from xslope.fileio import load_slope_data
        from xslope.mesh import get_material_polygons, build_mesh_from_polygons
        xlsx = os.path.join(_ROOT, "docs", "fem", "files", "xslope_griffiths1.xlsx")
        with contextlib.redirect_stdout(io.StringIO()):
            slope_data = load_slope_data(xlsx)
            mesh = build_mesh_from_polygons(get_material_polygons(slope_data),
                                            target_size=8.0, element_type="tri6")
            _MODEL["fd"] = fem.build_fem_data(slope_data, mesh)
    return _MODEL["fd"]


def _raises(fn, exc_type):
    """(raised, the exception or None)."""
    try:
        _quiet(fn)
    except exc_type as e:                   # noqa: PERF203
        return True, e
    except Exception as e:                  # a different type is a failure too
        return False, e
    return False, None


@contextlib.contextmanager
def _patched(name, value):
    orig = getattr(fem, name)
    setattr(fem, name, value)
    try:
        yield orig
    finally:
        setattr(fem, name, orig)


# ============================ 1. the corrector's catch ============================

#: A trial on the coarse Griffiths & Lane slope that runs past the first corrector
#: checkpoint (vp300), so the corrector is actually offered a state.
_F_CORR = 1.40


def _corrector_trial():
    return fem.solve_fem(_griffiths_coarse(), F=_F_CORR, max_iterations=1200,
                         max_iterations_ceiling=1200, fast_kernel=False)


def check_corrector_catch():
    print("\n1. the corrector's catch")

    def boom(exc):
        def _nr(*a, **k):
            raise exc
        return _nr

    with _patched("_solve_fem_newton", boom(NameError("name 'ks' is not defined"))):
        ok, e = _raises(_corrector_trial, NameError)
    check("a NameError inside the corrector propagates out of solve_fem", ok,
          f"{type(e).__name__ if e else 'nothing raised'}")

    with _patched("_solve_fem_newton", boom(KeyError("oob"))):
        ok, e = _raises(_corrector_trial, KeyError)
    check("so does a KeyError", ok, f"{type(e).__name__ if e else 'nothing raised'}")

    with _patched("_solve_fem_newton",
                  boom(np.linalg.LinAlgError("Singular matrix"))):
        sol = _quiet(_corrector_trial)
    atts = sol.get("corrector_attempts") or []
    check("a LinAlgError is a refusal, and the trial still gets its verdict",
          bool(atts) and all(not a.get("certified") for a in atts)
          and sol.get("verdict") is not None,
          f"{len(atts)} attempt(s), verdict {sol.get('verdict')}")
    check("  the refusal names LinAlgError",
          bool(atts) and all(a.get("refusal_type") == "LinAlgError"
                             and str(a.get("refusal", "")).startswith("LinAlgError")
                             for a in atts),
          "; ".join(f"{a.get('at')}: {a.get('refusal_type')}" for a in atts))

    with _patched("_solve_fem_newton", boom(RuntimeError("Factor is exactly singular"))):
        sol = _quiet(_corrector_trial)
    atts = sol.get("corrector_attempts") or []
    check("a singular factorization (RuntimeError) is a refusal naming it",
          bool(atts) and all(a.get("refusal_type") == "RuntimeError" for a in atts))

    # The hold test reads the same way.
    fd = _griffiths_coarse()
    state = {"u": np.zeros(1), "evp": []}
    with _patched("solve_fem", boom(AttributeError("'NoneType' has no attribute 'x'"))):
        ok, e = _raises(lambda: fem.corrector_hold_test(fd, 1.0, state), AttributeError)
    check("an AttributeError inside the hold test propagates", ok,
          f"{type(e).__name__ if e else 'nothing raised'}")
    with _patched("solve_fem", boom(FloatingPointError("overflow"))):
        out = _quiet(fem.corrector_hold_test, fd, 1.0, state)
    check("a FloatingPointError in the hold test is held=False, naming it",
          out.get("held") is False and out.get("error_type") == "FloatingPointError"
          and str(out.get("verdict", "")).startswith("FloatingPointError"),
          f"{out.get('held')}, {out.get('error_type')}")
    with _patched("solve_fem", boom(fem.FemInvariantError("seed from another model"))):
        ok, e = _raises(lambda: fem.corrector_hold_test(fd, 1.0, state),
                        fem.FemInvariantError)
    check("a FemInvariantError (a ValueError) is never read as a refusal", ok)


# ============================ 2. the other broad catches ==========================

def check_other_catches():
    print("\n2. the other broad catches")
    with tempfile.TemporaryDirectory() as d:
        stem = os.path.join(d, "model")
        with open(stem + "_fem_meta.json", "w") as f:
            f.write("{not json")
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            got = fem.import_fem_meta(stem)
        check("a corrupt meta sidecar is ignored, with a warning naming it",
              got is None and "could not be read" in buf.getvalue()
              and "model_fem_meta.json" in buf.getvalue(), buf.getvalue().strip())
        with open(stem + "_fem_meta.json", "w") as f:
            json.dump({"FS": 1.25}, f)
        check("a readable one is read", fem.import_fem_meta(stem) == {"FS": 1.25})

    check("an unrecognized unit system gives no unit",
          fem._ssrm_length_unit({"unit_system": "cubits"}) == ""
          and fem._ssrm_length_unit({"unit_system": "si"}) == "m")

    from scipy.sparse import csc_matrix
    K = csc_matrix(np.array([[1.0, 1.0], [1.0, 1.0]]))
    ok, e = _raises(lambda: fem._factorize_free_stiffness(K), RuntimeError)
    check("a singular stiffness still raises from the factorization", ok,
          f"{type(e).__name__ if e else 'nothing raised'}: {e}")

    fd = _griffiths_coarse()
    calls = [0]
    orig = fem.build_tri6_stiffness

    def bad(coords, E, nu):
        calls[0] += 1
        if calls[0] == 7:
            raise np.linalg.LinAlgError("negative Jacobian")
        return orig(coords, E, nu)

    with _patched("build_tri6_stiffness", bad):
        ok, e = _raises(lambda: fem.build_global_stiffness(
            fd["nodes"], fd["elements"], fd["element_types"],
            fd["element_materials"], fd["E_by_mat"], fd["nu_by_mat"], fd), ValueError)
    msg = str(e)
    check("an element whose stiffness cannot be built stops the assembly, naming "
          "the element and its material",
          ok and "element 6" in msg and "material" in msg
          and "negative Jacobian" in msg, msg)


# ======================= 3. no unreachable code, no unbound name ==================

def _unreachable(tree):
    """(function name, line) of every statement after a top-level return."""
    out = []
    for node in ast.walk(tree):
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            for i, st in enumerate(node.body[:-1]):
                if isinstance(st, (ast.Return, ast.Raise)):
                    out.append((node.name, node.body[i + 1].lineno))
                    break
    return out


def check_static():
    print("\n3. no unreachable code, no unbound name")
    with open(fem.__file__) as f:
        tree = ast.parse(f.read())
    dead = _unreachable(tree)
    check("no function body has statements after a return", not dead, f"{dead}")
    state = {"tn": np.zeros((2, 3)), "ts": np.zeros((2, 3)),
             "dt": np.zeros((2, 3)), "open": np.zeros((2, 3), dtype=bool)}
    ok, e = _raises(lambda: fem._nr_joint_slip([{"kind": "tie"}], state),
                    fem.FemInvariantError)
    check("_nr_joint_slip with a state and no joint group raises a plain error, "
          "not UnboundLocalError", ok, f"{type(e).__name__ if e else 'nothing'}")


# ============================ 4. the prepared model ===============================

def check_prepared():
    print("\n4. the prepared model and the continuation")
    fd = _griffiths_coarse()
    prep = _quiet(fem._prepare_fem_model, fd)
    kw = dict(F=1.0, max_iterations=20, max_iterations_ceiling=20, fast_kernel=False)
    sol = _quiet(fem.solve_fem, fd, _prepared=prep, **kw)
    check("a solve with the options the prepared model was built with runs",
          sol.get("iterations", 0) > 0)
    for opt, val in (("dt_scale", 0.5), ("tension_cutoff", True),
                     ("min_slip_depth", 2.0),
                     ("tension_cap_by_elem", np.full(len(fd["elements"]), 5.0)),
                     ("elastic_mask", np.ones(len(fd["elements"]), dtype=bool))):
        ok, e = _raises(lambda: fem.solve_fem(fd, _prepared=prep, **{opt: val}, **kw),
                        fem.FemInvariantError)
        check(f"a reuse with a different {opt} raises and names it",
              ok and opt in str(e), f"{e}")
    ok, e = _raises(lambda: fem.solve_fem(fd, _prepared=prep, min_slip_depth=0.0,
                                          **kw), fem.FemInvariantError)
    check("min_slip_depth = 0 is the same prepared model as none", not ok, f"{e}")

    prep2 = _quiet(fem._prepare_fem_model, fd)
    prep2["F_gravity"][3] += 1.0
    ok, e = _raises(lambda: fem.solve_fem(fd, _prepared=prep2, **kw),
                    fem.FemInvariantError)
    check("a reuse after a shared array changed raises", ok and "changed" in str(e),
          f"{e}")
    prep3 = _quiet(fem._prepare_fem_model, fd)
    prep3["gp_groups_static"][0]["B"][0, 0, 0] += 1.0
    ok, e = _raises(lambda: fem.solve_fem(fd, _prepared=prep3, **kw),
                    fem.FemInvariantError)
    check("  and so does one after a Gauss-point group's B changed", ok, f"{e}")

    # The continuation, on the stand-in solver (the bookkeeping, not the physics).
    common = dict(F_min=1.0, F_max=2.0, tolerance=0.01, failure_criterion="hybrid",
                  capture_failure_state=False)
    with _patched("solve_fem", _stand_in(lambda F: F <= 1.3, undecided_above=1.25)):
        r1 = _quiet(fem.solve_ssrm, fd, max_iterations=6000,
                    max_iterations_ceiling=6000, **common)
        check("the stand-in run is continuable",
              fem.ssrm_can_continue(r1) is not None, f"FS {r1.get('FS')}")
        edited = dict(fd)
        edited["c_by_mat"] = np.asarray(fd["c_by_mat"], dtype=float) * 1.1
        ok, e = _raises(lambda: fem.solve_ssrm(edited, resume=r1, max_iterations=20000),
                        ValueError)
        check("continuing on a model with an edited strength raises and names it",
              ok and "c_by_mat" in str(e) and "changed" in str(e), f"{e}")
        saved = fd["c_by_mat"]
        fd["c_by_mat"] = np.asarray(saved, dtype=float) * 0.9
        try:
            ok, e = _raises(lambda: fem.solve_ssrm(fd, resume=r1, max_iterations=20000),
                            ValueError)
        finally:
            fd["c_by_mat"] = saved
        check("  and so does continuing on the run's own model edited in place",
              ok and "c_by_mat" in str(e), f"{e}")
        rebuilt = {k: (np.array(v, copy=True) if isinstance(v, np.ndarray) else v)
                   for k, v in fd.items()}
        r2 = _quiet(fem.solve_ssrm, rebuilt, resume=r1, max_iterations=20000)
        check("continuing on an unchanged copy of the model runs",
              r2.get("resumed") is not None, f"FS {r2.get('FS')}")


def _stand_in(stands, undecided_above=None):
    """A solve_fem whose trial at F stands where ``stands(F)``; a trial that does
    not stand above ``undecided_above`` is undecided at the limit (continuable)."""
    def fake(fem_data, F=1.0, max_iterations=12000, max_iterations_ceiling=50000,
             _keep_resume=False, _resume_state=None, **_kw):
        limit = max(int(max_iterations), int(max_iterations_ceiling or 0))
        sol = dict(converged=False, stable=False, verdict="FAILED", u_ratio=3.0,
                   u_growth=0.5, max_displacement=0.01 * F, stop_reading=None,
                   creep_reading=None, plateau_iteration=None,
                   budget_extensions=0, iterations=50, exit_reason="diverging",
                   softened_1d_elements=np.zeros(0, dtype=bool),
                   _resume_state=None, F=F)
        if stands(F):
            sol.update(converged=True, stable=True, verdict="CONVERGED",
                       u_ratio=None, u_growth=None, iterations=200,
                       exit_reason="converged")
        elif undecided_above is not None and F > undecided_above and limit < 10000:
            sol.update(iterations=limit, verdict="AMBIGUOUS", u_growth=0.05,
                       exit_reason="inconclusive")
            if _keep_resume:
                sol["_resume_state"] = {"F": float(F), "iterations": limit,
                                        "iteration": limit - 1}
        return sol
    return fake


# ============================ 5. carried seeds ====================================

def check_seeds():
    print("\n5. carried seeds")
    fd = _griffiths_coarse()
    ok, e = _raises(lambda: fem.solve_fem(fd, F=1.0, max_iterations=10,
                                          fast_kernel=False,
                                          _softened_seed=np.ones(3, dtype=bool)),
                    fem.FemInvariantError)
    check("a post-peak set of the wrong length raises, naming both sizes",
          ok and "3 entries" in str(e) and "0 1D" in str(e), f"{e}")
    sol = _quiet(fem.solve_fem, fd, F=1.0, max_iterations=10, fast_kernel=False,
                 _softened_seed=np.zeros(0, dtype=bool))
    check("  an empty one is no seed", sol.get("iterations", 0) > 0)

    groups = [{"pairs": [(0, 0)] * 4}, {"pairs": [(1, 0)] * 2}]
    ok, e = _raises(lambda: fem._check_seed_state(np.zeros(9), [np.zeros((4, 4)),
                                                                np.zeros((2, 4))],
                                                  groups, 10, "the seed"),
                    fem.FemInvariantError)
    check("a Newton seed whose displacement is short raises, naming the sizes",
          ok and "(9,)" in str(e) and "(10,)" in str(e), f"{e}")
    ok, e = _raises(lambda: fem._check_seed_state(np.zeros(10), [np.zeros((4, 4))],
                                                  groups, 10, "the seed"),
                    fem.FemInvariantError)
    check("  and one with a missing Gauss-point group", ok and "1 group" in str(e),
          f"{e}")
    ok, e = _raises(lambda: fem._check_seed_state(np.zeros(10), [np.zeros((4, 4)),
                                                                 np.zeros((2, 4))],
                                                  groups, 10, "the seed"),
                    fem.FemInvariantError)
    check("  and a matching one passes", not ok and e is None, f"{e}")


# ============================ 6. monotonicity =====================================

def check_monotone():
    print("\n6. monotonicity")
    T = [{"F": 1.0, "stable": True, "exit_reason": "converged"},
         {"F": 1.5, "stable": False, "exit_reason": "diverging"},
         {"F": 1.25, "stable": True, "exit_reason": "converged"}]
    check("a monotone record reads as none", fem._ssrm_nonmonotone(T) is None)
    T2 = T + [{"F": 1.1, "stable": False, "exit_reason": "diverging"}]
    check("a standing trial above a failed one is the pair",
          fem._ssrm_nonmonotone(T2) == {"stood_F": 1.25, "failed_F": 1.1},
          f"{fem._ssrm_nonmonotone(T2)}")
    T3 = T + [{"F": 1.1, "stable": False, "exit_reason": "inconclusive"}]
    check("an undecided trial is not a failure", fem._ssrm_nonmonotone(T3) is None)

    fd = _griffiths_coarse()
    common = dict(F_min=1.0, F_max=2.0, tolerance=0.05, failure_criterion="hybrid",
                  capture_failure_state=False, max_iterations=6000,
                  max_iterations_ceiling=6000)
    # Stands below 0.9 and again between 1.9 and 2.1: F_min = 1 fails, F_max = 2
    # stands, both recorded (the pair is the highest standing F, 2.0, and the
    # lowest failed one, which the bisection below 1 then finds).
    buf = io.StringIO()
    with _patched("solve_fem", _stand_in(lambda F: F < 0.9 or 1.9 <= F < 2.1)):
        with contextlib.redirect_stdout(buf):
            r = fem.solve_ssrm(fd, **common)
    with _patched("solve_fem", _stand_in(lambda F: F < 0.9)):
        r0 = _quiet(fem.solve_ssrm, fd, **common)
    check("a run whose trials stood above a failure records the pair",
          (r.get("nonmonotone_trials") or {}).get("stood_F") == 2.0
          and (r.get("nonmonotone_trials") or {}).get("failed_F", 9.0) <= 1.0,
          f"{r.get('nonmonotone_trials')}")
    check("  and says so at debug level 0",
          "not monotone in F" in buf.getvalue())
    check("  and a monotone run carries no such key",
          "nonmonotone_trials" not in r0)


# ============================ 7. the in-situ state ================================

def check_init_state():
    print("\n7. the in-situ state")
    fd = _griffiths_coarse()
    real = fem._ssrm_displacement_limit

    def writing(*a, **k):
        st = k.get("_init_state")
        if st is not None:
            st["u"][0] += 1.0
        return real(*a, **k)

    with _patched("_ssrm_displacement_limit", writing):
        with _patched("solve_fem", _k0_stand_in()):
            ok, e = _raises(lambda: fem.solve_ssrm(
                fd, F_min=1.0, F_max=2.0, tolerance=0.05, k0=0.5,
                capture_failure_state=False, max_iterations=2000), fem.FemInvariantError)
    check("a trial that writes into the in-situ state is caught", ok, f"{e}")


def _k0_stand_in():
    """The stand-in, returning an in-situ state from the equilibration solve."""
    base = _stand_in(lambda F: F <= 1.3)

    def fake(fem_data, F=1.0, **kw):
        sol = base(fem_data, F=F, **kw)
        prep = kw.get("_prepared")
        sol.setdefault("plastic_elements", np.zeros(0, dtype=bool))
        sol.setdefault("unbalanced_force_ratio", 0.0)
        sol["_k0_state"] = {"u": np.zeros(prep["n_dof"]),
                            "evp": [np.zeros((len(g["pairs"]), 4))
                                    for g in prep["gp_groups_static"]]}
        return sol
    return fake


def main():
    print("=" * 72)
    print("FEM guard checks")
    print("=" * 72)
    check_corrector_catch()
    check_other_catches()
    check_static()
    check_prepared()
    check_seeds()
    check_monotone()
    check_init_state()
    print("\n" + "=" * 72)
    if FAILURES:
        print(f"FAILED ({len(FAILURES)}): " + ", ".join(FAILURES))
        return 1
    print("All FEM guard checks passed.")
    return 0


def run():
    """Failures as a list, for run_tests.py."""
    del FAILURES[:]
    main()
    return list(FAILURES)


if __name__ == "__main__":
    sys.exit(main())
