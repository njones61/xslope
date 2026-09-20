"""Unit checks for what makes a bracket edge ANSWERED: the corrector certification.

A strength-reduction bracket halves an interval by asking of each trial factor
whether the model stands there. A trial that reaches its sweep budget with nothing
to say has not answered, and the two tools that read a lock's trial record treat
it as undecided — which is right, and is why `tools/ssrm_trial_audit.py` exists.

But a trial can reach its budget and still be answered. On a jointed model the
Newton corrector is offered the viscoplastic loop's own state at its checkpoints,
including the one at the budget exit; where it reaches equilibrium inside the
force, yield and displacement gates, `solve_ssrm` records a `corrector` block on
that trial. That block says an INDEPENDENT driver found a static equilibrium from
this trial's state. It answers the standing question whatever the sweep count
preceding it, and rows of the corpus turn on exactly this: RS2-52's standing edge
is certified after 258 605 sweeps and RS2-53's at 500 000 — both past the 250 000
sweeps their tags budget — and RJ-4's after 186 859. Read on the sweep count alone
they would be trials nothing had ruled on.

The same corrector is offered the K0 IN-SITU EQUILIBRATION, which is a solve at
full strength before any trial and on a jointed model is where a row's sweeps go.
Its ladder is continued past the trials' three rungs (`_K0_CORRECTOR_CHECKPOINTS`),
because the equilibration's product is a state and not a verdict: certifying it
earlier changes what the step costs, not what a bracket reads.

What this file locks:

  1. THE RULE. A certified trial is decided at, above and below the ceiling. An
     uncertified one is decided only below it. The rule is unchanged for every
     other reason a trial can fail to answer — the undecided verdicts, and the
     `inconclusive` exit.

  2. THE THREE READERS AGREE. `tools/ssrm_trial_audit.trial_decided`,
     `tools/lock_edges._decided` and `run_tests._edge_reading` must give the same
     answer for the same trial. They are what decide, respectively, whether a
     lock is budget-bound, whether its edge pair may be written, and whether the
     suite's two-trial check passes — and a disagreement between them means a
     pair is written that the audit then reads as unsafe.

  3. WHAT TRAVELS INTO THE RECORD. `fem.ssrm_run_record` must carry the
     certification into a committed sidecar (or the readers above cannot see it)
     and must NOT carry the machine time with it (or a sidecar changes when the
     machine does).

  4. THE K0 STEP. Its ladder is the trials' ladder continued, it certifies a
     jointed model's in-situ state and says where, and on an unjointed model it
     is inert — the same equilibration, reached at the same sweep, with the
     switch on and off.

  5. A REFUSAL'S RESIDUAL. A failed one-shot corrector publishes the residual of
     the state it actually returned. It must agree with `nr_diag['oob']`, not the
     zero that `last_oob` was initialized to before the attempt.
"""
import os
import sys

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _ROOT)
sys.path.insert(0, os.path.join(_ROOT, "tools"))

FAILURES = []


def check(label, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {label}" + (f"  — {detail}" if detail else ""))
    if not ok:
        FAILURES.append(label)


#: The jointed model the K0 leg reads: a slab on a cohesionless plane at 20 degrees
#: whose friction angle is two degrees above it. The in-situ field has a slope to
#: redistribute against and the interface is close enough to its limit that the
#: sweep does not settle inside the first checkpoint, so the corrector is actually
#: asked for something — it certifies at `vp300`, where the plain sweep would have
#: run to 311.
K0_SLAB = dict(beta=20.0, HV=3.0, phi_j=22.0, cj=0.0, ts=3.0, s1d=3.0,
               E_void=100.0)

CERT = {"driver_of_record": "corrector", "checkpoint": "vp_exit",
        "vp_iterations": 250000, "nr_iterations": 6, "nr_force_evals": 26,
        "oob": 3.1e-11, "force_tol": 1e-3, "yield_violation": 0.0,
        "yield_tol": 1e-6}


def trial(F=1.0, verdict="CONVERGED", iterations=100, cert=False, exit_reason=None):
    t = {"F": F, "verdict": verdict, "iterations": iterations}
    if exit_reason is not None:
        t["exit_reason"] = exit_reason
    if cert:
        t["corrector"] = dict(CERT)
    return t


def check_rule():
    print("\n1. The rule")
    import ssrm_trial_audit as audit

    C = 250000
    check("a plain trial under the ceiling is decided",
          audit.trial_decided(trial(iterations=1221), C))
    check("a plain trial AT the ceiling is not",
          not audit.trial_decided(trial(iterations=C), C))
    check("a plain trial past the ceiling is not",
          not audit.trial_decided(trial(iterations=C + 6), C))
    check("a CERTIFIED trial past the ceiling IS decided",
          audit.trial_decided(trial(iterations=C + 6, cert=True), C),
          "the corpus case: 250 006 sweeps, certified at the budget exit")
    check("a certified trial under the ceiling is decided too",
          audit.trial_decided(trial(iterations=350, cert=True), C))
    check("a certified FAILED trial past the ceiling is decided",
          audit.trial_decided(trial(verdict="FAILED", iterations=C + 1, cert=True), C))

    # The certification lifts the CEILING test and nothing else.
    check("an undecided verdict stays undecided even certified",
          not audit.trial_decided(
              trial(verdict="STABLE_STUCK", iterations=10, cert=True), C),
          "STABLE_STUCK is not a verdict, certified or not")
    check("an inconclusive exit stays undecided even certified",
          not audit.trial_decided(
              trial(iterations=10, cert=True, exit_reason="inconclusive"), C))
    check("an empty corrector block is not a certification",
          not audit.trial_decided(
              {"F": 1.0, "verdict": "CONVERGED", "iterations": C, "corrector": {}}, C))
    check("a non-dict corrector value is not a certification",
          not audit.trial_decided(
              {"F": 1.0, "verdict": "CONVERGED", "iterations": C,
               "corrector": True}, C))
    check("JOINT_SETTLED still decides",
          audit.trial_decided(trial(verdict="JOINT_SETTLED", iterations=99), C))


def check_standing_verdicts():
    """A JOINT_SETTLED trial STANDS; it is not a refusal.

    The bracket a record ended on is read as the highest standing trial and the
    lowest refused one. Reading JOINT_SETTLED as a refusal puts the refused edge at
    the BOTTOM of the bracket on any row that has one, and the row is then reported
    as checking a pair its own record did not end on — which is what happened to
    RJ-4, whose F = 0.5 trial settles on its joints."""
    print("\n1b. Which verdicts stand")
    import ssrm_trial_audit as audit

    import json
    import tempfile
    rec = {"trials": [
        {"F": 0.5, "verdict": "JOINT_SETTLED", "iterations": 5001},
        {"F": 1.125, "verdict": "CONVERGED", "iterations": 1024},
        {"F": 1.30078125, "verdict": "CONVERGED", "iterations": 186870},
        {"F": 1.3203125, "verdict": "FAILED", "iterations": 75061},
        {"F": 3.0, "verdict": "FAILED", "iterations": 181},
    ]}
    fh = tempfile.NamedTemporaryFile("w", suffix=".json", delete=False)
    json.dump(rec, fh)
    fh.close()
    got = audit.audit_one({"max_iter": "250000"}, fh.name)
    os.unlink(fh.name)
    lo, hi = (got or {}).get("edges", (None, None))
    check("the standing edge is the highest trial that stood",
          lo == 1.30078125, f"got {lo}")
    check("the refused edge is the lowest trial that failed",
          hi == 1.3203125, f"got {hi}")
    check("every trial is counted decided", (got or {}).get("decided") == 5)


def check_readers_agree():
    print("\n2. The three readers agree")
    import ssrm_trial_audit as audit
    import lock_edges
    import run_tests

    C = 250000
    cases = [
        ("plain, under", trial(iterations=1221)),
        ("plain, at", trial(iterations=C)),
        ("plain, past", trial(iterations=C + 6)),
        ("certified, past", trial(iterations=C + 6, cert=True)),
        ("certified, under", trial(iterations=350, cert=True)),
        ("stuck, certified", trial(verdict="STABLE_STUCK", iterations=10, cert=True)),
        ("inconclusive, certified",
         trial(iterations=10, cert=True, exit_reason="inconclusive")),
    ]
    for label, t in cases:
        a = bool(audit.trial_decided(t, C))
        b = bool(lock_edges._decided(t, C))
        _, c = run_tests._edge_reading([t], t["F"], C)
        check(f"{label}: audit / lock_edges / suite agree",
              a == b == bool(c), f"{a} / {b} / {bool(c)}")


def check_record():
    print("\n3. What travels into a committed record")
    from xslope import fem

    full = dict(CERT)
    full["wall"] = 9.61
    full["attempts"] = [{"at": "vp300", "certified": False, "wall": 0.4}]
    result = {"trials": [{"F": 1.0, "verdict": "CONVERGED", "iterations": 250006,
                          "corrector": full, "corrector_attempts": [1, 2],
                          "wall": 12.0}]}
    rec = fem.ssrm_run_record(result)
    got = (rec.get("trials") or [{}])[0]
    cert = got.get("corrector")
    check("the certification reaches the record", isinstance(cert, dict) and bool(cert))
    check("its out-of-balance ratio travels", cert.get("oob") == CERT["oob"])
    check("its driver travels", cert.get("driver_of_record") == "corrector")
    check("the wall time does NOT travel", "wall" not in (cert or {}))
    check("the attempts log does NOT travel", "attempts" not in (cert or {}))
    check("the trial's own wall time does not travel", "wall" not in got)
    check("corrector_attempts does not travel", "corrector_attempts" not in got)

    # A refused or never-offered trial carries no certification at all, so the
    # key's PRESENCE is what the readers may rely on.
    plain = {"trials": [{"F": 1.0, "verdict": "FAILED", "iterations": 2061}]}
    got2 = (fem.ssrm_run_record(plain).get("trials") or [{}])[0]
    check("a trial with no corrector carries no corrector key",
          "corrector" not in got2)
    empty = {"trials": [{"F": 1.0, "verdict": "FAILED", "iterations": 2061,
                         "corrector": {}}]}
    got3 = (fem.ssrm_run_record(empty).get("trials") or [{}])[0]
    check("an empty corrector block is dropped rather than written",
          "corrector" not in got3)


def _k0_equilibration(fem_data, ladder_on=True):
    """One SSRM's in-situ equilibration: its record, and the ladder the step was
    handed.

    `trial_factors` keeps it to a single trial after the equilibration, which is
    the cheapest way to reach the step through the path that owns the switch.
    """
    import contextlib
    import io
    from xslope import fem

    seen = {"n": 0, "rungs": "not called"}
    orig = fem.solve_fem

    def wrapped(*a, **k):
        seen["n"] += 1
        if seen["n"] == 1:
            seen["rungs"] = k.get("_corrector_rungs")
        return orig(*a, **k)

    was = fem.K0_CORRECTOR_LADDER_ON
    fem.K0_CORRECTOR_LADDER_ON = ladder_on
    fem.solve_fem = wrapped
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            res = fem.solve_ssrm(fem_data, F_min=1.0, F_max=1.5, k0=1.0,
                                 trial_factors=[1.0], max_iterations=6000,
                                 max_iterations_ceiling=6000,
                                 capture_failure_state=False, debug_level=0)
    finally:
        fem.K0_CORRECTOR_LADDER_ON = was
        fem.solve_fem = orig
    return res.get("k0_equilibration") or {}, seen["rungs"]


def check_k0_step():
    """The equilibration is a solve like any other, and the corrector is offered
    it like any other — with the ladder continued, because the step's product is
    a state rather than a verdict (see fem._K0_CORRECTOR_CHECKPOINTS)."""
    print("\n4. The K0 in-situ equilibration")
    import warnings

    from xslope import fem

    n = len(fem._CORRECTOR_CHECKPOINTS)
    check("the K0 ladder is the trials' ladder, continued",
          fem._K0_CORRECTOR_CHECKPOINTS[:n] == fem._CORRECTOR_CHECKPOINTS
          and len(fem._K0_CORRECTOR_CHECKPOINTS) > n,
          f"{fem._K0_CORRECTOR_CHECKPOINTS}")
    check("the continued rungs are the ones a creeping equilibration reaches",
          min(fem._K0_CORRECTOR_CHECKPOINTS[n:]) >= 10000)

    # A JOINTED model: a slab on a plane at 20 degrees, where the K0 field has a
    # slope to redistribute against and the sweep does not settle at once.
    sys.path.insert(0, os.path.join(_ROOT, "test"))
    from joint_element_check import slab_model

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        _d, jointed, _g = slab_model(**K0_SLAB)
        eq, rungs = _k0_equilibration(jointed)
    check("the jointed step is handed the continued ladder",
          tuple(rungs or ()) == fem._K0_CORRECTOR_CHECKPOINTS, f"{rungs}")
    check("a jointed in-situ state is established", bool(eq.get("stable")),
          f"{eq.get('verdict')} in {eq.get('iterations')} sweep(s)")
    check("the corrector certified it, and the record says where",
          str(eq.get("certified_at") or "").startswith(("vp", "rule:", "gate:")),
          f"at {eq.get('certified_at')}, {eq.get('nr_iterations')} Newton "
          f"iteration(s)")
    _, off_rungs = _k0_equilibration(jointed, ladder_on=False)
    check("and with the switch off it runs the trials' three rungs",
          off_rungs is None, f"{off_rungs}")

    # An UNJOINTED model is never handed the continued ladder, so the switch
    # cannot reach it: the same equilibration, at the same sweep, either way.
    from xslope.fem import build_fem_data
    from xslope.fileio import load_slope_data
    from xslope.mesh import build_mesh_from_polygons, get_material_polygons

    base = os.path.join(_ROOT, 'docs', 'fem', 'files', 'xslope_griffiths1.xlsx')
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sd = load_slope_data(base)
        mesh = build_mesh_from_polygons(get_material_polygons(sd),
                                        target_size=2.0, element_type='tri6')
        plain = build_fem_data(sd, mesh)
        on, on_rungs = _k0_equilibration(plain, ladder_on=True)
        off, _ = _k0_equilibration(plain, ladder_on=False)
    check("an unjointed step is never handed the continued ladder",
          on_rungs is None, f"{on_rungs}")
    keys = ("iterations", "converged", "verdict", "certified_at", "n_plastic",
            "max_displacement", "unbalanced_force_ratio")
    check("so its equilibration is the same with the switch on and off",
          all(on.get(k) == off.get(k) for k in keys)
          and bool(on.get("converged")),
          f"{on.get('iterations')} sweep(s), max|u| = "
          f"{on.get('max_displacement', float('nan')):.6g}")


def check_refusal_residual():
    """A refused seeded increment reports that increment's nonzero residual."""
    print("\n5. A refused corrector's residual")
    import contextlib
    import io
    import math

    import run_tests as rt
    from xslope import fem

    page = os.path.join(_ROOT, "docs", "tutorials",
                        "fem03_block_wall_joints.md")
    tag = next(t for t in rt.parse_test_tags(page)
               if t.get("benchmark") == "FEM-3-sheet-jointed-ssrm")
    tag = dict(tag)
    tag["file"] = os.path.normpath(
        os.path.join(os.path.dirname(page), tag["file"]))
    with contextlib.redirect_stdout(io.StringIO()):
        fem_data, kwargs, _fmin, _fmax, tol = rt.build_fem_ssrm_case(tag)
    kwargs = dict(kwargs)
    kwargs.update(max_iterations=300, max_iterations_ceiling=300,
                  capture_failure_state=False)
    with contextlib.redirect_stdout(io.StringIO()):
        with rt._force_fast_kernel(fem, False):
            result = fem.solve_ssrm(
                fem_data, F_min=1.0, F_max=2.0, tolerance=tol,
                debug_level=1, trial_factors=[1.5546875], **kwargs)
    attempt = next(a for a in result["trials"][0]["corrector_attempts"]
                   if not a.get("certified"))
    public = attempt.get("oob")
    diagnostic = (attempt.get("nr_diag") or {}).get("oob")
    check("the fixture reaches a refused corrector",
          attempt.get("exit_reason") == "diverging",
          f"{attempt.get('at')}: {attempt.get('exit_reason')}")
    check("the refusal publishes a finite, nonzero out-of-balance",
          public is not None and math.isfinite(public) and public > 0.0,
          f"{public}")
    check("the public and inner diagnostic readings are the same attempt",
          public == diagnostic, f"{public} / {diagnostic}")


def main():
    print("=" * 72)
    print("Corrector certification checks")
    print("=" * 72)
    check_rule()
    check_standing_verdicts()
    check_readers_agree()
    check_record()
    check_k0_step()
    check_refusal_residual()
    print("\n" + "=" * 72)
    if FAILURES:
        print(f"FAILED ({len(FAILURES)}): " + ", ".join(FAILURES))
        return 1
    print("All corrector certification checks passed.")
    return 0


def run():
    """Failures as a list, for run_tests.py."""
    del FAILURES[:]
    main()
    return list(FAILURES)


if __name__ == '__main__':
    sys.exit(main())
