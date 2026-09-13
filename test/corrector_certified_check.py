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
preceding it, and three rows of the RS2 joint corpus turned on exactly this: their
standing edges were certified at 250 002, 250 006 and 250 003 sweeps against a
250 000 ceiling, and were being read as though nothing had ruled on them.

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


def main():
    print("=" * 72)
    print("Corrector certification checks")
    print("=" * 72)
    check_rule()
    check_readers_agree()
    check_record()
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
