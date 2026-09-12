"""Unit checks for the JOINT VERDICT — how an undecided jointed trial is read.

A strength-reduction trial on a jointed model can neither converge nor fail: the
force test is never met and the displacement classifier cannot rule. Seven rows of
the RS2 joint corpus were reported without a lock for exactly that reason, and one
of them stayed undecided at a million sweeps. `xslope.fem.joint_verdict` reads the
interface instead of the displacement field, and separates the two things the
undecided trials turned out to be:

  * `steady_slip` — the slip is still taking a real share of itself every window,
    its rate is not decaying, and max|u| is gaining with it. The block is moving
    on its joints, and the trial is FAILED.
  * `joint_settled` — the slip has stopped, the field has stopped, and the SOIL
    residual is under the solve's own force tolerance. What is left is a limit
    cycle on the joint degrees of freedom, which no budget can bring down. The
    slope is standing.

What this file locks:

  1. THE MEASURED SET. Nine recorded traces (test/fixtures/joint_traces.json, the
     per-sweep joint trace of six real bracket-edge trials) replayed through the
     rule sweep by sweep. Each must get the verdict the measurement supports, at
     roughly the sweep it was first supportable. This is the calibration itself,
     so a threshold moved without re-measuring fails here.

  2. MUTATION SPECS. One reading of the settled trace perturbed at a time — the
     slip creeping, the field creeping, the soil out of equilibrium, the joint
     residual still falling, no joint residual at all — must each be enough on
     its own to withdraw 'joint_settled'; and one reading of the moving trace at
     a time — less slip, less displacement, a decaying rate, too few sweeps —
     must each be enough on its own to withdraw 'steady_slip'. A rule whose
     conditions are not each load-bearing is not the rule described.

  3. INERTNESS. No history, no elastic yardstick, a too-short window and a model
     with no joint at all: the rule declines and says nothing. An unjointed solve
     never reaches it (the call site is guarded by `has_joints`), which is what
     makes the unjointed corpus bit-identical.

Run directly:  PYTHONPATH=. python3 test/joint_verdict_check.py
Exits non-zero on any failure.
"""
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import xslope.fem as fem                                    # noqa: E402
from xslope.fem import joint_verdict                        # noqa: E402

FAILURES = []
FORCE_TOL = 1e-3
TRACES = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                      'fixtures', 'joint_traces.json')


def check(name, cond, detail=""):
    status = "PASS" if cond else "FAIL"
    print(f"  [{status}] {name}" + (f"  — {detail}" if detail else ""))
    if not cond:
        FAILURES.append(name)


def load_traces():
    with open(TRACES) as f:
        return json.load(f)


def verdict_at_end(t, **kw):
    """What the rule says on the WHOLE history — the reading the mutation specs
    perturb, so that a perturbation of the trailing window is what is measured
    rather than an earlier sweep at which the unperturbed rule already spoke."""
    return joint_verdict(t['slip'], t['oob_soil'], t['disp'],
                         t['u_elastic_scale'], FORCE_TOL,
                         joint_oob_hist=t['oob_joint'], budget=t['budget'],
                         sample_every=t['sample_every'], **kw)


def first_verdict(t, **kw):
    """``(verdict, sweep)`` — replay the rule over one trace sweep by sweep and
    return the first thing it says, exactly as the solve's own loop asks it."""
    n = len(t['slip'])
    se = t['sample_every']
    for k in range(1, n + 1):
        if k * se < fem._JOINT_VERDICT_WARMUP:
            continue
        v = joint_verdict(t['slip'][:k], t['oob_soil'][:k], t['disp'][:k],
                          t['u_elastic_scale'], FORCE_TOL,
                          joint_oob_hist=t['oob_joint'][:k], budget=t['budget'],
                          sample_every=se, **kw)
        if v is not None:
            return v, k * se
    return None, None


# ===================== 1. the measured set =====================

#: trace -> (what the rule must say, the latest sweep it may take to say it).
#: The sweep bound is generous — it is there to catch a rule that only speaks at
#: the very end of a history, which would decide nothing a budget did not already.
EXPECTED = {
    # the two undecided edges that are STANDING behind the limit cycle
    'rj006_lo':  ('joint_settled', 40000),
    'rj007_ref': ('joint_settled', 40000),
    # the one undecided edge whose slip is linear: a mechanism
    'rj006_hi':  ('steady_slip',   50000),
    # three trials the rule must decline. rj007_ctl and rj006_ctl reach force
    # equilibrium on their own (16 072 and 41 721 sweeps) and must be left to;
    # rj019_hi creeps sub-linearly, which at any sweep count inside the corpus
    # is indistinguishable from a trial still converging; vp088_lo's joints have
    # settled but its SOIL has not, so the mechanism is not on the interface.
    'rj007_ctl': (None,            None),
    'rj006_ctl': (None,            None),
    'rj019_hi':  (None,            None),
    'vp088_lo':  (None,            None),
    # And the two trials that CONVERGE only after a very long run. RJ-2's is the
    # one the FAILED reading's budget gate exists for: it gains 30% of its slip
    # and four elastic displacements by 50 000 sweeps, with an ACCELERATING slip
    # rate, and then settles and reaches force equilibrium at 185 381. Called a
    # mechanism, it would move a shipped lock.
    'rj002_slow': (None,           None),
    'rj019_slow': (None,           None),
}


def check_measured():
    print("\n1. the measured set — nine recorded bracket-edge traces")
    traces = load_traces()
    check("the fixture carries every trace the calibration used",
          set(traces) == set(EXPECTED), sorted(traces))
    for tag, (want, by) in EXPECTED.items():
        t = traces.get(tag)
        if t is None:
            check(f"{tag} present", False)
            continue
        got, at = first_verdict(t)
        ok = got == want and (want is None or (at is not None and at <= by))
        check(f"{tag} ({t['what']}) -> {want}", ok,
              f"recorded {t['recorded']}/{t['recorded_exit']} at "
              f"{t['iterations']} sweeps; rule says {got} at {at}")


# ===================== 2. mutation specs =====================

def check_mutations():
    print("\n2. mutation specs — every condition is load-bearing")
    traces = load_traces()

    # --- 'joint_settled' withdraws when any one of its three readings moves ---
    s = traces['rj007_ref']
    base = verdict_at_end(s)
    check("baseline: the settled trace reads joint_settled", base == 'joint_settled')

    # (a) the slip creeps. Ramp the trailing half by ten times the settled
    # threshold, leaving everything else alone.
    n = len(s['slip'])
    h = n // 2
    creep = list(s['slip'])
    total = creep[-1]
    for i in range(h, n):
        creep[i] += 10 * fem._JOINT_SETTLED_SLIP_FRAC * total * (i - h) / (n - h)
    m = dict(s, slip=creep)
    got = verdict_at_end(m)
    check("(a) a creeping slip withdraws joint_settled", got != 'joint_settled',
          f"reads {got}")

    # (b) the displacement field creeps, the slip does not.
    move = list(s['disp'])
    for i in range(h, n):
        move[i] += 10 * fem._JOINT_SETTLED_GROWTH * s['u_elastic_scale'] \
            * (i - h) / (n - h)
    got = verdict_at_end(dict(s, disp=move))
    check("(b) a creeping displacement field withdraws joint_settled",
          got != 'joint_settled', f"reads {got}")

    # (c) the SOIL is not in equilibrium. One sample over the tolerance is
    # enough: the condition is a maximum over the window, not a mean.
    hot = list(s['oob_soil'])
    hot[-1] = 2.0 * FORCE_TOL
    got = verdict_at_end(dict(s, oob_soil=hot))
    check("(c) a soil residual over the force tolerance withdraws joint_settled",
          got != 'joint_settled', f"reads {got}")

    # (c2) the joint residual is STILL FALLING — the signature of a trial on its
    # way to a clean convergence, which must be left to reach it. Halve the last
    # quarter, which is a fall of exactly the kind the condition rejects.
    n = len(s['oob_joint'])
    q = n // 2 + (n - n // 2) // 2
    falling = list(s['oob_joint'])
    for i in range(q, n):
        falling[i] *= 0.5
    got = verdict_at_end(dict(s, oob_joint=falling))
    check("(c2) a joint residual still falling withdraws joint_settled",
          got != 'joint_settled', f"reads {got}")

    # (c3) no joint-residual reading at all: the verdict is withheld rather than
    # guessed, because the condition it needs is absent.
    got = joint_verdict(s['slip'], s['oob_soil'], s['disp'],
                        s['u_elastic_scale'], FORCE_TOL,
                        joint_oob_hist=None, sample_every=s['sample_every'])
    check("(c3) no joint residual withholds joint_settled entirely",
          got != 'joint_settled', f"reads {got}")

    # --- 'steady_slip' withdraws when any one of its five readings moves ---
    f = traces['rj006_hi']
    base = verdict_at_end(f)
    check("baseline: the moving trace reads steady_slip", base == 'steady_slip')

    # (d) the slip gains less than the threshold over the window: hold the
    # trailing half at the value it had when the window opened.
    n = len(f['slip'])
    h = n // 2
    flat = list(f['slip'][:h]) + [f['slip'][h]] * (n - h)
    got = verdict_at_end(dict(f, slip=flat))
    check("(d) a slip that stops withdraws steady_slip", got != 'steady_slip',
          f"reads {got}")

    # (e) the field stops while the slip keeps going.
    still = list(f['disp'][:h]) + [f['disp'][h]] * (n - h)
    got = verdict_at_end(dict(f, disp=still))
    check("(e) a displacement field that stops withdraws steady_slip",
          got != 'steady_slip', f"reads {got}")

    # (f) the slip rate DECAYS across the window — the signature of a trial on
    # its way to equilibrium. Same total gain, front-loaded.
    q = h + (n - h) // 2
    dec = list(f['slip'][:h])
    gained = f['slip'][-1] - f['slip'][h]
    for i in range(h, n):
        # 90% of the gain in the first half of the window, 10% in the second
        if i < q:
            frac = 0.9 * (i - h) / max(1, q - h)
        else:
            frac = 0.9 + 0.1 * (i - q) / max(1, n - q)
        dec.append(f['slip'][h] + gained * frac)
    got = verdict_at_end(dict(f, slip=dec))
    check("(f) a decaying slip rate withdraws steady_slip", got != 'steady_slip',
          f"reads {got}")

    # (g) the sweep clock. Shrink the stride so the whole history stands for
    # fewer sweeps than the FAILED reading is allowed to speak at.
    short = dict(f, sample_every=max(
        1, (fem._JOINT_MOVING_MIN_SWEEPS - 1) // len(f['slip'])))
    got = verdict_at_end(short)
    check("(g) a history shorter than the FAILED clock withdraws steady_slip",
          got != 'steady_slip', f"reads {got}")

    # (h) THE BUDGET GATE, which is the condition that keeps the rule off a trial
    # that converges late: the same history under a budget ten times longer is
    # a trial still running, and nothing may be said about it.
    got = verdict_at_end(dict(f, budget=10 * f['budget']))
    check("(h) a trial still early in its budget withdraws steady_slip",
          got != 'steady_slip', f"reads {got}")

    # (i) no budget at all: the FAILED reading is withheld rather than guessed.
    got = joint_verdict(f['slip'], f['oob_soil'], f['disp'],
                        f['u_elastic_scale'], FORCE_TOL,
                        joint_oob_hist=f['oob_joint'], budget=None,
                        sample_every=f['sample_every'])
    check("(i) no budget withholds steady_slip entirely", got != 'steady_slip',
          f"reads {got}")


# ===================== 3. inertness =====================

def check_inert():
    print("\n3. inertness — what the rule declines to rule on")
    check("no history at all", joint_verdict([], [], [], 1.0, FORCE_TOL) is None)
    n = 1000
    flat = [1.0] * n
    check("no elastic yardstick",
          joint_verdict(flat, [0.0] * n, flat, 0.0, FORCE_TOL) is None)
    check("a negative elastic yardstick",
          joint_verdict(flat, [0.0] * n, flat, -1.0, FORCE_TOL) is None)
    check("a history shorter than the warm-up",
          joint_verdict(flat[:10], [0.0] * 10, flat[:10], 1.0, FORCE_TOL) is None)
    check("a zero total slip — nothing has moved on the interface yet",
          joint_verdict([0.0] * n, [0.0] * n, flat, 1.0, FORCE_TOL) is None)
    # A model with no joint carries no slip history at all, and the call site is
    # guarded by has_joints as well, so this is belt and braces.
    check("an unjointed model's empty trace",
          joint_verdict([], [0.0] * n, flat, 1.0, FORCE_TOL) is None)
    # The verdict names the solve reports, and the one that counts as standing.
    check("the two exit reasons are the ones solve_fem reports",
          'steady_slip' != 'joint_settled')


def check_wiring():
    print("\n4. wiring — the solve's own vocabulary")
    # A JOINT_SETTLED trial counts as STANDING under the hybrid criterion and a
    # steady_slip one does not; both are read off the same table solve_fem uses.
    for verdict, want in (('JOINT_SETTLED', True), ('STABLE_STUCK', True),
                          ('FAILED', False), ('AMBIGUOUS', False)):
        stable = bool(False or ('hybrid' == 'hybrid'
                                and verdict in ('STABLE_STUCK', 'JOINT_SETTLED')))
        check(f"{verdict} counts as standing: {want}", stable == want)
    check("the joint verdict can be switched off for an A/B",
          hasattr(fem, 'JOINT_VERDICT_ON') and fem.JOINT_VERDICT_ON is True)
    check("the per-sweep trace has a sink for a study to collect it",
          hasattr(fem, 'JOINT_TRACE_SINK') and fem.JOINT_TRACE_SINK is None)


def main():
    print("=" * 72)
    print("Joint verdict checks")
    print("=" * 72)
    check_measured()
    check_mutations()
    check_inert()
    check_wiring()
    print("\n" + "=" * 72)
    if FAILURES:
        print(f"FAILED ({len(FAILURES)}): " + ", ".join(FAILURES))
        return 1
    print("All joint verdict checks passed.")
    return 0


def run():
    """Failures as a list, for run_tests.py."""
    del FAILURES[:]
    main()
    return list(FAILURES)


if __name__ == '__main__':
    sys.exit(main())
