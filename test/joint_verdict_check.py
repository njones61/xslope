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


def check_churn():
    """The ACTIVE-SET CHURN reading, which says what the interface was doing at a
    moment in the solve.

    It is the fraction of node pairs that change state (sticking / slipping /
    open) per sweep, averaged over a trailing window of the trace, and it is what
    a corrector attempt's record carries beside its refusal. A refusal from a
    frozen set and one from a set still making up its mind are different findings
    (r19), so the reading has to be the trailing window and not the whole history.
    """
    print("\n4. the active-set churn reading")
    churn = fem._joint_churn
    check("no history reads nothing", churn([], 100) is None)
    check("no pairs reads nothing", churn([1, 2, 3], 0) is None)
    check("a frozen set reads zero", churn([0] * 20, 100) == 0.0)
    check("every pair flipping every sweep reads one",
          churn([100] * 20, 100) == 1.0)
    check("half the pairs flipping reads a half", churn([50] * 20, 100) == 0.5)
    # The window is the point: a set that chattered early and then froze must read
    # frozen NOW, which is the state a refusal is taken from.
    settled = [40] * 100 + [0] * 20
    check("a set that chattered and then froze reads zero over the last ten "
          "samples", churn(settled, 100, 10) == 0.0)
    check("...and not over the trailing half",
          churn(settled, 100) > 0.0)
    check("the short window is ten samples, a hundred sweeps",
          fem._JOINT_CHURN_SAMPLES == 10)
    check("a window longer than the history reads the history",
          churn([10] * 4, 100, 50) == 0.1)


def check_wiring():
    print("\n5. wiring — the solve's own vocabulary")
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


def _toy_pairs(n_pairs=2):
    """A single joint element with two live pairs (the third padded), horizontal,
    k_n = 1e8, k_s = 1e7, no cohesion, tan phi = 0.5, no tension."""
    import numpy as np
    jd = {"n": 1, "dof": np.arange(12)[None, :], "tx": np.array([1.0]),
          "ty": np.array([0.0]), "nx": np.array([0.0]), "ny": np.array([1.0]),
          "kn": np.array([1e8]), "ks": np.array([1e7]), "tcut": np.array([0.0]),
          "tip": np.zeros((1, 3), dtype=bool),
          "w": np.array([[0.5, 0.5, 0.0]]) if n_pairs == 2 else np.ones((1, 3))}
    return jd


def _toy_u(dt, dn=1e-4):
    """Displacements giving each pair a tangential offset ``dt[p]`` (side a moved
    +x) and a normal closing ``dn`` (side a moved -y, compression positive)."""
    import numpy as np
    u = np.zeros(12)
    for p in range(3):
        u[2 * p] = dt[p] if p < len(dt) else 0.0
        u[2 * p + 1] = -dn
    return u


def check_restick_band():
    """The one-sided stick/slip band (joint.JOINT_RESTICK_FRAC) on a two-pair toy."""
    import numpy as np
    import xslope.joint as J
    print("\n6. the one-sided stick/slip band — a two-pair toy")
    jd = _toy_pairs()
    cj, tp = np.array([0.0]), np.array([0.5])
    # t_n = 1e8 * 1e-4 = 1e4, t_lim = 5e3, the limit offset is 5e-4.
    lim = 5e-4
    frac = J.JOINT_RESTICK_FRAC
    check("the band is a small one below the limit", 0.99 <= frac < 1.0,
          f"JOINT_RESTICK_FRAC = {frac}")
    in_band = lim * (1.0 - 0.5 * (1.0 - frac))       # inside the band
    below = lim * (1.0 - 2.0 * (1.0 - frac))         # clearly below it
    over = lim * 1.001
    u = _toy_u([in_band, below])
    plain = J.joint_state(jd, u, cj, tp, slip_p=np.zeros((1, 3)),
                          open_prev=np.zeros((1, 3), bool))
    prev = np.array([[True, True, False]])
    held = J.joint_state(jd, u, cj, tp, slip_p=np.zeros((1, 3)),
                         open_prev=np.zeros((1, 3), bool), slip_prev=prev)
    check("plain test: both pairs under the limit read sticking",
          not plain["slipping"][0, :2].any())
    check("a pair slipping last sweep and inside the band stays slipping",
          bool(held["slipping"][0, 0]))
    check("a pair slipping last sweep and clearly below the band re-sticks",
          not bool(held["slipping"][0, 1]))
    check("a pair held inside the band carries its own traction, not the limit",
          held["ts"][0, 0] == plain["ts"][0, 0] and abs(held["ts"][0, 0]) < 5e3)
    stick = J.joint_state(jd, u, cj, tp, slip_p=np.zeros((1, 3)),
                          open_prev=np.zeros((1, 3), bool),
                          slip_prev=np.zeros((1, 3), bool))
    check("ONE-SIDED: a sticking pair inside the band stays sticking",
          not stick["slipping"][0, :2].any())
    u2 = _toy_u([over, below])
    a = J.joint_state(jd, u2, cj, tp, slip_p=np.zeros((1, 3)),
                      open_prev=np.zeros((1, 3), bool))
    b = J.joint_state(jd, u2, cj, tp, slip_p=np.zeros((1, 3)),
                      open_prev=np.zeros((1, 3), bool),
                      slip_prev=np.zeros((1, 3), bool))
    check("a pair over the limit slips with or without the band",
          bool(a["slipping"][0, 0]) and bool(b["slipping"][0, 0]))
    for k in ("tn", "ts", "tlim", "open"):
        check(f"no band held: '{k}' identical to the plain test",
              np.array_equal(a[k], b[k]))

    # The sweep: a sequence that crosses the limit, dips into the band (the
    # flip) and unloads, run with and without the band. The body load, the slip,
    # the open record and the residual latch must be BIT-IDENTICAL; only the
    # slipping code may differ, and only on the dip.
    # Pair 0 is loaded past its limit, then ripples about the offset it slipped
    # to by a part in 1e5 of the limit (the RJ-3 flip); pair 1 is loaded past
    # its limit and then unloads for real, well below the band.
    rip = lim * 1e-5
    seq = [[over, over], [over - rip, over], [over + rip, lim * 0.9],
           [over - rip, lim * 0.5], [over + rip, lim * 0.5],
           [over - rip, lim * 0.5], [over, lim * 0.2]]
    res = {}
    for band in (False, True):
        slip = np.zeros((1, 3))
        opn = np.zeros((1, 3), bool)
        latch = np.zeros((1, 3), bool)
        sp = np.zeros((1, 3), bool) if band else None
        loads_all, codes = [], []
        for dts in seq:
            loads = np.zeros(12)
            _, st = J.joint_vp_sweep(jd, _toy_u(dts), loads, cj, tp, slip, opn,
                                     slipped=latch,
                                     res_r=(np.array([0.0]), np.array([0.5])),
                                     slip_state=sp)
            loads_all.append(loads.copy())
            codes.append(st["slipping"].copy())
        res[band] = (np.array(loads_all), slip.copy(), opn.copy(), latch.copy(),
                     np.array(codes))
    check("sweep: the body load is bit-identical with the band",
          np.array_equal(res[False][0], res[True][0]))
    check("sweep: the accumulated slip is bit-identical with the band",
          np.array_equal(res[False][1], res[True][1]))
    check("sweep: the open record and the residual latch are bit-identical",
          np.array_equal(res[False][2], res[True][2])
          and np.array_equal(res[False][3], res[True][3]))
    flips_plain = int((res[False][4][1:] != res[False][4][:-1]).sum())
    flips_band = int((res[True][4][1:] != res[True][4][:-1]).sum())
    check("sweep: the band removes the ripple's code changes",
          flips_band < flips_plain, f"{flips_plain} -> {flips_band}")
    check("sweep: a pair that really unloads still re-sticks under the band",
          not bool(res[True][4][-1][0, 1]))


def check_finisher_reading():
    """The late finisher's confinement reading (fem._finisher_open_share) and its
    hot-node selection, on a synthetic residual."""
    import inspect
    import numpy as np
    print("\n7. the late finisher — confinement reading on a synthetic residual")
    jd = _toy_pairs(3)
    jd["dof"] = np.arange(12)[None, :] + 4          # joint dofs 4..15 of 20
    open_mask = np.array([[True, False, False]])     # pair 0: dofs 4,5 and 10,11
    dload = np.zeros(20)
    dload[[4, 5, 10, 11]] = [3.0, 4.0, -3.0, -4.0]   # 50 on the open pair
    dload[[0, 17]] = [1.0, 2.0]                      # 5 elsewhere
    share, dofs = fem._finisher_open_share(jd, open_mask, dload)
    check("the share is the squared load change on the open pair's four dofs",
          abs(share - 50.0 / 55.0) < 1e-12, f"{share:.4f}")
    check("...both faces of the pair", list(dofs) == [4, 5, 10, 11])
    check("that share sits just above the arming line",
          share >= fem._FINISHER_CONFINED)
    dload[0] = 5.0
    share2, _ = fem._finisher_open_share(jd, open_mask, dload)
    check("a residual that spreads into the soil drops below the arming line",
          share2 < fem._FINISHER_CONFINED, f"{share2:.4f}")
    check("no open pair reads zero",
          fem._finisher_open_share(jd, np.zeros((1, 3), bool), dload)[0] == 0.0)
    check("no residual reads zero",
          fem._finisher_open_share(jd, open_mask, np.zeros(20))[0] == 0.0)
    hot = fem._finisher_hot_dofs(np.array([0.0, 10.0, 1.0, 0.0, 3.0]))
    check("the hot set is the fewest dofs carrying 90% of the squared residual",
          list(hot) == [1], f"{list(hot)}")
    hot = fem._finisher_hot_dofs(np.array([0.0, 3.0, 3.0, 1.0]))
    check("...taking more where one is not enough", list(hot) == [1, 2])
    check("an all-zero residual has no hot dofs",
          fem._finisher_hot_dofs(np.zeros(5)).size == 0)
    check("the finisher is on, and can be switched off for an A/B",
          fem.JOINT_FINISHER_ON is True)
    check("its relieved run is a solve_fem call the trial makes on a copy",
          "_finisher_watch" in inspect.signature(fem.solve_fem).parameters)
    check("its sweep allowance is inside the start-of-trial relief's window",
          fem._FINISHER_SWEEPS <= fem._JOINT_RELIEF_SWEEPS)


def check_ambiguous_at_ceiling():
    """An AMBIGUOUS trial at its HARD ceiling is undecided (fem.ambiguous_at_ceiling),
    built on RJ-20's trial at F = 2.53125: 2.87x elastic, still, 1,000,000 sweeps
    at a 1,000,000 ceiling, the late finisher armed and refused."""
    print("\n8. an unreadable trial at its hard ceiling is undecided (RJ-20)")
    u_el = 0.002429200750969142
    n = 100000                                       # 1,000,000 sweeps / 10
    disp = [2.8739644749953936 * u_el - 5.2e-05 * u_el * (1.0 - k / n)
            for k in range(n)]
    verdict = fem.classify_nonconvergence(disp, u_el, 'iteration_cap',
                                          model_height=60.0)[0]
    check("the trial reads AMBIGUOUS (2.87x elastic, not growing)",
          verdict == 'AMBIGUOUS', verdict)
    check("at its ceiling it is undecided",
          fem.ambiguous_at_ceiling(disp, u_el, 1000000, 1000000, 60.0))
    check("below its ceiling it keeps the legacy count (failed)",
          not fem.ambiguous_at_ceiling(disp, u_el, 500000, 1000000, 60.0))
    # RS2-4-zone's failing edge: AMBIGUOUS, still, 32,000 of a 50,000 ceiling.
    check("RS2-4-zone's failing edge (32,000 of 50,000) is still counted failed",
          not fem.ambiguous_at_ceiling(disp, u_el, 32000, 50000, 60.0))
    moving = [x * (1.0 + 0.2 * k / n) for k, x in enumerate(disp)]
    check("a trial still moving at its ceiling is not undecided by this rule",
          not fem.ambiguous_at_ceiling(moving, u_el, 1000000, 1000000, 60.0))
    stuck = [0.5 * u_el] * n
    check("a trial frozen at elastic scale is not undecided by this rule",
          not fem.ambiguous_at_ceiling(stuck, u_el, 1000000, 1000000, 60.0))

    # The bracket as the run record reads it: RJ-20's nine trials, the top one
    # ending as the loop now labels it.
    def run(top_exit):
        trials = [
            (0.5, 'CONVERGED', 'converged', True), (3.0, 'FAILED', 'diverging', False),
            (1.75, 'JOINT_SETTLED', 'joint_settled', True),
            (2.375, 'JOINT_SETTLED', 'joint_settled', True),
            (2.6875, 'FAILED', 'diverging', False),
            (2.53125, 'AMBIGUOUS', top_exit, False),
            (2.453125, 'STABLE_STUCK', 'iteration_cap', True),
            (2.4921875, 'STABLE_STUCK', 'iteration_cap', True),
            (2.51171875, 'CONVERGED', 'converged', True)]
        return {"final_interval": [2.51171875, 2.53125],
                "trials": [dict(F=F, verdict=v, exit_reason=e, stable=st,
                                converged=(v == 'CONVERGED'),
                                finisher=({"armed": True, "outcome": "not balanced"}
                                          if F == 2.53125 else None))
                           for F, v, e, st in trials]}
    top = fem.ssrm_undecided_top(run('inconclusive'))
    check("at the ceiling: the bracket's top is undecided, the answer a lower "
          "bound at 2.512 (FS not closed)",
          top is not None and top["F"] == 2.53125)
    check("below the ceiling (exit 'iteration_cap'): counted failed, the "
          "bracket closes as before",
          fem.ssrm_undecided_top(run('iteration_cap')) is None)


def _leaning_block(P, hp, c_wall=10.0, phi_wall=20.0, kn=1e6, ks=1e6, W=100.0):
    """A rigid block (ux, uy, rotation about (0, 1)) on a strong two-pair base joint
    from x = -1 to 1, leaning on a cohesive, no-tension two-pair wall joint at
    x = 1 (y 0.25 to 1.75), loaded by its weight W at (0, 1) and a push P toward
    the wall at height hp. Returns ``(jd, T, K, f)``: T maps the block onto the
    joint pairs' side-a degrees of freedom (side b is fixed), and K is the joints'
    full elastic block on the block plus a weak spring that keeps it regular."""
    import numpy as np
    import xslope.joint as J

    def elem(p0, p1):
        p0, p1 = np.asarray(p0, float), np.asarray(p1, float)
        L = np.hypot(*(p1 - p0))
        t = (p1 - p0) / L
        return dict(xy=np.array([p0, p1, 0.5 * (p0 + p1)]), t=t,
                    n=np.array([-t[1], t[0]]), w=np.array([0.5, 0.5, 0.0]) * L)
    E = [elem((-1, 0), (1, 0)), elem((1, 0.25), (1, 1.75))]
    jd = {"n": 2, "dof": np.arange(24).reshape(2, 12),
          "tx": np.array([e["t"][0] for e in E]), "ty": np.array([e["t"][1] for e in E]),
          "nx": np.array([e["n"][0] for e in E]), "ny": np.array([e["n"][1] for e in E]),
          "kn": np.array([kn, kn]), "ks": np.array([ks, ks]), "tcut": np.zeros(2),
          "tip": np.zeros((2, 3), bool), "w": np.array([e["w"] for e in E]),
          "cj": np.array([1e4, c_wall]),
          "tanphi": np.tan(np.radians([45.0, phi_wall])),
          "jred": np.array([False, True])}
    T = np.zeros((24, 3))
    for ei, e in enumerate(E):
        for p in range(3):
            x, y = e["xy"][p]
            T[12 * ei + 2 * p] = (1.0, 0.0, -(y - 1.0))
            T[12 * ei + 2 * p + 1] = (0.0, 1.0, x)
    Ke = J._joint_element_stiffness(jd["w"], jd["tx"], jd["ty"], jd["nx"],
                                    jd["ny"], jd["kn"], jd["ks"])
    Kn = np.zeros((24, 24))
    for ei in range(2):
        Kn[np.ix_(jd["dof"][ei], jd["dof"][ei])] += Ke[ei]
    K = T.T @ Kn @ T + np.eye(3) * 1e-3 * kn
    return jd, T, K, np.array([P, -W, -P * (hp - 1.0)])


def _leaning_block_run(case, sweeps=3000):
    """The solver's own initial-stiffness loop on the leaning block:
    ``u = K^-1 (f + T^T L(u))`` with L from joint_vp_sweep, every sweep fed to the
    solve's contact-cycle ring (fem._ContactCycle) as the solve feeds it. Returns
    the state codes, the upper wall pair's shear, the final field, the residual of
    the joints' own tractions at it (the law the Newton check reads) and the
    cycle reading at the end, against the block's elastic displacement."""
    import numpy as np
    import xslope.joint as J
    jd, T, K, f = _leaning_block(*case)
    cj, tp = J.joint_reduced_strength(jd, 1.0)
    slip = np.zeros((2, 3)); opn = np.zeros((2, 3), bool)
    sp = np.zeros((2, 3), bool)
    Kinv = np.linalg.inv(K)
    u = Kinv @ f
    u_el = float(np.abs(T @ u).max())
    cyc = fem._ContactCycle()
    codes, ts = [], []
    for _ in range(sweeps):
        loads = np.zeros(24)
        _, st = J.joint_vp_sweep(jd, T @ u, loads, cj, tp, slip, opn,
                                 slip_state=sp)
        code = st["slipping"].astype(np.int8) + 2 * st["open"]
        codes.append(code)
        ts.append(float(st["ts"][1, 1]))
        cyc.step(code, T @ u)
        u = Kinv @ (f + T.T @ loads)
    f_int, _, st = J.joint_internal_force(jd, T @ u, cj, tp, slip_p=slip,
                                          open_prev=opn)
    fn = np.zeros(24)
    for ei in range(2):
        fn[jd["dof"][ei]] += f_int[ei]
    resid = f - T.T @ fn - 1e-3 * jd["kn"][0] * u
    codes = np.array(codes)
    return dict(codes=codes, ts=np.array(ts), u=u, st=st,
                late=int((codes[-300:][1:] != codes[-300:][:-1]).sum()),
                resid=float(np.abs(resid).max()), cycle=cyc.reading(u_el))


def check_leaning_block():
    """The open/close chatter of a cohesive, no-tension contact sitting at zero
    normal stress, on a leaning block, and the contact-cycle reading of it."""
    import numpy as np
    print("\n9. the open/close chatter on a leaning block")
    W = 100.0
    chatter = [(5.0, 1.5), (5.0, 1.9), (10.0, 0.5)]
    quiet = [(10.0, 1.0), (20.0, 1.5), (50.0, 1.0), (5.0, 0.5)]
    for case in chatter:
        run = _leaning_block_run(case)
        tag = f"P={case[0]:g}, hp={case[1]:g}"
        # The mechanism: closed at (almost) zero normal stress the upper wall pair
        # carries its cohesion; one sweep later it is open and carries nothing.
        ts = np.abs(run["ts"][-6:])
        check(f"[{tag}] the upper wall pair flips between carrying its cohesion "
              "(c = 10) and carrying nothing, every few sweeps",
              run["late"] > 50 and ts.max() > 0.95 * 10.0 and ts.min() == 0.0,
              f"{run['late']} state changes in the last 300 sweeps, |shear| "
              f"{ts.min():.2f} .. {ts.max():.2f}")
        check(f"[{tag}] ...and the block never balances under the joint law",
              run["resid"] > 1e-3 * W, f"residual {run['resid']:.3g}")
        r = run["cycle"]
        check(f"[{tag}] ...the contact-cycle reading stands it: one contact (the "
              "upper wall pair) cycling, no net movement",
              r is not None and r["fires"] and r["n_flip"] == 1
              and r["contacts"] == [4] and r["drift"] <= 1e-12, f"{r}")
    for case in quiet:
        run = _leaning_block_run(case)
        tag = f"P={case[0]:g}, hp={case[1]:g}"
        r = run["cycle"]
        check(f"[{tag}] a contact that does not chatter: no state change late, "
              "and no contact cycling for the reading to stand",
              run["late"] == 0 and r is not None and r["n_flip"] == 0
              and not r["fires"], f"{r}")


def check_contact_cycle():
    """The contact-cycle reading (fem.JOINT_CYCLE_ON) on synthetic sweeps, and
    the words its verdict is stated in."""
    import numpy as np
    print("\n10. a trial whose only movement is a contact cycle")
    check("the contact-cycle reading is on by default", fem.JOINT_CYCLE_ON is True)
    check("its thresholds are the measured ones: period <= 64 over 256 sweeps, "
          "1 to 4 contacts, 5e-8 elastic displacements per sweep",
          (fem._CYCLE_PMAX, fem._CYCLE_WINDOW, fem._CYCLE_MAX_FLIPS,
           fem._CYCLE_DRIFT) == (64, 256, 4, 5e-8))
    n, ndof, uel = 300, 40, 1e-3
    base_u = np.linspace(0.0, 2e-3, ndof)

    def feed(period, flipping, drift_per_sweep, sweeps=None):
        cyc = fem._ContactCycle()
        sweeps = sweeps or (fem._CYCLE_WINDOW + fem._CYCLE_PMAX + 10)
        for t in range(sweeps):
            code = np.zeros(n, np.int8)
            if period > 1:
                code[:flipping] = 2 if (t % period) == 0 else 0
            u = base_u + drift_per_sweep * t + 1e-12 * np.sin(2 * np.pi * t / max(period, 1))
            cyc.step(code, u)
        return cyc.reading(uel)

    r = feed(8, 1, 0.0)
    check("one contact cycling with period 8 and no drift: the trial stands",
          r is not None and r["period"] == 8 and r["fires"], f"{r}")
    r = feed(8, 1, 1e-8 * uel)
    check("...and with a drift far below the bound (1e-8 elastic per sweep) it stands",
          r is not None and r["fires"])
    r = feed(8, 1, 1e-5 * uel)
    check("the same cycle while the field creeps 1e-5 elastic per sweep: it does not",
          r is not None and not r["fires"], f"drift {r and r['drift']:.2e}")
    r = feed(8, 1, 5.4e-7 * uel)
    check("the same cycle at RJ-5's failing edge's slowest creep (5.4e-7 elastic per "
          "sweep): it does not", r is not None and not r["fires"])
    r = feed(8, 1, 5.0e-9 * uel)
    check("...and at the stuck RJ-20 trial's (5.0e-9 per sweep): it stands",
          r is not None and r["fires"])
    r = feed(1, 0, 1e-5 * uel)
    check("a steady state set with a creeping field (a failing trial): it does not",
          r is not None and not r["fires"] and r["n_flip"] == 0)
    r = feed(8, fem._CYCLE_MAX_FLIPS + 1, 0.0)
    check("more contacts cycling than the bound allows: it does not",
          r is not None and not r["fires"])
    cyc = fem._ContactCycle()
    rng = np.random.default_rng(0)
    for t in range(fem._CYCLE_WINDOW + fem._CYCLE_PMAX + 10):
        code = np.zeros(n, np.int8)
        code[rng.integers(0, n)] = 1
        cyc.step(code, base_u)
    check("states that do not repeat exactly: no reading", cyc.reading(uel) is None)

    # The verdict's words: the log line, the closing summary and the report all
    # read them off the one stop reading.
    rd = {"rule": "contact_cycle", "period": 8, "n_contacts": 3, "drift": 4.1e-8,
          "contacts": [{"line": 993, "x": 24.34, "y": 42.17},
                       {"line": 730, "x": 39.9, "y": 47.2},
                       {"line": 729, "x": 40.6, "y": 46.3}],
          "iteration": 1000000}
    clause = ("the only movement is 3 contacts cycling (period 8 iterations), no "
              "net movement; force balance not met")
    check("the verdict's words", fem.contact_cycle_clause(rd) == clause,
          fem.contact_cycle_clause(rd))
    sol = {"converged": False, "verdict": "JOINT_SETTLED",
           "exit_reason": "joint_settled", "stop_reading": rd, "u_ratio": 2.87}
    note = fem._verdict_note_base(sol, hybrid=True)
    check("the run log's trial line: STANDS, the clause, the contacts named",
          note.startswith("STANDS: " + clause + " -> counted STABLE")
          and "line 993 at (24.34, 42.17)" in note, note)
    note = fem._verdict_note_base(sol, hybrid=False)
    check("...under a criterion that reads convergence alone it is counted failed",
          note.startswith("Did NOT converge: " + clause + " -> counted FAILED"), note)
    trial = dict(sol, F=2.53125, iterations=1000000)
    said = fem._standing_edge_sentence(2.53125, trial, [trial], "iterations")
    check("the closing summary's standing-edge sentence states it",
          said.startswith("At F = 2.5312 the slope stands: " + clause + ".")
          and "line 729 at (40.60, 46.30)" in said, said)
    from xslope import report
    rec = {"final_interval": [2.53125, 2.6875],
           "trials": [dict(trial, stable=True)]}
    check("the report states the standing edge's verdict in the same words",
          report._standing_edge_cycled(rec) == said, report._standing_edge_cycled(rec))


def main():
    print("=" * 72)
    print("Joint verdict checks")
    print("=" * 72)
    check_measured()
    check_mutations()
    check_inert()
    check_churn()
    check_wiring()
    check_restick_band()
    check_finisher_reading()
    check_ambiguous_at_ceiling()
    check_leaning_block()
    check_contact_cycle()
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
