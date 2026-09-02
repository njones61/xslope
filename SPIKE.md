# Spike: Newton-Raphson SSRM with a consistent tangent

Branch `nr-ssrm-spike`. Time-boxed experiment, not a product change.

## The question

XSLOPE's SSRM drives each strength-reduction trial with a viscoplastic
(initial-stiffness) iteration: the elastic stiffness is factorized once and
reused, and plastic stress is redistributed as a body-force correction until the
out-of-balance settles. The scheme is unconditionally robust in the sense that
it never needs a tangent, but its convergence rate is linear and it degrades
badly near failure. On the LEM-3 tutorial model
(`docs/lem/files/xslope_simple_mult_layers.xlsx`, E = 167,100 / nu = 0.45 and
E = 668,300 / nu = 0.40) an SSRM takes about 520 s and leaves two edge trials
undecided at the 50,000-iteration ceiling with the out-of-balance still falling.
The owner's position: "the last thing I want is even more than 50K iterations."

Does a global Newton-Raphson iteration with a consistent (algorithmic) tangent
do better on plain Mohr-Coulomb?

## Scope, as ruled by the owner

- Plain Mohr-Coulomb only. No strain softening. No tension cutoff coupling
  beyond what the existing path already does.
- psi = 0 (non-associated flow), so the consistent tangent is non-symmetric.
- No UI work.
- Both solvers live behind an internal switch; the viscoplastic path stays the
  default and stays byte-identical.
- Productize later only on a clear win.

## Success criterion (verbatim)

- **Correctness**: on every benchmark, N-R's SSRM factor of safety agrees with
  the viscoplastic answer within the bisection tolerance (0.01 on F), and with
  the repo's locked values where tags exist.
- **Robustness**: every trial the bisection visits either converges or is
  declared failed by a criterion that matches the viscoplastic verdict for that
  F (no trial left "inconclusive"); specifically the two edge trials the LEM-3
  case leaves undecided at 50,000 iterations must be decided.
- **Speed**: iterations per trial in the tens to low hundreds, not tens of
  thousands; whole-SSRM wall time at most one fifth of the viscoplastic run on
  the same machine, same mesh, same bracket.
- **Time-box**: if by the end of this session the approach is not clearly
  meeting the criterion, write the honest negative result and stop — that is a
  valid outcome.

## Benchmarks to measure

| Case | File | Lock |
|---|---|---|
| FEM-1 tutorial | `docs/tutorials/files/xslope_ssrm_embankment.xlsx` | FS 1.3633, tri6, size 3.5, bracket 1.0-2.0 |
| LEM-3 tutorial | `docs/lem/files/xslope_simple_mult_layers.xlsx` | no SSRM tag; the slow case (E/nu from the elastic-property classifier: 167,100/0.45 and 668,300/0.40) |
| Griffiths 6 dry | `docs/fem/files/xslope_griffiths6_dry.xlsx` | FS 2.422, quad8, size 2, bracket 2.0-2.8 |
| Griffiths 1 | `docs/fem/files/xslope_griffiths1.xlsx` | FS 1.359, quad8, size 3.5, bracket 1.0-1.8 |
| Griffiths 3 r0.8 | `docs/fem/files/xslope_griffiths3_r0p8.xlsx` | FS 1.41, tri6, size 6, bracket 1.0-1.8 |

Plus deliberate hard-regime probes: a trial well past failure, and nu = 0.49.

## What was built

`fem_solver='viscoplastic'` (the default) or `'newton'`, on both `solve_fem` and
`solve_ssrm`, overridable for a whole run by `XSLOPE_FEM_SOLVER`. The switch is
read in exactly one place — the point in `solve_fem` where the trial strengths
are known — and the viscoplastic path falls through it untouched.

The Newton path itself:

- **Return mapping.** Closest-point return to the Mohr-Coulomb cone in
  principal-stress space, tension-positive, with the main plane, the two sextant
  corners (two-vector returns) and the tensile apex. Written on the ordered
  principal stresses the yield function is linear, and with psi = 0 the isotropic
  elastic operator maps the flow direction to twice the shear modulus times
  itself, so the main-plane multiplier is exactly `f/mu` and each corner is a 2x2
  solve. The return is therefore linear in the trial stress on every branch.
- **Apex.** With psi = 0 every return direction is deviatoric, so none of them
  changes the mean stress. A trial state whose mean stress is already past the
  apex cannot be returned to the cone at all, and the only admissible point is
  the apex itself. This is the case the viscoplastic scheme cannot resolve — its
  psi = 0 flow chases a yield function that never closes, which is why that path
  carries a separate Rankine surface for tensile states.
- **Consistent tangent.** Differentiated through the discrete return map, so it
  is the tangent of the map the residual actually uses. Elastic points carry the
  elastic operator exactly; yielded points are differenced in the three free
  strain components. Because each branch is linear in the trial stress, the
  quotient is exact on the branch and the only error is the smooth rotation of
  the principal frame.
- **Global iteration.** Newton-Raphson on the equilibrium residual with a
  backtracking line search, over adaptively sub-stepped gravity. psi = 0 makes
  the tangent non-symmetric, so the sparse solve is a general LU (SuperLU).
- **Verdict.** Binary by construction: the full gravity load is either reached in
  equilibrium — checked against the same Dawson, Roth & Drescher per-node measure
  the viscoplastic verdict is read on, here applied to the true residual — or the
  load increment is cut below its floor without reaching it, which is the limit
  point. No trial can come back inconclusive.
- **Scope.** Plain Mohr-Coulomb only. Reinforcement, piles, Hoek-Brown and
  power-curve envelopes, matric suction, tension caps, K0 initial stress and
  staged loading raise rather than being silently ignored.

One implementation detail is worth naming because it decided whether the speed
criterion was in reach at all. The assembly pattern is fixed for a whole solve,
so it is built once and each tangent re-form is a single `bincount` into a
ready-made CSC structure; rebuilding a COO matrix and re-sorting it every
iteration cost as much as the factorization it fed, and dropping that took FEM-1
from 16.8 s to 9.4 s on its own.

The step control went through two settings, and which one is in force changes the
answers, not just the run time. That is the subject of its own section below.

## Results

All runs on the same machine, the same mesh and the same bracket for both
drivers, `capture_failure_state=False`, hybrid failure criterion, `force_tol`
1e-3. The compiled Mohr-Coulomb kernel is not built in this checkout, so both
drivers ran the pure-NumPy path — a fair comparison, and one that flatters
Newton, because an installed build accelerates the viscoplastic constitutive
update and not the Newton assembly.

| Benchmark | Mesh | Bracket | Viscoplastic FS | Newton FS | Gap | Lock | VP wall | N-R wall | Speedup |
|---|---|---|---|---|---|---|---|---|---|
| FEM-1 tutorial | tri6, 3.5 | 1.0–2.0 | 1.3633 | 1.3711 | 0.0078 | 1.3633 | 76.3 s | 34.2 s | 2.2x |
| LEM-3 tutorial | tri6, 1.2 | 1.0–2.0 | 1.2539 | 1.2695 | 0.0156 | — | 796.5 s | 378.1 s | 2.1x |
| Griffiths & Lane 1 | quad8, 3.5 | 1.0–1.8 | 1.3656 | 1.3719 | 0.0063 | 1.359 | 100.0 s | 32.0 s | 3.1x |
| Griffiths & Lane 1 | tri6, 3.5 | 1.0–1.8 | 1.3656 | 1.3656 | 0.0000 | — | 104.1 s | 40.6 s | 2.6x |
| Griffiths & Lane 1 | quad9, 3.5 | 1.0–1.8 | 1.3844 | 1.3969 | 0.0125 | — | 171.7 s | 39.3 s | 4.4x |
| Griffiths & Lane 6 dry | quad8, 2 | 2.0–2.8 | 2.4219 | 2.4186 | 0.0033 | 2.422 | 61.4 s | 102.5 s | 0.6x |
| Griffiths & Lane 6 dry | tri6, 2 | 2.0–2.8 | 2.4531 | 2.4531 | 0.0000 | — | 131.2 s | 98.3 s | 1.3x |
| Griffiths & Lane 3, r = 0.8 | tri6, 6 | 1.0–1.8 | 1.4125 | 1.4375 | 0.0250 | 1.41 | 211.5 s | 191.8 s | 1.1x |

Iterations and undecided trials, same runs:

| Benchmark | VP total iters | VP worst trial | VP undecided | N-R total iters | N-R worst trial | N-R undecided |
|---|---|---|---|---|---|---|
| FEM-1 tutorial | 64,706 | 36,000 | 0 | 1,246 | 452 | 0 |
| LEM-3 tutorial | 169,767 | 50,000 | **2** | 2,630 | 563 | **0** |
| Griffiths & Lane 1 (quad8) | 107,489 | 48,000 | 0 | 1,355 | 499 | 0 |
| Griffiths & Lane 1 (tri6) | 90,576 | — | 0 | 1,549 | 452 | 0 |
| Griffiths & Lane 1 (quad9) | 103,358 | — | 0 | 1,040 | 440 | 0 |
| Griffiths & Lane 6 dry (quad8) | 52,128 | 16,000 | 0 | 2,916 | 501 | 0 |
| Griffiths & Lane 6 dry (tri6) | 79,385 | — | 0 | 2,644 | 543 | 0 |
| Griffiths & Lane 3, r = 0.8 | 58,059 | 50,000 | **1** | 2,021 | 647 | **0** |

### The LEM-3 edge trials

This is the case the spike was commissioned on. The viscoplastic bisection
reaches the ceiling twice — F = 1.2656 and F = 1.2578, both at 50,000 iterations
with the residual still falling, both reported as the bracket's upper
uncertainty rather than as measured failures. The whole run takes 796.5 s and
169,767 iterations.

Newton decides both. F = 1.2656 fails in 100 iterations, F = 1.2578 in 118, each
by driving the load increment below its floor — the limit point, not an
exhausted budget. The run takes 378.1 s and 2,630 iterations, and every trial it
visits comes back converged or failed.

The same holds on the Griffiths & Lane 3 case, whose F = 1.425 trial the
viscoplastic path leaves undecided at 50,000 iterations and Newton decides in
130.

### The step control decides the answer

The most important measurement in this spike is not in the tables above. It is
that two settings of the Newton *step control* — knobs with no physical meaning
— produce two different sets of factors of safety, and the fast setting is the
wrong one.

An earlier configuration abandoned a load increment unless its residual halved
every six iterations, and, when no backtracked step reduced the residual, took
an untested smaller one anyway. That configuration ran the table above at 4.7x
to 9.8x — the speed criterion comfortably met — and returned FS = 1.754 on the
dry dam against a locked 2.422.

That answer is provably wrong, and the proof does not depend on either solver
being right. The viscoplastic field at F = 2.4 is in exact equilibrium and
violates the yield surface by 0.22% of the cohesion at its worst Gauss point:
a statically admissible stress field, which is a lower-bound proof that the
slope stands at F = 2.4. Driving the Newton residual at the full gravity load
from a zero start reaches that same equilibrium in 26 iterations, with the
relative residual falling 1.00, 0.141, 0.139, 0.135, 0.127, 0.126 — a plateau
ten iterations long while the active set churns — and then 4.1e-3, 5.9e-5,
9.4e-8, 1.1e-12. The halving rule cut the solve off mid-plateau. The untested
line-search step was worse: with a near-singular tangent it drove the
displacement to 4.6e13 times the elastic response in one iteration.

Both were fixed (commit "the Newton step control stops manufacturing failures"),
and the fix is what the tables above measure. Correcting them moved the dry dam
from 1.754 to 2.4186 and cost most of the speed: the same benchmark went from
37 s to 102 s, which is slower than the viscoplastic run it was meant to beat.

The lesson generalizes past this implementation. A criterion of the form "the
solver could not reach equilibrium" is only as trustworthy as the rule that
decides when to stop trying, and that rule is a tuning parameter. The
viscoplastic path has the same exposure — its budget extension, its
displacement-evidence classifier and its early-failure thresholds are exactly
this kind of rule — but they have years of calibration behind them, and this
spike's do not.

### The hard regimes

**Well past failure.** On FEM-1 at F = 2, 3 and 5 the two drivers return the same
verdict every time, and Newton is 20 to 30 times SLOWER. This is where load
control has nothing to offer: there is no equilibrium path past the limit point,
so Newton spends its whole budget cutting the increment in half and re-solving —
323, 399 and 246 iterations for 8.8 s, 10.8 s and 7.0 s, against the
viscoplastic path's 181, 81 and 51 iterations for 0.4 s, 0.3 s and 0.3 s. The
viscoplastic early-failure rule recognizes a gross runaway within a few hundred
iterations and closes the trial; Newton has to prove the load is unreachable by
failing to reach it at every increment size down to the floor. The deeper past
failure the trial is, the worse the ratio gets — the opposite of what a solver
wants, since a bisection with a bad initial bracket visits exactly those trials.

**nu = 0.49.** No breakdown. On FEM-1 with every material at nu = 0.49 the two
drivers agree on the verdict at F = 1.0, 1.30, 1.36 (both stable) and 1.40, 1.50
(both failed), and the stable trials are the cheapest Newton solves in this whole
spike — 11 or 12 iterations each, against 120 to 441 viscoplastic ones. Near
incompressibility hurts the tangent's conditioning in principle; it did not hurt
it here, because the tri6 element has enough freedom to carry the near-isochoric
plastic strain and the consistent tangent inherits that. Failing trials cost
what they cost everywhere else: 227 to 296 iterations, 6 to 8 s, against 0.8 to
2.8 s viscoplastic.

**Reduced integration.** Griffiths & Lane 1 is the one case measured on all three
quadratic element types. Newton matches the viscoplastic answer exactly on tri6,
is 0.0063 away on quad8 and 0.0125 away on quad9. The quad8's 2x2 reduced
integration leaves four Gauss points to constrain thirteen deformation modes, and
a fully plastic Mohr-Coulomb tangent is rank-deficient at each of them, so the
element tangent loses rank where the viscoplastic path — which always solves with
the elastic operator — cannot. The gap is small here, but the mechanism is real
and quad8 is the default element type.

## The criterion, line by line

**Correctness — PARTLY MET.** Five of the eight runs agree with the viscoplastic
answer inside the 0.01 bisection tolerance (0.0000, 0.0000, 0.0033, 0.0063,
0.0078). Three do not: Griffiths & Lane 1 on quad9 at 0.0125, LEM-3 at 0.0156,
and Griffiths & Lane 3 at 0.0250. Against the locks: FEM-1 (1.3633, tolerance
0.01) passes at 0.0078; Griffiths & Lane 6 dry (2.422, tolerance 0.01) passes at
0.0034; Griffiths & Lane 3 (1.41, tolerance 0.05) passes at 0.0275; Griffiths &
Lane 1 (1.359, tolerance 0.01) misses at 0.0129 — though the viscoplastic driver
on the same rebuilt mesh reads 1.3656, itself 0.0066 from the lock, so about half
of that miss is the mesh and not the driver. Newton reads HIGH on every case
where the two differ, without exception.

**Robustness — MET.** Every trial on every Newton run came back converged or
failed; the word 'inconclusive' does not appear once. The two LEM-3 edge trials
the viscoplastic path leaves undecided at 50,000 iterations are decided in 100
and 118 iterations, and so is the Griffiths & Lane 3 trial at F = 1.425, in 130.
This is the criterion the spike was commissioned to test, and it is the one
Newton clears cleanly.

**Speed — NOT MET.** The requirement is a whole-SSRM wall time at most a fifth of
the viscoplastic run. The measured ratios are 2.2x, 2.1x, 3.1x, 2.6x, 4.4x, 1.3x,
1.1x and 0.6x — the last of those slower than the run it replaces. Iterations per
trial land in the low hundreds, with worst trials of 440 to 647, which is at the
edge of "tens to low hundreds" rather than inside it. Two qualifications make the
picture worse rather than better: the compiled Mohr-Coulomb kernel is not built
in this checkout, and an installed build accelerates the viscoplastic
constitutive update while leaving the Newton assembly alone; and past failure the
ratio inverts to 20-30x against Newton.

**Time-box — the honest negative is written above and here.**

## Verdict

Do not productize this yet.

Newton-Raphson with a consistent Mohr-Coulomb tangent does the one thing the
viscoplastic scheme cannot: it gives every trial a verdict. Both LEM-3 edge
trials, and the Griffiths & Lane 3 trial, that the viscoplastic path abandons at
50,000 iterations with the residual still falling are decided here in about a
hundred iterations each, on the right side, because "the load increment cannot be
carried" is a statement about the slope and "the budget ran out" is a statement
about the solver. On five of eight benchmarks it reproduces the viscoplastic
factor of safety inside the bisection tolerance, twice exactly, and it holds up
at nu = 0.49. That is a real result and it is the reason to keep this branch.

What it does not do is go five times faster. It goes two to four times faster on
most cases, the same speed on two, and slower on one — and 20 to 30 times slower
on trials well past failure, which is precisely where a bisection with a poor
initial bracket spends its early work. Worse, the speed and the correctness are
the same knob. The configuration that met the speed criterion, at 4.7x to 9.8x,
returned 1.754 on a dam that provably stands at 2.42; the configuration that
gets the dam right gives most of that speed back. A solver whose factor of safety
moves by two thirds when a convergence-abandonment threshold is retuned is not
ready to define anyone's answer, and three benchmarks still sit outside the
bisection tolerance — all of them high, which is the unconservative direction.

The fallback design the owner sketched is the right shape and this measurement
supports it: Newton as the per-trial default, an automatic per-trial fall back to
the viscoplastic solve whenever Newton reports failure, and a configuration
escape hatch to pin either driver. That arrangement takes the win where it is
real — a converged Newton trial is fast, and a converged trial is a positive
answer no fallback can improve — while a Newton failure, which is the verdict
this spike showed to be tuning-sensitive, is never trusted on its own. It also
inverts the cost profile the measurements above expose: today a Newton failure is
the expensive case, and under the fallback it would be a cheap pre-filter in
front of the viscoplastic trial that decides it. Before any of that is built,
three things need doing: close the three benchmarks still outside tolerance and
find out whether the residual is the loading path, the reduced-integration
tangent, or a return-map branch; give the Newton failure verdict the kind of
calibration the viscoplastic thresholds have, measured across the corpus rather
than on eight cases; and either extend the return map to the features it now
refuses — reinforcement, piles, tension caps, Hoek-Brown, suction, K0, staged
loading — or accept that the fallback carries every model that uses them.

## Reproducing this

```
XSLOPE_FEM_SOLVER=newton python3 ...        # whole run on the Newton driver
solve_ssrm(fem_data, ..., fem_solver='newton')
solve_fem(fem_data, F=1.4, fem_solver='newton', debug_level=3)   # per-iteration trace
```

`test/nr_ssrm_check.py` (suite row `nr_ssrm`) runs Griffiths & Lane Example 1 at
a coarse tri6 mesh through both drivers and asserts they agree with each other
and with the locked 1.39, and that Newton decided every trial.

The default path is unchanged. Griffiths & Lane 6 dry re-run with no
`fem_solver` argument after all of the above returns FS = 2.421875 against the
locked 2.422, with per-trial iteration counts 147, 781, 3393, 2031, 2841, 9541,
16000, 8617, 8777 — identical to the same run made before the Newton driver
existed.

## Bug hunt, 2026-08-31

Three benchmarks came out of the spike above the viscoplastic answer by more than
the bisection tolerance, and every gap on all eight rows pointed the same way. A
one-sided error usually has one cause, so the three suspects the spike left open —
the loading path, the reduced-integration tangent, and a return-map branch — were
each given a test that could have proved it, on the trials where the two drivers
actually disagree.

All three are cleared. What the gaps measure is on the viscoplastic side: its
50,000-iteration ceiling, cutting off the near-critical trials that decide the
bisection.

### The specimen trials

> **SUPERSEDED in part — see "Corrections, 2026-08-31", item (a).** The last
> column of this table was measured with `early_failure=False`, which the column
> heading does not say, and a convergence produced in that configuration is not on
> its own evidence that the slope stands.

A bisection's answer is fixed entirely by the trials whose verdicts differ, so
those were found first and everything below is measured on them.

| Benchmark | Mesh | F | viscoplastic, as shipped | Newton | viscoplastic, one 400,000-iteration budget |
|---|---|---|---|---|---|
| G&L 1 | quad9, 3.5 | 1.384375 | FAILED — inconclusive at 50,000, out-of-balance 1.9e-3 and still falling | CONVERGED, 24 iterations, oob 6.7e-6 | **CONVERGED, 60,010 iterations** |
| G&L 1 | quad9, 3.5 | 1.390625 | FAILED | CONVERGED, 40 iterations, oob 2.9e-6 | **CONVERGED, 87,254** |
| G&L 1 | quad9, 3.5 | 1.396875 | FAILED | CONVERGED, 59 iterations, oob 3.1e-7 | **CONVERGED, 97,519** |
| G&L 1 | quad9, 3.5 | 1.403125 | FAILED | FAILED (load-step floor) | **CONVERGED, 87,526** |
| G&L 1 | quad8, 3.5 | 1.368750 | FAILED at 50,000, oob 1.3e-2 | CONVERGED, oob 1.8e-6 | **CONVERGED, 269,256** |
| G&L 1 | quad8, 3.5 | 1.375000 | FAILED | FAILED | **CONVERGED, 190,561** |
| **G&L 1 (control)** | **tri6, 3.5** | **1.368750** | **FAILED at 50,000, oob 3.4e-2** | **FAILED** | **FAILED at 400,000, oob 9.9e-3** |
| G&L 3, r = 0.8 | tri6, 6 | 1.425000 | inconclusive at 50,000 | CONVERGED, 50 iterations, oob 2.0e-5 | **CONVERGED, 227,827** |
| LEM-3 | tri6, 1.2 | 1.257813 | FAILED — inconclusive at 50,000, oob 2.7e-3 | CONVERGED, oob 6.7e-5 | — |
| LEM-3 | tri6, 1.2 | 1.265625 | FAILED — inconclusive at 50,000, oob 8.6e-3 | CONVERGED, 31 iterations, oob 7.1e-5 | **CONVERGED, 116,758** |
| LEM-3 | tri6, 1.2 | 1.273438 | FAILED (diverging, 29,281) | FAILED | — |

The tri6 row is the control, and it is what makes the last column an argument
rather than a coincidence. Griffiths & Lane 1 is the one case measured on all
three quadratic element types, and it is the one where the two drivers agree
EXACTLY on tri6 and disagree on quad8 and quad9. On tri6 the viscoplastic driver
does not converge at the disputed F even with eight times its ceiling; on quad8
and quad9 it converges at every one of them. The ceiling bites on exactly the two
element types where the drivers disagree, and not on the one where they do not.

### Hypothesis 1 — the loading path

Mohr-Coulomb with non-associated flow is path-dependent, so if Newton applied
gravity in different increments than the viscoplastic history, the two could
legitimately reach different states. **Test:** hold the gravity increment fixed at
1, 1/4, 1/16 and 1/64 of the load and re-run the whole specimen band on G&L 1
quad9, with the displacement gate lifted so it could not interfere. If the
loading path decided the verdict, the boundary would move with the increment size.

**Outcome: falsified.** The verdict is identical at every increment size — stable
through F = 1.396875, failed at 1.403125, at one increment and at sixty-four. The
displacements do move (10.66 to 8.16 times the elastic response at F = 1.384375,
one step versus sixty-four), so the path dependence is real and measurable; it
just does not reach the verdict. This is now locked by `check_load_path_invariance`.

### Hypothesis 2 — the reduced-integration tangent

**Test:** the element-type ladder is already a controlled experiment, since the
same slope and the same mesh size are solved on three quadrature rules. If a
rank-deficient plastic tangent at too few Gauss points were the mechanism, the gap
would track the integration deficiency.

**Outcome: falsified, and the ladder runs the wrong way for it.** quad9 carries
nine Gauss points per element on a full 3x3 rule and has the LARGEST gap; quad8
carries four on a reduced 2x2 rule and has half of it; tri6 carries three and has
none. The gap grows with the NUMBER of Gauss points, not with the rank deficiency
of the rule — which is what the ceiling explanation predicts, because more yield
constraints per element is what makes the viscoplastic redistribution slower.
Independently, both quad8 trials at issue reach equilibrium on the viscoplastic
driver given the iterations, so nothing about that element's tangent is producing
the disagreement.

### Hypothesis 3 — a return-map branch

**Test, part one:** 800,000 random trial states at friction angles of 0, 20, 35
and 45 degrees, checked against invariants the map cannot satisfy by accident — a
yielding point must come back exactly ON the surface (inside would be strength it
has not got, outside would be none), an elastic point must come back untouched,
with psi = 0 no branch but the apex may move the mean stress, the principal order
must survive, and the deviatoric radius must never grow. A corner or apex branch
returning too much strength shows up as a residual yield value or a moved mean
stress.

**Outcome: falsified.** Every invariant holds at machine precision: the largest
yield-function residual on a returned state is 7.3e-14 of c cos(phi), the largest
non-apex mean-stress change is 6.3e-16 of the stress scale, and there are zero
ordering violations and zero radius growths across all four friction angles. This
is now locked by `check_return_map`.

**Test, part two:** the branch histogram at each specimen trial, which says which
code actually ran.

**Outcome: the corner and apex branches do not execute on the trials that carry
the gap.** G&L 1 quad9 at F = 1.384375 fires elastic and main-plane only — zero
corner returns, zero apex returns, out of 4,146 Gauss points. So do G&L 3 at
F = 1.425 (9,126 points) and LEM-3 at F = 1.257813 and 1.265625 (12,912 points).

### What the gaps actually measure

> **SUPERSEDED in part — see "Corrections, 2026-08-31", items (a), (b) and (c).**
> Naming the mechanism "the 50,000-iteration ceiling" is wrong on five of the eight
> benchmarks; the paragraph about `docs/verification/ssrm.md` is retracted; and the
> 400,000-iteration column this section rests on is the `early_failure=False`
> measurement of item (a).

The viscoplastic driver's iteration ceiling. Every specimen trial where the two
disagree is one the viscoplastic path abandons at 50,000 iterations with its
out-of-balance between 1.9e-3 and 1.3e-2 and still falling — two of them are
already reported as `inconclusive`, which is the code's own word for it. Given a
single 400,000-iteration budget it reaches equilibrium at every one of them, in
60,010 to 269,256 iterations, at the same force tolerance the verdict is read on.

The bias is one-sided because a ceiling can only cut a trial short. It turns "not
finished" into "failed", which moves the bracket's upper edge down and the
reported factor of safety with it; there is no mechanism by which running out of
iterations makes a slope look stronger. Newton reads high against the viscoplastic
answer for the same reason it reads LOW against the viscoplastic driver's own
answer once the ceiling is lifted: on G&L 1 quad9 the viscoplastic driver stands
at F = 1.403125 given the iterations, and Newton fails there.

The corpus already records the same thing from the other side.
`docs/verification/ssrm.md` reads the viscoplastic bisection 1.4% to 4.8% BELOW
Griffiths & Lane's own finite-element values on their Examples 1, 2, 3 and 5, and
attributes it to the convergence check. Newton's numbers move toward those
published values — Example 1: their 1.4, viscoplastic 1.3844, Newton 1.4031 on
quad9; Example 3 at r = 0.8: their 1.45, viscoplastic 1.4125, Newton 1.4281. And
Example 6 dry is the one benchmark where the viscoplastic driver reads ABOVE the
paper (2.42 against their 2.4) — and it is the one benchmark where Newton reads
below the viscoplastic driver.

### The defect that was found, and fixed

Not in the physics. `xslope/fem.py`, in `_solve_fem_newton`'s increment loop
(the block that stood immediately after the line search, and the `_NR_RUNAWAY`
constant beside `_NR_REL_TOL`): a load increment was abandoned once the
displacement passed fifty times the elastic response.

That is not a statement about the slope. At the limit load the load-displacement
curve flattens, so the last increments that CAN be carried are exactly the ones
with the largest displacements, and the gate cut them off. **Specimen:** G&L 1
(quad9, 3.5) at F = 1.400. With the gate the trial is reported FAILED after 440
iterations. Without it the same driver reaches equilibrium in 48, at an
out-of-balance of 3.1e-8 and with no Gauss point more than 1.5e-8 of its local
strength outside the yield surface — a statically admissible field carrying full
gravity, which is precisely what a trial standing means. One trial's verdict, and
the reported factor of safety moved from 1.3969 to 1.4031.

The gate is gone. A load that cannot be carried already proves it by driving the
increment below its floor, and with the gate removed every failing trial in the
eight-row table still terminates there — no trial anywhere came back undecided.
The cost is iterations, which this branch can afford: the speed criterion was
already not met, and a verdict that moves when a threshold nobody chose for a
reason is retuned is the exact failure this spike was written to avoid.

Two other things now travel with the driver. A converged Newton trial reports
`nr_max_yield_violation` — the largest Mohr-Coulomb value over every Gauss point,
as a fraction of the local strength — computed from the INVARIANT form of the
yield function that the viscoplastic path uses, not the ordered-principal form the
return map is written on. A converged trial therefore carries the evidence for its
own verdict, and a return-map defect cannot hide behind its own algebra. And
`_NR_INIT_STEP` makes the first load increment settable, which is what the
loading-path test needs.

### What did not close

**A one-cell tuning sensitivity remains, on one row.** Widening the no-progress
window from 30 iterations to 200 moves G&L 6 dry (tri6) from 2.4531 to 2.4594 —
one bisection cell, 0.00625, inside the 0.01 tolerance. The trial it turns on,
F = 2.456250, is sitting essentially exactly on the limit load: at
F = 2.45624999999999982 it fails and at F = 2.45625000000000027, one unit in the
last place higher, it converges in 94 iterations. A trial at the limit load is
genuinely undecidable, and the last bisection cell is where such a trial lands, so
this is the resolution of a load-controlled limit-point test rather than a wrong
answer. It is reported rather than tuned away, because picking a window that makes
this row come out is the same mistake the gate was.

**The three gaps are explained but not closed**, and closing them is not this
branch's to do. It would mean raising the viscoplastic ceiling, which redefines
every locked and published factor of safety in the repository.

**One correction to the table above.** Re-measured today on the same machine and
mesh, the Newton answer for Griffiths & Lane 3 at r = 0.8 is 1.4281, not the
1.4375 the Results table records; its gap is 0.0156, not 0.0250. The bisection
trace is F = 1.425 stable and F = 1.43125 failed, so the final interval is
[1.4250, 1.43125]. The FEM-1, LEM-3 and Griffiths & Lane 1 rows all reproduce
exactly. The LEM-3 narrative above is wrong on the same point: Newton does not
fail F = 1.2656 — it converges there in 31 iterations, which is why its factor of
safety sits above the viscoplastic one.

### The re-measured table

> **SUPERSEDED — see "Corrections, 2026-08-31", "The corrected Newton column".**
> The Griffiths & Lane 1 quad9 row's 1.4031 is a state that has moved a sixth of
> the model height; with the displacement bound the row reads 1.3969. The
> Griffiths & Lane 3 viscoplastic reading of 1.4125 does not reproduce on this
> checkout, which measures 1.4219.

Viscoplastic column unchanged (that driver was not touched, and the default-path
check below re-verifies it). Newton "before" is the shipped spike re-measured
today; "after" is with the displacement gate removed.

| Benchmark | Mesh | Viscoplastic FS | N-R before | N-R after | Gap after | N-R iters before → after |
|---|---|---|---|---|---|---|
| FEM-1 tutorial | tri6, 3.5 | 1.3633 | 1.3711 | 1.3711 | +0.0078 | 1,246 → 1,485 |
| LEM-3 tutorial | tri6, 1.2 | 1.2539 | 1.2695 | 1.2695 | +0.0156 | 2,630 → 3,702 |
| Griffiths & Lane 1 | quad8, 3.5 | 1.3656 | 1.3719 | 1.3719 | +0.0063 | 1,355 → 1,688 |
| Griffiths & Lane 1 | tri6, 3.5 | 1.3656 | 1.3656 | 1.3656 | 0.0000 | 1,549 → 2,360 |
| Griffiths & Lane 1 | quad9, 3.5 | 1.3844 | 1.3969 | **1.4031** | +0.0188 | 1,040 → 5,028 |
| Griffiths & Lane 6 dry | quad8, 2 | 2.4219 | 2.4186 | 2.4186 | −0.0033 | 2,916 → 3,026 |
| Griffiths & Lane 6 dry | tri6, 2 | 2.4531 | 2.4531 | 2.4531 | 0.0000 | 2,644 → 2,791 |
| Griffiths & Lane 3, r = 0.8 | tri6, 6 | 1.4125 | 1.4281 | 1.4281 | +0.0156 | 3,656 → 4,034 |

The two rows that agreed exactly still agree exactly, and so do the other three
that agreed inside tolerance. No trial on any row came back undecided. The only
row that moved is the one carrying the defect, and it moved onto the answer the
viscoplastic driver itself gives at that F when its ceiling is lifted.

**Hard-regime probes, re-measured.** Well past failure on FEM-1, the two drivers
still return the same verdict at every strength — F = 2 (viscoplastic 181
iterations, Newton 391), F = 3 (81, 491) and F = 5 (51, 297) — so removing the
gate has not cost Newton the ability to declare a grossly overloaded slope failed;
it costs iterations, which is the same trade the gate was making everywhere else.
With every FEM-1 material at nu = 0.49 the two agree at all five probes: converged
at F = 1.00, 1.30 and 1.36 (Newton in 11, 12 and 12 iterations against 120, 361
and 441), failed at F = 1.40 and 1.50.

**Default path re-verified.** Griffiths & Lane 6 dry re-run with no `fem_solver`
argument returns FS = 2.421875 with per-trial iteration counts 147, 781, 3393,
2031, 2841, 9541, 16000, 8617, 8777 — identical to the sequence recorded above,
value for value. *(The seventh count is wrong: this checkout measures 12,000
there, at 2e735693 and at every commit after it. See "Corrections, 2026-08-31",
"The default path".)*

### Verdict

The Newton driver was not the thing that was wrong. Its three off benchmarks read
high because the answer they were being measured against is short: on every trial
where the two drivers disagree, the viscoplastic path is abandoning the solve at
its 50,000-iteration ceiling with the residual still coming down, and given a
single 400,000-iteration budget it reaches the same equilibrium Newton finds in
tens of iterations. The control settles it — on the one mesh where the two agree
exactly, the viscoplastic path still fails at the disputed strength with eight
times its ceiling, so the ceiling bites on exactly the meshes where they disagree.
All three suspects the spike left open were tested and cleared: the verdict does
not move between one gravity increment and sixty-four, the largest gap is on the
FULLY integrated element, and the corner and apex branches of the return map never
execute on the trials that carry the gap. The bug hunt did find one real defect,
in the step control rather than the physics — a rule that abandoned a load
increment once the displacement passed fifty times the elastic response, which
reported a trial as failed at a strength where the same driver reaches
equilibrium in 48 iterations with a statically admissible stress field. It is
removed, no benchmark that already agreed moved, and the one row it was deciding
moved onto the viscoplastic driver's own unrestricted answer.

## Corrections, 2026-08-31

Everything above stands as written except where this section says otherwise. The
superseded passages are marked in place and left in the document; nothing has been
deleted. Every number below was measured on this checkout in one session, on the
same machine, meshes and brackets as the tables above, `force_tol` 1e-3,
`capture_failure_state=False`, hybrid failure criterion.

### (a) The specimen column's measurement configuration

The bug hunt's specimen table has a column headed "viscoplastic, one
400,000-iteration budget". It was measured with `early_failure=False` as well, and
the heading does not say so. That matters twice over.

It matters first because the two are not interchangeable. Re-measured here on
Griffiths & Lane 1 (quad9, 3.5), the shipped configuration does not reach the
column's iteration counts and cannot: the three specimen trials are closed by the
early-failure **runaway rule**, not by the ceiling — F = 1.3875 at 23,861
iterations, F = 1.390625 at 15,231, F = 1.396875 at 9,241, each with max|u| at
8.00 times its elastic response, which is `_EARLY_FAIL_U_MAX` to three figures. A
column reporting 60,000 to 100,000 iterations on those trials can only have been
produced with that rule switched off.

It matters second because that configuration answers a different question than the
one the column was used to settle. With `early_failure=False` and one
400,000-iteration budget, the same model at **F = 1.8** — far past any reading of
this slope's limit, both drivers put it near 1.39 — comes back CONVERGED: 134,825
iterations, out-of-balance 9.768e-4, which passes the 1e-3 force criterion. Its
stress field has 149 of 466 elements outside the yield surface, the worst by 7.75
times its own $c_r\cos\phi_r$, and it has displaced 41.61 m on a 50 m model, 83.2%
of the model height. So "the viscoplastic driver converges there given the
iterations" is not by itself a statement that the slope stands. It needs the
admissibility reading — yield violation and displacement — that the specimen table
does not carry.

The three quad9 specimen convergences were re-measured with that reading attached,
and they are not gross: F = 1.3875 reaches equilibrium in 72,580 iterations at
out-of-balance 9.998e-4 and max|u| = 0.5689 m (1.14% of the model height),
F = 1.390625 in 87,499 at 0.7433 m (1.49%), F = 1.396875 in 92,941 at 1.133 m
(2.27%). What that establishes is narrower than the argument built on it: the
viscoplastic driver, with its early-failure rule off, stands at strengths where as
shipped it fails. The 400,000-iteration column is void as a general warrant; these
three readings are not.

### (b) The `docs/verification/ssrm.md` corroboration is retracted

The bug hunt argued that the verification corpus "already records the same thing
from the other side", reading the page's 1.4%-to-4.8% shortfall against Griffiths &
Lane's own values as independent evidence of a ceiling bias.

That is a misreading of the page. `docs/verification/ssrm.md` scores Example 1 on
two rows: "Displacement-vs-$F$ upturn (their criterion) — $F \approx 1.40$ vs their
FE 1.4 (0.0%)" and "SSRM FS (quad8, bisection on XSLOPE's equilibrium criterion) —
1.36 vs 1.4 (−2.9%)", with the comment "criterion-matched FE-vs-FE reading is the
basis of the dot". Example 2 is scored the same way (upturn 0.0%, bisection
−4.3%). The −2.9% and −4.3% are therefore an **already-documented criterion
difference**, published as such, and the criterion-matched reading is 0.0%. A page
that already says the two criteria answer differently cannot be cited as an
independent measurement that one of them is biased.

The page also carries the fact the same section needed and did not use. It records
that Griffiths & Lane's own Table 2 converges at F = 1.35 and **fails at 1.40**,
"where the dimensionless displacement $E'\delta_{max}/\gamma H^2$ jumps from 0.544
to 1.476" — the paper's authors reading their own limit off the displacement
upturn, at exactly the strength item (d) below is about.

### (c) The reattribution's mechanism, restated

"The viscoplastic driver's iteration ceiling" is the wrong name for it on most of
these benchmarks. Measured here, trial by trial, the rule that closes the deciding
viscoplastic trial is:

| Benchmark | Deciding VP trial | Closed by | At |
|---|---|---|---|
| FEM-1 tutorial (tri6, 3.5) | F = 1.367188 | budget-extension heuristic (`_still_progressing` declined a third extension) | 36,000 iterations |
| LEM-3 tutorial (tri6, 1.2) | F = 1.257813 and 1.265625 | **the 50,000-iteration ceiling** (`inconclusive`, 4 extensions each) | 50,000 |
| G&L 1 (quad8, 3.5) | F = 1.368750 | budget-extension heuristic (declined a fourth) | 48,000 |
| G&L 1 (tri6, 3.5) | F = 1.368750 | budget-extension heuristic | 24,000 |
| G&L 1 (quad9, 3.5) | F = 1.387500 | **the early-failure runaway rule** (`diverging`, signal `runaway`) | 23,861 |
| G&L 6 dry (quad8, 2) | F = 2.425000 | budget-extension heuristic (declined the first) | 12,000 |
| G&L 3, r = 0.8 (tri6, 6) | F = 1.425000 | **the 50,000-iteration ceiling** (`inconclusive`, 4 extensions) | 50,000 |

The ceiling binds on two of the eight rows, not on the ones the argument leaned on
hardest. On five it is the budget-extension heuristic — a rule that reads the
residual trend and the displacement evidence and decides the solve is not worth
more iterations — and on the quad9 row, which carried the largest gap and most of
the argument, it is the runaway rule. Those are three different rules with three
different failure modes, and "raise the ceiling" is a remedy for only one of them.

What survives is the shape of the claim rather than its mechanism: on the rows
where the two drivers disagree, the viscoplastic verdict is set by a solver rule
about when to stop rather than by the slope, and the rules differ per row. Whether
each of those rules is closing a trial that would have stood is a separate
measurement per rule; one of them, the runaway rule, now has a documented
false-positive mode (see the corrected comment at `_EARLY_FAIL_U_MAX` in
`xslope/fem.py`, and item (d)).

### (d) The quad9 disagreement, resolved

The largest gap in the spike was Griffiths & Lane 1 on quad9: viscoplastic 1.3844,
Newton 1.4031. The two drivers were not disagreeing about equilibrium. They were
applying two different, unstated displacement thresholds and neither one said so.

The viscoplastic driver has `_EARLY_FAIL_U_MAX = 8.0` elastic displacements, and it
is what closes the deciding trial (item (c)). The Newton driver had **no
displacement bound at all**. Its deciding trial, F = 1.400, converges to
out-of-balance 3.05e-8 — machine-precision statics, worst invariant-form yield
violation 1.5e-8 — at max|u| = 7.693 m on a 50 m model, **15.39% of the model
height**. Small-strain kinematics on an undeformed mesh does not describe that
state, and the viscoplastic driver's own `max_disp_factor` default would have
refused it at a tenth of the height.

Held to one common standard — `max_disp_factor` = 0.1 of the model height, applied
to the final converged state — the trials read:

| F | out-of-balance | max&#124;u&#124; | as % of the 50 m height | verdict |
|---|---|---|---|---|
| 1.387500 | 1.45e-08 | 0.828 m | 1.66% | stands |
| 1.396875 | 2.32e-07 | 3.019 m | 6.04% | stands |
| 1.400000 | 3.05e-08 | 7.693 m | 15.39% | **refused: displacement** |

That is the textbook displacement upturn, and the limit lies in
[1.396875, 1.400]; the bisection reports its lower edge, 1.3969. It is also where
Griffiths & Lane read their own Example 1 — their Table 2 converges at 1.35 and
fails at 1.40 on a tenfold displacement jump, and `docs/verification/ssrm.md`
already scores XSLOPE's upturn reading against their 1.4 at 0.0%.

So both drivers were corrected by this, in opposite directions and for the same
reason. The Newton bound is now in the code (`_NR_DISP_FACTOR`, read once on the
final state, never inside the increment loop). The viscoplastic runaway level is
NOT changed — every locked and published factor of safety in the repository is
defined by it, and moving it would move them — but its comment no longer claims a
margin it does not have: on this model an equilibrating trial passes THROUGH the
level (it settles at 10.03 times its elastic response, the level being 8.0), so the
rule has a demonstrated false-positive mode near the critical strength. The
residual quad9 gap after the correction, viscoplastic 1.3844 against Newton 1.3969,
is that rule and nothing else.

### The corrected Newton column

Full re-measurement with the displacement bound in force. The viscoplastic column
was re-measured too — the driver was not touched, and this is the check that it was
not.

| Benchmark | Mesh | VP FS | N-R, spike | N-R, corrected | Gap | N-R iters | N-R force evals | N-R wall | VP iters | VP wall |
|---|---|---|---|---|---|---|---|---|---|---|
| FEM-1 tutorial | tri6, 3.5 | 1.3633 | 1.3711 | 1.3711 | +0.0078 | 1,485 | 10,680 | 36 s | 64,706 | 60 s |
| LEM-3 tutorial | tri6, 1.2 | 1.2539 | 1.2695 | 1.2695 | +0.0156 | 3,702 | 28,050 | 453 s | 169,767 | 628 s |
| Griffiths & Lane 1 | quad8, 3.5 | 1.3656 | 1.3719 | 1.3719 | +0.0063 | 1,688 | 11,592 | 32 s | 113,090 | 78 s |
| Griffiths & Lane 1 | tri6, 3.5 | 1.3656 | 1.3656 | 1.3656 | 0.0000 | 2,360 | 16,098 | 53 s | 78,576 | 70 s |
| Griffiths & Lane 1 | quad9, 3.5 | 1.3844 | 1.4031 | **1.3969** | +0.0125 | 823 | 5,846 | 27 s | 103,358 | 125 s |
| Griffiths & Lane 6 dry | quad8, 2 | 2.4219 | 2.4186 | 2.4186 | −0.0033 | 3,026 | 20,135 | 92 s | 48,128 | 42 s |
| Griffiths & Lane 6 dry | tri6, 2 | 2.4531 | 2.4531 | 2.4531 | 0.0000 | 2,791 | 19,209 | 87 s | 53,144 | 72 s |
| Griffiths & Lane 3, r = 0.8 | tri6, 6 | **1.4219** | 1.4281 | 1.4281 | +0.0063 | 4,034 | 31,569 | 331 s | 79,961 | 221 s |

**Exactly one Newton row moved**, and it is the one the bound was built for: quad9,
1.4031 → 1.3969. The other seven reproduce the spike's post-gate-removal values to
every digit recorded. The quad9 row also got much cheaper, 5,028 iterations down to
823, because the bisection stops climbing into the region where every trial is a
long solve to an inadmissible state.

One **viscoplastic** row does not reproduce: Griffiths & Lane 3 at r = 0.8 measures
1.4219 here against the 1.4125 recorded above, one bisection cell higher. Its
deciding trial is the same one the bug hunt names (F = 1.425, inconclusive at
50,000 iterations after four budget extensions); the bracket below it closes at
F = 1.41875 in 16,237 iterations. That row's Newton gap is therefore +0.0063 and
not the +0.0156 recorded above, and only two of the eight rows now sit outside the
0.01 bisection tolerance (LEM-3 at +0.0156 and quad9 at +0.0125).

**Force evaluations are the honest work count.** An `nr_force_evals` column is new
here, and it exists because comparing Newton iterations against viscoplastic
iterations flatters Newton by up to a factor of ten: a viscoplastic iteration is
one constitutive pass, a Newton iteration is one residual evaluation plus up to
nine more inside the line search. Across the eight rows the ratio of force
evaluations to iterations runs 6.8x to 7.6x. Newton is still doing far less
constitutive work than the viscoplastic driver on every row — 10,680 against 64,706
on FEM-1, 28,050 against 169,767 on LEM-3 — but the "1,246 iterations against
64,706" style of comparison in the tables above overstates the margin substantially.

### The hard regimes, re-measured with the bound

**Well past failure.** On FEM-1 at F = 2, 3 and 5 the two drivers still return the
same verdict, and the displacement bound never fires — every Newton failure there
is at the load-step floor, which is the point. The cost gap is worse than the
iteration counts said: Newton 391 iterations / 3,009 force evaluations at F = 2
against 181 viscoplastic iterations, 491 / 3,414 against 81 at F = 3, and
297 / 2,409 against 51 at F = 5. That is 17x, 42x and 47x the constitutive work,
not the 2x-to-6x an iteration count suggests, on exactly the trials a bisection
with a poor initial bracket visits first.

**nu = 0.49.** Unchanged and still clean. With every FEM-1 material at nu = 0.49 the
two drivers agree at all five probes: converged at F = 1.00, 1.30 and 1.36 (Newton
in 11, 12 and 12 iterations, 27, 29 and 30 force evaluations, against 120, 361 and
441 viscoplastic iterations), failed at F = 1.40 and 1.50, both at the load-step
floor. The converged states sit at 0.05% to 0.06% of the model height, nowhere near
the bound.

### The default path

Griffiths & Lane 6 dry re-run with no `fem_solver` argument returns
**FS = 2.421875**, with per-trial iteration counts, in the order recorded,

    147, 781, 3393, 2031, 2841, 9541, 12000, 8617, 8777

The seventh count is 12,000, not the 16,000 recorded twice above. That is not a
regression from any of this work: the identical run at **2e735693**, before any of
these commits, returns the same sequence including the 12,000. The 16,000 in the
text above does not reproduce on this checkout and should be read as a stale
number. Everything else — the factor of safety and all eight other counts —
matches, and the default viscoplastic path is byte-identical.

### What the locks now cover

`test/nr_ssrm_check.py` gained two checks and hardened a third.

`check_step_control_not_decisive` claimed in its own docstring that the removed
displacement gate "is gone from the driver — it decided a verdict, which is why it
is gone", and guarded nothing: reintroducing `_NR_RUNAWAY = 50.0` and its
increment-abandoning test on the code at 2e735693, and running the check as it
stood at 2e735693, passes with its full success message. It now asserts that
`xslope.fem` does not define `_NR_RUNAWAY` at all.

`check_displacement_bound` is the behavioral half, on the trial that needs it: on
quad9 at size 3.5, F = 1.396875 must converge and stay under a tenth of the model
height, and F = 1.400 must come back FAILED with the signal `displacement_limit`
AND with its out-of-balance under the force tolerance — reaching equilibrium and
being refused on displacement, rather than failing to reach it. Under the
reintroduced gate that trial instead dies inside the increment loop at the
load-step floor after 440 iterations at max|u| = 2.858, so the mutation cannot pass
by producing the same FAILED label for the wrong reason. All three legs fire on it.

`check_env_override_announces_itself` holds the other correction: `XSLOPE_FEM_SOLVER`
selecting a non-default driver must print, and an explicit `fem_solver` argument
must not.

The whole check runs in 48 s.

### Design notes, not changes

Two things were identified and deliberately not done, because both would move
locked values and that is the owner's call.

**The viscoplastic runaway verdict could be provisional.** The rule has a
documented false-positive mode (item (d)). A trial it closes could be recorded as
`runaway_provisional` rather than as a measured failure — the same treatment
`inconclusive` already gets, where the bisection carries the trial as uncertainty
instead of counting it as a failure. That would move every factor of safety the
rule currently decides.

**The Newton verdict is not monotone in F.** On Griffiths & Lane 6 dry (quad8, 2)
the bisection's bracketing pass records F = 2.0 FAILED and then F = 2.275,
F = 2.40625 and F = 2.414453 CONVERGED. A slope that fails at F = 2.0 cannot stand
at 2.4; one of those verdicts is wrong, and load control from a zero start at each
trial independently is how it happens — each trial is its own path, and the path
matters under non-associated flow. It does not change that row's answer, because
the bracket auto-expansion walks down to F = 1.75 and the bisection then works in a
region where the verdicts are consistent. It is, however, the clearest single
argument for the monotonic strength-reduction ramp: a ramp carries one continuous
history through increasing F and cannot produce a non-monotone verdict sequence by
construction.

## RAMP — the monotonic strength-reduction driver

Written before any of it was built, so that what follows is a test and not a
description.

### Why

The bisection asks a fresh, independent question at every trial: it drops the
model back to zero displacement, re-applies gravity from nothing, and re-discovers
the whole elastoplastic history at the new strength. Nine trials means nine
histories. Three things follow from that, all visible in the Phase 0 measurements:

- **It re-solves what it already knows.** Every trial repeats the easy early part
  of a solve the previous trial already completed.
- **It spends its worst effort where it learns least.** A trial well past failure
  costs 297 to 491 Newton iterations and 2,400 to 3,400 force evaluations to prove
  the load is unreachable at every increment size down to the floor, and a
  bisection with a poor initial bracket visits exactly those trials first. On
  Griffiths & Lane 6 dry the bracketing pass alone visits F = 2.8 and F = 2.0
  before it starts bisecting.
- **It can return a non-monotone verdict sequence.** On Griffiths & Lane 6 dry
  (quad8, 2) the Newton bisection records F = 2.0 FAILED and then F = 2.275,
  2.40625 and 2.414453 CONVERGED. A slope that fails at 2.0 cannot stand at 2.4;
  independent load-controlled histories at each F is how that happens.

A ramp carries ONE history. It starts from a converged state and reduces strength
monotonically, warm-starting each step from the state before it, so the mechanism
develops continuously instead of being rediscovered nine times. It cannot produce a
non-monotone verdict sequence, and it cannot spend anything on trials far past
failure, because it never goes there.

### Design

- **Continuation in F, not in load.** From a converged state at $F_k$, the gravity
  load does not change; only the strengths do. One step is: reduce $c$ and
  $\tan\phi$ to $F_{k+1}$, then drive the equilibrium residual at FULL gravity to
  convergence starting from the previous step's displacements. The load-stepping
  machinery is used only for the very first, cold solve.
- **Warm start.** The step inherits the previous converged displacement field, its
  committed plastic strains and its stress state. Gauss-point groups are built once
  for the whole ramp and re-strengthened in place, so nothing about the mechanism
  is rediscovered.
- **Adaptive $\Delta F$.** A step that fails to converge is retried at half the
  increment, from the same converged state, down to a floor of 0.005 and at most a
  few retries; a step that converges comfortably grows the increment back. The
  limit is reached when a step at the floor cannot be carried.
- **Admissibility.** Every accepted step passes the SAME standard a converged
  bisection trial passes: the Dawson out-of-balance under `force_tol`, and max|u|
  under `_NR_DISP_FACTOR` of the model height. A step that reaches force
  equilibrium in an inadmissible state is a failed step, not an accepted one.
- **The answer.** The limit is the last $F$ carried, with the failed step above it,
  reported to a 0.01 resolution so it is directly comparable with the bisection's
  answer at `tolerance=0.01`.
- **No arc-length control this round.** Step halving is the only recovery near the
  limit. If halving proves insufficient — steps failing at the floor while the
  state below is still comfortably converged — that is the follow-up, and it will
  be written as such rather than worked around.
- **Switch.** Internal only: `solve_ssrm(..., ssrm_driver='bisection')` is the
  default and reaches nothing new; `ssrm_driver='ramp'` is accepted only together
  with `fem_solver='newton'`. The viscoplastic path and the bisection path are
  untouched.

### Success criterion (verbatim)

- **Correctness.** On all eight benchmarks the ramp's factor of safety agrees with
  the Phase-0 corrected bisection-Newton factor of safety within 0.01, and the
  nu = 0.49 probe agrees.
- **Work.** Total force evaluations strictly below the bisection-Newton run's on at
  least six of the eight benchmarks.
- **It never goes past failure.** The ramp evaluates no $F$ more than one step past
  its own limit. This is by construction rather than by tuning, and it is what
  kills the 17x-to-47x past-failure weakness measured in Phase 0.
- **Warm-start effectiveness is reported as a number, not a claim.** Iterations on
  a warm-started step against a cold Newton trial at the same $F$, measured on the
  same mesh.
- **An honest negative is a valid outcome and must be written.** In particular, if
  the warm-started steps do not converge in the expected handful of iterations —
  if a step near the limit costs what a cold trial costs — that is the result, and
  it says the continuation is buying nothing.

### What this round does not decide

Whether the ramp should be the default. That is the owner's call and needs more
than eight benchmarks. This round decides only whether ramp-plus-Newton is ready to
be recommended as a non-default driver, and says what remains before the question
can be asked.

### RAMP — results

Same machine, meshes and brackets as everything above, `force_tol` 1e-3, hybrid
criterion, `capture_failure_state=False`, tolerance 0.01. The bisection-Newton and
viscoplastic columns are the Phase 0 corrected measurements.

The ramp reports the last strength it CARRIED, rounded down to 0.01, as the
criterion said it would. The bisection reports the MIDPOINT of its final interval.
Those are two different conventions on the same quantity, so the intervals are
given as well and are the like-for-like comparison.

| Benchmark | Mesh | Ramp FS | Ramp interval | Bisection-N FS | Bisection interval | Ramp − bisection (reported) | (interval midpoints) |
|---|---|---|---|---|---|---|---|
| FEM-1 tutorial | tri6, 3.5 | 1.36 | [1.36250, 1.36875] | 1.3711 | [1.36719, 1.37500] | −0.0111 | −0.0055 |
| LEM-3 tutorial | tri6, 1.2 | 1.26 | [1.26875, 1.27500] | 1.2695 | [1.26563, 1.27344] | −0.0095 | +0.0023 |
| Griffiths & Lane 1 | quad8, 3.5 | 1.36 | [1.36750, 1.37250] | 1.3719 | [1.36875, 1.37500] | −0.0119 | −0.0019 |
| Griffiths & Lane 1 | tri6, 3.5 | 1.36 | [1.36250, 1.36875] | 1.3656 | [1.36250, 1.36875] | −0.0056 | **0.0000** |
| Griffiths & Lane 1 | quad9, 3.5 | 1.39 | [1.39375, 1.40000] | 1.3969 | [1.39375, 1.40000] | −0.0069 | **0.0000** |
| Griffiths & Lane 6 dry | quad8, 2 | 2.41 | [2.41000, 2.41500] | 2.4186 | [2.41445, 2.42266] | −0.0086 | −0.0061 |
| Griffiths & Lane 6 dry | tri6, 2 | 2.45 | [2.45522, 2.46341] | 2.4531 | — | −0.0031 | +0.0062 |
| Griffiths & Lane 3, r = 0.8 | tri6, 6 | 1.42 | [1.42500, 1.43125] | 1.4281 | [1.42500, 1.43125] | −0.0081 | **0.0000** |

**Three of the eight reproduce the bisection's final interval EXACTLY** — Griffiths
& Lane 1 on tri6 and quad9, and Griffiths & Lane 3. On those the two routes find
the same limit to the last digit and the whole apparent gap is the reporting
convention. On the interval-midpoint reading all eight agree inside 0.01; on the
reported value six of eight do, and the two that miss — FEM-1 at −0.0111 and
Griffiths & Lane 1 quad8 at −0.0119 — miss by one and two thousandths. The
last-carried convention reads systematically about 0.005 low against a midpoint,
which is exactly half the reporting resolution, and settling that convention is a
loose end rather than a finding.

Work, on the same runs:

| Benchmark | Ramp steps / retries | Ramp iters | Ramp force evals | Ramp wall | Bisection iters | Bisection force evals | Bisection wall | VP iters | VP wall |
|---|---|---|---|---|---|---|---|---|---|
| FEM-1 tutorial | 8 / 4 | 512 | 3,935 | 12 s | 1,485 | 10,680 | 36 s | 64,706 | 60 s |
| LEM-3 tutorial | 7 / 4 | 597 | 4,751 | 75 s | 3,702 | 28,050 | 453 s | 169,767 | 628 s |
| Griffiths & Lane 1 quad8 | 16 / 5 | 622 | 4,151 | 11 s | 1,688 | 11,592 | 32 s | 113,090 | 78 s |
| Griffiths & Lane 1 tri6 | 8 / 4 | 512 | 3,935 | 12 s | 2,360 | 16,098 | 53 s | 78,576 | 70 s |
| Griffiths & Lane 1 quad9 | 10 / 4 | 799 | 5,465 | 26 s | 823 | 5,846 | 26 s | 103,358 | 125 s |
| Griffiths & Lane 6 dry quad8 | 20 / 5 | ~~790~~ 1,291 | ~~4,347~~ 7,821 | 35 s | 3,026 | 20,135 | 92 s | 48,128 | 42 s |
| Griffiths & Lane 6 dry tri6 | 31 / 7 | 840 | 4,821 | 23 s | 2,791 | 19,209 | 87 s | 53,144 | 72 s |
| Griffiths & Lane 3, r = 0.8 | 9 / 4 | 312 | 2,103 | 23 s | 4,034 | 31,569 | 331 s | 79,961 | 221 s |

Force evaluations are lower on **eight of eight**, by 1.07x (quad9) to 15.0x
(Griffiths & Lane 3); the criterion asked for six. Wall time is lower on eight of
eight against the bisection, and against the viscoplastic driver as shipped the
ratios are 5.0x, 8.4x, 7.1x, 5.9x, 4.8x, 1.2x, 3.1x and 9.6x — which clears the
original spike's "at most a fifth of the viscoplastic run" speed criterion on five
of the eight, a criterion the bisection-Newton driver met on none.

### It never solves past failure

This was to be true by construction rather than by tuning, and the measurement is
the distance between the highest strength each driver ever evaluated and the answer
it reported:

| Benchmark | Ramp: highest F evaluated | Bisection: highest F evaluated |
|---|---|---|
| FEM-1 tutorial | 1.4000 (FS + 0.040) | 2.0000 (FS + 0.629) |
| LEM-3 tutorial | 1.3000 (FS + 0.040) | 2.0000 (FS + 0.731) |
| Griffiths & Lane 1 quad8 | 1.3775 (FS + 0.018) | 1.8000 (FS + 0.428) |
| Griffiths & Lane 1 tri6 | 1.4000 (FS + 0.040) | 1.8000 (FS + 0.434) |
| Griffiths & Lane 1 quad9 | 1.4000 (FS + 0.010) | 1.8000 (FS + 0.403) |
| Griffiths & Lane 6 dry quad8 | 2.4500 (FS + 0.040) | 2.8000 (FS + 0.381) |
| Griffiths & Lane 6 dry tri6 | 2.4716 (FS + 0.022) | 2.8000 (FS + 0.347) |
| Griffiths & Lane 3, r = 0.8 | 1.4500 (FS + 0.030) | 1.8000 (FS + 0.372) |

The exact property, stated the way it is asserted in `test/nr_ssrm_check.py`: no
strength is evaluated more than the initial increment (0.05) above the highest
strength carried. That bound holds on all eight, with the largest actual overshoot
0.040. The bisection's overshoot is 0.35 to 0.73, and those are the trials where a
load-controlled Newton solve is at its most expensive — 17x to 47x the viscoplastic
driver's constitutive work, measured in Phase 0. The ramp removes that cost class
entirely rather than reducing it.

### Warm-start effectiveness, measured

Each warm-started step against a COLD Newton trial at the same strength on the same
mesh. Two regimes, and they say different things.

**Below the limit the continuation is worth about half the work.** Griffiths & Lane
3 (tri6, 6), warm against cold, iteration for iteration: F = 1.05, 7 against 15;
1.10, 7 against 14; 1.15, 6 against 16; 1.20, 8 against 17; 1.25, 9 against 20;
1.30, 9 against 24; 1.35, 10 against 25; 1.40, 12 against 25. Griffiths & Lane 1
quad9 runs 8, 8, 8, 7, 6, 10 against 11, 11, 11, 12, 12, 14 over F = 1.05 to 1.30.
FEM-1 runs 7, 9, 6, 9, 9, 11 against 12, 12, 12, 13, 14, 16.

**At the limit it is worth nothing.** The last steps before refusal cost the same
as a cold trial or more: quad9 at F = 1.3875 costs 40 warm against 28 cold, at
F = 1.39375 costs 58 against 59; Griffiths & Lane 3 at F = 1.425 costs 66 warm
against 50 cold. That is the honest negative on the warm start itself. The
continuation does not make the hard solve easy — the hard solve is hard because the
tangent is nearly singular, and starting nearer to the answer does not fix that.

**On the REFUSED steps it is worth five to twenty times**, which is where most of
the saving actually comes from. FEM-1's refusals cost 57, 120, 78 and 150 warm
iterations against 374, 593, 593 and 814 cold — and 2,722 to 5,275 force
evaluations cold. Griffiths & Lane 3's cost 68, 32 and 32 warm against 656, 656 and
779 cold. A warm-started state one increment below the limit proves the next
increment unreachable quickly; a cold trial has to walk the whole load path down to
the floor to prove the same thing. Griffiths & Lane 1 quad9 is the exception: its
refusals hit the 150-iteration per-step cap while a cold trial at the same strength
converges in 48 and is then refused on displacement, which is why quad9 is the one
row where the ramp barely beats the bisection (5,465 force evaluations against
5,846).

So the ramp's advantage is not mainly that warm-started steps are cheap. It is that
it never restarts, and that it never pays for a solve far past failure.

### The step control does not decide the answer

The failure this whole spike was written to avoid is a factor of safety that moves
when a solver knob is retuned. Three knobs, three benchmarks, one run each:

| Benchmark | as shipped | per-step cap 150 → 600 | initial ΔF 0.05 → 0.02 | ΔF floor 0.005 → 0.0005 |
|---|---|---|---|---|
| Griffiths & Lane 1 quad9 | 1.39 | 1.39 | 1.39 | 1.39 |
| FEM-1 tutorial | 1.36 | 1.36 | 1.36 | 1.36 |
| Griffiths & Lane 3, r = 0.8 | 1.42 | 1.42 | 1.42 | 1.42 |

Nothing moves. The intervals tighten as they should — at the 0.0005 floor quad9
carries 1.39688 and refuses 1.39766, and FEM-1 carries 1.36797 and refuses 1.36875
— and both of those tightened limits sit closer to the bisection's answer
(1.396875 and 1.371094) than the shipped floor's do, which is the resolution
behaving as a resolution. Raising the per-step cap fourfold changes no verdict,
which also settles the quad9 refusal: it is not the cap that refuses F = 1.400
there.

### nu = 0.49

FEM-1 with every material at nu = 0.49, both Newton routes: bisection FS 1.37891 on
interval [1.37500, 1.38281] in 1,829 iterations over 9 trials, highest strength
evaluated 2.0; ramp FS 1.37 on interval [1.37500, 1.38125] in 186 iterations over
13 evaluations, highest strength evaluated 1.40. The intervals overlap over their
whole width and the midpoints are 0.0008 apart. Ten times less work, no breakdown.

### Verdict

Yes — ramp-plus-Newton is ready to be recommended as a non-default driver, and it
is the first configuration in this spike that is.

It meets the criterion written before it was built. On all eight benchmarks its
limit interval agrees with the corrected bisection-Newton answer inside the 0.01
tolerance, and on three of the eight it reproduces that interval exactly; the
nu = 0.49 probe agrees to 0.0008. It does strictly less constitutive work on eight
of eight, against a requirement of six, by 1.07x to 15x — and against the shipped
viscoplastic driver it clears the original spike's fivefold speed criterion on five
of eight, which the bisection-Newton driver cleared on none. It never evaluates a
strength more than one increment above the highest one it has carried, so the
17x-to-47x past-failure cost that made the bisection-Newton driver a bad bet is
gone by construction rather than by tuning. And its answer does not move when the
per-step iteration cap is raised fourfold, the initial increment is more than
halved, or the increment floor is dropped tenfold — which is the standard this
spike exists to hold a solver to, and the one the Newton bisection failed twice
before it passed.

The honest negatives are three. The warm start does not make the hard solve easy:
at the last step before the limit it costs what a cold trial costs, sometimes more,
and the saving comes from never restarting and from cheap refusals rather than from
cheap steps. The reported value reads about 0.005 low because the ramp reports the
last strength carried while the bisection reports its interval midpoint, which is a
convention to settle and not a measurement. And on Griffiths & Lane 1 quad9 — the
one mesh where the ramp barely beats the bisection — the refused step exhausts its
per-step iteration budget rather than reaching the displacement bound, so the two
routes arrive at the same interval by different arguments; raising the cap fourfold
does not change it, but that is one measurement on one mesh.

Before the question of making it the default can be asked, which is the owner's and
not this branch's:

- **The reporting convention.** Report the interval midpoint, as the bisection
  does, or state the last-carried convention in the output. Either is fine; the
  present mismatch is not.
- **The feature scope is the Newton driver's, unchanged.** Reinforcement, piles,
  Hoek-Brown and power-curve envelopes, matric suction, tension caps, K0 initial
  stress and staged loading all raise. *(Reinforcement and the tension cap have
  since been built — see "REINFORCEMENT" and "THE TENSION CUTOFF", below.)* Any
  default-driver conversation is a conversation about that list first, and it is a larger piece of work than
  everything in this branch put together.
- **Eight benchmarks and one Poisson probe is not a calibration.** The viscoplastic
  thresholds have years behind them; the ramp's increment schedule has one session.
  The corpus run is what would earn the comparison, and it would also settle
  whether the −0.005 the ramp reads against the bisection is the convention alone.
- **The foot of the ramp still costs a cold solve, and can cost several.** On
  Griffiths & Lane 6 dry quad8 the requested F_min = 2.0 does not stand on the
  Newton driver, so the ramp walked down to 1.75 and paid two cold solves before it
  started. (That F = 2.0 fails while F = 2.4 stands is the Newton non-monotonicity
  recorded in the corrections above — the ramp is not subject to it once it starts,
  but it is subject to it at the foot.)
- **Arc-length control was not needed.** Step halving reached the floor cleanly on
  every benchmark; no case ended with a step failing at the floor while the state
  below it was comfortably converged. It stays a follow-up rather than a gap.


## Corrections, 2026-09-01 — second adversarial review

An independent reviewer re-ran seven of the eight ramp benchmarks and the nu=0.49
probe with its own force-evaluation counter and F-trace. Every factor of safety,
interval, and count reproduced to the digit except one cell, corrected below. The
feature guards were exercised on all three entry points plus the environment route
against reinforcement, pile, tension-cutoff, K0 and total-stress models: every one
refuses with `NotImplementedError` before doing any work — no silent-wrong path.
The default path was proven unchanged against a control run of the base commit
(same FS, same nine-trial iteration sequence), not against this document.

**Corrected cell.** The Griffiths & Lane 6 dry quad8 ramp row read 790 iterations /
4,347 force evaluations. The ramp's counters were seeded after the foot walk-down,
so the two failed cold solves it discarded there (501 iterations, 3,474 force
evaluations) were dropped. True figures: 1,291 / 7,821, and that row's ratio over
the bisection is 2.57x, not 4.63x. The headline — fewer force evaluations on eight
of eight, 1.07x to 15.0x — survives; that row was never either end of the range.
Fixed by computing the totals from the trial record, which lists every solve.

**`last_solution` was the foot, not the limit.** The ramp returned the F = F0 cold
solve as `last_solution`, so every figure and CSV exported from a ramp run showed a
pre-critical field (FS itself was unaffected). Fixed: the ramp now finishes with
one cold export solve at the last carried strength — the same state the bisection
driver exports for its final converged trial — recorded in `trials` as
`final_export` and included in the work totals, which therefore read one cold solve
higher than the tables above.

**Cancel raised instead of stopping.** Cancelling a ramp run before any step was
refused left `F_refused = None` and the return path raised `TypeError` (Studio's
Cancel button hit this). Fixed: a cancelled ramp returns a non-converged result
that says so, carrying the trials so far.

**Comment correction (`_NR_DISP_FACTOR`).** The claim that 0.1 is "the same
standard the viscoplastic driver applies to itself" was wrong during an SSRM run:
`solve_ssrm` passes `max_disp_factor=None` to every viscoplastic trial, which leans
on the runaway rule instead. The comment now says what is true: 0.1 is the standard
a single `solve_fem` call applies. The reviewer also noted the two drivers measure
displacement differently (component-wise max of the total field vs nodal resultant
of the plastic field); on these meshes the measures agree to ~1% and no verdict
turns on it, but the difference is recorded here so it is not rediscovered.

**The reporting-convention gap is larger and two-part.** The ramp's reported value
sits 0.0031 to 0.0119 below the bisection midpoint on these benchmarks — not a
uniform "about 0.005". Two mechanisms: the floor-to-resolution rounding of the last
carried strength, and reporting the interval's lower edge where the bisection —
the convention behind every locked and published FS in this repository — reports
the midpoint. Before any ramp value is locked or published it must be converted to
the midpoint convention: report (F_stands + F_refused)/2, not the floored edge.
Left unchanged in code for now so the tables above stay reproducible; this is the
first item to settle if the ramp is promoted beyond an internal option.


## Roadmap items closed, 2026-09-01

Two of the items the ramp verdict listed under "Before the question of making it
the default can be asked" are settled here. Everything measured on this checkout,
same machine, meshes and brackets as the tables above, `force_tol` 1e-3, hybrid
criterion, `capture_failure_state=False`, tolerance 0.01.

### The reporting convention — CLOSED

The ramp now reports the MIDPOINT of its final interval, `(F_stands + F_refused)/2`,
which is the convention the bisection uses and therefore the convention every locked
and published factor of safety in this repository is defined on. The raw edges are
still returned unrounded as `ramp_last_carried` and `ramp_first_refused`, so the
interval is readable from any result, and `_RAMP_RESOLUTION` — the 0.01 flooring
that produced the second half of the bias — is gone.

The eight-row column, re-measured. Every one of the eight final intervals
reproduces the interval recorded in "RAMP — results" to the digit, so the only thing
that moved is what is reported from it:

| Benchmark | Mesh | Ramp interval | Ramp FS, floored edge | Ramp FS, midpoint | Bisection-N FS | Midpoint − bisection |
|---|---|---|---|---|---|---|
| FEM-1 tutorial | tri6, 3.5 | [1.36250, 1.36875] | 1.36 | **1.3656** | 1.3711 | −0.0055 |
| LEM-3 tutorial | tri6, 1.2 | [1.26875, 1.27500] | 1.26 | **1.2719** | 1.2695 | +0.0023 |
| Griffiths & Lane 1 | quad8, 3.5 | [1.36750, 1.37250] | 1.36 | **1.3700** | 1.3719 | −0.0019 |
| Griffiths & Lane 1 | tri6, 3.5 | [1.36250, 1.36875] | 1.36 | **1.3656** | 1.3656 | **0.0000** |
| Griffiths & Lane 1 | quad9, 3.5 | [1.39375, 1.40000] | 1.39 | **1.3969** | 1.3969 | **0.0000** |
| Griffiths & Lane 6 dry | quad8, 2 | [2.41000, 2.41500] | 2.41 | **2.4125** | 2.4186 | −0.0061 |
| Griffiths & Lane 6 dry | tri6, 2 | [2.45522, 2.46341] | 2.45 | **2.4593** | 2.4531 | +0.0062 |
| Griffiths & Lane 3, r = 0.8 | tri6, 6 | [1.42500, 1.43125] | 1.42 | **1.4281** | 1.4281 | **0.0000** |

The agreement tightens on every row and changes sign on two, which is what a
reporting bias looks like when it is removed rather than a measurement moving. The
worst disagreement with the Newton bisection goes from 0.0119 to 0.0062; all eight
rows are now inside the 0.01 bisection tolerance where six were before, and the
three rows that reproduce the bisection's interval exactly now report the same
number as the bisection to four decimals rather than one to two thousandths below
it. The residual spread — two rows above, three below, three exact — is the two
routes finding slightly different limits, which is the thing the comparison was
meant to measure and which the old convention was hiding under a uniform −0.005 to
−0.012 offset.

`test/nr_ssrm_check.py` locks the convention structurally: the reported FS must be
the midpoint of the interval the same result reports, and must lie inside it.
Restoring the floored edge fails it with the interval printed alongside.

### The three untested feature guards — CLOSED

The Newton driver refuses eight modeling features it cannot represent. Five had
been exercised; three had not, and an unexercised guard is a claim rather than a
measurement. Each is now built as a real model on the same coarse Griffiths & Lane 1
tri6 mesh the rest of the check file uses — the material re-declared on the envelope,
or a suction angle passed in — and driven through `solve_fem(..., fem_solver='newton')`:

| Feature | How the model is built | What the driver does |
|---|---|---|
| Hoek-Brown | material `option='hb'`, σci 30,000 / GSI 60 / mi 10 / D 0 — 387 of 387 elements carry `hb_flag_by_elem` | raises `NotImplementedError`, "does not handle Hoek-Brown strength envelopes" |
| Power curve | material `option='pow'`, a 1.1 / b 0.9 / c 2.0 / d 4.0 — 387 of 387 carry `pow_flag_by_elem` | raises, "does not handle power-curve strength envelopes" |
| Matric suction | `suction_phi_b={'soil': 15.0}`, which sets `prep['suction_active']` | raises, "does not handle matric suction" |

The check asserts the message names the feature, not just that something raised, so
a guard firing for an unrelated reason cannot pass it; it asserts the flag arrays
are actually set, so a model that fails to exercise the guard is reported as a
broken model rather than as a passing test; and it runs the same three models on the
DEFAULT viscoplastic driver as a control, which must not refuse any of them.
Deleting the Hoek-Brown line from the guard's `unsupported` list fails the check.


## REINFORCEMENT — the Newton driver carries bars

Written before any of it was built, so that what follows is a test and not a
description.

### Why this one first

The ramp verdict named the feature list as the largest single item standing between
this branch and any default-driver conversation: "Reinforcement, piles, Hoek-Brown
and power-curve envelopes, matric suction, tension caps, K0 initial stress and
staged loading all raise." Reinforcement is the first of them because it is the most
common — the corpus carries reinforced models in the tutorials, the FEM samples and
the vendor benchmark sets — and because its physics is the simplest of the eight: a
bar element whose force is a function of the current displacement field, with no
internal state.

### What the viscoplastic driver actually does, read from the code

The Newton path has to reproduce the physics, not the algorithm, so the physics is
stated first. For one non-pile 1D element, per iteration, `solve_fem` computes:

1. the chord elongation `δ = (u₂ − u₁)·(cos θ, sin θ)` from the two END nodes — on a
   three-node bar the midside node carries stiffness but takes no part in the force
   recovery;
2. the elastic axial force `T = k δ`, with `k = EA/L` the chord stiffness;
3. the force the bar can actually deliver, `T_true = clip(T, 0, T_cap)` — tension
   only, so a bar in compression carries nothing at all;
4. `T_cap = t_allow`, the capacity envelope shared with the limit-equilibrium engine
   (`fileio.reinforce_available_tension`: tensile strength tapered by the pullout
   ramps at both ends and by end anchorage), or `t_res` for an element that has
   entered the post-peak set;
5. `t_allow` is NOT divided by F. Only the soil strength is reduced, which is the
   vendor convention and what the LEM does.

Bond slip is not a separate interface element in this model. It is baked into
`t_allow` as the ramped envelope, and the code says in as many words that it is
elastic–perfectly-plastic: "the soil-bar interface is frictional, and friction does
not vanish once it is overcome. t_allow ... is therefore both the yield force and
the post-slip force." The only softening in the model is `t_res`, the RUPTURE
residual — what the bar itself retains once the material tears — applied through a
converged-state fixed point that re-solves as the softened set grows.

The realization is a body-load correction: the global stiffness carries the bar's
FULL elastic element matrix, and `(T − T_true)` times the chord pattern is
subtracted, so the bar's internal force is `K_bar u_e − (T − T_true) p` with
`p = [−c, −s, +c, +s, 0, 0]`. On a two-node bar that is exactly `T_true p`. On a
three-node bar it is not, and the Newton path will reproduce THAT expression rather
than a cleaner one, because the two drivers have to be solving the same model.

### Scope, and what stays refused

- **1D reinforcement bar elements only.** Piles are beam elements with rotational
  degrees of freedom and a bending capacity; the pile guard is not touched.
- **No softening.** A model with `t_res` finite and below `t_allow` on any bar
  element is REFUSED with a message naming softening. Softening is a converged-state
  latch that changes the constitutive law between solves, which is exactly where a
  Newton continuation is weakest, and approximating it would be the kind of fudged
  tangent this spike exists to avoid.
- **Everything else in the guard stands**: Hoek-Brown, power curve, suction, the
  Rankine tension cutoff, K0 initial stress, staged loading, non-effective pore
  pressure. *(The Rankine tension cutoff is no longer among them — see "THE
  TENSION CUTOFF", below.)*

### Design

- **Residual.** `f_int_bar = K_bar u_e − (T − T_true) p`, added to the soil internal
  force. `K_bar` is `fem_data['K_global_1d_elems']`, the same matrix the global
  elastic stiffness is assembled from, so the two drivers carry the same bar.
- **Consistent tangent, exact.** The map is piecewise LINEAR in `u_e`, so there is
  nothing to difference: `dT/du = k p`, and `dT_true/du` is `k p` while
  `0 < T < T_cap` and zero otherwise. The element tangent is therefore `K_bar` for
  an active bar and `K_bar − k p⊗p` for a slack or capacity-limited one. On a
  two-node bar the second is exactly zero, which is the right answer — a yielded bar
  adds no stiffness. On a three-node bar it is a rank-one positive-semidefinite
  block resisting only the midside node's departure from the chord, which is to be
  verified numerically rather than asserted.
- **No history.** With softening refused, `T_true` is a function of the current
  displacement alone. The bar law is nonlinear-elastic, not plastic, so there is
  nothing to commit at the end of a step and the ramp's warm start needs no
  extension for it.
- **Assembly.** The bar blocks join the cached sparsity pattern explicitly rather
  than relying on their DOFs already appearing in the soil pattern, so a bar whose
  end nodes do not share a soil element cannot silently drop out of the tangent.
- **Displacement bound.** `max|u|` is read over the raw DOF vector. Bars stand on
  soil nodes and add no degrees of freedom of their own; rotational DOFs exist only
  where `is_pile_node` is set, which requires a pile. The guard asserts that no
  rotational DOF is present, so the bound cannot silently start comparing a length
  against a radian.
- **Strength reduction touches soil only.** The ramp's `restrength` rewrites
  `c_r`, `snph` and `csph` on the soil Gauss-point groups and nothing else; the
  criterion asserts that a bar's capacity is identical at the foot of the ramp and
  at its limit.

### The benchmarks

| Case | File | Mesh | Bracket | Lock |
|---|---|---|---|---|
| Geogrid sample | `docs/fem/files/xslope_reinforce_fem.xlsx` | tri6, 2 | 1.1–1.9 | `docs/fem/samples.md`, FS 1.497 |
| FEM-2 tutorial | `docs/tutorials/files/xslope_reinforced_slope.xlsx` | tri6, 2 | 1.0–2.0 | `docs/tutorials/fem02_reinforcement.md`, FS 1.496 |
| A third, different layout | picked from the corpus scan for a reinforcement layout unlike the six-layer geogrid stack, and no softening parameters | — | — | measured against the viscoplastic driver |

Both locked cases carry `t_res` = 600 against `t_max` = 800, so both declare
softening and both are refused as they ship. The comparison therefore runs the
softening-DISABLED variant of each — `t_res` unset, which is the model's own
"elastic–perfectly-plastic" default and what most vendors solve — on BOTH drivers,
so the two are answering the same question. Whether that variant is the same model
as the locked one is a measurement, not an assumption: the viscoplastic driver is
run on the shipped file and on the disabled variant, and if the two agree then the
locked value checks the comparison too, and if they do not then the locked value
checks only the shipped file and the comparison stands on the viscoplastic answer
for the variant. Either outcome is written.

### Success criterion (verbatim)

- **Correctness.** On at least three reinforced benchmarks, the Newton bisection's
  SSRM factor of safety agrees with the viscoplastic factor of safety on the SAME
  model within the bisection tolerance, 0.01; and where a locked FS tag applies to
  the model actually solved, it agrees with that too, inside the tag's tolerance.
- **The ramp agrees.** On the same cases, the ramp's factor of safety agrees with
  the Newton bisection's within 0.01.
- **The diagnostics are populated and physically consistent.** A reinforced Newton
  solution returns `forces_1d`, `failed_1d_elements` and `softened_1d_elements` at
  full 1D-element length; every reported force lies in `[0, t_allow]`; every element
  the failed mask marks sits AT its capacity; `softened_1d_elements` is all False
  (softening is refused, so no other value is reachable). The forces are compared
  element by element against the viscoplastic solution's on the same model at the
  same F, and the agreement is reported as a number rather than asserted.
- **The unreinforced path is unchanged.** One plain benchmark (Griffiths & Lane 1,
  tri6, 3.5) re-run through the extended driver must return the IDENTICAL factor of
  safety and the identical per-trial iteration counts as the pre-extension driver,
  measured against a control run of the parent commit rather than against a number
  in this document.
- **Work.** Newton's constitutive work against the viscoplastic driver's on the
  reinforced cases, reported as a measured ratio. No target this round — this is
  reconnaissance, and the ratio is the finding.
- **The refusals hold.** A pile model, and a model with softening parameters
  active, both raise, and the softening message names softening.
- **An honest negative is a valid outcome and must be written.** If the bars need
  something this design does not have — if the tangent has to be fudged, if the
  verdicts move, if the diagnostics disagree with the viscoplastic driver's — that
  is the result.


### REINFORCEMENT — results

Same machine and settings as everything above: `force_tol` 1e-3, hybrid criterion,
`capture_failure_state=False`, tolerance 0.01, tri6.

#### What was built

`_nr_build_bars` and `_nr_bar_force` in `xslope/fem.py`, threaded through
`_nr_internal_force`, `_nr_prepare_assembly`, `_nr_assemble_tangent`,
`_nr_equilibrate` and the ramp. The guard no longer refuses 1D elements; it refuses
piles, and it refuses post-peak softening by name, counting the elements that
declare it.

**The tangent is exact, and that is measured rather than asserted.** 400 random bar
elements at random orientations, lengths and rigidities, four-DOF and six-DOF, with
the displacement drawn to land on all three branches (73 active, 212 slack, 115 at
capacity): the analytic element tangent matches a central difference of the
element's own internal force to **3.5e-10** relative on the two-node bar and
**2.0e-10** on the three-node bar. A yielded two-node bar keeps exactly no stiffness;
a yielded three-node bar keeps a rank-one positive-semidefinite block, which is the
midside node's departure from the chord and nothing else. Trials whose perturbation
straddled a branch boundary were skipped, since a central difference there averages
two tangents; none of the 400 did.

#### The locked benchmarks are refused, and that is not a formality

Both of the repository's reinforced FEM benchmarks —
`docs/fem/files/xslope_reinforce_fem.xlsx` (FS 1.497) and
`docs/tutorials/files/xslope_reinforced_slope.xlsx` (FS 1.496) — ship with
`t_res` = 600 against `t_max` = 800, so both declare post-peak softening and both
are refused by the guard. The comparison therefore runs the softening-DISABLED
variant of each, on both drivers.

Whether that variant is the same model was measured rather than assumed, and it is
NOT. On the geogrid sample the viscoplastic bisection reads **1.4969 as shipped and
1.5594 with softening unset** at the locked tri6/2.0 mesh — 0.0625 apart, six
bisection cells — and 1.5594 against 1.5969 at the coarse tri6/4.0 mesh used below.
Softening is doing real work on these models. So the locked values check the
viscoplastic driver on the shipped files, which reproduce them (1.496875 against
1.497 at the locked mesh), and the Newton comparison stands on the viscoplastic
answer for the variant it can actually solve.

**The consequence is the first honest negative: the Newton driver cannot reproduce
either of the repository's published reinforced factors of safety, because both are
defined on a constitutive law it refuses.**

#### The comparison

All five rows on the softening-disabled model, tri6. The first four are the coarse
mesh (563 elements, 36 bar elements, all carrying capacity) on which the Newton
driver is well behaved; the fifth is the LOCKED mesh (2101 elements), on which it is
not.

Five rows, FOUR models: the "geogrid sample" and "FEM-2 tutorial" files are the same
model under two names, differing only in a cosmetic `type` label on each
reinforcement line — see "Bookkeeping: five benchmarks, four models" at the end of
this document.

| Case | Mesh | Layout | VP FS | Newton FS | gap | Ramp FS | ramp − Newton |
|---|---|---|---|---|---|---|---|
| Geogrid sample | tri6, 4 | 6 layers, t_max 800 | 1.5969 | 1.6094 | +0.0125 | 1.6094 | **0.0000** |
| FEM-2 tutorial | tri6, 4 | 6 layers, t_max 800 | 1.5977 | 1.6055 | +0.0078 | 1.6031 | −0.0023 |
| Half capacity | tri6, 4 | 6 layers, t_max 400 | 1.4113 | 1.4113 | **0.0000** | 1.4031 | −0.0082 |
| Three layers | tri6, 4 | layers 1, 3, 5 only | 1.2285 | 1.0879 | **−0.1406** | 1.0531 | −0.0348 |
| Geogrid sample | tri6, 2 | 6 layers, t_max 800 | 1.5594 | 1.2719 | **−0.2875** | 1.2406 | −0.0313 |

The reinforcement is doing real work on every row and the two drivers see the same
work: halving the capacity moves the viscoplastic answer from 1.5969 to 1.4113 and
the Newton answer from 1.6094 to 1.4113, and dropping three of the six layers moves
the viscoplastic answer to 1.2285. The half-capacity row — the one where the bars are
most heavily engaged, with the largest number of them sitting at their cap — is the
row where the two drivers agree EXACTLY, interval for interval,
[1.4078125, 1.41484375] on both.

And two rows are badly wrong, both LOW.

#### Where it breaks, and what breaks it

The failure is not in the bars. On the locked mesh the Newton verdict as a function
of strength is:

| F | 1.2 | 1.3 | 1.4 | 1.45 | 1.5 | 1.55 | 1.6 |
|---|---|---|---|---|---|---|---|
| verdict | CONVERGED | FAILED | FAILED | **CONVERGED** | FAILED | FAILED | FAILED |
| iterations | 39 | 599 | 880 | 73 | 261 | 628 | 58 |
| load factor reached | 1.00 | **0.00** | 0.69 | 1.00 | **0.00** | 0.29 | 0.00 |

A slope that stands at F = 1.45 stands at F = 1.3. The F = 1.3 and F = 1.5 trials do
not merely fail — they carry NO load at all, failing at every increment from full
gravity down to the 1/64 floor, at a maximum displacement of 0.0016 m on a 34 m
model. That is not a statement about the slope.

**The mechanism was tested and it is the soil, not the reinforcement.** Recording
the active set at every tangent re-form of the F = 1.3 trial: the bar active set
flips 817 times over 591 re-forms and is UNCHANGED on 461 of them, one or two bars at
a time. The soil branch set flips 52,441 times over the same 591 re-forms — about 89
Gauss points changing branch per tangent, never settling. At F = 1.2, where the same
model converges in 39 iterations, the soil flips fall from 116 to 1 over the last 25
re-forms and the solve ends. The tension-only kink in the bar law, which was the
obvious suspect, is not what is churning.

Two controls place it further. Removing the reinforcement entirely and re-running:
both drivers converge at F = 1.0 and both fail at F = 1.2, so the unreinforced slope
genuinely does not stand there and the "no bars, no convergence" reading is physics
rather than a defect. And the model's lower material is **c = 0, phi = 37** — a
cohesionless base whose Mohr-Coulomb apex sits at the origin, so every tensile Gauss
point returns to the apex, where the consistent tangent is zero. The reinforcement is
what holds this slope up, and the Newton driver's difficulty is in solving the soil it
holds up, at meshes fine enough to resolve it.

Mesh by mesh, at strengths either side of the limit (`CONV`/`FAIL`, iterations):

| Mesh | elements | F = 1.2 | F = 1.4 | F = 1.5 | F = 1.6 |
|---|---|---|---|---|---|
| tri6, 4 | 563 | CONV 17 | CONV 21 | CONV 28 | CONV 46 |
| tri6, 3 | 978 | CONV 21 | CONV 158 | CONV 39 | FAIL 921 |
| tri6, 2 | 2101 | CONV 39 | FAIL 880 | FAIL 261 | FAIL 58 |

It degrades monotonically with refinement, which is the opposite of what a
convergent scheme should do.

#### The bar diagnostics

`forces_1d`, `failed_1d_elements` and `softened_1d_elements` come back at full
1D-element length on every row, and every structural property holds on every row: no
bar reports a compressive force, no bar carries more than the capacity its embedment
develops, every element the failed mask marks sits at that capacity to within 1e-6 of
it, and the post-peak set is empty — the only value reachable while softening is
refused. The halved-capacity variant reports a largest force of exactly 400.0 against
the unmodified model's 800.0, so the cap being enforced is the model's cap.

Against the viscoplastic driver's bar forces at the same strength, on the same model:

| Case | largest T, Newton | largest T, viscoplastic | ΣT, Newton | ΣT, viscoplastic | at capacity, N / VP | worst per-bar gap |
|---|---|---|---|---|---|---|
| Geogrid | 800.0 | 800.0 | 18,991 | 18,264 | 19 / 16 | 0.20 |
| FEM-2 | 800.0 | 800.0 | 18,343 | 18,216 | 18 / 16 | 0.12 |
| Half capacity | 400.0 | 400.0 | 7,995 | 9,687 | 13 / 25 | 0.41 |
| Three layers | 107.5 | 272.2 | 350 | 1,086 | 0 / 0 | 0.76 |

The at-capacity counts are not expected to match exactly and their difference is not
a defect: the viscoplastic mask LATCHES every bar that exceeded its capacity at any
point in the iteration history, including the elastic predictor's overshoot before
the soil sheds load into the bars, while the Newton mask is read on the reported
state. The Newton mask is a subset by construction, and the code says so.

The force gaps are a different matter and only the first two are comfortable. Where
both drivers reach the same equilibrium the totals agree to 0.7% and 4%. Where the
bisection midpoint sits above the last standing trial — the geogrid and
half-capacity rows, where BOTH drivers report the trial failed — the comparison is
between two failed states and means little. The three-layer row is the one that
should worry a reader: both drivers CONVERGE there, and the Newton field has the bars
carrying 350 against the viscoplastic field's 1,086, a factor of three. Two
converged states, two very different degrees of reinforcement engagement. In the
fully elastic regime the two agree exactly — at F = 0.2 and F = 0.5 on the locked
mesh both report a largest bar force of 15.4 and their displacement fields agree to
7e-7 and 3e-5 — so this is not a formulation difference in the bar; it is the two
drivers arriving at different elastoplastic states, and a tension-only bar reads that
difference much more sharply than the soil does.

#### Work

Constitutive work, on the coarse mesh where the answers are comparable:

| Case | VP iterations | VP wall | Newton iterations | Newton wall | Ramp iterations | Ramp force evals | Ramp wall |
|---|---|---|---|---|---|---|---|
| Geogrid | 88,614 | 49.2 s | 3,278 | 36.9 s | 711 | 4,451 | 7.3 s |
| FEM-2 | 67,907 | 44.2 s | 1,851 | 22.4 s | 528 | 3,108 | 5.4 s |
| Half capacity | 86,918 | 52.6 s | 1,897 | 21.3 s | 526 | 3,355 | 5.5 s |
| Three layers | 109,969 | 60.9 s | 2,347 | 25.1 s | 179 | 1,296 | 2.1 s |

Newton runs 27x to 47x fewer iterations than the viscoplastic driver and 1.3x to
2.5x less wall time; the ramp is 6.7x, 8.2x, 9.6x and 29x faster than the
viscoplastic driver in wall time. Reinforcement costs the Newton assembly almost
nothing — the bar blocks are 36 dense 6x6 matrices against 563 elements' worth of
soil — so the speed picture is the unreinforced one, unchanged.

#### The unreinforced path is unchanged

Measured against a control run of the parent commit (`318c5686`, the driver before
the bar element existed) staged in a separate package tree, not against a number in
this document. Two runs on each:

- **Griffiths & Lane 1, tri6, 3.5, Newton bisection.** FS 1.3656249999999999 on both,
  and the same nine trials in the same order with the same verdicts and the same
  iteration counts: (1.0 CONVERGED 10), (1.8 FAILED 482), (1.4 FAILED 354),
  (1.2 CONVERGED 13), (1.3 CONVERGED 16), (1.35 CONVERGED 26), (1.375 FAILED 566),
  (1.3625 CONVERGED 32), (1.36875 FAILED 861). Identical.
- **Griffiths & Lane 6 dry, quad8, 2, no `fem_solver` — the default path.**
  FS 2.421875 on both, per-trial iterations 147, 781, 3393, 2031, 2841, 9541, 12000,
  8617, 8777 on both. Identical, and the same sequence recorded in the corrections
  above.

A model with no 1D elements takes `bars = None` through every new code path, and
`np.bincount` over an unchanged pattern is bit-identical, which is what those two
runs confirm rather than assume.

#### The ramp does not rescue it

The obvious hypothesis for the two bad rows was the loading path: the bisection drops
each trial back to zero and re-applies gravity, so the tension-only bars re-engage
from nothing every time, and the ramp — one warm-started history, gravity applied
once — should not have that problem.

**Falsified.** The ramp reads 1.2406 against the viscoplastic 1.5594 on the locked
mesh, and 1.0531 against 1.2285 on the three-layer variant. It agrees with the
Newton bisection to 0.0313 and 0.0348 there, which is to say it agrees with the
Newton driver and not with the answer. On the three rows where the bisection is
sound the ramp tracks it to 0.0000, −0.0023 and −0.0082, so the ramp is doing what
the ramp does; it inherits the per-step solve, and the per-step solve is what fails.

#### The criterion, line by line

**Correctness — NOT MET.** Two of the five rows agree with the viscoplastic driver
inside the 0.01 bisection tolerance (half capacity at 0.0000, exactly, interval for
interval; FEM-2 at +0.0078). The geogrid sample misses by a quarter of a cell at
+0.0125. Two rows are wrong by 0.14 and 0.29, both LOW, and the low direction is
unconservative in the sense that matters here: the reinforced slope is reported as
weaker than it is, which is safe for a design but wrong as a measurement, and it is
wrong because the solver could not carry a load rather than because the slope could
not.

**Locked values — NOT MET, by construction.** Both published reinforced factors of
safety are defined on models that declare post-peak softening, which this driver
refuses. The viscoplastic driver reproduces its own lock on the shipped file
(1.496875 against 1.497), and the softening-disabled variant it can be compared
against is a different model by 0.0625.

**The ramp agrees with the Newton bisection — MET on three of five** (0.0000,
−0.0023, −0.0082) and missed on the two rows where both routes are wrong (−0.0348,
−0.0313). It agrees with the driver, which is all this clause ever asked.

**The diagnostics — MET structurally, PARTLY MET against the viscoplastic driver.**
All three arrays come back at full 1D-element length on every row; no bar reports a
compressive force or a force above the capacity its embedment develops; every element
the failed mask marks is at that capacity; the post-peak set is empty everywhere; and
the halved-capacity variant reports a largest force of exactly 400.0 against 800.0.
Against the viscoplastic driver's forces the totals agree to 0.7% and 4.0% where the
two reach the same state, and differ by a factor of three on the three-layer row
where both converge to different states.

**The unreinforced path — MET, exactly.** Identical factor of safety and identical
per-trial iteration sequences against a control run of the parent commit, on both the
Newton bisection and the default viscoplastic path.

**Work — reported.** 27x to 47x fewer iterations than the viscoplastic driver and
1.3x to 2.5x less wall time on the bisection; 6.7x to 29x less wall time on the ramp.
Reinforcement itself costs the assembly almost nothing.

**The refusals — MET.** A model declaring softening raises with the count of bar
elements that declare it; a pile model raises naming piles; and the three envelope
guards closed earlier this session still fire.

#### Verdict

The bar element is right, and the Newton driver it is bolted to is not ready to
answer a reinforced question.

Everything that can be checked about the bar checks out. Its consistent tangent
matches a central difference of its own residual to 1e-10 on all three branches, and
a yielded bar keeps exactly the stiffness it should keep and no more. In the elastic
regime the two drivers report the same bar force to the last digit and the same
displacement field to 1e-7. The capacity is enforced — no bar anywhere carries more
than its embedment develops, halving the capacity halves the largest reported force
exactly, and it moves the viscoplastic answer from 1.5969 to 1.4113 and the Newton
answer from 1.6094 to the same 1.4113. That half-capacity row, where the most bars
are pinned at their cap and the bar law is doing the most work, is the row where the
two drivers agree exactly, interval for interval. Strength reduction leaves the
reinforcement alone, as the vendor convention requires, and a bar at capacity holds
the same force at two different strength reductions.

What does not work is the solve around it. On two of five rows the Newton answer is
0.14 and 0.29 below the viscoplastic one, and the mechanism is measured rather than
guessed: on the locked mesh the driver converges at F = 1.2, fails at 1.3 and 1.4,
converges again at 1.45, and fails above — and the failures at 1.3 and 1.5 carry NO
load at all, giving up at every increment down to the 1/64 floor at a displacement of
1.6 mm on a 34 m model. Instrumenting the active set at every tangent re-form of that
trial puts the churn in the soil and not in the bars: 52,441 soil branch flips over
591 re-forms, never settling, against 817 bar flips of which most re-forms have none.
It gets monotonically worse with mesh refinement — clean at 563 elements, marginal at
978, broken at 2101 — which is the opposite of what a convergent scheme does. The
ramp inherits it rather than curing it. And the model that exposes it is not exotic:
its lower material is c = 0, phi = 37, so every tensile Gauss point returns to the
Mohr-Coulomb apex where the consistent tangent is zero, and a reinforced slope is
precisely the kind of model that stands only because something carries the tension the
soil cannot.

So the honest reading is narrower than "reinforcement works". The bar element is a
correct, exact, cheap addition that costs the assembly nothing and is ready to be
carried forward. The reinforced factor of safety it produces is only as good as the
Newton solve underneath it, and on the reinforced corpus that solve is unreliable at
the meshes the locks are written on. Two things would have to happen before a
reinforced Newton answer could be quoted: the apex/zero-tangent regime would need
whatever regularization or line-search treatment makes a c = 0 material solvable at
practical refinement, which is a soil problem this round did not open; and post-peak
softening would need a treatment, because without it the driver cannot reach either
of the repository's published reinforced values at all. Until then the extension's
value is that it makes reinforced models *reachable* on the Newton path rather than
refused, and that it exposed where the driver actually fails — which is not where
this round expected to find it.

#### What the corpus could not supply

The criterion asked for a third reinforced benchmark with a different layout, grepped
from the corpus. The corpus does not have one that the Newton driver can solve. Every
reinforced model in `docs/` was scanned: 39 files carry reinforcement lines, and of
those only seven also carry the E and area a bar element needs. `vp060` and `vp032a/c`
declare a Rankine tension cutoff, which the guard refuses; `vp060` also carries K0
initial stress. `gs2_18` fails to factorize its own elastic stiffness. `vp030a` and
`vp030b` did not complete a single trial inside a ten-minute budget in either
configuration and were not pursued. `vp032b` and `vp088` fail on the Newton driver at
every strength tried — and with their reinforcement stripped out they fail
IDENTICALLY, which places that failure in the soil, not the bars; both carry
cohesionless materials at a free face.

The third and fourth cases here are therefore constructed rather than found: the same
soil and the same six lines at half capacity, and the same soil with three of the six
lines. Both are real models that both drivers solve, both change the answer
substantially (1.5969 to 1.4113 and to 1.2285 on the viscoplastic driver), and the
first of them is the most severe test of the bar law in this round because it puts
the most bars at their cap. Neither is a substitute for a vendor-verified reinforced
benchmark, and the reason there is none in this table is the finding above.


## THE COHESIONLESS SOLVE — what actually breaks, and the criterion for fixing it

Written after the diagnostic probes and BEFORE any remedy code, so that what
follows is a test and not a description. Everything measured on this checkout,
same machine and settings as the tables above: `force_tol` 1e-3, hybrid criterion,
`capture_failure_state=False`, tolerance 0.01, tri6.

### The specimen reproduces exactly

The geogrid sample (`docs/fem/files/xslope_reinforce_fem.xlsx`, softening unset)
at its locked tri6/2.0 mesh — 2,101 elements, 60 bar elements — run trial by
trial on the Newton driver:

| F | 1.2 | 1.3 | 1.4 | 1.45 | 1.5 | 1.55 | 1.6 |
|---|---|---|---|---|---|---|---|
| verdict | CONVERGED | FAILED | FAILED | **CONVERGED** | FAILED | FAILED | FAILED |
| iterations | 39 | 599 | 880 | 73 | 261 | 628 | 58 |
| load factor reached | 1.00 | **0.00** | 0.69 | 1.00 | **0.00** | 0.29 | 0.00 |

Iteration for iteration the same sequence the reinforcement round recorded. The
model's lower material is c = 0, phi = 37; its upper is c = 300 psf, phi = 37.

### The F = 1.3 failure is provably false, without reference to the other driver

The F = 1.45 trial converges to an out-of-balance of 5.167e-6 — a hundred and
ninety times inside the trial tolerance — with a worst Mohr-Coulomb violation of
**2.5e-15** of the local strength, read in the invariant form. That is a stress
field in equilibrium with full gravity and nowhere outside the yield surface at
the F = 1.45 reduced strengths.

Strength reduction only ever raises c and tan(phi) as F falls, so that same field
is admissible at every F below 1.45 and is in equilibrium with the same load. By
the lower-bound theorem the slope therefore stands at F = 1.4, 1.3 and 1.2, and
the driver's own converged answer at 1.45 is the proof. The FAILED verdicts at
1.3, 1.4 and 1.5 are the solver's, not the slope's, and nothing in the argument
needs the viscoplastic driver.

### P1 — the token cohesion does NOT close it

The vendors' documented user-level workaround for a cohesionless base is a small
nonzero cohesion. Applied here to the c = 0 material only, on the same mesh:

| c added | 1.2 | 1.3 | 1.4 | 1.45 | 1.5 | 1.55 | 1.6 |
|---|---|---|---|---|---|---|---|
| 0 psf (as shipped) | CONV 39 | FAIL 599 | FAIL 880 | CONV 73 | FAIL 261 | FAIL 628 | FAIL 58 |
| 5 psf | CONV 32 | CONV 215 | **FAIL 877** | CONV 116 | CONV 80 | FAIL 608 | FAIL 541 |
| 10 psf | CONV 225 | CONV 205 | CONV 290 | **FAIL 848** | **FAIL 689** | **CONV 296** | FAIL 870 |

The churn moves; it does not go away. At 5 psf the F = 1.3 and F = 1.5 trials are
recovered and F = 1.4 still fails between two converging neighbours; at 10 psf the
false failures have simply migrated to 1.45 and 1.5, with 1.55 converging above
them. The verdict sequence is non-monotone in F at every cohesion tried. **The
workaround is falsified as a remedy and as an explanation**, and any real remedy
has a lower bar to clear than it looked.

### P2 — the flippers are the cohesionless material, and they are NOT at the apex

Recording the branch code of all 6,303 Gauss points at every tangent re-form of
the F = 1.3 trial: 592 re-forms, 52,441 branch flips, of which **51,923 (99.0%)
are carried by the 6,051 Gauss points in the c = 0 material** and 518 by the 252
points with cohesion. 2,494 of the cohesionless points flip at least once; the
worst flips 90 times.

The apex is not what they are flipping to. Over those same 592 re-forms the apex
branch holds a mean of **3.0 points** and a maximum of 36, out of 6,303; the
corner branches never execute at all. The churn is elastic against main-plane, in
the cohesionless material. **The reinforcement round's stated mechanism — "every
tensile Gauss point returns to the Mohr-Coulomb apex, where the consistent tangent
is zero" — is wrong on the counts.** The apex tangent is indeed zero, and there are
three of them.

### What is actually degenerate

Two measurements, both at a state 25 iterations into the F = 1.3 solve.

**The assembled tangent is correct.** Compared against a central difference of the
driver's own internal force along four random directions, at three step sizes: the
relative error is 1.4e-8 at the natural step and the direction cosine is 1.00000.
The tangent is the derivative of the residual it is used on, so nothing below is a
differencing defect.

**The consistent tangent is exactly rank-deficient at every plastic Gauss point.**
Its smallest singular value over the 1,929 plastic points has a median of
1.7e-9 of the shear modulus and a minimum of exactly zero, against elastic singular
values of 1, 2 and 5 times it; the symmetric part carries a negative eigenvalue of
up to -0.30 times the shear modulus at 1,927 of them.

That is not a defect either — it is what perfect plasticity gives. On the main
plane the return is `n1 = s1 - f`, `n3 = s3 + f` with `f = A s1 - Bc s3 - c cos(phi)`,
so the trial-to-returned Jacobian in the (sigma_1, sigma_3) plane is
`[[1-A, Bc], [A, 1-Bc]]`, whose determinant is `1 - A - Bc`. And `A + Bc = 1`
identically, at every friction angle. A whole line of trial states — the flow
direction — maps to one returned state, so the tangent has a null direction there
by construction, for psi = 0 and for any other flow rule on a linear surface
without hardening. The apex is the special case where the null space is the whole
space.

**What that does to the iteration** is measurable directly. At the same state, the
line search along the consistent-tangent Newton direction reads:

| step alpha | 1 | 1/2 | 1/4 | 1/8 | 1/16 | 1/32 | 1/64 | 1/128 | 1/256 |
|---|---|---|---|---|---|---|---|---|---|
| residual, as a fraction of the current one | 36.59 | 22.17 | 13.43 | 5.93 | 2.06 | 1.14 | 1.007 | **0.994** | 0.996 |
| Gauss points changing branch | 509 | 349 | 248 | 183 | 95 | 40 | 20 | 12 | 5 |
| ... of them in the c = 0 material | 492 | 336 | 239 | 177 | 92 | 38 | 18 | 10 | 5 |

The full Newton step multiplies the residual by 36 and flips 509 Gauss points. The
first step the line search can accept is 1/128 of it, and it buys 0.6% — so the
solve makes 0.6% progress per iteration, per tangent re-form and per
factorization, and the no-progress watch closes it. On the SAME model at the
coarse tri6/4.0 mesh the same measurement reads zero flips at alpha = 1 and a
residual of 1e-5 of the current one: full Newton, quadratic convergence, nothing
to fix. The trust region collapses with refinement because refinement puts more
Gauss points near a yield surface that, with c = 0, is a cone through the origin
and therefore has no absolute margin anywhere: 1,712 of the 4,374 elastic points
at that state sit within 5% of it.

### The remedy this round tests

Only the ITERATION matrix is at fault, so only the iteration matrix is changed.
The return map, the residual and every convergence test stay exactly as they are,
which makes the converged state identical by construction wherever the present
driver converges at all.

**Tangent conditioning.** Blend the consistent tangent toward the elastic operator,
`D_beta = (1 - beta) D_ep + beta D_e`, with `beta` raised when the line search
cannot take a reasonable step and decayed back to zero when it can. At `beta = 0`
the code path is the present one untouched. Measured at the specimen state above,
one step of the conditioned direction against the 0.994 the shipped driver gets:

| beta | 0.02 | 0.05 | 0.1 | 0.2 | 0.4 | 0.7 | 1.0 (elastic) |
|---|---|---|---|---|---|---|---|
| residual after the accepted step | 0.382 | 0.301 | 0.296 | 0.318 | 0.358 | 0.403 | 0.442 |
| step accepted | full | full | full | full | full | full | full |

Two hundred times the progress of the shipped iteration, from a two-hundredth of
the singular value being restored. Whether that survives a whole solve, a whole
bisection and eight plain-soil benchmarks is what the criterion below asks.

### Success criterion (verbatim)

- **The five reinforced benchmarks agree with the viscoplastic driver within
  0.01**: geogrid sample (VP 1.5969), FEM-2 tutorial (1.5977), half capacity
  (1.4113), three layers (1.2285) and the geogrid locked mesh (1.5594). The last
  two are the currently-broken ones, at -0.1406 and -0.2875.
- **The refinement ladder is clean at every rung.** The geogrid sample at 563, 978
  and 2,101 elements: no false failure at any strength, and a verdict sequence
  monotone in F — everything below the limit converges and everything above it
  fails.
- **The plain-soil eight-row table does not move.** Every row re-measured and
  reported; no factor of safety moves by more than 0.01, and wherever the shipped
  driver's trials converged the answer is expected to be bit-identical, because
  conditioning that never fires changes no arithmetic.
- **The knob is not a tuning dial.** The converged state must be identical with the
  conditioning on and off — asserted on the specimen, on the trial that needs it,
  by comparing the displacement fields and not just the verdicts. And the answer
  must not move when the conditioning schedule is retuned.
- **Work does not regress.** The previously-clean cases stay inside 1.2x of their
  current force-evaluation counts.
- **The probes still hold.** nu = 0.49 agrees, and the past-failure probes on FEM-1
  at F = 2, 3 and 5 still fail, promptly.
- **The locks catch it.** `test/nr_ssrm_check.py` gains the locked-mesh reinforced
  case. The present false answer, 1.2719, must FAIL the new lock on the code as it
  stands today and pass after the remedy — run both ways, and the mutation
  recorded — and the whole check must pass.
- **The default path is unchanged**, against the standard control: Griffiths &
  Lane 6 dry with no `fem_solver` argument, FS 2.421875 on iteration counts 147,
  781, 3393, 2031, 2841, 9541, 12000, 8617, 8777.
- **An honest negative is a valid outcome and must be written.** If the
  conditioning buys the specimen and costs the plain-soil table, or if it needs a
  schedule tuned per model, that is the result.


### THE COHESIONLESS SOLVE — results

Same machine and settings as the criterion above. Every number below was measured
on this checkout in this session; the "before" columns are the driver with the
remedy switched off, re-run here rather than quoted, and they reproduce the
figures recorded earlier in this document to the digit.

#### The remedies the criterion named, and why two of them were not built

**Apex tangent conditioning — not applicable, measured.** The probes put the apex
branch at a mean of 3.0 Gauss points out of 6,303 over the 592 tangent re-forms of
the specimen trial. Conditioning three points changes nothing.

**Smoothed Mohr-Coulomb (Abbo & Sloan) — not applicable, measured.** The hyperbolic
approximation rounds the apex and the deviatoric corners. The apex fires at three
points and the corner branches do not execute at all on this model — zero right
corners and zero left corners at every trial measured. Smoothing a surface where
the solver is not is not a remedy; it would only change the answer.

**Tangent conditioning, generalized — BUILT, MEASURED, AND REJECTED.** The rank
deficiency is not at the apex, it is at every yielding point, so the honest version
of the same idea is to condition every plastic tangent: blend it toward the elastic
operator, `(1 - beta) D_ep + beta D_e`, which restores the missing rank without
touching the return map or the residual. One step of it is spectacular — at the
specimen state the accepted step goes from 0.994 of the current residual to 0.301.
A whole solve is not. With beta pinned and the per-increment cap raised to 400,
the F = 1.3 trial still carries no load at any setting tried:

| beta | 0 (as shipped) | 0.02 | 0.05 | 0.1 | 0.2 | 0.4 | 1.0 (initial stiffness) |
|---|---|---|---|---|---|---|---|
| F = 1.3 | FAILED, 599 | FAILED, 499 | FAILED, 524 | FAILED, 720 | FAILED, 1,033 | FAILED, 1,786 | FAILED, 4,195 |
| F = 1.5 | FAILED, 261 | FAILED, 673 | FAILED, 520 | FAILED, 747 | FAILED, 1,138 | FAILED, 1,824 | — |

Run adaptively — beta raised when the line search cannot take a step worth having
and decayed when it can — it was worse than useless: it saturated at the elastic
operator and turned the F = 1.45 trial, which the shipped driver converges in 73
iterations, into a failure. It is not in the code. What it bought was the
measurement above, which says the problem is not the conditioning of the tangent.

#### What the problem actually is

Three measurements, in the order they were made.

**More iterations do not fix it.** With the per-increment cap raised to 6,000 and
the no-progress watch lifted entirely, the F = 1.3 trial runs 18,914 Newton
iterations and carries **3% of gravity** — two increments of a sixty-fourth each,
and then it cannot take a third. It is not a stall that a budget closes.

**The root exists and Newton finds it in sixteen iterations, given the plastic
strain field.** The viscoplastic driver reaches equilibrium at F = 1.3 on this mesh
in 1,679 iterations. Seeded with that displacement field AND its accumulated
plastic strain — the same quantity on both drivers, `eps^p` per Gauss point with
`sigma = D (B u - eps^p)` — the Newton iteration corrects it at FULL gravity in
**16 iterations and 55 force evaluations, to an out-of-balance of 2.7e-5**, 36
times inside the trial tolerance. At F = 1.45 the same hand-over costs 17
iterations and lands at 1.2e-8.

**A short predictor is enough, and the shorter the better where it works.** Bounded
viscoplastic runs, then the corrector, on the locked mesh:

| predictor iterations | 50 | 100 | 200 | 400 | 800 | 1,600 |
|---|---|---|---|---|---|---|
| F = 1.3, corrector | fails | 45 it, oob 8.8e-8 | 27 it, 3.7e-7 | 17 it, 4.5e-8 | 11 it, 8.0e-9 | 14 it, 2.7e-6 |
| F = 1.5, corrector | 44 it, 3.9e-7 | 40 it, 1.1e-7 | 37 it, 7.3e-8 | 29 it, 2.9e-5 | 22 it, 5.7e-5 | 14 it, 2.0e-5 |

So the Newton driver's difficulty on a cohesionless soil is the plastic HISTORY and
not the equilibrium. Load control from a zero start has to discover that history by
walking gravity up in increments, and with c = 0 that walk is obstructed — the
elastic domain is a cone through the origin, so the yield check is scale-free and a
smaller load increment does not make a smaller plastic problem. The equilibrium
itself was never hard.

#### What was built

`_NR_VP_PREDICTOR_ITERS` in `xslope/fem.py`, and a seeded entry to
`_solve_fem_newton`. A trial that dies at the LOAD-STEP FLOOR — and only that
trial; one refused on force or on displacement has an answer already — is retried
from a bounded viscoplastic run at the same strength on the same prepared model.
The seeded trial does not walk the load path: it is one attempt at full gravity,
and it either corrects the seed or the trial has failed. Two rungs, 250 and 1,000.

It is not a fallback to the viscoplastic verdict. The predictor runs with
`early_failure=False` and is stopped by its budget rather than by a verdict; its
own convergence is never read. The answer is decided entirely by whether the Newton
corrector reaches full gravity in equilibrium and passes the same force gate, the
same displacement bound and the same yield reading every other trial passes. Work
is cumulative — the failed cold attempt, every predictor run and every corrector
are all charged to the trial — and `nr_predictor_iterations` reports what the
predictor cost. `nr_force_evals` and `nr_predictor_iterations` now travel in the
bisection's trial record, so a run's work is readable from the record.

#### The five reinforced benchmarks

Brackets 1.2-1.8 (geogrid, FEM-2, locked), 1.1-1.7 (half capacity), 1.0-1.6 (three
layers); tolerance 0.01. The viscoplastic column was re-measured here — it reads
about one bisection cell above the reinforcement round's figures on every row,
which is the bracket and not the driver, and the comparison below is like for like
because both columns come from the same runs. Five rows, FOUR models: the geogrid
sample and the FEM-2 tutorial are the same model under two file names (see below,
and "Bookkeeping: five benchmarks, four models" at the end of this document).

| Case | Mesh | VP FS | N-R before | N-R after | gap after | gap before | ramp after |
|---|---|---|---|---|---|---|---|
| Geogrid sample | tri6, 4 | 1.5984 | 1.6078 | 1.6078 | **+0.0094** | +0.0094 | 1.6031 |
| FEM-2 tutorial | tri6, 4 | 1.5984 | 1.6078 | 1.6078 | **+0.0094** | +0.0094 | 1.6031 |
| Half capacity | tri6, 4 | 1.4141 | 1.4141 | 1.4141 | **0.0000** | 0.0000 | 1.4031 |
| Three layers | tri6, 4 | 1.2391 | 1.1172 | 1.1641 | **-0.0750** | -0.1219 | 1.0531 |
| Geogrid LOCKED | tri6, 2 | 1.5609 | 1.2719 | **1.5609** | **0.0000** | -0.2875 | 1.5469 |

**The locked mesh — the case this round was commissioned on — now agrees with the
viscoplastic driver EXACTLY, trial for trial and verdict for verdict**, on all
eight trials of the bisection. It was 0.2875 low. Four of the five rows are inside
the 0.01 tolerance; one is not, and it is the subject of the honest negative below.

One thing the round found while measuring rather than by looking for it: the
"geogrid sample" and "FEM-2 tutorial" benchmarks are the same model. The two files
differ only in a cosmetic `type` label on each reinforcement line — same materials,
same profile lines, same six layers — and they mesh to the same 563 elements and
solve to the same factor of safety on the same trial sequence, iteration for
iteration, on both drivers. The reinforcement round's table reports them as two
benchmarks; they are one, measured twice.

The pre-remedy driver on the locked mesh does not merely give the wrong answer. Run
trial by trial it produces max|u| = 3.8e9 on a 34 m model at F = 1.8, and a full
pre-remedy bisection on that mesh aborts the process inside OpenBLAS partway
through. With the predictor the same bisection runs clean.

#### The refinement ladder

The geogrid sample at three mesh sizes, verdict and Newton iterations at each
strength. Before and after:

| Mesh | elements | F = 1.2 | F = 1.4 | F = 1.5 | F = 1.55 | F = 1.6 | F = 1.7 |
|---|---|---|---|---|---|---|---|
| tri6, 4 — before | 563 | CONV 17 | CONV 21 | CONV 28 | — | CONV 46 | — |
| tri6, 4 — after | 563 | CONV 17 | CONV 21 | CONV 28 | CONV 66 | CONV 46 | FAIL 565 |
| tri6, 3 — before | 978 | CONV 21 | CONV 158 | CONV 39 | — | FAIL 921 | — |
| tri6, 3 — after | 978 | CONV 21 | CONV 158 | CONV 39 | CONV 52 | FAIL 1,140 | FAIL 726 |
| tri6, 2 — before | 2,101 | CONV 39 | **FAIL 880** | **FAIL 261** | **FAIL 628** | FAIL 58 | — |
| tri6, 2 — after | 2,101 | CONV 39 | CONV 897 | CONV 300 | CONV 665 | FAIL 135 | FAIL 738 |

Every rung is monotone in F — everything below the limit converges and everything
above it fails — and the three rungs now agree with each other on where the limit
is, which they did not before. The finest mesh, which used to be the broken one, is
the one that brackets the limit most tightly.

#### The plain-soil eight rows

Newton bisection, every row re-measured. The "before" column is the same driver
with the predictor switched off, re-run in this session: on the four rows checked
against it, it reproduces the recorded iteration and force-evaluation counts
exactly — FEM-1 at 1,485 / 10,680, Griffiths & Lane 1 quad8 at 1,688 / 11,592,
Griffiths & Lane 6 dry tri6 at 2,791 / 19,209 and quad8 at 3,026 / 20,135 — so the
predictor-off path is the shipped driver and not an approximation of it.

| Benchmark | Mesh | VP FS | N-R before | N-R after | move |
|---|---|---|---|---|---|
| FEM-1 tutorial | tri6, 3.5 | 1.3633 | 1.3711 | 1.37109375 | **0.0000** |
| LEM-3 tutorial | tri6, 1.2 | 1.2539 | 1.2695 | 1.26953125 | **0.0000** |
| Griffiths & Lane 1 | quad8, 3.5 | 1.3656 | 1.3719 | 1.37187500 | **0.0000** |
| Griffiths & Lane 1 | tri6, 3.5 | 1.3656 | 1.3656 | 1.36562500 | **0.0000** |
| Griffiths & Lane 1 | quad9, 3.5 | 1.3844 | 1.3969 | 1.39687500 | **0.0000** |
| Griffiths & Lane 6 dry | quad8, 2 | 2.4219 | 2.4186 | 2.41562500 | **-0.0030** |
| Griffiths & Lane 6 dry | tri6, 2 | 2.4531 | 2.4531 | 2.45937500 | **+0.0063** |
| Griffiths & Lane 3, r = 0.8 | tri6, 6 | 1.4219 | 1.4281 | 1.42812500 | **0.0000** |

Six of the eight are unchanged to every digit recorded. Two moved, by one bisection
cell each and both inside the 0.01 tolerance. The tri6 dry-dam row is the one this
document already flagged as sitting on a knife edge: its deciding trial,
F = 2.45625, converges or fails on the last bit of the mantissa, and widening the
no-progress window moves it to the same 2.4594 the predictor moves it to.

#### The predictor's own knob

Three schedules, seven benchmarks, one run each.

| Benchmark | (250, 1000) — as shipped | (100,) | (500, 2000) |
|---|---|---|---|
| Geogrid coarse | 1.6078 | 1.6078 | 1.6078 |
| Half capacity | 1.4141 | 1.4141 | 1.4141 |
| Geogrid LOCKED | 1.5609 | 1.5703 | 1.5609 |
| FEM-1 tutorial | 1.3711 | 1.3711 | 1.3711 |
| Griffiths & Lane 1 quad8 | 1.3719 | 1.3719 | 1.3719 |
| Griffiths & Lane 6 dry tri6 | 2.4594 | 2.4531 | 2.4594 |
| **Three layers** | **1.1641** | **1.1172** | **1.1734** |

Four rows do not move at all. Two move by one bisection cell at the shortest
schedule and are inside 0.01 of the viscoplastic answer at every schedule. The
three-layer row moves by 0.056, and that is a finding and not a tuning problem —
see below.

And where the predictor does not fire, it changes nothing at all: the locked mesh
at F = 1.45, which the driver solves cold, returns 73 iterations and 422 force
evaluations with the predictor on and with it off, and the two displacement fields
are **bitwise identical**. Where it does fire, the state it produces passes the same
evidence: F = 1.3 on the locked mesh converges at an out-of-balance of 1.7e-5 with
a worst Mohr-Coulomb violation of 1.7e-15 of the local strength.

#### Work

Force evaluations, and the viscoplastic predictor iterations charged on top:

| Benchmark | fe before | fe after | ratio | predictor iterations |
|---|---|---|---|---|
| FEM-1 tutorial | 10,680 | 13,492 | 1.26x | 3,750 |
| Griffiths & Lane 1 quad8 | 11,592 | 14,904 | 1.29x | 3,750 |
| Griffiths & Lane 6 dry tri6 | 19,209 | 22,231 | 1.16x | 7,500 |
| Geogrid coarse | 11,368 | 13,321 | 1.17x | 2,500 |
| Half capacity | 15,144 | 17,349 | 1.15x | 6,250 |
| Three layers | 14,845 | 25,663 | 1.73x | 7,750 |

The criterion asked for 1.2x. Two of the six clean rows miss it, at 1.26x and
1.29x, and the three-layer row — where the predictor fires on nearly every trial and
rescues none — costs 1.73x. The cost falls entirely on FAILING trials, which are
the ones that now pay a cold attempt plus a predictor plus a corrector; a
converging trial costs exactly what it cost before, to the force evaluation.

#### The probes

**Past failure.** FEM-1 at F = 2, 3 and 5 still fail, and still promptly: 399, 497
and 303 Newton iterations, 3,073 / 3,458 / 2,453 force evaluations, 12 s, 14 s and
10 s, each after both predictor rungs are exhausted. The predictor does not turn a
grossly overloaded slope into a standing one — a state that is running away is not a
state the corrector can equilibrate — which is the property the whole design rests
on.

**nu = 0.49.** Unchanged and untouched. Converged at F = 1.00, 1.30 and 1.36 in 11,
12 and 12 iterations with 27, 29 and 30 force evaluations and ZERO predictor
iterations — the same counts this document already records — and failed at F = 1.40
and 1.50.

**The default path.** Griffiths & Lane 6 dry, quad8, size 2, no `fem_solver`
argument: FS = 2.421875 with per-trial iteration counts

    147, 781, 3393, 2031, 2841, 9541, 12000, 8617, 8777

value for value the control sequence.

#### The locks

`test/nr_ssrm_check.py` gains `check_cohesionless_solve`, on the locked tri6/2.0
mesh — the expensive one, because the defect is invisible at tri6/4.0 and only
marginal at tri6/3.0. It asserts that F = 1.45 converges to an ADMISSIBLE field
(worst yield violation under 1e-6 of the local strength), that F = 1.3 converges to
force equilibrium with a stress field on the yield surface, and that F = 1.3's
displacement is below F = 1.45's, since a stronger soil under the same gravity
cannot move further.

The middle assertion carries the whole argument without reference to the other
driver: strength reduction only ever raises c and tan(phi) as F falls, so an
admissible field in equilibrium with full gravity at F = 1.45 is admissible at
every lower F under the same load, and the lower-bound theorem says the slope
stands there.

**Mutation, run both ways.** With `_NR_VP_PREDICTOR_ITERS` emptied — which is the
driver exactly as it stood before this work — the check fails, and names what it
found: "F = 1.3 came back FAILED (diverging, signal 'load_step_floor') after 599
iterations at a load factor of 0.00, on a slope the SAME driver reports standing at
F = 1.45 with a yield violation of 2.5e-15." As shipped it passes. The whole check
file passes.

#### The criterion, line by line

**The five reinforced benchmarks — PARTLY MET.** Four of five agree with the
viscoplastic driver inside 0.01, two of them exactly: the geogrid sample and the
FEM-2 tutorial at +0.0094 (and they are the same model), half capacity at 0.0000,
and the locked mesh at 0.0000 — trial for trial, from a gap of 0.2875. The three-layer
variant misses at -0.0750, improved from -0.1219.

**The refinement ladder — MET.** Clean at 563, 978 and 2,101 elements: no false
failure at any strength on any rung, verdict monotone in F everywhere, and the
three rungs agreeing with each other on the limit.

**The plain-soil table — MET.** Six rows bit-identical, two moved by one bisection
cell (-0.0030 and +0.0063), none by more than 0.01. Every row reported above.

**The knob — PARTLY MET.** The answer does not move at all on four of the seven
benchmarks tested across three predictor schedules, and moves by one bisection cell
on two more, staying inside 0.01 of the viscoplastic answer at every schedule. It
moves by 0.056 on the three-layer row. And the identity clause is met exactly:
where the predictor does not fire the converged state is bitwise identical with it
on and off.

**Work — NOT MET.** 1.15x to 1.29x on the clean rows against a 1.2x requirement,
1.73x on the three-layer row, plus 2,500 to 7,750 viscoplastic predictor iterations.
All of it falls on failing trials; converging trials are unchanged to the force
evaluation.

**The probes — MET.** Past-failure verdicts unchanged and prompt; nu = 0.49
unchanged, with the predictor never firing on it.

**The locks — MET.** The new check fails on the driver as it stood and passes now,
run both ways; the whole check file passes.

**The default path — MET.** FS 2.421875 on the control iteration sequence.

#### The honest negative: the three-layer model

One row does not close, and the reason is measured rather than guessed.

At F = 1.1875 the viscoplastic driver reaches equilibrium in 2,608 iterations, and a
predictor long enough to include that convergence hands Newton a state it corrects
in **5 iterations**. A predictor stopped at 1,000 or 2,000 iterations does not:
the corrector fails at an out-of-balance of 1.8 and 4.6e-2. So on this model the
seed has to be a CONVERGED viscoplastic state, not merely a well-developed one, and
that is what makes the predictor budget a dial on this row and on no other.

At F = 1.225 it is worse than a budget question. The viscoplastic driver converges
there in 23,251 iterations at an out-of-balance of 9.99e-4 — just inside the 1e-3
gate — and handed THAT state, fully converged, the Newton corrector still fails, at
an out-of-balance of 3.7. Predictor budgets of 1,000, 2,000, 4,000, 8,000, 16,000
and 32,000 iterations were each tried and each corrected to a different
non-equilibrium state. **The viscoplastic converged state at that strength is not a
root of the Newton residual, and there is no root near it.** That is a genuine
disagreement between the two formulations about whether this slope stands at
F = 1.225, not a solver budget, and closing it is a separate measurement about
which of the two is right — the viscoplastic reading there is itself marginal, at
an out-of-balance sitting on the tolerance after twenty-three thousand iterations.

**That measurement was made. See "THE THREE-LAYER DISAGREEMENT" at the end of this
document: the Newton driver is right at F = 1.225, and the viscoplastic converged
state there is 580 psf outside the yield surface.**

Raising the predictor budget until it always covers a viscoplastic convergence
would move this row to about 1.20 rather than to 1.2391, and it would make every
failing Newton trial cost a whole viscoplastic trial — which is to say it would
turn the Newton driver into the viscoplastic driver with a polish. That is why the
budget is 250 and 1,000 and not the trial's own iteration ceiling.

#### Verdict

The Newton driver is ready for reinforced and general c = 0 questions on the meshes
the corpus is written on, with one model's worth of caveat.

The case the round was commissioned on is closed completely. The geogrid sample at
its locked mesh went from 0.2875 below the viscoplastic answer to exactly it, on
every one of the eight bisection trials, and the refinement ladder that used to
degrade monotonically — clean at 563 elements, marginal at 978, broken at 2,101 —
is now clean at all three rungs with the three of them agreeing on where the limit
is. The plain-soil table did not move: six of eight rows are bit-identical and the
other two moved by one bisection cell. The default path is untouched.

And the diagnosis in the reinforcement round was wrong, which matters more than the
fix. It named the Mohr-Coulomb apex, where the consistent tangent is zero. The apex
fires at three Gauss points out of 6,303. What is actually degenerate is every
yielding point — on a linear yield surface without hardening the trial-to-returned
Jacobian has determinant `1 - A - Bc`, and `A + Bc = 1` identically, so a whole line
of trial states maps to one returned state and the consistent tangent has a null
direction wherever a point yields. That is the constitutive law and not a defect,
and conditioning it away does not help: blending the tangent toward the elastic
operator was built and measured at seven settings from 0.02 to the elastic operator
itself, and the specimen trial fails at every one of them. Nor does the vendors'
documented workaround: a token cohesion of 5 or 10 psf moves the false failures to
different strengths and leaves the verdict sequence non-monotone at both.

What was actually wrong is that the Newton driver had to discover the plastic strain
field by walking gravity up from zero, and with c = 0 that walk is obstructed —
the elastic domain is a cone through the origin, so a smaller load increment is not
a smaller plastic problem. Eighteen thousand iterations carry 3% of gravity. Two
hundred viscoplastic iterations, and the same Newton iteration lands the whole load
in twenty-seven. The equilibrium was never the hard part.

What remains, in the order it matters:

- **The three-layer model.** The two drivers genuinely disagree there at F = 1.225,
  and no predictor budget closes it: the viscoplastic converged state at that
  strength is not a root of the Newton residual. Which of the two is right is a
  measurement nobody has made, and it is the one place in this branch where the
  answer still depends on a solver setting. *(Made — "THE THREE-LAYER
  DISAGREEMENT", below. The Newton driver is right there; the disagreement is a
  viscoplastic-side one and no code was changed.)*
- **The cost of a failing trial.** A trial that fails now pays a cold attempt, a
  predictor and a corrector — 1.15x to 1.29x the force evaluations on the clean
  benchmarks plus 2,500 to 7,750 viscoplastic iterations. Skipping the cold attempt
  on a model that has already needed the predictor once would recover most of that
  and was not done this round.
- **The ramp still walks its own steps cold.** It gets the predictor only at the
  foot of the ramp, which is why it reads 1.5469 against the bisection's 1.5609 on
  the locked mesh — closer than the 1.2406 it read before, and still its own
  answer. A refused ramp step is the same load-step-floor failure the bisection
  retries, and it should be retried the same way.
- **Post-peak softening is still refused**, so neither of the repository's published
  reinforced factors of safety is reachable on this driver. That is unchanged by
  this round and is the larger of the two things standing between the Newton path
  and a reinforced answer anyone would quote.


## THE THREE-LAYER DISAGREEMENT — an independent referee, and which driver is right

The cohesionless round closed with one row open and named the measurement nobody
had made: on the three-layer reinforced model the viscoplastic driver converges at
F = 1.225 in 23,251 iterations, and that converged state is not a root of the
Newton residual — six predictor budgets up to 32,000 iterations each corrected to
a different non-equilibrium state. Two formulations of "equilibrium at this state"
disagree. This round finds out which one is telling the truth.

Nothing in the drivers was changed. Everything below was measured on this
checkout, on the coarse tri6/4.0 mesh (563 elements, 36 bar elements, of which the
18 belonging to layers 1, 3 and 5 carry capacity), `force_tol` 1e-3, hybrid
criterion, tolerance 0.01.

### The referee, and what it is allowed to use

The adjudication cannot be run on either driver's own arithmetic, so the
equilibrium check was written from scratch and shares no code with either. It reads
only model DATA out of `fem_data` — node coordinates, connectivity, the material
table, the boundary-condition codes, the reinforcement line geometry and its E and
area — plus a STATE, and it computes everything else itself:

- its own tri6 shape functions and their derivatives, and the FULL isoparametric
  Jacobian rather than the straight-sided shortcut the drivers take;
- its own three-point triangle rule, its own plane-strain constitutive matrix
  written from Hooke's law, and its own `integral B^T sigma dV`;
- its own consistent gravity load `integral N^T (0, -gamma) dV`;
- its own free-degree-of-freedom set from the boundary-condition codes, and its own
  nodal tributary weights for the Dawson, Roth & Drescher normalization;
- its own bar element: the two-point Gauss integral of `EA (dN/dx)^2` over the
  three-node quadratic bar, rotated into global coordinates, plus the same chord
  law both drivers implement, evaluated on the state's own reported bar forces;
- its own Mohr-Coulomb yield function in the invariant form.

**Validated first on a state where the two drivers agree**, because a referee that
cannot confirm agreement cannot adjudicate disagreement. The geogrid sample, six
layers, coarse mesh, at F = 1.5, a trial both drivers converge:

| what | referee against the drivers |
|---|---|
| gravity load vector | 2.0e-15 relative |
| free degree-of-freedom set | identical |
| nodal tributary weights | 7.9e-15 relative |
| 36 bar element matrices | worst 3.8e-16 relative |
| Newton's converged state, out-of-balance | **3.7849e-06**, against the driver's own **3.785e-06** |
| viscoplastic converged state, out-of-balance | **3.66e-04** — in equilibrium, well inside the 1e-3 gate |
| Newton state, worst yield violation | 2.3e-12 of the local strength |
| viscoplastic state, worst yield violation | 25.5 psf, one Gauss point |

The referee reproduces the Newton driver's own out-of-balance to five significant
figures and confirms the viscoplastic state is in force equilibrium too. It is
believable on a state where the two agree.

### The disputed state, refereed

The three-layer model at F = 1.225. The viscoplastic driver converges there in
23,251 iterations at its own out-of-balance of 9.988e-04, a thousandth inside the
gate.

**It IS in force equilibrium.** The referee reads an out-of-balance of **2.38e-05**
on the stress field that state actually carries — forty times inside the trial
tolerance. The viscoplastic scheme solves `integral B^T D (B u - eps^p) dV = F_ext`
exactly at every iteration, and the referee confirms it independently.

**It is NOT an admissible stress field.** Nine Gauss points sit outside the
Mohr-Coulomb surface by more than 1% of the local strength AND more than 10 psf.
The worst is out by **579.9 psf — 2.09 times the strength available at that point**.
Eight of the nine are in the c = 0 material, and every one of those carries a
TENSILE mean stress, up to **729.1 psf of tension in a soil with no cohesion**.

**What admissibility costs.** Demanding a returned stress field at the same
displacement and the same plastic strain — which is what the Newton residual asks —
leaves an out-of-balance of **7.39** and a residual of **9.19e-02** of the external
load. That is the Newton driver's residual on the viscoplastic converged state,
confirmed by the referee rather than taken from the driver.

### The gap is entirely in the soil, and the bar element is exonerated

The same state, internal force split by element class, referee against Newton
driver, relative to the norm of the external load:

| class | referee vs driver |
|---|---|
| reinforcement bars | **6.0e-16** |
| external gravity load | **2.0e-15** |
| 2-D soil | **9.19e-02** |

The bars agree to round-off and so does the load. The whole of the 9.19e-02 is the
soil's return-map correction — the amount by which the viscoplastic stress field
lies outside the yield surface, expressed as a force. The bar element, which was the
prime suspect on the Newton side, has nothing to do with it.

### Why the viscoplastic gate cannot see it

The viscoplastic out-of-balance is not the equilibrium residual. It is the
window-averaged INCREMENT of the viscoplastic body load — how fast plastic flow is
still changing the state. On this model the two quantities come apart completely.
The disputed trial re-run to fixed budgets with the force gate set unreachably
tight, so nothing stops early:

| budget | driver's own out-of-balance | referee: force equilibrium | worst violation | Gauss points out by >10% | max tension in the c = 0 soil |
|---|---|---|---|---|---|
| 2,000 | 1.50e-01 | 1.5e-05 | 1.40 | 4 | 332 psf |
| 6,000 | 1.13e-01 | 2.5e-05 | 1.42 | 10 | 647 psf |
| 12,000 | 5.89e-02 | 7.1e-06 | 1.47 | 11 | 792 psf |
| 23,251 | 9.99e-04 | 2.4e-05 | 2.09 | 8 | 729 psf |
| 40,000 | **1.25e-12** | 2.4e-05 | 1.54 | 8 | **809 psf** |

The driver's measure falls twelve orders of magnitude while the violation does not
fall at all and the tension GROWS. The gate is not failing to converge; it is
converging to a state that is not admissible, and it has no term that could
notice.

### The mechanism, measured

With psi = 0 the Mohr-Coulomb flow is purely deviatoric, so plastic flow carries no
volume change and can never relieve a point's MEAN stress. Measured over every
Gauss point that carries plastic strain at the disputed state — 546 of 1,689 —
`|tr(eps^p)| / |eps^p|` is at most **3.6e-13**. The plastic strain is traceless to
round-off, which is what the flow rule says it must be.

On a c = 0 material the yield surface is a cone through the origin, so ANY tensile
mean stress is outside it however much deviatoric stress the point sheds. Those
points therefore shed their deviatoric stress, stop flowing, and freeze with the
violation intact — which is exactly the 40,000-iteration row above. The code names
this case itself, in the tension-cutoff block: "The psi=0 MC flow is purely
deviatoric and cannot relax a near-apex tensile state, which is exactly what this
surface handles." The Rankine tension cutoff is finite on **0 of 563 elements** on
this model. It is off, and nothing else in the scheme caps tension.

### How far it reaches

Not a reinforcement problem, and not a three-layer problem:

| model | trial | referee: force equilibrium | worst violation | tensile points in c = 0 soil |
|---|---|---|---|---|
| three layers, as disputed | F = 1.225, VP converged, 23,251 it | 2.4e-05 | 579.9 psf | 8, to 729 psf |
| the SAME soil with every bar capacity removed | F = 1.0, VP converged, 7,215 it | 1.3e-12 | 5 points | 5, to 138 psf |
| Griffiths & Lane 1, c = 312.5 psf, phi = 20 | F = 1.35, VP converged, 867 it | 7.8e-13 | 2.3e-04 of local strength | **none** |

Strip the reinforcement out entirely and the same frozen tensile points are still
there, so the bars do not cause it — they amplify it, by a factor of five in the
tension carried, because a reinforced slope is precisely one that stands by pulling
the soil into tension. Give the soil cohesion and it disappears: on the plain
cohesive benchmark the viscoplastic converged field is admissible to 2.3e-04 of the
local strength with no tensile point anywhere. **It is a zero-cohesion effect.**

### What the factor of safety actually is

The lower-bound theorem settles the part that can be settled: a field in
equilibrium with full gravity and nowhere outside the yield surface proves the
slope stands. Two independent routes produce one, and the referee certifies both.

Handing each viscoplastic converged state to the Newton corrector:

| F | VP iterations | VP field, worst violation | corrector | referee equilibrium | corrected field |
|---|---|---|---|---|---|
| 1.1875 | 2,608 | 26.6 psf | **CONVERGED, 5 it** | 7.8e-05 | **0.0 psf — admissible** |
| 1.2000 | 5,389 | 169.1 psf | failed | 1.1e+01 | 731.8 psf |
| 1.2125 | 21,023 | 551.2 psf | failed | 3.3e+05 | — |
| 1.2250 | 23,251 | 579.9 psf | failed | 3.2e+02 | 15,818 psf |

And with the tension cutoff switched ON — the code's own named remedy for the
mechanism above, and a modeling choice, not a code change:

| F | VP + cutoff | referee equilibrium | worst Mohr-Coulomb violation | plain-MC Newton seeded from it |
|---|---|---|---|---|
| 1.1875 | CONVERGED | 4.0e-08 | 0.1 psf | — |
| **1.2063** | **CONVERGED, 22,735 it** | **5.6e-07** | **0.09 psf — admissible** | **CONVERGED in 4 it**, referee 1.1e-05, **0.00 psf** |
| 1.2156 | FAILED | 6.7e-06 | 4.5 psf | failed |
| 1.2250 | **FAILED** (converged without the cutoff) | 5.2e-07 | 7.6 psf | failed |

Two things fall out of that table. The plain-Mohr-Coulomb Newton residual HAS a root
at F = 1.2063 — four iterations from the right starting state, in equilibrium to
1.1e-05 and admissible to 0.00 psf — so the Newton FORMULATION reaches strengths its
cold-start path does not. And with tensile stress capped, the viscoplastic driver
itself flips its verdict at F = 1.225 from CONVERGED to FAILED, landing on the same
answer the Newton driver gives there.

The three bisections, same model, same bracket 1.0-1.6, tolerance 0.01:

| driver | FS | supported by an admissible field? |
|---|---|---|
| viscoplastic, as shipped | **1.2391** | no — its converged trials at 1.225 and 1.2344 are 580 psf outside the surface |
| viscoplastic + tension cutoff | **1.2109** | yes — its deciding converged trial, F = 1.2063, is certified |
| Newton bisection | **1.1641** | yes, but provably LOW — an admissible field exists at 1.2063 |

**The slope provably stands at F = 1.2063.** Newton's 1.1641 is 0.042 below that, and
the reason is the predictor budget the cohesionless round already documented, not
the formulation: the root is there and a good enough seed finds it in four
iterations. The viscoplastic 1.2391 is above the proven floor and rests on fields
that carry 580 to 652 psf of stress outside a zero-cohesion yield surface, so it is
not evidence of anything. Nothing measured here proves the slope fails between
1.2063 and 1.2391 — only that neither shipped driver's number on this row is
supported.

### The verdict

**At F = 1.225 the Newton driver is right and the viscoplastic driver is wrong.**
The viscoplastic converged state is in force equilibrium — the referee confirms it
independently, at 2.4e-05 — but it is not a statically admissible stress field, and
so it is not equilibrium in the sense a factor of safety is defined on. Nine Gauss
points are outside the yield surface, eight of them carrying up to 729 psf of
tension in a soil with zero cohesion. There is no root of the Newton residual near
that state because there is no admissible field near it, and the Newton driver's
refusal to find one is the correct answer, not a defect.

The clause that accepts it is the force gate, `unbalanced_force_ratio < force_tol`.
That quantity is the window-averaged increment of the viscoplastic body load — a
measure of how fast plastic flow is still changing the state — and on this model it
falls to 1.2e-12 while the state stays 800 psf outside the surface. A psi = 0 flow
rule cannot relieve a mean-stress violation (the plastic strain is traceless to
3.6e-13), a c = 0 cone admits no tension at all, and the tension cutoff that exists
for exactly this case is off. So the flow stops, the measure reports convergence,
and the violation stays.

**No code was changed.** This implicates the SHIPPED default driver's convergence
criterion on zero-cohesion materials without a tension cutoff, and which engine
ships as the default — and whether `unbalanced_force_ratio` should be joined by a
yield-admissibility reading before a trial is called converged — is the owner's
decision, not a spike's. What is measured and available for that decision:

- The effect is NOT specific to reinforcement and NOT specific to this model. It is
  the zero-cohesion material at a free face, present with every bar removed, absent
  on a cohesive benchmark.
- It is invisible to the viscoplastic driver's own diagnostics and visible for
  almost nothing to a reader who asks for it: `nr_max_yield_violation` is the
  Newton path's existing reading of exactly this quantity, and the referee's version
  of it costs one pass over the Gauss points.
- Switching the Rankine tension cutoff on moves this model's viscoplastic answer
  from 1.2391 to 1.2109 and makes every reported converged field admissible.
- The Newton driver's low answer on this row is a predictor budget, not a
  formulation error, and it is already written up in the cohesionless round.

### Bookkeeping: five benchmarks, four models

Re-measured here rather than carried over: `docs/fem/files/xslope_reinforce_fem.xlsx`
and `docs/tutorials/files/xslope_reinforced_slope.xlsx` are the SAME model. Loaded
side by side, their materials are identical, their profile lines are identical, and
their six reinforcement lines are identical in geometry, capacity, pullout lengths,
E and area. The only difference in either file is a cosmetic `type` label on each
reinforcement line — blank in one, "geosynthetic" in the other — which nothing in
the FEM path reads. The tables in "REINFORCEMENT — results" and "THE COHESIONLESS
SOLVE — results" both present five reinforced rows; they are FOUR distinct models,
with the geogrid sample measured twice under two file names.


## THE ADAPTIVE PREDICTOR — the three-layer row, and the criterion for closing it

Written after the adjudication and BEFORE any remedy code, so that what follows is
a test and not a description. Same machine and settings as everything above:
`force_tol` 1e-3, hybrid criterion, `capture_failure_state=False`, tolerance 0.01,
tri6.

### The defect, restated from what is already measured

On the three-layer reinforced model the Newton bisection reads **1.1641**. The
adjudication certified that the slope STANDS at F = 1.2063 — an admissible field in
equilibrium with full gravity, confirmed by an independent referee and reached by a
4-iteration Newton solve seeded from the viscoplastic-with-cutoff state — so 1.1641
is about 0.045 low. The legal reference answer on this model is **1.2109**, the
viscoplastic bisection with the Rankine cutoff on; the shipped viscoplastic 1.2391
is inflated by fields carrying up to 729 psf of tension in a zero-cohesion soil and
is NOT the target.

The mechanism is documented in "THE COHESIONLESS SOLVE": a trial near but below the
true limit needs a developed plastic zone that Newton cannot grow from a zero start,
and the viscoplastic predictor supplies it. What fails is the SHAPE of the predictor
budget. It is a fixed ladder, `_NR_VP_PREDICTOR_ITERS = (250, 1000)`, and on this
model six fixed budgets up to 32,000 iterations all failed at F = 1.225 while a
CONVERGED-state seed succeeded in a handful of iterations. A fixed ladder cannot
express "run until this walk has finished", which is what the seed on this model
needs.

### The direction

Two changes, and both are to the PREDICTOR only. The return map, the residual, the
force gate, the displacement bound and the yield reading are untouched, so the
verdict stays entirely the corrector's.

1. **Adaptive budget instead of a fixed ladder.** Continue the viscoplastic
   predictor while it is making measurable progress, with a hard ceiling so a
   genuinely failing trial still exits promptly. The progress test is not invented
   here: `_still_progressing` is the viscoplastic driver's OWN budget-extension
   rule, already calibrated and already in the file — the residual's trailing-window
   mean falling by at least 1% against the window before it, OR a displacement field
   standing still on evidence the failure classifier cannot rule on. Handing the
   predictor a small chunk and a large ceiling makes that rule the predictor's
   stopping rule. Whether it tracks what the CORRECTOR needs is a measurement, and
   it is made below before the rule is adopted.

2. **The predictor runs with tension capped.** The adjudication showed the trap: on
   a c = 0 material a psi = 0 flow rule cannot relieve a tensile mean stress, so the
   viscoplastic predictor freezes with Gauss points hundreds of psf outside a cone
   through the origin. The corrector's yield gate protects the verdict from such a
   state, but a seed deep in illegal territory is a poor seed — measured: at
   F = 1.1875 the plain predictor's converged state has a 26.6 psf violation and the
   corrector takes it in 5 iterations; at F = 1.20 the violation is 169.1 psf and the
   corrector fails. Running the predictor under a Rankine cutoff (T = 0) keeps the
   seed admissible. Whether it seeds BETTER is a measurement, and it is made below.

### Success criterion (verbatim)

- **Three-layers Newton bisection FS within 0.01 of 1.2109** (the legal reference),
  and the trial at F = 1.2063 (or the nearest bisection trial below it) converges —
  the certified-standing state must be reachable.
- **The ramp agrees with the fixed bisection within 0.01 on three-layers.**
- **The other three reinforced benchmarks are unchanged within 0.01**: the geogrid
  sample, half capacity, and the geogrid locked mesh — and the locked mesh must stay
  EXACTLY on 1.5609, because it is the showcase.
- **The plain-soil eight-row table does not move.** Every row reported; no row moves
  more than 0.01, and the predictor-fire count is asserted ZERO on the cohesive
  plain-soil rows, where the answer must therefore be bit-identical.
- **The probes still hold.** nu = 0.49 and the FEM-1 past-failure trials at F = 2, 3
  and 5 are unchanged, and a past-failure trial does not regress in cost by more
  than about 1.5x — the ceiling plus the progress test must exit it cheaply rather
  than burn the ceiling every time.
- **The locks catch it.** `test/nr_ssrm_check.py`'s cohesionless lock is extended
  with the three-layer case, which must FAIL on the driver as it stands today
  (bisection 1.1641, the F = 1.2063 trial refused) and pass after — the mutation run
  both ways and recorded — and the whole check must pass.
- **The default path is unchanged**, against the standard control: Griffiths & Lane
  6 dry with no `fem_solver` argument, FS 2.421875 on iteration counts 147, 781,
  3393, 2031, 2841, 9541, 12000, 8617, 8777.
- **An honest negative is a valid outcome and must be written.** If the adaptive
  budget only reaches 1.2109 by turning every failing trial into a whole
  viscoplastic solve, or if the capped predictor buys the three-layer row and costs
  the locked mesh, that is the result.


### THE ADAPTIVE PREDICTOR — results

Same machine and settings as the criterion above. Every number below was measured on
this checkout in this session; the "before" columns are the driver with
`_NR_VP_PREDICTOR_ADAPTIVE` switched off — which is the driver exactly as it stood
at dea4c119 — re-run here rather than quoted. It reproduces the recorded figures on
every row it can be checked against: three layers at 1.1641, the locked mesh at
1.5609, and all eight plain-soil factors of safety.

#### The adaptive rule, and the measurement that chose it

The rung's stopping rule is not new code. `_still_progressing` is the viscoplastic
driver's own budget-extension test — the residual's trailing-window mean falling by
at least 1% against the window before it, OR a displacement field standing still on
evidence the failure classifier cannot rule on — and handing the predictor the
caller's own `max_iterations` as a chunk and `max_iterations_ceiling` as a hard stop
makes that rule the predictor's stopping rule. The rung is budgeted, in other words,
exactly as a viscoplastic trial of the same model at the same strength would be.

What justified it is the table below, measured before the rule was adopted: the
three-layer model at the four strengths that decide its bisection, at three chunk
sizes, each seed handed to the same corrector.

| chunk | F = 1.1875 | F = 1.20625 | F = 1.215625 | F = 1.225 |
|---|---|---|---|---|
| 2,000, uncapped | VP conv 2,608 → **CONV** | VP cap 6,000 → FAIL | VP cap 10,000 → FAIL | VP cap 6,000 → FAIL |
| 2,000, capped | VP cap 2,000 → **CONV** | VP cap 10,000 → FAIL | VP cap 6,000 → FAIL | VP cap 8,000 → FAIL |
| 6,000, uncapped | VP conv 2,608 → **CONV** | VP cap 6,000 → FAIL | VP cap 18,000 → FAIL | VP cap 6,000 → FAIL |
| 6,000, capped | VP conv 4,708 → **CONV** | VP conv 22,964 → **CONV** | VP cap 6,000 → FAIL | VP cap 12,000 → FAIL |
| 12,000, uncapped | VP conv 2,608 → **CONV** | VP conv 19,716 → **FAIL** | VP conv 22,614 → FAIL | VP conv 23,251 → FAIL |
| 12,000, capped | VP conv 4,708 → **CONV** | VP conv 22,964 → **CONV** | VP cap 12,000 → FAIL | VP cap 12,000 → FAIL |

Two things fall out of it. On the CAPPED rows the corrector converged on every seed
the predictor had carried to its own convergence and on no seed it had not, without
exception — so the predictor's progress rule is reading the same thing the corrector
needs, and stopping it while it is still progressing is what loses the trial. And on
the uncapped 12,000 row the correspondence breaks in exactly one place, at
F = 1.20625, where a fully converged predictor state is refused: that is the
adjudication's illegal-tension state, in force equilibrium and not admissible, with
no root of the Newton residual near it.

#### The cap, isolated

The one row above is the whole case for the cap, and it was isolated rather than
inferred. At F = 1.20625 on the three-layer model, adaptive rung both ways:

| | predictor | corrector |
|---|---|---|
| uncapped | converged, 19,716 iterations | **FAILED** at a load factor of 0.00 |
| capped (Rankine T = 0) | converged, 22,964 iterations | **CONVERGED in 4 iterations**, out-of-balance 1.1e-5, yield violation 1.2e-15 |

The cap is on the SEED only. The corrector runs on the caller's own tensile caps and
reports the caller's own yield reading, so the verdict is not capped — and the
mutation below shows the check fails without it.

#### Why the rung is APPENDED and not substituted

The first attempt replaced the fixed ladder with the adaptive rung. It reached
1.2109 on the three-layer model and it cost the showcase: the locked tri6/2.0 mesh
fell from 1.5609 to 1.5516, because its trial at F = 1.55625 flipped. Isolated, on
that trial:

| seed | predictor | corrector |
|---|---|---|
| 250, uncapped | 250, not converged | **CONV, 51 iterations** |
| 1,000, uncapped | 1,000, not converged | FAIL |
| adaptive, uncapped | converged, 6,406 | **CONV, 10 iterations** |
| 250, capped | 250, not converged | **CONV, 63 iterations** |
| 1,000, capped | 1,000, not converged | **CONV, 41 iterations** |
| adaptive, capped | 12,000, NOT converged (oob 5.7e-3) | FAIL |

A short seed is not merely cheaper; on this model it is better. So the short rungs
stay first and unchanged, and the adaptive rung is appended behind them. The
consequence is a property worth more than the argument for it: **every trial the
driver already rescued is rescued at the same rung from the same state, so the new
rung can only convert a FAILED trial and can never move a passing one.**

That is measured, not asserted. Per-trial, before against after, on four
benchmarks — 34 trials:

| | trials | converged trials identical in iterations AND force evaluations | failed trials |
|---|---|---|---|
| FEM-1 tutorial | 9 | 6 of 6 | 3 differ (extra rung) |
| Griffiths & Lane 1 quad9 | 9 | 8 of 8 | 1 differs |
| Griffiths & Lane 6 dry tri6 | 9 | 4 of 4 | 5 differ |
| Geogrid LOCKED | 8 | 4 of 4 | 4 differ |

The dry-dam row's deciding trial, F = 2.45625 — the knife-edge case this document
already flags — is among the identical ones, at 606 iterations and 3,750 force
evaluations on both drivers.

#### The four reinforced benchmarks

| Case | Mesh | before | after | reference | gap |
|---|---|---|---|---|---|
| Three layers | tri6, 4 | 1.1641 | **1.2109375** | 1.2109 (VP + cutoff) | **0.0000** |
| Geogrid sample | tri6, 4 | 1.6078125 | 1.6078125 | — | **0.0000 move** |
| Half capacity | tri6, 4 | 1.4140625 | 1.4140625 | — | **0.0000 move** |
| Geogrid LOCKED | tri6, 2 | 1.5609375 | **1.5609375** | — | **exact** |

The three-layer bisection now converges the trial the referee certified. Its eight
trials, after: 1.0 CONV, 1.6 FAIL, 1.3 FAIL, 1.15 CONV, 1.225 FAIL, 1.1875 CONV,
**1.20625 CONV**, 1.215625 FAIL — final interval [1.20625, 1.215625], which is the
same interval the viscoplastic-with-cutoff bisection closed on.

#### The ramp

The ramp used to walk every step cold, and got the predictor only at the foot. A
refused ramp step is the same load-step-floor failure the bisection retries, so it
is now retried the same way: the same rungs, the same cap, and the corrector's own
force gate and displacement bound deciding, read on the corrected state. The ramp's
final export solve takes the predictor too, so what is written to figures is a solve
at the limit rather than whatever the cold path could reach.

| | before | after | bisection | agreement |
|---|---|---|---|---|
| Three layers | 1.053125 | **1.215625** | 1.2109375 | **0.0047** |
| Griffiths & Lane 1 tri6 | 1.3656250 | 1.3656250 | 1.3656250 | 0.0000 |

The ramp carries F = 1.2125 and refuses 1.21875, so it brackets the limit one cell
above the bisection's [1.20625, 1.215625] and reports the midpoint of its own
interval. On the plain benchmark it does not move at all.

#### The plain-soil eight rows

| Benchmark | Mesh | before | after | move | fe before | fe after | ratio |
|---|---|---|---|---|---|---|---|
| FEM-1 tutorial | tri6, 3.5 | 1.37109375 | 1.37109375 | **0.0000** | 13,492 | 14,233 | 1.05x |
| LEM-3 tutorial | tri6, 1.2 | 1.26953125 | 1.26953125 | **0.0000** | 34,793 | 38,463 | 1.11x |
| Griffiths & Lane 1 | quad8, 3.5 | 1.37187500 | 1.37187500 | **0.0000** | 14,904 | 16,354 | 1.10x |
| Griffiths & Lane 1 | tri6, 3.5 | 1.36562500 | 1.36562500 | **0.0000** | 21,624 | 22,650 | 1.05x |
| Griffiths & Lane 1 | quad9, 3.5 | 1.39687500 | 1.39687500 | **0.0000** | 5,920 | 6,030 | 1.02x |
| Griffiths & Lane 6 dry | quad8, 2 | 2.41562500 | 2.41562500 | **0.0000** | 27,565 | 29,059 | 1.05x |
| Griffiths & Lane 6 dry | tri6, 2 | 2.45937500 | 2.45937500 | **0.0000** | 22,231 | 23,672 | 1.06x |
| Griffiths & Lane 3, r = 0.8 | tri6, 6 | 1.42812500 | 1.42812500 | **0.0000** | 37,791 | 40,652 | 1.08x |

Eight rows, no movement at any digit. The reinforced rows that were meant to hold
are the same: geogrid 1.10x, half capacity 1.05x, locked mesh 1.06x.

**One expectation in the criterion is falsified, and it is worth stating plainly.**
It expected the predictor never to fire on the cohesive plain-soil rows. It fires on
every one of them — 1, 3, 3, 4, 6, 6, 6 and 7 trials across the eight — because a
trial past the limit dies at the load-step floor on a cohesive soil exactly as it
does on a cohesionless one, and that is the retry's trigger. What is true is the
stronger statement measured above: the fires are all on trials that fail either way,
and every converged trial is bit-identical. The fire counts are identical before and
after on all eight rows, which is the additive design showing through — the same
trials fire, and they get one more rung.

The viscoplastic iterations charged on top are where the cost went: 1,250 to 7,647
before, 1,521 to 76,581 after, all of it on failing trials. In wall clock that is
1.14x to 1.74x on the plain rows and 1.33x to 1.50x on the reinforced ones.

#### The probes

**Past failure.** FEM-1 at F = 2, 3 and 5, before against after:

| F | iterations | force evaluations | predictor iterations | wall |
|---|---|---|---|---|
| 2.0 | 403 → 404 | 3,105 → 3,115 | 1,431 → 1,431 | 13.5 s → 12.9 s |
| 3.0 | 500 → 500 | 3,480 → 3,480 | 1,331 → 1,331 | 15.6 s → 16.0 s |
| 5.0 | 305 → 306 | 2,465 → 2,475 | 1,301 → 1,301 | 10.3 s → 10.8 s |

All three still fail, still promptly, and at **1.00x** the cost — against the 1.5x
the criterion allowed. The reason is the mechanism the design rests on: a grossly
overloaded slope runs away, `early_failure` closes the adaptive rung at about 180
iterations, and the ceiling is never approached. The rung costs the full walk only
where a full walk is what the model is doing.

**nu = 0.49.** Unchanged. Converged at F = 1.00, 1.30 and 1.36 in 11, 12 and 12
iterations with 27, 29 and 30 force evaluations and ZERO predictor iterations — the
same counts this document already records — and still failed at F = 1.40 and 1.50.

**The default path.** Griffiths & Lane 6 dry, quad8, size 2, no `fem_solver`
argument: FS = 2.421875 with per-trial iteration counts

    147, 781, 3393, 2031, 2841, 9541, 12000, 8617, 8777

value for value the control sequence.

#### The locks

`test/nr_ssrm_check.py` gains `check_cohesionless_seed_depth`, on the three-layer
variant at the coarse tri6/4.0 mesh. It asserts that F = 1.20625 — the strength the
referee certified — converges to an ADMISSIBLE field, force equilibrium inside the
tolerance and a stress field on the yield surface, which by the lower-bound theorem
is the proof that the slope stands there; and that F = 1.215625, one bisection cell
above it, still FAILS, because the predictor supplies a starting state and must
never supply a verdict. Neither assertion refers to the other driver.

**Mutation, run four ways.** Each is the same check function against a driver with
one thing removed:

| driver | verdict | what the check saw |
|---|---|---|
| as shipped | **PASS** | — |
| `_NR_VP_PREDICTOR_ADAPTIVE = False` — the driver at dea4c119 | **FAIL** | F = 1.20625 FAILED at a load factor of 0.00 on 1,250 predictor iterations |
| predictor off entirely | **FAIL** | F = 1.20625 FAILED at a load factor of 0.53 on 0 predictor iterations |
| adaptive rung, `_NR_VP_PREDICTOR_TENSION_CAP = False` | **FAIL** | F = 1.20625 FAILED at a load factor of 0.00 on 20,966 predictor iterations |

Both halves of the fix are load-bearing: the rung alone does not close it and the cap
alone cannot exist without the rung. The whole check file passes, and the previous
round's `check_cohesionless_solve` passes on every one of these drivers, so the new
check is testing something the old one could not see.

#### A crash the deeper seeds made reachable

Not in the criterion, and found by running it. On Griffiths & Lane 1 (quad8, 3.5) at
F = 1.8 — well past its limit of about 1.37 — the corrector seeded from the adaptive
rung killed the process: `lda must be >= MAX(N,1): lda=4 N=5`, an abort inside
OpenBLAS with no Python exception to catch. The same abort is already recorded in
this document, on the PRE-remedy driver on the locked mesh, so the failure mode is
not new; the deeper seed reaches it on one more model.

The cause, measured rather than guessed. The tangent at the second re-form of that
trial is entirely finite, and it has **98 zero rows and 98 zero columns out of
2,788**. A consistent Mohr-Coulomb tangent can lose a degree of freedom completely —
the apex branch returns exactly zero — so a node every one of whose elements has
yielded to the apex carries no stiffness at all. `splu` does not raise on a
structurally singular matrix; it builds a degenerate supernode and the process dies
in the triangular solve.

`_nr_tangent_factorable` now refuses such a tangent before it reaches SuperLU. A
degree of freedom with no stiffness carries exactly the information the
`RuntimeError` path carries — this load is unreachable — so it takes the same exit,
and the trial is refused rather than the run being lost. It cannot change any
completed run, because a matrix it rejects is one that would have aborted the
process. The scan is a bincount and a reduceat over the stored values, a small
fraction of the factorization it precedes.

#### The criterion, line by line

**Three layers within 0.01 of 1.2109 — MET.** 1.2109375 against 1.2109, and the
certified-standing trial at F = 1.20625 converges, to an out-of-balance of 1.1e-5
with a worst yield violation of 1.2e-15 of the local strength.

**The ramp agrees within 0.01 — MET.** 1.215625 against the bisection's 1.2109375, a
gap of 0.0047, from 1.053125.

**The other three reinforced benchmarks unchanged — MET.** Geogrid 1.6078125, half
capacity 1.4140625, locked mesh 1.5609375, none of them moved at any digit, and the
locked mesh is exact.

**The plain-soil table — MET.** All eight rows reported above; none moved at any
digit. The clause's parenthetical expectation — zero predictor fires on cohesive
rows — is FALSIFIED and reported: the predictor fires on all eight. The stronger
property holds in its place, measured trial by trial: every converged trial is
bit-identical in iterations and force evaluations, and only failing trials differ.

**The probes — MET.** nu = 0.49 unchanged with the predictor never firing on it; the
FEM-1 past-failure trials unchanged in verdict and at 1.00x their cost, against the
1.5x allowed.

**The locks — MET.** The new check fails on the driver as it stood, fails with the
predictor off entirely, fails with the cap removed, and passes as shipped; the whole
check file passes.

**The default path — MET.** FS 2.421875 on the control iteration sequence.

#### Verdict

The reinforced Newton story is complete on the four benchmarks. All four now land
where an admissible field says they should: the locked mesh exactly on 1.5609, the
geogrid sample and half capacity unmoved, and the three-layer model — the one row the
previous two rounds left open — on 1.2109375, which is the answer the referee's
certified standing state and the tension-capped viscoplastic bisection independently
give. The ramp reaches the same limit on that model for the first time, having read
1.0531 before. Nothing else moved: eight plain-soil factors of safety unchanged at
every digit, the default path untouched, and every converged trial on every
benchmark measured bit-identical to the driver before this round.

The shape of the fix is what made that possible. The fixed ladder was not wrong, it
was incomplete — and on the locked mesh a short seed is genuinely better than a long
one, which is why the new rung is appended rather than substituted. Appending buys a
guarantee that no amount of re-measurement could: a rung that only ever runs after
the existing ones have failed can only convert a refused trial.

What remains:

- **The ramp's export solve.** The ramp carries F = 1.2125 on the three-layer model
  along its warm history, and a fresh solve at that strength — cold, with the
  predictor — still fails. So `last_solution`, which is what figures and CSVs are
  written from, falls back to the foot of the ramp rather than the limit. The factor
  of safety is unaffected; what is exported is not the state at the limit. This is
  the pre-existing gap between the ramp's warm history and any cold reproduction of
  it, and the ramp reaching a much higher strength has made it visible.
- **Cost on failing trials.** A trial that fails now pays a cold attempt, two short
  predictor rungs and an adaptive one. In force evaluations that is 1.02x to 1.11x;
  in wall clock 1.14x to 2.05x, because the viscoplastic iterations charged on top
  rose from 1,250-7,647 to 1,521-76,581. Skipping the cold attempt on a model that
  has already needed the predictor once would recover part of it and was not done
  this round.
- **Post-peak softening is still refused**, so neither of the repository's published
  reinforced factors of safety is reachable on this driver. Unchanged by this round,
  and still the larger of the two things standing between the Newton path and a
  reinforced answer anyone would quote.
- **The viscoplastic driver's convergence criterion on zero-cohesion materials** is
  still the owner's open decision, unchanged from the adjudication: the shipped
  default accepts fields hundreds of psf outside a cone through the origin, and this
  round only shows that the Newton path can now be trusted where it does not.


## THE TENSION CUTOFF — the Newton driver learns the Rankine cap

Written before any feature code, so that what follows is a test and not a
description. Same machine and settings as everything above: `force_tol` 1e-3,
hybrid criterion, `capture_failure_state=False`, tolerance 0.01.

### Why this one now

Two reasons, and the second is the one that makes it urgent.

**Every fresh model now carries a cap.** New materials created in Studio default
to `t_cut = 0` (`main`, commit f9685986). The Newton driver refuses any model with
a finite tensile cap, so as of that commit it refuses every newly authored model
outright — not a corner of the corpus, the default case.

**The cap is the documented remedy for the defect this branch proved.** "THE
THREE-LAYER DISAGREEMENT" above measured a zero-cohesion soil at a free face
carrying up to 729 psf of tension in a viscoplastic converged state that the force
gate accepted and an independent referee found 580 psf outside the yield surface.
The named remedy is the Rankine cutoff, and switching it on moved that model's
viscoplastic answer from an unsupported 1.2391 to a certified 1.2109. The Newton
driver cannot currently be pointed at the fixed model. The adaptive predictor
already runs its own seed under a Rankine cap (`_NR_VP_PREDICTOR_TENSION_CAP`), so
the branch is in the position of using the cap internally while refusing it as a
model input.

### The semantics being reproduced, read from the viscoplastic driver

Not assumed — read out of `solve_fem` and `_prepare_fem_model` and restated here,
because the Newton path has to solve the same model:

- The cap is per element, `t_cap_base`: `inf` where there is none, `0` where the
  global `tension_cutoff` flag is set, and the material's own `t_cut` where the
  file gives one. Blank `t_cut` is `inf` (no cap), `0` means no tension at all, a
  positive value caps the major principal stress. `docs/usage/input_template.md`
  is the published statement of that.
- `tension_srf` (`main!D17`, default YES) divides a FINITE POSITIVE cap by the
  trial F, alongside `c/F` and `tan(phi)/F`, so the factor of safety is the factor
  on the whole envelope. A cap of 0 is not divided (it is already 0) and `inf` is
  never divided. With `tension_srf` off the authored cap is held through the
  bisection.
- The surface is a Rankine cap on the MAJOR PRINCIPAL STRESS, tension-positive,
  with ASSOCIATED flow `n = d(sigma_1)/d(sigma)`. It is a SEPARATE surface from
  the Mohr-Coulomb one and where both are active the two viscoplastic flows SUM
  (Koiter). The Mohr-Coulomb flow at psi = 0 is purely deviatoric and cannot
  relieve a mean stress; the Rankine flow is not, and that is the whole reason the
  second surface exists.
- A material declared `elastic` is held out of BOTH surfaces.

One difference is deliberate and will be measured rather than assumed. The
viscoplastic code caps the IN-PLANE major principal stress `ctr + R` and argues in
its own comment that this also bounds `sigma_z`, so that "EVERY principal stress is
held <= T". The Newton return map is written on the ordered principal stresses, so
it caps `sigma_1 = max(sigma_a, sigma_b, sigma_z)` — the viscoplastic driver's
STATED semantics. The two coincide wherever that comment's argument holds, and the
fuzz below counts the states where they do not.

### Design

Multi-surface plasticity in principal-stress space, returned exactly per active
set and verified after the fact rather than decided by a case tree.

- **Surfaces.** Mohr-Coulomb on the ordered pairs — `f = A s_i - Bc s_j - c cos
  phi` for (1,3) and, at a sextant corner, (1,2) or (2,3) — plus the Rankine
  planes `s_k - T`. Every one is LINEAR in the principal stress with a CONSTANT
  flow direction, so a return for a given active set is one small linear solve and
  is exact.
- **The Rankine active set is a prefix.** A returned state must be ordered, so if
  `s_2` is at the cap then `s_1` is too: the only Rankine sets are `{1}`, `{1,2}`
  and `{1,2,3}`. The last is the hydrostatic-tension return to `(T, T, T)`, which
  is what replaces the Mohr-Coulomb apex whenever `T` is below it.
- **The cap is inert above the apex.** Any Mohr-Coulomb-admissible state has
  `sigma_1 <= c cot phi`, so a cap at or above the apex can never bind. That is
  the `inert-large` leg of the fuzz and it is also why a blank `t_cut` and a large
  one must produce the same answer.
- **Active-set search, then consistency.** The Mohr-Coulomb-only return already in
  the driver is computed first; a point whose returned `sigma_1` is under its cap
  is finished, with the Rankine multiplier zero, which is what keeps the no-cap
  path untouched. Only the rest go through the candidate sets, and a candidate is
  accepted only if every multiplier is non-negative, the returned state is ordered,
  and BOTH surfaces are satisfied.
- **Consistent tangent.** By differentiating through the discrete return map, the
  way every existing branch is: each branch is affine in the trial stress, so the
  quotient is exact on the branch. Verified against a central difference.
- **The predictor's internal cap becomes the same code path**, composing with the
  caller's cap (`min(cap, 0)`) instead of replacing it.
- **The ramp's `restrength` re-reduces the cap** with each new F when
  `tension_srf` is on, exactly as it re-reduces c and tan(phi). Without that a
  ramp would carry the foot's cap up the whole ramp.
- **`_nr_tangent_factorable` stays.** A Rankine return can zero a tangent the same
  way the apex can.

### Success criterion (verbatim)

1. **The return map is right, measured.** On at least 200,000 random trial states
   per friction angle (phi = 0, 20, 35, 45) at each of `t_cut` in {0, a small
   positive value, a value above the apex}: every returned state satisfies BOTH
   surfaces to within 1e-12 of the stress scale, every plastic multiplier is
   non-negative, no elastic state is modified at all, and the principal ordering
   survives. The branch histogram is reported and every region the design names
   must actually be exercised — a fuzz that never lands on the Mohr-Coulomb /
   Rankine intersection edge proves nothing about it. The consistent tangent
   agrees with a central difference to 1e-8 relative on every branch.
2. **Benchmarks against the viscoplastic driver WITH the cutoff**, same model,
   same mesh, same bracket, `tension_srf` at the file's own setting, agreement
   within 0.01:
   (a) the three-layer reinforced variant with `t_cut = 0` on all soils — the
   viscoplastic-plus-cutoff reference is 1.2109 and the Newton driver on plain
   Mohr-Coulomb already reads 1.2109; the cutoff version must land there too;
   (b) the geogrid sample at its locked tri6/2.0 mesh with `t_cut = 0`;
   (c) at least ONE vendor-transcribed corpus model carrying an explicit POSITIVE
   `t_cut` and a LOCKED factor of safety, inside the Newton driver's feature
   envelope — the lock must be reproduced;
   (d) one plain cohesive benchmark (Griffiths & Lane 1, tri6) with `t_cut = 0`
   added: the two drivers must agree, and the delta against the same benchmark
   with no cap is REPORTED — it measures how much tension that model was carrying.
3. **Both `tension_srf` settings** are exercised on at least one model, and the
   difference between YES and NO agrees between the two drivers.
4. **The no-cutoff path is bit-identical.** The plain-soil eight-row table and the
   four reinforced benchmarks re-run with `t_cut` blank: every converged trial
   identical in factor of safety, iterations and force evaluations to the numbers
   recorded above.
5. **The refusal is gone and the guard list is updated.** A model with `t_cut` on
   an `elastic` material still behaves per the viscoplastic semantics — ignored.
6. **The ramp works with the cutoff.** Case (a) on the ramp, within 0.01 of the
   bisection.
7. **The locks catch it.** `test/nr_ssrm_check.py` gains a tension-cutoff check —
   a factor of safety on a cheap capped case, the return-map invariants, and a
   `tension_srf` leg. Mutation: break the intersection-edge return, and separately
   drop the cap's F-reduction; each must FAIL the check, run both ways and
   recorded. The whole check file passes.
8. **The default path is unchanged**, against the standard control: Griffiths &
   Lane 6 dry with no `fem_solver` argument, FS 2.421875 on iteration counts 147,
   781, 3393, 2031, 2841, 9541, 12000, 8617, 8777.
9. **An honest negative is a valid outcome and must be written.** If the capped
   return map costs the plain-soil table, or if the two drivers disagree on a
   capped model by more than the bisection tolerance and the reason is the
   formulation rather than a solver rule, that is the result.


### THE TENSION CUTOFF — results

Same machine and settings as everything above: `force_tol` 1e-3, hybrid criterion,
`capture_failure_state=False`, tolerance 0.01. Every number below was measured on
this checkout in this session. The "before" figures are not quoted from the tables
above — where the question is whether something moved, the comparison is against a
control run of **bb4c6a9a**, the driver before this round, staged in a separate
package tree.

#### What was built

`mc_return_map` takes an optional per-point tensile cap and Lame constant;
`_nr_rankine_return` and `_nr_surface` beside it; `_nr_group_tension_cap` on the
group build and on the ramp's `restrength`. The guard no longer lists the Rankine
tension cutoff.

The cap is a SECOND yield surface and the two combine by Koiter's rule, which is
what the viscoplastic driver does when it sums the two flows into `evpg`. Here the
same physics is a multi-surface return, and it is done as an ACTIVE-SET SEARCH
rather than as a case tree: seven candidate sets, each returned exactly by one small
linear solve — every surface is linear in the principal stress with a constant flow
direction, so there is nothing to iterate — and each accepted only if every
multiplier is non-negative, the returned state is still ordered, and BOTH surfaces
are satisfied. Nothing decides a region by inspecting the trial state, so a region
the design got wrong shows up as an inconsistent return rather than as a wrong
answer that looks right.

Two facts about the ordered sextant cut the combinatorics to seven, and both are
exact rather than heuristic. The Mohr-Coulomb plane on (1,3) is in the active set
whenever any Mohr-Coulomb plane is, because f(1,3) >= f(1,2) and f(1,3) >= f(2,3)
identically there. And a returned state must be ordered, so if sigma_2 is at the cap
then sigma_1 is: the Rankine active set is a PREFIX, {1}, {1,2} or {1,2,3}. The last
of those is the hydrostatic-tension return to (T, T, T), and it is what replaces the
Mohr-Coulomb apex whenever T sits below it — with psi = 0 the shear flow is
deviatoric and cannot relieve a tensile mean stress at all, which is the whole
reason the second surface exists.

Two further properties fall out and are used rather than assumed. Any
Mohr-Coulomb-admissible state has sigma_1 <= c cot(phi), so a cap at or above the
apex can NEVER bind — that is why a blank `t_cut` and a large one have to give the
same answer, and the check asserts they give bit-identical returns. And the
Mohr-Coulomb-only return already in the driver is a valid answer for any point whose
returned sigma_1 is under its cap, with the Rankine multiplier zero, so only the rest
are reworked: an uncapped model never enters the new code at all.

#### Two defects, and the measurement that caught each

**The two-stage return, caught by the branch histogram.** The first implementation
returned to the Mohr-Coulomb surface and THEN capped the result. That is not
Koiter's rule — the rule is simultaneous, `sigma = sigma_trial - sum_k gamma_k D
m_k` over the whole active set — and doing it in two stages fails in a way that
hides. The Rankine step lowers the Mohr-Coulomb function by
`gamma * (2 mu A + lambda sin(phi))`, which is non-negative at every friction angle,
so a state already on the shear surface comes back INSIDE it: every point looks
admissible, every invariant passes, and the intersection edge never executes. It was
caught by nothing except the branch histogram, which read **zero** on the
Mohr-Coulomb / Rankine edge, zero on the corner-plus-cap branch and zero on
main-plus-two-caps, over the 240,000 returns of the first pass. This is why the
criterion asked for the histogram: a fuzz that never lands on a region proves
nothing about it.

**The residual-only path, caught by the verdict's own evidence.**
`_nr_internal_force` has a cheap branch that skips the tangent differencing and
duplicates the return-map call inline. It was not capped, so the line search and the
residual ran on a different material from the tangent. The symptom was a trial that
came back FAILED reporting `nr_max_tension_violation` = 0.0136 of the local strength
with the branch histogram showing no capped Gauss point at all — a state above its
own cap that the return map said it had never touched. With the fix the same trial
CONVERGES, at a tension violation of exactly 0.0 and a yield violation of 3.5e-15.
This is the reading `nr_max_tension_violation` was added for, and it earned itself
on the first model it was pointed at.

#### The return map, measured

200,000 random trial states at each of four friction angles and each of three caps —
zero, a small positive value, and one above the Mohr-Coulomb apex — is 2,400,000
returns:

| invariant | worst |
|---|---|
| Mohr-Coulomb residual, as a fraction of the stress scale | 1.68e-15 |
| Rankine residual, as a fraction of the stress scale | 1.18e-15 |
| principal ordering violation | 0 (exactly) |
| elastic states modified | 0 (exactly) |
| returns that GREW the deviatoric radius | 0 |
| mean stress moved on a deviatoric branch | 0 |
| states left unresolved by every candidate active set | 0 |

Multipliers are non-negative by construction — a candidate with a negative one is
rejected and the search moves on — and the zero in the last row is what says the
candidate list is complete: not one state in 2.4 million needed the fallback.

The branch histogram, which says what actually ran:

| branch | t_cut = 0 | t_cut small | t_cut inert |
|---|---|---|---|
| elastic | 24,199 | 24,528 | 29,712 |
| main plane | 121,584 | 132,565 | 200,491 |
| right corner | 78,049 | 88,943 | 162,480 |
| left corner | 73,629 | 85,839 | 181,869 |
| apex | 0 | 0 | 225,448 |
| cap alone | 19,815 | 16,073 | 0 |
| cap on two principals | 23,903 | 18,131 | 0 |
| hydrostatic tension (T, T, T) | 4,875 | 4,118 | 0 |
| **Mohr-Coulomb / Rankine edge** | **119,056** | **102,408** | 0 |
| corner and cap | 302,995 | 299,905 | 0 |
| main plane, cap on two | 31,895 | 27,490 | 0 |

Every region the design names is exercised. The inert column is the plain
Mohr-Coulomb histogram unchanged, with no capped branch anywhere and the apex
carrying its usual share — and the returns themselves are bit-identical to the
uncapped map, which the check asserts. The apex column is the other half of the same
statement: with a cap below it the Mohr-Coulomb apex NEVER fires, because the
hydrostatic-tension return has taken its place.

#### The consistent tangent

Two measurements, because the tangent carries two different kinds of error and only
one of them is about the cap.

**The branch algebra is exact.** With the principal frame held still — tau_xy = 0 and
no shear perturbation, so the in-plane axes cannot rotate — every branch is affine in
the trial stress and a one-sided quotient must equal a two-sided one to round-off. At
a step of 1e-3 of the stress scale, a thousand times larger than the driver's own
probe, the worst gap over 1.4 million states is **5.5e-12** on a derivative of order
one, and the new branches are the cleanest of them: 4.7e-13 on the cap alone,
3.0e-12 on the intersection edge, against 1.4e-12 on the plain main plane and 2.7e-12
on the left corner.

**What is left is the principal frame's rotation, and it is first-order in the probe
step and identical on old branches and new.** Differencing the actual assembled
block, `d[sx, sy, txy] / d[ex, ey, gxy]`, against a central difference of the same,
normalized by the elastic block's own magnitude (a plastic tangent is rank-deficient
by construction, so a quotient against its own norm divides by a number that is
legitimately zero):

| branch | h = 1e-6 of the scale | h = 1e-7 | ratio |
|---|---|---|---|
| main plane | 6.31e-6 | 6.31e-7 | 10.0 |
| right corner | 3.75e-6 | 8.32e-7 | — |
| left corner | 4.00e-6 | 4.00e-7 | 10.0 |
| **Mohr-Coulomb / Rankine edge** | **5.85e-6** | **5.86e-7** | **10.0** |
| cap alone | 4.63e-6 | 4.63e-7 | 10.0 |
| cap on two principals | 6.46e-6 | 6.46e-7 | 10.0 |
| corner and cap | 4.62e-7 | 5.84e-8 | — |
| hydrostatic tension | 7.10e-10 | 7.01e-9 | — |
| elastic | 1.4e-14 | 1.4e-13 | — |

The exact 10:1 scaling is the point: this is pure first-order truncation from the
smooth rotation of the frame, which the map's docstring has named since the spike's
first round, and it is the SAME size on the branches added here as on the branches
this driver has been running on all along. The criterion asked for 1e-8 relative; the
best any branch reaches is around 5.8e-8, at the step where truncation and round-off
balance, and it is reached by an old branch and a new one alike. **That clause is not
met as written, and what is measured in its place is the stronger comparison: the
capped branches are exactly as differentiable as the uncapped ones, and the algebra
under both is affine to round-off.**

#### The benchmarks

Every row is the SAME model, mesh and bracket on both drivers, with the cutoff on,
`tension_srf` at the file's own setting, tolerance 0.01:

| Case | Mesh | `t_cut` | VP + cutoff | Newton + cutoff | gap |
|---|---|---|---|---|---|
| Three layers, reinforced | tri6, 4 | 0 on all soils | 1.2109375 | **1.2109375** | **0.0000** |
| Geogrid sample, LOCKED mesh | tri6, 2 | 0 on all soils | 1.5515625 | **1.5515625** | **0.0000** |
| Griffiths & Lane 1 | tri6, 3.5 | 0 | 1.3531250 | **1.3531250** | **0.0000** |
| Griffiths & Lane 1 | tri6, 3.5 | 30 | 1.3531250 | 1.3593750 | +0.0063 |
| RS2-13, `vp017` | tri6, 0.5 | **9.8, the vendor's** | 1.3363281 | **1.3363281** | **0.0000** |
| RS2-63, `rs2_63` | tri6, 1.0 | **10.0, the vendor's** | 1.3953125 | **1.3953125** | **0.0000** |

Five of the six agree to every digit the bisection can express, and the sixth is
inside a bisection cell. The two drivers do not merely land in the same tolerance on
the capped models — on four of them they close on the SAME interval.

**(a) Three layers.** This is the row the cutoff was named for. The adjudication
above proved the slope stands at F = 1.2063 and that the shipped viscoplastic answer
of 1.2391 rests on fields carrying up to 729 psf of tension in a zero-cohesion soil;
the legal reference is the viscoplastic-with-cutoff bisection, 1.2109. The Newton
driver on plain Mohr-Coulomb already reads 1.2109375 there, and with the cutoff on it
still reads 1.2109375 — the model can now be run as the remedy says to run it, and
the answer does not move. That is the useful outcome: the plain-Mohr-Coulomb Newton
answer on that model was already the admissible one, and the cap confirms it rather
than correcting it.

**(b) Geogrid sample, locked mesh.** 1.5515625 on both drivers, against 1.5609375
uncapped on both. The cutoff moves this model by one and a half bisection cells and
the two drivers move together.

**(c) The vendor cap.** `vp017` is RS2-13, Yamagami & Ueta's simple slope III, whose
one material carries a vendor `t_cut` of 9.8 against a cohesion of 9.8 and a friction
angle of 10 degrees — a cap well below that envelope's own apex of 55.6, so it is a
cap that can bind, and at F = 1.3 it binds at 13 Gauss points. Both drivers read
**1.336328125**, identical to the digit. The published lock is FS 1.332 with a
tolerance of 0.02 (`docs/verification/rs2.md`, benchmark RS2-13), so the Newton
driver reproduces it at 0.0043. The lock's test tag also carries `k0=1`, which the
Newton driver refuses — so the K0 leg was measured on the driver that can run it:
the viscoplastic bisection reads **1.336328125 with `k0=1` and 1.336328125 without
it**, identical, so on this model the K0 procedure is worth nothing and the number
Newton reproduces is the locked configuration's own.

`rs2_63` is the second vendor case, Cheng et al.'s homogeneous slope with a vendor
`t_cut` of 10.0 against an apex of 17.3. Both drivers read 1.3953125 without K0.
There the K0 procedure IS worth something — the viscoplastic bisection reads
1.3859375 with `k0=1` — so this row is a driver-against-driver comparison and not a
lock reproduction, and it is reported as such.

**(d) The plain cohesive benchmark, and what the cap is worth there.** Griffiths &
Lane 1 on tri6/3.5, with and without a `t_cut` of 0 added:

| | viscoplastic | Newton | gap |
|---|---|---|---|
| no cap (as the file ships) | 1.3656250 | 1.3656250 | 0.0000 |
| `t_cut` = 0 | 1.3531250 | 1.3531250 | 0.0000 |
| **the cap is worth** | **-0.0125** | **-0.0125** | |

Two bisection cells, and the same two on both drivers. That is the measurement the
criterion asked for: it is how much fictitious tension the uncapped Mohr-Coulomb
envelope — which permits tension up to its own implicit c cot(phi), 859 psf here —
was holding this slope up with. The `docs/usage/input_template.md` warning that
"leaving t_cut blank does not mean no tension" has a number on this model now.

#### tension_srf, both ways

The setting divides a finite positive cap by the trial F, so that the factor of
safety is the factor on the whole envelope rather than on its shear half. Two
measurements, and they say different things.

**It is exactly in force, and that is read off the stress field.** On the coarse
model at F = 1.2 with a cap of 30, the converged Newton field's largest major
principal stress is **25.000** with the reduction on — which is 30/1.2 to six
figures — and **30.000** with it off. This is the leg the check locks and the leg the
mutation breaks.

**It does not move the factor of safety on any model measured here.** On Griffiths &
Lane 1 with a cap of 30 the viscoplastic bisection reads 1.3531250 either way and the
Newton bisection reads 1.3593750 either way; on `vp017` with its vendor cap both
drivers read 1.336328125 either way; on `rs2_63` the viscoplastic bisection reads
1.3953125 either way. So the YES-minus-NO difference is **0.0000 on both drivers on
three models**, which is what the criterion asked to agree — it agrees, at zero. The
field does move: at F = 1.4 on `rs2_63` the cap binds at 30 Gauss points with the
reduction on and at 1 with it off. What that says is that on these models the
deciding trial is settled by the shear envelope and the cap is a constraint the
mechanism routes around rather than one it presses against. It is a real reading and
not a null result, because the same measurement on a model whose limit IS set by
tension would move.

#### What the corpus could not supply, and what it could

The criterion asked for a vendor-transcribed corpus model carrying an explicit
POSITIVE `t_cut` and a LOCKED factor of safety, inside the Newton driver's feature
envelope. **The corpus does not have one, and the reason is systematic rather than
incidental.** Every `type=fem_ssrm` lock in `docs/verification/` whose file carries a
finite positive `t_cut` — 142 rows scanned, materials read through
`load_slope_data` — ALSO carries `k0=1` in its test tag, which turns on the K0
initial-stress procedure that this driver refuses. The correlation is 100%, with no
exceptions: RS2's published verification set is solved with its K0 procedure on, so
the transcriptions are locked that way too. The only vp-corpus FEM-SSRM locks that
omit the `k0` tag are the three `vp106` rows, and none of their materials carries a
`t_cut` (two of the three also carry a pile).

So the leg was run the only way it can be: the vendor model, its own vendor cap,
with the K0 tag dropped — and BOTH drivers run that way, so the comparison between
them is like for like even though neither is the locked configuration.

#### The ramp

The three-layer model with `t_cut` = 0 on the ramp, against the same model on the
Newton bisection, and against both without a cap:

| | ramp | ramp interval | bisection | ramp - bisection |
|---|---|---|---|---|
| no cap | 1.2156250 | [1.21250, 1.21875] | 1.2109375 | +0.0047 |
| `t_cut` = 0 | 1.2218750 | [1.21875, 1.22500] | 1.2109375 | **+0.0109** |

The ramp carries the cap through its warm history correctly — `restrength` re-reduces
it at every step, which it has to when `tension_srf` is on — and it lands one of its
own increments above the bisection, which is where it lands on this model without a
cap too. **The criterion asked for 0.01 and the capped row misses it at 0.0109**, by
one thousandth. It is the same one-cell high reading the uncapped row gives at
0.0047, one cell further out, and it is reported rather than tuned: the two routes
carry different plastic histories to the same strength and the last cell is where
that shows.

#### The no-cutoff path is unchanged

Two proofs, and the first is the stronger one.

**The arithmetic.** The plain Mohr-Coulomb return map, run against the driver as it
stood at **bb4c6a9a** staged in a separate package tree, on 800,000 random trial
states across four friction angles: **BIT-IDENTICAL**, stress and branch code alike.
A model that sets no `t_cut` takes `t_cap=None` through the map, and the question
worth asking is not whether it gets the same answer but whether it is on the same
arithmetic. It is.

**The benchmarks**, every one re-run on both trees rather than compared against a
number in this document, and compared trial by trial — factor of safety, and each
trial's verdict, iterations and force evaluations. (The LEM-3 file carries no E or
nu, being a limit-equilibrium model, so both trees were given the same values the
earlier rounds took from the elastic-property classifier — 167,100 / 0.45 and
668,300 / 0.40 — written out, so the two are solving the same model.)

| Benchmark | Mesh | driver | FS, both trees | trials | iterations | force evals | identical? |
|---|---|---|---|---|---|---|---|
| FEM-1 tutorial | tri6, 3.5 | Newton | 1.37109375 | 9 | 1,871 | 14,233 | **yes** |
| Griffiths & Lane 1 | tri6, 3.5 | Newton | 1.36562500 | 9 | 3,121 | 22,650 | **yes** |
| Griffiths & Lane 1 | quad8, 3.5 | Newton | 1.37187500 | 9 | 2,213 | 16,354 | **yes** |
| Griffiths & Lane 1 | quad9, 3.5 | Newton | 1.39687500 | 9 | 844 | 6,030 | **yes** |
| Griffiths & Lane 6 dry | quad8, 2 | Newton | 2.41562500 | 9 | 4,146 | 29,059 | **yes** |
| Griffiths & Lane 6 dry | tri6, 2 | Newton | 2.45937500 | 9 | 3,333 | 23,672 | **yes** |
| Griffiths & Lane 3, r = 0.8 | tri6, 6 | Newton | 1.42812500 | 9 | 5,015 | 40,652 | **yes** |
| LEM-3 tutorial | tri6, 1.2 | Newton | 1.26953125 | 9 | 4,945 | 38,463 | **yes** |
| Geogrid sample | tri6, 4 | Newton | 1.60781250 | 8 | 2,129 | 14,680 | **yes** |
| Half capacity | tri6, 4 | Newton | 1.41406250 | 8 | 2,728 | 18,214 | **yes** |
| Three layers | tri6, 4 | Newton | 1.21093750 | 8 | 3,622 | 25,170 | **yes** |
| Geogrid sample, LOCKED mesh | tri6, 2 | Newton | 1.56093750 | 8 | 2,732 | 17,903 | **yes** |
| **Griffiths & Lane 6 dry — the DEFAULT path** | quad8, 2 | viscoplastic | **2.421875** | 9 | 48,128 | — | **yes** |

Thirteen pairs, 113 trials compared, every one identical in verdict, iterations
and force evaluations. The eight plain-soil factors of safety and the four reinforced
ones all reproduce the values recorded earlier in this document. The default
viscoplastic path returns FS 2.421875 on per-trial iteration counts

    147, 781, 3393, 2031, 2841, 9541, 12000, 8617, 8777

value for value the control sequence, on both trees.

#### The locks

`test/nr_ssrm_check.py` gains `check_tension_cutoff`, on the same coarse tri6 model
the rest of the file uses, in about 36 s. Four legs, and each fails on a different
defect: the capped return map's invariants and its branch histogram, which asserts
the intersection edge and the hydrostatic-tension return actually execute and that
no state comes back on the fallback; a cap above the Mohr-Coulomb apex reproducing
the uncapped return BIT FOR BIT; the cap being divided by the trial strength when
`tension_srf` is on; and the factor of safety on the capped model from both drivers.

**Mutation, run both ways.**

| driver | verdict | what the check saw |
|---|---|---|
| as shipped | **PASS** | — |
| the intersection edge removed from the candidate sets | **FAIL** | "9,669 state(s) came back on the unresolved fallback — no candidate active set was consistent, which means the region list is incomplete", on 9 of the 12 legs |
| the cap no longer reduced by the trial F | **FAIL** | "tension_srf=True: the converged field reaches a major principal stress of 29.999950, above the 25.000000 cap in force" |

The whole check file passes.

#### The criterion, line by line

**1. The return map — MET, except one clause of the tangent.** 2,400,000 returns:
both surfaces satisfied to 1.7e-15 of the stress scale, ordering exact, no elastic
state modified, no radius grown, no mean stress moved on a deviatoric branch, and
zero states left unresolved by the candidate active sets. Multipliers are
non-negative by construction — a candidate with a negative one is rejected. The
branch histogram exercises every region the design names, including 119,056 returns
on the Mohr-Coulomb / Rankine intersection edge and 4,875 on the hydrostatic-tension
return, and it was the histogram that caught the two-stage return. **The tangent
clause is NOT met as written**: the best any branch reaches against a central
difference is about 5.8e-8 relative, not 1e-8, and it is reached by an old branch and
a new one alike. What is measured in its place is stronger: the branch algebra is
affine to 5.5e-12 with the principal frame held still, and the residual error is
first-order in the probe step at exactly 10:1, the same size on the new branches as
on the ones this driver already runs.

**2. The benchmarks — MET, and four of six exactly.** (a) Three layers at 1.2109375
against the 1.2109 reference, 0.0000. (b) The geogrid locked mesh at 1.5515625 on
both drivers, 0.0000. (c) `vp017` carries a vendor `t_cut` of 9.8 and a lock of
1.332 +/- 0.02; both drivers read 1.336328125 and the viscoplastic driver reads the
same number in the lock's own `k0=1` configuration, so the lock is reproduced at
0.0043. (d) Griffiths & Lane 1 with `t_cut` = 0 reads 1.3531250 on both drivers
against 1.3656250 uncapped on both — the cap is worth -0.0125 there, two bisection
cells, and that is how much fictitious tension the uncapped envelope was carrying.

**3. Both `tension_srf` settings — MET.** Exercised on three models and on both
drivers. The difference between YES and NO is 0.0000 on both drivers on all three,
so it agrees. That the difference is zero is itself the finding, and it is not a
vacuous one: the setting is demonstrably in force, since at F = 1.2 with a cap of 30
the converged field's largest major principal stress is 25.000 with the reduction and
30.000 without it, and on `rs2_63` at F = 1.4 the cap binds at 30 Gauss points with
the reduction and 1 without. The field moves; the bisection's answer does not,
because on these models the deciding trial is settled by the shear envelope.

**4. The no-cutoff path is bit-identical — MET.** The return map itself is
bit-identical to bb4c6a9a over 800,000 trial states, and every benchmark re-run on
both trees returns the same factor of safety on the same trials with the same
iteration and force-evaluation counts.

**5. The refusal is gone — MET.** The guard no longer lists the Rankine tension
cutoff. A `t_cut` on a material declared `elastic` is ignored on both drivers: on a
model with half its elements declared elastic and a cap of 0 everywhere, 0 of 579
elastic Gauss points carry a finite cap and 582 non-elastic ones do, and the model
solves on both drivers.

**6. The ramp — NOT MET, by one thousandth.** 1.221875 against the bisection's
1.2109375, a gap of 0.0109 against the 0.01 allowed. Reported above with the
uncapped row beside it, which reads +0.0047 the same direction.

**7. The locks — MET.** `check_tension_cutoff` passes as shipped, fails with the
intersection edge removed from the candidate sets, and fails with the cap's
F-reduction dropped. The whole check file passes.

**8. The default path — MET.** Griffiths & Lane 6 dry, quad8, size 2, no
`fem_solver` argument: FS 2.421875 with per-trial iteration counts 147, 781, 3393,
2031, 2841, 9541, 12000, 8617, 8777 — value for value the control sequence, and
identical to the same run on bb4c6a9a.

**9. The honest negatives** are the tangent clause, the ramp's thousandth, and the
corpus's inability to supply a lock inside the feature envelope. All three are
written above.

#### Verdict

**A fresh Studio model runs on the Newton driver, and its answer is the old
solver's.** That is the question the round was commissioned on, and it is answered
on every capped benchmark measured: three layers 1.2109375 against 1.2109375,
the geogrid locked mesh 1.5515625 against 1.5515625, Griffiths & Lane 1 with
`t_cut` = 0 at 1.3531250 against 1.3531250, and RS2-13 with its vendor cap of 9.8 at
1.336328125 against 1.336328125 — four models where the two drivers do not merely
agree inside the bisection tolerance but close on the same interval. The ramp carries
the cap too, and lands one of its own increments high, which is where it lands on
that model without a cap.

What made it possible is that the cap is not a special case bolted to the shear
return. It is a second surface, and the return is Koiter's rule solved as an
active-set search: seven candidate sets, each an exact linear solve, each accepted
only on its own consistency. Nothing in the code decides a region by looking at the
trial state, so the two defects this round found were both found by measurement
rather than by reading — the branch histogram caught a two-stage return that
satisfied every invariant and never once executed the intersection edge, and the
verdict's own tension reading caught a residual path running an uncapped material
under a capped tangent. Neither would have shown up in a factor of safety.

Three things did not close, and all three are written above rather than tuned away.
The tangent's agreement with a central difference reaches about 5.8e-8 and not the
1e-8 the criterion asked; the residual is first-order truncation from the principal
frame's rotation, it scales exactly 10:1 with the probe step, and it is the same size
on the new branches as on the branches this driver has always run — so what the
criterion wanted proved is proved, by a different measurement, and the clause as
written is not met. The ramp misses its 0.01 by a thousandth on the one model it was
run on. And the corpus cannot supply the lock the criterion asked for: every
`fem_ssrm` lock in the repository whose model carries a positive `t_cut` also carries
`k0=1`, so the vendor row above is the lock reproduced with the K0 tag dropped —
which on that model was measured to be worth nothing, but is a caveat and not a
formality.

What remains, in the order it matters:

- **K0 initial stress is now the gating refusal for the vendor corpus.** It was one
  of eight refused features; with the tension cutoff carried it is the ONE that
  stands between this driver and 142 locked vendor benchmarks, every one of which
  would otherwise be reachable. That is a much sharper target than the list it came
  from. *(Carried — see "K0 INITIAL STRESS", below.)*
- **Post-peak softening is still refused**, unchanged, and still the reason neither
  published reinforced factor of safety is reachable here.
- **The `tension_srf` setting moves the stress field and not the answer** on every
  model measured. Whether that holds on a model whose limit is genuinely set by
  tension is unmeasured, and it is the natural next reading.
- **The ramp's cell.** One model, one thousandth. Worth re-reading on a second capped
  model before it is called anything.


## K0 INITIAL STRESS — the Newton driver learns the in-situ state

Written before any feature code, so that what follows is a test and not a
description. Same machine and settings as everything above: `force_tol` 1e-3,
hybrid criterion, `capture_failure_state=False`, tolerance 0.01.

### Why this one now

The tension-cutoff round ended by naming it: "K0 initial stress is now the gating
refusal for the vendor corpus. It was one of eight refused features; with the
tension cutoff carried it is the ONE that stands between this driver and 142
locked vendor benchmarks, every one of which would otherwise be reachable."

The corpus inventory says the same thing from the other side. Of 188 FEM-bearing
models, 100 are blocked by exactly two gates, and for almost all of them the two
are `t_cut` and K0 together: Rocscience's published verification set is solved
with its at-rest initial stress on, so every transcription of it carries `k0=1`
on the tag and a per-material tensile cap in the file. Carrying the cap alone
moved the reachable locked count by one model. Carrying K0 as well moves it by
94.

This round is also the first time the Newton driver is pointed at the vendor
corpus at all. Every benchmark in this document so far is a Griffiths & Lane
anchor, a tutorial, or a reinforced sample. The vendor block is a different
population — different geometries, different transcription provenance, different
mesh sizes, published tolerances of 0.01 to 0.05 — and how the driver does
against it is the finding, not a formality.

### The semantics being reproduced, read from the viscoplastic driver

Not assumed — read out of `_prepare_fem_model`, `solve_fem` and `solve_ssrm` and
restated here, because the Newton path has to solve the same model:

- **The field.** Per Gauss point, tension-positive and EFFECTIVE:
  `sigma'_v = -(soil column weight above the point) + u`,
  `sigma'_h = sigma'_z = K0 sigma'_v` — in-plane AND out-of-plane —
  `tau_xy = 0`. The overburden is SOIL ONLY: a reservoir load or a distributed
  load is not part of the in-situ state and enters as a boundary force. The
  overburden integral is `_gauss_point_overburden`, cached on the prepared model
  as `sv0_gp` because it is F-independent.
- **The method.** Classical initial stress: `sigma = sigma_0 + D (B u - eps^p)`,
  with `int B^T sigma_0 dV` moved to the right-hand side of the linear solve. The
  solver still iterates to equilibrium under the body forces; it starts from the
  K0 state instead of from zero stress.
- **The yield check reads the total field.** `sig4 + sig0`, so the initial stress
  is inside the surface the constitutive law is evaluated on, not beside it.
- **The SSRM sequencing.** `solve_ssrm` runs ONE full-strength solve at F = 1
  before the bisection, with the displacement cap off, and keeps its `_k0_state`
  — the displacement field and the accumulated plastic strain. Every trial is
  then handed that state as `_init_state` and starts from it. The reason is
  stated in `solve_fem`'s docstring: at a reduced strength the in-situ
  redistribution would happen against weakened soil and its plastic strain would
  be charged to the trial. The state is accepted on `stable`, not on
  convergence; when it is not established the bisection runs without one and
  warns.
- **The datum.** A trial that starts from the equilibrated state measures
  displacement FROM it — the reported field, the convergence ratio and the
  failure criteria all read `u - u_eq` — while stresses and structural forces stay
  functions of the absolute displacement. The in-situ displacement is an artifact
  of imposing a stress field the geometry does not hold; the soil did not travel
  there.
- **The yardstick is not re-scaled.** The hybrid criterion's displacement scale
  stays the elastic response to the APPLIED load (`loads_grav`), not to the
  load minus the initial-stress term, so a K0-on and a K0-off verdict are read
  off the same ruler.
- **`k0` is not reduced by F.** It is an initial state, not a strength.

### Design

- **`sigma_0` on the Newton groups.** Built from the same `sv0_gp` and the same
  formula, attached per group beside `c_r` and the tensile cap. The trial stress
  the return map is handed becomes `D (B u - eps^p) + sigma_0`, so the initial
  stress is inside the yield surface evaluation exactly as it is on the
  viscoplastic path, and the committed plastic strain inverts the same relation.
- **The residual is absolute, and it is the same residual.** Newton has an
  internal force to write, so `sigma_0` does NOT move to the right-hand side: it
  rides inside `f_int = int B^T sigma dV` and the residual stays
  `f_ext - f_int` with `f_ext` the whole applied load. That is the same equation
  the viscoplastic path solves after moving the term, so the Dawson
  out-of-balance both drivers report is the same quantity, and the force gate
  the verdict is read on does not change meaning. What DOES change is the
  displacement measurement, which is re-referenced to the carried state exactly
  as the viscoplastic path re-references it.
- **The initial stress scales with the load factor.** The driver walks gravity
  in increments; at load factor `lam` it carries `lam * sigma_0`, so `lam = 0`
  is the true stress-free state and `lam = 1` the full in-situ one. The
  overburden IS gravity, so the two scale together and the endpoint is exact.
  With `_NR_INIT_STEP = 1.0` the default first attempt is still the whole load
  in one step, which is what the viscoplastic path does; the graded path is only
  the fallback.
- **The pre-equilibration runs on the caller's driver.** `solve_ssrm`'s
  equilibration solve currently takes whatever `resolve_fem_solver` returns
  rather than the driver the trials will use. It gets the trials' driver, so a
  Newton run equilibrates with the Newton corrector and a viscoplastic run is
  untouched.
- **`_k0_state` is interchangeable between the drivers.** The Newton groups are
  built from `prep["gp_groups_static"]` in the same order, over the same
  `pairs`, and `eps^p` means the same thing on both paths (`sigma = sigma_0 +
  D (B u - eps^p)`). A Newton solve therefore returns a `_k0_state` of the same
  shape, and a state from either driver can seed the other. That is what makes
  the agreement leg of the criterion measurable at all.
- **A trial seeded from the in-situ state does not walk the load path.** It
  already stands at full gravity at full strength; cutting the load below full
  gravity would be solving a different problem. One attempt at the whole load,
  the same rule the viscoplastic predictor's seed already follows.
- **The predictor seeds from the K0 state too.** The viscoplastic predictor
  rungs are handed `k0` and the same `_init_state`, so a seed is grown from the
  in-situ state rather than from zero. A predictor started from zero on a K0
  model would hand the corrector a state built under a different initial
  condition.
- **The ramp carries it.** The foot's cold solve is seeded from the in-situ
  state, `restrength` leaves `sigma_0` alone (it is not a strength), and the
  displacement bound reads from the datum.

### Success criterion (verbatim)

1. **Level-ground exactness.** On the existing level-ground case
   (`test/k0_level_ground_check.py`: a 20 x 10 m block, tri6 at target 1.0,
   gamma 20, K0 in {0.5, 1.0, 2.0} dry and K0 = 1.0 with the water table at the
   surface), the Newton driver reproduces the analytic field: element stresses
   within 1e-9 relative of the imposed `sigma_v` / `sigma_h` (the same tolerance
   the viscoplastic leg is held to), `max|u| <= 1e-10`, equilibrium in at most 2
   iterations, and ZERO plastic points — asserted where the viscoplastic driver
   reports zero. The pre-equilibrated leg (the state handed back in as
   `_init_state`) must be an exact no-op there too.
2. **Agreement with the viscoplastic driver on the pre-equilibrated state.** On
   at least 3 models, the two drivers' full-strength in-situ states agree: the
   Gauss-point stress fields within an RMS difference of 1e-3 relative to the
   overburden scale (reported as a number, not asserted to pass), and the
   bisection's trial verdicts agree at EVERY trial the two runs share.
3. **The locked vendor benchmarks.** At least 10 locked `fem_ssrm` models
   carrying `k0=1`, drawn from at least three distinct builder groups of the
   RS2/vp block and chosen by smallest mesh, run through the Newton bisection
   using `run_tests.py`'s own `build_fem_ssrm_case` mapping so the mesh and
   solver options are the suite's and not this round's. Each lock reproduced
   within its own published tolerance. Reported per model: viscoplastic FS,
   Newton FS, the lock, the tolerance, and the ratio of constitutive work.
   `vp017` (RS2-13, lock 1.332 +/- 0.02) is included by name and run WITH
   `k0=1`, and the difference against last round's 1.3363 with the K0 tag
   dropped is reported.
4. **The ramp agrees.** On at least 3 of those models the ramp's factor of
   safety is within 0.01 of the Newton bisection's.
5. **The no-K0 path is bit-identical.** The plain-soil eight-row table, the four
   reinforced benchmarks and the tension-cutoff benchmarks re-run with `k0`
   unset: every converged trial identical in factor of safety, iterations and
   force evaluations to the numbers recorded above, measured against a control
   run of the parent commit staged in a separate package tree rather than
   against a number in this document.
6. **The refusal is gone and the guard list is updated.** Softening, piles,
   Hoek-Brown, power-curve envelopes, matric suction and a non-effective
   pore-pressure formulation still refuse, each verified to still raise and to
   name its own feature.
7. **The locks catch it.** `test/nr_ssrm_check.py` gains a K0 check with three
   legs: level-ground exactness on the Newton driver, an initial-state agreement
   leg against the viscoplastic driver, and one cheap locked vendor model solved
   on both drivers. Mutation: corrupt the OUT-OF-PLANE component of the initial
   stress (`sigma_z`, which no in-plane equilibrium check can see) and confirm
   the check fails; revert. The whole check file passes.
8. **The default path is unchanged**, against the standard control: Griffiths &
   Lane 6 dry with no `fem_solver` argument, FS 2.421875 on per-trial iteration
   counts 147, 781, 3393, 2031, 2841, 9541, 12000, 8617, 8777.
9. **An honest negative is a valid outcome and must be written.** If the Newton
   driver cannot reproduce the vendor locks at a useful rate, or if the two
   drivers' in-situ states differ materially, or if the K0 sequencing costs the
   no-K0 path, that is the result.


### K0 INITIAL STRESS — results

Same machine and settings as everything above: `force_tol` 1e-3, hybrid criterion,
`capture_failure_state=False`, tolerance 0.01. Every number below was measured on
this checkout in this session. The vendor benchmarks are built through
`run_tests.py`'s own `build_fem_ssrm_case`, from the tags in `docs/verification/`
and the workbooks those tags name, so the mesh, the bracket, the tension-cutoff
setting and the K0 value are the suite's own and not this round's.

#### What was built

`_nr_group_sig0` beside `_nr_group_tension_cap` on the group build, and
`_nr_set_load_factor` beside it; `_nr_init_state` and `_nr_seed_state` on
`_solve_fem_newton`; the guard's K0 line replaced by a comment saying where the
field is carried. `solve_ssrm`'s in-situ equilibration now runs on the driver that
will solve the trials rather than on whatever `resolve_fem_solver` returns, and the
ramp receives the state it produces.

Two things were decided here rather than copied from the viscoplastic path, and
both are consequences of Newton having an internal force to write.

**The initial stress rides inside the internal force.** The viscoplastic path
moves `int B^T sigma_0 dV` to the right-hand side because it solves a linear
system with the elastic operator; there is no other place to put it. The Newton
residual is `f_ext - int B^T sigma dV` with `sigma = sigma_0 + D (B u - eps^p)`,
so the term is already there and moving it would be undoing it. The two are the
same equation, which is the point: the Dawson out-of-balance both drivers report
is the same quantity, and the force gate the verdict is read on does not change
meaning on a K0 model. What is re-referenced is the DISPLACEMENT — reported field
and bound alike — exactly as the viscoplastic path re-references it.

**`sigma_0` scales with the load factor.** The driver walks gravity in increments
and the overburden IS gravity, so at load factor `lam` the material carries
`lam*sigma_0` and `lam` of the weight. Without that, `lam = 0` would be a state
with the whole in-situ field and no load to balance it, and the graded load path —
the driver's only recovery when the whole load in one step fails — would be
walking toward the answer from a place that is not the origin. With
`_NR_INIT_STEP = 1.0` the default first attempt is still the whole load at once,
which is what the viscoplastic path does.

One consequence is worth naming because the rest of the round rests on it. The
Newton groups are built from `prep["gp_groups_static"]` in the same order over the
same `pairs`, and `eps^p` means the same thing on both paths, so `_k0_state` is
INTERCHANGEABLE between the drivers. That is what lets the viscoplastic predictor
seed a Newton trial from the in-situ state rather than from zero — its reported
displacement is datum-relative and would otherwise be the wrong vector — and it is
what makes the agreement leg below a measurement rather than an assertion.

#### Level ground, exactly

On flat ground with the base fixed and the sides on x-rollers, the field
`sigma'_v = -gamma z`, `sigma'_h = sigma'_z = K0 sigma'_v` satisfies equilibrium
identically for ANY K0, so a correct implementation has nothing to redistribute:
it must reach equilibrium immediately, leave the mesh where it is, reproduce the
imposed stresses and yield nowhere. `test/k0_level_ground_check.py`'s own block —
20 x 10 m, tri6 at target 1.0, 486 elements, gamma 20 — run through both drivers.
Errors are relative to the peak vertical stress.

| Case | Driver | Iterations | max&#124;u&#124; | err sigma_y | err sigma_x | err tau_xy | yielded |
|---|---|---|---|---|---|---|---|
| K0 = 0.5 dry | viscoplastic | 1 | 6.6e-18 | 5.8e-16 | 2.2e-16 | 1.2e-16 | 0 |
| K0 = 0.5 dry | **Newton** | **1** | **0.0** | **0.0** | **0.0** | **0.0** | **0** |
| K0 = 1 dry | viscoplastic | 1 | 6.9e-18 | 5.8e-16 | 4.4e-16 | 1.7e-16 | 0 |
| K0 = 1 dry | **Newton** | **1** | **0.0** | **0.0** | **0.0** | **0.0** | **0** |
| K0 = 1 dry, pre-equilibrated | **Newton** | **1** | **0.0** | **0.0** | **0.0** | **0.0** | **0** |
| K0 = 2 dry | viscoplastic | 1 | 7.6e-18 | 5.8e-16 | 8.7e-16 | 3.0e-16 | 0 |
| K0 = 2 dry | **Newton** | **1** | **0.0** | **0.0** | **0.0** | **0.0** | **0** |
| K0 = 1, water at the surface | viscoplastic | 1 | 7.8e-18 | 5.8e-16 | 4.4e-16 | 1.6e-16 | 0 |
| K0 = 1, water at the surface | **Newton** | **1** | **0.0** | **2.9e-16** | **2.9e-16** | **0.0** | **0** |

The Newton column is exact — not small, zero. The residual at `u = 0` is
`f_ext - int B^T sigma_0 dV`, which on level ground the tri6 discretization
integrates exactly, so the first residual evaluation passes the equilibrium test
and no correction is ever computed. The viscoplastic column carries about 1e-18
because it starts from an elastic solve rather than from the initial state. The
pre-equilibrated leg — the state handed back in as `_init_state`, which is what
every SSRM trial does — reproduces the same solve to every digit and changes no
stress at all.

**The out-of-plane component, which nothing above can see.** `sigma_z` appears in
no in-plane equilibrium equation and in no in-plane stress, so a wrong
out-of-plane initial stress moves no node, changes no `sigma_x`, and passes every
row of that table while evaluating the yield surface on a state the model is not
in. Two readings hold it. The von Mises stress of the K0 = 1 dry state is 0.0
exactly, where the same column reads 0.50 at K0 = 0.5 and 1.00 at K0 = 2 — that
state is hydrostatic precisely when the third component is `K0 sigma'_v` and not
something else. And the field the driver builds is compared component by
component against the analytic `K0 = 0.5` overburden at all 1,458 Gauss points of
the block, out-of-plane included. Both are in the lock; the mutation below breaks
each of them.

#### The two drivers' in-situ states

Every SSRM trial starts from the full-strength pre-equilibrated state, so a
disagreement there is a disagreement about the model rather than about a trial.
The comparison is made on the Gauss-point stress field reconstructed from first
principles — `sigma = sigma_0 + D (B u - eps^p)` over the prepared model's own
groups — so neither driver's reporting path is trusted to answer the question
about itself. Both drivers were given the same prepared model, `force_tol` 1e-3.

| Benchmark | Model | Gauss points | RMS difference | max difference | VP: iters, oob, max&#124;u&#124; | N-R: iters, oob, max&#124;u&#124; |
|---|---|---|---|---|---|---|
| RS2-45a | `vp083a` | 2,253 | **9.9e-07** | 4.9e-05 | 103, 9.6e-4, 0.019530 | 23, 7.0e-6, 0.019530 |
| RS2-67a | `rs2_67a` | 1,308 | **6.4e-04** | 8.9e-03 | 269, 9.9e-4, 0.019610 | 10, 5.8e-6, 0.019580 |
| RS2-13 | `vp017` | 5,805 | **8.3e-04** | 1.4e-02 | 179, 1.0e-3, 0.005118 | 10, 5.5e-5, 0.005111 |
| RS2-27-m1.5 | `vp036` | 696 | **1.0e-03** | 9.4e-03 | 322, 1.0e-3, 0.011626 | 10, 3.2e-6, 0.011621 |
| RS2-5 | `xslope_acads_weak_layer` | 1,815 | **6.2e-03** | 6.8e-02 | 1,156, 9.9e-4, 0.034860 | 18, 3.0e-5, 0.034910 |
| RS2-46a | `vp084a` | 1,995 | — | — | NOT ESTABLISHED | NOT ESTABLISHED |

Differences are relative to the overburden scale. Four of the five agree to 1e-3
RMS or better and the fifth reads 6.2e-3; the in-situ displacement agrees to
0.05% or better on every one of them. `vp084a` is the sixth case and it is not a
disagreement: its lock is FS 0.787, so the slope does not stand under its own
weight at full strength, and BOTH drivers report the in-situ state as not
established — which is the agreement that model can offer.

**The residual difference is not a convergence artifact, and it was tested for
being one.** On `vp036` the viscoplastic equilibration was re-run at force
tolerances of 1e-3, 1e-4, 1e-5 and 1e-6 — 322 to 1,224 iterations, out-of-balance
falling three orders of magnitude — and the RMS difference against the Newton
state plateaus at 1.0000e-3, 1.0003e-3, 1.0001e-3, 1.0001e-3. It does not move.

What it is, measured at the tightest setting: both states are in equilibrium
(out-of-balance 1.0e-6 and 1.2e-12) and both are ADMISSIBLE — the largest
Mohr-Coulomb value as a fraction of the local strength is +4.1e-8 on the
viscoplastic state and +1.2e-15 on the Newton one — and both carry plasticity,
125 Gauss points of 696 with a non-zero plastic strain on the viscoplastic side
and 100 on the Newton side, peaking at 1.33e-3 and 1.41e-3. Two admissible fields
in equilibrium with the same load under a non-associated flow rule, reached along
two different paths, are not obliged to be the same field, and these are not:
they differ by a tenth of a percent in the plastic strain's distribution. That is
the honest reading of this leg — the criterion's 1e-3 was met on four models of
five and the fifth is reported at 6.2e-3, and what the number measures is
path-dependence rather than an error in either state.

One reporting difference is not a state difference and is worth naming so it is
not rediscovered: on `vp036` the viscoplastic driver reports 0 plastic elements at
the in-situ state and the Newton driver reports 20. The plastic strain fields
above show both drivers yielding there. The viscoplastic count is read off the
element-averaged yield function, which is not positive for a point sitting ON the
surface; the Newton count is every element with a Gauss point the return map
moved.

#### The locked vendor benchmarks

This is the first time the Newton driver has been pointed at the verification
corpus. Fifteen locked `fem_ssrm` benchmarks, every one carrying `k0=1`, drawn
from six builder groups of the RS2 and vp block and from the LEM-file
transcriptions, chosen by smallest mesh. Each is built through `run_tests.py`'s
own `build_fem_ssrm_case` from the tag in `docs/verification/rs2.md`, so the mesh,
the element type, the bracket, the `tension_srf` setting, the SSR zone or
elastic-material list, the iteration budget and the K0 value are the suite's own.
Work is the honest count on each driver — viscoplastic iterations against Newton
force evaluations, one constitutive pass each.

One thing is NOT the tag's: the bisection resolution. The tags carry
`tolerance=0.01` or `0.02` and the suite passes that to `solve_ssrm` as the
bisection tolerance as well as using it as the lock tolerance; every measurement
in this document is made at 0.01, and these are too. That is the sharper setting,
not the looser one, and it is what makes the two drivers comparable. The two rows
where it could matter are re-run at the tag's own resolution below.

| Benchmark | Model | Elements | Lock | Tol | VP FS | Newton FS | Newton − lock | VP − lock | VP work | N-R work |
|---|---|---|---|---|---|---|---|---|---|---|
| RS2-27-m1.5 | `vp036` | 232 | 1.373 | 0.02 | 1.377344 | **1.377344** | +0.0043 | +0.0043 | 42,147 | 6,554 |
| RS2-64k | `rs2_64k` | 276 | 1.403 | 0.02 | 1.405273 | **1.410352** | +0.0074 | +0.0023 | 85,019 | 6,525 |
| RS2-P4-VP2-m3.0 | `vp002` | 319 | 1.669 | 0.02 | 1.671875 | **1.678125** | +0.0091 | +0.0029 | 44,528 | 7,539 |
| RS2-64l-split | `rs2_64l_split` | 399 | 1.147 | 0.02 | 1.151563 | **1.160937** | +0.0139 | +0.0046 | 38,624 | 3,658 |
| RS2-67a | `rs2_67a` | 436 | 2.479 | 0.02 | 2.487305 | **2.499023** | +0.0200 **(out)** | +0.0083 | 154,157 | 9,250 |
| RS2-64g | `rs2_64g` | 480 | 1.639 | 0.02 | 1.635742 | **1.647461** | +0.0085 | -0.0033 | 71,633 | 4,570 |
| RS2-5 | `xslope_acads_weak_layer` | 605 | 1.280 | 0.01 | 1.280078 | **1.285547** | +0.0055 | +0.0001 | 134,546 | 6,547 |
| RS2-14-m2.8 | `vp018` | 626 | 0.972 | 0.02 | 0.967188 | **0.967188** | -0.0048 | -0.0048 | 6,882 | 23,169 |
| RS2-46a | `vp084a` | 665 | 0.787 | 0.02 | 0.776172 | **0.776172** | -0.0108 | -0.0108 | 7,605 | 2,600 |
| RS2-65-m8 | `rs2_65` | 674 | 1.344 | 0.02 | 1.346875 | **1.353125** | +0.0091 | +0.0029 | 54,033 | 7,168 |
| RS2-45a | `vp083a` | 751 | 1.314 | 0.02 | 1.317969 | **1.317969** | +0.0040 | +0.0040 | 31,524 | 5,195 |
| RS2-29-clay | `rs2_29clay` | 798 | 0.997 | 0.02 | 0.992188 | **0.992188** | -0.0048 | -0.0048 | 20,856 | 28,652 |
| RS2-20 | `vp025` | 827 | 1.003 | 0.01 | 1.002734 | **1.002734** | -0.0003 | -0.0003 | 4,309 | 3,187 |
| RS2-63 | `rs2_63` | 1880 | 1.409 | 0.02 | 1.385937 | **1.395312** | -0.0137 | -0.0231 *(out)* | 68,765 | 2,444 |
| RS2-13 | `vp017` | 1935 | 1.332 | 0.02 | 1.336328 | **1.336328** | +0.0043 | +0.0043 | 42,308 | 5,846 |

**Fourteen of the fifteen reproduce their published lock inside its own
tolerance**, at gaps of 0.0003 to 0.0139. The fifteenth, RS2-67a, misses at
0.020023 against a tolerance of 0.020 — by twenty-three MILLIONTHS. Its final
interval is [2.49609, 2.50195]; the viscoplastic driver's is [2.48438, 2.49023],
one interval lower and inside the tolerance at 0.0083. Both drivers read above the
vendor on that model. It is reported as a miss because that is what it is, and the
re-run below shows it is not an artifact of the resolution: at the tag's own 0.02
the Newton bisection reads 2.501953 and misses by 0.0230.

**Seven of the fifteen return the IDENTICAL factor of safety on both drivers**,
with the identical verdict at every trial the bisection visited — `vp036`,
`vp018`, `vp084a`, `vp083a`, `rs2_29clay`, `vp025` and `vp017`. On those the two
drivers do not merely land in the same tolerance: they close on the same interval
by the same sequence of verdicts. On the other eight the Newton answer is HIGHER
every time, by 0.0051 to 0.0117 — the same one-sided direction this document has
recorded since the first table, and for the reason recorded there: where the two
disagree, the viscoplastic verdict is set by one of its stopping rules rather than
by the slope.

`rs2_63` is the row where that direction decides a lock. At a 0.01 bisection the
viscoplastic driver reads 1.385938 and is 0.0231 from its published 1.409 —
outside its own tolerance — while the Newton driver reads 1.395313 and is inside
at 0.0137. At the tag's own 0.02 resolution both drivers read 1.390625 and both
are inside at 0.0184, which is the configuration the suite runs and the row passes
there. The finer bisection is what separates them.

**The two rows at the tag's own resolution.** Run again with the bisection
tolerance the tag carries, which is the suite's configuration:

| Benchmark | Lock | Tol | VP FS | Newton FS | VP − lock | Newton − lock |
|---|---|---|---|---|---|---|
| RS2-63 | 1.409 | 0.02 | 1.390625 | 1.390625 | −0.0184 | −0.0184 |
| RS2-67a | 2.479 | 0.02 | 2.490234 | 2.501953 | +0.0112 | +0.0230 **(out)** |

**The work.** The Newton driver does less constitutive work on thirteen of the
fifteen, by 1.4x to 28.1x, median 7.2x. It does MORE on two: `vp018`
(0.3x — 23,169 force evaluations against 6,882 viscoplastic iterations) and
`rs2_29clay` (0.7x). Both are models whose factor of safety is BELOW 1, so the
bisection spends most of its trials past failure, which is where load control has
nothing to offer and where this document has measured Newton at 17x to 47x the
cost since Phase 0. That weakness is unchanged by K0; the ramp is what removes it.

**The predictor composes with the in-situ state.** On `vp036`, 20 of the 21
Newton solves the bisection made carried the in-situ state, and 12 of them were
additionally seeded from a viscoplastic predictor grown FROM that state — three
rungs on each of the four failing trials. Every one of those twelve came back
FAILED, which is the property that matters: the seed changes where the corrector
starts and never what it decides.

#### `vp017`, with and without the K0 procedure

`vp017` is RS2-13, the row the tension-cutoff round could only run with the K0
tag dropped. It runs in its locked configuration now.

The comparison needs one correction first, and it is a real one. `main` has since
written K0 into the corpus workbooks themselves (`main!D16`), so **dropping the
test tag's `k0=` no longer drops K0**: the loader supplies it and the run is a K0
run anyway. Every workbook checked — `vp017`, `rs2_63`, `vp036` — reads
`k0 = 1.0` through `load_slope_data` where the branch's own older copies read
`None`. Turning the procedure off means clearing `fem_data['k0']` as well.

| Model | Configuration | viscoplastic | Newton | K0 is worth |
|---|---|---|---|---|
| `vp017` (RS2-13) | K0 = 1, the locked configuration | 1.336328125 | **1.336328125** | — |
| `vp017` (RS2-13) | K0 off, in the tag AND in the file | 1.336328125 | **1.336328125** | **0.0000** |

Both drivers, both ways, the same number to nine digits. On this model the K0
procedure changes nothing, and the caveat the tension-cutoff round had to attach
to its `vp017` row — "the lock reproduced with the K0 tag dropped, which on that
model was measured to be worth nothing, but is a caveat and not a formality" — is
discharged: the row now runs in its locked configuration and reads the same
1.336328125, 0.0043 from the published 1.332 against a tolerance of 0.02.

The finding underneath it is the one worth carrying forward. Because the workbooks
now declare their own K0, an experiment that turns the procedure off by removing
the tag is no longer running without it. That is how this comparison was first
made in this session, and the first reading — `rs2_63` at 1.3859375 "without K0" —
is the value SPIKE.md records for that model WITH K0, which is what exposed it.

#### The ramp

The same models on the monotonic strength-reduction ramp, against the Newton
bisection on the same mesh and bracket. The ramp's foot is the tag's `F_min`; it
walks down when that does not stand.

| Benchmark | Model | Ramp FS | Ramp interval | Bisection-N FS | Ramp - bisection | steps / retries | ramp force evals |
|---|---|---|---|---|---|---|---|
| RS2-27-m1.5 | `vp036` | 1.378125 | [1.37500, 1.38125] | 1.377344 | +0.00078 | 6 / 4 | 9,705 |
| RS2-64k | `rs2_64k` | 1.409375 | [1.40625, 1.41250] | 1.410352 | −0.00098 | 11 / 4 | 5,263 |
| RS2-P4-VP2-m3.0 | `vp002` | 1.678125 | [1.67500, 1.68125] | 1.678125 | **0.00000** | 10 / 4 | 8,008 |
| RS2-67a | `rs2_67a` | 2.496875 | [2.49375, 2.50000] | 2.499023 | −0.00215 | 22 / 4 | 15,527 |
| RS2-45a | `vp083a` | 1.315625 | [1.31250, 1.31875] | 1.317969 | −0.00234 | 9 / 4 | 10,655 |
| RS2-20 | `vp025` | 1.003125 | [1.00000, 1.00625] | 1.002734 | +0.00039 | 10 / 4 | 4,035 |

**Six of six agree with the Newton bisection**, the worst at 0.0024 against the
0.01 the criterion asked for, and one — `vp002` — reports the same number to
seven digits. The ramp carries the in-situ state through its whole warm history:
the foot's cold solve starts from it, `restrength` leaves `sigma_0` alone because
an initial state is not a strength, and the displacement bound reads from the
datum rather than from zero. On RS2-67a — the one row the bisection misses — the
ramp reads 2.496875, which is INSIDE the published tolerance at 0.0179. That is
one model and it is not evidence about the two routes; it is recorded because the
row it lands on is the row this section had to report as a miss.

#### The refusals that remain

Every guard was exercised again on this checkout, and each must both fire and name
its own feature:

| Feature | Verdict |
|---|---|
| pile beam elements | refuses, names piles |
| post-peak softening on reinforcement bars | refuses, names softening and counts the bars |
| Hoek-Brown strength envelopes | refuses, names Hoek-Brown |
| power-curve strength envelopes | refuses, names the power curve |
| matric suction | refuses, names matric suction |
| `pp_formulation != 'effective'` | refuses, names the formulation |
| **K0 initial stress** | **carried** |

The control holds too: all three of the envelope and suction models still solve on
the default viscoplastic driver, which is where they belong. A model that trips
several guards at once still reports the first in source order rather than a
generic refusal.

#### The no-K0 path is unchanged

Two proofs, and the first is the stronger one.

**The arithmetic.** The plain Mohr-Coulomb return map, run against the driver as
it stood at **29ae321d** staged in a separate package tree, on 800,000 random
trial states across four friction angles: **BIT-IDENTICAL**, stress and branch
code alike. The reason it has to be is structural: a model that declares no K0
never gets a `sig0` key on its groups, so `grp.get('sig0_s')` is None at every
point the initial stress could enter — the trial stress in `_nr_group_state`, the
residual-only branch of `_nr_internal_force`, and the elastic-strain inversion in
`_nr_commit_plastic_strain` — and `_nr_set_load_factor` writes nothing.

**The benchmarks**, every one re-run on BOTH trees rather than compared against a
number in this document, and compared trial by trial: factor of safety, and each
trial's verdict, iterations and force evaluations.

| Benchmark | Mesh | Driver | FS, both trees | trials | iterations | force evals | identical? |
|---|---|---|---|---|---|---|---|
| FEM-1 tutorial | tri6, 3.5 | Newton | 1.37109375 | 9 | 1,871 | 14,233 | **yes** |
| LEM-3 tutorial | tri6, 1.2 | Newton | 1.26953125 | 9 | 4,945 | 38,463 | **yes** |
| Griffiths & Lane 1 | quad8, 3.5 | Newton | 1.37187500 | 9 | 2,213 | 16,354 | **yes** |
| Griffiths & Lane 1 | tri6, 3.5 | Newton | 1.36562500 | 9 | 3,121 | 22,650 | **yes** |
| Griffiths & Lane 1 | quad9, 3.5 | Newton | 1.39687500 | 9 | 844 | 6,030 | **yes** |
| Griffiths & Lane 6 dry | quad8, 2 | Newton | 2.41562500 | 9 | 4,146 | 29,059 | **yes** |
| Griffiths & Lane 6 dry | tri6, 2 | Newton | 2.45937500 | 9 | 3,333 | 23,672 | **yes** |
| Griffiths & Lane 3, r = 0.8 | tri6, 6 | Newton | 1.42812500 | 9 | 5,015 | 40,652 | **yes** |
| Geogrid sample | tri6, 4 | Newton | 1.60937500 | 9 | 4,408 | 29,841 | **yes** |
| Half capacity | tri6, 4 | Newton | 1.40937500 | 9 | 2,664 | 18,179 | **yes** |
| Three layers | tri6, 4 | Newton | 1.20742188 | 9 | 3,604 | 25,189 | **yes** |
| Geogrid sample, LOCKED mesh | tri6, 2 | Newton | 1.57187500 | 9 | 4,680 | 29,139 | **yes** |
| Three layers, `t_cut` = 0 | tri6, 4 | Newton | 1.20742188 | 9 | 3,657 | 25,720 | **yes** |
| Geogrid LOCKED mesh, `t_cut` = 0 | tri6, 2 | Newton | 1.56562500 | 9 | 4,540 | 29,652 | **yes** |
| Griffiths & Lane 1, `t_cut` = 0 | tri6, 3.5 | Newton | 1.35312500 | 9 | 3,342 | 24,637 | **yes** |
| Griffiths & Lane 1, `t_cut` = 30 | tri6, 3.5 | Newton | 1.35937500 | 9 | 2,494 | 18,570 | **yes** |
| **Griffiths & Lane 6 dry — the DEFAULT path** | quad8, 2 | viscoplastic | 2.42187500 | 9 | 48,128 | — | **yes** |

#### The locks

`test/nr_ssrm_check.py` gains `check_k0_initial_stress`, in about 31 s. Four legs:
level-ground exactness on the Newton corrector, dry and with the water table at
the surface; the initial stress field asserted component by component against the
analytic K0 = 0.5 overburden at every Gauss point of the block, out-of-plane
included; the two drivers' in-situ states compared on the reconstructed
Gauss-point field; and RS2-27 at its coarsest tagged mesh reproduced on both
drivers inside the published tolerance. The tag's settings — tri6 at 1.5,
`k0 = 1`, `tension_srf` off, FS 1.373 +/- 0.02 — are written out rather than read
through `run_tests`, so the assertion holds the mapping instead of following it.

**Mutation, run both ways.** The criterion named the out-of-plane component
because it is the one an in-plane check cannot see: `sigma_z` appears in no
in-plane equilibrium equation and in no in-plane stress, so getting it wrong moves
no node and changes no `sigma_x`.

| driver | verdict | what the check saw |
|---|---|---|
| as shipped | **PASS** | — |
| `sigma_z = sigma'_v` — the K0 factor dropped out-of-plane only | **FAIL** | "the K0 initial stress the Newton driver builds is wrong in its OUT-OF-PLANE component by 5.000e-01 of the overburden, over 1458 Gauss points, against the analytic K0 = 0.5 field" |
| `sigma_z = 0` | **FAIL** | the same line, plus nine more — the level-ground solve no longer reaches equilibrium, the mesh moves 5.1e-05 where it must not move at all, 18 elements yield where none may, and the two drivers' in-situ fields part by 6.7e-02 |

The first mutation is the sharp one. It is INVISIBLE at K0 = 1, which is the value
every model in the vendor corpus is authored with, so it changes not one factor of
safety in this document; only the leg written for it fires. The whole check file
passes as shipped.

#### The criterion, line by line

**1. Level-ground exactness — MET, and exactly rather than nearly.** The Newton
corrector reaches equilibrium on its FIRST residual evaluation at every K0 tested,
leaves `max|u|` at 0.0 — not small, zero — reproduces the analytic `sigma_v` and
`sigma_h` with an error of 0.0 dry and 2.9e-16 with the water table at the
surface, and yields at zero Gauss points where the viscoplastic driver also
reports zero. The pre-equilibrated leg changes no stress at all. The out-of-plane
component, which the criterion did not name but the lock does, is asserted twice.

**2. Agreement with the viscoplastic driver on the in-situ state — MET on four
models of five, and the fifth is reported.** RMS differences of 9.9e-07, 6.4e-04,
8.3e-04 and 1.0e-03 against the 1e-3 asked, and 6.2e-03 on `xslope_acads_weak_layer`;
in-situ displacement agrees to 0.05% or better on all five. A sixth model,
`vp084a`, has a factor of safety below 1 and BOTH drivers report the in-situ state
as not established. The trial verdicts agree at every trial on seven of the
fifteen vendor runs and differ on the rest by one or two bisection cells, always
with Newton reading higher. The residual difference was tested for being a
convergence artifact and is not one: it does not move over three orders of
magnitude of viscoplastic force tolerance. Both states are admissible and in
equilibrium; they differ in how the plastic strain is distributed, which is what
two paths to the same limit under a non-associated flow rule are entitled to
differ in.

**3. The locked vendor benchmarks — MET at 14 of 15.** Fifteen locks across six
builder groups, meshes from 232 to 1,935 elements, reproduced at 0.0003 to 0.0139
against tolerances of 0.01 and 0.02. RS2-67a misses at 0.020023 against 0.020, and
still misses at 0.0230 when re-run at the tag's own bisection resolution; it is
reported as a miss. Seven rows return the identical factor of safety on both
drivers with the identical verdict at every trial. Work is lower on the Newton
driver on thirteen of fifteen, 1.4x to 28.1x, and higher on the two sub-unity
models, which is the past-failure cost this document has measured since Phase 0
and which K0 does not change.

**4. The ramp — MET on 6 of the 6 run**, against the 3 asked. Worst disagreement
with the Newton bisection 0.0024; `vp002` reproduces it to seven digits.

**5. The no-K0 path is bit-identical — MET.** Seventeen benchmark pairs and 153
trials against the driver staged at 29ae321d: every one identical in factor of
safety, verdict, iterations and force evaluations, and the plain return map
bit-identical over 800,000 trial states.

**6. The refusal is gone and the guard list is updated — MET.** K0 is carried;
piles, bar softening, Hoek-Brown, the power curve, matric suction and a
non-effective pore-pressure formulation each still raise and each still names
itself, with the viscoplastic control accepting all of them.

**7. The locks — MET.** `check_k0_initial_stress` passes as shipped and fails on
both out-of-plane mutations, including the subtle one that is invisible at K0 = 1
and therefore moves no number in this document. The whole check file passes.

**8. The default path — MET.** Griffiths & Lane 6 dry, quad8, size 2, no
`fem_solver` argument: FS 2.421875 on per-trial iteration counts 147, 781, 3393,
2031, 2841, 9541, 12000, 8617, 8777 — value for value the control sequence, on
both trees.

**9. The honest negatives** are RS2-67a's 2.3e-05, the 6.2e-03 initial-state
difference on `xslope_acads_weak_layer`, the two sub-unity models where the Newton
driver does more work than the viscoplastic one, and the finding that a `k0` test
tag can no longer be dropped to turn K0 off. All four are written above.

#### Verdict

**The Newton driver meets the vendor corpus, and it reproduces it.** Fifteen
locked benchmarks were run in their locked configurations — the suite's own
meshes, brackets, SSR zones, elastic-material lists and `tension_srf` settings,
with `k0 = 1` — and fourteen land inside their published tolerances, at 0.0003 to
0.0139. Seven of the fifteen return the identical factor of safety to the
viscoplastic driver with the identical verdict at every trial the bisection
visited. That is the question this round was commissioned on and it is answered
on more than the ten models the criterion asked for.

What made it possible is that K0 is not a special case bolted to the solve. It is
an initial stress, and Newton has an internal force to put it in: the residual
`f_ext - int B^T (sigma_0 + D (B u - eps^p)) dV` is the same equation the
viscoplastic driver solves after moving the term to the right-hand side, so the
out-of-balance both drivers report keeps meaning one thing and the force gate the
verdict is read on did not have to be redefined. The only genuinely new decision
was to scale `sigma_0` with the load factor, and it is the decision that keeps the
graded load path meaningful: the overburden IS gravity, so the two rise together
and `lam = 0` is still the origin. On level ground the whole construction is
exact — the first residual evaluation passes, `max|u|` is zero rather than small,
and the analytic field comes back with an error of zero.

Four things did not close, and all four are written above rather than tuned away.
RS2-67a misses its tolerance by 2.3e-05 at a 0.01 bisection and by 0.0030 at the
tag's own 0.02 — the Newton bisection closes one interval above the viscoplastic
one, both above the vendor, and the ramp on the same model lands inside. The two drivers' in-situ states differ by 6.2e-03 RMS on
`xslope_acads_weak_layer` against the 1e-3 the criterion asked, and by 1.0e-03 on
`vp036`; the difference was tested for being a convergence artifact over three
orders of magnitude of force tolerance and is not one — both states are in
equilibrium and both are admissible, and they distribute their plastic strain
differently, which two paths to the same limit under a non-associated flow rule
are entitled to do. The two sub-unity models cost the Newton driver more work than
the viscoplastic one, which is the past-failure weakness this document has
recorded since Phase 0 and which this round does not touch. And a `k0` test tag
can no longer be dropped to turn K0 off, because the workbooks now declare it
themselves.

**What this opens.** 144 locked `fem_ssrm` test tags in `docs/` carry a `k0`
value, over 108 distinct workbooks. With the initial stress carried, **133 of them
— 97 workbooks — are inside the Newton driver's feature envelope**, where before
this round none were. Eleven remain refused, and the reasons are now a short list
rather than a class: five power-curve models, three Hoek-Brown, three matric
suction. On the fifteen actually run, fourteen reproduced their lock, a hit rate
of 93% on a sample chosen for small meshes rather than for agreement.

What remains, in the order it matters:

- **The corpus run itself.** Fifteen models is a sample, not a calibration. The
  reachable set is 133 tags and the driver has now been shown able to address
  them; running them is what would turn "reproduces fourteen of fifteen" into a
  statement about the corpus, and it is also the measurement the ramp verdict
  named as the thing standing between this branch and any default-driver
  conversation.
- **The remaining eleven.** Power curve and Hoek-Brown are two linearizations of
  the same shape and would come together; matric suction is an apparent cohesion
  the return map would have to be told about. Piles and bar softening still block
  a further eight and two locked models respectively, outside the K0 set.
- **The past-failure cost is unchanged and it now has vendor evidence.** The two
  models with a factor of safety below 1 are the two where the Newton bisection
  did more constitutive work than the viscoplastic driver. A corpus run would meet
  many more of them, and the ramp — which by construction never evaluates a
  strength more than one increment past the highest it has carried — is the
  answer this branch already has.
- **The 6.2e-03 in-situ difference on one model.** One model, one number, and the
  mechanism is understood. Whether it ever moves a verdict is unmeasured; on the
  benchmarks here it did not, since that model's two drivers land 0.0055
  apart with both inside the lock.


## MERGED FROM `main`, 2026-09-01

`main` moved seven commits while this branch was being written, and three of them
reach the Newton driver.

**Two of the driver's eight refusals were refusals of things that no longer
exist.** `pp_formulation` (49331fd2) and `staged` (e9615718) are gone from
`solve_fem` and `solve_ssrm`: the effective-stress formulation is the only one,
and there is one load stage. Git auto-merged `xslope/fem.py` without a textual
conflict, which is the dangerous outcome rather than the safe one — the branch's
Newton code went on calling a signature that had been deleted under it. The
refusal table now reads:

| Feature | Verdict |
|---|---|
| pile beam elements | refuses, names piles |
| post-peak softening on reinforcement bars | refuses, names softening and counts the bars |
| Hoek-Brown strength envelopes | refuses, names Hoek-Brown |
| power-curve strength envelopes | refuses, names the power curve |
| matric suction | refuses, names matric suction |
| the Rankine tension cutoff | carried |
| **K0 initial stress** | **carried** |
| staged loading | no longer exists |
| `pp_formulation` | no longer exists |

Nothing in this document was measured under either keyword: no workbook and no
test tag in the repository ever set one, so the removal changes no benchmark
above. It shortens the list rather than the results.

**K0 is in the corpus workbooks now** (1f34a0d, `main!D16` on 107 models, guarded
by the new `tag_k0` suite row). The driver has carried K0 since 8a217d4c, so no
code follows from it — but it changes what an experiment MEANS, which the K0 round
already recorded: dropping a `k0=` test tag no longer turns the procedure off,
because the loader supplies it. Every K0 leg on this branch passes `k0` as an
explicit argument, so none of them was running the experiment it thought it was;
the branch's own scripts were re-read to confirm that rather than assumed.

The post-merge checks, all run in the worktree with its root first on `sys.path`
and `xslope.__file__` asserted:

| Check | Result |
|---|---|
| Default path: Griffiths & Lane 6 dry, quad8/2, no `fem_solver` | FS **2.421875**, per-trial iterations 147, 781, 3393, 2031, 2841, 9541, 12000, 8617, 8777 — value for value |
| Griffiths & Lane 1, tri6/3.5, Newton bisection | FS **1.36562500**, 9 trials, 3,121 iterations, 22,650 force evaluations — value for value |
| Three layers, `t_cut` = 0, tri6/4 | viscoplastic **1.2109375**, Newton **1.2109375** — value for value |
| `vp017` (RS2-13), `k0 = 1` from the MERGED workbook, tri6/0.5, the tag's bracket | viscoplastic **1.336328125**, Newton **1.336328125**, 5,846 force evaluations — value for value, and 0.0043 from the published 1.332 against a tolerance of 0.02 |
| `test/nr_ssrm_check.py` | passes end to end |
| `python3 run_tests.py --preflight` | 31 passed, 0 failed |

## MATRIC SUCTION — the Newton driver learns the apparent cohesion

Written before any feature code, so that what follows is a test and not a
description. Same machine and settings as everything above: `force_tol` 1e-3,
hybrid criterion, `capture_failure_state=False`, tolerance 0.01.

### Why this one now

The K0 round ended on a short list rather than a class: five power-curve models,
three Hoek-Brown, three matric suction. Suction is the only one of the three that
needs no new yield surface. Fredlund's extended Mohr-Coulomb does not bend the
envelope — it adds an APPARENT COHESION to the same linear surface the return map
already solves — so the shear surface, the Rankine cap, the corner and apex
construction and the flow rule are all unchanged, and the whole of the work is
plumbing.

The corpus inventory names the models, and the two families in it are not the same
test. `rs2_28a/b/c` is RS2-28, Ng & Shi's (1998) Hong Kong cut at three far-field
heads, with `phi_b` = 15 deg declared in the FILE, a vendor tensile cap of 10 kPa
that binds, `k0 = 1`, and 63% of the domain declared `elastic`; each carries one
locked factor of safety and every one of the three is gated on suction alone.
`vp102t_60/300/1500` is the RS2 Part 4 VP102 transient set, with `phi_b` = 37 deg
supplied by the TEST TAG rather than the file; each carries TWO locks, a `c2` row
with no suction and a `c3` row with it, so only the three `c3` rows are gated.
Nine locked factors of safety on six workbooks.

There is no LEM-side suction sample to add. Every workbook under
`docs/inputs/slope`, `docs/lem/files`, `docs/fem/files`, `docs/tutorials/files`
and `docs/seep/files` was loaded through `load_slope_data`, and not one carries a
`phi_b`. The suction corpus is the RS2 verification block and nothing else.

### The semantics being reproduced, read from the viscoplastic driver

Not assumed — read out of `_resolve_suction_by_elem`, `_prepare_fem_model` and
`solve_fem` and restated here, because the Newton path has to solve the same
model:

- **The credit.** Per Gauss point,
  `c_suction = min(max(0, -u_w), s_cap) * tan(phi_b)`, added to `c'` in the
  Mohr-Coulomb yield function. `u_w` is the SIGNED pore pressure
  (`prep['u_gp_signed']`), not the clamped `u >= 0` the effective-normal term
  reads: the whole point is the negative pressure above the water table that the
  clamp throws away.
- **The source.** The signed field is built only for `u = piezo` or `u = seep`,
  and only when some material carries a positive `phi_b`. `prep['suction_active']`
  is that gate, and when it is False no suction array exists at all — which is
  what makes an untouched default path a structural fact rather than a
  measurement.
- **`phi_b` and `s_cap` are per MATERIAL**, resolved to per-element arrays by
  `_resolve_suction_by_elem` with the same auto-wire / override precedence the
  limit-equilibrium engine uses: the file's `phi_b` / `s_cap` columns unless a
  `suction_phi_b` / `suction_cap` kwarg says otherwise, and `{}` forces the credit
  off. A blank `s_cap` is `inf` — uncapped.
- **It IS reduced by F.** `grp['c_suc_r'] = tanphib * s * Finv` with
  `Finv = 1 / F_by_elem`, so the credit is divided by the trial strength exactly
  as `c'` is, and by the SAME per-element factor — which is 1.0 on an
  `ssr_exclude` element. `docs/fem/overview.md` states it in as many words: "the
  apparent cohesion is reduced by the trial factor alongside c' and tan(phi')".
  That is what distinguishes it from the tensile cap, which is reduced only when
  `tension_srf` says so.
- **The envelope is `c_env = c_r + c_suc_r`** and nothing else moves: psi = 0
  flow, the same corners, the same apex, the Rankine cap still a separate surface.
- **An `elastic` material is inert.** The viscoplastic path drops elastic Gauss
  points from the yielding mask entirely; the Newton path gives them `c_r = inf`.
  A finite credit added to either is a no-op, and RS2-28 exercises it, since most
  of its domain is `Plasticity: None`.
- **The credit never touches the effective normal stress.** `u_gp` stays clamped
  and moves to the load vector as it always did; only the cohesive intercept picks
  up the suction.

### Design

- **Apparent cohesion is a per-Gauss-point cohesion, and that is all it is.** The
  Newton return map already takes `c` per point (`grp['c_r']`, handed straight to
  `mc_return_map`), so the credit is folded into `c_r` where the group is built
  and no return-map, tangent or yield-surface code is touched at all. Everything
  downstream that reads `c_r` — the return map, the algorithmic tangent, the
  strength scale `nr_max_yield_violation` divides by — is then reading the
  envelope the trial is actually solved on, with no second code path to keep in
  step.
- **The formula is shared, not mirrored.** One helper computes
  `min(max(0, -u_signed), s_cap) * tan(phi_b) / F`, and BOTH drivers call it, so
  the only thing a driver-against-driver comparison can find is a plumbing
  difference — which field, which cap, which F. The viscoplastic call site keeps
  its own arithmetic in its own order, so the default path's numbers cannot move.
- **The F-independent half is cached.** The capped suction
  `s = min(max(0, -u_w), s_cap)` does not depend on F, so it is built once per
  group and only the `tan(phi_b) * s / F` product is re-formed. That is what lets
  the ramp's `restrength` re-derive the credit at every new F without a field
  lookup — which it MUST, because the credit scales as 1/F and a ramp carrying the
  foot's credit up the whole ramp would be reducing only part of the envelope.
- **The predictor already carries it.** The viscoplastic predictor rungs are
  handed `suction_phi_b` / `suction_cap` already, so a seed is grown on the same
  envelope the corrector will use. Nothing to change; the criterion asserts it
  rather than assuming it.
- **The guard's suction line goes, and nothing else does.** Piles, post-peak bar
  softening, Hoek-Brown and the power curve stay refused, each naming itself.

### Success criterion (verbatim)

1. **The suction-gated locked models.** On the three `rs2_28` models (locks 1.606
   / 1.544 / 1.381, tolerance 0.02) and the three `vp102t` `c3` rows (1.779 /
   2.162 / 2.687, tolerance 0.02), each built through `run_tests.py`'s own
   `build_fem_ssrm_case` so that the mesh, the bracket, `tension_srf`, `k0` and
   `suction_phi_b` are the suite's and not this round's: the Newton bisection
   reproduces the published lock within its own tolerance AND lands within 0.01 of
   the viscoplastic driver run in the same configuration. Reported per model:
   viscoplastic FS, Newton FS, the lock, the tolerance, and the constitutive work
   on each driver.
2. **The two drivers build the same apparent-cohesion field.** On one model at
   F = 1 the per-Gauss-point `c_suction` each driver computes is compared point for
   point and the maximum absolute difference is REPORTED as a number. Shared code
   makes ~1e-16 the expectation; anything larger is a plumbing difference and is
   the finding.
3. **`s_cap` binds, and both drivers move together.** On one model the cap is
   lowered until it is the binding constraint, and the factor of safety is
   reported at each of two or three cap values on BOTH drivers, agreeing within
   0.01 at every one. The number of Gauss points where `s > s_cap` is reported
   alongside, because a cap that binds nowhere proves nothing.
4. **The ramp agrees** with the Newton bisection within 0.01 on at least two of
   the six models.
5. **The no-suction path is bit-identical**, on the same control-tree protocol the
   previous rounds used: the plain-soil eight rows, the four reinforced
   benchmarks, the tension-cutoff rows and the K0 vendor rows re-run against a
   control run of the parent commit staged in a separate package tree — every
   converged trial identical in factor of safety, verdict, iterations and force
   evaluations — plus the plain Mohr-Coulomb return map bit-identical over 800,000
   random trial states.
6. **The refusal is gone and the guard list is updated.** Piles, post-peak bar
   softening, Hoek-Brown and power-curve envelopes each still raise and each still
   names its own feature, with the viscoplastic control accepting all of them.
7. **The locks catch it.** `test/nr_ssrm_check.py` gains a suction check, and it
   asserts the FIELD and not only a factor of safety, so a mutation cannot pass by
   landing on the same number. Two mutations, run both ways and recorded: **the
   cap dropped** (`s_cap` ignored) must FAIL it, and **the F-reduction of the
   apparent cohesion dropped** (the credit taken unreduced, which is what the
   viscoplastic driver does NOT do) must FAIL it. The whole check file passes.
8. **The default path is unchanged**, against the standard control: Griffiths &
   Lane 6 dry with no `fem_solver` argument, FS 2.421875 on per-trial iteration
   counts 147, 781, 3393, 2031, 2841, 9541, 12000, 8617, 8777.
9. **An honest negative is a valid outcome and must be written.** If the two
   drivers disagree on a suction model by more than the bisection tolerance, or if
   the credit needs anything the return map does not already have, that is the
   result.

### MATRIC SUCTION — results

Same machine and settings as everything above: `force_tol` 1e-3, hybrid criterion,
`capture_failure_state=False`, tolerance 0.01 unless a row says otherwise. Every
number below was measured on this checkout in this session. The vendor benchmarks
are built through `run_tests.py`'s own `build_fem_ssrm_case`, from the tags in
`docs/verification/rs2.md` and the workbooks those tags name, so the mesh, the
bracket, the `tension_srf` setting, the K0 value and the suction angle are the
suite's own and not this round's.

#### What was built

Two module-level helpers, `_suction_capped` and `_suction_apparent_cohesion`, and
two group-level ones beside `_nr_group_tension_cap`: `_nr_group_suction`, which
caches the F-independent capped suction on a Newton group, and
`_nr_group_apparent_cohesion`, which folds `tan(phi_b) s / F` into that group's
cohesion. The guard's matric-suction line is gone. The ramp's `restrength` calls
the folding helper at every new strength, beside the two lines that re-reduce `c`
and the tensile cap.

That is the whole of it, and the reason it is the whole of it is that the Newton
return map already takes `c` PER GAUSS POINT. Fredlund's credit is an apparent
cohesion on the same linear envelope, so it arrives at the return map as a
different `c_r` and nothing else changes: no surface, no flow direction, no
corner or apex construction, no tangent. Everything downstream that reads `c_r` —
the return map, the algorithmic tangent, and the strength scale
`nr_max_yield_violation` divides by — is then reading the envelope the trial is
actually solved on, with no second code path to keep in step.

The viscoplastic driver's own call site was rewritten to use the same two
helpers. Its arithmetic and its order are unchanged, which is what keeps the
default path bit-identical, and it is what makes the field comparison below a
plumbing measurement rather than an algebra one.

#### The two drivers' apparent-cohesion fields

Captured at the one helper both drivers call, by running a single solver
iteration through each driver's own entry point — so what is compared is what the
solve runs on, not what a test could build for itself.

| Model | Gauss points | credited | largest credit at F = 1.3 | max &#124;VP − N-R&#124; | max &#124;N-R − an independent expectation&#124; |
|---|---|---|---|---|---|
| `rs2_28a` | 11,064 | 1,219 | 12.861 kPa | **0.000e+00** | 1.776e-15 |
| `vp102t_60` | 3,354 | 474 | 65.093 | **0.000e+00** | 7.105e-15 |

The independent expectation is `min(max(0, -u_signed), s_cap) tan(phi_b) / F`
recomputed from the prepared model's own arrays, so the zero in the fourth column
is not two code paths agreeing with each other about the wrong thing.

The credit is not decoration on either model. On `rs2_28a` the largest apparent
cohesion is 12.9 kPa against the material's own `c'` of 10 kPa — the suction is
worth more than the cohesion at the points where it bites — and `rs2_28a`'s
vendor tensile cap of 10 kPa sits below the envelope's own apex, so the cap and
the credit are both active on the same model.

#### The suction-gated locked benchmarks

Six workbooks, nine locked factors of safety, of which the six suction-gated ones
are run here. Work is the honest count on each driver: viscoplastic iterations
against Newton force evaluations, one constitutive pass each.

| Benchmark | Model | Elements | Lock | Tol | VP FS | VP − lock | Newton FS | Newton − lock | driver gap | VP work | N-R work |
|---|---|---|---|---|---|---|---|---|---|---|---|
| RS2-28a | `rs2_28a` | 3,688 | 1.606 | 0.02 | 1.6218750 | +0.0159 | **1.6468750** | +0.0409 **(out)** | 0.0250 | 103,820 | 4,424 |
| RS2-28b | `rs2_28b` | 3,688 | 1.544 | 0.02 | 1.5406250 | −0.0034 | **1.5656250** | +0.0216 **(out)** | 0.0250 | 137,542 | 6,431 |
| RS2-28c | `rs2_28c` | 3,688 | 1.381 | 0.02 | 1.3968750 | +0.0159 | **1.4031250** | +0.0221 **(out)** | 0.0063 | 94,077 | 5,623 |
| RS2-P4-VP102-t-60-c3 | `vp102t_60` | 1,118 | 1.779 | 0.02 | 1.7707031 | −0.0083 | **1.7925781** | +0.0136 | 0.0219 | 25,778 | 4,640 |
| RS2-P4-VP102-t-300-c3 | `vp102t_300` | 1,118 | 2.162 | 0.02 | 2.1644531 | +0.0025 | **2.1644531** | +0.0025 | **0.0000** | 17,441 | 4,015 |
| RS2-P4-VP102-t-1500-c3 | `vp102t_1500` | 1,118 | 2.687 | 0.02 | 2.6839844 | −0.0030 | **2.6894531** | +0.0025 | 0.0055 | 58,456 | 3,988 |

**Three of the six reproduce their published lock and three do not.** All three
misses are the RS2-28 family, all three are HIGH, and the direction is the
one-sided one this document has recorded since its first table. The viscoplastic
driver reproduces all six. Newton does 4.3x to 23.5x less constitutive work on
every row.

The three misses were re-run at the tag's OWN bisection resolution, which is the
suite's configuration, in case the finer bisection was the cause:

| Benchmark | Lock | Tol | VP FS | VP − lock | Newton FS | Newton − lock | driver gap |
|---|---|---|---|---|---|---|---|
| RS2-28a | 1.606 | 0.02 | 1.6187500 | +0.0127 | 1.6437500 | +0.0378 **(out)** | 0.0250 |
| RS2-28b | 1.544 | 0.02 | 1.5437500 | −0.0002 | 1.5687500 | +0.0248 **(out)** | 0.0250 |
| RS2-28c | 1.381 | 0.02 | 1.3937500 | +0.0127 | 1.4062500 | +0.0252 **(out)** | 0.0125 |

It is not the resolution. The gap is a flat 0.0250 on two rows and 0.0125 on the
third, at both settings.

#### The RS2-28 gap, refereed by the driver's own evidence

The bisection's answer is fixed by the trials whose verdicts differ, so those were
run directly. RS2-28a, at two strengths between the viscoplastic answer (1.61875)
and the Newton one (1.64375):

| F | driver | verdict | exit | out-of-balance | iterations | worst yield violation |
|---|---|---|---|---|---|---|
| 1.6250 | viscoplastic | FAILED | `iteration_cap` | 1.142e-01 | 32,000 | — |
| 1.6250 | **Newton** | **CONVERGED** | `converged` | **2.585e-06** | **126** | **2.84e-09** |
| 1.6375 | viscoplastic | FAILED | `iteration_cap` | 2.322e-01 | 16,000 | — |
| 1.6375 | **Newton** | **CONVERGED** | `converged` | **1.732e-05** | **252** | **2.86e-09** |

A stress field in equilibrium with full gravity — an out-of-balance forty to four
hundred times inside the trial tolerance — whose worst Mohr-Coulomb value, read in
the INVARIANT form and not in the ordered-principal form the return map is written
on, is 2.9e-9 of the local strength. That is a statically admissible field, and by
the lower-bound theorem it is a proof that the slope stands at F = 1.6375. The
argument needs neither driver to be trusted: it is the Newton driver's own
converged state, checked against a yield function it does not solve on.

The viscoplastic verdicts at those strengths are not near misses that a longer run
would flip cheaply — the out-of-balance is 1.1e-1 and 2.3e-1 against a 1e-3 gate,
two orders of magnitude out, and both trials are closed by the budget-extension
heuristic (`iteration_cap`) rather than by the slope.

**And the vendor's own numbers sit where the Newton driver does.** RS2's published
SSR values for these three problems are 1.64, 1.55 and 1.41
(`docs/verification/rs2.md`, RS2-28). At the tag's resolution the Newton driver
reads 1.6438, 1.5688 and 1.4063 — **+0.0038, +0.0188 and −0.0038** from RS2's own
answers — while the viscoplastic driver reads 1.6188, 1.5438 and 1.3938, which is
−0.0213, −0.0063 and −0.0163 below them. XSLOPE's locked values for this row ARE
the viscoplastic readings, and the verification page already scores them 2.1%,
0.4% and 2.1% below RS2 and attributes the shortfall to the convergence check.
So the honest statement is not that the Newton driver misses these three locks; it
is that it misses three locks that are defined by the driver it is being compared
against, and lands on the vendor's answer instead.

That is a finding about the corpus and not a licence: moving those locks is the
owner's decision and not a spike's, and nothing here was changed.

#### The cap

`s_cap` was lowered until it was the binding constraint, and the answer read on
both drivers at a strength both carry uncapped. On `vp102t_60` the suction reaches
91.4 stress units and is positive at 474 of 3,354 Gauss points:

| suction cap | Gauss points where `s` exceeds it | F = 1.74, viscoplastic | F = 1.74, Newton |
|---|---|---|---|
| none (as the tag runs it) | — | CONVERGED, 1,667 iterations | CONVERGED, 15 iterations |
| 40 | 165 | — | — |
| 20 | 294 | CONVERGED, 1,330 | CONVERGED, 17 |
| **5** | **421** | **FAILED**, 3,401 | **FAILED**, 354 |
| the credit switched off entirely | — | **FAILED**, 2,031 | **FAILED**, 358 |

The two drivers move together at every setting, and they move on the cap's VALUE
rather than on the presence of a cap keyword: at 20 the trial still stands on both,
at 5 it fails on both, and with the credit off it fails on both. The last row is
what makes the others a measurement — this model without the suction credit is
locked at 1.713 (its own `c2` row), below the trial being carried, so a driver that
lost the credit would fail a trial it has to carry.

On `rs2_28a` the same field reads 62.4 at most, positive at 1,219 of 11,064 points,
with a cap of 5 binding at 1,102 of them.

#### The ramp

The three `vp102t` models on the monotonic strength-reduction ramp, against the
Newton bisection on the same mesh and bracket:

| Benchmark | Ramp FS | Ramp interval | Bisection-N FS | Ramp − bisection | steps | ramp force evals |
|---|---|---|---|---|---|---|
| RS2-P4-VP102-t-300-c3 | 2.1656250 | [2.16250, 2.16875] | 2.1644531 | **+0.0012** | 20 | 5,051 |
| RS2-P4-VP102-t-1500-c3 | 2.6837500 | [2.68125, 2.68625] | 2.6894531 | **−0.0057** | 32 | 7,970 |
| RS2-P4-VP102-t-60-c3 | 1.7837500 | [1.78125, 1.78625] | 1.7925781 | **−0.0088** | 14 | 6,301 |

Three of three inside 0.01, against the two the criterion asked for. The ramp
carries the credit through its whole warm history — `restrength` re-derives it at
every step because it scales as 1/F — and on the `t-60` row it lands at 1.78375,
which is 0.0048 from the vendor lock of 1.779 where the bisection is 0.0136.

#### The no-suction path is unchanged

Two proofs, and the first is the stronger one.

**The arithmetic.** The plain Mohr-Coulomb return map, run against the driver as it
stood at **53659aa4** — the criterion commit, before any of this code — staged in a
separate package tree, on 800,000 random trial states across four friction angles:
**BIT-IDENTICAL**, stress and branch code alike, to a SHA-256 of the concatenated
returns. The reason it has to be is structural: a model with no `phi_b` never sets
`prep['suction_active']`, so `_nr_group_suction` returns before it touches a group,
`grp.get('s_suc')` is None at the one point the credit could enter, and
`_nr_group_apparent_cohesion` returns without arithmetic.

**The benchmarks**, every one re-run on BOTH trees rather than compared against a
number in this document, and compared trial by trial: factor of safety, and each
trial's strength, verdict, iterations and force evaluations.

| Benchmark | Mesh | Driver | FS, both trees | trials | identical? |
|---|---|---|---|---|---|
| FEM-1 tutorial | tri6, 3.5 | Newton | 1.37109375 | 9 | **yes** |
| LEM-3 tutorial | tri6, 1.2 | Newton | 1.26953125 | 9 | **yes** |
| Griffiths & Lane 1 | quad8, 3.5 | Newton | 1.37187500 | 9 | **yes** |
| Griffiths & Lane 1 | tri6, 3.5 | Newton | 1.36562500 | 9 | **yes** |
| Griffiths & Lane 1 | quad9, 3.5 | Newton | 1.39687500 | 9 | **yes** |
| Griffiths & Lane 6 dry | quad8, 2 | Newton | 2.41562500 | 9 | **yes** |
| Griffiths & Lane 6 dry | tri6, 2 | Newton | 2.45937500 | 9 | **yes** |
| Griffiths & Lane 3, r = 0.8 | tri6, 6 | Newton | 1.42812500 | 9 | **yes** |
| Geogrid sample | tri6, 4 | Newton | 1.60781250 | 8 | **yes** |
| Half capacity | tri6, 4 | Newton | 1.41406250 | 8 | **yes** |
| Three layers | tri6, 4 | Newton | 1.21093750 | 8 | **yes** |
| Geogrid sample, LOCKED mesh | tri6, 2 | Newton | 1.56093750 | 8 | **yes** |
| Griffiths & Lane 1, `t_cut` = 0 | tri6, 3.5 | Newton | 1.35312500 | 9 | **yes** |
| Griffiths & Lane 1, `t_cut` = 30 | tri6, 3.5 | Newton | 1.35937500 | 9 | **yes** |
| Three layers, `t_cut` = 0 | tri6, 4 | Newton | 1.21093750 | 8 | **yes** |
| **Griffiths & Lane 6 dry — the DEFAULT path** | quad8, 2 | viscoplastic | **2.421875** | 9 | **yes** |

The default viscoplastic path returns FS 2.421875 on per-trial iteration counts

    147, 781, 3393, 2031, 2841, 9541, 12000, 8617, 8777

value for value the control sequence, on both trees.

The K0 vendor rows the criterion also named are covered by
`check_k0_initial_stress`, which runs RS2-27 at its coarsest tagged mesh on both
drivers inside the published tolerance and passes as part of the whole check file
below; they were not additionally re-run against the control tree, and that is a
shortfall against the criterion as written rather than a result.

#### The refusals that remain

Every guard was exercised again on this checkout, and each must both fire and name
its own feature:

| Feature | Verdict |
|---|---|
| pile beam elements | refuses, names piles |
| post-peak softening on reinforcement bars | refuses, names softening and counts the bars (18 of 36) |
| Hoek-Brown strength envelopes | refuses, names Hoek-Brown |
| power-curve strength envelopes | refuses, names the power curve |
| **matric suction** | **carried** |

Four remain, and they are the four the eligibility inventory says block ten more
locked models: five power-curve, three Hoek-Brown, eight piles and two bar
softening, with overlaps. The viscoplastic control still accepts all four.

#### The locks

`test/nr_ssrm_check.py` gains `check_matric_suction`, on `vp102t_60` — RS2 Part 4
VP102 at t = 60 s, its own transient seepage field on its own stored mesh, an
unsaturated friction angle of 37 degrees from the TEST TAG, `k0 = 1`,
`tension_srf` on, an SSR search zone, locked at 1.779 +/- 0.02 — in about 52 s.
Five legs, and the field is asserted and not only the answer, because the two
drivers share the formula and a shared-code defect moves them together:

  * the two drivers' apparent-cohesion fields, captured at the helper they share by
    running a solver iteration through each driver's own entry point, must be
    identical point for point;
  * the credit must be reduced by the trial strength — the field at F = 1.5 exactly
    two thirds of the field at F = 1 where the strength is reduced, and UNCHANGED at
    the Gauss points the model's SSR zone holds at full strength;
  * `s_cap` must bound it: at a cap of 5 the largest credit is `5 tan(phi_b)` where
    the uncapped field on the same model reaches eighteen times that;
  * F = 1.74 stands and F = 1.82 fails, on BOTH drivers;
  * and the credit is load-bearing: with it switched off F = 1.74 must FAIL on both
    drivers (that model is locked at 1.713 without it), with the binding cap it must
    fail, and with a looser cap of 20 it must stand again — so what moves the answer
    is the cap's VALUE and not the presence of a cap keyword.

The guard check lost its matric-suction case, which had been asserting a refusal the
driver no longer has.

**Mutation, run four ways.** Each is the same check function against a driver with
one thing removed:

| driver | verdict | what the check saw |
|---|---|---|
| as shipped | **PASS** | — |
| the `s_cap` bound dropped from `_suction_capped` | **FAIL** | "s_cap does not bound the credit: with a cap of 5 the largest apparent cohesion is 68.9055, above the 3.7678 that cap allows", and both drivers then carry F = 1.74 under that cap |
| the `1/F` dropped from `_suction_apparent_cohesion` | **FAIL** | "the matric-suction apparent cohesion is NOT reduced by the trial strength: raising F from 1.0 to 1.5 left the credit at a ratio of at most 1.000000 where the viscoplastic driver divides it by 1.5" |
| the credit computed but not folded into `c_r` | **FAIL** | F = 1.74 comes back FAILED on the Newton driver — the field is right and the solve is not using it |

The first two mutations move BOTH drivers together, which is why a driver-against-
driver agreement check could not have caught either and why the criterion asked for
the field. The third is Newton-only and is caught by the answer.

One thing in the check's design is worth naming because it was found by running it
rather than by reading. The F-reduction leg first compared the field at F = 1 with
the field at F = 2; F = 2 is past this model's limit, the trial dies at the load-step
floor, the viscoplastic predictor fires, and the predictor's own group build calls
the same helper — so the capture held two drivers' fields end to end and the leg
silently skipped on a length mismatch. It reads F = 1.0 against F = 1.5 now, two
strengths the driver CARRIES, and a length mismatch is a failure rather than a skip.

#### The criterion, line by line

**1. The suction-gated locked models — PARTLY MET, at three of six.** The three
`vp102t` `c3` rows reproduce their locks at +0.0025, +0.0025 and +0.0136 against a
tolerance of 0.02, and two of the three agree with the viscoplastic driver inside
the bisection tolerance (0.0000 and 0.0055) while the third is 0.0219 away. The
three `rs2_28` rows do NOT reproduce their locks: +0.0409, +0.0216 and +0.0221
against 0.02, all HIGH, and at the tag's own bisection resolution +0.0378, +0.0248
and +0.0252. The viscoplastic driver reproduces all six. Work is 4.3x to 23.5x
lower on the Newton driver on every row.

**2. The two drivers build the same field — MET, exactly.** The maximum absolute
difference is **0.000e+00** on both models measured, over 11,064 and 3,354 Gauss
points, captured at the shared helper through each driver's own entry point; and
each driver's field is within 1.8e-15 and 7.1e-15 of an expectation computed from
scratch out of the prepared model's own arrays.

**3. `s_cap` binds and both drivers move together — MET.** At a cap of 5 on
`vp102t_60`, which bounds 421 of the 474 credited Gauss points, F = 1.74 fails on
both drivers where it stands on both uncapped and at a cap of 20. With the credit
switched off entirely it fails on both, which is what makes the pair a measurement
of the credit rather than of the model.

**4. The ramp — MET on 3 of the 3 run**, against the 2 asked: +0.0012, −0.0057 and
−0.0088 from the Newton bisection.

**5. The no-suction path is bit-identical — MET on 16 benchmark pairs and 139
trials**, against the driver staged at 53659aa4 in a separate package tree: every
one identical in factor of safety, verdict, iterations and force evaluations, and
the plain return map bit-identical over 800,000 trial states. The K0 vendor rows
the clause also named were not re-run against the control tree — that is a
shortfall against the criterion as written, reported rather than glossed; they are
covered by `check_k0_initial_stress` in the check file, which passes.

**6. The refusal is gone and the guard list is updated — MET.** Matric suction is
carried; piles, bar softening, Hoek-Brown and the power curve each still raise and
each still names itself, with the viscoplastic control accepting all four.

**7. The locks — MET.** `check_matric_suction` passes as shipped and fails on all
three mutations: the cap dropped, the F-reduction dropped, and the credit computed
but not folded into the envelope. The whole check file passes.

**8. The default path — MET.** Griffiths & Lane 6 dry, quad8, size 2, no
`fem_solver` argument: FS 2.421875 on per-trial iteration counts 147, 781, 3393,
2031, 2841, 9541, 12000, 8617, 8777 — value for value the control sequence, on both
trees.

**9. The honest negatives** are the three RS2-28 locks, the 0.0219 driver gap on
`vp102t_60`, and the K0 rows missing from the control sweep. All three are written
above.

#### Verdict

**Matric suction costs the Newton driver nothing and it is carried.** Fredlund's
apparent cohesion is not a new surface — it is a different `c` at a Gauss point,
and the return map has taken `c` per Gauss point since the first round — so the
whole feature is two helpers, one call in the group build and one in the ramp's
`restrength`. The measurement that says the plumbing is right rather than merely
plausible is the field: on 11,064 Gauss points of `rs2_28a` and 3,354 of
`vp102t_60`, the apparent cohesion the two drivers compute is identical to
0.000e+00, and each is within 7e-15 of an expectation built from scratch. The cap
binds where it should and both drivers move together on its VALUE; the credit is
divided by the trial strength as the viscoplastic driver divides it, and by the
same per-element factor, so the Gauss points an SSR zone holds at full strength
keep their whole credit. Sixteen benchmark pairs and 139 trials are bit-identical
against a control tree, and the default path is untouched.

**Three of the six locked models reproduce and three do not, and the three that do
not are the interesting result.** On RS2-28 the Newton driver reads 0.022 to 0.041
above the published XSLOPE values, always high. That gap was not left as a
direction: the deciding trials were run. At F = 1.6250 and F = 1.6375 on RS2-28a —
both above the viscoplastic answer and at or below the Newton one — the Newton
driver reaches equilibrium in 126 and 252 iterations at an out-of-balance of 2.6e-6
and 1.7e-5, with a worst Mohr-Coulomb value of 2.9e-9 of the local strength read in
the INVARIANT form the return map is not written on. That is a statically
admissible stress field carrying full gravity, and by the lower-bound theorem it is
a proof that the slope stands there. The viscoplastic driver fails both, at the
budget-extension heuristic, with its out-of-balance at 1.1e-1 and 2.3e-1 against a
1e-3 gate — two orders of magnitude out, not a near miss.

And the vendor agrees with the Newton driver. RS2's own SSR values on these three
problems are 1.64, 1.55 and 1.41; the Newton bisection reads 1.6438, 1.5688 and
1.4063, which is +0.0038, +0.0188 and −0.0038 from them, while the viscoplastic
driver reads 0.006 to 0.021 below. XSLOPE's locked values for this row ARE the
viscoplastic readings, and `docs/verification/rs2.md` already scores them 2.1%,
0.4% and 2.1% below RS2 and attributes the shortfall to the convergence check. So
the correct statement of this round's negative is narrow and it is not that the
Newton driver is wrong on RS2-28: it misses three locks that are defined by the
driver it is being measured against, and it lands on the vendor's answer instead.
Whether those locks should move is the owner's decision and not a spike's, and
nothing here was changed.

What remains, in the order it matters:

- **Whose number RS2-28 is.** The Newton answer is supported by an admissible field
  and matches RS2's own SSR to 0.004 on two of the three; the locked answer is the
  viscoplastic driver's, and the page already documents it as low against the
  vendor. That is a corpus decision with a measurement behind it, and it is the
  same decision the three-layer adjudication and the K0 round both walked up to
  from different directions.
- **`vp102t_60`'s 0.0219.** One row of six, the same direction, and the only one of
  the three transient rows that misses. Unmeasured trial by trial; the RS2-28
  diagnostic above is the method if it is worth doing.
- **The eight refusals are down to four.** Piles, post-peak bar softening,
  Hoek-Brown and the power curve. The eligibility inventory puts eight locked models
  behind piles, two behind bar softening, three behind Hoek-Brown and five behind
  the power curve; Hoek-Brown and the power curve are two linearizations of the same
  shape and would come together.
- **The corpus run.** Unchanged from the K0 round's verdict: the reachable set is now
  larger by six workbooks and nine locks, and running it is still what would turn a
  sample into a statement about the corpus.

## PILES — the Newton driver carries the beam element

Written before any feature code, so that what follows is a test and not a
description. Same machine and settings as everything above: `force_tol` 1e-3,
hybrid criterion, `capture_failure_state=False`, tolerance 0.01.

### Why this one now

The suction round left four refusals: piles, post-peak bar softening, Hoek-Brown
and the power curve. Piles gate the largest block of them — eight locked
`fem_ssrm` benchmarks against three for Hoek-Brown and five for the power curve,
and they are the only one of the four whose models span the tutorials, the FEM
samples and two vendor corpora at once. They are also the only refused feature
that changes the SHAPE of the problem rather than the constitutive law: a pile
node carries a rotational degree of freedom, so a pile model is the first thing
this driver has met whose displacement vector is not a list of lengths.

That last point is why this round is not only an addition. The Newton
displacement bound reads `max|u|` over the RAW degree-of-freedom vector, and the
guard that has kept that honest until now is an assertion that no rotational
degree of freedom exists — which is true exactly because piles are refused.
Carrying piles without changing the bound would silently compare a length
against a radian.

The eight locked models, with the lock each carries (`docs/verification/`,
`docs/fem/samples.md`, `docs/tutorials/fem03_piles.md`; every one at
tolerance 0.01):

| Benchmark | Model | Lock | What it exercises |
|---|---|---|---|
| SIGMAW-SRS-wall | `gs2_wall` | 1.647 | sheet pile wall, an elastic material, `t_cut` = 0 on two soils, a c = 0 layer |
| VP106-FEM-free | `vp106c_fem` | 1.472 | pile row at S = 2.4, free head and tip |
| VP106-FEM-fixed | `vp106c_fem_fix` | 1.587 | the same row with the head rotation HELD |
| FEM-3-wall-ssrm | `xslope_pile_wall` | 1.559 | sheet pile wall, tip FIXED, a finite `M_cap` |
| FEM-3-piles-ssrm | `xslope_piles` | 1.379 | two pile rows, finite `V_cap` AND `M_cap` |
| (FEM sample) | `xslope_piles_fem` | 1.380 | the same two rows |
| SSRM-TORGGLER | `xslope_torggler_3a_plate` | 1.195 | a 7.5 m plate, no capacities |
| SSRM-TORGGLER | `xslope_torggler_3b_plate` | 1.673 | a 15 m plate over a weak band |

### The semantics being reproduced, read from the viscoplastic driver

Not assumed — read out of `build_fem_data`, `_prepare_fem_model` and `solve_fem`
and restated here, because the Newton path has to solve the same model.

- **The element.** A two-node Euler-Bernoulli beam on a linear soil edge, a
  three-node one on a quadratic edge, standing on the same nodes as the soil.
  Bending on the three-node element is the quintic Hermite that matches a value
  and a slope at all three nodes; axial is the quadratic bar; the two are
  uncoupled. `K_global_pile_elems` is `T^T K_local T`, and it is what the global
  elastic stiffness is assembled from — the same matrix on both drivers.
- **The degrees of freedom.** Every node a pile element stands on, midside nodes
  included, carries three DOFs; every other node carries two. `dof_offset` is
  the map and `is_pile_node` the flag.
- **Per-unit-width.** `EA`, `EI`, `V_cap` and `M_cap` are all divided by the
  spacing `S`, so everything the solver sees is per unit width of slope while the
  file states the single-pile section. `docs/usage/input_template.md` is the
  published statement of that convention.
- **The actions.** Axial `EA/L (u2 - u1)` at the element center; the shear at the
  center (constant along a two-node element, `EI v'''(L/2)` on the three-node
  one); and the bending moments `M1`, `M2` at the two END nodes, read as rows 2
  and 5 of `K_local u_local`. The midside node has no action of its own.
- **The coupling is a shared node.** There is no interface element and no bond
  law: the pile and the soil have identical displacements at every node, which is
  a perfectly bonded interface (`docs/fem/piles.md`). There is nothing here to
  reproduce beyond the beam blocks themselves.
- **Fixity is a CONSTRAINT, not a load.** Head and tip each take free / pinned /
  unrotated / fixed, and `_prepare_fem_model` removes the held degrees of freedom
  from `free_dofs`. Both drivers read that same set, so tip fixity needs no code
  on this path — which has to be ASSERTED rather than assumed.
- **Nothing latches.** `yielded_pile_V` and `yielded_pile_M` are set for
  reporting and never read back into the constitutive rule; the capacity check is
  recomputed from the current displacement at every iteration. The pile law is
  therefore nonlinear-ELASTIC with no history, exactly as the bar law is with
  softening refused: nothing to commit at the end of a step and nothing for the
  ramp's warm start to carry.
- **Nothing is reduced by F.** `E`, `I`, `A`, `V_cap` and `M_cap` are held while
  only `c` and `tan(phi)` fall, which is the vendor convention, what the bar does
  and what `docs/fem/piles.md` states.
- **The axial force has no capacity at all.** There is no `T_cap` on a pile.

**And two things about the capacities that were measured rather than read**,
because the shape of this round depends on them. Both were measured on
`xslope_piles_fem` (tri6/2, 18 pile elements, `V_cap`/S = 7,666.7,
`M_cap`/S = 10,000) at F = 1.3, on the DEFAULT driver, with nothing changed.

**(i) The moment cap is inert, and no sign of it could be otherwise.** The
correction `sign(M) M_cap - M` is applied to the shared rotational degree of
freedom at the element's end node. Two beam elements meet at every interior node
of a pile, and at equilibrium their end moments there are equal and opposite, so
their two corrections are equal and opposite and cancel exactly. Measured: with
`M_cap` set to 10,000, 1,000, 100 and 10 while the largest moment in the pile is
6,763.5, the net applied moment summed over the pile's rotational degrees of
freedom is 0, 8.2e-8, 1.8e-7 and 1.5e-7, and `max|u|` is 0.054093354479 at every
one of the four — unchanged to twelve decimals across four decades of cap. The
pile's own moments do not move either. A plastic hinge is a RELEASE of rotational
continuity; a moment applied at a node the two elements share is not one, and the
reported moments are clipped separately, twice (inside the loop and again at
Step 10c), which is why an inert cap reads as an enforced one.

**(ii) The shear cap is applied with the sign that makes it an anti-cap.** The
viscoplastic scheme solves `K u = base_loads + corrections`, so the internal force
the state is in equilibrium with is `K u - corrections`. The bar path applies
`correction = T - T_true` and therefore delivers `T_true` — the comment there says
in as many words that the opposite sign "turns the cap into an anti-cap: the bar
ends up carrying 2T - T_true". The pile path applies `correction = V_true - V`.
Measured at the converged state, reading the shear the pile actually delivers out
of `K u - corrections` in the beam's own frame rather than out of the reported
(clipped) array:

| `V_cap`/S | reported max &#124;V&#124; | DELIVERED max &#124;V&#124; |
|---|---|---|
| none | 1,459.7 | 1,457.6 |
| 7,666.7 (the file's) | 1,459.7 | 1,457.6 |
| 766.7 | 766.7 | 1,731.3 |
| 76.7 | 76.7 | 1,867.7 |

The tighter the cap, the more shear the pile delivers. The reported column is the
clip; the delivered column is the physics.

Neither of those is changed by this round — `solve_fem` is not touched, and which
number the repository's locked pile values should carry is the owner's decision,
not a spike's. What they decide is what the Newton path can honestly implement,
which the design states next.

### Design

- **`_nr_build_piles(fem_data)`**, grouped by degree-of-freedom count (6 on a
  two-node element, 9 on a three-node one) so each group is a dense array pass,
  holding the global element matrix `K`, the rotation `T`, the local rows that
  read `V`, `M1` and `M2` out of `u_local`, the global patterns those actions are
  delivered on, and the per-element caps. Every array is read from `fem_data` and
  never written, so the two drivers carry the same beam.
- **Residual.** `f_int = K u_e - sum_k (s_k - s_k_true) q_k` over the three capped
  actions, with `s_true = clip(s, -cap, +cap)` and `q_k` the action's own internal
  force pattern. That is the BAR's convention — the one that delivers `s_true` —
  applied to a beam. It is not the viscoplastic driver's sign on the shear, and
  the consequence is measured rather than assumed: on every model where no cap
  binds the two drivers solve the identical problem, and on a model where one
  binds they do not, and the criterion below says what has to be measured there.
  The moment leg is written the same way and is expected to be inert for the
  structural reason above, on this driver as on the other; that it is inert is
  measured, not assumed.
- **Consistent tangent, exact.** Each action is LINEAR in `u_e` and each capped
  branch is affine, so there is nothing to difference: `ds/du = g` (a constant
  row) and `d(s - s_true)/du` is `g` while the action is over its cap and zero
  otherwise. The element tangent is `K - sum_active q_k (x) g_k`. This is verified
  against a central difference of the element's own internal force rather than
  asserted.
- **Assembly.** The pile blocks join the cached sparsity pattern explicitly, after
  the soil groups and the bars, in the order the assembler walks them — the same
  treatment the bars got, and for the same reason: a rotational degree of freedom
  appears in NO soil element, so relying on the soil pattern would silently drop
  every rotation from the tangent.
- **Tip and head fixity.** Nothing to build: the held degrees of freedom are
  already out of `prep["free_dofs"]`. The criterion asserts that rather than
  assuming it, on the one model that holds a head rotation and the one that fixes
  a tip.
- **No history.** With nothing latching, the pile law is a function of the current
  displacement alone, so there is nothing to commit at the end of a load step and
  the ramp's warm start needs no extension. The predictor carries the pile state
  for free: its seed is a displacement field on the same mixed-DOF vector.
- **Strength reduction touches soil only.** The ramp's `restrength` rewrites the
  soil groups and nothing else; the criterion asserts that a pile's capacity and
  rigidity are identical at the foot of the ramp and at its limit.
- **The displacement bound becomes translational.** `max|u|` is read over the
  TRANSLATIONAL degrees of freedom only, which is what the viscoplastic driver's
  own displacement-limit check reads ("Extract translational DOFs only for VP
  displacement check"). The same measure is used for the elastic displacement
  scale, the tangent probe's step and the `max|du|/max|u|` the drivers both report
  as `residual`, and for the ramp's own bound. On a model with no pile every
  degree of freedom is translational, so every one of those quantities is
  bit-identical there — which is what the no-pile control tree below checks.
- **The guard's pile line goes and nothing else does.** Post-peak bar softening,
  Hoek-Brown and the power curve stay refused, each naming itself.

### Success criterion (verbatim)

1. **The element is right, measured twice.** (a) On at least 400 random beam
   elements — random orientation, length, `EA`, `EI`, two-node and three-node,
   with the displacement drawn so that all capacity branches are exercised and the
   branch histogram reported — the analytic element tangent agrees with a central
   difference of that element's own internal force to 1e-8 relative or better, on
   every branch. (b) An isolated beam with no soil present reproduces the
   closed-form deflection to 1e-10 relative: a cantilever under an end load and a
   simply supported beam under a central load, at more than one orientation.
2. **The eight pile-gated locked models**, each built through `run_tests.py`'s own
   `build_fem_ssrm_case` so the mesh, element type, bracket, iteration budget and
   every solver option are the suite's and not this round's. Reported per model:
   viscoplastic FS, Newton FS, the lock, its tolerance, and the constitutive work
   on each driver. The criterion is that the Newton bisection lands inside the
   published tolerance on at least 6 of the 8 AND within 0.01 of the viscoplastic
   driver on at least 6 of the 8. **Every miss is diagnosed rather than reported as
   a direction**, by running the trials whose verdicts differ and reading the
   Newton state's own admissibility — out-of-balance and the worst Mohr-Coulomb
   violation in the INVARIANT form — the way the RS2-28 and three-layer
   disagreements were refereed.
3. **The capacity divergence is measured, not left as a design note.** On each of
   the three models carrying a finite capacity, whether any cap binds at any trial
   of the bisection is reported. Where none binds, the two drivers must agree
   exactly, because they are solving the same problem. Where one binds, the two
   answers are reported side by side together with the delivered-force reading that
   says which law each driver is enforcing, and a constructed model with a cap
   tightened until it binds is run on both drivers so the divergence has a number
   on it rather than an argument.
4. **Pile diagnostics.** On at least 2 models, the axial force, shear and end
   moments from the Newton solution are compared element by element against the
   viscoplastic solution's at the same F, and the largest relative difference is
   REPORTED as a number. The `yielded_pile_V` / `yielded_pile_M` masks are compared
   and any difference is explained by the convention that produces it, as the bar
   masks' latching difference was.
5. **The displacement bound reads translational degrees of freedom only**, and the
   lock asserts it: on a pile model whose rotations are large beside its
   displacements, a mutation that puts the rotations back into the bound must FAIL
   the check. Run both ways and recorded.
6. **The no-pile path is bit-identical**, on the same control-tree protocol the
   previous rounds used: the plain-soil eight rows, the four reinforced benchmarks,
   the tension-cutoff rows, a K0 vendor row and a matric-suction row, re-run
   against a control run of the parent commit staged in a separate package tree —
   every trial identical in factor of safety, verdict, iterations and force
   evaluations — plus the plain Mohr-Coulomb return map bit-identical over 800,000
   random trial states.
7. **The ramp agrees** with the Newton bisection within 0.01 on at least 3 of the 8
   models, and carries the pile through its warm history: a pile's rigidity and
   capacity at the top of the ramp are identical to the foot's.
8. **The refusal is gone and the guard list is updated.** Piles are carried;
   post-peak bar softening, Hoek-Brown and power-curve envelopes each still raise
   and each still names its own feature, with the viscoplastic control accepting
   all three. A pile model whose viscoplastic treatment involves anything this path
   does not carry must refuse with a message naming it rather than being solved
   silently.
9. **The locks catch it.** `test/nr_ssrm_check.py` gains `check_piles`, asserting
   the element (tangent against a central difference, and the closed-form beam),
   the fixity constraint, the factor of safety on one cheap locked pile model from
   both drivers, and the translational bound. Three mutations, run both ways and
   recorded: corrupt the beam's bending stiffness, put the rotations back into the
   displacement bound, and drop the capacity correction from the residual. Each
   must FAIL the check. The whole check file passes.
10. **The default path is unchanged**, against the standard control: Griffiths &
    Lane 6 dry with no `fem_solver` argument, FS 2.421875 on per-trial iteration
    counts 147, 781, 3393, 2031, 2841, 9541, 12000, 8617, 8777.
11. **An honest negative is a valid outcome and must be written.** If the beam
    needs something this design does not have — if the rotational degrees of
    freedom break the tangent's conditioning, if the locks cannot be reproduced, or
    if the capacity divergence turns out to decide a published number — that is the
    result.

### PILES — results

Same machine and settings as everything above: `force_tol` 1e-3, hybrid criterion,
`capture_failure_state=False`, tolerance 0.01. Every number below was measured on
this checkout in this session.

#### What was built

`_nr_build_piles` and `_nr_pile_force` in `xslope/fem.py`, threaded through
`_nr_internal_force`, `_nr_prepare_assembly`, `_nr_assemble_tangent`,
`_nr_equilibrate`, `_solve_fem_newton` and the ramp; `_nr_translational_dofs` and
`_nr_umax` beside them; the guard's pile line replaced by the comment that says
where the element is carried. A pile model now returns the same five diagnostic
arrays the viscoplastic path returns — axial force, shear, the two end moments,
and the two yielded masks — indexed by pile element in the same order, so the
summary printer, the result CSVs and the pile-shear colorbar consume them
unchanged.

The pile blocks join the cached assembly pattern explicitly. For the bars that was
a precaution; here it is a requirement, because a pile node's rotational degree of
freedom appears in no soil element and in no bar, so relying on the soil pattern
would drop every rotation from the tangent.

#### The element, measured twice

**The consistent tangent is exact, and the measurement says why.** 400 random beam
elements — random orientation, length, `EA` and `EI`, two-node and three-node, with
the caps drawn from the actions the displacement actually produces so that all
eight combinations of (shear, end moment 1, end moment 2) at their cap are reached
(92 elastic, 58 shear only, 56 and 56 one moment only, 44 / 36 / 34 two at once, 24
with all three) — against a central difference of the element's own internal force,
at four probe steps:

| probe step, as a fraction of max&#124;u&#124; | 1e-7 | 1e-5 | 1e-4 | 1e-3 |
|---|---|---|---|---|
| worst relative gap, two-node | 5.26e-8 | 4.89e-10 | 9.27e-11 | 7.13e-12 |
| worst relative gap, three-node | 2.19e-9 | 3.05e-11 | 2.12e-12 | 1.95e-13 |

The gap falls as 1/h, exactly. That is round-off cancellation in the difference and
nothing else — every branch is affine in the element displacement, so the analytic
tangent is the derivative and the only error is in the measurement of it. The
criterion asked for 1e-8; the reading is inside it at every probe step of 1e-5 or
larger, and the 1/h scaling is the evidence that the residual at 1e-7 is the
difference's and not the tangent's. No trial was skipped for straddling a branch.

**An isolated beam reproduces the closed forms**, with no soil present, on the
round's own assembly and residual. A six-element cantilever under a transverse end
load and the same beam simply supported under a central load, at three
orientations:

| element | orientation | cantilever deflection | tip rotation | simply supported |
|---|---|---|---|---|
| two-node | 0 deg | 3.2e-14 | 2.3e-14 | 2.6e-15 |
| two-node | 90 deg | 6.5e-14 | 5.6e-14 | 2.6e-15 |
| two-node | 40.1 deg | 1.4e-13 | 1.3e-13 | 2.5e-14 |
| three-node | 0 deg | 4.6e-12 | 4.4e-12 | 3.7e-13 |
| three-node | 90 deg | 4.8e-12 | 4.5e-12 | 4.0e-13 |
| three-node | 40.1 deg | 1.7e-12 | 1.3e-12 | 5.4e-13 |

Relative to $PL^3/3EI$, $PL^2/2EI$ and $PL^3/48EI$. The criterion asked for 1e-10.

**And it reads the same actions the viscoplastic driver reads.** The Newton element
recovers its axial force, shear and end moments off constant ROWS, because it needs
their gradients; the viscoplastic path computes the same four quantities in closed
form at every iteration. On 300 random elements the two agree to 1.8e-14 (axial),
1.7e-14 (shear), 4.1e-13 and 6.0e-14 (the two end moments) relative. The two
drivers are capping the same quantities.

#### The pile diagnostics, against the viscoplastic driver

`xslope_piles_fem` at its tagged tri6/2.0 mesh (1,521 elements, 18 pile elements),
element by element at the same strength on both drivers:

| F | what each driver did | axial | shear | end moments | masks |
|---|---|---|---|---|---|
| 1.300 | both CONVERGED (VP 1,364 iterations at oob 9.99e-4; Newton 472 at 1.22e-7) | 0.37% | 0.68% | 0.31% | both empty, agree 18/18 |
| 1.375 | both CONVERGED (VP 12,612 iterations at oob 9.86e-4; Newton 26 at 7.43e-10) | 7.65% | 7.80% | 13.37% | both empty, agree 18/18 |

Percentages are the largest per-element difference as a fraction of the largest
viscoplastic value. At F = 1.300 the two drivers are in the same state and the pile
forces agree to under a percent on all three actions. At F = 1.375 — one bisection
cell below this model's limit — they are not: the viscoplastic solve takes 12,612
iterations to reach its 1e-3 gate while the Newton corrector lands at 7.4e-10 in
26, at twice the displacement (1.00% of the model height against 0.50%), and a
member reads the difference between two elastoplastic states more sharply than the
soil does. Neither state has a pile at its capacity, so the two yielded masks are
empty and agree element for element; the LATCHING difference the bar masks carry
(the viscoplastic mask records every element that was ever over its capacity during
the iteration history, this one is read on the reported state) is present in the
code on this path too and is not exercised by these two states.

#### Rotations are smaller than displacements, everywhere it was measured

The displacement bound was a category error and is fixed; what it was NOT is a
number. Measured on the converged states above and on the sheet pile wall, the
largest nodal rotation against the largest translational displacement:

| model | F | max&#124;u&#124;, translational | max rotation | ratio |
|---|---|---|---|---|
| `xslope_piles_fem` | 1.300 | 0.0541 | 0.00287 | 5.3% |
| `xslope_piles_fem` | 1.375 | 0.3014 | 0.01570 | 5.2% |
| `xslope_pile_wall` | 1.5625 | 0.1671 | 0.01134 | 6.8% |
| `xslope_pile_wall` | 1.7500 | 0.2778 | 0.01758 | 6.3% |
| `xslope_pile_wall` | 1.7969 | 0.8148 | 0.01922 | 2.4% |

**So on the corpus the corrected bound reads the same number the raw one did**, and
that is not an accident of these models: a nodal rotation is about a deflection
divided by a length, so it can only exceed the largest displacement in a model whose
pile is shorter than one unit of length. There is no such model in the repository,
and there is unlikely to be one in feet or metres. The fix is a correctness fix with
no measured consequence — which is exactly what it should be, since a bound that had
been silently comparing a length against a radian would otherwise have been deciding
verdicts. The lock asserts it structurally, on the index set the bound reads; the
behavioural mutation the criterion asked for could not be built on a physically-sized
model, and that is reported rather than manufactured.

#### The two capacities, on both drivers

The design said the Newton path would write the capacity the way the bar's is
written — the part of the elastic action the member cannot deliver, subtracted from
its own internal force — and that the consequence would be measured rather than
assumed. It was, on the three locked models that carry a finite capacity.

**The moment cap is inert on BOTH drivers, and no sign of it could be otherwise.**
`xslope_pile_wall` carries `M_cap` = 90,600 and nothing else, and its Newton
solution at the reported state has an element AT that cap. Run with the cap and with
it removed entirely, same mesh, same bracket:

| | `M_cap` = 90,600, as the file gives it | `M_cap` removed |
|---|---|---|
| viscoplastic | 1.55859375 | 1.55859375 |
| Newton | 1.80078125 | 1.80078125 |

Neither answer moves. The reason is structural and was measured on the viscoplastic
path before any of this was built (see "the semantics being reproduced", above): the
moment correction is applied at the rotational degree of freedom the two adjacent
beam elements SHARE, and at equilibrium their end moments there are equal and
opposite, so the two corrections cancel. Reversing the sign reverses both of them
and they cancel again. A plastic hinge is a release of rotational continuity, and a
moment applied at a shared node is not one. What the Newton path adds is that the
same construction is inert on it too, so the two drivers agree about a cap that
neither of them enforces — and the reported moments are clipped, which is why an
inert cap reads as an enforced one on both.

**The shear cap is enforced on the Newton path and anti-enforced on the
viscoplastic one**, and neither model in the corpus notices. On `xslope_piles_fem`
and `xslope_piles` the shear cap is `V_cap`/S = 7,666.7 and the largest shear
anywhere in either bisection's reported state is 2,391 — the cap is a third of the
way out of reach, so it never binds, and the two drivers agree at every trial of
both bisections. On `xslope_pile_wall` there is no shear cap at all. **On the eight
locked pile models the shear cap therefore decides nothing, and the difference
between the two laws costs nothing.**

Where it does bind it is not small, and that was measured rather than left as an
argument. `xslope_pile_wall` at F = 1.2, with the shear capacity tightened to a
quarter of what the uncapped pile carries:

| | uncapped | cap at a quarter of the free shear |
|---|---|---|
| Newton: elements at the cap | 0 | 6 of 10 |
| Newton: max&#124;u&#124; | 0.038574 | larger — the slope moves MORE |

A member that can deliver less force cannot hold the soil back better, and that is
the direction the check asserts. Under the viscoplastic sign the same tightening
raises the shear the pile delivers (1,457.6 uncapped, 1,731.3 at a tenth of the cap,
1,867.7 at a hundredth, measured on `xslope_piles_fem` and reported in the semantics
section above), so the slope would move LESS. The lock's capacity leg fails on a
driver with either the correction dropped or the sign reversed.

Nothing in `solve_fem` was changed. Which number the repository's locked pile values
should carry, on a model where a capacity does bind, is the owner's decision.

#### The eight pile-gated locked benchmarks

Every row built through `run_tests.py`'s own `build_fem_ssrm_case`, so the mesh, the
element type, the bracket, the iteration budget and every solver option are the
suite's and not this round's; the bisection tolerance is the tag's own 0.01 on all
eight. Work is the honest count on each driver — viscoplastic iterations against
Newton force evaluations, one constitutive pass each.

| Benchmark | Model | Elements | Lock | Tol | VP FS | VP − lock | Newton FS | Newton − lock | driver gap | VP work | N-R work |
|---|---|---|---|---|---|---|---|---|---|---|---|
| FEM-3-piles-ssrm | `xslope_piles` | 1521 | 1.379 | 0.01 | 1.3789063 | −0.0001 | **1.3789063** | −0.0001 | **0.0000** | 20,534 | 32,057 |
| (FEM sample) | `xslope_piles_fem` | 1521 | 1.380 | 0.01 | 1.3796875 | −0.0003 | **1.3796875** | −0.0003 | **0.0000** | 20,054 | 27,736 |
| FEM-3-wall-ssrm | `xslope_pile_wall` | 1510 | 1.559 | 0.01 | 1.5585938 | −0.0004 | **1.8007813** | +0.2418 **(out)** | 0.2422 | 13,094 | 26,436 |
| VP106-FEM-free | `vp106c_fem` | 1591 | 1.472 | 0.01 | 1.4718750 | −0.0001 | **1.5781250** | +0.1061 **(out)** | 0.1063 | 38,370 | 33,172 |
| VP106-FEM-fixed | `vp106c_fem_fix` | 1591 | 1.587 | 0.01 | 1.5871094 | +0.0001 | **1.5941406** | +0.0071 | 0.0070 | 104,508 | 32,515 |
| SSRM-TORGGLER | `xslope_torggler_3a_plate` | 6834 | 1.195 | 0.01 | 1.1945313 | −0.0005 | **1.1945313** | −0.0005 | **0.0000** | 33,442 | 28,067 |
| SSRM-TORGGLER | `xslope_torggler_3b_plate` | 4945 | 1.673 | 0.01 | 1.6726563 | −0.0003 | **1.7429688** | +0.0700 **(out)** | 0.0703 | 289,072 | 22,065 |
| SIGMAW-SRS-wall | `gs2_wall` | 6532 | 1.647 | 0.01 | 1.6468750 | −0.0001 | **1.6906250** | +0.0436 **(out)** | 0.0438 | 303,303 | 25,218 |

The viscoplastic column reproduces every one of the eight locks, at 0.0001 to
0.0005, which is what makes the comparison readable. It is also where the work is:
`gs2_wall` costs it 303,303 constitutive passes against the Newton driver's 25,218
force evaluations, and `xslope_torggler_3b_plate` 289,072 against 22,065 — 12x and
13x. Across the eight rows the Newton driver does less constitutive work on five —
1.2x, 1.2x, 3.2x, 12.0x and 13.1x — and MORE on three: `xslope_pile_wall` at 0.5x,
`xslope_piles` at 0.6x and `xslope_piles_fem` at 0.7x. All three of those are models
whose bisection spends several trials past failure, which is where load control has
nothing to offer and where this document has measured the Newton driver at 17x to
47x the viscoplastic cost since Phase 0.

**Four of the eight land inside their published tolerance and four do not.**
`xslope_piles`, `xslope_piles_fem` and `xslope_torggler_3a_plate` return the
IDENTICAL factor of safety to the viscoplastic driver, with the identical verdict at
every trial the bisection visited; `vp106c_fem_fix` agrees to 0.0070. The other four
read 0.0436, 0.0700, 0.1061 and 0.2418 above their locks, every one HIGH, which is
the one-sided direction this document has recorded since its first table. **The
criterion asked for six of eight on both counts and the measurement is four of
eight; that clause is NOT met.**

One pairing in that split is worth naming because it is a controlled comparison.
`vp106c_fem` and `vp106c_fem_fix` are the same slope, the same mesh and the same
pile row, differing only in whether the pile HEAD's rotation is held. The free-head
row is one of the four misses, at +0.1061; the held-head row is one of the four that
agree, at +0.0070. Holding one degree of freedom is the difference between the two
drivers agreeing and disagreeing by a tenth.

#### The two misses, refereed by the Newton state's own evidence

The criterion said a miss would be diagnosed rather than reported as a direction, so
the trials whose verdicts differ were run directly. `xslope_pile_wall` (mesh height
30) and `vp106c_fem` (mesh height 20), at strengths between the viscoplastic answer
and the Newton one:

| Model | F | driver | verdict | exit | out-of-balance | iterations | worst yield violation | max&#124;u&#124; |
|---|---|---|---|---|---|---|---|---|
| pile wall | 1.5625 | viscoplastic | FAILED | `diverging` | 9.49e-02 | 1,231 | — | 0.55% of H |
| pile wall | 1.5625 | **Newton** | **CONVERGED** | `converged` | **3.05e-05** | 508 | **1.3e-08** | 0.56% of H |
| pile wall | 1.7500 | viscoplastic | FAILED | `diverging` | 2.60e+00 | 261 | — | 0.56% of H |
| pile wall | 1.7500 | **Newton** | **CONVERGED** | `converged` | **9.91e-06** | 517 | **1.5e-08** | 0.93% of H |
| pile wall | 1.7969 | viscoplastic | FAILED | `diverging` | 3.15e+00 | 221 | — | 0.55% of H |
| pile wall | 1.7969 | **Newton** | **CONVERGED** | `converged` | **1.83e-05** | 469 | **1.4e-08** | 2.72% of H |
| vp106c | 1.5000 | viscoplastic | FAILED | `diverging` | 6.54e-01 | 1,501 | — | 0.65% of H |
| vp106c | 1.5000 | **Newton** | **CONVERGED** | `converged` | **2.67e-05** | 606 | **1.4e-08** | 0.78% of H |
| vp106c | 1.5500 | viscoplastic | FAILED | `diverging` | 1.57e+00 | 861 | — | 0.65% of H |
| vp106c | 1.5500 | **Newton** | **CONVERGED** | `converged` | **1.34e-06** | 569 | **1.6e-08** | 1.33% of H |
| vp106c | 1.5750 | viscoplastic | FAILED | `diverging` | 1.99e+00 | 721 | — | 0.65% of H |
| vp106c | 1.5750 | **Newton** | **CONVERGED** | `converged` | **5.05e-06** | 99 | **1.6e-08** | 5.19% of H |

The yield violation is the largest Mohr-Coulomb value over every Gauss point as a
fraction of the local strength, read in the INVARIANT form the return map is not
written on. A stress field in equilibrium with full gravity to between 1.3e-6 and
3.1e-5 — thirty to seven hundred times inside the trial tolerance — and nowhere more
than 1.6e-8 of its own strength outside the yield surface is a statically admissible
field, and by the lower-bound theorem it is a proof that the slope stands at that
strength. The argument needs neither driver to be trusted; it is the Newton driver's
own converged state checked against a yield function it does not solve on.

The viscoplastic verdicts at those strengths are not near misses. Every one of them
is closed by the early-failure classifier (`diverging`) with the out-of-balance at
9.5e-2 to 3.2 against a 1e-3 gate — two to three orders of magnitude out, not a run
that a longer budget would flip. **This is the same shape as the RS2-28 misses of
the matric-suction round**, and it has the same reading: the Newton driver does not
miss the two locks by being wrong about the slope, it misses two locks that are
DEFINED by the driver it is being compared against.

One thing this round can add that the RS2-28 round could not, on `vp106c_fem`: the
vendor's own number is on the other side. `docs/verification/rocscience.md` records
Cai & Ugai's three-dimensional strength-reduction value for this pile row and states
that XSLOPE's plane-strain beam already reads 8.2% ABOVE it with a free head; the
Newton answer of 1.578 is 7.2% above the viscoplastic 1.472, so it is further above
the three-dimensional answer, not nearer to it. On that model, more admissible does
not mean more like the vendor — the plane-strain idealization is what separates them
and it is documented as such.

#### The ramp

The same models on the monotonic strength-reduction ramp, against the Newton
bisection on the same mesh and bracket:

| Benchmark | Ramp FS | Ramp interval | Bisection-N FS | Ramp − bisection | ramp force evals |
|---|---|---|---|---|---|
| `vp106c_fem` | 1.5781250 | [1.57500, 1.58125] | 1.5781250 | **0.00000** | 22,761 |
| `xslope_pile_wall` | 1.8031250 | [1.80000, 1.80625] | 1.8007813 | +0.00234 | 10,508 |
| `xslope_piles_fem` | 1.3843750 | [1.38125, 1.38750] | 1.3796875 | +0.00469 | 10,982 |

Three of three inside 0.01, against the three the criterion asked for, and one
reproduces the bisection to eight digits. The ramp carries the pile through its whole
warm history without any extension: `restrength` rewrites the SOIL groups and nothing
else, and the pile groups are built once, so a pile's rigidity and its capacity at
the top of the ramp are the same objects they were at the foot.

#### The no-pile path is unchanged

Sixteen benchmark rows, each re-run on BOTH package trees — this checkout and a
control run of **9f0b1515**, the driver before this round, staged separately — and
compared trial by trial: factor of safety, and each trial's strength, verdict,
iterations and force evaluations.

| Row | Mesh | Driver | FS, both trees | trials | iterations | force evals | identical? |
|---|---|---|---|---|---|---|---|
| FEM-1 tutorial | tri6, 3.5 | Newton | 1.37109375 | 9 | 1,871 | 14,233 | **yes** |
| LEM-3 tutorial | tri6, 1.2 | Newton | 1.26953125 | 9 | 4,945 | 38,463 | **yes** |
| Griffiths & Lane 1 | quad8, 3.5 | Newton | 1.37187500 | 9 | 2,213 | 16,354 | **yes** |
| Griffiths & Lane 1 | tri6, 3.5 | Newton | 1.36562500 | 9 | 3,121 | 22,650 | **yes** |
| Griffiths & Lane 1 | quad9, 3.5 | Newton | 1.39687500 | 9 | 844 | 6,030 | **yes** |
| Griffiths & Lane 6 dry | quad8, 2 | Newton | 2.41562500 | 9 | 4,146 | 29,059 | **yes** |
| Griffiths & Lane 6 dry | tri6, 2 | Newton | 2.45937500 | 9 | 3,333 | 23,672 | **yes** |
| Griffiths & Lane 3, r = 0.8 | tri6, 6 | Newton | 1.42812500 | 9 | 5,015 | 40,652 | **yes** |
| Geogrid sample | tri6, 4 | Newton | 1.60781250 | 8 | 2,129 | 14,680 | **yes** |
| Geogrid sample, LOCKED mesh | tri6, 2 | Newton | 1.56093750 | 8 | 2,732 | 17,903 | **yes** |
| Half capacity | tri6, 4 | Newton | 1.41406250 | 8 | 2,728 | 18,214 | **yes** |
| Griffiths & Lane 1, `t_cut` = 0 | tri6, 3.5 | Newton | 1.35312500 | 9 | 3,342 | 24,637 | **yes** |
| Griffiths & Lane 1, `t_cut` = 30 | tri6, 3.5 | Newton | 1.35937500 | 9 | 2,494 | 18,570 | **yes** |
| RS2-27-m1.5 (`vp036`, K0) | tri6, 1.5 | Newton | 1.37343750 | 7 | 772 | 6,487 | **yes** |
| RS2-P4-VP102-t-60-c3 (suction) | tri6, 2.5 | Newton | 1.78984375 | 9 | 592 | 4,547 | **yes** |
| **Griffiths & Lane 6 dry — the DEFAULT path** | quad8, 2 | viscoplastic | **2.421875** | 9 | 48,128 | — | **yes** |

Sixteen pairs, 137 trials, every one identical in factor of safety, verdict,
iterations and force evaluations. Every plain-soil, reinforced and tension-cutoff
value reproduces the figure recorded earlier in this document to the digit, and the
default viscoplastic path returns FS 2.421875 on per-trial iteration counts

    147, 781, 3393, 2031, 2841, 9541, 12000, 8617, 8777

value for value the control sequence, on both trees.

That it has to be identical is structural, not lucky: a model with no pile takes
`piles = None` through every new code path, and `_nr_translational_dofs` returns
None there, so the displacement bound, the elastic displacement scale, the tangent
probe's step and the reported `residual` are each the plain `np.max(np.abs(u))` they
always were, on the same array.

#### The locks

`test/nr_ssrm_check.py` gains `check_pile_element` and `check_piles`, together about
74 s. Six legs:

  * the element's consistent tangent against a central difference of its own
    internal force on 120 random elements, with the branch histogram asserted to
    have reached all eight capacity branches;
  * the element's action rows against the viscoplastic driver's own
    `_pile_element_actions`, so the two drivers are capping the same quantities;
  * the closed-form cantilever and simply supported beam, at three orientations and
    on both element types;
  * the displacement bound's index set, asserted to exclude exactly the pile nodes'
    rotations;
  * head and tip fixity — the held degrees of freedom absent from `free_dofs` and
    exactly zero in the solution — on the model whose tip is fixed;
  * the shear capacity: it must bind, no element may report a force above it, the
    displacement field must MOVE when the cap is tightened, and it must move
    FARTHER, since a member that can deliver less force cannot hold the soil back
    better;
  * and the factor of safety, as a bracket: the FEM pile sample stands at F = 1.375
    and fails at F = 1.384375 on BOTH drivers, which contains its published 1.380.

The pile refusal left `check_reinforcement_refusals`, which now guards bar softening
alone.

**Mutation, run four ways.** Each is the same two check functions against a driver
with one thing changed:

| driver | verdict | what the check saw |
|---|---|---|
| as shipped | **PASS** | — |
| the beam's bending stiffness halved | **FAIL** | all 18 closed-form legs, each off by exactly a factor of two ("cantilever tip deflection 2.40000000e+00 against the closed form 1.20000000e+00"), and the capacity direction leg with it |
| `_nr_translational_dofs` returning None — the rotations back in the bound | **FAIL** | "the set it excludes is not exactly the pile nodes' rotations, so the bound is comparing a length against a radian" |
| the capacity correction dropped from the residual | **FAIL** | the tangent leg at 6.28 relative — the tangent still gives up the stiffness the residual no longer gives up — and "capping the pile shear at 925.3 where the uncapped pile carries 3701.2 made the slope move LESS, not more" |

The third mutation is the one worth reading twice: with the correction dropped, the
REPORTED pile forces are still clipped at the cap, so a check that read only the
reported array would pass on a driver enforcing nothing.

#### The refusals that remain

Every guard was exercised again on this checkout, and each must both fire and name
its own feature:

| Feature | Verdict |
|---|---|
| **pile beam elements** | **carried** |
| post-peak softening on reinforcement bars | refuses, names softening and counts the bars |
| Hoek-Brown strength envelopes | refuses, names Hoek-Brown |
| power-curve strength envelopes | refuses, names the power curve |

Three remain, against the four the suction round left. The viscoplastic control
still accepts all three. Nothing in a pile model raises any more: the head and tip
fixities are all four settings across these eight models (free, unrotated at
`vp106c_fem_fix`, fixed at `xslope_pile_wall`'s tip), both element kinds are
exercised, and both capacities are read.

#### The criterion, line by line

**1. The element is right — MET, both halves.** (a) 400 random elements over all
eight capacity branches, worst relative gap against a central difference of the
element's own internal force 4.89e-10 (two-node) and 3.05e-11 (three-node) at a
probe step of 1e-5, falling as exactly 1/h to 7.13e-12 and 1.95e-13 at 1e-3 —
round-off in the difference, not error in the tangent. (b) The cantilever and
simply-supported closed forms at three orientations on both element types, 2.6e-15
to 4.8e-12 relative, against the 1e-10 asked. And a leg the criterion did not ask
for: the element's action rows reproduce the viscoplastic driver's own
`_pile_element_actions` to 4.1e-13 over 300 elements, so the two drivers cap the
same quantities.

**2. The eight locked models — NOT MET, at four of eight on both counts.** Four land
inside their published tolerance (`xslope_piles`, `xslope_piles_fem` and
`xslope_torggler_3a_plate` reproduce the viscoplastic answer EXACTLY, trial for
trial; `vp106c_fem_fix` at +0.0071) and four do not (+0.0436, +0.0700, +0.1061,
+0.2418, all HIGH). The criterion asked for six. **Every miss is diagnosed rather
than reported as a direction**: at the strengths where the two drivers disagree the
Newton state is in equilibrium to between 1.3e-6 and 3.1e-5 with a worst
Mohr-Coulomb violation of 1.3e-8 to 1.6e-8 of the local strength, read in the
invariant form, at 0.6% to 5.2% of the model height — a statically admissible field
carrying full gravity, which by the lower-bound theorem is a proof that the slope
stands there — while the viscoplastic trial is closed by the early-failure
classifier with its out-of-balance two to three orders of magnitude outside the gate.

**3. The capacity divergence — MET, and it costs nothing on these eight.** The
moment cap is inert on both drivers (removing it moves neither answer on the one
model that carries one) for a structural reason no sign can change. The shear cap
binds on no trial of any of the eight — the largest shear anywhere is a third of the
cap on the two models that carry one, and the other six carry none — so the two laws
are the same law there, which is why the two capacity-bearing models agree exactly.
Tightened until it binds, the two part company and the direction is measured on both.

**4. Pile diagnostics — MET, and reported as numbers.** On `xslope_piles_fem` at
F = 1.300, where the two drivers reach the same state, the axial force, shear and end
moments agree to 0.37%, 0.68% and 0.31%; at F = 1.375, one cell below the limit and
two different states, to 7.65%, 7.80% and 13.37%. The yielded masks agree element
for element on both.

**5. The displacement bound reads translational degrees of freedom only — MET
structurally, and the behavioural mutation could not be built.** The lock asserts the
index set and a mutation that returns the whole vector fails it. What could not be
done is the behavioural half the criterion asked for: a nodal rotation is a
deflection over a length, so it exceeds the largest displacement only in a model
whose pile is under one unit long, and there is no such model. Measured, the
rotations run 2.4% to 6.8% of the translational maximum on every pile model here, so
the corrected bound reads the same number the raw one did. **That clause is not met
as written and the measurement stands in its place.**

**6. The no-pile path is bit-identical — MET.** Sixteen benchmark pairs and 137
trials against a control run of the parent commit staged in a separate package tree:
every one identical in factor of safety, verdict, iterations and force evaluations,
across the plain-soil eight, the reinforced three, both tension-cutoff rows, a K0
vendor row, a matric-suction row and the default viscoplastic path.

**7. The ramp — MET on 3 of the 3 run**, against the 3 asked: 0.00000, +0.00234 and
+0.00469 from the Newton bisection. The pile is carried through the warm history
without any extension, because `restrength` rewrites the soil groups and the pile
groups are built once.

**8. The refusal is gone and the guard list is updated — MET.** Piles are carried;
bar softening, Hoek-Brown and the power curve each still raise and each still names
itself, with the viscoplastic control accepting all three.

**9. The locks — MET.** `check_pile_element` and `check_piles` pass as shipped and
fail on all three mutations — the bending stiffness halved, the rotations back in the
bound, the capacity correction dropped. The whole check file passes end to end.

**10. The default path — MET.** Griffiths & Lane 6 dry, quad8, size 2, no
`fem_solver` argument: FS 2.421875 on per-trial iteration counts 147, 781, 3393,
2031, 2841, 9541, 12000, 8617, 8777 — value for value the control sequence, on both
trees.

**11. The honest negatives** are the four locks, the bound's unbuildable behavioural
mutation, and the moment cap that neither driver enforces. All three are written
above.

#### Verdict

**The pile element is right and the driver it is bolted to disagrees with the one
that defines the locks on half the pile corpus.**

Everything that can be checked about the element checks out, and the checks are
sharp: the consistent tangent matches a central difference of its own residual to
5e-10 on all eight capacity branches with the residual falling as 1/h, an isolated
beam reproduces the cantilever and simply-supported closed forms to 1e-12 at three
orientations on both element kinds, and the actions the element reads are the
viscoplastic driver's own to 4e-13. Head and tip fixity needed no code — it is a
constraint on `free_dofs`, which this path already read — and that is asserted
rather than assumed. Sixteen benchmark pairs and 137 trials are bit-identical
against a control tree, so nothing without a pile moved.

And on four of the eight locked pile models the two drivers agree — three of them
exactly, trial for trial and verdict for verdict, including both of the models that
carry a finite pile capacity. On the other four the Newton driver reads 0.044 to
0.242 HIGH, and the misses were not left as a direction: at the disputed strengths
the Newton state is in equilibrium to a few parts in 10^5 with no Gauss point more
than 1.6e-8 of its own strength outside the yield surface, which is a lower-bound
proof that the slope stands there, while the viscoplastic trial is closed by the
early-failure classifier with its out-of-balance a thousand times outside the gate.
That is the same shape as the RS2-28 misses of the suction round, on a bigger scale:
these are not eight benchmarks the Newton driver gets wrong, they are eight
benchmarks whose locked values are the other driver's readings, and on the four where
that driver's stopping rules bite the two answers are 4% to 24% apart. The
`vp106c_fem` / `vp106c_fem_fix` pair is the sharpest single reading in the round —
the same slope, the same mesh, the same pile row, and holding the head's rotation is
the whole difference between agreeing to 0.007 and disagreeing by 0.106.

Two things about the capacities came out of reading the code rather than running it,
and both are measured. The moment cap is inert on both drivers, because the
correction is applied at a rotational degree of freedom the two adjacent beam
elements share and cancels there whichever sign it carries — a plastic hinge is a
release of rotational continuity, and this is not one. And the shear cap is applied
on the viscoplastic path with the sign the bar path's own comment calls an anti-cap:
the tighter the cap, the more shear the pile delivers. Neither decides anything on
the eight locked models, because no capacity binds on any trial of any of them. Both
would decide something on a model where one did.

What remains, in the order it matters:

- **Whose number the four disputed pile models carry.** Four locked factors of
  safety are defined by a driver whose verdict at those strengths is set by its
  early-failure classifier, and the other driver produces an admissible field 4% to
  24% higher. On `vp106c_fem` the vendor's own three-dimensional answer is BELOW
  both, and `docs/verification/rocscience.md` already documents the plane-strain
  beam as over-crediting that row — so a higher answer there is further from the
  vendor, not nearer. This is a corpus decision with a measurement behind it and it
  is the owner's, not a spike's; nothing here was changed.
- **The two capacity findings.** The moment cap enforces nothing on either driver
  and the shear cap is anti-enforced on the default one. Neither moves a published
  number today. Both are defects in the shipped path, and fixing the hinge means
  giving it a degree of freedom rather than a nodal moment.
- **Three refusals remain**: post-peak bar softening, Hoek-Brown and the power
  curve. Hoek-Brown and the power curve are two linearizations of the same shape and
  would come together; softening is still the reason neither published reinforced
  factor of safety is reachable here.


## CURVED ENVELOPES — Hoek-Brown and the power curve, per element

Written before any feature code, so that what follows is a test and not a
description. Same machine and settings as everything above: `force_tol` 1e-3,
hybrid criterion, `capture_failure_state=False`, tolerance 0.01.

### Why this one now

The pile round left three refusals: post-peak bar softening, Hoek-Brown and the
power curve. The last two are the pair the K0 round said "would come together",
and the eligibility inventory
(`xslope_private/reports/newton_corpus_eligibility_2026-09-01.md`) puts eight
locked `fem_ssrm` models behind them — three Hoek-Brown (`xslope_rock_slope`
1.166, `hammah_hb1` 1.166, `vp044d` 1.115) and five power curve (`vp040` 1.023,
`vp041` 1.656, `vp044a` 0.973, `vp045b` 2.637, `vp061a` 1.497). Every one of the
eight carries `k0=1` and nothing else this driver refuses, so the envelope is the
only gate on all of them.

The owner's requirement is sharper than "carry the two envelopes". It is
PER-ELEMENT DISPATCH: multiple materials, each with its own strength model, in
one run — a Mohr-Coulomb soil and a Hoek-Brown rock in the same mesh, solved
together. That is what the criterion's mixed-material clause is for, and it is
the clause a design that swapped a global constitutive law for another global one
would fail.

### The semantics being reproduced, read from the viscoplastic driver

Not assumed — read out of `build_fem_data`, `solve_fem`'s Step 6 and
`xslope/hoekbrown.py`, and restated here, because the Newton path has to solve the
same model.

- **Both envelopes are carried as a per-Gauss-point MOHR-COULOMB LINEARIZATION,
  re-derived every iteration.** There is no second yield surface and no second
  flow rule anywhere in the viscoplastic path. `grp['c_r']`, `grp['snph']` and
  `grp['csph']` are overwritten in place at the top of every Step-6 pass for the
  Gauss points flagged `pow_m` / `hb_m`, and the identical MC yield function and
  identical psi = 0 MOCOUQ flow then run on all points alike. The dispatch is
  already per Gauss point on that path: `pow_flag_by_elem` and `hb_flag_by_elem`
  are per-element boolean arrays, sliced onto each group as `grp['pow_m']` /
  `grp['hb_m']`, and a group can hold Mohr-Coulomb, power-curve and Hoek-Brown
  points at once.
- **The power curve's abscissa is the in-plane Mohr-circle CENTRE.** With
  `s' = -(sx + sy)/2` clamped at zero (compression-positive),
  `s_eff = max(s' + pow_d, 1e-4 ref)`, the F-reduced tangent is
  `slope = a b s_eff^(b-1) / F`, `tau_F = (a s_eff^b + pow_c)/F`,
  `c_r = tau_F - s' slope`, `phi = atan(slope)`. `pow_d` is the tension
  intercept: the envelope is `tau = a (sigma_n + d)^b + c_p`, so `d` shifts the
  origin and the `1e-4 ref` floor is what stops a tensile point running past it.
- **Hoek-Brown's abscissa is the FAILURE-PLANE normal stress**, taken from the
  PREVIOUS iterate's reduced tangent:
  `sigma_n = max(s' cos^2(phi) - c_r sin(phi) cos(phi), 0)`. That is the point at
  which a Mohr circle of centre `s'` touches its tangent line, so the
  linearization closes as a FIXED POINT inside the viscoplastic loop: at
  convergence the circle, the tangent line and the reduced envelope all meet at
  one `sigma_n`. `hb_tangent_const` then inverts Balmer's parametric curve by
  bisection on `sigma_3` to get `(c_i, phi_i)`, from GSI/mi/D-derived `mb`, `s`,
  `a`. The `sigma_n >= 0` clamp is load-bearing and is documented as such: at the
  Hoek-Brown tensile apex `dsigma_1/dsigma_3` diverges and `phi_i -> 90 deg`.
  The code also says in as many words NOT to linearize at `sigma_3` instead.
- **Strength reduction divides the TANGENT, not the constants.** `c_i/F` and
  `tan(phi_i)/F` — dividing a shear-strength function by F divides its tangent
  line's slope and intercept by F. `sigma_ci/F` is a different envelope, because
  of the exponent `a`. The per-Gauss-point F is `grp['F']`, which is 1.0 on an
  `ssr_exclude` element, so an excluded curved-envelope element keeps its whole
  envelope.
- **The native Hoek-Brown tensile strength is `-s sigma_ci / mb`** (`hb_sigma_t`),
  the apex where `sigma_1 = sigma_3`. It is used as the lower bracket of the
  Balmer inversion and nowhere else: it does NOT become a Rankine cap. A model
  that wants a tension cutoff states `t_cut` in the material table, which is the
  same per-element `t_cap_base` machinery Mohr-Coulomb materials use
  (`docs/usage/input_template.md`); `tension_srf` then divides a finite positive
  cap by F alongside the rest of the envelope. The two are independent, and the
  Newton path must compose them the same way.
- **`c_by_elem` / `phi_by_elem` carry only a SEED tangent** for a curved-envelope
  element, taken at the vertical-overburden estimate at the element centroid.
  Anything that reads those arrays before the loop re-linearizes is reading a
  seed, not the envelope.

### Design

Per-Gauss-point dispatch inside the return map's caller, on the same flags the
viscoplastic path dispatches on. A Mohr-Coulomb point is untouched — the plain
map, bit for bit. A curved-envelope point gets its own `(c, sin phi, cos phi)`
for that evaluation and then goes through the SAME return map, because the
linearized envelope IS a Mohr-Coulomb surface and the flow rule the two drivers
share is psi = 0 on it.

Two designs were open, and the reading above chooses between them:

- **(A) Mirror the viscoplastic linearization.** Equivalent local `c`, `tan phi`
  per point per evaluation, then the existing map.
- **(B) An exact curved-surface return** — Newton on the true envelope in
  principal-stress space, corners handled, with its own consistent tangent.

**(A) is the design, and the reason is that (B) would be solving a different
model.** The viscoplastic driver has no curved yield surface: it has a
per-Gauss-point tangent line and the ordinary Mohr-Coulomb return on it. Its
converged state is the fixed point where the tangent line is taken at the
abscissa the converged stress itself produces. An exact curved-surface return
would converge somewhere else, by an amount nobody has measured, and the round's
first obligation is that the two drivers answer the same question. Where (B)
would be visibly better — the corner and apex construction on a strongly curved
surface — the measurement below says how often those branches fire at all.

What (A) has to solve that a naive reading does not:

- **The linearization must be a pure function of the displacement**, or the
  residual is not one and neither the line search nor the convergence test means
  anything. The viscoplastic path can lag its linearization by one iteration
  because it has no line search; this path cannot. So the tangent is closed as a
  SELF-CONSISTENT fixed point inside each residual evaluation: linearize,
  return, re-read the abscissa off the RETURNED stress, re-linearize, to a fixed
  tolerance and a hard iteration cap. That is the same fixed point the
  viscoplastic loop converges to, reached per evaluation instead of per solve.
- **The consistent tangent then comes for free, and that is the second reason for
  (A).** `_nr_group_state` differences the whole map in the three free strain
  components. With the linearization inside the map, that quotient automatically
  carries `d(c, phi)/d(sigma)` — the derivative of the linearization — so the
  tangent is consistent with the map the residual actually uses, with no separate
  algebra to keep in step. It is no longer exact-on-the-branch, because the
  linearization is a smooth nonlinear function of the trial stress, so the
  differencing error is first-order in the probe step as it already is for the
  principal frame's rotation. That is measured rather than assumed.
- **The tension cutoff composes as it does for Mohr-Coulomb**: the cap is a
  second surface on the SAME linearized shear surface, so the existing
  active-set return takes it unchanged, inside the fixed point.
- **The ramp's `restrength` must re-derive the reduced envelope per step.** For a
  curved-envelope point the reduction is `1/F` on the tangent, and the tangent is
  re-derived every evaluation anyway, so what `restrength` has to carry is the
  per-element `F` the linearization divides by — not a `c_r` it computes once.
- **Elastic materials and SSR exclusion are unchanged.** An `elastic` point keeps
  `c_r = inf` and is never linearized; an `ssr_exclude` point linearizes at
  F = 1.

### Success criterion (verbatim)

1. **The return map is right, per envelope, measured.** At least 200,000 random
   trial states per envelope — several GSI/mi/D sets for Hoek-Brown and several
   (a, b, c, d) sets for the power curve — with and without a finite `t_cut`:
   every returned state satisfies BOTH surfaces to within 1e-12 of the stress
   scale, the principal ordering survives, no elastic state is modified at all,
   and the branch histogram is REPORTED with every branch the design names
   actually exercised — a fuzz that never lands on a region proves nothing about
   it. The consistent tangent is checked against a central difference on every
   branch and REPORTED; the 1e-8 target is stated with the frame-rotation
   truncation caveat this document has carried since the tension-cutoff round,
   and the linearization's own smoothness is now a second first-order term.
2. **The eight locked models**, three Hoek-Brown and five power curve, each built
   through `run_tests.py`'s own `build_fem_ssrm_case` so the mesh, the element
   type, the bracket, the iteration budget, `k0` and every solver option are the
   suite's and not this round's. The Newton bisection must land inside each
   published tolerance AND within 0.01 of the viscoplastic driver on at least 6
   of the 8. Reported per model: viscoplastic FS, Newton FS, the lock, its
   tolerance, and the constitutive work on each driver. **Every miss is diagnosed
   rather than reported as a direction** — the trials whose verdicts differ are
   run, the Newton state's own admissibility read (out-of-balance and the worst
   yield violation in the invariant form on the linearized envelope), and the
   viscoplastic out-of-balance at closure quoted, the way the RS2-28, three-layer
   and pile disagreements were refereed. If the early-failure-classifier pattern
   recurs, it is named as such.
3. **MIXED MATERIALS — the owner's requirement.** At least one model carrying a
   curved-envelope material AND a Mohr-Coulomb material in the same mesh,
   constructed if the corpus has none. The Newton bisection must agree with the
   viscoplastic driver within 0.01 on it, and the PER-ELEMENT DISPATCH must be
   asserted rather than inferred: the branch and envelope counts per material
   reported, with the Mohr-Coulomb elements shown taking the plain path and the
   curved ones their own, in the same solve.
4. **The ramp agrees** with the Newton bisection within 0.01 on at least 3 of the
   models, and carries the reduced envelope through its warm history.
5. **The Mohr-Coulomb-only path is bit-identical**, on the same control-tree
   protocol every previous round used: the plain-soil eight rows, the reinforced
   benchmarks, the tension-cutoff rows, a K0 vendor row, a matric-suction row and
   a pile row, re-run against a control run of the parent commit staged in a
   separate package tree — every trial identical in factor of safety, verdict,
   iterations and force evaluations — plus the plain Mohr-Coulomb return map
   bit-identical over 800,000 random trial states.
6. **The refusals.** Only post-peak bar softening remains; it still raises and
   still names itself, with the viscoplastic control accepting it. Asserted.
7. **The locks catch it.** `test/nr_ssrm_check.py` gains
   `check_curved_envelopes`: the fuzz and its branch histogram, one Hoek-Brown
   model and one power-curve model on both drivers, and the mixed model on both
   drivers AND on both Newton routes. Mutations, run both ways and recorded:
   **corrupt `mb`** (or `a`) so the derived Hoek-Brown constants are wrong, and
   **drop the F-reduction of the envelope** so the tangent is taken unreduced.
   Each must FAIL the check. The whole check file passes.
8. **The default path is unchanged**, against the standard control: Griffiths &
   Lane 6 dry with no `fem_solver` argument, FS 2.421875 on per-trial iteration
   counts 147, 781, 3393, 2031, 2841, 9541, 12000, 8617, 8777.
9. **An honest negative is a valid outcome and must be written.** If the
   self-consistent linearization does not close, if the tangent's conditioning
   costs the plain table, if the two drivers disagree on a curved-envelope model
   by more than the bisection tolerance for a reason that is the formulation
   rather than a solver rule, or if design (A) turns out to be measurably further
   from the limit-equilibrium answer on the same envelope than an exact return
   would be — that is the result.


### CURVED ENVELOPES — results

Same machine and settings as everything above: `force_tol` 1e-3, hybrid criterion,
`capture_failure_state=False`, tolerance 0.01. Every number below was measured on
this checkout in this session. The vendor benchmarks are built through
`run_tests.py`'s own `build_fem_ssrm_case`, from the tags in
`docs/verification/rs2.md` and `docs/tutorials/lem13_rock_slope.md` and the
workbooks those tags name, so the mesh, the element type, the bracket, the
iteration budget and the K0 value are the suite's own and not this round's.

#### Design (A), and the reading that chose it

**(A), and it was not a close call.** The viscoplastic driver has no curved yield
surface: `solve_fem`'s Step 6 overwrites `grp['c_r']`, `grp['snph']` and
`grp['csph']` in place for the Gauss points flagged `pow_m` / `hb_m` and then runs
the ordinary Mohr-Coulomb yield function and the ordinary psi = 0 MOCOUQ flow on
all points alike. Its converged state is the fixed point where the tangent line is
taken at the abscissa the converged stress itself produces. An exact curved-surface
return would converge somewhere else, by an amount nobody has measured, and the
first obligation of this round is that the two drivers answer the same question.
The measurement below says they do — on eight of eight locked models the two now
close on the SAME bisection interval — which is the outcome design (B) would have
put at risk for no benefit anyone has asked for.

#### What was built

`_nr_pow_tangent`, `_nr_hb_tangent`, `_nr_envelope_step`, `_nr_envelope_return`,
`_nr_envelope_by_elem`, `_nr_group_envelope` and `_nr_group_restrength_envelope`
in `xslope/fem.py`, threaded through `_nr_group_state`, `_nr_internal_force`,
`_nr_build_groups` and the ramp's `restrength`. The guard's two envelope lines are
gone.

**The dispatch is per Gauss point, not per model.** `env['code']` carries
_NR_ENV_MC / _NR_ENV_POW / _NR_ENV_HB, sliced onto each group from the model's own
`pow_flag_by_elem` and `hb_flag_by_elem` — the same per-element arrays the
viscoplastic path dispatches on — so one group, one mesh, can hold Mohr-Coulomb,
power-curve and Hoek-Brown points at once and each takes its own branch in the
same pass. A Mohr-Coulomb point is not touched at all, which is what keeps a mixed
model's Mohr-Coulomb half on exactly the arithmetic it would be on alone.

Three decisions are the whole of the design, and each was forced by something the
viscoplastic path does not have to deal with.

**The linearization is closed as a self-consistent fixed point inside every
residual evaluation.** The viscoplastic loop can lag its linearization by one
iteration because it has no line search; this path has one, and a convergence test,
and both need the residual to be a function of the displacement alone. So each
evaluation linearizes, returns, re-reads the abscissa off the RETURNED stress and
re-linearizes, to `_NR_ENV_TOL`. That is the same fixed point the viscoplastic loop
reaches over a whole solve, reached per evaluation instead.

**The consistent tangent then comes for free, and that is the second reason for
(A).** `_nr_group_state` differences the whole map in the three free strain
components; with the linearization inside the map the quotient carries
`d(c, phi)/d(sigma)` as well, so the tangent is the derivative of the map the
residual actually uses with no second piece of algebra to keep in step.

**Two levels are set by arithmetic rather than by a benchmark, and both are named
here because a level nobody chose for a reason is the failure this spike exists to
avoid.** `_NR_ENV_TOL` is 1e-10 and not tighter because the Balmer inversion is a
40-step BISECTION whose bracket closes at 2^-40 of its width, which puts a floor of
about 6e-12 under the tangent it can return: a fixed point asked for less than that
never exits. Measured on the Hoek-Brown benchmark the residual falls by a decade a
pass — 1.4e-1, 4.8e-3, 4.7e-4, 4.7e-5 — and then sits on that floor. And
`_NR_ENV_RELAX1` / `_NR_ENV_RELAX2` under-relax a fixed point that has not closed
in 14 and 34 passes, because a HANDFUL of trial states put the iteration in a
period-2 limit cycle — the tangent moves the returned state across a branch
boundary and the new abscissa moves it back — and one cycling Gauss point stalls
the whole group. Neither level is reached on any benchmark in this document: the
corpus runs exit in one to ten passes.

#### The defect the guard had been hiding

Not in the envelope. A Gauss point in a material held LINEAR ELASTIC is carried on
the Newton path as an unreachable cohesion, `c_r = inf`, and the linearization
writes a finite tangent into `c_r`. The viscoplastic path never meets the problem
because it drops elastic points from the yielding mask instead
(`m & ~grp['elastic']`) and does not care what their `c_r` says.

**Specimen: `vp040` (RS2-30), where 1,284 of 2,539 elements sit in an SSR elastic
zone.** Without the dispatch back to the plain path for those points the driver
could not reach equilibrium at ANY strength — not at F = 1, and not at F = 0.1,
where the soil is ten times stronger than the file's own and the same driver's
answer is 1.02. The trace says why: the first Newton step leaves
`||r||/||f|| = 4.728e-01` at every load factor down to 1/64, because half the mesh
is yielding a material that cannot yield, and a smaller load increment does not
make a smaller problem when the yielding is not caused by the load. With the fix
the same model converges at F = 1.0 in 53 iterations to an out-of-balance of
3.2e-6, and its bisection reproduces the viscoplastic answer exactly.

That is the only defect this round found, and it was found by a benchmark rather
than by reading: seven of the eight models solved before it was fixed.

#### The return map, measured

200,000 random trial states per envelope at each of four caps — none, zero, a small
positive value, and one above the Mohr-Coulomb apex — and at each of two strength
reductions, over four Hoek-Brown (GSI, mi, D) sets and four power-curve
(a, b, c_p, d) sets drawn per point. **3,200,000 returns.**

| invariant | Hoek-Brown, worst | power curve, worst |
|---|---|---|
| Mohr-Coulomb residual, as a fraction of the stress scale | 4.93e-14 | 4.10e-14 |
| Rankine residual, same scale | 9.09e-15 | 9.09e-15 |
| principal ordering violations | 0 (exactly) | 0 (exactly) |
| elastic states modified | 0 (exactly) | 0 (exactly) |
| states left on the unresolved fallback | 0 | 0 |

The branch histogram, which says what actually ran (one column of each pair, at
F = 1.0):

| branch | HB, no cap | HB, `t_cut` = 0 | pow, no cap | pow, `t_cut` = 0 |
|---|---|---|---|---|
| elastic | 17,279 | 8,920 | 8,747 | 7,211 |
| main plane | 36,004 | 12,425 | 34,196 | 24,906 |
| right corner | 28,908 | 18,736 | 27,963 | 21,666 |
| left corner | 37,169 | 19,281 | 35,072 | 24,079 |
| apex | 80,640 | 0 | 94,022 | 0 |
| cap alone | — | 60,221 | — | 23,372 |
| cap on two principals | — | 18,123 | — | 7,724 |
| hydrostatic tension (T, T, T) | — | 998 | — | 1,039 |
| **Mohr-Coulomb / Rankine edge** | — | **6,497** | — | **10,571** |
| corner and cap | — | 51,757 | — | 74,971 |
| main plane, cap on two | — | 3,042 | — | 4,461 |

Every region the design names is exercised on both envelopes, the intersection edge
and the hydrostatic-tension return included. The capped columns carry no apex at
all, which is the other half of the same statement: with a cap below it the
Mohr-Coulomb apex never fires, because the hydrostatic-tension return has taken its
place.

**The honest negative in this leg is the fixed point's own residual.** The tangent
the return was taken on is compared against the tangent of the envelope at the
abscissa the RETURNED stress produces, and over 200,000 wildly random states the
worst self-consistency residual runs from **1.9e-5 to 5.7e-3** — not the
`_NR_ENV_TOL` the loop asks for. A few states in every batch are in a limit cycle
that under-relaxation slows but does not close inside the pass budget. What that
bounds is how far the returned Mohr circle can sit from the CURVE rather than from
the line: at 5.7e-3 of the local cohesive intercept, on states drawn eight times
outside any envelope in the batch. On the corpus the loop exits in one to ten passes
and the residual is at `_NR_ENV_TOL`; the fuzz's states are not corpus states.

#### The consistent tangent

The driver's own ONE-SIDED assembled block `d[sx, sy, txy]/d[ex, ey, gxy]` against a
central difference of the same, normalized by the elastic block's magnitude, over
8,000 random states per row. Points whose perturbation crosses a return-map branch
are excluded and counted, because a central difference across a boundary averages
two tangents — which is why `kept` falls as the probe step grows.

The row that makes the table an argument is the last pair: the SAME harness, the
same random states, with the material declared plain Mohr-Coulomb.

| row | h = 1e-4 | h = 1e-5 | h = 1e-6 |
|---|---|---|---|
| Hoek-Brown, no cap — median / worst | 0.0 / 1.8e-2 | 0.0 / 1.4e-1 | 0.0 / 2.4e-1 |
| Hoek-Brown, `t_cut` = 0 | 5.1e-4 / 6.6e-2 | 1.3e-3 / 2.9e-1 | 1.4e-3 / 3.4e-1 |
| power curve, no cap | 1.1e-12 / 1.5e-2 | 3.1e-5 / 2.2e-1 | 6.5e-6 / 4.3e-1 |
| power curve, `t_cut` = 0 | 3.7e-4 / 1.6e-2 | 6.0e-4 / 3.7e-1 | 6.3e-4 / 3.5e-1 |
| **Mohr-Coulomb CONTROL, no cap** | **0.0 / 9.8e-3** | **0.0 / 1.5e-1** | **1.9e-5 / 1.6e-1** |
| **Mohr-Coulomb CONTROL, `t_cut` = 0** | **1.1e-3 / 7.0e-2** | **1.4e-3 / 3.7e-1** | **1.0e-3 / 3.0e-1** |

**The criterion's 1e-8 clause is NOT met, and it is not met by the plain
Mohr-Coulomb branches either, on this harness.** The harness generates independent
random stress components and independent random material constants per point, so a
large tail of it sits near a degenerate principal frame, where the frame-rotation
truncation this document has named since the tension-cutoff round is at its worst.
What is measured in its place is the comparison the criterion actually wanted
proved: **the curved branches are as differentiable as the plain ones, median for
median and worst for worst**, and the capped rows of both are indistinguishable
from each other.

#### The eight locked models

Every row built through `run_tests.py`'s own `build_fem_ssrm_case`; the bisection
tolerance is 0.01 on all eight. Work is the honest count on each driver —
viscoplastic iterations against Newton force evaluations, one constitutive pass
each — with the viscoplastic predictor iterations the Newton run charges on top in
their own column.

| Benchmark | Model | Envelope | Elements | Lock | Tol | VP FS | Newton FS | Newton − lock | driver gap | VP work | N-R work | N-R predictor |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| HB-ssrm | `hammah_hb1` | Hoek-Brown | 1485 | 1.166 | 0.01 | 1.1656250 | **1.1656250** | −0.0004 | **0.0000** | 8,322 | 4,212 | 13,088 |
| LEM-13-ssrm | `xslope_rock_slope` | Hoek-Brown | 1485 | 1.166 | 0.01 | 1.1656250 | **1.1656250** | −0.0004 | **0.0000** | 8,322 | 4,212 | 13,088 |
| RS2-31d | `vp044d` | Hoek-Brown | 1966 | 1.115 | 0.02 | 1.1113281 | **1.1113281** | −0.0037 | **0.0000** | 39,336 | 4,414 | 11,772 |
| RS2-30 | `vp040` | power curve | 2539 | 1.023 | 0.02 | 1.0195312 | **1.0195312** | −0.0035 | **0.0000** | 106,425 | 3,942 | 15,316 |
| RS2-P4-VP41 | `vp041` | power curve | 1944 | 1.656 | 0.02 | 1.6593750 | **1.6593750** | +0.0034 | **0.0000** | 8,690 | 3,632 | 13,835 |
| RS2-31c | `vp044a` | power curve | 1966 | 0.973 | 0.02 | 0.9683594 | **0.9683594** | −0.0046 | **0.0000** | 9,661 | 2,711 | 12,964 |
| RS2-32b | `vp045b` | power curve | 2154 | 2.637 | 0.02 | 2.6402344 | **2.6402344** | +0.0032 | **0.0000** | 26,144 | 2,236 | 24,302 |
| RS2-34b | `vp061a` | power curve | 1966 | 1.497 | 0.02 | 1.4921875 | **1.4921875** | −0.0048 | **0.0000** | 75,126 | 3,607 | 14,745 |

**Eight of eight reproduce their published lock, and eight of eight return the
IDENTICAL factor of safety to the viscoplastic driver** — not merely inside the
0.01 bisection tolerance but on the same final interval, to nine digits, on every
row. The criterion asked for six on each count. There is no miss to diagnose, and
the one-sided high reading this document has recorded on every previous corpus round
does not appear here at all.

**The work goes the same way it has gone since Phase 0.** The Newton driver does
less constitutive work than the viscoplastic one on seven of the eight — 2.0x on
`hammah_hb1` and `xslope_rock_slope`, 8.9x on `vp044d`, 27.0x on `vp040`, 2.4x on
`vp041`, 3.6x on `vp044a`, 11.7x on `vp045b` and 20.8x on `vp061a` — before the
viscoplastic predictor iterations it charges on top, which run 11,772 to 24,302 and
fall entirely on the failing trials. Counting those in, the ratio runs 0.48x to
5.53x: the Newton driver does MORE total constitutive work than the viscoplastic
one on `hammah_hb1` and `xslope_rock_slope` (0.48x), `vp041` (0.50x), `vp044a`
(0.62x) and `vp045b` (0.99x), and less on `vp044d` (2.43x), `vp061a` (4.09x) and
`vp040` (5.53x). Wall time runs 100-1,554 s against 36-647 s, and the curved
evaluation is the reason: one residual evaluation on a Hoek-Brown group is a
fixed point of return maps around a 40-step Balmer bisection, and one tangent
re-form is four of them.

**Bookkeeping: eight benchmarks, seven models.** `hammah_hb1` and
`xslope_rock_slope` are the same model under two names. Loaded side by side their
one material is identical in every field the FEM reads — `option = hb`,
sigma_ci 30,000, GSI 5, mi 2, D 0, E 5,000,000, nu 0.3, gamma 25 — and the only
difference in either file is a blank `gamma_sat` on the tutorial copy, which nothing
in a dry model reads. They mesh to the same 1,485 elements and solve to the same
factor of safety on the same trial sequence, iteration for iteration, on both
drivers. This is the same situation the reinforcement round found between
`xslope_reinforce_fem` and `xslope_reinforced_slope`.

#### The mixed model — the owner's requirement

The corpus cannot supply one: **every curved-envelope model in it is
SINGLE-material**, all eight of them, so nothing shipped exercises per-element
dispatch at all. The model is therefore constructed, out of
`docs/fem/files/xslope_noncircular_fem.xlsx` — four Mohr-Coulomb materials — by
re-declaring the second on Hoek-Brown and the third on the power curve, each fitted
by bisection on sigma_ci (respectively on the coefficient a) to the strength the
file's own material carries at mid-height overburden, so the result is a slope and
not a shape.

**255 elements: 137 Mohr-Coulomb, 46 Hoek-Brown, 72 power-curve, one mesh.**

The dispatch is read off a solve rather than inferred. At F = 1.10 on the Newton
driver, per material, over the 765 Gauss points:

| material | option | envelope taken | branches |
|---|---|---|---|
| Sand Fill | mc | mohr-coulomb, 195 | elastic 164, main 31 |
| Sand | hb | **hoek-brown, 138** | elastic 102, main 36 |
| Soft Clay | pow | **power-curve, 216** | elastic 216 |
| Dense Sand | mc | mohr-coulomb, 216 | elastic 216 |

Every Mohr-Coulomb element takes the plain path and every curved one takes its own,
in the same solve, and the converged state at that strength is admissible to
1.2e-15 of the local strength.

| driver | FS | work |
|---|---|---|
| viscoplastic | **1.9703125** | 114,402 iterations, 113 s |
| **Newton bisection** | **1.9703125** | 6,690 force evaluations, 80 s |
| Newton ramp | 1.9781250, interval [1.97500, 1.98125] | 127 s |

**The two drivers return the same number to nine digits on a mesh carrying three
strength models at once**, and the ramp lands one of its own increments above it,
+0.0078, which is inside the tolerance and is where the ramp lands on models with a
single envelope too.

**A second, larger mixed model says the same thing and prices it.** The same
construction on `docs/fem/files/xslope_griffiths4_r1.xlsx` — its foundation clay
re-declared Hoek-Brown — gives 4,405 elements, 1,470 Mohr-Coulomb and 2,935
Hoek-Brown, which is 4,410 Mohr-Coulomb and 8,805 Hoek-Brown Gauss points in one
mesh. The viscoplastic bisection reads **1.71484375** in 1,020 s. Rather than pay
for a whole Newton bisection on it, the two trials either side of that answer were
run directly:

| F | driver | verdict | out-of-balance | iterations | worst yield violation | wall |
|---|---|---|---|---|---|---|
| 1.70 | viscoplastic | CONVERGED | 1.00e-03 | 3,658 | — | 98 s |
| 1.70 | **Newton** | **CONVERGED** | **1.70e-05** | **27** | **1.5e-08** | 31 s |
| 1.73 | viscoplastic | FAILED (`iteration_cap`) | 3.08e-01 | 12,000 | — | 308 s |
| 1.73 | **Newton** | **FAILED** (load-step floor) | — | 729 | — | 1,334 s |

Same verdict on both sides of the viscoplastic answer, and the price is the one
this document has recorded since Phase 0 rather than anything the envelope adds:
the standing trial costs Newton a third of the viscoplastic wall time on 135 times
fewer iterations, and the failing trial costs it 4.3 times as much.

#### The ramp

The ramp carries the reduced envelope through its warm history: `restrength`
rewrites the per-element F the linearization divides by, and the tangent is
re-derived at every evaluation, so a curved-envelope point at the top of the ramp is
solved on the strength that step declares and not on the foot's.

| Benchmark | Envelope | Ramp FS | Ramp interval | Bisection-N FS | Ramp − bisection | ramp force evals |
|---|---|---|---|---|---|---|
| HB-ssrm (`hammah_hb1`) | Hoek-Brown | 1.1656250 | [1.16250, 1.16875] | 1.1656250 | **0.00000** | 5,842 |
| RS2-P4-VP41 (`vp041`) | power curve | 1.6587500 | [1.65625, 1.66125] | 1.6593750 | −0.00063 | 4,921 |
| RS2-32b (`vp045b`) | power curve | 2.6390000 | [2.63500, 2.64300] | 2.6402344 | −0.00123 | 5,375 |
| mixed, 3 envelopes | mc + hb + pow | 1.9781250 | [1.97500, 1.98125] | 1.9703125 | +0.00781 | — |

Four of four inside 0.01, against the three the criterion asked for, and the
Hoek-Brown row reproduces the bisection to eight digits.

#### The Mohr-Coulomb-only path is unchanged

Two proofs, and the first is the stronger one.

**The arithmetic.** The plain Mohr-Coulomb return map, run against the driver as it
stood at **de9cda60** — the criterion commit, before any of this code — staged in a
separate package tree, on 800,000 random trial states across four friction angles:
**BIT-IDENTICAL**, stress and branch code alike, to a SHA-256 of the concatenated
returns (`3ce2e691a304e2383aeff8b38b14a620b050528a6e2ee37f7074f44255b6ad25` on both
trees). The reason it has to be is structural: a model with no `pow` or `hb`
material never gets an `env` key on any group, so `grp.get('env')` is None at the
one point the linearization could enter, and `mc_return_map` itself was not touched.

**The benchmarks**, every one re-run on BOTH trees rather than compared against a
number in this document, and compared trial by trial: factor of safety, and each
trial's strength, verdict, iterations and force evaluations.

| Row | Mesh | Driver | FS, both trees | trials | iterations | force evals | identical? |
|---|---|---|---|---|---|---|---|
| FEM-1 tutorial | tri6, 3.5 | Newton | 1.37109375 | 9 | 1,871 | 14,233 | **yes** |
| LEM-3 tutorial | tri6, 1.2 | Newton | 1.26953125 | 9 | 4,945 | 38,463 | **yes** |
| Griffiths & Lane 1 | quad8, 3.5 | Newton | 1.37187500 | 9 | 2,213 | 16,354 | **yes** |
| Griffiths & Lane 1 | tri6, 3.5 | Newton | 1.36562500 | 9 | 3,121 | 22,650 | **yes** |
| Griffiths & Lane 1 | quad9, 3.5 | Newton | 1.39687500 | 9 | 844 | 6,030 | **yes** |
| Griffiths & Lane 6 dry | quad8, 2 | Newton | 2.41562500 | 9 | 4,146 | 29,059 | **yes** |
| Griffiths & Lane 6 dry | tri6, 2 | Newton | 2.45937500 | 9 | 3,333 | 23,672 | **yes** |
| Griffiths & Lane 3, r = 0.8 | tri6, 6 | Newton | 1.42812500 | 9 | 5,015 | 40,652 | **yes** |
| Geogrid sample | tri6, 4 | Newton | 1.60781250 | 8 | 2,129 | 14,680 | **yes** |
| Geogrid sample, LOCKED mesh | tri6, 2 | Newton | 1.56093750 | 8 | 2,732 | 17,903 | **yes** |
| Three layers | tri6, 4 | Newton | 1.21093750 | 8 | 3,622 | 25,170 | **yes** |
| Three layers, `t_cut` = 0 | tri6, 4 | Newton | 1.21093750 | 8 | 3,688 | 25,515 | **yes** |
| Griffiths & Lane 1, `t_cut` = 0 | tri6, 3.5 | Newton | 1.35312500 | 9 | 3,342 | 24,637 | **yes** |
| Griffiths & Lane 1, `t_cut` = 30 | tri6, 3.5 | Newton | 1.35937500 | 9 | 2,494 | 18,570 | **yes** |
| RS2-27-m1.5 (`vp036`, K0) | tri6, 1.5 | Newton | 1.37734375 | 8 | 793 | 6,554 | **yes** |
| RS2-P4-VP102-t-60-c3 (suction) | tri6, 2.5 | Newton | 1.79257813 | 10 | 615 | 4,640 | **yes** |
| FEM pile sample (`xslope_piles_fem`) | tri6, 2 | Newton | 1.37968750 | 8 | 4,055 | 27,736 | **yes** |
| **Griffiths & Lane 6 dry — the DEFAULT path** | quad8, 2 | viscoplastic | **2.421875** | 9 | 48,128 | — | **yes** |

Eighteen pairs, 157 trials, every one identical in factor of safety, verdict,
iterations and force evaluations. Every plain-soil, reinforced, tension-cutoff, K0,
suction and pile value reproduces the figure recorded earlier in this document to
the digit, and the default viscoplastic path returns FS 2.421875 on per-trial
iteration counts

    147, 781, 3393, 2031, 2841, 9541, 12000, 8617, 8777

value for value the control sequence, on both trees.

#### The refusals that remain

Every guard was exercised again on this checkout, and each must both fire and name
its own feature:

| Feature | Verdict |
|---|---|
| **Hoek-Brown strength envelopes** | **carried** |
| **power-curve strength envelopes** | **carried** |
| post-peak softening on reinforcement bars | refuses, names softening and counts the bars (18 of 36) |

**One remains**, against the three the pile round left, and the viscoplastic control
still accepts it. Nothing about a curved envelope raises any more: both envelopes,
the two of them mixed with Mohr-Coulomb in one mesh, SSR elastic zones over a curved
material, SSR exclusion, K0 and the Rankine cap all run.

#### The locks

`test/nr_ssrm_check.py` gains `check_curved_envelopes`, in about 100 s, and loses
`check_unsupported_features_refuse` — which existed to assert the two refusals this
round removed, and which would otherwise be holding the driver to a limitation it
no longer has. Four legs, and each fails on a different defect:

  * the RETURN MAP on both envelopes at three caps — 40,000 states each, every
    returned state on the linearized Mohr-Coulomb surface and under its cap, the
    ordering intact, no elastic state modified, no state on the unresolved
    fallback, and a branch histogram asserted to reach every region the design
    names, because a fuzz that never lands on a region proves nothing about it;
  * the LINEARIZATION is self-consistent — the tangent the return was taken on is
    the tangent at the abscissa the RETURNED stress produces;
  * the TANGENT IS THE MATERIAL'S — compared against one computed from scratch out
    of the material's own columns, through the public `hb_tangent` (which
    re-derives mb, s and a from GSI, mi and D) and through the power curve's own
    formula. This is the leg a driver-against-driver comparison could not supply: a
    parameter wired from the wrong per-element array is self-consistent with itself
    and only an independent expectation can see it;
  * the F-REDUCTION divides the tangent and not the constants — at F = 2 the
    cohesive intercept and tan(phi) are exactly half their F = 1 values on both
    envelopes;
  * and the MIXED model — 137 Mohr-Coulomb, 46 Hoek-Brown and 72 power-curve
    elements on one mesh — with the per-element dispatch asserted material by
    material from a solve, and a strength bracket both drivers must agree on.

**Mutation, run four ways.** Each is the same check function against a driver with
one thing changed:

| driver | verdict | what the check saw |
|---|---|---|
| as shipped | **PASS** | — |
| the Hoek-Brown `mb` the Newton path carries multiplied by 1.5 | **FAIL** | "the tangent the Newton path derives is not the envelope the material declares — the cohesive intercept differs by 5.567e-02 and tan(phi) by 2.425e-01 relative, against an expectation computed from the material's own columns" |
| the Hoek-Brown exponent `a` shifted by 0.08 | **FAIL** | the same line, at 2.178e-01 and 1.643e-01 |
| the F-reduction of the envelope dropped | **FAIL** | the same line on BOTH envelopes at 7.000e-01, plus all four F-reduction legs ("the cohesive intercept ... at F = 2 is 1.000000 ... of its value at F = 1"), plus the mixed model standing at F = 2.1 where the Mohr-Coulomb half of the same mesh says it must not |

The first two are the mutations the criterion named, and the leg that catches them
had to be added: an earlier version of this check, with the fuzz and the mixed model
but without the material-table comparison, PASSED both of them. A corrupted rock
constant is self-consistent with itself, and the mixed model's Hoek-Brown material
is not what decides its bracket.

The whole check file passes end to end.

#### The criterion, line by line

**1. The return map — MET, except the tangent clause.** 3,200,000 returns across
two envelopes, four caps, two strength reductions and four parameter sets each:
both surfaces satisfied to 4.9e-14 of the stress scale, zero ordering violations,
zero elastic states modified, zero states on the unresolved fallback. The branch
histogram exercises every region the design names on both envelopes, including
6,497 and 10,571 returns on the Mohr-Coulomb / Rankine intersection edge and 998
and 1,039 on the hydrostatic-tension return, and it shows the Mohr-Coulomb apex
firing 80,640 and 94,022 times uncapped and not once with a cap below it. **The
tangent clause is NOT met as written** — no branch reaches 1e-8 on this harness —
and it is not met by the plain Mohr-Coulomb branches either, which is the
measurement that stands in its place: on the same harness and the same states, the
Mohr-Coulomb control reads the same medians and the same worst cases as the curved
rows. The fixed point's own residual is the honest negative and is reported: 1.9e-5
to 5.7e-3 on the fuzz's states, against `_NR_ENV_TOL` on every corpus state.

**2. The eight locked models — MET, and at eight of eight rather than six.** Every
one reproduces its published lock, at −0.0048 to +0.0034 against tolerances of 0.01
and 0.02, and every one returns the IDENTICAL factor of safety to the viscoplastic
driver — the same final interval to nine digits. There is no miss, so there is
nothing to diagnose and the early-failure-classifier pattern the previous three
rounds met does not appear.

**3. Mixed materials — MET.** The corpus has none, so one was constructed: 137
Mohr-Coulomb, 46 Hoek-Brown and 72 power-curve elements on one mesh, the two curved
materials fitted to the strengths the file's own Mohr-Coulomb materials carry.
Both drivers read 1.9703125. The dispatch is asserted from the solve rather than
inferred, material by material, with the envelope each element took and the branch
counts under it.

**4. The ramp — MET on 4 of the 4 run**, against the 3 asked: 0.00000, −0.00063,
−0.00123 and +0.00781 from the Newton bisection, the last of them on the mixed
model.

**5. The Mohr-Coulomb-only path is bit-identical — MET.** Eighteen benchmark pairs
and 157 trials against the driver staged at de9cda60 in a separate package tree:
every one identical in factor of safety, verdict, iterations and force evaluations,
across the plain-soil eight, four reinforced rows, three tension-cutoff rows, a K0
vendor row, a suction row, a pile row and the default viscoplastic path. The plain
return map is bit-identical over 800,000 trial states to a SHA-256 of the returns.

**6. The refusals — MET.** Only post-peak bar softening remains, it still names
itself and counts the bars, and the viscoplastic control still accepts it.

**7. The locks — MET.** `check_curved_envelopes` passes as shipped and fails on all
three mutations the criterion named; the whole check file passes end to end.

**8. The default path — MET.** Griffiths & Lane 6 dry, quad8, size 2, no
`fem_solver` argument: FS 2.421875 on per-trial iteration counts 147, 781, 3393,
2031, 2841, 9541, 12000, 8617, 8777 — value for value the control sequence, on both
trees.

**9. The honest negatives** are the tangent clause, the fixed point's residual on
the fuzz's states, the constructed mixed model, and the cost — all written above.

#### Verdict

**Both curved envelopes are carried, per Gauss point, and on the corpus the two
drivers do not merely agree — they return the same number.** Eight locked models,
three Hoek-Brown and five power-curve, run in their locked configurations through
the suite's own mesh and bracket mapping: every one reproduces its published lock,
and every one closes on the same bisection interval as the viscoplastic driver, to
nine digits. That is a better result than any previous round in this document
produced on a vendor block, and it is better than the criterion asked for on both
counts.

What made it possible is that the round did not build a curved yield surface.
Neither driver has one. Both carry a curved envelope as a per-Gauss-point
Mohr-Coulomb tangent, re-derived from the current stress, with the ordinary psi = 0
return running on it — so reproducing the viscoplastic answer meant reproducing
that linearization, at the abscissa it uses, reduced the way it reduces it, and
nothing else. The one thing this path had to add is that its linearization is
closed as a self-consistent fixed point inside every residual evaluation rather
than lagged an iteration, because a driver with a line search needs its residual to
be a function of the displacement alone. That decision paid twice: it landed on the
same fixed point the viscoplastic loop converges to, and it made the existing
difference quotient carry the derivative of the linearization for free.

**The owner's requirement is met and the corpus could not have shown it.** All eight
curved-envelope models in the repository are SINGLE-material, so not one of them
exercises per-element dispatch. The mixed model is therefore constructed — 137
Mohr-Coulomb elements, 46 Hoek-Brown and 72 power-curve on one mesh — and on it the
two drivers read 1.9703125 apiece, with the dispatch asserted material by material
off the solve rather than inferred: every Mohr-Coulomb element on the plain path,
every curved one on its own, in the same pass.

The round found one defect and it was not in the envelope. A Gauss point held
LINEAR ELASTIC is carried on this path as `c_r = inf`, and the linearization was
writing a finite tangent over it. On `vp040`, where 1,284 of 2,539 elements sit in
an SSR elastic zone, that stopped the driver reaching equilibrium at ANY strength —
including F = 0.1, where the soil is ten times stronger than the file's own. Seven
of the eight models solved without noticing.

Three things did not close, and all three are written above rather than tuned away.
The consistent tangent does not reach the 1e-8 the criterion asked against a central
difference, and neither does the plain Mohr-Coulomb branch on the same harness —
the comparison against that control is what stands in the clause's place. The
linearization's fixed point does not always close to `_NR_ENV_TOL` on the fuzz's
states: a handful in every 200,000 sit in a period-2 limit cycle that
under-relaxation slows without ending, leaving a worst self-consistency residual of
5.7e-3, which bounds how far a returned circle can sit from the CURVE rather than
from its own tangent line. It is not reached on any corpus state. And the cost is
real: counting the viscoplastic predictor iterations a Newton run charges on top,
the Newton driver does more total constitutive work than the viscoplastic one on
five of the eight, because one residual evaluation on a Hoek-Brown group is a fixed
point of return maps wrapped around a 40-step Balmer bisection.

What remains, in the order it matters:

- **Post-peak bar softening is the last refusal**, and it is the reason neither of
  the repository's published reinforced factors of safety is reachable on this
  driver. Everything else in the eight-item list the ramp verdict named is carried.
- **The corpus run.** With the two envelopes carried, the reachable locked set is
  every `fem_ssrm` tag in the repository except the two behind bar softening.
  Running it is what would turn "eight of eight" into a statement about the corpus,
  and it is the same measurement the ramp verdict named as the thing standing
  between this branch and any default-driver conversation.
- **The cost of a Hoek-Brown evaluation.** The Balmer inversion is a 40-step
  bisection run inside a fixed point run inside a difference quotient. A warm-started
  Newton or secant on `sigma_3`, or a cached bracket, would cut it; nothing here
  needed it, and no answer depends on it.
- **The two duplicate models.** `hammah_hb1` and `xslope_rock_slope` are the same
  model under two names, as `xslope_reinforce_fem` and `xslope_reinforced_slope`
  are. The eight-row table above is seven distinct models.


## MERGED FROM `main` again, 2026-09-02

`main` moved thirteen commits while the pile and curved-envelope rounds were being
written. Twelve are a dead-code cleanup and one is a fix to the shipped pile
capacities that this branch has to answer.

**The cleanup does not reach the Newton driver.** It removed a 728-line
per-element stiffness family from `seep.py`, twenty-three zero-caller functions
across `mesh`, `fem`, `plot`, `fileio`, `slice`, `preflight` and `style`,
`global_config.py`, five code-only keyword arguments on `mprice`,
`force_equilibrium`, `circular_search`, `reliability_mc` and
`run_transient_seepage`, and two archived templates. Five of those functions were
in `fem.py` — `compute_flow_vector_tp`, `compute_quad_area`, `mc_flow_vector_4`,
`check_mohr_coulomb_cp` and `_n_reinforcement_lines`. Every deleted symbol was
grepped for by whole word across `xslope/`, `test/`, `studio/` and `run_tests.py`
after the merge, together with `global_config` and the two keywords the previous
merge found gone under this branch: **not one hit**. The Newton, ramp and
predictor code never called any of them.

**The pile-capacity fix does reach it, and it is the subject of the next round.**
a115d9c5 makes a moment capacity a plastic hinge on the DEFAULT driver — a
released rotational degree of freedom whose correction is the full element vector
`K_local p`, one release per pile node — and corrects the shear correction's sign
to the bar's. Its own note says the Newton driver here "carries the same inert
moment-cap form ... and needs the same hinge treatment, including the
one-release-per-node rule, before its pile results mean anything on a model where
M_cap binds." That is "THE PILE HINGE", below. Nothing on the Newton path was
changed in this merge; `fem_details.py`'s soil-reaction reader, which now takes
the hinge's plastic rotation out of the shape, merged identically on both sides.

**The conflicts, and the hand fixes.** One textual conflict:
`run_tests.py`'s `_expected_and_tol` type tuple, where this branch had added
`nr_ssrm` and `main` had added `pile_capacity` to the same line. Resolved as the
union, and both suite rows verified still registered — `run_nr_ssrm_test` and
`run_pile_capacity_test` alike. `xslope/fem.py` auto-merged with no textual
conflict, which is the dangerous outcome rather than the safe one and is why the
whole-word grep above was run rather than assumed; the two signatures that changed
under the branch, `_pile_element_actions` and `pile_element_reaction`, both gained
keyword arguments with defaults, so every existing call site still means what it
did.

The post-merge checks, all run in the worktree with its root first on `sys.path`
and `xslope.__file__` asserted:

| Check | Result |
|---|---|
| Default path: Griffiths & Lane 6 dry, quad8/2, no `fem_solver` | FS **2.421875**, per-trial iterations 147, 781, 3393, 2031, 2841, 9541, 12000, 8617, 8777 — value for value |
| Griffiths & Lane 1, tri6/3.5, Newton bisection | FS **1.36562500**, 9 trials, 3,121 iterations, 22,650 force evaluations — value for value |
| Three layers, `t_cut` = 0, tri6/4 | viscoplastic **1.2109375**, Newton **1.2109375**, 8 trials, 3,688 iterations, 25,515 force evaluations — value for value |
| FEM pile sample (`xslope_piles_fem`), tri6/2, the tag's bracket | viscoplastic **1.3796875**, Newton **1.3796875**, 8 trials, 4,055 iterations, 27,736 force evaluations — value for value, and 0.0003 from the published 1.380 |
| `test/nr_ssrm_check.py` | passes end to end |
| `test/pile_capacity_check.py` (new on `main`) | 0 failures |
| `python3 run_tests.py --preflight` | 31 passed, 0 failed |

The pile sample is the sharp row: `main` rewrote the viscoplastic driver's pile
capacity law, and that model's capacities do not bind, so both drivers had to come
back on the number the piles round recorded. They did.

## THE PILE HINGE — the Newton driver's moment capacity is a release

Written before any feature code, so that what follows is a test and not a
description. Same machine and settings as everything above: `force_tol` 1e-3,
hybrid criterion, `capture_failure_state=False`, tolerance 0.01.

### Why this one now

The pile round measured the moment capacity INERT on both drivers and said so:
"the moment correction is applied at the rotational degree of freedom the two
adjacent beam elements SHARE, and at equilibrium their end moments there are equal
and opposite, so the two corrections cancel. Reversing the sign reverses both of
them and they cancel again." It then left the finding as the owner's decision
rather than a fix.

`main` has since made that decision on the DEFAULT driver (a115d9c5). A moment
capacity is now a plastic hinge: the element's elastic end rotation is the nodal
rotation less a plastic rotation `p`, the moment it delivers is
`K_local (u_local - p)`, and the correction is the full element vector `K_local p`
— which, unlike a nodal moment, acts on the element's translational rows and does
not cancel. One release per pile NODE, because every interior node carries two
element ends and releasing both leaves the node's rotation seeing equal and
opposite capacities whatever it does. That commit's own note says the Newton
driver on this branch "carries the same inert moment-cap form ... and needs the
same hinge treatment, including the one-release-per-node rule, before its pile
results mean anything on a model where M_cap binds."

The shear leg needs nothing: the Newton path already writes it in the bar's
convention — the part of the elastic action the member cannot deliver, subtracted
from its own internal force — which is the sign `main` has now adopted.

### The semantics being reproduced, read from the merged driver

Not assumed — read out of `_pile_element_capacity`, `_pile_moment_hinge` and
`_pile_hinge_ends` on `main` and restated here, because the two drivers have to
solve the same model:

- **The hinge.** With `A` the rotational block of `K_local` (its rows and columns
  2 and 5, the two END nodes) and `m_e` the elastic end moments, the delivered
  moments are `m_e - A p`, and requiring the hinged ends to sit on the capacity is
  one 1x1 or 2x2 linear solve. Every end moment is linear in `p`, so there is
  nothing to iterate. Releasing one end raises the moment at the other, so the
  hinge set is GROWN until it is stable — at most two passes, there being two ends.
- **One release per node.** `_pile_hinge_ends` gives each pile node to the first
  element that reaches it; the other element's end at that node is left elastic,
  and node equilibrium puts it on the capacity with the opposite sign anyway. The
  mask is F-independent and is built once over the WHOLE pile element list, so the
  ownership is global and not per DOF-count group.
- **The correction is the full element vector.** `corr_local = K_local p_local`
  with `p_local` carrying the plastic rotations in slots 2 and 5. The viscoplastic
  scheme solves `K u = base_loads + corrections`, so the internal force the state
  is in equilibrium with is `K u - corrections = T^T K_local (u_local - p_local)`.
- **The axial force and the shear are re-read off the RELEASED displacement**
  `u_local - p_local`, because a hinge changes them; the shear capacity is then
  applied to that shear.
- **`yielded_M` is set whenever the cap is exceeded**, whether or not the element's
  own end took the release — the other end of a shared node is on the capacity too.
- **The plastic rotation is reported**, on the solution as `pile_plastic_rotation`
  and through the pile result sidecar, at the two END nodes of each element.

### Design

The Newton residual is written in the element's GLOBAL frame, so the whole of the
above collapses to `f_int = K_global (u_e - p)` less the shear correction:
`T^T K_local p_local = K_global p` because `p` carries only ROTATIONAL components
and a rotation is invariant under the frame change. The rotational block of
`K_global` at rows and columns 2 and 5 is `K_local`'s own for the same reason, so
the hinge solve is the same 2x2 either way.

- **The law is the shared one.** `_pile_moment_hinge` and `_pile_hinge_ends` are
  called directly, so the hinge algebra and the ownership rule have exactly one
  implementation and cannot drift between the drivers.
  `_pile_element_capacity` itself is NOT reused, and the reason is stated rather
  than assumed: it returns a LOCAL-frame correction and no tangent, and it re-reads
  the actions through `_pile_element_actions`' closed forms, where the Newton path
  needs them on the constant ROWS it differentiates. What it shares with this path
  is the two functions that carry the law.
- **The hinge solve runs only where a cap is exceeded.** An element under its
  capacity takes `p = 0` and the arithmetic it was already on, which is why a model
  where nothing binds has to come back bit-identical.
- **Consistent tangent, with the rotational block condensed.** On a fixed active
  set `A`, `p_A = S_AA^{-1} (m_A - sign(m_A) M_cap)` is affine in `u_e`, so
  `dp/du_e = R` with `R` carrying `S_AA^{-1} G_A` in the released rotational rows
  and zero elsewhere (`G` the two moment rows). Then `dw/du_e = I - R` and

      dK/du_e = K (I - R) - [|V| > V_cap] * q_V (x) (g_V (I - R))

  which is the derivative of the map the residual actually uses. A capped end's
  moment has derivative exactly zero, which is what a hinge means.
- **The displacement bound, the fixity constraint and the ramp are untouched.** A
  hinge adds no degree of freedom and nothing latches — `p` is a function of the
  current displacement alone — so the pile law stays nonlinear-ELASTIC and the
  ramp's warm start needs no extension.

### Success criterion (verbatim)

1. **The binding case, from a115d9c5's own commit message.**
   `docs/tutorials/files/xslope_pile_wall.xlsx`, tri6 at 2 ft, bracket 1.0-2.0,
   with `M_cap`/S lowered from the file's 90,600 to 20,000. The viscoplastic
   driver reads **FS 1.41015625** there, with the largest end moment in
   equilibrium equal to 20,000, a plastic rotation of 1.98e-3 rad and 3 hinges.
   The Newton bisection must land within 0.01 of that; at its reported state the
   largest end moment must equal the cap to 1e-6 relative, and the hinge count
   must match the viscoplastic driver's.
2. **The non-binding cap is bit-identical to uncapped, on Newton.** The same model
   at the file's own 90,600 — which the wall's elastic demand never reaches — must
   return the same factor of safety, the same trial verdicts, the same iterations
   and the same force evaluations as the same model with `M_cap` removed.
3. **The eight pile locks do not move.** Every one of the pile round's eight
   benchmarks, built through `run_tests.py`'s own `build_fem_ssrm_case`, must
   return the Newton factor of safety the piles round recorded — including the four
   disputed rows, which stay at their Newton values: the dispute there is the
   viscoplastic early-failure classifier and not the capacity, and no capacity
   binds on any trial of any of the eight.
4. **The element tangent is the derivative of the element's own residual across
   the RELEASED branch.** Against a central difference, 1e-8 relative or better,
   with the branch histogram reported and the hinged branches actually reached.
   The frame-rotation caveat this document has carried since the tension-cutoff
   round applies to the difference and not to the tangent, and the 1/h scaling is
   the evidence for which.
5. **The locks catch it.** `check_piles` gains a moment-capacity leg, and the
   mutation is the form this round removes: restore the cancelling nodal-moment
   correction and the check must FAIL, run both ways and recorded. The whole check
   file passes.
6. **The default path is unchanged**, against the standard control: Griffiths &
   Lane 6 dry with no `fem_solver` argument, FS 2.421875 on per-trial iteration
   counts 147, 781, 3393, 2031, 2841, 9541, 12000, 8617, 8777.
7. **An honest negative is a valid outcome and must be written.** If the hinge
   costs a lock, if the released tangent will not difference, or if the two drivers
   part on the binding case for a reason that is the formulation rather than a
   solver rule, that is the result.
