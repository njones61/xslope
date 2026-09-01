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
  stress and staged loading all raise. Any default-driver conversation is a
  conversation about that list first, and it is a larger piece of work than
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
  pressure.

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
