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
