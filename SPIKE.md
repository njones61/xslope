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
value for value.

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
