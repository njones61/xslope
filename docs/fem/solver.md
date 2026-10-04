---
title: "Finite element solver — XSLOPE"
description: "Finite element trial settings, viscoplastic iteration, convergence tests and shear strength reduction in XSLOPE."
---

# Solver

## Run settings

In [XSLOPE Studio](../studio/analysis.md#finite-element-fem), build a mesh and choose
**Run FEM…** for a single strength-reduction trial or an SSRM search.

![Run FEM dialog](../studio/images/analysis_run_fem_dialog.png){width=818}

The table maps each dialog control to its input and API setting. **—** means no Excel
cell or no argument on that function; SSRM-only arguments belong to `solve_ssrm()`.

| Dialog control | Excel input | API setting | Meaning |
|---|---|---|---|
| Analysis | — | `solve_fem()` or `solve_ssrm()` | [Single trial](#the-solve_fem-function) or [strength reduction search](#methodology). |
| F (single) | — | `solve_fem(F=)` | [Strength reduction factor](#the-solve_fem-function), default 1.0. |
| F min (SSRM) | `main!D21` | `solve_ssrm(F_min=)` | [Starting bracket](#the-solve_ssrm-function), default 1.0. |
| F max (SSRM) | `main!D22` | `solve_ssrm(F_max=)` | [Starting bracket](#the-solve_ssrm-function), default 2.0. |
| Tolerance (SSRM) | — | `solve_ssrm(tolerance=)` | [Final bracket width](#the-solve_ssrm-function), default 0.01; not the displacement convergence tolerance. |
| Max iterations per trial | — | `max_iterations=` (both functions) | [Trial budget](#creep-trend), default 12000. |
| Iteration ceiling | — | `max_iterations_ceiling=` (both functions) | [Hard limit on budget extensions](#creep-trend), default 50000. |
| Accelerate convergence | — | `accelerate=` (both functions) | [Lengthen admissible iteration steps](#creep-trend); `False` runs the ordinary iteration. |
| Side BC | `main!D23` | `build_fem_data()` reads `slope_data['side_bc']`; no solve argument | [Side restraint](#side-boundary-restraint), rollers by default. |
| K0 initial stress | `main!D16` (blank disables) | `k0=` (both functions) | [At-rest initialization and equilibration](#in-situ-equilibration), off when unspecified. |
| K0 | `main!D16` | `k0=` (both functions) | [At-rest lateral stress coefficient](#in-situ-equilibration); the dialog value starts at 1.0 when enabled. |
| Reduce the tensile cap with F (Tension SRF) | `main!D17` | `tension_srf=` (both functions) | [Reduce a stated positive tensile cap](#tensile-strength-in-ssrm), on by default. |
| Failure criterion | — | `failure_criterion=` (both functions) | [Trial classification](#ssrm-failure-criteria); the API default is `hybrid`. |
| Ignore surficial (skin) failures | — | `min_slip_depth=None` disables (both functions) | [Depth filter](#surficial-skin-failures-and-the-minimum-slip-depth-filter), off by default. |
| Min slip depth | — | `min_slip_depth=` (both functions) | [Depth below ground](#surficial-skin-failures-and-the-minimum-slip-depth-filter), in model length units. |
| SSR exclusions… | `polygon` sheet SSR zones; no cell for the material picker | `solve_ssrm(ssr_exclude=)`; `ssr_zone=` for a search polygon | [Hold selected regions at full strength](#ssr-exclusion-zones). |
| Capture failure-state mechanism | — | `solve_ssrm(capture_failure_state=)` | [One at-failure solve after the search](#the-solve_ssrm-function), on by default. |
| Capture margin | — | `solve_ssrm(capture_margin=)` | [Capture factor above FS](#the-solve_ssrm-function), default 0.15 of FS. |
| Set capture iteration budget | — | `solve_ssrm(capture_max_iterations=None)` uses the automatic budget | [Override the capture budget](#the-solve_ssrm-function). |
| Capture max iterations | — | `solve_ssrm(capture_max_iterations=)` | [Capture ceiling](#the-solve_ssrm-function); automatic budget is `max(max_iterations, 3000)`. |

The [Studio page](../studio/analysis.md#finite-element-fem) describes model checks,
running, cancellation, continuation and result controls. For a transient seepage model,
the additional [Seepage time](../studio/analysis.md#seepage-time) control selects the
instant used for pore pressures.

### Side boundary restraint

Side BC applies to single trials and SSRM searches alike. Set it before
`build_fem_data()`; the solve uses the resulting restraints. Rollers prevent horizontal
movement and leave vertical movement free; Fixed restrains both. The
[boundary-condition section](overview.md#what-xslope-assigns-automatically) describes
the fixed base, side edges and loads, including the added restraint of fixed sides.

## Elastic-plastic behavior: the viscoplastic algorithm {#elastic-plastic-behavior-viscoplastic-algorithm}

The **viscoplastic algorithm** of
[Griffiths & Lane (1999)](https://doi.org/10.1680/geot.1999.49.3.387) and Smith & Griffiths (2004)
returns stress a Gauss point cannot carry as a body load built from accumulated viscoplastic
strains. The elastic stiffness matrix is assembled and factorized **once**, then reused by
back-substitution for every iteration of every strength reduction trial.

![fem_ov_viscoplastic_loop.png](images/fem_ov_viscoplastic_loop.png){width=700}

### The iteration {#viscoplastic-iteration-process}

Stress is carried in the 4-component plane-strain form of Smith & Griffiths (their nst = 4), with
$\sigma_z$ explicit so the algorithm can relax it through plastic $\varepsilon_z$. At each Gauss
point on each iteration:

>- Total in-plane strains come from the current displacements, $\{\varepsilon\} = [B]\{u_e\}$, and
>  the **elastic** strains are $\{\varepsilon\} - \{\varepsilon^{vp}\}$, with
>  $\varepsilon_z^{el} = -\varepsilon_z^{vp}$ (total $\varepsilon_z = 0$). Using the elastic strain
>  rather than the total strain accounts for the stress relief already taken by plastic flow.<br>
>- The stress $\{\sigma\} = [D_e^{4}]\{\varepsilon^{el}\}$ is reduced to the invariants
>  $\sigma_m$, $\bar{\sigma} = \sqrt{3J_2}$ and $\theta$, and the yield function is evaluated in
>  invariant form.<br>
>- Where $f > 0$, a viscoplastic strain increment
>  $\Delta\varepsilon^{vp} = f \cdot \partial Q/\partial\sigma \cdot \Delta t$ is accumulated, using
>  the non-associated plastic potential with dilation angle $\psi = 0$ (no plastic volume change).
>  Within about 0.7° of the Lode-angle corners ($|\sin\theta| > 0.49$) the $\theta$-dependence is
>  frozen at the corner value, keeping the flow direction finite where $\tan 3\theta$ becomes infinite.<br>
>- The accumulated strains form the body-load correction
>  $\{F\} \mathrel{+}= \sum_{e} \int [B]^T [D_e] \{\varepsilon^{vp}\} \, dA$, and the system is
>  re-solved with the existing factorization.

The **pseudo-time step** $\Delta t$ is a numerical parameter, not physical time, taken from Smith &
Griffiths' Program 6.1 as $\Delta t = 4(1+\nu)/(3E)$. The alternative Mohr-Coulomb stability bound
$4(1+\nu)(1-2\nu)/[E(1-2\nu+\sin^2\phi)]$ drives a sustained limit cycle at Gauss points
in mild effective tension beneath reservoir loading; the smaller value is in the stable regime. The
per-iteration displacement increment scales with $\Delta t$, so the convergence tolerance and the
failure criterion are calibrated jointly with it — see the warning on `dt_scale` below.

A **tension cutoff** runs as a second viscoplastic yield surface through the same mechanism; because
it mainly affects SSRM results rather than ordinary stress analyses, it is described under
[Tensile strength in the SSRM](#tensile-strength-in-ssrm).

### Convergence criterion

A displacement-change test alone is not sufficient: a slope past its critical strength reduction
factor can creep slowly enough that the per-iteration change drops below any tolerance while the
slope is in fact failing. XSLOPE therefore requires **two conditions simultaneously**:

1. **Displacement settled** — Smith & Griffiths' CHECON test, the maximum-norm relative change
   between iterations:

>>$\dfrac{\max_i |U_i^{(k+1)} - U_i^{(k)}|}{\max_i |U_i^{(k+1)}|} < \text{tol}$

2. **Force equilibrium** — the criterion of
   [Dawson, Roth & Drescher (1999)](https://doi.org/10.1680/geot.1999.49.6.835): every node's
   out-of-balance force, normalized by the gravitational body force acting on *that* node, below a
   tolerance:

>>$\displaystyle\max_i \dfrac{|\,\mathbf{r}_i\,|}{|\,\mathbf{f}^{\,grav}_i\,|} < \text{force\_tol}$

Because this is an initial-stress scheme, each solve satisfies
$\int B^T D(Bu - \varepsilon^{vp})\,dV = F_{ext}$ *exactly* using the previous iteration's plastic
strains. What is still out of balance is therefore the amount by which the viscoplastic body load is
**still changing**. When plastic flow ceases that increment decays to zero; when the slope is
failing, flow never ceases and the increment plateaus at a non-zero value that feeds displacement
indefinitely.

The test is **local**, which makes it reliable. The increment is non-zero only at nodes
adjacent to Gauss points that are still flowing, so material that merely sits in equilibrium — a
deeper foundation, a longer runout — contributes exactly zero and cannot shift the maximum. A global
norm ratio measures the mechanism against the weight of the entire mesh and offers no such
protection: enlarging the domain changes the denominator without changing the slope. The denominator is
a **lumped** tributary weight, $\sum_e \gamma_e A_e / n_e$ over the elements touching the node, not
the consistent nodal gravity load — the consistent load is exactly zero at a tri6 corner node and
slightly negative at a quad8 corner, which would make the ratio there meaningless.

The converse does **not** hold: a plateau above the tolerance is not by itself evidence of failure.
A slope can stand perfectly still while the residual stalls above an absolute threshold it never
reaches. The size of the gap between the two regimes is problem-dependent — several orders of
magnitude on a Hoek-Brown slope, about two on the Griffiths & Lane benchmark, and on a Mohr-Coulomb
slope with a non-associated flow rule it can close entirely. The
[hybrid criterion](#2-hybrid-hybrid-default), the default, checks the displacement field
directly in that case.

Because the tolerance is absolute, four things affect it:

>- **The yield-surface limit cycle (why `oob_window` exists).** A *one-iteration* increment does not
>  decay on a settled slope: Gauss points resting exactly on the yield surface flip their flow
>  direction on alternate iterations, producing a clean **period-2** oscillation in the viscoplastic
>  body load whose amplitude is proportional to $\Delta t$. Damping the timestep shrinks it but never
>  removes it, so with a one-iteration window a stable slope is reported as failing forever.
>  Averaging the increment over `oob_window` iterations (default 10) cancels the mode exactly while
>  leaving genuine plastic drift untouched. The verdict is insensitive to the width — 10, 50 and 200
>  agree on the same $F$ at the same iteration — so this rejects a specific numerical mode rather
>  than tuning a threshold.<br>
>- **Iteration count.** The test demands *actual* force equilibrium rather than a decayed rate, and
>  displacements settle long before the per-node maximum does. Budget roughly **3× the iterations** a
>  rate-based criterion needs; a ceiling set too low silently truncates a converging solve and reports
>  it as failure, biasing FS **low**.<br>
>- **Element size.** The ratio scales roughly as $1/h$ — the numerator is an internal-force residual
>  ($\sim\sigma h$), the denominator a body force ($\sim\gamma h^2$). A coarser mesh narrows the margin.<br>
>- **Timestep scale.** The residual is the *increment* of the viscoplastic body load and is
>  proportional to $\Delta t$. Shrinking `dt_scale` shrinks the residual without making the slope any
>  more stable, so a failing state can be driven under an absolute `force_tol` and reported as
>  converged. Leave `dt_scale` at 1.0, and never lower it to make a model converge.

**Defaults and budgets.** $\text{tol} = 10^{-3}$ (displacement) and
$\texttt{force\_tol} = 10^{-3}$ (force equilibrium, Dawson's published value), with a budget of
`max_iterations` = 12000 iterations per trial. 1500–4000 iterations is normal well below failure —
consistent with Griffiths & Lane's reported 792 just below their Example 1 failure point — but the
count climbs steeply as a trial approaches the critical factor and as the mesh is refined: the same
reinforced slope reaches equilibrium at $F = 1.25$ in 5,054 iterations at 2.5 ft element size and
16,242 at 1 ft.

#### Trials that reach the iteration limit {#creep-trend}

A trial that reaches `max_iterations`
without converging is classified by the trend of its movement over the second half of its iterations,
split into five equal blocks: how far max&#124;u&#124; moved in each block, and the ratio of each
block's movement to the one before. The ratio is a ratio of movements, not a count of iterations, so
two runs that cover the same ground at different paces are classified the same way.

| Trend over the window | What happens |
|---|---|
| **Dying away**: every block moved forward, none moved more than the one before, and the movement shrinks at a steady ratio below 0.9 a block | The [Newton corrector](#finishing-a-trial-with-the-newton-corrector) is asked for the balanced state from where the trial is. A state the corrector certifies (force, yield and, on a jointed model, the hold test) stands, and the trial's `stop_reading` records the reading and the corrector's result. Where the corrector does not certify a state, the trial is decided as below |
| **Holding steady or growing**: the window moved at least 0.02 elastic displacements at a ratio of 0.9 or more a block; on a jointed model also the joint slip gaining 2% of itself or more at a rate not falling below 0.9 of the half-window before | `exit_reason = 'not_slowing'`, `FAILED`: the slope is sliding |
| **Still** (under $10^{-4}$ elastic displacements over the window) or **unclear** | The [hybrid classifier](#2-hybrid-hybrid-default) decides, as for any trial stopped at the limit |

Below `max_iterations_ceiling` (default 50000) a movement still dying away, or not yet clear, is given
another `max_iterations` worth and read again at the new limit; a trial holding steady, growing or
still stops at the limit. The dying-away reading, with its corrector attempt, is also taken at every
block end once the window is full (and, on a jointed model, past 5,000 iterations), so a creeping
trial can be certified before its limit. A jointed trial holding steady is counted as sliding from the
last tenth of `max_iterations` on, and never earlier: see
[the joint verdict](#the-joint-verdict) for the standing trial that looks like a slide at a third of
its iterations. Every level is one the classifier or the joint verdict already uses: 0.9 is the joint
verdict's "not decaying", 0.02 elastic displacements the classifier's "still moving", $10^{-4}$ the
joint verdict's "stopped".

**Acceleration on jointed models.** On a model with joints the iteration is accelerated by default (`accelerate=None`; pass
`accelerate=False` for the ordinary iteration, `True` to accelerate a model without joints).
From the trial's last corrector checkpoint on, each iteration's step is multiplied by a factor
between 1 and 50 read from the last two steps (Irons and Tuck), with the plastic strains and
the joint slip scaled with it; a step that would change a joint pair's open or slipping state,
trip its residual latch or add dilation is taken at the ordinary length. The balanced state is
the ordinary iteration's; what changes is how many iterations reach it. Over the 32 jointed
verification rows the answers are the same, and the whole set runs about a fifth faster. The
K0 in-situ solve and the hold test always run the ordinary iteration, and the Log's opening
lines for a run show whether acceleration was on.

Whether a trial dying away is certified depends on how close to rest it has come. On the
[Griffiths and Lane Example 2 slope](../verification/ssrm.md#verification-griffiths2) at
$F = 1.34375$ no state is certified at the 300, 1,000 or 3,000 checkpoints or at the block ends
before 11,200 iterations; at 11,200, with each block of 1,600 iterations moving the slope about 0.71
times as far as the one before, the corrector certifies the state at 1.85 times the elastic
displacement and the trial stands. On [RS2-28a](../verification/rs2.md#rs2-28) at $F = 1.7$ the
movement does not die away: at the 16,000-iteration limit the slope has moved 13 times the elastic
displacement, each block still moving it 96% as far as the one before, no state along the way is
certified, and the trial is counted as sliding.

**Inconclusive trials.** A trial that reaches `max_iterations_ceiling` still dying away or with no
clear trend, and with its out-of-balance still falling (the mean over the last 500 iterations at
least 1% below the mean over the 500 before), is neither settled nor failed, and it is reported as `exit_reason = 'inconclusive'`. The
[Newton corrector](#finishing-a-trial-with-the-newton-corrector) is applied to such trials, since a
trial still improving at the ceiling is the kind a locally quadratic iteration can finish; with the
default driver an inconclusive trial is therefore rare, and remains only where the corrector also
cannot certify a state. The bisection does
not count it as a failure, because doing so would bias the factor of safety low. It continues below
the inconclusive $F$. `inconclusive` lists such trials and `note` names the last of them in a sentence,
which is also printed to the log. When an inconclusive trial is still the top of the final bracket,
the search found no failure, so it reports the factor of safety as at least the bracket's bottom,
the highest strength the slope was confirmed to stand at (`fs_is_lower_bound = True`), because a
midpoint would assume a failure that no trial showed. Raise the ceiling or Max iterations per trial to
go further.

**Continuing a run with a higher limit.** Where the iteration limit is what stopped the trial at
the top of the final bracket (undecided at the limit, or counted failed while still slowing), the
run can be continued from where its trials stopped:
`solve_ssrm(fem_data, resume=result, max_iterations=N)`, or **Continue with a higher limit…** in Studio. The search walks the path a
fresh search at the new limit walks, from the original bracket. A trial on it the earlier run
decided is reused; one the earlier run left unfinished goes on from the iteration it stopped at,
with its displacements, plastic strains, joint state and every history the stopping rules read,
and reaches the state a trial run straight through at the new limit reaches; any other trial is
solved as usual. A trial the path never reaches keeps its earlier record. `trials` holds each
trial once, a continued one with `resumed_from`, `result['resumed']` lists the trials reused,
continued and solved afresh, and the closing summary names the F the run was continued from and
the new limit, with the wall time of both parts together. Continuation is not offered for a top trial that ended on the yield
check, since a higher limit does not change that result. The trials' end states are held in memory for the session the run was made in
(`result['resumable']`) and are never saved, so a run read back from its files cannot be
continued.

A **no-progress plateau** — 1500 iterations without improving on the lowest out-of-balance value
seen by more than 1% — is recorded (`plateau_iteration`, `plateau_ratio`) and **does not stop the
solve**: a plateau alone does not end a trial; the convergence tests, movement rules,
yield gate and iteration limits still apply. A plateau describes the residual and does not show whether the slope fails,
for two reasons: the residual is **not monotone** — a reinforced slope
whose out-of-balance sits at twice `force_tol` around iteration 9,500 climbs back an order of
magnitude and then converges at 16,242 — and the iterations a trial needs **grow with mesh
refinement** while a fixed window does not, so on that model a 1500-iteration window was 30% of the
required work at 2.5 ft element size and 9% at 1 ft. Stopping on the plateau would report those
trials as failed, and the bisection would close on the false failure, biasing the factor of safety
low by 18% on the finest mesh. The cost of letting trials run is borne by trials that do fail: they spend
the whole budget unless the rule below closes them sooner, which is why `max_iterations` should be
set to what the model needs rather than left to absorb hopeless trials.

**A trial whose movement is clearly running away is declared failed early** — max|u| past 8 times
that trial's own elastic displacement and still growing, or the out-of-balance flat over the last
2000 iterations while the field gains a whole elastic displacement — and stops there with
`exit_reason = 'diverging'` instead of spending the rest of its budget (`early_failure=False` turns
it off; the at-failure capture solve does not use this rule, and ends instead at the distance limit
described under `capture_failure_state` below or at its own budget). Both thresholds are measured in
the trial's own elastic response, and both sit far outside the range occupied by trials that go on
to reach equilibrium — which near the critical factor grow past five times elastic with a flat
residual — so the rule catches only gross runaways. On the default driver the runaway signal
triggers a [corrector](#finishing-a-trial-with-the-newton-corrector) attempt before it ends
anything: a displacement ratio taken on an iterate is weaker evidence than an equilibrium the
corrector can certify, so the rule applies only where that attempt fails.

The **displacement limit** (`max_disp_factor`) is disabled on the default criterion, and
deliberately so: its yardstick is the height of the *mesh*, not of the *slope*, so it loosens as a
model is given a deeper foundation. The force-equilibrium test has no such dependence.

**Submerged boundaries** converge like any other problem under the effective-stress formulation
combined with consistent boundary-load integration: the submerged soil carries its buoyant weight,
the flooded surface skin is in compression, and sub-critical trials reach true equilibrium (the G&L
Example 6 dam at $F = 1$ settles in a handful of iterations). A useful check on any submerged model
is a single solve at $F = 1$: flooded ground at working strength must settle quickly with an
essentially elastic strain field, and if it does not, suspect the inputs — loads inconsistent with
boundary pore pressures — rather than the solver knobs. Quadratic **triangles** (tri6) are preferred
over quad8 for this problem class, because the 2×2 reduced-integration quad has a zero-energy
hourglass mode that persistent near-surface forcing can excite.

#### Classifying trials on a jointed model {#the-joint-verdict}

On a model with [interface (joint) elements](reinforcement.md#joints-without-reinforcement) the
force-equilibrium test above is measured almost entirely on the joints. The residual is the
increment of the viscoplastic body load, and on a jointed model that increment is the shear traction
the interface could not hold on the last sweep: measured on a jointed rock slope whose soil residual
had fallen to $3\times10^{-8}$, the same solve's joint residual was $3\times10^{-2}$ — six orders of
magnitude apart on one mesh, in one solve.

That matters because the joint residual has a mode of its own. A pair sitting at its slip limit
alternates between slipping and sticking as the stresses around it change, and the resulting
oscillation in the body load is **not** the period-2 flicker `oob_window` was built to cancel: its
period is tens of sweeps, so the ten-sweep average passes it straight through. On a settled jointed
slope the joint residual therefore oscillates in a band with a **flat mean** — between
$1.2\times10^{-2}$ and $4.7\times10^{-2}$ on that slope, its mean constant to 0.04% over the last
half of the solve — and no iteration budget brings it under any tolerance. Such a trial reaches its
ceiling with nothing wrong with the slope and nothing left for the solver to do.

A jointed trial that neither converges nor fails is therefore classified from the **interface**
rather than from the displacement field, over the trailing half of its sweep history, in two ratios that carry no
length, stiffness or mesh size: the slip gained over the window as a fraction of the slip already
there, and the displacement gained over the window in the trial's own elastic displacements.

| Evidence | `exit_reason` | Verdict | Effect on the bisection |
|---|---|---|---|
| In the last tenth of the budget, or at the limit: the [trend reading](#creep-trend) finds the slip gaining ≥ 2% of itself over the window at a rate **not decaying at all** (the last half-window's rate at least 0.9 of the half before), or max&#124;u&#124; moving at a pace that does not slow, **and** max&#124;u&#124; gaining ≥ 0.02 elastic displacements | `not_slowing` | `FAILED` | Failed — the slope is moving on its joints |
| The slip gains ≤ 0.01% of itself, max&#124;u&#124; gains ≤ 10⁻⁴ elastic displacements, the residual on the nodes carrying *no* joint is under `force_tol` across the window, **and** the joint residual has stopped falling (its window mean at least 0.85 of the previous window's) | `joint_settled` | `JOINT_SETTLED` | **Not** failed: the slope is standing, and the bracket moves up |
| At the hard iteration ceiling only: the whole set of contact states repeats exactly with a period of at most 64 iterations over the last 256, 1 to 4 contacts change state within a period, **and** the field returns to itself after every period to within 5×10⁻⁸ elastic displacements per iteration | `joint_settled` (stop reading `contact_cycle`) | `JOINT_SETTLED` | **Not** failed: the slope stands, force balance not met, and the cycling contacts are named |
| Anything else | unchanged | unchanged | The [hybrid classifier's](#2-hybrid-hybrid-default) verdict stands |

Both tests apply only to a trial that would otherwise spend its whole budget, and only on a
jointed model (the sliding test through the trend reading, which every model has) — a trial that converges, a trial already failing, and every model without a joint
reach neither.

Each condition is strict because a looser version misclassifies trials, as the following cases show.

**Timing of the sliding test.** The sliding test is applied only in the last tenth of the trial's
budget and does not end a trial early. There is no early reading that separates a jointed mechanism from a jointed trial that
converges late. One bracket-edge trial that reaches force equilibrium at 185,381 sweeps has, at
50,000 sweeps, gained 30% of its slip and 1.76 elastic displacements of movement with an
*accelerating* slip rate, and sits at 4.3 times its elastic response — every reading a runaway
produces. It then settles. The test therefore waits until a trial that could still converge would
have done so, and at the limit it replaces an `AMBIGUOUS` classification based on a displacement
ratio with a measurement of the joint slip.

**Rate threshold.** The rate test uses 0.9. A trial that converges after 203,000 sweeps and one
that never converges had gained 16.2% and 18.7% of their slip, and 0.375 and 0.488 elastic
displacements, read at 25,000 sweeps — indistinguishable. They are separated only by the rate, which
falls away steadily on the one that finishes (0.64 at 25,000 sweeps, 0.44 at 100,000, 0.21 at
200,000) and does not change at all on a real mechanism (1.0002). A creep that is only slowing is
left to the displacement classifier, which already classifies it as failed.

**Joint residual test.** The settled classification also requires the joint residual to have
stopped falling. Without that condition, one trial was classified as settled nine thousand sweeps
before it converged: the slip, the displacement field and the soil were already at rest, and only
the joint residual was still falling, by 35% per window.

A `JOINT_SETTLED` trial is **not** a converged trial: it never met the force tolerance and
`converged` stays `False`. It establishes what the bisection needs, that the slope stands at that
strength, with the interface, the displacement field and the soil all at rest. The per-trial record
carries the verdict and its `exit_reason` like any other, so both tests are recorded.

**A contact cycle.** A joint with cohesion and no tensile strength can hold a contact whose balance
sits at zero normal stress. Closed, it carries shear up to its cohesion and the faces part; open, it
carries nothing and the faces touch again. No exact balance exists under the joint law at such a
contact: held closed, the rest of the slope balances with the contact in tension, and held open, it
balances with the faces overlapping. The iteration then repeats the same few contact states without
end while the rest of the slope stays still, and the force tolerance is never met. At the hard
iteration ceiling, where such a trial would otherwise end undecided, the last iterations are tested
as in the table above, and the trial stands. The trial line in the log, the closing summary and the
report state it as "stands: the only movement is 3 contacts cycling (period 8 iterations), no net
movement; force balance not met", and name each cycling contact by its line and location. The
movement bound, 5×10⁻⁸ elastic displacements per iteration, sits a factor of ten above the
measured stuck trial (5×10⁻⁹) and a factor of ten below the slowest failing trial whose contacts
were also cycling (5.4×10⁻⁷).

### Finishing a trial with the Newton corrector

The viscoplastic iteration approaches equilibrium from outside the yield surface, and its
convergence rate is linear, so a trial near the critical strength can spend tens of thousands of
iterations still improving and still undecided. XSLOPE runs a second, locally quadratic iteration on
top of it. The viscoplastic loop drives the solve and builds the plastic history; at a short series
of checkpoints — 300, 1,000 and 3,000 viscoplastic passes — at every block end where the
[trend reading](#creep-trend) finds the movement dying away, and again wherever one of the stopping
rules above would end the trial, the current displacement field and plastic strains are handed to a
single bounded Newton-Raphson solve at full gravity and this trial's reduced strengths. That solve
uses the consistent tangent of the Mohr-Coulomb return map, and from a starting state that has
already followed the load path it reaches equilibrium in tens of iterations where the viscoplastic loop needs
thousands.

A state the corrector returns ends the trial as standing only when it passes three checks:

>- **Force equilibrium** — the per-node Dawson measure described above, below `force_tol`, read on
>  the true residual at full gravity.<br>
>- **Yield** — the largest Mohr-Coulomb violation anywhere in the mesh, as a fraction of the local
>  strength, at or below $10^{-6}$. A field can be in force balance and still carry stress the
>  material cannot hold.<br>
>- **Displacement** — the translational movement below a tenth of the model height, read on the deep
>  degrees of freedom wherever a `min_slip_depth` filter is in force.

On a jointed model a fourth check follows, the **hold test**: the viscoplastic sweep is restarted
from the certified state, with its joint history, and the certificate stands only if that sweep
classifies the state as standing within 3,000 more sweeps having moved it by no more than a hundredth
of an elastic displacement. A Newton equilibrium on a set of contacts can be a saddle — a state the
sweep leaves as soon as it is allowed to — and the hold test separates it from a state the slope
remains in. Models without joints never run it.

**When the corrector does not certify a state.** Where any of the checks fails, or the Newton solve does not
converge, the attempt is recorded and control returns to the viscoplastic loop with nothing about
its state changed — the corrector works on a copy of the displacement field and of every plastic
strain, so an attempt that fails leaves the continuing iteration unchanged. The loop then runs on
to its next checkpoint or to its own exit exactly as it would have. The corrector can therefore turn
a trial that the stopping rules would have ended into one certified as standing, but it cannot make
a trial fail.

A trial the corrector decides carries a record of how, under `corrector`: which checkpoint produced
the certified state, how many viscoplastic passes seeded it, the corrector's iterations and force
evaluations, the three readings against their limits, and every attempt made along the way.
`iterations` counts the seed's passes plus the corrector's, so a trial is charged for all the work
that produced it.

### The yield check

Every solved field carries an admissibility reading taken in invariant form: the largest
Mohr-Coulomb violation as a fraction of the local strength (`max_yield_violation`), the count of
Gauss points more than 1% of their strength outside the surface (`n_yield_above_1pct`), the same
pair for the Rankine tension surface, and where in the mesh the worst violation sits
(`max_yield_at`). It is one pass over the Gauss points and costs no solve.

The strength scale the violation is divided by, $c\cos\phi + |\sigma_m|\sin\phi$, carries an
**absolute floor** of $10^{-4}$ of the model's own overburden scale — the largest unit weight in the
model, $\gamma_{sat}$ included, times the mesh height, so the floor is in the model's stress units
whatever they are. The floor matters
because both terms of that scale vanish together in a cohesionless material near a free surface. A
Gauss point there can carry a fraction of a millipascal of numerical residue over a strength scale
of a few tens of micropascals and read as several times its own strength outside the surface, which
is a statement about the denominator rather than about the slope. The floor is orders of magnitude
below any stress that carries a mechanism, so it leaves a real violation exactly where it was.

That reading is also a check. A viscoplastic state that satisfies both convergence conditions but
sits more than $10^{-2}$ of the local strength outside the yield surface does not end the trial: it
is handed to the corrector, and where the corrector certifies an admissible field the trial stands
on that. The force test cannot see this on its own, because the viscoplastic scheme is in force
balance at every iteration and yield is what it relaxes.

**When the yield check ends a trial.** The gate arms only when both readings show that
the loop has stopped changing: the residual has gone flat — 1500 iterations without
improving the lowest out-of-balance value by more than 1% — **and** the trailing
displacements have stopped growing, as measured by the hybrid criterion's growth test.
Until then, a corrector refusal leaves the trial running through its checkpoints,
runaway watch and budget, and the yield reading is taken again on the next state.
Once the gate is armed, the corrector makes a final attempt. If it does not certify
an admissible state, the trial ends with `exit_reason = 'yield_gate'` and remains
**undecided**: the bisection handles it on its uncertainty path, not as a failure.
A state that passes the check ends the trial as it otherwise would.

The threshold is looser than the corrector's $10^{-6}$ because the two states are reached in
different ways. A Newton state solves the equations the reading is taken from and measures $10^{-8}$
or better. A viscoplastic state approaches the surface from outside along the relaxation and stops
when the *displacement* increment settles, so what is left of its yield violation is set by a
displacement tolerance and not by a yield one; holding it to the corrector's figure would reject
most of the verification states, which have simply not finished relaxing. $10^{-2}$ is where that
residual ends and unrelaxed yield begins, and it is also the fraction the reported Gauss-point count
is taken against, so the check and the reported count are consistent.

When a state fails the check, look first at the material at `max_yield_at`. In a material with
$c = 0$ the Mohr-Coulomb surface passes through the origin, so any tensile mean stress there is
outside it — a cohesionless soil carries no tensile capacity whether or not a Rankine cap is written
beside it, which is why declaring `t_cut = 0` on such a material changes nothing about the surface.
On a material with real cohesion, a tensile cap is a genuine restriction and declaring one is the
modeling answer; see [Tensile strength in the SSRM](#tensile-strength-in-ssrm).

### Choosing the driver

`fem_solver` names the per-trial driver on both `solve_fem()` and `solve_ssrm()`:

>- **`'auto'`** (the default) — the viscoplastic loop with the corrector and the yield check
>  described above.<br>
>- **`'viscoplastic'`** — that loop on its own: no corrector, no ladder, no yield check. A trial ends
>  on the stopping rules and on nothing else.

Setting `XSLOPE_FEM_SOLVER` selects the driver for a whole process. When the environment rather than
a call argument selects a non-default driver, one warning line is printed, because a shell variable
left from an earlier session otherwise changes every factor of safety in a run without any indication.

The at-failure capture solve that `solve_ssrm()` makes past the critical strength runs the
viscoplastic loop with the corrector off. That solve exists to let the failure mechanism develop for
the deformation figure, and a certified equilibrium there would replace the field the figure is
drawn to show.

### A jointed model {#jointed-model-solver-policy}

A trial on a model with [interface (joint) elements](reinforcement.md#joints-without-reinforcement)
costs tens or hundreds of thousands of sweeps where an unjointed one costs hundreds. The cause is
the iteration, not the joint law. The loop factorizes the elastic stiffness once and returns
everything the material cannot carry as a body load, so its error contracts by the ratio of the
stiffness the assembly holds to the stiffness the material has — and a pair at its Mohr-Coulomb
limit carries a traction pinned at that limit, so its tangential stiffness is zero while the matrix
still holds the full $k_s$. Where slipping and open pairs are all that holds a block up,
the contraction factor is within a ten-thousandth of one: measured at 0.99947 on one rock-toppling
bracket edge, which is 13,000 sweeps to take the error down by a factor of a thousand. A trial that
is *failing* is slow for the mirror-image reason — both residuals go flat within a few hundred
sweeps and the block then slides at a few billionths of an elastic displacement per sweep until it
has moved the eight that the runaway rule reads.

On a jointed model, `'auto'` passes the state to the corrector as on any other model. The
interface's consistent tangent carries the friction cross term
$\partial t_s/\partial\Delta_n = \pm k_n\tan\phi_j$ — non-symmetric, and only the Newton path's
general factorization ever sees it — and the corrector's interface law is seeded with the sweep's
own accumulated slip, dilational opening, residual branch and opening history rather than with a
pristine interface. The two of those that a Newton step can *grow* — a pair that reaches its limit
during the step moves onto its residual branch, and one that slides further rides further up its
asperities — are advanced by the step's own return map, so what the corrector solves is the problem
it was seeded with. It certifies on the same three checks it applies anywhere: one
rock-toppling bracket edge that the plain loop takes 185,381 sweeps to converge is certified at 300
sweeps, at a force residual of $3\times10^{-11}$ with no Gauss point outside its surface, and that
model's whole bracket closes at the same place as with the plain loop, at a thirteenth of the
wall-clock time. An attempt that fails changes nothing, here as elsewhere, so the corrector cannot
make a jointed trial fail.

`joint_newton=False` on `solve_fem()` and `solve_ssrm()` runs the viscoplastic loop alone on a
jointed model, deciding the trial on the stopping rules and on
[the joint verdict](#the-joint-verdict); `fem.JOINT_NEWTON_ON = False` does the same for a whole
process, and `fem_solver='viscoplastic'` turns the corrector off on every model.

A second switch, **off by default**, addresses the *failing* edge, which the corrector cannot help
because it has no equilibrium to certify:

>- **`joint_tangent='slip'`** — the interface relief. The shear stiffness of every slipping pair,
>  and both stiffnesses of every open one, are taken out of the assembled free matrix down to
>  `joint_tangent_factor` (default 0.01) of their elastic value, and the matrix is refactorized
>  whenever the slipping and open set moves. The traction limit, the slip return and the state a
>  trial converges to are untouched. The relief runs as a **predictor** — until its state settles or
>  its own sweep budget runs out — and then switches off and passes the trial to the plain loop,
>  which decides it on its own trace, because every level and rate in the stopping rules is
>  calibrated on that trace. It cuts one benchmark's 196,201-sweep failing edge to 20,991, but it
>  also moves a bracket that the corrector reproduces exactly, so it is an option rather than the
>  default. Inert on a model with no joint.

`fem_solver='newton'` on a jointed model runs from a cold start with no accumulated slip behind it,
and on the rock-joint benchmarks it diverges on trials the viscoplastic loop converges. A divergence
there is reported as a failed trial, not as an error, but the Newton driver cannot bracket those
models.

### Surficial (skin) failures and the minimum-slip-depth filter

On a purely frictional face ($c = 0$) the critical mechanism is a shallow slide running parallel to
the slope, with $FS = \tan\phi / \tan\beta$ — a result *independent of depth*, so the shallowest
surface governs. The per-node force-equilibrium criterion detects this "skin" faithfully, and
because it is the true global minimum the reported factor of safety can sit well below a deeper,
more conventional mechanism, and below published values that report the deeper one. This is
physically correct but often not the engineering question; it occurs mainly on the steep frictional
faces of embankment dams.

The optional **`min_slip_depth`** parameter — on `solve_fem()`/`solve_ssrm()` and on the LEM
searches, **off by default** — excludes any failure shallower than the given depth below the ground
surface. In the FEM it acts on the **failure verdict**, not on the strength: nodes shallower than
the cutoff are left out of the per-node out-of-balance maximum, so a shallow skin alone can no
longer cause the slope to be classified as failing, while a deep-seated mechanism still trips the
criterion through its deep nodes. The filter applies to both tests a trial is decided on: the
[corrector](#finishing-a-trial-with-the-newton-corrector) measures its displacement bound on the
same deep degrees of freedom, so a skin that is sliding while the mass beneath it stands cannot
cause an otherwise admissible state to be rejected either. Nothing is held at full strength and no
element is masked — the skin still yields, but it no longer decides the classification. It is a
**run option rather than a file setting**:
pass `min_slip_depth=` to a solve or to a `circular_search()` / `noncircular_search()` call, or set
**Min slip depth** in Studio's Run FEM dialog. A depth greater than the depth of the mesh is
rejected.

Steering the mechanism by **depth** is one of two ways to keep a competing failure out of the
answer; the other is to steer it by **region**, which is what
[SSR search areas and exclusion zones](#ssr-exclusion-zones) do. Use the depth filter when the
mechanism to exclude is defined by how shallow it is (a face-parallel skin, which no zone boundary
separates from the deep surface because both run through the same material); use a zone when it
belongs to an identifiable part of the model (a stiff foundation, a shell, a bench).

**Choosing a value.** Sweep the depth and look for a plateau. As the depth increases, the factor of safety
holds at the surficial-skin value while the cutoff is still inside the failing band, rises as the
cutoff clears the band, then **flattens onto a plateau**: the deep-seated factor of safety. Any
depth on the flat part returns the same FS, so the choice is robust. Run a handful of depths (say
5, 10, 15, 20, 25% of the slope height) and read the trend:

>- Still rising → the cutoff is inside the surficial band; go deeper.<br>
>- Flat → that value is the deep-seated FS. Report it.<br>
>- Still climbing at a large fraction of the slope height → you are past the real mechanism and are
>  excluding genuine failure; back off to where it plateaued.

A large gap between the filter-off value and the plateau means a surficial skin was governing the
unfiltered result; a small gap means the deep mechanism already governs and the filter can stay off.
On a low fill over soft ground the plateau can sit deep as a fraction of height — the embankment on
soft ground of [RS2-66](../verification/rs2.md#rs2-66), on a 4 m soft band, is on its plateau at a
4 m cutoff on a 10 m fill (1.131 under a 4, 6 or 8 m cutoff alike), while the 162 m Talbingo dam
of [RS2-4](../verification/rs2.md#rs2-4) is already on its plateau by 10 m (its 1.67
downstream-bench skin against a 1.82–1.83 plateau held flat from 10 m out to 30 m). The same
model on a 2 m band has no plateau there — 4 m returns 1.131, and 6 or 8 m return 1.206, a
different and higher value. The plateau must be found on each model; it does not carry over from a
similar one.
Set the same `min_slip_depth` in the LEM search and the SSRM run so both report the same mechanism.

### The `solve_fem()` function

`solve_fem()` takes a FEM data dictionary from `build_fem_data()` and an optional strength reduction
factor, assembles and factors the stiffness once, and runs the viscoplastic loop to convergence or
to its iteration budget:

```python
from xslope.fileio import load_slope_data
from xslope.mesh import build_mesh_from_polygons, get_material_polygons
from xslope.fem import build_fem_data, solve_fem

slope_data = load_slope_data("docs/fem/files/xslope_griffiths1.xlsx")

# quadratic elements are required for a trustworthy factor of safety
mesh = build_mesh_from_polygons(get_material_polygons(slope_data),
                                target_size=6, element_type='tri6')

fem_data = build_fem_data(slope_data, mesh)
solution = solve_fem(fem_data, F=1.0, debug_level=1)

if solution['converged']:
    print(f"Converged in {solution['iterations']} iterations")
    print(f"Max displacement: {solution['max_displacement']:.6f}")
else:
    print(f"No equilibrium after {solution['iterations']} iterations "
          f"({solution['exit_reason']}, verdict {solution['verdict']})")
```

Its principal arguments:

>- **`F`** (default 1.0): strength reduction factor, applied as $c_r = c/F$ and
>  $\tan\phi_r = \tan\phi/F$.<br>
>- **`max_iterations`** (default 12000) and **`tolerance`** (default $10^{-3}$): the iteration budget
>  and the CHECON displacement tolerance.<br>
>- **`max_iterations_ceiling`** (default 50000): hard stop on the extension a trial still dying
>  away, or with no clear trend, is given at `max_iterations` (see the
>  [trend reading](#creep-trend)). Reaching it with the
>  out-of-balance still falling gives `exit_reason = 'inconclusive'`.<br>
>- **`force_tol`** (default $10^{-3}$): the per-node force-equilibrium tolerance; with `oob_window`
>  (default 10) the averaging width that cancels the yield-surface limit cycle.<br>
>- **`failure_criterion`** (default `"hybrid"`): how a non-converged trial is judged — see
>  [SSRM failure criteria](#ssrm-failure-criteria).<br>
>- **`max_disp_factor`** (default 0.1, `None` to disable): displacement backstop as a fraction of
>  mesh height. The SSRM's default path disables it.<br>
>- **`early_exit`** (default `True`): watch the residual for the no-progress plateau described
>  above and report it. The plateau does not end the solve.<br>
>- **`fem_solver`** (default `'auto'`): the per-trial driver — see
>  [Choosing the driver](#choosing-the-driver).<br>
>- **`k0`**, **`min_slip_depth`**, **`tension_cutoff`**, **`elastic_mask`**,
>  **`suction_phi_b`** / **`suction_cap`**: the options described in their own
>  sections; all default to off or to what the input file declares.<br>
>- **`debug_level`** (default 0): 0 silent, 1 summary, 2 per-iteration.

The returned dictionary carries `converged` and `stable`, the verdict metadata (`verdict`,
`u_ratio`, `u_growth`, `exit_reason`), `iterations`, the nodal `displacements` and
`displacements_elastic`, element `stresses` and `strains`, `plastic_elements`, and the 1D structural
element forces — everything `plot_fem_results()` and `export_fem_solution()` need. It also carries
the [yield reading](#the-yield-check) (`max_yield_violation`, `n_yield_above_1pct`,
`max_yield_at`, `yield_flagged`) and, where a corrector decided the trial, the `corrector` record.

## Shear strength reduction method (SSRM)

The SSRM (Matsui & San, 1992; Griffiths & Lane, 1999) reduces soil strength until the
finite element system can no longer find equilibrium under the applied loads, without assuming
a failure surface. The reduction factor at that transition is the factor of safety, consistent
with the limit-equilibrium definition.

### Methodology

Each trial divides both strength components by the trial factor,

>>$c_r = \dfrac{c}{F}$<br>
$\tan \phi_r = \dfrac{\tan \phi}{F}$

reducing $\tan\phi$ rather than $\phi$ so the scheme stays well behaved as the friction angle
approaches zero. As $F$ rises, more Gauss points yield, displacements grow, and at some point the
viscoplastic iteration stops reaching equilibrium at all. `solve_ssrm()` brackets that transition
and bisects it.

![fem_ov_ssrm_sweep.png](images/fem_ov_ssrm_sweep.png){width=760}

The sweep above is the [Griffiths & Lane Example 1](../verification/ssrm.md#verification-griffiths1)
sample file solved at fixed strength reduction factors on a deliberately coarse mesh: trials below
the critical factor settle to equilibrium with small viscoplastic displacement, trials above it
never settle and their displacement runs away. The bisection locates that transition — here 1.36 on
this illustration mesh, against the paper's 1.4 and the finer meshes used for the verification
benchmark. The displacement of a failing trial depends on the iteration budget it was given, which is
why the bisection uses the trial's *verdict* rather than the displacement magnitude.

### In-situ equilibration

Under a **slope** the same field is not an equilibrium. There is no soil column beside the face to
balance the lateral stress there, so a substantial share of the weight is left out of balance and
has to redistribute — about a quarter of it on Griffiths & Lane Example 1.

That redistribution is part of **establishing the in-situ state** rather than of strength reduction,
and the SSRM runs the two as separate steps. A $K_0$ analysis begins with one
**full-strength equilibration solve**: the $K_0$ field settles against the real geometry at
unreduced strength, and every bisection trial then starts from the resulting stress state, with a
zero displacement datum, and reduces strength from there. On a jointed model the state carried
into each trial includes every joint's slip history — how far each pair has slid, whether it stands
open, and its residual and dilation state — so a joint that slid while the slope settled starts
the trial where the settling left it, not pushed back onto its limit as if it had never moved. The
equilibration is solved once and shared by all trials, and its outcome is returned in
`result['k0_equilibration']`.

If the two steps were run together, every trial would repeat the in-situ redistribution against soil
already weakened by $F$ and charge the displacement and plastic strain it produces to the trial.
On Example 1 at $K_0 = 1$, $F = 1.2$, that gives about three times the displacement the strength
reduction actually causes.

Displacements are reported relative to the equilibrated state, because the in-situ travel is an
artifact of imposing a stress field the geometry does not hold in equilibrium, not motion of the
slope. Stresses and structural forces — bar tensions, pile end forces — are functions of the
absolute displacement and are unaffected by where its zero is put.

A single `solve_fem()` call at full strength *is* an equilibration solve. At a reduced $F$ a single
call does both at once, which is the sequencing the SSRM avoids; use `solve_ssrm()` when the two
must be kept apart. If the equilibration does not come back stable, the slope does not stand at full
strength with that initial stress ($FS < 1$): XSLOPE warns, and the bisection proceeds without a
carried in-situ state and finds the sub-unity factor of safety.

For a non-converged trial, note that the displacement scale the
[hybrid criterion](#ssrm-failure-criteria) measures against is the elastic response to the
**applied** load, which is the same quantity with or without $K_0$. What $K_0$ changes is the zero —
a trial's displacement is counted from the equilibrated in-situ state, so it carries only the
movement the strength reduction causes.

### SSRM failure criteria

Four criteria are selectable through the `failure_criterion` argument of `solve_ssrm()`.

#### 1. Non-convergence (`"non_convergence"`)

The classical Griffiths & Lane (1999) approach: bisection on whether the viscoplastic iteration
converges. In XSLOPE "converges" means **true equilibrium** — both the CHECON displacement test and
the force-equilibrium test — so the bisection brackets the genuine boundary between states that
reach static equilibrium and states that creep indefinitely.

The force-equilibrium half is Dawson, Roth & Drescher's, not Griffiths & Lane's, whose own criterion
is the displacement test plus an iteration ceiling. In practice the displacement test alone almost
never discriminates: a slope creeping steadily past its critical factor produces a bounded
per-iteration change measured against a growing total, so the ratio decays and the test passes on
states that are plainly failing.

Validated against Griffiths & Lane Example 1 (FS ≈ 1.40 vs published 1.4), their Example 6 dam
without free surface (≈ 2.4–2.5 vs published ~2.4), and the geogrid-reinforced slope (≈ 1.65 vs the
limit-equilibrium Spencer value 1.59 on the same model).

#### 2. Hybrid (`"hybrid"`, default) {#2-hybrid-hybrid-default}

The same bisection, with one addition: **a trial that fails to reach equilibrium must also show
displacement evidence of failure before the bisection counts it as a failed slope.**
Non-convergence on its own describes the solver; the hybrid criterion also checks the displacement
of the slope.

Two signals are read from the trial's own iteration history, both measured against its **elastic
displacement** — the purely elastic response to the same loads, which every solve already computes:

>- **Scale.** $u_{ratio} = \max|u| / \max|u|_{elastic}$ at the end of the solve, which shows whether
>  the field is beyond elastic scale.<br>
>- **Growth.** The gain in $\max|u|$ over the last quarter of the history, in elastic displacements,
>  which shows whether it is still moving.

`max|u|` is sampled every 10 iterations from the value the CHECON test already computes, so the
instrumentation costs nothing measurable and no extra solves are needed.

| Evidence | Verdict | Effect on the bisection |
|---|---|---|
| Beyond elastic scale **and** growing (or the displacement cap was tripped) | `FAILED` | Failed — same as the default criterion |
| At elastic scale **and** frozen | `STABLE_STUCK` | **Not** failed: the bracket moves up |
| One signal without the other, or too little history | `AMBIGUOUS` | Failed — the default criterion's verdict stands |

On a model with joints these two signals are not enough, because the quantity that separates a
settled jointed slope from one sliding on its joints is the **slip**, which the displacement field
does not report: a rock slope gaining 0.4 elastic displacements of movement every 25,000 sweeps
is classified `AMBIGUOUS` here because its max&#124;u&#124; has not yet passed 1.5× an elastic
response set by loading stiff rock. Such a trial is classified from the interface instead — see
[the joint verdict](#the-joint-verdict).

Requiring *both* signals in each direction keeps this conservative: the hybrid overrides only
where the evidence is unambiguous, and every trial's verdict, $u_{ratio}$ and growth come back in
`result['trials']`, so an override is never silent.

**Full-budget history.** Both signals are calibrated on solves that ran to
their iteration ceiling. A slow runaway takes far longer to become visible in the displacement field
than a residual plateau takes to appear, so a truncated history can look frozen while the slope is
accelerating. A plateau alone does not truncate the history read by the classifier.
A trial can still end under the other stopping rules described above.

**Calibration.** The thresholds are $u_{ratio} \le 1.25$ for "at elastic scale", $u_{ratio} \ge 1.5$
for "beyond it", and a growth of 0.02 elastic displacements over the trailing window for "still
moving". They come from measured behavior: stable-but-stuck trials sit at **1.0–1.1×** elastic and
are frozen there whether the budget is 10,000 iterations or 80,000, while genuinely failing trials
reach **4–21×** and are still growing when the budget runs out. Both signals are ratios, so when the
elastic displacement comes back smaller than $10^{-6}$ of the model height — a level model whose
initial stress is already in equilibrium — the verdict is `AMBIGUOUS` rather than a ratio taken
against rounding noise. The only verdict that does not depend on the elastic displacement is a
trial stopped by the displacement limit, which is an absolute fraction of mesh height and is evidence in its own right.

**Why it is the default.** All 103 FEM benchmarks were solved under both criteria on the same mesh
with the same options. No row comes back **lower** under the hybrid, and all but a few come back
identical to the last digit, because on a healthy model every non-converged trial carries real
displacement evidence and the extra test gives the same result as non-convergence.

| Case | Non-convergence | Hybrid | What the hybrid changes |
|---|---|---|---|
| [Griffiths & Lane Example 1](../verification/ssrm.md#verification-griffiths1) | 1.372 | 1.372 | Nothing — the majority case. Every non-converged trial is beyond elastic scale and still growing, so the displacement test gives the same result as non-convergence on every trial |
| [RS2-62c](../verification/rs2.md#rs2-62) | 0.769 | 0.769 | Nothing, by the other route: the $F = 0.775$ trial that sets the bracket spends its whole budget still at elastic scale ($u_{ratio} = 1.23$) but still moving (growth 0.22), one signal without the other, so its verdict is `AMBIGUOUS` and non-convergence's failed verdict stands |
| [RS2-48](../verification/rs2.md#rs2-48) baseline geotextile wall | *no bracket* | 0.994 | The hybrid gives a result where non-convergence gives none. Under the vendor's $T = 0$ cap the trials are stationary rather than collapsing, so non-convergence has no failure side to bisect: it drives the auto-bracket to its floor and returns no factor of safety, while the hybrid brackets the same model |

Agreement is therefore the usual outcome, and it is reached two ways: on a healthy model the
displacement evidence confirms the non-converged trial has failed, and where the evidence is
inconclusive the hybrid defers rather than overriding.

Pass `failure_criterion="non_convergence"` for the classical Griffiths & Lane verdict; it remains
fully supported, and every criterion returns the same per-trial records.

#### 3. Displacement limit (`"displacement_limit"`)

Bisection on whether the maximum viscoplastic displacement exceeds `max_disp_factor` of the mesh
height within the iteration budget. A simple physical backstop, but its verdict is coupled to the
budget for any state that creeps slowly rather than racing.

#### 4. Displacement catastrophe (`"displacement_increase"`)

Sweeps $F$, locates the sharpest upturn of displacement versus $F$ (the evidence Griffiths & Lane
present as their Figs 2 and 18), and refines around it; related to the average-residual-displacement
criterion of [Sun, Wang & Zhang (2021)](https://doi.org/10.1007/s10064-021-02237-y), and like it
reads a **characteristic point** rather than the global maximum. The point is selected
automatically — after the coarse sweep, the node whose plastic displacement grew fastest between the
lowest and highest $F$ becomes the measurement point and the curve is re-read there — which keeps
the measurement on the mechanism rather than on any localized background deformation that grows at
*all* $F$. A specific point can be supplied through `char_point=(x, y)`.

#### Choosing a criterion

| Problem class | Criterion | Why |
|---|---|---|
| All slope problems, including submerged boundaries and reservoir loading | `hybrid` (default) | Bisection on true equilibrium, with a non-converged trial required to show displacement evidence before it counts as failed. Checked over the whole FEM benchmark set: it matches `non_convergence` on all but a few rows and never returns a lower factor of safety |
| Reproducing the classical Griffiths & Lane (1999) verdict, or a published result obtained that way | `non_convergence` | The same bisection without the displacement-evidence test. Expect agreement with the default except where a slow-but-converging trial is truncated at the no-progress exit |
| Evidence and reporting | `displacement_increase` | Produces the displacement-vs-$F$ curve; read the upturn at the automatically selected characteristic point |

FEM-SSRM and limit equilibrium are different formulations, and some difference in computed factors
of safety is expected; running both — as the verification suite does — is the strongest consistency
check available.

### The `solve_ssrm()` function

```python
from xslope.fem import solve_ssrm

result = solve_ssrm(fem_data, F_min=1.0, F_max=2.0, tolerance=0.05, debug_level=1)

if result['converged']:
    print(f"Factor of Safety: {result['FS']:.2f}")
    print(f"Final interval: {result['final_interval']}")
```

Its principal arguments:

>- **`F_min`** (1.0) and **`F_max`** (2.0): the starting bracket. If the slope fails at `F_min` or
>  stands at `F_max`, the bracket **auto-expands** in steps of **`f_adjust`** (0.25) until it is
>  valid, bounded by **`f_min_floor`** (0.1), **`f_max_ceiling`** (10.0) and **`max_expand`** (20
>  steps each way) — so a wrong guess still finds the factor of safety, and a good guess simply
>  brackets on the first try.<br>
>- **`tolerance`** (0.01): the bisection stops when the bracket is narrower than this. The reported
>  FS is the bracket midpoint (± tolerance/2), or the bracket's bottom as a lower bound
>  (`fs_is_lower_bound`) when the top is an undecided trial; the bracket is returned in
>  `final_interval`.<br>
>- **`grid`** (`None`): bisect over a **fixed global grid** of this step instead of halving the
>  supplied bracket. Because the failure threshold sits between two fixed grid points — a property of
>  the slope and mesh, not of the bracket — every starting bracket then converges to the same cell and
>  the reported FS is independent of the bracket, at the same $\log_2$ cost. Used by the
>  [reliability analysis](../reliability/fem.md) for reproducible results.<br>
>- **`failure_criterion`** (`"hybrid"`), **`convergence_tol`** ($10^{-3}$), **`force_tol`**
>  ($10^{-3}$), **`max_iterations`** (12000), **`max_iterations_ceiling`** (50000): passed to each
>  trial.<br>
>- **`max_disp_factor`** (0.1): the displacement-limit fraction. It is what the
>  `"displacement_limit"` criterion bisects on; the equilibrium-based criteria disable it in their
>  trials.<br>
>- **`dt_scale`** (1.0): multiplier on the viscoplastic pseudo-time step. **Do not lower it to make
>  a model converge** — it shrinks the residual without making the slope any more stable, and can push
>  a failing state under an absolute `force_tol`.<br>
>- **`fem_solver`** (`'auto'`): the per-trial driver, passed to every trial — see
>  [Choosing the driver](#choosing-the-driver). The at-failure capture solve below always runs the
>  viscoplastic loop with the corrector off.<br>
>- **`k0`**, **`min_slip_depth`**,
>  **`ssr_exclude`** / **`ssr_zone`**, **`tension_cutoff_by_material`** / **`tension_srf`**,
>  **`elastic_materials`**, **`suction_phi_b`** / **`suction_cap`**: as described in
>  their own sections.<br>
>- **`n_sweep`** (10): coarse sweep points for the `"displacement_increase"` criterion.<br>
>- **`capture_failure_state`** (`True`) and **`capture_margin`** (0.15): after the bracket resolves,
>  take one extra solve at $F = \text{FS}\times(1 + \text{capture\_margin})$ — with the displacement
>  backstop and early exit off and a generous ceiling — so the unconverged field develops the
>  **at-failure mechanism** (the rotational collapse: crest settlement, toe heave). Right at the
>  critical factor the collapse develops too slowly to become visible in a finite number of
>  iterations, hence the margin. The capture stops once the section has moved 20% of the mesh
>  height, and keeps the last state short of that distance. Past critical the slide moves at a
>  steady rate, so without the limit the size of the drawn displacement would be set by the
>  iteration ceiling alone; the shear band forms long before the section has moved that far. A
>  capture that reaches its ceiling first is kept as it stands. A slope beyond critical never
>  passes through equilibrium, and a reinforcement element drops to its residual only on an
>  equilibrium state, so the capture solve starts from the post-peak set the bracket's failed-edge trial shed to: a layer that gave way at
>  its residual is shown at that residual in the at-failure field (`failed_edge_softened` in the
>  result). The field is returned as `failure_solution` and changes nothing
>  else, so turning it off leaves the factor of safety, the bracket and `last_solution` untouched —
>  which is what the reliability and sensitivity analyses do, since they never draw the field.

The result dictionary carries `FS`, the last converged solution (`last_solution`), `final_interval`,
the per-trial records (`trials`), and — with the capture on — `failure_solution`. Passing
`last_solution` to `plot_fem_results()` shows the near-critical converged state; passing
`failure_solution` shows the developed collapse mechanism.

### SSR search areas and exclusion zones {#ssr-exclusion-zones}

A strength reduction run finds the weakest mechanism in the model, which is not always the
mechanism the analysis is about: a stiff foundation, a bench, a shell or a face skin can hold the
global minimum while the surface of interest lies elsewhere. Both XSLOPE and the vendor codes handle
that by naming **where** the reduction applies, in one of two complementary forms.

>- An **exclusion area** names the part of the model held at **full strength** — everything else is
>  reduced. Use it when the competing mechanism is the one you can point at: "not through the
>  foundation", "not through the downstream shell".<br>
>- A **search area** names the part that **is** reduced — everything else is held at full strength.
>  Use it when the mechanism of interest is the one you can point at: a corridor around a proposed
>  slip surface, one face of an embankment, a single tier of a wall.

Both constrain the answer, so the factor of safety they return is conditional on the constraint —
run the unconstrained case as well. The depth-based alternative, which steers the mechanism without
naming a region, is [`min_slip_depth`](#surficial-skin-failures-and-the-minimum-slip-depth-filter).

**By material name.** `ssr_exclude` takes a list of material names. At every trial factor, every
zone *not* named has its $c$ and $\tan\phi$ divided by $F$ as usual, but a named zone keeps its full
strength and the developing shear band is forced up and out of it. This reproduces RS2's
per-material **Apply_SSR** flag, presented in its interface as an "SSR Exclusion Area".

```python
result = solve_ssrm(fem_data, F_min=1.2, F_max=1.5, tolerance=0.02,
                    ssr_exclude=["Foundation lower"])
```

Names must match a material's `name` exactly; an unknown name raises `ValueError` rather than
silently reducing nothing. [RS2-P4-VP67](../verification/rs2.md#p4-vp67) works through the
constrained/unconstrained pair on a USACE end-of-construction embankment: an unconstrained SSRM of
1.076 on a deep foundation mechanism, against 1.303 with the foundation's lower zone excluded,
on the same toe-circle family RS2 reports at 1.33 under its own exclusion area. In Studio
this is the **SSR exclusions…** button in the Run FEM dialog.

**By polygon, on the input file.** A zone is drawn where the rest of the model is drawn: a row on
the [**polygon** sheet](../usage/input_template.md#ssr-zones) whose **Type** is one of three words.

| Type | Meaning |
|:-----|:--------|
| **`ssr reduce`** | **Search area** — reduce **only inside**. |
| **`ssr hold`** | **Exclusion, full strength** — never reduced inside, but can still yield. |
| **`ssr elastic`** | **Exclusion, elastic** — linear elastic inside, cannot yield at all. |

(Template version 20 encoded the same three as negative Material IDs, −1 / −2 / −3, and those files
still load unchanged.)

Several zones combine by one rule: **the reduced region is the union of the `ssr reduce` zones,
minus the union of the `ssr hold` and `ssr elastic` zones**, defaulting to the whole model when no
search area is drawn. Exclusions therefore always carve out — of a search area they sit inside, or
of the model as a whole — and an interior hole in a search area is drawn by putting an `ssr hold`
(or `ssr elastic`) polygon on top of it.

These rows are **analysis overlays, not geometry**. They are never meshed, never become material
regions and never generate slices; they may overlap one another and cross material boundaries
freely. Membership is decided element by element, by where each element's centroid falls, so a zone
can be added to a finished model without disturbing it — the mesh, the material assignment and the
factor of safety are unchanged unless the zone actually constrains something. The limit-equilibrium
solvers ignore the rows entirely.

**By polygon, at the run.** `solve_ssrm()` also takes an explicit **`ssr_zone`** polygon (a vertex
list) — the programmatic primitive, and what RS2's "SSR Search Area" maps onto when a vendor model
is imported. It has **one sense only, reduce inside**, and it takes precedence over the file's own
zones with a warning when both are present rather than quietly intersecting them. It classifies
elements by the same centroid test, so the same polygon written into the file and passed as the
argument give identical answers.

A vendor **exclusion** polygon passed through `ssr_zone` must therefore be converted to its
**complement within the model outline**; passing it as-drawn reduces exactly the wrong region.
[RS2-4](../verification/rs2.md#rs2-4) is the worked case: RS2 holds the whole downstream benched
shell of the Talbingo dam at full strength, and reproducing that answer through the argument means
passing the complementary ring — everything upstream of the shell — as the search area.

The choice between the two depends on the shape of the region. If it **is** a material zone,
`ssr_exclude` names it directly and the mesh follows the zone boundary exactly. If it cuts across
the materials, it has to be a polygon. On [RS2-64b](../verification/rs2.md#rs2-64) the two paths
give identical factors of safety on the same mesh, and the conforming mesh's own answer sits about
6% below the non-conforming one — the difference is the discretization, not the masking.

### Tensile strength in the SSRM {#tensile-strength-in-ssrm}

Mohr-Coulomb is not a compression-only criterion. Extended into the tensile quadrant, the straight
envelope closes on an apex at

>>$\sigma'_t = -\dfrac{c}{\tan\phi}$

so a Mohr-Coulomb material given no further treatment carries an **implicit tensile strength** of
$c/\tan\phi$. For $\phi = 0$ — an undrained `cp` material, or `mc` entered with $\phi = 0$ — the
apex is at infinity and the implicit tensile strength is unbounded. That number is an artifact of
fitting a straight line to compression tests and extrapolating backwards, not a measured property:
real soil cracks at a small fraction of it, and often at zero.

![fem_ov_tension_cutoff.png](images/fem_ov_tension_cutoff.png){width=900}

In an ordinary stress analysis this rarely surfaces, because under the effective-stress formulation
a slope at working strength is in compression nearly everywhere. It surfaces in the **SSRM**, and it
surfaces asymmetrically: reducing strength divides $c$ and $\tan\phi$ by the same factor, so the
apex $c/\tan\phi$ is *invariant under reduction* — the tensile strength a material has at $F = 1$ it
still has at $F = 3$ (left panel above). Where the mechanism has to open a tension zone to develop —
a steep entry cut at a crest, a vertical face, the head of a scarp — that capacity acts as a
structural member holding the cut shut, and the model reaches *genuine* equilibrium at strength
reductions the real slope would never survive. Nothing in the failure criterion flags it: the states
are converged, force-balanced and budget-independent. The factor of safety simply comes out high.

**The cap.** The mat sheet's **t_cut** column sets a per-material tensile strength $T$, applied as a
**Rankine cutoff** $F_t = \sigma_1' - T$ that caps the major (most-tensile) principal effective
stress — a second viscoplastic yield surface driven by the same damped mechanism as the
Mohr-Coulomb surface, the two combined by Koiter's rule where both are active. It layers on top of
the shear envelope and never alters it. `t_cut = 0` means the material carries no tension at all; a
positive value caps the major principal stress there. The column is read automatically by
`solve_fem()` and `solve_ssrm()`; a script can override it per element with `tension_cap_by_elem`,
per material with `tension_cutoff_by_material`, or globally with the `tension_cutoff` flag, which is
simply the $T = 0$ case applied everywhere.

**Blank cutoff.** Where the cell is blank the material's own envelope decides how much tension it
carries, and XSLOPE enforces the envelope's limit: every Mohr-Coulomb element runs a Rankine cap at
its own apex $c/\tan\phi$, the most tension the Mohr-Coulomb envelope allows. The working cap is the smaller of the two, so a stated $T$ below the apex governs and the
apex governs where nothing is stated. It takes no state away from the shear envelope — every state
it can act on is one Mohr-Coulomb already forbids — and what it adds is a **return path**. The
$\psi = 0$ flow is purely deviatoric: it can shrink a stress circle but cannot move the circle's
center, so a Gauss point pulled to mean tension is inadmissible with no direction of flow that
returns it, and the iteration cannot settle. The Rankine flow is volumetric at the biaxial apex and
supplies exactly that return.

A **cohesionless** material carries no tension whether its cell is blank or
`0`: its apex sits at the origin, so the two entries describe the same admissible set and give the
same answer. And at $\phi = 0$ — a `cp` material, or `mc` entered with $\phi = 0$ — the apex is
at infinity, so only a stated $T$ bounds the tension there. Power-curve and Hoek-Brown elements are
left to their own envelopes, which carry a tensile strength of their own.

**Reducing the cap with $F$.** `tension_srf` decides whether the cap shrinks with the trial factor.
With `tension_srf=True` (**the default**) the solver divides it, $T_r = T/F$, exactly as it divides
$c$ and $\tan\phi$, so the reported factor of safety is the factor by which the *whole* envelope,
shear and tensile, is reduced (right panel above). With `tension_srf=False` the cap is held at its
authored value through the whole bisection. This is RS2's `tensilestrength_SRF` switch, and matching
it matters when the target is an RS2 answer: on a tension-controlled mechanism the two settings do
not converge to the same factor of safety. The default is on because it only ever acts *where a cap
exists and is positive* — a model with no `t_cut` and no global cutoff has no $T$ to reduce, so
every cap-less run (including all the Griffiths & Lane anchors) is identical either way, and a
cutoff of $T = 0$ is left where it is for the same reason, since $0/F$ is $0$ at every trial factor.
The apex cap is never divided either: $(c/F)/(\tan\phi/F) = c/\tan\phi$ is the same number at every
trial factor, which is the invariance the left panel above shows. The switch is reachable three
ways: the `tension_srf` keyword, the **Tension SRF** cell on the main sheet, and the matching
checkbox in Studio's Run FEM dialog, which is dimmed on both of those kinds of model.

**Which convention to run.** XSLOPE's default is *the envelope's own limit* — the Griffiths & Lane
convention, since those analyses state no separate tensile strength and a plain Mohr-Coulomb
material allows tension up to its apex and no further. Every
[Griffiths & Lane anchor](../verification/ssrm.md) in the verification suite is locked under
it. RS2 and Plaxis cap tension as a matter of course, writing an explicit per-material tensile
strength into the model (in Rocscience's own published verification models it is almost always
$T = c$, well below the apex). Neither convention is wrong, but they are not interchangeable, and
the difference is largest exactly where the mechanism is tension-controlled. **When comparing
against RS2 or Plaxis, set `t_cut` from the vendor model rather than leaving it blank**, and match
the vendor's tension-SRF switch. XSLOPE's RS2 reader does the first half automatically:
`xslope.rs2.read_fez` maps each material's tensile strength onto `t_cut`.

[RS2-62](../verification/rs2.md#rs2-62) — Cheng, Lansivaara & Wei's three-layer slope with a soft
band — is a benchmark where the tension cap controls the answer. Its vendor model caps the three materials at
$T$ = 20 / 0 / 10 kPa and reduces them with the SRF. Run uncapped, the cap soil's implicit
$c/\tan\phi \approx 28$ kPa holds the crest entry cut shut and the model equilibrates to
$F \ge 1.3$; run with the vendor caps and the tension SRF, the band mechanism mobilizes as limit
equilibrium predicts and the factor of safety is 0.769, against RS2's 0.81 and Plaxis' 0.82.

### Fast kernel

The cost of an SSRM run is dominated by the per–Gauss-point constitutive update, evaluated for every
Gauss point on every iteration of every trial. XSLOPE ships a **compiled kernel** that runs this
update in C (via Cython) instead of NumPy, which shortens a typical Mohr-Coulomb SSRM solve by
roughly a third to a half.

It is **used automatically when it is available**: the `fast_kernel` argument of `solve_fem()`
defaults to `"auto"`, meaning the compiled kernel if it is built and the NumPy path if it is not.
(`solve_ssrm()` has no such argument; its trials inherit the same automatic choice.) The installer
builds compile the kernel; a `pip install` gets a pure-Python wheel and therefore the NumPy path.
Both give the same answers, checked over all 103 FEM benchmarks, so which one a machine has
affects wall-clock time and nothing else. You never need the kernel to run XSLOPE.

The pure-NumPy path is the reference implementation: every factor of safety in the verification
suite is computed with it, and the compiled kernel reproduces it bit-for-bit. `fast_kernel=False`
forces the NumPy path. `fast_kernel=True` *requires* the kernel but warns
and falls back to NumPy if it has not been built, so the flag is always safe to set.

The kernel handles the standard Mohr-Coulomb path, including the Rankine tension cutoff, the
matric-suction term and the $K_0$ in-situ stress. Curved-envelope materials (power-curve and
Hoek-Brown) and all 1D reinforcement and pile work stay on the NumPy path automatically — a model
that mixes them accelerates its Mohr-Coulomb groups and leaves the rest unchanged.

To build it locally, with Cython installed:

```bash
pip install Cython
python setup_kernel.py build_ext --inplace
```

This compiles `xslope/_fem_kernel` next to its `.pyx` source; only the `.pyx` is tracked in the
repository. Once built, the default `"auto"` setting picks it up with no code change.
