---
title: "Tutorial LEM-11 — Reliability"
description: "A submerged clay slope in XSLOPE whose factor of safety is 1.354 and whose probability of failure is 17%: the standard-deviation columns, the Taylor series and Monte Carlo run on the same model, and what halving one σ buys."
---

# Tutorial LEM-11 — Reliability

A 30 ft slope in undrained clay, standing under 40 ft of water, with a factor of
safety of 1.354. The clay's strength and unit weight were not measured exactly —
they carry standard deviations on the materials sheet — and running those
standard deviations through the analysis says the slope has roughly one chance in
six of failing. **The factor of safety does not change; what changes is what is
known about it.**


<div class="tut-glance" markdown>
<div class="tgt-row">
<div class="tgt-tile"><span class="tg-label">Analysis</span><p>Limit equilibrium</p></div>
<div class="tgt-tile"><span class="tg-label">Open &amp; run</span><p>~15 min</p></div>
</div>
<div class="tgm-obj" markdown>
**Objectives** — Read the standard-deviation columns a reliability analysis runs
on, estimate the reliability index and probability of failure two ways — the
Taylor series over 1 + 2N searches and a 10,000-realization Monte Carlo campaign
— on one model and one surface, measure which of the two uncertainties actually
drives the answer, and halve it.
</div>
<p><span class="tg-pill">one material</span><span class="tg-pill">undrained strength</span><span class="tg-pill">distributed load</span><span class="tg-pill">standard deviations</span><span class="tg-pill">Taylor series (TSPM)</span><span class="tg-pill">Monte Carlo</span><span class="tg-pill">variance Pareto</span><span class="tg-pill">circular search</span></p>
<div class="tgm-model" markdown>**Completed model** — [xslope_prob_submerged_KEY.xlsx](../lem/files/xslope_prob_submerged_KEY.xlsx), the same file used by [LEM Sample Problem 15](../lem/samples.md#15-reliability-analysis-submerged-slope)</div>
</div>

---

## The slope

One soil: an undrained clay, γ = 120 pcf, c = 400 psf, φ = 0, over a rigid base
20 ft below the toe. The face rises 30 ft at 1.5:1, and the whole section is
under water — the free surface stands at elevation 40, so there is 40 ft of water
on the flat in front of the toe and 10 ft over the crest. That water is modelled
as a distributed load on the ground surface rather than as a pore-pressure model,
which is why the clay's **u** is `none`: a total-stress φ = 0 analysis carries no
pore pressures on its slice bases, and the water's weight arrives as the load it
is. [LEM-4](lem04_water_in_the_slope.md) covers the pore-pressure inputs this
model deliberately does without.

What makes the slope a reliability problem is the two extra numbers beside those
strengths. The unit weight carries a standard deviation of 8 pcf and the cohesion
one of 100 psf — 7% and 25% of their own values.

### Opening the model

Download
[xslope_prob_submerged_KEY.xlsx](../lem/files/xslope_prob_submerged_KEY.xlsx)
and open it in Studio — **File → Open**. The Inputs plot draws the section, the
starting circle the file carries, and the water as two blocks of load arrows:

![The loaded model](images/lem11_inputs.png){width=1000}

The arrows are longest on the flat at elevation 0, where the water is 40 ft deep
and presses 2496 psf, and shortest over the crest at elevation 30, where 10 ft of
water presses 624 psf. Along the face between them the load tapers linearly with
the depth of water above each point.

### The standard deviations

Open **Materials** and switch to **List view**. Above the form, the row of
checkboxes labelled **Show parameters for:** decides which inputs the editor
shows; **Reliability** is the one that is off by default, and ticking it puts a
**± σ** box beside every parameter that can carry a standard deviation:

![The clay with its standard deviations](images/lem11_studio_materials.png)

The clay reads **γ** = `120` **± σ** `8` and **c** = `400` **± σ** `100`. Those
are the mat sheet's `s(g)` and `s(c)` columns, which the table view shows as two
columns in a **Standard Deviations** block at the right of the sheet and which
the list view puts next to the values they belong to. A blank or zero σ means the
parameter is treated as known exactly and is not varied.

A standard deviation is easier to judge as a fraction of its own value — the
**coefficient of variation**, σ/x — which is 8/120 = 6.7% for the unit weight and
100/400 = 25% for the cohesion. Both sit inside the ranges published for those
properties: 3–10% for unit weight, 15–50% for undrained shear strength. The
[Reliability Analysis](../reliability/index.md#parameter-uncertainty) page
tabulates the rest and gives the two rules of thumb for estimating a σ from a
high and a low estimate rather than from a data set.

---

## The deterministic run

Click **Run LEM…** and choose **Method** = `Spencer` and **Analysis** =
`Auto search`, with the slice count left at 40:

![The Run LEM dialog on the loaded model](images/lem11_studio_run_lem.png)

Click **Run**. The search deepens the file's circle onto the base of the model:

![Spencer on the most-likely values](images/lem11_solution.png){width=1000}

**FS = 1.354**, on a circle centered at (23.14, 47.36), tangent to the rigid base
at elevation −20, 141.27 ft of failure surface carrying 370,624 lb/ft of clay.
This is the answer every deterministic analysis in the previous ten tutorials
would stop at, and it is the number the rest of this page is about: not whether
1.354 is right, but how much of it is real.

---

## The Taylor series

**Reliability…** sits beside **Parametric…** on the **Run** menu and on the
toolbar. Where a parametric sweep answers deterministic what-ifs, this one turns
the σ columns into a reliability index. Open it and leave the **Method** on
`Taylor series (TSPM)`, with **LEM method** = `Spencer` and 40 slices:

![The Reliability dialog on the Taylor series](images/lem11_studio_reliability_taylor.png)

**Search the critical surface at the mean values** is ticked by default and
should stay ticked: the surface is found once, at the most-likely values, and the
uncertainty is then evaluated on it. Below the controls, **Standard deviations in
this file** lists what the run will actually vary — `mat:soil:c = 400 ± 100
(COV 25%)` and `mat:soil:gamma = 120 ± 8 (COV 7%)`. If that box says no standard
deviations are set, **Run** is disabled; a reliability analysis with nothing
uncertain in it has no question to answer.

Click **Run**. The result view carries the reliability statistics in its title
over the surface they were computed on:

![The Taylor-series result](images/lem11_taylor.png){width=1000}

The Taylor Series Probability Method evaluates the factor of safety at the
most-likely values and then once more at ±σ for each uncertain parameter — 1 + 2N
solves, five of them here — and combines the swings into a variance. The
per-parameter table is written to the Log pane:

| Parameter | MLV | σ | F⁺ | F⁻ | ΔF |
|---|:---:|:---:|:---:|:---:|:---:|
| soil γ | 120 | 8 | 1.189 | 1.572 | 0.383 |
| soil c | 400 | 100 | 1.692 | 1.015 | 0.677 |

The unit-weight row runs the other way from the strength row: heavier clay is a
larger driving moment, so **F⁺ is the lower** of its pair. Only the size of
each swing enters the variance, σ<sub>F</sub>² = Σ(ΔF/2)², which here gives
**σ<sub>F</sub> = 0.389** and **COV<sub>F</sub> = 0.287**. From the factor of
safety and its coefficient of variation comes the lognormal reliability index
**β<sub>LN</sub> = 0.935**, a **reliability of 82.52%** and a **probability of
failure of 17.48%**.

The plot's legend names F⁺ and F⁻ surfaces that cannot be found on it, and that
is the result rather than a drawing fault: all five searches returned the same
circle, center (23.14, 47.36) tangent at elevation −20, so the perturbation
surfaces lie exactly under the F<sub>MLV</sub> surface drawn over them. On a
single-material φ = 0 slope neither scaling the strength nor scaling the weight
moves the critical circle — both scale the resisting and driving terms of every
trial surface the same way. On a layered slope the perturbed surfaces separate,
and where they separate a long way the parameter is moving the mechanism and not
just the number.

---

## Monte Carlo

Open **Reliability…** again and change **Method** to `Monte Carlo`. The three
controls below the surface options come live:

![The Reliability dialog on Monte Carlo](images/lem11_studio_reliability_mc.png)

**MC samples** `10000` is the number of realizations, **MC seed** `20240117` is
the random seed — fixed rather than taken from the clock, so the run reproduces
exactly — and **MC distribution** `Normal` is how each σ is read: as the standard
deviation of a normal distribution about the material's own value, truncated at
zero so a draw can never hand the solver a negative strength. Leave all three and
**Run**.

Ten thousand solves take appreciably longer than five. What comes back is not an
estimate of the distribution but the distribution itself:

![The Monte Carlo FS histogram](images/lem11_mc.png){width=1000}

**Mean FS = 1.381, σ<sub>F</sub> = 0.408**, and the shaded tail left of the
FS = 1 line holds **1,671 of the 10,000 realizations — P<sub>f</sub> = 16.71%**,
counted rather than inferred. The reliability index comes in two conventions:
β<sub>normal</sub> = (F̄ − 1)/σ<sub>F</sub> = **0.932**, and the lognormal
β<sub>LN</sub> = **0.969** computed from the same sample moments as the Taylor
series uses.

Against 17.48% from five solves, 16.71% from ten thousand. The two estimators
agree to eight tenths of a percentage point on the answer that matters, which is
the ordinary result — the [reliability
overview](../reliability/index.md#when-to-use-monte-carlo-versus-the-taylor-series)
sets out the cases where they part company, and
[VP29](../verification/rocscience.md#vp29) reaches the same conclusion on
Duncan's LASH terminal, where the choice of σ inputs moves the probability of
failure further than the choice of estimator does.

The one visible disagreement is the mean: the Taylor series reports
F<sub>MLV</sub> = 1.354 and Monte Carlo a mean of 1.381, 2% higher. Both are
right about different things. F<sub>MLV</sub> is the factor of safety *at* the
mean inputs; the Monte Carlo mean is the *mean of* the factors of safety, and
those differ whenever F curves. The Taylor table above measures the curvature
directly, one parameter at a time, by averaging each parameter's own F⁺ and F⁻
against F<sub>MLV</sub>. The cohesion's pair averages to 1.354, returning
F<sub>MLV</sub> to seven digits, because with φ = 0 the factor of safety is
exactly linear in c and a symmetric perturbation of it cancels. The unit weight's
pair averages **+0.0266 above** F<sub>MLV</sub>, and the Monte Carlo mean sits
**+0.0268 above** it. **The whole gap is the curvature in γ**, which a
first-order method drops by construction and a sampling method does not.

---

## Reading β against P<sub>f</sub>

The reliability index and the probability of failure are two spellings of one
result. β is the distance from the mean factor of safety to failure, measured in
standard deviations of the factor of safety — on the lognormal scale the
[reliability overview](../reliability/index.md#reliability-equation) gives — and
P<sub>f</sub> is the area of that distribution below FS = 1, so
P<sub>f</sub> = Φ(−β) and the two always move together:

| Estimator | F | σ<sub>F</sub> | β | P<sub>f</sub> |
|---|:---:|:---:|:---:|:---:|
| Taylor series (TSPM) | 1.354 | 0.389 | 0.935 (lognormal) | 17.48% |
| Monte Carlo, 10,000 samples | 1.381 | 0.408 | 0.969 (lognormal) | 16.71% (counted) |

β is the more stable of the two to quote, because P<sub>f</sub> is exponentially
sensitive to it: β = 0.935 is 17.5%, β = 2 is 2.3%, β = 3 is 0.13%. The
consequence for a sampled estimate is a resolution limit. Each of the 10,000
realizations is worth 0.01% of P<sub>f</sub>, so the 16.71% above rests on 1,671
of them and is solid, while a P<sub>f</sub> near 0.05% would rest on five and
would move by a fifth of itself if one realization landed differently. The same
limit is visible from the other side on
[VP28](../verification/rocscience.md#vp28), where SLOPE/W's own 10⁴-sample
campaign reports a 0.04% probability of failure that is four realizations.
Resolving a small P<sub>f</sub> by counting means paying for the samples; the
Taylor series gets there by extrapolating a fitted tail from β, which costs
nothing and assumes the tail's shape.

---

## Which uncertainty is worth reducing

σ<sub>F</sub> = 0.389 is built from two contributions, and they are not equal.
The **Parametric…** dialog computes the split: set **Method** to `Spencer` and
**Plot type** to `Variance Pareto (σ)`, which is offered only for a model that
carries standard deviations. The sweep table below it stays empty — this plot
reads every σ-carrying material and ignores the table, as the note under it says:

![The Parametric dialog set to the variance Pareto](images/lem11_studio_parametric_variance.png)

Click **Run**. The bars are the two parameters, tallest first, with the running
total over them:

![Each parameter's share of the factor-of-safety variance](images/lem11_variance.png){width=800}

Each bar is that parameter's own (ΔF/2)² term as a share of σ<sub>F</sub>², which
is the same arithmetic the Taylor series already did — the cohesion carries
**75.8%** of the variance and the unit weight **24.2%**. Cohesion dominates on
the width of its uncertainty and not on the slope's sensitivity to it. Divide
each ΔF by the coefficient of variation that produced it and the unit weight is
the more powerful parameter: 0.383 over 6.7% is 0.057 of factor of safety per 1%
of γ, against 0.677 over 25% for 0.027 per 1% of c. Per unit of change γ moves
the answer twice as far — but c is asked over an interval nearly four times
wider, and the variance is what the interval is squared into.

Three quarters of the uncertainty in the answer is therefore in one number, and
that number is the one further site investigation would narrow. Another round of
undrained testing on this clay would buy a smaller s(c) and nothing else.

### Halve it

Open **Materials**, select the clay, and change the **± σ** beside **c** from
`100` to `50`. In the table view — and on the mat sheet — this is the second cell
of the two-column **Standard Deviations** block:

| s(g) | s(c) |
|:---:|:---:|
| 8 | 50 |

Nothing else on the row moves: **γ** is still `120`, its σ still `8`, and the
cohesion itself is still `400`:

![The clay with its cohesion σ halved](images/lem11_studio_materials_edit.png)

Click **OK** and run **Reliability…** again on the Taylor series, then a second
time on Monte Carlo, with everything else as before:

![The Monte Carlo histogram after the edit](images/lem11_mc_tightened.png){width=1000}

The histogram is the same shape drawn on a narrower base. Beside the run before
it:

| Run | F | σ<sub>F</sub> | β<sub>LN</sub> | P<sub>f</sub> |
|---|:---:|:---:|:---:|:---:|
| TSPM, s(c) = 100 | 1.354 | 0.389 | 0.935 | 17.48% |
| TSPM, s(c) = 50 | 1.354 | 0.256 | 1.526 | 6.36% |
| Monte Carlo, s(c) = 100 | 1.381 | 0.408 | 0.969 | 16.71% |
| Monte Carlo, s(c) = 50 | 1.381 | 0.270 | 1.567 | 5.34% |

**The factor of safety is unchanged at 1.354, to every digit.** Nothing about the
slope moved — the same clay at the same strength on the same critical circle. The
probability of failure fell from roughly one in six to roughly one in sixteen,
11.1 percentage points, and the failed count in the Monte Carlo campaign fell
from 1,671 realizations to 534. What was bought was not strength but knowledge of
it, and only a probabilistic analysis can be shown that purchase at all.

The unit weight's contribution is what remains. Its ΔF is 0.383 before the edit
and 0.383 after — untouched, because nothing about γ changed — and it now carries
56% of a smaller variance instead of 24% of a larger one. A second halving of
s(c) would buy correspondingly less; the parameter to attack next is the one the
Pareto puts first *after* the edit, which is why the measurement is worth
repeating rather than doing once.

---

## When to use which

Both estimators read the same two columns and neither needs an input the other
does not, so the choice is about cost and about what is being asked.

- **The Taylor series** costs 1 + 2N solves — five here, seconds — and returns β
  and a lognormal P<sub>f</sub>. It is a first-order method: it summarizes F by
  its slope at the mean values, and everything above shows what that discards
  (the curvature in γ) and what it does not (the ranking of the parameters, the
  variance, β itself). It is the default, and the engine behind both the variance
  Pareto and the finite-element reliability path.
- **Monte Carlo** costs 10⁴ solves and returns the whole distribution, with
  P<sub>f</sub> counted from the tail rather than fitted to it. It earns that
  cost where the first-order assumptions break: a coefficient of variation so
  large that x − σ is not a physical value and the Taylor series declines the
  analysis outright, a strongly non-linear factor-of-safety response, or a
  wanted empirical tail. [VP34](../verification/rocscience.md#vp34) is the
  worked case of the first — a fill whose friction angle carries a 124%
  coefficient of variation.

Neither randomizes the failure surface. The slip surface is a decision, not a
random variable: it is found once at the most-likely values and then held while
the parameters vary, which is what both engines above did and what the commercial
limit-equilibrium codes do in their probabilistic modes.

---

## Conclusion

This tutorial demonstrated:

- The **± σ** boxes the materials editor's **Reliability** toggle reveals, on
  this clay's `s(g)` = 8 and `s(c)` = 100 — coefficients of variation of 6.7% and
  25%.
- A deterministic Spencer search at **FS = 1.354** on a circle tangent to the
  rigid base, and the Taylor series over that same search: **σ<sub>F</sub> =
  0.389**, **β<sub>LN</sub> = 0.935**, **P<sub>f</sub> = 17.48%** from 1 + 2N =
  five solves.
- Monte Carlo on the same surface: **mean FS = 1.381**, **σ<sub>F</sub> =
  0.408**, **β<sub>LN</sub> = 0.969**, and **P<sub>f</sub> = 16.71%** counted as
  1,671 of 10,000 realizations — within a percentage point of the Taylor series.
- Where the two means differ: the +0.027 gap is the curvature in γ, matched to
  the digit by the average of the Taylor table's own γ perturbations.
- The variance Pareto measuring **75.8%** of σ<sub>F</sub>² onto the cohesion and
  **24.2%** onto the unit weight, and what halving the dominant σ does — P<sub>f</sub>
  from **17.48% to 6.36%** (TSPM) and **16.71% to 5.34%** (Monte Carlo) with the
  factor of safety fixed at 1.354.

**Where to go next:** the [tutorials index](index.md) lists the series.
[Reliability Analysis](../reliability/index.md) gives the theory and the equation
each index is built from, the [Taylor Series](../reliability/taylor.md) and
[Monte Carlo](../reliability/monte_carlo.md) pages document the two engines, and
[Reliability Analysis (FEM)](../reliability/fem.md) runs the same Taylor series
with a strength-reduction factor of safety in place of a limit-equilibrium
search.
