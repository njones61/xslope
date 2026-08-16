---
title: "Tutorial LEM-11 — Reliability"
description: "A submerged clay slope in XSLOPE whose factor of safety is 1.354 and whose probability of failure is 17%: the standard-deviation columns, the Taylor series and Monte Carlo run on the same model, and what halving one σ buys."
---

# Tutorial LEM-11 — Reliability

Every tutorial before this one ends in a single number: enter the best
estimate of each soil property and a deterministic analysis returns one factor
of safety. A stochastic analysis starts from **more** data, not less. A site
investigation yields not just a best estimate of each property but a sense of
its spread — the scatter of the tests, the range two engineers would defend —
and a standard deviation entered beside the value puts that knowledge into the
model instead of leaving it in the report. The analysis carries it through the
computation, and the answer comes back correspondingly richer: a distribution
of factors of safety, summarized by a **probability of failure** and a
**reliability index** rather than by one number.
[Reliability Analysis](../reliability/index.md) gives the theory behind both.

A probability is also the more meaningful number to act on. In design and
retrofitting, "FS = 1.35" says nothing about how likely the slope is to fail
or what a proposed fix buys — one slope at 1.35 can be far riskier than
another at the same value — while "one chance in six, and the retrofit takes
it to one in a hundred" is a statement a decision can rest on, and one that
can be weighed against cost and consequence.

This page's slope — 30 ft of undrained clay standing under 40 ft of water —
makes the difference concrete. The deterministic answer is FS = 1.354; the
stochastic answer is that the slope has roughly one chance in six of failing.
**The factor of safety does not change; what changes is what is known about
it.**


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
<div class="tgm-model" markdown>**Completed model** — [xslope_reliability.xlsx](../lem/files/xslope_reliability.xlsx) — essentially the same problem as [LEM Sample Problem 15](../lem/samples.md#15-reliability-analysis-submerged-slope), with the standing water entered as a piezometric line instead of hand-typed surface loads</div>
</div>

---

## The slope

One soil: an undrained clay, γ = 120 pcf, c = 400 psf, φ = 0, over a rigid base
20 ft below the toe. The face rises 30 ft at 1.5:1, and the whole section is
under water — the free surface stands at elevation 40, so there is 40 ft of water
on the flat in front of the toe and 10 ft over the crest. The water surface is
entered once, as a flat piezometric line at elevation 40, and the model's
**Water loads** option is `auto`: at solve time the engine turns the standing
water above the ground into the surface loads it is. The clay's **u** stays
`none`, because a total-stress φ = 0 analysis carries no pore pressures on its
slice bases — here the line defines the reservoir, not a pore-pressure field.
[LEM-4](lem04_water_in_the_slope.md) covers the pore-pressure inputs this model
deliberately does without.

What makes the slope a reliability problem is the two extra pieces of data
beside those strengths: the site investigation supports a standard deviation of
8 pcf on the unit weight and 100 psf on the cohesion — 7% and 25% of their own
values — and the sheet carries them where a deterministic model would have
carried nothing.

### Opening the model

Download
[xslope_reliability.xlsx](../lem/files/xslope_reliability.xlsx)
and open it in Studio — **File → Open**. The Inputs plot draws the section, the
starting circle the file carries, the piezometric line at elevation 40 with the
water-table symbol on it, and the surface loads the water derives — drawn in
water blue and labelled as derived, because they are computed from the line
rather than entered:

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

The **Taylor Series Probability Method** (TSPM) is the workhorse of
geotechnical reliability, popularized in U.S. practice by Duncan: it solves
the model at the most-likely values and again with each uncertain parameter
moved one standard deviation up and one down — 1 + 2N solves, five here,
seconds of work. It goes first because it is the cheapest honest estimate:
from nothing but the σ columns it returns β, the probability of failure, and
a per-parameter breakdown the rest of this page is built on. It takes two
shortcuts to get there: it treats the factor of safety as changing in a
straight line with each parameter over that range, and it assumes the shape
of the FS distribution instead of observing it. Whether either shortcut
matters on this slope is exactly what the Monte Carlo section answers, by
solving the same problem with no shortcuts at all.

**Reliability…** sits beside **Parametric…** on the **Run** menu and on the
toolbar. Where a parametric sweep answers deterministic what-ifs, this one turns
the σ columns into a reliability index. Open it and leave the **Method** on
`Taylor series (TSPM)`, with **LEM method** = `Spencer` and 40 slices:

![The Reliability dialog on the Taylor series](images/lem11_studio_reliability_taylor.png)

The search is not done once. The Taylor series solves the model 1 + 2N
times — once at the most-likely values (MLV, the values entered on the mat
sheet), then twice more for each uncertain parameter, at MLV + σ and at
MLV − σ — and with **Search for the critical
surface at each solve (MLV and MLV ± σ)** ticked, that is exactly what
happens: **every one of those solves is a fresh search**, each perturbed
model finding its own critical surface. Leave it ticked here. Unticked,
every solve evaluates the first circle on the circles sheet instead — the
same surface a `Single surface` run solves — which is the right mode when
the surface itself is prescribed, as on a published benchmark's circle, an
observed failure surface, or a bedding-controlled plane, where the question
is that surface's own reliability. The checkbox restates
itself when the engine changes, because the sampling engines answer the
question differently — watch its label in the Monte Carlo step below. Below the controls, **Standard deviations in
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

Monte Carlo answers the same question by brute honesty instead of
approximation. One realization works like this:

1. **Draw** a value for each uncertain parameter from its probability
   distribution — here, a unit weight from the bell-shaped curve centered at
   120 with σ = 8, and a cohesion from the one centered at 400 with σ = 100.
2. **Build and solve.** That pair of draws is one complete version of the
   slope; its factor of safety is computed the ordinary way.
3. **Repeat** — ten thousand times — and keep every answer.

Because every value was drawn from its parameter's own probability
distribution — likely values drawn often, unlikely ones rarely — **each
realization is an equally probable version of the slope**. The objective is
to explore an exhaustive set of these versions and determine what portion of
them fail. The pile of ten thousand factors of safety *is* the answer: its
histogram is the distribution of FS, its mean and spread are read off
directly, and the probability of failure is simply the fraction of versions
that came back below FS = 1 — **counted, not inferred**.

The contrast with the Taylor series is the whole lesson. TSPM solved the
model five times at carefully chosen points (the MLV and ±σ probes), then
*assumed* a shape for the answer to turn those five numbers into a
probability — so its probability of failure is computed from an assumed
curve, not observed. Monte Carlo solves ten thousand versions and assumes
nothing: whatever the distribution of FS really is, it emerges. The price is the count,
and on a fixed surface where each solve costs milliseconds, the price is
affordable.

Open **Reliability…** again and change **Method** to `Monte Carlo`. The
controls below the surface options come live:

![The Reliability dialog on Monte Carlo](images/lem11_studio_reliability_mc.png)

**MC samples** `10000` is the number of realizations — how many versions of
the slope will be built and solved. Ten thousand is the default because the
probability of failure is a count, and a count is only as sharp as the pile
it comes from: with roughly one realization in six failing here, 10,000 of
them pin P<sub>f</sub> to well under a percentage point. The occasion to
change it is a *rarer* failure — a P<sub>f</sub> of a couple percent puts
only a couple hundred failures in the pile, and pinning it down as tightly
takes several times more realizations. Raising the count costs only time;
lowering it saves little and costs precision. The cleanest habit is to let the campaign size itself: tick
the convergence stop described below, and this field becomes the cap rather
than the count.

**MC seed** `20240117` starts the random-number generator. The value itself
means nothing — any integer serves — but keeping it fixed means the same
input file draws the same ten thousand realizations every time, so the
results repeat to the digit. There is no reason to change it in ordinary
work. The one good use for a different seed is a scatter check: rerun with
two or three other seeds, and P<sub>f</sub> should wobble only in its last
digits. If a conclusion flips when the seed changes, it was resting on too
few realizations — the fix is more samples, not a better seed.

**MC distribution** `Normal` sets the shape of the spread each σ describes. `Normal` draws
from the symmetric bell curve centered on the material's own value with the σ
column as its spread — the same reading the Taylor series gives those
columns. The alternative, `Lognormal`, is a skewed relative of the bell curve
that can never produce a negative value, matched to the same mean and σ; it
is worth considering when a σ is so large that the bell would spill past
zero. That is not the case here — the largest COV on this model is the
cohesion's 25% — and either way every draw is cut off at zero, so a sample
can never hand the solver a negative strength.

**MC sampling** chooses how the draws are laid down: `Random` is independent
draws; `Latin hypercube` cuts each distribution into equal-probability bins
with one draw per bin, which the
[Monte Carlo page](../reliability/monte_carlo.md) measures at roughly three
times the information per realization on this model.

**Stop when P_f converges** ends the campaign early once the probability of
failure is known to within a stated percentage of its own value — ±5% of a
17% answer means trusted between about 16% and 18%, which this model reaches
at 7,600 realizations — with the samples field as the cap; the
[Monte Carlo page](../reliability/monte_carlo.md) gives the rule and its
convergence plot. Leave it unticked here, and **Run**.

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
the two differ whenever F does not change in a straight line with a
parameter. The Taylor table above shows which parameter is responsible.
Average each parameter's own F⁺ and F⁻ and compare against F<sub>MLV</sub>:
if the response were a straight line, moving the parameter up and down by the
same amount would move F by equal and opposite amounts, and the average would
land back on F<sub>MLV</sub>. The cohesion's pair does exactly that — 1.354
to seven digits, because with φ = 0 the factor of safety really is a straight
line in c. The unit weight's pair averages **+0.0266 above**
F<sub>MLV</sub> — its response bends — and the Monte Carlo mean sits
**+0.0268 above** it. **The whole gap is the bend in the γ response**, which
the straight-line shortcut cannot see and a sampling method feels on every
draw.

---

## The response surface

Monte Carlo's answer still carries a sampling wobble. The 16.71% rests on a
count — 1,671 failing realizations — and rerunning with different seeds moves
that count; at 10,000 realizations, P<sub>f</sub> can land anywhere in a band
about ±0.7 percentage points wide (the same statistical band the convergence
stop watches). The only way to shrink it is more solves.

The third engine takes a different route to the same answer. Instead of
solving the real model for every realization, it solves it a **handful** of
times at chosen combinations of the uncertain parameters — the mean values,
plus points one and two standard deviations out in each direction — and fits
a smooth curved formula through those answers. That formula becomes a
stand-in for the solver: hand it any unit weight and cohesion and it predicts
the factor of safety instantly. The Monte Carlo campaign then runs on the
*formula* — ten million realizations, each one arithmetic rather than a
solve — and the counting wobble all but disappears. What is left is a new
question, whether the formula faithfully reproduces the real model, and the
engine spends most of its real solves answering exactly that.

Open **Reliability…** once more and set **Method** to `Response surface
(RS)`; the sample count and convergence stop gray out, because realizations
on the formula are nearly free:

![The Reliability dialog on the response surface](images/lem11_studio_reliability_rs.png)

Click **Run**. On this model the engine makes 9 design solves (two uncertain
parameters), then 700 more to *check itself* — 500 across the whole
distribution and 200 in the failure region — and reports in about four
seconds: **P<sub>f</sub> = 16.6%**, β<sub>LN</sub> = 0.998, with the fit's
credentials beside the answer: R² = 0.998, an rms fit error of 0.018 on the
factor of safety, and 1 of the 500 checked draws on the wrong side of
FS = 1. Paired against thirty thousand real solves of the identical draws,
the surrogate's P<sub>f</sub> is 0.10 percentage points low — the sampling
error is gone and the *fit* error is what remains, measured. The
[response-surface section](../reliability/monte_carlo.md#sampling-a-fitted-response-surface)
draws the fitted surface, the sample cloud and the F = 1 boundary in one
figure, and states the gate that makes the engine refuse a model it cannot
fit honestly.

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

All three engines read the same σ columns and none needs an input the others do
not, so the choice is about cost and about what is being asked.

| Engine | Real solves here | What it returns | What it assumes |
|---|:---:|---|---|
| Taylor series | 5 | β, lognormal P<sub>f</sub>, the per-parameter table | F changes in a straight line near the means |
| Monte Carlo | 2,000–10,000 | the whole distribution; P<sub>f</sub> *counted* | only the input distributions |
| Response surface | ~710 | MC-grade statistics with sampling error removed | F is quadratic — and it *checks* |

- **The Taylor series** is the screen: five solves, seconds, and the engine
  behind both the variance Pareto and the finite-element reliability path. It
  summarizes F by how it changes near the mean values, and everything above
  shows what that discards (the bend in the γ response) and what it keeps
  (the ranking, the variance, β itself).
- **Monte Carlo** is the reference: nothing assumed about the shape of F, the
  tail counted rather than fitted. The convergence stop sets its cost to the
  resolution actually asked, and Latin hypercube stretches each solve about
  three-fold. It is the engine to reach for when the straight-line shortcut
  breaks — a coefficient of variation so large that x − σ is not a physical
  value ([VP34](../verification/rocscience.md#vp34), a fill whose friction
  angle carries a 124% COV), a response that bends strongly, or a wanted
  empirical tail — and it is the arbiter when the other two disagree.
- **The response surface** is the precision instrument: tail resolution no
  count can afford, for a few hundred real solves — *when its gate accepts*.
  It refuses exactly the models Monte Carlo exists for (VP34's gate finds a
  third of the surrogate's failures have no real solution, and refuses), and
  its β differs from Monte Carlo's in the third decimal because a quadratic
  slightly under-disperses σ<sub>F</sub>. Quote its P<sub>f</sub> with the
  gate credentials it prints beside it.

A working habit: screen with the Taylor series, decide with Monte Carlo at the
default convergence stop, and reach for the response surface when the tail
itself is the question — checking it against a Monte Carlo run once, since the
two sampling engines share their draws by construction.

No engine randomizes the failure surface — the slip surface is a decision, not
a random variable. The Taylor series searches at every one of its 1 + 2N
solves, each perturbation finding its own critical surface; the two sampling
engines find the surface once, at the most-likely values, and hold it across
their realizations, which is what the commercial limit-equilibrium codes do in
their probabilistic modes.

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
- Where the two means differ: the +0.027 gap comes from the bend in the γ
  response, matched to the digit by the average of the Taylor table's own γ
  perturbations.
- The variance Pareto measuring **75.8%** of σ<sub>F</sub>² onto the cohesion and
  **24.2%** onto the unit weight, and what halving the dominant σ does — P<sub>f</sub>
  from **17.48% to 6.36%** (TSPM) and **16.71% to 5.34%** (Monte Carlo) with the
  factor of safety fixed at 1.354.
- The response surface reaching the same answer from **709 real solves** —
  P<sub>f</sub> = 16.6% with its fit credentials printed beside it — and the
  three-engine habit: screen with the Taylor series, decide with Monte Carlo,
  bring in the response surface when the tail is the question.

**Where to go next:** the [tutorials index](index.md) lists the series.
[Reliability Analysis](../reliability/index.md) gives the theory and the equation
each index is built from, the [Taylor Series](../reliability/taylor.md) and
[Monte Carlo](../reliability/monte_carlo.md) pages document the two engines, and
[Reliability Analysis (FEM)](../reliability/fem.md) runs the same Taylor series
with a strength-reduction factor of safety in place of a limit-equilibrium
search.
