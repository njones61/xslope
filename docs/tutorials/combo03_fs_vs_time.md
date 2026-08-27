---
title: "Tutorial COMBO-3 — Factor of Safety vs Time"
description: "Two dams run for stability against every saved frame of a transient drawdown — SEEP-3's earth dam at drained strengths, where the lowest factor of safety falls 12 days before the reservoir stops moving and the critical face changes from downstream to upstream and back, and COMBO-2's Johnson Reservoir dam as a three-stage rapid drawdown at each instant, where the curve passes through COMBO-2's own answer at day 50 and stays below the drained curve for the rest of the run."
---

# Tutorial COMBO-3 — Factor of Safety vs Time

A steady seepage analysis gives a stability run one pore-pressure field, so the
run has one answer. A transient analysis gives it a *sequence* of fields, one per
saved instant, and the stability question can be asked of each of them in turn.
The result comes back as a curve — a factor of safety against time — and it
answers something a single analysis cannot: **when** the slope is at its weakest,
not only how weak.

[Tutorial SEEP-3](seep03_reservoir_drawdown.md) built a small earth dam with a
granular shell and a clay core and lowered its reservoir 16 m over 45 days. That
page stopped at the pore pressures. This one carries them into a limit
equilibrium analysis of the dam and runs it against every frame the march saved.
The analysis is offered both slopes, because over a run this long neither one
governs throughout: a drawdown attacks the upstream face, since lowering the pool
takes the water off it and leaves the water inside the embankment, but a full
reservoir loads that same face and leaves the downstream slope the weaker of the
two.

The one thing SEEP-3's dam does not carry is strength, because a seepage analysis
reads none. Part 1 assigns it, and everything else about the model comes
across unchanged.

A curve like that is only as good as the strengths under it, and drained
strengths are not the right ones for every dam. **Part 1** asks the stability
question of a dam whose zones both carry drained effective-stress envelopes, so
each instant is an ordinary effective-stress analysis of the water as it stands
then. **Part 2** asks the same question of
[COMBO-2](combo02_rapid_drawdown.md)'s Johnson Reservoir dam, whose compacted
clay core carries an undrained envelope, so each instant is a full three-stage
rapid drawdown instead — and the two curves that come out of that one march are
what a designer compares.

<div class="tut-glance" markdown>
<div class="tgt-row">
<div class="tgt-tile"><span class="tg-label">Analysis</span><p>Seepage + limit equilibrium (drained, and rapid drawdown)</p></div>
<div class="tgt-tile"><span class="tg-label">Open &amp; run</span><p>~45 min</p></div>
</div>
<div class="tgm-obj" markdown>
**Objectives** — Learn what a factor-of-safety-versus-time curve is and what it
adds to a single analysis; how a transient seepage field reaches a stability run
one instant at a time, and where that instant is chosen; how to place one
starting circle per face so a curve that spans a full reservoir and an empty one
can report either, and how to read which face governs at each instant; how the
saved times set the curve's resolution and how to choose them; how to read the
curve's minimum against the drawdown that caused it; and how to run a three-stage
rapid drawdown at every instant of a march and compare that curve against the
drained one.
</div>
<p><span class="tg-pill">two materials</span><span class="tg-pill">three materials</span><span class="tg-pill">transient seepage</span><span class="tg-pill">u = seep</span><span class="tg-pill">saved frames</span><span class="tg-pill">seepage time</span><span class="tg-pill">FS vs time</span><span class="tg-pill">rapid drawdown</span><span class="tg-pill">three-stage procedure</span><span class="tg-pill">Kc = 1 envelope</span><span class="tg-pill">d and ψ</span><span class="tg-pill">stage times</span><span class="tg-pill">Spencer</span><span class="tg-pill">circular search</span><span class="tg-pill">starting circles</span><span class="tg-pill">minimum slip depth</span><span class="tg-pill">automatic water loads</span></p>
<div class="tgm-model" markdown>
**Part 1 model** — [xslope_earth_dam_fs_time.xlsx](files/xslope_earth_dam_fs_time.xlsx),
SEEP-3's completed dam with a strength band on its materials table, a starting
circle on each face and a minimum slip depth. It carries no mesh and no solution,
so every run below is made from scratch

**Part 2 model** — [xslope_johnson_fs_time.xlsx](files/xslope_johnson_fs_time.xlsx),
COMBO-2's completed Johnson Reservoir dam: an undrained envelope on the core,
`u = seep` on all three zones, the pool schedule and both stage times. This copy
ships meshed and marched, so it opens with the twenty-one pore-pressure fields
Part 2 sweeps already on it; [COMBO-2](combo02_rapid_drawdown.md) builds both
</div>
</div>

---

## Part 1 — Drained strengths at every instant

SEEP-3's dam is the simpler of the two: both its zones carry drained
effective-stress envelopes, so every point of the curve below is an ordinary
effective-stress analysis run against the pore pressures the march computed for
that moment.

### The dam and the drawdown

![The dam, its two zones, and the properties that make them behave differently](images/seep03_problem_sketch.png){width=1000}

The section is 110 m long and 22 m tall at the crest, sitting on rock at
elevation 0, with both faces sloping at about 2.3:1. A clay **core** runs the
full height from the rock up to 4 m below the crest, about 17 m wide at its base
and 8 m wide at its top, and the granular **shell** is everything else. The
reservoir stands at elevation 18 against the upstream face, the tailwater at
elevation 2 at the downstream toe.

The pool is held at 18 for 2 days, drawn down 16 m to elevation 2 over the next
45 days, and held there for the remaining 253 days of a 300-day run. The two
zones drain at very different rates — the shell follows the pool down, the core
does not — and the head the core keeps after the pool has gone is the reason a
drawdown is a stability problem at all.
[SEEP-3](seep03_reservoir_drawdown.md) builds that model from nothing: the
storage properties, the boundary set and the `pool` time series. The section, the
zones and the schedule the reservoir follows are all unchanged here. What this
file does change is how often the march stops to save a frame, and the section
below on the march says why.

#### What the curve is made of

Each saved frame carries a pore pressure at every node of the mesh. A limit
equilibrium run with `u = seep` reads one such field at every slice base, so one
frame produces one factor of safety. Nineteen frames produce nineteen, and
plotting them against the frame times draws the curve. No input changes between the
points — the same section, the same strengths, the same solver at every one — so
everything the curve does is the pore pressures moving and the reservoir load
leaving the face.

Two consequences follow, and both matter further down. **The curve can only be
as fine as the save schedule**, because a time that is not a saved frame has no
field to read and is never interpolated between two that are. And **the critical
surface moves**: the mechanism that governs a dam under a full reservoir is not
the one that governs it half drained, and on this dam it is not even on the same
slope, so each instant gets its own search over the whole section.

### The strengths this page adds

The first thing to see is what the file carries beyond SEEP-3's seepage
properties. Open [xslope_earth_dam_fs_time.xlsx](files/xslope_earth_dam_fs_time.xlsx) with
**File → Open…**. The toolbar's mode strip reads LEM | Seepage | FEM; leave it on
**LEM** for now and click **Materials** in the Inputs tree. Set the **Show
parameters for:** toggles to **LEM** alone:

![The two zones with their strength band filled](images/combo03_studio_materials.png)

| mat | name | γ (kN/m³) | γsat (kN/m³) | option | c′ (kPa) | φ′ (deg) | u |
|:---:|---|:---:|:---:|:---:|:---:|:---:|:---:|
| 1 | `shell` | 20 | 21 | `mc` | 0 | 32 | `seep` |
| 2 | `core` | 19 | 20 | `mc` | 10 | 25 | `seep` |

**These are typical values for a granular shell and a compacted clay core, chosen
for the exercise rather than measured on this dam.** SEEP-3's file carries none
at all, and nothing below depends on their exact size — what the page measures is
how one model's factor of safety moves as the water in it moves.

**Two unit weights per row.** γ applies above the water table and γsat below it,
and the slicer splits each slice's weight where the two meet. On a drawdown that
is not bookkeeping: the shell the pool drains gets lighter as it drains, and the
driving weight falls with the resistance.

**Both strengths are drained.** c′ and φ′ are effective-stress parameters, so
each slice base needs a pore pressure to compute an effective normal stress from,
and that is what the last column supplies.

**`u` is `seep` on both rows.** That column decides where a slice base gets its
pore pressure, and `seep` sends it to a solved seepage field —
[COMBO-1](combo01_seepage_stability.md#the-column-that-connects-the-modes) covers
its four values and what each one costs. With a transient solution loaded, `seep`
means one frame of that solution, and which frame is a choice made at the run
rather than here.

Click **OK**.

### Two starting circles, one per face

The strengths say how strong the slope is; the circles say where the search
starts looking. Click **Circles**. The file carries two, one on each face of the
dam:

| Xo | Yo | Option | Depth |
|:---:|:---:|:---:|:---:|
| 7 | 56 | `Depth` | 0 |
| 103 | 59 | `Depth` | 0 |

Each one is drawn on the deep mechanism its face can make. The center sits out
beyond the heel or beyond the toe, **Depth** puts the bottom of the circle at
elevation 0 — an elevation rather than a distance, and 0 is the rock the dam
stands on — and the radius carries the trace across the whole slope. The
upstream circle runs from (0.80, 0.34) at the heel up to the crest edge at
(51.50, 22.00), 57.63 m of surface; the downstream one runs from the crest edge
at (57.04, 22.00) down to (109.19, 0.33) at the toe, 58.88 m.

**Two circles because the critical face changes with the pool.** A drawdown
attacks the upstream slope, and that is the mechanism this page is mostly about.
But a full reservoir *loads* the upstream slope, and under 18 m of water the
weaker side of this dam is the downstream one. A curve that runs from a full pool
to an empty one crosses both states, so both faces have to be reachable.

**Both are drawn deep because a shallower seed refines into a smaller
mechanism.** A search steps away from its starting circle and stops where the
steps no longer lower the factor of safety, so a circle drawn part way up the
slope tends to settle near where it started rather than travel out to the deep
surface.

Neither circle cuts the core, which follows from this dam's proportions rather
than from any rule. The core is 8 m wide across its top and its crown stands 4 m
below the crest, so a circle entering at the crest edge and running to the rock
passes over its shoulder and spends its length in the shell. On a dam whose core
is wider, or carried closer to the crest, the critical circle usually cuts it,
and the starting circle should be drawn to cut it too.

A circle sketched through the core still reaches a shell slide from a single
face: one drawn at (22, 52) with a radius of 51 m, crossing 14 m of core, refines
on day 35 to 1.338 at (9.12, 46.76), against 1.331 for the circle in the table.
The search walks it out of the core and down to the rock.

Under the table, the **Search window** group holds the limits that confine where
a searched surface may run. There are ten, and a blank field is a limit that is
not applied;
[the template page](../usage/input_template.md#search-window-optional) lists all
ten. One is set here:

| Limit | Value | What it does |
|---|:---:|---|
| Min slip depth | 8 | The surface must reach at least 8 m below the ground surface |

It sets a floor on the *size* of the mechanism, and the shell is why the model
takes one. A cohesionless soil has no length scale in it: with c′ = 0 the factor
of safety on a shallow surface parallel to the face does not depend on how deep
that surface is drawn, so an unconstrained search walks down to a surficial wedge
whose failure is raveling rather than a slide. 8 m on a 22 m dam keeps the search
on slides of the embankment rather than on crest slivers. This is Slide2's
minimum-depth filter, and the Run LEM dialog offers the same limit as **Ignore
surficial (skin) failures** for a run that wants to set it rather than read it
off the file.

The entry and exit ranges are deliberately blank. Filling them would confine the
trace to one slope, and which slope governs at a given instant is the thing this
page measures.

![Two starting circles and the one search limit the file sets](images/combo03_studio_circles.png)

Two circles with **Grid search** off raises `circles.multiple_without_grid` in
the **Model checks** panel on every run below, because a seeded search screens
each starting circle on one coarse grid and then refines only the best-scoring
one. That warning is expected here: both circles were placed deliberately, on the
mechanism each face actually makes, so the screen compares two real surfaces and
the refinement follows whichever face is weaker at that instant. Where the
critical face is not known in advance, ticking **Grid search (auto-seed the
circular search)** sweeps a grid of centers against a range of tangent elevations
derived from the geometry and refines every competing family instead.

Click **OK** on the Circles editor. The Inputs plot draws the model the runs
below are made on:

![The stability model: two zones, the reservoir load, and the two starting circles](images/combo03_inputs.png){width=1000}

The blue arrows on the upstream face are the reservoir pressing on the slope,
marked **derived** because nothing entered them: **Water loads** on the `main`
sheet is `auto`, so the engine measures the pool against the ground and turns
what stands above it into a distributed load. The drawdown takes that load away,
and it is re-derived at every instant of the curve from the pool as it stood at
*that* moment. The two red dashed arcs are the starting circles, one under each
face; both centers sit above the top of the frame.

### The mesh and the initial condition

Every pore pressure on this page comes from the seepage engine, so the mesh and
the full-pool steady solution come before any stability run. Switch the mode
strip to **Seepage** (`Ctrl+2`) and click **Run → Build Mesh…** Set **Element type** to **Linear triangles (tri3)**, leave **Auto-size from
geometry** ticked, and set **Size divisions** to `64`. Click **Build**: **614
nodes and 1,089 triangles**, the mesh
[SEEP-3 builds](seep03_reservoir_drawdown.md#building-the-mesh) and every number
on this page is computed on. Head is a scalar field and no strength reduction
runs here, so linear elements are sufficient; a finite element stability analysis
on the same mesh would need quadratic ones.

Click **Run → Run Seep…**, set **Run type** to **Steady**, leave **Convergence
tol** and **Max iterations** at `0.00010000` and `400`, and click **Run**. The
unconfined iteration converges in 59 sweeps to a total discharge of **0.16488
m³/day per m** of dam, with the head running from 2.000 to 18.000 m and the pore
pressure from −78.6 to 176.6 kPa. The reservoir boundary is bound to the `pool`
series, and a steady solve reads a series-bound value at t = 0 — which on this
schedule is full pool — so the Log opens with
`Series 'pool' read at t = 0 for the steady solve: reservoir value 18`.

![The full-pool steady solution](images/seep03_steady.png){width=1000}

Nearly the whole 16 m head drop happens inside the core, and the downstream slope
stays dry. [SEEP-3](seep03_reservoir_drawdown.md#the-steady-solution-at-full-pool)
reads this field in detail. Its role here is twofold: the march starts from it,
and the baseline analysis below is run against it.

### The baseline: the dam under a full reservoir

The curve needs a starting point, and that is the factor of safety with the
reservoir full and the steady solution just made under it. Switch back to
**LEM** (`Ctrl+1`) and click **Run → Run LEM…** **Method** opens on **Spencer**,
which satisfies both force and moment equilibrium and is the method behind every
factor of safety on this page. **Analysis** opens on **Auto search**, which finds
the run its own critical circle, and **Number of slices** at 40. Leave every
field alone.

Click **Run**. The Log reports the limit it read off the file, then the
refinement:

```text
Applying the file's search window: min_slip_depth.
Searching for the critical circular surface with SPENCER…
[🔁 iteration 1] center=(103.00, 59.00), FS=1.5372, grid=8.8500
[🔁 iteration 4] center=(103.00, 56.79), FS=1.5311, grid=2.2125
[🔁 iteration 7] center=(103.00, 56.79), FS=1.5311, grid=0.2766
[✅ converged] Iter=7, FS=1.5311 (ΔFS<0.0005) at (x=103.00, y=56.79, depth=0.00)
Critical FS = 1.531
Sliding mass = 5,580.9 kN/m over 57.93 m of failure surface
```

The first line names the one limit the file sets. What follows is the refinement
itself: the search screens both starting circles, takes the downstream one, and
walks its center down in halving steps until the factor of safety stops moving by
more than 0.0005.

**The factor of safety at full pool is 1.531, and it is on the downstream face.**
The critical circle is centered at (103.00, 56.79) and tangent to the rock at
elevation 0, breaking at the crest at (58.12, 22.00) and daylighting at the toe
at (109.16, 0.34): 5,580.9 kN/m of dam over 57.93 m of surface.

![Spencer at full pool: the critical circle is on the downstream face](images/combo03_solution_full.png){width=1000}

The reservoir is on the other side. Its blue arrows reach **176.6 kPa** at the
heel under 18 m of water and taper to zero at the waterline, and every one of
them presses the upstream slope into the dam. That is why the upstream face is
not critical here: 18 m of water is a large stabilizing load, and the downstream
slope — which carries the seepage the reservoir drives through the section and
none of its weight — is the weaker of the two. Under the failure surface the pale
blue band is the pore pressure on each slice base, read off the steady field, and
the green hatched band above it is the effective normal stress the strength is
computed from.

### The march

Now the pool comes down. Switch to **Seepage** and click **Run → Run Seep…**
again. Set **Run type** to
**Transient (time-dependent)** and click **Run**. **Convergence tol** and **Max
iterations** gray out — they belong to the steady solve, and the march sets its
own step size from how fast the field is moving.

This is the long run of the page; every other run here is seconds. The march
solves its initial condition first — the same field the steady run just produced
— and then prints a line per saved frame:

```text
Running transient seepage analysis…
  t=2: frame saved (steps so far=6, mass-balance closure=7.60e-04)
  t=5: frame saved (steps so far=18, mass-balance closure=2.40e-01)
  t=10: frame saved (steps so far=91, mass-balance closure=6.97e-03)
  t=15: frame saved (steps so far=175, mass-balance closure=7.88e-03)
  t=20: frame saved (steps so far=282, mass-balance closure=2.72e-02)
  t=25: frame saved (steps so far=417, mass-balance closure=1.95e-02)
  t=30: frame saved (steps so far=519, mass-balance closure=1.87e-02)
  t=35: frame saved (steps so far=626, mass-balance closure=7.96e-03)
  t=40: frame saved (steps so far=753, mass-balance closure=9.38e-03)
  t=47: frame saved (steps so far=917, mass-balance closure=2.97e-02)
  ...
  t=240: frame saved (steps so far=1689, mass-balance closure=1.67e-02)
  t=300: frame saved (steps so far=1697, mass-balance closure=1.01e-02)
Transient seepage complete — 19 saved frame(s).
```

**Nineteen frames**, at t = 0, 2, 5, 10, 15, 20, 25, 30, 35, 40, 47, 55, 65, 80,
100, 130, 180, 240 and 300. Those nineteen are the whole of what a curve can be
drawn through, so the schedule is a modeling decision made before the stability
question is asked, and this one is deliberately uneven: a frame every five days
through the drawdown, which ends on day 47, then widening steps through the
recovery. The answer moves fastest while the pool is falling and slowest once the
core is left draining on its own, and the frames follow that.

The **mass-balance closure** on each line is a diagnostic: the gap between the
stored-water change and the net inflow, as a fraction of the flow that has
passed so far. Early in the march that flow is still small, so the ratio at t = 5
reads large without meaning anything — the absolute gap there is a third of a
cubic meter per meter, no larger than at later frames, and the head field at
t = 5 moves by less than 0.04 m under a step size four times finer.

### One instant at a time

The march has left nineteen pore-pressure fields on the file, and the stability
run can read any one of them. Switch back to **LEM** and open
**Run → Run LEM…** once more. The form has grown
a group since the baseline run:

![Run LEM with the Seepage time group the march adds](images/combo03_studio_run_lem.png)

**Seepage time** is the group the march adds, and it names the instant this run
reads. Without it the choice would be invisible — it would fall to wherever the
results view's play bar happened to be sitting. Three ways to name an instant, in
the order they are reached for:

**Saved frame** lists the nineteen the march stored, so the run starts
immediately. **Frame shown in the results viewer** names whatever frame the play
bar is sitting on, and is offered only while a transient results tab is open.
**Another time (re-marches the solution)** accepts any instant inside the run
duration, and serves it by re-running the march with that time added to the
save schedule — a full re-solve, because a field blended between two frames
solves nothing and is never invented.

Under the list, **Save as the model's stability time** writes the choice to the
`tseep` sheet, which is what makes a scripted or headless re-run read the same
frame. Left off, the choice governs this run only. The dialog opens on the last
saved frame, t = 300, and the note under the fields says why: with the model's
stability time blank, the last frame is what any run that does not choose one
will read.

Under the circles warning the panel carries a note, the dependency stated back:

> Pore pressures come from the transient seepage solution, at t = 300 day. That
> frame is read into the model when the run starts, in place of any stored field.

Set **Saved frame** to `t = 35 day` and click **Run**:

```text
Applying the file's search window: min_slip_depth.
Searching for the critical circular surface with SPENCER…
[🔁 iteration 1] center=(7.00, 64.40), FS=1.3602, grid=9.5862
[🔁 iteration 5] center=(7.00, 57.21), FS=1.3319, grid=2.3965
[✅ converged] Iter=11, FS=1.3313 (ΔFS<0.0005) at (x=7.00, y=56.91, depth=0.00)
Critical FS = 1.331
Sliding mass = 5,722.8 kN/m over 58.01 m of failure surface
```

**1.331, and the circle has moved to the other side of the dam** — center
(7.00, 56.91) against the full-pool run's (103.00, 56.79). Two thirds of the
reservoir is gone, and with it the load that made the upstream face the safer
one.

That is one instant. Run LEM answers one at a time, which is what the **Seepage
time** group is for, and repeating it nineteen times would draw the curve by
hand.

#### The whole curve in one run

Click **Run → Parametric…** and set **Mode** to **Factor of safety vs time**:

![The Parametric dialog sweeping the march's saved frames](images/combo03_studio_parametric.png)

The parameter picker the other three modes use is gone, because nothing is
substituted here: every point solves *this* model against a different instant's
pore pressures. In its place is a **Saved frames** list holding the nineteen the
march stored, all ticked — **All** and **None** set the whole list, and unticking
samples a long march, since each frame is a full stability run. Leave **Method** on
**Spencer** and **Number of slices** at 40 and **Re-search the critical
surface at each step** ticked, which is the right setting here because the
mechanism moves. Leave **Grid search (auto-seed the circular search)** off: the
two circles the file carries already reach both faces.

Click **Run**. Nineteen searches, a few seconds each, stream their factors of
safety to the Log:

```text
Factor of safety vs time (lem): 19 instant(s) of the transient march, spencer, re-searching at each…
  t = 0 day      spencer   FS = 1.5311
  t = 2 day      spencer   FS = 1.5311
  t = 5 day      spencer   FS = 1.5132
  t = 10 day     spencer   FS = 1.4572
  ...
Lowest factor of safety 1.3313 at t = 35 day (19 instant(s), 0 without a result).
```

The run opens an **FS vs Time** tab with the curve, its lowest instant ringed and
labelled, and the reservoir schedule that drives it drawn faintly behind:

![The FS vs Time result tab](images/combo03_studio_fs_time.png)

**Day 35 reads 1.3313 here too, and day 0 reads the baseline's 1.5311** — the
sweep reads the same minimum slip depth off the circles sheet and searches from
the same two circles, so one dialog run and nineteen swept ones land on the same
numbers. An instant that produces no result comes back as a row carrying its
reason rather than as a gap in the curve.

From a script the same run is `xslope.sensitivity.fs_vs_time`:

```python
from xslope.sensitivity import fs_vs_time

ok, res = fs_vs_time(slope_data, solution, methods=("spencer",), search=True)
print(f"minimum FS = {res['min_fs']:.3f} at t = {res['critical_time']:g}")
```

Alongside the per-instant table it returns `critical_time` and `min_fs`, and each
row carries the `Xo`, `Yo` and `R` of the circle that instant settled on — which
is where the face each point came out on is read.
[Factor of safety versus time](../parametric/sensitivity.md#factor-of-safety-versus-time)
carries the mode in full.

### The curve

| t (day) | pool (m) | FS | center | R (m) | face |
|:---:|:---:|:---:|:---:|:---:|:---:|
| 0 | 18.00 | 1.5311 | (103.00, 56.79) | 56.79 | downstream |
| 2 | 18.00 | 1.5311 | (103.00, 56.79) | 56.79 | downstream |
| 5 | 16.93 | 1.5132 | (5.20, 65.90) | 65.88 | upstream |
| 10 | 15.16 | 1.4572 | (5.20, 65.90) | 65.88 | upstream |
| 15 | 13.38 | 1.4128 | (5.20, 65.90) | 65.88 | upstream |
| 20 | 11.60 | 1.3773 | (6.25, 60.58) | 60.58 | upstream |
| 25 | 9.82 | 1.3509 | (7.00, 57.81) | 57.75 | upstream |
| 30 | 8.04 | 1.3344 | (7.00, 56.91) | 56.91 | upstream |
| 35 | 6.27 | **1.3313** | (7.00, 56.91) | 56.91 | upstream |
| 40 | 4.49 | 1.3436 | (7.00, 56.91) | 56.91 | upstream |
| 47 | 2.00 | 1.3813 | (7.00, 57.16) | 57.16 | upstream |
| 55 | 2.00 | 1.4386 | (7.00, 57.16) | 57.16 | upstream |
| 65 | 2.00 | 1.4828 | (7.00, 57.16) | 57.16 | upstream |
| 80 | 2.00 | 1.5187 | (7.00, 57.16) | 57.16 | upstream |
| 100 | 2.00 | 1.5482 | (103.00, 56.79) | 56.79 | downstream |
| 130 | 2.00 | 1.5516 | (103.00, 56.79) | 56.79 | downstream |
| 180 | 2.00 | 1.5550 | (103.00, 56.79) | 56.79 | downstream |
| 240 | 2.00 | 1.5566 | (103.00, 56.79) | 56.79 | downstream |
| 300 | 2.00 | 1.5572 | (103.00, 56.79) | 56.79 | downstream |

Every circle above is tangent to the rock at elevation 0. The **face** column is
read off the center: the crest runs from x = 51 to x = 59, so a center left of it
is an upstream mechanism and a center right of it a downstream one.

![The factor of safety at every saved instant, over the pool that drives it](images/combo03_curve.png){width=1000}

**Day 2 repeats day 0 exactly**, at 1.5311 on the same downstream circle. The
pool has not moved yet, so the field has not either, and the two-day hold SEEP-3
entered to check the march checks the stability side as well.

**The lowest of the nineteen is 1.3313, on day 35** — 13% below the full-pool
1.5311, and **12 days before the drawdown ends**. The pool stands at elevation
6.27 there, with 4.3 m of the 16 m drawdown still to come.

**Day 47, the end of the drawdown, reads 1.3813** — 0.05 *above* the minimum, and
still on the upstream face. The instant the pool stops falling is not the instant
the dam is weakest.

**After the fall the curve climbs, then creeps.** From 1.3813 it rises to 1.4386
by day 55, 1.4828 by day 65 and 1.5187 by day 80, all on the upstream face; the
answer passes to the downstream mechanism at day 100 and then flattens: 1.5482,
1.5516, 1.5550, 1.5566, 1.5572 at days 100 through 300. The last three frames
together move it by two thousandths — the core draining, slice by slice, out of
the downstream shell, each step gaining less than the one before.

#### Which face governs, and when

The face column is the second reading the curve carries, and it is not the same
shape as the factor of safety.

**At full pool the downstream face governs, at 1.5311.** Eighteen meters of water
stand against the upstream slope and press it into the dam, and the searched
upstream mechanism at that instant is the safer of the two. The downstream slope
carries the seepage the reservoir drives through the section and none of its
weight.

**From the first frame of the drawdown the upstream face governs, and it holds
the answer for two thirds of the run.** By day 5 the pool has dropped about a
meter and the critical circle is already at (5.20, 65.90) on the upstream slope
at 1.5132; the answer stays there through the fall, the minimum and the whole
recovery, down to 1.3313 on day 35 and back up to 1.5187 on day 80.

**Only at day 100 does the downstream face take the answer back**, at 1.5482,
and it keeps it to the end: 1.5516, 1.5550, 1.5566, 1.5572. By then the upstream
slope has finished shedding the water the reservoir put into it and has climbed
past the downstream mechanism, which has barely moved all run.

So the reported curve is two mechanisms in sequence, not one surface moving. The
figure marks each point by the face it came out on, and the steep rise from day
35 to day 80 belongs entirely to the upstream mechanism recovering — the answer
does not change hands anywhere in it. The handover happens between day 80 and day
100, once the recovering upstream face crosses the downstream one at about 1.55.
A run given one starting circle on one face would draw a smooth curve through one
of those mechanisms and miss the other entirely.

**The mechanism also grows.** At full pool the critical circle carries 5,580.9
kN/m over 57.93 m of surface; the day-35 circle carries 5,722.8 kN/m over
58.01 m of the other slope. Both are tangent to the rock, and neither is the
surface the other run would have reported. That movement is why **Re-search the
critical surface at each step** matters: one circle solved at every instant would
report a drop, but on a mechanism that governs at one end of the run only.

#### The critical instant, drawn

![Spencer on day 35, near the end of the drawdown](images/combo03_solution_min.png){width=1000}

The failure surface is on the other slope from the full-pool figure, and the
reservoir is why. Its load has fallen from **176.6 kPa** at the heel to **61.5
kPa**, and it reaches only to elevation 6.3 on the face instead of 18 — both the
pressure at its base and the height it acts over are now about a third of what
they were. Inside the dam the contours are still crowded over the core, which is
holding head the pool no longer balances. The surface runs the full height of the
upstream slope, from the heel at (0.80, 0.34) to the crest edge at (51.94, 22.00),
tangent to the rock.

### Choosing the time steps

The curve's minimum sits on a saved frame, which raises the question the save
schedule always raises: whether the true minimum falls between two of them. Here
the schedule was built to answer it — five-day frames across the whole drawdown,
because that is the stretch where the answer moves — so the check is a short one.

Asking for a set of unsaved instants at once is `fs_vs_time` with `times=` and
the seepage data to re-march from, which serves the whole set with **one**
re-march before the first solve rather than one per instant:

```python
ok, fine = fs_vs_time(slope_data, solution, seep_data=seep_data, remarch=True,
                      times=[20, 22.5, 25, 27.5, 30, 32.5, 35, 37.5, 40],
                      methods=("spencer",))
```

Those nine halve the resolution the schedule already has across the dip, and four
of them name no saved frame:

| t (day) | pool (m) | FS |
|:---:|:---:|:---:|
| 20 | 11.60 | 1.3773 |
| 22.5 | 10.71 | 1.3620 |
| 25 | 9.82 | 1.3510 |
| 27.5 | 8.93 | 1.3409 |
| 30 | 8.04 | 1.3343 |
| 32.5 | 7.16 | 1.3323 |
| 35 | 6.27 | **1.3312** |
| 37.5 | 5.38 | 1.3360 |
| 40 | 4.49 | 1.3434 |

**The minimum stays on day 35, at 1.3312 against the saved grid's 1.3313** — one
ten-thousandth, on a curve whose total swing is 0.23, and well inside the
search's own convergence tolerance of 0.0005. Every instant governs on the
upstream face.

A denser schedule is what makes that check cheap to trust, and a coarse one is
what makes it necessary: SEEP-3's own twelve-frame schedule saves nothing between
day 30 and day 47, so a curve drawn through it would report a minimum of 1.3344
on day 30 and never evaluate the instant that is actually lowest. The time to spend a re-march is when the saved
frames are far enough apart that the minimum falls on the first or last of them,
or when the curve's neighbors on either side of the minimum are still falling
steeply toward it.

### Reading the curve

Two things happen to the upstream slope when the pool is lowered, and they do not
happen at the same speed.

**The stabilizing load leaves immediately.** The reservoir presses on the face,
and that pressure resists the slide directly. It tracks the pool exactly — 176.6
kPa at the heel on day 0, 61.5 on day 35, and 19.6 from day 47 on, where all that
stands against the heel is the 2 m of water left at tailwater level — because a
load is a boundary condition rather than something the soil has to release.

**The pore pressure inside leaves slowly.** The shell can follow the pool: its
vertical conductivity over its drainable porosity is 0.5 / 0.22 = 2.3 m/day
against the pool's 0.36 m/day. The core cannot, at 0.005 / 0.03 = 0.17 m/day,
and it holds the head the full reservoir put into it long after the reservoir has
gone — 8.6 m of it still standing on day 47, which
[SEEP-3 measures](seep03_reservoir_drawdown.md#reading-the-frames).

The upstream factor of safety falls while the first is ahead of the second and
rises once the second catches up. That crossover, not the end of the drawdown, is
what day 35 marks. Over the 12 days after it the pool gives up its last 4.3 m
while the shell behind the face keeps draining, and the strength recovers faster
than the remaining load falls — 1.3436 by day 40, 1.3813 by day 47, 1.5187 by day
80. After that the curve follows the core emptying, which takes months and acts
on the downstream shell.

#### The relationship to a rapid-drawdown check

[COMBO-2](combo02_rapid_drawdown.md) runs the other analysis of the same
phenomenon: the three-stage Duncan, Wright and Wong procedure, on a dam whose core
carries an undrained $K_c = 1$ envelope. That procedure names **two** instants —
the slope before the drawdown and the slope after it — and computes a single
factor of safety for the second from consolidation stresses read at the first.
Where a transient march supplies the water, those two instants are two frames of
it, named by `stage_1` and `stage_2` on the `tseep` sheet.

A curve is the sequence version of that pair. It cannot replace the drawdown
check, because a drawdown check is one three-stage analysis with undrained
strengths in it rather than a sequence of single-stage analyses — the two answer
different questions, and this model, carrying no $d$ / $\psi$ pair, cannot be
asked the first. What the curve adds is *which* instant the second stage should
be read at. On this dam a check staged at the end of the drawdown would be
reading day 47, at **1.3813**, with the upstream mechanism already 0.05 into its
recovery; the worst of the nineteen saved instants is day 35, at **1.3313**. The
curve locates it, and once located it can be named as a stage time and the
drawdown check run there.

---

## Part 2 — Rapid drawdown at every instant

Part 1's curve is a sequence of drained analyses, which is the right question for
a slope that sheds pore water about as fast as the pool falls. A compacted clay
core does not, and the section above says what that costs: the model carries no
$d$ / $\psi$ pair, so it cannot be asked the rapid drawdown question at all.

[COMBO-2](combo02_rapid_drawdown.md)'s Johnson Reservoir dam can. Its core
carries a $K_c = 1$ envelope, and COMBO-2 runs the three-stage Duncan, Wright and
Wong procedure on it three times — once per statement of where the water is —
reading 1.181, 1.195 and 1.016. The last of those comes from a transient march,
with stage 1 read at t = 0 and stage 2 at t = 50, the instant the pool reaches
its residual level. That is one drawdown at one chosen instant. Part 2 runs the
same procedure at **every** saved instant of that march, and then runs Part 1's
kind of curve on the same twenty-one frames so the two can be set side by side.

### The model, and what Part 2 reads off it

The march the sweep reads is the one COMBO-2 already solved, so Part 2 opens it
rather than producing it again. Download
[xslope_johnson_fs_time.xlsx](files/xslope_johnson_fs_time.xlsx) and open it with
**File → Open…**, leaving the toolbar's mode strip on **LEM**. The mesh and the
march travel with it — both are stored beside the workbook and Studio reads them
on open — so the dam arrives meshed and marched, and the Log names each result as
it is restored. In **LEM** the inputs view draws that mesh behind the zones,
which the workbook on its own does not carry.

![The Johnson Reservoir dam as the shipped file opens, the mesh behind it](images/combo03_rapid_inputs.png){width=1000}

The section is 750 ft long: a 100 ft foundation on rock at elevation 0 and an
80 ft embankment on it, a sand **shell** on both faces and a compacted-clay
**core** carried 40 ft into the foundation as a cutoff key. The reservoir stands
at elevation 160 and the tailwater at 100, and the drawdown lowers the pool 50 ft
to a residual reservoir at 110. [SEEP-2](seep02_johnson_dam.md) builds the
geometry; COMBO-2 adds everything the drawdown needs. Three of those additions
are what this part reads.

**The undrained envelope.** Click **Materials** and set the **Show parameters
for:** toggles to **LEM** alone. Only the core carries a $K_c = 1$ envelope —
$d = 250$ psf and $\psi = 14°$, under its drained $c' = 400$ psf and
$\phi' = 18°$ — and the shell and the foundation are blank on both columns, which
declares them free-draining through the drawdown.

| mat | name | γ (pcf) | c′ (psf) | φ′ (deg) | d (psf) | ψ (deg) | u |
|:---:|---|:---:|:---:|:---:|:---:|:---:|:---:|
| 1 | `shell` | 130 | 100 | 35 | — | — | `seep` |
| 2 | `core` | 125 | 400 | 18 | 250 | 14 | `seep` |
| 3 | `foundation` | 127 | 100 | 27 | — | — | `seep` |

[COMBO-2 carries the procedure itself](combo02_rapid_drawdown.md#the-three-stage-procedure)
— what each stage computes, and how the strength stage 2 uses is interpolated
between the undrained envelope and the drained one according to the stress ratio
each slice consolidated at.

**`u` is `seep` on all three rows**, so every slice base takes its pore pressure
from a solved seepage field. The file's two piezometric lines were deleted when
COMBO-2 replaced them, which is why the figure above draws water surfaces and
their derived loads and no piezometric line at all.

**The schedule and the stage times.** Switch the mode strip to **Seepage**
(`Ctrl+2`) and click **Transient** in the Inputs dock. The `pool` series holds at
elevation 160 for five days and falls to 110 over the following 45; the run lasts
500 days; and the save times are listed outright, a frame every five days from
the start of the hold to the end of the fall and then widening steps through the
recovery. **Stage 1 time** is `0` and **Stage 2 time** is `50`. Click **OK**.

That schedule is this file's own — [COMBO-2](combo02_rapid_drawdown.md)'s
workbook keeps a coarser one and a 1,000-day run — and the reason is the sweep
below. A curve resolves nothing between two frames, so the instants are packed
where the answer moves and spread where it does not, and the run stops at day 500
because the curve has flattened well before it.

Every point of the curve below reads its consolidation stresses at stage 1, and
stage 1 stays at t = 0 for the whole sweep. What moves from point to point is
stage 2.

### The mesh and the march it carries

The schedule states what the reservoir does; the saved frames hold what the dam
did about it. Still in **Seepage**, click the **Seep · Transient** tab. The play
bar steps through **twenty-one** instants — t = 0, 5, 10, 15, 20, 25, 30, 35, 40,
45, 50, 60, 70, 80, 100, 130, 170, 220, 300, 400 and 500 — and every pore pressure
the runs below read comes from one of them.

Those fields were solved on **2,080 nodes and 3,923 triangles**: linear triangles
auto-sized at 100 divisions across the 750 ft section, which puts the target
element size at 750/100 = 7.5 ft.
[COMBO-2 builds that mesh](combo02_rapid_drawdown.md#meshing-and-both-solves) with
**Run → Build Mesh…** and
[marches the pool on it](combo02_rapid_drawdown.md#marching-it) with
**Run → Run Seep…**, a run of several minutes on this dam. The file downloaded
above arrives past both, so nothing below waits on a seepage solve.

### Running the rapid drawdown sweep

Those twenty-one fields are what the sweep questions, one at a time. Switch back to
**LEM** (`Ctrl+1`) and click **Run → Parametric…**

![The Parametric dialog with the Rapid drawdown box ticked](images/combo03_rapid_parametric.png)

Set **Mode** to **Factor of safety vs time**. **Method** opens on **Spencer**,
which satisfies both force and moment equilibrium and is the method behind every
number in this part, and **Number of slices** opens at 40; leave both.

**Saved frames** lists the twenty-one this march stored, all ticked. Leave them so:
the whole march is what the comparison at the end of this part needs.

**Rapid drawdown at each time → ticked.** Ticking it changes what a ticked
instant means. Left off, each instant is one single-stage analysis of that
instant's water. Ticked, each instant becomes stage 2 of a three-stage drawdown
whose stage 1 is the march's initial state, so the curve answers how safe the
slope is if the pool falls from full to where it stands at that moment.

**Re-search the critical surface at each step** goes gray and is held on. A
drawdown's critical surface is not the drained one and it moves as the drawn-down
field changes, so there is nothing left for the toggle to decide: every point is
searched from the file's starting circle at (275, 235). That circle is on the
upstream face, which is where a drawdown is a problem: lowering the pool takes
the water off that slope and leaves it inside the embankment. Every instant of
the sweep interrogates the same face.

The note under the two boxes rewrites itself to match, and the **Model checks**
panel carries one warning, which repeats COMBO-2's: two of the three materials
carry no $d$ / $\psi$ and keep their drained strength through the drawdown.
Boundary set 2 raises nothing here, because
[COMBO-2 clears it](combo02_rapid_drawdown.md#clearing-boundary-set-2) before its
own march and this file arrives without it. A transient run solves boundary set 1
only, and both stages come out of that march.

Click **Run**. Twenty of the twenty-one instants are drawdowns and each is
searched in full, so the Log fills slowly; the table lands when the last search
finishes:

```text
Rapid drawdown vs time (lem): 21 instant(s) of the transient march, spencer, re-searching at each…
Rapid drawdown versus time (stage 1 at t = 0)
t (day)  stage 1  stage 2       stage 3      FS  governs      Xo      Yo       R
      0        —        —             —       —        —       —       —       —   a drawdown from t = 0 to t = 0 is not a drawdown …
      5   1.5109   1.4563  not required  1.4563        2  242.74  256.34  165.90
     10   1.5142   1.3496  not required  1.3496        2  244.06  253.05  164.14
     15   1.5248   1.2704  not required  1.2704        2  246.78  250.47  163.54
     20   1.5352   1.2036  not required  1.2036        2  249.59  244.90  159.60
     25   1.5414   1.1489  not required  1.1489        2  249.59  243.48  159.38
     30   1.5440   1.1051  not required  1.1051        2  249.59  243.48  159.76
     35   1.5453   1.0711  not required  1.0711        2  248.17  244.19  161.02
     40   1.5465   1.0457  not required  1.0457        2  247.47  244.90  162.02
     45   1.5525   1.0279        1.0278  1.0278        3  244.64  244.19  162.68
     50   1.5545   1.0163        1.0157  1.0157        3  243.93  244.90  163.64
     60   1.5518   1.0523        1.0522  1.0522        3  245.35  242.07  160.48
     70   1.5506   1.0743  not required  1.0743        2  245.31  246.38  164.30
     80   1.5512   1.0902  not required  1.0902        2  245.31  244.96  163.08
    100   1.5514   1.1159  not required  1.1159        2  245.31  243.53  161.82
    130   1.5514   1.1414  not required  1.1414        2  245.31  243.53  161.82
    170   1.5514   1.1612  not required  1.1612        2  245.31  243.53  161.82
    220   1.5514   1.1744  not required  1.1744        2  245.31  243.53  161.82
    300   1.5514   1.1848  not required  1.1848        2  245.31  243.53  161.82
    400   1.5514   1.1893  not required  1.1893        2  245.31  243.53  161.82
    500   1.5514   1.1912  not required  1.1912        2  245.31  243.53  161.82
Lowest factor of safety 1.0157 at t = 50 day (21 instant(s), 1 without a result).
```

The t = 0 row's message is cut at the ellipsis to fit the page. In the Log it runs
on: *"a drawdown from t = 0 to t = 0 is not a drawdown -- the pool has not fallen
yet; stage 2 must be a later instant than stage 1."* At t = 0 stage 2 would be
stage 1, a fall from the full pool to itself, so that instant comes back as a row
carrying its reason rather than as a point on the curve. Twenty of the twenty-one
frames are drawdowns; the first is the state they all fall from.

### The rapid drawdown curve

The run opens a **Drawdown vs Time** tab — the same tab a single-stage sweep
uses, renamed for the curve it is holding:

![The rapid drawdown factor of safety at every saved instant, with its three stages behind it](images/combo03_rapid_curve.png)

One curve is drawn, and it carries the rapid drawdown factor of safety at each
instant — the lower of stages 2 and 3. The three stages themselves stay in the
table above, where the governing one is named per row; a figure holding a line
per stage would draw three answers where only one of them is reported. Behind
the curve, on the right axis, runs the pool schedule that drives it. The lowest
instant is ringed and labelled with its own time, the red guide marks FS = 1 —
drawn here because the curve comes within its own swing of it — and the legend
counts the instant that produced no result, so the missing t = 0 point is never
a silent gap.

**Day 50 reads 1.0157, and that is COMBO-2's answer.** COMBO-2's transient route
reports **1.016** on a circle at (243.93, 244.90) with a radius of 163.64 ft, and
the day-50 row above carries the same factor of safety on the same circle to the
hundredth of a foot. Nothing about that is a check of one page against another:
at t = 50 this sweep runs precisely the rapid drawdown COMBO-2 ran — stage 1 at 0,
stage 2 at 50, the same march, the same starting circle — so the curve passes
through COMBO-2's answer at day 50 because that answer is one of its points.

**Stage 3 runs on three rows, all of them at the bottom of the dip.** Days 45, 50
and 60 read 3 in the **governs** column, and on each the margin is a single
ten-thousandth — 1.0278 against stage 2's 1.0279, 1.0157 against 1.0163, 1.0522
against 1.0523. Every other row reads `not required`, which means no core slice
came out weaker drained than undrained and stage 2 is the answer.
[COMBO-2's sweep of the core's intercept](combo02_rapid_drawdown.md#the-governing-stage-flip)
finds why those margins are so thin: the handover from stage 2 to stage 3
falls at $d = 223$ psf, 27 psf under the 250 psf this core carries.

**Day 5 costs 0.05 with no drawdown in it.** The pool holds at elevation 160
until day 5, so at that instant stage 2's water is stage 1's water — the same
field, the same reservoir load on the face. The reported 1.4563 still sits 0.055
below the 1.5109 stage 1 reads on the same circle, and the whole of that gap is
the undrained envelope being applied to the core. It is the price of the strength
substitution alone, measured with the water held still.

**The minimum, 1.0157, falls on day 50 — the instant the pool stops falling.**
That is not a rule. Two things race while the pool drops: the support the
reservoir gives the face keeps falling until the pool stops, and the pore
pressures in the core start dissipating from the moment the drop begins. On this
dam the drawdown is fast against a core at 0.001 ft/day, the pore pressures barely
move during the 45 days, and the worst instant is where the load is lowest — the
end. Part 1's dam drained fast enough for dissipation to overtake the load loss
mid-fall, and its minimum came 12 days *before* its drawdown ended. A slow enough
drawdown shows no dip at all. Which of those a given dam does is exactly what
the sweep is for; the instant cannot be named in advance.

The saved times bracket that minimum closely, which is what they were chosen for:
day 45 reads 1.0278 and day 60 reads 1.0522, so the fall is resolved to within
five days of its end and the curve turns at a frame rather than between two.

**The recovery climbs toward COMBO-2's two-steady answer.** From 1.0157 the curve
rises to 1.0902 by day 80, 1.1414 by day 130, and then slows: 1.1612, 1.1744,
1.1848, 1.1893, 1.1912 at days 170 through 500, the last two frames together
moving it by 0.0019. COMBO-2's **two-steady-solution** answer is **1.195**, and
the curve is still 0.004 under it at day 500 and still climbing. Set 2 solves the
dam after the pore pressures have fully caught up with the residual pool, which a
march of finite length only approaches: by day 500 the core has nearly, but not
quite, finished draining.

**Stage 1 moves, and the water is not why.** Its column reads 1.5514 on the six
rows from day 100 on, 1.5512 on day 80, 1.5545 on day 50, 1.5453 on day 35 and
1.5109 on day 5. Every one of those is the same physical state — full pool,
drained, the march's initial field — so the spread is the surface it was
evaluated on. Each row searches for the circle that minimizes the *drawdown*, and
reports stage 1 on whichever circle that turned out to be.

### The same march at drained strengths

The question Part 1 asks can be put to this dam too, on the same twenty-one frames,
by leaving the **Rapid drawdown at each time** box unticked. Open **Run → Parametric…** again, keep every
other field, untick **Rapid drawdown at each time**, and click **Run**.
**Re-search the critical surface at each step** comes back live and stays ticked,
which is the right setting for the same reason it was held on before.

Each point now runs about three times faster than its drawdown counterpart,
because it solves the section once rather than three times:

```text
Factor of safety vs time (lem): 21 instant(s) of the transient march, spencer, re-searching at each…
Factor of safety versus time
t (day)      FS      Xo      Yo       R
      0  1.5097  240.00  261.13  169.97
      5  1.5097  240.00  261.13  169.97
     10  1.3944  242.71  254.34  164.71
     15  1.3093  244.00  250.40  162.76
     20  1.2374  248.18  242.77  157.49
     25  1.1779  249.59  242.07  157.77
     30  1.1301  249.59  239.24  155.60
     35  1.0930  248.88  239.24  156.31
     40  1.0656  246.76  241.36  159.20
     45  1.0469  246.76  241.01  158.89
     50  1.0350  243.89  244.96  163.73
     60  1.0765  245.31  242.11  160.55
     70  1.1037  245.31  244.96  163.08
     80  1.1244  245.31  243.53  161.82
    100  1.1587  245.31  243.53  161.82
    130  1.1949  245.31  243.53  161.82
    170  1.2252  245.31  243.53  161.82
    220  1.2475  244.60  244.60  163.07
    300  1.2684  245.67  243.18  161.32
    400  1.2802  244.60  242.82  161.51
    500  1.2870  244.60  244.25  162.76
```

Every instant produces a result here, including t = 0, which is an ordinary
analysis of the dam at full pool rather than a fall from itself. **Day 0 and day 5
read 1.5097 on the same circle**, because the pool has not moved between them and
neither has the field.

Two of these rows are COMBO-2's drained bracket runs, arrived at from a different
direction. **Day 0's 1.5097** is its full-pool drained figure of **1.510**, and
**day 50's 1.0350** is its day-50 drained figure of **1.035** — the same frame,
the same strengths. Its third drained figure, the long-term **1.311**, is not on
this curve: day 500 reads 1.2870 and is still climbing, because a march of this
length leaves the core part drained.

Studio's result tab holds one curve at a time, so the two are set on one pair of
axes here, drawn for this page. Running the sweep in both modes and plotting the
two results together reproduces it:

![The rapid drawdown curve and the single-stage curve of the same march](images/combo03_rapid_compare.png){width=1000}

**The rapid drawdown curve is the lower of the two at every instant**, and the distance
between them is not constant:

| t (day) | drawdown | single-stage | difference |
|:---:|:---:|:---:|:---:|
| 5 | 1.4563 | 1.5097 | 0.0535 |
| 10 | 1.3496 | 1.3944 | 0.0448 |
| 15 | 1.2704 | 1.3093 | 0.0389 |
| 20 | 1.2036 | 1.2374 | 0.0338 |
| 25 | 1.1489 | 1.1779 | 0.0290 |
| 30 | 1.1051 | 1.1301 | 0.0250 |
| 35 | 1.0711 | 1.0930 | 0.0219 |
| 40 | 1.0457 | 1.0656 | 0.0199 |
| 45 | 1.0278 | 1.0469 | 0.0191 |
| 50 | **1.0157** | **1.0350** | **0.0193** |
| 60 | 1.0522 | 1.0765 | 0.0244 |
| 70 | 1.0743 | 1.1037 | 0.0294 |
| 80 | 1.0902 | 1.1244 | 0.0342 |
| 100 | 1.1159 | 1.1587 | 0.0429 |
| 130 | 1.1414 | 1.1949 | 0.0536 |
| 170 | 1.1612 | 1.2252 | 0.0640 |
| 220 | 1.1744 | 1.2475 | 0.0731 |
| 300 | 1.1848 | 1.2684 | 0.0835 |
| 400 | 1.1893 | 1.2802 | 0.0909 |
| 500 | 1.1912 | 1.2870 | 0.0958 |

**The gap is narrowest where the slope is weakest.** Across the bottom of the dip
— days 40 to 50 — the two analyses run 0.019 to 0.020 apart, about 1.9% of the
rapid drawdown answer; by day 500 they are 0.096 apart, 8.0% of it. The undrained
envelope changes the answer least at the moment the answer matters most.

The pore pressure in the core explains both halves of that. At day 50 the core
still holds the head the full reservoir put into it, so the
*drained* strength there is already cut by that pore pressure and lands close to
the undrained strength stage 2 assigns — the same reading
[COMBO-2 takes](combo02_rapid_drawdown.md#the-three-answers-and-what-brackets-them)
from its bracket table, where the same envelope moves one field by 0.12 and the
other by 0.02. Afterward the two analyses part company: the drained strength
recovers as the core drains, and the undrained strength cannot, because stage 2
computes it from stage-1 consolidation stresses that are the same at every point
of the curve. One curve climbs to 1.29 and the other levels off near 1.19.

### Which curve to read

The two curves answer two questions, and neither is the other's conservative
version.

**Through the fall and just after it, read the rapid drawdown curve.** Its number for
this dam is **1.016, on day 50**. During those 45 days the core cannot shed the
reservoir's head, so the strength it can mobilize is the undrained one, and a
single-stage drained analysis at the same instant credits it with a strength it
does not have. That the single-stage curve reads 1.0350 there — only 0.019
higher — is a property of this dam, not of the method. The gap is the
shear-induced pore pressure the undrained envelope adds on top of the retained
head the transient solution already carries, and on this core the two strengths
sit close together at the stresses the surface sees
([COMBO-2](combo02_rapid_drawdown.md#the-governing-stage-flip) measures the
crossover at 223 psf against the core's 250). A more contractive core — a lower
$d$, a lower $\psi$ — a surface that crosses more of it, or a faster drawdown
would open that gap well beyond 0.02, with the single-stage curve reading too
high the whole way down. The sweep is what measures the gap for the dam at hand;
it costs one run.

**Long after the fall, read the single-stage curve.** By day 500 the core has
largely drained, and 1.2870 is what the slope actually has. The rapid drawdown
curve's 1.1912 at that instant is still answering "what if the pool fell from full
to here",
which by then is not what happened. It describes the long-term drawdown state,
the case COMBO-2's two steady solutions compute, and it applies only if the
reservoir is filled and drawn down again.

For a zoned dam with a core that does not drain, that puts the design number on
the rapid drawdown curve's minimum and the long-term number on the single-stage curve's
plateau. The sweep is what locates the first of those, without anyone having to
name the instant in advance.

---

## Conclusion

This tutorial covered:

- What a factor-of-safety-versus-time curve is: one stability analysis per saved
  frame of a transient seepage solution, with no input changing between the
  points, so the axis is time and everything the curve does is the water moving.
- How a frame reaches a stability run — `u = seep` on the materials table sends
  every slice base to a solved field, and Run LEM's **Seepage time** group names
  which frame that is, with the option to write the choice back to the model;
  **Run → Parametric… → Factor of safety vs time** sweeps all of them at once.
- Where to put a starting circle when either face may govern: one per face, each
  drawn on the deep mechanism its slope can make — entering near the crest,
  tangent to the rock, daylighting past the toe — because a seeded search refines
  only the best-screening circle, and a shallower seed settles into a smaller
  mechanism near where it started.
- The curve on Part 1's dam — 1.531 at full pool on the **downstream** face, a
  minimum of **1.331 on day 35** on the **upstream** face, 1.381 still upstream at
  the end of the drawdown on day 47, and a climb to 1.557 by day 300 as the core
  empties.
- How to read which face governs and when: downstream while the reservoir stands,
  upstream from day 5 through day 80 — the fall, the minimum and the whole
  recovery — and downstream again from day 100 on, once the upstream slope has
  climbed past it. The reported curve is two mechanisms in sequence, and the steep
  rise from day 35 to day 80 belongs entirely to the upstream one.
- That the saved times set the curve's resolution, so the schedule is packed where
  the answer moves: five-day frames across the drawdown, widening afterward.
  Halving that spacing across the dip — nine instants, four of them new, one
  re-march — left the minimum on day 35 and moved it by 0.0001.
- Why that minimum falls mid-drawdown: the reservoir load tracks the pool exactly
  while the pore pressure behind the face lags, and the upstream slope is weakest
  where the gap between them is widest — 12 days before the pool stops moving.
- How to turn the same sweep into a rapid drawdown curve — **Rapid drawdown at each
  time** on the Parametric dialog makes every ticked instant stage 2 of a
  three-stage analysis whose stage 1 is the march's initial state, holds the
  re-search toggle on, and refuses the initial instant itself, which would be a
  fall from the full pool to itself.
- The rapid drawdown curve on COMBO-2's Johnson Reservoir dam — 1.4563 on day 5
  with the pool still standing, a minimum of **1.0157 on day 50** where the fall
  ends, and a recovery to 1.1912 by day 500. Day 50 is COMBO-2's own
  transient-route answer of 1.016, on the same circle, because that instant is the
  drawdown COMBO-2 ran; stage 3 governs there and on the two frames either side of
  it, by a ten-thousandth each time.
- What the same twenty-one frames give at drained strengths — 1.5097 at full pool,
  1.0350 on day 50 and 1.2870 by day 500 — and that the rapid drawdown curve is
  the lower of the two at every instant, by 0.019 at day 50 and 0.096 by day 500.
  The undrained envelope changes the answer least where the slope is weakest,
  because the retained pore pressure has already cut the drained strength there.
- Which curve to read: the rapid drawdown curve through the fall and just after it,
  where the core cannot shed the reservoir's head; the single-stage curve for the
  long-term condition, once it has.

**Where to go next:** [SEEP-3](seep03_reservoir_drawdown.md) builds the transient
seepage model of Part 1's dam from nothing and reads it as frames, head histories
and a water budget. [COMBO-2](combo02_rapid_drawdown.md) builds Part 2's model
and runs the three-stage rapid-drawdown procedure on it at one instant, supplying
its two water states three ways; [SEEP-2](seep02_johnson_dam.md) builds that dam.
[Factor of safety versus time](../parametric/sensitivity.md#factor-of-safety-versus-time)
carries the sweep mode itself — the times, the re-march rule, the per-instant row
contract and the FEM variant. [Transient seepage](../seep/transient.md) carries
the formulation and the boundary types, and
[Rapid drawdown analysis](../lem/rapid.md) the three-stage procedure. The
[tutorials index](index.md) lists the series.
