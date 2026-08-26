---
title: "Tutorial COMBO-3 — Factor of Safety vs Time"
description: "Two dams run for stability against every saved frame of a transient drawdown — SEEP-3's earth dam at drained strengths, where the lowest factor of safety falls 17 days before the reservoir stops moving, and COMBO-2's Johnson Reservoir dam as a three-stage rapid drawdown at each instant, where the curve passes through COMBO-2's own answer at day 50 and stays below the drained curve for the rest of the run."
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
equilibrium analysis of the upstream slope — the side a drawdown attacks, since
lowering the pool takes the water off that face and leaves the water inside the
embankment — and runs the analysis against every frame the march saved.

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
one instant at a time, and where that instant is chosen; how the saved times set
the curve's resolution and what a finer set of them costs; how to read the
curve's minimum against the drawdown that caused it; and how to run a three-stage
rapid drawdown at every instant of a march and compare that curve against the
drained one.
</div>
<p><span class="tg-pill">two materials</span><span class="tg-pill">three materials</span><span class="tg-pill">transient seepage</span><span class="tg-pill">u = seep</span><span class="tg-pill">saved frames</span><span class="tg-pill">seepage time</span><span class="tg-pill">FS vs time</span><span class="tg-pill">rapid drawdown</span><span class="tg-pill">three-stage procedure</span><span class="tg-pill">Kc = 1 envelope</span><span class="tg-pill">d and ψ</span><span class="tg-pill">stage times</span><span class="tg-pill">Spencer</span><span class="tg-pill">circular search</span><span class="tg-pill">search window</span><span class="tg-pill">minimum slip depth</span><span class="tg-pill">automatic water loads</span></p>
<div class="tgm-model" markdown>
**Part 1 model** — [xslope_earth_dam_fs_time.xlsx](files/xslope_earth_dam_fs_time.xlsx),
SEEP-3's completed dam with a strength band on its materials table, a starting
circle on the upstream face and a search window. It carries no mesh and no
solution, so every run below is made from scratch

**Part 2 model** — [xslope_johnson_fs_time.xlsx](files/xslope_johnson_fs_time.xlsx),
COMBO-2's completed Johnson Reservoir dam: an undrained envelope on the core,
`u = seep` on all three zones, the pool schedule and both stage times. This copy
ships meshed and marched, so it opens with the twelve pore-pressure fields Part 2
sweeps already on it; [COMBO-2](combo02_rapid_drawdown.md) builds both
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
45 days, and held there for the remaining 313 days of a 360-day run. The two
zones drain at very different rates — the shell follows the pool down, the core
does not — and the head the core keeps after the pool has gone is the reason a
drawdown is a stability problem at all.
[SEEP-3](seep03_reservoir_drawdown.md) builds that model from nothing: the
storage properties, the boundary set, the `pool` time series and the saved-frame
schedule. Nothing about the seepage side changes here.

#### What the curve is made of

Each saved frame carries a pore pressure at every node of the mesh. A limit
equilibrium run with `u = seep` reads one such field at every slice base, so one
frame produces one factor of safety. Twelve frames produce twelve, and plotting
them against the frame times draws the curve. No input changes between the
points — the same section, the same strengths, the same solver at every one — so
everything the curve does is the pore pressures moving and the reservoir load
leaving the face.

Two consequences follow, and both matter further down. **The curve can only be
as fine as the save schedule**, because a time that is not a saved frame has no
field to read and is never interpolated between two that are. And **the critical
surface moves**: the mechanism that governs a dam under a full reservoir is not
the one that governs it half drained, so each instant gets its own search.

### The strengths this page adds

The first thing to see is what the file carries beyond SEEP-3's seepage
properties. Open [xslope_earth_dam_fs_time.xlsx](files/xslope_earth_dam_fs_time.xlsx) with
**File → Open…**. The toolbar's mode strip reads LEM | Seepage | FEM; leave it on
**LEM** for now and click **Materials** in the Inputs tree. Set the **Show
parameters for:** toggles to **LEM** alone:

![The two zones with their strength band filled](images/combo03_studio_materials.png)

| mat | name | γ (kN/m³) | γsat (kN/m³) | option | c′ (kPa) | φ′ (deg) | u |
|:---:|---|:---:|:---:|:---:|:---:|:---:|:---:|
| 1 | `shell` | 20 | 21 | `mc` | 0 | 35 | `seep` |
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

### The surface and the search window

The strengths say how strong the slope is; the circles say where the search
looks for the failure. Click **Circles**. The file carries one starting circle, on the **upstream** face:

| Xo | Yo | Option | Depth |
|:---:|:---:|:---:|:---:|
| 25 | 44 | `Depth` | 0 |

Center at x = 25, mid-way along a face that runs from the heel at x = 0 to the
upstream crest edge at x = 51, and at y = 44, two dam heights above the toe.
**Depth** is an elevation rather than a distance, and 0 is the rock the dam
stands on. A search refines away from this circle, so its job is to point the
refinement at the right part of the section, not to be the answer.

The file also declares a **search window** — a set of limits confining where a
searched surface may run. There are ten, and they are edited in the **Search
window** group under this editor's table — the same cells as the **Search window**
block at the right of the `circles` worksheet. A blank field is a limit that is not
applied; [the template page](../usage/input_template.md#search-window-optional)
lists all ten. Three of them are set here, filling five of the group's fields:

| Limit | Value | What it does |
|---|:---:|---|
| Entry x min / max | 42 / 59 | The crest-side end of the trace must land between the full-pool waterline on the face and the downstream crest edge |
| Exit x min / max | 0 / 42 | The toe-side end must land on the upstream face below that waterline |
| Min slip depth | 8 | The surface must reach at least 8 m below the ground surface |

The first two hold the search on the **upstream** slope. Entry names the
crest-side end of the trace and exit the toe-side end, whichever way the slope
faces, and the reservoir meets the face at (42, 18) — so the pair reads as *break
above the full-pool waterline, daylight below it*. Without them the critical
surface at full pool is a toe circle on the **downstream** slope, which the
reservoir never touches and a rapid drawdown curve says nothing about.

The third limit sets a floor on the *size* of the mechanism, and the shell is why
it takes one. A cohesionless soil has no length scale in it: with c′ = 0 the
factor of safety on a shallow surface parallel to the face does not depend on how
deep that surface is drawn, so an unconstrained search walks down to a surficial
wedge whose failure is raveling rather than a slide. 8 m on a 22 m dam holds the
search on slides of the embankment. This is Slide2's minimum-depth filter, and
the Run LEM dialog offers the same limit as **Ignore surficial (skin) failures**
for a run that wants to set it rather than read it off the file.

The preview draws the window along with the circle: the two bars lying on the
ground surface are the exit range on the upstream face and the entry range
straddling the crest, meeting at the waterline.

![The starting circle and the search window that confines the search](images/combo03_studio_circles.png)

Click **OK** on the Circles editor. The Inputs plot draws the model the runs
below are made on:

![The stability model: two zones, the reservoir load, and the starting circle](images/combo03_inputs.png){width=1000}

The blue arrows on the upstream face are the reservoir pressing on the slope,
marked **derived** because nothing entered them: **Water loads** on the `main`
sheet is `auto`, so the engine measures the pool against the ground and turns
what stands above it into a distributed load. The drawdown takes that load away,
and it is re-derived at every instant of the curve from the pool as it stood at
*that* moment. The red dashed arc is the starting circle, whose
center at (25, 44) sits above the top of the frame.

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

### The baseline: the slope under a full reservoir

The curve needs a starting point, and that is the factor of safety with the
reservoir full and the steady solution just made under it. Switch back to
**LEM** (`Ctrl+1`) and click **Run → Run LEM…** Leave every field at its
default. **Method** opens on **Spencer**, which satisfies both force
and moment equilibrium and is the method behind every factor of safety on this
page. **Analysis** opens on **Auto search**, which finds the run its own critical
circle, and **Number of slices** at 40.

Click **Run**. The Log reports the window before it reports anything else:

```text
Applying the file's search window: entry_range, exit_range, min_slip_depth.
Searching for the critical circular surface with SPENCER…
[🔁 iteration 1] center=(25.00, 44.00), FS=1.8820, grid=5.7995
[🔁 iteration 6] center=(23.55, 44.00), FS=1.8394, grid=1.4499
[🔁 iteration 12] center=(23.55, 45.81), FS=1.8308, grid=0.3625
[✅ converged] Iter=15, FS=1.8308 (ΔFS<0.0005) at (x=23.55, y=45.90, depth=5.63)
Critical FS = 1.831
Sliding mass = 4,733.6 kN/m over 46.09 m of failure surface
```

That first line confirms the window reached the run. Three of its limits are
named, the two ranges and the depth floor, and a search that reported none of
them would be searching the whole section.

**The factor of safety at full pool is 1.831**, on a circle centered at
(23.55, 45.90) with its tangent at elevation 5.63 — a surface through the
upstream shell that stops well above the rock.

![Spencer at full pool, with the reservoir load on the face](images/combo03_solution_full.png){width=1000}

The blue arrows are the reservoir, **176.6 kPa** at the heel under 18 m of water
and tapering to zero at the waterline. The pale blue band under the failure
surface is the pore pressure on each slice base, read off the steady field, and
the green hatched band above it is the effective normal stress the strength is
computed from. Both are what the drawdown is about to change.

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
  t=15: frame saved (steps so far=283, mass-balance closure=4.66e-03)
  t=30: frame saved (steps so far=716, mass-balance closure=1.82e-02)
  t=47: frame saved (steps so far=1207, mass-balance closure=2.89e-02)
  t=60: frame saved (steps so far=1579, mass-balance closure=2.67e-02)
  t=80: frame saved (steps so far=1802, mass-balance closure=3.99e-02)
  t=120: frame saved (steps so far=1940, mass-balance closure=4.92e-02)
  t=180: frame saved (steps so far=2050, mass-balance closure=3.25e-02)
  t=240: frame saved (steps so far=2076, mass-balance closure=1.59e-02)
  t=300: frame saved (steps so far=2084, mass-balance closure=9.54e-03)
  t=360: frame saved (steps so far=2098, mass-balance closure=6.84e-03)
Transient seepage complete — 12 saved frame(s).
```

**Twelve frames**, at t = 0, 2, 15, 30, 47, 60, 80, 120, 180, 240, 300 and 360 —
the union of the save-interval grid, the extra save times and the schedule's own
breakpoints, which
[SEEP-3 assembles](seep03_reservoir_drawdown.md#running-the-transient-march). Those
twelve are the whole of what a curve can be drawn through, and choosing them is a
modeling decision made before the stability question is asked.

### One instant at a time

The march has left twelve pore-pressure fields on the file, and the stability
run can read any one of them. Switch back to **LEM** and open
**Run → Run LEM…** once more. The form has grown
a group since the baseline run:

![Run LEM with the Seepage time group the march adds](images/combo03_studio_run_lem.png)

**Seepage time** appears only when a transient solution is loaded, and it names
the instant this run reads. Without it the choice would be invisible — it would
fall to wherever the results view's play bar happened to be sitting. Three ways
to name an instant, in the order they are reached for:

**Saved frame** lists the twelve the march stored, so the run starts
immediately. **Frame shown in the results viewer** names whatever frame the play
bar is sitting on, and is offered only while a transient results tab is open.
**Another time (re-marches the solution)** accepts any instant inside the run
duration, and serves it by re-running the march with that time added to the
save schedule — a full re-solve, because a field blended between two frames
solves nothing and is never invented.

Under the list, **Save as the model's stability time** writes the choice to the
`tseep` sheet, which is what makes a scripted or headless re-run read the same
frame. Left off, the choice governs this run only. The dialog opens on the last
saved frame, t = 360, and the note under the fields says why: with the model's
stability time blank, the last frame is what any run that does not choose one
will read.

The **Model checks** panel carries one note, the dependency stated back:

> Pore pressures come from the transient seepage solution, at t = 360 day. That
> frame is read into the model when the run starts, in place of any stored field.

Set **Saved frame** to `t = 30 day` and click **Run**:

```text
Applying the file's search window: entry_range, exit_range, min_slip_depth.
[✅ converged] Iter=20, FS=1.4962 (ΔFS<0.0005) at (x=7.46, y=54.67, depth=0.00)
Critical FS = 1.496
Sliding mass = 5,635.3 kN/m over 57.41 m of failure surface
```

That is one instant. Run LEM answers one at a time, which is what the **Seepage
time** group is for, and repeating it twelve times would draw the curve by hand.

#### The whole curve in one run

Click **Run → Parametric…** and set **Mode** to **Factor of safety vs time**:

![The Parametric dialog sweeping the march's saved frames](images/combo03_studio_parametric.png)

The parameter picker the other three modes use is gone, because nothing is
substituted here: every point solves *this* model against a different instant's
pore pressures. In its place is a **Saved frames** list holding the twelve the
march stored, all ticked — **All** and **None** set the whole list, and unticking
samples a long march, since each frame is a full stability run. Leave **Method** on
**Spencer** and **Number of slices** at 40 and **Re-search the critical
surface at each step** ticked, which is the right setting here because the
mechanism moves.

Click **Run**. Twelve searches, a few seconds each, stream their factors of safety
to the Log:

```text
Factor of safety vs time (lem): 12 instant(s) of the transient march, spencer, re-searching at each…
  t = 0 day      spencer   FS = 1.8308
  t = 2 day      spencer   FS = 1.8308
  t = 15 day     spencer   FS = 1.6073
  t = 30 day     spencer   FS = 1.4962
  ...
Lowest factor of safety 1.4962 at t = 30 day (12 instant(s), 0 without a result).
```

The run opens an **FS vs Time** tab with the curve, its lowest instant ringed and
labelled, and the reservoir schedule that drives it drawn faintly behind:

![The FS vs Time result tab](images/combo03_studio_fs_time.png)

**Day 30 reads 1.4962 here too** — the sweep applies the same search window off
the circles sheet that Run LEM applies, so one dialog run and twelve swept ones
land on the same number. An instant that produces no result comes back as a row
carrying its reason rather than as a gap in the curve.

From a script the same run is `xslope.sensitivity.fs_vs_time`:

```python
from xslope.sensitivity import fs_vs_time

ok, res = fs_vs_time(slope_data, solution, methods=("spencer",), search=True)
print(f"minimum FS = {res['min_fs']:.3f} at t = {res['critical_time']:g}")
```

Alongside the per-instant table it returns `critical_time` and `min_fs`.
[Factor of safety versus time](../parametric/sensitivity.md#factor-of-safety-versus-time)
carries the mode in full.

### The curve

| t (day) | pool (m) | FS | center | R (m) | tangent elevation (m) |
|:---:|:---:|:---:|:---:|:---:|:---:|
| 0 | 18.00 | 1.8308 | (23.55, 45.90) | 40.27 | 5.63 |
| 2 | 18.00 | 1.8308 | (23.55, 45.90) | 40.27 | 5.63 |
| 15 | 13.38 | 1.6073 | (18.99, 45.97) | 42.14 | 3.83 |
| 30 | 8.04 | **1.4962** | (7.46, 54.67) | 54.67 | 0.00 |
| 47 | 2.00 | 1.5760 | (10.97, 37.61) | 37.61 | 0.00 |
| 60 | 2.00 | 1.6484 | (8.50, 49.57) | 49.57 | 0.00 |
| 80 | 2.00 | 1.7056 | (7.65, 53.72) | 53.72 | 0.00 |
| 120 | 2.00 | 1.7457 | (7.52, 54.34) | 54.34 | 0.00 |
| 180 | 2.00 | 1.7630 | (7.67, 54.31) | 54.31 | 0.00 |
| 240 | 2.00 | 1.7663 | (7.67, 54.31) | 54.31 | 0.00 |
| 300 | 2.00 | 1.7675 | (7.67, 54.31) | 54.31 | 0.00 |
| 360 | 2.00 | 1.7680 | (7.67, 54.31) | 54.31 | 0.00 |

![The factor of safety at every saved instant, over the pool that drives it](images/combo03_curve.png){width=1000}

**Day 2 repeats day 0 exactly**, at 1.8308. The pool has not moved yet, so the
field has not either, and the two-day hold SEEP-3 entered to check the march
checks the stability side as well.

**The lowest of the twelve is 1.4962, on day 30** — 18% below full pool, and
**17 days before the drawdown ends**. The pool stands at elevation 8.04 there,
with 6 m of the 16 m drawdown still to come.

**Day 47, the end of the drawdown, reads 1.5760** — 0.08 *above* the minimum. The
instant the pool stops falling is not the instant the slope is weakest.

**The recovery is slow and does not quite finish.** From 1.5760 the curve climbs
to 1.6484 by day 60, 1.7056 by day 80, and then flattens: 1.7457, 1.7630, 1.7663,
1.7675, 1.7680 at days 120 through 360. Each step gains less than the one before,
which is the same diffusive relaxation the seepage outflow shows, and the last
three frames together move the answer by less than two thousandths. At 1.7680 the
empty dam still sits below the full-pool 1.8308: the reservoir was holding this
slope up, and taking it away costs more than draining the shell gives back.

**The mechanism grows.** At full pool the critical circle stops at elevation 5.63,
carrying 4,733.6 kN/m over 46.09 m of surface. From day 30 onward every critical
circle is tangent to the rock at elevation 0, and the day-30 surface carries
5,635.3 kN/m over 57.41 m — a fifth more weight on a quarter more surface. That
movement is why `search=True` matters: the same circle solved at every instant
would report the drop but would report it on a mechanism that only governs at one
end of the run.

#### The critical instant, drawn

![Spencer on day 30, mid drawdown](images/combo03_solution_min.png){width=1000}

Against the full-pool figure the difference sits in two places. The reservoir
load has fallen from **176.6 kPa** at the heel to **78.9 kPa**, and it reaches
only to elevation 8 on the face instead of 18 — both the pressure at its base and
the height it acts over are now under half what they were. Inside the dam the
contours are still crowded over the core, whose water table stands near elevation
13 with the pool 5 m below it, and the failure surface runs the full height of
the upstream slope, from the heel to the crest and tangent to the rock.

### Choosing the time steps

The curve's minimum sits on a saved frame, which raises the question the save
schedule always raises: whether the true minimum falls between two of them.

Answering it costs a re-march, and the Run LEM dialog's **Another time** option is
where a single instant is asked for. Asking for a set of them at once is
`fs_vs_time` with `times=` and the seepage data to re-march from, which serves the
whole set with **one** re-march before the first solve rather than one per
instant:

```python
ok, fine = fs_vs_time(slope_data, solution, times=[15, 20, 25, 30, 35, 40, 47],
                      methods=("spencer",), seep_data=seep_data, remarch=True)
```

| t (day) | pool (m) | FS |
|:---:|:---:|:---:|
| 15 | 13.38 | 1.6073 |
| 20 | 11.60 | 1.5640 |
| 25 | 9.82 | 1.5132 |
| 30 | 8.04 | **1.4959** |
| 35 | 6.27 | 1.5010 |
| 40 | 4.49 | 1.5166 |
| 47 | 2.00 | 1.5760 |

Five of those seven name no saved frame, so the whole set costs one full re-march
— the longest run on this page, made a second time. **The minimum comes back at
day 30, at 1.4959 against the saved grid's 1.4962** — a difference of 0.0002, on
a curve whose total swing is 0.33. The twelve frames had already resolved the dip.

Two readings come out of that. **The dip is broad**: the four instants from day 25
to day 40 all lie within 0.021 of the minimum, so a save schedule does not have to
be fine to find it — which is generally true of a drawdown, where the pool moves
over weeks. And **a
refinement pass is a check rather than a habit**: it re-solves the whole march to
move the answer, here, by less than the search's own convergence tolerance. The
time to spend a re-march is when the saved frames are far enough apart that the
minimum falls on the first or last of them, or when the curve's neighbors on
either side of the minimum are still falling steeply toward it.

### Reading the curve

Two things happen to the upstream slope when the pool is lowered, and they do not
happen at the same speed.

**The stabilizing load leaves immediately.** The reservoir presses on the face,
and that pressure resists the slide directly. It tracks the pool exactly — 176.6
kPa at the heel on day 0, 131.2 on day 15, 78.9 on day 30, and 19.6 from day 47
on, where all that stands against the heel is the 2 m of water left at tailwater
level — because a load is a boundary condition rather than something the soil has
to release.

**The pore pressure inside leaves slowly.** The shell can follow the pool: its
vertical conductivity over its drainable porosity is 0.5 / 0.22 = 2.3 m/day
against the pool's 0.36 m/day. The core cannot, at 0.005 / 0.03 = 0.17 m/day,
and it holds the head the full reservoir put into it long after the reservoir has
gone — 8.6 m of it still standing on day 47, which
[SEEP-3 measures](seep03_reservoir_drawdown.md#reading-the-frames).

The factor of safety falls while the first is ahead of the second and rises once
the second catches up. That crossover, not the end of the drawdown, is what
day 30 marks. Over the 17 days after it the pool gives up its last 6 m while the
shell behind the face drains through the elevations the failure surface runs at —
SEEP-3's shell node at (30.7, 8.9) goes from a pressure head of +0.92 m to −2.37 m
across that span, so the water table passes below it — and the strength recovers
faster than the remaining load falls. After day 47 the curve simply follows the
core emptying, which takes months.

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
be read at. On this dam a check staged at the end of the drawdown would be reading day 47,
where the slope is at **1.5760**; the worst of the twelve saved instants is day
30, at **1.4962**. The curve locates it, and once located it can be named as a
stage time and the drawdown check run there.

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
kind of curve on the same twelve frames so the two can be set side by side.

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
1,000 days and saves a frame every 200; and six extra save times — 5, 35, 50, 80,
150 and 300 — put frames where the regular grid would miss them. **Stage 1 time**
is `0` and **Stage 2 time** is `50`. Click **OK**.

Every point of the curve below reads its consolidation stresses at stage 1, and
stage 1 stays at t = 0 for the whole sweep. What moves from point to point is
stage 2.

### The mesh and the march it carries

The schedule states what the reservoir does; the saved frames hold what the dam
did about it. Still in **Seepage**, click the **Seep · Transient** tab. The play
bar steps through **twelve** instants — t = 0, 5, 35, 50, 80, 150, 200, 300, 400,
600, 800 and 1000 — and every pore pressure the runs below read comes from one of
them.

Those fields were solved on **2,080 nodes and 3,923 triangles**: linear triangles
auto-sized at 100 divisions across the 750 ft section, which puts the target
element size at 750/100 = 7.5 ft.
[COMBO-2 builds that mesh](combo02_rapid_drawdown.md#meshing-and-both-solves) with
**Run → Build Mesh…** and
[marches the pool on it](combo02_rapid_drawdown.md#marching-it) with
**Run → Run Seep…**, a run of several minutes on this dam. The file downloaded
above arrives past both, so nothing below waits on a seepage solve.

### Running the rapid drawdown sweep

Those twelve fields are what the sweep questions, one at a time. Switch back to
**LEM** (`Ctrl+1`) and click **Run → Parametric…**

![The Parametric dialog with the Rapid drawdown box ticked](images/combo03_rapid_parametric.png)

Set **Mode** to **Factor of safety vs time**. **Method** opens on **Spencer**,
which satisfies both force and moment equilibrium and is the method behind every
number in this part, and **Number of slices** opens at 40; leave both.

**Saved frames** lists the twelve this march stored, all ticked. Leave them so:
the whole march is what the comparison at the end of this part needs.

**Rapid drawdown at each time → ticked.** Ticking it changes what a ticked
instant means. Left off, each instant is one single-stage analysis of that
instant's water. Ticked, each instant becomes stage 2 of a three-stage drawdown
whose stage 1 is the march's initial state, so the curve answers how safe the
slope is if the pool falls from full to where it stands at that moment.

**Re-search the critical surface at each step** goes gray and is held on. A
drawdown's critical surface is not the drained one and it moves as the drawn-down
field changes, so there is nothing left for the toggle to decide: every point is
searched from the file's starting circle at (275, 235).

The note under the two boxes rewrites itself to match, and the **Model checks**
panel carries one warning, which repeats COMBO-2's: two of the three materials
carry no $d$ / $\psi$ and keep their drained strength through the drawdown.
Boundary set 2 raises nothing here, because
[COMBO-2 clears it](combo02_rapid_drawdown.md#clearing-boundary-set-2) before its
own march and this file arrives without it. A transient run solves boundary set 1
only, and both stages come out of that march.

Click **Run**. Eleven of the twelve instants are drawdowns and each is searched
in full, so the Log fills slowly; the table lands when the last search finishes:

```text
Rapid drawdown vs time (lem): 12 instant(s) of the transient march, spencer, re-searching at each…
Rapid drawdown versus time (stage 1 at t = 0)
t (day)  stage 1  stage 2       stage 3      FS  governs      Xo      Yo       R
      0        —        —             —       —        —       —       —       —   a drawdown from t = 0 to t = 0 is not a drawdown …
      5   1.5109   1.4563  not required  1.4563        2  242.74  256.34  165.90
     35   1.5453   1.0712  not required  1.0712        2  248.17  244.19  161.02
     50   1.5545   1.0164        1.0158  1.0158        3  243.93  244.90  163.64
     80   1.5512   1.0904  not required  1.0904        2  245.31  244.96  163.08
    150   1.5514   1.1524  not required  1.1524        2  245.31  243.53  161.82
    200   1.5514   1.1701  not required  1.1701        2  245.31  243.53  161.82
    300   1.5514   1.1848  not required  1.1848        2  245.31  243.53  161.82
    400   1.5514   1.1892  not required  1.1892        2  245.31  243.53  161.82
    600   1.5514   1.1921  not required  1.1921        2  245.31  243.53  161.82
    800   1.5514   1.1931  not required  1.1931        2  245.31  243.53  161.82
   1000   1.5514   1.1936  not required  1.1936        2  245.31  243.53  161.82
Lowest factor of safety 1.0158 at t = 50 day (12 instant(s), 1 without a result).
```

The t = 0 row's message is cut at the ellipsis to fit the page. In the Log it runs
on: *"a drawdown from t = 0 to t = 0 is not a drawdown -- the pool has not fallen
yet; stage 2 must be a later instant than stage 1."* At t = 0 stage 2 would be
stage 1, a fall from the full pool to itself, so that instant comes back as a row
carrying its reason rather than as a point on the curve. Eleven of the twelve
frames are drawdowns; the first is the state they all fall from.

### The rapid drawdown curve

The run opens a **Drawdown vs Time** tab — the same tab a single-stage sweep
uses, renamed for the curve it is holding:

![The rapid drawdown factor of safety at every saved instant, with its three stages behind it](images/combo03_rapid_curve.png)

The heavy line is the rapid drawdown factor of safety, the lower of stages 2 and
3, and the two thin dashed lines behind it are stages 1 and 2. Stage 2 lies under
the reported curve nearly everywhere, because it governs nearly everywhere, and
the dotted stage-3 line is drawn only at the one instant stage 3 was required.
The legend counts the instant that produced no result, so the missing t = 0 point
is never a silent gap.

**Day 50 reads 1.0158, and that is COMBO-2's answer.** COMBO-2's transient route
reports **1.016** on a circle at (243.93, 244.90) with a radius of 163.64 ft, and
the day-50 row above carries the same factor of safety on the same circle to the
hundredth of a foot. Nothing about that is a check of one page against another:
at t = 50 this sweep runs precisely the rapid drawdown COMBO-2 ran — stage 1 at 0,
stage 2 at 50, the same march, the same starting circle — so the curve passes
through COMBO-2's answer at day 50 because that answer is one of its points.

**Day 50 is also the only row where stage 3 runs.** Its **governs** column reads
3, and the reported 1.0158 is stage 3's against stage 2's 1.0164 — six
ten-thousandths. Every other row reads `not required`, which means no core slice
came out weaker drained than undrained and stage 2 is the answer.
[COMBO-2's sweep of the core's intercept](combo02_rapid_drawdown.md#the-governing-stage-flip)
finds why the margin at day 50 is that thin: the handover from stage 2 to stage 3
falls at $d = 223$ psf, 27 psf under the 250 psf this core carries.

**Day 5 costs 0.05 with no drawdown in it.** The pool holds at elevation 160
until day 5, so at that instant stage 2's water is stage 1's water — the same
field, the same reservoir load on the face. The reported 1.4563 still sits 0.055
below the 1.5109 stage 1 reads on the same circle, and the whole of that gap is
the undrained envelope being applied to the core. It is the price of the strength
substitution alone, measured with the water held still.

**The minimum, 1.0158, falls on day 50 — the instant the pool stops falling.**
That is not a rule. Two things race while the pool drops: the support the
reservoir gives the face keeps falling until the pool stops, and the pore
pressures in the core start dissipating from the moment the drop begins. On this
dam the drawdown is fast against a core at 0.001 ft/day, the pore pressures barely
move during the 45 days, and the worst instant is where the load is lowest — the
end. Part 1's dam drained fast enough for dissipation to overtake the load loss
mid-fall, and its minimum came 17 days *before* its drawdown ended. A slow enough
drawdown shows no dip at all. Which of those a given dam does is exactly what
the sweep is for; the instant cannot be named in advance.

The dip here is resolved only as finely as the saved frames: the nearest either
side of day 50 are day 35 at 1.0712 and day 80 at 1.0904, so the true minimum
could sit between 35 and 50. Part 1's [refinement pass](#choosing-the-time-steps)
— extra saved times around the end of the fall — is the way to pin it.

**The recovery flattens onto COMBO-2's two-steady answer.** From 1.0158 the curve
climbs to 1.0904 by day 80, 1.1524 by day 150, and then slows: 1.1701, 1.1848,
1.1892, 1.1921, 1.1931, 1.1936 at days 200 through 1000, the last three frames
together moving the answer by 0.0044. That plateau at 1.194 is
COMBO-2's **two-steady-solution** answer of **1.195**, reached from the other
direction. Set 2 solves the dam after the pore pressures have fully caught up
with the residual pool, and a march run to day 1000 arrives at very nearly the
same field, so the same drawdown asked against it gives very nearly the same
number.

**Stage 1 moves, and the water is not why.** Its column reads 1.5514 on the seven
rows from day 150 on, 1.5512 on day 80, 1.5545 on day 50, 1.5453 on day 35 and
1.5109 on day 5. Every one of those is the same physical state — full pool,
drained, the march's initial field — so the spread is the surface it was
evaluated on. Each row searches for the circle that minimizes the *drawdown*, and
reports stage 1 on whichever circle that turned out to be.

### The same march at drained strengths

The question Part 1 asks can be put to this dam too, on the same twelve frames,
by leaving the **Rapid drawdown at each time** box unticked. Open **Run → Parametric…** again, keep every
other field, untick **Rapid drawdown at each time**, and click **Run**.
**Re-search the critical surface at each step** comes back live and stays ticked,
which is the right setting for the same reason it was held on before.

Each point now runs about three times faster than its drawdown counterpart,
because it solves the section once rather than three times:

```text
Factor of safety vs time (lem): 12 instant(s) of the transient march, spencer, re-searching at each…
Factor of safety versus time
t (day)      FS      Xo      Yo       R
      0  1.5097  240.00  261.13  169.97
      5  1.5097  240.00  261.13  169.97
     35  1.0931  248.88  239.24  156.31
     50  1.0351  243.89  244.96  163.73
     80  1.1246  245.31  243.53  161.82
    150  1.2116  245.31  243.53  161.82
    200  1.2401  245.31  243.53  161.82
    300  1.2683  245.67  243.18  161.32
    400  1.2801  244.60  242.82  161.51
    600  1.2916  243.89  242.82  161.83
    800  1.2978  241.40  248.51  168.03
   1000  1.3015  241.40  248.51  168.03
```

Every instant produces a result here, including t = 0, which is an ordinary
analysis of the dam at full pool rather than a fall from itself. **Day 0 and day 5
read 1.5097 on the same circle**, because the pool has not moved between them and
neither has the field.

Three of these rows are COMBO-2's drained bracket runs, arrived at from a
different direction. **Day 0's 1.5097** is its full-pool drained figure of
**1.510**. **Day 50's 1.0351** is its day-50 drained figure of **1.035** — the
same frame, the same strengths. And **day 1000's 1.3015** approaches its
long-term drained figure of **1.311**, the gap being that a march stopped at day
1000 has not quite finished draining the core.

Studio's result tab holds one curve at a time, so the two are set on one pair of
axes here, drawn for this page. Running the sweep in both modes and plotting the
two results together reproduces it:

![The rapid drawdown curve and the single-stage curve of the same march](images/combo03_rapid_compare.png){width=1000}

**The rapid drawdown curve is the lower of the two at every instant**, and the distance
between them is not constant:

| t (day) | drawdown | single-stage | difference |
|:---:|---:|---:|---:|
| 5 | 1.4563 | 1.5097 | 0.0535 |
| 35 | 1.0712 | 1.0931 | 0.0219 |
| 50 | **1.0158** | **1.0351** | **0.0193** |
| 80 | 1.0904 | 1.1246 | 0.0342 |
| 150 | 1.1524 | 1.2116 | 0.0592 |
| 200 | 1.1701 | 1.2401 | 0.0700 |
| 300 | 1.1848 | 1.2683 | 0.0835 |
| 400 | 1.1892 | 1.2801 | 0.0909 |
| 600 | 1.1921 | 1.2916 | 0.0995 |
| 800 | 1.1931 | 1.2978 | 0.1046 |
| 1000 | 1.1936 | 1.3015 | 0.1079 |

**The gap is narrowest where the slope is weakest.** At day 50 the two analyses
are 0.019 apart, 1.9% of the rapid drawdown answer; by day 1000 they are 0.108 apart,
9.0% of it. The undrained envelope changes the answer least at the moment the
answer matters most.

The pore pressure in the core explains both halves of that. At day 50 the core
still holds the head the full reservoir put into it, so the
*drained* strength there is already cut by that pore pressure and lands close to
the undrained strength stage 2 assigns — the same reading
[COMBO-2 takes](combo02_rapid_drawdown.md#the-three-answers-and-what-brackets-them)
from its bracket table, where the same envelope moves one field by 0.12 and the
other by 0.02. Afterward the two analyses part company: the drained strength
recovers as the core drains, and the undrained strength cannot, because stage 2
computes it from stage-1 consolidation stresses that are the same at every point
of the curve. One curve climbs to 1.30 and the other stops at 1.19.

### Which curve to read

The two curves answer two questions, and neither is the other's conservative
version.

**Through the fall and just after it, read the rapid drawdown curve.** Its number for
this dam is **1.016, on day 50**. During those 45 days the core cannot shed the
reservoir's head, so the strength it can mobilize is the undrained one, and a
single-stage drained analysis at the same instant credits it with a strength it
does not have. That the single-stage curve reads 1.0351 there — only 0.019
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

**Long after the fall, read the single-stage curve.** By day 1000 the core has
drained, and 1.3015 is what the slope actually has. The rapid drawdown curve's 1.1936
at that instant is still answering "what if the pool fell from full to here",
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
- Why the curve needs a search window: the critical surface moves as the pore
  pressures change, and the limits on the circles sheet keep the movement inside
  one mechanism instead of letting it jump to the downstream slope or down to a
  surficial skin.
- The curve on Part 1's dam — 1.831 at full pool, a minimum of **1.496 on day
  30**, 1.576 at the end of the drawdown on day 47, and a slow recovery to 1.768
  by day 360 that never regains the full-pool value.
- That the saved times set the curve's resolution and cost a full re-march to
  refine: seven instants across the dip, five of them new, moved the minimum by
  0.0002 and left it on day 30.
- Why that minimum falls mid-drawdown: the reservoir load tracks the pool exactly
  while the pore pressure behind the face lags, and the slope is weakest where
  the gap between them is widest — 17 days before the pool stops moving.
- How to turn the same sweep into a rapid drawdown curve — **Rapid drawdown at each
  time** on the Parametric dialog makes every ticked instant stage 2 of a
  three-stage analysis whose stage 1 is the march's initial state, holds the
  re-search toggle on, and refuses the initial instant itself, which would be a
  fall from the full pool to itself.
- The rapid drawdown curve on COMBO-2's Johnson Reservoir dam — 1.4563 on day 5 with
  the pool still standing, a minimum of **1.0158 on day 50** where the fall ends,
  and a recovery that flattens at 1.1936 by day 1000. Day 50 is COMBO-2's own
  transient-route answer of 1.016, on the same circle, because that instant is
  the drawdown COMBO-2 ran.
- What the same twelve frames give at drained strengths — 1.5097 at full pool,
  1.0351 on day 50 and 1.3015 by day 1000 — and that the rapid drawdown curve is the
  lower of the two at every instant, by 0.019 at day 50 and 0.108 by day 1000.
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
