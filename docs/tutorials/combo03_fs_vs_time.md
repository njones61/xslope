---
title: "Tutorial COMBO-3 — Factor of Safety vs Time"
description: "The cored earth dam of Tutorial SEEP-3 given strengths and run for stability against every saved frame of its drawdown — how a transient seepage field reaches a limit equilibrium run one instant at a time, how the saved times set the curve's resolution, and why the lowest factor of safety falls 17 days before the reservoir stops moving."
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
reads none. This page assigns it, and everything else about the model comes
across unchanged.

<div class="tut-glance" markdown>
<div class="tgt-row">
<div class="tgt-tile"><span class="tg-label">Analysis</span><p>Seepage + limit equilibrium</p></div>
<div class="tgt-tile"><span class="tg-label">Open &amp; run</span><p>~30 min</p></div>
</div>
<div class="tgm-obj" markdown>
**Objectives** — Learn what a factor-of-safety-versus-time curve is and what it
adds to a single analysis; how a transient seepage field reaches a stability run
one instant at a time, and where that instant is chosen; how the saved times set
the curve's resolution and what a finer set of them costs; and how to read the
curve's minimum against the drawdown that caused it.
</div>
<p><span class="tg-pill">two materials</span><span class="tg-pill">transient seepage</span><span class="tg-pill">u = seep</span><span class="tg-pill">saved frames</span><span class="tg-pill">seepage time</span><span class="tg-pill">FS vs time</span><span class="tg-pill">Spencer</span><span class="tg-pill">circular search</span><span class="tg-pill">search window</span><span class="tg-pill">minimum slip depth</span><span class="tg-pill">automatic water loads</span></p>
<div class="tgm-model" markdown>
**Model** — [xslope_earth_dam_fs_time.xlsx](files/xslope_earth_dam_fs_time.xlsx),
SEEP-3's completed dam with a strength band on its materials table, a starting
circle on the upstream face and a search window. It carries no mesh and no
solution, so every run below is made from scratch
</div>
</div>

---

## The dam and the drawdown

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

### What the curve is made of

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

---

## The strengths this page adds

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

---

## The surface and the search window

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
reservoir never touches and a drawdown curve says nothing about.

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

---

## The mesh and the initial condition

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

---

## The baseline: the slope under a full reservoir

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

---

## The march

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

---

## One instant at a time

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

### The whole curve in one run

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

---

## The curve

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

### The critical instant, drawn

![Spencer on day 30, mid drawdown](images/combo03_solution_min.png){width=1000}

Against the full-pool figure the difference sits in two places. The reservoir
load has fallen from **176.6 kPa** at the heel to **78.9 kPa**, and it reaches
only to elevation 8 on the face instead of 18 — both the pressure at its base and
the height it acts over are now under half what they were. Inside the dam the
contours are still crowded over the core, whose water table stands near elevation
13 with the pool 5 m below it, and the failure surface runs the full height of
the upstream slope, from the heel to the crest and tangent to the rock.

---

## Choosing the time steps

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

---

## Reading the curve

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

### The relationship to a rapid-drawdown check

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
- The curve on this dam — 1.831 at full pool, a minimum of **1.496 on day 30**,
  1.576 at the end of the drawdown on day 47, and a slow recovery to 1.768 by
  day 360 that never regains the full-pool value.
- That the saved times set the curve's resolution and cost a full re-march to
  refine: seven instants across the dip, five of them new, moved the minimum by
  0.0002 and left it on day 30.
- Why the minimum falls mid-drawdown: the reservoir load tracks the pool exactly
  while the pore pressure behind the face lags, and the slope is weakest where
  the gap between them is widest — 17 days before the pool stops moving.

**Where to go next:** [SEEP-3](seep03_reservoir_drawdown.md) builds this dam's
transient seepage model from nothing and reads it as frames, head histories and a
water budget. [COMBO-2](combo02_rapid_drawdown.md) runs the three-stage
rapid-drawdown procedure on a larger dam, supplying its two water states three
ways.
[Factor of safety versus time](../parametric/sensitivity.md#factor-of-safety-versus-time)
carries the sweep mode itself — the times, the re-march rule, the per-instant row
contract and the FEM variant. [Transient seepage](../seep/transient.md) carries
the formulation and the boundary types, and
[Rapid drawdown analysis](../lem/rapid.md) the three-stage procedure. The
[tutorials index](index.md) lists the series.
