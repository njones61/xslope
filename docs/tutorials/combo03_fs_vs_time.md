---
title: "Tutorial COMBO-3 — Factor of Safety vs Time"
description: "Two dams run for stability against every saved frame of a transient drawdown — SEEP-3's earth dam at drained strengths, where the lowest factor of safety falls 12 days before the reservoir stops moving and the critical face changes from downstream to upstream and back, and COMBO-2's Johnson Reservoir dam as a three-stage rapid drawdown at each instant, where the curve passes through COMBO-2's own answer at day 50 and stays below the drained curve for the rest of the run."
---

# Tutorial COMBO-3 — Factor of Safety vs Time

A steady seepage analysis produces one pore-pressure field, so a stability run on
it has one answer. A transient analysis produces a field at every saved instant,
and the stability question can be asked at each of them. Plotting those answers
against time shows *when* the slope is weakest, not just how weak.

In [Tutorial SEEP-3](seep03_reservoir_drawdown.md) we built a small earth dam
with a granular shell and a clay core and lowered its reservoir 16 m over
45 days, stopping at the pore pressures. In **Part 1** we add strengths and run
one drained analysis per saved frame, searching both faces — a falling pool
weakens the upstream face, but a full reservoir loads it and leaves the
downstream slope weaker. In **Part 2** we repeat the exercise as a rapid
drawdown at every instant on [COMBO-2](combo02_rapid_drawdown.md)'s Johnson
Reservoir dam, whose clay core needs an undrained treatment.

<div class="tut-glance" markdown>
<div class="tgt-row">
<div class="tgt-tile"><span class="tg-label">Analysis</span><p>Seepage + limit equilibrium (drained, and rapid drawdown)</p></div>
<div class="tgt-tile"><span class="tg-label">Open &amp; run</span><p>~45 min</p></div>
</div>
<div class="tgm-obj" markdown>
**Objectives** — Run a stability analysis at every saved instant of a transient
seepage march, choose the saved times and the starting circles, read which face
governs, and repeat the run as a rapid drawdown for a dam with a clay core.
</div>
<p><span class="tg-pill">two materials</span><span class="tg-pill">three materials</span><span class="tg-pill">transient seepage</span><span class="tg-pill">u = seep</span><span class="tg-pill">saved frames</span><span class="tg-pill">seepage time</span><span class="tg-pill">FS vs time</span><span class="tg-pill">rapid drawdown</span><span class="tg-pill">three-stage procedure</span><span class="tg-pill">Kc = 1 envelope</span><span class="tg-pill">d and ψ</span><span class="tg-pill">stage times</span><span class="tg-pill">Spencer</span><span class="tg-pill">circular search</span><span class="tg-pill">starting circles</span><span class="tg-pill">minimum slip depth</span><span class="tg-pill">automatic water loads</span></p>
<div class="tgm-model" markdown>
**Part 1 model** — [xslope_earth_dam_fs_time.xlsx](files/xslope_earth_dam_fs_time.xlsx),
SEEP-3's dam with strengths on its materials table, a starting circle on each
face and a minimum slip depth. It carries no mesh and no solution, so we make
every run below from scratch

**Part 2 model** — [xslope_johnson_fs_time.xlsx](files/xslope_johnson_fs_time.xlsx),
COMBO-2's Johnson Reservoir dam — an undrained envelope on the core, `u = seep`
on all three zones, the pool schedule and both stage times. It ships meshed and
marched, with the twenty-one pore-pressure fields we sweep in Part 2 already on
it; we build both in [COMBO-2](combo02_rapid_drawdown.md)
</div>
</div>

---

## Part 1 — Drained analysis at each time step

We use the SEEP-3 dam for Part 1. Both of its zones carry drained strengths
(c′ and φ′), so at each saved time we run an ordinary effective-stress analysis
on the pore pressures from that time step.

### The dam and the drawdown

![The dam, its two zones, and the properties that make them behave differently](images/seep03_problem_sketch.png){width=1000}

The section is 110 m long and 22 m tall at the crest, on rock at elevation 0,
both faces sloping at about 2.3:1. A clay **core** runs from the rock to 4 m
below the crest, 17 m wide at its base and 8 m at its top; the granular **shell**
is everything else. The reservoir stands at elevation 18, the tailwater at
elevation 2 at the toe.

The pool holds at 18 for 2 days, falls 16 m to elevation 2 over the next 45 days,
and holds there for the remaining 253 days of a 300-day run. The shell follows
the pool down and the core does not, and the head the core keeps after the pool
has gone is what makes a drawdown a stability problem.
We built that model in [SEEP-3](seep03_reservoir_drawdown.md); only the save
schedule differs here.

#### How the curve is built

Each saved frame carries a pore pressure at every node of the mesh. A run with
`u = seep` reads that field at every slice base, so each frame gives one factor
of safety. We compute one for each of the nineteen frames and plot them against
time. Between one point and the next, only two things change: the pore
pressures and the reservoir load on the face.

This has two consequences. **The curve can only be as fine as the save
schedule**, because a time that is not a saved frame has no field to read. And
**the critical surface
moves** — the mechanism that governs under a full reservoir is not the one that
governs half drained, and on this dam not even on the same slope, so each instant
gets its own search.

### Adding the soil strength parameters

This dam was built for a seepage analysis in the [SEEP-3](seep03_reservoir_drawdown.md)
tutorial, so its materials carry conductivities and storage properties but no
strengths. A stability run needs
c′, φ′ and a unit weight for each zone. The workbook for this tutorial,
[xslope_earth_dam_fs_time.xlsx](files/xslope_earth_dam_fs_time.xlsx), is the
SEEP-3 model with those values added. Open it with
**File → Open…**, leave the mode strip (LEM | Seepage | FEM) on **LEM**, click
**Materials** in the Inputs tree, and set the **Show parameters for:** toggles to
**LEM** alone:

![The two zones with their strength band filled](images/combo03_studio_materials.png)

| mat | name | γ (kN/m³) | γsat (kN/m³) | option | c′ (kPa) | φ′ (deg) | u |
| :---: | --- | :---: | :---: | --- | :---: | :---: | --- |
| 1 | `shell` | 20 | 21 | `mc` | 0 | 32 | `seep` |
| 2 | `core` | 19 | 20 | `mc` | 10 | 25 | `seep` |

**These are typical values for a granular shell and a compacted clay core, chosen
for the exercise rather than measured on this dam.**

**Two unit weights per row.** γ applies above the water table and γsat below it,
and the slicer splits each slice's weight where the two meet, so the shell gets
lighter as the pool drains it.

**Both strengths are drained.** c′ and φ′ are effective-stress parameters, so
each slice base needs a pore pressure to form an effective normal stress from.

**`u` is `seep` on both rows.** That column decides where a slice base gets its
pore pressure, and `seep` sends it to a solved seepage field;
[COMBO-1](combo01_seepage_stability.md#the-column-that-connects-the-modes) covers
its four values. With a transient solution loaded, `seep` means one frame of it,
chosen at the run rather than here.

Click **OK**.

### Starting circles

Click **Circles**. The file carries two, one on each face of the dam:

| Xo | Yo | Option | Depth |
| :---: | :---: | --- | :---: |
| 7 | 56 | `Depth` | 0 |
| 103 | 59 | `Depth` | 0 |

There is one circle for each face, drawn on the deep mechanism that face can
make: its center sits beyond the heel or the toe, and **Depth** puts the bottom of the circle at elevation 0 — an
elevation rather than a distance, and 0 is the rock. The upstream circle runs
from (0.80, 0.34) to the crest edge at (51.50, 22.00), 57.63 m of surface; the
downstream one from (57.04, 22.00) to (109.19, 0.33), 58.88 m.

**Two circles because the critical face changes with the pool.** A drawdown
weakens the upstream slope, but a full reservoir *loads* it, and under 18 m of
water the downstream side is the weaker of the two. A curve from a full pool to
an empty one crosses both states.

**Both circles are drawn deep because a shallower seed refines into a smaller
mechanism.** A search steps away from its starting circle and stops where the
steps no longer lower the factor of safety, so a seed part way up the slope
settles near where it started.

Normally the critical circle in a zoned dam passes through the core, and the
starting circle should be drawn to cut it. This dam is small and its core is
narrow — 8 m wide at the top, with its crown 4 m below the crest — so the
critical circles stay mainly in the shell, and neither starting circle here cuts
the core.

Under the table, the **Search window** group holds ten limits on where a searched
surface may run; a blank field is a limit that is not applied, and
[the template page](../usage/input_template.md#search-window-optional) lists all
ten. The file sets one:

| Limit | Value | What it does |
| --- | :---: | --- |
| Min slip depth | 8 | The surface must reach at least 8 m below the ground surface |

That limit sets a floor on the *size* of the mechanism. With c′ = 0 the factor of
safety on a shallow surface parallel to the face does not depend on how deep it
is drawn, so an unconstrained search settles on a surficial wedge that ravels
rather than slides. 8 m on a 22 m dam keeps the search on embankment slides; the
Run LEM dialog offers the same limit as **Ignore surficial (skin) failures**.

The entry and exit ranges are left blank on purpose. Filling them would confine
the trace to one slope, and which slope governs at a given instant is what we
are measuring.

![Two starting circles and the one search limit the file sets](images/combo03_studio_circles.png)

Two circles with **Grid search** off raises `circles.multiple_without_grid` in
the **Model checks** panel on every run below, because a seeded search screens
each starting circle on one coarse grid and refines only the best-scoring one.
That is what we want here, since both circles sit on a real mechanism. Where
the critical face is not known in advance, **Grid search (auto-seed the circular
search)** sweeps a grid of centers against tangent elevations from the geometry
instead.

Click **OK** on the Circles editor. The Inputs plot draws the model:

![The stability model: two zones, the reservoir load, and the two starting circles](images/combo03_inputs.png){width=1000}

The blue arrows on the upstream face are the reservoir pressing on the slope,
marked **derived** because nothing entered them. **Water loads** on the `main`
sheet is `auto`, so the engine turns whatever pool stands above the ground into a
distributed load, re-derived at every instant from the pool as it stood then. The
two red dashed arcs are the starting circles, both centers sitting above the
frame.

### Building the mesh and solving the initial condition

Every pore pressure on this page comes from the seepage engine, so we build the
mesh and solve the full-pool steady state first. Switch the mode strip to
**Seepage** (`Ctrl+2`) and click **Run → Build Mesh…** Set **Element type** to
**Linear triangles (tri3)**, leave **Auto-size from geometry** ticked, and set
**Size divisions** to `64`. Click **Build**: **614 nodes and 1,089 triangles**,
the mesh [we built in SEEP-3](seep03_reservoir_drawdown.md#building-the-mesh).
Head is a scalar field, so linear elements are enough.

Click **Run → Run Seep…**, set **Run type** to **Steady**, leave **Convergence
tol** and **Max iterations** at `0.00010000` and `400`, and click **Run**. The
unconfined iteration converges in 59 sweeps to a total discharge of **0.16488
m³/day per m** of dam, the head running from 2.000 to 18.000 m and the pore
pressure from −78.6 to 176.6 kPa. The reservoir boundary is bound to the `pool`
series, read at t = 0, so the Log opens with
`Series 'pool' read at t = 0 for the steady solve: reservoir value 18`.

![The full-pool steady solution](images/seep03_steady.png){width=1000}

Nearly the whole 16 m head drop happens inside the core and the downstream slope
stays dry; we read this field in detail
[in SEEP-3](seep03_reservoir_drawdown.md#the-steady-solution-at-full-pool). The
march starts from it, and we run the baseline below against it.

### Baseline analysis at full pool

The curve starts from the factor of safety at full pool. Switch back to **LEM**
(`Ctrl+1`) and click **Run → Run LEM…** **Method** opens on **Spencer**, which
satisfies both force and moment equilibrium and is the method behind every factor
of safety on this page; **Analysis** opens on **Auto search**, which finds the
run its own critical circle, and **Number of slices** on 40. Leave every field
alone.

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

The first line names the one limit the file sets. The search screens both
starting circles, takes the downstream one, and steps its center down in halving
increments until the factor of safety stops moving by more than 0.0005.

**The factor of safety at full pool is 1.531, on the downstream face.** The
critical circle is centered at (103.00, 56.79), tangent to the rock at elevation
0, breaking at the crest at (58.12, 22.00) and daylighting at the toe at
(109.16, 0.34), carrying 5,580.9 kN/m of dam over 57.93 m of surface.

![Spencer at full pool: the critical circle is on the downstream face](images/combo03_solution_full.png){width=1000}

The reservoir is on the other side. Its blue arrows reach **176.6 kPa** at the
heel under 18 m of water and taper to zero at the waterline; that stabilizing
load is why the upstream face is not critical here. Under the failure surface the
pale blue band is the pore pressure on each slice base and the green hatched band
above it the effective normal stress.

### Running the transient seepage analysis

Now the pool comes down. Switch to **Seepage**, click **Run → Run Seep…** again,
set **Run type** to **Transient (time-dependent)** and click **Run**.
**Convergence tol** and **Max iterations** gray out — they belong to the steady
solve, and the march sets its own step size from how fast the field is moving.

The march solves its initial condition first, the field the steady run just
produced, then prints a line per saved frame:

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
100, 130, 180, 240 and 300 — all a curve can be drawn through, which makes the
save schedule a modeling decision taken before the stability question is asked.
It is uneven on purpose — five-day frames through the drawdown and widening
steps after — because the answer moves fastest while the pool is falling.

The **mass-balance closure** on each line is the gap between the stored-water
change and the net inflow, as a fraction of the flow passed so far. That flow is
small early on, so the ratio at t = 5 reads large without meaning much: the
absolute gap there is a third of a cubic meter per meter, and the field moves by
less than 0.04 m under a step size four times finer.

### Stability at one time step

The march has left nineteen pore-pressure fields on the file. Back in **LEM**,
reopen **Run → Run LEM…**; the form has grown a group:

![Run LEM with the Seepage time group the march adds](images/combo03_studio_run_lem.png)

**Seepage time** names the instant this run reads, three ways. **Saved frame**
lists the nineteen the march stored. **Frame shown in the results viewer** takes
whatever frame the play bar is on, and appears only while a transient results tab
is open. **Another time (re-marches the solution)** accepts any instant in the
run duration and re-runs the march with that time added to the save schedule; a
field blended between two frames is never invented.

**Save as the model's stability time** writes the choice to the `tseep` sheet, so
a scripted or headless re-run reads the same frame; left off, it governs this run
only. The dialog opens on the last saved frame, t = 300, which is what any run
reads while the model's stability time is blank.

Under the circles warning, the panel restates the dependency:

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
one. Run LEM answers one instant at a time, so we would need nineteen runs to
draw the whole curve this way.

#### Sweeping all time steps

A parametric sweep does that repetition itself. Click **Run → Parametric…** and
set **Mode** to **Factor of safety vs time**:

![The Parametric dialog sweeping the march's saved frames](images/combo03_studio_parametric.png)

Leave **Method** on **Spencer** and **Number of slices** at 40. The parameter
picker the other three modes use is gone, because nothing is substituted — every
point solves *this* model against a different instant's pore pressures. In its
place is a **Saved frames** list holding the nineteen the march stored, all
ticked, and unticking samples a long march. Leave **Rapid drawdown at each time**
unticked and **Grid search (auto-seed the circular search)** off, and leave
**Re-search the critical surface at each step** ticked, because the mechanism
moves.

Click **Run**. Nineteen searches report their factors of safety to the Log:

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
the pool schedule drawn faintly behind:

![The FS vs Time result tab](images/combo03_studio_fs_time.png)

**Day 35 reads 1.3313 here too, and day 0 reads the baseline's 1.5311** — the
sweep reads the same minimum slip depth and searches from the same two circles as
the dialog run. An instant that produces no result comes back as a row carrying
its reason rather than as a gap in the curve.

<!-- test: file=files/xslope_earth_dam_fs_time.xlsx, type=fs_vs_time, method=spencer, element_type=tri3, size_divisions=64, num_slices=40, expected_first=1.5311, critical_time=35, min_fs=1.3313, tolerance=0.005, benchmark=COMBO-3-drained -->

From a script we make the same run with `xslope.sensitivity.fs_vs_time`:

```python
from xslope.sensitivity import fs_vs_time

ok, res = fs_vs_time(slope_data, solution, methods=("spencer",), search=True)
print(f"minimum FS = {res['min_fs']:.3f} at t = {res['critical_time']:g}")
```

It returns `critical_time` and `min_fs` alongside a per-instant table whose rows
carry the `Xo`, `Yo` and `R` of the circle each instant settled on.
[Factor of safety versus time](../parametric/sensitivity.md#factor-of-safety-versus-time)
documents the mode in full.

### The factor of safety vs time curve

| t (day) | pool (m) | FS | center | R (m) | face |
| :---: | :---: | :---: | :---: | :---: | --- |
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

Every circle above bottoms at, or within a tenth of a meter of, the rock at
elevation 0. The **face** column is read off the center: the crest runs from
x = 51 to x = 59, so a center left of it is an upstream mechanism and a center
right of it a downstream one.

![The factor of safety at every saved instant, over the pool that drives it](images/combo03_curve.png){width=1000}

**Day 2 repeats day 0 exactly**, at 1.5311 on the same downstream circle, because
the pool has not moved and neither has the field.

**The lowest of the nineteen is 1.3313, on day 35** — 13% below the full-pool
1.5311, and **12 days before the drawdown ends**, with the pool at elevation 6.27
and 4.3 m of the 16 m drawdown still to come.

**Day 47, the end of the drawdown, reads 1.3813** — 0.05 *above* the minimum and
still upstream. The instant the pool stops falling is not the instant the dam is
weakest.

**After the fall the curve climbs, then flattens.** From 1.3813 it rises to
1.4386 by day 55, 1.4828 by day 65 and 1.5187 by day 80, all upstream; the answer
passes to the downstream mechanism at day 100 and levels out at 1.5482, 1.5516,
1.5550, 1.5566 and 1.5572 through day 300.

#### Which face governs

The face column does not move with the factor of safety.

**At full pool the downstream face governs, at 1.5311**, because 18 m of water
press the upstream slope into the dam.

**From the first frame of the drawdown the upstream face governs, and it holds
the answer for two thirds of the run.** By day 5 the pool has dropped about a
meter and the critical circle is already at (5.20, 65.90) at 1.5132; it stays
upstream through the fall, the minimum and the recovery, down to 1.3313 on day 35
and back to 1.5187 on day 80.

**At day 100 the downstream face takes it back**, at 1.5482, and keeps it to the
end, once the upstream slope has climbed past a downstream mechanism that has
barely moved all run.

So the reported curve is two mechanisms in sequence, not one surface moving, and
the figure marks each point by the face it came out on. The handover falls
between day 80 and day 100, where the two cross at about 1.55. A run given one
starting circle on one face would draw a smooth curve through one mechanism and
miss the other.

**The mechanism also grows.** At full pool the critical circle carries 5,580.9
kN/m over 57.93 m of surface; the day-35 circle carries 5,722.8 kN/m over
58.01 m of the other slope. Neither is the surface the other run would have
reported, which is why **Re-search the critical surface at each step** matters.

#### The critical time step

![Spencer on day 35, near the end of the drawdown](images/combo03_solution_min.png){width=1000}

The failure surface is on the other slope from the full-pool figure. The
reservoir load has fallen from **176.6 kPa** at the heel to **61.5 kPa** and
reaches only to elevation 6.3 instead of 18. The contours are still crowded over
the core, which holds head the pool no longer balances. The surface runs from the
heel at (0.80, 0.34) to the crest edge at (51.94, 22.00), tangent to the rock.

### Choosing the time steps

The minimum sits on a saved frame, which leaves open whether the true minimum
falls between two of them. `fs_vs_time` takes a set of unsaved instants
through `times=`, with the seepage data to re-march from, and serves the whole
set with **one** re-march before the first solve rather than one per instant:

```python
ok, fine = fs_vs_time(slope_data, solution, seep_data=seep_data, remarch=True,
                      times=[20, 22.5, 25, 27.5, 30, 32.5, 35, 37.5, 40],
                      methods=("spencer",))
```

Those nine halve the frame spacing across the dip, and four of them fall between
saved frames:

| t (day) | pool (m) | FS |
| :---: | :---: | :---: |
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

A coarse schedule is what makes that check necessary. SEEP-3's own twelve-frame
schedule saves nothing between day 30 and day 47, so a curve drawn through it
would report a minimum of 1.3344 on day 30 and never evaluate the lowest instant.
Re-march when the minimum falls on the first or last saved frame, or when the
points either side of it are still falling steeply.

### Reading the curve

Two things happen to the upstream slope when the pool is lowered, at different
speeds.

**The stabilizing load leaves immediately.** The reservoir's pressure on the face
resists the slide directly and tracks the pool exactly — 176.6 kPa at the heel on
day 0, 61.5 on day 35, and 19.6 from day 47 on — because a load is a boundary
condition, not something the soil has to release.

**The pore pressure inside leaves slowly.** The shell can follow the pool, its
vertical conductivity over drainable porosity being 0.5 / 0.22 = 2.3 m/day
against the pool's 0.36 m/day; the core cannot, at 0.005 / 0.03 = 0.17 m/day, and
still holds 8.6 m of head on day 47, as we measured
[in SEEP-3](seep03_reservoir_drawdown.md#reading-the-frames).

The factor of safety falls while the load is ahead of the drainage and rises once
the drainage catches up, so day 35 marks that crossover rather than the end of
the drawdown. Over the next 12 days the pool gives up its last 4.3 m while the
shell keeps draining, and the strength recovers faster than the load falls —
1.3436 by day 40, 1.3813 by day 47, 1.5187 by day 80.

#### Relationship to a rapid drawdown analysis

In [COMBO-2](combo02_rapid_drawdown.md) we analyze the same phenomenon with the
three-stage Duncan, Wright and Wong procedure, on a dam whose core carries an
undrained $K_c = 1$ envelope. That analysis names **two** instants, the slope
before the drawdown and the slope after it, and computes one factor of safety for
the second from consolidation stresses read at the first — two frames of a
transient march, named by `stage_1` and `stage_2` on the `tseep` sheet.

A curve cannot replace that check, and this model carries no $d$ / $\psi$ pair to
run one, but it says *which* instant the second stage should be read at. A check
staged at the end of the drawdown would read day 47, at **1.3813**, already 0.05
into the recovery; the worst of the nineteen saved instants is day 35, at
**1.3313**.

---

## Part 2 — Rapid drawdown analysis at each time step

The curve in Part 1 is a sequence of drained analyses, which suits a slope that
sheds pore water about as fast as the pool falls. A compacted clay core does not,
and the Part 1 model carries no $d$ / $\psi$ pair, so we cannot ask it the rapid
drawdown question at all.

[COMBO-2](combo02_rapid_drawdown.md)'s Johnson Reservoir dam can be. Its core
carries a $K_c = 1$ envelope, and in COMBO-2 we ran the three-stage procedure on
it three times, once per statement of where the water is, reading 1.181, 1.195
and 1.016 — the last from a transient march, stage 1 at t = 0 and stage 2 at
t = 50. Here we run it at **every** saved instant, then run Part 1's kind of
curve on the same twenty-one frames.

### The model

Download [xslope_johnson_fs_time.xlsx](files/xslope_johnson_fs_time.xlsx) and
open it with **File → Open…**, leaving the mode strip on **LEM**. The mesh and
the march are stored beside the workbook and Studio reads them on open, so the
dam arrives meshed and marched. In **LEM** the inputs view draws that mesh behind
the zones.

![The Johnson Reservoir dam as the shipped file opens, the mesh behind it](images/combo03_rapid_inputs.png){width=1000}

The section is 750 ft long, a 100 ft foundation on rock at elevation 0 with an
80 ft embankment on it, a sand **shell** on both faces and a compacted-clay
**core** carried 40 ft into the foundation as a cutoff key. The reservoir stands
at elevation 160 and the tailwater at 100, and the drawdown lowers the pool 50 ft
to 110. Three of the additions we made in COMBO-2 to
[SEEP-2](seep02_johnson_dam.md)'s geometry matter here.

**The undrained envelope.** Click **Materials** and set the **Show parameters
for:** toggles to **LEM** alone. Only the core carries a $K_c = 1$ envelope —
$d = 250$ psf and $\psi = 14°$, under its drained $c' = 400$ psf and
$\phi' = 18°$ — and blanks on the shell and the foundation declare them
free-draining through the drawdown.

| mat | name | γ (pcf) | c′ (psf) | φ′ (deg) | d (psf) | ψ (deg) | u |
| :---: | --- | :---: | :---: | :---: | :---: | :---: | --- |
| 1 | `shell` | 130 | 100 | 35 | — | — | `seep` |
| 2 | `core` | 125 | 400 | 18 | 250 | 14 | `seep` |
| 3 | `foundation` | 127 | 100 | 27 | — | — | `seep` |

We cover the procedure itself
[in COMBO-2](combo02_rapid_drawdown.md#the-three-stage-procedure) — what each
stage computes, and how stage 2's strength is interpolated between the undrained
and drained envelopes by the stress ratio each slice consolidated at.

**`u` is `seep` on all three rows**, so every slice base takes its pore pressure
from a solved seepage field. In COMBO-2 we deleted the file's two piezometric
lines when we replaced them, which is why the figure above draws a water surface
and a derived load and no piezometric line.

**The schedule and the stage times.** Switch the mode strip to **Seepage**
(`Ctrl+2`) and click **Transient** in the Inputs dock. The `pool` series holds at
elevation 160 for five days and falls to 110 over the following 45, the run lasts
500 days, and the save times are listed outright, a frame every five days to the
end of the fall and widening steps after. **Stage 1 time** is `0` and **Stage 2
time** is `50`. Click **OK**.

That schedule is this file's own; the [COMBO-2](combo02_rapid_drawdown.md)
workbook keeps a coarser one and a 1,000-day run. Its instants are packed where
the answer moves, and the run stops at day 500 because the curve has flattened
before then. Every point below reads its consolidation stresses at stage 1, which
stays at t = 0 for the whole sweep; what moves is stage 2.

### The mesh and transient solution shipped with the file

Still in **Seepage**, click the **Seep · Transient** tab. The play bar steps
through **twenty-one** instants — t = 0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50,
60, 70, 80, 100, 130, 170, 220, 300, 400 and 500 — and every pore pressure the
runs below read comes from one of them.

Those fields were solved on **2,080 nodes and 3,923 triangles**, linear triangles
auto-sized at 100 divisions across the 750 ft section, a target element size of
750/100 = 7.5 ft. We built
[that mesh in COMBO-2](combo02_rapid_drawdown.md#meshing-and-both-solves) and
[marched the pool on it](combo02_rapid_drawdown.md#marching-it); the file
downloaded above arrives past both.

### Running the rapid drawdown sweep

Switch back to **LEM** (`Ctrl+1`) and click **Run → Parametric…**

![The Parametric dialog with the Rapid drawdown box ticked](images/combo03_rapid_parametric.png)

Set **Mode** to **Factor of safety vs time**. **Method** opens on **Spencer** and
**Number of slices** on 40; leave both. **Saved frames** lists the twenty-one this
march stored, all ticked — leave them, because the comparison at the end of this
part needs the whole march.

**Rapid drawdown at each time → ticked.** Left off, each instant is one
single-stage analysis of that instant's water. Ticked, each instant becomes stage
2 of a three-stage drawdown whose stage 1 is the march's initial state, so the
curve gives the factor of safety if the pool falls from full to where it stands
at that moment.

**Re-search the critical surface at each step** goes gray and is held on, because
a drawdown's critical surface is not the drained one and moves as the field does.
Every point is searched from the file's starting circle at (275, 235), on the
upstream face — lowering the pool takes the water off that slope and leaves it
inside the embankment.

The **Model checks** panel repeats COMBO-2's one warning: two of the three
materials carry no $d$ / $\psi$ and keep their drained strength through the
drawdown. Boundary set 2 raises nothing, because
[we cleared it in COMBO-2](combo02_rapid_drawdown.md#clearing-boundary-set-2)
before that march; a transient run solves boundary set 1 only.

Click **Run**. Twenty of the twenty-one instants are drawdowns, each searched in
full, so the table lands when the last search finishes:

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

The t = 0 row's message is cut at the ellipsis to fit the page. In the Log it
runs on: *"a drawdown from t = 0 to t = 0 is not a drawdown -- the pool has not
fallen yet; stage 2 must be a later instant than stage 1."* The first frame is
the state the other twenty fall from.

### The rapid drawdown curve

The run opens a **Drawdown vs Time** tab, the same tab a single-stage sweep
uses:

![The rapid drawdown factor of safety at every saved instant, over the pool schedule that drives it](images/combo03_rapid_curve.png)

One curve is drawn, carrying the rapid drawdown factor of safety at each instant,
the lower of stages 2 and 3; the three stages stay in the table above, where the
governing one is named per row. Behind it, on the right axis, runs the pool
schedule. The lowest instant is ringed and labeled, the red guide marks FS = 1,
and the legend counts the instant that produced no result.

**Day 50 reads 1.0157, the answer we got in COMBO-2.** The transient route there
reports **1.0158** on a circle at (243.93, 244.90) with a radius of 163.64 ft,
and the day-50 row lands on that circle at 1.0157. At t = 50 this sweep poses the
drawdown we posed in COMBO-2 — stage 1 at 0, stage 2 at 50, the same march and
the same starting circle — so the curve passes through that answer.

<!-- test: file=files/xslope_johnson_fs_time.xlsx, type=fs_vs_time, method=spencer, rapid=true, march=file, num_slices=40, expected_first=1.4563, critical_time=50, min_fs=1.0157, tolerance=0.005, benchmark=COMBO-3-rapid -->

**Stage 3 runs on three rows, all of them at the bottom of the dip.** Days 45, 50
and 60 read 3 in the **governs** column, by a single ten-thousandth on days 45
and 60 and six on day 50 — 1.0279 / 1.0278, 1.0163 / 1.0157, 1.0523 / 1.0522.
Every other row reads `not required`, meaning no core slice came out weaker
drained than undrained. The margins are thin because the handover falls at
$d = 223$ psf, 27 psf under the 250 psf this core carries.

**Day 5 costs 0.05 with no drawdown in it.** The pool holds at elevation 160
until day 5, so stage 2's water is stage 1's water. The reported
1.4563 still sits 0.055 below the 1.5109 stage 1 reads on the same circle, and
that gap is the undrained envelope alone, measured with the water held still.

**The minimum, 1.0157, falls on day 50, the instant the pool stops falling.**
That is not a rule. Two effects compete while the pool drops: the reservoir's
support keeps falling until the pool stops, and the core's pore pressures start
dissipating from the moment the drop begins. On this dam the drawdown is fast
against a core at 0.001 ft/day, so the pore pressures barely move and the worst
instant is the one with the least load. Part 1's dam drained fast enough for
dissipation to overtake the load loss mid-fall, and its minimum came 12 days
*before* its drawdown ended.

Day 45 reads 1.0278 and day 60 reads 1.0522, so the saved times bracket that
minimum and the curve turns at a frame rather than between two.

**The recovery climbs toward COMBO-2's two-steady answer.** From 1.0157 the curve
rises to 1.0902 by day 80 and 1.1414 by day 130, then slows to 1.1612, 1.1744,
1.1848, 1.1893 and 1.1912 through day 500, the last two frames together moving it
by 0.0019. COMBO-2's **two-steady-solution**
answer is **1.195**, and the curve is still 0.004 under it at day 500, because by
then the core has not quite finished draining.

**Stage 1 moves, and the water is not why.** Its column reads 1.5514 on the seven
rows from day 100 on, 1.5512 on day 80, 1.5545 on day 50, 1.5453 on day 35 and
1.5109 on day 5, all of them the same physical state. Each row searches for the
circle that minimizes the *drawdown* and reports stage 1 on that circle.

### The same time steps with drained strengths

We can run the same twenty-one frames at drained strengths. Open
**Run → Parametric…** again, untick **Rapid drawdown at each time**, keep every
other field, and click **Run**; **Re-search the critical surface at each step**
comes back live and stays ticked.

Each point now runs about three times faster, solving the section once rather
than three times:

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

Every instant produces a result here, including t = 0, an ordinary full-pool
analysis rather than a fall from itself. **Day 0 and day 5 read 1.5097 on the
same circle**, the pool not having moved between them.

Two of these rows repeat the drained bracket runs from COMBO-2: **day 0's
1.5097** is the full-pool drained figure of **1.510** there, and **day 50's
1.0350** the day-50 figure of **1.035**. The third of those, the long-term
**1.311**, is not on this curve — day 500 reads 1.2870 and is still climbing,
because a march of this length leaves the core part drained.

<!-- test: file=files/xslope_johnson_fs_time.xlsx, type=fs_vs_time, method=spencer, march=file, num_slices=40, times=0;50;500, expected=1.5097;1.0350;1.2870, tolerance=0.005, benchmark=COMBO-3-single -->

Studio's result tab holds one curve at a time, so the figure below plots both
results on one pair of axes:

![The rapid drawdown curve and the single-stage curve of the same march](images/combo03_rapid_compare.png){width=1000}

**The rapid drawdown curve runs below the single-stage curve at every instant**,
and the distance between them is not constant:

| t (day) | drawdown | single-stage | difference |
| :---: | :---: | :---: | :---: |
| 5 | 1.4563 | 1.5097 | 0.0534 |
| 10 | 1.3496 | 1.3944 | 0.0448 |
| 15 | 1.2704 | 1.3093 | 0.0389 |
| 20 | 1.2036 | 1.2374 | 0.0338 |
| 25 | 1.1489 | 1.1779 | 0.0290 |
| 30 | 1.1051 | 1.1301 | 0.0250 |
| 35 | 1.0711 | 1.0930 | 0.0219 |
| 40 | 1.0457 | 1.0656 | 0.0199 |
| 45 | 1.0278 | 1.0469 | 0.0191 |
| 50 | **1.0157** | **1.0350** | **0.0193** |
| 60 | 1.0522 | 1.0765 | 0.0243 |
| 70 | 1.0743 | 1.1037 | 0.0294 |
| 80 | 1.0902 | 1.1244 | 0.0342 |
| 100 | 1.1159 | 1.1587 | 0.0428 |
| 130 | 1.1414 | 1.1949 | 0.0535 |
| 170 | 1.1612 | 1.2252 | 0.0640 |
| 220 | 1.1744 | 1.2475 | 0.0731 |
| 300 | 1.1848 | 1.2684 | 0.0836 |
| 400 | 1.1893 | 1.2802 | 0.0909 |
| 500 | 1.1912 | 1.2870 | 0.0958 |

**The gap is narrowest where the slope is weakest.** Across the bottom of the dip
— days 40 to 50 — the two analyses run 0.019 apart, about 1.9% of the rapid
drawdown answer; by day 500 they are 0.096 apart, 8.0% of it.

The pore pressure in the core explains both halves. At day 50 the core still
holds the head the full reservoir put into it, so the *drained* strength there is
already cut by that pore pressure and lands close to the undrained strength stage
2 assigns, the same reading
[we take in COMBO-2](combo02_rapid_drawdown.md#the-three-answers-and-what-brackets-them)
from its bracket table, where the same envelope moves one field by 0.12 and the
other by 0.02. Afterward the two separate: the drained strength recovers as the
core drains, and the undrained one cannot, because stage 2 computes it from
stage-1 consolidation stresses that never change. One curve climbs to 1.29 and
the other levels off near 1.19.

### Which curve to use

The two curves answer two questions, and neither is the other's conservative
version.

**Through the fall and just after it, read the rapid drawdown curve.** Its number
for this dam is **1.016, on day 50**. During those 45 days the core cannot shed
the reservoir's head, so the strength it can mobilize is the undrained one, and a
single-stage drained analysis credits it with strength it does not have. That the
single-stage curve reads only 0.019 higher is a property of this dam — on this
core the two strengths sit close together at the stresses the surface sees. A
more contractive core or a faster drawdown would open that gap well beyond
0.02.

**Long after the fall, read the single-stage curve.** By day 500 the core has
largely drained and 1.2870 is what the slope actually has. The rapid drawdown
curve's 1.1912 at that instant still answers "what if the pool fell from full to
here", the long-term drawdown state COMBO-2's two steady solutions compute, which
applies only if the reservoir is filled and drawn down again.

For a zoned dam with a core that does not drain, the design number is the rapid
drawdown curve's minimum and the long-term number is the single-stage curve's
plateau.

---

## Conclusion

In this tutorial we covered:

- What a factor-of-safety-versus-time curve is — one stability analysis per saved
  frame of a transient seepage solution, with nothing but the water changing
  between the points.
- How a frame reaches a stability run — `u = seep` sends every slice base to a
  solved field, Run LEM's **Seepage time** group names which frame, and
  **Run → Parametric… → Factor of safety vs time** sweeps all of them at once.
- Where to put a starting circle when either face may govern — one per face, each
  on the deep mechanism its slope can make, because a seeded search refines only
  the best-screening circle.
- The curve on Part 1's dam — 1.531 at full pool on the **downstream** face, a
  minimum of **1.331 on day 35** on the **upstream** face, 1.381 at the end of the
  drawdown on day 47, and a climb to 1.557 by day 300. The minimum falls
  mid-drawdown, 12 days before the pool stops moving, because the reservoir load
  tracks the pool exactly while the pore pressure behind the face lags; halving
  the frame spacing across the dip moved it by 0.0001.
- How to turn the same sweep into a rapid drawdown curve — **Rapid drawdown at
  each time** makes every ticked instant stage 2 of a three-stage analysis whose
  stage 1 is the march's initial state.
- The rapid drawdown curve on COMBO-2's Johnson Reservoir dam — a minimum of
  **1.0157 on day 50** where the fall ends and a recovery to 1.1912 by day 500;
  the same frames at drained strengths give 1.5097, 1.0350 and 1.2870.
- Which curve to read — the rapid drawdown curve through the fall and just after
  it, the single-stage curve for the long-term condition.

**Where to go next:** in [SEEP-3](seep03_reservoir_drawdown.md) we build the
transient seepage model of Part 1's dam; in
[COMBO-2](combo02_rapid_drawdown.md) we build Part 2's and run the three-stage
procedure at one instant, on the dam we build in
[SEEP-2](seep02_johnson_dam.md).
[Factor of safety versus time](../parametric/sensitivity.md#factor-of-safety-versus-time)
carries the sweep mode itself, [Transient seepage](../seep/transient.md) the
formulation and the boundary types, and
[Rapid drawdown analysis](../lem/rapid.md) the three-stage procedure. The
[tutorials index](index.md) lists the series.
