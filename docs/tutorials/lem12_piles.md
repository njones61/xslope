---
title: "Tutorial LEM-12 — Piles"
description: "A 20 ft clay slope stabilized by two rows of drilled shafts in XSLOPE: the pile force left blank so Ito & Matsui computes it from the diameter and the spacing at every trial surface, the moment capacity that caps it, what the spacing is worth, and the shallow surface that slides over the pile row."
---

# Tutorial LEM-12 — Piles

This tutorial shows how to model stabilizing piles in a limit equilibrium
analysis. The example is a 20 ft slope in medium-stiff clay that stands at a
factor of safety of 1.149 on its own; two rows of 2 ft drilled shafts
through the face bring it to 1.842. The pile force is not entered anywhere
in the model. XSLOPE computes it instead — the lateral force the soil can
push onto each pile, from the pile diameter and the center-to-center
spacing, by the Ito & Matsui (1975) method — and recomputes it for every
trial surface the search tries. This is the second of
the two ways XSLOPE models a pile: the tieback wall tutorial
([LEM-9](lem09_tieback_wall.md)) entered its soldier pile the other way,
with the force stated directly, and this page compares the two routes.


<div class="tut-glance" markdown>
<div class="tgt-row">
<div class="tgt-tile"><span class="tg-label">Analysis</span><p>Limit equilibrium</p></div>
<div class="tgt-tile"><span class="tg-label">Open &amp; explore</span><p>~15 min</p></div>
</div>
<div class="tgm-obj" markdown>
**Objectives** — Learn how to model stabilizing piles: how XSLOPE computes the
lateral force by Ito & Matsui, what pile spacing and structural capacity do to
it, and how to check for the shallow surface that bypasses the row.
</div>
<p><span class="tg-pill">one material</span><span class="tg-pill">piles</span><span class="tg-pill">Ito &amp; Matsui</span><span class="tg-pill">pile spacing</span><span class="tg-pill">structural capacity</span><span class="tg-pill">specified pile force</span><span class="tg-pill">circular search</span><span class="tg-pill">grid seeding</span></p>
<div class="tgm-model" markdown>**Completed model** — [xslope_piles.xlsx](../lem/files/xslope_piles.xlsx), the same file used by [LEM Sample Problem 10](../lem/samples.md#10-slope-stabilized-with-piles)</div>
</div>

---

## The slope

The slope itself is kept simple so that everything interesting on this page
comes from the piles. It is a single soil — a medium-stiff clay with γ = 120 pcf,
c = 200 psf, and φ = 20° — over a rigid
base 10 ft below the toe. The face rises 20 ft at 1:1 from (0, 0) to (20, 20),
with level ground either side. There is no water — the clay's **u** is `none` —
so the section is analyzed dry and the strength on every slice base is the
Mohr-Coulomb pair above.

Two rows of drilled shafts stand vertically through the face, 2 ft in
diameter at 6 ft centers, each row running along the length of the slope and
taken down to the rigid base. The 2D section models each row per foot of
slope — a pair of isolated piles would be a three-dimensional problem and
does not belong in a tool like this. The upper row is at
x = 10, where the ground is at elevation 10; the lower at x = 5, where it is at
elevation 5.

### Opening the model

Download [xslope_piles.xlsx](../lem/files/xslope_piles.xlsx) and open it in
Studio — **File → Open**. The Inputs plot draws the section, the starting circle
the file carries, and the two pile rows as green bars:

![The loaded model](images/lem12_inputs.png){width=1000}

The piles carry no force label, because on this model the force is not an input.
Where a pile states its own force — the soldier pile of
[LEM-9](lem09_tieback_wall.md#running-the-analysis) does — the Inputs plot prints
it beside the bar.

### The pile rows

Everything XSLOPE needs to know about the two pile rows sits in two rows of
the piles sheet. Read them closely before running anything. Open
**Piles** in the **Inputs** tree and press **Table view**. Its columns match
the piles worksheet; with **Show parameters for:** set to **LEM**, the
columns only the finite element engine reads (`E`, `I`, `Area`, `Fixity`)
are hidden:

![The two pile rows as the file carries them](images/lem12_studio_piles_table.png)

Each row of the sheet describes one pile row. The columns after the
endpoints control how the pile force is computed, limited, and applied:

| Label | x1 | y1 | x2 | y2 | H | Appl | D | S | Vcap | Mcap |
|---|:---:|:---:|:---:|:---:|:---:|---|:---:|:---:|:---:|:---:|
| `lower row` | 5 | 5 | 5 | -10 | | Active | 2 | 6 | 46000 | 60000 |
| `upper row` | 10 | 10 | 10 | -10 | | Active | 2 | 6 | 46000 | 60000 |

**H is empty on both rows.** `H` is the pile force per foot of slope — the
lateral resistance a row contributes to the equilibrium equations. Leaving
it blank while giving a diameter `D` and a spacing `S` tells XSLOPE to
compute that force by the Ito & Matsui method; entering a number in `H`
uses that number instead.
The direction of the pile force is not entered: XSLOPE takes it
perpendicular to the pile's own axis, which for these vertical shafts means
horizontal — the usual case for a stabilizing pile. `Appl` is `Active`: the force enters the equilibrium equations as it
stands. The alternative, `Passive`, treats the force as a resistance and
divides it by the factor of safety, so the support carries the same margin
as the soil strength. (A file that leaves the cell blank is read as
active; files XSLOPE saves write the choice out.)

**Vcap and Mcap are properties of one shaft**, not of a foot of slope: 46,000 lb
of shear capacity and 60,000 ft·lb of moment capacity, consistent with a 2 ft
reinforced concrete section at f′<sub>c</sub> = 4,000 psi. They bound whatever
the soil calculation produces, and on this model they bind hard.

Press **List view** to read one pile row at a time, as a form in four
groups —
**Identity**, **Geometry**, **Capacity / design**, **Behavior** — with the
preview drawing the selected pile bold on the section:

![The piles editor on the lower row](images/lem12_studio_piles.png){width=1000}

The **Capacity / design** group puts the empty **H** directly above the **D** and
**S** it will be computed from, and the two capacity fields carry their own
per-shaft units in their labels: **Vcap (per element, lb)** and
**Mcap (per element, lb·ft)**.

---

## Running the analysis

With the inputs understood, run the search and see what the slope does
with the pile rows in place. Click **Run LEM…** and choose **Method** = `Spencer` and **Analysis** =
`Auto search`, with the slice count left at 40:

![The Run LEM dialog on the loaded model](images/lem12_studio_run_lem.png)

Click **Run**. The search walks the circle the file carries down onto the rigid
base, and the plot draws every circle it tried in gray with the one it kept in
red:

![The Spencer search with the pile rows in place](images/lem12_search.png){width=1000}

**FS = 1.842**, on a circle centered at (3.79, 26.69) with a radius of 31.79,
tangent at elevation −5.10. The solution plot draws that surface with its base
stresses and marks the two points where it cuts the piles:

![Spencer on the surface the search found](images/lem12_solution.png){width=1000}

The surface enters 13.5 ft out in front of the toe, exits 14.9 ft behind the
crest and carries 62,198 lb/ft of clay on 61.27 ft of base. The two red dots are
the pile crossings, at elevation −5.07 on the lower shaft and −4.48 on the upper
one — 10.07 ft and 14.48 ft of soil above the surface at each pile, which is the
depth the force is computed over.

The Log pane prints the computed forces with the answer — one line per pile
row, running from the soil force through the capacity that governed it to
the force each row applies:

```text
Critical FS = 1.842
Sliding mass = 62,198.1 lb/ft over 61.27 ft of failure surface
Pile row 'lower row': Ito & Matsui soil force = 44,178 lb per pile; bending governs (Mcap/Lm, Lm = 3.94) -> 15,244 lb per pile = 2,540.7 lb/ft applied
Pile row 'upper row': Ito & Matsui soil force = 81,729 lb per pile; bending governs (Mcap/Lm, Lm = 5.47) -> 10,962 lb per pile = 1,827.0 lb/ft applied
Total pile resistance = 4,367.6 lb/ft
```

The **Analysis Report** (**File → Generate Report…**) carries the same
numbers in full: the pile inputs echoed, with `H` reading *computed* rather
than a number, and the slice table's **H<sub>p</sub>** column — *pile
resistance mobilized at the slice base, per unit thickness* — placing
**2,540.7** and **1,827.0 lb/ft** on the two slices the rows cross
([Analysis Report](../studio/reports.md)).

Each piece of those Log lines gets its own section below: where the soil
force comes from, and what "bending governs" means for the shafts.

### What the two rows are worth

FS = 1.842 describes the reinforced slope, but it does not say what the
piles contributed. Measuring that takes two comparisons: first the slope
without any piles, then the same surface with each row removed in turn.

Taking the piles out and searching the same slope again — same clay, same
geometry, same starting circle, a separate search with its own critical surface —
gives the slope on its own:

![Spencer on the same slope with no piles](images/lem12_solution_nopiles.png){width=1000}

**FS = 1.149** on a shallower circle from the toe, tangent at elevation −0.04,
moving 21,166 lb/ft against the 62,198 on the surface found with the pile
rows in place. The piles raise the
factor of safety by 0.69, and they also change the failure mechanism: the surface the
search settles on with them present is deeper, longer and three times the mass,
because the shallower one now has two shafts across it.

To separate the two rows' contributions, hold the surface still. On the
critical circle from the search with both rows present — **one surface, no search, only which piles are present
changing** — each row can be removed on its own:

| Piles present | Spencer FS |
|---|:---:|
| Neither | 1.481 |
| Lower row only | 1.675 |
| Upper row only | 1.613 |
| Both | 1.842 |

The lower row is worth slightly more than the upper one on this surface even
though it is the shallower of the two, because its crossing sits nearer the toe
where the base is flatter and the force resolves more directly against sliding.

---

## How the force is computed

The console log and the report say what the piles delivered. This section
explains where those
numbers come from, because understanding the calculation is what makes the
rest of the page — the capacity limits, the spacing study — predictable
rather than mysterious.

Ito & Matsui treat the soil between two adjacent piles as squeezing plastically
through the gap between them, and derive from Mohr-Coulomb plasticity the lateral
pressure that flow puts on each pile. The pressure at depth z is
p(z) = c·A₁ + γz·A₂, with A₁ and A₂ arching coefficients that depend only on the
friction angle and on the ratio of the center-to-center spacing S to the clear
gap S − D. On this clay at φ = 20° with D = 2 and S = 6, they are **A₁ = 7.569
ft** and **A₂ = 4.755 ft**.

The force on one pile is that pressure integrated from the ground surface down to
the failure surface, and the per-foot-of-slope force is that divided by the
spacing. **Both depend on where the trial surface is**, so both are recomputed at
every surface a search evaluates — a deeper surface means a taller soil column
pushing on the pile. On the critical surface above:

| | Lower row | Upper row |
|---|:---:|:---:|
| Depth to the surface, z<sub>f</sub> (ft) | 10.07 | 14.48 |
| Soil force per pile, F<sub>pile</sub> (lb) | 44,178 | 81,729 |
| Soil force per foot of slope, F<sub>pile</sub>/S (lb/ft) | 7,363.0 | 13,621.5 |

The deeper the surface, the larger the force: the upper row, whose head
stands 5 ft higher on the
face, has 44% more soil above the surface and develops 85% more force per shaft.
[Piles and Concrete Piers](../lem/piles.md#ito-matsui-1975-theory)
gives the derivation, the coefficients in full and the φ = 0 and c = 0 special
cases; the method is applicable for S/D between about **2 and 8**, and this
model's S/D = 3 sits in the middle of that band.

These are the forces the soil can deliver, not the forces the run applied.
The Log lines show 2,540.7 and 1,827.0 lb/ft — a fraction of the 7,363.0
and 13,621.5 computed above. The difference is the structural capacity of
the shafts, which the next section covers.

---

## What the structural capacity does

The Ito & Matsui force is the soil's capacity to push, not the pile's capacity to
resist. `Vcap` and `Mcap` are the second half of the answer: the force actually
used is the lesser of the soil force, the shear capacity, and the force the
moment capacity permits, M<sub>cap</sub>/L<sub>m</sub>, where L<sub>m</sub> is the
arm from the centroid of the computed pressure distribution down to the failure
surface.

On the critical surface both shafts are governed by bending:

| | Lower row | Upper row |
|---|:---:|:---:|
| Soil force per pile (lb) | 44,178 | 81,729 |
| Shear limit, V<sub>cap</sub> (lb) | 46,000 | 46,000 |
| Moment arm, L<sub>m</sub> (ft) | 3.94 | 5.47 |
| Bending limit, M<sub>cap</sub>/L<sub>m</sub> (lb) | 15,244 | 10,962 |
| Force used, per pile (lb) | 15,244 | 10,962 |
| Force used, per foot of slope (lb/ft) | 2,540.7 | 1,827.0 |

The lower shaft's soil force is inside its shear capacity and still nearly three
times what its section can carry in bending; the upper shaft exceeds both — by
78% in shear and by a factor of seven and a half in bending. The caps cut the
delivered force by 65% and 87%, and the deeper pile — the one the soil pushes
hardest on — ends up delivering *less* than the shallower one, because its
pressure centroid sits further above the failure surface and the same moment
capacity buys a smaller force at a longer arm.

To see how much the caps matter, remove them. This too is a **single-surface
run on the search's own critical circle, with only the two capacity cells
changing**:

| Capacities given | ΣH (lb/ft) | Spencer FS |
|---|:---:|:---:|
| Vcap and Mcap (as shipped) | 4,367.6 | 1.842 |
| Mcap only | 4,367.6 | 1.842 |
| Vcap only | 15,029.7 | 4.207 |
| Neither | 20,984.6 | 12.739 |

Bending alone reproduces the shipped answer exactly, so the moment capacity is
doing all of the limiting on this surface. Left uncapped, the soil force would
report a factor of safety of 12.739 on a slope that stands at 1.149 without the
piles — an answer that asks the two shafts to carry three and seven and a half
times the moment their sections are designed for.

---

## What the spacing is worth

The capacities are fixed by the concrete section, but the spacing is a
choice — it is the main variable a designer adjusts once the diameter is
set. Closing the spacing might be expected to spread the load, each shaft
carrying less. That would be true if the shafts were sharing a fixed
demand — a required support divided among however many piles are there to
carry it. The Ito & Matsui force is not a demand. It is a supply: the
largest force the sliding soil can press onto one pile before it squeezes
through the gap and escapes downhill. Narrow the gap and escape gets
harder, so the soil loads each pile more, not less — the sweep below
measures one pile carrying 188,549 lb at 3 ft spacing, four times what the
same pile carries at 6 ft. The limiting cases say the same thing: spread
far apart, the shafts are isolated posts the slope flows around; pushed
together until they touch, they are a wall, and a wall resists more per
foot than any row of separated posts.

The analysis then works per foot of slope: one pile serves S feet of it,
so the force per foot is the per-pile force divided by the spacing. At 6 ft
that is 44,178 / 6 = 7,363 lb/ft; at 3 ft it is 188,549 / 3 = 62,850. One
knob moved — the spacing — and the per-foot force rose 8.5-fold: a factor
of 4.3 because the arching raised the per-pile force, and the remaining
factor of 2 because the divisor is the spacing itself.

Widening the spacing to 12 ft is one cell on each row — **S** from `6` to `12`,
with **H** still blank so the force is recomputed at the new spacing:

![The piles editor at 12 ft spacing](images/lem12_studio_piles_spacing.png)

Sweeping the spacing across the applicable band is a **held-surface study: one
circle, the search's own critical surface, and nothing changing but S**. The
diameter stays at 2 ft, so S/D runs from 1.5 to 8. The soil-force column is the
uncapped Ito & Matsui force, lower row then upper row; ΣH is what survives the
capacities and reaches the slices:

| S (ft) | S/D | A₁ (ft) | A₂ (ft) | Soil force per pile (lb) | ΣH used (lb/ft) | Spencer FS |
|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| 3 | 1.5 | 41.709 | 17.181 | 188,549 / 336,891 | 8,386.6 | 2.354 |
| 4 | 2.0 | 15.184 | 7.526 | 76,378 / 138,639 | 6,405.7 | 2.072 |
| 5 | 2.5 | 9.792 | 5.564 | 53,580 / 98,345 | 5,192.2 | 1.929 |
| 6 | 3.0 | 7.569 | 4.755 | 44,178 / 81,729 | 4,367.6 | 1.842 |
| 8 | 4.0 | 5.621 | 4.046 | 35,938 / 67,165 | 3,316.3 | 1.741 |
| 10 | 5.0 | 4.741 | 3.725 | 32,217 / 60,588 | 2,673.4 | 1.684 |
| 12 | 6.0 | 4.241 | 3.544 | 30,104 / 56,854 | 2,239.4 | 1.648 |
| 16 | 8.0 | 3.695 | 3.345 | 27,796 / 52,774 | 1,690.7 | 1.604 |

The arching coefficients fall steeply as the gap opens — A₁ loses nine tenths of
its value between S/D = 1.5 and 8 — and the force per foot of slope falls with
them, from 8,386.6 lb/ft to 1,690.7. The table touches both limits of
the applicable band. The top row, at S/D = 1.5, is below it: that close, the
row acts more like a continuous wall than a line of piles soil can flow
between. (A wall itself is modeled the other way — a pile line with its
resistance stated directly per foot of wall, from the wall's structural
capacity, the way the soldier-pile wall of
[LEM-9](lem09_tieback_wall.md#the-problem) enters its 5,900 lb/ft — because
a wall has no gaps for the soil to arch across.) The bottom row sits at the band's upper edge, S/D = 8; spaced any
wider, the arching the theory is built on fades away and the computed force
overestimates what the row can develop. XSLOPE says so before the run: enter
the 3 ft spacing and the **Model checks** panel warns that S/D = 1.5 is
below the S/D of about 2 to 8 the theory is applicable for, and the Log adds
a warning at solve time that the computed force is very large.

Widening the spacing does not only lower the number on this surface. Searching
the model at 12 ft — a separate search with its own critical surface — puts the
answer somewhere else entirely:

![Spencer at 12 ft spacing](images/lem12_solution_wide.png){width=1000}

**FS = 1.409**, and the surface is not the deep one the table's 1.648 was measured
on: it leaves the toe, reaches 10.1 ft below the ground at its deepest, exits
8.5 ft behind the crest, and carries 20,525 lb/ft against the deep surface's
62,198. At 6 ft spacing this shallower surface stood above 1.842 — it lost to the
deep circle, because the rows crossing it supplied enough force to keep it
up. At 12 ft they supply less, and it drops to 1.409 and becomes the
answer. Widening the spacing did not just lower the factor of safety on the
deep circle; it changed which surface governs. That is why the held sweep
is only for reading the trend, and any spacing actually being considered
gets its own search.
[VP106](../verification/rocscience.md#vp106) is this same sweep on Cai & Ugai's
pile-reinforced slope, at four spacings, with the computed forces checked against
Slide2 and against the paper.

---

## Giving the force yourself

Everything so far has used the computed path. The other way to enter a pile
is to state its force outright, as
[LEM-9's soldier pile](lem09_tieback_wall.md#the-problem) does — a number from
a p-y analysis (a lateral load–deflection model of the shaft), a structural
check or a published chart, entered per foot of slope.
Typing one into **H** turns the Ito & Matsui calculation off for that row; nothing
else changes.

A natural test is to enter the forces the automatic run computed — the
two values the report printed:

| H |
|:---:|
| 2540.7 |
| 1827.0 |

![The piles editor with the force stated on both rows](images/lem12_studio_piles_h.png)

**Solved on the search's own critical circle — one surface, no search, only the
H column changing — the two paths give the same answer**, 1.842 either way: a
stated force and a computed force enter the equilibrium equations identically.
The difference is not in the arithmetic but in what happens when the surface
moves.

Searching with the stated forces shows the difference:

![Spencer with the force stated instead of computed](images/lem12_solution_statedh.png){width=1000}

**FS = 1.752**, on a shallower circle from the toe — tangent at elevation −0.04,
14.5 ft below the ground at its deepest, carrying 32,929 lb/ft. Held on that
circle, the two ways of entering the pile no longer agree:

| Surface | Computed (H blank) | Stated (2,540.7 / 1,827.0) | No piles |
|---|:---:|:---:|:---:|
| The auto search's critical circle | 1.842 | 1.842 | 1.481 |
| The shallower circle above | 1.896 | 1.752 | 1.213 |

The difference between the two entries is simple. With Ito & Matsui, every
trial surface gets its own H, because the force depends on how deep the
surface crosses each row; a stated H never changes. The two searches can
therefore settle on different circles — and here they do. On this shallow
circle the computed path develops **5,178.6 lb/ft** (a shallower crossing
on the lower row, a shorter moment arm on the upper), while the stated
numbers stay frozen at the 4,367.6 that was right for the deep circle —
16% less than these piles would supply here. Under-credited, the surface
reads 1.752 instead of 1.896, and the stated-force search reports it as
the minimum. A stated force is only correct for the surface it was
computed from, and the surface that governs is not known until the search
finishes. [VP54](../verification/rocscience.md#vp54) is the verified case of the
stated-force route, a micro-pile row entered with its shear resistance given
directly.

Choosing between the two comes down to where the force should come from.
Leave `H` blank for the usual stabilizing case — a row of vertical shafts
with a known diameter and spacing inside the applicable S/D band — so every
trial surface is credited with the force it would actually receive. State
`H` when the force comes from a better source than the arching theory: a
p-y analysis, a load test, or a wall's own structural rating; when the
member is a wall or a battered pile (Ito & Matsui requires vertical piles
and a gap to arch across); or when reproducing a model that specifies the
force, as VP54 does. And whenever a stated force meets a search, check
that the surface the search returns matches the failure depth the stated
force was derived for — the comparison above is what happens when it does
not.

---

## Looking above the pile row

One question remains before the design can be trusted: did the search look
everywhere it should? Every search so far started from the circle on the
sheet, and that circle
reaches the rigid base. **Grid search (auto-seed the circular search)** ignores
the circles sheet and sweeps a grid of centers against a range of tangent
elevations instead ([LEM-10](lem10_global_minimum.md#grid-search) is where that
tool is built). Back in **Run LEM…**, tick it and leave everything else:

![The Run LEM dialog with grid search on](images/lem12_studio_run_lem_grid.png)

Run again:

![Spencer with grid seeding](images/lem12_solution_bypass.png){width=1000}

**FS = 1.702**, below the 1.842 the seeded search reported, on a surface that
never engages the piles. It daylights (cuts back up to the ground surface)
at (10.00, 10.00), the head of the upper shaft, runs back to the crest at
x = 24.3, and is 6.5 ft deep at its deepest,
moving 6,867 lb/ft. The depth from the ground surface to the failure surface at
the pile is zero, so the Ito & Matsui integral over that depth is zero: the whole
pile row contributes **0.04 lb/ft** to a surface that slides over the top of it.

This is an important check for any pile design. The piles hold the deep mechanism, and
holding it promotes whatever the next mechanism is; here the next one runs above
the pile heads and takes nothing from them. It is not a numerical artifact —
**Min slip depth** at 5 ft leaves it untouched, because 6.5 ft of clay is moving —
and the same slope searched without piles at all returns 1.149 from either
seeding, so it is adding the piles that makes the shallow surface critical.
The same behavior is
recorded on [VP54](../verification/rocscience.md#vp54), where a free search on a
micro-piled slope finds a circle exiting upslope of the pile row.

Raising the pile heads, lengthening the row upslope, or adding a row above the
crest are the design answers; each is a change to the model, and each has to be
searched again from a grid.

---

## Conclusion

This tutorial covered:

- A pile row whose **H** is blank, which computes the lateral force from the
  diameter and spacing by Ito & Matsui, per trial surface — with the chain
  printed one Log line per row.
- The piled search against the unpiled one, and each row's own share measured
  on a held circle.
- The structural capacities capping the soil force — the moment capacity
  governed both rows here.
- Spacing as supply: closer piles arch more and divide by less, and wider
  spacing can hand the slope back a different mechanism.
- A stated force against a computed one: identical on the surface it came
  from, stale one surface away.
- Checking for the shallow bypass — a surface that daylights at a pile head
  and takes almost nothing from the row it slides over.

**Where to go next:** the [tutorials index](index.md) lists the series.
[Piles and Concrete Piers](../lem/piles.md) derives the Ito & Matsui equations,
the capacity checks and the per-unit-width convention;
[Sample Problem 10](../lem/samples.md#10-slope-stabilized-with-piles) catalogues
this model with its factor of safety by every method;
[LEM-9](lem09_tieback_wall.md) enters a pile with its force stated, beside the
tieback anchors it carries; and [VP106](../verification/rocscience.md#vp106) and
[VP54](../verification/rocscience.md#vp54) check the computed and stated routes
against published solutions.
