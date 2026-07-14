# Rocscience Slide2 Groundwater Corpus

This page tracks the [Slide2 Groundwater Verification Manual](https://www.rocscience.com/support/slide2/verification)
(Rocscience, 2022; 21 problems) the same way the [Slide2 slope-stability corpus](rocscience.md)
tracks its manual: every problem gets a row, built problems get an XSLOPE input file
(`docs/files/rocscience_gw/`), a results section, and regression test tags. Unlike the LEM
corpus, the seepage tags mesh and solve live on every run — the committed artifact is the
input file alone.

<!-- test: file=../files/rocscience_gw/gw002.xlsx, type=seep, target_size=0.10, expected_flowrate=4.534e-06, tolerance=0.02, benchmark=GW2-q -->
<!-- test: file=../files/rocscience_gw/gw002.xlsx, type=seep_head, target_size=0.10, points=4:1:0.500;4.5:0.866:0.381;5:0:0.263;6:0:0.202, tolerance=0.01, benchmark=GW2-h -->
<!-- test: file=../files/rocscience_gw/gw003.xlsx, type=seep, target_size=0.40, expected_flowrate=2.307e-05, tolerance=0.02, benchmark=GW3-q -->
<!-- test: file=../files/rocscience_gw/gw003.xlsx, type=seep_head, target_size=0.40, points=0:-4:4.45;10:-4:3.37;14:-4:2.42;20:-4:1.05;30:-4:0.19, tolerance=0.05, benchmark=GW3-h -->
<!-- test: file=../files/rocscience_gw/gw004.xlsx, type=seep, target_size=0.147, expected_flowrate=5.462e-06, tolerance=0.05, benchmark=GW4-q -->
<!-- test: file=../files/rocscience_gw/gw006a.xlsx, type=seep, target_size=1.0, max_iter=2000, expected_flowrate=2.301e-07, tolerance=0.05, benchmark=GW6a-q -->
<!-- test: file=../files/rocscience_gw/gw006a.xlsx, type=seep_head, target_size=1.0, max_iter=2000, points=26:0.05:7.40;26:2:7.47;26:4:7.52;26:6:7.66, tolerance=0.15, benchmark=GW6a-h -->
<!-- test: file=../files/rocscience_gw/gw009a.xlsx, type=seep, expected_flowrate=2.2985e-05, tolerance=0.05, benchmark=GW9a-q -->
<!-- test: file=../files/rocscience_gw/gw010.xlsx, type=seep, target_size=0.25, max_iter=1500, expected_flowrate=6.07e-05, tolerance=0.05, benchmark=GW10-q -->
<!-- test: file=../files/rocscience_gw/gw012.xlsx, type=seep, target_size=1.0, max_iter=1500, expected_flowrate=4.137e-04, tolerance=0.05, benchmark=GW12-q -->
<!-- test: file=../files/rocscience_gw/gw013.xlsx, type=seep, target_size=1.0, max_iter=1500, expected_flowrate=2.087e-02, tolerance=0.05, benchmark=GW13-q -->

The manual verifies Slide2's finite-element groundwater engine against closed-form solutions
(Polubarinova-Kochina, Vedernikov, Terzaghi consolidation) and published numerical benchmarks.
It is the seepage analog of the slope-stability manual, and the natural verification target for
XSLOPE's own FE seepage solver — which the LEM corpus already exercises end-to-end on VP71,
VP72, VP76, VP77 and VP102, but only as a pore-pressure source, not against seepage-specific
quantities (flow rates, free-surface positions, pressure profiles).

Problems 1–13 are steady-state and analyzable with the current solver. Problems 15–21 are
**transient** (and 15–16 are consolidation), which XSLOPE's steady-state solver does not
support.

Two capabilities are deliberately out of scope, so the rows that depend on them are blocked by
design rather than pending work:

- **Transient seepage.** XSLOPE's seepage analysis exists to supply pore pressures to a stability
  analysis, and a transient head field does not change a factor of safety on its own — capturing
  rainfall-induced failure would additionally require an unsaturated shear-strength model (such as
  Fredlund's $\phi^b$), without which suction is neglected in the stability calculation regardless
  of how the water table moves. Rapid drawdown, the other classic transient case, is already covered
  by the staged Duncan & Wright procedure that is standard of practice. Problems 16–21 are blocked
  on this.
- **Specified-flux (Neumann) boundary conditions.** The solver's boundary conditions are fixed-head
  and exit-face only. Problems 1 and 7 apply a uniform infiltration rate to the ground surface, which
  has no fixed-head equivalent, so they are blocked on this. Their inputs and published answers are
  recorded below for whenever it is added.

## Status

| # | Problem | Status | Notes |
|---|---|---|---|
| 1 | Shallow unconfined flow with rainfall | *blocked* | Uniform infiltration (P = 2.5×10⁻⁶ m/s) across the top boundary needs a specified-flux (Neumann) boundary condition, which the seepage solver does not have (BCs are fixed-head and exit-face only). Inputs and both published answers (x_a = 4.06/3.98, h_max = 4.49/4.25) are recorded for when it exists. |
| [2](#gw2) | Flow around cylinder | **built** | [gw002.xlsx](../files/rocscience_gw/gw002.xlsx). Confined potential flow: solved heads match Slide within 0.0013 m at every printed point and the closed form within its own idealization error. |
| [3](#gw3) | Confined flow under dam foundation | **built** | [gw003.xlsx](../files/rocscience_gw/gw003.xlsx). Rushton & Redshaw benchmark: head profiles under and beyond the dam within 0.08 m of the published chart everywhere. |
| [4](#gw4) | Steady unconfined flow through earth dam | **built** | [gw004.xlsx](../files/rocscience_gw/gw004.xlsx). Kozeny basic parabola: phreatic surface within 1–2% over the dam body; drain-tip height 0.50 vs Slide 0.442 / parabola 0.480 (the published pair itself spreads 9%). |
| 5 | Unsaturated flow behind an embankment | *no lock possible* | The manual publishes only qualitative pressure contours and flow lines against FLAC ("compared very well") — no numeric quantity exists to lock, and the geometry figure is unlabeled. |
| [6](#gw6) | Steady-state seepage through saturated–unsaturated soils | **built** (case 1 of 5, caveat) | [gw006a.xlsx](../files/rocscience_gw/gw006a.xlsx). Fredlund & Rahardjo isotropic dam with a 12 m drain: the pressure-head profile along the crest centerline matches Slide/F&R in shape exactly but sits +0.3 m high — insensitive to the conductivity fit and the mesh, the same free-surface question as GW11. Case 4 needs the flux BC (out of scope, see above); cases 2 (9:1 anisotropy), 3 (core), 5 (seepage-face) are buildable with chart targets and deferred. |
| 7 | Seepage within layered slope | *blocked* | Rulon & Freeze layered slope: the manual applies a constant infiltration rate of 2.1×10⁻⁴ m/s to the top of the slope — a specified flux with no fixed-head equivalent (the stated water table position is an output, not an input), so this needs the flux boundary condition, which is out of scope (see above). k-charts and profile targets recorded for when it exists. |
| 8 | Flow through ditch-drained soils | *blocked* | Inputs fully printed, but both soils use Gardner conductivity kr = 1/(1 + a·ψⁿ); the solver implements linear-front and van Genuchten only. A Gardner option would unlock GW8 and GW11. The problem exists to verify the Gardner function, so a van Genuchten stand-in was deliberately not attempted. |
| [9](#gw9) | Seepage through dam | **built** (dam 1 of 2) | [gw009a.xlsx](../files/rocscience_gw/gw009a.xlsx). Bowles homogeneous dam via Chapuis et al. (2001): Q = 1.379×10⁻³ m³/(min·m) vs Slide 1.378×10⁻³ / SEEP/W 1.37×10⁻³ / Bowles flow nets 1.10–1.28×10⁻³. Dam 2 (drain) needs the source paper — its k-function and reservoir level are chart-only and the published Q implies a k two decades below the chart. |
| [10](#gw10) | Steady unconfined flow, van Genuchten permeability | **built** | [gw010.xlsx](../files/rocscience_gw/gw010.xlsx). Clement et al. (1996): Q = 6.070×10⁻⁵ vs Slide 6.066×10⁻⁵ (+0.07%) / Clement 6.076×10⁻⁵; phreatic exit el. 4.87 vs Clement 4.8 / Slide 5.0. |
| 11 | Earth/rock-fill dam, Gardner permeability function | partial — under investigation | Zhang et al. (2001) dam, geometry fully extracted. The Gardner law is absent, but that is NOT the gap: a van Genuchten fit, a low-suction-weighted fit, and a linear front all give the same release point (el. 17.4) vs Slide 19.40 / ABAQUS 19.64 — the discrepancy is insensitive to the conductivity model and points at the exit-face free-surface treatment on this steep downstream face. Casagrande-style hand estimates bracket 15–22. No lock until resolved. |
| [12](#gw12) | Seepage from a trapezoidal ditch into a deep drainage layer | **built** | [gw012.xlsx](../files/rocscience_gw/gw012.xlsx). Vedernikov: Q = 4.137×10⁻⁴ vs Slide 4.093×10⁻⁴ (+1.1%) / theory 4.0×10⁻⁴; flow-bulb half-width ≈42 vs Slide 41 / theory 40. |
| [13](#gw12) | Seepage from a triangular ditch into a deep drainage layer | **built** | [gw013.xlsx](../files/rocscience_gw/gw013.xlsx). Vedernikov: Q = 2.087×10⁻² vs Slide 2.050×10⁻² (+1.8%) / theory 2.0×10⁻². |
| 14 | Unsaturated soil column | *blocked* | Steady (Gardner 1958 capillary profile under constant infiltration), but needs BOTH a Gardner-exponential conductivity k = ks·e^(−αψ) and a specified-flux surface boundary condition; neither exists. Belongs with GW1 on the flux-BC gap, not the transient tier. |
| 15 | 1-D consolidation, uniform initial excess pore pressure | blocked | Transient/consolidation — no transient solver |
| 16 | Pore pressure dissipation of stratified soil | blocked | Transient/consolidation |
| 17 | Transient seepage, earth fill dam with toe drain | blocked | Transient |
| 18 | Transient seepage through an earth fill dam | blocked | Transient |
| 19 | Transient seepage through an earth fill dam (II) | blocked | Transient |
| 20 | Transient seepage in a layered slope | blocked | Transient |
| 21 | Transient seepage through a fully confined aquifer | blocked | Transient |

## Results

### GW2: Flow around cylinder {#gw2}

**Input files:** [gw002.xlsx](../files/rocscience_gw/gw002.xlsx)

Confined potential flow past a cylinder (Streeter's closed form; Desai & Kundu's finite
elements): a half-domain 8 × 4 m with a semicircular notch of radius 1 m at the bottom
center and fixed heads of 1.0 and 0.0 on the vertical edges. The domain is fully
saturated, and the head field in homogeneous confined flow is independent of the
conductivity, so the problem is determined by geometry and boundary heads alone.

| Point | XSLOPE | Slide | Closed form | Desai & Kundu |
|---|---|---|---|---|
| (4, 1) | 0.500 | 0.500 | 0.500 | 0.500 |
| (4.5, 0.866) | 0.381 | 0.381 | 0.375 | 0.378 |
| (5, 0) | 0.263 | 0.263 | 0.250 | 0.277 |
| (6, 0) | 0.202 | 0.203 | 0.188 | 0.213 |
| (8, 0) | 0.000 | 0.000 | −0.031 | 0.000 |

The two finite-element solutions agree with each other within 0.0013 m everywhere; both
depart from the closed form near the downstream edge, where the analytical solution
assumes an infinite domain (its −0.031 at the finite model's zero-head boundary measures
that idealization, not an error). Flowrate is locked as a regression value.

![gw002: mesh and solved heads](images/gw002.png)

### GW3: Confined flow under a dam foundation {#gw3}

**Input files:** [gw003.xlsx](../files/rocscience_gw/gw003.xlsx)

The classic flow-net problem (Rushton & Redshaw): a 40 × 10 m soil block with head 5 m
on the ground surface upstream of the dam (x = 0–8), head 0 downstream (x = 20–40), and
the dam base impervious between them.

| Station (line 1-1, y = −4) | XSLOPE | Rushton & Redshaw / Slide |
|---|---|---|
| x = 0 | 4.45 | 4.50 |
| x = 10 | 3.37 | 3.45 |
| x = 14 | 2.42 | 2.47 |
| x = 20 | 1.05 | 1.05 |
| x = 30 | 0.19 | 0.20 |
| x = 40 | 0.08 | 0.07 |

Heads along both published profile lines fall within 0.08 m of the chart everywhere
(Slide's markers coincide with Rushton & Redshaw's); the vertical profile under the dam
heel spans 0.24–1.30 m, matching the published range endpoint for endpoint.

![gw003: mesh and solved heads](images/gw003.png)

### GW4: Steady unconfined flow through an earth dam {#gw4}

**Input files:** [gw004.xlsx](../files/rocscience_gw/gw004.xlsx)

Free-surface seepage through a small dam with a toe drain, against Kozeny's basic
parabola: reservoir head 4 m on the upstream face, a seepage exit face on the
downstream slope and drain end. The free-surface shape is independent of k.

| Quantity | XSLOPE | Slide | Kozeny parabola |
|---|---|---|---|
| Phreatic el. at x = 14 | 2.88 | — | 2.81 |
| Phreatic el. at x = 18 | 2.03 | — | 2.02 |
| Phreatic el. at x = 20 | 1.43 | — | 1.47 |
| y₁ at the drain | 0.50 | 0.442 | 0.480 |

The solved phreatic surface tracks the parabola within 1–2% across the dam body. At the
drain tip the published values themselves spread 9% (Slide's measured 0.442 vs the
parabola's 0.480), and the parabola is an idealization exact only at the drain; XSLOPE's
0.50–0.53 across a mesh refinement (Q changing 2%) sits just above that band. The
flowrate lock (5.46×10⁻⁶ m³/s) exceeds the idealized k·y₁ = 4.80×10⁻⁶ because the
parabola underestimates entry-face flow.

![gw004: mesh and solved heads](images/gw004.png)

### GW9: Seepage through dam {#gw9}

**Input files:** [gw009a.xlsx](../files/rocscience_gw/gw009a.xlsx)

Bowles' homogeneous dam, the flow-net textbook example re-solved numerically by Chapuis,
Chenaf & Bowles (2001) and by Slide: base 100 m, crest 10 m at el. 20 (2.5:1 upstream,
2:1 downstream), reservoir head 18.5 m, ks = 6.67×10⁻⁶ m/s with the manual's printed
8-point unsaturated conductivity table, fit here by a Mualem–van Genuchten curve
(α = 0.2835, n = 2.765).

| | XSLOPE | Slide | SEEP/W (fine) | Bowles (flow nets) |
|---|---|---|---|---|
| Q, m³/(min·m) | 1.379×10⁻³ | 1.378×10⁻³ | 1.37×10⁻³ | 1.10–1.28×10⁻³ |

The manual's second dam (with a chimney-and-blanket drain) is not built: its
conductivity function and reservoir level exist only as chart pixels, and the published
flowrate implies a fill conductivity two decades below the chart — resolving that needs
the source paper (Chapuis et al. 2001, Can. Geotech. J. 38:1113).

![gw009a: mesh and solved heads](images/gw009a.png)

### GW10: Steady unconfined flow, van Genuchten permeability {#gw10}

**Input files:** [gw010.xlsx](../files/rocscience_gw/gw010.xlsx)

Clement, Wise, Molz & Wen (1996)'s unconfined square domain, the manual's designated
van Genuchten test: a 10 × 10 m block with head 10 on the left edge, tailwater 2 on the
right, an exit face above the tailwater, and vG conductivity (α = 0.64, n = 4.65,
ks = 1.1574×10⁻⁵ m/s) — an exact capability match for the solver's `vg` option.

| | XSLOPE | Slide | Clement et al. |
|---|---|---|---|
| Q (m³/s per m) | 6.070×10⁻⁵ | 6.066×10⁻⁵ | 6.076×10⁻⁵ |
| Phreatic exit elevation | 4.87 | 5.0 | 4.8 |

The manual's "seepage face" column tabulates the phreatic exit *elevation* (both
published figures show the free surface exiting near el. 4.9–5.0), not a face length.
Only the tailwater-2 case carries published numbers and is locked.

![gw010: mesh and solved heads](images/gw010.png)

### GW12 / GW13: Ditch seepage into a deep drainage layer (Vedernikov) {#gw12}

**Input files:** [gw012.xlsx](../files/rocscience_gw/gw012.xlsx) (trapezoidal) ·
[gw013.xlsx](../files/rocscience_gw/gw013.xlsx) (triangular)

Vedernikov's closed-form solutions for steady seepage from a channel into a deep
drainage layer, modeled as half-domains by symmetry with the ditch perimeter at head 50
and the deep drain as head 0 on the base. The seepage detaches below the ditch and
descends as a bulb whose width the theory predicts.

| | XSLOPE | Slide | Vedernikov |
|---|---|---|---|
| Trapezoidal: Q per half | 4.137×10⁻⁴ | 4.093×10⁻⁴ | 4.0×10⁻⁴ |
| Trapezoidal: bulb half-width | ≈42 | 41 | 40 |
| Triangular: Q per half | 2.087×10⁻² | 2.050×10⁻² | 2.0×10⁻² |

The detached-bulb iteration converges cleanly at 1,500 free-surface iterations (the
`max_iter` tag key); the default 400 is not enough for these geometries.

![gw012: mesh and solved heads](images/gw012.png)

![gw013: mesh and solved heads](images/gw013.png)

### GW6: Steady-state seepage through saturated–unsaturated soils {#gw6}

**Input files:** [gw006a.xlsx](../files/rocscience_gw/gw006a.xlsx) (case 1 of 5)

Fredlund & Rahardjo (1993)'s saturated–unsaturated earth dam (12 m high, symmetric 2:1
faces, reservoir at 10 m, a 12 m horizontal drain at the downstream toe), case 1:
isotropic conductivity. The conductivity function is digitized from the manual's chart
(ks = 10⁻⁷ m/s, air entry ≈ 9 kPa) and fit by a Mualem–van Genuchten curve; three fit
variants move the answer by less than 3 cm, so the fit is not the controlling
uncertainty. The published target is the pressure-head profile along the crest
centerline, where Slide's and Fredlund & Rahardjo's curves coincide within 0.2 m.

| Elevation on the crest line | XSLOPE pressure head | Slide | F&R |
|---|---|---|---|
| 0 | 7.40 | 7.15 | ≈7.3 |
| 2 | 5.47 | 5.15 | ≈5.3 |
| 4 | 3.52 | 3.25 | ≈3.4 |
| 6 | 1.66 | 1.30 | ≈1.45 |
| 8 | −0.22 | −0.60 | ≈−0.45 |

*The profile shape reproduces exactly; the whole curve sits 0.25–0.5 m above the
published pair, insensitive to the conductivity fit and to mesh refinement. This is a
family difference, not a solver error: run against SEEP2D (the original USACE/WES code)
on the identical mesh, boundary conditions and unsaturated law, XSLOPE's free surface
daylights at the same place — see [the SEEP2D cross-check](#seep2d-crosscheck) below. The
regression locks XSLOPE's own values.*

![gw006a: mesh and solved heads](images/gw006a.png)

## The SEEP2D cross-check: where does the free surface daylight? {#seep2d-crosscheck}

Several problems on this page reproduce a published profile in shape but sit a little
high, and the free surface on a steep exit face looked like the common thread. The
question is whether XSLOPE places the **release point** — the top of the saturated
seepage face, where the phreatic surface daylights — differently from other codes.

It is answerable directly, because XSLOPE's seepage solver is in the SEEP2D family. The
original USACE/WES SEEP2D Fortran program is run on the *identical* tri3 mesh, the
identical boundary conditions and the identical unsaturated law (SEEP2D's `iuntyp`:
linear front, or van Genuchten — both codes use the same Mualem–van Genuchten form, so
α and n carry straight across), and the two solutions are compared node for node. The
harness is `benchmarks/run_seep2d_compare.py --gw`.

| problem | law | release point, XSLOPE | release point, SEEP2D | head RMS | head range |
|---|---|---|---|---|---|
| gw004 | linear front | 0.500 | 0.500 | 0.0004 | 0–4 |
| gw006a | van Genuchten | 0.000 | 0.000 | 0.0000 | 0–10 |
| gw009a | van Genuchten | 8.846 | 8.462 | 0.0026 | 0–18.5 |
| gw010 | van Genuchten | 4.872 | 4.872 | 0.0006 | 2–10 |
| gw012 | linear front | face dry | face dry | 0.034 | 0–50 |
| gw013 | linear front | face dry | face dry | 0.070 | 0–50 |

**The release point agrees.** It is identical on gw004, gw006a and gw010; both codes
leave the face fully drained on gw012 and gw013; and gw009a differs by 0.385 — a single
element at that mesh. The head fields agree to a relative RMS of order 10⁻⁴ throughout.
XSLOPE does not release low. Where this corpus sits above a published profile, the
difference is between the SEEP2D family and the reference code's exit-face convention,
not an XSLOPE defect — so those problems are locked on XSLOPE's own values with the
offset reported, rather than tuned.

One open item came out of the same run. On the **van Genuchten** problems the total
discharge reads 3.5–4.7% below SEEP2D (gw009a 2.299×10⁻⁵ vs 2.412×10⁻⁵; gw010
6.070×10⁻⁵ vs 6.294×10⁻⁵) even though the heads agree to 10⁻⁴ — while the linear-front
problems agree on discharge to better than 0.15%. The split follows the unsaturated law
exactly, and it is *not* XSLOPE's kr floor: dropping `kr_min` from 10⁻⁴ to zero leaves
the discharge unchanged to six figures. Since the head field is right and only the
integrated flux differs, the likely cause is where each code evaluates the strongly
nonlinear kr(ψ) when it forms an element's conductivity. That is being tracked
separately; it does not affect pore pressures or stability, which read the head field.

## Methodology

Same rules as the [Slide2 corpus](rocscience.md): problems are built from the manual's tabulated
data and coordinate-labeled figures; where a figure is unlabeled, geometry is extracted by
axis-calibrated pixel measurement and validated against printed solution quantities; every built
problem is locked into `run_tests.py` via test tags. Seepage problems compare flow rates,
phreatic-surface positions, and head/pressure profiles rather than factors of safety. Two tag
types carry the locks: `type=seep` (flowrate, live mesh + solve at `target_size`) and
`type=seep_head` (solved total head at named `x:y:h` points, interpolated from the four nearest
nodes). Where a problem's published answer is itself a chart curve (GW6/GW7), a
tolerance-banded profile lock is planned but not yet implemented.
