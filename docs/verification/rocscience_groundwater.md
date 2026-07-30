# Rocscience Slide2 Groundwater Corpus

The [Slide2 Groundwater Verification Manual](https://www.rocscience.com/help/slide2/verification-theory/verification-manuals)
(Rocscience, 2022) verifies Slide2's finite-element groundwater engine against closed-form
solutions (Polubarinova-Kochina, Vedernikov, Terzaghi consolidation) and published numerical
benchmarks. This page compares XSLOPE's finite-element seepage solver against the same
21 problems. Each problem has a row in the summary table; each built problem has an XSLOPE
input file, a results section, and figures. The quantities compared are seepage-specific —
flow rates, free-surface positions, head and pressure profiles — rather than factors of
safety.

Full bibliographic details for the author-year citations on this page are on the
shared [References](references.md) page.

<!-- test: file=files/rocscience_gw/gw001.xlsx, type=seep, target_size=0.2, expected_flowrate=2.500e-05, tolerance=0.02, benchmark=GW1-q -->
<!-- test: file=files/rocscience_gw/gw001.xlsx, type=seep_head, target_size=0.2, points=2:2:4.052;4:2:4.150;6:2:4.030, tolerance=0.02, benchmark=GW1-h -->
<!-- test: file=files/rocscience_gw/gw002.xlsx, type=seep, target_size=0.10, expected_flowrate=4.534e-06, tolerance=0.02, benchmark=GW2-q -->
<!-- test: file=files/rocscience_gw/gw002.xlsx, type=seep_head, target_size=0.10, points=4:1:0.500;4.5:0.866:0.381;5:0:0.263;6:0:0.202, tolerance=0.01, benchmark=GW2-h -->
<!-- test: file=files/rocscience_gw/gw003.xlsx, type=seep, target_size=0.10, expected_flowrate=2.351e-05, tolerance=0.02, benchmark=GW3-q -->
<!-- test: file=files/rocscience_gw/gw003.xlsx, type=seep_head, target_size=0.10, points=0:-4:4.47;10:-4:3.40;14:-4:2.44;20:-4:1.05;30:-4:0.19, tolerance=0.05, benchmark=GW3-h -->
<!-- test: file=files/rocscience_gw/gw004.xlsx, type=seep, target_size=0.06, max_iter=2500, expected_flowrate=5.487e-08, tolerance=0.02, benchmark=GW4-q -->
<!-- test: file=files/rocscience_gw/gw004.xlsx, type=seep_head, target_size=0.06, max_iter=2500, points=6:1:3.906;12:1:3.183;16:1:2.441;20:0.5:1.356, tolerance=0.01, benchmark=GW4-h -->
<!-- test: file=files/rocscience_gw/gw005.xlsx, type=seep, target_size=0.5, expected_flowrate=8.167e-11, tolerance=0.05, benchmark=GW5-q -->
<!-- test: file=files/rocscience_gw/gw005.xlsx, type=seep_head, target_size=0.5, points=5:2:9.041;15:2:7.090;25:2:5.093;15:8:9.855;35:1:4.087, tolerance=0.05, benchmark=GW5-h -->
<!-- test: file=files/rocscience_gw/gw006a.xlsx, type=seep, target_size=1.0, max_iter=2000, expected_flowrate=2.437e-07, tolerance=0.05, benchmark=GW6a-q -->
<!-- test: file=files/rocscience_gw/gw006a.xlsx, type=seep_head, target_size=1.0, max_iter=2000, points=26:0.05:7.40;26:2:7.47;26:4:7.52;26:6:7.66, tolerance=0.15, benchmark=GW6a-h -->
<!-- test: file=files/rocscience_gw/gw006b.xlsx, type=seep, target_size=1.0, max_iter=2000, expected_flowrate=1.639e-06, tolerance=0.05, benchmark=GW6b-q -->
<!-- test: file=files/rocscience_gw/gw006b.xlsx, type=seep_head, target_size=1.0, max_iter=2000, points=26:0.05:6.573;26:2:6.738;26:4:7.191;26:6:7.789, tolerance=0.05, benchmark=GW6b-h -->
<!-- test: file=files/rocscience_gw/gw006c.xlsx, type=seep, target_size=1.0, max_iter=2000, expected_flowrate=2.502e-08, tolerance=0.05, benchmark=GW6c-q -->
<!-- test: file=files/rocscience_gw/gw006c.xlsx, type=seep_head, target_size=1.0, max_iter=2000, points=26:0.05:5.602;26:2:6.056;26:4:6.678;26:6:7.478, tolerance=0.05, benchmark=GW6c-h -->
<!-- test: file=files/rocscience_gw/gw006e.xlsx, type=seep, target_size=1.0, max_iter=2000, expected_flowrate=1.686e-07, tolerance=0.05, benchmark=GW6e-q -->
<!-- test: file=files/rocscience_gw/gw006e.xlsx, type=seep_head, target_size=1.0, max_iter=2000, points=26:0.05:8.337;26:2:8.348;26:4:8.386;26:6:8.446, tolerance=0.05, benchmark=GW6e-h -->
<!-- test: file=files/rocscience_gw/gw009a.xlsx, type=seep, expected_flowrate=2.2985e-05, tolerance=0.05, benchmark=GW9a-q -->
<!-- test: file=files/rocscience_gw/gw009b.xlsx, type=seep, expected_flowrate=4.2885e-06, tolerance=0.05, benchmark=GW9b-q -->
<!-- test: file=files/rocscience_gw/gw010.xlsx, type=seep, target_size=0.25, max_iter=1500, expected_flowrate=6.07e-05, tolerance=0.05, benchmark=GW10-q -->
<!-- test: file=files/rocscience_gw/gw012.xlsx, type=seep, target_size=1.0, max_iter=1500, expected_flowrate=4.137e-04, tolerance=0.05, benchmark=GW12-q -->
<!-- test: file=files/rocscience_gw/gw013.xlsx, type=seep, target_size=1.0, max_iter=1500, expected_flowrate=2.087e-02, tolerance=0.05, benchmark=GW13-q -->

Problems 1–13 are steady-state. Problems 15–21 are **transient** (15–16 are
consolidation), solved with XSLOPE's
[transient seepage solver](../seep/transient.md) — the first three against
closed-form (or recomputed series) solutions:

- **15** — Terzaghi 1-D consolidation (single and double drainage), against the Eq 17.3 series.
- **16** — Pyrah two-layer consolidation, against a recomputed two-layer eigenfunction series.
- **21** — transient flow in a fully confined aquifer, against Ferris' erfc solution.
- **18** — transient flow through the GW6 earth dam (reservoir raised 4 m → 10 m, no drain),
  against the digitized Fig 20.5 toe-slope profile.
- **17** — the same dam with a toe drain, against the Fig 19-5 total-head contours (qualitative).
- **19** — transient flow below a lined lagoon as the pond is filled, against the digitized
  Fig 21.9 top-boundary pressure profile.
- **20** — transient flow through the GW7 layered slope as rainfall switches on, against the
  digitized Fig 22.7 query-line profile.

Problems **17–20** add an unsaturated-conductivity curve — 17/18 also a reservoir boundary
applied only where the face is submerged, 19/20 a stepped pond head and a stepped rainfall
flux. Their published targets are figures rather than closed forms; 18, 19 and 20 publish
line or marker plots on labeled axes and are compared numerically against digitizations of
them, while 17 publishes contour plates and stays qualitative.

Slide2 defines the conductivity and water-content curves independently, and each of these
four models ships both tables. The relative permeability is always fitted to the vendor's own
conductivity table, never derived from its retention curve. Whether one van Genuchten pair can
also carry the retention curve depends on the two ranges: [GW18](#gw18)'s retention line runs
to 4000 kPa, a hundred times its conductivity scale, so it is carried separately through the
Gardner option's independent drainage band; on 17, 19 and 20 the two ranges are the same
order and one pair carries both, which shifts transient *timing* — up to a measured factor of
two on 20 — while leaving the steady field unchanged.

The transient problems compare the head and pressure field as it evolves, not a factor of safety: a
transient head field does not change an FS on its own (rainfall-triggered failure would
additionally need an unsaturated shear-strength model such as Fredlund's $\phi^b$), and rapid
drawdown — the classic transient stability case — is handled by the staged Duncan & Wright
procedure. These problems exercise the seepage engine's storage and time-stepping directly.

Problems 1, 7, 8 and 14 apply a uniform infiltration rate to the ground surface, which has no
fixed-head equivalent. These require a
[specified-flux (Neumann) boundary condition](../seep/overview.md#specified-flux-boundary-conditions-neumann),
which the solver supports. Three of them — **1** (Dupuit recharge mound), **7** (Rulon &
Freeze layered slope) and **8** (ditch-drained aquifer) — are built on that boundary. The fourth,
**14**, is blocked for a different reason: its published quantity is the Gardner
(1958) capillary profile, which is derived for the *exponential* conductivity law
$k = k_s\,e^{\alpha\psi}$ — a law XSLOPE does not implement (its `gard` option is the power form
$k_r = 1/(1 + a\,\psi^n)$ that SEEP/W and Slide carry, a different function). With no exponential
law and no tabulated Slide value to compare against (the manual prints only Fig. 14.3/14.4 charts),
no comparable published value remains beyond the one-dimensional through-flux, which is
verified to machine precision.

## Status

**Match to the published value**

| Symbol | Meaning |
|---|---|
| 🟢 | within 3% of the vendor and/or reference figure |
| 🟡 | 3–6% |
| 🔴 | more than 6% |
| 🟣 | in progress |
| <span class="nodata">⊘</span> | insufficient data or out of scope |

The dot scores the **match quality of what is locked**, not how much of a problem is built; the partial/blocked detail is in the row text.

<div class="corpus-summary match" markdown>

Status values follow the [shared vocabulary](rocscience.md)
used across this section (**built**, *covered*, *partial*, *planned*, *blocked*,
*no lock possible*, *not supported*).

| # | Match | Problem | Results | Notes |
|---:|:-:|---|---|---|
| [1](#gw1) | 🟢 | Shallow unconfined flow with rainfall | Crest x_a ≈ 4.1 vs Slide 4.06 (~0.05 m) · h_max 4.61 vs Slide 4.49 (~0.12 m) · Q = P·L = 2.5×10⁻⁵ locked | slightly-high free-surface family, above Haar's Dupuit closed form (SEEP2D cross-check) |
| [2](#gw2) | 🟢 | Flow around cylinder | Solved heads within 0.0013 m of Slide at every printed point · closed form matched within its own idealization error | |
| [3](#gw3) | 🟢 | Confined flow under dam foundation | Head profiles under and beyond the dam within 0.08 m of the published Rushton & Redshaw / Slide chart everywhere | |
| [4](#gw4) | 🔴 | Steady unconfined flow through earth dam | Phreatic surface within 0.02–0.06 m of the Kozeny basic parabola over the dam body · y₁ above the drain toe 0.401 vs RS2's own solve of this model 0.395 (+1.5%) · vs Slide 0.442 (−9.3%) · vs Eq 4.1 0.486 (−17.5%) · drain-face entry offset x₁ 0.272 vs 0.226 (**+20.4%**, and it sets the dot) | RS2's Table 4.1 for the file this model is built from is the governing pairing and publishes both quantities; y₁ is green on it, x₁ is not, and x₁ turns on an unsaturated curve the vendor file does not store |
| [5](#gw5) | 🟢 | Unsaturated flow behind an embankment | Q = 8.167×10⁻¹¹ vs the one-dimensional closed form *k·b·i* = 8.0×10⁻¹¹ (+2.1%; +1.8% on the finest mesh) · solved pressure head inside Fig 5-4's own 1 m colour bands at 46 of 49 grid points (worst miss 0.025 m) | **built**; chart-keyed target, locked on XSLOPE's own field |
| [6](#gw6) | <span class="nodata">⊘</span> | Steady-state seepage through saturated–unsaturated soils | Pressure head along line 1-1: cases 2 and 5 reproduce the Slide/F&R curve almost exactly · cases 1 and 3 sit ~0.3–0.5 m high (mesh- and fit-insensitive; the published Slide/Ref[1] themselves scatter ~1.5 m near the crest on case 3) | **built** (4 of 5 cases); chart-only target, locked on XSLOPE's own field |
| [7](#gw7) | 🟢 | Seepage within layered slope | Total head along the manual's own query line within 0.005 m rms / 0.013 m worst of the Fig 22.7 steady markers over 21 stations (≈1% of the profile's head range) · water table at the toe el 0.30 vs the stated Slide / Rulon & Freeze 0.3 m (0.00 m) · perched zone and slope-face spring reproduced · Q = q·L = 1.68×10⁻⁴ locked | **built**; problem 7's own figures are chart-only, so the numeric target comes from problem 22's Fig 22.7 steady frame |
| [8](#gw8) | 🟢 | Flow through ditch-drained soils | Flux boundary exact — total inflow = *q*·*L*, the confined response matching the closed form to six figures · water table within 0.004–0.006 m of the Fig 8.3/8.4 line over the whole span (worst 2.4% of the 0.25 m divide mound) · Fig 8.3's labeled pressure-head contours within 0.010 m rms over 14 stations | **built**; flux rate and Soil B's Gardner *a* taken from the vendor model where it disagrees with the printed tables |
| [9](#gw9) | 🟢 | Seepage through dam | Dam 1: Q = 1.379×10⁻³ vs Slide 1.378×10⁻³ m³/(min·m) (+0.1%) · dam 2: Q = 4.29×10⁻⁶ vs Slide 4.23×10⁻⁶ m³/(s·m) (+1.4%) | **built** (both dams); body k read from Bowles (1984) Fig E9-2b, not the Chapuis caption |
| [10](#gw10) | 🟢 | Steady unconfined flow, van Genuchten permeability | Q = 6.070×10⁻⁵ vs Slide 6.066×10⁻⁵ (+0.1%) · vs Clement 6.076×10⁻⁵ (−0.1%) · phreatic exit el. 4.87 vs Slide 5.0 (−0.13 m) · vs Clement 4.8 (+0.07 m) | |
| [11](#gw11) | 🔴 | Earth/rock-fill dam, Gardner permeability function | Free-surface release point el. 17.90 (17.8 ± 0.15 across the mesh sweep) vs Slide 19.40 (−1.50 m) · vs ABAQUS 19.64 (−1.74 m) · the free surface itself within 0.76 m rms / 1.23 m worst of Slide's own drawn line over 18 stations | **built** (case 1 of 2, discrepancy); unsaturated law, mesh size, element order, the *k*<sub>r</sub> floor and the exit-face extent are each measured and none moves the release point |
| [12](#gw12) | 🟢 | Seepage from a trapezoidal ditch into a deep drainage layer | Q = 4.137×10⁻⁴ vs Slide 4.093×10⁻⁴ (+1.1%) · vs Vedernikov theory 4.0×10⁻⁴ (+3.4%) · flow-bulb half-width ≈42 vs Slide 41 (+1) · vs theory 40 (+2) | |
| [13](#gw12) | 🟢 | Seepage from a triangular ditch into a deep drainage layer | Q = 2.087×10⁻² vs Slide 2.050×10⁻² (+1.8%) · vs Vedernikov theory 2.0×10⁻² (+4.4%) | |
| [14](#gw14) | <span class="nodata">⊘</span> | Unsaturated soil column | | *blocked* — the closed form assumes an exponential conductivity law XSLOPE does not implement |
| [15](#gw15) | 🟢 | 1-D consolidation, uniform initial excess pore pressure | Isochrones within ≈0.3% of u₀ of the Terzaghi Eq 17.3 closed form | **built** (both cases) |
| [16](#gw16) | 🟢 | Pore pressure dissipation of stratified soil | Within ≈0.3–0.5% of u₀ of the recomputed Pyrah 1996 two-layer eigen-series | **built** (3 cases) |
| [17](#gw17) | <span class="nodata">⊘</span> | Transient seepage, earth fill dam with toe drain | At 15 h XSLOPE's *h* = 7 front stands 1.1–1.8 m inside the upstream face against RS2's 1.5–3.3 m · the steady field reproduces the Fig 19-5 total-head contours — reservoir 10 → drain 0, phreatic drawn to the toe — and at 16383 h, where RS2 is essentially steady, XSLOPE's station heads read 6.54 / 6.92 / 4.37 / 3.22 against a steady 7.20 / 7.52 / 5.84 / 4.40, settling to within 0.01 m of it by ≈5×10⁴ h | **built** (both vendor stage times + steady); contour-only target, locked as a regression guard |
| [18](#gw18) | 🟢 | Transient seepage through an earth fill dam | Toe-slope total head within 0.072 m rms / 0.106 m worst of the digitized Fig 20.5 profile at t = 0.6 h and 0.115 m rms / 0.187 m worst at t = 19656 h, over 11 stations — both inside the ≈0.2 m the chart can be read to · XSLOPE's steady profile within 0.074 m rms / 0.122 m worst of Fig 20.5's 19656 h curve | **built** (both vendor stage times + steady); the elastic *S*<sub>s</sub> acts in the saturated zone only, as Slide's *m*<sub>v</sub> does — the term that sets the drainage time-scale of every transient row here |
| [19](#gw19) | 🟢 | Transient seepage below a lagoon | Near-steady (11340 min) pressure head along the top boundary within 0.045 m rms / 0.126 m worst of the digitized Fig 21.9 markers over 20 stations — under 1% of the driving head · early frames 0.07–0.35 m rms, XSLOPE's mound running slightly *ahead* of RS2's | **built** (all four frames); both vendor curves are reproduced rather than fitted (Gardner conductivity through its three points, capacity band from its straight retention line), the near-steady frame is locked against the vendor profile |
| [20](#gw20) | <span class="nodata">⊘</span> | Transient seepage in a layered slope | RS2 reaches the steady query-line profile of Fig 22.7 by 208 s and XSLOPE by 400 s: the two differ by 0.105 m rms along the line at 208 s and by 0.005 m once both are steady · XSLOPE's own solved heads at 4.6 / 31 / 208 s are locked as a regression guard | **built** (all three frames); the vendor-facing lock for this geometry is [GW7](#gw7)'s steady profile, and the transient frames carry a measured ≈2× timing offset |
| [21](#gw21) | 🟢 | Transient seepage through a fully confined aquifer | Within ≈0.02 ft of the Ferris erfc closed form at 600 hr | **built** (both cases) |

</div>

## Results

### GW1: Shallow unconfined flow with rainfall {#gw1}

**Input files:** [gw001.xlsx](files/rocscience_gw/gw001.xlsx)

Flow between two long parallel rivers 10 m apart (Haar 1990; the Dupuit–Forchheimer
recharge problem), the manual's first specified-flux case. A 10 × 5 m block (base
impermeable, ground surface at el 5) with river heads h₁ = 3.75 m on the left edge and
h₂ = 3.0 m on the right, a uniform rainfall P = 2.5×10⁻⁶ m/s infiltrating the top, and
k = 1×10⁻⁵ m/s. The recharge builds an **internal** free surface that mounds above both
river levels; there is no daylighting seepage face, so the top edge is declared both the
flux boundary and an (inactive) exit face — the device that puts the solver on the
unsaturated free-surface path the internal mound needs. The free-surface position is set
by mass balance, independent of k and of the unsaturated model.

| Quantity | XSLOPE | Slide | Haar eqs 1.2–1.3 |
|---|---|---|---|
| x_a (crest position) | ≈ 4.1 | 4.06 (~0.05 m) | 3.98 (+0.12 m) |
| h_max (crest elevation) | 4.61 | 4.49 (+0.12 m) | 4.25 (+0.36 m) |
| Q (m³/s per m) | 2.500×10⁻⁵ | — | P·L = 2.5×10⁻⁵ (0.0%) |

The mound crest lands within ~0.05 m of Slide in position and ~0.12 m in height — XSLOPE
sits a touch above Slide, the same slightly-high free-surface bias the
[SEEP2D cross-check](#seep2d-crosscheck) documents across this panel; both finite-element
solutions sit above Haar's Dupuit closed form, which neglects the vertical flow the FE
models resolve near the crest. The flowrate lock is Q = P·L = 2.5×10⁻⁵ (exact by
construction, since all the rain enters and drains to the rivers); a head regression at
three interior stations guards the mound shape, which the flux total alone does not.

![gw001: mesh and solved heads](images/gw001.png)

### GW2: Flow around cylinder {#gw2}

**Input files:** [gw002.xlsx](files/rocscience_gw/gw002.xlsx)

Confined potential flow past a cylinder (Streeter's closed form; Desai & Kundu's finite
elements): a half-domain 8 × 4 m with a semicircular notch of radius 1 m at the bottom
center and fixed heads of 1.0 and 0.0 on the vertical edges. The domain is fully
saturated, and the head field in homogeneous confined flow is independent of the
conductivity, so the problem is determined by geometry and boundary heads alone.

| Point | XSLOPE | Slide | Closed form | Desai & Kundu |
|---|---|---|---|---|
| (4, 1) | 0.500 | 0.500 (0.000 m) | 0.500 (0.000 m) | 0.500 (0.000 m) |
| (4.5, 0.866) | 0.381 | 0.381 (0.000 m) | 0.375 (+0.006 m) | 0.378 (+0.003 m) |
| (5, 0) | 0.263 | 0.263 (0.000 m) | 0.250 (+0.013 m) | 0.277 (−0.014 m) |
| (6, 0) | 0.202 | 0.203 (−0.001 m) | 0.188 (+0.014 m) | 0.213 (−0.011 m) |
| (8, 0) | 0.000 | 0.000 (0.000 m) | −0.031 (+0.031 m) | 0.000 (0.000 m) |

The two finite-element solutions agree with each other within 0.0013 m everywhere; both
depart from the closed form near the downstream edge, where the analytical solution
assumes an infinite domain (its −0.031 at the finite model's zero-head boundary measures
that idealization, not an error). Flowrate is locked as a regression value.

![gw002: mesh and solved heads](images/gw002.png)

### GW3: Confined flow under a dam foundation {#gw3}

**Input files:** [gw003.xlsx](files/rocscience_gw/gw003.xlsx)

The classic flow-net problem (Rushton & Redshaw): a 40 × 10 m soil block with head 5 m
on the ground surface upstream of the dam (x = 0–8), head 0 downstream (x = 20–40), and
the dam base impervious between them.

| Station (line 1-1, y = −4) | XSLOPE | Slide |
|---|---|---|
| x = 0 | 4.47 | 4.50 (−0.03 m) |
| x = 10 | 3.40 | 3.45 (−0.05 m) |
| x = 14 | 2.44 | 2.47 (−0.03 m) |
| x = 20 | 1.05 | 1.05 (0.00 m) |
| x = 30 | 0.19 | 0.20 (−0.01 m) |
| x = 40 | 0.08 | 0.07 (+0.01 m) |

| Profile | XSLOPE | Slide |
|---|---|---|
| line 2-2 (x = 20), head span | 0.17–1.30 | 0.24–1.30 |

Both tables compare against Slide, whose markers coincide with Rushton & Redshaw's own
flow-net values everywhere on the two published profile lines, so one column carries both
sources. Heads along those lines fall within 0.08 m of the chart everywhere. The line 2-2 span is quoted as a
band rather than a delta: its lower end sits on the (20, 0) singularity discussed below,
where neither solution has a converged point value to difference.

This problem is solved on a finer mesh than the rest of the panel (`target_size` 0.1
rather than 0.4) because of the singularity at (20, 0), where the impervious dam base
meets the zero-head boundary. The head gradient is unbounded at that corner, so the
solution converges slowly in its immediate neighbourhood — the top of line 2-2 sits
right on top of it — while the rest of the domain is well resolved at any mesh size.
A coarse mesh reports a plausible number there for the wrong reason.

![gw003: mesh and solved heads](images/gw003.png)

### GW4: Steady unconfined flow through an earth dam {#gw4}

**Input files:** [gw004.xlsx](files/rocscience_gw/gw004.xlsx)

Free-surface seepage through a small dam with a **trapezoidal toe drain**, against Kozeny's
basic parabola: reservoir head 4 m on the upstream face, and a seepage exit face over the dry
upstream face, the crest, the downstream slope and the drain's contact face, with the drain toe
pinned at total head 0. The free-surface shape is independent of *k*.

**The drain contact is inclined, not vertical.** The dam's mesh boundary in the vendor model
`groundwater #004.fez` is (0,0)–(10,5)–(12.5,5)–(23.4857,1.5016)–(21.900,0): the downstream
slope runs *past* the printed dimension chain to (23.4857, 1.5016), then falls back down-left
along the drain's contact face to the toe at (21.900, 0). The manual's own figures show the same
thing — Fig 4.1 draws the trapezoidal drain that face bounds, and Fig 4.2's right-hand end puts
the seepage-face markers along the inclined contact — and the printed 10.00 + 2.50 + 9.50 = 22.00
chain measures to the *base toe*, not to the rightmost point of the slope.

That geometry decides the whole row, because **y₁ is defined at the drain**: Fig 4.5 dimensions
it as the free-surface height directly above the toe, and x₁ as the horizontal offset from the
toe to the point where the free surface enters the drain face.

| Quantity | XSLOPE | RS2 (this model) | Slide | Eqs 4.1–4.2 |
|---|---|---|---|---|
| Phreatic el. at x = 14 | 2.878 | — | — | 2.814 (+0.064 m) |
| Phreatic el. at x = 18 | 2.024 | — | — | 2.007 (+0.017 m) |
| Phreatic el. at x = 20 | 1.427 | — | — | 1.444 (−0.017 m) |
| y₁, free surface above the toe | 0.401 | 0.395 (+1.5%) | 0.442 (−9.3%) | 0.486 (−17.5%) |
| x₁, drain-face entry offset (**sets the dot**) | 0.272 | 0.226 (+20.4%) | 0.227 (+19.8%) | 0.243 (+11.9%) |
| Q (m³/s per m) | 5.487×10⁻⁸ | — | — | k·y₁ = 4.86×10⁻⁸ (+12.9%) |

The governing comparison is **RS2's own Table 4.1**, because the file this problem is built from
is RS2's, and RS2 publishes its solve of it (y₁ 0.395, x₁ 0.226). Slide's 0.442 / 0.227 is the
same model before conversion; the two vendors are 11% apart from each other on y₁, and both sit
below the Eq 4.1 idealization. XSLOPE's y₁ lands 0.006 m above RS2's.

y₁ is **corner-singular** — the free surface is nearly vertical where it enters the drain — so it
needs a fine mesh, and the earlier reading of 0.50 came from a coarse one on the wrong geometry.
Across `target_size` 0.25 / 0.147 / 0.09 / 0.06 / 0.045 / 0.03 (1.4 k → 93 k nodes) y₁ reads
0.346 / 0.377 / 0.394 / 0.401 / 0.404 / 0.407 while Q moves 0.6%; the tags run at 0.06.

x₁ is the weaker of the two, and it is the one that sets the dot: RS2's Table 4.1 publishes
both quantities, so the worse of them governs. It is a horizontal offset measured on that same
near-vertical stretch, and it depends directly on the one input this model does not carry: the
vendor material
is `gw_mode 0 / stype General`, RS2's built-in Simple curve, whose parameters the file does not
store, so XSLOPE stands a linear front (k<sub>r</sub> = 10⁻³ at 0.25 m of suction) in for it.
That assumption is a measured lever on the answer — moving the front to 0.10 m gives x₁ = 0.207
and y₁ = 0.426, and to 1.0 m gives 0.345 and 0.413 — and the builder's 0.25 m is not fitted to
either published number.

![gw004: mesh and solved heads](images/gw004.png)

### GW5: Unsaturated flow behind an embankment {#gw5}

**Input file:** [gw005.xlsx](files/rocscience_gw/gw005.xlsx)

A permeability-contrast problem taken from the FLAC manual (Coetzee et al. 1995) and
extended by the Slide manual to two materials. The domain is a 30 m × 10 m block with a 2 m
downstream shelf running out to *x* = 40; a 2 m band at 4 ≤ *y* ≤ 6 spans the block at a
saturated conductivity of 10⁻¹³ m/s against 10⁻¹⁰ m/s in the host, a factor of 1000. Total
head is 10 m on the whole left face and 4 m on the *x* = 30 step face (*y* ∈ [2, 4]) plus the
*x* = 40 end (*y* ∈ [0, 2]); every other edge is impermeable, and the published model
declares no seepage face. Geometry, conductivities and boundary cards are read from the
vendor RS2 model the manual links from its own §5.5 (`groundwater #005.fez`), so none of
this row's inputs are digitized.

Both of the vendor's conductivity functions are two-point **constant** curves, so — despite
the section title — the published model is a constant-conductivity linear problem in which
the unsaturated law never enters. The XSLOPE file reproduces that with *k*ᵣ ≡ 1.

The 1000× band isolates the block above it, which stagnates between total head 9.79 and
10.00, while the head below falls essentially linearly from 10 to 4 across the 30 m. Two
checks pin the answer, one code-independent and one against the vendor.

**Flowrate, against the one-dimensional closed form.** Almost all the flow passes through
the 4 m of host below the band, so *k·b·i* = 10⁻¹⁰ × 4 × (6/30) = 8.0×10⁻¹¹ m³/s per m.
XSLOPE reads 8.194 / 8.167 / 8.141 ×10⁻¹¹ over a tri3 1.0 m → tri3 0.5 m → tri6 0.5 m
ladder — a 0.6% spread, and +1.8% on the closed form at the finest.

**Pressure head, against Figure 5-4.** The manual's own comparison is with FLAC and is
presented as contour plates rather than tables, but Fig 5-4 carries a numeric key: pressure
head 0 → 10 m in ten 1 m bands. Read off a 300 dpi render calibrated on the domain outline
(the *x* and *y* scales agree to 0.1%), XSLOPE's solved pressure head falls inside RS2's own
band at **46 of 49** grid points spanning the whole domain including the downstream shelf.
The three exceptions miss a band edge by 0.003, 0.008 and 0.025 m — an order of magnitude
inside the ±0.5 m the 1 m banding itself resolves.

The published target is a chart, so — as the methodology note sets out for a chart target —
the lock is XSLOPE's own flowrate and total-head field, with the Fig 5-4 band comparison
above standing as the vendor check.

![gw005: mesh and solved heads](images/gw005.png)

### GW6: Steady-state seepage through saturated–unsaturated soils {#gw6}

**Input files:** [gw006a.xlsx](files/rocscience_gw/gw006a.xlsx) (case 1, isotropic) /
[gw006b.xlsx](files/rocscience_gw/gw006b.xlsx) (case 2, 9:1 anisotropy) /
[gw006c.xlsx](files/rocscience_gw/gw006c.xlsx) (case 3, core) /
[gw006e.xlsx](files/rocscience_gw/gw006e.xlsx) (case 5, seepage face)

This manual problem runs the same 12 m dam through five cases. The published target in every
case is the pressure-head profile along **line 1-1** (the crest centerline, x = 26) — a chart
curve (Figs 6.6 / 6.9 / 6.14 / 6.23) with no tabulated value and no numeric key — so, as the
methodology note sets out for a chart target, each case locks XSLOPE's own flowrate and
total-head field. The case-2/3/5 material and boundary data were read verbatim from the
vendor RS2 groundwater models (`groundwater #006_02/03/05.slw` in the RS2 Groundwater zip): case 2's `condx:1 condy:0.111111`
(kₕ = 9e-7, k_v = 1e-7), case 3's 100×-lower core over the mesh's material-2 footprint (a
rectangle x ∈ [24, 28], y ∈ [0, 10]), and case 5's crest-plus-downstream-slope seepage face.

**Case 1 — isotropic dam with a 12 m horizontal drain.**

Fredlund & Rahardjo (1993)'s saturated–unsaturated earth dam (12 m high, symmetric 2:1
faces, reservoir at 10 m, a 12 m horizontal drain at the downstream toe), case 1:
isotropic conductivity. The conductivity function is digitized from the manual's chart
(ks = 10⁻⁷ m/s, air entry ≈ 9 kPa) and fit by a Mualem–van Genuchten curve; three fit
variants move the answer by less than 3 cm, so the fit is not the controlling
uncertainty. The published target is the pressure-head profile along the crest
centerline, where Slide's and Fredlund & Rahardjo's curves coincide within 0.2 m.

| Elevation on the crest line | XSLOPE pressure head | Slide | F&R |
|---|---|---|---|
| 0 | 7.40 | 7.15 (+0.25 m) | ≈7.3 (≈+0.10 m) |
| 2 | 5.47 | 5.15 (+0.32 m) | ≈5.3 (≈+0.17 m) |
| 4 | 3.52 | 3.25 (+0.27 m) | ≈3.4 (≈+0.12 m) |
| 6 | 1.66 | 1.30 (+0.36 m) | ≈1.45 (≈+0.21 m) |
| 8 | −0.22 | −0.60 (+0.38 m) | ≈−0.45 (≈+0.23 m) |

*The profile shape reproduces exactly; the whole curve sits 0.25–0.5 m above the
published pair, insensitive to the conductivity fit and to mesh refinement. This is a
family difference, not a solver error: run against SEEP2D (the original USACE/WES code)
on the identical mesh, boundary conditions and unsaturated law, XSLOPE's free surface
daylights at the same place — see [the SEEP2D cross-check](#seep2d-crosscheck) below. The
regression locks XSLOPE's own values.*

![gw006a: mesh and solved heads](images/gw006a.png)

**Case 2 — anisotropic dam (kₕ = 9 kᵥ) with the horizontal drain.** The horizontal
conductivity is nine times the vertical (vendor `condx:1 condy:0.111111`, so kₕ = 9×10⁻⁷,
kᵥ = 10⁻⁷ m/s), which spreads the flow and lowers the phreatic surface; the unsaturated
relative shape is identical to case 1, so the same Mualem–vG fit is reused. Here XSLOPE
reproduces the published curve **almost exactly**. Fredlund & Rahardjo's curve is
indistinguishable from Slide's on Fig 6.9 at the scale the chart can be read to, so the
comparison column carries the one reading:

| Elevation on line 1-1 | XSLOPE pressure head | Slide (Fig 6.9) |
|---|---|---|
| 0 | 6.52 | ≈6.5 (≈+0.02 m) |
| 2 | 4.74 | ≈4.7 (≈+0.04 m) |
| 4 | 3.19 | ≈3.2 (≈−0.01 m) |
| 6 | 1.79 | ≈1.85 (≈−0.06 m) |
| 8 | 0.42 | ≈0.4 (≈+0.02 m) |

Flowrate 1.639×10⁻⁶ m³/s per m (locked with the total-head field).

![gw006b: mesh and solved heads (9:1 anisotropy)](images/gw006b.png)

**Case 3 — isotropic dam with a low-permeability core and the horizontal drain.** A
rectangular central core (x ∈ [24, 28], y ∈ [0, 10], read from the vendor mesh's material-2
footprint) with saturated k = 10⁻⁹ m/s — 100× lower than the 10⁻⁷ shell — is tiled into the
dam as four non-overlapping polygons (three shell pieces + the core). The core forces almost
the whole head drop across its 4 m width (Fig 6.13's crowded contours), throttling the
flowrate to 2.502×10⁻⁸ m³/s per m. Along line 1-1 (now inside the core) XSLOPE reproduces the
profile shape and sits at the high end of the published scatter:

| Elevation on line 1-1 | XSLOPE pressure head | Slide (Fig 6.14) | Ref[1] |
|---|---|---|---|
| 0 | 5.55 | ≈5.9 (≈−0.35 m) | ≈5.8 (≈−0.25 m) |
| 2 | 4.06 | ≈3.9 (≈+0.16 m) | ≈3.9 (≈+0.16 m) |
| 4 | 2.68 | ≈2.1 (≈+0.58 m) | ≈2.1 (≈+0.58 m) |
| 6 | 1.48 | ≈0.4 (≈+1.08 m) | ≈0.7 (≈+0.78 m) |
| 8 | 0.28 | ≈−1.2 (≈+1.48 m) | ≈−0.3 (≈+0.58 m) |

*The published Slide and Ref[1] curves themselves diverge ~1.5 m near the crest; XSLOPE
tracks the shape and sits above both up high — the same +0.5 m free-surface family as case 1.
Locked at XSLOPE's own values.*

![gw006c: mesh and solved heads (low-k core)](images/gw006c.png)

**Case 5 — isotropic dam with a downstream seepage face (no drain).** The horizontal toe
drain is replaced by the "unknown boundary condition": the crest and the whole downstream
slope are a seepage face where the phreatic surface may daylight (vendor seepage-face nodes
(22,11)→(50,1)). Without the drain the phreatic surface rides higher. XSLOPE again reproduces
the published curve **almost exactly**. As on case 2, Fredlund & Rahardjo's curve and
Slide's coincide on Fig 6.23 within the chart's read precision, so one column carries both:

| Elevation on line 1-1 | XSLOPE pressure head | Slide (Fig 6.23) |
|---|---|---|
| 0 | 8.29 | ≈8.4 (≈−0.11 m) |
| 2 | 6.35 | ≈6.4 (≈−0.05 m) |
| 4 | 4.39 | ≈4.5 (≈−0.11 m) |
| 6 | 2.45 | ≈2.5 (≈−0.05 m) |
| 8 | 0.50 | ≈0.55 (≈−0.05 m) |

Flowrate 1.686×10⁻⁷ m³/s per m (locked with the total-head field).

![gw006e: mesh and solved heads (seepage face)](images/gw006e.png)

*Case 4 (steady-state infiltration, a 10⁻⁸ m/s flux over the whole dam surface) is not built:
a surface flux boundary combined with a downstream exit face does not converge.*

### GW7: Seepage within a layered slope {#gw7}

**Input file:** [gw007.xlsx](files/rocscience_gw/gw007.xlsx)

Rulon & Freeze's sandbox slope (after Fredlund & Rahardjo 1993): a medium-sand slope
holding a thin fine-sand lens. Geometry from the vendor RS2 model (a 2.4 × 1.0 m box, a
2:1 downstream slope from the crest at (1.6, 1.0) to the toe at (0, 0.2), and the fine
lens as the y = 0.6–0.7 band from the slope face to the right wall); conductivity curves
from the manual's Fig 7.2, fit here by Mualem–van Genuchten (medium ks = 1.4×10⁻³ m/s,
fine ks = 5.5×10⁻⁵ m/s). A uniform 2.1×10⁻⁴ m/s falls on the crest; the submerged toe is
a tailwater head of 0.3 m; the rest of the slope is a seepage face; base and right wall
are no-flow.

The infiltration rate exceeds the fine sand's saturated conductivity, so the recharge
cannot drain vertically through the lens: it **perches** a water table on top of the fine
band and sheds laterally to daylight as a slope-face spring — the physical result Rulon &
Freeze observed. XSLOPE reproduces it: the perched saturated zone above the lens (the
crowded head contours across the low-k band), the free surface daylighting on the slope,
and the main water table exiting at el 0.30 at the toe — the "0.3 m from the toe" the
manual reports.

| Quantity | XSLOPE | Slide (after Rulon & Freeze) | Mass balance on the flux boundary |
|---|---|---|---|
| water table at the toe | el 0.30 | 0.3 m stated (0.00 m) | — |
| Q, m³/s per m | 1.680×10⁻⁴ | — | q·L = 1.68×10⁻⁴ (0.0%) |

Problem 7's own published targets are all chart curves — the water table (Fig 7.4) and
the total-head profiles along lines 1-1 and 2-2 (Figs 7.7, 7.8) — with no tabulated
number. The conductivity fit shifts the unsaturated detail but not the flowrate or the
perched-table topology.

**The steady profile has a numeric target, in a later chapter.** Manual problem 22
([GW20](#gw20)) re-runs this identical slope transiently, and its Fig 22.7 plots total
head along a vertical query line at *x* = 1.6 — the crest break, from the base up to the
crest — with marker values for RS2 at 4.6, 31 and 208 s. The 208 s frame has reached
steady state, which is this problem. The chart also calibrates itself: its 4.6 s markers
at the base read 0.302–0.305 against the model's initial total head of 0.300, so the
digitization is good to ≈0.005 m. Along that query line XSLOPE's steady field matches
RS2's markers to **0.005 m rms and 0.013 m at worst over 21 stations** — about 1% of the
0.56 m head range along the profile — and puts the step across the fine lens in the right
place:

| $y$ on the query line ($x=1.6$) | 1.00 | 0.85 | 0.70 | 0.65 | 0.60 | 0.30 | 0.00 |
|---|---|---|---|---|---|---|---|
| XSLOPE total head | 0.868 | 0.850 | 0.827 | 0.743 | 0.632 | 0.614 | 0.607 |
| Fig 22.7, RS2 at 208 s | 0.862 | 0.845 | 0.832 | 0.730 | 0.628 | 0.610 | 0.605 |

The same chart's Ref [1] (Fredlund & Rahardjo) curve runs a little above both — ≈0.878 near
the crest and ≈0.617 near the base, against RS2's 0.862 and 0.605 — so XSLOPE's 0.868 and
0.607 sit between the two sources, closer to RS2. Five stations of the profile are locked
alongside the flowrate (Q = q·L, exact by construction on the flux boundary) and the
three-station regression that guards the solved field.

<!-- test: file=files/rocscience_gw/gw007.xlsx, type=seep, target_size=0.04, max_iter=1000, expected_flowrate=1.680e-04, tolerance=0.02, benchmark=GW7-q -->
<!-- test: file=files/rocscience_gw/gw007.xlsx, type=seep_head, target_size=0.04, max_iter=1000, points=1:0.1:0.518;2:0.1:0.642;2.2:0.3:0.657, tolerance=0.02, benchmark=GW7-h -->
<!-- test: file=files/rocscience_gw/gw007.xlsx, type=seep_head, target_size=0.04, max_iter=1000, points=1.6:0.05:0.608;1.6:0.3:0.614;1.6:0.6:0.632;1.6:0.8:0.844;1.6:1:0.868, tolerance=0.02, benchmark=GW7-q22 -->

![gw007: mesh and solved heads](images/gw007.png)

### GW8: Flow through ditch-drained soils {#gw8}

**Input file:** [gw008.xlsx](files/rocscience_gw/gw008.xlsx)

A ditch-drained two-layer aquifer after Gureghian (1981), and the corpus' exercise of the
[specified-flux boundary](#flux-crosscheck) — the whole problem is driven by rainfall
infiltration on the ground surface, so it cannot be posed at all without one. A half-drain
spacing of 1.0 m over an impermeable base at 0.5 m depth; a coarse Soil A in the lowest 0.1 m
(*k* = 1.111111×10⁻³ m/s, Gardner *a* = 1000, *n* = 4.5) beneath a finer Soil B
(*k* = 1.111111×10⁻⁴ m/s, *a* = 277.777, *n* = 4.2); a uniform infiltration of
4.4444×10⁻⁵ m/s on the top; the water-free ditch as a seepage face on the left wall; base and
right-hand symmetry edge no-flow.

**Two inputs are taken from the vendor model rather than the printed tables**, because the two
disagree. The manual's text prints the infiltration as 4.4×10⁻⁶ m/s and its Table 8.1 prints
Soil B's Gardner *a* as 2777.7; the Slide model converted to `groundwater #008.fez` carries
4.4444×10⁻⁵ m/s on all twenty top-boundary flux edges and *a* = 277.777. The flux card's units
are fixed by two problems whose printed value is not in dispute — GW1 carries 2.5×10⁻⁶ against
a printed *P* = 2.5×10⁻⁶ m/s, GW7 carries 2.1×10⁻⁴ against a printed 2.1×10⁻⁴ m/s — so the card
is the distributed flux in m/s. The vendor set is also the arithmetically clean one:
4.4444×10⁻⁵ m/s is 0.16 m/hr and *q*/*k*<sub>B</sub> = 0.400 exactly, where the printed value is
clean in nothing.

Both substitutions are measured, not assumed. A two-layer Dupuit calculation at the symmetry
divide gives a mound of 0.063 m from the printed flux and 0.240 m from the vendor's, and both
published figures draw the water table at ≈0.25 m there. Gardner *a* is discriminated on Fig
8.3's own labeled pressure-head contours, digitized at five stations: *a* = 277.777 reproduces
all fourteen station values to 0.010 m rms (0.019 m worst), while *a* = 2777.7 never reaches
−0.16 m of suction anywhere in the section and therefore cannot draw the −0.17 m and −0.20 m
contours the figure shows at all.

**The flux boundary behaves exactly.** Total applied inflow is *q*·*L* = 4.44440×10⁻⁵ at every
mesh size tested, from 243 to 14 867 nodes; the confined form of the same model (fix the base,
remove the seepage face) produces a head rise of 0.163998 m against the one-dimensional hand
calculation *q*·(0.4/*k*<sub>B</sub> + 0.1/*k*<sub>A</sub>) = 0.163998 m, agreeing to six
figures. Reversing the sign of *q* inverts the response antisymmetrically.

**The published water table and pressure-head contours are reproduced.**

| Station (*x*) | XSLOPE water table | Slide (after Gureghian), Figs 8.3, 8.4 |
|---|---|---|
| 0.02 | 0.045 | 0.041 (+0.004 m) |
| 0.25 | 0.130 | 0.125 (+0.005 m) |
| 0.50 | 0.201 | 0.197 (+0.005 m) |
| 0.75 | 0.240 | 0.234 (+0.006 m) |
| 0.98 | 0.253 | 0.247 (+0.006 m) |

| Quantity | XSLOPE | Slide (after Gureghian) |
|---|---|---|
| water table at the symmetry edge | 0.253 m | ≈0.25 m |
| pressure head at 14 stations on Fig 8.3's −0.10 to −0.20 m contours | within 0.010 m rms / 0.019 m worst of the drawn contours | — |
| total head, top of the symmetry edge | 0.347 m | between the figure's −0.14 and −0.17 m suction contours |

The solved water table sits a uniform 0.004–0.006 m above the digitized line over the whole
1.0 m span — 2.4% of the 0.25 m divide mound at its largest, the same slightly-high free-surface
bias the [SEEP2D cross-check](#seep2d-crosscheck) documents across this panel. The flowrate lock
(*q*·*L*, exact by construction) is joined by a five-station head lock; heads move less than
0.003 m across a 0.05 → 0.0125 m mesh refinement.

<!-- test: file=files/rocscience_gw/gw008.xlsx, type=seep, target_size=0.025, element_type=tri3, expected_flowrate=4.4444e-05, tolerance=0.02, benchmark=GW8-q -->
<!-- test: file=files/rocscience_gw/gw008.xlsx, type=seep_head, target_size=0.025, element_type=tri3, points=0.25:0.25:0.183;0.5:0.25:0.220;0.75:0.25:0.243;1:0.25:0.253;1:0.5:0.347, tolerance=0.01, benchmark=GW8-h -->

### GW9: Seepage through dam {#gw9}

**Input files:** [gw009a.xlsx](files/rocscience_gw/gw009a.xlsx) (dam 1) ·
[gw009b.xlsx](files/rocscience_gw/gw009b.xlsx) (dam 2, toe drain)

Bowles' homogeneous dam, the flow-net textbook example re-solved numerically by Chapuis,
Chenaf & Bowles (2001) and by Slide: base 100 m, crest 10 m at el. 20 (2.5:1 upstream,
2:1 downstream), reservoir head 18.5 m, ks = 6.67×10⁻⁶ m/s with the manual's printed
8-point unsaturated conductivity table, fit here by a Mualem–van Genuchten curve
(α = 0.2835, n = 2.765).

| | XSLOPE | Slide | SEEP/W (fine) | Bowles (flow nets) |
|---|---|---|---|---|
| Q, m³/(min·m) | 1.379×10⁻³ | 1.378×10⁻³ (+0.1%) | 1.37×10⁻³ (+0.7%) | 1.10–1.28×10⁻³ |

**Dam 2 — Bowles' dam with a toe drain** (Bowles 1984, Example 9-2 / Fig E9-2b, p. 248;
Slide manual §9.2, Fig 9.5; Chapuis et al. 2001, Fig 5). Base 190 m; crest 10 m wide at
el. 45; symmetric 2:1 faces (upstream and downstream horizontal runs 90 m each); reservoir
head 40 m. A coarse toe drain (ks = 1.0×10⁻⁴ m/s) fills the downstream-toe triangle
(100, 0)–(190, 0)–(145, 22.5) — base 90 m, apex at mid-height of the downstream slope. The
body's saturated conductivity is ks = 2.0×10⁻⁷ m/s, carrying the dam-1 unsaturated k(u) curve.

| | XSLOPE | Slide | SEEP/W (2328 el.) | Bowles (flow net) |
|---|---|---|---|---|
| Q, m³/(s·m) | 4.29×10⁻⁶ | 4.23×10⁻⁶ (+1.4%) | 4.23×10⁻⁶ (+1.4%) | 3.8×10⁻⁶ (+12.9%) |

XSLOPE matches the two numerical benchmarks to 1.4% and Bowles' flow net to its graphical
accuracy. Note the units: dam 2 is worked per **second** (Bowles solves it in cm/s), where
dam 1 was per minute.

*Provenance.* The body conductivity comes straight from Bowles (1984), Example 9-2, Fig E9-2b:
the figure prints k = 2×10⁻⁵ cm/s = 2.0×10⁻⁷ m/s, and its flow net (n_f/n_d = 1.9/4) gives
q = k·h·n_f/n_d = 2.0×10⁻⁷ × 40 × 0.475 = 3.8×10⁻⁶ m³/(s·m) as an independent hand-check.
(Bowles' printed answer, 3.8×10⁻² m³/(s·m), applies the cm→m factor the wrong way in that same
line; the conductivity, head, and flow net are unambiguous.) Two errata in the secondary
sources are resolved by the figure:

- **Body k.** The 2.0×10⁻⁶ m/s in the Chapuis et al. (2001) Fig 5 caption is one decade high —
  a −6/−7 exponent slip. Bowles' value is 2.0×10⁻⁷ m/s, and the Slide manual's own Fig 9.6
  chart draws the earth-dam curve at ≈2×10⁻⁷.
- **Flowrate units.** The published Q, tabulated as m³/(min·m) beside dam 1, is actually
  m³/(s·m). Bowles 3.8×10⁻⁶ and Chapuis's SEEP/W and Slide 4.23×10⁻⁶ are all per second, and
  agree because all three used k ≈ 2×10⁻⁷ m/s.

Run at the caption's 2.0×10⁻⁶ m/s, XSLOPE returns 3.97×10⁻⁵ m³/(s·m) — an order of magnitude
above the published value (Q is nearly linear in k), which makes the exponent slip visible. The
locked value is XSLOPE's own Q at Bowles' conductivity, 4.29×10⁻⁶ m³/(s·m).

*Vendor check.* The RS2 Groundwater Verification set ships only `groundwater #009_01.fez`
(dam 1, no drain); the Slide manual's `Groundwater#09_2.sli` is not in the distributed model
set. Neither is needed — Bowles (1984) fixes the conductivity directly.

![gw009a: mesh and solved heads](images/gw009a.png)

![gw009b: mesh and solved heads](images/gw009b.png)

### GW10: Steady unconfined flow, van Genuchten permeability {#gw10}

**Input files:** [gw010.xlsx](files/rocscience_gw/gw010.xlsx)

Clement, Wise, Molz & Wen (1996)'s unconfined square domain, the manual's designated
van Genuchten test: a 10 × 10 m block with head 10 on the left edge, tailwater 2 on the
right, an exit face above the tailwater, and vG conductivity (α = 0.64, n = 4.65,
ks = 1.1574×10⁻⁵ m/s) — an exact capability match for the solver's `vg` option.

| | XSLOPE | Slide | Clement et al. |
|---|---|---|---|
| Q (m³/s per m) | 6.070×10⁻⁵ | 6.066×10⁻⁵ (+0.1%) | 6.076×10⁻⁵ (−0.1%) |
| Phreatic exit elevation | 4.87 | 5.0 (−0.13 m) | 4.8 (+0.07 m) |

The manual's "seepage face" column tabulates the phreatic exit *elevation* (both
published figures show the free surface exiting near el. 4.9–5.0), not a face length.
Only the tailwater-2 case carries published numbers and is locked.

![gw010: mesh and solved heads](images/gw010.png)

### GW11: Earth/rock-fill dam, Gardner permeability function {#gw11}

**Input file:** [gw011.xlsx](files/rocscience_gw/gw011.xlsx)

A 45 m homogeneous earth/rock-fill dam after Zhang et al. (2001), and the manual's dedicated
test of the **Gardner** relative-conductivity function, $k_r = 1/(1 + a\,\psi^n)$, with the
published $a$ = 0.15 and $n$ = 6. Crest 17 m, upstream run 89.1 m, downstream run 76.9 m,
impermeable base; reservoir at el. 40 m, no tailwater; the whole downstream slope is an exit
face. The published quantity is the **release point** — the elevation at which the free surface
daylights on the downstream face.

| | release point |
|---|---|
| XSLOPE (Gardner `unsat=gard`, `target_size` = 1.0) | 17.90 m |
| Slide | 19.397 m (−1.50 m) |
| ABAQUS (Zhang et al.) | 19.64 m (−1.74 m) |

*XSLOPE releases about 1.5 m low, and the difference is unexplained.*

**The free-surface field itself agrees much more closely than the release point does.** Fig 11.2
draws Slide's own free surface across the whole dam; digitized at 18 stations from x = 80 to
x = 148 (calibrated independently — the drawn line meets the upstream face at el. 40.0 exactly
where the reservoir does), XSLOPE tracks it to **0.76 m rms and 1.23 m worst**, 1.7% and 2.7% of
the 45 m dam. XSLOPE sits ~0.4 m *above* Slide's line near the reservoir and crosses to ~1.2 m
below it through the downstream half, so the gap is a difference in the surface's slope across
the dam body that accumulates into the release point, not a discontinuity at the exit.

Four candidate causes have been measured and none of them is it:

| candidate | test | release point |
|---|---|---|
| unsaturated law | Gardner (the manual's own), van Genuchten, linear front | 17.4–17.9 m |
| mesh size | `target_size` 2.0 / 1.0 / 0.6 (1.3 k → 14.5 k nodes) | 17.79 / 17.90 / 17.76 m |
| element order | tri6 at the vendor's own element count (794 vs its 865), and tri6 at `target_size` 1.0 | 17.81 / 17.90 m |
| relative-conductivity floor | $k_{r,\min}$ = 10⁻⁴ / 10⁻⁶ / 10⁻⁸ / 10⁻¹² | 17.90 m at every value |
| exit-face extent | adding the dry upstream face and the crest, and pinning the vendor's total head 0 at the toe (183, 0) | 17.90 m |

The element-order row matters because the vendor meshes this — and only this — groundwater
problem with 6-noded quadratic triangles, where XSLOPE takes $k_r$ per Gauss point rather than
from the element-average pressure; running the same element order moves the answer by 0.09 m.
The $k_{r,\min}$ row is a controlled override of the solver's internal floor, not an input
setting: no public call site passes `kr_min`, and the floor turns out to be inert here because
it binds only from 6.4 m of suction upward, where the conveyance is already negligible. There is
also no trend with refinement, and at the finest mesh the next face node *above* the release
point sits at 18.06 m — still 1.3 m below Slide, so the gap is not a nodal-resolution artifact.

Two notes on reading the manual. Its text gives the downstream slope as 1:1.171, but the
printed dimensions on Fig. 11.1 give 76.9/45 = 1.709 — a digit transposition. The figure's
dimensions are used here (they also agree with a direct measurement of the plate), and they are
the *more favourable* read: the text's slope releases at 16.54 m, lower still. Separately,
$k_s$ is not printed for this case; because the dam is homogeneous the free-surface position is
independent of $k_s$, so it does not affect the published target, but it does set the discharge —
which is why the regression tag below locks XSLOPE's own flowrate rather than a published one.

GW11 is the single exit-face problem in this panel where XSLOPE releases low. The
[SEEP2D cross-check](#seep2d-crosscheck) below puts XSLOPE on the *same* release point as the
original USACE code on every other exit-face problem, so the difference is specific to this
problem rather than a family-wide convention difference. Case 2 of the manual's
problem (the zoned dam with a foundation and toe drain) is not built.

<!-- test: file=files/rocscience_gw/gw011.xlsx, type=seep, target_size=1.0, max_iter=2000, expected_flowrate=7.814e-07, tolerance=0.05, benchmark=GW11-q -->
<!-- test: file=files/rocscience_gw/gw011.xlsx, type=seep_head, target_size=1.0, max_iter=2000, points=60:20:39.09;95:20:33.76;120:20:27.84;140:15:21.87;155:10:15.55, tolerance=0.05, benchmark=GW11-h -->

### GW12 / GW13: Ditch seepage into a deep drainage layer (Vedernikov) {#gw12}

**Input files:** [gw012.xlsx](files/rocscience_gw/gw012.xlsx) (trapezoidal) ·
[gw013.xlsx](files/rocscience_gw/gw013.xlsx) (triangular)

Vedernikov's closed-form solutions for steady seepage from a channel into a deep
drainage layer, modeled as half-domains by symmetry with the ditch perimeter at head 50
and the deep drain as head 0 on the base. The seepage detaches below the ditch and
descends as a bulb whose width the theory predicts.

| | XSLOPE | Slide | Vedernikov |
|---|---|---|---|
| Trapezoidal: Q per half | 4.137×10⁻⁴ | 4.093×10⁻⁴ (+1.1%) | 4.0×10⁻⁴ (+3.4%) |
| Trapezoidal: bulb half-width | ≈42 | 41 (+1) | 40 (+2) |
| Triangular: Q per half | 2.087×10⁻² | 2.050×10⁻² (+1.8%) | 2.0×10⁻² (+4.4%) |

The detached-bulb iteration converges cleanly at 1,500 free-surface iterations (the
`max_iter` tag key); the default 400 is not enough for these geometries.

![gw012: mesh and solved heads](images/gw012.png)

![gw013: mesh and solved heads](images/gw013.png)

### GW14: Unsaturated soil column {#gw14}

The manual's steady capillary-profile problem, after
[Gardner (1958)](https://doi.org/10.1097/00010694-195804000-00006): a soil column of length
L = 1 m with ks = 10⁻⁷ m/s and α = 1 m⁻¹, carrying a steady vertical flux
ν = ±8.64×10⁻⁴ m/d = ±10⁻⁸ m/s — evaporation drawing water up, or infiltration pushing it
down — whose steady suction profile Gardner solves in closed form.

The problem is blocked because the published profile and XSLOPE's conductivity model are not
the same function. Gardner's analytical profile is derived for the **exponential**
conductivity law k = ks·e^(αψ), and the vendor RS2 model is configured to match it
(`conductivity: Gardner, α = 1`). XSLOPE does not implement that law: its `gard` option is
the power form kr = 1/(1 + a·ψⁿ), the Gardner-named function that SEEP/W and Slide carry,
which is a different function from the exponential one the closed form assumes. Fitting the
one to the other would change the quantity under test rather than verify it.

Nor is there a vendor number to fall back on. The manual publishes only the Fig 14.3/14.4
charts for this problem, with no tabulated Slide value, so there is nothing numeric to
compare a differently-shaped profile against. What remains lockable is a 1-D through-flux,
which the [flux cross-check](#flux-crosscheck) already verifies to machine precision.

## Transient problems {#transient}

The three transient problems below exercise XSLOPE's uncoupled
[transient seepage solver](../seep/transient.md). Each has a **closed-form** target (Terzaghi, Ferris) or a **recomputed
analytical series** (Pyrah's two-layer consolidation), so the lock is the analytical value
itself and the tolerance only absorbs the numerical (mesh + backward-Euler) error, which is
reported for each.

All three are modelled as **saturated** columns/strips. The excess pore pressure (or aquifer
head rise) is carried on an arbitrary datum offset — a constant baseline head $h_\text{ref}$
added to every node — so the pressure head stays positive everywhere and the storage is $S_s$
throughout. The offset cancels out of the excess head $h-h_\text{ref}$, which is the physical
quantity; it only selects the solver's saturated branch, where the governing equation is exactly
the linear diffusion $\partial h/\partial t = (K/S_s)\nabla^2 h = c_v\nabla^2 h$. The uniform
non-steady initial condition is set with a repeated-time **step series**: the boundary head
holds the initial value at $t=0$ (so the solver's $t=0$ steady solve produces the uniform IC)
and steps to the drained value for $t>0$. No initial-field input is needed — the `tseep`
sheet alone defines the transient model.

### GW15: Terzaghi 1-D consolidation {#gw15}

**Input files:** [gw015a.xlsx](files/rocscience_gw/gw015a.xlsx) (double drainage) ·
[gw015b.xlsx](files/rocscience_gw/gw015b.xlsx) (single drainage)

The manual's designated consolidation test: a 1 m clay column with a uniform initial excess
pore pressure $u_0$ dissipating by one-dimensional flow. Case 1 drains at **both** faces
(drainage path $H=0.5$ m); case 2 drains at the **top only** over an impermeable base
($H=1.0$ m). The excess pore pressure follows Terzaghi's (1943) Eq 17.3,

$$ \frac{u_e}{u_0}=\sum_{m=0}^{\infty}\frac{2}{M}\sin\!\left(M Z\right)e^{-M^2 T_v},
\qquad M=\tfrac{\pi}{2}(2m+1),\quad T_v=\frac{c_v t}{H^2}, $$

with the normalized depth $Z=z/H$ measured from a drained face. The material is
$m_v=0.01\ \text{kPa}^{-1}$, $k=10^{-5}$ m/s, giving $S_s=\gamma_w m_v=0.0981\ \text{m}^{-1}$ and
$c_v=k/S_s=1.02\times10^{-4}\ \text{m}^2/\text{s}$ (the published target is dimensionless, so
$c_v$ only sets the real-time scale).

XSLOPE reproduces the isochrones (below) at every interior sample point and time:

| Case | Drainage path $H$ | Δ vs Terzaghi Eq 17.3 (max over the sampled depths and times) |
|---|---|---|
| 1 — drained at both faces | 0.5 m | 0.73% of $u_0$ |
| 2 — drained at the top only | 1.0 m | 0.47% of $u_0$ |

Both maxima fall at the **earliest** sampled time, where the isochrone is steepest and a
uniform initial condition is hardest to resolve; by the last sampled time the two cases are
within 0.27% and 0.24% of $u_0$.

The tags lock the closed-form total head
$h_\text{ref}+u_0\,(u_e/u_0)$ at three depths and two time factors per case, tolerance 0.6 m.

![gw015: Terzaghi isochrones, analytical vs XSLOPE](images/gw015.png)

<!-- test: file=files/rocscience_gw/gw015a.xlsx, type=tseep_head, target_size=0.02, time=500, points=0.125:0.25:154.766;0.125:0.5:176.533;0.125:0.75:154.766, tolerance=0.6, benchmark=GW15a-t500 -->
<!-- test: file=files/rocscience_gw/gw015a.xlsx, type=tseep_head, target_size=0.02, time=1000, points=0.125:0.25:132.924;0.125:0.5:146.551;0.125:0.75:132.924, tolerance=0.6, benchmark=GW15a-t1000 -->
<!-- test: file=files/rocscience_gw/gw015b.xlsx, type=tseep_head, target_size=0.02, time=2000, points=0.125:0.25:170.955;0.125:0.5:154.766;0.125:0.75:129.887, tolerance=0.6, benchmark=GW15b-t2000 -->
<!-- test: file=files/rocscience_gw/gw015b.xlsx, type=tseep_head, target_size=0.02, time=4000, points=0.125:0.25:143.010;0.125:0.5:132.924;0.125:0.75:117.821, tolerance=0.6, benchmark=GW15b-t4000 -->

### GW16: Pore-pressure dissipation in stratified soil {#gw16}

**Input files:** [gw016a.xlsx](files/rocscience_gw/gw016a.xlsx) (uniform) ·
[gw016b.xlsx](files/rocscience_gw/gw016b.xlsx) (A over B) ·
[gw016c.xlsx](files/rocscience_gw/gw016c.xlsx) (B over A)

Pyrah's (1996) two-layer consolidation column: a 1 m column, drained at the top and
impermeable at the base, of two 0.5 m layers with **the same coefficient of consolidation**
$c_v=1$ but different conductivity and compressibility — Soil A ($k=1$, $m_v=1$) and Soil B
($k=10$, $m_v=10$), in consistent units with $\gamma_w=1$. Because $c_v$ is uniform, the
dissipation is shaped entirely by the $k/m_v$ contrast at the layer interface, which is Pyrah's
point: **layer order matters even when $c_v$ does not.** Case 1 is a uniform column (a Terzaghi
check); cases 2 and 3 swap the layer order.

The analytical solution is a **recomputed eigenfunction series** rather than a digitized chart.
With equal $c_v$ the excess pore pressure satisfies $\partial u_e/\partial t=\partial^2
u_e/\partial z^2$ in both layers, so $u_e=\sum_n c_n Z_n(z)\,e^{-\beta_n^2 t}$ where the spatial
eigenfunctions $Z_n$ satisfy $Z'(0)=0$ (impermeable base), $Z(L)=0$ (drained top), and
continuity of head **and of Darcy flux** ($k\,\mathrm{d}Z/\mathrm{d}z$) at the interface. For
two equal-thickness layers this collapses to the eigenvalue equation
$k_\text{bot}\tan^2(\beta/2)=k_\text{top}$; the coefficients $c_n$ come from the
$S_s$-weighted projection of the uniform initial condition (all integrals in closed form). The
builder recomputes this series directly.

The layer-order labels follow the **upper/lower** reading ("A/B" = A on top). The eigen-series
and the solver use the identical assignment, so the comparison is order-independent; only the
match to Pyrah's Fig 18-3/18-4 depends on the convention.

XSLOPE tracks the recomputed series across all three cases:

| Case | Layer order | Δ vs the recomputed eigenfunction series (max over the sampled depths and times) |
|---|---|---|
| 1 | uniform column | 0.19% of $u_0$ |
| 2 | Soil A over Soil B | 0.52% of $u_0$ |
| 3 | Soil B over Soil A | 0.47% of $u_0$ |

As in [GW15](#gw15) each maximum falls at the earliest sampled time; by the last one the three
cases are within 0.13%, 0.31% and 0.14% of $u_0$.

The interface kink and the strong effect of layer order are clear in
the isochrones: with the low-permeability Soil A **on top** (case 2) the underlying Soil B stays
near its initial pressure far longer, while high-permeability Soil B on top (case 3) drains the
upper layer quickly. The tags lock the closed-form head at three depths and two times per case,
tolerance 4–6 m (of $u_0=1000$), with small time steps (`max_head_change_frac=0.005`) since the
residual error at the interface is temporal.

![gw016: Pyrah two-layer isochrones, recomputed series vs XSLOPE](images/gw016.png)

<!-- test: file=files/rocscience_gw/gw016a.xlsx, type=tseep_head, target_size=0.02, time=0.2, max_head_change_frac=0.005, points=0.25:0.25:816.227;0.25:0.5:653.176;0.25:0.75:402.084, tolerance=4.0, benchmark=GW16a-t0.2 -->
<!-- test: file=files/rocscience_gw/gw016a.xlsx, type=tseep_head, target_size=0.02, time=0.5, max_head_change_frac=0.005, points=0.25:0.25:442.557;0.25:0.5:362.188;0.25:0.75:241.899, tolerance=4.0, benchmark=GW16a-t0.5 -->
<!-- test: file=files/rocscience_gw/gw016b.xlsx, type=tseep_head, target_size=0.02, time=0.2, max_head_change_frac=0.005, points=0.25:0.25:1046.604;0.25:0.5:1013.431;0.25:0.75:562.623, tolerance=6.0, benchmark=GW16b-t0.2 -->
<!-- test: file=files/rocscience_gw/gw016b.xlsx, type=tseep_head, target_size=0.02, time=0.5, max_head_change_frac=0.005, points=0.25:0.25:945.851;0.25:0.5:916.037;0.25:0.75:512.850, tolerance=6.0, benchmark=GW16b-t0.5 -->
<!-- test: file=files/rocscience_gw/gw016c.xlsx, type=tseep_head, target_size=0.02, time=0.2, max_head_change_frac=0.005, points=0.25:0.25:601.928;0.25:0.5:340.125;0.25:0.75:255.692, tolerance=6.0, benchmark=GW16c-t0.2 -->
<!-- test: file=files/rocscience_gw/gw016c.xlsx, type=tseep_head, target_size=0.02, time=0.5, max_head_change_frac=0.005, points=0.25:0.25:181.529;0.25:0.5:131.238;0.25:0.75:119.462, tolerance=6.0, benchmark=GW16c-t0.5 -->

### GW17: Transient seepage through an earth fill dam with a toe drain {#gw17}

**Input files:** [gw017.xlsx](files/rocscience_gw/gw017.xlsx)

The same dam as [GW18](#gw18) but **with a 12 m toe drain** — a 0.5 m deep high-k strip
under the downstream toe ($x\in[40,52]$, $y\in[-0.5,0]$) held at total head 0 (the vendor's
$t_x=0$ drain nodes). The drain draws the phreatic surface down so the downstream slope
stays largely unsaturated and the flow concentrates into the drain. The reservoir is raised
4 m → 10 m at $t=0$ as in GW18. Every hydraulic input comes from the vendor model
`groundwater #019`: $k_s=10^{-7}$ m/s for the fill and $10^{-2}$ m/s for the drain,
$S_s=\gamma_w m_v = 10\times0.003=0.03\ \text{m}^{-1}$ (drain $10\times0.002$),
$S_y=0.4-0.1=0.3$ from its water-content table, and a Mualem–van Genuchten fit
($\alpha=0.232\ \text{m}^{-1}$, $n=2.93$) to its 5-point conductivity table. Unlike
[GW18](#gw18) this model's retention range (100 kPa ≈ 10 m of head) is the same order as
its conductivity curve, so one van Genuchten pair can carry both, and it fits the
conductivity table twice as well as a Gardner law would (0.055 against 0.113 decades rms
over the 0–20 m of suction the dam reaches).

The published targets are total-head and pressure-head **contours** at 15 h and 16383 h
(Figs 19-4…19-7, vs FlexPDE and SEEP/W, Pentland et al. 2001) — chart-only, no tabulated
profile. Both vendor stage times are solved and locked, along with XSLOPE's own steady
field; the locks are regression guards on XSLOPE's own values.

**The early frame now lands where the vendor puts it.** Sampling the 15 h frame at five
elevations, XSLOPE's $h=7$ m contour stands **1.1 to 1.8 m** inside the upstream face,
against RS2's wetting ramp at **1.5 to 3.3 m** — the reservoir raise has barely entered
the fill in both. (The order-of-magnitude front position this row previously reported came
from a saturated conductivity that was 3600× too large; the vendor's $10^{-7}$ m/s is
already the m/s value, and converting it once to the m/hr the schedule runs in gives
$3.6\times10^{-4}$ m/hr.)

**The late frame is a timing difference, not a field difference.** XSLOPE's steady
total-head field reproduces Fig 19-5 exactly as described — reservoir head 10 drawn down
through the dam to the toe drain at total head 0, phreatic surface descending to the drain
— but it is not yet fully there at 16383 h, where RS2 essentially is. At 16383 h XSLOPE's
four station heads read 6.54 / 6.92 / 4.37 / 3.22 against a steady 7.20 / 7.52 / 5.84 /
4.40, and the field settles to within 0.01 m of steady by $\approx5\times10^{4}$ h. The two
stations still furthest from steady sit between the crest and the toe drain, where the drain
paces the last of the drawdown. [GW18](#gw18) measures the same approach against a
digitized profile rather than contours.

*Reading Fig 19-4: its colour ramp runs the opposite way to Fig 19-5's on the same page.
At the toe drain, held at total head 0, Fig 19-4 is red where Fig 19-5 is blue; just inside
the submerged upstream face, held at 10, Fig 19-4 is blue where Fig 19-5 is red. Read with
Fig 19-5's key, Fig 19-4 places the 15 h front at the wrong end of the dam.*

![gw017: steady total-head field vs Fig 19-5](images/gw017.png)

<!-- test: file=files/rocscience_gw/gw017.xlsx, type=tseep_head, target_size=1.5, time=15, max_head_change_frac=0.25, points=26:4:2.068;26:8:2.128;32:10:1.486;36:8:1.004, tolerance=0.15, benchmark=GW17-t15 -->
<!-- test: file=files/rocscience_gw/gw017.xlsx, type=tseep_head, target_size=1.5, time=16383, max_head_change_frac=0.25, points=26:4:6.543;26:8:6.917;32:10:4.365;36:8:3.223, tolerance=0.15, benchmark=GW17-t16383 -->
<!-- test: file=files/rocscience_gw/gw017.xlsx, type=tseep_head, target_size=1.5, time=200000, max_head_change_frac=0.25, points=26:4:7.199;26:8:7.517;32:10:5.838;36:8:4.402, tolerance=0.15, benchmark=GW17-tsteady -->

### GW18: Transient seepage through an earth fill dam {#gw18}

**Input files:** [gw018.xlsx](files/rocscience_gw/gw018.xlsx)

The Fredlund & Rahardjo (1993) 12 m earth-fill dam — the same body as the steady
[GW6](#gw6) family (base 52 m, 4 m crest, symmetric 2:1 faces) — with the reservoir
**quickly raised from 4 m to 10 m at $t=0$** and no drain. The phreatic surface rises
through the dam and daylights on the downstream (toe) slope, the free seepage face.
It adds a moving reservoir boundary and an unsaturated-conductivity fit on top of the storage
formulation the [consolidation problems](#gw15) lock.

Storage $S_s=\gamma_w m_v = 9.81\times0.002=0.0196\ \text{m}^{-1}$ (the vendor RS2 model
carries $m_v=0.002$; the manual *text* prints 0.003 — the model value is used, since it
is what produced RS2's own Fig 20.5 result). Saturated conductivity
$k_s=10^{-7}$ m/s, the plateau of the vendor's own conductivity table and of Fig 18.1's
`k (m/s)` axis, expressed as $3.6\times10^{-4}$ m/hr for the schedule.

**The two unsaturated curves are transcribed independently**, because the vendor ships them
that way and they live on completely different scales. `groundwater #020`'s user-defined
function carries a conductivity table — (0, 10⁻⁷), (5, 10⁻⁷), (100, 10⁻¹⁰), suction in kPa —
and, separately, a water-content table running (0, 0.70) to (**4000**, 0.40). That retention
curve spreads its 0.3 of water content over 408 m of head, a hundred times the scale of any
van Genuchten fit to the conductivity curve, so a single $(\alpha, n)$ cannot carry both.
XSLOPE's Gardner option keeps them apart: $k_r$ comes from the power law
($a=0.0404$, $n=4.15$, fitted to the conductivity table at 0.120 decades rms over 0–12 m of
suction, better than the 0.161 a Mualem–van Genuchten fit reaches), and the drainage
capacity from an independent linear band $S_y/|h_0|$ with $S_y=0.70-0.40=0.30$ and
$h_0=-4000/\gamma_w=-407.7$ m — which *is* the vendor's retention line, reproduced rather
than fitted.

The reservoir raise is a **submerged-only Dirichlet series** on the upstream face (the
$t=0$ steady solve at reservoir 4 sets the initial condition; every upstream node with
$y\le h(t)$ is then held at $h(t)$, nodes above becoming exit-face nodes). Linear tri3 mesh.

The published target is **total head sampled along the toe slope**, compared with Fig 20.5
at the vendor's own two stage times, $t=0.6$ h and $t=19656$ h. Fig 20.5 is a labeled
chart, so it is digitized at all eleven of its own $x$ stations; XSLOPE's steady profile is
solved and locked as a third frame.

| $x$ (m) | 28 | 30 | 32 | 34 | 36 | 38 | 40 | 42 | 44 | 46 | 48 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| Fig 20.5, $t=0.6$ h | 2.887 | 2.775 | 2.687 | 2.560 | 2.448 | 2.323 | 2.177 | 2.047 | 1.882 | 1.680 | 1.406 |
| XSLOPE, $t=0.6$ h | 2.840 | 2.750 | 2.606 | 2.492 | 2.378 | 2.217 | 2.092 | 1.968 | 1.785 | 1.617 | 1.386 |
| Fig 20.5, $t=19656$ h | 8.330 | 8.001 | 7.639 | 7.238 | 6.794 | 6.286 | 5.683 | 4.970 | 4.019 | 3.014 | 2.009 |
| XSLOPE, $t=19656$ h | 8.192 | 7.896 | 7.452 | 7.100 | 6.720 | 6.125 | 5.608 | 4.993 | 3.927 | 3.000 | 2.130 |
| XSLOPE, steady | 8.250 | 7.962 | 7.525 | 7.176 | 6.793 | 6.188 | 5.657 | 5.019 | 3.929 | 3.000 | 2.131 |

| frame | rms | worst |
|---|---|---|
| $t=0.6$ h vs Fig 20.5 at 0.6 h | 0.072 m | 0.106 m |
| $t=19656$ h vs Fig 20.5 at 19656 h | 0.115 m | **0.187 m** |
| XSLOPE steady vs Fig 20.5 at 19656 h | 0.074 m | 0.122 m |

**The two solutions agree on the steady profile and on the approach to it.** XSLOPE's
steady toe-slope profile lands within 0.12 m of Fig 20.5's 19656 h curve at every station,
and the 19656 h frame of the transient march itself is within 0.115 m rms / 0.187 m worst of
that curve — below the ≈0.2 m the chart can be read to, so at the figure's own resolution the
two codes are at the same place at the same time. By 19656 h XSLOPE has closed to 0.053 m rms
of its own steady profile (it is within 0.003 m of it by $\approx3\times10^{4}$ h), which is
the same "already steady" state Fig 20.5's 19656 h curve represents.

The **storage convention** is what sets that timing, here and on every transient row of this
page. XSLOPE applies the elastic specific storage $S_s$ in the saturated zone only and takes
the storage above the phreatic surface from the retention curve alone — the same convention
Slide and SEEP/W use. On this problem the two coefficients differ by a factor of about 27
($S_s=0.0196$ against a retention capacity $S_y/|h_0| = 7.4\times10^{-4}\ \text{m}^{-1}$),
because the vendor's retention line is nearly flat over the 0–12 m of suction the dam
reaches, so which of the two acts above the water table governs the drainage time-scale
outright. The 0.6 h frame is insensitive to it: that early, the raise has barely entered the
fill and the response is saturated.

![gw018: toe-slope total head, XSLOPE vs digitized Fig 20.5](images/gw018.png)

<!-- test: file=files/rocscience_gw/gw018.xlsx, type=tseep_head, target_size=1.5, time=0.6, max_head_change_frac=0.25, points=30:11:2.750;35:8.5:2.416;40:6:2.092;45:3.5:1.738, tolerance=0.15, benchmark=GW18-t0.6 -->
<!-- test: file=files/rocscience_gw/gw018.xlsx, type=tseep_head, target_size=1.5, time=19656, max_head_change_frac=0.25, points=30:11:7.896;35:8.5:6.848;40:6:5.608;45:3.5:3.654, tolerance=0.15, benchmark=GW18-t19656 -->
<!-- test: file=files/rocscience_gw/gw018.xlsx, type=tseep_head, target_size=1.5, time=60000, max_head_change_frac=0.25, points=30:11:7.962;35:8.5:6.923;40:6:5.657;45:3.5:3.656, tolerance=0.15, benchmark=GW18-tsteady -->

### GW19: Transient seepage below a lagoon {#gw19}

**Input files:** [gw019.xlsx](files/rocscience_gw/gw019.xlsx)

A lined lagoon leaking into a two-layer aquifer, modelled as a **half-model** (the lagoon
centerline at $x=0$ is a no-flow symmetry plane): 19 m wide × 10 m deep, a 1 m **soil
liner** across the top ($y\in[9,10]$, the lower-permeability layer) over 9 m of **soil**
($y\in[0,9]$). A 2 m-wide lagoon ($x\in[0,2]$ at the surface) is filled with 1 m of water
at $t=0$ and leaks down through the liner.

From the vendor RS2 model both materials carry $m_v=0.002$ and porosity 0.7; the soil
Custom $k(\psi)$ has $k_s=6\times10^{-4}$ m/min (down two decades by 10 m suction) and the
liner $k_s=3.54\times10^{-4}$ m/min with the same relative shape. Both of the vendor's
unsaturated curves are **reproduced rather than fitted**, the same split
[GW18](#gw18) uses: the conductivity table is three points, $k_r=1/0.1/0.01$ at
$0/5/10$ m of suction, which the Gardner law $k_r=1/(1+a\psi^n)$ passes through exactly at
$n=\ln 11/\ln 2=3.459$, $a=9/5^{n}=0.0344$; and the water-content table is the two points
$(0,0.7)\rightarrow(100\ \text{kPa},0.5)$, a straight line, which is exactly the linear
drainage band the Gardner option uses for moisture capacity — $S_y=0.2$ over
$h_0=-100/10=-10$ m. Storage $S_s=\gamma_w m_v=10\times0.002=0.02\ \text{m}^{-1}$. The far field (right edge,
$y\in[0,5]$) is held at total head 5 — the regional water table 5 m below the surface, which
is also the **initial condition**: the $t=0$ steady solve at a lagoon head of 5 (equal to the
far field) gives a uniform water table at el 5, then the lagoon steps to total head 11
(el 10 + 1 m ponded) for $t>0$. The base, the centerline and the top away from the lagoon are
no-flow. Report times are the vendor stage schedule: 73 / 416 / 792 / 11340 min.

The published target is **pressure head along the top boundary** (Fig 21.9, vs Ref [1]
Fredlund & Rahardjo) plus pressure-head contours at the four times. Fig 21.9 has no tabulated
companion, but it is a marker plot on labeled axes carrying two anchors that fix its
calibration independently of any reading: the far-field markers sit at −5.000 ± 0.005 m,
exactly the initial water table 5 m below the top boundary, and the lagoon markers at
+0.996, exactly the 1 m of ponded water. Digitized against those anchors it resolves to
≈0.005 m — an order of magnitude finer than the ≈0.2 m read-off precision Fig 20.5 supports
for [GW18](#gw18).

**The near-steady frame.** Against RS2's own markers at 20 stations along the top boundary,
XSLOPE's 11340 min profile agrees to **0.045 m rms and 0.126 m at worst** over a pressure-head
range of more than 5 m — under 1% of the driving head. Five top-boundary stations of that
frame are locked.

| $x$ (m) | 3 | 6 | 10 | 14 | 19 |
|---|---|---|---|---|---|
| XSLOPE pressure head | −0.37 | −1.83 | −2.73 | −3.54 | −4.15 |
| Fig 21.9 (RS2) | −0.48 (+0.10) | −1.79 (−0.05) | −2.76 (+0.03) | −3.57 (+0.04) | −4.18 (+0.02) |

**The early frames.** These track less closely than the near-steady one — 0.35 / 0.14 /
0.07 m rms at 73 / 416 / 792 min — and the residual has a consistent sign: XSLOPE's mound runs
slightly *ahead* of RS2's, by up to 0.61 m at 73 min (at *x* = 3 m, the station nearest the
lagoon) and by up to 0.27 and 0.19 m at 416 and 792 min. The lead closes monotonically as the
field fills.

The rate the mound fills at is set by the moisture capacity, which is why this model is built
on the Gardner split rather than on one Mualem–van Genuchten pair carrying both vendor curves.
A single $(\alpha, n)$ fitted to the conductivity table also has to supply the capacity, which
then follows the *conductivity* curve's shape instead of the vendor's retention line; on this
problem that costs a factor of two to four at every frame (0.62 / 0.47 / 0.27 / 0.089 m rms,
against 0.35 / 0.14 / 0.07 / 0.045 for the split). The split is available here for the same
reason it is on [GW18](#gw18) — the retention table is a straight line, so the band
$S_y/|h_0|$ *is* that line — with the added benefit that Gardner passes through this model's
three conductivity points exactly where the van Genuchten fit returned 0.095 and 0.0128 for
$k_r$ = 0.1 and 0.01. The residual lead is what remains once both curves are the vendor's own.

![gw019: pressure head along the top boundary as the lagoon fills](images/gw019.png)

<!-- test: file=files/rocscience_gw/gw019.xlsx, type=tseep_head, target_size=0.8, time=73, max_head_change_frac=0.25, points=1:8:6.088;3:8:5.417;1:5:5.033, tolerance=0.15, benchmark=GW19-t73 -->
<!-- test: file=files/rocscience_gw/gw019.xlsx, type=tseep_head, target_size=0.8, time=416, max_head_change_frac=0.25, points=1:8:7.706;3:8:6.960;1:5:5.924, tolerance=0.15, benchmark=GW19-t416 -->
<!-- test: file=files/rocscience_gw/gw019.xlsx, type=tseep_head, target_size=0.8, time=792, max_head_change_frac=0.25, points=1:8:8.113;3:8:7.416;1:5:6.464, tolerance=0.15, benchmark=GW19-t792 -->
<!-- test: file=files/rocscience_gw/gw019.xlsx, type=tseep_head, target_size=0.8, time=11340, max_head_change_frac=0.25, points=1:8:9.470;3:8:9.028;1:5:8.626, tolerance=0.15, benchmark=GW19-t11340 -->
<!-- test: file=files/rocscience_gw/gw019.xlsx, type=tseep_head, target_size=0.8, time=11340, max_head_change_frac=0.25, points=3:10:9.626;6:10:8.166;10:10:7.266;14:10:6.465;19:10:5.846, tolerance=0.15, benchmark=GW19-top-t11340 -->

### GW20: Transient seepage in a layered slope {#gw20}

**Input files:** [gw020.xlsx](files/rocscience_gw/gw020.xlsx)

The steady [GW7](#gw7) Rulon & Freeze sandbox (2.4 × 1.0 m, a medium-sand slope with a thin
fine-sand lens) re-run **transiently** as rainfall switches on at $t=0$. The geometry,
materials and boundary layout are GW7's exactly; the transient additions are storage
($S_s=\gamma_w m_v=9.81\times0.002=0.0196\ \text{m}^{-1}$, $S_y=0.7-0.5=0.2$ from the
vendor's water-content table) and a rainfall flux
that is **zero at $t=0$ and steps to the GW7 value $2.1\times10^{-4}$ m/s** for $t>0$.

The initial condition is the $t=0$ steady solve with no infiltration — only the tailwater
head 0.3 at the toe holds, so the water table sits at el 0.3 (0.1 m above the toe). When the
rainfall (above the fine sand's $k_s=5.5\times10^{-5}$ m/s) switches on, water perches on the
fine lens and the perched mound builds toward the steady GW7 result. Report times are the
vendor schedule: 4.6 / 31 / 208 s (the medium-sand diffusion time over 1 m is $\approx14$ s,
so the 208 s frame has essentially reached the steady perched state).

The published target is **total head along a query line**: Fig 22.6 draws it vertical at
$x=1.6$ m — the crest break — running the full 1 m from the base to the crest, and Fig 22.7
plots RS2 and Ref [1] against depth along it at the three times. That chart has no tabulated
companion but calibrates itself: its 4.6 s markers at the base of the line read 0.302–0.305
against the model's initial total head of 0.300, so the digitization is good to ≈0.005 m.

**RS2 is at steady state by 208 s; XSLOPE arrives around 400 s.** RS2's 208 s markers sit on
the steady field — XSLOPE's own steady solve of the identical slope ([GW7](#gw7)) reproduces
them to 0.005 m rms over 21 stations, which is where GW7's numeric lock now comes from.
XSLOPE's *transient* march is not quite there at 208 s: it is still 0.105 m rms below,
uniformly along the whole line (0.096 to 0.126 m at every station). It reaches the profile by
400 s and then holds it:

| frame | rms Δ vs the Fig 22.7 RS2 markers | worst |
|---|---|---|
| $t$ = 4.6 s | 0.037 m | 0.073 m |
| $t$ = 31 s | 0.155 m | 0.244 m |
| $t$ = 208 s | 0.105 m | 0.126 m |
| $t$ = 400 s | 0.004 m | 0.010 m |
| $t$ = 800 s and beyond | 0.005 m | 0.011 m |

So the end state is right and the approach to it is roughly twice as slow. Unlike
[GW18](#gw18), the residual here is not a storage-*convention* difference — the perching
mechanism runs almost entirely in the unsaturated zone, where the retention capacity governs
and the elastic $S_s$ never enters. It is the retention *shape*: both unsaturated curves come
from the vendor's own tables — the Mualem–van Genuchten pairs are fitted to its conductivity
tables and $S_y=0.7-0.5=0.2$ is the drop across its water-content table — but a single
$(\alpha, n)$ has to carry both, so the moisture capacity follows the shape of the
conductivity curve (which falls over ~1 m of suction here) rather than the vendor's own
retention line (which falls over 10 m). That can change the rate at which the perched mound
fills but not the state it fills to — which is exactly the shape of the measurement, 0.11 m
apart at 208 s and 0.005 m apart at steady state. XSLOPE's own solved heads at four interior
stations are locked at the three published report times as a regression guard on the
transient frames themselves; the vendor-facing lock for this geometry is GW7's.

![gw020: total head along the query line as rainfall perches on the lens](images/gw020.png)

<!-- test: file=files/rocscience_gw/gw020.xlsx, type=tseep_head, target_size=0.04, time=4.6, max_head_change_frac=0.25, points=2.2:0.95:0.357;2:0.85:0.305;2:0.75:0.300;1.6:0.72:0.300, tolerance=0.15, benchmark=GW20-t4.6 -->
<!-- test: file=files/rocscience_gw/gw020.xlsx, type=tseep_head, target_size=0.04, time=31, max_head_change_frac=0.25, points=2.2:0.95:0.490;2:0.85:0.412;2:0.75:0.375;1.6:0.72:0.337, tolerance=0.15, benchmark=GW20-t31 -->
<!-- test: file=files/rocscience_gw/gw020.xlsx, type=tseep_head, target_size=0.04, time=208, max_head_change_frac=0.25, points=2.2:0.95:0.809;2:0.85:0.780;2:0.75:0.768;1.6:0.72:0.715, tolerance=0.15, benchmark=GW20-t208 -->

### GW21: Transient flow in a fully confined aquifer {#gw21}

**Input files:** [gw021a.xlsx](files/rocscience_gw/gw021a.xlsx) (IC = 0) ·
[gw021b.xlsx](files/rocscience_gw/gw021b.xlsx) (IC = 5 ft)

A 100 ft × 5 ft fully confined, fully saturated aquifer (imperial units): $k=4$ ft/hr,
$m_v=0.1$, $\gamma_w=62.4$, giving $S_s=6.24\ \text{ft}^{-1}$ and diffusivity
$D=T/S=k/(\gamma_w m_v)=0.641\ \text{ft}^2/\text{hr}$. The head at the left face is stepped up
5 ft at $t=0$; the aquifer is long enough that the far end stays undisturbed at 600 hr (the
front reaches $\sqrt{4Dt}\approx39$ ft), so the head rise follows J. G. Ferris' semi-infinite
solution (Tao & Xi 2006),

$$ \Delta h(x,t)=\Delta H\,\operatorname{erfc}\!\left(\frac{x}{\sqrt{4Dt}}\right). $$

Case 1 starts from zero head; case 2 from a uniform 5 ft steady head and steps to 10 ft.
XSLOPE reproduces the erfc profile at 600 hr to within **0.015 ft** across the domain (below).
The tags lock the closed-form head at five stations, tolerance 0.05 ft.

![gw021: Ferris confined-aquifer profile, erfc vs XSLOPE](images/gw021.png)

<!-- test: file=files/rocscience_gw/gw021a.xlsx, type=tseep_head, target_size=0.8, time=600, points=10:2.5:103.592;20:2.5:102.354;30:2.5:101.397;40:2.5:100.746;50:2.5:100.357, tolerance=0.05, benchmark=GW21a -->
<!-- test: file=files/rocscience_gw/gw021b.xlsx, type=tseep_head, target_size=0.8, time=600, points=10:2.5:108.592;20:2.5:107.354;30:2.5:106.397;40:2.5:105.746;50:2.5:105.357, tolerance=0.05, benchmark=GW21b -->

## The SEEP2D cross-check: where does the free surface daylight? {#seep2d-crosscheck}

Several problems on this page reproduce a published profile in shape but sit a little
high, and what they share is a free surface on a steep exit face. The question is whether
XSLOPE places the **release point** — the top of the saturated seepage face, where the
phreatic surface daylights — differently from other codes.

That question is answerable directly, because XSLOPE's seepage solver is in the SEEP2D family. The
original USACE/WES SEEP2D Fortran program is run on the *identical* tri3 mesh, the
identical boundary conditions and the identical unsaturated law (SEEP2D's `iuntyp`:
linear front, or van Genuchten — both codes use the same Mualem–van Genuchten form, so
α and n carry straight across), and the two solutions are compared node for node. The
harness is `benchmarks/run_seep2d_compare.py --gw`.

| problem | law | release point, XSLOPE | release point, SEEP2D | head RMS | head range |
|---|---|---|---|---|---|
| gw004 | linear front | 0.500 | 0.500 (0.0) | 0.0004 | 0–4 |
| gw006a | van Genuchten | 0.000 | 0.000 (0.0) | 0.0000 | 0–10 |
| gw009a | van Genuchten | 8.846 | 8.462 (+0.384) | 0.0026 | 0–18.5 |
| gw010 | van Genuchten | 4.872 | 4.872 (0.000) | 0.0006 | 2–10 |
| gw012 | linear front | face dry | face dry | 0.034 | 0–50 |
| gw013 | linear front | face dry | face dry | 0.070 | 0–50 |

**The release point agrees.** It is identical on gw004, gw006a and gw010; both codes
leave the face fully drained on gw012 and gw013; and gw009a differs by 0.384 — a single
element at that mesh. The head fields agree to a relative RMS of order 10⁻⁴ throughout.
XSLOPE does not release low. Where this corpus sits above a published profile, the
difference is between the SEEP2D family and the reference code's exit-face convention,
not an XSLOPE defect, so those problems are locked on XSLOPE's own values with the
offset reported.

One difference remains between the two codes. On the **van Genuchten** problems the total
discharge reads 3.5–4.7% below SEEP2D (gw009a 2.299×10⁻⁵ vs 2.412×10⁻⁵; gw010
6.070×10⁻⁵ vs 6.294×10⁻⁵) even though the heads agree to 10⁻⁴ — while the linear-front
problems agree on discharge to better than 0.15%. The split follows the unsaturated law
exactly, and it is *not* XSLOPE's kr floor: dropping `kr_min` from 10⁻⁴ to zero leaves
the discharge unchanged to six figures. Since the head field is right and only the
integrated flux differs, the likely cause is where each code evaluates the strongly
nonlinear kr(ψ) when it forms an element's conductivity. It does not affect pore
pressures or stability, which read the head field.

### The specified-flux boundary, against SEEP2D {#flux-crosscheck}

The same harness verifies the
[specified-flux (Neumann) boundary](../seep/overview.md#specified-flux-boundary-conditions-neumann),
because SEEP2D supports one natively. It is entered as *flowrate cards* — a boundary segment
and a flux — and SEEP2D forms its own consistent nodal loads from them (`seep2d.f`, the
`nflcd` block: `fx(i) += ½·L·q` at each end of the segment, and it sets the flux flags itself).
So SEEP2D is handed the raw segment and flux, never XSLOPE's assembled vector, and the two
codes' assemblies are compared rather than one being fed the other's answer.

The test problem is a confined rectangle, 10 × 4, isotropic *k* = 2, with a fixed head on one
side and a uniform inflow *q* on the other, no-flow above and below. One-dimensional Darcy
gives the head exactly: *h*(*x*) = *q*(*L* − *x*)/*k*, with total inflow *q*·*H*.

| | XSLOPE | SEEP2D |
|---|---|---|
| max abs. error vs the exact solution | **2.5×10⁻¹⁹** | 5.0×10⁻⁹ |
| total inflow (exact: 1.760000×10⁻⁵) | **1.760000×10⁻⁵** | 1.759800×10⁻⁵ (+2.0×10⁻⁹) |
| max abs. difference, XSLOPE − SEEP2D | \-- | **5.0×10⁻⁹** (2.3×10⁻⁴ of the head range) |

*XSLOPE reproduces the closed form to machine precision, and matches SEEP2D to that code's own
iteration tolerance.*

Machine-precision agreement on the *head field* — not merely on the total — is the check that
matters. A linear head field only satisfies **A**·**h** = **f** when **f** is the *consistent*
load vector, so an incorrect distribution along an edge (½, ½, 0 across a quadratic edge, say)
would still sum to *q·L* and pass a total-flux check while breaking the solution. It does not.
The quadratic weights (⅙, ⅙, ⅔ at corner, corner and midside) are verified the same way, and a
zero flux reproduces the no-flux answer bitwise.

## Methodology

Same rules as the [Slide2 corpus](rocscience.md): problems are built from the manual's tabulated
data and coordinate-labeled figures; where a figure is unlabeled, geometry is extracted by
axis-calibrated pixel measurement and validated against printed solution quantities; every built
problem is locked into `run_tests.py` via test tags. Seepage problems compare flow rates,
phreatic-surface positions, and head/pressure profiles rather than factors of safety. Three tag
types carry the locks: `type=seep` (flowrate, live mesh + solve at `target_size`),
`type=seep_head` (solved total head at named `x:y:h` points, interpolated from the four nearest
nodes), and `type=tseep_head` — the transient sibling of `seep_head`, which meshes and samples
head the identical way but pulls the field from the frame of a transient solve at a given save
time `t` (`time=…`, an entry the solver lands on exactly). A `type=tseep_head` tag names a file
carrying a v18 `tseep` sheet; optional `dt_max` / `max_head_change_frac` / `theta` tune the
stepper.

Where a problem's published answer is itself a chart, how it is handled depends on what the
chart resolves. A contour plate with a numeric key gives a band per point, and the comparison
is reported as the fraction of sampled points falling inside the vendor's own band (GW5's
Fig 5-4); a plate without one supports only a qualitative reading (GW17's Fig 19-4). A line or
marker plot on labeled axes is digitized and compared numerically — and where the model itself
fixes a value the chart must show (an initial condition, a boundary head), that value calibrates
the digitization, so its precision is a measured quantity rather than an estimate (GW7's
Fig 22.7, GW18's Fig 20.5, GW19's Fig 21.9). In every chart-target case the tag locks XSLOPE's
own solved field; the vendor comparison is reported in the row text with its precision.
