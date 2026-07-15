# Rocscience RS2 (SSRM) Corpus

This page tracks the RS2 Slope Stability Verification Manual (Rocscience, Parts I–III,
68 problems) the way the [Slide2 corpus](rocscience.md) tracks its manual — but for the
**shear strength reduction (SSRM)** method against XSLOPE's FEM/SSRM solver rather than
limit equilibrium. The long-standing SSRM anchors (Griffiths & Lane 1999 and the feature
samples) live on the [SSRM benchmarks page](ssrm.md).

Where a problem shares its geometry with a built Slide2 problem, the SSRM analysis runs on
the **same corpus input file** — the extraction is already validated there. SSRM results
use the Griffiths elastic convention (E = 10⁵ kPa, or its psf equivalent on the imperial
problems; ν = 0.3, ψ = 0). SSRM factors are insensitive to these, so the corpus builder
fills them into any material that does not publish its own, and the LEM problems carry
them as inert values. SSRM factors are quoted at the tagged mesh size: SSRM drifts a
percent or two with refinement, so tolerances are honest rather than tight. A recurring
pattern on this corpus: fine-mesh SSRM finds shallow-skin mechanisms that published
coarse-mesh SSRM analyses miss or deliberately suppress with a "can't fail" elastic region
(#23) — where the published value depends on such an artifice rather than the mechanics,
the problem is recorded as not lockable rather than tuned to match.

The same theme sets how far a *mesh* can be trusted. Where the failure mechanism is
pinned by geometry — a weak seam, a bedrock contact — the SSRM factor barely moves with
refinement (#18 returns the same value at two mesh sizes). Where nothing pins it, the
shear band is free to keep localizing as the elements shrink, because Mohr-Coulomb
without a regularizing length scale has nothing to stop it, and the factor drifts without
reaching a plateau (#14, under r<sub>u</sub> = 0.5). Such a problem is reported with its
whole mesh sweep and locked at a pinned mesh as a *regression* test, not advertised as a
converged value.

The RS2 manual is unusually cheap to build against: a large fraction of its problems are
**SSRM renditions of the same problems as the Slide2 LEM manual**, so the geometry and
materials are already extracted, validated, and sitting in the corpus input files — often
the only new work per problem is an SSRM run and a tag. Problems 56–58 additionally carry
published FS values from **Z-Soil, PLAXIS, and GEO FEM**, giving multi-program SSRM
cross-bearings.

## Methodology

Same discipline as the [Slide2 corpus](rocscience.md): geometry from the manuals'
coordinate-labeled figures (or reused directly from the Slide2 corpus input files where the
problem is shared), results locked into `run_tests.py` via `fem_ssrm` test tags — the tag
type and runner already exist and currently guard the Griffiths & Lane anchors. SSRM runs
are expensive (~1 min each), so this corpus will lean on coarse meshes with honest
tolerances, the same trade documented for the FEM reliability regression.

Each figure in the problem details below has two panels: the **left** panel is the FEM
model (elements, materials, boundary conditions) and the **right** panel is the maximum
shear strain contours at the critical SRF.

## Status

### Part I (1–34)

| # | Problem | Status | XSLOPE file / results |
|---:|---|---|---|
| [1](#rs2-1) | Simple slope stability assessment | **built** | [xslope_acads_simple.xlsx](../lem/files/xslope_acads_simple.xlsx). SSRM 0.986 vs RS2 SSRM 0.99, Slide Bishop 0.987, ACADS referee 1.00. |
| [2](#rs2-2) | Non-homogeneous slope | **built** | [vp003.xlsx](../files/rocscience/vp003.xlsx). SSRM 1.375 vs RS2 1.36, Slide Spencer 1.375, referee 1.39. |
| [3](#rs2-3) | Non-homogeneous slope with seismic load (0.15g) | **built** | [vp004.xlsx](../files/rocscience/vp004.xlsx). SSRM 0.948 vs RS2 0.97, Slide Spencer 0.991, referee 1.00. |
| [4](#rs2-4) | Dry Talbingo dam | **built** | [vp005.xlsx](../files/rocscience/vp005.xlsx). SSRM finds the true global minimum — the steeper downstream-bench skin (tan45/tan30.9 = 1.669) — at 1.678; the published RS2 1.88 / Slide 1.948 / referee 1.95 are the gentler upstream face. |
| [5](#rs2-5) | Water table with weak seam | **built** | [xslope_acads_weak_layer.xlsx](../lem/files/xslope_acads_weak_layer.xlsx). SSRM 1.286 (mesh-stable: 1.280 at 1.2 m) vs RS2 1.26, Slide Spencer 1.258, referee 1.24–1.27. |
| [6](#rs2-6) | Slope with load and pore pressure by water table (ACADS 4) | **built** (caveat) | [vp009.xlsx](../files/rocscience/vp009.xlsx). SSRM 0.81, mesh-stable, vs ACADS survey mean 0.808 and referee 0.78 — but +18% above RS2's SSRM 0.69 and Slide2's MC-optimized LEM 0.68–0.71. |
| [7](#rs2-7) | Pore pressure by digitized total head grid (ACADS 5) | **built** | [vp010.xlsx](../files/rocscience/vp010.xlsx). SSRM 1.483 vs RS2 SSRM 1.48 (+0.2%), on the FE-seepage model XSLOPE built for Slide2 VP10. Slide2 LEM 1.498–1.501, Giam 1.53. |
| [8](#rs2-8) | Saint-Alban test embankment | *no lock possible* | The grid encodes measured construction-induced pressures (see the Slide2 corpus VP11 row); RS2 SSRM 0.96 vs Pilot 1.04 recorded. |
| [9](#rs2-9) | Cubzac-les-Ponts test embankment | *no lock possible* | Measured pore-pressure grid plus a "can't fail" elastic face layer; RS2 SSRM 1.31 vs Pilot 1.24 recorded. |
| [10](#rs2-10) | Simple slope II (Arai & Tagyo ex. 1) | **built** | [xslope_arai_tagyo.xlsx](../lem/files/xslope_arai_tagyo.xlsx). SSRM 1.43 vs RS2 SSRM 1.40 (+2%), mesh-converged; LEM locks Bishop 1.404 / Spencer 1.401. |
| [11](#rs2-11) | Layered slope (Arai & Tagyo ex. 2) | **built** | [vp015.xlsx](../files/rocscience/vp015.xlsx). SSRM 0.41 vs RS2 SSRM 0.39 and Greco/Kim pattern-search 0.39–0.43; LEM locks 0.419–0.422. |
| [12](#rs2-12) | Simple slope + water table (Arai & Tagyo ex. 3) | **built** | [vp016.xlsx](../files/rocscience/vp016.xlsx). SSRM 1.12 vs RS2 SSRM 1.09 (+2.3%); LEM locks Bishop 1.112 / Spencer 1.113. |
| [13](#rs2-13) | Simple slope III (Yamagami & Ueta) | **built** | [vp017.xlsx](../files/rocscience/vp017.xlsx). SSRM 1.38 vs RS2 SSRM 1.33 and Greco Spencer 1.33; LEM locks Bishop 1.342 / Spencer 1.340. |
| [14](#rs2-14) | Simple slope, pore pressure by r<sub>u</sub> | **built** (caveat) | [vp018.xlsx](../files/rocscience/vp018.xlsx). The SSRM factor does not become mesh-independent on this model; the tag pins 2.0 m (0.948) as a regression lock, against RS2 SSRM 0.98, Slide2 Spencer 1.01 and Baker 1.02. |
| [15](#rs2-15) | Layered slope II (Greco ex. 4 / Yamagami & Ueta) | **built** | [vp019.xlsx](../files/rocscience/vp019.xlsx). SSRM 1.39 vs RS2 SSRM 1.39, Slide2 Spencer 1.398, Greco 1.40–1.42; mesh-converged. |
| [16](#rs2-16) | Layered slope and water table with weak seam (Greco ex. 5 / Chen & Shao) | **built** | [vp020.xlsx](../files/rocscience/vp020.xlsx). SSRM 0.972 vs RS2 SSRM 1.02, Slide2 Spencer 1.093 circular / 1.007 noncircular, Greco 0.973–1.1; LEM locks 1.086–1.091. |
| [17](#rs2-17) | Slope with three pore pressure conditions (Fredlund & Krahn) | **built** (dry + r<sub>u</sub>) | [vp021a](../files/rocscience/vp021a.xlsx) / [vp021b](../files/rocscience/vp021b.xlsx). Dry SSRM 2.02 vs RS2 SSRM 2.0, Slide2 M-P 2.075, F&K 2.076; r<sub>u</sub> = 0.25 SSRM 1.704 against F&K 1.761–1.766. |
| [18](#rs2-18) | Three pore pressure conditions and a weak seam (Fredlund & Krahn) | **built** (dry + r<sub>u</sub>) | [vp022a](../files/rocscience/vp022a.xlsx) / [vp022b](../files/rocscience/vp022b.xlsx). Dry SSRM 1.34 vs RS2 SSRM 1.34, Slide2 Bishop 1.382; r<sub>u</sub> = 0.25 SSRM 1.054 against Slide2 1.124 and F&K 1.124. |
| [19](#rs2-19) | Undrained layered slope (Low 1989) | **built** (caveat) | [vp024.xlsx](../files/rocscience/vp024.xlsx). SSRM 1.51 at the tagged mesh (1.48 at 0.6 m) vs RS2 SSRM 1.41, Slide2 LEM 1.439, Low 1.44. |
| [20](#rs2-20) | Slope with vertical load (Prandtl's wedge) | **built** | [vp025.xlsx](../files/rocscience/vp025.xlsx). SSRM 1.00 vs Prandtl theory 1.0 and RS2 SSRM 1.0; Slide2 Spencer 1.051 on the specified surface. |
| [21](#rs2-21) | Bearing capacity test prism (Prandtl II) | **built** | [vp026.xlsx](../files/rocscience/vp026.xlsx). SSRM 1.011 (1.020 at 1.2 m), converging on theory 1.0; RS2 SSRM 1.01; Slide2 Spencer 0.941 on the specified surface. |
| [22](#rs2-22) | Layered slope with undulating bedrock | **built** (SSRM variant) | [vp027_fem.xlsx](../files/rocscience/vp027_fem.xlsx). SSRM 1.558 (1.536 at 1.8 m, converging) vs RS2 SSRM **1.52** — +2.5% at the tagged mesh, +1.1% at the finer one. |
| [23](#rs2-23) | Underwater slope with linearly varying cohesion | *no lock possible* | RS2's published SSRM (1.12) depends on a "can't fail" elastic region; this slope's anchor remains the LEM lock ([VP29](rocscience.md#vp29), Spencer 1.145 on Duncan's surface). |
| [24](#rs2-24) | Layered slope with geosynthetic reinforcement | **built** | [vp032a](../files/rocscience/vp032a.xlsx) / [vp032c](../files/rocscience/vp032c.xlsx). SSRM 1.196 (H=7) and 0.991 (H=8.75) vs RS2 SSRM 1.15 / 0.95 (+4%, within the corpus SSRM scatter). |
| [25](#rs2-25) | Syncrude tailings dyke (El-Ramly et al. 2003) | **built** (caveat) | [vp033.xlsx](../files/rocscience/vp033.xlsx). SSRM 1.20 (mesh-converged) vs RS2 SSRM 1.29, Slide2 Bishop 1.305, El-Ramly 1.31. |
| [26](#rs2-26) | Clarence Cannon dam (Wolff & Harr 1987) | **built** | [vp034.xlsx](../files/rocscience/vp034.xlsx). SSRM 2.28 vs RS2 SSRM 2.29 (−0.4%); Slide2 GLE 2.333 / Spencer 2.383, W&H 2.36, XSLOPE LEM M-P 2.384. |
| [27](#rs2-27) | Homogeneous slope, pore pressure by r<sub>u</sub> | *blocked* | The FEM has no r<sub>u</sub> option (task, with RS2 #14); RS2's own text cites Slide2 [VP36](rocscience.md#vp36) (Li & Lumb), not VP21. |
| 28 | FE analysis with groundwater and stress | *planned* | Slide2 VP38 family. |
| [29](#rs2-29) | Geosynthetic-reinforced embankment on soft soil (Tandjiria 2002) | **built** (sand case) | [vp039c.xlsx](../files/rocscience/vp039c.xlsx). SSRM 1.253 vs RS2 SSRM 1.25, Spencer 1.209, Tandjiria 1.219. The clay case is not locked. |
| [30](#rs2-30) | Homogeneous slope, power-curve strength (Perry 1993) | **built** | [vp040.xlsx](../files/rocscience/vp040.xlsx). SSRM 0.910 vs RS2 SRF 0.91 (exact); Slide2 Janbu 0.944, Perry 0.98. |
| [31](#rs2-31) | M-C vs power curve (Baker 2003 ex. 1) | **built** (all three halves) | [vp044a](../files/rocscience/vp044a.xlsx) / [b](../files/rocscience/vp044b.xlsx) / [c](../files/rocscience/vp044c.xlsx). M-C SSRM 1.546 / 0.966 vs RS2 1.53 / 0.98; power-curve SSRM 0.934. |
| [32](#rs2-32) | Heading mismatch — body is Baker's example 2 | **built** (both halves) | [vp045a](../files/rocscience/vp045a.xlsx) / [vp045b](../files/rocscience/vp045b.xlsx). M-C SSRM 2.852 vs RS2 2.83 (+0.8%); power-curve SSRM 2.732 vs RS2 2.74 (−0.3%), Slide2 Spencer 2.662. |
| [33](#rs2-33) | Homogeneous slope with tension crack and water table (P&D test slope 2) | **built** (caveat) | [vp056.xlsx](../files/rocscience/vp056.xlsx). SSRM 1.247 vs RS2 SSRM 1.28 and an eight-program LEM table spanning 1.03–1.32. |
| [34](#rs2-34) | M-C vs power curve III (Baker 2003 ex. 3, London clay) | **built** (both halves) | [vp061a](../files/rocscience/vp061a.xlsx) / [vp061b](../files/rocscience/vp061b.xlsx). M-C SSRM 1.369 vs RS2 1.38; power-curve SSRM 1.502 vs RS2 1.47 / Slide2 Spencer 1.47 / Baker 1.48 (+2.1%). |

### Part II (35–58)

| # | Problem | Status | XSLOPE file / results |
|---:|---|---|---|
| 35 | Submerged slope | *planned* | Slide2 [VP64](rocscience.md#vp64) family. |
| [36](#rs2-36) | Seepage analysis, homogeneous slope (D&W Fig 6.37) | **built** (both cases) | [vp071a](../files/rocscience/vp071a.xlsx) / [b](../files/rocscience/vp071b.xlsx). SSRM 1.118 on the FE-seepage model and 1.118 on the piezo approximation vs RS2 SSRM 1.12 / 1.12; referee 1.138/1.141; XSLOPE LEM locks 1.132. |
| [37](#rs2-37) | Embankment with layered foundation (D&W Fig 6.39) | *reported, no lock* | RS2's SSRM is the artesian downstream-toe slide (0.95 in its table, 1.1 in its own convergence graph); XSLOPE's SSRM finds the deep mechanism at 1.31. |
| 38 | Cohesionless embankment on saturated clay foundation | *planned* | Slide2 [VP73](rocscience.md#vp73)/[VP74](rocscience.md#vp74) family. |
| 39, 41, 43 | Earth embankment, infinite-slope mechanism (I–III) | *planned* | Slide2 [VP69](rocscience.md#vp69) family. |
| [40](#rs2-40) | Dam with impermeable foundation (D&W Fig 7.24) | **built** (piezo case) | [vp077b.xlsx](../files/rocscience/vp077b.xlsx). SSRM finds the saturated-toe skin at 1.126 (true global minimum, ~5% below the idealized toe infinite slope); RS2 SSRM 1.53 reports a deeper face. FE-seepage case blocked. |
| [42](#rs2-42) | James dike | **built** | [vp075.xlsx](../files/rocscience/vp075.xlsx). SSRM 1.249 vs RS2 SSRM 1.26 (−0.9%); Slide2 noncircular LEM 1.11–1.16, referee 1.17. |
| [44](#rs2-44) | Seepage analysis for an earth embankment (D&W Fig 14.20-a) | **built** | [vp082.xlsx](../files/rocscience/vp082.xlsx). SSRM 1.511 vs RS2 SSRM 1.51 (+0.1%); Slide2 LEM 1.532/1.541, referee 1.528–1.542. |
| [45](#rs2-45) | Varying undrained shear strength profiles (D&W Fig 14.20-b) | **built** (caveat) | [vp083a](../files/rocscience/vp083a.xlsx) / [b](../files/rocscience/vp083b.xlsx). SSRM 1.41 / 1.37 vs RS2 SSRM 1.32 / 1.32 (D&W referee 1.28–1.33). |
| [46](#rs2-46) | Varying undrained strength profiles II (D&W Fig 15.9, c<sub>u</sub> = 300 + c<sub>z</sub>·z) | **built** | [vp084a–d](../files/rocscience/vp084a.xlsx). SSRM 0.80 / 0.95 / 1.08 / 1.19 vs RS2 SSRM 0.78 / 0.93 / 1.05 / 1.15 (+2–3%); D&W 0.75 / 0.90 / 1.03 / 1.13. |
| [47](#rs2-47) | Purely cohesive slope, varying thickness (D&W Fig 14.3) | **built** (30-ft case) | [vp078.xlsx](../files/rocscience/vp078.xlsx). SSRM 1.08 vs RS2 SSRM 1.03; D&W referee 1.124–1.135. The 46.5- and 60-ft variants are deferred. |
| [48–54](#rs2-48) | Multi-tiered geotextile walls (Leshchinsky & Han 2004) | *blocked* | Slide2 [VP87](rocscience.md#vp87)–VP93. The FEM reinforcement tensile-capacity cap is not enforced in SSRM. Extraction and published tables are complete for all seven; the series unlocks together when the cap is fixed. |
| 55 | Five 1.8-m tiers | *blocked* | Slide2 VP94; belongs to the #48–54 family — blocked with the series. |
| [56](#rs2-56) | Homogeneous slope vs Z-Soil, PLAXIS, GEO FEM (Pruska 2003, H = 7 m, 5 cases) | **built** | [rs2_56a](../files/rocscience/rs2_56a.xlsx) / [b](../files/rocscience/rs2_56b.xlsx). All five within ±3.3% of RS2's M-C and inside the four-program band; locks bracket the family (0.667 / 2.131). Full tables in [the Pruska section](#pruska). |
| [57](#rs2-57) | Pruska H = 10.5 m, 6 cases | **built** | [rs2_57a](../files/rocscience/rs2_57a.xlsx) / [b](../files/rocscience/rs2_57b.xlsx). All six within ±3.6% of RS2's M-C; locks 0.449 / 1.411. Full tables in [the Pruska section](#pruska). |
| [58](#rs2-58) | Pruska H = 14 m, 6 cases | **built** (5 of 6) | [rs2_58a](../files/rocscience/rs2_58a.xlsx) / [b](../files/rocscience/rs2_58b.xlsx). Four within ±3.6%; case 5 reads 0.667 vs a published 0.72–0.75 cluster and is unlocked pending explanation; locks 0.342 / 1.057. |

### Part III (59–68)

| # | Problem | Status | XSLOPE file / results |
|---:|---|---|---|
| 59 | Three-layered soil slope | *planned* | Görög & Török (2007), vs Slide2 + PLAXIS. |
| [60](#rs2-60) | Generalized Hoek–Brown, homogeneous slope | **built** (LEM) | [rs2_60a.xlsx](../files/rocscience/rs2_60a.xlsx) / [b](../files/rocscience/rs2_60b.xlsx) / [c](../files/rocscience/rs2_60c.xlsx). Three slope angles from Li, Merifield & Lyamin (2008) at GSI = 70, the strong-rock end of the criterion. Bishop/Spencer 1.009 / 1.017 / 1.008 against Li's F = 1.0. SSRM is not locked on this problem. |
| 61 | Local and global minima, homogeneous slope | *planned* | |
| 62 | Three-layered slope with a soft band | *planned* | |
| 63 | Homogeneous slope assessment | *planned* | |
| 64 | Three homogeneous landslides | *planned* | |
| 65 | Tailings dam | *planned* | |
| 66 | Embankment basal stability | *planned* | |
| 67 | Earth dam under steady & transient unsaturated seepage | *blocked* | Transient — blocked on a transient solver. |
| 68 | Seismically loaded slopes | *planned* | |

---

## Problem details

### RS2-1: Simple slope stability assessment {#rs2-1}

Slide2 counterpart: [VP1](rocscience.md#vp1) (ACADS 1a).

**Input files:** [xslope_acads_simple.xlsx](../lem/files/xslope_acads_simple.xlsx)

| Method | XSLOPE | Published |
|---|---|---|
| SSRM | 0.986 | RS2 SSRM 0.99 |

*Cross-bearings: Slide2 Bishop 0.987 (LEM); ACADS referee 1.00.*

FS reads 0.967 at half the element size — SSRM values are quoted at the tagged mesh.

<!-- test: file=../lem/files/xslope_acads_simple.xlsx, type=fem_ssrm, expected_fs=0.986, element_type=tri6, target_size=0.9, tolerance=0.01, f_min=0.7, f_max=1.3, max_iter=4000, benchmark=RS2-1 -->

![RS2-1: FEM model (left) and maximum shear strain contours at the critical SRF (right)](images/RS2-1.png)

### RS2-2: Non-homogeneous slope {#rs2-2}

Slide2 counterpart: [VP3](rocscience.md#vp3).

**Input files:** [vp003.xlsx](../files/rocscience/vp003.xlsx)

| Method | XSLOPE | Published |
|---|---|---|
| SSRM | 1.375 | RS2 SSRM 1.36 |

*Cross-bearings: Slide2 Spencer 1.375 (LEM); ACADS referee 1.39.*

<!-- test: file=../files/rocscience/vp003.xlsx, type=fem_ssrm, expected_fs=1.375, element_type=tri6, target_size=0.9, tolerance=0.01, f_min=1.0, f_max=1.7, max_iter=4000, benchmark=RS2-2 -->

![RS2-2: FEM model (left) and maximum shear strain contours at the critical SRF (right)](images/RS2-2.png)

### RS2-3: Non-homogeneous slope with seismic load (0.15g) {#rs2-3}

Slide2 counterpart: [VP4](rocscience.md#vp4).

**Input files:** [vp004.xlsx](../files/rocscience/vp004.xlsx)

| Method | XSLOPE | Published |
|---|---|---|
| SSRM | 0.948 | RS2 SSRM 0.97 |

*Cross-bearings: Slide2 Spencer 0.991 (LEM); ACADS referee 1.00.*

k is entered negative per the FEM sign convention — this is a left-facing slope, so the
pseudo-static force acts in −x, while the LEM takes the magnitude and directs it from the
failure surface.

<!-- test: file=../files/rocscience/vp004.xlsx, type=fem_ssrm, expected_fs=0.948, element_type=tri6, target_size=0.9, tolerance=0.01, f_min=0.7, f_max=1.3, max_iter=4000, benchmark=RS2-3 -->

![RS2-3: FEM model (left) and maximum shear strain contours at the critical SRF (right)](images/RS2-3.png)

### RS2-4: Dry Talbingo dam {#rs2-4}

Slide2 counterpart: [VP5](rocscience.md#vp5).

**Input files:** [vp005.xlsx](../files/rocscience/vp005.xlsx)

| Method | XSLOPE | Published |
|---|---|---|
| SSRM | 1.678 | RS2 SSRM 1.88 (upstream face) |

For a dry cohesionless dam the critical mechanism is a surface-parallel (infinite-slope)
slide, FS = tan φ / tan β, which is *independent of depth* — so the steepest face governs. The
per-node SSRM criterion finds the true global minimum on the steeper **downstream** bench
(30.9°): tan 45° / tan 30.9° = 1.669, and the FEM returns 1.678 (+0.5%). The published values
report the gentler, constrained **upstream** face (27.2°, the end-of-construction problem):
tan 45° / tan 27.2° = 1.948, matched by RS2 SSRM 1.88, Slide2 1.948 (all LEM methods collapse
to tan φ / tan β on a cohesionless face), and the ACADS referee 1.95. Both are correct
infinite-slope answers; XSLOPE reports the more critical one. The seeded LEM search
([VP5](rocscience.md#vp5)) stays on the upstream circle in the input file and locks 1.955, so
the LEM and SSRM entries for this dam report different faces by construction, not a discrepancy.

*Closed-form check: across φ = 35–45° (c = 0 materials only) the SSRM tracks tan φ / tan 30.9° to 0.3%.*

<!-- test: file=../files/rocscience/vp005.xlsx, type=fem_ssrm, expected_fs=1.678, element_type=tri6, target_size=6.5, tolerance=0.01, f_min=1.5, f_max=2.3, max_iter=16000, benchmark=RS2-4 -->

![RS2-4: FEM model (left) and maximum shear strain contours at the critical SRF (right)](images/RS2-4.png)

### RS2-5: Water table with weak seam {#rs2-5}

Slide2 counterpart: **VP7** (inventory-only on the LEM page — no detail section to link).

**Input files:** [xslope_acads_weak_layer.xlsx](../lem/files/xslope_acads_weak_layer.xlsx)

| Method | XSLOPE | Published |
|---|---|---|
| SSRM | 1.286 | RS2 SSRM 1.26 |

*Mesh-stable: 1.280 at 1.2 m. Cross-bearings: Slide2 Spencer 1.258 (LEM); ACADS referee 1.24–1.27.*

<!-- test: file=../lem/files/xslope_acads_weak_layer.xlsx, type=fem_ssrm, expected_fs=1.286, element_type=tri6, target_size=2.0, tolerance=0.01, f_min=0.9, f_max=1.6, max_iter=4000, benchmark=RS2-5 -->

![RS2-5: FEM model (left) and maximum shear strain contours at the critical SRF (right)](images/RS2-5.png)

### RS2-6: Slope with load and pore pressure by water table (ACADS 4) {#rs2-6}

Slide2 counterpart: [VP9](rocscience.md#vp9). Built with a caveat.

**Input files:** [vp009.xlsx](../files/rocscience/vp009.xlsx)

| Method | XSLOPE | Published |
|---|---|---|
| SSRM | 0.81 (mesh-stable) | RS2 SSRM 0.69 |

*Cross-bearings: ACADS survey mean 0.808 and referee 0.78; Slide2's MC-optimized LEM 0.68–0.71; XSLOPE's own LEM locks 0.724.*

XSLOPE's SSRM lands on the ACADS survey mean but sits +18% above RS2's SSRM and Slide2's LEM.
The published values themselves span 0.68–0.81 on this thin-weak-seam problem; under the same
investigation as [#16](#rs2-16).

<!-- test: file=../files/rocscience/vp009.xlsx, type=fem_ssrm, expected_fs=0.812, element_type=tri6, target_size=1.3, tolerance=0.02, f_min=0.3, f_max=1.3, max_iter=4000, benchmark=RS2-6 -->

![RS2-6: FEM model (left) and maximum shear strain contours at the critical SRF (right)](images/RS2-6.png)

### RS2-7: Pore pressure by digitized total head grid (ACADS 5) {#rs2-7}

Slide2 counterpart: [VP10](rocscience.md#vp10).

**Input files:** [vp010.xlsx](../files/rocscience/vp010.xlsx)

| Method | XSLOPE | Published |
|---|---|---|
| SSRM | 1.483 | RS2 SSRM 1.48 (+0.2%) |

*Cross-bearings: Slide2 LEM 1.498–1.501; Giam 1.53.*

The SSRM runs on the FE-seepage model XSLOPE built for Slide2 [VP10](rocscience.md#vp10) (the
grid is a stand-in for the flow solution; sidecars are tri6 so the SSRM plasticity is not
volumetrically locked).

<!-- test: file=../files/rocscience/vp010.xlsx, type=fem_ssrm, expected_fs=1.483, tolerance=0.01, f_min=1.0, f_max=2.2, max_iter=4000, benchmark=RS2-7 -->

![RS2-7: FEM model (left) and maximum shear strain contours at the critical SRF (right)](images/RS2-7.png)

### RS2-8: Saint-Alban test embankment {#rs2-8}

Slide2 counterpart: **VP11** (inventory-only on the LEM page — no detail section to link).

| Method | XSLOPE | Published |
|---|---|---|
| SSRM | *no lock possible* | RS2 SSRM 0.96 |

*Cross-bearings: Pilot 1.04 recorded.*

The grid encodes measured construction-induced pressures (see the Slide2 corpus VP11 row), so
there is nothing here XSLOPE can reproduce as a lock.

### RS2-9: Cubzac-les-Ponts test embankment {#rs2-9}

Slide2 counterpart: **VP13** (inventory-only on the LEM page — no detail section to link).

| Method | XSLOPE | Published |
|---|---|---|
| SSRM | *no lock possible* | RS2 SSRM 1.31 |

*Cross-bearings: Pilot 1.24 recorded.*

Measured pore-pressure grid plus a "can't fail" elastic face layer suppressing the true face
failure (FS 1.11 per RS2's own text).

### RS2-10: Simple slope II (Arai & Tagyo ex. 1) {#rs2-10}

Slide2 counterpart: [VP14](rocscience.md#vp14) (Arai & Tagyo 1).

**Input files:** [xslope_arai_tagyo.xlsx](../lem/files/xslope_arai_tagyo.xlsx)

| Method | XSLOPE | Published |
|---|---|---|
| SSRM | 1.43 | RS2 SSRM 1.40 (+2%) |

*Mesh-converged (1.428→1.434 over a 2.9× size change). Cross-bearings: XSLOPE LEM locks Bishop 1.404 / Spencer 1.401 vs Slide2 1.409 / 1.406.*

<!-- test: file=../lem/files/xslope_arai_tagyo.xlsx, type=fem_ssrm, expected_fs=1.427, element_type=tri6, target_size=2.2, tolerance=0.02, f_min=1.2, f_max=1.7, max_iter=4000, benchmark=RS2-10 -->

![RS2-10: FEM model (left) and maximum shear strain contours at the critical SRF (right)](images/RS2-10.png)

### RS2-11: Layered slope (Arai & Tagyo ex. 2) {#rs2-11}

Slide2 counterpart: [VP15](rocscience.md#vp15).

**Input files:** [vp015.xlsx](../files/rocscience/vp015.xlsx)

| Method | XSLOPE | Published |
|---|---|---|
| SSRM | 0.41 | RS2 SSRM 0.39 |

*Cross-bearings: Greco/Kim pattern-search 0.39–0.43; XSLOPE LEM locks 0.419–0.422.*

<!-- test: file=../files/rocscience/vp015.xlsx, type=fem_ssrm, expected_fs=0.413, element_type=tri6, target_size=1.9, tolerance=0.02, f_min=0.25, f_max=0.65, max_iter=4000, benchmark=RS2-11 -->

![RS2-11: FEM model (left) and maximum shear strain contours at the critical SRF (right)](images/RS2-11.png)

### RS2-12: Simple slope + water table (Arai & Tagyo ex. 3) {#rs2-12}

Slide2 counterpart: [VP16](rocscience.md#vp16).

**Input files:** [vp016.xlsx](../files/rocscience/vp016.xlsx)

| Method | XSLOPE | Published |
|---|---|---|
| SSRM | 1.12 | RS2 SSRM 1.09 (+2.3%) |

*Cross-bearings: XSLOPE LEM locks Bishop 1.112 / Spencer 1.113.*

The FEM piezo pore pressure uses the vertical-distance convention, consistent with the LEM
slicer and the published analyses.

<!-- test: file=../files/rocscience/vp016.xlsx, type=fem_ssrm, expected_fs=1.115, element_type=tri6, target_size=1.3, tolerance=0.02, f_min=0.9, f_max=1.45, max_iter=4000, benchmark=RS2-12 -->

![RS2-12: FEM model (left) and maximum shear strain contours at the critical SRF (right)](images/RS2-12.png)

### RS2-13: Simple slope III (Yamagami & Ueta) {#rs2-13}

Slide2 counterpart: [VP17](rocscience.md#vp17).

**Input files:** [vp017.xlsx](../files/rocscience/vp017.xlsx)

| Method | XSLOPE | Published |
|---|---|---|
| SSRM | 1.38 | RS2 SSRM 1.33 |

*Cross-bearings: Greco Spencer 1.33; XSLOPE LEM locks Bishop 1.342 / Spencer 1.340 vs Y&U 1.348 / 1.339.*

<!-- test: file=../files/rocscience/vp017.xlsx, type=fem_ssrm, expected_fs=1.384, element_type=tri6, target_size=0.5, tolerance=0.02, f_min=1.1, f_max=1.65, max_iter=4000, benchmark=RS2-13 -->

![RS2-13: FEM model (left) and maximum shear strain contours at the critical SRF (right)](images/RS2-13.png)

### RS2-14: Simple slope, pore pressure by r<sub>u</sub> {#rs2-14}

Slide2 counterpart: [VP18](rocscience.md#vp18) (this problem is Slide2 VP18, not VP21). Built
with a caveat.

**Input files:** [vp018.xlsx](../files/rocscience/vp018.xlsx)

| Method | XSLOPE | Published |
|---|---|---|
| SSRM (regression lock at the 2.0 m mesh) | 0.948 | RS2 SSRM 0.98 |

*Cross-bearings: Slide2 Spencer 1.01; Baker 1.02; XSLOPE LEM locks Spencer 1.033 on the same file.*

The FEM now carries the r<sub>u</sub> option, but the SSRM factor on *this* model does not
become mesh-independent: 0.986 → 0.948 → 0.902 → 0.873 as the target size goes 2.8 → 2.0 →
1.4 → 1.0 m, with no plateau. The tag pins 2.0 m (0.948) as a regression lock, chosen
mid-sweep rather than at the coarse end that happens to sit on RS2's 0.98 — the honest
reading is a value between roughly 0.87 and 0.99, straddling RS2 SSRM 0.98, Slide2 Spencer
1.01 and Baker 1.02.

The drift is a property of the model, not the r<sub>u</sub> plumbing: run the same slope dry
and the same meshes converge (2.127 → 2.135, +0.4%). With r<sub>u</sub> = 0.5 half the
overburden is cancelled, leaving so little effective confinement that the shear band keeps
localizing as the elements shrink — unregularized Mohr-Coulomb has no length scale to stop
it, and a tension cutoff changes nothing (0.987/0.948/0.901/0.870). LEM locks Spencer 1.033
on the same file.

<!-- test: file=../files/rocscience/vp018.xlsx, type=fem_ssrm, expected_fs=0.948, element_type=tri6, target_size=2.0, tolerance=0.02, f_min=0.7, f_max=1.3, max_iter=4000, benchmark=RS2-14 -->

![RS2-14: FEM model (left) and maximum shear strain contours at the critical SRF (right)](images/RS2-14.png)

### RS2-15: Layered slope II (Greco ex. 4 / Yamagami & Ueta) {#rs2-15}

Slide2 counterpart: [VP19](rocscience.md#vp19).

**Input files:** [vp019.xlsx](../files/rocscience/vp019.xlsx)

| Method | XSLOPE | Published |
|---|---|---|
| SSRM | 1.39 | RS2 SSRM 1.39 |

*Mesh-converged (1.386→1.377 over a 1.7× size change). Cross-bearings: Slide2 Spencer 1.398; Greco 1.40–1.42.*

<!-- test: file=../files/rocscience/vp019.xlsx, type=fem_ssrm, expected_fs=1.386, element_type=tri6, target_size=4.33, tolerance=0.02, f_min=1.1, f_max=1.7, max_iter=4000, benchmark=RS2-15 -->

![RS2-15: FEM model (left) and maximum shear strain contours at the critical SRF (right)](images/RS2-15.png)

### RS2-16: Layered slope and water table with weak seam (Greco ex. 5 / Chen & Shao) {#rs2-16}

Slide2 counterpart: [VP20](rocscience.md#vp20).

**Input files:** [vp020.xlsx](../files/rocscience/vp020.xlsx)

| Method | XSLOPE | Published |
|---|---|---|
| SSRM | 0.972 | RS2 SSRM 1.02 |

*Mesh sweep: 0.983 at 4.0 m, 0.961 at 2.2 m. Cross-bearings: Slide2 Spencer 1.093 circular / 1.007 noncircular; Greco 0.973–1.1; XSLOPE LEM locks 1.086–1.091 on the same file.*

The SSRM used to fall through its bracket here: the model's base is an *inclined* polygon
boundary, and the FEM fixed displacements only along the nodes at the single lowest
elevation, so the body hung from one corner and never reached equilibrium at any F. Fixing
the whole bottom polyline (see [#22](#rs2-22)) resolved it.

<!-- test: file=../files/rocscience/vp020.xlsx, type=fem_ssrm, expected_fs=0.972, element_type=tri6, target_size=3.0, tolerance=0.02, f_min=0.8, f_max=1.4, max_iter=4000, benchmark=RS2-16 -->

![RS2-16: FEM model (left) and maximum shear strain contours at the critical SRF (right)](images/RS2-16.png)

### RS2-17: Slope with three pore pressure conditions (Fredlund & Krahn) {#rs2-17}

Slide2 counterpart: **VP21** (inventory-only on the LEM page — no detail section to link).
Built for the dry and r<sub>u</sub> cases.

**Input files:** [vp021a.xlsx](../files/rocscience/vp021a.xlsx),
[vp021b.xlsx](../files/rocscience/vp021b.xlsx)

| Method | XSLOPE | Published |
|---|---|---|
| SSRM (vp021a, dry) | 2.02 | RS2 SSRM 2.0 |
| SSRM (vp021b, r<sub>u</sub> = 0.25) | 1.704 | — (RS2's table records no SSRM for this sub-case) |

*Cross-bearings — dry: Slide2 M-P 2.075, F&K 2.076. r<sub>u</sub> = 0.25: F&K's LEM 1.761–1.766 and Slide2's 1.760–1.763; mesh-stable (1.696 at 2.0 m) — the usual few percent of SSRM-under-LEM.*

Because RS2's table does not record its own SSRM for the r<sub>u</sub> sub-case, that
cross-check is still open. The water-table case awaits the VP21 case-3 input file.

<!-- test: file=../files/rocscience/vp021a.xlsx, type=fem_ssrm, expected_fs=2.018, element_type=tri6, target_size=3.0, tolerance=0.02, f_min=1.6, f_max=2.5, max_iter=4000, benchmark=RS2-17 -->
<!-- test: file=../files/rocscience/vp021b.xlsx, type=fem_ssrm, expected_fs=1.704, element_type=tri6, target_size=3.0, tolerance=0.02, f_min=1.2, f_max=2.2, max_iter=4000, benchmark=RS2-17b -->

![RS2-17: dry case (vp021a) — FEM model (left) and maximum shear strain contours at the critical SRF (right)](images/RS2-17.png)

![RS2-17b: r<sub>u</sub> = 0.25 case (vp021b) — FEM model (left) and maximum shear strain contours at the critical SRF (right)](images/RS2-17b.png)

### RS2-18: Three pore pressure conditions and a weak seam (Fredlund & Krahn) {#rs2-18}

Slide2 counterpart: [VP22](rocscience.md#vp22). Built for the dry and r<sub>u</sub> cases.

**Input files:** [vp022a.xlsx](../files/rocscience/vp022a.xlsx),
[vp022b.xlsx](../files/rocscience/vp022b.xlsx)

| Method | XSLOPE | Published |
|---|---|---|
| SSRM (vp022a, dry) | 1.34 | RS2 SSRM 1.34 |
| SSRM (vp022b, r<sub>u</sub> = 0.25) | 1.054 | — (RS2's table records no SSRM for this sub-case) |

*Cross-bearings — dry: Slide2 Bishop 1.382. r<sub>u</sub> = 0.25: Slide2 1.124 and F&K 1.124.*

This one returns the *same* factor at 3.0 m and 2.0 m — the mechanism is pinned by the weak
seam, a geometric feature, so it cannot migrate with refinement. The contrast with
[#14](#rs2-14) is the point: there, nothing pins the band. Same open RS2 sub-case
cross-check as [#17](#rs2-17); water-table case likewise pending.

<!-- test: file=../files/rocscience/vp022a.xlsx, type=fem_ssrm, expected_fs=1.336, element_type=tri6, target_size=3.0, tolerance=0.02, f_min=1.0, f_max=1.7, max_iter=4000, benchmark=RS2-18 -->
<!-- test: file=../files/rocscience/vp022b.xlsx, type=fem_ssrm, expected_fs=1.054, element_type=tri6, target_size=3.0, tolerance=0.02, f_min=0.8, f_max=1.8, max_iter=4000, benchmark=RS2-18b -->

![RS2-18: dry case (vp022a) — FEM model (left) and maximum shear strain contours at the critical SRF (right)](images/RS2-18.png)

![RS2-18b: r<sub>u</sub> = 0.25 case (vp022b) — FEM model (left) and maximum shear strain contours at the critical SRF (right)](images/RS2-18b.png)

### RS2-19: Undrained layered slope (Low 1989) {#rs2-19}

Slide2 counterpart: [VP24](rocscience.md#vp24) (this problem is Slide2 VP24). Built with a
caveat.

**Input files:** [vp024.xlsx](../files/rocscience/vp024.xlsx)

| Method | XSLOPE | Published |
|---|---|---|
| SSRM | 1.51 at the tagged mesh (1.48 at 0.6 m) | RS2 SSRM 1.41 |

*Cross-bearings: Slide2 LEM 1.439; Low 1.44.*

The two SSRM values straddle the LEM from opposite sides on this φ = 0 slope, and the XSLOPE
factor drifts −2% with refinement; quoted at the tagged mesh per the page convention.

<!-- test: file=../files/rocscience/vp024.xlsx, type=fem_ssrm, expected_fs=1.507, element_type=tri6, target_size=1.0, tolerance=0.02, f_min=1.1, f_max=1.8, max_iter=4000, benchmark=RS2-19 -->

![RS2-19: FEM model (left) and maximum shear strain contours at the critical SRF (right)](images/RS2-19.png)

### RS2-20: Slope with vertical load (Prandtl's wedge) {#rs2-20}

Slide2 counterpart: [VP25](rocscience.md#vp25).

**Input files:** [vp025.xlsx](../files/rocscience/vp025.xlsx)

| Method | XSLOPE | Published |
|---|---|---|
| SSRM | 1.00 (mesh pair 1.011→1.003) | RS2 SSRM 1.0 |

*Cross-bearings: Prandtl theory 1.0; Slide2 Spencer reads 1.051 on the specified surface.*

<!-- test: file=../files/rocscience/vp025.xlsx, type=fem_ssrm, expected_fs=1.003, element_type=tri6, target_size=0.8, tolerance=0.01, f_min=0.5, f_max=1.6, max_iter=4000, benchmark=RS2-20 -->

![RS2-20: FEM model (left) and maximum shear strain contours at the critical SRF (right)](images/RS2-20.png)

### RS2-21: Bearing capacity test prism (Prandtl II) {#rs2-21}

Slide2 counterpart: **VP26** (inventory-only on the LEM page — no detail section to link).

**Input files:** [vp026.xlsx](../files/rocscience/vp026.xlsx)

| Method | XSLOPE | Published |
|---|---|---|
| SSRM | 1.011 (1.020 at 1.2 m) | RS2 SSRM 1.01 |

*Converging on Prandtl theory 1.0. Cross-bearings: Slide2 Spencer 0.941 on the specified surface.*

The input file was extracted for this problem.

<!-- test: file=../files/rocscience/vp026.xlsx, type=fem_ssrm, expected_fs=1.011, element_type=tri6, target_size=0.8, tolerance=0.01, f_min=0.5, f_max=1.6, max_iter=4000, benchmark=RS2-21 -->

![RS2-21: FEM model (left) and maximum shear strain contours at the critical SRF (right)](images/RS2-21.png)

### RS2-22: Layered slope with undulating bedrock {#rs2-22}

Slide2 counterpart: [VP27](rocscience.md#vp27). Built on an SSRM variant.

**Input files:** [vp027_fem.xlsx](../files/rocscience/vp027_fem.xlsx)

| Method | XSLOPE | Published |
|---|---|---|
| SSRM | 1.558 (1.536 at 1.8 m, converging) | RS2 SSRM 1.52 |

*+2.5% at the tagged mesh, +1.1% at the finer one.*

Two FEM gaps had to close first: displacements are now fixed along the whole bottom
*polyline* of the domain rather than only the nodes at the lowest elevation (an undulating
bedrock base used to hang from one corner and never reach equilibrium — this also unblocked
[#16](#rs2-16)), and the FEM now applies the phreatic-inclination (Hu) correction that this
problem specifies, matching the LEM's Type=phreatic path.

The SSRM runs on [vp027_fem.xlsx](../files/rocscience/vp027_fem.xlsx), which differs from the
LEM file in one respect: the published model caps the crest with a **zero-strength** layer
(c = 0, φ = 0). That is a limit-equilibrium device — dead weight riding above the failure
surface — and it has no continuum equivalent, because a null Mohr-Coulomb yield surface
stays null however much you divide it by, so the cap can never carry deviatoric stress and
the SSRM correctly reports no equilibrium at any F (it is not a solver defect). The variant
gives the cap the parent soil's strength; everything else is untouched. vp027's LEM locks
stand on the as-published file.

<!-- test: file=../files/rocscience/vp027_fem.xlsx, type=fem_ssrm, expected_fs=1.558, element_type=tri6, target_size=2.5, tolerance=0.02, f_min=1.2, f_max=1.9, max_iter=4000, benchmark=RS2-22 -->

![RS2-22: FEM model (left) and maximum shear strain contours at the critical SRF (right)](images/RS2-22.png)

### RS2-23: Underwater slope with linearly varying cohesion {#rs2-23}

Slide2 counterpart: [VP29](rocscience.md#vp29).

| Method | XSLOPE | Published |
|---|---|---|
| SSRM | *no lock possible* | RS2 SSRM 1.12 |

*Stand-ins for the two readings of the "can't fail" region give 0.87 and 0.92; the unsuppressed SSRM minimum is 0.21. The slope's anchor remains the LEM lock ([VP29](rocscience.md#vp29), Spencer 1.145 on Duncan's surface).*

RS2's published SSRM depends on a "can't fail" elastic region whose boundary its text and
figure draw differently; without the patch the true SSRM minimum is the shallow skin above
el. −20 (FS 0.21) that the artifice suppresses. The comparison would test where the patch is
drawn, not the mechanics.

### RS2-24: Layered slope with geosynthetic reinforcement {#rs2-24}

Slide2 counterpart: [VP32](rocscience.md#vp32).

**Input files:** [vp032a.xlsx](../files/rocscience/vp032a.xlsx),
[vp032c.xlsx](../files/rocscience/vp032c.xlsx)

| Method | XSLOPE | Published |
|---|---|---|
| SSRM (vp032a, H = 7) | 1.196 (mesh-stable 1.202 at 2.2 m) | RS2 SSRM 1.15 |
| SSRM (vp032c, H = 8.75) | 0.991 | RS2 SSRM 0.95 |

*+4%, within the corpus SSRM scatter; geotextile as an FEM truss with EA = 2×10³ kN/m stated as convention.*

RS2's fully labeled figures also supplied the geometry that unlocked Slide2's
[VP32](rocscience.md#vp32) — LEM locks on the three printed circles live there.

<!-- test: file=../files/rocscience/vp032a.xlsx, type=fem_ssrm, expected_fs=1.196, element_type=tri6, target_size=1.5, tolerance=0.02, f_min=0.9, f_max=1.6, max_iter=4000, benchmark=RS2-24a -->
<!-- test: file=../files/rocscience/vp032c.xlsx, type=fem_ssrm, expected_fs=0.991, element_type=tri6, target_size=2.2, tolerance=0.02, f_min=0.7, f_max=1.4, max_iter=4000, benchmark=RS2-24b -->

![RS2-24a: H = 7 case (vp032a, SSRM 1.196) — FEM model (left) and maximum shear strain contours at the critical SRF (right)](images/RS2-24a.png)

![RS2-24b: H = 8.75 case (vp032c, SSRM 0.991) — FEM model (left) and maximum shear strain contours at the critical SRF (right)](images/RS2-24b.png)

### RS2-25: Syncrude tailings dyke (El-Ramly et al. 2003) {#rs2-25}

Slide2 counterpart: [VP33](rocscience.md#vp33). Built with a caveat.

**Input files:** [vp033.xlsx](../files/rocscience/vp033.xlsx)

| Method | XSLOPE | Published |
|---|---|---|
| SSRM | 1.20 (mesh-converged) | RS2 SSRM 1.29 |

*Cross-bearings: Slide2 Bishop 1.305; El-Ramly 1.31; XSLOPE's own LEM 1.299 on Slide's circle and 1.253 on a free composite search.*

The weak presheared clay-shale rewards less-constrained searches (XSLOPE's own LEM: 1.299 on
Slide's circle, 1.253 free composite search, SSRM 1.20); locked at the tagged mesh.

<!-- test: file=../files/rocscience/vp033.xlsx, type=fem_ssrm, expected_fs=1.199, element_type=tri6, target_size=5.0, tolerance=0.02, f_min=0.9, f_max=1.8, max_iter=4000, benchmark=RS2-25 -->

![RS2-25: FEM model (left) and maximum shear strain contours at the critical SRF (right)](images/RS2-25.png)

### RS2-26: Clarence Cannon dam (Wolff & Harr 1987) {#rs2-26}

Slide2 counterpart: [VP34](rocscience.md#vp34).

**Input files:** [vp034.xlsx](../files/rocscience/vp034.xlsx)

| Method | XSLOPE | Published |
|---|---|---|
| SSRM | 2.28 | RS2 SSRM 2.29 (−0.4%) |

*Cross-bearings: Slide2 GLE 2.333 / Spencer 2.383; W&H 2.36; XSLOPE LEM M-P 2.384.*

Polygon zones with the chimney drain.

<!-- test: file=../files/rocscience/vp034.xlsx, type=fem_ssrm, expected_fs=2.281, element_type=tri6, target_size=15.0, tolerance=0.02, f_min=1.7, f_max=3.0, max_iter=4000, benchmark=RS2-26 -->

![RS2-26: FEM model (left) and maximum shear strain contours at the critical SRF (right)](images/RS2-26.png)

### RS2-27: Homogeneous slope, pore pressure by r<sub>u</sub> {#rs2-27}

Slide2 counterpart: [VP36](rocscience.md#vp36). *Blocked*: the FEM has no r<sub>u</sub> option
(task, with RS2 [#14](#rs2-14)); RS2's own text cites Slide2 [VP36](rocscience.md#vp36)
(Li & Lumb), not VP21. No published/XSLOPE FS pair to tabulate.

### RS2-29: Geosynthetic-reinforced embankment on soft soil (Tandjiria 2002) {#rs2-29}

Slide2 counterpart: [VP39](rocscience.md#vp39). The manual's §29/§30 headings are swapped.
Built for the sand case.

**Input files:** [vp039c.xlsx](../files/rocscience/vp039c.xlsx)

| Method | XSLOPE | Published |
|---|---|---|
| SSRM (sand case, vp039c) | 1.253 | RS2 SSRM 1.25 |

*Cross-bearings: Spencer 1.209; Tandjiria 1.219.*

The clay case (RS2 SSRM 0.99) is not locked: its water-filled tension crack has no FEM
representation (XSLOPE reads 1.05 without it).

<!-- test: file=../files/rocscience/vp039c.xlsx, type=fem_ssrm, expected_fs=1.253, element_type=tri6, target_size=0.7, tolerance=0.02, f_min=0.9, f_max=1.7, max_iter=4000, benchmark=RS2-29 -->

![RS2-29: sand case (vp039c) — FEM model (left) and maximum shear strain contours at the critical SRF (right)](images/RS2-29.png)

### RS2-30: Homogeneous slope, power-curve strength (Perry 1993) {#rs2-30}

Slide2 counterpart: [VP40](rocscience.md#vp40). Swapped heading (see [#29](#rs2-29)).

**Input files:** [vp040.xlsx](../files/rocscience/vp040.xlsx)

| Method | XSLOPE | Published |
|---|---|---|
| SSRM | 0.910 | RS2 SRF 0.91 (exact) |

*Cross-bearings: Slide2 Janbu 0.944; Perry 0.98.*

The FEM linearizes the reduced envelope at the current stress per iteration.

<!-- test: file=../files/rocscience/vp040.xlsx, type=fem_ssrm, expected_fs=0.910, element_type=tri6, target_size=2.0, tolerance=0.02, f_min=0.5, f_max=1.5, max_iter=4000, benchmark=RS2-30 -->

![RS2-30: FEM model (left) and maximum shear strain contours at the critical SRF (right)](images/RS2-30.png)

### RS2-31: M-C vs power curve (Baker 2003 ex. 1) {#rs2-31}

Slide2 counterpart: [VP44](rocscience.md#vp44). Built, all three halves.

**Input files:** [vp044a.xlsx](../files/rocscience/vp044a.xlsx),
[vp044b.xlsx](../files/rocscience/vp044b.xlsx),
[vp044c.xlsx](../files/rocscience/vp044c.xlsx)

| Method | XSLOPE | Published |
|---|---|---|
| SSRM (vp044b, M-C) | 1.546 | RS2 SSRM 1.53 |
| SSRM (vp044c, M-C) | 0.966 | RS2 SSRM 0.98 |
| SSRM (vp044a, power curve) | 0.934 | RS2 SSRM 1.11 |

*Cross-bearings on the power-curve case: Slide2's own published band is Janbu 0.92 / Spencer 0.96, Baker ~0.97.*

XSLOPE's power-curve SSRM sits squarely on Slide2's band; RS2's 1.11 sits 15% above Slide2's
LEM on the same problem and is the outlier.

<!-- test: file=../files/rocscience/vp044b.xlsx, type=fem_ssrm, expected_fs=1.546, element_type=tri6, target_size=0.5, tolerance=0.02, f_min=1.1, f_max=2.0, max_iter=4000, benchmark=RS2-31a -->
<!-- test: file=../files/rocscience/vp044c.xlsx, type=fem_ssrm, expected_fs=0.966, element_type=tri6, target_size=0.5, tolerance=0.02, f_min=0.6, f_max=1.4, max_iter=4000, benchmark=RS2-31b -->
<!-- test: file=../files/rocscience/vp044a.xlsx, type=fem_ssrm, expected_fs=0.934, element_type=tri6, target_size=0.5, tolerance=0.02, f_min=0.5, f_max=1.6, max_iter=4000, benchmark=RS2-31c -->

![RS2-31a: Mohr-Coulomb case (vp044b, SSRM 1.546) — FEM model (left) and maximum shear strain contours at the critical SRF (right)](images/RS2-31a.png)

![RS2-31b: Mohr-Coulomb case (vp044c, SSRM 0.966) — FEM model (left) and maximum shear strain contours at the critical SRF (right)](images/RS2-31b.png)

![RS2-31c: power-curve case (vp044a, SSRM 0.934) — FEM model (left) and maximum shear strain contours at the critical SRF (right)](images/RS2-31c.png)

### RS2-32: Heading mismatch — body is Baker's example 2 {#rs2-32}

Slide2 counterpart: [VP45](rocscience.md#vp45). Built, both halves.

**Input files:** [vp045a.xlsx](../files/rocscience/vp045a.xlsx),
[vp045b.xlsx](../files/rocscience/vp045b.xlsx)

| Method | XSLOPE | Published |
|---|---|---|
| SSRM (vp045a, M-C) | 2.852 | RS2 SSRM 2.83 (+0.8%) |
| SSRM (vp045b, power curve) | 2.732 | RS2 SSRM 2.74 (−0.3%) |

*Cross-bearings: Slide2 Spencer 2.662.*

<!-- test: file=../files/rocscience/vp045a.xlsx, type=fem_ssrm, expected_fs=2.852, element_type=tri6, target_size=1.0, tolerance=0.02, f_min=2.3, f_max=3.4, max_iter=4000, benchmark=RS2-32 -->
<!-- test: file=../files/rocscience/vp045b.xlsx, type=fem_ssrm, expected_fs=2.732, element_type=tri6, target_size=1.0, tolerance=0.02, f_min=1.8, f_max=3.6, max_iter=4000, benchmark=RS2-32b -->

![RS2-32: Mohr-Coulomb case (vp045a, SSRM 2.852) — FEM model (left) and maximum shear strain contours at the critical SRF (right)](images/RS2-32.png)

![RS2-32b: power-curve case (vp045b, SSRM 2.732) — FEM model (left) and maximum shear strain contours at the critical SRF (right)](images/RS2-32b.png)

### RS2-33: Homogeneous slope with tension crack and water table (P&D test slope 2) {#rs2-33}

Slide2 counterpart: [VP56](rocscience.md#vp56). Swapped heading. Built with a caveat.

**Input files:** [vp056.xlsx](../files/rocscience/vp056.xlsx)

| Method | XSLOPE | Published |
|---|---|---|
| SSRM | 1.247 | RS2 SSRM 1.28 |

*Cross-bearings: an eight-program LEM table spanning 1.03–1.32.*

The model's dry tension crack has no FEM representation, worth ~2–3% here.

<!-- test: file=../files/rocscience/vp056.xlsx, type=fem_ssrm, expected_fs=1.247, element_type=tri6, target_size=4.0, tolerance=0.02, f_min=0.9, f_max=1.7, max_iter=4000, benchmark=RS2-33 -->

![RS2-33: FEM model (left) and maximum shear strain contours at the critical SRF (right)](images/RS2-33.png)

### RS2-34: M-C vs power curve III (Baker 2003 ex. 3, London clay) {#rs2-34}

Slide2 counterpart: [VP61](rocscience.md#vp61). Built, both halves.

**Input files:** [vp061a.xlsx](../files/rocscience/vp061a.xlsx),
[vp061b.xlsx](../files/rocscience/vp061b.xlsx)

| Method | XSLOPE | Published |
|---|---|---|
| SSRM (vp061b, M-C) | 1.369 | RS2 SSRM 1.38 |
| SSRM (vp061a, power curve) | 1.502 | RS2 SSRM 1.47 (+2.1%) |

*Cross-bearings on the power-curve case: Slide2 Spencer 1.47; Baker 1.48.*

<!-- test: file=../files/rocscience/vp061b.xlsx, type=fem_ssrm, expected_fs=1.369, element_type=tri6, target_size=0.5, tolerance=0.02, f_min=1.0, f_max=1.9, max_iter=4000, benchmark=RS2-34 -->
<!-- test: file=../files/rocscience/vp061a.xlsx, type=fem_ssrm, expected_fs=1.502, element_type=tri6, target_size=0.5, tolerance=0.02, f_min=1.0, f_max=2.2, max_iter=4000, benchmark=RS2-34b -->

![RS2-34: Mohr-Coulomb case (vp061b, SSRM 1.369) — FEM model (left) and maximum shear strain contours at the critical SRF (right)](images/RS2-34.png)

![RS2-34b: power-curve case (vp061a, SSRM 1.502) — FEM model (left) and maximum shear strain contours at the critical SRF (right)](images/RS2-34b.png)

### RS2-36: Seepage analysis, homogeneous slope (D&W Fig 6.37) {#rs2-36}

Slide2 counterpart: [VP71](rocscience.md#vp71) (= Slide2 VP71, not
[VP70](rocscience.md#vp70)). Built, both cases.

**Input files:** [vp071a.xlsx](../files/rocscience/vp071a.xlsx),
[vp071b.xlsx](../files/rocscience/vp071b.xlsx)

| Method | XSLOPE | Published |
|---|---|---|
| SSRM (vp071a, FE seepage) | 1.118 | RS2 SSRM 1.12 |
| SSRM (vp071b, piezo approximation) | 1.118 | RS2 SSRM 1.12 |

*Cross-bearings: referee 1.138/1.141; XSLOPE LEM locks 1.132.*

The seep case runs on tri6 sidecars.

<!-- test: file=../files/rocscience/vp071a.xlsx, type=fem_ssrm, expected_fs=1.118, tolerance=0.01, f_min=0.7, f_max=1.6, max_iter=4000, benchmark=RS2-36a -->
<!-- test: file=../files/rocscience/vp071b.xlsx, type=fem_ssrm, expected_fs=1.118, element_type=tri6, target_size=4.4, tolerance=0.01, f_min=0.7, f_max=1.6, max_iter=4000, benchmark=RS2-36b -->

![RS2-36a: FE-seepage case (vp071a) — FEM model (left) and maximum shear strain contours at the critical SRF (right)](images/RS2-36a.png)

![RS2-36b: piezometric-line case (vp071b) — FEM model (left) and maximum shear strain contours at the critical SRF (right)](images/RS2-36b.png)

### RS2-37: Embankment with layered foundation (D&W Fig 6.39) {#rs2-37}

Slide2 counterpart: [VP72](rocscience.md#vp72). Reported, no lock.

| Method | XSLOPE | Published |
|---|---|---|
| SSRM | 1.31 (deep mechanism) | RS2 SSRM 0.95 in its table, 1.1 in its own convergence graph (artesian downstream-toe slide) |

*Cross-bearings: XSLOPE LEM on the tangent circle 1.339; Slide2 1.15/1.16; referee 1.11.*

The two programs are not finding the same mechanism. Reproducing the toe mechanism needs
toe-refined meshing — noted with the artesian-toe discussion in the Slide2
[VP72](rocscience.md#vp72) section.

### RS2-40: Dam with impermeable foundation (D&W Fig 7.24) {#rs2-40}

Slide2 counterpart: [VP77](rocscience.md#vp77). Built for the piezo case.

**Input files:** [vp077b.xlsx](../files/rocscience/vp077b.xlsx)

| Method | XSLOPE | Published |
|---|---|---|
| SSRM (vp077b, piezo case) | 1.126 | RS2 SSRM 1.53 |

As with the dry Talbingo dam (RS2-4), the cohesionless downstream shell fails as a
surface-parallel skin rather than a deep rotation. Here the piezometric line daylights at the
toe, so the governing mechanism is the *saturated* toe: the seepage-parallel infinite-slope
limit is (140 − 62.4)/140 × tan 38° / tan 20° = 1.190, and the per-node SSRM returns 1.126. The
per-node out-of-balance nodes sit exactly on the daylighting toe (x ∈ [1162, 1234]). The finite,
partially-saturated toe geometry softens the FEM ~5% below the idealized infinite slope, so this
skin's analytic anchor is looser than RS2-4's ±0.5% — a φ-sweep clincher (as in RS2-4) is the
natural way to tighten it. The published RS2 SSRM 1.53 reports a deeper mechanism on the face;
XSLOPE reports the more critical toe skin as the true global minimum.

*The FE-seepage case is blocked — tri6 seepage does not converge on the high-contrast thick core (the documented tri3/tri6 trade).*

<!-- test: file=../files/rocscience/vp077b.xlsx, type=fem_ssrm, expected_fs=1.126, element_type=tri6, target_size=12.4, tolerance=0.02, f_min=1.1, f_max=2.2, max_iter=16000, benchmark=RS2-40 -->

![RS2-40: piezometric case (vp077b) — FEM model (left) and maximum shear strain contours at the critical SRF (right)](images/RS2-40.png)

### RS2-42: James dike {#rs2-42}

Slide2 counterpart: [VP75](rocscience.md#vp75).

**Input files:** [vp075.xlsx](../files/rocscience/vp075.xlsx)

| Method | XSLOPE | Published |
|---|---|---|
| SSRM | 1.249 | RS2 SSRM 1.26 (−0.9%) |

*Cross-bearings: Slide2 noncircular LEM 1.11–1.16; referee 1.17.*

<!-- test: file=../files/rocscience/vp075.xlsx, type=fem_ssrm, expected_fs=1.249, element_type=tri6, target_size=1.85, tolerance=0.02, f_min=0.8, f_max=1.8, max_iter=4000, benchmark=RS2-42 -->

![RS2-42: FEM model (left) and maximum shear strain contours at the critical SRF (right)](images/RS2-42.png)

### RS2-44: Seepage analysis for an earth embankment (D&W Fig 14.20-a) {#rs2-44}

Slide2 counterpart: [VP82](rocscience.md#vp82) (= Slide2 VP82, not
[VP76](rocscience.md#vp76) — §39's body carries VP76).

**Input files:** [vp082.xlsx](../files/rocscience/vp082.xlsx)

| Method | XSLOPE | Published |
|---|---|---|
| SSRM | 1.511 | RS2 SSRM 1.51 (+0.1%) |

*Cross-bearings: Slide2 LEM 1.532/1.541; referee 1.528–1.542.*

<!-- test: file=../files/rocscience/vp082.xlsx, type=fem_ssrm, expected_fs=1.511, element_type=tri6, target_size=2.0, tolerance=0.02, f_min=1.0, f_max=2.1, max_iter=4000, benchmark=RS2-44 -->

![RS2-44: FEM model (left) and maximum shear strain contours at the critical SRF (right)](images/RS2-44.png)

### RS2-45: Varying undrained shear strength profiles (D&W Fig 14.20-b) {#rs2-45}

Slide2 counterpart: [VP83](rocscience.md#vp83). Built with a caveat.

**Input files:** [vp083a.xlsx](../files/rocscience/vp083a.xlsx),
[vp083b.xlsx](../files/rocscience/vp083b.xlsx)

| Method | XSLOPE | Published |
|---|---|---|
| SSRM (vp083a) | 1.41 | RS2 SSRM 1.32 |
| SSRM (vp083b) | 1.37 | RS2 SSRM 1.32 |

*Cross-bearings: D&W referee 1.28–1.33.*

XSLOPE's SSRM reads high on φ=0 foundations (the [RS2-19](#rs2-19) family, under
investigation); locked as regression values at the tagged mesh.

<!-- test: file=../files/rocscience/vp083a.xlsx, type=fem_ssrm, expected_fs=1.412, element_type=tri6, target_size=4.0, tolerance=0.02, f_min=0.9, f_max=1.9, max_iter=4000, benchmark=RS2-45a -->
<!-- test: file=../files/rocscience/vp083b.xlsx, type=fem_ssrm, expected_fs=1.365, element_type=tri6, target_size=4.0, tolerance=0.02, f_min=0.9, f_max=1.9, max_iter=4000, benchmark=RS2-45b -->

![RS2-45a: vp083a (SSRM 1.412) — FEM model (left) and maximum shear strain contours at the critical SRF (right)](images/RS2-45a.png)

![RS2-45b: vp083b (SSRM 1.365) — FEM model (left) and maximum shear strain contours at the critical SRF (right)](images/RS2-45b.png)

### RS2-46: Varying undrained strength profiles II (D&W Fig 15.9, c<sub>u</sub> = 300 + c<sub>z</sub>·z) {#rs2-46}

Slide2 counterpart: [VP84](rocscience.md#vp84).

**Input files:** [vp084a–d](../files/rocscience/vp084a.xlsx)

| Method | XSLOPE | Published |
|---|---|---|
| SSRM (vp084a) | 0.80 | RS2 SSRM 0.78 |
| SSRM (vp084b) | 0.95 | RS2 SSRM 0.93 |
| SSRM (vp084c) | 1.08 | RS2 SSRM 1.05 |
| SSRM (vp084d) | 1.19 | RS2 SSRM 1.15 |

*+2–3%, the φ=0 pattern. Cross-bearings: D&W 0.75 / 0.90 / 1.03 / 1.13 for the four cases in the same order.*

<!-- test: file=../files/rocscience/vp084a.xlsx, type=fem_ssrm, expected_fs=0.804, element_type=tri6, target_size=4.0, tolerance=0.02, f_min=0.4, f_max=1.3, max_iter=4000, benchmark=RS2-46a -->
<!-- test: file=../files/rocscience/vp084b.xlsx, type=fem_ssrm, expected_fs=0.947, element_type=tri6, target_size=4.0, tolerance=0.02, f_min=0.5, f_max=1.4, max_iter=4000, benchmark=RS2-46b -->
<!-- test: file=../files/rocscience/vp084c.xlsx, type=fem_ssrm, expected_fs=1.082, element_type=tri6, target_size=4.0, tolerance=0.02, f_min=0.6, f_max=1.5, max_iter=4000, benchmark=RS2-46c -->
<!-- test: file=../files/rocscience/vp084d.xlsx, type=fem_ssrm, expected_fs=1.188, element_type=tri6, target_size=4.0, tolerance=0.02, f_min=0.7, f_max=1.7, max_iter=4000, benchmark=RS2-46d -->

![RS2-46a: vp084a (SSRM 0.804) — FEM model (left) and maximum shear strain contours at the critical SRF (right)](images/RS2-46a.png)

![RS2-46b: vp084b (SSRM 0.947) — FEM model (left) and maximum shear strain contours at the critical SRF (right)](images/RS2-46b.png)

![RS2-46c: vp084c (SSRM 1.082) — FEM model (left) and maximum shear strain contours at the critical SRF (right)](images/RS2-46c.png)

![RS2-46d: vp084d (SSRM 1.188) — FEM model (left) and maximum shear strain contours at the critical SRF (right)](images/RS2-46d.png)

### RS2-47: Purely cohesive slope, varying thickness (D&W Fig 14.3) {#rs2-47}

Slide2 counterpart: [VP78](rocscience.md#vp78). Built for the 30-ft case.

**Input files:** [vp078.xlsx](../files/rocscience/vp078.xlsx)

| Method | XSLOPE | Published |
|---|---|---|
| SSRM (30-ft case) | 1.08 | RS2 SSRM 1.03 |

*Cross-bearings: D&W referee 1.124–1.135.*

The 46.5- and 60-ft variants (RS2 1.02 / 1.02) need deepened-base builds and are deferred.

<!-- test: file=../files/rocscience/vp078.xlsx, type=fem_ssrm, expected_fs=1.081, element_type=tri6, target_size=4.0, tolerance=0.02, f_min=0.6, f_max=1.6, max_iter=4000, benchmark=RS2-47 -->

![RS2-47: 30-ft case (vp078) — FEM model (left) and maximum shear strain contours at the critical SRF (right)](images/RS2-47.png)

### RS2-48–54: Multi-tiered geotextile walls (Leshchinsky & Han 2004) {#rs2-48}

Slide2 counterparts: [VP87](rocscience.md#vp87)–VP93 (one-for-one, verified; only VP87 has a
detail section on the LEM page). *Blocked*.

| Method | XSLOPE | Published |
|---|---|---|
| SSRM (baseline wall, reinforced — any T<sub>a</sub>) | 1.72 | ≈1.0 |
| SSRM (baseline wall, unreinforced) | 0.41 | — |

*The FEM reinforcement tensile-capacity cap is not enforced in SSRM, which is why the SSRM is provably T<sub>a</sub>-insensitive on the baseline wall.*

Extraction and published tables are complete for all seven; the series unlocks together when
the cap is fixed. Problem 55 (five 1.8-m tiers, Slide2 VP94) belongs to this family and is
blocked with the series.

### RS2-56: Homogeneous slope vs Z-Soil, PLAXIS, GEO FEM (Pruska 2003, H = 7 m, 5 cases) {#rs2-56}

New corpus files (no Slide2 counterpart). Built: all five cases run.

**Input files:** [rs2_56a.xlsx](../files/rocscience/rs2_56a.xlsx),
[rs2_56b.xlsx](../files/rocscience/rs2_56b.xlsx)

| Method | XSLOPE | Published |
|---|---|---|
| SSRM (rs2_56a — case 2, weakest; lock) | 0.667 | RS2 SSRM 0.67 |
| SSRM (rs2_56b — case 5, strongest; lock) | 2.131 | RS2 SSRM 2.14 |

*All five cases land within ±3.3% of RS2's M-C and inside the four-program band (Z-Soil, PLAXIS, GEO FEM); the two locks bracket the family. Full case-by-case tables — including the Z-Soil / PLAXIS / GEO FEM / Slide2 columns — are in [the Pruska cross-bearing section](#pruska).*

<!-- test: file=../files/rocscience/rs2_56a.xlsx, type=fem_ssrm, expected_fs=0.667, element_type=tri6, target_size=0.8, tolerance=0.02, f_min=0.32, f_max=1.12, max_iter=4000, benchmark=RS2-56a -->
<!-- test: file=../files/rocscience/rs2_56b.xlsx, type=fem_ssrm, expected_fs=2.131, element_type=tri6, target_size=0.8, tolerance=0.02, f_min=1.79, f_max=2.59, max_iter=4000, benchmark=RS2-56b -->

![RS2-56a: case 2, (γ, c, φ) = (18, 5, 10), the weakest of the five (SSRM 0.667) — FEM model (left) and maximum shear strain contours at the critical SRF (right)](images/RS2-56a.png)

![RS2-56b: case 5, (γ, c, φ) = (24, 20, 30), the strongest of the five (SSRM 2.131) — FEM model (left) and maximum shear strain contours at the critical SRF (right)](images/RS2-56b.png)

### RS2-57: Pruska H = 10.5 m, 6 cases {#rs2-57}

New corpus files. Built: all six cases run.

**Input files:** [rs2_57a.xlsx](../files/rocscience/rs2_57a.xlsx),
[rs2_57b.xlsx](../files/rocscience/rs2_57b.xlsx)

| Method | XSLOPE | Published |
|---|---|---|
| SSRM (rs2_57a — case 1, weakest; lock) | 0.449 | RS2 SSRM 0.44 |
| SSRM (rs2_57b — case 6, strongest; lock) | 1.411 | RS2 SSRM 1.42 |

*All six cases land within ±3.6% of RS2's M-C; the two locks bracket the family. Full case-by-case tables — including the Z-Soil / PLAXIS / GEO FEM / Slide2 columns — are in [the Pruska cross-bearing section](#pruska).*

<!-- test: file=../files/rocscience/rs2_57a.xlsx, type=fem_ssrm, expected_fs=0.449, element_type=tri6, target_size=0.8, tolerance=0.02, f_min=0.1, f_max=0.89, max_iter=4000, benchmark=RS2-57a -->
<!-- test: file=../files/rocscience/rs2_57b.xlsx, type=fem_ssrm, expected_fs=1.411, element_type=tri6, target_size=0.8, tolerance=0.02, f_min=1.07, f_max=1.87, max_iter=4000, benchmark=RS2-57b -->

![RS2-57a: case 1, (γ, c, φ) = (18, 5, 10), the weakest of the six (SSRM 0.449) — FEM model (left) and maximum shear strain contours at the critical SRF (right)](images/RS2-57a.png)

![RS2-57b: case 6, (γ, c, φ) = (24, 20, 30), the strongest of the six (SSRM 1.411) — FEM model (left) and maximum shear strain contours at the critical SRF (right)](images/RS2-57b.png)

### RS2-58: Pruska H = 14 m, 6 cases {#rs2-58}

New corpus files. Built (5 of 6).

**Input files:** [rs2_58a.xlsx](../files/rocscience/rs2_58a.xlsx),
[rs2_58b.xlsx](../files/rocscience/rs2_58b.xlsx)

| Method | XSLOPE | Published |
|---|---|---|
| SSRM (rs2_58a — case 1, weakest; lock) | 0.342 | RS2 SSRM 0.33 |
| SSRM (rs2_58b — case 6, strongest; lock) | 1.057 | RS2 SSRM 1.06 |
| SSRM (case 5, c = 5, φ = 30; unlocked) | 0.667 | RS2 SSRM 0.72 |

*Four of the six land within ±3.6%; the two locks bracket the family. Case 5 reads 0.667 against a tight published 0.72–0.75 cluster (RS2 0.72, Z-Soil 0.75, PLAXIS 0.74, GEO FEM 0.73, Slide2 0.73) and is reported unlocked pending explanation. Full case-by-case tables in [the Pruska cross-bearing section](#pruska).*

<!-- test: file=../files/rocscience/rs2_58a.xlsx, type=fem_ssrm, expected_fs=0.342, element_type=tri6, target_size=0.8, tolerance=0.02, f_min=0.1, f_max=0.78, max_iter=4000, benchmark=RS2-58a -->
<!-- test: file=../files/rocscience/rs2_58b.xlsx, type=fem_ssrm, expected_fs=1.057, element_type=tri6, target_size=0.8, tolerance=0.02, f_min=0.71, f_max=1.51, max_iter=4000, benchmark=RS2-58b -->

![RS2-58a: case 1, (γ, c, φ) = (18, 5, 10), the weakest of the six (SSRM 0.342) — FEM model (left) and maximum shear strain contours at the critical SRF (right)](images/RS2-58a.png)

![RS2-58b: case 6, (γ, c, φ) = (24, 20, 30), the strongest of the six (SSRM 1.057) — FEM model (left) and maximum shear strain contours at the critical SRF (right)](images/RS2-58b.png)

---

## The Pruska cross-bearing (#56–58) {#pruska}

Pruska (2003) analyzed three homogeneous slopes (H = 7, 10.5, and 14 m over an 8-m
foundation) with five or six material sets each in four SSRM programs. RS2 reproduces the
study; XSLOPE's SSRM joins it here as a fifth column — 16 of the 17 cases land within
±4% of RS2 and inside the four-program band. (The study's Drucker-Prager columns are not
comparable; XSLOPE, like Slide2, analyzes Mohr-Coulomb only. Elastic constants are the
paper's published E = 5,000 kPa and per-case ν, not the corpus's usual convention.)

**H = 7 m (#56):** cases (γ, c, φ) = (24,20,10), (18,5,10), (24,20,20), (18,5,20), (24,20,30)

| Case | XSLOPE | RS2 | Z-Soil | PLAXIS | GEO FEM | Slide2 LEM |
|---|---|---|---|---|---|---|
| 1 | 1.254 | 1.22 | 1.21 | 1.22 | 1.31 | 1.22 |
| 2 | 0.667 | 0.67 | 0.71 | 0.68 | 0.73 | 0.66 |
| 3 | 1.689 | 1.68 | 1.64 | 1.65 | 1.71 | 1.64 |
| 4 | 1.016 | 1.05 | 0.95 | 0.99 | 1.17 | 1.02 |
| 5 | 2.131 | 2.14 | 1.98 | 2.09 | 2.19 | 2.08 |

**H = 10.5 m (#57):** cases 1–6 = (18,5,10), (24,20,10), (18,5,20), (24,20,20), (18,5,30), (24,20,30)

| Case | XSLOPE | RS2 | Z-Soil | PLAXIS | GEO FEM | Slide2 LEM |
|---|---|---|---|---|---|---|
| 1 | 0.449 | 0.44 | 0.46 | 0.44 | 0.48 | 0.44 |
| 2 | 0.818 | 0.79 | 0.83 | 0.85 | 0.91 | 0.80 |
| 3 | 0.687 | 0.69 | 0.71 | 0.71 | 0.73 | 0.69 |
| 4 | 1.107 | 1.11 | 1.14 | 1.17 | 1.18 | 1.10 |
| 5 | 0.944 | 0.96 | 0.98 | 0.97 | 1.03 | 0.95 |
| 6 | 1.411 | 1.42 | 1.52 | 1.45 | 1.54 | 1.40 |

**H = 14 m (#58):** same six material cases as #57

| Case | XSLOPE | RS2 | Z-Soil | PLAXIS | GEO FEM | Slide2 LEM |
|---|---|---|---|---|---|---|
| 1 | 0.342 | 0.33 | 0.34 | 0.35 | 0.35 | 0.34 |
| 2 | 0.606 | 0.59 | 0.61 | 0.59 | 0.63 | 0.60 |
| 3 | 0.523 | 0.52 | 0.54 | 0.53 | 0.59 | 0.53 |
| 4 | 0.833 | 0.83 | 0.84 | 0.82 | 0.86 | 0.84 |
| 5 | 0.667 | 0.72 | 0.75 | 0.74 | 0.73 | 0.73 |
| 6 | 1.057 | 1.06 | 1.07 | 1.06 | 1.10 | 1.08 |

Case 5 of the H = 14 m slope is the one outlier (−7.4% against a tight published
cluster; the same materials at H = 10.5 m agree within 1.6%) — it is reported here and
excluded from the regression locks pending an explanation. Each slope's locks bracket
its family (the weakest and strongest case).

### RS2-60: Generalized Hoek-Brown, homogeneous slope (Li et al. 2008) {#rs2-60}

**Input files:** [rs2_60a.xlsx](../files/rocscience/rs2_60a.xlsx) (β = 15°) ·
[rs2_60b.xlsx](../files/rocscience/rs2_60b.xlsx) (β = 30°) ·
[rs2_60c.xlsx](../files/rocscience/rs2_60c.xlsx) (β = 45°)

A homogeneous rock slope at three angles, after

> Li, A.J., Merifield, R.S., & Lyamin, A.V. (2008). "Stability charts for rock slopes based
> on the Hoek-Brown failure criterion." *International Journal of Rock Mechanics and Mining
> Sciences* 45(5), 689–700.

GSI = 70, $m_i$ = 15, $D$ = 0, γ = 23 kN/m³, ν = 0.3, $H$ = 1 m. This is the companion to the
[Hammah benchmark](#hoek-brown) at the opposite end of the criterion: GSI = 70 is a strong,
lightly-jointed rock mass ($a$ = 0.501, essentially the classical exponent), where Hammah's
GSI = 5 is a badly broken one ($a$ = 0.619).

The manual does not state σ<sub>ci</sub>. Li's Table 1 tabulates the *critical* ratio
σ<sub>ci</sub>/(γH) — the value at which collapse has just occurred — so σ<sub>ci</sub> is
recovered by multiplying it by γH, and **the verification target is FS = 1.0 by
construction** rather than an independently computed factor of safety:

| case | β | σ<sub>ci</sub>/(γH) | σ<sub>ci</sub> |
|---|---:|---:|---:|
| a | 15° | 0.026 | 0.598 kPa |
| b | 30° | 0.075 | 1.725 kPa |
| c | 45° | 0.176 | 4.048 kPa |

Those magnitudes look wrong for rock and are the trap in this problem. $H$ = 1 m makes
γH = 23 kPa, and the critical ratio is *less than one*, so σ<sub>ci</sub> is a fraction of
γH — sub-kPa to a few kPa. The problem is normalized: only the ratio matters, and a 1 m
slope in 0.6 kPa rock is the same problem as a 100 m slope in 60 kPa rock. Entering
σ<sub>ci</sub> in MPa, as Hoek-Brown convention would invite, overstates the strength a
thousandfold and the slope becomes trivially stable.

**Factors of safety:**

| case | Bishop | Spencer | Li (limit analysis) | Slide2 Spencer |
|---|---|---|---|---|
| a (β = 15°) | 1.009 | 1.009 | 1.0 | 1.011 |
| b (β = 30°) | 1.015 | 1.017 | 1.0 | 0.992 |
| c (β = 45°) | 1.003 | 1.008 | 1.0 | 1.035 |

*All three land within 1.7% of limiting equilibrium, confirming the Hoek-Brown
implementation at high GSI. SSRM is not locked on this problem.*

Li's Table 1 prints its last block as β = 10°, but the body text and the charts (Fig. 5 is
β = 15°; no β = 10° chart exists) both say 15° — a typographical error in the paper. RS2 read
it as 15° as well: its Slide2 value for case a (1.011) reproduces Li's own F for that row.

<!-- test: file=../files/rocscience/rs2_60a.xlsx, type=circular_search, method=spencer, expected_fs=1.009, num_slices=40, benchmark=RS2-60a -->
<!-- test: file=../files/rocscience/rs2_60b.xlsx, type=circular_search, method=spencer, expected_fs=1.017, num_slices=40, benchmark=RS2-60b -->
<!-- test: file=../files/rocscience/rs2_60c.xlsx, type=circular_search, method=spencer, expected_fs=1.008, num_slices=40, benchmark=RS2-60c -->

## Hoek-Brown verification (Hammah et al. 2005) {#hoek-brown}

The `hb` strength option is verified end-to-end — LEM *and* SSRM — against Example 1 of the
Rocscience method paper that introduced Hoek-Brown shear-strength reduction. It is the
low-GSI counterpart to [RS2-60](#rs2-60) above, which exercises the same criterion at
GSI = 70:

> Hammah, R.E., Yacoub, T.E., Corkum, B., & Curran, J.H. (2005). "The shear strength
> reduction method for the generalized Hoek-Brown criterion." *Proc. 40th U.S. Symposium
> on Rock Mechanics (ARMA/USRMS)*, Anchorage, Paper 05-810.

A 10 m high homogeneous rock slope at 45° in a very weak rock mass: $\sigma_{ci}$ = 30 MPa,
GSI = 5, $m_i$ = 2, $D$ = 0, $\gamma$ = 25 kN/m³, $E$ = 5000 MPa, $\nu$ = 0.3. This is a
demanding test of the criterion rather than of the geometry — at GSI = 5 the envelope is
strongly curved, and the exponent $a$ = 0.619 is far from the $a$ = 0.5 special case.

**Derived constants** reproduce the paper's Table 1 exactly:

| | XSLOPE | Hammah et al. |
|---|---|---|
| $m_b$ | 0.0672 | 0.067 |
| $s$ | 2.605e-5 | 2.5e-5 |
| $a$ | 0.6192 | 0.619 |

**Factors of safety:**

| Method | XSLOPE | Hammah et al. |
|---|---|---|
| Bishop simplified | 1.150 | 1.153 |
| Spencer | 1.152 | 1.152 |
| Janbu corrected | 1.144 | — |
| Morgenstern-Price | 1.148 | — |
| SSRM | 1.159 | 1.15 (generalized Hoek-Brown *and* equivalent Mohr-Coulomb) |

The paper dimensions only the slope itself (10 m high, 10 m run) and leaves the foundation
depth and lateral extents unstated. The answer does not depend on them: foundation depths of
2, 4, 6 and 10 m all return Bishop 1.150 / Spencer 1.152, because the critical mechanism
exits at the toe. SSRM converges on the published value from above as the mesh refines
(1.165 at 1.0 m, 1.158 at 0.6 m), and is quoted here at the tagged 0.9 m mesh.

Corps of Engineers (1.191) and Lowe & Karafiath (1.166) both converge on this slope. They
are the two methods that struggle on *strong* rock masses, where the instantaneous friction
angle at low confinement exceeds ~55°; at GSI = 5 the envelope is weak enough that they are
well behaved. See the note in the [LEM overview](../lem/overview.md#hoek-brown-strength).

<!-- test: file=../files/rocscience/hammah_hb1.xlsx, type=circular_search, num_slices=40, fs_bishop=1.150, fs_spencer=1.152, fs_janbu=1.144, fs_mprice=1.148, benchmark=HB-lem -->
<!-- test: file=../files/rocscience/hammah_hb1.xlsx, type=fem_ssrm, expected_fs=1.159, element_type=tri6, target_size=0.9, tolerance=0.01, f_min=0.8, f_max=1.6, max_iter=4000, benchmark=HB-ssrm -->
