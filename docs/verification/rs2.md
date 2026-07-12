# Rocscience RS2 (SSR) Corpus

!!! note "Status: stub — corpus not yet started"
    This page will track the RS2 Slope Stability Verification Manual (Rocscience,
    Parts I–III, 68 problems) the way the [Slide2 corpus](rocscience.md) tracks its
    manual — but for the **shear strength reduction (SSR)** method against XSLOPE's
    FEM/SSRM solver, rather than limit equilibrium. Nothing is built yet; the
    current SSRM anchors (Griffiths & Lane 1999 and the feature samples) live on the
    [SSRM benchmarks page](ssrm.md).

The RS2 manual is unusually cheap to build against: a large fraction of its problems are
**SSR renditions of the same problems as the Slide2 LEM manual**, so the geometry and
materials are already extracted, validated, and sitting in the corpus input files — often
the only new work per problem is an SSRM run and a tag. Problems 56–58 additionally carry
published FS values from **Z-Soil, PLAXIS, and GEO FEM**, giving multi-program SSR
cross-bearings.

## Problem inventory

### Part I (1–34)

| # | Problem | Slide2 counterpart |
|---|---|---|
| 1 | Simple slope stability assessment | VP1 (ACADS 1a) |
| 2 | Non-homogeneous slope | VP3 |
| 3 | Non-homogeneous slope with seismic load | VP4 |
| 4 | Dry Talbingo dam | VP5 |
| 5 | Water table with weak seam | VP9 |
| 6 | Slope with load and pore pressure by water table | VP9/VP10 family |
| 7 | Pore pressure by digitized total head grid | VP11 family |
| 8 | Slope stability with a pore pressure grid | VP12 family |
| 9 | Pore pressure grid with two limit sets | VP13 family |
| 10 | Simple slope stability assessment II | VP14 (Arai & Tagyo 1) |
| 11 | Layered slope stability assessment | VP15 |
| 12 | Simple slope with water table | VP16 |
| 13 | Simple slope III | VP17 |
| 14 | Simple slope, pore pressure by r<sub>u</sub> | VP21 |
| 15 | Layered slope II | VP19/VP20 family |
| 16 | Layered slope and water table with weak seam | VP20 |
| 17 | Slope with three pore pressure conditions | VP21 family |
| 18 | Three pore pressure conditions and a weak seam | VP22 |
| 19 | Undrained layered slope | VP23/VP24 |
| 20 | Slope with vertical load (Prandtl's wedge) | VP25/VP26 |
| 22 | Layered slope with undulating bedrock | VP27 |
| 23 | Underwater slope with linearly varying cohesion | VP29 family |
| 24 | Layered slope with geosynthetic reinforcement | VP30–32 family |
| 25 | Syncrude tailings dyke, multiple phreatic surfaces | VP33 family |
| 26 | Clarence Cannon dam | VP34 family |
| 27 | Homogeneous slope, pore pressure by r<sub>u</sub> | VP21 |
| 28 | FE analysis with groundwater and stress | VP38 family |
| 29 | Power-curve strength criterion | VP44 (Baker) |
| 30 | Geosynthetic-reinforced embankment on soft soil | VP39 |
| 31 | Mohr-Coulomb vs power curve | VP44/VP45 |
| 32 | Tension crack and water table | VP56/VP57 family |
| 33–34 | Mohr-Coulomb vs power curve (II, III) | VP61 (Baker 2003) |

### Part II (35–58)

| # | Problem | Slide2 counterpart |
|---|---|---|
| 35 | Submerged slope | VP64 family |
| 36 | Seepage analysis, homogeneous slope | VP70 |
| 37 | Seepage analysis, embankment with layered foundation | **VP72** |
| 38 | Cohesionless embankment on saturated clay foundation | VP73/VP74 family |
| 39, 41, 43 | Earth embankment, infinite-slope mechanism (I–III) | VP69 family |
| 40 | Seepage analysis, dam with impermeable foundation | **VP77** |
| 42 | Planned cross-section of James dike | **VP75** |
| 44 | Seepage analysis for an earth embankment | VP76 |
| 45–46 | Varying undrained shear strength profiles (I, II) | VP83/VP84 |
| 47 | Purely cohesive slope with varying thickness | VP78 family |
| 48–54 | Multi-tiered walls (tiers, fill, length, type, foundation, seepage, surcharge) | VP85–VP94 |
| 56–58 | Homogeneous slope vs Z-Soil, PLAXIS, GEO FEM (I–III) | — (multi-program SSR cross-bearing) |

### Part III (59–68)

| # | Problem | Source |
|---|---|---|
| 59 | Three-layered soil slope | Görög & Török (2007), vs Slide2 + PLAXIS |
| 60 | Generalized Hoek–Brown, homogeneous slope | needs Hoek–Brown strength (not in XSLOPE) |
| 61 | Local and global minima, homogeneous slope | |
| 62 | Three-layered slope with a soft band | |
| 63 | Homogeneous slope assessment | |
| 64 | Three homogeneous landslides | |
| 65 | Tailings dam | |
| 66 | Embankment basal stability | |
| 67 | Earth dam under steady & transient unsaturated seepage | transient — blocked on a transient solver |
| 68 | Seismically loaded slopes | |

## Methodology

Same discipline as the [Slide2 corpus](rocscience.md): geometry from the manuals'
coordinate-labeled figures (or reused directly from the Slide2 corpus input files where the
problem is shared), results locked into `run_tests.py` via `fem_ssrm` test tags — the tag
type and runner already exist and currently guard the Griffiths & Lane anchors. SSRM runs
are expensive (~1 min each), so this corpus will lean on coarse meshes with honest
tolerances, the same trade documented for the FEM reliability regression.
