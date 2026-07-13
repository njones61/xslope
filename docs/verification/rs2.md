# Rocscience RS2 (SSR) Corpus

This page tracks the RS2 Slope Stability Verification Manual (Rocscience, Parts I–III,
68 problems) the way the [Slide2 corpus](rocscience.md) tracks its manual — but for the
**shear strength reduction (SSR)** method against XSLOPE's FEM/SSRM solver rather than
limit equilibrium. The long-standing SSRM anchors (Griffiths & Lane 1999 and the feature
samples) live on the [SSRM benchmarks page](ssrm.md).

Where a problem shares its geometry with a built Slide2 problem, the SSR analysis runs on
the **same corpus input file** — the extraction is already validated there. SSR results
use the Griffiths elastic convention (E = 10⁵ kPa, ν = 0.3, ψ = 0; SSR factors are
insensitive to these), stored as inert values in the shared files. SSR factors are quoted
at the tagged mesh size: SSRM drifts a percent or two with refinement, so tolerances are
honest rather than tight. A recurring pattern on this corpus: fine-mesh SSRM finds
shallow-skin mechanisms that published coarse-mesh SSR analyses miss or deliberately
suppress with "can't fail" elastic regions (#22, #23) — where the published value
depends on such an artifice rather than the mechanics, the problem is recorded as not
lockable rather than tuned to match.

<!-- test: file=../lem/files/xslope_acads_simple.xlsx, type=fem_ssrm, expected_fs=0.986, element_type=tri6, target_size=0.9, tolerance=0.01, f_min=0.7, f_max=1.3, max_iter=4000, benchmark=RS2-1 -->
<!-- test: file=../files/rocscience/vp003.xlsx, type=fem_ssrm, expected_fs=1.375, element_type=tri6, target_size=0.9, tolerance=0.01, f_min=1.0, f_max=1.7, max_iter=4000, benchmark=RS2-2 -->
<!-- test: file=../files/rocscience/vp004.xlsx, type=fem_ssrm, expected_fs=0.948, element_type=tri6, target_size=0.9, tolerance=0.01, f_min=0.7, f_max=1.3, max_iter=4000, benchmark=RS2-3 -->
<!-- test: file=../files/rocscience/vp005.xlsx, type=fem_ssrm, expected_fs=1.909, element_type=tri6, target_size=6.5, tolerance=0.01, f_min=1.5, f_max=2.3, max_iter=4000, benchmark=RS2-4 -->
<!-- test: file=../lem/files/xslope_acads_weak_layer.xlsx, type=fem_ssrm, expected_fs=1.286, element_type=tri6, target_size=2.0, tolerance=0.01, f_min=0.9, f_max=1.6, max_iter=4000, benchmark=RS2-5 -->
<!-- test: file=../lem/files/xslope_arai_tagyo.xlsx, type=fem_ssrm, expected_fs=1.427, element_type=tri6, target_size=2.2, tolerance=0.02, f_min=1.2, f_max=1.7, max_iter=4000, benchmark=RS2-10 -->
<!-- test: file=../files/rocscience/vp015.xlsx, type=fem_ssrm, expected_fs=0.413, element_type=tri6, target_size=1.9, tolerance=0.02, f_min=0.25, f_max=0.65, max_iter=4000, benchmark=RS2-11 -->
<!-- test: file=../files/rocscience/vp016.xlsx, type=fem_ssrm, expected_fs=1.115, element_type=tri6, target_size=1.3, tolerance=0.02, f_min=0.9, f_max=1.45, max_iter=4000, benchmark=RS2-12 -->
<!-- test: file=../files/rocscience/vp017.xlsx, type=fem_ssrm, expected_fs=1.384, element_type=tri6, target_size=0.5, tolerance=0.02, f_min=1.1, f_max=1.65, max_iter=4000, benchmark=RS2-13 -->
<!-- test: file=../files/rocscience/vp019.xlsx, type=fem_ssrm, expected_fs=1.386, element_type=tri6, target_size=4.33, tolerance=0.02, f_min=1.1, f_max=1.7, max_iter=4000, benchmark=RS2-15 -->
<!-- test: file=../files/rocscience/vp021a.xlsx, type=fem_ssrm, expected_fs=2.018, element_type=tri6, target_size=3.0, tolerance=0.02, f_min=1.6, f_max=2.5, max_iter=4000, benchmark=RS2-17 -->
<!-- test: file=../files/rocscience/vp022a.xlsx, type=fem_ssrm, expected_fs=1.336, element_type=tri6, target_size=3.0, tolerance=0.02, f_min=1.0, f_max=1.7, max_iter=4000, benchmark=RS2-18 -->
<!-- test: file=../files/rocscience/vp024.xlsx, type=fem_ssrm, expected_fs=1.507, element_type=tri6, target_size=1.0, tolerance=0.02, f_min=1.1, f_max=1.8, max_iter=4000, benchmark=RS2-19 -->
<!-- test: file=../files/rocscience/vp025.xlsx, type=fem_ssrm, expected_fs=1.003, element_type=tri6, target_size=0.8, tolerance=0.01, f_min=0.5, f_max=1.6, max_iter=4000, benchmark=RS2-20 -->
<!-- test: file=../files/rocscience/vp032a.xlsx, type=fem_ssrm, expected_fs=1.196, element_type=tri6, target_size=1.5, tolerance=0.02, f_min=0.9, f_max=1.6, max_iter=4000, benchmark=RS2-24a -->
<!-- test: file=../files/rocscience/vp032c.xlsx, type=fem_ssrm, expected_fs=0.991, element_type=tri6, target_size=2.2, tolerance=0.02, f_min=0.7, f_max=1.4, max_iter=4000, benchmark=RS2-24b -->
<!-- test: file=../files/rocscience/vp045a.xlsx, type=fem_ssrm, expected_fs=2.852, element_type=tri6, target_size=1.0, tolerance=0.02, f_min=2.3, f_max=3.4, max_iter=4000, benchmark=RS2-32 -->

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
| 1 | Simple slope stability assessment — **built**: SSRM 0.986 vs RS2 SSR 0.99, Slide Bishop 0.987, ACADS referee 1.00 (on [xslope_acads_simple.xlsx](../lem/files/xslope_acads_simple.xlsx); FS reads 0.967 at half the element size — SSR values are quoted at the tagged mesh) | VP1 (ACADS 1a) |
| 2 | Non-homogeneous slope — **built**: SSRM 1.375 vs RS2 1.36, Slide Spencer 1.375, referee 1.39 (on [vp003.xlsx](../files/rocscience/vp003.xlsx)) | VP3 |
| 3 | Non-homogeneous slope with seismic load (0.15g) — **built**: SSRM 0.948 vs RS2 0.97, Slide Spencer 0.991, referee 1.00 (on [vp004.xlsx](../files/rocscience/vp004.xlsx); k is entered negative per the FEM sign convention — this is a left-facing slope, so the pseudo-static force acts in −x, while the LEM takes the magnitude and directs it from the failure surface) | VP4 |
| 4 | Dry Talbingo dam — **built**: SSRM 1.909 vs RS2 1.88, Slide 1.948, referee 1.95 (on [vp005.xlsx](../files/rocscience/vp005.xlsx)) | VP5 |
| 5 | Water table with weak seam — **built**: SSRM 1.286 (mesh-stable: 1.280 at 1.2 m) vs RS2 1.26, Slide Spencer 1.258, referee 1.24–1.27 (on [xslope_acads_weak_layer.xlsx](../lem/files/xslope_acads_weak_layer.xlsx); this problem is Slide VP7) | VP7 |
| 6 | Slope with load and pore pressure by water table | VP9/VP10 family |
| 7 | Pore pressure by digitized total head grid | VP11 family |
| 8 | Slope stability with a pore pressure grid | VP12 family |
| 9 | Pore pressure grid with two limit sets | VP13 family |
| 10 | Simple slope II (Arai & Tagyo ex. 1) — **built**: SSRM 1.43 vs RS2 SSR 1.40 (+2%), mesh-converged (1.428→1.434 over a 2.9× size change); LEM locks Bishop 1.404 / Spencer 1.401 vs Slide2 1.409 / 1.406 (on [xslope_arai_tagyo.xlsx](../lem/files/xslope_arai_tagyo.xlsx)) | VP14 (Arai & Tagyo 1) |
| 11 | Layered slope (Arai & Tagyo ex. 2) — **built**: SSRM 0.41 vs RS2 SSR 0.39 and Greco/Kim pattern-search 0.39–0.43; LEM locks 0.419–0.422 (on [vp015.xlsx](../files/rocscience/vp015.xlsx)) | VP15 |
| 12 | Simple slope + water table (Arai & Tagyo ex. 3) — **built**: SSRM 1.12 vs RS2 SSR 1.09 (+2.3%); LEM locks Bishop 1.112 / Spencer 1.113 (on [vp016.xlsx](../files/rocscience/vp016.xlsx)). The FEM piezo pore pressure uses the vertical-distance convention, consistent with the LEM slicer and the published analyses | VP16 |
| 13 | Simple slope III (Yamagami & Ueta) — **built**: SSRM 1.38 vs RS2 SSR 1.33 and Greco Spencer 1.33; LEM locks Bishop 1.342 / Spencer 1.340 vs Y&U 1.348 / 1.339 (on [vp017.xlsx](../files/rocscience/vp017.xlsx)) | VP17 |
| 14 | Simple slope, pore pressure by r<sub>u</sub> — *blocked*: RS2 SSR 0.98 vs Slide2 Spencer 1.01 / Baker 1.02; XSLOPE's FEM has no r<sub>u</sub> pore-pressure option (the LEM does — vp018.xlsx locks Spencer 1.033). This problem is Slide2 VP18, not VP21 | VP18 |
| 15 | Layered slope II (Greco ex. 4 / Yamagami & Ueta) — **built**: SSRM 1.39 vs RS2 SSR 1.39, Slide2 Spencer 1.398, Greco 1.40–1.42; mesh-converged (1.386→1.377 over a 1.7× size change) (on [vp019.xlsx](../files/rocscience/vp019.xlsx)) | VP19 |
| 16 | Layered slope and water table with weak seam (Greco ex. 5 / Chen & Shao) — *blocked*: RS2 SSR 1.02 vs Slide2 Spencer 1.093 circular / 1.007 noncircular; XSLOPE's SSRM does not reach equilibrium on this polygon-zone + sloping-water-table model (under investigation; the LEM locks 1.086–1.091 on the same file) | VP20 |
| 17 | Slope with three pore pressure conditions (Fredlund & Krahn) — **built** (dry case): SSRM 2.02 vs RS2 SSR 2.0, Slide2 M-P 2.075, F&K 2.076 (on [vp021a.xlsx](../files/rocscience/vp021a.xlsx)). The r<sub>u</sub> case awaits FEM r<sub>u</sub> support; the water-table case awaits the VP21 case-3 input file | VP21 |
| 18 | Three pore pressure conditions and a weak seam (Fredlund & Krahn) — **built** (dry case): SSRM 1.34 vs RS2 SSR 1.34, Slide2 Bishop 1.382 (on [vp022a.xlsx](../files/rocscience/vp022a.xlsx)). Same r<sub>u</sub>/water-table status as #17 | VP22 |
| 19 | Undrained layered slope (Low 1989) — **built** (caveat): SSRM 1.51 at the tagged mesh (1.48 at 0.6 m) vs RS2 SSR 1.41, Slide2 LEM 1.439, Low 1.44 — the SSR values straddle the LEM from opposite sides on this φ = 0 slope, and the XSLOPE factor drifts −2% with refinement; quoted at the tagged mesh per the page convention. This problem is Slide2 VP24 | VP24 |
| 20 | Slope with vertical load (Prandtl's wedge) — **built**: SSRM 1.00 vs Prandtl theory 1.0 and RS2 SSR 1.0 (mesh pair 1.011→1.003); Slide2 Spencer reads 1.051 on the specified surface (on [vp025.xlsx](../files/rocscience/vp025.xlsx)) | VP25 |
| 21 | Bearing capacity test prism (Prandtl II) — planned: RS2 SSR 1.01 vs theory 1.0; the same setup as #20 on the VP26 prism, expected equally clean | VP26 |
| 22 | Layered slope with undulating bedrock — *blocked*: RS2 SSR 1.52. Two gaps: the FEM fixes displacements only along a flat base (an undulating bedrock bottom never reaches equilibrium), and the problem specifies phreatic-inclination-corrected pore pressures, which the FEM does not yet apply (the LEM does, via Type=phreatic — vp027's LEM locks stand) | VP27 |
| 23 | Underwater slope with linearly varying cohesion — *no lock possible*: RS2's published SSR (1.12) depends on a "can't fail" elastic region whose boundary its text and figure draw differently; stand-ins for the two readings give 0.87 and 0.92, and without the patch the true SSR minimum is the shallow skin above el. −20 (FS 0.21) that the artifice suppresses. The comparison would test where the patch is drawn, not the mechanics — this slope's anchor remains the LEM lock (VP29, Spencer 1.145 on Duncan's surface) | VP29 |
| 24 | Layered slope with geosynthetic reinforcement — **built**: SSRM 1.196 (H=7 case; mesh-stable 1.202 at 2.2 m) and 0.991 (H=8.75) vs RS2 SSR 1.15 / 0.95 (+4%, within the corpus SSR scatter; geotextile as an FEM truss with EA = 2×10³ kN/m stated as convention). RS2's fully labeled figures also supplied the geometry that unlocked Slide2's VP32 — LEM locks on the three printed circles live there (on [vp032a](../files/rocscience/vp032a.xlsx)/[c](../files/rocscience/vp032c.xlsx)) | VP32 |
| 25 | Syncrude tailings dyke, multiple phreatic surfaces | VP33 family |
| 26 | Clarence Cannon dam | VP34 family |
| 27 | Homogeneous slope, pore pressure by r<sub>u</sub> | VP21 |
| 28 | FE analysis with groundwater and stress | VP38 family |
| 29 | (heading mismatch — body is the Tandjiria reinforced embankment) — *blocked*: the FEM has no tension-crack input, and the problem's LEM values depend on the crack (without it SSRM reads 1.053 vs the crack-included LEM lock 0.968); waits on the FEM tension-crack gap | VP39 |
| 30 | (heading mismatch — body is the Perry power-curve slope) — *blocked*: RS2 SRF 0.91 vs Slide2 Janbu 0.944 / Perry 0.98; the FEM has no power-curve strength (one gap blocks #30/#31/#33/#34) | VP40 |
| 31 | Mohr-Coulomb vs power curve — *blocked* (power half): FEM power-curve strength missing; the MC halves live under #32 | VP44/VP45 |
| 32 | (heading mismatch — body is Baker's example 2) — **built** (MC half): SSRM 2.852 (2.861 at 1.6 m) vs RS2 SRF 2.83 (+0.8%); Baker's own MC 2.66, Slide2 Spencer 2.76 (on [vp045a.xlsx](../files/rocscience/vp045a.xlsx)). Power half blocked with #30 (RS2 SRF 2.74) | VP45 |
| 33–34 | Mohr-Coulomb vs power curve (II, III) — *blocked*: FEM power-curve strength missing (same gap as #30/#31) | VP61 (Baker 2003) |

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
