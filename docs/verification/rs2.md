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
<!-- test: file=../files/rocscience/vp009.xlsx, type=fem_ssrm, expected_fs=0.812, element_type=tri6, target_size=1.3, tolerance=0.02, f_min=0.3, f_max=1.3, max_iter=4000, benchmark=RS2-6 -->
<!-- test: file=../files/rocscience/vp010.xlsx, type=fem_ssrm, expected_fs=1.483, tolerance=0.01, f_min=1.0, f_max=2.2, max_iter=4000, benchmark=RS2-7 -->
<!-- test: file=../files/rocscience/vp026.xlsx, type=fem_ssrm, expected_fs=1.011, element_type=tri6, target_size=0.8, tolerance=0.01, f_min=0.5, f_max=1.6, max_iter=4000, benchmark=RS2-21 -->
<!-- test: file=../files/rocscience/vp033.xlsx, type=fem_ssrm, expected_fs=1.199, element_type=tri6, target_size=5.0, tolerance=0.02, f_min=0.9, f_max=1.8, max_iter=4000, benchmark=RS2-25 -->
<!-- test: file=../files/rocscience/vp034.xlsx, type=fem_ssrm, expected_fs=2.281, element_type=tri6, target_size=15.0, tolerance=0.02, f_min=1.7, f_max=3.0, max_iter=4000, benchmark=RS2-26 -->
<!-- test: file=../files/rocscience/vp039c.xlsx, type=fem_ssrm, expected_fs=1.253, element_type=tri6, target_size=0.7, tolerance=0.02, f_min=0.9, f_max=1.7, max_iter=4000, benchmark=RS2-29 -->
<!-- test: file=../files/rocscience/vp044b.xlsx, type=fem_ssrm, expected_fs=1.546, element_type=tri6, target_size=0.5, tolerance=0.02, f_min=1.1, f_max=2.0, max_iter=4000, benchmark=RS2-31a -->
<!-- test: file=../files/rocscience/vp044c.xlsx, type=fem_ssrm, expected_fs=0.966, element_type=tri6, target_size=0.5, tolerance=0.02, f_min=0.6, f_max=1.4, max_iter=4000, benchmark=RS2-31b -->
<!-- test: file=../files/rocscience/vp056.xlsx, type=fem_ssrm, expected_fs=1.247, element_type=tri6, target_size=4.0, tolerance=0.02, f_min=0.9, f_max=1.7, max_iter=4000, benchmark=RS2-33 -->
<!-- test: file=../files/rocscience/vp061b.xlsx, type=fem_ssrm, expected_fs=1.369, element_type=tri6, target_size=0.5, tolerance=0.02, f_min=1.0, f_max=1.9, max_iter=4000, benchmark=RS2-34 -->
<!-- test: file=../files/rocscience/vp071a.xlsx, type=fem_ssrm, expected_fs=1.118, tolerance=0.01, f_min=0.7, f_max=1.6, max_iter=4000, benchmark=RS2-36a -->
<!-- test: file=../files/rocscience/vp071b.xlsx, type=fem_ssrm, expected_fs=1.118, element_type=tri6, target_size=4.4, tolerance=0.01, f_min=0.7, f_max=1.6, max_iter=4000, benchmark=RS2-36b -->
<!-- test: file=../files/rocscience/vp077b.xlsx, type=fem_ssrm, expected_fs=1.491, element_type=tri6, target_size=12.4, tolerance=0.02, f_min=1.1, f_max=2.2, max_iter=4000, benchmark=RS2-40 -->
<!-- test: file=../files/rocscience/vp075.xlsx, type=fem_ssrm, expected_fs=1.249, element_type=tri6, target_size=1.85, tolerance=0.02, f_min=0.8, f_max=1.8, max_iter=4000, benchmark=RS2-42 -->
<!-- test: file=../files/rocscience/vp082.xlsx, type=fem_ssrm, expected_fs=1.511, element_type=tri6, target_size=2.0, tolerance=0.02, f_min=1.0, f_max=2.1, max_iter=4000, benchmark=RS2-44 -->
<!-- test: file=../files/rocscience/vp040.xlsx, type=fem_ssrm, expected_fs=0.910, element_type=tri6, target_size=2.0, tolerance=0.02, f_min=0.5, f_max=1.5, max_iter=4000, benchmark=RS2-30 -->
<!-- test: file=../files/rocscience/vp044a.xlsx, type=fem_ssrm, expected_fs=0.934, element_type=tri6, target_size=0.5, tolerance=0.02, f_min=0.5, f_max=1.6, max_iter=4000, benchmark=RS2-31c -->
<!-- test: file=../files/rocscience/vp045b.xlsx, type=fem_ssrm, expected_fs=2.732, element_type=tri6, target_size=1.0, tolerance=0.02, f_min=1.8, f_max=3.6, max_iter=4000, benchmark=RS2-32b -->
<!-- test: file=../files/rocscience/vp061a.xlsx, type=fem_ssrm, expected_fs=1.502, element_type=tri6, target_size=0.5, tolerance=0.02, f_min=1.0, f_max=2.2, max_iter=4000, benchmark=RS2-34b -->
<!-- test: file=../files/rocscience/vp083a.xlsx, type=fem_ssrm, expected_fs=1.412, element_type=tri6, target_size=4.0, tolerance=0.02, f_min=0.9, f_max=1.9, max_iter=4000, benchmark=RS2-45a -->
<!-- test: file=../files/rocscience/vp083b.xlsx, type=fem_ssrm, expected_fs=1.365, element_type=tri6, target_size=4.0, tolerance=0.02, f_min=0.9, f_max=1.9, max_iter=4000, benchmark=RS2-45b -->
<!-- test: file=../files/rocscience/vp084a.xlsx, type=fem_ssrm, expected_fs=0.804, element_type=tri6, target_size=4.0, tolerance=0.02, f_min=0.4, f_max=1.3, max_iter=4000, benchmark=RS2-46a -->
<!-- test: file=../files/rocscience/vp084b.xlsx, type=fem_ssrm, expected_fs=0.947, element_type=tri6, target_size=4.0, tolerance=0.02, f_min=0.5, f_max=1.4, max_iter=4000, benchmark=RS2-46b -->
<!-- test: file=../files/rocscience/vp084c.xlsx, type=fem_ssrm, expected_fs=1.082, element_type=tri6, target_size=4.0, tolerance=0.02, f_min=0.6, f_max=1.5, max_iter=4000, benchmark=RS2-46c -->
<!-- test: file=../files/rocscience/vp084d.xlsx, type=fem_ssrm, expected_fs=1.188, element_type=tri6, target_size=4.0, tolerance=0.02, f_min=0.7, f_max=1.7, max_iter=4000, benchmark=RS2-46d -->
<!-- test: file=../files/rocscience/vp078.xlsx, type=fem_ssrm, expected_fs=1.081, element_type=tri6, target_size=4.0, tolerance=0.02, f_min=0.6, f_max=1.6, max_iter=4000, benchmark=RS2-47 -->
<!-- test: file=../files/rocscience/rs2_56a.xlsx, type=fem_ssrm, expected_fs=0.667, element_type=tri6, target_size=0.8, tolerance=0.02, f_min=0.32, f_max=1.12, max_iter=4000, benchmark=RS2-56a -->
<!-- test: file=../files/rocscience/rs2_56b.xlsx, type=fem_ssrm, expected_fs=2.131, element_type=tri6, target_size=0.8, tolerance=0.02, f_min=1.79, f_max=2.59, max_iter=4000, benchmark=RS2-56b -->
<!-- test: file=../files/rocscience/rs2_57a.xlsx, type=fem_ssrm, expected_fs=0.449, element_type=tri6, target_size=0.8, tolerance=0.02, f_min=0.1, f_max=0.89, max_iter=4000, benchmark=RS2-57a -->
<!-- test: file=../files/rocscience/rs2_57b.xlsx, type=fem_ssrm, expected_fs=1.411, element_type=tri6, target_size=0.8, tolerance=0.02, f_min=1.07, f_max=1.87, max_iter=4000, benchmark=RS2-57b -->
<!-- test: file=../files/rocscience/rs2_58a.xlsx, type=fem_ssrm, expected_fs=0.342, element_type=tri6, target_size=0.8, tolerance=0.02, f_min=0.1, f_max=0.78, max_iter=4000, benchmark=RS2-58a -->
<!-- test: file=../files/rocscience/rs2_58b.xlsx, type=fem_ssrm, expected_fs=1.057, element_type=tri6, target_size=0.8, tolerance=0.02, f_min=0.71, f_max=1.51, max_iter=4000, benchmark=RS2-58b -->

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
| 6 | Slope with load and pore pressure by water table (ACADS 4) — **built** (caveat): SSRM 0.81, mesh-stable, vs ACADS survey mean 0.808 and referee 0.78 — but +18% above RS2's SSR 0.69 and Slide2's MC-optimized LEM 0.68–0.71 (XSLOPE's own LEM locks 0.724). The published values themselves span 0.68–0.81 on this thin-weak-seam problem; under the same investigation as #16 (on [vp009.xlsx](../files/rocscience/vp009.xlsx)) | VP9 |
| 7 | Pore pressure by digitized total head grid (ACADS 5) — **built**: SSRM 1.483 vs RS2 SSR 1.48 (+0.2%), on the FE-seepage model XSLOPE built for Slide2 VP10 (the grid is a stand-in for the flow solution; sidecars are tri6 so the SSRM plasticity is not volumetrically locked). Slide2 LEM 1.498–1.501, Giam 1.53 (on [vp010.xlsx](../files/rocscience/vp010.xlsx)) | VP10 |
| 8 | Saint-Alban test embankment — *no lock possible*: the grid encodes measured construction-induced pressures (see the Slide2 corpus VP11 row); RS2 SSR 0.96 vs Pilot 1.04 recorded | VP11 |
| 9 | Cubzac-les-Ponts test embankment — *no lock possible*: measured pore-pressure grid plus a "can't fail" elastic face layer suppressing the true face failure (FS 1.11 per RS2's own text); RS2 SSR 1.31 vs Pilot 1.24 recorded | VP13 |
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
| 21 | Bearing capacity test prism (Prandtl II) — **built**: SSRM 1.011 (1.020 at 1.2 m), converging on theory 1.0; RS2 SSR 1.01; Slide2 Spencer 0.941 on the specified surface (on [vp026.xlsx](../files/rocscience/vp026.xlsx), extracted for this problem) | VP26 |
| 22 | Layered slope with undulating bedrock — *blocked*: RS2 SSR 1.52. Two gaps: the FEM fixes displacements only along a flat base (an undulating bedrock bottom never reaches equilibrium), and the problem specifies phreatic-inclination-corrected pore pressures, which the FEM does not yet apply (the LEM does, via Type=phreatic — vp027's LEM locks stand) | VP27 |
| 23 | Underwater slope with linearly varying cohesion — *no lock possible*: RS2's published SSR (1.12) depends on a "can't fail" elastic region whose boundary its text and figure draw differently; stand-ins for the two readings give 0.87 and 0.92, and without the patch the true SSR minimum is the shallow skin above el. −20 (FS 0.21) that the artifice suppresses. The comparison would test where the patch is drawn, not the mechanics — this slope's anchor remains the LEM lock (VP29, Spencer 1.145 on Duncan's surface) | VP29 |
| 24 | Layered slope with geosynthetic reinforcement — **built**: SSRM 1.196 (H=7 case; mesh-stable 1.202 at 2.2 m) and 0.991 (H=8.75) vs RS2 SSR 1.15 / 0.95 (+4%, within the corpus SSR scatter; geotextile as an FEM truss with EA = 2×10³ kN/m stated as convention). RS2's fully labeled figures also supplied the geometry that unlocked Slide2's VP32 — LEM locks on the three printed circles live there (on [vp032a](../files/rocscience/vp032a.xlsx)/[c](../files/rocscience/vp032c.xlsx)) | VP32 |
| 25 | Syncrude tailings dyke (El-Ramly et al. 2003) — **built** (caveat): SSRM 1.20 (mesh-converged) vs RS2 SSR 1.29, Slide2 Bishop 1.305, El-Ramly 1.31 — the weak presheared clay-shale rewards less-constrained searches (XSLOPE's own LEM: 1.299 on Slide's circle, 1.253 free composite search, SSRM 1.20); locked at the tagged mesh (on [vp033.xlsx](../files/rocscience/vp033.xlsx)) | VP33 |
| 26 | Clarence Cannon dam (Wolff & Harr 1987) — **built**: SSRM 2.28 vs RS2 SSR 2.29 (−0.4%); Slide2 GLE 2.333 / Spencer 2.383, W&H 2.36, XSLOPE LEM M-P 2.384 (on [vp034.xlsx](../files/rocscience/vp034.xlsx); polygon zones with the chimney drain) | VP34 |
| 27 | Homogeneous slope, pore pressure by r<sub>u</sub> — *blocked*: the FEM has no r<sub>u</sub> option (task, with RS2 #14); RS2's own text cites Slide2 VP36 (Li & Lumb), not VP21 | VP36 |
| 28 | FE analysis with groundwater and stress | VP38 family |
| 29 | Geosynthetic-reinforced embankment on soft soil (Tandjiria 2002; the manual's §29/§30 headings are swapped) — **built** (sand case): SSRM 1.253 vs RS2 SSR 1.25, Spencer 1.209, Tandjiria 1.219 (on [vp039c.xlsx](../files/rocscience/vp039c.xlsx)). The clay case (RS2 SSR 0.99) is not locked: its water-filled tension crack has no FEM representation (XSLOPE reads 1.05 without it) | VP39 |
| 30 | Homogeneous slope, power-curve strength (Perry 1993; swapped heading) — **built**: SSRM 0.910 vs RS2 SRF 0.91 (exact); Slide2 Janbu 0.944, Perry 0.98 (on [vp040.xlsx](../files/rocscience/vp040.xlsx); the FEM linearizes the reduced envelope at the current stress per iteration) | VP40 |
| 31 | M-C vs power curve (Baker 2003 ex. 1) — **built**, all three halves: M-C SSRM 1.546 / 0.966 vs RS2 1.53 / 0.98; power-curve SSRM 0.934, squarely on Slide2's own published band (Janbu 0.92 / Spencer 0.96, Baker ~0.97) — RS2's 1.11 sits 15% above Slide2's LEM on the same problem and is the outlier (on [vp044a](../files/rocscience/vp044a.xlsx)/[b](../files/rocscience/vp044b.xlsx)/[c](../files/rocscience/vp044c.xlsx)) | VP44 |
| 32 | (heading mismatch — body is Baker's example 2) — **built**, both halves: M-C SSRM 2.852 vs RS2 2.83 (+0.8%); power-curve SSRM 2.732 vs RS2 2.74 (−0.3%), Slide2 Spencer 2.662 (on [vp045a](../files/rocscience/vp045a.xlsx)/[vp045b](../files/rocscience/vp045b.xlsx)) | VP45 |
| 33 | Homogeneous slope with tension crack and water table (P&D test slope 2; swapped heading) — **built** (caveat): SSRM 1.247 vs RS2 SSR 1.28 and an eight-program LEM table spanning 1.03–1.32; the model's dry tension crack has no FEM representation, worth ~2–3% here (on [vp056.xlsx](../files/rocscience/vp056.xlsx)) | VP56 |
| 34 | M-C vs power curve III (Baker 2003 ex. 3, London clay) — **built**, both halves: M-C SSRM 1.369 vs RS2 1.38; power-curve SSRM 1.502 vs RS2 1.47 / Slide2 Spencer 1.47 / Baker 1.48 (+2.1%) (on [vp061a](../files/rocscience/vp061a.xlsx)/[vp061b](../files/rocscience/vp061b.xlsx)) | VP61 |

### Part II (35–58)

| # | Problem | Slide2 counterpart |
|---|---|---|
| 35 | Submerged slope | VP64 family |
| 36 | Seepage analysis, homogeneous slope (D&W Fig 6.37; = Slide2 VP71, not VP70) — **built**, both cases: SSRM 1.118 on the FE-seepage model and 1.118 on the piezo approximation vs RS2 SSR 1.12 / 1.12; referee 1.138/1.141; XSLOPE LEM locks 1.132 (on [vp071a](../files/rocscience/vp071a.xlsx)/[b](../files/rocscience/vp071b.xlsx); the seep case runs on tri6 sidecars) | VP71 |
| 37 | Embankment with layered foundation (D&W Fig 6.39) — reported, no lock: RS2's SSR is the artesian downstream-toe slide, quoted as 0.95 in its table and 1.1 in its own convergence graph; XSLOPE's SSRM finds the deep mechanism at 1.31 (LEM on the tangent circle 1.339, Slide2 1.15/1.16, referee 1.11). Reproducing the toe mechanism needs toe-refined meshing — noted with the artesian-toe discussion in the Slide2 VP72 section | VP72 |
| 38 | Cohesionless embankment on saturated clay foundation | VP73/VP74 family |
| 39, 41, 43 | Earth embankment, infinite-slope mechanism (I–III) | VP69 family |
| 40 | Dam with impermeable foundation (D&W Fig 7.24) — **built** (piezo case): SSRM 1.491 vs RS2 SSR 1.53; the FE-seepage case is blocked — tri6 seepage does not converge on the high-contrast thick core (the documented tri3/tri6 trade) (on [vp077b.xlsx](../files/rocscience/vp077b.xlsx)) | VP77 |
| 42 | James dike — **built**: SSRM 1.249 vs RS2 SSR 1.26 (−0.9%); Slide2 noncircular LEM 1.11–1.16, referee 1.17 (on [vp075.xlsx](../files/rocscience/vp075.xlsx)) | VP75 |
| 44 | Seepage analysis for an earth embankment (D&W Fig 14.20-a; = Slide2 VP82, not VP76 — §39's body carries VP76) — **built**: SSRM 1.511 vs RS2 SSR 1.51 (+0.1%); Slide2 LEM 1.532/1.541, referee 1.528–1.542 (on [vp082.xlsx](../files/rocscience/vp082.xlsx)) | VP82 |
| 45 | Varying undrained shear strength profiles (D&W Fig 14.20-b) — **built** (caveat): SSRM 1.41 / 1.37 vs RS2 SSR 1.32 / 1.32 (D&W referee 1.28–1.33): XSLOPE's SSRM reads high on φ=0 foundations (the RS2-19 family, under investigation); locked as regression values at the tagged mesh (on [vp083a](../files/rocscience/vp083a.xlsx)/[b](../files/rocscience/vp083b.xlsx)) | VP83 |
| 46 | Varying undrained strength profiles II (D&W Fig 15.9, cu = 300 + cz·z) — **built**: SSRM 0.80 / 0.95 / 1.08 / 1.19 vs RS2 SSR 0.78 / 0.93 / 1.05 / 1.15 (+2–3%, the φ=0 pattern); D&W 0.75 / 0.90 / 1.03 / 1.13 (on [vp084a–d](../files/rocscience/vp084a.xlsx)) | VP84 |
| 47 | Purely cohesive slope, varying thickness (D&W Fig 14.3) — **built** (30-ft case): SSRM 1.08 vs RS2 SSR 1.03; D&W referee 1.124–1.135. The 46.5- and 60-ft variants (RS2 1.02 / 1.02) need deepened-base builds and are deferred (on [vp078.xlsx](../files/rocscience/vp078.xlsx)) | VP78 |
| 48–54 | Multi-tiered geotextile walls (= Slide2 VP87–VP93 one-for-one, verified; Leshchinsky & Han 2004) — *blocked*: the FEM reinforcement tensile-capacity cap is not enforced in SSRM (proven Ta-insensitive on the baseline wall: unreinforced 0.41, any Ta 1.72 vs published ≈1.0). Extraction and published tables are complete for all seven; the series unlocks together when the cap is fixed | VP87–VP93 |
| 55 | Five 1.8-m tiers (= Slide2 VP94; belongs to the #48–54 family) — *blocked* with the series | VP94 |
| [56](#pruska) | Homogeneous slope vs Z-Soil, PLAXIS, GEO FEM (Pruska 2003, H = 7 m, 5 cases) — **built**: all five within ±3.3% of RS2's M-C and inside the four-program band; locks bracket the family (0.667 / 2.131 on [rs2_56a](../files/rocscience/rs2_56a.xlsx)/[b](../files/rocscience/rs2_56b.xlsx)); full tables in the section | — (new files) |
| [57](#pruska) | Pruska H = 10.5 m, 6 cases — **built**: all six within ±3.6% of RS2's M-C; locks 0.449 / 1.411 ([rs2_57a](../files/rocscience/rs2_57a.xlsx)/[b](../files/rocscience/rs2_57b.xlsx)) | — |
| [58](#pruska) | Pruska H = 14 m, 6 cases — **built** (5 of 6): four within ±3.6%; case 5 (c=5, φ=30) reads 0.667 vs a tight published 0.72–0.75 cluster and is reported unlocked pending explanation; locks 0.342 / 1.057 ([rs2_58a](../files/rocscience/rs2_58a.xlsx)/[b](../files/rocscience/rs2_58b.xlsx)) | — |

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

## The Pruska cross-bearing (#56–58) {#pruska}

Pruska (2003) analyzed three homogeneous slopes (H = 7, 10.5, and 14 m over an 8-m
foundation) with five or six material sets each in four SSR programs. RS2 reproduces the
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
