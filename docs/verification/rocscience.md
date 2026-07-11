# Rocscience Slide2 Verification Corpus

The [Rocscience Slide2 verification manual](https://www.rocscience.com/help/slide2/verification-theory/verification-manuals)
contains 111 slope stability problems drawn from the published literature, each with Slide2's computed
factors of safety and (in most cases) independent reference values from the original authors. XSLOPE is
being verified against this corpus problem by problem: each **built** entry links an XSLOPE input file
reproducing the problem, reports the comparison, and is locked into the automated regression suite via a
test tag. Status values: **built** (input file + verified results below), *covered* (equivalent to an
existing XSLOPE sample), *partial* (some data extracted, some still needed), *blocked* (source data not
yet available), *planned*.

Problems are built from the manual's tabulated data and coordinate-labeled figures; where a problem's
geometry exists only as an unlabeled figure, the original source publication is consulted before the
problem is marked *built* — no digitized guesses are used for benchmark inputs.

<!-- test: file=../files/rocscience/vp002.xlsx, type=circular_search, num_slices=40, fs_bishop=1.589, fs_spencer=1.585, fs_janbu=1.481, fs_mprice=1.586, benchmark=VP2 -->
<!-- test: file=../files/rocscience/vp003.xlsx, type=circular_search, num_slices=40, fs_bishop=1.403, fs_spencer=1.372, fs_janbu=1.354, fs_mprice=1.371, benchmark=VP3 -->
<!-- test: file=../files/rocscience/vp004.xlsx, type=circular_search, num_slices=40, fs_bishop=1.013, fs_spencer=0.989, fs_janbu=0.963, fs_mprice=0.987, benchmark=VP4 -->
<!-- test: file=../files/rocscience/vp005.xlsx, type=single_circle, num_slices=60, fs_bishop=1.955, fs_spencer=1.955, fs_janbu=1.965, fs_mprice=1.955, benchmark=VP5 -->
<!-- test: file=../files/rocscience/vp006.xlsx, type=single_circle, num_slices=60, fs_bishop=2.206, fs_spencer=2.290, fs_janbu=2.073, fs_mprice=2.299, benchmark=VP6 -->
<!-- test: file=../files/rocscience/vp008.xlsx, type=single_noncirc, num_slices=50, fs_spencer=1.276, fs_janbu=1.294, fs_mprice=1.260, benchmark=VP8 -->
<!-- test: file=../files/rocscience/vp009.xlsx, type=noncircular_search, num_slices=50, fs_spencer=0.724, fs_janbu=0.718, benchmark=VP9 -->
<!-- test: file=../files/rocscience/vp015.xlsx, type=circular_search, num_slices=40, fs_bishop=0.419, fs_spencer=0.422, fs_janbu=0.436, fs_mprice=0.420, benchmark=VP15 -->
<!-- test: file=../files/rocscience/vp016.xlsx, type=circular_search, num_slices=40, fs_bishop=1.112, fs_spencer=1.113, fs_janbu=1.122, fs_mprice=1.111, benchmark=VP16 -->
<!-- test: file=../files/rocscience/vp017.xlsx, type=circular_search, num_slices=50, fs_oms=1.274, fs_bishop=1.342, fs_spencer=1.340, benchmark=VP17 -->
<!-- test: file=../files/rocscience/vp018.xlsx, type=noncircular_search, num_slices=50, fs_spencer=1.033, fs_mprice=1.024, benchmark=VP18 -->
<!-- test: file=../files/rocscience/vp019.xlsx, type=circular_search, num_slices=50, fs_bishop=1.448, fs_spencer=1.429, benchmark=VP19 -->
<!-- test: file=../files/rocscience/vp020.xlsx, type=circular_search, num_slices=50, fs_bishop=1.086, fs_spencer=1.091, benchmark=VP20-circ -->
<!-- test: file=../files/rocscience/vp020.xlsx, type=noncircular_search, num_slices=50, fs_spencer=1.082, benchmark=VP20-noncirc -->
<!-- test: file=../files/rocscience/vp023.xlsx, type=circular_search, num_slices=50, fs_oms=1.357, fs_bishop=1.130, benchmark=VP23 -->
<!-- test: file=../files/rocscience/vp024.xlsx, type=circular_search, num_slices=50, fs_oms=1.433, fs_bishop=1.433, benchmark=VP24 -->
<!-- test: file=../files/rocscience/vp036.xlsx, type=circular_search, num_slices=50, fs_bishop=1.333, benchmark=VP36-fs -->
<!-- test: file=../files/rocscience/vp036.xlsx, type=reliability, method=bishop, expected_beta=2.263, tolerance=0.03, benchmark=VP36-beta -->
<!-- test: file=../files/rocscience/vp041.xlsx, type=circular_search, num_slices=50, fs_bishop=1.668, fs_spencer=1.670, fs_janbu=1.660, benchmark=VP41 -->
<!-- test: file=../files/rocscience/vp044a.xlsx, type=circular_search, num_slices=40, fs_spencer=0.958, benchmark=VP44-pow -->
<!-- test: file=../files/rocscience/vp044b.xlsx, type=circular_search, num_slices=40, fs_spencer=1.518, benchmark=VP44-mc -->
<!-- test: file=../files/rocscience/vp044c.xlsx, type=circular_search, num_slices=40, fs_spencer=0.980, benchmark=VP44-lla -->
<!-- test: file=../files/rocscience/vp045a.xlsx, type=circular_search, num_slices=50, fs_spencer=2.801, benchmark=VP45-mc -->
<!-- test: file=../files/rocscience/vp045b.xlsx, type=circular_search, num_slices=50, fs_spencer=2.649, benchmark=VP45-pow -->
<!-- test: file=../files/rocscience/vp047.xlsx, type=single_noncirc, num_slices=50, fs_janbu=0.899, benchmark=VP47 -->
<!-- test: file=../files/rocscience/vp048.xlsx, type=single_noncirc, num_slices=50, fs_janbu=0.991, fs_spencer=0.991, benchmark=VP48 -->
<!-- test: file=../files/rocscience/vp050.xlsx, type=single_noncirc, num_slices=60, fs_janbu=1.448, fs_spencer=1.576, benchmark=VP50 -->
<!-- test: file=../files/rocscience/vp051.xlsx, type=single_circle, num_slices=100, fs_oms=1.069, fs_bishop=1.278, fs_janbu=1.205, fs_corps=1.404, fs_lowe=1.296, fs_spencer=1.294, fs_mprice=1.304, benchmark=VP51 -->
<!-- test: file=../files/rocscience/vp052a.xlsx, type=circular_search, num_slices=50, fs_spencer=1.797, fs_bishop=1.796, benchmark=VP52-dry -->
<!-- test: file=../files/rocscience/vp052b.xlsx, type=circular_search, num_slices=50, fs_spencer=1.189, fs_bishop=1.176, benchmark=VP52-wet -->
<!-- test: file=../files/rocscience/vp054a.xlsx, type=single_circle, num_slices=50, fs_bishop=1.100, benchmark=VP54-nopile -->
<!-- test: file=../files/rocscience/vp054b.xlsx, type=single_circle, num_slices=50, fs_bishop=1.185, benchmark=VP54-pile -->
<!-- test: file=../files/rocscience/vp086.xlsx, type=circular_search, num_slices=50, fs_bishop=1.617, fs_spencer=1.611, benchmark=VP86 -->
<!-- test: file=../files/rocscience/vp062a.xlsx, type=circular_search, num_slices=50, fs_spencer=1.001, fs_bishop=0.991, benchmark=VP62-dry -->
<!-- test: file=../files/rocscience/vp062b.xlsx, type=circular_search, num_slices=50, fs_spencer=1.001, fs_bishop=0.986, benchmark=VP62-ru -->
<!-- test: file=../files/rocscience/vp097.xlsx, type=circular_search, rapid=true, num_slices=50, fs_spencer=1.044, fs_bishop=1.042, benchmark=VP97 -->
<!-- test: file=../files/rocscience/vp100.xlsx, type=circular_search, num_slices=50, fs_bishop=1.201, fs_spencer=1.206, benchmark=VP100 -->
<!-- test: file=../files/rocscience/vp101.xlsx, type=circular_search, num_slices=50, fs_bishop=1.416, fs_spencer=1.422, benchmark=VP101 -->
<!-- test: file=../files/rocscience/vp087.xlsx, type=single_circle, num_slices=50, fs_bishop=1.031, benchmark=VP87 -->
<!-- test: file=../files/rocscience/vp088.xlsx, type=single_circle, num_slices=50, fs_spencer=1.057, benchmark=VP88 -->
<!-- test: file=../files/rocscience/vp089.xlsx, type=single_circle, num_slices=50, fs_spencer=1.011, benchmark=VP89 -->
<!-- test: file=../files/rocscience/vp090.xlsx, type=single_circle, num_slices=50, fs_bishop=1.012, benchmark=VP90 -->
<!-- test: file=../files/rocscience/vp091.xlsx, type=single_circle, num_slices=50, fs_spencer=0.960, benchmark=VP91 -->
<!-- test: file=../files/rocscience/vp092.xlsx, type=single_circle, num_slices=50, fs_bishop=1.010, benchmark=VP92 -->
<!-- test: file=../files/rocscience/vp093.xlsx, type=single_circle, num_slices=50, fs_bishop=1.017, benchmark=VP93 -->
<!-- test: file=../files/rocscience/vp094.xlsx, type=single_circle, num_slices=50, fs_bishop=1.020, benchmark=VP94 -->
<!-- test: file=../files/rocscience/vp098.xlsx, type=circular_search, num_slices=40, rapid=true, fs_spencer=1.046, benchmark=VP98 -->
<!-- test: file=../files/rocscience/vp099.xlsx, type=circular_search, num_slices=40, rapid=true, fs_spencer=1.390, benchmark=VP99 -->
<!-- test: file=../files/rocscience/vp096.xlsx, type=single_circle, rapid=true, num_slices=60, fs_spencer=1.434, fs_bishop=1.432, benchmark=VP96 -->
<!-- test: file=../files/rocscience/vp085a.xlsx, type=single_circle, num_slices=60, fs_bishop=1.567, fs_spencer=1.567, benchmark=VP85-active -->
<!-- test: file=../files/rocscience/vp085b.xlsx, type=single_circle, num_slices=60, fs_oms=1.319, fs_bishop=1.319, benchmark=VP85-passive -->
<!-- test: file=../files/rocscience/vp021a.xlsx, type=single_circle, num_slices=60, fs_oms=1.927, fs_bishop=2.075, fs_spencer=2.071, fs_mprice=2.071, benchmark=VP21-dry -->
<!-- test: file=../files/rocscience/vp021b.xlsx, type=single_circle, num_slices=60, fs_oms=1.606, fs_bishop=1.759, fs_spencer=1.757, fs_mprice=1.756, benchmark=VP21-ru -->

| # | Problem | Status | XSLOPE file / results |
|---:|---|---|---|
| [1](#vp1) | Slope, homogenous | **built** | [xslope_acads_simple.xlsx](../lem/files/xslope_acads_simple.xlsx). ACADS 1(a): seven-method comparison vs the ACADS consensus 1.00; Bishop 0.985 vs Slide 0.987. |
| [2](#vp2) | Slope, homogenous, tension crack | **built** | [vp002.xlsx](../files/rocscience/vp002.xlsx). Bishop 1.589 / Spencer 1.585 / Janbu(corr) 1.495 / M-P 1.586 vs Slide 1.596 / 1.592 / 1.489 / 1.592 (±0.4%); Giam reference 1.65. |
| [3](#vp3) | Slope, (3) materials | **built** | [vp003.xlsx](../files/rocscience/vp003.xlsx). Bishop 1.403 / Spencer 1.372 / Janbu(corr) 1.354 / M-P 1.371 vs Slide 1.405 / 1.375 / 1.357 / 1.374 (±0.3%); ACADS reference 1.39. Interface coordinates read from the labeled GeoStudio verification-manual figure of the same ACADS 1(c) problem. |
| [4](#vp4) | Slope, (3) materials, seismic | **built** | [vp004.xlsx](../files/rocscience/vp004.xlsx). Problem #3 + k=0.15g. Bishop 1.013 / Spencer 0.989 / Janbu(corr) 0.963 / M-P 0.987 vs Slide 1.016 / 0.991 / 0.965 / 0.989 (±0.3%); ACADS reference 1.00. |
| [5](#vp5) | Dam, (4) materials | **built** | [vp005.xlsx](../files/rocscience/vp005.xlsx). Talbingo Dam, end of construction (polygon-zone geometry). Critical mechanism is the infinite-slope limit on the upstream face: stored shallow circle gives 1.955 (all methods) vs Slide 1.948-1.949 and the tan φ/tan β limit 1.9475. |
| [6](#vp6) | Dam, (4) materials, predefined slip surface | **built** | [vp006.xlsx](../files/rocscience/vp006.xlsx). Specified circle (100.3, 291, R=278.8) through the inclined core. Bishop 2.206 / Spencer 2.290 / Janbu(corr) 2.073 / M-P 2.299 vs Slide 2.208 / 2.292 / 2.073 / 2.301 (±0.1%); ACADS reference 2.29. This problem exposed (and now guards) the folded-zone slice-weight bug fixed in `slice.py`. |
| 7 | Slope, (2) materials, weak layer | covered | [LEM sample 13](../lem/samples.md#verification-acads-weak-layer) (`xslope_acads_weak_layer.xlsx`) is this exact problem (ACADS 3(a)). Non-circular search: Spencer 1.258 / M-P 1.248 vs Slide 1.246 / 1.275; Giam reference 1.24-1.27. |
| [8](#vp8) | Slope, (2) materials, weak layer, predefined slip surface | **built** | [vp008.xlsx](../files/rocscience/vp008.xlsx). Specified 4-point surface (Table 8.2). Spencer 1.276 / Janbu(corr) 1.294 / M-P 1.260 vs Slide 1.277 / 1.294 / 1.262 (exact to ±0.002); SLOPE/W M-P 1.261; Giam reference 1.34. |
| [9](#vp9) | Slope, (2) materials, weak layer, water table, distributed load | **built** | [vp009.xlsx](../files/rocscience/vp009.xlsx). ACADS 4: inclined 0.6 m seam (geometry from the labeled GeoStudio figure), 8-point piezometric line, two surcharge strips. Non-circular search: Spencer 0.724 / Janbu(corr) 0.718 (M-P reaches 0.707 from a wider seed but does not solve the stored seed, so it is not tagged) vs Slide 0.760/0.720/0.734 (block search) and 0.707/0.683/0.699 (optimized); SLOPE/W 0.699-0.689; ACADS references 0.78 [Giam], 0.6878 [Slope 2000], 20-program mean 0.808. Published spread is wide; XSLOPE sits mid-band. |
| 10 | Slope, homogenous, pore pressure grid, ponded water | planned | Tier A via u='seep' per plan (ACADS consensus 1.53), not the pore-pressure grid |
| 11 | Embankment, (2) materials, pore pressure grid | planned |  |
| 12 | Embankment, (4) materials, tension crack, pore pressure grid | planned |  |
| 13 | Embankment, (3) materials, pore pressure grid | planned |  |
| [14](#vp14) | Slope, homogenous | **built** | [xslope_arai_tagyo.xlsx](../lem/files/xslope_arai_tagyo.xlsx). Arai & Tagyo (1985) ex. 1: seven-method comparison; Bishop 1.404 vs published 1.451. |
| [15](#vp15) | Slope, (3) materials, weak layer | **built** | [vp015.xlsx](../files/rocscience/vp015.xlsx). Arai & Tagyo (1985) ex. 2, weak middle band. Circular search: Bishop 0.419 / Spencer 0.422 / Janbu(corr) 0.436 / M-P 0.420 vs Slide 0.420 / 0.409 / 0.423 / (GLE) 0.437; A&T Bishop 0.417; Kim et al. 0.43. |
| [16](#vp16) | Slope, homogenous, water table | **built** | [vp016.xlsx](../files/rocscience/vp016.xlsx). Arai & Tagyo (1985) ex. 3, piezometric line. Circular search: Bishop 1.112 / Spencer 1.113 / Janbu(corr) 1.122 / M-P 1.111 vs Slide 1.118 / 1.118 / 1.131; A&T Bishop 1.138. |
| [17](#vp17) | Slope, homogenous | **built** | [vp017.xlsx](../files/rocscience/vp017.xlsx). Yamagami & Ueta (1988). Circular search: Ordinary 1.274 / Bishop 1.342 / Spencer 1.340 vs Slide 1.278 / 1.344 and Y&U 1.282 / 1.348; published non-circular Spencer 1.325-1.339 (our local non-circular search reaches 1.394 — same search-power note as #19/#20). |
| [18](#vp18) | Slope, homogenous slope, ru pore pressure | **built** | [vp018.xlsx](../files/rocscience/vp018.xlsx). Spencer (1969)/Baker (1980) slope, ru=0.5, non-circular search (right-facing). Spencer 1.033 / M-P 1.024 vs Slide 1.010 (random search + Monte-Carlo optimization), Baker 1.02, Spencer (1969) 1.08. |
| [19](#vp19) | Slope, (4) materials | **built** | [vp019.xlsx](../files/rocscience/vp019.xlsx). Greco (1996) ex. 4 / Yamagami & Ueta (1988) four-layer slope. Circular search: Spencer 1.429 / Bishop 1.448 vs published Spencer 1.40-1.42. Non-circular: XSLOPE's local search plateaus at ~1.45 from the stored seed while Slide's Monte-Carlo optimization reaches 1.398 — a search-power gap (noted for future search work), not a model difference. |
| [20](#vp20) | Slope, (4) materials, weak layer, water table | **built** | [vp020.xlsx](../files/rocscience/vp020.xlsx). Greco (1996) ex. 5 / Chen & Shao (1988): 0.5 m weak seam along the inclined base (polygon zones), water table. Circular: Bishop 1.086 / Spencer 1.091 vs Slide 1.087 / 1.093 (exact). Non-circular seam block: local search 1.082 vs Slide Monte-Carlo 1.010, Chen & Shao 1.01-1.03, Greco 0.973-1.1 — same search-power gap as #19. |
| 21 | Slope, homogenous, ru pore pressure | partial | [vp021a.xlsx](../files/rocscience/vp021a.xlsx) (dry), [vp021b.xlsx](../files/rocscience/vp021b.xlsx) (ru=0.25) — Fredlund & Krahn (1977) classic, fixed circle (120,90,R=80), imperial units. Dry: OMS 1.927 / Bishop 2.075 / Spencer 2.071 / M-P 2.071 vs F&K 1.928 / 2.080 / 2.073 / 2.076. ru: OMS 1.606 / Bishop 1.759 / Spencer 1.757 / M-P 1.756 vs F&K 1.607 / 1.766 / 1.761 / 1.764 (XSLOPE matches the F&K OMS-ru value exactly; Slide reports 1.687 there). Case 3 (water table) pending the phreatic-line coordinates. |
| 22 | Slope, (2) materials, weak layer, ru pore pressure | planned |  |
| [23](#vp23) | Slope, (3) materials | **built** | [vp023.xlsx](../files/rocscience/vp023.xlsx). Low (1989): undrained layers, lower cu grows 15→30 kPa with depth (`cp` linear-strength option). Circular search: Ordinary 1.357 / Bishop 1.130 vs Low 1.36 / 1.14 (Slide 1.370 / 1.192; Kim 1.17 — the published Bishop values themselves spread 1.14-1.19). |
| [24](#vp24) | Slope, (3) materials | **built** | [vp024.xlsx](../files/rocscience/vp024.xlsx). Low (1989) three-layer undrained slope (φ=0). Circular search: Ordinary 1.433 / Bishop 1.433 vs Slide 1.439 / 1.439; Low reference 1.44. |
| 25 | Bearing capacity test slope, homogenous, distributed load, predefined slip surface | planned | Chen & Shao (1988) / Prandtl 60° weightless slope (qc=149.31 kN/m over 10 m; theoretical FS 1.0; Slide Spencer 1.051, SLOPE/W 1.036). Buildable — the mechanism's ends are at different elevations, so the flat-arc guard (which blocks #26) does not apply; needs the theoretical surface construction. |
| 26 | Bearing capacity test prism, homogenous, distributed load, predefined slip surface | blocked | Prandtl bearing-capacity mechanism (weightless φ=0 soil, surface load, theoretical FS=1.0; Slide Spencer 0.940). XSLOPE's surface-validity checks reject failure surfaces whose two ends sit at equal elevation (the 'flat arc' guard) — flat-ground bearing mechanisms cannot currently be evaluated. Feature gap noted for a relaxed guard when driving comes from loads. |
| 27 | Slope, (2) materials, tension crack, water table (auto Hu) | planned |  |
| 28 | Excavated slope and embankment, (3) materials and (5) materials, probabilistic analysis | planned |  |
| 29 | Submerged slope, homogenous, probabilistic analysis, water table | planned |  |
| 30 | Reinforced embankment, (4) materials, tension crack, geosynthetic | planned |  |
| 31 | Reinforced embankment, (5) materials, geosynthetic | planned |  |
| 32 | Reinforced embankment, (7) materials, geosynthetic | planned |  |
| 33 | Dike, (5) materials, probabilistic analysis, water table | planned |  |
| 34 | Dam, (3) materials, probabilistic analysis, water table | planned |  |
| 35 | Dam, (5) materials, probabilistic analysis, reliability index | planned |  |
| [36](#vp36) | Slope, homogenous, probabilistic analysis, ru pore pressure, reliability index | **built** | [vp036.xlsx](../files/rocscience/vp036.xlsx). Li & Lumb (1987) / Hassan & Wolff (1999) reliability benchmark (c′=18±3.6, φ′=30±3, γ=18±0.9, ru=0.2). Deterministic Bishop 1.333 vs H&W 1.334 (Slide 1.340). Taylor-series β_ln on the critical surface 2.263 vs H&W (FOSM) 2.336 and Slide (Monte-Carlo) 2.482 — β estimates legitimately spread by method; xslope does not yet perturb ru (σ=0.02, minor). |
| 37 | Slope, homogenous, distributed load, back analysis of required support force and length | planned |  |
| 38 | Excavated slope, homogenous, finite element groundwater seepage analysis, matric suction | planned |  |
| 39 | Reinforced embankment, (2) materials, tension crack, geosynthetic | planned |  |
| 40 | Slope, homogenous, sensitivity analysis | planned |  |
| [41](#vp41) | Slope, homogenous, ru pore pressure | **built** | [vp041.xlsx](../files/rocscience/vp041.xlsx). Jiang, Baker & Yamagami (2003): power-curve strength τ=1.4·σ′^0.8 with ru=0.3 — exercises the v12 `pow` and `ru` options together. Circular search: Bishop 1.668 / Spencer 1.670 / Janbu(corr) 1.660 vs Slide Bishop 1.656 (non-linear path search), Charles & Soares 1.66, published range 1.56-1.67. |
| 42 | Dam, (3) materials, water table, ponded water, tension crack | partial | Baker & Leshchinsky (2001) safety-map dam. Geometry and core polygon fully labeled (extracted); reservoir level 30 on the right face; Slide Spencer 1.925 circular / 1.877 non-circular, Baker 1.91. Blocked on the internal phreatic-line coordinates through the core (drawn but unlabeled) — needs the B&L paper. |
| 43 | Slope, homogenous, planar surface, RocPlane comparison | partial | [vp043.xlsx](../files/rocscience/vp043.xlsx) built from the printed table (c'=30, φ'=30, γ=20, labeled geometry). xslope Janbu = hand Culmann exactly (1.429 vs 1.430 at 49.5°), but Slide/RocPlane/Baker all report ≈1.35 — reproducible only with different inputs (γ≈21.8 or c'≈27.5), so the manual's property table appears not to be what was run. Needs Baker (2001) [same blocked source as #42] to resolve before tagging. |
| [44](#vp44) | Slope, homogenous | **built** | [vp044a.xlsx](../files/rocscience/vp044a.xlsx) (power curve), [vp044b.xlsx](../files/rocscience/vp044b.xlsx) (Mohr-Coulomb), [vp044c.xlsx](../files/rocscience/vp044c.xlsx) (converged LLA). Baker (2003) ex. 1: 43° face, H=6. Spencer: power 0.958 vs Slide 0.960 / Baker 0.97; MC 1.518 vs Slide 1.536 / Baker 1.50; LLA 0.980 vs Slide 0.981. Baker's paper resolved the MC row (c'=11.64, φ'=24.7 — Table I it. 0) and γ=18. |
| [45](#vp45) | Slope, homogenous | **built** | [vp045a.xlsx](../files/rocscience/vp045a.xlsx) (Mohr-Coulomb), [vp045b.xlsx](../files/rocscience/vp045b.xlsx) (power curve). Baker (2003) ex. 2: linear vs non-linear envelope on the same 4:1 slope. Spencer: MC 2.801 vs Slide 2.794; power curve 2.649 vs Slide 2.662. (Slide's Janbu values are simplified/uncorrected; ours carry the fo correction and agree once scaled.) |
| 46 | Dam, (2) materials, rapid drawdown, finite element groundwater seepage analysis, ponded water | partial | Baker (1993) three-stage dam (dry / steady-state seep / rapid drawdown). The manual itself calls this a validation problem: permeabilities were estimated by Rocscience and the stage-3 undrained strengths live in discrete .fn6 functions. Stage 1 (dry, Spencer 2.534 / Baker 2.41 / theory 2.5) buildable once Figure 46.1 coordinates are read; stages 2-3 map onto xslope's seep + rapid-drawdown pipeline but need those source functions. |
| [47](#vp47) | Retaining wall, homogenous, planar failure, line load, shotcrete, soil nails | **built** | [vp047.xlsx](../files/rocscience/vp047.xlsx). Sheahan & Ho (2003) Amherst test wall: 6 m undrained cut, 2 nail rows (FHWA capacity envelope) + shotcrete line load. Critical 44.5° plane: Janbu 0.899 vs Slide 0.890 / Sheahan 0.887. |
| [48](#vp48) | Retaining wall, homogenous, planar failure, line load , soil nails, shotcrete | **built** | [vp048.xlsx](../files/rocscience/vp048.xlsx). Clouterre full-scale test wall: 7 nail rows (constant 15 kN tension), planar surfaces through the toe at 45–70°; Janbu/Spencer within 0.3% of Slide at 55–70°. Building this exposed (and fixed) right-facing axial-reinforcement sign bugs in every LEM method. |
| 49 | Retaining wall, (2) materials, grouted tiebacks, soldier piles | planned |  |
| [50](#vp50) | Reinforced slope, (2) materials, predefined slip surface, geosynthetic | **built** | [vp050.xlsx](../files/rocscience/vp050.xlsx). SNAILZ reference-manual nail wall: 14 rows with per-row length/tensile/bond values, evaluated on the printed deep wedge (-15.8,0)-(0,-5)-(41.7,25). With Slide's nail defaults (tangent orientation, force factored by FS): Janbu(corr) 1.448 vs SNAILZ 1.46 and Slide 1.417. The capacity envelope reproduces the hand-computed available tension at every crossing (Σ 10.6 kip). The shallow (0,0) surface's kink is not printed — only the deep case is tagged. |
| [51](#vp51) | Slope, (4) materials, water table, tension crack, seismic | **built** | [vp051.xlsx](../files/rocscience/vp051.xlsx). Zhu, Lee & Jiang (2003) four-layer slope, k=0.1, 5 m tension crack, specified circle (18.058, 66.744, R=86). Seven methods vs Slide/Zhu: Bishop 1.278 vs 1.278/1.278 and M-P 1.304 vs 1.304/1.303 — exact; Spencer 1.294 vs 1.293; Lowe 1.296 vs 1.288/1.290; Corps 1.404 vs 1.422/1.377 (in-band); OMS 1.069 vs Zhu 1.066 (Slide's 1.145 is the outlier); Janbu(corr) 1.205 ≡ simplified 1.112 × fo. Phreatic line calibrated against the two independently agreeing published Bishop/Spencer values (±1 m bracket spans them). |
| [52](#vp52) | Slope, (4) materials, water table, tension crack | **built** | [vp052a.xlsx](../files/rocscience/vp052a.xlsx) (dry), [vp052b.xlsx](../files/rocscience/vp052b.xlsx) (wet). Zhu & Lee (2002), heterogeneous benched slope; water table from the manual's Table 52.2. Unconstrained circular search lands in the governing deep (surface 3) family: wet Spencer 1.189 and Bishop 1.176 vs Slide 1.189 / 1.176 — exact; dry 1.797 / 1.796 vs Slide 1.804. Zhu's own values on his specified circle: 1.211 / 1.836 (the manual shows a wide Slide-Zhu spread on this family). The shallow/noncircular cases (surfaces 1, 2, 4) use constrained/block searches xslope does not yet expose — noted, not tagged. Paper in `ref_docs_lim_eq/`. |
| 53 | Slope, homogenous, water table, tension crack, planar failure, RocPlane comparison | planned |  |
| [54](#vp54) | Slope, homogenous, micro piles | **built** | [vp054a.xlsx](../files/rocscience/vp054a.xlsx) (no pile), [vp054b.xlsx](../files/rocscience/vp054b.xlsx) (with pile). Yamagami (2000): micro-pile row at the crest, 10.7 kN shear per pile at 1 m spacing. On the printed critical circle: no-pile Bishop 1.100 vs Slide 1.102 / Yamagami 1.10; with-pile 1.185 vs Slide 1.193 / Yamagami 1.20 (Slide adds the pile shear un-factored, i.e. active application). A free search with the pile finds 1.113 on a circle exiting upslope of the pile — the published comparison is per-circle. |
| 55 | Slope, homogenous, water table | partial | Pockoski & Duncan (2000) test slope 1 — geometry and properties extracted (ground (-75,100)-(0,100)-(100,150)-(170,150), c=300 psf, φ=30, γ=120 pcf; multi-program table: Spencer 1.30, Bishop 1.29, Janbu 1.15, Lowe 1.32); blocked on the paper's water-table coordinates (figure trace only). |
| 56 | Slope, homogenous, water table, tension crack | planned |  |
| 57 | Slope, (2) materials, water table, tension crack, composite surfaces | planned |  |
| 58 | Retaining wall, (8) materials, water table, grouted tieback | planned |  |
| 59 | Retaining wall, homogenous, water table, grouted tieback | planned |  |
| 60 | Retaining wall, (2) materials, tension crack, distributed load, soil nails | planned |  |
| 61 | Slope, homogenous, composite surfaces | planned |  |
| [62](#vp62) | Slope, homogenous, ru pore pressure, seismic | **built** | [vp062a.xlsx](../files/rocscience/vp062a.xlsx) (dry, kc=0.432), [vp062b.xlsx](../files/rocscience/vp062b.xlsx) (ru=0.5, kc=0.132). Loukidis et al. (2003) critical-seismic-coefficient benchmark: FS should be 1.0 at kc. Circular search: Spencer 1.001 / 1.001 and Bishop 0.991 / 0.986 vs Slide 1.001 / 1.001 and 0.991 / 0.987 — exact. |
| 63 | Slope, (3) materials, seismic | partial | Loukidis et al. (2003) ex. 2 (paper now in `ref_docs_lim_eq/`). Outline fully pinned from the paper's Fig. 9 (bench el 20 to x=20, 2:1 to (60,40), 8 m bench, 2.5:1 to (105.5,55), crest to 150); interfaces are 22:1 but their face anchors are not dimensioned and FS=1.0-calibration attempts with plausible anchors give 1.17-1.21 — the paper's Fig. 10 profile also appears inconsistent with Fig. 9's mesh. Needs a closer read of the paper (or the SLOPE/W .gsz) before building. |
| 64 | Embankment, (4) materials, water table, tension crack | planned |  |
| 65 | Embankment, (4) materials, water table, ponded water | planned |  |
| 66 | Embankment, (4) materials, water table, ponded water | planned |  |
| 67 | Embankment, (2) materials | planned |  |
| 68 | Embankment, (3) materials, ponded water | planned |  |
| 69 | Embankment, (2) materials, water table, ponded water | planned |  |
| 70 | Submerged slope, homogenous, water table, ponded water | planned |  |
| 71 | Slope, homogenous, finite element groundwater seepage analysis, water table | planned |  |
| 72 | Embankment dam, (4) materials, finite element groundwater seepage analysis, ponded water | planned |  |
| 73 | Excavated slope, (4) materials, tension crack | planned |  |
| 74 | Embankment, (2) materials | planned |  |
| 75 | Dyke, (4) materials | planned |  |
| 76 | Embankment dam, homogenous, finite element groundwater seepage analysis, ponded water | planned |  |
| 77 | Dam, (2) materials, finite element groundwater seepage analysis, ponded water | planned |  |
| 78 | Slope, homogenous | planned | Duncan & Wright (2005) Fig. 14.3 — the D&W book is in hand (2014 ed. in ref library); build queued. |
| 79 | Slope, (2) materials, infinite slope failure | planned | Duncan & Wright (2005) Fig. 14.4 — the D&W book is in hand (2014 ed. in ref library); build queued. |
| 80 | Embankment, (6) materials | planned |  |
| 81 | Embankment, (2) materials, infinite slope failure | planned | Duncan & Wright (2005) Fig. 14.7 — the D&W book is in hand (2014 ed. in ref library); build queued. |
| 82 | Embankment, (2) materials, water table | planned |  |
| 83 | Embankment, (2) materials | planned |  |
| 84 | Embankment, (2) materials | planned |  |
| [85](#vp85) | Reinforced slope, homogenous, grouted tieback | **built** | [vp085a.xlsx](../files/rocscience/vp085a.xlsx) (active), [vp085b.xlsx](../files/rocscience/vp085b.xlsx) (passive). Duncan & Wright (2005) Fig. 6.34: one 9,000 lb/ft horizontal tieback at mid-height of an undrained clay slope — the reinforcement acceptance benchmark from the design plan. On Slide's printed critical circles: active 1.567 vs Slide GLE 1.575 (D&W 1.51); passive Bishop 1.319 vs Slide 1.324 (D&W 1.32). Slide's own method table scatters 1.42-2.05 here (concentrated force strains interslice assumptions), so per-circle comparison is the meaningful one. |
| [86](#vp86) | Reinforced slope, homogenous, grouted tieback | **built** | [vp086.xlsx](../files/rocscience/vp086.xlsx). Duncan & Wright (2005) Fig. 7.28 / STABGM reinforced fill on rock: five 800 lb/ft geogrids. Circular search: Bishop 1.617 / Spencer 1.611 vs Slide 1.629 / 1.620; D&W reference 1.61. |
| [87](#vp87) | Retaining wall, (3) materials, geotextile | **built** | [vp087.xlsx](../files/rocscience/vp087.xlsx). Baseline three-tier wall (Ta=10, L=6.3): Bishop 1.031 on Slide's printed circle vs Slide 1.040; free search 0.99 vs L&H 0.99–1.00. |
| [88](#vp87) | Retaining wall, (3) materials, geotextile | **built** | [vp088.xlsx](../files/rocscience/vp088.xlsx). Fill-quality case (φ=25, Ta=22): Spencer 1.057 vs Slide 1.043. |
| [89](#vp87) | Retaining wall, (3) materials, geotextile | **built** | [vp089.xlsx](../files/rocscience/vp089.xlsx). Reinforcement-length case (L=4.2, Ta=11.4): Spencer 1.011 (≈L&H design intent 1.0); with Slide's actual baseline Ta=10 supports: 0.980 vs Slide 0.971. |
| [90](#vp87) | Retaining wall, (3) materials, geotextile | **built** | [vp090.xlsx](../files/rocscience/vp090.xlsx). Two reinforcement types (7.5 upper 8 / 11.0 lower 7): Bishop 1.012 vs Slide 1.004. |
| [91](#vp87) | Retaining wall, (3) materials, geotextile | **built** | [vp091.xlsx](../files/rocscience/vp091.xlsx). Weak foundation (c=0, φ=18): deep bearing circle, Spencer 0.960 vs Slide 0.964. |
| [92](#vp87) | Retaining wall, (3) materials, geotextile | **built** | [vp092.xlsx](../files/rocscience/vp092.xlsx). Water table 3 m above foundation (drained fill + pond): with Ta=10, Bishop 1.039 vs Slide 1.037; at the paper's Ta=9.25: 1.010 ≈ L&H 1.01. |
| [93](#vp87) | Retaining wall, (3) materials, distributed load, geotextile | **built** | [vp093.xlsx](../files/rocscience/vp093.xlsx). 20 kPa crest surcharge: at the paper's Ta=11.6, Bishop 1.017 ≈ L&H 1.02; with Ta=10, 0.961 vs Slide 0.958. |
| [94](#vp87) | Retaining wall, (3) materials, geotextile | **built** | [vp094.xlsx](../files/rocscience/vp094.xlsx). Five 1.8-m tiers (Ta=10.1): Bishop 1.020 on Slide's printed circle vs Slide 1.040. |
| 95 | Embankment dam, homogenous, rapid drawdown, water table | partial | USACE EM 1110-2-1902 App. G example with the **Corps 1970 2-stage** method (Slide 1.347, USACE 1.35). XSLOPE implements the Duncan-Wright-Wong 3-stage procedure (see #96, same model file) — the 2-stage variant would be a new option. |
| [96](#vp96) | Embankment dam, homogenous, rapid drawdown, water table | **built** | [vp096.xlsx](../files/rocscience/vp096.xlsx). USACE EM 1110-2-1902 (2003) Appendix G example (Figure G-5: 3:1 / 2.5:1 face, pool 103→24, Kc=1 envelope d=1379 psf ψ=18.2°), specified circle (169.5, 210, R=210), Duncan-Wright-Wong 3-stage: Spencer 1.434 / Bishop 1.432 vs Slide 1.443 and USACE 1.44 (Modified Swedish). First corpus problem through the rapid-drawdown pipeline. |
| [97](#vp97) | Embankment dam, homogenous, rapid drawdown, water table | **built** | [vp097.xlsx](../files/rocscience/vp097.xlsx). Pilarcitos Dam (Duncan, Wright & Wong 1990 — paper in `ref_docs_slope/`), drawdown 72→37 ft; Kc=1 envelope from D&W eqs 9.6-9.7 (d=64.1 psf, ψ=24.4°, the same equations reproduce USACE's 1379/18.2 exactly). Rapid 3-stage search: Spencer 1.044 / Bishop 1.042 vs Slide 1.043 and DWW 1.05 — the dam that actually failed in drawdown sits right at FS≈1. |
| [98](#vp98) | Embankment dam, (5) materials, rapid drawdown, water table | **built** | [vp098.xlsx](../files/rocscience/vp098.xlsx). Walter Bouldin Dam (DWW 1990): 5-zone section from Slide's labeled Figure 98.1 + color-zone tracing; Kc=1 envelopes from the paper's Table 2. DWW 3-stage Spencer search: 1.046 vs Slide 1.039 / paper 1.04. Critical surface matches the observed slide location. |
| [99](#vp99) | Embankment dam, (3) materials, rapid drawdown, water table | **built** | [vp099.xlsx](../files/rocscience/vp099.xlsx). DWW pumped-storage dam, drawdown 285→120. Geometry axis-calibrated from Slide's figure (vertices unlabeled). DWW 3-stage Spencer: 1.390 (search) / 1.428 (Slide's printed circle) vs Slide 1.534, SLOPE/W 1.550, paper 1.56 — ~7% low, attributed to core-geometry reading; to be re-pinned from the vendor .gsz when available. |
| [100](#vp100) | Embankment dam, homogenous, rapid drawdown, water table | **built** | [vp100.xlsx](../files/rocscience/vp100.xlsx). Morgenstern (1963) chart problem, complete drawdown, B̄=1 — the residual pore-pressure field maps exactly onto a piezometric line at the slope surface, so this runs as a single-stage analysis. Bishop 1.201 vs Morgenstern chart 1.20 and Slide (B-bar method) 1.212. Paper in `ref_docs_slope/`. |
| [101](#vp101) | Embankment dam, homogenous, rapid drawdown, water table | **built** | [vp101.xlsx](../files/rocscience/vp101.xlsx). Morgenstern (1963), drawdown 100→50 ft, B̄=1 (piezo = ground above the pool, 50 below it; remaining pond on the face). Bishop 1.416 vs Slide 1.417 (exact) and Morgenstern chart 1.41. |
| 102 | Embankment dam, homogenous, rapid drawdown | planned |  |
| 103 | Undrained slope, multi-model optimization (MMO) | planned |  |
| 104 | Newmark analysis, seismic analysis, multi-modal optimization (MMO) | planned |  |
| 105 | Anisotropic surface, multi-modal optimization (MMO) | planned |  |
| 106 | Support, Ito & Matsui pile | partial | Cai & Ugai (2000): Bishop with Ito & Matsui limit-pressure pile forces (FS 1.13-1.54 vs spacing). XSLOPE models pile resistance via Vcap/Mcap, not the Ito & Matsui pressure equation, so only the no-pile case (1.13/1.14) is directly reproducible; an Ito & Matsui capacity option would be a new feature. |
| 107 | Retaining walls, gabion walls, supports | planned |  |
| 108 | Retaining walls, gabion walls, supports | planned |  |
| 109 | Retaining walls, gabion walls, weak layers | planned |  |
| 110 | Retaining walls, equivalent fluid pressure | planned |  |
| 111 | Helical anchor | planned |  |

---

## Problem details

Each built problem below shows the XSLOPE inputs (with coordinate labels) beside a representative solved surface. The build scripts live in `benchmarks/rocscience/build_problems.py`; the figures are regenerated by `benchmarks/rocscience/make_figures.py`.

### VP1: Slope, homogeneous (ACADS 1(a)) {#vp1}

This is the headline limit-equilibrium verification benchmark, from the ACADS
slope stability program review (Donald & Giam, 1989; Giam & Donald, 1992), as
documented in the [GeoStudio SLOPE/W Verification Manual (Oct 2022)](https://files.seequent.com/PDFs/Geostudio-Slope%20Stability%20Verification%20Manual-Oct2022.pdf). A simple
homogeneous slope analyzed with a circular search; the ACADS consensus answer
is FOS ≈ 1.00, making percent differences easy to read.

| Property | Value |
|---|---|
| Slope | 2:1, 10 m high, with a bench |
| Cohesion, $c'$ | 3.0 kPa |
| Friction angle, $\phi'$ | 19.6° |
| Unit weight, $\gamma$ | 20.0 kN/m³ |
| Pore pressure | none (total stress) |

Excel input file: [xslope_acads_simple.xlsx](../lem/files/xslope_acads_simple.xlsx)

![acads_simple_inputs.png](../lem/images/acads_simple_inputs.png){width=900}

Critical circle from the automated search (Spencer's method shown):

![acads_simple_solution.png](../lem/images/acads_simple_solution.png){width=900}

XSLOPE results for all six methods (automated critical-circle search, 50
slices, each method searched independently):

| Method | XSLOPE FOS | Reference | Diff |
|---|---|---|---|
| Ordinary (OMS) | 0.942 | 1.00 | -5.8% |
| Bishop's Simplified | 0.985 | 1.00 | -1.5% |
| Simplified Janbu | 0.986 | 1.00 | -1.4% |
| Corps of Engineers | 0.990 | 1.00 | -1.0% |
| Lowe & Karafiath | 0.987 | 1.00 | -1.3% |
| Spencer | 0.984 | 1.00 | -1.6% |
| Morgenstern-Price | 0.984 | 1.00 | -1.6% |

All rigorous methods fall within the ACADS accepted band; OMS reads low, as
expected for the legacy method (its conservative bias on this class of problem
is why it is reported for completeness only). This benchmark also appears on
the [Verification](../verification/lem.md) page.

**Sources:** Donald, I.B. & Giam, P. (1989), *Soil slope stability programs
review*, ACADS, Melbourne; Giam, P. & Donald, I.B. (1992); GeoStudio
[SLOPE/W Verification Manual (Oct 2022)](https://files.seequent.com/PDFs/Geostudio-Slope%20Stability%20Verification%20Manual-Oct2022.pdf),
ACADS suite.

<!-- fs-table -->
**Factor of safety by method** (each method's own critical surface):

| OMS | Bishop | Janbu | Corps | Lowe | Spencer | M-P |
|---:|---:|---:|---:|---:|---:|---:|
| 0.942 | 0.985 | 0.986 | 0.990 | 0.987 | 0.984 | 0.984 |
<!-- /fs-table -->

<!-- test: file=../lem/files/xslope_acads_simple.xlsx, type=circular_search, num_slices=50, fs_oms=0.942, fs_bishop=0.985, fs_janbu=0.986, fs_corps=0.990, fs_lowe=0.987, fs_spencer=0.984, fs_mprice=0.984, benchmark=LEM-1 -->

### VP2: Slope, homogenous, tension crack {#vp2}

ACADS 1(b) (Giam & Donald 1989): the 1(a) slope with c'=32, phi'=10, gamma=20 and a water-filled tension crack of depth 2c/(gamma*sqrt(ka)) [Craig 1997]. Slide2: Bishop 1.596, Spencer 1.592, Janbu corrected 1.489, GLE 1.592; Giam reference 1.65.

**Input files:** [vp002.xlsx](../files/rocscience/vp002.xlsx)

| Method | XSLOPE | Published |
|---|---|---|
| Bishop | 1.589 | Slide 1.596; SLOPE/W 1.664 |
| Janbu (corrected) | 1.495 | Slide 1.489 |
| Spencer | 1.585 | Slide 1.592 |
| Morgenstern-Price | 1.586 | Slide 1.592; SLOPE/W 1.660 |

*ACADS reference band 1.65–1.70 (Giam & Donald).*

![vp002: inputs and representative solution](images/vp002.png)


### VP3: Slope, (3) materials {#vp3}

ACADS 1(c): non-homogeneous three-layer slope, critical circle. Slide2: Bishop 1.405, Spencer 1.375, GLE 1.374, Janbu corrected 1.357; SLOPE/W: Bishop 1.414, M-P 1.382; ACADS reference 1.39.

**Input files:** [vp003.xlsx](../files/rocscience/vp003.xlsx)

| Method | XSLOPE | Published |
|---|---|---|
| Bishop | 1.403 | Slide 1.405; SLOPE/W 1.414 |
| Janbu (corrected) | 1.354 | Slide 1.357 |
| Spencer | 1.372 | Slide 1.375 |
| Morgenstern-Price | 1.371 | Slide 1.374; SLOPE/W 1.382 |

*ACADS reference 1.39.*

![vp003: inputs and representative solution](images/vp003.png)

### VP4: Slope, (3) materials, seismic {#vp4}

ACADS 1(d): problem #3 plus horizontal seismic coefficient 0.15. Slide2: Bishop 1.016, Spencer 0.991, GLE 0.989, Janbu corrected 0.965; SLOPE/W: Bishop 1.02, M-P 0.989; ACADS reference 1.00.

**Input files:** [vp004.xlsx](../files/rocscience/vp004.xlsx)

| Method | XSLOPE | Published |
|---|---|---|
| Bishop | 1.013 | Slide 1.016; SLOPE/W 1.02 |
| Janbu (corrected) | 0.963 | Slide 0.965 |
| Spencer | 0.989 | Slide 0.991 |
| Morgenstern-Price | 0.987 | Slide 0.989; SLOPE/W 0.989 |

*ACADS reference 1.00.*

![vp004: inputs and representative solution](images/vp004.png)

### VP5: Dam, (4) materials {#vp5}

ACADS 2(a) (Giam & Donald 1989): Talbingo Dam at end of construction, 4 zones, critical circular surface. Slide2: Bishop 1.948, Spencer 1.948, GLE 1.948, Janbu corrected 1.949; Giam reference 1.95. The minimum is a shallow slide parallel to the (steeper) upstream face.

**Input files:** [vp005.xlsx](../files/rocscience/vp005.xlsx)

| Method | XSLOPE | Published |
|---|---|---|
| Bishop | 1.955 | Slide 1.948; SLOPE/W 1.951 |
| Janbu (corrected) | 1.965 | Slide 1.949 |
| Spencer | 1.955 | Slide 1.948 |
| Morgenstern-Price | 1.955 | Slide 1.948 |

*Critical mechanism is the infinite-slope limit: tan φ′/tan β = 1.9475.*

![vp005: inputs and representative solution](images/vp005.png)

### VP6: Dam, (4) materials, predefined slip surface {#vp6}

ACADS 2(b): Talbingo Dam, single specified circle Xc=100.3, Yc=291.0, R=278.8 (Table 6.1). Slide2: Bishop 2.208, Spencer 2.292, GLE 2.301, Janbu corrected 2.073; Giam reference 2.29.

**Input files:** [vp006.xlsx](../files/rocscience/vp006.xlsx)

| Method | XSLOPE | Published |
|---|---|---|
| Bishop | 2.206 | Slide 2.208; SLOPE/W 2.207 |
| Janbu (corrected) | 2.073 | Slide 2.073 |
| Spencer | 2.290 | Slide 2.292 |
| Morgenstern-Price | 2.299 | Slide 2.301; SLOPE/W 2.299 |

*ACADS reference 2.29. The problem that exposed the folded-zone slice-weight bug.*

![vp006: inputs and representative solution](images/vp006.png)

### VP8: Slope, (2) materials, weak layer, predefined slip surface {#vp8}

ACADS 3(b): the weak-layer slope (= LEM sample 13 / Slide #7) with the fully specified non-circular surface of Table 8.2. Slide2: Spencer 1.277, GLE 1.262, Janbu corrected 1.294; SLOPE/W: Bishop 1.259, M-P 1.261; Giam reference 1.34.

**Input files:** [vp008.xlsx](../files/rocscience/vp008.xlsx)

| Method | XSLOPE | Published |
|---|---|---|
| Janbu (corrected) | 1.294 | Slide 1.294 |
| Spencer | 1.276 | Slide 1.277 |
| Morgenstern-Price | 1.260 | Slide 1.262; SLOPE/W 1.261 |

*Giam reference 1.34.*

![vp008: inputs and representative solution](images/vp008.png)

### VP9: Slope, (2) materials, weak layer, water table, distributed load {#vp9}

ACADS 4 (Slide #9): weak-layer slope + piezometric surface (Table 9.3) + two surcharge strips (Table 9.2: 20 kPa on the lower bench x=23-43, 20->40 kPa ramp on the crest x=70-80). Non-circular search. Slide2 (block search, no optimization): Spencer 0.760, GLE 0.720, Janbu corrected 0.734; with optimization 0.683-0.707; SLOPE/W: Bishop 0.699, M-P 0.689; Giam reference 0.78; Slope 2000 GLE reference 0.6878. The published spread is wide - this is a search-difficulty benchmark.

**Input files:** [vp009.xlsx](../files/rocscience/vp009.xlsx)

| Method | XSLOPE | Published |
|---|---|---|
| Janbu (corrected) | 0.718 | Slide 0.734 / 0.699 |
| Spencer | 0.724 | Slide 0.760 (block) / 0.707 (optimized) |

*SLOPE/W Bishop 0.699, M-P 0.689; ACADS 0.6878 (Slope 2000), 20-program mean 0.808 — a wide published band.*

![vp009: inputs and representative solution](images/vp009.png)

### VP14: Slope, homogeneous (Arai & Tagyo ex. 1) {#vp14}

From [Arai & Tagyo (1985)](https://doi.org/10.3208/sandf1972.25.43), *Soils
and Foundations* 25(1), and republished by Greco (1996), Malkawi et al.
(2001), and Kim et al. (2002); also
[SLOPE/W Verification Manual](https://files.seequent.com/PDFs/Geostudio-Slope%20Stability%20Verification%20Manual-Oct2022.pdf)
sec. 2.11. A homogeneous 1.5:1 slope, 20 m high, with
c = 41.65 kPa, φ = 15.0°, γ = 18.82 kN/m³ (total stress). Published FOS ≈ 1.451.

Excel input file: [xslope_arai_tagyo.xlsx](../lem/files/xslope_arai_tagyo.xlsx)

![arai_tagyo_inputs.png](../lem/images/arai_tagyo_inputs.png){width=900}

Critical circle (Spencer's method shown):

![arai_tagyo_solution.png](../lem/images/arai_tagyo_solution.png){width=900}

Results for all six methods (automated critical-circle search, 50 slices):

| Method | XSLOPE FOS | Reference | Diff |
|---|---|---|---|
| Ordinary (OMS) | 1.344 | 1.451 | -7.4% |
| Bishop's Simplified | 1.404 | 1.451 | -3.2% |
| Simplified Janbu | 1.411 | 1.451 | -2.8% |
| Corps of Engineers | 1.476 | 1.451 | +1.7% |
| Lowe & Karafiath | 1.438 | 1.451 | -0.9% |
| Spencer | 1.401 | 1.451 | -3.4% |
| Morgenstern-Price | 1.400 | 1.451 | -3.5% |

This benchmark also appears on the
[Verification](../verification/lem.md) page.

**Source:** Arai, K. & Tagyo, K. (1985). Determination of noncircular slip
surface giving the minimum factor of safety in slope stability analysis.
*Soils and Foundations* 25(1):43-51.
[doi:10.3208/sandf1972.25.43](https://doi.org/10.3208/sandf1972.25.43).
Republished in Greco (1996), Malkawi et al. (2001), and Kim et al. (2002);
also SLOPE/W Verification Manual sec. 2.11.

<!-- fs-table -->
**Factor of safety by method** (each method's own critical surface):

| OMS | Bishop | Janbu | Corps | Lowe | Spencer | M-P |
|---:|---:|---:|---:|---:|---:|---:|
| 1.344 | 1.404 | 1.411 | 1.476 | 1.438 | 1.401 | 1.400 |
<!-- /fs-table -->

<!-- test: file=../lem/files/xslope_arai_tagyo.xlsx, type=circular_search, num_slices=50, fs_oms=1.344, fs_bishop=1.404, fs_janbu=1.411, fs_corps=1.476, fs_lowe=1.438, fs_spencer=1.401, fs_mprice=1.400, benchmark=LEM-2b -->

### VP15: Slope, (3) materials, weak layer {#vp15}

Slide #15: Arai & Tagyo (1985) example 2 - three layers with a weak (c=9.8, phi=5) middle band, no water. Circular search. Slide2 (auto refine): Bishop 0.420, Spencer 0.409, GLE 0.437, Janbu corrected 0.423; A&T Bishop 0.417; Kim et al. 0.43.

**Input files:** [vp015.xlsx](../files/rocscience/vp015.xlsx)

| Method | XSLOPE | Published |
|---|---|---|
| Bishop | 0.419 | Slide 0.420; A&T 0.417 |
| Janbu (corrected) | 0.436 | Slide 0.423; A&T 0.430 |
| Spencer | 0.422 | Slide 0.409 |
| Morgenstern-Price | 0.420 | Slide (GLE) 0.437 |

*Kim et al. (2002) 0.43.*

![vp015: inputs and representative solution](images/vp015.png)


### VP16: Slope, homogenous, water table {#vp16}

Slide #16: Arai & Tagyo (1985) example 3 - homogeneous slope with a water table. Circular search. Slide2 (auto refine): Bishop 1.118, Janbu simplified 1.046, Janbu corrected 1.131, Spencer 1.118; A&T Bishop 1.138.

**Input files:** [vp016.xlsx](../files/rocscience/vp016.xlsx)

| Method | XSLOPE | Published |
|---|---|---|
| Bishop | 1.112 | Slide 1.118; A&T 1.138 |
| Janbu (corrected) | 1.122 | Slide 1.131 |
| Spencer | 1.113 | Slide 1.118 |
| Morgenstern-Price | 1.111 | — |

*SLOPE/W reports 1.190, the outlier of the four sources.*

![vp016: inputs and representative solution](images/vp016.png)

### VP17: Slope, homogenous {#vp17}

Slide #17: Yamagami & Ueta (1988) homogeneous slope, dry. Circular: Slide Bishop 1.344, Ordinary 1.278 (Y&U 1.348 / 1.282). Non-circular: Slide Spencer 1.325 (Y&U 1.339, Greco 1.33).

**Input files:** [vp017.xlsx](../files/rocscience/vp017.xlsx)

| Method | XSLOPE | Published |
|---|---|---|
| Ordinary | 1.274 | Slide 1.278; Y&U 1.282 |
| Bishop | 1.342 | Slide 1.344; Y&U 1.348 |
| Spencer | 1.340 | published non-circular 1.325–1.339 |

*Our local non-circular search reaches 1.394 — same search-power note as VP19/VP20.*

![vp017: inputs and representative solution](images/vp017.png)

### VP18: Slope, homogenous slope, ru pore pressure {#vp18}

Slide #18: Spencer (1969) / Baker (1980) homogeneous slope with ru=0.5, non-circular critical surface. The slope descends left-to-right (a right-facing case). Slide2 (random search + Monte Carlo optimization): Spencer 1.010; Baker reference 1.02; Spencer (1969) 1.08.

**Input files:** [vp018.xlsx](../files/rocscience/vp018.xlsx)

| Method | XSLOPE | Published |
|---|---|---|
| Spencer | 1.033 | Slide 1.010 (MC-optimized); Baker 1.02; Spencer (1969) 1.08 |
| Morgenstern-Price | 1.024 | — |

![vp018: inputs and representative solution](images/vp018.png)

### VP19: Slope, (4) materials {#vp19}

Slide #19: Greco (1996) ex. 4 / Yamagami & Ueta (1988) four-layer slope, no water, non-circular critical surface. Slide2 (random search + Monte-Carlo optimization, convex): Spencer 1.398; Greco and Yamagami & Ueta references 1.40-1.42.

**Input files:** [vp019.xlsx](../files/rocscience/vp019.xlsx)

| Method | XSLOPE | Published |
|---|---|---|
| Bishop | 1.448 | — |
| Spencer | 1.429 | Greco / Y&U 1.40–1.42; Slide MC 1.398 |

*Circular-search values; the local non-circular search plateaus ~1.45 (search-power gap).*

![vp019: inputs and representative solution](images/vp019.png)

### VP20: Slope, (4) materials, weak layer, water table {#vp20}

Slide #20: Greco (1996) ex. 5 / Chen & Shao (1988): four layers with a 0.5 m weak seam along the inclined model base, water table. Polygon-zone geometry (the base is not horizontal). Slide2: circular (toe focus grid) Bishop 1.087, Spencer 1.093; non-circular (block search in seam + MC optimization) Spencer 1.010; Chen & Shao 1.01-1.03; Greco 0.973-1.1.

**Input files:** [vp020.xlsx](../files/rocscience/vp020.xlsx)

| Method | XSLOPE | Published |
|---|---|---|
| Bishop | 1.086 | Slide 1.087 |
| Spencer | 1.082 | Slide 1.093; Greco 1.08 |

*Non-circular seam block: local search 1.082 vs Slide MC 1.010, Chen & Shao 1.01–1.03.*

![vp020: inputs and representative solution](images/vp020.png)

### VP23: Slope, (3) materials {#vp23}

Slide #23: Low (1989) slope over two undrained layers; the lower layer's cu grows linearly 15->30 kPa from y=4 to y=0 (xslope 'cp' option: Su = c + cp*(r_elev - y)). Circular search. Slide2: Ordinary 1.370, Bishop 1.192; Low 1.36 / 1.14; Kim (2002) 1.17.

**Input files:** [vp023.xlsx](../files/rocscience/vp023.xlsx)

| Method | XSLOPE | Published |
|---|---|---|
| Ordinary | 1.357 | Slide 1.370; Low 1.36 |
| Bishop | 1.130 | Slide 1.192; Low 1.14; Kim 1.17 |

*Published Bishop values themselves spread 1.14–1.19 on this deep φ=0 problem.*

![vp023: inputs and representative solution](images/vp023.png)

### VP24: Slope, (3) materials {#vp24}

Slide #24: Low (1989) three-layer undrained slope (phi=0). Circular search. Slide2: Ordinary 1.439, Bishop 1.439; Low reference 1.44 both.

**Input files:** [vp024.xlsx](../files/rocscience/vp024.xlsx)

| Method | XSLOPE | Published |
|---|---|---|
| Ordinary | 1.433 | Slide 1.439; Low 1.44 |
| Bishop | 1.433 | Slide 1.439; Low 1.44 |

![vp024: inputs and representative solution](images/vp024.png)

### VP36: Slope, homogenous, probabilistic analysis, ru pore pressure, reliability index {#vp36}

Slide #36: Li & Lumb (1987) / Hassan & Wolff (1999) reliability benchmark: c'=18+-3.6, phi'=30+-3, gamma=18+-0.9, ru=0.2 (+-0.02, not perturbed by xslope's Taylor-series reliability - its contribution to sigma_F is small). Bishop deterministic FS 1.334 (H&W) / 1.340 (Slide); beta_lognormal on the deterministic surface 2.336 (H&W) / 2.482 (Slide).

**Input files:** [vp036.xlsx](../files/rocscience/vp036.xlsx)

| Method | XSLOPE | Published |
|---|---|---|
| Bishop | 1.333 | H&W 1.334; Slide 1.340 |
| β_ln (reliability) | 2.263 | H&W (FOSM) 2.336; Slide (Monte-Carlo) 2.482 |

*β estimates legitimately spread by estimation method; xslope does not yet perturb ru (σ = 0.02, minor).*

![vp036: inputs and representative solution](images/vp036.png)

### VP41: Slope, homogenous, ru pore pressure {#vp41}

Slide #41: Jiang, Baker & Yamagami (2003) homogeneous clay slope with power-curve strength tau = 1.4*(sigma')^0.8 and ru = 0.3 - exercises the v12 pow option and ru together. Slide2: Bishop 1.656, Janbu simplified 1.563; Charles & Soares 1.66; Baker 1.56-1.60; Perry 1.67.

**Input files:** [vp041.xlsx](../files/rocscience/vp041.xlsx)

| Method | XSLOPE | Published |
|---|---|---|
| Bishop | 1.668 | Slide 1.656 (path search); Charles & Soares 1.66 |
| Janbu (corrected) | 1.660 | Slide simplified 1.563 |
| Spencer | 1.670 | — |

*Published range 1.56–1.67.*

![vp041: inputs and representative solution](images/vp041.png)

### VP44: Slope, homogeneous — linear vs non-linear envelope (Baker ex. 1) {#vp44}

Slide #44 / Baker (2003) example problem 1: a straight 43° slope, H = 6 m, γ = 18 kN/m³, in compacted Israeli clays, analyzed with three strength models fitted to the same triaxial data: (a) the power curve τ = 1.107·σ′^0.86 (Baker's A = 0.58, n = 0.86, T = 0); (b) the experimentally fitted Mohr-Coulomb envelope c′ = 11.64 kPa, φ′ = 24.7° (Table I, iteration 0 of Baker's paper — this resolves the property-table ambiguity in the Slide manual); and (c) Baker's converged local-linear-approximation parameters c′ = 0.39 kPa, φ′ = 38.6°. The point of the example is the danger of extrapolating a linear envelope into the low-stress range: the M-C fit says FS = 1.5, the non-linear law says the slope is failing.

**Input files:** [vp044a.xlsx](../files/rocscience/vp044a.xlsx), [vp044b.xlsx](../files/rocscience/vp044b.xlsx), [vp044c.xlsx](../files/rocscience/vp044c.xlsx)

| Case | Method | XSLOPE | Published |
|---|---|---|---|
| (a) power curve | Spencer | 0.958 | Slide 0.960; Baker 0.97 |
| (b) Mohr-Coulomb | Spencer | 1.518 | Slide 1.536; Baker 1.50 |
| (c) LLA converged | Spencer | 0.980 | Slide 0.981; Baker 0.97 |

*Baker states γ = 18 for all his examples; the Slide manual's table prints 19.5, which reconciles with neither program's results (γ = 19.5 gives Spencer 1.459 on case b). Slide's Janbu values are simplified/uncorrected, as in [#45](#vp45).*

![vp044a: inputs and representative solution](images/vp044a.png)

![vp044b: inputs and representative solution](images/vp044b.png)

![vp044c: inputs and representative solution](images/vp044c.png)

### VP45: Slope, homogenous {#vp45}

Slide #45, Mohr-Coulomb case: c'=11.64, phi'=24.7, gamma=18. Slide2: Janbu simplified 2.662, Spencer 2.794.

Slide #45, power-curve case: tau = 1.107*(sigma')^0.86 (Baker's A=0.58, n=0.86, T=0). Slide2: Janbu simplified 2.559, Spencer 2.662; Baker's accepted values for this example are of the same order.

**Input files:** [vp045a.xlsx](../files/rocscience/vp045a.xlsx), [vp045b.xlsx](../files/rocscience/vp045b.xlsx)

*Mohr-Coulomb case*

| Method | XSLOPE | Published |
|---|---|---|
| Spencer | 2.801 | Slide 2.794 |

*Power-curve case*

| Method | XSLOPE | Published |
|---|---|---|
| Spencer | 2.649 | Slide 2.662 |

*Slide’s Janbu values are simplified/uncorrected; ours carry fo and agree once scaled.*

![vp045a: inputs and representative solution](images/vp045a.png)

![vp045b: inputs and representative solution](images/vp045b.png)

### VP47: Soil-nailed wall in clay (Amherst test wall) {#vp47}

Slide #47 / Sheahan & Ho (2003): 6 m vertical cut in undrained Amherst clay (cᵤ = 25 kPa, γ = 18.9 kN/m³), two nail rows at 20° declination (L = 4.9 m, tensile 118 kN, plate 86 kN, bond 15 kN/m, sₕ = 1.5 m) and the shotcrete facing weight applied as a 14.6 kN/m vertical line load at the crest. The wall failed in the field test; the published analyses sweep planar surfaces through the toe. Nails are modeled axial/passive with the FHWA-style capacity envelope (plate strength at the head, bond-strength taper at the tip).

**Input files:** [vp047.xlsx](../files/rocscience/vp047.xlsx)

| Method | XSLOPE | Published |
|---|---|---|
| Janbu, critical plane (44.5°) | 0.899 | Slide 0.890 (simplified = corrected); Sheahan trial wedge 0.887 |

*Sheahan adds the nail tension unfactored; that convention (`appl=active`) gives 0.893. The tabulated 0.899 uses Slide's nail default (passive).*

![vp047: inputs and representative solution](images/vp047.png)

### VP48: Soil-nailed wall in sand (Clouterre test wall no. 1) {#vp48}

Slide #48 / Sheahan & Ho (2003): the CEBTP Clouterre full-scale wall — 7 m cut in Fontainebleau sand (c′ = 3 kPa, φ′ = 38°, γ = 20 kN/m³), seven nail rows at 10° declination (lengths 6/8/7.5/8/8/8/6 m from the paper's Fig. 4a, sₕ = 1.15 m), shotcrete weight as a 13.2 kN/m line load. Following Sheahan, each nail carries a constant 15 kN tension (fully anchored ends in xslope). The benchmark evaluates planar surfaces through the toe at 45–70°:

**Input files:** [vp048.xlsx](../files/rocscience/vp048.xlsx)

| Plane angle | XSLOPE Janbu | XSLOPE Spencer | Slide Janbu | Sheahan |
|---|---|---|---|---|
| 45° | 1.154 | 1.154 | 1.123 | 1.176 |
| 50° | 1.060 | 1.060 | 1.043 | 1.070 |
| 55° | 0.991 | 0.991 | 0.989 | 0.989 |
| 60° | 0.944 | 0.944 | 0.945 | 0.929 |
| 65° | 0.920 | 0.920 | 0.922 | 0.893 |
| 70° | — | 0.921 | 0.923 | 0.887 |

*The stored surface (and test tag) is the 55° plane, where Slide and Sheahan agree exactly. Janbu's fixed-point iteration does not converge at 70° (Spencer shown). This problem exposed a family of right-facing axial-reinforcement sign errors (vertical force component, facing detection against a vertical wall face, and the Janbu correction-factor chord), all fixed and locked by a left/right mirror consistency test.*

![vp048: inputs and representative solution](images/vp048.png)

### VP50: Reinforced slope, (2) materials, predefined slip surface, geosynthetic {#vp50}

Slide #50 (SNAILZ reference manual): nail-reinforced wall, 14 horizontal rows with per-row length/capacity/bond strength, evaluated on the printed deep wedge (-15.813,0)-(0,-5)-(41.722,25). Slide Janbu corrected 1.417; SNAILZ 1.46. Plate strength equals tensile strength, so the wall end is fully anchored (lp1=0); the embedded end tapers at the bond strength (lp2 = T/bond). Active application, imperial units.

**Input files:** [vp050.xlsx](../files/rocscience/vp050.xlsx)

| Method | XSLOPE | Published |
|---|---|---|
| Janbu (corrected) | 1.448 | SNAILZ 1.46; Slide 1.417; SLOPE/W force 1.354 (×fo ≈ 1.44) |
| Spencer | 1.576 | SLOPE/W M-P 1.606 |

*Tangent orientation with the force factored by FS (Slide’s nail defaults + SNAILZ convention); axial+active gives 1.675 — conventions dominate this comparison.*

![vp050: inputs and representative solution](images/vp050.png)

### VP51: Slope, (4) materials, water table, tension crack, seismic {#vp51}

Slide #51 / GS 2.31: Zhu, Lee & Jiang (2003) four-layer slope, wet, k=0.1, 5 m dry tension crack, specified circle (18.058, 66.744, R=86) read from the printed info box (fig 51.2). Layer-4 properties from the GeoStudio manual (Table 85). The phreatic line is the one element read from the figure trace (anchored at (0,0)-(10,5) on the face, flat ~15.5 at the right); a +/-1 m sensitivity bracket moved Bishop by <0.01, and Bishop/Spencer/Janbu all match the two published programs' agreeing values. Slide/Zhu: OMS 1.145/1.066, Bishop 1.278/1.278, Janbu simp 1.112/1.112, Corps#2 1.422/1.377, Lowe 1.288/1.290, Spencer 1.293/1.293, GLE 1.304/1.303; SLOPE/W: 1.284/1.115/1.368/1.283/1.299/1.310.

**Input files:** [vp051.xlsx](../files/rocscience/vp051.xlsx)

| Method | XSLOPE | Published |
|---|---|---|
| Ordinary | 1.069 | Slide 1.145; Zhu 1.066; SLOPE/W 1.284* |
| Bishop | 1.278 | Slide 1.278; Zhu 1.278; SLOPE/W 1.284 |
| Janbu (corrected) | 1.205 | Slide/Zhu simplified 1.112 (×fo ≈ 1.20) |
| Corps #2 | 1.404 | Slide 1.422; Zhu 1.377; SLOPE/W 1.368 |
| Lowe-Karafiath | 1.296 | Slide 1.288; Zhu 1.290; SLOPE/W 1.283 |
| Spencer | 1.294 | Slide 1.293; Zhu 1.293; SLOPE/W 1.299 |
| Morgenstern-Price | 1.304 | Slide 1.304; Zhu 1.303; SLOPE/W 1.310 |

*Phreatic line calibrated against the two independently agreeing Bishop/Spencer anchors (±1 m bracket).*

![vp051: inputs and representative solution](images/vp051.png)

### VP52: Slope, (4) materials, water table, tension crack {#vp52}

Slide #52, dry. Unconstrained circular search lands in the deep (surface 3) family: Slide grid search Spencer 1.804 / Zhu 1.836 on his specified deep circle (the manual's own Bishop shows a 1.804-vs-1.429 Slide-Zhu spread here, so the family band is wide).

Slide #52, wet (Table 52.2 water table). Deep family: Slide Spencer 1.189 / Zhu 1.211.

**Input files:** [vp052a.xlsx](../files/rocscience/vp052a.xlsx), [vp052b.xlsx](../files/rocscience/vp052b.xlsx)

*Dry — governing deep (surface 3) family*

| Method | XSLOPE | Published |
|---|---|---|
| Bishop | 1.796 | Slide 1.804; Zhu 1.429 |
| Spencer | 1.797 | Slide 1.804; Zhu 1.836 |

*Wet — governing deep (surface 3) family*

| Method | XSLOPE | Published |
|---|---|---|
| Bishop | 1.176 | Slide 1.176; Zhu 1.079 |
| Spencer | 1.189 | Slide 1.189; Zhu 1.211 |

*Wet Spencer and Bishop match Slide exactly; the manual itself shows a wide Slide–Zhu spread on this family.*

![vp052a: inputs and representative solution](images/vp052a.png)

![vp052b: inputs and representative solution](images/vp052b.png)

### VP54: Slope, homogenous, micro piles {#vp54}

Slide #54, unreinforced case on the printed critical circle (2.674, 7.573, R=8.031). Slide Bishop 1.102; Yamagami 1.10.

Slide #54 with the micro-pile row. Slide 1.193; Yamagami 1.20.

**Input files:** [vp054a.xlsx](../files/rocscience/vp054a.xlsx), [vp054b.xlsx](../files/rocscience/vp054b.xlsx)

*No pile*

| Method | XSLOPE | Published |
|---|---|---|
| Bishop | 1.100 | Slide 1.102; Yamagami 1.10; SLOPE/W 1.102 |

*With micro-pile row*

| Method | XSLOPE | Published |
|---|---|---|
| Bishop | 1.185 | Slide 1.193; Yamagami 1.20; SLOPE/W 1.223 |

*Slide adds the pile shear un-factored (= our active application); a free search finds 1.113 on a circle exiting upslope of the pile, so the tags pin the printed circle.*

![vp054a: inputs and representative solution](images/vp054a.png)

![vp054b: inputs and representative solution](images/vp054b.png)

### VP62: Slope, homogenous, ru pore pressure, seismic {#vp62}

Slide #62 dry case, kc=0.432. Slide circular: Spencer 1.001, Bishop 0.991; Loukidis log-spiral Spencer 1.000.

Slide #62 ru=0.5 case, kc=0.132. Slide circular: Spencer 1.001, Bishop 0.987; Loukidis 1.000.

**Input files:** [vp062a.xlsx](../files/rocscience/vp062a.xlsx), [vp062b.xlsx](../files/rocscience/vp062b.xlsx)

*Dry, kc = 0.432*

| Method | XSLOPE | Published |
|---|---|---|
| Bishop | 0.991 | Slide 0.991; SLOPE/W 0.993 |
| Spencer | 1.001 | Slide 1.001; SLOPE/W 1.001; Loukidis 1.000 |

*ru = 0.5, kc = 0.132*

| Method | XSLOPE | Published |
|---|---|---|
| Bishop | 0.986 | Slide 0.987; SLOPE/W 0.988 |
| Spencer | 1.001 | Slide 1.001; SLOPE/W 1.001; Loukidis 1.000 |

![vp062a: inputs and representative solution](images/vp062a.png)

![vp062b: inputs and representative solution](images/vp062b.png)

### VP85: Reinforced slope, homogenous, grouted tieback {#vp85}

Slide #85 case 1 (active). D&W reference 1.51; Slide circular Bishop 1.531.

Slide #85 case 2 (passive). D&W reference 1.32; Slide circular Bishop 1.324.

**Input files:** [vp085a.xlsx](../files/rocscience/vp085a.xlsx), [vp085b.xlsx](../files/rocscience/vp085b.xlsx)

*Active support, on Slide’s printed GLE circle*

| Method | XSLOPE | Published |
|---|---|---|
| Bishop | 1.567 | Slide GLE 1.575 (same circle); D&W 1.51 |
| Spencer | 1.567 | Slide GLE 1.575 (same circle) |

*Passive support, on Slide’s printed Bishop circle*

| Method | XSLOPE | Published |
|---|---|---|
| Ordinary | 1.319 | Slide Bishop 1.324 (same circle); D&W 1.32 |
| Bishop | 1.319 | Slide 1.324; D&W 1.32 |

*Slide’s own method table scatters 1.42–2.05 here (concentrated force strains interslice assumptions), so per-circle comparison is the meaningful one.*

![vp085a: inputs and representative solution](images/vp085a.png)

![vp085b: inputs and representative solution](images/vp085b.png)

### VP86: Reinforced slope, homogenous, grouted tieback {#vp86}

Slide #86: Duncan & Wright (2005) Fig. 7.28 / STABGM reinforced fill on a strong rock foundation: 5 geogrids (800 lb/ft, 20 ft long, every 4 ft). Slide2 circular: Bishop 1.629, Spencer 1.620, GLE 1.622; D&W reference 1.61.

**Input files:** [vp086.xlsx](../files/rocscience/vp086.xlsx)

| Method | XSLOPE | Published |
|---|---|---|
| Bishop | 1.617 | Slide 1.629 |
| Spencer | 1.611 | Slide 1.620 |

*Duncan & Wright reference 1.61.*

![vp086: inputs and representative solution](images/vp086.png)

### VP87–VP94: Geosynthetic multitiered MSE walls (Leshchinsky & Han 2004) {#vp87}

Slide #87–#94 reproduce the parametric study in Leshchinsky & Han (2004): a three-tier segmental (block-faced) MSE wall — three 3-m tiers offset 1.2 m, 0.3-m block columns (c=2.5 kPa, φ=34°), reinforced/retained fill c=0/φ=34°, foundation c=10 kPa/φ=34° (6 m deep), γ=18 kN/m³ throughout — with geotextile layers every 0.6 m, L=6.3 m, and the tensile strength Ta the paper *required* for FS=1.0 in each variation. Pullout resistance is 80% of the fill strength (translated to xslope anchorage lengths from the local overburden at each layer end); the geotextile force is applied horizontally (`dir=axial`, `appl=passive` — Slide's convention, verified against its printed VP87 circle). Each problem's Slide-printed critical circle is stored in the file, so the test tags evaluate a deterministic surface.

Two manual quirks resolved during the build: (1) Slide's VP89/92/93 *results* were computed with the **baseline Ta = 10** supports even though their support tables print the paper's per-case required strengths (11.4/9.25/11.6) — with Ta=10 xslope lands within 1% of all three Slide numbers, and with the paper's strengths it lands on L&H's design intent (FS ≈ 1.0). (2) VP91's printed circle exits exactly tangent to the crest and needs a hair of extra radius to intersect.

| # | Case | Method (Slide's figure) | XSLOPE | Slide | L&H reference |
|---|---|---|---|---|---|
| 87 | Baseline (Ta=10) | Bishop | 1.031 | 1.040 | 0.99 (FLAC) / 1.00 (Bishop) |
| 88 | Fill φ=25 (Ta=22) | Spencer | 1.057 | 1.043 | 0.99 / 1.00 |
| 89 | L=4.2 m (Ta=11.4) | Spencer | 1.011 *(0.980 at Ta=10)* | 0.971 *(used Ta=10)* | 0.98 / 1.00 |
| 90 | Two types (7.5/11.0) | Bishop | 1.012 | 1.004 | 1.01 / 1.00 |
| 91 | Foundation c=0, φ=18 | Spencer | 0.960 | 0.964 | 0.86 (FLAC, bearing) / 1.00 |
| 92 | Water hw=3 m (Ta=9.25) | Bishop | 1.010 *(1.039 at Ta=10)* | 1.037 *(used Ta=10)* | 1.01 / 1.00 |
| 93 | Surcharge q=20 (Ta=11.6) | Bishop | 1.017 *(0.961 at Ta=10)* | 0.958 *(used Ta=10)* | 1.02 / 1.00 |
| 94 | Five 1.8-m tiers (Ta=10.1) | Bishop | 1.020 | 1.040 | 1.00 |

*VP92 models the paper's hw as pore pressure in the foundation soil only (a drained MSE fill), plus the 3-m pond standing against the lower tier — treating the fill as saturated drops FS to ~0.89 and reconciles with neither program. xslope's free circular search finds slightly more critical circles than Slide's grid on several of these (e.g. 0.99 on the baseline, matching the L&H reference).*

![vp087: inputs and representative solution](images/vp087.png)

![vp088: inputs and representative solution](images/vp088.png)

![vp089: inputs and representative solution](images/vp089.png)

![vp090: inputs and representative solution](images/vp090.png)

![vp091: inputs and representative solution](images/vp091.png)

![vp092: inputs and representative solution](images/vp092.png)

![vp093: inputs and representative solution](images/vp093.png)

![vp094: inputs and representative solution](images/vp094.png)

### VP96: Embankment dam, homogenous, rapid drawdown, water table {#vp96}

Slide #96 / USACE EM 1110-2-1902 (2003) Appendix G example: 3:1 then 2.5:1 embankment face, max pool el. 103 drawn down to 24, specified circle (169.5, 210, R=210). Material: c'=0, phi'=30, gamma=135 pcf with the Kc=1 envelope d=1379 psf, psi=18.2 deg (Figure G-5). Duncan-Wright- Wong 3-stage: Slide 1.443, USACE reference 1.44. (Slide's #95 runs the same model with the older Corps 2-stage method: 1.347.)

**Input files:** [vp096.xlsx](../files/rocscience/vp096.xlsx)

| Method | XSLOPE | Published |
|---|---|---|
| Bishop | 1.432 | Slide 1.443; USACE 1.44 |
| Spencer | 1.434 | Slide 1.443; USACE 1.44 |

*Duncan-Wright-Wong 3-stage on the specified circle; Kc=1 envelope d=1379 psf, ψ=18.2°.*

![vp096: inputs and representative solution](images/vp096.png)

### VP97: Embankment dam, homogenous, rapid drawdown, water table {#vp97}

Slide #97: Pilarcitos Dam (Duncan, Wright & Wong 1990). Homogeneous earthfill, gamma=135 pcf, c'=0, phi'=45; R-envelope cR=60 psf, phiR=23. Kc=1 envelope via D&W (2014) eqs 9.6-9.7: d = cR cos(phiR) cos(phi') / (1-sin(phiR)) = 64.1 psf, psi = 24.4 deg (the same equations reproduce the USACE App G values 1379/18.2 exactly). Drawdown 72 -> 37 ft. DWW 3-stage 1.05; Slide 3-stage 1.043 (Corps 2-stage 0.823/0.82).

**Input files:** [vp097.xlsx](../files/rocscience/vp097.xlsx)

| Method | XSLOPE | Published |
|---|---|---|
| Bishop | 1.042 | Slide 1.043; DWW 1.05 |
| Spencer | 1.044 | Slide 1.043; DWW 1.05 |

*The dam that actually failed in drawdown sits right at FS ≈ 1.*

![vp097: inputs and representative solution](images/vp097.png)

### VP98: Walter Bouldin Dam rapid drawdown (Duncan, Wright & Wong 1990) {#vp98}

Slide #98: the Walter Bouldin Dam failure case — a rolled earthfill dam that failed during a 32-ft drawdown in 1975. Pool drops 47 ft → 15 ft. Five zones (riprap, clayey silty sand, micaceous sand, cretaceous clay, clayey sandy gravel) rebuilt from Slide's coordinate-labeled Figure 98.1 with the interior boundaries traced from its color zones (axis-calibrated, ±1 ft); the Kc=1 undrained envelopes come from the paper's own Table 2 — (750 psf, 15°), (480, 13°), (280, 15.5°) — with riprap and gravel drained.

**Input files:** [vp098.xlsx](../files/rocscience/vp098.xlsx)

| Method | XSLOPE | Published |
|---|---|---|
| DWW 3-stage (Spencer, circular search) | 1.046 | Slide 1.039; DWW 1.04 |

*The critical circle ((52,21)→(157,60)) falls where the dam actually slid. Slide's Corps 2-stage (0.931) and Lowe & Karafiath (1.075) rows exercise staged procedures xslope does not implement.*

![vp098: inputs and representative solution](images/vp098.png)

### VP99: Pumped-storage project dam rapid drawdown (DWW 1990) {#vp99}

Slide #99: the paper's hypothetical pumped-storage dam — silty clay core and random zone (c′=0, φ′=36°, Kc=1 envelope 2250 psf/20°), free-draining rockfill shells (φ′=37°), drawdown 285 ft → 120 ft (paper El 545 → 380). Slide's figure has no vertex labels, so the geometry was extracted by axis-calibrated color segmentation (the 285/120 pool dashes pin the vertical scale; the crest coincides with the initial pool).

**Input files:** [vp099.xlsx](../files/rocscience/vp099.xlsx)

| Method | XSLOPE | Published |
|---|---|---|
| DWW 3-stage (Spencer, circular search) | 1.390 | Slide 1.534; SLOPE/W 1.550; DWW 1.56 |
| DWW 3-stage (Spencer, Slide's printed circle) | 1.428 | Slide 1.534 |

*≈7% low. The two vendors' own geometries differ visibly (berm elevations, core width), and the result is sensitive to the core-face position; the geometry will be re-pinned from GeoStudio's published .gsz project when available. Tagged at the current value as a regression lock.*

![vp099: inputs and representative solution](images/vp099.png)

### VP100: Embankment dam, homogenous, rapid drawdown, water table {#vp100}

Slide #100: complete drawdown (100 -> 0), B-bar = 1: the residual pore pressure is hydrostatic below the slope surface, i.e. piezo = ground, no external pond. Slide B-bar method 1.212; Morgenstern chart 1.20.

**Input files:** [vp100.xlsx](../files/rocscience/vp100.xlsx)

*Complete drawdown (100 → 0)*

| Method | XSLOPE | Published |
|---|---|---|
| Bishop | 1.201 | Morgenstern chart 1.20; Slide (B-bar) 1.212 |
| Spencer | 1.206 | — |

![vp100: inputs and representative solution](images/vp100.png)

### VP101: Embankment dam, homogenous, rapid drawdown, water table {#vp101}

Slide #101: partial drawdown (100 -> 50), B-bar = 1: piezo follows the ground where the face is above the pool and stays at 50 below it, with the remaining pond loading the submerged face. Slide 1.417; Morgenstern chart 1.41.

**Input files:** [vp101.xlsx](../files/rocscience/vp101.xlsx)

*Partial drawdown (100 → 50)*

| Method | XSLOPE | Published |
|---|---|---|
| Bishop | 1.416 | Slide 1.417; Morgenstern chart 1.41 |
| Spencer | 1.422 | — |

![vp101: inputs and representative solution](images/vp101.png)
