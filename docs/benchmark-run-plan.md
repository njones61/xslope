# xslope Verification & Validation — Benchmark Run Plan

Purpose: a complete, build-and-run summary of every benchmark case selected for
the Verification and Validation section of the JGGE Technical Note. Copy this
file into the xslope working folder, build an input file per case, run it, and
record the xslope result in the "Capture" blanks. Then compute the percent
difference and report values back for the manuscript tables.

**Verification philosophy (two tiers per mode):**
1. *Analytical anchor* — match a closed-form solution where one exists.
2. *Established-code / multiply-published cross-check* — agreement with tools
   practitioners already trust on realistic problems no closed form covers.

**Ground rules:**
- Pull exact node coordinates and properties from the cited source. Do not
  estimate geometry from figures.
- ACADS "correct" answers are consensus/averaged values: report agreement as
  "within the accepted band," not against a single true value.
- For SSRM: default **non_convergence** criterion (bisection on true
  viscoplastic equilibrium — CHECON + plastic-settled); for submerged/
  reservoir problems use **displacement_increase** (the displacement-vs-F
  upturn, the evidence G&L themselves publish, Figs 2/18) with
  tension_cutoff. The literal iteration-ceiling form of non-convergence is
  not transferable between codes (see docs/fem/overview.md, "Choosing a
  Failure Criterion").
- Record the units used per case (the sources mix SI and English).

---

## Case-to-table map

| Case ID | Mode | Feeds manuscript table | Tier |
|---|---|---|---|
| LEM-1 | Limit equilibrium (6 methods), circular | `tab:lem` | code/consensus cross-check |
| LEM-2 | Limit equilibrium, non-circular (weak layer) | `tab:lem` companion | code/consensus cross-check |
| LEM-2b | Limit equilibrium (6 methods), circular | `tab:lem` companion | published literature |
| LEM-3 (optional) | Limit equilibrium, method ordering | supplemental / text | published literature |
| SEEP-1 | FE seepage, confined analytical | `tab:seep` | analytical anchor (exact) |
| SEEP-1c | FE seepage, sheetpile cutoff | `tab:seep` | analytical anchor (Pavlovsky) |
| SEEP-1b | FE seepage, unconfined free surface | sample / text | analytical (Kozeny, ±4%) |
| SEEP-2 | FE seepage, code cross-check | `tab:seep` | established code (SEEP2D) |
| SSRM-1 | FE slope stability | `tab:ssrm` | published (G&L Ex. 1) |
| SSRM-2 | FE slope stability | `tab:ssrm` | published (G&L Ex. 6) |
| SSRM-3 (optional) | FE slope stability | supplemental | published (Rocscience tabulated) |

---

## Mode 1 — Limit Equilibrium (six methods)

Six methods per case: OMS, Bishop's Simplified, Simplified Janbu, Corps of
Engineers, Lowe & Karafiath, Spencer. (OMS and Simplified Janbu are reported for
completeness only; they are legacy/educational, not recommended for design.)

### LEM-1 — ACADS simple homogeneous slope (headline LEM table)
- **Source:** GeoStudio SLOPE/W Verification Manual (Oct 2022), ACADS suite;
  original Donald & Giam (1989), reviewed Giam & Donald (1992).
  PDF: https://files.seequent.com/PDFs/Geostudio-Slope%20Stability%20Verification%20Manual-Oct2022.pdf
- **Geometry/materials:** pull exact slope geometry and single-soil properties
  (c', phi', gamma) from the manual's ACADS simple-slope section. Circular search.
- **Reference FOS:** ACADS accepted band (~1.0 region); record the manual's
  per-method reference values.
- **What to run:** all six methods on the same critical circle search.
- **Capture:**

  | Method | xslope FOS | Reference FOS | Diff (%) |
  |---|---|---|---|
  | Ordinary (OMS) | 0.942 | 1.00 | -5.8 |
  | Bishop's Simplified | 0.987 | 1.00 | -1.3 |
  | Simplified Janbu | 0.992 | 1.00 | -0.8 |
  | Corps of Engineers | 0.991 | 1.00 | -0.9 |
  | Lowe & Karafiath | 0.987 | 1.00 | -1.3 |
  | Spencer | 0.986 | 1.00 | -1.4 |

  *xslope FOS = automated critical-circle search, 50 slices, each method searched
  independently. Corps of Engineers uses variant 2 (per-slice ground-parallel
  side forces; see note under LEM-2).*

### LEM-2 — ACADS weak-layer slope (non-circular)
- **Source:** GeoStudio SLOPE/W Verification Manual (Oct 2022), ACADS weak-layer
  case (sec. 2.7); ACADS suite, Donald & Giam (1989).
- **Geometry/materials:** 2:1 slope with a thin low-strength interlayer. Soil 1
  c'=28.5, phi'=20, gamma=18.84; Weak Layer c'=0, phi'=10, gamma=18.84; piezo
  line at the base of the weak band. Coordinates from `benchmarks/build_lem.py`
  (`build_acads_weak_layer`). This is the **non-circular search** test — the
  critical surface runs along the weak layer with a back scarp to the crest.
- **Reference FOS:** ACADS accepted band ≈ 1.26 (force/complete-equilibrium
  methods).
- **What to run:** non-circular critical-surface search; report the methods that
  apply to non-circular surfaces (Spencer, Corps of Engineers, Lowe & Karafiath,
  Simplified Janbu).
- **Capture:**

  | Method | xslope FOS | Reference FOS | Diff (%) |
  |---|---|---|---|
  | Spencer | 1.279 | ~1.26 | +1.5 |
  | Corps of Engineers | 1.355 | ~1.26 | +7.6 |
  | Lowe & Karafiath | 1.268 | ~1.26 | +0.6 |
  | Simplified Janbu | 1.278 | ~1.26 | +1.4 |

  *Corps of Engineers uses **variant 2** (side forces parallel to the ground
  surface at each slice top — the standard "Corps of Engineers #2" convention;
  USACE EM 1110-2-1902). The single-inclination "#1" convention (crest-to-toe
  chord) is fragile under an unconstrained non-circular search: it lets the
  search dive to a surface with a near-vertical toe segment where the fixed
  inclination collapses the force balance (Corps FS → 1.05 while Spencer rates
  that same surface 1.53). Variant 2 is xslope's default for this reason. Corps
  reads slightly high here (+7.6%), consistent with ground-parallel side forces
  on this geometry.*

### LEM-2b — Arai & Tagyo (1985) homogeneous slope (circular)
- **Source:** Arai & Tagyo (1985), Soils and Foundations 25(1):43-51; SLOPE/W
  manual sec. 2.11; re-published in Greco (1996), Malkawi et al. (2001), Kim et
  al. (2002).
- **Geometry/materials:** homogeneous 1.5:1 slope, 20 m high; single soil
  c=41.65, phi=15.0, gamma=18.82 (total stress). Coordinates from
  `benchmarks/build_lem.py` (`build_arai_tagyo`). Circular search.
- **Reference FOS:** published Spencer/Bishop ≈ 1.451 (records vary 1.40–1.45
  across re-publications).
- **What to run:** circular critical-surface search, all six methods.
- **Capture:**

  | Method | xslope FOS | Reference FOS | Diff (%) |
  |---|---|---|---|
  | Ordinary (OMS) | 1.344 | 1.451 | -7.4 |
  | Bishop's Simplified | 1.404 | 1.451 | -3.2 |
  | Simplified Janbu | 1.441 | 1.451 | -0.7 |
  | Corps of Engineers | 1.477 | 1.451 | +1.8 |
  | Lowe & Karafiath | 1.439 | 1.451 | -0.8 |
  | Spencer | 1.402 | 1.451 | -3.4 |

### LEM-3 (optional) — Fredlund & Krahn (1977) method ordering
- **Source:** Fredlund & Krahn (1977), Can. Geotech. J. 14(3):429-439.
- **Use:** confirms expected inter-method behavior (OMS low; Bishop and Spencer
  converge closely). Good for one sentence of interpretation; may live in text or
  Supplemental rather than its own table.
- **Capture:** per-method FOS vs the paper's tabulated values.

---

## Mode 2 — FE Seepage

### SEEP-1 — Confined radial flow (analytical anchor, saturated)
- **Source:** steady confined (saturated) flow; exact solution of Laplace's
  equation in polar coordinates (radial flow between two concentric equipotential
  arcs). Standard analytical solution (e.g. Bear, *Dynamics of Fluids in Porous
  Media*).
- **Problem:** a quarter-annulus confined aquifer (impervious base/cap implied by
  plane strain); inner arc r1=10 at head h1=30, outer arc r2=30 at head h2=10, the
  two straight radial edges are no-flow streamlines. Built via the **polygon**
  input (`benchmarks/build_seep.py::build_confined_radial`). Chosen as the
  analytical anchor instead of an unconfined dam because it is *exact* with no
  free-surface tangency: the discharge integral for a toe-drain Kozeny dam is
  sensitive (~±4%) to where the drain-start boundary lands relative to the point
  where the free surface meets the base.
- **Closed form:**
  - head:  `h(r) = h1 + (h2 - h1)*ln(r/r1)/ln(r2/r1)`
  - discharge (sector angle alpha):  `q = k*alpha*(h1 - h2)/ln(r2/r1)`
    = `(pi/2)*20/ln 3` = **28.596** (k=1).
- **Capture:**

  | Quantity | xslope | Analytical | Diff (%) |
  |---|---|---|---|
  | Discharge q (tri6) | 28.5961 | 28.596 | <0.01 |
  | Discharge q (tri3) | 28.6001 | 28.596 | +0.01 |
  | Max nodal head error (tri6) | 0.004 | 0 (exact) | 0.02 (of 20-unit drop) |

  *Mesh-converged: tri6 q = 28.5961 at both 2k and 6k nodes; quad8 identical. The
  only error is faceting of the curved equipotential arcs by the polygon boundary.*
  *Note: this case initially produced a singular system with tri3 elements, which
  exposed a real bug — gmsh inherits the winding of the input polygon ring, and
  the linear-triangle assembly silently skips clockwise elements (`area <= 0`).
  Fixed by normalizing all 2D elements to CCW in `build_mesh_from_polygons`
  (`xslope/mesh.py::ensure_ccw_elements`); CW polygon input now meshes correctly
  for every element type.*

### SEEP-1c — Partially penetrating sheetpile (confined, exact)
- **Source:** Pavlovsky's conformal-mapping solution for a single cutoff wall in a
  homogeneous confined stratum of finite thickness; Harr, *Groundwater and
  Seepage* (1962); Polubarinova-Kochina (1962).
- **Problem:** wall of depth s in a stratum of thickness T = 20, head loss
  H = 10 across it (25/15), k = 1, horizontal extent 4T each side (truncation
  error ~ exp(-pi*L/T) ~ 3e-6). Wall modeled with the V-notch crack idiom
  (w/T = 0.005), as in the clay-blanket sample. Built by
  `benchmarks/build_seep.py::build_sheetpile` for s/T = 0.5 and 0.75.
- **Closed form:** `q = k*H*K(lam')/(2*K(lam))`, `lam = sin(pi*s/(2T))`. At
  s/T = 1/2 the modulus is self-dual, so **q = k*H/2 exactly** (no elliptic
  integrals in the headline number). Second exact check: by antisymmetry the
  head on the wall plane below the tip is exactly (H1+H2)/2.
- **Capture** (tri6):

  | Case | xslope q | Exact q | Diff (%) | Head below tip |
  |---|---|---|---|---|
  | s/T = 0.5, 15k nodes | 5.0288 | 5.0000 | +0.58 | 20.0000 (exact) |
  | s/T = 0.5, 59k nodes | 5.0101 | 5.0000 | +0.20 | 20.0000 (exact) |
  | s/T = 0.75, 15k nodes | 3.4246 | 3.4032 | +0.63 | 20.0000 (exact) |
  | s/T = 0.75, 59k nodes | 3.4122 | 3.4032 | +0.27 | 20.0000 (exact) |

  *Error halves with mesh refinement (first order, set by the r^-1/2 tip
  singularity) and converges to the exact value from above; the residual is the
  finite notch width plus tip discretization. This case complements the radial
  anchor: it exercises a thin internal barrier with a singular tip — the
  geometry cutoff-wall analyses actually involve.*

  *Formula independently verified: a from-scratch finite-difference solution of
  the same BVP (5-point Laplacian on the antisymmetric half-domain, flux through
  the throat below the tip, sqrt(h)-Richardson extrapolation) reproduces the
  K(lam')/2K(lam) form factor to ~0.4-0.5% at s/T = 0.25, 0.5, and 0.75 — so the
  xslope comparison is anchored to a confirmed closed form, not a recalled one.*

### SEEP-1b — Kozeny dam (free-surface sample, retained for scrutiny)
- **Source:** Kozeny (1931) basic parabola + Casagrande entrance correction.
- **Problem:** homogeneous earth dam, horizontal toe drain, impervious base. Built
  with the upstream face as the **exact confocal-parabola equipotential** of the
  Kozeny flow net (`build_seep.py::build_kozeny_dam`), which removes the entrance
  error of a straight face. Kept as a seepage *sample* (see `docs/seep/samples.md`)
  and a free-surface check, not the headline analytical anchor.
- **Findings:** with the confocal face, xslope *brackets* the closed-form
  q = k*y0 = 4.0 within ±4% — drain at the focus gives +3.8%, drain at the free-
  surface/base vertex gives -4% — the spread being the drain-start tangency
  sensitivity, not a solver bias. The FE free surface tracks the analytical
  parabola `x = 90 - y^2/8` to ~1-2 length units. (A straight 2.5:1 face instead
  gives a converged -7%, the classic basic-parabola entrance over-prediction.)

### SEEP-2 — SEEP2D code cross-check (Johnson Reservoir)
- **Source:** SEEP2D (USACE/WES, Fred Tracy), via GMS. Author has access.
- **Problem:** the **Johnson Reservoir** zoned earth dam, the same geometry used
  in the worked example (permeable shell, low-permeability core, foundation;
  reservoir at elev. 160, tailwater at elev. 100). Run SEEP2D and xslope on the
  **same** geometry and BCs. (xslope can import SEEP2D input files directly, which
  simplifies matching.) This validates the very seepage solution that drives the
  worked-example stability analyses, so the demo rests on a verified field rather
  than an unchecked one.
- **What to compare:** nodal heads, total discharge (xslope gave 1.94 ft^3/day per
  unit length), phreatic surface location.
- **Capture:**

  | Quantity | xslope | SEEP2D | Diff (%) |
  |---|---|---|---|
  | Total discharge q | 1.9575 | 1.9603 | -0.14 |
  | Nodal heads (all 2604 nodes) | RMS dh = 0.105 ft | (60-ft head range) | 0.18 |

  *DONE — and stronger than originally planned: rather than comparing against a
  separate GMS run on a different mesh, the original USACE/WES SEEP2D Fortran
  code (ref_docs/ref_docs_seep) was compiled with gfortran and run on the EXACT
  same tri3 mesh topology, BCs, and material parameters, exported from xslope
  via `benchmarks/run_seep2d_compare.py` (RCM-renumbered for SEEP2D's banded
  solver). Identical-mesh comparison: q within 0.14%, nodal heads RMS 0.105 ft
  over a 60-ft range; max |dh| ~ 2 ft localized near the free surface, where
  the codes' unsaturated-kr details differ.*

---

## Mode 3 — FE Slope Stability (SSRM)

Elastic-perfectly-plastic Mohr-Coulomb, plane strain, viscoplastic algorithm
(Smith & Griffiths 4-component formulation with Lode-angle corner handling).
Failure = displacement catastrophe: sharp upturn of the converged displacement
vs the strength-reduction factor F (G&L Figs 2/18).

### SSRM-1 — Griffiths & Lane (1999) Example 1 (homogeneous)  [DONE]
- **Source:** Griffiths & Lane (1999), Geotechnique 49(3):387-403, Example 1
  (paper in `ref_docs/ref_docs_fem/`); Rocscience re-run:
  https://www.rocscience.com/assets/resources/learning/papers/Application-of-the-Finite-Element-Method-to-Slope-Stability.pdf
- **Geometry/materials:** homogeneous 2:1 slope (26.57 deg), D = 1 (no
  foundation), crest width 1.2H, slope run 2H; phi' = 20 deg, c'/(gamma*H) =
  0.05, nominal E and nu; left rollers, fixed base; psi = 0. Existing model:
  `docs/fem/files/xslope_griffiths1.xlsx` (H = 50, gamma = 125, c' =
  gamma*H*0.05 — same dimensionless group, English units).
- **Reference FOS:** G&L FE = 1.4 (their algorithm converges at 1.35 and
  fails at 1.40, Table 2); Bishop & Morgenstern (1960) chart = 1.380.
- **Capture (final, tol 0.01):** xslope FOS = **1.36** under the strict
  true-equilibrium criterion (CHECON + plastic-settled); the
  displacement-vs-F upturn sits at **F ~ 1.40**, matching G&L's reported
  1.4 (their tolerant convergence corresponds to the upturn) and within
  -1.5% of the B&M chart 1.380. All readings agree within +/-3%. Sweep:
  max|u| flat 0.17-0.21 through F = 1.35, then 0.29 / 0.65 / 1.16 / 2.39
  at F = 1.40 / 1.45 / 1.5 / 1.6.

### SSRM-2 — Griffiths & Lane (1999) Example 6 (earth dam, free surface)
- **Source:** G&L (1999) Example 6, Figs 16-21; original cross-section from
  Torres & Coffman (1997); Rocscience re-run appendix.
- **Geometry/materials (Fig. 16, metres):** two-sided embankment on a 7.3-m
  foundation layer; dam height H = 21.3 above foundation level; crest width
  7.3; upstream face 18 deg, downstream face 23 deg; dam base 124.4 with
  33.5 aprons each side. Homogeneous: c' = 13.8 kPa, phi' = 37 deg,
  gamma = 18.2 kN/m3 above and below the water table. Free surface from
  reservoir level (17.1 above foundation) at the upstream face down to the
  downstream toe at foundation level; pore pressure = gamma_w x vertical
  depth below free surface (G&L's simplification = xslope piezo line).
  Reservoir water load applied as normal stress on the submerged upstream
  face (xslope dload). Rollers on the foundation ends, fixed base.
- **Reference FOS (paper, Figs 18-19):** **with free surface (reservoir
  full) ~1.9; no free surface (before filling) ~2.4** — both confirmed by
  the companion limit-equilibrium runs (1.90 / 2.42). (Earlier note had
  these two swapped; failure is on the steeper downstream face, so the wet
  dam is the weaker case.)
- **What to run:** SSRM for both cases.
- **Capture** (quad8, displacement-catastrophe criterion):

  | Case | xslope FOS | G&L FOS | Diff (%) |
  |---|---|---|---|
  | Example 6 (full reservoir, free surface) | 2.10 | ~1.9 | +11 |
  | Example 6 (no free surface) | 2.45 | ~2.4 | +2 |

  *Notes: (1) The wet result is mesh-converged — identical FS (~2.1) at
  target_size 3.5 and 2.0. (2) The input model is independently validated:
  xslope's Spencer on the same file gives FOS = 1.915 vs G&L's limit-
  equilibrium companion value 1.90 (+0.8%), with the critical circle on the
  downstream face as published. (3) The relative reservoir effect matches:
  wet/dry = 0.82 (xslope) vs 0.79 (G&L). (4) The +6-10% FEM-high offset
  mirrors the pre-existing internal LEM-vs-SSRM gap on the Johnson
  Reservoir worked example (Spencer 1.26 vs SSRM ~1.40, ~+11%) — a
  consistent xslope-FEM-vs-LEM characteristic to discuss in the manuscript,
  not a per-case anomaly. (5) G&L's published FE values were obtained with
  their non-convergence criterion inside their code's numerical regime,
  which is not fully documented in the paper and was found to be
  non-transferable (see docs/fem/overview.md).*

### SSRM-3 (optional) — Rocscience tabulated cross-checks
Clean tabulated geometry + FOS, easy to reproduce; good as supplemental
robustness checks if room allows.
- **Their Ex. 1** (= Slide verification #1): homogeneous + foundation;
  gamma = 20.2, phi = 19.6, c = 3.0; Bishop 0.988 / Spencer 0.987 / GLE 0.987.
- **Their Ex. 2** (= Slide verification #3): three-layer; soils
  (c=0, phi=38), (c=5.3, phi=23), (c=7.2, phi=20), all gamma = 19.5; FOS ~1.38-1.41.

---

## Open decisions to resolve while running

1. **Headline LEM problem:** RESOLVED — LEM-1 (ACADS simple homogeneous, circular)
   anchors `tab:lem`. The non-circular test is the ACADS weak-layer case (LEM-2);
   Arai & Tagyo is now the homogeneous circular cross-check (LEM-2b). All three
   are built by `benchmarks/build_lem.py` and run by `benchmarks/run_lem.py`.
2. **Seepage cross-check geometry:** RESOLVED — use Johnson Reservoir for SEEP-2,
   the same geometry as the worked example. The SEEP2D cross-check then validates
   the worked-example seepage field directly.
3. **Units per case:** G&L/Arai & Tagyo are SI; Johnson Reservoir is English.
   Keep each case internally consistent and label units in the tables.
4. **LEM-vs-SSRM gap:** the worked example shows Spencer 1.26 vs SSRM 1.40 on
   Johnson Reservoir (~11%). If an internal consistency check is wanted, confirm
   the two methods are locating the same mechanism.

---

## Source quick list

- **SLOPE/W Verification Manual (Oct 2022)** — ACADS suite, Arai & Tagyo, Greco,
  rapid drawdown; dimensioned geometry + reference FOS.
  https://files.seequent.com/PDFs/Geostudio-Slope%20Stability%20Verification%20Manual-Oct2022.pdf
- **Griffiths & Lane (1999)** — inside.mines.edu/~vgriffit (Griffiths' faculty page).
- **Rocscience "Application of the FEM to Slope Stability"** — tabulated G&L re-runs.
  https://www.rocscience.com/assets/resources/learning/papers/Application-of-the-Finite-Element-Method-to-Slope-Stability.pdf
- **Arai & Tagyo (1985)**, Soils and Foundations 25(1):43-51.
- **Kozeny (1931)**; **Casagrande** entrance correction; **Polubarinova-Kochina (1962)**.
- **SEEP2D** (USACE/WES, Tracy) via GMS.
