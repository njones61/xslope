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
- For SSRM, replicate Griffiths & Lane's **non-convergence** failure criterion.
- Record the units used per case (the sources mix SI and English).

---

## Case-to-table map

| Case ID | Mode | Feeds manuscript table | Tier |
|---|---|---|---|
| LEM-1 | Limit equilibrium (6 methods), circular | `tab:lem` | code/consensus cross-check |
| LEM-2 | Limit equilibrium, non-circular (weak layer) | `tab:lem` companion | code/consensus cross-check |
| LEM-2b | Limit equilibrium (6 methods), circular | `tab:lem` companion | published literature |
| LEM-3 (optional) | Limit equilibrium, method ordering | supplemental / text | published literature |
| SEEP-1 | FE seepage, analytical | `tab:seep` | analytical anchor |
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
  | Corps of Engineers | 0.990 | 1.00 | -1.0 |
  | Lowe & Karafiath | 0.986 | 1.00 | -1.4 |
  | Spencer | 0.986 | 1.00 | -1.4 |

  *xslope FOS = automated critical-circle search, 20 slices, each method searched
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
  | Ordinary (OMS) | 1.340 | 1.451 | -7.7 |
  | Bishop's Simplified | 1.403 | 1.451 | -3.3 |
  | Simplified Janbu | 1.436 | 1.451 | -1.0 |
  | Corps of Engineers | 1.474 | 1.451 | +1.6 |
  | Lowe & Karafiath | 1.436 | 1.451 | -1.0 |
  | Spencer | 1.400 | 1.451 | -3.5 |

### LEM-3 (optional) — Fredlund & Krahn (1977) method ordering
- **Source:** Fredlund & Krahn (1977), Can. Geotech. J. 14(3):429-439.
- **Use:** confirms expected inter-method behavior (OMS low; Bishop and Spencer
  converge closely). Good for one sentence of interpretation; may live in text or
  Supplemental rather than its own table.
- **Capture:** per-method FOS vs the paper's tabulated values.

---

## Mode 2 — FE Seepage

### SEEP-1 — Kozeny / Casagrande analytical anchor
- **Source:** Kozeny (1931) basic parabola + Casagrande entrance correction;
  Polubarinova-Kochina (1962) for rigorous free-surface checks.
- **Problem:** homogeneous earth dam with horizontal toe drain (the classic
  unconfined free-surface case with a closed-form phreatic surface and discharge).
- **What to compare:**
  - Discharge `q` (per unit length) vs Kozeny closed form.
  - Phreatic surface position (entrance corrected) vs the basic parabola.
- **Capture:**

  | Quantity | xslope | Analytical | Diff (%) |
  |---|---|---|---|
  | Discharge q | ____ | ____ | ____ |
  | Phreatic exit / surface point | ____ | ____ | ____ |

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
  | Total discharge q | ____ | ____ | ____ |
  | Phreatic surface (sampled x) | ____ | ____ | ____ |
  | Head at check node(s) | ____ | ____ | ____ |

---

## Mode 3 — FE Slope Stability (SSRM)

Elastic-perfectly-plastic Mohr-Coulomb, plane strain, viscoplastic algorithm.
Failure = non-convergence at the critical strength reduction factor (match G&L).

### SSRM-1 — Griffiths & Lane (1999) Example 1 (homogeneous)  [already matched]
- **Source:** Griffiths & Lane (1999), Geotechnique 49(3):387-403; Rocscience
  re-run (Sec. 3) with tabulated geometry:
  https://www.rocscience.com/assets/resources/learning/papers/Application-of-the-Finite-Element-Method-to-Slope-Stability.pdf
- **Geometry/materials:** homogeneous slope; E = 100000, nu = 0.3, gamma = 20,
  phi = 20 deg, c = 10 (consistent SI), dimensionless group c/(gamma*H) = 0.05.
- **Reference FOS:** Bishop 1.40 / Griffiths FE 1.42.
- **Status:** xslope already reproduces ~1.42 in the code; confirm and record.
- **Capture:** xslope FOS = ____ (target ~1.42), Diff (%) = ____.

### SSRM-2 — Griffiths & Lane (1999) Example 6 (earth dam, free surface)
- **Source:** G&L (1999); Rocscience App. Ex. V tabulated geometry.
- **Geometry/materials:** two-sided earth dam with free surface; c = 13.8 kPa,
  phi = 37 deg, gamma = 18.2 kN/m3.
- **Reference FOS:** full reservoir FE ~2.45; reservoir empty ~1.85.
- **What to run:** SSRM with the seepage free surface for the full-reservoir case,
  and the empty/no-pool case.
- **Capture:**

  | Case | xslope FOS | G&L FOS | Diff (%) |
  |---|---|---|---|
  | Example 6 (full reservoir) | ____ | 2.45 | ____ |
  | Example 6 (empty) | ____ | 1.85 | ____ |

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
