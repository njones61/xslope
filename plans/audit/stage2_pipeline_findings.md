# Stage 2 Geometry & Slice Pipeline Audit — Findings (June 2026)

Multi-agent audit (22 agents, 12 confirmed / 0 refuted). Read-only; no fixes applied.
Every confirmed defect is LATENT on the shipped benchmarks (left-facing, ascending-X,
top-to-bottom-ordered polygons, non-negative loads, in-range mat_ids) — no published
number changes. Source of truth for Stage-2 fix decisions.

# Stage 2 Audit Synthesis — Geometry & Slice Pipeline

## 1. Verdict on Overall Pipeline Validity

The geometry and slice pipeline is **fundamentally sound but carries three validity-critical distributed-load defects that must be fixed before publication**, plus a cluster of bounded silent-correctness faults. The core integration machinery is verified correct: slice geometry (alpha, y_cb, dl, weight), pore-pressure resolution (piezo vertical-distance rule and seep FE interpolation), circle-option three-way equivalence, and the distributed-load resultant convention all reproduce hand calculations to machine precision, and the pipeline passes unit-invariance (rel ~1e-16) and slice-convergence (smooth, monotone) cleanly. The failures are concentrated in two places: (a) the distributed-load path under right-facing geometry and descending-X input order, where OMS/Bishop silently return wrong-signed or zero load contributions, and (b) base-material assignment, which depends on polygon iteration order rather than depth. None of these corrupt any *published benchmark number* — every shipped benchmark is left-facing, ascending-X, top-to-bottom-ordered, and uses non-negative loads — so all confirmed fixes are no-ops on the verified path. The defects are latent landmines for user-authored inputs, not errors in the shipped results.

## 2. Confirmed Findings, Ranked by Severity

### VALIDITY-CRITICAL

**F10 — Distributed-load moment arm not sign-corrected for right-facing slopes (OMS + Bishop)**
- Area: facing-symmetry | Location: `solve.py:232` (oms `a_dx`), `solve.py:333` (bishop `a_dx`); `a_dy` at 220/334
- Defect: `a_dx = d_x - Xo` is built from raw coordinates with no right-facing sign correction, while the weight term uses signed alpha. On a right-facing slope the D driving moment flips sign, so an upslope/crest surcharge that must LOWER FS instead RAISES it. **Bishop has the identical defect** (the seed named OMS only; the battery's claim that "Bishop is symmetric" is FALSE).
- Evidence: Mirror of `xslope_simple_embankment` with crest surcharge: OMS/Bishop orig (left) = **0.98036** vs mirror (right) = **1.84701** (~88% error); Spencer/Janbu/Corps/L-K symmetric to ~2e-5. D-moment `(1/R)(sum_Dx - sum_Dy)` = **+6093.40** (left) vs **−6093.75** (right). With no load both facings give **1.28078**; the surcharge correctly drops left to 0.98036 but wrongly raises right to 1.84701.
- Fix: `right_facing = slice_df['y_lt'].iloc[0] > slice_df['y_rt'].iloc[-1]`; negate **only** `a_dx` (and `a_dy`) when right_facing; leave `a_s/kw` and `a_t/T` untouched. Negating `a_dy` too is wrong (inclined case breaks to 1.0844). Apply identically to oms and bishop.
- Published-benchmark impact: **None** — all six shipped benchmarks verified left-facing; fix is a no-op.

**DLOAD-ORDER-DROP — Entire distributed load silently dropped on descending-X input order**
- Area: distributed-loads / property battery | Location: `slice.py:859-861` (and 867 for dloads2)
- Defect: `np.interp` requires monotonically increasing `xp`; when dload polyline points are entered descending-X, the interp returns 0 everywhere and the **entire load is silently applied as zero** — no error, no warning.
- Evidence: `xslope_submerged`, critical circle, ns=40: ascending-X `sum(dload)=174362.3`, Bishop FS=**1.1539**; reverse point order only → `sum(dload)=0.0` exactly, Bishop FS=**0.6006**. Deterministic at ns=20/30/40/60.
- Fix: Sort each dload line's (X, Normal) pairs by ascending X before building the interp lambda (`order = np.argsort(xs)`), at `slice.py:859` and `:867`; or sort in `load_slope_data`.
- Published-benchmark impact: **None** — all shipped dload tabs are ascending-X; fix is a no-op.

> Note: F10 and the battery's OMS-BISHOP-DLOAD-MOMENT-MIRROR are the **same defect** (moment-arm sign), independently confirmed by two harnesses (embankment and foundation). One fix closes both.

### WRONG-BUT-BOUNDED

**SLICE-1 — Base material assigned by polygon iteration order, not depth**
- Area: slice-geometry | Location: `slice.py:933-935`; parser `fileio.py:178-233`
- Defect: `base_material_idx = mat_index` is set unconditionally for the LAST h>0 polygon in iteration order. The parser appends polygons in spreadsheet column order with no elevation sort and no validation, so out-of-order blocks silently bind the wrong base c/phi/u → wrong FS. Geometry and weight stay correct.
- Evidence: `xslope_acads_weak_layer` reordered [top,bottom,weak]: Spencer **1.99919 → 1.67895** (−16%), OMS 1.91628 → 1.61974; weight identical to the last digit (5010.763119335172, diff 0.0); 8/39 slices change base material. Shapely point-in-polygon ground truth: original 0 mismatches, reordered 8.
- Fix: Select base material by deepest h>0 layer — track `min(overlap_bot)`, guard `if overlap_bot < best_bot: best_bot=overlap_bot; base_material_idx=mat_index` (init `best_bot=+inf`). Verified 0 mismatches on reordered case, identical on original. Optionally also sort polygons top-to-bottom at load and/or warn.
- Published-benchmark impact: **None** — fix is identical to current code on correctly-ordered (shipped) inputs.

**PP-1 — Off-mesh seep base points silently get u=0**
- Area: pore-pressure | Location: `mesh.py:2621-2622`, consumed at `slice.py:1131`
- Defect: `find_element_containing_point` returns −1 outside the mesh and `interpolate_at_point` silently returns 0.0; `max(0.0,...)` preserves the zero. A base below the phreatic surface but outside the mesh gets u=0 → over-predicted FS, no warning.
- Evidence: `xslope_earth_dam_rapid`: 0/40 base points outside on the benchmark circle (latent, not triggered). Constructed deep circle: 23/43 points off-mesh, all u=0; Bishop FS=**6.255** (fallback) vs **3.248** (snapped to nearest mesh value) — ~1.9x unconservative. Mitigating safeguard: the circular search's depth clamp coincides with mesh bottom (182), so the production search path cannot reach sub-mesh points; trigger requires a hand-entered deep circle or a mesh not spanning domain depth.
- Fix: Distinguish −1 from a true zero; warn or clamp to nearest boundary value. Purely additive on the benchmark (fallback never fires).
- Published-benchmark impact: **None**.

**DLOAD-1 — Negative (uplift/suction) loads silently dropped by `qC > 0` gate**
- Area: distributed-loads | Location: `slice.py:942, 952`
- Defect: Strict-positivity gate drops any load whose center intensity is negative → D=0 for that slice.
- Evidence: Uniform Normal=−100 over width 10: current D=**0.0** vs intended **−1000.0**.
- Fix: Change gate to `if qC != 0` (and `qC2 != 0` at 952), or drop the heuristic and always compute qL/qR (interp left/right=0 makes it safe).
- Published-benchmark impact: **None** — all shipped Normal values ≥ 0.

**fileio-02 — Profile mat_id has no upper-bound check; silent wrong-material fallback**
- Area: fileio | Location: `fileio.py:307-317` vs `219-223`; fallback `slice.py:911-914`
- Defect: Polygon parser validates mat_id against `len(materials)` and raises; profile parser only checks `<0`. An out-of-range profile mat_id is stored verbatim, then `generate_slices` silently falls back to polygon-order index → wrong strength/unit-weight, no error. Docstring (line 250) claims a count check that does not exist.
- Evidence: `xslope_simple_mult_layers`, deep polygon mat_id typo'd 0→99: runs with no error, slice c-set silently becomes {400, 800} (deep slices bound to c=800 instead of intended c=400).
- Fix: In profile parser raise on `mat_id >= len(materials)` and on `<0`, mirroring `_parse_polygon_sheet`; preserve the None case.
- Published-benchmark impact: **None** — all shipped mat_ids in range.

**fileio-03 — Blank unit weight silently parsed as gamma=0 (weightless soil)**
- Area: fileio | Location: `fileio.py:359-361`, guard at `378`
- Defect: `_num()` coerces blank/non-numeric to 0.0; the `c_to_x_empty` guard requires ALL of C:X blank, so a row with c/phi set but gamma blank passes → gamma=0 → zero soil weight for that layer → silently higher FS (non-conservative direction).
- Evidence: Row {g:NaN, c:500, f:30} accepted, gamma=0.0; `soil_weight += h*gamma*dx = 0` for that layer (slice.py:932). The all-blank case is still caught by slice.py:567.
- Fix: For mc-option LEM materials require gamma present and >0; distinguish blank from explicit-0 via `pd.isna`. Condition on LEM use so seep-only rows (legitimately gamma-less) are unaffected.
- Published-benchmark impact: **None**.

**fileio-01 — FEM-only (SSRM) inputs wrongly rejected ("must include circular or non-circular surface")**
- Area: fileio | Location: `fileio.py:1067-1075`
- Defect: Validation exempts seep-only models but not FEM-only. A pure-FEM input (mesh + E/nu materials, no surface) is rejected, though `fem.py` reads no surface data. Latent (all shipped FEM files also carry a surface; 0/8 rejected).
- Fix: Add `is_fem_only` exemption symmetric to `is_seepage_only`.
- Published-benchmark impact: **None**.

**SEARCH-1 — Non-circular search collapses onto degenerate near-vertical surface for rigorous methods**
- Area: search | Location: `search.py:226`
- Defect: FS is sharply sensitive to initial `movement_distance`; a Free endpoint can slide to near-zero slice width, producing a near-vertical base, and greedy descent converges to a spurious low minimum — affecting **rigorous Spencer/Janbu**, not just Corps v1. Returns as a converged success, not a loud failure.
- Evidence: `xslope_noncircular`, Spencer, ns=20: md∈{2,4,5,6,7}=**1.7390**; md=8=**0.5571** (3.12x drop). Degenerate first segment 5 deep over 0.543 wide → alpha=−83.8°, n_eff_min=−7578, Z∈[−5432,6732]. Threshold non-monotone (md=7.5 safe, 7.9/8/8.5 collapse, 9.0 recovers).
- Fix: Characterization only (per task constraint). Real fix is the search-admissibility filter co-landing with the bracketed root finder (solve.py:668-678).
- Published-benchmark impact: **None** — all callers use default md=4.0.

**SEARCH-2 — Free force-equilibrium search descends toward degenerate (negative-base-normal/tension) surfaces**
- Area: search | Location: `search.py:24`
- Defect: Free non-circular search drives Corps-v1 onto inadmissible surfaces (tensile base normal) and reports a spurious low minimum.
- Evidence: `xslope_acads_weak_layer`: Corps-v1 **1.0480** vs Spencer **1.2791** (−18%). Corps-v1 min: 3/32 negative base normals (min −124.6), 23/32 tension boundaries. Cache: 181/198 (91%) have ≥1 negative base normal, 198/198 (100%) have tension.
- Fix: Characterization only; defer admissibility filter to land with bracketed root finder. **No code change here.**
- Published-benchmark impact: **None**.

### COSMETIC

**SLICE-2** (`slice.py:1161-1184`) — non-circular alpha via finite difference blends segments at a vertex; latent (benchmark max error 5.75e-12; constructed thin slice width 0.015 gives 41.19° vs true 45°). Optional analytic-segment fix.
**FS-1** (`slice.py:178`) — `yi < Yo` filter drops above-center intersections, yielding opaque "got 1" rejection; consistent with the bottom-half-only modeled arc, never a wrong FS (true crossings y=5.65/24.35, only lower returned). Doc/message clarity only.
**fileio-04** (`fileio.py:701-707`) — malformed surface coords default to NaN rather than failing at load; but NaN fails closed downstream (GEOSException / "got 1"/"got 0"), never a numeric FS. (Partial: auditor's "(False, got 1.)" for endpoint NaN was actually a GEOSException.) Add finite-value validation at load.

## 3. Property/Invariant Battery + F10 Check

- **Unit invariance: PASS.** Scaling all stress/weight quantities by k left FS identical to machine precision for all 7 methods on both piezo (worst rel 2.7e-16) and dry (4.3e-16). Holds even with dloads present in as-loaded ascending orientation.
- **Slice-count convergence: PASS.** ns=20/50/100/200 converged smoothly and monotonically toward a limit for every method, no jumps (dry maxΔ=3.5e-3 OMS/Bishop/Spencer, 1.1e-2 Janbu; piezo maxΔ<7.1e-3).
- **Mirror symmetry: PASS for dry, piezo/seepage, and seismic (kw=0.15, rel 1.6e-5); FAILS only for distributed loads.** The recently-added Corps/L-K interslice-angle negation on right-facing slopes (solve.py:768-771, 828-831) is confirmed mirror-symmetric. The mirror failures are exactly the two validity-critical dload bugs (F10 sign + DLOAD-ORDER-DROP), both input-order-/facing-specific. The residual ~1e-5 on passing cases is benign arc-discretization noise (XSLOPE structurally rejects perfectly symmetric "flat arc" surfaces), not a sign error.
- **F10 focused check: CONFIRMED, and it also hits Bishop.** Faithful mirror (no-dload baseline symmetric, proving the mirror is physical: OMS 1.2443 both sides). With crest surcharge, left correctly drops (Δ=−0.0510) but right wrongly rises (OMS Δ=+0.1368, Bishop Δ=+0.2429); Spencer symmetric (theta flip). Seismic and tension arms already symmetric, confirming only `a_dx/a_dy` are defective. In-memory fix (negate `a_dx/a_dy` when right_facing) restores right-facing OMS to 1.1933, exactly matching left, with left-facing unchanged.

## 4. Prioritized Fix Queue for Norm

| # | Finding | Severity | Effort | Risk to published numbers |
|---|---------|----------|--------|---------------------------|
| 1 | **F10 / OMS-BISHOP-DLOAD-MOMENT-MIRROR** — sign-correct `a_dx`/`a_dy` for right-facing in oms() and bishop() | validity-critical | low | None (all benchmarks left-facing) |
| 2 | **DLOAD-ORDER-DROP** — sort dload (X,Normal) ascending before interp at slice.py:859/867 | validity-critical | low | None (all benchmarks ascending-X) |
| 3 | **SLICE-1** — base material by deepest overlap, not iteration order (+ optional load-time sort/warn) | wrong-but-bounded | low | None (no-op on ordered inputs) |
| 4 | **fileio-02** — profile mat_id upper-bound check; **fileio-03** — require positive gamma for LEM materials | wrong-but-bounded | low | None |
| 5 | **PP-1** — warn/clamp off-mesh seep base points instead of silent u=0 | wrong-but-bounded | low | None (additive) |
| 6 | **DLOAD-1** — `qC != 0` gate (loops in with #1/#2 dload work) | wrong-but-bounded | trivial | None |
| 7 | **fileio-01** — add `is_fem_only` validation exemption | wrong-but-bounded | trivial | None |
| 8 | **SEARCH-1 / SEARCH-2** — DEFER to the admissibility-filter + bracketed-root-finder work unit; do NOT add a standalone filter (per design constraint) | wrong-but-bounded | large | None now (characterization only) |
| 9 | **SLICE-2, FS-1, fileio-04** — cosmetic robustness/message/validation polish; batch opportunistically | cosmetic | trivial | None |

**No fix in this queue changes a published benchmark number** — every confirmed defect is latent on the shipped inputs (left-facing, ascending-X, top-to-bottom-ordered polygons, non-negative loads, in-range mat_ids, surface-bearing FEM files). Recommended gate for publication: land #1 and #2 (validity-critical, both low-effort no-ops on benchmarks) plus #3 before release; #4–#7 are low-risk hardening that should ship alongside; #8 stays on the existing planned work unit.
---

## 5. Batch Resolution (applied)

Stage-2 fix batch #1–#7 applied and verified; benchmarks (run_lem.py) unchanged
across all three LEM cases — every fix confirmed a no-op on the shipped path.

| # | Finding | Resolution | Verification |
|---|---------|-----------|--------------|
| 1 | **F10** | `oms()`+`bishop()`: negate **`a_dx` only** when `right_facing` | Search-based mirror, see note below |
| 2 | **DLOAD-ORDER-DROP** | sort each dload line by ascending X before `np.interp` (slice.py) | descending input now = ascending (FS 1.1543) |
| 3 | **SLICE-1** | base material by deepest overlap (`base_overlap_bot`, init +inf) | reordered polygons = original FS 1.2786 |
| 4a | **fileio-02** | mat_id range check moved to validation block (materials parsed there) | benchmarks load clean |
| 4b | **fileio-03** | require positive gamma for materials referenced by geometry (LEM only; seep-only exempt) | benchmarks load clean |
| 5 | **PP-1** | `interpolate_at_point(..., return_found=True)`; warn on off-mesh seep base point instead of silent u=0 | rapid-drawdown seep LEM: u 0–5628, Bishop 2.23, 0 false warnings |
| 6 | **DLOAD-1** | `qC != 0` / `qC2 != 0` gate (was `> 0`) | negative loads now retained |
| 7 | **fileio-01** | `is_mesh_analysis` exemption (mesh + no LEM surface = seep/FEM run) | FEM/seep templates no longer rejected |

**F10 — resolved to `a_dx` ONLY, not `a_dy`.** §3's focused check suggested "negate
`a_dx/a_dy`", but it only exercised a *flat-crest* surcharge (β≈0, so the `a_dy`/sin β
term is ~0 and cannot distinguish the two fixes). A search-based mirror harness with a
**sloping-face** surcharge (β≠0, exercises both arms) settles it: negating `a_dx` only
restores symmetry to <0.01% on *both* the flat-crest and sloping-face cases, whereas
negating `a_dy` too would break the sloping-face case. Physical reason: `generate_slices`
already sign-flips β for right-facing (slice.py:896), and `sum_Dy = Σ D·sinβ·(Yo−d_y)`
inherits that correction; `sum_Dx = Σ D·cosβ·(d_x−Xo)` does not (cosβ is even in β), so
only its geometric arm needs flipping. `a_s`/`a_t` left untouched as specified.

Mirror results (acads_simple, crest & sloping-face surcharge), left vs right-facing:
flat-crest OMS 0.36743/0.36741, Bishop 0.48968/0.48966; sloping-face OMS 0.73533/0.73527,
Bishop 0.76841/0.76841 — all <0.01%. No-surcharge baseline symmetric to 0.001%.

**PILE-FACING — RESOLVED. Real battered-pile right-facing bug, now fixed in all methods.**
The F10 investigation flagged the pile moment term `H·sin(θ_p)·(x_pile − Xo)` as a latent
right-facing analog. Settling "real vs. discretization" required peeling apart **three**
layered effects (this is why early numbers looked contradictory):

1. **An apparent ~16% asymmetry was MY harness error** — I mirrored with θ_p→180−θ_p,
   treating θ_p as a global bearing. The convention is **θ_p relative to the resisting
   direction** (horizontal component always opposes the slide; θ_p only tilts up/down), so
   the same θ_p is kept on both facings.
2. **A ~1.5% asymmetry at θ_p=0 was ALSO a harness artifact** — my geometry mirror left
   `ground_surface` in *descending*-x order, and the Ito & Matsui auto-H path uses
   `np.interp(x, gs_x, gs_y)` (slice.py:1063, 1100), which silently returns garbage for
   non-ascending `xp`. That perturbed z_f → H → the M_cap/L_m cap. With `ground_surface`
   re-sorted ascending (as `load_slope_data` always produces — verified) the θ_p=0
   asymmetry vanishes (**0.0004%**). Not a production bug; production ground_surface is
   ascending. (Optional hardening: sort before those `np.interp` calls.)
3. **The remaining battered-pile (θ_p≠0) asymmetry is REAL** — ~2.9% (OMS/Bishop) and ~4.8%
   (Spencer) under the clean (sorted-gs) harness. Root cause: the vertical pile force's
   moment uses `x_pile − Xo` in real coordinates while α runs in the orientation-normalized
   frame (OMS/Bishop), and Spencer's "reflect everything" block negated the *entire* H_pile,
   wrongly flipping the facing-independent **vertical** (upward) component. Janbu/Corps/Lowe
   are force-only (no pile moment) and were already symmetric (0.001%).

**Fixes applied & verified (clean harness → all ≤0.0006%; no-op on shipped benchmarks):**
- `oms()`/`bishop()`: flip the `(x_pile − Xo)` moment arm when `right_facing` (same family
  as the F10 a_dx fix).
- `spencer()`: drop H_pile from the reflect-everything negation; instead flip **only** the
  horizontal component (`H_cos_tp`) for right_facing, leaving the vertical component alone.
- Janbu/Corps/Lowe: no change needed (verified symmetric).
- Published-benchmark impact: **None** (all shipped slopes left-facing; shipped xslope_piles
  FS unchanged to 6 digits). Convention documented in `docs/lem/piles.md` §"Setting the
  Force Angle".
- Stage-3 leftover (optional, low priority): defensive ascending-sort before the pile-path
  `np.interp` ground lookups; permanent pile mirror-symmetry property test in the harness.

### Still deferred (unchanged)
- **SEARCH-1 / SEARCH-2 (#8)** — stay with the admissibility-filter + bracketed-root-finder
  work unit (Stage-2 seed item in plan_comprehensive_audit.md). No standalone filter.
- **SLICE-2 / FS-1 / fileio-04 (#9)** — cosmetic; not addressed in this batch.
