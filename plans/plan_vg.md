# van Genuchten Unsaturated Model — Design & Scope

Add **van Genuchten (vG)** as an optional per-material unsaturated relative-permeability
model alongside xslope's existing **linear front (lf)** model. The material table
gains a model-selector column with two values — `lf` (linear front, the current
behavior and default) and `vg` (van Genuchten). This document captures the
rationale, the design decision, and a full change scope (input, solver, plotting,
docs, validation) before implementation.

**Status:** 🟡 **Design / scoped** — approved to implement; template redesign by the
author in progress. Open naming/parameter decisions in §9.

---

## 1. Rationale

The unsaturated kr curve is, for xslope's use case, close to a "don't care" on the
physics — yet there are real **perception** and **interoperability** reasons to
support vG, and a technical wrinkle that makes it far cheaper to add than its
reputation suggests.

**Physics — why it barely matters for results:**

- **Stability ignores suction.** Pore pressures above the phreatic surface are
  clamped to zero in LEM and in the FEM SSRM, so the kr shape in the unsaturated
  zone never reaches the stability result. Two kr models that place the phreatic
  surface in the same location give identical strength behavior.
- **Steady-state.** xslope solves ∇·(k∇h)=0 — there is **no storage term**. Most of
  vG's instability reputation comes from *transient* Richards solvers (the steep
  `C(h)=dθ/dh` capacity term, mass lumping, sharp wetting fronts) — none of which
  exists here. Steady-state carries only the kr(h) nonlinearity in the stiffness
  matrix, the milder half.
- **Flow rate.** Mostly set by the saturated zone and the phreatic-surface position;
  the unsaturated kr nudges it but rarely dominates for slope problems (exception:
  genuinely unsaturated-driven cases — infiltration into a dry embankment, thin
  saturated zones — not the core target).

**Why do it anyway:**

- **Perception.** vG is the de-facto "standard" unsaturated model; its absence reads
  as a gap even where it doesn't change answers. A vG option removes that, and lets
  us position the linear front as a *deliberate, robust choice* rather than a
  limitation.
- **Interoperability.** SLOPE/W-SEEP/W, Slide, and GeoStudio define hydraulic
  functions with vG (or Fredlund–Xing). Native vG makes the planned `.gsz` / Slide
  importer **lossless** for hydraulic functions instead of needing a fitted
  vG→(kr0,h0) conversion to defend.
- **Soils database.** Carsel & Parrish (1988) — the standard unsaturated dataset —
  *is a table of vG α and n by texture class.* Native vG makes that whole soils-DB
  cluster map **1:1 with no conversion** (this flips decision **D5** in
  [`plan_soils_db.md`](plan_soils_db.md) from "defer unsaturated" to "enabled
  cleanly").

**The wrinkle that de-risks it — steady-state vG kr needs only α and n.**
vG–Mualem relative permeability:

$$S_e = \big[\,1+(\alpha|h|)^n\,\big]^{-m},\quad m = 1-\tfrac{1}{n},\qquad
k_r = S_e^{\,L}\big[\,1-(1-S_e^{1/m})^{m}\,\big]^{2},\quad L=\tfrac12$$

θ_r and θ_s cancel out of kr — they only scale water content / storage, which a
steady-state solve never uses. So a vG kr option needs **exactly two parameters
(α, n)** — the same count as the current `(kr0, h0)` — and is a drop-in alternative
kr *function*, not a schema explosion. Because the unsaturated zone is non-critical,
the one remaining steep spot (dkr/dh near h=0⁻) can be tamed with a **kr floor** and
near-saturation smoothing **without any accuracy cost**, since it doesn't touch the
pore pressures that matter.

---

## 2. Design decision

- **New per-material column: unsaturated model selector**, values `lf` (linear
  front, **default** = today's behavior) and `vg` (van Genuchten). Short header —
  see §9 D1 (working choice: **`unsat`**).
- **vG parameters: α and n only** (steady-state ⇒ no θ_r/θ_s needed). New columns,
  named to avoid the clash with the existing `alpha` (permeability-tensor
  orientation) — see §9 D2 (working choice: **`vg_a`**, **`vg_n`**). The existing
  `kr0`/`h0` columns remain the **lf** parameters; the `unsat` column selects which
  pair the solver reads.
- **Regularized vG kr**: a configurable floor `kr_min` (≈1e-4) and air-entry
  smoothing near saturation, justified by the non-criticality of the unsaturated
  zone (§9 D3).
- **lf stays the default and the documented recommended model** for stability work;
  vG is offered for compatibility, import fidelity, and user preference.
- **Confined (saturated) solves are unaffected** — kr only enters the unsaturated /
  unconfined path.

---

## 3. Input / template changes

| Area | Change |
| --- | --- |
| Template | Add the `unsat` model column (`lf`/`vg`) and the `vg_a` / `vg_n` parameter columns to the **mat** sheet. Author is redesigning the template. |
| Back-compat | Missing `unsat` ⇒ default `lf`; existing files (no new columns) load unchanged and behave exactly as today. |
| `fileio.load_slope_data` | Read the new columns into each material dict (`unsat`, `vg_a`, `vg_n`), with defaults (`unsat="lf"`, `vg_a`/`vg_n` = NaN/0 when lf). |
| `_blank_material` (`document.py`) | Add the new keys so New/empty projects and DXF import carry them. |
| `save_slope_data_to_xlsx` | Write the new columns back (round-trip fidelity). |
| Validation | When `unsat=vg`, require valid `vg_a>0`, `vg_n>1`; when `unsat=lf`, require valid `kr0`/`h0` (existing check at `seep.py:168`). Clear error messages naming the offending material. |
| Template sync + roundtrip | Update `run_tests.py` roundtrip key list + the docs-master/packaged template sync (the two-copy guard). Add at least one `vg` material to a roundtrip fixture. |

---

## 4. Solver changes (`seep.py`)

The kr evaluation is centralized, so the change is localized:

- **New kr function** `kr_vg_vec(p, vg_a, vg_n, kr_min)` — the Mualem–vG kr above,
  vectorized over `(n_elements, n_sample_points)`, with the `kr_min` floor and
  near-saturation smoothing. (Pressure head `p<0` in the unsaturated zone; `p>=0`
  ⇒ kr=1.)
- **Dispatch** — replace direct `kr_frontal_vec` calls with a `kr_relative_vec`
  that takes per-element `model` + both parameter pairs and selects:
  `np.where(model==VG, kr_vg_vec(p, a, n), kr_frontal_vec(p, kr0, h0))` (both
  branches computed vectorized, then selected — the assembly already factors
  `ke = factor(kr_avg)·ke_sat`, so averaging/scaling is model-agnostic once kr is
  computed).
- **Thread the model + (vg_a, vg_n) arrays** through the same paths that carry
  `kr0/h0` today: `build_seep_data` (build the per-material arrays), the
  `mat_props` array (§ lines ~270, 346), `solve_unsaturated`, `_assembly_data`,
  and the post-processing `compute_velocity` / `solve_flow_function_unsaturated`.
- **Legacy per-element kr path** — `kr_frontal`, `_kr_factor`, and the
  `*_stiffness_matrix_kr` functions: confirm which are still live (velocity / flow
  function vs. the vectorized assembly) and give them the same dispatch, or route
  them through `kr_relative_*`.
- **Convergence** — keep the existing hybrid head/flow-closure test; add
  under-relaxation of the kr update if vG needs it (the floor usually suffices).
- **No change** to `solve_confined`, geometry, BC handling, or element matrices.

---

## 5. Plotting / output changes

- **Phreatic surface** is the `p=0` (or head=elevation) contour — independent of the
  kr model, so `plot_seep_solution` phreatic extraction is unaffected.
- **Velocity / flow vectors / flow function** use kr — they must read the per-material
  model (covered by the §4 threading).
- **Seep solution CSV** (`export/import_seep_solution`) stores nodal head/results,
  not kr params — unaffected.
- **`plot_seep_data`** may optionally annotate the unsaturated model per zone (nice,
  not required).

---

## 6. Documentation changes

- **`docs/seep/overview.md`** — add an "Unsaturated zone modeling" section:
  - Position the **linear front** model as the deliberate, robust default, and
    explain *why it's sufficient* (suction neglected in stability; steady-state
    pore pressures below the phreatic surface insensitive to kr shape). *(This is
    the "defuse the perceived limitation" track — worth writing regardless.)*
  - Document the **vG** option, the (α, n)-only steady-state formulation, the
    regularization, and when to choose which.
- **`docs/usage/input_template.md`** — document the `unsat`/`vg_a`/`vg_n` columns
  (the kr0/h0 section already exists at line ~142).
- **`docs/api/seep.md`** — picks up the new functions via mkdocstrings.
- **Studio docs** — the Materials editor / Run Seep pages mention the unsat model
  selector (minor).
- Per project convention: docs describe current behavior; no changelog prose.

---

## 7. Validation

- **Equivalence on benchmarks** — re-run the seepage sample suite (`docs/seep/samples.md`,
  `run_tests.py --seep`) with materials left as `lf` ⇒ **byte-for-byte unchanged**
  (the dispatch must not perturb the lf path).
- **lf vs vg comparison** — on a representative unconfined problem (e.g. the earth
  dam), solve with `lf` and with `vg` params chosen to approximate the same air-entry,
  and confirm: (a) the **phreatic surface** matches closely, (b) the **flow rate**
  matches within a stated tolerance, (c) the **stability FS** (LEM and SSRM) is
  unchanged — the evidence for the docs claim.
- **vG correctness** — verify `kr_vg_vec` against the closed-form Mualem–vG kr for a
  range of α, n (unit test), including the floor/smoothing behavior.
- **Robustness** — confirm vG converges on a dry-ish problem within comparable
  iterations to lf; record any need for under-relaxation.
- **Round-trip** — a `vg` material survives load→save→load (add to the roundtrip
  fixtures).
- **Import** — once the `.gsz` importer lands, a vG material imports without
  conversion (cross-link to that work).

---

## 8. Connections to other plans

- **Soils database** ([`plan_soils_db.md`](plan_soils_db.md)) — native vG flips **D5**:
  the Carsel–Parrish (1988) α, n table maps directly onto `vg_a`/`vg_n` with no
  conversion, so the unsaturated cluster can be in scope rather than deferred.
- **Model importers** (plan_studio.md §14) — vG support makes `.gsz` / Slide hydraulic
  functions import losslessly.

---

## 9. Open questions / decisions

- **D1 — Model-column header.** Short header for the `lf`/`vg` selector.
  Candidates: **`unsat`** (clearest to an engineer — "unsaturated model"),
  `kr_fn` / `kr_mod` (precise — "relative-permeability function"), `usm`.
  *Working choice:* `unsat`.
- **D2 — vG parameter column names.** Must avoid the existing `alpha` (tensor
  orientation). Candidates: **`vg_a` / `vg_n`**, `a_vg` / `n_vg`,
  `alpha_vg` / `n_vg`. *Working choice:* `vg_a`, `vg_n`. (Reuse vs. dedicated
  columns: dedicated is cleaner; `kr0`/`h0` stay as the lf params.)
- **D3 — Regularization.** Value of the `kr_min` floor and whether to add explicit
  air-entry smoothing (modified-vG) vs. a bare floor. *Rec:* start with a floor
  (~1e-4) + simple clamp; add smoothing only if convergence needs it.
- **D4 — Future-proofing for transient/θ.** Reserve θ_r/θ_s columns now (cheap) in
  case transient seepage is ever added, or omit until needed? *Rec:* omit — add when
  transient is actually scoped.
- **D5 — Other models.** Brooks–Corey or Fredlund–Xing too (SLOPE/W supports
  several)? *Rec:* vG only for v1; revisit per importer demand.
- **D6 — Default for new/blank materials.** `lf` (preserve current behavior) — yes.

---

## 10. Phased roadmap

- **Phase 0 — Decisions + template.** Resolve §9 (esp. D1, D2); author redesigns the
  mat sheet with the new columns; confirm the live kr code paths (vectorized
  assembly vs. legacy per-element).
- **Phase 1 — Engine.** `kr_vg_vec` + `kr_relative_vec` dispatch; thread model +
  (vg_a, vg_n) through `build_seep_data`, `solve_unsaturated`, `_assembly_data`,
  velocity/flow-function; fileio read/write + `_blank_material`; validation.
- **Phase 2 — Validation.** Unit tests for vG kr; lf-unchanged regression;
  lf-vs-vg equivalence on the earth dam (phreatic / flow / FS); roundtrip fixture.
- **Phase 3 — Docs.** Seepage overview unsaturated section (both tracks),
  input-template columns, Studio mentions.
- **Phase 4 — Downstream.** Wire vG into the soils-DB unsaturated cluster (D5 there)
  and the `.gsz` importer when those land.

---

## 11. Risks / watch-items

- **lf regression.** The dispatch must leave the lf path bit-identical — guard with
  the unchanged-benchmark test.
- **Convergence on dry problems.** vG kr is steep near saturation; the floor +
  optional under-relaxation must keep the free-surface iteration stable. Non-critical
  accuracy means we can regularize freely.
- **Parameter confusion.** Two α's (tensor orientation vs. vG) is a footgun — the
  distinct `vg_a` header and clear docs matter.
- **Units of α.** vG α has units of 1/length (inverse air-entry head) — must be
  consistent with the model's length unit (English vs. metric), like every other
  dimensional input.
