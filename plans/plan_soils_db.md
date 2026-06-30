# Soils Database — Design & Strategy

A planned feature that lets a user pick **reasonable material properties** for a
problem from a curated, cited reference database — covering the LEM (strength),
seepage (hydraulic), and FEM (elastic) properties a single `xslope` material
carries, plus the reliability uncertainty columns. The big commercial packages
advertise a "soils database"; xslope's material model is unusually well-suited to
one because a single material row spans all three analysis types *and* carries
variability, and because there are recognized published databases for exactly
those pieces.

**Status:** 🟡 **Design / strategy** — not started. This document captures the
strategy and, critically, the **open decisions to resolve before any
implementation** (§7). Nothing is built yet.

---

## 1. Motivation & goals

**Goal.** Give users a fast, trustworthy starting point for material properties —
"this is a stiff clay (CL), give me reasonable γ, c, φ (and k, E, ν, and typical
COVs)" — that they then review and adjust. A *starting estimate*, never a
substitute for site-specific investigation.

**Why it fits xslope specifically.** An xslope material is not just c/φ. One row
holds (see `_blank_material`):

- **Strength (LEM):** `gamma`, `option` (`mc` / `cp` / `none`), `c`, `phi`, `cp`,
  `r_elev`, `psi` (Kc=1 envelope for rapid drawdown), `u` (pore-pressure option).
- **Reliability:** `sigma_gamma`, `sigma_c`, `sigma_phi`, `sigma_cp`, `sigma_d`,
  `sigma_psi` (standard deviations driving reliability analysis).
- **Seepage:** `k1`, `k2` (major/minor hydraulic conductivity), `alpha` (tensor
  orientation, usually 0 so k1=kx, k2=ky), `kr0`, `h0` (simplified unsaturated
  flow parameters).
- **FEM:** `E`, `nu`.

Each of those clusters has a recognized, citable reference source for *typical
values and ranges* (§5). The reliability cluster is the sleeper feature: most
commercial soils databases don't help with variability, but xslope already does
reliability analysis, so prefilling the σ-columns from published coefficients of
variation is genuinely differentiating.

**Non-goals (for now).** A user/company "save your own materials" library
(separate feature — §7 D11); a full unsaturated soil-water-characteristic-curve
fitting tool; automatic strength selection from in-situ test data (SPT/CPT
correlations) — possible later, out of v1.

---

## 2. What it is (three thin layers + one integration)

1. **Reference dataset** — a curated, versioned data file (likely JSON shipped in
   `xslope/resources/`) keyed by soil classification, where each entry holds a
   *representative value + range* (and where available a COV) for each property,
   each with a **citation**. This is the bulk of the real work — data curation,
   not code.
2. **Engine helpers** — an `xslope.soils` module: `soil_presets()` returns the
   data; `material_from_preset(class, drained=…, units=…)` returns a ready
   material dict. Engine-side so **notebooks/CLI benefit too** (the standard
   xslope pattern — every capability is added to the package first).
3. **Studio picker** — a **"Pick from soils database…"** button in the Materials
   editor that opens a browsable picker (filter by class; see value **and** range
   **and** source), filling the selected material row on choose. The user always
   edits afterward — it is a starting point, never a commitment.
4. **Assistant integration** — expose the dataset to the [Studio AI
   assistant](plan_studio.md) so "make layer 2 a soft clay" pulls consistent,
   cited values rather than improvising (cheap to add once the dataset exists).

---

## 3. How a material maps to the database

| xslope cluster | Fields | Database provides | Notes / mapping risk |
| --- | --- | --- | --- |
| Unit weight | `gamma` | typical + range by class | low risk; well-tabulated |
| Strength `mc` | `c`, `phi` | drained c′/φ′ by class | key by **class × drainage** (D4) |
| Strength `cp` | `c`, `cp`, `r_elev` | undrained su and su-increase rate | needs su / su-rate data, not c′/φ′; reference-elevation is site-specific (leave to user) |
| Rapid drawdown | `psi` | Kc=1 envelope friction angle | niche; likely defer |
| Reliability | `sigma_*` | COVs (→ σ = COV·value) | **Phoon–Kulhawy** COVs; differentiating |
| Seepage | `k1`, `k2`, `alpha` | saturated k range by class | order-of-magnitude k is textbook; default alpha=0, k1=k2 unless user sets anisotropy |
| Unsaturated | `kr0`, `h0` | (approx) from texture | xslope uses a **simplified** model, not van Genuchten — Carsel–Parrish vG params need conversion/approximation (D5) |
| FEM elastic | `E`, `nu` | typical E and ν by class | E spans a wide range (stress-level dependent) — present as range |

The two mapping subtleties to resolve are **drainage** (mc vs. cp / drained vs.
undrained — D4) and the **simplified unsaturated model** (D5).

---

## 4. Classification & granularity

Primary key is almost certainly the **Unified Soil Classification System (USCS)**
class (GW, GP, GM, GC, SW, SP, SM, SC, ML, CL, OL, MH, CH, OH, Pt) — the
classification practitioners already think in and the key most strength/k
correlations are tabulated against. The unsaturated cluster is the exception:
van Genuchten databases (Carsel–Parrish) are keyed by **USDA texture class**
(sand, loam, silty clay, …), so the unsaturated mapping needs a USCS↔texture
cross-walk (or its own texture key). Whether to also support AASHTO is an open
question (D1).

---

## 5. Data sources & scientific basis

| Cluster | Candidate source(s) | License posture |
| --- | --- | --- |
| γ, φ, c by class | NAVFAC DM-7.1/7.2; Das; Bowles; Duncan & Wright correlations | **NAVFAC = US-gov public domain** (preferred); textbooks = cite, re-derive, don't copy tables |
| Saturated k ranges | Das; Terzaghi–Peck–Mesri; Freeze & Cherry | order-of-magnitude k is common knowledge; cite |
| Unsaturated (vG) | **Carsel & Parrish (1988)** | journal — values are facts; cite, don't reproduce table formatting verbatim |
| E, ν ranges | NAVFAC DM-7; Bowles | as above |
| COVs (σ-columns) | **Phoon & Kulhawy (1999)** | journal — cite |

**Licensing principle (D6):** factual property *values* aren't copyrightable, but a
specific *table's selection/arrangement* can be. Lead with public-domain sources
(NAVFAC, USDA), and for journal-derived values, re-derive into our own schema with
clear citation rather than copying a paper's table. Get an explicit sign-off on the
sourcing before publishing the dataset.

---

## 6. Draft data schema (illustrative — not final)

Canonical units stored in SI; converted on apply (D3). Keyed by USCS class:

```json
{
  "version": 1,
  "units": "SI",
  "classes": {
    "CL": {
      "name": "Lean clay",
      "gamma":   { "typ": 19.0, "min": 17.0, "max": 21.0, "unit": "kN/m3", "src": "NAVFAC DM-7.1" },
      "drained": { "c": { "typ": 5,  "min": 0,  "max": 15, "unit": "kPa" },
                   "phi": { "typ": 28, "min": 25, "max": 32, "unit": "deg" } },
      "undrained": { "su": { "typ": 50, "min": 25, "max": 100, "unit": "kPa" } },
      "k":  { "typ": 1e-9, "min": 1e-11, "max": 1e-8, "unit": "m/s", "src": "Das" },
      "E":  { "typ": 15000, "min": 5000, "max": 40000, "unit": "kPa" },
      "nu": { "typ": 0.40, "min": 0.35, "max": 0.45 },
      "cov": { "phi": 0.12, "c": 0.30, "gamma": 0.05, "src": "Phoon & Kulhawy 1999" },
      "confidence": "established"
    }
  }
}
```

`material_from_preset("CL", drained=True, units="english")` → fills `gamma`, sets
`option="mc"`, `c`/`phi` from `drained`, `k1=k2=k.typ`, `E`, `nu`, and
`sigma_* = cov·value`; leaves `r_elev`, `psi`, geometry, etc. untouched.

---

## 7. Open questions / decisions to resolve (before implementing)

Each carries a recommendation to converge quickly.

- **D1 — Classification key.** USCS only, or also AASHTO / USDA texture?
  *Rec:* USCS primary, with a USDA-texture cross-walk used only for the unsaturated
  cluster. Defer AASHTO.
- **D2 — Value representation.** Single typical value, typical+range, or full
  distribution? *Rec:* typical + min/max range, plus an optional COV (the range is
  what engineers want to see; the COV feeds reliability).
- **D3 — Units.** Store both English & metric, or store SI canonical and convert?
  *Rec:* store SI canonical; convert on apply using the system implied by
  `gamma_water` (xslope is unit-agnostic, so this must be explicit).
- **D4 — Drainage / strength option.** How do entries map to `mc` vs. `cp`, and
  drained c′/φ′ vs. undrained su? *Rec:* key strength by **class × drainage**;
  drained→`option=mc` (c,φ), undrained→`option=cp` with su as `c` and a typical
  su-rate as `cp` (leave `r_elev` to the user).
- **D5 — Unsaturated model.** **Superseded by [`plan_vg.md`](plan_vg.md).** If xslope
  gains native van Genuchten support (planned), the Carsel–Parrish (1988) α, n table
  maps **1:1** onto the `vg_a`/`vg_n` columns with no conversion — so the unsaturated
  cluster becomes cleanly in-scope rather than deferred. *Rec (revised):* sequence the
  unsaturated cluster after vG lands; until then, omit it (no lossy kr0/h0 fitting).
- **D6 — Data sourcing & licensing.** Which references, and is reproducing their
  tables acceptable? *Rec:* lead with public-domain (NAVFAC, USDA); cite and
  re-derive journal values; get explicit sign-off before publishing.
- **D7 — Scope of v1.** Which clusters and how many classes? *Rec:* full USCS class
  list for **γ, c, φ (strength)** + **k ranges**, plus **Phoon–Kulhawy COVs**;
  add FEM E/ν and unsaturated in v2.
- **D8 — Surfacing.** Engine module + Studio picker, or one first? *Rec:* both,
  engine first (notebooks benefit), then the Studio picker.
- **D9 — Fill behavior.** Overwrite the whole row, or only known fields? Confirm
  before clobbering user-entered values? *Rec:* fill only the fields the preset
  knows, never silently overwrite a non-default user value (confirm), and record
  the provenance.
- **D10 — Confidence flags.** Mark which values are well-established vs. rough?
  *Rec:* yes — a per-entry `confidence` flag surfaced in the picker.
- **D11 — Personal/company library.** In scope, or separate? *Rec:* **separate,
  later feature** (save/reuse your own calibrated materials); out of v1.
- **D12 — Assistant integration.** Expose the dataset to the AI assistant? *Rec:*
  yes — cheap once the dataset exists; gives consistent cited values.
- **D13 — Disclaimer / framing.** How to present presets responsibly? *Rec:* always
  show the range + citation + "starting estimate" language; consider a one-time
  acknowledgment on first use.
- **D14 — Maintenance / versioning.** Where the dataset lives and how it evolves.
  *Rec:* versioned JSON in `xslope/resources/`, extensible, with a schema version
  field; keep the docs/packaged copies in sync like the template.

---

## 8. Phased roadmap (do not start until §7 is resolved)

- **Phase 0 — Decisions.** Resolve §7 (especially D4, D5, D6, D7). Confirm the
  `seep.py` unsaturated model and the strength-option semantics against the engine.
- **Phase 1 — Dataset + engine.** Curate the v1 dataset (D7 scope) into the JSON
  schema (§6); build `xslope.soils` (`soil_presets`, `material_from_preset`) with
  unit conversion and citations; round-trip test that a preset produces a valid
  material dict for each strength option.
- **Phase 2 — Studio picker.** A "Pick from soils database…" button in the
  Materials editor → a browsable picker showing value/range/source/confidence →
  fills the row (D9 behavior). Provenance shown; undoable like any edit.
- **Phase 3 — Breadth + assistant.** Add FEM E/ν, the unsaturated cluster (D5),
  and assistant access (D12). Optionally start the personal-library feature (D11).

---

## 9. Risks / watch-items

- **Professional responsibility.** Presets must read as starting estimates, not
  answers — framing (D13) is a feature, not a footnote. This is how the commercial
  tools handle it too.
- **Data quality & sourcing.** Garbage-in undermines trust; curation and citation
  (D6) are the make-or-break. Wide, honest ranges beat false precision.
- **Parameter-model mismatch.** The unsaturated (D5) and undrained-rate (`cp`,
  `r_elev`) mappings are where published correlations don't map 1:1 onto xslope's
  fields — scope carefully or defer.
- **Units.** Silent unit errors are dangerous in geotech; conversion (D3) must be
  explicit and tested in both systems.
