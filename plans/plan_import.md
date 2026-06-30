# Importing Models from Other Slope-Stability Software — Exploration & Plan

A high-value direction beyond DXF: **importing problems from other slope-stability
software** — GeoStudio **SLOPE/W**, Rocscience **Slide2**, and similar — mapping their
geometry, materials, piezometric/water conditions, and (where sensible) failure
surfaces onto `slope_data`, the same target the DXF importer and the Studio assistant
already populate. This would let users bring existing models into XSlope without
re-drawing them.

**Status:** 🟡 **Exploration — not scoped.** An initial scout (2026-06-27) confirmed
SLOPE/W `.gsz` is parseable and the strong first target. Next step is a feasibility
spike (§4). Extracted from `plan_studio.md` so it can grow on its own.

---

## 1. Goal

Let users open a model authored in another LEM/seepage package and get a populated
`slope_data` — geometry (material zones), materials, water/pore-pressure conditions,
and analysis definitions — plus a list of caveats for anything that doesn't map 1:1.
Surfaced in Studio as **File → Import → <format>**, mirroring the DXF import path
(engine-side parser → `slope_data` + caveats → live document → user reviews → Save As).

## 2. Target formats

- **GeoStudio SLOPE/W (`.gsz`)** — the first target (zipped XML; see findings below).
  Also the de-facto interchange format (Slide2 itself imports `.gsz`).
- **Rocscience Slide2 (`.slmd`)** — harder; undocumented/proprietary, samples not freely
  downloadable.
- Others (PLAXIS, Geostudio SEEP/W for seepage, etc.) — later, demand-driven.

## 3. Unknowns to resolve first

- **File formats & docs.** What does each package's native file look like, and is the
  format documented or reverse-engineerable? Some are XML-ish or text; others are
  opaque/binary or project bundles. SLOPE/W (`.gsz` — a zipped XML project) and Slide
  are the obvious first targets to investigate.
- **Sample files.** We need a corpus of real `.gsz` / Slide files (ideally with known
  answers) to develop and regression-test an importer.
- **Semantic mapping.** Their material models, pore-pressure definitions (Ru, piezo
  lines, pressure grids), and surface conventions don't map 1:1 onto `slope_data`;
  scope what's faithfully importable vs. flagged-for-review.
- **Per-importer module.** Likely one `xslope` importer per format (engine-side, so
  notebooks benefit too), surfaced in Studio as File → Import → <format>, each
  returning a populated `slope_data` + a list of caveats — mirroring the DXF path.

## 4. Next step — feasibility spike

Collect format documentation and a few sample files per package, then prototype a
parser for the most tractable one (likely SLOPE/W's zipped XML) to see how cleanly its
geometry/materials map across. Deliver a `.gsz → slope_data` proof of concept
(Regions→polygons, Materials, WaterItems) + a caveats list, on the PyGeoStudio samples.

## 5. Findings from an initial scout (2026-06-27)

- **SLOPE/W `.gsz` is the strong first target — confirmed parseable.** A `.gsz` is a
  ZIP holding the model as plain XML (`GSIData` root) + a `mesh_*.ply` + result CSVs.
  The XML maps closely onto `slope_data`: `<Geometries>` → **Points** `(X,Y)` / **Lines**
  / **Regions** (material zones → `polygons` + `mat_id`); `<Materials>` → `<Material>`
  (strength/hydraulic) → `materials`; `<WaterItems>` → pore-pressure/piezo;
  `<StabilityItems>` → slip surfaces; `<Analyses>` → analysis defs. (Verified by
  unzipping `Rapid drawdown.gsz` and walking the XML.)
- **Reference implementation:** [PyGeoStudio](https://github.com/MoiseRousseau/PyGeoStudio)
  — a Python `.gsz` reader/writer to study (or depend on). ⚠️ **No LICENSE file** as of
  the scout — check rights before reusing its code or redistributing its samples.
- **Sample files:** PyGeoStudio's `examples/GeoStudio_files/` has 5 real `.gsz` files,
  incl. `Rapid drawdown.gsz` and `Reinforcement with Anchors.gsz` (both relevant here).
  Official GeoStudio examples (Seequent/GeoSlope) exist but sit behind a Seequent-ID
  login. ⚠️ These are Seequent's example files — keep them in a **git-ignored** dev
  folder, don't commit, until licensing is clear.
- **Slide2 `.slmd`:** harder — samples ship with the install
  (`C:\Users\Public\Documents\Rocscience\Slide2 Examples`), not freely downloadable;
  format undocumented/proprietary. Notably Slide2 *imports* SLOPE/W `.gsz`, reinforcing
  `.gsz` as the de-facto interchange format and first target.
- **Recommendation:** start with `.gsz` → `slope_data` (Regions→polygons, Materials,
  WaterItems), using PyGeoStudio's samples + parser as the reference for a spike.

## 6. Connections to other work

- **van Genuchten support** ([`plan_vg.md`](plan_vg.md)) — native vG makes a `.gsz`/Slide
  importer **lossless** for hydraulic functions (their unsaturated curves are vG/Fredlund),
  instead of needing a fitted conversion to the linear-front model.
- **DXF importer** (`plan_studio.md` Phase 6) — the template: an engine-side parser →
  `slope_data` + caveats → live document → user fills placeholders → Save As. A `.gsz`
  importer should mirror that structure (and the Studio File → Import wizard).
- **Soils database** ([`plan_soils_db.md`](plan_soils_db.md)) — imported materials with
  missing properties could be backfilled from the soils-DB presets.
