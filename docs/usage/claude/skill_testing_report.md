# XSLOPE "build-from-sketch" skill — testing & improvement report

## Objective
Stress-test the XSLOPE Claude skill (builds an XSLOPE model from a pasted sketch + supplemental
text, then runs LEM / Seep / FEM analyses), fix accuracy and robustness problems, and make sure
the improvements ship to the XSlope Studio AI assistant.

## Method
- **Blind-test harness:** general-purpose Sonnet subagents (Opus orchestration), each restricted
  to the skill doc + one sketch + supplemental text, sandboxed so they cannot read reference/KEY
  files or search the filesystem. Each returns a built `.xlsx` + JSON results + friction notes.
- **Diagnosis (per the user's method):** for every off result, `load_slope_data(reference)` vs
  `load_slope_data(agent)` and diff the dicts — compare data, don't theorize. This overturned two
  of my early guesses (point-contact junctions, BC placement) and found the real cause.

## Core architectural change (validated earlier in the effort)
The skill no longer writes Excel cells directly. It builds the in-memory `slope_data` dict and
calls `save_slope_data_to_xlsx(slope_data, path, template=...)` — the exact inverse of
`load_slope_data`, so round-trip fidelity is guaranteed. Result: **mechanism 11/11 clean** across
Round 4 — zero build/save failures; every agent with complete data produced a valid, loadable,
plottable model.

## Test outcomes
- **Clean wins (complete data → correct answer):** F03 (FEM SSRM 1.703 vs 1.68), L09 (~1–2%).
- **Guardrails working:** agents correctly STOPPED and asked instead of fabricating data when a
  sketch lacked essential inputs (conductivities, property tables, pile specs).
- **Dropped as invalid sketch tests** (not drawings / under-specified): L15, L10, L17, L18, F02,
  F04 (all `plot_inputs` renders); L05 (no horizontal scale); S07 (no conductivity given).

## The one real engine-usage bug found (high value) — FIXED
**Overlapping zoned-dam polygons → wrong seepage** (J01, S05). Agents built a dam's shell and
foundation as simple shapes that *overlapped* the core. Point-in-polygon material lookup still
resolved correctly, but the **mesher bridged the high-k shell over the low-k core**, short-
circuiting the barrier → ~5× too much flow → wrong pore pressures → wrong LEM/FEM.
- Proven by the diff method: agent geometry → flow 9.516; carve the core out (non-overlapping) →
  **1.958, exactly the reference**.
- **Fix:** `_validate_polygons_no_overlap()` in fileio.py now raises `ValueError` (naming both
  Mat IDs and the overlap area) if any two zones overlap; touching (shared edges) is allowed. Doc
  rewritten: zones must **tile** (no gaps, no overlaps, shared edges match); a cored section is
  one concave polygon with a notch, or use profile lines (`build_polygons` tiles correctly).
  Verified: catches the J01/S05 files, no false positives on valid polygon references, roundtrip
  15/15, regression green.

## Skill-doc improvements written this session (all synced to the packaged copy)
1. **No-overlap rule + hard validation** [J01, S05] — see above.
2. **Max Depth / base** (guideline #4) [L01] — infer the base from the drawn bottom of a material
   zone (any drawing style, not an enumerated list); Max Depth applies to **profile input only**;
   with polygons the base is implicit in each zone's bottom edge, so building a bounded zone as a
   single polygon removes the ambiguity.
3. **Don't infer undimensioned coordinates** (Phase 1 step 3) [L05] — if a whole direction lacks
   scale (e.g. only vertical thicknesses, no horizontal run), **stop and ask**; don't assume the
   drawing is to scale. Offers the "if it's to scale, shall I infer?" prompt.
4. **Trust labels over the drawing** (step 3) [L09] — sketches are often not to scale; when a
   coordinate must be measured, calibrate the pixel scale off a feature whose true size is labeled.
5. **Non-circular = search, not single-surface** [L07] — the drawn surface is only a search seed;
   Free endpoints on the ground surface, Horiz interior points seeded ~0.1 above the weak-layer
   base, run `auto_search`. Verified: search converges to 1.739 (vs the agent's single-surface
   1.99). Base elevation for the bottom material is guess-or-prompt (doesn't affect the surface).
6. **`staged=True` for reservoir/pore-pressure FEM** [J01] — apply gravity first (dry), then add
   water loads and pore pressures: construction history (built, then filled). Confirmed live and
   functional in `solve_fem` (two-stage viscoplastic loop, carries displacement + strains forward).

## Studio packaging — confirmed
- `studio/ai/assistant.py._load_skill_text()` loads the docs master in a checkout, falls back to
  the wheel-shipped `xslope/resources/xslope_skill.md` for pip installs.
- `run_tests.py` has a `template_sync` guard asserting the two copies are byte-identical — the
  build fails if the packaged copy drifts.
- The Studio trailer overrides only the file-save step (mutate `slope_data` in memory, never write
  an .xlsx); the shared dict-building core now carries all the fixes above. No stale
  nesting/overlap/write_cells language remains in the Studio path.

## Net
- Mechanism: robust (dict + `save_slope_data_to_xlsx`), 11/11 clean.
- One real engine-usage bug found and fixed (overlapping zones → wrong seepage), now guarded.
- Five workflow/interpretation gaps closed in the skill doc (base depth, scale inference, label
  trust, non-circular search, staged FEM).
- Remaining accuracy "misses" all trace to invalid/under-specified sketch test cases, not to the
  skill mechanism.
