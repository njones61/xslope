# XSLOPE Graphical User Interface

A design + implementation plan for a cross-platform desktop application that wraps
the existing `xslope` Python package. The GUI lets a user open a problem from an
Excel file, view and edit the inputs graphically, run LEM / seepage / FEM analyses,
and view results — without writing code or running a notebook.

---

## 1. Goals & Constraints

**Goals**

- Cross-platform desktop app (macOS + Windows at minimum; Linux for free if we choose well).
- Open an existing filled-out Excel problem, render its geometry and inputs graphically.
- Edit every input through forms/dialogs; the graphics update automatically after each edit.
- Save / Save As back to the **same Excel format**; new problems start from the standard template.
- Run **LEM**, **Seepage**, and **FEM** analyses, plus a **meshing** step, each with a run-options dialog.
- Display results, including the two LEM display modes (search results vs. solution).
- Zoom/pan/select on the main graphics canvas; double-click an input to edit it.
- DXF import (to start a project) and DXF export (of the current display).
- Optional "smart editing" so coincident features move together.
- Display-options dialog: per-feature visibility toggles + color/line-style/size.

**Hard constraints**

- The computational engine **must remain the existing `xslope` Python package** — we wrap it, we do not re-implement it in another language.
- We keep the **existing Excel file format** for load and save.

---

## 2. Technology Stack (decided)

**PySide6 (Qt for Python) + embedded Matplotlib**, single Python process, calling `xslope` directly. Reasons tied to this codebase:

- **Reuse all plotting for free** — every low-level `plot_*` helper already takes an `ax`, so we embed a Matplotlib `FigureCanvasQTAgg` and call the existing plotting code unchanged. The single biggest lever.
- **No language/process boundary** — the GUI imports `xslope` and calls functions directly; `slope_data`/results (Shapely geometry, NumPy arrays, pandas DataFrames) pass by reference, so there is no serialization layer to build or maintain.
- **Zoom/pan + pick built in** — `NavigationToolbar2QT` gives zoom/pan/home; Matplotlib events (`button_press_event`, `pick_event`) give double-click-to-edit.
- **Rich forms/tables** — Qt's `QTableView`/`QFormLayout`/`QDockWidget` suit the materials/vertex tables and run-options dialogs.
- **Background work** — `QThread` workers run long solves (searches, SSRM) with progress/cancel.
- **Licensing & packaging** — PySide6 is LGPL (compatible with Apache-2.0; PyQt is GPL — avoid); PyInstaller/Briefcase produce native installers.

*Ruled out:* web/Electron/Tauri/PyWebView (forces rewriting the Matplotlib plotting layer in JS, or accepting image-based round-tripped interactivity, plus a process/serialization bridge — all cost, no benefit for a desktop-only goal); Tkinter (weak tables/forms); Streamlit/Dash (reactive model fights canvas editing). A hosted web version stays possible later by wrapping the same engine functions in FastAPI, without affecting the desktop app.

---

## 3. Application Name — **XSlope Studio** (decided)

The app is named **XSlope Studio** — it keeps the brand while "Studio" clearly marks the
desktop application as distinct from the `xslope` Python library (the engine import stays
`xslope`). The window title, about box, installer, and `console_script` use this name.

---

## 4. High-Level Architecture

Single Python process, document-view (MVC-ish) structure:

```
┌─────────────────────────────────────────────────────────────┐
│                     PySide6 Application                       │
│                                                               │
│  ┌─────────────┐   ┌──────────────────────┐  ┌────────────┐  │
│  │  Editors /  │   │   Canvas (Matplotlib  │  │  Results / │  │
│  │  Forms      │◄─►│   FigureCanvasQTAgg)  │  │  Log panes │  │
│  │  (dock)     │   │   zoom/pan/pick       │  │  (tabs)    │  │
│  └─────┬───────┘   └──────────┬───────────┘  └─────┬──────┘  │
│        │                      │                     │         │
│        ▼                      ▼                     ▼         │
│  ┌────────────────────────────────────────────────────────┐  │
│  │     ProjectDocument  (single source of truth)          │  │
│  │   slope_data dict  +  derived results  +  dirty flag    │  │
│  └───────────────────────────┬────────────────────────────┘  │
│                              │ direct function calls          │
│        ┌─────────────────────┴───────────────────┐           │
│        ▼               ▼               ▼          ▼           │
│   QThread workers for long-running engine calls              │
└────────┼───────────────┼───────────────┼──────────┼─────────┘
         ▼               ▼               ▼          ▼
   xslope.fileio   xslope.solve/   xslope.seep/  xslope.mesh
   (load/SAVE)     search/slice    fem           cad (DXF)
```

**Core idea:** the GUI is fundamentally a *structured editor over the `slope_data` dictionary*. `slope_data` is the single source of truth. Any edit mutates `slope_data` → emits a "changed" signal → the canvas re-renders the affected layers. Results (slice_df, fs_cache, seep/fem solutions) are derived artifacts cached on the document and shown in dedicated result views.

### Key components

- **`ProjectDocument`** — owns `slope_data`, the source file path, a dirty flag, and cached results. Provides typed accessors/mutators per input category and emits change signals. Owns **snapshot-based undo/redo** (deep-copies `slope_data` onto an undo stack on each edit).
- **`CanvasWidget`** — embeds a Matplotlib figure; renders the current view by calling existing `xslope.plot` helpers with the document's `slope_data` and a `StyleConfig`. Handles pan/zoom (toolbar) and pick events (double-click → open editor).
- **Editor panels/dialogs** — one per input category (see §6). Tables for tabular data, forms for scalars.
- **Run controllers** — assemble kwargs from run-options dialogs, launch engine calls on a worker thread, stream progress/log, store results on the document.
- **Result views** — tabbed canvases for the different display modes (see §7).
- **Style/Display manager** — holds a `StyleConfig` (visibility + style per layer) driving the renderer and the Display Options dialog.

---

## 5. Engine Integration — what exists vs. what we must add

The Explore of the codebase confirms a clean, callable API. Mapping of GUI actions → engine calls:

| GUI action | Engine call(s) |
| --- | --- |
| Open file | `load_slope_data(path)` → `slope_data` |
| Render inputs | `plot_inputs(...)` or low-level `plot_*` `ax` helpers |
| Generate slices (single surface) | `generate_slices(slope_data, circle=…, non_circ=…, num_slices=…)` |
| Run a method | `oms/bishop/janbu/spencer/corps/lowe/mprice(slice_df)` → `{'FS': …}` |
| Auto-search (circular) | `circular_search(slope_data, method, …)` → `(fs_cache, converged, search_path, circle_cache)` |
| Auto-search (non-circular) | `noncircular_search(slope_data, method, …)` |
| Reliability | `advanced.reliability(slope_data, method, circular=…)` |
| Rapid drawdown | `rapid=True` flag through the above |
| Plot LEM solution | `plot_solution(slope_data, slice_df, failure_surface, results)` |
| Plot search | `plot_circular_search_results(...)` / `plot_noncircular_search_results(...)` |
| Build mesh | `get_material_polygons(...)` + `build_mesh_from_polygons(...)`; `export_mesh_to_json` |
| Seepage | `build_seep_data(mesh, slope_data)` → `run_seepage_analysis(...)`; `plot_seep_*` |
| FEM | `build_fem_data(slope_data, mesh)` → `solve_fem` / `solve_ssrm`; `plot_fem_*` |
| DXF import | `cad.import_dxf(...)` / `read_dxf_polygons(...)` |
| DXF export | `save_dxf=True` on the active `plot_*` call, or `cad.export_dxf(...)` |

### Gaps to build **in the `xslope` package** (so notebooks/scripts benefit too)

1. **`save_slope_data_to_xlsx(slope_data, path, template=…)`** — ✅ **DONE** (`xslope/fileio.py`). Inverse of `load_slope_data`: writes the source input tables back into the template via the formatting-preserving `write_cells_to_xlsx`. `template=<blank>` creates a new file (New / Save As); `template=None` edits in place (Save). Verified by a load→save→load round-trip across 13 files spanning all input categories, wired into `run_tests.py` (`--roundtrip`).
2. **Embeddable, styleable plotting.** Add optional `ax=`/`fig=` and `style=` (`StyleConfig`) params to the high-level `plot_*` functions (`plot_inputs`, `plot_solution`, the search/seep/fem plots) so the GUI passes its embedded Matplotlib figure and applies per-layer styling, instead of each function creating its own figure with hardcoded styles. Keeps a single layout path — no GUI-side re-implementation. Start by reusing plots as-is (figure embedding first), then thread `style=` through incrementally.
3. **Progress/cancel callbacks (nice-to-have).** Searches and SSRM print to stdout today — capture that for the log pane initially; later thread an optional `progress_callback` through `circular_search` / `solve_ssrm` for a real progress bar and cancellation.

---

## 6. Input Editing

`slope_data` categories and their editors (each opens from the menu/side-panel tree, and via double-click on the corresponding canvas artist):

| Category | `slope_data` keys | Editor type |
| --- | --- | --- |
| Global / main | `gamma_water`, `tcrack_depth`, `tcrack_water`, `k_seismic`, `max_depth` | Form |
| Materials | `materials` (name, γ, option, c, φ, cp, r_elev, d, ψ, u, σ's, k1, k2, α, kr0, h0, E, ν) | Table (cols vary by analysis mode) |
| Profile lines | `profile_lines` (coords + mat_id) | Table per line + on-canvas vertex pick |
| Polygons | `polygons` (per-material zones) | Table / derived from profile |
| Piezometric lines | `piezo_line`, `piezo_line2` | Vertex table |
| Failure surfaces | `circular`, `circles` (Xo,Yo,Depth,R), `non_circ` (X,Y,Movement) | Form (circles) / vertex table (non-circ) |
| Distributed loads | `dloads`, `dloads2` | Table per load block |
| Reinforcement | `reinforce_lines` / `reinforcement_lines` | Table |
| Piles | `pile_lines` | Table |
| Seepage BC | `seepage_bc`, `seepage_bc2` (specified heads, exit face) | Table + on-canvas pick |

**Interaction model**

- A dockable **Inputs panel** (tree or accordion) lists categories; selecting one opens its editor.
- **Double-click on the canvas** maps the picked Matplotlib artist back to a `slope_data` object and opens the right editor (requires tagging artists with identifiers when rendering).
- Edits validate, mutate `slope_data`, mark the document dirty, and trigger a re-render of affected layers.
- **Smart editing (global toggle, later phase):** maintain a coincidence index so that moving a profile vertex also moves coincident dload/piezo/reinforcement/BC points. Implement as an opt-in pass over `slope_data` after a geometry edit.

---

## 7. Running Analyses & Results Display

### Run-options dialogs (kwargs map directly to engine signatures)

- **LEM:** method (`oms…mprice`), `num_slices`, analysis type (single surface / auto-search / reliability), surface type (circular / non-circular), `rapid` drawdown, tolerances (`fs_tol`, `tol`, `max_iter`).
- **Meshing:** `target_size` or size-divisions (auto), `element_type` (`tri3/tri6/quad4/quad8/quad9`), remesh toggle; writes `{stem}_mesh.json`.
- **Seepage:** convergence `tol`, BC set (1 / 2), variable to plot (head…), contour levels, flowlines/vectors.
- **FEM:** single vs. SSRM, `F` / `F_min` / `F_max`, `failure_criterion`, element type, mesh size, `tolerance`.

### Long-running execution

Run on a `QThread` worker; capture stdout into a **Log/Console pane**; show a progress bar and a **Cancel** button (cooperative cancellation once callbacks are threaded through). Results are stored on the document and handed to the appropriate result view.

### Result views (tabbed canvas area)

A **view selector / tab strip** above the canvas switches between modes, each its own embedded figure:

- **Inputs** — `plot_inputs` (mode = lem/seep/fem changes the material table columns).
- **LEM · Search** — `plot_circular_search_results` / `plot_noncircular_search_results` (all trial surfaces + critical + search path).
- **LEM · Solution** — `plot_solution` (slices, base stresses, thrust line) for the critical surface.
- **LEM · Reliability** — `plot_reliability_results`.
- **Seepage** — `plot_seep_data` (mesh+BC) and `plot_seep_solution` (head contours / phreatic / flowlines).
- **FEM** — `plot_fem_data` (mesh+BC+reinforcement), `plot_fem_results` (deformation / shear strain / displacement vectors), `plot_ssrm_convergence`.

This cleanly resolves the "two LEM display types" question: they're just two tabs (Search and Solution) populated from the same run.

---

## 8. Canvas: Zoom / Pan / Select / DXF

- **Zoom/pan:** the figure is rendered to a pixmap shown in a `QGraphicsView`, so the **whole figure** (axes, title, labels, margins) zooms/pans uniformly like an image viewer — wheel-zoom anchored at the cursor, drag-to-pan, Fit / 100% / +/− (not Matplotlib's axes-only data zoom). To keep text crisp instead of pixelating on zoom-in, the pixmap is **re-rasterized on demand**: the scene uses fixed logical units (inches × `BASE_DPI`) while the backing bitmap is redrawn at a DPI matched to the current zoom (× a supersample factor), debounced and capped to bound memory (see `studio/canvas.py`).
- **Select / edit:** enable `picker` on rendered artists, tag each with its `slope_data` reference; double-click → open editor (§6).
- **Display Options dialog:** a `StyleConfig` of per-layer `{visible, color, linestyle, linewidth/size, alpha, ...}` plus figure-level settings (background, font size, DPI); the renderer reads it. Editing re-renders live. (Drives the §5.2 styling work, and is persisted per §8a.)
- **DXF export:** "Export current view to DXF" calls the active `plot_*` with `save_dxf=True` (already supported across plot functions), or `cad.export_dxf(slope_data, path)` for a clean geometry export.
- **DXF import:** wizard → `cad.read_dxf_polygons` to list layers → map layers to materials → `cad.import_dxf(...)` writes into a template copy → `load_slope_data` opens it.

### 8a. Style persistence (per-project + global default)

Styling (colors, line styles, widths, visibility, fills, figure background/fonts) is
**not** stored in the Excel file — the input template format stays untouched. Instead it
is kept in JSON sidecars, resolved in three deep-merged tiers:

1. **Factory defaults** — shipped with the app (`studio/resources/style_factory.json`). Always present, so "Reset to factory" can never fail.
2. **User global default** — written by a **"Set as default"** action. Lives in the OS app-config dir via Qt's `QStandardPaths.AppConfigLocation` (macOS `~/Library/Application Support/XSlope Studio/style_default.json`, Windows `%APPDATA%\XSlope Studio\…`) — not a home-dir dotfile.
3. **Per-project sidecar** — `{stem}_style.json`, written **next to the `.xlsx`**, matching the existing `{stem}_mesh.json` / `{stem}_seep.csv` sidecar convention. So styling travels with the project when the xlsx + sidecar are shared, without altering the xlsx.

**Effective style = factory ⊕ global ⊕ project** (deep-merge). Each file stores only
**sparse overrides relative to factory** (not a full snapshot) plus a `version` field —
so files stay tiny, remain portable (a colleague sees your exact styling regardless of
*their* global default), and unset keys automatically pick up new factory values when the
app updates. The `StyleConfig` the renderer reads is the merged result; the plot path
never needs to know which tier a value came from.

Example `{stem}_style.json`:

```json
{
  "version": 1,
  "layers": {
    "failure_surface": { "color": "#cc0000", "linewidth": 2.0 },
    "slices":          { "facecolor": "#e8f0ff", "alpha": 0.3 },
    "piezo_line":      { "color": "#1f77b4", "linestyle": "--", "visible": true }
  },
  "figure": { "dpi": 300, "font_size": 9 }
}
```

**Lifecycle hooks:** Open → load `{stem}_style.json` if present, else global, else factory.
Display Options edits → live `StyleConfig` re-render. Save / Save As → also write
`{stem}_style.json` beside the xlsx (only if it differs from default, so unstyled projects
don't litter sidecars). "Set as default" → write current overrides to the global file.
"Reset to factory" / "Revert to default" → clear the relevant override tier.

---

## 9. File / Project Lifecycle

- **New** → copy standard template (`docs/inputs/input_template.xlsx`) into a working doc; user fills via forms; `save_slope_data_to_xlsx` writes it.
- **Open** → `load_slope_data`; render Inputs view. Auto-load `{stem}_mesh.json` / `{stem}_seep.csv` if present (already handled by `load_slope_data`), and `{stem}_style.json` if present (§8a).
- **Save / Save As** → `save_slope_data_to_xlsx` (§5.1), preserving the existing format; also writes the `{stem}_style.json` style sidecar when the style differs from default (§8a).
- **Recent files**, dirty-state prompt on close, and per-project sidecars (`{stem}_mesh.json`, `{stem}_seep.csv`, `{stem}_style.json`) all following the existing `{stem}_*` convention.

---

## 10. Suggested Project Structure

The GUI lives in a **`studio/`** subfolder inside the existing repo, separate from the engine package:

```
xslope/                # engine package (unchanged; now includes save_slope_data_to_xlsx)
studio/                # XSlope Studio desktop app
  __init__.py
  app.py              # QApplication entry point  → console_script "xslope-studio"
  document.py         # ProjectDocument (slope_data + results + signals)
  canvas.py           # Matplotlib canvas widget, pick/zoom/pan
  renderer.py         # StyleConfig + compose plot_* helpers
  views/              # inputs, lem_search, lem_solution, seep, fem result views
  editors/            # one module per input category (forms/tables)
  runners/            # lem, seep, fem, mesh run-controllers (QThread workers)
  dialogs/            # run-options, display-options, dxf import/export wizards
  style.py            # StyleConfig + 3-tier resolve/merge + sidecar I/O (§8a)
  resources/          # icons, bundled blank template, style_factory.json
```

- The app imports the engine as `import xslope` — no copy of the engine, single source of truth.
- Ship **both** distribution channels (decided): `pip install "xslope[gui]"` (pulls PySide6) for Python users, **and** packaged native installers (`.dmg` / `.msi`) for non-Python users.
- Entry point `xslope-studio` (console_script) for dev launch.

---

## 11. Phased Roadmap

**Phase 0 — Engine prerequisite** ✅ **DONE**
- `save_slope_data_to_xlsx(slope_data, path, template=…)` implemented in `xslope/fileio.py` and verified by a round-trip test across 13 representative files (all input categories). See §5.1. The round-trip check is wired into `run_tests.py` as a `roundtrip` test type (runs in the default suite; standalone via `python run_tests.py --roundtrip`).

**Phase 1 — Skeleton + read-only viewer** ✅ **DONE**
- `studio/` PySide6 app: shell, menus, dockable Inputs summary panel + Log pane, embedded Matplotlib canvas with zoom/pan + scroll-zoom, LEM/Seep/FEM mode selector, recent files.
- File → Open renders the Inputs view via `plot_inputs(fig=…)` (the `fig=` param added to `plot.py`). `ProjectDocument` holds `slope_data` with snapshot undo/redo. Packaged as the `xslope-studio` entry point (`gui` extra). Verified headlessly (offscreen).

**Phase 2 — Editing** ✅ **DONE**
- ✅ Modal editor framework (`studio/editors.py`): `Field` spec + reusable `_EditableTable`, `TableEditorDialog`, `FormEditorDialog`, `TabbedTableEditorDialog`; launched by **single-click** on an Inputs-tree category; edits go through `ProjectDocument.begin_edit/commit_edit` → re-render + dirty. Editors preserve unshown record keys.
- ✅ Editors for: global params, materials, circles (Option Depth/Radius/Intercept), non-circular surface, piezometric lines (2 tabs), distributed loads (2 sets), Head BC (2 sets: specified heads + exit face), piles (optional-float fields), reinforcement (rebuilds the derived display format via the extracted `build_reinforce_lines`), and profile-line geometry (master/detail, resyncs polygons/ground surface).
- ✅ Save (in place) / Save As (fresh file from bundled blank template); unsaved-changes prompt on close and before New/Open.
- ✅ Polygon-sheet geometry editor (for polygon-based files; profile-based files use the profile editor). Shares one `MatGeometryDialog` master/detail with the profile editor; the Inputs tree marks **Polygons** editable only when there are no profile lines.
- ✅ New (from template) — `ProjectDocument.new()` seeds a minimal valid in-memory skeleton (the `xslope_simple1` single-material slope + one circle, since a blank template won't `load`); edited via the forms, written out by Save As.

**Phase 3 — LEM analysis**
- Run-options dialog; worker-thread execution; Search + Solution + Reliability result tabs.
- Single-surface, auto-search, rapid drawdown, reliability.

**Phase 4 — Meshing + Seepage + FEM**
- Meshing dialog; Seepage run + result view; FEM single/SSRM + result views; progress/cancel.

**Phase 5 — Canvas selection + Display Options + style persistence**
- Pick/double-click-to-edit; `StyleConfig` + Display Options dialog (visibility/colors/styles).
- Style persistence (§8a): factory/global/project 3-tier resolve+merge, `{stem}_style.json` sidecar I/O, "Set as default" / "Reset to factory" actions.

**Phase 6 — DXF + Smart editing + polish**
- DXF import wizard and export-current-view; optional coincident smart-editing; undo/redo.

**Phase 7 — Packaging & distribution**
- PyInstaller or Briefcase native installers (`.dmg`/`.msi`); macOS code-signing/notarization; bundle gmsh for FEM; CI build matrix.

---

## 12. Decisions (all settled)

- **GUI stack:** PySide6 + embedded Matplotlib, single-process, calling `xslope` directly (§2).
- **Application name:** XSlope Studio.
- **Code location:** a `studio/` subfolder inside the existing repo (§10).
- **v1 scope:** LEM first — open / edit / save + LEM analysis and results; Seepage and FEM follow (Phase 4).
- **Distribution:** both — native `.dmg`/`.msi` installers *and* `pip install xslope[gui]`.
- **Undo/redo:** snapshot-based — deep-copy `slope_data` onto an undo stack per edit (the dicts are tiny, so this gives full undo/redo cheaply) (§4).
- **Canvas editing (v1):** forms/tables + double-click-to-edit only; the canvas is view + selection. Interactive drag of geometry is deferred to a later phase (§6, §8).
- **Plot embedding:** add optional `ax=`/`fig=` (and `style=`) params to the high-level `plot_*` functions so the GUI passes its embedded figure and applies styling — one layout path, no GUI-side re-implementation (§5).
- **Document model:** single project per window (no MDI) in v1.
- **Recompute:** explicit **Run** to (re)solve; only inputs re-render live on edit.
- **Units:** the engine is unit-agnostic; the `gamma_water` field implicitly sets the system (≈62.4 → English pcf/psf/ft, ≈9.81 → metric kN/m³/kPa/m). The GUI may *infer* that from `gamma_water` to show convenient unit labels, but performs **no enforcement or conversion** — the user is responsible for confirming units are consistent.
- **Packaging detail:** `studio/` ships inside the `xslope` distribution (added to `pyproject` packages + the `xslope-studio` entry point), enabled via the `gui` extra.

Nothing blocking remains; detail-level choices (icon/branding, menu layout) settle during implementation.

---

## 13. Key Risks / Watch-items

- **Save fidelity** — done and covered by the `--roundtrip` test; keep that test green as input categories evolve.
- **Long solves blocking the UI** — must run on worker threads from the start; FEM/SSRM also need cancellation.
- **Restyling existing plots** — the hardcoded styles in `plot_*` mean Display Options needs either added style kwargs or a GUI-side renderer; scope carefully.
- **Packaging heavy deps** — gmsh (FEM) and the scientific stack inflate installer size and complicate signing.
- **Coincident smart-editing** — genuinely tricky geometry bookkeeping; keep it opt-in and late.
```