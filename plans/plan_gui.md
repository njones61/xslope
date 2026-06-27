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

- **Zoom/pan:** the figure is rendered to a pixmap shown in a `QGraphicsView`, so the **whole figure** (axes, title, labels, margins) zooms/pans uniformly like an image viewer — wheel-zoom anchored at the cursor, drag-to-pan, a checkable **Zoom-to-box** tool (drag a rectangle to zoom into it, via `RubberBandDrag`), Fit / 100% / +/− (not Matplotlib's axes-only data zoom). To keep text crisp instead of pixelating on zoom-in, the pixmap is **re-rasterized on demand**: the scene uses fixed logical units (inches × `BASE_DPI`) while the backing bitmap is redrawn at a DPI matched to the current zoom (× a supersample factor), debounced and capped to bound memory (see `studio/canvas.py`).
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
- **Open** → `load_slope_data`; render Inputs view. Auto-load `{stem}_mesh.json` / `{stem}_seep.csv` if present (already handled by `load_slope_data`), and `{stem}_style.json` if present (§8a). **Auto-restore saved solutions:** if the mesh is present, rebuild seep/FEM data on it and read back any `{stem}_seep.csv` / `{stem}_seep2.csv` (both BC sets) and `{stem}_fem_nodes/elements.csv` sidecars (engine `import_seep_solution` / `import_fem_solution`, the inverses of the export functions), populating the Seep · Solution (per BC) and FEM · Results tabs immediately — no re-solve. The FEM SSRM factor of safety is restored from a `{stem}_fem_meta.json` sidecar (`import_fem_meta`), since it is not in the node/element CSVs. Best-effort: a sidecar whose node/element count no longer matches the mesh is skipped, not fatal.
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
- ✅ New — `ProjectDocument.new()` creates an **empty** project: every `slope_data` category present but empty, no geometry, blank canvas. The Inputs tree starts at zero and stays editable (Profile lines editable so the first line can be added); the user builds up via the editors (add a material, then a profile line — the profile editor rebuilds the ground surface), and Save As writes it through the template.

**Phase 3 — LEM analysis** ✅ **DONE**
- ✅ End-to-end solve framework. `RunLemDialog` (method / analysis / surface / num_slices / rapid / diagnostic) → `LemRunner` (`QThread`, `studio/runners.py`) runs the engine off the GUI thread; engine output streams to the Log pane live via the thread-safe stdout tee. Central area is a `QTabWidget` view strip (Inputs + lazily-added result tabs, cleared on Open/New); a run logs a banner + result.
- ✅ **Single surface** — circular *and* non-circular (`generate_slices` + `solve_selected`) → **LEM · Solution** tab via `plot_solution(fig=…)`.
- ✅ **Auto-search** — circular (`circular_search`) and non-circular (`noncircular_search`) → **LEM · Search** tab (`plot_circular_search_results` / `plot_noncircular_search_results`, both given `fig=`) showing all trial surfaces + critical + search path, plus the critical surface in the Solution tab. Search iteration progress streams to the log.
- ✅ Rapid drawdown flag wired through single/search; `fig=` added to all LEM plot functions, mirroring `plot_inputs`.
- ✅ **Reliability** (`advanced.reliability`) → **LEM · Reliability** tab (`plot_reliability_results`, given `fig=`) + Solution tab for the MLV surface. A determinate **progress bar** in the status bar tracks the `1 + 2N` searches via a `progress_callback` threaded through `reliability()` (engine-side, so notebooks benefit too); other runs show a busy bar. Surface-type choice is hidden unless the file has both circular and non-circular.
- ✅ **Cooperative cancel** — a `cancel_check` callable is threaded through `circular_search` / `noncircular_search` / `reliability` (engine-side; checked at iteration boundaries, raises `AnalysisCancelled`). The worker exposes `cancel()` (sets a `threading.Event`) and a Cancel button by the progress bar; cancelling aborts cleanly without killing the thread, emits `cancelled`, and leaves no result tab. Close also cancels an in-flight run.
- ✅ Remaining solver tolerances in the dialog (`fs_tol`, `tol`, `max_iter`). A **"Search tolerances"** group on `RunLemDialog`, enabled only for auto-search / reliability (the search-driven analyses); `tol` (geometric refinement) is greyed out for non-circular surfaces, which have no such parameter. Threaded through `LemRunner` to `circular_search` / `noncircular_search`, and through `reliability()` (engine-side — added `fs_tol`/`tol`/`max_iter` params that forward to its `1 + 2N` internal searches, so notebooks benefit too). Unset values fall through to each engine function's own default.

**Phase 4 — Meshing + Seepage + FEM** 🚧 **IN PROGRESS**
- ✅ **Build Mesh** — `BuildMeshDialog` (element type; target size, entered or auto-sized as slope-width / divisions per the drivers) → `MeshRunner` (`QThread`) builds via `get_material_polygons` + `build_mesh_from_polygons`, including reinforcement/pile constraint lines so the mesh also serves FEM. Result shows in a **Mesh** tab (`plot_mesh`, given `fig=`) and the mesh is stored on `slope_data['mesh']` (so it appears in the Inputs view) and written to the `{stem}_mesh.json` sidecar. Build progress shows a busy bar.
- ✅ **Mode-driven Run** — one Run action whose label/dispatch follow the LEM/Seep/FEM mode; Build Mesh shows only in Seep/FEM; Seep/FEM Run gated on a built mesh (`_update_run_actions`).
- ✅ **Run Seep** — `RunSeepDialog` (BC set 1/2, tol, plot variable head/u/v_mag/i_mag, contour levels, flowlines/vectors/fill/phreatic) → `SeepRunner` (`QThread`) runs `build_seep_data` → `run_seepage_analysis`; **Seep · Data** + **Seep · Solution** result tabs (`plot_seep_data` / `plot_seep_solution`, both given `fig=`); solution written to `{stem}_seep.csv` (`_seep2.csv` for BC 2). Convergence trace streams to the Log pane. **Each BC set keeps its own tab pair** (BC 2 → "Seep · Data 2" / "Seep · Solution 2"), so rapid-drawdown problems show BC 1 and BC 2 side by side and switch freely; results are held per BC in `doc.results["seep_solutions"][bc]`.
- ✅ **Display dock** — per-plot-type display options live in a context-sensitive "Display" dock under the Inputs tree (not in the run dialogs); its page follows the active result tab and re-renders the cached result live. Panels for every view that has options, each exposing the full set of display kwargs the underlying `plot_*` accepts: Inputs (material table + placement, legend column layout), LEM · Search (highlight critical), LEM · Solution (slice numbers, seep contours), Mesh (nodes / labels / padding), Seep · Data and FEM · Data (BC / nodes / labels / fill opacity; FEM also BC symbol size), Seep · Solution (variable, levels, fill opacity, vector scale, padding, flow lines / vectors / fill / phreatic), FEM · Results (shear strain / deformation / displacement vectors — the diagnostic stress/strain/yield/displace_mag types are omitted; deform %, mesh / reinforcement / element labels, and the displacement-vector options gated to that plot type). Styling the `plot_*` functions hardcode (colors, widths, fonts) is still deferred to the Phase 5 `StyleConfig`. Reliability has no plot options, so it shows the placeholder.
- ✅ **Image export** — a per-view **"Save…"** button on each canvas toolbar writes the current figure to PNG / PDF / SVG via `savefig`; PNG prompts for a DPI, vector formats skip it. Saves at the figure's true inch size so resolution is independent of on-screen zoom. Per-view (not per-panel) so Inputs / Reliability can export too. The Save button also offers **DXF export** of the current view.
- ✅ **Auto legend layout** — engine-side measure-based legend layout across all result-view plots (wide but neat, multi-column), with **legend column controls** exposed on every result-view Display panel (so notebooks benefit too).
- ✅ **Run FEM** — `RunFemDialog` (single trial / SSRM, F or F_min/F_max + tolerance + failure criterion) → `FemRunner` (`QThread`) runs `build_fem_data` → `solve_fem` / `solve_ssrm`; **FEM · Data** + **FEM · Results** tabs (`plot_fem_data` / `plot_fem_results`, both given `fig=`), display options (plot type, deform %) on the Display dock. SSRM reports FS and supports **cancel** (a `cancel_check` threaded through `solve_ssrm` + its helpers); solution exported via `export_fem_solution`.
- Decisions: Run is **mode-driven** with dynamic text ("Run LEM/Seep/FEM"); a mesh is **explicit** (Build Mesh first — Seep/FEM stay disabled until `slope_data['mesh']` is present), not auto-built.

**Phase 4 — done** (Meshing + Seepage + FEM all wired; display panels for all views; image export incl. DXF; auto legend layout). Remaining polish: progress granularity for FEM.

**File lifecycle (§9) updates landed:**
- ✅ **Sidecar reconciliation on save** — Save/Save As reconciles the `{stem}_mesh.json` / `{stem}_seep(.2).csv` / FEM sidecars against the current document so stale solution files don't outlive the inputs they came from.
- ✅ **Stale-result invalidation** — editing inputs invalidates a stale LEM solution, and a geometry edit that makes the mesh invalid clears the mesh (Seep/FEM re-gate on a rebuild).

**Phase 5 — Canvas selection + Display Options + style persistence** 🚧 **PARTIAL**
- ✅ **Pick / double-click-to-edit** on the Inputs canvas. Since the canvas is a rasterised pixmap (no Matplotlib pick events), the double-click is mapped viewport → scene → axes data coords (`studio/canvas.py`, `picked` signal) and hit-tested against the input geometry (`studio/picking.py`) to resolve the editor category — profile/polygon lines, piezo lines, distributed loads, reinforcement, piles, circles, non-circular points, with a material-zone interior falling back to the materials editor — then opens the same editor the Inputs tree uses. Inputs view only; result canvases stay view-only. Tolerance tracks zoom (~12 px in data units).
- ⬜ `StyleConfig` + Display Options dialog (per-layer visibility / colors / line-styles). *(Per-view Display dock panels exist from Phase 4, but the per-layer `StyleConfig` over colors/line-styles is not yet built.)*
- ⬜ Style persistence (§8a): factory/global/project 3-tier resolve+merge, `{stem}_style.json` sidecar I/O, "Set as default" / "Reset to factory" actions.
- ✅ **Canvas rendering polish (§8):** Fit frames the content bbox (not the whole figure) with a cushion; crisp text on Retina by reading the device-pixel ratio from the screen and matching render DPI to the fitted scale; autofit retries until the shown tab is laid out and re-fits on each (re)render.

**Phase 6 — DXF + Smart editing + polish** 🚧 **PARTIAL**
- ✅ **DXF export** — offered on each view's Save button (export current view).
- ✅ **DXF import** — File → Import DXF brings polygons into the live document.
- ⬜ DXF import *wizard* (layer→material mapping); optional coincident smart-editing; full undo/redo coverage.

**Phase 7 — Packaging & distribution** 🚧 **PARTIAL**
- ✅ **Custom app icon** — branded "X" app icon for Dock / taskbar.
- ⬜ PyInstaller or Briefcase native installers (`.dmg`/`.msi`); macOS code-signing/notarization; bundle gmsh for FEM; CI build matrix.

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

---

## 14. AI Chat Assistant (scoping — not yet built)

A dockable **chat panel** inside Studio that drives the app and the engine with
natural language and images, against a user-chosen model (Claude, OpenAI/GPT,
local Ollama, …). The aim is **far more than the existing `/xslope` skill**:
not just "build an input file from a sketch," but *any* interaction with the
project data and the `xslope` Python API — e.g. "vary the slope angle 20→30° and
plot FS vs angle," "add a piezo line 2 ft below the crest and re-run Spencer,"
"why did the SSRM not converge?".

### 14.1 Use cases (what "more than the skill" means)

- **Build / edit inputs** — from a sketch (vision) or prose: materials, geometry,
  piezo lines, loads, reinforcement, BCs.
- **Run analyses** — LEM / seepage / FEM / meshing, single or parametric sweeps.
- **Arbitrary scripting** — parameter studies, custom plots, batch over files,
  post-processing of `slice_df` / `fs_cache` / seep & FEM solutions, exporting
  CSVs — i.e. anything you'd do in a notebook with `import xslope`.
- **Explain / debug** — interpret results, diagnose non-convergence, suggest fixes.

### 14.2 This is an agent with code execution, not a chatbot

The decisive capability is a **Python execution tool** the model drives — a
persistent kernel/REPL with `xslope` imported and the **live document in scope**
(`slope_data`, `results`, the current file path, a Matplotlib figure sink). Most
requests reduce to "write and run a little Python, show me the figure/number,"
exactly like the CLI/notebook workflow the user already relies on. Everything
else (typed convenience tools, the skill knowledge) sits *around* that core.

Two ways the agent can touch the project, likely both:
1. **Live document tools** — call into `ProjectDocument` (`begin_edit`/`commit_edit`)
   and the existing run-controllers so edits re-render on the canvas and land on
   the undo stack; results flow into the existing result tabs. Tight, native feel.
2. **Code execution** — a sandboxed-ish `run_python(code)` with the engine and the
   document handed in, capturing stdout + any figures. This is what unlocks
   "anything a notebook can do" and is the main reuse of the skill's existing code
   patterns. Figures it produces are shown inline in the chat (and/or a result tab).

**Build = populate the in-memory project, not write a file (decided, Studio path
only).** The standalone `/xslope` skill ends by *saving an `.xlsx`* and **remains a
first-class feature unchanged** — the CLI/IDE, file-based workflow (surgical xlsx
writer and all) keeps working for people who want it. Inside Studio, the chat takes
a *different* build path: the agent **populates the live `slope_data` dict of a
fresh (or current) project** — which renders immediately on the canvas for review —
and the user persists it with **Save As**. So *for the Studio path* the surgical
xlsx writer isn't needed (the document is built directly); what Studio reuses from
the skill is its **schema knowledge** (the sheet/category layout, geometry rules,
examples) mapped onto `slope_data` keys and the §6 editor `apply` functions. File
write/reload stays available as a secondary path. Net: two coexisting front ends to
the same engine — the file-first standalone skill *and* the document-first Studio
assistant — not a replacement of one by the other.

### 14.3 Provider abstraction (bring-your-own-model)

Each provider differs in API shape and tool-call format (Anthropic Messages,
OpenAI function-calling/Responses, Ollama's OpenAI-compatible endpoint, etc.).
Options, roughly increasing in effort/control:

- **LiteLLM** (recommended starting point) — one Python interface over Anthropic /
  OpenAI / Ollama / many more, including tool-calling and vision; we implement the
  agent loop once. Optional `xslope[ai]` extra.
- **Provider SDK adapters** — thin in-house `LLMProvider` interface with per-vendor
  adapters; more code, no extra dep, full control of quirks.
- **Subprocess an existing agent** (e.g. the Claude Code CLI / Claude Agent SDK)
  pointed at the project dir — reuses a battle-tested agent + the skill verbatim,
  but is Claude-only and shells out, so it doesn't satisfy the multi-provider goal.
  Possible as a fast Claude-only spike to validate UX before building the native loop.

A **Settings** panel selects provider + model, stores API keys (OS keychain via
`keyring`, *not* plaintext in QSettings), and sets the Ollama base URL for local
models. Capability varies by model — tool use and vision aren't universal (esp.
small local models); the UI should degrade (e.g. disable image attach when the
selected model has no vision).

### 14.4 Tool surface (what the model can call)

- `run_python(code)` — **the core**; persistent namespace, engine + live document
  preloaded, returns stdout/result + captured figures.
- `get_inputs()` / `update_inputs(category, data)` — structured read/edit of
  `slope_data` (reuse the §6 editor `apply` logic so validation/resync is shared).
- `build_mesh` / `run_lem` / `run_seep` / `run_fem` — reuse the Phase 3/4 runners
  (already threaded + cancellable); results populate the existing tabs.
- `get_results()` — FS, flowrate, SSRM FS, convergence, etc.
- `show_plot(view)` / figure capture — surface plots in chat and/or the tab strip.
- `read_file` / `write_file` (incl. the xlsx) — so the file-based skill patterns
  still work and Studio reloads the document after.
- Image attachment (vision) for the sketch→inputs use case.

### 14.5 Reusing the existing skill

The standalone `/xslope` skill stays as-is (see §14.2) — this section is about how
*Studio* draws on it. `docs/usage/claude/xslope.md` (the skill body — template cell
layout, the surgical xlsx writer, per-sheet helpers, run snippets) becomes **domain
knowledge injected into the Studio system prompt**, not the limit of behavior, and
Studio leans on the schema parts rather than the file-writer (§14.2). Caveat:
it's written for a file-first agent and currently lives in `docs/` and is
repo-bound (the `SKILL.md` `!`cat …`` include only resolves inside the repo). To
use it from an installed Studio we must **package the prompt** (ship it under
`xslope/resources/` and locate it like the template) and keep the docs master and
packaged copy in sync — the same two-copy + `run_tests.py` guard pattern we just
set up for the template.

### 14.6 Architecture & UI

- **`ChatDock`** — a right-side `QDockWidget`: transcript (markdown + inline
  figures + collapsible "ran code" blocks), input box with image attach, model
  picker, Stop button.
- **`AssistantController`** — owns the conversation, runs the agent loop on a
  `QThread` (streaming tokens + tool calls back to the UI), dispatches tool calls
  to the document/engine/runners, and feeds results back to the model.
- **`providers/`** — the LiteLLM-or-adapter layer behind a small interface.
- **`tools/`** — the tool schemas + their Python implementations (the bridge to
  `ProjectDocument`, runners, and the kernel).
- **Kernel** — an in-process namespace (or a Jupyter kernel) for `run_python`;
  in-process is simplest but shares the GUI process (see safety).

### 14.7 Safety, trust & cost

- **Arbitrary code execution** is the point and the risk — same trust model as
  Claude Code: the agent can read/write files and run Python as the user. Need an
  **approval mode** (auto-run vs confirm-before-run, at least for `write_file`/
  `run_python`), a visible log of every action, and easy Stop/undo (the snapshot
  undo stack already helps for input edits).
- **In-process vs sandbox** — running in the GUI process is simplest but a bad
  snippet can hang/crash Studio; a subprocess/Jupyter kernel isolates it at the
  cost of plumbing. Start in-process behind confirm-to-run; revisit.
- **Network egress & keys** — prompts/images leave the machine for hosted models
  (Ollama stays local); make that explicit. Keys in the OS keychain. Surface
  token usage / rough cost.

### 14.8 Packaging

- Optional extra: `pip install "xslope[ai]"` pulls the provider lib (e.g.
  LiteLLM). Studio degrades gracefully (chat dock hidden/disabled) when it's absent.
- Bundle the skill prompt with the package (§14.5).
- Native installers would include the AI deps; keep them out of the base install.

### 14.9 Phased approach

- **A — Spike** ✅ **BUILT** — Claude-only (`anthropic` SDK, `claude-opus-4-8`,
  adaptive thinking) in `studio/ai/` (`kernel.py` persistent in-process namespace
  with `xslope` + the live `doc`/`slope_data`/`results`; `assistant.py` manual
  tool-use loop on a `QThread` with the single `run_python` tool, marshalling each
  call back to the GUI thread for a confirm dialog + document mutation + re-render)
  and `studio/chat_dock.py` (right-side dock: transcript with inline figures and
  "ran code" blocks, autonomy toggle confirm/auto, Stop). Skill prompt is the
  system context (`docs/usage/claude/xslope.md`, repo-bound for now). Optional
  `xslope[ai]` extra; the dock loads without `anthropic` and reports if absent.
  Verified end-to-end with a mocked client (text → tool_use → run_python on the
  live doc → result → final text). **Remaining for Phase A polish:** package the
  skill prompt (§14.5), streaming responses, API-key UX.
- **B — Multi-provider** ✅ **BUILT** — the agent loop now runs over **LiteLLM**
  (OpenAI message format), so the same loop drives Claude / OpenAI / local **Ollama**
  (free, no key). `studio/ai/config.py` holds the provider/model selection
  (`QSettings`) and API keys (OS keychain via `keyring`, QSettings fallback) and
  produces the `litellm.completion` kwargs; `studio/ai/settings_dialog.py` is the
  Settings dialog (provider, model, key, Ollama URL) opened from a **Settings…**
  button in the dock; the dock shows the active provider·model. `litellm.drop_params`
  smooths over per-provider param gaps. Verified end-to-end with a mocked
  `litellm.completion` (tool_calls round-trip) plus config/keyring round-trip.
  **Capability-aware UI** — `config.capabilities()` reports per-provider tool/vision
  support (None = depends on a local model); the dock shows a caption when the
  selected model can't / may not run code, and the Settings note states tool-use
  support per provider. **Claude prompt caching** is re-threaded: for Anthropic the
  system prompt goes out as a `cache_control: ephemeral` content block (plain string
  for other providers, which would 400 on list content), so the large skill prompt
  bills at cache-read rates on repeat turns. Refinement still open: adaptive thinking
  is Claude-specific and not yet re-threaded; per-model (not per-provider) capability
  detection.
- **C — Native tools + live document:** 🚧 **PARTIAL** — the live document, runners,
  and convenience helpers are reachable from the kernel: a preloaded **`run_lem()`**
  helper drives the LEM runner, a **`sensitivity()`** helper does parametric sweeps
  (writes a CSV + one plot to an accessible output folder), and assistant edits are
  **transactional** (resynced derived geometry, no auto-solve), invalidating stale
  mesh/solution like the GUI editors. Authoritative `slope_data` schemas are injected
  so the assistant builds the **live doc** (not an `.xlsx`); an empty project auto-opens
  so the first snippet works. **Still open:** dedicated structured input-edit / run /
  results *tools* (vs. `run_python` helpers), inline-figure refinements, richer
  confirm-to-run diffs.
- **D — Vision & polish:** 🚧 **PARTIAL** — **vision** is in (paste/drop images into the
  chat); provider coverage broadened beyond Claude / OpenAI / Ollama to **DeepSeek** and
  **Z.ai (GLM)** with per-model vision capability; chat UX polish (Enter-to-send,
  wrapping message blocks, "New chat" resets the kernel, actionable error messages).
  **Still open:** sketch→inputs ergonomics, cost meter, conversation save/restore.

### 14.10 Decisions (settled)

1. **Provider strategy** — **LiteLLM**: one dependency unifying Claude / OpenAI /
   Ollama (tool-calling + vision); the agent loop is written once behind a small
   interface. Optional `xslope[ai]` extra.
2. **Execution model** — **in-process kernel** for `run_python` (engine + live
   document directly in scope); a bad snippet can hang Studio, so it's gated behind
   the autonomy mode + a Stop button. Revisit isolation if it bites.
3. **Autonomy** — **user-selectable mode** (like Claude Code permissions): an
   *auto-run* mode and a *confirm-before-acting* mode (preview/diff + approve for
   code execution and file/document writes), chosen in Settings and switchable per
   session. Default to confirm.
4. **MVP target** — **Claude-only spike first** (Phase A): in-process `run_python`
   + file/document tools + skill prompt as system context + plain transcript;
   generalize to multi-provider (Phase B) once the UX is proven.
5. **Scope of "edit"** — the agent **populates / mutates the live `slope_data`
   document** (rendered live, persisted via Save As); file write + reload is a
   secondary path, not the primary one (see §14.2).