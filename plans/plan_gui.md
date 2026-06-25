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

## 2. Technology Stack — Trade Study

Two families are on the table: a **native-Python desktop** stack (PySide6 + embedded
Matplotlib) and a **web front-end** stack (Electron/Tauri/PyWebView + the Python engine
behind an API). Both can satisfy "cross-platform desktop." They differ sharply in how
they handle this codebase's two defining facts: (a) all visualization is Matplotlib, and
(b) `slope_data`/results are full of **Shapely geometry, NumPy arrays, and pandas
DataFrames**.

### 2a. The web / Electron path, evaluated fairly

Concrete architectures if we go web:

- **Electron** (Chromium + Node) shell + Python engine as a **sidecar process**, talking over local HTTP/WebSocket (FastAPI) or JSON-RPC over stdio.
- **Tauri** (Rust shell, system webview) + Python sidecar — much smaller binary than Electron, same sidecar/IPC model.
- **PyWebView** — a lightweight system-webview window with Python in the *main* process and a JS↔Python bridge (no Node, no second runtime). The most Python-friendly web option.

Within any of these, visualization has only two honest options, and both have a real cost:

1. **Render Matplotlib server-side to PNG/SVG and ship the image to the web view.** Reuses *all* existing plotting code — good. But the canvas becomes an **image**: smooth zoom/pan/pick at the vector level is lost. You bolt on image pan/zoom, re-render on zoom via a **round-trip to Python** (laggy), and implement "double-click to edit" by mapping click pixels back through the exported axes transform. Workable, but the interactive feel is second-rate and the plumbing is non-trivial.
2. **Rewrite the plotting in JS (Plotly/D3/canvas).** Great client-side interactivity — but you re-implement *and then forever maintain a parallel rendering layer* for slices, base stresses, thrust lines, mesh, head contours, flowlines, FEM deformation/shear-strain/vectors, and the material tables (~370 KB of Matplotlib code across `plot.py`/`plot_seep.py`/`plot_fem.py`). Two plotting code paths that must stay in sync as the engine evolves. Highest cost in the whole project.

**The serialization tax (applies to Electron/Tauri sidecar, less so PyWebView):**
`slope_data` contains `LineString`/`Polygon` (Shapely), the mesh holds NumPy arrays, and
results hold pandas `DataFrame`s. None of these are JSON-serializable as-is. A process
boundary forces a **custom encode/decode layer for every category**, maintained in
lockstep with the engine. In-process Python passes them by reference — that layer simply
does not exist.

**Where the web path genuinely wins** (worth being honest about):

- If you later want the **same app to run in a browser or as a hosted/multi-user cloud service**, a web UI is already most of the way there. *(The stated goal is a local desktop app to replace Colab — not a web service.)*
- **Richer UI styling/theming** via HTML/CSS, and easier custom/modern look.
- If the team were **JS/web-skilled** rather than Python-scientific (not the case here).

**Net:** the web path's costs (image-based interactivity *or* a duplicated plotting
layer, plus a serialization boundary) are paid up front and forever, against benefits
(browser/cloud reach, CSS theming) that the current goal doesn't ask for.

### 2b. Future-proofing note

Choosing native **does not burn the web bridge.** The engine API is clean and
function-oriented; if a hosted/browser version is ever wanted, the *same* functions can
be wrapped in a FastAPI service later without touching the desktop app. So we can take
the low-cost native path now and keep the web option open as an additive future project,
rather than pay the web tax today.

### Recommendation: **PySide6 (Qt for Python) + embedded Matplotlib**

For concrete reasons tied to this codebase:

| Requirement | Why PySide6 + Matplotlib wins |
| --- | --- |
| Same language as engine | The GUI runs in the same Python process and calls `xslope` functions directly — **no IPC, no serialization boundary, no second language to maintain.** |
| Reuse existing plots | Every `plot_*` low-level helper already takes an `ax`. We embed a Matplotlib `FigureCanvasQTAgg` and call the **existing plotting code unchanged.** This is the single biggest lever — we get most of the visualization "for free." |
| Zoom / pan | Matplotlib's `NavigationToolbar2QT` gives zoom/pan/home out of the box; custom pick (double-click-to-edit) via Matplotlib event callbacks (`button_press_event`, `pick_event`). |
| Rich forms & tables | Qt has first-class `QTableView`, `QFormLayout`, `QDockWidget` (dockable side panels), validators, file dialogs, etc. — ideal for the materials table, vertex tables, and run-options dialogs. |
| Packaging | Mature: PyInstaller (one-folder) or Briefcase (native `.dmg` / `.msi` installers). |
| Licensing | PySide6 is **LGPL**, compatible with the project's Apache-2.0 license (PyQt is GPL/commercial — avoid). |
| Background work | `QThread` workers for long solves (searches, SSRM), with progress + cancel, while keeping the UI responsive. |

### Alternatives considered (and why not)

- **Electron / Tauri + web front-end + Python backend.** The user is open to a different wrapper language, but this path forces us to **rewrite the entire plotting layer** (Matplotlib doesn't render in a browser; we'd move to Plotly/D3) and adds a cross-process bridge (HTTP/WebSocket) plus a packaging headache (bundling a Python sidecar). All cost, little benefit here.
- **Tkinter + Matplotlib.** Built-in and embeds Matplotlib, but dated widgets and weak table/validation support make the forms-heavy parts painful.
- **Streamlit / Dash / Panel (web app in a webview).** Fast to prototype, but the reactive-rerun model fights interactive canvas editing (pan/pick/double-click), and it isn't really a "desktop app."
- **Toga / BeeWare, Dear PyGui, Kivy.** Either immature for this use or require rewriting all visualization.

**Decision:** PySide6 + embedded Matplotlib, single-process, engine called directly. (Open to confirming — see §12.)

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

- **`ProjectDocument`** — owns `slope_data`, the source file path, a dirty flag, and cached results. Provides typed accessors/mutators per input category and emits change signals. Owns undo/redo (command stack) eventually.
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

1. **Full `save_slope_data_to_xlsx(slope_data, path, template=…)`** — ✅ **DONE** (in `xslope/fileio.py`). The inverse of `load_slope_data`: writes the **source input tables** (main params, materials, profile/polygon geometry, piezo, circles, non-circ, dloads, reinforce, piles, seep BC) back into the template's sheet/cell layout, building on the package's existing formatting-preserving `write_cells_to_xlsx`. Contract: pass `template=<blank template>` to create a file ("New"/"Save As"); pass `template=None` to edit an existing file in place ("Save").

   Key correctness details handled: writes geometry to **either** the `profile` sheet (when `profile_lines` present) **or** the `polygon` sheet (mutually exclusive, matching the loader); skips **derived** geometry (`ground_surface`, `domain_polygon`, `polygons` when profile-sourced); never writes **formula** cells (the row-6 material-name XLOOKUP) but does write the Mat IDs that feed them; collapses every circle to `Option="Depth"` (loader recomputes `R = Yo - Depth` identically); leaves the pile `qp/theta` column blank (auto-derived on load); reconstructs reinforcement from the raw `reinforcement_lines` (the round-trippable form); leaves empty `option`/`u` cells blank.

   **Verified** with a round-trip test (load → save → load is identity) across 13 representative files spanning all categories — circular & non-circular surfaces, profile & polygon geometry, reinforcement, piles, distributed loads (both sets), second piezo line, reliability sigmas, and seepage BCs (both sets). All pass.
2. **Optional: a thin restyle hook on the renderer.** Existing `plot_*` functions hardcode colors/styles. To support the Display-Options dialog without forking the plotting code, add optional style kwargs (color/linestyle/linewidth/visibility) to the low-level `ax` helpers, *or* introduce a GUI-side renderer that composes the low-level helpers and applies a `StyleConfig`. Start by reusing plots as-is; add styling incrementally.
3. **Progress/cancel callbacks (nice-to-have).** Long searches and SSRM currently print to stdout. We can capture stdout for the log pane initially; later, thread an optional `progress_callback` through `circular_search` / `solve_ssrm` for a real progress bar and cancellation.

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

- **Zoom/pan:** Matplotlib `NavigationToolbar2QT` (built in). Add fit-to-extent ("Home") and scroll-wheel zoom.
- **Select / edit:** enable `picker` on rendered artists, tag each with its `slope_data` reference; double-click → open editor (§6).
- **Display Options dialog:** a `StyleConfig` of per-layer `{visible, color, linestyle, linewidth/size}`; the renderer reads it. Toggling re-renders. (Drives the §5.2 styling work.)
- **DXF export:** "Export current view to DXF" calls the active `plot_*` with `save_dxf=True` (already supported across plot functions), or `cad.export_dxf(slope_data, path)` for a clean geometry export.
- **DXF import:** wizard → `cad.read_dxf_polygons` to list layers → map layers to materials → `cad.import_dxf(...)` writes into a template copy → `load_slope_data` opens it.

---

## 9. File / Project Lifecycle

- **New** → copy standard template (`docs/inputs/input_template.xlsx`) into a working doc; user fills via forms; `save_slope_data_to_xlsx` writes it.
- **Open** → `load_slope_data`; render Inputs view. Auto-load `{stem}_mesh.json` / `{stem}_seep.csv` if present (already handled by `load_slope_data`).
- **Save / Save As** → `save_slope_data_to_xlsx` (new function, §5.1). Preserve the existing format.
- **Recent files**, dirty-state prompt on close, and a per-project sidecar for mesh/seep/fem outputs following the existing `{stem}_*.json/.csv` convention.

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
  resources/          # icons, bundled blank template
```

- The app imports the engine as `import xslope` — no copy of the engine, single source of truth.
- Ship **both** distribution channels (decided): `pip install "xslope[gui]"` (pulls PySide6) for Python users, **and** packaged native installers (`.dmg` / `.msi`) for non-Python users.
- Entry point `xslope-studio` (console_script) for dev launch.

---

## 11. Phased Roadmap

**Phase 0 — Engine prerequisite** ✅ **DONE**
- `save_slope_data_to_xlsx(slope_data, path, template=…)` implemented in `xslope/fileio.py` and verified by a round-trip test across 13 representative files (all input categories). See §5.1. *(Recommended follow-up: fold the round-trip check into the repo's `run_tests.py` so it stays green.)*

**Phase 1 — Skeleton + read-only viewer**
- PySide6 app shell, menus, dockable panels, embedded Matplotlib canvas with zoom/pan.
- Open Excel → render Inputs view via `plot_inputs`. Recent files, log pane.

**Phase 2 — Editing**
- Editors for all input categories; mutate `slope_data` → auto re-render.
- Save / Save As / New (from template). Document dirty-state handling.

**Phase 3 — LEM analysis**
- Run-options dialog; worker-thread execution; Search + Solution + Reliability result tabs.
- Single-surface, auto-search, rapid drawdown, reliability.

**Phase 4 — Meshing + Seepage + FEM**
- Meshing dialog; Seepage run + result view; FEM single/SSRM + result views; progress/cancel.

**Phase 5 — Canvas selection + Display Options**
- Pick/double-click-to-edit; `StyleConfig` + Display Options dialog (visibility/colors/styles).

**Phase 6 — DXF + Smart editing + polish**
- DXF import wizard and export-current-view; optional coincident smart-editing; undo/redo.

**Phase 7 — Packaging & distribution**
- PyInstaller or Briefcase native installers (`.dmg`/`.msi`); macOS code-signing/notarization; bundle gmsh for FEM; CI build matrix.

---

## 12. Decisions

**Settled**

- **GUI stack:** **PySide6 (Qt for Python) + embedded Matplotlib**, single-process, calling `xslope` directly. Chosen for plotting reuse and engine symmetry (see §2 trade study). The web/Electron path was evaluated and set aside — its costs (image-based interactivity *or* a duplicated JS plotting layer, plus a serialization bridge) aren't justified by the desktop-only goal, and a FastAPI-wrapped web version remains possible later (§2b).
- **Application name:** **XSlope Studio**.
- **Code location:** a **`studio/`** subfolder inside the existing repo (see §10).
- **v1 scope:** **LEM first** — open / edit / save + LEM analysis and results in the first usable release; Seepage and FEM follow (Phase 4). The roadmap is already ordered this way.
- **Distribution:** **both** — native `.dmg`/`.msi` installers *and* `pip install xslope[gui]`.

**Still open**

1. Nothing blocking. Remaining choices (icon/branding, exact menu layout, undo/redo granularity) are detail-level and can be settled during implementation.

---

## 13. Key Risks / Watch-items

- **Save fidelity** — the Excel save must perfectly round-trip every category; this is the foundation and the main net-new engine code.
- **Long solves blocking the UI** — must run on worker threads from the start; FEM/SSRM also need cancellation.
- **Restyling existing plots** — the hardcoded styles in `plot_*` mean Display Options needs either added style kwargs or a GUI-side renderer; scope carefully.
- **Packaging heavy deps** — gmsh (FEM) and the scientific stack inflate installer size and complicate signing.
- **Coincident smart-editing** — genuinely tricky geometry bookkeeping; keep it opt-in and late.
```