# Running Analyses

Studio runs the same LEM, seepage, and FEM analyses as the `xslope` library, each
behind a run-options dialog. Long solves run on a background thread, so the window
stays responsive, with live output in the [Log pane](interface.md#the-log-pane)
and a Cancel button.

Pick the analysis type with the **Mode** selector, then click **Run** (its label
follows the mode: **Run LEM…**, **Run Seep…**, **Run FEM…**).

---

## Limit equilibrium (LEM)

In **LEM** mode, **Run LEM…** opens a dialog with:

![Run LEM dialog](images/analysis_run_lem_dialog.png)

- **Method** — OMS, Bishop, Janbu, Corps of Engineers, Lowe & Karafiath, Spencer,
  or Morgenstern–Price.
- **Analysis** — single surface, automated search, or reliability.
- **Surface** — circular or non-circular (shown only when the file has both).
- **Number of slices**, the **rapid drawdown** flag, and a **diagnostic** toggle.
- **Composite surfaces** — lets a circle deeper than the bottom of the model be
  truncated at it and run along the base between the two crossings (see
  [Composite Failure Surfaces](../lem/overview.md#composite-failure-surfaces)).
  Off by default, and only available for circular surfaces. Turn it on when the
  base of the model is a real impenetrable boundary — bedrock, or a weak seam
  resting on it — because the critical mechanism there follows the base and no
  ordinary circle can reach it. Leave it off when the bottom of the model is
  simply how deep you chose to look.
- **Grid search** — seeds the circular search from an automatic
  grid-and-tangent sweep instead of (only) the circles sheet (see
  [Grid Seeding](../lem/search.md#grid-seeding-global-search)). Off by default,
  available for circular auto-search and reliability. Turn it on to protect
  against the local-minimum trap of a single starting circle — a seed in the
  wrong family can converge 20% or more too high with no warning — or when you
  have no idea where the critical circle is (the circles sheet may even be
  empty). It reports the most critical surface anywhere in the model; leave it
  off to interrogate a specific mechanism with your own circles.
- **Search tolerances** (`fs_tol`, `tol`, `max_iter`) — enabled for the
  search-driven analyses.

The result depends on the analysis type:

| Analysis | Result tabs |
| --- | --- |
| **Single surface** | **LEM · Solution** — the surface with slices, base stresses, and thrust line. |
| **Automated search** | **LEM · Search** (all trial surfaces + critical + search path) *and* **LEM · Solution** (the critical surface). |
| **Reliability** | **LEM · Reliability** plus the **Solution** for the most-likely-value surface. A determinate progress bar tracks the `1 + 2N` searches. |

![LEM Search view](images/analysis_lem_search.png)

![LEM Solution view](images/analysis_lem_solution.png)

Both circular and non-circular surfaces are supported for single solves and
searches. Search iteration progress streams to the Log pane, and a search (or a
reliability run) can be **cancelled** from the status bar.

![LEM Reliability view](images/analysis_lem_reliability.png)

---

## Building a mesh

Seepage and FEM run on a finite-element mesh, which you build explicitly. In
**Seepage** or **FEM** mode, **Build Mesh** opens a dialog with:

![Build Mesh dialog](images/analysis_build_mesh_dialog.png)

- **Element type** — `tri3`, `tri6`, `quad4`, `quad8`, or `quad9`.
- **Target size** — entered directly, or auto-sized as the slope width divided by a
  number of divisions.

The mesh is built on a background thread (it includes reinforcement and pile
constraint lines, so it serves FEM too), shown in a **Mesh** tab, and written to a
`{stem}_mesh.json` sidecar. Seep/FEM **Run** stays disabled until a mesh exists; a
geometry edit that invalidates the mesh clears it and re-gates Run.

![Mesh view](images/analysis_mesh_view.png)

!!! note "Meshing needs gmsh"
    Mesh generation uses **gmsh**, which is installed by the `fem` extra
    (`pip install "xslope[gui,fem]"`). See [Installation](index.md#installation).

---

## Seepage

In **Seepage** mode, **Run Seep…** opens a dialog with the convergence tolerance,
the **BC set** (1 or 2), the **variable to plot** (head, pressure, velocity
magnitude, gradient magnitude), contour levels, and toggles for flow lines,
vectors, fill, and the phreatic surface.

![Run Seep dialog](images/analysis_run_seep_dialog.png)

The run produces two tabs — **Seep · Data** (mesh + boundary conditions) and
**Seep · Solution** (the chosen variable, contours, phreatic surface, flow lines /
vectors) — and writes the solution to `{stem}_seep.csv` (`_seep2.csv` for BC 2).
The convergence trace streams to the Log pane. Each BC set keeps its own tab pair,
so rapid-drawdown problems show BC 1 and BC 2 together.

![Seep Data view](images/analysis_seep_data.png)

![Seep Solution view](images/analysis_seep_solution.png)

---

## Finite element (FEM)

In **FEM** mode, **Run FEM…** offers a **single trial** or an **SSRM** run (the
Shear Strength Reduction Method), with `F` (or `F_min`/`F_max`), a tolerance, and
the failure criterion.

![Run FEM dialog](images/analysis_run_fem_dialog.png)

The run produces **FEM · Data** (mesh + boundary conditions + reinforcement) and
**FEM · Results** (deformation, shear strain, displacement vectors). An SSRM run
reports the factor of safety and can be **cancelled** mid-run. The solution is
exported alongside the model so it can be restored on the next Open without
re-solving.

![FEM Data view](images/analysis_fem_data.png)

![FEM Results view](images/analysis_fem_results.png)

---

## Display options per view

Each result view has its own **Display** panel (in the left dock) exposing the
options the underlying plot accepts — slice numbers and seep contours on the LEM
solution; nodes / labels / padding on the mesh; variable, levels, vector scale, and
flow-line toggles on the seep solution; plot type and deformation scale on FEM
results; legend column layout on every view. Changing an option re-renders the
**cached** result instantly — there's no re-solve. See
[The Display dock](interface.md#the-display-dock).

![A Display panel](images/analysis_display_panel.png)

---

## Exporting views: images and DXF

Every canvas has a **Save…** button that exports the **current view**:

![Save view dialog](images/analysis_save_dialog.png)

- **Image** — PNG, PDF, or SVG. PNG prompts for a DPI; vector formats don't. The
  figure is saved at its true inch size, so resolution is independent of the
  on-screen zoom.
- **DXF (rendered view)** — the drawn picture as a layered DXF, good for dropping
  into a CAD document. This is *lossy* as a re-import source; for that, use the
  structured geometry export below.

---

## DXF import and export

Studio exchanges model geometry with CAD through DXF.

**Export geometry** — **File → Export Geometry (DXF)…** writes the structured,
layered DXF: material zones on per-material layers, and profile lines, circles,
reinforcement, distributed loads, and piezo lines on reserved feature layers. This
is the clean companion the importer reads (distinct from the per-view *rendered*
DXF above).

**Import** — **File → Import DXF…** opens a wizard that lists every layer in the
drawing and lets you map each one to an input feature — material zone, profile
line, piezo line, distributed load, reinforcement, failure circles, or *ignore* —
with a material column for zones and profiles.

![DXF import wizard](images/analysis_dxf_import_wizard.png)

Defaults are seeded from xslope's own export layer names and the geometry kind, but
you can override anything, so a DXF drawn in external CAD with arbitrary layer names
maps just as well. Geometry populates the features; non-geometric properties (load
magnitudes, reinforcement strengths, material properties, circle depth) come in as
editable **placeholders** to fill in afterward. The import replaces the current
project (you're prompted to discard first) and is left unsaved so you can complete
it and Save As.

!!! info "DXF layer conventions"
    The layer naming and entity conventions are shared with the library's
    `xslope.cad` module — see [DXF Import/Export](../usage/dxf.md) for the full
    layer table and format details.

!!! note "DXF support needs ezdxf"
    Reading and writing DXF uses the **ezdxf** package (installed with the `gui`
    extra). If it's missing, the import/export actions show an actionable install
    message.

---

## GeoStudio (SLOPE/W) import and export

Studio also exchanges whole models — not just geometry — with GeoStudio SLOPE/W.

**Import** — **File → Import GeoStudio (SLOPE/W)…** reads a `.gsz`. Unlike DXF, there
is no mapping wizard to work through: a `.gsz` already knows what its geometry means, so
material zones, strengths, water conditions and the seismic coefficient all arrive
identified. The one prompt is **which analysis** to import, because a GeoStudio file
usually holds several over the same geometry — and they can differ in *materials*, not
just in slip surface, so the choice changes the model you get.

**Export** — **File → Export to GeoStudio (SLOPE/W)…** writes the current model out as a
`.gsz`. It needs a polygon-based model (material zones), since a profile-line model has
no regions to map onto.

Both directions replace nothing silently: whatever cannot cross the format boundary —
SLOPE/W's search definition, reinforcement, piles, loads, non-Mohr-Coulomb strengths — is
listed in a notes dialog and in the Log pane, so you know exactly what to re-create by
hand.

!!! info "Units and what survives the trip"
    A `.gsz` carries no unit-system field, so XSLOPE infers it from the unit weight of
    water and refuses to guess when it's neither metric nor imperial. See
    [GeoStudio Import/Export](../usage/geostudio.md) for the full mapping table, the
    per-analysis materials wrinkle, and the limits of export.
