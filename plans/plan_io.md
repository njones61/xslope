# Plan: File I/O and Interchange

xslope reads and writes several kinds of files. This document is the home for
all file-I/O planning — the Excel template, CAD geometry interchange (DXF/DWG),
and (eventually) native project files from other geotechnical software.

The polygon-geometry feature (see [`plan_polygons.md`](plan_polygons.md))
introduced the first new I/O surface beyond the Excel template (DXF import of
polygons). The CAD section below was originally drafted there and has been moved
here so the I/O concerns live in one place.


## 1. Module Architecture

Today all file handling lives in **`fileio.py`** (Excel template parsing). As we
add CAD and vendor-format support, the cleanest split is by *concern*:

| Module | Concern | Status |
|--------|---------|--------|
| `fileio.py` | Excel input template (xlsx) read/write | exists |
| `cad.py` | Neutral CAD geometry interchange (DXF; later DWG/DGN) | planned (§3) |
| `interop.py` | Native project files of other software packages | future stub (§4) |

**Why separate `cad.py` and `interop.py`?** CAD is "here are some polygons/lines"
— pure geometry, no analysis semantics. A vendor project file (SLOPE/W, Slide2,
GEO5, PLAXIS) is "here is a whole model in someone else's vocabulary" — materials,
methods, pore-pressure models, and results that must be *translated* onto xslope's
data structures, per vendor. Different problem, different module.

### 1.1 Possible future package layout

If file I/O grows past two or three formats, promote these flat modules into a
package so each vendor format gets its own file:

```
xslope/io/
    __init__.py      # re-exports the public API (import_*/export_*)
    template.py      # current fileio.py (xlsx)
    cad.py           # DXF / DWG / DGN
    interop/
        __init__.py  # uniform import_*/export_* across vendors
        slopew.py
        slide2.py
        geo5.py
```

**Not now.** Flat modules (`fileio.py`, `cad.py`, future `interop.py`) are fine
until the vendor list grows. Refactor to the package only when `interop` holds
more than ~two formats.


## 2. Template I/O (`fileio.py`)

**Status**: exists. `load_slope_data()` parses the xlsx template into `slope_data`.

**Planned work** (tracked in detail in `plan_polygons.md`):
- Parse the new `polygons` sheet → `slope_data['polygons']` (plan_polygons §3).
- Derive ground surface / domain polygon from polygons (plan_polygons §4).

**Writer** (done): a surgical, formula-safe xlsx writer now lives in
**`fileio.py`** as `write_cells_to_xlsx(filepath, updates)` plus the
`cell_ref(row, col)` helper. It edits cell XML in place, forces
`fullCalcOnLoad`, and drops the stale `calcChain.xml` so populated sheets don't
open blank in Excel (see the `xlsx calcChain/recalc` and `no overwrite formulas`
notes in project memory). The `poly_test/` builder scripts (`_build_seep_polys.py`,
`_build_lost_lake_poly.py`) now import it rather than carrying their own copies.
The CAD importer (§3) reuses the same function.

> Requires the `zip` CLI (used to replace single members inside the xlsx archive
> in place). Target only value/precedent cells — never cells holding formulas.


## 3. CAD Import / Export (`xslope/cad.py`)

Introduce a new module, **`xslope/cad.py`**, that converts between DXF files and
xslope input templates in **both directions**.

**Module name**: `import.py` is not an option — `import` is a Python reserved word
(`from xslope import import` is a syntax error). `cad.py` is preferred over
`dxf.py` because §3.5 anticipates DWG/DGN support, and over `import_export.py` for
brevity. The intended public API reads cleanly:

```python
from xslope.cad import import_dxf, export_dxf
```

The module has two halves:
- **Import** (`import_dxf`): read a DXF cross-section → polygons with materials →
  write the `polygons` sheet of an input template.
- **Export** (`export_dxf`): write an input template back out to DXF — not just
  polygons, but profile lines, failure surfaces, search circles, reinforcement,
  distributed loads, and piezo lines — each on its own named layer.

Round-tripping (export → import) is the basis for testing: a template exported to
DXF and re-imported should reproduce the same geometry.

### 3.1 Test Fixtures (already generated)

`poly_test/` contains DXF fixtures generated from the polygon-sheet xlsx problems,
plus two reusable generator scripts:

- **`_build_dxf_from_poly.py`** — *clean* fixtures: one closed `LWPOLYLINE` per
  zone, one DXF layer per material (layer name = material name, uppercased, spaces
  → underscores). Outputs `xslope_{sea_trench,earth_dam,johnson_res,levee,lost_lake}_poly.dxf`.
  The lost-lake file has 12 zones across 8 material layers; the levee/sea-trench
  files exercise the **multiple-polygons-share-one-layer** case.
- **`_build_dxf_messy.py`** — *defective* variants (all from sea-trench geometry)
  that exercise the importer's recovery paths:
  `xslope_sea_trench_messy_{unclosed,polyline,line_segments,arcs}.dxf`.

These are the canonical inputs for importer development and regression tests. Each
generator round-trips its output back through `ezdxf` as a built-in smoke test.
The clean-fixture logic in `_build_dxf_from_poly.py` is the seed for `export_dxf`'s
polygon path.

### 3.2 Import: DXF → input template

DXF is the Phase 1 format: open, well-documented, mature `ezdxf` support, and the
common exchange format for cross-sections (AutoCAD, MicroStation, etc.). Every
entity belongs to a **layer** (e.g. "CLAY_FILL", "BEDROCK", "SAND_LENS");
geotechnical drawings organize material zones by layer, so layer names drive
material assignment.

**Background — commercial patterns**: the commercial packages use one of two
*targeted* approaches: (1) import one feature at a time (separate command per
feature type), or (2) require features to follow a reserved layer-naming
convention and route entities by layer. Importing an arbitrary, unstructured DXF
of a full cross-section is not realistic in either case.

**Scope & strategy (decided)**:

- **Import is polygons-only.** `import_dxf` reads *material zones* and nothing
  else. Rationale: polygons are the only feature where CAD geometry saves real
  work — a material zone is a complex closed shape that is painful to type. Piezo
  lines (~5 points), distributed loads (a few values), and search circles (3
  numbers) are faster to enter directly in the template than to round-trip
  through a precisely-structured CAD file, so they stay as spreadsheet entry.
- **Per-feature primitive.** `import_dxf` is a thin wrapper over a
  `dxf_to_polygons(path, layers=...)` primitive. Keeping that primitive separate
  leaves the door open to add `dxf_to_piezo`, `dxf_to_dloads`, etc. later (the
  per-feature pattern) without restructuring — but we do **not** build them now.
- **Material identity = bare layer name.** Layer `CLAY_CORE` → a material named
  `CLAY_CORE`; the user confirms or remaps in the validation step (Step 3). This
  matches the fixtures we already generate. Reserved feature-layer names (used by
  export, §3.3) are the only ones treated specially; everything else is a
  material zone. (Collision risk — a material literally named `PIEZO` — is
  negligible and surfaces in the validation table anyway.)
- **Export stays full** (§3.3) and **shares one layer convention** with import,
  so any future per-feature importer is the exact inverse of what export writes.

**Workflow**: a standalone utility that reads a DXF, helps the user validate and
assign materials, and writes the `polygons` sheet.

**Step 1 — Read polygons (robustly).** Real exports are messy; the reader must
handle more than the happy path. The messy fixtures pin down each case:

- **Both polyline types** — read `LWPOLYLINE` *and* heavyweight `POLYLINE`. (`messy_polyline.dxf`)
- **Arc bulges** — LWPOLYLINE segments can carry a bulge (curved interface). Use
  `entity.flattening(sag)` to tessellate arcs into vertices; naively reading
  `get_points()` silently drops the curvature and corrupts the computed area. (`messy_arcs.dxf`)
- **Unclosed rings** — if a polyline isn't closed (no close flag, no repeated end
  vertex), auto-close it and emit a warning (per plan_polygons §3.2). (`messy_unclosed.dxf`)
- **Loose LINE segments** — a zone may be drawn as disconnected `LINE` entities
  never joined into a polyline. Attempt to stitch segments sharing endpoints into
  closed rings; if stitching is ambiguous, reject with a clear message. (`messy_line_segments.dxf`)

```python
def read_dxf_polygons(dxf_path, arc_sag=0.05):
    """Read closed zones from DXF as {coords, layer}, handling LWPOLYLINE,
    POLYLINE, arc bulges (flattened), and unclosed rings. Loose LINE entities
    are stitched into rings separately (or reported)."""
    import ezdxf
    doc = ezdxf.readfile(dxf_path)
    msp = doc.modelspace()
    polygons, warnings = [], []
    for e in msp.query("LWPOLYLINE POLYLINE"):
        coords = [(p.x, p.y) for p in e.flattening(arc_sag)]  # tessellates bulges
        if not e.is_closed:
            warnings.append(f"open polyline on layer '{e.dxf.layer}' auto-closed")
        polygons.append({"coords": coords, "layer": e.dxf.layer})
    # loose LINE entities → stitch shared endpoints into rings (or report)
    ...
    return polygons, warnings
```

> The original draft of this reader read only `LWPOLYLINE` via `get_points()` and
> `entity.closed`. Against the fixtures that would skip the `POLYLINE`/`LINE` files
> entirely, drop the arc bulges, and reject the unclosed file instead of
> recovering — hence the expanded reader above.

**Step 2 — Display for validation.** Plot all polygons with layer-name labels and
print a summary table so the user can confirm the import (correct zones extracted,
geometry sound, nothing missing or duplicated):

```
Poly ID │ Layer Name       │ Vertices │ Area
────────┼──────────────────┼──────────┼────────
1       │ CLAY_FILL        │ 6        │ 450.2
2       │ SAND_FOUNDATION  │ 8        │ 320.1
3       │ BEDROCK          │ 5        │ 890.0
```

**Step 3 — Material assignment.** List unique layer names; the user maps each to a
material ID. All polygons on a layer get that material — the natural mapping since
CAD layers correspond to material zones:

```
Layer Name       │ Assign Material ID
─────────────────┼───────────────────
CLAY_FILL        │ 1
SAND_FOUNDATION  │ 2
BEDROCK          │ 3
```

**Step 4 — Write the template.** Populate the `polygons` sheet (Poly ID, Mat ID,
vertices). Optionally seed the `mat` sheet with material names derived from the
layer names (e.g. layer "CLAY_FILL" → material named CLAY_FILL with placeholder
properties the user then fills in), giving a head start over blank rows. Reuse the
surgical, formula-safe xlsx writer now in `fileio.py` (§2):
`fileio.write_cells_to_xlsx`.

### 3.3 Export: input template → DXF

The reverse direction turns a complete xslope model into a layered CAD drawing —
useful for documentation, sharing with CAD users, and round-trip testing. Each
entity class goes on its own layer:

| Source data | DXF entity | Layer(s) | Reserved? |
|-------------|-----------|----------|-----------|
| Polygons (`polygons` sheet) | closed LWPOLYLINE | one per material (e.g. `CLAY_CORE`) | no — material name |
| Profile lines | LWPOLYLINE (open) | `PROFILE_<mat>` | yes |
| Failure surface(s) | LWPOLYLINE | `FAILURE_SURFACE` | yes |
| Search circles (centers + radii) | CIRCLE + POINT | `SEARCH_CIRCLES` | yes |
| Reinforcement lines | LINE | `REINFORCEMENT` | yes |
| Distributed loads | LWPOLYLINE | `DLOADS` | yes |
| Piezo / water table | LWPOLYLINE | `PIEZO` | yes |

**Shared layer convention.** The reserved layer names above are the single
vocabulary shared between export and import. On import (§3.2) any layer matching a
reserved name is a non-polygon feature; **every other layer is a material zone
named by the layer**. The current polygons-only importer simply *ignores* the
reserved feature layers; if per-feature importers are added later, they consume
exactly these names — making import the precise inverse of export. Layer naming:
uppercase, spaces/illegal DXF chars → underscores (DXF disallows `< > / \ " : ; ? * | = ``).

The polygon path is exactly the clean-fixture logic from `_build_dxf_from_poly.py`,
generalized to the other entity types.

```python
def export_dxf(slope_data, dxf_path, version="R2010"):
    """Write polygons, profile lines, failure surfaces, search circles,
    reinforcement, distributed loads, and piezo lines to a layered DXF."""
    ...
```

### 3.4 Implementation Order (CAD)

1. `export_dxf` polygon path — lift directly from `_build_dxf_from_poly.py`.
2. `dxf_to_polygons(path, layers=...)` primitive + `import_dxf` wrapper — read
   material zones from the clean fixtures, ignoring reserved feature layers;
   validation plot/table; write the `polygons` sheet. (Polygons-only — no
   per-feature importers; see §3.2 scope decision.)
3. `import_dxf` robustness — POLYLINE, bulge flattening, auto-close, LINE stitching
   (the messy fixtures).
4. `export_dxf` remaining entity types (profiles, failure surfaces, circles, etc.)
   on their reserved layers.
5. Round-trip test: export → re-import reproduces the polygon geometry (feature
   layers are written by export but ignored by the polygons-only import).

### 3.5 Future: DWG, DGN

DWG requires Autodesk libraries or the ODA File Converter. DGN is niche. Add based
on user demand. Naming the module `cad.py` (rather than `dxf.py`) leaves room for
these formats without a rename.


## 4. Interop: Other Software Packages (`xslope/interop.py`) — STUB

> **Status: not started.** Captured here so the design intent isn't lost.

Read (and possibly write) the *native project files* of other geotechnical
software, mapping their model semantics onto xslope's `slope_data`. Unlike CAD
(§3), this is more than geometry — materials, analysis methods, pore-pressure
options, and sometimes results must be translated per vendor.

Candidate formats (prioritize by user demand):

| Software | File(s) | Notes |
|----------|---------|-------|
| GeoStudio SLOPE/W | `.gsz` (zipped XML) | Common; XML is parseable without vendor libs |
| Rocscience Slide2 | `.slim` / `.sli` | Format may need reverse-engineering |
| GEO5 | GEO5 project files | DXF interchange already covered by §3 |
| PLAXIS | project / command files | Primarily FEM; relevant to xslope's FEM side |

**Open questions for this section:**
- Import-only first, or bidirectional? (Import is the higher-value direction.)
- How much semantic mapping vs. geometry-only? A geometry-only path could route
  through `cad.py` concepts; full model import needs a per-vendor translator.
- Where do *results* fit — import another tool's FS/failure surface for comparison?

**TODO**: pick a first target format (likely SLOPE/W `.gsz`, since it is zipped
XML and widely used), prototype a reader, and define the vendor→xslope material
and method mapping.


## 5. Open Questions

1. **Writer promotion** (§2): finalize the public signature for the xlsx writer in
   `fileio.py` and migrate the `poly_test/` scripts to use it.
2. **Package vs. flat modules** (§1.1): decide the trigger for refactoring into
   `xslope/io/` — proposed: when `interop` gains a second format.
3. **CAD round-trip fidelity**: arc flattening (`arc_sag`) is lossy; decide the
   default tolerance and whether to preserve true arcs on export for entities that
   are genuinely curved.
