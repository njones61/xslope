# DXF Import/Export

XSLOPE can exchange slope geometry with CAD software through the **DXF** format,
using the `xslope.cad` module. This is most useful for **polygon-based geometry**:
material zones are complex closed shapes that are tedious to type into the input
template but are exactly what a CAD drawing already contains.

There are two directions:

- **Import** — read material-zone polygons from a DXF and write them to the
  `polygons` sheet of an input template (`import_dxf`).
- **Export** — write a complete XSLOPE model out to a layered DXF for
  documentation, sharing, or editing in CAD (`export_dxf`).

!!! note "Import is polygons-only"
    Import reads **material zones** and nothing else. Material zones are the one
    feature where CAD geometry saves real effort. Piezometric lines, distributed
    loads, and search circles are quick to enter directly in the template, so they
    are not imported — enter them on their respective sheets after importing the
    geometry. (Export, by contrast, writes the *full* model.)

## The layer convention

XSLOPE organizes a DXF by **layer**. Each material zone is a closed polyline on a
layer **named after its material** (e.g. `CLAY_CORE`, `BEDROCK`, `SAND`). On
import, the layer name becomes the material name; all polylines on the same layer
get the same material. Other model features, when exported, go on a fixed set of
**reserved layer names**:

| Model feature | DXF entity | Layer | Reserved? |
|---------------|-----------|-------|-----------|
| Material zones | closed `LWPOLYLINE` | one per material (e.g. `CLAY_CORE`) | no — material name |
| Profile lines | open `LWPOLYLINE` | `PROFILE_<material>` | yes |
| Failure surface | `LWPOLYLINE` | `FAILURE_SURFACE` | yes |
| Search circles | `CIRCLE` + `POINT` | `SEARCH_CIRCLES` | yes |
| Reinforcement lines | `LINE` | `REINFORCEMENT` | yes |
| Distributed loads | `LWPOLYLINE` | `DLOADS` | yes |
| Piezometric / water table | `LWPOLYLINE` | `PIEZO` | yes |

Layer names are uppercased with spaces and illegal characters replaced by
underscores (DXF disallows `< > / \ " : ; ? * | = ` `` ` ``). On import, any layer
matching a reserved name is ignored (import is polygons-only); **every other layer
is treated as a material zone**.

![Placeholder: a slope cross-section in CAD with each material zone on its own named layer (CLAY_CORE, SAND, BEDROCK, …).](images/dxf_cad_layers.png)

*Placeholder — screenshot of a layered cross-section open in CAD.*

## Exporting a model to DXF

`export_dxf` writes whatever is present in the loaded model to a layered DXF:

```python
from xslope.fileio import load_slope_data
from xslope.cad import export_dxf

slope_data = load_slope_data("my_slope.xlsx")
export_dxf(slope_data, "my_slope.dxf")          # DXF version R2010 by default
```

Material zones are written as closed polylines on per-material layers; profile
lines, search circles, reinforcement, distributed loads, and piezometric lines are
written on their reserved layers (see the table above). The result opens directly
in any CAD program.

![Placeholder: the exported DXF opened in CAD, showing filled/outlined material zones and the reserved feature layers.](images/dxf_export_example.png)

*Placeholder — the exported `my_slope.dxf` in CAD.*

## Importing polygons from a DXF

### Step 1 — Organize the CAD drawing

Draw each material zone as a **closed polyline** on a **layer named after the
material**. All polylines on a layer are assigned to the same material. The
importer is tolerant of common real-world variations:

- both lightweight `LWPOLYLINE` and heavyweight `POLYLINE` entities,
- **arc segments** (polyline bulges) — tessellated to straight segments,
- **unclosed** polylines — closed automatically (with a warning),
- a zone drawn as **loose `LINE` segments** — stitched into a ring by shared
  endpoints.

### Step 2 — Review the extracted zones

Before importing, print a summary to confirm the extraction — that the right zones
were found, the geometry looks correct, and nothing is missing or duplicated:

```python
from xslope.cad import summarize_dxf

summarize_dxf("my_slope.dxf")
```

```
Poly ID │ Layer Name         │ Vertices │       Area
────────┼────────────────────┼──────────┼───────────
      1 │ EMBANKMENT         │        4 │    2000.00
      2 │ FOUNDATION         │        4 │    4250.00

Unique material layers (2): EMBANKMENT, FOUNDATION
```

![Placeholder: plot of the extracted material-zone polygons, colored by layer.](images/dxf_import_zones.png)

*Placeholder — extracted zones plotted by layer for visual validation.*

### Step 3 — Import into a template

`import_dxf` copies an input template, writes the extracted polygons to the
`polygons` sheet, and seeds the `mat` sheet with the layer names as material
names:

```python
from xslope.cad import import_dxf

result = import_dxf(
    "my_slope.dxf",                          # source DXF
    "docs/inputs/input_template.xlsx",       # template to copy
    "my_slope.xlsx",                         # output input file
)
print(result["layer_to_mat"])   # {'EMBANKMENT': 1, 'FOUNDATION': 2}
print(result["warnings"])       # e.g. ["open polyline on layer 'SAND' auto-closed"]
```

By default each unique layer (in first-appearance order) is assigned material IDs
`1, 2, 3, …`. To control the mapping, pass `material_map`:

```python
import_dxf(
    "my_slope.dxf", "docs/inputs/input_template.xlsx", "my_slope.xlsx",
    material_map={"EMBANKMENT": 1, "FOUNDATION": 2},
)
```

For finer control, the underlying reader is also available directly:

```python
from xslope.cad import dxf_to_polygons

polygons, warnings = dxf_to_polygons("my_slope.dxf")   # optional: layers=[...]
# polygons -> [{'coords': [(x, y), ...], 'layer': 'EMBANKMENT'}, ...]
```

### Step 4 — Fill in material properties

Import seeds only the material **names**. Open the generated `my_slope.xlsx`, go
to the `mat` sheet, and fill in the strength and other properties for each
material. The file is not analysis-ready until the material properties are
provided.

![Placeholder: the imported input file opened in XSLOPE via plot_inputs(), showing filled material zones and the hatched domain base.](images/dxf_import_result.png)

*Placeholder — the imported geometry rendered with `plot_inputs()`.*

## Sample files

A sample exported DXF (a two-material slope on a sloping base) is provided for
experimentation:

- [sample_polygons.dxf](files/sample_polygons.dxf) — two material zones
  (`EMBANKMENT`, `FOUNDATION`) on separate layers.

Try it:

```python
from xslope.cad import summarize_dxf, import_dxf
summarize_dxf("sample_polygons.dxf")
import_dxf("sample_polygons.dxf", "docs/inputs/input_template.xlsx", "sample.xlsx")
```

## Notes and limitations

- **Polygons only on import.** Other features (piezo, loads, circles) are entered
  in the template; export writes them, but import ignores their reserved layers.
- **Material properties are not in the DXF.** Import seeds names only; properties
  are filled in on the `mat` sheet afterward.
- **Round-trip.** Exporting a model and re-importing it reproduces the material-zone
  geometry exactly (feature layers are written by export but ignored on import).
- **Format.** DXF only for now. DWG requires Autodesk libraries or the ODA File
  Converter and may be added based on demand.
