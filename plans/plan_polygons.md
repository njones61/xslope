# Plan: Defining Slope Geometry with Polygons


## 1. Problem Overview

xslope currently uses **profile lines** to define slope geometry. Each profile line represents a material boundary, ordered top-to-bottom. The topmost profile is the ground surface, lower profiles define layer interfaces, and a horizontal `max_depth` parameter serves as the bottom boundary.

This approach works well for simple layered slopes but becomes limiting for:
- Irregular bedrock surfaces or dipping strata
- Lenticular deposits (lens-shaped inclusions)
- Complex fill geometries (e.g., dam cross-sections with multiple zones)
- Geometries imported from CAD software, which are typically drawn as closed regions

**Proposed change**: Allow users to define slope geometry using **polygons** as an alternative to profile lines. Each polygon represents a material zone (a closed region with assigned material properties). Users can choose either input method depending on their preference and the complexity of the geometry.


## 2. Current Architecture

### 2.1 Profile Lines Data Flow

```
Excel "profile" sheet
    ↓
load_slope_data()  (fileio.py:139-212)
    ├─→ max_depth          (Row B2 — horizontal bottom boundary)
    ├─→ profile_lines      (list of dicts: {'coords': [...], 'mat_id': int})
    └─→ build_ground_surface()  (fileio.py:24-97)
         └─→ ground_surface (LineString — highest points from all profiles)

slope_data
    ├─→ generate_slices()  (slice.py)
    │    ├─→ Ground surface: explicit from slope_data['ground_surface']
    │    ├─→ Slice x-coords: from profile points, loads, failure surface
    │    ├─→ Material per layer: profile_lines[i].mat_id → materials[mat_id]
    │    └─→ Layer heights: between adjacent profile lines at each x
    │
    ├─→ build_polygons()   (mesh.py:1235-1700)
    │    ├─→ Converts profile lines → closed polygons with mat_id
    │    ├─→ Bottom of lowest polygon = max_depth
    │    └─→ Used by mesh generation (seepage + FEM)
    │
    └─→ circular_search()  (search.py)
         ├─→ Depth constrained by: max_depth (scalar)
         └─→ Circle validity: must intersect ground_surface exactly twice
```

### 2.2 Key Data Structures

**Profile lines** (`slope_data['profile_lines']`):
```python
[
    {'coords': [(x1,y1), (x2,y2), ...], 'mat_id': 0},  # topmost = ground surface
    {'coords': [(x1,y1), (x2,y2), ...], 'mat_id': 1},  # second layer interface
    ...
]
```

**Polygons from `build_polygons()`** (used for mesh generation):
```python
[
    {'coords': [(x1,y1), (x2,y2), ...], 'mat_id': 0},  # closed polygon, material zone
    ...
]
```

**Materials** (`slope_data['materials']`):
```python
[
    {'name': 'Clay', 'gamma': 18.0, 'c': 10.0, 'phi': 25.0, 'E': 5000, ...},
    ...
]
```

### 2.3 How Materials Are Assigned to Slices

In `generate_slices()` (slice.py:799-845), for each slice:

1. Iterate over profile lines top-to-bottom
2. For each profile, compute the height $h$ of that layer within the slice (between this profile and the next lower one, clipped to the failure surface)
3. The lowermost layer with $h > 0$ determines the **base material** (controls $c$, $\varphi$, $u$)
4. All layers contribute to slice weight: $W = \sum h_j \cdot \gamma_j \cdot \Delta x$

### 2.4 How the Search Algorithm Uses Boundaries

In `circular_search()` (search.py):

- **Bottom boundary**: `max_depth` (scalar). Circles are clamped: `depth = max(depth, max_depth)`. The depth optimization loop never tests circles below `max_depth`.
- **Left/right boundaries**: Not explicitly checked. Instead, `get_sorted_intersections()` finds where the failure surface crosses the ground surface. If fewer than 2 intersection points are found (circle goes off the side), `generate_slices()` returns `success=False` and the circle is assigned `FS = 9999` (rejected).
- **Ground surface**: The circle must intersect `ground_surface` exactly twice to define valid entry/exit points.


## 3. Input: "polygons" Sheet

### 3.1 Sheet Layout

Each polygon is defined by a **Polygon ID**, a **Material ID**, and a series of $(x, y)$ vertices. Multiple polygons are listed sequentially — a new Polygon ID signals the start of a new polygon.

```
Row 1: Sheet title ("Polygon Geometry")
Row 2: Headers
Row 3+: Data

     A          B       C      D
  ──────────────────────────────────
2 │ Poly ID  │ Mat ID │  x   │  y
  ──────────────────────────────────
3 │ 1        │ 1      │ 0.0  │ 30.0
4 │          │        │ 20.0 │ 30.0
5 │          │        │ 40.0 │ 25.0
6 │          │        │ 80.0 │ 25.0
7 │          │        │ 80.0 │ 15.0
8 │          │        │ 0.0  │ 15.0
9 │ 2        │ 2      │ 0.0  │ 15.0
10│          │        │ 80.0 │ 15.0
11│          │        │ 80.0 │ 5.0
12│          │        │ 0.0  │ 10.0
13│          │        │      │
```

### 3.2 Column Definitions

| Column | Field | Description |
|--------|-------|-------------|
| A | Poly ID | Integer identifier for the polygon (new value = new polygon) |
| B | Mat ID | Material ID (1-based, references the `mat` sheet) |
| C | $x$ | X-coordinate of vertex |
| D | $y$ | Y-coordinate of vertex |

- Poly ID and Mat ID only need to appear on the **first row** of each polygon; subsequent rows for the same polygon leave columns A and B blank.
- Reading stops when both columns C and D are empty.
- Polygons are **automatically closed** — the last vertex is connected back to the first. A warning message is printed if the user did not explicitly close the loop.

### 3.3 Design Decisions

- **Winding order**: Users can input vertices in either CW or CCW order. Shapely normalizes to CCW internally, so no detection or enforcement is needed.
- **Nesting is automatic**: If polygon B is geometrically contained within polygon A (`A.contains(B)` in Shapely), B is nested inside A. No parent-child column needed — nesting is determined from geometry.
- **Polygon vs. profile choice**: The `polygons` sheet is an alternative to the `profile` sheet. If the `polygons` sheet exists and contains data, it is used; otherwise, profile lines are used. Both should not be specified simultaneously.
- **Material assignment**: Each polygon has a `mat_id` referencing the `mat` sheet, identical to how profile lines reference materials.

### 3.4 Data Structure

Store as `slope_data['polygons']` — a list of Shapely `Polygon` objects with associated material IDs:

```python
slope_data['polygons'] = [
    {'polygon': Polygon([(x1,y1), ...]), 'mat_id': 0},
    {'polygon': Polygon([(x1,y1), ...]), 'mat_id': 1},
    ...
]
```

### 3.5 Validation

- Each polygon must have at least 3 vertices
- Polygons must be valid (not self-intersecting): `polygon.is_valid`
- No overlapping polygons (same area claimed by two materials) — check with pairwise `intersection().area > tolerance`
- Warn on gaps between polygons (uncovered area within the bounding box)
- Mat ID must reference a valid material in the `mat` sheet


## 4. Derived Geometry

### 4.1 Domain Polygon

The **domain polygon** is the union of all input polygons:

```python
from shapely.ops import unary_union
domain_polygon = unary_union([p['polygon'] for p in slope_data['polygons']])
```

This defines the total extent of the slope model. Failure surfaces must stay within this boundary.

### 4.2 Ground Surface Extraction

The ground surface is the **upper boundary** of the domain polygon. With polygons, this is no longer an explicit input — it must be derived.

**Algorithm:**

1. Compute the domain polygon (union of all polygons)
2. Extract the exterior ring: `domain_polygon.exterior`
3. The ground surface is the portion of the exterior ring that forms the **top boundary** — the segments where $y$ is locally maximum

Practically, this means tracing the exterior from the top-left corner rightward along the upper boundary to the top-right corner. The left and right boundaries are vertical (or near-vertical) segments connecting the top to the bottom.

**Approach**: Extract the exterior coordinates, identify the leftmost and rightmost points at the top of the domain, and take the segment between them that has the highest elevation:

```python
exterior = list(domain_polygon.exterior.coords)
# Find the top-left and top-right transition points
# (where the boundary transitions from vertical side to top surface)
# Trace the top boundary between them → ground_surface LineString
```

This replaces the current `build_ground_surface()` function when polygons are used.

### 4.3 Replacing the Limiting Depth

With profile lines, `max_depth` is a horizontal line — the de-facto bottom boundary. With polygons, the bottom boundary is the **lower boundary of the domain polygon**, which can be irregular (dipping bedrock, stepped foundations, etc.).

**Key insight**: The bottom boundary does not need to be extracted as a separate entity. The domain polygon itself serves as the containment boundary. The failure surface must stay within the domain polygon — that's the only check needed.

`max_depth` is only associated with the profile lines sheet and is not used when polygons are provided. The domain polygon itself defines all boundaries.


## 5. LEM Implementation

### 5.1 Unified Code Path via `profiles_to_polygons()`

The existing `build_polygons()` function in `mesh.py` already converts profile lines to polygons. The cleanest architecture is:

1. If the user provides **polygons**: use them directly
2. If the user provides **profile lines**: convert to polygons using `build_polygons()` early in the pipeline

Then all downstream code (slice generation, search, etc.) works with polygons only — **one code path** instead of two.

This means `generate_slices()` needs a polygon-based implementation, and the profile-based layer height calculation becomes a legacy path (or is removed entirely once polygons are the canonical internal representation).

### 5.2 Polygon-Based Slice Generation

The current `generate_slices()` computes layer heights by iterating over profile lines. With polygons, the approach changes to **geometric intersection**:

For each slice between $x_l$ and $x_r$:

1. **Build a vertical strip**: `strip = box(x_l, -inf, x_r, +inf)`
2. **Intersect each material polygon with the strip**: `zone = polygon.intersection(strip)`
3. **Clip to between ground surface and failure surface**: intersect the zone with the region above the failure surface and below the ground surface
4. **Compute area and centroid** of each clipped zone → gives layer height, weight contribution, and centroid for $y_{cg}$

```python
for i, mat_poly in enumerate(slope_data['polygons']):
    zone = mat_poly['polygon'].intersection(strip)
    clipped = zone.intersection(slice_region)  # between ground and failure surface
    if not clipped.is_empty:
        area = clipped.area
        h = area / dx  # effective height
        gamma = materials[mat_poly['mat_id']]['gamma']
        weight += gamma * area
        # centroid for y_cg weighted average
```

**Base material**: The material of the polygon that contains the **center of the slice base** (the point on the failure surface at $x_c$). This determines $c$, $\varphi$, and pore pressure option — same role as the current "lowermost layer with $h > 0$" logic.

```python
base_point = Point(x_c, y_cb)
for mat_poly in slope_data['polygons']:
    if mat_poly['polygon'].contains(base_point):
        base_mat_id = mat_poly['mat_id']
        break
```

### 5.3 Slice Boundary Points (fixed x-coordinates)

Currently, fixed x-coordinates for slice boundaries come from profile line vertices, distributed load transitions, failure surface points, etc. With polygons:

- **Replace profile vertices with polygon vertices**: extract all x-coordinates from polygon exteriors that lie on or above the failure surface
- All other sources of fixed x-coordinates remain the same (loads, failure surface, piezo intersections)

### 5.4 Search Algorithm and Domain Boundaries

The current search algorithm constrains circles using `max_depth` (scalar) and ground surface intersection count. With polygons:

**Circle validation** uses the same existing logic:

1. **Ground surface intersection** (unchanged): `get_sorted_intersections()` finds entry/exit points. If fewer than 2 intersections → circle goes off the left or right side → rejected (`success=False`).

2. **Bottom boundary** (new): After confirming 2 valid ground surface intersections, check whether the failure surface stays within the domain polygon:
   ```python
   if not domain_polygon.contains(failure_surface):
       return False, "Failure surface extends outside domain"
   ```
   This replaces the scalar `depth >= max_depth` check. A circle that dips below an irregular bedrock surface gets rejected, just as one that goes below `max_depth` does today.

3. **Left/right boundaries**: Not a new problem. A circle that extends past the left or right edge of the domain will only intersect the ground surface once (or not at all), triggering the existing rejection in step 1. No new code needed.

**Performance note**: The `domain_polygon.contains(failure_surface)` check is more expensive than the scalar `depth >= max_depth` comparison. For the search grid with thousands of trial circles, this could matter. Optimizations:
- **Flat bottom fast path**: Detect once during initialization whether the bottom boundary of the domain polygon is horizontal (or approximately so). If it is, use a simple scalar depth check just like the current `max_depth` logic — no geometric containment needed. This applies to profile-line inputs (which always have a flat `max_depth` bottom) and to polygon inputs where the bottom happens to be flat. Only irregular bottom boundaries require the full containment check.
- **Quick reject**: Before the full containment check, compare the circle's lowest point against `domain_polygon.bounds[1]` (minimum $y$ of the bounding box). If the circle is above the bounding box minimum, it's likely contained — skip the expensive check.
- **Prepared geometry**: Use `shapely.prepared.prep(domain_polygon)` for faster repeated containment checks against many failure surfaces.

**Non-circular surface search**: Control points must stay within the domain polygon. Replace the current `y >= max_depth` constraint with `domain_polygon.contains(Point(x, y))` for each movable control point.


## 6. Seepage and FEM

### 6.1 Current Path

Both seepage and FEM currently use `build_polygons()` to convert profile lines into material zone polygons, then pass those polygons to the mesh generator.

### 6.2 Changes

When polygons are provided directly:

- **Skip `build_polygons()`** — the polygons are already in the required format
- Pass `slope_data['polygons']` directly to the mesh generator
- Ensure the polygon format matches what the mesher expects (list of dicts with `'coords'` and `'mat_id'`)

This is straightforward since the mesh generator already works with polygons — we're just removing the conversion step.


## 7. CAD Import / Export

CAD interchange (DXF import of polygons, and DXF export of the full model) is part
of the broader file-I/O design and is documented in **[`plan_io.md`](plan_io.md) §3**.

In brief: a new module **`xslope/cad.py`** provides `import_dxf` (DXF cross-section
→ `polygons` sheet) and `export_dxf` (template → layered DXF). The importer maps
DXF layers to materials and must handle real-world messiness (`POLYLINE` as well as
`LWPOLYLINE`, arc bulges, unclosed rings, loose `LINE` segments). Test fixtures —
clean and deliberately messy — already exist in `poly_test/` (generated by
`_build_dxf_from_poly.py` and `_build_dxf_messy.py`). See `plan_io.md` for the full
reader/writer design, the entity→layer mapping, and the CAD implementation order.


## 8. Special Cases

### 8.1 Nested Polygons (Embedded Zones)

When polygon B is inside polygon A (e.g., a sand lens within a clay deposit):

- **Detection**: Automatic via `A.contains(B)` — no user input needed for parent-child relationships
- **Slice generation**: When computing layer weights, both polygons intersect the slice strip. The inner polygon "overrides" the outer polygon in the overlap region. Compute: outer contribution = `A.intersection(strip).difference(B)`, inner contribution = `B.intersection(strip)`.
- **Base material**: The `contains(base_point)` check naturally returns the innermost polygon (check inner polygons first, ordered by area ascending).

### 8.2 Voids / Tunnels

A polygon with a hole (e.g., access tunnel in a dam):

- **Data model**: Shapely `Polygon` natively supports interior rings: `Polygon(exterior, [hole])`
- **LEM**: A void within the slice has zero weight and zero strength. The slice weight calculation would need to subtract the void area. This is conceptually simple but changes the assumption that slices are solid columns — would need careful handling of failure surfaces that pass through voids.
- **FEM**: Mesh generation can handle holes naturally (the hole is not meshed, and its boundary becomes a free surface).
- **Recommendation**: Defer void support to a future phase. The data model supports it (Shapely handles it), but the LEM implications need careful thought. For now, require all polygons to be solid (no interior rings).

### 8.3 Gaps Between Polygons

If polygons don't fully tile the cross-section, there are gaps with no assigned material. This is a modeling error — the user needs to fix their geometry. Detect during validation:

```python
domain = unary_union([p['polygon'] for p in polygons])
bbox = box(*domain.bounds)
gaps = bbox.difference(domain)
if gaps.area > tolerance:
    warn("Gaps detected between polygons — uncovered area has no material")
```


## 9. Migration Strategy

### 9.1 Backward Compatibility

- Existing input files with profile lines continue to work unchanged
- If both `profile` and `polygons` sheets contain data, raise an error (don't try to merge)
- `max_depth` is only used with profile lines; with polygons, the domain polygon defines all boundaries

### 9.2 Internal Representation

Convert profiles to polygons early, then use one code path:
- `build_polygons()` already exists and does this conversion
- All downstream code (slicing, search) works with polygons
- Profile-based slice generation becomes legacy and is removed once polygon path is verified


## 10. Implementation Order

1. **[done] `fileio.py`**: Unified convert-early representation. Geometry is always
   stored as `slope_data['polygons']` (list of `{'polygon': Polygon, 'mat_id': int}`,
   mat_id 0-based): the `polygons` sheet is used directly if present, otherwise
   profile lines are converted via `build_polygons()`. Profiles optional;
   both-specified raises an error.
2. **[done] `fileio.py`**: `build_ground_surface_from_polygons()` derives the
   ground surface (upper boundary of the polygon union) and the domain polygon
   (`slope_data['domain_polygon']`) from the unified polygons; `max_depth` becomes
   the domain's min elevation. Verified the union-derived ground surface is
   identical to the old `build_ground_surface()` for profile inputs.
3. **[done] `slice.py`**: Polygon-based slice generation. `generate_slices()`
   computes layer heights, weight, centre of gravity, and base material from each
   polygon's vertical extent at the slice centre (fast numpy interpolation via
   `_build_polygon_edges`), and takes slice-boundary breakpoints from polygon
   vertices. Equivalence-preserving (plan §12.3): produces bit-for-bit identical FS
   to the profile path on the LEM regression suite, and identical slices with
   `profile_lines` emptied (polygon-sheet input). The Ito & Matsui pile auto-H
   feature (`intersect_pile_with_materials`) was also migrated to polygons, so
   `generate_slices()` no longer reads `profile_lines` at all.
4. **[done] domain containment**: `generate_slices()` rejects any failure surface
   (circular or non-circular) that leaves the domain polygon, via a cached
   `prep(domain).covers(clipped_surface)` check. Flat-bottomed domains (all
   profile-line inputs, and polygon inputs with a horizontal base) take a fast path
   that skips the geometric test — the scalar `max_depth` clamp already bounds them
   — so existing analyses are unchanged (LEM regression identical). Only irregular
   bottoms (e.g. dipping bedrock) incur the containment check. `max_depth` is now
   the domain's minimum elevation for every input.
5. **[done] `mesh.py`**: `get_material_polygons(slope_data, reinf_lines)` returns
   mesh-ready polygons — the stored `slope_data['polygons']` for polygon inputs
   (skipping `build_polygons()`), or `build_polygons()` for profile inputs (with
   reinforcement integration). `main_seep`/`main_fem`/`main_mesh` use it, so FEM and
   seepage run on polygon-sheet inputs (verified end-to-end on the seep fixtures).
6. **[done] `plot.py`**: filled material-zone polygons via `plot_polygons_on_ax`,
   hatched domain base via `plot_domain_base`, and a shared `plot_base_geometry`
   used by `plot_inputs`/`plot_solution`/search-results plots; `compute_ylim` and
   `get_plot_elements_bounds` handle polygon inputs.
7. **Testing**: Verify polygon-based results match profile-based results for equivalent geometries
8. **[done] CAD import/export** (`xslope/cad.py`): see [`plan_io.md`](plan_io.md) §3.
   `export_dxf` writes a layered DXF; `dxf_to_polygons`/`import_dxf` read material
   zones (robust to LWPOLYLINE/POLYLINE, arc bulges, unclosed rings, loose LINE
   segments) and write the `polygons` sheet, seeding material names. Verified against
   all `poly_test/` fixtures and a full export→read round-trip (areas match exactly).
9. **Documentation**: Input template, sample problems, migration guide


## 11. Documentation Updates

- **Input template** (`docs/usage/input_template.md`): Add `polygons` sheet section with column definitions, examples, and guidance on when to use polygons vs. profiles
- **LEM overview** (`docs/lem/overview.md`): Add note on polygon geometry support; explain how ground surface and domain boundary are derived
- **Sample problems**: Add at least one polygon-based example (same slope as an existing profile-based example, to show equivalence)
- **FEM/Seepage**: Note that polygons can be used directly without profile-to-polygon conversion
- **CAD import**: New section documenting the DXF import workflow


## 12. Open Questions

1. **Performance of polygon-based slicing**: The intersection approach (Shapely `intersection()` per slice per polygon) may be slower than the current profile-based arithmetic. Need to benchmark. If too slow, consider caching polygon-strip intersections or using a spatial index (`STRtree`).

2. **Pore pressure with polygons**: Currently, pore pressure option (`'none'`, `'piezo'`, `'seep'`) is per-material. This doesn't change with polygons — the base material still controls pore pressure. But if the user has a seepage mesh, the mesh was generated from the same polygons, so integration should be seamless.

3. **Profile-to-polygon equivalence**: The existing `build_polygons()` makes certain assumptions (vertical side boundaries, `max_depth` as bottom). Need to verify that converting profiles to polygons and running through the polygon path gives identical results to the current profile path.
