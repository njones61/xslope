# XSLOPE Input Template

## Overview

The XSLOPE input template is the primary means of defining slope stability problems in xslope. It is an Excel workbook that contains all necessary information about the slope geometry, material properties, loading conditions, boundary conditions, and analysis parameters. The template is designed to support three main types of analysis:

- **Limit Equilibrium Method (LEM)**: Classical slope stability analysis using methods like Bishop, Spencer, Janbu, and others
- **Seepage Analysis**: Finite element groundwater flow analysis for steady-state and transient conditions
- **Finite Element Method (FEM)**: Stress-deformation analysis including the Shear Strength Reduction Method (SSRM)

The template uses a structured format with multiple worksheets (tabs), each dedicated to a specific aspect of the problem definition. This organization makes it easy to prepare complex slope stability analyses while maintaining clarity and avoiding errors.

Once analyses have been run, a project is more than this workbook: meshes and solutions
are written to sidecar files named after it. See
[Project Packaging](packaging.md) for what the full set is and how to move it as one
file.

## Download

A template for the Excel file can be downloaded here:

[input_template.xlsx](../inputs/input_template.xlsx)

The input file can be modified using any spreadsheet software such as Microsoft Excel, LibreOffice Calc, or Google Sheets.

## Loading the Template

The template is loaded into xslope using the `load_slope_data()` function from the `fileio` module:

```python
from xslope.fileio import load_slope_data

# Load the input template
slope_data = load_slope_data("path/to/your/input_template.xlsx")
```

This function reads all worksheets, validates the data, and returns a dictionary containing all parsed information. The resulting `slope_data` dictionary is then used by various xslope modules for analysis.

Loading checks that the file is well formed, not that it is complete — a
half-built model must always open. Whether a model carries everything the analysis
you are about to run actually needs is checked separately, when the run starts; see
[Preflight Input Checks](preflight.md).

## Template Structure

The template consists of 14 worksheets, each serving a specific purpose. Different worksheets are used by different analysis types: Limit Equilibrium Method (LEM), seepage analysis (SEEP), and Finite Element Method (FEM).

| Sheet Name | Description | LEM | SEEP | FEM |
|------------|-------------|:---:|:----:|:---:|
| [**main**](#worksheet-main) | Global parameters and instructions | X   | X    | X   |
| [**mat**](#worksheet-mat) | Material properties including strength, permeability, and stiffness | X   | X    | X   |
| [**profile**](#worksheet-profile) | XY coordinates of profile lines defining slope geometry | X   | X    | X   |
| [**polygon**](#worksheet-polygon) | Material zones defined as closed polygons (alternative to **profile**) | X   | X    | X   |
| [**piezo**](#worksheet-piezo) | Piezometric lines for pore pressure calculations | X   |      | X   |
| [**circles**](#worksheet-circles) | Circular failure surface definitions | X   |      |     |
| [**non-circ**](#worksheet-non-circ) | Non-circular failure surface coordinates | X   |      |     |
| [**dloads**](#worksheet-dloads) | Distributed surface loads | X   |      | X   |
| [**reinforce**](#worksheet-reinforce) | Soil reinforcement elements (anchors, nails, geosynthetics) | X   |      | X   |
| [**joints**](#worksheet-joints) | Slip surfaces with no reinforcement in them (rock joints, block contacts) |     |      | X   |
| [**piles**](#worksheet-piles) | Pile and concrete pier support elements | X   |      | X   |
| [**lloads**](#worksheet-lloads) | Line loads (concentrated forces on the ground surface) | X   |      |     |
| [**seep bc**](#worksheet-seep-bc) | Seepage analysis boundary conditions |     | X    |     |
| [**tseep**](#worksheet-tseep) | Transient seepage time series and run controls |     | X    |     |

The following sections describe each worksheet in detail, including the data structure and how it is used in analysis.

To check your geometry visually before running an analysis, open the file in [XSLOPE Studio](../studio/index.md) or
call `plot_inputs()` from Python — both draw the profile lines, material zones, piezometric lines, loads,
reinforcement, and boundary conditions to scale.

---

## Worksheet: main

![sheet_main.png](images/sheet_main.png)

The **main** worksheet provides global parameters that apply to all analyses and serves as the instruction page for the template. This tab contains:

- **Template version**: Tracks the template format for compatibility. The current version is **27**. xslope rejects a file whose version is newer than the installed version supports, so an older install cannot misread a newer template — and older files load unchanged, with anything a previous template lacked simply staying at its default.
- **Units** (`SI` or `Imperial`): declares the unit system for the model. Selecting a system fixes the unit weight of water to its standard value (**9.81 kN/m³** for SI, **62.4 pcf** for Imperial) and records the system with the model. XSLOPE is unit-agnostic and never converts your numbers — the declaration simply keeps the model's units explicit and self-consistent (SI = m, kPa, kN/m³; Imperial = ft, psf, pcf). If you leave this blank, xslope **infers** the system from the unit weight of water you enter (≈9.81 → SI, ≈62.4 → Imperial).
- **Time** (`sec`, `min`, `hr`, `day`, or `yr`): declares the time unit for every time-bearing quantity — hydraulic conductivity (length/time), specified flux, and the transient-seepage series and durations on the **tseep** sheet. Because xslope never converts, this one declared time unit governs them all together. Unlike the unit system, the time unit is **never inferred or guessed** (a wrong time label is worse than none), so it applies only when you set it here. Leave it blank for a static model with no time-bearing inputs; the **tseep** sheet requires it to be set.
- **Unit weight of water** (γw) — **[F/L³]**: used in pore pressure calculations. When you select a unit system, this cell is auto-filled with the canonical value, but you may **override** it — a value you type takes precedence (e.g. ≈10.05 kN/m³ or 64 pcf for seawater), and xslope warns at load time if your value differs from the canonical one by more than about 2%. With the Units selector blank, the value you enter here is what determines the inferred system.
- **Tension crack parameters**: Depth **[L]** and water-surface elevation **[L]** within tension cracks at the top of the failure surface
- **Seismic coefficient** (kh) — **[–]**: Horizontal seismic acceleration coefficient (dimensionless) for pseudo-static earthquake analysis

### Run options

The cells below the seismic coefficient record **how this model is meant to be
analyzed**, so the analysis travels with the file instead of living only in a dialog.
Every one of them is optional, and **blank always means unspecified** — the solver or
the Studio dialog uses its own default. A choice you make in a Studio dialog always takes
precedence over the value in the file.

- **LEM method** (`oms`, `janbu`, `bishop`, `corps`, `lowe`, `spencer`, `mprice`, or `all`) — **[–]**: the limit-equilibrium method this model is built for. `all` means "run every method" and is used by the batch drivers; the Studio run dialog, which solves one method at a time, leaves its own default in place for it. An unrecognized value is an error at load time, never a silent fallback to some other method.
- **Number of slices** — **[–]**: the slice count for limit-equilibrium runs. Minimum 2; blank uses the solver default.
- **K0 initial stress (FEM)** — **[–]**: at-rest lateral earth pressure coefficient for the FEM's **initial stress state**. Blank starts from zero stress and switches gravity on in one step, so the initial lateral stress is whatever elasticity produces — σ<sub>h</sub> = ν/(1−ν)·σ<sub>v</sub>, about 0.43·σ<sub>v</sub> at ν = 0.3, a number set by the stiffness rather than by the soil. Enter a value and the initial stress is built from the overburden instead, σ<sub>h</sub> = K0·σ<sub>v</sub> both in-plane and out-of-plane, and then iterated to equilibrium — in an SSRM run as a separate full-strength step, so that establishing the in-situ state is not charged to the strength reduction. Normally consolidated sand sits near Jaky's 1 − sin φ′; compacted fills and overconsolidated clays run at 1.0 and above. Set **1.0** to reproduce an RS2 result, whose models are authored with an isotropic K = 1 field stress. The effect is small on a cohesive slope and worth several percent on a reinforced or near-cohesionless one, always in the direction of a higher factor of safety. A value here must be positive; leave the cell blank for the gravity turn-on. See [K0 initial stress](../fem/overview.md#k0-initial-stress) for the formulation, how to choose a value, and the measured sensitivity.
- **Tension SRF (FEM)** (`YES` or `NO`) — **[–]**: whether the tensile-strength cap (`t_cut`) is reduced along with c and tan φ during a strength reduction. `YES` (what the template ships with) makes the factor of safety the factor on the whole strength envelope, shear and tensile. `NO` holds each cap at its authored value through the bisection. The setting changes nothing unless some material declares a cutoff **above zero**: a blank `t_cut` leaves the material its envelope's own apex, which strength reduction never moves, and a cutoff of `0` — a soil that carries no tension — is still `0` however it is divided. Studio's Run FEM dialog dims the matching checkbox on both kinds of model.
- **Mesh element type** (`tri3`, `tri6`, `quad4`, `quad8`, or `quad9`) — **[–]**: the element type the Build Mesh dialog opens on. Use quadratic elements (`tri6`, `quad8`, `quad9`) for every stress analysis; linear elements are for seepage only, because in plasticity they lock volumetrically and overestimate the factor of safety.
- **Mesh target size** — **[L]**: the target element size the Build Mesh dialog opens on. Setting it also turns auto-sizing off, because a size in the file is taken as the intended size.
- **1D element size** — **[L]**: the target element size along the model's 1D members — the pile and reinforcement lines — where it should differ from the mesh target size above. Blank means the members follow the mesh target size; a value here must be positive.
- **SSRM F min** / **SSRM F max** — **[–]**: the strength-reduction bracket the SSRM search starts from. F min must be less than F max.
- **Water loads** (`auto` or `manual`) — **[–]**: who supplies the weight of water standing on the slope. `auto` — what a new file carries — hands it to the engine, which derives the ponded-water load at solve time from the model's own water definition: the seepage boundary conditions wherever a seepage analysis is defined, otherwise the piezometric line. The two dloads sheets then carry **non-water** loads only — a surcharge, a footing, traffic — and the derived load is drawn on every plot in its own color, so that a load you did not enter is clearly visible. `manual` leaves the water load to you on the dloads sheets. In Studio the same selector sits in the **Global parameters** dialog. A blank cell means `auto` in a version-22-or-newer file and `manual` in every older one, and that is a correctness requirement rather than a preference: an older file already carries its water load typed in, and deriving a second one under it would count the reservoir twice. Use `manual` for a deliberate exception, and for a vendor-faithful transcription whose point is to reproduce another program's input exactly. See [Automatic water loads](preflight.md#automatic-water-loads).
- **Surface family** (`circular` or `non-circular`) — **[–]**: which failure-surface family a model that defines **both** a circular surface (the [circles](#worksheet-circles) sheet) and a non-circular one (the [non-circ](#worksheet-non-circ) sheet) actually means. Blank — what the template ships, and the normal state — means "whichever family the model defines", which is unambiguous for all but the rare file carrying both; there the circular surface is used and the non-circular one is ignored, and the [model checks](preflight.md) report this. Fill the cell (or answer the **Surface** selector in Studio's Run LEM dialog, which writes your choice here) and that answer travels with the file: the run, the plots and the next session all read the same one. A family named here is only honored when the model actually defines it, so a value left behind by a surface you later deleted can never select a surface that is not there.
- **Side BC** (`rollers` or `fixed`) — **[–]**: how the FEM restrains the left and right truncation boundaries of the model. `rollers` (the default, and what every file that leaves this blank gets) fixes horizontal movement and leaves the face free to settle vertically. `fixed` clamps both directions. Rollers are the usual choice: a truncation boundary is an artificial cut through ground that continues, and the ground on the other side of it does not hold the face up. `fixed` adds shear restraint that ground does not have, and stiffens a domain truncated close to the slope; use it to reproduce a program that fixes its side boundaries. The bottom boundary is fully fixed either way.

**Dimensional notation.** The field descriptions on this page tag each quantity with its dimensions in
brackets, built from three base dimensions — **[L]** length, **[F]** force, **[t]** time — plus the derived
forms the page uses: stress **[F/L²]**, unit weight **[F/L³]**, strength gradient **[F/L²/L]** (stress per
unit change in depth — dimensionally the same as [F/L³], written this way because it means something
different), hydraulic conductivity and flux **[L/t]**, specific storage **[1/L]**, moment **[F·L]**,
angle **[deg]**, and dimensionless **[–]**. XSLOPE never
converts, so a bracket means the same in either system; under a declaration it resolves to concrete units —
SI with **Time = day** gives [L] = m, [F/L²] = kPa, [F/L³] = kN/m³, [L/t] = m/day, while Imperial gives ft,
psf, pcf, ft/day.

These global parameters are accessed throughout the analysis. For example, the unit weight of water is used in 
computing pore pressures from piezometric lines, and the seismic coefficient is used to add horizontal inertial 
forces to each slice in limit equilibrium calculations. The tension crack parameters allow the simulation of a tension
crack at the top of the slope for cohesive soils to reduce the likelihood that negative normal forces develop along the 
face, which are unconservative in the limit equilibrium method. Filling the crack with water adds an extra level of 
conservatism as this applies a driving force to the failure surface.

The run options are read the same way everywhere: the Studio dialogs open on them, and a
library caller gets them on the loaded `slope_data` (`lem_method`, `num_slices`, `k0`,
`tension_srf`, `element_type`, `target_size`, `ssrm_f_min`, `ssrm_f_max`, `side_bc`) as
`None` wherever the file leaves a cell blank.

The declared **Units** and **Time** selectors are carried through to the display layer as well: they set the unit labels shown on generated plots and a
units-provenance line (for example `# units: SI, time: day`) written at the top of exported seepage and FEM result
files. When nothing is declared, those labels and the provenance line are omitted.

---

## Worksheet: mat

The **mat** worksheet has one row per material, starting in row 11. The number in
the **mat** column is the material ID that the **profile** and **polygon** sheets
refer to. Add rows below the last one for more materials.

The columns fall into three bands, each shown below. A color legend on the sheet
marks each column **LEM only**, **LEM & FEM**, or **FEM only**, and the cells a
material's options don't use are grayed out.

### Shear strength and stiffness

![sheet_mat1.png](images/sheet_mat1.png)

| Column | Units | What it is |
| --- | --- | --- |
| **γ** | F/L³ | Unit weight above the water table |
| **γsat** | F/L³ | Unit weight below the water table. Blank uses γ everywhere |
| **option** | | Strength model: `mc`, `cp`, `pow`, `hb`, or `elastic` (table below) |
| **c**, **φ** | F/L², deg | Cohesion and friction angle |
| **c/p**, **r-elev** | F/L²/L, L | Strength gain per unit depth, and the elevation it starts from (`cp`) |
| **d**, **psi** | F/L², deg | Kc = 1 envelope for [rapid drawdown](../lem/rapid.md) |
| **t_cut** | F/L² | Tensile strength cap (FEM only) |
| **E**, **nu** | F/L², – | Young's modulus and Poisson's ratio (FEM only) |
| **u** | | Pore pressure source: `piezo`, `seep`, `ru`, or `none` |
| **ru** | – | Pore pressure ratio, read when **u** = `ru` |
| **phi_b**, **s_cap** | deg, F/L² | Strength from suction above the water table, and its cap |
| **pow_a** … **pow_d** | | Power-curve envelope (`pow`) |
| **hb_sci** … **hb_d** | F/L², –, –, – | Hoek-Brown inputs (`hb`) |

**Unit weights.** Enter total unit weights, never buoyant ones. XSLOPE takes the
water off through the pore pressure, so a buoyant weight removes it twice. The
model has one water table, from the seepage solution if there is one and from the
piezometric line otherwise. Below it, a material with a γsat weighs γsat,
whatever its **u** option.

**Strength options.** Each option reads its own columns:

| option | Strength | Columns read |
| --- | --- | --- |
| `mc` | Mohr-Coulomb: τ = c + σ′ tan φ | c, φ |
| `cp` | Undrained strength increasing below an elevation: S<sub>u</sub> = c + c/p × (r-elev − y) | c, c/p, r-elev |
| `pow` | Power curve: τ = pow_a (σ′ + pow_d)<sup>pow_b</sup> + pow_c | pow_a, pow_b, pow_c, pow_d |
| `hb` | Generalized Hoek-Brown, for rock | hb_sci, hb_gsi, hb_mi, hb_d |
| `elastic` | Can't fail. Search surfaces may run along it but not through it | none |

`mc` with φ = 0 is a constant undrained strength. Enter **hb_sci** in the
model's stress units, not MPa (30 MPa is 30,000 kPa). The
[LEM overview](../lem/overview.md) explains each envelope, including
[Hoek-Brown](../lem/overview.md#hoek-brown-strength).

**Pore pressure (u).**

- `piezo`: from the line on the **piezo** sheet.
- `seep`: from a seepage solution; see
  [Seepage and slope stability](../seep/seep_slope.md).
- `ru`: u = ru × the vertical stress of the soil column.
- `none`: no pore pressure. Use it with undrained (total-stress) strengths.

Any other entry is rejected when the file loads. A model imported from RS2 gets
its **u** values from RS2's water settings; a pore-pressure grid or a groundwater
solution has no counterpart and arrives as zero pore pressure, flagged in the
import notes.

**Tensile cap (t_cut).** Only the FEM reads it; the LEM uses a
[tension crack](#worksheet-main) instead. `0` means no tension. Blank lets an `mc`
material carry the tension its envelope implies, c / tan φ, and leaves a `cp`
material with no limit at all, which can hold a steep crest together in a
strength reduction and raise the factor of safety. Set it for cohesive soils, and
always when matching an RS2 or Plaxis model. See
[Tensile strength in the SSRM](../fem/solver.md#tensile-strength-in-ssrm).

**Suction (phi_b, s_cap).** A positive **phi_b** adds strength s × tan φ<sup>b</sup>
from the suction s above the water table. It is read only for `mc`, `pow`, and
`hb` materials with **u** = `piezo` or `seep`. With `piezo`, the suction grows
without limit above the line, so **always set s_cap**. See
[Matric suction](../lem/overview.md#matric-suction-apparent-cohesion-above-the-water-table).

**Stiffness.** The FEM needs E > 0 and 0 ≤ ν < 0.5 for every material.

### Standard deviations

![sheet_mat2.png](images/sheet_mat2.png)

**σ(γ)**, **σ(c)**, **σ(φ)**, **σ(c/p)**, **σ(d)**, and **σ(ψ)** are the standard
deviations used by the reliability analysis, each in the units of its property.

### Seepage

![sheet_mat3.png](images/sheet_mat3.png)

| Column | Units | What it is |
| --- | --- | --- |
| **k1**, **k2** | L/t | Major and minor hydraulic conductivity. Equal values are isotropic |
| **alpha** | deg | Angle of the k1 direction from horizontal, usually 0 |
| **unsat** | | Unsaturated model: `lf` (default), `vg`, or `gard` |
| **kr0**, **h0** | –, L | Linear-front parameters (`lf`) |
| **a**, **n** | | van Genuchten α and n (`vg`), or the Gardner parameters (`gard`) |
| **l** | – | van Genuchten pore-connectivity exponent. Blank is 0.5 |
| **Ss** | 1/L | Specific storage (transient only) |
| **Sy** | – | Specific yield (transient, unconfined only) |

Leave **Ss** and **Sy** blank for steady-state flow. With a
[**tseep** sheet](#worksheet-tseep), every material needs **Ss**, and **Sy** is
needed when the model has an exit face. The
[seepage overview](../seep/overview.md) defines all of these, with typical
[van Genuchten values](../seep/overview.md#van-genuchten-model) by soil type.

---

## Worksheet: profile

![sheet_profile.png](images/sheet_profile.png)

The **profile** worksheet defines the slope geometry using XY coordinates of profile lines are the primary input for 
all types of analysis (LEM, SEEP, and FEM). 
Each profile line the top of a soil layer or profile and all of the soil below that line and above all of the lower 
profile lines is assumed to consist of the material associated with the profile line. The material id listed for each profile line references one of the materials in the 
material properties table. The material name in row 6 is found by using the material ids in row 5 to look up the 
name in the second column of the materials table.

To illustrate how profile lines can be used to define the geometry of a slope, consider the following slope 
with three 
layers:

![profile_lines.png](images/profile_lines.png)

The profile lines should always be drawn in the order of increasing depth, from top to bottom and the XY coordinates 
defining the line should be listed from left to right. All profile-line coordinates are lengths **[L]**. In the example above, the top profile line has three points, 
the next line has three points, and the last line has two points. Each line should have at least two points. The 
bottom of the slope is defined by the Max Depth parameter **[L]** (an elevation) at the top of the **profile** worksheet. This defines a 
horizontal base to the problem. During a limit equilibrium analysis using an automated search algorithm, the failure 
surface is not allowed to go below this depth. Thus, it can be thought of as a bedrock surface.

A profile-line model needs the bottom elevation in cell B2. A workbook that leaves it blank still opens,
and the model checks stop a run until it is entered.

Each profile line also has an optional **Size** cell — a target finite-element size **[L]** for
the material zone that line builds, used only when a mesh is generated. Leave it blank and the
zone meshes at the global target size; set it and elements inside that zone are driven down to
the value given, growing smoothly back to the global size outside it. Use it to resolve a thin
or critical layer without refining the whole model. It has no effect on a limit-equilibrium run,
which does not mesh.

The template includes tables for 15 profile lines, organized horizontally. However, you can copy additional tables to the right as needed. There is no limit to the number of profile lines that can be defined. Furthermore, each table includes 20 rows of XY coordinates, but you can add as many rows as needed.

During analysis, xslope uses these lines to:

1. Construct the ground surface by finding the highest elevation at each x-coordinate
2. Determine slice geometry when a failure surface is intersected with the profile
3. Assign material properties to slices based on which layer they fall within
4. Build polygons for finite element meshing in seepage or FEM analysis

---

## Worksheet: polygon

The **polygon** worksheet is an **alternative to the profile worksheet** for defining slope
geometry. Instead of describing each soil layer as a left-to-right line, each material zone is
described as a **closed polygon** — a self-contained region with an assigned material. Use whichever
method best fits the geometry; do **not** fill in both the **profile** and **polygon** worksheets in
the same file.

![sheet_polygon.png](images/sheet_polygon.png)

Polygons are well suited to geometries that are awkward to express as stacked profile lines, such as:

- Irregular or dipping bedrock surfaces
- Lenticular deposits (lens-shaped inclusions) embedded within another material
- Complex fill geometries such as zoned dam cross-sections
- Geometry imported from CAD software, which is typically drawn as closed regions

Each polygon is defined in its own table, organized horizontally just like the profile tables. For
each polygon you provide:

- A **Type**, which says what kind of region the polygon is. `material` (the default, and what a
  blank cell means) is a soil zone; the three `ssr` values are [SSR zone](#ssr-zones) analysis
  overlays; `refine` is a [mesh refinement region](#refine-regions); `joints` is a
  [joint region](#joint-regions). An unrecognized word is rejected at load time, never read as a
  material zone.
- A **Material ID**, which references one of the materials in the **mat** worksheet. As with the
  profile sheet, the name below it is filled in automatically. Both cells gray out for any Type
  other than `material`, because an overlay has no material.
- An optional **Size** — a target finite-element size **[L]** inside the polygon, used only when a
  mesh is generated. Blank means the global target size. Size is available on **every** polygon and
  layers on top of the Type rather than replacing it: a material zone may carry one, an SSR zone may
  carry one, and a `refine` polygon carries nothing else. See
  [Local mesh sizes](#refine-regions).
- The polygon **vertices** as XY coordinates **[L]**, one vertex per row.

A few rules govern how the vertices are interpreted:

- **Winding order does not matter** — vertices may be listed clockwise or counter-clockwise.
- **Polygons are closed automatically** — the last vertex is connected back to the first, so you do
  not need to repeat the starting point.
- **Nesting is detected from the geometry** — if one polygon lies entirely inside another (for example
  a sand lens within a clay deposit), the inner polygon overrides the outer one in the overlap region.
  No parent/child bookkeeping is required.
- Each polygon must have at least three vertices, and the polygons together should tile the cross
  section without overlaps or gaps.

The template includes tables for 15 polygons, organized horizontally, and you can copy additional
tables to the right as needed. Each table includes room for many vertices, and you can add as many
rows as required.

Unlike the profile worksheet, the **polygon** worksheet does **not** use a Max Depth parameter. The
overall extent of the model — including the ground surface and the bottom boundary — is derived from
the **union of all polygons** (the domain polygon). The ground surface is taken from the upper edge of
that union, and during an automated limit equilibrium search the failure surface is constrained to
stay within the domain polygon, which can therefore represent an irregular bedrock surface directly.

### SSR zones {#ssr-zones}

A polygon whose **Type** begins with `ssr` is not a material zone at all. It is an **SSR zone** — an
analysis overlay that tells the finite-element [strength reduction method](../fem/solver.md#ssr-exclusion-zones)
which part of the model to weaken.

| Type | Meaning |
|:-----|:--------|
| **`ssr reduce`** | **Search area.** Strength reduction applies **only inside** this polygon. |
| **`ssr hold`** | **Exclusion, full strength.** Inside is never reduced, but it can still yield. |
| **`ssr elastic`** | **Exclusion, elastic.** Inside is linear elastic and cannot yield at all. |

Any other word is rejected at load time rather than ignored, so a typo can never quietly run the
model unconstrained.

**Overlays.** An SSR zone is never meshed, never becomes a material region and
never generates slices. It changes nothing about the model except which elements the strength reduction
touches, and it is classified element by element, by where each element's center falls. A zone can
therefore be added to a finished model safely: the mesh and the factor of safety are untouched unless the
zone actually constrains something. Zones may overlap each other and may cross material boundaries
freely — the no-overlap rule above applies to material zones only. The limit-equilibrium solvers ignore
them entirely.

**How several zones combine.** There is one rule:

> The reduced region is the union of the `ssr reduce` zones, minus the union of the `ssr hold` and
> `ssr elastic` zones. With no `ssr reduce` zone drawn, it is the **whole model** minus those exclusions.

So exclusions always carve out — of a search area they sit inside, or of the model as a whole — and an
interior hole in a search area is drawn simply by putting an `ssr hold` (or `ssr elastic`) polygon on
top of it. An `ssr elastic` zone additionally makes its elements linear elastic, the same treatment as a
material whose strength option is `elastic`, but addressed by outline instead of by material.

**When to use one.** Use a search area when the mechanism you want to measure is confined to part of the
model and a competing one elsewhere would otherwise take over — a stiff foundation held at full strength
forces the failure up into the fill above it. This is RS2's "SSR Search Area" and "SSR Exclusion Area",
and both senses are drawn directly: mark the region you want reduced `ssr reduce`, or mark the region you
want held `ssr hold`. Both forms, and when each is the natural one to use, are covered under
[SSR search areas and exclusion zones](../fem/solver.md#ssr-exclusion-zones).

A zone drawn on the polygon sheet applies to every run of the model. A search-area polygon passed
explicitly to a run takes precedence over the file's zones (XSLOPE warns when both are present, rather
than quietly intersecting the two).

**Zones on a profile-sheet model.** Because zones are not geometry, the "do not fill in both sheets"
rule does not apply to them: a model whose geometry lives on the **profile** sheet may still carry SSR
zone rows on the **polygon** sheet.

### Local mesh sizes and refine regions {#refine-regions}

The **Size** cell on a polygon (and the matching cell on a [profile line](#worksheet-profile)) sets a
target finite-element size **[L]** inside that region. Blank means the global target size — the value
in the Build Mesh dialog, or **Mesh target size** on the [main](#worksheet-main) sheet. Set it and the
elements inside the region are driven down to the value given, growing smoothly back to the global size
outside it so there is no abrupt jump in element size at the boundary.

Size is **independent of Type**. A material zone with a Size is that soil zone, meshed finer. An SSR
zone with a Size is the same analysis overlay, with the elements it selects resolved more finely. And a
polygon whose Type is **`refine`** is nothing but a Size: it carries no material, never becomes a mesh
region, never generates slices and is invisible to every solver — the only thing it does is make the
mesh finer where it is drawn. A `refine` polygon must therefore carry a Size; one without does nothing,
and the [preflight check](preflight.md) reports this before a mesh-based analysis runs.

Use a refine region to resolve something the geometry does not mark out on its own — the ground under a
footing, the zone a slip surface is expected to pass through, or the tip of a cutoff wall. Use a
material or profile-line Size when the thing to resolve *is* a layer. Neither affects a
limit-equilibrium run, which does not mesh, and neither changes the mesh at all when left blank.

A Size only ever **refines**. A value at or above the global target size cannot make the mesh coarser
there; xslope warns at mesh time rather than leaving a setting with no effect to pass unnoticed.

Refine regions may overlap anything, including each other and material boundaries; where several apply,
the smallest size applies.

### Joint regions {#joint-regions}

A polygon whose **Type** is `joints` is an **authoring region**: it marks where a generated joint network exists.
Like a refine region it carries no material, is never meshed as a region, never generates slices and is invisible
to every solver — and unlike one it does not change the mesh either. The only thing that reads it is a
[generated joint set](#joint-sets) that names it as its region, and a region no set names has no effect.

Draw one when the ground a set exists in is not a material zone: one block of an outcrop, the rock behind a wall,
the part of a unit above a bench. A set that *is* confined to a material needs no polygon at all — it can name the
material directly.

The region's **name** is the block header over its column — the `Polygon #3` cell, which you can type over. A set
refers to it by that name (`poly:North block`) or by its number among the joint regions (`poly:1`). A header left
as the template wrote it is the block's number rather than a name, so such a region is referred to by number.

Two cells on a joint region do nothing, and the [model checks](preflight.md) report both rather than letting either
pass unnoticed: a **Material ID**, because the region is not a soil zone (if the polygon was meant to *be* that zone, its
Type is the cell to change), and a **Size**, because nothing about the region is meshed — draw a `refine` polygon
over the same ground to refine the mesh there.

---

## Worksheet: piezo

![sheet_piezo.png](images/sheet_piezo.png)

The **piezo** worksheet defines up to two water lines used to compute pore water pressures in limit equilibrium
analysis. A line is used by every material on the **mat** sheet whose pore pressure option is set to "piezo", so it
applies to soils analyzed with effective-stress (drained) strengths. Below the line, pore pressures are positive;
above it, they are zero.

Each line has a **Type** cell that states what the line represents, which determines how the pressure head is
computed:

- **piezo** (the default — a blank cell means the same): the line is a true *piezometric line*, i.e., the locus of
  pressure heads measured by piezometers at the slip surface. The pore pressure at a slice base is the full static
  head, u = γw · Δz, where Δz is the vertical distance from the slice base up to the line (upper figure below).

- **phreatic**: the line is the *phreatic surface* (the water table) in a slope with inclined, roughly
  surface-parallel steady seepage. In that flow regime the equipotentials are not vertical, so the head at depth is
  less than the static column and the static head is reduced by cos²θ, where θ is the inclination of the line
  segment directly above the slice base (lower figure below). This is the same correction as Slide2's "Hu: Auto"
  and the phreatic-surface option in XSTABL, and it matters most on steep phreatic lines — on a flat line the two
  types are identical.

![piezo.png](images/piezo.png)

Choosing between them: if the line came from piezometer readings, a flow net, or a seepage analysis, it is a
piezometric line — use **piezo**. If all you have is the water table position, **phreatic** is the appropriate
shortcut. For full rigor on complex problems, run a finite element seepage analysis and set the material pore
pressure option to "seep" instead; the pore pressures then come from the seepage solution directly and neither line
type is needed.

A line applies only over its own horizontal extent — its pressures are not extended past either end. That makes a
line that stops short of the section a legitimate model, for a water body that stops there too (a reservoir on one
side of a dam and no tailwater on the other). It does mean the line must reach everywhere it is *read*: if a slice
base, mesh node, or Gauss point in a material set to "piezo" falls beyond the end of the line, the analysis stops
with an error rather than treating that point as having no pore pressure. To model dry ground past the end of a
line, carry the line on at an elevation below the section. For the same reason, a material set to "piezo" in a file
whose piezo worksheet is empty is an error too; a material with no water takes **none**.

The worksheet provides space for two lines, which supports rapid drawdown analysis:

- **Piezometric Line 1 (columns A–B)**: steady-state or initial condition
- **Piezometric Line 2 (columns D–E)**: post-drawdown condition (only required for rapid drawdown)

Each line has its own Type cell and requires at least two XY coordinate pairs **[L]**, ordered from left to right. The table
is formatted for 20 rows, but coordinates can be entered beyond the bottom of the table as needed.

A file with no Type cell treats its lines as **piezo**.

---

## Worksheet: circles

![sheet_circles.png](images/sheet_circles.png)

The **circles** worksheet defines circular failure surfaces for limit equilibrium analysis. Circular surfaces are 
the most common assumption in slope stability analysis and are required for methods like Bishop's Simplified Method 
and Spencer's Method. XSLOPE supports up to 10 circular failure surfaces, each of which can be analyzed 
individually or used as starting points when searching for a critical failure surface with a minimum factor of 
safety using an automated search algorithm. A seepage or finite element model can leave this sheet empty: only a limit equilibrium run needs a surface, here or on the **non-circ** sheet.

![circles.png](images/circles.png){width=1000px}

Each row in the circles table specifies one circular failure surface with the following 
parameters:

- **Xo, Yo** — **[L]**: Coordinates of the circle center
- **Option**: Method for defining circle size - "Depth", "Radius", or "Intercept"
>>- **Depth**: Specify depth below ground surface at center location
>>- **Radius**: Directly specify circle radius
>>- **Intercept**: Specify a point (Xi, Yi) where circle should pass through
- **Depth, R, Xi, Yi** — **[L]**: Associated values depending on option selected

During a limit equilibrium (LEM) analysis, XSLOPE performs the following steps:

1. Constructs the circular arc geometry from the parameters
2. Finds intersection points with the ground surface
3. Divides the arc into slices
4. Assigns material properties to each slice based on its position

A common problem in limit equilibrium analysis is finding the critical 
failure surface. The automated search algorithm sometimes converges to a location corresponding to a local minimum 
of the factor of safety, but this location may not correspond to the global minimum. To find the critical surface, 
it is common practice to start the search at multiple locations and then analyze the results to identify the 
critical surface. XSLOPE automates this process by testing each of the defined circle locations when performing an 
automated search and then continuing to iterate from the location with the lowest factor of safety. When defining multiple circles, a good strategy is to start define one circle passing through the toe of the slope 
(for steep slopes) and one circle at the base of each soil layer.

### Search window (optional)

The **Search window** block at the right of the sheet
confines the automated circular search to a chosen region of the geometry. This is the
counterpart to the "search limits" other programs offer, and its purpose is to make the
search settle on a particular **local** minimum rather than the global one — a benched
slope, for example, has several competing mechanisms and you may want the factor of
safety of one of them specifically.

Every limit is optional and independent. Leave a cell blank and that limit is not
applied; leave all ten blank and the search is exactly the unconstrained search
described above.

| Cell | Meaning |
|---|---|
| **Entry X min / max** — **[L]** | X range the failure surface's crest-side (higher-ground) endpoint must fall in. |
| **Exit X min / max** — **[L]** | X range the toe-side (lower-ground) endpoint must fall in. |
| **Center box X min / X max / Y min / Y max** — **[L]** | Rectangle the circle centers are confined to. The refining grid stays inside it, so the search cannot walk out. |
| **Max tangent depth** — **[L]** | The lowest **elevation** the circle's bottom (its tangent point) may reach. |
| **Min slip depth** — **[L]** | Minimum depth below the ground surface a surface must reach. Rejects shallow surficial "skin" mechanisms, whose factor of safety is depth-independent on a cohesionless face and would otherwise win. |

The same ten cells are edited in XSLOPE Studio by the **Search window** group under the
circles table (see [Editing inputs](../studio/editing.md)), which draws the entry and exit
ranges and the center box on its preview.

A range is applied only when **both** of its ends are filled; XSLOPE does not supply a
missing end. Entry and exit ranges and the tangent-depth limit **reject** a trial surface
that violates them rather than clamping it, so the reported minimum lies within the window; the center box confines the grid
itself. Both ends of every range must be increasing (min ≤ max), which is checked when the
file loads.

---

## Worksheet: non-circ

![sheet_noncirc.png](images/sheet_noncirc.png)

The **non-circ** worksheet allows definition of arbitrary non-circular failure surfaces. Some slopes include thin 
layers with especially weak soils. In such cases, a failure surface where much of the surface is confined to the 
weak layer can be more critical than a circular failure surface. A non-circular failure surface is defined by a set 
of XY points, listed from left to right. The table is formatted for 20 rows, but extra rows can be added below the table if needed. Generally, the leftmost point is the entry point and the rightmost point is 
the exit point and both should correspond to the ground suface. 

![noncircular.png](images/noncircular.png)

The worksheet contains three columns:

- **X, Y** — **[L]**: Coordinates of points along the failure surface, defined sequentially from left to right
- **Movement**: Direction of movement constraint at each point (used in automated search) — a selector, not a numeric value

Non-circular surfaces can only be used with LEM methods that support non-circular failure surfaces. For the methods 
supported in XSLOPE, the following table defines which methods support non-circular failure surfaces:

| Method                                         |         Non-circular         |
|------------------------------------------------|:----------------------------:|
| Ordinary Method of Slices                      |              No              |
| Bishop's Simplified Method                     |              No              |
| Janbu Method                                   |             Yes              |
| Force Equilibrium - Corps of Engineer's Method |             Yes              |
| Force Equilibrium - Lowe & Karafiath Method    |             Yes              |
| Spencer's Method                               |             Yes              |

For the movement option, the following values are supported:

- **Free**: No movement constraint
- **Horiz**: Point moves in the horizontal direction only
- **Fixed**: Fixed movement constraint at each point

For the **Free** option, if the point is the first or last point in the list, the movement constraint is applied 
such that the point moves left or right along the ground surface. For interior points, the movement constraint is 
applied such that the point moves in the direction of tangent of the failure surface at that point. This process is 
described in more detail in the [Automated Search Algorithms](../lem/search.md) section.

---

## Worksheet: dloads

![sheet_dloads.png](images/sheet_dloads.png)

The **dloads** and **dloads (2)** worksheets define distributed surface loads applied to the slope. These represent surcharge loads 
such as traffic, buildings, stockpiled materials, or other surface loading. They are also used with submerged slopes 
to represent the force of the water on the slope. During limit equilibrium analysis, distributed loads are applied to the top of each slice, which affects either or both of the driving and resisting forces depending on the slope angle and load orientation. The **dloads** sheet defines loads used in a normal slope stability analysis or the first stage of a rapid drawdown analysis. The **dloads (2)** sheet defines loads used in the second stage of a rapid drawdown analysis, and has the same layout:

![sheet_dloads2.png](images/sheet_dloads2.png)

Each worksheet is formatted for 4 distributed loads, but additional loads can be added by copying and pasting more tables to the right. Each table is formatted for up to 20 rows, but additional rows can be added below the end of table if necessary.

Each distributed load is defined by a **Direction** and a series of points with:

- **X, Y** — **[L]**: Coordinates of points along the load distribution line, ordered from left to right
- **Normal** — **[F/L²]**: Load intensity (force per unit area) along the line

At least two points are required to define each load block.

**Direction** (`normal` or `vertical`) sets which way the load pushes. `normal` — the default, and
what a blank cell means — applies it **perpendicular to the loaded line**. That is right for
anything that acts as a pressure on a surface, and it is what water does, so ponded water and
reservoir loads always use it. `vertical` applies the same intensity **straight down**, which is what
a dead-weight surcharge does: a stockpile, a fill, or an equipment load presses down under gravity
whatever the ground beneath it is doing.

The distinction matters only on an **inclined** loaded surface, where it can be large. On ground
falling at angle β the perpendicular reading carries a horizontal component of tan β times the load —
on a 6° crest that is 11% of the surcharge, pushed sideways into the hill, which is a real restraining
force the surcharge does not actually apply. Both readings carry the same total force; they differ
only in where it points. On level ground the two are identical. The direction is set per load block,
so a model may mix them — a pond on the face and a stockpile on the crest — and both sheets carry the
option independently, so the two stages of a rapid drawdown may differ.

The direction applies to limit-equilibrium slices and to the finite-element edge tractions alike. The load distribution typically follows the ground 
surface. The points should be listed in order from left to right. For example, consider the following load distribution:

![dist_loads.png](images/dist_loads.png)

**Point order.** Left to right is the convention, but a line entered right to left is not an error and does not
change the answer: a load line that runs monotonically in decreasing X is re-oriented to increasing X when the file
is read, and a `normal` load's push direction is taken from the **geometry** — perpendicular to the line, into the
material — rather than from the order the points were typed in. A line whose X genuinely turns back on itself (an
overhang, or a small backstep at a toe) has no increasing-X form and is kept exactly as authored.

The forces are defined using a unit width perpendicular to the plane of the slope. For the water load, 
the force at each point would be the unit weight of water multiplied by the height of water above the point in 
question. For the example shown above, the distributed load would consist of three points, the first two points 
having the normal force and the last point having a normal force of zero. The surcharge load on the right would be 
defined by two points with a normal force for each.

---

## Worksheet: reinforce

![sheet_reinforce.png](images/sheet_reinforce.png)

The **reinforce** worksheet defines soil reinforcement elements such as soil nails, ground anchors, geosynthetic 
reinforcement, or tiebacks. These elements provide additional resistance to sliding by mobilizing tensile force
where the failure surface crosses them. Each reinforcement object is represented as a straight line defined by the XY coordinates of 
the endpoints. Each line also has a set of properties that define the type of support, the strength of the
reinforcement, and the anchorage at each end. The [Reinforcement Overview](../reinforcement/overview.md#the-columns)
describes every column, and [Reinforcement Types](../reinforcement/types.md) gives the entries for each kind of
support, with typical values.

![reinforce.png](images/reinforce.png)

The template is formatted for up to 30 reinforcement lines (rows 3-32), but additional rows can be added to the
table as needed.

---

## Worksheet: joints

![sheet_joints.png](images/sheet_joints.png)

The **joints** worksheet defines slip surfaces that carry no reinforcement: a rock joint, a bedding plane, the
contact between a facing block and the one beneath it, the back of a retaining wall against the soil behind it.
Each is a straight line between two points, and the finite element mesh is **split** along it, so the material on
one side is a separate body from the material on the other and the two can slide along it and part across it. The
reinforce sheet's `Joint` column does the same thing to a line that *does* carry a sheet; this worksheet is for the
lines that do not.

The template is formatted for up to 30 joint lines (rows 3-32), but additional rows can be added as needed. Every
column on this sheet is **FEM only** — the limit equilibrium engines do not read it.

Each joint line is defined by:

- **Geometry**:<br>
>>Label: Name used in error messages, summaries, and plots (optional). A [generated
set](#joint-sets) names its rows `bed-01`, `bed-02`, … after the set they came from<br>
>>x1, y1 **[L]**: Start point coordinates<br>
>>x2, y2 **[L]**: End point coordinates<br>
- **Strength**:<br>
>>c **[F/L²]**: Joint cohesion. Blank is 0, which with a stated `phi` is the purely frictional contact a dry joint
usually is.<br>
>>phi **[degrees]**: Joint friction angle. **Required** — `c` alone defaults to zero, so a blank `phi` would be a
frictionless surface the two sides slide along under any load at all.<br>
>>c_res **[F/L²]**, phi_res **[degrees]**: Residual cohesion and friction angle — what the joint keeps once it has
slipped. A rough surface shears through its asperities the first time it reaches its limit and does not rebuild
them, so the drop is instantaneous and permanent: the pair stays on the residual branch for the rest of the run
even where it later closes or unloads. Blank is the peak value, which is a joint with no residual branch at all.
Neither may exceed its peak.<br>
>>dil **[degrees]**: Dilation angle. A rough joint rides up on its asperities as it slides, so a unit of slip opens
it by $\tan$`dil`. Where the material around it holds it closed that opening builds normal traction against `kn`,
and with it shear strength; where the sliding block is free to lift, it lifts. Blank is 0, a joint that slides
without opening.<br>
>>t_cut **[F/L²]**: Tension cutoff. The joint **opens** when the normal traction passes it — carrying neither shear
nor normal traction — and closes again when the two faces come back into contact. Blank is 0, which is a contact
that holds no tension at all.<br>
- **Stiffness**:<br>
>>kn **[F/L³]**: Normal stiffness of the interface. Leave blank to derive it from the softer adjacent soil,
$E_{adj}$ over a virtual thickness of 0.1 × the 1D element size.<br>
>>ks **[F/L³]**: Shear stiffness of the interface. Leave blank to derive it the same way, from $G_{adj}$.<br>
>>Jred: **No** holds the joint at full strength through a strength reduction — for a construction detail rather
than the design question. Blank or **Yes** reduces `c` and $\tan$`phi` by the trial factor along with the soil's,
which is the usual choice.<br>

Joint lines may **meet** — at a T, at a crossing, at a corner, end to end — and where they do the split copies the
shared node once per wedge of material around it, so each element keeps the material on its own side of every line
through the point. That is what lets a facing column with a joint on its back face, a joint under its base and a
joint between every course of blocks slide and tip the way a stack of blocks does, instead of behaving as one solid.
What two joint lines may not do is lie **on** one another along part of their length, or run along the outer
boundary of the section, where there is material on one side only and nothing for the other face of the joint to be;
preflight reports both as errors, naming the lines.

The interface element, its constitutive law, the derived stiffnesses and what the results show are the same as for a
jointed reinforcement line: see
[Joints without reinforcement](../reinforcement/fem.md#joints-without-reinforcement).

### Generated sets {#joint-sets}

A jointed rock mass is rarely described one line at a time, so a whole **set** of them can be generated at once —
from Studio's [Build network](../studio/editing.md#build-network) dialog, or from `xslope.joints` in a script. A
generated line is an ordinary row on this sheet: it can be read, edited by hand, saved and reloaded like a typed
one, and nothing downstream can tell the two apart.

What marks it is its **Label**, which is the set's name and the row's place in the set — `bed-01`, `bed-02`, … —
and that name is what the plots, the messages and the report show. The description the rows were generated from
is not stored: the joint lines are the input. To change a network, remove it (**Remove set** in Studio's joints
editor removes every row of the set the selected row belongs to) and build another. A label of any other
shape — `base joint` — is a name you chose, and it belongs to no set.

---

## Worksheet: piles

![sheet_piles.png](images/sheet_piles.png)

The **piles** worksheet defines piles, piers and walls that resist a slide in shear and bending. Each is a straight
line between its two end points, vertical or battered. The end points can be entered in either order, and the higher
end is the head. The [Piles and Walls Overview](../piles/overview.md#the-columns) describes every column, and
[Pile and Wall Types](../piles/types.md) gives the entries for each kind of pile and wall, with typical values.

![Pile Example](images/pile_example.png){width=800px}

The template is formatted for up to 20 piles, but additional rows can be added to the table as needed.

---

## Worksheet: lloads

![sheet_lloads.png](images/sheet_lloads.png)

The **lloads** worksheet defines line loads — concentrated forces applied at a point on the ground surface, per
unit width of slope. A typical use is the weight of a soil nail wall's shotcrete facing in a limit equilibrium
model, applied on the crest just behind the top of the face; in a finite element model, draw the facing as a
polygon instead
([A soil nail and a shotcrete facing](../piles/types.md#a-soil-nail-and-a-shotcrete-facing)). (For loads spread
over an area of the surface, use the **dloads** worksheet instead.)

The template is formatted for up to 20 line loads (rows 3-22). Each line load is defined by:

- **Label**: Name used in error messages, summaries, and plots (optional)
- **x, y** — **[L]**: Coordinates of the point where the load acts. The point must lie on (or within a small tolerance of)
  the ground surface.
- **P** — **[F/L]** (per unit width of slope): Force magnitude
- **Angle** — **[deg]**: Direction of the force measured from horizontal in degrees. Leave blank for the default of **−90**
  (straight down — a weight).

During limit equilibrium analysis, the load is applied to the slice whose top boundary contains the point, entering
the equilibrium equations as force components with a real moment arm (analogous to the pile force terms, but applied
at the top of the slice rather than at the failure surface). The finite element analysis applies it as a force at
the nearest mesh node ([Loads](../fem/overview.md#distributed-loads)).

---

## Worksheet: seep bc

![sheet_seepbc.png](images/sheet_seepbc.png)

The **seep bc** and **seep bc (2)** worksheets define boundary conditions for finite element seepage analysis. 
Boundary conditions 
specify where water enters or exits the domain and the magnitude of hydraulic head on the boundary. There are three 
types of boundary conditions: specified head, specified flux, and exit face. **Specified head** boundaries correspond to free water on 
the face of the slope and the magnitude of the head is the height of water above the datum defined for the problem. 
**Exit face** boundaries conditions are used for unconfined problems are applied to the "downstream" side of the slope 
where water exits the slope. In the unconfined seepage solution, the phreatic surface intersects the exit face at 
some point (exit point) that is determined as part of an iterative solution process. For points on the exit face below 
the exit point, a head = elevation (zero pressure) condition is applied. For points on the exit face above the exit 
point, the head is determined by the pore pressure equation. 

Each of the five **Head/Flux BC** tables carries a **type** cell (a dropdown reading
`head`, `reservoir`, or `flux`) above its value cell. The type determines how the value
directly to its right is interpreted:

| type | value is | units |
|---|---|---|
| `head` | specified head — a **plain Dirichlet** total head held at every node of the polyline, at all times (may be a suction/negative-pressure head) | **[L]** |
| `reservoir` | a reservoir/tailwater **level** — a **submerged-only** Dirichlet: nodes at or below the level are held at it, nodes above it become seepage (exit) faces | **[L]** |
| `flux` | specified flux — the **normal Darcy velocity**, **positive into the domain** | **[L/t]** |

Both `head` and `reservoir` are total heads (height of water above the datum) and, for a
level drawn at or below the applied value, behave identically; they differ only in what
happens to a node the water level leaves above the line (a `head` holds it; a `reservoir`
releases it to an exit face). Use `reservoir` for a free-water face that rises or falls
(dam pool, tailwater, pond) and `head` for a drained face or an imposed head that must be
enforced regardless of elevation.

A **specified flux** (Neumann) boundary prescribes the *rate* at which water crosses the
boundary instead of the head on it. The usual use is rainfall infiltration or recharge applied
to the ground surface, where the water table position is a result of the analysis and so
cannot be entered as a head. Two points to be careful about:

- The flux is a flow per unit *area* of boundary, not a total discharge spread over the
  boundary — a rate of 1×10⁻⁶ m/s applied along a 40 m surface admits 4×10⁻⁵ m³/s per meter of
  slope, not 1×10⁻⁶.
- The sign matters and is easy to get backwards. **Positive is inflow** (infiltration,
  recharge). A negative value withdraws water.

Unlike a head boundary, whose points are simply pinned, a flux boundary's points define a
**polyline whose edges carry the load**, so it needs at least two points and they must run
along the boundary of the mesh. Zero flux is the default condition on any boundary that is not
otherwise specified, so a flux boundary only needs to be defined where the flux is non-zero. A
model with no specified-head boundary and no exit face anywhere is singular — head would be
determined only up to an additive constant — and xslope does not solve it.

**Time-varying boundaries (transient).** A head or flux **value** cell may hold the *name* of a
time series defined on the [**tseep** sheet](#worksheet-tseep) instead of a number. When it does,
that boundary becomes time-dependent — its head (for a `head` type) or its flux (for a `flux`
type) follows the named series through the run. A series is simply a curve of numbers versus time;
whether those numbers mean head or flux is decided by the block's own **type** cell, so one series
can drive several boundaries. The name must match a series header on the tseep sheet exactly —
xslope reports an error listing the available series names if it does not — and a name is only
valid when a tseep sheet is present. A constant boundary takes a plain number. A series on a `head` block is a plain Dirichlet held at $h(t)$ at every node of the
polyline at all times; a series on a `reservoir` block is the rising/falling level. For a
series-driven **reservoir** boundary, draw the polyline over the *full* face the water can ever
cover: at each time step only the currently submerged nodes are held at the series level, while
nodes above the waterline automatically become seepage-face (exit face) nodes so the
still-saturated soil behind a freshly exposed face can drain through it — see
[head types](../seep/transient.md#head-types-head-and-reservoir).

For a typical unconfined problem, there is one upstream specified head boundary condition and a single downstream 
exit face. For confined problems, there is typically one upstream and one downstream specified head boundary 
condition. Additional boundary conditions can be defined to represent the water level in an 
excavation, infiltration on the ground surface, etc. The sheet is formatted for one exit face and up to 5 head or flux boundaries. Additional 
boundary 
conditions can be added by copying and pasting more tables to the right. Each table is formatted for up to 20 rows, 
but additional rows can be added below the end of table if necessary.

For most analyses, only the main **seep bc** sheet is used. However, the **seep bc (2)** sheet is used for rapid drawdown 
analysis where a second seepage solution is used to calculate the pore pressures corresponding to the drawdown 
condition. It has the same layout as the main **seep bc** sheet:

![sheet_seepbc2.png](images/sheet_seepbc2.png)

During seepage analysis, xslope:

1. Builds a finite element mesh from the **profile** geometry
2. Applies specified head values at nodes on the specified head boundaries
3. Integrates any specified flux over the boundary edges it covers, into consistent nodal inflows
4. Applies exit face conditions where water exits
5. Solves the Laplace equation (∇·(k∇h) = 0) for hydraulic head at all nodes
6. Computes pore pressures (u = γw(h - y)) and flow velocities

For coupled seepage-LEM analysis, the computed pore pressures can be exported and later interpolated onto slice bases using the "seep" option in the **mat** worksheet.

See the [Seepage Analysis](../seep/overview.md) section for more details and the [samples](../seep/samples.md) 
section for sample input files illustrating various boundary conditions

---

## Worksheet: tseep

![sheet_tseep.png](images/sheet_tseep.png)

The **tseep** worksheet defines a **transient** (time-dependent) seepage analysis: a set of time
series that drive time-varying boundary conditions, together with the controls for the transient
run. It is used only for seepage. **A filled-in tseep sheet is what makes a seepage analysis
transient** — leave the sheet empty (as it is in the blank template) and seepage is the
steady-state analysis described in the [seep bc](#worksheet-seep-bc) section.

This section documents the tseep sheet as an *input*. How the transient analysis is formulated and
solved — the storage term, the time-stepping scheme, the initial and boundary conditions, and the
saved output frames — is described in [Transient Seepage](../seep/transient.md).

**Time-series table** (left side of the sheet). The first column, headed **time** **[t]**, lists times in
increasing order, in the model's declared [time unit](#worksheet-main). Each remaining column is a
**named series**: the column header is the series name (the template ships with the default names
`t1`…`t5`, but any short name works — in the image above two have been renamed `reservoir` and
`rain`), and the cells beneath hold that series' values at the listed times. A
[seep bc](#worksheet-seep-bc) value cell references a series by this name, and the referencing
block's **type** decides whether the numbers are read as heads **[L]** or as fluxes **[L/t]** (matching
the head/flux BC convention above) — so one series can drive several boundaries.

Each series is interpolated **linearly** between its values, and held constant before its first
time and after its last. A **blank cell** inside a series means "no breakpoint here" — the series
interpolates straight through it — so each series needs values only at its own breakpoints, and
independent series with different breakpoints can share the single time column (the `rain` series
in the image is blank at times where only `reservoir` changes). A **step change** is entered by
repeating a time on two consecutive rows with different values.

**Run controls** (right side of the sheet):

- **duration** — **[t]**: the total time to simulate. It may extend past the last time in the table; series
  values are held constant beyond their last breakpoint.
- **save_interval** — **[t]**: the spacing between saved output frames. Leave blank for an automatic default
  (roughly 50 frames over the duration).
- **save_times** — **[t]**: an optional list of extra individual times, entered down the column, at which a
  frame is additionally saved — useful for capturing specific instants, such as the times at which
  a published verification result reports its values.
- **stage_1**, **stage_2** — **[t]**: optional times that couple a transient run to a
  [rapid-drawdown](../lem/rapid.md#transient-solution) stability analysis — the two instants (for example the steady
  full-reservoir state at `stage_1` and the drawn-down state at `stage_2`) whose pore pressures
  supply the two stages of the drawdown calculation. Leave both blank when not doing rapid
  drawdown; if you set one you must set the other, with `stage_1` earlier than `stage_2`.
- **stability_time** — **[t]**: the optional instant that a single-time
  [LEM or FEM stability analysis](../seep/transient.md#stability-time) with `u = seep` reads its
  pore pressures from. Like the stage times it selects a frame out of the transient solution rather than
  changing it, and the stepper is forced to land on it. **Leave it blank and a stability
  run reads the last saved frame** — usually the drained end state. A run started from Studio can
  name a different instant for that run alone; storing it here makes a scripted or headless
  re-run read the same frame.

A tseep sheet requires the **Time** unit on the [main](#worksheet-main) sheet to be set, and every
material to carry a specific storage `Ss` (with `Sy` as well on unconfined models) — see
[Worksheet: mat](#worksheet-mat). Any save time, stage time, or series breakpoint that falls beyond
the run **duration** draws a load-time warning, since it would never be reached — and so does a
**stability_time** beyond it.

---

## Notes

- All variables must use one consistent unit system throughout the template. Declare it with the
  **Units** selector on the [main](#worksheet-main) sheet — **SI** (m, kPa, kN/m³) or **Imperial**
  (ft, psf, pcf). XSLOPE never converts your numbers; the declaration fixes the unit weight of
  water and keeps the model's units explicit. Time-bearing quantities (permeability, flux,
  transient series) additionally follow the single **Time** unit declared there.
- Angles are specified in degrees **[deg]**
- The template is designed to be flexible - you need not fill in all worksheets for every analysis
- Always check your geometry visually after entering data — open the file in XSLOPE Studio or call `plot_inputs()`
- Templates can be saved and reused for parametric studies or similar projects
