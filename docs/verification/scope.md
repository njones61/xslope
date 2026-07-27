# Vendor Feature Coverage

Every vendor program XSLOPE is verified against — Rocscience Slide2 and RS2, and
Seequent GeoStudio — carries features XSLOPE does not implement. This page is the
list of them, with the reason for each. It exists so that the scope of the package
is stated rather than inferred, and so that a *blocked* or *not supported* row on a
corpus page can point at a named ruling instead of re-arguing it.

Everything listed here is a deliberate position, not a backlog. Features XSLOPE
**does** support are documented in the [Input Template](../usage/input_template.md)
reference and in the analysis sections; where a supported capability exists but a
particular benchmark file does not yet carry a vendor's value for it, that is recorded
row by row on the corpus pages, not here.

---

## What was surveyed {#surveyed}

The list is derived from a field-by-field sweep of every vendor model file in the
verification archives. Each file was parsed and every input token — key/value records
in the Rocscience formats, XML elements and attributes in a GeoStudio `.gsz` —
enumerated mechanically, then clustered into engineering features and checked against
what the XSLOPE code actually does.

| Archive | Files | What it is |
|---|---:|---|
| RS2 (native) | 156 | RS2 slope-stability verification models |
| RS2 (Slide2 import) | 108 | the same archive's "(Slide2 Import)" models — exact geometry for 52 Slide2 problems |
| RS2 groundwater | 34 | RS2 groundwater verification models |
| Slide2 | 69 | Slide2 tutorial and hand-calculation models |
| SLOPE/W verification | 43 | the GeoStudio SLOPE/W verification corpus |
| SLOPE/W examples | 53 | the GeoStudio 2025.2 SLOPE/W example library |
| SEEP/W examples | 32 | the GeoStudio 2025.2 SEEP/W example library |
| **Total** | **495** | plus the Slide2, RS2 and GeoStudio verification manuals, read for feature enumeration |

That sweep produced 103 engineering features. Roughly two-thirds are supported and
carried by the input template; the remainder are on this page.

## Carried versus selected {#carried-versus-selected}

**A usage count below is the number of files that carry a feature *non-trivially* — not
the number of files whose schema mentions it.** The distinction matters more than it
sounds. RS2 writes its entire material-property dictionary into every model file
regardless of which constitutive model the file actually selects, so cap models,
hardening parameters, thermal properties and unsaturated curves appear in all 298 RS2
files while being selected by none. GeoStudio likewise writes SIGMA/W and QUAKE/W
fields into every `.gsz`. About four fifths of the raw fields in the RS2 archives are
that kind of boilerplate.

Those fields are still listed here rather than dropped, because a feature that appears
in every file and is selected in none looks exactly like a coverage gap until someone
checks. Rows of that kind are marked *out of scope* and say so in the note.

## Status vocabulary {#status}

Corpus pages use a [shared status vocabulary](rocscience.md); in that vocabulary every
row on this page is *not supported* — a deliberate exclusion rather than an unbuilt
item of work. Within that, the **Ruling** column separates two kinds:

- **gap** — a capability XSLOPE does not have. It could be built; nothing in the
  benchmark corpus currently demands it, so it is documented instead. Where a benchmark
  row *is* blocked by one, the row names it.
- **out of scope** — a standing exclusion: a discipline XSLOPE is not (dynamic
  analysis, coupled consolidation, soil-cover climate modelling), or vendor schema
  boilerplate that no model selects.

Each table gives the vendor feature and what it is, how many vendor files select it,
which benchmarks it affects, and what to do instead in XSLOPE where an alternative
exists.

!!! note "Documented, not built"
    Some vendor features look like gaps and turn out to have zero demonstrated need. The
    **vertical seismic coefficient** is the exemplar: all three vendor formats carry it,
    all three XSLOPE importers detect and report it, and every one of the 495 files in
    the archives sets it to zero. Prefabricated **wick drains** are the same story — the
    field is present in all 298 RS2 files and selected in none. Neither is built, and
    both are listed here so that finding the field in a vendor file does not read as an
    oversight.

---

## Geometry and meshing {#geometry}

| Vendor feature | Vendor usage | Benchmarks affected | XSLOPE position | Ruling |
|---|---|---|---|---|
| **Mapped / structured meshing regions** — quadrilateral mesh patterns imposed on a region rather than free triangulation | RS2 groundwater 8; most GeoStudio files | none | XSLOPE meshes unstructured only. This affects discretization quality, never the input physics; the seepage and FEM benchmarks are mesh-converged at their tagged sizes. | gap |
| **Axisymmetric, plan-view and 1-D analysis modes** — solving a radial or one-dimensional idealization instead of a plane-strain section | 2 RS2 groundwater; 14 GeoStudio | none currently blocked | Slope stability is inherently plane-strain, so the gap is seepage-only: the seepage solver accepts planar problems, and there is no axisymmetric element formulation. Confined problems in polar coordinates can still be posed as plan-view planar models — the [confined radial flow anchor](seep.md#verification-confined-radial) is exactly that. | gap |
| **Borehole / RSLog soil-profile import** — building layer polygons from borehole logs | 63 Slide2 | none | A data-preparation convenience that produces the same polygons the `polygon` and `profile` sheets already accept. Not a physics input. | gap |
| **Weak-layer entity** — a thin seam declared as a special object so the search engine treats it preferentially | 2 Slide2 | none | Model a weak seam as an ordinary thin polygon. The vendor entity only steers the search, and search steering is already handled by the search-window arguments. | gap |

## Strength models {#strength}

| Vendor feature | Vendor usage | Benchmarks affected | XSLOPE position | Ruling |
|---|---|---|---|---|
| **Residual / brittle post-peak soil strength** — a second strength pair the material drops to after peak is passed | 110 RS2-native models | every RS2 SSRM benchmark is solved at peak strength only; no row is locked against a residual-strength answer | XSLOPE has no strain softening for soil. Reinforcement lines do soften (a residual tensile capacity, `Tres`), but the continuum does not. The RS2 importer reports the dropped residual values. This is the largest single gap by vendor usage, and closing it is a solver project, not a template field. | gap |
| **Anisotropic strength** — cohesion and friction that vary with the orientation of the slice base relative to a bedding or discontinuity dip | 2 RS2-native, 5 RS2-import, 7 Slide2, 6 GeoStudio | [Slide2 VP105](rocscience.md) and [GeoStudio §2.47](geostudio.md) are blocked on it | Strength in XSLOPE is isotropic: one strength option per material, applied the same way at every base angle. All three importers flag an anisotropic material as *not the same strength* rather than silently importing the horizontal values. | gap |
| **Shear–normal (discrete) strength function** — a tabulated τ against σ_n envelope | 7 Slide2, 3 GeoStudio | GeoStudio "Material Model Shear Normal" | For a smoothly curved envelope, use the power-curve option (`option='pow'`), which fits the same shape analytically. A general tabulated envelope would need a new material option in both the slice and FEM strength paths. | gap |
| **SHANSEP** — undrained strength normalized by effective overburden and OCR | 3 Slide2, 3 GeoStudio | Slide2 Tutorial 31; GeoStudio "SHANSEP" | The nearest supported analog is the c/p option (`option='cp'`), which grows undrained strength linearly with depth below a datum. SHANSEP proper needs stress history and effective overburden tracked per slice. | gap |
| **Compound strength** — a base material overlaid by a second envelope, whichever governs | 2 GeoStudio | GeoStudio "Material Model Compound Strength" | One strength option per material is mutually exclusive by design. Combining two envelopes on one material is a solver change. | gap |
| **Drucker–Prager** | 17 RS2-native models | none locked against a Drucker–Prager answer | The RS2 importer maps it to Mohr–Coulomb and says so. Adding it to the viscoplastic core is a solver project no benchmark currently demands. | gap |
| **Jointed-rock strength models** — Barton–Bandis and joint-network rock-mass models | 10 RS2-native, 5 RS2-import, 2 GeoStudio | RS2 §32 family; GeoStudio "Rock Joint and Rock Mass Shear Strength Models" | Rock masses are modelled with the generalized Hoek–Brown option (`option='hb'`), which is a rock-mass envelope, not a discrete-joint model. The RS2 importer reports jointed materials as unmapped. Distinct from the interface *elements* listed under [reinforcement and structures](#reinforcement). | gap |
| **Unit weight as a function of depth or position** | 2 GeoStudio | GeoStudio "Unit Weight Definition", "Spatial Variation of Soil Properties" | Each material carries one moist and one saturated unit weight. For a step change with depth, split the layer into two polygons with different materials. | gap |
| **Cap, hardening and critical-state models** — Cam-Clay, hardening soil, non-linear elastic | schema present in every RS2 and most GeoStudio files; **selected in none** | none | These are the boilerplate case described [above](#carried-versus-selected): RS2 writes the whole property dictionary for every material, GeoStudio writes SIGMA/W fields into every project. Listed so that finding them in a file is not mistaken for a gap. | out of scope |
| **Anisotropic (transversely isotropic) elasticity** | schema present in most files; **selected in none** — every RS2 file declares linear elastic | none | Boilerplate, as above. | out of scope |
| **Dilation angle ψ** — a non-zero plastic potential angle in the FEM | 8 RS2-native models | RS2 §62 family | FEM dilation is fixed at zero (non-associated flow), the conservative strength-reduction convention. Note that the `psi` column on the `mat` sheet is the Duncan & Wright rapid-drawdown parameter, not a dilation angle — the two must not be conflated. | out of scope |

## Pore pressure and water {#pore-pressure}

| Vendor feature | Vendor usage | Benchmarks affected | XSLOPE position | Ruling |
|---|---|---|---|---|
| **Several independent piezometric surfaces assigned per material** | 8 GeoStudio | GeoStudio "Probabilistic – Syncrude Dyke", "Staged rapid drawdown – CE Example", "Pore-Water Pressure Definition for a Levee Stability" | The `piezo` sheet carries two lines, and the second is reserved as the drawn-down pairing of the first for rapid drawdown — it is not a free second surface. Where per-material water genuinely differs, run FE seepage and use `u='seep'`, which gives a full field rather than a set of lines. The GeoStudio importer reports that it took the first surface. | gap |
| **B-bar / excess pore pressure from loading** — Δu = B̄·Δσ generated by construction | 2 Slide2, 5 GeoStudio | Slide2 Tutorial 12; GeoStudio "Excess Pore-Water Pressure using Skempton's B-bar", "Submerged Slope with Excess Pore Water Pressure" | Requires a construction history to compute Δσ against, and XSLOPE has no [staged construction](#analysis). Transient seepage covers dissipation of an *imposed* excess pressure, which is how the consolidation benchmarks are built. | gap |
| **Air-entry / critical-suction threshold** — the vendor form of the suction cap, expressed as an air-entry value | 74 RS2 files, 11 SLOPE/W verification models | RS2 §38 family; GeoStudio rows using a capped suction | Suction credited to strength is capped by `s_cap`, a ceiling on the suction *stress* — not the vendor's air-entry threshold on the suction itself. The two agree until suction exceeds the air-entry value; both importers say so explicitly. | gap |
| **Pore-pressure grid / spatial pore-pressure function** — pore pressure supplied as an interpolated grid or a 3-D function | 23 Slide2, 7 GeoStudio | Slide2 Tutorial 05; GeoStudio "Lanester Embankment", "Pore-Water Pressures Defined Using a Spatial Function" | Use FE seepage: `u='seep'` with a companion nodal pore-pressure file is the equivalent capability and is strictly more general, since the field satisfies the flow equation. Where a published grid is a *measured* construction-induced field with no flow solution behind it, no seepage run can regenerate it and the corpus rows say so. The Slide2 importer flags the grid. | gap |
| **Pore-air pressure** — a non-zero u_a in the unsaturated effective-stress terms | 2 GeoStudio | GeoStudio "Effect of Pore-Air Pressure on Stability" | The matric-suction strength model assumes u_a = 0 explicitly, so the suction term reduces to (−u_w)·tan φ_b. | gap |

## Seepage {#seepage}

| Vendor feature | Vendor usage | Benchmarks affected | XSLOPE position | Ruling |
|---|---|---|---|---|
| **Fredlund–Xing and user-tabulated relative-conductivity curves** | a minority of the unsaturated files; the archives otherwise use van Genuchten or Gardner forms | none blocked on the curve shape alone | Three relative-conductivity models are supported: linear-front (`lf`), van Genuchten (`vg`) and Gardner (`gard`). Fredlund–Xing and free-form tabulated curves are not. Note the Gardner option is the power form; an *exponential* Gardner law is a different curve, and the one groundwater row that needs it is blocked with that reason recorded. | gap |
| **Volumetric water-content function (SWCC)** | 21 RS2 groundwater and 45 GeoStudio models define an explicit curve | the transient groundwater rows | Standing ruling: the SWCC is diagnostic only. Transient storage uses specific storage and specific yield plus the capacity term implied by the relative-conductivity model, not a tabulated θ(ψ). Where a vendor's steeper SWCC makes its wetting front lag XSLOPE's, the affected benchmark rows carry that caveat. | out of scope |
| **Climate and free-drainage boundary conditions** — land–climate interaction (evaporation, snow, solar radiation), unit-gradient outlets, head-versus-volume | 10 GeoStudio; RS2's climate function blocks are boilerplate present in all 298 files and selected in none | [SEEP/W "Infiltration into multi-layered system"](geostudio.md#seepw-t06); GeoStudio "Evaporation from the Wilson Soil Column", "Soil Cover Modeling" | The boundary-condition set is total head, pressure head, reservoir, specified flux, potential seepage face, and time-varying head or flux. Atmospheric coupling is a soil-cover discipline outside slope stability. A unit-gradient outlet has no equivalent — an exit face clamps the boundary to zero pressure rather than letting it drain freely. | out of scope |
| **Discharge section** — flow reported across an arbitrary internal query line | 4 RS2 groundwater, 1 Slide2 | groundwater rows whose published answer is a sectional flow rate | The seepage solver reports the total boundary flow rate, which is what the built groundwater rows lock against. An internal section is a post-processing query rather than a model input; if a row is ever locked on one, this becomes a reporting need. | gap |
| **Thermal and freeze–thaw, gas phase, solute transport, density-dependent flow** | 5 SEEP/W examples; RS2's thermal blocks are boilerplate in all 298 files | none | Outside the scope of the package. | out of scope |

## Loads {#loads}

| Vendor feature | Vendor usage | Benchmarks affected | XSLOPE position | Ruling |
|---|---|---|---|---|
| **Distributed load applied at an angle** — a surface pressure not perpendicular to the surface it acts on | 5 RS2-native, 4 RS2-import, 6 Slide2 | RS2 rows whose angled loads the importer declines to map | A distributed load is perpendicular to its own drawn line by construction; there is no angle field. An inclined surface load can be entered as an equivalent line load on the `lloads` sheet, which does take a direction. The RS2 importer reports angled loads as not imported. | gap |
| **Vertical seismic coefficient** | **zero** — every file in all 495 carries k_v = 0 | none | Pseudo-static loading uses a horizontal coefficient only. All three importers detect a non-zero vertical coefficient and report it rather than dropping it silently. See [Documented, not built](#status). | gap |
| **Seismic time history, Newmark displacement, dynamic analysis** | 11 RS2 groundwater, 1 Slide2, 1 GeoStudio | GeoStudio "Newmark Deformation Analysis" | Pseudo-static analysis is the whole seismic story in XSLOPE. Dynamic response and permanent-displacement analysis are separate disciplines. | out of scope |
| **Prescribed (non-zero) displacements as a boundary condition** | 2 RS2-native, 5 RS2-import | RS2 §24 and §32 families | FEM boundary conditions are generated from the mesh geometry — fixed base, roller sides — and the boundary-condition types are free, fixed, roller and applied force. There is no prescribed-displacement type. | gap |

## Reinforcement and structures {#reinforcement}

| Vendor feature | Vendor usage | Benchmarks affected | XSLOPE position | Ruling |
|---|---|---|---|---|
| **Capacity reduction factors** — creep, installation damage, durability and chemical factors applied to a nominal geosynthetic strength | 4 Slide2, 5 GeoStudio | Slide2 Tutorial 37; GeoStudio "Reinforcement with Anchors" / "…with Geosynthetics" | The `reinforce` sheet takes the **already-derated** long-term design strength as `Tmax`. The vendor stores the derivation chain; multiply the factors through before entering. This is arithmetic rather than physics, but it is a transcription trap worth naming. | gap |
| **Interface-dependent pullout** — bond resistance growing with normal stress through an adhesion and interface friction angle | 4 Slide2, 5 GeoStudio | GeoStudio "Snailz – Reinforced Slope", "Reinforcement with Geosynthetics" | The FEM has an opt-in Coulomb bond-slip model for reinforcement lines. The LEM applies a constant pullout rate over the anchorage length instead, which the GeoStudio importer flags as an approximation of the vendor's stress-dependent form. | gap |
| **Reinforcement force distribution and shear option** — force applied perpendicular to the surface, or distributed along the bar rather than concentrated | 5 GeoStudio | GeoStudio "Reinforcement with Anchors", "Snailz – Reinforced Slope" | Reinforcement direction is `tangent` (along the slip surface) or `axial` (along the bar). There is no perpendicular option and no distributed-along-the-bar application. The GeoStudio importer flags both. | gap |
| **Reinforcement sets and patterns** — one object that expands into N parallel lines by count, spacing and length | 2 Slide2, 3 GeoStudio | GeoStudio "Pockoski and Duncan – Tie-Back Wall", "Maccaferri 37 m High MSE Wall" | Enter the expanded lines. The pattern object is a drawing generator, not a physics input; the GeoStudio importer names it in its unrecognized-entity list. | gap |
| **Liners, beams and shotcrete** — surface structural elements | 10 RS2-native, 5 RS2-import | RS2 §32 family | The pile element (a six-degree-of-freedom beam) is the nearest analog, but it is a discrete pile, not a surface liner. The RS2 importer reports liners as not imported. | gap |
| **Joint and interface elements** — zero-thickness elements with normal and shear stiffness and their own strength | 10 RS2-native, 5 RS2-import | RS2 §32 family | There is no interface element in the FEM: adjacent materials share nodes and are perfectly bonded. This is distinct from reinforcement bond-slip, which models load transfer between a bar and the surrounding soil, not slip along a continuum interface. | gap |
| **Prefabricated vertical (wick) drains** | **zero** — field present in all 298 RS2 files, selected in none | none | See [Documented, not built](#status). | out of scope |

## Analysis machinery {#analysis}

| Vendor feature | Vendor usage | Benchmarks affected | XSLOPE position | Ruling |
|---|---|---|---|---|
| **Sarma method / non-vertical slices** | 2 Slide2 | the nine problems of the Slide2 Sarma verification manual | Slices are vertical in every method XSLOPE implements. Sarma's inclined-slice formulation is a method-sized addition; the importer reports it rather than substituting another method. | gap |
| **Morgenstern–Price side-force functions beyond constant and half-sine** — trapezoidal and user-defined shapes | 17 Slide2, 2 GeoStudio | GeoStudio "Chen and Shao Frictionless Slope", "Gravity retaining wall" | `mprice()` accepts `f_type='constant'` (which reproduces Spencer) and `f_type='half_sine'` (the default); anything else raises. In practice the factor of safety is weakly sensitive to the shape, and the corpus rows lock against the half-sine. | gap |
| **Metaheuristic search optimizers** — cuckoo search, particle swarm, simulated annealing, multi-modal optimization | most Slide2 files declare them as defaults; 9 GeoStudio select one | GeoStudio "Search Techniques Block Specified"; Slide2 Tutorial 38 | Circular searches use adaptive grid refinement and non-circular searches use coordinate descent. This is a difference in critical-surface search *quality*, not in physics — a different optimizer against the same objective. Individual modes of a multi-modal problem are reachable by constraining the entry, exit and tangent-depth windows. | gap |
| **Tension crack specified as a line or an angle** | 4 GeoStudio (10 more use the supported constant-depth form) | GeoStudio rows using the Angle option | The `main` sheet carries a constant crack depth and the fraction of it filled with water. An arbitrary crack polyline, or a crack angle from which the depth is derived, is not supported: the GeoStudio importer takes the deepest point of a crack surface and declines the Angle option outright. | gap |
| **K₀ / in-situ field stress** — an initial stress state specified independently of the gravity solve | every RS2 file; 3 GeoStudio | every RS2 strength-reduction benchmark, in principle | XSLOPE has no K₀ input. Initial stresses emerge from the elastic gravity solve, which gives K₀ = ν/(1−ν). RS2's gravity turn-on instead produces an isotropic K₀ = 1 initial state. The difference is real, but it does not survive into a fully plastic critical strength-reduction factor: the affected rows lock within about 1%. | gap |
| **Staged construction** — multiple stages with element birth and death, staged loads and staged water | 22 RS2-import, 14 RS2 groundwater, 3 Slide2, 23 GeoStudio | Slide2 Tutorial 31; GeoStudio "Staged rapid drawdown" family, "Multi-stage Pseudostatic Analysis" | Model each stage as its own input file. The only staging XSLOPE has is the three-stage Duncan & Wright rapid-drawdown procedure, which is limited to the LEM, and an optional dry-then-flooded gravity sequence within a single FEM solve. Generic staged construction is a combined solver and template project. | gap |
| **Consolidation, coupled flow–deformation, creep** | 11 RS2-native, 64 RS2-import, all 34 RS2 groundwater, 1 GeoStudio | none — the one-dimensional consolidation rows are reachable as uncoupled transient seepage | Coupled consolidation is outside scope. Note that the viscoplastic "creep" pseudo-time inside the FEM is a stress-return algorithm, not physical creep. | out of scope |
| **Springs and user-specified restraints** | 10 RS2 groundwater | none locked | FEM restraints are generated from the mesh geometry; there is no user-defined boundary-condition entity and no spring element. | gap |

## Probabilistic analysis and design standards {#probabilistic}

| Vendor feature | Vendor usage | Benchmarks affected | XSLOPE position | Ruling |
|---|---|---|---|---|
| **Correlated random variables** — a correlation matrix among the sampled inputs | 1 Slide2 tutorial populates it (the section is a default in all 69 files); 3 GeoStudio | GeoStudio "Basic Probabilistic Stability Analysis", "Cannon Dam #2"; Slide2 Tutorial 27 | Monte Carlo draws each parameter independently; there is no covariance structure. The one coupling that is enforced is physical: saturated and moist unit weight move together. | gap |
| **Point-estimate (PEM) sampling** | offered by RS2; the archives' probabilistic models use Monte Carlo | none | Reliability is available three ways — Taylor series, Monte Carlo, and a strength-reduction (FEM) variant — none of them point-estimate. The Taylor-series method serves the same purpose at similar cost. | gap |
| **Spatial variability / random fields** — a correlation length and covariance function applied to a material property field | 2 Slide2, 2 GeoStudio | Slide2 Tutorials 33/34; GeoStudio "Spatial Variation of Soil Properties" | Standing scope ruling: uncertainty is treated parameter-by-parameter, not as a spatially correlated field. | out of scope |
| **Design standards and partial factors** — Eurocode 7 design approaches, favourable/unfavourable factor sets | 1 Slide2, 7 GeoStudio | Slide2 Tutorial 21; GeoStudio "Eurocode design Case 1/2/3" | There is no partial-factor layer, and the reported factor of safety is unfactored — the Slide2 importer warns when a design standard was selected. The alternative works and is verified: bake the factors into the authored material and load values and run an ordinary analysis. Both Eurocode rows in the SLOPE/W corpus are built that way and reproduce GeoStudio's Overdesign Factor to within 0.1% — see [§2.45](geostudio.md#gs-2-45) and [§2.46](geostudio.md#gs-2-46). | gap |

## Workflow and bookkeeping {#workflow}

| Vendor feature | Vendor usage | Benchmarks affected | XSLOPE position | Ruling |
|---|---|---|---|---|
| **Several analyses or scenarios in one file** | every multi-scenario Slide2 project and multi-analysis `.gsz` | none — the importers ask which one to take | One input workbook is one model. The Slide2 and GeoStudio importers expose the scenario and analysis lists and import the one you name, which also forces the choice to be deliberate: GeoStudio assigns materials *per analysis*, so two analyses over one geometry can be different soils. | gap |
| **Parent/child analysis chaining** — one analysis taking its initial or pore-pressure conditions from another | 23 GeoStudio | every SEEP/W-parented SLOPE/W row | The equivalent is the seepage companion file: `u='seep'` reads a solved nodal pore-pressure field. XSLOPE carries the field, not the chain that produced it — and the GeoStudio importer reads a solved parent field directly. | gap |
| **Display and solver-bookkeeping fields** — colours, hatches, fonts, page layout, output routing, version stamps | every file in every archive | none | Not model inputs. The one exception is deliberate: the GeoStudio importer reads material colours so that imported models plot with the same palette they had in GeoStudio. | out of scope |

---

## Related reading {#related}

- [Verification and Validation](index.md) — how the benchmark corpus is organized.
- [Rocscience Slide2 Corpus](rocscience.md) — the shared status vocabulary, and the
  problem-by-problem record.
- [Input Template](../usage/input_template.md) — what a shipped model file carries.
- [GeoStudio Import/Export](../usage/geostudio.md) — what an import brings across and
  what it reports as unmapped.
