# GeoStudio (SLOPE/W) Import/Export

XSLOPE reads and writes GeoStudio SLOPE/W project files (`.gsz`), so a model built in
GeoStudio can be brought into XSLOPE without redrawing it, and an XSLOPE model can be
handed to a GeoStudio user.

A `.gsz` is a ZIP holding the model as a single XML document. Unlike a DXF — which is
just lines on layers, and needs a person to state what each layer *means* — a `.gsz`
identifies its regions, materials and water conditions, so the import needs no mapping
wizard. The only choice is **which analysis** to import.

## Importing a GeoStudio model

=== "Studio"

    **File → Import GeoStudio (SLOPE/W)…**, choose the `.gsz`, and pick an analysis.
    The model lands on the canvas unsaved, with a list of anything that didn't come
    across. Review it, fill any gaps, then **Save As**.

=== "Script"

    ```python
    from xslope.geostudio import read_gsz, list_analyses, gsz_to_slope_data

    gsz = read_gsz("Cannon Dam.gsz")
    for a in list_analyses(gsz):
        print(a["id"], a["name"], a["method"])

    slope_data, caveats = gsz_to_slope_data(gsz, analysis_id=1)
    for c in caveats:
        print("note:", c)
    ```

    Or go straight to an XSLOPE input file:

    ```python
    from xslope.geostudio import import_gsz
    from xslope.fileio import default_template_path

    caveats = import_gsz("Cannon Dam.gsz", default_template_path(), "cannon.xlsx")
    ```

### Choosing the analysis {#why-the-analysis-matters}

A `.gsz` usually holds several analyses over one geometry, and they differ in more than
the slip surface. GeoStudio assigns **materials per analysis**, not per region — so the
same slope can be one soil in the first analysis and a different soil in the second. The
choice changes the model, so XSLOPE asks for it rather than guessing.

### What comes across

| GeoStudio | XSLOPE |
|---|---|
| Regions | Material zones (`polygons`) |
| Region→material assignment (per analysis) | `mat_id` |
| Mohr-Coulomb materials (c′, φ′) | `option='mc'` with c′, φ′ |
| **Undrained materials** (φ = 0, strength entered as cohesion) | `option='mc'` with c = Su, φ = 0 |
| **SHANSEP** with a constant ratio (its older name is S=f(overburden)) | `option='ratio'` with c/p = the ratio and c = 0, see below |
| **S=f(depth)** (undrained strength rising with depth below the top of its layer) | `option='cp'` measured from that top, see below |
| **Spatial Mohr-Coulomb** with a cohesion that is a straight line in elevation, on φ = 0 | `option='cp'`, see below |
| **Bedrock** (impenetrable) | `option='elastic'` — an impenetrable material, see below |
| Material color | The zone color in Studio |
| Piezometric surface, and the materials assigned to use it | Piezo line, on the materials that use it; Type `phreatic` only where the analysis applies SLOPE/W's phreatic correction. A line that stops short of the model's edges is continued with no water beyond its ends, see below |
| **Ru**, per material | `u = ru`, with the same value |
| **Water above the ground surface** | **A distributed load** — see below |
| **A parent SEEP/W analysis's pore-pressure field** | The `seep` option, on SEEP/W's own mesh — see below |
| **A staged rapid drawdown**: two piezometric surfaces and the drawdown strengths | XSLOPE's rapid drawdown: piezometric lines 1 and 2, and d and ψ — see below |
| **Surcharges** (a wedge of fill above the ground — the surcharge pressure is its *unit weight*) | `dloads`, varying with the depth of fill |
| **Reinforcement** — nails, geosynthetics, anchors, tiebacks | `reinforcement_lines` — see below |
| **Line loads** (magnitude + trend/plunge) | `line_loads` |
| **Tension crack** (+ percent filled with water) | `tcrack_depth`, `tcrack_water` |
| **Unit weight above the water table**, where the material switches it on | `gamma` above the water table and `gamma_sat` below — see below |
| Unit weight of water | `gamma_water` |
| Horizontal seismic coefficient | `k_seismic` |
| SLOPE/W's critical circle, if the file was solved and it fits | A single trial circle |

The ground surface and domain are derived from the zones, as for any polygon-based model.

### Materials

Mohr-Coulomb materials carry straight over as `option='mc'` with their c′ and φ′.
GeoStudio's undrained (φ = 0) materials store the undrained strength in the cohesion
field, and import as `option='mc'` with c = S<sub>u</sub> and φ = 0.

**Strength from the vertical effective stress.** SLOPE/W's SHANSEP model with a
constant ratio gives S<sub>u</sub> = (τ/σ′<sub>v</sub>) σ′<sub>v</sub>. It imports as the
[`ratio` option](../lem/overview.md#undrained-strength-ratio), with c/p set to the
ratio. SHANSEP has no minimum strength, so c is set to 0. A SHANSEP material often still
carries a cohesion and friction angle from an earlier model; SLOPE/W does not use them,
and neither does the import. Both programs take σ′<sub>v</sub> as the weight of the soil
column above the slice base, plus surcharges and ponded water, minus the pore pressure.
Neither counts point loads. On SLOPE/W's own long-term SHANSEP example, the two
strengths agree within 0.02 kPa at every column of SLOPE/W's critical surface.

Two other SHANSEP forms have no exact match:

- **A ratio of S × OCR<sup>m</sup>, with OCR varying with elevation.** XSLOPE has one
  ratio per material, so the import uses the lowest, at the smallest OCR. Its strength
  is then at or below SLOPE/W's, and the notes give the largest difference.
- **A strength read from a curve of σ′<sub>v</sub>.** No single ratio is that curve, so
  the material is imported with no strength, and the notes ask you to set one.

SLOPE/W often pairs SHANSEP with B-bar, to hold σ′<sub>v</sub> at its value before a
fill is placed. XSLOPE has no B-bar (see below), so it counts the fill's weight in
σ′<sub>v</sub>, and its strength there is higher than SLOPE/W's. The notes say so.

**Strength rising with depth.** SLOPE/W's S=f(depth) gives
S<sub>u</sub> = S<sub>top</sub> + r d, where d is the depth below the top of the
material's own layer, measured at each slice. It imports as the
[`cp` option](../lem/overview.md), with c = S<sub>top</sub>, c/p = r, and r_elev at the
top of the layer. Where that top is level, the strength is the same as SLOPE/W's.
Where it is not, r_elev is its highest point, so XSLOPE's strength is never above
SLOPE/W's, and the notes give the largest difference. S=f(depth) can also cap the
strength at a maximum, which XSLOPE does not have. The import names any field on these
materials that it does not read.

**Strength that varies in space.** SLOPE/W's Spatial Mohr-Coulomb model lets the
cohesion, friction angle and unit weight each follow a function of x, a function of y,
or a function over x and y. One case is exactly XSLOPE's `cp` option: a cohesion that is
a straight line in elevation, rising with depth, on a friction angle of zero. It imports
as `cp`, with c the cohesion at the line's upper point, r_elev that point's elevation and
c/p the rate of increase. SLOPE/W holds the strength at its end values beyond the line's
two ends, and so does `cp` above r_elev. Below the line's lower end `cp` keeps rising
while SLOPE/W's strength stays constant, and the notes say so when the material reaches
that deep. On SLOPE/W's own example (a foundation whose cohesion rises from 5.4 kPa at
y = 9 m to 27.52 kPa at y = 0), the imported model matches SLOPE/W's factor of safety to
within 0.2% on its critical circle.

XSLOPE has one value per material for every other case. A cohesion or friction angle
function imports as its lowest value, so XSLOPE's strength is at or below SLOPE/W's, and
a unit weight function imports as the average of its points. The notes name the
material, the function, its range and the value used.

SLOPE/W's S=f(datum), and every other strength model XSLOPE does not have, are imported
with whatever cohesion and friction angle they carry, and the notes say that this is not
the same strength.

SLOPE/W can also scale any strength by a function of the slice base's inclination (an
anisotropy function). XSLOPE has no such scaling, so the import uses the unscaled
strength and says so.

**Unit weight above the water table.** A SLOPE/W material's unit weight is its
saturated weight, and it applies everywhere unless the material switches on a
separate weight above the water table:

| SLOPE/W | XSLOPE |
|---|---|
| No switch (a weight above the water table may be stored, and is not used) | γ = the unit weight, throughout |
| A constant weight above the water table | γ = that weight, γ<sub>sat</sub> = the unit weight |
| A weight above the water table computed from the volumetric water content | γ = the unit weight, throughout, and the notes say so |

XSLOPE splits the two weights at its water table, the piezometric line or the top of a
seepage field's positive pressure, where SLOPE/W splits them where the pore pressure turns
negative; the two are the same water table unless it is perched. On SLOPE/W's own Unit
Weight Definition example (16 above, 20 below), the split puts the weight of the critical
slip surface's sliding mass within 0.1% of SLOPE/W's. XSLOPE has no weight that follows
the water content, so in the third case its slices above the water table are heavier than
SLOPE/W's.

GeoStudio's **Bedrock** is impenetrable — slip surfaces cannot enter it. XSLOPE's
**`elastic`** option is the exact analog: an impenetrable zone that cannot fail and that
a slip surface may ride along but never cut into. So a bedrock region is imported as an
elastic material, keeping its unit weight; its c and φ are dropped because it cannot
fail. That reproduces SLOPE/W exactly, and the material round-trips back to Bedrock on
export. Importing bedrock as ordinary soil would let trial surfaces cut straight
through it.

Matric suction is reported rather than converted. SLOPE/W parameterizes unsaturated
(negative pore-pressure) shear strength as a **suction ceiling on the piezometric
surface** — a maximum suction, and a cap on how much of it counts — not as a
per-material φᵇ. That is a threshold on the pore pressure, semantically different from
XSLOPE's Fredlund `phi_b` / `s_cap` apparent-cohesion pair, so the import reports it and
leaves the choice to you: set `phi_b` (and `s_cap`) by hand if you need unsaturated
strength credited.

### Water and pore pressure

A piezometric surface's points belong to the **analysis that owns it**, not to the
shared geometry table, and the import resolves them that way (export writes them back
under the same convention). Resolved against the wrong list, a water table can double
back on itself or land meters off, and the model still solves, with pore pressures that
change the factor of safety by several percent.

SLOPE/W's **phreatic correction** is a setting on the analysis. Where it is on, the head
under the line is reduced by cos² of the line's slope, and the line imports as Type
`phreatic`, which does the same. Where it is off, which is the default, SLOPE/W uses the
plain head under the line, and so does the import. One difference remains: SLOPE/W leaves
the head uncorrected where the line runs along the ground surface, and XSLOPE corrects it
there too. The caveats name that stretch when a model has one. Export writes a `phreatic`
line back as SLOPE/W's switch.

**A line that stops short of the model.** SLOPE/W applies no pore pressure beyond the
ends of a piezometric surface. Its Arai and Tagyo example shows this: the water table
starts partway up the slope, and SLOPE/W's saved pore pressure is zero at every slice
base to the left of that point. XSLOPE needs a line that spans every slice it reads, so
the import continues the line straight down to the model's base and along the base to
the model's edge. That gives the same zero pore pressure, and the notes state where the
line was continued. Continuing the line level instead would put water under the toe
that SLOPE/W does not have: on that example it puts XSLOPE's factor of safety 11% above
SLOPE/W's (the median over SLOPE/W's 251 saved circles), against 0.2% below with the
import's rule. If the water does continue past the line's end in
your model, edit the line.

**Ru** imports as XSLOPE's `ru` option, with each material's own value. Both programs
take the pore pressure as ru times the weight of the soil above the slice base. SLOPE/W
also counts a surcharge in that weight and XSLOPE does not, and the caveats report it
when a model has both.

Ponded water needs no conversion at all. GeoStudio stores no
ponded-water object: where the water rises **above** the ground surface, SLOPE/W simply
carries its weight implicitly, from the water surface itself. XSLOPE models it the same
way — with **Water loads** set to `auto` on the main sheet, the engine derives the
ponded-water load from the model's own water definition at solve time (see
[Automatic water loads](preflight.md#automatic-water-loads)). So the import carries the
**piezometric surface**, sets the mode to `auto`, and writes **no water load**: the two
programs hold the reservoir implicitly, and the dloads sheet of an imported model carries
non-water loads only.

The load itself is unchanged — γ_w × depth, normal to the ground and tapering to zero at
the waterline. Only the source of the load differs. The import caveats state this: *water
loads derive automatically from the imported water definition*. Do not add the reservoir
to the dloads sheet as well; that is the double count preflight warns about.

Where the water comes from a **SEEP/W field**, there is no water surface to carry —
**the file does not state where the water is** — and SLOPE/W still loads the submerged
face. That reservoir has to be recovered from the head field itself, which the import
does (described next), and because nothing downstream can re-derive a load from an
imported field, those models import with **Water loads = `manual`** and the recovered
load written on the dloads sheet. The caveat names the cell. It is the one GeoStudio
import that carries a water load explicitly, because its water surface exists only as a
computed result.

### Pore pressure from a SEEP/W analysis

A SLOPE/W analysis can take its pore pressure from a **parent SEEP/W run** rather than a
piezometric line — a whole finite-element head field, which is the point of using SEEP/W
at all (a piezo line *means* hydrostatic pressure beneath it, and a drawdown or a
transient seepage problem is where that assumption fails).

XSLOPE imports that field. A solved `.gsz` carries the transferred pore pressure in the
stability analysis's **own** result folder, on the mesh SEEP/W solved on, so nothing has
to be re-solved: the mesh becomes XSLOPE's mesh and the nodal pressures become its
`seep_u`, with every material set to the `seep` option. The file must have been saved
**solved** — an unsolved one holds no field, and XSLOPE cannot run SEEP/W's analysis for
it. The import reports this.

Such an analysis records **no water surface anywhere**, and yet SLOPE/W still puts the
weight of the reservoir on the submerged face. It derives that reservoir from the head
field itself, and so does XSLOPE: at any point on the ground surface, water stands to
elevation `y + u/γ_w`.

That single rule reproduces both cases exactly. Under a reservoir the ground is a
specified-head boundary, so `u = γ_w (H − y)` and the implied level comes out at the
reservoir level *H*, to the millimeter, at every submerged point. On a drained or
seepage face `u = 0`, the implied level collapses onto the ground, and no load is
produced. On GeoStudio's own rapid-drawdown example this puts water on the face at the
full-reservoir step and none at any drawn-down step — which is exactly what SLOPE/W's
own per-slice surcharge forces do. In that example the load matters more than the
pressures: of the 13% of factor of safety lost without this treatment, most comes from
the missing reservoir rather than the missing pore pressures.

This is the one import that writes a water load, and it sets **Water loads** to
`manual` so the engine does not derive a second one. XSLOPE's automatic mode reads a
piezometric line or seepage head boundaries; an imported head field is neither, so a
model left on `auto` here would simply have no reservoir. Keep the mode on `manual`
unless you replace the field with a water definition the engine can read.

!!! warning "Importing a transient analysis"
    A SEEP/W run that marches in time saves a result set per time step, and SLOPE/W
    solves the slope at **every one** — the factor of safety in a drawdown problem can
    swing by 50% across them. An XSLOPE model is a single state in time, so the import
    selects one step.

    It selects the **governing** step: the one where SLOPE/W's own factor of safety is
    lowest. The choice, the
    full list of steps and their factors of safety are all reported, and `step='NNN'`
    takes any other.

### Staged rapid drawdown {#staged-rapid-drawdown}

SLOPE/W's staged rapid drawdown uses the three-stage method of Duncan, Wright, and Wong
(1990), the same method as XSLOPE's [rapid drawdown](../lem/rapid.md), and it states the
problem with the same three inputs, so it imports whole:

| SLOPE/W | XSLOPE |
|---|---|
| The higher of the analysis's two piezometric surfaces, the water before drawdown | Piezometric line 1 |
| The lower surface, the water after drawdown | Piezometric line 2 |
| Drawdown strengths c<sub>R</sub> and φ<sub>R</sub> on a Mohr-Coulomb soil that does not drain | d and ψ, the K<sub>c</sub> = 1 envelope |
| A soil with no drawdown strengths | A free-draining material, with d and ψ blank |

SLOPE/W's drawdown strengths are an R envelope drawn tangent to the Mohr circles. SLOPE/W
turns them into the K<sub>c</sub> = 1 envelope with Duncan, Wright, and Brandon's (2014)
equations, and the import uses the same ones (see
[When the lab reports an R envelope](../lem/rapid.md#r-envelope)):

>>$d = c_R \dfrac{\cos \phi_R \cos \phi'}{1 - \sin \phi_R} \qquad \psi = \arctan \dfrac{\sin \phi_R \cos \phi'}{1 - \sin \phi_R}$

An XSLOPE input file does not record which analysis to run, so there is no rapid drawdown
switch to carry. Run the imported model as a rapid drawdown: in Studio, **Run → Run LEM…**
with Rapid drawdown ticked; in a script, `rapid=True`. Run as an ordinary analysis, it is
the slope before drawdown. The critical circle imported from a solved file carries
SLOPE/W's rapid drawdown factor of safety, so the two programs can be compared on it; the
[verification page](../verification/geostudio.md#gs-2-40) does this for SLOPE/W's four
rapid drawdown examples.

The notes report anything that is not simply two surfaces, one above the other: a third
surface, two surfaces that cross, a soil assigned to the lowered surface, drawdown
strengths on a soil that is not Mohr-Coulomb, and a staged rapid drawdown whose water is
not two piezometric surfaces. XSLOPE needs both d and ψ, and reads a zero in either as a
blank, so a soil whose c<sub>R</sub> or φ<sub>R</sub> is zero is named too.

### Tension cracks

A tension crack imports as `tcrack_depth`, with its percent filled with water as
`tcrack_water`. GeoStudio switches a crack **off by omitting the option**, while keeping
the crack line's geometry in the file. The import honors the option, not the leftover
geometry, so a crack the analysis no longer uses is not restored. If those leftover
points were read, they would put a water-filled crack into models that have none — invisible in
a φ = 0 undrained analysis, but worth a few percent in a c′ = 0 one.

### Reinforcement

A GeoStudio reinforcement carries a capacity, a plate capacity, a pullout resistance and
an out-of-plane spacing, but **no direction**: the direction is implied by the
reinforcement type. XSLOPE stores it explicitly, so the import has to supply it, and the
mapping was chosen by comparison with *SLOPE/W's own factors of safety*: reinforcement
acts **along the bar** (`dir='axial'`). On the SLOPE/W verification manual's
geosynthetic-reinforced embankment ([§2.24](../verification/geostudio.md#gs-2-24)) this
reproduces SLOPE/W to within 0.3% for the clay fill and 0.6% for the sand fill.

GeoStudio writes a setting only when it differs from the default, so a missing setting
takes GeoStudio's default:

| GeoStudio | XSLOPE |
|---|---|
| No type | An anchor, GeoStudio's default type |
| F of S Dependent: No (the default) | `appl='active'`, a known load |
| F of S Dependent: Yes | `appl='passive'`, the force divided by the factor of safety |
| Force distribution: Distributed (the default) | Applied where the bar crosses the slip surface, and the notes say so |

**Anchors.** An anchor develops its pullout resistance only over its **bond length**, the
grouted length at its buried end; the free length in front of it carries no bond. So an
anchor's capacity imports as the smaller of its tensile capacity and its pullout
resistance times its bond length, developed over the bond length from the buried end.
Where a slip surface crosses the free length, that is SLOPE/W's anchor force exactly. On
SLOPE/W's Reinforcement with Anchors example, a pullout resistance of 200 kN/m over a
3 m bond, with the anchors 2 m apart, gives 300 kN per meter of wall, the force in
SLOPE/W's own slice results, and the imported model matches SLOPE/W's factor of safety
to within 0.3% on its critical circles. A nail is
bonded over its whole length, so its bond length plays no part.

An anchor or nail can also carry a **shear force**, and GeoStudio has pile and
user-defined reinforcement types. XSLOPE's reinforcement lines carry tension only, so
these are not imported, and the notes name each one. Leaving out an anchor's shear force
lowers the factor of safety: by 5.5% on SLOPE/W's example with shear.

Everything is converted to XSLOPE's **per unit width** convention by dividing through by
the spacing, and GeoStudio's pullout resistance per unit length becomes the bond length
that develops full capacity (`lp = t_max / rate`), which is the same law written the
other way round. The plate capacity becomes the end capacity at whichever end sits at the
face. A file from an earlier GeoStudio 2025 build states the capacity as
`<TensileCapacity>` and gives no plate capacity; GeoStudio 2025.2 writes the same model
with both capacities equal to it, so the import reads it that way and says so.

Where a geosynthetic defines its pullout through **interface adhesion and friction**
rather than a constant resistance, XSLOPE's constant-rate bond length is an approximation
of a stress-dependent law, and the import reports this.

### What does not come across, and why

Import **reports every one of these** — it never drops something that changes the answer
without telling you. Read the caveats it returns.

- **The search definition.** SLOPE/W searches with an entry/exit range or a
  grid-and-radius; neither is an XSLOPE search, and a wrong surface is worse than none.
  If the file was saved *solved*, XSLOPE imports **SLOPE/W's critical circle** instead —
  which makes the model complete, and lets you compare the two programs on identical
  geometry. A search with a moment axis saves the axis in place of each circle's
  center; the import then recovers the critical circle's center from SLOPE/W's saved
  slices, and the other circles of that search are not used. Run a search afterwards to find XSLOPE's own critical surface. If that circle
  will not build on the model (SLOPE/W may have scored a *composite* surface, truncated
  along a bedrock boundary), XSLOPE takes the most critical of SLOPE/W's circles that
  does build, or imports no surface if none does, and the caveats report either case. A
  block, fully specified, or automatic-search surface is not a circle, so no surface is
  imported from one, or from an analysis that takes its surface from such a parent.
  An *optimized* search ends on a non-circular surface, which SLOPE/W saves with a
  circle fitted to it. XSLOPE imports the most critical circle of the search before the
  optimization instead, and the caveats give SLOPE/W's lower, optimized factor of safety.
- **Reinforcement *sets*** — GeoStudio's out-of-plane staging groups. The reinforcement
  itself imports; the grouping does not.
- **Partial factors.** SLOPE/W can solve an analysis with a set of partial factors (a
  Eurocode 7 design approach), dividing the strengths and multiplying the loads by them.
  XSLOPE applies none, so the model imports with the strengths and loads the file states,
  and its factor of safety is not SLOPE/W's factored one. The notes name the set and every
  factor in it that is not 1.
- **Inclined tension cracks.** XSLOPE's crack is vertical. Where GeoStudio's crack line
  sits at a varying depth, the deepest is taken (conservative) and the difference reported.
- **Spatial pore-pressure functions and B-bar.** A spatial function imports as zero.
  XSLOPE has no B-bar, so where SLOPE/W adds B-bar's excess pore pressure to a
  piezometric line, the line imports and the excess pressure does not. The caveats
  report both. The spatial function is one-way: export *writes* one (it is how a
  seepage field crosses), but import does not read one back, so an XSLOPE model
  exported and re-imported returns dry, with a caveat reporting it.
- **The SEEP/W analysis itself** — its mesh is imported, but its conductivity functions
  and boundary conditions are not, so XSLOPE cannot *re-solve* the seepage problem, only
  read the answer SEEP/W already computed. An **unsolved** SEEP/W-fed file therefore
  arrives with no pore pressure at all, and is flagged.
- **Anything XSLOPE does not recognize.** Unrecognized GeoStudio elements are reported
  by name rather than ignored, so a feature added in a later GeoStudio version cannot
  change a factor of safety without notice.

## Unit systems

The same `.gsz` schema holds both metric and imperial models, and GeoStudio declares
which in the project's coordinate settings. XSLOPE reads that and records it on the
imported project's **unit-system declaration** (the same `SI` / `Imperial` setting the
template carries) — but there is nothing to convert. XSLOPE never converts your numbers:
the declaration only fixes the unit weight of water, drives the unit labels on the Studio
forms and plot axes/colorbars, and keeps the model self-consistent (SI = m, kPa, kN/m³;
Imperial = ft, psf, pcf). An import builds a whole new project whose geometry, strengths
and unit weight of water all come from the same file, so they are consistent with each
other by construction; the declaration just names which system they are in. Where a file
omits the attribute, the system is inferred from the unit weight of water: ≈9.807 → SI;
≈62.4 → Imperial.

## Exporting to GeoStudio

=== "Studio"

    **File → Export to GeoStudio (SLOPE/W)…** writes the current model to a `.gsz`.

=== "Script"

    ```python
    from xslope.geostudio import export_gsz
    from xslope.fileio import load_slope_data

    slope_data = load_slope_data("my_model.xlsx")
    caveats = export_gsz(slope_data, "my_model.gsz")
    ```

Material zones become regions, materials become SLOPE/W materials (see
[Strength options on export](#strength-options-on-export)), a piezo line becomes a
piezometric surface (its points written under the
owning analysis, the mirror image of how the import reads them; a `phreatic` line
switches on SLOPE/W's phreatic correction), a material with a γ<sub>sat</sub> different
from its γ is written with SLOPE/W's unit weight (γ<sub>sat</sub>) and its separate unit
weight above the water table (γ), a rapid drawdown model
becomes SLOPE/W's staged rapid drawdown (see
[Exporting a rapid drawdown](#exporting-a-rapid-drawdown)), and **distributed loads
become surcharges** — GeoStudio's surcharge is a wedge of fill, so an XSLOPE pressure
profile is written as fill whose *depth* reproduces the pressure exactly. A **seepage
field** is written as SLOPE/W's spatial pore-pressure function, described below. Export
requires a **polygon-based** model: a profile-line model has no regions to map onto, and
is rejected rather than approximated.

!!! warning "Do not re-create the ponded water"
    Ponded water is the one distributed load that is deliberately **not** written, and it
    must not be added by hand. GeoStudio has no ponded-water object: it *derives* the
    reservoir, and the pressure it puts on the slope, from the **piezometric surface** —
    which is written. Re-creating it as a surcharge would count the water twice. The export
    caveats report this when it happens.

    This is the mirror image of the import, and of XSLOPE's own automatic mode. A model
    with **Water loads** on `auto` has no water block to leave out in the first place:
    its piezometric line goes across, and GeoStudio derives the reservoir from it just
    as XSLOPE's engine does. A `manual` model's water block is recognized and skipped.

A `.gsz` still cannot carry everything XSLOPE models. Failure surfaces, piles, and the
material properties `t_cut`, `phi_b` and `s_cap` are **not written**: SLOPE/W has no
material encoding for them, and inventing one would corrupt the file. Each omission is
reported as a caveat so you know what to re-create on the GeoStudio side.

### Strength options on export

| XSLOPE | GeoStudio |
|---|---|
| `mc` | Mohr-Coulomb, with c′ and φ′ |
| `elastic` | Bedrock (impenetrable) |
| `ratio` | SHANSEP with a constant ratio, c/p |
| `cp` | S=f(depth) |
| `pow`, `hb` | Mohr-Coulomb, with the material's c and φ, which is not the same strength |

SHANSEP has no minimum strength, so a `ratio` material's c is not written. Where
c/p × σ′<sub>v</sub> is below c, SLOPE/W's strength is lower than XSLOPE's, and the
caveats say so.

S=f(depth) measures depth from the top of the material's own layer. A `cp` material
whose top is level and at or below r_elev therefore crosses exactly, with the strength
at that top and the same rate. A `cp` material that reaches above r_elev, where XSLOPE
holds the strength at c, is written to match XSLOPE below r_elev. Above it, SLOPE/W's
strength falls below c. A top that is not level is written from its highest point. Each
case is reported with the size of the difference.

GeoStudio has a Hoek-Brown model, but the export does not write it yet. GeoStudio has
no power-curve model as such. Its Shear-Normal function could hold the curve as a
table of points, which the export does not write.

### Exporting a rapid drawdown

A model with both piezometric lines and a d and ψ on a Mohr-Coulomb material is written
as SLOPE/W's staged rapid drawdown, the way SLOPE/W's own examples write one: the rapid
drawdown switch on the analysis, line 1 and line 2 as its two piezometric surfaces, and
each d and ψ as the drawdown strengths c<sub>R</sub> and φ<sub>R</sub> that SLOPE/W turns
back into the same d and ψ. The XSLOPE file does not say whether a model is to be run as
a rapid drawdown, so a model that carries both is always exported as one. To solve it in
GeoStudio as an ordinary analysis, turn its rapid drawdown off there.

SLOPE/W checks that a soil's c<sub>R</sub> is above its c′ and its φ<sub>R</sub> below its
φ′. XSLOPE has no such rule, so an envelope that fails the check is written as it is, and
the notes name the soil. Distributed loads after drawdown (the `dloads (2)` sheet) are not
written, because SLOPE/W applies one set of surcharges in every stage; the lowered pool
needs no load of its own, since GeoStudio derives it from the lower surface. A model whose
water is a seepage field cannot become a staged rapid drawdown, which takes its water as
two piezometric surfaces, so its d and ψ are left out, and the notes say so.

### Exporting a seepage field

Where the model's pore pressure is a finite-element **seepage field**, the *pressures*
cross and the *seepage model* does not.

SLOPE/W accepts pore pressure as a **spatial function**: a set of discrete points, each
carrying a pressure head, which it interpolates between. That is a data input rather than
a seepage analysis, and it matches what XSLOPE has. So the export writes
one point per node of XSLOPE's seepage mesh, at full resolution, carrying `u / γ_w`, and
points the analysis at it — the intent being that a GeoStudio user opening the file
solves the slope on the same pore pressures XSLOPE solved on. Every tag XSLOPE writes is
taken from GeoStudio's own sample and verification files, including models that carry a
points-only spatial function, and the export is checked against that schema; it has not
been tested in a GeoStudio installation, so check an exported coupled model carefully the
first time it is opened in GeoStudio.

What does **not** cross is the seepage problem: there is no SEEP/W analysis in the file,
no mesh, no conductivity functions and no boundary conditions, so the field cannot be
re-solved or altered on the GeoStudio side — only used. Re-creating the seepage analysis
in SEEP/W is the way to get that back. XSLOPE writes no SEEP/W analysis rather than a
plausible-looking one: an invented mesh with invented properties and a solution it never
computed would be indistinguishable, in the file, from a real one.

A piezometric line is still *not* offered as a substitute. It **means** hydrostatic
pressure beneath it, which is the very assumption a seepage analysis is run to avoid, so
XSLOPE does not substitute one for the other. Where a model defines both, the seepage
field takes precedence and the export reports this.

GeoStudio applies a spatial function to **every
material in the analysis**, so a model where only some materials
draw from the seepage field arrives with all of them doing so. And ponded water — written
as a surcharge here, as everywhere — **must be kept**, which is the opposite of the
piezometric-surface case above: a spatial function gives GeoStudio no water surface, so it
derives no reservoir of its own, and deleting the surcharge would lose the weight of the
water entirely.

Suction (negative pore pressure above the phreatic surface) is written as it stands.
Neither program credits it by default — SLOPE/W only through a φᵇ, which is not written,
and XSLOPE clamps the effective-normal pore pressure at zero — so the two agree unless you
set a φᵇ in GeoStudio. The count of negative points is reported.

The export reports all of this: how many points were written, that they are the pressures
and not the model, the suction count, and the ponded-water rule.

### How export is checked

A round trip through XSLOPE's own reader would not test the export, because the reader and
writer share the same reading of the schema and would agree with each other even where
GeoStudio rejects the file. Export is therefore checked against GeoStudio's own documents,
in both directions:

- **Every tag XSLOPE writes is one GeoStudio writes.** This catches an invented tag.
- **Every tag GeoStudio *always* writes, XSLOPE writes too** — unless it is on an explicit
  list of things XSLOPE deliberately omits, each with a reason. This direction is the more
  important one: the first check cannot catch an **omission**, and an omission is worse, because the file still opens and still looks
  right. A file missing its `ComputedPhysics` block, for instance, draws perfectly and
  names and colors its materials correctly while leaving their strengths unreachable,
  with no warning.

Both directions are checked in the test suite, so a regression fails locally rather than in
GeoStudio. `tools/gsz_export_diff.py` regenerates them: it imports a real `.gsz`, exports
it straight back out, and diffs XSLOPE's XML against the vendor's.

## Notes and limitations

- **Stability analyses only.** A `.gsz` can hold SEEP/W, SIGMA/W and other GeoStudio
  analyses. XSLOPE imports the SLOPE/W ones, and reads a SEEP/W analysis's **solved
  pore-pressure field** where a stability analysis depends on it — but it does not import
  a seepage *problem* to re-solve, and writes no analysis but SLOPE/W.
- Only one piezometric surface is imported, except from a staged rapid drawdown, which
  takes two; any other analysis that defines several uses the first, for every material
  that uses one. SLOPE/W can assign a different surface to each material, and the notes
  say when a file defines more than one.
- Import replaces the current project. In Studio you are prompted to discard unsaved
  changes first, exactly as when opening a file; the imported model is then left unsaved
  so you can review it and **Save As**.
- No third-party package is needed — `.gsz` is a ZIP of XML and a binary PLY mesh, both
  handled by the standard library.
