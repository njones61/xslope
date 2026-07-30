# GeoStudio (SLOPE/W) Import/Export

XSLOPE reads and writes GeoStudio SLOPE/W project files (`.gsz`), so a model built in
GeoStudio can be brought into XSLOPE without redrawing it, and an XSLOPE model can be
handed to a GeoStudio user.

A `.gsz` is a ZIP holding the model as a single XML document. Unlike a DXF — which is
just lines on layers, and needs a human to say what each layer *means* — a `.gsz` is
semantically complete: its regions, materials and water conditions already know what
they are. So the import needs no mapping wizard. The one thing it does ask is **which
analysis** to import.

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

### Why the analysis matters

A `.gsz` usually holds several analyses over one geometry, and they differ in more than
the slip surface. GeoStudio assigns **materials per analysis**, not per region — so the
same slope can be one soil in the first analysis and a different soil in the second. The
choice is not cosmetic, which is why XSLOPE makes you make it rather than guessing.

### What comes across

| GeoStudio | XSLOPE |
|---|---|
| Regions | Material zones (`polygons`) |
| Region→material assignment (per analysis) | `mat_id` |
| Mohr-Coulomb materials (c′, φ′) | `option='mc'` with c′, φ′ |
| **Undrained materials** (φ = 0, strength entered as cohesion) | `option='mc'` with c = Su, φ = 0 |
| **Bedrock** (impenetrable) | `option='elastic'` — an impenetrable material, see below |
| Material colour | The zone colour in Studio |
| Piezometric surface, and the materials assigned to use it | Piezo line, on the materials that use it |
| **Water above the ground surface** | **A distributed load** — see below |
| **A parent SEEP/W analysis's pore-pressure field** | The `seep` option, on SEEP/W's own mesh — see below |
| **Surcharges** (a wedge of fill above the ground — the surcharge pressure is its *unit weight*) | `dloads`, varying with the depth of fill |
| **Reinforcement** — nails, geosynthetics, anchors, tiebacks | `reinforcement_lines` — see below |
| **Line loads** (magnitude + trend/plunge) | `line_loads` |
| **Tension crack** (+ percent filled with water) | `tcrack_depth`, `tcrack_water` |
| Unit weight of water | `gamma_water` |
| Horizontal seismic coefficient | `k_seismic` |
| SLOPE/W's critical circle, if the file was solved and it fits | A single trial circle |

The ground surface and domain are derived from the zones, as for any polygon-based model.

### Materials

Mohr-Coulomb materials carry straight over as `option='mc'` with their c′ and φ′.
GeoStudio's undrained (φ = 0) materials store the undrained strength in the cohesion
field, and import as `option='mc'` with c = S<sub>u</sub> and φ = 0.

GeoStudio's **Bedrock** is impenetrable — slip surfaces cannot enter it. XSLOPE's
**`elastic`** option is the exact analog: an impenetrable zone that cannot fail and that
a slip surface may ride along but never cut into. So a bedrock region is imported as an
elastic material, keeping its unit weight; its c and φ are dropped because it cannot
fail. That reproduces SLOPE/W exactly, and the material round-trips back to Bedrock on
export. Importing bedrock as ordinary soil would let trial surfaces cut straight
through it.

Matric suction is reported rather than converted. SLOPE/W parameterises unsaturated
(negative pore-pressure) shear strength as a **suction ceiling on the piezometric
surface** — a maximum suction, and a cap on how much of it counts — not as a
per-material φᵇ. That is a threshold on the pore pressure, semantically different from
XSLOPE's Fredlund `phi_b` / `s_cap` apparent-cohesion pair, so the import reports it and
leaves the choice to you: set `phi_b` (and `s_cap`) by hand if you need unsaturated
strength credited.

### Water and pore pressure

A piezometric surface's points belong to the **analysis that owns it**, not to the
shared geometry table, and the import resolves them that way (export writes them back
under the same convention). The distinction is silent but not small: resolved against
the wrong list, a water table can double back on itself or land metres off and the model
still solves — at pore pressures several percent wrong in factor of safety.

Ponded water is synthesised on import, because GeoStudio stores no ponded-water object:
where the water rises **above** the ground surface, SLOPE/W simply carries its weight
implicitly. XSLOPE needs that weight to exist explicitly, so the import **creates a
distributed load** of γ_w × depth, normal to the ground and tapering to zero at the
waterline. Without it, the weight of the reservoir would silently vanish.

Where the water comes from a piezometric line, the waterline is drawn. Where it comes
from a SEEP/W field, **nothing in the file says where the water is** — and SLOPE/W still
loads the submerged face, so the water surface has to be recovered from the head field
itself, as described next.

### Pore pressure from a SEEP/W analysis

A SLOPE/W analysis can take its pore pressure from a **parent SEEP/W run** rather than a
piezometric line — a whole finite-element head field, which is the point of using SEEP/W
at all (a piezo line *means* hydrostatic pressure beneath it, and a drawdown or a
transient seepage problem is precisely where that assumption fails).

XSLOPE imports that field. A solved `.gsz` carries the transferred pore pressure in the
stability analysis's **own** result folder, on the mesh SEEP/W solved on, so nothing has
to be re-solved: the mesh becomes XSLOPE's mesh and the nodal pressures become its
`seep_u`, with every material set to the `seep` option. The file must have been saved
**solved** — an unsolved one holds no field, and XSLOPE cannot run SEEP/W's analysis for
it. That is said loudly rather than assumed away.

Such an analysis records **no water surface anywhere**, and yet SLOPE/W still puts the
weight of the reservoir on the submerged face. It derives that reservoir from the head
field itself, and so does XSLOPE: at any point on the ground surface, water stands to
elevation `y + u/γ_w`.

That single rule reproduces both cases exactly. Under a reservoir the ground is a
specified-head boundary, so `u = γ_w (H − y)` and the implied level comes out at the
reservoir level *H*, to the millimetre, at every submerged point. On a drained or
seepage face `u = 0`, the implied level collapses onto the ground, and no load is
produced. On GeoStudio's own rapid-drawdown example this puts water on the face at the
full-reservoir step and none at any drawn-down step — which is exactly what SLOPE/W's
own per-slice surcharge forces do. The load is worth more than the pressures it
accompanies: of the 13% of factor of safety that example loses without this treatment,
**most is the missing reservoir, not the missing pore pressures.**

!!! warning "A transient analysis is many models, and XSLOPE is one"
    A SEEP/W run that marches in time saves a result set per time step, and SLOPE/W
    solves the slope at **every one** — the factor of safety in a drawdown problem can
    swing by 50% across them. An XSLOPE model is a single state in time, so the import
    must choose.

    It chooses the **governing** step: the one where SLOPE/W's own factor of safety is
    lowest, which is the condition a drawdown analysis exists to find. The choice, the
    full list of steps and their factors of safety are all reported, and `step='NNN'`
    takes any other.

### Tension cracks

A tension crack imports as `tcrack_depth`, with its percent filled with water as
`tcrack_water`. GeoStudio switches a crack **off by omitting the option**, while keeping
the crack line's geometry in the file. The import honors the option, not the leftover
geometry, so a crack the analysis no longer uses is not resurrected. Read naively, those
leftover points would put a water-filled crack into models that have none — invisible in
a φ = 0 undrained analysis, but worth a few percent in a c′ = 0 one.

### Reinforcement

A GeoStudio reinforcement carries a capacity, a plate capacity, a pullout resistance and
an out-of-plane spacing — but **no direction and no dependence on the factor of safety**.
Both are implied by the reinforcement type. XSLOPE stores them explicitly, so the import
has to supply them, and the mapping is *measured against SLOPE/W's own factors of safety*
rather than reasoned about: reinforcement acts **along the bar** (`dir='axial'`), as a
known load (`appl='active'`). On the manual's reinforced-embankment benchmark this
reproduces SLOPE/W to within 0.3%, where applying the force tangent to the slip surface
does not converge at all.

Everything is converted to XSLOPE's **per unit width** convention by dividing through by
the spacing, and GeoStudio's pullout resistance per unit length becomes the bond length
that develops full capacity (`lp = t_max / rate`), which is the same law written the
other way round. The plate capacity becomes the end capacity at whichever end sits at the
face.

Where a geosynthetic defines its pullout through **interface adhesion and friction**
rather than a constant resistance, XSLOPE's constant-rate bond length is an approximation
of a stress-dependent law, and the import says so.

### What does not come across, and why

Import **reports every one of these** — it never drops something that changes the answer
without telling you. Read the caveats it returns.

- **The search definition.** SLOPE/W searches with an entry/exit range or a
  grid-and-radius; neither is an XSLOPE search, and a wrong surface is worse than none.
  If the file was saved *solved*, XSLOPE imports **SLOPE/W's critical circle** instead —
  which makes the model complete, and lets you compare the two programs on identical
  geometry. Run a search afterwards to find XSLOPE's own critical surface. If that circle
  will not build on the model (SLOPE/W may have scored a *composite* surface, truncated
  along a bedrock boundary), no surface is imported and you are told so.
- **Reinforcement *sets*** — GeoStudio's out-of-plane staging groups. The reinforcement
  itself imports; the grouping does not.
- **Inclined tension cracks.** XSLOPE's crack is vertical. Where GeoStudio's crack line
  sits at a varying depth, the deepest is taken (conservative) and the difference reported.
- **Vertical seismic coefficient.** XSLOPE models only the horizontal one.
- **Pore pressure that is neither a piezometric surface nor a SEEP/W field.** Ru and
  spatial pore-pressure functions import as zero, and say so. The traffic here is
  one-way: export *writes* a spatial function (it is how a seepage field crosses), but
  import does not read one back, so an XSLOPE model exported and re-imported returns
  dry — loudly, not quietly.
- **The SEEP/W analysis itself** — its mesh is imported, but its conductivity functions
  and boundary conditions are not, so XSLOPE cannot *re-solve* the seepage problem, only
  read the answer SEEP/W already computed. An **unsolved** SEEP/W-fed file therefore
  arrives with no pore pressure at all, and is flagged.
- **Anything XSLOPE has never seen.** Unrecognised GeoStudio elements are reported by
  name rather than ignored, so a feature added in a future version cannot silently change
  someone's factor of safety.

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

Material zones become regions, materials become Mohr-Coulomb materials with their
strengths, a piezo line becomes a piezometric surface (its points written under the
owning analysis, the mirror image of how the import reads them), and **distributed loads
become surcharges** — GeoStudio's surcharge is a wedge of fill, so an XSLOPE pressure
profile is written as fill whose *depth* reproduces the pressure exactly. A **seepage
field** is written as SLOPE/W's spatial pore-pressure function, described below. Export
requires a **polygon-based** model: a profile-line model has no regions to map onto, and
is rejected rather than approximated.

!!! warning "Do not re-create the ponded water"
    Ponded water is the one distributed load that is deliberately **not** written, and it
    must not be added by hand. GeoStudio has no ponded-water object: it *derives* the
    reservoir, and the pressure it puts on the slope, from the **piezometric surface** —
    which is written. Re-creating it as a surcharge would count the water twice. Export
    says so explicitly when it happens.

A `.gsz` still cannot carry everything XSLOPE models. An **elastic** material *is* written
— as GeoStudio's impenetrable **Bedrock** — but failure surfaces, piles, other
non-Mohr-Coulomb strengths (c/φ, power curve, Hoek-Brown), and the v16/v17 material
features `t_cut`, `phi_b` and `s_cap` are **not written**: SLOPE/W has no material
encoding for them, and inventing one would corrupt the file. Each omission is reported as
a caveat so you know what to re-create on the GeoStudio side.

### Exporting a seepage field

Where the model's pore pressure is a finite-element **seepage field**, the *pressures*
cross and the *seepage model* does not — and the difference is worth being precise about.

SLOPE/W accepts pore pressure as a **spatial function**: a set of discrete points, each
carrying a pressure head, which it interpolates between. That is a data input, not a
seepage analysis, and it is exactly the shape of what XSLOPE has. So the export writes
one point per node of XSLOPE's seepage mesh, at full resolution, carrying `u / γ_w`, and
points the analysis at it. A GeoStudio user opening the file solves the slope on the same
pore pressures XSLOPE solved on. (The schema for this comes from GeoStudio's own sample
and verification files, which is where every tag XSLOPE writes comes from.)

What does **not** cross is the seepage problem: there is no SEEP/W analysis in the file,
no mesh, no conductivity functions and no boundary conditions, so the field cannot be
re-solved or altered on the GeoStudio side — only used. Re-creating the seepage analysis
in SEEP/W is the way to get that back. XSLOPE writes no SEEP/W analysis rather than a
plausible-looking one: an invented mesh with invented properties and a solution it never
computed would be indistinguishable, in the file, from a real one.

A piezometric line is still *not* offered as a substitute. It **means** hydrostatic
pressure beneath it, which is the very assumption a seepage analysis is run to avoid, so
XSLOPE will not quietly swap one for the other. Where a model defines both, the seepage
field wins and the export says so.

Two things to know about the result. GeoStudio applies a spatial function to **every
material in the analysis**, not to a chosen subset, so a model where only some materials
draw from the seepage field arrives with all of them doing so. And ponded water — written
as a surcharge here, as everywhere — **must be kept**, which is the opposite of the
piezometric-surface case above: a spatial function gives GeoStudio no water surface, so it
derives no reservoir of its own, and deleting the surcharge would lose the weight of the
water entirely.

Suction (negative pore pressure above the phreatic surface) is written as it stands.
Neither program credits it by default — SLOPE/W only through a φᵇ, which is not written,
and XSLOPE clamps the effective-normal pore pressure at zero — so the two agree unless you
set a φᵇ in GeoStudio. The count of negative points is reported so the question is at
least asked.

The export reports all of this: how many points were written, that they are the pressures
and not the model, the suction count, and the ponded-water rule.

### How export is checked

A round-trip through XSLOPE's own reader proves **nothing** — reader and writer share the
same idea of the schema, agree perfectly with each other, and pass while GeoStudio quietly
refuses half the file. So export is checked against GeoStudio's own documents instead, in
both directions:

- **Every tag XSLOPE writes is one GeoStudio writes.** This catches an invented tag.
- **Every tag GeoStudio *always* writes, XSLOPE writes too** — unless it is on an explicit
  list of things XSLOPE consciously omits, each with a reason. This is the direction that
  matters, and the one that is easy to forget: the first check cannot catch an
  **omission**, and an omission is worse, because the file still opens and still looks
  right. A file missing its `ComputedPhysics` block, for instance, draws perfectly and
  names and colours its materials correctly while leaving their strengths unreachable,
  saying nothing.

Both directions are locked in the test suite, so a regression fails locally rather than in
GeoStudio. `tools/gsz_export_diff.py` regenerates them: it imports a real `.gsz`, exports
it straight back out, and diffs XSLOPE's XML against the vendor's.

## Notes and limitations

- **Stability analyses only.** A `.gsz` can hold SEEP/W, SIGMA/W and other GeoStudio
  analyses. XSLOPE imports the SLOPE/W ones, and reads a SEEP/W analysis's **solved
  pore-pressure field** where a stability analysis depends on it — but it does not import
  a seepage *problem* to re-solve, and writes no analysis but SLOPE/W.
- Only one piezometric surface is imported; a file defining several uses the first.
- Import replaces the current project. In Studio you are prompted to discard unsaved
  changes first, exactly as when opening a file; the imported model is then left unsaved
  so you can review it and **Save As**.
- No third-party package is needed — `.gsz` is a ZIP of XML and a binary PLY mesh, both
  handled by the standard library.
