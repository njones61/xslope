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
| Regions (`Points` + `PointIDs`) | Material zones (`polygons`) |
| Region→material assignment (per analysis) | `mat_id` |
| Mohr-Coulomb materials (`CohesionPrime`, `PhiPrime`) | `option='mc'` with c′, φ′ |
| **Undrained materials** (`UndrainedPhiZero`, strength in `Cohesion`) | `option='mc'` with c = Su, φ = 0 |
| **Bedrock** (impenetrable) | **Excluded from the domain** — see below |
| Material colour | The zone colour in Studio |
| Piezometric surface + `MaterialUsesPiezs` | Piezo line, on the materials that use it |
| **Piezo line above the ground surface** | **A distributed load** — see below |
| **Surcharges** (a wedge of fill above the ground — `Pressure` is its *unit weight*) | `dloads`, varying with the depth of fill |
| **Reinforcement** — nails, geosynthetics, anchors, tiebacks | `reinforcement_lines` — see below |
| **Line loads** (`Value` + trend/plunge) | `line_loads` |
| **Tension crack** (+ `PctFilledWithWater`) | `tcrack_depth`, `tcrack_water` |
| Unit weight of water | `gamma_water` |
| Horizontal seismic coefficient | `k_seismic` |
| SLOPE/W's critical circle, if the file was solved and it fits | A single trial circle |

The ground surface and domain are derived from the zones, as for any polygon-based model.

!!! info "Bedrock becomes geometry, not a material"
    GeoStudio's **Bedrock** is impenetrable — slip surfaces cannot enter it. XSLOPE has no
    such material: its impenetrable boundary is the **domain** itself. So bedrock regions
    are *excluded* from the imported model, and the domain floor lands on the top of the
    bedrock. That reproduces SLOPE/W exactly. Importing bedrock as ordinary soil would let
    trial surfaces cut straight through it.

!!! info "Ponded water is synthesised"
    GeoStudio stores no ponded-water object: where the piezometric line rises **above** the
    ground surface, SLOPE/W simply carries the water's weight implicitly. XSLOPE needs that
    weight to exist explicitly, so the import **creates a distributed load** of γ_w × depth,
    normal to the ground and tapering to zero at the waterline. Without it, the weight of
    the reservoir would silently vanish.

!!! info "Reinforcement: GeoStudio implies the direction, XSLOPE states it"
    A GeoStudio reinforcement carries a capacity, a plate capacity, a pullout resistance
    and an out-of-plane spacing — but **no direction and no F-of-S dependence**. Both are
    implied by the `Type`. XSLOPE stores them explicitly, so the import has to supply them,
    and the mapping is *measured against SLOPE/W's own factors of safety* rather than
    reasoned about: reinforcement acts **along the bar** (`dir='axial'`), as a known load
    (`appl='active'`). On the manual's reinforced-embankment benchmark this reproduces
    SLOPE/W to within 0.3%, where applying the force tangent to the slip surface does not
    converge at all.

    Everything is converted to XSLOPE's **per unit width** convention by dividing through
    by the spacing, and GeoStudio's "pullout resistance per unit length" becomes the bond
    length that develops full capacity (`lp = t_max / rate`), which is the same law written
    the other way round. The plate capacity becomes the end capacity at whichever end sits
    at the face.

    Where a geosynthetic defines its pullout through **interface adhesion and friction**
    rather than a constant resistance, XSLOPE's constant-rate bond length is an
    approximation of a stress-dependent law, and the import says so.

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
- **Pore pressure other than a piezometric surface.** Ru, pressure grids and
  finite-element pore pressures import as zero.
- **Anything XSLOPE has never seen.** Unrecognised GeoStudio elements are reported by
  name rather than ignored, so a feature added in a future version cannot silently change
  someone's factor of safety.

!!! info "Units need no conversion"
    The same `.gsz` schema holds both metric and imperial models, and GeoStudio declares
    which in `<Coordinates><EngCoords UnitSystem="…">`. XSLOPE reads that and **reports
    it**, but there is nothing to convert: XSLOPE has no unit setting anywhere. It is
    unit-agnostic — consistent units in, the same units out. An import builds a whole new
    project whose geometry, strengths and unit weight of water all come from the same
    file, so they are consistent with each other by construction. The caveat exists only
    to tell you which system the numbers you are about to read are in. (Where a file omits
    the attribute, the system is inferred from the unit weight of water: ≈9.807 → kN/m³,
    m, kPa; ≈62.4 → lb/ft³, ft, psf.)

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

Material zones become regions, materials become Mohr-Coulomb materials, and a piezo line
becomes a piezometric surface. Export requires a **polygon-based** model: a profile-line
model has no regions to map onto, and is rejected rather than approximated.

A `.gsz` cannot carry everything XSLOPE models. Failure surfaces, reinforcement, piles,
distributed and line loads, and non-Mohr-Coulomb strengths are **not written**, and each
is reported as a caveat so you know what to re-create on the GeoStudio side.

!!! note "Check exported files in GeoStudio"
    Export is checked two ways: the file round-trips through XSLOPE's own reader with no
    loss, and every XML tag it writes is one that GeoStudio's own files use (a schema
    conformance test, so a tag in the wrong place fails the suite). A round-trip through
    our own reader alone proves nothing — reader and writer can share the same wrong
    assumption — which is why the conformance check exists. Still, GeoStudio is the only
    authority on its own format: open an exported file there before relying on it.

## Notes and limitations

- **SLOPE/W only, for now.** A `.gsz` can also hold SEEP/W, SIGMA/W and other GeoStudio
  analyses; XSLOPE currently reads only the SLOPE/W ones. This is a gap, not a limit of
  the format — see below.
- Only one piezometric surface is imported; a file defining several uses the first.
- Import replaces the current project. In Studio you are prompted to discard unsaved
  changes first, exactly as when opening a file; the imported model is then left unsaved
  so you can review it and **Save As**.
- No third-party package is needed — `.gsz` is a ZIP of XML, both handled by the
  standard library.

### SEEP/W is a natural next step

Nothing about a SEEP/W analysis is unreadable — it sits in the same XML, and a solved
`.gsz` carries its mesh and nodal pore pressures. Those map onto XSLOPE's own seepage
inputs directly: hydraulic conductivities and boundary conditions onto `seepage_bc`, and
a solved head field onto the `seep` pore-pressure option. It simply isn't implemented
yet. If you have SEEP/W models you want to bring over, that is worth knowing — the work
is scoped, not speculative.
