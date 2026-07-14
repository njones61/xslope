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
| Materials — Mohr-Coulomb `UnitWeight`, `CohesionPrime`, `PhiPrime` | Materials — γ, c′, φ′ |
| Region→material assignment (per analysis) | `mat_id` |
| Piezometric surface | Piezo line, with materials set to `piezo` |
| Unit weight of water | `gamma_water` |
| Seismic coefficient | `k_seismic` |
| SLOPE/W's critical circle, if the file was solved | A single trial circle |

The ground surface and domain are derived from the zones, as they are for any
polygon-based model.

### What does not, and why

Import reports these rather than dropping them quietly — read the caveats it returns.

- **The search definition.** SLOPE/W searches with an entry/exit range or a
  grid-and-radius; neither is an XSLOPE search, and a wrong surface is worse than none.
  If the file was saved *solved*, XSLOPE imports **SLOPE/W's critical circle** instead —
  which makes the model complete, and lets you compare the two programs on identical
  geometry. Run a search afterwards to find XSLOPE's own critical surface.
- **Non-Mohr-Coulomb strengths.** Anything else (undrained, Hoek-Brown, curved
  envelopes) is imported with its c and φ and flagged. Check it before solving.
- **Pore pressure other than a piezometric surface.** Ru, pressure grids and
  finite-element pore pressures have no direct mapping and import as zero.
- **Reinforcement, piles, and loads.**

!!! warning "Units are implicit"
    A `.gsz` carries no unit-system field — the numbers are just numbers, and the format
    holds both metric and imperial models. XSLOPE infers the system from the unit weight
    of water (≈9.807 → kN/m³, m, kPa; ≈62.4 → lb/ft³, ft, psf) and **refuses to guess**
    on anything else, because silently assuming the wrong system would rescale the whole
    model. The unit system is reported as a caveat on every import — confirm it matches
    your template.

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

!!! note "Export is verified against XSLOPE's reader, not against GeoStudio"
    The exported file round-trips faithfully through XSLOPE's own `.gsz` reader —
    geometry, materials, water and seismic coefficient all survive unchanged. Whether
    GeoStudio itself opens the file has not been verified against the application.
    Treat export as a convenience, and check the result in GeoStudio before relying on it.

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
