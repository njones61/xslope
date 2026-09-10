# docs_line_loads — recorded assistant session

Captured by `tools/assistant_sessions.py` against the live Studio assistant, offscreen.

- Project: `docs/lem/files/xslope_simple_mult_layers.xlsx`
- Provider / model: `anthropic` / `claude-opus-5`
- Confirm-before-running: off (unattended capture)

## Turn 1 — Does xslope support line loads?

```text
You: Does xslope support line loads?

Ran code:
    docs('line loads')

Output:
    - 2. Loads on the Crest — Tutorial LEM-2 — Loads on the Crest  [9.07]
      https://xslope.readthedocs.io/en/latest/tutorials/lem02_loads_on_the_crest/#tutorial-lem-2-loads-on-the-crest
      …m-obj" markdown> **Objectives** — Learn how to apply loads to a slope model: how to enter distributed and line loads and set their direction, what a seismic coefficient does, how each load moves the factor of safety and the critical surface, and how to sweep strength against a target factor of safety. </div> <p><span c…
    - Input Template — Worksheet: main  [9.02]
      https://xslope.readthedocs.io/en/latest/usage/input_template/#worksheet-main
      …et/tab as appropriate. Template version: Units: Imperial mat Material properties. One per set per profile line Time: day profile Profile lines describing the slope geometry Unit weight of water: polygon Polygons describing the slope geometry Tension crack depth: piezo Piezometric line to define pore pressures Depth of…
    - 2. Loads on the Crest — The same force as a line load  [8.96]
      https://xslope.readthedocs.io/en/latest/tutorials/lem02_loads_on_the_crest/#the-same-force-as-a-line-load
      A **line load** is a concentrated force per unit width, acting at one point on the ground surface: a footing, or the weight of a wall facing. The surcharge above amounts to 750 psf × 10 ft = **7500 lb/ft**. Now we put all of it on the single point at the middle of the strip — one row, in the columns the `lloads` worksh…
    - Overview — Advanced Loading Conditions  [8.82]
      https://xslope.readthedocs.io/en/latest/lem/overview/#advanced-loading-conditions
      …stands: the seepage boundary conditions wherever a seepage analysis is defined, otherwise the piezometric line. The dloads sheets then carry non-water loads only, and the derived load is drawn on every plot in its own colour so a load nobody typed is more visible rather than less. On `manual` — what every file written…
    - Editing Inputs — Double-click to edit  [8.72]
      https://xslope.readthedocs.io/en/latest/studio/editing/#double-click-to-edit
      …he click to the drawing's coordinates and hit-tests the geometry, so this works on: - profile and polygon lines (opens the geometry editor at that line), - the max-depth line, piezometric lines, distributed loads, line loads, reinforcement, piles, - failure circles and non-circular points, - seepage boundary conditions…
    - Input Template — Worksheet: lloads  [8.7]
      https://xslope.readthedocs.io/en/latest/usage/input_template/#worksheet-lloads
      The **lloads** worksheet defines line loads — concentrated forces applied at a point on the ground surface, per unit width of slope. A typical use is the self-weight of a facing element, such as the shotcrete plate of a soil nail wall, applied as a point load on the wall face. (For loads spread over an area of the surf…
    - Overview — Distributed loads  [8.64]
      https://xslope.readthedocs.io/en/latest/fem/overview/#distributed-loads
      Force boundary conditions in XSLOPE come from the **dloads** sheets: line loads given as a sequence of coordinates with load intensities (force per unit length), shared with the limit-equilibrium solvers, which convert them to a resultant on each slice. Hydrostatic pressure on a submerged face need not be entered at al…
    - 2. Loads on the Crest — Conclusion  [8.61]
      https://xslope.readthedocs.io/en/latest/tutorials/lem02_loads_on_the_crest/#conclusion
      This tutorial covered: - Distributed loads along the ground surface and line loads at a point, each with a direction. - The difference direction makes — a normal pressure against a vertical dead weight — and the slopes where it matters. - The seismic coefficient as a second kind of demand, applied to every slice and al…

Assistant: Yes — line loads are a first-class input.

**What they are:** a concentrated force per unit width applied at a single point on the ground surface — a footing, or the self-weight of a facing element (the classic case is the shotcrete plate of a soil nail wall). They live on the `lloads` worksheet ([Input Template — Worksheet: lloads](https://xslope.readthedocs.io/en/latest/usage/input_template/#worksheet-lloads)).

**Fields** (`slope_data['line_loads'][i]`):

| Key | Meaning |
|:---|:---|
| `x`, `y` | the point of application — must lie **on** the ground surface (a small mismatch is snapped, a large one refused) |
| `P` | magnitude, force per unit width, always **positive** |
| `angle` | degrees from horizontal, default `-90` (straight down) |
| `label` | optional name |

**Who reads them:** every LEM method converts the force to a resultant on the slice containing the point; the FEM applies it as a concentrated force at the nearest mesh node ([lem/overview — Advanced Loading Conditions](https://xslope.readthedocs.io/en/latest/lem/overview/#advanced-loading-conditions)). In Studio you can double-click a line load on the canvas to edit it ([studio/editing](https://xslope.readthedocs.io/en/latest/studio/editing/#double-click-to-edit)).

**Worked example:** Tutorial LEM-2 takes a 750 psf × 10 ft crest surcharge and re-applies the same 7,500 lb/ft as a single line load, comparing the two — [Loads on the Crest — the same force as a line load](https://xslope.readthedocs.io/en/latest/tutorials/lem02_loads_on_the_crest/#the-same-force-as-a-line-load).

For a load spread over an area of surface, use `dloads` instead. Want me to add one to the current model? I'd need the point, the magnitude, and the direction.
```

Tokens: input 39,024 · cached_input 37,115 · output 724 · calls 2 · 10.0s

## Session total

- Tokens: input 39,024 · cached_input 37,115 · output 724 · calls 2
- Wall time: 14.5s
