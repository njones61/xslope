# docs_line_loads — recorded assistant session

Captured by `tools/assistant_sessions.py` against the live Studio assistant, offscreen.

- Project: `docs/lem/files/xslope_simple_mult_layers.xlsx`
- Provider / model: `anthropic` / `claude-opus-5`
- Confirm-before-running: off (unattended capture)

## Turn 1 — Does xslope support line loads?

````text
You: Does xslope support line loads?

Ran code:
    r = docs("line load")
    for row in r:
        print(row['score'], row['page'], '|', row['heading'], '|', row['url'], row['anchor'])
        print('   ', row['snippet'][:250].replace('\n',' '))

Output:
    - 2. Loads on the Crest — The same force as a line load  [12.23]
      https://xslope.readthedocs.io/en/latest/tutorials/lem02_loads_on_the_crest/#the-same-force-as-a-line-load
      A **line load** is a concentrated force per unit width, acting at one point on the ground surface: a footing, or the weight of a wall facing. The surcharge above amounts to 750 psf × 10 ft = **7500 lb/ft**. Now we put all of it on the single point at the middle of the strip — one row, in the columns the `lloads` worksh…
    - Ordinary Method of Slices — Complete Formulation  [8.13]
      https://xslope.readthedocs.io/en/latest/lem/oms/#complete-formulation
      …ace <br> $\theta_p$ = angle of pile force from horizontal (positive = counterclockwise/upward) <br> $L$ = line load at point $f$ on the top of the slice <br> $\delta$ = angle of the line load from horizontal (default $-90°$ = straight down) <br> The rest of the forces are the same as before. The **distributed load** re…
    - Input Template — Worksheet: lloads  [8.13]
      https://xslope.readthedocs.io/en/latest/usage/input_template/#worksheet-lloads
      The **lloads** worksheet defines line loads — concentrated forces applied at a point on the ground surface, per unit width of slope. A typical use is the self-weight of a facing element, such as the shotcrete plate of a soil nail wall, applied as a point load on the wall face. (For loads spread over an area of the surf…
    - Overview — Advanced Loading Conditions  [7.95]
      https://xslope.readthedocs.io/en/latest/lem/overview/#advanced-loading-conditions
      …stands: the seepage boundary conditions wherever a seepage analysis is defined, otherwise the piezometric line. The dloads sheets then carry non-water loads only, and the derived load is drawn on every plot in its own colour so a load nobody typed is more visible rather than less. On `manual` — what every file written…
    - Janbu Method — Complete Formulation with Extra Forces  [7.94]
      https://xslope.readthedocs.io/en/latest/lem/janbu/#complete-formulation-with-extra-forces
      …gle $\psi$ from horizontal ($\psi = \alpha$ for tangent/flexible reinforcement, the default; $\psi$ = the line's own inclination for axial/rigid reinforcement) <br> $T$ = tension crack water force <br> $H$ = pile/pier force at point $e$ on the failure surface <br> $\theta_p$ = angle of pile force from horizontal (posit…
    - Overview — Distributed loads  [7.77]
      https://xslope.readthedocs.io/en/latest/fem/overview/#distributed-loads
      Force boundary conditions in XSLOPE come from the **dloads** sheets: line loads given as a sequence of coordinates with load intensities (force per unit length), shared with the limit-equilibrium solvers, which convert them to a resultant on each slice. Hydrostatic pressure on a submerged face need not be entered at al…
    - Bishop's Simplified Method — Complete Formulation  [7.7]
      https://xslope.readthedocs.io/en/latest/lem/bishop/#complete-formulation
      …gle $\psi$ from horizontal ($\psi = \alpha$ for tangent/flexible reinforcement, the default; $\psi$ = the line's own inclination for axial/rigid reinforcement) <br> $T$ = tension crack water force <br> $H$ = pile/pier force at point $e$ on the failure surface <br> $\theta_p$ = angle of pile force from horizontal (posit…
    - 2. Loads on the Crest — 2. Loads on the Crest  [7.66]
      https://xslope.readthedocs.io/en/latest/tutorials/lem02_loads_on_the_crest/
      --- title: "Tutorial LEM-2 — Loads on the Crest" description: "Put a surcharge on the crest of a slope in XSLOPE — as a distributed load, as a line load, pushed normal to the surface or straight down — add a seismic coefficient, and solve for the cohesion that would carry it all at a target factor of safety." ---
    12.23 2. Loads on the Crest | The same force as a line load | https://xslope.readthedocs.io/en/latest/tutorials/lem02_loads_on_the_crest/ the-same-force-as-a-line-load
        A **line load** is a concentrated force per unit width, acting at one point on the ground surface: a footing, or the weight of a wall facing. The surcharge above amounts to 750 psf × 10 ft = **7500 lb/ft**. Now we put all of it on the single point at
    8.13 Ordinary Method of Slices | Complete Formulation | https://xslope.readthedocs.io/en/latest/lem/oms/ complete-formulation
        …ace <br> $\theta_p$ = angle of pile force from horizontal (positive = counterclockwise/upward) <br> $L$ = line load at point $f$ on the top of the slice <br> $\delta$ = angle of the line load from horizontal (default $-90°$ = straight down) <br> The
    8.13 Input Template | Worksheet: lloads | https://xslope.readthedocs.io/en/latest/usage/input_template/ worksheet-lloads
        The **lloads** worksheet defines line loads — concentrated forces applied at a point on the ground surface, per unit width of slope. A typical use is the self-weight of a facing element, such as the shotcrete plate of a soil nail wall, applied as a p
    7.95 Overview | Advanced Loading Conditions | https://xslope.readthedocs.io/en/latest/lem/overview/ advanced-loading-conditions
        …stands: the seepage boundary conditions wherever a seepage analysis is defined, otherwise the piezometric line. The dloads sheets then carry non-water loads only, and the derived load is drawn on every plot in its own colour so a load nobody typed i
    7.94 Janbu Method | Complete Formulation with Extra Forces | https://xslope.readthedocs.io/en/latest/lem/janbu/ complete-formulation-with-extra-forces
        …gle $\psi$ from horizontal ($\psi = \alpha$ for tangent/flexible reinforcement, the default; $\psi$ = the line's own inclination for axial/rigid reinforcement) <br> $T$ = tension crack water force <br> $H$ = pile/pier force at point $e$ on the failu
    7.77 Overview | Distributed loads | https://xslope.readthedocs.io/en/latest/fem/overview/ distributed-loads
        Force boundary conditions in XSLOPE come from the **dloads** sheets: line loads given as a sequence of coordinates with load intensities (force per unit length), shared with the limit-equilibrium solvers, which convert them to a resultant on each sli
    7.7 Bishop's Simplified Method | Complete Formulation | https://xslope.readthedocs.io/en/latest/lem/bishop/ complete-formulation
        …gle $\psi$ from horizontal ($\psi = \alpha$ for tangent/flexible reinforcement, the default; $\psi$ = the line's own inclination for axial/rigid reinforcement) <br> $T$ = tension crack water force <br> $H$ = pile/pier force at point $e$ on the failu
    7.66 2. Loads on the Crest | 2. Loads on the Crest | https://xslope.readthedocs.io/en/latest/tutorials/lem02_loads_on_the_crest/ 
        --- title: "Tutorial LEM-2 — Loads on the Crest" description: "Put a surcharge on the crest of a slope in XSLOPE — as a distributed load, as a line load, pushed normal to the surface or straight down — add a seismic coefficient, and solve for the coh

Assistant: Yes. Line loads are a first-class input — the `lloads` worksheet, `slope_data['line_loads']` in memory.

**What it is:** a concentrated force per unit width of slope applied *at a point on the ground surface* — a footing, or the self-weight of a facing element such as a soil nail wall's shotcrete plate. (Loads spread over an area of the surface go on the `dloads` sheets instead.)

**Record:**

```python
slope_data['line_loads'].append(
    {'x': 40.0, 'y': 20.0, 'P': 7500.0, 'angle': -90, 'label': 'footing'}
)
```

| Field | Meaning |
|:--|:--|
| `x`, `y` | the application point — must lie **on** the ground surface (the loader snaps a small mismatch, refuses a big one), so take it from the profile, not a sketch dimension |
| `P` | magnitude, force per unit width (lb/ft here, kN/m in SI). Must be **positive** |
| `angle` | degrees from horizontal, default `-90` = straight down |
| `label` | free text for plots/report |

**Who reads it:** every LEM method carries it — in the slice equations it is the term `L` at point `f` on top of the slice, acting at angle `delta` from horizontal (see the OMS complete formulation, https://xslope.readthedocs.io/en/latest/lem/oms/#complete-formulation). The FEM applies it as a concentrated force on the nearest mesh node.

**Worked example:** Tutorial LEM-2 takes a 750 psf × 10 ft crest strip = 7500 lb/ft and re-applies the identical resultant as a single line load, comparing the two — https://xslope.readthedocs.io/en/latest/tutorials/lem02_loads_on_the_crest/#the-same-force-as-a-line-load. Field reference: https://xslope.readthedocs.io/en/latest/usage/input_template/#worksheet-lloads

Your current model has no line loads. Want me to add one?
````

Tokens: input 40,323 · cached_input 18,723 · output 839 · calls 2 · 11.8s

## Session total

- Tokens: input 40,323 · cached_input 18,723 · output 839 · calls 2
- Wall time: 16.3s
