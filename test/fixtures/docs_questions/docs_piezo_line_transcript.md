# docs_piezo_line — recorded assistant session

Captured by `tools/assistant_sessions.py` against the live Studio assistant, offscreen.

- Project: `docs/lem/files/xslope_simple_mult_layers.xlsx`
- Provider / model: `anthropic` / `claude-opus-5`
- Confirm-before-running: off (unattended capture)

## Turn 1 — How do I add a piezometric line?

````text
You: How do I add a piezometric line?

Ran code:
    r = docs("piezometric line pore pressure piezo", limit=8)
    for x in r:
        print(x['score'], x['page'], '|', x['heading'], '|', x['url'], x['anchor'])
        print('   ', x['snippet'][:300].replace('\n',' '))

Output:
    - Input Template — Worksheet: piezo  [18.3]
      https://xslope.readthedocs.io/en/latest/usage/input_template/#worksheet-piezo
      …e pressure head is computed: - **piezo** (the default — a blank cell means the same): the line is a true *piezometric line*, i.e., the locus of pressure heads measured by piezometers at the slip surface. The pore pressure at a slice base is the full static head, u = γw · Δz, where Δz is the vertical distance from the s…
    - Seep-Slope Integration — Pore-pressure options  [17.8]
      https://xslope.readthedocs.io/en/latest/seep/seep_slope/#pore-pressure-options
      …----------| | none | No pore pressures — dry soils, or a total stress analysis | | piezo | Derived from a piezometric line defined on the `piezo` sheet | | seep | Interpolated from the finite-element seepage solution | | ru | A pore-pressure ratio applied to the vertical total stress | Materials that are not set to `se…
    - Overview — Pore pressures  [17.6]
      https://xslope.readthedocs.io/en/latest/lem/overview/#pore-pressures
      …to supply the pore pressure $u$, selected per material through the **u** column of the input template. **Piezometric line** (`u = piezo`, the piezo sheet's `Type` = *piezo* or blank). The line is a true piezometric line — the locus of pressure heads a set of piezometers would read at the slip surface. The pore pressur…
    - 4. Water in the Slope — 3. The piezo worksheet  [17.01]
      https://xslope.readthedocs.io/en/latest/tutorials/lem04_water_in_the_slope/#3-the-piezo-worksheet
      …for the **Type**, and enter (or copy-paste) the eight coordinates from the table above, as shown below: A piezometric line is the model's water table, given as its own polyline: the pore pressure at a point on the failure surface is γ<sub>w</sub> times the vertical distance below the line, and zero above it. The sheet'…
    - Overview — Pore pressures  [16.72]
      https://xslope.readthedocs.io/en/latest/fem/overview/#pore-pressure-options
      …-----------------------| | `none` | none | $u = 0$; the yield check is a total-stress check | | `piezo` | piezometric line | $u = \gamma_w (z_{piezo} - y_{gp})$ from the line elevation above the point | | `ru` | pore-pressure ratio | $u = r_u\,\sigma_v$, with $\sigma_v$ the weight of the soil column above the point | |…
    - 5. A Weak Layer, Non-Circular — 3. The piezo worksheet  [16.53]
      https://xslope.readthedocs.io/en/latest/tutorials/lem05_weak_layer_noncircular/#3-the-piezo-worksheet
      …` for the **Type**, and enter (or copy-paste) the two coordinates from the table above, as shown below: A piezometric line is a water table given as its own polyline: pore pressure at a point on the failure surface is γ<sub>w</sub> times the vertical distance below it. It spans the full section because the water table…
    - 4. Water in the Slope — Check its work  [15.66]
      https://xslope.readthedocs.io/en/latest/tutorials/lem04_water_in_the_slope/#check-its-work
      …ure source reads u = 0 everywhere. If any material came back with `none`, say: *"All three soils read the piezometric line — set every pore pressure option to piezo."* - **One piezometric line, eight points, spanning x = 0 to 320.** A line that stops short of the surface's reach is an error at run time, not a warning.…
    - Sample Problems — 16. Saturated vs. Moist Unit Weight (γ_sat)  [15.3]
      https://xslope.readthedocs.io/en/latest/lem/samples/#16-saturated-vs-moist-unit-weight-_sat
      …ight — γ_sat below the water table, γ above it. The water table comes from whichever the model defines: a piezometric line on the **piezo** sheet, or the phreatic surface of a finite-element seepage solution. It belongs to the problem rather than to any one material, so it splits the weight whether or not the material'…
    18.3 Input Template | Worksheet: piezo | https://xslope.readthedocs.io/en/latest/usage/input_template/ worksheet-piezo
        …e pressure head is computed: - **piezo** (the default — a blank cell means the same): the line is a true *piezometric line*, i.e., the locus of pressure heads measured by piezometers at the slip surface. The pore pressure at a slice base is the full static head, u = γw · Δz, where Δz is the vertica
    17.8 Seep-Slope Integration | Pore-pressure options | https://xslope.readthedocs.io/en/latest/seep/seep_slope/ pore-pressure-options
        …----------| | none | No pore pressures — dry soils, or a total stress analysis | | piezo | Derived from a piezometric line defined on the `piezo` sheet | | seep | Interpolated from the finite-element seepage solution | | ru | A pore-pressure ratio applied to the vertical total stress | Materials th
    17.6 Overview | Pore pressures | https://xslope.readthedocs.io/en/latest/lem/overview/ pore-pressures
        …to supply the pore pressure $u$, selected per material through the **u** column of the input template. **Piezometric line** (`u = piezo`, the piezo sheet's `Type` = *piezo* or blank). The line is a true piezometric line — the locus of pressure heads a set of piezometers would read at the slip surfa
    17.01 4. Water in the Slope | 3. The piezo worksheet | https://xslope.readthedocs.io/en/latest/tutorials/lem04_water_in_the_slope/ 3-the-piezo-worksheet
        …for the **Type**, and enter (or copy-paste) the eight coordinates from the table above, as shown below: A piezometric line is the model's water table, given as its own polyline: the pore pressure at a point on the failure surface is γ<sub>w</sub> times the vertical distance below the line, and zero
    16.72 Overview | Pore pressures | https://xslope.readthedocs.io/en/latest/fem/overview/ pore-pressure-options
        …-----------------------| | `none` | none | $u = 0$; the yield check is a total-stress check | | `piezo` | piezometric line | $u = \gamma_w (z_{piezo} - y_{gp})$ from the line elevation above the point | | `ru` | pore-pressure ratio | $u = r_u\,\sigma_v$, with $\sigma_v$ the weight of the soil colum
    16.53 5. A Weak Layer, Non-Circular | 3. The piezo worksheet | https://xslope.readthedocs.io/en/latest/tutorials/lem05_weak_layer_noncircular/ 3-the-piezo-worksheet
        …` for the **Type**, and enter (or copy-paste) the two coordinates from the table above, as shown below: A piezometric line is a water table given as its own polyline: pore pressure at a point on the failure surface is γ<sub>w</sub> times the vertical distance below it. It spans the full section bec
    15.66 4. Water in the Slope | Check its work | https://xslope.readthedocs.io/en/latest/tutorials/lem04_water_in_the_slope/ check-its-work
        …ure source reads u = 0 everywhere. If any material came back with `none`, say: *"All three soils read the piezometric line — set every pore pressure option to piezo."* - **One piezometric line, eight points, spanning x = 0 to 320.** A line that stops short of the surface's reach is an error at run 
    15.3 Sample Problems | 16. Saturated vs. Moist Unit Weight (γ_sat) | https://xslope.readthedocs.io/en/latest/lem/samples/ 16-saturated-vs-moist-unit-weight-_sat
        …ight — γ_sat below the water table, γ above it. The water table comes from whichever the model defines: a piezometric line on the **piezo** sheet, or the phreatic surface of a finite-element seepage solution. It belongs to the problem rather than to any one material, so it splits the weight whether

Assistant: Two things have to be set — the line, and the materials that read it. Setting only one is the usual mistake.

## 1. The line

In Studio it's the **piezo** worksheet/tab: a polyline of (x, y) points, left to right, spanning the **full width** of the section (yours is x = -30 to 90 — a line that stops short of the failure surface's reach is a run-time error, not a warning). The sheet's **Type** column picks how head is computed:

| Type | Meaning |
|:--|:--|
| `piezo` (or blank, the default) | true piezometric line — u = γw · Δz, the full vertical distance from the line down to the slice base |
| `phreatic` | that head reduced by cos²(local inclination of the line) — the sloping-water-table correction |

## 2. Each material's pore-pressure option

On the materials sheet set the **u** column to `piezo` for every material that should feel it. Materials left at `none` read u = 0 no matter what line exists. (Other options: `seep` = FE seepage field, `ru` = pore-pressure ratio.)

Your two soils are both φ = 0 at `u = none`, i.e. a total-stress analysis — adding a piezo line changes nothing until you switch them to `piezo`, and in a φ = 0 undrained run you normally *shouldn't*.

## Doing it from here

If you'd rather I set it, it's a small edit to `slope_data`:

```python
slope_data['piezo_line'] = [(-30, 12), (10, 12), (40, 4), (90, 4)]
slope_data['piezo_phreatic'] = False        # True = phreatic correction
for m in slope_data['materials']:
    m['u'] = 'piezo'
```

Two related rules worth knowing: use **total** unit weights, never buoyant ones (set `gamma_sat` and the weight is split at the line automatically); and if the water stands *above* the ground surface anywhere, `water_loads = 'auto'` derives the hydrostatic surface load itself — don't add a `dload` for the pool.

Docs: [Input Template — piezo worksheet](https://xslope.readthedocs.io/en/latest/usage/input_template/#worksheet-piezo) and [LEM Overview — Pore pressures](https://xslope.readthedocs.io/en/latest/lem/overview/#pore-pressures); worked walkthrough in [Tutorial 4 — Water in the Slope](https://xslope.readthedocs.io/en/latest/tutorials/lem04_water_in_the_slope/#3-the-piezo-worksheet).

Tell me the water elevations you want and I'll add the line to this model.
````

Tokens: input 40,471 · cached_input 37,016 · output 1,229 · calls 2 · 17.5s

## Session total

- Tokens: input 40,471 · cached_input 37,016 · output 1,229 · calls 2
- Wall time: 19.5s
