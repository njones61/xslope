# docs_piezo_line — recorded assistant session

Captured by `tools/assistant_sessions.py` against the live Studio assistant, offscreen.

- Project: `docs/lem/files/xslope_simple_mult_layers.xlsx`
- Provider / model: `anthropic` / `claude-opus-5`
- Confirm-before-running: off (unattended capture)

## Turn 1 — How do I add a piezometric line?

````text
You: How do I add a piezometric line?

Ran code:
    r = docs('piezometric line')

Output:
    - 2. Rapid Drawdown — Clearing the piezometric lines  [13.24]
      https://xslope.readthedocs.io/en/latest/tutorials/combo02_rapid_drawdown/#clearing-the-piezometric-lines
      …s replace both of the jobs they did. Where a seepage analysis is defined, its head boundaries outrank the piezometric line for the water loads, and the `u` column set to `seep` outranks it for pore pressure — a line answers only for a model that has neither. Left on the file, the pair would be read by nothing, so delet…
    - 4. Water in the Slope — 2. The piezometric line  [13.23]
      https://xslope.readthedocs.io/en/latest/tutorials/lem04_water_in_the_slope/#2-the-piezometric-line
      Click **Piezometric lines**. The editor holds the two lines on tabs — the water table on **Line 1**, and **Line 2 (rapid drawdown)**, which stays empty here. Press **Add row** eight times and enter the points from the table above, left to right. The preview draws the line on the section as you type, which is where a mi…
    - Seep-Slope Integration — Piezometric line vs. seepage-derived pore pressures  [13.23]
      https://xslope.readthedocs.io/en/latest/seep/seep_slope/#piezometric-line-vs-seepage-derived-pore-pressures
      A piezometric line is a drawn estimate of the phreatic surface; a seepage solution honours the flow field and the conductivities. The right-facing slope below (a three-layer silt / sand / clay profile) is solved both ways on the **same single circle** ($X_0 = 175$, $Y_0 = 100$), so the pore-pressure model is the only d…
    - Seep-Slope Integration — Part a — piezometric line  [13.2]
      https://xslope.readthedocs.io/en/latest/seep/seep_slope/#part-a-piezometric-line
      All three materials use the `piezo` option, with the surface entered on the **piezo** sheet. Input file: [xslope_rface_PIEZO_KEY.xlsx](files/xslope_rface_PIEZO_KEY.xlsx) {width=900} {width=900} Spencer's method gives **FS = 1.94**.
    - Rocscience Slide2 Corpus — 🟢 VP76: Homogeneous dam, FE seepage vs. piezometric line (D&W Fig. 7.19)  [13.04]
      https://xslope.readthedocs.io/en/latest/verification/rocscience/#vp76
      …the entire upstream face. As in VP71, pore pressures are modeled two ways — an FE seepage analysis and a piezometric line — and the critical circle is found by free search. **Input files:** [vp076a.xlsx](files/rocscience/vp076a.xlsx) (FE seepage), [vp076b.xlsx](files/rocscience/vp076b.xlsx) (piezometric line) | Case |…
    - Rocscience Slide2 Corpus — 🟢 VP71: FE seepage vs. piezometric line, same slope (D&W Figs. 6.37–6.38)  [12.99]
      https://xslope.readthedocs.io/en/latest/verification/rocscience/#vp71
      …ope is solved two ways — pore pressures from a finite-element seepage analysis, and pore pressures from a piezometric line — and the two must agree. Case 1 runs XSLOPE's own FE seepage solver (specified heads of 40 and 75 on the two vertical boundaries, the ground surface an exit face) and feeds the nodal pore pressure…
    - Rocscience Slide2 Corpus — 🟢 VP69: Steady-seepage dam with a piezometric line (USACE example F-6)  [12.89]
      https://xslope.readthedocs.io/en/latest/verification/rocscience/#vp69
      …a granular foundation (c' = 0, φ' = 35°, γ = 125 pcf) under steady seepage. Pore pressures come from the piezometric line, which leaves the pool surface at el. 100, drops to the chimney drain, follows it to the tailwater at el. 22.5 and runs out flat to the downstream face; USACE takes u as the vertical distance to th…
    - 2. Rapid Drawdown — Part 1 — Two piezometric lines  [12.55]
      https://xslope.readthedocs.io/en/latest/tutorials/combo02_rapid_drawdown/#part-1-two-piezometric-lines
      The first way to describe the two water states is also the simplest: a piezometric line for each, one at full pool and one drawn down, sketched from the piezometers or from judgment, with no seepage run behind either. The file already carries both lines, so we look at them first and then run the procedure on the pair.

Assistant: Two pieces are needed: the **line** itself, and the materials that **read** it.

## In Studio (GUI)

Click **Piezometric lines** in the inputs dock. The editor has two tabs — **Line 1** (the water table) and **Line 2 (rapid drawdown)**, which stays empty unless you're staging a drawdown. Add rows and type the points left to right; the preview draws them on the section as you go. Then, in the **Materials** table, set the `u` column to `piezo` for every material whose pore pressure should come from that line. (Tutorial: [Water in the Slope §2](https://xslope.readthedocs.io/en/latest/tutorials/lem04_water_in_the_slope/#2-the-piezometric-line))

## Or have me do it

```python
slope_data['piezo_line'] = [(-30, 12), (10, 12), (30, 6), (90, 6)]
slope_data['piezo_phreatic'] = False
for m in slope_data['materials']:
    m['u'] = 'piezo'
```

Three things worth knowing:

| Item | Behavior |
|:--|:--|
| `piezo_phreatic` | `False` (default) = true piezometric line, u = γ_w × vertical distance below it. `True` = phreatic surface, that head reduced by cos²(local line inclination). |
| Ponded water | Where the line rises above the ground, the pool load is added automatically as long as `water_loads = 'auto'` (it is, in your model). Do **not** also enter a `dloads` block — that counts the pool twice. |
| Precedence | A material with `u = 'seep'`, or a defined seepage BC, outranks the line. A line answers only for a model with neither. |

Your model currently has `u = 'none'` on both materials and no line, so it's running fully dry.

If you tell me the water level(s) — where the line starts on the left, where it exits the face, and its tailwater elevation — I'll enter it and set the materials to `piezo`. I won't run anything unless you ask.
````

Tokens: input 39,262 · cached_input 37,116 · output 837 · calls 2 · 12.5s

## Session total

- Tokens: input 39,262 · cached_input 37,116 · output 837 · calls 2
- Wall time: 14.7s
