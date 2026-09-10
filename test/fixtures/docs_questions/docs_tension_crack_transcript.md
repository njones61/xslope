# docs_tension_crack — recorded assistant session

Captured by `tools/assistant_sessions.py` against the live Studio assistant, offscreen.

- Project: `docs/lem/files/xslope_simple_mult_layers.xlsx`
- Provider / model: `anthropic` / `claude-opus-5`
- Confirm-before-running: off (unattended capture)

## Turn 1 — Does xslope have a tension crack?

```text
You: Does xslope have a tension crack?

Ran code:
    docs('tension crack')

Output:
    - GeoStudio (SLOPE/W) Corpus — 🟢 2.2 — ACADS Tension Crack  [16.28]
      https://xslope.readthedocs.io/en/latest/verification/geostudio/#gs-2-2
      ACADS problem 1(b): a homogeneous slope with a water-filled tension crack, verifying tension-crack handling in the limit-equilibrium solution. **Input:** [vp002.xlsx](files/rocscience/vp002.xlsx) · **Rocscience detail:** [VP2](rocscience.md#vp2) | Method | XSLOPE | SLOPE/W | Slide | ACADS band | |---|---|---|---|---| |…
    - Sample Problems — 14. Tension Crack  [16.24]
      https://xslope.readthedocs.io/en/latest/lem/samples/#14-tension-crack
      A slope whose upper layer has cohesion, so an unmodified analysis produces non-physical **tension at the crest** (and an inverted line of thrust) that unconservatively raises the factor of safety. The remedy is a **tension crack** at the top of the slope, whose depth follows $d_{crack} = \dfrac{2 c_d}{\gamma}\tan\!\lef…
    - GeoStudio Import/Export — Tension cracks  [16.09]
      https://xslope.readthedocs.io/en/latest/usage/geostudio/#tension-cracks
      A tension crack imports as `tcrack_depth`, with its percent filled with water as `tcrack_water`. GeoStudio switches a crack **off by omitting the option**, while keeping the crack line's geometry in the file. The import honors the option, not the leftover geometry, so a crack the analysis no longer uses is not resurrec…
    - Rocscience RS2 (SSRM) Corpus — 🟢 RS2-33: Homogeneous slope with tension crack and water table (P&D test slope 2)  [15.99]
      https://xslope.readthedocs.io/en/latest/verification/rs2/#rs2-33
      …| Eight-program LEM table | |---|---|---|---| | SSRM | 1.269 | 1.28 (−0.9%) | 1.03–1.32 | The model's dry tension crack has no FEM representation, worth ~2–3% here.
    - Rocscience Slide2 Corpus — 🟢 VP2: Slope, homogenous, tension crack  [15.84]
      https://xslope.readthedocs.io/en/latest/verification/rocscience/#vp2
      ACADS 1(b) (Giam & Donald 1989): the 1(a) slope with c'=32, phi'=10, gamma=20 and a water-filled tension crack of depth 2c/(gamma*sqrt(ka)) [Craig 1997]. **Input files:** [vp002.xlsx](files/rocscience/vp002.xlsx) | Method | XSLOPE | Slide | SLOPE/W | Giam | |---|---|---|---|---| | Bishop | 1.589 | 1.596 (−0.4%) | 1.664…
    - Rocscience Slide2 Corpus — 🟢 VP60: Soil-nailed wall with tension crack and surcharges  [15.61]
      https://xslope.readthedocs.io/en/latest/verification/rocscience/#vp60
      …sf, φ = 0) carrying a 250-psf crest surcharge plus a 500-psf strip over the first 7.3 ft, with a dry 7-ft tension crack and five passive nail rows at 15° (25,918 lb tensile at 5-ft spacing, bond 1,508 lb/ft). Evaluated on Slide's printed critical circle, truncated by the crack at its printed endpoint; at that geometry…
    - Rocscience RS2 (SSRM) Corpus — The clay case: RS2's tension crack is geometry  [15.43]
      https://xslope.readthedocs.io/en/latest/verification/rs2/#the-clay-case-rs2s-tension-crack-is-geometry
      …kment as a separate file (`slope stability #029_clay`, published SSR **0.99**), and that model states the tension crack **geometrically** rather than constitutively. The crest is cut 2.06 m down from the 9.0 m embankment crest, exactly the theoretical crack depth 2c/γ = 2 × 20 / 19.4 = 2.06 m, and the removed wedge's w…
    - Rocscience Slide2 Corpus — 🟢 VP52: Slope, (4) materials, water table, tension crack  [15.32]
      https://xslope.readthedocs.io/en/latest/verification/rocscience/#vp52
      Slide #52, dry and wet (Table 52.2 water table). An unconstrained circular search lands in the deep (surface 3) family, where Slide runs a grid search and Zhu a specified deep circle, and the published values spread widely. **Input files:** [vp052a.xlsx](files/rocscience/vp052a.xlsx), [vp052b.xlsx](files/rocscience/vp0…

Assistant: Yes — tension cracks are a first-class input.

**How you specify one** — two scalars on the main sheet, which you can set directly on `slope_data`:

| Key | Meaning |
|---|---|
| `tcrack_depth` | depth of the crack below the crest ground surface (length units) |
| `tcrack_water` | fraction (percent) of the crack filled with water |

The crack truncates the trial surface at the crest end, removes the shear resistance over that depth, and — when filled — applies the hydrostatic thrust of the water in the crack as a driving force.

**Why it exists:** in a cohesive upper layer, an unmodified LEM solution develops non-physical *tension* at the crest (and an inverted line of thrust), which unconservatively raises FS. The documented crack depth is

> d_crack = 2 c_d / gamma × tan(45 + phi_d/2)  ... (the sample page gives the full form)

See `lem/samples/#14-tension-crack` (https://xslope.readthedocs.io/en/latest/lem/samples/#14-tension-crack).

**Verified against other codes:**
- ACADS problem 1(b), water-filled crack — Bishop FS 1.589 in XSLOPE vs 1.596 Slide, 1.664 SLOPE/W: https://xslope.readthedocs.io/en/latest/verification/rocscience/#vp2 (and the SLOPE/W row at `verification/geostudio/#gs-2-2`)
- Four-material slope with water table and crack: `verification/rocscience/#vp52`
- A caution from the FEM corpus: a dry tension crack has **no FEM/SSRM representation** — worth ~2–3% on the RS2-33 comparison (https://xslope.readthedocs.io/en/latest/verification/rs2/#rs2-33). It's an LEM feature.

GeoStudio imports also carry it: the crack comes in as `tcrack_depth`/`tcrack_water`, honoring the option flag rather than leftover crack geometry (`usage/geostudio/#tension-cracks`).

Your current model has no crack set. If you want one on this two-layer embankment (c = 400 psf, phi = 0, gamma = 130 pcf), say the depth and water fill and I'll add it.
```

Tokens: input 39,360 · cached_input 37,117 · output 838 · calls 2 · 12.6s

## Session total

- Tokens: input 39,360 · cached_input 37,117 · output 838 · calls 2
- Wall time: 14.6s
