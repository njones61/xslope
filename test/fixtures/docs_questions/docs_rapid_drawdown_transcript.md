# docs_rapid_drawdown — recorded assistant session

Captured by `tools/assistant_sessions.py` against the live Studio assistant, offscreen.

- Project: `docs/lem/files/xslope_simple_mult_layers.xlsx`
- Provider / model: `anthropic` / `claude-opus-5`
- Confirm-before-running: off (unattended capture)

## Turn 1 — Can xslope do rapid drawdown?

````text
You: Can xslope do rapid drawdown?

Ran code:
    docs('rapid drawdown')
    corpus_index('rapid drawdown')

Output:
    - 2. Rapid Drawdown — 2. Rapid Drawdown  [15.24]
      https://xslope.readthedocs.io/en/latest/tutorials/combo02_rapid_drawdown/
      --- title: "Tutorial COMBO-2 — Rapid Drawdown" description: "The Johnson Reservoir dam with its pool lowered 50 ft, from elevation 160 to a residual pool at 110 — the Duncan, Wright and Wong three-stage procedure run from three statements of where the water is: a sketched piezometric pair, two steady seepage solutions,…
    - Overview — Rapid Drawdown Analysis  [15.16]
      https://xslope.readthedocs.io/en/latest/lem/overview/#rapid-drawdown-analysis
      Rapid drawdown analysis represents a specialized application that can use any of the other methods as its foundation. This approach is specifically designed for dam and levee analysis where water level changes create unique stability challenges. The method accounts for undrained conditions and provides multi-stage anal…
    - Rapid Drawdown Analysis — Rapid Drawdown Analysis  [15.11]
      https://xslope.readthedocs.io/en/latest/lem/rapid/#rapid-drawdown-analysis
      When designing a dam or levee, it is crucial to consider the potential for rapid drawdown of the water level. This can occur due to various factors such as dam failure, sluice gate operation, or other emergency situations. When the water level drops rapidly, it can lead to significant changes in the forces acting on th…
    - 2. Rapid Drawdown — What rapid drawdown is  [15.09]
      https://xslope.readthedocs.io/en/latest/tutorials/combo02_rapid_drawdown/#what-rapid-drawdown-is
      …on, $t$ the time over which the pool falls and $D$ the longest path pore water has to travel to get out. [Rapid drawdown analysis](../lem/rapid.md#when-does-rapid-drawdown-apply) carries the rubric and the $c_v$ table: above $T = 3$ the material drains as fast as the pool falls and the drawdown is not rapid for it; bel…
    - Rapid Drawdown Analysis — When Does Rapid Drawdown Apply?  [15.09]
      https://xslope.readthedocs.io/en/latest/lem/rapid/#when-does-rapid-drawdown-apply
      Rapid drawdown occurs when the pool is lowered rapidly enough that pore water pressures in some of the soils cannot dissipate quickly enough to maintain stability. To determine when rapid drawdown may apply, Duncan, et al. (1992) suggested using the dimensionless time factor, T, from consolidation theory. T is defined…
    - GeoStudio (SLOPE/W) Corpus — 🟢 2.43 — Rapid Drawdown – Pilarcitos Dam  [15.06]
      https://xslope.readthedocs.io/en/latest/verification/geostudio/#gs-2-43
      A three-stage rapid-drawdown analysis of Pilarcitos Dam — a homogeneous earthfill embankment (Duncan, Wright & Wong 1990) drawn down from 72 to 37 ft, a documented drawdown-failure case. It verifies XSLOPE's staged rapid-drawdown procedure. **Input:** [vp097.xlsx](files/rocscience/vp097.xlsx) · **Rocscience detail:** […
    - Sample Problems — 12. Rapid Drawdown (Johnson Reservoir Dam)  [15.03]
      https://xslope.readthedocs.io/en/latest/lem/samples/#12-rapid-drawdown-johnson-reservoir-dam
      Built and run step by step in [COMBO-2](../tutorials/combo02_rapid_drawdown.md). {width=900} **Factor of safety by method** (each method's own critical surface): | OMS | Bishop | Janbu | Corps | Lowe | Spencer | M-P | |---:|---:|---:|---:|---:|---:|---:| | 1.247 | 1.439 | 1.366 | 1.719 | 1.548 | 1.498 | 1.510 |
    - GeoStudio (SLOPE/W) Corpus — 🟢 2.40 — Rapid Drawdown – Walter Bouldin Dam  [15.02]
      https://xslope.readthedocs.io/en/latest/verification/geostudio/#gs-2-40
      The Walter Bouldin Dam is a rolled earthfill dam that failed during a rapid drawdown in 1975; this case verifies XSLOPE's Duncan-Wright-Wong three-stage rapid drawdown procedure. **Input:** [vp098.xlsx](files/rocscience/vp098.xlsx) · **Rocscience detail:** [VP98](rocscience.md#vp98) | Method | XSLOPE | SLOPE/W | DWW |…
    - 16. Saturated vs. Moist Unit Weight (γ_sat)  [rapid-drawdown]
      https://xslope.readthedocs.io/en/latest/lem/samples/#16-saturated-vs-moist-unit-weight-_sat
    - 12. Rapid Drawdown (Johnson Reservoir Dam)  [rapid-drawdown]
      https://xslope.readthedocs.io/en/latest/lem/samples/#12-rapid-drawdown-johnson-reservoir-dam
    - Tutorial COMBO-2 — Rapid Drawdown  [rapid-drawdown]
      https://xslope.readthedocs.io/en/latest/tutorials/combo02_rapid_drawdown/#tutorial-combo-2-rapid-drawdown
    - The dam  [rapid-drawdown]
      https://xslope.readthedocs.io/en/latest/tutorials/combo02_rapid_drawdown/#the-dam
    - SEEPW-T03 — Rapid drawdown  [transient-seepage]
      https://xslope.readthedocs.io/en/latest/verification/geostudio/#seepw-t03
    - VP96: Embankment dam, homogenous, rapid drawdown, water table  [rapid-drawdown]
      https://xslope.readthedocs.io/en/latest/verification/rocscience/#vp96
    - VP97: Embankment dam, homogenous, rapid drawdown, water table  [rapid-drawdown]
      https://xslope.readthedocs.io/en/latest/verification/rocscience/#vp97
    - VP98: Walter Bouldin Dam rapid drawdown (Duncan, Wright & Wong 1990)  [rapid-drawdown]
      https://xslope.readthedocs.io/en/latest/verification/rocscience/#vp98
    - VP99: Pumped-storage project dam rapid drawdown (DWW 1990)  [rapid-drawdown]
      https://xslope.readthedocs.io/en/latest/verification/rocscience/#vp99
    - VP100: Embankment dam, homogenous, rapid drawdown, water table  [homogeneous]
      https://xslope.readthedocs.io/en/latest/verification/rocscience/#vp100

Assistant: Yes — it's a first-class capability, using the **Duncan, Wright & Wong (1990) three-stage procedure** (as documented in Duncan & Wright's text).

## How it works

Three stages on the same trial surface:

1. **Stage 1 — before drawdown.** Drained (effective-stress) analysis with the full pool. This gives the consolidation stresses on the base of each slice.
2. **Stage 2 — immediately after drawdown.** Undrained strengths derived from those stage-1 consolidation stresses, with the pool now at the drawn-down level (the low pool still counts as a surface load and as pore pressure below it).
3. **Stage 3 — check.** Where the drained strength at stage-2 stresses is *lower* than the undrained strength, the drained value governs and is substituted. The reported FS is stage 2/3.

Page: `lem/rapid/#rapid-drawdown-analysis`, plus the applicability rubric (the consolidation time factor T = c_v·t/D², with T > 3 meaning the soil drains as fast as the pool falls, so drawdown is *not* rapid for it) at `lem/rapid/#when-does-rapid-drawdown-apply`.

## What you supply

| Input | Meaning |
|:--|:--|
| `piezo_line` | stage 1 — the full pool |
| `piezo_line2` | stage 2 — the drawn-down pool |
| `dloads` / `dloads2` | stage-specific surface loads (leave to `water_loads='auto'` for the pool itself) |
| `seepage_bc` / `seepage_bc2` | the seepage alternative — two steady solves instead of sketched lines |
| materials | `c`, `phi` **and** `cp`, `r_elev`, `d`, `psi` (the R-envelope parameters for the undrained stage) |

Run it with:

```python
res = run_lem(rapid=True, search=True)
```

A transient version exists too: with a `tseep` sheet, `run_lem(rapid=True, seep_time='model')` stages the two drawdown instants the sheet names.

## Worked examples and verification

- **Tutorial COMBO-2 — Johnson Reservoir Dam**, pool lowered 50 ft from el. 160 to 110, built three ways (sketched piezometric pair, two steady seepage solves, transient): `tutorials/combo02_rapid_drawdown/`. Its results table (`lem/samples/#12-rapid-drawdown-johnson-reservoir-dam`) gives Spencer 1.498, Bishop 1.439, OMS 1.247.
- **Pilarcitos Dam** (drawn down 72 → 37 ft, a documented drawdown failure): `verification/geostudio/#gs-2-43`.
- **Walter Bouldin Dam**, which actually failed in a 1975 drawdown, compared against SLOPE/W and the original DWW numbers: `verification/geostudio/#gs-2-40`.

Your current model has no piezometric lines at all, so it isn't set up for a drawdown run — want me to check whether that's something you'd like to add?
````

Tokens: input 40,123 · cached_input 37,119 · output 1,130 · calls 2 · 17.2s

## Session total

- Tokens: input 40,123 · cached_input 37,119 · output 1,130 · calls 2
- Wall time: 19.5s
