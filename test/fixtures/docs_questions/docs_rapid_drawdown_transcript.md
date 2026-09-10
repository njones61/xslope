# docs_rapid_drawdown — recorded assistant session

Captured by `tools/assistant_sessions.py` against the live Studio assistant, offscreen.

- Project: `docs/lem/files/xslope_simple_mult_layers.xlsx`
- Provider / model: `anthropic` / `claude-opus-5`
- Confirm-before-running: off (unattended capture)

## Turn 1 — Can xslope do rapid drawdown?

````text
You: Can xslope do rapid drawdown?

Ran code:
    r = docs('rapid drawdown', limit=8)
    for row in r:
        print(row['page'], '|', row['heading'], '|', row['url'], row['anchor'])
    print('---CORPUS---')
    for row in corpus_index('rapid drawdown'):
        print(row)

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
    2. Rapid Drawdown | 2. Rapid Drawdown | https://xslope.readthedocs.io/en/latest/tutorials/combo02_rapid_drawdown/ 
    Overview | Rapid Drawdown Analysis | https://xslope.readthedocs.io/en/latest/lem/overview/ rapid-drawdown-analysis
    Rapid Drawdown Analysis | Rapid Drawdown Analysis | https://xslope.readthedocs.io/en/latest/lem/rapid/ rapid-drawdown-analysis
    2. Rapid Drawdown | What rapid drawdown is | https://xslope.readthedocs.io/en/latest/tutorials/combo02_rapid_drawdown/ what-rapid-drawdown-is
    Rapid Drawdown Analysis | When Does Rapid Drawdown Apply? | https://xslope.readthedocs.io/en/latest/lem/rapid/ when-does-rapid-drawdown-apply
    GeoStudio (SLOPE/W) Corpus | 🟢 2.43 — Rapid Drawdown – Pilarcitos Dam | https://xslope.readthedocs.io/en/latest/verification/geostudio/ gs-2-43
    Sample Problems | 12. Rapid Drawdown (Johnson Reservoir Dam) | https://xslope.readthedocs.io/en/latest/lem/samples/ 12-rapid-drawdown-johnson-reservoir-dam
    GeoStudio (SLOPE/W) Corpus | 🟢 2.40 — Rapid Drawdown – Walter Bouldin Dam | https://xslope.readthedocs.io/en/latest/verification/geostudio/ gs-2-40
    ---CORPUS---
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
    {'topic': 'rapid-drawdown', 'title': '16. Saturated vs. Moist Unit Weight (γ_sat)', 'url': 'https://xslope.readthedocs.io/en/latest/lem/samples/#16-saturated-vs-moist-unit-weight-_sat'}
    {'topic': 'rapid-drawdown', 'title': '12. Rapid Drawdown (Johnson Reservoir Dam)', 'url': 'https://xslope.readthedocs.io/en/latest/lem/samples/#12-rapid-drawdown-johnson-reservoir-dam'}
    {'topic': 'rapid-drawdown', 'title': 'Tutorial COMBO-2 — Rapid Drawdown', 'url': 'https://xslope.readthedocs.io/en/latest/tutorials/combo02_rapid_drawdown/#tutorial-combo-2-rapid-drawdown'}
    {'topic': 'rapid-drawdown', 'title': 'The dam', 'url': 'https://xslope.readthedocs.io/en/latest/tutorials/combo02_rapid_drawdown/#the-dam'}
    {'topic': 'transient-seepage', 'title': 'SEEPW-T03 — Rapid drawdown', 'url': 'https://xslope.readthedocs.io/en/latest/verification/geostudio/#seepw-t03'}
    {'topic': 'rapid-drawdown', 'title': 'VP96: Embankment dam, homogenous, rapid drawdown, water table', 'url': 'https://xslope.readthedocs.io/en/latest/verification/rocscience/#vp96'}
    {'topic': 'rapid-drawdown', 'title': 'VP97: Embankment dam, homogenous, rapid drawdown, water table', 'url': 'https://xslope.readthedocs.io/en/latest/verification/rocscience/#vp97'}
    {'topic': 'rapid-drawdown', 'title': 'VP98: Walter Bouldin Dam rapid drawdown (Duncan, Wright & Wong 1990)', 'url': 'https://xslope.readthedocs.io/en/latest/verification/rocscience/#vp98'}
    {'topic': 'rapid-drawdown', 'title': 'VP99: Pumped-storage project dam rapid drawdown (DWW 1990)', 'url': 'https://xslope.readthedocs.io/en/latest/verification/rocscience/#vp99'}
    {'topic': 'homogeneous', 'title': 'VP100: Embankment dam, homogenous, rapid drawdown, water table', 'url': 'https://xslope.readthedocs.io/en/latest/verification/rocscience/#vp100'}

Assistant: Yes — it's a first-class feature, not a workaround.

## What it implements

The **Duncan, Wright & Wong (1990) three-stage procedure**, run on top of any of the standard methods (OMS, Bishop, Janbu, Corps, Lowe, Spencer, M-P). Stage 1 solves the pre-drawdown (full pool) condition to get consolidation stresses on the failure surface; stage 2 re-solves with the pool lowered using undrained strengths derived from those stresses; stage 3 takes the lower of drained/undrained strength slice by slice. Full write-up: [`lem/rapid/`](https://xslope.readthedocs.io/en/latest/lem/rapid/#rapid-drawdown-analysis), including the *T = c_v·t / D²* screening rubric for **when** drawdown is actually "rapid" ([when-does-rapid-drawdown-apply](https://xslope.readthedocs.io/en/latest/lem/rapid/#when-does-rapid-drawdown-apply)).

## What you supply

| Input | Meaning |
|:--|:--|
| `piezo_line` / `piezo_line2` | pool before / after drawdown |
| `dloads` / `dloads2` | stage-1 / stage-2 surface loads (skip on `water_loads='auto'`) |
| `seepage_bc` / `seepage_bc2` | the seepage alternative — two steady solves instead of sketched lines |
| material `c`, `phi`, `cp`, `r_elev`, `d`, `psi` | drained envelope plus the R-envelope for the undrained stage |

Three ways to state where the water is — a sketched piezometric pair, two steady seepage solutions, or a transient march — all shown in [Tutorial COMBO-2, Johnson Reservoir Dam](https://xslope.readthedocs.io/en/latest/tutorials/combo02_rapid_drawdown/) (pool dropped 160 → 110 ft; Spencer FS = 1.498, [sample 12](https://xslope.readthedocs.io/en/latest/lem/samples/#12-rapid-drawdown-johnson-reservoir-dam)).

## How you run it

```python
res = run_lem(rapid=True, search=True)          # sketched piezo pair
res = run_lem(rapid=True, seep_time='model')    # transient: uses tseep stage_1/stage_2
```

## Verification

Two documented drawdown failures are in the corpus with published comparisons: [Walter Bouldin Dam](https://xslope.readthedocs.io/en/latest/verification/geostudio/#gs-2-40) (1975 failure) and [Pilarcitos Dam](https://xslope.readthedocs.io/en/latest/verification/geostudio/#gs-2-43), both against SLOPE/W and the original DWW numbers, plus Rocscience cases [VP96–VP100](https://xslope.readthedocs.io/en/latest/verification/rocscience/#vp98).

Your current model (`xslope_simple_mult_layers`) has no piezometric lines at all, so it isn't set up for drawdown — say the word if you'd like me to add the two pool levels.
````

Tokens: input 41,626 · cached_input 37,019 · output 1,647 · calls 2 · 23.7s

## Session total

- Tokens: input 41,626 · cached_input 37,019 · output 1,647 · calls 2
- Wall time: 26.0s
