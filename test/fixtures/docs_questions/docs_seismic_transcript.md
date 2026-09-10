# docs_seismic — recorded assistant session

Captured by `tools/assistant_sessions.py` against the live Studio assistant, offscreen.

- Project: `docs/lem/files/xslope_simple_mult_layers.xlsx`
- Provider / model: `anthropic` / `claude-opus-5`
- Confirm-before-running: off (unattended capture)

## Turn 1 — Can it do seismic loading?

```text
You: Can it do seismic loading?

Ran code:
    docs('seismic loading k_seismic pseudo-static')

Output:
    - GeoStudio (SLOPE/W) Corpus — 🟢 2.4 — ACADS Non-Homogeneous + Seismic  [27.92]
      https://xslope.readthedocs.io/en/latest/verification/geostudio/#gs-2-4
      ACADS 1(d): the three-material non-homogeneous slope of §2.3 with a horizontal pseudo-static seismic coefficient of 0.15 added, verifying seismic loading over layered strengths. **Input:** [vp004.xlsx](files/rocscience/vp004.xlsx) · **Rocscience detail:** [VP4](rocscience.md#vp4) | Method | XSLOPE | SLOPE/W | ACADS | |…
    - Overview — Seismic forces  [25.83]
      https://xslope.readthedocs.io/en/latest/fem/overview/#seismic-forces
      Seismic loading uses the pseudo-static method in both the limit-equilibrium and finite element solvers: a constant horizontal acceleration, expressed as a fraction $k$ of gravity, applied to the whole soil mass as an additional body force, >>$b_{x,seismic} = k \gamma, \qquad b_{y,seismic} = 0$ integrated to nodal force…
    - Overview — Equilibrium  [23.23]
      https://xslope.readthedocs.io/en/latest/fem/overview/#equilibrium
      …onents and $b_x$, $b_y$ are body forces — gravity, $b_x = 0$ and $b_y = -\gamma$, plus the pseudo-static [seismic](#seismic-forces) term when one is applied.
    - Rocscience RS2 (SSRM) Corpus — 🟢 RS2-3: Non-homogeneous slope with seismic load (0.15g)  [22.43]
      https://xslope.readthedocs.io/en/latest/verification/rs2/#rs2-3
      …| 1.00 (−5.2%) | k is entered negative per the FEM sign convention — this is a left-facing slope, so the pseudo-static force acts in −x, while the LEM takes the magnitude and directs it from the failure surface.
    - Published Problems — 🟢 FHWA Example E7 — seismic loading on the E4 wall  [21.82]
      https://xslope.readthedocs.io/en/latest/verification/published/#fhwa-e7
      …cceleration within the reinforced zone out to k<sub>av</sub> = 0.211 g, and the model carries that as its seismic coefficient. Internal Step 4 reduces $F^{*}$ to 80 percent of its static value for seismic loading, so every level's Delta is taken on 0.8F\*. The static resistances the reduction applies to are E4's own, c…
    - 2. Loads on the Crest — A second kind of demand  [19.97]
      https://xslope.readthedocs.io/en/latest/tutorials/lem02_loads_on_the_crest/#a-second-kind-of-demand
      A **seismic coefficient** applies a horizontal force of k × W to every slice, at the slice's center of gravity — the pseudo-static way of asking what an earthquake would do. It is a single global number: - **Studio** — open **Global parameters** and set **Seismic coefficient k**. - **Excel** — `main!D13`. - **Assistant…
    - Overview — Advanced Loading Conditions  [18.93]
      https://xslope.readthedocs.io/en/latest/lem/overview/#advanced-loading-conditions
      …ificantly affect the stability calculations. These advanced loading conditions include distributed loads, seismic loads, reinforcement forces, and tension cracks. **Distributed Loads**: Surface loads such as traffic, construction equipment, or other surcharges can be represented as distributed forces acting on the top…
    - GeoStudio (SLOPE/W) Corpus — 🟢 2.39 — Loukidis – Seismic Coefficient #2  [18.74]
      https://xslope.readthedocs.io/en/latest/verification/geostudio/#gs-2-39
      …ides, and a strong φ = 45° base (c = 15, γ = 19). The target is not a factor of safety but the **critical seismic coefficient** k꜀, the horizontal pseudo-static coefficient at which the searched minimum FS = 1. **Input:** [vp063.xlsx](files/rocscience/vp063.xlsx) · **Rocscience detail:** [VP63](rocscience.md#vp63) · **…

Assistant: Yes — pseudo-static seismic loading, in both the LEM and the FEM.

**How it works.** A single global horizontal coefficient `k` (a fraction of g) applies a force `k·W` to every slice at its center of gravity in the LEM; in the FEM it enters as an extra body force `bx = k·γ`, `by = 0` (`fem/overview/#seismic-forces`). There is no vertical component and no time-history option — it is the classic pseudo-static method, described for the LEM at `lem/overview/#advanced-loading-conditions`.

**Where you set it.** One number, not a per-material or per-slice input:

- Studio → **Global parameters** → *Seismic coefficient k*
- Excel → `main!D13`
- Here → `slope_data['k_seismic'] = 0.15`

Sign convention: the LEM takes the magnitude and directs the force away from the failure surface (driving), so a positive `k` is what you normally want. The FEM is direction-sensitive — for a left-facing slope `k` is entered negative, since the pseudo-static force must act in −x (`verification/rs2/#rs2-3`).

**Verified against.** ACADS problem 1(d) — a three-layer non-homogeneous slope with k = 0.15 — is in the corpus at `verification/geostudio/#gs-2-4`, compared to SLOPE/W and the ACADS reference; the same problem run by SSRM is `verification/rs2/#rs2-3`. There is also a critical-seismic-coefficient study (find the k that drives FS to 1.0) at `verification/geostudio/#gs-2-39` — that one is a `design_sweep({'global': 'k_seismic'}, ..., target_fs=1.0)` job.

Your model currently has no seismic coefficient set. Want me to add one and run it?
```

Tokens: input 39,173 · cached_input 37,115 · output 694 · calls 2 · 10.8s

## Session total

- Tokens: input 39,173 · cached_input 37,115 · output 694 · calls 2
- Wall time: 13.2s
