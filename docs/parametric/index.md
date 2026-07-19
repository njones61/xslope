# Parametric Studies

Vary the inputs and watch the answer move. XSLOPE groups three closely related studies
under one **Parametric** umbrella — one Studio button, one API family, one parameter
grammar:

- **[Sensitivity](sensitivity.md)** — sweep one or more parameters and see how the factor of
  safety responds. The results feed a family of plots: a tornado, scaled-sensitivity bars, a
  spider plot, and — when the model carries standard deviations — a variance-contribution
  Pareto and a Monte Carlo rank-correlation chart.
- **[Design](design.md)** — sweep one parameter across an explicit range and find the value
  at which FS meets a target (FS = 1.5, say).
- **[Back-Analysis](back_analysis.md)** — the same single-parameter sweep framed as a
  failure investigation: a slide has occurred, so FS = 1.0 is known, and the study
  back-calculates the parameter value consistent with the observed failure.

This is the geotechnical staple — Duncan & Wright present exactly these charts (FS vs
parameter, and tornado diagrams comparing several parameters at their low/high bounds) —
and half of slope-stability judgment is knowing *which* parameter matters on a given slope.

Sensitivity is deliberately distinct from [reliability analysis](../reliability/taylor.md): the
Taylor-series method perturbs parameters by ±σ to estimate the *distribution* of FS,
which requires the parameters to be independent. A sensitivity sweep asserts nothing
statistical — it simply evaluates the model across a range — which also makes it the
right tool for **correlated fit coefficients** (the power-curve `A` and `b`, for
example) that the reliability method must not treat as independent. The Hoek-Brown
inputs are in the same category: `hb_sci`, `hb_gsi`, `hb_mi` and `hb_d` are all
sweepable here, but they carry no standard-deviation columns and reliability rejects
them, because $m_b$, $s$ and $a$ are all *derived* from GSI, $m_i$ and $D$ — perturbing
them independently would be meaningless.

Sweeps are configured entirely through the API — sensitivity describes an analysis you
run, not a property of the model, so nothing is added to the Excel input template. The same
engine drives the point-and-click
[Parametric study dialog in Studio](../studio/analysis.md#parametric-study)
and the recipes in the [`/xslope` Claude Code skill](../usage/claude/index.md), so a study
set up one way reads the same the others.

## Analysis modes: LEM, FEM, seepage

Every sweep runs through one of three engine `mode=` values: `'lem'` (limit
equilibrium, output FS — the default), `'fem'` (a full finite-element SSRM solve per
point, output FS), or `'seep'` (a rebuilt-and-solved seepage problem per point, output
the **total discharge q**). `mode='fem'` and `mode='seep'` both need a mesh in
`slope_data['mesh']`. Not every function in the family takes a `mode=` argument —
the two plots that reuse the reliability machinery directly are LEM-only:

| Function | `lem` | `fem` | `seep` |
|---|:---:|:---:|:---:|
| `sensitivity()` | yes | yes | yes |
| `tornado()` | yes | yes | yes |
| `design()` | yes | yes | yes |
| `back_analysis()` | yes | yes | yes |
| `scaled_sensitivity()` | yes | yes | yes |
| `variance_contribution()` | yes | — | — |
| `mc_rank_correlation()` | yes | — | — |

`variance_contribution()` and `mc_rank_correlation()` reuse the LEM Taylor-series and
Monte Carlo [reliability](../reliability/index.md) machinery directly (`xslope.advanced.reliability`
and `reliability_mc`), which is a limit-equilibrium path only — see
[Monte Carlo in xslope](../reliability/monte_carlo.md#monte-carlo-in-xslope) for why. They take
no `mode=` argument at all rather than accepting one they cannot honor.

## Addressing a parameter

A sweep needs an unambiguous name for the thing being varied. Parameter references are
strings of the form `"kind:name:field"`:

| kind | name | fields | example |
|---|---|---|---|
| `mat` | material name | strength fields valid for the material's `option` (`c`, `phi` for `mc`; `c`, `cp` for `cp`), plus `gamma`, `gamma_sat`, `ru`, `d`, `psi` | `"mat:Clay:c"` |
| `reinforce` | line label | `t_max`, `t_res`, `lp1`, `lp2`, `tend1`, `tend2`, `spacing` | `"reinforce:Row 2:t_max"` |
| `piles` | pile label | `H`, `theta`, `D`, `S`, `V_cap`, `M_cap` | `"piles:Pile 1:H"` |
| `global` | — | `k_seismic`, `tcrack_depth`, `tcrack_water` | `"global:k_seismic"` |
| `seep` | material name | `k1`, `k2`, `alpha`, `kr0`, `h0` | `"seep:Sand:k1"` |
| `geom` | — | `piezo:dy` — shift the piezometric line vertically (the value is a *delta*) | `"geom:piezo:dy"` |

Every reference is validated against the loaded model before anything runs: an unknown
kind, an unmatched name, or a field the material's strength option does not use raises
immediately, naming what was given and what exists. Names resolve case-insensitively;
duplicated names are an error rather than a guess. Sweeping `c` on a `cp` material is
rejected for the same reason the reliability module rejects it — the model would not
read the value, and the sweep would silently report zero sensitivity.

Note that `gamma` and `gamma_sat` are the same soil weighed two ways, so sweeping
`gamma` moves `gamma_sat` by the same absolute delta (the same coupling the reliability
module applies); `gamma_sat` remains separately addressable when that is what you mean.

## Discovering and specifying parameters

Rather than hand-write every reference, `list_params(slope_data)` enumerates every sweepable
parameter in the loaded model — each material's option-aware strength and general fields,
plus the global `k_seismic` — as plain dicts. It is the menu a GUI parameter-picker (or an
assistant driving the API) chooses from, so a reference is never guessed:

```python
from xslope.sensitivity import list_params

for p in list_params(slope_data):
    print(p['ref'], p['value'], p['sigma'])
    # mat:Soil:c        3.0   1.8
    # mat:Soil:phi      19.6  2.744
    # mat:Soil:gamma    20.0  1.2
    # global:k_seismic  0.0   None
```

Each entry carries `ref` (the canonical string), `kind`, `name`, `index` (the 1-based
material index), `field`, `value` (the current value, or `None` if unset), `sigma` (the
reliability standard deviation — `sigma_c`, `sigma_phi`, … — if the model carries a non-zero
one, so a picker can offer a one-click ±σ range), and a short `label`. Blank or zero-valued
fields are still listed, so a design study can target them with explicit bounds.

Anywhere a `"kind:name:field"` string is accepted, the entry points also accept the
equivalent **dict** or **tuple** — often what a GUI or an assistant naturally produces. All
forms resolve to the same setter and validate identically:

```python
design(slope_data, {"material": "Soil", "property": "c"}, low=6, high=18)  # dict, by name
design(slope_data, {"material": 1, "property": "phi"}, low=15, high=25)     # 1-based index
design(slope_data, {"global": "k_seismic"}, low=0.0, high=0.3)             # a global
design(slope_data, ("mat", "Soil", "c"), low=6, high=18)                    # tuple form
```

The dict accepts `ref` (passed straight through), `material`/`name` and `property`/`field`,
or a `global` key; a material may be named or given by 1-based index.
