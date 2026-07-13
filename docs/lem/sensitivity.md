# Sensitivity Analysis

Vary one input over a range and report how the factor of safety moves. This is the
geotechnical staple — Duncan & Wright present exactly these charts (FS vs parameter, and
tornado diagrams comparing several parameters at their low/high bounds) — and half of
slope-stability judgment is knowing *which* parameter matters on a given slope.

Sensitivity is deliberately distinct from [reliability analysis](reliability.md): the
Taylor-series method perturbs parameters by ±σ to estimate the *distribution* of FS,
which requires the parameters to be independent. A sensitivity sweep asserts nothing
statistical — it simply evaluates the model across a range — which also makes it the
right tool for **correlated fit coefficients** (the power-curve `A` and `b`, for
example) that the reliability method must not treat as independent.

Sweeps are configured entirely through the API — sensitivity describes an analysis you
run, not a property of the model, so nothing is added to the Excel input template.

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

## Running a sweep

```python
from xslope import load_slope_data
from xslope.sensitivity import sensitivity
from xslope.plot import plot_sensitivity

slope_data = load_slope_data("docs/lem/files/xslope_acads_simple.xlsx")

success, result = sensitivity(
    slope_data,
    param="mat:Soil:c",
    rel_range=0.5, n=9,          # ±50% about the base value, 9 points
    # values=[...],              # ...or give the values explicitly
    methods=("bishop", "spencer"),
    search=True,                 # re-search the critical surface per point
)
plot_sensitivity(result['df'], target_fs=1.2)
```

`search=True` is the default and the honest setting: **the critical surface moves as
parameters change**, and re-solving a fixed surface silently understates sensitivity.
`search=False` re-solves the stored surface (`circles[0]` or the non-circular polyline)
instead — roughly fifty times faster, and correct when the question is "given this
surface" (prescribed-surface benchmarks, for example).

The result is a tidy long-format DataFrame — one row per (value × method), with the
unmodified model included as a flagged `is_base` row:

| column | meaning |
|---|---|
| `param` | canonical reference, e.g. `"mat:Soil:c"` |
| `value`, `rel` | the swept value and `value / base_value` |
| `is_base` | this row is the unmodified model |
| `method`, `fs` | solver name and factor of safety |
| `success`, `msg` | per-point outcome — a failed point is a row, not an exception, so a sweep that breaks at value 7 of 9 still reports 1–6 (finding where things break is often the point) |
| `Xo`, `Yo`, `R` | the critical circle per point (searched sweeps), so a jump of the critical surface is visible in the data — `plot_sensitivity` draws jumped points open |

## Sweeping anything else: `modify=`

For changes that are not a single stored scalar — geometry above all — pass a callable
instead of a parameter reference. The callable receives a copy of the model and the
swept value, and owns whatever consistency its edit requires (rebuilding the material
polygons and ground surface after moving profile points, keeping reinforcement anchored
to a moved face):

```python
import math
from shapely.geometry import Polygon
from xslope.fileio import build_ground_surface_from_polygons
from xslope.mesh import build_polygons

def set_slope_angle(sd, beta_deg):
    """Rotate the slope face about the toe to beta_deg, then rebuild the
    polygon geometry (slice weights come from the polygons, not the raw
    profile lines — see the Slope Design page for this pattern)."""
    prof = sd['profile_lines'][0]['coords']
    x_toe, y_toe = prof[1]
    _, y_top = prof[2]
    prof[2] = (x_toe + (y_top - y_toe) / math.tan(math.radians(beta_deg)), y_top)
    polys = [{'polygon': Polygon(p['coords']), 'mat_id': p['mat_id']}
             for p in build_polygons(slope_data={'profile_lines': sd['profile_lines'],
                                                 'max_depth': sd.get('max_depth')})]
    sd['polygons'] = polys
    sd['ground_surface'], sd['domain_polygon'] = build_ground_surface_from_polygons(polys)
    return sd

success, result = sensitivity(slope_data, modify=set_slope_angle,
                              label="slope angle (deg)",
                              values=[20, 22.5, 25, 27.5, 30],
                              methods=("bishop",))
```

Built-in references and `modify=` callables are one code path — every `param` reference
resolves internally to exactly the setter signature `modify=` accepts — and the engine
validates the modified model at every point (polygon validity, ground surface present)
precisely because setters may be user-written: a broken edit becomes a `success=False`
row naming what broke, never a silently inconsistent answer.

## Tornado diagrams

To compare several parameters at once, `tornado()` evaluates each at its low and high
bound (two solves per parameter plus one shared base case) and `plot_tornado` draws the
Duncan-style summary — horizontal bars of the FS swing, sorted so the parameter that
matters most sits on top:

```python
from xslope.sensitivity import tornado
from xslope.plot import plot_tornado

success, result = tornado(
    slope_data,
    ["mat:Soil:c", "mat:Soil:phi", "mat:Soil:gamma"],
    rel_range=0.25,              # ±25% per parameter...
    # bounds={"mat:Soil:c": (2.0, 5.0)},   # ...or explicit bounds per ref
    method="bishop",
)
plot_tornado(result)
```

For the shipped ACADS sample (a weak c–φ soil), the tornado ranks φ far ahead of c and
γ — the ±25% φ band alone swings FS from 0.77 to 1.21 across FS = 1.

## Worked sample

The sweep below reproduces the base case of the ACADS simple slope sample
([xslope_acads_simple.xlsx](files/xslope_acads_simple.xlsx), the same file used
throughout these pages) and brackets it over ±50% of the cohesion, re-searching the
critical surface at every point. The regression suite locks the end points: at c = 1.5
the critical circle has retreated toward the face (the surface-jump columns show the
radius growing as cohesion falls), and at c = 4.5 the slope crosses FS = 1.

<!-- test: file=files/xslope_acads_simple.xlsx, type=sensitivity, param=mat:Soil:c, method=bishop, num_slices=30, n=3, rel_range=0.5, expected_base=0.985, expected_low=0.882, expected_high=1.073, tolerance=0.01, benchmark=SENS-acads-c -->

| c (kPa) | Bishop FS |
|---|---|
| 1.5 (−50%) | 0.882 |
| 3.0 (base) | 0.985 |
| 4.5 (+50%) | 1.073 |

A `geom:piezo:dy` sweep on any of the water-table samples answers the most common
geometry question — "how sensitive is this slope to the water table?" — without writing
a setter: the reference shifts the piezometric line vertically by the swept delta.

For a published benchmark of the sweep itself, see
[verification problem VP40](../verification/rocscience.md#vp40) — Perry (1993)'s
power-curve sensitivity study, where the swept ΔFS matches Slide's published curves
within about a percent at every range endpoint.
