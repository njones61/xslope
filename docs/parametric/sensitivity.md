# Sensitivity

See [Addressing a parameter](index.md#addressing-a-parameter) and
[Discovering and specifying parameters](index.md#discovering-and-specifying-parameters) for the
parameter grammar `sensitivity()` and every other function on this page shares with
[Design](design.md) and [Back-Analysis](back_analysis.md).

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

**A searched sweep honours the model's own search window.** If the circles sheet declares
entry and exit ranges, a center box, a maximum tangent depth or a minimum slip depth, every
point is searched inside those limits — exactly as Studio's Run LEM path reads them, so a
windowed model gives the same surface family from a script as from the interface. Pass
`search_opts={...}` to set or override any limit per call (a circular-search keyword wins
over the file), or `use_file_window=False` to search unconstrained regardless of what the
file declares. The same two settings are available on `design`, `back_analysis`, `tornado`,
`scaled_sensitivity`, `fs_vs_time`, `reliability` and `reliability_mc`, and they mean the same
thing on each — every path that searches on the model's behalf reads the window the same way.

The window matters because a sweep answers "how far does FS move when this input moves", and
it answers by re-searching at every point. Left unconstrained on a slope with competing
minima — a benched face, an embankment on a soft foundation — one point can settle in a deep
foundation circle where its neighbour found a shallow one, and the step between them reads as
sensitivity to the parameter when it is really a change in what was measured. The `Xo`/`Yo`/`R`
columns are there to make such a jump visible when no window is declared; declaring one keeps
it from happening.

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

`plot_sensitivity` draws one line per method with FS = 1 (and an optional `target_fs`) as
guide lines, marks the unmodified model as a labelled **base case** entry in the legend — so
the black square in the plot reads as `base case (value, FS = …)` rather than an unexplained
point — and draws any point where the critical surface jumped as an open circle.

`mode=` selects the engine that evaluates each swept point (see
[Analysis modes: LEM, FEM, seepage](index.md#analysis-modes-lem-fem-seepage)). `mode='lem'`
(the default) runs the limit-equilibrium methods and reports a factor of safety. `mode='fem'` runs a
full SSRM solve at every point — the output is still a factor of safety, but each step is
a complete finite-element analysis (expect minutes per step; run it in the background) and
the model must carry a mesh in `slope_data['mesh']`. `mode='seep'` rebuilds and solves the
seepage problem at every point and reports the **total discharge q** through the section
(the same flow rate the seepage solver returns), and likewise needs a mesh. The result's
value column stays `fs` for a consistent schema; `result['output']` (`'FS'` or `'q'`) and
`result['output_label']` say what it really is, so a discharge sweep is never mislabelled a
factor of safety.

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
resolves internally to exactly the setter signature `modify=` accepts.

### What is checked, and when

A sweep validates in two stages.

**Once, before the first point,** the base model goes through the full
[preflight input check](../usage/preflight.md) for whichever engine the sweep will run.
A model that cannot be analysed at all is refused there, with the field named — rather
than failing identically at every one of nine points.

**Then, per point,** the substituted value is re-checked against only the rules that read
the field being swept, and the modified model against the geometry (polygon validity,
ground surface present) because setters may be user-written. A value that would produce a
wrong answer — a unit weight at or below zero, a negative cohesion, a friction angle at or
above 90° — makes that point a `success=False` row **carrying the reason as its `msg`**,
and the sweep carries on. Points before and after it are unaffected.

Finally, each point's result is screened before it is recorded: a factor of safety of 9999
is the search's no-admissible-surface flag rather than an answer, and a negative one is not
a result either. Both become failed rows with a stated reason instead of numbers that would
dominate every plot and derivative built from the sweep.

Pass `check_inputs=False` to skip both validation stages. The result screen still runs — a
sentinel recorded as a success is a defect in the record, not a strictness setting.

## Factor of safety versus time

A [transient seepage](../seep/transient.md) run produces a *sequence* of pore-pressure
fields, one per saved instant. `fs_vs_time` runs the stability analysis against each of
them in turn and tabulates the result — the coupled output a drawdown or an infiltration
study is actually run for, and the answer to *when* the slope is at its weakest rather
than only *how weak*:

```python
from xslope.sensitivity import fs_vs_time
from xslope.plot import plot_sensitivity

success, result = fs_vs_time(
    slope_data, transient_solution,   # from run_transient_seepage / import_transient_solution
    methods=("spencer",),
    # times=[0, 21600, 86400],        # ...or a chosen set instead of every saved frame
)
if success:
    print(f"lowest FS = {result['min_fs']:.3f} at t = {result['critical_time']:g}")
    plot_sensitivity(result['df'])
```

It is a sweep like every other mode on this page — run the model N times and tabulate —
and it returns the same tidy long-format DataFrame, so `plot_sensitivity` draws it with no
special handling. One difference is worth naming: **no input is modified at any step.**
Each point solves the same model against a different *computed* field, so there is no
substituted value to validate and no base case to compare against. The axis is time, and
the `param` column reads `time`.

Alongside `df` the result carries `critical_time` and `min_fs` — the instant of the lowest
factor of safety and its value — plus `critical` (the same pair per method when several
ran), `times`, `n_failed`, and `solution`. `xslope.plot.plot_fs_vs_time` draws the result
directly: the curve, its lowest instant annotated, and the model's drawdown schedule behind
it. In Studio the same run is the Parametric dialog's
[Factor of safety vs time](../studio/analysis.md#factor-of-safety-vs-time) mode, where the
instants are ticked off the loaded march and the curve opens in an **FS vs Time** tab.

When the march is done the run prints its own table — one aligned row per evaluated instant,
carrying the factor of safety, the circle it was found on, and (on a drawdown curve) the
three stage values and which stage governed. An instant that produced nothing carries its
reason on its row. The same rows come back as `result['table']` (a list of dicts) and the
printed block as `result['table_text']`, so a report or a notebook reads what the console
showed rather than rebuilding it. `print_table=False` keeps the rows and prints nothing.

### Each instant as a rapid drawdown

`rapid=True` turns every instant into a full three-stage
[Duncan, Wright and Brandon drawdown](../lem/rapid.md) rather than a single-stage
analysis:

```python
ok, result = fs_vs_time(slope_data, transient_solution, rapid=True)
```

Stage 1 is the march's **initial** state — the `tseep` sheet's `stage_1` (normally t = 0 at
full pool), or the earliest saved frame where the sheet names none — and stage 2 is the frame
at *t*, with stage 3 re-checking the drawn-down section against drained strengths where they
are the lower. The reported value is the drawdown's own, the lower of stages 2 and 3, so the
curve answers *how safe is this slope if the pool falls to where it stands at t* — asked
against pore pressures the march computed rather than a single assumed drawn-down state. The
rows carry `stage1_FS`, `stage2_FS`, `stage3_FS`, `stage3_run` and `governs` beside `fs`, and
`plot_fs_vs_time` draws the three stages as thin dashed lines behind the reported curve.

Two consequences follow from the construction and are worth stating:

* every point is an **auto search from the model's starting circle**, never the stored
  circle, so `search` is not consulted on this branch — a drawdown's critical surface is not
  the drained one, and it moves with the drawn-down field.
* the stage values on a row are read **on that row's critical surface**, which is the
  drawdown's, so `stage1_FS` varies slightly along the curve. It is the full-pool factor of
  safety of the surface the drawdown settled on, not of the surface a full-pool search would
  have found.

The instant stage 1 itself is read at cannot be a drawdown — a fall from the pool to itself —
and comes back as a `success=False` row saying so. `mode='fem'` is refused: the three-stage
procedure is a limit-equilibrium construction with no SSRM equivalent.

**The instant is never interpolated.** A time that names no saved frame is served by
re-marching with that instant injected into the saved schedule — pass `seep_data=` to
allow it, and the whole set of missing times is served by *one* re-march before the first
solve. Without `seep_data` such a time becomes a `success=False` row saying so. A field
blended between two frames is not a solution of anything, which is why the mode declines
to invent one.

Two settings matter more here than in a value sweep:

* `search=True` (the default) re-searches the critical surface at every instant. That is
  the right default because the critical surface *moves* as the pore pressures change,
  which is the phenomenon being studied.
* the **search window** keeps the curve on one mechanism. `fs_vs_time` reads the circles
  sheet's search limits (entry and exit ranges, center box, tangent depth, minimum slip
  depth) from the model, exactly as Studio's Run LEM path does, so a windowed model gives
  the same family from a script; `search_opts=` overrides it per call and
  `use_file_window=False` ignores it. Without a window, a curve can jump between competing
  minima from one instant to the next, and the jump reads as a change in the slope rather
  than a change in which surface was measured.
* `search_opts={'seed': 'grid'}` — Studio's **Grid search (auto-seed the circular search)**
  box on the same dialog — seeds every instant from a geometry-derived sweep of centers and
  tangent elevations instead of the circles sheet alone, which is what a curve needs when
  the critical mechanism sits in a different part of the section at one end of the run than
  at the other (a dam governed by its downstream face at full pool and by its upstream face
  during a drawdown).

`mode='fem'` runs a full SSRM solve per instant instead (needs a mesh in
`slope_data['mesh']`); `mode='seep'` is refused, because the seepage solution is this
run's input rather than its output.

The per-instant contract is the same as a value sweep's: the base model is preflighted
once before the first solve, and an instant that produces no result is a `success=False`
row **carrying its reason**, never a gap. That matters more for a curve than for a sweep —
a dropped instant moves the curve's minimum to a healthier value with nothing in the
record to say a point went missing.

A worked comparison against a vendor-published curve — both drawdown rates of the
GeoStudio rapid-drawdown example, 11 instants each — is on the
[GeoStudio verification page](../verification/geostudio.md#seepw-t03).

!!! note "Two different questions"
    The default curve and the `rapid=True` curve are different analyses of the same
    physical problem. The default is a sequence of single-stage analyses, one per instant,
    each reading the pore pressures computed for that moment — *how safe is the slope in
    the state it is in at t*. The drawdown curve is a sequence of three-stage
    [rapid drawdowns](../lem/rapid.md), each reading two instants through undrained
    strength envelopes — *how safe is the slope if the pool falls from full to where it
    stands at t*. Neither substitutes for the other, and on a model that carries `d` and ψ
    the drawdown curve is the lower of the two everywhere the undrained envelope bites.

## Sensitivity plots

A sweep across several parameters feeds a family of views. Which one to read depends on the
question: *which parameter matters most* (tornado, scaled bars), *how the answer moves across
the range* (spider), or *where the uncertainty comes from* (variance Pareto, Monte Carlo
rank). Each is a function in `xslope.plot`; all are reachable from the
[Studio **Parametric** dialog](../studio/analysis.md#parametric-study)'s plot-type selector
and from the assistant's `parametric_sweep` call. The
[tornado](#tornado-diagrams) — the most common — is covered in its own section below.

### Tornado diagrams

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

`plot_tornado` draws the base-case FS as a labelled vertical reference line and stacks the
bars widest-on-top by default — the classic Duncan ordering that gives the diagram its name.
Pass `widest_on_top=False` to invert the stack (widest at the bottom); the parameter is kept
for programmatic callers, and Studio deliberately exposes no toggle for it.

When you have already run *full* per-parameter sweeps — a GUI that draws an FS-vs-value curve
per parameter for click-through, for instance — `tornado_from_sweeps()` assembles the same
diagram from those DataFrames with no extra solves, since `plot_tornado` reads each
parameter's lowest- and highest-value FS:

```python
from xslope.sensitivity import sensitivity, tornado_from_sweeps

sweeps = {ref: sensitivity(slope_data, param=ref, rel_range=0.25, n=5,
                           methods=("bishop",))[1]['df']
          for ref in ("mat:Soil:c", "mat:Soil:phi", "mat:Soil:gamma")}
result = tornado_from_sweeps(sweeps, method="bishop")   # {'df', 'base_fs', 'method'}
plot_tornado(result)
```

### Scaled-sensitivity bars

`scaled_sensitivity()` estimates the local derivative $\partial F/\partial p$ for each
parameter with a **central difference at ±1%** (relative) about the base case, then reports
it under three scalings so parameters with different units compare on one axis. `plot_scaled_sensitivity()`
draws one vertical bar per parameter — height is the magnitude, color the sign (green: FS
*rises* with the parameter; red: FS *falls*):

- **`elasticity`** (the default) — the dimensionless $\dfrac{\partial F}{\partial p}\cdot\dfrac{p}{F}$,
  i.e. the percent change in FS per percent change in the parameter.
- **`per_1pct`** — the change in FS for a 1% change in the parameter.
- **`per_sigma`** — the change in FS for a one-σ change (offered only for parameters that
  carry a standard deviation).

```python
from xslope.sensitivity import scaled_sensitivity
from xslope.plot import plot_scaled_sensitivity

ok, result = scaled_sensitivity(slope_data,
                                ["mat:Soil:c", "mat:Soil:phi", "mat:Soil:gamma"],
                                method="bishop")
plot_scaled_sensitivity(result, scaling="elasticity")
```

![Scaled-sensitivity bars](images/parametric_scaled.png){width=800}

For the ACADS sample the friction angle dominates on the elasticity scaling — a percent
change in φ moves FS several times more than a percent change in c or γ. The estimator (the
central-difference derivative and each scaling) is documented in the `scaled_sensitivity`
docstring.

### Spider plot

`plot_spider()` overlays every parameter's FS-vs-value curve on one **normalized** x-axis
(percent change from the base value), with a black base-case marker at the origin and an
FS = 1 guide. The steepness of each line *is* its sensitivity, and a curving line reveals
nonlinearity a single tornado bar would hide. It reuses the per-parameter sweeps you already
ran:

```python
from xslope.sensitivity import sensitivity
from xslope.plot import plot_spider

sweeps = {ref: sensitivity(slope_data, param=ref, rel_range=0.3, n=7,
                           methods=("bishop",))[1]["df"]
          for ref in ("mat:Soil:c", "mat:Soil:phi", "mat:Soil:gamma")}
plot_spider(sweeps)
```

![Spider plot](images/parametric_spider.png){width=800}

### Variance-contribution Pareto

When the model carries standard deviations, `variance_contribution()` asks *which
uncertainties actually drive the scatter in FS*. It reuses the Taylor-series
[reliability](../reliability/taylor.md) machinery — each parameter's variance term is
$\left(\dfrac{\partial F}{\partial p}\,\sigma_p\right)^2$, exactly the $(\Delta F/2)^2$ the
TSPM already computes — normalizes each to a percent of $\mathrm{Var}(F)$, and
`plot_variance_pareto()` draws them descending with a cumulative line:

```python
from xslope.sensitivity import variance_contribution
from xslope.plot import plot_variance_pareto

ok, result = variance_contribution(slope_data, method="bishop")
plot_variance_pareto(result)
```

![Variance-contribution Pareto](images/parametric_variance.png){width=800}

Note how this reorders the parameters versus the elasticity bars: c and φ contribute almost
equally to $\mathrm{Var}(F)$ even though φ has the larger elasticity, because c's standard
deviation is a much larger fraction of its mean. A scaled bar answers *how strong is the
response*; the Pareto answers *how much does this parameter's uncertainty matter* — the
per-σ scaling above is the bridge between the two.

### Monte Carlo rank correlation

`mc_rank_correlation()` runs a Monte Carlo [reliability](../reliability/monte_carlo.md) campaign and computes
the **Spearman rank correlation** between each sampled input and the resulting FS.
`plot_mc_rank_correlation()` draws them signed and sorted:

```python
from xslope.sensitivity import mc_rank_correlation
from xslope.plot import plot_mc_rank_correlation

ok, result = mc_rank_correlation(slope_data, method="bishop", n_samples=10000)
plot_mc_rank_correlation(result)
```

![Monte Carlo rank correlation](images/parametric_mc_rank.png){width=800}

This is a **global** sensitivity measure: unlike the scaled bars — a derivative *at the base
case* — it ranks parameters over the whole sampled distribution, so it captures a
parameter's spread and any nonlinearity in its effect. The local and global measures can
disagree, and the disagreement is itself informative; read the two together.

## Worked sample

The sweep below reproduces the base case of the ACADS simple slope sample
([xslope_acads_simple.xlsx](../lem/files/xslope_acads_simple.xlsx), the same file used
throughout these pages) and brackets it over ±50% of the cohesion, re-searching the
critical surface at every point. The regression suite locks the end points: at c = 1.5
the critical circle has retreated toward the face (the surface-jump columns show the
radius growing as cohesion falls), and at c = 4.5 the slope crosses FS = 1.

<!-- test: file=../lem/files/xslope_acads_simple.xlsx, type=sensitivity, param=mat:Soil:c, method=bishop, num_slices=30, n=3, rel_range=0.5, expected_base=0.985, expected_low=0.882, expected_high=1.073, tolerance=0.01, benchmark=SENS-acads-c -->

| c (kPa) | Bishop FS |
|---|---|
| 1.5 (−50%) | 0.882 |
| 3.0 (base) | 0.985 |
| 4.5 (+50%) | 1.073 |

A `geom:piezo:dy` sweep on any of the water-table samples answers the most common
geometry question — "how sensitive is this slope to the water table?" — without writing
a setter: the reference shifts the piezometric line vertically by the swept delta.

For a published benchmark of the sweep itself, see
[verification problem VP40](../verification/rocscience.md#vp40) — [Perry (1993)](https://doi.org/10.1144/GSL.QJEGH.1994.027.P3.04)'s
power-curve sensitivity study, where the swept ΔFS matches Slide's published curves
within about a percent at every range endpoint.
