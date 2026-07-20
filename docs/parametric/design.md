# Design Mode

Where a sweep asks *how much does FS move*, a **design study** asks the inverse — *what value
of this parameter gives FS = 1.5?* `design()` runs a fixed number of evenly spaced solves
across an explicit `[low, high]` range and linearly interpolates the parameter value where the
FS curve crosses the target. It is the deterministic-design staple: "vary the undrained
strength between X and Y and find where FS reaches the design factor." It shares the same
parameter grammar as [Sensitivity](sensitivity.md) — see
[Addressing a parameter](index.md#addressing-a-parameter) — and the same three
[engine modes](index.md#analysis-modes-lem-fem-seepage): `mode='lem'` (the default) and
`mode='fem'` both target a factor of safety, while `mode='seep'` sweeps the seepage problem
and targets the **total discharge q** instead — `target_fs` and the crossing/bracketing logic
below apply identically, just against q rather than FS.

In XSLOPE Studio this runs behind the **Parametric…** button (its Design mode) — see
[Studio: Parametric study](../studio/analysis.md#parametric-study) for the dialog.

## Running a design study

```python
from xslope.sensitivity import design
from xslope.plot import plot_sensitivity

success, result = design(
    slope_data,
    param="mat:Soil:c",          # or {"material": "Soil", "property": "c"}, or a tuple
    low=6, high=18, steps=7,     # 7 evenly spaced solves across [6, 18]
    target_fs=1.5,
    method="bishop",             # one method — a design study locates a single curve
    num_slices=30,
)
print(result['message'])
plot_sensitivity(result['df'], target_fs=result['target_fs'])
```

`design()` returns `(success, result)`; on success `result` carries the sweep DataFrame
(`result['df']`, exactly the sensitivity shape above) plus a summary:

| field | meaning |
|---|---|
| `crossing` | interpolated parameter value at FS = `target_fs`, or `None` if the target is not reached |
| `crossings` | every crossing found — a non-monotonic curve can cross twice; `crossing` is the first |
| `bracketed` | `True` only when the target is crossed *inside* `[low, high]` |
| `fs_range` | `(min FS, max FS)` over the successful sweep points |
| `direction` | `'increasing'` / `'decreasing'` / `'non-monotonic'` trend of FS vs the parameter |
| `extend` | on a miss, `'above {high}'` or `'below {low}'` — which way to widen the range; else `None` |
| `message` | one-line human-readable summary |
| `param`, `target_fs`, `base_value`, `runtime` | canonical ref, the target, the parameter's current value, wall-clock seconds |

`design()` also takes `progress_callback` and `cancel_check` hooks (a per-point progress
callback and a cooperative cancel), which is how Studio streams a progress bar and a Cancel
button over a background sweep; a plain data-in/data-out call leaves both `None`.

## Sweeping anything else: `modify=`

`design()` shares the sensitivity engine's escape hatch, so the quantity you design need not
be a single stored scalar. When it is **geometry** — a slope angle, a berm width — pass a
callable in place of a parameter reference: `modify` is a `(slope_data, value) -> slope_data`
setter and `label` names the swept axis, and the two are mutually exclusive with `param`
(exactly one per call). `low`/`high`/`steps` sweep the callable's value the same way, and the
crossing, bracketing, and `extend`/`direction` reporting are identical. The callable owns
whatever consistency its edit requires — here, rebuilding the material polygons and ground
surface after it moves a profile point, since slice weights come from the polygons, not the
raw profile lines:

```python
import math
from shapely.geometry import Polygon
from xslope.fileio import build_ground_surface_from_polygons
from xslope.mesh import build_polygons

def set_slope_angle(sd, beta_deg):
    """Rotate the slope face about the toe to beta_deg, then rebuild the polygon
    geometry (see the Slope Design page for this rebuild pattern)."""
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

success, result = design(
    slope_data, modify=set_slope_angle, label="slope angle (deg)",
    low=12, high=25, steps=8,        # flatten the face from 25° toward 12°...
    target_fs=1.5, method="bishop",  # ...and find the angle that reaches FS = 1.5
)
print(result['message'])
# FS = 1.5 at slope angle (deg) = 16.74 (interpolated between solves).
```

On the shipped ACADS slope this brackets the target — FS climbs from about 1.04 at 25° to
2.02 at 12° as the face flattens, crossing FS = 1.5 at a slope angle of **16.74°**. A
built-in `param` reference and a `modify=` callable are **one code path** — every reference
resolves internally to exactly the setter signature `modify=` takes — so the answer is
identical whichever way you name the swept axis, and the engine validates the modified model
at each step (polygon validity, ground surface present) precisely because a setter may be
user-written: a broken edit becomes an honest `success=False` sweep point, never a silently
inconsistent crossing.

## Honest about misses

**The engine never extrapolates.** A crossing is reported *only* when
the target is bracketed by two actual solves. If the swept range never reaches the target,
`bracketed` is `False`, `crossing` is `None`, and `extend` names the direction to widen the
range — the study reports that it fell short rather than projecting a value past the last
solve.

## Worked example

On the [ACADS simple slope](../lem/files/xslope_acads_simple.xlsx) used throughout
these pages (base cohesion c = 3 kPa, base Bishop FS = 0.985), sweeping the cohesion from 6 to
18 kPa in seven steps and asking where FS reaches 1.5:

```
FS = 1.5 at mat:Soil:c = 13.4 (interpolated between solves).
```

`result['crossing']` is 13.40, `result['bracketed']` is `True`, and `result['direction']` is
`'increasing'`; the sweep spans FS = 1.154 (c = 6) to FS = 1.693 (c = 18). Re-solving at the
interpolated c = 13.40 returns FS = 1.500 — the linear interpolation between the c = 12 and
c = 14 solves lands the target to within 0.01%.

Ask the same question over a range that never reaches the target — say c from 3 to 9 — and the
study declines to guess:

```
FS = 1.5 is not reached for mat:Soil:c in [3, 9] — FS spans [0.985, 1.302].
Extend the range above 9 to bracket it.
```

Here `result['bracketed']` is `False`, `result['crossing']` is `None`, and `result['extend']`
is `'above 9'`.
