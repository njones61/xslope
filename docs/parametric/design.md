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
