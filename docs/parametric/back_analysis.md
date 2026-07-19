# Back-Analysis Mode

Back-analysis is the same single-parameter sweep as [Design](design.md), inverted for a *failure
investigation*. A slide has already occurred, so the factor of safety at the moment of
failure is known to be exactly 1.0; the unknown is a strength (or pore-pressure, or loading)
parameter, and the study back-calculates the value **consistent with the observed failure** —
the mobilized shear strength implied by the slide, most commonly. `back_analysis()` is
`design()` with `target_fs` defaulting to 1.0 and the result worded for the forensic reading.
It takes the same parameter grammar (see
[Addressing a parameter](index.md#addressing-a-parameter)) and the same `mode=`/`fem_opts=`/
`seep_opts=` engine controls as `design()` — in `mode='seep'` the back-calculated quantity is
whatever discharge q is consistent with the observed condition, rather than a strength value.

```python
from xslope.sensitivity import back_analysis

success, result = back_analysis(
    slope_data,
    param="mat:Soil:c",          # the strength parameter to back-calculate
    low=1.0, high=6.0, steps=11, # sweep the plausible range...
    # target_fs=1.0,             # ...to the known failure condition (the default)
    method="bishop",
)
print(result['message'])
# Back-analysis: mat:Soil:c = 3.25 gives FS = 1 (the value consistent with the
# observed failure).
```

`result['crossing']` is the back-calculated value; `result['study']` is `'back_analysis'`;
every other field carries the same meaning as [`design()`](design.md#running-a-design-study).
The same [*never-extrapolate* discipline](design.md#honest-about-misses) applies — if the
swept range never reaches FS = 1.0, `bracketed` is `False` and `extend` says which way to
widen it, rather than guessing a value past the last solve.
