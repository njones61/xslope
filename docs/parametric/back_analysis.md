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

In XSLOPE Studio this runs behind the **Parametric…** button (its Back-Analysis mode) — see
[Studio: Parametric study](../studio/analysis.md#parametric-study) for the dialog.

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

## Inherited `modify=` callable

Because `back_analysis()` is `design()` with a forensic target, it inherits the same
[`modify=` escape hatch](design.md#sweeping-anything-else-modify) unchanged: pass a
`(slope_data, value) -> slope_data` callable and a `label` in place of a `param` reference
(exactly one per call) whenever the unknown is not a single stored scalar. A common forensic
use is a **water-table elevation** — the phreatic surface at the moment of a slide is rarely
recorded, so sweep it and back-calculate the level consistent with FS = 1.0:

```python
def set_water_table(sd, elev):
    """Set the piezometric line to a horizontal elevation `elev`."""
    sd['piezo_line'] = [(x, elev) for x, _ in sd['piezo_line']]
    return sd

success, result = back_analysis(
    slope_data, modify=set_water_table, label="water-table elevation (m)",
    low=8.0, high=16.0, steps=9,   # sweep the plausible phreatic range...
    # target_fs=1.0,               # ...to the known failure condition (the default)
    method="bishop",
)
print(result['message'])
# Back-analysis: water-table elevation (m) = 13.4 gives FS = 1 (the value
# consistent with the observed failure).
```

`result['crossing']` is the back-calculated elevation; every other field carries the same
meaning as [`design()`](design.md#running-a-design-study), and the same never-extrapolate
discipline applies whichever way the swept axis is named.
