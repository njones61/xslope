"""Where the automatic water load is derived, and how often (``xslope.search``).

With automatic water loads on (``main!D23`` = auto) the weight of the standing
water is not something the user typed: it is measured, at solve time, from the
model's own water definition against the ground surface. The slicer does that
measurement defensively at the top of every ``generate_slices`` call -- which is
correct, and which under a search means once per TRIAL SURFACE. A circular
search slices several thousand trials to report one factor of safety, and every
one of those derivations returns the same blocks, because a search varies the
circle and nothing else. Neither the water line nor the ground surface is
something it touches.

So the derivation is hoisted to the entry of the search and the derived model is
what every trial is sliced against. That is a performance change with a
correctness claim attached, and this module pins both halves of it:

  * **once per search, not once per trial.** The count is asserted against a run
    that really did slice hundreds of surfaces, so a regression to per-trial
    derivation cannot hide behind a small search.
  * **the same answer.** A model sliced after the hoist must produce the identical
    slice table, load for load and weight for weight, as the same model sliced
    before it -- and the identical factor of safety.
  * **nothing is cached between runs.** The reuse is scoped to one search, over a
    model that cannot change while it runs. A model whose water definition has
    been edited -- which is every second in Studio -- must be re-derived, and the
    load it gets must be the new one. A memo keyed on the model would answer the
    second run with the first run's reservoir, and that is the failure this
    check exists to make loud.
  * **manual mode is untouched.** A model that supplies its own water loads
    derives nothing, at either place.

File-light: it loads two shipped samples and runs short searches
(``max_iter`` clipped, few slices) -- the counts and the identities are what is
under test, not the criticals.

Run directly:  PYTHONPATH=. python3 test/water_hoist_check.py
"""

import warnings
from pathlib import Path

import numpy as np

warnings.filterwarnings('ignore')

from xslope import water as WATER
from xslope.fileio import load_slope_data
from xslope.slice import generate_slices

ROOT = Path(__file__).resolve().parent.parent / 'docs'
#: An embankment against a reservoir: automatic water loads, real derived blocks.
AUTO = ROOT / 'lem' / 'files' / 'xslope_earth_dam_up.xlsx'
#: An automatic-water model that also carries a non-circular starting surface.
AUTO_NONCIRC = ROOT / 'verification' / 'files' / 'rocscience' / 'vp034.xlsx'


class _Counter:
    """Counts (and passes through) every call to the derivation itself.

    ``with_water_loads`` reaches ``derive_water_loads`` as a module global, so
    replacing the module attribute intercepts the real derivation wherever it is
    reached from -- the slicer's defensive call included. That is the point: a
    reintroduced per-trial derivation would be counted here no matter which
    function made it.
    """

    def __init__(self):
        self.n = 0
        self._orig = WATER.derive_water_loads

    def __enter__(self):
        def counted(*a, **kw):
            self.n += 1
            return self._orig(*a, **kw)
        WATER.derive_water_loads = counted
        return self

    def __exit__(self, *exc):
        WATER.derive_water_loads = self._orig
        return False


def _auto_model():
    sd = load_slope_data(str(AUTO))
    return sd


def _short_circular(sd, **kw):
    from xslope.search import circular_search
    kw.setdefault('num_slices', 12)
    kw.setdefault('max_iter', 3)
    return circular_search(sd, 'bishop', **kw)


#: The most derivations a whole search may make. A hoisted search makes one pass
#: over the model (one call per stage) and the entry preflight makes its own,
#: reading the same water definition for the double-count rule -- a small fixed
#: budget either way. What matters is that it is FIXED: per-trial derivation runs
#: to two per trial, which is hundreds.
MAX_PER_SEARCH = 6


# ---------------------------------------------------------------------------

def check_circular_search_derives_a_fixed_number_of_times(failures):
    """A search's derivation count is a small constant, not a count of trials."""
    sd = _auto_model()
    if WATER.water_loads_mode(sd) != 'auto':
        failures.append(f"{AUTO.name} is not an automatic-water model; this check "
                        f"would pass vacuously")
        return
    seen = []
    for max_iter in (2, 6):
        with _Counter() as c:
            fs_cache, _conv, _path, circle_cache = _short_circular(
                _auto_model(), max_iter=max_iter)
        seen.append((max_iter, c.n, len(circle_cache), bool(fs_cache)))

    for max_iter, n, trials, got in seen:
        if not got:
            failures.append(f"max_iter={max_iter}: the search produced no surfaces, "
                            f"so nothing was sliced")
        if trials < 100:
            failures.append(f"max_iter={max_iter}: only {trials} trials -- too few "
                            f"for the count to mean anything")
        if n > MAX_PER_SEARCH:
            failures.append(f"max_iter={max_iter}: a search over {trials} trials "
                            f"derived the water load {n} times; a hoisted search "
                            f"derives a fixed few (<= {MAX_PER_SEARCH}), and a count "
                            f"near {2 * trials} means it is back to once per trial")
    if seen[0][1] != seen[1][1]:
        failures.append(f"the derivation count tracks the search's length: "
                        f"{seen[0][1]} at {seen[0][2]} trials vs {seen[1][1]} at "
                        f"{seen[1][2]}. It must not depend on how long the search runs")


def check_noncircular_search_derives_a_fixed_number_of_times(failures):
    """The same contract on the non-circular search."""
    from xslope.search import noncircular_search

    sd = load_slope_data(str(AUTO_NONCIRC))
    if WATER.water_loads_mode(sd) != 'auto' or not sd.get('non_circ'):
        failures.append(f"{AUTO_NONCIRC.name} no longer carries automatic water "
                        f"loads and a non-circular surface; this check would pass "
                        f"vacuously")
        return
    with _Counter() as c:
        try:
            fs_cache, _conv, path = noncircular_search(
                sd, 'janbu', num_slices=15, max_iter=3, diagnostic=False)
        except Exception as e:                                 # noqa: BLE001
            failures.append(f"the non-circular search failed: {e}")
            return
    if c.n > MAX_PER_SEARCH:
        failures.append(f"a non-circular search over {len(path or [])} steps derived "
                        f"the water load {c.n} times; a hoisted search derives a "
                        f"fixed few (<= {MAX_PER_SEARCH})")
    if not fs_cache:
        failures.append("the non-circular search produced no surfaces")


def check_hoisting_changes_nothing(failures):
    """The slice table and the factor of safety are the same, to the last digit,
    whether the model was derived by the slicer or handed in already derived."""
    from xslope import solve

    sd = _auto_model()
    circle = (sd.get('circles') or [None])[0]
    if circle is None:
        failures.append(f"{AUTO.name} carries no circle to slice")
        return
    ok_a, res_a = generate_slices(sd, circle=circle, num_slices=30,
                                  check_inputs=False)
    ok_b, res_b = generate_slices(WATER.with_water_loads(sd), circle=circle,
                                  num_slices=30, check_inputs=False)
    if not (ok_a and ok_b):
        failures.append(f"slicing failed: raw={res_a if not ok_a else 'ok'}, "
                        f"hoisted={res_b if not ok_b else 'ok'}")
        return
    df_a, df_b = res_a[0], res_b[0]
    if list(df_a.columns) != list(df_b.columns) or len(df_a) != len(df_b):
        failures.append("the hoisted model produced a differently shaped slice table")
        return
    for col in df_a.columns:
        try:
            a = np.asarray(df_a[col].values, dtype=float)
            b = np.asarray(df_b[col].values, dtype=float)
        except (TypeError, ValueError):
            continue
        if not np.array_equal(a, b, equal_nan=True):
            worst = float(np.nanmax(np.abs(a - b)))
            failures.append(f"column '{col}' differs after the hoist "
                            f"(largest difference {worst:g}); the derived load must "
                            f"be bit-identical, not merely close")
    ok1, r1 = solve.bishop(df_a)
    ok2, r2 = solve.bishop(df_b)
    if not (ok1 and ok2):
        failures.append("bishop refused one of the two slice tables")
    elif repr(r1['FS']) != repr(r2['FS']):
        failures.append(f"FS moved with the hoist: {r1['FS']!r} vs {r2['FS']!r}")


def check_a_changed_water_state_is_re_derived(failures):
    """Nothing survives the call. Change what the water weighs and the next
    derivation must see it -- a memo keyed on the model would serve the old one,
    which is what Studio's live editing would hit first."""
    sd = _auto_model()
    first = WATER.with_water_loads(sd)
    r1 = (first.get('water_derived') or {}).get(1, {}).get('resultant')

    heavier = dict(sd)
    heavier['gamma_water'] = float(sd['gamma_water']) * 2.0
    # a fresh model, not one handed back its own earlier derivation
    for k in ('dloads_derived', 'dloads2_derived', 'water_derived'):
        heavier.pop(k, None)

    with _Counter() as c:
        second = WATER.with_water_loads(heavier)
    r2 = (second.get('water_derived') or {}).get(1, {}).get('resultant')
    if c.n == 0:
        failures.append("the edited model derived nothing at all -- it was served "
                        "a cached answer")
    if r1 is None or r2 is None:
        failures.append(f"no resultant was reported: before={r1!r}, after={r2!r}")
    elif abs(r2 - 2.0 * r1) > 1e-6 * abs(r1):
        failures.append(f"doubling the unit weight of water took the derived load "
                        f"from {r1!r} to {r2!r}, not to twice the first -- the "
                        f"second reading is not of the edited model")


def check_manual_mode_derives_nothing(failures):
    """A model that supplies its own water loads derives nothing, anywhere."""
    sd = _auto_model()
    sd['water_loads'] = 'manual'
    with _Counter() as c:
        fs_cache, _conv, _path, _cc = _short_circular(sd)
    if c.n != 0:
        failures.append(f"a manual-mode search derived the water load {c.n} times; "
                        f"manual mode must derive nothing at all")
    if not fs_cache:
        failures.append("the manual-mode search produced no surfaces")


CHECKS = [
    check_circular_search_derives_a_fixed_number_of_times,
    check_noncircular_search_derives_a_fixed_number_of_times,
    check_hoisting_changes_nothing,
    check_a_changed_water_state_is_re_derived,
    check_manual_mode_derives_nothing,
]


def run():
    """Returns a list of failure descriptions (empty on success)."""
    failures = []
    for chk in CHECKS:
        try:
            chk(failures)
        except Exception as e:                                 # noqa: BLE001
            failures.append(f"{chk.__name__} raised {type(e).__name__}: {e}")
    return failures


if __name__ == '__main__':
    import sys
    bad = run()
    for f in bad:
        print(f"FAIL  {f}")
    print("water_hoist_check: "
          + (f"{len(bad)} failure(s)" if bad else f"{len(CHECKS)} checks passed"))
    sys.exit(1 if bad else 0)
