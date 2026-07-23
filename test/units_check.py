"""Unit checks for xslope.units — the unit-system single source of truth.

Phase 1 of the units plan is "declare + label, never convert", and its hard
contract is *fix-first with no behavior change*: gamma_water's canonical values now
live in exactly one place, and every scattered silent fallback (9.81 / 62.4 / 9.807)
is replaced by a loud, no-default read. This file locks that contract at the unit
level:

  1. GAMMA_W carries the two canonical constants and nothing drifts.
  2. require_gamma_water returns the model's value when present and raises loudly —
     never substituting a canonical constant — when it is genuinely missing. That
     loud read is the mechanism every solver/export/plot site now shares.
  3. infer_unit_system reproduces the elastic_props gamma-magnitude heuristic
     (max gamma > 50 -> imperial) and returns None when there is no signal.
  4. labels() honors the plan's key contract, including empty strings for an
     undeclared (None) system and for time-bearing keys with no declared time unit.
  5. The two heaviest consumers — build_seep_data and build_fem_data — actually fire
     the loud behavior when a model omits gamma_water (no silent 9.81/62.4 that could
     flip the unit system mid-pipeline).

Run directly:  PYTHONPATH=. python3 test/units_check.py
"""

import numpy as np

from xslope import units
from xslope.units import GAMMA_W, require_gamma_water, infer_unit_system, labels


def _expect_raises(fn, needle, label):
    """Return [] if fn() raises with `needle` in the message, else a one-item failure."""
    try:
        fn()
    except Exception as e:  # noqa: BLE001 — the point is to catch whatever it raises
        if needle in str(e):
            return []
        return [f"{label}: raised but message lacked {needle!r}: {e!r}"]
    return [f"{label}: expected a raise mentioning {needle!r}, none occurred"]


def check_constants():
    fails = []
    if GAMMA_W != {"si": 9.81, "imperial": 62.4}:
        fails.append(f"GAMMA_W drifted: {GAMMA_W!r}")
    return fails


def check_require_gamma_water():
    fails = []
    # Present -> returned as a float, unchanged.
    for val in (9.81, 62.4, 10.05, 64):
        got = require_gamma_water({"gamma_water": val})
        if got != float(val):
            fails.append(f"require_gamma_water returned {got!r} for {val!r}")
    # Works for any mapping key name style (fem_data / water blocks share the key).
    if require_gamma_water({"gamma_water": 9.81, "extra": 1}) != 9.81:
        fails.append("require_gamma_water mishandled a richer mapping")
    # Missing key -> raise; None value -> raise; both mention gamma_water.
    fails += _expect_raises(lambda: require_gamma_water({}),
                            "gamma_water", "require missing-key")
    fails += _expect_raises(lambda: require_gamma_water({"gamma_water": None}),
                            "gamma_water", "require None-value")
    # The context phrase is woven into the message.
    fails += _expect_raises(lambda: require_gamma_water({}, "seepage analysis"),
                            "seepage analysis", "require context")
    return fails


def check_infer_unit_system():
    fails = []
    cases = [
        ([{"gamma": 18.0}, {"gamma": 20.0}], "si"),          # metric band
        ([{"gamma": 9.81}], "si"),
        ([{"gamma": 120.0}], "imperial"),                     # imperial band
        ([{"gamma": 18.0}, {"gamma": 125.0}], "imperial"),   # max wins
        ([], None),                                           # no materials
        ([{"gamma": 0.0}], None),                             # no positive gamma
        ([{"gamma": None}, {"c": 5}], None),                  # unset gammas
        ([{"gamma": float("nan")}], None),                    # NaN is not a signal
        ([{"gamma": 50.0}], "si"),                            # threshold is strict >
        ([{"gamma": 50.001}], "imperial"),
    ]
    for materials, expected in cases:
        got = infer_unit_system(materials)
        if got != expected:
            fails.append(f"infer_unit_system({materials!r}) = {got!r}, expected {expected!r}")
    return fails


def check_labels():
    fails = []
    keys = {"length", "stress", "unit_weight", "force_per_len", "k", "flowrate", "time"}

    # None -> every key present and empty (legacy/undeclared degrades to bare labels).
    none_lbl = labels(None)
    if set(none_lbl) != keys:
        fails.append(f"labels(None) keys = {set(none_lbl)!r}")
    if any(v != "" for v in none_lbl.values()):
        fails.append(f"labels(None) not all empty: {none_lbl!r}")

    # SI, no time unit: length/stress/unit_weight/force_per_len set; time-bearing empty.
    si = labels("si")
    for k, want in (("length", "m"), ("stress", "kPa"),
                    ("unit_weight", "kN/m³"), ("force_per_len", "kN/m")):
        if si[k] != want:
            fails.append(f"labels('si')[{k!r}] = {si[k]!r}, expected {want!r}")
    for k in ("k", "flowrate", "time"):
        if si[k] != "":
            fails.append(f"labels('si')[{k!r}] should be empty w/o time_unit, got {si[k]!r}")

    # Imperial, no time unit.
    imp = labels("imperial")
    for k, want in (("length", "ft"), ("stress", "psf"),
                    ("unit_weight", "pcf"), ("force_per_len", "lb/ft")):
        if imp[k] != want:
            fails.append(f"labels('imperial')[{k!r}] = {imp[k]!r}, expected {want!r}")

    # With a declared time unit the time-bearing keys fill in.
    si_day = labels("si", "day")
    for k, want in (("time", "day"), ("k", "m/day"), ("flowrate", "m³/day per m")):
        if si_day[k] != want:
            fails.append(f"labels('si','day')[{k!r}] = {si_day[k]!r}, expected {want!r}")
    imp_s = labels("imperial", "s")
    for k, want in (("time", "s"), ("k", "ft/s"), ("flowrate", "ft³/s per ft")):
        if imp_s[k] != want:
            fails.append(f"labels('imperial','s')[{k!r}] = {imp_s[k]!r}, expected {want!r}")

    # Case-insensitive system token.
    if labels("SI") != labels("si") or labels("Imperial") != labels("imperial"):
        fails.append("labels() is not case-insensitive on unit_system")

    # None system ignores time_unit entirely (stays all-empty).
    if any(v != "" for v in labels(None, "day").values()):
        fails.append("labels(None, 'day') should stay all-empty")

    # Unknown system is a loud error, not a silent guess.
    fails += _expect_raises(lambda: labels("metric"), "Unknown unit system",
                            "labels bogus system")
    return fails


def _tiny_mesh():
    """A one-triangle mesh — enough for build_seep_data / the FEM piezo path to reach
    their gamma_water read without any heavy assembly."""
    return {
        "nodes": np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]]),
        "elements": np.array([[0, 1, 2]]),
        "element_types": np.array([3]),
        "element_materials": np.array([1]),
    }


def check_seep_entry_point_raises():
    """build_seep_data must fire the loud read (not silently assume 9.81) when a
    model omits gamma_water."""
    from xslope.seep import build_seep_data
    slope_data = {"materials": [{"name": "m1", "k1": 1.0, "k2": 1.0}]}  # no gamma_water
    return _expect_raises(lambda: build_seep_data(_tiny_mesh(), slope_data),
                          "gamma_water", "build_seep_data missing-gamma")


def check_fem_entry_point_raises():
    """build_fem_data must fire the loud read on the piezo pore-pressure path when a
    model omits gamma_water (previously a silent 62.4/9.81 that could flip units)."""
    from xslope.fem import build_fem_data
    slope_data = {
        "materials": [{"name": "m1", "u": "piezo", "c": 10.0, "phi": 20.0,
                       "gamma": 18.0, "E": 1.0e5, "nu": 0.3}],
        "piezo_line": [(0.0, 0.5), (1.0, 0.5)],
        # deliberately no gamma_water
    }
    return _expect_raises(lambda: build_fem_data(slope_data, mesh=_tiny_mesh()),
                          "gamma_water", "build_fem_data missing-gamma")


def main():
    print("xslope.units unit checks:")
    checks = [
        ("constants", check_constants),
        ("require_gamma_water", check_require_gamma_water),
        ("infer_unit_system", check_infer_unit_system),
        ("labels", check_labels),
        ("seep entry point", check_seep_entry_point_raises),
        ("fem entry point", check_fem_entry_point_raises),
    ]
    failures = []
    for name, fn in checks:
        fs = fn()
        status = "ok" if not fs else f"FAIL ({len(fs)})"
        print(f"  {name:22s} {status}")
        failures += fs
    if failures:
        print("\nFAILURES:")
        for f in failures:
            print(f"  - {f}")
        raise SystemExit(1)
    print("\nAll xslope.units checks passed.")


if __name__ == "__main__":
    main()
