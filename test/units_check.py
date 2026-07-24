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

Phase 2 adds the DECLARATION plumbing (unit_system / time_unit on slope_data) and
locks its contract here too:

  6. infer_system_from_gamma_water bands: ~9.81 AND GeoStudio's 9.807 -> SI, ~62.4 ->
     Imperial, everything else (incl. a seawater override) -> None.
  7. normalize_unit_system folds every importer vocabulary (metric/si/imperial/unknown)
     onto si/imperial/None -- the one-liner each import site persists through.
  8. units_check is a warnings-only cross-check: silent for an undeclared (None) model;
     fires on a divergent gamma_water, materials that look like the other system, and an
     E magnitude out of the declared system's band; silent when all three agree.
  9. Loading a real <= v17 template infers a unit_system but NEVER a time_unit, and a
     real .s2d import infers its unit_system from the file's bare gamma_water.

Run directly:  PYTHONPATH=. python3 test/units_check.py
"""

import os

import numpy as np

from xslope import units
from xslope.units import (GAMMA_W, require_gamma_water, infer_unit_system, labels,
                          infer_system_from_gamma_water, normalize_unit_system,
                          units_check)

_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


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


def check_infer_system_from_gamma_water():
    fails = []
    cases = [
        (9.81, "si"),        # xslope SI canonical
        (9.807, "si"),       # GeoStudio vendor canonical -- must land in the SI band
        (9.80, "si"),
        (9.72, "si"),        # just inside the 0.15 band
        (62.4, "imperial"),  # Imperial canonical
        (61.5, "imperial"),  # just inside the 1.0 band
        (63.3, "imperial"),
        (10.05, None),       # seawater override -> off-band, unlabeled (not mislabeled)
        (64.0, None),        # seawater Imperial -> off-band
        (5.0, None),         # nonsense magnitude
        (30.0, None),        # between the bands
        (None, None),        # no value
        (float("nan"), None),
        (float("inf"), None),
        ("not a number", None),
    ]
    for g, expected in cases:
        got = infer_system_from_gamma_water(g)
        if got != expected:
            fails.append(f"infer_system_from_gamma_water({g!r}) = {got!r}, expected {expected!r}")
    return fails


def check_normalize_unit_system():
    fails = []
    # Every vocabulary each importer speaks -> the three slope_data values.
    cases = [
        ("si", "si"), ("SI", "si"), ("Si", "si"),
        ("metric", "si"), ("Metric", "si"), ("  metric ", "si"),  # gsz + slide2
        ("imperial", "imperial"), ("Imperial", "imperial"),       # gsz + rs2 + slide2
        ("unknown", None),                                        # rs2 no-tag path
        ("", None), ("english", None), (None, None),
    ]
    for token, expected in cases:
        got = normalize_unit_system(token)
        if got != expected:
            fails.append(f"normalize_unit_system({token!r}) = {got!r}, expected {expected!r}")
    return fails


def check_units_check():
    fails = []

    # Undeclared (None) system -> completely silent (a legacy model is unchanged).
    for none_model in ({"unit_system": None, "gamma_water": 999.0,
                        "materials": [{"gamma": 999.0, "E": 1e12}]},
                       {"gamma_water": 62.4}):  # key absent entirely
        w = units_check(none_model)
        if w != []:
            fails.append(f"units_check on undeclared model should be silent, got {w!r}")

    def _fires(model, needle, label):
        w = units_check(model)
        if not any(needle in m for m in w):
            fails.append(f"{label}: expected a warning containing {needle!r}, got {w!r}")

    def _silent_on(model, needle, label):
        w = units_check(model)
        if any(needle in m for m in w):
            fails.append(f"{label}: expected NO warning containing {needle!r}, got {w!r}")

    # 1. gamma_water divergence (>2%): 11.0 vs SI 9.81 is ~12% off.
    _fires({"unit_system": "si", "gamma_water": 11.0}, "gamma_water",
           "gamma_water divergent")
    # canonical value is silent.
    _silent_on({"unit_system": "si", "gamma_water": 9.81}, "gamma_water",
               "gamma_water canonical")
    # Seawater override (~2.4% off) is legal but DOES warn (names both values).
    _fires({"unit_system": "si", "gamma_water": 10.05}, "gamma_water",
           "gamma_water seawater warns-but-legal")

    # 2. Materials look like the other system.
    _fires({"unit_system": "si", "materials": [{"gamma": 120.0}, {"gamma": 125.0}]},
           "look", "materials mismatch")
    # Consistent materials are silent on the mismatch check.
    _silent_on({"unit_system": "si", "materials": [{"gamma": 18.0}, {"gamma": 20.0}]},
               "look", "materials consistent")

    # 3. E magnitude out of band.
    _fires({"unit_system": "si", "materials": [{"gamma": 18.0, "E": 5.0e8}]},
           "Young's modulus", "E too large for SI")
    _fires({"unit_system": "imperial", "materials": [{"gamma": 120.0, "E": 8000.0}]},
           "Young's modulus", "E too small for Imperial")
    # In-band E is silent.
    _silent_on({"unit_system": "si", "materials": [{"gamma": 18.0, "E": 1.0e5}]},
               "Young's modulus", "E in SI band")
    _silent_on({"unit_system": "imperial", "materials": [{"gamma": 120.0, "E": 2.0e6}]},
               "Young's modulus", "E in Imperial band")

    # A fully-consistent declared model is entirely silent.
    good = {"unit_system": "si", "gamma_water": 9.81,
            "materials": [{"gamma": 18.0, "E": 1.0e5}, {"gamma": 20.0, "E": 1.5e5}]}
    if units_check(good) != []:
        fails.append(f"units_check on a consistent SI model should be silent, "
                     f"got {units_check(good)!r}")

    # normalize is applied inside: an importer's "metric" declaration is checked as SI.
    _fires({"unit_system": "metric", "gamma_water": 11.0}, "gamma_water",
           "units_check normalizes 'metric'")
    return fails


def check_loader_time_never_inferred():
    """A real <= v17 template infers a unit_system from its gamma_water but NEVER a
    time_unit (plan law: time is never inferred for template files)."""
    from xslope.fileio import load_slope_data
    fails = []
    path = os.path.join(_REPO_ROOT, "docs", "lem", "files", "xslope_acads_simple.xlsx")
    if not os.path.exists(path):
        return [f"loader fixture missing: {path}"]
    import warnings as _w
    with _w.catch_warnings():
        _w.simplefilter("ignore")  # unit-label notes are expected; we assert on values
        sd = load_slope_data(path)
    if "unit_system" not in sd or "time_unit" not in sd:
        fails.append("loader did not add unit_system/time_unit keys")
        return fails
    if sd["time_unit"] is not None:
        fails.append(f"loader inferred a time_unit ({sd['time_unit']!r}); it must be None "
                     f"for a template file")
    # Whatever it inferred must equal the gamma_water-band inference (labels only).
    expected = infer_system_from_gamma_water(sd.get("gamma_water"))
    if sd["unit_system"] != expected:
        fails.append(f"loader unit_system {sd['unit_system']!r} != gamma_water inference "
                     f"{expected!r}")
    if sd["unit_system"] not in ("si", "imperial"):
        fails.append(f"expected a real template to infer si/imperial, got "
                     f"{sd['unit_system']!r} (gamma_water={sd.get('gamma_water')!r})")
    return fails


def check_s2d_import_unit_system():
    """A real .s2d import infers its unit_system from the file's bare gamma_water."""
    from xslope.seep import import_seep2d
    fails = []
    path = os.path.join(_REPO_ROOT, "docs", "inputs", "seep", "seep2d", "lface.s2d")
    if not os.path.exists(path):
        return [f"s2d fixture missing: {path}"]
    mesh = import_seep2d(path)
    if "unit_system" not in mesh:
        fails.append("import_seep2d did not add a unit_system key")
        return fails
    expected = infer_system_from_gamma_water(mesh.get("unit_weight"))
    if mesh["unit_system"] != expected:
        fails.append(f"s2d unit_system {mesh['unit_system']!r} != gamma_water inference "
                     f"{expected!r} (unit_weight={mesh.get('unit_weight')!r})")
    return fails


def check_importer_persistence_mapping():
    """Persistence for the gsz/rs2/slide2 importers, at the function level.

    No cheap real .gsz/.fez/.sli fixture exists in the repo (the gsz suite authors a
    synthetic file; rs2/slide2 corpora live outside it), so per the plan we unit-test the
    persistence FUNCTION against the exact detection vocabulary each importer feeds it.
    The end-to-end import path is exercised by run_tests.py --gsz / --rs2 / --slide2."""
    fails = []
    # geostudio._detect_units emits "metric"/"imperial" (raises otherwise) -> si/imperial.
    # rs2._detect_units emits "metric"/"imperial"/"unknown" -> si/imperial/None.
    # slide2 reads md["units"], "metric"/"imperial" -> si/imperial.
    for token, expected, who in (
        ("metric", "si", "geostudio/slide2"),
        ("imperial", "imperial", "geostudio/rs2/slide2"),
        ("unknown", None, "rs2 no-tag"),
    ):
        got = normalize_unit_system(token)
        if got != expected:
            fails.append(f"{who}: normalize_unit_system({token!r}) = {got!r}, "
                         f"expected {expected!r}")
    return fails


# --- v18 template: selector semantics, tseep parse, BC series refs, Ss/Sy rules ----
# These exercise the fileio Phase-3 code end to end (write a real v18 workbook, reload
# it) rather than at the unit level, because the contract IS the file round-trip. They
# need the v18 master template and a seep corpus base file; missing fixtures skip.
_V18_TEMPLATE = os.path.join(_REPO_ROOT, "docs", "inputs", "input_template.xlsx")
_V18_BASE = os.path.join(_REPO_ROOT, "docs", "seep", "files", "xslope_earth_dam1.xlsx")


def _v18_fixtures_present():
    return os.path.exists(_V18_TEMPLATE) and os.path.exists(_V18_BASE)


def _load_silent(path):
    import warnings as _w
    from xslope.fileio import load_slope_data
    with _w.catch_warnings():
        _w.simplefilter("ignore")
        return load_slope_data(path)


def _make_v18(mutate=None, raw_cells=None):
    """Build a v18 workbook from the seep base file: load it, apply `mutate(sd)`, save
    through the v18 master template, then optionally overwrite raw cells (to simulate
    hand edits the writer would not itself produce, e.g. a blank gamma_w). Returns the
    temp path (caller deletes)."""
    import tempfile
    from xslope.fileio import save_slope_data_to_xlsx, write_cells_to_xlsx
    sd = _load_silent(_V18_BASE)
    if mutate:
        mutate(sd)
    tmp = tempfile.mktemp(suffix=".xlsx")
    save_slope_data_to_xlsx(sd, tmp, template=_V18_TEMPLATE)
    if raw_cells:
        write_cells_to_xlsx(tmp, raw_cells)
    return tmp


def _expect_load_error(path, needle, label):
    import warnings as _w
    from xslope.fileio import load_slope_data
    try:
        with _w.catch_warnings():
            _w.simplefilter("ignore")
            load_slope_data(path)
    except Exception as e:  # noqa: BLE001
        if needle in str(e):
            return []
        return [f"{label}: raised but message lacked {needle!r}: {e!r}"]
    return [f"{label}: expected a raise mentioning {needle!r}, none occurred"]


def check_v18_selector_semantics():
    """v18 Units selector drives gamma_w: blank-D10 autofills the canonical value; a
    filled D10 wins and diverging from canonical warns; the writer never leaks the
    template's stock Time default into a model with no time_unit."""
    import os as _os
    import warnings as _w
    if not _v18_fixtures_present():
        return []
    fails = []

    def _si_no_tseep(sd):
        sd["unit_system"] = "si"
        sd["time_unit"] = None
        for m in sd["materials"]:
            m["Ss"] = None
            m["Sy"] = None

    # 1. Selector = SI, D10 blank -> canonical 9.81 autofill.
    tmp = _make_v18(_si_no_tseep, raw_cells={"main": {"D10": None}})
    try:
        sd = _load_silent(tmp)
        if sd["gamma_water"] != GAMMA_W["si"]:
            fails.append(f"blank-D10 autofill: gamma_water {sd['gamma_water']!r} != 9.81")
        if sd["unit_system"] != "si":
            fails.append(f"blank-D10 autofill: unit_system {sd['unit_system']!r} != 'si'")
    finally:
        _os.remove(tmp)

    # 2. Selector = SI, D10 = 11.0 (>2% off) -> D10 wins AND a divergence warning fires.
    tmp = _make_v18(_si_no_tseep, raw_cells={"main": {"D10": 11.0}})
    try:
        with _w.catch_warnings(record=True) as rec:
            _w.simplefilter("always")
            sd = _load_silent_recording(tmp)
        if sd["gamma_water"] != 11.0:
            fails.append(f"filled-D10 override: gamma_water {sd['gamma_water']!r} != 11.0")
        msgs = [str(x.message) for x in rec]
        if not any("gamma_water" in m and "off" in m for m in msgs):
            fails.append(f"filled-D10 override: no divergence warning, got {msgs!r}")
    finally:
        _os.remove(tmp)

    # 3. No time_unit -> the writer must NOT inherit the template's stock "day".
    tmp = _make_v18(_si_no_tseep)
    try:
        sd = _load_silent(tmp)
        if sd["time_unit"] is not None:
            fails.append(f"time_unit leak: got {sd['time_unit']!r}, expected None")
    finally:
        _os.remove(tmp)
    return fails


def _load_silent_recording(path):
    """load_slope_data with warnings NOT suppressed (the caller records them)."""
    from xslope.fileio import load_slope_data
    return load_slope_data(path)


def _transient_mutation(sd, *, time_unit="day", ss=1e-5, sy=0.15, series="res"):
    sd["unit_system"] = "imperial"
    sd["time_unit"] = time_unit
    for m in sd["materials"]:
        m["Ss"] = ss
        m["Sy"] = sy
    sd["tseep"] = {
        "times": [0.0, 10.0, 30.0, 100.0],
        "series": {series: [18.0, 15.0, None, 12.0], "tail": [0.0, None, 1e-6, 2e-6]},
        "duration": 120.0, "save_interval": 5.0, "save_times": [7.5, 42.0],
        "stage_1": 0.0, "stage_2": 100.0,
    }
    if sd["seepage_bc"].get("specified_heads"):
        sd["seepage_bc"]["specified_heads"][0]["head"] = series


def check_v18_tseep_roundtrip():
    """A full tseep block (series table w/ gaps, controls, save_times, stages) plus a
    series-bound head BC round-trips through the writer and loader unchanged."""
    import os as _os
    if not _v18_fixtures_present():
        return []
    fails = []
    tmp = _make_v18(_transient_mutation)
    try:
        sd = _load_silent(tmp)
        ts = sd.get("tseep")
        if ts is None:
            return ["tseep round-trip: reloaded model has no 'tseep' key"]
        want = {
            "times": [0.0, 10.0, 30.0, 100.0],
            "series": {"res": [18.0, 15.0, None, 12.0], "tail": [0.0, None, 1e-6, 2e-6]},
            "duration": 120.0, "save_interval": 5.0, "save_times": [7.5, 42.0],
            "stage_1": 0.0, "stage_2": 100.0,
        }
        for k, v in want.items():
            if ts.get(k) != v:
                fails.append(f"tseep[{k!r}] round-trip: {ts.get(k)!r} != {v!r}")
        head0 = sd["seepage_bc"]["specified_heads"][0]["head"]
        if head0 != "res":
            fails.append(f"series-bound BC: head {head0!r} != 'res'")
        if sd["materials"][0].get("Ss") != 1e-5 or sd["materials"][0].get("Sy") != 0.15:
            fails.append("Ss/Sy did not round-trip on materials")
    finally:
        _os.remove(tmp)
    return fails


def check_v18_tseep_absent_is_legacy():
    """A v18 file with NO transient data adds no 'tseep' key (bit-identical legacy)."""
    import os as _os
    if not _v18_fixtures_present():
        return []
    fails = []

    def _no_tseep(sd):
        sd["unit_system"] = "imperial"
        sd["time_unit"] = None
        for m in sd["materials"]:
            m["Ss"] = None
            m["Sy"] = None

    tmp = _make_v18(_no_tseep)
    try:
        sd = _load_silent(tmp)
        if "tseep" in sd:
            fails.append(f"empty tseep sheet produced a key: {sd['tseep']!r}")
    finally:
        _os.remove(tmp)
    return fails


def check_v18_validation_errors():
    """The four load-time validation rules fire with clear, specific messages."""
    import os as _os
    from xslope.fileio import write_cells_to_xlsx
    if not _v18_fixtures_present():
        return []
    fails = []

    # tseep present but no time_unit declared.
    tmp = _make_v18(lambda sd: _transient_mutation(sd, time_unit=None))
    fails += _expect_load_error(tmp, "Time unit", "tseep-requires-time_unit")
    _os.remove(tmp)

    # Ss required for every material with a tseep sheet.
    tmp = _make_v18(lambda sd: _transient_mutation(sd, ss=None))
    fails += _expect_load_error(tmp, "Ss", "Ss-required")
    _os.remove(tmp)

    # Sy required when unconfined (the base file has an exit-face BC).
    base = _load_silent(_V18_BASE)
    unconfined = bool(base["seepage_bc"].get("exit_face")
                      or base["seepage_bc2"].get("exit_face"))
    if unconfined:
        tmp = _make_v18(lambda sd: _transient_mutation(sd, sy=None))
        fails += _expect_load_error(tmp, "Sy", "Sy-required-unconfined")
        _os.remove(tmp)

    # A BC bound to an undefined series name -> hard error.
    def _bad_series(sd):
        _transient_mutation(sd)
        sd["seepage_bc"]["specified_heads"][0]["head"] = "does_not_exist"
    tmp = _make_v18(_bad_series)
    fails += _expect_load_error(tmp, "not defined", "bad-series-name")
    _os.remove(tmp)

    # A non-numeric BC value with NO tseep sheet -> error (nothing to bind to).
    def _no_tseep(sd):
        sd["unit_system"] = "imperial"
        sd["time_unit"] = None
        for m in sd["materials"]:
            m["Ss"] = None
            m["Sy"] = None
    tmp = _make_v18(_no_tseep)
    write_cells_to_xlsx(tmp, {"seep bc": {"F3": "myseries"}})  # F3 = first block value cell
    fails += _expect_load_error(tmp, "time-series", "string-BC-without-tseep")
    _os.remove(tmp)
    return fails


def check_elastic_props_declared_system():
    """The elastic_props classifier honors a declared unit system when given, and is
    byte-identical to the gamma heuristic when undeclared (feedback: wire the classifier
    to declared units; heuristic stays as the fallback)."""
    import sys
    fails = []
    bench = os.path.join(_REPO_ROOT, "benchmarks", "rocscience")
    if not os.path.exists(os.path.join(bench, "elastic_props.py")):
        return [f"elastic_props fixture missing: {bench}"]
    if bench not in sys.path:
        sys.path.insert(0, bench)
    import elastic_props as ep

    si_mat = {"name": "a", "gamma": 18.0, "c": 60.0, "phi": 25.0}     # heuristic -> SI
    imp_mat = {"name": "b", "gamma": 120.0, "c": 1250.0, "phi": 25.0}  # heuristic -> Imperial

    # Undeclared: identical to passing the heuristic flag explicitly.
    if ep.classify(si_mat, ep.is_imperial([si_mat])) != ep.classify(si_mat, False):
        fails.append("elastic_props: undeclared SI classify drifted from heuristic")

    # A declared system OVERRIDES the passed flag.
    if ep.classify(imp_mat, True, declared_system="si") != ep.classify(imp_mat, False):
        fails.append("elastic_props: declared_system='si' did not override imperial=True")
    if ep.classify(si_mat, False, declared_system="imperial") != ep.classify(si_mat, True):
        fails.append("elastic_props: declared_system='imperial' did not override imperial=False")

    # 'metric' folds to SI; an unrecognized token falls back to the passed flag.
    if ep.classify(imp_mat, True, declared_system="metric") != ep.classify(imp_mat, False):
        fails.append("elastic_props: declared_system='metric' should behave as SI")
    if ep.classify(si_mat, True, declared_system="bogus") != ep.classify(si_mat, True):
        fails.append("elastic_props: an unrecognized declared_system should fall back to the flag")

    # assign_elastic_props: undeclared == declared-consistent (idempotent, unchanged builds).
    import copy
    m1, m2 = copy.deepcopy([si_mat]), copy.deepcopy([si_mat])
    r1 = [(x[1], x[2], x[3]) for x in ep.assign_elastic_props(m1)]
    r2 = [(x[1], x[2], x[3]) for x in ep.assign_elastic_props(m2, declared_system="si")]
    if r1 != r2:
        fails.append(f"elastic_props: assign undeclared {r1!r} != declared-consistent {r2!r}")
    return fails


def main():
    print("xslope.units unit checks:")
    checks = [
        ("constants", check_constants),
        ("require_gamma_water", check_require_gamma_water),
        ("infer_unit_system", check_infer_unit_system),
        ("labels", check_labels),
        ("seep entry point", check_seep_entry_point_raises),
        ("fem entry point", check_fem_entry_point_raises),
        ("infer_system_from_gamma_water", check_infer_system_from_gamma_water),
        ("normalize_unit_system", check_normalize_unit_system),
        ("units_check", check_units_check),
        ("loader time never inferred", check_loader_time_never_inferred),
        ("s2d import unit_system", check_s2d_import_unit_system),
        ("importer persistence mapping", check_importer_persistence_mapping),
        ("v18 selector semantics", check_v18_selector_semantics),
        ("v18 tseep round-trip", check_v18_tseep_roundtrip),
        ("v18 tseep absent = legacy", check_v18_tseep_absent_is_legacy),
        ("v18 validation errors", check_v18_validation_errors),
        ("elastic_props declared system", check_elastic_props_declared_system),
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
