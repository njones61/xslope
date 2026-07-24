# Copyright 2025 Norman L. Jones
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Unit-system single source of truth for xslope.

xslope is deliberately unit-agnostic: the solver core never converts, it consumes
the numbers as entered and returns results in the same system. What this module
adds is the small amount of unit *bookkeeping* that must be consistent everywhere:

- :data:`GAMMA_W` -- the canonical unit weight of water in each declared system.
  This is the ONLY place the 9.81 / 62.4 constants live; every other module imports
  them from here. The values are used for *seeding* a new project and for *detection*
  (inferring a system from a gamma_water magnitude), never as a silent runtime
  fallback.
- :func:`require_gamma_water` -- read a model's ``gamma_water`` with NO silent
  default. A missing value raises loudly rather than being replaced by a hardcoded
  constant, because a substituted constant of the wrong system rescales every pore
  pressure in the model -- a far worse failure than a hard stop (the "u-option
  lesson": a silent default that can flip the unit system is worse than an error).
- :func:`labels` -- the display unit strings for a declared system (length, stress,
  unit weight, ...). Consumers arrive in a later phase; it is defined and tested now
  so the contract is fixed.
- :func:`infer_unit_system` -- the gamma-magnitude heuristic (from material unit
  weights), for cross-check warnings and undeclared vendor imports ONLY, never to
  silently set behavior.
- :func:`infer_system_from_gamma_water` -- the sibling heuristic that reads the ONE
  quantity physics fixes (the unit weight of *water*) against the SI/Imperial
  magnitude bands. Used to label pre-selector (<= v17) template files and s2d imports,
  where gamma_water is the only declaration available. Returns ``None`` (unlabeled)
  for anything off-band rather than guessing.
- :func:`normalize_unit_system` -- fold every vocabulary the importers already speak
  (``"metric"`` / ``"Metric"`` / ``"si"`` / ``"imperial"`` / ``"unknown"``) onto the
  canonical ``"si"`` / ``"imperial"`` / ``None`` that ``slope_data["unit_system"]``
  carries, so persistence is one line at every import site.
- :func:`units_check` -- the load/import-time sanity checks (warnings, never errors):
  gamma_water vs the declared system's canonical value, the material-gamma heuristic
  vs the declared system, and a soft Young's-modulus magnitude band. Skipped entirely
  when the model declares no system (returns ``[]``).

The declaration itself is a ``unit_system`` ("si" / "imperial" / ``None``) and
``time_unit`` (str / ``None``) pair on ``slope_data``. Phase 2 introduces the code
paths that SET them (the loader's blank-inference for <= v17 files, and importer
persistence) and CONSUME them (``units_check``, and FEM/seep export provenance),
without the v18 template selector (Phase 3) or the Studio/plot labels (Phase 4).
``time_unit`` is NEVER inferred for template files -- a wrong time label is worse than
none (the published min-vs-sec mislabel lesson) -- so it stays ``None`` until the v18
cell or a vendor import genuinely declares it for a quantity we ingest.
"""

import math

#: Canonical unit weight of water per declared unit system -- the ONLY place these
#: two constants live. ``"si"`` is 9.81 kN/m3 (Norm, 2026-07-23); ``"imperial"`` is
#: 62.4 lb/ft3. Read a *model's* gamma_water with :func:`require_gamma_water`; these
#: canonical values are for seeding and magnitude detection, not runtime fallback.
GAMMA_W = {"si": 9.81, "imperial": 62.4}


def require_gamma_water(data, where=None):
    """Return ``data["gamma_water"]``, or raise a clear error when it is absent.

    Reads the unit weight of water with **no silent default**. Every solver, export,
    and plotting-math site that needs gamma_w goes through here, so that a genuinely
    missing value fails loudly instead of being quietly replaced by a hardcoded 9.81
    or 62.4. Substituting a canonical value of the wrong system would rescale every
    pore pressure in the model -- the exact silent unit-flip the units plan exists to
    prevent.

    Parameters
    ----------
    data : mapping
        Any mapping that carries a ``gamma_water`` key -- ``slope_data``,
        ``fem_data``, a parsed vendor water block, and so on.
    where : str, optional
        A short context phrase woven into the error message (e.g.
        ``"seepage analysis"``), naming the site that needed the value.

    Raises
    ------
    ValueError
        If ``gamma_water`` is missing or ``None``.
    """
    g = data.get("gamma_water")
    if g is None:
        ctx = f" ({where})" if where else ""
        raise ValueError(
            f"gamma_water (unit weight of water) is required but missing{ctx}. "
            f"There is deliberately no default: substituting a canonical value "
            f"({GAMMA_W['si']} SI / {GAMMA_W['imperial']} Imperial) could silently "
            f"flip the model's unit system and rescale every pore pressure. Set "
            f"gamma_water on the model."
        )
    return float(g)


def infer_unit_system(materials):
    """Guess the unit system from material unit weights, or ``None`` if there is no
    signal.

    The two systems' unit weights sit in far-apart bands -- soil gamma is roughly
    15-25 kN/m3 in SI but 90-160 pcf in Imperial -- so a single threshold on the
    largest material gamma separates them cleanly: ``max(gamma) > 50`` is Imperial,
    otherwise SI. Promoted from ``benchmarks/rocscience/elastic_props.is_imperial``.

    For CROSS-CHECK WARNINGS and UNDECLARED vendor imports ONLY. Never call this to
    silently *set* a model's behavior: a declared unit system always wins, and an
    inferred one is only ever the basis of a warning (e.g. "declared SI, but the
    gammas look Imperial").

    Parameters
    ----------
    materials : iterable of mapping
        Material dictionaries, each optionally carrying a numeric ``"gamma"``.

    Returns
    -------
    {"si", "imperial", None}
        ``None`` when no positive, finite gamma is present (no signal to infer from).
    """
    gammas = []
    for m in materials or []:
        try:
            g = float(m.get("gamma"))
        except (TypeError, ValueError):
            continue
        if math.isnan(g) or math.isinf(g) or g <= 0:
            continue
        gammas.append(g)
    if not gammas:
        return None
    return "imperial" if max(gammas) > 50.0 else "si"


#: Half-width of the magnitude band around each canonical unit weight of water that
#: :func:`infer_system_from_gamma_water` accepts as a match. The SI band
#: (9.81 +/- 0.15 -> 9.66..9.96) also captures GeoStudio's vendor canonical 9.807, so
#: a metric .gsz whose gamma_w round-trips as 9.807 still labels SI. The Imperial band
#: (62.4 +/- 1.0) is wider in absolute terms, matching the bands the GeoStudio importer
#: already uses. Both are tight enough that a seawater/brine override (~10.05 / ~64)
#: falls OFF-band and labels ``None`` -- the entered value is still used for physics;
#: it is only left unlabeled.
_GAMMA_W_BAND = {"si": 0.15, "imperial": 1.0}


def infer_system_from_gamma_water(gamma_water):
    """Guess the unit system from the unit weight of *water*, or ``None`` if off-band.

    The unit weight of water is the one quantity physics pins, so its magnitude is a
    reliable label for a file that declares no system: ~9.81 (or GeoStudio's 9.807)
    is SI, ~62.4 is Imperial. This is the blank-inference the <= v17 loader applies and
    the detection the s2d importer uses -- neither of which carries an explicit
    selector. Anything outside both bands (including a seawater override) returns
    ``None``: the value is still used as entered, it is simply left unlabeled rather
    than mislabeled.

    Parameters
    ----------
    gamma_water : float or None
        A model's unit weight of water. ``None`` / non-numeric returns ``None``.

    Returns
    -------
    {"si", "imperial", None}
    """
    if gamma_water is None:
        return None
    try:
        g = float(gamma_water)
    except (TypeError, ValueError):
        return None
    if math.isnan(g) or math.isinf(g):
        return None
    for system, canonical in GAMMA_W.items():
        if abs(g - canonical) < _GAMMA_W_BAND[system]:
            return system
    return None


def normalize_unit_system(value):
    """Fold an importer/vendor unit-system token onto ``"si"`` / ``"imperial"`` / ``None``.

    ``slope_data["unit_system"]`` speaks exactly three values, but the importers detect
    in their own vendor vocabulary -- GeoStudio and Slide2 say ``"metric"``, RS2 says
    ``"unknown"`` when its RFCunits string carries no tag. This maps them all so each
    import site persists with a single call. Case-insensitive; ``"metric"`` -> ``"si"``;
    any unrecognized token (``"unknown"``, ``""``, ``None``) -> ``None`` so an
    undetermined system stays undeclared rather than guessed.

    Parameters
    ----------
    value : str or None
        A unit-system token in any casing: ``"si"``, ``"metric"``, ``"imperial"``,
        ``"unknown"``, ``None``, ...

    Returns
    -------
    {"si", "imperial", None}
    """
    if value is None:
        return None
    v = str(value).strip().lower()
    if v in ("si", "metric"):
        return "si"
    if v == "imperial":
        return "imperial"
    return None


#: Plausible Young's-modulus magnitude band per declared system, for the soft
#: :func:`units_check` warning. Derived from the docs' per-system E table
#: (``docs/fem/overview.md``: SI ~2,000-10,000,000 kPa across soft clay to soft rock;
#: Imperial ~41,800-208,800,000 psf), widened with margin so only a value that is
#: clearly the WRONG system's magnitude trips it (e.g. an E entered in kPa on an
#: Imperial model, or vice versa). A soft cross-check, never a hard rule.
_E_BOUNDS = {"si": (1_000.0, 15_000_000.0),
             "imperial": (20_000.0, 300_000_000.0)}

#: Fractional divergence of a model's gamma_water from its declared canonical value
#: that :func:`units_check` will warn about (~2%). The seawater/brine override
#: (~10.05 kN/m3 / ~64 pcf, ~2.4% high) diverges past this and DOES draw a warning --
#: that is legal, the warning just names both values so the override is a conscious
#: choice, not a silent unit flip.
_GAMMA_W_DIVERGENCE = 0.02


def units_check(slope_data):
    """Load/import-time unit sanity checks -- a list of warning strings, never a raise.

    Returns the messages so the caller surfaces them its own way: the loader emits
    them through the :mod:`warnings` module, the importers append them to their caveat
    list. Every check is a WARNING, not an error -- a mismatched unit label never
    blocks a model from loading. When the model declares no system
    (``unit_system`` is ``None``) the function returns ``[]`` immediately: with nothing
    declared there is nothing to cross-check, and an undeclared (legacy) model must be
    exactly as silent as it is today.

    The three checks:

    1. **gamma_water vs the declared system.** More than ~2% off the canonical value
       (:data:`GAMMA_W`) draws a warning naming BOTH values. The seawater/brine
       override is legal and still trips this -- the warning documents the choice.
    2. **material gammas vs the declared system.** :func:`infer_unit_system` on the
       materials disagreeing with the declaration (e.g. "declared SI, gammas look
       Imperial") draws a warning; the bands are far apart, so false positives are
       unlikely.
    3. **Young's-modulus magnitude.** Any material whose ``E`` is set but lands
       outside the declared system's plausible band (:data:`_E_BOUNDS`) draws a soft
       warning -- the classic "E entered in the wrong system's units" mistake.

    Parameters
    ----------
    slope_data : mapping
        A model dict carrying at least ``unit_system``; optionally ``gamma_water`` and
        ``materials``.

    Returns
    -------
    list of str
        Warning messages, in check order. Empty when ``unit_system`` is ``None`` or no
        check fired.
    """
    system = normalize_unit_system(slope_data.get("unit_system"))
    if system is None:
        return []

    warnings_out = []
    canonical = GAMMA_W[system]

    # 1. gamma_water vs the declared system's canonical value.
    g = slope_data.get("gamma_water")
    if g is not None:
        try:
            gw = float(g)
        except (TypeError, ValueError):
            gw = None
        if gw is not None and not (math.isnan(gw) or math.isinf(gw)) and canonical > 0:
            if abs(gw - canonical) / canonical > _GAMMA_W_DIVERGENCE:
                warnings_out.append(
                    f"declared unit system is {system!r}, whose canonical unit weight "
                    f"of water is {canonical:g}, but this model uses gamma_water = "
                    f"{gw:g} ({abs(gw - canonical) / canonical * 100:.1f}% off). If "
                    f"that is a deliberate seawater/brine value this is fine; otherwise "
                    f"check the unit system.")

    # 2. Material unit weights vs the declared system.
    materials = slope_data.get("materials") or []
    inferred = infer_unit_system(materials)
    if inferred is not None and inferred != system:
        warnings_out.append(
            f"declared unit system is {system!r}, but the material unit weights look "
            f"{inferred!r} (max gamma {'> 50' if inferred == 'imperial' else '<= 50'}). "
            f"Check the unit system -- xslope does not convert, so a mislabel means the "
            f"numbers are being read in the wrong system.")

    # 3. Young's-modulus magnitude (soft cross-check).
    lo, hi = _E_BOUNDS[system]
    off_band = []
    for m in materials:
        try:
            E = float(m.get("E"))
        except (TypeError, ValueError):
            continue
        if math.isnan(E) or math.isinf(E) or E <= 0:
            continue
        if E < lo or E > hi:
            off_band.append((str(m.get("name", "?")), E))
    if off_band:
        shown = ", ".join(f"{name} (E={E:g})" for name, E in off_band[:4])
        more = "" if len(off_band) <= 4 else f", and {len(off_band) - 4} more"
        unit = "kPa" if system == "si" else "psf"
        warnings_out.append(
            f"declared unit system is {system!r} (E in {unit}), but "
            f"{'these materials have' if len(off_band) > 1 else 'this material has'} a "
            f"Young's modulus outside the plausible {system} range "
            f"[{lo:g}, {hi:g}] {unit}: {shown}{more}. Check that E was entered in "
            f"{unit}, not the other system's units.")

    return warnings_out


#: Human-readable system token for the ``# units:`` export provenance comment --
#: the same "SI" / "Imperial" wording the v18 template selector shows, so a results
#: CSV and the file it came from name the system identically.
_SYSTEM_DISPLAY = {"si": "SI", "imperial": "Imperial"}


def units_comment_line(unit_system, time_unit=None):
    """A ``# units: ...`` provenance comment for a results CSV, or ``""`` undeclared.

    Written at the TOP of the seep / FEM result CSVs so a saved solution records the
    unit system it was computed in. Returns the empty string when the model declares no
    system (:func:`normalize_unit_system` of ``unit_system`` is ``None``) -- an
    undeclared model's exports stay byte-identical to today. When a system is declared
    the line is ``"# units: SI\\n"``, gaining ``", time: <t>"`` only when a time unit is
    also declared::

        # units: SI, time: day

    The leading ``#`` makes it an ordinary CSV comment: readers skip it with
    ``pandas.read_csv(..., comment="#")``, exactly as they already skip the trailing
    ``# Total Flowrate:`` footer.

    Parameters
    ----------
    unit_system : {"si", "imperial", None} or a vendor token
        Passed through :func:`normalize_unit_system`, so ``"metric"`` etc. are accepted.
    time_unit : str, optional
        A time token (``"day"``, ``"sec"``, ...); omitted from the line when falsy.

    Returns
    -------
    str
        ``""`` when undeclared, else a single newline-terminated comment line.
    """
    system = normalize_unit_system(unit_system)
    if system is None:
        return ""
    line = f"# units: {_SYSTEM_DISPLAY[system]}"
    t = str(time_unit).strip() if time_unit else ""
    if t:
        line += f", time: {t}"
    return line + "\n"


#: The keys every :func:`labels` result carries, in a stable order.
_LABEL_KEYS = ("length", "stress", "unit_weight", "force_per_len",
               "k", "flowrate", "time")


def labels(unit_system, time_unit=None):
    """Display unit strings for a declared unit system.

    Parameters
    ----------
    unit_system : {"si", "imperial", None}
        The declared system (case-insensitive), or ``None`` for an undeclared
        (legacy) model.
    time_unit : str, optional
        A time token (``"s"``, ``"sec"``, ``"min"``, ``"hr"``, ``"day"``, ...) or
        ``None``. Never inferred -- an unlabeled time base stays unlabeled.

    Returns
    -------
    dict
        A dict with these string keys, ready to append as ``f" ({lbl['stress']})"``:

        ============= ================= =========================
        key           SI                Imperial
        ============= ================= =========================
        length        ``"m"``           ``"ft"``
        stress        ``"kPa"``         ``"psf"``
        unit_weight   ``"kN/m³"``       ``"pcf"``
        force_per_len ``"kN/m"``        ``"lb/ft"``
        k             ``"m/<t>"``       ``"ft/<t>"``
        flowrate      ``"m³/<t> per m"``  ``"ft³/<t> per ft"``
        time          ``"<t>"``         ``"<t>"``
        ============= ================= =========================

        When ``unit_system`` is ``None`` every value is an empty string, so a legacy
        (undeclared) model degrades to today's bare, unit-less labels. The
        time-bearing keys (``k``, ``flowrate``, ``time``) are also empty when
        ``time_unit`` is ``None``: the time unit is never guessed, so an unlabeled
        time base stays unlabeled rather than inventing "m/?" (the min-vs-sec
        mislabel lesson).

    Raises
    ------
    ValueError
        If ``unit_system`` is neither ``None`` nor a recognized system.
    """
    out = {k: "" for k in _LABEL_KEYS}
    if unit_system is None:
        return out
    key = str(unit_system).strip().lower()
    if key not in GAMMA_W:
        raise ValueError(
            f"Unknown unit system {unit_system!r}; expected 'si', 'imperial', or None."
        )
    if key == "si":
        length, out["stress"], out["unit_weight"], out["force_per_len"] = (
            "m", "kPa", "kN/m³", "kN/m")
    else:
        length, out["stress"], out["unit_weight"], out["force_per_len"] = (
            "ft", "psf", "pcf", "lb/ft")
    out["length"] = length
    if time_unit:
        t = str(time_unit).strip()
        out["time"] = t
        out["k"] = f"{length}/{t}"
        out["flowrate"] = f"{length}³/{t} per {length}"
    return out
