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
- :func:`infer_unit_system` -- the gamma-magnitude heuristic, for cross-check
  warnings and undeclared vendor imports ONLY, never to silently set behavior.

The declaration itself (a ``unit_system`` / ``time_unit`` on ``slope_data``, the
template selector, importer persistence, and the labels applied to Studio forms and
plots) lands in the later phases of the units plan; this module is the fix-first
foundation they build on.
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
