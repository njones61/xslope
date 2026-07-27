"""Unit-system and elastic-property (E, nu) bookkeeping for the benchmark builders.

The two belong together: which E a material gets depends on the unit system the file
is written in, and both are things a builder inherits from whatever donor file it
loaded rather than declaring for itself. :func:`resolve_unit_system` fixes the label,
:func:`assign_elastic_props` fixes the constants.

The benchmark builders copy a material dict from a base template and never set E/nu, so
every corpus material inherited a single unit-blind E = 100,000 with nu = 0.3. For the
IMPERIAL problems that means E = 100,000 *psf* (~4.8 MPa) — below the docs' own softest
recommendation (soft clay, 41,800 psf) for what are rockfill dam shells. It does not change
the factor of safety (a perfectly-plastic collapse load is independent of the elastic
constants, and this was verified: FS is invariant to E across a 100x sweep), but it produces
absurd displacements (max|u| = 324 ft on vp077b) and misleads anyone reading the deformation
output.

E and nu are picked from the table in docs/fem/overview.md by soil type, which we infer from
(c, phi) and the strength option. Values are the MIDPOINT of each published range.

    Soil Type      E [kPa]              nu
    Soft Clay      2,000 -  15,000      0.40-0.50
    Medium Clay    15,000 -  50,000     0.35-0.45
    Stiff Clay     50,000 - 200,000     0.20-0.40
    Loose Sand     10,000 -  25,000     0.25-0.35
    Medium Sand    25,000 -  75,000     0.30-0.40
    Dense Sand     75,000 - 200,000     0.25-0.35
    Gravel         100,000 - 500,000    0.15-0.30
    Rock Fill      50,000 - 300,000     0.20-0.35
    Soft Rock      1,000,000 - 10,000,000  0.15-0.30
"""

import math

KPA_TO_PSF = 20.8854342

# (name, E_kPa, nu) — E is the midpoint of the published range
_SOFT_CLAY   = ('Soft Clay',    8_000.0,   0.45)
_MEDIUM_CLAY = ('Medium Clay',  32_000.0,  0.40)
_STIFF_CLAY  = ('Stiff Clay',   125_000.0, 0.30)
_LOOSE_SAND  = ('Loose Sand',   17_000.0,  0.30)
_MEDIUM_SAND = ('Medium Sand',  50_000.0,  0.35)
_DENSE_SAND  = ('Dense Sand',   137_000.0, 0.30)
_ROCK_FILL   = ('Rock Fill',    175_000.0, 0.28)
_SOFT_ROCK   = ('Soft Rock',    5_500_000.0, 0.22)


def finite(x):
    """A finite float, or 0.0 for None / NaN / inf / non-numeric — so an unset or
    NaN cell (fileio reads a blank as NaN) reads as 'no value', not a real number."""
    try:
        v = float(x)
    except (TypeError, ValueError):
        return 0.0
    return 0.0 if (math.isnan(v) or math.isinf(v)) else v


def is_imperial(materials):
    """Unit system from unit weight: metric gamma ~ 15-25 kN/m3, imperial ~ 90-160 pcf.
    The two bands are far apart, so a single threshold is safe."""
    gammas = [finite(m.get('gamma')) for m in materials]
    gammas = [g for g in gammas if g > 0]
    if not gammas:
        return False
    return max(gammas) > 50.0


def _declared_imperial(declared_system):
    """Resolve an explicit unit-system declaration to the ``imperial`` bool, or ``None``
    when nothing is declared.

    ``'si'`` / ``'metric'`` -> ``False``, ``'imperial'`` -> ``True``, anything else
    (including ``None``) -> ``None``. This lets a caller that KNOWS the model's unit
    system (``slope_data['unit_system']``, from the v18 selector or a vendor import)
    override the gamma-magnitude heuristic below — closing the "wire the classifier to
    declared units" follow-up (feedback_fem_E_nu_units). When it returns ``None`` the
    heuristic (:func:`is_imperial`) still decides, so behavior is unchanged for the
    undeclared benchmark builds that call this today.
    """
    if declared_system is None:
        return None
    d = str(declared_system).strip().lower()
    if d in ('si', 'metric'):
        return False
    if d == 'imperial':
        return True
    return None


def classify(mat, imperial=False, declared_system=None):
    """Return (soil_type, E, nu) in the model's own units.

    ``declared_system`` (``'si'`` / ``'imperial'`` / ``None``), when given, OVERRIDES the
    ``imperial`` flag with the model's declared unit system — use it instead of the
    gamma-magnitude heuristic when the system is known. When ``None`` (the default) the
    passed ``imperial`` flag is honored unchanged.
    """
    decl = _declared_imperial(declared_system)
    if decl is not None:
        imperial = decl
    opt = str(mat.get('option', 'mc')).strip().lower()
    c = float(mat.get('c', 0) or 0)
    phi = float(mat.get('phi', 0) or 0)

    # Rock strength models (Hoek-Brown) -> rock, not soil.
    if opt in ('hb', 'hoek', 'hoek-brown'):
        soil = _SOFT_ROCK
    else:
        # c is in the model's own units; compare against a kPa threshold converted if needed.
        c_kpa = c / KPA_TO_PSF if imperial else c
        if phi >= 1.0 and c_kpa < 1.0:
            # purely frictional
            soil = _ROCK_FILL if phi >= 40.0 else (_DENSE_SAND if phi >= 33.0 else
                                                   (_MEDIUM_SAND if phi >= 28.0 else _LOOSE_SAND))
        elif phi < 1.0:
            # purely cohesive (undrained clay)
            soil = (_SOFT_CLAY if c_kpa < 25.0 else
                    (_MEDIUM_CLAY if c_kpa < 75.0 else _STIFF_CLAY))
        else:
            # c-phi soil: stiffness driven by the cohesive component
            soil = (_MEDIUM_CLAY if c_kpa < 50.0 else _STIFF_CLAY)

    name, E_kpa, nu = soil
    E = E_kpa * KPA_TO_PSF if imperial else E_kpa
    return name, round(E, -2), nu


# The unit-blind value every material inherited from the base template. ONLY these get
# reassigned — anything else was set deliberately by a builder and must be left alone.
# (rs2_56/57/58 carry E=5000 with nu=0.35 vs 0.30 as a POISSON-RATIO STUDY PAIR; overwriting
# their nu would destroy the benchmark. hammah_hb1 carries a deliberate E=5e6 for HB rock.)
INHERITED_DEFAULT = (100_000.0, 0.3)


def resolve_unit_system(slope_data):
    """Label a builder's model from ITS OWN numbers and store it on ``slope_data``.

    Almost every corpus builder starts from ``load_slope_data(<donor file>)`` and then
    rewrites the geometry and properties -- including switching the model to English
    units by setting ``gamma_water = 62.4`` and pcf unit weights. What it does NOT
    rewrite is ``unit_system``, which arrives from the metric donor as ``'si'`` and,
    since the v18 writer persists that label into the Units selector, silently
    relabels an English-unit file as SI. This re-derives the label from the model as
    it now stands, so a builder never has to remember to update it:

      1. the unit weight of WATER, the one quantity physics pins (9.81 / 62.4);
      2. failing that (an off-band or absent gamma_w), the material unit weights.

    Both heuristics come from :mod:`xslope.units`, so the label matches what the
    loader would infer for the same file. Returns the resolved system (or ``None``
    when neither heuristic has a signal, in which case ``unit_system`` is left as-is).
    """
    from xslope.units import infer_system_from_gamma_water, infer_unit_system
    resolved = (infer_system_from_gamma_water(slope_data.get('gamma_water'))
                or infer_unit_system(slope_data.get('materials', [])))
    if resolved:
        slope_data['unit_system'] = resolved
    return resolved


def assign_elastic_props(materials, force=False, declared_system=None):
    """Set E and nu in place on every material that does NOT already carry a
    deliberately-set, non-default value — i.e. one whose E is unset (None / 0 / NaN)
    or still equal to the inherited unit-blind default. Materials that carry a real,
    non-default E/nu (vendor .fez models, the rs2_56/57/58 Poisson study pair, HB
    rock) are left untouched unless force=True.

    ``declared_system`` (``'si'`` / ``'imperial'`` / ``None``): when a caller knows the
    model's unit system (e.g. ``slope_data['unit_system']``), pass it to use it instead
    of the gamma-magnitude heuristic. When ``None`` (the default) the heuristic
    :func:`is_imperial` decides exactly as before, so existing builds are unchanged.

    Idempotent: once a material has been given a physical soil-type E it reads as
    deliberate on the next call and is preserved, so a plain rebuild reproduces the
    same file.

    Returns a list of (name, soil_type_or_status, E, nu, changed)."""
    decl = _declared_imperial(declared_system)
    imperial = is_imperial(materials) if decl is None else decl
    out = []
    for m in materials:
        E_old = finite(m.get('E'))
        nu_old = finite(m.get('nu'))
        unset = E_old <= 0.0
        inherited = (round(E_old, 3), round(nu_old, 3)) == INHERITED_DEFAULT
        deliberate = not (unset or inherited)
        soil, E, nu = classify(m, imperial, declared_system=declared_system)
        if deliberate and not force:
            out.append((str(m.get('name', '?')), 'KEPT (set deliberately)', E_old, nu_old, False))
            continue
        m['E'] = float(E)
        m['nu'] = float(nu)
        out.append((str(m.get('name', '?')), soil, E, nu, True))
    return out


#: Back-compat alias — run_tests.run_fem_elastic_units_test imports this name.
_finite = finite
