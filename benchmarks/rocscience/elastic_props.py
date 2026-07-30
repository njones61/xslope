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

E and nu are picked from the table in docs/fem/overview.md by soil type, which we infer
from (c, phi) and the strength option. That soil-type table and the classifier that
reads it live in :func:`xslope.units.classify_elastic` — the vendor importers need the
same last-resort fallback for a file that states no elastic pair, so the table is in the
shipped package and this module is the builder-side face of it.
"""

import math

from xslope.units import KPA_TO_PSF, classify_elastic  # noqa: F401  (KPA_TO_PSF re-export)


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


#: Builder-side alias for the shipped soil-type classifier. It carries the table (see
#: :func:`xslope.units.classify_elastic`); this module only names it, so a builder and
#: a vendor import that both fall back to a soil type land on the same constants.
classify = classify_elastic


# The unit-blind value every material inherited from the base template. ONLY these get
# reassigned — anything else was set deliberately by a builder and must be left alone.
# (hammah_hb1 carries a deliberate E=5e6 for HB rock.)
#
# This pair is a SENTINEL, not a measurement, and 100000/0.3 is also a perfectly ordinary
# thing for a vendor model to specify — RS2 #067 specifies exactly it. Reading the sentinel
# off the value alone therefore cannot tell "nobody set this" from "the vendor set this to
# the same number", and the classifier used to overwrite the vendor's own constants on all
# six rs2_67 files for that reason. The fix is not a different number: it is that whoever
# DID set a value says so, via ``pinned`` below.
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


def assign_elastic_props(materials, force=False, declared_system=None, pinned=None):
    """Set E and nu in place on every material that does NOT already carry a
    deliberately-set, non-default value — i.e. one whose E is unset (None / 0 / NaN)
    or still equal to the inherited unit-blind default. Materials that carry a real,
    non-default E/nu (vendor .fez models, HB rock) are left untouched unless force=True.

    ``pinned`` is the set of material NAMES a caller has already spoken for — pass
    ``vendor_tcut.apply_vendor_e_nu``'s return value. Those are skipped whatever their
    values are, which is what keeps a vendor pair that happens to equal
    :data:`INHERITED_DEFAULT` from being read as "nobody set this" and reclassified.

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
    pinned = set(pinned or ())
    out = []
    for m in materials:
        name = str(m.get('name', '?'))
        E_old = finite(m.get('E'))
        nu_old = finite(m.get('nu'))
        if name.strip() in pinned and not force:
            out.append((name, 'KEPT (vendor model)', E_old, nu_old, False))
            continue
        unset = E_old <= 0.0
        inherited = (round(E_old, 3), round(nu_old, 3)) == INHERITED_DEFAULT
        deliberate = not (unset or inherited)
        soil, E, nu = classify(m, imperial, declared_system=declared_system)
        if deliberate and not force:
            out.append((name, 'KEPT (set deliberately)', E_old, nu_old, False))
            continue
        m['E'] = float(E)
        m['nu'] = float(nu)
        out.append((name, soil, E, nu, True))
    return out


#: Back-compat alias — run_tests.run_fem_elastic_units_test imports this name.
_finite = finite
