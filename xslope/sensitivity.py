"""Sensitivity analysis: vary one input over a range and report how the answer moves.

The geotechnical staple (Duncan & Wright's FS-vs-parameter charts and tornado
diagrams), distinct from reliability(): reliability perturbs parameters by
+/-sigma to estimate a distribution; sensitivity sweeps a user range and
asserts nothing statistical — which also makes it the designated tool for
correlated fit coefficients (power-curve A and b) that MFOSM must not touch.

Two ways to say what varies, mutually exclusive per call:

- ``param``: a validated parameter reference "kind:name:field" (or tuple),
  e.g. "mat:Clay:c", "reinforce:Row 2:t_max", "piles:Pile 1:H",
  "global:k_seismic", "geom:piezo:dy". Every ref resolves to a setter.
- ``modify``: a user callable ``(slope_data, value) -> slope_data`` plus a
  ``label`` — the escape hatch that makes any change sweepable (geometry
  especially; see main_design.py's set_slope_angle for the archetype). The
  engine treats built-in refs and user setters identically.

Sweep inputs are API-only by design — sensitivity describes an analysis you
run, not a property of the model, so nothing here reads or writes template
columns.
"""

import copy
import time

import numpy as np
import pandas as pd

__all__ = ['sensitivity', 'tornado', 'design', 'back_analysis', 'set_param',
           'resolve_param', 'list_params', 'tornado_from_sweeps',
           'scaled_sensitivity', 'variance_contribution', 'mc_rank_correlation',
           'fs_vs_time']


# ---------------------------------------------------------------------------
# Parameter addressing
# ---------------------------------------------------------------------------

# fields addressable on every material regardless of strength option
_MAT_GENERAL_FIELDS = ('gamma', 'gamma_sat', 'ru', 'd', 'psi')
# seep fields live on the material rows too; settable here so a later coupled
# seep+lem analysis reuses the same refs
_SEEP_FIELDS = ('k1', 'k2', 'alpha', 'kr0', 'h0')
_REINFORCE_FIELDS = ('t_max', 't_res', 'lp1', 'lp2', 'tend1', 'tend2', 'spacing')
# plan-name -> slope_data key for pile fields
_PILE_FIELDS = {'H': 'H', 'theta': 'theta_p', 'D': 'D_pile', 'S': 'S',
                'V_cap': 'V_cap', 'M_cap': 'M_cap'}
_GLOBAL_FIELDS = ('k_seismic', 'tcrack_depth', 'tcrack_water')


def _mat_strength_fields(material, mat_name):
    """Strength fields addressable on this material, per its option — the same
    option-awareness reliability() enforces (sweeping c on a cp material is the
    same class of silent error as u='peizo'). Power-curve coefficients and
    Hoek-Brown constants are addressable HERE and deliberately not in
    reliability(): the pow_* are correlated fit coefficients and mb/s/a are all
    derived from GSI/mi/D, so neither set takes an independent standard
    deviation. Sensitivity is the designated tool for those."""
    if material.get('option') == 'pow':
        return ('pow_a', 'pow_b', 'pow_c', 'pow_d')
    if material.get('option') == 'hb':
        return ('hb_sci', 'hb_gsi', 'hb_mi', 'hb_d')
    from .advanced import _strength_param_mapping
    mapping = _strength_param_mapping(material, mat_name)   # raises on unknown option
    return tuple(k for k in mapping.keys() if k != 'gamma')


def _find_by_name(items, name, kind, name_key='name'):
    """Case-insensitive unique lookup by name, or by 1-based index when ``name`` is
    an integer. Raises naming what was given and what exists."""
    # An integer (not the bool subclass) addresses the item by 1-based position —
    # the "material index or name" the AI/GUI specs allow.
    if isinstance(name, int) and not isinstance(name, bool):
        if 1 <= name <= len(items):
            return name - 1, items[name - 1]
        raise ValueError(f"{kind} index {name} out of range (1..{len(items)}).")
    want = str(name).strip().lower()
    hits = [(i, it) for i, it in enumerate(items)
            if str(it.get(name_key, '')).strip().lower() == want]
    if len(hits) == 1:
        return hits[0]
    have = [str(it.get(name_key, '')) or f'<unnamed #{i+1}>' for i, it in enumerate(items)]
    if not hits:
        raise ValueError(f"No {kind} named '{name}'. Available: {have}")
    raise ValueError(f"{len(hits)} {kind}s share the name '{name}' — rename them "
                     f"so the reference is unambiguous. Available: {have}")


def _set_material_field(sd, idx, field, value, couple_gamma=True):
    """THE single mutation path for material fields, shared with reliability()'s
    _perturbed_slope_data. gamma and gamma_sat are the same soil weighed two
    ways (correlation ~1), so setting gamma moves gamma_sat by the same
    absolute delta; gamma_sat stays separately addressable for direct sweeps."""
    m = sd['materials'][idx]
    if couple_gamma and field == 'gamma' and m.get('gamma_sat') is not None:
        delta = value - m['gamma']
        m['gamma_sat'] = m['gamma_sat'] + delta
    m[field] = value
    return sd


def _copy_for_edit(slope_data):
    """Copy slope_data one level deep plus fresh copies of every list-of-dicts
    the setters touch. Geometry objects (shapely) are shared — setters that
    change geometry (modify=) are responsible for rebuilding them. The finite-
    element ``mesh`` is shared by reference: seep/fem sweeps rebuild their data
    from the (copied) materials on the SAME mesh, and no setter mutates it."""
    sd = slope_data.copy()
    for key in ('materials', 'reinforcement_lines', 'reinforce_lines',
                'pile_lines', 'dloads'):
        if sd.get(key):
            sd[key] = [copy.deepcopy(it) for it in sd[key]]
    # Seepage BC dicts (specified_heads / exit_face / fluxes) are addressable by
    # the 'seep_bc' setter, so give each its own deep copy before an edit.
    for key in ('seepage_bc', 'seepage_bc2'):
        if sd.get(key):
            sd[key] = copy.deepcopy(sd[key])
    if sd.get('piezo_line'):
        sd['piezo_line'] = list(sd['piezo_line'])
    if sd.get('piezo_line2'):
        sd['piezo_line2'] = list(sd['piezo_line2'])
    return sd


def _dict_to_ref(d):
    """Normalize a plain-dict parameter spec to a (kind, name, field) tuple the
    rest of resolve_param understands. Accepts the shapes an LLM or GUI naturally
    writes:
        {'ref': 'mat:Clay:c'}                      -> passed straight through
        {'material': 'Clay', 'property': 'c'}      -> ('mat', 'Clay', 'c')
        {'material': 2, 'property': 'phi'}         -> ('mat', 2, 'phi')  (1-based index)
        {'global': 'k_seismic'}                    -> ('global', 'k_seismic')
        {'kind': 'seep', 'name': 'Clay', 'field': 'k1'}
        {'seep_bc': {'set': 1, 'head_index': 0}}   -> ('seep_bc', 1, 0)
    """
    d = {str(k).lower(): v for k, v in d.items()}
    if d.get('ref'):
        return d['ref']
    field = d.get('field') or d.get('property') or d.get('prop')
    if 'seep_bc' in d:
        spec = d['seep_bc']
        if isinstance(spec, dict):
            spec = {str(k).lower(): v for k, v in spec.items()}
            return ('seep_bc', spec.get('set', 1), spec.get('head_index', 0))
        # bare value addresses set 1, that head index
        return ('seep_bc', 1, spec)
    if 'global' in d:
        g = d['global']
        return ('global', field if g is True else g)
    kind = str(d.get('kind') or 'mat').lower()
    if kind == 'global':
        return ('global', field)
    name = d.get('name') or d.get('material') or d.get('mat') or d.get(kind)
    return (kind, name, field)


def resolve_param(slope_data, ref):
    """Validate a parameter reference against THIS model and resolve it.

    Returns (canonical_ref, setter, base_value) where setter is
    ``(slope_data, value) -> slope_data`` operating on a fresh copy.
    Raises ValueError naming what was given and what exists on any miss.
    """
    if isinstance(ref, dict):
        ref = _dict_to_ref(ref)
    if isinstance(ref, (tuple, list)):
        # keep native types (a material index stays an int for _find_by_name)
        parts = list(ref)
    else:
        parts = str(ref).split(':')
    if len(parts) == 2:
        kind, name, field = parts[0], '', parts[1]
    elif len(parts) == 3:
        kind, name, field = parts
    else:
        raise ValueError(f"Parameter ref '{ref}' is not 'kind:name:field' "
                         f"(or 'global:field' / 'geom:piezo:dy').")
    kind = str(kind).strip().lower()

    if kind == 'mat':
        idx, mat = _find_by_name(slope_data['materials'], name, 'material')
        mat_name = mat.get('name', f'Material_{idx+1}')
        allowed = _mat_strength_fields(mat, mat_name) + _MAT_GENERAL_FIELDS
        if field not in allowed:
            raise ValueError(
                f"Field '{field}' is not addressable on material '{mat_name}' "
                f"(option='{mat.get('option')}'). Allowed: {sorted(allowed)}")
        base = mat.get(field)
        if base is None:
            raise ValueError(f"Material '{mat_name}' has no value for '{field}'.")
        canonical = f"mat:{mat_name}:{field}"

        def setter(sd, value, _idx=idx, _field=field):
            sd = _copy_for_edit(sd)
            return _set_material_field(sd, _idx, _field, value)
        return canonical, setter, float(base)

    if kind == 'seep':
        idx, mat = _find_by_name(slope_data['materials'], name, 'material')
        mat_name = mat.get('name', f'Material_{idx+1}')
        if field not in _SEEP_FIELDS:
            raise ValueError(f"Field '{field}' is not a seep property. "
                             f"Allowed: {sorted(_SEEP_FIELDS)}")
        base = mat.get(field)
        if base is None:
            raise ValueError(f"Material '{mat_name}' has no value for '{field}'.")
        canonical = f"seep:{mat_name}:{field}"

        def setter(sd, value, _idx=idx, _field=field):
            sd = _copy_for_edit(sd)
            sd['materials'][_idx][_field] = value
            return sd
        return canonical, setter, float(base)

    if kind == 'seep_bc':
        # A specified-head boundary value: name = BC set (1 or 2), field = the
        # 0-based index into that set's 'specified_heads' list.
        try:
            bc_set = int(name)
            head_idx = int(field)
        except (TypeError, ValueError):
            raise ValueError(f"seep_bc ref must be 'seep_bc:<set>:<head_index>' "
                             f"(e.g. seep_bc:1:0); got set={name!r}, index={field!r}.")
        bc_key = 'seepage_bc2' if bc_set == 2 else 'seepage_bc'
        heads = (slope_data.get(bc_key) or {}).get('specified_heads') or []
        if not heads:
            raise ValueError(f"BC set {bc_set} has no specified-head boundaries "
                             f"('{bc_key}'.specified_heads is empty).")
        if not (0 <= head_idx < len(heads)):
            raise ValueError(f"seep_bc head index {head_idx} out of range "
                             f"(BC set {bc_set} has {len(heads)} specified head(s)).")
        base = heads[head_idx].get('head')
        if base is None:
            raise ValueError(f"BC set {bc_set} head #{head_idx} has no 'head' value.")
        canonical = f"seep_bc:{bc_set}:{head_idx}"

        def setter(sd, value, _key=bc_key, _i=head_idx):
            sd = _copy_for_edit(sd)
            sd[_key]['specified_heads'][_i]['head'] = value
            return sd
        return canonical, setter, float(base)

    if kind == 'reinforce':
        lines = slope_data.get('reinforcement_lines') or slope_data.get('reinforce_lines')
        if not lines:
            raise ValueError("The model has no reinforcement lines.")
        idx, line = _find_by_name(lines, name, 'reinforcement line', name_key='label')
        if field not in _REINFORCE_FIELDS:
            raise ValueError(f"Field '{field}' is not addressable on a reinforcement "
                             f"line. Allowed: {sorted(_REINFORCE_FIELDS)}")
        canonical = f"reinforce:{line.get('label')}:{field}"

        def setter(sd, value, _idx=idx, _field=field):
            sd = _copy_for_edit(sd)
            for key in ('reinforcement_lines', 'reinforce_lines'):
                if sd.get(key):
                    sd[key][_idx][_field] = value
            return sd
        return canonical, setter, float(line[field])

    if kind == 'piles':
        piles = slope_data.get('pile_lines')
        if not piles:
            raise ValueError("The model has no pile lines.")
        idx, pile = _find_by_name(piles, name, 'pile line', name_key='label')
        key = _PILE_FIELDS.get(field, field if field in _PILE_FIELDS.values() else None)
        if key is None:
            raise ValueError(f"Field '{field}' is not addressable on a pile line. "
                             f"Allowed: {sorted(_PILE_FIELDS)}")
        base = pile.get(key)
        if base is None:
            raise ValueError(f"Pile '{pile.get('label')}' has no value for '{field}' "
                             f"(it may be auto-computed — e.g. H=None means Ito & Matsui).")
        canonical = f"piles:{pile.get('label')}:{field}"

        def setter(sd, value, _idx=idx, _key=key):
            sd = _copy_for_edit(sd)
            sd['pile_lines'][_idx][_key] = value
            return sd
        return canonical, setter, float(base)

    if kind == 'global':
        if field not in _GLOBAL_FIELDS:
            raise ValueError(f"'{field}' is not a global parameter. "
                             f"Allowed: {sorted(_GLOBAL_FIELDS)}")
        base = slope_data.get(field, 0.0) or 0.0
        canonical = f"global:{field}"

        def setter(sd, value, _field=field):
            sd = _copy_for_edit(sd)
            sd[_field] = value
            return sd
        return canonical, setter, float(base)

    if kind == 'geom':
        # first named geometry transform: vertical water-table shift. The value
        # is a DELTA (dy), so the base value is 0 by construction.
        if (name, field) != ('piezo', 'dy'):
            raise ValueError(f"Unknown geometry transform 'geom:{name}:{field}'. "
                             f"Available: geom:piezo:dy (more arrive over time); "
                             f"for anything else write a modify= setter.")
        if not slope_data.get('piezo_line'):
            raise ValueError("geom:piezo:dy: the model has no piezometric line.")

        def setter(sd, value):
            sd = _copy_for_edit(sd)
            sd['piezo_line'] = [(x, y + value) for x, y in sd['piezo_line']]
            if sd.get('piezo_line2'):
                sd['piezo_line2'] = [(x, y + value) for x, y in sd['piezo_line2']]
            return sd
        return "geom:piezo:dy", setter, 0.0

    raise ValueError(f"Unknown parameter kind '{kind}'. "
                     f"Known: mat, reinforce, piles, global, seep, seep_bc, geom.")


def set_param(slope_data, ref, value):
    """Return a copy of slope_data with the referenced parameter set to value.
    The shared mutation path for sensitivity(), tornado(), and any back-analysis
    built on the same addressing."""
    _, setter, _ = resolve_param(slope_data, ref)
    return setter(slope_data, value)


def _resolve_sweep_spec(slope_data, param, modify, label):
    """Resolve the exclusive ``param`` / ``modify`` sweep contract into a single
    ``(canonical, setter, base_value)`` triple — the ONE code path a built-in
    reference and a user setter both travel, shared by sensitivity() and design()
    (so the docs' promise that "built-in references and modify= callables are one
    code path" holds for both, not by duplication but by this common resolver).

    Exactly one of:
      * ``param`` — a validated "kind:name:field" reference, resolved (and validated
        against THIS model) by resolve_param to a ``(slope_data, value) -> slope_data``
        setter with a numeric base value; or
      * ``modify`` — a user callable of that same signature, which requires ``label``
        (the swept axis's name); its base value is NaN (there is no stored scalar).

    Raises ValueError on any contract violation or invalid reference, so a caller
    returns a clean ``(False, message)`` rather than crashing."""
    if (param is None) == (modify is None):
        raise ValueError("Provide exactly one of param= or modify=.")
    if modify is not None and not label:
        raise ValueError("modify= requires label= (the df's param string).")
    if param is not None:
        return resolve_param(slope_data, param)
    return str(label), modify, float('nan')


# ---------------------------------------------------------------------------
# The sweep engine
# ---------------------------------------------------------------------------

# Which preflight analysis a sweep mode's per-step re-check is evaluated against.
_PREFLIGHT_BASE = {'lem': 'lem', 'fem': 'fem', 'seep': 'seep'}

# Above this the value in the 'fs' column is the search's no-admissible-surface
# sentinel (fs_fail = 9999), not a factor of safety. It is recorded like any other
# result and flows straight into plot_tornado, design()'s crossing interpolation and
# scaled_sensitivity's derivative, where it dominates every axis it touches.
_FS_SENTINEL = 100.0


def _geometry_baseline(slope_data):
    """Which material polygons are ALREADY invalid before any edit is made.

    ``_validate_model`` used to report every invalid polygon as "invalid after the
    edit", which on a model whose geometry was invalid to begin with blames the
    sweep for a defect it inherited -- and does so at every single point, so the
    whole study reads as a per-step failure. Seventy-five files in the corpus carry
    a ring shapely calls invalid and still solve correctly, so the honest question
    per step is not *is this polygon valid* but *did this edit break it*.
    """
    base = set()
    for i, p in enumerate(slope_data.get('polygons') or []):
        poly = p.get('polygon')
        if poly is not None and not poly.is_valid:
            base.add(i)
    return frozenset(base)


def _validate_model(sd, baseline=frozenset(), field=None, mode='lem',
                    selection=None):
    """Per-step precondition check of a modified model.

    Two halves, and both are needed:

    * **Geometry.** Setters -- user-written ``modify=`` callables above all -- are
      not trusted to leave the model coherent. A polygon the edit invalidated, or a
      ground surface it removed, is a failed point. A polygon that was invalid
      before the edit (``baseline``) is not: that is the model's own condition and
      the base-model check is where it belongs.
    * **Field ranges.** The substituted value itself is checked, against exactly
      the preflight rules that read the field being swept
      (:func:`xslope.preflight.rules_for_field`) -- a filtered subset, not the whole
      registry, so a tornado over 8 parameters x 9 points stays cheap. Without this
      a sweep that steps a unit weight to zero records the solver's 9999 sentinel
      as a success, and one that steps it negative records a negative factor of
      safety as a success.

    Returns the reason string for a point that must be skipped, or ``None``.
    """
    for i, p in enumerate(sd.get('polygons') or []):
        if i in baseline:
            continue
        poly = p.get('polygon')
        if poly is not None and not poly.is_valid:
            return f"material polygon (mat_id {p.get('mat_id')}) is invalid after the edit"
    gs = sd.get('ground_surface')
    if gs is None or gs.is_empty:
        return "model has no ground surface after the edit"
    if not field:
        return None
    from .preflight import preflight
    analysis = _PREFLIGHT_BASE.get(mode, 'lem')
    report = preflight(sd, analysis, selection, fields=[field])
    errs = report.errors
    return errs[0].message if errs else None


def _screen_point(rows, mode):
    """Post-step result screen: refuse to record a sentinel as a success.

    A point can fail without the solver saying so. ``circular_search`` scores an
    inadmissible trial ``fs_fail = 9999`` and returns it like any other answer, and
    a model with a negative strength solves to a negative factor of safety through
    the success path. Both are recorded with ``success=True`` today and both then
    dominate every plot and every derivative built from the sweep. This is where
    they become the ``success=False`` row the docstring already promises.

    Rows are screened in place and returned. Mode ``seep`` reports a discharge
    rather than a factor of safety, so only its non-finite values are screened --
    the sign of q is a direction, not an error.
    """
    for r in rows:
        if not r.get('success'):
            continue
        fs = r.get('fs')
        try:
            bad = fs is None or not np.isfinite(float(fs))
        except (TypeError, ValueError):
            bad = True
        if bad:
            r['success'] = False
            r['msg'] = r.get('msg') or 'the solve returned no usable value'
            continue
        if mode == 'seep':
            continue
        fs = float(fs)
        if fs >= _FS_SENTINEL:
            r['success'] = False
            r['msg'] = (f"factor of safety came back as {fs:g}, which is the "
                        f"search's no-admissible-surface sentinel rather than a "
                        f"result")
        elif fs < 0:
            r['success'] = False
            r['msg'] = (f"factor of safety came back as {fs:g}. A negative factor "
                        f"of safety is not a result -- the model at this value is "
                        f"non-physical")
    return rows


def _run_lem_point(sd, methods, search, num_slices, search_opts=None):
    """Evaluate one model. Returns list of dicts (one per method).

    ``search_opts`` is an optional dict of extra ``circular_search`` keyword
    arguments -- the search-WINDOW constraints (``entry_range``, ``exit_range``,
    ``center_box``, ``tangent_depth``, ``min_slip_depth``) above all. Callers
    build it with :func:`_file_search_window` so that a model carrying a window
    on its circles sheet is searched inside it here exactly as it is on Studio's
    Run LEM path.

    Only the circular branch takes the window, which is the same line Studio
    draws: ``circular_search`` is the function that accepts entry/exit/centre
    limits, and a non-circular search has no notion of them. The one limit
    ``noncircular_search`` does understand -- ``min_slip_depth`` -- is forwarded
    when it is there, and everything else is dropped rather than raising.
    """
    from .slice import generate_slices
    from .solve import solve_selected

    rows = []
    kw = dict(search_opts or {})
    circular = bool(sd.get('circular', True))
    if search:
        from .search import circular_search, noncircular_search
        for method in methods:
            try:
                if circular:
                    out = circular_search(sd, method, num_slices=num_slices, **kw)
                    best = out[0][0]
                    fs = best.get('FS')
                    R = (best['Yo'] - best['Depth']
                         if best.get('Yo') is not None and best.get('Depth') is not None
                         else np.nan)
                    rows.append({'method': method, 'fs': fs, 'success': fs is not None,
                                 'msg': '', 'Xo': best.get('Xo'), 'Yo': best.get('Yo'),
                                 'R': R})
                else:
                    nc_kw = {k: v for k, v in kw.items() if k == 'min_slip_depth'}
                    out = noncircular_search(sd, method, num_slices=num_slices,
                                             **nc_kw)
                    best = out[0][0]
                    fs = best.get('FS')
                    rows.append({'method': method, 'fs': fs, 'success': fs is not None,
                                 'msg': '', 'Xo': np.nan, 'Yo': np.nan, 'R': np.nan})
            except Exception as e:                        # noqa: BLE001
                rows.append({'method': method, 'fs': np.nan, 'success': False,
                             'msg': f'search failed: {e}', 'Xo': np.nan,
                             'Yo': np.nan, 'R': np.nan})
        return rows

    # fixed-surface evaluation (fast; right for "given this surface" questions)
    # check_inputs=False: a swept value is a deliberate perturbation of the model,
    # not a user mistake, and refusing here would abort the whole sweep. The swept
    # value IS validated -- once against the whole registry before the sweep starts,
    # and per step against the rules that read the swept field (_validate_model) --
    # so what reaches this call has already been screened, and screened into a
    # success=False ROW rather than an exception.
    try:
        if circular:
            circ = sd['circles'][0]
            ok, res = generate_slices(sd, circle=circ, num_slices=num_slices,
                                      check_inputs=False)
        else:
            ok, res = generate_slices(sd, non_circ=sd['non_circ'],
                                      num_slices=num_slices, check_inputs=False)
        if not ok:
            raise RuntimeError(res)
        df = res[0]
    except Exception as e:                                # noqa: BLE001
        return [{'method': m, 'fs': np.nan, 'success': False,
                 'msg': f'generate_slices failed: {e}',
                 'Xo': np.nan, 'Yo': np.nan, 'R': np.nan} for m in methods]
    for method in methods:
        # The solver call is inside the try for the same reason the slicing is: a
        # point that raises must become one failed ROW. Left bare, a TypeError here
        # escaped sensitivity() entirely and destroyed every point already computed
        # -- a sweep dying at value 7 of 9 returned nothing at all, contradicting
        # this function's own contract.
        try:
            r = solve_selected(method, df)
        except Exception as e:                            # noqa: BLE001
            rows.append({'method': method, 'fs': np.nan, 'success': False,
                         'msg': f'{method} failed: {e}', 'Xo': np.nan,
                         'Yo': np.nan, 'R': np.nan})
            continue
        if isinstance(r, str):
            rows.append({'method': method, 'fs': np.nan, 'success': False, 'msg': r,
                         'Xo': np.nan, 'Yo': np.nan, 'R': np.nan})
        else:
            rows.append({'method': method, 'fs': r['FS'], 'success': True, 'msg': '',
                         'Xo': (sd['circles'][0]['Xo'] if circular else np.nan),
                         'Yo': (sd['circles'][0]['Yo'] if circular else np.nan),
                         'R': (sd['circles'][0]['R'] if circular else np.nan)})
    return rows


# output quantity per engine mode: (short name, axis/label). LEM and FEM both
# report a factor of safety; a seepage sweep reports total discharge q.
_OUTPUT_BY_MODE = {
    'lem': ('FS', 'Factor of Safety'),
    'fem': ('FS', 'Factor of Safety'),
    'seep': ('q', 'Total discharge, q'),
}


def _fail_rows(method, msg):
    return [{'method': method, 'fs': np.nan, 'success': False, 'msg': msg,
             'Xo': np.nan, 'Yo': np.nan, 'R': np.nan}]


def _run_fem_point(sd, fem_opts, cancel_check=None):
    """Evaluate one model with the SSRM (xslope.fem). Output quantity is FS, so
    the row keeps the 'fs' column like the LEM path. The mesh must already live in
    ``sd['mesh']`` (a sweep runs every point on the SAME mesh — only the material
    fields change); build_fem_data re-derives the element properties from the
    edited materials each time."""
    from .fem import build_fem_data, solve_ssrm
    from .search import AnalysisCancelled

    mesh = sd.get('mesh')
    if mesh is None:
        return _fail_rows('ssrm', "FEM sweep needs a mesh in slope_data['mesh'] — "
                                  "build one first.")
    opts = fem_opts or {}
    try:
        fem_data = build_fem_data(sd, mesh)
        result = solve_ssrm(
            fem_data,
            F_min=opts.get('F_min', 1.0), F_max=opts.get('F_max', 2.0),
            tolerance=opts.get('tolerance', 0.01),
            failure_criterion=opts.get('failure_criterion', 'non_convergence'),
            min_slip_depth=opts.get('min_slip_depth'),
            debug_level=opts.get('debug_level', 0),
            # Sensitivity reads only result['FS'] per sweep point, never the
            # at-failure field — skip the extra non-converging capture solve.
            capture_failure_state=False,
            cancel_check=cancel_check)
    except AnalysisCancelled:
        raise                                          # cancel propagates to the sweep
    except Exception as e:                             # noqa: BLE001
        return _fail_rows('ssrm', f'SSRM failed: {e}')
    if not result.get('converged', False):
        return _fail_rows('ssrm', f"SSRM did not converge: "
                                  f"{result.get('error', 'unknown error')}")
    fs = result.get('FS')
    return [{'method': 'ssrm', 'fs': fs, 'success': fs is not None, 'msg': '',
             'Xo': np.nan, 'Yo': np.nan, 'R': np.nan}]


def _run_seep_point(sd, seep_opts):
    """Evaluate one model with the seepage solver (xslope.seep). Output quantity
    is the TOTAL DISCHARGE q — the same 'flowrate' the seep regression suite locks
    — carried in the 'fs' column so the DataFrame schema is mode-independent. The
    mesh must already live in ``sd['mesh']``."""
    from .seep import build_seep_data, run_seepage_analysis

    mesh = sd.get('mesh')
    if mesh is None:
        return _fail_rows('seep', "Seep sweep needs a mesh in slope_data['mesh'] — "
                                  "build one first.")
    opts = seep_opts or {}
    bc = opts.get('bc', 1)
    tol = opts.get('tol', 1e-4)
    try:
        seep_data = build_seep_data(mesh, sd, seep_bc=bc)
        solution = run_seepage_analysis(seep_data, tol=tol,
                                        max_iter=int(opts.get('max_iter', 400)))
    except Exception as e:                             # noqa: BLE001
        return _fail_rows('seep', f'seepage solve failed: {e}')
    if solution is None or not solution.get('converged', True):
        return _fail_rows('seep', 'seepage solution did not converge (q unreliable)')
    q = solution.get('flowrate')
    return [{'method': 'seep', 'fs': q, 'success': q is not None, 'msg': '',
             'Xo': np.nan, 'Yo': np.nan, 'R': np.nan}]


def _run_point(sd, mode, methods, search, num_slices, fem_opts, seep_opts,
               cancel_check=None, search_opts=None):
    """Evaluate one model for the active engine mode."""
    if mode == 'fem':
        return _run_fem_point(sd, fem_opts, cancel_check=cancel_check)
    if mode == 'seep':
        return _run_seep_point(sd, seep_opts)
    return _run_lem_point(sd, methods, search, num_slices, search_opts=search_opts)


def _sweep_selection(slope_data, canonical, values, mode, search):
    """What the sweep's rules are evaluated against: the run, not the file.

    A sweep is described entirely by the call (sensitivity.py:19-22), so every
    sensitivity rule reads this rather than the model. It is built once and used
    twice -- for the base-model preflight and for each per-step re-check -- so the
    two can never be asked different questions.
    """
    return {'param': canonical, 'values': list(values), 'mode': mode,
            'search': bool(search), 'base': _PREFLIGHT_BASE.get(mode, 'lem'),
            'surface': 'circular' if slope_data.get('circular', True)
                       else 'noncircular'}


def _sweep_gate(slope_data, selection, check_inputs=True):
    """The one full preflight of the BASE model, before any point is evaluated.

    Every point of a sweep is the same model with one number changed, so a defect
    in the base model is a defect in all of them -- and reporting it once at the
    door is both cheaper and far clearer than reporting it n times as a per-step
    failure. Returns the refusal message, or ``None``.
    """
    if not check_inputs:
        return None
    from .preflight import preflight
    report = preflight(slope_data, 'sensitivity', selection)
    errs = report.errors
    if not errs:
        return None
    if len(errs) == 1:
        return errs[0].message
    return (f"This model cannot be swept ({len(errs)} problem(s) found):\n"
            + "\n".join(f"  {i + 1}. {e.message}" for i, e in enumerate(errs)))


def sensitivity(slope_data, param=None, modify=None, label=None, values=None,
                rel_range=0.5, n=9, mode='lem', analysis=None, methods=('spencer',),
                search=True, num_slices=40, fem_opts=None, seep_opts=None,
                search_opts=None, use_file_window=True,
                debug_level=0, progress_callback=None, cancel_check=None,
                check_inputs=True):
    """Sweep one input; report FS (and the critical surface) per point.

    Parameters:
        slope_data: model dict (never modified).
        param: parameter reference "kind:name:field" — see resolve_param.
        modify: callable (slope_data, value) -> slope_data; exclusive with
            param; requires label. The callable owns geometry consistency
            (rebuild polygons/ground surface if it moves them — see
            main_design.rebuild_geometry) and any dependent-feature coupling.
        label: the df's param string when modify is used.
        values: iterable of swept values. Default: base*(1 +/- rel_range),
            n points (requires a nonzero base value).
        mode: which engine evaluates each point — 'lem' (default; limit
            equilibrium, output FS), 'fem' (each point is a full SSRM solve via
            xslope.fem, output FS — minutes per point, needs slope_data['mesh']),
            or 'seep' (each point rebuilds and solves the seepage problem, output
            = TOTAL DISCHARGE q; also needs slope_data['mesh']).
        analysis: deprecated alias for ``mode`` (kept for back-compatibility).
        methods: LEM method names, any subset of the seven (mode='lem' only).
        fem_opts: dict of SSRM knobs for mode='fem' (F_min, F_max, tolerance,
            failure_criterion, min_slip_depth) — mirrors solve_ssrm's defaults,
            except failure_criterion, which stays pinned at 'non_convergence'
            here (solve_ssrm's own default is now 'hybrid') so a sweep's shape is
            not silently redefined; pass it explicitly to opt in.
        seep_opts: dict for mode='seep' (bc set, tol) — mirrors the seep runner.
        search_opts: extra ``circular_search`` keyword arguments applied at every
            searched point — the search WINDOW above all (``entry_range`` /
            ``exit_range`` / ``center_box`` / ``tangent_depth`` /
            ``min_slip_depth``), which is what keeps a sweep on one mechanism
            instead of letting it jump families between values.
        use_file_window: fold the model's own circles-sheet search window
            (v19, ``slope_data['search_window']``) into ``search_opts``, exactly
            as Studio's Run LEM path does — so a windowed model gets the same
            surface family from a sweep as from a single run, and the curve
            reports how FS moves rather than which minimum was measured.
            Explicit ``search_opts`` win. Set False to search unconstrained
            regardless of what the file declares.
        search: re-search the critical surface per point (default — the
            critical surface MOVES as parameters change and a fixed-surface
            sweep silently understates sensitivity) vs re-solve the stored
            surface (~50x faster; right for prescribed-surface questions).
        num_slices: slices per evaluation.
        progress_callback: optional callable(done, total, label) invoked once per
            swept point (the base case is point 1), for a GUI progress bar.
        cancel_check: optional callable() -> bool; checked before every point and
            if it returns True an xslope.search.AnalysisCancelled is raised so a
            background runner can abort the sweep cleanly. Never set for a plain
            data-in/data-out call.
        check_inputs: run the input checks (default True). Validation happens in
            two stages: ONE full preflight of the base model before the first
            point, then a per-step re-check of only the rules that read the field
            being swept. A step whose value carries an ERROR is SKIPPED WITH A
            STATED REASON — a success=False row carrying the rule's own message —
            and never stops the sweep. Pass False to bypass both stages; the
            post-step result screen still runs, because a sentinel recorded as a
            success is a defect in the record, not a strictness setting.

    Returns:
        (success, result): result['df'] is a tidy long-format DataFrame (one
        row per value x method; the unmodified model is the flagged is_base
        row), plus 'param', 'base_value', 'runtime'. A failed point is a
        success=False ROW, not an exception — a sweep that dies at value 7 of
        9 should still report 1-6.
    """
    if analysis is not None:                           # deprecated alias for mode=
        mode = analysis
    if mode not in _OUTPUT_BY_MODE:
        return False, (f"mode='{mode}' is not recognised "
                       f"(choose 'lem', 'fem', or 'seep').")
    output, output_label = _OUTPUT_BY_MODE[mode]
    if mode in ('fem', 'seep') and slope_data.get('mesh') is None:
        return False, (f"mode='{mode}' needs a finite-element mesh in "
                       f"slope_data['mesh'] — build one before sweeping.")
    t0 = time.perf_counter()
    try:
        canonical, setter, base_value = _resolve_sweep_spec(
            slope_data, param, modify, label)
    except (ValueError, KeyError) as e:
        return False, str(e)

    if values is None:
        if not np.isfinite(base_value) or base_value == 0:
            return False, ("values= is required when the base value is 0 or undefined "
                           "(a relative range about it is meaningless).")
        values = np.linspace(base_value * (1 - rel_range),
                             base_value * (1 + rel_range), n)
    values = np.asarray(list(values), dtype=float)

    # Stage 1 of the two-stage contract: one full preflight of the BASE model,
    # analysis-appropriate for the engine this sweep will run. A model that cannot
    # be analysed at all is refused here, once, instead of failing at every point.
    selection = _sweep_selection(slope_data, canonical, values, mode, search)
    gate = _sweep_gate(slope_data, selection, check_inputs=check_inputs)
    if gate:
        return False, gate
    # The geometry the sweep INHERITED, so a per-step check blames only the edit.
    baseline = _geometry_baseline(slope_data)
    swept_field = selection['param'].split(':')[-1] if modify is None else None

    methods = (methods,) if isinstance(methods, str) else tuple(methods)
    rows = []

    def add_rows(value, is_base, point_rows):
        for pr in point_rows:
            rows.append({'param': canonical, 'value': value,
                         'rel': (value / base_value if base_value not in (0.0,)
                                 and np.isfinite(base_value) else np.nan),
                         'is_base': is_base, 'analysis': mode, **pr})

    total = 1 + len(values)
    done = [0]

    def _tick(label):
        done[0] += 1
        if progress_callback is not None:
            try:
                progress_callback(done[0], total, label)
            except Exception:                             # noqa: BLE001
                pass

    def _check_cancel():
        if cancel_check is not None and cancel_check():
            from .search import AnalysisCancelled
            raise AnalysisCancelled()

    # error-row method labels match what a successful point of this mode carries
    fail_methods = {'fem': ('ssrm',), 'seep': ('seep',)}.get(mode, methods)

    # The search window is read from the BASE model, once, and applied at every
    # point: a sweep varies a strength or a load, never the circles sheet, so the
    # window is a property of the study rather than of the step. Read here for the
    # same reason fs_vs_time reads it -- a searched sweep that ignores the file's
    # window can report a different surface FAMILY at each value, and the jump in
    # FS then reads as sensitivity to the parameter instead of a change in what
    # was measured.
    kw = dict(search_opts or {})
    if use_file_window:
        for k, v in _file_search_window(slope_data, already=kw).items():
            kw[k] = v

    def _point(sd):
        # Stage 3: the post-step result screen. A point that produced a sentinel or
        # a non-physical answer is recorded as a failure with its reason, never as
        # a success carrying the number.
        return _screen_point(_run_point(sd, mode, methods, search, num_slices,
                                        fem_opts, seep_opts,
                                        cancel_check=cancel_check,
                                        search_opts=kw), mode)

    # base case first: the unmodified model
    _check_cancel()
    add_rows(base_value, True, _point(slope_data))
    _tick(f"{canonical} (base)")

    for v in values:
        _check_cancel()
        try:
            sd = setter(_copy_for_edit(slope_data), v) if modify is not None \
                else setter(slope_data, v)
        except Exception as e:                            # noqa: BLE001
            add_rows(v, False, [{'method': m, 'fs': np.nan, 'success': False,
                                 'msg': f'setter failed: {e}', 'Xo': np.nan,
                                 'Yo': np.nan, 'R': np.nan} for m in fail_methods])
            _tick(f"{canonical} = {v:g}")
            continue
        # Stage 2: the per-step re-check, over ONLY the rules that read the swept
        # field. A step carrying an ERROR is skipped with the rule's own reason.
        err = _validate_model(sd, baseline,
                              field=swept_field if check_inputs else None,
                              mode=mode, selection=selection)
        if err:
            add_rows(v, False, [{'method': m, 'fs': np.nan, 'success': False,
                                 'msg': err, 'Xo': np.nan, 'Yo': np.nan,
                                 'R': np.nan} for m in fail_methods])
            _tick(f"{canonical} = {v:g}")
            continue
        if debug_level > 0:
            print(f"sensitivity: {canonical} = {v:g}")
        add_rows(v, False, _point(sd))
        _tick(f"{canonical} = {v:g}")

    df = pd.DataFrame(rows, columns=['param', 'value', 'rel', 'is_base', 'analysis',
                                     'method', 'fs', 'success', 'msg',
                                     'Xo', 'Yo', 'R'])
    # output-quantity metadata travels WITH the df (plots read it) and in the
    # result dict. 'fs' stays the value column for API compatibility even when it
    # carries q; 'output'/'output_label' say what it really is.
    df['output'] = output
    df['output_label'] = output_label
    return True, {'df': df, 'param': canonical, 'base_value': base_value,
                  'mode': mode, 'output': output, 'output_label': output_label,
                  'runtime': time.perf_counter() - t0}


def tornado(slope_data, params, rel_range=0.25, bounds=None, mode='lem',
            analysis=None, method='spencer', search=True, num_slices=40,
            fem_opts=None, seep_opts=None, search_opts=None,
            use_file_window=True):
    """Duncan-style tornado: the output quantity at the low/high bound of each parameter.

    Parameters:
        params: list of parameter refs.
        rel_range: default bounds base*(1 -/+ rel_range) per parameter.
        bounds: optional {ref: (low, high)} overriding rel_range per ref.
        mode: engine mode 'lem' (FS) / 'fem' (FS, SSRM) / 'seep' (discharge q).
        method: one LEM method (a tornado mixes parameters, not methods; mode='lem').
        search_opts / use_file_window: the search window every bar is measured
            inside — see :func:`sensitivity`. A tornado compares parameters
            against each other, so it matters more here than anywhere that every
            bar is measured on the same surface family.

    Returns (success, result): result['df'] has one row per (param, bound)
    plus the shared base row; result feeds plot_tornado. This is MFOSM's exact
    evaluation pattern with a range instead of sigma."""
    if analysis is not None:
        mode = analysis
    frames = []
    base_fs = None
    for i, ref in enumerate(params):
        try:
            canonical, _, base_value = resolve_param(slope_data, ref)
        except (ValueError, KeyError) as e:
            return False, str(e)
        if bounds and ref in bounds:
            lo, hi = bounds[ref]
        else:
            if base_value == 0:
                return False, (f"{canonical}: base value is 0 — give explicit "
                               f"bounds for this ref via bounds=.")
            lo, hi = base_value * (1 - rel_range), base_value * (1 + rel_range)
        ok, res = sensitivity(slope_data, param=ref, values=[lo, hi],
                              mode=mode, methods=(method,), search=search,
                              num_slices=num_slices, fem_opts=fem_opts,
                              seep_opts=seep_opts, search_opts=search_opts,
                              use_file_window=use_file_window)
        if not ok:
            return False, res
        df = res['df']
        if i == 0:
            base_fs = df.loc[df['is_base'], 'fs'].iloc[0]
            frames.append(df)                    # keep the base row once
        else:
            frames.append(df.loc[~df['is_base']])
    out = pd.concat(frames, ignore_index=True)
    output, output_label = _OUTPUT_BY_MODE.get(mode, _OUTPUT_BY_MODE['lem'])
    return True, {'df': out, 'base_fs': base_fs, 'method': method,
                  'mode': mode, 'output': output, 'output_label': output_label}


def tornado_from_sweeps(sweeps, base_fs=None, method=None):
    """Assemble a tornado result (the dict plot_tornado consumes) from full
    per-parameter sensitivity() sweeps.

    Where tornado() re-solves each parameter's two endpoints, this reuses sweeps a
    caller already ran — e.g. a GUI that runs a full FS-vs-value curve per
    parameter for click-through and wants the tornado for free. plot_tornado reads
    each parameter's lowest- and highest-value FS, so a full sweep yields the same
    bar as a 2-point one.

    Parameters:
        sweeps: dict {canonical_ref: df} or a sequence of sensitivity DataFrames,
            each carrying the base row plus one parameter's swept points.
        base_fs: base-case FS for the reference line; taken from the sweeps'
            is_base row if omitted.
        method: the LEM method label (for the plot's axis title).

    Returns {'df', 'base_fs', 'method'} — the shape tornado() returns.
    """
    dfs = list(sweeps.values()) if isinstance(sweeps, dict) else list(sweeps)
    if not dfs:
        return {'df': pd.DataFrame(), 'base_fs': base_fs, 'method': method}
    combined = pd.concat(dfs, ignore_index=True)
    if base_fs is None:
        b = combined.loc[combined['is_base'] & combined['success'], 'fs']
        base_fs = float(b.iloc[0]) if len(b) else None
    # carry the output-quantity labels the sweep dfs recorded, so plot_tornado
    # titles/axes match the swept quantity (FS or q)
    output = combined['output'].iloc[0] if 'output' in combined and len(combined) else 'FS'
    output_label = (combined['output_label'].iloc[0]
                    if 'output_label' in combined and len(combined)
                    else 'Factor of Safety')
    return {'df': combined, 'base_fs': base_fs, 'method': method,
            'output': output, 'output_label': output_label}


def _target_crossings(values, fs, target):
    """Parameter values where the FS curve crosses ``target``, linearly
    interpolated between adjacent solves. Handles a non-monotonic curve (returns
    every crossing) and endpoints exactly on target."""
    xs = []
    n = len(values)
    for i in range(n - 1):
        f0, f1 = fs[i], fs[i + 1]
        if not (np.isfinite(f0) and np.isfinite(f1)):
            continue
        g0, g1 = f0 - target, f1 - target
        if g0 == 0.0:
            xs.append(values[i])
        elif (g0 < 0) != (g1 < 0) and g1 != 0.0:
            t = g0 / (g0 - g1)
            xs.append(values[i] + t * (values[i + 1] - values[i]))
    if n and np.isfinite(fs[-1]) and (fs[-1] - target) == 0.0:
        xs.append(values[-1])
    uniq = []
    for x in xs:
        if not any(abs(x - u) < 1e-9 for u in uniq):
            uniq.append(float(x))
    return uniq


def design(slope_data, param=None, low=None, high=None, steps=11, target_fs=1.5,
           mode='lem', analysis=None, method='spencer', search=True, num_slices=40,
           fem_opts=None, seep_opts=None, search_opts=None, use_file_window=True,
           modify=None, label=None,
           progress_callback=None, cancel_check=None, debug_level=0,
           check_inputs=True):
    """Design sweep: vary ONE input from ``low`` to ``high`` and find the value at
    which the output quantity meets ``target_fs``.

    The deterministic-design staple — "vary the undrained strength between X and Y
    and find where FS = 1.5". Runs ``steps`` evenly spaced solves across
    [low, high] (re-searching the critical surface at each step by default), then
    linearly interpolates the parameter value where the output curve crosses
    ``target_fs``. Honest about misses: if the target is never reached inside the
    swept range, ``bracketed`` is False and ``extend``/``message`` say which way to
    widen the range. In mode='seep' the swept quantity is total discharge q and
    ``target_fs`` is the target q; the crossing logic is identical.

    The swept axis is named exactly the way sensitivity() names it — exclusively
    either ``param`` OR ``modify`` (they share one resolver, ``_resolve_sweep_spec``,
    so a built-in reference and a user setter are one code path here too):

        param: parameter reference — a "kind:name:field" string (e.g. "mat:Clay:c",
            "global:k_seismic"), a (kind, name, field) tuple (name may be a 1-based
            material index), or a dict {'material': name|index, 'property': field}
            / {'global': field}. See resolve_param for the full grammar.
        modify: a user callable ``(slope_data, value) -> slope_data`` (exclusive with
            param; requires ``label``) — the escape hatch for designing anything that
            is not a single stored scalar, geometry above all (a slope angle, a berm
            width; see main_design.set_slope_angle for the archetype). low/high/steps
            sweep the callable's VALUE axis identically, and the crossing/bracketing/
            extend semantics are unchanged; a callable that leaves the model invalid
            at a point becomes an honest ``success=False`` sweep point, never a
            silently inconsistent crossing.
        label: the swept-axis name when ``modify`` is used (the df's ``param`` string).
        low, high: inclusive bounds of the swept value (required).
        steps: number of solves across [low, high] (>= 2).
        target_fs: the design output value to locate (default 1.5 — a factor of
            safety in mode lem/fem, a discharge q in mode seep).
        mode: engine mode 'lem' (FS) / 'fem' (FS, SSRM) / 'seep' (discharge q).
        method: one LEM method name (mode='lem' only).
        fem_opts / seep_opts: engine knobs forwarded to sensitivity() for the
            'fem' / 'seep' modes (see sensitivity()).
        search_opts / use_file_window: the search window every step is measured
            inside — see :func:`sensitivity`. The crossing this function reports
            is only a design value if every step measured the same mechanism.
        search: re-search the critical surface at each step (default; correct — the
            critical surface moves as the parameter changes) or re-solve the stored
            surface (faster, but understates the movement). mode='lem' only.
        num_slices: slices per evaluation.
        progress_callback: optional callable(done, total, label) for a GUI.
        cancel_check: optional callable() -> bool; True mid-sweep raises
            xslope.search.AnalysisCancelled. Never set for a plain data call.

    Returns:
        (success, result). On success ``result`` carries:
          'df'         — the sensitivity DataFrame (output vs value, one method).
          'param'      — canonical parameter ref.
          'target_fs'  — the target output value.
          'crossing'   — interpolated parameter value at output = target, or None.
          'crossings'  — every crossing found (list; usually one).
          'bracketed'  — True iff the target is crossed inside [low, high].
          'fs_range'   — (min, max) output over the successful sweep points.
          'direction'  — 'increasing' / 'decreasing' / 'non-monotonic' trend of the
                         output vs the parameter (None if undetermined).
          'extend'     — when not bracketed, 'above {high}' or 'below {low}': which
                         way to widen the range to bracket the target; else None.
          'message'    — one-line human-readable summary.
          'output'/'output_label' — 'FS'/'q' and the axis label for the quantity.
          'base_value' — the parameter's current (unmodified) value.
          'runtime'    — wall-clock seconds.
    """
    if analysis is not None:
        mode = analysis
    if low is None or high is None:
        return False, "design() requires low= and high= (the swept-value bounds)."
    if int(steps) < 2:
        return False, "steps must be >= 2."
    values = np.linspace(float(low), float(high), int(steps))
    # param vs modify is validated by sensitivity() via the shared _resolve_sweep_spec,
    # so design and sensitivity enforce the identical exclusive contract (one path).
    ok, res = sensitivity(slope_data, param=param, modify=modify, label=label,
                          values=values, mode=mode,
                          methods=(method,), search=search, num_slices=num_slices,
                          fem_opts=fem_opts, seep_opts=seep_opts,
                          search_opts=search_opts, use_file_window=use_file_window,
                          progress_callback=progress_callback,
                          cancel_check=cancel_check, debug_level=debug_level,
                          check_inputs=check_inputs)
    if not ok:
        return False, res

    df = res['df']
    canonical = res['param']
    out_name = res.get('output', 'FS')                 # 'FS' or 'q'
    target = float(target_fs)
    swept = df.loc[~df['is_base'] & df['success']].sort_values('value')
    vals = swept['value'].to_numpy(dtype=float)
    fs = swept['fs'].to_numpy(dtype=float)
    crossings = _target_crossings(vals, fs, target)
    bracketed = len(crossings) > 0
    fs_min = float(np.min(fs)) if len(fs) else float('nan')
    fs_max = float(np.max(fs)) if len(fs) else float('nan')

    direction = None
    if len(fs) >= 2:
        diffs = np.diff(fs)
        if np.all(diffs >= -1e-9):
            direction = 'increasing'
        elif np.all(diffs <= 1e-9):
            direction = 'decreasing'
        else:
            direction = 'non-monotonic'

    extend = None
    if not bracketed and len(fs):
        rises = (direction == 'increasing') or (direction != 'decreasing'
                                                and fs[-1] >= fs[0])
        need_higher = target > fs_max
        if need_higher:
            extend = f"above {high:g}" if rises else f"below {low:g}"
        else:                                  # target below the whole swept range
            extend = f"below {low:g}" if rises else f"above {high:g}"

    if bracketed:
        crossing = float(crossings[0])
        message = (f"{out_name} = {target:g} at {canonical} = {crossing:.4g} "
                   f"(interpolated between solves).")
        if len(crossings) > 1:
            message += f"  [{len(crossings)} crossings; first reported]"
    else:
        crossing = None
        message = (f"{out_name} = {target:g} is not reached for {canonical} in "
                   f"[{low:g}, {high:g}] — {out_name} spans "
                   f"[{fs_min:.3g}, {fs_max:.3g}]."
                   + (f" Extend the range {extend} to bracket it." if extend else ""))

    return True, {'df': df, 'param': canonical, 'target_fs': target,
                  'crossing': crossing, 'crossings': [float(c) for c in crossings],
                  'bracketed': bracketed, 'fs_range': (fs_min, fs_max),
                  'direction': direction, 'extend': extend, 'message': message,
                  'mode': mode, 'output': res.get('output', 'FS'),
                  'output_label': res.get('output_label', 'Factor of Safety'),
                  'base_value': res['base_value'], 'runtime': res['runtime']}


def back_analysis(slope_data, param=None, low=None, high=None, steps=11, target_fs=1.0,
                  mode='lem', analysis=None, method='spencer', search=True,
                  num_slices=40, fem_opts=None, seep_opts=None, search_opts=None,
                  use_file_window=True, modify=None, label=None,
                  progress_callback=None, cancel_check=None, debug_level=0,
                  check_inputs=True):
    """Forensic back-analysis: find the parameter value that makes the slope
    limiting (FS = 1.0 by default).

    This is :func:`design` framed for a failure investigation. A slide has
    occurred, so the factor of safety at the time of failure is known to be 1.0;
    the unknown is a strength (or pore-pressure, or loading) parameter, and the
    back-analysis inverts the problem — *what value of this parameter is
    consistent with the observed failure?* The mechanics are identical to a design
    sweep; only the target and the interpretation differ, so this is a thin wrapper
    over :func:`design` with ``target_fs`` defaulting to 1.0.

    The returned ``crossing`` is the **back-calculated** value (e.g. the mobilized
    shear strength implied by the failure). ``bracketed`` False means the swept
    range never reaches FS = 1.0 — widen it per ``extend``. All other fields carry
    the same meaning as :func:`design`; ``result['study']`` is set to
    ``'back_analysis'`` so a caller can label the plot accordingly.

    Being a thin wrapper over :func:`design`, back-analysis **inherits the same
    exclusive** ``param`` **/** ``modify`` **contract** unchanged: pass a
    ``(slope_data, value) -> slope_data`` callable plus ``label`` to back-calculate
    anything that is not a stored scalar (a water-table elevation, say — sweep the
    phreatic surface and find the level consistent with the observed failure).
    """
    ok, res = design(slope_data, param=param, low=low, high=high, steps=steps,
                     target_fs=target_fs, mode=mode, analysis=analysis, method=method,
                     search=search, num_slices=num_slices, fem_opts=fem_opts,
                     seep_opts=seep_opts, search_opts=search_opts,
                     use_file_window=use_file_window, modify=modify, label=label,
                     progress_callback=progress_callback, cancel_check=cancel_check,
                     debug_level=debug_level, check_inputs=check_inputs)
    if not ok:
        return False, res
    res['study'] = 'back_analysis'
    out = res.get('output', 'FS')
    if res.get('bracketed') and res.get('crossing') is not None:
        res['message'] = (f"Back-analysis: {res['param']} = {res['crossing']:.4g} "
                          f"gives {out} = {target_fs:g} (the value consistent with "
                          f"the observed failure).")
    return True, res


# ---------------------------------------------------------------------------
# Factor of safety versus time
# ---------------------------------------------------------------------------

def _file_search_window(slope_data, already=()):
    """The circles-sheet search window (v19) as ``circular_search`` kwargs.

    The same reading Studio's Run-LEM path applies, so a windowed model gets the
    same surface family from a script as from the interface. A limit is forwarded
    only when the file supplies BOTH ends of it -- a half-declared range is not a
    window, and inventing the missing end would constrain the search in a
    direction nobody asked for. Anything the caller set (``already``) wins.
    """
    sw = (slope_data or {}).get('search_window') or {}
    if not sw:
        return {}
    kw = {}

    def pair(name, lo, hi):
        if name in already:
            return
        if sw.get(lo) is not None and sw.get(hi) is not None:
            kw[name] = (float(sw[lo]), float(sw[hi]))

    pair('entry_range', 'entry_x_min', 'entry_x_max')
    pair('exit_range', 'exit_x_min', 'exit_x_max')
    box = ('center_box_x_min', 'center_box_y_min',
           'center_box_x_max', 'center_box_y_max')
    if 'center_box' not in already and all(sw.get(k) is not None for k in box):
        kw['center_box'] = tuple(float(sw[k]) for k in box)
    if 'tangent_depth' not in already and sw.get('max_tangent_depth') is not None:
        gs = (slope_data or {}).get('ground_surface')
        try:
            y_top = max(y for _x, y in gs.coords)
        except Exception:                                  # noqa: BLE001
            y_top = None
        if y_top is not None and y_top > float(sw['max_tangent_depth']):
            kw['tangent_depth'] = (float(sw['max_tangent_depth']), float(y_top))
    if 'min_slip_depth' not in already and sw.get('min_slip_depth') is not None:
        kw['min_slip_depth'] = float(sw['min_slip_depth'])
    return kw


def _time_selection(slope_data, times, mode, search):
    """What an FS-vs-time run's preflight rules are evaluated against.

    ``param`` is None: no input is substituted at any point -- every step solves
    the SAME model against a different computed pore-pressure field -- so the
    swept-parameter rules have nothing to read and correctly stay silent. The
    engine mode, the surface family and the times are still declared, because the
    engine-appropriate rules do read those.
    """
    return {'param': None, 'values': list(times), 'times': list(times),
            'mode': mode, 'search': bool(search),
            'base': _PREFLIGHT_BASE.get(mode, 'lem'),
            'surface': 'circular' if slope_data.get('circular', True)
                       else 'noncircular'}


def fs_vs_time(slope_data, transient_solution, times=None, mode='lem',
               methods=('spencer',), search=True, num_slices=40, fem_opts=None,
               search_opts=None, use_file_window=True, seep_data=None,
               remarch=True, progress_callback=None,
               cancel_check=None, check_inputs=True, debug_level=0):
    """Factor of safety at every instant of a transient seepage solution.

    The coupled-analysis output the vendors publish: a transient march produces a
    sequence of pore-pressure fields, and this runs the chosen stability analysis
    against each one in turn and tabulates the result. It is a sweep like every
    other mode here -- run the model N times and report -- with one difference
    that makes it the cleanest instance of the pattern: **no input is modified at
    any step.** Each point solves the same model against a different computed
    field, so there is no substituted value to validate and no base case to
    compare against; the axis is time.

    The instant is never interpolated. A requested time that names no saved frame
    is served by re-marching with that instant injected into ``save_times``
    (:func:`xslope.seep.remarch_for_times`) when ``seep_data`` is supplied, and is
    otherwise a ``success=False`` row carrying that reason -- a field blended
    between two frames is not a solution of anything.

    Parameters:
        slope_data: model dict (never modified).
        transient_solution: the march to read, from
            ``xslope.seep.run_transient_seepage`` or
            ``import_transient_solution``.
        times: instants to evaluate. Default: every saved frame, which is the
            whole published curve. A list restricts (or extends -- see
            ``seep_data``) the set.
        mode: which engine evaluates each instant -- 'lem' (default) or 'fem'
            (a full SSRM solve per frame; needs ``slope_data['mesh']``). Mode
            'seep' is refused: the seepage solution is the INPUT here, not the
            output.
        methods: LEM method names, any subset of the seven (mode='lem' only).
            Every method is evaluated at every instant, so a multi-method run
            reports one curve per method from a single set of frames.
        search: re-search the critical surface at each instant (default -- and
            the right default here, because the critical surface MOVES as the
            pore pressures change, which is the whole phenomenon) or re-solve the
            stored surface.
        num_slices: slices per evaluation.
        fem_opts: SSRM knobs for mode='fem' (see :func:`sensitivity`).
        search_opts: extra ``circular_search`` keyword arguments -- the search
            WINDOW above all (``entry_range`` / ``exit_range`` / ``center_box`` /
            ``tangent_depth`` / ``min_slip_depth``), which is how a curve is kept
            on one mechanism instead of jumping families between instants.
        use_file_window: fold the model's own circles-sheet search window
            (v19, ``slope_data['search_window']``) into ``search_opts``, exactly
            as Studio's Run LEM path does. Explicit ``search_opts`` win. Set
            False to search unconstrained regardless of what the file declares.
        seep_data, remarch: with both supplied, times naming no saved frame are
            served by ONE re-march with all of them injected, before the first
            solve -- never one re-march per instant.
        progress_callback: optional callable(done, total, label), once per
            instant.
        cancel_check: optional callable() -> bool; True raises
            ``xslope.search.AnalysisCancelled``.
        check_inputs: run the input checks (default True). Same two-stage
            contract as :func:`sensitivity`: ONE full preflight of the base model
            for the engine this run will use, then a per-instant re-check. An
            instant that cannot be evaluated is a ``success=False`` ROW carrying
            its reason and never stops the run. The post-step result screen runs
            either way.

    Returns:
        (success, result). ``result['df']`` is the per-time table in the same
        tidy long format every other mode returns (one row per time x method,
        with ``value`` holding the time and ``param`` = 'time'), so
        ``xslope.plot.plot_sensitivity`` draws it unchanged. Plus:

          'times'          -- the instants evaluated, in order.
          'critical_time'  -- the time of the lowest factor of safety.
          'min_fs'         -- that lowest factor of safety.
          'critical'       -- {method: {'time', 'fs'}} when several methods ran.
          'remarched'      -- True if a re-march was needed to serve the times.
          'solution'       -- the (possibly re-marched) solution the frames came
                              from, so a caller keeps it instead of paying twice.
          'n_failed'       -- instants that produced no result.
          'mode' / 'output' / 'output_label' / 'runtime'.

    A rapid-drawdown analysis is a different question about the same problem and
    is deliberately not available here: this curve is a sequence of single-stage
    analyses against successive fields, where a Duncan-Wright-Brandon drawdown is
    one three-stage analysis reading two of them.
    """
    from .seep import (has_transient_frame, remarch_for_times,
                       select_transient_frame_u, transient_frame_index)
    from .water import with_water_loads

    if mode not in ('lem', 'fem'):
        return False, (f"mode={mode!r} is not available for an FS-vs-time run "
                       f"(choose 'lem' or 'fem'). A seepage sweep reports "
                       f"discharge, and the seepage solution is this run's INPUT.")
    if transient_solution is None:
        return False, ("fs_vs_time needs a transient seepage solution -- run or "
                       "import the march first (xslope.seep.run_transient_seepage "
                       "/ import_transient_solution).")
    if not (transient_solution.get('times') or []):
        return False, "the transient solution carries no saved frames."
    if mode == 'fem' and slope_data.get('mesh') is None:
        return False, ("mode='fem' needs a finite-element mesh in "
                       "slope_data['mesh'] -- build one before running.")
    output, output_label = _OUTPUT_BY_MODE[mode]
    t0 = time.perf_counter()

    if times is None:
        wanted = [float(t) for t in transient_solution['times']]
    else:
        wanted = sorted({float(t) for t in times})
    if not wanted:
        return False, "times= is empty; there is nothing to evaluate."

    # A requested instant with no frame is made into a computed state, once, for
    # the whole set -- never interpolated, and never one re-march per time.
    remarched = False
    missing = [t for t in wanted if not has_transient_frame(transient_solution, t)]
    if missing and remarch and seep_data is not None:
        try:
            transient_solution = remarch_for_times(
                seep_data, slope_data, missing)
            remarched = True
        except Exception as e:                             # noqa: BLE001
            return False, f"re-march for the requested times failed: {e}"

    methods = (methods,) if isinstance(methods, str) else tuple(methods)
    selection = _time_selection(slope_data, wanted, mode, search)
    # The base model is gated WITH a frame in it. Every point of this run carries
    # one, so a model whose materials read u = seep is not missing a pore-pressure
    # field -- gating the bare model would refuse the whole run for the one
    # condition the run itself removes at every step.
    gate_sd = _copy_for_edit(slope_data)
    try:
        select_transient_frame_u(gate_sd, transient_solution, time=wanted[0])
    except (ValueError, KeyError, IndexError):
        gate_sd = slope_data                   # no frame there: gate what we have
    gate = _sweep_gate(gate_sd, selection, check_inputs=check_inputs)
    if gate:
        return False, gate
    baseline = _geometry_baseline(slope_data)

    kw = dict(search_opts or {})
    if use_file_window:
        for k, v in _file_search_window(slope_data, already=kw).items():
            kw[k] = v

    fail_methods = ('ssrm',) if mode == 'fem' else methods
    rows = []
    total, done = len(wanted), [0]

    def _tick(label):
        done[0] += 1
        if progress_callback is not None:
            try:
                progress_callback(done[0], total, label)
            except Exception:                              # noqa: BLE001
                pass

    def _check_cancel():
        if cancel_check is not None and cancel_check():
            from .search import AnalysisCancelled
            raise AnalysisCancelled()

    def add_rows(t, point_rows):
        for pr in point_rows:
            rows.append({'param': 'time', 'value': float(t), 'rel': np.nan,
                         'is_base': False, 'analysis': mode, **pr})

    def fail(t, msg):
        add_rows(t, [{'method': m, 'fs': np.nan, 'success': False, 'msg': msg,
                      'Xo': np.nan, 'Yo': np.nan, 'R': np.nan}
                     for m in fail_methods])

    for t in wanted:
        _check_cancel()
        # The frame is placed on a COPY, and any water load derived for an earlier
        # instant is dropped with it: the pool pressing on the slope belongs to the
        # moment being solved, and a load carried over from t = 0 would load a
        # drained slope with a full reservoir.
        sd = _copy_for_edit(slope_data)
        for key in ('dloads_derived', 'dloads2_derived', 'water_derived'):
            sd.pop(key, None)
        try:
            transient_frame_index(transient_solution, t)
        except ValueError as e:
            fail(t, str(e))
            _tick(f"t = {t:g}")
            continue
        try:
            select_transient_frame_u(sd, transient_solution, time=t)
        except Exception as e:                             # noqa: BLE001
            fail(t, f"could not read the frame at t = {t:g}: {e}")
            _tick(f"t = {t:g}")
            continue
        # Derive the automatic water load ONCE per instant. The slicer derives it
        # too, defensively, but it does so per call -- and a searched point calls it
        # once per trial surface, which on a reservoir model is most of the run
        # time. Deriving here is the same load (with_water_loads is idempotent and
        # reads the frame's own time, just set above), computed a few hundred times
        # less often. In manual mode it returns the model by identity.
        sd = with_water_loads(sd)
        err = _validate_model(sd, baseline, field=None, mode=mode,
                              selection=selection)
        if err:
            fail(t, err)
            _tick(f"t = {t:g}")
            continue
        if debug_level > 0:
            print(f"fs_vs_time: t = {t:g}")
        add_rows(t, _screen_point(
            _run_point(sd, mode, methods, search, num_slices, fem_opts, None,
                       cancel_check=cancel_check, search_opts=kw), mode))
        _tick(f"t = {t:g}")

    df = pd.DataFrame(rows, columns=['param', 'value', 'rel', 'is_base', 'analysis',
                                     'method', 'fs', 'success', 'msg',
                                     'Xo', 'Yo', 'R'])
    df['output'] = output
    df['output_label'] = output_label

    critical = {}
    for m, g in df.groupby('method'):
        ok_rows = g.loc[g['success']]
        if len(ok_rows):
            i = ok_rows['fs'].astype(float).idxmin()
            critical[str(m)] = {'time': float(df.at[i, 'value']),
                                'fs': float(df.at[i, 'fs'])}
    lead = fail_methods[0] if fail_methods else None
    head = critical.get(str(lead)) or (next(iter(critical.values()))
                                       if critical else None)
    return True, {'df': df, 'param': 'time', 'times': wanted,
                  'critical_time': head['time'] if head else None,
                  'min_fs': head['fs'] if head else None,
                  'critical': critical, 'remarched': remarched,
                  'solution': transient_solution,
                  'n_failed': int((~df['success']).sum()),
                  'mode': mode, 'output': output, 'output_label': output_label,
                  'runtime': time.perf_counter() - t0}


# ---------------------------------------------------------------------------
# Scaled sensitivity, variance contribution, and rank correlation
#
# These three feed the Parametric-study plot family alongside the tornado and the
# per-parameter spider. The first is a LOCAL measure (a derivative at the base
# case); the last two are its statistical cousins that reuse the reliability
# machinery in advanced.py rather than reimplementing it.
# ---------------------------------------------------------------------------

def _sigma_by_ref(slope_data, mode):
    """Map every sweepable ref to its reliability sigma (or None), from list_params —
    so the per-sigma scaling and any sigma-gated plot know which refs carry one."""
    out = {}
    try:
        for e in list_params(slope_data, mode='seep' if mode == 'seep' else 'lem'):
            out[e['ref']] = e.get('sigma')
    except Exception:                                     # noqa: BLE001
        pass
    return out


def scaled_sensitivity(slope_data, params, method='spencer', search=True,
                       num_slices=40, mode='lem', rel_step=0.01, fem_opts=None,
                       seep_opts=None, search_opts=None, use_file_window=True,
                       progress_callback=None, cancel_check=None):
    """Local scaled-sensitivity coefficients — the vertical-bar cousin of the
    tornado, made comparable across parameters with different units.

    For each parameter *p* with base value *p0* and base output *F0*, the derivative
    dF/dp is estimated with a CENTRAL DIFFERENCE at +/- ``rel_step`` (relative,
    default 1%)::

        dF/dp  ~=  [ F(p0*(1+h)) - F(p0*(1-h)) ] / (2*h*p0),   h = rel_step

    and three scalings are reported so bars of different-unit parameters can be
    compared on one axis:

      * ``elasticity``  = (dF/dp) * (p0/F0)   — dimensionless; the percent change in
        F per percent change in p. The DEFAULT the plot shows.
      * ``per_1pct``    = (dF/dp) * (0.01*p0) — the change in F for a 1% change in p.
      * ``per_sigma``   = (dF/dp) * sigma_p   — the change in F for a one-sigma change
        in p; present only where the model carries sigma_p (reliability columns).

    A parameter whose base value is 0 or undefined has no relative step and is
    skipped (recorded in ``result['skipped']``) rather than dividing by zero.

    Returns (success, result). ``result['bars']`` is a list of dicts (one per
    analyzable parameter), each carrying every scaling, ``dfdp``, ``sign``
    (+1 if F rises with p, -1 if it falls), ``sigma`` and ``base_value``, sorted by
    ``|elasticity|`` descending. Also: ``fs_base``, ``rel_step``, ``method``,
    ``mode``, ``output``/``output_label``, and ``has_sigma``.
    """
    params = [params] if isinstance(params, (str, dict, tuple)) else list(params)
    sigma_by_ref = _sigma_by_ref(slope_data, mode)
    h = float(rel_step)
    bars, skipped, fs_base = [], [], None
    for ref in params:
        try:
            canonical, _, base_value = resolve_param(slope_data, ref)
        except (ValueError, KeyError) as e:
            return False, str(e)
        if not np.isfinite(base_value) or base_value == 0:
            skipped.append({'param': canonical, 'reason':
                            'base value is 0 or undefined — no relative step'})
            continue
        lo, hi = base_value * (1 - h), base_value * (1 + h)
        ok, res = sensitivity(slope_data, param=ref, values=[lo, hi], mode=mode,
                              methods=(method,), search=search, num_slices=num_slices,
                              fem_opts=fem_opts, seep_opts=seep_opts,
                              search_opts=search_opts,
                              use_file_window=use_file_window,
                              progress_callback=progress_callback,
                              cancel_check=cancel_check)
        if not ok:
            return False, res
        df = res['df']
        b = df.loc[df['is_base'] & df['success'], 'fs']
        f0 = float(b.iloc[0]) if len(b) else np.nan
        swept = df.loc[~df['is_base'] & df['success']].sort_values('value')
        if len(swept) < 2 or not np.isfinite(f0) or f0 == 0:
            skipped.append({'param': canonical, 'reason':
                            'base or perturbed solve failed — no derivative'})
            continue
        if fs_base is None:
            fs_base = f0
        f_minus = float(swept['fs'].iloc[0])
        f_plus = float(swept['fs'].iloc[-1])
        dfdp = (f_plus - f_minus) / (hi - lo)
        sigma = sigma_by_ref.get(canonical)
        bars.append({
            'param': canonical, 'base_value': base_value, 'fs_base': f0,
            'dfdp': dfdp, 'elasticity': dfdp * base_value / f0,
            'per_1pct': dfdp * 0.01 * base_value,
            'per_sigma': (dfdp * sigma if sigma else None),
            'sigma': (float(sigma) if sigma else None),
            'sign': 1 if dfdp >= 0 else -1,
        })
    bars.sort(key=lambda d: abs(d['elasticity']), reverse=True)
    output, output_label = _OUTPUT_BY_MODE.get(mode, _OUTPUT_BY_MODE['lem'])
    return True, {'bars': bars, 'skipped': skipped, 'fs_base': fs_base,
                  'method': method, 'mode': mode, 'rel_step': h,
                  'output': output, 'output_label': output_label,
                  'has_sigma': any(b.get('per_sigma') is not None for b in bars)}


def variance_contribution(slope_data, method='spencer', circular=True, rapid=False,
                          search=True, progress_callback=None, cancel_check=None,
                          fs_tol=None, tol=None, max_iter=None, composite=False,
                          seed='circles'):
    """Per-parameter contribution to Var(FS), reusing the Taylor-series reliability
    machinery in :func:`xslope.advanced.reliability` (not reimplemented).

    The TSPM writes sigma_F^2 = SUM_i (delta_F_i / 2)^2, where delta_F_i is the
    F+ minus F- swing when parameter *i* moves by +/- sigma. Each term is exactly
    that parameter's variance contribution (dF/dp * sigma_p)^2; this function runs
    the reliability analysis, extracts those terms, normalizes each to a percent of
    the total variance, sorts them descending and accumulates them — the input to a
    variance-contribution Pareto. Requires the model to carry standard deviations
    (reliability's own precondition).

    Returns (success, result). ``result['bars']`` is one dict per uncertain
    parameter (``label``, ``param``, ``variance``, ``pct`` = % of Var(FS),
    ``cumulative`` = running % after sorting, ``delta_F``, ``sigma``), plus
    ``sigma_F``, ``F_MLV``, ``COV_F`` and ``method``.
    """
    from .advanced import reliability
    ok, res = reliability(slope_data, method, rapid=rapid, circular=circular,
                          debug_level=-1, search=search,
                          progress_callback=progress_callback, cancel_check=cancel_check,
                          fs_tol=fs_tol, tol=tol, max_iter=max_iter,
                          composite=composite, seed=seed)
    if not ok:
        return False, res
    bars = []
    for p in res.get('param_info', []):
        dF = p.get('delta_F')
        if dF is None:
            continue
        bars.append({
            'param': f"mat:{p['material_name']}:{p['param']}",
            'material_id': p['material_id'], 'field': p['param'],
            'label': f"{p['material_name']} · {p['param']}",
            'variance': (dF / 2.0) ** 2, 'delta_F': dF, 'sigma': p['std'],
        })
    total = sum(b['variance'] for b in bars)
    for b in bars:
        b['pct'] = 100.0 * b['variance'] / total if total else 0.0
    bars.sort(key=lambda b: b['variance'], reverse=True)
    cum = 0.0
    for b in bars:
        cum += b['pct']
        b['cumulative'] = cum
    return True, {'bars': bars, 'sigma_F': res.get('sigma_F'),
                  'F_MLV': res.get('F_MLV'), 'COV_F': res.get('COV_F'),
                  'method': method}


def mc_rank_correlation(slope_data, method='bishop', n_samples=None, rng_seed=None,
                        distribution='normal', circular=True, rapid=False,
                        search=True, num_slices=40, progress_callback=None,
                        cancel_check=None):
    """Spearman rank correlation between each sampled input and FS, from a Monte
    Carlo reliability run (:func:`xslope.advanced.reliability_mc`).

    Where :func:`scaled_sensitivity` is a LOCAL measure (a derivative at the base
    case), this is a GLOBAL one: it ranks parameters by how strongly they drive FS
    across the whole sampled distribution, capturing nonlinearity and the parameter's
    own spread. The two answer different questions and are meant to be read together.

    A single reliability_mc campaign is run; its per-realization input matrix
    (``param_samples``) and FS vector (``fs_samples``) are correlated
    parameter-by-parameter with :func:`scipy.stats.spearmanr` over the admissible
    realizations. Requires standard deviations (reliability_mc's precondition) and,
    at 10^4 default samples, is the most expensive plot in the family.

    Returns (success, result). ``result['bars']`` is one dict per parameter
    (``label``, ``param``, ``rho`` in [-1, 1], ``sign``, ``sigma``), sorted by
    ``|rho|`` descending, plus ``method``, ``n_valid``, ``n_samples`` and
    ``distribution``.
    """
    from .advanced import reliability_mc
    kw = {}
    if n_samples is not None:
        kw['n_samples'] = int(n_samples)
    if rng_seed is not None:
        kw['rng_seed'] = int(rng_seed)
    ok, res = reliability_mc(slope_data, method, circular=circular, rapid=rapid,
                             distribution=distribution, search=search,
                             num_slices=num_slices, progress_callback=progress_callback,
                             cancel_check=cancel_check, debug_level=-1, **kw)
    if not ok:
        return False, res
    from scipy.stats import spearmanr
    samples = np.asarray(res.get('param_samples'), dtype=float)
    fs = np.asarray(res.get('fs_samples'), dtype=float)
    param_info = res.get('param_info', [])
    valid = np.isfinite(fs)
    fs_v = fs[valid]
    X = samples[valid] if samples.ndim == 2 else samples.reshape(len(fs), -1)[valid]
    bars = []
    for j, p in enumerate(param_info):
        col = X[:, j] if X.size else np.array([])
        if fs_v.size < 3 or col.size == 0 or float(np.std(col)) == 0.0:
            rho = float('nan')
        else:
            rho = float(spearmanr(col, fs_v).correlation)
        bars.append({
            'param': f"mat:{p['material_name']}:{p['param']}",
            'material_id': p['material_id'], 'field': p['param'],
            'label': f"{p['material_name']} · {p['param']}",
            'rho': rho, 'sign': 1 if (np.isfinite(rho) and rho >= 0) else -1,
            'sigma': p['std'],
        })
    bars.sort(key=lambda b: abs(b['rho']) if np.isfinite(b['rho']) else -1.0,
              reverse=True)
    return True, {'bars': bars, 'method': method, 'n_valid': res.get('n_valid'),
                  'n_samples': res.get('n_samples'),
                  'distribution': res.get('distribution')}


def _list_seep_params(slope_data):
    """The natural seepage sweep set: each material's hydraulic-conductivity and
    unsaturated fields, plus every specified-head boundary value (seep_bc refs).
    A seepage answer (discharge q) is driven by k and the boundary heads, so those
    are what a seep sensitivity/design study varies."""
    out = []
    for i, mat in enumerate(slope_data.get('materials') or []):
        name = mat.get('name') or f'Material_{i + 1}'
        for field in _SEEP_FIELDS:
            if field not in mat:
                continue
            val = mat.get(field)
            out.append({
                'ref': f'seep:{name}:{field}', 'kind': 'seep', 'name': name,
                'index': i + 1, 'field': field,
                'value': (float(val) if isinstance(val, (int, float))
                          and not isinstance(val, bool) else None),
                'sigma': None, 'label': f'{name} · {field}',
                'available': True, 'reason': '',
            })
    # Specified-head boundaries: one entry per head in each BC set present.
    for bc_set, key in ((1, 'seepage_bc'), (2, 'seepage_bc2')):
        heads = (slope_data.get(key) or {}).get('specified_heads') or []
        for j, bc in enumerate(heads):
            val = bc.get('head')
            out.append({
                'ref': f'seep_bc:{bc_set}:{j}', 'kind': 'seep_bc',
                'name': f'BC set {bc_set}', 'index': j,
                'field': f'head[{j}]',
                'value': (float(val) if isinstance(val, (int, float))
                          and not isinstance(val, bool) else None),
                'sigma': None, 'label': f'BC{bc_set} · head #{j} ({val:g})'
                if isinstance(val, (int, float)) else f'BC{bc_set} · head #{j}',
                'available': True, 'reason': '',
            })
    return out


def _unsweepable_reason(mat, mat_name, field, value):
    """Why this material field cannot be swept on this model, or ``''``.

    The same conditions the ``sensitivity.*`` preflight rules report, asked of one
    row so a picker can dim it. Offering a row and then refusing the run is the
    pattern the dim-with-reason design exists to remove.
    """
    if field == 'ru' and str(mat.get('u') or '').strip().lower() != 'ru':
        return (f"'{mat_name}' takes its pore pressure from u = "
                f"'{str(mat.get('u') or '').strip() or 'blank'}', so ru is never "
                f"read and sweeping it does not move the factor of safety.")
    if field in ('d', 'psi'):
        return (f"{field} is read only by the rapid drawdown analysis, which a "
                f"parametric study does not run.")
    if value is None:
        return (f"'{mat_name}' has no value for {field}, so there is no base value "
                f"to build a relative range around. Sweep it with explicit values=.")
    return ''


def list_params(slope_data, mode='lem'):
    """Enumerate every sweepable parameter in this model as plain dicts — the menu
    a GUI parameter-picker or an LLM uses to drive sensitivity()/design()/tornado().

    In the default mode ('lem'; 'fem' behaves the same), covers every material's
    numeric strength and general fields (option-aware, the same set resolve_param
    accepts) plus the global k_seismic. In mode='seep' the menu switches to the
    seepage-relevant set: each material's hydraulic fields (k1, k2, alpha, kr0, h0)
    plus every specified-head boundary value (seep_bc refs). Blank/zero-valued
    fields are still listed (value None/0) so a design sweep with explicit bounds
    can target them.

    Each entry:
        'ref'    — canonical "kind:name:field" string accepted by every entry point
        'kind'   — 'mat', 'seep', 'seep_bc', or 'global'
        'name'   — material / BC-set name (None for globals)
        'index'  — 1-based material index (None for globals)
        'field'  — property key
        'value'  — current value (float, or None if unset)
        'sigma'  — reliability std-dev for the field (sigma_c, sigma_phi, …) if the
                   model carries a non-zero one, else None — lets a GUI offer a
                   one-click +/-sigma range
        'label'  — short human label, e.g. 'Clay · c'
        'available' — False when this row cannot actually be swept on THIS model,
                   with 'reason' saying why. A picker should dim such a row rather
                   than offering it and refusing the run afterwards: the menu and
                   the run are answering the same question, so they must not
                   disagree. Two cases carry it — a field with no stored value (a
                   relative range about nothing), and a field the material's own
                   options make inert (ru where u is not 'ru'; d and psi, which
                   only a rapid drawdown analysis reads).
    """
    if mode == 'seep':
        return _list_seep_params(slope_data)
    out = []
    for i, mat in enumerate(slope_data.get('materials') or []):
        name = mat.get('name') or f'Material_{i + 1}'
        try:
            strength = _mat_strength_fields(mat, name)
        except Exception:                                 # noqa: BLE001
            strength = ()
        seen = set()
        for field in tuple(strength) + _MAT_GENERAL_FIELDS:
            if field in seen or field not in mat:
                continue
            seen.add(field)
            val = mat.get(field)
            sig = mat.get('sigma_' + field)
            value = (float(val) if isinstance(val, (int, float))
                     and not isinstance(val, bool) else None)
            reason = _unsweepable_reason(mat, name, field, value)
            entry = {
                'ref': f'mat:{name}:{field}', 'kind': 'mat', 'name': name,
                'index': i + 1, 'field': field,
                'value': value,
                'sigma': (float(sig) if isinstance(sig, (int, float))
                          and not isinstance(sig, bool) and sig else None),
                'label': f'{name} · {field}',
                'available': not reason, 'reason': reason,
            }
            out.append(entry)
    out.append({'ref': 'global:k_seismic', 'kind': 'global', 'name': None,
                'index': None, 'field': 'k_seismic',
                'value': float(slope_data.get('k_seismic') or 0.0),
                'sigma': None, 'label': 'k_seismic (global)',
                'available': True, 'reason': ''})
    return out
