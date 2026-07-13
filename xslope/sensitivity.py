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

__all__ = ['sensitivity', 'tornado', 'set_param', 'resolve_param']


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
    same class of silent error as u='peizo'). Power-curve coefficients are
    addressable HERE and deliberately not in reliability(): they are correlated
    fit coefficients, and sensitivity is the designated tool for those."""
    if material.get('option') == 'pow':
        return ('pow_a', 'pow_b', 'pow_c', 'pow_d')
    from .advanced import _strength_param_mapping
    mapping = _strength_param_mapping(material, mat_name)   # raises on unknown option
    return tuple(k for k in mapping.keys() if k != 'gamma')


def _find_by_name(items, name, kind, name_key='name'):
    """Case-insensitive unique lookup; raises naming what was given and what exists."""
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
    change geometry (modify=) are responsible for rebuilding them."""
    sd = slope_data.copy()
    for key in ('materials', 'reinforcement_lines', 'reinforce_lines',
                'pile_lines', 'dloads'):
        if sd.get(key):
            sd[key] = [copy.deepcopy(it) for it in sd[key]]
    if sd.get('piezo_line'):
        sd['piezo_line'] = list(sd['piezo_line'])
    if sd.get('piezo_line2'):
        sd['piezo_line2'] = list(sd['piezo_line2'])
    return sd


def resolve_param(slope_data, ref):
    """Validate a parameter reference against THIS model and resolve it.

    Returns (canonical_ref, setter, base_value) where setter is
    ``(slope_data, value) -> slope_data`` operating on a fresh copy.
    Raises ValueError naming what was given and what exists on any miss.
    """
    if isinstance(ref, (tuple, list)):
        parts = [str(p) for p in ref]
    else:
        parts = str(ref).split(':')
    if len(parts) == 2:
        kind, name, field = parts[0], '', parts[1]
    elif len(parts) == 3:
        kind, name, field = parts
    else:
        raise ValueError(f"Parameter ref '{ref}' is not 'kind:name:field' "
                         f"(or 'global:field' / 'geom:piezo:dy').")
    kind = kind.strip().lower()

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
                     f"Known: mat, reinforce, piles, global, seep, geom.")


def set_param(slope_data, ref, value):
    """Return a copy of slope_data with the referenced parameter set to value.
    The shared mutation path for sensitivity(), tornado(), and any back-analysis
    built on the same addressing."""
    _, setter, _ = resolve_param(slope_data, ref)
    return setter(slope_data, value)


# ---------------------------------------------------------------------------
# The sweep engine
# ---------------------------------------------------------------------------

def _validate_model(sd):
    """Engine-side sanity check of a modified model. Setters (especially
    user-written modify= callables) are not trusted to keep the model
    coherent; a failure here becomes a success=False row, not a crash."""
    if sd.get('polygons'):
        for p in sd['polygons']:
            poly = p.get('polygon')
            if poly is not None and not poly.is_valid:
                return f"material polygon (mat_id {p.get('mat_id')}) is invalid after the edit"
    gs = sd.get('ground_surface')
    if gs is None or gs.is_empty:
        return "model has no ground surface after the edit"
    return None


def _run_lem_point(sd, methods, search, num_slices):
    """Evaluate one model. Returns list of dicts (one per method)."""
    from .slice import generate_slices
    from .solve import solve_selected

    rows = []
    circular = bool(sd.get('circular', True))
    if search:
        from .search import circular_search, noncircular_search
        for method in methods:
            try:
                if circular:
                    out = circular_search(sd, method, num_slices=num_slices)
                    best = out[0][0]
                    fs = best.get('FS')
                    R = (best['Yo'] - best['Depth']
                         if best.get('Yo') is not None and best.get('Depth') is not None
                         else np.nan)
                    rows.append({'method': method, 'fs': fs, 'success': fs is not None,
                                 'msg': '', 'Xo': best.get('Xo'), 'Yo': best.get('Yo'),
                                 'R': R})
                else:
                    out = noncircular_search(sd, method, num_slices=num_slices)
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
    try:
        if circular:
            circ = sd['circles'][0]
            ok, res = generate_slices(sd, circle=circ, num_slices=num_slices)
        else:
            ok, res = generate_slices(sd, non_circ=sd['non_circ'], num_slices=num_slices)
        if not ok:
            raise RuntimeError(res)
        df = res[0]
    except Exception as e:                                # noqa: BLE001
        return [{'method': m, 'fs': np.nan, 'success': False,
                 'msg': f'generate_slices failed: {e}',
                 'Xo': np.nan, 'Yo': np.nan, 'R': np.nan} for m in methods]
    for method in methods:
        r = solve_selected(method, df)
        if isinstance(r, str):
            rows.append({'method': method, 'fs': np.nan, 'success': False, 'msg': r,
                         'Xo': np.nan, 'Yo': np.nan, 'R': np.nan})
        else:
            rows.append({'method': method, 'fs': r['FS'], 'success': True, 'msg': '',
                         'Xo': (sd['circles'][0]['Xo'] if circular else np.nan),
                         'Yo': (sd['circles'][0]['Yo'] if circular else np.nan),
                         'R': (sd['circles'][0]['R'] if circular else np.nan)})
    return rows


def sensitivity(slope_data, param=None, modify=None, label=None, values=None,
                rel_range=0.5, n=9, analysis='lem', methods=('spencer',),
                search=True, num_slices=40, debug_level=0):
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
        analysis: 'lem' (Phase 1; 'fem' and 'seep' arrive later).
        methods: LEM method names, any subset of the seven.
        search: re-search the critical surface per point (default — the
            critical surface MOVES as parameters change and a fixed-surface
            sweep silently understates sensitivity) vs re-solve the stored
            surface (~50x faster; right for prescribed-surface questions).
        num_slices: slices per evaluation.

    Returns:
        (success, result): result['df'] is a tidy long-format DataFrame (one
        row per value x method; the unmodified model is the flagged is_base
        row), plus 'param', 'base_value', 'runtime'. A failed point is a
        success=False ROW, not an exception — a sweep that dies at value 7 of
        9 should still report 1-6.
    """
    if analysis != 'lem':
        return False, f"analysis='{analysis}' is not implemented yet (Phase 1 is LEM)."
    if (param is None) == (modify is None):
        return False, "Provide exactly one of param= or modify=."
    if modify is not None and not label:
        return False, "modify= requires label= (the df's param string)."

    t0 = time.perf_counter()
    if param is not None:
        try:
            canonical, setter, base_value = resolve_param(slope_data, param)
        except (ValueError, KeyError) as e:
            return False, str(e)
    else:
        canonical, setter, base_value = str(label), modify, np.nan

    if values is None:
        if not np.isfinite(base_value) or base_value == 0:
            return False, ("values= is required when the base value is 0 or undefined "
                           "(a relative range about it is meaningless).")
        values = np.linspace(base_value * (1 - rel_range),
                             base_value * (1 + rel_range), n)
    values = np.asarray(list(values), dtype=float)

    methods = (methods,) if isinstance(methods, str) else tuple(methods)
    rows = []

    def add_rows(value, is_base, point_rows):
        for pr in point_rows:
            rows.append({'param': canonical, 'value': value,
                         'rel': (value / base_value if base_value not in (0.0,)
                                 and np.isfinite(base_value) else np.nan),
                         'is_base': is_base, 'analysis': analysis, **pr})

    # base case first: the unmodified model
    add_rows(base_value, True, _run_lem_point(slope_data, methods, search, num_slices))

    for v in values:
        try:
            sd = setter(_copy_for_edit(slope_data), v) if modify is not None \
                else setter(slope_data, v)
        except Exception as e:                            # noqa: BLE001
            add_rows(v, False, [{'method': m, 'fs': np.nan, 'success': False,
                                 'msg': f'setter failed: {e}', 'Xo': np.nan,
                                 'Yo': np.nan, 'R': np.nan} for m in methods])
            continue
        err = _validate_model(sd)
        if err:
            add_rows(v, False, [{'method': m, 'fs': np.nan, 'success': False,
                                 'msg': err, 'Xo': np.nan, 'Yo': np.nan,
                                 'R': np.nan} for m in methods])
            continue
        if debug_level > 0:
            print(f"sensitivity: {canonical} = {v:g}")
        add_rows(v, False, _run_lem_point(sd, methods, search, num_slices))

    df = pd.DataFrame(rows, columns=['param', 'value', 'rel', 'is_base', 'analysis',
                                     'method', 'fs', 'success', 'msg',
                                     'Xo', 'Yo', 'R'])
    return True, {'df': df, 'param': canonical, 'base_value': base_value,
                  'runtime': time.perf_counter() - t0}


def tornado(slope_data, params, rel_range=0.25, bounds=None, analysis='lem',
            method='spencer', search=True, num_slices=40):
    """Duncan-style tornado: FS at the low/high bound of each parameter.

    Parameters:
        params: list of parameter refs.
        rel_range: default bounds base*(1 -/+ rel_range) per parameter.
        bounds: optional {ref: (low, high)} overriding rel_range per ref.
        method: one LEM method (a tornado mixes parameters, not methods).

    Returns (success, result): result['df'] has one row per (param, bound)
    plus the shared base row; result feeds plot_tornado. This is MFOSM's exact
    evaluation pattern with a range instead of sigma."""
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
                              analysis=analysis, methods=(method,),
                              search=search, num_slices=num_slices)
        if not ok:
            return False, res
        df = res['df']
        if i == 0:
            base_fs = df.loc[df['is_base'], 'fs'].iloc[0]
            frames.append(df)                    # keep the base row once
        else:
            frames.append(df.loc[~df['is_base']])
    out = pd.concat(frames, ignore_index=True)
    return True, {'df': out, 'base_fs': base_fs, 'method': method}
