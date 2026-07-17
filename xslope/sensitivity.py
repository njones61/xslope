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

__all__ = ['sensitivity', 'tornado', 'design', 'set_param', 'resolve_param',
           'list_params', 'tornado_from_sweeps']


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
               cancel_check=None):
    """Evaluate one model for the active engine mode."""
    if mode == 'fem':
        return _run_fem_point(sd, fem_opts, cancel_check=cancel_check)
    if mode == 'seep':
        return _run_seep_point(sd, seep_opts)
    return _run_lem_point(sd, methods, search, num_slices)


def sensitivity(slope_data, param=None, modify=None, label=None, values=None,
                rel_range=0.5, n=9, mode='lem', analysis=None, methods=('spencer',),
                search=True, num_slices=40, fem_opts=None, seep_opts=None,
                debug_level=0, progress_callback=None, cancel_check=None):
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
            failure_criterion, min_slip_depth) — mirrors solve_ssrm's defaults.
        seep_opts: dict for mode='seep' (bc set, tol) — mirrors the seep runner.
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

    def _point(sd):
        return _run_point(sd, mode, methods, search, num_slices,
                          fem_opts, seep_opts, cancel_check=cancel_check)

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
        err = _validate_model(sd)
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
            fem_opts=None, seep_opts=None):
    """Duncan-style tornado: the output quantity at the low/high bound of each parameter.

    Parameters:
        params: list of parameter refs.
        rel_range: default bounds base*(1 -/+ rel_range) per parameter.
        bounds: optional {ref: (low, high)} overriding rel_range per ref.
        mode: engine mode 'lem' (FS) / 'fem' (FS, SSRM) / 'seep' (discharge q).
        method: one LEM method (a tornado mixes parameters, not methods; mode='lem').

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
                              seep_opts=seep_opts)
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


def design(slope_data, param, low, high, steps=11, target_fs=1.5,
           mode='lem', analysis=None, method='spencer', search=True, num_slices=40,
           fem_opts=None, seep_opts=None,
           progress_callback=None, cancel_check=None, debug_level=0):
    """Design sweep: vary ONE parameter from ``low`` to ``high`` and find the value
    at which the output quantity meets ``target_fs``.

    The deterministic-design staple — "vary the undrained strength between X and Y
    and find where FS = 1.5". Runs ``steps`` evenly spaced solves across
    [low, high] (re-searching the critical surface at each step by default), then
    linearly interpolates the parameter value where the output curve crosses
    ``target_fs``. Honest about misses: if the target is never reached inside the
    swept range, ``bracketed`` is False and ``extend``/``message`` say which way to
    widen the range. In mode='seep' the swept quantity is total discharge q and
    ``target_fs`` is the target q; the crossing logic is identical.

    Parameters:
        param: parameter reference — a "kind:name:field" string (e.g. "mat:Clay:c",
            "global:k_seismic"), a (kind, name, field) tuple (name may be a 1-based
            material index), or a dict {'material': name|index, 'property': field}
            / {'global': field}. See resolve_param for the full grammar.
        low, high: inclusive bounds of the swept parameter value.
        steps: number of solves across [low, high] (>= 2).
        target_fs: the design output value to locate (default 1.5 — a factor of
            safety in mode lem/fem, a discharge q in mode seep).
        mode: engine mode 'lem' (FS) / 'fem' (FS, SSRM) / 'seep' (discharge q).
        method: one LEM method name (mode='lem' only).
        fem_opts / seep_opts: engine knobs forwarded to sensitivity() for the
            'fem' / 'seep' modes (see sensitivity()).
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
    if int(steps) < 2:
        return False, "steps must be >= 2."
    values = np.linspace(float(low), float(high), int(steps))
    ok, res = sensitivity(slope_data, param=param, values=values, mode=mode,
                          methods=(method,), search=search, num_slices=num_slices,
                          fem_opts=fem_opts, seep_opts=seep_opts,
                          progress_callback=progress_callback,
                          cancel_check=cancel_check, debug_level=debug_level)
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
            })
    return out


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
            out.append({
                'ref': f'mat:{name}:{field}', 'kind': 'mat', 'name': name,
                'index': i + 1, 'field': field,
                'value': (float(val) if isinstance(val, (int, float))
                          and not isinstance(val, bool) else None),
                'sigma': (float(sig) if isinstance(sig, (int, float))
                          and not isinstance(sig, bool) and sig else None),
                'label': f'{name} · {field}',
            })
    out.append({'ref': 'global:k_seismic', 'kind': 'global', 'name': None,
                'index': None, 'field': 'k_seismic',
                'value': float(slope_data.get('k_seismic') or 0.0),
                'sigma': None, 'label': 'k_seismic (global)'})
    return out
