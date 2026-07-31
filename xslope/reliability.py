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

"""The reliability analysis family (Taylor Series Probability Method and Monte
Carlo), plus the variance-decomposition helpers the Parametric study reuses.

Public entry points:

* :func:`reliability` -- the front door; ``engine='taylor'`` (default) or ``'mc'``.
* :func:`reliability_taylor` -- Taylor Series Probability Method (limit equilibrium).
* :func:`reliability_mc` -- Monte Carlo reliability on a fixed surface.
* :func:`reliability_fem` -- Taylor series driven by the finite-element SSRM.

Historically this code lived in :mod:`xslope.advanced`, which now re-exports every
name here for backward compatibility, so ``from xslope.advanced import reliability``
(and ``reliability_mc`` / ``reliability_fem`` / the private helpers) keeps working.
"""

import time

import numpy as np
import pandas as pd
from scipy.stats import norm
from tabulate import tabulate

# The refusal messages live in preflight, next to the rules that report the same
# conditions BEFORE a run starts, so the sentence a dialog shows beforehand and the
# sentence an engine refuses with afterwards are one string with one definition.
from .preflight import (elastic_sigma_message, no_sigmas_message,
                        sigma_exceeds_mean_message, unsupported_option_message)

#: A factor of safety at or above this is the search's no-admissible-surface
#: sentinel (``fs_fail = 9999``) or an equally non-physical answer, never a real
#: result. Recording one as a success puts the sentinel straight into sigma_F.
FS_SENTINEL = 100.0

#: Below this the coefficient of variation of F has rounded out of the answer
#: entirely: beta comes back infinite and the run prints "Reliability: 100.00%".
#: The engines used to test for EXACT zero, which this does not reach.
COV_FLOOR = 1e-9


def reliability(slope_data, method='bishop', *args, engine='taylor', **kwargs):
    """Front door to the reliability family: pick the analysis engine, then hand off.

    This is the one public entry point Studio, the notebooks and the regression
    suite call. ``engine`` selects which reliability method runs and every other
    argument is forwarded to that engine unchanged, so a plain
    ``reliability(slope_data, method, ...)`` call is exactly equivalent to
    ``reliability_taylor(slope_data, method, ...)`` — the default — and all the
    established Taylor-series benchmarks run through the unchanged path.

    Parameters
    ----------
    slope_data : dict
        The usual slope-data dictionary.
    method : str
        The limit-equilibrium solver name ('oms', 'bishop', 'janbu', 'corps',
        'lowe', 'spencer', 'mprice'). Both engines take it as their first analysis
        argument. (Passed as the second positional/keyword, exactly as before.)
    engine : {'taylor', 'mc'}, keyword-only
        Which reliability engine to run:

        * ``'taylor'`` (default) -> :func:`reliability_taylor`, the Taylor Series
          Probability Method (1 + 2N limit-equilibrium searches).
        * ``'mc'`` -> :func:`reliability_mc`, a Monte Carlo sampling campaign on a
          fixed failure surface.

        Aliases are accepted case-insensitively ('tspm' for Taylor;
        'monte_carlo' / 'montecarlo' for MC).
    *args, **kwargs
        Forwarded verbatim to the selected engine. See each engine's signature for
        the accepted options (they overlap but are not identical — e.g. Monte Carlo
        adds ``n_samples`` / ``rng_seed`` / ``distribution``).

    Returns
    -------
    (success, result) : tuple
        Whatever the selected engine returns.

    See Also
    --------
    reliability_taylor : Taylor Series Probability Method (limit equilibrium).
    reliability_mc : Monte Carlo reliability on a fixed surface (limit equilibrium).
    reliability_fem : Taylor series with the finite-element SSRM factor of safety.
    """
    key = str(engine).lower().replace('-', '_')
    if key in ('taylor', 'tspm', 'taylor_series'):
        return reliability_taylor(slope_data, method, *args, **kwargs)
    if key in ('mc', 'monte_carlo', 'montecarlo'):
        return reliability_mc(slope_data, method, *args, **kwargs)
    raise ValueError(
        f"reliability(): unknown engine={engine!r}. Use 'taylor' (default, the "
        f"Taylor Series Probability Method) or 'mc' (Monte Carlo).")


def _reliability_gate(slope_data, selection, base='lem', check_inputs=True):
    """The full base-model preflight a reliability run starts with.

    A reliability analysis is 1 + 2N solves, and every one of them uses the same
    model: a defect in it is a defect in all of them. Checking once at the door
    costs one pass over the registry and saves the whole campaign -- the Taylor
    path used to discover an unsupported strength option at the perturbation step,
    *after* the critical-surface search, and to raise rather than return, breaking
    the ``(success, result)`` contract its callers handle.

    Returns a message string when the run must be refused, else ``None``. The
    reliability-specific rules and every rule of the base analysis are evaluated
    together, so what is reported here is what a Run dialog would already have
    dimmed the button over.
    """
    if not check_inputs:
        return None
    from .preflight import preflight
    sel = dict(selection or {})
    sel['base'] = base
    report = preflight(slope_data, 'reliability', sel)
    errs = report.errors
    if not errs:
        return None
    if len(errs) == 1:
        return errs[0].message
    return ("This model cannot be run as a reliability analysis "
            f"({len(errs)} problem(s) found):\n"
            + "\n".join(f"  {i + 1}. {e.message}" for i, e in enumerate(errs)))


def _cov_floor_reason(cov):
    """Why a coefficient of variation this small is no variability at all.

    The engines used to test ``COV_F == 0`` exactly, so a standard deviation that
    moved the factor of safety by less than its own last digit passed the test and
    came out as beta = infinity, printed as "Reliability: 100.00%" -- the most
    confident answer the program can give, produced by having nothing to say.
    """
    return (f"Reliability: the coefficient of variation of the factor of safety is "
            f"{cov:.3g}, which is indistinguishable from zero. The standard "
            f"deviations provided do not move the factor of safety, so the "
            f"reliability index is infinite and the probability of failure would be "
            f"reported as exactly zero. Check that the sigma columns are set on the "
            f"materials the failure surface actually crosses.")


def _sentinel_reason(label, value):
    """Why a factor of safety this large is a sentinel rather than an answer.

    ``circular_search`` scores an inadmissible trial ``fs_fail = 9999`` and returns
    it like any other result. Fed into the Taylor combination it dominates sigma_F
    completely, and the reliability index that comes out is arithmetic performed on
    a flag. Nothing downstream filters it, so this is where it stops.
    """
    return (f"Reliability: {label} came back as {value:.4g}, which is the search's "
            f"no-admissible-surface sentinel rather than a factor of safety. The "
            f"perturbed model has no analyzable failure surface, so the reliability "
            f"index cannot be formed from it.")


def reliability_taylor(slope_data, method, rapid=False, circular=True, debug_level=0,
                progress_callback=None, cancel_check=None,
                fs_tol=None, tol=None, max_iter=None, composite=False, seed='circles',
                search=True, search_opts=None, use_file_window=True,
                check_inputs=True):
    """
    Performs reliability analysis using the Taylor Series Probability Method (TSPM).

    Parameters:
        slope_data : dict
            Dictionary containing slope geometry, materials, and other input data
        method : str
            The limit equilibrium method name to use ('oms', 'bishop', 'janbu', 'spencer', etc.)
        rapid : bool, optional
            If True, performs rapid drawdown analysis (default: False)
        circular : bool, optional
            If True, uses circular search; if False, uses noncircular search (default: True)
        debug_level : int, optional
            Debug output level: 0=basic, 1=intermediate, 2=detailed (default: 0)
        progress_callback : callable, optional
            Called as ``progress_callback(done, total, label)`` to report progress.
            The analysis runs ``1 + 2N`` searches (one critical-surface search plus
            ``F+``/``F-`` per uncertain parameter ``N``); ``total`` is None until the
            parameter count is known. Exceptions raised by the callback are ignored.
        fs_tol, tol, max_iter : float/int, optional
            Search-convergence tolerances forwarded to the internal
            ``circular_search`` / ``noncircular_search`` calls (all ``1 + 2N`` of
            them). Any left as ``None`` use that search function's own default.
            ``tol`` only applies to circular search (noncircular has no ``tol``).
        search_opts : dict, optional
            Extra search keyword arguments applied to every one of the ``1 + 2N``
            searches — the search WINDOW above all (``entry_range`` /
            ``exit_range`` / ``center_box`` / ``tangent_depth`` /
            ``min_slip_depth``). What keeps the ``F+`` and ``F-`` searches on the
            same mechanism as the most-likely-value search they are differenced
            against.
        use_file_window : bool, optional
            Fold the model's own circles-sheet search window
            (``slope_data['search_window']``) into ``search_opts``, exactly as
            Studio's Run LEM path and a parametric sweep do — default True, so a
            windowed model gets the same surface family here as anywhere else.
            Explicit ``search_opts`` win per limit. Set False to search
            unconstrained regardless of what the file declares.

    Returns:
        tuple: (success, result) where result contains reliability analysis results

    See Also:
        reliability : the front door — ``engine='taylor'`` routes here.
        reliability_mc : the Monte Carlo counterpart (sampling on a fixed surface).
        reliability_fem : the same Taylor-series method with an FEM SSRM factor of safety.
    """

    # Only forward tolerances the user actually set; circular search also takes tol.
    _search_kwargs = {}
    if fs_tol is not None:
        _search_kwargs['fs_tol'] = fs_tol
    if max_iter is not None:
        _search_kwargs['max_iter'] = max_iter
    _circ_kwargs = dict(_search_kwargs)
    if tol is not None:
        _circ_kwargs['tol'] = tol
    if composite:
        _circ_kwargs['composite'] = True
    if seed != 'circles':
        _circ_kwargs['seed'] = seed

    # The search window, read once from the model and applied to all 1 + 2N
    # searches. A Taylor run DIFFERENCES the F+ and F- searches against each
    # other, so it matters more here than in a single run that every one of them
    # settles in the same surface family: a perturbation that jumps to a
    # competing minimum contributes its FAMILY GAP to the variance and is read as
    # sensitivity to the parameter. The circular branch takes the whole window;
    # the non-circular one takes the single limit its search understands.
    from .search import file_search_window, noncircular_search_opts
    _win = dict(search_opts or {})
    if use_file_window:
        _win.update(file_search_window(slope_data, already=_win))
    _circ_kwargs.update(_win)
    _search_kwargs.update(noncircular_search_opts(_win))

    def _progress(done, total, label):
        if progress_callback is not None:
            try:
                progress_callback(done, total, label)
            except Exception:
                pass

    start_time = time.time()

    # One full preflight of the base model before the first of 1 + 2N solves. This
    # subsumes the old "at least one standard deviation" test and adds every other
    # condition that would have surfaced later (or, for an unsupported strength
    # option, as a raise) -- see _reliability_gate.
    gate = _reliability_gate(
        slope_data,
        {'engine': 'taylor', 'rapid': rapid, 'search': search,
         'surface': 'circular' if circular else 'noncircular'},
        check_inputs=check_inputs)
    if gate:
        return False, gate

    # Import search functions and solve module here to avoid circular import
    from .search import circular_search, noncircular_search
    from . import solve

    if debug_level >= 1:
        print("=== RELIABILITY ANALYSIS ===")
        print(f"Method: {method}")
        print(f"Rapid drawdown: {rapid}")
        print(f"Circular search: {circular}")
    
    # search=False evaluates the SPECIFIED surface (circles[0] or non_circ) for
    # F_MLV and every perturbation instead of re-searching — the right mode when
    # a benchmark prescribes the slip surface (e.g. Duncan's LASH terminal,
    # corpus VP29), and immune to search pathologies on submerged slopes.
    def _solve_fixed(sd_):
        from .slice import generate_slices
        from . import solve as _solve
        # check_inputs=False: the Taylor perturbation deliberately moves c, phi and
        # gamma off their most-likely values, so a perturbed input is not a user
        # mistake and must not refuse the run.
        if circular:
            ok_, res_ = generate_slices(sd_, circle=sd_['circles'][0],
                                        num_slices=40, composite=composite,
                                        check_inputs=False)
        else:
            ok_, res_ = generate_slices(sd_, non_circ=sd_['non_circ'], num_slices=40,
                                        check_inputs=False)
        if not ok_:
            return None
        df_, surf_ = res_
        ok2, r_ = getattr(_solve, method)(df_)
        if not ok2:
            return None
        return [{"FS": r_['FS'], "slices": df_, "failure_surface": surf_,
                 "solver_result": r_}]

    # Step 1: Find the critical failure surface using search
    _progress(0, None, "Searching for the critical surface…")
    if not search:
        fs_cache = _solve_fixed(slope_data)
        converged = True
        if fs_cache is None:
            return False, "Fixed-surface evaluation failed at the most likely values"
    elif circular:
        if debug_level >= 1:
            print("Performing circular search...")
        fs_cache, converged, search_path, circle_cache = circular_search(slope_data, method, rapid=rapid, cancel_check=cancel_check, **_circ_kwargs)
    else:
        if debug_level >= 1:
            print("Performing noncircular search...")
        fs_cache, converged, search_path = noncircular_search(slope_data, method, rapid=rapid, cancel_check=cancel_check, **_search_kwargs)
    
    if not fs_cache:
        return False, "Search failed - no results found"
    
    if not converged and debug_level >= 1:
        print("Warning: Search did not fully converge - results may be less reliable")
    
    # Get the critical (minimum FS) result
    critical_result = fs_cache[0]  # First item has minimum FS
    F_MLV = critical_result["FS"]
    if F_MLV is None or not np.isfinite(F_MLV) or F_MLV >= FS_SENTINEL:
        return False, _sentinel_reason("the factor of safety at the most likely "
                                       "values (F_MLV)", float(F_MLV or np.nan))
    critical_slices = critical_result["slices"]
    critical_surface = critical_result["failure_surface"]

    if debug_level >= 1:
        print(f"Critical factor of safety (F_MLV): {F_MLV:.4f}")
    
    # Store the fs_cache for plotting
    reliability_fs_cache = [{"name": "MLV", "result": critical_result}]
    
    # Step 2: Identify parameters with standard deviations
    materials = slope_data['materials']
    
    # Find parameters that have standard deviations
    param_info = []
    
    for i, material in enumerate(materials):
        mat_name = material.get('name', f'Material_{i+1}')

        param_mappings = _strength_param_mapping(material, mat_name)

        for param, std_key in param_mappings.items():
            if std_key in material and material[std_key] > 0:
                param_info.append({
                    'material_id': i + 1,  # Use 1-based index
                    'material_name': mat_name,
                    'param': param,
                    'mlv': material[param],
                    'std': material[std_key]
                })
    
    if debug_level >= 1:
        print(f"Found {len(param_info)} parameters with standard deviations:")
        for p in param_info:
            print(f"  Material {p['material_id']}: {p['param']} = {p['mlv']:.3f} ± σ={p['std']:.3f}")
    
    # Validate up front: a strength/weight parameter cannot be reduced below zero.
    # If MLV - sigma < 0 for any parameter, the "minus sigma" perturbation is
    # non-physical, the search finds no admissible surface (returns the fs_fail
    # sentinel), and the reliability index comes out as garbage. Abort with a clear
    # message before running the expensive perturbation searches.
    invalid = [p for p in param_info if (p['mlv'] - p['std']) < 0]
    if invalid:
        details = "; ".join(
            f"material {p['material_id']} {p['param']} (mean={p['mlv']:.3g}, sigma={p['std']:.3g})"
            for p in invalid)
        return False, sigma_exceeds_mean_message(details)

    # Step 3: Calculate F+ and F- for each parameter using TSPM
    total_steps = 1 + 2 * len(param_info)   # critical search + F+/F- per parameter
    _progress(1, total_steps, "Critical surface found")
    delta_F_values = []

    for i, param in enumerate(param_info):
        from .search import _check_cancel
        _check_cancel(cancel_check)
        if debug_level >= 1:
            print(f"\nProcessing parameter {i+1}/{len(param_info)}: Material {param['material_id']}, {param['param']}")
        
        # Create modified slope_data copies
        # Shared with the FEM path: shifts the target parameter by +/- sigma AND
        # keeps gamma_sat coupled to gamma (same soil weighed two ways). The old
        # inline copy here perturbed only 'gamma', so for materials with
        # gamma_sat defined the unit-weight derivative silently evaluated to
        # ZERO wherever the slice weight came from gamma_sat (any submerged
        # mass) — found by VP29, where Duncan's gamma term is a fifth of sigma_F.
        slope_data_plus = _perturbed_slope_data(slope_data, materials, param, +1)
        slope_data_minus = _perturbed_slope_data(slope_data, materials, param, -1)
        
        # Calculate F+ and F-
        if not search:
            fs_cache_plus = _solve_fixed(slope_data_plus)
            fs_cache_minus = _solve_fixed(slope_data_minus)
        elif circular:
            fs_cache_plus, _, _, _ = circular_search(slope_data_plus, method, rapid=rapid, cancel_check=cancel_check, **_circ_kwargs)
            fs_cache_minus, _, _, _ = circular_search(slope_data_minus, method, rapid=rapid, cancel_check=cancel_check, **_circ_kwargs)
        else:
            fs_cache_plus, _, _ = noncircular_search(slope_data_plus, method, rapid=rapid, cancel_check=cancel_check, **_search_kwargs)
            fs_cache_minus, _, _ = noncircular_search(slope_data_minus, method, rapid=rapid, cancel_check=cancel_check, **_search_kwargs)
        
        if not fs_cache_plus or not fs_cache_minus:
            return False, f"Failed to calculate F+ or F- for parameter {param['param']}"
        
        F_plus = fs_cache_plus[0]["FS"]
        F_minus = fs_cache_minus[0]["FS"]

        # The post-step sentinel screen. A perturbed model with no admissible
        # surface comes back as fs_fail = 9999 and is otherwise indistinguishable
        # from a result; squared into sigma_F it swamps every other term.
        for label, value in ((f"F+ for material {param['material_id']} "
                              f"{param['param']}", F_plus),
                             (f"F- for material {param['material_id']} "
                              f"{param['param']}", F_minus)):
            if value is None or not np.isfinite(value) or value >= FS_SENTINEL:
                return False, _sentinel_reason(label, float(value or np.nan))

        # Store results for plotting
        reliability_fs_cache.append({
            "name": f"{param['param']}+",
            "result": fs_cache_plus[0]
        })
        reliability_fs_cache.append({
            "name": f"{param['param']}-",
            "result": fs_cache_minus[0]
        })
        
        delta_F = abs(F_plus - F_minus)
        delta_F_values.append(delta_F)
        
        param['F_plus'] = F_plus
        param['F_minus'] = F_minus
        param['delta_F'] = delta_F

        _progress(1 + 2 * (i + 1), total_steps,
                  f"Parameter {i + 1}/{len(param_info)}: "
                  f"mat {param['material_id']} {param['param']}")

        if debug_level >= 1:
            print(f"  F+ = {F_plus:.4f}, F- = {F_minus:.4f}, ΔF = {delta_F:.4f}")
    
    # Step 4: Calculate sigma_F and COV_F
    sigma_F = np.sqrt(sum([(df / 2)**2 for df in delta_F_values]))
    COV_F = sigma_F / F_MLV
    
    # Step 5: Calculate reliability index and probability of failure
    if COV_F < COV_FLOOR:
        return False, _cov_floor_reason(COV_F)
    
    beta_ln = np.log(F_MLV / np.sqrt(1 + COV_F**2)) / np.sqrt(np.log(1 + COV_F**2))
    reliability = norm.cdf(beta_ln)
    prob_failure = 1 - reliability
    
    if debug_level >= 1:
        print(f"\nσ_F = {sigma_F:.4f}")
        print(f"COV_F = {COV_F:.4f}")
        print(f"β_ln = {beta_ln:.4f}")
        print(f"Reliability = {reliability*100:.2f}%")
        print(f"Probability of failure = {prob_failure*100:.2f}%")
    
    # Print summary table
    if debug_level >= 0:
        print("\n=== RELIABILITY ANALYSIS RESULTS ===")
        
        # Parameter table
        table_data = []
        for param in param_info:
            table_data.append([
                f"Mat {param['material_id']} {param['param']}",
                f"{param['mlv']:.3f}",
                f"{param['std']:.3f}",
                f"{param['mlv'] + param['std']:.3f}",
                f"{param['mlv'] - param['std']:.3f}",
                f"{param['F_plus']:.3f}",
                f"{param['F_minus']:.3f}",
                f"{param['delta_F']:.3f}"
            ])
        
        headers = ["Parameter", "MLV", "σ", "MLV+σ", "MLV-σ", "F+", "F-", "ΔF"]
        colalign = ["left", "center", "center", "center", "center", "center", "center", "center"]
        print(tabulate(table_data, headers=headers, tablefmt="grid", colalign=colalign))
        
        # Summary statistics
        print(f"\nSummary Statistics:")
        print(f"F_MLV: {F_MLV:.3f}")
        print(f"σ_F: {sigma_F:.3f}")
        print(f"COV_F: {COV_F:.3f}")
        print(f"β_ln: {beta_ln:.3f}")
        print(f"Reliability: {reliability*100:.2f}%")
        print(f"Probability of failure: {prob_failure*100:.2f}%")
    
    # Prepare results
    result = {
        'method': f'{method}_reliability',
        'F_MLV': F_MLV,
        'sigma_F': sigma_F,
        'COV_F': COV_F,
        'beta_ln': beta_ln,
        'reliability': reliability,
        'prob_failure': prob_failure,
        'param_info': param_info,
        'fs_cache': reliability_fs_cache,
        'critical_surface': critical_surface,
        'critical_slices': critical_slices
    }

    elapsed = time.time() - start_time
    print(f"\nReliability analysis completed in {elapsed:.2f} seconds.")

    return True, result


def _strength_param_mapping(material, mat_name):
    """Map the strength parameters a material's model actually uses to their sigma keys.

    Perturbing a field the model does not use (e.g. phi on a cp material) would do
    nothing and silently drop that material's strength uncertainty from COV_F. gamma
    applies to every model. A blank option carries no strength parameters -- legal on
    seep-only material rows; generate_slices raises if one reaches a failure surface.

    An ``elastic`` material is the same kind of legal-but-empty row. It cannot fail and
    a slip surface cannot enter it, so it has NOTHING to perturb -- and a deterministic
    elastic zone (bedrock beneath a probabilistic model) is an ordinary, legitimate part
    of a probabilistic model. It is only an elastic material that itself carries a
    standard deviation that has to be refused, because there is no way to put that
    uncertainty through the solver.

    Raises:
        ValueError: on a strength model with no defined perturbation set, rather than
            silently falling back to Mohr-Coulomb's; and on a genuinely probabilistic
            elastic material, rather than silently ignoring its standard deviation.
    """
    option = material.get('option')
    mapping = {'gamma': 'sigma_gamma'}
    if option == 'mc':
        mapping.update({'c': 'sigma_c', 'phi': 'sigma_phi'})
    elif option == 'cp':
        mapping.update({'c': 'sigma_c', 'cp': 'sigma_cp'})
    elif option == 'elastic':
        stated = sorted(k for k in ('sigma_gamma', 'sigma_c', 'sigma_phi', 'sigma_cp')
                        if (material.get(k) or 0) > 0)
        if stated:
            raise ValueError(elastic_sigma_message(mat_name, stated))
        return {}
    elif option:
        raise ValueError(unsupported_option_message(mat_name, option))
    return mapping


def _reliability_param_info(materials):
    """Build the list of uncertain strength parameters to perturb for TSPM, shared
    by the LEM and FEM reliability paths. Perturbs only the parameters a material's
    strength model uses — see :func:`_strength_param_mapping`. Returns (param_info,
    error) — error is a message string if no sigmas are set or a mean-minus-sigma
    would go negative (non-physical), else None."""
    has_std = any(
        m.get('sigma_gamma', 0) != 0 or m.get('sigma_c', 0) != 0 or
        m.get('sigma_phi', 0) != 0 or m.get('sigma_cp', 0) != 0
        for m in materials)
    if not has_std:
        return None, no_sigmas_message()

    param_info = []
    for i, material in enumerate(materials):
        mat_name = material.get('name', f'Material_{i+1}')
        param_mappings = _strength_param_mapping(material, mat_name)
        for param, std_key in param_mappings.items():
            if std_key in material and material[std_key] > 0:
                param_info.append({
                    'material_id': i + 1, 'material_name': mat_name, 'param': param,
                    'mlv': material[param], 'std': material[std_key]})

    invalid = [p for p in param_info if (p['mlv'] - p['std']) < 0]
    if invalid:
        details = "; ".join(
            f"material {p['material_id']} {p['param']} (mean={p['mlv']:.3g}, sigma={p['std']:.3g})"
            for p in invalid)
        return None, sigma_exceeds_mean_message(details)
    return param_info, None


def _perturbed_slope_data(slope_data, materials, param, sign):
    """A shallow copy of slope_data whose materials are copied and the one target
    parameter shifted by ``sign * std`` (sign +1 -> MLV+sigma, -1 -> MLV-sigma).

    gamma and gamma_sat are the same soil weighed two ways
    (gamma_sat - gamma = n*gamma_w*(1 - S_r), correlation ~1), so perturbing
    gamma shifts gamma_sat by the SAME ABSOLUTE delta — there is deliberately no
    independent sigma_gamma_sat, which could otherwise produce moist soil
    heavier than saturated soil inside the FS derivative."""
    from .sensitivity import _set_material_field
    sd = slope_data.copy()
    sd['materials'] = [m.copy() for m in materials]
    idx = param['material_id'] - 1
    if idx < len(sd['materials']):
        # single shared mutation path with sensitivity()'s set_param — the
        # gamma/gamma_sat coupling lives there and only there
        _set_material_field(sd, idx, param['param'],
                            param['mlv'] + sign * param['std'])
    return sd


def _finalize_reliability(F_MLV, param_info, delta_F_values, method_label, debug_level=0):
    """TSPM combination shared by the LEM and FEM paths: sigma_F from the parameter
    delta_Fs, COV_F, lognormal beta, reliability and probability of failure, plus a
    printed summary table. Returns a result dict, or an error string."""
    sigma_F = np.sqrt(sum((df / 2) ** 2 for df in delta_F_values))
    COV_F = sigma_F / F_MLV if F_MLV else 0.0
    if COV_F < COV_FLOOR:
        return _cov_floor_reason(COV_F)
    beta_ln = np.log(F_MLV / np.sqrt(1 + COV_F ** 2)) / np.sqrt(np.log(1 + COV_F ** 2))
    reliability = float(norm.cdf(beta_ln))
    prob_failure = 1 - reliability

    if debug_level >= 0:
        print("\n=== RELIABILITY ANALYSIS RESULTS ===")
        table_data = [[
            f"Mat {p['material_id']} {p['param']}", f"{p['mlv']:.3f}", f"{p['std']:.3f}",
            f"{p['mlv'] + p['std']:.3f}", f"{p['mlv'] - p['std']:.3f}",
            f"{p['F_plus']:.3f}", f"{p['F_minus']:.3f}", f"{p['delta_F']:.3f}"]
            for p in param_info]
        headers = ["Parameter", "MLV", "σ", "MLV+σ", "MLV-σ", "F+", "F-", "ΔF"]
        print(tabulate(table_data, headers=headers, tablefmt="grid",
                       colalign=["left"] + ["center"] * 7))
        print(f"\nF_MLV: {F_MLV:.3f}\nσ_F: {sigma_F:.3f}\nCOV_F: {COV_F:.3f}\n"
              f"β_ln: {beta_ln:.3f}\nReliability: {reliability*100:.2f}%\n"
              f"Probability of failure: {prob_failure*100:.2f}%")

    return {
        'method': method_label, 'F_MLV': F_MLV, 'sigma_F': sigma_F, 'COV_F': COV_F,
        'beta_ln': beta_ln, 'reliability': reliability, 'prob_failure': prob_failure,
        'param_info': param_info,
    }


def reliability_fem(slope_data, mesh=None, F_min=0.5, F_max=2.0, element_type='tri6',
                    target_size=None, tolerance=0.001, failure_criterion='non_convergence',
                    max_iterations=3000, debug_level=0, progress_callback=None,
                    cancel_check=None, check_inputs=True):
    """Reliability analysis (Taylor Series Probability Method) using the FEM SSRM
    solver for the factor of safety instead of a limit-equilibrium search.

    Same method and math as ``reliability()`` — F_MLV at the most-likely values,
    then F+/F- for each uncertain strength parameter (± its standard deviation) —
    but each F comes from ``solve_ssrm``. The bracket auto-expands, so the shifted
    perturbation runs bracket robustly without hand-tuned F_min/F_max.

    ``tolerance`` (the SSRM grid/precision) defaults TIGHTER here than for a single
    solve (0.001 vs 0.01). TSPM combines ~1+2N factors of safety, and the reliability
    index is sensitive to F_MLV when COV_F is small (dβ/dF ≈ 1/(F·√(ln(1+COV²))), ≈ 9
    at COV 0.1), so an imprecise FS would visibly move β/reliability. Each SSRM here
    runs on a FIXED global grid (``grid=tolerance``; see ``solve_ssrm(grid=...)``):
    every F_MLV and perturbation lands on the same grid cell regardless of the
    F_min/F_max bracket, so the reliability is reproducible to every decimal for a
    given mesh — not just to ±tolerance/2. Results still depend on the mesh, as with
    any FE analysis.

    Perturbs the same strength parameters as the LEM path (c, phi for mc; c, cp for
    cp; gamma for both). E and nu are NOT perturbed: a sensitivity check shows E has
    no effect on the FS (halving and doubling give an identical FS) and nu only
    ~1% over its full plausible range (negligible at a realistic sigma).

    All ``1 + 2N`` trials share ONE mesh (built here if not supplied), rebuilding
    only the material→element mapping per perturbation.

    Returns (success, result) with the same reliability keys as ``reliability()``
    plus ``mlv_solution`` (the SSRM result at the most-likely values) and ``mesh``.
    """
    from .fem import build_fem_data, solve_ssrm
    from .mesh import (get_material_polygons, build_mesh_from_polygons, extract_point_constraints,
                       extract_constraint_line_geometry)
    from .search import _check_cancel

    def _progress(done, total, label):
        if progress_callback is not None:
            try:
                progress_callback(done, total, label)
            except Exception:
                pass

    start_time = time.time()
    materials = slope_data['materials']

    gate = _reliability_gate(slope_data, {'engine': 'taylor',
                                          'mesh': mesh if mesh is not None
                                          else slope_data.get('mesh')},
                             base='fem', check_inputs=check_inputs)
    if gate:
        return False, gate

    param_info, err = _reliability_param_info(materials)
    if err:
        return False, err

    # One mesh for every trial (reuse an attached mesh, else build from geometry
    # like the FEM solve path). Perturbations rebuild only build_fem_data.
    if mesh is None:
        mesh = slope_data.get('mesh')
    if mesh is None:
        if target_size is None:
            xs = [x for x, _ in slope_data['ground_surface'].coords]
            target_size = (max(xs) - min(xs)) / 100
        constraint_lines, _n_reinf, _n_pile = extract_constraint_line_geometry(slope_data)
        polygons = get_material_polygons(slope_data, reinf_lines=constraint_lines)
        mesh = build_mesh_from_polygons(polygons, target_size=target_size,
                                        element_type=element_type, lines=constraint_lines,
                                        point_constraints=extract_point_constraints(slope_data))

    # grid=tolerance: bisect each SSRM on a fixed global grid so every factor of
    # safety (F_MLV and all perturbations) is independent of the F_min/F_max bracket
    # — the reliability is then reproducible to every decimal for a given mesh, not
    # just to +/- tolerance/2. See solve_ssrm(grid=...).
    # capture_failure_state is passed PER CALL (default off, see _fs): only the
    # F_MLV solve's solution is rendered (Studio's FEM Results view via
    # mlv_solution), so only it captures the at-failure mechanism for the
    # deformation/vector/strain panels; the 2N perturbation solves use only
    # res['FS'] and skip the extra non-converging solve (it would multiply the cost
    # for nothing rendered).
    ssrm_kw = dict(tolerance=tolerance, grid=tolerance,
                   failure_criterion=failure_criterion,
                   max_iterations=max_iterations, debug_level=max(0, debug_level - 1))

    def _fs(sd, fmin, fmax, capture=False):
        res = solve_ssrm(build_fem_data(sd, mesh), F_min=fmin, F_max=fmax,
                         cancel_check=cancel_check,
                         capture_failure_state=capture, **ssrm_kw)
        if not res.get('converged'):
            return None, None, res.get('error', 'SSRM did not converge')
        return res['FS'], res, None

    if debug_level >= 1:
        print("=== RELIABILITY ANALYSIS (FEM / SSRM) ===")

    total_steps = 1 + 2 * len(param_info)
    _progress(0, total_steps, "Solving SSRM at most-likely values…")
    # capture=True here only: mlv_solution is the one SSRM result rendered (Studio's
    # FEM Results deformation/vector/strain panels), so it captures the at-failure
    # mechanism like a normal single SSRM solve. FS/reliability are unaffected.
    F_MLV, mlv_solution, err = _fs(slope_data, F_min, F_max, capture=True)
    if err:
        return False, f"Reliability (FEM): the most-likely-values solve failed — {err}"
    if debug_level >= 1:
        print(f"F_MLV = {F_MLV:.4f}")
    _progress(1, total_steps, f"F_MLV = {F_MLV:.3f}")

    # Centre the perturbation brackets on F_MLV so every bisection stays short. The
    # window must hold the LARGEST single-parameter F+/F- shift (a dominant, high-COV
    # parameter can move the FS by ~COV·F ≈ 0.3-0.5), so keep it generous: with the
    # tight tolerance the width costs only ~1 log2 step either way, but a too-narrow
    # window would trip the auto-expansion. The expansion is still a safety net for
    # anything beyond this.
    fmin_p = max(0.1, F_MLV - 0.5)
    fmax_p = F_MLV + 0.5

    delta_F_values = []
    for i, param in enumerate(param_info):
        _check_cancel(cancel_check)
        if debug_level >= 1:
            print(f"\nParameter {i+1}/{len(param_info)}: "
                  f"material {param['material_id']} {param['param']}")
        F_plus, _sp, e1 = _fs(_perturbed_slope_data(slope_data, materials, param, +1), fmin_p, fmax_p)
        F_minus, _sm, e2 = _fs(_perturbed_slope_data(slope_data, materials, param, -1), fmin_p, fmax_p)
        if e1 or e2:
            return False, (f"Reliability (FEM): perturbation solve failed for material "
                           f"{param['material_id']} {param['param']} — {e1 or e2}")
        delta_F = abs(F_plus - F_minus)
        delta_F_values.append(delta_F)
        param.update(F_plus=F_plus, F_minus=F_minus, delta_F=delta_F)
        if debug_level >= 1:
            print(f"  F+ = {F_plus:.4f}, F- = {F_minus:.4f}, ΔF = {delta_F:.4f}")
        _progress(1 + 2 * (i + 1), total_steps,
                  f"Parameter {i+1}/{len(param_info)}: mat {param['material_id']} {param['param']}")

    result = _finalize_reliability(F_MLV, param_info, delta_F_values,
                                   method_label='fem_reliability', debug_level=debug_level)
    if isinstance(result, str):
        return False, result
    result['mlv_solution'] = mlv_solution
    result['mesh'] = mesh

    print(f"\nReliability (FEM) analysis completed in {time.time() - start_time:.2f} seconds.")
    return True, result


# Monte Carlo reliability defaults. MC_DEFAULT_SEED is a FIXED CONSTANT (never
# time-based): a given input file must reproduce the same beta and probability of
# failure on every run so the regression suite can lock them. Override per-call
# with reliability_mc(rng_seed=..., n_samples=...).
MC_DEFAULT_SEED = 20240117
MC_DEFAULT_SAMPLES = 10000

# Physical floor each sampled parameter is truncated to (a strength or unit weight
# cannot go negative). This truncation IS the phi>=0 bound that governs the
# high-COV problems (e.g. VP34's Phase I fill, COV 124%): a plain normal would
# draw negative friction angles, and clamping them at zero is exactly how the
# published Monte Carlo treatments handle that bound.
_MC_FLOOR = {'gamma': 1e-6, 'c': 0.0, 'phi': 0.0, 'cp': 0.0}


def _mc_param_info(materials):
    """Uncertain strength/weight parameters for Monte Carlo, as a list of dicts
    (material_id, param, mlv, std). Same parameter set as the Taylor-series path
    (:func:`_reliability_param_info`) but WITHOUT its mean-minus-sigma>=0 rejection:
    Monte Carlo handles COV>100% by truncating samples at the physical floor, so a
    high-COV parameter that TSPM must decline is admissible here. Returns
    (param_info, error); error is a message when no sigmas are set."""
    has_std = any(
        m.get('sigma_gamma', 0) != 0 or m.get('sigma_c', 0) != 0 or
        m.get('sigma_phi', 0) != 0 or m.get('sigma_cp', 0) != 0
        for m in materials)
    if not has_std:
        return None, no_sigmas_message()
    param_info = []
    for i, material in enumerate(materials):
        mat_name = material.get('name', f'Material_{i+1}')
        for param, std_key in _strength_param_mapping(material, mat_name).items():
            if std_key in material and material[std_key] > 0:
                param_info.append({'material_id': i + 1, 'material_name': mat_name,
                                   'param': param, 'mlv': material[param],
                                   'std': material[std_key]})
    return param_info, None


def _mc_sampled_slope_data(slope_data, materials, param_info, values):
    """A shallow copy of slope_data whose materials are copied and every uncertain
    parameter set to its sampled value for this trial. Uses the same
    ``_set_material_field`` mutation path as the TSPM perturbation, so the
    gamma/gamma_sat coupling is applied identically."""
    from .sensitivity import _set_material_field
    sd = slope_data.copy()
    sd['materials'] = [m.copy() for m in materials]
    for p, v in zip(param_info, values):
        idx = p['material_id'] - 1
        if idx < len(sd['materials']):
            _set_material_field(sd, idx, p['param'], float(v))
    return sd


def reliability_mc(slope_data, method, rapid=False, circular=True, debug_level=0,
                   n_samples=MC_DEFAULT_SAMPLES, rng_seed=MC_DEFAULT_SEED,
                   distribution='normal', search=True, num_slices=40,
                   progress_callback=None, cancel_check=None,
                   fs_tol=None, tol=None, max_iter=None, composite=False,
                   seed='circles', search_opts=None, use_file_window=True,
                   check_inputs=True):
    """Monte Carlo reliability analysis — the sampling counterpart to the
    Taylor-series :func:`reliability`.

    Draws ``n_samples`` independent realizations of every uncertain material
    parameter from the standard deviations in the mat sheet (s(g), s(c), s(f),
    s(c/p)), evaluates the factor of safety of each realization on a FIXED failure
    surface, and reports the sample statistics.

    **Limit-equilibrium only.** This is deliberately an LEM path: a Monte Carlo
    campaign needs 10^4 factor-of-safety evaluations, which is affordable with a
    limit-equilibrium solve but not with the finite-element SSRM. FEM reliability
    stays on the Taylor series (:func:`reliability_fem`, 1+2N solves). ``method``
    must therefore name an LEM solver ('bishop', 'spencer', ...).

    **The failure surface is never randomized.** The slip surface is a decision
    variable, not a random variable, so its geometry is held fixed across all
    realizations (the given surface with ``search=False``, or the most-likely-values
    critical surface with ``search=True``). This fixed-surface campaign is the
    analogue of Slide2's "Global Minimum" probabilistic method and matches what
    every published benchmark did (a prescribed or deterministic-critical surface).
    See the Reliability documentation for the surface-treatment discussion and the
    optional per-realization-minimum mode.

    Parameters
    ----------
    method : str
        Limit-equilibrium method name ('bishop', 'spencer', ...).
    n_samples : int
        Number of realizations (default 10000).
    rng_seed : int
        Seed for ``numpy.random.default_rng`` — a fixed constant by default, so the
        result is bit-reproducible. Pass a different int to see sampling scatter.
    distribution : {'normal', 'lognormal'}
        Per-parameter input distribution. Default 'normal' (mean = MLV, std =
        sigma), the same interpretation the Taylor-series method places on the
        sigma columns. 'lognormal' matches the same mean and std (method of
        moments) and requires a positive mean. Every parameter is truncated at its
        physical floor (>= 0), which is how the friction-angle >= 0 bound is handled
        on high-COV problems.
    search : bool
        With ``search=False`` the specified surface (``circles[0]`` or ``non_circ``)
        is evaluated for every realization — the right mode when a benchmark
        prescribes the slip surface (VP28, VP34). With ``search=True`` the critical
        surface is found once at the most-likely values and then held fixed across
        all realizations (a full re-search per realization is not performed —
        prohibitive at 10^4 samples). Noncircular ``search=True`` is not supported.
    num_slices : int
        Slice count for each evaluation (default 40, matching ``reliability``'s
        fixed-surface path so the MC mean FS and the TSPM F_MLV are comparable).
    search_opts : dict, optional
        Extra search keyword arguments for the ``search=True`` critical-surface
        search — the search WINDOW above all (``entry_range`` / ``exit_range`` /
        ``center_box`` / ``tangent_depth`` / ``min_slip_depth``). That one search
        fixes the surface every realization is evaluated on, so it decides which
        mechanism the whole campaign describes.
    use_file_window : bool, optional
        Fold the model's own circles-sheet search window
        (``slope_data['search_window']``) into ``search_opts``, exactly as
        Studio's Run LEM path and a parametric sweep do — default True. Explicit
        ``search_opts`` win per limit. Set False to search unconstrained
        regardless of what the file declares. Ignored when ``search=False``,
        which searches for nothing.

    Returns
    -------
    (success, result) : tuple
        On success ``result`` carries ``mean_FS``, ``sigma_F``, ``COV_F``, both
        reliability-index conventions (``beta_normal`` = (mean-1)/sigma_F,
        ``beta_ln`` = lognormal from the sample moments), the empirical probability
        of failure ``pf_empirical`` (fraction of realizations with FS<1), the
        distribution-fitted ``pf_normal`` / ``pf_lognormal``, the raw ``fs_samples``,
        and the sampled input matrix ``param_samples`` (n_samples x n_params, column
        order matching ``param_info``) so a global-sensitivity measure such as a
        Spearman rank correlation of each input against FS can reuse the same draws.

    See Also
    --------
    reliability : the front door — ``engine='mc'`` routes here.
    reliability_taylor : the Taylor-series counterpart (1 + 2N searches).
    """
    from .search import circular_search, noncircular_search, _check_cancel
    from .slice import generate_slices
    from . import solve

    def _progress(done, total, label):
        if progress_callback is not None:
            try:
                progress_callback(done, total, label)
            except Exception:
                pass

    start_time = time.time()
    materials = slope_data['materials']

    gate = _reliability_gate(
        slope_data,
        {'engine': 'mc', 'rapid': rapid, 'search': search, 'n_samples': n_samples,
         'surface': 'circular' if circular else 'noncircular'},
        check_inputs=check_inputs)
    if gate:
        return False, gate

    param_info, err = _mc_param_info(materials)
    if err:
        return False, err

    _LEM_METHODS = {'oms', 'bishop', 'janbu', 'corps', 'lowe', 'spencer', 'mprice'}
    if method not in _LEM_METHODS:
        return False, (f"Monte Carlo reliability is limit-equilibrium only; '{method}' "
                       f"is not an LEM solver. Use one of {sorted(_LEM_METHODS)}. FEM "
                       "reliability uses the Taylor series (reliability_fem).")
    solver = getattr(solve, method)

    # Only forward tolerances the caller actually set.
    _search_kwargs = {}
    if fs_tol is not None:
        _search_kwargs['fs_tol'] = fs_tol
    if max_iter is not None:
        _search_kwargs['max_iter'] = max_iter
    _circ_kwargs = dict(_search_kwargs)
    if tol is not None:
        _circ_kwargs['tol'] = tol
    if composite:
        _circ_kwargs['composite'] = True
    if seed != 'circles':
        _circ_kwargs['seed'] = seed

    # The model's own search window, read the same way every other searching path
    # reads it. The one search this engine runs fixes the surface all n_samples
    # realizations are evaluated on, so the window decides which mechanism the
    # reported distribution belongs to. Only the circular branch searches here at
    # all -- search=True is refused for noncircular below -- so there is nothing
    # to hand the non-circular subset to.
    _win = dict(search_opts or {})
    if use_file_window:
        from .search import file_search_window
        _win.update(file_search_window(slope_data, already=_win))
    _circ_kwargs.update(_win)

    # ---- Resolve the fixed evaluation surface -----------------------------
    fixed_circle = None
    fixed_noncirc = None
    if not search:
        if circular:
            circs = slope_data.get('circles')
            if not circs:
                return False, "Monte Carlo (search=False): no circle is specified in the input."
            fixed_circle = circs[0]
        else:
            fixed_noncirc = slope_data.get('non_circ')
            if not fixed_noncirc:
                return False, "Monte Carlo (search=False): no non-circular surface is specified."
    else:
        if not circular:
            return False, ("Monte Carlo with search=True is only supported for circular "
                           "surfaces; specify the surface and use search=False for "
                           "noncircular problems.")
        if debug_level >= 1:
            print("Finding the critical circle at the most-likely values…")
        fs_cache, _conv, _path, ccache = circular_search(
            slope_data, method, rapid=rapid, cancel_check=cancel_check, **_circ_kwargs)
        if not fs_cache:
            return False, "Monte Carlo: critical-surface search failed."
        crit = (ccache[0] if ccache else fs_cache[0])
        fixed_circle = {'Xo': crit['Xo'], 'Yo': crit['Yo'], 'R': crit['R'],
                        'Depth': crit.get('Depth')}

    def _eval(sd):
        # check_inputs=False: every Monte Carlo realization is a sampled model, so
        # the same reasoning as the Taylor path applies.
        if circular:
            ok, res = generate_slices(sd, circle=fixed_circle, num_slices=num_slices,
                                      composite=composite, check_inputs=False)
        else:
            ok, res = generate_slices(sd, non_circ=fixed_noncirc,
                                      num_slices=num_slices, check_inputs=False)
        if not ok:
            return None
        ok2, r = solver(res[0])
        if not ok2:
            return None
        fs = r.get('FS')
        if fs is None or not np.isfinite(fs):
            return None
        return float(fs)

    F_MLV = _eval(slope_data)
    if F_MLV is None:
        return False, "Monte Carlo: evaluation at the most-likely values failed."

    if debug_level >= 1:
        print("=== MONTE CARLO RELIABILITY ANALYSIS ===")
        print(f"Method: {method} | samples: {n_samples} | seed: {rng_seed} | "
              f"distribution: {distribution}")
        print(f"F at most-likely values: {F_MLV:.4f}")

    # ---- Draw the sample matrix (n_samples x n_params) --------------------
    rng = np.random.default_rng(rng_seed)
    npar = len(param_info)
    sample_matrix = np.empty((n_samples, npar))
    for j, p in enumerate(param_info):
        mlv, std = p['mlv'], p['std']
        if distribution == 'lognormal':
            if mlv <= 0:
                return False, (f"lognormal distribution requires a positive mean "
                               f"(material {p['material_id']} {p['param']} mean={mlv}).")
            s_ln = np.sqrt(np.log(1.0 + (std / mlv) ** 2))
            m_ln = np.log(mlv) - 0.5 * s_ln ** 2
            col = rng.lognormal(m_ln, s_ln, n_samples)
        elif distribution == 'normal':
            col = rng.normal(mlv, std, n_samples)
        else:
            return False, f"Unknown distribution '{distribution}'. Use 'normal' or 'lognormal'."
        col = np.maximum(col, _MC_FLOOR.get(p['param'], 0.0))
        if p['param'] == 'phi':
            col = np.minimum(col, 89.0)
        sample_matrix[:, j] = col

    # ---- Evaluate every realization on the fixed surface ------------------
    fs_vals = np.empty(n_samples)
    valid = np.ones(n_samples, dtype=bool)
    report_every = max(1, n_samples // 20)
    for k in range(n_samples):
        if k % report_every == 0:
            _check_cancel(cancel_check)
            _progress(k, n_samples, f"Monte Carlo sample {k}/{n_samples}")
        sd_k = _mc_sampled_slope_data(slope_data, materials, param_info, sample_matrix[k])
        fk = _eval(sd_k)
        if fk is None:
            valid[k] = False
            fs_vals[k] = np.nan
        else:
            fs_vals[k] = fk
    _progress(n_samples, n_samples, "Monte Carlo complete")

    fs_ok = fs_vals[valid]
    n_valid = int(fs_ok.size)
    n_invalid = n_samples - n_valid
    if n_valid < 2:
        return False, ("Monte Carlo: fewer than two admissible realizations — the "
                       "sampled surface is non-analyzable for almost every draw.")

    mean_FS = float(np.mean(fs_ok))
    sigma_F = float(np.std(fs_ok, ddof=1))
    COV_F = sigma_F / mean_FS if mean_FS else 0.0

    # Empirical probability of failure over the admissible realizations. An
    # inadmissible draw (no analyzable surface) is not counted as FS<1 — it is
    # reported separately so a large invalid fraction is never hidden inside PF.
    n_fail = int(np.count_nonzero(fs_ok < 1.0))
    pf_empirical = n_fail / n_valid

    beta_normal = (mean_FS - 1.0) / sigma_F if sigma_F > 0 else float('inf')
    pf_normal = float(norm.cdf(-beta_normal))
    if COV_F > 0:
        beta_ln = float(np.log(mean_FS / np.sqrt(1 + COV_F ** 2)) /
                        np.sqrt(np.log(1 + COV_F ** 2)))
    else:
        beta_ln = float('inf')
    pf_lognormal = float(norm.cdf(-beta_ln))

    if debug_level >= 0:
        print("\n=== MONTE CARLO RELIABILITY RESULTS ===")
        table = [[f"Mat {p['material_id']} {p['param']}", f"{p['mlv']:.3f}",
                  f"{p['std']:.3f}", f"{p['std'] / p['mlv'] * 100:.1f}%" if p['mlv'] else "—"]
                 for p in param_info]
        print(tabulate(table, headers=["Parameter", "MLV", "σ", "COV"],
                       tablefmt="grid", colalign=["left", "center", "center", "center"]))
        print(f"\nSamples: {n_samples} (valid {n_valid}, invalid {n_invalid}) | "
              f"seed {rng_seed} | {distribution}")
        print(f"Mean FS: {mean_FS:.4f}   σ_F: {sigma_F:.4f}   COV_F: {COV_F:.4f}")
        print(f"β (normal): {beta_normal:.4f}   β (lognormal): {beta_ln:.4f}")
        print(f"PF empirical: {pf_empirical*100:.3f}%   "
              f"PF normal: {pf_normal*100:.3f}%   PF lognormal: {pf_lognormal*100:.3f}%")

    result = {
        'method': f'{method}_reliability_mc',
        'F_MLV': F_MLV,
        'mean_FS': mean_FS,
        'sigma_F': sigma_F,
        'COV_F': COV_F,
        'beta_ln': beta_ln,
        'beta_normal': beta_normal,
        'reliability': 1.0 - pf_empirical,
        'prob_failure': pf_empirical,
        'pf_empirical': pf_empirical,
        'pf_normal': pf_normal,
        'pf_lognormal': pf_lognormal,
        'n_samples': n_samples,
        'n_valid': n_valid,
        'n_invalid': n_invalid,
        'rng_seed': rng_seed,
        'distribution': distribution,
        'param_info': param_info,
        'fs_samples': fs_vals,
        # The sampled input matrix (n_samples x n_params), column j is param_info[j].
        # Carried so a global-sensitivity measure (e.g. Spearman rank correlation of
        # each input against FS — see sensitivity.mc_rank_correlation) can be computed
        # from the same realizations, without re-sampling.
        'param_samples': sample_matrix,
    }

    print(f"\nMonte Carlo reliability analysis completed in "
          f"{time.time() - start_time:.2f} seconds.")
    return True, result
