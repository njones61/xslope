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

* :func:`reliability` -- the front door; ``engine='taylor'`` (default), ``'mc'`` or ``'rs'``.
* :func:`reliability_taylor` -- Taylor Series Probability Method (limit equilibrium).
* :func:`reliability_mc` -- Monte Carlo reliability on a fixed surface.
* :func:`reliability_rs` -- response-surface reliability on a fixed surface: a
  quadratic surrogate fitted to a few dozen real solves, gated against held-out
  real solves, then sampled millions of times.
* :func:`reliability_fem` -- Taylor series driven by the finite-element SSRM.

Historically this code lived in :mod:`xslope.advanced`, which now re-exports every
name here for backward compatibility, so ``from xslope.advanced import reliability``
(and ``reliability_mc`` / ``reliability_fem`` / the private helpers) keeps working.
"""

import time

import math
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
        'lowe', 'spencer', 'mprice'). Every engine takes it as its first analysis
        argument. (Passed as the second positional/keyword, exactly as before.)
    engine : {'taylor', 'mc', 'rs'}, keyword-only
        Which reliability engine to run:

        * ``'taylor'`` (default) -> :func:`reliability_taylor`, the Taylor Series
          Probability Method (1 + 2N limit-equilibrium searches).
        * ``'mc'`` -> :func:`reliability_mc`, a Monte Carlo sampling campaign on a
          fixed failure surface.
        * ``'rs'`` -> :func:`reliability_rs`, a response-surface campaign: a
          quadratic surrogate fitted to a few dozen real solves, sampled millions
          of times, and checked against held-out real solves before it is used.

        Aliases are accepted case-insensitively ('tspm' for Taylor;
        'monte_carlo' / 'montecarlo' for MC; 'response_surface' / 'rsm' for RS).
    *args, **kwargs
        Forwarded verbatim to the selected engine. See each engine's signature for
        the accepted options (they overlap but are not identical — e.g. Monte Carlo
        adds ``n_samples`` / ``rng_seed`` / ``distribution``, and the response
        surface adds ``n_surrogate`` / ``n_gate`` / ``alpha``).

    Returns
    -------
    (success, result) : tuple
        Whatever the selected engine returns. A response-surface run whose
        surrogate fails its accuracy gate returns ``(False, message)`` like any
        other refusal — the caller decides whether to fall back to Monte Carlo.

    See Also
    --------
    reliability_taylor : Taylor Series Probability Method (limit equilibrium).
    reliability_mc : Monte Carlo reliability on a fixed surface (limit equilibrium).
    reliability_rs : Response-surface reliability on a fixed surface (limit equilibrium).
    reliability_fem : Taylor series with the finite-element SSRM factor of safety.
    """
    key = str(engine).lower().replace('-', '_')
    if key in ('taylor', 'tspm', 'taylor_series'):
        return reliability_taylor(slope_data, method, *args, **kwargs)
    if key in ('mc', 'monte_carlo', 'montecarlo'):
        return reliability_mc(slope_data, method, *args, **kwargs)
    if key in ('rs', 'rsm', 'response_surface', 'responsesurface'):
        return reliability_rs(slope_data, method, *args, **kwargs)
    raise ValueError(
        f"reliability(): unknown engine={engine!r}. Use 'taylor' (default, the "
        f"Taylor Series Probability Method), 'mc' (Monte Carlo) or 'rs' "
        f"(response surface).")


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

    # Center the perturbation brackets on F_MLV so every bisection stays short. The
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


def _draw_sample_matrix(param_info, n, rng, distribution, sampling='random'):
    """Draw ``n`` realizations of every uncertain parameter — the ONE place the
    input distributions are defined.

    Returns ``(matrix, error)``: an ``n x len(param_info)`` array whose column j is
    ``param_info[j]``, or ``(None, message)`` for an unusable distribution.

    Both sampling engines draw here, which is the point: a response surface that
    was sampled from anything other than the distribution the real path draws from
    would answer a different question, however well the surrogate fitted. That
    includes the truncations — every draw is floored at :data:`_MC_FLOOR` (a
    strength or unit weight cannot go negative) and a friction angle is additionally
    capped at 89 degrees, so the samples the surrogate sees are the samples the
    solver would have been handed.

    Drawing is column-by-column and stateless apart from ``rng``, so a long campaign
    may draw in chunks (``_draw_sample_matrix(..., chunk, rng, ...)`` repeatedly)
    without changing the distribution — only the ordering of the stream.

    ``sampling='lhs'`` draws by Latin hypercube instead: each marginal is
    stratified into ``n`` equal-probability bins with one draw per bin
    (scipy's ``qmc.LatinHypercube``, seeded from ``rng`` — deterministic like
    everything else), pushed through the same inverse CDFs and the same
    truncations. Stratification reduces the sampling variance of the
    estimates; chunked draws stratify per chunk, which is valid but weaker
    than one full-campaign hypercube. NOTE: LHS draws are not independent, so
    the convergence stop's i.i.d. confidence half-width OVERSTATES the
    uncertainty under LHS — the stop errs conservative (runs longer than
    strictly needed), never the other way.
    """
    npar = len(param_info)
    if sampling not in ('random', 'lhs'):
        return None, f"Unknown sampling '{sampling}'. Use 'random' or 'lhs'."
    u = None
    if sampling == 'lhs':
        from scipy.stats import qmc
        u = qmc.LatinHypercube(d=npar, seed=rng).random(n)
    out = np.empty((n, npar))
    for j, p in enumerate(param_info):
        mlv, std = p['mlv'], p['std']
        if distribution == 'lognormal':
            if mlv <= 0:
                return None, (f"lognormal distribution requires a positive mean "
                              f"(material {p['material_id']} {p['param']} mean={mlv}).")
            s_ln = np.sqrt(np.log(1.0 + (std / mlv) ** 2))
            m_ln = np.log(mlv) - 0.5 * s_ln ** 2
            if u is None:
                col = rng.lognormal(m_ln, s_ln, n)
            else:
                col = np.exp(m_ln + s_ln * norm.ppf(u[:, j]))
        elif distribution == 'normal':
            if u is None:
                col = rng.normal(mlv, std, n)
            else:
                col = norm.ppf(u[:, j], loc=mlv, scale=std)
        else:
            return None, f"Unknown distribution '{distribution}'. Use 'normal' or 'lognormal'."
        out[:, j] = _truncate_samples(col, p['param'])
    return out, None


def _truncate_samples(col, param):
    """Apply the physical floor (and the friction-angle cap) to a column of draws.

    Split out from :func:`_draw_sample_matrix` because the response surface has to
    truncate its DESIGN points the same way: a design point at MLV - 2 sigma that
    lands below the floor is solved at the floor, and the surrogate is fitted
    against the coordinate that was actually solved.
    """
    col = np.maximum(col, _MC_FLOOR.get(param, 0.0))
    if param == 'phi':
        col = np.minimum(col, 89.0)
    return col


def _fixed_surface_context(slope_data, method, label, rapid=False, circular=True,
                           search=True, num_slices=40, composite=False, seed='circles',
                           fs_tol=None, tol=None, max_iter=None, search_opts=None,
                           use_file_window=True, cancel_check=None, debug_level=0):
    """Everything a fixed-surface sampling campaign needs before it can sample.

    Resolves the one failure surface every realization is evaluated on, builds the
    slice template that makes a realization cost milliseconds, evaluates the factor
    of safety at the most-likely values, and works out what the solver can be seeded
    with. Returns ``(True, ctx)`` or ``(False, message)``; ``ctx`` carries

    ``evaluate``
        ``sd -> float | None``, the factor of safety of a sampled slope-data dict on
        the fixed surface (None when the model has no admissible solution).
    ``F_MLV``
        The factor of safety at the most-likely values, evaluated WITHOUT a solver
        seed so it is the same number the engine has always reported.
    ``circle`` / ``non_circ``
        The fixed surface, for the record.
    ``template``
        The slice template, or None when the model declined the shortcut.

    ``label`` names the engine in the refusal messages ("Monte Carlo", "Response
    surface"), which is the only thing that differs between the two callers.
    """
    from .search import circular_search
    from .slice import generate_slices, update_slice_materials
    from . import solve

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
    # reads it. The one search this engine runs fixes the surface all realizations
    # are evaluated on, so the window decides which mechanism the reported
    # distribution belongs to. Only the circular branch searches here at all --
    # search=True is refused for noncircular below -- so there is nothing to hand
    # the non-circular subset to.
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
                return False, f"{label} (search=False): no circle is specified in the input."
            fixed_circle = circs[0]
        else:
            fixed_noncirc = slope_data.get('non_circ')
            if not fixed_noncirc:
                return False, f"{label} (search=False): no non-circular surface is specified."
    else:
        if not circular:
            return False, (f"{label} with search=True is only supported for circular "
                           f"surfaces; specify the surface and use search=False for "
                           f"noncircular problems.")
        if debug_level >= 1:
            print("Finding the critical circle at the most-likely values…")
        fs_cache, _conv, _path, ccache = circular_search(
            slope_data, method, rapid=rapid, cancel_check=cancel_check, **_circ_kwargs)
        if not fs_cache:
            return False, f"{label}: critical-surface search failed."
        # The circle every realization is evaluated on is the CRITICAL one. It is
        # taken from circle_cache rather than fs_cache because only circle_cache
        # carries the radius — and by lowest FS rather than by position, because
        # circle_cache is in the order the circles were TRIED: its first entry is
        # the search's first trial circle, which is not an answer to anything.
        crit = min(ccache, key=lambda c: c['FS']) if ccache else fs_cache[0]
        fixed_circle = {'Xo': crit['Xo'], 'Yo': crit['Yo'], 'R': crit.get('R'),
                        'Depth': crit.get('Depth')}

    # The slice-template fast path: the fixed surface means the geometry of the
    # slice table never changes across realizations — only the material
    # properties multiplied onto it. Generate once with the per-band breakdown
    # and rewrite the material-derived columns per realization; realizations
    # fall back to a full rebuild whenever the model declines the shortcut
    # (stress-dependent envelopes, ru pore pressure, suction), and the template
    # is verified against a rebuilt table at the most-likely values before it
    # is trusted at all.
    _template = None
    if not rapid:
        try:
            if circular:
                okt, rest = generate_slices(slope_data, circle=fixed_circle,
                                            num_slices=num_slices,
                                            composite=composite,
                                            check_inputs=False,
                                            material_breakdown=True)
            else:
                okt, rest = generate_slices(slope_data, non_circ=fixed_noncirc,
                                            num_slices=num_slices,
                                            check_inputs=False,
                                            material_breakdown=True)
            if okt:
                cand = rest[0]
                probe = cand.copy()
                update_slice_materials(probe, slope_data['materials'])
                same = all(
                    np.allclose(probe[col].astype(float),
                                cand[col].astype(float),
                                rtol=0, atol=1e-9, equal_nan=True)
                    for col in ('w', 'kw', 'y_cg', 'c', 'phi', 'c1', 'phi1',
                                'd', 'psi'))
                if same:
                    _template = cand
        except (ValueError, KeyError):
            _template = None

    _seed_kw = {}

    def _eval(sd):
        # check_inputs=False: every realization is a sampled model, so the same
        # reasoning as the Taylor path applies.
        if _template is not None:
            df = _template.copy()
            try:
                update_slice_materials(df, sd['materials'])
            except ValueError:
                return _eval_rebuild(sd)
            ok2, r = solver(df, **_seed_kw)
            if not ok2:
                return None
            fs = r.get('FS')
            if fs is None or not np.isfinite(fs):
                return None
            return float(fs)
        return _eval_rebuild(sd)

    def _eval_rebuild(sd):
        if circular:
            ok, res = generate_slices(sd, circle=fixed_circle, num_slices=num_slices,
                                      composite=composite, check_inputs=False)
        else:
            ok, res = generate_slices(sd, non_circ=fixed_noncirc,
                                      num_slices=num_slices, check_inputs=False)
        if not ok:
            return None
        ok2, r = solver(res[0], **_seed_kw)
        if not ok2:
            return None
        fs = r.get('FS')
        if fs is None or not np.isfinite(fs):
            return None
        return float(fs)

    F_MLV = _eval(slope_data)
    if F_MLV is None:
        return False, f"{label}: evaluation at the most-likely values failed."

    # Seed every realization's iterative solve from the most-likely-values
    # solution: each sampled model's root sits within a few percent of it, so
    # Newton (Spencer) and the fixed-point loop (Bishop) start almost converged.
    # Methods whose signatures take no seed run exactly as before.
    import inspect as _inspect
    _mlv_theta = None
    _solver_params = set(_inspect.signature(solver).parameters)
    if 'theta_seed' in _solver_params and _template is not None:
        # One extra solve at the most-likely values fetches Spencer's theta for
        # the seed; the F_MLV above stays the seedless, historical evaluation.
        try:
            _df0 = _template.copy()
            update_slice_materials(_df0, slope_data['materials'])
            _ok0, _r0 = solver(_df0)
            if _ok0 and isinstance(_r0, dict):
                _mlv_theta = _r0.get('theta')
        except ValueError:
            pass
    if 'fs_seed' in _solver_params:
        _seed_kw['fs_seed'] = F_MLV
    if 'theta_seed' in _solver_params and _mlv_theta is not None:
        _seed_kw['theta_seed'] = _mlv_theta

    return True, {'evaluate': _eval, 'F_MLV': F_MLV, 'circle': fixed_circle,
                  'non_circ': fixed_noncirc, 'template': _template}


#: The limit-equilibrium solvers a sampling engine will run. Both Monte Carlo and
#: the response surface need thousands to millions of factor-of-safety evaluations,
#: which is a limit-equilibrium proposition only.
_LEM_METHODS = {'oms', 'bishop', 'janbu', 'corps', 'lowe', 'spencer', 'mprice'}


def reliability_mc(slope_data, method, rapid=False, circular=True, debug_level=0,
                   n_samples=MC_DEFAULT_SAMPLES, rng_seed=MC_DEFAULT_SEED,
                   distribution='normal', search=True, num_slices=40,
                   progress_callback=None, cancel_check=None,
                   fs_tol=None, tol=None, max_iter=None, composite=False,
                   seed='circles', search_opts=None, use_file_window=True,
                   check_inputs=True, converge_rel=None, converge_check=100,
                   converge_min=500, sampling='lhs'):
    """Monte Carlo reliability analysis — the sampling counterpart to the
    Taylor-series :func:`reliability`.

    Draws ``n_samples`` realizations of every uncertain material parameter
    (Latin hypercube by default; ``sampling='random'`` for independent draws)
    from the standard deviations in the mat sheet (s(g), s(c), s(f),
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
    converge_rel : float, optional
        Statistical-convergence stopping tolerance on the empirical probability
        of failure, RELATIVE to its own value (0.10 = stop when P_f is known to
        ±10% of itself at 95% confidence). A relative rule self-scales: an
        absolute half-percentage-point on a P_f of 17% and of 2% are entirely
        different amounts of knowledge, and the sample count the rule demands
        grows as (1−p)/p, roughly ten times more realizations at 2% than at
        17%. When set, the campaign checks every ``converge_check``
        realizations and stops once 1.96·sqrt(p(1−p)/n) ≤ converge_rel·p;
        ``n_samples`` becomes the cap. The rule never fires before
        ``converge_min`` realizations or before 10 failures have been observed
        (below that the normal approximation behind the half-width is not
        valid, which protects rare-event problems from stopping on false
        confidence). Because the whole sample matrix is drawn up front from the
        fixed seed, a converged run is exactly the first ``n_used`` realizations
        of the full run — bit-reproducible, like everything else here. Default
        None: run all ``n_samples``, the previous behavior unchanged.
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
    from .search import _check_cancel

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

    if method not in _LEM_METHODS:
        return False, (f"Monte Carlo reliability is limit-equilibrium only; '{method}' "
                       f"is not an LEM solver. Use one of {sorted(_LEM_METHODS)}. FEM "
                       "reliability uses the Taylor series (reliability_fem).")

    okc, ctx = _fixed_surface_context(
        slope_data, method, "Monte Carlo", rapid=rapid, circular=circular,
        search=search, num_slices=num_slices, composite=composite, seed=seed,
        fs_tol=fs_tol, tol=tol, max_iter=max_iter, search_opts=search_opts,
        use_file_window=use_file_window, cancel_check=cancel_check,
        debug_level=debug_level)
    if not okc:
        return False, ctx
    _eval = ctx['evaluate']
    F_MLV = ctx['F_MLV']

    if debug_level >= 1:
        print("=== MONTE CARLO RELIABILITY ANALYSIS ===")
        print(f"Method: {method} | samples: {n_samples} | seed: {rng_seed} | "
              f"distribution: {distribution}")
        print(f"F at most-likely values: {F_MLV:.4f}")

    # ---- Draw the sample matrix (n_samples x n_params) --------------------
    rng = np.random.default_rng(rng_seed)
    sample_matrix, err = _draw_sample_matrix(param_info, n_samples, rng,
                                             distribution, sampling=sampling)
    if err:
        return False, err

    # ---- Evaluate every realization on the fixed surface ------------------
    fs_vals = np.empty(n_samples)
    valid = np.ones(n_samples, dtype=bool)
    report_every = max(1, n_samples // 20)
    n_used = n_samples
    pf_trace = [] if converge_rel is not None else None
    pf_ci_half = None
    run_fail = run_valid = 0
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
            run_valid += 1
            if fk < 1.0:
                run_fail += 1
        if (converge_rel is not None and (k + 1) % converge_check == 0
                and run_valid > 0):
            p = run_fail / run_valid
            half = 1.96 * math.sqrt(p * (1.0 - p) / run_valid)
            pf_trace.append((run_valid, p, half))
            # The half-width is a normal approximation; below ~10 observed
            # failures it is not credible, so a rare-event campaign keeps
            # sampling to the cap rather than stopping on false confidence.
            if (k + 1 >= converge_min and run_fail >= 10 and p > 0
                    and half <= converge_rel * p):
                n_used = k + 1
                pf_ci_half = half
                break
    _progress(n_used, n_used, "Monte Carlo complete")
    if n_used < n_samples:
        fs_vals = fs_vals[:n_used]
        valid = valid[:n_used]

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
        print(f"\nSamples: {n_used} of {n_samples} requested "
              f"(valid {n_valid}, invalid {n_invalid}) | seed {rng_seed} | {distribution}")
        if converge_rel is not None:
            if n_used < n_samples:
                print(f"Converged at n = {n_used}: P_f {pf_empirical*100:.2f}% "
                      f"known to ±{pf_ci_half/pf_empirical*100:.1f}% of itself "
                      f"(target ±{converge_rel*100:.0f}%, 95% CI)")
            else:
                print(f"Convergence target ±{converge_rel*100:.0f}% of P_f not "
                      f"met within the {n_samples}-sample cap.")
        print(f"Mean FS: {mean_FS:.4f}   σ_F: {sigma_F:.4f}   COV_F: {COV_F:.4f}")
        print(f"β (normal): {beta_normal:.4f}   β (lognormal): {beta_ln:.4f}")
        print(f"PF empirical: {pf_empirical*100:.3f}%   "
              f"PF normal: {pf_normal*100:.3f}%   PF lognormal: {pf_lognormal*100:.3f}%")

    result = {
        'n_requested': n_samples,
        'n_used': n_used,
        'pf_ci_half': pf_ci_half,
        'pf_trace': pf_trace,
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


# ---------------------------------------------------------------------------
# Response surface
# ---------------------------------------------------------------------------

#: Realizations the surrogate is sampled at. Ten million is not a cost decision —
#: evaluating a polynomial that many times is a few seconds of vector arithmetic —
#: it is a resolution decision: at this count the sampling half-width on a
#: probability of failure of 17% is +/-0.02% absolute, which is far below the
#: surrogate's own fit error, so the answer's uncertainty is the FIT and nothing
#: else. Reaching the same sampling resolution with real solves would take ~53
#: million of them.
RS_DEFAULT_SURROGATE = 10_000_000

#: Realizations evaluated per vectorized chunk. The chunk is drawn, standardized,
#: expanded into the quadratic basis and reduced to running totals, so peak memory
#: is chunk x (number of polynomial terms) doubles rather than n_surrogate of them.
RS_DEFAULT_CHUNK = 1_000_000

#: Held-out REAL solves the fitted surrogate is measured against before any of its
#: answers are believed. Drawn from the same distributions as the surrogate sampling
#: and never used in the fit.
RS_GATE_SAMPLES = 500

#: Axial (star) distance of the central composite design, in standard deviations.
#: The corners sit at +/-1 sigma; the axial points reach further out so the
#: curvature the quadratic has to carry is fitted over the range the campaign
#: actually samples rather than extrapolated into it.
RS_DEFAULT_ALPHA = 2.0

#: The surrogate is accepted when its root-mean-square error over the held-out
#: gate draws is within this fraction of the real spread of the factor of safety
#: over those same draws, AND its coefficient of determination is at least
#: :data:`RS_R2_MIN`, AND at most :data:`RS_PF_DISAGREE_MAX` of the gate draws are
#: put on the wrong side of F = 1. All three are measured, never assumed.
#:
#: The two fit thresholds are calibrated against what the fit error costs the
#: ANSWER, measured by solving tens of thousands of realizations for real and
#: predicting the same draws with the surrogate: at an axial distance of 2 sigma
#: the submerged-slope model (two parameters, 50,000 paired realizations) fits to
#: 0.042 of the real spread with R^2 = 0.9982 and lands its probability of failure
#: within 0.5% of its own value, and Duncan's LASH terminal (VP29, 20,000 paired
#: realizations) fits to 0.016 with R^2 = 0.9997 for the same 0.5%. A 10,000-sample
#: Monte Carlo of either carries a 95% sampling half-width of about +/-4% of P_f, so
#: a surrogate admitted at these thresholds is well inside the noise of the campaign
#: it replaces. A tighter R^2 bar (0.999) would refuse the first of those two
#: surrogates while it was answering to 0.5%.
RS_RMS_FRACTION = 0.05
RS_R2_MIN = 0.995

#: Fraction of the gate draws the surrogate may put on the wrong side of F = 1.
#: This is the fit error measured in the units the answer is reported in: each
#: disagreement is a realization the surrogate counts as a failure and the solver
#: does not, or the reverse, so the disagreement rate BOUNDS the error in the
#: probability of failure (the two directions partly cancel in the count itself —
#: measured, the net error is under half the disagreement rate). The two models
#: above disagree on 0.2-0.5% of their draws.
RS_PF_DISAGREE_MAX = 0.02

#: A gate draw the real pipeline cannot solve has no surrogate counterpart — the
#: polynomial answers everywhere. A few such draws are excluded and reported; more
#: than this fraction means a material part of the sampled population is
#: inadmissible, which a surrogate defined over the whole space cannot represent,
#: and the run is refused rather than quietly answering for the admissible part.
RS_MAX_GATE_INVALID = 0.02

#: Realizations from the FAILURE REGION the surrogate is additionally checked on:
#: draws the surrogate itself counted as F < 1, taken from the sampling pass and
#: solved for real. A gate drawn uniformly from the input distributions barely
#: visits that region — it is a few percent of the population — and the region is
#: the whole of the answer, so it is measured separately. On VP34 a uniform gate
#: found 1.4% of its draws inadmissible while 36% of the realizations the surrogate
#: was counting as failures had no real solution at all: inadmissibility is
#: concentrated exactly where P_f is counted.
RS_TAIL_SAMPLES = 200

#: The failure-region checks. Above either fraction the run is refused: the
#: surrogate is either counting realizations the model cannot solve, or classifying
#: them differently from the solver where the classification IS the answer.
RS_TAIL_INVALID_MAX = 0.05
RS_TAIL_DISAGREE_MAX = 0.10

#: Fewest failure-region draws that make those rates worth testing. Below this the
#: numbers are reported and the refusal is not applied — a model with almost no
#: predicted failures has almost no counted mass to get wrong.
RS_TAIL_MIN = 20

#: Surrogate realizations retained for plotting. The histogram is drawn from a
#: fixed-stride subsample so the result dict stays a few hundred kilobytes instead
#: of hundreds of megabytes; every reported statistic still comes from the full
#: n_surrogate draws.
RS_PLOT_SAMPLES = 10_000

#: Upper bound on the design size. The factorial part of a central composite design
#: doubles with every uncertain parameter, so a model with many sigmas would spend
#: its afternoon solving corners. Refused rather than run.
RS_MAX_DESIGN = 4096


def _ccd_design(d, alpha):
    """Coded coordinates of the central composite design in ``d`` parameters.

    Rows are, in order: the ``2^d`` factorial corners at +/-1, the ``2d`` axial
    points at +/-``alpha``, and the centre. Coded units are standard deviations
    about the most-likely values.
    """
    import itertools
    corners = np.array(list(itertools.product((-1.0, 1.0), repeat=d)), dtype=float)
    axial = np.zeros((2 * d, d))
    for j in range(d):
        axial[2 * j, j] = -alpha
        axial[2 * j + 1, j] = alpha
    return np.vstack([corners, axial, np.zeros((1, d))])


def _quad_terms(z):
    """Full quadratic basis of the standardized coordinates ``z`` (n x d).

    Columns are the constant, the ``d`` linear terms, then the ``d(d+1)/2``
    quadratic terms (squares and cross products) in row-major order — the same
    order the fitted coefficient vector is stored in.
    """
    n, d = z.shape
    cols = [np.ones(n)]
    for j in range(d):
        cols.append(z[:, j])
    for j in range(d):
        for k in range(j, d):
            cols.append(z[:, j] * z[:, k])
    return np.column_stack(cols)


def _quad_term_labels(param_info):
    """Human-readable names for the columns of :func:`_quad_terms`."""
    names = [f"m{p['material_id']}.{p['param']}" for p in param_info]
    labels = ['1'] + list(names)
    for j, a in enumerate(names):
        for b in names[j:]:
            labels.append(f"{a}*{b}" if a != b else f"{a}^2")
    return labels


def reliability_rs(slope_data, method, rapid=False, circular=True, debug_level=0,
                   n_surrogate=RS_DEFAULT_SURROGATE, n_gate=RS_GATE_SAMPLES,
                   alpha=RS_DEFAULT_ALPHA, chunk=RS_DEFAULT_CHUNK,
                   rng_seed=MC_DEFAULT_SEED, distribution='normal', search=True,
                   num_slices=40, progress_callback=None, cancel_check=None,
                   fs_tol=None, tol=None, max_iter=None, composite=False,
                   seed='circles', search_opts=None, use_file_window=True,
                   check_inputs=True, sampling='lhs'):
    """Response-surface reliability analysis — a Monte Carlo campaign whose factor
    of safety comes from a fitted surrogate instead of a solve, with the surrogate
    measured against real solves before any of its answers are used.

    The engine samples the same distributions as :func:`reliability_mc`, on the same
    fixed failure surface, and reports the same statistics. What changes is where
    the factor of safety comes from:

    1. **Design.** A central composite design about the most-likely values — the
       ``2^d`` factorial corners at +/-1 sigma, ``2d`` axial points at
       +/-``alpha`` sigma, and the centre — each solved with the real pipeline.
    2. **Fit.** A full quadratic in the ``d`` uncertain parameters
       (``1 + d + d(d+1)/2`` coefficients), least squares.
    3. **Gate.** ``n_gate`` further realizations, drawn from the sampling
       distributions and never used in the fit, solved for real and compared with
       the surrogate. The run is REFUSED unless the root-mean-square error is
       within :data:`RS_RMS_FRACTION` of the real spread over those draws, the
       coefficient of determination is at least :data:`RS_R2_MIN`, and no more
       than :data:`RS_PF_DISAGREE_MAX` of the draws are put on the wrong side of
       F = 1.
    4. **Sampling.** ``n_surrogate`` realizations (default ten million) evaluated
       through the polynomial in vectorized chunks.

    The arithmetic that motivates it: the 95% half-width on an empirical
    probability of failure is ``1.96*sqrt(p(1-p)/n)``, so resolving a P_f near 17%
    to +/-0.01% absolute needs ~53 million real solves — days of them. The
    surrogate reaches that resolution for a few dozen real solves plus the gate,
    and moves the error budget from sampling noise, which the report can no longer
    see, to fit error, which the gate measures and prints.

    Parameters
    ----------
    n_surrogate : int
        Realizations drawn from the surrogate (default ten million).
    n_gate : int
        Held-out real solves the surrogate is measured against (default 500).
    alpha : float
        Axial distance of the design, in standard deviations (default 2.0).
    chunk : int
        Surrogate realizations evaluated per vectorized pass. Reduced
        automatically for a model with many parameters, whose quadratic basis is
        wider.
    rng_seed : int
        Seed for the surrogate draws. The gate draws come from an independent
        stream of the same seed, so both are reproducible and the gate is never
        fitted.

    Everything else is interpreted exactly as in :func:`reliability_mc`
    (``distribution``, ``search``, ``num_slices``, ``search_opts``,
    ``use_file_window``, the solver tolerances).

    Returns
    -------
    (success, result) : tuple
        On success ``result`` carries the Monte Carlo result keys — ``mean_FS``,
        ``sigma_F``, ``COV_F``, ``beta_normal``, ``beta_ln``, ``pf_empirical``,
        ``pf_normal``, ``pf_lognormal``, ``F_MLV``, ``param_info`` — so a
        histogram or a rank correlation reads it unchanged, plus the surrogate's
        credentials: ``n_design`` and ``n_real_solves``, ``gate_rms``,
        ``gate_max_error``, ``gate_r2``, ``gate_sigma``, ``gate_pf_disagree``,
        ``rs_coefficients`` and ``rs_terms``. ``fs_samples`` and ``param_samples`` are a fixed-stride
        SUBSAMPLE of the surrogate draws (``fs_samples_stride`` realizations
        apart, :data:`RS_PLOT_SAMPLES` of them) so plotting stays light; every
        statistic reported comes from all ``n_surrogate``.

        A surrogate that fails the gate returns ``(False, message)`` naming the
        measured errors.

    See Also
    --------
    reliability : the front door — ``engine='rs'`` routes here.
    reliability_mc : the same campaign with every realization solved for real.
    """
    from .search import _check_cancel

    def _progress(done, total, label):
        if progress_callback is not None:
            try:
                progress_callback(done, total, label)
            except Exception:
                pass

    start_time = time.time()
    materials = slope_data['materials']

    gate_msg = _reliability_gate(
        slope_data,
        {'engine': 'rs', 'rapid': rapid, 'search': search,
         'surface': 'circular' if circular else 'noncircular'},
        check_inputs=check_inputs)
    if gate_msg:
        return False, gate_msg

    param_info, err = _mc_param_info(materials)
    if err:
        return False, err

    if method not in _LEM_METHODS:
        return False, (f"Response-surface reliability is limit-equilibrium only; "
                       f"'{method}' is not an LEM solver. Use one of "
                       f"{sorted(_LEM_METHODS)}. FEM reliability uses the Taylor "
                       f"series (reliability_fem).")

    d = len(param_info)
    n_design = 2 ** d + 2 * d + 1
    if n_design > RS_MAX_DESIGN:
        return False, (
            f"Response surface: this model carries {d} uncertain parameters, so the "
            f"central composite design would be {n_design} real solves "
            f"(2^{d} corners + {2 * d} axial + 1 centre), past the {RS_MAX_DESIGN} "
            f"limit. Monte Carlo samples the same distributions without a design "
            f"that doubles with every parameter.")
    n_terms = 1 + d + d * (d + 1) // 2

    okc, ctx = _fixed_surface_context(
        slope_data, method, "Response surface", rapid=rapid, circular=circular,
        search=search, num_slices=num_slices, composite=composite, seed=seed,
        fs_tol=fs_tol, tol=tol, max_iter=max_iter, search_opts=search_opts,
        use_file_window=use_file_window, cancel_check=cancel_check,
        debug_level=debug_level)
    if not okc:
        return False, ctx
    _eval = ctx['evaluate']
    F_MLV = ctx['F_MLV']

    mlv = np.array([p['mlv'] for p in param_info], dtype=float)
    std = np.array([p['std'] for p in param_info], dtype=float)

    if debug_level >= 1:
        print("=== RESPONSE-SURFACE RELIABILITY ANALYSIS ===")
        print(f"Method: {method} | parameters: {d} | design: {n_design} solves | "
              f"surrogate draws: {n_surrogate} | seed: {rng_seed} | {distribution}")
        print(f"F at most-likely values: {F_MLV:.4f}")

    # ---- 1. Design: solve the central composite design for real -----------
    coded = _ccd_design(d, float(alpha))
    # Physical coordinates, truncated exactly as a sampled realization is: a design
    # point below a parameter's physical floor is SOLVED at the floor, and the fit
    # is then carried out against the coordinate that was actually solved rather
    # than the nominal one that was not.
    phys = mlv + coded * std
    for j, p in enumerate(param_info):
        phys[:, j] = _truncate_samples(phys[:, j], p['param'])

    # Peak memory during sampling is one chunk of the quadratic basis; a model with
    # many parameters has a wider basis, so the chunk narrows to hold the product
    # roughly constant. Settled here so the progress total counts real passes.
    chunk = int(max(50_000, min(int(chunk), 12_000_000 // max(1, n_terms))))
    n_chunks = max(1, -(-int(n_surrogate) // chunk))
    total_steps = n_design + int(n_gate) + n_chunks + RS_TAIL_SAMPLES
    fs_design = np.empty(n_design)
    for i in range(n_design):
        if i % 8 == 0:
            _check_cancel(cancel_check)
            _progress(i, total_steps, f"Design solve {i + 1}/{n_design}")
        sd_i = _mc_sampled_slope_data(slope_data, materials, param_info, phys[i])
        fi = _eval(sd_i)
        if fi is None:
            where = ", ".join(
                f"material {p['material_id']} {p['param']}={phys[i, j]:.4g}"
                for j, p in enumerate(param_info))
            return False, (
                f"Response surface: the design point {i + 1} of {n_design} has no "
                f"analyzable solution ({where}). The surrogate is fitted to the "
                f"whole design, so a design point that cannot be solved leaves the "
                f"fit undefined over part of the sampled range. Monte Carlo can "
                f"report on such a model — it excludes the inadmissible draws and "
                f"counts them.")
        fs_design[i] = fi

    # ---- 2. Fit the quadratic ---------------------------------------------
    z_design = (phys - mlv) / std
    X = _quad_terms(z_design)
    coef, _res, rank, _sv = np.linalg.lstsq(X, fs_design, rcond=None)

    def _surrogate(x_phys):
        """Factor of safety of physical parameter rows through the fitted quadratic."""
        return _quad_terms((x_phys - mlv) / std) @ coef

    # ---- 3. The gate: held-out real solves ---------------------------------
    # An independent stream of the same seed: reproducible, and never drawn from
    # the stream the surrogate is sampled with, so no gate point can also be a
    # sampled point by construction.
    gate_rng = np.random.default_rng([int(rng_seed), 1])
    gate_x, err = _draw_sample_matrix(param_info, int(n_gate), gate_rng, distribution)
    if err:
        return False, err
    gate_real = np.empty(int(n_gate))
    gate_ok = np.ones(int(n_gate), dtype=bool)
    for i in range(int(n_gate)):
        if i % 25 == 0:
            _check_cancel(cancel_check)
            _progress(n_design + i, total_steps, f"Gate solve {i + 1}/{n_gate}")
        sd_i = _mc_sampled_slope_data(slope_data, materials, param_info, gate_x[i])
        fi = _eval(sd_i)
        if fi is None:
            gate_ok[i] = False
            gate_real[i] = np.nan
        else:
            gate_real[i] = fi
    n_gate_valid = int(np.count_nonzero(gate_ok))
    n_gate_invalid = int(n_gate) - n_gate_valid
    if n_gate_valid < max(20, n_terms + 1):
        return False, (
            f"Response surface: only {n_gate_valid} of {n_gate} gate realizations "
            f"could be solved, which is too few to measure the surrogate against. "
            f"A surrogate that cannot be checked is not used.")
    if n_gate_invalid > RS_MAX_GATE_INVALID * int(n_gate):
        return False, (
            f"Response surface: {n_gate_invalid} of {n_gate} gate realizations "
            f"({n_gate_invalid / n_gate * 100:.1f}%) have no analyzable solution, "
            f"past the {RS_MAX_GATE_INVALID * 100:.0f}% limit. The surrogate returns "
            f"a factor of safety everywhere, so it would answer for a population "
            f"that part of the real model cannot. Use Monte Carlo, which excludes "
            f"and counts the inadmissible draws.")

    gate_pred = _surrogate(gate_x[gate_ok])
    gate_truth = gate_real[gate_ok]
    resid = gate_pred - gate_truth
    gate_rms = float(np.sqrt(np.mean(resid ** 2)))
    gate_max = float(np.max(np.abs(resid)))
    gate_sigma = float(np.std(gate_truth, ddof=1))
    ss_tot = float(np.sum((gate_truth - np.mean(gate_truth)) ** 2))
    gate_r2 = float(1.0 - np.sum(resid ** 2) / ss_tot) if ss_tot > 0 else float('nan')
    # The fit error in the units the answer is reported in: a gate draw the
    # surrogate and the solver disagree about is a realization that would be
    # counted into the probability of failure by one and not the other.
    n_disagree = int(np.count_nonzero((gate_pred < 1.0) != (gate_truth < 1.0)))
    pf_disagree = n_disagree / n_gate_valid

    rms_limit = RS_RMS_FRACTION * gate_sigma
    if not (gate_rms <= rms_limit and gate_r2 >= RS_R2_MIN
            and pf_disagree <= RS_PF_DISAGREE_MAX):
        return False, (
            f"Response surface: the fitted surrogate did not pass its accuracy gate "
            f"on {n_gate_valid} held-out real solves. RMS error {gate_rms:.5f} "
            f"(limit {rms_limit:.5f} = {RS_RMS_FRACTION:.0%} of the real spread "
            f"sigma = {gate_sigma:.4f}), maximum error {gate_max:.5f}, R^2 = "
            f"{gate_r2:.6f} (minimum {RS_R2_MIN}), and {n_disagree} of "
            f"{n_gate_valid} draws ({pf_disagree * 100:.1f}%, limit "
            f"{RS_PF_DISAGREE_MAX * 100:.0f}%) put on the wrong side of F = 1. A "
            f"quadratic in {d} parameter(s) does not describe this model's factor "
            f"of safety well enough for its tail to be trusted; run Monte Carlo, "
            f"which evaluates every realization.")

    # ---- 4. Sample the surrogate ------------------------------------------
    stride = max(1, int(n_surrogate) // RS_PLOT_SAMPLES)
    rng = np.random.default_rng(rng_seed)
    n_left = int(n_surrogate)
    offset = 0
    n_fail = 0
    sum_d = 0.0
    sum_d2 = 0.0
    keep_fs = []
    keep_x = []
    # The first realizations the surrogate itself counts as failures, kept for the
    # failure-region gate below. They are the head of an i.i.d. stream, so they are
    # a fair sample of the region the probability of failure is a count of.
    tail_x = []
    n_tail_want = int(RS_TAIL_SAMPLES)
    step = n_design + int(n_gate)
    while n_left > 0:
        _check_cancel(cancel_check)
        m = min(chunk, n_left)
        xs, err = _draw_sample_matrix(param_info, m, rng, distribution,
                                      sampling=sampling)
        if err:
            return False, err
        fs = _surrogate(xs)
        failed = fs < 1.0
        n_fail += int(np.count_nonzero(failed))
        if len(tail_x) < n_tail_want and failed.any():
            tail_x.extend(xs[failed][:n_tail_want - len(tail_x)])
        dev = fs - F_MLV
        sum_d += float(dev.sum())
        sum_d2 += float((dev ** 2).sum())
        i0 = (-offset) % stride
        if i0 < m:
            keep_fs.append(fs[i0::stride])
            keep_x.append(xs[i0::stride])
        offset += m
        n_left -= m
        step += 1
        _progress(step, total_steps,
                  f"Surrogate realizations {offset:,}/{n_surrogate:,}")

    n_tot = int(n_surrogate)
    mean_FS = F_MLV + sum_d / n_tot
    var = (sum_d2 - n_tot * (mean_FS - F_MLV) ** 2) / (n_tot - 1)
    sigma_F = float(np.sqrt(max(var, 0.0)))
    COV_F = sigma_F / mean_FS if mean_FS else 0.0
    pf_empirical = n_fail / n_tot

    beta_normal = (mean_FS - 1.0) / sigma_F if sigma_F > 0 else float('inf')
    pf_normal = float(norm.cdf(-beta_normal))
    if COV_F > 0:
        beta_ln = float(np.log(mean_FS / np.sqrt(1 + COV_F ** 2)) /
                        np.sqrt(np.log(1 + COV_F ** 2)))
    else:
        beta_ln = float('inf')
    pf_lognormal = float(norm.cdf(-beta_ln))

    fs_samples = np.concatenate(keep_fs) if keep_fs else np.empty(0)
    param_samples = (np.vstack(keep_x) if keep_x
                     else np.empty((0, d)))

    # ---- 5. The failure-region gate ---------------------------------------
    # Everything above is a statement about a region a uniform gate hardly
    # visits. These are the realizations the surrogate COUNTED, solved for real.
    tail_n = len(tail_x)
    tail_invalid = tail_disagree = None
    tail_rms = None
    if tail_n:
        tail_arr = np.asarray(tail_x)
        tail_real = np.empty(tail_n)
        tail_ok = np.ones(tail_n, dtype=bool)
        for i in range(tail_n):
            if i % 25 == 0:
                _check_cancel(cancel_check)
                _progress(n_design + int(n_gate) + n_chunks + i, total_steps,
                          f"Failure-region check {i + 1}/{tail_n}")
            sd_i = _mc_sampled_slope_data(slope_data, materials, param_info, tail_arr[i])
            fi = _eval(sd_i)
            if fi is None:
                tail_ok[i] = False
                tail_real[i] = np.nan
            else:
                tail_real[i] = fi
        n_tail_invalid = int(tail_n - np.count_nonzero(tail_ok))
        tail_invalid = n_tail_invalid / tail_n
        n_tail_valid = int(np.count_nonzero(tail_ok))
        if n_tail_valid:
            tail_pred = _surrogate(tail_arr[tail_ok])
            tail_disagree = float(np.count_nonzero(tail_real[tail_ok] >= 1.0)
                                  / n_tail_valid)
            tail_rms = float(np.sqrt(np.mean((tail_pred - tail_real[tail_ok]) ** 2)))
        if tail_n >= RS_TAIL_MIN:
            if tail_invalid > RS_TAIL_INVALID_MAX:
                return False, (
                    f"Response surface: {n_tail_invalid} of {tail_n} realizations the "
                    f"surrogate counts as failures ({tail_invalid * 100:.0f}%, limit "
                    f"{RS_TAIL_INVALID_MAX * 100:.0f}%) have no analyzable solution at "
                    f"all. The probability of failure would be a count of realizations "
                    f"the model cannot solve, and the surrogate cannot tell them from "
                    f"the ones it can. Use Monte Carlo, which excludes and counts the "
                    f"inadmissible draws — note that its probability of failure is then "
                    f"a fraction of the ADMISSIBLE realizations.")
            if tail_disagree is not None and tail_disagree > RS_TAIL_DISAGREE_MAX:
                return False, (
                    f"Response surface: the solver puts {tail_disagree * 100:.0f}% of "
                    f"the realizations the surrogate counts as failures back above "
                    f"F = 1 (limit {RS_TAIL_DISAGREE_MAX * 100:.0f}%), with an RMS "
                    f"error of {tail_rms:.5f} there against {gate_rms:.5f} over the "
                    f"population as a whole. The surrogate is least accurate exactly "
                    f"where the answer is counted; run Monte Carlo.")

    if debug_level >= 0:
        print("\n=== RESPONSE-SURFACE RELIABILITY RESULTS ===")
        table = [[f"Mat {p['material_id']} {p['param']}", f"{p['mlv']:.3f}",
                  f"{p['std']:.3f}", f"{p['std'] / p['mlv'] * 100:.1f}%" if p['mlv'] else "—"]
                 for p in param_info]
        print(tabulate(table, headers=["Parameter", "MLV", "σ", "COV"],
                       tablefmt="grid", colalign=["left", "center", "center", "center"]))
        print(f"\nDesign: {n_design} real solves (2^{d} corners at ±1σ, {2 * d} axial "
              f"at ±{alpha:g}σ, 1 centre) | quadratic terms: {n_terms} | rank {rank}")
        print(f"Gate: {n_gate_valid} held-out real solves | R² = {gate_r2:.6f} | "
              f"RMS = {gate_rms:.5f} (≤ {rms_limit:.5f}) | max = {gate_max:.5f} | "
              f"real σ over the gate draws = {gate_sigma:.4f}")
        print(f"Gate F<1 disagreement: {n_disagree} of {n_gate_valid} draws "
              f"({pf_disagree * 100:.2f}%, limit {RS_PF_DISAGREE_MAX * 100:.0f}%) — "
              f"the measured bound on the error in P_f")
        if tail_n:
            print(f"Failure region: {tail_n} counted failures solved for real | "
                  f"inadmissible {tail_invalid * 100:.1f}% "
                  f"(limit {RS_TAIL_INVALID_MAX * 100:.0f}%) | solver disputes "
                  f"{tail_disagree * 100:.1f}% of them "
                  f"(limit {RS_TAIL_DISAGREE_MAX * 100:.0f}%) | RMS there "
                  f"{tail_rms:.5f}")
        else:
            print("Failure region: the surrogate predicts no failures, so there is "
                  "nothing to check there.")
        print(f"Surrogate: {n_tot:,} realizations | seed {rng_seed} | {distribution} "
              f"| real solves in total: {n_design + int(n_gate) + tail_n + 1}")
        print(f"Mean FS: {mean_FS:.4f}   σ_F: {sigma_F:.4f}   COV_F: {COV_F:.4f}")
        print(f"β (normal): {beta_normal:.4f}   β (lognormal): {beta_ln:.4f}")
        print(f"PF empirical: {pf_empirical * 100:.3f}%   "
              f"PF normal: {pf_normal * 100:.3f}%   PF lognormal: {pf_lognormal * 100:.3f}%")

    result = {
        'method': f'{method}_reliability_rs',
        'engine': 'rs',
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
        'n_surrogate': n_tot,
        'n_samples': n_tot,
        'n_used': n_tot,
        'n_valid': n_tot,
        'n_invalid': 0,
        'rng_seed': rng_seed,
        'distribution': distribution,
        'param_info': param_info,
        # The design and its verdict — the credentials of every number above.
        'n_design': n_design,
        'n_gate': int(n_gate),
        'n_gate_valid': n_gate_valid,
        'n_gate_invalid': n_gate_invalid,
        'n_real_solves': n_design + int(n_gate) + tail_n + 1,
        'alpha': float(alpha),
        'gate_rms': gate_rms,
        'gate_max_error': gate_max,
        'gate_r2': gate_r2,
        'gate_sigma': gate_sigma,
        'gate_rms_limit': rms_limit,
        'gate_pf_disagree': pf_disagree,
        'gate_pf_disagree_n': n_disagree,
        # The failure-region check: how many of the realizations the surrogate
        # counted as failures the solver cannot solve, and how many it puts back
        # above F = 1.
        'tail_n': tail_n,
        'tail_invalid': tail_invalid,
        'tail_disagree': tail_disagree,
        'tail_rms': tail_rms,
        'rs_coefficients': coef,
        'rs_terms': _quad_term_labels(param_info),
        # A fixed-stride subsample of the surrogate draws, for plotting only: the
        # histogram reads fs_samples and a rank correlation reads both, while every
        # statistic above comes from all n_surrogate realizations.
        'fs_samples': fs_samples,
        'param_samples': param_samples,
        'fs_samples_stride': stride,
    }

    print(f"\nResponse-surface reliability analysis completed in "
          f"{time.time() - start_time:.2f} seconds.")
    return True, result
