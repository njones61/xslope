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

import time

import numpy as np
import pandas as pd
from scipy.stats import norm
from tabulate import tabulate


def validate_rapid_drawdown(slope_data):
    """
    Validates that slope_data has the required inputs for rapid drawdown analysis.
    Prints warnings for missing optional inputs and raises ValueError for missing required inputs.
    """
    materials = slope_data['materials']

    # Hard error: at least one material must have non-zero d or psi
    has_d_psi = any(m.get('d', 0) != 0 or m.get('psi', 0) != 0 for m in materials)
    if not has_d_psi:
        raise ValueError("Rapid drawdown requires at least one material with non-zero d or psi values. Check your input template.")

    # Warning: no second set of distributed loads
    if not slope_data.get('dloads2'):
        print("[WARNING] Rapid drawdown: no second set of distributed loads (dloads2) found.")

    # Warning: piezo method selected but no second piezo line
    has_piezo = any(m.get('u') == 'piezo' for m in materials)
    if has_piezo and not slope_data.get('piezo_line2'):
        print("[WARNING] Rapid drawdown: piezo method selected but no second piezo line found.")

    # Warning: seep method selected but no second seep solution
    has_seep = any(m.get('u') == 'seep' for m in materials)
    if has_seep and 'seep_u2' not in slope_data:
        print("[WARNING] Rapid drawdown: seep method selected but no second seep solution found.")


def rapid_drawdown(df, method_name, debug_level=1):
    """
    Performs rapid drawdown analysis using a three-stage approach.
    
    Parameters:
        df: pandas.DataFrame
            Slice data with all required columns including rapid drawdown specific data:
            - c, phi: current strength parameters
            - c1, phi1: original strength parameters (for stage 3)
            - d, psi: rapid drawdown parameters for low-K materials
            - u: pore pressure (stage 1)
            - u2: pore pressure for lowered pool (stage 2)
            - dload, d_x, d_y: distributed loads (stage 1)
            - dload2, d_x2, d_y2: distributed loads for lowered pool (stage 2)
        method_name: str
            The method name to use ('oms', 'bishop', 'spencer', etc.)
        debug_level: int
            0: no output, 1: print FS at each stage, >1: detailed debug info
    
    Returns:
        Tuple(bool, dict): (True, result_dict) or (False, error_message)
    """
    # Curved-envelope materials are incompatible with the staged strength
    # overrides of the drawdown procedure (both mutate per-slice c/phi); refuse
    # clearly rather than silently letting one clobber the other.
    for _flag, _opt, _name in (('pow_flag', 'pow', 'power-curve'),
                               ('hb_flag', 'hb', 'Hoek-Brown')):
        if _flag in df.columns and bool(df[_flag].any()):
            return False, (f"One or more slices use the {_name} strength option "
                           f"(option='{_opt}'), which is not supported in rapid "
                           "drawdown analysis.")

    # Import solve module and get the method function
    from . import solve
    method_func = getattr(solve, method_name)

    # Work on a copy: the analysis overwrites strength and load columns (Stage 2
    # swaps in the drawdown pore pressures / loads and undrained strengths), so do
    # not mutate the caller's slice DataFrame. We DO write the winning stage's
    # per-slice results (n_eff, u, c, phi) back to it at the end, so the caller —
    # and the plots/search cache that read it — see the rapid-drawdown stresses
    # rather than the stale Stage-0 values.
    caller_df = df
    df = df.copy()

    # Validate that d and psi parameters are present for at least some slices
    if (df['d'] == 0).all() and (df['psi'] == 0).all():
        return False, "Rapid drawdown requires d and psi parameters for low-K materials. All values are zero — check your input template."

    if debug_level >= 1:
        print("=== RAPID DRAWDOWN ANALYSIS ===")
    
    # Stage 1: Pre-drawdown conditions
    if debug_level >= 1:
        print("Stage 1: Pre-drawdown conditions...")
    
    # Use original conditions (c, phi, u, dload, d_x, d_y)
    success, result_stage1 = method_func(df)
    if not success:
        return False, f"Stage 1 failed: {result_stage1}"
    
    stage1_FS = result_stage1['FS']
    if debug_level >= 1:
        print(f"Stage 1 FS = {stage1_FS:.4f}")

    # The three-stage procedure presumes the slope is stable BEFORE drawdown: stage 1
    # exists only to supply consolidation stresses. If FS < 1 the mobilized shear
    # tau_fc = (1/FS)(c' + sigma'_fc tan phi') exceeds the failure envelope, so the
    # consolidation stress state lies above it and K1 > Kf. Since K1 rises
    # monotonically with tau_fc (hence with 1/FS) and equals Kf exactly at FS = 1,
    # `stage1_FS < 1` is precisely the condition `K1 > Kf`, uniformly over slices.
    #
    # Eq (5) then EXTRAPOLATES beyond the Kc=Kf envelope -- which the source defines as
    # the physical extreme -- driving tau_ff below it and eventually negative, where the
    # max(0, tau_ff) clamp below turns it into a silent zero-strength slice. A search
    # calling this on trial surfaces would let such a surface win on a fictitious FS~0.
    # The negative-stress fallback does NOT catch this: sigma'_3c stays positive.
    if stage1_FS < 1.0:
        return False, (
            f"Rapid drawdown requires a slope that is stable before drawdown, but the "
            f"Stage 1 (full pool) factor of safety is {stage1_FS:.4f} < 1. The "
            f"consolidation stresses are undefined because the stress state lies above "
            f"the failure envelope (K1 > Kf)."
        )

    # Calculate consolidation stresses for each slice
    # N_eff should be available from the method function
    if 'n_eff' not in df.columns:
        return False, "Stage 1 did not compute n_eff values"
    
    # Calculate sigma_fc and tau_fc for each slice
    sigma_fc = df['n_eff'] / df['dl']  # Equation (2)
    tau_fc = (1.0 / stage1_FS) * (df['c'] + sigma_fc * np.tan(np.radians(df['phi'])))  # Equation (3)
    
    if debug_level >= 2:
        print("Stage 1 consolidation stresses:")
        for i in range(len(df)):
            print(f"  Slice {i+1}: sigma_fc = {sigma_fc.iloc[i]:.2f}, tau_fc = {tau_fc.iloc[i]:.2f}")
    
    # Stage 2: Post-drawdown conditions with undrained strengths
    if debug_level >= 1:
        print("Stage 2: Post-drawdown conditions with undrained strengths...")
   
    # Update pore pressures and distributed loads for stage 2
    df['u'] = df['u2']
    df['dload'] = df['dload2']
    df['d_x'] = df['d_x2']
    df['d_y'] = df['d_y2']
    
    # Process each slice for undrained strength calculation
    for i in range(len(df)):
        # Check if this slice has low-K material (d or psi are not zero)
        d_val = df.iloc[i]['d']
        psi_val = df.iloc[i]['psi']

        if d_val > 0 or psi_val > 0:
            # Low-K material - calculate undrained strength
            if debug_level >= 2:
                print(f"Processing low-K material for slice {i+1}")
            
            # Get consolidation stresses for this slice
            sigma_fc_i = sigma_fc.iloc[i]
            tau_fc_i = tau_fc.iloc[i]
            phi_deg = df.iloc[i]['phi1']  # Use original phi for calculations
            c_val = df.iloc[i]['c1']      # Use original c for calculations
            
            phi_rad = np.radians(phi_deg)
            cos_phi = np.cos(phi_rad)
            sin_phi = np.sin(phi_rad)

            # tau_ff for the two envelopes: d-psi (Kc=1) and c'-phi' (Kc=Kf).
            tau_ff_k1 = d_val + sigma_fc_i * np.tan(np.radians(psi_val))  # d-psi curve (Kc=1)
            tau_ff_kf = c_val + sigma_fc_i * np.tan(phi_rad)             # c-phi curve (Kc=Kf)

            # Fall back to the lower of the two curves (the doc's negative-stress rule)
            # whenever the K1/Kf interpolation is ill-conditioned: cos(phi) ~ 0, the Kf
            # denominator factor ~ 0, or a non-positive minor principal stress sigma'_3c
            # on either envelope (eqs 7, 8). Previously these cases hit `continue` and the
            # slice silently kept its drained strength in the undrained Stage-2 solve.
            #
            # sigma'_3c (eq 7) is also the DENOMINATOR of K1 (eq 4), so the test must
            # exclude zero, not just negatives -- the source says "negative (or zero)".
            # At exactly zero K1 is +inf and tau_ff is NaN, which max(0.0, NaN) would
            # silently return as 0.0.
            kf_first = sigma_fc_i - c_val * cos_phi   # Kf denominator factor (eq 6)
            use_fallback = abs(cos_phi) < 1e-12 or abs(kf_first) < 1e-12
            if not use_fallback:
                sigma3_k1 = sigma_fc_i + tau_fc_i * (sin_phi - 1) / cos_phi          # eq (7)
                sigma3_kf = kf_first * (1 - sin_phi) / (cos_phi ** 2)                # eq (8)
                use_fallback = sigma3_k1 <= 0 or sigma3_kf <= 0

            if use_fallback:
                tau_ff = min(tau_ff_k1, tau_ff_kf)
            else:
                K1 = (sigma_fc_i + tau_fc_i * (sin_phi + 1) / cos_phi) / \
                     (sigma_fc_i + tau_fc_i * (sin_phi - 1) / cos_phi)               # eq (4)
                Kf = ((sigma_fc_i + c_val * cos_phi) * (1 + sin_phi)) / \
                     (kf_first * (1 - sin_phi))                                       # eq (6)
                if abs(Kf - 1) < 1e-12:
                    tau_ff = tau_ff_k1
                else:
                    tau_ff = ((Kf - K1) * tau_ff_k1 + (K1 - 1) * tau_ff_kf) / (Kf - 1)  # eq (5)

            # The Kc=1 and Kc=Kf envelopes bound the physically possible states, so
            # eq (5) interpolates and never extrapolates. The Stage-1 FS >= 1 guard above
            # already assures K1 <= Kf; a non-finite tau_ff would mean that reasoning
            # failed, and must not be laundered into 0.0 by the clamp below.
            if not np.isfinite(tau_ff):
                return False, (
                    f"Rapid drawdown: undrained strength is not finite for slice {i+1} "
                    f"(sigma'_fc={sigma_fc_i:.4g}, tau_fc={tau_fc_i:.4g}). The K1/Kf "
                    "interpolation is degenerate."
                )

            tau_ff = max(0.0, tau_ff)   # undrained shear strength cannot be negative

            if debug_level >= 2:
                print(f"  Slice {i+1}: tau_ff_k1={tau_ff_k1:.3f}, tau_ff_kf={tau_ff_kf:.3f}, "
                      f"fallback={use_fallback} -> tau_ff={tau_ff:.3f}")

            # Set undrained strength parameters
            df.iloc[i, df.columns.get_loc('c')] = float(tau_ff)
            df.iloc[i, df.columns.get_loc('phi')] = 0.0
        else:
            # High-K material - keep original c and phi
            if debug_level >= 2:
                print(f"Slice {i+1}: High-K material, keeping original c and phi")
    
    # Calculate Stage 2 FS
    success, result_stage2 = method_func(df)
    if not success:
        return False, f"Stage 2 failed: {result_stage2}"
    
    stage2_FS = result_stage2['FS']
    if debug_level >= 1:
        print(f"Stage 2 FS = {stage2_FS:.4f}")

    stage2_state = df.copy()   # per-slice n_eff/u/c/phi for the Stage-2 result

    # Stage 3: Check drained strengths
    if debug_level >= 1:
        print("Stage 3: Checking drained strengths...")
    
    # Check if any low-K slices need drained strength
    need_stage3 = False
    
    for i in range(len(df)):
        d_val = df.iloc[i]['d']
        psi_val = df.iloc[i]['psi']
        
        if d_val > 0 or psi_val > 0:
            # This is a low-K material slice
            if 'n_eff' not in df.columns:
                return False, "Stage 2 did not compute n_eff values"
            
            # Calculate drained strength using equations (9) and (10)
            sigma_prime = df.iloc[i]['n_eff'] / df.iloc[i]['dl']  # Equation (9)
            tau_drained = df.iloc[i]['c1'] + sigma_prime * np.tan(np.radians(df.iloc[i]['phi1']))  # Equation (10)
            
            # Compare with undrained strength (current c value)
            tau_undrained = df.iloc[i]['c']
            
            if debug_level >= 2:
                print(f"Slice {i+1}: tau_drained = {tau_drained:.4f}, tau_undrained = {tau_undrained:.4f}")
            
            if tau_drained < tau_undrained:
                # Use drained strength
                df.iloc[i, df.columns.get_loc('c')] = float(df.iloc[i]['c1'])
                df.iloc[i, df.columns.get_loc('phi')] = float(df.iloc[i]['phi1'])
                need_stage3 = True
                
                if debug_level >= 2:
                    print(f"  Using drained strength for slice {i+1}")
    
    if need_stage3:
        if debug_level >= 1:
            print("Stage 3: Recalculating FS with drained strengths...")
        
        success, result_stage3 = method_func(df)
        if not success:
            return False, f"Stage 3 failed: {result_stage3}"
        
        stage3_FS = result_stage3['FS']
        if debug_level >= 1:
            print(f"Stage 3 FS = {stage3_FS:.4f}")
    else:
        stage3_FS = stage2_FS
        result_stage3 = result_stage2
        if debug_level >= 1:
            print("Stage 3: No drained strength adjustments needed")

    stage3_state = df.copy()   # per-slice state after Stage 3 (== Stage 2 if skipped)

    # Final FS is the lower of Stage 2 and Stage 3
    if stage2_FS < stage3_FS:
        final_FS = stage2_FS
        result = result_stage2
        winning_state = stage2_state
    else:
        final_FS = stage3_FS
        result = result_stage3
        winning_state = stage3_state

    # Hand the winning stage's per-slice results back to the caller's DataFrame so
    # the cached/plotted slices carry real rapid-drawdown stresses (not Stage-0).
    for col in ("n_eff", "u", "c", "phi"):
        if col in winning_state.columns and col in caller_df.columns:
            caller_df[col] = winning_state[col].values
    
    if debug_level >= 1:
        print(f"Final rapid drawdown FS = {final_FS:.4f}")
        print("=== END RAPID DRAWDOWN ANALYSIS ===")
    
    # Append stage FS to result
    result['stage1_FS'] = stage1_FS
    result['stage2_FS'] = stage2_FS
    result['stage3_FS'] = stage3_FS

    return True, result


def reliability(slope_data, method, rapid=False, circular=True, debug_level=0,
                progress_callback=None, cancel_check=None,
                fs_tol=None, tol=None, max_iter=None, composite=False, seed='circles',
                search=True):
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

    Returns:
        tuple: (success, result) where result contains reliability analysis results
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

    def _progress(done, total, label):
        if progress_callback is not None:
            try:
                progress_callback(done, total, label)
            except Exception:
                pass

    start_time = time.time()

    # Validate that at least one material has non-zero standard deviations
    has_std = any(
        m.get('sigma_gamma', 0) != 0 or m.get('sigma_c', 0) != 0 or
        m.get('sigma_phi', 0) != 0 or m.get('sigma_cp', 0) != 0
        for m in slope_data['materials']
    )
    if not has_std:
        return False, ("Reliability analysis requires standard deviations for at least one "
                        "material property (columns L-Q in the mat sheet). None were provided.")

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
        if circular:
            ok_, res_ = generate_slices(sd_, circle=sd_['circles'][0],
                                        num_slices=40, composite=composite)
        else:
            ok_, res_ = generate_slices(sd_, non_circ=sd_['non_circ'], num_slices=40)
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
        return False, ("Reliability: the standard deviation exceeds the mean (COV > 100%) for "
                       f"{details}. mean - sigma is negative, which is non-physical. Reduce the "
                       "standard deviation(s) so mean - sigma >= 0.")

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
    if COV_F == 0:
        return False, "COV_F is zero - no parameter variability"
    
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

    Raises:
        ValueError: on a strength model with no defined perturbation set, rather than
            silently falling back to Mohr-Coulomb's.
    """
    option = material.get('option')
    mapping = {'gamma': 'sigma_gamma'}
    if option == 'mc':
        mapping.update({'c': 'sigma_c', 'phi': 'sigma_phi'})
    elif option == 'cp':
        mapping.update({'c': 'sigma_c', 'cp': 'sigma_cp'})
    elif option:
        raise ValueError(
            f"Reliability analysis does not support the strength option "
            f"option='{option}' on material '{mat_name}'. Supported: mc, cp."
        )
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
        return None, ("Reliability analysis requires standard deviations for at least one "
                      "material property (columns L-Q in the mat sheet). None were provided.")

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
        return None, ("Reliability: the standard deviation exceeds the mean (COV > 100%) for "
                      f"{details}. mean - sigma is negative, which is non-physical. Reduce the "
                      "standard deviation(s) so mean - sigma >= 0.")
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
    if COV_F == 0:
        return "COV_F is zero - no parameter variability"
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
                    cancel_check=None):
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
    ssrm_kw = dict(tolerance=tolerance, grid=tolerance,
                   failure_criterion=failure_criterion,
                   max_iterations=max_iterations, debug_level=max(0, debug_level - 1))

    def _fs(sd, fmin, fmax):
        res = solve_ssrm(build_fem_data(sd, mesh), F_min=fmin, F_max=fmax,
                         cancel_check=cancel_check, **ssrm_kw)
        if not res.get('converged'):
            return None, None, res.get('error', 'SSRM did not converge')
        return res['FS'], res, None

    if debug_level >= 1:
        print("=== RELIABILITY ANALYSIS (FEM / SSRM) ===")

    total_steps = 1 + 2 * len(param_info)
    _progress(0, total_steps, "Solving SSRM at most-likely values…")
    F_MLV, mlv_solution, err = _fs(slope_data, F_min, F_max)
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