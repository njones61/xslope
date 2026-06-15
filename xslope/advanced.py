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
    
    # Import solve module and get the method function
    from . import solve
    method_func = getattr(solve, method_name)

    # Work on a copy: the analysis overwrites strength and load columns (Stage 2
    # swaps in the drawdown pore pressures / loads and undrained strengths), so do
    # not mutate the caller's slice DataFrame.
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
            # denominator factor ~ 0, or a negative minor principal stress sigma'_3c on
            # either envelope (eqs 7, 8). Previously these cases hit `continue` and the
            # slice silently kept its drained strength in the undrained Stage-2 solve.
            kf_first = sigma_fc_i - c_val * cos_phi   # Kf denominator factor (eq 6)
            use_fallback = abs(cos_phi) < 1e-12 or abs(kf_first) < 1e-12
            if not use_fallback:
                sigma3_k1 = sigma_fc_i + tau_fc_i * (sin_phi - 1) / cos_phi          # eq (7)
                sigma3_kf = kf_first * (1 - sin_phi) / (cos_phi ** 2)                # eq (8)
                use_fallback = sigma3_k1 < 0 or sigma3_kf < 0

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
    
    # Final FS is the lower of Stage 2 and Stage 3
    if stage2_FS < stage3_FS:
        final_FS = stage2_FS
        result = result_stage2
    else:
        final_FS = stage3_FS
        result = result_stage3
    
    if debug_level >= 1:
        print(f"Final rapid drawdown FS = {final_FS:.4f}")
        print("=== END RAPID DRAWDOWN ANALYSIS ===")
    
    # Append stage FS to result
    result['stage1_FS'] = stage1_FS
    result['stage2_FS'] = stage2_FS
    result['stage3_FS'] = stage3_FS

    return True, result


def reliability(slope_data, method, rapid=False, circular=True, debug_level=0):
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
    
    Returns:
        tuple: (success, result) where result contains reliability analysis results
    """
    
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
    
    # Step 1: Find the critical failure surface using search
    if circular:
        if debug_level >= 1:
            print("Performing circular search...")
        fs_cache, converged, search_path, circle_cache = circular_search(slope_data, method, rapid=rapid)
    else:
        if debug_level >= 1:
            print("Performing noncircular search...")
        fs_cache, converged, search_path = noncircular_search(slope_data, method, rapid=rapid)
    
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

        # Perturb only the strength parameters the material's strength model actually
        # uses: 'mc' (Mohr-Coulomb) uses c and phi; 'cp' (undrained strength
        # Su = c + cp*max(0, r_elev - y)) uses both the base c and the rate cp. gamma
        # applies to both models. Perturbing a field the model does not use (e.g. phi on
        # a cp material) would do nothing and silently drop that material's strength
        # uncertainty from COV_F.
        if material.get('option') == 'cp':
            param_mappings = {'gamma': 'sigma_gamma', 'c': 'sigma_c', 'cp': 'sigma_cp'}
        else:
            param_mappings = {'gamma': 'sigma_gamma', 'c': 'sigma_c', 'phi': 'sigma_phi'}

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
    delta_F_values = []

    for i, param in enumerate(param_info):
        if debug_level >= 1:
            print(f"\nProcessing parameter {i+1}/{len(param_info)}: Material {param['material_id']}, {param['param']}")
        
        # Create modified slope_data copies
        slope_data_plus = slope_data.copy()
        slope_data_minus = slope_data.copy()
        slope_data_plus['materials'] = [mat.copy() for mat in materials]
        slope_data_minus['materials'] = [mat.copy() for mat in materials]
        
        # Find the material and modify the parameter (use 0-based index)
        mat_index = param['material_id'] - 1
        if mat_index < len(slope_data_plus['materials']):
            slope_data_plus['materials'][mat_index][param['param']] = param['mlv'] + param['std']
        
        if mat_index < len(slope_data_minus['materials']):
            slope_data_minus['materials'][mat_index][param['param']] = param['mlv'] - param['std']
        
        # Calculate F+ and F-
        if circular:
            fs_cache_plus, _, _, _ = circular_search(slope_data_plus, method, rapid=rapid)
            fs_cache_minus, _, _, _ = circular_search(slope_data_minus, method, rapid=rapid)
        else:
            fs_cache_plus, _, _ = noncircular_search(slope_data_plus, method, rapid=rapid)
            fs_cache_minus, _, _ = noncircular_search(slope_data_minus, method, rapid=rapid)
        
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