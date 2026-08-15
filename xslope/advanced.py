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

import numpy as np


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
    if 'beta2' in df.columns:
        # The stage-2 load carries its OWN inclination: the two dloads sheets declare
        # their directions independently, so a stage-1 vertical surcharge beside a
        # stage-2 normal pool (or the reverse) must not inherit the other's beta.
        df['beta'] = df['beta2']
    
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


# ---------------------------------------------------------------------------
# The reliability family moved to :mod:`xslope.reliability`. It is re-exported
# here so every historical import path keeps working unchanged, e.g.
#   from xslope.advanced import reliability, reliability_mc, reliability_fem
# Prefer importing these from xslope.reliability in new code.
# ---------------------------------------------------------------------------
from .reliability import (  # noqa: F401  (re-exported for backward compatibility)
    reliability,
    reliability_taylor,
    reliability_mc,
    reliability_rs,
    reliability_fem,
    MC_DEFAULT_SEED,
    MC_DEFAULT_SAMPLES,
    _MC_FLOOR,
    _strength_param_mapping,
    _reliability_param_info,
    _perturbed_slope_data,
    _finalize_reliability,
    _mc_param_info,
    _mc_sampled_slope_data,
)
