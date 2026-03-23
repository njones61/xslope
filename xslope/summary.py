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

"""Summary printing functions for LEM analysis results."""

import numpy as np


def print_ito_matsui_summary(slope_data, slice_df):
    """Print Ito & Matsui pile force summary for the critical failure surface.

    For each pile where H was auto-computed, prints force details, coefficients,
    and structural capacity check results (if V_cap/M_cap are provided).
    """
    if slice_df is None or 'h_pile' not in slice_df.columns:
        return
    pile_slices = slice_df[slice_df['h_pile'] > 0]
    if pile_slices.empty or not slope_data.get('pile_lines'):
        return

    from .ito_matsui import (ito_matsui_coefficients,
                             intersect_pile_with_materials,
                             compute_ito_matsui_force,
                             compute_ito_matsui_force_and_moment_arm)

    for pile_info, (_, ps) in zip(slope_data['pile_lines'], pile_slices.iterrows()):
        if pile_info.get('H') is not None:
            continue
        if not pile_info.get('D_pile') or not pile_info.get('S'):
            continue

        D_p, S_p = pile_info['D_pile'], pile_info['S']
        D1 = S_p - D_p
        gs_coords = np.array(slope_data['ground_surface'].coords)
        y_gnd = np.interp(ps['x_pile'], gs_coords[:, 0], gs_coords[:, 1])
        depth = y_gnd - ps['y_pile']
        segments = intersect_pile_with_materials(
            ps['x_pile'], y_gnd, ps['y_pile'],
            slope_data['profile_lines'], slope_data['materials'])
        H_val, F_pile = compute_ito_matsui_force(D_p, S_p, segments)
        A1, A2 = ito_matsui_coefficients(D1, D_p, segments[0]['phi']) if segments else (0, 0)

        label = pile_info.get('label', f'x={ps["x_pile"]:.2f}')
        print(f'  === Ito & Matsui Summary ({label}) ===')
        print(f'  Pile diameter (D)          = {D_p:.1f}')
        print(f'  Pile spacing (S)           = {S_p:.1f}')
        print(f'  Clear spacing (D1 = S - D) = {D1:.1f}')
        print(f'  Depth to failure surface   = {depth:.1f}')
        print(f'  Coefficients: A1 = {A1:.3f}, A2 = {A2:.3f}')
        print(f'  Force per pile (F_pile)    = {F_pile:.0f}')
        print(f'  Force per unit width (H)   = {H_val:.1f}')

        # Structural capacity check
        V_cap = pile_info.get('V_cap')
        M_cap = pile_info.get('M_cap')
        if V_cap is not None or M_cap is not None:
            _, _, L_m = compute_ito_matsui_force_and_moment_arm(D_p, S_p, segments)
            print(f'  --- Structural Capacity Check ---')
            F_capped = F_pile
            if V_cap is not None:
                governs = F_pile > V_cap
                print(f'  V_cap = {V_cap:.0f}  {"[GOVERNS]" if governs else "[OK]"}')
                F_capped = min(F_capped, V_cap)
            if M_cap is not None:
                F_from_moment = M_cap / L_m if L_m > 0 else float('inf')
                governs = F_pile > F_from_moment
                print(f'  M_cap = {M_cap:.0f}, L_m = {L_m:.2f}, '
                      f'M_cap/L_m = {F_from_moment:.0f}  {"[GOVERNS]" if governs else "[OK]"}')
                if L_m > 0:
                    F_capped = min(F_capped, F_from_moment)
            if F_capped < F_pile:
                H_capped = F_capped / S_p
                print(f'  Capped F_pile = {F_capped:.0f}, Capped H = {H_capped:.1f}')


def print_rapid_drawdown_summary(results):
    """Print rapid drawdown stage results for the critical surface."""
    if results and 'stage1_FS' in results:
        print(f"=== RAPID DRAWDOWN SUMMARY (Critical Surface) ===")
        print(f"Stage 1 FS = {results['stage1_FS']:.4f}")
        print(f"Stage 2 FS = {results['stage2_FS']:.4f}")
        print(f"Stage 3 FS = {results['stage3_FS']:.4f}")
        print(f"Final rapid drawdown FS = {results['FS']:.4f}")


def print_no_solution_warning():
    """Print warning when no valid solution was found (FS=9999 for all surfaces)."""
    from .slice import _ito_matsui_warned
    print(f"[WARNING] No valid solution found (FS=9999 for all surfaces).")
    if _ito_matsui_warned:
        print(f"    Ito & Matsui pile forces exceed driving forces for all trial surfaces.")
        print(f"    This may mean the piles are sufficient to stabilize the slope. Note that")
        print(f"    Ito & Matsui computes an upper bound on soil resistance and does not check")
        print(f"    pile structural capacity. The available resistance should be taken as the")
        print(f"    lesser of the Ito & Matsui soil resistance and the pile shear/moment")
        print(f"    capacity. If the structural capacity governs, specify H directly in the")
        print(f"    piles sheet as H = min(structural capacity, Ito & Matsui) / S.")
    else:
        print(f"    Check input parameters.")
