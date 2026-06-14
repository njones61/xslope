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

from math import sin, cos, tan, radians, atan, atan2, degrees

import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar, root_scalar, newton
from shapely.geometry import LineString, Point
from tabulate import tabulate

from .advanced import rapid_drawdown

def solve_selected(method_name, slice_df, rapid=False):
    """
    Executes a specified limit equilibrium solution method and displays results.

    Parameters
    ----------
    method_name : str
        Name of the solution method function to call. Must be one of:
        'oms', 'bishop', 'janbu', 'spencer', 'corps_engineers', 'lowe_karafiath'
    slice_df : pandas.DataFrame
        Slice dataframe containing all required columns for the specified method
        (see individual method documentation for column requirements)
    rapid : bool, optional
        If True, performs rapid drawdown analysis using the specified method.
        Default is False.

    Returns
    -------
    dict or str
        If successful: dictionary containing method results (includes 'FS' and method-specific parameters)
        If failed: error message string

    Notes
    -----
    This function automatically prints the factor of safety and method-specific
    parameters to the console. For methods with additional parameters:
    - Spencer: displays theta (interslice force angle)
    - Janbu: displays fo (correction factor)
    - Corps of Engineers: displays theta
    """

    func = globals()[method_name]

    if rapid:
        success, result = rapid_drawdown(slice_df, method_name)
    else:
        success, result = func(slice_df)
    if not success:
        print(f'Error: {result}')
        return result

    if func == oms:
        print(f'OMS: FS={result["FS"]:.3f}')
    elif func == bishop:
        print(f'Bishop: FS={result["FS"]:.3f}')
    elif func == spencer:
        print(f'Spencer: FS={result["FS"]:.3f}, theta={result["theta"]:.2f}')
    elif func == janbu:
        print(f'Janbu Corrected FS={result["FS"]:.3f}, fo={result["fo"]:.2f}')
    elif func == corps_engineers:
        print(f'Corps Engineers: FS={result["FS"]:.3f}, theta={result["theta"]:.2f}')
    elif func == lowe_karafiath:
        print(f'Lowe & Karafiath: FS={result["FS"]:.3f}')

    return result

def solve_all(slice_df, rapid=False):
    """
    Executes all available limit equilibrium solution methods sequentially.

    Runs six different limit equilibrium methods on the provided slice dataframe
    and displays the factor of safety for each method. This is useful for comparing
    results across multiple solution approaches.

    Parameters
    ----------
    slice_df : pandas.DataFrame
        Slice dataframe containing all required columns for all methods.
        Must include: 'alpha', 'phi', 'c', 'w', 'u', 'dl', 'dload', 'd_x', 'd_y',
        'beta', 'kw', 't', 'y_t', 'p', 'x_c', 'y_cg', and additional columns
        required for specific methods (e.g., 'r', 'xo', 'yo' for circular methods).

    Returns
    -------
    None
        Results are printed to console for each method.

    Notes
    -----
    Methods executed in order:
    1. Ordinary Method of Slices (OMS)
    2. Bishop's Simplified Method
    3. Janbu's Simplified Method
    4. Corps of Engineers Method
    5. Lowe & Karafiath Method
    6. Spencer's Method

    If any method fails, an error message is displayed but execution continues
    with the remaining methods.
    """
    solve_selected('oms', slice_df, rapid=rapid)
    solve_selected('bishop', slice_df, rapid=rapid)
    solve_selected('janbu', slice_df, rapid=rapid)
    solve_selected('corps_engineers', slice_df, rapid=rapid)
    solve_selected('lowe_karafiath', slice_df, rapid=rapid)
    solve_selected('spencer', slice_df, rapid=rapid)

def oms(slice_df, debug=False):
    """
    Computes FS by direct application of Equation 9 (Ordinary Method of Slices).

    Inputs
    ------
    slice_df : pandas.DataFrame
        Must contain exactly these columns (length = n slices):
          'alpha'   (deg)   = base inclination αᵢ
          'phi'     (deg)   = friction angle φᵢ
          'c'             = cohesion cᵢ
          'w'             = slice weight Wᵢ
          'u'             = pore pressure force/unit‐length on base, uᵢ
          'dl'            = base length Δℓᵢ
          'dload'         = resultant distributed load Dᵢ
          'd_x','d_y'     = centroid (x,y) at which Dᵢ acts
          'beta'   (deg)   = distributed-load inclination βᵢ
          'kw'            = seismic horizontal kWᵢ
          't'             = tension‐crack horizontal Tᵢ  (zero except one slice)
          'y_t'           = y‐loc of Tᵢ's line of action (zero except that one slice)
          'p'             = reinforcement uplift pᵢ (zero if none)
          'x_c','y_cg'    = slice‐centroid (x,y) for seismic moment arm
          'h_pile'        = pile/pier force on base (zero if none)
          'theta_p' (RAD)  = pile force inclination from horizontal — stored and
                            consumed in RADIANS, unlike every other angle here
          'r'             = radius of circular failure surface
          'xo','yo'       = x,y coordinates of circle center

    Returns
    -------
    (bool, dict_or_str)
      • If success: (True, {'method':'oms', 'FS': <computed value>})
      • If denominator → 0 or other fatal error: (False, "<error message>")


    """
    if 'r' not in slice_df.columns:
        return False, "Circle is required for OMS method."

    # 1) Unpack circle‐center and radius as single values
    Xo = slice_df['xo'].iloc[0]    # Xoᵢ (x-coordinate of circle center)
    Yo = slice_df['yo'].iloc[0]    # Yoᵢ (y-coordinate of circle center)
    R  = slice_df['r'].iloc[0]     # Rᵢ (radius of circular failure surface)

    # 2) Pull arrays directly from slice_df
    alpha_deg = slice_df['alpha'].values    # αᵢ in degrees
    phi_deg   = slice_df['phi'].values      # φᵢ in degrees
    c     = slice_df['c'].values        # cᵢ
    W     = slice_df['w'].values        # Wᵢ
    u     = slice_df['u'].values        # uᵢ (pore‐force per unit length)
    dl     = slice_df['dl'].values       # Δℓᵢ
    D     = slice_df['dload'].values        # Dᵢ
    d_x    = slice_df['d_x'].values      # d_{x,i}
    d_y    = slice_df['d_y'].values      # d_{y,i}
    beta_deg  = slice_df['beta'].values     # βᵢ in degrees
    kw    = slice_df['kw'].values       # kWᵢ
    T     = slice_df['t'].values        # Tᵢ (zero except one slice)
    y_t    = slice_df['y_t'].values      # y_{t,i} (zero except one slice)
    P     = slice_df['p'].values        # pᵢ
    x_c    = slice_df['x_c'].values      # x_{c,i}
    y_cg = slice_df['y_cg'].values       # y_{cg,i} coordinate of slice centroid

    # Pile forces (backward compatible — zeros if no pile data)
    n = len(slice_df)
    H_pile  = slice_df['h_pile'].values  if 'h_pile'  in slice_df.columns else np.zeros(n)
    theta_p = slice_df['theta_p'].values if 'theta_p' in slice_df.columns else np.zeros(n)
    y_pile  = slice_df['y_pile'].values  if 'y_pile'  in slice_df.columns else np.zeros(n)
    x_pile  = slice_df['x_pile'].values  if 'x_pile'  in slice_df.columns else np.zeros(n)

    # 3) Convert angles to radians
    alpha = np.radians(alpha_deg)   # αᵢ [rad]
    phi   = np.radians(phi_deg)     # φᵢ [rad]
    beta  = np.radians(beta_deg)    # βᵢ [rad]

    # 4) Precompute sines/cosines
    sin_alpha = np.sin(alpha)          # sin(αᵢ)
    cos_alpha = np.cos(alpha)          # cos(αᵢ)
    sin_ab    = np.sin(alpha - beta)   # sin(αᵢ−βᵢ)
    cos_ab    = np.cos(alpha - beta)   # cos(αᵢ−βᵢ)
    tan_phi   = np.tan(phi)            # tan(φᵢ)

    # ————————————————————————————————————————————————————————
    # 5) Build the NUMERATOR (Eqn 8) = Σᵢ [  cᵢ·Δℓᵢ
    #       + (Wᵢ·cosαᵢ + Dᵢ·cos(αᵢ−βᵢ) − kWᵢ·sinαᵢ − Tᵢ·sinαᵢ − uᵢ·Δℓᵢ )·tanφᵢ ]
    #


    # N′ᵢ = Wᵢ·cosαᵢ + Dᵢ·cos(αᵢ−βᵢ) − kWᵢ·sinαᵢ − Tᵢ·sinαᵢ − uᵢ·Δℓᵢ + Hᵢ·sin(αᵢ−θₚᵢ)
    N_eff = (
        W * cos_alpha
      + D * cos_ab
      - kw * sin_alpha
      - T * sin_alpha
      - (u * dl)
      + H_pile * np.sin(alpha - theta_p)
    )

    #   Σ  Dᵢ·sinβᵢ·(Yo - d_{y,i}) 
    a_dy = Yo - d_y
    sum_Dy = np.sum(D * np.sin(beta) * a_dy)

    numerator = np.sum(c * dl + N_eff * tan_phi)

    # ————————————————————————————————————————————————————————
    # 6) Build each piece of the DENOMINATOR exactly as Eqn 8:

    #  (A) = Σ [ Wᵢ · sinαᵢ ]
    sum_W = np.sum(W * sin_alpha)

    #  (B) = Σ  Dᵢ·cosβᵢ·(Xo - d_{x,i})
    # The vertical D-component (D·cosβ) acts at horizontal arm (d_x - Xo). That arm
    # flips sign when the slope faces right, and cosβ (unlike sinβ) is NOT corrected
    # by the β sign-flip applied in generate_slices, so flip it here. The horizontal
    # D-component term (a_dy / sum_Dy) is already sign-correct via that β-flip and is
    # left untouched; likewise the seismic (a_s) and tension-crack (a_t) arms.
    right_facing = slice_df['y_lt'].iat[0] > slice_df['y_rt'].iat[-1]
    a_dx = d_x - Xo
    if right_facing:
        a_dx = -a_dx
    sum_Dx = np.sum(D * np.cos(beta) * a_dx)

    #  (C) = Σ [ kWᵢ * (Yo - y_{cg,i}) ]
    a_s = Yo - y_cg
    sum_kw = np.sum(kw * a_s)

    #  (D) = Σ [ Tᵢ * (Yo - y_{t,i}) ]
    a_t = Yo - y_t
    sum_T = np.sum(T * a_t)

    # Pile resisting moment about circle center:
    # M_pile = Σ [ H·cos(θₚ)·(Yo - y_pile) + H·sin(θₚ)·(x_pile - Xo) ]
    # The vertical-component arm (x_pile - Xo) is in real coordinates while alpha runs in
    # the orientation-normalized frame, so it must flip on right-facing slopes (same
    # correction as the distributed-load a_dx arm). Only matters for battered piles
    # (theta_p != 0); for theta_p = 0 this term is zero.
    pile_x_arm = -(x_pile - Xo) if right_facing else (x_pile - Xo)
    sum_pile_moment = np.sum(
        H_pile * np.cos(theta_p) * (Yo - y_pile)
      + H_pile * np.sin(theta_p) * pile_x_arm
    )

    # Put them together with their 1/R factors:
    # P, sum_Dy, and pile moment are known resisting forces, not factored by FS
    denominator = sum_W + (1.0 / R) * (sum_Dx + sum_kw + sum_T) - np.sum(P) - (1.0 / R) * sum_Dy - (1.0 / R) * sum_pile_moment

    # 7) Finally compute FS = (numerator)/(denominator)
    if denominator <= 0:
        return False, "Net driving moment is non-positive (resisting forces exceed driving forces)."
    FS = numerator / denominator

    # 8) Store effective normal forces in the DataFrame
    slice_df['n_eff'] = N_eff

    if debug==True:
        print(f'numerator = {numerator:.4f}')
        print(f'denominator = {denominator:.4f}')
        print(f'Sum_W = {sum_W:.4f}')
        print(f'Sum_Dx = {sum_Dx:.4f}')
        print(f'Sum_Dy = {sum_Dy:.4f}')
        print(f'Sum_kw = {sum_kw:.4f}')
        print(f'Sum_T = {sum_T:.4f}')
        if np.any(H_pile > 0):
            print(f'sum_pile_moment = {sum_pile_moment:.4f}')
        print('N_eff =', np.array2string(N_eff, precision=4, separator=', '))

    # 9) Return success and the FS
    return True, {'method': 'oms', 'FS': FS}

def bishop(slice_df, debug=False, tol=1e-6, max_iter=100):
    """
    Computes FS using the complete Bishop's Simplified Method (Equation 10) and computes N_eff (Equation 8).
    Requires circular slip surface and full input data structure consistent with OMS.

    Parameters:
        slice_df : pandas.DataFrame with required columns (see OMS spec)
        debug : bool, if True prints diagnostic info
        tol : float, convergence tolerance
        max_iter : int, maximum iteration steps

    Returns:
        (bool, dict | str): (True, {'method': 'bishop', 'FS': value}) or (False, error message)
    """

    if 'r' not in slice_df.columns:
        return False, "Circle is required for Bishop method."

    # 1) Unpack circle‐center and radius as single values
    Xo = slice_df['xo'].iloc[0]    # Xoᵢ (x-coordinate of circle center)
    Yo = slice_df['yo'].iloc[0]    # Yoᵢ (y-coordinate of circle center)
    R  = slice_df['r'].iloc[0]     # Rᵢ (radius of circular failure surface)

    # Load input arrays
    alpha = np.radians(slice_df['alpha'].values)
    phi   = np.radians(slice_df['phi'].values)
    c     = slice_df['c'].values
    W     = slice_df['w'].values
    u     = slice_df['u'].values
    dl    = slice_df['dl'].values
    D     = slice_df['dload'].values
    d_x   = slice_df['d_x'].values
    d_y   = slice_df['d_y'].values
    beta  = np.radians(slice_df['beta'].values)
    kw    = slice_df['kw'].values
    T     = slice_df['t'].values
    y_t   = slice_df['y_t'].values
    P     = slice_df['p'].values
    x_c   = slice_df['x_c'].values
    y_cg  = slice_df['y_cg'].values

    # Pile forces (backward compatible — zeros if no pile data)
    n = len(slice_df)
    H_pile  = slice_df['h_pile'].values  if 'h_pile'  in slice_df.columns else np.zeros(n)
    theta_p = slice_df['theta_p'].values if 'theta_p' in slice_df.columns else np.zeros(n)
    y_pile  = slice_df['y_pile'].values  if 'y_pile'  in slice_df.columns else np.zeros(n)
    x_pile  = slice_df['x_pile'].values  if 'x_pile'  in slice_df.columns else np.zeros(n)

    # Trigonometric terms
    sin_alpha = np.sin(alpha)
    cos_alpha = np.cos(alpha)
    tan_phi   = np.tan(phi)
    sin_beta  = np.sin(beta)
    cos_beta  = np.cos(beta)

    # Moment arms. The vertical-D-component arm (a_dx) flips sign on right-facing
    # slopes and is not corrected by the β sign-flip in generate_slices (cosβ is
    # even in β), so flip it here. a_dy (horizontal-D-component) is already correct
    # via that β-flip; a_s/a_t track the sliding-direction convention. See oms().
    right_facing = slice_df['y_lt'].iat[0] > slice_df['y_rt'].iat[-1]
    a_dx = d_x - Xo
    if right_facing:
        a_dx = -a_dx
    a_dy = Yo - d_y
    a_s  = Yo - y_cg
    a_t  = Yo - y_t

    # Pile resisting moment about circle center (same term as OMS) — flip the
    # vertical-component arm (x_pile - Xo) on right-facing slopes, as for a_dx above.
    pile_x_arm = -(x_pile - Xo) if right_facing else (x_pile - Xo)
    sum_pile_moment = np.sum(
        H_pile * np.cos(theta_p) * (Yo - y_pile)
      + H_pile * np.sin(theta_p) * pile_x_arm
    )

    # Denominator (moment equilibrium) — pile moment is a known resisting force, not factored by FS
    sum_W = np.sum(W * sin_alpha)
    sum_Dx = np.sum(D * cos_beta * a_dx)
    sum_Dy = np.sum(D * sin_beta * a_dy)
    sum_kw = np.sum(kw * a_s)
    sum_T = np.sum(T * a_t)
    denominator = sum_W + (1.0 / R) * (sum_Dx + sum_kw + sum_T) - np.sum(P) - (1.0 / R) * sum_Dy - (1.0 / R) * sum_pile_moment

    if denominator <= 0:
        return False, "Net driving moment is non-positive (resisting forces exceed driving forces)."

    # Vertical component of pile force (upward, reduces net downward force on slice)
    H_sin_tp = H_pile * np.sin(theta_p)

    # Iterative solution
    F = 1.0
    for _ in range(max_iter):
        # Compute N_eff — H·sin(θp) enters vertical equilibrium
        num_N = (
            W + D * cos_beta - P * sin_alpha - H_sin_tp
            - u * dl * cos_alpha
            - (c * dl * sin_alpha) / F
        )
        denom_N = cos_alpha + (sin_alpha * tan_phi) / F
        N_eff = num_N / denom_N

        # Numerator for FS — same H·sin(θp) term in the weight group
        shear = (
            c * dl * cos_alpha
            + (W + D * cos_beta - P * sin_alpha - H_sin_tp - u * dl * cos_alpha) * tan_phi
        )
        numerator = np.sum(shear / denom_N)
        F_new = numerator / denominator

        if abs(F_new - F) < tol:
            slice_df['n_eff'] = N_eff
            if debug:
                print(f"FS = {F_new:.6f}")
                print(f"Numerator = {numerator:.6f}")
                print(f"Denominator = {denominator:.6f}")
                if np.any(H_pile > 0):
                    print(f"sum_pile_moment = {sum_pile_moment:.6f}")
                print("N_eff =", np.array2string(N_eff, precision=4, separator=', '))
            return True, {'method': 'bishop', 'FS': F_new}

        F = F_new

    return False, "Bishop method did not converge within the maximum number of iterations."

def janbu(slice_df, debug=False, tol=1e-6, max_iter=100):
    """
    Computes FS using Janbu's Simplified Method with the f0 correction factor.

    Janbu's Simplified Method satisfies overall HORIZONTAL FORCE equilibrium of
    the sliding mass with the inter-slice shear forces neglected. The effective
    base normal force is obtained from VERTICAL equilibrium of each slice — the
    same m_alpha = cosα + sinα·tanφ/F relation used by Bishop's method — which
    makes the factor of safety iterative in F. The base (uncorrected) factor of
    safety is then multiplied by Janbu's empirical correction factor
    f0(d/L, soil type), which compensates for neglecting the inter-slice shear.

    This is the standard Janbu Simplified formulation (matching e.g. GeoStudio
    SLOPE/W's F_f at lambda=0). It is NOT the Ordinary/Fellenius normal force —
    using N = W·cosα − u·dl with a ΣW·sinα denominator (an earlier xslope error)
    collapses the base FS onto the Ordinary Method of Slices.

    Implements the full slice force budget: distributed loads (dload at the load
    inclination beta), seismic force (kw, horizontal), tension-crack water force
    (t, horizontal), line reinforcement (p), and pile forces (h_pile/theta_p).

    Parameters:
        slice_df : pandas.DataFrame with required columns (see OMS spec)
        debug : bool, if True prints diagnostic info
        tol : float, convergence tolerance on F
        max_iter : int, maximum fixed-point iterations

    Returns:
        (bool, dict | str): (True, {'method': 'janbu', 'FS': corrected_value,
                            'fo': correction_factor, 'FS_base': uncorrected_value})
                           or (False, error message)
    """

    # Load input arrays
    alpha = np.radians(slice_df['alpha'].values)
    phi = np.radians(slice_df['phi'].values)
    c = slice_df['c'].values
    W = slice_df['w'].values
    u = slice_df['u'].values
    dl = slice_df['dl'].values
    D = slice_df['dload'].values
    beta = np.radians(slice_df['beta'].values)
    kw = slice_df['kw'].values
    T = slice_df['t'].values
    P = slice_df['p'].values

    # Pile forces (backward compatible — zeros if no pile data)
    n = len(slice_df)
    H_pile  = slice_df['h_pile'].values  if 'h_pile'  in slice_df.columns else np.zeros(n)
    theta_p = slice_df['theta_p'].values if 'theta_p' in slice_df.columns else np.zeros(n)

    # Trigonometric terms
    sin_alpha = np.sin(alpha)
    cos_alpha = np.cos(alpha)
    tan_phi = np.tan(phi)
    cos_beta = np.cos(beta)
    sin_beta = np.sin(beta)

    # Pile force components: vertical H·sin(θp) into vertical equilibrium,
    # horizontal H·cos(θp) into the horizontal force balance.
    H_sin_tp = H_pile * np.sin(theta_p)
    H_cos_tp = H_pile * np.cos(theta_p)

    # External horizontal driving forces (independent of F): seismic (driving),
    # distributed-load horizontal component (driving), tension-crack water
    # (driving), reinforcement and pile horizontal components (resisting).
    horiz_ext = (np.sum(kw) + np.sum(D * sin_beta) + np.sum(T)
                 - np.sum(P) - np.sum(H_cos_tp))

    # Iterate F: the base normal depends on F through m_alpha, exactly as Bishop.
    F = 1.0
    N_eff = None
    converged = False
    for _ in range(max_iter):
        m_alpha = cos_alpha + sin_alpha * tan_phi / F
        # Effective base normal from VERTICAL slice equilibrium (Bishop's N').
        num_N = (W + D * cos_beta - P * sin_alpha - H_sin_tp
                 - u * dl * cos_alpha - (c * dl * sin_alpha) / F)
        N_eff = num_N / m_alpha
        N_total = N_eff + u * dl
        # Horizontal force equilibrium: resisting shear (projected horizontally)
        # over driving (base-normal horizontal component + external horizontal).
        resisting = np.sum((c * dl + N_eff * tan_phi) * cos_alpha)
        driving = np.sum(N_total * sin_alpha) + horiz_ext
        if driving <= 0:
            return False, "Net horizontal driving force is non-positive (resisting forces exceed driving forces)."
        F_new = resisting / driving
        if abs(F_new - F) < tol:
            F = F_new
            converged = True
            break
        F = F_new

    if not converged:
        return False, "Janbu method did not converge within the maximum number of iterations."

    FS_base = F

    # === Compute Janbu correction factor ===

    # Get failure surface endpoints
    x_l = slice_df['x_l'].iloc[0]  # leftmost x
    y_lt = slice_df['y_lt'].iloc[0]  # leftmost top y
    x_r = slice_df['x_r'].iloc[-1]  # rightmost x
    y_rt = slice_df['y_rt'].iloc[-1]  # rightmost top y

    # Length of failure surface (straight line approximation)
    L = np.hypot(x_r - x_l, y_rt - y_lt)

    # Calculate perpendicular distance from each slice center to failure surface line
    x0 = slice_df['x_c'].values
    y0 = slice_df['y_cb'].values

    # Distance from point to line formula: |ax + by + c| / sqrt(a² + b²)
    # Line equation: (y_rt - y_lt)x - (x_r - x_l)y + (x_r * y_lt - y_rt * x_l) = 0
    numerator_dist = np.abs((y_rt - y_lt) * x0 - (x_r - x_l) * y0 + x_r * y_lt - y_rt * x_l)
    dists = numerator_dist / L
    d = np.max(dists)  # maximum perpendicular distance

    dL_ratio = d / L

    # Determine b1 factor by overall soil type (Janbu's three correction-factor
    # curves). The classification is global by design — Janbu's chart is read
    # once for the whole surface — using a small tolerance rather than exact
    # float equality so that, e.g., a φ stored as 1e-12 still reads as φ=0.
    phi_sum = slice_df['phi'].sum()
    c_sum = slice_df['c'].sum()

    if phi_sum < 1e-9:    # c-only soil (undrained, φ = 0)
        b1 = 0.67
    elif c_sum < 1e-9:    # φ-only soil (no cohesion)
        b1 = 0.31
    else:                 # c-φ soil
        b1 = 0.50

    # Correction factor. The polynomial 1 + b1·(d/L − 1.4·(d/L)²) is a fit to
    # Janbu's correction-factor charts and is only valid up to the chart domain;
    # it peaks at d/L = 1/2.8 ≈ 0.357 and turns over (eventually < 1) beyond
    # that, which is unphysical. d/L is geometrically unbounded, so clamp it to
    # the peak: f0 saturates at its maximum rather than decreasing for very deep
    # surfaces. Typical slip surfaces have d/L well below 0.357 (the repo
    # benchmarks are ~0.13–0.20), so this only guards pathological geometries.
    DL_PEAK = 1.0 / 2.8
    dL_eff = min(dL_ratio, DL_PEAK)
    if debug and dL_ratio > DL_PEAK:
        print(f"Janbu f0: d/L = {dL_ratio:.3f} exceeds chart peak {DL_PEAK:.3f}; clamped.")
    fo = 1 + b1 * (dL_eff - 1.4 * dL_eff ** 2)

    # Final corrected factor of safety
    FS = FS_base * fo

    # Store effective normal forces in DataFrame
    slice_df['n_eff'] = N_eff

    if debug:
        print(f"FS_base (uncorrected) = {FS_base:.6f}")
        print(f"d/L ratio = {dL_ratio:.4f}")
        print(f"b1 factor = {b1:.2f}")
        print(f"fo correction = {fo:.4f}")
        print(f"FS_corrected = {FS:.6f}")
        if np.any(H_pile > 0):
            print(f"sum_pile_horizontal = {np.sum(H_cos_tp):.6f}")
        print("N_eff =", np.array2string(N_eff, precision=4, separator=', '))

    return True, {
        'method': 'janbu',
        'FS': FS,
        'fo': fo,
        'FS_base': FS_base,
    }


def force_equilibrium(slice_df, theta_list, fs_guess=1.5, tol=1e-6, max_iter=50, debug=False, right_facing=False):
    """
    Limit‐equilibrium by force equilibrium in X & Y with variable interslice angles.

    Parameters:
        slice_df (pd.DataFrame): must contain columns
            'alpha' (slice base inclination, degrees),
            'phi'   (slice friction angle, degrees),
            'c'     (cohesion),
            'dl'    (slice base length),
            'w'     (slice weight),
            'u'     (pore force per unit length),
            'dload' (distributed load),
            'beta'  (distributed load inclination, degrees),
            'kw'    (seismic force),
            't'     (tension crack water force),
            'p'     (reinforcement force),
            'h_pile' (pile force, optional), 'theta_p' (pile inclination, RADIANS, optional)
        theta_list (array-like): slice‐boundary force inclinations (degrees),
                                 length must be n+1 if there are n slices
        fs_guess (float): initial guess for factor of safety
        tol (float): convergence tolerance on residual
        max_iter (int): maximum number of Newton (secant) iterations
        debug (bool): print residuals during iteration

    Returns:
        (bool, dict or str):
           - If converged: (True, {'method':'force_equilibrium','FS':<value>})
           - If failed:   (False, "error message")
    """
    import numpy as np

    n = len(slice_df)
    if len(theta_list) != n+1:
        return False, f"theta_list length ({len(theta_list)}) must be n+1 ({n+1})"

    # extract and convert to radians
    alpha   = np.radians(slice_df['alpha'].values)
    phi     = np.radians(slice_df['phi'].values)
    c       = slice_df['c'].values
    w       = slice_df['w'].values
    u       = slice_df['u'].values
    dl      = slice_df['dl'].values
    D       = slice_df['dload'].values
    beta    = np.radians(slice_df['beta'].values)
    kw      = slice_df['kw'].values
    T       = slice_df['t'].values
    P       = slice_df['p'].values
    theta   = np.radians(np.asarray(theta_list))

    # Pile forces (backward compatible — zeros if no pile data)
    H_pile  = slice_df['h_pile'].values  if 'h_pile'  in slice_df.columns else np.zeros(n)
    theta_p = slice_df['theta_p'].values if 'theta_p' in slice_df.columns else np.zeros(n)

    N = np.zeros(n)  # normal forces on slice bases
    Z = np.zeros(n+1)  # interslice forces, Z[0] = 0 by definition (no force entering leftmost slice)

    def residual(FS):
        """Return the right‐side interslice force Z[n] for a given FS."""
        c_m       = c / FS
        tan_phi_m = np.tan(phi) / FS
        Z[:] = 0.0  # reset Z for each call
        for i in range(n):
            ca, sa = np.cos(alpha[i]), np.sin(alpha[i])
            cb, sb = np.cos(beta[i]), np.sin(beta[i])

            # Matrix A coefficients from equations (6) and (7)
            A = np.array([
                [tan_phi_m[i]*ca - sa,   -np.cos(theta[i+1])],
                [tan_phi_m[i]*sa + ca,   -np.sin(theta[i+1])]
            ])

            # Vector b from equations (6) and (7)
            # Pile: subtract H·cos(θp) from b0 (horizontal), subtract H·sin(θp) from b1 (vertical)
            b0 = (
                -c_m[i]*dl[i]*ca
                - P[i]*ca
                + u[i]*dl[i]*sa
                - Z[i]*np.cos(theta[i])
                - D[i]*sb
                + kw[i]
                + T[i]
                - H_pile[i]*np.cos(theta_p[i])
            )
            b1 = (
                -c_m[i]*dl[i]*sa
                - P[i]*sa
                - u[i]*dl[i]*ca
                + w[i]
                - Z[i]*np.sin(theta[i])
                + D[i]*cb
                - H_pile[i]*np.sin(theta_p[i])
            )
            
            N_i, Z_ip1 = np.linalg.solve(A, np.array([b0, b1]))
            Z[i+1] = Z_ip1
            N[i] = N_i  # store normal force on slice base
        return Z[n]

    if debug:
        r0 = residual(fs_guess)
        print(f"FS_guess={fs_guess:.6f} → residual={r0:.4g}")

    # Root-find FS such that the right-end interslice force Z[n] = 0, by Newton-
    # secant from fs_guess. (A bracketed brentq solver was tried to also catch the
    # genuinely-low-FS surfaces where the single-guess secant diverges, but it
    # resurrects non-physical over-strength roots near FS->0 on surfaces where the
    # secant correctly fails, which the admissibility check below does not reliably
    # separate from legitimate low-FS solutions. Deferred — see SEARCH/F6 item in
    # plans/plan_comprehensive_audit.md.)
    try:
        FS_opt = newton(residual, fs_guess, tol=tol, maxiter=max_iter)
    except Exception as e:
        return False, f"force_equilibrium failed to converge: {e}"

    # A non-positive factor of safety is unphysical: the secant has converged onto a
    # spurious root (it happens on some steep right-facing surfaces). Reject it.
    if not np.isfinite(FS_opt) or FS_opt <= 0:
        return False, f"force_equilibrium: non-physical factor of safety ({FS_opt})"

    # Re-evaluate at the root so N and Z reflect FS_opt (newton's last call may differ).
    residual(FS_opt)

    # Admissibility guard. A free search can drive the force-equilibrium solver onto
    # grossly non-physical surfaces (a large fraction of slices in base tension, or
    # pervasive interslice tension) and report a spurious low FS. Reject those by
    # EXTENT only — a few negative base normals or a non-monotonic thrust line occur
    # in valid solutions and are NOT rejected (valid benchmark criticals run ~0-4%
    # negative normals, <=20% interslice tension). See the SEARCH/F6 item in
    # plans/plan_comprehensive_audit.md.
    frac_N_neg = float(np.mean(N < 0)) if n else 0.0
    # Interslice "tension" sign is convention-dependent: on right-facing slopes the
    # caller negates theta_list, which flips the sign of Z, so there tension is Z>0.
    Z_int = Z[1:n]
    z_tension = (Z_int > 0) if right_facing else (Z_int < 0)
    frac_Z_tension = float(np.mean(z_tension)) if n > 1 else 0.0
    if frac_N_neg > 0.5 or frac_Z_tension > 0.5:
        return False, (
            "force_equilibrium: inadmissible solution "
            f"({100*frac_N_neg:.0f}% of base normals in tension, "
            f"{100*frac_Z_tension:.0f}% interslice tension)")

    slice_df['n_eff'] = N  # store effective normal forces in slice_df
    slice_df['z'] = Z[:-1]  # store interslice forces in slice_df, adjust length to n slices

    if debug:
        r_opt = residual(FS_opt)
        print(f" Converged FS = {FS_opt:.6f}, residual = {r_opt:.4g}")

    return True, {'FS': FS_opt}

def corps_engineers(slice_df, variant=2, debug=False):
    """
    Corps of Engineers (Modified Swedish) force-equilibrium solver.

    The interslice forces are assumed parallel at an inclination set by one of the
    two standard Corps conventions (USACE EM 1110-2-1902). Both are offered by
    commercial codes (e.g. SLOPE/W "Corps of Engineers #1/#2"):

      variant=1: a single constant θ = inclination of the line from the bottom to
                 the top of the failure surface (the crest-to-toe chord).
      variant=2 (default): θ at each slice boundary parallel to the GROUND SURFACE
                 (slice top).

    Variant 2 is the default because it is robust under an unconstrained
    search (variant 1 uses one fixed θ everywhere and can return a spuriously
    low FS on surfaces with steep segments, which a self-search exploits).
    Note that variant 2 is NOT the conservative choice: the ground-parallel
    convention systematically gives the HIGHEST (least conservative) factor of
    safety among XSLOPE's methods — typically a few percent to ~15% above
    Spencer — and the USACE manual itself notes the "average embankment slope"
    assumption can be unconservative. Report a rigorous method (Spencer) for
    design; use Corps for comparison.

    Both Corps conventions and Lowe & Karafiath set the inter-slice inclination
    from the signed ground/base slope and negate it on right-facing slopes, so
    the force-equilibrium factor of safety is mirror-symmetric.

    Parameters:
        slice_df (pd.DataFrame): Must include ['x_l','y_lt','x_r','y_rt'] plus all
                           columns required by force_equilibrium.
        variant (int): 1 = crest-to-toe chord, 2 = per-slice ground-parallel.

    Returns:
        Tuple(bool, dict or str): Whatever force_equilibrium returns.
    """
    n = len(slice_df)

    if variant == 1:
        # endpoints of the slip surface -> one constant θ
        x0, y0 = slice_df['x_l'].iat[0], slice_df['y_lt'].iat[0]
        x1, y1 = slice_df['x_r'].iat[-1], slice_df['y_rt'].iat[-1]
        dx = x1 - x0
        dy = y1 - y0
        if abs(dx) < 1e-12:
            theta_deg = 90.0
        else:
            theta_deg = np.degrees(np.arctan2(dy, dx))  # signed crest-to-toe chord
        theta_list = np.full(n+1, theta_deg)
        theta_out = theta_deg
    else:
        # per-boundary θ parallel to the ground surface (slice top)
        x_l = slice_df['x_l'].values
        y_lt = slice_df['y_lt'].values
        x_r = slice_df['x_r'].values
        y_rt = slice_df['y_rt'].values
        slope_top = (y_rt - y_lt) / (x_r - x_l)   # signed ground-surface slope per slice
        theta_list = np.zeros(n+1)
        for j in range(n+1):
            if j == 0:
                s = slope_top[0]
            elif j == n:
                s = slope_top[-1]
            else:
                s = 0.5*(slope_top[j-1] + slope_top[j])
            theta_list[j] = np.degrees(np.arctan(s))
        theta_out = float(np.mean(theta_list))

    # The force-equilibrium engine marches slices left→right with a hard-coded
    # sliding sense (Z[0]=0 at the left). On a right-facing slope the signed
    # ground/chord slope carries the wrong sign for that march, so the
    # inter-slice inclinations are negated — the same convention lowe_karafiath
    # uses. This is a no-op for left-facing slopes (the common case and every
    # shipped benchmark). Without it, force-equilibrium FS is asymmetric under
    # mirroring and can violate the expected method ordering on right-facing
    # geometries.
    right_facing = slice_df['y_lt'].iat[0] > slice_df['y_rt'].iat[-1]
    if right_facing:
        theta_list = -theta_list
        theta_out = -theta_out

    slice_df['theta'] = theta_list[:-1]  # store theta in slice_df. Adjust length to n slices.

    # delegate to your force_equilibrium solver
    success, results = force_equilibrium(slice_df, theta_list, debug=debug, right_facing=right_facing)
    if not success:
        return success, results
    else:
        results['method'] = 'corps_engineers'  # append method
        results['theta'] = theta_out           # append theta (mean for variant 2)
        return success, results

def lowe_karafiath(slice_df, debug=False):
    """
    Lowe-Karafiath limit equilibrium: variable interslice inclinations equal to
    the average of the top‐and bottom‐surface slopes of the two adjacent slices
    at each boundary.
    """
    n = len(slice_df)

    # grab boundary coords
    x_l = slice_df['x_l'].values
    y_lt = slice_df['y_lt'].values
    y_lb = slice_df['y_lb'].values
    x_r = slice_df['x_r'].values
    y_rt = slice_df['y_rt'].values
    y_rb = slice_df['y_rb'].values

    # determine facing
    right_facing = (y_lt[0] > y_rt[-1])

    # precompute each slice's top & bottom slopes
    widths   = (x_r - x_l)
    slope_top    = (y_rt - y_lt) / widths
    slope_bottom = (y_rb - y_lb) / widths

    # build θ_list for j=0..n
    if debug:
        print("boundary slopes (top/bottom) avg, θ_list:")  # header for debug list

    theta_list = np.zeros(n+1)
    for j in range(n+1):
        if j == 0:
            st = slope_top[0]
            sb = slope_bottom[0]
        elif j == n:
            st = slope_top[-1]
            sb = slope_bottom[-1]
        else:
            st = 0.5*(slope_top[j-1] + slope_top[j])
            sb = 0.5*(slope_bottom[j-1] + slope_bottom[j])

        avg_slope = 0.5*(st + sb)
        theta = np.degrees(np.arctan(avg_slope))

        # sign convention
        if right_facing:
            theta_list[j] =  -theta
        else:
            theta_list[j] = theta

        if debug:
            print(f"  j={j:2d}: st={st:.3f}, sb={sb:.3f}, θ={theta:.3f}°")

    slice_df['theta'] = theta_list[:-1]  # store theta in slice_df. Adjust length to n slices.

    # call your force_equilibrium solver
    success, results = force_equilibrium(slice_df, theta_list, debug=debug, right_facing=right_facing)
    if not success:
        return success, results
    else:
        results['method'] = 'lowe_karafiath'  # append method
        return success, results

def spencer(slice_df, tol=1e-4, max_iter = 100, debug_level=0):
    """
    Spencer's Method using Steve G. Wright's formulation from the UTEXAS v2  user manual.
    

    Parameters:
        slice_df (pd.DataFrame): must contain columns
            'alpha' (slice base inclination, degrees),
            'phi'   (slice friction angle, degrees),
            'c'     (cohesion),
            'dl'    (slice base length),
            'w'     (slice weight),
            'u'     (pore force per unit length),
            'dload' (distributed load),
            'beta'  (distributed load inclination, degrees),
            'kw'    (seismic force),
            't'     (tension crack water force),
            'p'     (reinforcement force),
            'h_pile' (pile force, optional), 'theta_p' (pile inclination, RADIANS, optional)

    Returns:
        float: FS where FS_force = FS_moment
        float: beta (degrees)
        bool: converged flag
    """

    alpha = np.radians(slice_df['alpha'].values)  # slice base inclination, degrees
    phi   = np.radians(slice_df['phi'].values)  # slice friction angle, degrees  
    c     = slice_df['c'].values  # cohesion
    dx    = slice_df['dx'].values  # slice width
    dl    = slice_df['dl'].values  # slice base length
    W     = slice_df['w'].values  # slice weight
    u     = slice_df['u'].values  # pore presssure
    x_c   = slice_df['x_c'].values  # center of base x-coordinate
    y_cb  = slice_df['y_cb'].values  # center of base y-coordinate
    y_lb   = slice_df['y_lb'].values  # left side base y-coordinate
    y_rb   = slice_df['y_rb'].values  # right side base y-coordinate
    P     = slice_df['dload'].values  # distributed load resultant 
    beta  = np.radians(slice_df['beta'].values)  # distributed load inclination, degrees
    kw    = slice_df['kw'].values  # seismic force
    V     = slice_df['t'].values  # tension crack water force
    y_v   = slice_df['y_t'].values  # tension crack water force y-coordinate
    R     = slice_df['p'].values  # reinforcement force

    # Pile forces (backward compatible — zeros if no pile data)
    n_slices = len(slice_df)
    H_pile    = slice_df['h_pile'].values  if 'h_pile'  in slice_df.columns else np.zeros(n_slices)
    theta_pile = slice_df['theta_p'].values if 'theta_p' in slice_df.columns else np.zeros(n_slices)
    x_pile    = slice_df['x_pile'].values  if 'x_pile'  in slice_df.columns else np.zeros(n_slices)
    y_pile    = slice_df['y_pile'].values  if 'y_pile'  in slice_df.columns else np.zeros(n_slices)

    # For now, we assume that reinforcement is flexible and therefore is parallel to the failure surface
    # at the bottom of the slice. Therefore, the psi value used in the derivation is set to alpha, 
    # and the point of action is the center of the base of the slice.
    psi = alpha  # psi is the angle of the reinforcement force from the horizontal
    y_r = y_cb  # y_r is the y-coordinate of the point of action of the reinforcement
    x_r = x_c  # x_r is the x-coordinate of the point of action of the reinforcement

    # use variable names to match the derivation.
    x_p = slice_df['d_x'].values  # distributed load x-coordinate
    y_p = slice_df['d_y'].values  # distributed load y-coordinate
    y_k = slice_df['y_cg'].values  # seismic force y-coordinate
    x_b = x_c  # center of base x-coordinate
    y_b = y_cb  # center of base y-coordinate


    tan_p = np.tan(phi)  # tan(phi)

    y_ct = slice_df['y_ct'].values
    right_facing = (y_ct[0] > y_ct[-1])
    # If right facing, swap angles and strengths. For most methods, you can use the normal angle conventions
    # and get the right answer. But for Spencer, due to the way that the moment equation is written,
    # you need to swap the angles and strengths if the slope is right facing.
    if right_facing:
        alpha = -alpha
        beta = -beta
        psi = -psi
        R = -R
        c = -c
        kw = -kw
        tan_p = -tan_p

    # pre-compute the trigonometric functions
    cos_a = np.cos(alpha)  # cos(alpha)
    sin_a = np.sin(alpha)  # sin(alpha)
    #tan_p = np.tan(phi)  # tan(phi)    # moved above
    cos_b = np.cos(beta)  # cos(beta)
    sin_b = np.sin(beta)  # sin(beta)
    sin_psi = np.sin(psi)  # sin(psi)
    cos_psi = np.cos(psi)  # cos(psi)

    # Pile force components. On right-facing slopes only the HORIZONTAL component
    # flips with the facing (it reverses with the sliding direction); the vertical
    # (upward) component is facing-independent — up is up. Flipping the whole H_pile
    # (as the other reflected quantities above) would wrongly negate the vertical
    # component and break battered (theta_p != 0) piles.
    H_cos_tp = H_pile * np.cos(theta_pile)  # horizontal component
    H_sin_tp = H_pile * np.sin(theta_pile)  # vertical component (upward)
    if right_facing:
        H_cos_tp = -H_cos_tp

    Fh = - kw - V + P * sin_b + R * cos_psi + H_cos_tp       # Equation (1) + pile horizontal
    Fv = - W - P * cos_b + R * sin_psi + H_sin_tp        # Equation (2) + pile vertical (upward)
    Mo = - P * sin_b * (y_p - y_b) - P * cos_b * (x_p - x_b) \
        + kw * (y_k - y_b) + V * (y_v - y_b) - R * cos_psi * (y_r - y_b) + R * sin_psi * (x_r - x_b) \
        - H_cos_tp * (y_pile - y_b) + H_sin_tp * (x_pile - x_b)  # Equation (3) + pile moment
    
    # ========== BEGIN SOLUTION ==========

    def compute_Q_and_yQ(F, theta_rad):
        """Compute Q and y_Q for given F and theta values."""
        # Equation (24): m_alpha
        ma = 1 / (np.cos(alpha - theta_rad) + np.sin(alpha - theta_rad) * tan_p / F)

        # Equation (23): Q
        Q = (- Fv * sin_a - Fh * cos_a - (c / F) * dl + (Fv * cos_a - Fh * sin_a + u * dl) * tan_p / F) * ma

        # Equation (26): y_Q with numerical safeguard
        # Add small epsilon to prevent divide-by-zero when Q * cos(theta) is very small
        Q_cos_theta = Q * np.cos(theta_rad)
        eps = 1e-10
        safe_denom = np.where(np.abs(Q_cos_theta) < eps, eps * np.sign(Q_cos_theta + eps), Q_cos_theta)
        y_q = y_b + Mo / safe_denom

        return Q, y_q
    
    def compute_residuals(F, theta_rad):
        """Compute residuals R1 and R2 for given F and theta values."""
        Q, y_q = compute_Q_and_yQ(F, theta_rad)

        # Equation (27): R1 = sum(Q)
        R1 = np.sum(Q)

        # Equation (28): R2 = sum(Q * (x_b * sin(theta) - y_Q * cos(theta)))
        R2 = np.sum(Q * (x_b * np.sin(theta_rad) - y_q * np.cos(theta_rad)))

        return R1, R2, Q, y_q


    def compute_derivatives(F, theta_rad, Q, y_q):

        """Compute all derivatives needed for Newton's method."""
        # Precompute trigonometric terms
        cos_alpha_theta = np.cos(alpha - theta_rad)
        sin_alpha_theta = np.sin(alpha - theta_rad)
        cos_theta = np.cos(theta_rad)
        sin_theta = np.sin(theta_rad)
        
        # Constants for Q expression (Equations 45-49)
        C1 = -Fv * sin_a - Fh * cos_a
        C2 = -c * dl + (Fv * cos_a - Fh * sin_a + u * dl) * tan_p
        C3 = cos_alpha_theta
        C4 = sin_alpha_theta * tan_p
        
        # Denominator for Q
        denom_Q = C3 + C4 / F
        
        # First-order partial derivatives of Q (Equations 50-51)
        dQ_dF = (-1 / denom_Q**2) * ((denom_Q * C2 / F**2) - (C1 + C2 / F) * C4 / F**2)
        
        dC3_dtheta = sin_alpha_theta  # Equation (55)
        dC4_dtheta = -cos_alpha_theta * tan_p  # Equation (56)
        dQ_dtheta = (-1 / denom_Q**2) * (C1 + C2 / F) * (dC3_dtheta + dC4_dtheta / F)

        # Partial derivatives of y_Q (Equations 59-60)
        # Add numerical safeguard to prevent divide-by-zero when Q * cos_theta is very small
        Q_cos_theta = Q * cos_theta
        eps = 1e-10  # Small epsilon to prevent exact zero division

        # Use np.where to handle element-wise operations safely
        # Where |Q * cos_theta| < eps, set derivatives to a large value to signal ill-conditioning
        safe_denom = np.where(np.abs(Q_cos_theta) < eps, eps * np.sign(Q_cos_theta + eps), Q_cos_theta)

        dyQ_dF = (-1 / safe_denom**2) * Mo * dQ_dF * cos_theta
        dyQ_dtheta = (-1 / safe_denom**2) * Mo * (dQ_dtheta * cos_theta - Q * sin_theta)
        
        # First-order partial derivatives of R1 (Equations 35-36)
        dR1_dF = np.sum(dQ_dF)
        dR1_dtheta = np.sum(dQ_dtheta)

        # First-order partial derivatives of R2 (Equations 40-41)
        dR2_dF = np.sum(dQ_dF * (x_b * sin_theta - y_q * cos_theta)) - np.sum(Q * dyQ_dF * cos_theta)
        dR2_dtheta = np.sum(dQ_dtheta * (x_b * sin_theta - y_q * cos_theta)) + np.sum(Q * (x_b * cos_theta + y_q * sin_theta - dyQ_dtheta * cos_theta))

        return dR1_dF, dR1_dtheta, dR2_dF, dR2_dtheta

    
    # Compute safe theta bounds to avoid m_alpha singularity.
    # m_alpha = 1/(cos(α-θ) + sin(α-θ)·tan(φ)/F)
    # Rewrite denominator as √(1+k²)·cos(α-θ-arctan(k)), k=tan_p/F.
    # Positive when |α - θ - arctan(k)| < π/2, giving per-slice θ bounds.
    def safe_theta_bounds(F_val, margin_deg=5):
        margin = np.radians(margin_deg)
        k = tan_p / F_val
        offset = alpha - np.arctan(k)  # per-slice offset
        lo = np.max(offset) - np.pi/2 + margin
        hi = np.min(offset) + np.pi/2 - margin
        if lo >= hi:
            lo, hi = -np.pi/2, np.pi/2  # fallback if no valid range
        return max(lo, -np.pi/2), min(hi, np.pi/2)

    def newton_solve(F0, theta0_rad, max_iter, tol, debug_level):
        """Run Newton iteration with backtracking line search. Returns (converged, F, theta_rad, iteration)."""
        theta_lo, theta_hi = safe_theta_bounds(F0)
        if theta0_rad < theta_lo or theta0_rad > theta_hi:
            theta0_rad = (theta_lo + theta_hi) / 2

        F = F0
        theta_rad = theta0_rad
        stagnation_count = 0

        for iteration in range(max_iter):
            R1, R2, Q, y_q = compute_residuals(F, theta_rad)
            residual_norm = R1**2 + R2**2

            if debug_level >= 1:
                if iteration == 0:
                    print(f"Iteration {1} - Initial: F = {F:.3f}, theta = {np.degrees(theta_rad):.3f}°, R1 = {R1:.6e}, R2 = {R2:.6e}")
                else:
                    print(f"Iteration {iteration + 1} - Updated: F = {F:.3f}, theta = {np.degrees(theta_rad):.3f}°, R1 = {R1:.6e}, R2 = {R2:.6e}")

            if abs(R1) < tol and abs(R2) < tol:
                if debug_level >= 1:
                    print(f"Converged in {iteration + 1} iterations, R1 = {R1:.6e}, R2 = {R2:.6e}")
                return True, F, theta_rad, iteration

            dR1_dF, dR1_dtheta, dR2_dF, dR2_dtheta = compute_derivatives(F, theta_rad, Q, y_q)

            J = np.array([[dR1_dF, dR1_dtheta],
                          [dR2_dF, dR2_dtheta]])

            try:
                cond_num = np.linalg.cond(J)
                if cond_num > 1e12:
                    return False, F, theta_rad, iteration
            except:
                return False, F, theta_rad, iteration

            try:
                delta_solution = np.linalg.solve(J, np.array([-R1, -R2]))
                delta_F = delta_solution[0]
                delta_theta = delta_solution[1]
            except np.linalg.LinAlgError:
                return False, F, theta_rad, iteration

            if debug_level >= 1:
                print(f"          Newton: delta_F = {delta_F:.3f}, delta_theta = {np.degrees(delta_theta):.3f}°, {delta_theta: .3f} (rad)")

            # Backtracking line search: find step size that reduces residual norm
            step = 1.0
            for _ in range(20):
                F_trial = F + step * delta_F
                theta_trial = theta_rad + step * delta_theta

                if F_trial <= 0:
                    step *= 0.5
                    continue

                theta_lo, theta_hi = safe_theta_bounds(F_trial)
                theta_trial = np.clip(theta_trial, theta_lo, theta_hi)

                R1_trial, R2_trial, _, _ = compute_residuals(F_trial, theta_trial)
                trial_norm = R1_trial**2 + R2_trial**2

                if trial_norm < residual_norm:
                    break
                step *= 0.5

            if debug_level >= 1 and step < 1.0:
                print(f"          Line search: step = {step:.4f}")

            # Detect stagnation (line search found no improvement)
            if step < 1e-6:
                stagnation_count += 1
                if stagnation_count >= 5:
                    if debug_level >= 1:
                        print(f"          Stagnation detected, exiting Newton solver early")
                    return False, F, theta_rad, iteration
            else:
                stagnation_count = 0

            F = F + step * delta_F
            theta_rad = theta_rad + step * delta_theta

            if F <= 0:
                F = 0.1

            theta_lo, theta_hi = safe_theta_bounds(F)
            theta_rad = np.clip(theta_rad, theta_lo, theta_hi)

        return False, F, theta_rad, max_iter - 1

    def max_tension_ratio(F_val, theta_val):
        """Compute max base tension as a multiple of cohesive capacity (c*dl) per slice.
        Tension beyond what cohesion can sustain is physically anomalous."""
        _, _, Q_val, _ = compute_residuals(F_val, theta_val)
        N_eff_val = -Fv * cos_a + Fh * sin_a + Q_val * np.sin(alpha - theta_val) - u * dl
        tension_mask = N_eff_val < 0
        if not np.any(tension_mask):
            return 0.0
        cohesive_capacity = np.abs(c) * dl  # use abs(c) in case of right_facing sign flip
        safe_capacity = np.where(cohesive_capacity > 0, cohesive_capacity, 1.0)
        ratios = np.where(
            tension_mask & (cohesive_capacity > 0),
            np.abs(N_eff_val) / safe_capacity,
            0.0
        )
        return np.max(ratios)

    MAX_TENSION_RATIO = 2.0  # reject if tension exceeds 2x cohesive capacity on any slice
    best_candidate = None  # (F, theta_rad, tension_ratio) — fallback if no clean solution

    def accept_or_save(F_val, theta_val, label=""):
        """Accept solution if tension magnitude is reasonable; otherwise save as fallback."""
        nonlocal best_candidate
        tr = max_tension_ratio(F_val, theta_val)
        if tr <= MAX_TENSION_RATIO:
            return True
        if best_candidate is None or tr < best_candidate[2]:
            best_candidate = (F_val, theta_val, tr)
        if debug_level >= 1:
            print(f"  {label}F={F_val:.3f}, θ={np.degrees(theta_val):.1f}° → "
                  f"max tension {tr:.1f}x cohesive capacity — continuing search")
        return False

    # Strategy 1: Newton with default initial guess
    F0 = 1.5
    theta0_rad = np.radians(-8.0) if right_facing else np.radians(8)

    converged, F, theta_rad, iteration = newton_solve(F0, theta0_rad, max_iter, tol, debug_level)
    if converged and not accept_or_save(F, theta_rad, "Newton default: "):
        converged = False

    # Strategy 2: Newton with Bishop's FS and multiple theta starts
    if not converged:
        try:
            success_bishop, result_bishop = bishop(slice_df)
            if success_bishop:
                F0_bishop = result_bishop['FS']
                if debug_level >= 1:
                    print(f"Retrying with Bishop FS = {F0_bishop:.3f} as initial guess")
                theta_tries = [theta0_rad, -theta0_rad, 0.0, np.radians(15), np.radians(-15)]
                for theta_try in theta_tries:
                    conv_b, F_b, theta_b, _ = newton_solve(F0_bishop, theta_try, max_iter, tol, debug_level)
                    if conv_b:
                        if accept_or_save(F_b, theta_b, "Newton Bishop: "):
                            converged, F, theta_rad = True, F_b, theta_b
                            break
        except Exception:
            pass  # Bishop may fail for non-circular surfaces

    # Strategy 3: scipy.optimize.root with multiple starts (ordered by |theta|)
    if not converged:
        from scipy.optimize import root as scipy_root

        def residual_vec(x):
            if x[0] <= 0.01:
                return np.array([1e10, 1e10])
            R1_v, R2_v, _, _ = compute_residuals(x[0], x[1])
            if not (np.isfinite(R1_v) and np.isfinite(R2_v)):
                return np.array([1e10, 1e10])
            return np.array([R1_v, R2_v])

        def jac_mat(x):
            if x[0] <= 0.01:
                return np.eye(2) * 1e10
            R1_v, R2_v, Q_v, yq_v = compute_residuals(x[0], x[1])
            if not (np.isfinite(R1_v) and np.isfinite(R2_v)):
                return np.eye(2) * 1e10
            d = compute_derivatives(x[0], x[1], Q_v, yq_v)
            return np.array([[d[0], d[1]], [d[2], d[3]]])

        F_starts = [1.0, 1.5, 2.0, 0.5, 3.0]
        try:
            sb, rb = bishop(slice_df)
            if sb:
                F_starts.insert(0, rb['FS'])
        except Exception:
            pass

        # Order by |theta| to prefer physically reasonable solutions first
        theta_starts = np.radians(np.array([0, 5, -5, 10, -10, 15, -15, 20, -20, 25, -25], dtype=float))

        for F_try in F_starts:
            for theta_try in theta_starts:
                try:
                    sol = scipy_root(residual_vec, [F_try, float(theta_try)],
                                     jac=jac_mat, method='hybr')
                    if sol.success and sol.x[0] > 0.01:
                        R1_c, R2_c, _, _ = compute_residuals(sol.x[0], sol.x[1])
                        if abs(R1_c) < tol and abs(R2_c) < tol:
                            if accept_or_save(sol.x[0], sol.x[1], "scipy: "):
                                converged = True
                                F = sol.x[0]
                                theta_rad = sol.x[1]
                                if debug_level >= 1:
                                    print(f"Converged via scipy root: F={F:.3f}, theta={np.degrees(theta_rad):.3f}°")
                                break
                except Exception:
                    continue
            if converged:
                break

    if not converged:
        if best_candidate is not None:
            return False, (f"Spencer's method: only solutions with anomalous base tension found "
                           f"({best_candidate[2]:.1f}x cohesive capacity, "
                           f"θ={np.degrees(best_candidate[1]):.1f}°)")
        return False, "Spencer's method did not converge within the maximum number of iterations."

    # Final computation of Q and y_q
    R1, R2, Q, y_q = compute_residuals(F, theta_rad)

    if debug_level >= 2:
        ma = 1 / (np.cos(alpha - theta_rad) + np.sin(alpha - theta_rad) * tan_p / F)
        slice_df['ma'] = ma
        slice_df['Q'] = Q
        slice_df['y_q'] = y_q
        slice_df['Fh'] = Fh
        slice_df['Fv'] = Fv
        slice_df['Mo'] = Mo
        # Print F and theta to 12 decimal places
        print(f"F = {F:.12f}, theta = {np.degrees(theta_rad):.12f}°")
        # Report the residuals
        print(f"R1 = {R1:.6e}, R2 = {R2:.6e}")
        # Debug print values per slice
        for i in range(len(Q)):
            print(f"Slice {i+1}: ma = {ma[i]:.3f}, Q = {Q[i]:.1f}, y_q = {y_q[i]:.2f}, Fh = {Fh[i]:.1f}, Fv = {Fv[i]:.1f}, Mo = {Mo[i]:.2f}")

    
    # Convert theta to degrees for output
    theta_opt = np.degrees(theta_rad)
    
    # ========== END SOLUTION ==========

    # Store theta in df
    slice_df['theta'] = theta_opt

    # --- Compute N_eff using Equation (18) ---
    N_eff = - Fv * cos_a + Fh * sin_a + Q * np.sin(alpha - theta_rad) - u * dl
    slice_df['n_eff'] = N_eff

    # --- Compute interslice forces Z using Equation (67) ---
    n = len(Q)
    Z = np.zeros(n+1)
    for i in range(n):
        Z[i+1] = Z[i] - Q[i] 
    slice_df['z'] = Z[:-1]        # Z_i acting on slice i's left face
 

    # --- Compute line of thrust using Equation (69) ---
    yt_l = np.zeros(n)  # the y-coordinate of the line of thrust on the left side of the slice.
    yt_r = np.zeros(n)  # the y-coordinate of the line of thrust on the right side of the slice.
    yt_l[0] = y_lb[0]  
    sin_theta = np.sin(theta_rad)
    cos_theta = np.cos(theta_rad)
    for i in range(n):
        if i == n - 1:
            yt_r[i] = y_rb[i]
        else:
            yt_r[i] = y_b[i] - ((Mo[i] - Z[i] * sin_theta * dx[i] / 2 - Z[i+1] * sin_theta * dx[i] / 2 - Z[i] * cos_theta * (yt_l[i] - y_b[i])) / (Z[i+1] * cos_theta))
            yt_l[i+1] = yt_r[i]
    slice_df['yt_l'] = yt_l
    slice_df['yt_r'] = yt_r
    
    # --- Return results ---
    results = {}
    results['method'] = 'spencer'
    results['FS'] = F
    results['theta'] = theta_opt


    return True, results






