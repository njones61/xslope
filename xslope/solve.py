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
from .hoekbrown import hb_tangent

def solve_selected(method_name, slice_df, rapid=False):
    """
    Executes a specified limit equilibrium solution method and displays results.

    Parameters
    ----------
    method_name : str
        Name of the solution method function to call. Must be one of:
        'oms', 'bishop', 'janbu', 'spencer', 'corps', 'lowe',
        'mprice'
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
        for w in result.get('warnings', []):
            print(f'  Spencer admissibility warning: {w}')
    elif func == janbu:
        print(f'Janbu Corrected FS={result["FS"]:.3f}, fo={result["fo"]:.2f}')
    elif func == corps:
        print(f'Corps Engineers: FS={result["FS"]:.3f}, theta={result["theta"]:.2f}')
        for w in result.get('warnings', []):
            print(f'  Corps admissibility warning: {w}')
    elif func == lowe:
        print(f'Lowe & Karafiath: FS={result["FS"]:.3f}')
        for w in result.get('warnings', []):
            print(f'  Lowe & Karafiath admissibility warning: {w}')
    elif func == mprice:
        print(f'Morgenstern-Price ({result["f_type"]}): FS={result["FS"]:.3f}, '
              f'lambda={result["lambda"]:.3f}')
        for w in result.get('warnings', []):
            print(f'  Morgenstern-Price admissibility warning: {w}')

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
    7. Morgenstern-Price Method

    If any method fails, an error message is displayed but execution continues
    with the remaining methods.
    """
    solve_selected('oms', slice_df, rapid=rapid)
    solve_selected('bishop', slice_df, rapid=rapid)
    solve_selected('janbu', slice_df, rapid=rapid)
    solve_selected('corps', slice_df, rapid=rapid)
    solve_selected('lowe', slice_df, rapid=rapid)
    solve_selected('spencer', slice_df, rapid=rapid)
    solve_selected('mprice', slice_df, rapid=rapid)


def _moment_arms(slice_df, Xo, Yo, alpha, right_facing):
    """Per-slice moment arms about the center of rotation (Xo, Yo), for the two
    moment methods (OMS, Bishop).

        xr   horizontal arm of the slice weight        ( = R·sinα  on a circle )
        a_S  moment arm of the base shear              ( = R       on a circle )
        a_N  perpendicular offset of the base normal   ( = 0       on a circle )

    On a true circular arc the base normal points straight at the center, so it
    exerts no moment, and every base is at radius R — which is what lets the
    textbook OMS/Bishop equations factor R out and drop the normal entirely.
    Neither holds on a COMPOSITE surface (a circle truncated at bedrock, see
    slice.CompositeSurface): the run along the floor is not at radius R and its
    normal misses the center. These are the general arms — SLOPE/W's per-slice R
    and its normal offset f — and they collapse identically to (R·sinα, R, 0) when
    the surface is a circle, so no circular result moves.

    The x-arms are taken in the orientation-normalized frame that alpha lives in
    (mirrored on right-facing slopes), the same correction as a_dx / pile_x_arm.
    """
    xr = slice_df['x_c'].values - Xo
    if right_facing:
        xr = -xr
    yr = slice_df['y_cb'].values - Yo
    sa, ca = np.sin(alpha), np.cos(alpha)
    a_S = xr * sa - yr * ca
    a_N = xr * ca + yr * sa
    return xr, a_S, a_N


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
    if 'r' not in slice_df.columns or slice_df['r'].iloc[0] is None:
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

    # v12 support/load terms (all zero when the columns are absent -> pre-v12
    # behavior). Axial reinforcement mirrors the pile-force pattern (a known
    # force at a point with a real direction); line loads mirror the
    # distributed-load pattern (an equivalent resultant at an inclination in the
    # dload convention); passive terms are factored by FS where supported.
    def _c12(name):
        return slice_df[name].values.astype(float) if name in slice_df.columns else np.zeros(n)
    P_pt  = _c12('p_pt')
    pa_cx = _c12('pa_cx'); pa_cy = _c12('pa_cy')
    pa_mx = _c12('pa_mx'); pa_my = _c12('pa_my')
    pp_cx = _c12('pp_cx'); pp_cy = _c12('pp_cy')
    pp_mx = _c12('pp_mx'); pp_my = _c12('pp_my')
    H_pas = _c12('h_pile_pas')
    LL    = _c12('lload')
    ll_b  = np.radians(_c12('ll_beta'))
    ll_x  = _c12('ll_x'); ll_y = _c12('ll_y')
    has_passive = bool(np.any(P_pt) or np.any(pp_cx) or np.any(pp_cy) or np.any(H_pas))
    # Axial components are stored normalized to cos(psi) > 0; on a right-facing
    # slope the tension on the sliding mass is the NEGATION of that direction,
    # so its vertical component flips sign. Only the vertical-equilibrium copies
    # flip: the horizontal component is consumed as a resisting magnitude, and
    # the moment terms flip their own x-part where used.
    _rf = slice_df['y_lb'].iat[0] > slice_df['y_rb'].iat[-1]
    pa_cy_v = -pa_cy if _rf else pa_cy
    pp_cy_v = -pp_cy if _rf else pp_cy

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
      + (H_pile - H_pas) * np.sin(alpha - theta_p)
      + (sin_alpha * pa_cx - cos_alpha * pa_cy_v)   # axial reinforcement, normal comp
      + LL * np.cos(alpha - ll_b)                 # line load (dload pattern)
    )

    #   Σ  Dᵢ·sinβᵢ·(Yo - d_{y,i}) 
    a_dy = Yo - d_y
    sum_Dy = np.sum(D * np.sin(beta) * a_dy)

    # ————————————————————————————————————————————————————————
    # 6) Build each piece of the DENOMINATOR exactly as Eqn 8, as MOMENTS about the
    # center of rotation. (The classic form divides the moment terms by R instead;
    # the two are the same equation, but only the moment form survives a composite
    # surface, where the base is not everywhere at radius R.)

    #  (B) = Σ  Dᵢ·cosβᵢ·(Xo - d_{x,i})
    # The vertical D-component (D·cosβ) acts at horizontal arm (d_x - Xo). That arm
    # flips sign when the slope faces right, and cosβ (unlike sinβ) is NOT corrected
    # by the β sign-flip applied in generate_slices, so flip it here. The horizontal
    # D-component term (a_dy / sum_Dy) is already sign-correct via that β-flip and is
    # left untouched; likewise the seismic (a_s) and tension-crack (a_t) arms.
    # An explicit facing from generate_slices (flat-arc / override case) wins;
    # otherwise fall back to the historical geometric test, byte-identical.
    right_facing = slice_df.attrs.get(
        'right_facing', slice_df['y_lb'].iat[0] > slice_df['y_rb'].iat[-1])

    # Moment arms about (Xo, Yo). a_S == R and a_N == 0 on a true circle.
    xr, a_S, a_N = _moment_arms(slice_df, Xo, Yo, alpha, right_facing)

    #  (A) = Σ [ Wᵢ · sinαᵢ ] · R  =  Σ [ Wᵢ · xrᵢ ]
    sum_W = np.sum(W * xr)

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
    h_act = H_pile - H_pas
    sum_pile_moment = np.sum(
        h_act * np.cos(theta_p) * (Yo - y_pile)
      + h_act * np.sin(theta_p) * pile_x_arm
    )
    sum_pile_moment_pas = np.sum(
        H_pas * np.cos(theta_p) * (Yo - y_pile)
      + H_pas * np.sin(theta_p) * pile_x_arm
    )

    # Axial reinforcement moment about the circle center (pile pattern: the
    # x-arm flips on right-facing slopes). M = cx*(Yo - y_r) + cy*(x_r - Xo),
    # linearized through the slice sums (mx = sum cy*x_r, my = sum cx*y_r).
    # (facing-invariant: on right-facing slopes the stored axial components and
    # the resisting-rotation sense both negate, so the raw expression stands)
    ra_x_part = (pa_mx - Xo * pa_cy)
    rp_x_part = (pp_mx - Xo * pp_cy)
    sum_reinf_moment_a = np.sum(pa_cx * Yo - pa_my + ra_x_part)
    sum_reinf_moment_p = np.sum(pp_cx * Yo - pp_my + rp_x_part)

    # Line-load moment terms mirror the distributed-load terms exactly.
    ll_x_arm = -(ll_x - Xo) if right_facing else (ll_x - Xo)
    sum_LLx = np.sum(LL * np.cos(ll_b) * ll_x_arm)
    sum_LLy = np.sum(LL * np.sin(ll_b) * (Yo - ll_y))

    # Put them together. Every term is now a moment about (Xo, Yo): the support/load
    # terms already were (that is why the classic form divides them by R), and the
    # slice terms pick up their arms here. Active P, sum_Dy, and the pile/axial
    # moments are known resisting forces, not factored by FS; line loads are known
    # applied loads. The base normal joins the driving side through a_N — zero for
    # every slice of a circle, non-zero along a composite surface's floor run.
    denominator = (sum_W + sum_Dx + sum_kw + sum_T + sum_LLx
                   - np.sum(P * a_S) - (sum_Dy + sum_LLy)
                   - sum_pile_moment - sum_reinf_moment_a
                   - np.sum((N_eff + u * dl) * a_N))

    # Passive support (ultimate capacities mobilized with the soil) joins the
    # RESISTING side: numerator gains the capacity moments. OMS is non-iterative,
    # so passive components cannot enter N_eff; their resisting moment is the
    # entire contribution (consistent with OMS's moment-only basis).
    numerator = (np.sum((c * dl + N_eff * tan_phi) * a_S) + np.sum(P_pt * a_S)
                 + sum_reinf_moment_p + sum_pile_moment_pas)

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

    if 'r' not in slice_df.columns or slice_df['r'].iloc[0] is None:
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

    # v12 support/load terms (all zero when the columns are absent)
    def _c12(name):
        return slice_df[name].values.astype(float) if name in slice_df.columns else np.zeros(n)
    P_pt  = _c12('p_pt')
    pa_cx = _c12('pa_cx'); pa_cy = _c12('pa_cy')
    pa_mx = _c12('pa_mx'); pa_my = _c12('pa_my')
    pp_cx = _c12('pp_cx'); pp_cy = _c12('pp_cy')
    pp_mx = _c12('pp_mx'); pp_my = _c12('pp_my')
    H_pas = _c12('h_pile_pas')
    LL    = _c12('lload')
    ll_b  = np.radians(_c12('ll_beta'))
    ll_x  = _c12('ll_x'); ll_y = _c12('ll_y')
    # Vertical-equilibrium copies of the axial components: flip on right-facing
    # slopes (see oms(); the moment terms below flip their own x-part instead).
    _rf = slice_df['y_lb'].iat[0] > slice_df['y_rb'].iat[-1]
    pa_cy_v = -pa_cy if _rf else pa_cy
    pp_cy_v = -pp_cy if _rf else pp_cy

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
    # Explicit facing (flat-arc / override) wins; else historical test, byte-identical.
    right_facing = slice_df.attrs.get(
        'right_facing', slice_df['y_lb'].iat[0] > slice_df['y_rb'].iat[-1])
    a_dx = d_x - Xo
    if right_facing:
        a_dx = -a_dx
    a_dy = Yo - d_y
    a_s  = Yo - y_cg
    a_t  = Yo - y_t

    # Slice moment arms about (Xo, Yo): a_S == R and a_N == 0 on a true circle, so
    # this reproduces the classic R-normalized Bishop equation exactly. A composite
    # surface (circle truncated at bedrock) needs the general arms. See _moment_arms.
    xr, a_S, a_N = _moment_arms(slice_df, Xo, Yo, alpha, right_facing)

    # Pile resisting moment about circle center (same term as OMS) — flip the
    # vertical-component arm (x_pile - Xo) on right-facing slopes, as for a_dx above.
    pile_x_arm = -(x_pile - Xo) if right_facing else (x_pile - Xo)
    h_act = H_pile - H_pas
    sum_pile_moment = np.sum(
        h_act * np.cos(theta_p) * (Yo - y_pile)
      + h_act * np.sin(theta_p) * pile_x_arm
    )
    sum_pile_moment_pas = np.sum(
        H_pas * np.cos(theta_p) * (Yo - y_pile)
      + H_pas * np.sin(theta_p) * pile_x_arm
    )

    # Axial reinforcement and line-load moment terms (see oms())
    # (facing-invariant: on right-facing slopes the stored axial components and
    # the resisting-rotation sense both negate, so the raw expression stands)
    ra_x_part = (pa_mx - Xo * pa_cy)
    rp_x_part = (pp_mx - Xo * pp_cy)
    sum_reinf_moment_a = np.sum(pa_cx * Yo - pa_my + ra_x_part)
    sum_reinf_moment_p = np.sum(pp_cx * Yo - pp_my + rp_x_part)
    ll_x_arm = -(ll_x - Xo) if right_facing else (ll_x - Xo)
    sum_LLx = np.sum(LL * np.cos(ll_b) * ll_x_arm)
    sum_LLy = np.sum(LL * np.sin(ll_b) * (Yo - ll_y))

    # Driving moment about (Xo, Yo) — active pile/axial moments are known resisting
    # forces, not factored by FS; line loads are known applied loads. Everything here
    # is FS-independent; the base-normal moment is not, and is added in the loop.
    sum_W = np.sum(W * xr)
    sum_Dx = np.sum(D * cos_beta * a_dx)
    sum_Dy = np.sum(D * sin_beta * a_dy)
    sum_kw = np.sum(kw * a_s)
    sum_T = np.sum(T * a_t)
    denom_fixed = (sum_W + sum_Dx + sum_kw + sum_T + sum_LLx
                   - np.sum(P * a_S) - (sum_Dy + sum_LLy)
                   - sum_pile_moment - sum_reinf_moment_a)

    # Passive capacity moments join the resisting side (the numerator of F)
    res_passive = np.sum(P_pt * a_S) + sum_reinf_moment_p + sum_pile_moment_pas

    # Vertical component of pile force (upward, reduces net downward force on slice)
    H_sin_tp = H_pile * np.sin(theta_p)

    # Iterative solution
    F = 1.0
    for _ in range(max_iter):
        # Compute N_eff — H·sin(θp) enters vertical equilibrium
        # Vertical equilibrium: known loads at full value; passive support
        # mobilizes with the soil, so its vertical components carry 1/F.
        vert_known = (W + D * cos_beta + LL * np.cos(ll_b)
                      - P * sin_alpha - pa_cy_v - (H_pile - H_pas) * np.sin(theta_p))
        vert_pas = (P_pt * sin_alpha + pp_cy_v + H_pas * np.sin(theta_p))
        num_N = (
            vert_known - vert_pas / F
            - u * dl * cos_alpha
            - (c * dl * sin_alpha) / F
        )
        denom_N = cos_alpha + (sin_alpha * tan_phi) / F
        N_eff = num_N / denom_N

        # The base normal exerts a moment only where it misses the center of
        # rotation (a_N = 0 on every slice of a true arc, non-zero along a composite
        # surface's floor run). It depends on N_eff, so it lives inside the loop.
        denominator = denom_fixed - np.sum((N_eff + u * dl) * a_N)
        if denominator <= 0:
            return False, "Net driving moment is non-positive (resisting forces exceed driving forces)."

        # Resisting moment: the mobilized base shear (c·dl + N'·tanφ) on its own arm
        numerator = np.sum((c * dl + N_eff * tan_phi) * a_S) + res_passive
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

    # v12 support/load terms (all zero when the columns are absent)
    def _c12(name):
        return slice_df[name].values.astype(float) if name in slice_df.columns else np.zeros(n)
    P_pt  = _c12('p_pt')
    pa_cx = _c12('pa_cx'); pa_cy = _c12('pa_cy')
    pp_cx = _c12('pp_cx'); pp_cy = _c12('pp_cy')
    H_pas = _c12('h_pile_pas')
    LL    = _c12('lload')
    ll_b  = np.radians(_c12('ll_beta'))
    # Vertical-equilibrium copies of the axial components: flip on right-facing
    # slopes (see oms()).
    _rf = slice_df['y_lb'].iat[0] > slice_df['y_rb'].iat[-1]
    pa_cy_v = -pa_cy if _rf else pa_cy
    pp_cy_v = -pp_cy if _rf else pp_cy
    H_sin_act = (H_pile - H_pas) * np.sin(theta_p)
    H_cos_act = (H_pile - H_pas) * np.cos(theta_p)
    H_sin_pas = H_pas * np.sin(theta_p)
    H_cos_pas = H_pas * np.cos(theta_p)

    # External horizontal driving forces (independent of F): seismic (driving),
    # distributed-load horizontal component (driving), tension-crack water
    # (driving), reinforcement and pile horizontal components (resisting).
    #
    # Reinforcement is applied parallel to the slice base (the flexible assumption,
    # psi = alpha), so its HORIZONTAL component is P*cos(alpha) -- matching the
    # P*sin(alpha) vertical component already used in num_N below, and matching
    # spencer's `R * cos_psi` with psi = alpha. Summing the bare magnitude here
    # over-credited the reinforcement and raised FS non-conservatively.
    horiz_ext = (np.sum(kw) + np.sum(D * sin_beta) + np.sum(T) + np.sum(LL * np.sin(ll_b))
                 - np.sum(P * cos_alpha) - np.sum(pa_cx) - np.sum(H_cos_act))

    # Iterate F: the base normal depends on F through m_alpha, exactly as Bishop.
    F = 1.0
    N_eff = None
    converged = False
    for _ in range(max_iter):
        m_alpha = cos_alpha + sin_alpha * tan_phi / F
        # Effective base normal from VERTICAL slice equilibrium (Bishop's N').
        num_N = (W + D * cos_beta + LL * np.cos(ll_b) - P * sin_alpha - pa_cy_v - H_sin_act
                 - (P_pt * sin_alpha + pp_cy_v + H_sin_pas) / F
                 - u * dl * cos_alpha - (c * dl * sin_alpha) / F)
        N_eff = num_N / m_alpha
        N_total = N_eff + u * dl
        # Horizontal force equilibrium: resisting shear (projected horizontally)
        # over driving (base-normal horizontal component + external horizontal).
        resisting = np.sum((c * dl + N_eff * tan_phi) * cos_alpha)
        # Passive support mobilizes with the soil: its horizontal components
        # reduce the driving side factored by 1/F.
        driving = (np.sum(N_total * sin_alpha) + horiz_ext
                   - (np.sum(P_pt * cos_alpha) + np.sum(pp_cx) + np.sum(H_cos_pas)) / F)
        if driving <= 0:
            return False, "Net horizontal driving force is non-positive (resisting forces exceed driving forces)."
        F_new = resisting / driving
        if abs(F_new - F) < tol:
            F = F_new
            converged = True
            break
        # Damped update: the passive-support 1/F terms on the driving side can
        # make the plain fixed-point map non-contractive (oscillating between
        # two values); the half-step converges to the same fixed point.
        F = 0.5 * (F + F_new)

    if not converged:
        return False, "Janbu method did not converge within the maximum number of iterations."

    FS_base = F

    # === Compute Janbu correction factor ===

    # Get failure surface endpoints
    x_l = slice_df['x_l'].iloc[0]   # leftmost x
    y_lt = slice_df['y_lb'].iloc[0]   # left surface endpoint (base corner; the slice
    x_r = slice_df['x_r'].iloc[-1]  # rightmost x        top sits at crest elevation
    y_rt = slice_df['y_rb'].iloc[-1]  # right surface endpoint  against a vertical face)

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


def _equilibrium_march(alpha, phi, c, w, u, dl, D, beta, kw, T, P, H_pile, theta_p, theta, FS,
                       P_pt=None, pa_cx=None, pa_cy=None, pp_cx=None, pp_cy=None,
                       H_pas=None, LL=None, ll_b=None):
    """Shared fixed-FS force-equilibrium march (left -> right over slices).

    The per-slice engine extracted from `force_equilibrium` so it can be reused:
    the force-equilibrium methods (Corps of Engineers, Lowe-Karafiath — via
    `force_equilibrium`) and, in future, Morgenstern-Price all march on it. It does
    NOT root-find FS — it evaluates ONE trial `FS` and returns the resulting base
    normals and interslice forces. Callers form their own residuals from the
    output: force closure is `Z[n]`; M-P additionally accumulates a moment residual
    from the returned `N`.

    Parameters (per-slice numpy arrays of length n; angles already in RADIANS):
        alpha, phi, c, w, u, dl  base inclination, friction, cohesion, weight,
                                 pore pressure/length, base length
        D, beta                  distributed surface load magnitude and inclination
        kw                       horizontal seismic force
        T                        tension-crack water force (nonzero only on the
                                 cracked slice)
        P                        reinforcement force (parallel to base)
        H_pile, theta_p          pile force magnitude and inclination (theta_p in
                                 RADIANS, like everything else here)
        theta                    interslice-force inclination at each BOUNDARY,
                                 length n+1 (radians)
        FS                       trial factor of safety
    Returns:
        (N, Z): N = base normal force per slice (length n, effective normal);
                Z = interslice resultant at each boundary (length n+1, Z[0] = 0).

    SIGN / RIGHT-FACING CONVENTION  ── do not forget ─────────────────────────────
    This march uses every per-slice quantity (alpha, beta, kw, T, P, ...) EXACTLY
    as passed; it performs NO internal sign flips for slope facing. Right-facing
    slopes are handled by the CALLER negating the whole `theta_list` before the
    call (see `corps` / `lowe`). Negating theta flips the sign
    of Z, which is why the tension guard back in `force_equilibrium` keys its
    Z<0 / Z>0 test on the `right_facing` flag.

    This DIFFERS from `spencer()`, which instead flips alpha, beta, psi, R, c, kw,
    V, tan_p and the pile horizontal component INTERNALLY (its `if right_facing:`
    block). Morgenstern-Price is built on THIS march, so its moment accumulator
    must follow the caller-negates-theta convention used here, NOT Spencer's
    internal-flip convention. The `f(x)=1 == Spencer` validation must therefore be
    run on a right-facing geometry specifically, to confirm the two reconcile.
    ──────────────────────────────────────────────────────────────────────────────
    """
    n = len(alpha)

    # Fully VECTORIZED march. Per slice the 2x2 gives, in closed form,
    #   Z[i+1] = (a11·b1 − a21·b0)/det   with   b0 = b0base − Z[i]·ct[i],
    #                                            b1 = b1base − Z[i]·st[i],
    # which is AFFINE in Z[i]:  Z[i+1] = p[i] + q[i]·Z[i]  (a LINEAR recurrence). A
    # linear recurrence has a closed form via cumulative products, so the whole march
    # is np.cumprod/np.cumsum — no Python loop (~2x faster than the per-slice loop,
    # and it speeds up Corps/L-K too). The multiplier q stays ~1 for physical slip
    # surfaces (cumprod magnitude ~1 even at 120 slices), so it is numerically stable;
    # we fall back to the scalar recurrence if cumprod ever goes non-finite.
    ca = np.cos(alpha); sa = np.sin(alpha)
    cb = np.cos(beta);  sb = np.sin(beta)
    ct = np.cos(theta); st = np.sin(theta)            # interslice angle, length n+1
    cp = np.cos(theta_p); sp = np.sin(theta_p)
    c_m = c / FS
    tan_phi_m = np.tan(phi) / FS

    # v12 support/load terms (zeros when not supplied). Axial reinforcement
    # mirrors the pile terms; passive support carries 1/FS (mobilizes with the
    # soil, like c_m); line loads mirror the distributed-load terms.
    z = np.zeros(n)
    P_pt  = z if P_pt  is None else P_pt
    pa_cx = z if pa_cx is None else pa_cx
    pa_cy = z if pa_cy is None else pa_cy
    pp_cx = z if pp_cx is None else pp_cx
    pp_cy = z if pp_cy is None else pp_cy
    H_pas = z if H_pas is None else H_pas
    LL    = z if LL    is None else LL
    ll_b  = z if ll_b  is None else ll_b
    h_act_cp = (H_pile - H_pas) * cp
    h_act_sp = (H_pile - H_pas) * sp
    a11 = tan_phi_m * ca - sa                          # A first column (N coefficients)
    a21 = tan_phi_m * sa + ca
    a12 = -ct[1:]; a22 = -st[1:]                       # A second column at boundary i+1
    det = a11 * a22 - a12 * a21
    b0base = (-c_m*dl*ca - P*ca - pa_cx + u*dl*sa - D*sb - LL*np.sin(ll_b) + kw + T
              - h_act_cp - (P_pt*ca + pp_cx + H_pas*cp) / FS)
    b1base = (-c_m*dl*sa - P*sa - pa_cy - u*dl*ca + w + D*cb + LL*np.cos(ll_b)
              - h_act_sp - (P_pt*sa + pp_cy + H_pas*sp) / FS)
    cti = ct[:-1]; sti = st[:-1]                       # interslice angle at boundary i

    # Z[i+1] = p[i] + q[i]·Z[i]
    p = (a11 * b1base - a21 * b0base) / det
    q = (a21 * cti - a11 * sti) / det
    Z = np.empty(n + 1)
    Z[0] = 0.0
    Pc = np.empty(n + 1)
    Pc[0] = 1.0
    np.cumprod(q, out=Pc[1:])                          # Pc[i] = prod_{j<i} q[j]
    if np.all(np.isfinite(Pc)) and np.all(Pc[1:] != 0.0):
        Z[1:] = Pc[1:] * np.cumsum(p / Pc[1:])
    else:                                              # pathological geometry → scalar scan
        for i in range(n):
            Z[i + 1] = p[i] + q[i] * Z[i]

    b0 = b0base - Z[:-1] * cti
    b1 = b1base - Z[:-1] * sti
    N = (b0 * a22 - a12 * b1) / det
    return N, Z


def _moment_about_origin(N, FS, *, x_c, y_cg, y_cb, alpha, c, phi, dl, u, w,
                         D, beta, d_x, d_y, kw, V, y_t, P, H_pile, theta_p,
                         x_pile, y_pile, P_pt=None, pa_mx=None, pa_my=None,
                         pp_mx=None, pp_my=None, H_pas=None, LL=None, ll_b=None,
                         ll_x=None, ll_y=None):
    """Moment of the WHOLE sliding mass about the origin (0,0), CCW positive.

    This is the Morgenstern-Price moment residual (plan §4a). Its root in `FS` at a
    fixed interslice-force distribution is `F_m`. The moment of a force with global
    components `(F_x, F_y)` acting at `(x, y)` is `M = x·F_y − y·F_x`.

    The interslice forces are INTERNAL: on the whole mass they cancel exactly (equal
    and opposite on each shared boundary), regardless of `λ` or the thrust line, so
    they do not appear here — that is what lets us sum about the origin without
    knowing the line of thrust.

    All angle arrays are in RADIANS. `N` is the EFFECTIVE base normal from the
    march; the total base reaction for its moment arm is `N + u·dl` (effective +
    pore uplift, plan §4a rows 2/2b), and the mobilized base shear is
    `S = (c·dl + N·tanφ)/FS`.

    Per-slice terms (summed), each straight from the §4a table:
        weight W        →  −W·x_c
        normal N′+U     →  (N + u·dl)·(x_c·cosα + y_cb·sinα)
        base shear S    →  S·(x_c·sinα − y_cb·cosα)
        surface load D  →  −D·(d_x·cosβ + d_y·sinβ)
        seismic kW      →  kW·y_cg
        tension water V →  V·y_t
        reinforcement R →  R·(x_c·sinα − y_cb·cosα)
        pile H          →  H·(x_pile·sinθ_p − y_pile·cosθ_p)

    SIGN/RIGHT-FACING: convention-agnostic — it uses whatever arrays it is given.
    Facing is handled by the caller (see `_mp_march`); this keeps the math one place.
    """
    sa, ca = np.sin(alpha), np.cos(alpha)
    sb, cb = np.sin(beta), np.cos(beta)
    S = (c * dl + N * np.tan(phi)) / FS          # mobilized base shear (N is effective)
    M = (
        -w * x_c                                  # weight (0,−W) at (x_c, y_cg)
        + (N + u * dl) * (x_c * ca + y_cb * sa)   # total normal N′+U at (x_c, y_cb)
        + S * (x_c * sa - y_cb * ca)              # base shear at (x_c, y_cb)
        - D * (d_x * cb + d_y * sb)               # surface load at (d_x, d_y)
        + kw * y_cg                               # seismic (horizontal) at y_cg
        + V * y_t                                 # tension-crack water (horiz) at y_t
        + P * (x_c * sa - y_cb * ca)              # reinforcement (∥ base) at (x_c, y_cb)
        + H_pile * (x_pile * np.sin(theta_p) - y_pile * np.cos(theta_p))  # pile
    )
    Msum = float(np.sum(M))
    # v12 terms (zeros when not supplied): axial reinforcement (components at the
    # crossing point, moment about origin = mx - my), passive support factored by
    # FS, and line loads (distributed-load pattern at their own point).
    z = np.zeros_like(x_c)
    P_pt  = z if P_pt  is None else P_pt
    pa_mx = z if pa_mx is None else pa_mx
    pa_my = z if pa_my is None else pa_my
    pp_mx = z if pp_mx is None else pp_mx
    pp_my = z if pp_my is None else pp_my
    H_pas = z if H_pas is None else H_pas
    LL    = z if LL    is None else LL
    ll_b  = z if ll_b  is None else ll_b
    ll_x  = z if ll_x  is None else ll_x
    ll_y  = z if ll_y  is None else ll_y
    Msum += float(np.sum(
        (pa_mx - pa_my)                                     # axial active
        + (pp_mx - pp_my) / FS                              # axial passive
        + P_pt * (x_c * sa - y_cb * ca) / FS                # tangent passive
        - H_pas * (x_pile * np.sin(theta_p) - y_pile * np.cos(theta_p)) * (1.0 - 1.0 / FS)
        - LL * (ll_x * np.cos(ll_b) + ll_y * np.sin(ll_b))  # line load
    ))
    return Msum


def _mp_extract(slice_df, right_facing):
    """Pull every array the M-P march needs out of the DataFrame ONCE, apply the
    right-facing flips, and return a plain dict of numpy arrays.

    Hoisting this out of the per-iteration residual avoids re-reading ~25 pandas
    columns on every Newton step — `DataFrame.__getitem__` otherwise dominates the
    M-P solve (it is read ~13× per solve, once per residual evaluation).

    RIGHT-FACING flips mirror spencer()'s internal set (alpha, beta, phi=tan_p, c, P,
    kw, V) plus the pile `theta_p -> π-θp` (flips cos θp, not sin θp) so that f(x)=1
    reduces to Spencer — validated by the f(x)=1 == Spencer gate (S3b).
    """
    n = len(slice_df)
    g = lambda col: slice_df[col].values.astype(float)
    A = dict(
        alpha=np.radians(slice_df['alpha'].values), phi=np.radians(slice_df['phi'].values),
        c=g('c'), w=g('w'), u=g('u'), dl=g('dl'), D=g('dload'),
        beta=np.radians(slice_df['beta'].values), kw=g('kw'), V=g('t'), P=g('p'),
        H_pile=g('h_pile') if 'h_pile' in slice_df.columns else np.zeros(n),
        theta_p=g('theta_p') if 'theta_p' in slice_df.columns else np.zeros(n),
        x_c=g('x_c'), y_cg=g('y_cg'), y_cb=g('y_cb'),
        d_x=g('d_x'), d_y=g('d_y'), y_t=g('y_t'),
        x_pile=g('x_pile') if 'x_pile' in slice_df.columns else np.zeros(n),
        y_pile=g('y_pile') if 'y_pile' in slice_df.columns else np.zeros(n),
        # v12 support/load terms
        P_pt=g('p_pt') if 'p_pt' in slice_df.columns else np.zeros(n),
        pa_cx=g('pa_cx') if 'pa_cx' in slice_df.columns else np.zeros(n),
        pa_cy=g('pa_cy') if 'pa_cy' in slice_df.columns else np.zeros(n),
        pa_mx=g('pa_mx') if 'pa_mx' in slice_df.columns else np.zeros(n),
        pa_my=g('pa_my') if 'pa_my' in slice_df.columns else np.zeros(n),
        pp_cx=g('pp_cx') if 'pp_cx' in slice_df.columns else np.zeros(n),
        pp_cy=g('pp_cy') if 'pp_cy' in slice_df.columns else np.zeros(n),
        pp_mx=g('pp_mx') if 'pp_mx' in slice_df.columns else np.zeros(n),
        pp_my=g('pp_my') if 'pp_my' in slice_df.columns else np.zeros(n),
        H_pas=g('h_pile_pas') if 'h_pile_pas' in slice_df.columns else np.zeros(n),
        LL=g('lload') if 'lload' in slice_df.columns else np.zeros(n),
        ll_b=np.radians(g('ll_beta')) if 'll_beta' in slice_df.columns else np.zeros(n),
        ll_x=g('ll_x') if 'll_x' in slice_df.columns else np.zeros(n),
        ll_y=g('ll_y') if 'll_y' in slice_df.columns else np.zeros(n),
    )
    if right_facing:
        for k in ('alpha', 'beta', 'phi', 'c', 'P', 'kw', 'V'):
            A[k] = -A[k]
        A['theta_p'] = np.pi - A['theta_p']
        # Axial reinforcement enters in real components at real positions and
        # the stored direction (cos(psi) > 0) is the NEGATION of the real
        # tension on a right-facing slope: negate the whole vector (all sums),
        # as in spencer().
        for k in ('pa_cx', 'pa_cy', 'pa_mx', 'pa_my',
                  'pp_cx', 'pp_cy', 'pp_mx', 'pp_my'):
            A[k] = -A[k]
        # Line loads mirror beta's flip
        A['ll_b'] = -A['ll_b']
        # tangent-passive mirrors P
        A['P_pt'] = -A['P_pt']
    return A


def _mp_residuals(A, f_vals, lam, FS):
    """`(N, Z, force_res, moment_res)` from PRE-EXTRACTED arrays `A` (`_mp_extract`),
    at interslice scale λ and trial FS. No DataFrame access — this is what the
    Approach-A/B inner loops call every iteration.

    `θ_j = arctan(λ·f_vals_j)` (M-P's `tan θ = λ·f`); `force_res = Z[n]` → 0 at force
    equilibrium; `moment_res` → 0 at moment equilibrium. With `f_vals ≡ 1` this is
    Spencer (constant θ) — the make-or-break regression.
    """
    theta = np.arctan(lam * f_vals)
    N, Z = _equilibrium_march(A['alpha'], A['phi'], A['c'], A['w'], A['u'], A['dl'],
                              A['D'], A['beta'], A['kw'], A['V'], A['P'],
                              A['H_pile'], A['theta_p'], theta, FS,
                              P_pt=A['P_pt'], pa_cx=A['pa_cx'], pa_cy=A['pa_cy'],
                              pp_cx=A['pp_cx'], pp_cy=A['pp_cy'], H_pas=A['H_pas'],
                              LL=A['LL'], ll_b=A['ll_b'])
    force_res = Z[-1]
    moment_res = _moment_about_origin(
        N, FS, x_c=A['x_c'], y_cg=A['y_cg'], y_cb=A['y_cb'], alpha=A['alpha'],
        c=A['c'], phi=A['phi'], dl=A['dl'], u=A['u'], w=A['w'], D=A['D'], beta=A['beta'],
        d_x=A['d_x'], d_y=A['d_y'], kw=A['kw'], V=A['V'], y_t=A['y_t'], P=A['P'],
        H_pile=A['H_pile'], theta_p=A['theta_p'], x_pile=A['x_pile'], y_pile=A['y_pile'],
        P_pt=A['P_pt'], pa_mx=A['pa_mx'], pa_my=A['pa_my'], pp_mx=A['pp_mx'],
        pp_my=A['pp_my'], H_pas=A['H_pas'], LL=A['LL'], ll_b=A['ll_b'],
        ll_x=A['ll_x'], ll_y=A['ll_y'])
    return N, Z, force_res, moment_res


def _mp_march(slice_df, lam, f_vals, FS, right_facing=False):
    """Morgenstern-Price march for a single `(λ, FS)`: extract arrays then evaluate.

    Thin wrapper over `_mp_extract` + `_mp_residuals`, kept for direct/standalone
    callers (tests, gates). The solver itself extracts ONCE and calls
    `_mp_residuals` in the loop. Returns `(N, Z, force_res, moment_res)`.
    """
    n = len(slice_df)
    f_vals = np.asarray(f_vals, dtype=float)
    if len(f_vals) != n + 1:
        raise ValueError(f"f_vals length ({len(f_vals)}) must be n+1 ({n+1})")
    return _mp_residuals(_mp_extract(slice_df, right_facing), f_vals, lam, FS)


def _admissibility_warnings(c, N_eff, Z, y_lt=None, y_lb=None, yt_l=None):
    """Report-only admissibility screen shared by spencer/mprice/corps/lowe.

    Returns a list of Duncan & Wright admissibility notes for an ALREADY-ACCEPTED
    solution — it never affects FS, convergence, or acceptance. The tension guard
    in each solver normalizes base tension by cohesive capacity, so a c=0 slice
    can carry unbounded base tension without tripping it (VP30: -71 kN/m on the
    cohesionless crack-face sliver scores 0.0 and passes). These report the
    signatures that guard cannot see. Deliberately narrow — every accepted
    solution already satisfies the guard on cohesive slices, so re-flagging their
    tension would be noise:
      - base tension only on COHESIONLESS slices (the guard's blind spot);
      - interslice tension only against a clearly compressive field (when every
        Z is negative the sign convention itself is ambiguous, e.g. right-facing
        nailed walls, and no verdict is offered);
      - thrust line outside the slice on >10% of interior boundaries, with
        NaN/inf ratios counting as outside — skipped when no thrust line is
        supplied (corps/lowe expose no line of thrust).

    Z must arrive in the physical convention (tension < 0); callers whose march
    stores the opposite sign negate it before passing it in.
    """
    warns = []
    c = np.asarray(c, dtype=float)
    N_eff = np.asarray(N_eff, dtype=float)
    Z = np.asarray(Z, dtype=float)
    n_scale = float(np.max(np.abs(N_eff))) if len(N_eff) else 0.0
    cohesionless = np.abs(c) <= 1e-9
    bad_n = np.flatnonzero(cohesionless & (N_eff < -0.01 * n_scale))
    if bad_n.size:
        worst = int(bad_n[np.argmin(N_eff[bad_n])])
        warns.append(f"base tension on {bad_n.size} cohesionless slice(s), "
                     f"worst N' = {float(N_eff[worst]):.1f} at slice {worst + 1}")
    Z_int = Z[1:-1]
    if Z_int.size:
        z_min = float(np.min(Z_int))
        z_max = float(np.max(Z_int))
        if z_min < 0 and z_max > 0 and -z_min > 0.10 * z_max:
            warns.append(f"interslice tension (min Z = {z_min:.1f} vs max "
                         f"compression {z_max:.1f})")
    if yt_l is not None and y_lt is not None and y_lb is not None:
        y_lt_arr = np.asarray(y_lt, dtype=float)
        y_lb_arr = np.asarray(y_lb, dtype=float)
        yt_l_arr = np.asarray(yt_l, dtype=float)
        h_bnd = y_lt_arr[1:] - y_lb_arr[1:]          # interior boundaries 1..n-1
        tall = h_bnd > 1e-9
        if np.any(tall):
            t_ratio = (yt_l_arr[1:][tall] - y_lb_arr[1:][tall]) / h_bnd[tall]
            inside = np.isfinite(t_ratio) & (t_ratio >= -0.05) & (t_ratio <= 1.05)
            frac_out = 1.0 - float(np.mean(inside))
            if frac_out > 0.10:
                warns.append(f"line of thrust outside the slice on {frac_out:.0%} "
                             f"of boundaries")
    return warns


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

    # v12 support/load terms (all zero when the columns are absent)
    def _c12(name):
        return slice_df[name].values.astype(float) if name in slice_df.columns else np.zeros(n)
    P_pt  = _c12('p_pt')
    pa_cx = _c12('pa_cx'); pa_cy = _c12('pa_cy')
    pp_cx = _c12('pp_cx'); pp_cy = _c12('pp_cy')
    H_pas = _c12('h_pile_pas')
    LL    = _c12('lload')
    ll_b  = np.radians(_c12('ll_beta'))
    # The march mixes the orientation-normalized frame (alpha, theta) with the
    # stored axial components (normalized to cos(psi) > 0). The horizontal
    # component enters as a resisting magnitude (facing-agnostic), but the real
    # tension on a right-facing slope is the negation of the stored direction,
    # so its VERTICAL component flips (same correction as janbu()).
    if right_facing:
        pa_cy = -pa_cy
        pp_cy = -pp_cy

    N = np.zeros(n)  # normal forces on slice bases (filled by each march)
    Z = np.zeros(n+1)  # interslice forces, Z[0] = 0 by definition (no force entering leftmost slice)

    def residual(FS):
        """Right-side interslice force Z[n] for a given FS (force-closure residual).

        Delegates the per-slice march to the shared `_equilibrium_march`, then
        copies its output into the persistent N/Z arrays (via slice assignment) so
        the admissibility check and the slice_df writes after the root-find see the
        values at the converged FS. Pure extraction — arithmetic is unchanged.
        """
        N_march, Z_march = _equilibrium_march(
            alpha, phi, c, w, u, dl, D, beta, kw, T, P, H_pile, theta_p, theta, FS,
            P_pt=P_pt, pa_cx=pa_cx, pa_cy=pa_cy, pp_cx=pp_cx, pp_cy=pp_cy,
            H_pas=H_pas, LL=LL, ll_b=ll_b)
        N[:] = N_march
        Z[:] = Z_march
        return Z[n]

    if debug:
        r0 = residual(fs_guess)
        print(f"FS_guess={fs_guess:.6f} → residual={r0:.4g}")

    # Root-find FS such that the right-end interslice force Z[n] = 0, by Newton-
    # secant from fs_guess. (A bracketed brentq solver was tried to also catch the
    # genuinely-low-FS surfaces where the single-guess secant diverges, but it
    # resurrects non-physical over-strength roots near FS->0 on surfaces where the
    # secant correctly fails, which the admissibility check below does not reliably
    # separate from legitimate low-FS solutions. Deferred.)
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
    # negative normals, <=20% interslice tension).
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

    # Report-only admissibility screen (corps/lowe inherit it). The parallel-force
    # march exposes no line of thrust, so only the base-tension and interslice-
    # tension signatures apply; the third is skipped. On right-facing slopes the
    # caller negates theta_list, flipping Z's sign, so tension there is Z>0 — pass
    # the physical-convention Z (tension < 0) the helper expects.
    Z_phys = -Z if right_facing else Z
    warns = _admissibility_warnings(c, N, Z_phys)

    if debug:
        r_opt = residual(FS_opt)
        print(f" Converged FS = {FS_opt:.6f}, residual = {r_opt:.4g}")

    return True, {'FS': FS_opt, 'warnings': warns}

def corps(slice_df, variant=2, debug=False):
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
        x0, y0 = slice_df['x_l'].iat[0], slice_df['y_lb'].iat[0]
        x1, y1 = slice_df['x_r'].iat[-1], slice_df['y_rb'].iat[-1]
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
    # inter-slice inclinations are negated — the same convention lowe
    # uses. This is a no-op for left-facing slopes (the common case and every
    # shipped benchmark). Without it, force-equilibrium FS is asymmetric under
    # mirroring and can violate the expected method ordering on right-facing
    # geometries.
    # Explicit facing (flat-arc / override) wins; else historical test, byte-identical.
    right_facing = slice_df.attrs.get(
        'right_facing', slice_df['y_lb'].iat[0] > slice_df['y_rb'].iat[-1])
    if right_facing:
        theta_list = -theta_list
        theta_out = -theta_out

    slice_df['theta'] = theta_list[:-1]  # store theta in slice_df. Adjust length to n slices.

    # delegate to your force_equilibrium solver
    success, results = force_equilibrium(slice_df, theta_list, debug=debug, right_facing=right_facing)
    if not success:
        return success, results
    else:
        results['method'] = 'corps'  # append method
        results['theta'] = theta_out           # append theta (mean for variant 2)
        return success, results

def lowe(slice_df, debug=False):
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

    # determine facing (explicit flat-arc / override wins; else historical test)
    right_facing = slice_df.attrs.get('right_facing', (y_lb[0] > y_rb[-1]))

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
        results['method'] = 'lowe'  # append method
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
            'y_lt'  (left-top y, used with 'y_lb' by the report-only
                     admissibility screen for the thrust-line check),
            'h_pile' (pile force, optional), 'theta_p' (pile inclination, RADIANS, optional)

    Returns:
        float: FS where FS_force = FS_moment
        float: beta (degrees)
        bool: converged flag
        The results dict also carries 'warnings': a list of admissibility
        notes (base tension on cohesionless slices, interslice tension,
        thrust line outside the slices) — empty when the solution is clean.
        Warnings never affect FS or convergence.
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

    # v12 support/load terms (all zero when the columns are absent)
    def _c12(name):
        return slice_df[name].values.astype(float) if name in slice_df.columns else np.zeros(n_slices)
    pa_cx = _c12('pa_cx'); pa_cy = _c12('pa_cy')
    pa_mx = _c12('pa_mx'); pa_my = _c12('pa_my')
    pp_cx = _c12('pp_cx'); pp_cy = _c12('pp_cy')
    pp_mx = _c12('pp_mx'); pp_my = _c12('pp_my')
    P_pt  = _c12('p_pt')
    H_pas = _c12('h_pile_pas')
    LL    = _c12('lload')
    ll_b  = np.radians(_c12('ll_beta'))
    ll_x  = _c12('ll_x'); ll_y = _c12('ll_y')

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

    # Explicit facing (flat-arc / override) wins; else historical y_cb test. On a
    # flat arc y_cb[0] and y_cb[-1] are near-equal, so this comparison is
    # otherwise noise-driven — see generate_slices / _resolve_right_facing.
    right_facing = slice_df.attrs.get('right_facing', (y_cb[0] > y_cb[-1]))
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
        V = -V   # tension-crack water thrust: flips with the sliding direction, like kw
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
        # Axial reinforcement enters in REAL components at real positions, and
        # the stored direction (cos(psi) > 0) is the NEGATION of the real
        # tension on a right-facing slope: negate the whole vector (all four
        # sums; verified against the left-facing mirror oracle). Line loads
        # mirror beta.
        pa_cx = -pa_cx; pa_cy = -pa_cy
        pa_mx = -pa_mx; pa_my = -pa_my
        pp_cx = -pp_cx; pp_cy = -pp_cy
        pp_mx = -pp_mx; pp_my = -pp_my
        P_pt = -P_pt          # tangent passive mirrors R
        ll_b = -ll_b

    # Passive pile portion: mobilized (divided by F), so it is removed from the
    # constant pile components and carried in the F-dependent passive arrays.
    frac_pas = np.divide(H_pas, H_pile, out=np.zeros_like(H_pas), where=H_pile != 0)
    H_cos_pas = H_cos_tp * frac_pas
    H_sin_pas = H_sin_tp * frac_pas
    H_cos_tp = H_cos_tp - H_cos_pas
    H_sin_tp = H_sin_tp - H_sin_pas

    Fh = (- kw - V + P * sin_b + R * cos_psi + H_cos_tp      # Equation (1) + pile horizontal
          + pa_cx + LL * np.sin(ll_b))                         # + axial reinf + line load
    Fv = (- W - P * cos_b + R * sin_psi + H_sin_tp             # Equation (2) + pile vertical (upward)
          + pa_cy - LL * np.cos(ll_b))                         # + axial reinf + line load

    # Passive (ultimate) support components: enter every equation divided by F.
    # Tangent-passive acts parallel to the base at the base center (like R);
    # axial-passive at its crossing point; passive piles per their components.
    Fh_pas = P_pt * cos_psi + pp_cx + H_cos_pas
    Fv_pas = P_pt * sin_psi + pp_cy + H_sin_pas
    Mo_pas = (- (pp_my - y_b * pp_cx) + (pp_mx - x_b * pp_cy)
              - H_cos_pas * (y_pile - y_b) + H_sin_pas * (x_pile - x_b))
    Mo = - P * sin_b * (y_p - y_b) - P * cos_b * (x_p - x_b) \
        + kw * (y_k - y_b) + V * (y_v - y_b) - R * cos_psi * (y_r - y_b) + R * sin_psi * (x_r - x_b) \
        - H_cos_tp * (y_pile - y_b) + H_sin_tp * (x_pile - x_b) \
        - (pa_my - y_b * pa_cx) + (pa_mx - x_b * pa_cy) \
        - LL * np.sin(ll_b) * (ll_y - y_b) - LL * np.cos(ll_b) * (ll_x - x_b)  # Eq (3) + pile/axial/line-load moments
    
    # ========== BEGIN SOLUTION ==========

    def compute_Q_and_yQ(F, theta_rad):
        """Compute Q and y_Q for given F and theta values. Passive support
        mobilizes with the soil: its components carry 1/F."""
        Fv_e = Fv + Fv_pas / F
        Fh_e = Fh + Fh_pas / F
        # Equation (24): m_alpha
        ma = 1 / (np.cos(alpha - theta_rad) + np.sin(alpha - theta_rad) * tan_p / F)

        # Equation (23): Q
        Q = (- Fv_e * sin_a - Fh_e * cos_a - (c / F) * dl + (Fv_e * cos_a - Fh_e * sin_a + u * dl) * tan_p / F) * ma

        # Equation (26): y_Q with numerical safeguard
        # Add small epsilon to prevent divide-by-zero when Q * cos(theta) is very small
        Q_cos_theta = Q * np.cos(theta_rad)
        eps = 1e-10
        safe_denom = np.where(np.abs(Q_cos_theta) < eps, eps * np.sign(Q_cos_theta + eps), Q_cos_theta)
        y_q = y_b + (Mo + Mo_pas / F) / safe_denom

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
        
        # Constants for Q expression (Equations 45-49), generalized: passive
        # support makes C1, C2 (and Mo) functions of F —
        #   C1(F) = C1_0 + C1_p/F,  C2(F) = C2_0 + C2_p/F,  Mo(F) = Mo + Mo_pas/F
        # so N(F) = C1 + C2/F and dN/dF = -(C1_p + C2_0 + 2*C2_p/F)/F^2, written
        # below as -(C1_p + C2 + C2_p/F)/F^2 using C2_0 = C2 - C2_p/F. With no
        # passive elements (C1_p = C2_p = 0) every expression reduces to the
        # original closed forms exactly.
        Fv_e = Fv + Fv_pas / F
        Fh_e = Fh + Fh_pas / F
        C1 = -Fv_e * sin_a - Fh_e * cos_a
        C2 = -c * dl + (Fv_e * cos_a - Fh_e * sin_a + u * dl) * tan_p
        C1_p = -Fv_pas * sin_a - Fh_pas * cos_a
        C2_p = (Fv_pas * cos_a - Fh_pas * sin_a) * tan_p
        C3 = cos_alpha_theta
        C4 = sin_alpha_theta * tan_p
        
        # Denominator for Q
        denom_Q = C3 + C4 / F
        
        # First-order partial derivatives of Q (Equations 50-51, generalized)
        dN_dF = -(C1_p + C2 + C2_p / F) / F**2
        dQ_dF = (1 / denom_Q**2) * (dN_dF * denom_Q + (C1 + C2 / F) * C4 / F**2)
        
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

        Mo_e = Mo + Mo_pas / F
        dyQ_dF = ((-Mo_pas / F**2) / safe_denom
                  + (-1 / safe_denom**2) * Mo_e * dQ_dF * cos_theta)
        dyQ_dtheta = (-1 / safe_denom**2) * Mo_e * (dQ_dtheta * cos_theta - Q * sin_theta)
        
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

    # --- Compute N_eff using Equation (18), with passive components at 1/F ---
    Fv_fin = Fv + Fv_pas / F
    Fh_fin = Fh + Fh_pas / F
    Mo_fin = Mo + Mo_pas / F
    N_eff = - Fv_fin * cos_a + Fh_fin * sin_a + Q * np.sin(alpha - theta_rad) - u * dl
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
            yt_r[i] = y_b[i] - ((Mo_fin[i] - Z[i] * sin_theta * dx[i] / 2 - Z[i+1] * sin_theta * dx[i] / 2 - Z[i] * cos_theta * (yt_l[i] - y_b[i])) / (Z[i+1] * cos_theta))
            yt_l[i+1] = yt_r[i]
    slice_df['yt_l'] = yt_l
    slice_df['yt_r'] = yt_r

    # --- Admissibility screen (report-only; acceptance above is unchanged) ---
    # max_tension_ratio normalizes base tension by cohesive capacity, so a c=0
    # slice can carry unlimited base tension without tripping it (VP30: -71 kN/m
    # on the cohesionless crack-face sliver scores 0.0 and passes). The shared
    # helper reports the Duncan & Wright signatures that guard cannot see; here
    # spencer keeps its own thrust line (yt_l), so all three fire. Z is already
    # in the physical convention (tension < 0) for both facings.
    warns = _admissibility_warnings(
        c, N_eff, Z,
        y_lt=slice_df['y_lt'].values,
        y_lb=slice_df['y_lb'].values,
        yt_l=yt_l)
    if warns and debug_level >= 1:
        print("Spencer admissibility warnings: " + "; ".join(warns))

    # --- Return results ---
    results = {}
    results['method'] = 'spencer'
    results['FS'] = F
    results['theta'] = theta_opt
    results['warnings'] = warns

    return True, results


def _mp_f_vals(slice_df, f_type):
    """Interslice-force shape function f evaluated at each slice BOUNDARY.

    The boundaries are x_b = [x_l[0], x_r[0], ..., x_r[n-1]] (length n+1), and the
    normalized position is s_j = (x_b[j] - x_L) / (x_R - x_L) in [0, 1]. Only the
    SHAPE of f matters (tan θ = λ·f, so scaling f just rescales λ), so each is
    normalized to a unit peak by convention.

    v1 library:
        'constant'  -> f = 1            (recovers Spencer; the regression case)
        'half_sine' -> f = sin(π s)     (textbook / GeoStudio default; f=0 at ends)
    """
    x_l = slice_df['x_l'].values.astype(float)
    x_r = slice_df['x_r'].values.astype(float)
    x_b = np.concatenate([[x_l[0]], x_r])           # boundary x-coords, length n+1
    span = x_b[-1] - x_b[0]
    s = (x_b - x_b[0]) / span if span != 0 else np.zeros_like(x_b)
    if f_type == 'constant':
        return np.ones_like(s)
    elif f_type == 'half_sine':
        return np.sin(np.pi * s)
    else:
        raise ValueError(f"unknown f_type {f_type!r} (use 'constant' or 'half_sine')")


def _mp_line_of_thrust(slice_df, Z, theta_rad, right_facing, FS=1.0):
    """Line of thrust for Morgenstern-Price (post-process diagnostic; plan §4 opt b).

    Generalizes Spencer's per-slice moment recurrence (its eqs 68-69) to a
    per-BOUNDARY interslice angle. Taking moments on each slice about its base
    center fixes where the right interslice force acts, given the left: the weight
    and base reaction act at the base center (thin-slice x_c = x_b → zero arm), so
    only the external-load moment `Mo` and the interslice `Z`/θ terms remain. With a
    constant θ this reduces EXACTLY to Spencer's recurrence (the f(x)=1 == Spencer
    check). Stores 'yt_l'/'yt_r' in slice_df.

    `theta_rad` = per-boundary interslice angle (length n+1, radians); `Z` = the
    interslice resultant at each boundary (length n+1) from the march. The same
    right-facing reflection the march/Spencer use is applied to the load terms.
    """
    n = len(slice_df)
    alpha = np.radians(slice_df['alpha'].values)
    beta  = np.radians(slice_df['beta'].values)
    P   = slice_df['dload'].values            # distributed load resultant
    kw  = slice_df['kw'].values
    V   = slice_df['t'].values                # tension-crack water
    R   = slice_df['p'].values                # reinforcement
    dx  = slice_df['dx'].values
    x_b = slice_df['x_c'].values; y_b = slice_df['y_cb'].values
    y_lb = slice_df['y_lb'].values; y_rb = slice_df['y_rb'].values
    x_p = slice_df['d_x'].values; y_p = slice_df['d_y'].values   # surface-load point
    y_k = slice_df['y_cg'].values; y_v = slice_df['y_t'].values
    H_pile     = slice_df['h_pile'].values  if 'h_pile'  in slice_df.columns else np.zeros(n)
    theta_pile = slice_df['theta_p'].values if 'theta_p' in slice_df.columns else np.zeros(n)
    x_pile = slice_df['x_pile'].values if 'x_pile' in slice_df.columns else np.zeros(n)
    y_pile = slice_df['y_pile'].values if 'y_pile' in slice_df.columns else np.zeros(n)
    psi = alpha; y_r = y_b; x_r = x_b           # reinforcement ∥ base, acts at base center

    # v12 support/load terms (zeros when the columns are absent)
    def _c12(name):
        return slice_df[name].values.astype(float) if name in slice_df.columns else np.zeros(n)
    P_pt  = _c12('p_pt')
    pa_cx = _c12('pa_cx'); pa_cy = _c12('pa_cy')
    pa_mx = _c12('pa_mx'); pa_my = _c12('pa_my')
    pp_cx = _c12('pp_cx'); pp_cy = _c12('pp_cy')
    pp_mx = _c12('pp_mx'); pp_my = _c12('pp_my')
    H_pas = _c12('h_pile_pas')
    LL    = _c12('lload')
    ll_b  = np.radians(_c12('ll_beta'))
    ll_x  = _c12('ll_x'); ll_y = _c12('ll_y')

    if right_facing:                            # same reflection as spencer()/_mp_march
        beta = -beta; psi = -psi; R = -R; kw = -kw; V = -V
        # axial: negate the whole force vector, as in spencer()
        pa_cx = -pa_cx; pa_cy = -pa_cy; pa_mx = -pa_mx; pa_my = -pa_my
        pp_cx = -pp_cx; pp_cy = -pp_cy; pp_mx = -pp_mx; pp_my = -pp_my
        P_pt = -P_pt                            # tangent passive mirrors R
        ll_b = -ll_b                            # line loads mirror beta
    sin_b, cos_b = np.sin(beta), np.cos(beta)
    sin_psi, cos_psi = np.sin(psi), np.cos(psi)
    frac_pas = np.divide(H_pas, H_pile, out=np.zeros_like(H_pas), where=H_pile != 0)
    H_cos_tp = H_pile * np.cos(theta_pile)
    H_sin_tp = H_pile * np.sin(theta_pile)
    if right_facing:
        H_cos_tp = -H_cos_tp                    # horizontal pile component flips, vertical not
    H_cos_pas = H_cos_tp * frac_pas
    H_sin_pas = H_sin_tp * frac_pas
    H_cos_tp = H_cos_tp - H_cos_pas
    H_sin_tp = H_sin_tp - H_sin_pas

    # Mo: moment of the external loads about the base center (Spencer eq 3 + pile
    # + v12 axial reinforcement, passive support factored by FS, and line loads).
    Mo = (- P * sin_b * (y_p - y_b) - P * cos_b * (x_p - x_b)
          + kw * (y_k - y_b) + V * (y_v - y_b)
          - R * cos_psi * (y_r - y_b) + R * sin_psi * (x_r - x_b)
          - H_cos_tp * (y_pile - y_b) + H_sin_tp * (x_pile - x_b)
          - (pa_my - y_b * pa_cx) + (pa_mx - x_b * pa_cy)
          + (- (pp_my - y_b * pp_cx) + (pp_mx - x_b * pp_cy)
             - P_pt * cos_psi * (y_r - y_b) + P_pt * sin_psi * (x_r - x_b)
             - H_cos_pas * (y_pile - y_b) + H_sin_pas * (x_pile - x_b)) / FS
          - LL * np.sin(ll_b) * (ll_y - y_b) - LL * np.cos(ll_b) * (ll_x - x_b))

    st, ct = np.sin(theta_rad), np.cos(theta_rad)   # per boundary (length n+1)
    yt_l = np.zeros(n)
    yt_r = np.zeros(n)
    yt_l[0] = y_lb[0]
    for i in range(n):
        if i == n - 1:
            yt_r[i] = y_rb[i]
        else:
            denom = Z[i + 1] * ct[i + 1]
            if denom == 0:
                yt_r[i] = y_b[i]
            else:
                yt_r[i] = y_b[i] - ((Mo[i]
                                     - Z[i] * st[i] * dx[i] / 2
                                     - Z[i + 1] * st[i + 1] * dx[i] / 2
                                     - Z[i] * ct[i] * (yt_l[i] - y_b[i])) / denom)
            yt_l[i + 1] = yt_r[i]
    slice_df['yt_l'] = yt_l
    slice_df['yt_r'] = yt_r


def mprice(slice_df, f_type='half_sine', fs_guess=1.5, tol=1e-6,
                      max_iter=50, lambda_bracket=(-1.5, 1.5), solver='auto',
                      debug_level=0):
    """Morgenstern-Price complete-equilibrium method (Approach A: F_f / F_m crossing).

    M-P satisfies BOTH force and moment equilibrium with a VARIABLE interslice-force
    inclination: ``tan(θ_j) = λ · f(x_j)``, where ``f`` is a prescribed shape
    function (see `_mp_f_vals`) and ``λ`` a scalar. Spencer is the special case
    ``f ≡ 1`` (then ``θ = arctan λ`` is constant).

    Approach A (Fredlund-Krahn / GLE style): for each ``λ`` find
        ``F_f(λ)`` = the FS root of the force residual  (Z[n] = 0), and
        ``F_m(λ)`` = the FS root of the moment residual (moment about origin = 0),
    then find ``λ*`` where ``F_f = F_m``. The crossing FS is the M-P factor of
    safety; ``λ*`` fixes the interslice-force distribution. Both residuals come from
    `_mp_march`.

    Parameters:
        f_type: 'constant' (== Spencer regression) or 'half_sine' (design default).
        lambda_bracket: search interval for ``λ`` (sign change of ``F_f − F_m``).
    Returns:
        (True, {'method','FS','lambda','f_type','theta'(deg, per boundary, n+1)})
        or (False, message).

    Right-facing slopes are handled: `_mp_extract` mirrors the per-slice arrays, and
    ``f ≡ 1`` reproduces `spencer()` to machine precision on both facings.
    """
    from scipy.optimize import brentq

    n = len(slice_df)
    if n < 2:
        return False, "Morgenstern-Price needs at least 2 slices."

    # Facing detection mirrors spencer(); _mp_march applies the right-facing flip set.
    # Explicit facing (flat-arc / override) wins; else the historical y_cb test.
    _ycb = slice_df['y_cb'].values
    right_facing = bool(slice_df.attrs.get('right_facing', _ycb[0] > _ycb[-1]))

    try:
        f_vals = _mp_f_vals(slice_df, f_type)
    except ValueError as e:
        return False, str(e)

    # Extract all per-slice arrays ONCE (with right-facing flips); the inner Newton
    # loops call _mp_residuals(A, ...) — no per-iteration DataFrame access.
    A = _mp_extract(slice_df, right_facing)

    # Seed FS from Bishop where possible (as Spencer does), else fs_guess.
    seed = fs_guess
    try:
        ok_b, res_b = bishop(slice_df)
        if ok_b:
            seed = res_b['FS']
    except Exception:
        pass

    def F_f(lam, s):
        """Force-equilibrium FS at this λ: the FS root of `force_res` (= Z[n]). The
        force residual is monotonic in FS, so this root is unique and well-behaved —
        no branch-jumping (unlike the moment-only FS, which is multivalued)."""
        return newton(lambda FS: _mp_residuals(A, f_vals, lam, FS)[2],
                      s, tol=tol, maxiter=max_iter)

    def h(lam):
        """Moment residual evaluated AT the force-equilibrium FS for this λ:
        ``h(λ) = moment_res(F_f(λ), λ)``. This is the key to a robust Approach A.
        `F_f(λ)` is smooth and single-valued, so `h` is a SMOOTH scalar function of λ
        alone whose root is the M-P solution (where force AND moment both close) — it
        sidesteps the multivalued moment-only FS curve `F_m(λ)` (and its asymptote)
        entirely. (This curve search is the slow-but-robust reference path; the
        Approach-B Newton, S5, is the fast path.)"""
        return _mp_residuals(A, f_vals, lam, F_f(lam, seed))[3]

    lo, hi = lambda_bracket

    def approach_A():
        """Robust reference (slow): scan h(λ) on a grid, brentq each sign change,
        and take the crossing nearest λ=0 (the physical one). ~61 h-evals."""
        grid = np.linspace(lo, hi, 61)
        hv = []
        for lam in grid:
            try:
                hv.append(h(lam))
            except Exception:
                hv.append(np.nan)
        hv = np.array(hv)
        crossings = []
        for i in range(len(grid) - 1):
            a, b = hv[i], hv[i + 1]
            if np.isfinite(a) and np.isfinite(b) and a * b <= 0 and not (a == 0 and b == 0):
                try:
                    crossings.append(brentq(h, grid[i], grid[i + 1], xtol=1e-9, maxiter=100))
                except Exception:
                    pass
        return min(crossings, key=abs) if crossings else None

    def approach_B():
        """Fast path (Approach B): 2-D Newton (scipy `hybr`) on the SCALED
        `(force_res, moment_res)` residuals over `(FS, λ)`, seeded at `(Bishop-FS,
        0)`. Both residuals come from ONE `_mp_march` per evaluation, so this is ~3×
        faster than the h(λ) secant (which nests a full `F_f` root-find inside every
        step). Scaling the two residuals to O(1) (they differ by ~1e3×) is what makes
        `hybr` well-conditioned; seeding near the physical solution keeps it off the
        multivalued moment-only branch. Returns λ* or None (→ fall back to A)."""
        from scipy.optimize import root
        f0 = _mp_residuals(A, f_vals, 0.0, seed)
        sf = abs(f0[2]) or 1.0       # force-residual scale at the seed
        sm = abs(f0[3]) or 1.0       # moment-residual scale at the seed

        def R(x):
            F, lam = x
            if not (np.isfinite(F) and np.isfinite(lam)) or F <= 0.05:
                return [1e3, 1e3]
            try:
                _, _, fr, mr = _mp_residuals(A, f_vals, lam, F)
            except Exception:
                return [1e3, 1e3]
            return [fr / sf, mr / sm]

        try:
            sol = root(R, [seed, 0.0], method='hybr')
        except Exception:
            return None
        F_b, lam_b = sol.x
        if not (sol.success and F_b > 0.05 and (lo - 0.5) <= lam_b <= (hi + 0.5)):
            return None
        rb = R(sol.x)                # confirm it's a true root, not a merit-fn stall
        if abs(rb[0]) < 1e-4 and abs(rb[1]) < 1e-4:
            return float(lam_b), float(F_b)   # λ* and its FS (skips a redundant re-solve)
        return None

    FS = None
    if solver == 'A':
        lam_star = approach_A()
    else:                       # 'auto' (default) / 'B': fast 2-D Newton, fall back to A
        b = approach_B()
        if b is not None:
            lam_star, FS = b           # B gives FS directly — no redundant F_f re-solve
        elif solver == 'auto':
            lam_star = approach_A()
        else:
            lam_star = None
    if lam_star is None:
        return False, ("Morgenstern-Price: no F_f/F_m crossing found in "
                       f"λ ∈ [{lo}, {hi}].")

    # Recover N, Z, θ at λ* (FS already known from Approach B; else solve it via A's F_f).
    if FS is None:
        FS = F_f(lam_star, seed)
    N, Z, force_res, moment_res = _mp_residuals(A, f_vals, lam_star, FS)
    theta_deg = np.degrees(np.arctan(lam_star * f_vals))   # per boundary, length n+1

    # Admissibility guard. M-P's λ-optimisation actively seeks the lowest FS, so —
    # unlike the fixed-θ methods — it can exploit surfaces that balance only with
    # partly TENSILE interslice forces (non-physical: soil cannot pull slices apart),
    # giving a spuriously low FS on a surface Spencer's solver simply fails to
    # converge on. We reject those. force_equilibrium's own note is that valid
    # criticals run ≤20% interslice tension; our physical M-P criticals are ≤18%,
    # while spurious ones (e.g. simple_embankment) hit ~33%, so a 30% INTERSLICE cap
    # cleanly separates them. Base-normal tension keeps the looser 50% extent
    # backstop (it is not what the λ-optimisation exploits). Interslice tension is
    # Z<0 here: M-P keeps the un-negated θ convention for BOTH facings (unlike
    # force_equilibrium, which negates theta_list for right-facing and flips Z).
    frac_N_neg = float(np.mean(N < 0))
    Z_int = Z[1:-1]
    frac_Z_tension = float(np.mean(Z_int < 0)) if len(Z_int) else 0.0
    if frac_N_neg > 0.5 or frac_Z_tension > 0.30:
        return False, ("Morgenstern-Price: inadmissible solution "
                       f"({100*frac_N_neg:.0f}% of base normals in tension, "
                       f"{100*frac_Z_tension:.0f}% interslice tension)")

    slice_df['n_eff'] = N
    slice_df['z'] = Z[:-1]
    slice_df['theta'] = theta_deg[1:]   # right-boundary interslice angle per slice
    _mp_line_of_thrust(slice_df, Z, np.arctan(lam_star * f_vals), right_facing, FS=FS)

    # Report-only admissibility screen (never affects FS/λ or acceptance). M-P
    # keeps the un-negated θ convention for both facings, so its Z is already
    # physical (tension < 0); the thrust line just stored gives the third
    # signature. On VP30 this fires interslice tension and thrust-outside — the
    # inadmissibility the λ-optimisation exploits, silent until now.
    warns = _admissibility_warnings(
        slice_df['c'].values, N, Z,
        y_lt=slice_df['y_lt'].values,
        y_lb=slice_df['y_lb'].values,
        yt_l=slice_df['yt_l'].values)

    if debug_level >= 1:
        print(f"Morgenstern-Price ({f_type}): FS = {FS:.5f}, λ = {lam_star:.5f}, "
              f"force_res = {force_res:.2e}, moment_res = {moment_res:.2e}")
        for w in warns:
            print(f"  M-P admissibility warning: {w}")

    return True, {
        'method': 'mprice',
        'FS': FS,
        'lambda': lam_star,
        'f_type': f_type,
        'theta': theta_deg,
        'warnings': warns,
    }


# =====================================================================
# Nonlinear strength coupling — options 'pow' (v12) and 'hb' (v14)
#
# Both envelopes are curved, so the shear strength depends on the base normal
# stress, which itself depends on FS. Each is handled the same way: linearize
# the envelope at the current sigma'_n into an instantaneous Mohr-Coulomb
# tangent (c_i, phi_i), solve, recompute sigma'_n, repeat until FS is
# stationary (_with_nonlinear_strength).
# =====================================================================

def _hb_update_strength(slice_df, sigma_n):
    """Re-linearize the generalized Hoek-Brown envelope at the current base
    normal stress into an instantaneous tangent (c_i, phi_i) for every
    hb-flagged slice (Balmer transformation; see xslope.hoekbrown)."""
    m = slice_df['hb_flag'].values.astype(bool)
    if not m.any():
        return
    s = np.maximum(np.asarray(sigma_n, dtype=float)[m], 0.0)
    c_i, phi_i = hb_tangent(s,
                            slice_df['hb_sci'].values[m],
                            slice_df['hb_gsi'].values[m],
                            slice_df['hb_mi'].values[m],
                            slice_df['hb_d'].values[m])
    c_vals = slice_df['c'].values.copy()
    phi_vals = slice_df['phi'].values.copy()
    c_vals[m] = c_i
    phi_vals[m] = phi_i
    slice_df['c'] = c_vals
    slice_df['phi'] = phi_vals


def _pow_update_strength(slice_df, sigma_n):
    """Re-linearize the power-curve envelope tau = a*(sigma'_n + d)^b + c_p at
    the current base normal stress into an instantaneous-tangent (c_i, phi_i)
    for every pow-flagged slice: phi_i = arctan(a*b*(sigma+d)^(b-1)) and c_i =
    tau(sigma) - sigma*tan(phi_i). sigma is clamped at zero (no tensile normal
    stress) and (sigma + d) is floored at a small fraction of the mean stress so
    b < 1 curves cannot produce an unbounded tangent slope at the origin."""
    m = slice_df['pow_flag'].values.astype(bool)
    if not m.any():
        return
    a = slice_df['pow_a'].values[m]
    b = slice_df['pow_b'].values[m]
    cp = slice_df['pow_c'].values[m]
    d = slice_df['pow_d'].values[m]
    s = np.maximum(np.asarray(sigma_n, dtype=float)[m], 0.0)
    ref = max(1.0, float(np.mean(s)) if len(s) else 1.0)
    s_eff = np.maximum(s + d, 1e-4 * ref)
    slope = a * b * np.power(s_eff, b - 1.0)
    tau = a * np.power(s_eff, b) + cp
    c_vals = slice_df['c'].values.copy()
    phi_vals = slice_df['phi'].values.copy()
    c_vals[m] = tau - s * slope
    phi_vals[m] = np.degrees(np.arctan(slope))
    slice_df['c'] = c_vals
    slice_df['phi'] = phi_vals


def _has_nonlinear(slice_df):
    """True if any slice sits on a curved envelope ('pow' or 'hb')."""
    for flag in ('pow_flag', 'hb_flag'):
        if flag in slice_df.columns and slice_df[flag].any():
            return True
    return False


def _with_nonlinear_strength(method):
    """Outer fixed-point wrapper giving every solution method curved-envelope
    support ('pow' and 'hb'): solve with the current tangent strengths, update
    sigma'_n from the method's n_eff, re-tangent, repeat until FS is stationary.
    slice_df arrives seeded with tangent strengths at a no-FS normal-stress
    estimate (slice.py), so the first inner solve is already close. Damped (50%)
    sigma updates prevent oscillation on strongly curved envelopes."""
    import functools

    @functools.wraps(method)
    def wrapped(slice_df, *args, **kwargs):
        if not _has_nonlinear(slice_df):
            return method(slice_df, *args, **kwargs)
        dl = slice_df['dl'].values.astype(float)
        sigma = None
        fs_prev = None
        for _ in range(40):
            ok, res = method(slice_df, *args, **kwargs)
            if not ok:
                return ok, res
            n_eff = slice_df['n_eff'].values.astype(float)
            sigma_new = np.where(dl > 0, n_eff / np.where(dl > 0, dl, 1.0), 0.0)
            sigma = sigma_new if sigma is None else 0.5 * sigma + 0.5 * sigma_new
            if 'pow_flag' in slice_df.columns:
                _pow_update_strength(slice_df, sigma)
            if 'hb_flag' in slice_df.columns:
                _hb_update_strength(slice_df, sigma)
            if fs_prev is not None and abs(res['FS'] - fs_prev) < 1e-4:
                return ok, res
            fs_prev = res['FS']
        return False, ("Nonlinear strength iteration did not converge in 40 outer "
                       "iterations (envelope tangent vs normal stress).")
    return wrapped


# Rebind the public solution methods with the nonlinear-strength outer loop.
# Import-time rebinding covers both `solve.oms(...)` attribute lookups and
# `from xslope.solve import oms` (imports resolve after module execution).
oms = _with_nonlinear_strength(oms)
bishop = _with_nonlinear_strength(bishop)
janbu = _with_nonlinear_strength(janbu)
spencer = _with_nonlinear_strength(spencer)
corps = _with_nonlinear_strength(corps)
lowe = _with_nonlinear_strength(lowe)
mprice = _with_nonlinear_strength(mprice)
