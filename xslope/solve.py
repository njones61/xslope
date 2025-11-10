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

def solve_all(slice_df):
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
    solve_selected('oms', slice_df)
    solve_selected('bishop', slice_df)
    solve_selected('janbu', slice_df)
    solve_selected('corps_engineers', slice_df)
    solve_selected('lowe_karafiath', slice_df)
    solve_selected('spencer', slice_df)

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
          'd'             = resultant distributed load Dᵢ
          'd_x','d_y'     = centroid (x,y) at which Dᵢ acts
          'beta'   (deg)   = top slope βᵢ
          'kw'            = seismic horizontal kWᵢ
          't'             = tension‐crack horizontal Tᵢ  (zero except one slice)
          'y_t'           = y‐loc of Tᵢ's line of action (zero except that one slice)
          'p'             = reinforcement uplift pᵢ (zero if none)
          'x_c','y_cg'    = slice‐centroid (x,y) for seismic moment arm
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
    # 5) Build the NUMERATOR = Σᵢ [  cᵢ·Δℓᵢ
    #                               + (Wᵢ·cosαᵢ + Dᵢ·cos(αᵢ−βᵢ) − kWᵢ·sinαᵢ − Tᵢ·sinαᵢ − uᵢ·Δℓᵢ )·tanφᵢ
    #                               + pᵢ  ] + Σ  Dᵢ·sinβᵢ·(Yo - d_{y,i}) 
    #


    # N′ᵢ = Wᵢ·cosαᵢ + Dᵢ·cos(αᵢ−βᵢ) − kWᵢ·sinαᵢ − Tᵢ·sinαᵢ − uᵢ·Δℓᵢ
    N_eff = (
        W * cos_alpha
      + D * cos_ab
      - kw * sin_alpha
      - T * sin_alpha
      - (u * dl)
    )  

    #   Σ  Dᵢ·sinβᵢ·(Yo - d_{y,i}) 
    a_dy = Yo - d_y
    sum_Dy = np.sum(D * np.sin(beta) * a_dy)

    numerator = np.sum(c * dl + N_eff * tan_phi + P)+ (1.0 / R) * sum_Dy

    # ————————————————————————————————————————————————————————
    # 6) Build each piece of the DENOMINATOR exactly as Eqn 9:

    #  (A) = Σ [ Wᵢ · sinαᵢ ]
    sum_W = np.sum(W * sin_alpha)

    #  (B) = Σ  Dᵢ·cosβᵢ·(Xo - d_{x,i}) 
    a_dx = d_x - Xo
    sum_Dx = np.sum(D * np.cos(beta) * a_dx)

    #  (C) = Σ [ kWᵢ * (Yo - y_{cg,i}) ]
    a_s = Yo - y_cg
    sum_kw = np.sum(kw * a_s)

    #  (D) = Σ [ Tᵢ * (Yo - y_{t,i}) ]
    a_t = Yo - y_t
    sum_T = np.sum(T * a_t)

    # Put them together with their 1/R factors:
    denominator = sum_W + (1.0 / R) * (sum_Dx + sum_kw + sum_T)

    # 7) Finally compute FS = (numerator)/(denominator)
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

    # Trigonometric terms
    sin_alpha = np.sin(alpha)
    cos_alpha = np.cos(alpha)
    tan_phi   = np.tan(phi)
    sin_beta  = np.sin(beta)
    cos_beta  = np.cos(beta)

    # Moment arms
    a_dx = d_x - Xo
    a_dy = Yo - d_y
    a_s  = Yo - y_cg
    a_t  = Yo - y_t

    # Denominator (moment equilibrium)
    sum_W = np.sum(W * sin_alpha)
    sum_Dx = np.sum(D * cos_beta * a_dx)
    sum_Dy = np.sum(D * sin_beta * a_dy)
    sum_kw = np.sum(kw * a_s)
    sum_T = np.sum(T * a_t)
    denominator = sum_W + (1.0 / R) * (sum_Dx + sum_kw + sum_T)

    # Iterative solution
    F = 1.0
    for _ in range(max_iter):
        # Compute N_eff from Equation (8)
        num_N = (
            W + D * cos_beta - P * sin_alpha
            - u * dl * cos_alpha
            - (c * dl * sin_alpha) / F
        )
        denom_N = cos_alpha + (sin_alpha * tan_phi) / F
        N_eff = num_N / denom_N

        # Numerator for FS from Equation (10)
        shear = (
            c * dl * cos_alpha
            + (W + D * cos_beta - P * sin_alpha - u * dl * cos_alpha) * tan_phi
            + P
        )
        numerator = np.sum(shear / denom_N) + (1.0 / R) * sum_Dy
        F_new = numerator / denominator

        if abs(F_new - F) < tol:
            slice_df['n_eff'] = N_eff
            if debug:
                print(f"FS = {F_new:.6f}")
                print(f"Numerator = {numerator:.6f}")
                print(f"Denominator = {denominator:.6f}")
                print("N_eff =", np.array2string(N_eff, precision=4, separator=', '))
            return True, {'method': 'bishop', 'FS': F_new}

        F = F_new

    return False, "Bishop method did not converge within the maximum number of iterations."

def janbu(slice_df, debug=False):
    """
    Computes FS using Janbu's Simplified Method with correction factor (Equation 7).

    Implements the complete formulation including distributed loads, seismic forces,
    reinforcement, and tension crack water forces. Applies Janbu correction factor
    based on d/L ratio and soil type.

    Parameters:
        slice_df : pandas.DataFrame with required columns (see OMS spec)
        debug : bool, if True prints diagnostic info

    Returns:
        (bool, dict | str): (True, {'method': 'janbu_simplified', 'FS': value, 'fo': correction_factor})
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

    # Trigonometric terms
    sin_alpha = np.sin(alpha)
    cos_alpha = np.cos(alpha)
    tan_phi = np.tan(phi)
    sin_beta_alpha = np.sin(beta - alpha)
    cos_beta_alpha = np.cos(beta - alpha)

    # Effective normal forces (Equation 10)
    N_eff = W * cos_alpha - kw * sin_alpha + D * cos_beta_alpha - T * sin_alpha - u * dl

    # Numerator: resisting forces (shear resistance)
    numerator = np.sum(c * dl + N_eff * tan_phi + P)

    # Denominator: driving forces parallel to base (Equation 6)
    denominator = np.sum(W * sin_alpha + kw * cos_alpha - D * sin_beta_alpha + T * cos_alpha)

    # Base factor of safety (Equation 7)
    if abs(denominator) < 1e-12:
        return False, "Division by zero in Janbu method: driving forces sum to zero"

    FS_base = numerator / denominator

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

    # Determine b1 factor based on soil type
    phi_sum = slice_df['phi'].sum()
    c_sum = slice_df['c'].sum()

    if phi_sum == 0:  # c-only soil (undrained, φ = 0)
        b1 = 0.67
    elif c_sum == 0:  # φ-only soil (no cohesion)
        b1 = 0.31
    else:  # c-φ soil
        b1 = 0.50

    # Correction factor
    fo = 1 + b1 * (dL_ratio - 1.4 * dL_ratio ** 2)

    # Final corrected factor of safety
    FS = FS_base * fo

    # Store effective normal forces in DataFrame
    slice_df['n_eff'] = N_eff

    if debug:
        print(f"FS_base = {FS_base:.6f}")
        print(f"d/L ratio = {dL_ratio:.4f}")
        print(f"b1 factor = {b1:.2f}")
        print(f"fo correction = {fo:.4f}")
        print(f"FS_corrected = {FS:.6f}")
        print(f"Numerator = {numerator:.6f}")
        print(f"Denominator = {denominator:.6f}")
        print("N_eff =", np.array2string(N_eff, precision=4, separator=', '))

    return True, {
        'method': 'janbu',
        'FS': FS,
        'fo': fo
    }


def force_equilibrium(slice_df, theta_list, fs_guess=1.5, tol=1e-6, max_iter=50, debug=False):
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
            'd'     (distributed load),
            'beta'  (distributed load inclination, degrees),
            'kw'    (seismic force),
            't'     (tension crack water force),
            'p'     (reinforcement force)
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
            b0 = (
                -c_m[i]*dl[i]*ca 
                - P[i]*ca 
                + u[i]*dl[i]*sa 
                - Z[i]*np.cos(theta[i]) 
                - D[i]*sb 
                + kw[i] 
                + T[i]
            )
            b1 = (
                -c_m[i]*dl[i]*sa 
                - P[i]*sa 
                - u[i]*dl[i]*ca 
                + w[i] 
                - Z[i]*np.sin(theta[i]) 
                + D[i]*cb
            )
            
            N_i, Z_ip1 = np.linalg.solve(A, np.array([b0, b1]))
            Z[i+1] = Z_ip1
            N[i] = N_i  # store normal force on slice base
        return Z[n]

    if debug:
        r0 = residual(fs_guess)
        print(f"FS_guess={fs_guess:.6f} → residual={r0:.4g}")

    # use Newton‐secant (no derivative) with single initial guess
    try:
        FS_opt = newton(residual, fs_guess, tol=tol, maxiter=max_iter)
    except Exception as e:
        return False, f"force_equilibrium failed to converge: {e}"

    slice_df['n_eff'] = N  # store effective normal forces in slice_df
    slice_df['z'] = Z[:-1]  # store interslice forces in slice_df, adjust length to n slices

    if debug:
        r_opt = residual(FS_opt)
        print(f" Converged FS = {FS_opt:.6f}, residual = {r_opt:.4g}")

    return True, {'FS': FS_opt}

def corps_engineers(slice_df, debug=False):
    """
    Corps of Engineers style force equilibrium solver.

    1. Computes a single θ from the slope between
       (x_l[0], y_lt[0]) and (x_r[-1], y_rt[-1]).
    2. Builds a constant θ array of length n+1.
    3. Calls force_equilibrium(slice_df, theta_array).

    Parameters:
        slice_df (pd.DataFrame): Must include at least ['x_l','y_lt','x_r','y_rt']
                           plus all the columns required by force_equilibrium:
                           ['alpha','phi','c','dl','w','u','dx'].

    Returns:
        Tuple(bool, dict or str): Whatever force_equilibrium returns.
    """
    # endpoints of the slip surface
    x0, y0 = slice_df['x_l'].iat[0], slice_df['y_lt'].iat[0]
    x1, y1 = slice_df['x_r'].iat[-1], slice_df['y_rt'].iat[-1]

    # compute positive slope‐angle
    dx = x1 - x0
    dy = y1 - y0
    if abs(dx) < 1e-12:
        theta_deg = 90.0
    else:
        theta_deg = abs(np.degrees(np.arctan2(dy, dx)))

    # one theta per slice boundary
    n = len(slice_df)
    theta_list = np.full(n+1, theta_deg)

    slice_df['theta'] = theta_list[:-1]  # store theta in slice_df. Adjust length to n slices.

    # delegate to your force_equilibrium solver
    success, results = force_equilibrium(slice_df, theta_list, debug=debug)
    if not success:
        return success, results
    else:
        results['method'] = 'corps_engineers'  # append method
        results['theta'] = theta_deg           # append theta
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
    success, results = force_equilibrium(slice_df, theta_list, debug=debug)
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
            'd'     (distributed load),
            'beta'  (distributed load inclination, degrees),
            'kw'    (seismic force),
            't'     (tension crack water force),
            'p'     (reinforcement force)

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

    Fh = - kw - V + P * sin_b + R * cos_psi       # Equation (1)
    Fv = - W - P * cos_b + R * sin_psi        # Equation (2)
    Mo = - P * sin_b * (y_p - y_b) - P * cos_b * (x_p - x_b) \
        + kw * (y_k - y_b) + V * (y_v - y_b) - R * cos_psi * (y_r - y_b) + R * sin_psi * (x_r - x_b) # Equation (3)
    
    # ========== BEGIN SOLUTION ==========
    
    def compute_Q_and_yQ(F, theta_rad):
        """Compute Q and y_Q for given F and theta values."""
        # Equation (24): m_alpha
        ma = 1 / (np.cos(alpha - theta_rad) + np.sin(alpha - theta_rad) * tan_p / F)
        
        # Equation (23): Q
        Q = (- Fv * sin_a - Fh * cos_a - (c / F) * dl + (Fv * cos_a - Fh * sin_a + u * dl) * tan_p / F) * ma
        
        # Equation (26): y_Q
        y_q = y_b + Mo / (Q * np.cos(theta_rad))
        
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
        dyQ_dF = (-1 / (Q * cos_theta)**2) * Mo * dQ_dF * cos_theta
        dyQ_dtheta = (-1 / (Q * cos_theta)**2) * Mo * (dQ_dtheta * cos_theta - Q * sin_theta)
        
        # First-order partial derivatives of R1 (Equations 35-36)
        dR1_dF = np.sum(dQ_dF)
        dR1_dtheta = np.sum(dQ_dtheta)
        
        # First-order partial derivatives of R2 (Equations 40-41)
        dR2_dF = np.sum(dQ_dF * (x_b * sin_theta - y_q * cos_theta)) - np.sum(Q * dyQ_dF * cos_theta)
        dR2_dtheta = np.sum(dQ_dtheta * (x_b * sin_theta - y_q * cos_theta)) + np.sum(Q * (x_b * cos_theta + y_q * sin_theta - dyQ_dtheta * cos_theta))
        
        return dR1_dF, dR1_dtheta, dR2_dF, dR2_dtheta

    
    # Initial guesses
    F0 = 1.5
    if right_facing:
        theta0_rad = np.radians(-8.0) 
    else:
        theta0_rad = np.radians(8) 
    
    # Newton iteration
    F = F0
    theta_rad = theta0_rad
    
    for iteration in range(max_iter):
        # Compute residuals
        R1, R2, Q, y_q = compute_residuals(F, theta_rad)
        
        if debug_level >= 1:
            if iteration == 0:
                print(f"Iteration {1} - Initial: F = {F:.3f}, theta = {np.degrees(theta_rad):.3f}°, R1 = {R1:.6e}, R2 = {R2:.6e}")
            else:
                print(f"Iteration {iteration + 1} - Updated: F = {F:.3f}, theta = {np.degrees(theta_rad):.3f}°, R1 = {R1:.6e}, R2 = {R2:.6e}")
        
        # Check convergence
        if abs(R1) < tol and abs(R2) < tol:
            if debug_level >= 1:
                print(f"Converged in {iteration + 1} iterations, R1 = {R1:.6e}, R2 = {R2:.6e}")
            break
        
        # Compute derivatives
        dR1_dF, dR1_dtheta, dR2_dF, dR2_dtheta = compute_derivatives(F, theta_rad, Q, y_q)
        
        # Basic Newton method (Equations 31-32)
        # Build Jacobian matrix
        J = np.array([[dR1_dF, dR1_dtheta], 
                      [dR2_dF, dR2_dtheta]])
        
        # Check condition number for numerical stability
        try:
            cond_num = np.linalg.cond(J)
            if cond_num > 1e12:
                return False, f"Ill-conditioned Jacobian matrix (condition number: {cond_num:.2e})"
        except:
            return False, "Unable to compute Jacobian condition number"
        
        # Solve using matrix form for better numerical stability
        try:
            delta_solution = np.linalg.solve(J, np.array([-R1, -R2]))
            delta_F = delta_solution[0]
            delta_theta = delta_solution[1]
        except np.linalg.LinAlgError:
            return False, "Singular Jacobian matrix in Newton iteration"

        if debug_level >= 1:
            print(f"          Newton: delta_F = {delta_F:.3f}, delta_theta = {np.degrees(delta_theta):.3f}°, {delta_theta: .3f} (rad)")

        # Add step size control to prevent large jumps
        max_delta_F = 0.5  # Maximum allowed change in F per iteration
        max_delta_theta = np.radians(20)  # Maximum allowed change in theta per iteration (20 degrees)
        
        # Apply step size limiting
        if abs(delta_F) > max_delta_F:
            delta_F = np.sign(delta_F) * max_delta_F
            if debug_level >= 1:
                print(f"          Step limited: delta_F clamped to {delta_F:.3f}")
                
        if abs(delta_theta) > max_delta_theta:
            delta_theta = np.sign(delta_theta) * max_delta_theta
            if debug_level >= 1:
                print(f"          Step limited: delta_theta clamped to {np.degrees(delta_theta):.3f}°")

        # Update values
        F += delta_F
        theta_rad += delta_theta
        
        # Ensure F stays positive
        if F <= 0:
            F = 0.1
        
        # Limit theta to reasonable range
        theta_rad = np.clip(theta_rad, -np.pi/2, np.pi/2)
    
    # Check if we converged
    if iteration >= max_iter - 1:
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


def spencer_OLD(df, tol=1e-4, max_iter = 100, debug_level=2):
    """
    Spencer's Method using Steve G. Wright's formulation from the UTEXAS v2  user manual.
    

    Parameters:
        df (pd.DataFrame): must contain columns
            'alpha' (slice base inclination, degrees),
            'phi'   (slice friction angle, degrees),
            'c'     (cohesion),
            'dl'    (slice base length),
            'w'     (slice weight),
            'u'     (pore force per unit length),
            'd'     (distributed load),
            'beta'  (distributed load inclination, degrees),
            'kw'    (seismic force),
            't'     (tension crack water force),
            'p'     (reinforcement force)

    Returns:
        float: FS where FS_force = FS_moment
        float: beta (degrees)
        bool: converged flag
    """

    alpha = np.radians(df['alpha'].values)  # slice base inclination, degrees
    phi   = np.radians(df['phi'].values)  # slice friction angle, degrees  
    c     = df['c'].values  # cohesion
    dx    = df['dx'].values  # slice width
    dl    = df['dl'].values  # slice base length
    W     = df['w'].values  # slice weight
    u     = df['u'].values  # pore presssure
    x_c   = df['x_c'].values  # center of base x-coordinate
    y_cb  = df['y_cb'].values  # center of base y-coordinate
    y_lb   = df['y_lb'].values  # left side base y-coordinate
    y_rb   = df['y_rb'].values  # right side base y-coordinate
    P     = df['dload'].values  # distributed load resultant 
    beta  = np.radians(df['beta'].values)  # distributed load inclination, degrees
    kw    = df['kw'].values  # seismic force
    V     = df['t'].values  # tension crack water force
    y_v   = df['y_t'].values  # tension crack water force y-coordinate
    R     = df['p'].values  # reinforcement force

    # For now, we assume that reinforcement is flexible and therefore is parallel to the failure surface
    # at the bottom of the slice. Therefore, the psi value used in the derivation is set to alpha, 
    # and the point of action is the center of the base of the slice.
    psi = alpha  # psi is the angle of the reinforcement force from the horizontal
    y_r = y_cb  # y_r is the y-coordinate of the point of action of the reinforcement
    x_r = x_c  # x_r is the x-coordinate of the point of action of the reinforcement

    # use variable names to match the derivation.
    x_p = df['d_x'].values  # distributed load x-coordinate
    y_p = df['d_y'].values  # distributed load y-coordinate
    y_k = df['y_cg'].values  # seismic force y-coordinate
    x_b = x_c  # center of base x-coordinate
    y_b = y_cb  # center of base y-coordinate

    # pre-compute the trigonometric functions
    cos_a = np.cos(alpha)  # cos(alpha)
    sin_a = np.sin(alpha)  # sin(alpha)
    tan_p = np.tan(phi)  # tan(phi)
    cos_b = np.cos(beta)  # cos(beta)
    sin_b = np.sin(beta)  # sin(beta)
    sin_psi = np.sin(psi)  # sin(psi)
    cos_psi = np.cos(psi)  # cos(psi)

    Fh = - kw - V + P * sin_b + R * cos_psi       # Equation (1)
    Fv = - W - P * cos_b + R * sin_psi        # Equation (2)
    Mo = - P * sin_b * (y_p - y_b) - P * cos_b * (x_p - x_b) \
        + kw * (y_k - y_b) + V * (y_v - y_b) - R * cos_psi * (y_r - y_b) + R * sin_psi * (x_r - x_b) # Equation (3)

    def compute_Q(F, theta_rad):
        ma = 1 / (np.cos(alpha - theta_rad) + np.sin(alpha - theta_rad) * tan_p / F)  # Equation (24)
        Q = (- Fv * sin_a - Fh * cos_a - (c / F) * dl + (Fv * cos_a - Fh * sin_a + u * dl) * tan_p / F) * ma     # Equation (23)
        y_q = y_b + Mo / (Q * np.cos(theta_rad))   # Equation (26)
        return Q, y_q

    fs_min = 0.01
    fs_max = 20.0

    def fs_force(theta_rad):
        def residual(F):
            Q, y_q = compute_Q(F, theta_rad)
            return Q.sum()  # Equation (15)
        result = minimize_scalar(lambda F: abs(residual(F)), bounds=(fs_min, fs_max), method='bounded', options={'xatol': tol})
        return result.x

    def fs_moment(theta_rad):
        def residual(F):
            Q, y_q = compute_Q(F, theta_rad)
            return np.sum(Q * (x_b * np.sin(theta_rad) - y_q * np.cos(theta_rad)))  # Equation (16)
        result = minimize_scalar(lambda F: abs(residual(F)), bounds=(fs_min, fs_max), method='bounded', options={'xatol': tol})
        return result.x

    def fs_difference(theta_deg):
        theta_rad = np.radians(theta_deg)
        Ff = fs_force(theta_rad)
        Fm = fs_moment(theta_rad)
        return Ff - Fm

    # Robust theta root-finding with multiple strategies
    theta_opt = None
    convergence_error = None
    
    # Strategy 1: Try multiple starting points for Newton's method

    newton_starting_points = [3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
    
    # Pre-evaluate fs_difference for all starting points and sort by absolute value
    starting_point_evaluations = []
    for theta_guess in newton_starting_points:
        try:
            fs_diff = fs_difference(theta_guess)
            starting_point_evaluations.append((theta_guess, abs(fs_diff), fs_diff))
        except Exception as e:
            if debug_level >= 1:
                print(f"Failed to evaluate fs_difference at {theta_guess:.1f} deg: {e}")
            continue
    
    # Sort by absolute value of fs_difference (smallest first)
    starting_point_evaluations.sort(key=lambda x: x[1])
    
    if debug_level >= 1:
        print("Starting points sorted by |fs_difference|:")
        for theta_guess, abs_fs_diff, fs_diff in starting_point_evaluations:
            print(f"  {theta_guess:.1f}°: |fs_diff| = {abs_fs_diff:.6f}, fs_diff = {fs_diff:.6f}")
    
    for theta_guess, abs_fs_diff, fs_diff in starting_point_evaluations:
        try:
            if debug_level >= 1:
                print(f"Trying Newton's method with initial guess {theta_guess:.1f} deg (|fs_diff| = {abs_fs_diff:.6f})")
            theta_candidate = newton(fs_difference, x0=theta_guess, tol=tol, maxiter=max_iter)
            
            # Check if the solution is valid
            if (abs(theta_candidate) <= 59 and 
                abs(fs_difference(theta_candidate)) <= 0.01 and
                fs_force(np.radians(theta_candidate)) < fs_max - 1e-3):
                theta_opt = theta_candidate
                if debug_level >= 1:
                    print(f"Newton's method succeeded with starting point {theta_guess:.1f} deg")
                break
        except Exception as e:
            if debug_level >= 1:
                print(f"Newton's method failed with starting point {theta_guess:.1f} deg: {e}")
            continue
    
    # Strategy 2: If Newton's method failed, try adaptive grid search
    if theta_opt is None:
        if debug_level >= 1:
            print("Newton's method failed for all starting points, trying adaptive grid search...")
        
        # First, do a coarse sweep to identify promising regions
        theta_coarse = np.linspace(-60, 60, 121)  # More points for better resolution
        fs_diff_coarse = []
        
        for theta in theta_coarse:
            try:
                fs_diff_coarse.append(fs_difference(theta))
            except Exception:
                fs_diff_coarse.append(np.nan)
        
        fs_diff_coarse = np.array(fs_diff_coarse)
        
        # Find regions where sign changes occur
        sign_changes = []
        for i in range(len(fs_diff_coarse) - 1):
            if (not np.isnan(fs_diff_coarse[i]) and 
                not np.isnan(fs_diff_coarse[i+1]) and
                fs_diff_coarse[i] * fs_diff_coarse[i+1] < 0):
                sign_changes.append((theta_coarse[i], theta_coarse[i+1]))
        
        # Try root_scalar on each bracket
        for bracket in sign_changes:
            try:
                if debug_level >= 1:
                    print(f"Trying root_scalar with bracket {bracket}")
                sol = root_scalar(fs_difference, bracket=bracket, method='brentq', xtol=tol)
                theta_candidate = sol.root
                
                # Check if the solution is valid
                if (abs(theta_candidate) <= 59 and 
                    abs(fs_difference(theta_candidate)) <= 0.01 and
                    fs_force(np.radians(theta_candidate)) < fs_max - 1e-3):
                    theta_opt = theta_candidate
                    if debug_level >= 1:
                        print(f"root_scalar succeeded with bracket {bracket}")
                    break
            except Exception as e:
                if debug_level >= 1:
                    print(f"root_scalar failed with bracket {bracket}: {e}")
                continue
    
    # Strategy 3: If still no solution, try global optimization
    if theta_opt is None:
        if debug_level >= 1:
            print("All root-finding methods failed, trying global optimization...")
        
        try:
            # Use minimize_scalar to find the minimum of |fs_difference|
            result = minimize_scalar(
                lambda theta: abs(fs_difference(theta)), 
                bounds=(-60, 60), 
                method='bounded', 
                options={'xatol': tol}
            )
            
            if result.success and abs(fs_difference(result.x)) <= 0.01:
                theta_opt = result.x
                if debug_level >= 1:
                    print(f"Global optimization succeeded with theta = {theta_opt:.6f} deg")
            else:
                convergence_error = f"Global optimization failed: {result.message}"
                
        except Exception as e:
            convergence_error = f"Global optimization failed: {e}"
    
    # Check if we found a solution
    if theta_opt is None:
        if convergence_error:
            return False, f"Spencer's method failed to converge: {convergence_error}"
        else:
            return False, "Spencer's method: No valid solution found with any method."

    theta_rad = np.radians(theta_opt)
    FS_force = fs_force(theta_rad)
    FS_moment = fs_moment(theta_rad)

    df['theta'] = theta_opt  # store theta in df.

    # --- Compute N_eff ---
    Q, y_q = compute_Q(FS_force, theta_rad)
    N_eff = - Fv * cos_a + Fh * sin_a + Q * np.sin(alpha - theta_rad) - u * dl   # Equation (18)

    # ---  compute interslice forces Z  ---
    n = len(Q)
    Z = np.zeros(n+1)
    for i in range(n):
        Z[i+1] = Z[i] - Q[i] 

    # --- store back into df ---
    df['z']     = Z[:-1]        # Z_i acting on slice i's left face
    df['n_eff'] = N_eff

    # --- compute line of thrust ---
    yt_l = np.zeros(n)  # the y-coordinate of the line of thrust on the left side of the slice.
    yt_r = np.zeros(n)  # the y-coordinate of the line of thrust on the right side of the slice.
    yt_l[0] = y_lb[0]  
    sin_theta = np.sin(theta_rad)
    cos_theta = np.cos(theta_rad)
    for i in range(n):
        if i == n - 1:
            yt_r[i] = y_rb[i]
        else:
            yt_r[i] = y_b[i] - ((Mo[i] - Z[i] * sin_theta * dx[i] / 2 - Z[i+1] * sin_theta * dx[i] / 2 - Z[i] * cos_theta * (yt_l[i] - y_b[i])) / (Z[i+1] * cos_theta))  # Equation (30)
            yt_l[i+1] = yt_r[i]
    df['yt_l'] = yt_l
    df['yt_r'] = yt_r
    
    # --- Check convergence ---
    converged = abs(FS_force - FS_moment) < tol
    if not converged:
        return False, "Spencer's method did not converge within the maximum number of iterations."
    else:
        results = {}
        results['method'] = 'spencer'
        results['FS'] = FS_force
        results['theta'] = theta_opt

        # debug print values per slice
        if debug_level >= 2:
            for i in range(len(Q)):
                print(f"Slice {i+1}: Q = {Q[i]:.1f}, y_q = {y_q[i]:.2f}, Fh = {Fh[i]:.1f}, Fv = {Fv[i]:.1f}, Mo = {Mo[i]:.2f}")

        return True, results




