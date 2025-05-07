
import numpy as np
import pandas as pd
from shapely.geometry import LineString, Point
from math import sin, cos, tan, radians, atan, atan2, degrees
from scipy.optimize import minimize_scalar, root_scalar
from tabulate import tabulate


def oms(df, circular=True):
    """
    Computes the Factor of Safety (FS) using the Ordinary Method of Slices (OMS).
    This method works on circular failure surfaces only.

    Parameters:
        df (pd.DataFrame): DataFrame containing slice data with columns:
            - 'alpha': base angle (degrees)
            - 'w': total slice weight
            - 'c': cohesion
            - 'phi': friction angle (degrees)
            - 'u': pore pressure
            - 'dl': base length
            - 'shear_reinf' (optional): reinforcement force opposing sliding (FL)
            - 'normal_reinf' (optional): reinforcement force contributing to base normal force (FT)

    Returns:
        dict:
            - float: Computed Factor of Safety (FS)
            - np.ndarray: Normal force on the base of each slice
    """

    if not circular:
        return False, 'OMS method is only applicable for circular failure surfaces.'

    alpha_rad = np.radians(df['alpha'])

    cos_alpha = np.cos(alpha_rad)
    sin_alpha = np.sin(alpha_rad)
    cos2_alpha = cos_alpha ** 2

    W = df['w'].values
    shear_reinf = df.get('shear_reinf', 0).values
    c = df['c'].values
    phi = np.radians(df['phi']).values
    u = df['u'].values
    dl = df['dl'].values

    N = W * cos_alpha - u * dl * cos2_alpha
    numerator = c * dl + N * np.tan(phi)
    denominator = W * sin_alpha - shear_reinf
    FS = numerator.sum() / denominator.sum() if denominator.sum() != 0 else float('inf')

    results = {}
    results['method'] = 'oms'
    results['FS'] = FS
    return True, results

def bishop(df, circular=True, tol=1e-6, max_iter=100):
    """
    Computes the Factor of Safety (FS) using Bishop's Simplified Method.
    This method works on circular failure surfaces only.
    It iterates on the factor of safety until convergence is achieved.

    Parameters:
        df (pd.DataFrame): DataFrame containing slice data with necessary columns.
        tol (float, optional): Convergence tolerance. Default is 1e-6.
        max_iter (int, optional): Maximum number of iterations. Default is 100.

    Returns:
        tuple:
            - float: Computed Factor of Safety (FS)
            - np.ndarray: Normal force on the base of each slice
            - bool: Whether the solution converged
    """

    if not circular:
        return False, 'Bishop method is only applicable for circular failure surfaces.'

    alpha_rad = np.radians(df['alpha'])
    cos_alpha = np.cos(alpha_rad)
    sin_alpha = np.sin(alpha_rad)
    cos2_alpha = cos_alpha**2
    tan_phi = np.tan(np.radians(df['phi']))

    W = df['w']
    shear_reinf = df.get('shear_reinf', 0)
    c = df['c']
    dl = df['dl']
    u = df['u']

    # Right-hand side: sum of W * sin(alpha)
    denominator = (W * sin_alpha - shear_reinf).sum()

    # Start iteration with an initial guess
    converge = False
    F_guess = 1.0
    N = W * cos_alpha - u * dl * cos2_alpha
    num = c * dl + N * tan_phi
    for _ in range(max_iter):
        denom = cos_alpha + (sin_alpha * tan_phi) / F_guess
        terms = num / denom
        F_calc = terms.sum() / denominator
        if abs(F_calc - F_guess) < tol:
            converge = True
            break
        F_guess = F_calc


    if not converge:
        return False, 'Bishop method did not converge within the maximum number of iterations.'
    else:
        results = {}
        results['method'] = 'bishop'
        results['FS'] = F_calc
        return True, results


def janbu_corrected(df, circular=True, tol=1e-6, max_iter=100):
    """
    Computes the Factor of Safety (FS) using Janbu's Simplified Method with correction factor.
    Applies the Janbu correction based on d/L ratio and soil type.

    Parameters:
        df (pd.DataFrame): DataFrame containing slice data.
        tol (float, optional): Convergence tolerance. Default is 1e-6.
        max_iter (int, optional): Maximum number of iterations. Default is 100.

    Returns:
        float: Corrected Factor of Safety (FS)
        float: Correction factor (fo)
        bool: Convergence flag
    """
    alpha_rad = np.radians(df['alpha'])
    phi_rad = np.radians(df['phi'])
    tan_phi = np.tan(phi_rad)

    W = df['w']
    shear_reinf = df.get('shear_reinf', 0)
    normal_reinf = df.get('normal_reinf', 0)
    c = df['c'].values
    dl = df['dl'].values
    u = df['u'].values

    sin_alpha = np.sin(alpha_rad)
    cos_alpha = np.cos(alpha_rad)

    # === Compute Janbu correction factor ===
    x_l = df['x_l'].iloc[0]
    y_lt = df['y_lt'].iloc[0]
    x_r = df['x_r'].iloc[-1]
    y_rt = df['y_rt'].iloc[-1]

    L = np.hypot(x_r - x_l, y_rt - y_lt)
    x0 = df['x_c'].values
    y0 = df['y_cb'].values
    numerator = np.abs((y_rt - y_lt) * x0 - (x_r - x_l) * y0 + x_r * y_lt - y_rt * x_l)
    dists = numerator / L
    d = np.max(dists)
    dL_ratio = d / L

    phi_sum = df['phi'].sum()
    c_sum = df['c'].sum()
    if phi_sum == 0:
        b1 = 0.69
    elif c_sum == 0:
        b1 = 0.31
    else:
        b1 = 0.50

    fo = 1 + b1 * (dL_ratio - 1.4 * dL_ratio**2)

    # === Iterative FS solver ===
    F = 1.0  # initial guess
    converged = False
    for _ in range(max_iter):
        N = W * cos_alpha + normal_reinf
        S = c * dl + (N - u * dl) * tan_phi / F
        T = W * sin_alpha - shear_reinf

        F_new = S.sum() / T.sum()

        if abs(F_new - F) < tol:
            converged = True
            break

        F = F_new

    # === Return solution ===
    FS = F * fo

    if not converged:
        return False, 'Janbu-Corrected method did not converge within the maximum number of iterations.'
    else:
        results = {}
        results['method'] = 'janbu_corrected'
        results['FS'] = FS
        results['fo'] = fo
        return True, results

    return results

def spencer(df, circular=True):
    """
    Spencer's Method using Steve G. Wright's formulation.
    Solves for FS_force and FS_moment independently using the Wright Q equation.

    Parameters:
        df (pd.DataFrame): Must include:
            'alpha', 'phi', 'c', 'w', 'u', 'dl', 'x_c', 'y_cb'

    Returns:
        float: FS where FS_force = FS_moment
        float: beta (degrees)
        bool: converged flag
    """

    beta_bounds = (-60, 60)
    tol = 1e-6
    max_iter = 100

    alpha = np.radians(df['alpha'].values)
    phi = np.radians(df['phi'].values)
    c = df['c'].values
    dx = df['dx'].values
    w = df['w'].values
    u = df['u'].values
    x_c = df['x_c'].values
    y_cb = df['y_cb'].values

    R = 120  # For circular case

    def compute_Q(F, theta_rad):
        theta_diff = alpha - theta_rad
        sec_alpha = 1 / np.cos(alpha)
        term1 = w * np.sin(alpha)
        term2 = (c / F) * dx * sec_alpha
        term3 = (w * np.cos(alpha) - u * dx * sec_alpha) * (np.tan(phi) / F)
        numerator = term1 - term2 - term3
        denominator = np.cos(theta_diff) * (1 + (np.tan(theta_diff) * np.tan(phi)) / F)
        Q = numerator / denominator
        return Q

    fs_min = 0.01
    fs_max = 20.0

    def fs_force(theta_rad):
        def residual(F):
            Q = compute_Q(F, theta_rad)
            return Q.sum()
        result = minimize_scalar(lambda F: abs(residual(F)), bounds=(fs_min, fs_max), method='bounded', options={'xatol': tol})
        return result.x

    def fs_moment(theta_rad):
        def residual(F):
            Q = compute_Q(F, theta_rad)
            if circular:
                theta_diff = alpha - theta_rad
                return np.sum(Q * np.cos(theta_diff))
            else:
                return np.sum(-Q * x_c * np.sin(theta_rad) + Q * y_cb * np.cos(theta_rad))
        result = minimize_scalar(lambda F: abs(residual(F)), bounds=(fs_min, fs_max), method='bounded', options={'xatol': tol})
        return result.x

    def fs_difference(theta_deg):
        theta_rad = radians(theta_deg)
        Ff = fs_force(theta_rad)
        Fm = fs_moment(theta_rad)
        return abs(Ff - Fm)

    result = minimize_scalar(fs_difference, bounds=beta_bounds, method='bounded', options={'xatol': tol})
    theta_opt = result.x
    theta_rad = radians(theta_opt)
    FS_force = fs_force(theta_rad)
    FS_moment = fs_moment(theta_rad)

    # Check convergence
    converged = abs(FS_force - FS_moment) < tol
    if not converged:
        return False, "Spencer's method did not converge within the maximum number of iterations."
    else:
        results = {}
        results['method'] = 'spencer'
        results['FS'] = FS_force
        results['theta'] = theta_opt

        # debug print values per slice
        Q = compute_Q(FS_force, theta_rad)
        for i in range(len(Q)):
            print(f"Slice {i}: Q = {Q[i]:.3f}, alpha = {degrees(alpha[i]):.2f}, phi = {degrees(phi[i]):.2f}, c = {c[i]:.2f}, w = {w[i]:.2f}, u = {u[i]:.2f}")

        return True, results


import numpy as np


import numpy as np

import numpy as np
from shapely.geometry import LineString

import numpy as np
from shapely.geometry import LineString

import numpy as np

def compute_line_of_thrust(df, FS, theta_deg, debug=False):
    """
    Compute the line of thrust following the Slope Tools methodology,
    with optional debug output.

    Parameters
    ----------
    df : pd.DataFrame
        Must include columns:
            'alpha' (deg), 'phi' (deg), 'c', 'w', 'u',
            'dl', 'dx', 'x_l', 'x_r', 'y_cb'
    FS : float
        Factor of safety from Spencer’s solution.
    theta_deg : float
        Common side‐force angle (deg) from Spencer’s solution.
    debug : bool, default False
        If True, assemble and export a per‐slice Excel debug table.

    Returns
    -------
    dict
        'N', 'Z', 'X', 'E', 'delta_y', 'thrust_line'
    """
    # unpack
    n       = len(df)
    alpha   = np.radians(df['alpha'].values)
    phi     = np.radians(df['phi'].values)
    c       = df['c'].values
    w       = df['w'].values
    u       = df['u'].values
    dl      = df['dl'].values
    dx      = df['dx'].values
    theta   = np.radians(theta_deg)

    # mobilized properties
    c_m       = c / FS
    tan_phi_m = np.tan(phi) / FS

    # allocate
    N       = np.zeros(n)
    Z       = np.zeros(n + 1)  # side‐force magnitudes
    Z[0]    = 0.0

    # 1) sweep left→right for N[i], Z[i+1]
    for i in range(n):
        A = np.array([
            [ tan_phi_m[i]*np.cos(alpha[i]) - np.sin(alpha[i]),  -np.cos(theta) ],
            [ tan_phi_m[i]*np.sin(alpha[i]) + np.cos(alpha[i]),  -np.sin(theta) ]
        ])
        b = np.array([
            -c_m[i]*np.cos(alpha[i]) + u[i]*dl[i]*tan_phi_m[i]*np.cos(alpha[i]) - Z[i]*np.cos(theta),
            -c_m[i]*np.sin(alpha[i]) + u[i]*dl[i]*tan_phi_m[i]*np.sin(alpha[i]) + w[i]   - Z[i]*np.sin(theta)
        ])
        sol       = np.linalg.solve(A, b)
        N[i]      = sol[0]
        Z[i+1]    = sol[1]

    # 2) decompose side‐forces
    X = Z * np.cos(theta)
    E = Z * np.sin(theta)

    # 3) moment‐arm sweep for Δy
    delta_y = np.zeros(n+1)
    for i in range(n):
        num = E[i]*delta_y[i] + (X[i] + X[i+1])*(dx[i]/2)
        delta_y[i+1] = num / E[i+1] if E[i+1] != 0 else np.nan

    # 4) build a Shapely LineString
    from shapely.geometry import LineString
    x_bound = np.empty(n+1)
    x_bound[0]  = df['x_l'].iloc[0]
    x_bound[1:] = df['x_r'].values

    ycb      = df['y_cb'].values
    y_bound  = np.array([ ycb[min(i, n-1)] + delta_y[i] for i in range(n+1) ])
    thrust_line = LineString(np.column_stack([x_bound, y_bound]))

    # --- debug table & Excel export ---
    if debug:
        import pandas as pd
        rows = []
        for i in range(n):
            # recompute A, b for residuals
            A = np.array([
                [ tan_phi_m[i]*np.cos(alpha[i]) - np.sin(alpha[i]),  -np.cos(theta) ],
                [ tan_phi_m[i]*np.sin(alpha[i]) + np.cos(alpha[i]),  -np.sin(theta) ]
            ])
            b = np.array([
                -c_m[i]*np.cos(alpha[i]) + u[i]*dl[i]*tan_phi_m[i]*np.cos(alpha[i]) - Z[i]*np.cos(theta),
                -c_m[i]*np.sin(alpha[i]) + u[i]*dl[i]*tan_phi_m[i]*np.sin(alpha[i]) + w[i]   - Z[i]*np.sin(theta)
            ])
            # force‐balance residuals
            fx_res = A[0,0]*N[i]     + A[0,1]*Z[i+1] - b[0]
            fy_res = A[1,0]*N[i]     + A[1,1]*Z[i+1] - b[1]
            # moment‐balance residual
            m_res  = E[i+1]*delta_y[i+1] - (E[i]*delta_y[i] + (X[i]+X[i+1])*(dx[i]/2))

            rows.append({
                'slice':           i,
                'alpha (deg)':     np.degrees(alpha[i]),
                "N'":              N[i],     # effective normal force
                'dl':              dl[i],    # base length
                'dx':              dx[i],    # slice width
                'u':               u[i],      # pore pressure
                'w':               w[i],        # weight
                'S':               c_m[i] + N[i]*tan_phi_m[i],  # mobilized shear strength
                'Z_left':          Z[i],
                'Z_right':         Z[i+1],
                'X_left':          X[i],
                'X_right':         X[i+1],
                'E_left':          E[i],
                'E_right':         E[i+1],
                'delta_y_left':    delta_y[i],
                'delta_y_right':   delta_y[i+1],
                'fx_residual':     fx_res,
                'fy_residual':     fy_res,
                'moment_residual': m_res
            })

        df_debug = pd.DataFrame(rows)
        df_debug.to_excel("thrust_calc_resuls.xlsx", index=False)
        print("Debug table written to thrust_calc_resuls.xlsx")

    return {
        'N': N,
        'Z': Z,
        'X': X,
        'E': E,
        'delta_y': delta_y,
        'thrust_line': thrust_line
    }