
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


def extract_spencer_Q(df, FS, theta_deg, debug=False):
    """
    Recomputes Q_i values from Spencer method after FS and theta are known,
    and enforces consistent sign convention based on the dominant direction
    (per Spencer's assumption of parallel and consistently oriented interslice forces).

    Parameters:
        df (pd.DataFrame): Must contain 'alpha', 'phi', 'c', 'dx', 'w', 'u'
        FS (float): Factor of safety
        theta_deg (float): Spencer interslice force inclination in degrees

    Returns:
        np.ndarray: Q values per slice (adjusted to have consistent direction)
    """
    import numpy as np

    theta_rad = np.radians(theta_deg)
    alpha = np.radians(df['alpha'].values)
    phi = np.radians(df['phi'].values)
    c = df['c'].values
    dx = df['dx'].values
    w = df['w'].values
    u = df['u'].values

    sec_alpha = 1 / np.cos(alpha)
    theta_diff = alpha - theta_rad

    term1 = w * np.sin(alpha)
    term2 = (c / FS) * dx * sec_alpha
    term3 = (w * np.cos(alpha) - u * dx * sec_alpha) * (np.tan(phi) / FS)

    numerator = term1 - term2 - term3
    denominator = np.cos(theta_diff) * (1 + (np.tan(theta_diff) * np.tan(phi)) / FS)

    Q = numerator / denominator

    # for debugging: print values per slice
    print('POST Spencer Q values:')
    for i in range(len(Q)):
        print(f"Slice {i}: Q = {Q[i]:.3f}, alpha = {degrees(alpha[i]):.2f}, phi = {degrees(phi[i]):.2f}, c = {c[i]:.2f}, w = {w[i]:.2f}, u = {u[i]:.2f}")

    return Q


def compute_line_of_thrust_spencer(df, FS, theta_deg, debug=False):
    """
    Computes the line of thrust for Spencer's method of slices by:
      1. Extracting Q_i side-force resultants from the Spencer solution.
      2. Summing moments about each slice center to locate thrust.
      3. Computing the effective normal force N' on each slice base.

    Parameters:
        df (pd.DataFrame): Must contain columns 'alpha', 'dx', 'w', 'u', 'y_cb'.
        FS (float): Factor of safety from Spencer's method.
        theta_deg (float): Interslice force inclination (degrees).
        debug (bool): If True, export detailed debug DataFrame to Excel.

    Returns:
        y_thrust (np.ndarray): Thrust-line elevations at each slice center (length = n).
        N_prime (np.ndarray): Effective normal force on each slice base (length = n).
    """
    # --- 1) Extract Q_i from Spencer's solution ---
    Q = extract_spencer_Q(df, FS, theta_deg, debug=debug)  # length n
    n = len(Q)

    # --- 2) Solve for Z_i (parallel interslice forces) ---
    Z = np.zeros(n+1)
    for i in range(n):
        Z[i+1] = Z[i] - Q[i]
    Z[-1] = 0.0  # enforce boundary condition

    # --- 3) Decompose Z into components along/normal to base ---
    theta = radians(theta_deg)
    X = Z * sin(theta)   # vertical side force component
    E = Z * cos(theta)   # horizontal side force component

    # --- 4) Moment equilibrium sweep about slice centers ---
    dx = df['dx'].values       # slice base widths
    y_lb = df['y_lb'].values  # y at left bottom of base
    y_rb = df['y_rb'].values  # y at right bottom of base
    y_cb = df['y_cb'].values  # y at slice center
    tol = 1e-8

    # initialize moment arm above base center
    dy = np.zeros(n+1)
    dy[0] = 0.0  # no moment arm at left boundary

    # For each slice, enforce:
    #   -E[i]*dy[i] - X[i]*(dx[i]/2) + E[i+1]*dy[i+1] - X[i+1]*(dx[i]/2) = 0
    # => dy[i+1] = (E[i]*dy[i] + X[i]*(dx[i]/2) + X[i+1]*(dx[i]/2)) / E[i+1]
    y_thrust = np.zeros(n+1)
    y_thrust[0] = y_lb[0]
    for i in range(n):
        denom = E[i+1]
        if abs(denom) > tol:
            dy[i+1] = (E[i]*dy[i] + X[i]*(dx[i]/2) + X[i+1]*(dx[i]/2)) / denom
        else:
            dy[i+1] = 0.0
        y_thrust[i+1] = y_cb[i] + dy[i+1]

    # thrust elevations at slice edges (one more than number of slices)
    x_edges = list(df['x_l'].values) + [df['x_r'].values[-1]]
    from shapely.geometry import LineString
    line_of_thrust = LineString(zip(x_edges, y_thrust))

    # --- 5) Compute effective normal force N' per slice ---
    alpha = np.radians(df['alpha'].values)
    sin_a = np.sin(alpha)
    cos_a = np.cos(alpha)
    c_m = df['c'].values/FS
    phi_rad = np.radians(df['phi'].values)
    tan_phi_m = np.tan(phi_rad) / FS
    dl = df['dl'].values
    w = df['w'].values
    u = df['u'].values
    N_prime = np.zeros(n)
    sigma_prime = np.zeros(n)
    num = np.zeros(n)
    denom = np.zeros(n)
    for i in range(n):
        num[i] = - c_m[i] * dl[i] * sin_a[i] - u[i] * dl[i] * cos_a[i] + w[i] - X[i] + X[i+1]
        denom[i] = tan_phi_m[i] * sin_a[i] + cos_a[i]
        N_prime[i] = num[i] / denom[i] if abs(denom[i]) > tol else 0.0
        sigma_prime[i] = N_prime[i] / dl[i]

    # --- Debug export if needed ---
    if debug:
        import pandas as pd
        rows = []
        for i in range(n):
            rows.append({
                'slice': i,
                'N_prime': N_prime[i],
                'num': num[i],
                'denom': denom[i],
                'sigma_prime': sigma_prime[i],
                'dl': dl[i],
                'dx': dx[i],
                'alpha_rad': alpha[i],
                'sin_a': sin_a[i],
                'cos_a': cos_a[i],
                'y_cb': y_cb[i],
                'y_lb': y_lb[i],
                'y_rb': y_rb[i],
                'u': u[i],
                'w': w[i],
                'c_m': c_m[i],
                'tan_phi_m': tan_phi_m[i],
                'Q': Q[i],
                'Z_left': Z[i],
                'Z_right': Z[i+1],
                'X_left': X[i],
                'X_right': X[i + 1],
                'E_left': E[i],
                'E_right': E[i+1],
                'dy_left': dy[i],
                'dy_right': dy[i+1],
                'y_left': y_thrust[i],
                'y_right': y_thrust[i+1],
            })
        df_debug = pd.DataFrame(rows)
        df_debug.to_excel("thrust_calc_results.xlsx", index=False)
        print("Debug table written to thrust_calc_results.xlsx")

    return line_of_thrust, sigma_prime
