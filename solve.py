
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

    N_eff = W * cos_alpha - u * dl * cos2_alpha
    numerator = c * dl + N_eff * np.tan(phi)
    denominator = W * sin_alpha - shear_reinf
    FS = numerator.sum() / denominator.sum() if denominator.sum() != 0 else float('inf')

    df['n_eff'] = N_eff  # store effective normal forces in df

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
    N_eff = W * cos_alpha - u * dl * cos2_alpha
    num = c * dl + N_eff * tan_phi
    for _ in range(max_iter):
        denom = cos_alpha + (sin_alpha * tan_phi) / F_guess
        terms = num / denom
        F_calc = terms.sum() / denominator
        if abs(F_calc - F_guess) < tol:
            converge = True
            break
        F_guess = F_calc

    df['n_eff'] = N_eff  # store effective normal forces in df

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
    c = df['c'].values
    dl = df['dl'].values
    u = df['u'].values

    sin_alpha = np.sin(alpha_rad)
    cos_alpha = np.cos(alpha_rad)

    # === Calculate base FS ===

    N_eff = W * cos_alpha - u * dl
    F = sum(c * dl + N_eff * tan_phi) / sum(W * sin_alpha)

    df['n_eff'] = N_eff  # store effective normal forces in df

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

    # === Return solution ===
    FS = F * fo

    results = {}
    results['method'] = 'janbu_corrected'
    results['FS'] = FS
    results['fo'] = fo
    return True, results


def force_equilibrium_OLD(df, theta_list, circular=True,
                      fs_min=0.5, fs_max=4.0, tol=1e-6, max_iter=50, debug=False):
    """
    Limit‐equilibrium by force equilibrium in X & Y with variable interslice angles.

    Parameters:
        df (pd.DataFrame): must contain columns
            'alpha' (slice base inclination, degrees),
            'phi'   (slice friction angle, degrees),
            'c'     (cohesion),
            'dl'    (slice base length),
            'w'     (slice weight),
            'u'     (pore force per unit length)
        theta_list (array-like): slice-boundary force inclinations (degrees),
                                 length must be n+1 if there are n slices
        circular (bool): unused but kept for API consistency
        fs_min, fs_max (float): bracket for FS search
        tol (float): convergence tolerance on residual
        max_iter (int): maximum bisection steps
        debug (bool): print residuals during iteration

    Returns:
        (bool, dict or str):
           - If converged: (True, {'method':'force_equilibrium','FS':<value>})
           - If failed:   (False, "error message")
    """
    n = len(df)
    if len(theta_list) != n+1:
        return False, f"theta_list length ({len(theta_list)}) must be n+1 ({n+1})"

    # extract and convert to radians
    alpha   = np.radians(df['alpha'].values)
    phi     = np.radians(df['phi'].values)
    c       = df['c'].values
    w       = df['w'].values
    u       = df['u'].values
    dl      = df['dl'].values
    theta   = np.radians(np.asarray(theta_list))

    def residual(FS):
        """Return the last interslice force Z[n] for a given FS."""
        c_m       = c / FS
        tan_phi_m = np.tan(phi) / FS
        Z = np.zeros(n+1)
        # Z[0] = 0 by definition (no force entering leftmost slice)
        for i in range(n):
            ca, sa = np.cos(alpha[i]), np.sin(alpha[i])
            # build coefficient matrix: [N_i, Z_{i+1}]
            A = np.array([
                [tan_phi_m[i]*ca - sa,   -np.cos(theta[i+1])],
                [tan_phi_m[i]*sa + ca,   -np.sin(theta[i+1])]
            ])
            # RHS vector
            b0 = -c_m[i]*dl[i]*ca + u[i]*dl[i]*sa - Z[i]*np.cos(theta[i])
            b1 = -c_m[i]*dl[i]*sa - u[i]*dl[i]*ca + w[i]    - Z[i]*np.sin(theta[i])
            N_i, Z_ip1 = np.linalg.solve(A, np.array([b0, b1]))
            Z[i+1] = Z_ip1
        return Z[n] # side force on right side of last slice should be zero if forces balance.

    # bracket check
    r_lo = residual(fs_min)
    r_hi = residual(fs_max)
    if debug:
        print(f"FS={fs_min:.4f} -> residual={r_lo:.4g}; FS={fs_max:.4f} -> residual={r_hi:.4g}")
    if r_lo * r_hi > 0:
        return False, (f"Cannot bracket root in [{fs_min}, {fs_max}]: "
                       f"residuals are {r_lo:.4g}, {r_hi:.4g}")

    # bisection
    lo, hi = fs_min, fs_max
    for it in range(max_iter):
        mid = 0.5*(lo + hi)
        r_mid = residual(mid)
        if debug:
            print(f" iter {it:2d}: FS={mid:.6f}, residual={r_mid:.4g}")
        if abs(r_mid) < tol:
            FS = mid
            break
        if r_lo * r_mid < 0:
            hi, r_hi = mid, r_mid
        else:
            lo, r_lo = mid, r_mid
    else:
        return False, f"force_equilibrium did not converge after {max_iter} iterations"

    return True, {'FS': FS}

from scipy.optimize import newton

def force_equilibrium(df, theta_list, fs_guess=1.5, tol=1e-6, max_iter=50, debug=False):
    """
    Limit‐equilibrium by force equilibrium in X & Y with variable interslice angles.

    Parameters:
        df (pd.DataFrame): must contain columns
            'alpha' (slice base inclination, degrees),
            'phi'   (slice friction angle, degrees),
            'c'     (cohesion),
            'dl'    (slice base length),
            'w'     (slice weight),
            'u'     (pore force per unit length)
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

    n = len(df)
    if len(theta_list) != n+1:
        return False, f"theta_list length ({len(theta_list)}) must be n+1 ({n+1})"

    # extract and convert to radians
    alpha   = np.radians(df['alpha'].values)
    phi     = np.radians(df['phi'].values)
    c       = df['c'].values
    w       = df['w'].values
    u       = df['u'].values
    dl      = df['dl'].values
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
            A = np.array([
                [tan_phi_m[i]*ca - sa,   -np.cos(theta[i+1])],
                [tan_phi_m[i]*sa + ca,   -np.sin(theta[i+1])]
            ])
            b0 = -c_m[i]*dl[i]*ca + u[i]*dl[i]*sa - Z[i]*np.cos(theta[i])
            b1 = -c_m[i]*dl[i]*sa - u[i]*dl[i]*ca + w[i]    - Z[i]*np.sin(theta[i])
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

    df['n_eff'] = N  # store effective normal forces in df
    df['z'] = Z[:-1]  # store interslice forces in df, adjust length to n slices

    if debug:
        r_opt = residual(FS_opt)
        print(f" Converged FS = {FS_opt:.6f}, residual = {r_opt:.4g}")

    return True, {'FS': FS_opt}

def corps_engineers(df, circular=True, debug=True):
    """
    Corps‐of‐Engineers style force‐equilibrium solver.

    1. Computes a single θ from the slope between
       (x_l[0], y_lt[0]) and (x_r[-1], y_rt[-1]).
    2. Builds a constant θ‐array of length n+1.
    3. Calls force_equilibrium(df, theta_array, circular).

    Parameters:
        df (pd.DataFrame): Must include at least ['x_l','y_lt','x_r','y_rt']
                           plus all the columns required by force_equilibrium:
                           ['alpha','phi','c','dl','w','u','dx'].
        circular (bool): Passed through to force_equilibrium (unused).

    Returns:
        Tuple(bool, dict or str): Whatever force_equilibrium returns.
    """
    # endpoints of the slip surface
    x0, y0 = df['x_l'].iat[0], df['y_lt'].iat[0]
    x1, y1 = df['x_r'].iat[-1], df['y_rt'].iat[-1]

    # compute positive slope‐angle
    dx = x1 - x0
    dy = y1 - y0
    if abs(dx) < 1e-12:
        theta_deg = 90.0
    else:
        theta_deg = abs(np.degrees(np.arctan2(dy, dx)))

    # one theta per slice boundary
    n = len(df)
    theta_list = np.full(n+1, theta_deg)

    df['theta'] = theta_list[:-1]  # store theta in df. Adjust length to n slices.

    # delegate to your force_equilibrium solver
    success, results = force_equilibrium(df, theta_list, debug=debug)
    if not success:
        return success, results
    else:
        results['method'] = 'corps_engineers'  # append method
        results['theta'] = theta_deg           # append theta
        return success, results

def lowe_karafiath(df, circular=True, debug=True):
    """
    Lowe-Karafiath limit equilibrium: variable interslice inclinations equal to
    the average of the top‐and bottom‐surface slopes of the two adjacent slices
    at each boundary.
    """
    n = len(df)

    # grab boundary coords
    x_l = df['x_l'].values
    y_lt = df['y_lt'].values
    y_lb = df['y_lb'].values
    x_r = df['x_r'].values
    y_rt = df['y_rt'].values
    y_rb = df['y_rb'].values

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

    df['theta'] = theta_list[:-1]  # store theta in df. Adjust length to n slices.

    # call your force_equilibrium solver
    success, results = force_equilibrium(df, theta_list, debug=debug)
    if not success:
        return success, results
    else:
        results['method'] = 'lowe_karafiath'  # append method
        return success, results

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
    dl = df['dl'].values
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

    df['theta'] = theta_opt  # store theta in df.

    # --- NEW: compute interslice forces Z and effective normals n_eff ---
    Q = compute_Q(FS_force, theta_rad)
    n = len(Q)
    Z = np.zeros(n+1)
    # Q_i = Z_i - Z_{i+1}  ⇒  Z_{i+1} = Z_i - Q_i
    for i in range(n):
        Z[i+1] = Z[i] - Q[i]

    # horizontal component X = Z * sin(θ)
    X = Z * np.sin(theta_rad)

    # --- compute effective normal N' per slice from your doc:
    # N' = [ -c_m·dl·sinα  - u·dl·cosα  + W  - X_i  + X_{i+1} ]
    #      -----------------------------------------------
    #          tan(φ_m)·sinα  +  cosα

    # mobilized cohesion and friction
    c_m      = c / FS_force
    tan_phi_m = np.tan(phi) / FS_force

    # X runs from 0..n, and X[i] is the horizontal side‐force on the left face of slice i
    X_i   = X[:-1]
    X_ip1 = X[1:]

    num = (
        -c_m * dl * np.sin(alpha)
        - u   * dl * np.cos(alpha)
        + w
        - X_i
        + X_ip1
    )
    denom = tan_phi_m * np.sin(alpha) + np.cos(alpha)

    N_eff = num / denom

    # store back into df
    df['z']     = Z[:-1]        # Z_i acting on slice i’s left face
    df['n_eff'] = N_eff

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

def compute_line_of_thrust(df, FS, debug=False):
    """
    Compute the line of thrust via dual-sweep moment equilibrium,
    with optional debug output. This should only be called after
    Spencer's method has been run to obtain the factor of safety (FS).
    It does not work for methods that do not satisfy complete force
    and moment equilibrium.

    Parameters
    ----------
    df : pd.DataFrame
        Must include columns:
            'alpha', 'phi', 'c', 'w', 'u',
            'dl', 'dx', 'x_l', 'x_r', 'y_lb', 'y_rb'
    FS : float
        Factor of safety from Spencer’s solution.
    theta_deg : float
        Side-force inclination (deg).
    debug : bool
        If True, export a per-slice debug table to Excel.

    Returns
    -------
    dict
        'sigma_eff'   : N_eff / dl per slice,
        'delta_y'     : averaged moment arm at each boundary (n+1),
        'thrust_line' : Shapely LineString of the averaged y's.
    """
    # 1) unpack & force equilibrium
    n      = len(df)
    alpha  = np.radians(df['alpha'].values)
    phi    = np.radians(df['phi'].values)
    theta  = np.zeros(n+1)
    theta[:-1]   = np.radians(df['theta'].values)
    c      = df['c'].values
    w      = df['w'].values
    u      = df['u'].values
    dl     = df['dl'].values
    dx     = df['dx'].values
    y_lb   = df['y_lb'].values
    y_rb   = df['y_rb'].values

    c_m       = c / FS
    tan_phi_m = np.tan(phi) / FS

    N_eff = df['n_eff'].values
    Z     = np.zeros(n+1)
    Z[:-1]  = df['z'].values

    tol = 1e-8

    # 1) decompose side-forces
    X = Z * np.sin(theta)
    E = Z * np.cos(theta)

    # 2a) left-to-right moment sweep about each slice's lower-right
    delta_y_L = np.zeros(n+1)  # Moment arms to side forces relative to slice lower-right.
    y_L       = np.zeros(n+1)  # Absolute y-coordinates of the thrust line on all slice boundaries based on left sweep
    delta_y_L[0] = 0.0         # First slice pivot (starting from the left)
    y_L[0]      = y_lb[0]      # First slice left side y-coordinate

    for i in range(n-1):  # Loop from 0 to n-2. Last slice = n-1 and the right-side moment arm on that slice is fixed.

        arm_left = y_L[i] - y_rb[i] #  left-side moment arm E[i] (pivot at y_rb[i])
        num = (
            E[i]*arm_left
            + X[i]*dx[i]
            - w[i]*(dx[i]/2)
            + (N_eff[i] + u[i]*dl[i])*(dl[i]/2)
        )
        delta_y_L[i+1]  = num / E[i+1] if abs(E[i+1]) > tol else 0  # right-side moment arm
        y_L[i+1]  = y_rb[i] + delta_y_L[i+1]    # absolute y value on right side

    delta_y_L[n] = 0.0         # Last (right-most) slice pivot
    y_L[n]       = y_rb[n-1]   # Last (right-most) slice right side y-coordinate

    # 2b) right-to-left moment sweep about each slice's lower-left
    delta_y_R = np.zeros(n+1)   # Moment arms to side forces relative to slice lower-left.
    y_R       = np.zeros(n+1)   # Absolute y-coordinates of the thrust line on all slice boundaries based on right sweep
    delta_y_R[n] = 0.0          # First slice pivot (starting from the right)
    y_R[n]      = y_rb[n-1]     # First slice left side y-coordinate

    for i in range(n-1, 0, -1):  # Loop from n-1 to 1. Last slice = 0 and the left-side moment arm on that slice is fixed.

        arm_right = y_R[i+1] - y_lb[i] #  right-side moment arm for E[i+1] (pivot at y_lb[i])
        num = (
            E[i+1]*arm_right
            - X[i+1]*dx[i]
            - w[i]*(dx[i]/2)
            + (N_eff[i] + u[i]*dl[i])*(dl[i]/2)
        )

        delta_y_R[i]  = num / E[i] if abs(E[i]) > tol else 0
        y_R[i]     = y_lb[i] + delta_y_R[i]

    delta_y_R[0] = 0.0       # Last (left-most) slice pivot
    y_R[0] = y_lb[0]         # Last (left-most) slice left side y-coordinate

    # 2c) average both sweeps
    y_bound = 0.5 * (y_L + y_R)

    # 3) build LineString
    x_bound = np.empty(n+1)
    x_bound[0]  = df['x_l'].iat[0]
    x_bound[1:] = df['x_r'].values
    thrust_line = LineString(np.column_stack([x_bound, y_bound]))

    # 4) debug export
    if debug:

        rows = []
        for i in range(n):

            # recompute arm_left and arm_right here for debug
            arm_left  = y_L[i]   - y_rb[i]
            arm_right = y_R[i+1] - y_lb[i]

            if i==n-1:
                m_res_L = 0
            else:
                m_res_L = (
                    - E[i]*arm_left
                    - X[i]*dx[i]
                    + E[i + 1] * delta_y_L[i + 1]
                    + w[i]*(dx[i]/2)
                    - (N_eff[i] + u[i]*dl[i])*(dl[i]/2)
                )
            if i==0:
                m_res_R = 0
            else:
                m_res_R = (
                    - E[i] * delta_y_R[i]
                    + E[i+1]*arm_right
                    - X[i+1]*dx[i]
                    - w[i]*(dx[i]/2)
                    + (N_eff[i] + u[i]*dl[i])*(dl[i]/2)
                )

            rows.append({
                'slice':           i,
                'N_eff':           N_eff[i],
                'dl':              dl[i],
                'dx':              dx[i],
                'y_lb':            y_lb[i],
                'y_rb':            y_rb[i],
                'u':               u[i],
                'w':               w[i],
                'c_m':             c_m[i],
                'tan_phi_m':       tan_phi_m[i],
                'Z_left':          Z[i],
                'Z_right':         Z[i+1],
                'X_left':          X[i],
                'X_right':         X[i+1],
                'E_left':          E[i],
                'E_right':         E[i+1],
                'delta_y_L': delta_y_L[i],
                'delta_y_R': delta_y_R[i],
                'y_L_abs_left':    y_L[i],
                'y_L_abs_right':   y_L[i+1],
                'y_R_abs_left':    y_R[i],
                'y_R_abs_right':   y_R[i+1],
                'y_ave_left':      y_bound[i],
                'y_ave_right':     y_bound[i+1],
                'mom_res_L':       m_res_L,
                'mom_res_R':       m_res_R
            })

        df_debug = pd.DataFrame(rows)
        df_debug.to_excel("thrust_calc_results.xlsx", index=False)
        print("Debug table written to thrust_calc_results.xlsx")

    # 6) return results
    return thrust_line

