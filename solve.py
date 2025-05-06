
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



def sweep_thrust_left(df, Q, theta_deg, debug=False):
    """
    Line of thrust calculation using Spencer's method with left-to-right sweep.
    Thrust line is calculated as vertical distance above the base of each slice.
    """

    theta = np.radians(theta_deg)
    x_c = df['x_c'].values   # x-coordinates of slice centers
    x_l = df['x_l'].values   # x-coordinates of slice left edges
    x_r = df['x_r'].values   # x-coordinates of slice right edges
    y_cb = df['y_cb'].values # y-coordinates of slice centers at bottom
    y_lb = df['y_lb'].values # y-coordinates of slice left edges at bottom
    y_lt = df['y_lt'].values # y-coordinates of slice left edges at top
    y_rb = df['y_rb'].values # y-coordinates of slice right edges at bottom
    y_rt = df['y_rt'].values # y-coordinates of slice right edges at top

    w = df['w'].values        # slice weights

    n = len(Q)                   # number of slices
    y_thrust = [None] * (n + 1)  # For n slices, we need n+1 y-thrust points

    Z_left_x = 0.0
    Z_left_y = 0.0
    y_thrust[0] = (y_lt[0] - y_lb[0]) / 3.0  # start thrust at 1/3 height of leftmost slice (in case it is tension crack)

    # if debug:
    #     print(f"{'Slice':>5} | {'Q':>8} | {'Qx':>8} | {'Qy':>8} | {'ZLx':>8} | {'ZLy':>8} | {'ZRx':>8} | {'ZRy':>8} | {'dx_L':>8} | {'dy_L':>8} | {'dx_R':>8} | {'dy_R':>8} | {'y_L':>10} | {'y_R':>10}")
    #     print("-" * 150)

    debug_rows = []

    for i in range(n):
        if i == n - 1: # last slice, assign y_thrust to 1/3 height of rightmost slice boundary and break
            y_thrust[i + 1] = (y_rt[i] - y_rb[i]) / 3.0
            break

        xc = x_c[i]
        xl = x_l[i]
        xr = x_r[i]
        ycb = y_cb[i]
        ylt = y_lt[i]
        ylb = y_lb[i]
        yrt = y_rt[i]
        yrb = y_rb[i]

        Qx = Q[i] * np.cos(theta)
        Qy = Q[i] * np.sin(theta)

        Z_right_x = Qx - Z_left_x
        Z_right_y = Qy - Z_left_y


        # For all slices, to find the moment arm (line of thrust) on the right side of slice,
        # sum moments around the right corner of the slice.
        # Moment equation: (CCW rotation = positive)
        # M = Z_left_y*(xr-xl) - Z_left_x*(y_thrust_left+(ylb-yrb) + Z_right_y*0 + Z_right_x*(dy_right) + Qx*(ycb-yrb) - Qy*(xr-xc) = 0
        # Rearranging gives:
        # dy_right = (- Z_left_y*(xr-xl) + Z_left_x*(y_thrust_left+(ylb-yrb)) - Z_right_y*0 - Qx*(ycb-yrb) + Qy*(xr-xc)) / Z_right_x

        y_thrust_left = y_thrust[i]
        if abs(Z_right_x) > 1e-6:
            dy_right = (- Z_left_y*(xr-xl) + Z_left_x*(y_thrust_left+(ylb-yrb)) - Qx*(ycb-yrb) + Qy*(xr-xc)) / Z_right_x
        else:
            dy_right = 0.0
        y_thrust[i + 1] = dy_right

        if debug:
            debug_rows.append({
                'Slice': i,
                'x_c': xc,
                'x_l': xl,
                'x_r': xr,
                'y_cb': ycb,
                'y_lt': ylt,
                'y_lb': ylb,
                'y_rt': yrt,
                'y_rb': yrb,
                'Q': Q[i],
                'Qx': Qx,
                'Qy': Qy,
                'Z_left_x': Z_left_x,
                'Z_left_y': Z_left_y,
                'Z_right_x': Z_right_x,
                'Z_right_y': Z_right_y,
                'dy_right': dy_right,
                'y_thrust[i]': y_thrust[i],
                'y_thrust[i+1]': y_thrust[i + 1],
            })

        Z_left_x = -Z_right_x
        Z_left_y = -Z_right_y

        if debug:
            debug_df = pd.DataFrame(debug_rows)
            debug_df.to_excel('sweep_left_debug.xlsx', index=False)

    return y_thrust

def sweep_thrust_right(df, Q, theta_deg):
    """
    Performs a right-to-left sweep to compute the line of thrust based on interslice force transmission
    using Spencer's assumptions. Assumes all interslice forces are parallel at angle theta.

    Parameters:
        df (pd.DataFrame): Must include 'x_c', 'x_l', 'x_r', 'y_cb', 'y_rb', 'w'
        Q (np.ndarray): Net interslice force resultants for each slice (length n)
        theta_deg (float): Spencer's interslice force inclination in degrees

    Returns:
        y_thrust (list): y-locations of interslice forces (length = n+1)
    """
    import numpy as np

    theta = np.radians(theta_deg)
    x_c = df['x_c'].values
    x_l = df['x_l'].values
    x_r = df['x_r'].values
    y_cb = df['y_cb'].values
    y_rb = df['y_rb'].values
    w = df['w'].values

    n = len(Q)
    y_thrust = [None] * (n + 1)

    # Initialize: no force entering from the right
    Z_right_x = 0.0
    Z_right_y = 0.0
    y_thrust[-1] = y_rb[-1]  # start thrust at right base elevation

    for i in reversed(range(n)):

        if i == 0:
            y_thrust[0] = y_cb[0]  # anchor leftmost point to base center
            break

        Wi = w[i]
        yb = y_cb[i]
        xc = x_c[i]
        xl = x_l[i]
        xr = x_r[i]

        # Spencer Q force components
        Qx = Q[i] * np.cos(theta)
        Qy = Q[i] * np.sin(theta)

        # Left-side force = Q - right
        Z_left_x = Qx - Z_right_x
        Z_left_y = Qy - Z_right_y

        # Moment arm for Z_right
        dy_right = y_thrust[i + 1] - yb if y_thrust[i + 1] is not None else 0.0
        dx_right = xr - xc

        # Moment from right-side interslice force
        M_right = Z_right_x * dy_right - Z_right_y * dx_right

        # Moment from weight = 0 (acts through x_c)
        M_w = 0.0

        # Moment arm for Z_left_y
        dx_left = xc - xl

        # Compute vertical offset to balance moment
        if abs(Z_left_x) > 1e-6:
            dy_left = (M_right + M_w + Z_left_y * dx_left) / Z_left_x
        else:
            dy_left = 0.0

        y_thrust[i] = yb + dy_left

        # Flip direction to pass to next slice
        Z_right_x = -Z_left_x
        Z_right_y = -Z_left_y

    return y_thrust



def compute_line_of_thrust_sweep(df, Q, theta_deg):
    """
    Wrapper function that runs both left-to-right and right-to-left thrust sweeps,
    and returns all results for inspection and comparison.

    Parameters:
        df (pd.DataFrame): Must include 'x_c', 'y_cb', 'y_cg', 'w'
        Q (np.ndarray): Spencer interslice force resultants (length = n)
        theta_deg (float): Interslice force inclination angle in degrees

    Returns:
        tuple: (y_left, y_right, y_avg)
            y_left: thrust line from left-to-right sweep (length = n+1)
            y_right: thrust line from right-to-left sweep (length = n+1)
            y_avg: element-wise average of y_left and y_right
    """

    y_left = sweep_thrust_left(df, Q, theta_deg, debug=True)
    y_right = sweep_thrust_right(df, Q, theta_deg)
    y_avg = [
        (yl + yr) / 2 if yl is not None and yr is not None else None
        for yl, yr in zip(y_left, y_right)
    ]
    return y_left, y_right, y_avg


### THIS IS THE DEEPSEEK SOLUTION ###

import numpy as np
import pandas as pd
from shapely.geometry import LineString
from math import sin, cos, tan, radians, degrees
import matplotlib.pyplot as plt  # For optional plotting

import numpy as np
import pandas as pd
from shapely.geometry import LineString
from math import sin, cos, tan, radians, degrees

import numpy as np
import pandas as pd
from shapely.geometry import LineString
from math import sin, cos, tan, radians, degrees


def compute_line_of_thrust(df, spencer_results, debug=False, excel_path=None):
    """
    Computes the line of thrust using Spencer's method with proper force equilibrium.

    Args:
        df: DataFrame containing slice data (alpha, phi, c, w, u, dl, x_c, y_cb, y_ct, x_l, x_r)
        spencer_results: Dictionary with 'FS' (factor of safety) and 'theta' (interslice angle in degrees)
        debug: If True, returns detailed debug information
        excel_path: Optional path to save debug Excel file

    Returns:
        tuple: (success, result_dict)
            success: Boolean indicating if calculation succeeded
            result_dict: Dictionary containing:
                - 'thrust_line': LineString geometry
                - 'normal_forces': Array of normal forces
                - 'interslice_forces': Array of interslice forces
                - 'residuals': Dictionary of max residuals
                - 'debug_df': DataFrame with detailed calculations (if debug=True)
    """
    try:

        # Initialize variables
        FS = float(spencer_results['FS'])
        theta = radians(float(spencer_results['theta']))
        n_slices = len(df)

        Z = np.zeros(n_slices + 1)  # Interslice forces (Z[0] and Z[-1] are boundaries)
        N = np.zeros(n_slices)  # Normal forces
        y_thrust = np.zeros(n_slices + 1)  # Thrust line y-coordinates

        # Boundary conditions
        Z[0] = 0.0  # No force on left boundary
        Z[-1] = 0.0  # No force on right boundary
        y_thrust[0] = float(df.iloc[0]['y_ct'])  # Left boundary starts at ground surface
        y_thrust[-1] = float(df.iloc[-1]['y_ct'])  # Right boundary ends at ground surface

        # Precompute trigonometric terms
        alpha = np.radians(df['alpha'].values)
        phi = np.radians(df['phi'].values)
        cos_a = np.cos(alpha)
        sin_a = np.sin(alpha)
        tan_phi = np.tan(phi)
        sin_theta_a = np.sin(theta - alpha)
        cos_theta_a = np.cos(theta - alpha)

        # Debug storage
        debug_data = []

        # Forward pass: Solve for interslice forces (Z) and normal forces (N)
        for i in range(n_slices):
            W = df.iloc[i]['w']
            U = df.iloc[i]['u'] * df.iloc[i]['dl']
            c_dl = df.iloc[i]['c'] * df.iloc[i]['dl']

            # Build and solve the 2x2 system for N and Z[i+1]
            A = np.array([
                [tan_phi[i] / FS, -sin_theta_a[i]],  # F_parallel equation
                [1.0, -cos_theta_a[i]]  # F_perpendicular equation
            ])
            B = np.array([
                W * sin_a[i] - Z[i] * sin_theta_a[i] - c_dl / FS,  # F_parallel
                W * cos_a[i] - U - Z[i] * cos_theta_a[i]  # F_perpendicular
            ])

            try:
                N[i], Z[i + 1] = np.linalg.solve(A, B)
            except np.linalg.LinAlgError:
                return False, f"Singular matrix in slice {i + 1}, cannot solve equilibrium"

            if debug:
                # PROPER force equilibrium verification in x-y coordinates
                S = (c_dl + N[i] * tan_phi[i]) / FS
                S_x = S * cos_a[i]  # Shear force x-component
                S_y = S * sin_a[i]  # Shear force y-component
                N_x = N[i] * sin_a[i]  # Normal force x-component
                N_y = N[i] * cos_a[i]  # Normal force y-component

                # Interslice force components
                E_left = Z[i] * cos(theta)
                X_left = Z[i] * sin(theta)
                E_right = -Z[i + 1] * cos(theta)  # Negative because acting on opposite face
                X_right = -Z[i + 1] * sin(theta)

                # Sum of forces in x-direction
                Fx_residual = E_left + E_right + N_x + S_x

                # Sum of forces in y-direction
                Fy_residual = X_left + X_right + N_y + S_y - W[i] + U * sin_a[i]  # U acts normal to base

                debug_data.append({
                    'Slice': i + 1,
                    'N': N[i],
                    'Z_left': Z[i],
                    'Z_right': Z[i + 1],
                    'S': S,
                    'Fx_Residual': Fx_residual,
                    'Fy_Residual': Fy_residual,
                    'Z_ratio': Z[i] / Z[i + 1] if Z[i + 1] != 0 else np.inf
                })

        # Backward pass: Compute thrust line positions
        for i in range(n_slices - 1, -1, -1):
            y_base = df.iloc[i]['y_cb']
            x_center = df.iloc[i]['x_c']

            # Moment arms (relative to base center)
            arm_E_right = y_thrust[i + 1] - y_base
            arm_X_right = df.iloc[i]['x_r'] - x_center
            arm_X_left = df.iloc[i]['x_l'] - x_center

            # Moment equilibrium equation (sum M = 0)
            numerator = (Z[i + 1] * cos(theta) * arm_E_right -
                         Z[i + 1] * sin(theta) * arm_X_right +
                         Z[i] * sin(theta) * arm_X_left)

            if abs(Z[i] * cos(theta)) > 1e-10:
                y_thrust[i] = y_base + numerator / (Z[i] * cos(theta))
            else:
                y_thrust[i] = y_base + (df.iloc[i]['y_ct'] - y_base) / 2

        # Build thrust line coordinates
        thrust_points = []

        # Left boundary
        thrust_points.append((df.iloc[0]['x_l'], y_thrust[0]))

        # Points between slices
        for i in range(n_slices):
            if i < n_slices - 1:
                x_pos = (df.iloc[i]['x_r'] + df.iloc[i + 1]['x_l']) / 2
            else:
                x_pos = df.iloc[i]['x_r']  # Right boundary
            thrust_points.append((x_pos, y_thrust[i + 1]))

        # Right boundary
        thrust_points.append((df.iloc[-1]['x_r'], y_thrust[-1]))

        # Prepare results
        result = {
            'thrust_line': LineString(thrust_points),
            'normal_forces': N,
            'interslice_forces': Z,
            'residuals': {
                'max_Fx': max(abs(d['Fx_Residual']) for d in debug_data) if debug else None,
                'max_Fy': max(abs(d['Fy_Residual']) for d in debug_data) if debug else None,
                'max_Z_ratio': max(abs(d['Z_ratio']) for d in debug_data) if debug else None
            }
        }

        # Debug output
        if debug:
            debug_df = pd.DataFrame(debug_data)
            result['debug_df'] = debug_df

            print("\n=== EQUILIBRIUM VERIFICATION ===")
            print(f"Max Fx residual: {result['residuals']['max_Fx']:.2e}")
            print(f"Max Fy residual: {result['residuals']['max_Fy']:.2e}")
            print(f"Max Z left/right ratio: {result['residuals']['max_Z_ratio']:.2f}")

            if excel_path:
                with pd.ExcelWriter(excel_path) as writer:
                    # Slice calculations
                    debug_df.to_excel(writer, sheet_name='Slice_Calculations', index=False)

                    # Summary
                    summary = pd.DataFrame({
                        'Parameter': ['FS', 'Theta (deg)', 'Max Fx Error', 'Max Fy Error', 'Mean N', 'Mean |Z|'],
                        'Value': [
                            FS,
                            degrees(theta),
                            result['residuals']['max_Fx'],
                            result['residuals']['max_Fy'],
                            np.mean(N),
                            np.mean(np.abs(Z[1:-1]))  # Exclude boundaries
                        ]
                    })
                    summary.to_excel(writer, sheet_name='Summary', index=False)

                    # Thrust line coordinates
                    thrust_df = pd.DataFrame({
                        'x': [p[0] for p in thrust_points],
                        'y': [p[1] for p in thrust_points],
                        'Location': ['Left Boundary'] +
                                    [f'Between {i + 1}-{i + 2}' for i in range(n_slices)] +
                                    ['Right Boundary']
                    })
                    thrust_df.to_excel(writer, sheet_name='Thrust_Line', index=False)

        return True, result

    except Exception as e:
        return False, f"Error in calculation: {str(e)}"