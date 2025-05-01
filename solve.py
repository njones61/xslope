
import numpy as np
from math import sin, cos, tan, radians, atan, atan2, degrees
from scipy.optimize import minimize_scalar, root_scalar


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

    # Enforce consistent sign (parallel and consistently oriented interslice forces)
    Q_nonzero = Q[np.abs(Q) > 1e-6]
    if len(Q_nonzero) > 0:
        sign = np.sign(Q_nonzero[0])
        Q = sign * np.abs(Q)

    # for debugging: print values per slice
    print('POST Spencer Q values:')
    for i in range(len(Q)):
        print(f"Slice {i}: Q = {Q[i]:.3f}, alpha = {degrees(alpha[i]):.2f}, phi = {degrees(phi[i]):.2f}, c = {c[i]:.2f}, w = {w[i]:.2f}, u = {u[i]:.2f}")

    return Q

def compute_thrust_line_left_sweep(df, Q, theta_deg):
    """
    Computes the y-location of the line of thrust using vector-based moment equilibrium (left to right),
    honoring the sign of Q_i.

    Parameters:
        df (pd.DataFrame): Must include 'x_l', 'x_r', 'x_c', 'y_cb', 'y_cg', 'w'
        Q (np.ndarray): Interslice force magnitudes (length n), signed
        theta_deg (float): Spencer's interslice force inclination (degrees)

    Returns:
        list: y-locations of interslice forces at each internal boundary (length = n+1)
    """
    import numpy as np

    theta = np.radians(theta_deg)
    x_l = df['x_l'].values
    x_r = df['x_r'].values
    x_c = df['x_c'].values
    y_cb = df['y_cb'].values
    y_cg = df['y_cg'].values

    n = len(Q)
    y_thrust = [None] * (n + 1)
    y_thrust[0] = y_cb[0]

    for i in range(n):
        x_cg = x_c[i]
        y_cg_i = y_cg[i]
        y_i = y_thrust[i]
        x_i = x_l[i]
        Q_i = Q[i]
        Q_ip1 = Q[i + 1] if i < n - 1 else 0
        x_ip1 = x_r[i]

        R_vec_i = Q_i * np.array([np.cos(theta), np.sin(theta)])
        r_vec_i = np.array([x_i - x_cg, y_i - y_cg_i])
        M_i = np.cross(r_vec_i, R_vec_i)

        dx = x_ip1 - x_cg
        R_ip1_x = Q_ip1 * np.cos(theta)
        R_ip1_y = Q_ip1 * np.sin(theta)

        if abs(R_ip1_x) > 1e-3:
            y_ip1 = y_cg_i + (M_i - dx * R_ip1_y) / R_ip1_x
        else:
            y_ip1 = y_cg_i

        y_thrust[i + 1] = y_ip1

    return y_thrust

def compute_thrust_line_right_sweep(df, Q, theta_deg):
    """
    Computes the y-location of the line of thrust using vector-based moment equilibrium (right to left),
    honoring the sign of Q_i.

    Parameters:
        df (pd.DataFrame): Must include 'x_l', 'x_r', 'x_c', 'y_cb', 'y_cg', 'w'
        Q (np.ndarray): Interslice force magnitudes (length n), signed
        theta_deg (float): Spencer's interslice force inclination (degrees)

    Returns:
        list: y-locations of interslice forces at each internal boundary (length = n+1)
    """
    import numpy as np

    theta = np.radians(theta_deg)
    x_l = df['x_l'].values
    x_r = df['x_r'].values
    x_c = df['x_c'].values
    y_cb = df['y_cb'].values
    y_cg = df['y_cg'].values

    n = len(Q)
    y_thrust = [None] * (n + 1)
    y_thrust[-1] = y_cb[-1]  # start at right end

    for i in reversed(range(n)):
        x_cg = x_c[i]
        y_cg_i = y_cg[i]
        y_ip1 = y_thrust[i + 1]
        x_ip1 = x_r[i]
        x_i = x_l[i]
        Q_i = Q[i]
        Q_ip1 = Q[i + 1] if i < n - 1 else 0

        # Moment from known force on right
        R_vec_ip1 = Q_ip1 * np.array([np.cos(theta), np.sin(theta)])
        r_vec_ip1 = np.array([x_ip1 - x_cg, y_ip1 - y_cg_i])
        M_ip1 = np.cross(r_vec_ip1, R_vec_ip1)

        dx = x_i - x_cg
        R_i_x = Q_i * np.cos(theta)
        R_i_y = Q_i * np.sin(theta)

        if abs(R_i_x) > 1e-3:
            y_i = y_cg_i + (M_ip1 - dx * R_i_y) / R_i_x
        else:
            y_i = y_cg_i

        y_thrust[i] = y_i

    return y_thrust

def compute_thrust_line_spencer(df, Q, theta_deg):
    """
    Computes the line of thrust using both left and right recursive moment sweeps,
    returning left, right, and averaged vertical thrust line positions.

    Parameters:
        df (pd.DataFrame): Must contain 'x_c', 'y_cb', 'y_cg', 'w', and geometry
        Q (np.ndarray): Interslice force magnitudes from Spencer solution (length n)
        theta_deg (float): Spencer's constant interslice force angle (degrees)

    Returns:
        y_left: list of y-values from left-to-right sweep
        y_right: list of y-values from right-to-left sweep
        y_avg: list of averaged y-values at each interslice interface
    """
    y_left = compute_thrust_line_left_sweep(df, Q, theta_deg)
    y_right = compute_thrust_line_right_sweep(df, Q, theta_deg)

    #y_debug = debug_thrust_line_left_sweep(df, Q, theta_deg)

    y_avg = [
        (yl + yr) / 2 if yl is not None and yr is not None else None
        for yl, yr in zip(y_left, y_right)
    ]

    return y_left, y_right, y_avg

def compute_thrust_line_by_resultant(df, FS):
    """
    Computes the line of thrust by calculating the point of application of the resultant force
    on each slice using local moment balance with y_cg (corrected for dload).

    Parameters:
        df (pd.DataFrame): Slice DataFrame. Must include:
            - 'w', 'alpha', 'phi', 'c', 'dx', 'u', 'y_cb', 'y_cg'
        FS (float): Factor of safety

    Returns:
        list: y-locations (one per slice) where the resultant force acts on the base
    """
    import numpy as np

    w = df['w'].values
    alpha = np.radians(df['alpha'].values)
    phi = np.radians(df['phi'].values)
    c = df['c'].values
    dx = df['dx'].values
    u = df['u'].values
    y_cb = df['y_cb'].values
    y_cg = df['y_cg'].values  # now includes distributed load

    y_thrust = []

    for i in range(len(df)):
        Wi = w[i]
        α = alpha[i]
        ϕ = phi[i]
        ci = c[i]
        dxi = dx[i]
        ui = u[i]
        y_base = y_cb[i]
        y_cg_i = y_cg[i]

        # Weight components
        Wx = Wi * np.sin(α)
        Wy = Wi * np.cos(α)

        # Base resistance (reduced by FS)
        sec_α = 1 / np.cos(α)
        N_base = (Wy - ui * dxi * sec_α) / FS
        S_base = (ci * dxi + N_base * np.tan(ϕ)) / FS

        # Net force vector
        Rx = Wx - S_base
        Ry = Wy - N_base
        R = np.hypot(Rx, Ry)

        # Find point on base where R must act to match moment of W at y_cg
        if R > 1e-3:
            M_w = Wi * (y_cg_i - y_base)  # moment about base center
            y_R = y_base + M_w / R
        else:
            y_R = y_base  # fallback

        y_thrust.append(y_R)

    return y_thrust

def compute_thrust_line_projected_to_base(df, FS):
    """
    Computes the location of the thrust line by balancing moments about the base,
    locating the thrust point along the actual slice base (not just vertical offset).

    Parameters:
        df (pd.DataFrame): Must include:
            'w', 'alpha', 'phi', 'c', 'dx', 'u', 'x_l', 'x_r', 'y_lb', 'y_rb', 'y_cg'
        FS (float): Factor of safety

    Returns:
        list of (x, y): Coordinates on the base where the resultant acts per slice
    """
    import numpy as np

    w = df['w'].values
    alpha = np.radians(df['alpha'].values)
    phi = np.radians(df['phi'].values)
    c = df['c'].values
    dx = df['dx'].values
    u = df['u'].values
    y_cg = df['y_cg'].values
    x_l = df['x_l'].values
    x_r = df['x_r'].values
    y_lb = df['y_lb'].values
    y_rb = df['y_rb'].values

    thrust_points = []

    for i in range(len(df)):
        Wi = w[i]
        α = alpha[i]
        ϕ = phi[i]
        ci = c[i]
        dxi = dx[i]
        ui = u[i]
        ycg = y_cg[i]

        # Geometry of base
        x1, y1 = x_l[i], y_lb[i]
        x2, y2 = x_r[i], y_rb[i]
        dx_base = x2 - x1
        dy_base = y2 - y1
        L = np.hypot(dx_base, dy_base)

        if L == 0:
            thrust_points.append((x1, y1))
            continue

        # Unit vectors
        tx = dx_base / L
        ty = dy_base / L
        nx = -ty
        ny = tx

        # Weight components
        Wx = Wi * np.sin(α)
        Wy = Wi * np.cos(α)

        # Base normal and shear resistance (reduced by FS)
        sec_α = 1 / np.cos(α)
        N_base = (Wy - ui * dxi * sec_α) / FS
        S_base = (ci * dxi + N_base * np.tan(ϕ)) / FS

        # Base reaction vector
        R_base_x = N_base * nx + S_base * tx
        R_base_y = N_base * ny + S_base * ty

        R_mag = np.hypot(R_base_x, R_base_y)
        if R_mag < 1e-6:
            thrust_points.append((0.5 * (x1 + x2), 0.5 * (y1 + y2)))
            continue

        # Moment from weight about (x1, y1)
        rw_x = 0.5 * (x1 + x2) - x1
        rw_y = ycg - y1
        Mw = Wx * rw_y - Wy * rw_x

        # Lever arm distance along R direction
        R_unit_x = R_base_x / R_mag
        R_unit_y = R_base_y / R_mag
        lever = Mw / R_mag

        # Thrust location
        px = x1 + lever * R_unit_y
        py = y1 - lever * R_unit_x

        thrust_points.append((px, py))

    return thrust_points

def compute_thrust_line_by_intersection(df, FS):
    """
    Computes the thrust line by constructing the resultant force vector for each slice
    and finding its intersection with the slice base line.

    Parameters:
        df (pd.DataFrame): Must include:
            'w', 'alpha', 'phi', 'c', 'dx', 'u', 'x_l', 'x_r', 'y_lb', 'y_rb', 'y_cg', 'x_c'
        FS (float): Factor of safety

    Returns:
        list of (x, y): Coordinates of the intersection point on the base per slice
    """
    import numpy as np
    from shapely.geometry import LineString, Point

    w = df['w'].values
    alpha = np.radians(df['alpha'].values)
    phi = np.radians(df['phi'].values)
    c = df['c'].values
    dx = df['dx'].values
    u = df['u'].values
    y_cg = df['y_cg'].values
    x_c = df['x_c'].values
    x_l = df['x_l'].values
    x_r = df['x_r'].values
    y_lb = df['y_lb'].values
    y_rb = df['y_rb'].values

    thrust_points = []

    for i in range(len(df)):
        Wi = w[i]
        α = alpha[i]
        ϕ = phi[i]
        ci = c[i]
        dxi = dx[i]
        ui = u[i]
        xcg = x_c[i]
        ycg = y_cg[i]

        # Base line
        x1, y1 = x_l[i], y_lb[i]
        x2, y2 = x_r[i], y_rb[i]
        base_line = LineString([(x1, y1), (x2, y2)])

        # Weight components
        Wx = Wi * np.sin(α)
        Wy = Wi * np.cos(α)

        # Base resistance (normal + shear) reduced by FS
        sec_α = 1 / np.cos(α)
        N_base = (Wy - ui * dxi * sec_α) / FS
        S_base = (ci * dxi + N_base * np.tan(ϕ)) / FS

        # Net resultant force vector
        Rx = Wx - S_base
        Ry = Wy - N_base
        R_mag = np.hypot(Rx, Ry)

        if R_mag < 1e-6:
            thrust_points.append((xcg, ycg))
            continue

        # Project force vector from CG to intersect base
        scale = 1000
        x_end = xcg + scale * Rx / R_mag
        y_end = ycg + scale * Ry / R_mag
        thrust_vector_line = LineString([(xcg, ycg), (x_end, y_end)])

        intersection = thrust_vector_line.intersection(base_line)
        if isinstance(intersection, Point):
            thrust_points.append((intersection.x, intersection.y))
        else:
            thrust_points.append(((x1 + x2) / 2, (y1 + y2) / 2))  # fallback

    return thrust_points


def compute_thrust_line_stabl5(df, Q, theta_deg):
    """
    Computes the line of thrust using the STABL5 approach:
    Starting from the first slice, recursively determines the y-location where each interslice force
    must act to satisfy moment equilibrium, assuming known Q and constant theta.

    Parameters:
        df (pd.DataFrame): Must include 'x_l', 'x_r', 'x_c', 'y_cb', 'y_cg', 'w'
        Q (np.ndarray): Interslice force magnitudes (length n)
        theta_deg (float): Spencer's constant interslice force inclination (degrees)

    Returns:
        list: y-locations (length = n+1) of interslice forces (line of thrust)
    """
    import numpy as np

    theta = np.radians(theta_deg)
    x_l = df['x_l'].values
    x_r = df['x_r'].values
    x_c = df['x_c'].values
    y_cb = df['y_cb'].values
    y_cg = df['y_cg'].values
    w = df['w'].values

    n = len(Q)
    y_thrust = [None] * (n + 1)
    y_thrust[0] = y_cb[0]  # start thrust at base center of slice 0

    for i in range(n):
        Wi = w[i]
        xi = x_l[i]
        xip1 = x_r[i]
        xc = x_c[i]
        ycg = y_cg[i]
        yb = y_cb[i]

        Q_i = Q[i]
        Q_ip1 = Q[i + 1] if i < n - 1 else 0

        R_i_x = Q_i * np.cos(theta)
        R_ip1_x = Q_ip1 * np.cos(theta)

        M_left = R_i_x * (y_thrust[i] - yb)
        M_w = Wi * (ycg - yb)

        denom = R_ip1_x
        if abs(denom) > 1e-3:
            y_ip1 = (M_left + M_w) / denom + yb
        else:
            print(f"Warning: slice {i}: R_ip1_x is nearly zero, falling back to base y.")
            y_ip1 = yb

        y_thrust[i + 1] = y_ip1

    return y_thrust