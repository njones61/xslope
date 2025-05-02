## MISC THINGS WE TRIED BUT DIDN'T WORK


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


## OLD STUFF

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


### VERSIONS AAS OF 8:30 pm on Thurday

def sweep_thrust_left(df, Q, theta_deg):
    """
    Performs a left-to-right sweep to compute the line of thrust based on interslice force transmission
    using Spencer's assumptions. Assumes all interslice forces are parallel at angle theta.

    Parameters:
        df (pd.DataFrame): Must include 'x_c', 'x_l', 'x_r', 'y_cb', 'y_lb', 'w'
        Q (np.ndarray): Net interslice force resultants for each slice (length n)
        theta_deg (float): Spencer interslice force inclination in degrees

    Returns:
        y_thrust (list): y-locations of interslice forces (length = n+1)
    """
    import numpy as np

    theta = np.radians(theta_deg)
    x_c = df['x_c'].values
    x_l = df['x_l'].values
    x_r = df['x_r'].values
    y_cb = df['y_cb'].values
    y_lb = df['y_lb'].values
    w = df['w'].values

    n = len(Q)
    y_thrust = [None] * (n + 1)

    # Initialize: no force entering from the left
    Z_left_x = 0.0
    Z_left_y = 0.0
    y_thrust[0] = y_lb[0]  # start thrust at left base elevation

    for i in range(n):

        if i == n - 1:
            y_thrust[i + 1] = y_cb[i]  # anchor final point to base center
            break

        Wi = w[i]
        yb = y_cb[i]
        xc = x_c[i]
        xl = x_l[i]
        xr = x_r[i]

        # Spencer Q force components
        Qx = Q[i] * np.cos(theta)
        Qy = Q[i] * np.sin(theta)

        # Right-side force = Q - left
        Z_right_x = Qx - Z_left_x
        Z_right_y = Qy - Z_left_y

        # Moment arm for Z_left
        dy_left = y_thrust[i] - yb if y_thrust[i] is not None else 0.0
        dx_left = xc - xl

        # Moment from incoming interslice force
        M_left = Z_left_x * dy_left - Z_left_y * dx_left

        # Moment from weight = 0 (acts through x_c)
        M_w = 0.0

        # Moment arm for Z_right_y
        dx_right = xr - xc

        # Compute vertical offset to balance moment
        if abs(Z_right_x) > 1e-6:
            dy_right = (M_left + M_w + Z_right_y * dx_right) / Z_right_x
        else:
            dy_right = 0.0

        y_thrust[i + 1] = yb + dy_right

        # Flip force direction to pass to next slice
        Z_left_x = -Z_right_x
        Z_left_y = -Z_right_y

    return y_thrust

def sweep_thrust_left_debug(df, Q, theta_deg):
    """
    Debug version of sweep_thrust_left. Prints intermediate moment terms per slice.
    """
    import numpy as np

    theta = np.radians(theta_deg)
    x_c = df['x_c'].values
    x_l = df['x_l'].values
    x_r = df['x_r'].values
    y_cb = df['y_cb'].values
    y_lb = df['y_lb'].values
    w = df['w'].values

    n = len(Q)
    y_thrust = [None] * (n + 1)

    Z_left_x = 0.0
    Z_left_y = 0.0
    y_thrust[0] = y_lb[0]

    print(f"{'Slice':>5} | {'ZLx':>8} | {'ZLy':>8} | {'ZRy':>8} | {'dx_R':>8} | {'dy_L':>8} | {'M_L':>10} | {'Δy_R':>10} | {'y_R':>10}")
    print("-" * 90)

    for i in range(n):
        if i == n - 1:
            y_thrust[i + 1] = y_cb[i]
            break

        Wi = w[i]
        yb = y_cb[i]
        xc = x_c[i]
        xl = x_l[i]
        xr = x_r[i]

        Qx = Q[i] * np.cos(theta)
        Qy = Q[i] * np.sin(theta)

        Z_right_x = Qx - Z_left_x
        Z_right_y = Qy - Z_left_y

        dy_left = y_thrust[i] - yb if y_thrust[i] is not None else 0.0
        dx_left = xc - xl
        M_left = Z_left_x * dy_left - Z_left_y * dx_left

        M_w = 0.0
        dx_right = xr - xc

        if abs(Z_right_x) > 1e-6:
            dy_right = (M_left + M_w + Z_right_y * dx_right) / Z_right_x
        else:
            dy_right = 0.0

        y_next = yb + dy_right
        y_thrust[i + 1] = y_next

        print(f"{i:5d} | {Z_left_x:8.2f} | {Z_left_y:8.2f} | {Z_right_y:8.2f} | {dx_right:8.2f} | {dy_left:8.2f} | {M_left:10.2f} | {dy_right:10.2f} | {y_next:10.2f}")

        Z_left_x = -Z_right_x
        Z_left_y = -Z_right_y

    return y_thrust

def sweep_thrust_right(df, Q, theta_deg):
    """
    Performs a right-to-left sweep to compute the line of thrust based on interslice force transmission
    using Spencer's assumptions. Assumes all interslice forces are parallel at angle theta.

    Parameters:
        df (pd.DataFrame): Must include 'x_c', 'y_cb', 'y_cg', 'w'
        Q (np.ndarray): Net interslice force resultants for each slice (length n)
        theta_deg (float): Spencer interslice force inclination in degrees

    Returns:
        y_thrust (list): y-locations of interslice forces (length = n+1)
    """
    import numpy as np

    theta = np.radians(theta_deg)
    x_c = df['x_c'].values
    y_cb = df['y_cb'].values
    y_cg = df['y_cg'].values
    w = df['w'].values

    n = len(Q)
    y_thrust = [None] * (n + 1)

    # Initialize: no force entering right side of last slice
    Z_right_x = 0.0
    Z_right_y = 0.0
    y_thrust[-1] = y_cb[-1]  # base center of last slice

    for i in reversed(range(n)):

        if i == 0:  # First slice
            y_thrust[i] = y_cb[i] # Set final thrust point at base of first slice
            break

        Wi = w[i]
        yb = y_cb[i]
        ycg = y_cg[i]

        Qx = Q[i] * np.cos(theta)
        Qy = Q[i] * np.sin(theta)

        # Left-side force = Q - right-side force
        Z_left_x = Qx - Z_right_x
        Z_left_y = Qy - Z_right_y

        # Moment from right interslice force (about base center)
        M_right = Z_right_x * (y_thrust[i + 1] - yb) if y_thrust[i + 1] is not None else 0.0

        # Moment from weight
        x_cg = x_c[i]  # placeholder: weight acts at center of base, need to adjust later
        x_base = x_c[i]  # center of base
        M_w = Wi * (x_cg - x_base)  # should be zero for now

        # Solve for thrust location on left
        denom = Z_left_x
        if abs(denom) > 1e-6:
            y_left = yb + (M_right + M_w) / Z_left_x
        else:
            y_left = yb

        y_thrust[i] = y_left

        # Propagate force to next slice
        Z_right_x = -Z_left_x
        Z_right_y = -Z_left_y

    return y_thrust