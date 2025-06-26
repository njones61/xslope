import numpy as np
import pandas as pd
from shapely.geometry import LineString, Point
from math import sin, cos, tan, radians, atan, atan2, degrees
from scipy.optimize import minimize_scalar, root_scalar, newton
from tabulate import tabulate


def oms(df, debug=False):
    """
    Computes FS by direct application of Equation 9 (Ordinary Method of Slices).

    Inputs
    ------
    df : pandas.DataFrame
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

    NOTES
    -----
    Implements exactly:
      FS
      = [ Σ { cᵢ·Δℓᵢ
            + [ Wᵢ·cosαᵢ + Dᵢ·cos(αᵢ−βᵢ) − kWᵢ·sinαᵢ − Tᵢ·sinαᵢ − uᵢ·Δℓᵢ ]·tanφᵢ
            + pᵢ } ]
        / [  Σ(Wᵢ·sinαᵢ)
           + (1/R)·Σ[ Dᵢ·cosβᵢ·(Xo - d_{x,i})  −  Dᵢ·sinβᵢ·(Yo - d_{y,i}) ]
           + (1/R)·Σ[ kWᵢ·(Yo - y_{cg,i}) ]
           + (1/R)·Σ[ Tᵢ·(Yo - y_{t,i}) ]  ].

    """
    if 'r' not in df.columns:
        return False, "Circle is required for OMS method."

    # 1) Unpack circle‐center and radius as single values
    Xo = df['xo'].iloc[0]    # Xoᵢ (x-coordinate of circle center)
    Yo = df['yo'].iloc[0]    # Yoᵢ (y-coordinate of circle center)
    R  = df['r'].iloc[0]     # Rᵢ (radius of circular failure surface)

    # 2) Pull arrays directly from df
    alpha_deg = df['alpha'].values    # αᵢ in degrees
    phi_deg   = df['phi'].values      # φᵢ in degrees
    c     = df['c'].values        # cᵢ
    W     = df['w'].values        # Wᵢ
    u     = df['u'].values        # uᵢ (pore‐force per unit length)
    dl     = df['dl'].values       # Δℓᵢ
    D     = df['dload'].values        # Dᵢ
    d_x    = df['d_x'].values      # d_{x,i}
    d_y    = df['d_y'].values      # d_{y,i}
    beta_deg  = df['beta'].values     # βᵢ in degrees
    kw    = df['kw'].values       # kWᵢ
    T     = df['t'].values        # Tᵢ (zero except one slice)
    y_t    = df['y_t'].values      # y_{t,i} (zero except one slice)
    P     = df['p'].values        # pᵢ
    x_c    = df['x_c'].values      # x_{c,i}
    y_cg = df['y_cg'].values       # y_{cg,i} coordinate of slice centroid

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
    #                               + pᵢ  ]
    #

    N_eff = (
        W * cos_alpha
      + D * cos_ab
      - kw * sin_alpha
      - T * sin_alpha
      - (u * dl)
    )  # N′ᵢ = Wᵢ·cosαᵢ + Dᵢ·cos(αᵢ−βᵢ) − kWᵢ·sinαᵢ − Tᵢ·sinαᵢ − uᵢ·Δℓᵢ

    numerator = np.sum(c * dl + N_eff * tan_phi + P)

    # ————————————————————————————————————————————————————————
    # 6) Build each piece of the DENOMINATOR exactly as Eqn 9:

    #  (A) = Σ [ Wᵢ · sinαᵢ ]
    sum_W = np.sum(W * sin_alpha)

    #  (B) = Σ [ Dᵢ·cosβᵢ·(Xo - d_{x,i})  −  Dᵢ·sinβᵢ·(Yo - d_{y,i}) ]
    a_dx = d_x - Xo
    a_dy = Yo - d_y
    sum_Dx = np.sum(D * np.cos(beta) * a_dx)
    sum_Dy = np.sum(D * np.sin(beta) * a_dy)

    #  (C) = Σ [ kWᵢ * (Yo - y_{cg,i}) ]
    a_s = Yo - y_cg
    sum_kw = np.sum(kw * a_s)

    #  (D) = Σ [ Tᵢ * (Yo - y_{t,i}) ]
    a_t = Yo - y_t
    sum_T = np.sum(T * a_t)

    # Put them together with their 1/R factors:
    denominator = sum_W + (1.0 / R) * (sum_Dx - sum_Dy + sum_kw + sum_T)

    # 7) Finally compute FS = (numerator)/(denominator)
    FS = numerator / denominator

    # 8) Store effective normal forces in the DataFrame
    df['n_eff'] = N_eff

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

def bishop(df, debug=False, tol=1e-6, max_iter=100):
    """
    Computes FS using the complete Bishop's Simplified Method (Equation 10) and computes N_eff (Equation 8).
    Requires circular slip surface and full input data structure consistent with OMS.

    Parameters:
        df : pandas.DataFrame with required columns (see OMS spec)
        debug : bool, if True prints diagnostic info
        tol : float, convergence tolerance
        max_iter : int, maximum iteration steps

    Returns:
        (bool, dict | str): (True, {'method': 'bishop', 'FS': value}) or (False, error message)
    """

    if 'r' not in df.columns:
        return False, "Circle is required for Bishop method."

    # 1) Unpack circle‐center and radius as single values
    Xo = df['xo'].iloc[0]    # Xoᵢ (x-coordinate of circle center)
    Yo = df['yo'].iloc[0]    # Yoᵢ (y-coordinate of circle center)
    R  = df['r'].iloc[0]     # Rᵢ (radius of circular failure surface)

    # Load input arrays
    alpha = np.radians(df['alpha'].values)
    phi   = np.radians(df['phi'].values)
    c     = df['c'].values
    W     = df['w'].values
    u     = df['u'].values
    dl    = df['dl'].values
    D     = df['dload'].values
    d_x   = df['d_x'].values
    d_y   = df['d_y'].values
    beta  = np.radians(df['beta'].values)
    kw    = df['kw'].values
    T     = df['t'].values
    y_t   = df['y_t'].values
    P     = df['p'].values
    x_c   = df['x_c'].values
    y_cg  = df['y_cg'].values

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
    denominator = sum_W + (1.0 / R) * (sum_Dx - sum_Dy + sum_kw + sum_T)

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
        numer_slice = shear / denom_N
        F_new = np.sum(numer_slice) / denominator

        if abs(F_new - F) < tol:
            df['n_eff'] = N_eff
            if debug:
                print(f"FS = {F_new:.6f}")
                print(f"Numerator = {np.sum(numer_slice):.6f}")
                print(f"Denominator = {denominator:.6f}")
                print("N_eff =", np.array2string(N_eff, precision=4, separator=', '))
            return True, {'method': 'bishop', 'FS': F_new}

        F = F_new

    return False, "Bishop method did not converge within the maximum number of iterations."

def janbu(df, debug=False):
    """
    Computes FS using Janbu's Simplified Method with correction factor (Equation 7).

    Implements the complete formulation including distributed loads, seismic forces,
    reinforcement, and tension crack water forces. Applies Janbu correction factor
    based on d/L ratio and soil type.

    Parameters:
        df : pandas.DataFrame with required columns (see OMS spec)
        debug : bool, if True prints diagnostic info

    Returns:
        (bool, dict | str): (True, {'method': 'janbu_simplified', 'FS': value, 'fo': correction_factor})
                           or (False, error message)
    """

    # Load input arrays
    alpha = np.radians(df['alpha'].values)
    phi = np.radians(df['phi'].values)
    c = df['c'].values
    W = df['w'].values
    u = df['u'].values
    dl = df['dl'].values
    D = df['dload'].values
    beta = np.radians(df['beta'].values)
    kw = df['kw'].values
    T = df['t'].values
    P = df['p'].values

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
    x_l = df['x_l'].iloc[0]  # leftmost x
    y_lt = df['y_lt'].iloc[0]  # leftmost top y
    x_r = df['x_r'].iloc[-1]  # rightmost x
    y_rt = df['y_rt'].iloc[-1]  # rightmost top y

    # Length of failure surface (straight line approximation)
    L = np.hypot(x_r - x_l, y_rt - y_lt)

    # Calculate perpendicular distance from each slice center to failure surface line
    x0 = df['x_c'].values
    y0 = df['y_cb'].values

    # Distance from point to line formula: |ax + by + c| / sqrt(a² + b²)
    # Line equation: (y_rt - y_lt)x - (x_r - x_l)y + (x_r * y_lt - y_rt * x_l) = 0
    numerator_dist = np.abs((y_rt - y_lt) * x0 - (x_r - x_l) * y0 + x_r * y_lt - y_rt * x_l)
    dists = numerator_dist / L
    d = np.max(dists)  # maximum perpendicular distance

    dL_ratio = d / L

    # Determine b1 factor based on soil type
    phi_sum = df['phi'].sum()
    c_sum = df['c'].sum()

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
    df['n_eff'] = N_eff

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
    D       = df['dload'].values
    beta    = np.radians(df['beta'].values)
    kw      = df['kw'].values
    T       = df['t'].values
    P       = df['p'].values
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

    df['n_eff'] = N  # store effective normal forces in df
    df['z'] = Z[:-1]  # store interslice forces in df, adjust length to n slices

    if debug:
        r_opt = residual(FS_opt)
        print(f" Converged FS = {FS_opt:.6f}, residual = {r_opt:.4g}")

    return True, {'FS': FS_opt}

def corps_engineers(df, debug=False):
    """
    Corps of Engineers style force equilibrium solver.

    1. Computes a single θ from the slope between
       (x_l[0], y_lt[0]) and (x_r[-1], y_rt[-1]).
    2. Builds a constant θ array of length n+1.
    3. Calls force_equilibrium(df, theta_array).

    Parameters:
        df (pd.DataFrame): Must include at least ['x_l','y_lt','x_r','y_rt']
                           plus all the columns required by force_equilibrium:
                           ['alpha','phi','c','dl','w','u','dx'].

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

def lowe_karafiath(df, debug=False):
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

def spencer(df, tol=1e-6, debug=False):
    """
    Spencer's Method using Steve G. Wright's formulation.
    Solves for FS_force and FS_moment independently using the Wright Q equation.

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

    # Check if the 'r' column is present in df. If so, circular = True.
    circular = 'r' in df.columns

    tol = 1e-6
    max_iter = 100

    alpha = np.radians(df['alpha'].values)
    phi   = np.radians(df['phi'].values)
    c     = df['c'].values
    dx    = df['dx'].values
    dl    = df['dl'].values
    w     = df['w'].values
    u     = df['u'].values
    x_c   = df['x_c'].values
    y_cb  = df['y_cb'].values
    D     = df['dload'].values
    beta  = np.radians(df['beta'].values)
    kw    = df['kw'].values
    T     = df['t'].values
    P     = df['p'].values
    cos_a = np.cos(alpha)
    sin_a = np.sin(alpha)
    tan_p = np.tan(phi)

    def compute_Q(F, theta_rad):
        term1 = w * sin_a + kw * cos_a + T * cos_a - P + D * np.sin(alpha - beta)
        term2 = (c / F) * dl
        term3 = (w * cos_a + D * np.cos(alpha - beta) - kw * sin_a - T * sin_a - u * dl) * tan_p / F
        numerator = term1 - term2 - term3
        theta_diff = alpha - theta_rad
        denominator = np.cos(theta_diff) * (1 + (np.tan(theta_diff) * tan_p) / F)
        Q = numerator / denominator
        return Q

    ## STABL5 Formulation
    # def compute_Q(F, theta_rad):
    #     Ca = c * dl
    #     Ua = u * dl
    #     S1 = Ca + tan_p * (w * cos_a  - kw * sin_a - Ua + D * np.cos(alpha - beta))
    #     S2 =  D * np.sin(alpha - beta) - w * sin_a - kw * cos_a 
    #     theta_diff = alpha - theta_rad
    #     S3 = tan_p * np.tan(theta_diff)
    #     Q = (S1/F + S2) / (np.cos(theta_diff) * (1 + S3/F))
    #     return Q

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
        theta_rad = np.radians(theta_deg)
        Ff = fs_force(theta_rad)
        Fm = fs_moment(theta_rad)
        return Ff - Fm

    # First try Newton's method to find theta_opt
    try:
        theta_guess = 10.0  # This seems to work well for most cases
        if debug:
            print(f"Trying Newton's method for theta root-finding with initial guess {theta_guess:.1f} deg")
        theta_opt = newton(fs_difference, x0=theta_guess, tol=tol, maxiter=max_iter)
        # Check if the solution is valid
        if (
            abs(theta_opt) > 59 or
            abs(fs_difference(theta_opt)) > 0.01 or
            fs_force(np.radians(theta_opt)) >= fs_max - 1e-3
        ):
            raise ValueError("Newton's method converged to an invalid solution, trying sweep...")
    except Exception:
        # If Newton's method fails, try sweeping theta to find a bracket where fs_difference changes sign,
        #  then use root_scalar to find theta_opt
        if debug: 
            print("Newton's method failed or gave invalid result, sweeping theta to find a bracket...")
        theta_range = np.linspace(-60, 60, 49)
        fs_diff_values = [fs_difference(theta) for theta in theta_range]
        bracket = None
        for i in range(len(fs_diff_values) - 1):
            if fs_diff_values[i] * fs_diff_values[i+1] < 0:
                bracket = [theta_range[i], theta_range[i+1]]
                break
        if bracket is None:
            return False, "Spencer's method: No sign change found in fs_difference over theta range."
        try:
            if debug:
                print(f"Trying root_scalar with bracket {bracket}")
            sol = root_scalar(fs_difference, bracket=bracket, method='brentq', xtol=tol)
            theta_opt = sol.root
        except Exception as e:
            return False, f"Spencer's method failed to converge: {e}"

    theta_rad = np.radians(theta_opt)
    FS_force = fs_force(theta_rad)
    FS_moment = fs_moment(theta_rad)

    df['theta'] = theta_opt  # store theta in df.

    # Simplified method for computing N_eff
    Q = compute_Q(FS_force, theta_rad)
    N_eff = w * cos_a + D * np.cos(alpha - beta) + Q * np.sin(alpha - theta_rad) - kw * sin_a - T * sin_a - u * dl 

    # ---  compute interslice forces Z  ---
    n = len(Q)
    Z = np.zeros(n+1)
    for i in range(n):
        Z[i+1] = Z[i] - Q[i] 

    # store back into df
    df['z']     = Z[:-1]        # Z_i acting on slice i's left face
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
        if debug:
            for i in range(len(Q)):
                print(f"Slice {i}: Q = {Q[i]:.3f}, alpha = {np.degrees(alpha[i]):.2f}, phi = {np.degrees(phi[i]):.2f}, c = {c[i]:.2f}, w = {w[i]:.2f}, u = {u[i]:.2f}")

        return True, results


def compute_line_of_thrust_sweep(df, FS, debug=True):
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
            'dl', 'dx', 'x_l', 'x_r', 'y_lb', 'y_rb',
            'd', 'beta', 'kw', 't', 'p',
            'theta'
    FS : float
        Factor of safety from Spencer's solution.
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
    phi    = np.radians(df['phi'].values)
    theta  = np.zeros(n+1)  # interslice force inclination
    theta[:-1]   = np.radians(df['theta'].values)  # interslice force inclination
    c      = df['c'].values  # cohesion
    w      = df['w'].values  # slice weight
    u      = df['u'].values  # pore pressure
    dl     = df['dl'].values  # slice base length
    dx     = df['dx'].values  # slice width
    x_r    = df['x_r'].values  # right side x-coordinate
    x_l    = df['x_l'].values  # left side x-coordinate
    x_c    = df['x_c'].values  # center of base x-coordinate
    y_lb   = df['y_lb'].values  # left side base y-coordinate
    y_rb   = df['y_rb'].values  # right side base y-coordinate
    y_cb   = df['y_cb'].values  # center base y-coordinate
    D      = df['dload'].values     # distributed load
    d_x    = df['d_x'].values   # distributed load x-coordinate
    d_y    = df['d_y'].values   # distributed load y-coordinate
    beta   = np.radians(df['beta'].values)  # distributed load inclination
    kw     = df['kw'].values  # seismic force
    y_cg   = df['y_cg'].values  # seismic force y-coordinate
    T      = df['t'].values  # tension crack water force
    y_t    = df['y_t'].values  # tension crack water force y-coordinate

    c_m       = c / FS    # mobilized cohesion
    tan_phi_m = np.tan(phi) / FS  # mobilized friction angle

    N_eff = df['n_eff'].values  # effective normal force
    Z     = np.zeros(n+1)  # interslice force
    Z[:-1]  = df['z'].values  # interslice force on left face

    # Determine if the slope is left-facing or right-facing
    right_facing = (y_lb[0] > y_rb[-1])
    if right_facing:  # if right-facing, flip the interslice force and theta
        Z = -Z
        theta = -theta
        beta = -beta

    tol = 1e-8

    # 1) decompose side-forces
    X = Z * np.sin(theta)
    E = Z * np.cos(theta)

    # 2a) left-to-right moment sweep about each slice's lower-right
    delta_y_L    = np.zeros(n+1)  # Moment arms to side forces relative to slice lower-right.
    y_L          = np.zeros(n+1)  # Absolute y-coordinates of the thrust line on all slice boundaries based on left sweep
    delta_y_L[0] = 0.0         # First slice pivot (starting from the left)
    y_L[0]       = y_lb[0]      # First slice left side y-coordinate

    for i in range(n-1):  # Loop from 0 to n-2. Last slice = n-1 and the right-side moment arm on that slice is fixed.

        arm_left = y_L[i] - y_rb[i] #  left-side moment arm E[i] (pivot at y_rb[i])
        N = N_eff[i] + u[i] * dl[i]  #  normal force
        num = (
            E[i] * arm_left  
            + X[i] * dx[i]
            - w[i] * dx[i] / 2
            + N * dl[i] / 2
            - D[i] * np.cos(beta[i]) * (x_r[i] - d_x[i])  # D*cos(β)*a_dx relative to right corner
            + D[i] * np.sin(beta[i]) * (d_y[i] - y_rb[i])  # D*sin(β)*a_dy relative to right corner
            - kw[i] * (y_cg[i] - y_rb[i])  # kW*a_k relative to right corner
            - T[i] * (y_t[i] - y_rb[i])    # T*a_t relative to right corner
        )
        if D[i] != 0:
            print(f"LEFT SWEEP:for slice {i}, a_dx = {x_r[i] - d_x[i]:.2f}, a_dy = {d_y[i] - y_rb[i]:.2f}")
            print(f"beta = {beta[i]:.2f}, beta(deg) = {np.degrees(beta[i]):.2f}")
        delta_y_L[i+1]  = num / E[i+1] if abs(E[i+1]) > tol else 0  # right-side moment arm
        y_L[i+1]  = y_rb[i] + delta_y_L[i+1]    # absolute y value on right side

    delta_y_L[n] = 0.0         # Last (right-most) slice pivot
    y_L[n]       = y_rb[n-1]   # Last (right-most) slice right side y-coordinate

    # 2b) right-to-left moment sweep about each slice's lower-left
    delta_y_R    = np.zeros(n+1)   # Moment arms to side forces relative to slice lower-left.
    y_R          = np.zeros(n+1)   # Absolute y-coordinates of the thrust line on all slice boundaries based on right sweep
    delta_y_R[n] = 0.0          # First slice pivot (starting from the right)
    y_R[n]       = y_rb[n-1]     # First slice left side y-coordinate

    for i in range(n-1, 0, -1):  # Loop from n-1 to 1. Last slice = 0 and the left-side moment arm on that slice is fixed.

        arm_right = y_R[i+1] - y_lb[i] #  right-side moment arm for E[i+1] (pivot at y_lb[i])
        N = N_eff[i] + u[i] * dl[i]  #  normal force
        num = (
            E[i+1] * arm_right
            - X[i+1] * dx[i]
            - w[i] * dx[i] / 2
            + N * dl[i] / 2
            - D[i] * np.cos(beta[i]) * (d_x[i] - x_l[i])  # D*cos(β)*a_dx relative to left corner
            - D[i] * np.sin(beta[i]) * (d_y[i] - y_lb[i]) # D*sin(β)*a_dy relative to left corner
            + kw[i] * (y_cg[i] - y_lb[i])  # kW*a_k relative to left corner
            + T[i] * (y_t[i] - y_lb[i])    # T*a_t relative to left corner
        )
        if D[i] != 0:
            print(f"Right SWEEP:for slice {i}, a_dx = {d_x[i] - x_l[i]:.2f}, a_dy = {d_y[i] - y_lb[i]:.2f}")
            print(f"beta = {beta[i]:.2f}, beta(deg) = {np.degrees(beta[i]):.2f}")
        delta_y_R[i] = num / E[i] if abs(E[i]) > tol else 0
        y_R[i] = y_lb[i] + delta_y_R[i]

    delta_y_R[0] = 0.0       # Last (left-most) slice pivot
    y_R[0] = y_lb[0]         # Last (left-most) slice left side y-coordinate

    # 2c) average both sweeps
    y_bound = 0.5 * (y_L + y_R)

    # 3) build LineString
    x_bound = np.empty(n+1)
    x_bound[0] = df['x_l'].iat[0]
    x_bound[1:] = df['x_r'].values
    thrust_line = LineString(np.column_stack([x_bound, y_bound]))

    # 4) debug export  (OUT OF DATE!!!)
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
                    + D[i] * np.cos(beta[i]) * (x_r[i] - d_x[i])  # D*cos(β)*a_dx relative to right corner
                    - D[i] * np.sin(beta[i]) * (d_y[i] - y_rb[i])  # D*sin(β)*a_dy relative to right corner
                    + kw[i] * (y_cg[i] - y_rb[i])  # kW*a_k relative to right corner
                    + T[i] * (y_t[i] - y_rb[i])    # T*a_t relative to right corner
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
                    - D[i] * np.cos(beta[i]) * (d_x[i] - x_l[i])  # D*cos(β)*a_dx relative to left corner
                    - D[i] * np.sin(beta[i]) * (d_y[i] - y_lb[i]) # D*sin(β)*a_dy relative to left corner
                    + kw[i] * (y_cg[i] - y_lb[i])  # kW*a_k relative to left corner
                    + T[i] * (y_t[i] - y_lb[i])    # T*a_t relative to left corner
                )

            rows.append({
                'slice':           i,
                'N_eff':           N_eff[i],
                'dl':              dl[i],
                'dx':              dx[i],
                'x_c':             x_c[i],
                'x_l':             x_l[i],
                'x_r':             x_r[i],
                'y_lb':            y_lb[i],
                'y_rb':            y_rb[i],
                'y_cb':            y_cb[i],
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
                'D':               D[i],
                'beta':            beta[i],
                'd_x':             d_x[i],
                'd_y':             d_y[i],
                'kw':              kw[i],
                'y_cg':            y_cg[i],
                'y_t':             y_t[i],
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

def compute_line_of_thrust_center(df, FS, debug=True):
    """
    Compute the line of thrust using the "center of base of slice" method.
    This implements equation 8 from the line_of_thrust.md documentation:
    
    Δy_{i+1} = (E_i * Δy_i + X_i * Δx/2 + X_{i+1} * Δx/2 - D * cos(β) * d_ax + D * sin(β) * d_ay - kW * a_k - T * a_t) / E_{i+1}
    
    This method sums moments about the center of the base of each slice.
    Only the side forces create moments since W, S, and N all go through the base center.
    Additional forces (D, kW, T) are included with their respective moment arms.
    
    Parameters
    ----------
    df : pd.DataFrame
        Must include columns:
            'alpha', 'phi', 'c', 'w', 'u',
            'dl', 'dx', 'x_l', 'x_r', 'y_lb', 'y_rb',
            'dload', 'beta', 'kw', 't', 'p',
            'theta', 'n_eff', 'z', 'd_x', 'd_y', 'y_cg', 'y_t'
    FS : float
        Factor of safety from Spencer's solution.
    debug : bool
        If True, export a per-slice debug table to Excel.

    Returns
    -------
    dict
        'sigma_eff'   : N_eff / dl per slice,
        'delta_y'     : moment arms at each boundary (n+1),
        'thrust_line' : Shapely LineString of the y coordinates.
    """
    # 1) unpack & force equilibrium
    n      = len(df)
    theta  = np.zeros(n+1)  # interslice force inclination
    theta[:-1]   = np.radians(df['theta'].values)  # interslice force inclination
    dx     = df['dx'].values  # slice width
    x_c    = df['x_c'].values  # center of base x-coordinate
    x_l    = df['x_l'].values  # left side x-coordinate
    x_r    = df['x_r'].values  # right side x-coordinate
    y_cb   = df['y_cb'].values  # center base y-coordinate
    y_lb   = df['y_lb'].values  # left side base y-coordinate
    y_rb   = df['y_rb'].values  # right side base y-coordinate
    D      = df['dload'].values     # distributed load
    d_x    = df['d_x'].values   # distributed load x-coordinate
    d_y    = df['d_y'].values   # distributed load y-coordinate
    beta   = np.radians(df['beta'].values)  # distributed load inclination
    kw     = df['kw'].values  # seismic force
    y_cg   = df['y_cg'].values  # seismic force y-coordinate
    T      = df['t'].values  # tension crack water force
    y_t    = df['y_t'].values  # tension crack water force y-coordinate

    Z     = np.zeros(n+1)  # interslice force
    Z[:-1]  = df['z'].values  # interslice force on left face

    # Determine if the slope is left-facing or right-facing
    right_facing = (y_lb[0] > y_rb[-1])
    if right_facing:  # if right-facing, flip the interslice force and theta
        Z = -Z
        theta = -theta
        beta = -beta

    tol = 1e-8

    # 1) decompose side-forces
    X = Z * np.sin(theta)
    E = Z * np.cos(theta)

    # 2) Compute moment arms using center of base method with additional forces
    # Start from left side and sweep to right
    delta_y = np.zeros(n+1)  # Moment arms to side forces relative to center of base
    y_thrust = np.zeros(n+1)  # Absolute y-coordinates of the thrust line
    
    

    y_thrust[0] = y_lb[0]  # First slice left side y-coordinate
    
    for i in range(n):
        # For each slice, solve for the right side moment arm using equation 8
        # Δy_{i+1} = (E_i * Δy_i + X_i * Δx/2 + X_{i+1} * Δx/2 - D * cos(β) * a_dx + D * sin(β) * a_dy - kW * a_k - T * a_t) / E_{i+1}
        
        if i < n-1:  # Not the last slice

            if i == 0:
                delta_yi = 0.0     # For first slice, left side force is zero
            else:
                delta_yi = y_thrust[i] - y_cb[i]  # need to recompute this for each slice as y_cb changes

            # Calculate moment arms for additional forces relative to center of base
            # a_dx = horizontal distance from point d to pivot point (center of base)
            a_dx = x_c[i] - d_x[i]
            
            # a_dy = vertical distance from point d to pivot point (center of base)
            a_dy = d_y[i] - y_cb[i]
            
            # a_k = vertical distance from c.g. to pivot point (center of base)
            a_k = y_cg[i] - y_cb[i]
            
            # a_t = vertical distance from tension crack to pivot point (center of base)
            a_t =y_t[i]- y_cb[i]
            
            numerator = (E[i] * delta_yi + 
                        X[i] * dx[i] / 2 + 
                        X[i+1] * dx[i] / 2 -
                        D[i] * np.cos(beta[i]) * a_dx +
                        D[i] * np.sin(beta[i]) * a_dy -
                        kw[i] * a_k -
                        T[i] * a_t)
            
            if abs(E[i+1]) > tol:
                delta_yip1 = numerator / E[i+1]
            else:
                delta_yip1 = 0.0
                
            # Convert to absolute y coordinate
            y_thrust[i+1] = y_cb[i] + delta_yip1
        else:
            # For the last slice, pass through the right side y-coordinate
            y_thrust[i+1] = y_rb[i]

    # 3) build LineString
    x_bound = np.empty(n+1)
    x_bound[0] = df['x_l'].iat[0]
    x_bound[1:] = df['x_r'].values
    thrust_line = LineString(np.column_stack([x_bound, y_thrust]))

    # 4) debug export
    if debug:
        rows = []
        for i in range(n):
            # Calculate moment residual for verification
            if i < n-1:
                # Calculate moment arms for additional forces
                a_dx = x_c[i] - d_x[i]
                a_dy = d_y[i] - y_cb[i]
                a_k = y_cg[i] - y_cb[i]
                a_t = y_t[i] - y_cb[i] 

                delta_yi = y_thrust[i] - y_cb[i]
                delta_yip1 = y_thrust[i+1] - y_cb[i]
                
                m_res = (-E[i] * delta_yi  - 
                        X[i] * dx[i] / 2 + 
                        E[i+1] * delta_yip1 - 
                        X[i+1] * dx[i] / 2 +
                        D[i] * np.cos(beta[i]) * a_dx -
                        D[i] * np.sin(beta[i]) * a_dy +
                        kw[i] * a_k +
                        T[i] * a_t)
            else:
                m_res = 0.0

            rows.append({
                'slice':           i,
                'dx':              dx[i],
                'x_c':             x_c[i],
                'x_l':             x_l[i],
                'x_r':             x_r[i],
                'y_lb':            y_lb[i],
                'y_rb':            y_rb[i],
                'y_cb':            y_cb[i],
                'Z_left':          Z[i],
                'Z_right':         Z[i+1],
                'X_left':          X[i],
                'X_right':         X[i+1],
                'E_left':          E[i],
                'E_right':         E[i+1],
                'D':               D[i],
                'beta':            beta[i],
                'd_x':             d_x[i],
                'd_y':             d_y[i],
                'kw':              kw[i],
                'y_cg':            y_cg[i],
                'y_t':             y_t[i],
                'a_dx':            a_dx if i < n-1 else 0.0,
                'a_dy':            a_dy if i < n-1 else 0.0,
                'a_k':             a_k if i < n-1 else 0.0,
                'a_t':             a_t if i < n-1 else 0.0,
                'delta_y_left':    delta_yi,
                'delta_y_right':   delta_yip1,
                'y_thrust_left':   y_thrust[i],
                'y_thrust_right':  y_thrust[i+1],
                'moment_residual': m_res
            })

        df_debug = pd.DataFrame(rows)
        df_debug.to_excel("thrust_calc_results_center_base.xlsx", index=False)
        print("Debug table written to thrust_calc_results_center_base.xlsx")

    # 5) return results
    return thrust_line

def spencer2(df, tol=1e-4, max_iter = 100, debug=True):
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
    P     = df['dload'].values  # distributed load resultant 
    beta  = np.radians(df['beta'].values)  # distributed load inclination, degrees
    kw    = df['kw'].values  # seismic force
    V     = df['t'].values  # tension crack water force
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

    Fh = - kw + P * sin_b + R * cos_psi       # Equation (1)
    Fv = - W - P * cos_b + R * sin_psi        # Equation (2)
    Mo = - P * sin_b * (y_p - y_b) - P * cos_b * (x_p - x_b) \
        + kw * (y_k - y_b) - R * cos_psi * (y_r - y_b) + R * sin_psi * (x_r - x_b) # Equation (3)

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

    # First try Newton's method to find theta_opt
    try:
        theta_guess = 6.948  # This seems to work well for most cases
        if debug:
            print(f"Trying Newton's method for theta root-finding with initial guess {theta_guess:.1f} deg")
        theta_opt = newton(fs_difference, x0=theta_guess, tol=tol, maxiter=max_iter)
        # Check if the solution is valid
        if (
            abs(theta_opt) > 59 or
            abs(fs_difference(theta_opt)) > 0.01 or
            fs_force(np.radians(theta_opt)) >= fs_max - 1e-3
        ):
            raise ValueError("Newton's method converged to an invalid solution, trying sweep...")
    except Exception:
        # If Newton's method fails, try sweeping theta to find a bracket where fs_difference changes sign,
        #  then use root_scalar to find theta_opt
        if debug: 
            print("Newton's method failed or gave invalid result, sweeping theta to find a bracket...")
        theta_range = np.linspace(-60, 60, 49)
        fs_diff_values = [fs_difference(theta) for theta in theta_range]
        bracket = None
        for i in range(len(fs_diff_values) - 1):
            if fs_diff_values[i] * fs_diff_values[i+1] < 0:
                bracket = [theta_range[i], theta_range[i+1]]
                break
        if bracket is None:
            return False, "Spencer's method: No sign change found in fs_difference over theta range."
        try:
            if debug:
                print(f"Trying root_scalar with bracket {bracket}")
            sol = root_scalar(fs_difference, bracket=bracket, method='brentq', xtol=tol)
            theta_opt = sol.root
        except Exception as e:
            return False, f"Spencer's method failed to converge: {e}"

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
        yt_r[i] = y_b[i] - ((Mo[i] - Z[i] * sin_theta * dx[i] / 2 - Z[i+1] * sin_theta * dx[i] / 2 - Z[i] * cos_theta * (yt_l[i] - y_b[i])) / (Z[i+1] * cos_theta))  # Equation (30)
        if i < n-1:
            yt_l[i+1] = yt_r[i]
    df['yt_l'] = yt_l
    df['yt_r'] = yt_r

    print("Hey Norm!")
    
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
        if debug:
            for i in range(len(Q)):
                print(f"Slice {i+1}: Q = {Q[i]:.1f}, y_q = {y_q[i]:.2f}, dx = {dx[i]:.1f}, P = {P[i]:.1f}, sin_b = {sin_b[i]:.3f}, cos_b = {cos_b[i]:.3f}, Fh = {Fh[i]:.1f}, Fv = {Fv[i]:.1f}, Mo = {Mo[i]:.2f}, alpha = {np.degrees(alpha[i]):.2f}, beta = {np.degrees(beta[i]):.2f}, c = {c[i]:.2f}, w = {W[i]:.2f}, u = {u[i]:.2f}")

        return True, results