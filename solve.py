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
            df['n_eff'] = N_eff
            if debug:
                print(f"FS = {F_new:.6f}")
                print(f"Numerator = {numerator:.6f}")
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

def spencer(df, tol=1e-4, max_iter = 100, debug_level=2):
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
        
        # Second derivatives of C3 and C4 (Equations 57-58)
        d2C3_dtheta2 = -cos_alpha_theta
        d2C4_dtheta2 = -sin_alpha_theta * tan_p
        
        # Second-order partial derivatives of Q (Equations 52-54)
        d2Q_dF2 = (1 / denom_Q**3) * (
            denom_Q * (2 * denom_Q * C2 / F**3 - 2 * (C1 + C2 / F) * C4 / F**3) -
            2 * C4 / F**2 * (denom_Q * C2 / F**2 - (C1 + C2 / F) * C4 / F**2)
        )
        
        d2Q_dFdtheta = (-1 / denom_Q**3) * (
            denom_Q * (C2 / F**2 * (dC3_dtheta + dC4_dtheta / F) - (C1 + C2 / F) * dC4_dtheta / F**2) -
            2 * (dC3_dtheta + dC4_dtheta / F) * (denom_Q * C2 / F**2 - (C1 + C2 / F) * C4 / F**2)
        )
        
        d2Q_dtheta2 = (-1 / denom_Q**3) * (
            denom_Q * (C1 + C2 / F) * (d2C3_dtheta2 + d2C4_dtheta2 / F) -
            2 * (C1 + C2 / F) * (dC3_dtheta + dC4_dtheta / F)**2
        )
        
        # Second-order partial derivatives of y_Q (Equations 61-63)
        dyQ_dF = (-1 / (Q * cos_theta)**2) * Mo * dQ_dF * cos_theta
        dyQ_dtheta = (-1 / (Q * cos_theta)**2) * Mo * (dQ_dtheta * cos_theta - Q * sin_theta)
        
        d2yQ_dF2 = (-1 / (Q**2 * cos_theta)) * Mo * (d2Q_dF2 - 2 / Q * dQ_dF**2)
        d2yQ_dFdtheta = (-1 / (Q**2 * cos_theta)) * Mo * (d2Q_dFdtheta + dQ_dF * np.tan(theta_rad) - 2 * dQ_dF * dQ_dtheta / Q)
        d2yQ_dtheta2 = (-1 / (Q**2 * cos_theta)) * Mo * (2 * d2Q_dtheta2 * np.tan(theta_rad) - d2Q_dtheta2 + Q + 2 / Q * (dQ_dtheta - Q * np.tan(theta_rad))**2)
        
        # Second-order partial derivatives of R1 (Equations 37-39)
        d2R1_dF2 = np.sum(d2Q_dF2)
        d2R1_dFdtheta = np.sum(d2Q_dFdtheta)
        d2R1_dtheta2 = np.sum(d2Q_dtheta2)
        
        # Second-order partial derivatives of R2 (Equations 42-44)
        d2R2_dF2 = np.sum(d2Q_dF2 * (x_b * sin_theta - y_q * cos_theta)) - 2 * np.sum(dQ_dF * dyQ_dF * cos_theta) - np.sum(Q * d2yQ_dF2 * cos_theta)
        d2R2_dFdtheta = np.sum(d2Q_dFdtheta * (x_b * sin_theta - y_q * cos_theta)) + np.sum(dQ_dF * (x_b * cos_theta + y_q * sin_theta - dyQ_dtheta * cos_theta)) - np.sum(dQ_dtheta * dyQ_dF * cos_theta) - np.sum(Q * (d2yQ_dFdtheta * cos_theta - dyQ_dF * sin_theta))
        d2R2_dtheta2 = np.sum(d2Q_dtheta2 * (x_b * sin_theta - y_q * cos_theta)) + 2 * np.sum(dQ_dtheta * (x_b * cos_theta + y_q * sin_theta - dyQ_dtheta * cos_theta)) - np.sum(Q * (x_b * sin_theta - y_q * cos_theta - 2 * dyQ_dtheta * sin_theta + d2yQ_dtheta2 * cos_theta))
        
        return dR1_dF, dR1_dtheta, dR2_dF, dR2_dtheta, d2R1_dF2, d2R1_dFdtheta, d2R1_dtheta2, d2R2_dF2, d2R2_dFdtheta, d2R2_dtheta2

    
    # Initial guesses
    F0 = 1.5
    theta0_rad = np.radians(-8.0) 
    
    # Newton iteration
    F = F0
    theta_rad = theta0_rad
    extended = False
    
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
        
        # Compute  derivatives
        dR1_dF, dR1_dtheta, dR2_dF, dR2_dtheta, d2R1_dF2, d2R1_dFdtheta, d2R1_dtheta2, d2R2_dF2, d2R2_dFdtheta, d2R2_dtheta2 = compute_derivatives(F, theta_rad, Q, y_q)
        
        # Basic Newton method (Equations 31-32)
        # denominator = dR1_dF * dR2_dtheta - dR1_dtheta * dR2_dF  # this seems wrong
        denominator = dR1_dtheta * dR2_dF - dR1_dF * dR2_dtheta
        
        if abs(denominator) < 1e-12:
            return False, "Singular Jacobian matrix in Newton iteration"
        
        delta_F = (R1 * dR2_dtheta - R2 * dR1_dtheta) / denominator
        delta_theta = (R2 * dR1_dF - R1 * dR2_dF) / denominator

        if debug_level >= 1 and not extended:
            print(f"          Basic: delta_F = {delta_F:.3f}, delta_theta = {np.degrees(delta_theta):.3f}°, {delta_theta: .3f} (rad)")
        
        # Check if we should switch to extended Newton method
        if abs(delta_F) < 0.5 and abs(delta_theta) < 0.15:

            if debug_level >= 1 and not extended:
                print(f"*** Switching to extended Newton method ***")
            extended = True
            
            # Build coefficient matrix for extended Newton
            A11 = dR1_dF + 0.5 * delta_F * d2R1_dF2 + 0.5 * delta_theta * d2R1_dFdtheta
            A12 = dR1_dtheta + 0.5 * delta_F * d2R1_dFdtheta + 0.5 * delta_theta * d2R1_dtheta2
            A21 = dR2_dF + 0.5 * delta_F * d2R2_dF2 + 0.5 * delta_theta * d2R2_dFdtheta
            A22 = dR2_dtheta + 0.5 * delta_F * d2R2_dFdtheta + 0.5 * delta_theta * d2R2_dtheta2
            
            # Solve extended Newton system
            A = np.array([[A11, A12], [A21, A22]])
            b = np.array([-R1, -R2])
            
            try:
                delta_solution = np.linalg.solve(A, b)
                delta_F = delta_solution[0]
                delta_theta = delta_solution[1]
            except np.linalg.LinAlgError:
                # Fall back to basic Newton if extended method fails
                pass

        if debug_level >= 1:
            print(f"          Extended: delta_F = {delta_F:.3f}, delta_theta = {np.degrees(delta_theta):.3f}°, {delta_theta: .3f} (rad)")

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
    

    ###### DEBUG RIGHT_FACING SLOPE ######
    print(f"F = {F:.12f}, theta = {np.degrees(theta_rad):.12f}°")
    # print("THETA_RAD CONVERTED TO 8.89552236670, F = 1.936856063097")
    # theta_rad = np.radians(8.895522366703)  # Hardwired to override the solution.
    # F = 1.936856063097          # Hardwired to override the solution.
    R1, R2, Q, y_q = compute_residuals(F, theta_rad)
    print(f"R1 = {R1:.6e}, R2 = {R2:.6e}")



    # Final computation of Q and y_q
    Q, y_q = compute_Q_and_yQ(F, theta_rad)

    if debug_level >= 2: 
        ma = 1 / (np.cos(alpha - theta_rad) + np.sin(alpha - theta_rad) * tan_p / F)
        df['ma'] = ma
        df['Q'] = Q
        df['y_q'] = y_q
        df['Fh'] = Fh
        df['Fv'] = Fv
        df['Mo'] = Mo
        # Debug print values per slice
        for i in range(len(Q)):
            print(f"Slice {i+1}: ma = {ma[i]:.3f}, Q = {Q[i]:.1f}, y_q = {y_q[i]:.2f}, Fh = {Fh[i]:.1f}, Fv = {Fv[i]:.1f}, Mo = {Mo[i]:.2f}")

    
    # Convert theta to degrees for output
    theta_opt = np.degrees(theta_rad)
    
    # ========== END SOLUTION ==========

    # Store theta in df
    df['theta'] = theta_opt

    # --- Compute N_eff using Equation (18) ---
    N_eff = - Fv * cos_a + Fh * sin_a + Q * np.sin(alpha - theta_rad) - u * dl
    df['n_eff'] = N_eff

    # --- Compute interslice forces Z using Equation (67) ---
    n = len(Q)
    Z = np.zeros(n+1)
    for i in range(n):
        Z[i+1] = Z[i] - Q[i] 
    df['z'] = Z[:-1]        # Z_i acting on slice i's left face
 

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
    df['yt_l'] = yt_l
    df['yt_r'] = yt_r
    
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




def rapid_drawdown(df, method_func, debug_level=1):
    """
    Performs rapid drawdown analysis using a three-stage approach.
    
    Parameters:
        df : pandas.DataFrame
            Slice data with all required columns including rapid drawdown specific data:
            - c, phi: current strength parameters
            - c1, phi1: original strength parameters (for stage 3)
            - d, psi: rapid drawdown parameters for low-K materials
            - u: pore pressure (stage 1)
            - u2: pore pressure for lowered pool (stage 2)
            - dload, d_x, d_y: distributed loads (stage 1)
            - dload2, d_x2, d_y2: distributed loads for lowered pool (stage 2)
        method_func : function
            The method function to use (oms, bishop, spencer, etc.)
        debug_level : int
            0: no output, 1: print FS at each stage, >1: detailed debug info
    
    Returns:
        Tuple(bool, dict): (True, result_dict) or (False, error_message)
    """
    
    if debug_level >= 1:
        print("=== RAPID DRAWDOWN ANALYSIS ===")
    
    # Stage 1: Pre-drawdown conditions
    if debug_level >= 1:
        print("Stage 1: Pre-drawdown conditions...")
    
    # Use original conditions (c, phi, u, dload, d_x, d_y)
    success, result_stage1 = method_func(df)
    if not success:
        return False, f"Stage 1 failed: {result_stage1}"
    
    stage1_FS = result_stage1['FS']
    if debug_level >= 1:
        print(f"Stage 1 FS = {stage1_FS:.4f}")
    
    # Calculate consolidation stresses for each slice
    # N_eff should be available from the method function
    if 'n_eff' not in df.columns:
        return False, "Stage 1 did not compute n_eff values"
    
    # Calculate sigma_fc and tau_fc for each slice
    sigma_fc = df['n_eff'] / df['dl']  # Equation (2)
    tau_fc = (1.0 / stage1_FS) * (df['c'] + sigma_fc * np.tan(np.radians(df['phi'])))  # Equation (3)
    
    if debug_level >= 2:
        print("Stage 1 consolidation stresses:")
        for i in range(len(df)):
            print(f"  Slice {i+1}: sigma_fc = {sigma_fc.iloc[i]:.2f}, tau_fc = {tau_fc.iloc[i]:.2f}")
    
    # Stage 2: Post-drawdown conditions with undrained strengths
    if debug_level >= 1:
        print("Stage 2: Post-drawdown conditions with undrained strengths...")
   
    # Update pore pressures and distributed loads for stage 2
    df['u'] = df['u2']
    df['dload'] = df['dload2']
    df['d_x'] = df['d_x2']
    df['d_y'] = df['d_y2']
    
    # Process each slice for undrained strength calculation
    for i in range(len(df)):
        # Check if this slice has low-K material (d and psi are not zero)
        d_val = df.iloc[i]['d']
        psi_val = df.iloc[i]['psi']
        
        if d_val > 0 and psi_val > 0:
            # Low-K material - calculate undrained strength
            if debug_level >= 2:
                print(f"Processing low-K material for slice {i+1}")
            
            # Get consolidation stresses for this slice
            sigma_fc_i = sigma_fc.iloc[i]
            tau_fc_i = tau_fc.iloc[i]
            phi_deg = df.iloc[i]['phi1']  # Use original phi for calculations
            c_val = df.iloc[i]['c1']      # Use original c for calculations
            
            # Calculate K1 using equation (4)
            phi_rad = np.radians(phi_deg)
            if abs(np.cos(phi_rad)) < 1e-12:
                if debug_level >= 2:
                    print(f"  Warning: cos(phi) near zero for slice {i+1}, skipping K1 calculation")
                continue
            
            K1 = (sigma_fc_i + tau_fc_i * (np.sin(phi_rad) + 1) / np.cos(phi_rad)) / \
                 (sigma_fc_i + tau_fc_i * (np.sin(phi_rad) - 1) / np.cos(phi_rad))
            
            if debug_level >= 2:
                print(f"  K1 = {K1:.4f}")
            
            # Calculate Kf using equation (6)
            if abs(sigma_fc_i - c_val * np.cos(phi_rad)) < 1e-12:
                if debug_level >= 2:
                    print(f"  Warning: denominator near zero for Kf calculation in slice {i+1}")
                continue
            
            Kf = ((sigma_fc_i + c_val * np.cos(phi_rad)) * (1 + np.sin(phi_rad))) / \
                 ((sigma_fc_i - c_val * np.cos(phi_rad)) * (1 - np.sin(phi_rad)))
            
            if debug_level >= 2:
                print(f"  Kf = {Kf:.4f}")
            
            # Check for negative stresses using equations (7) and (8)
            sigma3_k1 = sigma_fc_i + tau_fc_i * (np.sin(phi_rad) - 1) / np.cos(phi_rad)  # Equation (7)
            sigma3_kf = (sigma_fc_i - c_val * np.cos(phi_rad)) * (1 - np.sin(phi_rad)) / (np.cos(phi_rad)**2)  # Equation (8)
            
            if debug_level >= 2:
                print(f"  sigma3_k1 = {sigma3_k1:.4f}, sigma3_kf = {sigma3_kf:.4f}")
            
            # Calculate tau_ff values for both curves
            tau_ff_k1 = d_val + sigma_fc_i * np.tan(np.radians(psi_val))  # d-psi curve
            tau_ff_kf = c_val + sigma_fc_i * np.tan(phi_rad)  # c-phi curve
            
            if debug_level >= 2:
                print(f"  tau_ff_k1 = {tau_ff_k1:.4f}, tau_ff_kf = {tau_ff_kf:.4f}")
            
            # Determine which tau_ff to use
            if sigma3_k1 < 0 or sigma3_kf < 0:
                # Use the lower of the two curves
                tau_ff = min(tau_ff_k1, tau_ff_kf)
                if debug_level >= 2:
                    print(f"  Negative stress detected, using lower curve: tau_ff = {tau_ff:.4f}")
            else:
                # Interpolate using equation (5)
                if abs(Kf - 1) < 1e-12:
                    tau_ff = tau_ff_k1
                else:
                    tau_ff = ((Kf - K1) * tau_ff_k1 + (K1 - 1) * tau_ff_kf) / (Kf - 1)
                
                if debug_level >= 2:
                    print(f"  Interpolated tau_ff = {tau_ff:.4f}")
            
            # Set undrained strength parameters
            df.iloc[i, df.columns.get_loc('c')] = float(tau_ff)
            df.iloc[i, df.columns.get_loc('phi')] = 0.0
            
            if debug_level >= 2:
                print(f"  Set c = {tau_ff:.4f}, phi = 0.0 for slice {i+1}")
        else:
            # High-K material - keep original c and phi
            if debug_level >= 2:
                print(f"Slice {i+1}: High-K material, keeping original c and phi")
    
    # Calculate Stage 2 FS
    success, result_stage2 = method_func(df)
    if not success:
        return False, f"Stage 2 failed: {result_stage2}"
    
    stage2_FS = result_stage2['FS']
    if debug_level >= 1:
        print(f"Stage 2 FS = {stage2_FS:.4f}")
    
    # Stage 3: Check drained strengths
    if debug_level >= 1:
        print("Stage 3: Checking drained strengths...")
    
    # Check if any low-K slices need drained strength
    need_stage3 = False
    
    for i in range(len(df)):
        d_val = df.iloc[i]['d']
        psi_val = df.iloc[i]['psi']
        
        if d_val > 0 and psi_val > 0:
            # This is a low-K material slice
            if 'n_eff' not in df.columns:
                return False, "Stage 2 did not compute n_eff values"
            
            # Calculate drained strength using equations (9) and (10)
            sigma_prime = df.iloc[i]['n_eff'] / df.iloc[i]['dl']  # Equation (9)
            tau_drained = df.iloc[i]['c1'] + sigma_prime * np.tan(np.radians(df.iloc[i]['phi1']))  # Equation (10)
            
            # Compare with undrained strength (current c value)
            tau_undrained = df.iloc[i]['c']
            
            if debug_level >= 2:
                print(f"Slice {i+1}: tau_drained = {tau_drained:.4f}, tau_undrained = {tau_undrained:.4f}")
            
            if tau_drained < tau_undrained:
                # Use drained strength
                df.iloc[i, df.columns.get_loc('c')] = float(df.iloc[i]['c1'])
                df.iloc[i, df.columns.get_loc('phi')] = float(df.iloc[i]['phi1'])
                need_stage3 = True
                
                if debug_level >= 2:
                    print(f"  Using drained strength for slice {i+1}")
    
    if need_stage3:
        if debug_level >= 1:
            print("Stage 3: Recalculating FS with drained strengths...")
        
        success, result_stage3 = method_func(df)
        if not success:
            return False, f"Stage 3 failed: {result_stage3}"
        
        stage3_FS = result_stage3['FS']
        if debug_level >= 1:
            print(f"Stage 3 FS = {stage3_FS:.4f}")
    else:
        stage3_FS = stage2_FS
        if debug_level >= 1:
            print("Stage 3: No drained strength adjustments needed")
    
    # Final FS is the lower of Stage 2 and Stage 3
    if stage2_FS < stage3_FS:
        final_FS = stage2_FS
        result = result_stage2
    else:
        final_FS = stage3_FS
        result = result_stage3
    
    if debug_level >= 1:
        print(f"Final rapid drawdown FS = {final_FS:.4f}")
        print("=== END RAPID DRAWDOWN ANALYSIS ===")
    
    # Append stage FS to result
    result['stage1_FS'] = stage1_FS
    result['stage2_FS'] = stage2_FS
    result['stage3_FS'] = stage3_FS

    return True, result