### ARCHIVE OF OLD METHODS - DOES NOT INCLUDE COMPLETE EQUATIONS ###



def oms_OLD(df, circle=None, circular=True, debug=True):
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
    c = df['c'].values
    phi = np.radians(df['phi']).values
    u = df['u'].values
    dl = df['dl'].values

    N_eff = W * cos_alpha - u * dl * cos2_alpha
    numerator = c * dl + N_eff * np.tan(phi)
    denominator = W * sin_alpha
    FS = numerator.sum() / denominator.sum() if denominator.sum() != 0 else float('inf')

    df['n_eff'] = N_eff  # store effective normal forces in df

    if debug==True:
        print(f'numerator = {numerator.sum():.4f}')
        print(f'denominator = {denominator.sum():.4f}')
        print('N_eff =', np.array2string(N_eff.values, precision=4, separator=', '))

    results = {}
    results['method'] = 'oms'
    results['FS'] = FS
    return True, results

def bishop_OLD(df, circle, circular=True, tol=1e-6, max_iter=100, debug=True):
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

    if circular == False:
        return False, 'Bishop method is only applicable for circular failure surfaces.'

    alpha_rad = np.radians(df['alpha'])
    cos_alpha = np.cos(alpha_rad)
    sin_alpha = np.sin(alpha_rad)
    cos2_alpha = cos_alpha**2
    tan_phi = np.tan(np.radians(df['phi']))

    W = df['w']
    c = df['c']
    dl = df['dl']
    u = df['u']

    # Right-hand side: sum of W * sin(alpha)
    denominator = (W * sin_alpha).sum()

    # Start iteration with an initial guess
    converge = False
    F_guess = 1.0
    N_eff = W * cos_alpha - u * dl * cos2_alpha
    num = c * dl + N_eff * tan_phi
    for _ in range(max_iter):

        num_N = (
            W - u * dl * cos_alpha
            - (c * dl * sin_alpha) / F_guess
        )
        denom_N = cos_alpha + (sin_alpha * tan_phi) / F_guess
        N_eff = num_N / denom_N

        numerator = (c * dl + N_eff * tan_phi).sum()
        F_calc = numerator / denominator
        if abs(F_calc - F_guess) < tol:
            converge = True
            break
        F_guess = F_calc

    df['n_eff'] = N_eff  # store effective normal forces in df

    if debug:
        print(f"FS = {F_calc:.6f}")
        print(f"Numerator = {numerator:.6f}")
        print(f"Denominator = {denominator:.6f}")
        print("N_eff =", np.array2string(N_eff.to_numpy(), precision=4, separator=', '))

    if not converge:
        return False, 'Bishop method did not converge within the maximum number of iterations.'
    else:
        results = {}
        results['method'] = 'bishop'
        results['FS'] = F_calc
        return True, results

def janbu_corrected_OLD(df, circle=None, circular=True, tol=1e-6, max_iter=100):
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



def force_equilibrium_OLD(df, theta_list, fs_guess=1.5, tol=1e-6, max_iter=50, debug=False):
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