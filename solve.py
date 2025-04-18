
import numpy as np
from math import sin, cos, tan, radians, atan2, degrees

def oms(df):
    """
    Computes the Factor of Safety (FS) using the Ordinary Method of Slices (OMS).

    Parameters:
        df (pd.DataFrame): DataFrame containing slice data with columns:
            - 'alpha': base angle (degrees)
            - 'w': total slice weight
            - 'c': cohesion
            - 'phi': friction angle (degrees)
            - 'u': pore pressure
            - 'dl': base length

    Returns:
        tuple:
            - float: Computed Factor of Safety (FS)
            - np.ndarray: Normal force on the base of each slice
    """

    alpha_rad = np.radians(df['alpha'])

    cos_alpha = np.cos(alpha_rad)
    sin_alpha = np.sin(alpha_rad)
    cos2_alpha = cos_alpha ** 2

    W = df['w']
    c = df['c']
    phi = np.radians(df['phi'])
    u = df['u']
    dl = df['dl']

    N = W * cos_alpha - u * dl * cos2_alpha
    numerator = c * dl + N * np.tan(phi)
    denominator = W * sin_alpha
    FS = numerator.sum() / denominator.sum() if denominator.sum() != 0 else float('inf')

    return FS, N


def bishops(df, tol=1e-6, max_iter=100):
    """
    Computes the Factor of Safety (FS) using Bishop's Simplified Method.

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

    return F_calc, N, converge



def spencer(df, tol=1e-6, max_iter=100):
    """
    Computes the Factor of Safety (FS) using Spencer's Method.

    Parameters:
        df (pd.DataFrame): DataFrame containing slice data.
        tol (float, optional): Convergence tolerance. Default is 1e-6.
        max_iter (int, optional): Maximum number of iterations. Default is 100.

    Returns:
        tuple:
            - float: Computed Factor of Safety (FS)
            - float: Inter-slice force inclination angle (degrees)
    """

    alpha_rad = np.radians(df['alpha'])
    phi_rad = np.radians(df['phi'])
    tan_phi = np.tan(phi_rad)

    W = df['w'].values
    c = df['c'].values
    dl = df['dl'].values
    u = df['u'].values

    sin_alpha = np.sin(alpha_rad)
    cos_alpha = np.cos(alpha_rad)
    cos2_alpha = cos_alpha ** 2

    # Initial guesses
    F = 1.0
    beta = 0.0  # in radians

    for _ in range(max_iter):
        sin_beta = np.sin(beta)
        cos_beta = np.cos(beta)

        # Compute resisting force per slice
        denom = cos_alpha * cos_beta + sin_alpha * sin_beta * np.tan(phi_rad)
        num = c * dl * cos_beta + (W - u * dl) * (cos_alpha * cos_beta - sin_alpha * sin_beta) * np.tan(phi_rad)
        F_new = num.sum() / (W * sin_alpha).sum()

        # Update beta using ratio of residual horizontal and vertical forces
        sin_beta_new = (W * sin_alpha - num / F_new * sin_alpha * tan_phi) / W
        beta_new = np.arcsin(np.clip(sin_beta_new.mean(), -1.0, 1.0))  # averaged beta

        if abs(F_new - F) < tol and abs(beta_new - beta) < tol:
            return F_new, degrees(beta_new)

        F = F_new
        beta = beta_new

    return F, degrees(beta)


def janbu_simple(df, tol=1e-6, max_iter=100):
    """
    Computes the Factor of Safety (FS) using Janbu's Simplified Method.

    Parameters:
        df (pd.DataFrame): DataFrame containing slice data.
        tol (float, optional): Convergence tolerance. Default is 1e-6.
        max_iter (int, optional): Maximum number of iterations. Default is 100.

    Returns:
        float: Computed Factor of Safety (FS)
    """

    alpha_rad = np.radians(df['alpha'])
    phi_rad = np.radians(df['phi'])
    tan_phi = np.tan(phi_rad)

    W = df['w'].values
    c = df['c'].values
    dl = df['dl'].values
    u = df['u'].values

    sin_alpha = np.sin(alpha_rad)
    cos_alpha = np.cos(alpha_rad)

    F = 1.0  # initial guess

    for _ in range(max_iter):
        N = W * cos_alpha
        S = c * dl + (N - u * dl) * tan_phi / F
        T = W * sin_alpha

        F_new = S.sum() / T.sum()

        if abs(F_new - F) < tol:
            return F_new

        F = F_new

    return F


def janbu_corrected(df, tol=1e-6, max_iter=100):
    """
    Computes the Factor of Safety (FS) using Janbu's Corrected Method.

    Parameters:
        df (pd.DataFrame): DataFrame containing slice data.
        tol (float, optional): Convergence tolerance. Default is 1e-6.
        max_iter (int, optional): Maximum number of iterations. Default is 100.

    Returns:
        tuple:
            - float: Computed Factor of Safety (FS)
            - float: Horizontal force ratio (lambda)
            - bool: Whether the solution converged
    """

    alpha_rad = np.radians(df['alpha'])
    phi_rad = np.radians(df['phi'])
    tan_phi = np.tan(phi_rad)

    W = df['w'].values
    c = df['c'].values
    dl = df['dl'].values
    u = df['u'].values

    sin_alpha = np.sin(alpha_rad)
    cos_alpha = np.cos(alpha_rad)

    F = 1.0
    lambda_ = 0.0
    converged = False

    for _ in range(max_iter):
        m = 1 + lambda_ * tan_phi / F

        N = W * cos_alpha
        S = c * dl + (N - u * dl) * tan_phi / F
        R = S / m
        T = W * sin_alpha

        F_new = R.sum() / T.sum()
        lambda_new = (R * sin_alpha).sum() / (R * cos_alpha).sum()

        if abs(F_new - F) < tol and abs(lambda_new - lambda_) < tol:
            converged = True
            break

        F = F_new
        lambda_ = lambda_new

    return F, lambda_, converged


def morgenstern_price(df, function=lambda x: 1.0, tol=1e-6, max_iter=100):
    """
    Computes the Factor of Safety (FS) using the Morgenstern-Price Method.

    Parameters:
        df (pd.DataFrame): DataFrame containing slice data.
        function (callable, optional): A callable defining the inter-slice force function ψ(x),
            where x is normalized from 0 to 1 across slices. Default is constant function (1.0).
        tol (float, optional): Convergence tolerance. Default is 1e-6.
        max_iter (int, optional): Maximum number of iterations. Default is 100.

    Returns:
        tuple:
            - float: Computed Factor of Safety (FS)
            - float: Inter-slice force ratio (lambda)
            - bool: Whether the solution converged
    """
    alpha_rad = np.radians(df['alpha'])
    phi_rad = np.radians(df['phi'])
    tan_phi = np.tan(phi_rad)

    W = df['w'].values
    c = df['c'].values
    dl = df['dl'].values
    u = df['u'].values

    sin_alpha = np.sin(alpha_rad)
    cos_alpha = np.cos(alpha_rad)

    # Inter-slice function shape (normalized to [0, 1])
    n = len(df)
    x_norm = np.linspace(0, 1, n)
    psi = np.array([function(xi) for xi in x_norm])

    F = 1.0
    lam = 0.0
    converged = False

    for _ in range(max_iter):
        m = 1 + lam * psi * tan_phi / F

        N = W * cos_alpha
        S = c * dl + (N - u * dl) * tan_phi / F
        R = S / m
        T = W * sin_alpha

        F_new = R.sum() / T.sum()
        num = (R * psi * sin_alpha).sum()
        den = (R * psi * cos_alpha).sum()
        lam_new = num / den if den != 0 else 0.0

        if abs(F_new - F) < tol and abs(lam_new - lam) < tol:
            converged = True
            break

        F = F_new
        lam = lam_new

    return F, lam, converged