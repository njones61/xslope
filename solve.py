
import numpy as np
from math import sin, cos, tan, radians, atan2, degrees

def oms(df):
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