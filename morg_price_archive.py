import numpy as np
from math import sin, cos, tan, radians, atan2, degrees
from scipy.optimize import minimize_scalar, root_scalar

#######################################################################################
## None of these work. I put them here as an archive in case we want to use them later.
#######################################################################################

def morgenstern_price(df, psi_function=lambda x: 1.0, tol=1e-3, max_iter=100):
    """
    Full Morgenstern-Price Method with slice-wise recursion.
    Enforces both force and moment equilibrium.

    Parameters:
        df (pd.DataFrame): Slice data with required fields.
        psi_function (callable): Defines shape of interslice force function ψ(x)
        tol (float): Convergence tolerance
        max_iter (int): Max number of iterations

    Returns:
        float: Factor of Safety
        float: Lambda (scaling factor for interslice function)
        bool: Converged flag
    """
    alpha = np.radians(df['alpha'].values)
    phi = np.radians(df['phi'].values)
    c = df['c'].values
    dl = df['dl'].values
    w = df['w'].values
    u = df['u'].values

    n = len(df)
    x_norm = np.linspace(0, 1, n)
    psi = np.array([psi_function(xi) for xi in x_norm])

    sin_alpha = np.sin(alpha)
    cos_alpha = np.cos(alpha)
    tan_phi = np.tan(phi)

    F = 1.0
    lam = 0.0
    converged = False

    for _ in range(max_iter):
        E = np.zeros(n + 1)  # interslice shear forces
        X = np.zeros(n + 1)  # interslice normal forces
        N_base = np.zeros(n)
        T_base = np.zeros(n)

        # Forward recursion for X and E
        for i in range(n):
            N = w[i] * cos_alpha[i] - u[i] * dl[i] * cos_alpha[i] + X[i] * psi[i] * sin_alpha[i]
            S = c[i] * dl[i] + (N) * tan_phi[i] / F
            R = S  # resisting shear
            T = w[i] * sin_alpha[i] + E[i]
            X[i + 1] = X[i] + lam * psi[i] * (T - R) * cos_alpha[i] / (1 + lam * psi[i] * tan_phi[i] / F)
            E[i + 1] = E[i] + lam * psi[i] * (T - R) * sin_alpha[i] / (1 + lam * psi[i] * tan_phi[i] / F)

            N_base[i] = N
            T_base[i] = T

        # Recalculate F using total resistance and driving
        S_total = c * dl + N_base * tan_phi / F
        R_total = S_total
        T_total = w * sin_alpha + E[:-1]
        F_new = R_total.sum() / T_total.sum()

        # Recalculate lambda from force direction matching
        num = (R_total * psi * sin_alpha).sum()
        den = (R_total * psi * cos_alpha).sum()
        lam_new = num / den if den != 0 else lam

        if abs(F_new - F) < tol and abs(lam_new - lam) < tol:
            converged = True
            break

        F = F_new
        lam = lam_new

        # F = 0.7 * F + 0.3 * F_new
        # lam_new = np.clip(lam_new, lam - 0.05, lam + 0.05)
        # lam = lam_new

    return F, lam, converged


def morgenstern_price_full(df, psi_function=lambda x: np.sin(np.pi * x), tol=1e-6, max_iter=100):
    """
    Full Morgenstern-Price Method with slice-wise recursion.
    Enforces both force and moment equilibrium.

    Parameters:
        df (pd.DataFrame): Slice data with required fields.
        psi_function (callable): Defines shape of interslice force function ψ(x)
        tol (float): Convergence tolerance
        max_iter (int): Max number of iterations

    Returns:
        float: Factor of Safety
        float: Lambda (scaling factor for interslice function)
        bool: Converged flag
    """
    alpha = np.radians(df['alpha'].values)
    phi = np.radians(df['phi'].values)
    c = df['c'].values
    dl = df['dl'].values
    w = df['w'].values
    u = df['u'].values

    n = len(df)
    x_norm = np.linspace(0, 1, n)
    psi = np.array([psi_function(xi) for xi in x_norm])

    sin_alpha = np.sin(alpha)
    cos_alpha = np.cos(alpha)
    tan_phi = np.tan(phi)

    F = 1.0
    lam = 0.0
    converged = False

    fs_history = []
    lam_history = []

    for _ in range(max_iter):
        E = np.zeros(n + 1)  # interslice shear forces
        X = np.zeros(n + 1)  # interslice normal forces
        N_base = np.zeros(n)
        T_base = np.zeros(n)

        # Forward recursion for X and E
        for i in range(n):
            N = w[i] * cos_alpha[i] - u[i] * dl[i] * cos_alpha[i] + X[i] * psi[i] * sin_alpha[i]
            S = c[i] * dl[i] + (N) * tan_phi[i] / F
            R = S  # resisting shear
            T = w[i] * sin_alpha[i]
            X[i + 1] = X[i] + lam * psi[i] * (T - R) * cos_alpha[i] / (1 + lam * psi[i] * tan_phi[i] / F)
            E[i + 1] = E[i] + lam * psi[i] * (T - R) * sin_alpha[i] / (1 + lam * psi[i] * tan_phi[i] / F)
            N_base[i] = N
            T_base[i] = T

        # Recalculate F using total resistance and driving
        S_total = c * dl + N_base * tan_phi / F
        R_total = S_total
        T_total = w * sin_alpha + E[:-1]
        F_new = R_total.sum() / T_total.sum()

        # Recalculate lambda from force direction matching
        num = (R_total * psi * sin_alpha).sum()
        den = (R_total * psi * cos_alpha).sum()
        lam_new = num / den if den != 0 else lam

        if abs(F_new - F) < tol and abs(lam_new - lam) < tol:
            converged = True
            break

        fs_history.append(F_new)
        lam_history.append(lam_new)
        F = 0.7 * F + 0.3 * F_new
        lam_new = np.clip(lam_new, lam - 0.05, lam + 0.05)
        lam = lam_new

    import matplotlib.pyplot as plt
    plt.figure(figsize=(8, 4))
    plt.plot(fs_history, label='FS')
    plt.plot(lam_history, label='Lambda')
    plt.xlabel('Iteration')
    plt.ylabel('Value')
    plt.title('Morgenstern-Price Convergence History')
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()

    return F, lam, converged

def morgenstern_price_system(df, psi_function=lambda x: 1.0, tol=1e-6, max_iter=100):
    """
    Full Morgenstern-Price method using system-of-equations approach.
    Solves global force and moment equilibrium using slice-by-slice interaction.

    Parameters:
        df (pd.DataFrame): Slice data with columns 'alpha', 'phi', 'c', 'w', 'u', 'dl'
        psi_function (callable): ψ(x) function for interslice force distribution
        tol (float): Convergence tolerance
        max_iter (int): Maximum iterations

    Returns:
        float: Factor of Safety
        float: lambda (interslice scaling factor)
        bool: convergence status
    """
    alpha = np.radians(df['alpha'].values)
    phi = np.radians(df['phi'].values)
    c = df['c'].values
    l = df['dl'].values
    w = df['w'].values
    u = df['u'].values

    sin_alpha = np.sin(alpha)
    cos_alpha = np.cos(alpha)
    tan_phi = np.tan(phi)

    n = len(df)
    x_norm = np.linspace(0, 1, n)
    psi = np.array([psi_function(xi) for xi in x_norm])

    def residual_system(F):
        # Construct tridiagonal system Ax = b for interslice forces
        A = np.zeros((n, n))
        b = np.zeros(n)

        for i in range(n):
            m_i = 1 + psi[i] * tan_phi[i] / F
            k_i = psi[i] * tan_phi[i] / F
            N = w[i] * cos_alpha[i] - u[i] * l[i]
            S = c[i] * l[i] + N * tan_phi[i] / F
            T = w[i] * sin_alpha[i]

            A[i, i] = m_i
            if i > 0:
                A[i, i - 1] = -k_i

            b[i] = T - S

        try:
            X = np.linalg.solve(A, b)
        except np.linalg.LinAlgError:
            return 1e6  # non-invertible system

        residual = np.sum(X * psi * sin_alpha)  # net shear
        return abs(residual)

    fa = residual_system(0.01)
    fb = residual_system(10.0)
    print("Residual at FS = 0.01:", fa)
    print("Residual at FS = 10.0:", fb)

    if fa * fb > 0:
        print("❌ Cannot solve: residual does not change sign over the bracket.")
        return None, None, False

    result = root_scalar(residual_system, bracket=[0.01, 10.0], method='brentq', xtol=tol)
    F_opt = result.root
    converged = result.converged
    lam_opt = 1.0  # in this simplified version, lambda is not iterated explicitly

    return F_opt, lam_opt, converged