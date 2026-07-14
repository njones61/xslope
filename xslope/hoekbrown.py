"""Generalized Hoek-Brown strength for slope stability (template v14, option 'hb').

Hoek-Brown is written in PRINCIPAL stresses:

    sigma1' = sigma3' + sigma_ci * (mb * sigma3'/sigma_ci + s)^a

but limit equilibrium and Mohr-Coulomb plasticity both need shear strength as a
function of the NORMAL stress on the failure plane, tau(sigma_n'). The two are
connected by Balmer's (1952) transformation, which gives the point on the Mohr
envelope that is tangent to the Mohr circle through (sigma3', sigma1'):

    dsigma1/dsigma3 = 1 + a * mb * (mb*sigma3'/sigma_ci + s)^(a-1)
    sigma_n' = sigma3' + (sigma1' - sigma3') / (dsigma1/dsigma3 + 1)
    tau      = (sigma_n' - sigma3') * sqrt(dsigma1/dsigma3)

That is a PARAMETRIC curve in sigma3', so getting tau at a prescribed sigma_n'
means inverting sigma_n'(sigma3') numerically. sigma_n' is strictly increasing in
sigma3', so a bracketed bisection is robust and needs no derivatives.

The instantaneous friction angle follows analytically:

    tan(phi_i) = dtau/dsigma_n = (dsigma1/dsigma3 - 1) / (2 * sqrt(dsigma1/dsigma3))

and the instantaneous cohesion is the intercept c_i = tau - sigma_n * tan(phi_i).
This (c_i, phi_i) tangent is exactly what the rest of XSLOPE consumes, so
Hoek-Brown rides the same nonlinear-strength machinery as the power curve: the
LEM re-linearizes at the converged base normal stress (solve._with_nonlinear_
strength) and the FEM re-linearizes each viscoplastic iteration.

Strength reduction: dividing the SHEAR strength by F divides both c_i and
tan(phi_i) by F, so SSR needs no special handling — reduce the tangent, not the
Hoek-Brown constants (sigma_ci/F is NOT the same thing, because of the exponent).
"""

import numpy as np


def hb_constants(gsi, mi, d):
    """Rock-mass constants (mb, s, a) from GSI, the intact constant mi, and the
    disturbance factor D — Hoek, Carranza-Torres & Corkum (2002). These are
    DERIVED, never entered: nobody has mb/s/a to hand, but GSI/mi/D come off a
    field log."""
    gsi = np.asarray(gsi, dtype=float)
    mi = np.asarray(mi, dtype=float)
    d = np.clip(np.asarray(d, dtype=float), 0.0, 1.0)
    mb = mi * np.exp((gsi - 100.0) / (28.0 - 14.0 * d))
    s = np.exp((gsi - 100.0) / (9.0 - 3.0 * d))
    a = 0.5 + (np.exp(-gsi / 15.0) - np.exp(-20.0 / 3.0)) / 6.0
    return mb, s, a


def hb_sigma_t(sci, mb, s):
    """Hoek-Brown tensile strength: the apex, where sigma1' = sigma3'. sigma3'
    can never go below this."""
    with np.errstate(divide='ignore', invalid='ignore'):
        st = -s * sci / mb
    return np.where(np.asarray(mb) > 0, st, 0.0)


def _balmer(sig3, sci, mb, s, a):
    """(sigma_n, tau, dsig1_dsig3) at a given minor principal stress."""
    sig3 = np.asarray(sig3, dtype=float)
    base = np.maximum(mb * sig3 / sci + s, 0.0)
    sig1 = sig3 + sci * base ** a
    # d(sigma1)/d(sigma3); >= 1 always, so the sqrt below is safe
    dsig = 1.0 + a * mb * base ** (a - 1.0)
    dsig = np.maximum(dsig, 1.0)
    sig_n = sig3 + (sig1 - sig3) / (dsig + 1.0)
    tau = (sig_n - sig3) * np.sqrt(dsig)
    return sig_n, tau, dsig


def hb_tangent(sigma_n, sci, gsi, mi, d, iters=60):
    """Instantaneous Mohr-Coulomb tangent (c_i, phi_i in DEGREES) of the
    generalized Hoek-Brown envelope at normal stress ``sigma_n``.

    Inverts Balmer's parametric curve by bisection on sigma3' (sigma_n is
    monotone increasing in sigma3', so the bracket is guaranteed). Vectorized
    over sigma_n; the rock constants broadcast.
    """
    mb, s, a = hb_constants(gsi, mi, d)
    return hb_tangent_const(sigma_n, sci, mb, s, a, iters=iters)




def hb_tangent_const(sigma_n, sci, mb, s, a, iters=60):
    """hb_tangent with the rock constants already derived — for hot loops (the
    FEM re-linearizes at every Gauss point on every viscoplastic iteration, and
    re-deriving mb/s/a there would repeat the same exponentials thousands of
    times).

    THE sigma_n >= 0 CLAMP BELOW IS LOAD-BEARING, not defensive tidying. It bounds
    the search away from the Hoek-Brown tensile apex, where the envelope turns
    vertical: dsigma1/dsigma3 diverges there, so tan(phi_i) -> infinity and
    phi_i -> 90 degrees (NaN exactly at the apex). Gauss points next to a slope face
    routinely sit in tension, and letting them reach the apex hands them a
    near-vertical envelope — effectively infinite strength. Clamping at zero instead
    bounds phi_i at its zero-confinement value (~60 degrees for a typical rock mass):
    steep, but finite and physical.
    """
    sigma_n = np.maximum(np.asarray(sigma_n, dtype=float), 0.0)
    sci = np.asarray(sci, dtype=float)

    lo = np.broadcast_to(hb_sigma_t(sci, mb, s), sigma_n.shape).astype(float).copy()
    # Upper bracket: sigma_n(sigma3) >= sigma3, so sigma3 = sigma_n always
    # over-shoots. Pad so the root is strictly interior.
    hi = np.maximum(sigma_n, lo + 1e-9) + 1.0

    for _ in range(iters):
        mid = 0.5 * (lo + hi)
        sn, _tau, _dd = _balmer(mid, sci, mb, s, a)
        too_high = sn > sigma_n
        hi = np.where(too_high, mid, hi)
        lo = np.where(too_high, lo, mid)

    sig3 = 0.5 * (lo + hi)
    _sn, tau, dsig = _balmer(sig3, sci, mb, s, a)

    tan_phi = (dsig - 1.0) / (2.0 * np.sqrt(dsig))
    phi = np.degrees(np.arctan(tan_phi))
    c = tau - sigma_n * tan_phi
    return np.maximum(c, 0.0), phi
