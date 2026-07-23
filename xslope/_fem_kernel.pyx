# cython: language_level=3
# cython: boundscheck=False
# cython: wraparound=False
# cython: cdivision=True
# cython: initializedcheck=False
"""
Compiled inner kernel for the SSRM viscoplastic iteration loop (fem.solve_fem).

This is a PROTOTYPE (Cython sandbox). It reproduces the Mohr-Coulomb Step-6
constitutive Gauss-point update -- INCLUDING the Rankine/tension cutoff (second
yield surface) and the matric-suction apparent-cohesion term -- with explicit
per-Gauss-point loops in place of the vectorised einsum / np.add.at reference.

It covers ONLY the standard Mohr-Coulomb material path. Power-curve and
Hoek-Brown Gauss-point groups (which re-linearise their envelope every
iteration) are NOT handled here; the caller keeps those groups on the NumPy
reference. 1D truss / pile corrections live outside the group loop and are
untouched.

The scatter into `loads` is performed in Gauss-point-major, dof-within-GP order,
identical to the reference `np.add.at(loads, dof.ravel(), contrib.ravel())`, so
the accumulation order is preserved. The only source of drift vs the reference
is floating-point re-association inside the per-GP matrix-vector products
(einsum/BLAS vs explicit C sums).

Reference: fem.solve_fem, the Step-6 block (commit 9090c46 vectorised it,
5f65f70 added prepared-model reuse).

Building (optional; the pure-NumPy path is the permanent default and oracle)
--------------------------------------------------------------------------
This module is NOT built by the normal `pip install`. It is opt-in and built
locally by anyone who has Cython installed:

    python setup_kernel.py build_ext --inplace

That compiles `xslope/_fem_kernel.<ext>.so` next to this source. `fem.solve_fem`
imports it only when called with `fast_kernel=True`; if it is not built it warns
and falls back to the NumPy reference. The generated `.c` and `.so` are build
artifacts (git-ignored) — only this `.pyx` is tracked.
"""

from libc.math cimport sqrt, sin, cos, asin, tan, fabs

# sqrt(3.0), hoisted (matches sq3 = np.sqrt(3.0) in the reference loop)
cdef double SQ3 = 1.7320508075688772935274463415058723

# Max supported dofs-per-element (quad8 = 16). tri3=6, quad4=8, tri6=12, quad8=16.
DEF MAX_NDOF = 32


def mc_step6(double[::1] u,
             double[::1] loads,
             double[:, :, ::1] B,
             double[:, :, ::1] D4,
             double[::1] w,
             Py_ssize_t[:, ::1] dof,
             double[::1] dt_r,
             double[::1] c_r,
             double[::1] c_suc,
             double[::1] snph,
             double[::1] csph,
             double[::1] t_cap,
             double[::1] u_gp,
             unsigned char[::1] elastic,
             double[:, ::1] evp,
             double dt,
             int has_cap,
             int has_elastic):
    """Run one iteration's Mohr-Coulomb Step-6 update for ONE Gauss-point group.

    Mutates `evp` (G,4) in place and scatters the viscoplastic body-load
    correction into `loads`. Returns the number of Mohr-Coulomb yielding Gauss
    points (n_yielding), matching the reference's np.count_nonzero(m).

    Shapes: B (G,3,ndof) tension-positive strain-displacement; D4 (G,4,4)
    4-component plane-strain constitutive; w (G,) integration weights; dof
    (G,ndof) global dof indices (intp); dt_r (G,) per-GP Rankine pseudo-timestep;
    c_r/snph/csph (G,) F-reduced cohesion / sin(phi_r) / cos(phi_r); c_suc (G,)
    F-reduced suction apparent cohesion (zeros if inactive); t_cap (G,) per-GP
    tensile cap (inf where none); u_gp (G,) effective-normal pore pressure added
    to sx,sy,sz (already zeroed under the 'effective' formulation); elastic (G,)
    uint8 pure-elastic flag.
    """
    cdef Py_ssize_t G = B.shape[0]
    cdef Py_ssize_t ndof = B.shape[2]
    cdef Py_ssize_t g, j
    cdef int n_yielding = 0

    cdef double ue[MAX_NDOF]

    cdef double e0, e1, e2
    cdef double ep0, ep1, ep2, ep3
    cdef double x0, x1, x2, x3
    cdef double s0, s1, s2, s3
    cdef double ug, sx, sy, txy, sz
    cdef double sigm, dsbar, dxv, dyv, dzv
    cdef double ds3, sine, theta, snth, csth
    cdef double c_env, f
    cdef double dsb, xj2, inv
    cdef double a2_0, a2_1, a2_2, a2_3
    cdef double a3_0, a3_1, a3_2, a3_3
    cdef double cs3, tn3, K, Kp, C2, C3
    cdef double flow_0, flow_1, flow_2, flow_3, fac
    cdef double cap, ctr, Rc, s1v, Rf, half, Ft
    cdef double sd0, sd1, sd2, contrib, wg
    cdef int corner, is_elastic

    for g in range(G):
        # ---- gather element displacement and compute strain eps = B u_e ----
        for j in range(ndof):
            ue[j] = u[dof[g, j]]
        e0 = 0.0
        e1 = 0.0
        e2 = 0.0
        for j in range(ndof):
            e0 += B[g, 0, j] * ue[j]
            e1 += B[g, 1, j] * ue[j]
            e2 += B[g, 2, j] * ue[j]

        # ---- effective stress sig4 = D4 (eps4 - evp) (tension-positive) ----
        ep0 = evp[g, 0]
        ep1 = evp[g, 1]
        ep2 = evp[g, 2]
        ep3 = evp[g, 3]
        x0 = e0 - ep0
        x1 = e1 - ep1
        x2 = e2 - ep2
        x3 = -ep3
        s0 = D4[g, 0, 0] * x0 + D4[g, 0, 1] * x1 + D4[g, 0, 2] * x2 + D4[g, 0, 3] * x3
        s1 = D4[g, 1, 0] * x0 + D4[g, 1, 1] * x1 + D4[g, 1, 2] * x2 + D4[g, 1, 3] * x3
        s2 = D4[g, 2, 0] * x0 + D4[g, 2, 1] * x1 + D4[g, 2, 2] * x2 + D4[g, 2, 3] * x3
        s3 = D4[g, 3, 0] * x0 + D4[g, 3, 1] * x1 + D4[g, 3, 2] * x2 + D4[g, 3, 3] * x3

        # sig_eff: add pore pressure to sx, sy, sz (indices 0,1,3); NOT txy (2)
        ug = u_gp[g]
        sx = s0 + ug
        sy = s1 + ug
        txy = s2
        sz = s3 + ug

        # ---- MC invariants ----
        sigm = (sx + sy + sz) / 3.0
        dsbar = sqrt(((sx - sy) * (sx - sy) + (sy - sz) * (sy - sz)
                      + (sz - sx) * (sz - sx) + 6.0 * txy * txy) / 2.0)
        dxv = sx - sigm
        dyv = sy - sigm
        dzv = sz - sigm

        if dsbar > 1e-10:
            ds3 = dsbar * dsbar * dsbar
            sine = -13.5 * (dxv * dyv * dzv - dzv * txy * txy) / ds3
        else:
            sine = 0.0
        if sine > 1.0:
            sine = 1.0
        elif sine < -1.0:
            sine = -1.0
        theta = asin(sine) / 3.0
        snth = sin(theta)
        csth = cos(theta)

        c_env = c_r[g] + c_suc[g]
        f = (sigm * snph[g]
             + dsbar * (csth / SQ3 - snth * snph[g] / 3.0)
             - c_env * csph[g])

        is_elastic = has_elastic and elastic[g]

        # ---- Mohr-Coulomb yield + viscoplastic flow (psi = 0, corner freeze) ----
        if (f > 0.0) and (dsbar > 1e-20) and (not is_elastic):
            n_yielding += 1
            dsb = dsbar
            xj2 = dsb * dsb / 3.0
            inv = 3.0 / (2.0 * dsb)
            a2_0 = inv * dxv
            a2_1 = inv * dyv
            a2_2 = inv * (2.0 * txy)
            a2_3 = inv * dzv
            a3_0 = dxv * dxv + txy * txy - (2.0 / 3.0) * xj2
            a3_1 = dyv * dyv + txy * txy - (2.0 / 3.0) * xj2
            a3_2 = -2.0 * dzv * txy
            a3_3 = dzv * dzv - (2.0 / 3.0) * xj2
            corner = fabs(snth) > 0.49
            if corner:
                C2 = 0.5
                C3 = 0.0
            else:
                cs3 = cos(3.0 * theta)
                tn3 = tan(3.0 * theta)
                K = csth / SQ3
                Kp = -snth / SQ3
                C2 = K - Kp * tn3
                C3 = -4.5 * Kp / (cs3 * dsb * dsb)
            flow_0 = C2 * a2_0 + C3 * a3_0
            flow_1 = C2 * a2_1 + C3 * a3_1
            flow_2 = C2 * a2_2 + C3 * a3_2
            flow_3 = C2 * a2_3 + C3 * a3_3
            fac = f * dt
            ep0 += fac * flow_0
            ep1 += fac * flow_1
            ep2 += fac * flow_2
            ep3 += fac * flow_3

        # ---- Rankine tension cutoff (second yield surface, Koiter sum) ----
        if has_cap:
            cap = t_cap[g]
            ctr = 0.5 * (sx + sy)
            Rc = sqrt((0.5 * (sx - sy)) * (0.5 * (sx - sy)) + txy * txy)
            s1v = ctr + Rc
            # s1v > cap is never true where cap == inf (capless GP)
            if (s1v > cap) and (not is_elastic):
                Rf = Rc if Rc > 1e-10 else 1e-10
                half = 0.5 * (sx - sy)
                Ft = (s1v - cap) * dt_r[g]
                ep0 += Ft * (0.5 + 0.5 * half / Rf)
                ep1 += Ft * (0.5 - 0.5 * half / Rf)
                ep2 += Ft * (txy / Rf)

        # ---- write back accumulated viscoplastic strain ----
        evp[g, 0] = ep0
        evp[g, 1] = ep1
        evp[g, 2] = ep2
        evp[g, 3] = ep3

        # ---- body-load correction: B^T (D4 evp)[:3] * w, scattered to dofs ----
        sd0 = D4[g, 0, 0] * ep0 + D4[g, 0, 1] * ep1 + D4[g, 0, 2] * ep2 + D4[g, 0, 3] * ep3
        sd1 = D4[g, 1, 0] * ep0 + D4[g, 1, 1] * ep1 + D4[g, 1, 2] * ep2 + D4[g, 1, 3] * ep3
        sd2 = D4[g, 2, 0] * ep0 + D4[g, 2, 1] * ep1 + D4[g, 2, 2] * ep2 + D4[g, 2, 3] * ep3
        wg = w[g]
        for j in range(ndof):
            contrib = (B[g, 0, j] * sd0 + B[g, 1, j] * sd1 + B[g, 2, j] * sd2) * wg
            loads[dof[g, j]] += contrib

    return n_yielding
