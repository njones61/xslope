"""A reported yielded element is one whose stress is at its strength.

``plastic_elements`` (the element table's ``plastic`` column, ``plastic_fraction``,
the yielded-element counts and the outlines on the stress plot) was read as
``f > 1e-8`` on the element-averaged stress on the viscoplastic path, and as "a
Gauss point on a plastic branch of the return map" on the Newton corrector's.
Both miss the points the solve holds at their strength: the viscoplastic iteration
relaxes an over-stressed point back onto the yield surface, so it ends a hair
either side of f = 0, and an average of points on a convex surface lies inside
it; the corrector, started from an equilibrated seed, returns most of them on the
elastic branch. On the two-layer clay sample (tri6 at 2 ft, phi = 0) the flag read
0%, 8.5% and 19.5% at F = 1.0, 1.125 and 1.219, where 13%, 22% and 40% of the
elements carry a stress level of 0.97 or more. The stress plot read
``yield_function > 0`` instead and outlined nothing at all.

An element is now flagged when one of its Gauss points has a stress level
SL = tau / tau_f of at least 1 - 1% (``xslope.fem._AT_YIELD_TOL``), where
tau = (s1 - s3)/2 and tau_f = c cos(phi) - s sin(phi) on the reduced envelope,
or its major principal stress has reached the tension cap (to the model's
round-off floor). The code's Mohr-Coulomb yield function is exactly f = tau -
tau_f, so the test is f >= -0.01 tau_f. Reporting only: no verdict reads the flag,
so no factor of safety moves.

Legs:

1. The point rule against an independent stress level, computed from the
   principal stresses by an eigen-decomposition of the 3x3 tensor rather than
   from the invariants the code uses: random states, phi = 0, 20 and 35 degrees,
   and exact states built at SL = 1.02, 1.0, 0.995 (flagged) and 0.985, 0.9
   (not). A point past the Mohr-Coulomb apex is flagged; a point with c = inf
   (a material held elastic on the Newton path) never is. The cap: s1 = T, and
   s1 within the floor below T, flagged; s1 two floors below T not.
2. A level layer of phi = 0 clay under its own weight, sides on rollers, where
   the answer is known in closed form. Elastic plane strain gives
   tau = gamma z (1 - K)/2 with K = nu/(1 - nu), so the clay reaches c at the
   depth z* = 2c / (gamma (1 - K)), and below it confined flow holds tau = c.
   Every element wholly deeper than z* + 0.5 m must be flagged and every element
   wholly shallower than z* - 0.5 m must not, on the default driver, on the
   plain viscoplastic one, and with the Newton corrector called at sweep 10 so it
   is the driver of record; ``plastic_fraction`` is the flag's mean; the stress
   plot outlines exactly the flagged elements. Before the change the deep test
   failed on all three (4, 4 and 6 of 117 deep elements unflagged), and the stress
   plot outlined by ``yield_function > 0``.
3. The same layer with every element held linear elastic (``elastic_mask``):
   its deep elements lie outside the envelope (yield_function > 0) and none is
   flagged.

Run directly:  PYTHONPATH=. python3 test/plastic_flag_check.py
"""

import contextlib
import copy
import io
import math
import warnings
from pathlib import Path

import numpy as np
from shapely.geometry import LineString, Polygon

warnings.filterwarnings("ignore")

BASE_FILE = (Path(__file__).resolve().parents[1] / "docs" / "fem" / "files"
             / "xslope_griffiths1.xlsx")
W, H = 12.0, 8.0          # level layer, m
GAMMA, C, NU = 20.0, 20.0, 0.3
TARGET = 1.0              # tri6 target size, m
MARGIN = 0.5              # m either side of z* left out of the depth comparison


def _quiet(fn, *a, **k):
    with contextlib.redirect_stdout(io.StringIO()):
        return fn(*a, **k)


# --------------------------------------------------------------------- leg 1
def _sl_independent(sig, c, phi):
    """Stress level from the principal stresses of the full 3x3 tensor."""
    sx, sy, txy, sz = sig
    ev = np.linalg.eigvalsh(np.array([[sx, txy, 0.0], [txy, sy, 0.0], [0.0, 0.0, sz]]))
    s1, s3 = ev[-1], ev[0]
    tau, s = 0.5 * (s1 - s3), 0.5 * (s1 + s3)
    tau_f = c * math.cos(phi) - s * math.sin(phi)
    return tau, tau_f, s1


def _state_at(sl, c, phi, s_ctr, rot):
    """A plane-strain stress (tension-positive) with in-plane centre s_ctr and stress
    level sl, sigma_z set to the centre so the in-plane circle is the largest."""
    tau_f = c * math.cos(phi) - s_ctr * math.sin(phi)
    r = sl * tau_f
    sx = s_ctr + r * math.cos(2 * rot)
    sy = s_ctr - r * math.cos(2 * rot)
    txy = r * math.sin(2 * rot)
    return np.array([sx, sy, txy, s_ctr])


def _leg_points(failures):
    from xslope.fem import _at_yield_points, _AT_YIELD_TOL
    floor = 0.5
    rng = np.random.default_rng(544)
    n_bad = 0
    for phi_deg in (0.0, 20.0, 35.0):
        phi = math.radians(phi_deg)
        c = 30.0
        sig = rng.normal(size=(4000, 4)) * 60.0 - np.array([80.0, 80.0, 0.0, 80.0])
        got = _at_yield_points(sig, np.full(len(sig), c), np.full(len(sig), math.sin(phi)),
                               np.full(len(sig), math.cos(phi)), floor=floor)
        for k, sg in enumerate(sig):
            tau, tau_f, _ = _sl_independent(sg, c, phi)
            tf = max(tau_f, floor)
            margin = tau - (1.0 - _AT_YIELD_TOL) * tf
            if abs(margin) < 1e-7 * max(1.0, tf):
                continue                       # on the line: either answer is right
            if bool(got[k]) != (margin > 0):
                n_bad += 1
        # exact states
        for sl, want in ((1.02, True), (1.0, True), (0.995, True),
                         (0.985, False), (0.9, False)):
            for s_ctr in (-150.0, -20.0):
                st = _state_at(sl, c, phi, s_ctr, rot=0.37)
                got1 = bool(_at_yield_points(st[None, :], [c], [math.sin(phi)],
                                             [math.cos(phi)], floor=floor)[0])
                if got1 != want:
                    failures.append(f"point rule: SL = {sl} at phi = {phi_deg:g}, centre "
                                    f"{s_ctr:g} read {'at' if got1 else 'inside'} yield")
    if n_bad:
        failures.append(f"point rule: {n_bad} random states disagree with the stress "
                        "level from the principal stresses")
    # past the apex (all-round tension in a c-phi soil): f >= 0, flagged
    phi = math.radians(30.0)
    st = np.array([[80.0, 70.0, 0.0, 75.0]])          # c/tan(phi) = 52 < s
    if not _at_yield_points(st, [30.0], [math.sin(phi)], [math.cos(phi)], floor=floor)[0]:
        failures.append("point rule: a state past the Mohr-Coulomb apex is not flagged")
    # a material held elastic (c = inf) is never at yield, however far out
    st = np.array([[-10.0, -900.0, 50.0, -300.0]])
    if _at_yield_points(st, [np.inf], [0.0], [1.0], floor=floor)[0]:
        failures.append("point rule: a point with c = inf is flagged")
    # the tension cap, T = 0, the shear well inside (c large)
    for s1, want in ((0.0, True), (1e-12, True), (-1e-12, True), (-0.5 * floor, True),
                     (-2.0 * floor, False)):
        st = np.array([[s1, -100.0, 0.0, -50.0]])
        got1 = bool(_at_yield_points(st, [500.0], [0.0], [1.0], t_cap=[0.0],
                                     floor=floor)[0])
        if got1 != want:
            failures.append(f"cap rule: s1 = {s1:g} with T = 0 read "
                            f"{'at' if got1 else 'below'} the cap")
    print(f"  point rule: random states {'agree' if not n_bad else 'DISAGREE'}; exact "
          "states, apex, elastic and cap cases checked")


# --------------------------------------------------------------------- legs 2-3
def _layer():
    from xslope.fem import build_fem_data
    from xslope.fileio import load_slope_data
    from xslope.mesh import build_mesh_from_polygons, get_material_polygons
    d = copy.deepcopy(_quiet(load_slope_data, str(BASE_FILE)))
    d["unit_system"] = "metric"
    poly = Polygon([(0.0, H), (W, H), (W, 0.0), (0.0, 0.0)])
    d["profile_lines"] = []
    d["polygons"] = [{"polygon": poly, "mat_id": 0}]
    d["domain_polygon"] = poly
    d["ground_surface"] = LineString([(0.0, H), (W, H)])
    d["max_depth"] = 0.0
    d["circles"], d["non_circ"], d["piezo_line"] = [], [], []
    d["piezo_phreatic"] = False
    d["materials"] = d["materials"][:1]
    d["materials"][0].update(dict(name="clay", gamma=GAMMA, gamma_sat=None, option="mc",
                                  c=C, phi=0.0, E=20000.0, nu=NU, t_cut=None,
                                  u="none", ru=0.0))
    mesh = _quiet(build_mesh_from_polygons, get_material_polygons(d),
                  target_size=TARGET, element_type="tri6")
    return _quiet(build_fem_data, d, mesh)


def _leg_layer(failures):
    from xslope.fem import solve_fem
    from xslope.plot_fem import _yielded_elements
    fd = _layer()
    nodes, els = np.asarray(fd["nodes"]), np.asarray(fd["elements"])
    ys = nodes[els[:, :3], 1]
    z_top, z_bot = H - ys.max(axis=1), H - ys.min(axis=1)
    K = NU / (1.0 - NU)
    z_star = 2.0 * C / (GAMMA * (1.0 - K))
    deep, shallow = z_top > z_star + MARGIN, z_bot < z_star - MARGIN
    if deep.sum() < 20 or shallow.sum() < 20:
        failures.append(f"layer: too few elements either side of z* = {z_star:.2f} m "
                        f"({deep.sum()} deep, {shallow.sum()} shallow)")
        return
    for solver, rungs in (("auto", None), ("viscoplastic", None), ("auto", (10,))):
        sol = _quiet(solve_fem, fd, F=1.0, max_iterations=4000, fem_solver=solver,
                     _corrector_rungs=rungs)
        flag = np.asarray(sol["plastic_elements"], dtype=bool)
        drv = (sol.get("corrector") or {}).get("driver_of_record", "viscoplastic")
        if rungs and drv != "corrector":
            failures.append(f"layer: the corrector called at sweep {rungs[0]} is not the "
                            f"driver of record ({drv}), so the Newton path's flag is "
                            "not exercised")
        miss, extra = int(np.sum(deep & ~flag)), int(np.sum(shallow & flag))
        old_vp = int(np.sum(np.asarray(sol["yield_function"])[deep] > 1e-8))
        print(f"  layer, {solver:12s} ({drv}, {sol.get('verdict')}): z* = {z_star:.2f} m; "
              f"{int(deep.sum())} elements below z*+{MARGIN:g}, {miss} unflagged; "
              f"{int(shallow.sum())} above z*-{MARGIN:g}, {extra} flagged; "
              f"yield_function > 1e-8 on {old_vp} of the deep ones")
        if not sol.get("converged"):
            failures.append(f"layer ({solver}): the solve did not converge")
        if miss:
            failures.append(f"layer ({solver}): {miss} of {int(deep.sum())} elements below "
                            f"z* = {z_star:.2f} m (+{MARGIN:g}) are at the strength of the "
                            "clay and not flagged")
        if extra:
            failures.append(f"layer ({solver}): {extra} elements above z* = {z_star:.2f} m "
                            f"(-{MARGIN:g}) are flagged; they carry at most "
                            f"{(z_star - MARGIN) / z_star:.2f} of their strength")
        if abs(float(sol["plastic_fraction"]) - float(flag.mean())) > 1e-12:
            failures.append(f"layer ({solver}): plastic_fraction "
                            f"{sol['plastic_fraction']:.4f} is not the flag's share "
                            f"{flag.mean():.4f}")
        if not np.array_equal(_yielded_elements(sol, len(els)), flag):
            failures.append(f"layer ({solver}): the stress plot does not outline the "
                            "flagged elements")

    # every element held linear elastic: outside the envelope, never flagged
    sol = _quiet(solve_fem, fd, F=1.0, max_iterations=4000,
                 elastic_mask=np.ones(len(els), dtype=bool))
    flag = np.asarray(sol["plastic_elements"], dtype=bool)
    n_out = int(np.sum(np.asarray(sol["yield_function"])[deep] > 0))
    print(f"  layer, all elastic: {n_out} of {int(deep.sum())} deep elements outside "
          f"the envelope, {int(flag.sum())} flagged")
    if n_out < deep.sum() // 2:
        failures.append(f"elastic layer: only {n_out} deep elements lie outside the "
                        "envelope, so the leg does not test the mask")
    if flag.any():
        failures.append(f"elastic layer: {int(flag.sum())} elements held linear elastic "
                        "are flagged as yielded")


def run():
    """Returns a list of failure strings (empty = pass)."""
    failures = []
    print("Plastic flag check:")
    _leg_points(failures)
    _leg_layer(failures)
    return failures


def main():
    import xslope
    print(f"xslope from {xslope.__file__}")
    failures = run()
    if failures:
        print("\nFAILURES:")
        for f in failures:
            print(f"  - {f}")
        raise SystemExit(1)
    print("\nA yielded element is one with a Gauss point at its strength, on both "
          "drivers, and never one held elastic.")


if __name__ == "__main__":
    main()
