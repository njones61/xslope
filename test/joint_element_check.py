"""The interface (joint) element, against four closed forms.

A jointed constraint line is a slip surface: the mesh is split along it and a
pair of zero-thickness interface elements spans each station, upper soil to bar
and bar to lower soil. Their law is Goodman's — elastic tractions
``t_n = k_n delta_n`` and ``t_s = k_s delta_t``, a Mohr-Coulomb limit
``|t_s| <= c_j + t_n tan phi_j``, a tension cutoff that opens the joint, and
perfectly plastic slip — integrated at the element's own nodes.

Four rows, each building its model in memory and solving it on the reference
(viscoplastic, NumPy) path:

  1. **Goodman direct shear.** A slab on a horizontal joint, driven by a
     horizontal seismic coefficient so the shear on the interface is known
     exactly. Below the limit the joint carries the whole applied shear and
     ``t_s = k_s delta_t`` at every node pair; at the limit every slipping pair
     sits EXACTLY on ``c_j + t_n tan phi_j``; and the seismic coefficient at
     which the slab stops standing brackets the closed form
     ``k = c_j / sigma_n + tan phi_j``. Run again at ``k_n`` x 100, the nodal
     tractions along one element must still agree with each other — the check on
     the integration rule, since Gauss quadrature on a zero-thickness interface
     oscillates at high normal stiffness and nodal (Lobatto) integration does
     not.
  2. **Block on an inclined plane.** The same slab on a plane at beta. It stands
     when ``tan phi_j > tan beta`` and does not when it is less, and the
     strength reduction returns ``FS = tan phi_j / tan beta``.
  3. **Pullout of a sheet.** A sheet buried in uniform soil, both ends free. The
     tension it can develop at a station is the interface shear integrated from
     the free end, and that has to be the SAME envelope
     ``fileio.reinforce_available_tension`` applies as a capacity —
     ``2 (a + sigma'_v tan delta) s`` per unit width — at every embedment.
  4. **The infinite-slope form.** A slab on a plane at beta with a COHESIVE
     interface, so the cohesion term the third row reads as an adhesion is read
     again as a strength: ``FS = (c_j + gamma H cos^2 beta tan phi_j) /
     (gamma H sin beta cos beta)``.

Three more legs read what the rows do not reach: the tension cutoff and the
tied end, on an imposed displacement field; the stiffness the SOFTER adjacent
soil sets on each side of a material crossing, and the line's own ``kn`` / ``ks``
overriding it; and that a model with no jointed line builds exactly the
``fem_data`` it always did. A last leg sweeps the joint stiffness over two
orders of magnitude on rows 2 and 4, on the strict failure criterion (see
:func:`_leg_stiffness` for why that is the one that answers the question).

**The end columns.** Rows 1, 2 and 4 need a slab that can move, and a model
truncated on rollers cannot give one: the slab's ends would be held. Each is
therefore built with a column of near-weightless, very soft ELASTIC material at
each end — empty space, in the way an excavation is modelled — so the slab has
two free faces. The jointed line runs the whole width and its ends are buried in
those columns, which is also what keeps the two soil faces' shared node at each
end of the line (the crack tip the mesher leaves there) out of the mechanism.

Run directly:  PYTHONPATH=. python3 test/joint_element_check.py
"""

import contextlib
import copy
import io
import math
import os
import sys
import time
import warnings

warnings.filterwarnings('ignore')

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

import numpy as np
from shapely.geometry import LineString, Polygon

from xslope.fem import build_fem_data, solve_fem, solve_ssrm
from xslope.fileio import (build_reinforce_lines, load_slope_data,
                           reinforce_available_tension)
from xslope.joint import (_joint_element_stiffness, joint_kinematics,
                          joint_state, mesh_has_joints, tie_internal_force)
from xslope.mesh import (add_intersection_points_to_polygons,
                         build_mesh_from_polygons,
                         extract_constraint_line_geometry,
                         get_material_polygons)

#: A shipped FEM model, used only for the boilerplate every slope_data carries
#: (units, gamma_water, solver options). Every row replaces the geometry, the
#: materials and the reinforcement wholesale.
BASE_FILE = os.path.join(_ROOT, 'docs', 'fem', 'files', 'xslope_griffiths1.xlsx')

GAMMA = 20.0            # kN/m3, the soil in every row
E_SOIL = 30000.0        # kPa
NU_SOIL = 0.3

_base_sd = None


def _base():
    global _base_sd
    if _base_sd is None:
        with contextlib.redirect_stdout(io.StringIO()):
            _base_sd = load_slope_data(BASE_FILE)
    return copy.deepcopy(_base_sd)


def _material(name, **kw):
    # The soil is deliberately far stronger than the interface in every row, so
    # the joint is the only thing that can fail and the closed form is a
    # statement about the joint alone. A soil that can yield beside the slab's
    # end faces puts its own mechanism into the answer.
    m = dict(name=name, gamma=GAMMA, gamma_sat=None, option='mc', c=100.0,
             phi=45.0, E=E_SOIL, nu=NU_SOIL, t_cut=None, u='none', ru=0.0)
    m.update(kw)
    return m


def _finish(d, polys_and_ids, domain, ground, y_bottom, materials, line,
            cj, phi_j, ts, s1d, kn=None, ks=None, jred='yes', tend1=0.0):
    """Common tail: install the geometry and the jointed line, mesh, build."""
    d['unit_system'] = 'metric'
    d['gamma_water'] = 9.81
    d['profile_lines'] = []
    d['polygons'] = [{'polygon': Polygon(r), 'mat_id': i} for r, i in polys_and_ids]
    d['domain_polygon'] = domain
    d['ground_surface'] = ground
    d['max_depth'] = y_bottom
    d['circles'] = []
    d['non_circ'] = []
    d['piezo_line'] = []
    d['piezo_phreatic'] = False
    d['materials'] = materials
    d['reinforcement_lines'] = [dict(
        label='joint', x1=line[0][0], y1=line[0][1], x2=line[1][0], y2=line[1][1],
        t_max=1.0e6, t_res=float('nan'), lp1=0.0, lp2=0.0, tend1=0.0, tend2=0.0,
        E=2.0e5, area=1.0e-5, spacing=None, adhesion=cj, delta=phi_j,
        kn=kn, ks=ks, jred=jred)]
    d['reinforcement_lines'][0]['tend1'] = tend1
    d['reinforce_lines'] = build_reinforce_lines(d['reinforcement_lines'])
    d['pile_lines'] = []
    lines, _n_r, _n_p = extract_constraint_line_geometry(d)
    polys = get_material_polygons(d, reinf_lines=lines)
    with contextlib.redirect_stdout(io.StringIO()):
        joint_lines = {0: {'tend1': tend1}} if tend1 else [0]
        mesh = build_mesh_from_polygons(polys, target_size=ts, element_type='tri6',
                                        lines=lines, element_size_1d=s1d,
                                        joint_lines=joint_lines)
        fem_data = build_fem_data(d, mesh)
    return d, mesh, fem_data


# --------------------------------------------------------------------------
# The slab model (rows 1, 2, 4): a soil block whose upper part rides on a
# jointed plane, with a soft, near-weightless elastic column at each end so the
# slab has two free faces and the line's crack tips sit clear of the mechanism.
# --------------------------------------------------------------------------

def slab_model(beta=20.0, HV=3.0, X1=36.0, XA=6.0, XB=30.0, YB=-6.0, VD=1.0,
               ts=1.5, s1d=1.0, phi_j=30.0, cj=0.0, E_void=0.1,
               k_seismic=0.0, kn=None, ks=None, jred='yes', k_scale=1.0):
    """A slab of vertical thickness ``HV`` on a jointed plane at ``beta``.

    The plane runs the full width; the slab is the soil between ``XA`` and
    ``XB``; the end columns run from ``VD`` below the plane to the ground
    surface.
    """
    T = math.tan(math.radians(beta))
    y = lambda x: x * T
    d = _base()
    soil = [(0.0, YB), (X1, YB), (X1, y(X1) - VD), (XB, y(XB) - VD),
            (XB, y(XB) + HV), (XA, y(XA) + HV), (XA, y(XA) - VD), (0.0, -VD)]
    vl = [(0.0, -VD), (XA, y(XA) - VD), (XA, y(XA) + HV), (0.0, HV)]
    vr = [(XB, y(XB) - VD), (X1, y(X1) - VD), (X1, y(X1) + HV), (XB, y(XB) + HV)]
    line = [(0.0, 0.0), (X1, y(X1))]
    rings = add_intersection_points_to_polygons([soil, vl, vr], [line])
    mats = [_material('soil'),
            _material('void', option='elastic', c=0.0, phi=0.0, E=E_void,
                      gamma=0.001, nu=0.2)]
    ids = [(rings[0], 0), (rings[1], 1), (rings[2], 1)]
    domain = Polygon([(0.0, HV), (X1, y(X1) + HV), (X1, YB), (0.0, YB)])
    ground = LineString([(0.0, HV), (X1, y(X1) + HV)])
    d, mesh, fem_data = _finish(d, ids, domain, ground, YB, mats, line,
                                cj, phi_j, ts, s1d, kn=kn, ks=ks, jred=jred)
    fem_data['k_seismic'] = k_seismic
    if k_scale != 1.0:
        jd = fem_data['joint_data']
        jd['kn'] = jd['kn'] * k_scale
        jd['ks'] = jd['ks'] * k_scale
        jd['K'] = _joint_element_stiffness(jd['w'], jd['tx'], jd['ty'],
                                           jd['nx'], jd['ny'], jd['kn'], jd['ks'])
    geom = dict(T=T, beta=beta, HV=HV, XA=XA, XB=XB, X1=X1, cj=cj, phi_j=phi_j,
                sigma_n=GAMMA * HV * math.cos(math.radians(beta)) ** 2,
                W=GAMMA * HV * (XB - XA))
    return d, fem_data, geom


def _slab_masks(fem_data, geom, margin=4.0):
    """Which joint elements lie well inside the slab, and which are the upper
    joint of their pair."""
    jd = fem_data['joint_data']
    x = fem_data['nodes'][jd['conn'][:, 0], 0]
    inside = (x > geom['XA'] + margin) & (x < geom['XB'] - margin)
    return inside, jd['side'] == 1


def _solve(fem_data, F=1.0, max_iterations=8000):
    with contextlib.redirect_stdout(io.StringIO()):
        return solve_fem(fem_data, F=F, max_iterations=max_iterations,
                         fast_kernel=False)


def _ssrm(fem_data, F_min=1.0, F_max=2.5, tolerance=0.005, **kw):
    # The iteration budget is capped well below the shipped ceiling. A penalty
    # interface makes the viscoplastic iteration crawl at a factor just under
    # the critical one -- the correction is a fixed force against a stiffness
    # the joint dominates -- and an uncapped trial there spends tens of
    # thousands of sweeps to decide what the next bisection step decides for it.
    kw.setdefault('max_iterations', 3000)
    kw.setdefault('max_iterations_ceiling', 6000)
    with contextlib.redirect_stdout(io.StringIO()):
        return solve_ssrm(fem_data, F_min=F_min, F_max=F_max,
                          tolerance=tolerance, **kw)


# --------------------------------------------------------------------------
# Row 1 — Goodman direct shear
# --------------------------------------------------------------------------

ROW1 = dict(beta=0.0, HV=3.0, X1=36.0, XA=6.0, XB=30.0, YB=-6.0, VD=1.0,
            ts=1.5, s1d=1.5, phi_j=20.0, cj=10.0, E_void=100.0)


def _leg_goodman(failures, results):
    """Direct shear: the elastic branch, the yield surface, and the slip load."""
    sigma = GAMMA * ROW1['HV']
    tanphi = math.tan(math.radians(ROW1['phi_j']))
    k_pred = ROW1['cj'] / sigma + tanphi
    d, fem_data, geom = slab_model(**ROW1)
    jd = fem_data['joint_data']
    inside, upper = _slab_masks(fem_data, geom)
    w = jd['w']
    live = inside[:, None] & (w > 0.0)

    # --- the elastic branch, at two loads well below the limit -------------
    line = []
    for k in (0.10, 0.20):
        fem_data['k_seismic'] = k
        sol = _solve(fem_data)
        if not sol['converged']:
            failures.append(f"row 1: the slab does not stand at k = {k:g}, well "
                            f"below the closed-form slip coefficient {k_pred:.4f}")
        dt, dn = joint_kinematics(jd, sol['displacements'])
        ts_res = abs(float((sol['joint_ts'][upper] * w[upper]).sum()))
        applied = k * geom['W']
        err = abs(ts_res - applied) / applied
        if err > 0.01:
            failures.append(f"row 1: at k = {k:g} the joint carries {ts_res:.2f} "
                            f"of the {applied:.2f} applied shear ({100*err:.2f}% "
                            f"adrift); the interface is not in equilibrium with "
                            f"the load")
        # t_s = k_s delta_t, node pair by node pair, on the state as solved.
        want = jd['ks'][:, None] * dt
        scale = max(float(np.max(np.abs(sol['joint_ts'][live]))), 1e-12)
        slope_err = float(np.max(np.abs(sol['joint_ts'][live] - want[live]))) / scale
        if slope_err > 1e-9:
            failures.append(f"row 1: at k = {k:g} the shear traction is not "
                            f"k_s times the relative tangential displacement "
                            f"(worst {slope_err:.2e} of the peak traction)")
        line.append((float(np.mean(np.abs(dt[live]))),
                     float(np.mean(np.abs(sol['joint_ts'][live])))))
        results.append(f"row 1  k={k:.2f}  joint shear {ts_res:.2f} kN/m vs "
                       f"applied {applied:.2f}  ({100*err:+.3f}%)")

    # The elastic slope, read as a secant between the two loads.
    (dt1, ts1), (dt2, ts2) = line
    k_s_measured = (ts2 - ts1) / (dt2 - dt1)
    k_s_stated = float(np.median(jd['ks'][inside]))
    if abs(k_s_measured - k_s_stated) / k_s_stated > 1e-6:
        failures.append(f"row 1: the elastic shear slope reads "
                        f"{k_s_measured:.1f} against the element's own k_s "
                        f"{k_s_stated:.1f}")
    results.append(f"row 1  elastic slope {k_s_measured:.1f} kPa/m vs k_s "
                   f"{k_s_stated:.1f} kPa/m")

    # --- the yield surface, pointwise, on a state with the joint at its limit
    fem_data['k_seismic'] = 0.5
    sol = _solve(fem_data)
    slipping = sol['joint_slipping'] | (
        np.abs(sol['joint_ts']) >= sol['joint_tlim'] - 1e-9)
    slipping &= (w > 0.0)
    n_slip = int(np.count_nonzero(slipping))
    if n_slip == 0:
        failures.append("row 1: no joint pair reaches its limit at k = 0.5, so "
                        "the yield surface is never exercised")
    else:
        want = ROW1['cj'] + np.maximum(sol['joint_tn'], 0.0) * tanphi
        rel = np.abs(np.abs(sol['joint_ts'][slipping]) - want[slipping]) / \
            np.maximum(want[slipping], 1e-12)
        worst = float(np.max(rel))
        if worst > 1e-6:
            failures.append(f"row 1: a slipping pair sits {worst:.2e} off "
                            f"c_j + t_n tan phi_j; the shear traction is not on "
                            f"the Mohr-Coulomb surface")
        results.append(f"row 1  {n_slip} pairs at the limit, worst departure "
                       f"from c_j + t_n tan phi_j {worst:.2e}")
    over = np.abs(sol['joint_ts']) - sol['joint_tlim']
    if float(np.max(over[w > 0.0])) > 1e-6:
        failures.append("row 1: a joint pair carries more shear than its "
                        "Mohr-Coulomb limit allows")

    # --- the slip load, read off the two branches -------------------------
    # The direct-shear curve is a straight elastic branch, t_s = k sigma_n,
    # meeting a plateau at the Mohr-Coulomb limit. Their intersection IS the
    # closed form, and reading it off the tractions rather than off whether the
    # slab still stands keeps the answer clear of the iteration budget: at a
    # coefficient just under the limit a penalty interface creeps, and how many
    # sweeps that creep is allowed decides a convergence bracket without saying
    # anything about the joint.
    fem_data['k_seismic'] = 0.90
    sol_hi = _solve(fem_data)
    plateau = float(np.mean(np.abs(sol_hi['joint_ts'][live])))
    k_slip = plateau / sigma
    err = abs(k_slip - k_pred) / k_pred
    if err > 0.02:
        failures.append(f"row 1: the shear traction plateaus at {plateau:.4f}, "
                        f"a slip coefficient of {k_slip:.4f} against the closed "
                        f"form c_j/sigma_n + tan phi_j = {k_pred:.4f} "
                        f"({100*err:+.2f}%)")
    results.append(f"row 1  plateau {plateau:.4f} kPa -> slip at k = "
                   f"{k_slip:.4f} vs closed form {k_pred:.4f} "
                   f"({100*(k_slip-k_pred)/k_pred:+.2f}%)")

    # --- the integration rule, at a hundred times the normal stiffness -----
    # Gauss quadrature on a zero-thickness interface oscillates MORE the stiffer
    # the joint is; nodal integration leaves the node pairs uncoupled, so the
    # traction profile is whatever the soil beside it delivers and does not
    # change with k_n. What is read is therefore the SPREAD of the three nodal
    # tractions within an element at the default stiffness and at a hundred
    # times it: the second must not exceed the first.
    def _spread(k_factor):
        _d, fd2, g2 = slab_model(**dict(ROW1, k_seismic=0.10))
        jd2 = fd2['joint_data']
        if k_factor != 1.0:
            jd2['kn'] = jd2['kn'] * k_factor
            jd2['K'] = _joint_element_stiffness(jd2['w'], jd2['tx'], jd2['ty'],
                                                jd2['nx'], jd2['ny'],
                                                jd2['kn'], jd2['ks'])
        sol2 = _solve(fd2)
        inside2, _u2 = _slab_masks(fd2, g2)
        tn2 = sol2['joint_tn']
        worst = 0.0
        for i in np.flatnonzero(inside2 & (jd2['n_pairs'] >= 3)):
            v = tn2[i, :3]
            worst = max(worst, float((v.max() - v.min())
                                     / max(abs(v.mean()), 1e-12)))
        return worst

    s1 = _spread(1.0)
    s100 = _spread(100.0)
    if s100 > s1 + 0.005:
        failures.append(f"row 1: the nodal tractions within one joint element "
                        f"spread {100*s1:.2f}% at the default normal stiffness "
                        f"and {100*s100:.2f}% at a hundred times it — the "
                        f"profile is oscillating with the penalty")
    if s100 > 0.05:
        failures.append(f"row 1: the nodal tractions within one joint element "
                        f"spread {100*s100:.2f}% on a uniform state")
    results.append(f"row 1  within-element traction spread {100*s1:.3f}% at "
                   f"k_n, {100*s100:.3f}% at k_n x 100")


# --------------------------------------------------------------------------
# Row 2 — block on an inclined plane
# --------------------------------------------------------------------------

ROW2 = dict(beta=20.0, HV=3.0, X1=36.0, XA=6.0, XB=30.0, YB=-6.0, VD=1.0,
            ts=1.5, s1d=1.0, phi_j=30.0, cj=0.0, E_void=1.0)


def _leg_incline(failures, results, k_scale=1.0, quiet=False, criterion=None):
    """Block on an inclined plane: it stands iff tan phi_j > tan beta, and the
    strength reduction returns tan phi_j / tan beta."""
    T = math.tan(math.radians(ROW2['beta']))
    expected = math.tan(math.radians(ROW2['phi_j'])) / T

    if k_scale == 1.0:
        # Stands with the joint stronger than the plane, and does not with it
        # weaker — the two sides of tan phi_j = tan beta.
        for phi_j, should_stand in ((ROW2['phi_j'], True),
                                    (0.6 * ROW2['beta'], False)):
            _d, fd, g = slab_model(**dict(ROW2, phi_j=phi_j))
            sol = _solve(fd)
            if bool(sol['converged']) != should_stand:
                failures.append(
                    f"row 2: with tan phi_j = {math.tan(math.radians(phi_j)):.4f} "
                    f"against tan beta = {T:.4f} the block "
                    f"{'must' if should_stand else 'must not'} stand, and it "
                    f"{'did not' if should_stand else 'did'}")
            results.append(f"row 2  phi_j = {phi_j:g} deg: "
                           f"{'stands' if sol['converged'] else 'slides'}")

    _d, fd, g = slab_model(**dict(ROW2, k_scale=k_scale))
    _kw = {} if criterion is None else {'failure_criterion': criterion}
    res = _ssrm(fd, F_min=1.0, F_max=2.5, tolerance=0.005, **_kw)
    FS = res.get('FS')
    if FS is None:
        failures.append("row 2: the strength reduction returned no factor of "
                        "safety")
        return None
    err = (FS - expected) / expected
    if k_scale == 1.0 and abs(err) > 0.03:
        failures.append(f"row 2: the strength reduction returns FS = {FS:.4f} "
                        f"against tan phi_j / tan beta = {expected:.4f} "
                        f"({100*err:+.2f}%)")
    if not quiet:
        results.append(f"row 2  SSRM FS = {FS:.4f} vs tan phi_j / tan beta = "
                       f"{expected:.4f}  ({100*err:+.2f}%)")
    return FS


# --------------------------------------------------------------------------
# Row 3 — pullout of a sheet
# --------------------------------------------------------------------------

ROW3 = dict(W=20.0, H=10.0, y_sheet=5.0, x1=5.0, x2=15.0, adhesion=5.0,
            delta=25.0, ts=1.0, s1d=0.25)


def _pullout_model(ts=None, s1d=None, tend1=0.0):
    p = ROW3
    d = _base()
    ring = [(0.0, p['H']), (p['W'], p['H']), (p['W'], 0.0), (0.0, 0.0)]
    domain = Polygon(ring)
    ground = LineString([(0.0, p['H']), (p['W'], p['H'])])
    line = [(p['x1'], p['y_sheet']), (p['x2'], p['y_sheet'])]
    rings = add_intersection_points_to_polygons([ring], [line])
    return _finish(d, [(rings[0], 0)], domain, ground, 0.0,
                   [_material('soil', c=1.0e5, phi=35.0)], line,
                   p['adhesion'], p['delta'], ts or p['ts'], s1d or p['s1d'],
                   tend1=tend1)


def _leg_pullout(failures, results):
    """The joints must produce the pullout envelope the capacity law applies."""
    p = ROW3
    d, mesh, fem_data = _pullout_model()
    jd = fem_data['joint_data']
    sol = _solve(fem_data)
    if not sol['converged']:
        failures.append("row 3: the buried sheet model did not reach "
                        "equilibrium under its own weight")
    upper = jd['side'] == 1
    w = jd['w']
    sigma_v = GAMMA * (p['H'] - p['y_sheet'])

    # The normal traction on a horizontal plane IS the vertical stress there.
    tn_mean = float((sol['joint_tn'][upper] * w[upper]).sum() / w[upper].sum())
    if abs(tn_mean - sigma_v) / sigma_v > 0.05:
        failures.append(f"row 3: the joint reads a mean normal traction of "
                        f"{tn_mean:.2f} where the overburden is {sigma_v:.2f}")
    results.append(f"row 3  mean normal traction {tn_mean:.2f} kPa vs "
                   f"overburden {sigma_v:.2f} kPa")

    # The capacity the two faces can develop between end 1 and a station, from
    # the interface limit the solve reports, against the envelope the capacity
    # law applies at the same station.
    x_a = fem_data['nodes'][jd['conn'][:, 0], 0]
    x_b = fem_data['nodes'][jd['conn'][:, 1], 0]
    x_mid = 0.5 * (x_a + x_b)
    length = p['x2'] - p['x1']
    line_row = d['reinforcement_lines'][0]
    # Embedments are taken from end 1 and no further than mid-span, where the
    # envelope's near-end branch is the one that governs and the comparison is
    # against the same integral the joints are developing.
    worst = 0.0
    for frac in (0.25, 0.375, 0.5):
        s = frac * length
        take = x_mid <= p['x1'] + s
        got = float((sol['joint_tlim'][take] * w[take]).sum())
        want = reinforce_available_tension(
            s, length - s, line_row['t_max'], 0.0, 0.0, 0.0, 0.0,
            pullout=line_row.get('_pullout_profile'))
        rel = abs(got - want) / want
        worst = max(worst, rel)
        if rel > 0.01:
            failures.append(
                f"row 3: at an embedment of {s:g} m the joints develop "
                f"{got:.2f} kN/m against the capacity envelope's {want:.2f} "
                f"kN/m ({100*rel:+.2f}%)")
        results.append(f"row 3  embedment {s:4.1f} m: joints {got:8.2f} kN/m, "
                       f"envelope {want:8.2f} kN/m ({100*(got-want)/want:+.2f}%)")

    # The bar on a jointed line has no bond-slip cap: its only limit is rupture.
    t_allow = np.asarray(fem_data['t_allow_by_1d_elem'])
    if t_allow.size and not np.allclose(t_allow, line_row['t_max']):
        failures.append("row 3: a bar on a jointed line still carries the "
                        "pullout envelope as its cap; on a jointed line the "
                        "grip is the joints' and the bar's only limit is Tmax")
    return worst


# --------------------------------------------------------------------------
# Row 4 — two-layer slope, the infinite-slope closed form
# --------------------------------------------------------------------------

ROW4 = dict(beta=18.0, HV=3.0, X1=36.0, XA=6.0, XB=30.0, YB=-6.0, VD=1.0,
            ts=1.5, s1d=1.0, phi_j=22.0, cj=6.0, E_void=1.0)


def _row4_expected():
    b = math.radians(ROW4['beta'])
    H = ROW4['HV']
    tanphi = math.tan(math.radians(ROW4['phi_j']))
    return ((ROW4['cj'] + GAMMA * H * math.cos(b) ** 2 * tanphi)
            / (GAMMA * H * math.sin(b) * math.cos(b)))


def _leg_infinite(failures, results, k_scale=1.0, quiet=False, criterion=None):
    expected = _row4_expected()
    _d, fd, g = slab_model(**dict(ROW4, k_scale=k_scale))
    _kw = {} if criterion is None else {'failure_criterion': criterion}
    res = _ssrm(fd, F_min=1.0, F_max=2.5, tolerance=0.005, **_kw)
    FS = res.get('FS')
    if FS is None:
        failures.append("row 4: the strength reduction returned no factor of "
                        "safety")
        return None
    err = (FS - expected) / expected
    if k_scale == 1.0 and abs(err) > 0.03:
        failures.append(f"row 4: the strength reduction returns FS = {FS:.4f} "
                        f"against the infinite-slope form {expected:.4f} "
                        f"({100*err:+.2f}%)")
    if not quiet:
        results.append(f"row 4  SSRM FS = {FS:.4f} vs infinite slope "
                       f"{expected:.4f}  ({100*err:+.2f}%)")
    return FS


# --------------------------------------------------------------------------
# The stiffness default, and the untouched path
# --------------------------------------------------------------------------

def _leg_stiffness(failures, results):
    """k_n and k_s at the default and an order of magnitude either side.

    Read on the STRICT failure criterion. The shipped hybrid one asks a
    non-converged trial for displacement evidence before calling it a failed
    slope, and the slip a viscoplastic sweep puts into a joint is the excess
    traction divided by k_s: a ten-times stiffer interface travels ten times
    less per sweep, so the same failing trial reaches the same residual having
    moved half as far, and the evidence test reads it as standing. The residual
    itself is not deceived — it is the same 2.8 at both stiffnesses on the same
    trial — so the criterion that reads the residual alone is the one that says
    whether the DEFAULT carries a result.
    """
    table = []
    for name, leg in (('row 2', _leg_incline), ('row 4', _leg_infinite)):
        row = {}
        for scale in (0.1, 1.0, 10.0):
            row[scale] = leg(failures, results, k_scale=scale, quiet=True,
                             criterion='non_convergence')
        table.append((name, row))
        vals = [v for v in row.values() if v is not None]
        if len(vals) == 3:
            spread = max(vals) - min(vals)
            if spread > 0.05:
                failures.append(
                    f"{name}: the factor of safety moves {spread:.4f} over two "
                    f"orders of magnitude of joint stiffness — the default is "
                    f"not a numerical penalty but a result")
            results.append(f"{name}  stiffness x0.1 / x1 / x10: "
                           + " / ".join(f"{row[s]:.4f}" for s in (0.1, 1.0, 10.0))
                           + f"  spread {spread:.4f}")
    return table


def _leg_untouched(failures, results):
    """A model with no jointed line builds exactly the fem_data it always did."""
    d, mesh, fem_data = _pullout_model()
    if not mesh_has_joints(mesh):
        failures.append("identity: the jointed model carries no joint elements")
    if 'joint_data' not in fem_data:
        failures.append("identity: the jointed model's fem_data carries no "
                        "joint_data")

    d2 = _base()
    p = ROW3
    ring = [(0.0, p['H']), (p['W'], p['H']), (p['W'], 0.0), (0.0, 0.0)]
    d2['unit_system'] = 'metric'
    d2['gamma_water'] = 9.81
    d2['profile_lines'] = []
    d2['polygons'] = [{'polygon': Polygon(ring), 'mat_id': 0}]
    d2['domain_polygon'] = Polygon(ring)
    d2['ground_surface'] = LineString([(0.0, p['H']), (p['W'], p['H'])])
    d2['max_depth'] = 0.0
    d2['circles'] = []
    d2['non_circ'] = []
    d2['piezo_line'] = []
    d2['piezo_phreatic'] = False
    d2['materials'] = [_material('soil', c=1.0e5, phi=35.0)]
    d2['reinforcement_lines'] = []
    d2['reinforce_lines'] = []
    d2['pile_lines'] = []
    with contextlib.redirect_stdout(io.StringIO()):
        mesh2 = build_mesh_from_polygons(get_material_polygons(d2),
                                         target_size=p['ts'],
                                         element_type='tri6')
        fem2 = build_fem_data(d2, mesh2)
    if mesh_has_joints(mesh2):
        failures.append("identity: an unjointed mesh reports joint elements")
    if 'joint_data' in fem2:
        failures.append("identity: an unjointed model's fem_data carries a "
                        "joint_data key, so it is not the dictionary it was "
                        "before joints existed")
    sol = _solve(fem2)
    for key in ('joint_tn', 'joint_ts', 'joint_slip', 'tie_forces'):
        if np.asarray(sol[key]).size:
            failures.append(f"identity: an unjointed solve reports a non-empty "
                            f"'{key}'")
    results.append("identity  an unjointed model writes no joint_data and "
                   "reports empty joint fields")


# --------------------------------------------------------------------------
# The tension cutoff and the tied end
# --------------------------------------------------------------------------

def _leg_opening_and_ties(failures, results):
    """A joint carrying tension opens; a tied end is plastic at its capacity.

    Neither is reached by the four closed-form rows — every one of them keeps
    the interface in compression and leaves both ends free — so both are read
    here on the real arrays, under a displacement field imposed on the buried
    sheet's joints and ties.
    """
    d, mesh, fem_data = _pullout_model(ts=2.0, s1d=1.0, tend1=25.0)
    jd = fem_data['joint_data']
    n_dof = int(fem_data['n_dof_total'])
    delta = 1.0e-3

    def _field(sign):
        """Move the UPPER SOIL face by ``sign * delta`` along the joint normal.

        Only that face: a bar node is side b of the upper joint and side a of
        the lower one, so moving every side-a node would move both faces of the
        upper joint together and leave it with no relative displacement at all.
        """
        u = np.zeros(n_dof)
        for i in np.flatnonzero(jd['side'] == 1):
            for pair in range(int(jd['n_pairs'][i])):
                dx, dy = jd['dof'][i, 2 * pair], jd['dof'][i, 2 * pair + 1]
                u[dx] = sign * delta * jd['nx'][i]
                u[dy] = sign * delta * jd['ny'][i]
        return u

    cj, tanphi = jd['cj'], jd['tanphi']
    st_open = joint_state(jd, _field(+1.0), cj, tanphi)
    # The upper joints are the ones the imposed field opens; a tip pair shares
    # its two nodes with the pair below it and never opens.
    live = (jd['w'] > 0.0) & (jd['side'] == 1)[:, None]
    interior = live & ~jd['tip']
    if not bool(st_open['open'][interior].all()):
        failures.append("opening: a joint pulled into tension past its cutoff "
                        "did not open")
    if float(np.max(np.abs(st_open['tn'][interior]))) > 0.0 or \
            float(np.max(np.abs(st_open['ts'][interior]))) > 0.0:
        failures.append("opening: an open joint still carries traction")

    st_shut = joint_state(jd, _field(-1.0), cj, tanphi)
    if bool(st_shut['open'][interior].any()):
        failures.append("opening: a joint pressed into contact is reported open")
    want = jd['kn'][:, None] * delta
    err = float(np.max(np.abs(st_shut['tn'][interior] - np.broadcast_to(
        want, st_shut['tn'].shape)[interior]))) / float(np.max(want))
    if err > 1e-12:
        failures.append(f"opening: a closed joint's normal traction is not "
                        f"k_n delta_n (worst {err:.2e})")
    results.append(f"opening  {int(interior.sum())} pairs open in tension and "
                   f"carry k_n delta_n in compression")

    td = jd.get('ties')
    if td is None or td['n'] == 0:
        failures.append("ties: a jointed line with a stated end anchorage "
                        "produced no tie")
        return
    line = d['reinforcement_lines'][0]
    le = float(np.median(fem_data['elem_length_1d'][
        np.asarray(fem_data['elem_length_1d']) > 0]))
    want_k = line['E'] * line['area'] / le
    if abs(float(td['k'][0]) - want_k) / want_k > 1e-6:
        failures.append(f"ties: the tie stiffness is {float(td['k'][0]):.4g}, "
                        f"not EA over one 1D element length ({want_k:.4g})")
    for factor, at_cap in ((0.5, False), (5.0, True)):
        u = np.zeros(n_dof)
        d_rel = factor * td['cap'][0] / td['k'][0]
        u[td['dof'][0, 0]] = d_rel
        f, _Ke, cap_flag = tie_internal_force(td, u)
        mag = float(np.hypot(f[0, 0], f[0, 1]))
        want = min(td['k'][0] * d_rel, td['cap'][0])
        if abs(mag - want) / want > 1e-9 or bool(cap_flag[0]) != at_cap:
            failures.append(f"ties: at {factor:g} times the capacity "
                            f"displacement the tie delivers {mag:.4f} rather "
                            f"than {want:.4f}")
    results.append(f"ties  {td['n']} tie(s), stiffness {float(td['k'][0]):.4g} "
                   f"kN/m, plastic at {float(td['cap'][0]):.4g} kN/m")


# --------------------------------------------------------------------------
# The derived stiffness across a material boundary
# --------------------------------------------------------------------------

def _leg_crossing_stiffness(failures, results):
    """The default k_n and k_s are the SOFTER adjacent soil's, side by side.

    A jointed line crossing a material boundary has a different soil under each
    half of it, so the derived stiffness has to change at the crossing. Built
    with a stiff soil on the left of x = 10 and a soft one on the right, and
    read straight off ``fem_data`` — no solve is needed to say what the arrays
    hold.
    """
    p = ROW3
    d = _base()
    x_cut = 0.5 * (p['x1'] + p['x2'])
    left = [(0.0, p['H']), (x_cut, p['H']), (x_cut, 0.0), (0.0, 0.0)]
    right = [(x_cut, p['H']), (p['W'], p['H']), (p['W'], 0.0), (x_cut, 0.0)]
    line = [(p['x1'], p['y_sheet']), (p['x2'], p['y_sheet'])]
    rings = add_intersection_points_to_polygons([left, right], [line])
    E_left, E_right = 60000.0, 15000.0
    mats = [_material('stiff', c=1.0e5, phi=35.0, E=E_left),
            _material('soft', c=1.0e5, phi=35.0, E=E_right)]
    _d, _mesh, fem_data = _finish(
        d, [(rings[0], 0), (rings[1], 1)], Polygon(left[:1] + right),
        LineString([(0.0, p['H']), (p['W'], p['H'])]), 0.0, mats, line,
        p['adhesion'], p['delta'], 1.0, 1.0)
    jd = fem_data['joint_data']
    x_mid = 0.5 * (fem_data['nodes'][jd['conn'][:, 0], 0]
                   + fem_data['nodes'][jd['conn'][:, 1], 0])
    nu = NU_SOIL
    # One element either side of the crossing stands on the crossing node and
    # reads BOTH soils, which is the rule working; the comparison is taken clear
    # of it.
    for tag, sel, E in (('left of the crossing', x_mid < x_cut - 1.5, E_left),
                        ('right of the crossing', x_mid > x_cut + 1.5, E_right)):
        if not np.any(sel):
            failures.append(f"crossing: no joint element {tag}")
            continue
        want_kn = E / (0.1 * jd['L'][sel])
        want_ks = (E / (2.0 * (1.0 + nu))) / (0.1 * jd['L'][sel])
        if not np.allclose(jd['kn'][sel], want_kn, rtol=1e-9):
            failures.append(f"crossing: the derived normal stiffness {tag} is "
                            f"not E / (0.1 x the element length) on that soil")
        if not np.allclose(jd['ks'][sel], want_ks, rtol=1e-9):
            failures.append(f"crossing: the derived shear stiffness {tag} is "
                            f"not G / (0.1 x the element length) on that soil")
        results.append(f"crossing  {tag}: k_n {float(jd['kn'][sel][0]):.0f}, "
                       f"k_s {float(jd['ks'][sel][0]):.0f} kPa/m (E = {E:g})")
    if float(np.max(jd['kn'][x_mid < x_cut - 1.5])) <= \
            float(np.max(jd['kn'][x_mid > x_cut + 1.5])):
        failures.append("crossing: the two sides of the material boundary carry "
                        "the same joint stiffness, so the softer-soil rule is "
                        "not being read per element")

    # A stated k_n / k_s on the line wins over the derived default.
    _d2, _m2, fd2 = _finish(
        _base(), [(rings[0], 0), (rings[1], 1)], Polygon(left[:1] + right),
        LineString([(0.0, p['H']), (p['W'], p['H'])]), 0.0, mats, line,
        p['adhesion'], p['delta'], 1.0, 1.0, kn=1234.0, ks=567.0)
    jd2 = fd2['joint_data']
    if not (np.allclose(jd2['kn'], 1234.0) and np.allclose(jd2['ks'], 567.0)):
        failures.append("crossing: a stated kn / ks on the line does not "
                        "override the derived default")
    results.append("crossing  a stated kn / ks on the line overrides the "
                   "derived default")


def run():
    """Returns a list of failure strings (empty = pass)."""
    failures, results = [], []
    t0 = time.time()
    _leg_goodman(failures, results)
    _leg_incline(failures, results)
    _leg_pullout(failures, results)
    _leg_infinite(failures, results)
    _leg_opening_and_ties(failures, results)
    _leg_crossing_stiffness(failures, results)
    _leg_untouched(failures, results)
    _leg_stiffness(failures, results)
    print("Interface (joint) element check "
          f"({time.time() - t0:.0f} s):")
    for line in results:
        print("  " + line)
    return failures


def main():
    failures = run()
    if failures:
        print("\nFAILURES:")
        for f in failures:
            print(f"  - {f}")
        raise SystemExit(1)
    print("\nThe interface element reproduces Goodman direct shear, the block "
          "on a plane, the pullout envelope and the infinite-slope form, and "
          "the joint stiffness default carries no result.")


if __name__ == '__main__':
    main()
