"""Flow trace: particles that follow a steady seepage velocity field like meteors.

Each particle is lit for a fixed span of display time along its path: it
brightens, then dims and goes out, trailing a short streak of its last few
positions. Streaks are long where the flow is fast and short where it is slow.

The work splits in two. :func:`build_flow_trace` does everything once: it reads
the solved nodal velocity, finds the inflow boundaries, places release points
along them by equal flux, drops area particles, traces every path, and fixes
each particle's phase and lit window. :meth:`FlowTrace.frame` then turns a
display time into line segments and colors with nothing but interpolation along
the stored paths. Every particle's state is a function of ``(t - phase) mod P``,
so the frame at ``t = P`` is the frame at ``t = 0``: the loop has no seam.

Rendering (figure, background, animation, export) is
:func:`xslope.plot_seep.animate_flow_trace`.
"""

import logging

import numpy as np
import matplotlib.tri as mtri

logger = logging.getLogger(__name__)

SPEED_SCALES = ("true", "sqrt", "log")

# Decades of speed the "log" scale spreads over before it floors: a particle
# 10**-_LOG_DECADES slower than the reference speed moves at the floor speed.
_LOG_DECADES = 3.0
_LOG_FLOOR = 0.05


# ---------------------------------------------------------------------------
# The field: linear sub-triangles, nodal velocity, boundary edges
# ---------------------------------------------------------------------------

def _linear_subtriangles(elements, element_types):
    """Split every element into linear triangles over its own nodes.

    tri6 makes the four corner/mid-side triangles (as plot_seep contours them);
    quad4 two triangles; quad8 four corner triangles plus the inner mid-side
    quad as two; quad9 eight triangles around its center node. Every node that
    carries a solved value is a vertex, so the velocity is interpolated from all
    of them, and the sub-triangles tile the element exactly.
    """
    tris, owner = [], []
    for e, (en, et) in enumerate(zip(elements, element_types)):
        en = list(en)
        if et == 3:
            t = [en[:3]]
        elif et == 6:
            t = [[en[0], en[3], en[5]], [en[3], en[1], en[4]],
                 [en[5], en[4], en[2]], [en[3], en[4], en[5]]]
        elif et == 4:
            t = [[en[0], en[1], en[2]], [en[0], en[2], en[3]]]
        elif et == 8:
            t = [[en[0], en[4], en[7]], [en[4], en[1], en[5]],
                 [en[5], en[2], en[6]], [en[6], en[3], en[7]],
                 [en[4], en[5], en[6]], [en[4], en[6], en[7]]]
        elif et == 9:
            c = en[8]
            ring = [en[0], en[4], en[1], en[5], en[2], en[6], en[3], en[7]]
            t = [[ring[i], ring[(i + 1) % 8], c] for i in range(8)]
        else:
            raise ValueError(f"element type {et} is not supported")
        tris.extend(t)
        owner.extend([e] * len(t))
    return np.asarray(tris, dtype=int), np.asarray(owner, dtype=int)


def _nodal_porosity(porosity, n_nodes, elements, element_types, element_materials):
    """Porosity at each node: a scalar, or one value per material (1-based mats)."""
    if np.isscalar(porosity):
        return np.full(n_nodes, float(porosity))
    by_mat = np.asarray(porosity, dtype=float)
    tot = np.zeros(n_nodes)
    cnt = np.zeros(n_nodes)
    for en, et, m in zip(elements, element_types, element_materials):
        en = np.asarray(en[:et])
        tot[en] += by_mat[int(m) - 1]
        cnt[en] += 1
    with np.errstate(invalid="ignore"):
        out = tot / cnt
    out[cnt == 0] = np.nan
    return out


class TraceField:
    """A seepage solution's velocity on a linear triangulation, ready to trace."""

    def __init__(self, seep_data, solution, porosity=None):
        nodes = np.asarray(seep_data["nodes"], dtype=float)
        elements = seep_data["elements"]
        et = seep_data.get("element_types")
        if et is None:
            et = np.full(len(elements), 3)
        et = np.asarray(et)
        mats = seep_data.get("element_materials")
        if mats is None:
            mats = np.ones(len(elements), dtype=int)

        tris, owner = _linear_subtriangles(elements, et)
        # Drop zero-area slivers (they break the barycentric inverse).
        p = nodes[tris]
        area2 = ((p[:, 1, 0] - p[:, 0, 0]) * (p[:, 2, 1] - p[:, 0, 1])
                 - (p[:, 2, 0] - p[:, 0, 0]) * (p[:, 1, 1] - p[:, 0, 1]))
        # Orient every triangle counterclockwise so the outward normals below agree.
        flip = area2 < 0
        tris[flip] = tris[flip][:, [0, 2, 1]]
        area2 = np.abs(area2)
        keep = area2 > 1e-14 * max(np.ptp(nodes[:, 0]), np.ptp(nodes[:, 1])) ** 2
        tris, owner, area2 = tris[keep], owner[keep], area2[keep]

        self.nodes = nodes
        self.tris = tris
        self.tri_material = np.asarray(mats)[owner]
        self.area = 0.5 * area2
        self.h = np.sqrt(area2)          # a length per triangle (edge of the same area)
        self.triang = mtri.Triangulation(nodes[:, 0], nodes[:, 1], tris)
        self.finder = self.triang.get_trifinder()

        # Barycentric inverse: lambda_{1,2} = M @ (x - x0), lambda_0 = 1 - sum.
        p = nodes[tris]
        a = p[:, 1] - p[:, 0]
        b = p[:, 2] - p[:, 0]
        det = a[:, 0] * b[:, 1] - a[:, 1] * b[:, 0]
        self._x0 = p[:, 0]
        self._inv = np.stack([np.stack([b[:, 1], -b[:, 0]], -1),
                              np.stack([-a[:, 1], a[:, 0]], -1)], 1) / det[:, None, None]

        # Velocity. The solve's nodal velocity is the Darcy flux (element values
        # averaged to nodes). The particle speed is the seepage velocity q/n when a
        # porosity is given; otherwise the Darcy velocity, and the caller says so.
        darcy = np.asarray(solution["velocity"], dtype=float)
        self.darcy = darcy
        if porosity is not None:
            n = _nodal_porosity(porosity, len(nodes), elements, et, mats)
            if np.any(~np.isfinite(n[np.unique(tris)])) or np.any(n[np.unique(tris)] <= 0):
                raise ValueError("porosity must be positive on every node of the mesh")
            with np.errstate(invalid="ignore", divide="ignore"):
                self.vel = darcy / n[:, None]
            self.velocity_kind = "seepage"
        else:
            self.vel = darcy
            self.velocity_kind = "darcy"

        # Where water flows. An unconfined solve has a dry zone above the phreatic
        # surface whose velocity is the kr-scaled trickle; area particles go only
        # where the pressure head is not negative. A confined solve is saturated.
        head = np.asarray(solution["head"], dtype=float)
        pc = (head - nodes[:, 1])[tris].mean(axis=1)
        try:
            from .plot_seep import flownet_has_phreatic
            unconfined = bool(flownet_has_phreatic(seep_data, solution))
        except Exception:
            unconfined = bool(np.any(np.asarray(seep_data.get("bc_type", [])) == 2))
        self.unconfined = unconfined
        self.wet = pc >= 0 if unconfined else np.ones(len(tris), bool)

        self._boundary(seep_data)
        self.bnd_inflow = self._edge_inflow(seep_data, solution, elements, et)

        wet_nodes = np.unique(tris[self.wet]) if self.wet.any() else np.unique(tris)
        speed = np.linalg.norm(self.vel[wet_nodes], axis=1)
        speed = speed[np.isfinite(speed)]
        self.v_ref = float(np.percentile(speed, 90)) if speed.size else 1.0
        if not self.v_ref > 0:
            self.v_ref = float(speed.max()) if speed.size and speed.max() > 0 else 1.0
        x0, y0 = nodes[np.unique(tris)].min(axis=0)
        x1, y1 = nodes[np.unique(tris)].max(axis=0)
        self.bbox = (x0, x1, y0, y1)
        self.width = max(x1 - x0, y1 - y0)

    # -- boundary ------------------------------------------------------------
    def _boundary(self, seep_data):
        """Boundary edges of the triangulation, their outward normals, and which
        are open (a particle crossing it leaves the mesh) or no-flow walls."""
        tris = self.tris
        nb = self.triang.neighbors
        ti, k = np.nonzero(nb == -1)
        # matplotlib: neighbors[t, j] is across the edge from point j to point j+1.
        a = tris[ti, k]
        b = tris[ti, (k + 1) % 3]
        pa, pb = self.nodes[a], self.nodes[b]
        d = pb - pa
        L = np.hypot(d[:, 0], d[:, 1])
        # Counterclockwise triangles: the outward normal is the edge direction
        # turned clockwise.
        nrm = np.stack([d[:, 1], -d[:, 0]], -1) / L[:, None]
        self.bnd_a, self.bnd_b = a, b
        self.bnd_pa, self.bnd_pb, self.bnd_len, self.bnd_n = pa, pb, L, nrm
        self.bnd_tri = ti

        bc_type = seep_data.get("bc_type")
        flux = seep_data.get("flux_nodal")
        if bc_type is None:
            self.bnd_open = np.ones(len(a), bool)
            return
        bc_type = np.asarray(bc_type)
        is_open = bc_type != 0
        if flux is not None:
            is_open = is_open | (np.asarray(flux) != 0)
        self.bnd_open = is_open[a] & is_open[b]

    _ELEMENT_EDGES = {
        3: [(0, None, 1), (1, None, 2), (2, None, 0)],
        6: [(0, 3, 1), (1, 4, 2), (2, 5, 0)],
        4: [(0, None, 1), (1, None, 2), (2, None, 3), (3, None, 0)],
        8: [(0, 4, 1), (1, 5, 2), (2, 6, 3), (3, 7, 0)],
        9: [(0, 4, 1), (1, 5, 2), (2, 6, 3), (3, 7, 0)],
    }

    def _edge_inflow(self, seep_data, solution, elements, element_types):
        """Inflow density (flow per boundary length) at both ends of each
        boundary sub-edge, from the solve's nodal flows; None when the solution
        carries none (seeding then uses the nodal velocity)."""
        q = solution.get("q")
        bc = seep_data.get("bc_type")
        if q is None or bc is None or len(q) != len(self.nodes):
            return None
        bc = np.asarray(bc)
        inflow = np.where(bc != 0, np.asarray(q, dtype=float), 0.0)
        flux = seep_data.get("flux_nodal")
        if flux is not None:
            inflow = np.where(bc != 0, inflow, np.asarray(flux, dtype=float))
        is_open = (bc != 0) | (inflow != 0)
        # Boundary edges of the elements themselves (corner pair seen once).
        seen = {}
        for en, et in zip(elements, element_types):
            for a, m, b in self._ELEMENT_EDGES[int(et)]:
                key = (min(en[a], en[b]), max(en[a], en[b]))
                mid = en[m] if m is not None else None
                seen[key] = None if key in seen else (en[a], mid, en[b])
        edges = [v for v in seen.values() if v is not None]
        edges = [e for e in edges if all(is_open[x] for x in e if x is not None)]
        # Each node's flow over the boundary length its shape function covers,
        # integral(N_i ds): L/2 at each end of a linear edge; L/6 at the ends
        # and 2L/3 at the middle of a quadratic one. That is the row-sum lumped
        # boundary mass matrix, so the density is exact for a uniform flux on any
        # edge grading, and a quadratic edge's 1/6-2/3-1/6 loads come out flat.
        wN = np.zeros(len(self.nodes))
        for a, m, b in edges:
            if m is None:
                L = np.hypot(*(self.nodes[b] - self.nodes[a]))
                wN[a] += L / 2
                wN[b] += L / 2
            else:
                L = (np.hypot(*(self.nodes[m] - self.nodes[a]))
                     + np.hypot(*(self.nodes[b] - self.nodes[m])))
                wN[a] += L / 6
                wN[b] += L / 6
                wN[m] += 2 * L / 3
        rho = np.zeros(len(self.nodes))
        on = wN > 0
        rho[on] = inflow[on] / wN[on]
        self.inflow_total_q = float(np.sum(np.maximum(inflow[on], 0.0)))
        # Density at both ends of every boundary sub-edge (0 off the open edges).
        open_sub = set()
        for a, m, b in edges:
            if m is None:
                open_sub.add((min(a, b), max(a, b)))
            else:
                open_sub.add((min(a, m), max(a, m)))
                open_sub.add((min(m, b), max(m, b)))
        mask = np.array([(min(a, b), max(a, b)) in open_sub
                         for a, b in zip(self.bnd_a, self.bnd_b)])
        return np.column_stack([np.where(mask, rho[self.bnd_a], 0.0),
                                np.where(mask, rho[self.bnd_b], 0.0)])

    # -- evaluation ----------------------------------------------------------
    def locate(self, xy):
        return self.finder(xy[:, 0], xy[:, 1])

    def velocity_at(self, xy, ti=None):
        """Linear interpolation of the particle velocity; (nan, nan) outside."""
        xy = np.atleast_2d(xy)
        if ti is None:
            ti = self.locate(xy)
        out = np.full((len(xy), 2), np.nan)
        ok = ti >= 0
        if ok.any():
            t = ti[ok]
            l12 = np.einsum("nij,nj->ni", self._inv[t], xy[ok] - self._x0[t])
            lam = np.column_stack([1.0 - l12.sum(1), l12])
            out[ok] = np.einsum("ni,nij->nj", lam, self.vel[self.tris[t]])
        return out

    def display_speed_factor(self, speed, speed_scale):
        """g(|v| / v_ref): display speed relative to the reference particle's."""
        r = np.maximum(speed, 0.0) / self.v_ref
        if speed_scale == "true":
            return r
        if speed_scale == "sqrt":
            return np.sqrt(r)
        if speed_scale == "log":
            with np.errstate(divide="ignore"):
                g = 1.0 + np.log10(np.maximum(r, 1e-300)) / _LOG_DECADES
            return np.maximum(g, _LOG_FLOOR)
        raise ValueError(f"speed_scale must be one of {SPEED_SCALES}, got {speed_scale!r}")


# ---------------------------------------------------------------------------
# Seeding
# ---------------------------------------------------------------------------

def inflow_release_points(field, n, sub=8):
    """``n`` release points along the inflow boundaries, an equal share of the
    total inflow apart.

    The inflow density along the boundary comes from the solve's own nodal flows
    when the solution has them (``q``, the boundary reactions, whose positive sum
    is the reported flowrate, plus any specified-flux loads), each divided by the
    boundary length its node covers (see ``TraceField._edge_inflow``); without
    ``q`` it is the normal component of the nodal Darcy velocity. The density is
    linear along each boundary sub-edge, integrated over ``sub`` pieces with
    negative (outflow) parts clipped. Release points sit at the midpoints of
    ``n`` equal slices of the cumulative inflow, so each carries 1/n of it, like
    the flow channels of a flow net.

    Returns (points nudged just inside the mesh, the inflow each one carries, the
    total inflow found on the boundary).
    """
    if n <= 0:
        return np.zeros((0, 2)), 0.0, 0.0
    if field.bnd_inflow is not None:
        qa, qb = field.bnd_inflow[:, 0], field.bnd_inflow[:, 1]
    else:
        qa = -(field.darcy[field.bnd_a] * field.bnd_n).sum(1)   # inward flux density
        qb = -(field.darcy[field.bnd_b] * field.bnd_n).sum(1)
    cand = field.bnd_open & ((qa > 0) | (qb > 0))
    idx = np.nonzero(cand)[0]
    s_ = (np.arange(sub) + 0.5) / sub
    dens = np.maximum(qa[:, None] * (1 - s_) + qb[:, None] * s_, 0.0)
    piece_all = dens * (field.bnd_len[:, None] / sub)
    if idx.size == 0:
        return np.zeros((0, 2)), 0.0, 0.0
    # Walk the edges head to tail (counterclockwise triangles orient the whole
    # boundary one way round), so a slice that straddles two edges is one
    # contiguous stretch of boundary and the spacing stays equal across it.
    idx = _chain_order(field, idx)
    flat = piece_all[idx].ravel()
    cum = np.concatenate([[0.0], np.cumsum(flat)])
    total = cum[-1]
    if not total > 0:
        return np.zeros((0, 2)), 0.0, 0.0
    targets = (np.arange(n) + 0.5) / n * total
    j = np.clip(np.searchsorted(cum, targets, side="right") - 1, 0, flat.size - 1)
    frac = np.where(flat[j] > 0, (targets - cum[j]) / np.where(flat[j] > 0, flat[j], 1), 0.5)
    e = idx[j // sub]
    sp = ((j % sub) + np.clip(frac, 0, 1)) / sub
    pts = field.bnd_pa[e] + (field.bnd_pb[e] - field.bnd_pa[e]) * sp[:, None]
    # Nudge inside: along the inward normal, by a small part of the edge's triangle.
    inside = pts - field.bnd_n[e] * (1e-4 * field.h[field.bnd_tri[e]])[:, None]
    for _ in range(4):
        miss = field.locate(inside) < 0
        if not miss.any():
            break
        inside[miss] -= field.bnd_n[e[miss]] * (1e-2 * field.h[field.bnd_tri[e[miss]]])[:, None]
    return inside, total / n, total


def _chain_order(field, idx):
    """``idx`` (boundary edge indices) reordered into head-to-tail chains."""
    a, b = field.bnd_a[idx], field.bnd_b[idx]
    by_start = {int(n): k for k, n in enumerate(a)}
    ends = set(int(n) for n in b)
    used = np.zeros(len(idx), bool)
    order = []
    # Chain starts first (an edge no other candidate leads into), then loops.
    starts = [k for k in range(len(idx)) if int(a[k]) not in ends]
    for k0 in starts + list(range(len(idx))):
        k = k0
        while k is not None and not used[k]:
            used[k] = True
            order.append(k)
            k = by_start.get(int(b[k]))
    return idx[np.asarray(order, int)]


def area_points(field, n, rng):
    """``n`` points uniform by area over the wetted (flowing) part of the mesh."""
    w = field.area * field.wet
    if not w.sum() > 0:
        w = field.area.copy()
    t = rng.choice(len(w), size=n, p=w / w.sum())
    r1, r2 = rng.random(n), rng.random(n)
    fl = r1 + r2 > 1
    r1[fl], r2[fl] = 1 - r1[fl], 1 - r2[fl]
    p = field.nodes[field.tris[t]]
    return p[:, 0] + r1[:, None] * (p[:, 1] - p[:, 0]) + r2[:, None] * (p[:, 2] - p[:, 0])


# ---------------------------------------------------------------------------
# Tracing
# ---------------------------------------------------------------------------

def _first_crossing(field, p0, p1):
    """First boundary edge each segment p0->p1 crosses: (edge index or -1, t)."""
    n = len(p0)
    best_e = np.full(n, -1)
    best_t = np.full(n, np.inf)
    A, B = field.bnd_pa, field.bnd_pb
    # Only edges near the segment matter; chunk to bound memory.
    for s in range(0, n, 256):
        q0, q1 = p0[s:s + 256], p1[s:s + 256]
        r = q1 - q0                                   # (m, 2)
        e = B - A                                     # (E, 2)
        denom = r[:, None, 0] * e[None, :, 1] - r[:, None, 1] * e[None, :, 0]
        w = A[None] - q0[:, None]                     # (m, E, 2)
        with np.errstate(divide="ignore", invalid="ignore"):
            t = (w[..., 0] * e[None, :, 1] - w[..., 1] * e[None, :, 0]) / denom
            u = (w[..., 0] * r[:, None, 1] - w[..., 1] * r[:, None, 0]) / denom
        hit = (t >= -1e-9) & (t <= 1 + 1e-9) & (u >= -1e-9) & (u <= 1 + 1e-9) & (denom != 0)
        t = np.where(hit, t, np.inf)
        k = np.argmin(t, axis=1)
        tk = t[np.arange(len(k)), k]
        ok = np.isfinite(tk)
        best_e[s:s + 256][ok] = k[ok]
        best_t[s:s + 256][ok] = tk[ok]
    return best_e, best_t


def _nearest_edge(field, pts):
    """Index of the boundary edge nearest each point."""
    A, B = field.bnd_pa, field.bnd_pb
    d = B - A
    L2 = np.maximum((d ** 2).sum(1), 1e-300)
    out = np.empty(len(pts), int)
    for i, p in enumerate(pts):
        t = np.clip(((p - A) * d).sum(1) / L2, 0, 1)
        c = A + d * t[:, None]
        out[i] = np.argmin(((c - p) ** 2).sum(1))
    return out


def trace_paths(field, starts, speed_scale="sqrt", u_ref=1.0, max_tau=np.inf,
                ds_frac=0.3, max_steps=20000):
    """Trace every start point forward through the field, all at once.

    Integration is by arc length with a midpoint (RK2) step of ``ds_frac`` times
    the local triangle size along the interpolated velocity direction. Each step
    records the display time it took, ``ds / (u_ref * g(|v|/v_ref))``, so a path
    can later be sampled at any display time. A path ends when it crosses an open
    boundary edge (it left the mesh: an exit face, a drain, a head boundary), when
    its display time passes ``max_tau``, or when it stalls (no progress against a
    no-flow wall, or a zero velocity). On a no-flow wall the step's normal part is
    dropped and the particle slides along it.

    Returns a dict: ``xy`` (list of (m, 2) arrays), ``tau`` (list of (m,)),
    ``exited`` (bool array).
    """
    starts = np.asarray(starts, dtype=float)
    n = len(starts)
    pos = starts.copy()
    tau = np.zeros(n)
    alive = field.locate(pos) >= 0
    exited = np.zeros(n, bool)
    stall = np.zeros(n, int)
    rec_i, rec_xy, rec_tau = [np.arange(n)], [pos.copy()], [tau.copy()]
    for _ in range(max_steps):
        idx = np.nonzero(alive)[0]
        if idx.size == 0:
            break
        x = pos[idx]
        ti = field.locate(x)
        v0 = field.velocity_at(x, ti)
        sp0 = np.hypot(v0[:, 0], v0[:, 1])
        bad = ~(sp0 > 0) | ~np.isfinite(sp0) | (ti < 0)
        h = field.h[np.maximum(ti, 0)]
        ds = ds_frac * h
        with np.errstate(invalid="ignore", divide="ignore"):
            d0 = v0 / sp0[:, None]
        xm = x + 0.5 * ds[:, None] * d0
        tm = field.locate(xm)
        vm = field.velocity_at(xm, tm)
        spm = np.hypot(vm[:, 0], vm[:, 1])
        use0 = (tm < 0) | ~(spm > 0)
        vm[use0], spm[use0] = v0[use0], sp0[use0]
        with np.errstate(invalid="ignore", divide="ignore"):
            dm = vm / spm[:, None]
        xn = x + ds[:, None] * dm
        dtau = ds / (u_ref * field.display_speed_factor(spm, speed_scale))
        tn = field.locate(xn)

        out = (tn < 0) & ~bad
        step_len = ds.copy()
        if out.any():
            o = np.nonzero(out)[0]
            e, t = _first_crossing(field, x[o], xn[o])
            no_hit = e < 0
            # Left the mesh without crossing an edge we can find (a corner): treat
            # as a stall at the last good point.
            e_ok = ~no_hit
            leave = np.zeros(len(o), bool)
            leave[e_ok] = field.bnd_open[e[e_ok]]
            # Leaving through an open edge: end at the crossing point.
            lo = o[leave]
            frac = np.clip(t[leave], 0, 1)
            xn[lo] = x[lo] + (xn[lo] - x[lo]) * frac[:, None]
            dtau[lo] *= frac
            step_len[lo] *= frac
            exited[idx[lo]] = True
            # Hitting a wall: keep the tangential part, stay just inside.
            wo = o[e_ok & ~leave]
            if wo.size:
                ew = e[e_ok & ~leave]
                nrm = field.bnd_n[ew]
                d = xn[wo] - x[wo]
                d_t = d - (d * nrm).sum(1)[:, None] * nrm
                cand = x[wo] + d_t - nrm * (1e-6 * field.h[field.bnd_tri[ew]])[:, None]
                ok = field.locate(cand) >= 0
                # Sliding off the wall's end (a corner): out through an open edge
                # there means the particle leaves the mesh.
                if (~ok).any():
                    bad_c = np.nonzero(~ok)[0]
                    e2, _t2 = _first_crossing(field, x[wo[bad_c]], cand[bad_c])
                    near = np.where(e2 >= 0, e2, _nearest_edge(field, cand[bad_c]))
                    gone = np.zeros(len(wo), bool)
                    gone[bad_c] = field.bnd_open[near]
                    exited[idx[wo[gone]]] = True
                cand[~ok] = x[wo][~ok]
                moved = np.hypot(*(cand - x[wo]).T)
                frac = np.where(ds[wo] > 0, moved / ds[wo], 0)
                xn[wo] = cand
                dtau[wo] *= np.maximum(frac, 1e-3)
                step_len[wo] = moved
            # No crossing found: the step left through a corner. It leaves the
            # mesh if the nearest boundary edge is open, else it stalls in place.
            st = o[no_hit]
            if st.size:
                near = _nearest_edge(field, x[st])
                gone = field.bnd_open[near]
                exited[idx[st[gone]]] = True
                xn[st] = x[st]
                step_len[st] = np.where(gone, ds[st], 0.0)
                dtau[st[~gone]] = 0.0
                dtau[st[gone]] = 0.0
        progressed = step_len > 1e-3 * ds
        stall[idx] = np.where(progressed, 0, stall[idx] + 1)
        stall[idx[bad]] = 99
        good = ~bad
        g = idx[good]
        pos[g] = xn[good]
        tau[g] = tau[g] + np.where(np.isfinite(dtau[good]), dtau[good], np.inf)
        rec_i.append(g)
        rec_xy.append(pos[g].copy())
        rec_tau.append(tau[g].copy())
        done = exited[idx] | (stall[idx] >= 5) | ~np.isfinite(tau[idx]) | (tau[idx] > max_tau)
        alive[idx[done]] = False

    I = np.concatenate(rec_i)
    XY = np.concatenate(rec_xy)
    T = np.concatenate(rec_tau)
    order = np.argsort(I, kind="stable")
    I, XY, T = I[order], XY[order], T[order]
    cuts = np.searchsorted(I, np.arange(n + 1))
    xy_list, tau_list = [], []
    for k in range(n):
        xs, ts = XY[cuts[k]:cuts[k + 1]], T[cuts[k]:cuts[k + 1]]
        fin = np.isfinite(ts)
        xs, ts = xs[fin], ts[fin]
        # Keep display time strictly increasing (a stalled repeat adds nothing).
        if ts.size > 1:
            keep = np.concatenate([[True], np.diff(ts) > 0])
            xs, ts = xs[keep], ts[keep]
        xy_list.append(xs)
        tau_list.append(ts)
    return {"xy": xy_list, "tau": tau_list, "exited": exited}


# ---------------------------------------------------------------------------
# The trace: particles, schedule, frames
# ---------------------------------------------------------------------------

def _brightness(u, rise=0.25):
    """Meteor brightness over its lit span, u in [0, 1): up, then a slower fade."""
    u = np.asarray(u, dtype=float)
    up = np.clip(u / rise, 0, 1)
    down = np.clip((1 - u) / (1 - rise), 0, 1) ** 1.5
    return np.where(u < rise, up, down)


class FlowTrace:
    """Precomputed particles and their paths; :meth:`frame` draws any instant."""

    def __init__(self, field, paths, path_id, phase, s0, kind, loop_seconds,
                 lit_seconds, tail_seconds, n_tail, info):
        self.field = field
        self.P = float(loop_seconds)
        self.L = float(lit_seconds)
        self.tail = float(tail_seconds)
        self.n_tail = int(n_tail)
        self.path_id = np.asarray(path_id, int)
        self.phase = np.asarray(phase, float)
        self.s0 = np.asarray(s0, float)
        self.kind = np.asarray(kind)
        self.info = info
        # Concatenate the paths into one monotonic time axis: each path's display
        # time is offset so the ranges never overlap, and a vectorized searchsorted
        # finds any particle's position.
        taus, xys = paths["tau"], paths["xy"]
        dur = np.array([t[-1] if t.size else 0.0 for t in taus])
        self.path_end = dur
        offs = np.concatenate([[0.0], np.cumsum(dur + 1.0)[:-1]])
        self._off = offs
        self._T = np.concatenate([t + o for t, o in zip(taus, offs)])
        self._XY = np.concatenate(xys)
        cnt = np.array([t.size for t in taus])
        self._first = np.concatenate([[0], np.cumsum(cnt)[:-1]])
        self._last = self._first + cnt - 1

    def _position(self, pid, tq):
        """Positions of paths ``pid`` at display times ``tq`` (clamped to the path)."""
        tq = np.clip(tq, 0.0, self.path_end[pid])
        g = tq + self._off[pid]
        j = np.searchsorted(self._T, g, side="right")
        j = np.clip(j, self._first[pid] + 1, np.maximum(self._last[pid], self._first[pid] + 1))
        j = np.minimum(j, len(self._T) - 1)
        i = j - 1
        t0, t1 = self._T[i], self._T[j]
        with np.errstate(invalid="ignore", divide="ignore"):
            w = np.where(t1 > t0, (g - t0) / (t1 - t0), 0.0)
        w = np.clip(w, 0, 1)
        single = self._last[pid] == self._first[pid]
        w[single] = 0.0
        i = np.where(single, self._first[pid], i)
        j = np.where(single, self._first[pid], j)
        return self._XY[i] + (self._XY[j] - self._XY[i]) * w[:, None]

    def frame(self, t):
        """Segments (m, 2, 2), alpha (m,) and tail rank (m,) at display time ``t``.

        Rank 0 is the segment at the particle's head. Everything depends on
        ``(t - phase) mod P``, so ``frame(t + P)`` is ``frame(t)``.
        """
        age = np.mod(t - self.phase, self.P)
        lit = age < self.L
        if not lit.any():
            return np.zeros((0, 2, 2)), np.zeros(0), np.zeros(0, int)
        k = np.nonzero(lit)[0]
        pid = self.path_id[k]
        head = self.s0[k] + age[k]
        end = self.path_end[pid]
        nt = self.n_tail
        dt = self.tail / nt
        # Gone once the whole tail has drained into the path's end.
        show = head - nt * dt < end
        k, pid, head, end = k[show], pid[show], head[show], end[show]
        a = _brightness(age[k] / self.L)
        steps = np.arange(nt + 1)
        tq = head[:, None] - steps[None, :] * dt              # (m, nt+1)
        tq = np.minimum(tq, end[:, None])
        # An area particle has no path before its drop: its tail grows from there.
        tq = np.maximum(tq, 0.0)
        P = self._position(np.repeat(pid, nt + 1), tq.ravel()).reshape(len(k), nt + 1, 2)
        seg = np.stack([P[:, :-1], P[:, 1:]], axis=2).reshape(-1, 2, 2)
        fade = (1.0 - steps[:-1] / nt)                         # head 1 -> tail 1/nt
        alpha = (a[:, None] * fade[None, :]).ravel()
        rank = np.tile(steps[:-1], len(k))
        # Zero-length pieces (a tail clamped at its start or end) draw nothing.
        L2 = ((seg[:, 1] - seg[:, 0]) ** 2).sum(1)
        keep = (alpha > 0.003) & (L2 > 0)
        # Keep a lone dot visible for a particle that has not moved yet: a head
        # segment of zero length is drawn as a tiny cap by round line ends.
        head_only = (rank == 0) & (alpha > 0.003)
        keep |= head_only
        return seg[keep], alpha[keep], rank[keep]


def build_flow_trace(seep_data, solution, n_particles=1500, inflow_share=0.35,
                     speed_scale="sqrt", loop_seconds=10.0, lit_seconds=2.0,
                     tail_seconds=0.5, n_tail=8, streak_length=0.04,
                     porosity=None, seed=0, ds_frac=0.3, field=None):
    """Precompute a looping flow trace of a steady seepage solution.

    Parameters
    ----------
    n_particles : about how many meteors are lit at any moment.
    inflow_share : the part of them released from the inflow boundaries (0..1);
        the rest are area particles dropped at random over the wetted domain.
    speed_scale : "true" (display time proportional to travel time), "sqrt" or
        "log" (speeds compressed so slow zones still move; paths are unchanged).
    loop_seconds : P, the loop period in display seconds.
    lit_seconds : how long each meteor is lit (<= P).
    tail_seconds, n_tail : the streak is the last ``tail_seconds`` of the path,
        drawn as ``n_tail`` segments fading toward the tail.
    streak_length : the streak of a particle at the reference speed (the 90th
        percentile of the wetted nodes' speed), as a fraction of the model width.
        This fixes the time scale for every speed_scale.
    porosity : None (trace the Darcy velocity), a scalar, or one value per
        material; given, particles move at the seepage velocity q / n.
    """
    if speed_scale not in SPEED_SCALES:
        raise ValueError(f"speed_scale must be one of {SPEED_SCALES}, got {speed_scale!r}")
    if not 0 <= inflow_share <= 1:
        raise ValueError("inflow_share must be between 0 and 1")
    P = float(loop_seconds)
    # The lit span divides the loop: m spans per loop, each exactly P/m.
    m = max(1, int(round(P / min(float(lit_seconds), P))))
    L = P / m
    rng = np.random.default_rng(seed)
    if field is None:
        field = TraceField(seep_data, solution, porosity=porosity)
    u_ref = streak_length * field.width / tail_seconds    # model length / display s

    n_in = int(round(n_particles * inflow_share))
    n_area = int(n_particles) - n_in
    paths = {"xy": [], "tau": [], "exited": []}
    path_id, phase, s0, kind = [], [], [], []
    info = {"velocity_kind": field.velocity_kind, "speed_scale": speed_scale,
            "v_ref": field.v_ref, "u_ref": u_ref,
            # One display second stands for this much model time (true scale).
            "time_per_display_second": u_ref / field.v_ref}

    # Inflow: one release point per lit meteor wanted, each emitting m = P/L
    # particles per loop (own random phase and lit start each), so about one of
    # them is lit at any moment.
    if n_in:
        pts, q_each, q_total = inflow_release_points(field, n_in)
        info.update(n_release=len(pts), inflow_total=q_total,
                    flowrate=solution.get("flowrate"))
        if len(pts):
            # A trip need not be traced past the point where nothing can be lit:
            # cap the display time at many loops (true speed: deep slow trips).
            res = trace_paths(field, pts, speed_scale, u_ref, max_tau=200.0 * P,
                              ds_frac=ds_frac)
            per = m
            n_short = 0
            for r in range(len(pts)):
                xy, tt = res["xy"][r], res["tau"][r]
                if tt.size < 2:
                    n_short += 1
                    continue
                D = tt[-1]
                # A trip that ended inside the mesh (stalled or capped) keeps its
                # lit window, tail included, short of the end: nothing sits still.
                room = D - L if res["exited"][r] else D - L - tail_seconds
                if not res["exited"][r] and room <= 0:
                    n_short += 1
                    continue
                pidx = len(paths["tau"])
                paths["xy"].append(xy)
                paths["tau"].append(tt)
                paths["exited"].append(bool(res["exited"][r]))
                for _ in range(per):
                    path_id.append(pidx)
                    phase.append(rng.random() * P)
                    # The lit span starts at a random time along the trip.
                    s0.append(rng.random() * max(room, 0.0))
                    kind.append("inflow")
            info.update(inflow_exited=int(np.sum(res["exited"])),
                        inflow_dropped=n_short)

    # Area: each slot lives L_a = P/m, m lives per loop, each life from a fresh
    # random drop; the lives tile the loop, so the slot is always lit.
    if n_area:
        La = L
        drops = area_points(field, n_area * m, rng)
        res = trace_paths(field, drops, speed_scale, u_ref, max_tau=La + tail_seconds,
                          ds_frac=ds_frac)
        ph = rng.random(n_area) * P
        for j in range(n_area):
            for life in range(m):
                r = j * m + life
                xy, tt = res["xy"][r], res["tau"][r]
                if tt.size == 0:
                    continue
                pidx = len(paths["tau"])
                paths["xy"].append(xy)
                paths["tau"].append(tt)
                paths["exited"].append(bool(res["exited"][r]))
                path_id.append(pidx)
                phase.append(np.mod(ph[j] + life * La, P))
                s0.append(0.0)
                kind.append("area")
        info.update(area_slots=n_area, area_lives=m)
    info["lit_seconds"] = L
    if not paths["tau"]:
        raise ValueError("no particle could be traced: the field has no inflow and no wetted area")
    return FlowTrace(field, paths, path_id, phase, s0, kind, P, L, tail_seconds,
                     n_tail, info)
