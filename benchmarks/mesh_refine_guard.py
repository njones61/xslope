"""Guard for feature-aware auto mesh refinement in ``build_mesh_from_polygons``.

The mesher accepts an optional ``refine_factor`` (with an optional ``refine_features``
list): when set, gmsh native size fields (Distance + Threshold per feature, composed
with a Min background field) drive the local element size down to
``target_size/refine_factor`` near model features — reinforcement/pile lines, crack /
notch tips, and thin material zones — growing smoothly back to ``target_size`` away
from them. With ``refine_factor=None`` (the default) no feature band is created and
the background size field is the requested size everywhere.

This guard asserts the contract on tiny, hand-checkable fixtures (no corpus files):

  1. INVARIANCE (byte-identical OFF). ``refine_factor=None`` produces a mesh whose
     node / element arrays are byte-identical to omitting the argument entirely, twice.

  2. REFINE NEAR A LINE. On a 10x10 square with one interior reinforcement line,
     ``refine_factor`` makes the mean element edge NEAR the line materially smaller
     than the target size, while elements FAR from the line stay near the target size
     (refinement is local, not global).

  3. DETERMINISM. The refined mesh is bit-identical when built twice from the same
     inputs (gmsh is seeded/ordered) — the regression locks depend on this.

  4. FEATURE DETECTION. ``detect_crack_tips`` finds the tip of a V-notch — a slit the
     material wraps around — and ignores ordinary corners, convex or not (the corpus
     cases are in test/quad_mesh_check.py); ``detect_thin_zones`` flags a slender band
     as a whole-thin zone sized to fit >= 3 elements across, and leaves a thick block
     alone.

  5. VALIDATION. ``refine_factor <= 1`` and an unknown ``refine_features`` entry both
     raise ValueError (they do not silently no-op).

  7. HIGH-CONTRAST INTERFACES (opt-in, seepage-specific). ``detect_interface_edges``
     flags a material boundary whose two sides differ in k1 by >= 100x and ignores a
     below-threshold jump; ``refine_features=['interfaces']`` + ``material_k`` refines
     locally on that boundary, deterministically. It is OFF the default set (a factor
     with ``refine_features`` omitted is byte-identical to OFF on a plain two-material
     block) and no-ops when ``material_k`` is not supplied.

  8. TERMINATION (the hang trap). Refinement on a section carrying embedded constraint
     lines must TERMINATE, at every factor the dialog can produce (1.05-8), not just
     produce a good mesh. Each build runs in a subprocess under a generous timeout,
     because the failure mode is a spin inside gmsh's C++ edge recovery that no
     in-process timer can interrupt.

Run from the repo root:  PYTHONPATH=. python3 benchmarks/mesh_refine_guard.py
Exits non-zero on any failure.
"""
import hashlib
import math
import os
import sys
import time

import numpy as np

from xslope.mesh import (build_mesh_from_polygons, detect_crack_tips,
                         detect_thin_zones, detect_interface_edges)

SQUARE = [(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0)]
REINF_LINE = [(2.0, 5.0), (8.0, 5.0)]     # interior, horizontal at y = 5

# Two stacked material zones sharing a horizontal boundary at y = 5. The lower zone
# (region 0) is ~100x less permeable than the upper (region 1) — a clay core / pervious
# shell in miniature. 'interfaces' refinement must concentrate on that shared edge only.
IFACE_LOW = [(0.0, 0.0), (10.0, 0.0), (10.0, 5.0), (0.0, 5.0)]     # region 0, low k
IFACE_HIGH = [(0.0, 5.0), (10.0, 5.0), (10.0, 10.0), (0.0, 10.0)]  # region 1, high k
IFACE_K = {0: 1.0e-7, 1: 1.0e-5}          # exactly 100x contrast across y = 5

# Wall-rooted, inclined, long constraint line — a VP60 soil nail in miniature. Rooted
# on the left vertical face at (0, 20), diving in at 15 deg for 16 units. The geo-kernel
# embed cannot recover such a line (its nodes float inside 2D triangles -> orphans ->
# singular stiffness); the OCC-fragment fallback must make it conform.
WALL_DOMAIN = [(0.0, 0.0), (20.0, 0.0), (20.0, 24.0), (0.0, 24.0)]
_c15, _s15 = math.cos(math.radians(15)), math.sin(math.radians(15))
WALL_LINE = [(0.0, 20.0), (16.0 * _c15, 20.0 - 16.0 * _s15)]

# THE HANG TRAP fixture (check 8). A thin inclined facing band against a bulk zone,
# with reinforcement rows running out of the band into the bulk — the reinforced-wall
# idiom, in miniature. Two things about it are what made the mesher spin:
#
#   * TRAP_LINES[0]'s first segment (0, 0)-(2, 0) lies ON the edge the two zones share,
#     so a constraint curve and a zone boundary want the same locus; and
#   * TRAP_LINES[1] starts inside the band and crosses into the bulk, so its segments
#     belong to different zones and no single surface contains the line.
#
# Either one, meshed with an embedded curve in a surface that does not contain it,
# puts gmsh into an unbounded edge-recovery retry ("Splitting those edges and trying
# again - level N"): a hang, not a slow mesh. It needs a refinement band to show up —
# without one the constraint curves carry too few nodes for the recovery to have work
# to fail at — and it is not monotonic in the factor, so no cap on the factor avoids it.
TRAP_THIN = [(0.0, 0.0), (10.0, 8.0), (12.0, 8.0), (2.0, 0.0)]
TRAP_BULK = [(2.0, 0.0), (12.0, 8.0), (30.0, 8.0), (30.0, -4.0),
             (-10.0, -4.0), (-10.0, 0.0), (0.0, 0.0)]
TRAP_LINES = [[(0.0, 0.0), (2.0, 0.0), (8.0, 0.0)],
              [(1.0, 0.4), (6.0, 0.4), (14.0, 0.4)]]
# Every factor the dialog can produce, spanning its 3.0 default. gmsh must finish on
# all of them; TRAP_TIMEOUT is generous because the assertion is termination, not speed
# (the builds land in well under a second — the pre-fix hang ran past seven minutes).
TRAP_FACTORS = (1.05, 1.2, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0)
TRAP_TIMEOUT = 120.0


def _orphan_1d(mesh):
    """Count 1D-element corner nodes not shared with any 2D element (should be 0 for a
    conforming reinforcement mesh)."""
    e1d = mesh.get("elements_1d")
    if e1d is None or len(e1d) == 0:
        return 0
    used = set(int(n) for n in np.unique(np.asarray(mesh["elements"])))
    return len({int(nd) for e in e1d for nd in e[:2] if int(nd) not in used})


def _trap_polygons():
    """The trap fixture's zones, prepared the way every caller prepares them: with the
    constraint-line/edge intersections inserted as polygon vertices. That integration
    is what makes 'the zone that contains this segment' well defined, so the guard
    exercises the same invariant the mesher relies on rather than a private one."""
    from xslope.mesh import add_intersection_points_to_polygons
    polys = [{"coords": list(TRAP_THIN), "mat_id": 0},
             {"coords": list(TRAP_BULK), "mat_id": 1}]
    return add_intersection_points_to_polygons(polys, TRAP_LINES)


def _trap_build(factor, element_type="tri6"):
    """One trap build. Runs in the CHILD process (see _trap_run)."""
    return build_mesh_from_polygons(_trap_polygons(), target_size=2.0,
                                    element_type=element_type,
                                    lines=TRAP_LINES, refine_factor=float(factor))


def _trap_run(factor, element_type="tri6"):
    """Build the trap fixture at ``factor`` in a SUBPROCESS under a timeout.

    A subprocess is the point: the failure this guards against is non-termination
    inside gmsh's C++ mesher, which no Python-level timer can interrupt. Run in-process
    it would hang the suite instead of failing it. Returns (status, detail) where
    status is 'ok', 'hang', or 'error'."""
    import subprocess
    cmd = [sys.executable, os.path.abspath(__file__), "--hang-trap",
           str(factor), element_type]
    env = dict(os.environ)
    env["PYTHONPATH"] = (os.getcwd() + os.pathsep + env.get("PYTHONPATH", "")).rstrip(os.pathsep)
    try:
        p = subprocess.run(cmd, capture_output=True, text=True,
                           timeout=TRAP_TIMEOUT, env=env)
    except subprocess.TimeoutExpired:
        return "hang", f"no mesh after {TRAP_TIMEOUT:g}s"
    for ln in reversed(p.stdout.splitlines()):
        if ln.startswith("TRAP_OK "):
            _tag, nodes, orph, secs = ln.split()
            return "ok", (int(nodes), int(orph), float(secs))
    tail = (p.stderr or p.stdout or "").strip().splitlines()
    return "error", (tail[-1][:160] if tail else "no output")


def _digest(mesh):
    h = hashlib.sha256()
    for k in ("nodes", "elements", "element_types", "element_materials"):
        h.update(np.ascontiguousarray(np.asarray(mesh[k])).tobytes())
    return h.hexdigest()[:16]


def _mean_edge(mesh):
    """Mean tri edge length and element centroids (linear tri3 corners)."""
    nodes, elems = mesh["nodes"], mesh["elements"]
    sizes, cents = [], []
    for e in elems:
        pts = nodes[e[:3]]
        sizes.append(np.mean([np.linalg.norm(pts[i] - pts[(i + 1) % 3]) for i in range(3)]))
        cents.append(pts.mean(axis=0))
    return np.array(sizes), np.array(cents)


def _dist_to_line(cents, a, b):
    a, b = np.array(a), np.array(b)
    ab = b - a
    out = []
    for p in cents:
        t = np.clip(np.dot(p - a, ab) / np.dot(ab, ab), 0.0, 1.0)
        out.append(np.linalg.norm(p - (a + t * ab)))
    return np.array(out)


def main():
    failures = []
    ts = 1.0
    poly = [{"coords": SQUARE, "mat_id": 0}]

    # (1) INVARIANCE: refine_factor=None is byte-identical to omitting the argument.
    m_omit = build_mesh_from_polygons(poly, target_size=ts, element_type="tri3",
                                      lines=[REINF_LINE])
    m_none = build_mesh_from_polygons(poly, target_size=ts, element_type="tri3",
                                      lines=[REINF_LINE], refine_factor=None)
    if _digest(m_omit) != _digest(m_none):
        failures.append("refine_factor=None differs from omitting the argument")
    else:
        # and it is itself reproducible
        m_omit2 = build_mesh_from_polygons(poly, target_size=ts, element_type="tri3",
                                           lines=[REINF_LINE])
        if _digest(m_omit) != _digest(m_omit2):
            failures.append("OFF mesh is not reproducible run-to-run")
        else:
            print(f"[invariance] refine_factor=None is byte-identical to OFF, twice "
                  f"({_digest(m_omit)}, {len(m_omit['nodes'])} nodes)")

    # (2) REFINE NEAR A LINE: local, not global.
    factor = 4.0
    m_ref = build_mesh_from_polygons(poly, target_size=ts, element_type="tri3",
                                     lines=[REINF_LINE], refine_factor=factor,
                                     refine_features=["reinforcement"])
    sizes, cents = _mean_edge(m_ref)
    d = _dist_to_line(cents, REINF_LINE[0], REINF_LINE[1])
    near = sizes[d < 0.5]
    far = sizes[d > 4.0]
    local = ts / factor
    if near.size == 0 or far.size == 0:
        failures.append("could not sample near/far element populations")
    elif near.min() >= 0.5 * ts:
        # the field drives the smallest near-line elements toward local = ts/factor
        failures.append(f"no genuinely small elements near the line: min {near.min():.3f} "
                        f"(target {ts}, local {local})")
    elif far.mean() <= 1.5 * near.mean():
        failures.append(f"far elements not coarse relative to near: far {far.mean():.3f} "
                        f"vs near {near.mean():.3f} (refinement not local)")
    else:
        print(f"[refine] near-line min edge {near.min():.3f} / mean {near.mean():.3f} "
              f"(local target {local}), far mean edge {far.mean():.3f} (target {ts}) — "
              f"refinement is local")

    # (3) DETERMINISM: refined mesh identical when rebuilt.
    m_ref2 = build_mesh_from_polygons(poly, target_size=ts, element_type="tri3",
                                      lines=[REINF_LINE], refine_factor=factor,
                                      refine_features=["reinforcement"])
    if _digest(m_ref) != _digest(m_ref2):
        failures.append(f"refined mesh not deterministic: {_digest(m_ref)} vs {_digest(m_ref2)}")
    else:
        print(f"[determinism] refined mesh reproducible ({_digest(m_ref)}, "
              f"{len(m_ref['nodes'])} nodes vs {len(m_omit['nodes'])} OFF)")

    # (4) FEATURE DETECTION (pure geometry, no gmsh).
    # A square with a narrow V-notch cut into the top; the tip is (5, 6).
    notch = [(0.0, 0.0), (10.0, 0.0), (10.0, 10.0),
             (5.05, 10.0), (5.0, 6.0), (4.95, 10.0), (0.0, 10.0)]
    tips = detect_crack_tips([notch])
    if tips != [(5.0, 6.0)]:
        failures.append(f"crack-tip detection wrong: {tips} (expected [(5.0, 6.0)])")
    elif detect_crack_tips([SQUARE]):
        failures.append("crack-tip detection false-positived on a plain square")
    else:
        print(f"[detect] V-notch tip found at {tips[0]}, plain corners ignored")

    # A slender band (20 x 0.4) is whole-thin; a thick block (10 x 10) is not.
    band = [(0.0, 0.0), (20.0, 0.0), (20.0, 0.4), (0.0, 0.4)]
    zones = detect_thin_zones([band, SQUARE], target_size=1.0)
    whole = [z for z in zones if z["kind"] == "whole" and z["poly_index"] == 0]
    if not whole:
        failures.append(f"thin band not flagged whole-thin: {zones}")
    elif any(z["poly_index"] == 1 for z in zones):
        failures.append("thick 10x10 block wrongly flagged thin")
    elif not (0.10 <= whole[0]["size"] <= 0.15):   # 0.4/3 ~ 0.133
        failures.append(f"thin-zone size {whole[0]['size']:.3f} not ~width/3")
    else:
        print(f"[detect] slender band -> whole-thin size {whole[0]['size']:.3f} "
              f"(>= 3 across 0.4 m), thick block ignored")

    # (5) VALIDATION: bad inputs raise, never silently no-op.
    try:
        build_mesh_from_polygons(poly, target_size=ts, element_type="tri3",
                                 lines=[REINF_LINE], refine_factor=1.0)
        failures.append("refine_factor=1.0 did not raise")
    except ValueError:
        pass
    try:
        build_mesh_from_polygons(poly, target_size=ts, element_type="tri3",
                                 lines=[REINF_LINE], refine_factor=3.0,
                                 refine_features=["bogus"])
        failures.append("unknown refine_features entry did not raise")
    except ValueError:
        print("[guard] refine_factor<=1 and unknown feature both rejected")

    # (6) WALL-ROOTED LINE CONFORMITY. A long inclined constraint line rooted on a
    # vertical face (a VP60 soil nail in miniature) must conform into the 2D mesh —
    # zero orphan 1D nodes — via the OCC-fragment fallback, deterministically. The
    # interior horizontal line (which the primary embed handles) must stay conforming
    # too, so the fallback logic never disturbs a line that already worked.
    wall_poly = [{"coords": WALL_DOMAIN, "mat_id": 0}]
    m_wall = build_mesh_from_polygons(wall_poly, target_size=1.5, element_type="tri6",
                                      lines=[WALL_LINE])
    m_wall2 = build_mesh_from_polygons(wall_poly, target_size=1.5, element_type="tri6",
                                       lines=[WALL_LINE])
    orph = _orphan_1d(m_wall)
    if orph:
        failures.append(f"wall-rooted line left {orph} orphan 1D nodes (did not conform)")
    elif _digest(m_wall) != _digest(m_wall2):
        failures.append(f"wall-rooted conforming mesh not deterministic: "
                        f"{_digest(m_wall)} vs {_digest(m_wall2)}")
    elif _orphan_1d(build_mesh_from_polygons(
            [{"coords": SQUARE, "mat_id": 0}], target_size=ts, element_type="tri6",
            lines=[REINF_LINE])) != 0:
        failures.append("interior reinforcement line regressed to orphan 1D nodes")
    else:
        print(f"[conform] wall-rooted inclined nail conforms (0 orphan 1D nodes, "
              f"{len(m_wall['nodes'])} nodes, deterministic {_digest(m_wall)}); "
              f"interior line still conforms")

    # (7) HIGH-CONTRAST MATERIAL INTERFACES (opt-in, seepage-specific). The pure
    # detector flags a >= 100x k-jump boundary and ignores a below-threshold one and
    # the outer edges. With 'interfaces' selected + material_k, refinement concentrates
    # on the shared boundary (near-interface elements materially smaller than far),
    # deterministically. It is NOT in the default set: a factor with refine_features
    # omitted must leave this two-material toy byte-identical to OFF (no line / crack /
    # thin feature here), and selecting 'interfaces' without material_k must no-op.
    iface_polys = [{"coords": IFACE_LOW, "mat_id": 0}, {"coords": IFACE_HIGH, "mat_id": 1}]
    region_ids = [0, 1]
    edges = detect_interface_edges([IFACE_LOW, IFACE_HIGH], region_ids, IFACE_K)
    edges_below = detect_interface_edges([IFACE_LOW, IFACE_HIGH], region_ids,
                                         {0: 1.0e-6, 1: 1.0e-5})   # only 10x
    if edges != [((0.0, 5.0), (10.0, 5.0))]:
        failures.append(f"interface detection wrong: {edges} "
                        f"(expected the single shared edge at y=5)")
    elif edges_below:
        failures.append(f"interface detector fired on a 10x (< 100x) contrast: {edges_below}")
    else:
        m_iface_off = build_mesh_from_polygons(iface_polys, target_size=ts,
                                               element_type="tri3")
        # default set (no 'interfaces') + a factor: byte-identical to OFF on this toy
        m_iface_def = build_mesh_from_polygons(iface_polys, target_size=ts,
                                               element_type="tri3", refine_factor=4.0)
        # 'interfaces' selected but no material_k: nothing to refine -> OFF
        m_iface_nok = build_mesh_from_polygons(iface_polys, target_size=ts,
                                               element_type="tri3", refine_factor=4.0,
                                               refine_features=["interfaces"])
        m_iface = build_mesh_from_polygons(iface_polys, target_size=ts, element_type="tri3",
                                           refine_factor=4.0, refine_features=["interfaces"],
                                           material_k=IFACE_K)
        m_iface2 = build_mesh_from_polygons(iface_polys, target_size=ts, element_type="tri3",
                                            refine_factor=4.0, refine_features=["interfaces"],
                                            material_k=IFACE_K)
        if _digest(m_iface_def) != _digest(m_iface_off):
            failures.append("'interfaces' leaked into the default feature set "
                            "(refine_features omitted changed the two-material mesh)")
        elif _digest(m_iface_nok) != _digest(m_iface_off):
            failures.append("'interfaces' without material_k did not no-op to OFF")
        elif _digest(m_iface) != _digest(m_iface2):
            failures.append(f"interface-refined mesh not deterministic: "
                            f"{_digest(m_iface)} vs {_digest(m_iface2)}")
        else:
            sizes, cents = _mean_edge(m_iface)
            d = np.abs(cents[:, 1] - 5.0)          # distance to the interface at y = 5
            near = sizes[d < 0.5]
            far = sizes[d > 4.0]
            local = ts / 4.0
            if near.size == 0 or far.size == 0:
                failures.append("interface: could not sample near/far element populations")
            elif near.min() >= 0.5 * ts:
                failures.append(f"interface: no genuinely small elements near the boundary: "
                                f"min {near.min():.3f} (target {ts}, local {local})")
            elif far.mean() <= 1.5 * near.mean():
                failures.append(f"interface: far elements not coarse relative to near: "
                                f"far {far.mean():.3f} vs near {near.mean():.3f} (not local)")
            else:
                print(f"[interface] >=100x k-jump edge detected (10x ignored); near-boundary "
                      f"min edge {near.min():.3f} / mean {near.mean():.3f} (local {local}), "
                      f"far mean {far.mean():.3f} (target {ts}) — local, deterministic, opt-in "
                      f"(off the default set, no-ops without material_k)")

    # (8) THE HANG TRAP: refinement on a section with embedded constraint lines must
    # TERMINATE, at every factor the dialog can produce. This is the one check whose
    # assertion is not about the mesh that comes back but about one coming back at all:
    # the regression it locks made "Refine near features" — a tick box whose factor
    # defaults to 3 — freeze Studio outright on any reinforced model. Because it is
    # not monotonic in the factor, the whole range is swept rather than sampled.
    trap_rows, trap_bad = [], []
    for _f in TRAP_FACTORS:
        status, detail = _trap_run(_f)
        trap_rows.append((_f, status, detail))
        if status != "ok":
            trap_bad.append(f"refine_factor={_f:g}: {detail}")
        elif detail[1]:
            trap_bad.append(f"refine_factor={_f:g}: {detail[1]} orphan 1D node(s)")
    if trap_bad:
        failures.append("constraint lines + refinement did not terminate cleanly — "
                        + "; ".join(trap_bad))
    else:
        _worst = max(r[2][2] for r in trap_rows)
        _sizes = [r[2][0] for r in trap_rows]
        print(f"[hang trap] {len(trap_rows)} factors "
              f"({TRAP_FACTORS[0]:g}-{TRAP_FACTORS[-1]:g}) on a thin band + bulk zone "
              f"with a boundary-riding and a zone-crossing reinforcement line: all "
              f"terminate (worst {_worst:.2f}s of {TRAP_TIMEOUT:g}s allowed), "
              f"{_sizes[0]}-{_sizes[-1]} nodes, no orphan 1D nodes")

    if failures:
        print("\nFAILED:")
        for f in failures:
            print("  -", f)
        return 1
    print("\nOK: refinement is byte-identical OFF, refines locally near features, is "
          "deterministic, detects V-notch tips and thin bands, and validates its inputs.")
    return 0


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "--hang-trap":
        # Child of _trap_run: build once, report, exit. Never reached in a normal run.
        _t0 = time.time()
        _m = _trap_build(sys.argv[2], sys.argv[3] if len(sys.argv) > 3 else "tri6")
        print(f"TRAP_OK {len(_m['nodes'])} {_orphan_1d(_m)} {time.time() - _t0:.2f}")
        sys.exit(0)
    sys.exit(main())
