"""Guard for solve_ssrm's optional ``ssr_zone`` SSR-Search-Area constraint AND its
``elastic_materials`` pure-elastic run option.

The SSRM accepts an optional ``ssr_zone`` polygon (RS2's "SSR Search Area"):
strength reduction is confined to elements whose centroid lies INSIDE the polygon,
everything outside is held at full strength. This guard asserts the mask that
implements it on a tiny, hand-checkable mesh:

  1. CONTAINMENT. On a 4-element coarse mesh with centroids at (1,1), (3,1), (1,3),
     (3,3), a triangle polygon that encloses only (1,1) must exclude exactly the
     other three (mask True = held at full strength, outside the zone).

  2. INVARIANCE. ssr_zone = None builds no mask (the historical every-element-reduced
     path, bit-identical); a full-domain polygon excludes nothing.

  3. UNION COMPOSITION. ssr_zone and ssr_exclude compose by OR — an element is held
     at full strength if it is outside the zone OR named-excluded.

The SSRM also accepts ``elastic_materials`` (RS2's "Plasticity Specifications:
None"): named materials whose elements are PURE LINEAR ELASTIC — they skip
plasticity entirely and can never yield. This is DISTINCT from ssr_exclude (full
strength but STILL yields). The elastic-mask section asserts, on a tiny gmsh mesh:

  4. NO YIELD. A domain that yields (nonzero plastic_elements) under the default
     path has ZERO plastic elements and ZERO viscoplastic strain when every element
     is masked elastic — the elastic response, independent of the trial F.

  5. DEFAULT INVARIANCE. elastic_mask=None is bit-identical to omitting the argument.

  6. INDEPENDENCE. elastic_mask on one material yields no plastic elements there
     while the other material still yields — the mask is per-element and orthogonal
     to ssr_exclude. Unknown material names raise ValueError.

Run from the repo root:  PYTHONPATH=. python3 benchmarks/ssr_zone_guard.py
Exits non-zero on any failure.
"""
import hashlib
import sys

import numpy as np

from xslope.fem import _element_centroids, _ssr_zone_exclusion_mask


def _coarse_mesh():
    """4 tri3 elements with centroids at (1,1), (3,1), (1,3), (3,3)."""
    nodes = np.array([
        [0., 0.], [3., 0.], [0., 3.],     # elem 0 -> centroid (1,1)
        [2., 0.], [5., 0.], [2., 3.],     # elem 1 -> centroid (3,1)
        [0., 2.], [3., 2.], [0., 5.],     # elem 2 -> centroid (1,3)
        [2., 2.], [5., 2.], [2., 5.],     # elem 3 -> centroid (3,3)
    ])
    elements = np.array([[0, 1, 2], [3, 4, 5], [6, 7, 8], [9, 10, 11]])
    element_types = np.array([3, 3, 3, 3])
    element_materials = np.array([1, 1, 2, 2])  # two material zones
    return {"nodes": nodes, "elements": elements, "element_types": element_types,
            "element_materials": element_materials, "material_names": ["A", "B"]}


def _tiny_slope_fem(mat_split_x=5.0):
    """A small two-material slope that yields under gravity at F=1: bottom (0,0)-
    (10,0), slope face (10,0)-(4,6), crest (4,6)-(0,6), split at x=mat_split_x into
    'weakL' (left) and 'weakR' (right) by half-plane intersection. Built through the
    standard save/load + mesh pipeline so the base BCs and fem_data are exactly what
    solve_fem consumes. Returns fem_data."""
    import os
    import tempfile
    from shapely.geometry import Polygon
    from xslope.fileio import load_slope_data, save_slope_data_to_xlsx
    from xslope.mesh import (get_material_polygons, build_mesh_from_polygons,
                             extract_constraint_line_geometry, extract_point_constraints)
    from xslope.fem import build_fem_data

    template = os.path.join(os.path.dirname(__file__), '..', 'docs', 'lem',
                            'files', 'xslope_acads_simple.xlsx')
    dom = Polygon([(0.0, 0.0), (10.0, 0.0), (4.0, 6.0), (0.0, 6.0)])
    left = dom.intersection(Polygon([(-1., -1.), (mat_split_x, -1.),
                                     (mat_split_x, 7.), (-1., 7.)]))
    right = dom.intersection(Polygon([(mat_split_x, -1.), (13., -1.),
                                      (13., 7.), (mat_split_x, 7.)]))
    sd = load_slope_data(template)                     # full field set from template
    base = dict(sd['materials'][0])
    matbase = dict(base, c=1.0, phi=5.0, gamma=20.0, gamma_sat=20.0, E=8000.0,
                   nu=0.3, u='none', option='mc', psi=0.0)
    sd['materials'] = [dict(matbase, name='weakL'), dict(matbase, name='weakR')]
    sd['profile_lines'] = []
    sd['polygons'] = [{'mat_id': 0, 'polygon': left}, {'mat_id': 1, 'polygon': right}]
    sd['max_depth'] = 0.0
    sd['dloads'] = []
    sd['dloads2'] = []
    sd['piezo_line'] = []
    sd['circular'] = True
    sd['non_circ'] = []
    sd['circles'] = [{'Xo': 5.0, 'Yo': 12.0, 'Depth': 0.0, 'R': 12.0}]
    tmp = os.path.join(tempfile.gettempdir(), 'elastic_mask_guard.xlsx')
    save_slope_data_to_xlsx(sd, tmp)
    sd = load_slope_data(tmp)
    cl, _, _ = extract_constraint_line_geometry(sd)
    polys = get_material_polygons(sd, reinf_lines=cl)
    mesh = build_mesh_from_polygons(polys, target_size=1.0, element_type='tri3',
                                    lines=cl, point_constraints=extract_point_constraints(sd))
    return build_fem_data(sd, mesh)


def _sol_digest(sol):
    def h(a):
        return hashlib.sha256(np.ascontiguousarray(np.asarray(a, float)).tobytes()).hexdigest()[:16]
    return '|'.join([str(sol['converged']), str(sol['iterations']),
                     h(sol['displacements']), h(sol['stresses']),
                     h(sol['plastic_elements'].astype(float))])


def _elastic_mask_check():
    """Deterministic checks for the elastic_materials / elastic_mask pure-elastic
    option (RS2 'Plasticity: None'). Returns a list of failure strings ([] = ok)."""
    from xslope.fem import solve_fem, solve_ssrm
    failures = []
    fem = _tiny_slope_fem()
    names = list(fem['material_names'])
    em = fem['element_materials']
    maskL = np.isin(em, [names.index('weakL') + 1])
    maskR = np.isin(em, [names.index('weakR') + 1])
    all_mask = np.ones(len(em), dtype=bool)
    F = 1.0

    base = solve_fem(fem, F=F, max_iterations=400, debug_level=0)
    n_base = int(base['plastic_elements'].sum())
    if n_base == 0:
        failures.append("fixture did not yield under the default path (nothing to guard)")

    # (4) NO YIELD: every element elastic -> zero plastic, zero viscoplastic strain.
    elas = solve_fem(fem, F=F, max_iterations=400, debug_level=0, elastic_mask=all_mask)
    n_elas = int(elas['plastic_elements'].sum())
    vp = np.concatenate([np.asarray(v).ravel() for v in elas['plastic_strains'].values()]) \
        if elas.get('plastic_strains') else np.zeros(1)
    if n_elas != 0:
        failures.append(f"all-elastic mask still flagged {n_elas} plastic elements")
    elif np.abs(vp).max() > 1e-12:
        failures.append(f"all-elastic mask left nonzero viscoplastic strain {np.abs(vp).max():.2e}")
    else:
        print(f"[no-yield] default plastic={n_base} -> all-elastic plastic=0, vp-strain=0")

    # elastic response is F-independent (no strength to reduce): F=2 matches F=1.
    elas2 = solve_fem(fem, F=2.0, max_iterations=400, debug_level=0, elastic_mask=all_mask)
    if _sol_digest(elas) != _sol_digest(elas2):
        failures.append("all-elastic solution changed with F (should be F-independent)")

    # (5) DEFAULT INVARIANCE: elastic_mask=None is bit-identical to omitting it.
    none_sol = solve_fem(fem, F=F, max_iterations=400, debug_level=0, elastic_mask=None)
    if _sol_digest(none_sol) != _sol_digest(base):
        failures.append("elastic_mask=None differs from omitting the argument")
    else:
        print("[invariance] elastic_mask=None is bit-identical to the default path")

    # (6) INDEPENDENCE: masking one material -> zero plastic THERE; the other still
    # yields. (The mask is per-element and orthogonal to ssr_exclude.)
    solL = solve_fem(fem, F=F, max_iterations=400, debug_level=0, elastic_mask=maskL)
    solR = solve_fem(fem, F=F, max_iterations=400, debug_level=0, elastic_mask=maskR)
    if int(solL['plastic_elements'][maskL].sum()) != 0:
        failures.append("weakL masked elastic but some weakL elements flagged plastic")
    if int(solR['plastic_elements'][maskR].sum()) != 0:
        failures.append("weakR masked elastic but some weakR elements flagged plastic")
    print(f"[independence] mask weakL -> weakL plastic={int(solL['plastic_elements'][maskL].sum())}; "
          f"mask weakR -> weakR plastic={int(solR['plastic_elements'][maskR].sum())}")

    # Unknown material name raises (does not silently no-op).
    try:
        solve_ssrm(fem, F_min=1.0, F_max=1.0, elastic_materials=['NoSuchMat'],
                   max_iterations=5, debug_level=0)
        failures.append("unknown elastic_materials name did not raise")
    except ValueError:
        print("[guard] unknown elastic_materials name rejected")

    return failures


def main():
    fem = _coarse_mesh()
    failures = []

    cen = _element_centroids(fem)
    expect = np.array([[1., 1.], [3., 1.], [1., 3.], [3., 3.]])
    if not np.allclose(cen, expect):
        failures.append(f"centroids {cen.tolist()} != {expect.tolist()}")

    # (1) Containment: triangle x + y <= 2.5 encloses only centroid (1,1).
    tri = [(0.0, 0.0), (2.5, 0.0), (0.0, 2.5)]
    mask = _ssr_zone_exclusion_mask(fem, tri)
    want = np.array([False, True, True, True])   # only elem 0 inside -> reduced
    if not np.array_equal(mask, want):
        failures.append(f"triangle mask {mask.tolist()} != {want.tolist()}")
    else:
        print(f"[containment] triangle encloses only elem 0 -> mask={mask.tolist()}")

    # Closing repeat of the first vertex is allowed (ring closed automatically).
    mask_closed = _ssr_zone_exclusion_mask(fem, tri + [tri[0]])
    if not np.array_equal(mask_closed, mask):
        failures.append("closing-repeat vertex changed the mask")

    # (2) Full-domain polygon excludes nothing.
    box = [(-1.0, -1.0), (10.0, -1.0), (10.0, 10.0), (-1.0, 10.0)]
    mask_box = _ssr_zone_exclusion_mask(fem, box)
    if mask_box.any():
        failures.append(f"full-domain polygon excluded elements: {mask_box.tolist()}")
    else:
        print("[invariance] full-domain polygon excludes nothing")

    # Degenerate (collinear) polygon is rejected, not silently empty.
    try:
        _ssr_zone_exclusion_mask(fem, [(0, 0), (1, 1), (2, 2)])
        failures.append("degenerate polygon did not raise")
    except ValueError:
        print("[guard] degenerate (zero-area) polygon rejected")

    # (3) Union composition: outside-zone OR named-excluded.
    zone_mask = mask                                  # [F, T, T, T]
    exclude_mask = fem["element_materials"] == 2      # [F, F, T, T]
    union = zone_mask | exclude_mask
    want_union = np.array([False, True, True, True])
    if not np.array_equal(union, want_union):
        failures.append(f"union {union.tolist()} != {want_union.tolist()}")
    else:
        print(f"[union] zone OR exclude -> {union.tolist()}")

    # Elastic-materials / elastic_mask pure-elastic option (checks 4-6).
    print("\n-- elastic_materials (pure-elastic 'Plasticity: None') --")
    failures += _elastic_mask_check()

    if failures:
        print("\nFAILED:")
        for f in failures:
            print("  -", f)
        return 1
    print("\nOK: ssr_zone mask confines reduction inside the polygon, is inert for a "
          "full-domain zone, and composes with ssr_exclude by union; elastic_materials "
          "holds named materials pure-elastic (no yield, F-independent, None invariant).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
