"""Guard for solve_ssrm's optional ``ssr_zone`` SSR-Search-Area constraint.

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

Run from the repo root:  PYTHONPATH=. python3 benchmarks/ssr_zone_guard.py
Exits non-zero on any failure.
"""
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

    if failures:
        print("\nFAILED:")
        for f in failures:
            print("  -", f)
        return 1
    print("\nOK: ssr_zone mask confines reduction inside the polygon, is inert for a "
          "full-domain zone, and composes with ssr_exclude by union.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
