# Copyright 2025 Norman L. Jones
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Divergence fence for the optional compiled Mohr-Coulomb SSRM kernel.

The pure-NumPy Step-6 path in ``fem.solve_fem`` is the oracle: every locked
factor of safety is defined by it. The compiled kernel (``fast_kernel=True``,
built with ``setup_kernel.py``) must reproduce it bit-for-bit. This module solves
three small, coarse cases -- one plain Mohr-Coulomb, one with the Rankine tension
cutoff, one with matric suction -- both ways and fails on:

  * any factor-of-safety disagreement between the two paths, or
  * any displacement-field max-abs difference above 1e-8 after a fixed number of
    viscoplastic iterations at fixed F (the sensitive kernel-identity signal;
    bit-faithful kernels land near 1e-14).

ROLE: this is the drift fence between the compiled kernel and the reference path,
and it is what makes the suite's fast-first-with-fallback mode sound. Without it,
a reference-only physics edit (fem.py's Step-6 block changed, the kernel left
untouched) would still pass every lock, because fast-first would happily run the
unchanged, now-stale kernel and never notice the two paths had diverged. This
gate catches that disagreement directly, on every run, instead of relying on an
unrelated lock to surface it.

LIMIT: its cases are small and well-conditioned -- they verify agreement, not
knife-edge behavior. RS2-62c-class divergence (see xslope/_fem_kernel.pyx's
KNOWN LIMIT: a near-bifurcation mechanism where floating-point re-association
flips a bisection verdict) is caught by the fast-first-with-fallback protocol's
reference re-check, not by this gate staying green. This module is a required
companion to fast-first-with-fallback -- remove it only if fast-first is also
removed.

It is CHEAP by construction (coarse meshes, short iteration budgets) so it can run
as a routine CI gate in well under two minutes. When the compiled kernel is not
built it reports "unavailable" and the caller skips it -- the NumPy path needs no
cross-check against itself.

Run standalone:  PYTHONPATH=. python3 benchmarks/kernel_xcheck.py
Exits non-zero on any disagreement; exits 0 (with a notice) when the kernel is
not built.
"""
import os
import sys
import warnings

import numpy as np

warnings.filterwarnings("ignore")

# Field-difference ceiling. Bit-faithful kernels sit ~1e-14 after a few hundred
# iterations; 1e-8 is a wide fence that still catches any real divergence.
FIELD_TOL = 1e-8

# Base directory: repo root is this file's parent's parent.
_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# (name, file relative to repo root, element_type, target_size, solve_fem extras)
_CASES = [
    ("MC tri6", "docs/fem/files/xslope_griffiths1.xlsx", "tri6", 12.0, {}),
    ("tension-cutoff tri6", "docs/fem/files/xslope_griffiths1.xlsx", "tri6", 12.0,
     {"tension_cutoff": True}),
    ("suction tri6", "docs/fem/files/xslope_griffiths5_0p7.xlsx", "tri6", 11.0,
     {"suction_phi_b": {"soil": 15.0}, "suction_cap": 200.0}),
]


def kernel_available():
    """True when the compiled kernel is importable (built via setup_kernel.py)."""
    try:
        from xslope import _fem_kernel  # noqa: F401
        return True
    except ImportError:
        return False


def files_present():
    """True when the case input files are on disk (absent in engine-only clones)."""
    return all(os.path.exists(os.path.join(_ROOT, c[1])) for c in _CASES)


def _build(file_rel, element_type, target_size):
    import xslope.fem as fem
    from xslope.fileio import load_slope_data
    from xslope.mesh import (get_material_polygons, build_mesh_from_polygons,
                             extract_constraint_line_geometry, extract_point_constraints)
    sd = load_slope_data(os.path.join(_ROOT, file_rel))
    lines, _nr, _np = extract_constraint_line_geometry(sd)
    polys = get_material_polygons(sd, reinf_lines=lines)
    mesh = build_mesh_from_polygons(polys, target_size=target_size,
                                    element_type=element_type, lines=lines,
                                    point_constraints=extract_point_constraints(sd))
    return fem.build_fem_data(sd, mesh)


def check():
    """Cross-check every case both paths. Returns a list of failure strings (empty
    on full agreement). Assumes kernel_available() and files_present() are True."""
    import xslope.fem as fem
    failures = []

    for name, file_rel, et, ts, extra in _CASES:
        fd = _build(file_rel, et, ts)

        # (1) Kernel-identity: fixed F, fixed iteration count, deterministic. Both
        # paths must integrate the same viscoplastic field to ~machine precision.
        common = dict(max_iterations=250, tolerance=1e-12, force_tol=1e-30,
                      max_disp_factor=None, early_exit=False)
        s_ref = fem.solve_fem(fd, F=1.3, fast_kernel=False, **common, **extra)
        s_fast = fem.solve_fem(fd, F=1.3, fast_kernel=True, **common, **extra)
        d_ref = np.asarray(s_ref["displacements"])
        d_fast = np.asarray(s_fast["displacements"])
        dmax = float(np.max(np.abs(d_ref - d_fast)))
        if not (dmax <= FIELD_TOL):
            failures.append(f"{name}: field max|du|={dmax:.2e} > {FIELD_TOL:.0e}")

        # (2) FS identity: a short SSRM both ways must agree exactly. Narrow window
        # + loose tolerance keeps the bisection to a couple of trials (cheap).
        ssrm = dict(F_min=1.2, F_max=1.6, tolerance=0.05, max_iterations=4000,
                    debug_level=0)
        with _force_kernel(fem, False):
            r_ref = fem.solve_ssrm(fd, **ssrm, **extra)
        with _force_kernel(fem, True):
            r_fast = fem.solve_ssrm(fd, **ssrm, **extra)
        fs_ref, fs_fast = r_ref.get("FS"), r_fast.get("FS")
        if fs_ref is None or fs_fast is None:
            failures.append(f"{name}: SSRM did not return an FS "
                            f"(ref={fs_ref}, fast={fs_fast})")
        elif fs_ref != fs_fast:
            failures.append(f"{name}: FS_ref={fs_ref:.6f} != FS_fast={fs_fast:.6f} "
                            f"(dFS={abs(fs_ref - fs_fast):.2e})")

    return failures


class _force_kernel:
    """Context manager: force fem.solve_fem's fast_kernel flag for the nested
    solve_ssrm trials (solve_ssrm calls solve_fem internally without the flag)."""
    def __init__(self, fem_mod, on):
        self._fem = fem_mod
        self._on = on
        self._orig = fem_mod.solve_fem

    def __enter__(self):
        orig, on = self._orig, self._on

        def _wrap(*a, **k):
            k["fast_kernel"] = on
            return orig(*a, **k)
        self._fem.solve_fem = _wrap
        return self

    def __exit__(self, *exc):
        self._fem.solve_fem = self._orig
        return False


def main():
    if not kernel_available():
        print("kernel_xcheck: compiled xslope._fem_kernel not built; nothing to "
              "cross-check (build it with `python setup_kernel.py build_ext "
              "--inplace`). SKIP.")
        return 0
    if not files_present():
        print("kernel_xcheck: case input files not present (engine-only clone). SKIP.")
        return 0
    failures = check()
    if failures:
        print("kernel_xcheck: FAIL")
        for f in failures:
            print("  -", f)
        return 1
    print(f"kernel_xcheck: PASS ({len(_CASES)} cases, both paths agree; "
          f"field tol {FIELD_TOL:.0e})")
    return 0


if __name__ == "__main__":
    sys.path.insert(0, _ROOT)
    raise SystemExit(main())
