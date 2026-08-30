"""Build the SOLVED companion of COMBO-1's model — the same workbook, with its
mesh, its seepage field and its strength-reduction solution beside it.

``docs/tutorials/combo01_seepage_stability.md`` ships
``files/xslope_johnson_res.xlsx`` carrying no mesh and no solution, because the
page's whole subject is making the three runs. A reader who only wants to LOOK at
the finished model — or a later tutorial that opens it to work on something else —
should not have to sit through a strength reduction first. This writes that second
copy:

    xslope_johnson_res_solved.xlsx              the model, byte-identical inputs
    xslope_johnson_res_solved_mesh.json         the one mesh all three runs share
    xslope_johnson_res_solved_seep.csv          the steady field, BC set 1
    xslope_johnson_res_solved_fem_nodes.csv     the converged strength-reduction field
    xslope_johnson_res_solved_fem_elements.csv
    xslope_johnson_res_solved_fem_meta.json     the run record: bracket, tolerance,
                                                criterion, trials, FS
    xslope_johnson_res_solved_fem_failure_*      the at-failure mechanism

Those are exactly the names ``studio.main_window._load_solution_sidecars`` looks
for beside a workbook, so opening the solved copy in Studio attaches the mesh, the
seepage solution and the finite element solution with no solver run. The limit
equilibrium answer is NOT among them: no sidecar carries an LEM solution, so the
Spencer search is the one run a reader of the solved file still makes. It is
checked here all the same, because a file that ships as "the solved COMBO-1 model"
has to be the model whose published numbers are COMBO-1's.

THE PAGE'S OWN TAGS ARE THE SETTINGS. Element type, mesh divisions, method, slice
count, bracket and tolerance are read from the ``<!-- test: -->`` tags on
``combo01_seepage_stability.md``, and each run has to land on the number that tag
locks, within that tag's own tolerance, or nothing is written:

    seepage flow rate   1.925 ft^3/day/ft  +/- 0.005
    Spencer (u = seep)  FS = 1.248         +/- 0.005
    Spencer (u = none)  FS = 1.618         +/- 0.005   (the page's contrast run)
    SSRM                FS = 1.2305        +/- 0.01

Young's modulus and Poisson's ratio are the file's own — the materials table
carries E and nu on all three rows (shell 2.0e6 / 0.20, core 2.0e5 / 0.40,
foundation 1.0e6 / 0.35) — so nothing about the model is stated a second time here.

Deterministic: the same mesh, the same boundary conditions, the same bracket. The
strength reduction runs on the pure-NumPy reference kernel, the path the locks are
defined on. Roughly 15 minutes, nearly all of it the strength reduction.

    PYTHONPATH=. python3 tools/build_combo01_solved.py            # solve and write
    PYTHONPATH=. python3 tools/build_combo01_solved.py --check    # solve, write nothing
"""
import argparse
import contextlib
import io
import os
import shutil
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import run_tests as RT                                                   # noqa: E402
from xslope.fileio import load_slope_data                                # noqa: E402
from xslope.package import FEM_SOLUTION_SIDECARS, MESH_SIDECAR, SEEP_SIDECARS  # noqa: E402

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TUT_DIR = os.path.join(REPO_ROOT, "docs", "tutorials")
PAGE = os.path.join(TUT_DIR, "combo01_seepage_stability.md")
SOURCE = os.path.join(TUT_DIR, "files", "xslope_johnson_res.xlsx")
SOLVED = os.path.join(TUT_DIR, "files", "xslope_johnson_res_solved.xlsx")

#: The seepage dialog's defaults, which the page leaves alone and names in its
#: text: "Leave Convergence tol at 0.0001 and Max iterations at 400."
SEEP_SETTINGS = dict(tol=1e-4, max_iter=400)


def _quiet(fn, *a, **k):
    with contextlib.redirect_stdout(io.StringIO()):
        return fn(*a, **k)


def tags():
    """The page's four locks, keyed by the run each one describes."""
    found = {}
    for t in RT.parse_test_tags(PAGE):
        kind = t.get("type")
        if kind == "seep":
            found["seep"] = t
        elif kind == "circular_search":
            found["lem_none" if t.get("u_option") == "none" else "lem"] = t
        elif kind == "fem_ssrm":
            found["ssrm"] = t
    missing = {"seep", "lem", "lem_none", "ssrm"} - set(found)
    if missing:
        raise SystemExit(f"{os.path.basename(PAGE)} carries no tag for {sorted(missing)}")
    return found


def check(name, got, tag, key="expected_fs"):
    """One published number against the tag that locks it."""
    want, tol = float(tag[key]), float(tag["tolerance"])
    ok = abs(got - want) <= tol
    print(f"  {name:<22} {got:.4f}   tag {want:g} +/- {tol:g}   "
          f"{'ok' if ok else 'MISMATCH'}")
    if not ok:
        raise SystemExit(f"{name} came out {got:.4f} against the page's {want:g} "
                         f"+/- {tol:g} — nothing written")


def clear_sidecars(stem):
    """Remove any companion of a previous build, so the solve starts from the
    workbook alone and cannot read last run's mesh or field back in."""
    for suffix in (MESH_SIDECAR,) + SEEP_SIDECARS + FEM_SOLUTION_SIDECARS:
        path = f"{stem}{suffix}"
        if os.path.exists(path):
            os.unlink(path)


def build_mesh(slope_data, tag):
    """The mesh Build Mesh makes at the page's settings — quadratic triangles at
    the section width over the tag's size divisions. The same call
    ``studio.runners.MeshWorker.build`` makes: this model declares no constraint
    line, no size region and no thin zone, so the dialog's own defaults and this
    reduce to one mesh."""
    from xslope.mesh import (build_mesh_from_polygons,
                             extract_constraint_line_geometry,
                             extract_point_constraints, extract_size_regions,
                             get_material_polygons)
    lines, _n_reinf, _n_pile = extract_constraint_line_geometry(slope_data)
    polygons = get_material_polygons(slope_data, reinf_lines=lines)
    return _quiet(build_mesh_from_polygons, polygons,
                  target_size=RT._tag_target_size(tag, slope_data),
                  element_type=tag["element_type"], lines=lines or None,
                  element_size_1d=slope_data.get("element_size_1d"),
                  point_constraints=extract_point_constraints(slope_data),
                  size_regions=extract_size_regions(slope_data))


def solve_seep(mesh, slope_data):
    from xslope.seep import (apply_steady_stability_field, build_seep_data,
                             run_seepage_analysis)
    seep_data = _quiet(build_seep_data, mesh, slope_data, seep_bc=1)
    solution = _quiet(run_seepage_analysis, seep_data, **SEEP_SETTINGS)
    if not solution.get("converged", True):
        raise SystemExit("the steady seepage solve did not converge — nothing written")
    # The field belongs to the model the moment it exists: this is what puts it
    # where a u = 'seep' stability run reads it, and it is the same call Studio's
    # seepage runner makes.
    apply_steady_stability_field(slope_data, solution, bc=1, verbose=False)
    return seep_data, solution


def spencer(slope_data, tag):
    """The page's search: Spencer, Auto search, the tag's slice count."""
    from xslope.search import circular_search
    fs_cache, _converged, _path, _circles = _quiet(
        circular_search, slope_data, tag["method"],
        num_slices=int(tag.get("num_slices", 30)))
    if not fs_cache or fs_cache[0]["FS"] >= 9999:
        raise SystemExit("the circular search found no valid surface")
    return fs_cache[0]["FS"], len(fs_cache)


def solve_ssrm(slope_data, mesh, tag):
    """The strength reduction at the page's bracket, on the seepage run's own mesh,
    with the at-failure mechanism captured (Studio's default, and what the
    ``_fem_failure_*`` companions carry)."""
    import xslope.fem as fem
    fem_data = _quiet(fem.build_fem_data, slope_data, mesh)
    options = {"F_min": float(tag["f_min"]), "F_max": float(tag["f_max"]),
               "tolerance": float(tag["tolerance"])}
    with RT._force_fast_kernel(fem, False):
        result = _quiet(fem.solve_ssrm, fem_data, debug_level=0,
                        capture_failure_state=True, **options)
    if not result.get("converged"):
        raise SystemExit(f"SSRM did not converge: {result.get('error')}")
    return fem_data, result, options


def write(stem, mesh, seep_data, seep_solution, fem_data, result, options):
    """Every companion, written the way the writer that owns each one writes it."""
    from xslope.fem import export_fem_solution, ssrm_run_record
    from xslope.mesh import export_mesh_to_json
    from xslope.seep import export_seep_solution

    export_mesh_to_json(mesh, f"{stem}{MESH_SIDECAR}")
    _quiet(export_seep_solution, seep_data, seep_solution, f"{stem}_seep.csv")

    field = result["last_solution"]
    # The meta Studio's own writer would put beside the run: the run's record of
    # itself (bracket, tolerance, criterion, final interval, per-trial verdicts),
    # then the three facts the record does not name — the factor of safety, what
    # was run, and the strength reduction factor the WRITTEN FIELD was solved at.
    # F is not FS, and _restore_fem_sidecar titles the panel from it.
    meta = ssrm_run_record(result, fem_data, options)
    meta["FS"] = result["FS"]
    meta["analysis"] = "ssrm"
    meta["F"] = field.get("F")
    _quiet(export_fem_solution, fem_data, field, stem, meta=meta,
           failure_solution=result.get("failure_solution"))


def manifest(stem):
    """What was written, and how much of it."""
    folder, base = os.path.dirname(stem), os.path.basename(stem)
    names = sorted(n for n in os.listdir(folder)
                   if n == base + ".xlsx" or n.startswith(base + "_"))
    total = 0
    print("\nWrote:")
    for name in names:
        size = os.path.getsize(os.path.join(folder, name))
        total += size
        print(f"  {name:<50} {size / 1024:>9.1f} KB")
    print(f"  {'TOTAL':<50} {total / 1024 ** 2:>9.2f} MB")
    return total


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--check", action="store_true",
                    help="run everything and verify the locks, write nothing")
    args = ap.parse_args(argv)

    t0 = time.time()
    tag = tags()
    stem = os.path.splitext(SOLVED)[0]

    if not args.check:
        clear_sidecars(stem)
        shutil.copyfile(SOURCE, SOLVED)     # the same model, byte for byte
        with open(SOURCE, "rb") as a, open(SOLVED, "rb") as b:
            if a.read() != b.read():
                raise SystemExit("the copy is not byte-identical to the source")
        print(f"Copied {os.path.basename(SOURCE)} -> {os.path.basename(SOLVED)}")

    # Solved from the copy where there is one, so the file that ships is the file
    # the numbers came from; its companions were just cleared, so it is read bare.
    model = SOLVED if (not args.check and os.path.exists(SOLVED)) else SOURCE
    slope_data = _quiet(load_slope_data, model)

    mesh = build_mesh(slope_data, tag["seep"])
    slope_data["mesh"] = mesh               # the field and the FEM run read it here
    print(f"Mesh: {tag['seep']['element_type']}, "
          f"{RT._tag_target_size(tag['seep'], slope_data):g} ft target, "
          f"{len(mesh['nodes'])} nodes, {len(mesh['elements'])} elements")

    seep_data, seep_solution = solve_seep(mesh, slope_data)
    print("\nAgainst the page:")
    check("seepage flow rate", seep_solution["flowrate"], tag["seep"],
          key="expected_flowrate")

    fs_lem, n_circles = spencer(slope_data, tag["lem"])
    check("Spencer, u = seep", fs_lem, tag["lem"])

    # The page's contrast run: the same mesh and the same solved field, with the
    # materials' pore-pressure column turned off. It writes nothing — it is a check
    # that the shipped model is the one whose 1.618 the page publishes.
    dry = _quiet(load_slope_data, model)
    dry["mesh"] = mesh
    for m in dry.get("materials", []):
        m["u"] = "none"
    fs_dry, _n = spencer(dry, tag["lem_none"])
    check("Spencer, u = none", fs_dry, tag["lem_none"])

    print(f"\nStrength reduction (reference kernel, bracket "
          f"[{tag['ssrm']['f_min']:g}, {tag['ssrm']['f_max']:g}])…")
    t1 = time.time()
    fem_data, result, options = solve_ssrm(slope_data, mesh, tag["ssrm"])
    lo, hi = result.get("final_interval", (float("nan"),) * 2)
    print(f"  {len(result['trials'])} trials, final bracket [{lo:.4f}, {hi:.4f}], "
          f"{time.time() - t1:.0f}s")
    check("SSRM", result["FS"], tag["ssrm"])

    if args.check:
        print(f"\nEvery published number reproduces. Nothing written (--check). "
              f"{time.time() - t0:.0f}s")
        return 0

    write(stem, mesh, seep_data, seep_solution, fem_data, result, options)
    manifest(stem)
    print(f"\n{time.time() - t0:.0f}s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
