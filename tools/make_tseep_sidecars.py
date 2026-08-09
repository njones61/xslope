"""Build the transient seepage companions for the two docs transient samples.

A model that has been solved once keeps its results beside it, and a report of an
already-solved model reads them back rather than marching again. The two
transient samples shipped figures only, so a report of either had to re-solve
them or say nothing. This writes what a solved transient model carries:

    {stem}_mesh.json        the mesh the march was run on
    {stem}_tseep.csv        time, node_id, head, u -- one block per saved frame
    {stem}_tseep_meta.json  the frame ledger, the boundary rates, the mass
                            balance, the stage times and the provenance

The mesh is the one the sample's shipped figures are drawn on
(docs/seep/images/tseep_history.py and tools/make_earth_dam_tseep_figures.py
carry the same characteristic lengths), so a report of the model shows the same
run the documentation does.

Deterministic: same mesh, same march, same frames. Re-running reproduces every
file it writes.

    PYTHONPATH=. python3 tools/make_tseep_sidecars.py            # both
    PYTHONPATH=. python3 tools/make_tseep_sidecars.py johnson    # one
"""
import contextlib
import io
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

from xslope.fileio import load_slope_data                                # noqa: E402
from xslope.mesh import (build_mesh_from_polygons, export_mesh_to_json,  # noqa: E402
                         get_material_polygons)
from xslope.seep import (build_seep_data, build_tseep_data,              # noqa: E402
                         export_transient_solution, run_transient_seepage)

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FILES = os.path.join(REPO_ROOT, "docs", "seep", "files")

#: The characteristic length each sample is meshed at — the one its shipped
#: figures were drawn on.
SAMPLES = {
    "earth_dam": ("xslope_earth_dam_tseep.xlsx", 110.0 / 64.0),
    "johnson": ("xslope_johnson_res_tseep.xlsx", 10.0),
}


def _quiet(fn, *a, **k):
    with contextlib.redirect_stdout(io.StringIO()):
        return fn(*a, **k)


def build(name):
    xlsx, target = SAMPLES[name]
    path = os.path.join(FILES, xlsx)
    stem = os.path.splitext(path)[0]

    slope_data = load_slope_data(path)
    mesh = _quiet(build_mesh_from_polygons, get_material_polygons(slope_data),
                  target, "tri3")
    mesh_path = f"{stem}_mesh.json"
    export_mesh_to_json(mesh, mesh_path)

    seep_data = _quiet(build_seep_data, mesh, slope_data)
    solution = _quiet(run_transient_seepage, seep_data,
                      build_tseep_data(slope_data), verbose=False)
    csv_path, meta_path = _quiet(
        export_transient_solution, seep_data, solution, stem,
        input_file=os.path.basename(path),
        mesh_file=os.path.basename(mesh_path))

    print(f"{name}: {len(mesh['nodes'])} nodes, {len(solution['frames'])} frames, "
          f"converged={solution['converged']}, "
          f"closure={solution['mass_balance']['final_closure']:.3e}")
    for written in (mesh_path, csv_path, meta_path):
        print(f"  wrote {os.path.relpath(written, REPO_ROOT)}")


def main(argv):
    for name in (argv[1:] or list(SAMPLES)):
        build(name)


if __name__ == "__main__":
    main(sys.argv)
