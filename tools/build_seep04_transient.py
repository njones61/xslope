"""Build Tutorial SEEP-4's time-varying rain model — the workbook, the march, and
the companions a solved model carries.

SEEP-4 solves the same dam twice under steady rain: dry weather, then a constant
1e-8 m/s falling on it forever. This builds the third run the page closes on, in
which the rain starts, holds and stops:

    docs/tutorials/files/xslope_dam_infiltration_storm.xlsx
    docs/tutorials/files/xslope_dam_infiltration_storm_mesh.json
    docs/tutorials/files/xslope_dam_infiltration_storm_tseep.csv
    docs/tutorials/files/xslope_dam_infiltration_storm_tseep_meta.json

It starts from the page's own completed steady file, so the section, the soil, the
reservoir, the drain and the three rain blocks are the ones the page builds — the
only edits are the ones the page's transient section asks the reader to make:

  * STORAGE. A transient run needs ``Ss`` and ``Sy`` on every material, and the
    source problem is steady and publishes neither. They are chosen from the tables
    in docs/seep/transient.md for the dam fill this problem describes: a compacted
    silt, ``Ss = 3e-4`` /m (the stiff-clay / dense-silt band) and ``Sy = 0.15`` (the
    silt band) — a plausible drainable porosity for such a fill, picked for a
    readable response time. The page says so: chosen, not measured. (``Sy`` is not
    corroborated by the material's van Genuchten pair, because XSLOPE *defines*
    theta_s - theta_r := ``Sy`` and the material carries no independent retention
    data; the table is the whole justification.)

  * TWO SERIES, NOT ONE. The three flux blocks do not carry one rate — the crest
    takes the vertical rain rate and the two 2:1 faces take it times cos θ = 2/√5,
    which is the projection the page spends a section establishing. A single series
    driving all three would apply the crest rate to the faces and put 10% more
    water into the dam than fell on it. So the schedule is written twice: ``storm``
    is the rain itself (the crest block) and ``storm_face`` is the same curve
    scaled by 2/√5 (both face blocks). One series drives two boundaries, which is
    what the tseep sheet is for.

  * THE SCHEDULE. Times are in SECONDS, because the model's declared Time unit is
    sec and k is in m/s; the day figures are the readable form the page quotes.
    Thirty dry days first (a model already at its steady dry answer must not
    drift), a 30-day ramp to the page's own 1e-8 m/s, a 140-day hold, a 30-day
    fall back to nothing, then 370 days of draining.

The march is run on the page's own mesh — tri3 at a 1.0 m target size, the 473
nodes and 833 triangles the page reports — and asserted against the head the page
locks with its ``tseep_head`` tag, so a change that moves the run fails here before
it reaches the page.

Deterministic: same mesh, same march, same frames. Re-running reproduces every
file it writes.

    PYTHONPATH=. python3 tools/build_seep04_transient.py
"""
import contextlib
import io
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np                                                       # noqa: E402

from xslope.fileio import (default_template_path, load_slope_data,       # noqa: E402
                           save_slope_data_to_xlsx)
from xslope.mesh import (build_mesh_from_polygons, export_mesh_to_json,  # noqa: E402
                         extract_size_regions, get_material_polygons)
from xslope.seep import (build_seep_data, build_tseep_data,              # noqa: E402
                         export_transient_solution, run_transient_seepage,
                         transient_frame_index)

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FILES = os.path.join(REPO_ROOT, "docs", "tutorials", "files")
BASE = "xslope_dam_infiltration.xlsx"
OUT = "xslope_dam_infiltration_storm.xlsx"

#: Seconds in a day. The model's Time unit is sec (its k is in m/s), so the whole
#: schedule is written in seconds; days are how the page reads it.
DAY = 86400.0

#: Storage for the compacted dam fill, from the tables in docs/seep/transient.md.
#: Ss is a 1/length and the model is in meters.
STORAGE = dict(Ss=3.0e-4, Sy=0.15)

#: The vertical rain rate at the crest of the storm — the page's own steady rate.
PEAK = 1.0e-8
#: cos(theta) for a 2:1 face: 2/sqrt(5). The page derives it; it is repeated here
#: rather than imported so the two series can never drift apart.
COS_FACE = 0.894427191

#: The storm, in days: dry, dry, full rain, full rain, dry. Linear between.
SCHEDULE_DAYS = [0.0, 30.0, 60.0, 200.0, 230.0]
SCHEDULE_RATE = [0.0, 0.0, PEAK, PEAK, 0.0]
DURATION_DAYS = 600.0
SAVE_INTERVAL_DAYS = 25.0

#: The mesh the page builds and every run on this model uses.
TARGET_SIZE = 1.0
ELEMENT_TYPE = "tri3"

#: The page's lock: total head at three interior stations in the peak frame, day
#: 200. Literals, matched to the tag in seep04_dam_infiltration.md.
LOCK_TIME_DAYS = 200.0
LOCK_POINTS = [(26.0, 8.0, 8.4242),      # crest centerline, just above the water table
               (26.0, 4.0, 8.0568),      # crest centerline, deep in the saturated zone
               (34.0, 2.0, 5.3144)]      # downstream, between the crest and the drain
LOCK_TOL = 0.05


def _quiet(fn, *a, **k):
    with contextlib.redirect_stdout(io.StringIO()):
        return fn(*a, **k)


def storm_model():
    """The page's completed steady model with the storm written into it: storage on
    the material, the two series on the tseep sheet, and the series names in the
    three flux values."""
    sd = _quiet(load_slope_data, os.path.join(FILES, BASE))
    for m in sd["materials"]:
        m.update(STORAGE)

    times = [d * DAY for d in SCHEDULE_DAYS]
    sd["tseep"] = dict(
        times=times,
        series={"storm": list(SCHEDULE_RATE),
                "storm_face": [q * COS_FACE for q in SCHEDULE_RATE]},
        duration=DURATION_DAYS * DAY,
        save_interval=SAVE_INTERVAL_DAYS * DAY,
        save_times=[],
    )

    # The three flux blocks, in the order the page enters them: the upstream face
    # above the waterline, the crest, the downstream face. The crest takes the
    # vertical rain; both faces take its projection.
    fluxes = sd["seepage_bc"]["specified_fluxes"]
    if len(fluxes) != 3:
        raise SystemExit(f"{BASE} carries {len(fluxes)} flux blocks, expected 3")
    for block, name in zip(fluxes, ("storm_face", "storm", "storm_face")):
        block["flux"] = name
    return sd


def march(path):
    """Mesh the saved workbook and run its transient march. Returns
    (mesh, mesh_path, seep_data, solution)."""
    sd = _quiet(load_slope_data, path)
    mesh = _quiet(build_mesh_from_polygons, get_material_polygons(sd),
                  TARGET_SIZE, ELEMENT_TYPE,
                  size_regions=extract_size_regions(sd))
    mesh_path = f"{os.path.splitext(path)[0]}_mesh.json"
    export_mesh_to_json(mesh, mesh_path)
    seep_data = _quiet(build_seep_data, mesh, sd)
    solution = _quiet(run_transient_seepage, seep_data, build_tseep_data(sd),
                      verbose=False)
    return mesh, mesh_path, seep_data, solution


def head_sampler(seep_data, head):
    """Total head at an (x, y) station — inverse-distance over the four nearest
    nodes, the same sampling run_tests' tseep_head rows use."""
    nodes = seep_data["nodes"]

    def at(xq, yq):
        d2 = (nodes[:, 0] - xq) ** 2 + (nodes[:, 1] - yq) ** 2
        idx = np.argsort(d2)[:4]
        w = 1.0 / np.maximum(d2[idx], 1e-12)
        return float(np.sum(w * head[idx]) / np.sum(w))
    return at


def check_lock(seep_data, solution):
    """Assert the three heads the page's tseep_head tag locks. Returns the sampled
    values so the manifest can print them."""
    fi = transient_frame_index(solution, LOCK_TIME_DAYS * DAY)
    at = head_sampler(seep_data,
                      np.asarray(solution["frames"][fi]["head"], dtype=float))
    got = []
    for x, y, want in LOCK_POINTS:
        h = at(x, y)
        got.append(h)
        if abs(h - want) > LOCK_TOL:
            raise SystemExit(
                f"LOCK FAILED at ({x:g}, {y:g}) on day {LOCK_TIME_DAYS:g}: "
                f"the page locks h = {want:.4f} m, this run gives {h:.4f} m "
                f"(tolerance {LOCK_TOL:g}). Either the physics moved or the page's "
                f"tag is stale — decide which before editing either.")
    return got


def build():
    out = os.path.join(FILES, OUT)
    save_slope_data_to_xlsx(storm_model(), out, template=default_template_path())
    mesh, mesh_path, seep_data, solution = march(out)
    csv_path, meta_path = _quiet(
        export_transient_solution, seep_data, solution,
        os.path.splitext(out)[0],
        input_file=os.path.basename(out),
        mesh_file=os.path.basename(mesh_path))
    locked = check_lock(seep_data, solution)

    frames = solution["frames"]
    pi = transient_frame_index(solution, LOCK_TIME_DAYS * DAY)
    first, peak, last = frames[0], frames[pi], frames[-1]
    # Mass balance is reported at the PEAK, not at the end: by day 600 the section
    # has given back everything it took in, so both the stored change and the net
    # cumulative inflow are near zero and their ratio is meaningless.
    mb = solution["mass_balance"]["per_frame"][pi]

    print("SEEP-4 storm model")
    print(f"  mesh          {len(mesh['nodes'])} nodes, "
          f"{len(mesh['elements'])} {ELEMENT_TYPE} elements at {TARGET_SIZE:g} m")
    print(f"  march         {len(frames)} saved frames over "
          f"{DURATION_DAYS:g} days, converged={solution['converged']}")
    print(f"  mass balance  at day {LOCK_TIME_DAYS:g}: stored "
          f"{mb['stored_change']:.4f} m3/m against net inflow "
          f"{mb['cumulative_inflow']:.4f} m3/m "
          f"({100 * mb['closure']:.1f}% apart)")
    print("  drain outflow (m3/s per m)")
    for label, fr in (("day 0, dry", first),
                      (f"day {LOCK_TIME_DAYS:g}, peak", peak),
                      (f"day {DURATION_DAYS:g}, drained", last)):
        print(f"    {label:22s} {fr['outflow']:.4e}   "
              f"(in {fr['inflow']:.4e})")
    print(f"  locked heads at day {LOCK_TIME_DAYS:g}")
    for (x, y, want), h in zip(LOCK_POINTS, locked):
        print(f"    ({x:g}, {y:g})  {h:.4f} m   (page locks {want:.4f}, "
              f"tol {LOCK_TOL:g})")
    print("  files")
    for written in (out, mesh_path, csv_path, meta_path):
        print(f"    {os.path.relpath(written, REPO_ROOT)}")


if __name__ == "__main__":
    build()
