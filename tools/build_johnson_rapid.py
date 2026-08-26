"""Build the Johnson Reservoir rapid-drawdown workbooks from one parameter block.

Four files come out of the single set of constants below, so the rapid-drawdown
family cannot drift apart:

  ``docs/lem/files/xslope_johnson_res_rapid.xlsx``
      The worked example of ``docs/lem/rapid.md`` ("Rapid Drawdown from a Transient
      Solution"): the transient route alone, no second boundary set, no piezometric
      lines.

  ``docs/tutorials/files/xslope_johnson_rapid_start.xlsx``
      The starter file of tutorial COMBO-2. Strengths, the undrained core's d/psi,
      ``Water loads: auto`` and the coarse **piezometric pair** — and nothing else:
      no seepage boundary sets and no transient schedule, because building those is
      what the tutorial does. Every material takes ``u = piezo``.

  ``docs/tutorials/files/xslope_johnson_rapid.xlsx``
      The completed COMBO-2 model: the same section with **both** seepage boundary
      sets (full pool and drawn down), the transient pool schedule with its
      ``stage_1`` / ``stage_2`` times, and ``u = seep`` on every material. The
      piezometric pair is deleted: the boundary sets outrank it for water loads and
      ``u = seep`` outranks it for pore pressure, so it would feed nothing. No
      ``_seep.csv`` / ``_seep2.csv`` companions ship with it: a steady run writes
      them, and a transient run stages its two frames in memory.

  ``docs/tutorials/files/xslope_johnson_fs_time.xlsx`` (+ companions)
      The model tutorial COMBO-3's Part 2 opens: the completed model with boundary
      set 2 removed, shipped with the mesh and the march already on it as
      ``_mesh.json``, ``_tseep.csv`` and ``_tseep_meta.json``. Part 2 sweeps the
      march rather than producing it, so the reader opens a dam that is already
      meshed and marched and goes straight to the sweep. It carries a base name of
      its own because a companion mesh and march sitting next to
      ``xslope_johnson_rapid.xlsx`` would load themselves when COMBO-2 opens that
      file, and COMBO-2 has the reader build both. Set 2 is dropped because a
      transient march reads neither of its stages from it, and left on the file it
      raises ``rapid.stage2_bc_ignored`` on every run of the sweep.

Everything is built deterministically from the committed transient **seepage**
sample ``docs/seep/files/xslope_johnson_res_tseep.xlsx``, so the cross-section,
zones, conductivities, storage properties, saved-frame schedule and stage times
can never diverge from the seepage side of the same dam. That file is never
modified. The one value the tutorial pair overrides is the level the reservoir is
lowered to: the tutorial's drawdown stops at a residual pool 10 ft deep, while the
worked example of ``docs/lem/rapid.md`` keeps the base file's total drawdown to
the tailwater datum.

Run:  PYTHONPATH=. python3 tools/build_johnson_rapid.py           # all four
      PYTHONPATH=. python3 tools/build_johnson_rapid.py fs_time   # the solved set

The solved set carries the transient march, so it takes a few minutes; the other
three are written in seconds.
"""

from __future__ import annotations

import contextlib
import copy
import io
import os
import sys

from xslope.fileio import (load_slope_data, save_slope_data_to_xlsx,
                           default_template_path)

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SEEP_FILES = os.path.join(REPO_ROOT, "docs", "seep", "files")
LEM_FILES = os.path.join(REPO_ROOT, "docs", "lem", "files")
TUT_FILES = os.path.join(REPO_ROOT, "docs", "tutorials", "files")

BASE = os.path.join(SEEP_FILES, "xslope_johnson_res_tseep.xlsx")

# --- the one parameter block ------------------------------------------------

#: Undrained core, Kc = 1 (CU) envelope, below the drained c' = 400 / phi' = 18.
#: The free-draining sand shell and the silty-sand foundation carry no d/psi and
#: are analyzed drained through the drawdown.
CORE_D, CORE_PSI = 250.0, 14.0

#: The critical upstream circle, located by a rapid-drawdown circular search over
#: the transient-staged fields. Yo - Depth = R keeps the tangent elevation
#: consistent. It toes near the upstream base, daylights just past the crest, and
#: cuts twelve of forty slices through the core.
CIRCLE = dict(Xo=275.0, Yo=235.0, R=160.0, Depth=75.0)

#: The pool pair. The full pool is the reservoir at elevation 160 that the base
#: file's `pool` series starts from. The drawdown stops 10 ft above the
#: elevation-100 tailwater datum: a residual pool at elevation 110 stands against
#: the upstream toe at the end of it, which is the level the second seepage
#: boundary set states as a steady problem and the level the tutorial's `pool`
#: series ramps down to.
POOL_FULL, POOL_DOWN = 160.0, 110.0

#: The drawn-down steady boundary set: the residual pool along the upstream
#: foreshore and up the face to elevation 110, the tailwater at 100, and the
#: downstream slope as an exit face. Ten feet of head still cross the section, so
#: this is an unconfined flow problem with a phreatic surface of its own, stated
#: exactly the way set 1 states the full-pool one.
BC2 = {
    "specified_heads": [
        {"head": POOL_DOWN,
         "coords": [(0.0, 100.0), (200.0, 100.0), (220.0, 110.0)], "kind": "head"},
        {"head": 100.0, "coords": [(550.0, 100.0), (750.0, 100.0)], "kind": "head"},
    ],
    "specified_fluxes": [],
    "exit_face": [(380.0, 180.0), (550.0, 100.0)],
}

#: The tutorial's reservoir schedule: full pool held for five days, then lowered
#: to the residual pool over the following 45. Only the last breakpoint differs
#: from the base seepage sample's series, whose times and stage instants are kept.
POOL_SCHEDULE = [POOL_FULL, POOL_FULL, POOL_DOWN]

#: The coarse piezometric pair. These are five- and six-point polylines a designer
#: would sketch from piezometer readings, not traces of a solved field: Line 1
#: mimics the full-pool phreatic surface of boundary set 1, Line 2 the drawn-down
#: one of boundary set 2 — the residual pool at 110 upstream, falling through the
#: core and running out to tailwater at 100. Both were read off the solved fields
#: and rounded. Line 1 leaves the upstream shell AT the reservoir level and Line 2
#: at the residual pool level, so the pools they describe are the ones the two
#: seepage boundary sets state; a line that came off the shell a foot low would put
#: two different pools on one model. Both tails sit on the ground surface from
#: x = 550 out, so neither invents a pond on the downstream foreshore.
PIEZO_1 = [(0.0, 160.0), (360.0, 160.0), (410.0, 120.0), (550.0, 100.0),
           (750.0, 100.0)]
PIEZO_2 = [(0.0, 110.0), (220.0, 110.0), (360.0, 107.0), (410.0, 103.0),
           (550.0, 100.0), (750.0, 100.0)]

#: The mesh the solved file ships, stated the way Studio's Build Mesh dialog states
#: it: linear triangles, auto-sized from the geometry at 100 divisions across the
#: 750 ft section. The tutorials quote the 2,080 nodes and 3,923 triangles that come
#: out of it. Neither refinement path acts on this section — it carries no
#: constraint line, no size region, and no zone thin enough for the thin-zone pass —
#: so the dialog's defaults reduce to a plain target size here.
MESH_ELEMENT_TYPE = "tri3"
MESH_SIZE_DIVISIONS = 100


def _tutorial_schedule(sd):
    """The base file's transient schedule with the pool ramped to the residual
    level instead of to the tailwater datum. The saved-frame schedule and both
    stage times are the base file's own."""
    ts = copy.deepcopy(sd["tseep"])
    ts["series"]["pool"] = list(POOL_SCHEDULE)
    return ts


def _base():
    """The seepage sample with the two rapid-drawdown additions every file shares:
    the undrained core and the critical upstream circle."""
    sd = load_slope_data(BASE)
    sd["materials"][1]["d"] = CORE_D
    sd["materials"][1]["psi"] = CORE_PSI
    # max_depth is left at the base file's value: the fixed circle carries its own
    # tangent Depth, and overriding max_depth would clip the domain through the core.
    sd["circular"] = True
    sd["circles"] = [dict(CIRCLE)]
    return sd


def _save(sd, path):
    save_slope_data_to_xlsx(sd, path, template=default_template_path())
    print(f"wrote {os.path.relpath(path, REPO_ROOT)}")
    return path


def build_lem_worked_example():
    """``docs/lem/rapid.md``'s transient worked example — the base file plus the
    undrained core and the circle, and nothing else."""
    return _save(_base(), os.path.join(LEM_FILES, "xslope_johnson_res_rapid.xlsx"))


def build_tutorial_completed():
    """COMBO-2's completed model: both boundary sets, the schedule, u = seep.

    No piezometric lines. Once the boundary sets exist they state where the water
    stands for every load, and ``u = seep`` states every pore pressure, so a line
    left on the file would be read by nothing (``water.water_line_for_stage``:
    seepage head boundaries wherever a seepage analysis is defined, otherwise the
    piezometric line). The starter file keeps the pair, which is what the tutorial's
    first run is made from.

    Boundary set 2 IS kept. Part 2 of the tutorial runs the two-steady route on it,
    and Part 3 has the reader clear it before the transient march — the same pattern
    the piezometric pair follows one part earlier, each input deleted at the point
    the next route stops reading it."""
    sd = _base()
    sd["tseep"] = _tutorial_schedule(sd)
    sd["seepage_bc2"] = copy.deepcopy(BC2)
    sd["has_seepage_bc2"] = True
    sd["piezo_line"] = []
    sd["piezo_line2"] = []
    for m in sd["materials"]:
        m["u"] = "seep"
    return _save(sd, os.path.join(TUT_FILES, "xslope_johnson_rapid.xlsx"))


def build_tutorial_starter():
    """COMBO-2's starter: the piezometric pair only — no boundary sets, no
    schedule, u = piezo."""
    sd = _base()
    sd["piezo_line"] = list(PIEZO_1)
    sd["piezo_line2"] = list(PIEZO_2)
    for m in sd["materials"]:
        m["u"] = "piezo"
    empty = {"specified_heads": [], "specified_fluxes": [], "exit_face": []}
    sd["seepage_bc"] = copy.deepcopy(empty)
    sd["seepage_bc2"] = copy.deepcopy(empty)
    sd["has_seepage_bc2"] = False
    sd["tseep"] = None
    return _save(sd, os.path.join(TUT_FILES, "xslope_johnson_rapid_start.xlsx"))


def _quiet(fn, *a, **k):
    with contextlib.redirect_stdout(io.StringIO()):
        return fn(*a, **k)


def build_tutorial_fs_time():
    """COMBO-3 Part 2's model: the completed model without boundary set 2, meshed
    and marched.

    Set 2 goes because a transient rapid drawdown reads stage 1 and stage 2 out of
    the march, so a second boundary set states nothing the run consults and the
    checks say so on every sweep. COMBO-2's Part 3 has the reader clear it by hand
    at this same point; the solved file arrives with it already cleared.

    The mesh and the march are the ones COMBO-2 builds — tri3, auto-sized at 100
    divisions, and the schedule's own twelve saved frames — so Part 2's numbers are
    computed on exactly the state COMBO-2 leaves behind."""
    from xslope.mesh import (build_mesh_from_polygons, export_mesh_to_json,
                             get_material_polygons)
    from xslope.seep import (build_seep_data, build_tseep_data,
                             export_transient_solution, run_transient_seepage)

    sd = _base()
    sd["tseep"] = _tutorial_schedule(sd)
    sd["seepage_bc2"] = {"specified_heads": [], "specified_fluxes": [],
                         "exit_face": []}
    sd["has_seepage_bc2"] = False
    sd["piezo_line"] = []
    sd["piezo_line2"] = []
    for m in sd["materials"]:
        m["u"] = "seep"
    path = _save(sd, os.path.join(TUT_FILES, "xslope_johnson_fs_time.xlsx"))
    stem = os.path.splitext(path)[0]

    # Re-read what was written, so the companions are built on the model as the file
    # states it rather than on the dict that produced it.
    sd = load_slope_data(path)
    xs = [x for x, _ in sd["ground_surface"].coords]
    target = (max(xs) - min(xs)) / MESH_SIZE_DIVISIONS
    mesh = _quiet(build_mesh_from_polygons, get_material_polygons(sd), target,
                  MESH_ELEMENT_TYPE)
    mesh_path = f"{stem}_mesh.json"
    export_mesh_to_json(mesh, mesh_path)

    seep_data = _quiet(build_seep_data, mesh, sd)
    solution = _quiet(run_transient_seepage, seep_data, build_tseep_data(sd),
                      verbose=False)
    csv_path, meta_path = _quiet(
        export_transient_solution, seep_data, solution, stem,
        input_file=os.path.basename(path),
        mesh_file=os.path.basename(mesh_path))

    print(f"  mesh: {len(mesh['nodes'])} nodes, {len(mesh['elements'])} elements "
          f"(target size {target:g})")
    print(f"  march: {len(solution['frames'])} frames, "
          f"converged={solution['converged']}, "
          f"closure={solution['mass_balance']['final_closure']:.3e}")
    for written in (mesh_path, csv_path, meta_path):
        print(f"  wrote {os.path.relpath(written, REPO_ROOT)} "
              f"({os.path.getsize(written) / 1024:.0f} KB)")
    return path


def build(which=None):
    if which in ("fs_time", "solved"):
        return [build_tutorial_fs_time()]
    return [build_lem_worked_example(),
            build_tutorial_starter(),
            build_tutorial_completed(),
            build_tutorial_fs_time()]


if __name__ == "__main__":
    build(sys.argv[1] if len(sys.argv) > 1 else None)
