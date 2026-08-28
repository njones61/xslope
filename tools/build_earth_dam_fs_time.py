"""Build Tutorial COMBO-3 Part 1's pair of models — the cored earth dam of Tutorial
SEEP-3, meshed and marched, without the strengths a stability run reads
(``xslope_earth_dam_fs_time_start.xlsx``, the file the page opens) and with them
(``xslope_earth_dam_fs_time.xlsx``, the file it ends on).

COMBO-3 puts a factor-of-safety-versus-time curve on the drawdown SEEP-3 solved,
and it needs no new seepage input to do it. The build starts from
``build_earth_dam_tseep.seep03_model()`` — the same section, the same two zones,
the same storage, the same boundary set — and adds only what a limit equilibrium
run reads: a starting circle on each face, a minimum slip depth that keeps the
search off surficial slivers, and a denser saved-frame schedule for the curve to
be drawn through.

Both files arrive SOLVED. The mesh and the whole transient march ship beside each
workbook as ``_mesh.json``, ``_tseep.csv`` and ``_tseep_meta.json``, so the reader
opens a dam whose nineteen pore-pressure fields already exist and goes straight to
the stability runs. Building a transient seepage model is Tutorial SEEP-3's
subject, and COMBO-3 does not repeat it. The march is the same for both files:
seepage reads conductivity and storage, and none of the five cells the reader
fills changes a pore pressure.

The two workbooks differ in ten cells. The start file leaves ``gamma``,
``gamma_sat``, ``c``, ``phi`` and ``u`` blank on both zones — the whole strength
band, unit weights included — which is what the page's materials step has the
reader type in; the completed file carries them, and the page's tagged results
are read off it. Only the strength model, ``option = mc``, ships filled on both
files: it is a property of the material, not a value the exercise measures.
Blank is the honest state for a value nobody has entered — the workbook still
loads, and the Model checks panel reports the gap: ``mat.gamma_nonpositive`` on
the missing unit weight, ``seep_field.no_consumer`` on the pore-pressure option
left at its default of ``none`` while a solved seepage field sits unread.

The strengths are typical values for a granular shell and a compacted clay core,
chosen for the exercise rather than measured on this dam; SEEP-3's file carries
none, because a seepage analysis reads none. The tutorial page states this in the
same words.

Run:  PYTHONPATH=. python3 tools/build_earth_dam_fs_time.py

The build solves the whole march once and writes both file sets from it, so it
takes a couple of minutes.
"""

from __future__ import annotations

import contextlib
import io
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from build_earth_dam_tseep import (REPO_ROOT, TUTORIAL_FILES, _write,  # noqa: E402
                                   seep03_model)

OUT = "xslope_earth_dam_fs_time.xlsx"
OUT_START = "xslope_earth_dam_fs_time_start.xlsx"

#: The five cells per zone the page has the reader type, in the order the
#: Materials editor lists them (``gamma`` and ``gamma_sat`` sit left of ``option``,
#: ``c`` and ``phi`` right of it, ``u`` further right still). They are left blank
#: on the start file and filled on the completed one; everything else on the two
#: mat rows — the name and the strength model, ``option = mc`` — is identical.
TYPED_IN = ("gamma", "gamma_sat", "c", "phi", "u")

#: The limit equilibrium band the page reads out, indexed to the base file's
#: material order. Drained effective-stress strengths on both zones, and
#: ``u = seep`` on both so every slice base reads its pore pressure from the solved
#: field rather than from a sketched line — which is the whole point of a curve
#: driven by a march. ``gamma`` is the moist unit weight above the water table and
#: ``gamma_sat`` the saturated one below it; the slicer splits each slice's weight
#: at the water table, so a drawdown that drains the shell lightens it.
STRENGTHS = [
    dict(option="mc", gamma=20.0, gamma_sat=21.0, c=0.0, phi=32.0, u="seep"),
    dict(option="mc", gamma=19.0, gamma_sat=20.0, c=10.0, phi=25.0, u="seep"),
]

#: Two starting circles, one per face, each placed on the deep mechanism its
#: slope can make: entering near the crest, tangent to the rock (``Depth`` is an
#: ELEVATION, and 0 is the rock the dam sits on), and daylighting near the toe.
#: The upstream face runs from the heel at x = 0 to the upstream crest edge at
#: x = 51, the downstream face from the downstream crest edge at x = 59 to the toe
#: at x = 110, so a circle spanning either face has its center out beyond the toe
#: and a radius near the span. Both faces are offered because both govern over a
#: run this long: a drawdown attacks the upstream slope, but a full reservoir
#: loads it, and under a full pool the weaker mechanism is the downstream one.
#: Which one wins at a given instant is left to the search.
#:
#: Both centers are one metre off the round value a hand would reach for, and that
#: is deliberate on the upstream side: a circle centered at (7, 55) with R 55
#: passes exactly through the crest-edge vertex at (51, 22), and a starting circle
#: that meets the ground surface at one of its own vertices has a degenerate trace.
CIRCLES = [dict(Xo=7.0, Yo=56.0, Depth=0.0, R=56.0),      # upstream face
           dict(Xo=103.0, Yo=59.0, Depth=0.0, R=59.0)]    # downstream face

#: The one search limit the file keeps. ``min_slip_depth`` puts a floor on the
#: SIZE of the mechanism. The shell is cohesionless, so nothing in it sets a
#: scale: the shallower a trial surface is drawn, the lower it scores, and an
#: unfiltered search walks down to a surficial wedge whose factor of safety is the
#: infinite-slope value and whose failure is raveling rather than a slide. 8 m
#: below the ground surface, on a 22 m dam, holds the search on slides of the
#: embankment rather than on crest slivers. No entry or exit range is declared:
#: which face governs is what the curve measures, so confining the trace to one of
#: them would answer the question in the input.
SEARCH_WINDOW = dict(min_slip_depth=8.0)

#: The saved-frame schedule this page sweeps, replacing SEEP-3's coarser one. A
#: curve is only as fine as the frames under it — a time that is not saved has no
#: field to read and is never interpolated — so the instants are packed where the
#: answer moves and spread where it does not. Every five days through the
#: drawdown, which ends on day 47; then widening steps through the recovery, when
#: the core is draining and each frame gains less than the one before. SEEP-3's
#: own file keeps its schedule: this is a stability question asked of the same
#: march, not a change to the seepage tutorial.
#:
#: ``save_interval`` is set to the full duration so the interval grid contributes
#: only the last frame and the list below states the schedule outright. The
#: stepper lands exactly on every saved time, and the widest gap between two of
#: them is 60 days, so the interval never caps a step.
DURATION = 300.0
SAVE_TIMES = [2.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 47.0,
              55.0, 65.0, 80.0, 100.0, 130.0, 180.0, 240.0, 300.0]

#: The mesh the file ships, stated the way Studio's Build Mesh dialog states it:
#: linear triangles, auto-sized from the geometry at 64 divisions across the 110 m
#: section. Head is a scalar field, so linear elements resolve it. This is SEEP-3's
#: own mesh, so the fields under COMBO-3's curve are the fields that tutorial
#: solves.
MESH_ELEMENT_TYPE = "tri3"
MESH_SIZE_DIVISIONS = 64


def _model(strengths=True):
    """SEEP-3's completed model with both starting surfaces, the search floor and
    the denser saved-frame schedule.

    ``strengths=False`` returns the same model with the five cells the reader
    types left unset. They are None rather than 0.0 on purpose: the writer puts a
    blank cell down for None, and a cohesionless material and an unfilled one are
    not the same model.
    """
    sd = seep03_model()
    for mat, upd in zip(sd["materials"], STRENGTHS):
        mat.update(upd)
        if not strengths:
            mat.update({k: None for k in TYPED_IN})
    sd["tseep"] = dict(sd["tseep"], duration=DURATION,
                       save_interval=DURATION, save_times=list(SAVE_TIMES))
    sd["circles"] = [dict(c) for c in CIRCLES]
    sd["circular"] = True
    sd["search_window"] = dict(SEARCH_WINDOW)
    return sd


def _quiet(fn, *a, **k):
    with contextlib.redirect_stdout(io.StringIO()):
        return fn(*a, **k)


def _companions(path, mesh, seep_data, solution, export_mesh_to_json,
                export_transient_solution):
    """The mesh and the march, written under one workbook's own base name.

    Companion discovery is by base name — ``{base}_mesh.json``, ``{base}_tseep.csv``
    and ``{base}_tseep_meta.json`` beside the ``.xlsx`` — so each workbook needs its
    own set, and each set names its own workbook in the meta file's provenance.
    """
    stem = os.path.splitext(path)[0]
    mesh_path = f"{stem}_mesh.json"
    export_mesh_to_json(mesh, mesh_path)
    csv_path, meta_path = _quiet(
        export_transient_solution, seep_data, solution, stem,
        input_file=os.path.basename(path),
        mesh_file=os.path.basename(mesh_path))
    for written in (mesh_path, csv_path, meta_path):
        print(f"  wrote {os.path.relpath(written, REPO_ROOT)} "
              f"({os.path.getsize(written) / 1024:.0f} KB)")


def build():
    """Both workbooks, and the mesh and march that ship beside each.

    The march is solved once. Seepage reads conductivity and storage, so the three
    strength cells that separate the two files cannot move a pore pressure, and
    solving twice would only produce the same nineteen frames again.
    """
    from xslope.fileio import load_slope_data
    from xslope.mesh import (build_mesh_from_polygons, export_mesh_to_json,
                             extract_size_regions, get_material_polygons)
    from xslope.seep import (build_seep_data, build_tseep_data,
                             export_transient_solution, run_transient_seepage)

    path = _write(_model(), TUTORIAL_FILES, OUT)
    start_path = _write(_model(strengths=False), TUTORIAL_FILES, OUT_START)

    # Re-read what was written, so the companions are built on the model as the
    # file states it rather than on the dict that produced it.
    sd = load_slope_data(path)
    xs = [x for x, _ in sd["ground_surface"].coords]
    target = (max(xs) - min(xs)) / MESH_SIZE_DIVISIONS
    mesh = _quiet(build_mesh_from_polygons, get_material_polygons(sd), target,
                  MESH_ELEMENT_TYPE, size_regions=extract_size_regions(sd))
    seep_data = _quiet(build_seep_data, mesh, sd, seep_bc=1)
    solution = _quiet(run_transient_seepage, seep_data, build_tseep_data(sd),
                      verbose=False)
    print(f"  mesh: {len(mesh['nodes'])} nodes, {len(mesh['elements'])} elements "
          f"(target size {target:g})")
    print(f"  march: {len(solution['frames'])} frames, "
          f"converged={solution['converged']}, "
          f"closure={solution['mass_balance']['final_closure']:.3e}")

    for p in (path, start_path):
        _companions(p, mesh, seep_data, solution, export_mesh_to_json,
                    export_transient_solution)
    return path


if __name__ == "__main__":
    build()
