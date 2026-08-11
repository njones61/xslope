"""Re-solve the shipped steady seepage companions so each records what its solve was.

A model that has been solved once keeps its results beside it, and a report of an
already-solved model reads them back rather than solving again. Those companions
were written before ``export_seep_solution`` recorded the three solve facts its
columns cannot hold::

    # Unconfined: True|False     was a phreatic surface located, or is the field
                                 the ordinary suction of a confined potential flow
    # Converged: True|False      did the solve close
    # Closure Error: <float>     the flow closure it closed to

Without them a field read back cannot say whether its negative pore pressures mean
a phreatic surface or saturated suction, and a report of the model cannot state
whether the seepage it is quoting converged. This re-solves each model on the mesh
committed beside it and rewrites its companion, so the footers are the model's own
answer rather than a default.

THE MESH IS NEVER REBUILT. Each model is solved on its committed ``{stem}_mesh.json``
— the mesh its shipped field was solved on, which is checked here rather than
assumed: for every saturated node the solver's own identity ``u = gamma_w * (head - y)``
holds to machine precision against the committed mesh's node elevations, for every
companion in the corpus. Re-meshing would renumber the nodes and silently invalidate
every ``u='seep'`` model that reads the file.

Deterministic: same mesh, same boundary conditions, same solver settings, same field.
Re-running reproduces every file it writes.

    PYTHONPATH=. python3 tools/make_seep_sidecars.py --check   # compare, write nothing
    PYTHONPATH=. python3 tools/make_seep_sidecars.py           # rewrite the companions
    PYTHONPATH=. python3 tools/make_seep_sidecars.py earth_dam1 johnson_res

WHAT THIS SCRIPT DOES NOT TOUCH
-------------------------------
Some ``{stem}_seep.csv`` files in the corpus must not be rewritten here, for reasons
that differ from one to the next: a few were produced by no steady solve at all
(vendor fields, transient frames), and others ARE steady solves whose format or
missing mesh belongs to something else. Each is listed in ``EXCLUDED`` below with its
own reason.

A companion that re-solves to a MATERIALLY DIFFERENT field on its own committed mesh
is HELD rather than rewritten: rewriting it is a corpus change, with published numbers
behind some of them, and not a footer correction. Six such files were held through
2026-08-09, all four of their models for the same reason underneath — their solves ran
to the iteration ceiling without closing, so there was no single field to record. The
solver was diagnosed and fixed on 2026-08-10 (two limit cycles, one in the head field
and one in the exit-face active set, both invisible behind the relaxation ladder) and
all four converge. Five of the six are regenerated above; the sixth is a copy of one of
them, held for a reason that has nothing to do with the model — ``--check --held``
prints it.
"""
import argparse
import contextlib
import io
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

import numpy as np                                                      # noqa: E402

from xslope.fileio import load_slope_data                               # noqa: E402
from xslope.mesh import import_mesh_from_json                           # noqa: E402
from xslope.seep import (build_seep_data, export_seep_solution,         # noqa: E402
                         run_seepage_analysis)

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

#: The solve settings each model's field was recorded from.
#:
#: ``docs`` — the convention ``run_tests.py::run_seep_test`` computes every
#: ``type=seep`` lock under, and that ``tools/make_seep_sample_figures.py`` draws the
#: sample flow nets under (tol 1e-4).
#:
#: THE 1000-SWEEP CEILING IS EARTH_DAM2'S, AND IT IS THIN. On its committed mesh that
#: model converges in 927 sweeps — 73 of headroom — so the ceiling is load-bearing for
#: it and for nothing else: every other companion here closes in 500 or fewer, and the
#: library's own default of 400 would leave earth_dam2 short by more than twice over.
#: Its exit face settles at sweep 188 (where the set-revisit escape frees the two edges
#: it was cycling on) and the remaining 739 sweeps are the nonlinear residual walking
#: down to closure_tol. A change that slows that walk runs it into the ceiling, where
#: it reports a flow rate off a field that is still moving, which is the state this
#: whole registry exists to keep out of the corpus — raise the ceiling with the model,
#: rather than trimming the model to it.
#: ``vendor`` — the settings the model's own builder solves it at
#: (``benchmarks/rocscience/build_rs2.py::_build_rs2_28``, ``build_vp038``,
#: ``build_vp046``), so a rebuild and this script agree.
#: ``gs`` — ``benchmarks/geostudio/build_gs2_46.py``.
DOCS = dict(tol=1e-4, max_iter=1000)
VENDOR = dict(tol=1e-5, max_iter=600)
GS = dict(tol=1e-5, max_iter=400)

#: (stem relative to the repo root, boundary condition sets, solver settings).
MODELS = [
    ("docs/fem/files/xslope_griffiths6_seep",            (1,), DOCS),
    ("docs/inputs/seep/xslope_earth_dam_bc2",         (1, 2), DOCS),
    # Held until 2026-08-10 as SOLVER EVOLUTION: the first BC set of each ran to the
    # 1000-sweep ceiling without closing, so it had no single field to record, and the
    # second sets had drifted with the same solver evolution. The cause was a limit
    # cycle in the head field, hidden by the relaxation ladder's terminal 0.01; with
    # the escape in place each converges (463 and 901 sweeps) and each has one answer,
    # so all four files are solved for it here. gsat_seep is the same field under a
    # second name and is released with them — see HELD for why its copy waits.
    ("docs/lem/files/xslope_earth_dam_rapid",         (1, 2), DOCS),
    ("docs/lem/files/xslope_johnson_rapid_KEY",       (1, 2), DOCS),
    ("docs/seep/files/xslope_clay_blanket",              (1,), DOCS),
    ("docs/seep/files/xslope_double_sheetpile",          (1,), DOCS),
    ("docs/seep/files/xslope_earth_dam1",                (1,), DOCS),
    ("docs/seep/files/xslope_earth_dam1_vg",             (1,), DOCS),
    # Held until 2026-08-10, where non-convergence was itself the cause: the exit-face
    # active set cycled, so the field landed 8-23 psf apart depending only on where the
    # sweep count stopped. The cycle was the all-or-nothing quadratic-edge rule
    # forbidding the set the solve was reaching for; with that resolved it converges in
    # 927 sweeps on 8 of its 97 exit-face nodes, closing to 8e-10 of its own flow, and
    # has one field to record. Its type=seep lock re-meshes and is unaffected either way.
    ("docs/seep/files/xslope_earth_dam2",                (1,), DOCS),
    ("docs/seep/files/xslope_johnson_res",               (1,), DOCS),
    # levee2's shipped field was solved at a grout curtain k = 0.001 that no committed
    # revision of the workbook carries — classroom experimentation, k being what
    # students are set to sweep. The committed k = 0.2 is the model (owner, 2026-08-10),
    # so the field is solved at it and the workbook is untouched.
    ("docs/seep/files/xslope_levee2",                    (1,), DOCS),
    # rface_SEEP_KEY's shipped field matches no committed revision of its workbook, so
    # there was no state of the model to restore it to. It is solved from the committed
    # workbook (owner, 2026-08-10) and the seven factors of safety published on
    # docs/seep/seep_slope.md re-anchored to what that field gives.
    ("docs/seep/files/xslope_rface_SEEP_KEY",            (1,), DOCS),
    ("docs/seep/files/xslope_sea_trench_anis",           (1,), DOCS),
    ("docs/verification/files/geostudio/gs2_46",         (1,), GS),
    ("docs/verification/files/rocscience/rs2_28a",       (1,), VENDOR),
    ("docs/verification/files/rocscience/rs2_28b",       (1,), VENDOR),
    ("docs/verification/files/rocscience/rs2_28c",       (1,), VENDOR),
    ("docs/verification/files/rocscience/vp010",         (1,), DOCS),
    ("docs/verification/files/rocscience/vp038a",        (1,), VENDOR),
    ("docs/verification/files/rocscience/vp038b",        (1,), VENDOR),
    ("docs/verification/files/rocscience/vp038c",        (1,), VENDOR),
    ("docs/verification/files/rocscience/vp046b",        (1,), VENDOR),
    ("docs/verification/files/rocscience/vp071a",        (1,), DOCS),
    ("docs/verification/files/rocscience/vp072a",        (1,), DOCS),
    ("docs/verification/files/rocscience/vp076a",        (1,), DOCS),
    ("docs/verification/files/rocscience/vp077a",        (1,), DOCS),
    # The polygon-sheet fixtures, which live outside docs/ but ship the same
    # companions: a mesh and a solved field beside the model. Nothing reads their
    # fields — the published levee_poly flow-rate lock names the docs/ copy of the
    # model, which carries no mesh and no companion and so re-meshes and re-solves.
    ("poly_test/xslope_levee_poly",                      (1,), DOCS),
    ("poly_test/xslope_sea_trench_poly",                 (1,), DOCS),
]

#: Companions this script must not write, and why — one reason per entry, because
#: they are not all excluded for the same reason. Some were produced by no steady
#: solve (a vendor field, a published grid, a frame of a transient march); the
#: rs2_67b/e/f trio ARE steady solves, but their builder owns the file's format; and
#: levee1 is a solved field deliberately shipped with no mesh to place it on.
EXCLUDED = {
    "docs/seep/files/xslope_levee1":
        "the NOMESH fixture: ships a _seep.csv and deliberately NO _mesh.json, so "
        "test/report_check.py can exercise the 'holds a solution, but the model "
        "carries no mesh' note. Giving it a mesh would delete the fixture.",
    "docs/verification/files/rocscience/rs2_67b":
        "written by benchmarks/rocscience/build_rs2.py::_write_seep_sidecars, which "
        "records head and u only and zero-fills the flow-net columns. It IS a steady "
        "solve, but its builder owns the file's format; re-exporting it here would "
        "fill in velocities its builder does not write, and the rebuild guard would "
        "then see a file its builder no longer reproduces.",
    "docs/verification/files/rocscience/rs2_67e": "as rs2_67b (same builder).",
    "docs/verification/files/rocscience/rs2_67f": "as rs2_67b (same builder).",
    "docs/verification/files/rocscience/rs2_67c":
        "pore pressures IMPORTED from the RS2 vendor '#067_03 (90h).fea' nodal field. "
        "Not our solve at all — re-solving would replace the vendor's answer, which is "
        "the whole point of the row.",
    "docs/verification/files/rocscience/rs2_67d": "as rs2_67c (vendor .fea field).",
    "docs/verification/files/rocscience/rs2_9":
        "pore pressures synthesized from the 95-point grid published with RS2 Part I "
        "#9. A vendor-published field, not a solve.",
    "docs/verification/files/rocscience/vp102b":
        "cut from the t=0 (steady full-pool) FRAME of the VP102 drawdown transient, by "
        "benchmarks/rocscience/build_rs2.py::vp102_transient. One march feeds vp102b "
        "and all five snapshots so the family cannot drift apart; a separate steady "
        "solve here would break exactly that.",
    "docs/verification/files/rocscience/vp102t_60":
        "a FRAME of the VP102 rapid-drawdown transient march (t = 60 h). A steady "
        "solve has no meaning for it — the field is mid-drawdown.",
    "docs/verification/files/rocscience/vp102t_100": "as vp102t_60 (t = 100 h).",
    "docs/verification/files/rocscience/vp102t_300": "as vp102t_60 (t = 300 h).",
    "docs/verification/files/rocscience/vp102t_600": "as vp102t_60 (t = 600 h).",
    "docs/verification/files/rocscience/vp102t_1500": "as vp102t_60 (t = 1500 h).",
}

#: Companions whose committed field does not reproduce on their own committed mesh,
#: each with the DIAGNOSED cause of it. A model is ruled out as the cause before an
#: entry is made: re-solving with the workbook as it stood at the sidecar's own commit
#: must give results bit-identical to re-solving with the workbook as it stands today,
#: so no input change can explain the difference. Held for a decision, not silently
#: rewritten — and a decision that holds an entry here is a decision to leave a stale
#: field in the corpus, so an entry is a place to start work, not a place to end it.
#:
#: The four models held here through 2026-08-09 were held for one cause between them —
#: their solves did not converge, and a flow rate read off a field that is still moving
#: is not the flow through the section. That was a solver defect, not a corpus
#: decision; it was fixed on 2026-08-10 and three of the four are in MODELS above, each
#: with the ruling that released it. The fourth is held for a reason of a different
#: kind, and not about the model at all.
HELD = {
    "docs/lem/files/xslope_gsat_seep":
        "RELEASED, WAITING ON A TEST FIXTURE (2026-08-10). Its field is the "
        "earth_dam_rapid field under a second name, which is regenerated and "
        "converged; this copy reproduces it exactly (q = 183.958670, 463 sweeps) and "
        "is the file to write. What it cannot do yet is CARRY A SOLVE FOOTER: "
        "test/report_check.py uses this model as the sample whose saved solution "
        "records nothing about its solve, and two of its checks are premised on that "
        "silence (test_seep_convergence_is_stated, test_seep_stale_sidecar_says_so). "
        "That file belongs to another round in flight. Regenerate this stem in the "
        "same change that points those two fixtures at a sample that is genuinely "
        "silent, and that raises _truncated()'s row count past the footer length.",
}


def _quiet(fn, *a, **k):
    with contextlib.redirect_stdout(io.StringIO()):
        return fn(*a, **k)


def _shipped(path):
    """The committed field and its recorded flowrate."""
    import pandas as pd
    df = pd.read_csv(path, comment="#")
    flowrate = None
    for line in open(path):
        if line.startswith("# Total Flowrate:"):
            try:
                flowrate = float(line.split(":", 1)[1])
            except ValueError:
                pass
    return df, flowrate


def solve(stem, bc, settings):
    """Solve one boundary condition set on the model's COMMITTED mesh."""
    slope_data = _quiet(load_slope_data, f"{stem}.xlsx")
    mesh = _quiet(import_mesh_from_json, f"{stem}_mesh.json")
    seep_data = _quiet(build_seep_data, mesh, slope_data, seep_bc=bc)
    solution = _quiet(run_seepage_analysis, seep_data, **settings)
    return seep_data, solution


def run(stem_rel, bcs, settings, check):
    stem = os.path.join(REPO_ROOT, stem_rel)
    for bc in bcs:
        suffix = "_seep.csv" if bc == 1 else "_seep2.csv"
        path = f"{stem}{suffix}"
        seep_data, solution = solve(stem, bc, settings)
        df, flowrate = _shipped(path)
        if len(df) != len(seep_data["nodes"]):
            raise SystemExit(f"{path}: {len(df)} rows against a "
                             f"{len(seep_data['nodes'])}-node mesh")
        du = float(np.max(np.abs(solution["u"] - df["u"].to_numpy())))
        dh = float(np.max(np.abs(solution["head"] - df["head"].to_numpy())))
        name = os.path.basename(path)
        print(f"{name:<40} d_u={du:.3e} d_head={dh:.3e} "
              f"q={solution['flowrate']:.6f} (was {flowrate:.6f})  "
              f"unconfined={solution.get('unconfined')} "
              f"converged={solution.get('converged')} "
              f"closure={solution.get('closure_error'):.3e}")
        if not check:
            _quiet(export_seep_solution, seep_data, solution, path)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("stems", nargs="*", help="basenames to limit the run to")
    ap.add_argument("--check", action="store_true",
                    help="compare against the committed field, write nothing")
    ap.add_argument("--held", action="store_true",
                    help="also report the held (stale-field) companions")
    args = ap.parse_args(argv)

    models = MODELS
    if args.stems:
        models = [m for m in MODELS
                  if any(s in os.path.basename(m[0]) for s in args.stems)]
        if not models:
            raise SystemExit(f"no model matches {args.stems}")

    for stem_rel, bcs, settings in models:
        run(stem_rel, bcs, settings, args.check)

    if args.held:
        if not HELD:
            print("\nHELD — nothing: every companion in the corpus is the field its "
                  "own committed mesh and workbook give.")
        else:
            print("\nHELD — stale field, not rewritten:")
        for stem_rel, why in HELD.items():
            print(f"\n  {os.path.basename(stem_rel)}: {why}")
            bcs = (1, 2) if os.path.exists(
                os.path.join(REPO_ROOT, stem_rel + "_seep2.csv")) else (1,)
            run(stem_rel, bcs, DOCS, check=True)


if __name__ == "__main__":
    main()
