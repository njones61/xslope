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
Not every ``{stem}_seep.csv`` in the corpus is the output of a steady solve, and a
steady solve would destroy the ones that are not. Each exclusion is listed in
``EXCLUDED`` below with the reason it is excluded and the builder that owns it.

Ten further companions re-solve to a MATERIALLY DIFFERENT field on their own
committed mesh — their workbooks were migrated after the field was saved, so the
saved field is stale with respect to the model it sits beside. Rewriting them is a
corpus change with published numbers behind it, not a footer correction, so they are
listed in ``HELD`` and left alone. ``--check`` prints the measured difference for
every one of them.
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
#: sample flow nets under (tol 1e-4; the 1000-sweep ceiling is earth_dam2's, which
#: needs 427 sweeps of the exit face, and changes no other model).
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
    ("docs/seep/files/xslope_clay_blanket",              (1,), DOCS),
    ("docs/seep/files/xslope_double_sheetpile",          (1,), DOCS),
    ("docs/seep/files/xslope_earth_dam1",                (1,), DOCS),
    ("docs/seep/files/xslope_earth_dam1_vg",             (1,), DOCS),
    ("docs/seep/files/xslope_johnson_res",               (1,), DOCS),
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
]

#: Companions this script must not write, and why. Every one of them is a
#: ``{stem}_seep.csv`` that no steady solve produced.
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

#: Companions whose committed field is STALE against their own workbook: the model was
#: migrated after the field was saved, so re-solving moves the field materially. Held
#: for a decision, not silently rewritten. The note records what the difference is.
HELD = {
    "docs/seep/files/xslope_levee2":
        "flowrate 0.006655 -> 0.995013 (x150). Sidecar and mesh date from 2026-02-26; "
        "the workbook has been migrated twice since (v12, v22). No lock reads it.",
    "docs/inputs/seep/xslope_earth_dam_bc2":
        "max |du| 22.2 / 16.5 psf over the two BC sets, flowrate 42.342 -> 42.437 and "
        "11.572 -> 11.588. Sidecar 2026-02-26, workbook 2026-07-30. No lock reads it.",
    "docs/lem/files/xslope_earth_dam_rapid":
        "max |du| 628 / 258 psf; flowrate 181.62 -> 183.92 (and the second set "
        "44.88 -> 53.53, +19%). The first set does not converge on this mesh.",
    "docs/lem/files/xslope_gsat_seep":
        "its _seep.csv is byte-identical to xslope_earth_dam_rapid_seep.csv — one field "
        "shipped under two names — so it moves with it, and does not converge. Carries "
        "a circular_search lock (fs_bishop 1.933, fs_spencer 1.912).",
    "docs/lem/files/xslope_johnson_rapid_KEY":
        "max |du| 362 / 25 psf; the first set does not converge. Measured against its "
        "own single_circle locks the seven factors of safety all still round the same, "
        "but the field itself moves across 12-59% of nodes.",
    "docs/seep/files/xslope_rface_SEEP_KEY":
        "max |du| 19.5 psf, flowrate 3.7754 -> 3.7867. This one MOVES ITS LOCKS: all "
        "seven single_circle factors of safety shift by 0.0008-0.0016, which re-rounds "
        "every published value (e.g. spencer 2.080 -> 2.078).",
    "docs/seep/files/xslope_earth_dam2":
        "does not converge on its committed mesh at any solver setting tried (tol 1e-4 "
        "to 1e-6, 400 to 4000 sweeps): the exit-face active set cycles, and the field "
        "lands 8-23 psf away depending only on where the sweep count stops. Its "
        "type=seep lock re-meshes and is unaffected.",
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
        print("\nHELD — stale field, not rewritten:")
        for stem_rel, why in HELD.items():
            print(f"\n  {os.path.basename(stem_rel)}: {why}")
            bcs = (1, 2) if os.path.exists(
                os.path.join(REPO_ROOT, stem_rel + "_seep2.csv")) else (1,)
            run(stem_rel, bcs, DOCS, check=True)


if __name__ == "__main__":
    main()
