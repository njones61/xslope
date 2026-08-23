"""Build Tutorial FEM-2's file PAIR (the reinforced slope, LEM against FEM).

The model is the six-layer geogrid slope that Tutorial LEM-8 builds and that
[LEM Sample Problem 9](../docs/lem/samples.md#9-reinforced-slope) locks Spencer's
answer for.  Its finite element counterpart is
[FEM Sample Problem 1](../docs/fem/samples.md).  This builder DERIVES both
tutorial files from those two committed models rather than restating a number
from either:

    docs/lem/files/xslope_reinforce.xlsx        (the section, soil, surcharge,
                                                 circles and the six geogrid
                                                 lines with their pullout data)
    docs/fem/files/xslope_reinforce_fem.xlsx    (the elastic constants, and the
                                                 three FEM-only reinforcement
                                                 columns Tres / E / Area)
      -> docs/tutorials/files/xslope_reinforced_slope_start.xlsx
      -> docs/tutorials/files/xslope_reinforced_slope.xlsx

so the model a reader opens can never drift from the one the two sample pages
measure.  Nothing here hand-edits an xlsx, and no input is the tutorial's own:
every value comes from one of the two base models.

The starter is the LEM-8 model plus the elastic constants
--------------------------------------------------------
The page's arc is: recall the limit-equilibrium answer, then give the same
section to the finite element engine and watch what the reinforcement does that
the slices could not express.  So the starter carries **everything limit
equilibrium needs and everything the continuum needs**, and is missing exactly
the two things the page teaches a reader to supply:

* **no mesh** — neither file has a mesh sidecar, and the starter declares no
  element type or target size, so the reader builds the mesh from the Build Mesh
  dialog.  (The completed file declares tri6 at the tutorial's target size, which
  is what makes that dialog open on the tutorial's own mesh.)
* **no FEM reinforcement data** — Tres, E and Area are blank on the starter's six
  reinforcement rows.  A blank Tres is not zero and not a defect: the loader
  reads it as NaN, meaning *no post-peak drop*, and the finite element engine
  softens only where ``t_res`` is finite (``fileio.py``, the ``t_res`` comment on
  the reinforcement loader).  So the starter's own state is a legitimate
  elastic-perfectly-plastic model, which is the page's first run; the reader
  types the residual capacity in afterwards for the second.

Everything else on the starter is the completed file's, so the only inputs the
reader adds are the mesh and the three reinforcement columns.

Run:  PYTHONPATH=. python3 tools/build_reinforced_slope_tutorial.py         # both
      PYTHONPATH=. python3 tools/build_reinforced_slope_tutorial.py start   # one
"""

from __future__ import annotations

import math
import os
import sys

from xslope.fileio import (load_slope_data, save_slope_data_to_xlsx,
                           default_template_path)

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LEM_BASE = os.path.join(REPO_ROOT, "docs", "lem", "files", "xslope_reinforce.xlsx")
FEM_BASE = os.path.join(REPO_ROOT, "docs", "fem", "files", "xslope_reinforce_fem.xlsx")
TUTORIAL_FILES = os.path.join(REPO_ROOT, "docs", "tutorials", "files")

START_OUT = "xslope_reinforced_slope_start.xlsx"
DONE_OUT = "xslope_reinforced_slope.xlsx"

#: The mesh the tutorial solves on, declared on the completed file's main sheet.
#: tri6 because a factor of safety needs quadratic elements (docs/fem/overview.md,
#: "Element type and volumetric locking"), and 2 ft because that is the target size
#: the FEM sample page's own test tag builds this model at.
ELEMENT_TYPE = "tri6"
TARGET_SIZE = 2.0

#: The reinforcement columns the finite element engine reads and the limit
#: equilibrium engine does not: the residual (post-peak) capacity, the bar's
#: axial modulus and its area per unit width.  Blank on the starter, taken from
#: the FEM sample model on the completed file.
FEM_REINFORCE_KEYS = ("t_res", "E", "area")

#: The geometry and capacity terms that must agree between the two base models
#: before either can be called "the same slope".  Checked, not assumed.
SHARED_REINFORCE_KEYS = ("x1", "y1", "x2", "y2", "t_max", "lp1", "lp2",
                         "adhesion", "delta",
                         "tend1", "tend2", "spacing")


def _same(a, b, tol=1e-9):
    if isinstance(a, float) and isinstance(b, float):
        if a != a and b != b:          # both unset (NaN): the same blank
            return True
        return math.isclose(a, b, rel_tol=tol, abs_tol=tol)
    return a == b


def _bases():
    """The two committed models, with the FEM sample's mesh dropped.

    ``load_slope_data`` picks up ``xslope_reinforce_fem_mesh.json`` from beside
    the FEM sample; neither tutorial file has a sidecar, so the in-memory mesh
    goes too and no run in this builder can inherit it.
    """
    lem = load_slope_data(LEM_BASE)
    fem = load_slope_data(FEM_BASE)
    fem.pop("mesh", None)

    # The two files must describe one slope. Anything else and the tutorial's
    # "the same model, now in the finite element engine" claim is not true.
    if str(lem["ground_surface"]) != str(fem["ground_surface"]):
        raise SystemExit("the two base models do not share a ground surface")
    lr, fr = lem["reinforcement_lines"], fem["reinforcement_lines"]
    if len(lr) != len(fr):
        raise SystemExit("the two base models carry different reinforcement counts")
    for i, (a, b) in enumerate(zip(lr, fr), start=1):
        for k in SHARED_REINFORCE_KEYS:
            if not _same(a[k], b[k]):
                raise SystemExit(f"reinforcement line {i}: {k} differs between "
                                 f"the base models ({a[k]!r} vs {b[k]!r})")
    for a, b in zip(lem["materials"], fem["materials"]):
        for k in ("name", "gamma", "option", "c", "phi"):
            if not _same(a[k], b[k]):
                raise SystemExit(f"material {a['name']}: {k} differs between the "
                                 f"base models ({a[k]!r} vs {b[k]!r})")
    return lem, fem


def _model():
    """The tutorial model: the LEM-8 slope carrying the FEM sample's elastic
    constants, with the FEM-only reinforcement columns still blank."""
    lem, fem = _bases()
    for m, fm in zip(lem["materials"], fem["materials"]):
        m["E"], m["nu"] = fm["E"], fm["nu"]
    for r in lem["reinforcement_lines"]:
        r["t_res"] = float("nan")
        r["E"] = float("nan")
        r["area"] = float("nan")
    return lem, fem


def _write(sd, filename):
    out = os.path.join(TUTORIAL_FILES, filename)
    save_slope_data_to_xlsx(sd, out, template=default_template_path())
    print(f"wrote {os.path.relpath(out, REPO_ROOT)}")
    return out


def build_start():
    """The starter: the reinforced slope with elastic constants and no mesh, and
    with Tres / E / Area left blank on every reinforcement row.

    Blank rather than zero on all three.  ``save_slope_data_to_xlsx`` writes an
    unset Tres as an empty cell, and E and Area go out as the NaN they came in
    as, which the writer's number formatter also puts out empty — so the reader
    opens three genuinely empty columns and the loader reads them back as "not
    given" rather than as "zero", which for Tres would mean brittle rupture.
    """
    sd, _fem = _model()
    return _write(sd, START_OUT)


def build_done():
    """The completed file: the same model with the three reinforcement columns
    filled from the FEM sample, and the mesh it is solved on declared."""
    sd, fem = _model()
    for r, fr in zip(sd["reinforcement_lines"], fem["reinforcement_lines"]):
        for k in FEM_REINFORCE_KEYS:
            r[k] = fr[k]
    sd["element_type"] = ELEMENT_TYPE
    sd["target_size"] = TARGET_SIZE
    return _write(sd, DONE_OUT)


TARGETS = {"start": build_start, "done": build_done}


def main(argv):
    for name in argv[1:] or list(TARGETS):
        if name not in TARGETS:
            raise SystemExit(f"unknown target {name!r}; "
                             f"choices: {list(TARGETS)}")
        TARGETS[name]()


if __name__ == "__main__":
    main(sys.argv)
