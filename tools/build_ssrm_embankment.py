"""Build Tutorial FEM-1's file PAIR (the strength-reduction embankment).

The model is Griffiths & Lane (1999) Example 1 — the homogeneous 2:1 embankment
at c/gamma-H = 0.05 with phi = 20 degrees that is the base SSRM benchmark.  That
is the verification case at
[griffiths1](../docs/verification/ssrm.md#verification-griffiths1), and this
builder DERIVES both tutorial files from the committed corpus file rather than
restating its numbers:

    docs/fem/files/xslope_griffiths1.xlsx                (the base)
      -> docs/tutorials/files/xslope_ssrm_embankment_start.xlsx
      -> docs/tutorials/files/xslope_ssrm_embankment.xlsx

so the section, the soil and the elastic constants a reader types can never drift
from the ones the verification page measures.  Nothing here hand-edits an xlsx.

What the tutorial pair changes from the base, and why
----------------------------------------------------
* **The starter carries no elastic constants.**  E and nu are cleared to blank
  cells — not zeroed, which would be a stiffness of zero rather than an unfilled
  one.  The reader's arc is limit equilibrium first (which needs neither), then
  the two cells the finite element engine does need, then a mesh, then the
  strength reduction.  Everything else on the starter is the completed file's,
  so the only thing the reader adds before meshing is the pair.
* **Neither file carries a mesh sidecar.**  The base ships one
  (``xslope_griffiths1_mesh.json``, the quad8 benchmark mesh); the tutorial
  builds its own from the Build Mesh dialog, and a file that already had one
  would skip the step the page teaches.  Both tutorial files are therefore the
  workbook alone, as the SEEP-4 pair is.
* **The completed file declares the mesh it was solved on** — tri6 at a target
  element size of 3.5 — on the main sheet, so Studio's Build Mesh dialog opens on
  the tutorial's own mesh instead of on the auto size.  The base leaves both
  blank.  The strength-reduction bracket is deliberately left blank in both
  files, so the Run FEM dialog opens on its own defaults (F in [1.0, 2.0]),
  which is the bracket the page walks.

The geometry, the material (gamma = 125 pcf, c = 312.5 psf, phi = 20 degrees,
option ``mc``), the starting circle, the English units and — in the completed
file — the elastic pair E = 2,088,500 psf and nu = 0.3 are the base file's,
verbatim.

Run:  PYTHONPATH=. python3 tools/build_ssrm_embankment.py           # both
      PYTHONPATH=. python3 tools/build_ssrm_embankment.py start     # one
"""

from __future__ import annotations

import os
import sys

from xslope.fileio import (load_slope_data, save_slope_data_to_xlsx,
                           default_template_path)

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BASE = os.path.join(REPO_ROOT, "docs", "fem", "files", "xslope_griffiths1.xlsx")
TUTORIAL_FILES = os.path.join(REPO_ROOT, "docs", "tutorials", "files")

START_OUT = "xslope_ssrm_embankment_start.xlsx"
DONE_OUT = "xslope_ssrm_embankment.xlsx"

#: The mesh the tutorial solves on, declared on the completed file's main sheet.
#: tri6 because a factor of safety needs quadratic elements (docs/fem/overview.md,
#: "Element type and volumetric locking"), and 3.5 because that is the target size
#: the verification page's own element-type tag builds every quadratic type at.
ELEMENT_TYPE = "tri6"
TARGET_SIZE = 3.5


def _model():
    """The base model, loaded once, with its committed mesh dropped.

    ``load_slope_data`` picks up ``xslope_griffiths1_mesh.json`` from beside the
    base file; neither tutorial file has a sidecar, so the in-memory mesh goes
    too and no run in this builder can accidentally inherit the benchmark's
    quad8 discretization.
    """
    sd = load_slope_data(BASE)
    sd.pop("mesh", None)
    return sd


def _write(sd, filename):
    out = os.path.join(TUTORIAL_FILES, filename)
    save_slope_data_to_xlsx(sd, out, template=default_template_path())
    print(f"wrote {os.path.relpath(out, REPO_ROOT)}")
    return out


def build_start():
    """The starter: the embankment with no elastic constants.

    E and nu are set to ``None`` rather than to a number, which the material
    writer puts out as an empty cell — an unfilled property, not a zero one.  The
    file still carries a failure surface, so it loads on the ordinary analysis
    path and the limit-equilibrium search the tutorial opens with runs on it as
    delivered.
    """
    sd = _model()
    for m in sd["materials"]:
        m["E"] = None
        m["nu"] = None
    return _write(sd, START_OUT)


def build_done():
    """The completed file: the base's elastic pair, plus the mesh it is solved on."""
    sd = _model()
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
