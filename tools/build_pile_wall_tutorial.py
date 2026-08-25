"""Build Tutorial FEM-3's wall PAIR: a continuous sheet pile wall on the pile
tutorial's own slope.

FEM-3's first half puts two rows of drilled shafts through both engines.  Its
second half is the case the finite element engine is right for — a member that
really is continuous out of plane, whose section constants are already per foot
of wall and whose stiffness is therefore not smeared over a spacing.  Keeping
that member on the SAME slope as the pile rows is what makes the two halves
readable against each other: identical geometry, identical clay, identical mesh
settings, one structural line instead of two.

    docs/lem/files/xslope_piles.xlsx            the pile tutorial's model
      -> docs/tutorials/files/xslope_pile_wall_start.xlsx   both rows removed
      -> docs/tutorials/files/xslope_pile_wall.xlsx         the wall row added

The starter is the bare slope; the completed file is the same slope with the one
row the page's reader enters.  Both are written at the current template version
and are sidecar-free — no ``_mesh.json`` and no ``_seep.csv`` — because the page
builds the mesh in Studio and the model carries no seepage.

The wall
--------
A PZ-27 steel sheet pile section, driven from the crest of the face at (10, 10)
down to the rigid base at (10, -10), stated per foot of wall:

    E     4.176e9 psf     29,000 ksi
    I     0.00888 ft^4/ft   184.20 in^4/ft
    Area  0.0551 ft^2/ft      7.94 in^2/ft
    Mcap  90,600 lb*ft/ft   Fy 36 ksi x S 30.2 in^3/ft, the moment at first yield
    S     1

``S`` = 1 passes those constants into the beam elements unchanged: a wall has no
out-of-plane gap, so there is nothing to smear.  ``D`` is blank because there is
no circular shaft to derive a section from and no gap for an arching theory to
act in, ``Vcap`` is blank because the section's shear rating is not part of what
the page checks, and ``H`` is blank because a strength reduction run never reads
a stated pile force.  The row ships with its tip free; the page has the reader
set it to fixed and re-run.

Run:  PYTHONPATH=. python3 tools/build_pile_wall_tutorial.py         # both
      PYTHONPATH=. python3 tools/build_pile_wall_tutorial.py start   # one
"""

from __future__ import annotations

import os
import sys
import warnings

warnings.filterwarnings("ignore")

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO_ROOT)

from xslope.fileio import (default_template_path, load_slope_data,  # noqa: E402
                           save_slope_data_to_xlsx)

PILES_BASE = os.path.join(REPO_ROOT, "docs", "lem", "files", "xslope_piles.xlsx")
TUTORIAL_FILES = os.path.join(REPO_ROOT, "docs", "tutorials", "files")

START_OUT = "xslope_pile_wall_start.xlsx"
DONE_OUT = "xslope_pile_wall.xlsx"

#: The mesh the page's numbers are measured on, declared in the files so both
#: models open their Build Mesh dialog on it: quadratic triangles (the only
#: element the strength reduction path runs on) at 2 ft, the pile model's own
#: target size.
ELEMENT_TYPE = "tri6"
TARGET_SIZE = 2.0

#: The strength reduction bracket the page runs every trial over.  1.6 would sit
#: below the socketed wall's answer, so the files declare 2.0 and the dialog
#: opens there.
SSRM_F_MIN, SSRM_F_MAX = 1.0, 2.0

#: The PZ-27 sheet pile row, per foot of wall.  See the module docstring for
#: where each number comes from.
WALL_ROW = {
    "x1": 10.0, "y1": 10.0, "x2": 10.0, "y2": -10.0,
    "H": None, "theta_p": 0.0, "D_pile": None, "S": 1.0,
    "E": 4.176e9, "I": 0.00888, "area": 0.0551,
    "V_cap": None, "M_cap": 90600.0,
    "head_fixity": "free", "tip_fixity": "pinned",
    "appl": "active", "label": "sheet pile wall",
}


def _base():
    """The pile tutorial's slope with both pile rows removed and the mesh and
    strength reduction settings the page runs at declared on the main sheet."""
    sd = load_slope_data(PILES_BASE)
    if len(sd["pile_lines"]) != 2:
        raise SystemExit(f"{os.path.basename(PILES_BASE)} carries "
                         f"{len(sd['pile_lines'])} pile rows, expected 2")
    sd.pop("mesh", None)
    sd.pop("seep_u", None)
    sd["pile_lines"] = []
    sd["element_type"] = ELEMENT_TYPE
    sd["target_size"] = TARGET_SIZE
    sd["ssrm_f_min"] = SSRM_F_MIN
    sd["ssrm_f_max"] = SSRM_F_MAX
    return sd


def _write(sd, filename):
    out = os.path.join(TUTORIAL_FILES, filename)
    save_slope_data_to_xlsx(sd, out, template=default_template_path())
    for sidecar in (f"{os.path.splitext(out)[0]}_mesh.json",
                    f"{os.path.splitext(out)[0]}_seep.csv"):
        if os.path.exists(sidecar):
            raise SystemExit(f"tutorial downloads are sidecar-free, but "
                             f"{os.path.basename(sidecar)} sits beside {filename}")
    print(f"wrote {os.path.relpath(out, REPO_ROOT)}")
    return out


def build_start():
    """The starter: the slope on its own, with no structural line in it.

    This is the page's first run of the second half — mesh it, run the strength
    reduction, and read what the slope holds without the wall.
    """
    return _write(_base(), START_OUT)


def build_done():
    """The completed file: the same slope with the one wall row the reader
    enters, tip free."""
    sd = _base()
    sd["pile_lines"] = [dict(WALL_ROW)]
    return _write(sd, DONE_OUT)


TARGETS = {"start": build_start, "done": build_done}


def main(argv):
    for name in argv[1:] or list(TARGETS):
        if name not in TARGETS:
            raise SystemExit(f"unknown target {name!r}; choices: {list(TARGETS)}")
        TARGETS[name]()


if __name__ == "__main__":
    main(sys.argv)
