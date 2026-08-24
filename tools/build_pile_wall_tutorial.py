"""Build Tutorial FEM-3's wall PAIR (the continuous sheet pile wall, in the FEM).

FEM-3 compares what the two engines do with a row of piles.  Its second half is
the case the finite element engine is right for: a member that really is
continuous out of plane, where a plane-strain beam's EA and EI are already per
metre and nothing has to be smeared.  The model is the SIGMA/W example
transcribed for the verification corpus:

    docs/verification/files/geostudio/gs2_wall.xlsx   the slope + the wall
    docs/verification/files/geostudio/gs2_wall_none.xlsx   the slope alone
      -> docs/tutorials/files/xslope_pile_wall_start.xlsx   (no wall row)
      -> docs/tutorials/files/xslope_pile_wall.xlsx         (with the wall row)

and this builder DERIVES the tutorial pair from the corpus pair rather than
restating the model.  ``benchmarks/geostudio/build_gs2_wall.py`` stays the
authority on what the model is; this file only re-expresses it in the form the
campaign's Download rule requires.  Every value a reader sees came out of the
committed corpus workbook.

Why a tutorial copy at all
--------------------------
Two reasons, both structural:

* **The corpus has no starter.**  ``gs2_wall_none.xlsx`` is a benchmark in its
  own right, locked at its own factor of safety; the tutorial needs a file that
  is the *completed* model minus exactly the one thing the page teaches (the
  wall row), so that "add one row and re-run" is literally what the reader does.
  The starter is therefore cut from ``gs2_wall.xlsx``, and this builder asserts
  that what remains is the same model ``gs2_wall_none.xlsx`` carries.
* **The corpus files carry sidecars.**  Both ship a committed mesh and a
  committed seepage field, because a ``u = 'seep'`` model has to run its
  stability solve on the mesh its nodal field was written on.  Tutorial
  downloads are sidecar-free (the reader builds the mesh, runs the seepage and
  runs the SSRM), so the tutorial copies must not sit beside a ``_mesh.json`` or
  a ``_seep.csv``.

The tutorial copies are written at the current template version; the corpus pair
is older, and rebuilding the corpus is the corpus builder's business, not this
one's.

What the starter leaves out, and what it keeps
----------------------------------------------
Out: the single ``piles`` sheet row.  That is the page's whole teaching content
for this half — *which* cells a continuous member fills (x1 y1 x2 y2, E, Area,
I, S = 1) and which it correctly leaves blank (D, Vcap, Mcap), and the
re-mesh / re-seep / re-run loop that adding a constraint line forces.

Kept: everything else, including the mesh declaration (tri6 at 1.5 m with the
0.3 m local override across the weak clay band).  The mesh is not what this page
teaches and the reader builds it on *both* models, so both files open their
Build Mesh dialog on the settings the page's numbers were measured at.  This is
the one place the pair differs from FEM-2's, whose starter withholds the mesh
because meshing is part of what FEM-2 teaches.

Run:  PYTHONPATH=. python3 tools/build_pile_wall_tutorial.py         # both
      PYTHONPATH=. python3 tools/build_pile_wall_tutorial.py start   # one
"""

from __future__ import annotations

import math
import os
import sys
import warnings

warnings.filterwarnings("ignore")

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO_ROOT)

from xslope.fileio import (default_template_path, load_slope_data,  # noqa: E402
                           save_slope_data_to_xlsx)

CORPUS = os.path.join(REPO_ROOT, "docs", "verification", "files", "geostudio")
WALL_BASE = os.path.join(CORPUS, "gs2_wall.xlsx")
NONE_BASE = os.path.join(CORPUS, "gs2_wall_none.xlsx")
TUTORIAL_FILES = os.path.join(REPO_ROOT, "docs", "tutorials", "files")

START_OUT = "xslope_pile_wall_start.xlsx"
DONE_OUT = "xslope_pile_wall.xlsx"

#: The mesh both models are solved on, and the one the corpus builder emits:
#: tri6 (quadratic elements are the only ones the SSRM path is verified on) at
#: 1.5 m globally, with the weak clay band's own 0.3 m size override riding
#: through on its profile line.
ELEMENT_TYPE = "tri6"
TARGET_SIZE = 1.5

#: The model fields the starter must share with the corpus's own no-wall model.
#: Checked rather than assumed: if cutting the wall row out of gs2_wall does not
#: land on gs2_wall_none, the two halves of the page are not the same slope.
SHARED_KEYS = ("unit_system", "gamma_water", "max_depth", "k_seismic",
               "tcrack_depth")
MATERIAL_KEYS = ("name", "option", "gamma", "gamma_sat", "c", "phi", "E", "nu",
                 "u", "t_cut", "k1", "k2", "alpha", "kr0", "h0")


def _same(a, b, tol=1e-9):
    if isinstance(a, float) and isinstance(b, float):
        if a != a and b != b:          # both unset (NaN): the same blank
            return True
        return math.isclose(a, b, rel_tol=tol, abs_tol=tol)
    return a == b


def _base():
    """The corpus wall model, with its committed mesh dropped.

    ``load_slope_data`` picks up ``gs2_wall_mesh.json`` and ``gs2_wall_seep.csv``
    from beside the corpus workbook.  Neither tutorial file has a sidecar, so the
    in-memory mesh and the nodal pore pressure field go too, and nothing written
    here can inherit either.
    """
    sd = load_slope_data(WALL_BASE)
    sd.pop("mesh", None)
    sd.pop("seep_u", None)
    if len(sd["pile_lines"]) != 1:
        raise SystemExit(f"{os.path.basename(WALL_BASE)} carries "
                         f"{len(sd['pile_lines'])} pile rows, expected 1")
    sd["element_type"] = ELEMENT_TYPE
    sd["target_size"] = TARGET_SIZE
    return sd


def _check_starter_is_the_corpus_no_wall_model(sd):
    """The starter, cut from the wall model, must be the corpus's no-wall model."""
    none = load_slope_data(NONE_BASE)
    none.pop("mesh", None)
    none.pop("seep_u", None)
    if none["pile_lines"]:
        raise SystemExit(f"{os.path.basename(NONE_BASE)} carries a pile row")
    if str(sd["ground_surface"]) != str(none["ground_surface"]):
        raise SystemExit("the corpus pair does not share a ground surface")
    for k in SHARED_KEYS:
        if not _same(sd[k], none[k]):
            raise SystemExit(f"{k} differs across the corpus pair "
                             f"({sd[k]!r} vs {none[k]!r})")
    if len(sd["materials"]) != len(none["materials"]):
        raise SystemExit("the corpus pair carries different material counts")
    for a, b in zip(sd["materials"], none["materials"]):
        for k in MATERIAL_KEYS:
            if not _same(a.get(k), b.get(k)):
                raise SystemExit(f"material {a['name']}: {k} differs across the "
                                 f"corpus pair ({a.get(k)!r} vs {b.get(k)!r})")
    if len(sd["profile_lines"]) != len(none["profile_lines"]):
        raise SystemExit("the corpus pair carries different profile line counts")
    for a, b in zip(sd["profile_lines"], none["profile_lines"]):
        if a["mat_id"] != b["mat_id"] or list(a["coords"]) != list(b["coords"]) \
                or not _same(a["size"], b["size"]):
            raise SystemExit("a profile line differs across the corpus pair")
    if sd["seepage_bc"] != none["seepage_bc"]:
        raise SystemExit("the seepage boundary conditions differ across the "
                         "corpus pair")


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
    """The starter: the whole slope — geometry, the three zones including the
    elastic OC soil, the tensile cut-off on both Mohr-Coulomb zones, the seepage
    boundary conditions, the weak-clay size override — and no wall.

    This is also the page's first run: mesh it, solve the seepage on that mesh,
    run the SSRM, and read the slope's own factor of safety.
    """
    sd = _base()
    sd["pile_lines"] = []
    _check_starter_is_the_corpus_no_wall_model(sd)
    return _write(sd, START_OUT)


def build_done():
    """The completed file: the same model with the one wall row the reader adds.

    The row is the corpus model's, unchanged — a vertical member from (38, 10) to
    (38, 1) with E = 2x10^8 kPa, Area = 0.02 m2/m and I = 0.0005 m4/m at S = 1,
    and D, Vcap and Mcap left blank.  S = 1 because the wall is continuous out of
    plane, so the section constants are already per metre and XSLOPE's spacing
    divisor must not reduce them; D is blank because there is no gap between
    members for an arching theory to act in, so no Ito & Matsui limit pressure
    applies; the capacities are blank because the corpus transcription carries no
    structural rating for the section.
    """
    sd = _base()
    return _write(sd, DONE_OUT)


TARGETS = {"start": build_start, "done": build_done}


def main(argv):
    for name in argv[1:] or list(TARGETS):
        if name not in TARGETS:
            raise SystemExit(f"unknown target {name!r}; choices: {list(TARGETS)}")
        TARGETS[name]()


if __name__ == "__main__":
    main(sys.argv)
