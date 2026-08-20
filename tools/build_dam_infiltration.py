"""Build Tutorial SEEP-4's file PAIR (infiltration and flux boundaries).

The model is the Fredlund & Rahardjo (1993) earth dam of the Rocscience
groundwater manual's Problem 6 — 12 m high, 4 m crest, symmetric 2:1 faces,
52 m base, reservoir at 10 m, a 12 m horizontal drain at the downstream toe —
under steady rainfall infiltration.  That is verification case 4 of
[GW6](../docs/verification/rocscience_groundwater.md#gw6), and this builder
DERIVES both tutorial files from the committed corpus file rather than
restating its numbers:

    docs/verification/files/rocscience_gw/gw006d.xlsx   (the base)
      -> docs/tutorials/files/xslope_dam_infiltration_start.xlsx
      -> docs/tutorials/files/xslope_dam_infiltration.xlsx

so the geometry, the material and the boundary values a reader types can never
drift from the ones the verification page measures.  gw006d itself is written by
``benchmarks/rocscience/build_groundwater.py``; nothing here hand-edits an xlsx.

The base's own dry-weather sibling, gw006a (case 1), carries an IDENTICAL
material and geometry and differs only by having no flux boundaries, so the
tutorial's rain / no-rain contrast is the two corpus cases side by side.

What the tutorial pair changes from the base, and why
----------------------------------------------------
* **The placeholder failure circle is dropped.**  Every groundwater-corpus file
  carries one so the workbook has a surface in it; these two are seepage models
  and nothing solves a factor of safety on them.  Dropping it is also what makes
  the STARTER load only as an editor file: with no surface, no mesh and no
  boundary conditions, ``load_slope_data`` refuses it as not yet runnable, which
  is exactly the state the tutorial starts the reader in.
* **The starter's ``seep bc`` sheet is emptied.**  Head, fluxes and exit face all
  go, so the reader enters every boundary condition on the page — the whole
  point of the tutorial — rather than correcting a partly filled sheet.

* **The infiltration covers the whole exposed surface**, waterline (20, 10) over
  the crest to the downstream toe (52, 0), where the corpus file stops one vendor
  mesh edge short at each end — (22, 11) to (50, 1).  Those vendor extents are
  tied to the nodes of the vendor's own mesh, and the verification file transcribes
  them exactly because reproducing a published answer means reproducing the model
  that produced it.  A tutorial teaches the modelling instead: rain falls on the
  surface, so the boundary is drawn on the GEOMETRY, corner to corner, and the
  solver resolves what happens where it meets the reservoir head and the toe
  drain (the specified head wins at a shared node — see docs/seep/overview.md,
  "Specified flux (Neumann)").  The rates are the base file's own, not restated
  here, so the two files can only ever differ in extent.

Everything else — the profile line, the material (k = 1e-7 m/s isotropic, the
Mualem–van Genuchten fit α = 0.2452 / n = 2.5739 shared by all five GW6 cases),
SI units, and, in the completed file, the reservoir head and the toe-drain exit
face — is the base file's, verbatim.

The infiltration goes in as three blocks rather than one because the rain rate
and the normal flux a boundary takes are not the same number.  The vendor applies
a VERTICAL 1e-8 m/s over the exposed surface; XSLOPE's flux boundary condition
takes the normal Darcy velocity (see docs/seep/overview.md, "Specified flux
(Neumann)"), which on a 2:1 face is 1e-8 x 2/sqrt(5) = 8.94427191e-9 and across
the horizontal crest is the full 1e-8.  The three blocks are contiguous — each
one starts where the last ended — so the coverage has no gap; they are three
entries only because a flux entry carries one rate.

Run:  PYTHONPATH=. python3 tools/build_dam_infiltration.py           # both
      PYTHONPATH=. python3 tools/build_dam_infiltration.py start     # one
"""

from __future__ import annotations

import os
import sys

from xslope.fileio import (load_slope_data, save_slope_data_to_xlsx,
                           default_template_path)

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BASE = os.path.join(REPO_ROOT, "docs", "verification", "files",
                    "rocscience_gw", "gw006d.xlsx")
TUTORIAL_FILES = os.path.join(REPO_ROOT, "docs", "tutorials", "files")

START_OUT = "xslope_dam_infiltration_start.xlsx"
DONE_OUT = "xslope_dam_infiltration.xlsx"

#: The corpus file's LEM placeholder circle, dropped from both tutorial files.
_EMPTY_SURFACES = dict(circles=[], non_circ=[])

#: The exposed surface of the dam, waterline to downstream toe, as the vertices the
#: rain boundary is drawn through: up the upstream face above the reservoir, across
#: the crest, down the downstream face.  Consecutive pairs become the flux blocks.
FLUX_VERTICES = [(20.0, 10.0), (24.0, 12.0), (28.0, 12.0), (52.0, 0.0)]


def _model():
    """The base model with its placeholder failure surface removed.

    Both tutorial files come through here, so the soil and the section are the
    same object in both, and both are the corpus file's.
    """
    sd = load_slope_data(BASE)
    sd.update(_EMPTY_SURFACES)
    return sd


def _write(sd, filename):
    out = os.path.join(TUTORIAL_FILES, filename)
    save_slope_data_to_xlsx(sd, out, template=default_template_path())
    print(f"wrote {os.path.relpath(out, REPO_ROOT)}")
    return out


def build_start():
    """The starter: the dam, its soil and its units, with no boundary at all.

    Not a partly filled sheet — the reservoir head, the three infiltration
    polylines and the toe drain are all the reader's to enter, and until they
    are the file has nothing to solve.
    """
    sd = _model()
    sd["seepage_bc"] = {"specified_heads": [], "specified_fluxes": [],
                        "exit_face": []}
    sd["seepage_bc2"] = {"specified_heads": [], "specified_fluxes": [],
                         "exit_face": []}
    sd["has_seepage_bc2"] = False
    return _write(sd, START_OUT)


def _rates(base_fluxes):
    """The base file's two infiltration rates, read off the base file itself.

    A horizontal block takes the vertical rain unprojected, so the crest rate is
    the one whose block is level; a sloping block takes the rain's component
    normal to it, which is smaller.  Reading them this way rather than by list
    position means a change to the corpus file's rates travels here, and a change
    to its block ORDER does not break the build silently.
    """
    level = [f["flux"] for f in base_fluxes if f["coords"][0][1] == f["coords"][-1][1]]
    sloping = [f["flux"] for f in base_fluxes if f["coords"][0][1] != f["coords"][-1][1]]
    if not level or not sloping:
        raise SystemExit(f"{os.path.basename(BASE)} no longer carries one level "
                         "and one sloping infiltration block; the tutorial "
                         "builder cannot read its rates")
    return max(sloping), max(level)


def _full_extent_fluxes(base_fluxes, vertices=FLUX_VERTICES):
    """The rain redrawn over the whole exposed surface, at the base file's rates.

    One block per segment of ``vertices``, each taking the level rate where the
    segment is horizontal and the sloping rate where it is not.  The blocks are
    contiguous by construction, so the coverage is gapless from the waterline to
    the toe.
    """
    sloping, level = _rates(base_fluxes)
    blocks = []
    for a, b in zip(vertices, vertices[1:]):
        blocks.append({"flux": level if a[1] == b[1] else sloping,
                       "coords": [a, b]})
    return blocks


def build_done():
    """The completed file: the base's boundaries, with the rain redrawn on the
    geometry rather than on the vendor's mesh edges."""
    sd = _model()
    bc = dict(sd["seepage_bc"])
    bc["specified_fluxes"] = _full_extent_fluxes(bc["specified_fluxes"])
    sd["seepage_bc"] = bc
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
