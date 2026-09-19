"""Build Tutorial FEM-3's files: the block wall, and part 3's four sheet models.

FEM-3 puts slip joints into a structure.  Its models are built here rather than
derived from the verification corpus: the wall is an ordinary segmental wall
rather than a published case, and part 3's pairs exist to be run twice each with
one cell changed.

    docs/tutorials/files/xslope_block_wall_start.xlsx      the section and the
                                                           four materials, with
                                                           the joints sheet empty
    docs/tutorials/files/xslope_block_wall.xlsx            part 1: the wall on
                                                           its seven joints
    docs/tutorials/files/xslope_block_wall_grid.xlsx       part 2: three geogrid
                                                           layers tied into the
                                                           blocks
    docs/tutorials/files/xslope_base_geotextile_bonded.xlsx
    docs/tutorials/files/xslope_base_geotextile_jointed.xlsx
    docs/tutorials/files/xslope_liner_bonded.xlsx
    docs/tutorials/files/xslope_liner_jointed.xlsx         part 3: two pairs,
                                                           each differing by the
                                                           reinforce sheet's
                                                           Joint column alone

Nothing here hand-edits an xlsx: each file is written through the package writer
at the current template version, then re-loaded and put through the finite
element model checks.

The wall
--------
Six 0.6 m courses of 1.2 m deep modular block — 3.6 m of face — on a 3.2 m
foundation, with a 6 m reinforced fill zone behind the blocks and a 2:1 backfill
slope rising 2 m above the crest, on a 24 m section.

The block depth is 1.2 m and not the 0.6 m of a smaller unit because 1.2 m is
what lets the wall stand on its own: the same wall on 0.6 m deep block has no
equilibrium at any sweep budget, and part 1 needs a wall that stands before part
2 can say what the geogrid adds to it.  The face is vertical: a real segmental
wall is battered, and a battered wall has a stepped back face, which turns seven
joint lines into twelve without changing what the page teaches.

Seven joint lines: one under the base of the block column, one on its back face
against the fill, and one between each pair of courses.  The blocks are declared
elastic, so the strength reduction divides the joints' strength and the soil's
and has nothing in the block itself to divide.

The geogrid of part 2 is three layers on the 0.6 m, 1.8 m and 3.0 m course
lines, 3.0 m long, tied into the block column at 40 kN/m.  Each is a reinforce
line with ``Joint = Yes``: a sheet whose front end stops on the back face of the
blocks touches that joint line, and a bonded bar standing on a node the mesh
split has copied once per wedge of material has no defined side to attach to, so
both the model checks and the mesher refuse the bonded version of this layout by
name.  That is why the bonded-against-jointed comparison of part 3 is made on
models that carry no other joint line.

Part 3's two pairs
------------------
A base geotextile under an embankment on soft clay, where the failure surface
crosses the sheet, and a smooth geomembrane liner under the same embankment,
where the fill can slide out along it.  Each is built twice, bonded and jointed,
and the two files of a pair differ in one cell.

Run:  PYTHONPATH=. python3 tools/build_block_wall.py           # all seven
      PYTHONPATH=. python3 tools/build_block_wall.py wall      # one
"""

from __future__ import annotations

import copy
import os
import sys

from shapely.geometry import Polygon

from xslope.fileio import (default_template_path, load_slope_data,
                           save_slope_data_to_xlsx)
from xslope.preflight import preflight

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(REPO, "docs", "tutorials", "files")
BASE_WB = os.path.join(OUT, "xslope_ssrm_embankment.xlsx")

# ---- the wall ---------------------------------------------------------------
WALL_X0 = 8.0                       # front face of the block column
BLOCK_DEPTH = 1.2
WALL_X1 = WALL_X0 + BLOCK_DEPTH     # back face
COURSE = 0.6
N_COURSE = 6
WALL_TOP = COURSE * N_COURSE        # 3.6
REINF_ZONE = 6.0                    # depth of the reinforced fill behind the wall
REINF_BACK = WALL_X1 + REINF_ZONE   # back of the reinforced fill zone
BACKSLOPE_RISE = 2.0                # the backfill slope rises this far above the crest
CREST_X = WALL_X1 + 2.0 * BACKSLOPE_RISE    # 2:1
CREST_Y = WALL_TOP + BACKSLOPE_RISE
DOMAIN_X1 = 24.0
FOUND_Y0 = -3.2

WALL_TARGET_SIZE = 0.8              # m, global
BLOCK_SIZE = 0.3                    # m, on each block polygon

PHI_BASE, PHI_BACK, PHI_COURSE = 34.0, 30.0, 35.0

SHEET_Y = (0.6, 1.8, 3.0)           # the course lines the layers sit on
SHEET_LEN = 3.0
SHEET_TMAX = 40.0
SHEET_TIE = 40.0                    # connection capacity at the block column

# ---- part 3's embankment ----------------------------------------------------
EMB_H = 5.0
EMB_TOE_L, EMB_CREST_L, EMB_CREST_R, EMB_TOE_R = 14.0, 24.0, 36.0, 46.0
EMB_X1 = 60.0
EMB_TARGET_SIZE = 1.2


def _base():
    return load_slope_data(BASE_WB)


def mat(template, name, gamma, E, nu, c=0.0, phi=0.0, option="mc"):
    m = copy.deepcopy(template)
    m.update(name=name, gamma=gamma, gamma_sat=None, option=option,
             c=c, phi=phi, E=E, nu=nu, u="none", ru=0.0,
             t_cut=(0.0 if option == "mc" else None))
    for k in ("cp", "r_elev", "d", "psi", "sigma_gamma", "sigma_c", "sigma_phi",
              "sigma_cp", "sigma_d", "sigma_psi"):
        m[k] = 0.0
    return m


def poly(coords, mat_id, size=None):
    return {"polygon": Polygon(coords), "mat_id": mat_id, "size": size}


def joint(label, p1, p2, phi, c=0.0):
    nan = float("nan")
    return {"label": label, "x1": p1[0], "y1": p1[1], "x2": p2[0], "y2": p2[1],
            "c": c, "phi": phi, "c_res": nan, "phi_res": nan, "dil": nan,
            "t_cut": 0.0, "kn": nan, "ks": nan, "jred": ""}


def sheet(label, p1, p2, t_max, adhesion, delta, tend1=0.0, tend2=0.0,
          E=1.0e6, area=1.0e-3, is_joint=""):
    nan = float("nan")
    return {"label": label, "x1": p1[0], "y1": p1[1], "x2": p2[0], "y2": p2[1],
            "t_max": t_max, "t_res": nan, "lp1": 0.0, "lp2": 0.0,
            "E": E, "area": area, "type": "Geosynthetic", "dir": "Tangent",
            "appl": "Active", "tend1": tend1, "tend2": tend2, "spacing": 1.0,
            "joint": is_joint, "kn": nan, "ks": nan, "jred": "",
            "adhesion": adhesion, "delta": delta}


def model(materials, polygons, circle, target_size, joint_lines=(),
          reinforcement_lines=(), max_depth=0.0):
    sd = _base()
    sd["materials"] = list(materials)
    sd["profile_lines"] = []
    sd["polygons"] = list(polygons)
    sd["joint_lines"] = list(joint_lines)
    sd["reinforcement_lines"] = list(reinforcement_lines)
    for key in ("reinforce_lines", "joint_zones", "ssr_zones", "refine_zones",
                "non_circ", "dloads", "dloads2", "dload_dirs", "dload2_dirs",
                "pile_lines", "line_loads", "piezo_line", "piezo_line2"):
        sd[key] = []
    sd["circles"] = [dict(circle)]
    sd["element_type"] = "tri6"
    sd["target_size"] = target_size
    sd["element_size_1d"] = None
    sd["unit_system"] = "si"
    sd["gamma_water"] = 9.81
    sd["max_depth"] = max_depth
    sd["water_loads"] = "manual"
    sd["mesh"] = None
    sd["k0"] = None
    sd["tension_srf"] = None
    sd["ssrm_f_min"] = None
    sd["ssrm_f_max"] = None
    return sd


# ---- the wall ---------------------------------------------------------------

def wall_materials(template):
    return [
        mat(template, "foundation", gamma=20.0, c=15.0, phi=30.0,
            E=40000.0, nu=0.3),
        mat(template, "reinforced fill", gamma=20.0, c=0.0, phi=36.0,
            E=50000.0, nu=0.3),
        mat(template, "retained fill", gamma=19.0, c=5.0, phi=28.0,
            E=25000.0, nu=0.3),
        mat(template, "block", gamma=23.0, E=1.0e7, nu=0.2, option="elastic"),
    ]


def wall_polygons():
    polys = [
        poly([(0.0, FOUND_Y0), (DOMAIN_X1, FOUND_Y0), (DOMAIN_X1, 0.0),
              (0.0, 0.0)], 0),
        poly([(WALL_X1, 0.0), (REINF_BACK, 0.0), (REINF_BACK, WALL_TOP),
              (WALL_X1, WALL_TOP)], 1),
        poly([(WALL_X1, WALL_TOP), (CREST_X, CREST_Y), (DOMAIN_X1, CREST_Y),
              (DOMAIN_X1, 0.0), (REINF_BACK, 0.0), (REINF_BACK, WALL_TOP)], 2),
    ]
    # One polygon per course, each with a local Size: a 0.6 m course cannot be
    # resolved at the global target size, and a local override costs far fewer
    # nodes than refining the whole 24 m section would.
    for i in range(N_COURSE):
        y0 = i * COURSE
        polys.append(poly([(WALL_X0, y0), (WALL_X1, y0),
                           (WALL_X1, y0 + COURSE), (WALL_X0, y0 + COURSE)],
                          3, size=BLOCK_SIZE))
    return polys


def wall_joints():
    j = [joint("base", (WALL_X0, 0.0), (WALL_X1, 0.0), PHI_BASE),
         joint("back face", (WALL_X1, 0.0), (WALL_X1, WALL_TOP), PHI_BACK)]
    for i in range(1, N_COURSE):
        y = i * COURSE
        j.append(joint(f"course-{i:02d}", (WALL_X0, y), (WALL_X1, y),
                       PHI_COURSE))
    return j


def wall_sheets():
    return [sheet(f"grid-{n + 1:02d}", (WALL_X1, y), (WALL_X1 + SHEET_LEN, y),
                  t_max=SHEET_TMAX, adhesion=1.0, delta=30.0,
                  tend1=SHEET_TIE, tend2=0.0, is_joint="Yes")
            for n, y in enumerate(SHEET_Y)]


WALL_CIRCLE = {"Xo": 10.6, "Yo": 11.2, "Depth": 0.0, "R": 11.2}


def _wall(joint_lines=(), sheets=()):
    sd = _base()
    return model(wall_materials(sd["materials"][0]), wall_polygons(),
                 WALL_CIRCLE, WALL_TARGET_SIZE, joint_lines=joint_lines,
                 reinforcement_lines=sheets, max_depth=FOUND_Y0)


def build_start():
    """The section, the four materials and the six block polygons, with nothing
    on the joints sheet.  The reader types the seven lines into it."""
    return _wall()


def build_wall():
    """Part 1: the wall standing on its own seven contacts."""
    return _wall(joint_lines=wall_joints())


def build_wall_grid():
    """Part 2: the same wall with three geogrid layers tied into the blocks."""
    return _wall(joint_lines=wall_joints(), sheets=wall_sheets())


# ---- part 3 -----------------------------------------------------------------

#: Cohesion of the part 3 embankment fill, kPa. See build_base_geotextile.
EMB_FILL_C = 5.0


def emb_polygons(clay_y0):
    return [
        poly([(0.0, clay_y0), (EMB_X1, clay_y0), (EMB_X1, 0.0), (0.0, 0.0)], 0),
        poly([(EMB_TOE_L, 0.0), (EMB_TOE_R, 0.0),
              (EMB_CREST_R, EMB_H), (EMB_CREST_L, EMB_H)], 1),
    ]


def build_base_geotextile(is_joint):
    """An embankment on soft clay over a base geotextile.  The critical surface
    cuts up through the fill and crosses the sheet, so the sheet carries tension
    across it either way it is modeled.

    The fill carries EMB_FILL_C of cohesion: a cohesionless fill at 34 degrees
    on 2:1 slopes fails on its own face at tan(34)/tan(26.6) = 1.35, shallower
    than anything through the clay, and that face governed every run of the
    first cut of these models without touching the sheet."""
    sd = _base()
    t = sd["materials"][0]
    mats = [mat(t, "soft clay", gamma=17.0, c=20.0, phi=0.0, E=8000.0, nu=0.35),
            mat(t, "embankment fill", gamma=20.0, c=EMB_FILL_C, phi=34.0,
                E=25000.0, nu=0.3)]
    line = [sheet("base geotextile", (EMB_TOE_L, 0.0), (EMB_TOE_R, 0.0),
                  t_max=100.0, adhesion=5.0, delta=20.0, E=2.0e6,
                  is_joint=is_joint)]
    return model(mats, emb_polygons(-4.0),
                 {"Xo": 30.0, "Yo": 15.0, "Depth": -4.0, "R": 19.0},
                 EMB_TARGET_SIZE, reinforcement_lines=line, max_depth=-4.0)


def build_liner(is_joint):
    """The same embankment on a smooth geomembrane liner over a firm
    foundation.  The interface is the weakest thing in the section, so the fill
    can slide out along the liner — which a bonded bar cannot represent."""
    sd = _base()
    t = sd["materials"][0]
    mats = [mat(t, "foundation", gamma=20.0, c=30.0, phi=32.0,
                E=50000.0, nu=0.3),
            mat(t, "embankment fill", gamma=20.0, c=EMB_FILL_C, phi=34.0,
                E=25000.0, nu=0.3)]
    line = [sheet("liner", (EMB_TOE_L, 0.0), (EMB_TOE_R, 0.0),
                  t_max=50.0, adhesion=0.5, delta=10.0, E=2.0e6,
                  is_joint=is_joint)]
    # The same 4 m foundation as the geotextile pair: one section, one sketch.
    return model(mats, emb_polygons(-4.0),
                 {"Xo": 30.0, "Yo": 15.0, "Depth": -4.0, "R": 19.0},
                 EMB_TARGET_SIZE, reinforcement_lines=line, max_depth=-4.0)


BUILDS = {
    "start": ("xslope_block_wall_start.xlsx", build_start),
    "wall": ("xslope_block_wall.xlsx", build_wall),
    "grid": ("xslope_block_wall_grid.xlsx", build_wall_grid),
    "sheet_bonded": ("xslope_base_geotextile_bonded.xlsx",
                     lambda: build_base_geotextile("No")),
    "sheet_jointed": ("xslope_base_geotextile_jointed.xlsx",
                      lambda: build_base_geotextile("Yes")),
    "liner_bonded": ("xslope_liner_bonded.xlsx", lambda: build_liner("No")),
    "liner_jointed": ("xslope_liner_jointed.xlsx", lambda: build_liner("Yes")),
}


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    names = argv or list(BUILDS)
    print("template:", default_template_path())
    print(f"wall: {N_COURSE} courses of {COURSE:g} m on {BLOCK_DEPTH:g} m deep "
          f"block, face x = {WALL_X0:g} to {WALL_X1:g}, top y = {WALL_TOP:g}")
    bad = 0
    for name in names:
        if name not in BUILDS:
            print(f"unknown model {name!r}; known: {', '.join(BUILDS)}")
            return 2
        fname, fn = BUILDS[name]
        path = os.path.join(OUT, fname)
        save_slope_data_to_xlsx(fn(), path)
        back = load_slope_data(path)
        rep = preflight(back, "fem")
        errs = [f for f in rep.findings if f.severity == "error"]
        warns = [f for f in rep.findings if f.severity == "warning"]
        print(f"\n{fname}")
        print(f"  joint lines {len(back['joint_lines'])}, "
              f"reinforce lines {len(back['reinforcement_lines'])}, "
              f"polygons {len(back['polygons'])}, "
              f"materials {len(back['materials'])}")
        for f in errs:
            print(f"  ERROR   {f.rule_id}: {f.message}")
        for f in warns:
            print(f"  warning {f.rule_id}: {f.message}")
        bad += len(errs)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.path.insert(0, REPO)
    sys.exit(main())
