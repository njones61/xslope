"""Build Tutorial LEM-13's rock-slope PAIR (the Hoek-Brown models).

LEM-13 runs the generalized Hoek-Brown strength option at the two ends of the
criterion, and both models are DERIVED from the committed verification corpus
rather than restated here:

    docs/verification/files/rocscience/hammah_hb1.xlsx
      -> docs/tutorials/files/xslope_rock_slope.xlsx        (Part A)

    docs/verification/files/rocscience/rs2_60c.xlsx
      -> docs/tutorials/files/xslope_rock_slope_li.xlsx     (Part B)

so the section, the four Hoek-Brown inputs and the unit weight a reader reads off
the page can never drift from the ones
[the verification page](../docs/verification/rs2.md#hoek-brown) measures. Nothing
here hand-edits an xlsx; both files are written through the package writer at the
current template version.

Part A — Hammah, Yacoub, Corkum & Curran (2005), ARMA/USRMS 05-810, Example 1.
A 10 m, 45 degree slope in a badly broken rock mass: sigma_ci = 30 MPa entered as
30,000 kPa, GSI = 5, mi = 2, D = 0, gamma = 25 kN/m3, with E = 5,000 MPa and
nu = 0.3 for the strength-reduction half of the page. GSI = 5 puts the exponent
a at 0.619, far from the classical 0.5, which is what makes the case a test of
the criterion rather than of the geometry. Locked at Bishop 1.150 / Spencer 1.152
in the LEM and SSRM 1.159 on tri6 elements at a 0.9 m target size.

Part B — Li, Merifield & Lyamin (2008), IJRMMS 45, 689-700, the beta = 45 degree
case. The same criterion at GSI = 70 and mi = 15: a strong, lightly jointed rock
mass whose exponent a is 0.501. The problem is NORMALIZED — H = 1 m, so gamma*H
is 23 kPa and the critical strength ratio puts sigma_ci at 4.37 kPa. Those
magnitudes carry the page's units trap, since Hoek-Brown convention invites MPa.
Locked at Spencer 1.035.

What the tutorial files change from the corpus, and why
------------------------------------------------------
* **The corpus's own starting circle is kept.** Both pages are open-and-run: the
  reader downloads a completed model and never places a circle, so the seed's
  only job is to put the search on the answer the verification page measures. It
  does that and the standing Xo-mid-slope / Yo-at-twice-the-height / Depth-circle
  seed does not, quite: on Li's 1 m section that seed refines onto a marginally
  different circle and reports 1.033 against the corpus's 1.035, and on Hammah's
  it reports Bishop 1.151 against the corpus's 1.150. Both are the same
  mechanism at the same daylight point, but a tutorial that printed 1.033 beside
  a verification page printing 1.035 would read as a disagreement between two
  pages about one model.
* **gamma_sat is cleared.** Neither model has a water table, and a saturated
  unit weight that can never apply raises a model check on a page whose runs are
  meant to open clean.
* **Part A declares the finite element run it is solved on** — tri6 at a 0.9 m
  target size, the strength-reduction bracket 0.8 to 1.6, and K0 = 1 — so
  Studio's Build Mesh and Run FEM dialogs open on the tutorial's own run instead
  of on the defaults. Those are the settings the SSRM lock tag names. Neither
  file carries a mesh sidecar: the page builds the mesh in Studio.
* **Part B declares no finite element settings.** SSRM is not locked on Li's
  problem, and the page runs it in the limit equilibrium engine only.

Run:  PYTHONPATH=. python3 tools/build_rock_slope.py          # both
      PYTHONPATH=. python3 tools/build_rock_slope.py li       # one
"""

from __future__ import annotations

import os
import sys
import warnings

warnings.filterwarnings("ignore")

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)

from xslope.fileio import (default_template_path, load_slope_data,  # noqa: E402
                           save_slope_data_to_xlsx)

CORPUS = os.path.join(REPO_ROOT, "docs", "verification", "files", "rocscience")
TUTORIAL_FILES = os.path.join(REPO_ROOT, "docs", "tutorials", "files")

HAMMAH_BASE = os.path.join(CORPUS, "hammah_hb1.xlsx")
LI_BASE = os.path.join(CORPUS, "rs2_60c.xlsx")

HAMMAH_OUT = "xslope_rock_slope.xlsx"
LI_OUT = "xslope_rock_slope_li.xlsx"

#: Part A's finite element run, from the SSRM lock tag at
#: docs/verification/rs2.md (element_type=tri6, target_size=0.9, f_min=0.8,
#: f_max=1.6, k0=1). tri6 because a factor of safety needs quadratic elements.
ELEMENT_TYPE = "tri6"
TARGET_SIZE = 0.9
SSRM_F_MIN, SSRM_F_MAX = 0.8, 1.6
K0 = 1.0


def _base(path):
    """The corpus model with the one change both tutorial files share: no saturated
    unit weight on a section that has no water in it."""
    sd = load_slope_data(path)
    sd.pop("mesh", None)
    sd.pop("seep_u", None)
    if len(sd["materials"]) != 1:
        raise SystemExit(f"{os.path.basename(path)} carries "
                         f"{len(sd['materials'])} materials, expected 1")
    mat = dict(sd["materials"][0])
    if str(mat.get("option", "")).strip().lower() != "hb":
        raise SystemExit(f"{os.path.basename(path)} is not a Hoek-Brown model")
    mat["gamma_sat"] = None
    sd["materials"] = [mat]
    sd["circular"] = True
    if len(sd["circles"]) != 1:
        raise SystemExit(f"{os.path.basename(path)} carries "
                         f"{len(sd['circles'])} starting circles, expected 1")
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


def build_hammah():
    """Part A: the weak-rock slope, carrying the finite element run as well as the
    limit equilibrium one, since it is the page's two-engine comparison."""
    sd = _base(HAMMAH_BASE)
    sd["element_type"] = ELEMENT_TYPE
    sd["target_size"] = TARGET_SIZE
    sd["ssrm_f_min"] = SSRM_F_MIN
    sd["ssrm_f_max"] = SSRM_F_MAX
    sd["k0"] = K0
    return _write(sd, HAMMAH_OUT)


def build_li():
    """Part B: the strong-rock slope, limit equilibrium only."""
    return _write(_base(LI_BASE), LI_OUT)


TARGETS = {"hammah": build_hammah, "li": build_li}


def main(argv):
    for name in argv[1:] or list(TARGETS):
        if name not in TARGETS:
            raise SystemExit(f"unknown target {name!r}; choices: {list(TARGETS)}")
        TARGETS[name]()


if __name__ == "__main__":
    main(sys.argv)
