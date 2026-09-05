"""Build Tutorial LEM-13's rock-slope model (the Hoek-Brown tutorial file).

LEM-13 runs the generalized Hoek-Brown strength option on one section, and that
model is DERIVED from the committed verification corpus rather than restated here:

    docs/verification/files/rocscience/hammah_hb1.xlsx
      -> docs/tutorials/files/xslope_rock_slope.xlsx

so the section, the four Hoek-Brown inputs and the unit weight a reader reads off
the page can never drift from the ones
[the verification page](../docs/verification/rs2.md#hoek-brown) measures. Nothing
here hand-edits an xlsx; the file is written through the package writer at the
current template version.

The problem — Hammah, Yacoub, Corkum & Curran (2005), ARMA/USRMS 05-810,
Example 1.
A 10 m, 45 degree slope in a badly broken rock mass: sigma_ci = 30 MPa entered as
30,000 kPa, GSI = 5, mi = 2, D = 0, gamma = 25 kN/m3, with E = 5,000 MPa and
nu = 0.3 for the strength-reduction half of the page. GSI = 5 puts the exponent
a at 0.619, far from the classical 0.5, which is what makes the case a test of
the criterion rather than of the geometry. Locked at Bishop 1.150 / Spencer 1.152
in the LEM and SSRM 1.153 on tri6 elements at a 0.9 m target size.

What the tutorial file changes from the corpus, and why
-------------------------------------------------------
* **The corpus's own starting circle is kept.** The page is open-and-run: the
  reader downloads a completed model and never places a circle, so the seed's only
  job is to put the search on the answer the verification page measures. It does
  that and the standing Xo-mid-slope / Yo-at-twice-the-height / Depth-circle seed
  does not, quite — that seed reports Bishop 1.151 against the corpus's 1.150, the
  same mechanism at the same daylight point but a marginally different circle, and
  a tutorial printing 1.151 beside a verification page printing 1.150 would read as
  two pages disagreeing about one model.
* **gamma_sat is cleared.** The model has no water table, and a saturated
  unit weight that can never apply raises a model check on a page whose runs are
  meant to open clean.
* **The tensile cutoff is declared at zero**, the value the rest of the shipped
  files carry, so Run FEM does not open on an unbounded-tension warning. The rock's
  own Hoek-Brown tensile strength is 11.6 kPa, so the cap has next to nothing to
  trim and the strength reduction answers the same.
* **The file declares the finite element run it is solved on** — tri6 at a 0.9 m
  target size, the strength-reduction bracket 0.8 to 1.6, and K0 = 1 — so
  Studio's Build Mesh and Run FEM dialogs open on the tutorial's own run instead
  of on the defaults. Those are the settings the SSRM lock tag names. The file
  carries no mesh sidecar: the page builds the mesh in Studio.

Run:  PYTHONPATH=. python3 tools/build_rock_slope.py
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
HAMMAH_OUT = "xslope_rock_slope.xlsx"

#: Part A's finite element run, from the SSRM lock tag at
#: docs/verification/rs2.md (element_type=tri6, target_size=0.9, f_min=0.8,
#: f_max=1.6, k0=1). tri6 because a factor of safety needs quadratic elements.
ELEMENT_TYPE = "tri6"
TARGET_SIZE = 0.9
SSRM_F_MIN, SSRM_F_MAX = 0.8, 1.6
K0 = 1.0

#: The tensile cutoff every shipped tutorial file declares. The rock's own
#: Hoek-Brown tensile strength is s*sigma_ci/mb = 11.6 kPa, so a zero cap trims
#: little; declaring it keeps Run FEM's unbounded-tension check quiet on a page
#: whose runs are meant to open clean.
TENSILE_CUTOFF = 0.0


def _base(path):
    """The corpus model with the two changes the tutorial file makes: no saturated
    unit weight on a section that has no water in it, and a declared tensile
    cutoff."""
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
    mat["t_cut"] = TENSILE_CUTOFF
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
    """The weak-rock slope, carrying the finite element run as well as the limit
    equilibrium one, since the page is a two-engine comparison."""
    sd = _base(HAMMAH_BASE)
    sd["element_type"] = ELEMENT_TYPE
    sd["target_size"] = TARGET_SIZE
    sd["ssrm_f_min"] = SSRM_F_MIN
    sd["ssrm_f_max"] = SSRM_F_MAX
    sd["k0"] = K0
    return _write(sd, HAMMAH_OUT)


if __name__ == "__main__":
    build_hammah()
