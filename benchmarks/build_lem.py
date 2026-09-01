"""Build the LEM benchmark input files (all metric: m, kN/m3, kPa).

Geometry and properties pulled from the GeoStudio SLOPE/W Verification Manual
(Oct 2022): ACADS Simple Slope (sec 2.1), ACADS Weak Layer (sec 2.7), and the
Arai & Tagyo Homogeneous Slope (sec 2.11). Coordinates read directly from the
manual's dimensioned geometry figures.

Run from the repo root:  PYTHONPATH=. python3 benchmarks/build_lem.py
"""

from benchmarks._xlsx_writer import (
    new_file, write_cells_to_xlsx, cell_ref,
    main_cells, material_cells, profile_line_cells, circle_cells,
    noncirc_cells, piezo_cells,
)
from benchmarks.rocscience.vendor_tcut import VENDOR_T_CUT, VENDOR_E_NU
from benchmarks.tag_k0 import tag_k0

TEMPLATE = "docs/inputs/input_template.xlsx"
OUTDIR = "docs/lem/files"

# These three files double as RS2 SSRM benchmarks (RS2-1 / RS2-5 / RS2-10 on
# docs/verification/rs2.md), so they carry the vendor .fez tensile caps. t_cut is a
# FEM-only field — the LEM solvers never read it — so the LEM locks here are
# untouched. The values live with every other RS2 cap in vendor_tcut.VENDOR_T_CUT.
_T = VENDOR_T_CUT

# All three are also SSRM rows, and RS2 authors its models at an isotropic at-rest
# field stress (Kx = Kz = 1), so each file DECLARES K0 = 1 in main!D16 rather than
# leaving the convention on the doc page's test tag. tag_k0 reads the value off that
# tag (benchmarks/tag_k0.py) so the cell and the tag cannot drift apart.


def _vendor_e_nu(dst_name, mat_name):
    """(E, nu) transcribed from the vendor .fez these three files reproduce.

    They are RS2-1 / RS2-5 / RS2-10, so RS2's own model specifies their elastic
    constants and the published SSR was solved with them — a specified input, which
    transcribes verbatim rather than being replaced by a soil-type classifier value
    (see vendor_tcut.VENDOR_E_NU). Written here explicitly because these builders
    write cells directly and never pass through the RS2 save wrapper.
    """
    nu, E = VENDOR_E_NU[dst_name][mat_name]
    return E, nu


# ---------------------------------------------------------------------------
# LEM-1 : ACADS Simple Slope (SLOPE/W 2.1) — circular, single soil
#   Geometry: bench (0,25)-(30,25); 2:1 face (30,25)->(50,35); crest
#   (50,35)-(90,35); base elevation 20.  10 m high.
#   The flats run out to x = 0 and x = 90 rather than stopping at the
#   dimensioned 20 / 70: the domain must extend well past the slope face on
#   both sides so a searched circle never runs off the end of the ground
#   surface. The face itself is unchanged, so the ACADS answer is unchanged.
#   Soil c=3.0, phi=19.6, gamma=20.0 (total stress).  ACADS FOS = 1.00.
#   This file is also the reliability/sensitivity sample (docs/lem/design.md,
#   docs/parametric/*), which reads the strength standard deviations below.
# ---------------------------------------------------------------------------
def build_acads_simple():
    dst = f"{OUTDIR}/xslope_acads_simple.xlsx"
    new_file(dst, TEMPLATE)
    u = {}
    _E_AS = _vendor_e_nu("xslope_acads_simple.xlsx", "Soil")
    u['main'] = main_cells(gamma_w=9.81, template=TEMPLATE,
                           k0=tag_k0("xslope_acads_simple.xlsx"))
    u['mat'] = material_cells(1, "Soil", 20.0, "mc", c=3.0, phi=19.6, u="none", template=TEMPLATE,
                              t_cut=_T["xslope_acads_simple.xlsx"]["Soil"],
                              E=_E_AS[0], nu=_E_AS[1],
                              # Std devs published with the design sample
                              # (docs/lem/design.md): the reliability locks read them.
                              sigma_gamma=1.2, sigma_c=1.8, sigma_phi=2.744)
    prof = {'B2': 20}  # base elevation
    prof.update(profile_line_cells(1, 1, [(0, 25), (30, 25), (50, 35), (90, 35)]))
    u['profile'] = prof
    circ = {}
    circ.update(circle_cells(1, 40, 45, option="Intercept", xi=30, yi=25))  # toe circle
    circ.update(circle_cells(2, 42, 48, option="Depth", depth=20))          # base-tangent
    u['circles'] = circ
    write_cells_to_xlsx(dst, {k: v for k, v in u.items() if v})
    print("built", dst)
    return dst


# ---------------------------------------------------------------------------
# LEM-2 : ACADS Weak Layer (SLOPE/W 2.7) — non-circular, weak interlayer
#   Ground (20,27.75)-(43,27.75)-(67.5,40)-(84,40); 2:1 face.
#   Weak band: top y=27, base y=26.5 (Soil 2) spanning x=20..84.
#   Soil 1 above and below the band; base y=20; piezo at y=26.5.
#   Soil1 c=28.5 phi=20 g=18.84 ; Soil2 c=0 phi=10 g=18.84.  ACADS FOS = 1.26.
# ---------------------------------------------------------------------------
def build_acads_weak_layer():
    dst = f"{OUTDIR}/xslope_acads_weak_layer.xlsx"
    new_file(dst, TEMPLATE)
    u = {}
    u['main'] = main_cells(gamma_w=9.81, template=TEMPLATE,
                           k0=tag_k0("xslope_acads_weak_layer.xlsx"))
    mat = {}
    _t = _T["xslope_acads_weak_layer.xlsx"]
    _e1 = _vendor_e_nu("xslope_acads_weak_layer.xlsx", "Soil 1")
    _e2 = _vendor_e_nu("xslope_acads_weak_layer.xlsx", "Weak Layer")
    mat.update(material_cells(1, "Soil 1", 18.84, "mc", c=28.5, phi=20.0, u="piezo", template=TEMPLATE,
                              t_cut=_t["Soil 1"],
                              E=_e1[0], nu=_e1[1]))
    mat.update(material_cells(2, "Weak Layer", 18.84, "mc", c=0.0, phi=10.0, u="piezo", template=TEMPLATE,
                              t_cut=_t["Weak Layer"],
                              E=_e2[0], nu=_e2[1]))
    u['mat'] = mat
    prof = {'B2': 20}
    prof.update(profile_line_cells(1, 1, [(20, 27.75), (43, 27.75), (67.5, 40), (84, 40)]))
    prof.update(profile_line_cells(2, 2, [(20, 27.0), (84, 27.0)]))   # top of weak layer
    prof.update(profile_line_cells(3, 1, [(20, 26.5), (84, 26.5)]))   # base of weak layer
    u['profile'] = prof
    u['piezo'] = piezo_cells([(20, 26.5), (84, 26.5)])  # at base of weak layer
    # Non-circular block: exit at toe bench, slide along weak layer (y~26.55), back scarp to crest.
    # The 'Horiz' points are seeded just above the base of the weak layer (base y=26.5,
    # top y=27.0): the search only moves them horizontally, so this seed elevation is the
    # sliding plane. Placing the slip near the base of a weak layer is standard practice and
    # matches the ACADS/SLOPE/W reference (FOS ~1.26); the layer-center seed reads ~1.5% high.
    u['non-circ'] = noncirc_cells([
        (38, 27.75, "Free"),    # exit on bench near toe
        (44, 26.55, "Horiz"),   # toe-side, near base of weak layer
        (66, 26.55, "Horiz"),   # crest-side, near base of weak layer
        (75, 40.0, "Free"),     # entry on crest
    ])
    write_cells_to_xlsx(dst, {k: v for k, v in u.items() if v})
    print("built", dst)
    return dst


# ---------------------------------------------------------------------------
# LEM-2b : Arai & Tagyo Homogeneous Slope (SLOPE/W 2.11) — circular, single soil
#   Ground (-30,15)-(18,15)-(48,35)-(100,35); 1.5:1 face; base y=0.  20 m high.
#   As in LEM-1 the flats are run out well past the dimensioned 0 / 66 so a
#   searched circle cannot leave the ground surface; the face is unchanged.
#   Soil c=41.65, phi=15.0, gamma=18.82 (total stress).  Published FOS = 1.451.
# ---------------------------------------------------------------------------
def build_arai_tagyo():
    dst = f"{OUTDIR}/xslope_arai_tagyo.xlsx"
    new_file(dst, TEMPLATE)
    u = {}
    _E_AT = _vendor_e_nu("xslope_arai_tagyo.xlsx", "Soil")
    u['main'] = main_cells(gamma_w=9.81, template=TEMPLATE,
                           k0=tag_k0("xslope_arai_tagyo.xlsx"))
    u['mat'] = material_cells(1, "Soil", 18.82, "mc", c=41.65, phi=15.0, u="none", template=TEMPLATE,
                              t_cut=_T["xslope_arai_tagyo.xlsx"]["Soil"],
                              E=_E_AT[0], nu=_E_AT[1])
    prof = {'B2': 0}
    prof.update(profile_line_cells(1, 1, [(-30, 15), (18, 15), (48, 35), (100, 35)]))
    u['profile'] = prof
    circ = {}
    circ.update(circle_cells(1, 33, 55, option="Intercept", xi=18, yi=15))  # toe circle
    circ.update(circle_cells(2, 35, 55, option="Depth", depth=0))           # base-tangent
    u['circles'] = circ
    write_cells_to_xlsx(dst, {k: v for k, v in u.items() if v})
    print("built", dst)
    return dst


if __name__ == "__main__":
    build_acads_simple()
    build_acads_weak_layer()
    build_arai_tagyo()
    print("\nLEM benchmark input files built.")
