"""Build the FE-seepage benchmark input files.

SEEP-1 — Kozeny / Casagrande analytical anchor: a homogeneous earth dam on an
impervious base with a horizontal toe drain. The phreatic surface is a parabola
with focus at the upstream end of the drain, and the discharge has the closed
form  q = k * y0,  y0 = sqrt(d^2 + h^2) - d, where h is the reservoir depth and
d is the horizontal distance from the focus to the point where the reservoir
water surface meets the upstream face.

Geometry (units arbitrary; k carries them):
  outline  (0,0)-(60,24)-(70,24)-(118,0)   upstream 2.5:1, crest 10 wide, ds 2:1
  base     impervious at y=0
  reservoir depth h=20 -> water meets upstream face at B=(50,20)
  toe drain (exit face) (98,0)-(118,0)      focus F=(98,0)
  => d=48, h=20, y0 = sqrt(48^2+20^2)-48 = 52-48 = 4.0,  q = 4.0*k

Run from the repo root:  PYTHONPATH=. python3 benchmarks/build_seep.py
"""
from benchmarks._xlsx_writer import (
    new_file, write_cells_to_xlsx,
    main_cells, material_cells, polygon_cells, profile_line_cells, seep_bc_cells,
)

TEMPLATE = "docs/inputs/input_template.xlsx"
OUTDIR = "docs/seep/files"

# --- Kozeny exact-solution geometry -----------------------------------------
# focus (drain tip) at x=XF, base y=0; phreatic focal distance y0; reservoir
# depth h. The EXACT upstream equipotential (head=h) is the confocal parabola
#   x = XF - (d - y^2/(2*y0)) ... derived below as x = (XF-d) + y^2/(2*b)
# Choose d=48, h=20  ->  y0 = sqrt(d^2+h^2)-d = 52-48 = 4,  q = k*y0 = 4.
XF = 88.0          # focus (upstream tip of horizontal toe drain)
H = 20.0           # reservoir depth
D = 48.0           # focus-to-entry horizontal distance -> entry A at x=XF-D=40
Y0 = (D**2 + H**2) ** 0.5 - D          # = 4.0
KOZENY = dict(k=1.0, h=H, d=D, y0=Y0, q=1.0 * Y0)

# Upstream equipotential (head=H) is the parabola confocal with the phreatic
# parabola, passing through the entry A=(XF-D, H) and the base. Solving the
# orthogonal confocal family gives  x_model = (XF - B) + y^2/(2B)  with B chosen
# so the curve passes through A; for d=48,h=20 this is B=100, i.e.
#   x = 38 + y^2/200,  from (38,0) at the base to (40,20) at the water surface.
def _upstream_parabola(n=11):
    pts = []
    for i in range(n):
        y = H * i / (n - 1)
        x = (XF - 50.0) + y * y / 200.0   # = 38 + y^2/200
        pts.append((round(x, 4), round(y, 4)))
    return pts


def build_kozeny_dam():
    """Homogeneous dam, horizontal toe drain, impervious base. The upstream face
    is built (via the polygon input) as the EXACT confocal-parabola equipotential
    of the Kozeny flow net, so the FE discharge matches q=k*y0 to within
    discretization error instead of the ~7% basic-parabola entrance error that a
    straight 2.5:1 face produces."""
    dst = f"{OUTDIR}/xslope_kozeny_dam.xlsx"
    new_file(dst, TEMPLATE)
    up = _upstream_parabola()             # (38,0) ... (40,20)
    # The free surface terminates at the parabola VERTEX (downstream of the
    # focus by y0/2). The drain (head=0 exit face) must start there; the base
    # upstream of the vertex is impervious (saturated, positive pressure).
    xv = XF + Y0 / 2.0                     # vertex x = 90 (free surface meets base)
    # closed boundary ring: parabolic face, freeboard, crest, ds 2:1 face, toe,
    # base back to heel.
    ring = up + [(40.0, 24.0), (50.0, 24.0), (98.0, 0.0), up[0]]
    u = {}
    u['main'] = main_cells(gamma_w=9.81)
    u['mat'] = material_cells(9, 1, "Dam", 20.0, "mc", c=0.0, phi=30.0, u="seep",
                              k1=1.0, k2=1.0, alpha=0.0, kr0=0.01, h0=-1.0)
    u['polygon'] = polygon_cells(1, 1, ring)
    u['profile'] = {'B2': 0}              # max depth (impervious base) = 0
    u['seep bc'] = seep_bc_cells(
        exit_face=[(xv, 0.0), (98.0, 0.0)],     # toe drain, from free-surface exit
        head1=H, head1_pts=up,                  # reservoir on parabolic face
    )
    write_cells_to_xlsx(dst, {k: v for k, v in u.items() if v})
    print("built", dst, "| y0 =", Y0, "| q_analytical =", KOZENY['q'])
    return dst


# --- SEEP-1 : confined radial flow (saturated analytical anchor) -------------
# Quarter annulus, inner radius R1 at head H1, outer radius R2 at head H2; the
# two straight radial edges are no-flow streamlines. Steady confined (saturated)
# flow has the EXACT closed form (Laplace in polar, r-only):
#   head:  h(r) = H1 + (H2-H1)*ln(r/R1)/ln(R2/R1)
#   discharge (sector angle alpha):  q = k*alpha*(H1-H2)/ln(R2/R1)
# No free surface, no exit face -> a clean confined verification with no
# tangency sensitivity. Geometry/formula are exact; only the arc discretization
# introduces (small) error.
import math

R1, R2 = 10.0, 30.0
H1, H2 = 30.0, 10.0
ALPHA = math.pi / 2.0                  # quarter sector
CONFINED = dict(k=1.0, r1=R1, r2=R2, h1=H1, h2=H2, alpha=ALPHA,
                q=1.0 * ALPHA * (H1 - H2) / math.log(R2 / R1))


def _arc(radius, n=36):
    """Quarter-circle arc points from angle 0 to pi/2 (CCW)."""
    return [(round(radius * math.cos(ALPHA * i / n), 5),
             round(radius * math.sin(ALPHA * i / n), 5)) for i in range(n + 1)]


def build_confined_radial():
    dst = f"{OUTDIR}/xslope_confined_radial.xlsx"
    new_file(dst, TEMPLATE)
    inner = _arc(R1)                    # (10,0) -> (0,10)
    outer = _arc(R2)                    # (30,0) -> (0,30)
    # CCW ring: x-axis out, outer arc, y-axis in, inner arc back.
    ring = ([(R1, 0.0)] + [(R2, 0.0)] + outer[1:] +
            [(0.0, R1)] + list(reversed(inner))[1:])
    u = {}
    u['main'] = main_cells(gamma_w=9.81)
    u['mat'] = material_cells(9, 1, "Soil", 20.0, "mc", c=0.0, phi=30.0, u="seep",
                              k1=1.0, k2=1.0, alpha=0.0)
    u['polygon'] = polygon_cells(1, 1, ring)
    u['profile'] = {'B2': 0}
    u['seep bc'] = seep_bc_cells(
        head1=H1, head1_pts=inner,         # inner arc (high head)
        head2=H2, head2_pts=outer,         # outer arc (low head)
    )
    write_cells_to_xlsx(dst, {k: v for k, v in u.items() if v})
    print("built", dst, "| q_analytical =", round(CONFINED['q'], 4))
    return dst


# --- SEEP-1c : partially penetrating sheetpile (confined, exact) -------------
# Pavlovsky's conformal-mapping solution (Harr 1962; Polubarinova-Kochina 1962)
# for a single cutoff wall of depth s in a homogeneous confined stratum of
# thickness T with head loss H across the wall:
#   q = k*H*K(lam')/(2*K(lam)),  lam = sin(pi*s/(2T)),  lam' = cos(pi*s/(2T))
# (K = complete elliptic integral of the first kind). At s/T = 1/2 the modulus
# is self-dual (lam = lam') so q = k*H/2 EXACTLY. A second exact check: by
# antisymmetry the head on the wall plane below the tip is exactly (H1+H2)/2.
T_SP = 20.0      # stratum thickness
L_SP = 80.0      # horizontal extent each side of the wall (4T: truncation
                 # error ~ exp(-pi*L/T) ~ 3e-6)
H1_SP, H2_SP = 25.0, 15.0
W_SP = 0.1       # notch width at the surface (V-notch crack idiom, w/T = 0.005)


def sheetpile_q_exact(s_over_T, k=1.0, H=H1_SP - H2_SP):
    from scipy.special import ellipk   # takes the PARAMETER m = modulus^2
    import math
    lam2 = math.sin(math.pi * s_over_T / 2.0) ** 2
    return k * H * ellipk(1.0 - lam2) / (2.0 * ellipk(lam2))


def build_sheetpile(s_over_T=0.5):
    sfx = str(int(round(100 * s_over_T)))
    dst = f"{OUTDIR}/xslope_sheetpile_s{sfx}.xlsx"
    new_file(dst, TEMPLATE)
    s = s_over_T * T_SP
    u = {}
    u['main'] = main_cells(gamma_w=9.81)
    u['mat'] = material_cells(9, 1, "Sand", 20.0, "mc", c=0.0, phi=30.0, u="seep",
                              k1=1.0, k2=1.0, alpha=0.0, kr0=0.01, h0=-1.0)
    # ground surface with V-notch for the wall (same idiom as the clay-blanket
    # sample: profile dips to the wall tip and back up, leaving a crack)
    prof = {'B2': 0}
    prof.update(profile_line_cells(1, 1, [
        (-L_SP, T_SP), (-W_SP / 2, T_SP), (0.0, T_SP - s), (W_SP / 2, T_SP),
        (L_SP, T_SP)]))
    u['profile'] = prof
    u['seep bc'] = seep_bc_cells(
        head1=H1_SP, head1_pts=[(-L_SP, T_SP), (-W_SP / 2, T_SP)],
        head2=H2_SP, head2_pts=[(W_SP / 2, T_SP), (L_SP, T_SP)],
    )
    write_cells_to_xlsx(dst, {k: v for k, v in u.items() if v})
    print("built", dst, f"| s/T={s_over_T} | q_exact={sheetpile_q_exact(s_over_T):.4f}")
    return dst


if __name__ == "__main__":
    build_kozeny_dam()
    build_confined_radial()
    build_sheetpile(0.5)
    build_sheetpile(0.75)
    print("\nSEEP benchmark input files built.")
