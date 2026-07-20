"""Build the FE-seepage benchmark input files.

SEEP-1 anchors are exact confined solutions (no free-surface tangency issues):
confined radial flow (quarter annulus) and Pavlovsky's partially penetrating
sheetpile. The unconfined solver is verified code-to-code against SEEP2D on
the Johnson Reservoir problem (see benchmarks/run_seep2d_compare.py).

Run from the repo root:  PYTHONPATH=. python3 benchmarks/build_seep.py
"""
from benchmarks._xlsx_writer import (
    new_file, write_cells_to_xlsx,
    main_cells, material_cells, polygon_cells, profile_line_cells, seep_bc_cells,
)

TEMPLATE = "docs/inputs/input_template.xlsx"
OUTDIR = "docs/seep/files"

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
    u['mat'] = material_cells(1, "Soil", 20.0, "mc", c=0.0, phi=30.0, u="seep",
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
# Heads chosen so the section is physically consistent as a vertical
# cross-section: downstream head equals the stratum top (zero ponding, p >= 0
# everywhere); upstream ponded 10 above it. Only the difference H = 10 enters
# the confined solution.
H1_SP, H2_SP = 30.0, 20.0
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
    u['mat'] = material_cells(1, "Sand", 20.0, "mc", c=0.0, phi=30.0, u="seep",
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
    build_confined_radial()
    build_sheetpile(0.5)
    build_sheetpile(0.75)
    print("\nSEEP benchmark input files built.")
