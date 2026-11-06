"""Regenerate the two figures of "Slices on a Circle" in docs/lem/overview.md.

circle_slices_width_angle.png
    One circle meeting the crest of a 20 ft, 2V:1H slope vertically, cut into
    twelve slices two ways: equal widths, and equal angles. The last slice is
    shaded, with the arc it spans.

circle_slice_base.png
    One steep slice on a circle, drawn large: the arc between its corners, the
    chord, the middle of the chord where the base forces act, the middle of the
    arc where the pore pressure is read, and the point on the circle directly
    below the slice's center line.

    PYTHONPATH=. python3 docs/lem/images/circle_slicing.py
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Arc, FancyArrowPatch, Polygon

HERE = os.path.dirname(os.path.abspath(__file__))
SOIL = "#f7c4a0"
SHADE = "#e8735a"
WATER = "#1f6fbf"
DPI = 200

# The slope and circle of test/circle_arc_check.py (STEEP): tangent to the firm
# layer at the toe's elevation, meeting the crest vertically.
GROUND = [(-8.0, 0.0), (0.0, 0.0), (10.0, 20.0), (37.0, 20.0)]
XO, YO, R = 2.7317, 20.0049, 20.0049
N = 12


def ground_y(x):
    gx, gy = zip(*GROUND)
    return np.interp(x, gx, gy)


def arc_y(x):
    return YO - np.sqrt(np.maximum(R * R - (x - XO) ** 2, 0.0))


def theta(x):
    return np.degrees(np.arcsin(np.clip((x - XO) / R, -1.0, 1.0)))


def entry_exit():
    """Where the circle meets the ground: the face near the toe, and the crest."""
    xs = np.linspace(0.0, 10.0, 200001)
    d = ground_y(xs) - arc_y(xs)
    x_in = xs[np.argmax(d > 0)]
    x_out = XO + np.sqrt(R * R - (YO - 20.0) ** 2)
    return float(x_in), float(x_out)


def draw_section(ax, bounds, title):
    x0, x1 = bounds[0], bounds[-1]
    soil = list(GROUND) + [(GROUND[-1][0], -3.0), (GROUND[0][0], -3.0)]
    ax.add_patch(Polygon(soil, closed=True, fc=SOIL, ec="none", zorder=0))
    gx, gy = zip(*GROUND)
    ax.plot(gx, gy, color="black", lw=1.2, zorder=3)
    # the last slice, shaded
    xs = np.linspace(bounds[-2], bounds[-1], 200)
    ax.fill_between(xs, arc_y(xs), ground_y(xs), color=SHADE, alpha=0.55, lw=0, zorder=1)
    for x in bounds[1:-1]:
        ax.plot([x, x], [arc_y(x), ground_y(x)], ls="--", color="black", lw=0.9, zorder=2)
    xs = np.linspace(x0, x1, 600)
    ax.plot(xs, arc_y(xs), color="black", lw=2.4, zorder=4)
    # the arc the last slice spans, in the shade's own color
    t = np.radians(np.linspace(theta(bounds[-2]), theta(bounds[-1]), 200))
    ax.plot(XO + R * np.sin(t), YO - R * np.cos(t), color="#c0392b", lw=4.2,
            solid_capstyle="butt", zorder=5)
    tm = np.radians(0.5 * (theta(bounds[-2]) + theta(bounds[-1])))
    span = theta(bounds[-1]) - theta(bounds[-2])
    ax.annotate("last slice: %.0f° of arc" % span,
                xy=(XO + R * np.sin(tm) + 0.3, YO - R * np.cos(tm)), xytext=(25.0, 9.0),
                fontsize=12.5, ha="left", va="center", color="#c0392b",
                arrowprops=dict(arrowstyle="-", color="#c0392b", lw=0.9,
                                shrinkA=3, shrinkB=9))
    ax.text(-7.0, 22.6, title, fontsize=13.5, ha="left", va="bottom", weight="bold")
    ax.set_xlim(-8.0, 37.0)
    ax.set_ylim(-3.0, 24.5)
    ax.set_aspect("equal")
    ax.axis("off")


def split(fixed, length, cut):
    """Slice boundaries the way the slicer places them: each stretch between fixed
    points takes its share of the N slices by ``length`` and is cut by ``cut``."""
    lens = [length(a, b) for a, b in zip(fixed[:-1], fixed[1:])]
    out = [fixed[0]]
    for (a, b), L in zip(zip(fixed[:-1], fixed[1:]), lens):
        n = max(1, int(round(L / sum(lens) * N)))
        out += list(cut(a, b, n)[1:])
    return np.array(out)


def figure_width_angle():
    x_in, x_out = entry_exit()
    # the crest is a ground vertex above the surface, so it is a slice boundary
    fixed = [x_in, GROUND[2][0], x_out]
    widths = split(fixed, lambda a, b: b - a, lambda a, b, n: np.linspace(a, b, n + 1))

    def by_angle(a, b, n):
        xs = XO + R * np.sin(np.radians(np.linspace(theta(a), theta(b), n + 1)))
        xs[0], xs[-1] = a, b
        return xs
    angles = split(fixed, lambda a, b: theta(b) - theta(a), by_angle)
    for name, xs in (("widths", widths), ("angles", angles)):
        print("   %s: %d slices, arcs %s deg" % (name, len(xs) - 1,
              ", ".join("%.1f" % (theta(b) - theta(a)) for a, b in zip(xs[:-1], xs[1:]))))
    fig, axes = plt.subplots(2, 1, figsize=(7.4, 8.6))
    draw_section(axes[0], widths, "Equal widths")
    draw_section(axes[1], angles, "Equal angles")
    fig.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.01, hspace=0.04)
    out = os.path.join(HERE, "circle_slices_width_angle.png")
    fig.savefig(out, dpi=DPI, bbox_inches="tight", pad_inches=0.05)
    plt.close(fig)
    print("wrote", out)


def figure_base():
    """One slice from 42 to 82 degrees of a circle of radius 10, centered at O,
    drawn zoomed: the radius through M runs off the panel toward O. C, M, P and the arc length
    are labeled beside the slice with pointer lines that stop short of what they
    point at; the corner letters, the angle and the two forces are labeled in place."""
    r = 10.0
    t1, t2 = np.radians(42.0), np.radians(82.0)
    tm = 0.5 * (t1 + t2)
    on = lambda t: np.array([r * np.sin(t), -r * np.cos(t)])
    A, B, M = on(t1), on(t2), on(tm)
    C = 0.5 * (A + B)
    tp = np.arcsin(C[0] / r)
    P = on(tp)
    top = 2.6
    x_lo, x_hi, y_lo, y_hi = 1.4, 17.4, -8.9, 3.4
    right = 11.9                    # where the labels on the right begin
    units_per_inch = 1.75
    fig, ax = plt.subplots(figsize=((x_hi - x_lo) / units_per_inch,
                                    (y_hi - y_lo) / units_per_inch))

    ts = np.linspace(t1, t2, 200)
    arc = np.array([on(t) for t in ts])
    poly = [(A[0], top)] + [tuple(p) for p in arc] + [(B[0], top)]
    ax.add_patch(Polygon(poly, closed=True, fc=SOIL, ec="none", zorder=0))
    ax.plot([A[0] - 1.0, B[0] + 1.0], [top, top], color="black", lw=1.2)
    for x, y in ((A[0], A[1]), (B[0], B[1])):
        ax.plot([x, x], [y, top], ls="--", color="black", lw=0.9)
    tt = np.linspace(np.radians(36), np.radians(88), 300)
    ax.plot(r * np.sin(tt), -r * np.cos(tt), color="0.7", lw=1.0, zorder=1)

    # the radius through M, toward the center O, clipped at the panel's left edge
    def toward_o(p, color, lw):
        q = p * (x_lo + 0.1) / p[0]
        ax.plot([p[0], q[0]], [p[1], q[1]], color=color, lw=lw, zorder=1)
        return q
    q = toward_o(M, "0.55", 0.9)
    # along the line, just above it, near its far end
    u = (M - q) / np.linalg.norm(M - q)
    nrm = np.array([-u[1], u[0]]) if u[0] > 0 else np.array([u[1], -u[0]])
    at = q + 0.25 * u + 0.16 * nrm
    ax.text(at[0], at[1], "toward O", fontsize=11, color="0.4", ha="left", va="bottom",
            rotation=np.degrees(np.arctan2(u[1], u[0])), rotation_mode="anchor")

    # the base: arc and chord
    ax.plot(arc[:, 0], arc[:, 1], color="black", lw=2.6, zorder=4)
    ax.plot([A[0], B[0]], [A[1], B[1]], color="black", lw=1.0, ls=(0, (6, 3)), zorder=4)
    # the slice's center line and the weight on it
    ax.plot([C[0], C[0]], [P[1] - 0.6, top + 0.3], ls=":", color="0.2", lw=1.1)
    mid = 0.5 * (top + C[1])         # halfway down the slice, on its center line
    ax.add_patch(FancyArrowPatch((C[0], mid + 1.0), (C[0], mid - 1.0),
                                 arrowstyle="-|>", mutation_scale=14, color="black", lw=1.4))
    ax.text(C[0] + 0.22, mid, "W", fontsize=13.5, ha="left", va="center")
    # the base normal at C, along the radius toward O, labeled below the arrow
    d = -C / np.linalg.norm(C)
    up = np.array([-d[1], d[0]])
    if up[1] < 0:
        up = -up
    n_len = 1.25
    tip = C + n_len * d
    ax.add_patch(FancyArrowPatch(tuple(C), tuple(tip), arrowstyle="-|>", shrinkA=5,
                                 mutation_scale=14, color="black", lw=1.4, zorder=6))
    npos = C + 0.55 * n_len * d + 0.48 * up
    ax.text(npos[0], npos[1], "N", fontsize=13.5, ha="center", va="center")
    # alpha, the chord's inclination, drawn at A
    ax.plot([A[0], A[0] + 2.0], [A[1], A[1]], color="0.3", lw=0.9)
    ax.add_patch(Arc(tuple(A), 2.8, 2.8, theta1=0, theta2=np.degrees(tm), color="0.3",
                     lw=0.9))
    ax.text(A[0] + 1.62, A[1] + 0.55, "α", fontsize=14, ha="left", va="center")
    # the tangent at P, flatter than the chord
    tan = np.array([np.cos(tp), np.sin(tp)])
    ax.plot([P[0] - 2.4 * tan[0], P[0] + 2.4 * tan[0]],
            [P[1] - 2.4 * tan[1], P[1] + 2.4 * tan[1]], color="0.35", lw=1.3, ls="--",
            zorder=5)
    # the three points
    ax.plot(*C, "o", ms=8, mfc="white", mec="black", mew=1.8, zorder=7)
    ax.plot(*M, "o", ms=8, color=WATER, zorder=7)
    ax.plot(*P, "o", ms=7, mfc="white", mec="0.4", mew=1.5, zorder=7)
    ax.text(A[0] - 0.47, A[1] + 0.03, "A", fontsize=13.5, ha="center", va="center")
    ax.text(B[0] + 0.32, B[1] + 0.3, "B", fontsize=13.5, ha="left", va="center")

    def label(target, xy_text, text, color, ha="left"):
        ax.annotate(text, xy=tuple(target), xytext=xy_text, fontsize=12, ha=ha,
                    va="center", color=color, linespacing=1.2,
                    arrowprops=dict(arrowstyle="-", color=color, lw=0.9,
                                    shrinkA=4, shrinkB=8))

    label(C, (A[0] - 0.55, C[1] - 0.6), "C, the middle of the chord AB:\n"
          "the base forces act here", "black", ha="right")
    label(on(t1 + 0.8 * (t2 - t1)), (right, -1.6), "Δℓ, the length of the arc AB",
          "black")
    label(M, (right, -3.9), "M, the middle of the arc AB:\nthe pore pressure is read here",
          WATER)
    label(P, (right, -6.3), "P, the point on the circle directly\nbelow the slice's center line\n"
          "(dashed line: the circle's tangent at P)", "0.35")
    ax.set_xlim(x_lo, x_hi)
    ax.set_ylim(y_lo, y_hi)
    ax.set_aspect("equal")
    ax.axis("off")
    out = os.path.join(HERE, "circle_slice_base.png")
    fig.savefig(out, dpi=DPI, bbox_inches="tight", pad_inches=0.05)
    plt.close(fig)
    print("wrote", out)


if __name__ == "__main__":
    figure_width_angle()
    figure_base()
