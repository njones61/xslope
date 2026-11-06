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
                xy=(XO + R * np.sin(tm) + 0.3, YO - R * np.cos(tm)), xytext=(26.0, 9.0),
                fontsize=10.5, ha="left", va="center", color="#c0392b",
                arrowprops=dict(arrowstyle="-", color="#c0392b", lw=0.9))
    ax.text(-7.0, 22.6, title, fontsize=12, ha="left", va="bottom", weight="bold")
    ax.set_xlim(-8.0, 37.0)
    ax.set_ylim(-3.0, 24.5)
    ax.set_aspect("equal")
    ax.axis("off")


def figure_width_angle():
    x_in, x_out = entry_exit()
    widths = np.linspace(x_in, x_out, N + 1)
    t = np.radians(np.linspace(theta(x_in), theta(x_out), N + 1))
    angles = XO + R * np.sin(t)
    angles[0], angles[-1] = x_in, x_out
    fig, axes = plt.subplots(2, 1, figsize=(8.6, 9.6))
    draw_section(axes[0], widths, "Equal widths")
    draw_section(axes[1], angles, "Equal angles")
    fig.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.01, hspace=0.04)
    out = os.path.join(HERE, "circle_slices_width_angle.png")
    fig.savefig(out, dpi=DPI, bbox_inches="tight", pad_inches=0.05)
    plt.close(fig)
    print("wrote", out)


def figure_base():
    """One slice from 45 to 80 degrees of a circle of radius 10, centered at O."""
    r = 10.0
    t1, t2 = np.radians(45.0), np.radians(80.0)
    tm = 0.5 * (t1 + t2)
    on = lambda t: np.array([r * np.sin(t), -r * np.cos(t)])
    A, B, M = on(t1), on(t2), on(tm)
    C = 0.5 * (A + B)
    tp = np.arcsin(C[0] / r)
    P = on(tp)
    top = 3.0
    lab = 12.2                      # the column the labels stand in

    fig, ax = plt.subplots(figsize=(8.4, 6.6))
    ts = np.linspace(t1, t2, 200)
    arc = np.array([on(t) for t in ts])
    poly = [(A[0], top)] + [tuple(p) for p in arc] + [(B[0], top)]
    ax.add_patch(Polygon(poly, closed=True, fc=SOIL, ec="none", zorder=0))
    ax.plot([A[0] - 1.2, B[0] + 1.2], [top, top], color="black", lw=1.2)
    for x, y in ((A[0], A[1]), (B[0], B[1])):
        ax.plot([x, x], [y, top], ls="--", color="black", lw=0.9)
    # the circle beyond the slice, faint, and the radii to A, M and B
    tt = np.linspace(np.radians(25), np.radians(88), 300)
    ax.plot(r * np.sin(tt), -r * np.cos(tt), color="0.7", lw=1.0, zorder=1)
    for p in (A, M, B):
        ax.plot([0, p[0]], [0, p[1]], color="0.7", lw=0.8, zorder=1)
    ax.plot([0], [0], "o", color="0.35", ms=5)
    ax.text(-0.3, 0.3, "O", fontsize=12, ha="right", va="bottom", color="0.3")
    # the base: arc and chord
    ax.plot(arc[:, 0], arc[:, 1], color="black", lw=2.6, zorder=4)
    ax.plot([A[0], B[0]], [A[1], B[1]], color="black", lw=1.0, ls=(0, (6, 3)), zorder=4)
    # the slice's center line, where the weight acts
    ax.plot([C[0], C[0]], [P[1] - 0.5, top + 0.25], ls=":", color="0.2", lw=1.1)
    ax.add_patch(FancyArrowPatch((C[0], top - 0.3), (C[0], top - 2.4),
                                 arrowstyle="-|>", mutation_scale=14, color="black", lw=1.4))
    ax.text(C[0] + 0.18, top - 1.3, "W", fontsize=12, ha="left", va="center")
    # the base normal at C, along the radius toward O
    d = -C / np.linalg.norm(C)
    ax.add_patch(FancyArrowPatch(tuple(C), tuple(C + 1.8 * d), arrowstyle="-|>",
                                 mutation_scale=14, color="black", lw=1.4, zorder=6))
    ax.text(*(C + 1.8 * d + np.array([0.0, 0.42])), "N", fontsize=12, ha="center",
            va="center")
    # alpha, the chord's inclination, drawn at A
    ax.plot([A[0], A[0] + 2.2], [A[1], A[1]], color="0.3", lw=0.9)
    ax.add_patch(Arc(tuple(A), 3.0, 3.0, theta1=0, theta2=np.degrees(tm), color="0.3",
                     lw=0.9))
    ax.text(A[0] + 1.75, A[1] + 0.62, "α", fontsize=13, ha="left", va="center")
    # the tangent at P, flatter than the chord
    tan = np.array([np.cos(tp), np.sin(tp)])
    ax.plot([P[0] - 2.6 * tan[0], P[0] + 2.6 * tan[0]],
            [P[1] - 2.6 * tan[1], P[1] + 2.6 * tan[1]], color="0.35", lw=1.3, ls="--",
            zorder=5)
    # the three points
    ax.plot(*C, "o", ms=8, mfc="white", mec="black", mew=1.8, zorder=7)
    ax.plot(*M, "o", ms=8, color=WATER, zorder=7)
    ax.plot(*P, "o", ms=7, mfc="white", mec="0.45", mew=1.5, zorder=7)
    ax.text(A[0] - 0.4, A[1] - 0.15, "A", fontsize=12, ha="center", va="center")
    ax.text(B[0] + 0.25, B[1] + 0.35, "B", fontsize=12, ha="left", va="center")

    def label(xy, y, text, color, x=lab):
        ax.annotate(text, xy=tuple(xy), xytext=(x, y), fontsize=10.5, ha="left",
                    va="center", color=color,
                    arrowprops=dict(arrowstyle="-", color=color, lw=0.8,
                                    shrinkA=2, shrinkB=5))

    label(C, 1.3, "C, the middle of the chord AB:\nthe base forces act here", "black",
          x=0.4)
    label(M, -3.6, "M, the middle of the arc AB:\nthe pore pressure is read here", WATER)
    label(P, -6.2, "P, the circle directly below the\nslice's center line, with its\n"
          "tangent (dashed)", "0.4")
    ta = t1 + 0.78 * (t2 - t1)
    ax.text(*(on(ta) * (1.0 + 0.075)), "Δℓ", fontsize=13, ha="left", va="center")
    ax.set_xlim(-1.0, lab + 7.2)
    ax.set_ylim(-10.0, 4.0)
    ax.set_aspect("equal")
    ax.axis("off")
    out = os.path.join(HERE, "circle_slice_base.png")
    fig.savefig(out, dpi=DPI, bbox_inches="tight", pad_inches=0.05)
    plt.close(fig)
    print("wrote", out)


if __name__ == "__main__":
    figure_width_angle()
    figure_base()
