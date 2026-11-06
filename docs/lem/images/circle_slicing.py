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
                arrowprops=dict(arrowstyle="-", color="#c0392b", lw=0.9,
                                shrinkA=3, shrinkB=9))
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
    """One slice from 42 to 82 degrees of a circle of radius 10, centered at O,
    drawn zoomed: the radii run off the panel toward O. Each point carries its
    letter; the key beside the slice says what each letter is."""
    r = 10.0
    t1, t2 = np.radians(42.0), np.radians(82.0)
    tm = 0.5 * (t1 + t2)
    on = lambda t: np.array([r * np.sin(t), -r * np.cos(t)])
    A, B, M = on(t1), on(t2), on(tm)
    C = 0.5 * (A + B)
    tp = np.arcsin(C[0] / r)
    P = on(tp)
    top = 2.6
    x_lo, x_hi, y_lo, y_hi = 4.2, 11.6, -8.9, 3.4

    fig, ax = plt.subplots(figsize=(9.6, 6.9))
    ts = np.linspace(t1, t2, 200)
    arc = np.array([on(t) for t in ts])
    poly = [(A[0], top)] + [tuple(p) for p in arc] + [(B[0], top)]
    ax.add_patch(Polygon(poly, closed=True, fc=SOIL, ec="none", zorder=0))
    ax.plot([A[0] - 1.0, B[0] + 1.0], [top, top], color="black", lw=1.2)
    for x, y in ((A[0], A[1]), (B[0], B[1])):
        ax.plot([x, x], [y, top], ls="--", color="black", lw=0.9)
    # the circle beyond the slice, faint
    tt = np.linspace(np.radians(30), np.radians(88), 300)
    ax.plot(r * np.sin(tt), -r * np.cos(tt), color="0.7", lw=1.0, zorder=1)

    # radii from A, M and B toward the center O, clipped at the panel's left edge
    def toward_o(p, color, lw):
        k = (p[0] - x_lo) / p[0]
        q = p * (1.0 - k)
        ax.plot([p[0], q[0]], [p[1], q[1]], color=color, lw=lw, zorder=1)
        return q
    for p in (A, B):
        toward_o(p, "0.75", 0.8)
    q = toward_o(M, "0.55", 0.9)
    ax.text(q[0] + 0.15, q[1] + 0.35, "to the center O", fontsize=10, color="0.4",
            ha="left", va="bottom", rotation=np.degrees(np.arctan2(-M[1], -M[0])) - 180,
            rotation_mode="anchor")
    # the base: arc and chord
    ax.plot(arc[:, 0], arc[:, 1], color="black", lw=2.6, zorder=4)
    ax.plot([A[0], B[0]], [A[1], B[1]], color="black", lw=1.0, ls=(0, (6, 3)), zorder=4)
    # the slice's center line, where the weight acts
    ax.plot([C[0], C[0]], [P[1] - 0.6, top + 0.3], ls=":", color="0.2", lw=1.1)
    ax.add_patch(FancyArrowPatch((C[0], top - 0.35), (C[0], top - 2.4),
                                 arrowstyle="-|>", mutation_scale=14, color="black", lw=1.4))
    ax.text(C[0] + 0.2, top - 1.35, "W", fontsize=12, ha="left", va="center")
    # the base normal at C, along the radius toward O
    d = -C / np.linalg.norm(C)
    tip = C + 1.7 * d
    ax.add_patch(FancyArrowPatch(tuple(C), tuple(tip), arrowstyle="-|>", shrinkA=5,
                                 mutation_scale=14, color="black", lw=1.4, zorder=6))
    ax.text(tip[0] + 0.42, tip[1] - 0.42, "N", fontsize=12, ha="center", va="center")
    # alpha, the chord's inclination, drawn at A
    ax.plot([A[0], A[0] + 2.0], [A[1], A[1]], color="0.3", lw=0.9)
    ax.add_patch(Arc(tuple(A), 2.8, 2.8, theta1=0, theta2=np.degrees(tm), color="0.3",
                     lw=0.9))
    ax.text(A[0] + 1.62, A[1] + 0.55, "α", fontsize=13, ha="left", va="center")
    # the tangent at P, flatter than the chord
    tan = np.array([np.cos(tp), np.sin(tp)])
    ax.plot([P[0] - 2.4 * tan[0], P[0] + 2.4 * tan[0]],
            [P[1] - 2.4 * tan[1], P[1] + 2.4 * tan[1]], color="0.35", lw=1.3, ls="--",
            zorder=5)
    # the three points and their letters
    ax.plot(*C, "o", ms=8, mfc="white", mec="black", mew=1.8, zorder=7)
    ax.plot(*M, "o", ms=8, color=WATER, zorder=7)
    ax.plot(*P, "o", ms=7, mfc="white", mec="0.4", mew=1.5, zorder=7)
    ax.text(C[0] - 0.42, C[1] - 0.38, "C", fontsize=12, ha="center", va="center",
            zorder=8)
    ax.text(M[0] + 0.48, M[1] + 0.05, "M", fontsize=12, ha="left", va="center",
            color=WATER, zorder=8)
    ax.text(P[0] + 0.48, P[1] - 0.22, "P", fontsize=12, ha="left", va="center",
            color="0.35", zorder=8)
    ax.text(A[0] - 0.42, A[1] - 0.2, "A", fontsize=12, ha="center", va="center")
    ax.text(B[0] + 0.32, B[1] + 0.3, "B", fontsize=12, ha="left", va="center")
    tb = t1 + 0.8 * (t2 - t1)
    ax.text(*(on(tb) * 1.07), "Δℓ", fontsize=13, ha="left", va="center")

    # the key, beside the slice
    key = [("C", "the middle of the chord AB, on the slice's center line\n"
                 "and on the radius through M: the base forces act here", "black"),
           ("M", "the middle of the arc AB: the pore pressure is read here", WATER),
           ("P", "the circle directly below the slice's center line;\n"
                 "dashed, the circle's tangent there", "0.35"),
           ("Δℓ", "the length of the arc AB", "black"),
           ("α", "the inclination of the chord AB, which is\n"
                 "parallel to the circle's tangent at M", "black")]
    y = 1.6
    for sym, text, color in key:
        ax.text(x_hi + 0.5, y, sym, fontsize=12, ha="left", va="top", color=color,
                weight="bold")
        ax.text(x_hi + 1.45, y, text, fontsize=10.5, ha="left", va="top", color=color,
                linespacing=1.25)
        y -= 0.62 + 0.37 * (1 + text.count("\n"))
    ax.set_xlim(x_lo, x_hi + 9.6)
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
