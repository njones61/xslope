"""Generate the XSlope Studio app icon (run once; commit the PNG).

Draws a geotechnical slope-stability motif — a layered slope with a circular
failure surface and a few slices — on a rounded-square tile. Renders to
``icon.png`` (1024², the source the app loads via QIcon and macOS uses for the
dock). Re-run after edits: ``python studio/resources/make_icon.py``.
"""

from __future__ import annotations

import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyBboxPatch, Polygon
from matplotlib.path import Path
from matplotlib.patches import PathPatch

S = 1024
BLUE = "#1f4e79"      # tile background (deep blue)
SKY = "#2d6cb0"
SOIL = "#e7d6a8"      # upper soil (tan)
SOIL2 = "#c9a86a"     # foundation (darker tan)
ARC = "#e8552d"       # failure surface (orange-red)
LINE = "#143a5c"      # ground/profile line


def _clip(ax, patch_xy, artists):
    clip = Polygon(patch_xy, closed=True, transform=ax.transData)
    for a in artists:
        a.set_clip_path(clip)


def main():
    here = os.path.dirname(__file__)
    out = os.path.join(here, "icon.png")

    fig = plt.figure(figsize=(S / 100, S / 100), dpi=100)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, S)
    ax.set_ylim(0, S)
    ax.axis("off")

    # Rounded-square tile.
    pad = 36
    tile = FancyBboxPatch(
        (pad, pad), S - 2 * pad, S - 2 * pad,
        boxstyle="round,pad=0,rounding_size=190",
        linewidth=0, facecolor=BLUE, zorder=0)
    ax.add_patch(tile)
    tile_xy = tile.get_path().vertices  # for clipping children to the rounded tile
    # FancyBboxPatch path is in axes-ish coords; rebuild a simple rounded clip by
    # using the patch itself as the clip path instead.

    tile_clipped = []     # clip to the rounded tile
    soil_clipped = []     # clip to the soil body (so the slip surface stays inside)

    # Ground/profile line: left bench -> slope up -> crest bench.
    gx = [150, 430, 700, 880]
    gy = [430, 430, 690, 690]
    # Upper soil layer (tan) down to the foundation top.
    top = list(zip(gx, gy))
    upper = Polygon(top + [(880, 300), (150, 300)], closed=True,
                    facecolor=SOIL, edgecolor="none", zorder=2)
    ax.add_patch(upper); tile_clipped.append(upper)
    # Foundation layer (darker) below.
    found = Polygon([(150, 300), (880, 300), (880, 175), (150, 175)], closed=True,
                    facecolor=SOIL2, edgecolor="none", zorder=2)
    ax.add_patch(found); tile_clipped.append(found)

    # Circular failure surface: center above the slope, passing through the toe.
    toe = np.array([430, 430])
    cx, cy = 560, 980
    R = np.hypot(cx - toe[0], cy - toe[1])
    th = np.linspace(np.deg2rad(202), np.deg2rad(320), 240)
    ax_arc = cx + R * np.cos(th)
    ay_arc = cy + R * np.sin(th)
    arc, = ax.plot(ax_arc, ay_arc, color=ARC, lw=26, solid_capstyle="round", zorder=4)
    soil_clipped.append(arc)

    # A few slice lines from the ground down to the arc, for the LEM motif.
    for xs in np.linspace(380, 700, 6):
        ytop = np.interp(xs, gx, gy)
        k = np.argmin(np.abs(ax_arc - xs))
        yb = ay_arc[k]
        if yb < ytop:
            ln, = ax.plot([xs, xs], [yb, ytop], color=LINE, lw=5, alpha=0.5, zorder=3)
            soil_clipped.append(ln)

    # Ground line on top.
    gl, = ax.plot(gx, gy, color=LINE, lw=14, solid_capstyle="round",
                  solid_joinstyle="round", zorder=5)
    tile_clipped.append(gl)

    # Brand mark: a bold "X" in the empty upper-left blue field.
    xmark = ax.text(355, 720, "X", ha="center", va="center",
                    fontsize=240, fontweight="bold", color="#f5f1e6",
                    family="DejaVu Sans", zorder=6)
    tile_clipped.append(xmark)

    for a in tile_clipped:
        a.set_clip_path(tile)
    for a in soil_clipped:        # slip surface + slices stay inside the slope body
        a.set_clip_path(upper)

    fig.savefig(out, transparent=True)
    plt.close(fig)
    print("wrote", out)


if __name__ == "__main__":
    main()
