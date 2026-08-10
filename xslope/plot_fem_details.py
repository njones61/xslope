# Copyright 2025 Norman L. Jones
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Per-member detail figures for the FEM's reinforcement bars and piles.

These are profile plots, not geometry: the horizontal axis of a reinforcement
figure is position along the bar and the vertical axis of a pile figure is depth
along the pile, so neither takes the equal-aspect framing the geometry figures
use. Every series comes from :mod:`xslope.fem_details`, which is also what the
dialog's CSV export writes — the picture and the numbers are the same data.

Legends appear only where a panel genuinely carries more than one series whose
identity is not already in the axis label: the reinforcement force panel (the
mobilized force against its capacity envelope) and the pile soil-reaction panel
(mobilized against the Ito & Matsui limit). The others are single-series and are
labelled by their axis.
"""

import numpy as np

# Colours, kept in one place so the two figures read as one family.
C_FORCE = "#1f5fa9"        # mobilized quantity
C_ENVELOPE = "#333333"     # declared capacity
C_BAND = "#d95f0e"         # failure band
C_PEAK = "#c0392b"         # peak / maximum marker
C_SOFT = "#8e44ad"         # softened to residual
C_PULL = "#000000"         # pulled out
C_LIMIT = "#7f8c8d"        # limiting resistance envelope
GRID = dict(alpha=0.25, linewidth=0.6)


#: The drawing band each stacked panel of a member detail figure gets, as a
#: fraction of the figure's WIDTH. A profile is read for where along the member
#: something happens, so its drawing area is wide and shallow; the bond panel
#: gets half the force panel's band, carrying one curve and no legend. Only the
#: bands are stated — the decoration around them (title, axis labels, tick
#: labels) is measured off the drawn figure, so the same bands hold at any font
#: size and the figure ends up as tall as its content and no taller.
DETAIL_BANDS = (0.20, 0.10)


def _fit_stacked_panels(fig, bands):
    """Resize ``fig`` so its stacked panels get ``bands`` of drawing height.

    The decoration — everything that is not axes — is measured from the figure
    as it stands, then the height is set to that plus the requested bands. A
    detail figure printed at text width then costs a page a strip rather than a
    third of a sheet, which is what lets every member have one.
    """
    fig.tight_layout()
    width, height = fig.get_size_inches()
    drawn = sum(ax.get_position().height for ax in fig.axes) * height
    decoration = max(height - drawn, 0.0)
    wanted = decoration + width * sum(bands)
    if wanted > 0 and abs(wanted - height) > 0.01:
        fig.set_size_inches(width, wanted)
        fig.tight_layout()


def _axis_label(base, unit):
    return f"{base} ({unit})" if unit else base


def _thin_ticks(ax, axis="x", nbins=4):
    """Keep a narrow profile panel's tick labels from colliding."""
    from matplotlib.ticker import MaxNLocator
    loc = MaxNLocator(nbins=nbins, prune=None)
    (ax.xaxis if axis == "x" else ax.yaxis).set_major_locator(loc)


def _panel_obstacles(ax, renderer):
    """``(segments, boxes)`` in display pixels: every curve drawn in this panel,
    and every label and legend already standing in it.

    What a peak annotation has to miss. The curves are read off the artists
    rather than off the profile arrays, so a panel that grows a series gets it
    for free.
    """
    from matplotlib.text import Text
    segments, boxes = [], []
    for line in ax.lines:
        pts = line.get_xydata()
        if pts is None or len(pts) < 2:
            continue
        px = ax.transData.transform(pts)
        segments += list(zip(px[:-1], px[1:]))
    for artist in ax.texts:
        try:
            boxes.append(Text.get_window_extent(artist, renderer))
        except Exception:
            pass
    legend = ax.get_legend()
    if legend is not None:
        try:
            boxes.append(legend.get_window_extent(renderer))
        except Exception:
            pass
    return segments, boxes


def _annotate_inside(ax, xy, text, color, fontsize=8.5, fontweight="bold",
                     leader=False):
    """Annotate a point where the panel has room for it.

    A peak sits on the curve it is the peak OF, and on a member well inside its
    capacity that curve runs along the envelope's own descent — so an offset
    chosen from the point's quadrant alone lands the label on the dashes it is
    meant to be read against. The offset is solved instead: a ring of candidates
    scored on what each would cover — a curve, another label, the legend, the
    edge of the panel — and the cheapest taken. Same rule as the model figure's
    coordinate labels (:func:`xslope.plot.plot_coordinate_labels`), and the same
    box/segment test underneath it.

    Which is why the failure band's own label is placed through here too. Set
    against the band's left edge, it ran off the right of the panel whenever the
    band sat on the last element of the line — the label a reader needs is the
    one that is fully on the page, and where that is depends on the panel.
    """
    import math
    from matplotlib.font_manager import FontProperties
    from matplotlib.transforms import Bbox
    from .plot import _box_crosses_segment

    fig = ax.figure
    try:
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
    except Exception:
        renderer = None
    if renderer is None:
        ax.annotate(text, xy=xy, xytext=(8, 8), textcoords="offset points",
                    color=color, fontsize=fontsize, fontweight=fontweight,
                    ha="left", va="bottom", zorder=8)
        return

    prop = FontProperties(size=fontsize, weight=fontweight)
    # Measured line by line: a two-line label is as wide as its widest line and
    # as tall as both, which is not what a single measurement of the whole
    # string reports.
    rows = str(text).split("\n")
    sizes = [renderer.get_text_width_height_descent(r, prop, False) for r in rows]
    w = max(s[0] for s in sizes)
    h = sum(s[1] for s in sizes) * 1.2
    scale = fig.dpi / 72.0
    px, py = ax.transData.transform(xy)
    frame = ax.get_window_extent(renderer)
    segments, boxes = _panel_obstacles(ax, renderer)

    angles = [math.radians(a) for a in range(0, 360, 15)]
    radii = [fontsize * 1.4 * f for f in (0.8, 1.4, 2.0, 2.8, 3.8, 5.2)]
    best = None
    for r in radii:
        for a in angles:
            dx, dy = r * math.cos(a), r * math.sin(a)
            ha = "left" if dx > 2 else ("right" if dx < -2 else "center")
            va = "bottom" if dy > 2 else ("top" if dy < -2 else "center")
            ox, oy = px + dx * scale, py + dy * scale
            bx = ox - (w if ha == "right" else w / 2 if ha == "center" else 0)
            by = oy - (h if va == "top" else h / 2 if va == "center" else 0)
            bb = Bbox.from_bounds(bx - 2, by - 2, w + 4, h + 4)
            area = max(bb.width * bb.height, 1.0)
            # Outside the panel is unreadable; over a label or the legend hides
            # something already said; over a curve hides the profile itself.
            inside = (max(0.0, min(bb.x1, frame.x1) - max(bb.x0, frame.x0))
                      * max(0.0, min(bb.y1, frame.y1) - max(bb.y0, frame.y0)))
            cost = 9000.0 * (area - inside) / area
            cost += 5000.0 * sum(bb.overlaps(ob) for ob in boxes)
            cost += 300.0 * min(sum(_box_crosses_segment(bb, *s)
                                    for s in segments), 4)
            cost += 4.0 * r / (fontsize * 1.4)
            if best is None or cost < best[0]:
                best = (cost, dx, dy, ha, va)

    _cost, dx, dy, ha, va = best
    # A label the solver had to push well clear of its point needs to say which
    # point it belongs to. Drawn only when it is far enough for that to be a
    # question — a label sitting against the thing it names needs no line to it.
    arrow = None
    if leader and math.hypot(dx, dy) > 2.0 * h / scale:
        arrow = dict(arrowstyle="-", color=color, linewidth=0.7,
                     shrinkA=1.0, shrinkB=1.0, alpha=0.8)
    # On the same white backing the panel's other labels carry: a dense profile
    # has places where every offset is over SOMETHING, and the label that lands
    # on one is read rather than dissolved into the grid behind it.
    ax.annotate(text, xy=xy, xytext=(dx, dy), textcoords="offset points",
                color=color, fontsize=fontsize, fontweight=fontweight, ha=ha, va=va,
                bbox=dict(boxstyle="round,pad=0.15", facecolor="white",
                          edgecolor="none", alpha=0.8),
                arrowprops=arrow, zorder=8)


def _title(profile):
    util = profile.get("peak_utilization")
    bits = [profile["label"]]
    # Which field the profile was read from, disclosed the way the result plots
    # disclose it. Only the at-failure state is named: the converged field is
    # what a plain title has always meant.
    if profile.get("field_state") == "failure":
        bits.append("at failure")
    if util is not None and np.isfinite(util):
        bits.append(f"peak {util:.0%}")
    bits.append(profile.get("status", ""))
    return " — ".join(b for b in bits if b)


# --------------------------------------------------------------------------
# reinforcement
# --------------------------------------------------------------------------

def plot_reinforcement_detail(profile, fig=None, show_bond=True):
    """Draw the detail figure for one reinforcement line onto ``fig``.

    Upper panel: mobilized axial force along the bar, over the dashed capacity
    envelope — the pullout ramp developing from each free end over its
    development length, the tensile plateau at Tmax in the middle, and the step
    to the connection capacity at an anchored end. The failure band the
    mechanism field puts on this bar is shaded, and elements that softened or
    pulled out are marked where they occur.

    The point of greatest utilization is ringed where it is a point, and where
    the line holds that utilization over a stretch — which is the usual case,
    the force being capped by a flat envelope — the whole stretch is drawn
    instead, every sample on it ringed and the extent stated.

    Lower panel (``show_bond``): the bond transfer rate implied by the force
    gradient, dT/ds — the force the ground hands the bar per unit of its length.
    There is no companion slip series because the formulation has no slip degree
    of freedom; the bar sits on the continuum's own nodes.
    """
    import matplotlib.pyplot as plt

    if fig is None:
        fig = plt.figure(figsize=(9.5, 6.0))
    fig.clear()
    u = profile.get("units", {}) or {}

    has_bond = show_bond and len(profile.get("bond_s", [])) > 0
    if has_bond:
        axes = fig.subplots(2, 1, sharex=True,
                            gridspec_kw={"height_ratios": list(DETAIL_BANDS)})
        ax, ax_b = axes
    else:
        ax = fig.subplots(1, 1)
        ax_b = None

    s, T = profile["s"], profile["T"]

    # Failure band first, so everything else draws over it.
    lo, hi = profile.get("band_lo"), profile.get("band_hi")
    if lo is not None and hi is not None:
        if hi - lo < 1e-9:                 # a single element: a rule, not a band
            for a in (ax, ax_b):
                if a is not None:
                    a.axvline(lo, color=C_BAND, linewidth=1.2, linestyle=(0, (4, 3)),
                              zorder=1)
        else:
            for a in (ax, ax_b):
                if a is not None:
                    a.axvspan(lo, hi, color=C_BAND, alpha=0.13, zorder=0, linewidth=0)

    # Capacity envelope.
    if profile.get("env_s") is not None:
        ax.plot(profile["env_s"], profile["env_T"], linestyle="--", linewidth=1.4,
                color=C_ENVELOPE, label="Capacity envelope", zorder=3,
                gid="DETAIL_CAPACITY")
    else:
        ax.step(s, profile["t_cap"], where="mid", linestyle="--", linewidth=1.4,
                color=C_ENVELOPE, label="Element capacity", zorder=3,
                gid="DETAIL_CAPACITY")

    # Residual capacity, only where the line actually softens somewhere.
    t_res = np.asarray(profile.get("t_res", []), dtype=float)
    if len(t_res) and np.any(np.isfinite(t_res) & (t_res > 0)):
        ax.step(s, t_res, where="mid", linestyle=":", linewidth=1.1,
                color=C_SOFT, label="Residual capacity", zorder=3)

    # Mobilized force.
    ax.plot(s, T, "-o", color=C_FORCE, linewidth=1.8, markersize=3.5,
            label="Mobilized force", zorder=5, gid="DETAIL_PROFILE")

    soft_s = profile.get("softened_s", [])
    if len(soft_s):
        ax.plot(soft_s, np.interp(soft_s, s, T), "s", color=C_SOFT, markersize=6,
                label="Softened", zorder=6)
    pull_s = profile.get("pullout_s", [])
    if len(pull_s):
        ax.plot(pull_s, np.zeros(len(pull_s)), "x", color=C_PULL, markersize=8,
                markeredgewidth=2, label="Pulled out", zorder=6)

    ax.set_ylabel(_axis_label("Axial force", u.get("force")), fontsize=9)
    ax.set_title(_title(profile), fontsize=11)
    ax.grid(True, **GRID)
    ax.set_xlim(0.0, profile.get("length") or (s[-1] if len(s) else 1.0))
    ax.set_ylim(bottom=0.0)
    ax.legend(loc="best", fontsize=8, framealpha=0.85)

    # Where the line stands at its greatest utilization: one point, or the
    # stretch it holds that utilization over. A stretch is drawn as a stretch —
    # the run of curve it covers, thickened, with every sample on it ringed. The
    # single ring the figure used to carry sat on the first sample of that run
    # and marked it as the one place the reader should look at.
    span = profile.get("peak_span")
    tied = np.asarray(profile.get("peak_indices", []), dtype=int)
    if span is not None and len(tied):
        # Run by run: a sample between two tied ones that is NOT itself at the
        # maximum breaks the highlight, so what is thickened is the curve the
        # line is actually at capacity along and never a chord across a dip.
        for run in np.split(tied, np.where(np.diff(tied) > 1)[0] + 1):
            if len(run) > 1:
                ax.plot(s[run], T[run], "-", color=C_PEAK, linewidth=4.5,
                        alpha=0.30, solid_capstyle="butt", zorder=6,
                        gid="DETAIL_PEAK_SPAN")
        ax.plot(s[tied], T[tied], "o", color=C_PEAK, markersize=6,
                markerfacecolor="none", markeredgewidth=1.5, zorder=7)
    elif profile.get("peak_s") is not None:
        ps, pt = profile["peak_s"], profile["peak_T"]
        ax.plot([ps], [pt], "o", color=C_PEAK, markersize=7,
                markerfacecolor="none", markeredgewidth=1.8, zorder=7)
    if ax_b is not None:
        ax_b.axhline(0.0, color="0.6", linewidth=0.8)
        ax_b.plot(profile["bond_s"], profile["bond_q"], "-", color=C_FORCE,
                  linewidth=1.4)
        ax_b.set_ylabel(_axis_label("Bond transfer\ndT/ds", u.get("line_load")),
                        fontsize=9)
        ax_b.grid(True, **GRID)
        ax_b.set_xlabel(_axis_label("Position along line", u.get("length")))
    else:
        ax.set_xlabel(_axis_label("Position along line", u.get("length")))

    _fit_stacked_panels(fig, DETAIL_BANDS if has_bond else DETAIL_BANDS[:1])

    # The labels go on last, on the final layout: each is placed against the
    # curves, the legend and the labels already standing, and all of those move
    # when the panels are sized to their bands. The band's label goes first, so
    # the one naming the utilization is placed knowing where it ended up.
    if lo is not None:
        band_x = 0.5 * (lo + hi) if hi is not None else lo
        _annotate_inside(ax, (band_x, ax.get_ylim()[1]), "failure band", C_BAND,
                         fontsize=8, fontweight="normal")

    if span is not None and len(tied):
        mid = 0.5 * (span[0] + span[1])
        _annotate_inside(
            ax, (mid, float(np.interp(mid, s, T))),
            f"{profile['peak_utilization']:.0%} of capacity\n"
            f"from {span[0]:,.2f} to {span[1]:,.2f}"
            f"{(' ' + u['length']) if u.get('length') else ''}", C_PEAK)
    elif profile.get("peak_s") is not None:
        _annotate_inside(
            ax, (profile["peak_s"], profile["peak_T"]),
            f"{profile['peak_T']:,.0f}"
            f"{(' ' + u['force']) if u.get('force') else ''}"
            f"  ({profile['peak_utilization']:.0%})", C_PEAK)
    return fig


# --------------------------------------------------------------------------
# piles
# --------------------------------------------------------------------------

def plot_pile_detail(profile, fig=None):
    """Draw the detail figure for one pile onto ``fig``.

    Four panels sharing one depth axis, pile head at the top: lateral
    displacement, shear, bending moment, and the lateral soil reaction with the
    Ito & Matsui limiting resistance dashed beside it. The maximum-moment depth
    is marked, and the depth at which the failure band crosses the pile is ruled
    across all four so the profiles can be read against the mechanism.

    Capacity lines are drawn only for the capacities the model declares. Shear
    and moment capacities come from the ``Vcap`` and ``Mcap`` inputs; there is no
    fallback, because the pile inputs carry force capacities and not section
    properties, and a capacity computed from an assumed section would be an
    invention rather than a reading of the model.
    """
    import matplotlib.pyplot as plt

    if fig is None:
        fig = plt.figure(figsize=(11.0, 6.0))
    fig.clear()
    u = profile.get("units", {}) or {}
    axes = fig.subplots(1, 4, sharey=True)
    ax_u, ax_v, ax_m, ax_p = axes

    depth_max = profile.get("length") or 1.0

    def _panel(ax, x, y, xlabel, cap=None, cap_label=None):
        ax.axvline(0.0, color="0.6", linewidth=0.8, zorder=1)
        ax.plot(x, y, "-", color=C_FORCE, linewidth=1.7, zorder=4,
                gid="DETAIL_PROFILE")
        if cap is not None and np.isfinite(cap) and cap > 0:
            lim = max(np.max(np.abs(x)) if len(x) else 0.0, 0.0)
            if cap <= 3.0 * max(lim, 1e-12):     # in range: worth drawing
                from matplotlib.transforms import blended_transform_factory
                for sgn in (-1.0, 1.0):
                    ax.axvline(sgn * cap, color=C_ENVELOPE, linestyle="--",
                               linewidth=1.1, zorder=3)
                ax.text(cap, 0.02, f"{cap_label} ",
                        transform=blended_transform_factory(ax.transData,
                                                            ax.transAxes),
                        color=C_ENVELOPE, fontsize=8, va="bottom", ha="right",
                        rotation=90, zorder=8)
        ax.set_xlabel(xlabel, fontsize=9)
        ax.grid(True, **GRID)
        ax.tick_params(labelsize=8)
        _thin_ticks(ax, "x", 4)
        for lbl in ax.get_xticklabels():
            lbl.set_rotation(30)
            lbl.set_horizontalalignment("right")

    _panel(ax_u, profile["u_lateral"], profile["node_depth"],
           _axis_label("Lateral displacement", u.get("length")))
    _panel(ax_v, profile["shear"], profile["elem_depth"],
           _axis_label("Shear", u.get("force")), profile.get("V_cap"), "Vcap")
    _panel(ax_m, profile["moment"], profile["moment_depth"],
           _axis_label("Moment", u.get("moment")), profile.get("M_cap"), "Mcap")

    # Soil reaction, with the limiting-resistance envelope beside it.
    ax_p.axvline(0.0, color="0.6", linewidth=0.8, zorder=1)
    ax_p.plot(profile["reaction"], profile["reaction_depth"], "-", color=C_FORCE,
              linewidth=1.7, label="Mobilized", zorder=4)
    if profile.get("limit_p") is not None:
        lp = np.asarray(profile["limit_p"], dtype=float)
        ld = np.asarray(profile["limit_depth"], dtype=float)
        mob = float(np.max(np.abs(profile["reaction"]))) if len(profile["reaction"]) else 0.0
        lo_lim = float(np.nanmin(lp)) if np.any(np.isfinite(lp)) else np.inf
        # The limiting resistance grows with depth, and for a pile well inside
        # its working range it sits far above anything mobilized. Letting it set
        # the scale would flatten the mobilized profile onto the axis, so the
        # panel is scaled to the mobilized profile and the limit is drawn over
        # it, running off the sides wherever it exceeds that. The peak
        # mobilization is stated in the panel either way, which is the reading
        # the comparison exists to give and does not depend on the scale.
        half = 1.3 * mob if mob > 0 else (1.1 * lo_lim if np.isfinite(lo_lim) else 0.0)
        if half > 0:
            ax_p.set_xlim(-half, half)
        if lo_lim <= half:
            ax_p.plot(lp, ld, "--", color=C_LIMIT, linewidth=1.2,
                      label="Ito & Matsui limit", zorder=3)
            ax_p.plot(-lp, ld, "--", color=C_LIMIT, linewidth=1.2, zorder=3)
            ax_p.legend(loc="lower right", fontsize=7.5, framealpha=0.85)
        ratio = profile.get("reaction_ratio")
        if ratio is not None:
            ax_p.text(0.5, 0.985, f"peak {ratio:.0%} of limit",
                      transform=ax_p.transAxes, fontsize=7.5, color=C_LIMIT,
                      ha="center", va="top", zorder=8,
                      bbox=dict(facecolor="white", edgecolor="none", alpha=0.75,
                                pad=1))
    ax_p.set_xlabel(_axis_label("Soil reaction", u.get("line_load")), fontsize=9)
    ax_p.grid(True, **GRID)
    ax_p.tick_params(labelsize=8)
    _thin_ticks(ax_p, "x", 4)
    for lbl in ax_p.get_xticklabels():
        lbl.set_rotation(30)
        lbl.set_horizontalalignment("right")

    ax_u.set_ylabel(_axis_label("Depth below pile head", u.get("length")))
    ax_u.set_ylim(depth_max, 0.0)

    # Depth rules that belong to every panel.
    band = profile.get("band_depth")
    for ax in axes:
        if band is not None:
            ax.axhline(band, color=C_BAND, linestyle=(0, (4, 3)), linewidth=1.2,
                       zorder=2)
        if profile.get("max_moment_depth") is not None:
            ax.axhline(profile["max_moment_depth"], color=C_PEAK, linewidth=0.8,
                       alpha=0.35, zorder=2)
    if profile.get("max_moment") is not None:
        mm, md = profile["max_moment"], profile["max_moment_depth"]
        ax_m.plot([mm], [md], "o", color=C_PEAK, markersize=6,
                  markerfacecolor="none", markeredgewidth=1.6, zorder=7)

    fig.suptitle(_title(profile), fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.95))

    # The labels go on last, on the final layout: the moment panel is one of
    # four sharing a depth axis and is narrow, so where a label fits is a
    # question about the drawn panel and not about the moment. The band's label
    # is placed by the same solver, on the same terms.
    if band is not None:
        _annotate_inside(ax_u, (ax_u.get_xlim()[0], band), "failure band",
                         C_BAND, fontsize=8, fontweight="normal")
    if profile.get("max_moment") is not None:
        _annotate_inside(ax_m, (profile["max_moment"],
                                profile["max_moment_depth"]),
                         f"Mmax {profile['max_moment']:,.0f}\n"
                         f"at {profile['max_moment_depth']:,.2f}"
                         f"{(' ' + u['length']) if u.get('length') else ''}",
                         C_PEAK, fontsize=8.0)
    return fig


# --------------------------------------------------------------------------
# where the members are
# --------------------------------------------------------------------------

#: How far past the members the locator looks, as a fraction of the larger side
#: of their own bounding box. A map of the bars alone is a map of nothing: the
#: reader is trying to place them on the slope, so the frame carries enough of
#: the section around them to be recognisable, and then stops — the whole
#: section drawn small would put the members back where they were unreadable.
MAP_CONTEXT = 0.55

#: The uniform cushion the frame keeps around everything it draws, as a fraction
#: of its larger dimension — the geometry figures' own (:mod:`xslope.plot_fem`).
MAP_PAD = 0.035

#: The tallest the frame is allowed to be, as a fraction of its width. Piles run
#: down rather than along, so their own bounding box is tall and narrow, and an
#: equal-aspect frame around it printed at text width is a figure taller than the
#: page it locates them on. Past this the frame is widened rather than the
#: drawing squashed — it takes in more of the section, which is what a locator is
#: for — and never past the section's own extent.
MAP_MAX_ASPECT = 0.6

C_MAP_MEMBER = "#7f8c8d"   # a member of the kind being mapped
C_MAP_OTHER = "#bdc3c7"    # a member of the other kind, drawn for context


def plot_member_map(fem_data, slope_data=None, kind="reinforcement",
                    highlight=None, ax=None, fig=None, labels=True):
    """Draw where one kind of member sits in the section.

    The section's outline, every member the analysis solved, and — where
    ``highlight`` names one by its index — that member picked out from the rest.
    ``labels`` names each member of ``kind`` where it lies, so a row of the
    forces table and a member on the slope are the same member. A highlighted
    member is named either way — a panel too narrow to carry every name still
    has to say which one it is picking out.

    One drawing serves two readers: the report prints it once above the member
    subsection's detail figures, and the Studio details dialog shows it beside
    the member list, highlighting whichever member is selected. They are the same
    picture because they are the same call.
    """
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection
    from .fem_details import member_lines
    from .plot_fem import _get_mesh_boundary

    if ax is None:
        if fig is None:
            fig = plt.figure(figsize=(6.5, 3.5))
        fig.clear()
        ax = fig.subplots(1, 1)
    fig = ax.figure

    members = member_lines(fem_data, slope_data, kind)
    others = member_lines(fem_data, slope_data,
                          "reinforcement" if kind == "pile" else "pile")

    nodes = np.asarray(fem_data.get("nodes", np.zeros((0, 2))), dtype=float)
    segs = [[nodes[a], nodes[b]] for a, b in _get_mesh_boundary(fem_data)]
    if segs:
        ax.add_collection(LineCollection(segs, colors="0.55", linewidths=1.0,
                                         alpha=0.8, zorder=1, gid="MESH"))

    for member in others:
        (x1, y1), (x2, y2) = member["ends"]
        ax.plot([x1, x2], [y1, y2], "-", color=C_MAP_OTHER, linewidth=2.0,
                zorder=2)
    for member in members:
        (x1, y1), (x2, y2) = member["ends"]
        picked = highlight is not None and member["index"] == highlight
        ax.plot([x1, x2], [y1, y2], "-",
                color=C_PEAK if picked else C_MAP_MEMBER,
                linewidth=3.4 if picked else 2.2, zorder=4 if picked else 3,
                solid_capstyle="round")

    # The frame: the members, with enough of the section around them to place
    # them, never more of it than the section has.
    xs = [v for m in members for v in (m["ends"][0][0], m["ends"][1][0])]
    ys = [v for m in members for v in (m["ends"][0][1], m["ends"][1][1])]
    if xs and len(nodes):
        grow = MAP_CONTEXT * max(max(xs) - min(xs), max(ys) - min(ys), 1e-9)
        mx0, mx1 = float(nodes[:, 0].min()), float(nodes[:, 0].max())
        my0, my1 = float(nodes[:, 1].min()), float(nodes[:, 1].max())
        x0, x1 = max(min(xs) - grow, mx0), min(max(xs) + grow, mx1)
        y0, y1 = max(min(ys) - grow, my0), min(max(ys) + grow, my1)
        # A frame taller than MAP_MAX_ASPECT of its width is widened until it is
        # not, as far as the section reaches.
        wanted = (y1 - y0) / MAP_MAX_ASPECT
        if wanted > x1 - x0:
            spare = 0.5 * (wanted - (x1 - x0))
            x0, x1 = max(x0 - spare, mx0), min(x1 + spare, mx1)
        pad = MAP_PAD * max(x1 - x0, y1 - y0, 1e-9)
        ax.set_xlim(x0 - pad, x1 + pad)
        ax.set_ylim(y0 - pad, y1 + pad)
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    for side in ax.spines.values():
        side.set_visible(False)

    # Names last, on the drawn frame: each goes where the panel has room for it,
    # by the same solver the detail figures' labels use — the members are drawn
    # as lines, so a name never lands on one.
    # A picked-out member is named whether or not the rest are: a map drawn in a
    # panel too narrow for six names still has to say which one is picked out.
    named = [m for m in members
             if labels or (highlight is not None and m["index"] == highlight)]
    if named:
        fig.tight_layout()
        for member in named:
            (x1, y1), (x2, y2) = member["ends"]
            picked = highlight is not None and member["index"] == highlight
            _annotate_inside(ax, (0.5 * (x1 + x2), 0.5 * (y1 + y2)),
                             member["label"],
                             C_PEAK if picked else "#333333", fontsize=8,
                             fontweight="bold" if picked else "normal",
                             leader=True)
    return fig


def plot_detail(profile, fig=None, **kwargs):
    """Dispatch to the figure builder for this profile's member kind."""
    if profile.get("kind") == "pile":
        return plot_pile_detail(profile, fig=fig)
    return plot_reinforcement_detail(profile, fig=fig, **kwargs)
