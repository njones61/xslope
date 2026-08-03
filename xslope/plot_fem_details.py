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


def _axis_label(base, unit):
    return f"{base} ({unit})" if unit else base


def _thin_ticks(ax, axis="x", nbins=4):
    """Keep a narrow profile panel's tick labels from colliding."""
    from matplotlib.ticker import MaxNLocator
    loc = MaxNLocator(nbins=nbins, prune=None)
    (ax.xaxis if axis == "x" else ax.yaxis).set_major_locator(loc)


def _annotate_inside(ax, xy, text, color):
    """Annotate a point, offset towards whichever side of the axes has room, so
    the label stays inside the panel instead of running off its edge."""
    x0, x1 = ax.get_xlim()
    y0, y1 = ax.get_ylim()
    fx = (xy[0] - x0) / (x1 - x0) if x1 != x0 else 0.5
    fy = (xy[1] - y0) / (y1 - y0) if y1 != y0 else 0.5
    dx, ha = (-8, "right") if fx > 0.55 else (8, "left")
    dy, va = (-8, "top") if fy > 0.55 else (8, "bottom")
    ax.annotate(text, xy=xy, xytext=(dx, dy), textcoords="offset points",
                color=color, fontsize=8.5, fontweight="bold", ha=ha, va=va,
                zorder=8)


def _title(profile):
    util = profile.get("peak_utilization")
    bits = [profile["label"]]
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
    mechanism field puts on this bar is shaded, the peak-utilization point is
    annotated, and elements that softened or pulled out are marked where they
    occur.

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
                            gridspec_kw={"height_ratios": [3, 1]})
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
                color=C_ENVELOPE, label="Capacity envelope", zorder=3)
    else:
        ax.step(s, profile["t_cap"], where="mid", linestyle="--", linewidth=1.4,
                color=C_ENVELOPE, label="Element capacity", zorder=3)

    # Residual capacity, only where the line actually softens somewhere.
    t_res = np.asarray(profile.get("t_res", []), dtype=float)
    if len(t_res) and np.any(np.isfinite(t_res) & (t_res > 0)):
        ax.step(s, t_res, where="mid", linestyle=":", linewidth=1.1,
                color=C_SOFT, label="Residual capacity", zorder=3)

    # Mobilized force.
    ax.plot(s, T, "-o", color=C_FORCE, linewidth=1.8, markersize=3.5,
            label="Mobilized force", zorder=5)

    soft_s = profile.get("softened_s", [])
    if len(soft_s):
        ax.plot(soft_s, np.interp(soft_s, s, T), "s", color=C_SOFT, markersize=6,
                label="Softened", zorder=6)
    pull_s = profile.get("pullout_s", [])
    if len(pull_s):
        ax.plot(pull_s, np.zeros(len(pull_s)), "x", color=C_PULL, markersize=8,
                markeredgewidth=2, label="Pulled out", zorder=6)

    ax.set_ylabel(_axis_label("Axial force", u.get("force")))
    ax.set_title(_title(profile), fontsize=11)
    ax.grid(True, **GRID)
    ax.set_xlim(0.0, profile.get("length") or (s[-1] if len(s) else 1.0))
    ax.set_ylim(bottom=0.0)
    ax.legend(loc="best", fontsize=8, framealpha=0.85)

    # Peak utilization, after the limits are settled so the label can be placed
    # against them.
    if profile.get("peak_s") is not None:
        ps, pt = profile["peak_s"], profile["peak_T"]
        ax.plot([ps], [pt], "o", color=C_PEAK, markersize=7,
                markerfacecolor="none", markeredgewidth=1.8, zorder=7)
        _annotate_inside(ax, (ps, pt),
                         f"{pt:,.0f}{(' ' + u['force']) if u.get('force') else ''}"
                         f"  ({profile['peak_utilization']:.0%})", C_PEAK)
    if lo is not None:
        from matplotlib.transforms import blended_transform_factory
        ax.text(lo, 0.97, " failure band",
                transform=blended_transform_factory(ax.transData, ax.transAxes),
                color=C_BAND, fontsize=8, va="top", ha="left", zorder=8,
                bbox=dict(facecolor="white", edgecolor="none", alpha=0.75, pad=1))

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

    fig.tight_layout()
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
        ax.plot(x, y, "-", color=C_FORCE, linewidth=1.7, zorder=4)
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
    if band is not None:
        ax_u.annotate("failure band", xy=(0.03, band),
                      xycoords=ax_u.get_yaxis_transform(), color=C_BAND,
                      fontsize=8, va="bottom", ha="left", zorder=8)
    if profile.get("max_moment") is not None:
        mm, md = profile["max_moment"], profile["max_moment_depth"]
        ax_m.plot([mm], [md], "o", color=C_PEAK, markersize=6,
                  markerfacecolor="none", markeredgewidth=1.6, zorder=7)
        _annotate_inside(ax_m, (mm, md),
                         f"Mmax {mm:,.0f}\nat {md:,.2f}"
                         f"{(' ' + u['length']) if u.get('length') else ''}", C_PEAK)

    fig.suptitle(_title(profile), fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    return fig


def plot_detail(profile, fig=None, **kwargs):
    """Dispatch to the figure builder for this profile's member kind."""
    if profile.get("kind") == "pile":
        return plot_pile_detail(profile, fig=fig)
    return plot_reinforcement_detail(profile, fig=fig, **kwargs)
