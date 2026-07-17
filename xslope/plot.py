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

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.path import Path
from shapely.geometry import LineString
import math

from .slice import generate_failure_surface, domain_lower_envelope

# Configure matplotlib for better text rendering
plt.rcParams.update({
    "text.usetex": False,
    "font.family": "sans-serif",
    "font.size": 10
})

# Consistent color for materials (Tableau tab10)
def get_material_color(idx):
    cmap = plt.get_cmap("tab10")
    # Use the colormap callable rather than relying on the (typing-unknown) `.colors` attribute.
    return cmap(idx % getattr(cmap, "N", 10))


# Header/usage colors mirrored from studio.editors (LEM red, seepage green, FEM
# blue) so a material's strength/kr plots read as the same color family as the
# columns they visualize.
_STRENGTH_COLOR = "#c00000"   # LEM
_KR_COLOR = "#2e7d32"         # seepage
_REF_COLOR = "#0432ff"        # reference elevation / annotation


def _material_hint(ax, msg, title=""):
    """A blank axes with a centered hint — shown when a strength/kr model's
    parameters are missing, so an incomplete material draws a prompt rather than a
    broken plot."""
    ax.text(0.5, 0.5, msg, ha="center", va="center", transform=ax.transAxes,
            fontsize=10, color="#555555", wrap=True)
    ax.set_xticks([])
    ax.set_yticks([])
    if title:
        ax.set_title(title)
    return ax


def plot_material_strength(ax, material, n=200, sigma_max=100.0):
    """Draw one material's shear-strength model into ``ax``, EXACTLY as the solver
    evaluates it — no re-derivation. The four ``option`` values each get the plot
    that confirms their model:

      - ``mc``  Mohr-Coulomb line  τ = c + σ′·tan(φ)  in τ–σ′ space (slice.py).
      - ``cp``  undrained strength profile  Sᵤ = c + cp·max(0, r_elev − y),
                plotted vs elevation y with r_elev annotated (slice.py:1841).
      - ``pow`` power curve  τ = a·(σ′+d)^b + cₚ  over a σ′ range, matching the
                envelope solve._pow_update_strength re-linearizes (solve.py:2335).
      - ``hb``  generalized Hoek-Brown envelope in τ–σₙ form, traced point-by-point
                through hoekbrown.hb_tangent — the same instantaneous tangent the
                LEM/FEM consume via solve._hb_update_strength (solve.py:2314).

    A blank ``option`` (valid for seep-only materials) or missing parameters draw a
    centered hint instead of a broken plot. ``sigma_max`` sets the σ′/σₙ range for
    the mc/pow lines; the hb range is derived from σci and the cp range from the
    c/cp gradient so each is self-scaling. Pure — takes only an Axes and a material
    dict, so studio panels and static figures share it.
    """
    option = str(material.get("option", "") or "").strip().lower()

    def g(key):
        try:
            return float(material.get(key, 0) or 0)
        except (TypeError, ValueError):
            return 0.0

    if option == "mc":
        c, phi = g("c"), g("phi")
        if c == 0 and phi == 0:
            return _material_hint(ax, "enter c and φ", "Mohr–Coulomb")
        s = np.linspace(0.0, sigma_max, n)
        tau = c + s * np.tan(np.radians(phi))
        ax.plot(s, tau, color=_STRENGTH_COLOR, lw=2)
        ax.set_xlabel("σ′  (effective normal stress)")
        ax.set_ylabel("τ  (shear strength)")
        ax.set_title(f"Mohr–Coulomb   (c={c:g}, φ={phi:g}°)")
        ax.set_xlim(left=0)
        ax.set_ylim(bottom=0)
        ax.grid(True, alpha=0.3)

    elif option == "cp":
        c, cp, r_elev = g("c"), g("cp"), g("r_elev")
        if c == 0 and cp == 0:
            return _material_hint(ax, "enter c (and cp)", "c-with-depth")
        # Elevation window self-scaled from the gradient: show enough below r_elev
        # for Sᵤ to change by ~c, plus a short constant band above. No geometry is
        # available to a pure material plot, so the gradient sets the scale.
        if cp != 0:
            span = min(max(max(abs(c), 1.0) / abs(cp), 0.5), 50.0)
        else:
            span = 5.0
        y = np.linspace(r_elev - span, r_elev + 0.2 * span, n)
        su = c + cp * np.maximum(0.0, r_elev - y)
        ax.plot(su, y, color=_STRENGTH_COLOR, lw=2)
        ax.axhline(r_elev, color=_REF_COLOR, ls="--", lw=1)
        import matplotlib.transforms as _mt
        tr = _mt.blended_transform_factory(ax.transAxes, ax.transData)
        ax.text(0.03, r_elev, f"r_elev = {r_elev:g}", transform=tr,
                va="bottom", ha="left", color=_REF_COLOR, fontsize=9)
        ax.set_xlabel("Sᵤ  (undrained strength)")
        ax.set_ylabel("elevation, y")
        ax.set_title(f"c-with-depth   (c={c:g}, cp={cp:g})")
        ax.set_xlim(left=0)
        ax.grid(True, alpha=0.3)

    elif option == "pow":
        a, b, cp_c, d = g("pow_a"), g("pow_b"), g("pow_c"), g("pow_d")
        if a == 0:
            return _material_hint(ax, "enter pow_a … pow_d", "Power curve")
        s = np.linspace(0.0, sigma_max, n)
        # Mirror solve._pow_update_strength's evaluation: (σ′+d) floored so a b<1
        # curve stays finite near the origin.
        s_eff = np.maximum(s + d, 1e-4 * max(1.0, sigma_max))
        tau = a * np.power(s_eff, b) + cp_c
        ax.plot(s, tau, color=_STRENGTH_COLOR, lw=2)
        ax.set_xlabel("σ′  (effective normal stress)")
        ax.set_ylabel("τ  (shear strength)")
        ax.set_title(f"Power curve   τ = {a:g}·(σ′+{d:g})^{b:g} + {cp_c:g}")
        ax.set_xlim(left=0)
        ax.set_ylim(bottom=0)
        ax.grid(True, alpha=0.3)

    elif option == "hb":
        sci, gsi, mi, dfac = g("hb_sci"), g("hb_gsi"), g("hb_mi"), g("hb_d")
        if sci <= 0 or gsi <= 0 or mi <= 0:
            return _material_hint(ax, "enter hb_sci, hb_gsi, hb_mi", "Hoek–Brown")
        from .hoekbrown import hb_tangent
        # Range derived from σci so both stiff rock and weak rock mass show the
        # envelope's curvature; the low-stress end (curvature) is what matters for
        # slopes. Evaluated through the SAME hb_tangent the solver linearizes with:
        # τ(σₙ) = c_i + σₙ·tan(φ_i).
        smax = max(0.2 * sci, 1.0)
        s = np.linspace(0.0, smax, n)
        c_i, phi_i = hb_tangent(s, sci, gsi, mi, dfac)
        tau = c_i + s * np.tan(np.radians(phi_i))
        ax.plot(s, tau, color=_STRENGTH_COLOR, lw=2)
        ax.set_xlabel("σₙ  (normal stress)")
        ax.set_ylabel("τ  (shear strength)")
        ax.set_title(f"Hoek–Brown   (σci={sci:g}, GSI={gsi:g}, mi={mi:g})")
        ax.set_xlim(left=0)
        ax.set_ylim(bottom=0)
        ax.grid(True, alpha=0.3)

    else:
        return _material_hint(ax, "no strength model\n(blank option)")

    return ax


def plot_material_kr(ax, material, n=200):
    """Draw one material's unsaturated relative-conductivity curve kr vs matric
    suction into ``ax``, using seep.py's OWN kr functions (imported, never
    re-derived), dispatched on the material's ``unsat`` model:

      - ``lf``   linear front — seep.kr_frontal_vec (params kr0, h0).
      - ``vg``   van Genuchten–Mualem — seep.kr_vg_vec (params vg_a, vg_n).
      - ``gard`` Gardner power form — seep.kr_gardner_vec (params vg_a=a, vg_n=n).

    kr is on a log axis (it spans decades) with the suction range self-scaled to
    each model so the full wet→dry decline is framed. Missing/invalid parameters
    draw a centered hint. Pure — Axes + material dict only.
    """
    from .seep import kr_frontal_vec, kr_vg_vec, kr_gardner_vec
    kr_min = 1e-4
    unsat = str(material.get("unsat", "lf") or "lf").strip().lower()

    def g(key):
        try:
            return float(material.get(key, 0) or 0)
        except (TypeError, ValueError):
            return 0.0

    if unsat == "lf":
        kr0, h0 = g("kr0"), g("h0")
        if not (kr0 > 0 and h0 < 0):
            return _material_hint(ax, "enter kr0 > 0 and h0 < 0", "Linear front (lf)")
        smax = 1.3 * abs(h0)
        psi = np.linspace(0.0, smax, n)
        kr = kr_frontal_vec(-psi, kr0, h0)
        title = f"Linear front   (kr0={kr0:g}, h0={h0:g})"
    elif unsat in ("vg", "gard"):
        a, nn = g("vg_a"), g("vg_n")
        if unsat == "vg":
            if not (a > 0 and nn > 1):
                return _material_hint(ax, "enter vg_a > 0 and vg_n > 1",
                                      "van Genuchten (vg)")
            smax = 20.0 / a
            psi = np.linspace(0.0, smax, n)
            kr = kr_vg_vec(-psi, a, nn)
            title = f"van Genuchten   (α={a:g}, n={nn:g})"
        else:
            if not (a > 0 and nn > 0):
                return _material_hint(ax, "enter vg_a > 0 and vg_n > 0",
                                      "Gardner (gard)")
            # Suction at which kr reaches the floor: a·ψⁿ = 1/kr_min − 1.
            smax = 1.1 * (max(1.0 / kr_min - 1.0, 1.0) / a) ** (1.0 / nn)
            psi = np.linspace(0.0, smax, n)
            kr = kr_gardner_vec(-psi, a, nn)
            title = f"Gardner   (a={a:g}, n={nn:g})"
    else:
        return _material_hint(ax, "no unsaturated model")

    ax.plot(psi, kr, color=_KR_COLOR, lw=2)
    ax.set_yscale("log")
    ax.set_ylim(top=1.5)
    # Frame the actual decline: trim the flat floor tail (kr monotone-decreasing,
    # so the last sampled value is the asymptotic floor — kr0 for lf, kr_min for
    # vg/gard) so a sharp curve doesn't sit in a sea of dead space.
    floor = float(kr[-1])
    reached = np.where(kr <= floor * 1.02 + 1e-12)[0]
    right = psi[reached[0]] * 1.1 if len(reached) else psi[-1]
    ax.set_xlim(0, right if right > 0 else psi[-1])
    ax.set_xlabel("matric suction, ψ")
    ax.set_ylabel("relative conductivity, kr")
    ax.set_title(title)
    ax.grid(True, which="both", alpha=0.3)
    return ax

def get_dload_legend_handler():
    """
    Creates and returns a custom legend entry for distributed loads.
    Returns a tuple of (handler_class, dummy_patch) for use in matplotlib legends.
    """
    # Create a line with built-in arrow marker
    dummy_line = Line2D([0.0, 1.0], [0, 0],  # Two points to define line
                       color='purple', 
                       alpha=0.7, 
                       linewidth=2,
                       marker='>',  # Built-in right arrow marker
                       markersize=6,  # Smaller marker size
                       markerfacecolor='purple',
                       markeredgecolor='purple',
                       drawstyle='steps-post',  # Draw line then marker
                       solid_capstyle='butt')
    
    return None, dummy_line


def plot_profile_lines(ax, profile_lines, materials=None, labels=False, style=None):
    """
    Plots the profile lines for each material in the slope.

    Parameters:
        ax: matplotlib Axes object
        profile_lines: List of profile line dicts, each with 'coords' and 'mat_id' keys
        materials: List of material dictionaries (optional, for color mapping)
        labels: If True, add index labels to each profile line (default: False)
        style: optional style sheet (see xslope.style); None → defaults. Line color
            comes from the material; the profile_line feature owns width/linestyle.

    Returns:
        None
    """
    from .style import resolve_style, material_style, feature_style
    style = resolve_style(style)
    fs = feature_style(style, "profile_line")
    for i, line in enumerate(profile_lines):
        coords = line['coords']
        xs, ys = zip(*coords)

        # Material index from mat_id (already 0-based); fall back to the line index.
        mat_idx = line.get('mat_id') if (materials and line.get('mat_id') is not None) else i
        if not (materials and 0 <= mat_idx < len(materials)):
            mat_idx = i
        color = material_style(style, mat_idx)["color"]

        ax.plot(xs, ys, color=color, linewidth=fs.get("linewidth", 1.0),
                linestyle=fs.get("linestyle", "-"),
                label=f'Profile {i+1}', gid=f'PROFILE_{i+1}')

        if labels:
            _add_profile_index_label(ax, coords, i + 1, color)


def _add_profile_index_label(ax, line, index, color):
    """
    Adds an index label to a profile line, positioned on a suitable segment.

    Parameters:
        ax: matplotlib Axes object
        line: List of (x, y) coordinates for the profile line
        index: The index number to display (1-based)
        color: Color for the label text

    Returns:
        None
    """
    if len(line) < 2:
        return

    # Build list of segments with their properties
    segments = []
    for j in range(len(line) - 1):
        x1, y1 = line[j]
        x2, y2 = line[j + 1]
        dx = x2 - x1
        dy = y2 - y1
        length = np.sqrt(dx**2 + dy**2)

        if length < 1e-9:
            continue

        # Calculate how horizontal the segment is (1.0 = perfectly horizontal)
        horizontalness = abs(dx) / length

        # Midpoint of the segment
        mid_x = (x1 + x2) / 2
        mid_y = (y1 + y2) / 2

        segments.append({
            'length': length,
            'horizontalness': horizontalness,
            'mid_x': mid_x,
            'mid_y': mid_y,
            'position': j  # segment index in the line
        })

    if not segments:
        return

    # Calculate the total line length and find segments in the middle third
    total_length = sum(s['length'] for s in segments)

    # Score segments: prefer longer, more horizontal segments near the middle
    # Avoid very short segments (less than 5% of total length)
    min_length_threshold = 0.05 * total_length

    n_segments = len(segments)
    best_segment = None
    best_score = -1

    for seg in segments:
        # Skip very short segments
        if seg['length'] < min_length_threshold:
            continue

        # Calculate how close to the middle this segment is (0-1 scale, 1 = center)
        position_ratio = (seg['position'] + 0.5) / n_segments
        middle_score = 1.0 - 2.0 * abs(position_ratio - 0.5)  # 1.0 at center, 0.0 at ends

        # Score: weight length, horizontalness, and middle position
        # Length is normalized by total length
        length_score = seg['length'] / total_length

        # Combined score: prioritize horizontalness, then middle position, then length
        score = (seg['horizontalness'] * 2.0 +
                 middle_score * 1.5 +
                 length_score * 1.0)

        if score > best_score:
            best_score = score
            best_segment = seg

    # Fallback: if no segment passed the threshold, use the longest one
    if best_segment is None:
        best_segment = max(segments, key=lambda s: s['length'])

    # Place the label at the midpoint of the chosen segment
    ax.text(
        best_segment['mid_x'],
        best_segment['mid_y'],
        str(index),
        fontsize=7,
        color=color,
        fontfamily='monospace',
        ha='center',
        va='center',
        bbox=dict(
            boxstyle='round,pad=0.3',
            facecolor='white',
            edgecolor=color,
            linewidth=0.8
        ),
        zorder=10
    )


def plot_polygons_on_ax(ax, polygons, materials=None, labels=False, style=None):
    """
    Fill the material-zone polygons on an existing Axes (the polygon analog of
    plot_profile_lines). Used by plot_inputs when geometry is defined by polygons
    rather than profile lines.

    Parameters:
        ax: matplotlib Axes object
        polygons: list of dicts with a 'polygon' (shapely Polygon) and 'mat_id' key
        materials: list of material dicts (for color + legend names), optional
        labels: if True, place the material number at each polygon centroid
        style: optional style sheet (see xslope.style); None → defaults. Fill
            color/alpha/hatch come from the material; the polygon feature owns the
            edge line width.
    """
    from .style import resolve_style, material_style, feature_style
    style = resolve_style(style)
    lw = feature_style(style, "polygon").get("linewidth", 1.0)
    seen_labels = set()
    for i, poly in enumerate(polygons):
        geom = poly['polygon'] if isinstance(poly, dict) else poly
        mat_idx = poly.get('mat_id') if isinstance(poly, dict) else i
        if mat_idx is None:
            mat_idx = i
        ms = material_style(style, mat_idx)
        color = ms["color"]

        # Legend label: material name if available, once per material.
        mat_name = None
        if materials and 0 <= mat_idx < len(materials):
            item = materials[mat_idx]
            mat_name = item.get('name') if isinstance(item, dict) else item
        label = mat_name if mat_name else f'Material {mat_idx + 1}'
        legend_label = label if label not in seen_labels else None
        if legend_label:
            seen_labels.add(label)

        xs, ys = geom.exterior.xy
        ax.fill(xs, ys, color=color, alpha=ms.get("alpha", 0.6),
                hatch=ms.get("hatch"), label=legend_label, gid=label)
        ax.plot(xs, ys, color=color, linewidth=lw, gid=label)
        # Draw any interior rings (holes) as outlines.
        for ring in geom.interiors:
            rx, ry = ring.xy
            ax.plot(rx, ry, color=color, linewidth=lw, linestyle='--', gid=label)

        if labels:
            c = geom.representative_point()
            ax.text(
                c.x, c.y, str(mat_idx + 1),
                fontsize=7, color='black', fontfamily='monospace',
                ha='center', va='center',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                          edgecolor=color, linewidth=0.8),
                zorder=10,
            )


def plot_max_depth(ax, profile_lines, max_depth, style=None):
    """
    Plots a horizontal line representing the maximum depth limit with hash marks.

    Parameters:
        ax: matplotlib Axes object
        profile_lines: List of profile line dicts, each with 'coords' key containing coordinate tuples
        max_depth: Maximum allowed depth for analysis

    Returns:
        None
    """
    if max_depth is None:
        return
    from .style import resolve_style, feature_style
    fs = feature_style(resolve_style(style), "max_depth")
    color = fs.get("color", "black")
    lw = fs.get("linewidth", 1.5)
    ls = fs.get("linestyle", "-")
    x_vals = [x for line in profile_lines for x, _ in line['coords']]
    x_min = min(x_vals)
    x_max = max(x_vals)
    ax.hlines(max_depth, x_min, x_max, colors=color, linewidth=lw, linestyles=ls,
              label='Max Depth', gid='MAX_DEPTH')

    x_diff = x_max - x_min
    spacing = x_diff / 100
    length = x_diff / 80

    angle_rad = np.radians(60)
    dx = length * np.cos(angle_rad)
    dy = length * np.sin(angle_rad)
    x_hashes = np.arange(x_min, x_max, spacing)[1:]
    for x in x_hashes:
        ax.plot([x, x - dx], [max_depth, max_depth - dy], color=color, linewidth=1, gid='MAX_DEPTH')


def plot_domain_base(ax, domain_polygon, label='Max Depth', style=None):
    """Draw the domain's lower boundary as a hatched base line (the polygon
    analog of plot_max_depth). Works for a flat bottom — reproducing the old
    horizontal hatched 'Max Depth' line — and for an irregular/sloping bottom,
    where the hatch marks follow the base."""
    if domain_polygon is None:
        return
    from .style import resolve_style, feature_style
    fs = feature_style(resolve_style(style), "max_depth")
    color = fs.get("color", "black")
    lw = fs.get("linewidth", 1.5)
    ls = fs.get("linestyle", "-")
    base = domain_lower_envelope(domain_polygon)
    if len(base) < 2:
        return
    bx = np.array([p[0] for p in base])
    by = np.array([p[1] for p in base])
    ax.plot(bx, by, color=color, linewidth=lw, linestyle=ls, label=label, gid='MAX_DEPTH')

    x_min, x_max = bx[0], bx[-1]
    x_diff = x_max - x_min
    if x_diff <= 0:
        return
    spacing = x_diff / 100
    length = x_diff / 80
    angle_rad = np.radians(60)
    dx = length * np.cos(angle_rad)
    dy = length * np.sin(angle_rad)
    for x in np.arange(x_min, x_max, spacing)[1:]:
        y = np.interp(x, bx, by)
        ax.plot([x, x - dx], [y, y - dy], color=color, linewidth=1, gid='MAX_DEPTH')


def plot_base_geometry(ax, slope_data, labels=False, style=None):
    """Plot the base geometry on an Axes: profile lines + max-depth line when the
    input uses profile lines, or filled material polygons + a hatched domain base
    when it uses polygons. Shared by plot_inputs and all the solution/search plots
    so both geometry types render consistently everywhere. `style` (see
    xslope.style) styles the material/profile/polygon layers; None → defaults."""
    if slope_data.get('profile_lines'):
        plot_profile_lines(ax, slope_data['profile_lines'],
                           materials=slope_data.get('materials'), labels=labels,
                           style=style)
        plot_max_depth(ax, slope_data['profile_lines'], slope_data['max_depth'], style=style)
    elif slope_data.get('polygons'):
        plot_polygons_on_ax(ax, slope_data['polygons'],
                            materials=slope_data.get('materials'), style=style)
        plot_domain_base(ax, slope_data.get('domain_polygon'), style=style)


def plot_coordinate_labels(ax, slope_data, fontsize=7, arrows=False, style=None):
    """Annotate the geometry's vertex coordinates, verification-manual style.

    Labels every unique vertex of the profile lines (or material-zone polygon
    exteriors) with "(x, y)". Vertices shared by several lines — ground points
    repeated on each layer's profile line, coincident polygon corners — are
    labelled once. Values print with %g so integers stay clean.

    Parameters:
        ax: matplotlib Axes object
        slope_data: dict with 'profile_lines' or 'polygons'
        fontsize: label font size in points (default 7)
        arrows: if True, tie each label to its vertex with a thin gray leader
            line and allow labels to be pushed well clear of dense clusters;
            if False (default) labels stay adjacent to their vertices and no
            leaders are drawn
        style: reserved for future style-sheet control (unused)
    """
    seen = set()
    points = []

    def _add(x, y):
        key = (round(float(x), 6), round(float(y), 6))
        if key not in seen:
            seen.add(key)
            points.append(key)

    if slope_data.get('profile_lines'):
        for line in slope_data['profile_lines']:
            for x, y in line['coords']:
                _add(x, y)
        md = slope_data.get('max_depth')
        if md is not None and slope_data['profile_lines']:
            xs = [x for line in slope_data['profile_lines'] for x, _ in line['coords']]
            _add(min(xs), md)
            _add(max(xs), md)
    elif slope_data.get('polygons'):
        for poly in slope_data['polygons']:
            coords = list(poly['polygon'].exterior.coords)
            if len(coords) > 1 and coords[0] == coords[-1]:
                coords = coords[:-1]
            for x, y in coords:
                _add(x, y)

    if not points:
        return
    # Placement is collision-aware with leader lines. Label boxes are measured
    # manually (text metrics + transData) rather than via the annotation
    # extent API, which reports stale positions for artists that have not
    # been drawn. Each label walks a ring of candidate offsets at increasing
    # radius until its box clears every already-placed label, and a thin gray
    # leader ties displaced labels back to their vertex so dense clusters (a
    # dam crest with several nearly coincident corners) stay attributable.
    # Right/left-edge points only get inward candidates.
    import math
    from matplotlib.font_manager import FontProperties
    from matplotlib.transforms import Bbox

    xs = [p[0] for p in points]
    x_min, x_max = min(xs), max(xs)
    span = max(x_max - x_min, 1e-9)
    fig = ax.figure
    renderer = None
    try:
        fig.canvas.draw()          # settle limits; obtain a renderer
        renderer = fig.canvas.get_renderer()
    except Exception:
        pass

    prop = FontProperties(size=fontsize)
    scale = fig.dpi / 72.0         # offset points -> pixels

    def candidates(near_left, near_right):
        angles = [55, 125, 20, 160, 90, 35, 145, -35, -145, 70, 110, -90,
                  0, 180, -15, -165]
        if near_right:
            angles = [a for a in angles if math.cos(math.radians(a)) < 0.3]
        if near_left:
            angles = [a for a in angles if math.cos(math.radians(a)) > -0.3]
        radii = (9, 20, 32, 46, 62, 80, 100, 122) if arrows else (6, 12, 20)
        for r in radii:
            for a in angles:
                yield (r * math.cos(math.radians(a)), r * math.sin(math.radians(a)))

    placed = []
    for x, y in sorted(points, key=lambda p: (p[0], -p[1])):
        near_right = (x_max - x) < 0.02 * span
        near_left = (x - x_min) < 0.02 * span
        label = f"({x:g}, {y:g})"

        dx_best, dy_best, ha_best, va_best, bb_best, n_best = 6, 6, "left", "bottom", None, None
        if renderer is not None:
            w, h, _ = renderer.get_text_width_height_descent(label, prop, False)
            px, py = ax.transData.transform((x, y))
            for dx, dy in candidates(near_left, near_right):
                ha = "left" if dx > 2 else ("right" if dx < -2 else "center")
                va = "bottom" if dy > 2 else ("top" if dy < -2 else "center")
                ox, oy = px + dx * scale, py + dy * scale
                x0 = ox - (w if ha == "right" else w / 2 if ha == "center" else 0)
                y0 = oy - (h if va == "top" else h / 2 if va == "center" else 0)
                bb = Bbox.from_bounds(x0 - 2, y0 - 2, w + 4, h + 4)
                n = sum(bb.overlaps(pb) for pb in placed)
                if n_best is None or n < n_best:
                    dx_best, dy_best, ha_best, va_best, bb_best, n_best = dx, dy, ha, va, bb, n
                if n == 0:
                    break
            if bb_best is not None:
                placed.append(bb_best)

        kwargs = {}
        # Leaders only where a collision actually displaced the label beyond
        # the near ring — labels sitting beside their vertex need no arrow.
        displaced = (dx_best ** 2 + dy_best ** 2) ** 0.5 > 13
        if arrows and displaced:
            kwargs['arrowprops'] = dict(arrowstyle="-", color="0.45",
                                        linewidth=0.5, shrinkA=1, shrinkB=1)
        ax.annotate(label, (x, y), textcoords="offset points",
                    xytext=(dx_best, dy_best), fontsize=fontsize, color="black",
                    ha=ha_best, va=va_best,
                    bbox=dict(boxstyle="round,pad=0.15", facecolor="white",
                              edgecolor="none", alpha=0.8),
                    zorder=6, gid="COORD_LABEL", **kwargs)


def plot_failure_surface(ax, failure_surface):
    """
    Plots the failure surface as a black line.

    Parameters:
        ax: matplotlib Axes object
        failure_surface: Shapely LineString representing the failure surface

    Returns:
        None
    """
    if failure_surface:
        x_clip, y_clip = zip(*failure_surface.coords)
        ax.plot(x_clip, y_clip, 'k-', linewidth=2, label="Failure Surface", gid="FAILURE_SURFACE")

def plot_slices(ax, slice_df, fill=True):
    """
    Plots the slices used in the analysis.

    Parameters:
        ax: matplotlib Axes object
        slice_df: DataFrame containing slice data
        fill: Boolean indicating whether to fill the slices with color

    Returns:
        None
    """
    if slice_df is not None:
        for _, row in slice_df.iterrows():
            if fill:
                xs = [row['x_l'], row['x_l'], row['x_r'], row['x_r'], row['x_l']]
                ys = [row['y_lb'], row['y_lt'], row['y_rt'], row['y_rb'], row['y_lb']]
                ax.plot(xs, ys, 'r-', gid='SLICES')
                ax.fill(xs, ys, color='red', alpha=0.1, gid='SLICES')
            else:
                ax.plot([row['x_l'], row['x_l']], [row['y_lb'], row['y_lt']], 'k-', linewidth=0.5, gid='SLICES')
                ax.plot([row['x_r'], row['x_r']], [row['y_rb'], row['y_rt']], 'k-', linewidth=0.5, gid='SLICES')

def plot_slice_numbers(ax, slice_df):
    """
    Plots the slice number in the middle of each slice at the middle height.
    Numbers are 1-indexed.

    Parameters:
        ax: matplotlib Axes object
        slice_df: DataFrame containing slice data

    Returns:
        None
    """
    if slice_df is not None:
        for _, row in slice_df.iterrows():
            # Calculate middle x-coordinate of the slice
            x_middle = row['x_c']
            
            # Calculate middle height of the slice
            y_middle = (row['y_cb'] + row['y_ct']) / 2
            
            # Plot the slice number (1-indexed)
            slice_number = int(row['slice #'])
            ax.text(x_middle, y_middle, str(slice_number), 
                   ha='center', va='center', fontsize=8, fontweight='bold',
                   bbox=dict(boxstyle="round,pad=0.2", facecolor='white', alpha=0.8))

def plot_piezo_line(ax, slope_data, style=None):
    """
    Plots the piezometric line(s) with markers at their midpoints.

    Parameters:
        ax: matplotlib Axes object
        data: Dictionary containing plot data with 'piezo_line' and optionally 'piezo_line2'
        style: optional style sheet (see xslope.style); None → defaults.

    Returns:
        None
    """
    from .style import resolve_style, feature_style
    style = resolve_style(style)

    def _plot_touching_v_marker(ax, x, y, color, markersize=8, extra_gap_points=0.0):
        """
        Place an inverted triangle marker so its tip visually touches the line at (x, y).
        We do this in display coordinates (points/pixels) so it scales consistently.
        """
        from matplotlib.markers import MarkerStyle
        from matplotlib.transforms import offset_copy

        # Compute the distance (in "marker units") from the marker origin to the tip.
        # Matplotlib scales marker vertices by `markersize` (in points) for Line2D.
        ms = MarkerStyle("v")
        path = ms.get_path().transformed(ms.get_transform())
        verts = np.asarray(path.vertices)
        min_y = float(verts[:, 1].min())  # tip is the lowest y
        tip_offset_points = (-min_y) * float(markersize) + float(extra_gap_points)

        # Offset the marker center upward in point units so the tip lands at (x, y).
        trans = offset_copy(ax.transData, fig=ax.figure, x=0.0, y=tip_offset_points, units="points")
        ax.plot([x], [y], marker="v", color=color, markersize=markersize, linestyle="None", transform=trans)

    def plot_single_piezo_line(ax, piezo_line, color, label, linewidth=2, linestyle='-'):
        """Internal function to plot a single piezometric line"""
        if not piezo_line:
            return

        piezo_xs, piezo_ys = zip(*piezo_line)
        ax.plot(piezo_xs, piezo_ys, color=color, linewidth=linewidth,
                linestyle=linestyle, label=label, gid='PIEZO')
        
        # Find middle x-coordinate and corresponding y value
        if len(piezo_xs) > 1:
            # Sort by x to ensure monotonic input for interpolation
            pairs = sorted(zip(piezo_xs, piezo_ys), key=lambda p: p[0])
            sx, sy = zip(*pairs)
            x_min, x_max = min(sx), max(sx)
            mid_x = (x_min + x_max) / 2
            mid_y = float(np.interp(mid_x, sx, sy))
            # Slight negative gap so the marker visually "touches" the line (not floating above it)
            _plot_touching_v_marker(ax, mid_x, mid_y, color=color, markersize=8, extra_gap_points=2.0)
    
    # Plot both piezometric lines
    f1 = feature_style(style, "piezo_line")
    f2 = feature_style(style, "piezo_line2")
    plot_single_piezo_line(ax, slope_data.get('piezo_line'), f1.get('color', 'b'),
                           "Piezometric Line", f1.get('linewidth', 2), f1.get('linestyle', '-'))
    plot_single_piezo_line(ax, slope_data.get('piezo_line2'), f2.get('color', 'skyblue'),
                           "Piezometric Line 2", f2.get('linewidth', 2), f2.get('linestyle', '-'))

def plot_seepage_bc_lines(ax, slope_data, style=None):
    """
    Plots seep boundary-condition lines for seep-only workflows.

    Plots the primary seepage BCs (seepage_bc) and, if present, the second set
    (seepage_bc2) with different colors to distinguish them. `style` (see
    xslope.style) styles BC set 1; set 2 keeps its distinct default colors so a
    rapid-drawdown pair stays visually separable.
    """
    from .style import resolve_style, feature_style
    style = resolve_style(style)
    def _plot_touching_v_marker(ax, x, y, color, markersize=8, extra_gap_points=2.0):
        """Place an inverted triangle so its tip visually sits on the line at (x, y)."""
        from matplotlib.markers import MarkerStyle
        from matplotlib.transforms import offset_copy

        ms = MarkerStyle("v")
        path = ms.get_path().transformed(ms.get_transform())
        verts = np.asarray(path.vertices)
        min_y = float(verts[:, 1].min())
        tip_offset_points = (-min_y) * float(markersize) + float(extra_gap_points)
        trans = offset_copy(ax.transData, fig=ax.figure, x=0.0, y=tip_offset_points, units="points")
        ax.plot([x], [y], marker="v", color=color, markersize=markersize, linestyle="None", transform=trans)

    def _plot_one_bc_set(ax, seepage_bc, geom_width, x_min_geom, x_max_geom,
                         head_line_color, water_level_color, exit_face_color, label_suffix="",
                         head_lw=3, head_ls="--", water_lw=2, exit_lw=3, exit_ls="--",
                         flux_color="darkgreen", flux_lw=3, flux_ls="-."):
        """Plot a single set of seepage boundary conditions."""
        specified_heads = seepage_bc.get("specified_heads") or []
        specified_fluxes = seepage_bc.get("specified_fluxes") or []
        exit_face = seepage_bc.get("exit_face") or []

        for i, sh in enumerate(specified_heads):
            coords = sh.get("coords") or []
            if len(coords) < 2:
                continue

            xs, ys = zip(*coords)
            ax.plot(
                xs, ys,
                color=head_line_color, linewidth=head_lw, linestyle=head_ls,
                label=f"Specified Head Line{label_suffix}" if i == 0 else "",
            )

            head_val = sh.get("head", None)
            if head_val is None:
                continue

            if isinstance(head_val, (list, tuple, np.ndarray)):
                if len(head_val) != len(coords):
                    continue
                heads = [float(h) for h in head_val]
            else:
                try:
                    head_scalar = float(head_val)
                except (TypeError, ValueError):
                    continue
                heads = [head_scalar] * len(coords)

            tol = 1e-6
            is_vertical = (max(xs) - min(xs)) <= tol
            if is_vertical:
                x0 = float(xs[0])
                y_head = float(heads[0])
                seg_len = 0.04 * geom_width
                gap = 0.01 * geom_width
                is_right = x0 >= 0.5 * (x_min_geom + x_max_geom)
                if is_right:
                    wl_xs = [x0 + gap, x0 + gap + seg_len]
                else:
                    wl_xs = [x0 - gap - seg_len, x0 - gap]
                wl_ys = [y_head, y_head]
            else:
                wl_xs = list(xs)
                wl_ys = heads

            ax.plot(
                wl_xs, wl_ys,
                color=water_level_color, linewidth=water_lw, linestyle="-",
                label=f"Specified Head Water Level{label_suffix}" if i == 0 else "",
            )

            if len(wl_xs) > 1:
                try:
                    pairs = sorted(zip(wl_xs, wl_ys), key=lambda p: p[0])
                    sx, sy = zip(*pairs)
                    mid_x = 0.5 * (min(sx) + max(sx))
                    mid_y = float(np.interp(mid_x, sx, sy))
                    _plot_touching_v_marker(ax, mid_x, mid_y, color=water_level_color, markersize=8, extra_gap_points=2.0)
                except Exception:
                    pass

        for i, sf in enumerate(specified_fluxes):
            coords = sf.get("coords") or []
            if len(coords) < 2:
                continue
            fx_xs, fx_ys = zip(*coords)
            ax.plot(
                fx_xs, fx_ys,
                color=flux_color, linewidth=flux_lw, linestyle=flux_ls,
                label=f"Specified Flux{label_suffix}" if i == 0 else "",
            )
            try:
                q_val = float(sf.get("flux"))
            except (TypeError, ValueError):
                continue
            mid = len(coords) // 2
            if len(coords) % 2 == 0:
                lx = 0.5 * (coords[mid - 1][0] + coords[mid][0])
                ly = 0.5 * (coords[mid - 1][1] + coords[mid][1])
            else:
                lx, ly = coords[mid]
            ax.annotate(
                f"q = {q_val:g}", xy=(lx, ly), xytext=(0, 6),
                textcoords="offset points", ha="center", va="bottom",
                color=flux_color, fontsize=8,
            )

        if len(exit_face) >= 2:
            ex_xs, ex_ys = zip(*exit_face)
            ax.plot(
                ex_xs, ex_ys,
                color=exit_face_color, linewidth=exit_lw, linestyle=exit_ls,
                label=f"Exit Face{label_suffix}",
            )

    # Geometry x-extent (used for vertical-head-line derived segment length / side)
    x_vals = []
    for line in slope_data.get("profile_lines", []):
        try:
            xs_line, _ = zip(*line['coords'])
            x_vals.extend(xs_line)
        except Exception:
            pass
    gs = slope_data.get("ground_surface", None)
    if not x_vals and gs is not None and hasattr(gs, "coords"):
        x_vals.extend([x for x, _ in gs.coords])
    x_min_geom = min(x_vals) if x_vals else 0.0
    x_max_geom = max(x_vals) if x_vals else 1.0
    geom_width = max(1e-9, x_max_geom - x_min_geom)

    # Plot primary BCs
    seepage_bc = slope_data.get("seepage_bc") or {}
    has_bc2 = slope_data.get("has_seepage_bc2", False)
    label_suffix = " (BC 1)" if has_bc2 else ""
    fh = feature_style(style, "seep_bc")
    fw = feature_style(style, "seep_water_level")
    fe = feature_style(style, "seep_exit_face")
    ff = feature_style(style, "seep_flux")
    _plot_one_bc_set(ax, seepage_bc, geom_width, x_min_geom, x_max_geom,
                     head_line_color=fh.get("color", "darkblue"),
                     water_level_color=fw.get("color", "lightskyblue"),
                     exit_face_color=fe.get("color", "red"), label_suffix=label_suffix,
                     head_lw=fh.get("linewidth", 3), head_ls=fh.get("linestyle", "--"),
                     water_lw=fw.get("linewidth", 2),
                     exit_lw=fe.get("linewidth", 3), exit_ls=fe.get("linestyle", "--"),
                     flux_color=ff.get("color", "darkgreen"),
                     flux_lw=ff.get("linewidth", 3), flux_ls=ff.get("linestyle", "-."))

    # Plot second set of BCs if present
    if has_bc2:
        seepage_bc2 = slope_data.get("seepage_bc2") or {}
        _plot_one_bc_set(ax, seepage_bc2, geom_width, x_min_geom, x_max_geom,
                         head_line_color="steelblue", water_level_color="powderblue",
                         exit_face_color="orangered", label_suffix=" (BC 2)",
                         flux_color="seagreen")

def plot_tcrack_surface(ax, slope_data, style=None):
    """
    Plots the tension crack surface as a thin dashed red line, clipped to max_depth.

    Parameters:
        ax: matplotlib Axes object
        slope_data: Dictionary containing tcrack_surface and max_depth

    Returns:
        None
    """
    tcrack_surface = slope_data.get('tcrack_surface')
    if tcrack_surface is None:
        return

    from .style import resolve_style, feature_style
    fs = feature_style(resolve_style(style), "tcrack")
    color = fs.get('color', 'red')
    linestyle = fs.get('linestyle', ':')
    linewidth = fs.get('linewidth', 1.5)

    max_depth = slope_data.get('max_depth')
    if max_depth is None:
        # No clipping needed
        x_vals, y_vals = tcrack_surface.xy
        ax.plot(x_vals, y_vals, linestyle=linestyle, color=color, linewidth=linewidth, label='Tension Crack Depth', gid='TENSION_CRACK')
        return

    # Get coordinates and clip to max_depth with interpolation
    coords = list(tcrack_surface.coords)
    x_clipped = []
    y_clipped = []

    for i in range(len(coords)):
        x1, y1 = coords[i]

        if y1 >= max_depth:
            # Point is above max_depth, include it
            x_clipped.append(x1)
            y_clipped.append(y1)

        # Check if segment crosses max_depth (need to interpolate)
        if i < len(coords) - 1:
            x2, y2 = coords[i + 1]
            # Check if segment crosses max_depth
            if (y1 < max_depth and y2 >= max_depth) or (y1 >= max_depth and y2 < max_depth):
                # Interpolate to find crossing point
                t = (max_depth - y1) / (y2 - y1)
                x_cross = x1 + t * (x2 - x1)
                x_clipped.append(x_cross)
                y_clipped.append(max_depth)

    if x_clipped:
        ax.plot(x_clipped, y_clipped, linestyle=linestyle, color=color, linewidth=linewidth, label='Tension Crack Depth', gid='TENSION_CRACK')


def plot_tcrack_water_force(ax, slice_df, slope_data):
    """
    Plots the triangular water pressure distribution on the tension crack face.

    The water in the tension crack creates a triangular pressure distribution
    acting horizontally on the side of the top slice. Pressure is zero at the
    water surface and maximum (gamma_w * water_depth) at the bottom.
    The triangle is drawn on the outside of the slice, with arrows pointing
    toward the slice to show force direction. The base of the triangle is
    scaled to equal the water depth.

    Parameters:
        ax: matplotlib Axes object
        slice_df: DataFrame containing slice data with 't' and 'y_t' columns
        slope_data: Dictionary containing slope data including tcrack_water

    Returns:
        None
    """
    tcrack_water = slope_data.get('tcrack_water', 0)
    if tcrack_water <= 0:
        return

    # Find the slice with the tension crack force
    t_forces = slice_df['t'].abs()
    if t_forces.max() == 0:
        return

    # Get the slice with the tension crack force
    tcrack_slice_idx = t_forces.idxmax()
    tcrack_slice = slice_df.loc[tcrack_slice_idx]

    t_force = tcrack_slice['t']
    y_rb = tcrack_slice['y_rb']
    y_rt = tcrack_slice['y_rt']

    # Determine if right-facing or left-facing based on sign of t
    # Negative t means right-facing (force acts to the right, on left side of first slice)
    # Positive t means left-facing (force acts to the left, on right side of last slice)
    right_facing = t_force < 0

    if right_facing:
        # Water on left side of slice, triangle extends left (outside), arrows point right (into slice)
        x_base = tcrack_slice['x_l']
        triangle_direction = -1  # triangle extends left (outside the slice)
        arrow_direction = 1      # arrows point right (into the slice)
    else:
        # Water on right side of slice, triangle extends right (outside), arrows point left (into slice)
        x_base = tcrack_slice['x_r']
        triangle_direction = 1   # triangle extends right (outside the slice)
        arrow_direction = -1     # arrows point left (into the slice)

    # Water surface is at ground level, bottom of water is at y_rb
    y_water_top = y_rt  # top of water (at ground surface)
    y_water_bottom = y_rb  # bottom of water (at failure surface)
    water_depth = y_water_top - y_water_bottom

    if water_depth <= 0:
        return

    # Scale so that the base of the triangle equals the water depth
    max_length = tcrack_water  # base of triangle = water depth

    # Arrow head dimensions (same style as distributed loads)
    head_length = max_length / 8
    head_width = head_length * 0.8

    # Draw triangular pressure distribution (on outside of slice)
    num_arrows = 5
    y_positions = np.linspace(y_water_bottom, y_water_top, num_arrows + 1)[:-1]

    for y_pos in y_positions:
        # Arrow length proportional to depth (0 at top, max_length at bottom)
        depth_from_surface = y_water_top - y_pos
        arrow_length = max_length * (depth_from_surface / water_depth)

        if arrow_length < 0.1:
            continue  # Skip very short arrows

        # Arrow starts from outside (triangle edge) and points toward slice
        x_start = x_base + arrow_length * triangle_direction
        dx = -arrow_length * triangle_direction  # direction toward slice

        # Draw arrow using same style as distributed loads
        if head_length > arrow_length:
            # Draw a simple line without arrowhead for short arrows
            ax.plot([x_start, x_base], [y_pos, y_pos],
                   color='blue', linewidth=2, alpha=0.7)
        else:
            ax.arrow(x_start, y_pos, dx, 0,
                    head_width=head_width, head_length=head_length,
                    fc='blue', ec='blue', alpha=0.7,
                    length_includes_head=True)

    # Draw the triangular outline (pressure diagram) on outside of slice
    triangle_x = [x_base, x_base + max_length * triangle_direction, x_base]
    triangle_y = [y_water_top, y_water_bottom, y_water_bottom]
    ax.fill(triangle_x, triangle_y, color='lightblue', alpha=0.3, edgecolor='blue', linewidth=1)


def plot_dloads(ax, slope_data, style=None):
    """
    Plots distributed loads as arrows along the surface.
    """
    from .style import resolve_style, feature_style
    style = resolve_style(style)
    gamma_w = slope_data['gamma_water']
    ground_surface = slope_data['ground_surface']

    def plot_single_dload_set(ax, dloads, color, label, linewidth=1.5):
        """Internal function to plot a single set of distributed loads"""
        if not dloads:
            return
            
        # find the max horizontal length of the ground surface
        max_horizontal_length_ground = 0
        for pt in ground_surface.coords:
            max_horizontal_length_ground = max(max_horizontal_length_ground, pt[0])

        arrow_spacing = max_horizontal_length_ground / 60

        # find the max dload value
        max_dload = 0
        for line in dloads:
            max_dload = max(max_dload, max(pt['Normal'] for pt in line))

        arrow_height = max_dload / gamma_w
        # Head size: nominally proportional to the arrow height, but capped by
        # the arrow spacing so the heads of tall arrows (deep ponded water)
        # cannot grow into each other and merge into a solid band.
        head_length = min(arrow_height / 12, 0.75 * arrow_spacing)
        head_width = 0.8 * head_length
        
        # Find the maximum load value for scaling
        max_load = 0
        for line in dloads:
            max_load = max(max_load, max(pt['Normal'] for pt in line))
        
        for line in dloads:
            if len(line) < 2:
                continue
                
            xs = [pt['X'] for pt in line]
            ys = [pt['Y'] for pt in line]
            ns = [pt['Normal'] for pt in line]
            
            # Process line segments
            for i in range(len(line) - 1):
                x1, y1, n1 = xs[i], ys[i], ns[i]
                x2, y2, n2 = xs[i+1], ys[i+1], ns[i+1]
                
                # Calculate segment direction (perpendicular to this segment)
                dx = x2 - x1
                dy = y2 - y1
                segment_length = np.sqrt(dx**2 + dy**2)
                
                if segment_length == 0:
                    continue
                    
                # Normalize the segment direction
                dx_norm = dx / segment_length
                dy_norm = dy / segment_length
                
                # Perpendicular direction (rotate 90 degrees CCW)
                perp_dx = -dy_norm
                perp_dy = dx_norm
                
                # Generate arrows along this segment
                dx_abs = abs(x2 - x1)
                num_arrows = max(1, int(round(dx_abs / arrow_spacing)))
                if dx_abs == 0:
                    t_values = np.array([0.0, 1.0])
                else:
                    t_values = np.linspace(0, 1, num_arrows + 1)
                
                # Store arrow top points for connecting line
                top_xs = []
                top_ys = []
                
                # Add start point if it's the first segment and load is zero
                if i == 0 and n1 == 0:
                    top_xs.append(x1)
                    top_ys.append(y1)
                
                for t in t_values:
                    # Interpolate position along segment
                    x = x1 + t * dx
                    y = y1 + t * dy
                    
                    # Interpolate load value
                    n = n1 + t * (n2 - n1)
                    
                    # Scale arrow height based on equivalent water depth
                    if max_load > 0:
                        water_depth = n / gamma_w
                        arrow_height = water_depth  # Direct water depth, not scaled relative to max
                    else:
                        arrow_height = 0
                    
                    # For very small arrows, just store surface point for connecting line
                    if arrow_height < 0.5:
                        top_xs.append(x)
                        top_ys.append(y)
                        continue
            
                    
                    # Calculate arrow start point (above surface)
                    arrow_start_x = x + perp_dx * arrow_height
                    arrow_start_y = y + perp_dy * arrow_height
                    
                    # Store points for connecting line
                    top_xs.append(arrow_start_x)
                    top_ys.append(arrow_start_y)
                    
                    # Draw arrow - extend all the way to surface point
                    arrow_length = np.sqrt((x - arrow_start_x)**2 + (y - arrow_start_y)**2)
                    if head_length > arrow_length:
                        # Draw a simple line without arrowhead
                        ax.plot([arrow_start_x, x], [arrow_start_y, y], 
                               color=color, linewidth=2, alpha=0.7, gid='DLOADS')
                    else:
                        # Draw arrow with head
                        ax.arrow(arrow_start_x, arrow_start_y, 
                                x - arrow_start_x, y - arrow_start_y,
                                head_width=head_width, head_length=head_length, 
                                fc=color, ec=color, alpha=0.7,
                                length_includes_head=True, gid='DLOADS')
                
                # Add end point if it's the last segment and load is zero
                if i == len(line) - 2 and n2 == 0:
                    top_xs.append(x2)
                    top_ys.append(y2)
                
                # Draw connecting line at arrow tops
                if top_xs:
                    ax.plot(top_xs, top_ys, color=color, linewidth=linewidth, alpha=0.8, gid='DLOADS')

            # Draw the surface line itself
            ax.plot(xs, ys, color=color, linewidth=linewidth, alpha=0.8, label=label, gid='DLOADS')

    df1 = feature_style(style, "dloads")
    df2 = feature_style(style, "dloads2")
    dloads = slope_data['dloads']
    dloads2 = slope_data.get('dloads2', [])
    plot_single_dload_set(ax, dloads, df1.get('color', 'purple'), 'Distributed Load',
                          df1.get('linewidth', 1.5))
    plot_single_dload_set(ax, dloads2, df2.get('color', 'orange'), 'Distributed Load 2',
                          df2.get('linewidth', 1.5))

def plot_circles(ax, slope_data, style=None):
    """
    Plots starting circles with center markers and arrows.

    Parameters:
        ax (matplotlib axis): The plotting axis
        slope_data (dict): Slope data dictionary containing circles

    Returns:
        None
    """
    from .style import resolve_style, feature_style
    fs = feature_style(resolve_style(style), "failure_surface")
    c_color = fs.get('color', 'red')
    c_ls = fs.get('linestyle', '--')
    c_lw = fs.get('linewidth', 1.5)

    circles = slope_data['circles']
    tcrack_depth = slope_data.get('tcrack_depth', 0)

    for i, circle in enumerate(circles):
        Xo = circle['Xo']
        Yo = circle['Yo']
        R = circle['R']
        # theta = np.linspace(0, 2 * np.pi, 100)
        # x_circle = Xo + R * np.cos(theta)
        # y_circle = Yo + R * np.sin(theta)
        # ax.plot(x_circle, y_circle, 'r--', label='Circle')

        # Plot the portion of the circle in the slope (clipped to tension crack if present)
        ground_surface = slope_data['ground_surface']
        success, result = generate_failure_surface(ground_surface, circular=True, circle=circle, tcrack_depth=tcrack_depth)
        if not success:
            print(f"Warning: Circle {i+1} (Xo={Xo:.2f}, Yo={Yo:.2f}, R={R:.2f}) could not be plotted: {result}")
            continue
        # result = (x_min, x_max, y_left, y_right, clipped_surface)
        x_min, x_max, y_left, y_right, clipped_surface = result
        if not isinstance(clipped_surface, LineString):
            clipped_surface = LineString(clipped_surface)
        x_clip, y_clip = zip(*clipped_surface.coords)
        ax.plot(x_clip, y_clip, color=c_color, linestyle=c_ls, linewidth=c_lw,
                label="Circle", gid='CIRCLES')

        # Center marker — annotation layer (add_artist), not data: a center far
        # above the section must not inflate the autoscaled view. It still
        # draws whenever the equal-aspect view reaches it, and clips otherwise.
        center_marker = Line2D([Xo], [Yo], marker='+', color=c_color,
                               linestyle='None', markersize=10, gid='CIRCLES')
        center_marker.set_in_layout(False)  # keep tight bbox from reserving
        ax.add_artist(center_marker)        # space for an off-view center

        # Arrow direction: point from center to midpoint of failure surface
        mid_idx = len(x_clip) // 2
        x_mid = x_clip[mid_idx]
        y_mid = y_clip[mid_idx]

        dx = x_mid - Xo
        dy = y_mid - Yo

        # Normalize direction vector
        length = np.hypot(dx, dy)
        if length != 0:
            dx /= length
            dy /= length

        # Draw arrow with pixel-based head size
        ann = ax.annotate('',
                    xy=(Xo + dx * R, Yo + dy * R),  # arrow tip
                    xytext=(Xo, Yo),                 # arrow start
                    arrowprops=dict(
                        arrowstyle='-|>',
                        color=c_color,
                        lw=1.0,            # shaft width in points
                        mutation_scale=20  # head size in points
                    ))
        # Annotations default to clip_on=False; with the center allowed outside
        # the view, the shaft must stop at the axes edge and the annotation's
        # full (unclipped) extent must stay out of tight-bbox layout math.
        ann.arrow_patch.set_clip_box(ax.bbox)
        ann.set_in_layout(False)

def plot_non_circ(ax, non_circ, style=None):
    """
    Plots a non-circular failure surface.

    Parameters:
        ax: matplotlib Axes object
        non_circ: List of coordinates representing the non-circular failure surface

    Returns:
        None
    """
    if not non_circ or len(non_circ) == 0:
        return
    from .style import resolve_style, feature_style
    fs = feature_style(resolve_style(style), "failure_surface")
    # Handle both dict format {'X': x, 'Y': y} and tuple format (x, y)
    if isinstance(non_circ[0], dict):
        xs = [p['X'] for p in non_circ]
        ys = [p['Y'] for p in non_circ]
    else:
        xs, ys = zip(*non_circ)
    ax.plot(xs, ys, color=fs.get('color', 'red'), linestyle=fs.get('linestyle', '--'),
            linewidth=fs.get('linewidth', 1.5), label='Non-Circular Surface')

def plot_lem_material_table(ax, materials, xloc=0.6, yloc=0.7):
    """
    Adds a limit equilibrium material properties table to the plot.

    Displays soil properties for limit equilibrium analysis including unit weight (γ),
    cohesion (c), friction angle (φ), and optionally dilation angle (d) and
    dilatancy angle (ψ). Supports both Mohr-Coulomb (mc) and constant-phi (cp) options.

    Parameters:
        ax: matplotlib Axes object to add the table to
        materials: List of material property dictionaries with keys:
            - 'name': Material name (str)
            - 'gamma': Unit weight (float)
            - 'option': Material model - 'mc' or 'cp' (str)
            - 'c': Cohesion for mc option (float)
            - 'phi': Friction angle for mc option (float)
            - 'cp': Constant phi for cp option (float)
            - 'r_elev': Reference elevation for cp option (float)
            - 'd': Dilation angle, optional (float)
            - 'psi': Dilatancy angle, optional (float)
        xloc: x-location of table bottom-left corner in axes coordinates (0-1, default: 0.6)
        yloc: y-location of table bottom-left corner in axes coordinates (0-1, default: 0.7)

    Returns:
        None
    """
    if not materials:
        return

    # Check if any materials have non-zero d and psi values
    has_d_psi = any(mat.get('d', 0) > 0 or mat.get('psi', 0) > 0 for mat in materials)

    # Check material options
    options = set(mat['option'] for mat in materials)

    # Decide column headers
    if options == {'mc'}:
        if has_d_psi:
            col_labels = ["Mat", "Name", "γ", "c", "φ", "d", "ψ"]
        else:
            col_labels = ["Mat", "Name", "γ", "c", "φ"]
    elif options == {'cp'}:
        if has_d_psi:
            col_labels = ["Mat", "Name", "γ", "cp", "rₑ", "d", "ψ"]
        else:
            col_labels = ["Mat", "Name", "γ", "cp", "rₑ"]
    else:
        if has_d_psi:
            col_labels = ["Mat", "Name", "γ", "c / cp", "φ / rₑ", "d", "ψ"]
        else:
            col_labels = ["Mat", "Name", "γ", "c / cp", "φ / rₑ"]

    # Build table rows
    table_data = []
    for idx, mat in enumerate(materials):
        name = mat['name']
        gamma = mat['gamma']
        option = mat['option']
        
        if option == 'mc':
            c = mat['c']
            phi = mat['phi']
            if has_d_psi:
                d = mat.get('d', 0)
                psi = mat.get('psi', 0)
                d_str = f"{d:.1f}" if d > 0 or psi > 0 else "-"
                psi_str = f"{psi:.1f}" if d > 0 or psi > 0 else "-"
                row = [idx+1, name, f"{gamma:.1f}", f"{c:.1f}", f"{phi:.1f}", d_str, psi_str]
            else:
                row = [idx+1, name, f"{gamma:.1f}", f"{c:.1f}", f"{phi:.1f}"]
        elif option == 'cp':
            cp = mat['cp']
            r_elev = mat['r_elev']
            if has_d_psi:
                d = mat.get('d', 0)
                psi = mat.get('psi', 0)
                d_str = f"{d:.1f}" if d > 0 or psi > 0 else "-"
                psi_str = f"{psi:.1f}" if d > 0 or psi > 0 else "-"
                row = [idx+1, name, f"{gamma:.1f}", f"{cp:.2f}", f"{r_elev:.1f}", d_str, psi_str]
            else:
                row = [idx+1, name, f"{gamma:.1f}", f"{cp:.2f}", f"{r_elev:.1f}"]
        else:
            if has_d_psi:
                d = mat.get('d', 0)
                psi = mat.get('psi', 0)
                d_str = f"{d:.1f}" if d > 0 or psi > 0 else "-"
                psi_str = f"{psi:.1f}" if d > 0 or psi > 0 else "-"
                row = [idx+1, name, f"{gamma:.1f}", "-", "-", d_str, psi_str]
            else:
                row = [idx+1, name, f"{gamma:.1f}", "-", "-"]
        table_data.append(row)

    # Adjust table width based on number of columns
    table_width = 0.25 if has_d_psi else 0.2

    # Choose table height based on number of materials (uniform across table types)
    num_rows = max(1, len(materials))
    table_height = 0.06 + 0.035 * num_rows  # header + per-row estimate
    table_height = min(0.35, table_height)  # cap to avoid overflows for many rows

    # Add the table
    table = ax.table(cellText=table_data,
                     colLabels=col_labels,
                     loc='upper right',
                     colLoc='center',
                     cellLoc='center',
                     bbox=[xloc, yloc, table_width, table_height])
    table.auto_set_font_size(False)
    table.set_fontsize(8)
    # Auto layout based on content (shared method for all table types)
    auto_size_table_to_content(ax, table, col_labels, table_data, table_width, table_height)

def plot_seep_material_table(ax, seep_data, xloc=0.6, yloc=0.7):
    """
    Adds a seep material properties table to the plot.

    Displays hydraulic properties for seep analysis including hydraulic conductivities
    (k₁, k₂), anisotropy angle, and unsaturated flow parameters (kr₀, h₀).

    Parameters:
        ax: matplotlib Axes object to add the table to
        seep_data: Dictionary containing seep material properties with keys:
            - 'k1_by_mat': List of primary hydraulic conductivity values (float)
            - 'k2_by_mat': List of secondary hydraulic conductivity values (float)
            - 'angle_by_mat': List of anisotropy angles in degrees (float)
            - 'kr0_by_mat': List of relative permeability at residual saturation (float)
            - 'h0_by_mat': List of pressure head parameters (float)
            - 'material_names': List of material names (str), optional
        xloc: x-location of table bottom-left corner in axes coordinates (0-1, default: 0.6)
        yloc: y-location of table bottom-left corner in axes coordinates (0-1, default: 0.7)

    Returns:
        None
    """
    k1_by_mat = seep_data.get("k1_by_mat")
    k2_by_mat = seep_data.get("k2_by_mat")
    angle_by_mat = seep_data.get("angle_by_mat")
    kr0_by_mat = seep_data.get("kr0_by_mat")
    h0_by_mat = seep_data.get("h0_by_mat")
    material_names = seep_data.get("material_names", [])
    if k1_by_mat is None or len(k1_by_mat) == 0:
        return
    col_labels = ["Mat", "Name", "k₁", "k₂", "Angle", "kr₀", "h₀"]
    table_data = []
    for idx in range(len(k1_by_mat)):
        k1 = k1_by_mat[idx]
        k2 = k2_by_mat[idx] if k2_by_mat is not None else 0.0
        angle = angle_by_mat[idx] if angle_by_mat is not None else 0.0
        kr0 = kr0_by_mat[idx] if kr0_by_mat is not None else 0.0
        h0 = h0_by_mat[idx] if h0_by_mat is not None else 0.0
        material_name = material_names[idx] if idx < len(material_names) else f"Material {idx+1}"
        row = [idx + 1, material_name, f"{k1:.3f}", f"{k2:.3f}", f"{angle:.1f}", f"{kr0:.4f}", f"{h0:.2f}"]
        table_data.append(row)
    # Dimensions
    num_rows = max(1, len(k1_by_mat))
    table_width = 0.45
    table_height = 0.10 + 0.06 * num_rows
    table_height = min(0.50, table_height)
    table = ax.table(cellText=table_data,
                     colLabels=col_labels,
                     loc='upper right',
                     colLoc='center',
                     cellLoc='center',
                     bbox=[xloc, yloc, table_width, table_height])
    table.auto_set_font_size(False)
    table.set_fontsize(8)
    auto_size_table_to_content(ax, table, col_labels, table_data, table_width, table_height)

def plot_fem_material_table(ax, fem_data, xloc=0.6, yloc=0.7, width=0.6, height=None):
    """
    Adds a finite element material properties table to the plot.

    Displays material properties for FEM analysis including unit weight (γ), cohesion (c),
    friction angle (φ), Young's modulus (E), and Poisson's ratio (ν).

    Parameters:
        ax: matplotlib Axes object to add the table to
        fem_data: Dictionary containing FEM material properties with keys:
            - 'c_by_mat': List of cohesion values (float)
            - 'phi_by_mat': List of friction angle values in degrees (float)
            - 'E_by_mat': List of Young's modulus values (float)
            - 'nu_by_mat': List of Poisson's ratio values (float)
            - 'gamma_by_mat': List of unit weight values (float)
            - 'material_names': List of material names (str), optional
        xloc: x-location of table bottom-left corner in axes coordinates (0-1, default: 0.6)
        yloc: y-location of table bottom-left corner in axes coordinates (0-1, default: 0.7)
        width: Table width in axes coordinates (0-1, default: 0.6)
        height: Table height in axes coordinates (0-1, default: auto-calculated)

    Returns:
        None
    """
    c_by_mat = fem_data.get("c_by_mat")
    phi_by_mat = fem_data.get("phi_by_mat")
    E_by_mat = fem_data.get("E_by_mat")
    nu_by_mat = fem_data.get("nu_by_mat")
    gamma_by_mat = fem_data.get("gamma_by_mat")
    material_names = fem_data.get("material_names", [])
    if c_by_mat is None or len(c_by_mat) == 0:
        return
    col_labels = ["Mat", "Name", "γ", "c", "φ", "E", "ν"]
    table_data = []
    for idx in range(len(c_by_mat)):
        c = c_by_mat[idx]
        phi = phi_by_mat[idx] if phi_by_mat is not None else 0.0
        E = E_by_mat[idx] if E_by_mat is not None else 0.0
        nu = nu_by_mat[idx] if nu_by_mat is not None else 0.0
        gamma = gamma_by_mat[idx] if gamma_by_mat is not None else 0.0
        material_name = material_names[idx] if idx < len(material_names) else f"Material {idx+1}"
        row = [idx + 1, material_name, f"{gamma:.1f}", f"{c:.1f}", f"{phi:.1f}", f"{E:.0f}", f"{nu:.2f}"]
        table_data.append(row)
    if height is None:
        num_rows = max(1, len(c_by_mat))
        height = 0.06 + 0.035 * num_rows
        height = min(0.32, height)
    table = ax.table(cellText=table_data,
                     colLabels=col_labels,
                     loc='upper right',
                     colLoc='center',
                     cellLoc='center',
                     bbox=[xloc, yloc, width, height])
    table.auto_set_font_size(False)
    table.set_fontsize(8)
    auto_size_table_to_content(ax, table, col_labels, table_data, width, height)
def auto_size_table_to_content(ax, table, col_labels, table_data, table_width, table_height, min_row_frac=0.02, row_pad=1.35, col_min_frac=0.08, col_max_frac=0.15):
    """
    Automatically adjusts table column widths and row heights based on content.

    Measures text extents using the matplotlib renderer and sets column widths proportional
    to content while enforcing minimum and maximum constraints. The "Name" column gets
    more space (18-30%) while numeric columns are constrained to prevent excessive whitespace.
    Row heights are uniform and based on the tallest content in each row.

    Parameters:
        ax: matplotlib Axes object containing the table
        table: matplotlib Table object to be resized
        col_labels: List of column header labels (str)
        table_data: List of lists containing table cell data
        table_width: Total table width in axes coordinates (0-1)
        table_height: Total table height in axes coordinates (0-1)
        min_row_frac: Minimum row height as fraction of axes height (default: 0.02)
        row_pad: Padding factor applied to measured row heights (default: 1.35)
        col_min_frac: Minimum column width as fraction of table width for numeric columns (default: 0.08)
        col_max_frac: Maximum column width as fraction of table width for numeric columns (default: 0.15)

    Returns:
        None

    Notes:
        - The "Name" column is automatically left-aligned and gets 18-30% of table width
        - Numeric columns are center-aligned and constrained to 8-15% of table width
        - All three material table types (LEM, SEEP, FEM) use the same sizing parameters
    """
    # Force draw to get a valid renderer
    try:
        ax.figure.canvas.draw()
        renderer = ax.figure.canvas.get_renderer()
    except Exception:
        renderer = None

    ncols = len(col_labels)
    nrows = len(table_data) + 1  # include header
    # Measure text widths per column in pixels
    widths_px = [1.0] * ncols
    if renderer is not None:
        for c in range(ncols):
            max_w = 1.0
            for r in range(nrows):
                cell = table[(r, c)]
                text = cell.get_text()
                try:
                    bbox = text.get_window_extent(renderer=renderer)
                    max_w = max(max_w, bbox.width)
                except Exception:
                    pass
            widths_px[c] = max_w
    total_w = sum(widths_px) if sum(widths_px) > 0 else float(ncols)
    col_fracs = [w / total_w for w in widths_px]
    # Clamp extreme column widths to keep numeric columns from becoming too wide
    clamped = []
    for i, frac in enumerate(col_fracs):
        label = str(col_labels[i]).lower()
        min_frac = col_min_frac
        max_frac = col_max_frac
        if label == "name":
            min_frac = 0.18
            max_frac = 0.30
        clamped.append(min(max(frac, min_frac), max_frac))
    # Re-normalize to sum to 1.0
    s = sum(clamped)
    if s > 0:
        col_fracs = [c / s for c in clamped]
    # Compute per-row pixel heights based on text extents, convert to axes fraction
    axes_h_px = None
    if renderer is not None:
        try:
            axes_h_px = ax.get_window_extent(renderer=renderer).height
        except Exception:
            axes_h_px = None
    # Fallback axes height if needed (avoid division by zero)
    if not axes_h_px or axes_h_px <= 0:
        axes_h_px = 800.0  # arbitrary but reasonable default
    row_heights_frac = []
    for r in range(nrows):
        max_h_px = 1.0
        if renderer is not None:
            for c in range(ncols):
                try:
                    bbox = table[(r, c)].get_text().get_window_extent(renderer=renderer)
                    max_h_px = max(max_h_px, bbox.height)
                except Exception:
                    pass
        # padding factor to provide breathing room around text
        padded_px = max_h_px * row_pad
        # Convert to axes fraction with minimum clamp
        rh = max(padded_px / axes_h_px, min_row_frac)
        row_heights_frac.append(rh)

    # Apply column widths and per-row heights
    for r in range(nrows):
        for c in range(ncols):
            cell = table[(r, c)]
            cell.set_width(table_width * col_fracs[c])
            cell.set_height(row_heights_frac[r])
            # Left-align the "Name" column if present
            label = str(col_labels[c]).lower()
            if label == "name":
                cell.get_text().set_ha('left')

def plot_base_stresses(ax, slice_df, scale_frac=0.5, alpha=0.3):
    """
    Plots base normal stresses for each slice as bars.

    Parameters:
        ax: matplotlib Axes object
        slice_df: DataFrame containing slice data
        scale_frac: Fraction of plot height for bar scaling
        alpha: Transparency for bars

    Returns:
        None
    """
    u = slice_df['u'].values  # pore pressure (stress)
    dl = slice_df['dl'].values
    with np.errstate(divide='ignore', invalid='ignore'):
        # convert effective normal force to stress; a degenerate dl→0 yields
        # inf/nan, filtered out below rather than warning.
        n_eff = slice_df['n_eff'].values / dl
    heights = slice_df['y_ct'] - slice_df['y_cb']
    max_ht = heights.max() if not heights.empty else 1.0
    max_bar_len = max_ht * scale_frac

    # Guard against degenerate slices (dl→0 gives inf/nan n_eff) and an all-zero
    # stress field: use only finite values and never divide by zero below.
    finite_neff = n_eff[np.isfinite(n_eff)]
    max_stress = float(np.max(np.abs(finite_neff))) if finite_neff.size else 1.0
    finite_u = u[np.isfinite(u)]
    max_u = float(np.max(np.abs(finite_u))) if finite_u.size else 1.0
    # Scale BOTH the effective-stress and pore-pressure bars by a shared reference
    # so neither overruns the plot. Using max_stress alone blows up the pore bars
    # in rapid drawdown, where pore pressure ≫ the (tiny) effective normal stress.
    ref = max(max_stress, max_u)
    if ref == 0:
        ref = 1.0

    for i, (index, row) in enumerate(slice_df.iterrows()):
        if i >= len(n_eff):
            break

        x1, y1 = row['x_l'], row['y_lb']
        x2, y2 = row['x_r'], row['y_rb']

        stress = n_eff[i]
        pore = u[i]

        if not np.isfinite(stress):
            continue                 # degenerate slice (dl→0); nothing to draw
        if not np.isfinite(pore):
            pore = 0.0

        dx = x2 - x1
        dy = y2 - y1
        length = np.hypot(dx, dy)
        if length == 0:
            continue

        nx = -dy / length
        ny = dx / length

        # --- Normal stress trapezoid ---
        bar_len = (abs(stress) / ref) * max_bar_len
        direction = -np.sign(stress)

        x1_top = x1 + direction * bar_len * nx
        y1_top = y1 + direction * bar_len * ny
        x2_top = x2 + direction * bar_len * nx
        y2_top = y2 + direction * bar_len * ny

        poly_x = [x1, x2, x2_top, x1_top]
        poly_y = [y1, y2, y2_top, y1_top]

        ax.fill(poly_x, poly_y, facecolor='none', edgecolor='red' if stress <= 0 else 'limegreen', hatch='.....',
                linewidth=1, gid='EFF_NORMAL_STRESS')

        # --- Pore pressure trapezoid ---
        u_len = (pore / ref) * max_bar_len
        u_dir = -1  # always into the base

        ux1_top = x1 + u_dir * u_len * nx
        uy1_top = y1 + u_dir * u_len * ny
        ux2_top = x2 + u_dir * u_len * nx
        uy2_top = y2 + u_dir * u_len * ny

        poly_ux = [x1, x2, ux2_top, ux1_top]
        poly_uy = [y1, y2, uy2_top, uy1_top]

        ax.fill(poly_ux, poly_uy, color='blue', alpha=alpha, edgecolor='k', linewidth=1, gid='PORE_PRESSURE')


def plot_thrust_line_from_df(ax, slice_df,
                            color: str = 'red',
                            linestyle: str = '--',
                            linewidth: float = 1,
                            label: str = 'Line of Thrust'):
    """
    Plots the line of thrust from the slice dataframe.

    Parameters:
        ax: matplotlib Axes object
        slice_df: DataFrame containing slice data with 'yt_l' and 'yt_r' columns
        color: Color of the line
        linestyle: Style of the line
        linewidth: Width of the line
        label: Label for the line in the legend

    Returns:
        None
    """
    # Check if required columns exist
    if 'yt_l' not in slice_df.columns or 'yt_r' not in slice_df.columns:
        return
    
    # Create thrust line coordinates from slice data
    thrust_xs = []
    thrust_ys = []
    
    for _, row in slice_df.iterrows():
        # Add left point of current slice
        thrust_xs.append(row['x_l'])
        thrust_ys.append(row['yt_l'])
        
        # Add right point of current slice (same as left point of next slice)
        thrust_xs.append(row['x_r'])
        thrust_ys.append(row['yt_r'])
    
    # Plot the thrust line
    ax.plot(thrust_xs, thrust_ys,
            color=color,
            linestyle=linestyle,
            linewidth=linewidth,
            label=label,
            gid='LINE_OF_THRUST')

def compute_ylim(data, slice_df, scale_frac=0.5, pad_fraction=0.1):
    """
    Computes y-limits for plotting based on slice data.

    Parameters:
        data: Input data
        slice_df: pandas.DataFrame with slice data, must have 'y_lt' and 'y_lb' for stress‐bar sizing
        scale_frac: fraction of max slice height used when drawing stress bars
        pad_fraction: fraction of total range to pad above/below finally

    Returns:
        (y_min, y_max) suitable for ax.set_ylim(...)
    """
    import numpy as np

    y_vals = []

    # 1) collect all profile line elevations
    for line in data.get('profile_lines', []):
        coords = line['coords']
        if hasattr(coords, "xy"):
            _, ys = coords.xy
        else:
            _, ys = zip(*coords)
        y_vals.extend(ys)

    # 2) explicitly include the deepest allowed depth
    if "max_depth" in data and data["max_depth"] is not None:
        y_vals.append(data["max_depth"])

    # 2b) include the full vertical extent of the domain polygon. This is the
    # geometry source for polygon inputs (where profile_lines/max_depth are empty),
    # and is a no-op for profile inputs (the domain spans ground..max_depth).
    domain = data.get('domain_polygon')
    if domain is not None:
        _, dy0, _, dy1 = domain.bounds
        y_vals.extend([dy0, dy1])

    # 2c) line-load arrows extend above the ground at their application points;
    # their tails (and the label riding on them) must stay inside the view
    y_vals.extend(y for _, y in _line_load_tails(data))

    if not y_vals:
        return 0.0, 1.0

    y_min = min(y_vals)
    y_max = max(y_vals)

    # 3) ensure the largest stress bar will fit
    #    stress‐bar length = scale_frac * slice height
    heights = slice_df["y_lt"] - slice_df["y_lb"]
    if not heights.empty:
        max_bar = heights.max() * scale_frac
        y_min -= max_bar
        y_max += max_bar

    # 4) account for distributed loads extending above ground surface
    gamma_w = data.get('gamma_water', 62.4)
    for dloads in [data.get('dloads', []), data.get('dloads2', [])]:
        if dloads:
            for line in dloads:
                for pt in line:
                    # dload arrows extend above surface by load/gamma_w (water depth equivalent)
                    load = pt.get('Normal', 0)
                    if load > 0:
                        y_max = max(y_max, pt.get('Y', 0) + load / gamma_w)

    # 5) add a final small pad
    pad = (y_max - y_min) * pad_fraction
    return y_min - pad, y_max + pad

# ========== FOR PLOTTING INPUT DATA  =========

def plot_reinforcement_lines(ax, slope_data, style=None):
    """
    Plots the reinforcement lines from slope_data.
    
    Parameters:
        ax: matplotlib Axes object
        slope_data: Dictionary containing slope data with 'reinforce_lines' key
        
    Returns:
        None
    """
    if 'reinforce_lines' not in slope_data or not slope_data['reinforce_lines']:
        return

    from .style import resolve_style, feature_style
    rfs = feature_style(resolve_style(style), "reinforcement")
    tension_points_plotted = False  # Track if tension points have been added to legend

    for i, line in enumerate(slope_data['reinforce_lines']):
        # Extract x and y coordinates from the line points
        xs = [point['X'] for point in line]
        ys = [point['Y'] for point in line]

        # Plot the reinforcement line with a distinctive style
        ax.plot(xs, ys, color=rfs.get('color', 'darkgray'),
                linewidth=rfs.get('linewidth', 3), linestyle=rfs.get('linestyle', '-'),
                alpha=rfs.get('alpha', 0.8), label='Reinforcement Line' if i == 0 else "")
        
        # Add markers at each point to show tension values
        for j, point in enumerate(line):
            tension = point.get('T', 0.0)
            if tension > 0:
                # Use smaller marker size proportional to tension (normalized)
                max_tension = max(p.get('T', 0.0) for p in line)
                marker_size = 10 + 15 * (tension / max_tension) if max_tension > 0 else 10
                ax.scatter(point['X'], point['Y'], s=marker_size, 
                          color='red', alpha=0.7, zorder=5,
                          label='Tension Points' if not tension_points_plotted else "")
                tension_points_plotted = True


def _line_load_tails(slope_data):
    """Arrow-tail points for the line-load glyphs: the head sits on the point
    of application, the tail 6% of the model span opposite the force
    direction. Both plot_line_loads (drawing) and compute_ylim (view limits)
    use these, so the arrow and its label are always inside the axes."""
    import numpy as np
    loads = slope_data.get('line_loads') or []
    if not loads:
        return []
    gs = slope_data.get('ground_surface')
    if gs is not None and not gs.is_empty:
        xs = [p[0] for p in gs.coords]
        span = max(xs) - min(xs)
    else:
        span = 100.0
    alen = 0.06 * span
    tails = []
    for ll in loads:
        ang = np.radians(ll.get('angle', -90.0))
        tails.append((ll['x'] - np.cos(ang) * alen,
                      ll['y'] - np.sin(ang) * alen))
    return tails


def plot_line_loads(ax, slope_data, style=None):
    """
    Plots line loads (v12 'lloads') as arrows at their points of application on
    the ground surface. The arrow points IN the direction the force acts (a
    straight-down load draws as a downward arrow ending at the point).

    Parameters:
        ax: matplotlib Axes object
        slope_data: Dictionary containing slope data with 'line_loads' key
        style: optional style sheet (see xslope.style); None -> defaults.
    """
    loads = slope_data.get('line_loads') or []
    if not loads:
        return
    tails = _line_load_tails(slope_data)

    for i, (ll, (tx, ty)) in enumerate(zip(loads, tails)):
        # tail offset opposite the force direction so the head lands on the point
        ax.annotate('', xy=(ll['x'], ll['y']), xytext=(tx, ty),
                    arrowprops=dict(arrowstyle='-|>', color='purple', lw=2),
                    annotation_clip=False)
        ax.annotate(f"L={ll['P']:.0f}", (tx, ty),
                    textcoords="offset points", xytext=(4, 4),
                    fontsize=8, color='purple', fontweight='bold')
        # annotations don't autoscale: an invisible data point at the tail
        # makes the axes grow to keep the arrow (and its label) in view
        ax.plot([tx], [ty], linestyle='None')
        if i == 0:
            ax.plot([], [], color='purple', lw=2, label='Line Load')


def plot_piles(ax, slope_data, slice_df=None, style=None):
    """
    Plots pile lines from slope_data and optionally marks failure surface intersections.

    Parameters:
        ax: matplotlib Axes object
        slope_data: Dictionary containing slope data with 'pile_lines' key
        slice_df: Optional DataFrame — if provided, marks pile-failure surface intersection points
        style: optional style sheet (see xslope.style); None → defaults. Piles are
            structural, so color + width only (always solid).
    """
    if 'pile_lines' not in slope_data or not slope_data['pile_lines']:
        return

    from .style import resolve_style, feature_style
    pf = feature_style(resolve_style(style), "piles")
    pcolor = pf.get('color', 'green')
    plw = pf.get('linewidth', 4)

    for i, pile in enumerate(slope_data['pile_lines']):
        xs = [pile['x1'], pile['x2']]
        ys = [pile['y1'], pile['y2']]
        ax.plot(xs, ys, color=pcolor, linewidth=plw, linestyle='-',
                alpha=0.9, solid_capstyle='butt',
                label='Pile' if i == 0 else "")
        # Annotate with H value
        if pile.get('H') is not None:
            mid_x = (pile['x1'] + pile['x2']) / 2
            mid_y = (pile['y1'] + pile['y2']) / 2
            ax.annotate(f"H={pile['H']:.0f}", (mid_x, mid_y),
                        textcoords="offset points", xytext=(8, 0),
                        fontsize=8, color=pcolor, fontweight='bold')

    # Mark failure surface intersection points from slice_df
    if slice_df is not None and 'h_pile' in slice_df.columns:
        pile_slices = slice_df[slice_df['h_pile'] > 0]
        if not pile_slices.empty:
            ax.scatter(pile_slices['x_pile'], pile_slices['y_pile'],
                       marker='o', s=40, color='red', zorder=6,
                       label='Pile-Surface Intersection')


# Figure-fraction y the legend's bottom is pinned to (consistent across plots).
_LEGEND_BOTTOM = 0.03
# Compact, consistent text sizes across all plots.
_TITLE_FONTSIZE = 11
_TICK_FONTSIZE = 8
_LEGEND_FONTSIZE = 8
# Tight legend packing (matplotlib's defaults — columnspacing 2.0, handlelength 2.0,
# handletextpad 0.8 — are wide and waste row width, forcing long-label legends onto
# extra rows). Applied to BOTH the column-fit trial and the final legend so the
# measured width matches what's drawn.
_LEGEND_PACK = dict(columnspacing=1.2, handlelength=1.6, handletextpad=0.5)


def _fit_legend_ncol(ax, fig, handles, labels, anchor):
    """Choose a legend column count: the fewest rows whose balanced columns fit.

    The legend is anchored to (and centered on) the whole FIGURE width — not the
    axes box, which equal-aspect box-adjust can shrink — so candidates are measured
    against the figure width with the same tight packing the final legend uses
    (_LEGEND_PACK). We try 1 row, then 2, … and take the first whose balanced
    column count (ceil(n/rows)) fits; e.g. a 12-item legend lays out 6×2. Falls
    back to two rows when no renderer is available (a backend without one)."""
    n = len(labels)
    if n <= 1:
        return max(1, n)
    try:
        fig.canvas.draw()                      # ensure a renderer + axes layout
        renderer = fig.canvas.get_renderer()
        budget = fig.bbox.width                # frameless, centered: may fill the width
        for rows in range(1, n + 1):
            ncol = math.ceil(n / rows)
            leg = ax.legend(handles=handles, labels=labels, loc="upper center",
                            bbox_to_anchor=anchor, ncol=ncol,
                            fontsize=_LEGEND_FONTSIZE, **_LEGEND_PACK)
            fits = leg.get_window_extent(renderer).width <= budget
            leg.remove()
            if fits:
                return ncol                    # fewest rows that fit width-wise
        return 1
    except Exception:
        # No renderer: two rows is the common neat layout for these legends.
        return max(1, math.ceil(n / 2))


def _legend_below(ax, fig, anchor=(0.5, -0.12), handles=None, labels=None,
                  legend_ncol="auto", show_legend=True, **kw):
    """Draw a horizontal legend below the axes. With legend_ncol="auto" the column
    count is auto-fit — as wide as fits the axes width, with the fewest, neatly-
    balanced rows (see _fit_legend_ncol); pass an int to force a column count.
    Pass explicit handles/labels, else they come from the axes. Extra kwargs
    (e.g. frameon) pass through to ax.legend. show_legend=False skips the legend
    itself but still styles the title/ticks and reclaims its reserved bottom room."""
    # Title + tick sizes are set here (the one call every geometry plot ends with)
    # whether or not the legend is drawn.
    ax.tick_params(labelsize=_TICK_FONTSIZE)
    if ax.get_title():
        ax.title.set_fontsize(_TITLE_FONTSIZE)
    if not show_legend:
        try:
            fig.set_layout_engine("none")
        except Exception:
            try:
                fig.set_tight_layout(False)
            except Exception:
                pass
        top = 0.94
        if ax.get_title():                   # title-aware top; modest bottom for ticks
            try:
                fig.canvas.draw()
                th = ax.title.get_window_extent(fig.canvas.get_renderer()).height / fig.bbox.height
                top = max(0.80, min(0.94, 1.0 - th - 0.03))
            except Exception:
                pass
        fig.subplots_adjust(top=top, bottom=0.11)
        return None
    if handles is None:
        handles, labels = ax.get_legend_handles_labels()
    if labels is None:                       # handles given without explicit labels
        labels = [h.get_label() for h in handles]
    # De-duplicate by label, preserving order, and drop matplotlib's "_"-hidden
    # labels — so e.g. N failure circles all labeled "Circle" yield one entry.
    if handles:
        seen, dh, dl = set(), [], []
        for h, l in zip(handles, labels):
            if l in seen or (isinstance(l, str) and l.startswith("_")):
                continue
            seen.add(l); dh.append(h); dl.append(l)
        handles, labels = dh, dl
    if not handles:
        return None
    kw.setdefault("frameon", False)          # frameless by default; toggle via legend_frame
    kw.setdefault("fontsize", _LEGEND_FONTSIZE)
    for k, v in _LEGEND_PACK.items():        # tight packing; must match _fit_legend_ncol
        kw.setdefault(k, v)
    if legend_ncol == "auto":
        ncol = _fit_legend_ncol(ax, fig, handles, labels, anchor)
    else:
        ncol = max(1, int(legend_ncol))
    # Pin the legend to a FIXED position at the figure bottom (figure coordinates),
    # not relative to the axes — otherwise the gap below the axes scales with the
    # axes height, so plots with a taller data range (e.g. inputs showing circle
    # centers high up) push the legend far down while compact plots hug it. Figure-
    # anchoring makes the placement identical across every plot type.
    gap = 0.035                              # fixed figure-fraction gap below the axis
    # Anchor the legend by its TOP (loc="upper center"); its exact y is set below,
    # once layout tells us where the axes box actually sits.
    leg = ax.legend(handles=handles, labels=labels, loc="upper center",
                    bbox_to_anchor=(0.5, _LEGEND_BOTTOM), bbox_transform=fig.transFigure,
                    ncol=ncol, **kw)
    # Manual margins must survive the next canvas draw — an active tight/constrained
    # layout engine would repack the title against the top edge (clipping it) and
    # undo the reserved bottom. Disable it, then set top (title room) and bottom
    # (legend room) explicitly.
    try:
        fig.set_layout_engine("none")
    except Exception:
        try:
            fig.set_tight_layout(False)
        except Exception:
            pass
    try:
        fig.canvas.draw()
        leg_h = leg.get_window_extent().height / fig.bbox.height
    except Exception:
        n_rows = max(1, math.ceil(len(labels) / max(1, ncol)))
        leg_h = 0.045 * n_rows
    # Reserve a bottom margin sized for the legend (as if the axes filled down to it)
    # and a title-aware top, then lay out.
    bottom = min(0.55, _LEGEND_BOTTOM + leg_h + gap)
    top = 0.94
    if ax.get_title():                       # reserve enough top for the (maybe multi-line) title
        try:
            th = ax.title.get_window_extent(fig.canvas.get_renderer()).height / fig.bbox.height
            top = max(0.80, min(0.94, 1.0 - th - 0.03))
        except Exception:
            pass
    fig.subplots_adjust(top=top, bottom=bottom)
    # Pin the legend just under the REAL axes-box bottom so it hugs the x-axis on
    # every plot. For datalim plots the box fills the figure, so its bottom sits at
    # the reserved margin (unchanged placement). For box-adjust plots the box can be
    # short and float high (a wide-thin seep domain) — the legend follows it up
    # instead of stranding at the figure floor. (The title rides the box top the
    # same way, so title + plot + legend stay grouped.)
    try:
        fig.canvas.draw()
        ax.apply_aspect()                    # finalize the box-adjust position now
        # Use the DRAWN box extent (reflects box-adjust shrink + any colorbar
        # divider), converted to a figure fraction — get_position() can lag the
        # aspect on the first draw.
        y0 = ax.get_window_extent().y0 / fig.bbox.height
    except Exception:
        y0 = bottom
    anchor_y = max(_LEGEND_BOTTOM + leg_h, y0 - gap)
    leg.set_bbox_to_anchor((0.5, anchor_y), transform=fig.transFigure)
    return leg


def material_legend_handles(materials, style=None, alpha=None):
    """Filled-patch legend handles for material zones — the single consistent
    material-zone swatch used across the inputs / mesh / seep / fem plots, so a zone
    looks the same (filled patch in its style color) in every view. One patch per
    material, labeled by material name. ``alpha`` overrides the per-material style
    alpha when given (e.g. a busier result view that wants lighter zones)."""
    import matplotlib.patches as mpatches
    from .style import resolve_style, material_style
    style = resolve_style(style)
    handles = []
    for i, m in enumerate(materials or []):
        ms = material_style(style, i)
        a = ms.get("alpha", 0.6) if alpha is None else alpha
        name = (m.get("name") if isinstance(m, dict) else m) or f"Material {i + 1}"
        handles.append(mpatches.Patch(facecolor=ms["color"], alpha=a,
                                      edgecolor="none", label=name))
    return handles


def plot_inputs(
    slope_data,
    title="Slope Geometry and Inputs",
    figsize=(12, 7),
    mat_table=False,
    save_png=False,
    save_dxf=False,
    dpi=300,
    mode="lem",
    tab_loc="top",
    legend_ncol="auto",
    legend_frame=False,
    show_title=True,
    show_legend=True,
    label_coordinates=False,
    coord_label_size=7,
    coord_arrows=False,
    fig=None,
    style=None,
):
    """
    Creates a plot showing the slope geometry and input parameters.

    Parameters:
        slope_data: Dictionary containing plot data
        title: Title for the plot
        figsize: Tuple of (width, height) in inches for the plot
        mat_table: Controls material table display. Can be:
            - True: Use tab_loc for positioning (default)
            - False: Don't show material table
            - 'auto': Use tab_loc for positioning
            - String: Specific location from valid placements (see tab_loc)
        save_png: If True, save plot as PNG file (default: False)
        dpi: Resolution for saved PNG file (default: 300)
        mode: Which material properties table to display:
            - "lem": Limit equilibrium materials (γ, c, φ, optional d/ψ)
            - "seep": Seepage properties (k₁, k₂, Angle, kr₀, h₀)
            - "fem": FEM properties (γ, c, φ, E, ν)
        tab_loc: Table placement when mat_table is True or 'auto'. Valid options:
            - "upper left": Top-left corner of plot area
            - "upper right": Top-right corner of plot area
            - "upper center": Top-center of plot area
            - "lower left": Bottom-left corner of plot area
            - "lower right": Bottom-right corner of plot area
            - "lower center": Bottom-center of plot area
            - "center left": Middle-left of plot area
            - "center right": Middle-right of plot area
            - "center": Center of plot area
            - "top": Above plot area, horizontally centered
        legend_ncol: Legend column count. Use "auto" (default) to lay it out
            automatically — as wide as fits the axes, with the fewest, neatly
            balanced rows — or pass an int to force a specific column count.
        label_coordinates: If True, annotate every unique profile-line /
            polygon vertex with its "(x, y)" coordinates, in the style of the
            vendor verification manuals (default: False).
        coord_label_size: Font size in points for the coordinate labels
            (default: 7).
        coord_arrows: If True, labels in dense clusters are pushed clear and
            tied back to their vertices with thin leader lines; by default
            (False) labels stay adjacent to their vertices with no leaders.
        fig: Optional existing Matplotlib Figure to draw into (used for embedding in a
            GUI canvas). When None (default) a new pyplot figure is created and shown;
            when provided, the figure is cleared and reused and plt.show() is skipped.

    Returns:
        The Matplotlib Figure that was drawn into.
    """
    own_fig = fig is None
    if own_fig:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig.clear()
        ax = fig.add_subplot(111)

    from .style import resolve_style
    style = resolve_style(style)

    # Plot mesh in background if available
    mesh = slope_data.get('mesh')
    if mesh is not None:
        from matplotlib.collections import LineCollection
        m_nodes = mesh["nodes"]
        m_elements = mesh["elements"]
        m_etypes = mesh["element_types"]
        lines = []
        for elem, etype in zip(m_elements, m_etypes):
            if etype in (3, 6):  # tri3 / tri6 – corner edges only
                edges = [(elem[0], elem[1]), (elem[1], elem[2]), (elem[2], elem[0])]
            elif etype in (4, 8, 9):  # quad4 / quad8 / quad9
                edges = [(elem[0], elem[1]), (elem[1], elem[2]),
                         (elem[2], elem[3]), (elem[3], elem[0])]
            else:
                continue
            for n0, n1 in edges:
                lines.append(m_nodes[[n0, n1]])
        if lines:
            from .style import feature_style as _fs
            mfs = _fs(style, "mesh")
            lc = LineCollection(lines, colors=mfs.get("color", "gray"),
                                alpha=mfs.get("alpha", 0.25),
                                linewidths=mfs.get("linewidth", 0.5),
                                linestyles=mfs.get("linestyle", "-"))
            ax.add_collection(lc)

    # Plot geometry: profile lines if provided (drawn as before), otherwise the
    # material-zone polygons.
    plot_base_geometry(ax, slope_data, labels=True, style=style)
    if mode == "fem" or (mode == "lem" and any(m.get('u') == 'piezo' for m in slope_data.get('materials', []))):
        plot_piezo_line(ax, slope_data, style=style)
    if mode == "seep":
        plot_seepage_bc_lines(ax, slope_data, style=style)
    if mode != "seep":
        plot_dloads(ax, slope_data, style=style)
    plot_tcrack_surface(ax, slope_data, style=style)
    plot_reinforcement_lines(ax, slope_data, style=style)
    plot_piles(ax, slope_data, style=style)
    plot_line_loads(ax, slope_data, style=style)

    if mode == "lem":
        if slope_data['circular']:
            plot_circles(ax, slope_data, style=style)
        elif slope_data.get('non_circ') and len(slope_data['non_circ']) > 0:
            plot_non_circ(ax, slope_data['non_circ'], style=style)

    # Seismic coefficient annotation
    k_seismic = slope_data.get('k_seismic', 0.0)
    if k_seismic and mode in ("lem", "fem"):
        if mode == "fem":
            # FEM: sign of k directly gives direction (+x → right, -x → left)
            arrow = "\u2192" if k_seismic > 0 else "\u2190"
        else:
            # LEM: k acts toward the toe (downslope), infer from ground surface
            gs = slope_data.get('ground_surface')
            if gs is not None and not gs.is_empty:
                coords = list(gs.coords)
                y_left = coords[0][1]
                y_right = coords[-1][1]
                y_peak = max(c[1] for c in coords)
                # Dam detection: both ends are substantially lower than the peak
                threshold = 0.3 * (y_peak - min(y_left, y_right))
                if (y_peak - y_left) > threshold and (y_peak - y_right) > threshold:
                    arrow = "\u2194"  # ↔  both faces
                elif y_left > y_right:
                    arrow = "\u2192"  # →  toe on right
                else:
                    arrow = "\u2190"  # ←  toe on left
            else:
                arrow = ""
        k_text = f"k = {abs(k_seismic):g} {arrow}".strip()
        ax.text(0.98, 0.97, k_text, transform=ax.transAxes,
                fontsize=10, fontweight='bold', ha='right', va='top',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='lightyellow',
                          edgecolor='orange', linewidth=1.0, alpha=0.9))

    # Handle material table display
    if mat_table:
        # Helpers to adapt slope_data materials into formats expected by table functions
        def _build_seep_data():
            materials = slope_data.get('materials', [])
            return {
                "k1_by_mat": [m.get('k1', 0.0) for m in materials],
                "k2_by_mat": [m.get('k2', 0.0) for m in materials],
                "angle_by_mat": [m.get('alpha', 0.0) for m in materials],
                "kr0_by_mat": [m.get('kr0', 0.0) for m in materials],
                "h0_by_mat": [m.get('h0', 0.0) for m in materials],
                "material_names": [m.get('name', f"Material {i+1}") for i, m in enumerate(materials)],
            }

        def _build_fem_data():
            materials = slope_data.get('materials', [])
            return {
                "c_by_mat": [m.get('c', 0.0) for m in materials],
                "phi_by_mat": [m.get('phi', 0.0) for m in materials],
                "E_by_mat": [m.get('E', 0.0) for m in materials],
                "nu_by_mat": [m.get('nu', 0.0) for m in materials],
                "gamma_by_mat": [m.get('gamma', 0.0) for m in materials],
                "material_names": [m.get('name', f"Material {i+1}") for i, m in enumerate(materials)],
            }

        def _estimate_table_dims():
            """Estimate table dimensions based on mode and materials."""
            materials = slope_data.get('materials', [])
            num_rows = max(1, len(materials))

            if mode == "lem":
                has_d_psi = any(mat.get('d', 0) > 0 or mat.get('psi', 0) > 0 for mat in materials)
                width = 0.25 if has_d_psi else 0.2
                height = min(0.35, 0.06 + 0.035 * num_rows)
            elif mode == "fem":
                width = 0.60
                height = min(0.32, 0.06 + 0.035 * num_rows)
            elif mode == "seep":
                width = 0.45
                height = min(0.50, 0.10 + 0.06 * num_rows)
            else:
                raise ValueError(f"Unknown mode '{mode}'. Expected one of: 'lem', 'seep', 'fem'.")

            return width, height

        def _plot_table(ax, xloc, yloc):
            """Plot the appropriate material table based on mode."""
            if mode == "lem":
                plot_lem_material_table(ax, slope_data['materials'], xloc=xloc, yloc=yloc)
            elif mode == "seep":
                plot_seep_material_table(ax, _build_seep_data(), xloc=xloc, yloc=yloc)
            elif mode == "fem":
                width, height = _estimate_table_dims()
                plot_fem_material_table(ax, _build_fem_data(), xloc=xloc, yloc=yloc, width=width, height=height)

        def _calculate_position(location, margin=0.03):
            """Calculate xloc, yloc for a given location string."""
            width, height = _estimate_table_dims()

            # Location map with default coordinates
            position_map = {
                'upper left': (margin, max(0.0, 1.0 - margin - height)),
                'upper right': (max(0.0, 1.0 - width - margin), max(0.0, 1.0 - margin - height)),
                'upper center': (0.35, 0.70),
                'lower left': (0.05, 0.05),
                'lower right': (0.70, 0.05),
                'lower center': (0.35, 0.05),
                'center left': (0.05, 0.35),
                'center right': (0.70, 0.35),
                'center': (0.35, 0.35),
                'top': ((1.0 - width) / 2.0, 1.16)
            }

            return position_map.get(location, position_map['upper right'])

        # Determine which location to use
        placement = mat_table if isinstance(mat_table, str) and mat_table != 'auto' else tab_loc

        # Validate placement
        valid_placements = ['upper left', 'upper right', 'upper center', 'lower left',
                          'lower right', 'lower center', 'center left', 'center right', 'center', 'top']
        if placement not in valid_placements:
            raise ValueError(f"Unknown placement '{placement}'. Expected one of: {', '.join(valid_placements)}.")

        # Calculate position and plot table
        xloc, yloc = _calculate_position(placement)
        _plot_table(ax, xloc, yloc)

        # Adjust y-limits to prevent table overlap with plot data
        if placement in ("upper left", "upper right", "upper center"):
            _, height = _estimate_table_dims()
            margin = 0.03
            bottom_fraction = max(0.0, 1.0 - margin - height)

            y_min_curr, y_max_curr = ax.get_ylim()
            y_range = y_max_curr - y_min_curr
            if y_range > 0:
                elem_bounds = get_plot_elements_bounds(ax, slope_data)
                if elem_bounds:
                    y_top = max(b[3] for b in elem_bounds)
                    y_norm = (y_top - y_min_curr) / y_range
                    if y_norm >= bottom_fraction and bottom_fraction > 0:
                        y_max_new = y_min_curr + (y_top - y_min_curr) / bottom_fraction
                        ax.set_ylim(y_min_curr, y_max_new)

    ax.set_aspect('equal', adjustable='datalim')  # ✅ Equal aspect

    # Add a bit of headroom so plotted lines/markers don't touch the top border
    y0, y1 = ax.get_ylim()
    if y1 > y0:
        pad = 0.05 * (y1 - y0)
        ax.set_ylim(y0, y1 + pad)
    ax.grid(False)

    # Coordinate labels go on AFTER the aspect and limits are final: their
    # collision layout measures label boxes in display pixels, so placing them
    # earlier (before the equal-aspect rescale) invalidates the measurements.
    if label_coordinates:
        plot_coordinate_labels(ax, slope_data, fontsize=coord_label_size,
                               arrows=coord_arrows, style=style)

    # Get legend handles and labels
    handles, labels = ax.get_legend_handles_labels()

    # Add distributed load to legend if present
    if slope_data['dloads']:
        handler_class, dummy_line = get_dload_legend_handler()
        handles.append(dummy_line)
        labels.append('Distributed Load')

    if show_title:
        ax.set_title(title)
    fig.tight_layout()
    # Single shared legend recipe: below the axes, auto-columns, frameless by
    # default (frame toggled per-view via legend_frame). Reserves bottom margin.
    _legend_below(ax, fig, handles=handles, labels=labels,
                  legend_ncol=legend_ncol, frameon=legend_frame, show_legend=show_legend)

    base_name = 'plot_' + title.lower().replace(' ', '_').replace(':', '').replace(',', '')
    if save_png:
        fig.savefig(base_name + '.png', dpi=dpi, bbox_inches='tight')
    if save_dxf:
        from .cad import axes_to_dxf
        axes_to_dxf(ax, base_name + '.dxf')

    if own_fig:
        plt.show()
    return fig

# ========== Main Plotting Function =========

def plot_solution(slope_data, slice_df, failure_surface, results, figsize=(12, 7), slice_numbers=False, seep_contours=True, save_png=False, save_dxf=False, dpi=300, legend_ncol="auto", legend_frame=False, show_title=True, show_legend=True, fig=None, style=None):
    """
    Plots the full solution including slices, numbers, thrust line, and base stresses.

    Parameters:
        data: Input data
        slice_df: DataFrame containing slice data
        failure_surface: Failure surface geometry
        results: Solution results
        figsize: Tuple of (width, height) in inches for the plot
        fig: Optional existing Matplotlib Figure to draw into (used for embedding in a
            GUI canvas). When None (default) a new pyplot figure is created and shown;
            when provided, the figure is cleared and reused and plt.show() is skipped.

    Returns:
        The Matplotlib Figure that was drawn into.
    """
    own_fig = fig is None
    if own_fig:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig.clear()
        ax = fig.add_subplot(111)
    ax.grid(False)

    plot_base_geometry(ax, slope_data, style=style)
    plot_slices(ax, slice_df, fill=False)
    plot_failure_surface(ax, failure_surface)
    if any(m.get('u') == 'piezo' for m in slope_data.get('materials', [])):
        plot_piezo_line(ax, slope_data, style=style)

    # Seep overlays: head contours and phreatic surface when any material uses seep
    has_seep = any(m.get('u') == 'seep' for m in slope_data.get('materials', []))
    mesh = slope_data.get('mesh')
    seep_u = slope_data.get('seep_u')
    if seep_contours and has_seep and mesh is not None and seep_u is not None:
        import matplotlib.tri as mtri
        m_nodes = mesh['nodes']
        m_elements = mesh['elements']
        m_etypes = mesh.get('element_types', np.full(len(m_elements), 3))
        gamma_w = slope_data.get('gamma_water', 62.4)
        head = seep_u / gamma_w + m_nodes[:, 1]

        # Build triangulation for contouring (subdivide higher-order elements)
        all_tris = []
        for idx, elem in enumerate(m_elements):
            etype = m_etypes[idx]
            if etype == 3:
                all_tris.append(elem[:3])
            elif etype == 6:
                all_tris.append([elem[0], elem[3], elem[5]])
                all_tris.append([elem[3], elem[1], elem[4]])
                all_tris.append([elem[5], elem[4], elem[2]])
                all_tris.append([elem[3], elem[4], elem[5]])
            elif etype in (4, 8, 9):
                all_tris.append([elem[0], elem[1], elem[2]])
                all_tris.append([elem[0], elem[2], elem[3]])
        if all_tris:
            triang = mtri.Triangulation(m_nodes[:, 0], m_nodes[:, 1], all_tris)
            # Head contours
            levels = np.linspace(np.min(head), np.max(head), 20)
            ax.tricontour(triang, head, levels=levels, colors='k', linewidths=0.5, alpha=0.5)
            # Phreatic surface (u = 0)
            if np.min(seep_u) < 0:
                cs_phreatic = ax.tricontour(triang, seep_u, levels=[0], colors='black', linewidths=2.0, alpha=0.5)
                # Place inverted triangle marker at the midpoint of the phreatic contour
                from matplotlib.markers import MarkerStyle
                from matplotlib.transforms import offset_copy
                for seg in cs_phreatic.allsegs[0]:
                    if len(seg) > 1:
                        # Compute cumulative arc length along the segment
                        diffs = np.diff(seg, axis=0)
                        arc = np.concatenate([[0], np.cumsum(np.hypot(diffs[:, 0], diffs[:, 1]))])
                        mid_arc = arc[-1] / 2.0
                        mid_x = float(np.interp(mid_arc, arc, seg[:, 0]))
                        mid_y = float(np.interp(mid_arc, arc, seg[:, 1]))
                        # Offset marker so tip touches the line (same technique as piezo line)
                        ms = MarkerStyle("v")
                        path = ms.get_path().transformed(ms.get_transform())
                        tip_offset = (-float(np.asarray(path.vertices)[:, 1].min())) * 8.0 + 2.0
                        trans = offset_copy(ax.transData, fig=ax.figure, x=0.0, y=tip_offset, units="points")
                        ax.plot([mid_x], [mid_y], marker="v", color="black", markersize=8,
                                linestyle="None", transform=trans, alpha=0.5)
                        break  # marker on the longest/first segment only

    plot_dloads(ax, slope_data, style=style)
    plot_tcrack_surface(ax, slope_data, style=style)
    plot_tcrack_water_force(ax, slice_df, slope_data)
    plot_reinforcement_lines(ax, slope_data, style=style)
    plot_piles(ax, slope_data, slice_df=slice_df, style=style)
    plot_line_loads(ax, slope_data, style=style)
    if slice_numbers:
        plot_slice_numbers(ax, slice_df)
    # plot_material_table(ax, data['materials'], xloc=0.75) # Adjust this so that it fits with the legend

    alpha = 0.3
    if results['method'] in ('spencer', 'mprice'):
        # M-P draws the thrust line too once it is computed; plot_thrust_line_from_df
        # skips gracefully while 'yt_l'/'yt_r' are absent (the M-P thrust line is a
        # deferred post-process diagnostic).
        plot_thrust_line_from_df(ax, slice_df)

    plot_base_stresses(ax, slice_df, alpha=alpha)

    import matplotlib.patches as mpatches
    normal_patch = mpatches.Patch(facecolor='none', edgecolor='green', hatch='.....', label="Eff Normal Stress (σ')")
    pore_patch = mpatches.Patch(color='blue', alpha=alpha, label='Pore Pressure (u)')

    # Get legend handles and labels
    handles, labels = ax.get_legend_handles_labels()
    handles.extend([normal_patch, pore_patch])
    labels.extend(["Eff Normal Stress (σ')", 'Pore Pressure (u)'])
    
    # Add distributed load to legend if present
    if slope_data['dloads']:
        handler_class, dummy_line = get_dload_legend_handler()
        handles.append(dummy_line)
        labels.append('Distributed Load')
    
    ax.set_aspect('equal', adjustable='datalim')

    fs = results['FS']
    method = results['method']
    if method == 'oms':
        title = f'OMS: FS = {fs:.3f}'
    elif method == 'bishop':
        title = f'Bishop: FS = {fs:.3f}'
    elif method == 'spencer':
        theta = results['theta']
        title = f'Spencer: FS = {fs:.3f}, θ = {theta:.2f}°'
    elif method == 'janbu':
        fo = results['fo']
        title = f'Janbu-Corrected: FS = {fs:.3f}, fo = {fo:.2f}'
    elif method == 'corps':
        theta = results['theta']
        title = f'Corps Engineers: FS = {fs:.3f}, θ = {theta:.2f}°'
    elif method == 'lowe':
        title = f'Lowe & Karafiath: FS = {fs:.3f}'
    elif method == 'mprice':
        title = (f"Morgenstern-Price ({results.get('f_type', '')}): "
                 f"FS = {fs:.3f}, λ = {results['lambda']:.3f}")
    else:
        title = f'{method}: FS = {fs:.3f}'
    if show_title:
        ax.set_title(title)

    # zoom y‐axis to just cover the slope and depth, with a little breathing room (thrust line can be outside)
    ymin, ymax = compute_ylim(slope_data, slice_df, pad_fraction=0.05)
    ax.set_ylim(ymin, ymax)

    fig.tight_layout()
    _legend_below(ax, fig, handles=handles, labels=labels,
                  legend_ncol=legend_ncol, frameon=legend_frame, show_legend=show_legend)

    base_name = 'plot_' + title.lower().replace(' ', '_').replace(':', '').replace(',', '').replace('°', 'deg')
    if save_png:
        fig.savefig(base_name + '.png', dpi=dpi, bbox_inches='tight')
    if save_dxf:
        from .cad import axes_to_dxf
        axes_to_dxf(ax, base_name + '.dxf')

    if own_fig:
        plt.show()
    return fig

# ========== Functions for Search Results =========

def plot_failure_surfaces(ax, fs_cache):
    """
    Plots all failure surfaces from the factor of safety cache.

    Parameters:
        ax: matplotlib Axes object
        fs_cache: List of dictionaries containing failure surface data and FS values

    Returns:
        None
    """
    for i, result in reversed(list(enumerate(fs_cache))):
        surface = result['failure_surface']
        if surface is None or surface.is_empty:
            continue
        x, y = zip(*surface.coords)
        color = 'red' if i == 0 else 'gray'
        lw = 2 if i == 0 else 1
        ax.plot(x, y, color=color, linestyle='-', linewidth=lw, alpha=1.0 if i == 0 else 0.6,
                gid='CRITICAL_SURFACE' if i == 0 else 'TESTED_SURFACES')

def plot_circle_centers(ax, fs_cache):
    """
    Plots the centers of circular failure surfaces.

    Parameters:
        ax: matplotlib Axes object
        fs_cache: List of dictionaries containing circle center data

    Returns:
        None
    """
    for result in fs_cache:
        ax.plot(result['Xo'], result['Yo'], 'ko', markersize=3, alpha=0.6, gid='CIRCLE_CENTERS')

def plot_search_path(ax, search_path):
    """
    Plots the search path used to find the critical failure surface.

    Parameters:
        ax: matplotlib Axes object
        search_path: List of dictionaries containing search path coordinates

    Returns:
        None
    """
    if len(search_path) < 2:
        return  # need at least two points to draw an arrow

    for i in range(len(search_path) - 1):
        start = search_path[i]
        end = search_path[i + 1]
        dx = end['x'] - start['x']
        dy = end['y'] - start['y']
        ax.arrow(start['x'], start['y'], dx, dy,
                 head_width=1, head_length=2, fc='green', ec='green', length_includes_head=True,
                 gid='SEARCH_PATH')

def plot_circular_search_results(slope_data, fs_cache, search_path=None, circle_cache=None, highlight_fs=True, figsize=(12, 7), save_png=False, save_dxf=False, dpi=300, legend_ncol="auto", legend_frame=False, show_title=True, show_legend=True, fig=None, style=None):
    """
    Creates a plot showing the results of a circular failure surface search.

    Parameters:
        slope_data: Dictionary containing plot data
        fs_cache: List of dictionaries containing failure surface data and FS values
        search_path: List of dictionaries containing search path coordinates
        circle_cache: List of dictionaries containing all tested circles (for plotting)
        highlight_fs: Boolean indicating whether to highlight the critical failure surface
        figsize: Tuple of (width, height) in inches for the plot
        fig: Optional existing Matplotlib Figure to draw into (used for embedding in a
            GUI canvas). When None (default) a new pyplot figure is created and shown;
            when provided, the figure is cleared and reused and plt.show() is skipped.

    Returns:
        The Matplotlib Figure that was drawn into.
    """
    own_fig = fig is None
    if own_fig:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig.clear()
        ax = fig.add_subplot(111)

    plot_base_geometry(ax, slope_data, style=style)
    if any(m.get('u') == 'piezo' for m in slope_data.get('materials', [])):
        plot_piezo_line(ax, slope_data, style=style)
    plot_dloads(ax, slope_data, style=style)
    plot_tcrack_surface(ax, slope_data, style=style)

    # Plot all tested circles from circle_cache (light gray)
    if circle_cache:
        first_plotted = True
        for result in circle_cache:
            surface = result.get('failure_surface')
            if surface is None or surface.is_empty:
                continue
            x, y = zip(*surface.coords)
            label = 'Tested Circle' if first_plotted else None
            ax.plot(x, y, color='gray', linestyle='-', linewidth=0.5, alpha=0.5, label=label, gid='TESTED_SURFACES')
            first_plotted = False

    # Plot only the critical circle from fs_cache (red)
    if fs_cache:
        critical = fs_cache[0]
        surface = critical.get('failure_surface')
        if surface is not None and not surface.is_empty:
            x, y = zip(*surface.coords)
            ax.plot(x, y, color='red', linestyle='-', linewidth=2, label='Critical Circle', gid='CRITICAL_SURFACE')
        # Plot critical circle center
        ax.plot(critical['Xo'], critical['Yo'], 'ro', markersize=5, gid='CIRCLE_CENTERS')

    # Plot all circle centers from fs_cache
    plot_circle_centers(ax, fs_cache)

    if search_path:
        plot_search_path(ax, search_path)

    ax.set_aspect('equal', adjustable='datalim')
    ax.grid(False)
    if show_title and highlight_fs and fs_cache:
        critical_fs = fs_cache[0]['FS']
        ax.set_title(f"Critical Factor of Safety = {critical_fs:.3f}")

    fig.tight_layout()
    _legend_below(ax, fig, legend_ncol=legend_ncol, frameon=legend_frame, show_legend=show_legend)

    if save_png:
        fig.savefig('plot_circular_search_results.png', dpi=dpi, bbox_inches='tight')
    if save_dxf:
        from .cad import axes_to_dxf
        axes_to_dxf(ax, 'plot_circular_search_results.dxf')

    if own_fig:
        plt.show()
    return fig

def plot_noncircular_search_results(slope_data, fs_cache, search_path=None, highlight_fs=True, figsize=(12, 7), save_png=False, save_dxf=False, dpi=300, legend_ncol="auto", legend_frame=False, show_title=True, show_legend=True, fig=None, style=None):
    """
    Creates a plot showing the results of a non-circular failure surface search.

    Parameters:
        slope_data: Dictionary containing plot data
        fs_cache: List of dictionaries containing failure surface data and FS values
        search_path: List of dictionaries containing search path coordinates
        highlight_fs: Boolean indicating whether to highlight the critical failure surface
        figsize: Tuple of (width, height) in inches for the plot
        fig: Optional existing Matplotlib Figure to draw into (used for embedding in a
            GUI canvas). When None (default) a new pyplot figure is created and shown;
            when provided, the figure is cleared and reused and plt.show() is skipped.

    Returns:
        The Matplotlib Figure that was drawn into.
    """
    own_fig = fig is None
    if own_fig:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig.clear()
        ax = fig.add_subplot(111)

    # Plot basic profile elements
    plot_base_geometry(ax, slope_data, style=style)
    if any(m.get('u') == 'piezo' for m in slope_data.get('materials', [])):
        plot_piezo_line(ax, slope_data, style=style)
    plot_dloads(ax, slope_data, style=style)
    plot_tcrack_surface(ax, slope_data, style=style)

    # Plot all failure surfaces from cache
    first_tested = True
    for i, result in reversed(list(enumerate(fs_cache))):
        surface = result['failure_surface']
        if surface is None or surface.is_empty:
            continue
        x, y = zip(*surface.coords)
        if i == 0:
            # Critical surface
            ax.plot(x, y, color='red', linestyle='-', linewidth=2, alpha=1.0, label='Critical Surface', gid='CRITICAL_SURFACE')
        else:
            # Tested surfaces
            label = 'Tested Surface' if first_tested else None
            ax.plot(x, y, color='gray', linestyle='-', linewidth=1, alpha=0.6, label=label, gid='TESTED_SURFACES')
            first_tested = False

    # Plot search path if provided
    if search_path:
        for i in range(len(search_path) - 1):
            start = search_path[i]
            end = search_path[i + 1]
            # For non-circular search, we need to plot the movement of each point
            start_points = np.array(start['points'])
            end_points = np.array(end['points'])
            
            # Plot arrows for each moving point
            for j in range(len(start_points)):
                dx = end_points[j, 0] - start_points[j, 0]
                dy = end_points[j, 1] - start_points[j, 1]
                if abs(dx) > 1e-6 or abs(dy) > 1e-6:  # Only plot if point moved
                    ax.arrow(start_points[j, 0], start_points[j, 1], dx, dy,
                            head_width=1, head_length=2, fc='green', ec='green',
                            length_includes_head=True, alpha=0.6, gid='SEARCH_PATH')

    ax.set_aspect('equal', adjustable='datalim')
    ax.grid(False)
    if show_title and highlight_fs and fs_cache:
        critical_fs = fs_cache[0]['FS']
        ax.set_title(f"Critical Factor of Safety = {critical_fs:.3f}")

    fig.tight_layout()
    _legend_below(ax, fig, legend_ncol=legend_ncol, frameon=legend_frame, show_legend=show_legend)

    if save_png:
        fig.savefig('plot_noncircular_search_results.png', dpi=dpi, bbox_inches='tight')
    if save_dxf:
        from .cad import axes_to_dxf
        axes_to_dxf(ax, 'plot_noncircular_search_results.dxf')

    if own_fig:
        plt.show()
    return fig

def plot_reliability_results(slope_data, reliability_data, figsize=(12, 7), save_png=False, save_dxf=False, dpi=300, legend_ncol="auto", legend_frame=False, show_title=True, show_legend=True, fig=None, style=None):
    """
    Creates a plot showing the results of reliability analysis.

    Parameters:
        slope_data: Dictionary containing plot data
        reliability_data: Dictionary containing reliability analysis results
        figsize: Tuple of (width, height) in inches for the plot
        fig: Optional existing Matplotlib Figure to draw into (used for embedding in a
            GUI canvas). When None (default) a new pyplot figure is created and shown;
            when provided, the figure is cleared and reused and plt.show() is skipped.

    Returns:
        The Matplotlib Figure that was drawn into.
    """
    own_fig = fig is None
    if own_fig:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig.clear()
        ax = fig.add_subplot(111)

    # Plot basic slope elements (same as other search functions)
    plot_base_geometry(ax, slope_data, style=style)
    if any(m.get('u') == 'piezo' for m in slope_data.get('materials', [])):
        plot_piezo_line(ax, slope_data, style=style)
    plot_dloads(ax, slope_data, style=style)
    plot_tcrack_surface(ax, slope_data, style=style)

    # Plot reliability-specific failure surfaces
    fs_cache = reliability_data['fs_cache']
    
    # Plot all failure surfaces
    added_plus_legend = False
    added_minus_legend = False
    for i, fs_data in enumerate(fs_cache):
        result = fs_data['result']
        name = fs_data['name']
        failure_surface = result['failure_surface']

        # Convert failure surface to coordinates
        if hasattr(failure_surface, 'coords'):
            coords = list(failure_surface.coords)
        else:
            coords = failure_surface

        x_coords = [pt[0] for pt in coords]
        y_coords = [pt[1] for pt in coords]

        # Color and styling based on surface type
        if name == "MLV":
            # Highlight the MLV (critical) surface in red
            ax.plot(x_coords, y_coords, color='red', linewidth=3,
                   label=f'$F_{{MLV}}$ Surface (FS={result["FS"]:.3f})', zorder=10)
        else:
            # Other surfaces in different colors
            if '+' in name:
                color = 'blue'
                alpha = 0.7
                label = '$F^+$ surfaces' if not added_plus_legend else None
                added_plus_legend = True
            else:  # '-' in name
                color = 'green'
                alpha = 0.7
                label = '$F^-$ surfaces' if not added_minus_legend else None
                added_minus_legend = True

            ax.plot(x_coords, y_coords, color=color, linewidth=1.5,
                   alpha=alpha, label=label, zorder=5)



    # Standard finalization
    ax.set_aspect('equal', adjustable='datalim')
    ax.grid(False)

    # Title with reliability statistics using mathtext
    F_MLV = reliability_data['F_MLV']
    sigma_F = reliability_data['sigma_F']
    COV_F = reliability_data['COV_F']
    reliability = reliability_data['reliability']
    prob_failure = reliability_data['prob_failure']

    if show_title:
        ax.set_title(f"Reliability Analysis Results\n"
                    f"$F_{{MLV}}$ = {F_MLV:.3f}, $\\sigma_F$ = {sigma_F:.3f}, "
                    f"$COV_F$ = {COV_F:.3f}\n"
                    f"Reliability = {reliability*100:.2f}%, $P_f$ = {prob_failure*100:.2f}%")

    fig.tight_layout()
    _legend_below(ax, fig, legend_ncol=legend_ncol, frameon=legend_frame, show_legend=show_legend)

    if save_png:
        filename = 'plot_reliability_results.png'
        fig.savefig(filename, dpi=dpi, bbox_inches='tight')

    if save_dxf:
        from .cad import axes_to_dxf
        axes_to_dxf(ax, 'plot_reliability_results.dxf')

    if own_fig:
        plt.show()
    return fig

def plot_mesh(mesh, materials=None, figsize=(12, 7), pad_frac=0.05, show_nodes=True, label_elements=False, label_nodes=False, save_png=False, save_dxf=False, dpi=300, legend_ncol="auto", legend_frame=False, show_title=True, show_legend=True, fig=None, style=None):
    """
    Plot the finite element mesh with material regions.

    Parameters:
        mesh: Mesh dictionary with 'nodes', 'elements', 'element_types', and 'element_materials' keys
        materials: Optional list of material dictionaries for legend labels
        figsize: Figure size tuple
        pad_frac: Fraction of mesh size to use for padding around plot
        show_nodes: If True, plot points at node locations
        label_elements: If True, label each element with its number at its centroid
        label_nodes: If True, label each node with its number
        fig: Optional existing Matplotlib Figure to draw into (used for embedding in a
            GUI canvas). When None (default) a new pyplot figure is created and shown;
            when provided, the figure is cleared and reused and plt.show() is skipped.
    """
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch
    from matplotlib.collections import PolyCollection
    import numpy as np

    nodes = mesh["nodes"]
    elements = mesh["elements"]
    element_types = mesh["element_types"]
    mat_ids = mesh["element_materials"]

    own_fig = fig is None
    if own_fig:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig.clear()
        ax = fig.add_subplot(111)
    
    # Group elements by material ID
    material_elements = {}
    for i, (element, elem_type, mid) in enumerate(zip(elements, element_types, mat_ids)):
        if mid not in material_elements:
            material_elements[mid] = []
        
        # Only process 2D elements (skip 1D elements which have elem_type 2)
        if elem_type == 2:  # Skip 1D elements
            continue
        
        # Use corner nodes to define element boundary (no subdivision needed)
        if elem_type in [3, 6]:  # Triangular elements (linear or quadratic)
            element_coords = [nodes[element[0]], nodes[element[1]], nodes[element[2]]]
        elif elem_type in [4, 8, 9]:  # Quadrilateral elements (linear or quadratic)
            element_coords = [nodes[element[0]], nodes[element[1]], nodes[element[2]], nodes[element[3]]]
        else:
            continue  # Skip unknown element types
        
        material_elements[mid].append(element_coords)
    
    legend_elements = []
    
    # Plot 1D elements FIRST (bottom layer) if present in mesh
    if "elements_1d" in mesh and "element_types_1d" in mesh and "element_materials_1d" in mesh:
        elements_1d = mesh["elements_1d"]
        element_types_1d = mesh["element_types_1d"]
        mat_ids_1d = mesh["element_materials_1d"]
        
        # Group 1D elements by material ID
        material_lines = {}
        for i, (element_1d, elem_type_1d, mid_1d) in enumerate(zip(elements_1d, element_types_1d, mat_ids_1d)):
            if mid_1d not in material_lines:
                material_lines[mid_1d] = []
            
            # Get line coordinates based on actual number of nodes
            # elem_type_1d contains the number of nodes (2 for linear, 3 for quadratic)
            if elem_type_1d == 2:  # Linear 1D element (2 nodes)
                # Skip zero-padded elements
                if element_1d[1] != 0:  # Valid second node
                    line_coords = [nodes[element_1d[0]], nodes[element_1d[1]]]
                else:
                    continue  # Skip invalid element
            elif elem_type_1d == 3:  # Quadratic 1D element (3 nodes)
                # For visualization, connect all three nodes or just endpoints
                line_coords = [nodes[element_1d[0]], nodes[element_1d[1]], nodes[element_1d[2]]]
            else:
                continue  # Skip unknown 1D element types
            
            material_lines[mid_1d].append(line_coords)
        
        # Plot 1D elements with distinctive style
        for mid_1d, lines_list in material_lines.items():
            for line_coords in lines_list:
                xs = [coord[0] for coord in line_coords]
                ys = [coord[1] for coord in line_coords]
                ax.plot(xs, ys, color='red', linewidth=3, alpha=0.8, solid_capstyle='round')
        
        # Add 1D elements to legend
        if material_lines:
            legend_elements.append(plt.Line2D([0], [0], color='red', linewidth=3, 
                                            alpha=0.8, label='1D Elements'))
    
    # Material colors (style overrides → palette default). Mesh material IDs are
    # 1-based (gmsh); the style sheet keys by 0-based mat_id, so map mid-1 — this
    # also aligns the zone colors with the Inputs view.
    from .style import resolve_style, material_style
    _st = resolve_style(style)
    def _mat_color(mid):
        return material_style(_st, int(mid) - 1)["color"]

    # Plot 2D elements SECOND (middle layer)
    for mid, elements_list in material_elements.items():
        # Create polygon collection for this material
        poly_collection = PolyCollection(elements_list, gid='MESH',
                                       facecolor=_mat_color(mid),
                                       edgecolor='k',
                                       alpha=0.4,
                                       linewidth=0.5)
        ax.add_collection(poly_collection)
        
        # Add to legend
        if materials and mid <= len(materials) and materials[mid-1].get('name'):
            label = materials[mid-1]['name']  # Convert to 0-based indexing
        else:
            label = f'Material {mid}'
        
        legend_elements.append(Patch(facecolor=_mat_color(mid),
                                   edgecolor='k',
                                   alpha=0.4,
                                   label=label))
    
    # Label 2D elements if requested
    if label_elements:
        for idx, (element, element_type) in enumerate(zip(elements, element_types)):
            # Calculate element centroid based on element type
            if element_type == 3:  # 3-node triangle
                element_coords = nodes[element[:3]]
            elif element_type == 6:  # 6-node triangle
                element_coords = nodes[element[:6]]
            elif element_type == 4:  # 4-node quad
                element_coords = nodes[element[:4]]
            elif element_type == 8:  # 8-node quad
                element_coords = nodes[element[:8]]
            elif element_type == 9:  # 9-node quad
                element_coords = nodes[element[:9]]
            else:
                continue  # Skip unknown element types
            
            centroid = np.mean(element_coords, axis=0)
            ax.text(centroid[0], centroid[1], str(idx+1),
                    ha='center', va='center', fontsize=6, color='black', alpha=0.7,
                    zorder=12)
    
    # Label 1D elements if requested (with different color)
    if label_elements and "elements_1d" in mesh:
        elements_1d = mesh["elements_1d"]
        element_types_1d = mesh["element_types_1d"]
        
        for idx, (element_1d, elem_type_1d) in enumerate(zip(elements_1d, element_types_1d)):
            # Skip zero-padded elements
            if elem_type_1d == 2 and element_1d[1] != 0:  # Linear 1D element
                # Calculate midpoint of line element
                coord1 = nodes[element_1d[0]]
                coord2 = nodes[element_1d[1]]
                midpoint = (coord1 + coord2) / 2
                ax.text(midpoint[0], midpoint[1], f"1D{idx+1}",
                        ha='center', va='center', fontsize=6, color='black', alpha=0.9,
                        zorder=13)
            elif elem_type_1d == 3 and element_1d[2] != 0:  # Quadratic 1D element
                # Use middle node as label position (if it exists)
                midpoint = nodes[element_1d[1]]
                ax.text(midpoint[0], midpoint[1], f"1D{idx+1}",
                        ha='center', va='center', fontsize=6, color='black', alpha=0.9,
                        zorder=13)
    
    # Plot nodes LAST (top layer) if requested
    if show_nodes:
        # Plot all nodes - if meshing is correct, all nodes should be used
        ax.plot(nodes[:, 0], nodes[:, 1], 'k.', markersize=2)
        # Add to legend
        legend_elements.append(plt.Line2D([0], [0], marker='o', color='w', 
                                        markerfacecolor='k', markersize=6, 
                                        label=f'Nodes ({len(nodes)})', linestyle='None'))
    
    # Label nodes if requested
    if label_nodes:
        # Label all nodes
        for i, (x, y) in enumerate(nodes):
            ax.text(x + 0.5, y + 0.5, str(i+1), fontsize=6, color='blue', alpha=0.7,
                    ha='left', va='bottom', zorder=14)
    
    ax.set_aspect('equal', adjustable='datalim')
    if show_title:
        ax.set_title("Finite Element Mesh with Material Regions (Triangles and Quads)")

    # Add cushion
    x_min, x_max = nodes[:, 0].min(), nodes[:, 0].max()
    y_min, y_max = nodes[:, 1].min(), nodes[:, 1].max()
    x_pad = (x_max - x_min) * pad_frac
    y_pad = (y_max - y_min) * pad_frac
    ax.set_xlim(x_min - x_pad, x_max + x_pad)
    ax.set_ylim(y_min - y_pad, y_max + y_pad)

    fig.tight_layout()
    # Add legend if we have materials (after tight_layout so its reserved bottom
    # margin isn't clobbered)
    if legend_elements:
        _legend_below(ax, fig, handles=legend_elements,
                      legend_ncol=legend_ncol, frameon=legend_frame, show_legend=show_legend)

    if save_png:
        fig.savefig('plot_mesh.png', dpi=dpi, bbox_inches='tight')
    if save_dxf:
        from .cad import axes_to_dxf
        axes_to_dxf(ax, 'plot_mesh.dxf')

    if own_fig:
        plt.show()
    return fig


def plot_polygons(
    polygons,
    materials=None,
    nodes=False,
    legend=True,
    title="Material Zone Polygons",
    figsize=(10, 6),
    save_png=False,
    save_dxf=False,
    dpi=300,
):
    """
    Plot all material zone polygons in a single figure.
    
    Parameters:
        polygons: List of polygon coordinate lists or dicts with "coords"/"mat_id"
        materials: Optional list of material dicts (with key "name") or list of material
            name strings. If provided, the material name will be used in the legend.
        nodes: If True, plot each polygon vertex as a dot.
        legend: If True, show the legend.
        title: Plot title
        figsize: Matplotlib figure size tuple, e.g. (10, 6)
    """
    import matplotlib.pyplot as plt
    
    fig, ax = plt.subplots(figsize=figsize)
    for i, polygon in enumerate(polygons):
        coords = polygon.get("coords", []) if isinstance(polygon, dict) else polygon
        xs = [x for x, y in coords]
        ys = [y for x, y in coords]
        mat_idx = polygon.get("mat_id") if isinstance(polygon, dict) else i
        if mat_idx is None:
            mat_idx = i
        mat_name = None
        if materials is not None and 0 <= mat_idx < len(materials):
            item = materials[mat_idx]
            if isinstance(item, dict):
                mat_name = item.get("name", None)
            elif isinstance(item, str):
                mat_name = item
        label = mat_name if mat_name else f"Material {mat_idx}"
        ax.fill(xs, ys, color=get_material_color(mat_idx), alpha=0.6, label=label, gid=label)
        ax.plot(xs, ys, color=get_material_color(mat_idx), linewidth=1, gid=label)
        if nodes:
            # Avoid legend clutter by not adding a label here.
            ax.scatter(xs, ys, color='k', s=30, marker='o', zorder=3)
    ax.set_xlabel('X Coordinate')
    ax.set_ylabel('Y Coordinate')
    ax.set_title(title)
    if legend:
        ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_aspect('equal', adjustable='datalim')
    plt.tight_layout()
    
    base_name = 'plot_' + title.lower().replace(' ', '_').replace(':', '').replace(',', '')
    if save_png:
        plt.savefig(base_name + '.png', dpi=dpi, bbox_inches='tight')
    if save_dxf:
        from .cad import axes_to_dxf
        axes_to_dxf(ax, base_name + '.dxf')

    plt.show()


def plot_polygons_separately(polygons, materials=None, save_png=False, dpi=300):
    """
    Plot each polygon in a separate matplotlib frame (subplot), with vertices as round dots.
    
    Parameters:
        polygons: List of polygon coordinate lists or dicts with "coords"/"mat_id"
        materials: Optional list of material dicts (with key "name") or list of material
            name strings. If provided, the material name will be included in each subplot title.
    """
    import matplotlib.pyplot as plt
    
    n = len(polygons)
    fig, axes = plt.subplots(n, 1, figsize=(8, 3 * n), squeeze=False)
    for i, polygon in enumerate(polygons):
        coords = polygon.get("coords", []) if isinstance(polygon, dict) else polygon
        xs = [x for x, y in coords]
        ys = [y for x, y in coords]
        ax = axes[i, 0]
        mat_idx = polygon.get("mat_id") if isinstance(polygon, dict) else i
        if mat_idx is None:
            mat_idx = i
        ax.fill(xs, ys, color=get_material_color(mat_idx), alpha=0.6, label=f'Material {mat_idx}')
        ax.plot(xs, ys, color=get_material_color(mat_idx), linewidth=1)
        ax.scatter(xs, ys, color='k', s=30, marker='o', zorder=3, label='Vertices')
        ax.set_xlabel('X Coordinate')
        ax.set_ylabel('Y Coordinate')
        mat_name = None
        if materials is not None and 0 <= mat_idx < len(materials):
            item = materials[mat_idx]
            if isinstance(item, dict):
                mat_name = item.get("name", None)
            elif isinstance(item, str):
                mat_name = item
        if mat_name:
            ax.set_title(f'Material {mat_idx}: {mat_name}')
        else:
            ax.set_title(f'Material {mat_idx}')
        ax.grid(True, alpha=0.3)
        ax.set_aspect('equal', adjustable='datalim')
        # Intentionally no legend: these plots are typically used for debugging geometry,
        # and legends can obscure key vertices/edges.
    plt.tight_layout()
    
    if save_png:
        filename = 'plot_polygons_separately.png'
        plt.savefig(filename, dpi=dpi, bbox_inches='tight')
    
    plt.show()


def find_best_table_position(ax, materials, plot_elements_bounds):
    """
    Find the best position for the material table to avoid overlaps.
    
    Parameters:
        ax: matplotlib Axes object
        materials: List of materials to determine table size
        plot_elements_bounds: List of (x_min, x_max, y_min, y_max) for existing elements
        
    Returns:
        (xloc, yloc) coordinates for table placement
    """
    # Calculate table size based on number of materials and columns
    num_materials = len(materials)
    has_d_psi = any(mat.get('d', 0) > 0 or mat.get('psi', 0) > 0 for mat in materials)
    table_height = 0.05 + 0.025 * num_materials  # Height per row
    table_width = 0.25 if has_d_psi else 0.2
    
    # Define candidate positions (priority order) - with margins from borders
    candidates = [
        (0.05, 0.70),  # upper left
        (0.70, 0.70),  # upper right  
        (0.05, 0.05),  # lower left
        (0.70, 0.05),  # lower right
        (0.35, 0.70),  # upper center
        (0.35, 0.05),  # lower center
        (0.05, 0.35),  # center left
        (0.70, 0.35),  # center right
        (0.35, 0.35),  # center
    ]
    
    # Check each candidate position for overlaps
    for xloc, yloc in candidates:
        table_bounds = (xloc, xloc + table_width, yloc - table_height, yloc)
        
        # Check if table overlaps with any plot elements
        overlap = False
        for elem_bounds in plot_elements_bounds:
            elem_x_min, elem_x_max, elem_y_min, elem_y_max = elem_bounds
            table_x_min, table_x_max, table_y_min, table_y_max = table_bounds
            
            # Check for overlap
            if not (table_x_max < elem_x_min or table_x_min > elem_x_max or
                   table_y_max < elem_y_min or table_y_min > elem_y_max):
                overlap = True
                break
        
        if not overlap:
            return xloc, yloc
    
    # If all positions have overlap, return the first candidate
    return candidates[0]


def get_plot_elements_bounds(ax, slope_data):
    """
    Get bounding boxes of existing plot elements to avoid overlaps.
    
    Parameters:
        ax: matplotlib Axes object
        slope_data: Dictionary containing slope data
        
    Returns:
        List of (x_min, x_max, y_min, y_max) tuples for plot elements
    """
    bounds = []
    
    # Get axis limits
    x_min, x_max = ax.get_xlim()
    y_min, y_max = ax.get_ylim()
    
    # Profile lines bounds
    if slope_data.get('profile_lines'):
        for line in slope_data['profile_lines']:
            if line:
                coords = line['coords']
                xs = [p[0] for p in coords]
                ys = [p[1] for p in coords]
                bounds.append((min(xs), max(xs), min(ys), max(ys)))
    elif slope_data.get('domain_polygon') is not None:
        # Polygon input: use the domain polygon's bounding box.
        dx0, dy0, dx1, dy1 = slope_data['domain_polygon'].bounds
        bounds.append((dx0, dx1, dy0, dy1))

    # Distributed loads bounds
    if 'dloads' in slope_data and slope_data['dloads']:
        for dload_set in slope_data['dloads']:
            if dload_set:
                xs = [p['X'] for p in dload_set]
                ys = [p['Y'] for p in dload_set]
                bounds.append((min(xs), max(xs), min(ys), max(ys)))
    
    # Reinforcement lines bounds
    if 'reinforce_lines' in slope_data and slope_data['reinforce_lines']:
        for line in slope_data['reinforce_lines']:
            if line:
                xs = [p['X'] for p in line]
                ys = [p['Y'] for p in line]
                bounds.append((min(xs), max(xs), min(ys), max(ys)))
    
    return bounds

def plot_sensitivity(df, target_fs=None, figsize=(8, 5), save_png=False,
                     dpi=300, fig=None, style=None):
    """Plot a sensitivity() sweep: the output quantity vs parameter value, one
    line per method.

    The base case is marked, and — when the output is a factor of safety — FS = 1
    gets a horizontal guide (plus an optional target). If the sweep re-searched
    the critical surface, points where the surface JUMPED between neighbours (a
    different failure mode taking over) are drawn open — the Xo/Yo/R columns make
    the jump detectable in the data, not just suspected.

    The y-axis / title follow the swept quantity: 'Factor of Safety' for an LEM or
    FEM sweep, 'Total discharge, q' for a seepage sweep (the df's 'output' /
    'output_label' columns carry it; absent columns default to FS for backward
    compatibility, so old sweeps plot exactly as before).

    Parameters:
        df: the DataFrame from sensitivity() (result['df']).
        target_fs: optional design target to draw as a second guide line (an FS in
            LEM/FEM, a discharge q in a seepage sweep).
    """
    import numpy as np
    if fig is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        ax = fig.add_subplot(111)

    param = df['param'].iloc[0]
    # output-quantity awareness (backward-compatible defaults)
    output = df['output'].iloc[0] if 'output' in df.columns and len(df) else 'FS'
    output_label = (df['output_label'].iloc[0]
                    if 'output_label' in df.columns and len(df)
                    else 'Factor of Safety')
    is_fs = output == 'FS'
    for method, g in df.groupby('method'):
        pts = g.loc[~g['is_base'] & g['success']].sort_values('value')
        ax.plot(pts['value'], pts['fs'], marker='o', label=method)
        # surface-jump detection: a step in (Xo, Yo, R) far larger than its
        # neighbours' median step means a different critical surface took over
        if pts[['Xo', 'Yo', 'R']].notna().all().all() and len(pts) >= 3:
            geo = pts[['Xo', 'Yo', 'R']].to_numpy()
            steps = np.linalg.norm(np.diff(geo, axis=0), axis=1)
            med = np.median(steps)
            if med > 0:
                for k in np.nonzero(steps > 5 * med)[0]:
                    ax.plot(pts['value'].iloc[k + 1], pts['fs'].iloc[k + 1],
                            marker='o', mfc='white', mec='C3', ms=10, zorder=5,
                            linestyle='None')
        base = g.loc[g['is_base'] & g['success']]
        if not base.empty and np.isfinite(base['value'].iloc[0]):
            ax.plot(base['value'].iloc[0], base['fs'].iloc[0], 's', color='k',
                    ms=8, zorder=6,
                    label=f"base case ({base['value'].iloc[0]:g}, "
                          f"{output} = {base['fs'].iloc[0]:.3g})")
    if is_fs:                                    # FS = 1 guide only for a safety factor
        ax.axhline(1.0, color='r', linestyle='--', linewidth=0.8, label='FS = 1')
    if target_fs is not None:
        ax.axhline(target_fs, color='gray', linestyle='--', linewidth=0.8,
                   label=f'{output} = {target_fs:g}')
    ax.set_xlabel(param)
    ax.set_ylabel(output_label)
    ax.set_title(f'Sensitivity: {output} vs {param}')
    ax.legend()
    ax.grid(True, alpha=0.4)
    fig.tight_layout()
    if save_png:
        fig.savefig(f"sensitivity_{param.replace(':', '_')}.png", dpi=dpi,
                    bbox_inches='tight')
    return fig


def plot_tornado(result, figsize=(8, 5), save_png=False, dpi=300, fig=None,
                 style=None, widest_on_top=True):
    """Duncan-style tornado diagram from tornado(): horizontal bars of the output
    swing between each parameter's low and high bound, sorted by span, with the
    base-case value as the vertical reference.

    The swept quantity follows result['output'] / result['output_label'] — a
    factor of safety for an LEM/FEM tornado, total discharge q for a seepage one
    (absent keys default to FS, so an old result plots as before).

    Parameters:
        result: the dict from tornado() (carries 'df', 'base_fs', 'method', and
            optionally 'output' / 'output_label').
        widest_on_top: the classic tornado stacking (default). False inverts,
            widest at the bottom.
    """
    import numpy as np
    df, base_fs = result['df'], result['base_fs']
    output = result.get('output', 'FS')
    output_label = result.get('output_label', 'Factor of Safety')
    # FS reads best fixed-point; a tiny discharge q needs significant figures.
    vfmt = (lambda x: f'{x:.2f}') if output == 'FS' else (lambda x: f'{x:.3g}')
    if fig is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        ax = fig.add_subplot(111)

    bars = []
    for param, g in df.loc[~df['is_base']].groupby('param'):
        ok = g.loc[g['success']].sort_values('value')
        if ok.empty:
            continue
        lo_fs, hi_fs = ok['fs'].iloc[0], ok['fs'].iloc[-1]
        bars.append((param, lo_fs, hi_fs, abs(hi_fs - lo_fs)))
    bars.sort(key=lambda b: b[3], reverse=not widest_on_top)  # barh draws k=0 at the bottom

    for k, (param, lo_fs, hi_fs, _span) in enumerate(bars):
        left, width = min(lo_fs, hi_fs), abs(hi_fs - lo_fs)
        ax.barh(k, width, left=left, height=0.55, color='C0', alpha=0.75)
        ax.annotate(vfmt(lo_fs), (lo_fs, k), textcoords='offset points',
                    xytext=(-6, 0), ha='right', va='center', fontsize=8)
        ax.annotate(vfmt(hi_fs), (hi_fs, k), textcoords='offset points',
                    xytext=(6, 0), ha='left', va='center', fontsize=8)
    ax.set_yticks(range(len(bars)))
    ax.set_yticklabels([b[0] for b in bars])
    # Explicit limits, not margins(): the Studio canvas re-fits axes after render,
    # which drops margin padding and lets the widest bar's low/high value labels
    # collide with the tick labels at the spine. Fixed limits survive the re-fit.
    if bars:
        ends = [b[1] for b in bars] + [b[2] for b in bars]
        if base_fs is not None and np.isfinite(base_fs):
            ends.append(base_fs)
        lo_min, hi_max = min(ends), max(ends)
        span = (hi_max - lo_min) or max(abs(hi_max), 1.0)
        ax.set_xlim(lo_min - 0.15 * span, hi_max + 0.15 * span)
    if base_fs is not None and np.isfinite(base_fs):
        ax.axvline(base_fs, color='k', linewidth=1.0,
                   label=f'base {output} = {vfmt(base_fs)}')
        ax.legend()
    # LEM method label only makes sense for a factor of safety.
    method = result.get('method', '')
    xlabel = output_label + (f" ({method})" if method and output == 'FS' else "")
    ax.set_xlabel(xlabel)
    ax.set_title(f'Tornado: {output} swing per parameter')
    ax.grid(True, axis='x', alpha=0.4)
    fig.tight_layout()
    if save_png:
        fig.savefig('tornado.png', dpi=dpi, bbox_inches='tight')
    return fig
