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
from .units import require_gamma_water


def declared_unit_labels(data):
    """The :func:`xslope.units.labels` dict for a ``slope_data`` / ``seep_data`` /
    ``fem_data`` mapping's declared unit system, or ``None`` when it declares none.

    The single place the plotters resolve unit labels. Returning ``None`` (not the
    all-empty dict) lets a call site guard with a plain truthiness check and keeps an
    undeclared model's plots byte/pixel-identical to today: no axis unit, no colorbar
    unit, no flowrate unit. Time-bearing labels (k, flowrate) stay empty until a time
    unit is also declared, so they are never guessed."""
    from .units import labels, normalize_unit_system
    system = normalize_unit_system((data or {}).get("unit_system"))
    if system is None:
        return None
    return labels(system, (data or {}).get("time_unit"))

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

    elif option == "elastic":
        # Elastic / infinite strength (v16): no shear envelope — it cannot fail.
        return _material_hint(ax, "elastic — no strength model\n(cannot fail)",
                              "Elastic")

    else:
        return _material_hint(ax, "no strength model\n(blank option)")

    return ax


def _kr_curve(material, n=200):
    """One material's unsaturated relative-conductivity curve, or why it has none.

    Returns ``(psi, kr, title)`` — matric suction, relative conductivity, and the
    model with its parameters — for a material whose unsaturated model is
    assigned and parameterized, and ``(None, message, title)`` for one that is
    not, where ``message`` is what has to be entered for the curve to exist.

    The numbers come from seep.py's OWN kr functions (imported, never
    re-derived), dispatched on the material's ``unsat`` model:

      - ``lf``   linear front — seep.kr_frontal_vec (params kr0, h0).
      - ``vg``   van Genuchten–Mualem — seep.kr_vg_vec (params vg_a, vg_n).
      - ``gard`` Gardner power form — seep.kr_gardner_vec (params vg_a=a, vg_n=n).

    The suction range is self-scaled to each model so the full wet→dry decline is
    sampled. Pure — a material dict in, arrays out.
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
            return None, "enter kr0 > 0 and h0 < 0", "Linear front (lf)"
        smax = 1.3 * abs(h0)
        psi = np.linspace(0.0, smax, n)
        kr = kr_frontal_vec(-psi, kr0, h0)
        title = f"Linear front   (kr0={kr0:g}, h0={h0:g})"
    elif unsat == "vg":
        a, nn = g("vg_a"), g("vg_n")
        if not (a > 0 and nn > 1):
            return None, "enter vg_a > 0 and vg_n > 1", "van Genuchten (vg)"
        smax = 20.0 / a
        psi = np.linspace(0.0, smax, n)
        kr = kr_vg_vec(-psi, a, nn)
        title = f"van Genuchten   (α={a:g}, n={nn:g})"
    elif unsat == "gard":
        a, nn = g("vg_a"), g("vg_n")
        if not (a > 0 and nn > 0):
            return None, "enter vg_a > 0 and vg_n > 0", "Gardner (gard)"
        # Suction at which kr reaches the floor: a·ψⁿ = 1/kr_min − 1.
        smax = 1.1 * (max(1.0 / kr_min - 1.0, 1.0) / a) ** (1.0 / nn)
        psi = np.linspace(0.0, smax, n)
        kr = kr_gardner_vec(-psi, a, nn)
        title = f"Gardner   (a={a:g}, n={nn:g})"
    else:
        return None, "no unsaturated model", ""
    return psi, kr, title


#: The unsaturated model whose kr is a straight ramp: kr falls linearly from 1 at
#: h = 0 to kr0 at h = h0. It spans no decades — it is a line — and it is the one
#: model of the three that a logarithmic ordinate misrepresents.
_KR_LINEAR_MODELS = ("lf",)


def _kr_model(material):
    """A material's assigned unsaturated model, normalized: ``lf``, ``vg``,
    ``gard``, or whatever else was entered. Unassigned reads as the linear front,
    the template's default and what :func:`_kr_curve` dispatches on."""
    try:
        return str((material or {}).get("unsat", "lf") or "lf").strip().lower()
    except AttributeError:
        return "lf"


def _kr_yscale(materials):
    """The ordinate scale the kr curves of ``materials`` are read on: ``"linear"``
    or ``"log"``.

    The scale follows the models assigned, not the sampled values — the model is
    what the curve means, and a shape read off its own numbers is read off the
    thing in question. The linear front is a straight ramp from kr = 1 to kr0
    over one head interval; drawn against a logarithmic ordinate it becomes a
    decaying curve falling off a cliff, and the front — the whole content of the
    model — is unreadable. van Genuchten and Gardner genuinely fall through
    decades, and a linear ordinate flattens all of that against zero.

    So: every material drawn on the linear front, and the ordinate is linear over
    0 to 1, where the ramp reads as the straight front it is. Any material drawn
    on a model that spans decades, and the ordinate is logarithmic for the whole
    axes — one axes carries one scale, and a linear front drawn on it is still
    the same curve, correctly plotted.
    """
    models = [_kr_model(m) for m in materials or []]
    if models and all(m in _KR_LINEAR_MODELS for m in models):
        return "linear"
    return "log"


def _apply_kr_yscale(ax, scale):
    """Set ``ax``'s kr ordinate to ``scale``, with the limits and gridding that
    scale is read with.

    Linear: 0 to just past 1, the full range kr is defined over, with grid lines
    on the major ticks. Log: bounded above just past 1 and left to autoscale
    below onto whatever floor the curves reach, with the minor decade lines shown
    — on a log axis the decade subdivisions are how a reader places a value.
    """
    if scale == "linear":
        ax.set_yscale("linear")
        ax.set_ylim(0, 1.05)
        ax.grid(True, which="major", alpha=0.3)
    else:
        ax.set_yscale("log")
        ax.set_ylim(top=1.5)
        ax.grid(True, which="both", alpha=0.3)
    return ax


def _kr_extent(psi, kr):
    """Where a kr curve's x axis should end: the suction at which it reaches its
    floor, and no further.

    kr is monotone-decreasing, so the last sampled value is the asymptotic floor
    — kr0 for the linear front, the solver's kr_min for the others. Framing on
    where that floor is first reached keeps a sharp curve out of a sea of dead
    space.
    """
    floor = float(kr[-1])
    reached = np.where(kr <= floor * 1.02 + 1e-12)[0]
    right = psi[reached[0]] * 1.1 if len(reached) else psi[-1]
    return float(right if right > 0 else psi[-1])


def material_kr_curves(materials, n=200):
    """``[(index, material, psi, kr, model), ...]`` for every material that carries
    a drawable unsaturated curve.

    ``index`` is the material's 1-based position in ``materials`` — the number the
    report and the model figure key it by. A material assigned no unsaturated
    model, or one whose parameters are missing, is simply absent: it is analyzed
    saturated, and there is no curve to draw for it.
    """
    out = []
    for i, material in enumerate(materials or []):
        psi, kr, title = _kr_curve(material, n)
        if psi is not None:
            out.append((i + 1, material, psi, kr, title))
    return out


def plot_material_kr(ax, material, n=200):
    """Draw one material's unsaturated relative-conductivity curve kr vs matric
    suction into ``ax``.

    The ordinate scale is chosen by :func:`_kr_yscale` from the model the
    material is assigned — linear for the linear front, logarithmic for a model
    that spans decades — so the editor draws a curve on the same axis the report
    draws it on. The suction range is framed by :func:`_kr_extent`. Missing or
    invalid parameters draw a centered hint. Pure — Axes + material dict only.
    """
    psi, kr, title = _kr_curve(material, n)
    if psi is None:
        return _material_hint(ax, kr, title)

    ax.plot(psi, kr, color=_KR_COLOR, lw=2)
    _apply_kr_yscale(ax, _kr_yscale([material]))
    ax.set_xlim(0, _kr_extent(psi, kr))
    ax.set_xlabel("matric suction, ψ")
    ax.set_ylabel("relative conductivity, kr")
    ax.set_title(title)
    return ax


#: Line styles the combined kr plot cycles through. Two materials given the same
#: unsaturated parameters draw the same curve, and a solid line laid over a solid
#: line hides it; a dashed one laid over a solid one does not.
_KR_DASHES = ["-", "--", "-.", ":"]


def plot_material_kr_set(materials, n=200, fig=None, figsize=(8, 5), style=None,
                         show_title=True, show_legend=True, unit_labels=None,
                         legend_ncol="auto", legend_frame=False,
                         abscissa="suction"):
    """Every material's unsaturated conductivity curve on one axes.

    The same curves the material editor draws one at a time, together: what the
    flow solver reduces each material's conductivity to above the phreatic
    surface. Each curve takes its material's own zone color, so a curve is found
    in the legend and again in the model figure.

    ``abscissa`` is the convention the same curves are read in:

      - ``"suction"`` (default): matric suction ψ, positive and increasing to
        the right as the material dries.
      - ``"head"``: pressure head h = −ψ, negative through the unsaturated range
        and running to zero at the right-hand edge, where the material is
        saturated and kr = 1.

    They are one set of numbers under two sign conventions — the curve of a
    material is the same curve in both, mirrored about the ordinate — so the two
    are drawn from the same :func:`material_kr_curves` evaluation, in the same
    order, and a material keeps its color, its dash pattern, and its ordinate
    scale across both.

    That scale comes from :func:`_kr_yscale`, off the models the drawn materials
    are assigned: linear when they are all the linear front, logarithmic as soon
    as one of them spans decades.

    Materials with no unsaturated model are left out — they are analyzed
    saturated and have no curve — and a set with none of them draws nothing and
    returns None, so a caller can ask for the figure and be told there isn't one.

    ``unit_labels`` is the model's declared units, from
    :func:`declared_unit_labels`; both abscissae are a head, and carry the
    length unit.
    """
    curves = material_kr_curves(materials, n)
    if not curves:
        return None
    head = str(abscissa).strip().lower() == "head"

    own_fig = fig is None
    if own_fig:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig.clear()
        ax = fig.add_subplot(111)

    from .style import resolve_style, material_style
    st = resolve_style(style)

    right = 0.0
    for k, (index, material, psi, kr, _model) in enumerate(curves):
        name = str(material.get("name") or "").strip() or f"Material {index}"
        # Color by the material's own index and dash by position in the set —
        # both read off the same curve list for either abscissa, so the pair of
        # figures keys a material identically.
        ax.plot(-psi if head else psi, kr, label=name,
                color=material_style(st, index - 1)["color"],
                ls=_KR_DASHES[k % len(_KR_DASHES)], lw=2)
        right = max(right, _kr_extent(psi, kr))

    # One scale for the axes, chosen from the models on it and from nothing about
    # the abscissa — the pair is the same curves under two sign conventions, and
    # a pair read on two ordinates does not mirror.
    _apply_kr_yscale(ax, _kr_yscale([m for _i, m, *_rest in curves]))
    if right > 0:
        ax.set_xlim(-right, 0) if head else ax.set_xlim(0, right)
    length = (unit_labels or {}).get("length") or ""
    unit = f" ({length})" if length else ""
    ax.set_xlabel(("pressure head, h" if head else "matric suction, ψ") + unit)
    ax.set_ylabel("relative conductivity, $k_r$")
    if show_title:
        ax.set_title("Unsaturated relative conductivity")
    _legend_below(ax, fig, legend_ncol=legend_ncol, frameon=legend_frame,
                  show_legend=show_legend)
    if own_fig:
        plt.show()
    return fig

def _derived_water_blocks(slope_data):
    """The derived water blocks a model would be drawn with (empty in manual mode).

    Asked by the legend and the framing, which must account for a load the user
    never typed; it goes through the same derivation the plot and the solver use.
    """
    from .water import with_water_loads, derived_blocks
    sd = with_water_loads(slope_data)
    return derived_blocks(sd, 1) + derived_blocks(sd, 2)


def _derived_water_color(style=None):
    """The stage-1 derived water-load colour, for the legend key."""
    from .style import resolve_style, feature_style
    return feature_style(resolve_style(style), "dloads_derived").get(
        "color", "royalblue")


def get_dload_legend_handler(color='purple'):
    """
    Creates and returns a custom legend entry for distributed loads.
    Returns a tuple of (handler_class, dummy_patch) for use in matplotlib legends.

    ``color`` selects which load set the key stands for: the user's own blocks
    (purple, the default) or the engine's derived water loads (water-blue).
    """
    # Create a line with built-in arrow marker
    dummy_line = Line2D([0.0, 1.0], [0, 0],  # Two points to define line
                       color=color, 
                       alpha=0.7, 
                       linewidth=2,
                       marker='>',  # Built-in right arrow marker
                       markersize=6,  # Smaller marker size
                       markerfacecolor=color,
                       markeredgecolor=color,
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

        # The legend entry carries the material the line bounds: a figure of
        # numbered profile lines and a table of named materials otherwise leave a
        # reader to work out which is which, and nothing on either says.
        name = ""
        if materials and 0 <= mat_idx < len(materials):
            name = str((materials[mat_idx] or {}).get("name") or "").strip()
        label = f'Profile {i+1}' + (f' ({name})' if name else '')

        ax.plot(xs, ys, color=color, linewidth=fs.get("linewidth", 1.0),
                linestyle=fs.get("linestyle", "-"),
                label=label, gid=f'PROFILE_{i+1}')

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


# SSR zone kinds -> (style feature name, legend label). The labels are the
# canonical display triad, identical to the template's row-6 echo formulas and to
# the loader/Studio wording.
SSR_ZONE_FEATURES = {
    'reduce': ('ssr_zone_reduce', 'SSR reduce'),
    'hold': ('ssr_zone_hold', 'SSR hold'),
    'hold_elastic': ('ssr_zone_elastic', 'SSR elastic'),
}


def plot_ssr_zones(ax, slope_data, style=None):
    """Draw the v20 SSR zone overlays: a dashed boundary over a low-alpha wash.

    SSR zones are analysis overlays, not geometry — they are never meshed and never
    become material regions — so they are drawn as annotation: above the material
    fills, below every line and label (zorder 1.5, between the Patch and Line2D
    defaults). One legend entry per KIND, using the canonical codes "SSR reduce" /
    "SSR hold" / "SSR elastic"; repeated rows of the same kind share it.

    Parameters:
        ax: matplotlib Axes object
        slope_data: dict; reads 'ssr_zones' (absent/empty → nothing is drawn)
        style: optional style sheet (see xslope.style); None → defaults
    """
    zones = slope_data.get('ssr_zones') or []
    if not zones:
        return

    from matplotlib.colors import to_rgba
    from matplotlib.patches import Polygon as MplPolygon
    from .style import resolve_style, feature_style
    style = resolve_style(style)

    seen = set()
    for zone in zones:
        kind = str(zone.get('kind', '')).strip()
        feature, label = SSR_ZONE_FEATURES.get(kind, (None, None))
        if feature is None:
            continue
        coords = [(float(x), float(y)) for x, y in (zone.get('polygon') or [])]
        if len(coords) < 3:
            continue
        fs = feature_style(style, feature)
        color = fs.get('color', 'black')
        legend_label = label if label not in seen else None
        seen.add(label)
        ax.add_patch(MplPolygon(
            coords, closed=True,
            facecolor=to_rgba(color, fs.get('alpha', 0.10)),
            edgecolor=color, linestyle=fs.get('linestyle', '--'),
            linewidth=fs.get('linewidth', 1.6),
            zorder=1.5, label=legend_label, gid=f'SSR_ZONE_{kind.upper()}'))


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


def _label_geometry(slope_data):
    """``(points, segments)`` — the vertices a coordinate label names and the
    lines they are corners of.

    ``points`` is every unique vertex of the profile lines (or of the material
    zone polygon exteriors), plus the two ends of the max-depth line where one
    is drawn. ``segments`` is the same geometry as ``((x0, y0), (x1, y1))``
    pairs: what the labels have to stay off, and what the corner each label
    names is a corner OF.
    """
    seen = set()
    points = []
    segments = []

    def _add(x, y):
        key = (round(float(x), 6), round(float(y), 6))
        if key not in seen:
            seen.add(key)
            points.append(key)
        return key

    def _run(coords):
        keys = [_add(x, y) for x, y in coords]
        for a, b in zip(keys, keys[1:]):
            if a != b:
                segments.append((a, b))

    if slope_data.get('profile_lines'):
        for line in slope_data['profile_lines']:
            _run(line['coords'])
        md = slope_data.get('max_depth')
        if md is not None:
            xs = [x for line in slope_data['profile_lines'] for x, _ in line['coords']]
            _run([(min(xs), md), (max(xs), md)])
    elif slope_data.get('polygons'):
        for poly in slope_data['polygons']:
            coords = list(poly['polygon'].exterior.coords)
            if len(coords) > 1 and coords[0] == coords[-1]:
                coords = coords[:-1]
            _run(coords + coords[:1])
    return points, segments


def _box_crosses_segment(box, p0, p1):
    """Does the segment ``p0``-``p1`` pass through the rectangle ``box``?

    Liang-Barsky clip: the segment is inside for some parameter in [0, 1] iff
    the entering parameter never passes the leaving one.
    """
    x0, y0 = p0
    x1, y1 = p1
    dx, dy = x1 - x0, y1 - y0
    t0, t1 = 0.0, 1.0
    for p, q in ((-dx, x0 - box.x0), (dx, box.x1 - x0),
                 (-dy, y0 - box.y0), (dy, box.y1 - y0)):
        if p == 0:
            if q < 0:
                return False
        else:
            r = q / p
            if p < 0:
                if r > t1:
                    return False
                t0 = max(t0, r)
            else:
                if r < t0:
                    return False
                t1 = min(t1, r)
    return t0 <= t1


def _segments_cross(a, b):
    """Do two open segments, each ``((x0, y0), (x1, y1))``, properly cross?

    Used on leader lines, where a shared endpoint is not a crossing — two labels
    hung off the same vertex are not confusable — so the test is on strict sign
    changes of both orientations.
    """
    (ax0, ay0), (ax1, ay1) = a
    (bx0, by0), (bx1, by1) = b

    def side(px, py, qx, qy, rx, ry):
        return (qx - px) * (ry - py) - (qy - py) * (rx - px)

    d1 = side(ax0, ay0, ax1, ay1, bx0, by0)
    d2 = side(ax0, ay0, ax1, ay1, bx1, by1)
    d3 = side(bx0, by0, bx1, by1, ax0, ay0)
    d4 = side(bx0, by0, bx1, by1, ax1, ay1)
    return (d1 * d2 < 0) and (d3 * d4 < 0)


def _outward_angle(here, neighbours):
    """The direction, in radians, of the widest empty wedge at a vertex.

    A coordinate label belongs in the space the section does NOT occupy, and at
    a corner that space is the largest gap between the edges meeting there: at a
    crest corner it points up and away from the fill, at the toe of a slope it
    points down and out, and on a straight run it points square off the line.
    Read off the geometry itself, so no model needs a placement rule of its own.
    """
    angles = sorted(math.atan2(ny - here[1], nx - here[0])
                    for nx, ny in neighbours)
    if not angles:
        return math.pi / 2                       # nothing to lean away from: up
    if len(angles) == 1:
        return angles[0] + math.pi               # straight back down the edge
    gaps = [(angles[(i + 1) % len(angles)] - a
             + (2 * math.pi if i + 1 == len(angles) else 0.0), a)
            for i, a in enumerate(angles)]
    width, start = max(gaps)
    return start + width / 2


def plot_coordinate_labels(ax, slope_data, fontsize=7, arrows=False, style=None):
    """Annotate the geometry's vertex coordinates, verification-manual style.

    Labels every unique vertex of the profile lines (or material-zone polygon
    exteriors) with "(x, y)". Vertices shared by several lines — ground points
    repeated on each layer's profile line, coincident polygon corners — are
    labelled once. Values print with %g so integers stay clean.

    Placement is solved, not tabulated, and it is solved NEAR. A label belongs
    beside the point it names — a character's width off it, the least a box can
    stand clear and still be a separate mark — and only company pushes it out.
    So the candidate offsets are a ring of directions at a series of distances
    walked from the tightest outward, and the search takes the first distance at
    which some direction is unobjectionable: no box over another label, none off
    the edge of the axes, none across a geometry line, none nearer another
    labelled corner than its own, no leader across another's. Among the
    directions that clear at that distance it takes the one pointing into the
    widest empty wedge at the corner (:func:`_outward_angle` — the direction the
    section is not). Nothing a reader can work around is traded against distance,
    which is why an uncrowded vertex gets its label at arm's length rather than
    at the end of the compromise a whole-figure cost would have reached.

    Only where the whole ring is objectionable — a crest cluster with four
    vertices inside one label width — does the cost of a position get weighed:
    an overlap against a departure from the frame against a misattribution
    against a crossed line, with distance priced as a square so the least-bad
    position is also the nearest one that is that good. Every term is measured
    off the figure being drawn, so a crowded dam crest and a bare slope are laid
    out by the same rule and neither needs a number of its own.

    Parameters:
        ax: matplotlib Axes object
        slope_data: dict with 'profile_lines' or 'polygons'
        fontsize: label font size in points (default 7)
        arrows: if True, tie EVERY label to its vertex with a thin gray leader
            line; if False (default) only the labels a collision pushed clear of
            their vertex get one, since a label sitting beside its point needs no
            leader and one that had to travel is unattributable without it.
            Either way a leader stops half a character short of the vertex it
            points at rather than meeting it, the way a callout is drawn
        style: reserved for future style-sheet control (unused)
    """
    from matplotlib.font_manager import FontProperties
    from matplotlib.transforms import Bbox

    points, segments = _label_geometry(slope_data)
    if not points:
        return

    fig = ax.figure
    try:
        fig.canvas.draw()          # settle limits; obtain a renderer
        renderer = fig.canvas.get_renderer()
    except Exception:
        renderer = None
    if renderer is None:
        # No renderer, no measurements: fall back to a fixed offset rather than
        # laying out blind. Nothing is drawn worse than it was before.
        for x, y in points:
            ax.annotate(f"({x:g}, {y:g})", (x, y), textcoords="offset points",
                        xytext=(6, 6), fontsize=fontsize, color="black",
                        ha="left", va="bottom",
                        bbox=dict(boxstyle="round,pad=0.15", facecolor="white",
                                  edgecolor="none", alpha=0.8),
                        zorder=6, gid="COORD_LABEL")
        return

    prop = FontProperties(size=fontsize)
    scale = fig.dpi / 72.0         # offset points -> pixels
    line_pt = fontsize * 1.4       # one line of this text, in points

    # Everything below is measured in display pixels: the label boxes, the
    # geometry lines they must not cross, and the frame they must stay inside.
    pixels = {p: tuple(ax.transData.transform(p)) for p in points}
    seg_px = [(pixels[a], pixels[b]) for a, b in segments]
    frame = ax.get_window_extent(renderer)

    # Text already on the axes is an obstacle like any placed label: the profile
    # numbers sit on the lines whose ends these labels name, and a coordinate
    # printed over one hides the layer it identifies.
    from matplotlib.text import Text as _Text
    standing = []
    for artist in list(ax.texts):
        try:
            standing.append(_Text.get_window_extent(artist, renderer))
        except Exception:
            pass

    neighbours = {p: [] for p in points}
    for a, b in segments:
        neighbours[a].append(b)
        neighbours[b].append(a)
    outward = {p: _outward_angle(p, neighbours[p]) for p in points}

    # One character of this text, in points: the least a label can stand off the
    # point it names and still read as a mark of its own rather than as part of
    # the drawing. Measured from the font in use, so it is a character's width at
    # whatever size the caller asked for.
    char_pt = renderer.get_text_width_height_descent("0", prop, False)[0] / scale

    # A leader stops just short of the corner it points at instead of running
    # into it: a callout line that meets its vertex reads as part of the section
    # rather than as a line drawn to it, and at a crest the two are the same
    # ink. The gap is half a character — the same measurement the tightest ring
    # is made of, so it is a hair beside the type at any font size, and being in
    # points it stays that hair at any dpi rather than growing with the pixels.
    leader_gap = char_pt * 0.75

    # What the label end of a leader costs it: the rounded box the coordinate is
    # printed in (its pad is in font units) and the point the arrow is already
    # pulled back from that box. Both come off the run before the vertex-end gap
    # does, which is how a leader knows whether anything is left of it to draw.
    BOX_PAD = 0.15 * fontsize
    SHRINK_A = 1.0

    # The ring: sixteen directions at seven distances, walked from the tightest
    # outward — a character's width off the point, then a line of text further
    # each time. Both are read off the type, so a bigger font simply spreads the
    # same layout further.
    angles = [math.radians(a) for a in range(0, 360, 360 // 16)]
    radii = [char_pt + line_pt * f for f in (0.0, 0.4, 0.9, 1.6, 2.6, 3.8, 5.2)]

    # The breathing room a label keeps around its text, in pixels: what it wants
    # clear of its neighbours is a little more than the ink itself.
    PAD = 2.0

    def box_at(px, py, w, h, dx, dy):
        """``(ha, va, clear box, text box)``, in pixels, for an offset of
        ``(dx, dy)`` points.

        The two boxes are the same rectangle at two sizes. What the label needs
        kept clear is the text plus its breathing room; what a reader SEES, and
        so what says which corner the coordinate belongs beside, is the text
        alone.
        """
        ha = "left" if dx > 2 else ("right" if dx < -2 else "center")
        va = "bottom" if dy > 2 else ("top" if dy < -2 else "center")
        ox, oy = px + dx * scale, py + dy * scale
        x0 = ox - (w if ha == "right" else w / 2 if ha == "center" else 0)
        y0 = oy - (h if va == "top" else h / 2 if va == "center" else 0)
        return (ha, va,
                Bbox.from_bounds(x0 - PAD, y0 - PAD, w + 2 * PAD, h + 2 * PAD),
                Bbox.from_bounds(x0, y0, w, h))

    # The most constrained points are placed first — the ones with the most
    # company inside a label's own reach — so a crest cluster picks its wedges
    # before an isolated toe point that could go anywhere claims one.
    reach = radii[-1] * scale

    def crowding(p):
        px, py = pixels[p]
        return sum(1 for q in points
                   if q != p and abs(pixels[q][0] - px) < reach
                   and abs(pixels[q][1] - py) < reach)

    order = sorted(points, key=lambda p: (-crowding(p), p[0], -p[1]))

    def leader_of(point, bb):
        """The run a reader's eye makes: the point, to the middle of its box."""
        return (pixels[point], (bb.x0 + bb.width / 2, bb.y0 + bb.height / 2))

    def _gap(bb, px, py):
        """How far a point lies from a box, in pixels; zero inside it."""
        dx = max(bb.x0 - px, 0.0, px - bb.x1)
        dy = max(bb.y0 - py, 0.0, py - bb.y1)
        return math.hypot(dx, dy)

    def misread(point, bb):
        """Is this box nearer some OTHER labelled vertex than its own?

        A label with no leader is read as belonging to whatever it sits beside,
        so along a stack of vertices a label-height apart — the three right-hand
        ends of a dam's profile lines, one above the next — every label can land
        beside its neighbour's point without one leader meeting another. Nothing
        in the drawing then says the coordinates are shifted by one. The label
        has to be nearest what it names, or be tied to it.
        """
        own = _gap(bb, *pixels[point])
        return any(_gap(bb, *pixels[q]) < own for q in points if q != point)

    def offframe(tb):
        """Is any of the label's text outside the axes?

        Measured on the text and not on its breathing room: a coordinate may sit
        right up against the frame, and on a section whose corners ARE the frame
        that is the only place several of them can go.
        """
        return not (frame.x0 <= tb.x0 and tb.x1 <= frame.x1
                    and frame.y0 <= tb.y0 and tb.y1 <= frame.y1)

    def readability(point, bb, tb, r, a, occupied):
        """What a box at this offset costs a READER, crossings aside.

        ``occupied`` is the text already on the axes and the boxes of the other
        labels this one has to live with. Every term is an area or a count
        measured on the figure being drawn, and the scale of each says how much
        worse it is than the next: a label under another label is unreadable, a
        label off the frame is not there at all, and a label over a line merely
        hides the section.
        """
        area = max(bb.width * bb.height, 1.0)
        hits = [max(0.0, min(bb.x1, pb.x1) - max(bb.x0, pb.x0))
                * max(0.0, min(bb.y1, pb.y1) - max(bb.y0, pb.y0))
                for pb in occupied]
        cost = sum(5000.0 + 2000.0 * over / area for over in hits if over > 0)
        outside = (area - max(0.0, min(bb.x1, frame.x1) - max(bb.x0, frame.x0))
                   * max(0.0, min(bb.y1, frame.y1) - max(bb.y0, frame.y0)))
        cost += 9000.0 * max(0.0, outside) / area
        cost += 90.0 * min(sum(_box_crosses_segment(bb, *s) for s in seg_px), 3)
        # Sitting nearer another vertex than its own: a coordinate read against
        # the wrong corner. Priced under an overlap, which hides a label
        # outright, and over everything a reader can still work around.
        cost += 2500.0 * misread(point, tb)
        # Distance from the point it names, priced as a square: among positions
        # that are equally readable the nearest is wanted, and a label two lines
        # out is more than twice as adrift as one a line out, not twice.
        cost += 20.0 * (r / line_pt) ** 2
        # Departure from the direction that corner has room in.
        cost += 12.0 * (1.0 - math.cos(a - outward[point]))
        return cost

    #: What one leader crossing another costs. Two leaders that cross swap the
    #: labels for a reader — the coordinate at the end of a line is read back
    #: along the WRONG line — so no arrangement that reads correctly is ever
    #: worth trading for one that does not, whatever else it gains. The weight
    #: is above the sum of every other term rather than merely large, which
    #: makes non-crossing a constraint the search cannot outbid and leaves it a
    #: cost rather than a hard filter, so a geometry that admits no crossing-free
    #: layout still gets the layout with the fewest.
    CROSSING = 1.0e6

    def objectionable(point, bb, tb, lead, occupied, leaders):
        """Is there anything wrong with this position that a reader cannot work
        around?

        The five defects the placement exists to avoid, as a yes-or-no: the box
        hides another label, it leaves the axes, it lies across a line of the
        section, it has come to rest nearer some other labelled corner than the
        one it names, or its leader crosses another label's. A position with none
        of these is as good as a position gets, so the search stops at the first
        distance that offers one instead of buying a soft preference with the
        reader's eye-travel.
        """
        if any(bb.overlaps(pb) for pb in occupied):
            return True
        if offframe(tb):
            return True
        if any(_box_crosses_segment(bb, *s) for s in seg_px):
            return True
        if misread(point, tb):
            return True
        return any(_segments_cross(lead, other) for other in leaders)

    # --- place: nearest position that is unobjectionable, else cheapest ------
    placed = list(standing)
    layout = {}                    # point -> (dx, dy, ha, va, bb, tb,
                                   #           r, a, lead, w, h)
    leaders = []
    for point in order:
        x, y = point
        label = f"({x:g}, {y:g})"
        w, h, _ = renderer.get_text_width_height_descent(label, prop, False)
        px, py = pixels[point]

        best = None                # least bad position that is ON the axes
        anywhere = None            # least bad position at all
        hugs = None                # nearest position with nothing wrong with it
        for r in radii:
            here = None
            for a in angles:
                dx, dy = r * math.cos(a), r * math.sin(a)
                ha, va, bb, tb = box_at(px, py, w, h, dx, dy)
                lead = leader_of(point, bb)
                cost = (readability(point, bb, tb, r, a, placed)
                        + CROSSING * sum(_segments_cross(lead, other)
                                         for other in leaders))
                cand = (cost, dx, dy, ha, va, bb, tb, r, a, lead)
                if anywhere is None or cost < anywhere[0]:
                    anywhere = cand
                # A label off the frame is not on the figure at all, so no
                # position on it is ever traded for one off it — the frame is a
                # filter, not a term.
                if not offframe(tb):
                    if best is None or cost < best[0]:
                        best = cand
                    if not objectionable(point, bb, tb, lead, placed, leaders):
                        if here is None or cost < here[0]:
                            here = cand
            # This distance is clear somewhere: take the direction with the most
            # room and stop, rather than walking further for a softer preference.
            if here is not None:
                hugs = here
                break

        _cost, dx, dy, ha, va, bb, tb, r, a, lead = (
            hugs if hugs is not None else best if best is not None else anywhere)
        placed.append(bb)
        leaders.append(lead)
        layout[point] = (dx, dy, ha, va, bb, tb, r, a, lead, w, h)

    # --- and then untangle what placing them one at a time could not see -----
    #
    # Each label is placed against the ones already down, so a point reached
    # late can find every direction it wants taken by a leader an earlier point
    # laid across it. Exchanging the two offsets is the fix a draughtsman makes:
    # the pair that crossed is the pair whose labels are on each other's side,
    # and swapping them puts each label back on its own. The exchange is kept
    # only where it costs the figure nothing else, and the pass runs to a fixed
    # point over pairs taken in a fixed order, so the drawing is the same every
    # time it is made.
    def crossing_pairs():
        pts = list(order)
        return [(p, q) for i, p in enumerate(pts) for q in pts[i + 1:]
                if _segments_cross(layout[p][8], layout[q][8])]

    for _sweep in range(len(order) + 1):
        pairs = crossing_pairs()
        if not pairs:
            break
        improved = False
        for p, q in pairs:
            trial = {}
            for here, there in ((p, q), (q, p)):
                w, h = layout[here][9], layout[here][10]
                r2, a2 = layout[there][6], layout[there][7]
                dx, dy = r2 * math.cos(a2), r2 * math.sin(a2)
                px, py = pixels[here]
                ha, va, bb, tb = box_at(px, py, w, h, dx, dy)
                trial[here] = (dx, dy, ha, va, bb, tb, r2, a2,
                               leader_of(here, bb), w, h)
            # An exchange that walks a label off the axes is no exchange: the
            # crossing it undoes costs a reader less than a label that is not on
            # the figure.
            if any(offframe(trial[t][5]) and not offframe(layout[t][5])
                   for t in (p, q)):
                continue
            swapped = dict(layout)
            swapped.update(trial)
            others = [layout[o][4] for o in order if o not in (p, q)]
            pts = list(order)

            def total(state):
                out = 0.0
                for point, mate in ((p, q), (q, p)):
                    _dx, _dy, _ha, _va, bb, tb, r, a, _lead, _w, _h = state[point]
                    out += readability(point, bb, tb, r, a,
                                       standing + others + [state[mate][4]])
                out += CROSSING * sum(
                    _segments_cross(state[i][8], state[j][8])
                    for n, i in enumerate(pts) for j in pts[n + 1:])
                return out

            if total(swapped) < total(layout) - 1e-9:
                layout = swapped
                improved = True
                break
        if not improved:
            break

    # --- draw ---------------------------------------------------------------
    for point in order:
        x, y = point
        dx, dy, ha, va, _bb, tb, r, _a, _lead, _w, _h = layout[point]
        kwargs = {}
        # A label at the tight offset sits beside its point, a character's width
        # off it, and reads as its own. Every other label on the figure moved for
        # a reason a reader cannot see, so it carries a leader — as does one that
        # came to rest nearer some other vertex than the one it names. The figure
        # then says of every coordinate either that it is beside its corner or
        # which corner it belongs to, and never neither.
        tied = misread(point, tb)
        # The open run between the box and the vertex, in points, and the gap
        # taken off its far end. A label whose box is nearly against its point
        # has less run than the gap wants, so the gap is never allowed to eat
        # the line it is trimming: the leader shortens with it.
        run = _gap(tb, *pixels[point]) / scale - BOX_PAD - SHRINK_A
        gap = max(0.0, min(leader_gap, run / 2.0))
        # And where trimming would leave a stub shorter than the gap itself,
        # there is nothing worth drawing: the coordinate is already against its
        # corner, which is what a leader would have been drawn to say. A label
        # the reader could otherwise misattribute keeps its line regardless —
        # short, but the only thing on the figure saying which corner is whose —
        # as does every label when the caller asked for leaders throughout.
        if arrows or tied or (r > radii[0] + 1e-9 and run - gap >= gap):
            kwargs['arrowprops'] = dict(arrowstyle="-", color="0.45",
                                        linewidth=0.5, shrinkA=SHRINK_A,
                                        shrinkB=gap)
        ax.annotate(f"({x:g}, {y:g})", (x, y), textcoords="offset points",
                    xytext=(dx, dy), fontsize=fontsize, color="black",
                    ha=ha, va=va,
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

#: Slice-number labels: the size one is drawn at when there is room for it, the
#: smallest it is ever reduced to, and the rounded box's padding in matplotlib's
#: own units (a fraction of the font size, either side of the text).
#:
#: The floor is what still reads AS PRINTED, not as rendered. A report figure is
#: drawn at :data:`xslope.report.FIGURE_SIZE` and placed at the text width of a
#: portrait page, about seven tenths of it, so type in the figure reaches the
#: page seven tenths of the size it was set at — a floor chosen against the
#: rendered figure put three-point numbers on the page. The size a key opens on
#: is matplotlib's own default text size, which is what the rest of the figure's
#: type is set at: a slice number should read like the axis it sits over, not
#: smaller than it, and a key with room for that size should take it.
#:
#: The floor binds, and where it binds it is the number that is kept and the
#: clearance that is given up. A slice a sixth the width of the ones either side
#: — which is what a material boundary landing beside a slice boundary leaves —
#: has no room at its own mid-height for a number anybody could read, at any size
#: down to the floor, and a key printed at four points to separate that one pair
#: keys nothing at all.
SLICE_LABEL_PT = 10.0
SLICE_LABEL_MIN_PT = 6.0
SLICE_LABEL_PAD = 0.2


def slice_label_size(ax, slice_df, base=SLICE_LABEL_PT):
    """The size every slice number is set at.

    Every number stands at the middle of its own slice, and they are all one
    size, so the only thing left to decide is what that size is. What a number
    has to clear is its NEIGHBOURS, not the walls of its own slice. On a curved
    surface the middles climb: over a rising limb two labels a third of a label
    width apart across the page are a whole label height apart up it, and neither
    touches the other. Sizing them to fit between the slice boundaries ignored
    that and drove a forty-slice key to its floor for a collision that was not
    happening — three-point numbers on the printed page.

    So the collision is the measurement, taken on the axes as laid out. A label
    box is ``text width + 2·pad·size`` by ``text height + 2·pad·size`` and both
    are linear in the size, so a pair of labels clears at any size up to the
    larger of what their horizontal and their vertical separation allows, and the
    tightest pair on the figure governs. The whole key comes down to that size
    together: a key whose numbers are not all the same size, or not all at the
    same height in their slices, is read as saying something about the slices.

    Where there is no renderer to measure with — a figure that has never been
    drawn — ``base`` is returned and matplotlib does what it always did.
    """
    if slice_df is None or not len(slice_df):
        return base
    label = max((str(int(n)) for n in slice_df['slice #'].values), key=len)
    try:
        probe = ax.text(0, 0, label, fontsize=base, fontweight='bold')
        extent = probe.get_window_extent(ax.figure.canvas.get_renderer())
        probe.remove()
    except Exception:
        return base
    to_pt = 72.0 / ax.figure.dpi                          # display pixels -> points
    # The widest label's box at ``base``, in points. The widest is used for all
    # of them: one size is drawn, so the size has to suit the largest number.
    box_w = extent.width * to_pt + 2 * SLICE_LABEL_PAD * base
    box_h = extent.height * to_pt + 2 * SLICE_LABEL_PAD * base
    if box_w <= 0 or box_h <= 0:
        return base

    where = ax.transData.transform(np.column_stack([
        slice_df['x_c'].values.astype(float),
        slice_label_height(slice_df)]))
    dx = np.abs(where[:, None, 0] - where[None, :, 0]) * to_pt
    dy = np.abs(where[:, None, 1] - where[None, :, 1]) * to_pt
    # Each pair clears at any size up to whichever separation gives it more room.
    room = np.maximum(dx / box_w, dy / box_h) * base
    np.fill_diagonal(room, np.inf)
    return float(min(base, max(SLICE_LABEL_MIN_PT, float(room.min()))))


def slice_label_height(slice_df):
    """Where every slice number stands: the mid-height of its own slice, at the
    slice's own center line.

    Halfway between the base and the top at ``x_c``, and nothing else — no
    collision offset, no alternating step. The numbers then follow the mass in
    one band, and a reader looking for slice 20 knows where in slice 20 to look.
    """
    return (slice_df['y_cb'].values.astype(float)
            + slice_df['y_ct'].values.astype(float)) / 2.0


def plot_slice_numbers(ax, slice_df, fontsize=None):
    """
    Plots the slice number in the middle of each slice at the middle height.
    Numbers are 1-indexed.

    Parameters:
        ax: matplotlib Axes object
        slice_df: DataFrame containing slice data
        fontsize: label size in points; None sizes them all, together, so that no
            number is printed over its neighbour at that height
            (:func:`slice_label_size`), so a hundred-slice surface reads as well
            as a fifteen-slice one.

    Returns:
        None
    """
    if slice_df is not None:
        if fontsize is None:
            fontsize = slice_label_size(ax, slice_df)
        heights = slice_label_height(slice_df)
        for i, (_, row) in enumerate(slice_df.iterrows()):
            # Plot the slice number (1-indexed) at the middle of the slice, at
            # the middle of its height.
            slice_number = int(row['slice #'])
            ax.text(float(row['x_c']), float(heights[i]), str(slice_number),
                   ha='center', va='center', fontsize=fontsize, fontweight='bold',
                   bbox=dict(boxstyle=f"round,pad={SLICE_LABEL_PAD}",
                             facecolor='white', alpha=0.8),
                   zorder=6, gid='SLICE_NUMBER')


#: What the sliced mass DRAWS, by the gid each of those artists carries. A
#: slice-key figure is framed on these and on nothing else: the ground surface,
#: the material zones and the loads are context that routinely runs a long way
#: past the slices, and a frame that included them would not be a key to
#: anything.
SLICED_MASS_GIDS = ("SLICES", "FAILURE_SURFACE", "EFF_NORMAL_STRESS",
                    "PORE_PRESSURE")


def sliced_mass_bounds(ax, gids=SLICED_MASS_GIDS):
    """``(x0, x1, y0, y1)`` around everything the sliced mass has drawn, or None.

    Read off the artists themselves rather than recomputed from the slice table,
    so the base-stress bars — which stand off the base by a fraction of the slice
    height and are the reason a frame on the slice corners alone clips — are
    inside the box by construction.
    """
    want = set(gids)
    xs, ys = [], []
    for artist in ax.get_children():
        if artist.get_gid() not in want:
            continue
        if hasattr(artist, "get_xydata"):
            xy = np.asarray(artist.get_xydata(), dtype=float)
        elif hasattr(artist, "get_path"):
            xy = np.asarray(artist.get_path().vertices, dtype=float)
        else:
            continue
        if xy.size == 0:
            continue
        good = np.isfinite(xy).all(axis=1)
        if not good.any():
            continue
        xs.append(xy[good, 0])
        ys.append(xy[good, 1])
    if not xs:
        return None
    x = np.concatenate(xs)
    y = np.concatenate(ys)
    return float(x.min()), float(x.max()), float(y.min()), float(y.max())

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
            draw_water_level_symbol(ax, mid_x, mid_y, color=color, markersize=8, extra_gap_points=2.0)
    
    # Plot both piezometric lines
    f1 = feature_style(style, "piezo_line")
    f2 = feature_style(style, "piezo_line2")
    plot_single_piezo_line(ax, slope_data.get('piezo_line'), f1.get('color', 'b'),
                           "Piezometric Line", f1.get('linewidth', 2), f1.get('linestyle', '-'))
    plot_single_piezo_line(ax, slope_data.get('piezo_line2'), f2.get('color', 'skyblue'),
                           "Piezometric Line 2", f2.get('linewidth', 2), f2.get('linestyle', '-'))


# --- Seep boundary-condition color families -------------------------------------
# Each seep BC (set × type) is drawn in ONE hue: the boundary polyline in the DARK
# shade, its water-level line(s) + apex-down symbol in the LIGHTER shade. A boundary
# and a level that coincide (e.g. a tailwater head drawn along the ground at its own
# elevation) therefore still read as one feature, while the distinct hues keep the
# families apart when a rapid-drawdown pair (two BC sets) is shown together.
#
# Set-1 head defers to the style sheet (seep_bc / seep_water_level → navy / pale
# blue). The reservoir family is a DIFFERENT blue — both are water — separated from
# the navy head by hue (cyan-leaning azure vs indigo) and by lightness/saturation so
# the two never confuse. Set-2 head carries the non-blue distinction (rose), chosen
# to avoid the page's existing vocabulary — exit face (red / orangered), specified
# flux (green), dload (purple):
#   set-1 reservoir  → azure  (dark cerulean boundary / light azure levels)
#   set-2 head       → rose   (the constant-steady rapid-drawdown second set)
# Set 2 is the constant-steady rapid-drawdown set and NEVER carries a reservoir or a
# time-varying value (fileio rejects both at load time), so it has no reservoir hue.
_SEEP_BC_RESERVOIR = ("#0277bd", "#4fc3f7")   # set-1 reservoir  (dark azure / light azure)
_SEEP_BC_SET2_HEAD = ("#9c2a6e", "#e39ec8")   # set-2 head       (dark / light)


def _eval_bc_series_at(tseep, name, t):
    """Reservoir/head series value at time ``t``: linear between the series' own
    non-blank breakpoints, held constant beyond the ends (mirrors the solver's
    ``seep._eval_series``). Returns None when the series is unavailable. Shared by
    the inputs BC rendering and the solution BC-level overlay so both read a series
    exactly as the solver does."""
    if not tseep:
        return None
    times = tseep.get("times") or []
    vals = (tseep.get("series") or {}).get(name)
    if not times or vals is None:
        return None
    ta, va = [], []
    for tt, vv in zip(times, vals):
        if vv is not None:
            ta.append(float(tt)); va.append(float(vv))
    if not ta:
        return None
    if t <= ta[0]:
        return va[0]
    if t >= ta[-1]:
        return va[-1]
    j = int(np.searchsorted(ta, t, side="right")) - 1
    j = min(max(j, 0), len(ta) - 2)
    if ta[j + 1] == ta[j]:
        return va[j + 1]
    return va[j] + (va[j + 1] - va[j]) * (t - ta[j]) / (ta[j + 1] - ta[j])


def _reservoir_surface_x(coords, h):
    """(x_upstream, x_meets_face) for a horizontal water surface at level ``h``
    against a boundary polyline: the surface spans from the polyline's upstream
    (min-x) extent out to where the face first rises above ``h``, so it never draws
    into the embankment body."""
    pxs = [c[0] for c in coords]
    x_up = min(pxs)
    x_face = max(pxs)  # h at or above the face top: cover the whole extent
    for (x0, y0), (x1, y1) in zip(coords[:-1], coords[1:]):
        if y0 == y1 == h:
            x_face = max(x0, x1)                       # flat run exactly at h
        elif min(y0, y1) <= h <= max(y0, y1) and y1 != y0:
            x_face = x0 + (h - y0) / (y1 - y0) * (x1 - x0)
            break
    return x_up, x_face


def seep_bc_level_color(style, set_no, kind):
    """Light water-level shade for a seep BC's (set, kind) — the lighter half of the
    boundary/level color pair (see the family notes above). Set-1 head defers to the
    style sheet (``seep_water_level``); set-1 reservoir is light azure; set-2 head is
    light rose. Set 2 never carries a reservoir (fileio enforces it), so a (2,
    reservoir) request falls back to the set-2 head shade rather than inventing a
    hue that no valid file can produce."""
    from .style import resolve_style, feature_style
    kind = str(kind or "head").strip().lower()
    if int(set_no) == 1 and kind != "reservoir":
        return feature_style(resolve_style(style), "seep_water_level").get(
            "color", "lightskyblue")
    if int(set_no) == 1:
        return _SEEP_BC_RESERVOIR[1]
    return _SEEP_BC_SET2_HEAD[1]


def draw_water_level_symbol(ax, x, y, color, markersize=8, extra_gap_points=2.0):
    """Place an inverted-triangle water-surface marker so its TIP visually sits on the
    waterline at (x, y). Positioned in display (point) units via an offset transform,
    so it renders at a consistent size on any domain and its tip always touches the
    line. This is the ONE water-level symbol in the package — the inputs BC rendering,
    the piezometric-line marker, and the solution-frame water-level overlay all call
    it, so the symbol is identical everywhere."""
    from matplotlib.markers import MarkerStyle
    from matplotlib.transforms import offset_copy

    # Matplotlib scales marker vertices by `markersize` (points) for a Line2D marker;
    # the "v" tip is the lowest vertex, so offset the marker center up by that much
    # (plus a small gap) to land the tip exactly on (x, y).
    ms = MarkerStyle("v")
    path = ms.get_path().transformed(ms.get_transform())
    verts = np.asarray(path.vertices)
    min_y = float(verts[:, 1].min())
    tip_offset_points = (-min_y) * float(markersize) + float(extra_gap_points)
    trans = offset_copy(ax.transData, fig=ax.figure, x=0.0, y=tip_offset_points,
                        units="points")
    ax.plot([x], [y], marker="v", color=color, markersize=markersize,
            linestyle="None", transform=trans)


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

    def _plot_one_bc_set(ax, seepage_bc, geom_width, x_min_geom, x_max_geom,
                         head_line_color, water_level_color, exit_face_color, label_suffix="",
                         head_lw=3, head_ls="--", water_lw=2, exit_lw=3, exit_ls="--",
                         flux_color="darkgreen", flux_lw=3, flux_ls="-.",
                         reservoir_color=_SEEP_BC_RESERVOIR[0], reservoir_water_color=_SEEP_BC_RESERVOIR[1],
                         tseep=None):
        """Plot a single set of seepage boundary conditions.

        A head block carries a ``kind`` ("head" or "reservoir") and a value that is
        either a constant (a number) or the name of a tseep time series (v18). A
        reservoir (submerged-only) boundary is drawn distinctly from a fixed head,
        and a series-bound boundary draws TWO waterlines -- the t = 0 level and the
        t = end level -- each with the standard apex-down water symbol. A constant
        head renders exactly as before.
        """
        specified_heads = seepage_bc.get("specified_heads") or []
        specified_fluxes = seepage_bc.get("specified_fluxes") or []
        exit_face = seepage_bc.get("exit_face") or []

        # Each distinct legend label is emitted once (robust to block ordering:
        # a reservoir block ahead of a head block must not suppress the head label).
        _labeled = set()
        def _once(text):
            if text in _labeled:
                return ""
            _labeled.add(text)
            return text

        # t = end of the run (series held constant beyond it) for the second waterline
        t_end = None
        if tseep:
            t_end = tseep.get("duration")
            if t_end is None:
                _tt = tseep.get("times") or []
                t_end = max(_tt) if _tt else 0.0

        for i, sh in enumerate(specified_heads):
            coords = sh.get("coords") or []
            if len(coords) < 2:
                continue

            xs, ys = zip(*coords)
            is_reservoir = str(sh.get("kind", "head")).lower() == "reservoir"
            line_color = reservoir_color if is_reservoir else head_line_color
            line_label = ("Reservoir Boundary" if is_reservoir
                          else "Specified Head Line") + label_suffix
            ax.plot(
                xs, ys,
                color=line_color, linewidth=head_lw, linestyle=head_ls,
                label=_once(line_label),
            )

            head_val = sh.get("head", None)
            if head_val is None:
                continue

            # v18: a string value names a tseep series -> a time-varying boundary.
            # Draw its level at t = 0 and at t = end, each a horizontal reservoir
            # surface with the standard apex-down water symbol (tip on the line).
            if isinstance(head_val, str):
                h0 = _eval_bc_series_at(tseep, head_val, 0.0)
                h1 = _eval_bc_series_at(tseep, head_val, float(t_end)) if t_end is not None else None
                wl_color = reservoir_water_color if is_reservoir else water_level_color
                wl_label = ("Reservoir Level (t = 0, t = end)" if is_reservoir
                            else "Head Level (t = 0, t = end)") + label_suffix
                for hlev, tlab in ((h0, "t = 0"), (h1, "t = end")):
                    if hlev is None:
                        continue
                    x_up, x_face = _reservoir_surface_x(coords, float(hlev))
                    mid_x = 0.5 * (x_up + x_face)
                    ax.plot([x_up, x_face], [hlev, hlev], color=wl_color,
                            linewidth=water_lw, linestyle="-", label=_once(wl_label))
                    draw_water_level_symbol(ax, mid_x, float(hlev),
                                            color=wl_color, markersize=8,
                                            extra_gap_points=2.0)
                    # label above the water symbol so it stays clear even when the
                    # drawn-down surface is a short segment against a steep face
                    ax.annotate(tlab, xy=(mid_x, hlev), xytext=(0, 11),
                                textcoords="offset points", ha="center", va="bottom",
                                color=wl_color, fontsize=8)
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

            # A constant reservoir (numeric value on a submerged-only face) still
            # belongs to the azure reservoir family: draw its level in reservoir_water_color so a
            # boundary/level pair stays one hue. The common plain-head case is
            # unchanged (water_level_color, "Specified Head Water Level").
            const_wl_color = reservoir_water_color if is_reservoir else water_level_color
            const_wl_label = (f"Reservoir Water Level{label_suffix}" if is_reservoir
                              else f"Specified Head Water Level{label_suffix}")
            ax.plot(
                wl_xs, wl_ys,
                color=const_wl_color, linewidth=water_lw, linestyle="-",
                label=_once(const_wl_label),
            )

            if len(wl_xs) > 1:
                try:
                    pairs = sorted(zip(wl_xs, wl_ys), key=lambda p: p[0])
                    sx, sy = zip(*pairs)
                    mid_x = 0.5 * (min(sx) + max(sx))
                    mid_y = float(np.interp(mid_x, sx, sy))
                    draw_water_level_symbol(ax, mid_x, mid_y, color=const_wl_color, markersize=8, extra_gap_points=2.0)
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
    tseep = slope_data.get("tseep")
    _plot_one_bc_set(ax, seepage_bc, geom_width, x_min_geom, x_max_geom,
                     head_line_color=fh.get("color", "darkblue"),
                     water_level_color=fw.get("color", "lightskyblue"),
                     exit_face_color=fe.get("color", "red"), label_suffix=label_suffix,
                     head_lw=fh.get("linewidth", 3), head_ls=fh.get("linestyle", "--"),
                     water_lw=fw.get("linewidth", 2),
                     exit_lw=fe.get("linewidth", 3), exit_ls=fe.get("linestyle", "--"),
                     flux_color=ff.get("color", "darkgreen"),
                     flux_lw=ff.get("linewidth", 3), flux_ls=ff.get("linestyle", "-."),
                     tseep=tseep)

    # Plot second set of BCs if present. Set 2 is the constant-steady rapid-drawdown
    # set: it never carries a reservoir type or a time-varying (tseep series) value
    # (fileio rejects both at load time), so it draws in its own single hue family —
    # a rose head pair (dark boundary / light level) distinct from set 1's navy head
    # and azure reservoir. It therefore needs no reservoir-color override and no tseep.
    if has_bc2:
        seepage_bc2 = slope_data.get("seepage_bc2") or {}
        _plot_one_bc_set(ax, seepage_bc2, geom_width, x_min_geom, x_max_geom,
                         head_line_color=_SEEP_BC_SET2_HEAD[0],
                         water_level_color=_SEEP_BC_SET2_HEAD[1],
                         exit_face_color="orangered", label_suffix=" (BC 2)",
                         flux_color="seagreen")


def plot_derived_water_lines(ax, slope_data, style=None):
    """Plot the water surface a model's seepage head/reservoir boundaries imply.

    A head or reservoir block traced along the ground surface IS a statement of
    where a pool stands, and the engine turns exactly that statement into the
    ponded-water distributed load it applies. A plot that draws the loads must
    therefore draw their source, or the reader sees a water load with no water.

    The line comes from :func:`xslope.water.water_line_for_stage` — the one
    derivation the solver's automatic water loads, the vendor importers and the
    preflight remedy all read — so the drawn pool and the applied load can never
    be different pools. Only the runs that stand ABOVE the ground are drawn: the
    derived line lies on the ground everywhere else, and tracing that would draw a
    second ground surface. Each run is carried one sample past its ends so it
    meets the ground at the shoreline, and marked with the package's one
    water-level symbol.

    Piezometric lines are :func:`plot_piezo_line`'s; this layer covers only the
    boundary-condition source.

    Returns the runs it drew, as ``[[(x, y), ...], ...]`` (empty when the model
    states no pool).
    """
    from .water import water_line_for_stage, _y_on

    ground = slope_data.get('ground_surface')
    if ground is None or ground.is_empty:
        return []
    xs_g = [c[0] for c in ground.coords]
    tol = max(1e-12, 1e-6 * (max(xs_g) - min(xs_g)))

    drawn = []
    for stage, key in ((1, 'seepage_bc'), (2, 'seepage_bc2')):
        bc = slope_data.get(key) or {}
        heads = bc.get('specified_heads') or []
        if not heads:
            continue
        pts = water_line_for_stage(slope_data, stage=stage).get('points') or []
        if len(pts) < 2:
            continue

        wet = [(_y_on(ground, float(x)) is not None
                and y - _y_on(ground, float(x)) > tol) for x, y in pts]
        runs, i = [], 0
        while i < len(pts):
            if not wet[i]:
                i += 1
                continue
            j = i
            while j + 1 < len(pts) and wet[j + 1]:
                j += 1
            runs.append(pts[max(0, i - 1):min(len(pts), j + 2)])
            i = j + 1
        if not runs:
            continue

        kinds = {str(b.get('kind') or 'head').strip().lower() for b in heads}
        color = seep_bc_level_color(style, stage,
                                    'reservoir' if 'reservoir' in kinds else 'head')
        label = "Water Surface" if stage == 1 else "Water Surface 2"
        for run in runs:
            rx, ry = zip(*run)
            ax.plot(rx, ry, color=color, linewidth=2, linestyle='-',
                    label=label, gid='WATER_LINE')
            label = None                       # one legend key per stage
            mid = len(run) // 2
            draw_water_level_symbol(ax, run[mid][0], run[mid][1], color=color,
                                    markersize=8, extra_gap_points=2.0)
            drawn.append(list(run))
    return drawn


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

    Four sets, not two: the user's own blocks from the two dloads sheets (purple /
    orange) and, in automatic water-load mode, the ponded-water loads the engine
    derived for each stage (water-blue). The derived sets are drawn from the same
    derivation the solver uses, so what a plot shows and what an analysis carries
    are the same load — automatic must never mean invisible, and a load the user
    did not type is exactly the kind that has to be MORE visible, not less.
    """
    from .style import resolve_style, feature_style
    from .water import with_water_loads, DERIVED_KEYS
    style = resolve_style(style)
    slope_data = with_water_loads(slope_data)
    gamma_w = slope_data['gamma_water']
    ground_surface = slope_data['ground_surface']

    def plot_single_dload_set(ax, dloads, color, label, linewidth=1.5, dirs=None):
        """Internal function to plot a single set of distributed loads.

        `dirs` is the parallel per-block direction list ('normal' | 'vertical', v21).
        A vertical block is DRAWN vertical — the arrows stand straight up off the
        loaded surface instead of leaning with it — because the direction is the whole
        point of the option and an input plot that draws both the same way would hide
        the one thing the user set."""
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
        
        for _li, line in enumerate(dloads):
            if len(line) < 2:
                continue

            vertical = str((dirs or [])[_li] if dirs and _li < len(dirs)
                           else 'normal').lower() == 'vertical'
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
                
                # Perpendicular direction (rotate 90 degrees CCW), or straight up
                # for a vertical (gravity-surcharge) block.
                if vertical:
                    perp_dx, perp_dy = 0.0, 1.0
                else:
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
                          df1.get('linewidth', 1.5),
                          dirs=slope_data.get('dload_dirs'))
    plot_single_dload_set(ax, dloads2, df2.get('color', 'orange'), 'Distributed Load 2',
                          df2.get('linewidth', 1.5),
                          dirs=slope_data.get('dload2_dirs'))
    dfd1 = feature_style(style, "dloads_derived")
    dfd2 = feature_style(style, "dloads2_derived")
    plot_single_dload_set(ax, slope_data.get(DERIVED_KEYS[1]) or [],
                          dfd1.get('color', 'royalblue'), 'Water Load (derived)',
                          dfd1.get('linewidth', 1.5))
    plot_single_dload_set(ax, slope_data.get(DERIVED_KEYS[2]) or [],
                          dfd2.get('color', 'mediumturquoise'),
                          'Water Load 2 (derived)', dfd2.get('linewidth', 1.5))

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

def _table_value(value):
    """A material property as a number for a comparison, with an unset one as 0.

    A blank mat-sheet cell edited in Studio now reaches the model as ``None``
    rather than as an invented ``0.0``, so anything that compares a property has to
    survive one. Absent reads as zero here, which is what the comparisons meant all
    along -- ``mat.get('d', 0) > 0`` was already asking "is a dilation entered?".
    """
    try:
        f = float(value)
    except (TypeError, ValueError):
        return 0.0
    return 0.0 if f != f else f


def _table_cell(value, fmt="{:.1f}"):
    """A material property formatted for a table cell, or ``-`` when it is unset.

    The counterpart of :func:`_table_value` for display: an absent property prints
    as a dash rather than as ``0.0``, because a table that shows a unit weight of
    0.0 for a material the user has not filled in is stating something false.
    """
    try:
        f = float(value)
    except (TypeError, ValueError):
        return "-"
    return "-" if f != f else fmt.format(f)


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
    has_d_psi = any(_table_value(mat.get('d')) > 0 or _table_value(mat.get('psi')) > 0
                    for mat in materials)

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
        gamma_str = _table_cell(mat.get('gamma'))
        option = mat['option']
        d = _table_value(mat.get('d'))
        psi = _table_value(mat.get('psi'))
        d_str = f"{d:.1f}" if d > 0 or psi > 0 else "-"
        psi_str = f"{psi:.1f}" if d > 0 or psi > 0 else "-"
        tail = [d_str, psi_str] if has_d_psi else []

        if option == 'mc':
            cells = [_table_cell(mat.get('c')), _table_cell(mat.get('phi'))]
        elif option == 'cp':
            cells = [_table_cell(mat.get('cp'), "{:.2f}"),
                     _table_cell(mat.get('r_elev'))]
        else:
            cells = ["-", "-"]
        table_data.append([idx + 1, name, gamma_str] + cells + tail)

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

    # 4) account for distributed loads extending above ground surface — the derived
    # water loads among them, or a tall reservoir nobody typed would be clipped out
    # of the very view that exists to show it.
    from .water import with_water_loads as _with_water
    data = _with_water(data)
    _dload_sets = [data.get('dloads', []), data.get('dloads2', []),
                   data.get('dloads_derived', []), data.get('dloads2_derived', [])]
    if any(_dload_sets):
        # Only a model that actually has distributed loads needs gamma_w here; read
        # it loudly (no silent default that could flip the unit system) rather than
        # forcing every load-free plot to carry a gamma_water.
        gamma_w = require_gamma_water(data, "distributed-load plot extents")
        for dloads in _dload_sets:
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


def adaptive_colorbar_ticks(fig, cbar, steps=(2, 5, 10), min_ticks=2,
                            label_heights=1.8):
    """Thin a vertical colorbar's ticks to its DRAWN height so labels never stack.

    Allow at most one tick per ``label_heights`` label-heights of the colorbar's
    rendered height (a legibility margin; the aspect enters only through the measured
    height, not any per-figure constant), with ``steps=(2, 5, 10)`` restricting the
    tick VALUES to nice integers/round decimals. A short/wide bar thins to a few
    readable ticks; a full-height bar keeps finer labeling. Shared by the FEM result
    plots and the seepage flow nets (originally the plot_seep sheetpile fix,
    a35e308), so both get uncrowded, round-valued colorbars from one rule.

    Measures the colorbar axes' window extent after a draw, so it works for a
    make_axes_locatable cax, a constrained-layout colorbar, and a manually placed
    fractional-height cax alike.
    """
    from matplotlib.ticker import MaxNLocator, ScalarFormatter
    from matplotlib.font_manager import FontProperties
    try:
        fig.canvas.draw()
        h_px = cbar.ax.get_window_extent().height
    except Exception:
        return                                  # no renderer: leave default ticks
    h_pts = h_px * 72.0 / float(fig.dpi)
    label_pts = FontProperties(size=plt.rcParams["ytick.labelsize"]).get_size_in_points()
    max_ticks = max(min_ticks, int(h_pts // (label_heights * label_pts)))
    cbar.locator = MaxNLocator(nbins=max_ticks, steps=list(steps))
    # Short labels for near-zero ranges: factor a shared power out of the tick
    # labels (a "×10⁻⁴" header) instead of stacking seven-decimal strings, which
    # smear at panel size — the tiny viscoplastic strain fields (e.g. RS2-32b,
    # max ~1e-4) are the motivating case. Applied via cbar.formatter (update_ticks
    # re-applies it; setting the axis formatter directly would be clobbered). A
    # no-op for normal ranges, so it is safe to apply unconditionally.
    fmt = ScalarFormatter(useMathText=True)
    fmt.set_powerlimits((-3, 4))                # offset only outside 1e-3..1e4
    cbar.formatter = fmt
    cbar.update_ticks()


# Mesh element-edge linewidth: calibrated so a mesh whose elements render at
# EDGE_REF_PTS points across draws at BASE_LW (today's 0.5 pt weight — a sparse
# teaching mesh), scaling linearly with rendered edge size and clamped to a
# hairline floor / a modest ceiling. Fine meshes (small rendered edges) drop to
# the floor so they read as a mesh instead of fusing into a solid fill; coarse
# meshes stay at today's weight. Shared by the plot.py mesh overlays (plot_mesh
# and the plot_inputs mesh background); plot_fem.py applies the identical rule.
_EDGE_LW_BASE = 0.5       # weight at the reference (teaching-mesh) edge size
_EDGE_LW_REF_PTS = 28.0   # rendered edge length (points) that maps to _EDGE_LW_BASE
_EDGE_LW_MIN = 0.25       # hairline floor (dense meshes)
_EDGE_LW_MAX = 1.0        # ceiling (very coarse meshes)


def adaptive_edge_linewidth(ax, fig, segments,
                            base_lw=_EDGE_LW_BASE, ref_pts=_EDGE_LW_REF_PTS,
                            lw_min=_EDGE_LW_MIN, lw_max=_EDGE_LW_MAX, sample=4000):
    """Line width (points) for a mesh's element edges, scaled to their RENDERED
    size so a dense mesh thins to legible hairlines instead of fusing into a solid
    black fill (the "nearly solid" disease at high element counts).

    ``segments`` is an iterable of edges, each a 2-point ``[(x0, y0), (x1, y1)]`` in
    DATA coordinates. The MEDIAN edge length is measured in display points (via the
    axes' data->display transform, after a draw so the transform is final) and
    mapped linearly ``base_lw * median_pts / ref_pts``, clamped to ``[lw_min,
    lw_max]``. Calibrated (ref_pts) so a sparse teaching mesh keeps ~``base_lw``
    (today's 0.5 pt) while a fine mesh drops toward the floor. Returns ``base_lw``
    if there is no renderer or no measurable edge, so behaviour degrades safely."""
    import numpy as _np
    try:
        segs = list(segments)
        if not segs:
            return base_lw
        if len(segs) > sample:                 # subsample: the median is stable
            idx = _np.linspace(0, len(segs) - 1, sample).astype(int)
            segs = [segs[i] for i in idx]
        pts = _np.asarray([[s[0][0], s[0][1], s[1][0], s[1][1]] for s in segs],
                          dtype=float)
        fig.canvas.draw()
        tr = ax.transData
        p0 = tr.transform(pts[:, 0:2])
        p1 = tr.transform(pts[:, 2:4])
        px = _np.hypot(p1[:, 0] - p0[:, 0], p1[:, 1] - p0[:, 1])
        med_px = float(_np.median(px[px > 0])) if _np.any(px > 0) else 0.0
        if med_px <= 0:
            return base_lw
        med_pts = med_px * 72.0 / float(fig.dpi)
        return float(min(lw_max, max(lw_min, base_lw * med_pts / ref_pts)))
    except Exception:
        return base_lw


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

    def _xdecor_band_frac():
        """Figure-fraction height of the x-axis decoration band that hangs BELOW
        the axes spine — the tick labels and, if present, the x-axis label. The
        legend must clear THIS band, not just the spine box, or its top row lands
        on the tick numbers (the crowding Norm flagged on the FEM Inputs panel).
        Measured in pixels, which are font-fixed and so invariant to wherever
        subplots_adjust later puts the box."""
        try:
            r = fig.canvas.get_renderer()
            spine_y0 = ax.get_window_extent().y0
            lo = spine_y0
            for t in ax.get_xticklabels():
                if not t.get_text():
                    continue
                bb = t.get_window_extent(r)
                if bb.y0 < spine_y0:              # only labels hanging below the spine
                    lo = min(lo, bb.y0)
            xlab = ax.xaxis.get_label()
            if xlab.get_text():
                bb = xlab.get_window_extent(r)
                if bb.y0 < spine_y0:
                    lo = min(lo, bb.y0)
            return max(0.0, (spine_y0 - lo) / fig.bbox.height)
        except Exception:
            return 0.0

    # The band between the spine bottom and the lowest tick label / axis label.
    # It is what the legend has to be pushed clear of.
    band = _xdecor_band_frac()
    # Reserve a bottom margin sized for legend + the x-decoration band that sits
    # above it (as if the axes filled down to it), plus a title-aware top.
    bottom = min(0.55, _LEGEND_BOTTOM + leg_h + gap + band)
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
        band = _xdecor_band_frac()           # re-measure after the final layout
    except Exception:
        y0 = bottom
    # Pin the legend a fixed gap below the x-decoration band (tick labels + axis
    # label), which itself hangs `band` below the spine bottom (y0) — so the
    # legend clears the tick numbers instead of landing on them.
    anchor_y = max(_LEGEND_BOTTOM + leg_h, y0 - band - gap)
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
    frame="fill",
    pad_frac=0.035,
    fig=None,
    style=None,
    show_mesh=None,
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
        mode: Which engine's view of the model to draw, and which material
            properties table to display with it:
            - "lem": Limit equilibrium materials (γ, c, φ, optional d/ψ)
            - "seep": Seepage properties (k₁, k₂, Angle, kr₀, h₀)
            - "fem": FEM properties (γ, c, φ, E, ν)
            - "shared": the SECTION every engine shares — its geometry and its
              material zones, and nothing else. Everything an analysis applies
              to the section or reads off it is suppressed: no trial circles or
              non-circular surfaces (they belong with the search that produced
              them), no background mesh (it belongs with the mesh figure), no
              reinforcement lines or piles, no loads — neither the distributed
              loads nor the line loads, nor the ponded-water loads derived from
              a water surface — and no water lines, neither a piezometric line
              nor the pool a head boundary states. Each of those is read by a
              particular analysis: a piezometric line sets pore pressure, a
              pool becomes a load, a bar carries tension. Each stability
              engine's own model figure draws the ones it reads; drawn here too
              they were the same lines twice, one page apart, on a model whose
              only feature is its water or its members. This is the plot the
              Analysis Report's Project Definition uses. Its material table, if
              asked for, is the LEM one.
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
        frame: How the panel is framed around its content:
            - "fill" (default): equal aspect with adjustable="datalim" — the data
              limits are padded out to fill the figure aspect. Best for the
              interactive studio canvas, which sizes the figure to the viewport
              and wants the geometry to fill it.
            - "content": equal aspect with adjustable="box" — the axes box shrinks
              to the data's TRUE proportions and a single uniform cushion
              (``pad_frac`` × the larger domain dimension, in data units, on both
              axes) is placed around the content, so the geometry is never drawn
              steeper than reality, the visual margins are equal, and nothing
              touches the frame. For a wide-thin domain the excess figure height
              becomes an outer margin that the caller's ``bbox_inches="tight"``
              crops — so "content" is meant for figure files saved tight-bbox
              (the corpus figure generators), not for the fill-the-viewport canvas.
        pad_frac: Cushion for frame="content", as a fraction of the larger domain
            dimension (default 0.035 ≈ 3.5%). Ignored for frame="fill".
        fig: Optional existing Matplotlib Figure to draw into (used for embedding in a
            GUI canvas). When None (default) a new pyplot figure is created and shown;
            when provided, the figure is cleared and reused and plt.show() is skipped.
        show_mesh: Whether an engine mode draws the analysis mesh behind the
            model. None (default) leaves it to the mode — every engine view
            draws it, the shared model does not. False suppresses it, which is
            what a caller that gives the mesh a figure of its own wants: drawn
            twice, the underlay is a grid over the zones it is there to show.

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

    # Plot mesh in background if available. The shared-model plot leaves it out:
    # the mesh is an engine's discretisation of the model, not the model. A
    # caller that draws the mesh a second time, in a figure of its own, says so
    # with ``show_mesh=False`` and gets the model without it.
    _mesh_bg_lc = None
    _mesh_bg_segments = None
    mesh = (slope_data.get('mesh')
            if mode != "shared" and show_mesh is not False else None)
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
            # Finalize the edge width after layout (below): a dense background mesh
            # must thin to a hairline, not fuse into a solid fill under the geometry.
            _mesh_bg_lc = lc
            _mesh_bg_segments = lines

    # Plot geometry: profile lines if provided (drawn as before), otherwise the
    # material-zone polygons.
    plot_base_geometry(ax, slope_data, labels=True, style=style)
    # SSR zone overlays are an FEM concept (they constrain the strength reduction),
    # so they are drawn only on the FEM input plot — where a reader is looking at
    # the model the SSRM will run.
    if mode == "fem":
        plot_ssr_zones(ax, slope_data, style=style)
    # A defined piezometric line is the model's water table whether or not a
    # material's u option consumes it: with u = "none" it still splits each
    # slice's weight into saturated below / moist above, so it changes the
    # answer and must be visible. It is a stability input — it sets the pore
    # pressure on a slice base and at a node — so it is drawn where the engine
    # that reads it is documented, and not on the shared section. Seep mode
    # draws its own water surfaces.
    if mode not in ("seep", "shared"):
        plot_piezo_line(ax, slope_data, style=style)
        # The water surface the head/reservoir boundaries state, drawn beside
        # the ponded-water load derived from it: a load with no water above it
        # is a load with no source.
        plot_derived_water_lines(ax, slope_data, style=style)
    if mode == "seep":
        plot_seepage_bc_lines(ax, slope_data, style=style)
    plot_tcrack_surface(ax, slope_data, style=style)
    # The members and the loads belong to the engine that carries them, not to
    # the section: each engine's model figure draws its own, and the shared
    # model leaves them out (see the "shared" entry under ``mode``). A seepage
    # analysis applies no load at all, so its view draws none either.
    if mode not in ("shared", "seep"):
        plot_dloads(ax, slope_data, style=style)
        plot_line_loads(ax, slope_data, style=style)
    if mode != "shared":
        plot_reinforcement_lines(ax, slope_data, style=style)
        plot_piles(ax, slope_data, style=style)

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

            if mode in ("lem", "shared"):
                has_d_psi = any(_table_value(mat.get('d')) > 0
                                or _table_value(mat.get('psi')) > 0
                                for mat in materials)
                width = 0.25 if has_d_psi else 0.2
                height = min(0.35, 0.06 + 0.035 * num_rows)
            elif mode == "fem":
                width = 0.60
                height = min(0.32, 0.06 + 0.035 * num_rows)
            elif mode == "seep":
                width = 0.45
                height = min(0.50, 0.10 + 0.06 * num_rows)
            else:
                raise ValueError(f"Unknown mode '{mode}'. Expected one of: "
                                 f"'lem', 'seep', 'fem', 'shared'.")

            return width, height

        def _plot_table(ax, xloc, yloc):
            """Plot the appropriate material table based on mode."""
            if mode in ("lem", "shared"):
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

    if frame == "content":
        # Frame the panel to its CONTENT: box-adjust equal aspect (the axes box
        # shrinks to the data's true proportions instead of padding the data out
        # to fill an arbitrary figure aspect — no interior dead-space slab), with a
        # single uniform cushion in DATA units on BOTH axes so the visual margins
        # are equal under equal aspect and the geometry never touches the frame.
        bb = ax.dataLim
        cx0, cx1 = bb.intervalx
        cy0, cy1 = bb.intervaly
        if cx1 > cx0 and cy1 > cy0:
            pad = pad_frac * max(cx1 - cx0, cy1 - cy0)
            ax.set_xlim(cx0 - pad, cx1 + pad)
            ax.set_ylim(cy0 - pad, cy1 + pad)
        ax.set_aspect('equal', adjustable='box')
    else:
        ax.set_aspect('equal', adjustable='datalim')  # ✅ Equal aspect
        # Add a bit of headroom so plotted lines/markers don't touch the top border
        y0, y1 = ax.get_ylim()
        if y1 > y0:
            pad = 0.05 * (y1 - y0)
            ax.set_ylim(y0, y1 + pad)
    ax.grid(False)

    # Get legend handles and labels
    handles, labels = ax.get_legend_handles_labels()

    # Add distributed load to legend if present
    if slope_data['dloads']:
        handler_class, dummy_line = get_dload_legend_handler()
        handles.append(dummy_line)
        labels.append('Distributed Load')
    # And the derived water loads, which are their own entry. A model in automatic
    # mode may carry no user dloads at all and still have a reservoir drawn on it,
    # so a legend keyed only on slope_data['dloads'] would leave the one load the
    # user cannot see in the sheets out of the plot that exists to show it.
    if _derived_water_blocks(slope_data):
        handles.append(get_dload_legend_handler(color=_derived_water_color(style))[1])
        labels.append('Water Load (derived)')

    if show_title:
        ax.set_title(title)
    # Axis length units, only when the model declares a unit system (units plan
    # phase 4). Undeclared models get no axis label — pixel-identical to today.
    _ul = declared_unit_labels(slope_data)
    if _ul and _ul.get("length"):
        ax.set_xlabel(f"x ({_ul['length']})")
        ax.set_ylabel(f"y ({_ul['length']})")
    fig.tight_layout()
    # Single shared legend recipe: below the axes, auto-columns, frameless by
    # default (frame toggled per-view via legend_frame). Reserves bottom margin.
    _legend_below(ax, fig, handles=handles, labels=labels,
                  legend_ncol=legend_ncol, frameon=legend_frame, show_legend=show_legend)

    # Density-aware background-mesh edge width, now that the layout (and so the
    # data->display scale) is final.
    if _mesh_bg_lc is not None and _mesh_bg_segments:
        _mesh_bg_lc.set_linewidth(adaptive_edge_linewidth(ax, fig, _mesh_bg_segments))

    # Coordinate labels go on last, after tight_layout and the reserved legend
    # margin have settled the axes box. Their placement is solved in display
    # pixels, so laying them out before the layout is final solves the puzzle at
    # one scale and prints it at another: the points move together as the axes
    # shrink, and a clearance measured on the wide box closes on the narrow one.
    if label_coordinates:
        plot_coordinate_labels(ax, slope_data, fontsize=coord_label_size,
                               arrows=coord_arrows, style=style)

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

def plot_solution(slope_data, slice_df, failure_surface, results, figsize=(12, 7), slice_numbers=False, seep_contours=True, save_png=False, save_dxf=False, dpi=300, legend_ncol="auto", legend_frame=False, show_title=True, show_legend=True, fig=None, style=None, frame="fill", pad_frac=0.035):
    """
    Plots the full solution including slices, numbers, thrust line, and base stresses.

    Parameters:
        data: Input data
        slice_df: DataFrame containing slice data
        failure_surface: Failure surface geometry
        results: Solution results
        figsize: Tuple of (width, height) in inches for the plot
        slice_numbers: Label each slice with its 1-indexed number, sized to fit
            inside the slice it names.
        frame: "fill" (default) frames the whole model, as the results view does.
            "slices" frames the SLICED MASS — the slices, the failure surface and
            the base-stress bars — with a uniform cushion, which is what makes a
            slice-key figure readable beside its slice table.
        pad_frac: Cushion for frame="slices", as a fraction of the larger of the
            sliced mass's two dimensions. Ignored for frame="fill".
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
    # Drawn whenever the model defines one — a piezometric line read only for
    # saturated unit weights (u = "none") still moves the factor of safety.
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
        gamma_w = require_gamma_water(slope_data, "seep head contours")
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
    # Slice numbers go on last of all, once the frame and the layout are final —
    # they are sized against the slice widths as they will actually print.
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
    # And the derived water loads, which are their own entry. A model in automatic
    # mode may carry no user dloads at all and still have a reservoir drawn on it,
    # so a legend keyed only on slope_data['dloads'] would leave the one load the
    # user cannot see in the sheets out of the plot that exists to show it.
    if _derived_water_blocks(slope_data):
        handles.append(get_dload_legend_handler(color=_derived_water_color(style))[1])
        labels.append('Water Load (derived)')
    
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

    if frame == "slices":
        # Frame the sliced mass: box-adjust equal aspect (the axes box takes the
        # data's true proportions instead of padding the mass out to fill the
        # figure) with one uniform cushion in DATA units on both axes, so the
        # margins read equal and nothing touches the frame.
        bounds = sliced_mass_bounds(ax)
        if bounds is not None:
            bx0, bx1, by0, by1 = bounds
            if bx1 > bx0 and by1 > by0:
                pad = pad_frac * max(bx1 - bx0, by1 - by0)
                ax.set_xlim(bx0 - pad, bx1 + pad)
                ax.set_ylim(by0 - pad, by1 + pad)
                ax.set_aspect('equal', adjustable='box')

    # Axis length units, only when the model declares a unit system (units plan
    # phase 4); undeclared models get no axis label — pixel-identical to today.
    _ul = declared_unit_labels(slope_data)
    if _ul and _ul.get("length"):
        ax.set_xlabel(f"x ({_ul['length']})")
        ax.set_ylabel(f"y ({_ul['length']})")

    fig.tight_layout()
    _legend_below(ax, fig, handles=handles, labels=labels,
                  legend_ncol=legend_ncol, frameon=legend_frame, show_legend=show_legend)

    # The labels are measured against the axes as laid out, so they are the last
    # thing on the plot: the frame, the aspect and the legend's reserved margin
    # are all settled by here and a slice is exactly as wide as it will print.
    if slice_numbers:
        plot_slice_numbers(ax, slice_df)

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


def plot_reliability_histogram(result, figsize=(9, 5.5), bins='auto', show_fits=True,
                               save_png=False, dpi=300, show_title=True,
                               show_legend=True, fig=None, style=None):
    """Histogram of the Monte Carlo factor-of-safety samples from
    :func:`xslope.reliability.reliability_mc`.

    Draws the distribution of FS over the admissible realizations, marks the
    FS = 1 failure line and the sample mean, shades the FS < 1 (failure) tail, and
    — when ``show_fits`` and the spread is well conditioned — overlays the fitted
    normal and lognormal PDFs (method of moments on the sample mean/σ). The title
    carries the mean FS, σ_F, both reliability-index conventions (β normal and β
    lognormal) and the empirical probability of failure.

    ``result`` is the dict returned by ``reliability_mc`` (it reads ``fs_samples``
    and the summary statistics), or a bare array/sequence of FS samples.
    """
    own_fig = fig is None
    if own_fig:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig.clear()
        ax = fig.add_subplot(111)

    # Accept either the reliability_mc result dict or a bare sample array.
    if isinstance(result, dict):
        samples = np.asarray(result.get('fs_samples', []), dtype=float)
    else:
        samples = np.asarray(result, dtype=float)
        result = {}
    fs = samples[np.isfinite(samples)]

    if fs.size < 2:
        ax.text(0.5, 0.5, "Not enough admissible Monte Carlo samples to plot.",
                ha='center', va='center', transform=ax.transAxes)
        ax.set_axis_off()
        if own_fig:
            plt.show()
        return fig

    # Prefer the engine's reported statistics; fall back to the samples themselves
    # so a bare array still plots sensibly.
    mean_FS = float(result.get('mean_FS', np.mean(fs)))
    sigma_F = float(result.get('sigma_F', np.std(fs, ddof=1)))
    beta_normal = result.get('beta_normal')
    beta_ln = result.get('beta_ln')
    pf = result.get('pf_empirical', result.get('prob_failure'))
    if pf is None:
        pf = float(np.count_nonzero(fs < 1.0)) / fs.size

    ax.hist(fs, bins=bins, density=True, color='#4c78a8', alpha=0.75,
            edgecolor='white', linewidth=0.4, label='FS samples')

    # Shade the FS < 1 (failure) tail and mark the failure line + mean.
    ax.axvspan(min(fs.min(), 1.0), 1.0, color='#d62728', alpha=0.06, zorder=0)
    ax.axvline(1.0, color='#d62728', linewidth=1.6, linestyle='-',
               label='FS = 1 (failure)')
    ax.axvline(mean_FS, color='black', linewidth=1.3, linestyle='--',
               label=f'mean FS = {mean_FS:.3f}')

    # Fitted normal / lognormal overlays (method of moments), guarded against a
    # degenerate spread or a non-positive mean (lognormal only).
    if show_fits and sigma_F > 0:
        from scipy.stats import norm as _norm, lognorm as _lognorm
        xs = np.linspace(fs.min(), fs.max(), 400)
        ax.plot(xs, _norm.pdf(xs, mean_FS, sigma_F), color='#f58518', lw=1.8,
                label='normal fit')
        if mean_FS > 0:
            s_ln = np.sqrt(np.log(1.0 + (sigma_F / mean_FS) ** 2))
            scale = mean_FS / np.sqrt(1.0 + (sigma_F / mean_FS) ** 2)
            ax.plot(xs, _lognorm.pdf(xs, s_ln, scale=scale), color='#54a24b',
                    lw=1.8, linestyle='--', label='lognormal fit')

    ax.set_xlabel('Factor of Safety')
    ax.set_ylabel('Probability density')
    ax.grid(True, axis='y', alpha=0.3)

    if show_title:
        bits = [f"mean FS = {mean_FS:.3f}", f"$\\sigma_F$ = {sigma_F:.3f}"]
        beta_bits = []
        if beta_normal is not None and np.isfinite(beta_normal):
            beta_bits.append(f"$\\beta_{{normal}}$ = {beta_normal:.3f}")
        if beta_ln is not None and np.isfinite(beta_ln):
            beta_bits.append(f"$\\beta_{{ln}}$ = {beta_ln:.3f}")
        line2 = ", ".join(beta_bits + [f"$P_f$ = {pf * 100:.2f}%"])
        n_txt = ""
        if result.get('n_valid') is not None:
            n_txt = f"   (n = {result.get('n_valid')})"
        ax.set_title("Monte Carlo Reliability — FS distribution\n"
                     + ", ".join(bits) + n_txt + "\n" + line2)

    if show_legend:
        ax.legend(loc='best', fontsize=8, framealpha=0.9)

    fig.tight_layout()
    if save_png:
        fig.savefig('plot_reliability_histogram.png', dpi=dpi, bbox_inches='tight')
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

    # Plot 2D elements SECOND (middle layer). Their edge linewidth is set to a
    # placeholder here and finalized (adaptive_edge_linewidth) after layout, once
    # the data->display transform is known — so a dense mesh thins to a legible
    # hairline instead of fusing into a solid black fill.
    mesh_edge_collections = []
    mesh_edge_segments = []
    for mid, elements_list in material_elements.items():
        for ec in elements_list:                 # corner-edge segments for width calc
            for a in range(len(ec)):
                mesh_edge_segments.append([tuple(ec[a]), tuple(ec[(a + 1) % len(ec)])])
        # Create polygon collection for this material
        poly_collection = PolyCollection(elements_list, gid='MESH',
                                       facecolor=_mat_color(mid),
                                       edgecolor='k',
                                       alpha=0.4,
                                       linewidth=_EDGE_LW_BASE)
        ax.add_collection(poly_collection)
        mesh_edge_collections.append(poly_collection)

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

    # Density-aware element-edge width, now that the layout (and so the
    # data->display scale) is final.
    if mesh_edge_collections and mesh_edge_segments:
        _lw = adaptive_edge_linewidth(ax, fig, mesh_edge_segments)
        for _pc in mesh_edge_collections:
            _pc.set_linewidth(_lw)

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

    # Distributed loads bounds — the derived water loads included, so the material
    # table is not placed on top of a reservoir the user did not type.
    for dload_set in (list(slope_data.get('dloads') or [])
                      + _derived_water_blocks(slope_data)):
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


def _bar_chart_height_in(n_bars, per_bar_in=0.5, margin_in=1.4, min_in=2.4,
                         max_in=None):
    """Content-proportionate figure height for a categorical bar chart (tornado,
    and similarly-shaped sweep plots): a fixed allowance per bar (row pitch) plus
    fixed margins for the title/axis labels/legend, so a 4-bar chart stays
    compact and a 20-bar chart grows naturally.

    This is deliberately independent of any *viewport* height a caller (Studio's
    MplCanvas) might otherwise stretch the figure to: a sparse categorical bar
    chart isn't a dense drawing like a slope cross-section, so filling an
    arbitrary window height balloons the bars (in pixels) far past the fixed
    point-size of their labels — wrong proportions at any capture size. Width
    still follows the caller's figsize/figure (Studio still fits it to the
    viewport); only height is content-driven, here and in Studio alike.
    """
    h = margin_in + max(int(n_bars), 1) * per_bar_in
    h = max(h, min_in)
    if max_in is not None:
        h = min(h, max_in)
    return h


def plot_tornado(result, figsize=(8, 5), save_png=False, dpi=300, fig=None,
                 style=None, widest_on_top=True):
    """Duncan-style tornado diagram from tornado(): horizontal bars of the output
    swing between each parameter's low and high bound, sorted by span, with the
    base-case value as the vertical reference.

    The swept quantity follows result['output'] / result['output_label'] — a
    factor of safety for an LEM/FEM tornado, total discharge q for a seepage one
    (absent keys default to FS, so an old result plots as before).

    The figure HEIGHT is content-proportionate (see ``_bar_chart_height_in``) —
    a fixed allowance per bar plus fixed margins — not stretched to fill an
    arbitrary figsize/viewport height, so bar thickness stays in normal
    proportion to the tick/annotation text no matter how tall the window is.
    Width still follows ``figsize`` (standalone) or the caller's figure
    (Studio: the viewport width).

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

    # Build the bar list FIRST — both branches below need the bar count to size
    # the figure height content-proportionately.
    bars = []
    for param, g in df.loc[~df['is_base']].groupby('param'):
        ok = g.loc[g['success']].sort_values('value')
        if ok.empty:
            continue
        lo_fs, hi_fs = ok['fs'].iloc[0], ok['fs'].iloc[-1]
        bars.append((param, lo_fs, hi_fs, abs(hi_fs - lo_fs)))
    bars.sort(key=lambda b: b[3], reverse=not widest_on_top)  # barh draws k=0 at the bottom

    height_in = _bar_chart_height_in(len(bars))
    if fig is None:
        fig, ax = plt.subplots(figsize=(figsize[0], height_in))
    else:
        ax = fig.add_subplot(111)
        w_in, _ = fig.get_size_inches()
        fig.set_size_inches(w_in, height_in, forward=False)

    # Fat bars, tight row pitch — a categorical bar chart should read as mostly
    # bars, not gaps.
    bar_height = 0.78
    for k, (param, lo_fs, hi_fs, _span) in enumerate(bars):
        left, width = min(lo_fs, hi_fs), abs(hi_fs - lo_fs)
        ax.barh(k, width, left=left, height=bar_height, color='C0', alpha=0.75)
        ax.annotate(vfmt(lo_fs), (lo_fs, k), textcoords='offset points',
                    xytext=(-6, 0), ha='right', va='center', fontsize=8)
        ax.annotate(vfmt(hi_fs), (hi_fs, k), textcoords='offset points',
                    xytext=(6, 0), ha='left', va='center', fontsize=8)
    ax.set_yticks(range(len(bars)))
    ax.set_yticklabels([b[0] for b in bars])
    # Snug y-limits (not the default ~5% margins()) so the axis band is mostly
    # bars, matching the fattened bar_height above.
    if bars:
        pad = (1 - bar_height) / 2 + 0.1
        ax.set_ylim(-0.5 - pad, len(bars) - 0.5 + pad)
    # Explicit x limits, not margins(): the Studio canvas re-fits axes after render,
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


# ---------------------------------------------------------------------------
# Parametric-study plots (the tornado's companions): scaled-sensitivity bars,
# a variance-contribution Pareto, a spider plot, and Monte Carlo rank bars.
# ---------------------------------------------------------------------------

# Sign colors shared by the scaled bars and the rank bars: FS RISES with the
# parameter (stabilizing) vs FALLS. Matplotlib tab:green / tab:red, so they read
# as the same family as the rest of the palette.
_SIGN_POS_COLOR = '#2ca02c'   # FS increases with the parameter
_SIGN_NEG_COLOR = '#d62728'   # FS decreases with the parameter

# Scaled-sensitivity y-axis labels per scaling key.
_SCALED_META = {
    'elasticity': ('elasticity', 'FS elasticity  |∂F/∂p · p/F|'),
    'per_1pct':   ('per_1pct',   '|ΔFS| per 1% change in parameter'),
    'per_sigma':  ('per_sigma',  '|ΔFS| per 1σ change in parameter'),
}


def _param_short_label(bar_or_ref):
    """Compact parameter label for a bar/curve. Prefers an explicit 'label' the
    engine set (e.g. 'Soil · phi'), else derives one from the canonical ref."""
    if isinstance(bar_or_ref, dict):
        if bar_or_ref.get('label'):
            return bar_or_ref['label']
        ref = str(bar_or_ref.get('param', ''))
    else:
        ref = str(bar_or_ref)
    parts = ref.split(':')
    if len(parts) == 3:
        return f"{parts[1]} · {parts[2]}"
    return parts[-1] or ref


def _sign_legend(ax, output='FS'):
    """Two proxy handles explaining the bar colors (sign of the response)."""
    from matplotlib.patches import Patch
    up = Patch(facecolor=_SIGN_POS_COLOR, label=f'{output} increases with parameter')
    dn = Patch(facecolor=_SIGN_NEG_COLOR, label=f'{output} decreases with parameter')
    ax.legend(handles=[up, dn], fontsize=8, loc='best')


def plot_scaled_sensitivity(result, scaling='elasticity', figsize=(8, 5),
                            save_png=False, dpi=300, fig=None, style=None):
    """Scaled-sensitivity bars from ``scaled_sensitivity()`` — a vertical bar per
    parameter whose HEIGHT is the magnitude of the (scaled) local sensitivity and
    whose COLOR is the sign of the response (green: FS rises with the parameter;
    red: FS falls). Made comparable across parameters with different units.

    ``scaling`` picks which of the three coefficients the height shows:
      * ``'elasticity'`` (default) — dimensionless ∂F/∂p·p/F.
      * ``'per_1pct'``  — ΔFS for a 1% change in the parameter.
      * ``'per_sigma'`` — ΔFS for a one-σ change (only parameters carrying a σ are
        drawn; the others have no per-σ coefficient).

    The derivative behind every scaling is a central difference at ±1% (relative)
    about the base case (see ``scaled_sensitivity``). Bars are sorted by magnitude,
    largest at the left.
    """
    if fig is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        ax = fig.add_subplot(111)
    key, ylabel = _SCALED_META.get(scaling, _SCALED_META['elasticity'])
    output = result.get('output', 'FS')
    bars = [b for b in result.get('bars', []) if b.get(key) is not None]
    bars.sort(key=lambda b: abs(b[key]), reverse=True)
    if not bars:
        ax.text(0.5, 0.5, f"No parameters carry a '{scaling}' coefficient.",
                ha='center', va='center', transform=ax.transAxes)
        ax.set_axis_off()
        return fig
    heights = [abs(b[key]) for b in bars]
    colors = [_SIGN_POS_COLOR if b[key] >= 0 else _SIGN_NEG_COLOR for b in bars]
    labels = [_param_short_label(b) for b in bars]
    x = range(len(bars))
    ax.bar(x, heights, color=colors, alpha=0.85, width=0.6)
    for xi, h in zip(x, heights):
        ax.annotate(f'{h:.3g}', (xi, h), textcoords='offset points',
                    xytext=(0, 3), ha='center', va='bottom', fontsize=8)
    ax.set_xticks(list(x))
    ax.set_xticklabels(labels, rotation=20, ha='right')
    ax.set_ylabel(ylabel)
    ax.set_ylim(0, max(heights) * 1.18)
    scal_word = {'elasticity': 'elasticity', 'per_1pct': 'per 1%',
                 'per_sigma': 'per σ'}.get(scaling, scaling)
    fbase = result.get('fs_base')
    sub = f"   (base {output} = {fbase:.3g})" if fbase is not None else ""
    ax.set_title(f'Scaled sensitivity — {scal_word}{sub}')
    _sign_legend(ax, output=output)
    ax.grid(True, axis='y', alpha=0.4)
    fig.tight_layout()
    if save_png:
        fig.savefig(f'scaled_sensitivity_{scaling}.png', dpi=dpi, bbox_inches='tight')
    return fig


def plot_variance_pareto(result, figsize=(8, 5), save_png=False, dpi=300, fig=None,
                         style=None):
    """Variance-contribution Pareto from ``variance_contribution()``: each uncertain
    parameter's share of Var(FS) as a descending bar, with the cumulative share as
    an overlaid line — the standard way to read *which* uncertainties actually drive
    the reliability. Bars are labelled '% of FS variance'; the terms are the
    Taylor-series (dF/dp·σ_p)² contributions, so they sum to 100%.
    """
    if fig is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        ax = fig.add_subplot(111)
    bars = list(result.get('bars', []))
    if not bars:
        ax.text(0.5, 0.5, "No uncertain parameters (no σ set).",
                ha='center', va='center', transform=ax.transAxes)
        ax.set_axis_off()
        return fig
    pct = [b['pct'] for b in bars]
    cum = [b['cumulative'] for b in bars]
    labels = [_param_short_label(b) for b in bars]
    x = list(range(len(bars)))
    ax.bar(x, pct, color='C0', alpha=0.8, width=0.6, label='% of Var(FS)')
    for xi, p in zip(x, pct):
        ax.annotate(f'{p:.1f}%', (xi, p), textcoords='offset points',
                    xytext=(0, 3), ha='center', va='bottom', fontsize=8)
    ax2 = ax.twinx()
    ax2.plot(x, cum, color=_SIGN_NEG_COLOR, marker='o', lw=1.5, label='cumulative')
    ax2.set_ylim(0, 105)
    ax2.set_ylabel('Cumulative %')
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=20, ha='right')
    ax.set_ylabel('% of FS variance')
    ax.set_ylim(0, max(pct) * 1.18 if pct else 1)
    sig = result.get('sigma_F')
    cov = result.get('COV_F')
    bits = []
    if sig is not None:
        bits.append(f'σ_F = {sig:.3g}')
    if cov is not None:
        bits.append(f'COV_F = {cov:.3g}')
    sub = ('   (' + ', '.join(bits) + ')') if bits else ''
    ax.set_title(f'Variance contribution to FS (Taylor series){sub}')
    ax.grid(True, axis='y', alpha=0.4)
    fig.tight_layout()
    if save_png:
        fig.savefig('variance_pareto.png', dpi=dpi, bbox_inches='tight')
    return fig


def plot_spider(sweeps, x_mode='pct', figsize=(8, 5), save_png=False, dpi=300,
                fig=None, style=None, target_fs=None):
    """Spider (radial-style) plot: FS versus each parameter swept across its range,
    all curves sharing one normalized x-axis so a family of parameters can be read
    on a single chart. One labelled curve per parameter, a black base-case marker at
    the origin (the unperturbed model), and an FS = 1 guide.

    ``sweeps`` is a ``{ref: DataFrame}`` mapping of ``sensitivity()`` sweeps (or a
    result dict carrying a ``'sweeps'`` key). ``x_mode='pct'`` (default) plots the
    percent change of each parameter from its base value; ``'value'`` plots the raw
    swept value (only sensible when the parameters share units).
    """
    if isinstance(sweeps, dict) and 'sweeps' in sweeps:
        sweeps = sweeps['sweeps']
    items = list(sweeps.items()) if isinstance(sweeps, dict) else list(sweeps)
    if fig is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        ax = fig.add_subplot(111)
    base_fs = None
    output = 'FS'
    output_label = 'Factor of Safety'
    for ref, df in items:
        if df is None or len(df) == 0:
            continue
        if 'output' in df.columns:
            output = df['output'].iloc[0]
            output_label = df['output_label'].iloc[0]
        base_row = df.loc[df['is_base'] & df['success']]
        if base_fs is None and len(base_row):
            base_fs = float(base_row['fs'].iloc[0])
        pts = df.loc[~df['is_base'] & df['success']].sort_values('value')
        if pts.empty:
            continue
        if x_mode == 'pct':
            if 'rel' in pts.columns and pts['rel'].notna().all():
                xv = (pts['rel'].to_numpy(dtype=float) - 1.0) * 100.0
            elif len(base_row):
                bv = float(base_row['value'].iloc[0])
                xv = (pts['value'].to_numpy(dtype=float) / bv - 1.0) * 100.0 \
                    if bv else pts['value'].to_numpy(dtype=float)
            else:
                xv = pts['value'].to_numpy(dtype=float)
        else:
            xv = pts['value'].to_numpy(dtype=float)
        ax.plot(xv, pts['fs'], marker='o', ms=4, label=_param_short_label(ref))
    if base_fs is not None:
        bx = 0.0 if x_mode == 'pct' else np.nan
        if np.isfinite(bx):
            ax.plot(bx, base_fs, 's', color='k', ms=9, zorder=6,
                    label=f'base case ({output} = {base_fs:.3g})')
    if output == 'FS':
        ax.axhline(1.0, color='r', linestyle='--', linewidth=0.8, label='FS = 1')
    if target_fs is not None:
        ax.axhline(target_fs, color='gray', linestyle='--', linewidth=0.8,
                   label=f'{output} = {target_fs:g}')
    ax.set_xlabel('Change from base case (%)' if x_mode == 'pct'
                  else 'Parameter value')
    ax.set_ylabel(output_label)
    ax.set_title(f'Spider: {output} vs each parameter')
    ax.legend(fontsize=8, ncol=2 if len(items) > 6 else 1)
    ax.grid(True, alpha=0.4)
    fig.tight_layout()
    if save_png:
        fig.savefig('spider.png', dpi=dpi, bbox_inches='tight')
    return fig


def plot_mc_rank_correlation(result, figsize=(8, 5), save_png=False, dpi=300,
                             fig=None, style=None):
    """Monte Carlo rank-correlation bars from ``mc_rank_correlation()``: the Spearman
    correlation between each sampled input and FS, one horizontal bar per parameter,
    signed (green +, red −) and sorted by magnitude with the strongest on top.

    This is a GLOBAL sensitivity measure — it ranks parameters over the whole
    sampled distribution, so it captures a parameter's spread and any nonlinearity —
    and is meant to be read alongside the LOCAL scaled-sensitivity bars, which are a
    derivative at the base case. The two can disagree, and the disagreement is
    informative.
    """
    if fig is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        ax = fig.add_subplot(111)
    bars = [b for b in result.get('bars', [])
            if b.get('rho') is not None and np.isfinite(b['rho'])]
    # barh draws index 0 at the bottom, so sort ascending -> strongest on top.
    bars = sorted(bars, key=lambda b: abs(b['rho']))
    if not bars:
        ax.text(0.5, 0.5, "No admissible Monte Carlo correlations.",
                ha='center', va='center', transform=ax.transAxes)
        ax.set_axis_off()
        return fig
    y = list(range(len(bars)))
    vals = [b['rho'] for b in bars]
    colors = [_SIGN_POS_COLOR if v >= 0 else _SIGN_NEG_COLOR for v in vals]
    ax.barh(y, vals, color=colors, alpha=0.85, height=0.6)
    for yi, v in zip(y, vals):
        ax.annotate(f'{v:+.2f}', (v, yi), textcoords='offset points',
                    xytext=(6 if v >= 0 else -6, 0),
                    ha='left' if v >= 0 else 'right', va='center', fontsize=8)
    ax.axvline(0.0, color='k', linewidth=0.8)
    ax.set_xlim(-1.05, 1.05)
    ax.set_yticks(y)
    ax.set_yticklabels([_param_short_label(b) for b in bars])
    ax.set_xlabel('Spearman rank correlation with FS')
    n_valid = result.get('n_valid')
    n_samples = result.get('n_samples')
    sub = (f'   ({n_valid}/{n_samples} valid samples)'
           if n_valid is not None and n_samples is not None else '')
    ax.set_title(f'Monte Carlo rank correlation — global sensitivity{sub}')
    ax.grid(True, axis='x', alpha=0.4)
    fig.tight_layout()
    if save_png:
        fig.savefig('mc_rank_correlation.png', dpi=dpi, bbox_inches='tight')
    return fig
