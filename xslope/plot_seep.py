import logging
import math

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.tri as tri
from matplotlib.collections import LineCollection, PatchCollection
from matplotlib.patches import Polygon
from matplotlib.ticker import MaxNLocator
import numpy as np

from . import colormaps as _colormaps  # noqa: F401  (registers the BGYR ramp by name)

logger = logging.getLogger(__name__)

# Below this peak-to-peak spread the stream function phi is treated as flat (no
# through-flow → no flow net). Guards the flow-line contour so a degenerate transient
# storage-release frame can't raise "Contour levels must be increasing". A real
# (steady or through-flowing transient) frame has phi range orders of magnitude larger.
_PHI_FLAT_TOL = 1e-9


def flownet_has_phreatic(seep_data, solution):
    """Whether :func:`plot_seep_solution` draws a phreatic surface for this
    solution (with its default ``phreatic=True``).

    A phreatic surface is the p = 0 contour of an UNCONFINED solve, where the water
    table is located as part of the answer. A confined solve is a single saturated
    Laplace solve whose head is a potential, not a water level: negative pore
    pressure is routine there and the p = 0 contour is an artifact. A solve that
    does not record which it was is read as unconfined, the plotter's long-standing
    default. And an unconfined field that never goes into suction has no p = 0
    contour to draw.

    Exported so a caption or a paragraph can name what the figure carries instead
    of asserting a line the figure may not contain.
    """
    if not (solution or {}).get("unconfined", True):
        return False
    u = (solution or {}).get("u")
    if u is None:
        return False
    u = np.asarray(u, dtype=float)
    return bool(u.size and np.min(u) < 0)


def flownet_has_flowlines(seep_data, solution):
    """Whether :func:`plot_seep_solution` draws flow LINES for this solution (with
    its default ``flowlines=True``, plotting head).

    Flow lines are contours of the stream function, spaced by the flow they carry,
    so they need three things: a stream function with real range, a total flowrate
    to space them by, and the conductivities the spacing is computed from. A
    solution read back from a file that records no flowrate — or a pore-pressure
    field imported from another program, which carries no stream function at all —
    has no flow lines, and the figure is head contours alone.

    Exported so a caption or a paragraph can name what the figure carries.
    """
    if (seep_data or {}).get("k1_by_mat") is None:
        return False
    solution = solution or {}
    if solution.get("flowrate") is None:
        return False
    phi = solution.get("phi")
    if phi is None:
        return False
    phi = np.asarray(phi, dtype=float)
    if not phi.size:
        return False
    # RELATIVE to phi's own magnitude, not an absolute floor: phi carries the units
    # of the flowrate, so a low-conductivity problem has a small phi range that is
    # nonetheless the whole flow net (Rocscience GW#5 runs at k = 1e-10 m/s and spans
    # phi 0 -> 8e-11, well under a fixed 1e-9, yet is perfectly resolved). A
    # degenerate frame has phi exactly constant, which fails the relative test at any
    # scale.
    scale = float(np.max(np.abs(phi)))
    spread = float(np.ptp(phi))
    return bool(spread > max(_PHI_FLAT_TOL * scale, 0.0) and spread > 0.0)


def flownet_base_material(seep_data, solution, levels=20):
    """The material a flow net should be scaled to, as a 1-based id.

    ``plot_seep_solution`` draws ``q·(levels - 1) / (k_base·Δh)`` flow channels,
    so the base material decides how dense the net is and nothing else: a zone
    whose conductivity is far above the one carrying the through-flow leaves a
    net with no flow lines in it, and one far below asks for hundreds.

    The choice is the one the shipped sample figures were built on — "the zone
    that controls the through-flow, so the net is readable"
    (:file:`tools/make_seep_sample_figures.py`) — stated as the arithmetic it
    always was: take the zone whose channel count lands nearest a readable one.
    A flow net reads as curvilinear squares, so a readable net has about half as
    many flow channels as it has potential drops, and the target follows the
    contour count rather than being a number of its own. Nearest is measured on
    the ratio, because the candidates differ by orders of magnitude.

    Reproduces every base material the sample figures name: the low-conductivity
    core of each zoned dam, the foundation under Johnson's near-cutoff core, and
    the single through-flow layer of the one-material problems.

    Returns 1 where the flow rate or the head drop cannot be read, which is what
    the parameter has always defaulted to.
    """
    k1 = (seep_data or {}).get("k1_by_mat")
    if k1 is None or len(k1) == 0:
        return 1
    k2 = (seep_data or {}).get("k2_by_mat")
    head = (solution or {}).get("head")
    q = (solution or {}).get("flowrate")
    if head is None or q is None:
        return 1
    hdrop = float(np.max(head) - np.min(head))
    q = float(q)
    if not (hdrop > 0 and q > 0):
        return 1

    drops = max(int(levels) - 1, 1)
    target = drops / 2.0
    best, best_miss = 1, None
    for i in range(len(k1)):
        # The equivalent conductivity of an anisotropic zone, the same
        # sqrt(k1·k2) the flow-line count itself is computed from.
        k = float(np.sqrt(float(k1[i]) * float(k2[i] if k2 is not None else k1[i])))
        if not k > 0:
            continue
        channels = q * drops / (k * hdrop)
        if not channels > 0:
            continue
        miss = abs(math.log(channels / target))
        if best_miss is None or miss < best_miss:
            best, best_miss = i + 1, miss
    return best


def _draw_seep_bc_levels(ax, seep_data, solution, style):
    """Overlay each head/reservoir boundary's instantaneous water level on a seep
    solution frame (the ``show_bc_levels`` reading aid).

    Reuses the SAME series evaluator, reservoir-surface clip, AND water-level symbol
    the inputs BC rendering uses (``xslope.plot``), so a level is read and drawn
    exactly as the inputs plot draws it: a series value is evaluated at this frame's
    time (``solution['time']``, or 0 for a steady solution), a constant value is drawn
    flat. Each level is a thin waterline in the boundary's LIGHT shade plus the shared
    apex-down water symbol (``draw_water_level_symbol`` — tip on the line, sized in
    points so it is identical to the inputs symbol on any domain). A reservoir
    (submerged-only) face is clipped to where the level meets the face; a plain head is
    drawn across its boundary's x-extent. Returns the legend handles for the drawn
    level types (empty when nothing was drawn).

    Draws nothing (returns []) when the seep_data carries no BC geometry — e.g. a
    solution reconstructed from JSON — so the overlay degrades gracefully.
    """
    from .plot import (_eval_bc_series_at, _reservoir_surface_x,
                       seep_bc_level_color, draw_water_level_symbol)

    seepage_bc = seep_data.get("seepage_bc") or {}
    heads = seepage_bc.get("specified_heads") or []
    if not heads:
        return []
    set_no = int(seep_data.get("seep_bc", 1) or 1)
    tseep = seep_data.get("tseep")
    # Frame time: transient frames carry it; a steady solution has none → t = 0, at
    # which a series holds its first value and a constant is itself.
    t = solution.get("time")
    t = 0.0 if t is None else float(t)

    drawn = {}   # kind -> light color, for one legend entry per drawn type
    for h in heads:
        coords = h.get("coords") or []
        if len(coords) < 2:
            continue
        kind = str(h.get("kind", "head")).strip().lower()
        val = h.get("head")
        if isinstance(val, str):
            level = _eval_bc_series_at(tseep, val, t)
        else:
            try:
                level = float(val)
            except (TypeError, ValueError):
                level = None
        if level is None:
            continue

        color = seep_bc_level_color(style, set_no, kind)
        if kind == "reservoir":
            x_up, x_face = _reservoir_surface_x(coords, float(level))
        else:
            xs = [c[0] for c in coords]
            x_up, x_face = min(xs), max(xs)
        if x_face <= x_up:                       # degenerate (level meets a corner)
            continue
        ax.plot([x_up, x_face], [level, level], color=color, lw=2.0,
                linestyle="-", zorder=6, gid="BC_LEVEL")
        draw_water_level_symbol(ax, 0.5 * (x_up + x_face), float(level),
                                color=color, markersize=8, extra_gap_points=2.0)
        drawn.setdefault(kind, color)

    handles = []
    for kind, color in drawn.items():
        label = "Reservoir level" if kind == "reservoir" else "Head level"
        handles.append(plt.Line2D([0], [0], color=color, lw=2.0, label=label))
    return handles


def plot_seep_data(seep_data, figsize=(12, 7), show_nodes=False, show_bc=False, label_elements=False, label_nodes=False, alpha=0.6, save_png=False, save_dxf=False, dpi=300, legend_ncol="auto", legend_frame=False, show_title=True, show_legend=True, fig=None, style=None):
    """
    Plots a mesh colored by material zone.
    Supports both triangular and quadrilateral elements.

    Args:
        seep_data: Dictionary containing seep data from import_seep2d
        show_nodes: If True, plot node points
        show_bc: If True, plot boundary condition nodes
        label_elements: If True, label each element with its number at its centroid
        label_nodes: If True, label each node with its number just above and to the right
        fig: Optional existing Matplotlib Figure to draw into (embeds in a GUI canvas);
            when provided the figure is cleared/reused and plt.show() is skipped.
    """

    # Extract data from seep_data
    nodes = seep_data["nodes"]
    elements = seep_data["elements"]
    element_materials = seep_data["element_materials"]
    element_types = seep_data.get("element_types", None)
    bc_type = seep_data["bc_type"]

    own_fig = fig is None
    if own_fig:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig.clear()
        ax = fig.add_subplot(111)
    materials = np.unique(element_materials)

    # Material colors (style overrides → palette default). Mesh material IDs are
    # 1-based (gmsh); the style sheet / inputs key by 0-based mat_id, so map mat-1
    # — this also aligns the zone colors with the Inputs view.
    from .style import resolve_style, material_style, feature_style
    _st = resolve_style(style)
    mat_to_color = {mat: material_style(_st, int(mat) - 1)["color"] for mat in materials}

    # If element_types is not provided, assume all triangles (backward compatibility)
    if element_types is None:
        element_types = np.full(len(elements), 3)

    # Batch polygons and edge lines by material for efficient rendering
    mat_fill_polys = {mat: [] for mat in materials}
    edge_segments = []

    for idx, element_nodes in enumerate(elements):
        element_type = element_types[idx]
        mat = element_materials[idx]

        if element_type == 3:  # Linear triangle
            mat_fill_polys[mat].append(nodes[element_nodes[:3]])

        elif element_type == 6:  # Quadratic triangle - subdivide into 4 sub-triangles
            n0, n1, n2 = nodes[element_nodes[0]], nodes[element_nodes[1]], nodes[element_nodes[2]]
            n3, n4, n5 = nodes[element_nodes[3]], nodes[element_nodes[4]], nodes[element_nodes[5]]
            mat_fill_polys[mat].extend([
                np.array([n0, n3, n5]),
                np.array([n3, n1, n4]),
                np.array([n5, n4, n2]),
                np.array([n3, n4, n5]),
            ])
            edge_segments.extend([[n0, n1], [n1, n2], [n2, n0]])

        elif element_type == 4:  # Linear quadrilateral
            mat_fill_polys[mat].append(nodes[element_nodes[:4]])

        elif element_type == 8:  # Quadratic quadrilateral - subdivide into 4 sub-quads
            n0, n1, n2, n3 = nodes[element_nodes[0]], nodes[element_nodes[1]], nodes[element_nodes[2]], nodes[element_nodes[3]]
            n4, n5, n6, n7 = nodes[element_nodes[4]], nodes[element_nodes[5]], nodes[element_nodes[6]], nodes[element_nodes[7]]
            center = np.array([(n0[0]+n1[0]+n2[0]+n3[0]+n4[0]+n5[0]+n6[0]+n7[0]) / 8,
                               (n0[1]+n1[1]+n2[1]+n3[1]+n4[1]+n5[1]+n6[1]+n7[1]) / 8])
            mat_fill_polys[mat].extend([
                np.array([n0, n4, center, n7]),
                np.array([n4, n1, n5, center]),
                np.array([center, n5, n2, n6]),
                np.array([n7, center, n6, n3]),
            ])
            edge_segments.extend([[n0, n1], [n1, n2], [n2, n3], [n3, n0]])

        elif element_type == 9:  # 9-node quadrilateral
            n0, n1, n2, n3 = nodes[element_nodes[0]], nodes[element_nodes[1]], nodes[element_nodes[2]], nodes[element_nodes[3]]
            n4, n5, n6, n7 = nodes[element_nodes[4]], nodes[element_nodes[5]], nodes[element_nodes[6]], nodes[element_nodes[7]]
            center = nodes[element_nodes[8]]
            mat_fill_polys[mat].extend([
                np.array([n0, n4, center, n7]),
                np.array([n4, n1, n5, center]),
                np.array([center, n5, n2, n6]),
                np.array([n7, center, n6, n3]),
            ])
            edge_segments.extend([[n0, n1], [n1, n2], [n2, n3], [n3, n0]])

    # Render filled polygons as batched PatchCollections (one per material)
    for mat in materials:
        polys = mat_fill_polys[mat]
        if not polys:
            continue
        has_edge = any(element_types[i] in (3, 4) for i, m in enumerate(element_materials) if m == mat)
        has_no_edge = any(element_types[i] in (6, 8, 9) for i, m in enumerate(element_materials) if m == mat)
        color = mat_to_color[mat]

        if has_edge and not has_no_edge:
            patch_list = [Polygon(p) for p in polys]
            pc = PatchCollection(patch_list, facecolor=color, edgecolor='k', linewidth=0.5, alpha=alpha, gid='MESH_FILL')
            ax.add_collection(pc)
        elif has_no_edge and not has_edge:
            patch_list = [Polygon(p) for p in polys]
            pc = PatchCollection(patch_list, facecolor=color, edgecolor='none', alpha=alpha, gid='MESH_FILL')
            ax.add_collection(pc)
        else:
            linear_polys = []
            sub_polys = []
            sub_idx = 0
            for i in range(len(elements)):
                if element_materials[i] != mat:
                    continue
                et = element_types[i]
                if et in (3, 4):
                    linear_polys.append(polys[sub_idx]); sub_idx += 1
                elif et in (6, 8, 9):
                    sub_polys.extend(polys[sub_idx:sub_idx+4]); sub_idx += 4
            if linear_polys:
                pc = PatchCollection([Polygon(p) for p in linear_polys], facecolor=color, edgecolor='k', linewidth=0.5, alpha=alpha, gid='MESH_FILL')
                ax.add_collection(pc)
            if sub_polys:
                pc = PatchCollection([Polygon(p) for p in sub_polys], facecolor=color, edgecolor='none', alpha=alpha, gid='MESH_FILL')
                ax.add_collection(pc)

    # Render outer boundary edges of quadratic elements as a single LineCollection
    if edge_segments:
        lc = LineCollection(edge_segments, colors='k', linewidths=0.5, gid='MESH')
        ax.add_collection(lc)

    # Label element numbers at centroids if requested
    if label_elements:
        for idx, element_nodes in enumerate(elements):
            element_type = element_types[idx]
            if element_type == 3:
                element_coords = nodes[element_nodes[:3]]
            elif element_type == 4:
                element_coords = nodes[element_nodes[:4]]
            elif element_type == 6:
                element_coords = nodes[element_nodes[:6]]
            elif element_type == 8:
                element_coords = nodes[element_nodes[:8]]
            else:
                element_coords = nodes[element_nodes[:9]]
            centroid = np.mean(element_coords, axis=0)
            ax.text(centroid[0], centroid[1], str(idx+1),
                    ha='center', va='center', fontsize=6, color='black', alpha=0.4,
                    zorder=10)

    if show_nodes:
        ax.plot(nodes[:, 0], nodes[:, 1], 'k.', markersize=0.75, gid='MESH_NODES')

    # Label node numbers if requested
    if label_nodes:
        for i, (x, y) in enumerate(nodes):
            ax.text(x + 0.5, y + 0.5, str(i+1), fontsize=6, color='blue', alpha=0.7,
                    ha='left', va='bottom', zorder=11)

    # Get material names if available
    material_names = seep_data.get("material_names", [])
    
    legend_handles = []
    for mat in materials:
        # Use material name if available, otherwise use "Material {mat}"
        if material_names and mat <= len(material_names):
            label = material_names[mat - 1]  # Convert to 0-based index
        else:
            label = f"Material {mat}"
        
        legend_handles.append(
            mpatches.Patch(facecolor=mat_to_color[mat], alpha=alpha,
                           edgecolor="none", label=label)
        )

    if show_bc:
        bc1 = nodes[bc_type == 1]
        bc2 = nodes[bc_type == 2]
        if len(bc1) > 0:
            h1, = ax.plot(bc1[:, 0], bc1[:, 1], 'bs', label="Fixed Head (bc_type=1)", gid='SEEP_FIXED_HEAD')
            legend_handles.append(h1)
        if len(bc2) > 0:
            h2, = ax.plot(bc2[:, 0], bc2[:, 1], 'ro', label="Exit Face (bc_type=2)", gid='SEEP_EXIT_FACE')
            legend_handles.append(h2)
        # Flux nodes are NOT Dirichlet, so they carry bc_type 0 and have to be found from
        # the assembled nodal loads instead. Inflow and outflow get opposite-pointing
        # markers: a flipped sign is the classic way to get a flux boundary wrong, and
        # this makes it visible rather than something you infer from the head field.
        flux_nodal = seep_data.get("flux_nodal")
        if flux_nodal is not None:
            fn = np.asarray(flux_nodal, dtype=float)
            fs = feature_style(_st, "seep_flux")
            fcolor = fs.get("color", "darkgreen")
            for mask, marker, what, gid in (
                    (fn > 0, '^', 'in', 'SEEP_FLUX_IN'),
                    (fn < 0, 'v', 'out', 'SEEP_FLUX_OUT')):
                pts = nodes[mask]
                if len(pts) > 0:
                    hf, = ax.plot(pts[:, 0], pts[:, 1], marker=marker, linestyle='none',
                                  color=fcolor, markeredgecolor='k', markeredgewidth=0.4,
                                  label=f"Specified Flux ({what}flow)", gid=gid)
                    legend_handles.append(hf)

    from .plot import _legend_below
    # Box-adjust (not datalim): seepage domains are typically wide and thin, so
    # filling the figure via datalim would squish the mesh into a tiny band.
    ax.set_aspect("equal")

    # Add a bit of headroom so the mesh/BC markers don't touch the top border
    y0, y1 = ax.get_ylim()
    if y1 > y0:
        pad = 0.05 * (y1 - y0)
        ax.set_ylim(y0, y1 + pad)

    # Count element types for title. element_types holds the NODE COUNT, so the
    # quadratic families are 6 (tri6) and 8/9 (quad8/quad9) — counting only 3 and 4
    # reported "0 triangles" on every higher-order mesh.
    num_triangles = int(np.sum(np.isin(element_types, (3, 6))))
    num_quads = int(np.sum(np.isin(element_types, (4, 8, 9))))
    if num_triangles > 0 and num_quads > 0:
        title = f"Finite Element Mesh with Material Zones ({num_triangles} triangles, {num_quads} quads)"
    elif num_quads > 0:
        title = f"Finite Element Mesh with Material Zones ({num_quads} quadrilaterals)"
    else:
        title = f"Finite Element Mesh with Material Zones ({num_triangles} triangles)"

    if show_title:
        ax.set_title(title)
    fig.tight_layout()
    # Single combined legend below the plot, after tight_layout so the reserved
    # bottom margin (for multi-row legends) isn't clobbered.
    _legend_below(ax, fig, handles=legend_handles,
                  legend_ncol=legend_ncol, frameon=legend_frame, show_legend=show_legend)

    base_name = 'plot_' + title.lower().replace(' ', '_').replace(':', '').replace(',', '').replace('(', '').replace(')', '')
    if save_png:
        fig.savefig(base_name + '.png', dpi=dpi, bbox_inches='tight')
    if save_dxf:
        from .cad import axes_to_dxf
        axes_to_dxf(ax, base_name + '.dxf')

    if own_fig:
        plt.show()
    return fig


def plot_seep_solution(seep_data, solution, figsize=(12, 7), levels=20, base_mat=1, fill_contours=True, phreatic=True, alpha=0.4, pad_frac=0.05, mesh=True, variable="head", vectors=False, vector_scale=0.05, flowlines=True, cmap="Spectral_r", cbar_shrink=0.8, vmin=None, vmax=None, save_png=False, save_dxf=False, dpi=300, legend_ncol="auto", legend_frame=False, show_title=True, show_legend=True, show_bc_levels=False, fig=None, style=None):
    """
    Plot seep analysis results including head contours, flow lines, and phreatic surface.
    
    This function visualizes the results of a seep analysis by plotting contours of various
    nodal variables (head, pore pressure, velocity magnitude, or gradient magnitude). When
    plotting head, flow lines are also overlaid. The plot properly handles mesh aspect ratios
    and supports both linear and quadratic triangular and quadrilateral elements.
    
    Parameters:
    -----------
    seep_data : dict
        Dictionary containing seep mesh data from import_seep2d. Required keys include:
        'nodes', 'elements', 'element_materials', 'element_types' (optional), and
        'k1_by_mat' (optional, for flowline calculation).
    solution : dict
        Dictionary containing solution results from run_seepage_analysis. Required keys include:
        'head' (array of total head values at nodes), 'velocity' (array of velocity vectors),
        'gradient' (array of hydraulic gradient vectors). Optional keys: 'phi' (stream function),
        'flowrate' (total flow rate), 'u' (pore pressure), 'v_mag' (velocity magnitude),
        'i_mag' (gradient magnitude).
    figsize : tuple of float, optional
        Figure size in inches (width, height). Default is (14, 6).
    levels : int, optional
        Number of contour levels to plot. Default is 20.
    base_mat : int, optional
        Material ID (1-based) used to compute hydraulic conductivity for flow function
        calculation. Default is 1. Only used when variable="head".
    fill_contours : bool, optional
        If True, shows filled contours with color map. If False, only black solid
        contour lines are shown. Default is True.
    phreatic : bool, optional
        If True, plots the phreatic surface (where pressure head = 0) as a thick red line.
        Default is True. Only plotted if pore pressure is negative somewhere in the domain.
    alpha : float, optional
        Transparency level (0-1) for material zone fill colors. Default is 0.4.
    pad_frac : float, optional
        Fraction of mesh extent to add as padding around the plot boundaries. Default is 0.05.
    mesh : bool, optional
        If True, overlays element edges in light gray. Default is True.
    variable : str, optional
        Nodal variable to contour. Options: "head" (default), "u" (pore pressure),
        "v_mag" (velocity magnitude), "i_mag" (gradient magnitude). When "head" is selected,
        flow lines can be overlaid if flowlines=True. Other variables do not include flow lines.
    vectors : bool, optional
        If True, plots velocity vectors as arrows at each node. Default is False.
    vector_scale : float, optional
        Scale factor for vector lengths. Maximum vector length will be x_range * vector_scale,
        where x_range is the x-extent of the mesh. Default is 0.05.
    flowlines : bool, optional
        If True and variable="head", overlays flow lines (stream function contours) on the plot.
        Default is True. Only applicable when variable="head". A flow net requires
        divergence-free through-flow, so it exists for a steady solution but NOT for a
        transient (storage-release) frame, whose flow is sourced from released storage
        and has no stream function; request velocity vectors to read the instantaneous
        flow direction of a transient frame instead. (A degenerate-phi transient frame
        for which flow lines are requested simply draws none — the ptp guard below —
        rather than raising.)
    vmin, vmax : float, optional
        Fixed lower/upper bounds for the contour levels and the colorbar. Default None
        auto-scales each render to its own data (byte-identical to before). Passing an
        explicit range pins a frame to a shared scale — e.g. a transient series passes
        the full-pool (t = 0) head range to every panel so the colorbars match and the
        drawdown reads as one continuous story instead of each frame re-normalizing.
    show_bc_levels : bool, optional
        If True, overlay each specified-head / reservoir boundary's WATER LEVEL for
        this frame — a thin waterline in the boundary's light shade plus an apex-down
        water symbol (tip on the line), clipped to the reservoir side of the face. For
        a transient frame the level is the boundary's series evaluated at the frame's
        time (``solution['time']``), so the pool visibly drops through a playback; for
        a steady solution (or a constant boundary) it is the constant level. This is a
        reading aid, NOT the dashed boundary-location line drawn by the inputs view.
        Default is False, which keeps every committed steady figure byte-identical.
        Requires the BC geometry carried in ``seep_data`` (``build_seep_data`` adds it);
        a seep_data without it simply draws nothing.

    Returns:
    --------
    None
        Displays the plot using matplotlib.pyplot.show().
    
    Notes:
    ------
    - The function automatically subdivides quadratic elements (tri6, quad8, quad9) for
      proper visualization and contouring.
    - Flowlines are only plotted when variable="head" and if 'phi' and 'flowrate' are present
      in solution and 'k1_by_mat' is present in seep_data.
    - The plot includes a colorbar for contours when fill_contours=True.
    - The title includes flowrate information if available in the solution dictionary and
      variable="head".
    """
    # Validate variable parameter
    valid_variables = ["head", "u", "v_mag", "i_mag"]
    if variable not in valid_variables:
        raise ValueError(f"variable must be one of {valid_variables}, got '{variable}'")
    
    # Extract data from seep_data and solution
    nodes = seep_data["nodes"]
    elements = seep_data["elements"]
    element_materials = seep_data["element_materials"]
    element_types = seep_data.get("element_types", None)  # New field for element types
    k1_by_mat = seep_data.get("k1_by_mat")  # Use .get() in case it's not present
    k2_by_mat = seep_data.get("k2_by_mat")
    
    # Extract the variable to plot
    if variable not in solution:
        raise ValueError(f"Variable '{variable}' not found in solution dictionary. Available keys: {list(solution.keys())}")
    contour_data = solution[variable]
    
    # Extract head and flowline-related data (only needed for head plots)
    head = solution["head"]
    phi = solution.get("phi")
    flowrate = solution.get("flowrate")
    
    # Determine if we should plot flow lines (only for head and if flowlines=True)
    plot_flowlines = (variable == "head" and flowlines)


    # No layout engine: the colorbar goes in an explicit divider axes (below), so
    # _legend_below can set manual margins like every other plot. (A constrained-
    # layout colorbar is incompatible with those margins — switching engines once
    # a gridspec colorbar exists raises RuntimeError.)
    own_fig = fig is None
    if own_fig:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig.clear()
        try:
            fig.set_layout_engine("none")
        except Exception:
            pass
        ax = fig.add_subplot(111)

    # If element_types is not provided, assume all triangles (backward compatibility)
    if element_types is None:
        element_types = np.full(len(elements), 3)
    
    # Count element types
    tri3_count = np.sum(element_types == 3)
    tri6_count = np.sum(element_types == 6) 
    quad4_count = np.sum(element_types == 4)
    quad8_count = np.sum(element_types == 8)
    quad9_count = np.sum(element_types == 9)
    
    logger.debug("Plotting %s linear triangles, %s quadratic triangles, "
                 "%s linear quads, %s 8-node quads, %s 9-node quads",
                 tri3_count, tri6_count, quad4_count, quad8_count, quad9_count)

    # Plot material zones first (if element_materials provided)
    if element_materials is not None:
        materials = np.unique(element_materials)

        # Material colors (style overrides → palette default). Mesh material IDs are
        # 1-based (gmsh); the style sheet / inputs key by 0-based mat_id, so map mat-1
        # — this also aligns the zone colors with the Inputs view.
        from .style import resolve_style, material_style
        _st = resolve_style(style)
        mat_to_color = {mat: material_style(_st, int(mat) - 1)["color"] for mat in materials}

        # Batch polygons by material for efficient rendering
        mat_fill_polys = {mat: [] for mat in materials}

        for idx, element_nodes in enumerate(elements):
            element_type = element_types[idx]
            mat = element_materials[idx]

            if element_type == 3:  # Linear triangle
                mat_fill_polys[mat].append(nodes[element_nodes[:3]])

            elif element_type == 6:  # Quadratic triangle - subdivide into 4 sub-triangles
                n0, n1, n2 = nodes[element_nodes[0]], nodes[element_nodes[1]], nodes[element_nodes[2]]
                n3, n4, n5 = nodes[element_nodes[3]], nodes[element_nodes[4]], nodes[element_nodes[5]]
                mat_fill_polys[mat].extend([
                    np.array([n0, n3, n5]),
                    np.array([n3, n1, n4]),
                    np.array([n5, n4, n2]),
                    np.array([n3, n4, n5]),
                ])

            elif element_type == 4:  # Linear quadrilateral
                mat_fill_polys[mat].append(nodes[element_nodes[:4]])

            elif element_type == 8:  # Quadratic quadrilateral - subdivide into 4 sub-quads
                n0, n1, n2, n3 = nodes[element_nodes[0]], nodes[element_nodes[1]], nodes[element_nodes[2]], nodes[element_nodes[3]]
                n4, n5, n6, n7 = nodes[element_nodes[4]], nodes[element_nodes[5]], nodes[element_nodes[6]], nodes[element_nodes[7]]
                center = np.array([(n0[0]+n1[0]+n2[0]+n3[0]+n4[0]+n5[0]+n6[0]+n7[0]) / 8,
                                   (n0[1]+n1[1]+n2[1]+n3[1]+n4[1]+n5[1]+n6[1]+n7[1]) / 8])
                mat_fill_polys[mat].extend([
                    np.array([n0, n4, center, n7]),
                    np.array([n4, n1, n5, center]),
                    np.array([center, n5, n2, n6]),
                    np.array([n7, center, n6, n3]),
                ])

            elif element_type == 9:  # 9-node quadrilateral
                n0, n1, n2, n3 = nodes[element_nodes[0]], nodes[element_nodes[1]], nodes[element_nodes[2]], nodes[element_nodes[3]]
                n4, n5, n6, n7 = nodes[element_nodes[4]], nodes[element_nodes[5]], nodes[element_nodes[6]], nodes[element_nodes[7]]
                center = nodes[element_nodes[8]]
                mat_fill_polys[mat].extend([
                    np.array([n0, n4, center, n7]),
                    np.array([n4, n1, n5, center]),
                    np.array([center, n5, n2, n6]),
                    np.array([n7, center, n6, n3]),
                ])

        # Render as batched PatchCollections (one per material)
        for mat in materials:
            polys = mat_fill_polys[mat]
            if polys:
                patch_list = [Polygon(p) for p in polys]
                pc = PatchCollection(patch_list, facecolor=mat_to_color[mat], edgecolor='none', alpha=alpha, gid='ZONE_FILL')
                ax.add_collection(pc)

    # Set up contour levels
    # Contour/colorbar range. Default (vmin/vmax None) auto-scales to THIS frame's data
    # — byte-identical to before. Passing an explicit range pins every panel of a
    # series to ONE scale (e.g. the full-pool t = 0 head range), so the colorbars match
    # across panels and a drawdown reads as one continuous story — late frames fade
    # toward uniform rather than re-normalizing into a bullseye. The frame's own data is
    # a subset of a series-wide range, so nothing is clipped.
    vmin = float(np.min(contour_data)) if vmin is None else float(vmin)
    vmax = float(np.max(contour_data)) if vmax is None else float(vmax)
    contour_levels = np.linspace(vmin, vmax, levels)

    # For contouring, subdivide tri6 elements into 4 subtriangles
    all_triangles_for_contouring = []
    for idx, element_nodes in enumerate(elements):
        element_type = element_types[idx]
        if element_type == 3:  # Linear triangular elements
            all_triangles_for_contouring.append(element_nodes[:3])
        elif element_type == 6:  # Quadratic triangular elements
            # Standard GMSH tri6 ordering: 3 = edge 0-1; 4 = edge 1-2; 5 = edge 2-0
            # Create 4 subtriangles: 0-3-5, 3-1-4, 5-4-2, 3-4-5
            subtriangles = [
                [element_nodes[0], element_nodes[3], element_nodes[5]],  # 0-3-5 (corner at 0)
                [element_nodes[3], element_nodes[1], element_nodes[4]],  # 3-1-4 (corner at 1)
                [element_nodes[5], element_nodes[4], element_nodes[2]],  # 5-4-2 (corner at 2)
                [element_nodes[3], element_nodes[4], element_nodes[5]]   # 3-4-5 (center)
            ]
            all_triangles_for_contouring.extend(subtriangles)
        elif element_type in [4, 8, 9]:  # Quadrilateral elements
            tri1 = [element_nodes[0], element_nodes[1], element_nodes[2]]
            tri2 = [element_nodes[0], element_nodes[2], element_nodes[3]]
            all_triangles_for_contouring.extend([tri1, tri2])
    triang = tri.Triangulation(nodes[:, 0], nodes[:, 1], all_triangles_for_contouring)

    # Variable labels for colorbar and title
    variable_labels = {
        "head": "Total Head",
        "u": "Pore Pressure",
        "v_mag": "Velocity Magnitude",
        "i_mag": "Hydraulic Gradient Magnitude"
    }
    variable_label = variable_labels[variable]

    # Declared unit labels (units plan phase 4), or None when the model declares no
    # system. Drives the flowrate-title unit (below) and the colorbar unit; None keeps
    # every string byte-identical to today. Head is a length, pore pressure a stress,
    # velocity magnitude the same dimension as k (length/time); hydraulic gradient is
    # dimensionless, so it takes no unit.
    from .plot import declared_unit_labels
    _unit_labels = declared_unit_labels(seep_data)
    _seep_var_unit_key = {"head": "length", "u": "stress", "v_mag": "k"}

    # Filled contours (only if fill_contours=True). The colorbar itself is added at
    # the very end (after _legend_below lays out the axes) so it can be sized to
    # cbar_shrink × the plot height and tracks the box-adjusted plot box.
    contourf = None
    if fill_contours:
        # Fill opacity drives the contour-fill wash (0 = clean line flow net, the
        # default; raise it for a colored head field). The colorbar still shows.
        contourf = ax.tricontourf(triang, contour_data, levels=contour_levels, cmap=cmap, vmin=vmin, vmax=vmax, alpha=alpha)
        contourf.set_gid('CONTOUR_FILL')

    # Solid lines for contours. linestyles is pinned solid so the equipotentials
    # read as a crisp flow net regardless of sign (matplotlib would otherwise dash
    # any negative level via rcParam contour.negative_linestyle — e.g. for signed
    # variables like pore pressure).
    _cs = ax.tricontour(triang, contour_data, levels=contour_levels, colors="k",
                        linewidths=0.5, linestyles="solid")
    _cs.set_gid('HEAD_CONTOURS')

    # Phreatic surface (pressure head = 0)
    # Only meaningful for UNCONFINED flow. A confined solve is a single saturated
    # Laplace solve in which the head is a potential, not a water level: negative
    # pore pressure is routine (any node above the boundary head has u < 0) and the
    # p = 0 contour is an artifact, not a water table. Drawing it there produces a
    # figure that contradicts its own flow net -- flow lines correctly fill the whole
    # saturated domain while the bogus "phreatic surface" implies most of it is dry.
    #
    # The condition is flownet_has_phreatic, so a caller asking what this figure
    # carries gets the answer from the function that decides it.
    has_phreatic = False
    if phreatic and flownet_has_phreatic(seep_data, solution):
        elevation = nodes[:, 1]  # y-coordinate is elevation
        pressure_head = head - elevation
        _csp = ax.tricontour(triang, pressure_head, levels=[0], colors="black", linewidths=2.0)
        _csp.set_gid('PHREATIC')
        has_phreatic = True

    # A stream function (and therefore flow lines) exists only for divergence-free
    # flow. Under transient storage exchange a pure storage-release frame (boundary
    # inflow = 0, no through-flow) has a flat phi, for which np.linspace collapses to
    # equal contour levels and tricontour raises "Contour levels must be increasing".
    # That is PHYSICS, not an error: no through-flow -> no flow net. Skip the flow-line
    # overlay for such a frame so the shared renderer never crashes (which, swallowed
    # by a caller's try/except, would freeze an animation on the last good frame). A
    # steady solve always has a head drop, so phi has range and this never triggers —
    # the steady figure is unchanged.
    #
    # The test, the flowrate the lines are spaced by and the conductivities the
    # spacing is computed from are all flownet_has_flowlines, so a caller asking what
    # this figure carries gets the answer from the function that decides it. The one
    # answer also drives the legend and the transient subtitle below: a "Flow line"
    # key on a figure that drew none is the same false statement in a smaller font.
    has_flowlines = flownet_has_flowlines(seep_data, solution)
    if plot_flowlines and has_flowlines:
        # Compute head drop for flowline calculation
        hdrop = vmax - vmin
        if base_mat > len(k1_by_mat):
            print(f"Warning: base_mat={base_mat} is larger than number of materials ({len(k1_by_mat)}). Using material 1.")
            base_mat = 1
        elif base_mat < 1:
            print(f"Warning: base_mat={base_mat} is less than 1. Using material 1.")
            base_mat = 1
        base_k1 = k1_by_mat[base_mat - 1]
        base_k2 = k2_by_mat[base_mat - 1] if k2_by_mat is not None else base_k1
        # For anisotropic media, the equivalent conductivity is sqrt(k1*k2)
        # This ensures the flow net cells have the correct aspect ratio sqrt(k1/k2)
        base_k = np.sqrt(base_k1 * base_k2)
        ne = levels - 1
        nf = (flowrate * ne) / (base_k * hdrop)
        phi_levels = max(round(nf) + 1, 2)
        logger.debug("Computed nf: %.2f, using %s φ contours (flowrate=%.3f, base k=%.4f, head drop=%.3f)",
                     nf, phi_levels, flowrate, base_k, hdrop)
        phi_contours = np.linspace(np.min(phi), np.max(phi), phi_levels)
        _csf = ax.tricontour(triang, phi, levels=phi_contours, colors="blue", linewidths=0.7, linestyles="solid")
        _csf.set_gid('FLOWLINES')

    # Plot velocity vectors if requested
    if vectors:
        velocity = solution.get("velocity")
        if velocity is not None:
            # Calculate x_range for scaling
            x_min_vec = nodes[:, 0].min()
            x_max_vec = nodes[:, 0].max()
            x_range = x_max_vec - x_min_vec
            max_vector_length = x_range * vector_scale
            
            # Get velocity magnitude
            v_mag = solution.get("v_mag")
            if v_mag is None:
                # Calculate v_mag if not available
                v_mag = np.linalg.norm(velocity, axis=1)
            
            # Find maximum velocity magnitude
            max_v_mag = np.max(v_mag)
            
            # Scale vectors: if max_v_mag > 0, scale so max vector has length max_vector_length
            if max_v_mag > 0:
                scale_factor = max_vector_length / max_v_mag
                velocity_scaled = velocity * scale_factor
            else:
                velocity_scaled = velocity
            
            # Plot vectors using quiver
            _q = ax.quiver(nodes[:, 0], nodes[:, 1], velocity_scaled[:, 0], velocity_scaled[:, 1],
                     angles='xy', scale_units='xy', scale=1, width=0.002, headwidth=2.5,
                     headlength=3, headaxislength=2.5, color='black', alpha=0.7)
            _q.set_gid('VELOCITY')

    # Plot element edges if requested
    if mesh:
        mesh_segments = []
        for element, elem_type in zip(elements, element_types if element_types is not None else [3]*len(elements)):
            if elem_type in (3, 6):
                # Triangle: corner edges 0-1, 1-2, 2-0
                mesh_segments.extend([
                    [nodes[element[0]], nodes[element[1]]],
                    [nodes[element[1]], nodes[element[2]]],
                    [nodes[element[2]], nodes[element[0]]],
                ])
            elif elem_type in (4, 8, 9):
                # Quad: corner edges 0-1, 1-2, 2-3, 3-0
                mesh_segments.extend([
                    [nodes[element[0]], nodes[element[1]]],
                    [nodes[element[1]], nodes[element[2]]],
                    [nodes[element[2]], nodes[element[3]]],
                    [nodes[element[3]], nodes[element[0]]],
                ])
        if mesh_segments:
            lc = LineCollection(mesh_segments, colors='darkgray', linewidths=0.5, alpha=0.7, gid='MESH')
            ax.add_collection(lc)

    # Plot the mesh boundary
    try:
        boundary = get_ordered_mesh_boundary(nodes, elements, element_types)
        ax.plot(boundary[:, 0], boundary[:, 1], color="black", linewidth=1.0, label="Mesh Boundary", gid='MESH_BOUNDARY')
    except Exception as e:
        print(f"Warning: Could not plot mesh boundary: {e}")

    # Add cushion around the mesh
    x_min, x_max = nodes[:, 0].min(), nodes[:, 0].max()
    y_min, y_max = nodes[:, 1].min(), nodes[:, 1].max()
    x_pad = (x_max - x_min) * pad_frac
    y_pad = (y_max - y_min) * pad_frac
    ax.set_xlim(x_min - x_pad, x_max + x_pad)
    ax.set_ylim(y_min - y_pad, y_max + y_pad)

    # Opt-in BC water levels: each head/reservoir boundary's instantaneous level for
    # this frame (the pool drops through a transient playback). Drawn after the limits
    # so the symbol can be sized to the domain; its legend handles are collected here
    # and appended below. Default off keeps every committed steady figure byte-stable.
    bc_level_handles = []
    if show_bc_levels:
        try:
            bc_level_handles = _draw_seep_bc_levels(ax, seep_data, solution, style)
        except Exception as e:
            print(f"Warning: Could not draw BC water levels: {e}")

    # Build title based on variable
    def _q_fmt(v):
        # Fixed 3-decimal for O(0.1)+ flowrates (keeps the mesh-vs-exact deviation
        # legible, e.g. 5.029); 3 significant figures below so tiny flows never
        # collapse to "0.000".
        return f"{v:.3f}" if abs(v) >= 0.1 else f"{v:.3g}"
    # Append the flowrate unit ("m³/day per m") only when the model declares a unit
    # system AND a time unit; undeclared/time-less models stay unlabeled (byte-
    # identical to today).
    q_unit = f" {_unit_labels['flowrate']}" if (_unit_labels and _unit_labels.get("flowrate")) else ""
    # Transient frames carry separate boundary inflow/outflow (they differ under
    # storage exchange, so a single "Total Flowrate" would misrepresent the frame)
    # and a frame time. These keys are absent on steady solutions, so the steady
    # title is byte-identical.
    inflow = solution.get("inflow")
    outflow = solution.get("outflow")
    frame_time = solution.get("time")
    # A transient frame carries separate boundary inflow/outflow (both keys absent on
    # a steady solution). Its title is a compact two-line variant: a short main line
    # (identity + frame time — NOT a narration of the display options, which the
    # legend already names) plus a smaller second line carrying the per-frame
    # mass-balance numbers. Flow lines do NOT apply to a transient (storage-release)
    # solution — a flow net requires divergence-free through-flow, which a draining
    # frame is not — so Studio does not request them and velocity vectors read the
    # instantaneous flow direction instead. The honest "no through-flow" note is kept
    # ONLY for a direct caller that explicitly requested flow lines on such a frame
    # (phi degenerate): it explains why none were drawn. A steady solution has neither
    # inflow/outflow key, so it takes the byte-identical legacy branch below.
    subtitle = None
    if inflow is not None and outflow is not None:
        title = "Seepage Solution"
        if frame_time is not None:
            t_unit = f" {_unit_labels['time']}" if (_unit_labels and _unit_labels.get("time")) else ""
            title += f" — t = {frame_time:g}{t_unit}"
        if plot_flowlines and not has_flowlines:
            subtitle = "no through-flow — flow lines undefined"
        else:
            subtitle = f"Inflow {_q_fmt(inflow)} / Outflow {_q_fmt(outflow)}{q_unit}"
    else:
        # Steady solution: a compact two-line title mirroring the transient variant
        # above — a fixed "Seepage Solution" main line plus a smaller second line for
        # the total flowrate. The flow-net elements the old title used to narrate
        # ("Flow Net: Total Head Contours and Flowlines with Phreatic Surface", the
        # confined variant simply dropping the phreatic clause) are now named by the
        # legend, so the main line carries only the identity. A non-head field (pore
        # pressure, velocity, gradient) is identified by the colorbar, so it takes the
        # bare main line with no second line.
        title = "Seepage Solution"
        if variable == "head" and flowrate is not None:
            subtitle = f"Total Flowrate: {_q_fmt(flowrate)}{q_unit}"
        if frame_time is not None:
            t_unit = f" {_unit_labels['time']}" if (_unit_labels and _unit_labels.get("time")) else ""
            title += f" — t = {frame_time:g}{t_unit}"
    if show_title:
        if subtitle:
            # Two-line title: main line at the normal size, second line ~0.8x beneath
            # it. Raise the title with extra pad so the smaller line sits between it
            # and the axes (top→bottom: main, subtitle, plot) and tight_layout reserves
            # the room — so nothing clips at the (narrow) Studio viewport width.
            ax.set_title(title)
            base_fs = ax.title.get_fontsize()
            ax.set_title(title, pad=base_fs * 1.7 + 6.0)
            ax.annotate(subtitle, xy=(0.5, 1.0), xycoords="axes fraction",
                        xytext=(0.0, 4.0), textcoords="offset points",
                        ha="center", va="bottom", fontsize=base_fs * 0.82,
                        annotation_clip=False)
        else:
            ax.set_title(title)

    # Equal aspect (1:1) so the flow net reads as curvilinear squares. Box-adjust
    # (not datalim) because seepage domains are characteristically wide and thin
    # (blankets, dams) — datalim would expand the y-limits to fill the figure and
    # squish the data into a tiny band; box-adjust fits it snugly as a wide strip.
    ax.set_aspect("equal")

    # Legend (below the axes, frameless by default) for the flow-net features that
    # were actually drawn — material zones plus phreatic / contour / flow lines —
    # alongside the colorbar. Swatches use a fixed readable alpha so the material
    # color key stays visible even when the fill opacity is 0 (clean flow net).
    from .plot import _legend_below
    mat_names = seep_data.get("material_names", [])
    leg_handles = []
    # Material swatches belong in the legend only when the material zone color is
    # what the fill actually shows (fill_contours=False). With fill_contours=True the
    # fill IS the variable field (its key is the colorbar), so a material-colored
    # swatch would name a color the eye can't find on the plot — the legend must
    # describe what is drawn.
    if element_materials is not None and not fill_contours:
        for mat in materials:
            nm = mat_names[mat - 1] if (mat_names and mat <= len(mat_names)) else f"Material {mat}"
            leg_handles.append(mpatches.Patch(facecolor=mat_to_color[mat], alpha=0.6,
                                              edgecolor="none", label=nm))
    if has_phreatic:
        leg_handles.append(plt.Line2D([0], [0], color="black", lw=2.0, label="Phreatic surface"))
    # Sentence case (capitalize() → "Total head contour") to match "Flow line".
    leg_handles.append(plt.Line2D([0], [0], color="black", lw=0.5,
                                  label=f"{variable_label.capitalize()} contour"))
    if plot_flowlines and has_flowlines:
        leg_handles.append(plt.Line2D([0], [0], color="blue", lw=0.7, label="Flow line"))
    if vectors:
        leg_handles.append(plt.Line2D([0], [0], color="black", lw=0, marker=r"$\rightarrow$",
                                      markersize=10, label="Velocity"))
    # BC water-level entries (one per drawn type), when the overlay is on.
    leg_handles.extend(bc_level_handles)
    # Tighten the left/right margins so the (wide-thin) domain fills the width like
    # the data/inputs plots. Safe now that the colorbar is deferred to after the
    # legend: tight_layout sees no gridspec colorbar to choke on. Reserve right-side
    # room only when a colorbar will actually be drawn.
    try:
        fig.tight_layout()
        if fill_contours:
            fig.subplots_adjust(right=min(fig.subplotpars.right, 0.90))
    except Exception:
        pass
    _legend_below(ax, fig, handles=leg_handles, legend_ncol=legend_ncol, frameon=legend_frame, show_legend=show_legend)

    # Colorbar last: a manual axes to the right of the (now laid-out) plot, its
    # height cbar_shrink × the plot height and centered on it. Manual placement (vs
    # a gridspec/divider colorbar) needs no layout engine — which would clash with
    # the manual margins above — and honors the "Colorbar size" control.
    if contourf is not None:
        try:
            fig.canvas.draw()
            ax.apply_aspect()
        except Exception:
            pass
        pos = ax.get_position()
        shrink = min(1.0, max(0.1, cbar_shrink))
        ch = pos.height * shrink
        cy = pos.y0 + (pos.height - ch) / 2.0
        cax = fig.add_axes([pos.x1 + 0.012, cy, 0.014, ch])
        # Build the colorbar from a full-opacity mappable so it stays a solid color
        # key even when the fill wash is faint/off (fill opacity 0 = clean flow net).
        from matplotlib.cm import ScalarMappable
        from matplotlib.colors import Normalize
        sm = ScalarMappable(norm=Normalize(vmin=vmin, vmax=vmax), cmap=cmap)
        # Append the field's unit to the colorbar label only (not the legend/title
        # copies of variable_label), and only when declared; undeclared = today's text.
        _cb_key = _seep_var_unit_key.get(variable)
        _cb_unit = _unit_labels.get(_cb_key) if (_unit_labels and _cb_key) else None
        cbar_label = f"{variable_label} ({_cb_unit})" if _cb_unit else variable_label
        cbar = fig.colorbar(sm, cax=cax, label=cbar_label)
        # Tick density adapts to the colorbar's drawn height so the labels never
        # stack (a short/wide domain thins to a few readable ticks; a tall bar keeps
        # finer labeling). This rule — originally added here for the sheetpile flow
        # net — now lives in the shared adaptive_colorbar_ticks helper so the FEM
        # result plots reuse it; behaviour is unchanged for these flow nets.
        from .plot import adaptive_colorbar_ticks
        adaptive_colorbar_ticks(fig, cbar)


    base_name = 'plot_' + title.lower().replace(' ', '_').replace(':', '').replace(',', '').replace('—', '').replace('(', '').replace(')', '')
    if save_png:
        fig.savefig(base_name + '.png', dpi=dpi, bbox_inches='tight')
    if save_dxf:
        from .cad import axes_to_dxf
        axes_to_dxf(ax, base_name + '.dxf')

    if own_fig:
        plt.show()
    return fig


    # plot_seep_material_table has been moved to xslope/plot.py


def get_ordered_mesh_boundary(nodes, elements, element_types=None):
    """
    Extracts the outer boundary of the mesh and returns it as an ordered array of points.
    Supports both triangular and quadrilateral elements.

    Returns:
        np.ndarray of shape (N, 2): boundary coordinates in order (closed loop)
    """
    import numpy as np
    from collections import defaultdict, deque

    # If element_types is not provided, assume all triangles (backward compatibility)
    if element_types is None:
        element_types = np.full(len(elements), 3)

    # Step 1: Count all edges
    edge_count = defaultdict(int)
    edge_to_nodes = {}

    for i, element_nodes in enumerate(elements):
        element_type = element_types[i]
        
        if element_type == 3:
            # Triangle: 3 edges
            for j in range(3):
                a, b = sorted((element_nodes[j], element_nodes[(j + 1) % 3]))
                edge_count[(a, b)] += 1
                edge_to_nodes[(a, b)] = (element_nodes[j], element_nodes[(j + 1) % 3])  # preserve direction
        elif element_type == 4:
            # Quadrilateral: 4 edges
            for j in range(4):
                a, b = sorted((element_nodes[j], element_nodes[(j + 1) % 4]))
                edge_count[(a, b)] += 1
                edge_to_nodes[(a, b)] = (element_nodes[j], element_nodes[(j + 1) % 4])  # preserve direction
        elif element_type == 6:
            # 6-node triangle: 3 edges (use only corner nodes 0,1,2)
            for j in range(3):
                a, b = sorted((element_nodes[j], element_nodes[(j + 1) % 3]))
                edge_count[(a, b)] += 1
                edge_to_nodes[(a, b)] = (element_nodes[j], element_nodes[(j + 1) % 3])  # preserve direction
        elif element_type in [8, 9]:
            # Higher-order quadrilaterals: 4 edges (use only corner nodes 0,1,2,3)
            for j in range(4):
                a, b = sorted((element_nodes[j], element_nodes[(j + 1) % 4]))
                edge_count[(a, b)] += 1
                edge_to_nodes[(a, b)] = (element_nodes[j], element_nodes[(j + 1) % 4])  # preserve direction

    # Step 2: Keep only boundary edges (appear once)
    boundary_edges = [edge_to_nodes[e] for e, count in edge_count.items() if count == 1]

    if not boundary_edges:
        raise ValueError("No boundary edges found.")

    # Step 3: Build adjacency for boundary walk
    adj = defaultdict(list)
    for a, b in boundary_edges:
        adj[a].append(b)
        adj[b].append(a)

    # Step 4: Walk all boundary segments
    all_boundary_nodes = []
    remaining_edges = set(boundary_edges)
    
    while remaining_edges:
        # Start a new boundary segment
        start_edge = remaining_edges.pop()
        start_node = start_edge[0]
        current_node = start_edge[1]
        
        segment = [start_node, current_node]
        remaining_edges.discard((current_node, start_node))  # Remove reverse edge if present
        
        # Walk this segment until we can't continue
        while True:
            # Find next edge from current node
            next_edge = None
            for edge in remaining_edges:
                if edge[0] == current_node:
                    next_edge = edge
                    break
                elif edge[1] == current_node:
                    next_edge = (edge[1], edge[0])  # Reverse the edge
                    break
            
            if next_edge is None:
                break
                
            next_node = next_edge[1]
            segment.append(next_node)
            remaining_edges.discard(next_edge)
            remaining_edges.discard((next_node, current_node))  # Remove reverse edge if present
            current_node = next_node
            
            # Check if we've closed the loop
            if current_node == start_node:
                break
        
        all_boundary_nodes.extend(segment)
    
    # If we have multiple segments, we need to handle them properly
    # For now, just return the first complete segment
    if all_boundary_nodes:
        # Ensure the boundary is closed
        if all_boundary_nodes[0] != all_boundary_nodes[-1]:
            all_boundary_nodes.append(all_boundary_nodes[0])
        return nodes[all_boundary_nodes]
    else:
        raise ValueError("No boundary nodes found.")