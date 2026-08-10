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

import warnings

import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Polygon

from . import colormaps as _colormaps  # noqa: F401  (registers the BGYR ramp by name)
from .plot import adaptive_colorbar_ticks, declared_unit_labels


def _fem_cbar_label(fem_data, base, unit_key):
    """``"<base> (<unit>)"`` when ``fem_data`` declares a unit system and supplies a
    non-empty label for ``unit_key`` (e.g. ``"stress"``, ``"force_per_len"``,
    ``"length"``); the bare ``base`` otherwise.

    Units plan phase 4: labels the FEM stress / force / displacement colorbars where a
    system is declared. Undeclared models get the exact same string as today, so those
    figures stay pixel-identical."""
    labels = declared_unit_labels(fem_data)
    if labels:
        unit = labels.get(unit_key)
        if unit:
            return f"{base} ({unit})"
    return base

# Median rendered element edge (device px) below which two full interleaved grids
# (original + deformed) tangle; below it the original mesh collapses to its domain
# boundary outline so the deformed grid alone carries the deformation.
_DENSE_TANGLE_EDGE_PX = 9.0


def _fs_title(base, F, fs=None, at_failure=False):
    """Result-panel title that keeps the SSRM factor of safety and the rendered
    viscoplastic state unambiguous.

    Normal (converged) panels are rendered at the LAST CONVERGED strength-reduction
    factor (``solution['F']``); the reported factor of safety is the SSRM bracket
    midpoint (``result['FS']``), which at a wide tolerance can round differently. When
    both are known and differ at two-decimal display, the title names both; otherwise
    it keeps the simple ``F = X.XX`` form (or just ``base`` when no F is available).

    ``at_failure`` marks a panel rendering the UNCONVERGED at-failure field (captured
    a margin beyond critical). It leads with the factor of safety and stops there —
    "at Failure" already discloses the field is the unconverged trial state, so a
    parenthetical naming the trial F would only repeat that disclosure as noise.
    """
    if at_failure and fs is not None:
        return f"{base}  FS = {fs:.2f}"
    if F is None:
        return base
    if fs is None or f"{fs:.2f}" == f"{F:.2f}":
        return f"{base}  F={F:.2f}"
    return f"{base}  FS = {fs:.2f} (rendered at last converged F = {F:.2f})"


#: What the viscoplastic maximum shear strain field is called — the name the
#: documentation uses for it, and therefore the one the panel caption, the
#: sentence citing it, the colorbar and the results view's own list all use. It
#: had four: "Maximum shear strain" (which is the name of a DIFFERENT column,
#: the total-strain one), "vp_shear_strain", "VP Max Shear Strain" and "Shear
#: strain".
SHEAR_STRAIN_LABEL = "Viscoplastic shear strain"


def _vector_cbar_label(disp_elastic):
    """The displacement-vector colorbar's label, naming the field it colors: the
    viscoplastic magnitude where the elastic part is known and can be taken out,
    the total magnitude otherwise. It read "VP Displacement Magnitude" either
    way."""
    return ("Viscoplastic displacement magnitude, |u|" if disp_elastic is not None
            else "Displacement magnitude, |u|")


def deformation_scale(fem_data, field, deform_percent=15):
    """The exaggeration the deformed-mesh panel is drawn at.

    The displacements a slope settles under gravity are invisible at true scale,
    so the deformed grid is drawn at a multiplier chosen to make the largest one
    ``deform_percent`` of the mesh height. That multiplier is a fact about the
    figure a reader needs — a shape drawn at 54x is not the shape of the slope —
    and it is derived HERE so the number the report prints and the number the
    figure is drawn at cannot be two numbers.

    ``field`` is the solve_fem field the panel renders (the at-failure snapshot
    where one is being drawn), and the scale is measured on its viscoplastic
    displacement where the elastic part is known, matching plot_deformed_mesh.
    """
    nodes = fem_data["nodes"]
    disp = field.get("displacements", np.zeros(2 * len(nodes)))
    disp_elastic = field.get("displacements_elastic", None)
    if disp_elastic is not None:
        disp = disp - disp_elastic
    u_arr, v_arr = _extract_uv(disp, fem_data)
    max_disp = np.max(np.sqrt(u_arr**2 + v_arr**2))
    mesh_height = np.max(nodes[:, 1]) - np.min(nodes[:, 1])
    if max_disp <= 1e-30:
        return 1.0
    return max(1.0, (mesh_height * deform_percent / 100) / max_disp)


def displacement_magnitude(fem_data, field):
    """Nodal viscoplastic displacement magnitudes for one solve_fem field.

    The quantity the displacement-vector panel draws an arrow of, and the one the
    deformed grid is scaled by — read here once so a caller pinning two panels to
    a shared scale measures what those panels measure.
    """
    nodes = fem_data["nodes"]
    disp = (field or {}).get("displacements", np.zeros(2 * len(nodes)))
    disp_elastic = (field or {}).get("displacements_elastic", None)
    if disp_elastic is not None:
        disp = disp - disp_elastic
    u_arr, v_arr = _extract_uv(disp, fem_data)
    return np.sqrt(u_arr ** 2 + v_arr ** 2)


def shear_strain_field(fem_data, field):
    """The per-element array the viscoplastic shear strain panel contours, with
    the same fallback that panel makes — total shear strain where a solution
    carries no viscoplastic one. ``None`` where neither is there."""
    values = (field or {}).get("vp_shear_strain", None)
    if values is not None:
        return np.asarray(values, dtype=float)
    strains = (field or {}).get("strains", None)
    if strains is not None and np.asarray(strains).shape[1:] and \
            np.asarray(strains).shape[1] >= 4:
        return np.asarray(strains, dtype=float)[:, 3]
    return None


def shared_panel_scales(fem_data, fields, deform_percent=15):
    """One color range, one exaggeration and one arrow scale for a set of fields.

    A section drawn at two states — the mechanism at failure and the last
    converged trial — is two pictures of ONE run, and each auto-scaled to its own
    field is two pictures of two. Held together, a strain color means the same
    strain in both, the deformed grids are exaggerated by the same multiplier, and
    an arrow of a given length means the same displacement: the difference
    between the pictures is then the difference between the states. This is the
    rule the paired seepage sets are drawn under, applied to the field states.

    Returns ``{"vmin", "vmax", "deform_scale", "vector_max"}``, each ``None``
    where the set gives nothing to hold: a variable one field does not carry, or
    one flat across the set, leaves the panels that draw it to scale themselves,
    which is what a run drawn at one state wants.

    The exaggeration is the SMALLER of the two the fields ask for — drawn at the
    larger, the field that moved further would leave the section.
    """
    fields = [f for f in (fields or []) if f]
    out = {"vmin": None, "vmax": None, "deform_scale": None, "vector_max": None}
    if len(fields) < 2:
        return out

    lo, hi = [], []
    for field in fields:
        values = shear_strain_field(fem_data, field)
        if values is None:
            lo, hi = [], []
            break
        nodal = _nodal_average(fem_data, values)
        finite = nodal[np.isfinite(nodal)]
        if not finite.size:
            lo, hi = [], []
            break
        lo.append(float(finite.min()))
        hi.append(float(finite.max()))
    if lo and max(hi) > min(lo):
        out["vmin"], out["vmax"] = min(lo), max(hi)

    scales = [deformation_scale(fem_data, field, deform_percent)
              for field in fields]
    if scales:
        out["deform_scale"] = min(scales)

    peaks = [float(displacement_magnitude(fem_data, field).max())
             for field in fields]
    if peaks and max(peaks) > 0:
        out["vector_max"] = max(peaks)
    return out


def _place_deform_legend(ax, show_legend=True):
    """Draw the deformed-mesh legend (Original / Deformed, plus any reinforcement
    entries) INSIDE the deformation axes, in the empty corner above the slope
    profile. A compact semi-opaque frame keeps it legible over the mesh while
    consuming zero layout height — replacing the old full-width legend band that
    ate a strip between panels."""
    if not show_legend:
        return
    handles, labels = ax.get_legend_handles_labels()
    if not handles:
        return
    leg = ax.legend(handles, labels, loc='upper right', fontsize=9, frameon=True,
                    framealpha=0.9, borderpad=0.4, handlelength=1.6, labelspacing=0.3)
    # The plotted grids may be sub-point hairlines on a dense mesh; give the legend
    # swatches a fixed readable weight so Original/Deformed stay identifiable.
    for lh in getattr(leg, "legend_handles", getattr(leg, "legendHandles", [])):
        try:
            lh.set_linewidth(2.0)
        except Exception:
            pass


def _place_stacked_cbars(fig, ax, specs):
    """Stack full-height colorbars to the right of ``ax`` via make_axes_locatable —
    innermost first — each in its own appended slot so their ticks and rotated labels
    never collide (the bug when a field colorbar and a reinforcement-force colorbar
    shared one slot). ``specs`` is a list of (mappable, label); a None mappable
    reserves an invisible slot of the same width (used to keep sibling panels
    x-aligned). The first slot hugs the axes; later ones get a wide pad to clear the
    previous bar's tick+axis labels. Returns the list of Colorbar objects (None where
    a slot was reserved but empty). Self-adjusting — no overlap at any figure aspect."""
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    div = make_axes_locatable(ax)
    cbars = []
    for k, (mappable, label) in enumerate(specs):
        pad = 0.15 if k == 0 else 0.75
        cax = div.append_axes("right", size="3%", pad=pad)
        if mappable is None:
            cax.set_axis_off()
            cbars.append(None)
            continue
        cb = fig.colorbar(mappable, cax=cax)
        cb.set_label(label, rotation=270, labelpad=15)
        cbars.append(cb)
    return cbars


def _extract_uv(disp, fem_data):
    """Extract per-node u,v displacements from a mixed-DOF displacement vector."""
    dof_offset = fem_data.get("dof_offset", None)
    if dof_offset is not None:
        n_nodes = len(fem_data["nodes"])
        u = np.array([disp[dof_offset[i]] for i in range(n_nodes)])
        v = np.array([disp[dof_offset[i] + 1] for i in range(n_nodes)])
    else:
        u = disp[0::2]
        v = disp[1::2]
    return u, v


def plot_fem_data(fem_data, figsize=(12, 7), show_nodes=False, show_bc=True,
                  label_elements=False, label_nodes=False, alpha=0.6, bc_symbol_size=0.03, save_png=False, save_dxf=False, dpi=300, legend_ncol="auto", legend_frame=False, show_title=True, show_legend=True, fig=None, style=None):
    """
    Plots a FEM mesh colored by material zone with boundary conditions displayed.

    Args:
        fem_data: Dictionary containing FEM data from build_fem_data
        figsize: Figure size
        show_nodes: If True, plot node points
        show_bc: If True, plot boundary condition symbols
        label_elements: If True, label each element with its number at its centroid
        label_nodes: If True, label each node with its number just above and to the right
        alpha: Transparency for element faces
        bc_symbol_size: Size factor for boundary condition symbols (as fraction of mesh size)
    """
    from matplotlib.collections import PatchCollection

    # Extract data from fem_data
    nodes = fem_data["nodes"]
    elements = fem_data["elements"]
    element_materials = fem_data["element_materials"]
    element_types = fem_data.get("element_types", None)
    bc_type = fem_data["bc_type"]
    bc_values = fem_data["bc_values"]

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
    from .style import resolve_style, material_style
    _st = resolve_style(style)
    mat_to_color = {mat: material_style(_st, int(mat) - 1)["color"] for mat in materials}

    # If element_types is not provided, assume all triangles (backward compatibility)
    if element_types is None:
        element_types = np.full(len(elements), 3)

    # Batch polygons and edge lines by material for efficient rendering
    # Key: material -> list of polygon vertex arrays
    mat_fill_polys = {mat: [] for mat in materials}
    edge_segments = []  # for outer boundaries of quadratic elements

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
            # Outer boundary edges
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
            # All linear elements — draw with edges
            patch_list = [Polygon(p) for p in polys]
            pc = PatchCollection(patch_list, facecolor=color, edgecolor='k', linewidth=0.5, alpha=alpha, gid='MESH_FILL')
            ax.add_collection(pc)
        elif has_no_edge and not has_edge:
            # All quadratic sub-polys — no edges on fills
            patch_list = [Polygon(p) for p in polys]
            pc = PatchCollection(patch_list, facecolor=color, edgecolor='none', alpha=alpha, gid='MESH_FILL')
            ax.add_collection(pc)
        else:
            # Mixed — separate linear (with edges) and quadratic sub-polys (no edges)
            linear_polys = []
            sub_polys = []
            sub_idx = 0
            for i in range(len(elements)):
                if element_materials[i] != mat:
                    continue
                et = element_types[i]
                if et == 3:
                    linear_polys.append(polys[sub_idx]); sub_idx += 1
                elif et == 4:
                    linear_polys.append(polys[sub_idx]); sub_idx += 1
                elif et == 6:
                    sub_polys.extend(polys[sub_idx:sub_idx+4]); sub_idx += 4
                elif et in (8, 9):
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
        ax.plot(nodes[:, 0], nodes[:, 1], 'k.', markersize=2, gid='MESH_NODES')

    # Label node numbers if requested
    if label_nodes:
        for i, (x, y) in enumerate(nodes):
            ax.text(x + 0.5, y + 0.5, str(i+1), fontsize=6, color='blue', alpha=0.7,
                    ha='left', va='bottom', zorder=11)

    # Get material names if available
    material_names = fem_data.get("material_names", [])

    legend_handles = []
    for mat in materials:
        if material_names and mat <= len(material_names):
            label = material_names[mat - 1]
        else:
            label = f"Material {mat}"
        legend_handles.append(
            patches.Patch(facecolor=mat_to_color[mat], alpha=alpha,
                          edgecolor="none", label=label)
        )

    # Plot 1D elements (reinforcement truss + pile beam) using LineCollection
    elements_1d = fem_data.get("elements_1d", np.array([]).reshape(0, 3))
    pile_elem_mask = fem_data.get("pile_elem_mask", np.zeros(len(elements_1d), dtype=bool))
    n_reinf_plotted = 0
    n_pile_plotted = 0
    if len(elements_1d) > 0:
        reinf_segs = []
        pile_segs = []
        for elem_idx in range(len(elements_1d)):
            elem_nodes_1d = elements_1d[elem_idx]
            seg = [nodes[elem_nodes_1d[0]], nodes[elem_nodes_1d[1]]]
            if pile_elem_mask[elem_idx]:
                pile_segs.append(seg)
                n_pile_plotted += 1
            else:
                reinf_segs.append(seg)
                n_reinf_plotted += 1
        if reinf_segs:
            lc = LineCollection(reinf_segs, colors='red', linewidths=2.5, zorder=5, gid='REINFORCEMENT')
            ax.add_collection(lc)
            legend_handles.append(
                plt.Line2D([0], [0], color='red', lw=2.5, label=f'Reinforcement ({n_reinf_plotted} elements)')
            )
        if pile_segs:
            lc = LineCollection(pile_segs, colors='green', linewidths=3.5, zorder=5, gid='PILES')
            ax.add_collection(lc)
            legend_handles.append(
                plt.Line2D([0], [0], color='green', lw=3.5, label=f'Pile ({n_pile_plotted} elements)')
            )

    # Plot boundary conditions
    if show_bc:
        saved_roller_x = fem_data.get("roller_x_nodes", set())
        _plot_boundary_conditions(ax, nodes, bc_type, bc_values, legend_handles, bc_symbol_size, saved_roller_x)

    from .plot import _legend_below
    # Adjust plot limits to accommodate force arrows
    x_min, x_max = nodes[:, 0].min(), nodes[:, 0].max()
    y_min, y_max = nodes[:, 1].min(), nodes[:, 1].max()
    
    # Add extra space for force arrows if they exist
    force_nodes = np.where(bc_type == 4)[0]
    if len(force_nodes) > 0:
        # Find the extent of force arrows
        mesh_size = min(x_max - x_min, y_max - y_min)
        symbol_size = mesh_size * bc_symbol_size
        
        # Add padding for force arrows (they extend outward from nodes)
        y_padding = symbol_size * 4  # Extra space above for upward arrows
        x_padding = (x_max - x_min) * 0.05  # Standard padding
        y_padding_bottom = (y_max - y_min) * 0.05
    else:
        # Standard padding
        x_padding = (x_max - x_min) * 0.05
        y_padding = (y_max - y_min) * 0.05
        y_padding_bottom = y_padding
    
    ax.set_xlim(x_min - x_padding, x_max + x_padding)
    ax.set_ylim(y_min - y_padding_bottom, y_max + y_padding)
    # Box-adjust (the default) keeps the requested x/y limits and shrinks the axes
    # box to a snug wide strip — matching plot_seep_data and the FEM result plots.
    # (adjustable="datalim" would instead expand the data range to fill the axes,
    # which overrides the limits set above and makes matplotlib log a "Ignoring
    # fixed limits…" warning on every redraw.)
    ax.set_aspect("equal")
    
    # Count element types for title
    num_tri = np.sum((element_types == 3) | (element_types == 6))
    num_quad = np.sum((element_types == 4) | (element_types == 8) | (element_types == 9))
    num_1d = len(elements_1d)
    parts = []
    if num_tri > 0:
        parts.append(f"{num_tri} triangles")
    if num_quad > 0:
        parts.append(f"{num_quad} quads")
    if n_reinf_plotted > 0:
        parts.append(f"{n_reinf_plotted} reinforcement")
    if n_pile_plotted > 0:
        parts.append(f"{n_pile_plotted} pile")
    title = f"FEM Mesh with Material Zones ({', '.join(parts)})"
    
    if show_title:
        ax.set_title(title)
    fig.tight_layout()
    # Combined legend below the plot, after tight_layout so its reserved bottom
    # margin (for multi-row legends) isn't clobbered.
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


def _plot_boundary_conditions(ax, nodes, bc_type, bc_values, legend_handles, bc_symbol_size=0.03, saved_roller_x=None):
    """
    Plot boundary condition symbols on the mesh.

    BC types, as build_fem_data assigns them:
    0 = free (do nothing)
    1 = fixed, u = v = 0 (small triangle below node)
    2 = x roller, u = 0 and v free (small circle + line, left/right sides)
    3 = y roller — never assigned: the side restraint is 'rollers' (code 2) or
        'fixed' (code 1), and no other step writes a 3
    4 = specified force (vector arrow)

    The legend says what each symbol HOLDS, in words, not the code behind it: a
    reader of the mesh figure is being told which displacements the solution was
    found under, and 'bc_type=3' beside a symbol drawn on the code-2 nodes told
    them the wrong axis in the code's own vocabulary.
    """
    from matplotlib.collections import PatchCollection

    # Get mesh bounds for symbol sizing
    x_min, x_max = nodes[:, 0].min(), nodes[:, 0].max()
    y_min, y_max = nodes[:, 1].min(), nodes[:, 1].max()
    mesh_size = min(x_max - x_min, y_max - y_min)
    symbol_size = mesh_size * bc_symbol_size

    # Fixed boundary conditions (type 1) - triangle below node
    fixed_nodes = np.where(bc_type == 1)[0]
    if len(fixed_nodes) > 0:
        triangle_height = symbol_size
        triangle_width = symbol_size * 0.8
        tri_patches = []
        for node_idx in fixed_nodes:
            x, y = nodes[node_idx]
            tri_patches.append(patches.Polygon([
                [x - triangle_width/2, y - triangle_height],
                [x + triangle_width/2, y - triangle_height],
                [x, y]
            ], closed=True))
        pc = PatchCollection(tri_patches, facecolor='none', edgecolor='red', linewidth=1.5)
        ax.add_collection(pc)

        legend_handles.append(
            plt.Line2D([0], [0], marker='^', color='red', linestyle='None',
                      markersize=8, label='Fixed')
        )

    # X-roller boundary conditions (type 2) - circle + line on left/right sides
    x_roller_nodes = np.where(bc_type == 2)[0]
    if saved_roller_x:
        x_roller_nodes = np.unique(np.concatenate([x_roller_nodes, np.array(sorted(saved_roller_x), dtype=int)]))
    if len(x_roller_nodes) > 0:
        circle_radius = symbol_size * 0.4
        line_length = symbol_size
        x_mid = (x_min + x_max) / 2
        circle_patches = []
        roller_line_segs = []
        for node_idx in x_roller_nodes:
            x, y = nodes[node_idx]
            is_left_side = x < x_mid
            if is_left_side:
                cx = x - circle_radius
                lx = cx - circle_radius
            else:
                cx = x + circle_radius
                lx = cx + circle_radius
            circle_patches.append(patches.Circle((cx, y), circle_radius))
            roller_line_segs.append([[lx, y - line_length/2], [lx, y + line_length/2]])

        pc = PatchCollection(circle_patches, facecolor='none', edgecolor='blue', linewidth=1)
        ax.add_collection(pc)
        lc = LineCollection(roller_line_segs, colors='blue', linewidths=1)
        ax.add_collection(lc)

        legend_handles.append(
            plt.Line2D([0], [0], marker='o', color='blue', linestyle='None',
                      markersize=6, markerfacecolor='none', markeredgewidth=1,
                      label='Roller, free vertically')
        )

    # Specified force boundary conditions (type 4) - vector arrows
    force_nodes = np.where(bc_type == 4)[0]
    if len(force_nodes) > 0:
        force_magnitudes = np.array([np.sqrt(bc_values[ni][0]**2 + bc_values[ni][1]**2) for ni in force_nodes])
        max_force = force_magnitudes.max() if len(force_magnitudes) > 0 else 0
        if max_force > 0:
            scale = symbol_size * 3 / max_force
            gap = symbol_size * 0.5
            arrow_segs = []
            tip_xs, tip_ys = [], []

            for node_idx in force_nodes:
                x, y = nodes[node_idx]
                fx, fy = bc_values[node_idx]
                scaled_fx = fx * scale
                scaled_fy = fy * scale
                mag = np.sqrt(scaled_fx**2 + scaled_fy**2)
                if mag == 0:
                    continue
                ux, uy = scaled_fx / mag, scaled_fy / mag
                tail_x = x - scaled_fx
                tail_y = y - scaled_fy
                tip_x = x - ux * gap
                tip_y = y - uy * gap
                arrow_segs.append([[tail_x, tail_y], [tip_x, tip_y]])
                tip_xs.append(tip_x)
                tip_ys.append(tip_y)

            if arrow_segs:
                lc = LineCollection(arrow_segs, colors='green', linewidths=1.2)
                ax.add_collection(lc)
                ax.plot(tip_xs, tip_ys, marker='v', color='green', markersize=3, linestyle='None')

        legend_handles.append(
            plt.Line2D([0], [0], marker=r'$\rightarrow$', color='green', linestyle='None',
                      markersize=12, label='Applied force')
        )

def plot_fem_results(fem_data, solution, plot_type=['deformation', 'shear_strain', 'displace_vector'],
                    deform_percent=15, show_mesh=True, show_reinforcement=True, figsize=(12, 8), label_elements=False,
                    plot_nodes=False, plot_elements=False, plot_boundary=True, displacement_tolerance=0.5,
                    scale_vectors=True, cmap=None, cbar_shrink=None, save_png=False, save_dxf=False, dpi=300, legend_ncol="auto", legend_frame=False, show_title=True, show_legend=True, fig=None,
                    mesh_on_fields=False, fs=None, failure_solution=None,
                    show_original='outline', deformed_color='k', deform_scale=None,
                    field_state=None, strain_state=None, color_by_magnitude=False, vector_cmap='viridis',
                    vmin=None, vmax=None, vector_max=None):
    """
    Plot FEM results with various visualization options.

    Parameters:
        fem_data: FEM data dictionary from build_fem_data
        solution: FEM solution dictionary from solve_fem
        plot_type: Comma-separated plot types. Valid types:
            'deformation' - deformed mesh overlay
            'displace_mag' - displacement magnitude contours
            'displace_vector' - displacement vectors at corner nodes
            'stress' - von Mises stress contours
            'strain' - equivalent strain contours
            'shear_strain' - viscoplastic max shear strain contours
            'yield' - Mohr-Coulomb yield function contours
        deform_percent: Target deformation as percentage of mesh height (default 15).
        show_mesh: Show mesh lines where the mesh IS the content — the deformation
            panel's original-vs-deformed grid (and the displace_vector panel's edge
            context). It does NOT overlay edges on the filled-field contour panels
            (shear_strain / stress / strain / displace_mag / yield); those edges
            muddy the fill the way they did the seepage flow nets, so they are opt-in
            via ``mesh_on_fields`` instead. This keeps the default three-panel figure
            (deformation / shear_strain / displace_vector) clean: the deformation
            grid stays, the shear-strain band reads uncluttered.
        mesh_on_fields: Overlay light element edges on the filled-field contour
            panels (shear_strain / stress / strain / displace_mag / yield). Off by
            default (the fill is the content); turn on to inspect element boundaries
            against the field.
        show_reinforcement: Show reinforcement elements
        figsize: Figure size (width, height)
        label_elements: Show element ID labels at centroids
        plot_nodes: For displace_vector, show dots at node locations
        plot_elements: For displace_vector, show all element edges
        plot_boundary: For displace_vector, show boundary edges only (default)
        displacement_tolerance: Fraction of max displacement below which vectors are hidden
        scale_vectors: For displace_vector, auto-scale vectors for visibility
        color_by_magnitude: For displace_vector, color each arrow by |u| with a real
            colorbar instead of solid black. Default False keeps the paper-figure
            convention (black arrows, no field colorbar) unchanged.
        vector_cmap: Colormap for displace_vector's color_by_magnitude (default
            'viridis').
        cmap: Color ramp for the shear-strain contours (matplotlib colormap name).
            None keeps the default ('coolwarm').
        cbar_shrink: Colorbar length as a fraction of the axes height (0–1).
            None keeps the automatic size (depends on the number of panels).
        save_png: Save figure to PNG file
        dpi: Resolution for saved PNG
        fs: Optional SSRM factor of safety (the bracket-midpoint result['FS']). When
            given and it differs at display rounding from the last-converged F the
            field was rendered at (solution['F']), the panel titles name both — e.g.
            "FS = 0.45 (rendered at last converged F = 0.43)". When omitted or equal
            at two decimals, titles keep the simple "F = X.XX" form.
        failure_solution: Optional at-failure (unconverged) solve_fem field captured by
            solve_ssrm (result['failure_solution']). When given, the deformation and
            displacement-vector panels render it — the runaway rotational mechanism
            Griffiths & Lane plot — instead of the sub-critical last-converged
            settlement, and the deformation title reads "…at Failure" leading with FS.
            None (default) falls back to ``solution`` for those panels.
        field_state: Which field EVERY result panel renders, when ``failure_solution``
            is given — deformation, displace_vector, AND the filled-contour panels
            (displace_mag / stress / strain / shear_strain / yield) all follow the
            SAME selection, so a multi-panel results figure never mixes states (the
            strain band, the displacement mechanism, and the deformed mesh always
            trace the same field). 'failure' (default) renders every panel from the
            at-failure (unconverged) field, titled "...at Failure  FS = X". 'converged'
            renders every panel from the last-converged ``solution`` instead — the
            deformation panel gets its own auto-scale exaggeration on that field, and
            every panel keeps the established dual-title convention ("FS = X (rendered
            at last converged F = Y)"). With no ``failure_solution``, both selections
            are identical (there's only one field to render) — automatic fallback.
        strain_state: Pre-generalization alias for ``field_state`` (it originally
            governed only the filled-contour panels, before the deformation and
            displace_vector panels were brought under the same switch). Kept for
            backward compatibility; pass either kwarg — ``field_state`` wins if both
            are given.
        show_original: Original-mesh reference on the deformation panel — a tri-state
            passed to plot_deformed_mesh: 'outline' (default) a dashed light boundary
            outline at every density, 'mesh' the full light-gray grid (outline when
            dense), False none.
        deformed_color: Deformed-grid color on the deformation panel (default 'k'
            paper black; 'blue' for the former house style).
        deform_scale: Explicit deformation multiplier. None (default) auto-computes
            it so the rendered field's max displacement is ``deform_percent`` of the
            mesh height; a value overrides the auto-computation.
        vmin, vmax: Contour range for the shear-strain panel, and ``vector_max``
            the |u| the longest displacement arrow stands for. All three are None
            by default, which scales every panel to the field it draws — what a run
            drawn at ONE state wants. A caller drawing the SAME run at two states
            resolves one set of them from both fields (:func:`shared_panel_scales`)
            and passes it to both calls, so a color, an exaggeration and an arrow
            length mean the same thing in both pictures.
    """

    nodes = fem_data["nodes"]
    elements = fem_data["elements"]
    element_types = fem_data["element_types"]
    displacements = solution.get("displacements", np.zeros(2 * len(nodes)))

    # Carry the SSRM factor of safety into the solution the panels see, so titles can
    # name BOTH it and the last-converged F the field was rendered at. Shallow copy:
    # shares the arrays, never mutates the caller's dict (result['last_solution']).
    if fs is not None:
        solution = {**solution, "_ssrm_fs": fs}

    # ``field_state`` is the SINGLE switch governing EVERY result panel — deformation,
    # displace_vector, and the filled-contour panels (displace_mag / stress / strain /
    # shear_strain / yield) alike — so a multi-panel results figure never mixes states.
    # 'failure' (default) puts every panel on the at-failure (unconverged) field
    # solve_ssrm captured; 'converged' puts every panel on the last-converged
    # ``solution`` instead (its own auto-scale exaggeration, and the established dual-
    # title convention). ``strain_state`` is the pre-generalization kwarg name (it used
    # to govern only the contour panels); kept as a back-compat alias.
    if field_state is None:
        field_state = strain_state if strain_state is not None else 'failure'
    if field_state not in ('failure', 'converged'):
        raise ValueError(f"Unknown field_state: '{field_state}'. Valid: 'failure', 'converged'.")

    # The at-failure (unconverged) mechanism solve_ssrm captured, tagged so every
    # panel's title can disclose it ("_at_failure") and carry the SSRM FS. Absent a
    # captured field there's nothing to switch to — ``solution`` (the converged field)
    # serves as both states, so field_state becomes a no-op (automatic fallback).
    if failure_solution is not None:
        failure_field = {**failure_solution, "_at_failure": True}
        if fs is not None:
            failure_field["_ssrm_fs"] = fs
    else:
        failure_field = solution

    # ONE selection drives every panel: deformation and displace_vector read
    # deform_field directly; the filled-contour panels read contour_field, which is
    # simply the SAME field (kept as a separate name for readability at the call
    # sites below) — so the strain band, the displacement mechanism, and the deformed
    # mesh always trace the same F.
    deform_field = failure_field if field_state == 'failure' else solution
    contour_field = deform_field

    # Accept a single string or a list of strings
    if isinstance(plot_type, str):
        plot_types = [plot_type.strip().lower()]
    else:
        plot_types = [pt.strip().lower() for pt in plot_type]
    valid_types = ['displace_mag', 'displace_vector', 'deformation', 'stress', 'strain', 'shear_strain', 'yield']
    
    # Validate plot types
    for pt in plot_types:
        if pt not in valid_types:
            raise ValueError(f"Unknown plot_type: '{pt}'. Valid types: {valid_types}")
    
    # Auto-calculate deformation scale so max displacement is deform_percent of the mesh
    # height, measured on the field the deformation panel actually renders (the at-
    # failure field when present) so the exaggeration is honest. An explicit deform_scale
    # kwarg overrides. Use VP displacement if available (matches plot_deformed_mesh).
    if deform_scale is None:
        deform_scale = deformation_scale(fem_data, deform_field, deform_percent)
    
    # Create subplots based on number of plot types.
    n_plots = len(plot_types)

    # Mesh bounds and a single uniform cushion, in DATA UNITS, applied to BOTH axes.
    # With equal aspect this gives every panel identical visual breathing room on all
    # four sides so the content never welds to the frame; the height derivation and
    # the colorbar (which spans the padded axes box) both inherit it. A few percent
    # of the larger domain dimension — self-adjusting, not a fixed data-unit constant.
    nodes = fem_data["nodes"]
    x_min, x_max = np.min(nodes[:, 0]), np.max(nodes[:, 0])
    y_min, y_max = np.min(nodes[:, 1]), np.max(nodes[:, 1])
    pad = 0.035 * max(x_max - x_min, y_max - y_min)
    if pad <= 0:
        pad = 1.0
    x_margin = y_margin = pad

    # Single-panel plots (the Studio case: one result at a time) fill the space —
    # the real colorbar is placed with make_axes_locatable so it tracks the plot-box
    # height (mirroring plot_seep_solution). The default multi-panel figure stacks
    # equal-aspect strips (deformation / shear_strain / displace_vector) and gives
    # each the SAME make_axes_locatable treatment: the field bar is full-height and
    # every panel keeps an identical-width colorbar slot so the panels stay x-aligned.
    # This "deferred colorbar" layout applies whenever every requested panel is one
    # whose internal colorbar we can suppress; any other field type falls back to the
    # legacy inline colorbars so nothing else regresses.
    single = n_plots == 1
    _deferrable = {'deformation', 'shear_strain', 'displace_vector'}
    defer_cbars = (not single) and all(pt in _deferrable for pt in plot_types)

    own_fig = fig is None
    if not own_fig:
        # Embedded: reuse the caller's figure (GUI canvas). Build the same panel
        # layout on it via fig.subplots instead of creating a new pyplot figure.
        fig.clear()
        if not single and not defer_cbars:
            # Legacy multi-panel uses constrained layout; the deferred-colorbar path
            # uses make_axes_locatable + tight_layout instead (the two don't mix).
            try:
                fig.set_layout_engine("constrained")
            except Exception:
                pass

    if single:
        # With equal aspect, a wide/short slope only fills a thin band of a tall
        # figure, leaving the colorbar towering over the actual plot. Size the
        # figure height to the data aspect ratio so the image fills the figure
        # and the colorbar matches its height.
        data_w = (x_max - x_min) + 2 * x_margin
        data_h = (y_max - y_min) + 2 * y_margin
        if data_w > 0:
            single_height = figsize[0] * (data_h / data_w)
            # Clamp to a sensible range so very flat/steep slopes stay readable
            single_height = float(np.clip(single_height, 2.0, figsize[1]))
        else:
            single_height = figsize[1]
        if own_fig:
            fig, ax = plt.subplots(figsize=(figsize[0], single_height))
        else:
            ax = fig.add_subplot(111)
        axes = [ax]
    else:
        # Content derives the frame: each stacked equal-aspect strip renders at ~the
        # axes width (figure width minus L/R margins and the colorbar slot) times the
        # padded domain aspect, plus a fixed slab for its title/ticks. Summing gives a
        # figure height that packs the panels as compact, identically sized rows with
        # no stretching to fill an over-tall grid cell — the same content-derives-the-
        # frame principle used for the single-panel path. Self-adjusting to the slope's
        # own aspect; the 1.4"/0.6" allowances are structural chrome, not per-figure.
        data_w = (x_max - x_min) + 2 * x_margin
        data_h = (y_max - y_min) + 2 * y_margin
        aspect = (data_h / data_w) if data_w > 0 else 0.4
        axes_w = max(1.0, figsize[0] - 1.4)      # width minus margins + colorbar slot
        panel_h = axes_w * aspect
        total_height = float(np.clip(n_plots * (panel_h + 0.6), 3.0, 60.0))
        if defer_cbars:
            # Plain subplots + tight_layout so the make_axes_locatable colorbars
            # (placed after the loop) compose cleanly; constrained layout fights them.
            if own_fig:
                fig, axes = plt.subplots(n_plots, 1, figsize=(figsize[0], total_height))
            else:
                axes = fig.subplots(n_plots, 1)
        else:
            # Legacy path (non-default panel combos): constrained layout keeps the
            # inline colorbars x-aligned.
            if own_fig:
                fig, axes = plt.subplots(n_plots, 1,
                                         figsize=(figsize[0], total_height),
                                         layout='constrained')
            else:
                axes = fig.subplots(n_plots, 1)
        if not isinstance(axes, (list, np.ndarray)):
            axes = [axes]
        elif isinstance(axes, np.ndarray):
            axes = list(axes)


    # For a single-panel plot, the colorbar is deferred and placed manually
    # (seep-style) after layout so it matches the plot-box height; capture the
    # contour mappable + its label here.
    single_mappable = None
    single_cbar_label = None
    # Reinforcement/pile force colorbar specs deferred out of the shear-strain panel
    # (populated only when reinforcement forces are present) so they can be placed
    # beside the field colorbar without collision.
    reinf_cbar_specs = []
    # displace_vector's OWN deferred mappable (only set when color_by_magnitude is
    # on) — placed on its own panel below, the same mappable/label-return
    # convention as the shear-strain field above.
    vector_mappable = None
    vector_cbar_label = _vector_cbar_label(
        deform_field.get("displacements_elastic", None))

    # In the single-panel case AND the deferred multi-panel case the sub-plotters
    # suppress their own (real or dummy) colorbar; plot_fem_results places the
    # make_axes_locatable colorbars uniformly after the loop.
    defer_panel_cbar = single or defer_cbars

    # Plot each type
    for i, pt in enumerate(plot_types):
        ax = axes[i]

        # Calculate colorbar parameters based on number of plots
        if n_plots == 1:
            cb_shrink = 0.8
            cbar_labelpad = 20
        elif n_plots == 2:
            cb_shrink = 0.7  # Slightly larger than before
            cbar_labelpad = 15
        else:  # 3 or more plots
            cb_shrink = 0.5  # Slightly larger than before
            cbar_labelpad = 12
        # Explicit override from the caller (the Studio colorbar-size control).
        if cbar_shrink is not None:
            cb_shrink = cbar_shrink

        # Filled-field contour panels take ``mesh_on_fields`` (opt-in edge overlay,
        # default off) so the fill reads clean; the deformation panel keeps
        # ``show_mesh`` because there the grid IS the content, and displace_vector
        # keeps it for its edge/boundary context.
        if pt == 'displace_mag':
            plot_displacement_contours(ax, fem_data, contour_field, mesh_on_fields, show_reinforcement,
                                     cbar_shrink=cb_shrink, cbar_labelpad=cbar_labelpad, label_elements=label_elements)
        elif pt == 'displace_vector':
            vector_mappable = plot_displacement_vectors(ax, fem_data, deform_field, show_mesh, show_reinforcement,
                                    cbar_shrink=cb_shrink, cbar_labelpad=cbar_labelpad, label_elements=label_elements,
                                    plot_nodes=plot_nodes, plot_elements=plot_elements, plot_boundary=plot_boundary,
                                    displacement_tolerance=displacement_tolerance, scale_vectors=scale_vectors,
                                    color_by_magnitude=color_by_magnitude, vector_cmap=vector_cmap,
                                    vector_max=vector_max,
                                    single_panel=defer_panel_cbar)
        elif pt == 'deformation':
            plot_deformed_mesh(ax, fem_data, deform_field, deform_scale,
                             show_original=show_original, deformed_color=deformed_color,
                             show_reinforcement=show_reinforcement,
                             cbar_shrink=cb_shrink, cbar_labelpad=cbar_labelpad,
                             label_elements=label_elements, single_panel=defer_panel_cbar,
                             at_failure=deform_field.get("_at_failure", False))
        elif pt == 'stress':
            plot_stress_contours(ax, fem_data, contour_field, mesh_on_fields, show_reinforcement,
                               cbar_shrink=cb_shrink, cbar_labelpad=cbar_labelpad, label_elements=label_elements)
        elif pt == 'strain':
            plot_strain_contours(ax, fem_data, contour_field, mesh_on_fields, show_reinforcement,
                               cbar_shrink=cb_shrink, cbar_labelpad=cbar_labelpad, label_elements=label_elements)
        elif pt == 'shear_strain':
            single_mappable, reinf_cbar_specs = plot_shear_strain_contours(
                ax, fem_data, contour_field, mesh_on_fields, show_reinforcement,
                cbar_shrink=cb_shrink, cbar_labelpad=cbar_labelpad, label_elements=label_elements,
                cmap=cmap, single_panel=defer_panel_cbar, vmin=vmin, vmax=vmax)
            single_cbar_label = SHEAR_STRAIN_LABEL
        elif pt == 'yield':
            plot_yield_function_contours(ax, fem_data, contour_field, mesh_on_fields, show_reinforcement,
                                        cbar_shrink=cb_shrink, cbar_labelpad=cbar_labelpad, label_elements=label_elements)

        # Set consistent axis limits for all plots (including single plots)
        ax.set_xlim(x_min - x_margin, x_max + x_margin)
        ax.set_ylim(y_min - y_margin, y_max + y_margin)
        ax.set_aspect('equal')
        if not show_title:
            ax.set_title("")

    # Single-panel layout (the Studio case: one result shown at a time). When
    # there's a colorbar, attach it with make_axes_locatable so it tracks the
    # (equal-aspect, wide/short) plot's real height instead of towering over it,
    # then a single tight_layout gives symmetric margins AND reserves room for the
    # colorbar's tick + axis labels so nothing is clipped. No hand-tuned margins.
    if single:
        ax = axes[0]
        # Field colorbar (inner) plus any deferred reinforcement/pile force colorbars
        # (outer, separated) — each a full-height slot, no collision. At most one of
        # single_mappable (shear_strain) / vector_mappable (displace_vector) is ever
        # set here since n_plots == 1, so appending the vector spec never collides.
        specs = ([(single_mappable, single_cbar_label)] if single_mappable is not None
                 else []) + list(reinf_cbar_specs)
        if vector_mappable is not None:
            specs = specs + [(vector_mappable, vector_cbar_label)]
        cbars = _place_stacked_cbars(fig, ax, specs) if specs else []
        # A single-panel deformation view carries its Original/Deformed legend inside
        # the axes too (empty corner above the profile), so Studio matches the stack.
        if plot_types[0] == 'deformation':
            _place_deform_legend(axes[0], show_legend)
        try:
            fig.tight_layout()
        except Exception:
            pass
        # Height-adaptive, round-valued ticks (after layout, so each cax has its
        # final drawn height) — the full-height caxes carry them cleanly.
        for cb in cbars:
            if cb is not None:
                adaptive_colorbar_ticks(fig, cb)

    elif defer_cbars:
        # Deferred multi-panel layout: give every stacked panel the SAME set of
        # colorbar slots via make_axes_locatable so the panels stay x-aligned; put the
        # real (full-height) field colorbar — and any deferred reinforcement/pile
        # force bars — on the shear-strain panel, the displace_vector magnitude
        # colorbar (only when color_by_magnitude is on) on ITS OWN panel, and reserve
        # identical invisible slots on whichever panels have none. Slot count is
        # equalized to the WIDER of the two panels' own spec lists (normally just the
        # shear-strain field, unaffected) so every panel stays the same width either
        # way. Each cax tracks its panel's padded, equal-aspect box, so the bars span
        # the frame including the cushion.
        field_idx = plot_types.index('shear_strain') if 'shear_strain' in plot_types else None
        field_specs = ([(single_mappable, single_cbar_label)] if single_mappable is not None
                       else []) + list(reinf_cbar_specs)
        vector_idx = plot_types.index('displace_vector') if 'displace_vector' in plot_types else None
        vector_specs = [(vector_mappable, vector_cbar_label)] if vector_mappable is not None else []
        n_slots = max(len(field_specs), len(vector_specs))
        field_cbars = []
        vector_cbars = []
        for j, ax_j in enumerate(axes):
            if j == field_idx and field_specs:
                padded = field_specs + [(None, "")] * (n_slots - len(field_specs))
                field_cbars = _place_stacked_cbars(fig, ax_j, padded)
            elif j == vector_idx and vector_specs:
                padded = vector_specs + [(None, "")] * (n_slots - len(vector_specs))
                vector_cbars = _place_stacked_cbars(fig, ax_j, padded)
            elif n_slots:
                # Match the widest panel's reserved width with invisible slots.
                _place_stacked_cbars(fig, ax_j, [(None, "")] * n_slots)
        # Original/Deformed legend inside the deformation panel (zero layout height),
        # replacing the old full-width legend band between panels.
        if plot_types[0] == 'deformation':
            _place_deform_legend(axes[0], show_legend)
        try:
            fig.tight_layout()
        except Exception:
            pass
        for cb in field_cbars + vector_cbars:
            if cb is not None:
                adaptive_colorbar_ticks(fig, cb)

    else:
        # Legacy multi-panel (non-default combos): inline colorbars already drawn;
        # just drop the deformation legend inside its panel.
        if plot_types and plot_types[0] == 'deformation':
            _place_deform_legend(axes[0], show_legend)

    if save_png:
        fig.savefig('fem_results.png', dpi=dpi, bbox_inches='tight')
    if save_dxf:
        from .cad import axes_to_dxf
        # One DXF per panel (each plot type), since the figure is multi-panel.
        for i, pt in enumerate(plot_types):
            if i < len(axes):
                axes_to_dxf(axes[i], f'fem_results_{pt}.dxf')

    if own_fig:
        plt.show()
    
    # Return appropriate values
    if n_plots == 1:
        return fig, axes[0]
    else:
        return fig, axes


def plot_displacement_contours(ax, fem_data, solution, show_mesh=True, show_reinforcement=True,
                              cbar_shrink=0.8, cbar_labelpad=20, label_elements=False):
    """
    Plot total displacement magnitude as filled contours using the viridis colormap.
    """
    nodes = fem_data["nodes"]
    elements = fem_data["elements"]
    element_types = fem_data["element_types"]
    displacements = solution.get("displacements", np.zeros(2 * len(nodes)))
    
    # Calculate displacement magnitudes
    u, v = _extract_uv(displacements, fem_data)
    disp_mag = np.sqrt(u**2 + v**2)
    
    # Create triangulation for contouring
    triangles = []
    for i, elem in enumerate(elements):
        elem_type = element_types[i]
        if elem_type == 3:  # Triangle
            triangles.append([elem[0], elem[1], elem[2]])
        elif elem_type == 4:  # Quad - split into triangles
            triangles.append([elem[0], elem[1], elem[2]])
            triangles.append([elem[0], elem[2], elem[3]])
        elif elem_type == 6:  # 6-node triangle - use corner nodes
            triangles.append([elem[0], elem[1], elem[2]])
        elif elem_type in [8, 9]:  # 8-node or 9-node quad - use corner nodes
            triangles.append([elem[0], elem[1], elem[2]])
            triangles.append([elem[0], elem[2], elem[3]])
    
    if triangles:
        triangles = np.array(triangles)
        
        # Create contour plot
        tcf = ax.tricontourf(nodes[:, 0], nodes[:, 1], triangles, disp_mag, gid='DISPLACEMENT_CONTOURS',
                         
                           levels=20, cmap='viridis', alpha=0.8)
        
        # Colorbar
        cbar = ax.figure.colorbar(tcf, ax=ax)
        cbar.set_label(_fem_cbar_label(fem_data, 'Displacement Magnitude', 'length'),
                       rotation=270, labelpad=cbar_labelpad)
        adaptive_colorbar_ticks(ax.figure, cbar)
    
    # Plot mesh
    if show_mesh:
        plot_mesh_lines(ax, fem_data, color='black', alpha=0.3, linewidth=0.5)
    
    # Plot reinforcement
    if show_reinforcement and 'elements_1d' in fem_data:
        plot_reinforcement_lines(ax, fem_data, solution)
    
    # Add element labels if requested
    if label_elements:
        _add_element_labels(ax, fem_data)
    
    ax.set_aspect('equal')
    title = 'Displacement Magnitude Contours'
    if solution.get("_at_failure", False):
        title += ' at Failure'
    ax.set_title(title)


def _get_mesh_boundary(fem_data):
    """
    Compute the boundary edges of the mesh.
    
    Returns:
        boundary_edges: List of (node1, node2) tuples representing boundary edges
    """
    nodes = fem_data["nodes"]
    elements = fem_data["elements"]
    element_types = fem_data["element_types"]
    
    # Count how many times each edge appears
    edge_count = {}
    
    for i, elem in enumerate(elements):
        elem_type = element_types[i]
        
        # Define edges for each element type
        if elem_type == 3:  # Triangle
            edges = [(elem[0], elem[1]), (elem[1], elem[2]), (elem[2], elem[0])]
        elif elem_type == 4:  # Quadrilateral
            edges = [(elem[0], elem[1]), (elem[1], elem[2]), (elem[2], elem[3]), (elem[3], elem[0])]
        elif elem_type == 6:  # 6-node triangle - use corner nodes
            edges = [(elem[0], elem[1]), (elem[1], elem[2]), (elem[2], elem[0])]
        elif elem_type in [8, 9]:  # 8-node or 9-node quad - use corner nodes
            edges = [(elem[0], elem[1]), (elem[1], elem[2]), (elem[2], elem[3]), (elem[3], elem[0])]
        else:
            continue
        
        # Count each edge (both directions)
        for edge in edges:
            # Normalize edge direction (smaller node first)
            normalized_edge = tuple(sorted(edge))
            edge_count[normalized_edge] = edge_count.get(normalized_edge, 0) + 1
    
    # Boundary edges appear only once
    boundary_edges = [edge for edge, count in edge_count.items() if count == 1]
    
    return boundary_edges


# Reference cutoff used SOLELY to freeze the quiver length-scale (see the autoscale
# note in plot_displacement_vectors below). Matches plot_fem_results' own default
# displacement_tolerance, so a call at that default resolves the identical scale a
# plain scale=None autoscale would — default renders are unaffected — while every
# OTHER cutoff reuses that same frozen value instead of computing its own.
_VECTOR_SCALE_REF_TOLERANCE = 0.5

#: The fraction of the section's width the arrow standing for a pinned pair's
#: peak displacement is drawn at (the seepage figures' vector_scale default,
#: plot_seep.py). Pinned panels are drawn in absolute data units off this one
#: number, so a drawn length maps to one displacement in every panel of a pair.
_VECTOR_PIN_FRACTION = 0.05


def plot_displacement_vectors(ax, fem_data, solution, show_mesh=True, show_reinforcement=True,
                             cbar_shrink=0.8, cbar_labelpad=20, label_elements=False,
                             plot_nodes=False, plot_elements=False, plot_boundary=True,
                             displacement_tolerance=1e-6, scale_vectors=True, single_panel=False,
                             color_by_magnitude=False, vector_cmap='viridis',
                             vector_max=None):
    """
    Plot displacement vectors at corner nodes of each element.

    If viscoplastic solution data is available (displacements_elastic in solution),
    plots VP displacement (total - elastic) to show the failure mechanism rather
    than the gravity settlement. Otherwise plots total displacement.

    Parameters:
        ax: Matplotlib axes
        fem_data: FEM data dictionary
        solution: FEM solution dictionary
        show_mesh: Show mesh lines or boundary
        show_reinforcement: Show reinforcement elements
        label_elements: Show element ID labels
        plot_nodes: If True, show dots at all node locations
        plot_elements: If True, show all element edges
        plot_boundary: If True, show only boundary edges (default)
        displacement_tolerance: Fraction of max displacement below which vectors are hidden
        scale_vectors: If True (default), auto-scale vectors for visibility
        color_by_magnitude: If True, color each arrow by its |u| instead of solid
            black, and return the mappable (a real colorbar, labeled "VP Displacement
            Magnitude, |u|") instead of drawing/skipping the blank alignment bar —
            same mappable/label-return convention as plot_shear_strain_contours, so
            plot_fem_results places it with the same _place_stacked_cbars call.
            Default False keeps today's solid-black rendering (and blank alignment
            colorbar) bit-for-bit.
        vector_cmap: Colormap for color_by_magnitude (default 'viridis').
        vector_max: The |u| the LONGEST arrow of the whole SET of panels stands
            for. None (default) lets matplotlib scale this field's own arrows,
            which is what a field drawn alone wants. Given, the panel is drawn in
            absolute data units: vector_max maps to a fixed fraction of the
            section's width (_VECTOR_PIN_FRACTION) and every arrow is |u| times
            that one length-per-displacement, so a drawn length means the same
            displacement in every panel a caller pins to one number (see
            shared_panel_scales). This is the rule the paired seepage velocity
            panels are drawn under (plot_seep.py).

    Returns:
        mappable: The colored Quiver artist when color_by_magnitude is True (for the
            caller to place a colorbar from), else None.
    """
    nodes = fem_data["nodes"]
    elements = fem_data["elements"]
    element_types = fem_data["element_types"]
    displacements = solution.get("displacements", np.zeros(2 * len(nodes)))

    # Use VP displacement (total - elastic) if available, to show failure mechanism
    # This removes the gravity settlement and shows only plastic deformation
    disp_elastic = solution.get("displacements_elastic", None)
    if disp_elastic is not None:
        disp_vp = displacements - disp_elastic
        u, v = _extract_uv(disp_vp, fem_data)
    else:
        u, v = _extract_uv(displacements, fem_data)
    disp_mag = np.sqrt(u**2 + v**2)
    max_disp_mag = np.max(disp_mag)

    if max_disp_mag < 1e-30:
        print("Warning: No VP displacements to plot")
        return None

    # Collect corner nodes only (avoid mid-side nodes of quad8/tri6)
    corner_nodes = set()
    for i, elem in enumerate(elements):
        et = element_types[i]
        if et == 8 or et == 9:
            for j in range(4):
                corner_nodes.add(elem[j])
        elif et == 6:
            for j in range(3):
                corner_nodes.add(elem[j])
        else:
            for j in range(et):
                corner_nodes.add(elem[j])

    # Absolute threshold
    abs_tol = displacement_tolerance * max_disp_mag

    # Element edges (full light-gray mesh) and the boundary outline (black) are
    # independent context layers — either, both, or neither. (show_mesh is kept
    # for API compatibility but no longer gates these; the caller drives them
    # directly via plot_elements / plot_boundary.)
    if plot_elements:
        plot_mesh_lines(ax, fem_data, color='lightgray', alpha=0.5, linewidth=0.5)
    if plot_boundary:
        boundary_edges = _get_mesh_boundary(fem_data)
        for edge in boundary_edges:
            x_coords = [nodes[edge[0], 0], nodes[edge[1], 0]]
            y_coords = [nodes[edge[0], 1], nodes[edge[1], 1]]
            ax.plot(x_coords, y_coords, 'k-', alpha=0.7, linewidth=1.0)

    # Plot small vectors at corner nodes
    corner_list = sorted(corner_nodes)
    cx = nodes[corner_list, 0]
    cy = nodes[corner_list, 1]
    cu = u[corner_list]
    cv = v[corner_list]
    cmag = disp_mag[corner_list]

    mask = cmag > abs_tol

    if np.sum(mask) == 0:
        print("Warning: All displacements below tolerance")
        return None

    # scale_vectors=True: let Matplotlib auto-size the arrows so they are visible
    # (relative magnitudes preserved). scale_vectors=False: draw each arrow at its
    # true displacement magnitude in data units (may be very small for plastic VP
    # displacements), useful for reading actual displacement sizes.
    quiver_style = dict(angles='xy', width=0.002, headwidth=3, headlength=4,
                       headaxislength=3, pivot='tail')
    pinned = (vector_max is not None and float(vector_max) > 0)
    if scale_vectors and pinned:
        # One scale in ABSOLUTE data units for every panel pinned to this
        # vector_max — the rule the paired seepage velocity panels are drawn
        # under (plot_seep.py): the pair's peak displacement is drawn at a fixed
        # fraction of the section's width, and every arrow in every panel is |u|
        # times that same length-per-displacement. Nothing here depends on this
        # panel's own arrow population, so two panels of one pair cannot resolve
        # two scales — which is what multiplying each panel's OWN matplotlib
        # autoscale did (the autoscale is 1.8 * mean(|u|) * max(10, sqrt(N)) over
        # THIS field's visible arrows, a different number in each panel).
        x_range = float(np.max(nodes[:, 0]) - np.min(nodes[:, 0]))
        pin_factor = (x_range * _VECTOR_PIN_FRACTION) / float(vector_max)
        cu = cu * pin_factor
        cv = cv * pin_factor
        scale_kwargs = {"scale_units": "xy", "scale": 1.0}
    elif scale_vectors:
        # mpl's scale=None autoscale is 1.8 * mean(|u|) * max(10, sqrt(N)) over
        # whatever ARRAY is actually handed to quiver() — so deleting the
        # below-cutoff entries BEFORE calling quiver (today's mask[...] indexing)
        # changes both N and the mean, and hence the scale, every time
        # displacement_tolerance changes (Studio's "Vector cutoff" spinner
        # rescales the survivors live as you drag it). Fix: resolve the autoscale
        # once from the FIXED reference cutoff's own field (an invisible temporary
        # quiver, immediately drawn and removed) and reuse that frozen number for
        # the actually-visible arrows, so arrow-length-per-|u| no longer depends on
        # the live cutoff. At the reference cutoff itself — today's default,
        # _VECTOR_SCALE_REF_TOLERANCE — the visible mask IS the reference mask, so
        # this is a no-op short-circuit straight back to today's scale=None call.
        ref_mask = cmag > _VECTOR_SCALE_REF_TOLERANCE * max_disp_mag
        if not np.any(ref_mask):
            ref_mask = mask
        if np.array_equal(ref_mask, mask):
            scale_kwargs = {"scale": None}
        else:
            _ref_q = ax.quiver(cx[ref_mask], cy[ref_mask], cu[ref_mask], cv[ref_mask],
                              alpha=0, scale=None, **quiver_style)
            # Quiver works its autoscale out lazily, at draw time. A figure with
            # no interactive canvas behind it — which is every figure the report
            # renders into — never draws here, and the number came back None; ask
            # the artist for it directly, and fall back to the draw for a backend
            # where that private step is not there.
            try:
                _ref_q._init()
            except Exception:
                pass
            if _ref_q.scale is None:
                ax.figure.canvas.draw()
            scale = _ref_q.scale
            _ref_q.remove()
            if scale is None:
                # Nothing to hold this panel to; mpl scales it as it always did.
                scale_kwargs = {"scale": None}
            else:
                scale_kwargs = {"scale": scale}
    else:
        scale_kwargs = {"scale_units": "xy", "scale": 1.0}

    mappable = None
    if color_by_magnitude:
        from matplotlib.colors import Normalize
        # The colorbar spans the same |u| the arrow lengths do, so a pinned pair
        # is one scale in both readings of the field.
        top = max(float(vector_max), max_disp_mag) if vector_max else max_disp_mag
        _q = ax.quiver(cx[mask], cy[mask], cu[mask], cv[mask], cmag[mask],
                  gid='DISPLACE_VECTORS', cmap=vector_cmap,
                  norm=Normalize(vmin=0.0, vmax=top), alpha=0.9,
                  **quiver_style, **scale_kwargs)
        mappable = _q
    else:
        _q = ax.quiver(cx[mask], cy[mask], cu[mask], cv[mask], gid='DISPLACE_VECTORS',
                  color='black', alpha=0.7, **quiver_style, **scale_kwargs)

    # Plot node dots if requested
    if plot_nodes:
        ax.plot(nodes[:, 0], nodes[:, 1], 'k.', markersize=1, alpha=0.4, gid='MESH_NODES')

    # Plot reinforcement
    if show_reinforcement and 'elements_1d' in fem_data:
        plot_reinforcement_lines(ax, fem_data, solution)

    # Add element labels if requested
    if label_elements:
        _add_element_labels(ax, fem_data)

    # Colorbar. Three cases:
    #   color_by_magnitude, single_panel: defer — return the mappable so the caller
    #     places a REAL full-height colorbar via _place_stacked_cbars, exactly the
    #     mappable/label-return convention plot_shear_strain_contours uses.
    #   color_by_magnitude, not single_panel (legacy inline-colorbar figures): draw
    #     the real colorbar inline now, same styling as the other filled-field panels.
    #   not color_by_magnitude: today's behavior — a blank alignment colorbar so
    #     stacked panels stay x-aligned, skipped entirely for a single-panel plot
    #     (nothing to align to). Bit-for-bit unchanged from before this option existed.
    if color_by_magnitude:
        if not single_panel:
            cbar = ax.figure.colorbar(mappable, ax=ax, shrink=cbar_shrink)
            cbar.set_label(
                _fem_cbar_label(fem_data, _vector_cbar_label(disp_elastic),
                                'length'),
                rotation=270, labelpad=cbar_labelpad)
            adaptive_colorbar_ticks(ax.figure, cbar)
    elif not single_panel:
        dummy_data = np.array([[0, 1]])
        dummy_im = ax.imshow(dummy_data, cmap='viridis', alpha=0)
        cbar = ax.figure.colorbar(dummy_im, ax=ax, shrink=cbar_shrink)
        cbar.set_label('', color='white')
        cbar.set_ticks([])
        cbar.set_ticklabels([])
        cbar.outline.set_color('white')
        cbar.outline.set_linewidth(0)

    F = solution.get("F", None)
    # One name for the panel — the documentation's, and the results view's — with
    # the field it drew as a qualifier, not as a second name for the same thing.
    title = ('Displacement Vectors (viscoplastic)' if disp_elastic is not None
             else 'Displacement Vectors')
    # at_failure when the panel is drawing the captured unconverged field (routed here
    # by plot_fem_results): leads with FS; "at Failure" already discloses the state.
    title = _fs_title(title, F, solution.get("_ssrm_fs"),
                      at_failure=solution.get("_at_failure", False))
    ax.set_title(title, fontsize=12, pad=15)

    return mappable if single_panel else None


def plot_stress_contours(ax, fem_data, solution, show_mesh=True, show_reinforcement=True,
                        cbar_shrink=0.8, cbar_labelpad=20, label_elements=False):
    """
    Plot von Mises stress contours with yielding elements highlighted.
    """
    nodes = fem_data["nodes"]
    elements = fem_data["elements"]
    element_types = fem_data["element_types"]
    stresses = solution.get("stresses", np.zeros((len(elements), 4)))
    
    # Use yield function to determine plastic elements for consistency
    # If yield_function is available, use it; otherwise fall back to plastic_elements
    yield_function = solution.get("yield_function", None)
    if yield_function is not None:
        plastic_elements = yield_function > 0  # F > 0 means yielding
    else:
        plastic_elements = solution.get("plastic_elements", np.zeros(len(elements), dtype=bool))
    
    # Extract von Mises stresses
    von_mises = stresses[:, 3]  # 4th column is von Mises stress
    
    # Create element patches with color based on stress
    patches_list = []
    stress_values = []
    
    for i, elem in enumerate(elements):
        elem_type = element_types[i]
        if elem_type == 3:  # Triangle
            coords = nodes[elem[:3]]
            patch = Polygon(coords, closed=True)
            patches_list.append(patch)
            stress_values.append(von_mises[i])
        elif elem_type == 4:  # Quadrilateral
            coords = nodes[elem[:4]]
            patch = Polygon(coords, closed=True)
            patches_list.append(patch)
            stress_values.append(von_mises[i])
        elif elem_type == 6:  # 6-node triangle - use corner nodes
            coords = nodes[elem[:3]]
            patch = Polygon(coords, closed=True)
            patches_list.append(patch)
            stress_values.append(von_mises[i])
        elif elem_type in [8, 9]:  # 8-node or 9-node quad - use corner nodes
            coords = nodes[elem[:4]]
            patch = Polygon(coords, closed=True)
            patches_list.append(patch)
            stress_values.append(von_mises[i])
    
    if patches_list:
        from matplotlib.collections import PatchCollection
        
        # Create patch collection
        p = PatchCollection(patches_list, alpha=0.8, edgecolors='none', gid='STRESS_CONTOURS')
        p.set_array(np.array(stress_values))
        p.set_cmap('plasma')
        ax.add_collection(p)
        
        # Colorbar
        cbar = ax.figure.colorbar(p, ax=ax)
        cbar.set_label(_fem_cbar_label(fem_data, 'von Mises Stress', 'stress'),
                       rotation=270, labelpad=cbar_labelpad)
        adaptive_colorbar_ticks(ax.figure, cbar)
    
    # Highlight plastic elements with thick boundary
    if np.any(plastic_elements):
        for i, elem in enumerate(elements):
            if plastic_elements[i]:
                elem_type = element_types[i]
                if elem_type == 3:  # Triangle
                    coords = nodes[elem[:3]]
                    coords = np.vstack([coords, coords[0]])  # Close the polygon
                    ax.plot(coords[:, 0], coords[:, 1], 'r-', linewidth=2, alpha=0.8)
                elif elem_type == 4:  # Quadrilateral
                    coords = nodes[elem[:4]]
                    coords = np.vstack([coords, coords[0]])  # Close the polygon
                    ax.plot(coords[:, 0], coords[:, 1], 'r-', linewidth=2, alpha=0.8)
                elif elem_type == 6:  # 6-node triangle - use corner nodes
                    coords = nodes[elem[:3]]
                    coords = np.vstack([coords, coords[0]])  # Close the polygon
                    ax.plot(coords[:, 0], coords[:, 1], 'r-', linewidth=2, alpha=0.8)
                elif elem_type in [8, 9]:  # 8-node or 9-node quad - use corner nodes
                    coords = nodes[elem[:4]]
                    coords = np.vstack([coords, coords[0]])  # Close the polygon
                    ax.plot(coords[:, 0], coords[:, 1], 'r-', linewidth=2, alpha=0.8)
    
    # Plot mesh
    if show_mesh:
        plot_mesh_lines(ax, fem_data, color='gray', alpha=0.3, linewidth=0.3)
    
    # Plot reinforcement with force visualization
    if show_reinforcement and 'elements_1d' in fem_data:
        plot_reinforcement_forces(ax, fem_data, solution)
    
    # Add element labels if requested
    if label_elements:
        _add_element_labels(ax, fem_data)
    
    ax.set_aspect('equal')
    title = 'von Mises Stress (Red outline = Yielding/Plastic Elements)'
    if solution.get("_at_failure", False):
        title += ' at Failure'
    ax.set_title(title)


def plot_deformed_mesh(ax, fem_data, solution, deform_scale=1.0,
                       show_original='outline', deformed_color='k', show_reinforcement=True,
                       cbar_shrink=0.8, cbar_labelpad=20, label_elements=False,
                       single_panel=False, at_failure=False, show_mesh=None):
    """
    Plot deformed mesh overlay on original mesh.

    Draws the deformed mesh in ``deformed_color`` (default paper black 'k') over an
    original-mesh reference governed by ``show_original``. Uses VP displacement
    (total - elastic) when available to show the failure mechanism rather than gravity
    settlement. The deform_scale parameter amplifies the displacements for visibility.

    Parameters:
        deform_scale: Displacement multiplier for visibility.
        show_original: Original (undeformed) reference — a tri-state:
            'outline' (default) draws a light DASHED boundary outline only, at every
                mesh density (Griffiths & Lane convention; Norm's Variant E look);
            'mesh' draws the full light-gray original grid, collapsing to the dashed
                outline when the mesh is dense enough that two grids would tangle;
            False/None draws no original reference.
        deformed_color: Color of the deformed grid (default 'k'; pass 'blue' for the
            former house style).
        at_failure: When True the title reads "…at Failure" and leads with the SSRM
            factor of safety (the field is the unconverged at-failure state).
        show_mesh: Legacy boolean alias for ``show_original`` (True→'mesh',
            False→off); prefer ``show_original``. None (default) leaves show_original
            in effect.
    """
    # Legacy alias: an explicit show_mesh bool maps onto the tri-state so older
    # callers keep working (True == the former full-grid-or-outline behavior).
    if show_mesh is not None:
        show_original = 'mesh' if show_mesh else False

    nodes = fem_data["nodes"]
    elements = fem_data["elements"]
    element_types = fem_data["element_types"]
    displacements = solution.get("displacements", np.zeros(2 * len(nodes)))

    # Use VP displacement (total - elastic) if available, to show failure mechanism
    disp_elastic = solution.get("displacements_elastic", None)
    if disp_elastic is not None:
        disp = displacements - disp_elastic
    else:
        disp = displacements

    # Calculate deformed node positions
    u, v = _extract_uv(disp, fem_data)
    nodes_deformed = nodes + deform_scale * np.column_stack([u, v])
    
    # Mesh-line hierarchy with a device-pixel floor so nothing renders sub-pixel
    # (which smears each line into an antialiased haze). The deformed grid is the
    # story: crisp, full opacity, at the floored adaptive width, in deformed_color.
    # The original is a light reference governed by show_original (below).
    floor_pt = _mesh_pixel_floor_pt(ax)
    lw = _adaptive_mesh_linewidth(ax, fem_data, floor_pt)
    _, n_across = _mesh_edge_stats(fem_data)
    tangle = _rendered_edge_px(ax, n_across) < _DENSE_TANGLE_EDGE_PX

    fem_data_deformed = fem_data.copy()
    fem_data_deformed["nodes"] = nodes_deformed

    # Original (undeformed) reference — tri-state show_original:
    #   'mesh'    : full light-gray grid, collapsing to the dashed boundary outline
    #               when the mesh is dense enough that two interleaved grids tangle
    #               (median rendered edge below a few device px).
    #   'outline' : the dashed light boundary outline ONLY, at every density — the
    #               Griffiths & Lane convention (the deformed grid carries the
    #               deformation, the outline marks the undeformed extent). Default.
    #   False/None: no original reference.
    if show_original == 'mesh':
        if tangle:
            _plot_boundary_outline(ax, fem_data, color='0.55', alpha=0.6, linewidth=lw,
                                   linestyle='--', label='Original (outline)')
        else:
            plot_mesh_lines(ax, fem_data, color='lightgray', alpha=0.5,
                            linewidth=lw, label='Original')
    elif show_original:      # 'outline' (default; any truthy value that isn't 'mesh')
        _plot_boundary_outline(ax, fem_data, color='0.55', alpha=0.7,
                               linewidth=max(lw, floor_pt), linestyle='--',
                               label='Original (outline)')

    # Plot deformed mesh
    plot_mesh_lines(ax, fem_data_deformed, color=deformed_color, alpha=1.0,
                    linewidth=lw, label='Deformed')
    
    # Plot members in both original and deformed configurations. The two
    # configurations are told apart by COLOR, on both kinds of member: the
    # original is the muted reference (gray reinforcement, the pile's own green)
    # and the deformed is drawn in red, the same red the deformed grid's story
    # is told in. Both pile overlays were green and lay one on top of the other.
    if show_reinforcement and 'elements_1d' in fem_data:
        plot_reinforcement_lines(ax, fem_data, solution, color='gray', alpha=0.5,
                                 linewidth=2, label='Original Reinforcement',
                                 pile_color=PILE_COLOR)
        plot_reinforcement_lines(ax, fem_data_deformed, solution, color='red',
                                 alpha=0.8, linewidth=2,
                                 label='Deformed Reinforcement',
                                 pile_color=PILE_DEFORMED_COLOR)
    
    # Add element labels if requested
    if label_elements:
        _add_element_labels(ax, fem_data_deformed)  # Label on deformed mesh
    
    # Add a dummy colorbar to maintain consistent spacing with other plots so the
    # x-axis alignment stays consistent across stacked subplots. Skipped for a
    # single-panel plot (the Studio case), where there's nothing to align to and
    # the invisible colorbar would just steal the right margin.
    if not single_panel:
        dummy_data = np.array([[0, 1]])
        dummy_im = ax.imshow(dummy_data, cmap='viridis', alpha=0)
        cbar = ax.figure.colorbar(dummy_im, ax=ax, shrink=cbar_shrink)
        cbar.set_label('Deformation Scale', rotation=270, labelpad=cbar_labelpad, color='white')
        cbar.set_ticks([])  # Remove tick marks
        cbar.set_ticklabels([])  # Remove tick labels

        # Make the colorbar completely invisible by setting colors to background
        cbar.outline.set_color('white')  # Make the border invisible
        cbar.outline.set_linewidth(0)    # Remove the border line

    # Note: Axis limits will be set by the calling function for consistent multi-plot alignment
    # When used as a standalone plot, matplotlib will auto-scale appropriately
    F = solution.get("F", None)
    disp_label = 'Viscoplastic Deformation' if disp_elastic is not None else 'Mesh Deformation'
    if at_failure:
        disp_label += ' at Failure'
    scale_str = f'{deform_scale:.0f}' if deform_scale >= 10 else f'{deform_scale:.1f}'
    base = f'{disp_label} (Scale = {scale_str}x)'
    # The at-failure field is the UNCONVERGED state, solved a margin beyond critical to
    # develop the mechanism; _fs_title(at_failure=...) leads with FS — "at Failure"
    # already carries the disclosure, so no trial-F clause.
    title = _fs_title(base, F, solution.get("_ssrm_fs"), at_failure=at_failure)
    ax.set_title(title, fontsize=12, pad=15)


def _add_element_labels(ax, fem_data):
    """
    Add element ID labels at element centers.
    """
    nodes = fem_data["nodes"]
    elements = fem_data["elements"]
    element_types = fem_data["element_types"]
    
    for i, elem in enumerate(elements):
        elem_type = element_types[i]
        
        # Get element nodes for centroid calculation
        if elem_type == 3:  # Triangle
            elem_nodes = nodes[elem[:3]]
        elif elem_type == 4:  # Quad
            elem_nodes = nodes[elem[:4]]
        elif elem_type == 6:  # 6-node triangle - use corner nodes
            elem_nodes = nodes[elem[:3]]
        elif elem_type in [8, 9]:  # 8 or 9-node quad - use corner nodes
            elem_nodes = nodes[elem[:4]]
        else:
            continue
            
        # Calculate centroid
        centroid = np.mean(elem_nodes, axis=0)
        
        # Add label (1-based indexing for display)
        ax.text(centroid[0], centroid[1], str(i+1),
                ha='center', va='center', fontsize=6, 
                color='darkblue', alpha=0.7, zorder=100)


def _mesh_edge_stats(fem_data):
    """(median element edge length in data units, elements-across-the-domain-width).

    Measures actual corner-edge lengths rather than inferring from element count, so
    it stays honest on refined / mixed-density meshes. Sampled for speed on big
    meshes."""
    nodes = fem_data["nodes"]
    elements = fem_data["elements"]
    element_types = fem_data.get("element_types", None)
    if element_types is None or len(elements) == 0:
        return 0.0, 20.0
    step = max(1, len(elements) // 2000)
    lengths = []
    for i in range(0, len(elements), step):
        et = element_types[i]
        nc = 3 if et in (3, 6) else 4 if et in (4, 8, 9) else int(et)
        pts = nodes[elements[i][:nc]]
        d = np.hypot(*(pts - np.roll(pts, -1, axis=0)).T)
        lengths.append(np.median(d))
    med_edge = float(np.median(lengths)) if lengths else 0.0
    x_w = float(nodes[:, 0].max() - nodes[:, 0].min())
    n_across = (x_w / med_edge) if med_edge > 0 else 20.0
    return med_edge, max(n_across, 1.0)


def _rendered_edge_px(ax, n_across):
    """Median element edge length in device pixels at the axes' current drawn width
    — the on-screen element size that decides whether two interleaved grids tangle."""
    fig = ax.figure
    try:
        axes_w_px = float(ax.get_window_extent().width)
    except Exception:
        axes_w_px = 0.0
    if axes_w_px <= 1.0:
        axes_w_px = fig.get_size_inches()[0] * fig.dpi * 0.8
    return axes_w_px / max(n_across, 1.0)


def _mesh_pixel_floor_pt(ax):
    """Linewidth floor in points that renders as one full device pixel at the
    figure's dpi (72 pt/inch ÷ dpi). Lines never go below this, so they stay crisp
    instead of smearing into a sub-pixel antialiased haze; density beyond the floor
    is carried by alpha / the boundary-outline fallback, never by thinner geometry."""
    dpi = ax.figure.dpi or 100.0
    return 72.0 / dpi


def _adaptive_mesh_linewidth(ax, fem_data, floor_pt=None):
    """Mesh-line weight from the RENDERED element size, so a dense mesh reads as a
    legible hairline grid instead of fusing into a solid blob while a coarse teaching
    mesh keeps ~its usual weight. The on-screen element size (points) = axes width
    (points) / elements-across; the line is a fixed fraction of that. Calibrated so a
    ~20-element-across mesh renders at ~1.0 pt, then clamped to [pixel floor, 1.5 pt]
    — the floor keeps lines a crisp full device pixel (no sub-pixel blur). Uses the
    axes' current drawn width (tracks the Studio canvas) and the mesh's own edge
    stats (robust to refinement)."""
    med_edge, n_across = _mesh_edge_stats(fem_data)
    if floor_pt is None:
        floor_pt = _mesh_pixel_floor_pt(ax)
    if med_edge <= 0:
        return max(1.0, floor_pt)
    fig = ax.figure
    try:
        axes_w_px = float(ax.get_window_extent().width)
    except Exception:
        axes_w_px = 0.0
    if axes_w_px <= 1.0:
        axes_w_px = fig.get_size_inches()[0] * fig.dpi * 0.8
    axes_w_pts = axes_w_px * 72.0 / fig.dpi
    elem_pts = axes_w_pts / n_across
    lw = 0.030 * elem_pts
    return float(np.clip(lw, floor_pt, 1.5))


def plot_mesh_lines(ax, fem_data, color='black', alpha=1.0, linewidth=None, label=None):
    """
    Plot mesh element boundaries.

    linewidth=None (the default) sizes the line adaptively from the rendered element
    size (:func:`_adaptive_mesh_linewidth`) so dense meshes stay a legible hairline
    grid; pass an explicit value to force a fixed weight.
    """
    nodes = fem_data["nodes"]
    elements = fem_data["elements"]
    element_types = fem_data["element_types"]

    if linewidth is None:
        linewidth = _adaptive_mesh_linewidth(ax, fem_data)

    lines = []
    for i, elem in enumerate(elements):
        elem_type = element_types[i]
        if elem_type == 3:  # Triangle
            # Add triangle edges
            edges = [(elem[0], elem[1]), (elem[1], elem[2]), (elem[2], elem[0])]
        elif elem_type == 4:  # Quadrilateral
            # Add quad edges
            edges = [(elem[0], elem[1]), (elem[1], elem[2]), (elem[2], elem[3]), (elem[3], elem[0])]
        elif elem_type == 6:  # 6-node triangle - use corner nodes
            edges = [(elem[0], elem[1]), (elem[1], elem[2]), (elem[2], elem[0])]
        elif elem_type in [8, 9]:  # 8-node or 9-node quad - use corner nodes
            edges = [(elem[0], elem[1]), (elem[1], elem[2]), (elem[2], elem[3]), (elem[3], elem[0])]
        else:
            continue
        
        for edge in edges:
            line_coords = nodes[[edge[0], edge[1]]]
            lines.append(line_coords)
    
    if lines:
        lc = LineCollection(lines, colors=color, alpha=alpha, linewidths=linewidth, label=label, gid='MESH')
        ax.add_collection(lc)


def _plot_boundary_outline(ax, fem_data, color='0.55', alpha=0.6, linewidth=1.0,
                           linestyle='-', label=None):
    """Draw only the mesh's outer boundary (domain outline), no interior element
    edges — the stand-in for the original grid so it never tangles with the deformed
    grid drawn on top. ``linestyle`` defaults to solid; the deformation panel passes
    '--' for the dashed undeformed-extent outline (Griffiths & Lane convention)."""
    nodes = fem_data["nodes"]
    segs = [[nodes[a], nodes[b]] for a, b in _get_mesh_boundary(fem_data)]
    if segs:
        lc = LineCollection(segs, colors=color, alpha=alpha, linewidths=linewidth,
                            linestyles=linestyle, label=label, gid='MESH')
        ax.add_collection(lc)


#: What a pile is drawn in wherever its geometry is shown: the package-wide pile
#: color (xslope.style's "piles"), so a pile looks the same on every figure.
PILE_COLOR = 'green'

#: What a pile is drawn in where a figure shows it DEFORMED beside its original
#: position. The deformation panel already draws the deformed reinforcement in
#: red over the original in gray, so red on that panel means "this is the
#: deformed configuration" and not "this is a reinforcement line" — the
#: reinforcement's own geometry color is gray. Original and deformed piles were
#: both drawn in green and could not be told apart at any exaggeration.
PILE_DEFORMED_COLOR = 'red'


def plot_reinforcement_lines(ax, fem_data, solution, color='red', alpha=1.0,
                             linewidth=2, label=None, pile_color=PILE_COLOR):
    """
    Plot reinforcement and pile elements as lines with distinct colors.

    ``color`` is the reinforcement color; ``pile_color`` the pile color, which
    defaults to the package-wide green. A caller drawing the same members twice
    — original and deformed — passes a different color for each configuration on
    BOTH kinds, or the two overlays are indistinguishable.
    """
    if 'elements_1d' not in fem_data:
        return

    nodes = fem_data["nodes"]
    elements_1d = fem_data["elements_1d"]
    element_types_1d = fem_data["element_types_1d"]
    pile_elem_mask = fem_data.get("pile_elem_mask", np.zeros(len(elements_1d), dtype=bool))

    reinf_lines = []
    pile_lines = []
    for i, elem in enumerate(elements_1d):
        elem_type = element_types_1d[i]
        if elem_type >= 2:
            line_coords = nodes[elem[:2]]
            if pile_elem_mask[i]:
                pile_lines.append(line_coords)
            else:
                reinf_lines.append(line_coords)

    if reinf_lines:
        lc = LineCollection(reinf_lines, colors=color, alpha=alpha, linewidths=linewidth, label=label, gid='REINFORCEMENT')
        ax.add_collection(lc)
    if pile_lines:
        pile_label = label.replace('Reinforcement', 'Pile') if label and 'Reinforcement' in label else None
        lc = LineCollection(pile_lines, colors=pile_color, alpha=alpha, linewidths=linewidth + 1, label=pile_label, gid='PILES')
        ax.add_collection(lc)


#: The legend entries that state something about a member's FORCE. Every one of
#: them is a claim read off a solved array, so none may appear over a solution
#: that carries no such array — see :func:`_member_forces_in_solution`. Named
#: once here so ``test_a_member_overlay_claims_no_force_it_was_not_given`` can
#: pin their absence without restating the strings.
MEMBER_FORCE_LEGEND_LABELS = (
    'Inactive (no tension)',
    'At residual (Tres)',
    'Pulled out',
)

#: What a member drawn as geometry alone is called in the legend: the kind of
#: member it is, and nothing about what it carries.
REINFORCEMENT_GEOMETRY_LABEL = 'Reinforcement line'
PILE_GEOMETRY_LABEL = 'Pile line'

#: The neutral pair a geometry-only member is drawn in — an outline and a core,
#: the same weights the classified draws use, in one color that ranks nothing.
_GEOMETRY_MEMBER_COLOR = '#BBBBBB'


def _elem_flags(solution, key, n):
    """A per-1D-element boolean array from ``solution``, length ``n``.

    A solution that never recorded the flag, or recorded it at some other
    length, yields all-False: the state is unknown, and an unknown state is one
    the figure does not mark.
    """
    values = (solution or {}).get(key)
    if values is None:
        return np.zeros(n, dtype=bool)
    values = np.asarray(values, dtype=bool)
    return values if values.shape == (n,) else np.zeros(n, dtype=bool)


def _member_forces_in_solution(solution, key):
    """Whether ``solution`` actually carries the solved force array ``key``.

    A solution reloaded from companion files written before the member results
    existed has the nodal and element fields and nothing else. Substituting
    zeros for the array that is missing does not produce a neutral drawing: it
    produces a confident one. Every bar classifies as carrying no tension and
    the overlay prints "Inactive (no tension)" over six members whose force was
    never computed. This is the predicate
    :func:`xslope.report._member_forces_recorded` applies to the same arrays
    before it prints a forces table.
    """
    values = solution.get(key)
    return values is not None and len(values) > 0


def plot_reinforcement_forces(ax, fem_data, solution, draw_cbar=True):
    """
    Plot reinforcement elements colored by force level.

    Color scheme:
    - Blue to green to yellow to red: 0 to Tmax (tension force ramp)
    - Magenta: element has softened to its residual capacity Tres
    - Black: element has pulled out — softened with no residual left, carrying
      nothing
    - Green: element carrying no tension (inactive or in compression)

    The two marked states are read from ``softened_1d_elements``, the array the
    solver sets when an element DROPS to its residual capacity. An element that
    has reached its peak capacity and is holding there is a different thing: it
    is at Tmax, not at Tres, and it rides the force ramp at the color its force
    earns. Classifying from ``failed_1d_elements`` instead — the yield latch,
    which the solver sets for reporting the moment an element reaches its cap —
    put every bar that was holding its full capacity under a legend entry
    reading "At residual", which is the one state those bars were not in.

    Where the solution carries no force array for a kind of member, that kind is
    drawn as GEOMETRY ONLY: one neutral color, a legend entry naming the kind,
    no classification and no colorbar. The force colors above are read off a
    solved array, and a solution that never recorded one cannot be colored by it.
    The two kinds are judged separately — a run that solved reinforcement forces
    and no pile forces colors the bars and draws the piles plain.

    draw_cbar=True (default) draws the force colorbar(s) inline on ``ax`` with
    ``colorbar(ax=ax)`` — fine when it is the only colorbar. When a field colorbar
    (e.g. VP shear strain) is also present, pass draw_cbar=False: the force
    colorbar(s) are NOT drawn and their (ScalarMappable, label) specs are returned so
    the caller can lay them out beside the field bar without collision (both via
    make_axes_locatable). Always returns a list of (mappable, label) specs.
    """
    cbar_specs = []
    if 'elements_1d' not in fem_data:
        return cbar_specs

    from matplotlib.colors import LinearSegmentedColormap
    import matplotlib.cm as cm

    nodes = fem_data["nodes"]
    elements_1d = fem_data["elements_1d"]
    reinf_forces_known = _member_forces_in_solution(solution, "forces_1d")
    forces_1d = solution.get("forces_1d", np.zeros(len(elements_1d)))
    t_allow = fem_data.get("t_allow_by_1d_elem", np.ones(len(elements_1d)))
    t_res = fem_data.get("t_res_by_1d_elem", np.zeros(len(elements_1d)))
    softened_1d = _elem_flags(solution, "softened_1d_elements", len(elements_1d))

    # Find global Tmax (max of all t_allow values)
    t_max_global = t_allow.max() if len(t_allow) > 0 else 1.0

    # Custom colormap: blue -> white -> red (coolwarm style)
    force_cmap = LinearSegmentedColormap.from_list(
        'force_ramp', ['#2166ac', '#f7f7f7', '#d73027'], N=256)

    # Classify and draw each element
    normal_lines = []
    normal_colors = []
    tres_lines = []
    pullout_lines = []
    inactive_lines = []

    pile_elem_mask = fem_data.get("pile_elem_mask", np.zeros(len(elements_1d), dtype=bool))
    pile_force_lines = []
    pile_force_colors = []
    forces_pile_lateral = solution.get("forces_pile_lateral", np.array([]))
    pile_forces_known = _member_forces_in_solution(solution, "forces_pile_lateral")

    # Members whose force the solution never recorded: drawn where they are, in
    # one neutral color, named by kind and classified as nothing.
    reinf_geometry_lines = []
    pile_geometry_lines = []

    # Build pile element index mapping: global 1d index -> pile force index
    pile_force_idx = 0

    for i in range(len(elements_1d)):
        elem = elements_1d[i]
        coords = nodes[elem[:2]]

        if pile_elem_mask[i]:
            # Pile element — color by lateral (shear) force, where one was solved.
            if not pile_forces_known:
                pile_geometry_lines.append(coords)
            elif pile_force_idx < len(forces_pile_lateral):
                pile_force_lines.append(coords)
                pile_force_colors.append(abs(forces_pile_lateral[pile_force_idx]))
            pile_force_idx += 1
            continue

        if not reinf_forces_known:
            reinf_geometry_lines.append(coords)
            continue

        force = forces_1d[i]
        # Softened: the element dropped off its peak onto whatever residual it
        # was given. With no residual left and no force in it, that drop was a
        # pullout — the same three-part definition
        # :func:`xslope.fem_details.reinforcement_profile` marks a pullout by,
        # so the overlay and the member detail figure mark the same elements.
        is_softened = softened_1d[i]

        if is_softened and t_res[i] < 1e-6 and force < 1e-6:
            pullout_lines.append(coords)
        elif is_softened:
            tres_lines.append(coords)
        elif force > 1e-6:
            ratio = min(force / t_max_global, 1.0) if t_max_global > 0 else 0.0
            normal_lines.append(coords)
            normal_colors.append(force_cmap(ratio))
        else:
            inactive_lines.append(coords)

    # Draw the members whose force was never solved: geometry, neutral, named by
    # kind. No ramp, no state, no colorbar — there is nothing to color them by.
    for geometry_lines, geometry_label in ((reinf_geometry_lines, REINFORCEMENT_GEOMETRY_LABEL),
                                           (pile_geometry_lines, PILE_GEOMETRY_LABEL)):
        if not geometry_lines:
            continue
        lc_outline = LineCollection(geometry_lines, colors='black', linewidths=4.5,
                                    alpha=0.9, zorder=3.9)
        ax.add_collection(lc_outline)
        lc = LineCollection(geometry_lines, colors=_GEOMETRY_MEMBER_COLOR, linewidths=3,
                            alpha=0.9, zorder=4)
        ax.add_collection(lc)
        ax.plot([], [], '-', color=_GEOMETRY_MEMBER_COLOR, linewidth=3, alpha=0.9,
                label=geometry_label)

    # Draw inactive elements (green, solid)
    if inactive_lines:
        lc_outline = LineCollection(inactive_lines, colors='black', linewidths=4.5, alpha=0.9, zorder=3.9)
        ax.add_collection(lc_outline)
        lc = LineCollection(inactive_lines, colors='#00CC00', linewidths=3, alpha=0.9, zorder=4)
        ax.add_collection(lc)
        ax.plot([], [], '-', color='#00CC00', linewidth=3, alpha=0.9, label='Inactive (no tension)')

    # Draw normal tension elements (force-colored)
    if normal_lines:
        lc_outline = LineCollection(normal_lines, colors='black', linewidths=4.5, alpha=0.9, zorder=4.9)
        ax.add_collection(lc_outline)
        lc = LineCollection(normal_lines, colors=normal_colors, linewidths=3, alpha=0.9, zorder=5)
        ax.add_collection(lc)

        # Force colorbar — draw inline, or hand back the spec for collision-free
        # placement beside the field colorbar.
        sm = cm.ScalarMappable(cmap=force_cmap, norm=plt.Normalize(0, t_max_global))
        sm.set_array([])
        # Reinforcement/pile forces are per unit width of slope (force/length). Label
        # once so both the inline colorbar and the deferred stacked-colorbar spec agree.
        reinf_force_label = _fem_cbar_label(fem_data, 'Reinforcement Force', 'force_per_len')
        cbar_specs.append((sm, reinf_force_label))
        if draw_cbar:
            cbar = ax.figure.colorbar(sm, ax=ax, shrink=0.6, pad=0.02)
            cbar.set_label(reinf_force_label, rotation=270, labelpad=15, fontsize=10)

    # Draw elements that softened onto their residual capacity (magenta)
    if tres_lines:
        lc_outline = LineCollection(tres_lines, colors='black', linewidths=4.5, alpha=0.9, zorder=5.9)
        ax.add_collection(lc_outline)
        lc = LineCollection(tres_lines, colors='magenta', linewidths=3, alpha=0.9, zorder=6)
        ax.add_collection(lc)
        ax.plot([], [], '-', color='magenta', linewidth=3, label='At residual (Tres)')

    # Draw pulled-out elements (black, solid)
    if pullout_lines:
        lc_outline = LineCollection(pullout_lines, colors='black', linewidths=4.5, alpha=0.9, zorder=5.9)
        ax.add_collection(lc_outline)
        lc = LineCollection(pullout_lines, colors='black', linewidths=3, alpha=0.9, zorder=6)
        ax.add_collection(lc)
        ax.plot([], [], '-', color='black', linewidth=3, alpha=0.9, label='Pulled out')

    # Draw pile elements colored by lateral (shear) force
    if pile_force_lines:
        from matplotlib.colors import Normalize
        max_lateral = max(pile_force_colors) if pile_force_colors else 1.0
        pile_cmap = plt.cm.Greens
        pile_norm = Normalize(vmin=0, vmax=max_lateral if max_lateral > 0 else 1.0)
        colors = [pile_cmap(pile_norm(v)) for v in pile_force_colors]
        lc_outline = LineCollection(pile_force_lines, colors='black', linewidths=5, alpha=0.9, zorder=5.9)
        ax.add_collection(lc_outline)
        lc = LineCollection(pile_force_lines, colors=colors, linewidths=3.5, alpha=0.9, zorder=6)
        ax.add_collection(lc)
        sm = cm.ScalarMappable(cmap=pile_cmap, norm=pile_norm)
        sm.set_array([])
        pile_force_label = _fem_cbar_label(fem_data, 'Pile Shear Force', 'force_per_len')
        cbar_specs.append((sm, pile_force_label))
        if draw_cbar:
            cbar = ax.figure.colorbar(sm, ax=ax, shrink=0.6, pad=0.02)
            cbar.set_label(pile_force_label, rotation=270, labelpad=15, fontsize=10)

    # Add legend if any special states exist
    handles, labels = ax.get_legend_handles_labels()
    if handles:
        ax.legend(loc='lower right', fontsize=9, framealpha=0.9)

    return cbar_specs


def plot_reinforcement_force_profiles(fem_data, solution, figsize=(12, 8), save_png=False, dpi=300):
    """
    Plot axial force profiles along each reinforcement line as subplots.
    """
    if 'elements_1d' not in fem_data:
        print("No reinforcement elements found")
        return None, None
    
    nodes = fem_data["nodes"]
    elements_1d = fem_data["elements_1d"]
    element_materials_1d = fem_data["element_materials_1d"]
    forces_1d = solution.get("forces_1d", np.zeros(len(elements_1d)))
    t_allow = fem_data.get("t_allow_by_1d_elem", np.ones(len(elements_1d)))
    t_res = fem_data.get("t_res_by_1d_elem", np.zeros(len(elements_1d)))
    failed_1d = solution.get("failed_1d_elements", np.zeros(len(elements_1d), dtype=bool))
    
    # Group elements by reinforcement line (material ID)
    unique_lines = np.unique(element_materials_1d)
    n_lines = len(unique_lines)
    
    if n_lines == 0:
        print("No reinforcement lines found")
        return None, None
    
    # Create subplot layout
    if n_lines <= 3:
        fig, axes = plt.subplots(n_lines, 1, figsize=figsize, squeeze=False)
        axes = axes.flatten()
    else:
        rows = int(np.ceil(n_lines / 2))
        fig, axes = plt.subplots(rows, 2, figsize=figsize, squeeze=False)
        axes = axes.flatten()
    
    for line_idx, line_id in enumerate(unique_lines):
        ax = axes[line_idx]
        
        # Get elements for this line
        line_elements = np.where(element_materials_1d == line_id)[0]
        
        if len(line_elements) == 0:
            continue
        
        # Get element positions along the line
        positions = []
        forces = []
        t_allow_line = []
        t_res_line = []
        failed_line = []
        
        for elem_idx in line_elements:
            elem = elements_1d[elem_idx]
            # Use midpoint of element
            mid_point = 0.5 * (nodes[elem[0]] + nodes[elem[1]])
            # Distance along line (simplified - use x-coordinate)
            positions.append(mid_point[0])
            forces.append(forces_1d[elem_idx])
            t_allow_line.append(t_allow[elem_idx])
            t_res_line.append(t_res[elem_idx])
            failed_line.append(failed_1d[elem_idx])
        
        # Sort by position
        sorted_indices = np.argsort(positions)
        positions = np.array(positions)[sorted_indices]
        forces = np.array(forces)[sorted_indices]
        t_allow_line = np.array(t_allow_line)[sorted_indices]
        t_res_line = np.array(t_res_line)[sorted_indices]
        failed_line = np.array(failed_line)[sorted_indices]
        
        # Plot force profile
        ax.plot(positions, forces, 'b-o', linewidth=2, markersize=6, label='Tensile Force')
        ax.plot(positions, t_allow_line, 'g--', linewidth=1, label='Allowable Force')
        
        if np.any(t_res_line > 0):
            ax.plot(positions, t_res_line, 'orange', linestyle='--', linewidth=1, label='Residual Force')
        
        # Mark failed elements
        if np.any(failed_line):
            failed_positions = positions[failed_line]
            failed_forces = forces[failed_line]
            ax.scatter(failed_positions, failed_forces, color='red', s=100, marker='x', 
                      linewidth=3, label='Failed Elements', zorder=10)
        
        # Formatting. Reinforcement forces are per unit width of slope (force/length);
        # label with the declared force_per_len unit, or the bare 'Force' (byte-for-byte
        # as before) when the model declares no unit system.
        ax.set_xlabel('Position along line')
        ax.set_ylabel(_fem_cbar_label(fem_data, 'Force', 'force_per_len'))
        ax.set_title(f'Reinforcement Line {line_id} Force Profile')
        ax.grid(True, alpha=0.3)
        ax.legend()
        
        # Set y-limits to show all relevant values
        max_val = max(np.max(np.abs(forces)), np.max(t_allow_line))
        if max_val > 0:
            ax.set_ylim([-max_val * 0.1, max_val * 1.1])
    
    # Hide unused subplots
    for i in range(n_lines, len(axes)):
        axes[i].set_visible(False)
    
    plt.tight_layout()
    
    if save_png:
        filename = 'plot_reinforcement_force_profiles.png'
        plt.savefig(filename, dpi=dpi, bbox_inches='tight')
    
    return fig, axes


#: What the search figure draws each trial as. A trial the section STOOD under
#: moves the lower end of the interval up; one it did not moves the upper end
#: down. Two marks, so the closing is readable in grayscale as well as in color.
_TRIAL_STYLE = {
    True: ("#2e7d32", "o", "the section stood"),
    False: ("#c62828", "v", "the section did not stand"),
}


def ssrm_trials(record):
    """The per-trial record a strength reduction run kept, oldest first.

    ``record`` is :func:`xslope.fem.solve_ssrm`'s own result, or the meta sidecar
    a saved run was written with (:func:`xslope.fem.ssrm_run_record` carries the
    trials into it unchanged). Trials with no factor to place are dropped: a
    record is only as good as the numbers in it.
    """
    out = []
    for trial in ((record or {}).get("trials") or []):
        try:
            F = float(trial.get("F"))
        except (TypeError, ValueError):
            continue
        out.append({"F": F,
                    "role": str(trial.get("role") or "bisect"),
                    # 'stable' is the verdict the bisection acted on; a record
                    # written before the hybrid criterion existed carries only
                    # 'converged', which was then the same thing.
                    "stood": bool(trial.get("stable",
                                            trial.get("converged", False)))})
    return out


def ssrm_has_convergence_history(record):
    """True when a run's record carries enough trials to draw the search closing.

    The question asked before :func:`plot_ssrm_convergence` is called, so a caller
    counts the figures it will get right. A displacement catastrophe run keeps no
    trial record at all, and a run restored from a file written before the record
    was persisted keeps none either.
    """
    return len(ssrm_trials(record)) >= 2


def ssrm_interval_history(record):
    """``[(F_lo, F_hi), …]`` — the interval after each trial, in trial order.

    The search is replayed from the roles the trials were solved under, which is
    what :func:`xslope.fem.solve_ssrm` moves its own ends by: a lower-bound trial
    sets the low end whether or not it stood (the bracketing loop lowers it until
    it does), an upper-bound trial sets the high end, and each bisection trial
    moves the end its verdict belongs to. An end not yet established is ``None``.
    """
    lo = hi = None
    history = []
    for trial in ssrm_trials(record):
        if trial["role"] == "lower":
            lo = trial["F"]
        elif trial["role"] == "upper":
            hi = trial["F"]
        elif trial["stood"]:
            lo = trial["F"]
        else:
            hi = trial["F"]
        history.append((lo, hi))
    return history


def plot_ssrm_convergence(record, fs=None, tolerance=None, figsize=(7.4, 4.4),
                          fig=None, show_title=True, show_legend=True,
                          save_png=False, dpi=300):
    """How the strength reduction search closed on the factor of safety.

    One axes: every trial the run solved, at the factor it was solved at, marked
    by whether the section stood under it; the interval the search had narrowed
    to after each of them, shaded behind them; and the factor of safety the run
    reported, ruled across. Read left to right it is the search closing.

    ``record`` is :func:`xslope.fem.solve_ssrm`'s result or the meta sidecar a
    saved run was written with. ``fs`` and ``tolerance`` default to the record's
    own. Ask :func:`ssrm_has_convergence_history` first: a run that kept no trial
    record — a displacement catastrophe run keeps none — cannot be drawn, and
    raises rather than producing an empty axes.

    This read ``F_history`` and ``convergence_history``, which solve_ssrm has
    never emitted under any criterion, so it drew nothing on every run there has
    ever been.

    Returns ``(fig, ax)``.
    """
    trials = ssrm_trials(record)
    if len(trials) < 2:
        raise ValueError("this run kept no trial record to draw: "
                         f"{len(trials)} trial(s) with a factor on them")
    intervals = ssrm_interval_history(record)
    if fs is None:
        fs = (record or {}).get("FS")
    if tolerance is None:
        tolerance = (record or {}).get("tolerance")

    own_fig = fig is None
    if own_fig:
        fig = plt.figure(figsize=figsize)
    else:
        fig.clear()
    ax = fig.subplots(1, 1)

    steps = list(range(1, len(trials) + 1))

    # The interval, behind everything: a band from the low end to the high end,
    # held flat across the trial that produced it (step='post'), so its width at
    # any trial is the width the search had reached by then.
    band = [(n, lo, hi) for n, (lo, hi) in zip(steps, intervals)
            if lo is not None and hi is not None and hi > lo]
    if band:
        bx = [n for n, _lo, _hi in band] + [band[-1][0] + 0.5]
        blo = [lo for _n, lo, _hi in band] + [band[-1][1]]
        bhi = [hi for _n, _lo, hi in band] + [band[-1][2]]
        ax.fill_between(bx, blo, bhi, step="post", color="#7f8c9a", alpha=0.18,
                        linewidth=0, label="the interval still open")

    # The trials themselves, in the order they were solved.
    ax.plot(steps, [t["F"] for t in trials], "-", color="0.55", lw=1.0,
            zorder=2)
    for stood in (True, False):
        color, marker, label = _TRIAL_STYLE[stood]
        xs = [n for n, t in zip(steps, trials) if t["stood"] is stood]
        ys = [t["F"] for t in trials if t["stood"] is stood]
        if xs:
            ax.plot(xs, ys, marker, color=color, ms=6, ls="none", label=label,
                    zorder=3)

    if fs is not None:
        note = f"FS = {float(fs):.3f}"
        if tolerance is not None:
            note += f" ± {float(tolerance) / 2:g}"
        ax.axhline(float(fs), color="#1f4e79", ls="--", lw=1.6, zorder=4,
                   label=note)

    ax.set_xlabel("Trial, in the order it was solved")
    ax.set_ylabel("Strength reduction factor $F$")
    ax.set_xlim(0.5, len(trials) + 0.5)
    ax.set_xticks(steps)
    ax.grid(alpha=0.25)
    if show_title:
        ax.set_title("Strength reduction search", fontsize=11)
    if show_legend:
        ax.legend(loc="best", fontsize=8.5, framealpha=0.9)

    try:
        fig.tight_layout()
    except Exception:
        pass
    if save_png:
        fig.savefig('plot_ssrm_convergence.png', dpi=dpi, bbox_inches='tight')
    return fig, ax


def plot_strain_contours(ax, fem_data, solution, show_mesh=True, show_reinforcement=True,
                        cbar_shrink=0.8, cbar_labelpad=20, label_elements=False):
    """
    Plot von Mises equivalent strain contours computed from total strains.
    """
    nodes = fem_data["nodes"]
    elements = fem_data["elements"]
    element_types = fem_data["element_types"]
    strains = solution.get("strains", np.zeros((len(elements), 4)))
    
    if strains.shape[1] < 3:
        print("Warning: Strain data not available or incomplete")
        return
    
    # Calculate equivalent strain (von Mises equivalent strain)
    # For plane strain: equiv_strain = sqrt(2/3) * sqrt(eps_x^2 + eps_y^2 + eps_x*eps_y + 3/4*gamma_xy^2)
    eps_x = strains[:, 0]
    eps_y = strains[:, 1]
    gamma_xy = strains[:, 2]
    
    equiv_strain = np.sqrt((2/3) * (eps_x**2 + eps_y**2 + eps_x*eps_y + 0.75*gamma_xy**2))
    
    # Plot contours
    _plot_element_contours(ax, fem_data, equiv_strain, 'Equivalent Strain', 
                          show_mesh, show_reinforcement, cbar_shrink, cbar_labelpad, label_elements)


def plot_shear_strain_contours(ax, fem_data, solution, show_mesh=True, show_reinforcement=True,
                              cbar_shrink=0.8, cbar_labelpad=20, label_elements=False, cmap=None,
                              single_panel=False, vmin=None, vmax=None):
    """
    Plot viscoplastic max shear strain contours.

    Uses accumulated viscoplastic strains from the solution (vp_shear_strain key).
    Falls back to total shear strain if VP data is not available.

    When ``single_panel`` is True the inline colorbar is suppressed and the contour
    mappable is returned so the caller can place the colorbar manually (sized to the
    plot box). Returns the mappable (or None).

    ``vmin`` / ``vmax`` pin the contour range, for a caller drawing this field at
    more than one state on one scale (:func:`shared_panel_scales`).
    """
    nodes = fem_data["nodes"]
    elements = fem_data["elements"]
    element_types = fem_data["element_types"]

    vp_shear_strain = solution.get("vp_shear_strain", None)
    if vp_shear_strain is None:
        # Fallback to total shear strain if VP not available
        strains = solution.get("strains", np.zeros((len(elements), 4)))
        if strains.shape[1] >= 4:
            vp_shear_strain = strains[:, 3]
        else:
            print("Warning: Shear strain data not available")
            return None, []

    # show_mesh draws the element edges over the contours (reinforcement is drawn
    # separately below with force-based coloring, so it stays False here).
    mappable = _plot_nodal_contours(ax, fem_data, vp_shear_strain, SHEAR_STRAIN_LABEL,
                        show_mesh, False, cbar_shrink, cbar_labelpad,
                        colormap=cmap or 'coolwarm', label_elements=label_elements,
                        draw_cbar=not single_panel, vmin=vmin, vmax=vmax)

    # Draw reinforcement with force-based coloring. When the strain colorbar is
    # deferred (single_panel — placed by plot_fem_results via make_axes_locatable),
    # defer the force colorbar the same way and hand its spec back, so the two bars
    # get separate full-height slots instead of colliding in one.
    reinf_cbar_specs = []
    if show_reinforcement and 'elements_1d' in fem_data:
        reinf_cbar_specs = plot_reinforcement_forces(
            ax, fem_data, solution, draw_cbar=not single_panel)

    F = solution.get("F", None)
    at_failure = solution.get("_at_failure", False)
    title = SHEAR_STRAIN_LABEL
    if at_failure:
        title += ' at Failure'
    # at_failure when plot_fem_results routed this panel onto the at-failure field
    # (field_state='failure', the default): leads with FS, matching the
    # deformation/displace_vector panels so the figure tells one story.
    title = _fs_title(title, F, solution.get("_ssrm_fs"), at_failure=at_failure)
    ax.set_title(title, fontsize=12, pad=15)
    return mappable, reinf_cbar_specs


def plot_yield_function_contours(ax, fem_data, solution, show_mesh=True, show_reinforcement=True, 
                                cbar_shrink=0.8, cbar_labelpad=20, label_elements=False):
    """
    Plot yield function values (Mohr-Coulomb failure criterion).
    Positive values indicate yielding/failure, negative values indicate elastic state.
    """
    nodes = fem_data["nodes"]
    elements = fem_data["elements"]
    element_types = fem_data["element_types"]
    yield_function = solution.get("yield_function", None)
    
    if yield_function is None:
        print("Warning: Yield function data not available in solution")
        # Create dummy data
        yield_function = np.zeros(len(elements))
    
    # Create custom colormap for yield function visualization
    # Strong blue for very negative (very safe), white near zero, red for positive (yielding)
    from matplotlib.colors import LinearSegmentedColormap
    
    # Define color transitions for yield function
    # F < 0: shades of blue/green (elastic/safe)
    # F = 0: white/light gray (critical)
    # F > 0: shades of red (yielding/plastic)
    colors_below = ['#0000FF', '#0066FF', '#00AAFF', '#00DDDD', '#CCCCCC']  # Blue to gray
    colors_above = ['#FFCCCC', '#FF9999', '#FF6666', '#FF3333', '#FF0000', '#CC0000']  # Light red to dark red
    
    # Create custom colormap with sharp transition at F=0
    n_bins = 256
    n_below = int(n_bins * 0.7)  # 70% for negative values
    n_above = n_bins - n_below   # 30% for positive values
    
    from matplotlib.colors import ListedColormap
    colors_below_interp = plt.cm.Blues_r(np.linspace(0.2, 0.9, n_below))
    colors_above_interp = plt.cm.Reds(np.linspace(0.3, 1.0, n_above))
    colors_all = np.vstack([colors_below_interp, colors_above_interp])
    cmap_yield = ListedColormap(colors_all)
    
    # Set visualization bounds - asymmetric to focus on near-yield region
    vmin = -200  # Cap negative values for better contrast
    vmax = 50    # Positive values are more important
    
    # Plot each element as a colored patch
    from matplotlib.collections import PatchCollection
    from matplotlib.patches import Polygon
    patches_list = []
    values_list = []
    
    for i, elem in enumerate(elements):
        elem_type = element_types[i]
        if elem_type == 3:  # Triangle
            coords = nodes[elem[:3]]
        elif elem_type == 4:  # Quad
            coords = nodes[elem[:4]]
        elif elem_type == 6:  # 6-node triangle - use corner nodes
            coords = nodes[elem[:3]]
        elif elem_type in [8, 9]:  # 8 or 9-node quad - use corner nodes
            coords = nodes[elem[:4]]
        else:
            continue
            
        patch = Polygon(coords, closed=True)
        patches_list.append(patch)
        # Clip values for visualization
        values_list.append(np.clip(yield_function[i], vmin, vmax))
    
    if patches_list:
        # Element edges are baked into the patch collection; gate them on show_mesh
        # (the caller's mesh_on_fields opt-in) so the yield field reads clean by
        # default like the other filled-field panels.
        p = PatchCollection(patches_list, alpha=0.9,
                            edgecolors='gray' if show_mesh else 'none', linewidths=0.3)
        p.set_array(np.array(values_list))
        p.set_cmap(cmap_yield)
        p.set_clim(vmin, vmax)
        ax.add_collection(p)
        
        # Add colorbar with custom ticks
        cbar = ax.figure.colorbar(p, ax=ax, shrink=cbar_shrink)
        cbar.set_label('Yield Function F', rotation=270, labelpad=cbar_labelpad)
        
        # Set custom ticks to highlight key values
        tick_values = [-200, -100, -50, -20, -10, -5, 0, 5, 10, 20, 50]
        tick_labels = ['-200', '-100', '-50', '-20', '-10', '-5', '0', '5', '10', '20', '50']
        # Filter ticks to those within bounds
        valid_ticks = [(v, l) for v, l in zip(tick_values, tick_labels) if vmin <= v <= vmax]
        if valid_ticks:
            tick_values, tick_labels = zip(*valid_ticks)
            cbar.set_ticks(tick_values)
            cbar.set_ticklabels(tick_labels)
        
        # Add a line at F=0
        cbar.ax.axhline(y=0, color='black', linewidth=2)
    
    # Add yield function values as text on elements (if requested or for yielding elements)
    for i, elem in enumerate(elements):
        elem_type = element_types[i]
        
        # Get element centroid
        if elem_type == 3:  # Triangle
            elem_nodes = nodes[elem[:3]]
        elif elem_type == 4:  # Quad
            elem_nodes = nodes[elem[:4]]
        elif elem_type == 6:  # 6-node triangle - use corner nodes
            elem_nodes = nodes[elem[:3]]
        elif elem_type in [8, 9]:  # 8 or 9-node quad - use corner nodes
            elem_nodes = nodes[elem[:4]]
        else:
            continue
            
        centroid = np.mean(elem_nodes, axis=0)
        
        # Show values for elements that are close to yielding or already yielding
        # or if label_elements is True
        f_val = yield_function[i]
        
        if label_elements or f_val > -50:  # Show if requested or if close to yielding
            # Format the number based on magnitude
            if abs(f_val) < 10:
                text = f'{f_val:.1f}'
            else:
                text = f'{f_val:.0f}'
            
            # Choose text color based on value
            if f_val > 0:
                color = 'white'  # White on red background
                fontweight = 'bold'
            elif f_val > -10:
                color = 'black'  # Black on light background
                fontweight = 'normal'
            else:
                color = 'white'  # White on blue background
                fontweight = 'normal'
            
            # Only show for elements near yield or if explicitly requested
            if label_elements or f_val > -30:
                ax.text(centroid[0], centroid[1], text,
                       ha='center', va='center', fontsize=5,
                       color=color, fontweight=fontweight, alpha=0.8)
    
    # Highlight yielding elements with thick red border
    for i, elem in enumerate(elements):
        if yield_function[i] > 0:
            elem_type = element_types[i]
            if elem_type == 3:  # Triangle
                coords = nodes[elem[:3]]
            elif elem_type == 4:  # Quad
                coords = nodes[elem[:4]]
            elif elem_type == 6:  # 6-node triangle - use corner nodes
                coords = nodes[elem[:3]]
            elif elem_type in [8, 9]:  # 8 or 9-node quad - use corner nodes
                coords = nodes[elem[:4]]
            else:
                continue
            
            # Close the polygon
            coords = np.vstack([coords, coords[0]])
            ax.plot(coords[:, 0], coords[:, 1], 'k-', linewidth=2.5, alpha=1.0)  # Black border for yielding elements
    
    # Add reinforcement if requested
    if show_reinforcement and 'elements_1d' in fem_data:
        plot_reinforcement_lines(ax, fem_data, solution)
    
    # Add title indicating yield state
    title = 'Yield Function (Red: F>0 Yielding/Plastic, Blue: F<0 Elastic)'
    if solution.get("_at_failure", False):
        title += ' at Failure'
    ax.set_title(title, fontsize=12, pad=15)
    
    # Add statistics to the plot
    n_yielding = np.sum(yield_function > 0)
    n_total = len(yield_function)
    n_critical = np.sum((yield_function > -10) & (yield_function <= 0))  # Near yielding
    
    stats_text = f'Yielding: {n_yielding}/{n_total} elements\n'
    stats_text += f'Critical (F>-10): {n_critical} elements\n'
    stats_text += f'Max F: {np.max(yield_function):.1f}\n'
    stats_text += f'Min F: {np.min(yield_function):.1f}'
    
    ax.text(0.02, 0.98, stats_text,
            transform=ax.transAxes, fontsize=9, verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.7))


def _plot_element_contours(ax, fem_data, values, label, show_mesh=True, show_reinforcement=True,
                          cbar_shrink=0.8, cbar_labelpad=20, label_elements=False, colormap='viridis'):
    """
    Plot element-based scalar data as colored patches (one color per element).
    """
    nodes = fem_data["nodes"]
    elements = fem_data["elements"]
    element_types = fem_data["element_types"]
    
    # For element-based values, we need to interpolate to nodes or use a different approach
    # Let's use a simpler approach: plot each element as a colored patch
    
    # Create contour plot by directly coloring elements
    if np.max(values) > np.min(values):  # Only plot if there's variation
        # Normalize values for colormap
        vmin, vmax = np.min(values), np.max(values)
        norm = plt.Normalize(vmin=vmin, vmax=vmax)
        cmap = plt.get_cmap(colormap)
        
        # Plot each element as colored patch
        for i, elem in enumerate(elements):
            elem_type = element_types[i]
            color = cmap(norm(values[i]))
            
            if elem_type == 3:  # Triangle
                coords = nodes[elem[:3]]
                triangle = plt.Polygon(coords, facecolor=color, edgecolor='none', alpha=0.8)
                ax.add_patch(triangle)
            elif elem_type == 4:  # Quad
                coords = nodes[elem[:4]]
                quad = plt.Polygon(coords, facecolor=color, edgecolor='none', alpha=0.8)
                ax.add_patch(quad)
            elif elem_type == 6:  # 6-node triangle - use corner nodes
                coords = nodes[elem[:3]]
                triangle = plt.Polygon(coords, facecolor=color, edgecolor='none', alpha=0.8)
                ax.add_patch(triangle)
            elif elem_type in [8, 9]:  # 8-node or 9-node quad - use corner nodes
                coords = nodes[elem[:4]]
                quad = plt.Polygon(coords, facecolor=color, edgecolor='none', alpha=0.8)
                ax.add_patch(quad)
        
        # Create colorbar using a ScalarMappable
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        cbar = ax.figure.colorbar(sm, ax=ax, pad=0.05)
        cbar.set_label(label, rotation=270, labelpad=cbar_labelpad)
        adaptive_colorbar_ticks(ax.figure, cbar)
    else:
        # Uniform values - just color all elements the same
        for i, elem in enumerate(elements):
            elem_type = element_types[i]
            if elem_type == 3:  # Triangle
                coords = nodes[elem[:3]]
                triangle = plt.Polygon(coords, facecolor='lightblue', edgecolor='none', alpha=0.7)
                ax.add_patch(triangle)
            elif elem_type == 4:  # Quad
                coords = nodes[elem[:4]]
                quad = plt.Polygon(coords, facecolor='lightblue', edgecolor='none', alpha=0.7)
                ax.add_patch(quad)
            elif elem_type == 6:  # 6-node triangle - use corner nodes
                coords = nodes[elem[:3]]
                triangle = plt.Polygon(coords, facecolor='lightblue', edgecolor='none', alpha=0.7)
                ax.add_patch(triangle)
            elif elem_type in [8, 9]:  # 8-node or 9-node quad - use corner nodes
                coords = nodes[elem[:4]]
                quad = plt.Polygon(coords, facecolor='lightblue', edgecolor='none', alpha=0.7)
                ax.add_patch(quad)
    
    # Overlay mesh if requested
    if show_mesh:
        for i, elem in enumerate(elements):
            elem_type = element_types[i]
            if elem_type == 3:  # Triangle
                coords = nodes[elem[:3]]
                triangle = plt.Polygon(coords, fill=False, edgecolor='black', linewidth=0.5, alpha=0.7)
                ax.add_patch(triangle)
            elif elem_type == 4:  # Quad
                coords = nodes[elem[:4]]
                quad = plt.Polygon(coords, fill=False, edgecolor='black', linewidth=0.5, alpha=0.7)
                ax.add_patch(quad)
            elif elem_type == 6:  # 6-node triangle - use corner nodes
                coords = nodes[elem[:3]]
                triangle = plt.Polygon(coords, fill=False, edgecolor='black', linewidth=0.5, alpha=0.7)
                ax.add_patch(triangle)
            elif elem_type in [8, 9]:  # 8-node or 9-node quad - use corner nodes
                coords = nodes[elem[:4]]
                quad = plt.Polygon(coords, fill=False, edgecolor='black', linewidth=0.5, alpha=0.7)
                ax.add_patch(quad)
    
    # Add reinforcement if requested
    if show_reinforcement:
        elements_1d = fem_data.get("elements_1d", np.array([]).reshape(0, 3))
        if len(elements_1d) > 0:
            for elem in elements_1d:
                if len(elem) >= 2:
                    x_coords = [nodes[elem[0], 0], nodes[elem[1], 0]]
                    y_coords = [nodes[elem[0], 1], nodes[elem[1], 1]]
                    ax.plot(x_coords, y_coords, 'r-', linewidth=2, alpha=0.8)
    
    # Add element labels if requested
    if label_elements:
        _add_element_labels(ax, fem_data)
    
    ax.set_aspect('equal')


def _nodal_average(fem_data, element_values):
    """Element values averaged onto the nodes — the array the filled-contour
    panels actually contour.

    Split out of :func:`_plot_nodal_contours` so a caller pinning two panels to
    one contour range measures the range on the SAME numbers the panels draw:
    averaging onto nodes narrows an element field, and a range taken off the
    element array would wash the fill out.
    """
    nodes = fem_data["nodes"]
    elements = fem_data["elements"]
    element_types = fem_data["element_types"]

    nodal_values = np.zeros(len(nodes))
    node_counts = np.zeros(len(nodes))

    for i, elem in enumerate(elements):
        elem_type = element_types[i]
        elem_nodes = elem[:elem_type] if elem_type <= len(elem) else elem

        # Add this element's value to all its nodes
        for node_id in elem_nodes:
            if node_id < len(nodes):
                nodal_values[node_id] += element_values[i]
                node_counts[node_id] += 1

    # Average values at nodes (avoid division by zero)
    valid_nodes = node_counts > 0
    nodal_values[valid_nodes] /= node_counts[valid_nodes]
    return nodal_values


def _plot_nodal_contours(ax, fem_data, element_values, label, show_mesh=True, show_reinforcement=True,
                        cbar_shrink=0.8, cbar_labelpad=20, colormap='viridis', label_elements=False,
                        draw_cbar=True, vmin=None, vmax=None):
    """
    Plot smooth filled contours by averaging element values to nodes and triangulating.

    Returns the contour mappable (or None if the field was uniform / empty) so a
    caller can place the colorbar itself; when ``draw_cbar`` is True (default) the
    colorbar is drawn inline as before.

    ``vmin`` / ``vmax`` pin the contour levels to a range the caller holds several
    panels to (see :func:`shared_panel_scales`); None (default) scales each panel
    to its own field.
    """
    nodes = fem_data["nodes"]
    elements = fem_data["elements"]
    element_types = fem_data["element_types"]

    nodal_values = _nodal_average(fem_data, element_values)

    # Create triangulation for smooth contouring
    triangles = []
    for i, elem in enumerate(elements):
        elem_type = element_types[i]
        if elem_type == 3:  # Triangle
            triangles.append([elem[0], elem[1], elem[2]])
        elif elem_type == 4:  # Quad - split into triangles
            triangles.append([elem[0], elem[1], elem[2]])
            triangles.append([elem[0], elem[2], elem[3]])
        elif elem_type == 6:  # 6-node triangle - use corner nodes
            triangles.append([elem[0], elem[1], elem[2]])
        elif elem_type in [8, 9]:  # 8-node or 9-node quad - use corner nodes
            triangles.append([elem[0], elem[1], elem[2]])
            triangles.append([elem[0], elem[2], elem[3]])
    
    if not triangles:
        print("No valid elements for contouring")
        return None

    import matplotlib.tri as tri
    triangles = np.array(triangles)
    
    # Create triangulation
    triang = tri.Triangulation(nodes[:, 0], nodes[:, 1], triangles)
    
    # Create smooth contour plot
    mappable = None
    # The range the levels span: the caller's where it pinned one, this field's
    # otherwise. A pinned range is used even where this field is flat inside it —
    # that IS the fact the pair is drawn to show.
    low = np.min(nodal_values) if vmin is None else float(vmin)
    high = np.max(nodal_values) if vmax is None else float(vmax)
    # The colorbar's out-of-range arrows are drawn only where this field really
    # does run past the pinned range — a range that spans the fields it was
    # resolved from clips nothing, and an arrow on it would say it had.
    extend = 'neither'
    if vmin is not None or vmax is not None:
        below = float(np.min(nodal_values)) < low - 1e-12
        above = float(np.max(nodal_values)) > high + 1e-12
        extend = ('both' if below and above else
                  'min' if below else 'max' if above else 'neither')
    if high > low:  # Only plot if there's variation
        levels = np.linspace(low, high, 20)
        cs = ax.tricontourf(triang, nodal_values, levels=levels, cmap=colormap,
                            extend=extend)
        # DXF layer named after the plotted quantity (e.g. "Viscoplastic shear strain").
        cs.set_gid((label or 'CONTOURS').upper().replace(' ', '_') + '_CONTOURS')
        mappable = cs

        # Add colorbar (unless the caller will place it itself — single-panel case).
        # Full frame-height bar (constrained-layout native: no shrink) with
        # height-adaptive, round-valued ticks so the labels never stack.
        if draw_cbar:
            cbar = ax.figure.colorbar(cs, ax=ax, pad=0.02)
            cbar.set_label(label, rotation=270, labelpad=20)
            adaptive_colorbar_ticks(ax.figure, cbar)
    else:
        # Uniform values - just color all elements the same
        uniform_color = plt.get_cmap(colormap)(0.5)
        for triangle_nodes in triangles:
            coords = nodes[triangle_nodes]
            triangle = plt.Polygon(coords, facecolor=uniform_color, edgecolor='none', alpha=0.8)
            ax.add_patch(triangle)

    # Overlay mesh element edges if requested (light gray, matching the element
    # edges drawn on the deformation / displacement-vector plots).
    if show_mesh:
        for i, elem in enumerate(elements):
            elem_type = element_types[i]
            if elem_type == 3:  # Triangle
                coords = nodes[elem[:3]]
                triangle = plt.Polygon(coords, fill=False, edgecolor='lightgray', linewidth=0.5, alpha=0.6)
                ax.add_patch(triangle)
            elif elem_type == 4:  # Quad
                coords = nodes[elem[:4]]
                quad = plt.Polygon(coords, fill=False, edgecolor='lightgray', linewidth=0.5, alpha=0.6)
                ax.add_patch(quad)
            elif elem_type == 6:  # 6-node triangle - use corner nodes
                coords = nodes[elem[:3]]
                triangle = plt.Polygon(coords, fill=False, edgecolor='lightgray', linewidth=0.5, alpha=0.6)
                ax.add_patch(triangle)
            elif elem_type in [8, 9]:  # 8-node or 9-node quad - use corner nodes
                coords = nodes[elem[:4]]
                quad = plt.Polygon(coords, fill=False, edgecolor='lightgray', linewidth=0.5, alpha=0.6)
                ax.add_patch(quad)

    # Add reinforcement if requested
    if show_reinforcement:
        elements_1d = fem_data.get("elements_1d", np.array([]).reshape(0, 3))
        if len(elements_1d) > 0:
            for elem in elements_1d:
                if len(elem) >= 2:
                    x_coords = [nodes[elem[0], 0], nodes[elem[1], 0]]
                    y_coords = [nodes[elem[0], 1], nodes[elem[1], 1]]
                    ax.plot(x_coords, y_coords, 'r-', linewidth=2, alpha=0.8)

    # Add element labels if requested
    if label_elements:
        _add_element_labels(ax, fem_data)

    ax.set_aspect('equal')
    return mappable