import numpy as np
import gmsh
import matplotlib.pyplot as plt
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import reverse_cuthill_mckee

def build_tri_mesh(points, target_size):

    gmsh.initialize()
    gmsh.model.add("paving_mesh")

    # Add points to Gmsh with target size
    point_tags = [gmsh.model.geo.addPoint(x, y, 0, target_size) for x, y in points]

    # Create lines between consecutive points (wrap last to first)
    line_tags = [
        gmsh.model.geo.addLine(point_tags[i], point_tags[(i + 1) % len(point_tags)])
        for i in range(len(point_tags))
    ]

    # Create loop and surface
    loop = gmsh.model.geo.addCurveLoop(line_tags)
    gmsh.model.geo.addPlaneSurface([loop])

    gmsh.model.geo.synchronize()
    gmsh.model.mesh.generate(2)

    # Extract nodes and triangles
    node_tags, node_coords, _ = gmsh.model.mesh.getNodes()
    nodes = np.array(node_coords).reshape(-1, 3)[:, :2]
    elements = gmsh.model.mesh.getElementsByType(2)[1].reshape(-1, 3) - 1  # zero-based

    gmsh.finalize()
    return nodes, elements

def build_tri_mesh_with_regions(polygons, region_ids, target_size):
    """
    Parameters:
        polygons    : List of lists of (x, y) tuples
        region_ids  : List of material IDs (one per polygon)
        target_size : Desired triangle size

    Returns:
        nodes       : np.ndarray of node coordinates (n_nodes, 2)
        elements    : np.ndarray of triangle vertex indices (n_elements, 3)
        mat_ids     : np.ndarray of material ID for each triangle (n_elements,)
    """
    import gmsh
    import numpy as np

    gmsh.initialize()
    gmsh.model.add("multi_region_mesh")
    point_map = {}  # maps (x, y) to Gmsh point tag
    point_tag_counter = 1

    def add_point(x, y):
        nonlocal point_tag_counter
        key = (x, y)
        if key not in point_map:
            tag = gmsh.model.geo.addPoint(x, y, 0, target_size)
            point_map[key] = tag
        return point_map[key]

    surface_tags = []
    for poly_pts, region_id in zip(polygons, region_ids):
        pt_tags = [add_point(x, y) for x, y in poly_pts]
        line_tags = [gmsh.model.geo.addLine(pt_tags[i], pt_tags[(i + 1) % len(pt_tags)])
                     for i in range(len(pt_tags))]
        loop = gmsh.model.geo.addCurveLoop(line_tags)
        surface = gmsh.model.geo.addPlaneSurface([loop])
        surface_tags.append(surface)
        gmsh.model.addPhysicalGroup(2, [surface], region_id)

    gmsh.model.geo.synchronize()
    gmsh.model.mesh.generate(2)

    # Get nodes
    node_tags, coords, _ = gmsh.model.mesh.getNodes()
    nodes = np.array(coords).reshape(-1, 3)[:, :2]

    # Get triangle elements and material IDs
    elements = []
    mat_ids = []

    for region_id in region_ids:
        phys_group = gmsh.model.getEntitiesForPhysicalGroup(2, region_id)
        for surface in phys_group:
            elem_type, elem_tags, node_tags = gmsh.model.mesh.getElements(2, surface)
            tri_conn = np.array(node_tags[0]).reshape(-1, 3) - 1  # zero-based
            elements.append(tri_conn)
            mat_ids.extend([region_id] * len(tri_conn))

    elements = np.vstack(elements)
    mat_ids = np.array(mat_ids)

    gmsh.finalize()
    return nodes, elements, mat_ids


def reorder_mesh(nodes, elements):
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import reverse_cuthill_mckee

    n_nodes = nodes.shape[0]
    rows, cols = [], []

    for tri in elements:
        for i in range(3):
            for j in range(i + 1, 3):
                rows.append(tri[i])
                cols.append(tri[j])
                rows.append(tri[j])
                cols.append(tri[i])

    adj = coo_matrix((np.ones(len(rows)), (rows, cols)), shape=(n_nodes, n_nodes))

    def bandwidth(el):
        if len(el) == 0:
            return 0
        return max(abs(int(i) - int(j)) for tri in el for i in tri for j in tri)

    bw_before = bandwidth(elements)

    # RCM reordering
    perm = reverse_cuthill_mckee(adj.tocsr())
    inverse_perm = np.zeros_like(perm)
    inverse_perm[perm] = np.arange(len(perm))

    new_nodes = nodes[perm]
    new_elements = np.array([[inverse_perm[i] for i in tri] for tri in elements])

    # Bandwidth after
    bw_after = bandwidth(new_elements)

    print(f"Bandwidth before reordering: {bw_before}")
    print(f"Bandwidth after reordering:  {bw_after}")

    return new_nodes, new_elements


def plot_mesh(nodes, elements):
    for tri in elements:
        tri_pts = nodes[np.append(tri, tri[0])]
        xs, ys = tri_pts[:, 0], tri_pts[:, 1]
        plt.plot(xs, ys, 'k-')

    plt.gca().set_aspect('equal')
    plt.show()


# Simple plot showing material regions
def plot_mesh_with_materials(nodes, elements, mat_ids):
    import matplotlib.pyplot as plt
    from matplotlib.cm import tab10

    for tri, mid in zip(elements, mat_ids):
        pts = nodes[np.append(tri, tri[0])]
        xs, ys = pts[:, 0], pts[:, 1]
        plt.fill(xs, ys, color=tab10(mid % 10), edgecolor='k', alpha=0.6)

    plt.gca().set_aspect('equal')
    plt.title("Two-Region Mesh with Material IDs")
    plt.show()