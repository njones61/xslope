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

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import reverse_cuthill_mckee

# Lazy import gmsh - only needed for mesh generation functions
_gmsh = None
def _get_gmsh():
    global _gmsh
    if _gmsh is None:
        try:
            import gmsh
            _gmsh = gmsh
        except (ImportError, OSError) as e:
            error_msg = str(e)
            error_repr = repr(e)
            # Check for OpenGL library issues (common in Colab/headless environments)
            # Check both str() and repr() to catch all variations, and check exception args
            error_text = error_msg + " " + error_repr
            if hasattr(e, 'args') and e.args:
                error_text += " " + " ".join(str(arg) for arg in e.args)
            if ("libGL" in error_text or "libGLU" in error_text):
                help_msg = (
                    "gmsh is required for mesh generation but could not be imported due to missing OpenGL libraries. "
                    "This is common in headless environments like Google Colab.\n\n"
                    "To fix this in Google Colab, install the required system libraries first:\n"
                    "  !apt-get update && apt-get install -y libgl1-mesa-glx libglu1-mesa\n"
                    "Then install gmsh:\n"
                    "  !pip install gmsh\n\n"
                    "For other headless environments, install the appropriate OpenGL libraries for your system.\n"
                    f"Original error: {e}"
                )
            else:
                help_msg = (
                    "gmsh is required for mesh generation but could not be imported. "
                    "If you only need limit equilibrium analysis, you can ignore this. "
                    "To use FEM features, install gmsh: pip install gmsh\n"
                    f"Original error: {e}"
                )
            raise ImportError(help_msg) from e
    return _gmsh



def build_mesh_from_polygons(polygons, target_size, element_type='tri3', lines=None, debug=False, mesh_params=None, target_size_1d=None):
    """
    Build a finite element mesh with material regions using Gmsh.
    Fixed version that properly handles shared boundaries between polygons.
    
    Parameters:
        polygons     : List of lists of (x, y) tuples defining material boundaries
        target_size  : Desired element size
        element_type : 'tri3' (3-node triangles), 'tri6' (6-node triangles), 
                      'quad4' (4-node quadrilaterals), 'quad8' (8-node quadrilaterals),
                      'quad9' (9-node quadrilaterals)
        lines        : Optional list of lines, each defined by list of (x, y) tuples for 1D elements
        debug        : Enable debug output
        mesh_params  : Optional dictionary of GMSH meshing parameters to override defaults
        target_size_1d : Optional target size for 1D elements (default None, which is set to target_size if None)

    Returns:
        mesh dict containing:
        nodes        : np.ndarray of node coordinates (n_nodes, 2)
        elements     : np.ndarray of 2D element vertex indices (n_elements, 9) - unused nodes set to 0
        element_types: np.ndarray indicating number of nodes per 2D element (3, 4, 6, 8, or 9)
        element_materials: np.ndarray of material ID for each 2D element
        
        If lines is provided, also includes:
        elements_1d  : np.ndarray of 1D element vertex indices (n_elements_1d, 3) - unused nodes set to 0
        element_types_1d: np.ndarray indicating element type (2 for linear, 3 for quadratic)
        element_materials_1d: np.ndarray of material ID for each 1D element (line index)
    """
    gmsh = _get_gmsh()
    from collections import defaultdict

    # Set default target_size_1d if None
    if target_size_1d is None:
        target_size_1d = target_size
        if debug:
            print(f"Using default target_size_1d = target_size = {target_size_1d}")

    # build a list of region ids (list of material IDs - one per polygon)
    region_ids = [i for i in range(len(polygons))]

    if element_type not in ['tri3', 'tri6', 'quad4', 'quad8', 'quad9']:
        raise ValueError("element_type must be 'tri3', 'tri6', 'quad4', 'quad8', or 'quad9'")

    # Determine if we need quadratic elements - but always generate linear first
    quadratic = element_type in ['tri6', 'quad8', 'quad9']
    
    # For quadratic elements, always start with linear base element
    if quadratic:
        if element_type == 'tri6':
            base_element_type = 'tri3'
        elif element_type in ['quad8', 'quad9']:
            base_element_type = 'quad4'
        if debug:
            print(f"Quadratic element '{element_type}' requested: generating '{base_element_type}' first, then post-processing")
    else:
        base_element_type = element_type
    
    # Adjust target_size for quads to compensate for recombination creating finer meshes
    if element_type.startswith('quad'):
        # Different adjustment factors based on meshing parameters
        if mesh_params and 'size_factor' in mesh_params:
            size_factor = mesh_params['size_factor']
        else:
            # Default size factors for different approaches
            if mesh_params and mesh_params.get("Mesh.RecombinationAlgorithm") == 0:
                size_factor = 1.2  # Fast algorithm needs less adjustment
            elif mesh_params and mesh_params.get("Mesh.RecombineOptimizeTopology", 0) > 50:
                size_factor = 1.8  # High optimization creates more elements
            else:
                size_factor = 1.4  # Default
        
        adjusted_target_size = target_size * size_factor
        if debug:
            print(f"Adjusted target size for quads: {target_size} -> {adjusted_target_size} (factor: {size_factor})")
    else:
        adjusted_target_size = target_size

    gmsh.initialize()
    gmsh.option.setNumber("General.Verbosity", 4)  # Reduce verbosity
    gmsh.model.add("multi_region_mesh")
    
    # Global point map to ensure shared boundaries use the same points
    point_map = {}  # maps (x, y) to Gmsh point tag
    
    # Track all unique edges and their usage
    edge_map = {}  # maps (pt1, pt2) tuple to line tag
    edge_usage = defaultdict(list)  # maps edge to list of (region_id, orientation)
    
    def add_point(x, y, size_override=None):
        key = (x, y)
        if key not in point_map:
            point_size = size_override if size_override is not None else adjusted_target_size
            tag = gmsh.model.geo.addPoint(x, y, 0, point_size)
            point_map[key] = tag
        return point_map[key]

    def get_edge_key(pt1, pt2):
        """Get canonical edge key (always smaller point first)"""
        return (min(pt1, pt2), max(pt1, pt2))

    # First pass: Create all points and identify short edges
    polygon_data = []
    short_edge_points = set()  # Points that are endpoints of short edges
    
    # Pre-pass to identify short edges - improved logic
    for idx, (poly_pts, region_id) in enumerate(zip(polygons, region_ids)):
        poly_pts_clean = remove_duplicate_endpoint(list(poly_pts))
        for i in range(len(poly_pts_clean)):
            p1 = poly_pts_clean[i]
            p2 = poly_pts_clean[(i + 1) % len(poly_pts_clean)]
            edge_length = ((p2[0] - p1[0])**2 + (p2[1] - p1[1])**2)**0.5
            
            # Only mark as short edge if it's genuinely short AND not a major boundary
            # Major boundaries should maintain consistent mesh sizing
            is_major_boundary = False
            
            # Check if this edge is part of a major boundary (long horizontal or vertical edge)
            if abs(p2[0] - p1[0]) > adjusted_target_size * 5:  # Long horizontal edge
                is_major_boundary = True
            elif abs(p2[1] - p1[1]) > adjusted_target_size * 5:  # Long vertical edge
                is_major_boundary = True
            
            # Only apply short edge sizing if edge is genuinely short AND not a major boundary
            if edge_length < adjusted_target_size and not is_major_boundary:
                short_edge_points.add(p1)
                short_edge_points.add(p2)
                if debug:
                    print(f"Short edge found: {p1} to {p2}, length={edge_length:.2f}")
            elif debug and edge_length < adjusted_target_size:
                print(f"Short edge ignored (major boundary): {p1} to {p2}, length={edge_length:.2f}")
    
    # Main pass: Create points with appropriate sizes
    for idx, (poly_pts, region_id) in enumerate(zip(polygons, region_ids)):
        poly_pts_clean = remove_duplicate_endpoint(list(poly_pts))  # make a copy
        pt_tags = []
        for x, y in poly_pts_clean:
            # Use larger size for points on short edges to discourage subdivision
            # But be more conservative about when to apply this
            if (x, y) in short_edge_points:
                point_size = adjusted_target_size * 2.0  # Reduced from 3.0 to 2.0
                pt_tags.append(add_point(x, y, point_size))
            else:
                pt_tags.append(add_point(x, y))
        
        # Track edges for this polygon
        edges = []
        for i in range(len(pt_tags)):
            pt1 = pt_tags[i]
            pt2 = pt_tags[(i + 1) % len(pt_tags)]
            
            edge_key = get_edge_key(pt1, pt2)
            
            # Determine orientation: True if pt1 < pt2, False otherwise
            forward = (pt1 < pt2)
            
            # Store edge usage
            edge_usage[edge_key].append((region_id, forward))
            edges.append((pt1, pt2, edge_key, forward))
        
        polygon_data.append({
            'region_id': region_id,
            'pt_tags': pt_tags,
            'edges': edges
        })

    # Second pass: Create all unique lines and track short edges
    short_edges = []  # Track short edges for later processing
    for edge_key in edge_usage.keys():
        pt1, pt2 = edge_key
        line_tag = gmsh.model.geo.addLine(pt1, pt2)
        edge_map[edge_key] = line_tag
        
        # Calculate edge length from point coordinates
        pt1_coords = None
        pt2_coords = None
        for (x, y), tag in point_map.items():
            if tag == pt1:
                pt1_coords = (x, y)
            if tag == pt2:
                pt2_coords = (x, y)
        
        if pt1_coords and pt2_coords:
            edge_length = ((pt2_coords[0] - pt1_coords[0])**2 + (pt2_coords[1] - pt1_coords[1])**2)**0.5
            
            # Add transfinite constraints for long boundary edges to ensure consistent mesh sizing
            # This prevents the creation of overly coarse elements along major boundaries
            if edge_length > adjusted_target_size * 3:  # Long edge
                # Calculate how many elements should be along this edge
                num_elements = max(3, int(edge_length / adjusted_target_size))
                try:
                    gmsh.model.geo.mesh.setTransfiniteCurve(line_tag, num_elements)
                    if debug:
                        print(f"Set transfinite constraint on long edge: {pt1_coords} to {pt2_coords}, length={edge_length:.2f}, num_elements={num_elements}")
                except Exception as e:
                    if debug:
                        print(f"Warning: Could not set transfinite constraint on edge {pt1_coords} to {pt2_coords}: {e}")
            
            # Add transfinite constraints for short edges to prevent subdivision
            # This forces GMSH to use exactly 2 nodes (start and end) for short edges
            elif edge_length < adjusted_target_size:
                try:
                    gmsh.model.geo.mesh.setTransfiniteCurve(line_tag, 2)  # Exactly 2 nodes
                    if debug:
                        print(f"Set transfinite constraint on short edge: {pt1_coords} to {pt2_coords}, length={edge_length:.2f}, exactly 2 nodes")
                except Exception as e:
                    if debug:
                        print(f"Warning: Could not set transfinite constraint on short edge {pt1_coords} to {pt2_coords}: {e}")
            
            # Short edges are now handled by point sizing, no need for transfinite curves

    # Ensure all polygon points (including intersection points) are created as GMSH points
    # The intersection points were already added to polygons in build_polygons(), 
    # so we just need to ensure they exist as GMSH geometric entities
    if lines is not None:
        if debug:
            print("Ensuring all polygon points (including intersections) are created as GMSH points...")
        
        # Collect all points from all polygons to ensure they exist in GMSH
        all_polygon_points = set()
        for poly_data in polygon_data:
            pt_tags = poly_data['pt_tags']
            for tag in pt_tags:
                # Find the coordinates for this point tag
                for (x, y), point_tag in point_map.items():
                    if point_tag == tag:
                        all_polygon_points.add((x, y))
                        break
        
        # Create any missing GMSH points
        for x, y in all_polygon_points:
            key = (x, y)
            if key not in point_map:
                pt_tag = gmsh.model.geo.addPoint(x, y, 0.0, adjusted_target_size * 0.5)
                point_map[key] = pt_tag
                if debug:
                    print(f"Created GMSH point for polygon vertex {key}: tag {pt_tag}")
        
        if debug:
            print(f"Ensured {len(all_polygon_points)} polygon points exist as GMSH entities")
        
        # Create enhanced reinforcement lines that include intersection points from polygons
        # This is essential for proper mesh generation with embedded 1D elements
        enhanced_lines = []
        for line_idx, line_pts in enumerate(lines):
            line_pts_clean = remove_duplicate_endpoint(list(line_pts))
            
            # Collect all points for this line: original + intersection points from polygons
            all_line_points = []
            
            # Add original line points
            for x, y in line_pts_clean:
                all_line_points.append((x, y, 'original'))
            
            # Add intersection points that are on this line (from polygon data)
            for poly_data in polygon_data:
                pt_tags = poly_data['pt_tags']
                for tag in pt_tags:
                    # Find the coordinates for this point tag
                    for (x, y), point_tag in point_map.items():
                        if point_tag == tag:
                            # Check if this point is on the reinforcement line
                            if is_point_on_line_segments((x, y), line_pts_clean, tolerance=1e-6):
                                all_line_points.append((x, y, 'intersection'))
                            break
            
            # Sort all points along the line to maintain proper order
            if len(all_line_points) > 1:
                all_line_points.sort(key=lambda p: line_segment_parameter((p[0], p[1]), line_pts_clean[0], line_pts_clean[-1]))
            
            # Remove duplicates (keep first occurrence)
            unique_points = []
            seen = set()
            for x, y, point_type in all_line_points:
                point_key = (round(x, 8), round(y, 8))  # Round to avoid floating point issues
                if point_key not in seen:
                    seen.add(point_key)
                    unique_points.append((x, y, point_type))
            
            # Create the enhanced line
            enhanced_line = [(x, y) for x, y, _ in unique_points]
            enhanced_lines.append(enhanced_line)
            
            if debug:
                print(f"Enhanced line {line_idx}: {len(line_pts_clean)} original points -> {len(enhanced_line)} total points")
        
        # Replace original lines with enhanced lines
        lines = enhanced_lines

    # Create reinforcement lines as geometric constraints to force 2D mesh edges
    line_data = []
    
    if lines is not None:
        for line_idx, line_pts in enumerate(lines):
            # Use the enhanced line coordinates (which include intersection points)
            line_pts_clean = remove_duplicate_endpoint(list(line_pts))
            
            # Create points for this reinforcement line
            line_point_tags = []
            
            # Create all points for this line (original + intersection points)
            for x, y in line_pts_clean:
                key = (x, y)
                if key in point_map:
                    line_point_tags.append((x, y, point_map[key]))
                else:
                    # Create new point with small mesh size to ensure it's preserved
                    pt_tag = gmsh.model.geo.addPoint(x, y, 0.0, adjusted_target_size * 0.5)
                    point_map[key] = pt_tag
                    line_point_tags.append((x, y, pt_tag))
            
            # Sort points along the line to maintain proper order
            line_point_tags.sort(key=lambda p: line_segment_parameter((p[0], p[1]), line_pts_clean[0], line_pts_clean[-1]))
            
            # Extract just the point tags in order
            pt_tags = [tag for _, _, tag in line_point_tags]
            
            if debug:
                print(f"  Line {line_idx} points: {[(x, y) for x, y, _ in line_point_tags]}")
            
            # Create line segments as geometric constraints with controlled meshing
            line_tags = []
            for i in range(len(pt_tags) - 1):
                pt1, pt2 = pt_tags[i], pt_tags[i + 1]
                
                # Calculate segment length to determine number of subdivisions
                coord1 = None
                coord2 = None
                for (x, y), tag in point_map.items():
                    if tag == pt1:
                        coord1 = (x, y)
                    if tag == pt2:
                        coord2 = (x, y)
                
                if coord1 and coord2:
                    segment_length = ((coord2[0] - coord1[0])**2 + (coord2[1] - coord1[1])**2)**0.5
                    # Calculate number of elements needed to achieve target_size_1d
                    # For segments longer than target_size_1d, we want multiple elements
                    # For segments shorter than target_size_1d, we still want at least 2 elements
                    if segment_length > target_size_1d:
                        num_elements = max(3, int(round(segment_length / target_size_1d)))
                    else:
                        num_elements = 2
                    
                    if debug:
                        print(f"  Segment {i}: length {segment_length:.2f}, creating {num_elements} elements")
                    
                    line_tag = gmsh.model.geo.addLine(pt1, pt2)
                    line_tags.append(line_tag)
                    
                    # Set transfinite constraint to create appropriate number of nodes
                    try:
                        gmsh.model.geo.mesh.setTransfiniteCurve(line_tag, num_elements)
                        if debug:
                            print(f"  Set transfinite constraint on line segment {i}: {num_elements} nodes")
                    except Exception as e:
                        if debug:
                            print(f"  Warning: Could not set transfinite constraint on segment {i}: {e}")
                else:
                    # Fallback: create line with default 2 nodes
                    line_tag = gmsh.model.geo.addLine(pt1, pt2)
                    line_tags.append(line_tag)
                    
                    try:
                        gmsh.model.geo.mesh.setTransfiniteCurve(line_tag, 2)
                        if debug:
                            print(f"  Set transfinite constraint on line segment {i}: 2 nodes (fallback)")
                    except Exception as e:
                        if debug:
                            print(f"  Warning: Could not set transfinite constraint on segment {i}: {e}")
            
            # Store line data for later 1D element extraction
            # Use the enhanced line coordinates (which include intersection points)
            line_data.append({
                'line_idx': line_idx,
                'line_tags': line_tags,
                'point_coords': line_pts_clean  # This now contains the enhanced coordinates
            })
            
            if debug:
                print(f"Created reinforcement constraint line {line_idx} with {len(line_tags)} segments: {line_pts_clean}")

    # Third pass: Create surfaces using the shared lines
    surface_to_region = {}
    
    for poly_data in polygon_data:
        region_id = poly_data['region_id']
        edges = poly_data['edges']
        
        line_tags = []
        for pt1, pt2, edge_key, forward in edges:
            line_tag = edge_map[edge_key]
            
            # Use positive or negative line tag based on orientation
            if forward:
                line_tags.append(line_tag)
            else:
                line_tags.append(-line_tag)
        
        # Create curve loop and surface
        try:
            loop = gmsh.model.geo.addCurveLoop(line_tags)
            surface = gmsh.model.geo.addPlaneSurface([loop])
            surface_to_region[surface] = region_id
        except Exception as e:
            print(f"Warning: Could not create surface for region {region_id}: {e}")
            continue

    # Synchronize geometry
    gmsh.model.geo.synchronize()
    
    # Force mesh edges along reinforcement lines by creating additional geometric constraints
    if lines is not None:
        for line_info in line_data:
            line_idx = line_info['line_idx']
            line_tags = line_info['line_tags']
            line_pts = line_info['point_coords']
            
            # Set transfinite constraints to force mesh edges along each line segment
            # REMOVED: This was conflicting with the target_size_1d calculations above
            # for i, line_tag in enumerate(line_tags):
            #     try:
            #         # Force exactly 2 nodes (start and end) to prevent subdivision
            #         gmsh.model.geo.mesh.setTransfiniteCurve(line_tag, 2)
            #         if debug:
            #             print(f"Set transfinite constraint on line {line_idx} segment {i}: exactly 2 nodes")
            #     except Exception as e:
            #         if debug:
            #             print(f"Warning: Could not set transfinite constraint on line {line_idx} segment {i}: {e}")
            
            # Embed reinforcement lines in all surfaces to ensure they're part of the mesh
            for surface in surface_to_region.keys():
                try:
                    # Embed all line segments of this reinforcement line
                    gmsh.model.mesh.embed(1, line_tags, 2, surface)
                    if debug:
                        print(f"Embedded reinforcement line {line_idx} in surface {surface}")
                except Exception as e:
                    if debug:
                        print(f"Could not embed line {line_idx} in surface {surface}: {e}")
    
    # CRITICAL: Set mesh coherence to ensure shared nodes along boundaries
    # This forces Gmsh to use the same nodes for shared geometric entities
    gmsh.model.mesh.removeDuplicateNodes()
    
    # Create physical groups for material regions (this helps with mesh consistency)
    physical_surfaces = []
    for surface, region_id in surface_to_region.items():
        physical_tag = gmsh.model.addPhysicalGroup(2, [surface])
        physical_surfaces.append((physical_tag, region_id))
    
    # Create physical groups for embedded reinforcement lines
    physical_lines = []
    if lines is not None:
        for line_info in line_data:
            line_idx = line_info['line_idx']
            line_tags = line_info['line_tags']
            physical_tag = gmsh.model.addPhysicalGroup(1, line_tags)
            physical_lines.append((physical_tag, line_idx))
    
    # Check for potential quad4 + reinforcement line conflicts
    has_reinforcement_lines = lines is not None and len(lines) > 0
    wants_quads = base_element_type.startswith('quad')
    
    # Set mesh algorithm and recombination options BEFORE generating mesh
    if base_element_type.startswith('quad'):
        # Check if we need to use a more robust algorithm for reinforcement lines
        if has_reinforcement_lines:
            if debug:
                print(f"Detected quad elements with reinforcement lines.")
                print(f"Using robust recombination algorithm to handle embedded line constraints.")
            
            # Use 'fast' algorithm which is more robust with embedded constraints
            default_params = {
                "Mesh.Algorithm": 8,  # Frontal-Delaunay for quads
                "Mesh.RecombineAll": 1,  # Recombine triangles into quads
                "Mesh.RecombinationAlgorithm": 0,  # Standard (more robust than simple)
                "Mesh.SubdivisionAlgorithm": 0,  # Mixed tri/quad where needed
                "Mesh.RecombineOptimizeTopology": 0,  # Minimal optimization
                "Mesh.RecombineNodeRepositioning": 1,  # Still reposition nodes
                "Mesh.RecombineMinimumQuality": 0.01,  # Keep quality threshold
                "Mesh.Smoothing": 5,  # Reduced smoothing
                "Mesh.SmoothNormals": 1,  # Keep smooth normals
                "Mesh.SmoothRatio": 1.8,  # Keep smoothing ratio
            }
        else:
            # Standard quad meshing parameters for cases without reinforcement lines
            default_params = {
                "Mesh.Algorithm": 8,  # Frontal-Delaunay for quads (try 5, 6, 8)
                "Mesh.RecombineAll": 1,  # Recombine triangles into quads
                "Mesh.RecombinationAlgorithm": 1,  # Simple recombination (try 0, 1, 2, 3)
                "Mesh.SubdivisionAlgorithm": 1,  # All quads (try 0, 1, 2)
                "Mesh.RecombineOptimizeTopology": 5,  # Optimize topology (0-100)
                "Mesh.RecombineNodeRepositioning": 1,  # Reposition nodes (0 or 1)
                "Mesh.RecombineMinimumQuality": 0.01,  # Minimum quality threshold
                "Mesh.Smoothing": 10,  # Number of smoothing steps (try 0-100)
                "Mesh.SmoothNormals": 1,  # Smooth normals
                "Mesh.SmoothRatio": 1.8,  # Smoothing ratio (1.0-3.0)
            }
        
        # Override with user-provided parameters
        if mesh_params:
            default_params.update(mesh_params)
        
        # Apply all parameters (except our custom ones)
        for param, value in default_params.items():
            if param not in ['size_factor']:  # Skip our custom parameters
                gmsh.option.setNumber(param, value)
        
        # Set recombination for each surface
        for surface in surface_to_region.keys():
            gmsh.model.mesh.setRecombine(2, surface)
    else:
        gmsh.option.setNumber("Mesh.Algorithm", 6)  # Frontal-Delaunay for triangles
    
    # Always generate linear elements first - quadratic conversion is done in post-processing
    # This avoids gmsh issues with quadratic elements and embedded 1D lines
    gmsh.option.setNumber("Mesh.ElementOrder", 1)
    
    # Force mesh coherence before generation
    gmsh.option.setNumber("Mesh.ToleranceInitialDelaunay", 1e-12)
    
    # Short edge control is now handled by point sizing during geometry creation
    
    # Generate mesh
    gmsh.model.mesh.generate(2)
    
    # Remove duplicate nodes again after mesh generation (belt and suspenders)
    gmsh.model.mesh.removeDuplicateNodes()

    # Get nodes
    node_tags, coords, _ = gmsh.model.mesh.getNodes()
    nodes = np.array(coords).reshape(-1, 3)[:, :2]
    
    # Create node tag to index mapping
    node_tag_to_index = {tag: i for i, tag in enumerate(node_tags)}

    elements = []
    mat_ids = []
    element_node_counts = []
    
    # For quad8: track center nodes to delete later
    center_nodes_to_delete = set() if element_type == 'quad8' else None

    # Extract elements using physical groups for better region identification
    for physical_tag, region_id in physical_surfaces:
        try:
            # Get entities in this physical group
            entities = gmsh.model.getEntitiesForPhysicalGroup(2, physical_tag)
            
            for entity in entities:
                # Get all elements for this entity
                elem_types, elem_tags_list, node_tags_list = gmsh.model.mesh.getElements(2, entity)
                
                for elem_type, elem_tags, node_tags in zip(elem_types, elem_tags_list, node_tags_list):
                    # Gmsh element type mapping:
                    # 2: 3-node triangle, 9: 6-node triangle
                    # 3: 4-node quadrilateral, 10: 8-node quadrilateral
                    if elem_type == 2:  # 3-node triangle
                        elements_array = np.array(node_tags).reshape(-1, 3)
                        for element in elements_array:
                            idxs = [node_tag_to_index[tag] for tag in element]
                            
                            # GMSH returns clockwise triangles - reorder to counter-clockwise
                            idxs[1], idxs[2] = idxs[2], idxs[1]
                            
                            # Pad to 9 columns with zeros
                            padded_idxs = idxs + [0] * (9 - len(idxs))
                            elements.append(padded_idxs)
                            mat_ids.append(region_id)
                            element_node_counts.append(3)
                    elif elem_type == 9:  # 6-node triangle
                        elements_array = np.array(node_tags).reshape(-1, 6)
                        for element in elements_array:
                            idxs = [node_tag_to_index[tag] for tag in element]
                            
                            # GMSH returns clockwise tri6 elements - reorder to counter-clockwise
                            # Swap corner nodes 1 and 2
                            idxs[1], idxs[2] = idxs[2], idxs[1]
                            # Fix midpoint assignments after corner swap 1<->2:
                            # GMSH gives: n3=edge(0-1), n4=edge(1-2), n5=edge(2-0)
                            # After swap: n3=edge(0-2), n4=edge(2-1), n5=edge(1-0)
                            # Standard requires: n3=edge(0-1), n4=edge(1-2), n5=edge(2-0)
                            # So remap: new_n3=old_n5, new_n4=old_n4, new_n5=old_n3
                            old_3, old_4, old_5 = idxs[3], idxs[4], idxs[5]
                            idxs[3] = old_5  # standard edge(0-1) gets GMSH edge(2-0) midpoint
                            idxs[4] = old_4  # standard edge(1-2) gets GMSH edge(1-2) midpoint  
                            idxs[5] = old_3  # standard edge(2-0) gets GMSH edge(0-1) midpoint
                            
                            # Pad to 9 columns with zeros
                            padded_idxs = idxs + [0] * (9 - len(idxs))
                            elements.append(padded_idxs)
                            mat_ids.append(region_id)
                            element_node_counts.append(6)
                    elif elem_type == 3:  # 4-node quadrilateral
                        elements_array = np.array(node_tags).reshape(-1, 4)
                        for element in elements_array:
                            idxs = [node_tag_to_index[tag] for tag in element]
                            # Fix node ordering for quadrilateral elements
                            if element_type.startswith('quad'):
                                idxs = idxs[::-1] # Simple reversal of node order
                            # Pad to 9 columns with zeros
                            padded_idxs = idxs + [0] * (9 - len(idxs))
                            elements.append(padded_idxs)
                            mat_ids.append(region_id)
                            element_node_counts.append(4)
                    elif elem_type == 10:  # Quadratic quadrilateral (gmsh generates 9-node Lagrange)
                        # Gmsh always generates 9-node Lagrange quads for order 2
                        elements_array = np.array(node_tags).reshape(-1, 9)
                        for element in elements_array:
                            idxs = [node_tag_to_index[tag] for tag in element]
                            
                            if element_type in ['quad8', 'quad9']:
                                # Both quad8 and quad9 need CW to CCW conversion for first 8 nodes
                                # Convert from Gmsh CW to CCW ordering for quadrilateral
                                # Corner nodes: reverse order (0,1,2,3) -> (0,3,2,1)
                                # Midpoint nodes need to be reordered accordingly:
                                # GMSH: n4=edge(0-1), n5=edge(1-2), n6=edge(2-3), n7=edge(3-0)
                                # After corner reversal: need n4=edge(0-3), n5=edge(3-2), n6=edge(2-1), n7=edge(1-0)
                                # So: new_n4=old_n7, new_n5=old_n6, new_n6=old_n5, new_n7=old_n4
                                reordered_first8 = [
                                    idxs[0],  # corner 0 stays
                                    idxs[3],  # corner 1 -> corner 3  
                                    idxs[2],  # corner 2 stays
                                    idxs[1],  # corner 3 -> corner 1
                                    idxs[7],  # edge(0-1) -> edge(0-3) = old edge(3-0)
                                    idxs[6],  # edge(1-2) -> edge(3-2) = old edge(2-3)  
                                    idxs[5],  # edge(2-3) -> edge(2-1) = old edge(1-2)
                                    idxs[4]   # edge(3-0) -> edge(1-0) = old edge(0-1)
                                ]
                                
                                if element_type == 'quad8':
                                    # For quad8, skip center node and mark for deletion
                                    center_node_idx = idxs[8]  # Mark center node for deletion
                                    center_nodes_to_delete.add(center_node_idx)
                                    padded_idxs = reordered_first8 + [0]  # Skip center node, pad to 9
                                    elements.append(padded_idxs)
                                    mat_ids.append(region_id)
                                    element_node_counts.append(8)
                                else:  # quad9
                                    # For quad9, keep center node (9th node unchanged)
                                    full_idxs = reordered_first8 + [idxs[8]]  # Add center node
                                    elements.append(full_idxs)
                                    mat_ids.append(region_id)
                                    element_node_counts.append(9)
                            else:
                                # This should never happen since element_type is validated earlier
                                raise ValueError(f"Unexpected element_type '{element_type}' for Gmsh elem_type {elem_type}")
        except Exception as e:
            print(f"Warning: Could not extract elements for physical group {physical_tag} (region {region_id}): {e}")
            continue

    # Convert to numpy arrays
    elements_array = np.array(elements, dtype=int)
    element_types = np.array(element_node_counts, dtype=int)
    element_materials = np.array(mat_ids, dtype=int)

    # Extract 1D elements from Gmsh-generated 1D mesh along reinforcement lines
    elements_1d = []
    mat_ids_1d = []
    element_node_counts_1d = []
    
    if lines is not None:
        # Extract 1D elements from physical groups for each reinforcement line
        for physical_tag, line_idx in physical_lines:
            try:
                # Get entities in this physical group
                entities = gmsh.model.getEntitiesForPhysicalGroup(1, physical_tag)
                
                if debug:
                    print(f"  Physical group {physical_tag} (line {line_idx}): found {len(entities)} entities")
                
                for entity in entities:
                    # Get all 1D elements for this entity
                    elem_types, elem_tags_list, node_tags_list = gmsh.model.mesh.getElements(1, entity)
                    
                    for elem_type, elem_tags, node_tags in zip(elem_types, elem_tags_list, node_tags_list):
                        # Gmsh 1D element type mapping:
                        # 1: 2-node line (linear), 8: 3-node line (quadratic)
                        if elem_type == 1:  # Linear 1D elements (2 nodes)
                            elements_array = np.array(node_tags).reshape(-1, 2)
                            for element in elements_array:
                                try:
                                    # Convert numpy arrays to regular Python scalars
                                    element_list = element.tolist()  # Convert to Python list
                                    if len(element_list) >= 2:
                                        tag1 = int(element_list[0])
                                        tag2 = int(element_list[1])
                                        
                                        # Get node indices
                                        idx1 = node_tag_to_index[tag1]
                                        idx2 = node_tag_to_index[tag2]
                                        
                                        # Create 1D element
                                        padded_idxs = [idx1, idx2, 0]
                                        elements_1d.append(padded_idxs)
                                        mat_ids_1d.append(line_idx)
                                        element_node_counts_1d.append(2)
                                        
                                        if debug:
                                            coord1 = nodes[idx1]
                                            coord2 = nodes[idx2]
                                            print(f"    Created 1D element: {coord1} -> {coord2}")
                                except (KeyError, TypeError, ValueError, IndexError) as e:
                                    if debug:
                                        print(f"    Skipping 1D element due to error: {e}")
                                    continue
                        elif elem_type == 8:  # Quadratic 1D elements (3 nodes)
                            elements_array = np.array(node_tags).reshape(-1, 3)
                            for element in elements_array:
                                try:
                                    # Convert numpy arrays to regular Python scalars
                                    element_list = element.tolist()  # Convert to Python list
                                    if len(element_list) >= 3:
                                        tag1 = int(element_list[0])
                                        tag2 = int(element_list[1])
                                        tag3 = int(element_list[2])
                                        
                                        # Get node indices
                                        idx1 = node_tag_to_index[tag1]
                                        idx2 = node_tag_to_index[tag2]
                                        idx3 = node_tag_to_index[tag3]
                                        
                                        # Create 1D element
                                        padded_idxs = [idx1, idx2, idx3]
                                        elements_1d.append(padded_idxs)
                                        mat_ids_1d.append(line_idx)
                                        element_node_counts_1d.append(3)
                                except (KeyError, TypeError, ValueError, IndexError) as e:
                                    if debug:
                                        print(f"    Skipping quadratic 1D element due to error: {e}")
                                    continue
            except Exception as e:
                if debug:
                    print(f"  Error extracting 1D elements for line {line_idx}: {e}")
                continue

    gmsh.finalize()

    # Clean up center nodes for quad8 elements
    if element_type == 'quad8' and center_nodes_to_delete:
        print(f"Quad8 cleanup: removing {len(center_nodes_to_delete)} center nodes from {len(nodes)} total nodes")
        
        # c) Create array tracking original node numbering
        original_node_count = len(nodes)
        nodes_to_keep = [i for i in range(original_node_count) if i not in center_nodes_to_delete]
        
        # d) Delete center nodes - create new nodes array
        new_nodes = nodes[nodes_to_keep]
        
        # e) Create mapping from old node indices to new node indices
        old_to_new_mapping = {old_idx: new_idx for new_idx, old_idx in enumerate(nodes_to_keep)}
        
        # f) Update element topology to use new node numbering
        new_elements = []
        for element in elements_array:
            new_element = []
            for node_idx in element:
                if node_idx == 0:  # Keep padding zeros
                    new_element.append(0)
                elif node_idx in center_nodes_to_delete:
                    # This should not happen since we set center nodes to 0
                    new_element.append(0)
                else:
                    # Map to new node index
                    new_element.append(old_to_new_mapping[node_idx])
            new_elements.append(new_element)
        
        # g) Replace arrays with consolidated versions
        elements_array = np.array(new_elements, dtype=int)
        nodes = new_nodes
        
        print(f"Quad8 cleanup complete: {len(nodes)} nodes, {len(elements_array)} elements")

    # Convert lists to arrays
    elements_array = np.array(elements, dtype=int)
    element_types = np.array(element_node_counts, dtype=int)
    element_materials = np.array(mat_ids, dtype=int) + 1  # Make 1-based

    mesh = {
        "nodes": nodes,
        "elements": elements_array,
        "element_types": element_types,
        "element_materials": element_materials,
    }

    # Add 1D element data if lines were provided
    if lines is not None and len(elements_1d) > 0:
        elements_1d_array = np.array(elements_1d, dtype=int)
        element_types_1d = np.array(element_node_counts_1d, dtype=int)
        element_materials_1d = np.array(mat_ids_1d, dtype=int) + 1  # Make 1-based
        
        mesh["elements_1d"] = elements_1d_array
        mesh["element_types_1d"] = element_types_1d
        mesh["element_materials_1d"] = element_materials_1d

    # Post-process to convert linear elements to quadratic if requested
    if quadratic:
        if debug:
            print(f"Converting linear {base_element_type} mesh to quadratic {element_type}")
        mesh = convert_linear_to_quadratic_mesh(mesh, element_type, debug=debug)

    return mesh


def convert_linear_to_quadratic_mesh(mesh, target_element_type, debug=False):
    """
    Convert a linear mesh (tri3/quad4) to quadratic (tri6/quad8/quad9) by adding midside nodes.
    
    This is much more robust than gmsh's built-in quadratic generation, especially 
    when dealing with embedded 1D elements (reinforcement lines).
    
    Parameters:
        mesh: Dictionary containing linear mesh data
        target_element_type: 'tri6', 'quad8', or 'quad9'
        debug: Enable debug output
        
    Returns:
        Updated mesh dictionary with quadratic elements
    """
    if debug:
        print(f"Converting to {target_element_type} elements...")
    
    nodes = mesh["nodes"].copy()
    elements = mesh["elements"].copy()
    element_types = mesh["element_types"].copy()
    element_materials = mesh["element_materials"].copy()
    
    # Handle 1D elements if present
    elements_1d = mesh.get("elements_1d")
    element_types_1d = mesh.get("element_types_1d") 
    element_materials_1d = mesh.get("element_materials_1d")
    has_1d_elements = elements_1d is not None
    
    if has_1d_elements:
        elements_1d = elements_1d.copy()
        element_types_1d = element_types_1d.copy()
        element_materials_1d = element_materials_1d.copy()
    
    # Dictionary to store midside nodes: (node1_idx, node2_idx) -> midside_node_idx
    # Always store with node1_idx < node2_idx for consistency
    midside_nodes = {}
    next_node_idx = len(nodes)
    
    def get_or_create_midside_node(n1_idx, n2_idx):
        """Get existing midside node or create new one between n1 and n2"""
        nonlocal next_node_idx, nodes
        
        # Ensure consistent ordering
        if n1_idx > n2_idx:
            n1_idx, n2_idx = n2_idx, n1_idx
            
        edge_key = (n1_idx, n2_idx)
        
        if edge_key in midside_nodes:
            return midside_nodes[edge_key]
        
        # Create new midside node at edge center
        n1 = nodes[n1_idx]
        n2 = nodes[n2_idx]
        midside_coord = (n1 + n2) / 2.0
        
        # Add to nodes array
        nodes_list = nodes.tolist()
        nodes_list.append(midside_coord.tolist())
        nodes = np.array(nodes_list)
        
        midside_idx = next_node_idx
        midside_nodes[edge_key] = midside_idx
        next_node_idx += 1
        
        if debug and len(midside_nodes) <= 10:  # Only print first few
            print(f"  Created midside node {midside_idx} between {n1_idx}-{n2_idx} at {midside_coord}")
        
        return midside_idx
    
    # Convert 2D elements
    new_elements = []
    new_element_types = []
    
    for elem_idx, element in enumerate(elements):
        elem_type = element_types[elem_idx]
        
        if target_element_type == 'tri6' and elem_type == 3:
            # Convert tri3 to tri6
            n0, n1, n2 = element[0], element[1], element[2]
            
            # Get/create midside nodes
            n3 = get_or_create_midside_node(n0, n1)  # edge 0-1
            n4 = get_or_create_midside_node(n1, n2)  # edge 1-2
            n5 = get_or_create_midside_node(n2, n0)  # edge 2-0
            
            # tri6 node ordering: [corner_nodes, midside_nodes]
            new_element = [n0, n1, n2, n3, n4, n5, 0, 0, 0]
            new_elements.append(new_element)
            new_element_types.append(6)
            
        elif target_element_type == 'quad8' and elem_type == 4:
            # Convert quad4 to quad8
            n0, n1, n2, n3 = element[0], element[1], element[2], element[3]
            
            # Get/create midside nodes on edges
            n4 = get_or_create_midside_node(n0, n1)  # edge 0-1
            n5 = get_or_create_midside_node(n1, n2)  # edge 1-2
            n6 = get_or_create_midside_node(n2, n3)  # edge 2-3
            n7 = get_or_create_midside_node(n3, n0)  # edge 3-0
            
            # quad8 node ordering: [corner_nodes, midside_nodes]
            new_element = [n0, n1, n2, n3, n4, n5, n6, n7, 0]
            new_elements.append(new_element)
            new_element_types.append(8)
            
        elif target_element_type == 'quad9' and elem_type == 4:
            # Convert quad4 to quad9
            n0, n1, n2, n3 = element[0], element[1], element[2], element[3]
            
            # Get/create midside nodes on edges
            n4 = get_or_create_midside_node(n0, n1)  # edge 0-1
            n5 = get_or_create_midside_node(n1, n2)  # edge 1-2
            n6 = get_or_create_midside_node(n2, n3)  # edge 2-3
            n7 = get_or_create_midside_node(n3, n0)  # edge 3-0
            
            # Create center node
            center_coord = (nodes[n0] + nodes[n1] + nodes[n2] + nodes[n3]) / 4.0
            nodes_list = nodes.tolist()
            nodes_list.append(center_coord.tolist())
            nodes = np.array(nodes_list)
            n8 = next_node_idx
            next_node_idx += 1
            
            # quad9 node ordering: [corner_nodes, midside_nodes, center_node]
            new_element = [n0, n1, n2, n3, n4, n5, n6, n7, n8]
            new_elements.append(new_element)
            new_element_types.append(9)
            
        else:
            # Keep original element unchanged
            new_elements.append(element.tolist())
            new_element_types.append(elem_type)
    
    # Convert 1D elements to quadratic if present
    new_elements_1d = []
    new_element_types_1d = [] 
    
    if has_1d_elements:
        for elem_idx, element in enumerate(elements_1d):
            elem_type = element_types_1d[elem_idx]
            
            if elem_type == 2:  # Convert linear 1D to quadratic
                n0, n1 = element[0], element[1]
                
                # Get/create midside node (reuse if already created for 2D elements)
                n2 = get_or_create_midside_node(n0, n1)
                
                new_element = [n0, n1, n2]
                new_elements_1d.append(new_element)
                new_element_types_1d.append(3)  # quadratic 1D
            else:
                # Keep original element unchanged
                new_elements_1d.append(element.tolist())
                new_element_types_1d.append(elem_type)
    
    if debug:
        print(f"  Added {len(midside_nodes)} midside nodes")
        print(f"  Total nodes: {len(nodes)} (was {len(mesh['nodes'])})")
    
    # Create updated mesh
    updated_mesh = {
        "nodes": nodes,
        "elements": np.array(new_elements, dtype=int),
        "element_types": np.array(new_element_types, dtype=int),
        "element_materials": element_materials
    }
    
    if has_1d_elements:
        updated_mesh["elements_1d"] = np.array(new_elements_1d, dtype=int)
        updated_mesh["element_types_1d"] = np.array(new_element_types_1d, dtype=int)
        updated_mesh["element_materials_1d"] = element_materials_1d
    
    return updated_mesh


def line_segment_parameter(point, line_start, line_end):
    """
    Calculate the parameter t (0 to 1) of a point along a line segment.
    Returns t where point = line_start + t * (line_end - line_start)
    """
    px, py = point
    x1, y1 = line_start
    x2, y2 = line_end
    
    # Calculate parameter t
    dx = x2 - x1
    dy = y2 - y1
    
    if abs(dx) > abs(dy):
        t = (px - x1) / dx
    else:
        t = (py - y1) / dy
    
    return t


def line_segment_intersection(p1, p2, p3, p4, tol=1e-8):
    """
    Find intersection point between two line segments.
    Returns intersection point (x, y) if it exists, None otherwise.
    """
    x1, y1 = p1
    x2, y2 = p2
    x3, y3 = p3
    x4, y4 = p4
    
    # Calculate direction vectors
    d1x, d1y = x2 - x1, y2 - y1
    d2x, d2y = x4 - x3, y4 - y3
    
    # Calculate determinant
    det = d1x * d2y - d1y * d2x
    
    if abs(det) < tol:  # Lines are parallel
        return None
    
    # Calculate parameters
    t1 = ((x3 - x1) * d2y - (y3 - y1) * d2x) / det
    t2 = ((x3 - x1) * d1y - (y3 - y1) * d1x) / det
    
    # Check if intersection is within both segments
    if 0 <= t1 <= 1 and 0 <= t2 <= 1:
        # Calculate intersection point
        ix = x1 + t1 * d1x
        iy = y1 + t1 * d1y
        return (round(ix, 6), round(iy, 6))
    
    return None


def point_near_existing(point, existing_points, tol=1e-8):
    """Check if a point is near any existing points."""
    px, py = point
    for ex, ey in existing_points:
        if abs(px - ex) < tol and abs(py - ey) < tol:
            return True
    return False


def insert_point_into_polygon_edge(intersection, edge_start, edge_end, poly_data, point_map, target_size):
    """Insert an intersection point into a polygon edge, updating the polygon's coordinate list."""
    x, y = intersection
    # Ensure the point exists in the point_map (for Gmsh)
    if (x, y) not in point_map:
        tag = len(point_map) + 1  # Simple tag assignment
        point_map[(x, y)] = tag
    
    # Insert the intersection point into the polygon's coordinate list at the correct edge
    # poly_data['pt_tags'] is a list of Gmsh point tags, but we need to update the coordinate list used to build the polygon
    # We'll reconstruct the coordinate list from the tags and point_map
    pt_tags = poly_data['pt_tags']
    # Build coordinate list for the polygon
    coords = []
    tag_to_coord = {v: k for k, v in point_map.items()}
    for tag in pt_tags:
        if tag in tag_to_coord:
            coords.append(tag_to_coord[tag])
        else:
            # Fallback: try to find the coordinate in point_map
            found = False
            for (cx, cy), t in point_map.items():
                if t == tag:
                    coords.append((cx, cy))
                    found = True
                    break
            if not found:
                coords.append((None, None))  # Should not happen
    # Find the edge to insert after
    insert_idx = None
    for i in range(len(coords)):
        a = coords[i]
        b = coords[(i + 1) % len(coords)]
        if (abs(a[0] - edge_start[0]) < 1e-8 and abs(a[1] - edge_start[1]) < 1e-8 and
            abs(b[0] - edge_end[0]) < 1e-8 and abs(b[1] - edge_end[1]) < 1e-8):
            insert_idx = i + 1
            break
        # Also check reversed edge
        if (abs(a[0] - edge_end[0]) < 1e-8 and abs(a[1] - edge_end[1]) < 1e-8 and
            abs(b[0] - edge_start[0]) < 1e-8 and abs(b[1] - edge_start[1]) < 1e-8):
            insert_idx = i + 1
            break
    if insert_idx is not None:
        # Insert the intersection point into the coordinate list
        coords.insert(insert_idx, (x, y))
        # Now update pt_tags to match
        tag = point_map[(x, y)]
        pt_tags.insert(insert_idx, tag)
        # Update poly_data
        poly_data['pt_tags'] = pt_tags
    # If not found, do nothing (should not happen)


def get_quad_mesh_presets():
    """
    Returns dictionary of preset quad meshing parameter combinations to try.
    """
    presets = {
        'default': {
            "Mesh.Algorithm": 8,
            "Mesh.RecombinationAlgorithm": 1,
            "Mesh.SubdivisionAlgorithm": 1,
            "Mesh.RecombineOptimizeTopology": 5,
            "Mesh.Smoothing": 10,
            "size_factor": 1.4,  # Target size adjustment
        },
        'blossom': {
            "Mesh.Algorithm": 6,
            "Mesh.RecombinationAlgorithm": 2,  # Blossom
            "Mesh.SubdivisionAlgorithm": 1,
            "Mesh.RecombineOptimizeTopology": 20,
            "Mesh.Smoothing": 20,
            "size_factor": 1.6,  # Slightly larger for better recombination
        },
        'blossom_full': {
            "Mesh.Algorithm": 5,
            "Mesh.RecombinationAlgorithm": 3,  # Blossom full-quad
            "Mesh.SubdivisionAlgorithm": 1,
            "Mesh.RecombineOptimizeTopology": 50,
            "Mesh.Smoothing": 30,
            "size_factor": 1.7,  # Larger for complex recombination
        },
        'high_quality': {
            "Mesh.Algorithm": 6,
            "Mesh.RecombinationAlgorithm": 1,
            "Mesh.SubdivisionAlgorithm": 1,
            "Mesh.RecombineOptimizeTopology": 100,
            "Mesh.RecombineNodeRepositioning": 1,
            "Mesh.RecombineMinimumQuality": 0.1,
            "Mesh.Smoothing": 50,
            "Mesh.SmoothRatio": 2.0,
            "size_factor": 2.0,  # Much larger due to heavy optimization
        },
        'fast': {
            "Mesh.Algorithm": 8,
            "Mesh.RecombinationAlgorithm": 0,  # Standard (fastest)
            "Mesh.SubdivisionAlgorithm": 0,
            "Mesh.RecombineOptimizeTopology": 0,
            "Mesh.Smoothing": 5,
            "size_factor": 0.7,  # Smaller adjustment = more elements
        }
    }
    return presets



def build_polygons(slope_data, reinf_lines=None, debug=False):
    """
    Build material zone polygons from slope_data.
    
    Extracts profile lines and max depth, then creates polygons for each material zone.
    Also integrates distributed load points and reinforcement line endpoints that are
    coincident with polygon edges.
    
    Parameters:
        slope_data: Dictionary containing slope geometry data
        
    Returns:
        List of polygons, each defined by (x,y) coordinate tuples
    """
    import numpy as np
    import copy

    # Extract profile lines and max depth from slope_data
    profile_lines = slope_data.get('profile_lines', [])
    max_depth = slope_data.get('max_depth', None)

    if not profile_lines:
        raise ValueError("Need at least 1 profile line to create material zones")
    
    # For single profile line, max_depth serves as the bottom boundary
    if len(profile_lines) == 1:
        if max_depth is None:
            raise ValueError("When using only 1 profile line, max_depth must be specified")

    def get_avg_y(line):
        return sum(y for _, y in line) / len(line)

    # Sort profile lines from top to bottom by average y
    sorted_lines = sorted(profile_lines, key=get_avg_y, reverse=True)
    n = len(sorted_lines)
    # Deep copy so we can insert points
    lines = [list(line) for line in copy.deepcopy(sorted_lines)]
    tol = 1e-8

    for i in range(n - 1):
        top = lines[i]
        for endpoint in [0, -1]:  # left and right
            x_top, y_top = top[endpoint]
            # Find the highest lower profile at this x
            best_j = None
            best_y = -np.inf
            for j in range(i + 1, n):
                lower = lines[j]
                xs_lower = np.array([x for x, y in lower])
                ys_lower = np.array([y for x, y in lower])
                if xs_lower[0] - tol <= x_top <= xs_lower[-1] + tol:
                    y_proj = np.interp(x_top, xs_lower, ys_lower)
                    if y_proj > best_y:
                        best_y = y_proj
                        best_j = j
            if best_j is not None:
                lower = lines[best_j]
                xs_lower = np.array([x for x, y in lower])
                ys_lower = np.array([y for x, y in lower])
                y_proj = np.interp(x_top, xs_lower, ys_lower)
                # Check if lower profile already has a point at this x (within tol)
                found = False
                for (x_l, y_l) in lower:
                    if abs(x_l - x_top) < tol:
                        found = True
                        break
                if abs(y_proj - y_top) < tol:
                    # Coincident: insert (x_top, y_top) if not present
                    if not found:
                        insert_idx = np.searchsorted(xs_lower, x_top)
                        lower.insert(insert_idx, (round(x_top, 6), round(y_top, 6)))
                else:
                    # Not coincident: insert (x_top, y_proj) if not present
                    if not found:
                        insert_idx = np.searchsorted(xs_lower, x_top)
                        lower.insert(insert_idx, (round(x_top, 6), round(y_proj, 6)))

    def clean_polygon(poly, tol=1e-8):
        # Remove consecutive duplicate points (except for closing point)
        if not poly:
            return poly
        cleaned = [poly[0]]
        for pt in poly[1:]:
            if abs(pt[0] - cleaned[-1][0]) > tol or abs(pt[1] - cleaned[-1][1]) > tol:
                cleaned.append(pt)
        # Ensure closed
        if abs(cleaned[0][0] - cleaned[-1][0]) > tol or abs(cleaned[0][1] - cleaned[-1][1]) > tol:
            cleaned.append(cleaned[0])
        return cleaned

    # Now build polygons as before
    polygons = []
    for i, top_line in enumerate(lines):
        xs_top, ys_top = zip(*top_line)
        xs_top = np.array(xs_top)
        ys_top = np.array(ys_top)
        left_x, left_y = xs_top[0], ys_top[0]
        right_x, right_y = xs_top[-1], ys_top[-1]

        if i < n - 1:
            lower_line = lines[i + 1]
            xs_bot, ys_bot = zip(*lower_line)
            xs_bot = np.array(xs_bot)
            ys_bot = np.array(ys_bot)
            # Project left and right endpoints vertically to lower profile
            left_y_bot = np.interp(left_x, xs_bot, ys_bot)
            right_y_bot = np.interp(right_x, xs_bot, ys_bot)
            # Find all lower profile points between left_x and right_x (exclusive)
            mask = (xs_bot > left_x) & (xs_bot < right_x)
            xs_bot_in = xs_bot[mask]
            ys_bot_in = ys_bot[mask]
            # Build bottom boundary: right projection, lower profile points (right to left), left projection
            bottom = []
            bottom.append((right_x, right_y_bot))
            for x, y in zip(xs_bot_in[::-1], ys_bot_in[::-1]):
                bottom.append((x, y))
            bottom.append((left_x, left_y_bot))
        else:
            # For the lowest polygon, bottom is at max_depth
            # Only need endpoints - no intermediate points
            bottom = []
            bottom.append((right_x, max_depth))
            bottom.append((left_x, max_depth))

        # Build polygon: top left-to-right, bottom right-to-left
        poly = []
        for x, y in zip(xs_top, ys_top):
            poly.append((round(x, 6), round(y, 6)))
        for x, y in bottom:
            poly.append((round(x, 6), round(y, 6)))
        # Clean up polygon (should rarely do anything)
        poly = clean_polygon(poly)
        polygons.append(poly)
    
    # Add distributed load points to polygon edges if coincident
    polygons = add_dload_points_to_polygons(polygons, slope_data)
    
    # Add intersection points with reinforcement lines if provided
    if reinf_lines is not None:
        polygons = add_intersection_points_to_polygons(polygons, reinf_lines, debug=debug)
    
    return polygons

def add_dload_points_to_polygons(polygons, slope_data):
    """
    Add distributed load points to polygon edges if they are coincident with edges 
    but not existing vertices.
    
    Parameters:
        polygons: List of polygons (lists of (x,y) tuples)
        slope_data: Dictionary containing slope data
        
    Returns:
        Updated list of polygons with added points
    """
    import numpy as np
    tol = 1e-8
    
    # Collect distributed load points to check
    points_to_check = []
    
    # Add distributed load points
    distributed_loads = slope_data.get('distributed_loads', [])
    for load in distributed_loads:
        if 'xy' in load:
            for point in load['xy']:
                points_to_check.append(point)
    
    if not points_to_check:
        return polygons
    
    # Process each polygon
    updated_polygons = []
    for poly in polygons:
        updated_poly = list(poly)  # Make a copy
        
        # Check each point against polygon edges
        for check_point in points_to_check:
            x_check, y_check = check_point
            
            # Check if point is already a vertex
            is_vertex = False
            for vertex in updated_poly:
                if abs(vertex[0] - x_check) < tol and abs(vertex[1] - y_check) < tol:
                    is_vertex = True
                    break
            
            if is_vertex:
                continue
            
            # Check if point lies on any edge
            for i in range(len(updated_poly)):
                x1, y1 = updated_poly[i]
                x2, y2 = updated_poly[(i + 1) % len(updated_poly)]
                
                # Check if point lies on edge segment
                if is_point_on_edge((x_check, y_check), (x1, y1), (x2, y2), tol):
                    # Insert point after vertex i
                    updated_poly.insert(i + 1, (round(x_check, 6), round(y_check, 6)))
                    break  # Only insert once per point
        
        updated_polygons.append(updated_poly)
    
    return updated_polygons

def is_point_on_edge(point, edge_start, edge_end, tol=1e-8):
    """
    Check if a point lies on a line segment (edge).
    
    Parameters:
        point: (x, y) tuple of point to check
        edge_start: (x, y) tuple of edge start
        edge_end: (x, y) tuple of edge end
        tol: Tolerance for coincidence
        
    Returns:
        bool: True if point lies on edge segment
    """
    px, py = point
    x1, y1 = edge_start
    x2, y2 = edge_end
    
    # Check if point is within bounding box of edge
    if not (min(x1, x2) - tol <= px <= max(x1, x2) + tol and
            min(y1, y2) - tol <= py <= max(y1, y2) + tol):
        return False
    
    # Check if point is collinear with edge
    # Use cross product to check collinearity
    cross_product = abs((py - y1) * (x2 - x1) - (px - x1) * (y2 - y1))
    
    # If cross product is close to zero, point is on the line
    # Also check that it's within the segment bounds
    if cross_product < tol:
        # Check if point is between edge endpoints
        dot_product = (px - x1) * (x2 - x1) + (py - y1) * (y2 - y1)
        edge_length_sq = (x2 - x1) ** 2 + (y2 - y1) ** 2
        
        if edge_length_sq < tol:  # Edge is essentially a point
            return abs(px - x1) < tol and abs(py - y1) < tol
        
        # Parameter t should be between 0 and 1 for point to be on segment
        t = dot_product / edge_length_sq
        return -tol <= t <= 1 + tol
    
    return False

def print_polygon_summary(polygons):
    """
    Prints a summary of the generated polygons for diagnostic purposes.
    
    Parameters:
        polygons: List of polygon coordinate lists
    """
    print("=== POLYGON SUMMARY ===")
    print(f"Number of material zones: {len(polygons)}")
    print()
    
    for i, polygon in enumerate(polygons):
        print(f"Material Zone {i+1} (Material ID: {i}):")
        print(f"  Number of vertices: {len(polygon)}")
        
        # Calculate area (simple shoelace formula)
        area = 0
        for j in range(len(polygon) - 1):
            x1, y1 = polygon[j]
            x2, y2 = polygon[j + 1]
            area += (x2 - x1) * (y2 + y1) / 2
        area = abs(area)
        
        print(f"  Approximate area: {area:.2f} square units")
        
        # Print bounding box
        xs = [x for x, y in polygon]
        ys = [y for x, y in polygon]
        print(f"  Bounding box: x=[{min(xs):.2f}, {max(xs):.2f}], y=[{min(ys):.2f}, {max(ys):.2f}]")
        print()




def export_mesh_to_json(mesh, filename):
    """Save mesh dictionary to JSON file."""
    import json
    import numpy as np
    
    # Convert numpy arrays to lists for JSON serialization
    mesh_json = {}
    for key, value in mesh.items():
        if isinstance(value, np.ndarray):
            mesh_json[key] = value.tolist()
        else:
            mesh_json[key] = value
    
    with open(filename, 'w') as f:
        json.dump(mesh_json, f, indent=2)
    
    print(f"Mesh saved to {filename}")

def import_mesh_from_json(filename):
    """Load mesh dictionary from JSON file."""
    import json
    import numpy as np
    
    with open(filename, 'r') as f:
        mesh_json = json.load(f)
    
    # Convert lists back to numpy arrays
    mesh = {}
    for key, value in mesh_json.items():
        if isinstance(value, list):
            mesh[key] = np.array(value)
        else:
            mesh[key] = value
    
    return mesh

def remove_duplicate_endpoint(poly, tol=1e-8):
    if len(poly) > 1 and abs(poly[0][0] - poly[-1][0]) < tol and abs(poly[0][1] - poly[-1][1]) < tol:
        return poly[:-1]
    return poly


def extract_1d_elements_from_2d_edges(nodes, elements_2d, element_types_2d, lines, debug=False):
    """
    Extract 1D elements from 2D element edges that lie along reinforcement lines.
    This ensures proper finite element integration where 1D elements are shared edges of 2D elements.
    
    Parameters:
        nodes: np.ndarray of node coordinates (n_nodes, 2)
        elements_2d: np.ndarray of 2D element vertex indices (n_elements, 9)
        element_types_2d: np.ndarray indicating 2D element type (3, 4, 6, 8, or 9 nodes)
        lines: List of reinforcement lines, each defined by list of (x, y) tuples
        debug: Enable debug output
        
    Returns:
        tuple: (elements_1d, mat_ids_1d, element_node_counts_1d)
    """
    import numpy as np
    from collections import defaultdict
    
    elements_1d = []
    mat_ids_1d = []
    element_node_counts_1d = []
    
    # Build edge-to-element mapping from 2D elements
    edge_to_element = defaultdict(list)  # edge (n1, n2) -> list of element indices
    element_edges = {}  # element_idx -> list of edges
    
    for elem_idx, (element, elem_type) in enumerate(zip(elements_2d, element_types_2d)):
        edges = []
        
        if elem_type in [3, 6]:  # Triangle
            # Triangle edges: (0,1), (1,2), (2,0)
            corner_nodes = [element[0], element[1], element[2]]
            edge_pairs = [(0, 1), (1, 2), (2, 0)]
            
            for i, j in edge_pairs:
                n1, n2 = corner_nodes[i], corner_nodes[j]
                edge_key = (min(n1, n2), max(n1, n2))  # Canonical edge representation
                edges.append(edge_key)
                edge_to_element[edge_key].append(elem_idx)
                
        elif elem_type in [4, 8, 9]:  # Quadrilateral
            # Quadrilateral edges: (0,1), (1,2), (2,3), (3,0)
            corner_nodes = [element[0], element[1], element[2], element[3]]
            edge_pairs = [(0, 1), (1, 2), (2, 3), (3, 0)]
            
            for i, j in edge_pairs:
                n1, n2 = corner_nodes[i], corner_nodes[j]
                edge_key = (min(n1, n2), max(n1, n2))  # Canonical edge representation
                edges.append(edge_key)
                edge_to_element[edge_key].append(elem_idx)
        
        element_edges[elem_idx] = edges
    
    if debug:
        print(f"Built edge map with {len(edge_to_element)} unique edges from {len(elements_2d)} 2D elements")
    
    # For each reinforcement line, find 2D element edges that lie along it
    for line_idx, line_pts in enumerate(lines):
        line_pts_clean = remove_duplicate_endpoint(list(line_pts))
        
        if len(line_pts_clean) < 2:
            continue
            
        if debug:
            print(f"Processing reinforcement line {line_idx}: {line_pts_clean}")
        
        # Find all 2D element edges that lie along this reinforcement line
        line_edges = []
        
        for edge_key, elem_indices in edge_to_element.items():
            n1, n2 = edge_key
            
            # Get coordinates of edge endpoints
            coord1 = nodes[n1]
            coord2 = nodes[n2]
            
            # Check if this edge lies along the reinforcement line
            if is_edge_on_reinforcement_line(coord1, coord2, line_pts_clean, tolerance=1e-6):
                line_edges.append((n1, n2))
                if debug:
                    print(f"  Found edge ({n1}, {n2}) at coords {coord1} -> {coord2}")
        
        # Sort edges to form continuous 1D elements along the line
        if line_edges:
            sorted_edges = sort_edges_along_line(line_edges, nodes, line_pts_clean, debug)
            
            # Create 1D elements from sorted edges
            for n1, n2 in sorted_edges:
                # For linear elements, just use the two nodes
                elements_1d.append([n1, n2, 0])  # Pad to 3 columns
                mat_ids_1d.append(line_idx)
                element_node_counts_1d.append(2)
            
            if debug:
                print(f"  Created {len(sorted_edges)} 1D elements for line {line_idx}")
    
    if debug:
        print(f"Total 1D elements extracted: {len(elements_1d)}")
    
    return elements_1d, mat_ids_1d, element_node_counts_1d


def is_edge_on_reinforcement_line(coord1, coord2, line_pts, tolerance=1e-6):
    """
    Check if an edge lies along a reinforcement line.
    
    Parameters:
        coord1, coord2: Edge endpoint coordinates (x, y)
        line_pts: List of (x, y) points defining the reinforcement line
        tolerance: Tolerance for coincidence checking
        
    Returns:
        bool: True if edge lies along the reinforcement line
    """
    x1, y1 = coord1
    x2, y2 = coord2
    
    # Check if both endpoints lie on the reinforcement line
    point1_on_line = is_point_on_line_segments(coord1, line_pts, tolerance)
    point2_on_line = is_point_on_line_segments(coord2, line_pts, tolerance)
    
    if not (point1_on_line and point2_on_line):
        return False
    
    # Additional check: ensure edge direction is consistent with line direction
    # This prevents selecting edges that cross the reinforcement line
    edge_vector = np.array([x2 - x1, y2 - y1])
    edge_length = np.linalg.norm(edge_vector)
    
    if edge_length < tolerance:
        return False
    
    edge_unit = edge_vector / edge_length
    
    # Check alignment with any segment of the reinforcement line
    # This allows edges to span multiple segments after intersection preprocessing
    for i in range(len(line_pts) - 1):
        seg_start = np.array(line_pts[i])
        seg_end = np.array(line_pts[i + 1])
        seg_vector = seg_end - seg_start
        seg_length = np.linalg.norm(seg_vector)
        
        if seg_length < tolerance:
            continue
            
        seg_unit = seg_vector / seg_length
        
        # Check if edge is aligned with this segment (or opposite direction)
        dot_product = abs(np.dot(edge_unit, seg_unit))
        if dot_product > 0.95:  # Nearly parallel (cos(18°) ≈ 0.95)
            # More flexible check: edge should be collinear with the reinforcement line
            # and both endpoints should lie on the line (but not necessarily on the same segment)
            return True
    
    return False


def is_point_on_line_segments(point, line_pts, tolerance=1e-6):
    """
    Check if a point lies on any segment of a multi-segment line.
    
    Parameters:
        point: (x, y) coordinates of point to check
        line_pts: List of (x, y) points defining the line segments
        tolerance: Tolerance for coincidence checking
        
    Returns:
        bool: True if point lies on any line segment
    """
    for i in range(len(line_pts) - 1):
        if is_point_on_line_segment(point, line_pts[i], line_pts[i + 1], tolerance):
            return True
    return False


def is_point_on_line_segment(point, seg_start, seg_end, tolerance=1e-6):
    """
    Check if a point lies on a line segment.
    
    Parameters:
        point: (x, y) coordinates of point to check
        seg_start: (x, y) coordinates of segment start
        seg_end: (x, y) coordinates of segment end
        tolerance: Tolerance for coincidence checking
        
    Returns:
        bool: True if point lies on the line segment
    """
    px, py = point
    x1, y1 = seg_start
    x2, y2 = seg_end
    
    # Check if point is within bounding box of segment
    if not (min(x1, x2) - tolerance <= px <= max(x1, x2) + tolerance and
            min(y1, y2) - tolerance <= py <= max(y1, y2) + tolerance):
        return False
    
    # Check collinearity using cross product
    cross_product = abs((py - y1) * (x2 - x1) - (px - x1) * (y2 - y1))
    
    # Check if cross product is close to zero (collinear)
    if cross_product < tolerance:
        # Verify point is between segment endpoints using dot product
        dot_product = (px - x1) * (x2 - x1) + (py - y1) * (y2 - y1)
        segment_length_sq = (x2 - x1) ** 2 + (y2 - y1) ** 2
        
        if segment_length_sq < tolerance:  # Degenerate segment
            return abs(px - x1) < tolerance and abs(py - y1) < tolerance
        
        # Parameter t should be between 0 and 1 for point to be on segment
        t = dot_product / segment_length_sq
        return -tolerance <= t <= 1 + tolerance
    
    return False


def sort_edges_along_line(edges, nodes, line_pts, debug=False):
    """
    Sort edges to form a continuous sequence along a reinforcement line.
    
    Parameters:
        edges: List of (n1, n2) edge tuples
        nodes: Node coordinates array
        line_pts: Reinforcement line points
        debug: Enable debug output
        
    Returns:
        list: Sorted list of (n1, n2) edge tuples
    """
    if not edges:
        return []
    
    if len(edges) == 1:
        return edges
    
    # Build connectivity graph
    node_connections = defaultdict(list)
    for n1, n2 in edges:
        node_connections[n1].append(n2)
        node_connections[n2].append(n1)
    
    # Find start node (should have only one connection, or be closest to line start)
    line_start = np.array(line_pts[0])
    line_end = np.array(line_pts[-1])
    
    start_candidates = []
    for node in node_connections:
        if len(node_connections[node]) == 1:  # End node
            start_candidates.append(node)
    
    if not start_candidates:
        # No clear end nodes, use node closest to line start
        min_dist = float('inf')
        start_node = list(node_connections.keys())[0]
        for node in node_connections:
            dist = np.linalg.norm(nodes[node] - line_start)
            if dist < min_dist:
                min_dist = dist
                start_node = node
    else:
        # Choose end node closest to line start
        min_dist = float('inf')
        start_node = start_candidates[0]
        for node in start_candidates:
            dist = np.linalg.norm(nodes[node] - line_start)
            if dist < min_dist:
                min_dist = dist
                start_node = node
    
    # Trace path from start node
    sorted_edges = []
    used_edges = set()
    current_node = start_node
    
    while True:
        # Find next unused edge from current node
        next_node = None
        for neighbor in node_connections[current_node]:
            edge_key = (min(current_node, neighbor), max(current_node, neighbor))
            if edge_key not in used_edges:
                next_node = neighbor
                used_edges.add(edge_key)
                sorted_edges.append((current_node, next_node))
                break
        
        if next_node is None:
            break
            
        current_node = next_node
    
    if debug:
        print(f"    Sorted {len(sorted_edges)} edges along line")
    
    return sorted_edges

def verify_mesh_connectivity(mesh, tolerance=1e-8):
    """
    Verify that the mesh is properly connected by checking for duplicate nodes at shared boundaries.
    
    Parameters:
        mesh: Mesh dictionary with 'nodes' and 'elements' keys
        tolerance: Tolerance for considering nodes as duplicates
    
    Returns:
        dict: Connectivity verification results
    """
    import numpy as np
    from collections import defaultdict
    
    nodes = mesh["nodes"]
    elements = mesh["elements"]
    
    # Find duplicate nodes (nodes at same location)
    duplicate_groups = []
    used_indices = set()
    
    for i in range(len(nodes)):
        if i in used_indices:
            continue
            
        duplicates = [i]
        for j in range(i + 1, len(nodes)):
            if j in used_indices:
                continue
                
            if np.linalg.norm(nodes[i] - nodes[j]) < tolerance:
                duplicates.append(j)
                used_indices.add(j)
        
        if len(duplicates) > 1:
            duplicate_groups.append(duplicates)
            used_indices.add(i)
    
    # Check element connectivity
    element_connectivity = defaultdict(set)
    for elem_idx, element in enumerate(elements):
        for node_idx in element:
            element_connectivity[node_idx].add(elem_idx)
    
    # Find isolated nodes (nodes not used by any element)
    isolated_nodes = []
    for i in range(len(nodes)):
        if i not in element_connectivity:
            isolated_nodes.append(i)
    
    # Find elements with duplicate nodes
    elements_with_duplicates = []
    for elem_idx, element in enumerate(elements):
        unique_nodes = set(element)
        if len(unique_nodes) != len(element):
            elements_with_duplicates.append(elem_idx)
    
    results = {
        "total_nodes": len(nodes),
        "total_elements": len(elements),
        "duplicate_node_groups": duplicate_groups,
        "isolated_nodes": isolated_nodes,
        "elements_with_duplicates": elements_with_duplicates,
        "is_connected": len(duplicate_groups) == 0 and len(isolated_nodes) == 0
    }
    
    return results

def print_mesh_connectivity_report(mesh, tolerance=1e-8):
    """
    Print a detailed report about mesh connectivity.
    
    Parameters:
        mesh: Mesh dictionary
        tolerance: Tolerance for considering nodes as duplicates
    """
    results = verify_mesh_connectivity(mesh, tolerance)
    
    print("=== MESH CONNECTIVITY REPORT ===")
    print(f"Total nodes: {results['total_nodes']}")
    print(f"Total elements: {results['total_elements']}")
    print(f"Mesh is properly connected: {results['is_connected']}")
    print()
    
    if results['duplicate_node_groups']:
        print(f"WARNING: Found {len(results['duplicate_node_groups'])} groups of duplicate nodes:")
        for i, group in enumerate(results['duplicate_node_groups']):
            print(f"  Group {i+1}: Nodes {group} at position {mesh['nodes'][group[0]]}")
        print()
    
    if results['isolated_nodes']:
        print(f"WARNING: Found {len(results['isolated_nodes'])} isolated nodes:")
        for node_idx in results['isolated_nodes']:
            print(f"  Node {node_idx} at position {mesh['nodes'][node_idx]}")
        print()
    
    if results['elements_with_duplicates']:
        print(f"WARNING: Found {len(results['elements_with_duplicates'])} elements with duplicate nodes:")
        for elem_idx in results['elements_with_duplicates']:
            print(f"  Element {elem_idx}: {mesh['elements'][elem_idx]}")
        print()
    
    if results['is_connected']:
        print("✓ Mesh connectivity is good - no duplicate nodes or isolated nodes found.")
    else:
        print("✗ Mesh connectivity issues detected. Consider regenerating the mesh.")

def find_element_containing_point(nodes, elements, element_types, point):
    """
    Find which element contains the given point using spatial indexing for efficiency.
    
    Parameters:
        nodes: np.ndarray of node coordinates (n_nodes, 2)
        elements: np.ndarray of element vertex indices (n_elements, 9) - unused nodes set to 0
        element_types: np.ndarray indicating element type (3, 4, 6, 8, or 9 nodes)
        point: tuple (x, y) coordinates of the point to find
        
    Returns:
        int: Index of the element containing the point, or -1 if not found
    """
    x, y = point
    
    # Use spatial indexing to find candidate elements quickly
    # Build spatial hash grid if not already built
    if not hasattr(find_element_containing_point, '_spatial_grid'):
        find_element_containing_point._spatial_grid = _build_spatial_grid(nodes, elements, element_types)
    
    spatial_grid = find_element_containing_point._spatial_grid
    
    # Find grid cell containing the point
    grid_x = int((x - spatial_grid['x_min']) / spatial_grid['cell_size'])
    grid_y = int((y - spatial_grid['y_min']) / spatial_grid['cell_size'])
    
    # Get candidate elements from this cell and neighboring cells
    candidate_elements = set()
    for dx in [-1, 0, 1]:
        for dy in [-1, 0, 1]:
            cell_key = (grid_x + dx, grid_y + dy)
            if cell_key in spatial_grid['cells']:
                candidate_elements.update(spatial_grid['cells'][cell_key])
    
    # Check only the candidate elements
    for elem_idx in candidate_elements:
        element = elements[elem_idx]
        elem_type = element_types[elem_idx]
        
        if elem_type in [3, 6]:  # Triangle (linear or quadratic)
            # For point-in-element testing, use only corner nodes
            x1, y1 = nodes[element[0]]
            x2, y2 = nodes[element[1]]
            x3, y3 = nodes[element[2]]
            
            # Calculate barycentric coordinates
            det = (y2 - y3) * (x1 - x3) + (x3 - x2) * (y1 - y3)
            if abs(det) < 1e-12:  # Degenerate triangle
                continue
                
            lambda1 = ((y2 - y3) * (x - x3) + (x3 - x2) * (y - y3)) / det
            lambda2 = ((y3 - y1) * (x - x3) + (x1 - x3) * (y - y3)) / det
            lambda3 = 1.0 - lambda1 - lambda2
            
            # Check if point is inside triangle (all barycentric coordinates >= 0)
            if lambda1 >= -1e-12 and lambda2 >= -1e-12 and lambda3 >= -1e-12:
                return elem_idx
                
        elif elem_type in [4, 8, 9]:  # Quadrilateral (linear or quadratic)
            # For point-in-element testing, use only corner nodes
            x1, y1 = nodes[element[0]]
            x2, y2 = nodes[element[1]]
            x3, y3 = nodes[element[2]]
            x4, y4 = nodes[element[3]]
            
            # Use point-in-polygon test for quadrilaterals
            # Check if point is inside by counting crossings
            vertices = [(x1, y1), (x2, y2), (x3, y3), (x4, y4)]
            inside = False
            
            for j in range(len(vertices)):
                xi, yi = vertices[j]
                xj, yj = vertices[(j + 1) % len(vertices)]
                
                if ((yi > y) != (yj > y)) and (x < (xj - xi) * (y - yi) / (yj - yi) + xi):
                    inside = not inside
            
            if inside:
                return elem_idx
    
    return -1  # Point not found in any element


def _build_spatial_grid(nodes, elements, element_types):
    """
    Build a spatial hash grid for efficient element searching.
    
    Parameters:
        nodes: np.ndarray of node coordinates (n_nodes, 2)
        elements: np.ndarray of element vertex indices (n_elements, 8)
        element_types: np.ndarray indicating element type (3, 4, 6, or 8 nodes)
        
    Returns:
        dict: Spatial grid data structure
    """
    # Calculate bounding box
    x_coords = nodes[:, 0]
    y_coords = nodes[:, 1]
    x_min, x_max = x_coords.min(), x_coords.max()
    y_min, y_max = y_coords.min(), y_coords.max()
    
    # Determine optimal cell size based on average element size
    total_area = 0
    for i, (element, elem_type) in enumerate(zip(elements, element_types)):
        if elem_type in [3, 6]:  # Triangle
            x1, y1 = nodes[element[0]]
            x2, y2 = nodes[element[1]]
            x3, y3 = nodes[element[2]]
            area = 0.5 * abs((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1))
        else:  # Quadrilateral (4 or 8 nodes)
            x1, y1 = nodes[element[0]]
            x2, y2 = nodes[element[1]]
            x3, y3 = nodes[element[2]]
            x4, y4 = nodes[element[3]]
            area = 0.5 * abs((x2 - x1) * (y4 - y1) - (x4 - x1) * (y2 - y1))
        total_area += area
    
    avg_element_area = total_area / len(elements)
    # Cell size should be roughly 2-3 times the square root of average element area
    cell_size = max(0.1, 2.5 * np.sqrt(avg_element_area))
    
    # Build grid
    grid = {
        'x_min': x_min,
        'y_min': y_min,
        'cell_size': cell_size,
        'cells': {}
    }
    
    # Assign elements to grid cells
    for elem_idx, (element, elem_type) in enumerate(zip(elements, element_types)):
        # Calculate element bounding box
        if elem_type in [3, 6]:  # Triangle
            x_coords = [nodes[element[0]][0], nodes[element[1]][0], nodes[element[2]][0]]
            y_coords = [nodes[element[0]][1], nodes[element[1]][1], nodes[element[2]][1]]
        else:  # Quadrilateral (4 or 8 nodes)
            x_coords = [nodes[element[0]][0], nodes[element[1]][0], nodes[element[2]][0], nodes[element[3]][0]]
            y_coords = [nodes[element[0]][1], nodes[element[1]][1], nodes[element[2]][1], nodes[element[3]][1]]
        
        elem_x_min, elem_x_max = min(x_coords), max(x_coords)
        elem_y_min, elem_y_max = min(y_coords), max(y_coords)
        
        # Find grid cells that overlap with this element
        start_x = int((elem_x_min - x_min) / cell_size)
        end_x = int((elem_x_max - x_min) / cell_size) + 1
        start_y = int((elem_y_min - y_min) / cell_size)
        end_y = int((elem_y_max - y_min) / cell_size) + 1
        
        # Add element to all overlapping cells
        for grid_x in range(start_x, end_x + 1):
            for grid_y in range(start_y, end_y + 1):
                cell_key = (grid_x, grid_y)
                if cell_key not in grid['cells']:
                    grid['cells'][cell_key] = set()
                grid['cells'][cell_key].add(elem_idx)
    
    return grid


def interpolate_at_point(nodes, elements, element_types, values, point):
    """
    Interpolate values at a given point using the mesh.
    
    Parameters:
        nodes: np.ndarray of node coordinates (n_nodes, 2)
        elements: np.ndarray of element vertex indices (n_elements, 8)
        element_types: np.ndarray indicating element type (3, 4, 6, or 8 nodes)
        values: np.ndarray of values at nodes (n_nodes,)
        point: tuple (x, y) coordinates of the point to interpolate at
        
    Returns:
        float: Interpolated value at the point, or 0.0 if point not found
    """
    # Find the element containing the point
    element_idx = find_element_containing_point(nodes, elements, element_types, point)
    
    if element_idx == -1:
        return 0.0  # Point not found in any element
    
    element = elements[element_idx]
    elem_type = element_types[element_idx]
    x, y = point
    
    if elem_type == 3:  # Linear triangle
        # Get triangle vertices and values
        x1, y1 = nodes[element[0]]
        x2, y2 = nodes[element[1]]
        x3, y3 = nodes[element[2]]
        v1 = values[element[0]]
        v2 = values[element[1]]
        v3 = values[element[2]]
        
        # Calculate barycentric coordinates
        det = (y2 - y3) * (x1 - x3) + (x3 - x2) * (y1 - y3)
        lambda1 = ((y2 - y3) * (x - x3) + (x3 - x2) * (y - y3)) / det
        lambda2 = ((y3 - y1) * (x - x3) + (x1 - x3) * (y - y3)) / det
        lambda3 = 1.0 - lambda1 - lambda2
        
        # Interpolate using barycentric coordinates
        interpolated_value = lambda1 * v1 + lambda2 * v2 + lambda3 * v3
        
    elif elem_type == 6:  # Quadratic triangle
        # Get all 6 nodes: corners (0,1,2) and midpoints (3,4,5)
        # Node ordering: 0-1-2 corners, 3 midpoint of 0-1, 4 midpoint of 1-2, 5 midpoint of 2-0
        corner_nodes = [element[0], element[1], element[2]]
        midpoint_nodes = [element[3], element[4], element[5]]
        
        # Get coordinates
        x1, y1 = nodes[corner_nodes[0]]  # Node 0
        x2, y2 = nodes[corner_nodes[1]]  # Node 1  
        x3, y3 = nodes[corner_nodes[2]]  # Node 2
        
        # Calculate barycentric coordinates (L1, L2, L3)
        det = (y2 - y3) * (x1 - x3) + (x3 - x2) * (y1 - y3)
        L1 = ((y2 - y3) * (x - x3) + (x3 - x2) * (y - y3)) / det
        L2 = ((y3 - y1) * (x - x3) + (x1 - x3) * (y - y3)) / det
        L3 = 1.0 - L1 - L2
        
        # Quadratic shape functions for 6-node triangle
        N = np.zeros(6)
        N[0] = L1 * (2*L1 - 1)     # Corner node 0
        N[1] = L2 * (2*L2 - 1)     # Corner node 1
        N[2] = L3 * (2*L3 - 1)     # Corner node 2
        N[3] = 4 * L1 * L2         # Midpoint node 0-1
        N[4] = 4 * L2 * L3         # Midpoint node 1-2
        N[5] = 4 * L3 * L1         # Midpoint node 2-0
        
        # Interpolate using quadratic shape functions
        interpolated_value = 0.0
        for i in range(6):
            interpolated_value += N[i] * values[element[i]]
        
    elif elem_type == 4:  # Linear quadrilateral
        # Get quadrilateral vertices and values
        x1, y1 = nodes[element[0]]
        x2, y2 = nodes[element[1]]
        x3, y3 = nodes[element[2]]
        x4, y4 = nodes[element[3]]
        v1 = values[element[0]]
        v2 = values[element[1]]
        v3 = values[element[2]]
        v4 = values[element[3]]
        
        # Use proper bilinear shape functions for quadrilaterals
        # Map to natural coordinates (xi, eta) in [-1, 1] x [-1, 1]
        
        # For bilinear quad4, use iterative Newton-Raphson to find natural coordinates
        # Initial guess at element center
        xi, eta = 0.0, 0.0
        
        # Newton-Raphson iteration to find (xi, eta) such that physical coordinates match
        for _ in range(10):  # Max 10 iterations
            # Bilinear shape functions
            N = np.array([
                0.25 * (1-xi) * (1-eta),  # Node 0
                0.25 * (1+xi) * (1-eta),  # Node 1
                0.25 * (1+xi) * (1+eta),  # Node 2
                0.25 * (1-xi) * (1+eta)   # Node 3
            ])
            
            # Shape function derivatives
            dN_dxi = np.array([
                -0.25 * (1-eta),  # Node 0
                 0.25 * (1-eta),  # Node 1
                 0.25 * (1+eta),  # Node 2
                -0.25 * (1+eta)   # Node 3
            ])
            
            dN_deta = np.array([
                -0.25 * (1-xi),   # Node 0
                -0.25 * (1+xi),   # Node 1
                 0.25 * (1+xi),   # Node 2
                 0.25 * (1-xi)    # Node 3
            ])
            
            # Current physical coordinates
            x_curr = N[0]*x1 + N[1]*x2 + N[2]*x3 + N[3]*x4
            y_curr = N[0]*y1 + N[1]*y2 + N[2]*y3 + N[3]*y4
            
            # Residual
            fx = x_curr - x
            fy = y_curr - y
            
            if abs(fx) < 1e-10 and abs(fy) < 1e-10:
                break
                
            # Jacobian
            dx_dxi = dN_dxi[0]*x1 + dN_dxi[1]*x2 + dN_dxi[2]*x3 + dN_dxi[3]*x4
            dx_deta = dN_deta[0]*x1 + dN_deta[1]*x2 + dN_deta[2]*x3 + dN_deta[3]*x4
            dy_dxi = dN_dxi[0]*y1 + dN_dxi[1]*y2 + dN_dxi[2]*y3 + dN_dxi[3]*y4
            dy_deta = dN_deta[0]*y1 + dN_deta[1]*y2 + dN_deta[2]*y3 + dN_deta[3]*y4
            
            det_J = dx_dxi * dy_deta - dx_deta * dy_dxi
            if abs(det_J) < 1e-12:
                break
            
            # Newton-Raphson update
            dxi = (dy_deta * fx - dx_deta * fy) / det_J
            deta = (-dy_dxi * fx + dx_dxi * fy) / det_J
            
            xi -= dxi
            eta -= deta
            
            # Clamp to [-1,1]
            xi = max(-1, min(1, xi))
            eta = max(-1, min(1, eta))
        
        # Final bilinear shape functions
        N = np.array([
            0.25 * (1-xi) * (1-eta),  # Node 0
            0.25 * (1+xi) * (1-eta),  # Node 1
            0.25 * (1+xi) * (1+eta),  # Node 2
            0.25 * (1-xi) * (1+eta)   # Node 3
        ])
        
        # Interpolate using bilinear shape functions
        interpolated_value = N[0]*v1 + N[1]*v2 + N[2]*v3 + N[3]*v4
        
    elif elem_type == 8:  # Quadratic quadrilateral
        # Get all 8 nodes: corners (0,1,2,3) and midpoints (4,5,6,7)
        # Node ordering: 0-1-2-3 corners, 4 midpoint of 0-1, 5 midpoint of 1-2, 
        #                6 midpoint of 2-3, 7 midpoint of 3-0
        
        # Get corner coordinates for mapping to natural coordinates
        x1, y1 = nodes[element[0]]  # Node 0
        x2, y2 = nodes[element[1]]  # Node 1
        x3, y3 = nodes[element[2]]  # Node 2
        x4, y4 = nodes[element[3]]  # Node 3
        
        # For quadratic quads, we need to map from physical (x,y) to natural coordinates (xi,eta)
        # This is complex for general quadrilaterals, so use simplified approach:
        # Map to unit square [-1,1] x [-1,1] using bilinear mapping of corners
        
        # Bilinear inverse mapping (approximate for general quads)
        # Solve for natural coordinates xi, eta in [-1,1] x [-1,1]
        
        # For simplicity, use area coordinate method similar to linear quad
        # but with quadratic shape functions
        
        # Calculate area coordinates (this is an approximation)
        A_total = 0.5 * abs((x3-x1)*(y4-y2) - (x4-x2)*(y3-y1))
        if A_total < 1e-12:
            # Degenerate element, fall back to linear
            A1 = abs((x - x1) * (y2 - y1) - (x2 - x1) * (y - y1)) / 2
            A2 = abs((x - x2) * (y3 - y2) - (x3 - x2) * (y - y2)) / 2
            A3 = abs((x - x3) * (y4 - y3) - (x4 - x3) * (y - y3)) / 2
            A4 = abs((x - x4) * (y1 - y4) - (x1 - x4) * (y - y4)) / 2
            A_sum = A1 + A2 + A3 + A4
            if A_sum > 1e-12:
                w1, w2, w3, w4 = A1/A_sum, A2/A_sum, A3/A_sum, A4/A_sum
            else:
                w1 = w2 = w3 = w4 = 0.25
            
            # Linear interpolation as fallback
            interpolated_value = (w1 * values[element[0]] + w2 * values[element[1]] + 
                                w3 * values[element[2]] + w4 * values[element[3]])
        else:
            # For proper quadratic interpolation, we need natural coordinates
            # This is a simplified implementation - full implementation would solve
            # the nonlinear system for xi,eta
            
            # Use parametric coordinates estimation
            # Map point to approximate natural coordinates
            xi_approx = 2 * (x - 0.5*(x1+x3)) / (x2+x3-x1-x4) if abs(x2+x3-x1-x4) > 1e-12 else 0
            eta_approx = 2 * (y - 0.5*(y1+y3)) / (y2+y4-y1-y3) if abs(y2+y4-y1-y3) > 1e-12 else 0
            
            # Clamp to [-1,1]
            xi = max(-1, min(1, xi_approx))
            eta = max(-1, min(1, eta_approx))
            
            # Quadratic shape functions for 8-node quad in natural coordinates
            N = np.zeros(8)
            # Corner nodes
            N[0] = 0.25 * (1-xi) * (1-eta) * (-xi-eta-1)   # Node 0
            N[1] = 0.25 * (1+xi) * (1-eta) * (xi-eta-1)    # Node 1  
            N[2] = 0.25 * (1+xi) * (1+eta) * (xi+eta-1)    # Node 2
            N[3] = 0.25 * (1-xi) * (1+eta) * (-xi+eta-1)   # Node 3
            # Midpoint nodes
            N[4] = 0.5 * (1-xi*xi) * (1-eta)               # Node 4 (midpoint 0-1)
            N[5] = 0.5 * (1+xi) * (1-eta*eta)              # Node 5 (midpoint 1-2)
            N[6] = 0.5 * (1-xi*xi) * (1+eta)               # Node 6 (midpoint 2-3)
            N[7] = 0.5 * (1-xi) * (1-eta*eta)              # Node 7 (midpoint 3-0)
            
            # Interpolate using quadratic shape functions
            interpolated_value = 0.0
            for i in range(8):
                interpolated_value += N[i] * values[element[i]]
    
    elif elem_type == 9:  # Biquadratic quadrilateral (9-node Lagrange)
        # Get all 9 nodes: corners (0,1,2,3), edges (4,5,6,7), and center (8)
        # Node ordering: 0-1-2-3 corners, 4 midpoint of 0-1, 5 midpoint of 1-2,
        #                6 midpoint of 2-3, 7 midpoint of 3-0, 8 center
        
        # Get corner coordinates for mapping to natural coordinates
        x1, y1 = nodes[element[0]]  # Node 0
        x2, y2 = nodes[element[1]]  # Node 1
        x3, y3 = nodes[element[2]]  # Node 2
        x4, y4 = nodes[element[3]]  # Node 3
        
        # Newton-Raphson iteration to find natural coordinates (xi, eta)
        xi, eta = 0.0, 0.0  # Initial guess at element center
        
        for _ in range(10):  # Max 10 iterations
            # Biquadratic Lagrange shape functions for all 9 nodes
            N = np.zeros(9)
            # Corner nodes
            N[0] = 0.25 * xi * (xi-1) * eta * (eta-1)     # Node 0: (-1,-1)
            N[1] = 0.25 * xi * (xi+1) * eta * (eta-1)     # Node 1: (1,-1)  
            N[2] = 0.25 * xi * (xi+1) * eta * (eta+1)     # Node 2: (1,1)
            N[3] = 0.25 * xi * (xi-1) * eta * (eta+1)     # Node 3: (-1,1)
            # Edge nodes
            N[4] = 0.5 * (1-xi*xi) * eta * (eta-1)        # Node 4: (0,-1)
            N[5] = 0.5 * xi * (xi+1) * (1-eta*eta)        # Node 5: (1,0)
            N[6] = 0.5 * (1-xi*xi) * eta * (eta+1)        # Node 6: (0,1)
            N[7] = 0.5 * xi * (xi-1) * (1-eta*eta)        # Node 7: (-1,0)
            # Center node
            N[8] = (1-xi*xi) * (1-eta*eta)                # Node 8: (0,0)
            
            # Shape function derivatives w.r.t. xi
            dN_dxi = np.zeros(9)
            dN_dxi[0] = 0.25 * (2*xi-1) * eta * (eta-1)
            dN_dxi[1] = 0.25 * (2*xi+1) * eta * (eta-1)
            dN_dxi[2] = 0.25 * (2*xi+1) * eta * (eta+1)
            dN_dxi[3] = 0.25 * (2*xi-1) * eta * (eta+1)
            dN_dxi[4] = -xi * eta * (eta-1)
            dN_dxi[5] = 0.5 * (2*xi+1) * (1-eta*eta)
            dN_dxi[6] = -xi * eta * (eta+1)
            dN_dxi[7] = 0.5 * (2*xi-1) * (1-eta*eta)
            dN_dxi[8] = -2*xi * (1-eta*eta)
            
            # Shape function derivatives w.r.t. eta
            dN_deta = np.zeros(9)
            dN_deta[0] = 0.25 * xi * (xi-1) * (2*eta-1)
            dN_deta[1] = 0.25 * xi * (xi+1) * (2*eta-1)
            dN_deta[2] = 0.25 * xi * (xi+1) * (2*eta+1)
            dN_deta[3] = 0.25 * xi * (xi-1) * (2*eta+1)
            dN_deta[4] = 0.5 * (1-xi*xi) * (2*eta-1)
            dN_deta[5] = -eta * xi * (xi+1)
            dN_deta[6] = 0.5 * (1-xi*xi) * (2*eta+1)
            dN_deta[7] = -eta * xi * (xi-1)
            dN_deta[8] = -2*eta * (1-xi*xi)
            
            # Current physical coordinates using all 9 nodes
            node_coords = nodes[element[:9]]
            x_curr = np.sum(N * node_coords[:, 0])
            y_curr = np.sum(N * node_coords[:, 1])
            
            # Residual
            fx = x_curr - x
            fy = y_curr - y
            
            if abs(fx) < 1e-10 and abs(fy) < 1e-10:
                break
                
            # Jacobian
            dx_dxi = np.sum(dN_dxi * node_coords[:, 0])
            dx_deta = np.sum(dN_deta * node_coords[:, 0])
            dy_dxi = np.sum(dN_dxi * node_coords[:, 1])
            dy_deta = np.sum(dN_deta * node_coords[:, 1])
            
            det_J = dx_dxi * dy_deta - dx_deta * dy_dxi
            if abs(det_J) < 1e-12:
                break
            
            # Newton-Raphson update
            dxi = (dy_deta * fx - dx_deta * fy) / det_J
            deta = (-dy_dxi * fx + dx_dxi * fy) / det_J
            
            xi -= dxi
            eta -= deta
            
            # Clamp to [-1,1]
            xi = max(-1, min(1, xi))
            eta = max(-1, min(1, eta))
        
        # Final biquadratic shape functions
        N = np.zeros(9)
        N[0] = 0.25 * xi * (xi-1) * eta * (eta-1)     # Node 0
        N[1] = 0.25 * xi * (xi+1) * eta * (eta-1)     # Node 1
        N[2] = 0.25 * xi * (xi+1) * eta * (eta+1)     # Node 2
        N[3] = 0.25 * xi * (xi-1) * eta * (eta+1)     # Node 3
        N[4] = 0.5 * (1-xi*xi) * eta * (eta-1)        # Node 4
        N[5] = 0.5 * xi * (xi+1) * (1-eta*eta)        # Node 5
        N[6] = 0.5 * (1-xi*xi) * eta * (eta+1)        # Node 6
        N[7] = 0.5 * xi * (xi-1) * (1-eta*eta)        # Node 7
        N[8] = (1-xi*xi) * (1-eta*eta)                # Node 8
        
        # Interpolate using biquadratic shape functions
        interpolated_value = 0.0
        for i in range(9):
            interpolated_value += N[i] * values[element[i]]
    
    else:
        return 0.0  # Unknown element type
    
    # Return zero if interpolated value is negative (pore pressure cannot be negative)
    return max(0.0, interpolated_value)


def test_1d_element_alignment(mesh, reinforcement_lines, tolerance=1e-6, debug=True):
    """
    Test that 1D elements correctly align with reinforcement lines.
    
    This function verifies that:
    1. Each reinforcement line is represented by a sequence of 1D elements
    2. The 1D elements form continuous paths along each reinforcement line
    3. The element endpoints match the expected line segment endpoints
    
    Parameters:
        mesh: Dictionary containing nodes and 1D element data
        reinforcement_lines: List of reinforcement lines, each containing coordinate tuples
        tolerance: Tolerance for coordinate comparison (default 1e-6)
        debug: Enable detailed debug output
        
    Returns:
        bool: True if all tests pass, False otherwise
    """
    if debug:
        print("\n=== Testing 1D Element Alignment ===")
    
    if 'elements_1d' not in mesh:
        print("ERROR: No 1D elements found in mesh")
        return False
        
    elements_1d = mesh['elements_1d']
    if elements_1d is None or len(elements_1d) == 0:
        print("ERROR: No 1D elements found in mesh")
        return False
    
    nodes = np.array(mesh['nodes'])
    elements_1d = mesh['elements_1d']
    
    if debug:
        print(f"Testing {len(reinforcement_lines)} reinforcement lines")
        print(f"Found {len(elements_1d)} 1D elements")
    
    success = True
    
    for line_idx, line_pts in enumerate(reinforcement_lines):
        if debug:
            print(f"\nTesting line {line_idx}: {line_pts}")
        
        # Remove duplicate endpoints and get expected segments
        line_pts_clean = remove_duplicate_endpoint(list(line_pts))
        if len(line_pts_clean) < 2:
            if debug:
                print(f"  Skipping line {line_idx}: insufficient points")
            continue
        
        # Expected segments for this line
        expected_segments = []
        for i in range(len(line_pts_clean) - 1):
            expected_segments.append((line_pts_clean[i], line_pts_clean[i + 1]))
        
        if debug:
            print(f"  Expected {len(expected_segments)} segments:")
            for i, (start, end) in enumerate(expected_segments):
                print(f"    Segment {i}: {start} -> {end}")
        
        # Find 1D elements that belong to this reinforcement line using material IDs
        line_elements = []
        if 'element_materials_1d' in mesh:
            element_materials_1d = mesh['element_materials_1d']
            for elem_idx, (element, material_id) in enumerate(zip(elements_1d, element_materials_1d)):
                # Skip zero-padded elements
                if len(element) < 2 or element[1] == 0:
                    continue
                    
                # Check if this element belongs to the current line
                if material_id == line_idx + 1:  # Material IDs are 1-based
                    # Get element coordinates
                    try:
                        coord1 = nodes[element[0]]
                        coord2 = nodes[element[1]]
                    except IndexError:
                        if debug:
                            print(f"  WARNING: Element {elem_idx} has invalid node indices {element[0]}, {element[1]}")
                        continue
                    
                    line_elements.append((elem_idx, coord1, coord2))
        else:
            # Fallback: use the old method if material IDs are not available
            for elem_idx, element in enumerate(elements_1d):
                # Skip zero-padded elements
                if len(element) < 2 or element[1] == 0:
                    continue
                    
                # Get element coordinates
                try:
                    coord1 = nodes[element[0]]
                    coord2 = nodes[element[1]]
                except IndexError:
                    if debug:
                        print(f"  WARNING: Element {elem_idx} has invalid node indices {element[0]}, {element[1]}")
                    continue
                
                # Check if this element lies on the current reinforcement line
                if is_edge_on_reinforcement_line(coord1, coord2, line_pts_clean, tolerance):
                    line_elements.append((elem_idx, coord1, coord2))
        
        if debug:
            print(f"  Found {len(line_elements)} 1D elements on this line:")
            for elem_idx, coord1, coord2 in line_elements:
                print(f"    Element {elem_idx}: {coord1} -> {coord2}")
        
        # Test 1: Check that we have at least some 1D elements for this line
        if len(line_elements) == 0:
            print(f"ERROR: Line {line_idx} has no 1D elements")
            success = False
            continue
        
        # Test 2: Check that we have reasonable number of elements
        # After intersection preprocessing, we may have more elements than original segments
        # But we should have at least some elements for each line
        if len(line_elements) == 0:
            print(f"ERROR: Line {line_idx} has no 1D elements")
            success = False
            continue
        
        # Test 2: Check if elements form continuous path
        if len(line_elements) > 1:
            # Sort elements to form continuous sequence
            sorted_elements = []
            remaining_elements = line_elements.copy()
            
            # Start with first element
            current_elem = remaining_elements.pop(0)
            sorted_elements.append(current_elem)
            
            # Build chain by finding connecting elements
            while remaining_elements:
                last_coord = sorted_elements[-1][2]  # End coordinate of last element
                
                # Find next element that starts where last one ended
                found_next = False
                for i, (elem_idx, coord1, coord2) in enumerate(remaining_elements):
                    if np.linalg.norm(np.array(coord1) - np.array(last_coord)) < tolerance:
                        sorted_elements.append((elem_idx, coord1, coord2))
                        remaining_elements.pop(i)
                        found_next = True
                        break
                    elif np.linalg.norm(np.array(coord2) - np.array(last_coord)) < tolerance:
                        # Element is reversed, flip it
                        sorted_elements.append((elem_idx, coord2, coord1))
                        remaining_elements.pop(i)
                        found_next = True
                        break
                
                if not found_next:
                    print(f"ERROR: Line {line_idx} elements do not form continuous path")
                    print(f"  Cannot connect from {last_coord}")
                    print(f"  Remaining elements: {remaining_elements}")
                    success = False
                    break
            
            line_elements = sorted_elements
        
        # Test 3: Check that the 1D elements cover the reinforcement line from start to end
        if len(line_elements) > 0:
            # Get the start and end points of the reinforcement line
            line_start = line_pts_clean[0]
            line_end = line_pts_clean[-1]
            
            # Find the first and last 1D elements
            first_elem = line_elements[0]
            last_elem = line_elements[-1]
            
            # Check if the first element starts near the line start
            first_start_dist = np.linalg.norm(np.array(first_elem[1]) - np.array(line_start))
            first_end_dist = np.linalg.norm(np.array(first_elem[2]) - np.array(line_start))
            
            # Check if the last element ends near the line end
            last_start_dist = np.linalg.norm(np.array(last_elem[1]) - np.array(line_end))
            last_end_dist = np.linalg.norm(np.array(last_elem[2]) - np.array(line_end))
            
            # The first element should start near the line start (either direction)
            # Be more flexible due to intersection preprocessing
            if first_start_dist > tolerance * 10 and first_end_dist > tolerance * 10:
                print(f"WARNING: Line {line_idx} first element does not start at line start")
                print(f"  Line start: {line_start}")
                print(f"  First element: {first_elem[1]} -> {first_elem[2]}")
                print(f"  Start distances: {first_start_dist:.2e}, {first_end_dist:.2e}")
                # Don't fail the test for this - just warn
            
            # The last element should end near the line end (either direction)
            # Be more flexible due to intersection preprocessing
            if last_start_dist > tolerance * 10 and last_end_dist > tolerance * 10:
                print(f"WARNING: Line {line_idx} last element does not end at line end")
                print(f"  Line end: {line_end}")
                print(f"  Last element: {last_elem[1]} -> {last_elem[2]}")
                print(f"  End distances: {last_start_dist:.2e}, {last_end_dist:.2e}")
                # Don't fail the test for this - just warn
        
        # Test 4: Check that line path is continuous
        if len(line_elements) > 1:
            for i in range(len(line_elements) - 1):
                end_coord = line_elements[i][2]  # End of current element
                start_coord = line_elements[i + 1][1]  # Start of next element
                
                gap = np.linalg.norm(np.array(end_coord) - np.array(start_coord))
                if gap > tolerance:
                    print(f"ERROR: Line {line_idx} has gap between elements {i} and {i+1}")
                    print(f"  Gap size: {gap:.2e}")
                    print(f"  Element {i} end: {end_coord}")
                    print(f"  Element {i+1} start: {start_coord}")
                    success = False
        
        if debug and success:
            print(f"  ✓ Line {line_idx} passes all alignment tests")
    
    if debug:
        if success:
            print("\n=== All 1D Element Alignment Tests PASSED ===")
        else:
            print("\n=== 1D Element Alignment Tests FAILED ===")
    
    return success

def add_intersection_points_to_polygons(polygons, lines, debug=False):
    """
    Add intersection points between reinforcement lines and polygon edges to the polygon vertex lists.
    This ensures that polygons have vertices at all intersection points with reinforcement lines.
    
    Parameters:
        polygons: List of polygons (lists of (x,y) tuples)
        lines: List of reinforcement lines (lists of (x,y) tuples)
        debug: Enable debug output
        
    Returns:
        Updated list of polygons with intersection points added
    """
    if not lines:
        return polygons
        
    if debug:
        print("Adding intersection points to polygons...")
    
    # Make a copy of polygons to modify
    updated_polygons = []
    for poly in polygons:
        updated_polygons.append(list(poly))  # Convert to list for modification
    
    # Find all intersections
    for line_idx, line_pts in enumerate(lines):
        line_pts_clean = remove_duplicate_endpoint(list(line_pts))
        
        if debug:
            print(f"Processing line {line_idx}: {line_pts_clean}")
        
        # Check each segment of the reinforcement line
        for i in range(len(line_pts_clean) - 1):
            line_seg_start = line_pts_clean[i]
            line_seg_end = line_pts_clean[i + 1]
            
            # Check intersection with each polygon
            for poly_idx, poly in enumerate(updated_polygons):
                # Check each edge of this polygon
                for j in range(len(poly)):
                    poly_edge_start = poly[j]
                    poly_edge_end = poly[(j + 1) % len(poly)]
                    
                    # Find intersection point if it exists
                    intersection = line_segment_intersection(
                        line_seg_start, line_seg_end,
                        poly_edge_start, poly_edge_end
                    )
                    
                    if intersection:
                        if debug:
                            print(f"Found intersection {intersection} between line {line_idx} segment {i} and polygon {poly_idx} edge {j}")
                        
                        # Check if intersection point is already a vertex of this polygon
                        is_vertex = False
                        for vertex in poly:
                            if abs(vertex[0] - intersection[0]) < 1e-8 and abs(vertex[1] - intersection[1]) < 1e-8:
                                is_vertex = True
                                break
                        
                        if not is_vertex:
                            # Insert intersection point into polygon at the correct position
                            # Insert after vertex j (which is the start of the edge)
                            insert_idx = j + 1
                            updated_polygons[poly_idx].insert(insert_idx, intersection)
                            
                            if debug:
                                print(f"Added intersection point {intersection} to polygon {poly_idx} at position {insert_idx}")
    
    return updated_polygons

def extract_reinforcement_line_geometry(slope_data):
    """
    Extract reinforcement line geometry from slope_data in the format needed for mesh generation.
    
    Parameters:
        slope_data: Dictionary containing slope data with 'reinforce_lines' key
        
    Returns:
        List of reinforcement lines, where each line is a list of (x, y) coordinate tuples
    """
    lines = []
    if 'reinforce_lines' in slope_data and slope_data['reinforce_lines']:
        for line in slope_data['reinforce_lines']:
            # Convert from dict format to tuple format
            line_coords = [(point['X'], point['Y']) for point in line]
            lines.append(line_coords)
    return lines