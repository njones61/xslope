from pathlib import Path
from xslope.fem import build_fem_data, solve_fem, solve_ssrm, print_reinforcement_summary
from xslope.fileio import load_slope_data
from xslope.mesh import build_polygons, build_mesh_from_polygons, export_mesh_to_json, extract_reinforcement_line_geometry
from xslope.plot import plot_inputs
from xslope.plot_fem import plot_fem_results, plot_fem_data

# Load Griffiths & Lane (1999) Example 1: homogeneous slope, c/γH=0.05, φ=20°
input_file = "docs/fem/files/xslope_reinforce_fem.xlsx"
slope_data = load_slope_data(input_file)

plot_inputs(slope_data, mode='fem', tab_loc='top')

input_path = Path(input_file)
target_size = 2 # Desired target element size for meshing (adjust as needed for finer/coarser mesh)
element_type = 'tri6' # @param ["quad4","quad8","tri3","tri6"]
remesh = True # @param {"type":"boolean"}

# Use existing mesh from slope_data if available, otherwise build a new one
if slope_data.get("mesh") is not None and not remesh:
    print("Using existing mesh file.")
    mesh = slope_data["mesh"]
else:
    print("No existing mesh found in slope_data or remeshing enabled, building new mesh from profile line data.")
    reinf_geom = extract_reinforcement_line_geometry(slope_data)
    polygons = build_polygons(slope_data, reinf_lines=reinf_geom)
    print(f"Building mesh with {len(polygons)} polygons, {len(reinf_geom)} reinforcement lines.")
    mesh = build_mesh_from_polygons(polygons, target_size=target_size, element_type=element_type, lines=reinf_geom)
    mesh_file = input_path.parent / f"{input_path.stem}_mesh.json"
    export_mesh_to_json(mesh, mesh_file)

# Build FEM data
fem_data = build_fem_data(slope_data, mesh)
plot_fem_data(fem_data, figsize=(14, 7), show_nodes=True, show_bc=True,
              label_elements=False, label_nodes=False, save_png=True)

analysis_type = "ssrm" # @param ["single","ssrm"]
failure_criterion = "non_convergence" # @param ["non_convergence","displacement_limit","displacement_increase","unbalanced_force"]

F = 1.55     # Initial guess for Factor of Safety (used for single analysis) - adjust as needed
F_min=1.2   # Minimum FS for SSRM search (adjust as needed)
F_max=1.8   # Maximum FS for SSRM search (adjust as needed)

if analysis_type == "single":
    solution = solve_fem(fem_data, F=F, debug_level=2)
    print(f"  Converged: {solution['converged']}, Iterations: {solution['iterations']}")
    print_reinforcement_summary(fem_data, solution)
    plot_fem_results(fem_data, solution, plot_type=['deformation', 'shear_strain', 'displace_vector'], save_png=True)
elif analysis_type == "ssrm":
    result = solve_ssrm(fem_data, F_min=F_min, F_max=F_max, tolerance=0.05, debug_level=1,
                        failure_criterion=failure_criterion)
    if result.get("converged", False):
        print(f"\nFactor of Safety: {result['FS']:.2f}")
        print(f"Method: {result.get('method', 'Unknown')}")
        print_reinforcement_summary(fem_data, result['last_solution'])
        plot_fem_results(fem_data, result['last_solution'],
                         plot_type=['deformation', 'shear_strain', 'displace_vector'], save_png=True)
    else:
        print(f"SSRM failed: {result.get('error', 'Unknown error')}")

