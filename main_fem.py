from pathlib import Path
from xslope.fem import build_fem_data, solve_fem, solve_ssrm, print_reinforcement_summary, print_pile_summary, print_detailed_element_summary, export_fem_solution
from xslope.fileio import load_slope_data
from xslope.mesh import get_material_polygons, build_mesh_from_polygons, export_mesh_to_json, extract_constraint_line_geometry
from xslope.plot import plot_inputs
from xslope.plot_fem import plot_fem_results, plot_fem_data

input_file = "../ce544/docs/unit2/12_fem/files/xslope_earth_dam_fem_KEY.xlsx"
slope_data = load_slope_data(input_file)

plot_inputs(slope_data, mode='fem', tab_loc='top', save_png=True)

input_path = Path(input_file)
auto_size = True # @param {"type":"boolean"}
size_divisions = 50 # @param {"type":"number"}
target_size = 1 # Desired target element size for meshing (adjust as needed for finer/coarser mesh)
element_type = 'tri6' # @param ["quad4","quad8","tri3","tri6"]
remesh = False # @param {"type":"boolean"}

# Use existing mesh from slope_data if available, otherwise build a new one
if slope_data.get("mesh") is not None and not remesh:
    print("Using existing mesh file.")
    mesh = slope_data["mesh"]
else:
    print("No existing mesh found in slope_data or remeshing enabled, building new mesh from profile line data.")
    constraint_lines, n_reinf, n_pile = extract_constraint_line_geometry(slope_data)
    polygons = get_material_polygons(slope_data, reinf_lines=constraint_lines)
    print(f"Building mesh with {len(polygons)} polygons, {n_reinf} reinforcement lines, {n_pile} pile lines.")

    if auto_size:
        # find the x-range of the ground_surface and use it to set the target size
        x_range = [min(x for x, _ in slope_data['ground_surface'].coords), max(x for x, _ in slope_data['ground_surface'].coords)]
        target_size = (x_range[1] - x_range[0]) / size_divisions
        print(f"Auto-calculated target element size: {target_size:.3f}")

    mesh = build_mesh_from_polygons(polygons, target_size=target_size, element_type=element_type, lines=constraint_lines)
    mesh_file = input_path.parent / f"{input_path.stem}_mesh.json"
    export_mesh_to_json(mesh, mesh_file)

# Build FEM data
fem_data = build_fem_data(slope_data, mesh)
plot_fem_data(fem_data, figsize=(14, 7), show_nodes=True, show_bc=True,
              label_elements=False, label_nodes=False, save_png=True)

analysis_type = "ssrm" # @param ["single","ssrm"]
failure_criterion = "non_convergence" # @param ["non_convergence","displacement_limit","displacement_increase","unbalanced_force"]
deform_percent = 15 # @param {"type":"number"} for plotting deformation results - percentage of slope height

F = 1.14     # Initial guess for Factor of Safety (used for single analysis) - adjust as needed
F_min=1.4  # Minimum FS for SSRM search (adjust as needed)
F_max=1.6   # Maximum FS for SSRM search (adjust as needed)

if analysis_type == "single":
    solution = solve_fem(fem_data, F=F, debug_level=2)
    print(f"  Converged: {solution['converged']}, Iterations: {solution['iterations']}")
    print_reinforcement_summary(fem_data, solution)
    print_pile_summary(fem_data, solution)
    print_detailed_element_summary(fem_data, solution)
    plot_fem_results(fem_data, solution, plot_type=['deformation', 'shear_strain', 'displace_vector'], deform_percent=deform_percent, save_png=True)
    export_fem_solution(fem_data, solution, input_path.parent / input_path.stem)
elif analysis_type == "ssrm":
    result = solve_ssrm(fem_data, F_min=F_min, F_max=F_max, tolerance=0.05, debug_level=1,
                        failure_criterion=failure_criterion)
    if result.get("converged", False):
        print(f"\nFactor of Safety: {result['FS']:.2f}")
        print(f"Method: {result.get('method', 'Unknown')}")
        print_reinforcement_summary(fem_data, result['last_solution'])
        print_pile_summary(fem_data, result['last_solution'])
        print_detailed_element_summary(fem_data, result['last_solution'])
        plot_fem_results(fem_data, result['last_solution'],
                         plot_type=['deformation', 'shear_strain', 'displace_vector'], deform_percent=deform_percent, save_png=True)
        export_fem_solution(fem_data, result['last_solution'], input_path.parent / input_path.stem)
    else:
        print(f"SSRM failed: {result.get('error', 'Unknown error')}")
