from pathlib import Path
from xslope.fem import build_fem_data, solve_fem, solve_ssrm
from xslope.fileio import load_slope_data
from xslope.mesh import build_polygons, build_mesh_from_polygons, export_mesh_to_json, import_mesh_from_json
from xslope.plot import plot_inputs
from xslope.plot_fem import plot_fem_results, plot_fem_data

# Load Griffiths & Lane (1999) Example 1: homogeneous slope, c/γH=0.05, φ=20°
input_file = "docs/fem/files/xslope_griffiths1.xlsx"
slope_data = load_slope_data(input_file)

plot_inputs(slope_data, mode='fem', tab_loc='top')

# Check to see if the mesh file already exists
input_path = Path(input_file)
mesh_file = input_path.parent / f"{input_path.stem}_mesh.json"

re_mesh = False
target_size = 2
element_type='quad8'
save_mesh = True

# Build the mesh if it doesn't exist or if re_mesh is True
if re_mesh or not mesh_file.exists():
    polygons = build_polygons(slope_data)
    print(f"Building mesh with {len(polygons)} polygons.")
    mesh = build_mesh_from_polygons(polygons, target_size=target_size, element_type=element_type)
    if save_mesh:
        export_mesh_to_json(mesh, mesh_file)
else:
    mesh = import_mesh_from_json(mesh_file)

# Build FEM data
fem_data = build_fem_data(slope_data, mesh)
plot_fem_data(fem_data, figsize=(14, 7), show_nodes=True, show_bc=True,
              label_elements=False, label_nodes=False)

analysis_type = "single" # @param ["single","ssrm"]
F = 1.25
F_min=0.5
F_max=1.5

if analysis_type == "single":
    solution = solve_fem(fem_data, F=F, debug_level=2)
    print(f"  Converged: {solution['converged']}, Iterations: {solution['iterations']}")
    plot_fem_results(fem_data, solution, plot_type=['deformation', 'shear_strain', 'displace_vector'], save_png=True)
elif analysis_type == "ssrm":
    result = solve_ssrm(fem_data, F_min=F_min, F_max=F_max, tolerance=0.05, debug_level=1)
    if result.get("converged", False):
        print(f"\nFactor of Safety: {result['FS']:.2f}")
        print(f"Method: {result.get('method', 'Unknown')}")
        plot_fem_results(fem_data, result['last_solution'],
                         plot_type=['deformation', 'shear_strain', 'displace_vector'], save_png=True)
    else:
        print(f"SSRM failed: {result.get('error', 'Unknown error')}")

