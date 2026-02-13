from xslope.fem import build_fem_data, solve_fem, solve_ssrm
from xslope.fileio import load_slope_data
from xslope.mesh import build_polygons, build_mesh_from_polygons
from xslope.plot import plot_inputs
from xslope.plot_fem import plot_fem_results, plot_fem_data

# Load Griffiths & Lane (1999) Example 1: homogeneous slope, c/γH=0.05, φ=20°
slope_data = load_slope_data("docs/fem/files/xslope_griffiths1_load.xlsx")

plot_inputs(slope_data, mode='fem', tab_loc='top')

# Build mesh
polygons = build_polygons(slope_data)
mesh = build_mesh_from_polygons(polygons, target_size=2, element_type='quad8')

# Build FEM data
fem_data = build_fem_data(slope_data, mesh)
plot_fem_data(fem_data, figsize=(14, 7), show_nodes=True, show_bc=True,
              label_elements=False, label_nodes=False)

analysis_type = "single" # @param ["single","ssrm"]
F = 1.36
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

