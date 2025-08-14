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

"""
Main script for testing finite element slope stability analysis.

This script demonstrates the use of the FEM implementation for slope stability
analysis using the Shear Strength Reduction Method (SSRM). It can work with
existing XSLOPE input files or create simple test cases.
"""

from fem import build_fem_data, solve_fem, solve_ssrm
from plot_fem import plot_fem_results, plot_reinforcement_force_profiles, plot_ssrm_convergence
from fileio import load_slope_data
from mesh import build_mesh_from_polygons
import matplotlib.pyplot as plt
import numpy as np
import sys
import os


def create_simple_slope_data():
    """
    Create a simple slope geometry for testing FEM implementation.
    
    This creates a basic left-facing slope with a single material and
    simple boundary conditions for testing purposes.
    """
    
    # Define a simple slope geometry
    # Ground surface from x=0 to x=20, with slope from y=10 to y=5
    ground_surface = [
        (0, 10),
        (10, 10),
        (20, 5)
    ]
    
    # Define slope boundaries - single material region
    slope_boundary = [
        (0, 0),    # Bottom left
        (20, 0),   # Bottom right  
        (20, 5),   # Top right (toe of slope)
        (10, 10),  # Top of slope
        (0, 10),   # Top left
        (0, 0)     # Back to start
    ]
    
    # Material properties
    materials = [
        {
            "name": "Clay",
            "strength_option": "mc",
            "c": 20000,      # Pa (20 kPa)
            "phi": 20,       # degrees
            "E": 50e6,       # Pa (50 MPa)
            "nu": 0.3,       # Poisson's ratio
            "gamma": 18000,  # N/m3 (18 kN/m3)
            "pp_option": "none"
        }
    ]
    
    # Simple slope data structure
    slope_data = {
        "ground_surface": ground_surface,
        "profile_lines": [slope_boundary],
        "materials": materials,
        "gamma_water": 9810,  # N/m3
        "k_seismic": 0.0,
        "max_depth": 0.0,
        "distributed_loads": [],
        "reinforcement_lines": []
    }
    
    return slope_data


def create_simple_slope_with_reinforcement():
    """
    Create a simple slope with reinforcement for testing.
    """
    
    # Start with basic slope
    slope_data = create_simple_slope_data()
    
    # Add reinforcement lines
    slope_data["reinforcement_lines"] = [
        {
            "x1": 2, "y1": 8,    # Start point
            "x2": 18, "y2": 2,   # End point
            "t_max": 50000,      # N (50 kN max tensile force)
            "t_res": 25000,      # N (25 kN residual force)
            "lp1": 2.0,          # m (pullout length left)
            "lp2": 2.0,          # m (pullout length right)
            "E": 200e9,          # Pa (200 GPa - steel)
            "area": 1e-4         # m2 (1 cm2)
        },
        {
            "x1": 4, "y1": 9,
            "x2": 16, "y2": 3,
            "t_max": 40000,
            "t_res": 20000,
            "lp1": 2.0,
            "lp2": 2.0,
            "E": 200e9,
            "area": 8e-5
        }
    ]
    
    return slope_data


def test_simple_fem_analysis():
    """
    Test basic FEM analysis without SSRM.
    """
    print("=== Testing Basic FEM Analysis ===")
    
    # Create simple slope
    slope_data = create_simple_slope_data()
    
    # Generate mesh
    polygons = [slope_data["profile_lines"][0]]
    mesh = build_mesh_from_polygons(polygons, target_size=2.0, element_type='tri3')
    
    print(f"Mesh created: {len(mesh['nodes'])} nodes, {len(mesh['elements'])} elements")
    
    # Build FEM data
    fem_data = build_fem_data(slope_data, mesh)
    
    print(f"FEM data built successfully")
    print(f"Boundary conditions: {np.sum(fem_data['bc_type'] != 0)} constrained nodes")
    
    # Solve FEM with F=1.0 (no strength reduction)
    print("\nSolving FEM with F=1.0...")
    solution = solve_fem(fem_data, F=1.0, debug_level=1)
    
    if solution["converged"]:
        print(f"✓ FEM solution converged in {solution['iterations']} iterations")
        max_disp = np.max(np.abs(solution["displacements"]))
        print(f"Maximum displacement: {max_disp:.6f} m")
        
        n_plastic = np.sum(solution["plastic_elements"])
        n_total = len(solution["plastic_elements"])
        print(f"Plastic elements: {n_plastic}/{n_total} ({100*n_plastic/n_total:.1f}%)")
        
        # Plot results
        fig, ax = plot_fem_results(fem_data, solution, plot_type='displacement')
        plt.savefig('fem_test_displacement.png', dpi=150, bbox_inches='tight')
        print("Displacement plot saved as 'fem_test_displacement.png'")
        
        fig, ax = plot_fem_results(fem_data, solution, plot_type='stress')
        plt.savefig('fem_test_stress.png', dpi=150, bbox_inches='tight')
        print("Stress plot saved as 'fem_test_stress.png'")
        
        plt.show()
        
        return True
    else:
        print("✗ FEM solution failed to converge")
        print(f"Error: {solution.get('error', 'Unknown error')}")
        return False


def test_fem_with_reinforcement():
    """
    Test FEM analysis with reinforcement elements.
    """
    print("\n=== Testing FEM with Reinforcement ===")
    
    # Create slope with reinforcement
    slope_data = create_simple_slope_with_reinforcement()
    
    # Generate mesh including reinforcement lines
    polygons = [slope_data["profile_lines"][0]]
    lines = []
    for reinf in slope_data["reinforcement_lines"]:
        lines.append([(reinf["x1"], reinf["y1"]), (reinf["x2"], reinf["y2"])])
    
    mesh = build_mesh_from_polygons(polygons, target_size=1.5, element_type='tri3', lines=lines)
    
    print(f"Mesh created: {len(mesh['nodes'])} nodes, {len(mesh['elements'])} 2D elements")
    if 'elements_1d' in mesh:
        print(f"Reinforcement: {len(mesh['elements_1d'])} 1D elements")
    
    # Build FEM data
    fem_data = build_fem_data(slope_data, mesh)
    
    # Solve FEM
    print("\nSolving FEM with reinforcement...")
    solution = solve_fem(fem_data, F=1.0, debug_level=1)
    
    if solution["converged"]:
        print(f"✓ FEM solution converged in {solution['iterations']} iterations")
        
        # Reinforcement force analysis
        forces_1d = solution.get("forces_1d", [])
        if len(forces_1d) > 0:
            max_force = np.max(np.abs(forces_1d))
            print(f"Maximum reinforcement force: {max_force:.1f} N")
        
        # Plot results
        fig, ax = plot_fem_results(fem_data, solution, plot_type='stress')
        plt.savefig('fem_test_reinforced_stress.png', dpi=150, bbox_inches='tight')
        print("Stress plot saved as 'fem_test_reinforced_stress.png'")
        
        # Plot reinforcement force profiles
        if len(forces_1d) > 0:
            fig, axes = plot_reinforcement_force_profiles(fem_data, solution)
            plt.savefig('fem_test_reinforcement_forces.png', dpi=150, bbox_inches='tight')
            print("Reinforcement force plot saved as 'fem_test_reinforcement_forces.png'")
        
        plt.show()
        
        return True
    else:
        print("✗ FEM solution failed to converge")
        return False


def test_ssrm_analysis():
    """
    Test SSRM analysis for factor of safety calculation.
    """
    print("\n=== Testing SSRM Analysis ===")
    
    # Create simple slope (potentially unstable)
    slope_data = create_simple_slope_data()
    
    # Make slope more unstable for testing
    slope_data["materials"][0]["c"] = 10000  # Reduce cohesion to 10 kPa
    slope_data["materials"][0]["phi"] = 15   # Reduce friction angle to 15°
    
    # Generate mesh
    polygons = [slope_data["profile_lines"][0]]
    mesh = build_mesh_from_polygons(polygons, target_size=2.0, element_type='tri3')
    
    # Build FEM data
    fem_data = build_fem_data(slope_data, mesh)
    
    # Run SSRM analysis
    print("\nRunning SSRM analysis...")
    ssrm_solution = solve_ssrm(fem_data, debug_level=1)
    
    if ssrm_solution["converged"]:
        FS = ssrm_solution["FS"]
        print(f"✓ SSRM analysis completed successfully")
        print(f"Factor of Safety: {FS:.3f}")
        print(f"SSRM iterations: {ssrm_solution['iterations_ssrm']}")
        
        # Plot convergence history
        fig, axes = plot_ssrm_convergence(ssrm_solution)
        plt.savefig('fem_test_ssrm_convergence.png', dpi=150, bbox_inches='tight')
        print("SSRM convergence plot saved as 'fem_test_ssrm_convergence.png'")
        
        # Plot final solution at critical F
        if "last_solution" in ssrm_solution:
            fig, ax = plot_fem_results(fem_data, ssrm_solution["last_solution"], plot_type='stress')
            ax.set_title(f'Final Solution at Critical FS = {FS:.3f}')
            plt.savefig('fem_test_ssrm_final.png', dpi=150, bbox_inches='tight')
            print("Final SSRM solution plot saved as 'fem_test_ssrm_final.png'")
        
        plt.show()
        
        return True
    else:
        print("✗ SSRM analysis failed")
        print(f"Error: {ssrm_solution.get('error', 'Unknown error')}")
        return False


def test_with_input_file(input_file):
    """
    Test FEM analysis using an existing XSLOPE input file.
    """
    print(f"\n=== Testing with Input File: {input_file} ===")
    
    if not os.path.exists(input_file):
        print(f"Input file not found: {input_file}")
        return False
    
    try:
        # Load slope data
        slope_data = load_slope_data(input_file)
        print("✓ Input file loaded successfully")
        
        # Check if mesh exists
        if 'mesh' not in slope_data or slope_data['mesh'] is None:
            print("No mesh found in slope_data, generating mesh...")
            # Extract polygons from profile lines  
            polygons = slope_data.get('profile_lines', [])
            if not polygons:
                print("No profile lines found for mesh generation")
                return False
            
            # Generate mesh
            mesh = build_mesh_from_polygons(polygons, target_size=2.0, element_type='tri3')
            slope_data['mesh'] = mesh
            print(f"✓ Mesh generated: {len(mesh['nodes'])} nodes, {len(mesh['elements'])} elements")
        
        # Build FEM data
        fem_data = build_fem_data(slope_data)
        print("✓ FEM data structure built")
        
        # Run SSRM analysis
        print("\nRunning SSRM analysis on input file...")
        ssrm_solution = solve_ssrm(fem_data, debug_level=1)
        
        if ssrm_solution["converged"]:
            FS = ssrm_solution["FS"]
            print(f"✓ SSRM analysis completed")
            print(f"Factor of Safety: {FS:.3f}")
            
            # Save results
            base_name = os.path.splitext(os.path.basename(input_file))[0]
            
            # Plot final solution
            if "last_solution" in ssrm_solution:
                fig, ax = plot_fem_results(fem_data, ssrm_solution["last_solution"], plot_type='stress')
                ax.set_title(f'{base_name}: FS = {FS:.3f}')
                plt.savefig(f'fem_{base_name}_solution.png', dpi=150, bbox_inches='tight')
                print(f"Solution plot saved as 'fem_{base_name}_solution.png'")
            
            # Plot convergence
            fig, axes = plot_ssrm_convergence(ssrm_solution)
            plt.savefig(f'fem_{base_name}_convergence.png', dpi=150, bbox_inches='tight')
            print(f"Convergence plot saved as 'fem_{base_name}_convergence.png'")
            
            plt.show()
            return True
        else:
            print("✗ SSRM analysis failed")
            return False
            
    except Exception as e:
        print(f"✗ Error processing input file: {e}")
        return False


def main():
    """
    Main function to run FEM tests.
    """
    print("XSLOPE Finite Element Analysis Test Suite")
    print("==========================================")
    
    # Check if input file was provided
    if len(sys.argv) > 1:
        input_file = sys.argv[1]
        success = test_with_input_file(input_file)
    else:
        # Run test suite with synthetic examples
        print("Running synthetic test cases...")
        
        # Test 1: Basic FEM
        success1 = test_simple_fem_analysis()
        
        # Test 2: FEM with reinforcement
        success2 = test_fem_with_reinforcement()
        
        # Test 3: SSRM analysis
        success3 = test_ssrm_analysis()
        
        success = success1 and success2 and success3
        
        if success:
            print("\n" + "="*50)
            print("✓ All tests completed successfully!")
            print("Check the generated PNG files for visualization results.")
        else:
            print("\n" + "="*50)  
            print("✗ Some tests failed. Check the error messages above.")
    
    return 0 if success else 1


if __name__ == "__main__":
    exit(main())