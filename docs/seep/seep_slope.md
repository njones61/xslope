# Seepage - Slope Stability Integration

The seepage analysis tools in XSLOPE can be used independently to simulation flow through porous media for a broad variety of applications. However, the seepage analysis tools can also be used in conjunction with the slope stability analysis tools. Since both seepage and slope stability analysis use the same input file and the same approach to defining material properties and model geometry, integrating the two analysis tools is straightforward. The same Excel input template is used for both analyses. The seepage analysis is performed first and the seepage solution, inlcuding pore pressures at nodes is saved to a set of files. Next, the slope stability analysis is performed and the nodal pore pressures are used in the slope stability analysis to calculate effective stresses for selected soils. Each step of this process is described in detail below.

## Saving the Seepage Solution

After performing the seepage analysis, the seepage solution should be saved to a set of two files. The first is a JSON file containing the 2D finite element mesh including both nodal coordinates and element topology. The second file is a CSV file containing the seepage solution, including hydraulic heads, pressure head, pore pressures, and  total head at each node. The seepage solution files should be saved using a specific naming convention. If the input file is named `myfile.xlsx`, the seepage files should be exported as `myfile_mesh.json` and `myfile_seep.csv`. In other words, the prefix of the input file should be used as the prefix for the seepage files and they should end with the suffix `_mesh.json` and `_seep.csv`.  The seepage solution files are saved in the same folder as the Excel input template. The following code snippet shows how to save the seepage solution files:

```python
from xslope.mesh import export_mesh_to_json
from xslope.seep import export_seep_solution

path = "docs/inputs/slope/xslope_lface.xlsx"
path_prefix = path.replace(".xlsx", "")  
slope_data = load_slope_data(path)

# CODE HERE TO PERFORM SEEPAGE ANALYSIS

export_mesh_to_json(mesh, path_prefix + "_mesh.json")
export_seep_solution(seep_data, solution, path_prefix + "_seep.csv")
```

## Rapid Drawdown Analysis

The input template is set up to allow the specification of two sets of boundary conditions for the seepage analysis. 
These are stored in the **seep bc** sheet and the **seep bc (2)** sheet. The first set is used for a typical seepage 
analysis, but when conducting a rapid drawdown analysis, the first set of boundary conditions can be used to simulate 
the full reservoir condition and the second set of boundary conditions can be used to simulate the drawn-down 
condition. Both of the solution can then be imported to define pore pressures for both the full reservoir and 
drawn-down conditions when performing slope stability analysis. When performing the seepage analysis, both solutions 
should be found using the same mesh geometry. Then a single mesh file is exported in JSON format and the two 
solutions are exported in CSV format. The second seepage solution should be save with a `_seep2.csv` extension. The 
following code snippet shows how to export the seepage solution files for the case with rapid drawdown:

```python
from xslope.mesh import export_mesh_to_json
from xslope.seep import export_seep_solution

path = "docs/inputs/slope/xslope_lface.xlsx"
path_prefix = path.replace(".xlsx", "")  
slope_data = load_slope_data(path)

# CODE HERE TO PERFORM SEEPAGE ANALYSIS

export_mesh_to_json(mesh, path_prefix + "_mesh.json")
export_seep_solution(seep_data, solution, path_prefix + "_seep.csv")
export_seep_solution(seep_data, solution2, path_prefix + "_seep2.csv")
```

## Selecting the SEEP Option for Pore Pressure Calculation

When performing a slope stability analysis, a set of shear strength parameters and a pore pressure option are 
entered for each material. The pore pressure options are as follows:

| Option | Description |
|------- |-------------|
| none   | No pore pressures. This is appropriate for dry soils or total stress analysis |
| piezo  | Pore presssure are derived from a piezometric line defined on the `piezo` tab |
| seep   | Pore pressures are derived from the results of a 2D finite element seepage analysis |

 For the `seep` option, the user first performs a seepage analysis as described above and exports the resulting finite 
 element mesh and solution file using the specified naming convention. For a rapid drawdown analysis, one mesh file 
 and two solution files are exported.

When XSLOPE reads an Excel input file, it automatically checks to see if a set of seepage files exists for the 
simulation. It does this by checking if the seepage files exist in the same location and with 
the same prefix as the Excel input file and 
ending with the suffix `_mesh.json`, `_seep.csv`, and `_seep2.csv`. If the seepage files exist, XSLOPE imports them 
and stores them in the slope_data dictionary. These pore pressures are then used to calculate pore pressures at the 
bottom of each slice in the slope stability analysis. 

### Pore Pressure Calculation

The primary output of seepage analysis is the hydraulic head distribution throughout the analysis domain. However, slope stability analysis requires pore pressure values rather than hydraulic heads. The conversion between these quantities follows the fundamental relationship:

>>$u = \gamma_w (h - z)$

where:
- $u$ is the pore water pressure
- $\gamma_w$ is the unit weight of water (typically 9.81 kN/m³ or 62.4 lbf/ft³)
- $h$ is the hydraulic head from seepage analysis
- $z$ is the elevation coordinate

This relationship automatically accounts for the effect of elevation on pore pressure while preserving the hydraulic gradients computed by the seepage analysis. Positive pore pressures indicate groundwater under pressure (below the phreatic surface), while zero or negative values indicate conditions above the groundwater table.

### Conservative Treatment of Negative Pore Pressures

In certain situations, the seepage analysis may compute negative pore pressures (tensions) in portions of the domain, particularly in unsaturated zones or regions with strong capillary effects. While these negative pore pressures may be physically realistic, their inclusion in slope stability analysis can lead to unconservative results by artificially increasing effective stresses and apparent slope stability.

XSLOPE implements a conservative approach by automatically setting any computed negative pore pressures to zero before transferring them to slope stability calculations. This treatment ensures that:

**Conservative Assessment:** Slope stability calculations do not benefit from potentially unreliable tensile strength in pore water
**Physical Realism:** Avoids dependence on soil-water tension effects that may not persist under loading conditions
**Robustness:** Eliminates sensitivity to uncertain unsaturated soil parameters that may not be well-characterized

### Interpolation to Slice Locations

The integration of seepage analysis results with limit equilibrium slope stability calculations requires interpolation of pore pressures from the seepage mesh nodes to the slice centroids used in limit equilibrium analysis. XSLOPE implements a robust interpolation scheme that preserves the accuracy of the computed pore pressure field while ensuring computational efficiency.

The interpolation process follows these key steps:

**Spatial Search:** For each slice centroid, the system identifies the seepage analysis element containing that point using efficient spatial search algorithms.

**Shape Function Evaluation:** The pore pressure at the slice centroid is computed using the finite element shape functions and the nodal pore pressure values of the containing element:

>>$u(x,y) = \sum_{i=1}^{n} N_i(x,y) \cdot u_i$

where $N_i$ are the element shape functions and $u_i$ are the nodal pore pressures.