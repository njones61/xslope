# XSLOPE Input Template

## Overview

The XSLOPE input template is the primary means of defining slope stability problems in xslope. It is an Excel workbook that contains all necessary information about the slope geometry, material properties, loading conditions, boundary conditions, and analysis parameters. The template is designed to support three main types of analysis:

- **Limit Equilibrium Method (LEM)**: Classical slope stability analysis using methods like Bishop, Spencer, Janbu, and others
- **Seepage Analysis**: Finite element groundwater flow analysis for steady-state and transient conditions
- **Finite Element Method (FEM)**: Stress-deformation analysis including the Shear Strength Reduction Method (SSRM)

The template uses a structured format with multiple worksheets (tabs), each dedicated to a specific aspect of the problem definition. This organization makes it easy to prepare complex slope stability analyses while maintaining clarity and avoiding errors.

## Download

A template for the Excel file can be downloaded here:

[input_template.xlsx](../inputs/input_template.xlsx)

The input file can be modified using any spreadsheet software such as Microsoft Excel, LibreOffice Calc, or Google Sheets.

## Loading the Template

The template is loaded into xslope using the `load_slope_data()` function from the `fileio` module:

```python
from xslope.fileio import load_slope_data

# Load the input template
slope_data = load_slope_data("path/to/your/input_template.xlsx")
```

This function reads all worksheets, validates the data, and returns a dictionary containing all parsed information. The resulting `slope_data` dictionary is then used by various xslope modules for analysis.

## Template Structure

The template consists of 10 worksheets, each serving a specific purpose. Different worksheets are used by different analysis types: Limit Equilibrium Method (LEM), seepage analysis (SEEP), and Finite Element Method (FEM).

| Sheet Name | Description | LEM | SEEP | FEM |
|------------|-------------|:---:|:----:|:---:|
| **main** | Global parameters and instructions | X   | X    | X   |
| **plot** | Auto-generated geometry preview | X   | X    | X   |
| **profile** | XY coordinates of profile lines defining slope geometry | X   | X    | X   |
| **mat** | Material properties including strength, permeability, and stiffness | X   | X    | X   |
| **piezo** | Piezometric lines for pore pressure calculations | X   |      | X   |
| **circles** | Circular failure surface definitions | X   |      |     |
| **non-circ** | Non-circular failure surface coordinates | X   |      |     |
| **dloads** | Distributed surface loads | X   |      | X   |
| **reinforce** | Soil reinforcement elements (anchors, nails, geosynthetics) | X   |      | X   |
| **seep bc** | Seepage analysis boundary conditions |     | X    |     |

The following sections describe each tab in detail, including the data structure and how it is used in analysis.

---

## Worksheet: main

![main tab screenshot](screenshots/main.png)

The **main** worksheet provides global parameters that apply to all analyses and serves as the instruction page for the template. This tab contains:

- **Unit weight of water** (γw): Used in pore pressure calculations and buoyancy effects
- **Tension crack parameters**: Depth and water level within tension cracks at the top of the failure surface
- **Seismic coefficient** (kh): Horizontal seismic acceleration coefficient for pseudo-static earthquake analysis
- **Template version**: Tracks template format for compatibility

These global parameters are accessed throughout the analysis. For example, the unit weight of water is used in computing pore pressures from piezometric lines, and the seismic coefficient is used to add horizontal inertial forces to each slice in limit equilibrium calculations. The tension crack parameters allow modeling of the common condition where tensile stresses at the crest of a slope cause cracking, which may fill with water and reduce stability.

---

## Worksheet: plot

![plot tab screenshot](screenshots/plot.png)

The **plot** worksheet contains an auto-generated visual preview of the slope geometry based on the inputs in other tabs. This plot updates automatically when you modify the **profile**, **piezo**, **circles**, or other geometry-related worksheets (note: auto-update requires Excel recalculation).

This worksheet is purely for quality control and visualization within the Excel environment. It allows you to quickly verify that:

- Profile lines are correctly positioned and form a reasonable slope geometry
- Piezometric lines are within the slope boundaries
- Circular failure surfaces intersect the ground surface at reasonable locations

The plot is not used by xslope during analysis - it exists solely to help you validate your inputs before running calculations. When working with complex multi-layer geometries or multiple water tables, this visual check can catch data entry errors early.

---

## Worksheet: profile

![profile tab screenshot](screenshots/profile.png)

The **profile** worksheet defines the slope geometry using XY coordinates of profile lines. Each profile line represents the interface between two soil layers or a ground surface. The template accommodates up to 10 profile lines, organized in three blocks of tables.

Each table has paired X and Y columns where you enter coordinates sequentially from left to right (for horizontal layers) or in any consistent direction (for irregular boundaries). At least two points are required per line, but you can use up to 15 points per profile line to capture complex geometries.

The topmost profile line defines the ground surface. Additional profile lines below define layer boundaries. During analysis, xslope uses these lines to:

1. Construct the ground surface by finding the highest elevation at each x-coordinate
2. Determine slice geometry when a failure surface is intersected with the profile
3. Assign material properties to slices based on which layer they fall within
4. Build polygons for finite element meshing in seepage or FEM analysis

Profile lines must be continuous and cover the full width of the analysis domain. The order matters: the first material in the **mat** worksheet corresponds to the first profile line, the second material to the second profile line, and so on.

---

## Worksheet: mat

![mat tab screenshot](screenshots/mat.png)

The **mat** worksheet defines material properties for each soil layer. Each row corresponds to one profile line from the **profile** worksheet, establishing a one-to-one correspondence between geometry and properties. The template supports up to 10 materials with comprehensive property definitions for strength, permeability, and stiffness.

**Strength Properties** (for LEM analysis):

- **c** (cohesion) and **φ** (friction angle): Mohr-Coulomb shear strength parameters
- **Option**: Specifies how strength varies with depth (constant, linearly increasing, etc.)
- **cp** and **r-elev**: Peak cohesion and reference elevation for variable strength profiles

**Pore Pressure Options** (column U):

- **piezo**: Use piezometric line from **piezo** worksheet
- **seep**: Interpolate from seepage analysis solution (requires mesh and solution files in cells L22:L24)
- **none**: No pore pressure
- **ru**: Pore pressure ratio

**Variability** (for reliability analysis):

- **σ(γ)**, **σ(c)**, **σ(φ)**, etc.: Standard deviations for probabilistic analysis

**Permeability** (for seepage analysis):

- **k1**, **k2**: Major and minor hydraulic conductivity (can be anisotropic)
- **alpha**: Orientation angle of permeability tensor
- **kr0**, **h0**: Unsaturated flow parameters (relative conductivity and suction head)

**Stiffness** (for FEM analysis):

- **E**: Young's modulus
- **ν**: Poisson's ratio

When using the "seep" pore pressure option, you must specify in cells L22:L24:

- L22: Mesh filename (JSON format from `export_mesh_to_json()`)
- L23: Solution filename for condition 1 (CSV format from `export_seep_solution()`)
- L24: Solution filename for condition 2 (optional, for rapid drawdown analysis)

---

## Worksheet: piezo

![piezo tab screenshot](screenshots/piezo.png)

The **piezo** worksheet defines piezometric lines for calculating pore water pressures in limit equilibrium analysis. A piezometric line (also called a phreatic surface or water table) represents the elevation at which pore pressure equals atmospheric pressure. Below this line, pore pressures are positive (hydrostatic); above it, they are zero (or negative if suction is considered).

The tab provides space for two piezometric lines (columns A-B and D-E), which is useful for rapid drawdown analysis:

- **First line (A-B)**: Steady-state or initial condition water table
- **Second line (D-E)**: Drawdown condition water table (optional)

Each piezometric line requires at least two XY coordinate pairs and can accommodate up to 14 points. The coordinates define the line, and xslope interpolates pore pressures at each slice base using the vertical distance from the slice base to the piezometric line (u = γw × h, where h is the vertical distance).

For complex problems with perched water tables or confined aquifers, piezometric lines may not provide sufficient flexibility. In those cases, use the "seep" option in the **mat** worksheet to couple with a full finite element seepage analysis.

---

## Worksheet: circles

![circles tab screenshot](screenshots/circles.png)

The **circles** worksheet defines circular failure surfaces for limit equilibrium analysis. Circular surfaces are the most common assumption in slope stability analysis and are required for methods like Bishop's Simplified Method.

Each row specifies one circular failure surface with the following parameters:

- **Xo, Yo**: Coordinates of the circle center
- **Option**: Method for defining circle size - "Depth", "Radius", or "Intercept"
>>- **Depth**: Specify depth below ground surface at center location
>>- **Radius**: Directly specify circle radius
>>- **Intercept**: Specify a point (Xi, Yi) where circle should pass through
- **Depth, R, Xi, Yi**: Associated values depending on option selected

The template supports up to 10 circles, which can be analyzed individually or used as starting points for automated search algorithms. During analysis, xslope:

1. Constructs the circular arc geometry from the parameters
2. Finds intersection points with the ground surface
3. Divides the arc into slices
4. Assigns material properties to each slice based on its position

For simple problems, manually defining a few trial circles is often sufficient. For complex geometries or finding the critical failure surface, use the `circular_search()` function which employs adaptive grid refinement to systematically explore the space of potential failure surfaces.

---

## Worksheet: non-circ

![non-circ tab screenshot](screenshots/non-circ.png)

The **non-circ** worksheet allows definition of arbitrary non-circular failure surfaces. This is essential for problems where circular surfaces are inappropriate, such as:

- Slopes with weak layers or interfaces that control failure geometry
- Wedge failures in rock slopes
- Composite surfaces combining circular and planar segments
- Surfaces identified from FEM analysis or optimization algorithms

The worksheet contains three columns:

- **X, Y**: Coordinates of points along the failure surface, defined sequentially from entry to exit
- **Movement**: Direction of movement constraint at each point (used in advanced analyses)

Non-circular surfaces require methods like Spencer's method, Morgenstern-Price, or the General Limit Equilibrium (GLE) method. Bishop's method, which assumes circular geometry, cannot be used with non-circular surfaces.

When using non-circular surfaces, ensure that:

- The surface starts and ends outside the slope (intersecting the ground surface or boundaries)
- Points are ordered consistently (left to right or right to left)
- The surface is continuous (no gaps between points)

You can generate non-circular surfaces manually, from field observations, or by exporting results from automated search algorithms like `noncircular_search()`.

---

## Worksheet: dloads

![dloads tab screenshot](screenshots/dloads.png)

The **dloads** worksheet defines distributed surface loads applied to the slope. These represent surcharge loads such as traffic, buildings, stockpiled materials, or other surface loading. The worksheet supports up to 8 distributed loads, organized in two sets of 4 blocks each.

Each distributed load is defined by a series of points with:

- **X, Y**: Coordinates of points along the load distribution line
- **Normal**: Normal stress (force per unit area) acting perpendicular to the line

The loads are applied as additional normal stresses on slices that intersect the load distribution line. At least two points are required to define each load block. The load distribution can follow the ground surface or can be defined at any orientation.

**Two load sets** are provided (rows 4-13 and rows 17-26) to support rapid drawdown or staged loading analysis:

- **First set**: Initial or steady-state loading condition
- **Second set**: Modified loading condition (for drawdown or staged construction)

During limit equilibrium analysis, distributed loads increase the normal force on the slice base, which affects both the driving and resisting forces depending on the slope angle and load orientation. Vertical loads on horizontal surfaces increase resistance, while loads on steep slopes may increase driving forces.

---

## Worksheet: reinforce

![reinforce tab screenshot](screenshots/reinforce.png)

The **reinforce** worksheet defines soil reinforcement elements such as soil nails, rock anchors, geosynthetic reinforcement, or tiebacks. These elements provide additional resistance to sliding by mobilizing tensile forces along the failure surface. The template accommodates up to 20 reinforcement lines (rows 3-22).

Each reinforcement element is defined by:

- **Geometry**:<br>
>>x1, y1: Start point coordinates <br>
>>x2, y2: End point coordinates<br>
- **Strength Properties**:<br>
>>Tmax: Maximum tensile force that can be mobilized<br>
>>Tres: Residual tensile force (for post-peak behavior)<br>
- **Bond Properties**:<br>
>>Lp1: Pullout bond length at start end<br>
>>Lp2: Pullout bond length at end end<br>
- **Stiffness** (for FEM analysis):<br>
>>E: Elastic modulus of reinforcement<br>
>>Area: Cross-sectional area<br>

The pullout lengths (Lp1, Lp2) control how tensile force develops along the reinforcement:

- If Lp = 0, the end is fully anchored (immediate maximum tension)
- If Lp > 0, tension develops linearly over the pullout length
- If the total line length < Lp1 + Lp2, only partial tension is mobilized

During limit equilibrium analysis, xslope:

1. Identifies reinforcement lines that intersect the failure surface
2. Calculates the mobilized tensile force based on pullout mechanics
3. Resolves the force into components tangent and normal to the failure surface
4. Adds the resisting component to the factor of safety calculation

Reinforcement is particularly important for analyzing stabilized slopes, MSE walls, and anchored retaining structures.

---

## Worksheet: seep bc

![seep bc tab screenshot](screenshots/seep bc.png)

The **seep bc** worksheet defines boundary conditions for finite element seepage analysis. Boundary conditions specify where water enters or exits the domain and control the groundwater flow solution. The worksheet provides two complete sets of boundary conditions to support rapid drawdown analysis.

**Boundary Condition Types:**

1. **Specified Head Boundaries** (columns B-C, E-F, H-I):

>   - Define locations where hydraulic head is known (e.g., reservoir levels, constant head boundaries)
>   - Each specified head region requires:
>>  - Head value (in the header row)
>>  - XY coordinates of points along the boundary (minimum 2 points)
>   - Up to three specified head regions can be defined per set

2. **Exit Face Boundary** (columns B-C, rows 16-23):

   - Defines where water exits the domain (e.g., downstream toe, seepage face)
   - Special boundary condition where head = elevation (zero pressure)
   - Requires at least 2 XY coordinate pairs

**Two Boundary Condition Sets:**

- **Set 1** (rows 3-23): Initial or steady-state condition (e.g., full reservoir)
- **Set 2** (rows 26-46): Modified condition (e.g., rapid drawdown, reservoir empty)

During seepage analysis, xslope:

1. Builds a finite element mesh from the **profile** geometry
2. Applies specified head values at nodes on the specified head boundaries
3. Applies exit face conditions where water exits
4. Solves the Laplace equation (∇·(k∇h) = 0) for hydraulic head at all nodes
5. Computes pore pressures (u = γw(h - y)) and flow velocities

For coupled seepage-LEM analysis, the computed pore pressures can be exported and later interpolated onto slice bases using the "seep" option in the **mat** worksheet.

---

## Notes

- All coordinates should use consistent units throughout the template (typically feet or meters)
- Angles are specified in degrees
- The template is designed to be flexible - you need not fill in all worksheets for every analysis
- For simple LEM analyses, you only need: main, profile, mat, piezo, and circles (or non-circ)
- For seepage analyses, you only need: main, profile, mat, and seep bc
- Always check the **plot** worksheet after entering data to visually verify your geometry
- Templates can be saved and reused for parametric studies or similar projects