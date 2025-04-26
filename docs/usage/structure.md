# Code Structure

The slopetools python package is designed to be modular and easy to use. The following sections outline the main components of the package and how they can be used to perform slope stability analysis. Each of these components can be used independently or in combination to achieve the desired results. The following sections provide a brief overview of each component, along with example code snippets to illustrate their usage. The features are described in more detail in the other sections of the documentation.

## Loading the Input File

The first step in using slopetools is to load the input file. The input file contains all the necessary information about the slope, including the geometry, material properties, and loading conditions. The input file is an Excel file with a specific structure. A template for the Excel file can be downloaded here:

[input_template.xlsx](../input_template.xlsx)

The input file is designed to be easy to use and can be modified using any spreadsheet software. The input file is divided into several tabs or sheets, each of which corresponds to a specific aspect of the slope stability analysis. The strucutre of each sheet is designed to be intuitive and easy to understand. The sheets are as follows:

<div class="wrapped-table">
  <table>
    <thead>
      <tr>
        <th>Sheet Name</th>
        <th>Description</th>
      </tr>
    </thead>
    <tbody>
      <tr>
        <td>main</td>
        <td>A brief set of instructions and global variables including the unit weight of water, and tension crack parameters.</td>
      </tr>
      <tr>
        <td>plot</td>
        <td>A plot of the slope geometry based on the inputs in profile, piezo, and other sheets. This is intended to provide a quick visual check of the inputs.</td>
      </tr>
      <tr>
        <td>profile</td>
        <td>A set of tables for inputing the XY coordinates of up to 10 profile lines.</td>
      </tr>
      <tr>
        <td>mat</td>
        <td>A set of tables for inputting the material properties of up to 10 materials. This includes unit weight and shear strength properties.</td>
      </tr>
      <tr>
        <td>piezo</td>
        <td>A table for inputting the XY coordinates of piezometric line used to calculate pore pressures.</td>
      </tr>
      <tr>
        <td>circles</td>
        <td>A table for inputting the geometry of up to ten candidate circular failure surfaces.</td>
      </tr>
      <tr>
        <td>non-circ</td>
        <td>A table for inputting the XY coordinates of a non-circular slip surface.</td>
      </tr>
      <tr>
        <td>dloads</td>
        <td>A set of tables for inputting up to 8 distributed loads.</td>
      </tr>
      <tr>
        <td>reinforce</td>
        <td>A set of tables for inputting up to 12 soil reinforcement lines.</td>
      </tr>
    </tbody>
  </table>
</div>

The table can be loaded into Python using the 'load_globals' function. This function loads the input file and 
returns a dictionary containing the data from each sheet. The data can then be accessed using the sheet name as the 
key. For example:

```python
from slopetools import load_globals

data = load_globals("input_template.xlsx")
profile_lines = data["profile_lines"]
materials = data["materials"]
piezo_line = data["piezo_line"]
gamma_w = data["gamma_water"]
circle = data["circles"][0]  # or whichever one you want
non_circ = data["non_circ"]
dloads = data["dloads"]
max_depth = data["max_depth"]
reinforce_lines = data["reinforce_lines"]
```

The format of the data in each category is as follows: 

```python
#: [lb/ft^3] Unit weight of water
gamma_water = 62.4

#: [ft] Depth of the crack
tcrack_depth = 0.0

#: [ft] Water level in the crack
tcrack_water = 0.0

#: List of profile lines, each a list of (x, y) tuples.
profile_lines = []
# Example: profile_lines = [[(x1, y1), (x2, y2)], [(x3, y3), (x4, y4)]]

#: List of material property dictionaries.
materials = []
# Example: materials = [{'gamma': 120, 'c': 30, 'phi': 20, 'piezo': 0.5, 'sigma_gamma': 0.1, 'sigma_c': 0.1, 'sigma_phi': 0.1}]

#: List of (x, y) coordinates defining the piezometric line.
piezo_line = []
# Example: piezo_line = [(x1, y1), (x2, y2)]

#: [ft] Maximum depth for circular failure surfaces.
max_depth = 0.0

#: List of circular failure surface dictionaries.
circles = []
# Example: circles = [{'Xo': 120, 'Yo': 80, 'Option': "Depth", 'Depth': -10, 'Xi': 5, 'Yi': 5}]

#: List of points in a non-circular surface with movement options.
non_circ = []
# Example: non_circ = [{'X': 120, 'Y': 80, 'Movement': "Free"}, {'X': 130, 'Y': 90, 'Movement': "Horiz"}]

#: List of distributed load dictionaries.
dloads = []
# Example: dloads = [{'X': 120, 'Y': 80, 'Normal': 100}, {'X': 130, 'Y': 90, 'Normal': 150}]

#: List of reinforcement line dictionaries.
reinforce_lines = []
# Example: reinforce_lines = [{'X': 120, 'Y': 80, 'FL': 10, 'FT': 20}, {'X': 130, 'Y': 90, 'FL': 15, 'FT': 25}]

```

## Building Slice Data

Once the data are loaded, the next step is generally to build the slice data using the circular or non-circular 
failure surface, profile lines, material properties, etc. Before building the slice data, you need to define the 
ground surface geometry. The ground surface is defined as the uppermost surface of the slope, defined by the profile 
line data. This is done using the `build_ground_surface` function in the `utils` module. This is how you call the function:

```python
from slice import build_ground_surface

ground_surface = build_ground_surface(profile_lines)
```

The slice data are then constructed using the `generate_slices` function in the `slice` module, as follows:

```python
from slice import generate_slices

df, failure_surface = generate_slices(
    profile_lines=profile_lines,
    materials=materials,
    ground_surface=ground_surface,
    circle=circle,
    num_slices=20,
    gamma_w=gamma_w,
    piezo_line=piezo_line,
    dloads=dloads,
)
```

The function returns a dataframe containing the slice data and the failure surface object. The slice data is a 
pandas dataframe including columns for properties associated with each slice, such as the slice number, width, 
weight, etc. The failure surface object contains the geometry of the failure surface, clipped such that includes 
only the portion that is inside the slope.

## Solving for Factor of Safety

Once the slice data are generated, the next step is to solve for the factor of safety. This is done using one of the 
solution methods available in the `solve` module. The `solve` module includes several different methods for solving 
for the factor of safety, including the following:

| Method Name | Description                                      |
| ----------- |--------------------------------------------------|
| oms       | Ordinary Method of Slices                        |
| bishop    | Bishop's Simplified Procedure                    |
| spencer  | Spencer's Method  (complete equilibrium)         |
| janbu     | Janbu's Method    (force equilibrium)            |
| morgenstern_price | Morgenstern-Price Method  (complete equilibrium) |

Each of these methods is described in more detail in the [Solution Techniques](../methods) section of this documentation.

Here is an example of how to call the bishop method:

```python
FS, N, converge = bishop(df)
if not converge:
    print("Bishop's method did not converge.")
```

The `bishop` function takes the slice data dataframe as input and returns the factor of safety, the normal forces on 
the bottom each slice, and a convergence flag. The convergence flag indicates whether the method converged to a solution. 

## Automated Search

--- UNDER CONSTRUCTION ---

## Plotting

After finding a solution, the results can be plotted using the `plot` module. The `plot` module includes several 
functions for rendering the slope geometry, failure surface, and slice data using matplotlib. Here is an example of 
how to call the `plot_slices` function:

```python
from plot import plot_slices

plot_slices(profile_lines, df, piezo_line=piezo_line, failure_surface=failure_surface, fs=FS, dloads=dloads, max_depth=max_depth)
```


The `plot` module includes several functions for plotting the results of the slope stability analysis. The

## Advanced Tools
