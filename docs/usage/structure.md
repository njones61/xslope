# Code Structure

The **slopetools** python package is designed to be modular and easy to use. The following sections outline the main components of the package and how they can be used to perform slope stability analysis. Each of these components can be used independently or in combination to achieve the desired results. The following sections provide a brief overview of each component, along with example code snippets to illustrate their usage. The features are described in more detail in the other sections of the documentation.

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

success, result = generate_slices(data, circle, num_slices=20)
if success:
    df, failure_surface = result
else:
    print(result)
```

The function returns a dataframe containing the slice data and the failure surface object. The slice data is a 
pandas dataframe including columns for properties associated with each slice, such as the slice number, width, 
weight, etc. The failure surface object contains the geometry of the failure surface, clipped such that includes 
only the portion that is inside the slope.

## Solving for Factor of Safety

Once the slice data are generated, the next step is to solve for the factor of safety. This is done using one of the 
solution methods available in the `solve` module. The `solve` module includes several different methods for solving 
for the factor of safety, including the following:

| Method Name     | Description                                     |
|-----------------|-------------------------------------------------|
| oms             | Ordinary Method of Slices                       |
| bishop          | Bishop's Simplified Procedure                   |
| spencer         | Spencer's Method  (complete equilibrium)        |
| janbu_corrected | Janbu's Corrected Method    (force equilibrium) |


Each of these methods is described in more detail in the [Solution Techniques](../methods) section of this documentation.

Here is an example of how to call the bishop method:

```python
success, result = spencer(df, circular=True)
if not success:
    print(f'Error: {result}')
else:
    print(f'Spencer: FS={result["FS"]:.3f}, theta={result["theta"]:.2f} degrees')
```

Each solver function has the same set of inputs and outputs. The inputs are the slice df and a flag to indicate whether is it is circular surface or non-circular. Each function returns two arguments. The first is a success flag. It success=False, the results argument contains a message describing how the method failed. If success=True, results is a dictionary containing the factor of safety, and values specific to the method. For example, with Spencer's method, the dictionary contains the side force inclincation, theta as shown in the example above. 

## Automated Search

--- UNDER CONSTRUCTION ---

## Plotting

After finding a solution, the results can be plotted using the `plot` module. The `plot` module includes several 
functions for rendering the slope geometry, failure surface, and slice data using matplotlib. Here is an example of 
how to call the `plot_inputs` function after loading the input file. 

```python
from fileio import load_globals
from plot import plot_inputs

data = load_globals("docs/input_template.xlsx")

plot_inputs(data)
```

To plot the solution for a specific method, you can use the `plot_solution` function. Here is an example:

```python
from plot import plot_solution

success, result = spencer(df, circular=True)
if not success:
    print(f'Error: {result}')
else:
    plot_solution(data, df, failure_surface, result)
```

The `plot` module includes several functions for plotting the results of the slope stability analysis. The

## Advanced Tools
