# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

xslope is a Python package for limit equilibrium slope stability analysis. It provides multiple solution methods for slope stability calculations including:

- Ordinary Method of Slices (OMS)
- Bishop's Simplified Method
- Janbu Method
- Spencer's Method
- Corps of Engineers Method
- Lowe & Karafiath Method
- Rapid Drawdown Analysis

The package also includes seepage analysis capabilities and integrates with finite element mesh analysis.

## Key Commands

### Running the Main Analysis
```bash
python main.py
```
Solves one circle on a shipped sample with Spencer, compares every LEM method
on the same slices, then searches for the critical circle.

### Documentation
```bash
# Build documentation (requires mkdocs)
mkdocs build
mkdocs serve  # For local development server
```

### Driver Scripts
Each is a plain demonstration script pointed at a shipped input file; edit the
path at the top to run your own model. All of them open matplotlib windows, so
run them from a desktop session (or set `MPLBACKEND=Agg` to skip the figures).
```bash
python main_lem.py            # LEM: inputs plot, auto search, critical surface
python main_search.py         # Circular search on its own, printing each iteration
python main_reliability.py    # Taylor-series reliability from the mat-sheet std devs
python main_fem.py            # FEM / SSRM run and results plots
python main_seep.py           # Seepage: mesh, solve, flow net, export the .csv
python main_seep2d.py         # Seepage from a legacy SEEP2D .s2d file
python main_design.py         # Slope-angle sweep to a target FS or reliability index
python main_mesh.py           # Mesh a reinforced model and check 1D element alignment
python main_fileio_test.py    # Load a workbook and dump every slope_data key
python main_pickle.py         # Cache a parsed model to a pickle and read it back
```

## Code Architecture

### Core Modules

- **`fileio.py`**: Handles Excel input file parsing and data loading. Contains `load_slope_data()` function for reading input templates.
- **`slice.py`**: Generates analysis slices from slope geometry. Main function `generate_slices()` creates the slice dataframe used by all solution methods.
- **`solve.py`**: Contains all limit equilibrium solution methods (oms, bishop, janbu, spencer, etc.). Each method returns `(success, result)` tuple.
- **`search.py`**: Implements automated search algorithms for finding critical failure surfaces, including `circular_search()` with adaptive grid refinement.
- **`plot.py`**: Visualization functions for slopes, slices, and results. Contains `plot_solution()` and `plot_inputs()`.
- **`seep.py`**: Seepage analysis integration with finite element mesh data.
- **`mesh.py`**: Finite element mesh handling and interpolation functions.

### Data Flow

1. Input data is loaded from Excel templates via `load_slope_data()`
2. Slice geometry is generated using `generate_slices()` 
3. Solution methods in `solve.py` process the slice dataframe
4. Results are visualized using functions in `plot.py`
5. Search algorithms in `search.py` can automate finding critical surfaces

### Input Files

- Editable master template: `docs/inputs/input_template.xlsx` (version in `main!D5`); bundled wheel copy `xslope/resources/input_template.xlsx` must stay byte-identical (checked by `run_tests.py`)
- Sample/benchmark input files: `docs/inputs/slope/`
- Companion files are auto-discovered by naming convention next to the .xlsx: `{base}_mesh.json`, `{base}_seep.csv`, `{base}_seep2.csv` (rapid drawdown)

### Key Data Structures

- **slice_df**: Main pandas DataFrame containing slice properties (alpha, phi, c, w, u, dl, etc.)
- **slope_data**: Dictionary containing all input parameters including ground surface, materials, circles, etc.
- **failure_surface**: Geometry object representing the failure surface
- **results**: Dictionary containing FS (factor of safety) and method-specific parameters

### Solution Method Integration

All solution methods follow the same pattern:
```python
success, result = method_name(slice_df)
if success:
    print(f"Method: FS={result['FS']:.3f}")
```

The `solve_all()` function in `main.py` demonstrates running all methods sequentially.

### Seepage Integration

Seepage analysis uses finite element mesh data to determine pore pressures:
- Mesh imported from JSON files using `import_mesh_from_json()`
- Pore pressures interpolated at slice locations
- Integrated into limit equilibrium analysis through the 'u' column in slice_df

## Development Notes

- The package uses shapely for geometric operations
- Pandas DataFrames are central to the slice-based analysis approach
- All angles are handled in degrees internally
- Error handling follows (success, result) pattern throughout
- Visualization uses matplotlib with custom styling