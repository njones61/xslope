# Sample Problems - Seepage Analysis

The following examples illustrate how to use XSLOPE to perform seepage analysis. The problems feature both saturated and unsaturated conditions. Each of the Excel input files below can be used with the following notebook which has been set up specifically for running seepage analyses:

<a href="https://colab.research.google.com/github/njones61/xslope/blob/main/notebooks/xslope_seep.ipynb" target="_"><img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open In Colab"/></a>

These problems feature standalone seepage analyses. For instructions on to run an integrated seepage analysis with slope stability analysis, see the [Integrated Seeapage and Slope Stability Analysis](seep_slope.md) page.

### 1. Sheetpile with Clay Blanket

This is a saturated problem with a partially penetrating sheetpile and a clay blanket. It should have an upstream head BC = 13m up 
to the tip of the blanket and a downstream head BC = 10m. The profile line should follow the edge of the sheetpile (down 
and then back up) with a small gap to ensure that there is a crack in the resulting mesh. 

![clay_blanket_prob.png](images/clay_blanket.png)

The following Excel file illustrates how the inputs should be structured. Since this is a fully saturated problem, 
the kr0 and h0 material parameters are ignored. 

Excel input file: [xslope_clay_blanket.xlsx](files/xslope_clay_blanket.xlsx)

The solution should look something like this:

![clay_blanket_solution.png](images/clay_blanket_solution.png){width=1200px}

<!-- test: file=files/xslope_clay_blanket.xlsx, type=seep, expected_flowrate=40.062, tolerance=0.05 -->

### 2. Sea Trench

This is another saturated problem representing the excavation of a trench in a harbor supported by a parallel set of sheetpile walls. The sheetpiles pass through an upper silt layer down to a lower permeability silty clay layer.

![sea_trench.png](images/sea_trench.png)

The properties of the soil layers are as follows:

| Soil Layer | K1  | K2  |
|:----------:|:---:|:---:|
|    Silt    | 0.5 | 0.5 |
| Silty Clay | 0.1 | 0.1 |

Since this is a fully saturated problem, the kr0 and h0 material parameters are ignored. The problem set up requires 3 profile lines: 1 at the top of the silt layer on the left side, 1 at the top of the silt layer on the right, and 1 at the top of the silty clay layer that goes all the way from the left side to the right side of the problem. This profile line includes a small gap at the location of each sheetpile penetration to create a no-flow boundary along the edge of the sheetpile. The following Excel file illustrates how the inputs should be structured.

Excel input file: [xslope_sea_trench.xlsx](files/xslope_sea_trench.xlsx)

Solution:

![sea_trench_solution.png](images/sea_trench_solution.png){width=1000px}

<!-- test: file=files/xslope_sea_trench.xlsx, type=seep, expected_flowrate=4.347, tolerance=0.05 -->

### 3. Earth Dam with Core

The following diagram illustrates a simple earth dam with a clay core and a granular shell:

![earth_dam1.png](images/earth_dam1.png)

This problem requires an upstream head BC, a small downstream head BC, and a downstream exit face BC from the crest 
of the dam down to the tailwater. To build the input file, the following list of coordinates can be used:

![earth_dam1_pts.png](images/earth_dam1_pts.png)

In this case, the solution is partially saturated, so the kr0 and h0 parameters must be specfied for each material. 
The following Excel file contains a complete set of inputs for this problem:

[xslope_earth_dam1.xlsx](files/xslope_earth_dam1.xlsx)

The solution should look something like this:

![earth_dam1_solution.png](images/earth_dam1_solution.png){width=1200px}

<!-- test: file=files/xslope_earth_dam1.xlsx, type=seep, expected_flowrate=42.342, tolerance=0.05 -->

### 4. Johnson Reservoir

This is another earth dam problem with a shell, a core, and a foundation. 

![johnson_res.png](images/johnson_res.png)

In this case, there is an upstream head BC = 160 ft and a downstream head BC = 100 ft on the flat part of the 
downstream foundation. The entire back side of the dam is an exit face BC. Again, this is a partially saturated problem.

The following file illustrates how to prepare the inputs:

[xslope_johnson_res.xlsx](files/xslope_johnson_res.xlsx)

The solution should look something like this:

![johnson_res_solution.png](images/johnson_res_solution.png){width=1200px}

<!-- test: file=files/xslope_johnson_res.xlsx, type=seep, expected_flowrate=1.953, tolerance=0.05 -->

### 5. Earth Dam with Core and Filter

This problem has the following cross-section:

![earth_dam2.png](images/earth_dam2.png)

In this case, there is a single upstream head BC = 60ft and the entire backside of the dam is an exit face BC. 

The following Excel file contains the problem inputs:

[xslope_earth_dam2.xlsx](files/xslope_earth_dam2.xlsx)

The solution should look something like this:

![earth_dam2_solution.png](images/earth_dam2_solution.png){width=1200px}

<!-- test: file=files/xslope_earth_dam2.xlsx, type=seep, expected_flowrate=1.269, tolerance=0.05 -->

### 6. Levee with Grouted Foundation

The following problem represents a levee underlain by a foundation with a grout curtain. 

![levee.png](images/levee.png)

The material properties of the soil layers are as follows:

|  Soil Layer   | K1 [m/day] | K2 [m/day] | $\alpha$   | kr0    | h0 |
|:-------------:|:----------:|:----------:|:----------:|:----------:|:----------:|
|     Levee     |    0.5     |    0.2     |    0       |   0.001 | -1 |
| Grout Curtain |    0.2     |    0.2     |    0       |   0.001 | -1 |
|  Foundation   |     2      |     1      |    0       |   0.001 | -1 |

The coordinate geometry is shown here:

![levee_coords.png](images/levee_coords.png)

The following file illustrates how to prepare the inputs. Unlike the other
samples, this one defines the geometry with the **`polygon` sheet** — each
material zone (levee, grout curtain, foundation) is entered as a closed polygon
rather than as profile lines:

[xslope_levee_poly.xlsx](files/xslope_levee_poly.xlsx)

Inputs plotted with the XSLOPE plot_inputs() function:

![levee_inputs.png](images/levee_inputs.png){width=1200px}

!!! note
    The finite-element mesh is shown overlaid on the material zones above because
    this model has already been run — the mesh is generated during the seepage
    solve and stored with the model. `plot_inputs()` draws the mesh whenever one is
    present; for a model that has not yet been meshed, only the polygon geometry is
    shown.

Solution:

![levee_results.png](images/levee_results.png){width=1200px}

<!-- test: file=files/xslope_levee_poly.xlsx, type=seep, expected_flowrate=1.431, tolerance=0.05 -->

### 7. Kozeny Dam (verification — retained for scrutiny)

This is a homogeneous earth dam on an impervious base with a horizontal toe drain,
built so the upstream face follows the **exact confocal-parabola equipotential** of
the Kozeny flow net. It is a verification problem: the discharge has the closed
form `q = k*y0`, with `y0 = sqrt(d^2 + h^2) - d`. For this geometry (focus-to-entry
`d = 48`, reservoir depth `h = 20`), `y0 = 4` so the analytical `q = 4.0` (k = 1).

The input file is built by `benchmarks/build_seep.py::build_kozeny_dam`:

[xslope_kozeny_dam.xlsx](files/xslope_kozeny_dam.xlsx)

Note: the FE discharge for a horizontal-toe-drain dam is **sensitive** to where the
drain-start boundary lands relative to the point where the free surface meets the
base (the parabola vertex). xslope brackets the closed form within about ±4%
(drain at the focus over-predicts ~+4%, drain at the vertex under-predicts ~-4%);
the default build puts the drain at the vertex and yields ~3.81 on the standard
test mesh. This sample is retained so that sensitivity can be studied further; the
exact (confined) seepage anchor is the confined-radial case in
`benchmarks/build_seep.py`.

![kozeny_dam_solution.png](images/kozeny_dam_solution.png){width=1100}

<!-- test: file=files/xslope_kozeny_dam.xlsx, type=seep, expected_flowrate=3.81, tolerance=0.05, benchmark=SEEP-1b -->

### 8. Confined Radial Flow (verification — exact)

A quarter-annulus confined flow problem: inner arc (r = 10) at head 30, outer arc
(r = 30) at head 10, straight radial edges as no-flow streamlines. Steady saturated
flow has the exact solution `q = k*(pi/2)*(h1-h2)/ln(r2/r1) = 28.596` (k = 1) and a
logarithmic head profile. Built via the **polygon** input by
`benchmarks/build_seep.py::build_confined_radial`; this is the analytical anchor for
the seepage verification suite.

[xslope_confined_radial.xlsx](files/xslope_confined_radial.xlsx)

![confined_radial_solution.png](images/confined_radial_solution.png){width=800}

<!-- test: file=files/xslope_confined_radial.xlsx, type=seep, expected_flowrate=28.596, element_type=tri6, target_size=2.0, tolerance=0.01, benchmark=SEEP-1 -->

### 9. Partially Penetrating Sheetpile (verification — exact)

A single sheetpile cutoff of depth s = 10 in a homogeneous confined stratum of
thickness T = 20, head loss H = 10 across the wall, k = 1. Pavlovsky's exact
conformal-mapping solution gives `q = k*H*K(lam')/(2*K(lam))` with
`lam = sin(pi*s/(2T))`; at s/T = 1/2 the modulus is self-dual so **q = k*H/2 = 5.0
exactly**. A second exact check: the head on the wall plane below the tip is
exactly (H1+H2)/2 by antisymmetry. The wall uses the same V-notch crack idiom as
the clay-blanket sample. Built by `benchmarks/build_seep.py::build_sheetpile`
(s/T = 0.75 companion file also available, exact q = 3.4032).

[xslope_sheetpile_s50.xlsx](files/xslope_sheetpile_s50.xlsx)

The flow net (head contours and flowlines) for the s/T = 0.5 case:

![sheetpile_s50_solution.png](images/sheetpile_s50_solution.png){width=1100}

<!-- test: file=files/xslope_sheetpile_s50.xlsx, type=seep, expected_flowrate=5.0, element_type=tri6, target_size=1.0, tolerance=0.01, benchmark=SEEP-1c -->