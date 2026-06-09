# Sample Problems - Limit Equilibrium Method

The following examples illustrate how to use XSLOPE to perform limit equilibrium slope stability analysis. Each of the Excel input files below can be uploaded and used with the following Google Colab notebook which has been set up specifically for running slope stability analyses:

<a href="https://colab.research.google.com/github/njones61/xslope/blob/main/notebooks/xslope_lem.ipynb" target="_"><img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open In Colab"/></a>

The notebook allows the user to select a variety of analysis options using simple form inputs and then runs the analysis using the selected method and plots the results.

### 1. Simple Embankment

This problem features a simple slope with a single material. 

![simple_embankment.png](sample_images/simple_embankment.png){width=700}

Excel input file: [xslope_simple_embankment.xlsx](files/xslope_simple_embankment.xlsx)

Inputs plotted with the XSLOPE plot_inputs() function:

![simple_embankment_inputs1.png](sample_images/simple_embankment_inputs1.png)

Solution (critical surface and factor of safety). The green bars on the base of each slice represent the effective 
stress on the base of the slice. The red bars correspond to tension at the base of the slice. The red dashed line 
represents the line of thrust computed using Spencer's method.

![simple_embankment_results1.png](sample_images/simple_embankment_results1.png){width=700}

<!-- test: file=files/xslope_simple_embankment.xlsx, type=circular_search, method=spencer, expected_fs=1.281, num_slices=40 -->

Here is copy of the input file with the following variations/changes:

a) Distributed load on top of slope. q = 750 psf<br>
b) Tension crack. Depth = 3 ft. <br>
c) Tension crack filled with water.<br>
d) Submerged by 10 ft depth of water (distributed load)

Excel input file: [xslope_simple_embankment_mods.xlsx](files/xslope_simple_embankment_mods.xlsx){width=700}

Inputs:

![simple_embankment_inputs2.png](sample_images/simple_embankment_inputs2.png){width=700}

Solution (critical surface and factor of safety):

![simple_embankment_results2.png](sample_images/simple_embankment_results2.png){width=700}

<!-- test: file=files/xslope_simple_embankment_mods.xlsx, type=circular_search, method=spencer, expected_fs=0.986, num_slices=30 -->

### 2. Simple Slope with Foundation

This problem involves a uniform material extending below the toe of the slope. 

![simple_foundation.png](sample_images/simple_foundation.png){width=700}

Excel input file: [xslope_simple_foundation.xlsx](files/xslope_simple_foundation.xlsx) 

Inputs plotted with the XSLOPE plot_inputs() function:

![simple_foundation_inputs.png](sample_images/simple_foundation_inputs.png){width=700}

Solution (critical surface and factor of safety):

![simple_foundation_results.png](sample_images/simple_foundation_results.png){width=700}

<!-- test: file=files/xslope_simple_foundation.xlsx, type=circular_search, method=spencer, expected_fs=0.964, num_slices=30 -->

### 3. Simple Slope with Multiple Layers

This problem involves a simple slope with multiple layers of material. 

![simple_mult_layers.png](sample_images/simple_mult_layers.png){width=700}

Excel input file: [xslope_simple_mult_layers.xlsx](files/xslope_simple_mult_layers.xlsx)

Inputs plotted with the XSLOPE plot_inputs() function. Note that in this case we use two starting circles - one at 
the base each of each of the two materials - to ensure that the search algorithm finds the critical surface 
corresponding to a global and not a local minimum. 

![simple_mult_layers_inputs.png](sample_images/simple_mult_layers_inputs.png){width=900}

Search results. Each gray line represent each circle used in the search. The dots represent the center of the 
circles used in the nine-point search algorithm, and the green arrows represent the path of grid centers taken to 
reach the critical surface. The red circle represents the critical surface.

![simple_mult_layers_search_results.png](sample_images/simple_mult_layers_search_results.png){width=900}

Solution (critical surface and factor of safety):

![simple_mult_layers_results.png](sample_images/simple_mult_layers_results.png){width=900}

<!-- test: file=files/xslope_simple_mult_layers.xlsx, type=circular_search, method=spencer, expected_fs=1.244, num_slices=30 -->

### 4. Submerged Slope

This problem features a slope submerged by 10 ft of water. 

![submerged_slope.png](sample_images/submerged_slope.png){width=600}

The submerged slope is analyzed by applying a distributed load over the entire slope based on the unit weight of 
water (62.4 lb/ft3) and the depth of the water at a particular point on the slope. 

Excel input file: [xslope_submerged.xlsx](files/xslope_submerged.xlsx){width=900}

Inputs plotted with the XSLOPE plot_inputs() function:

![submerged_slope_inputs.png](sample_images/submerged_slope_inputs.png){width=900}

Solution (critical surface and factor of safety):

![submerged_slope_results.png](sample_images/submerged_slope_results.png){width=900}

<!-- test: file=files/xslope_submerged.xlsx, type=circular_search, method=spencer, expected_fs=1.154, num_slices=30 -->

### 5. Slope with Multiple Materials and Piezometric Line

This problem features three layers of material with an effective stress analysis where pore pressures are derives 
from a piezometric line. 

![method_slices_problem.png](sample_images/method_slices_problem.png){width=900}

This problem is featured as part of a graduate course on slope stability analysis (CE 544 - Slope Stability Analysis)
at Brigham Young University. The problem used in two exercises to illustrate how to solve limit equilibrium slope 
stability problems using the method of slices and an Excel spreadsheet. The problem descriptions are here:

[Ordinary Method of Slices Exercise](https://byu-ce544.readthedocs.io/en/latest/unit2/04_limiteq2/limiteq2_class/)<br>
[Bishop Simplified Procedure Homework](https://byu-ce544.readthedocs.io/en/latest/unit2/04_limiteq2/limiteq2_hw/)

In these exercises, a single circular surface was analyzed. The following Excel input file illustrates the problem:

Excel input file: [xslope_method_slices_problem.xlsx](files/xslope_method_slices_problem.xlsx)

Inputs plotted with the XSLOPE plot_inputs() function:

![method_slices_problem_inputs.png](sample_images/method_slices_problem_inputs.png){width=900}

Here is the solution for just the starting circle (to match the problem description) using Bishop's simplified procedure:

![method_slices_problem_results.png](sample_images/method_slices_problem_results.png){width=900}

<!-- test: file=files/xslope_method_slices_problem.xlsx, type=single_circle, method=bishop, expected_fs=1.576, num_slices=30 -->

Here is the Excel input file with multiple starting circles for a global search for the critical surface:

Excel input file: [xslope_method_slices_problem2.xlsx](files/xslope_method_slices_problem2.xlsx)

Inputs plotted with the XSLOPE plot_inputs() function:

![method_slices_problem_inputs2.png](sample_images/method_slices_problem_inputs2.png){width=900}

Sarch results. This problem is a good example of the search path and the large number of circles that are sometimes 
tested in the search algorithm. In this case, the critical surface is isolated to sloughing of the 2nd layer.

![method_slices_problem_search_results2.png](sample_images/method_slices_problem_search_results2.png){width=900}

Solution (critical surface and factor of safety):

![method_slices_problem_results2.png](sample_images/method_slices_problem_results2.png){width=900}

<!-- test: file=files/xslope_method_slices_problem2.xlsx, type=circular_search, method=spencer, expected_fs=0.770, num_slices=30 -->

### 6. Slope with Eight Layers

This problem features a slope with eight soil layers. This problem was featured in the user manual for the UTEXASED 
slope stability analysis software developed by  at the University of Texas at Austin by Stephen G. Wright. If 
features a series of alternating layers, some of which are analyzed with an effective stress analysis and a 
piezometric line, and some of which are analyzed using a total stress analysis. We will assume that the base (max 
depth) is 10 ft below the top of the bottom material.

![eight_layers.png](sample_images/eight_layers.png){width=900}

To find the critical surface and the global minimum factor of safety, we must use a circle starting at the base of 
each layer. The following Excel input file illustrates the problem.

Excel input file: [xslope_eight_layers.xlsx](files/xslope_eight_layers.xlsx)

Inputs plotted with the XSLOPE plot_inputs() function:

![eight_layers_inputs.png](sample_images/eight_layers_inputs.png){width=900}

Search results:

![eight_layers_search_results.png](sample_images/eight_layers_search_results.png){width=900}

Solution (critical surface and factor of safety):

![eight_layers_results.png](sample_images/eight_layers_results.png){width=900}

<!-- test: file=files/xslope_eight_layers.xlsx, type=circular_search, method=spencer, expected_fs=1.189, num_slices=30 -->

### 7. Non-Circular Failure Surface

This problem features a thin weak layer in the foundation of a slope. In such cases, a non-circular failure surface 
constrained to fit in the weak layer often corresponds to the critical failure surface. This can be modeled with 
non-circular options in XSLOPE. This problem is also featured in the user manual for the UTEXASED slope stability 
analysis software developed by Stephen G. Wright at the University of Texas at Austin.

![noncircular.png](sample_images/noncircular.png){width=900}

The non-circular failure surface is modeled with the following Excel input file. The failure surface is defined by 
four points. The first and last point are assigned the "Free" option, which causes them to be automatically 
calculated based on the slope geometry. The two middle points are assigned the "Horiz" option, which causes them to 
be moved horizontally inside the weak layer.

Excel input file: [xslope_noncircular.xlsx](files/xslope_noncircular.xlsx)

Inputs plotted with the XSLOPE plot_inputs() function:

![noncircular_inputs.png](sample_images/noncircular_inputs.png){width=900}

Search results:

![noncircular_search_results.png](sample_images/noncircular_search_results.png){width=900}

!!! note
    The search algorithm for non-circular failure surfaces is highly sensitive to the starting location. It the 
    angle of the wedge at the toe of the slope is too steep, there will be tension at the toe of the slope and the search 
    will fail to find a correct solution.

Solution (critical surface and factor of safety):

![noncircular_results.png](sample_images/noncircular_results.png){width=900}

<!-- test: file=files/xslope_noncircular.xlsx, type=noncircular_search, method=spencer, expected_fs=1.739, num_slices=30 -->

### 8. Earth Dam 

This problem features a dam with a shell and a clay core on top of a foundation with a clay layer and a sand layer. 
This problem was featured on page 121 of Shear Strength and Slope Stability - Second Edition by Duncan, Wright, and 
Brandon. 

![earth_dam.png](sample_images/earth_dam.png)

The material properties are as follows:

|  Mat   | c' (psf) | $\phi$' (degrees) | γ (pcf) |
|:------:|:--------:|:-----------------:|:-------:|
| Shell  |    0     |        34         |   125   |
|  Core  |   100    |        26         |   122   |
|  Clay  |    0     |        24         |   123   |
|  Sand  |    0     |        32         |   127   |

**Upstream side of the dam**

First, we will analyze the upstream side. This is accomplished by defining starting circles on the upstream side of the 
dam. The following Excel input file illustrates the problem.

Excel input file: [xslope_earth_dam_up.xlsx](files/xslope_earth_dam_up.xlsx)

Inputs plotted with the XSLOPE plot_inputs() function:

![earth_dam_up_inputs.png](sample_images/earth_dam_up_inputs.png){width=900}

Search results:

![earth_dam_up_search_results.png](sample_images/earth_dam_up_search_results.png){width=900}

Solution (critical surface and factor of safety):

![earth_dam_up_results.png](sample_images/earth_dam_up_results.png){width=900}

<!-- test: file=files/xslope_earth_dam_up.xlsx, type=circular_search, method=spencer, expected_fs=1.800, num_slices=30 -->

**Downstream side of the dam**

Next, we will analyze the other side of the dam by defining starting circles on the downstream 
side of the dam. 

Excel input file: [xslope_earth_dam_down.xlsx](files/xslope_earth_dam_down.xlsx)

Inputs plotted with the XSLOPE plot_inputs() function:

![earth_dam_down_inputs.png](sample_images/earth_dam_down_inputs.png){width=900}

Search results:

![earth_dam_down_search_results.png](sample_images/earth_dam_down_search_results.png){width=900}

Solution (critical surface and factor of safety):

![earth_dam_down_results.png](sample_images/earth_dam_down_results.png){width=900}

<!-- test: file=files/xslope_earth_dam_down.xlsx, type=circular_search, method=spencer, expected_fs=1.558, num_slices=30 -->

### 9. Reinforced Slope

This problem features an engineered slope with six layers of geogrid reinforcement. This problem was featured in the 
user manual for the UTEXASED slope stability analysis software developed by Stephen G. Wright at the University of Texas 
at Austin.

![reinforce.png](sample_images/reinforce.png){width=900}

A 240 psf surcharge is applied along the slope crest. For each line of reinforcement, the full tensile force develops over a length of 4 ft. The toe of the slope corresponds
to (0, 0) and the top of the slope corresponds to (30, 24). 

The following Excel input file illustrates the problem. The soil reinforcement is entered in the "reinforce" sheet.

Excel input file: [xslope_reinforce.xlsx](files/xslope_reinforce.xlsx)

Inputs plotted with the XSLOPE plot_inputs() function:

![reinforce_inputs.png](sample_images/reinforce_inputs.png){width=900}

Solution (critical surface and factor of safety):

![reinforce_results.png](sample_images/reinforce_results.png){width=900}

<!-- test: file=files/xslope_reinforce.xlsx, type=circular_search, method=spencer, expected_fs=1.589, num_slices=30 -->

!!! note
    The solution for this problem found by XSLOPE is not the same as the solution found by UTEXASED. The difference
    is due to the fact that XSLOPE assumes the reinforcement is flexible and the force from the reinforcement is
    therefore parallel to the base of the slope. UTEXASED assumes the reinforcement is rigid and the force from the
    reinforcement is in the direction of the reinforcement line. The flexible assumption is more conservative. The
    UTEXASED solution for this problem is FS = 1.646.

### 10. Slope Stabilized with Piles

This problem features a 1:1 slope in a medium-stiff clay stabilized by two rows of drilled shafts.

![piles_inputs.png](sample_images/piles_inputs.png){width=900}

Excel input file: [xslope_piles.xlsx](files/xslope_piles.xlsx)

| Property | Value |
|----------|-------|
| Cohesion, $c$ | 200 psf |
| Friction angle, $\phi$ | 20 degrees |
| Unit weight, $\gamma$ | 120 pcf |
| Pile diameter, $D$ | 2.0 ft |
| Pile spacing, $S$ | 6.0 ft |
| $V_{\text{cap}}$ | 46,000 lb |
| $M_{\text{cap}}$ | 60,000 ft·lb |

#### Results Without Piles (FS = 1.15)

![piles_results_no_pile.png](sample_images/piles_results_no_pile.png){width=900}

#### Results With Piles (FS = 1.85)

![piles_results.png](sample_images/piles_results.png){width=900}

The two rows of piles increase the factor of safety from 1.15 to 1.85.

#### Ito & Matsui Summary

The pile force $H$ is not specified directly in the input file. Instead, XSLOPE auto-computes $H$ using the
Ito & Matsui (1975) method, which models the plastic flow of soil between adjacent piles to determine the lateral
resistance. Because $H$ is computed for each trial failure surface during the search, the pile resistance varies
with the depth of the failure surface at the pile location.

Structural capacity limits ($V_{\text{cap}}$ = 46,000 lb, $M_{\text{cap}}$ = 60,000 ft·lb) are specified for
each pile, consistent with a 2-ft diameter reinforced concrete section ($f'_c$ = 4000 psi). For the critical
failure surface, the Ito & Matsui soil forces far exceed the structural capacity, and the moment capacity
controls:

```text
  === Ito & Matsui Summary (Pile 1) ===
  Pile diameter (D)          = 2.0
  Pile spacing (S)           = 6.0
  Clear spacing (D1 = S - D) = 4.0
  Depth to failure surface   = 9.5
  Coefficients: A1 = 7.569, A2 = 4.755
  Force per pile (F_pile)    = 39810
  Force per unit width (H)   = 6635.1
  --- Structural Capacity Check ---
  V_cap = 46000  (F_pile within shear capacity)
  M_cap = 60000, L_m = 3.72, F_limit = M_cap/L_m = 16139  (F_pile exceeds moment capacity)
  Controlled by moment (M_cap/L_m)
  F_pile: 39810 -> 16139 (capped)
  H:      6635.1 -> 2689.8 (capped)

  === Ito & Matsui Summary (Pile 2) ===
  Pile diameter (D)          = 2.0
  Pile spacing (S)           = 6.0
  Clear spacing (D1 = S - D) = 4.0
  Depth to failure surface   = 13.9
  Coefficients: A1 = 7.569, A2 = 4.755
  Force per pile (F_pile)    = 76447
  Force per unit width (H)   = 12741.2
  --- Structural Capacity Check ---
  V_cap = 46000  (F_pile exceeds shear capacity)
  M_cap = 60000, L_m = 5.28, F_limit = M_cap/L_m = 11356  (F_pile exceeds moment capacity)
  Controlled by moment (M_cap/L_m)
  F_pile: 76447 -> 11356 (capped)
  H:      12741.2 -> 1892.6 (capped)
```

The Ito & Matsui soil forces (39,810 and 76,447 lb per pile) represent the theoretical upper bound on what the soil
can push onto the pile. These greatly exceed both the shear and moment capacities. After capping, the effective
pile forces are reduced by 59% and 85% respectively, with the moment capacity ($M_{\text{cap}} / L_m$) controlling
in both cases. Without the capacity checks, the LEM would overestimate the pile resistance and produce an
unconservatively high factor of safety.

#### LEM vs. FEM Comparison

The corresponding FEM analysis of this problem (see [FEM Samples](../fem/samples.md), Problem 3) gives
FS = 1.32 with piles — significantly lower than the LEM result of 1.85. This difference arises because
the LEM applies the Ito & Matsui force (even after capping) as a concentrated load at the failure surface,
while the FEM beam elements only develop as much resistance as the global deformation pattern naturally
produces. The FEM result is generally considered more realistic for pile-stabilized slopes.

<!-- test: file=files/xslope_piles.xlsx, type=circular_search, method=spencer, expected_fs=1.85, num_slices=40 -->

### 11. Slope on a Sloping (Irregular) Bottom

This problem demonstrates **polygon-based geometry** with an irregular lower boundary.
Instead of profile lines and a horizontal `max_depth`, the cross-section is defined
directly as two material-zone polygons — an embankment over a foundation — whose shared
base dips from left to right (from elevation −15 on the left to −5 on the right). The
failure surface is constrained to stay within the domain polygon: trial circles that would
dip below the sloping base are automatically rejected during the search, so the critical
surface follows the dipping bedrock rather than a fictitious flat cutoff.

Excel input file: [xslope_sloping_bottom.xlsx](files/xslope_sloping_bottom.xlsx)

<!-- test: file=files/xslope_sloping_bottom.xlsx, type=circular_search, method=spencer, expected_fs=1.244, num_slices=30 -->
