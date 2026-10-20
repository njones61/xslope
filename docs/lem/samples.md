# Sample Problems - Limit Equilibrium Method

These problems cover common limit equilibrium setups: a slope with a foundation, a
submerged slope, layered soils, an earth dam, a tension crack, saturated unit weights,
and a pile-stabilized slope. Each one links its Excel input file. Open the file in
[XSLOPE Studio](../studio/index.md) with **File → Open** and run it, or load it in
Python with `load_slope_data()`. The
[limit equilibrium tutorials](../tutorials/index.md#limit-equilibrium) build models
like these step by step, and comparisons with published solutions and commercial
programs are in [Verification and Validation](../verification/index.md).

For each problem below, the solution figure shows the critical surface and factor of safety for Spencer's method, and most problems carry a **Factor of safety by method** table with the result for every method. On a solution figure, the green bars on the base of each slice are the effective stress there, the red bars are tension, and the red dashed line is the line of thrust computed with Spencer's method. The tables follow these conventions:

- Each value is that method's **own** critical surface — every method runs its own search, so the surfaces (and therefore the factors of safety) are not identical between methods.
- The methods differ by how much equilibrium they satisfy: OMS is the most approximate (and usually the most conservative), while Bishop, Janbu (corrected), Spencer, the Corps of Engineers method, and Lowe-Karafiath each enforce more of the force/moment balance. The Corps and Lowe-Karafiath force-equilibrium methods are sensitive to the assumed interslice-force inclination and can fall above the rigorous Spencer value.
- For **purely cohesive** soils ($\phi = 0$), the methods are theoretically identical on any given surface. Small differences in those tables therefore come from each method's search settling on a slightly different critical surface, not from the methods themselves.

### 1. Simple Slope with Foundation

This problem involves a uniform material extending below the toe of the slope. 

![Simple slope with a foundation](sample_images/simple_foundation_problem_sketch.png){width=1000}

Excel input file: [xslope_simple_foundation.xlsx](files/xslope_simple_foundation.xlsx) 

Inputs plotted with the XSLOPE plot_inputs() function:

![simple_foundation_inputs.png](sample_images/simple_foundation_inputs.png){width=1000}

Solution (critical surface and factor of safety):

![simple_foundation_results.png](sample_images/simple_foundation_results.png){width=1000}

<!-- fs-table -->
**Factor of safety by method** (each method's own critical surface):

| OMS | Bishop | Janbu | Corps | Lowe | Spencer | M-P |
|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| 0.964 | 0.964 | 1.029 | 1.120 | 1.041 | 0.964 | 0.964 |
<!-- /fs-table -->

<!-- test: file=files/xslope_simple_foundation.xlsx, type=circular_search, num_slices=40, fs_oms=0.964, fs_bishop=0.964, fs_janbu=1.029, fs_corps=1.120, fs_lowe=1.041, fs_spencer=0.964, fs_mprice=0.964 -->

### 2. Submerged Slope

This problem features a slope submerged by 10 ft of water. 

![Submerged slope](sample_images/submerged_slope_problem_sketch.png){width=1000}

The submerged slope is analyzed by applying a distributed load over the entire slope based on the unit weight of 
water (62.4 lb/ft3) and the depth of the water at a particular point on the slope. 

Excel input file: [xslope_submerged.xlsx](files/xslope_submerged.xlsx){width=900}

Inputs plotted with the XSLOPE plot_inputs() function:

![submerged_slope_inputs.png](sample_images/submerged_slope_inputs.png){width=1000}

Solution (critical surface and factor of safety):

![submerged_slope_results.png](sample_images/submerged_slope_results.png){width=1000}

<!-- fs-table -->
**Factor of safety by method** (each method's own critical surface):

| OMS | Bishop | Janbu | Corps | Lowe | Spencer | M-P |
|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| 1.154 | 1.154 | 1.248 | 2.011 | 1.861 | 1.154 | 1.154 |
<!-- /fs-table -->

<!-- test: file=files/xslope_submerged.xlsx, type=circular_search, num_slices=40, fs_oms=1.154, fs_bishop=1.154, fs_janbu=1.248, fs_corps=2.011, fs_lowe=1.861, fs_spencer=1.154, fs_mprice=1.154 -->

### 3. Slope with Eight Layers

This problem features a slope with eight soil layers. This problem was featured in the user manual for the UTEXASED 
slope stability analysis software developed by Stephen G. Wright at the University of Texas at Austin. It 
has a series of alternating layers, some of which are analyzed with an effective stress analysis and a 
piezometric line, and some of which are analyzed using a total stress analysis. We will assume that the base (max 
depth) is 10 ft below the top of the bottom material.

![Slope with eight layers](sample_images/eight_layers_problem_sketch.png){width=1000}

In the input file, the slope face rises 20 ft over a 45-ft run (2.25H:1V), and the water
table / piezometric line is horizontal at 2 ft below the toe-level ground surface
(elevation −2, inside the top foundation layer).

To find the critical surface and the global minimum factor of safety, we must use a circle starting at the base of 
each layer. The following Excel input file illustrates the problem.

Excel input file: [xslope_eight_layers.xlsx](files/xslope_eight_layers.xlsx)

Inputs plotted with the XSLOPE plot_inputs() function:

![eight_layers_inputs.png](sample_images/eight_layers_inputs.png){width=1000}

Search results:

![eight_layers_search_results.png](sample_images/eight_layers_search_results.png){width=1000}

Solution (critical surface and factor of safety):

![eight_layers_results.png](sample_images/eight_layers_results.png){width=1000}

<!-- fs-table -->
**Factor of safety by method** (each method's own critical surface):

| OMS | Bishop | Janbu | Corps | Lowe | Spencer | M-P |
|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| 0.805 | 1.154 | 1.160 | 1.240 | 1.060 | 1.189 | 1.170 |
<!-- /fs-table -->

<!-- test: file=files/xslope_eight_layers.xlsx, type=circular_search, num_slices=40, fs_oms=0.805, fs_bishop=1.154, fs_janbu=1.160, fs_corps=1.240, fs_lowe=1.060, fs_spencer=1.189, fs_mprice=1.170 -->

A grid search on this model (`seed='grid'`, which sweeps a coarse grid of centers and
tangent depths before refining) reaches circles that the search from the file's circles
never visits, and on one of them the **force-equilibrium** equation has no root that
corresponds to a physically possible sliding mass. Its ten roots there range from 0.12 to
2.55, while Bishop, Spencer and Morgenstern-Price all give about 9.1 on the same slices.
The Corps of Engineers method rejects that surface with a message giving the reason
instead of reporting 0.12 (see
[Force Equilibrium Methods](force_eq.md#solving-for-the-factor-of-safety)), and its
grid-seeded search reports 1.240, the same as the table above.

### 4. Earth Dam (Duncan, Wright, and Brandon)

This problem shows how a piezometric line can overstate the factor of safety where
water flows upward through a foundation. The downstream slope of an earth dam is
analyzed twice: once with pore pressures from a piezometric line and once with pore
pressures from a finite-element seepage solution. Under the downstream toe, water
rises through the foundation clay, and a piezometric line cannot represent that flow.
On the critical circle at the toe, Spencer's factor of safety is 1.323 with the
piezometric line and 0.967 with the seepage solution.

The problem is Example 2 in Section 6.13.6 of *Soil Strength and Slope Stability*,
2nd edition, by Duncan, Wright, and Brandon (2014), pages 121–123. The book gives the
elevations, the face slopes, the material properties, and the permeabilities; the
other dimensions are scaled from its Figures 6.41 and 6.42.
[Slide2's verification problem 72](../verification/rocscience.md#vp72) is the same
example, on Slide's own drawing of the section.

![Earth dam with a clay core](sample_images/earth_dam_problem_sketch.png){width=1000}

The book gives one unit weight per zone, γ, and the downstream analyses use it. Only
the seepage solution uses the permeabilities, k. The saturated unit weights, γ_sat,
are for the upstream analysis at the end of this section.

|  Mat   | c' (psf) | $\phi$' (degrees) | γ (pcf) | γ_sat (pcf) | k (ft/min) |
|:------:|:--------:|:-----------------:|:-------:|:-----------:|:----------:|
| Shell  |    0     |        34         |   125   |     130     |  1 × 10⁻²  |
|  Core  |   100    |        26         |   122   |     127     |  1 × 10⁻⁶  |
|  Clay  |    0     |        24         |   123   |     128     |  1 × 10⁻⁵  |
|  Sand  |    0     |        32         |   127   |     132     |  1 × 10⁻³  |

The piezometric line is the zero-pressure line of the book's seepage solution (its
Figure 6.42), traced from the figure. It stays at the reservoir level up to the core,
drops through the core, and meets the downstream face 1.5 ft above the toe.

**Downstream side of the dam**

The sand at the bottom of the foundation is a hundred times more permeable than the
clay above it, so reservoir water reaches the downstream toe through the sand still at
a high head, and its only way out is up through the clay. Upward flow raises the pore
pressure above that of still water at the same depth, which lowers the effective stress
and the strength of the shell and the clay, neither of which has cohesion. A piezometric
line gives the pressure of still water, so there it gives too little.

As in the book, each set of pore pressures gets two minimum factors of safety, with
each method running its own search:

- the overall minimum, a shallow circle at the toe, from a grid of starting circles,
  with the circle leaving the ground downstream of x = 500 ft:
  `circular_search(slope_data, method, seed='grid', exit_range=(500, 638))`
- the minimum for circles tangent to El. 197, the bottom of the clay:
  `circular_search(slope_data, method, tangent_depth=(197, 197))`

*Pore pressures from the piezometric line*

Excel input file: [xslope_earth_dam_down.xlsx](files/xslope_earth_dam_down.xlsx)

Inputs plotted with the XSLOPE plot_inputs() function:

![earth_dam_down_inputs.png](sample_images/earth_dam_down_inputs.png){width=1000}

Toe circle, shown close up:

![earth_dam_down_toe_results.png](sample_images/earth_dam_down_toe_results.png){width=922}

<!-- fs-table -->
**Factor of safety by method** (each method's own critical surface):

| OMS | Bishop | Janbu | Corps | Lowe | Spencer | M-P |
|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| n/a\* | 1.327 | 1.224 | 1.381 | 1.327 | 1.323 | 1.324 |

\* OMS is not reported for this problem's toe circle. The circle is small, and its last slices, where it rises steeply to the ground beyond the toe, have bases inclined up to 57°. There the Fellenius base normal force, $W\cos\alpha + D\cos(\alpha-\beta) - u\,\Delta\ell$, goes negative, on 6 of the 40 slices of the surface it settles on, so the shear resistance it computes there is meaningless and the factor of safety it reports (0.845) sits far below every other method's. See the [OMS method note](oms.md).
<!-- /fs-table -->

<!-- test: file=files/xslope_earth_dam_down.xlsx, type=circular_search, num_slices=40, seed=grid, exit_range=500;638, na=oms, fs_bishop=1.327, fs_janbu=1.224, fs_corps=1.381, fs_lowe=1.327, fs_spencer=1.323, fs_mprice=1.324 -->

Circle tangent to El. 197:

![earth_dam_down_deep_results.png](sample_images/earth_dam_down_deep_results.png){width=1000}

<!-- fs-table -->
**Factor of safety by method** (each method's own critical surface):

| OMS | Bishop | Janbu | Corps | Lowe | Spencer | M-P |
|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| 1.394 | 1.572 | 1.537 | 1.609 | 1.579 | 1.568 | 1.570 |
<!-- /fs-table -->

<!-- test: file=files/xslope_earth_dam_down.xlsx, type=circular_search, num_slices=40, tangent_depth=197;197, fs_oms=1.394, fs_bishop=1.572, fs_janbu=1.537, fs_corps=1.609, fs_lowe=1.579, fs_spencer=1.568, fs_mprice=1.570 -->

*Pore pressures from the seepage solution*

The seepage solution uses the permeabilities above, with the reservoir at El. 302
upstream, and the downstream face and the ground beyond the toe free to discharge. The
vertical ends of the section and its base are impermeable, so the water that flows
through the sand has to rise through the clay. The mesh has 3.5 ft triangles.

Excel input file: [xslope_earth_dam_down_seep.xlsx](files/xslope_earth_dam_down_seep.xlsx)
(the seepage mesh and solution are bundled alongside it).

Toe circle, shown close up:

![earth_dam_down_seep_toe_results.png](sample_images/earth_dam_down_seep_toe_results.png){width=874}

<!-- fs-table -->
**Factor of safety by method** (each method's own critical surface):

| OMS | Bishop | Janbu | Corps | Lowe | Spencer | M-P |
|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| n/a\* | 0.938 | 0.849 | 1.044 | 0.942 | 0.967 | 0.966 |

\* OMS is not reported for this problem's toe circle. The circle is small, and its last slices, where it rises steeply to the ground beyond the toe, have bases inclined up to 58°. There the Fellenius base normal force, $W\cos\alpha + D\cos(\alpha-\beta) - u\,\Delta\ell$, goes negative, on 10 of the 40 slices of the surface it settles on, so the shear resistance it computes there is meaningless and the factor of safety it reports (0.355) sits far below every other method's. See the [OMS method note](oms.md).
<!-- /fs-table -->

<!-- test: file=files/xslope_earth_dam_down_seep.xlsx, type=circular_search, num_slices=40, seed=grid, exit_range=500;638, na=oms, fs_bishop=0.938, fs_janbu=0.849, fs_corps=1.044, fs_lowe=0.942, fs_spencer=0.967, fs_mprice=0.966 -->

Beneath the toe, the head at the top of the sand is 245 ft, 18 ft above the ground, and
it is lost through the 30 ft of clay at an upward gradient of 0.54 to 0.62. The
gradient at which upward flow would carry the clay's whole buoyant weight is
$(\gamma - \gamma_w)/\gamma_w = (123 - 62.4)/62.4 = 0.97$. Along the base of the
circle beyond the toe, the seepage pore pressures are 1.65 to 1.74 times the still-water
pressure the piezometric line gives there. Over the whole base, the water force pushing
up on the slide mass is 4,400 lb/ft from the seepage solution and 3,160 lb/ft from the
piezometric line; the soil above the base weighs 6,920 lb/ft.

Circle tangent to El. 197:

![earth_dam_down_seep_deep_results.png](sample_images/earth_dam_down_seep_deep_results.png){width=1000}

<!-- fs-table -->
**Factor of safety by method** (each method's own critical surface):

| OMS | Bishop | Janbu | Corps | Lowe | Spencer | M-P |
|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| 1.250 | 1.406 | 1.391 | 1.431 | 1.412 | 1.405 | 1.405 |
<!-- /fs-table -->

<!-- test: file=files/xslope_earth_dam_down_seep.xlsx, type=circular_search, num_slices=40, tangent_depth=197;197, fs_oms=1.250, fs_bishop=1.406, fs_janbu=1.391, fs_corps=1.431, fs_lowe=1.412, fs_spencer=1.405, fs_mprice=1.405 -->

This circle crosses the same upward flow, but most of its base lies under the
embankment, where the weight of the fill holds the effective stress at 78 percent or
more of the piezometric line's value wherever at least 27 ft of fill lies above the
foundation. Its factor of safety (Spencer) falls by about a tenth, from 1.568 to
1.405, where the toe circle's falls by more than a quarter.

**Comparison with the book**

Spencer's factors of safety from the tables above, beside the book's Table 6.8. The book
does not say which method it used.

| Pore pressures | Toe circle, book | Toe circle, XSLOPE | Tangent to El. 197, book | Tangent to El. 197, XSLOPE |
|:--|:---:|:---:|:---:|:---:|
| Piezometric line | 1.30 | 1.323 | 1.57 | 1.568 |
| Seepage solution | 1.11 | 0.967 | 1.37 | 1.405 |

With the piezometric line, both circles agree with the book. With the seepage solution,
the circle tangent to El. 197 is within 3 percent of the book's, and the toe circle is
lower. That circle is small, about 11 ft in radius, and sits where the upward flow
concentrates at the toe, so its factor of safety depends on how finely the mesh resolves
the toe. It falls as the mesh is refined and settles at about the 3.5 ft element size
used here. On both circles the seepage solution gives the lower factor of safety, and the book
concludes that where the flow is not predominantly horizontal, the pore pressures
should come from a seepage analysis.

**Upstream side of the dam**

The book analyzes only the downstream slope. The upstream slope is analyzed here under
the full reservoir, with the same piezometric line and starting circles on the upstream
side of the dam. Each zone takes its saturated unit weight, γ_sat, 5 pcf above γ, below
the piezometric line.

Excel input file: [xslope_earth_dam_up.xlsx](files/xslope_earth_dam_up.xlsx)

Inputs plotted with the XSLOPE plot_inputs() function:

![earth_dam_up_inputs.png](sample_images/earth_dam_up_inputs.png){width=1000}

Search results:

![earth_dam_up_search_results.png](sample_images/earth_dam_up_search_results.png){width=1000}

Solution (critical surface and factor of safety):

![earth_dam_up_results.png](sample_images/earth_dam_up_results.png){width=1000}

<!-- fs-table -->
**Factor of safety by method** (each method's own critical surface):

| OMS | Bishop | Janbu | Corps | Lowe | Spencer | M-P |
|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| n/a\* | 1.828 | 1.752 | 2.076 | 2.025 | 1.812 | 1.807 |

\* OMS is not reported for this problem. Its base normal force is the Fellenius value $W\cos\alpha + D\cos(\alpha-\beta) - u\,\Delta\ell$, which under a full reservoir goes negative on the deepest slices — a quarter of them on the surface it settles on here — so the shear resistance it computes there is meaningless and the factor of safety it reports (0.897) sits far below every other method's. See the [OMS method note](oms.md).
<!-- /fs-table -->

<!-- test: file=files/xslope_earth_dam_up.xlsx, type=circular_search, num_slices=40, fs_bishop=1.828, fs_janbu=1.752, fs_corps=2.076, fs_lowe=2.025, fs_spencer=1.812, fs_mprice=1.807 -->

### 5. Tension Crack

This problem is a slope whose upper layer has cohesion, so an analysis without a
crack produces non-physical **tension at the crest** (and an inverted line of thrust) that
unconservatively raises the factor of safety. The remedy is a **tension crack** at
the top of the slope, whose depth follows
$d_{crack} = \dfrac{2 c_d}{\gamma}\tan\!\left(45 + \dfrac{\phi_d}{2}\right)$ with the
mobilized strengths $c_d = c/F$, $\tan\phi_d = \tan\phi / F$. Because the crack
depth depends on $F$, it is iterated to convergence; this model carries the
converged depth (`tcrack_depth` = 4.5 ft on the **main** sheet), at which the crest
tension just vanishes. The crack is taken as full of water (`tcrack_water` = 4.5 ft). The complete-equilibrium methods agree (Spencer and
Morgenstern-Price both 1.414, matching Bishop).

![Slope with a tension crack](sample_images/tension_problem_sketch.png){width=1000}

Excel input file: [xslope_tension_KEY.xlsx](files/xslope_tension_KEY.xlsx)

Inputs plotted with the XSLOPE plot_inputs() function:

![tension_inputs1.png](sample_images/tension_inputs1.png){width=1000}

Solution (critical surface with the tension crack, Spencer's method):

![tension_results1.png](sample_images/tension_results1.png){width=1000}

<!-- fs-table -->
**Factor of safety by method** (each method's own critical surface):

| OMS | Bishop | Janbu | Corps | Lowe | Spencer | M-P |
|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| 1.413 | 1.414 | 1.448 | 1.673 | 1.544 | 1.414 | 1.414 |
<!-- /fs-table -->

<!-- test: file=files/xslope_tension_KEY.xlsx, type=circular_search, num_slices=40, fs_oms=1.413, fs_bishop=1.414, fs_janbu=1.448, fs_corps=1.673, fs_lowe=1.544, fs_spencer=1.414, fs_mprice=1.414 -->

### 6. Saturated vs. Moist Unit Weight (γ_sat)

Soil below the water table weighs more than the same soil above it. On the **mat**
sheet, `gamma` is the moist unit weight and `gamma_sat` the saturated unit weight;
when a material carries both, the water table splits each slice's weight — γ_sat
below the water table, γ above it. The water table comes from whichever the model
defines: a piezometric line on the **piezo** sheet, or the phreatic surface of a
finite-element seepage solution. The water table is part of the model rather than of
any one material, so it splits the weight whether or not the material's pore-pressure
option also takes pore pressures from it.

This problem features a slope in undrained clay ($S_u = 600$ psf, $\phi = 0$) with
an internal water table: γ = 120 pcf above the water table and γ_sat = 127 pcf
below. The strength is a total-stress strength, so pore pressure never enters the
analysis (`u = none`) and the piezometric line serves only to locate the water
table. Spencer's method gives $F = 1.420$.

![Clay slope with an internal water table](sample_images/gsat_problem_sketch.png){width=1000}

Excel input file: [xslope_gsat_sidecar.xlsx](files/xslope_gsat_sidecar.xlsx)

Solution (critical surface and factor of safety, Spencer's method):

![gsat_sidecar_results.png](sample_images/gsat_sidecar_results.png){width=1000}

<!-- test: file=files/xslope_gsat_sidecar.xlsx, type=circular_search, num_slices=40, fs_bishop=1.420, fs_spencer=1.420, fs_janbu=1.518 -->

The same slope built as two material zones split at the water table — 120 pcf clay
above, 127 pcf clay below — reproduces that factor of safety:
[xslope_gsat_zoned.xlsx](files/xslope_gsat_zoned.xlsx).

<!-- test: file=files/xslope_gsat_zoned.xlsx, type=gsat_pair, file2=files/xslope_gsat_sidecar.xlsx -->
<!-- test: file=files/xslope_gsat_zoned.xlsx, type=circular_search, num_slices=40, fs_spencer=1.420 -->

A model can also supply γ_sat with no water table anywhere
([xslope_gsat_nowater.xlsx](files/xslope_gsat_nowater.xlsx)); with no elevation
below which the saturated weight applies, the slope is analyzed with γ throughout.

<!-- test: file=files/xslope_gsat_nowater.xlsx, type=circular_search, num_slices=40, fs_bishop=1.468, fs_spencer=1.468, fs_janbu=1.565 -->

The same slope and water table with a drained strength ($c' = 200$ psf,
$\phi' = 25°$) and `u = piezo`: the piezometric line now supplies the pore
pressures as well as the unit-weight split.

Excel input file: [xslope_gsat_piezo.xlsx](files/xslope_gsat_piezo.xlsx)

Solution (critical surface and factor of safety, Spencer's method):

![gsat_piezo_results.png](sample_images/gsat_piezo_results.png){width=1000}

<!-- test: file=files/xslope_gsat_piezo.xlsx, type=circular_search, num_slices=40, fs_bishop=1.538, fs_spencer=1.539, fs_janbu=1.503 -->

The water table can also come from a seepage solution. This is the upstream slope
of the earth dam of [Problem 4](#4-earth-dam-duncan-wright-and-brandon), with the same saturated unit weights
and `u = seep`: pore pressures are interpolated from the finite-element
seepage solution, and the phreatic surface of that same solution — the heavy
contour carrying the water-table marker in the figure below — is the water table
that splits the weights. When a model carries both a seepage solution and a
piezometric line, the seepage solution sets the split.

Excel input file: [xslope_gsat_seep.xlsx](files/xslope_gsat_seep.xlsx) (the seepage
mesh and solution are bundled alongside it).

Solution (critical surface and factor of safety, with the seepage head contours):

![gsat_seep_results.png](sample_images/gsat_seep_results.png){width=1000}

<!-- test: file=files/xslope_gsat_seep.xlsx, type=circular_search, num_slices=40, fs_bishop=1.925, fs_spencer=1.906 -->

The same dam analyzed for [rapid drawdown](rapid.md), with the reservoir at
El. 302 ft before drawdown and El. 250 ft after. The premise of *rapid* drawdown is
that the pool falls faster than the low-permeability zones can drain, so the soil
stays saturated while the pore pressures fall. The slice weights use the
pre-drawdown water table in all three stages, while the pore pressures follow
the staged piezometric lines.

Excel input file: [xslope_gsat_rapid.xlsx](files/xslope_gsat_rapid.xlsx)

Solution (governing rapid-drawdown surface and factor of safety):

![gsat_rapid_results.png](sample_images/gsat_rapid_results.png){width=1000}

<!-- test: file=files/xslope_gsat_rapid.xlsx, type=circular_search, num_slices=40, rapid=true, fs_bishop=1.095, fs_spencer=1.094 -->

### 7. Pile-Stabilized Slope (Hassiotis et al. 1997)

This is a published pile-stabilization benchmark. It checks that XSLOPE's built-in
Ito & Matsui force reproduces the force used in the source's design. The slope is
homogeneous, 13.7 m high at 30°, with $c = 23.94$ kPa, $\phi = 10°$ and
$\gamma = 19.63$ kN/m³, and is dry. One row of 1.0 m piles at 2.5 m centers is
placed 13.7 m horizontally from the toe, and in a second case 23.1 m from the toe.
The published clear-to-center spacing ratio $D_2/D_1 = 0.6$ is reproduced exactly:
a 1.0 m pile at 2.5 m centers leaves a 1.5 m clear opening, and $1.5/2.5 = 0.6$.

![Pile-stabilized slope with the two pile stations](sample_images/hassiotis_problem_sketch.png){width=1000}

Inputs, with both pile stations drawn together:

![hassiotis_inputs.png](sample_images/hassiotis_inputs.png){width=1000}

Excel input files:
[xslope_hassiotis.xlsx](files/xslope_hassiotis.xlsx) (unreinforced),
[xslope_hassiotis_p1.xlsx](files/xslope_hassiotis_p1.xlsx) (row 13.7 m from the toe),
[xslope_hassiotis_p2.xlsx](files/xslope_hassiotis_p2.xlsx) (row 23.1 m from the toe).

| Property | Value |
|----------|-------|
| Slope height, $H$ | 13.7 m |
| Slope angle, $\beta$ | 30 degrees |
| Cohesion, $c$ | 23.94 kPa |
| Friction angle, $\phi$ | 10 degrees |
| Unit weight, $\gamma$ | 19.63 kN/m³ |
| Pile diameter, $D$ | 1.0 m |
| Pile spacing, $S$ | 2.5 m |
| Pile length, $L$ | 17 m |
| $V_{\text{cap}}$, $M_{\text{cap}}$ | not specified |

The piles carry no structural capacity limits. The source specifies none, and its
factors of safety are the full soil force, so a cap would change the quantity
being compared. $H$ is left blank, so the Ito & Matsui force is computed for every
trial surface from $D$, $S$ and the soil above that surface at the pile.

#### Unreinforced (FS = 1.105)

![hassiotis_results.png](sample_images/hassiotis_results.png){width=1000}

Bishop's method gives 1.106 against the 1.12 Hull & Poulos (1999) report with the
same method (−1.2%). Hassiotis et al. report 1.08 and Ausilio et al. (2001) 1.11,
both by the friction-circle method, which XSLOPE does not implement.

<!-- test: file=files/xslope_hassiotis.xlsx, type=circular_search, num_slices=50, fs_oms=1.056, fs_bishop=1.106, fs_janbu=1.105, fs_corps=1.167, fs_lowe=1.135, fs_spencer=1.105, fs_mprice=1.104, benchmark=LEM-HASSIOTIS -->

#### Pile row 13.7 m from the toe (FS = 1.855)

![hassiotis_p1_results.png](sample_images/hassiotis_p1_results.png){width=1000}

Bishop's method gives 1.859 against 1.82 (Hassiotis et al., friction circle, +2.1%).
The row raises the factor of safety by 68%.

#### Pile row 23.1 m from the toe (FS = 1.281)

![hassiotis_p2_results.png](sample_images/hassiotis_p2_results.png){width=1000}

Bishop's method gives 1.289 against 1.64 (Hassiotis et al., −21%). The row sits
0.6 m short of the crest, so every surface that reaches it crosses it within a few
meters of the pile head, where the soil column above the surface — and with it the
Ito & Matsui force — is small. Moving the entry limit 2 m further behind the crest
moves this factor of safety to 1.48, so for this row the result depends on where the
search is allowed to start.

#### Search limits

Both pile files declare a search window on their circles sheet: the surface
daylights within a few meters of the toe (exit 25–32 m), enters behind the crest
(entry 54–75 m) and keeps its lowest point above the pile tip. Without it the
search returns a deep surface that passes *below* the pile tip and collects no pile
force at all: 1.327 for the 13.7 m row, on a deeper mechanism of the unreinforced
slope that the piles do not cross. The published comparisons are for the mechanism
through the row, and their searches are restricted in the same way.

#### Ito & Matsui summary

On the critical surface for the 13.7 m row, the surface crosses the pile 5.72 m
below its head. At $D = 1.0$, $S = 2.5$ and $\phi = 10°$ the coefficients are
$A_1 = 3.2298$ and $A_2 = 1.5695$, giving a pressure of 77.3 kN/m at the head
rising to 253.4 kN/m at the surface, a force of 945.4 kN **per pile**, and

$$H = \frac{F_{\text{pile}}}{S} = \frac{945.4}{2.5} = 378.2 \ \text{kN/m of slope}$$

which is the per-unit-width value the slice equations apply, horizontally, at the
pile–surface intersection.

The same coefficients reproduce the force used in the source's design. Hassiotis et
al. state 72.4 kN/m at the pile head and a 561.7 kN/m overburden term at 17 m
depth; XSLOPE's $c\,A_1 = 77.3$ kN/m (+6.8%) and $\gamma A_2 \cdot 17 = 523.8$ kN/m
(−6.7%). The two departures are in opposite directions and largely cancel in the
integral: over the published 6.56 m depth to the slip surface, the published
trapezoid gives 1185.7 kN per pile and XSLOPE's exact integration of the same law
gives 1170.2 kN, a difference of 1.3%.

#### Force direction

XSLOPE applies the pile reaction horizontally for a vertical pile and does not
divide it by the computed factor of safety (`Appl = active`). Hull & Poulos note
that the plastic-deformation theory derives the force horizontally, which is the
direction used here. Hassiotis et al. apply it parallel to the slip surface and
also leave it unfactored, so their 1.82 is the value this case is compared with.
Their 1.45 comes from a boundary-element shear and moment divided by the global
factor of safety. That is a different force model, and the value is quoted only to
show the range of published results.

<!-- fs-table -->
**Factor of safety by method** (each method's own critical surface, 13.7 m row):

| OMS | Bishop | Janbu | Corps | Lowe | Spencer | M-P |
|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| 1.822 | 1.859 | 1.788 | 1.885 | 1.871 | 1.855 | 1.855 |
<!-- /fs-table -->

<!-- test: file=files/xslope_hassiotis_p1.xlsx, type=circular_search, num_slices=50, entry_range=54.0;75.0, exit_range=25.0;32.0, tangent_depth=11.01;33.7, fs_oms=1.822, fs_bishop=1.859, fs_janbu=1.788, fs_corps=1.885, fs_lowe=1.871, fs_spencer=1.855, fs_mprice=1.855, benchmark=LEM-HASSIOTIS -->
<!-- test: file=files/xslope_hassiotis_p2.xlsx, type=circular_search, num_slices=50, entry_range=54.0;75.0, exit_range=25.0;32.0, tangent_depth=16.44;33.7, fs_oms=1.244, fs_bishop=1.289, fs_janbu=1.373, fs_corps=1.385, fs_lowe=1.340, fs_spencer=1.281, fs_mprice=1.289, benchmark=LEM-HASSIOTIS -->
