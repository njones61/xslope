# XSLOPE

**xslope** is a comprehensive Python package for geotechnical slope stability and seepage analysis. Available on [PyPI](https://pypi.org/project/xslope/), it provides integrated tools for limit equilibrium slope stability analysis, finite element seepage analysis, and finite element-based slope stability analysis.

The package uses an Excel-based input template system, making it accessible to practitioners familiar with spreadsheet workflows while leveraging the power of Python for computational analysis and visualization.

## Limit Equilibrium Slope Stability Analysis

The limit equilibrium method (LEM) represents the fundamental approach to slope stability analysis, evaluating stability by examining the equilibrium of forces acting on a potential failure mass. A slope remains stable when the resisting forces, primarily the shear strength of the soil along the failure surface, exceed the driving forces such as weight and other destabilizing influences. The factor of safety is expressed as the ratio of available shear strength to the shear stress required for equilibrium.

![xslope_results_single.png](lem/images/xslope_results_single.png)

xslope employs the method of slices, a numerical technique that divides the potential failure mass into a series of vertical slices. Rather than analyzing the entire mass as a single unit, each slice is examined individually, with the overall stability determined by summing the forces and moments acting on all slices. This approach allows for the analysis of complex geometries, varying soil conditions, and multiple loading scenarios including distributed surface loads, seismic forces, reinforcement, and tension cracks.

The package implements six distinct limit equilibrium methods, each with different assumptions about inter-slice forces and equilibrium conditions:

- **Ordinary Method of Slices (OMS)** - The simplest approach, satisfying only moment equilibrium with no iteration required. Best for preliminary analysis where speed is prioritized over accuracy.
- **Simplified Janbu Method** - Satisfies horizontal force equilibrium with an empirical correction factor. Suitable for both circular and non-circular failure surfaces.
- **Bishop's Simplified Method** - The most widely used method, satisfying moment and vertical force equilibrium. Provides good accuracy for circular failure surfaces with reasonable computational efficiency.
- **Corps of Engineers Method** - A force equilibrium approach assuming horizontal inter-slice forces. Applicable to any failure surface geometry.
- **Lowe & Karafiath Method** - Similar to Corps Engineers but with inter-slice forces oriented at the average of the slope and failure surface angles. Particularly effective for seismic loading analysis.
- **Spencer's Method** - The most rigorous approach, satisfying complete force and moment equilibrium simultaneously. Generally considered the most accurate method available in xslope.

Beyond single surface analysis, xslope includes automated search algorithms that systematically evaluate thousands of candidate failure surfaces to identify the critical surface with the minimum factor of safety. These search capabilities use adaptive grid refinement to efficiently locate critical surfaces for both circular and non-circular geometries. The package also supports rapid drawdown analysis for dams and levees, as well as reliability analysis using Monte Carlo methods to compute probability of failure based on input parameter uncertainties.

## Finite Element Seepage Analysis

Pore water pressures play a critical role in slope stability, yet they are rarely uniform or static in natural slopes. Traditional approaches such as estimating pore pressures using depth below a piezometric line often fail to capture the complex groundwater flow patterns that develop in heterogeneous soil profiles with varying permeabilities and complex boundary conditions.

xslope addresses this challenge by providing comprehensive finite element seepage analysis capabilities that solve the complete groundwater flow equation throughout the slope domain. The package generates finite element meshes directly from the slope geometry defined in the Excel template, computes spatially varying pore pressure fields that accurately reflect site-specific hydrogeological conditions, and seamlessly integrates these results into slope stability calculations. Both saturated and unsaturated flow problems can be simulated, with the system automatically selecting the appropriate solution algorithm based on boundary conditions.

![johnson_res_solution.png](seepage/images/johnson_res_solution.png)

Beyond slope stability applications, xslope can be used as a standalone 2D groundwater flow analysis package for a wide range of seepage problems. This includes analysis of earth dams, excavations, foundation seepage, flow under sheetpiles, and other civil engineering applications where groundwater flow patterns need to be characterized. The Excel template-based workflow makes it straightforward to set up complex seepage problems with multiple material zones and boundary conditions.

![clay_blanket_solution.png](seepage/images/clay_blanket_solution.png)

## Finite Element Slope Stability Analysis

!!! note "Work in Progress"
    The finite element slope stability analysis module is currently under development.

The finite element method (FEM) provides a powerful numerical approach to slope stability analysis that overcomes fundamental limitations of traditional limit equilibrium methods. While limit equilibrium approaches require engineers to assume a failure surface geometry and then verify equilibrium, FEM allows potential failure mechanisms to emerge naturally through rigorous stress-strain analysis. Rather than imposing kinematic constraints through assumed failure surfaces, FEM solves the complete stress-strain problem throughout the slope domain, capturing complex stress redistribution as soil elements progressively reach failure and allowing failure zones to develop naturally without prior assumptions about their geometry or location.

xslope implements finite element slope stability analysis using an elastic-perfectly plastic constitutive model with Mohr-Coulomb failure criterion. The Shear Strength Reduction Method (SSRM) is employed to determine the factor of safety by systematically reducing soil shear strength parameters until the slope can no longer maintain equilibrium. This approach provides a rigorous and theoretically sound alternative to limit equilibrium methods while maintaining a consistent definition of factor of safety.

The finite element implementation seamlessly integrates with xslope's existing seepage analysis capabilities, enabling coupled stress-seepage analysis for slopes under varying groundwater conditions. Soil reinforcement systems including geotextiles, soil nails, and ground anchors are modeled using embedded truss elements that capture both material tensile strength and pullout resistance. Seismic loading is incorporated through the pseudo-static method, representing earthquake forces as equivalent static body forces throughout the slope mass.

All finite element analysis capabilities use the same Excel input template system as the limit equilibrium and seepage modules, with additional material properties (Young's modulus, Poisson's ratio) and reinforcement parameters (cross-sectional area, elastic modulus, pullout lengths) specified in dedicated columns of the existing tables. This unified input framework ensures consistency across all analysis types while minimizing the learning curve for users already familiar with xslope's limit equilibrium capabilities.
