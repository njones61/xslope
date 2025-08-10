# The Limit Equilibrium Method

## Introduction

The Limit Equilibrium Method represents the fundamental approach to slope stability analysis, evaluating the stability of slopes by examining the equilibrium of forces acting on a potential failure mass. This method operates on the principle that a slope remains stable when the resisting forces, primarily the shear strength of the soil, exceed or equal the driving forces such as weight and other destabilizing influences.

The limit equilibrium method operates under several key assumptions that form the foundation of the analysis. First, it assumes the soil mass is on the verge of failure, existing at what we call limiting equilibrium. The factor of safety (F) serves as our primary measure, defined as the ratio of available shear strength to required shear strength. The method also assumes that the failure surface is known or can be reasonably assumed, and that the soil behaves as a rigid-plastic material.

The factor of safety is expressed mathematically as:

> $F = \dfrac{\text{Available Shear Strength}}{\text{Required Shear Strength}} = \dfrac{\sum \text{Resisting Forces}}{\sum \text{Driving Forces}}$

### Shear Strength

The shear strength of soil follows the well-established Mohr-Coulomb failure criterion, which provides the theoretical basis for most slope stability analyses:

> $s = c' + \sigma' \tan \phi'$

Where $s$ represents the shear strength, $c'$ is the effective cohesion, $\sigma'$ is the effective normal stress, and $\phi'$ is the effective angle of internal friction. This relationship forms the cornerstone of limit equilibrium analysis and appears in various forms throughout the different methods.

### Pore Water Pressure

![pore_pressure_effects.png](images/pore_pressure_effects.png)

Pore water pressure effects play a critical role in accurate slope stability analysis. These pressures reduce the effective normal stress acting on the failure surface, which in turn decreases the available shear strength. Engineers can include these effects using either total or effective stress parameters, depending on the analysis requirements and available data.

The effect of pore water pressure on slope stability can be analyzed using several approaches:

**Effective Stress Analysis**: The most common approach uses effective stress parameters ($c'$ and $\phi'$) and explicitly includes pore water pressure in the calculations. This method provides the most accurate representation of soil behavior under drained conditions.

**Total Stress Analysis**: For undrained conditions, total stress parameters ($c_u$ and $\phi_u$) can be used, with pore water pressure effects implicitly included in the soil strength parameters.

**Seepage Analysis**: For complex groundwater conditions, separate seepage analysis can be performed to determine pore water pressure distributions, which are then used as input to the stability analysis.

## Method of Slices

### General Approach

The Method of Slices represents a numerical technique that divides the potential failure mass into a series of vertical slices for analysis. Rather than analyzing the entire mass as a single unit, each slice is examined individually, with the overall stability determined by summing the forces and moments acting on all slices. This approach allows us to handle complex geometries and varying soil conditions that would be impossible to analyze using simpler methods.

The fundamental concept behind the method of slices is that the failure mass can be discretized into a finite number of vertical slices, each with its own geometry and soil properties. By analyzing the equilibrium of each slice individually and then combining the results, we can determine the overall stability of the slope. This discretization approach makes it possible to handle irregular failure surfaces, varying soil conditions, and complex loading scenarios.

### Basic Slice Geometry

![slice_basic.png](images/slice_basic.png)

For each slice, we consider the fundamental forces that govern its behavior. The weight of the slice ($W$) acts through its center of gravity, while the normal force ($N$) and shear force ($S$) act on the base of the slice. The inclination angle of the slice base ($\alpha$) and the length of the slice base ($\Delta \ell$) define the geometry that determines how these forces interact.

The geometry of each slice is defined by several key parameters. The slice width ($\Delta x$) represents the horizontal distance between slice boundaries, while the slice height varies from the ground surface to the failure surface. The base length ($\Delta \ell$) is calculated as $\Delta x \sec \alpha$, where $\alpha$ is the inclination angle of the slice base. This geometric relationship is crucial for accurate force calculations.

### Mathematical Formulation

The mathematical analysis of each slice begins with the fundamental equilibrium equations. For a typical slice, we must satisfy both force and moment equilibrium conditions. The force equilibrium equations can be written as:

**Horizontal Force Equilibrium:**
> $\sum F_x = 0: S \cos \alpha - N \sin \alpha + E_{i+1} - E_i + kW = 0$

**Vertical Force Equilibrium:**
> $\sum F_y = 0: S \sin \alpha + N \cos \alpha - W + X_{i+1} - X_i = 0$

Where $kW$ represents the seismic force (if applicable) and the interslice forces $E_i$, $E_{i+1}$, $X_i$, and $X_{i+1}$ represent the horizontal and vertical components of the forces acting on the left and right sides of the slice.

The moment equilibrium equation, typically taken about the center of rotation for circular surfaces, can be expressed as:

> $\sum M_o = 0: W \cdot x_w - N \cdot x_n - S \cdot y_s + \text{interslice moments} = 0$

Where $x_w$, $x_n$, and $y_s$ represent the moment arms of the respective forces.

### Shear Strength and Mobilization

The shear force acting on the base of each slice is directly related to the soil's shear strength and the factor of safety. The mobilized shear strength represents the portion of the available strength that is actually resisting sliding:

> $\tau_m = \dfrac{c' + \sigma' \tan \phi'}{F}$

Where $\tau_m$ is the mobilized shear stress, $c'$ is the effective cohesion, $\sigma'$ is the effective normal stress, $\phi'$ is the effective friction angle, and $F$ is the factor of safety.

The total shear force on the slice base is then:

> $S = \tau_m \Delta \ell = \dfrac{c' \Delta \ell + N' \tan \phi'}{F}$

Where $N' = N - u \Delta \ell$ is the effective normal force, accounting for pore water pressure effects.

The general factor of safety equation for the method of slices captures the essence of this approach:

> $F = \dfrac{\sum \left[c' \Delta \ell + (N - u \Delta \ell) \tan \phi'\right]}{\sum W \sin \alpha}$

### Interslice Forces

![slice_interslice_forces.png](images/slice_interslice_forces.png)

The treatment of interslice forces represents one of the most critical aspects of the method of slices. These forces arise from the interaction between adjacent slices and must be properly accounted for to maintain equilibrium. The magnitude and direction of these forces depend on the specific method being used and the assumptions made about their behavior.

In the most general case, the interslice forces can be represented as:

> $E_i = Z_i \cos \theta_i$
> $X_i = Z_i \sin \theta_i$

Where $Z_i$ is the magnitude of the interslice force and $\theta_i$ is its inclination angle. Different methods make different assumptions about these parameters:

## Advanced Loading Conditions

![slice_advanced_forces.png](images/slice_advanced_forces.png)

Modern slope stability analysis often includes additional forces that act on individual slices. These forces can significantly affect the stability calculations and must be properly incorporated into the analysis.

**Distributed Loads**: Surface loads such as traffic, construction equipment, or other surcharges can be represented as distributed forces acting on the top of slices. These loads contribute to both the driving forces and the normal forces on the slice base.

**Seismic Loads**
![seismic_loading_analysis.png](images/seismic_loading_analysis.png)

For projects in seismically active regions, pseudo-static analysis can include seismic effects through a horizontal acceleration coefficient $k$. The resulting seismic force $kW$ acts through the center of gravity of each slice, providing a simplified but effective approach to seismic slope stability analysis.

**Reinforcement**

![reinforcement_forces.png](images/reinforcement_forces.png)

Modern slope stabilization often includes geosynthetic or structural reinforcement, which the limit equilibrium methods can accommodate. These forces provide tensile resistance to sliding, can create bending moments that affect stability, and contribute interface shear resistance along the reinforcement elements. For the limit equilibrium methods in XSLOPE, the reinforcement is assumed to be flexible and therefore acts in a direction parallel to the base of the slice in a direciton that resists shear. More comprehensize treatment of reinforcement is included in the [finite element method](../fem/overview.md).

**Tension Cracks**

In the upper part of the slope, the cohesion of a soil can be greater than the driving forces. Since soils can generally not withstand tension, this is unconservative. To address this problem, a tension crack can be added to a user-specified depth and the tension crack forms the upper boundary of the slices and not cohesive forces are allowed on the slice (crack) boundary. It is also possible to assume that the crack fills with water, providing a small force driving failure as an extra measure of conservative analysis.

## Limit Equilibrium Methods Supported in XSLOPE

The following limit equilibrium methods are supported in XSLOPE. Each method has a page describing in detail the equations used, solution technique, etc. 

**Ordinary Method of Slices (OMS):** The Ordinary Method of Slices provides the simplest approach to slope stability analysis, offering speed and simplicity that make it ideal for preliminary analysis and quick screening. Since it requires no iteration, results can be obtained rapidly. However, these advantages come with significant limitations: the method only satisfies moment equilibrium, can produce unconservative results, and completely neglects interslice forces. This makes OMS most suitable for preliminary analysis and simple geometries where quick results are more important than high accuracy. [Documentation](oms.md)

**Simplified Janbu Method:** Janbu's method offers a different approach by satisfying force equilibrium rather than moment equilibrium. This makes it suitable for circular or non-circular failure surfaces where moment equilibrium might be less critical. The method includes a correction factor to compensate for the missing moment equilibrium, though this factor is empirical and may not always provide accurate compensation. Janbu's method works well for general analysis and non-circular surfaces, providing conservative estimates that may be appropriate for preliminary design. [Documentation](janbu.md)

**Bishop's Simplified Procedure:** Bishop's method represents a significant improvement over OMS and Janbu by satisfying both moment and vertical force equilibrium. This approach provides more accurate results while maintaining reasonable computational efficiency. However, the requirement for circular surfaces limits its applicability, and the iterative solution process increases computational time. Bishop's method finds its niche in projects with circular failure surfaces and moderate accuracy requirements.
 [Documentation](janbu.md)

**Corps of Engineer's Method:** Force equilibrium method that assumes all side forces are horizontal. Satisfies horizontal and vertical force equilibrim but not moment equilibrium. Can be used with both circular and non-circular failure surfaces. Requires and iterative solution. [Documentation](force_eq.md)

**Force Equilibrium Methods:** The General Force Equilibrium Method provides a framework for explicit treatment of interslice forces using magnitude and angle parameters. This method works on any failure surface and offers good accuracy when interslice forces are important. However, it only satisfies force equilibrium, requires iterative solution, and depends on assumptions about interslice force angles. This method is most valuable when explicit interslice force treatment is needed and when complex geometries require more sophisticated analysis than simple methods can provide. [Documentation](force_eq.md)

The **Corps Engineers method** applies the force equilibrium framework with a constant interslice angle assumption, providing good accuracy while maintaining reasonable computational efficiency. This approach works well for complex geometries and provides a good balance between accuracy and speed. The method is particularly suitable for USACE projects and general engineering practice where consistency and reliability are important. However, it still only satisfies force equilibrium and requires iterative solution.

The **Lowe-Karafiath method** represents the most sophisticated force equilibrium approach, using variable interslice angles based on local slice geometry. This provides high accuracy for complex geometries and irregular failure surfaces where local variations significantly affect stability. However, this accuracy comes with increased complexity compared to the Corps Engineers method, and the method still only satisfies force equilibrium. Lowe-Karafiath is most valuable for projects requiring high accuracy on irregular failure surfaces.

**Spencer's Methos:** Spencer's method represents the most sophisticated approach available in XSLOPE, satisfying both force and moment equilibrium simultaneously. This comprehensive approach provides the highest accuracy and can handle both circular and non-circular failure surfaces. However, this accuracy comes at a cost: the method is the most complex to implement, requires iterative solution, and is computationally intensive. Spencer's method is generally considered the best and most accurate of the methods supported in XSLOPE. [Documentation](spencer.md)

## Automated Search for the Critical Factor of Safety 

For each of the solution methods, the factor of safety can be computed either for a single failure surface or XSLOPE can perform an exhaustive search where a large number of candidate failure surfaces are considered until the surface with the minimum or critical factor of safety is found. This process is described in more detail in the [search documentation](search.md).

## Rapid Drawdown Analysis

Rapid drawdown analysis represents a specialized application that can use any of the other methods as its foundation. This approach is specifically designed for dam and levee analysis where water level changes create unique stability challenges. The method accounts for undrained conditions and provides multi-stage analysis that captures the complex behavior of soils during rapid water level changes. However, this specialization limits its applicability to specific scenarios, and the multi-stage approach increases computational complexity. In XSLOPE, rapid drawdown analysis can be peformed with any of the supported limit equilibrium methods. [Documentation](rapid.md)

## Reliability Analysis

XSLOPE includes an option to perform a reliability analysis with any of the supported limit equilibrium methods. Rather than finding a single factor of safety, selected inputs are perturbed and the critical factor of safety is computed for each combination of inputs allowing the computation of a probability of failure. [Documentation](reliability.md)

