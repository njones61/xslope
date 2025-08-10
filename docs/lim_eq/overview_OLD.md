# The Limit Equilibrium Method

## Introduction

The Limit Equilibrium Method represents the fundamental approach to slope stability analysis, evaluating the stability of slopes by examining the equilibrium of forces acting on a potential failure mass. This method operates on the principle that a slope remains stable when the resisting forces, primarily the shear strength of the soil, exceed or equal the driving forces such as weight and other destabilizing influences.

## Theoretical Foundation

### Basic Concept

The limit equilibrium method operates under several key assumptions that form the foundation of the analysis. First, it assumes the soil mass is on the verge of failure, existing at what we call limiting equilibrium. The factor of safety (F) serves as our primary measure, defined as the ratio of available shear strength to required shear strength. The method also assumes that the failure surface is known or can be reasonably assumed, and that the soil behaves as a rigid-plastic material.

The factor of safety is expressed mathematically as:

> $F = \dfrac{\text{Available Shear Strength}}{\text{Required Shear Strength}} = \dfrac{\sum \text{Resisting Forces}}{\sum \text{Driving Forces}}$

### Shear Strength

The shear strength of soil follows the well-established Mohr-Coulomb failure criterion, which provides the theoretical basis for most slope stability analyses:

> $s = c' + \sigma' \tan \phi'$

Where $s$ represents the shear strength, $c'$ is the effective cohesion, $\sigma'$ is the effective normal stress, and $\phi'$ is the effective angle of internal friction. This relationship forms the cornerstone of limit equilibrium analysis and appears in various forms throughout the different methods.

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

- **Simple methods** (OMS, Bishop, Janbu) neglect interslice forces entirely or make simplifying assumptions
- **Advanced methods** (Spencer, Corps Engineers, Lowe-Karafiath) explicitly calculate these forces using various approaches

### Solution Process and Numerical Techniques

The analysis of the entire failure mass proceeds differently depending on the method being used. Each method employs specific numerical approaches to solve for the factor of safety, ranging from direct solutions to iterative techniques.

**Direct Solution Methods**: Simple methods like OMS and Janbu can often be solved directly for the factor of safety, as they involve fewer unknowns and simpler equations. These methods typically use a single equation that can be rearranged to solve directly for F.

**Iterative Solution Methods**: More sophisticated methods like Bishop's, Spencer's, and the force equilibrium methods require iterative solution techniques. These methods typically involve:

1. Assuming an initial value for the factor of safety
2. Solving the equilibrium equations for all slices
3. Checking whether the overall equilibrium is satisfied
4. Adjusting the factor of safety and repeating until convergence

For force equilibrium methods, the analysis proceeds systematically from one end of the failure surface to the other. For a typical left-facing slope, the analysis begins at the leftmost slice and progresses to the right. Each slice introduces new unknowns that must be solved using the equilibrium equations, with the calculated interslice force from one slice serving as input for the next slice.

For moment equilibrium methods like Bishop's, the solution involves iterating on the factor of safety until the calculated normal forces satisfy both moment and vertical force equilibrium for all slices simultaneously.

For complete equilibrium methods like Spencer's, the solution involves iterating on both the factor of safety and the interslice force angle until all equilibrium conditions are satisfied across the entire failure mass.

This progressive approach creates a system of equations that can be solved using various numerical techniques, depending on the complexity of the method and the assumptions made.

## Advanced Considerations

### Pore Water Pressure

![pore_pressure_effects.png](images/pore_pressure_effects.png)

Pore water pressure effects play a critical role in accurate slope stability analysis. These pressures reduce the effective normal stress acting on the failure surface, which in turn decreases the available shear strength. Engineers can include these effects using either total or effective stress parameters, depending on the analysis requirements and available data.

The effect of pore water pressure on slope stability can be analyzed using several approaches:

**Effective Stress Analysis**: The most common approach uses effective stress parameters ($c'$ and $\phi'$) and explicitly includes pore water pressure in the calculations. This method provides the most accurate representation of soil behavior under drained conditions.

**Total Stress Analysis**: For undrained conditions, total stress parameters ($c_u$ and $\phi_u$) can be used, with pore water pressure effects implicitly included in the soil strength parameters.

**Seepage Analysis**: For complex groundwater conditions, separate seepage analysis can be performed to determine pore water pressure distributions, which are then used as input to the stability analysis.

### Seismic Loading

![seismic_loading_analysis.png](images/seismic_loading_analysis.png)

For projects in seismically active regions, pseudo-static analysis can include seismic effects through a horizontal acceleration coefficient $k$. The resulting seismic force $kW$ acts through the center of gravity of each slice, providing a simplified but effective approach to seismic slope stability analysis.

The pseudo-static approach represents a simplified method for including seismic effects in slope stability analysis. More sophisticated approaches include:

**Response Spectrum Analysis**: Uses earthquake response spectra to determine the maximum seismic forces acting on the slope.

**Time History Analysis**: Analyzes the response of the slope to specific earthquake time histories, providing the most detailed representation of seismic effects.

**Newmark Sliding Block Analysis**: Estimates permanent displacement of slopes during earthquakes, providing additional insight into seismic performance.

### Reinforcement

![reinforcement_forces.png](images/reinforcement_forces.png)

Modern slope stabilization often includes geosynthetic or structural reinforcement, which the limit equilibrium methods can accommodate. These forces provide tensile resistance to sliding, can create bending moments that affect stability, and contribute interface shear resistance along the reinforcement elements.

The analysis of reinforced slopes requires consideration of several factors:

**Tensile Resistance**: The tensile strength of reinforcement provides direct resistance to sliding forces, improving the overall factor of safety.

**Interface Shear**: The interaction between reinforcement and soil provides additional shear resistance along the reinforcement surface.

**Bending Effects**: For structural reinforcement, bending moments can affect the stability analysis and must be properly accounted for in the calculations.

**Anchorage**: The anchorage of reinforcement at the top and bottom of the slope must be sufficient to develop the full tensile capacity of the reinforcement.

## Classification of Methods

Limit equilibrium methods can be systematically classified based on which equilibrium equations they satisfy, providing engineers with a framework for selecting the most appropriate analysis technique for their specific project.

### Moment Equilibrium Only

The Ordinary Method of Slices (OMS) represents the simplest approach, satisfying only moment equilibrium while neglecting interslice forces entirely. Bishop's Simplified Method improves upon this by including interslice normal forces and satisfying both moment and vertical force equilibrium, though it still neglects horizontal force equilibrium.

### Force Equilibrium Only

Janbu's Simplified Method takes a different approach, satisfying force equilibrium but not moment equilibrium. This method includes a correction factor to compensate for the missing moment equilibrium, making it suitable for non-circular failure surfaces where moment equilibrium might be less critical.

### Complete Equilibrium

Spencer's Method represents the most sophisticated approach, satisfying both force and moment equilibrium simultaneously. This method can handle both circular and non-circular failure surfaces and provides the highest accuracy, though at the cost of increased computational complexity.

### Force Equilibrium with Interslice Forces

The General Force Equilibrium Method provides a framework for explicit treatment of interslice forces using magnitude and angle parameters. The Corps Engineers Method applies this framework with a constant interslice angle θ based on overall slope geometry, while the Lowe-Karafiath Method uses variable interslice angles based on local slice geometry.

## Key Assumptions and Limitations

The limit equilibrium method relies on several simplifying assumptions that make the analysis tractable while introducing certain limitations. The rigid body assumption treats each slice as moving as a rigid unit, which simplifies the analysis but may not capture the actual deformation behavior of the soil. We also assume that the failure surface is planar and that the same factor of safety applies to all slices, which may not always be realistic for complex soil conditions.

These assumptions lead to several important limitations that engineers must consider. The method provides no information about slope deformation, requiring engineers to specify potential failure surfaces rather than predicting them. The simplified soil behavior model doesn't account for progressive failure, and most implementations assume plane strain conditions, limiting the analysis to two dimensions.

## Advanced Considerations

### Pore Water Pressure

![pore_pressure_effects.png](images/pore_pressure_effects.png)

Pore water pressure effects play a critical role in accurate slope stability analysis. These pressures reduce the effective normal stress acting on the failure surface, which in turn decreases the available shear strength. Engineers can include these effects using either total or effective stress parameters, depending on the analysis requirements and available data.

The effect of pore water pressure on slope stability can be analyzed using several approaches:

**Effective Stress Analysis**: The most common approach uses effective stress parameters ($c'$ and $\phi'$) and explicitly includes pore water pressure in the calculations. This method provides the most accurate representation of soil behavior under drained conditions.

**Total Stress Analysis**: For undrained conditions, total stress parameters ($c_u$ and $\phi_u$) can be used, with pore water pressure effects implicitly included in the soil strength parameters.

**Seepage Analysis**: For complex groundwater conditions, separate seepage analysis can be performed to determine pore water pressure distributions, which are then used as input to the stability analysis.

### Seismic Loading

![seismic_loading_analysis.png](images/seismic_loading_analysis.png)

For projects in seismically active regions, pseudo-static analysis can include seismic effects through a horizontal acceleration coefficient $k$. The resulting seismic force $kW$ acts through the center of gravity of each slice, providing a simplified but effective approach to seismic slope stability analysis.

The pseudo-static approach represents a simplified method for including seismic effects in slope stability analysis. More sophisticated approaches include:

**Response Spectrum Analysis**: Uses earthquake response spectra to determine the maximum seismic forces acting on the slope.

**Time History Analysis**: Analyzes the response of the slope to specific earthquake time histories, providing the most detailed representation of seismic effects.

**Newmark Sliding Block Analysis**: Estimates permanent displacement of slopes during earthquakes, providing additional insight into seismic performance.

### Reinforcement

![reinforcement_forces.png](images/reinforcement_forces.png)

Modern slope stabilization often includes geosynthetic or structural reinforcement, which the limit equilibrium methods can accommodate. These forces provide tensile resistance to sliding, can create bending moments that affect stability, and contribute interface shear resistance along the reinforcement elements.

The analysis of reinforced slopes requires consideration of several factors:

**Tensile Resistance**: The tensile strength of reinforcement provides direct resistance to sliding forces, improving the overall factor of safety.

**Interface Shear**: The interaction between reinforcement and soil provides additional shear resistance along the reinforcement surface.

**Bending Effects**: For structural reinforcement, bending moments can affect the stability analysis and must be properly accounted for in the calculations.

**Anchorage**: The anchorage of reinforcement at the top and bottom of the slope must be sufficient to develop the full tensile capacity of the reinforcement.

### Additional Loading Conditions

![slice_advanced_forces.png](images/slice_advanced_forces.png)

Modern slope stability analysis often includes additional forces that act on individual slices. These forces can significantly affect the stability calculations and must be properly incorporated into the analysis.

**Distributed Loads**: Surface loads such as traffic, construction equipment, or other surcharges can be represented as distributed forces acting on the top of slices. These loads contribute to both the driving forces and the normal forces on the slice base.

**Point Loads**: Concentrated loads from structures, equipment, or other sources can be represented as point forces acting at specific locations on the slope.

**Hydrostatic Pressures**: Water pressures acting on submerged portions of the slope or from tension cracks can significantly affect stability and must be properly included in the analysis.

**Thermal Effects**: For slopes in extreme environments, thermal effects can cause changes in soil properties and must be considered in the analysis.

## Comparison of Implemented Methods

The following table summarizes the key characteristics of each method implemented in XSLOPE, providing engineers with a quick reference for method selection:

| Method | Equilibrium Satisfied | Interslice Forces | Circular Surfaces | Non-Circular Surfaces | Iteration Required | Accuracy | Use Case |
|--------|----------------------|-------------------|-------------------|----------------------|-------------------|----------|----------|
| **OMS** | Moment only | Neglected | ✓ | ✓ | No | Low | Quick screening |
| **Bishop** | Moment + Vertical | Horizontal only | ✓ | ✗ | Yes | Medium | Circular surfaces |
| **Janbu** | Force only | Neglected | ✓ | ✓ | No | Low-Medium | General analysis |
| **Spencer** | Force + Moment | Parallel forces | ✓ | ✓ | Yes | High | Comprehensive analysis |
| **General Force Equilibrium** | Force only | Explicit θ | ✓ | ✓ | Yes | Medium-High | General force equilibrium |
| **Corps Engineers** | Force only | Constant θ | ✓ | ✓ | Yes | Medium-High | USACE standard |
| **Lowe-Karafiath** | Force only | Variable θ | ✓ | ✓ | Yes | Medium-High | Variable interslice angles |
| **Rapid Drawdown** | Varies by method | Varies by method | ✓ | ✓ | Varies | High | Dam/levee analysis |

### Equations and Unknowns

Each method employs different mathematical approaches to solve for the factor of safety. The following table provides insight into the computational complexity and solution methods:

| Method | Primary Equation | Unknowns | Solution Method |
|--------|------------------|-----------|-----------------|
| **OMS** | $F = \dfrac{\sum (c' \Delta \ell + N' \tan \phi')}{\sum W \sin \alpha}$ | $F$ | Direct |
| **Bishop** | $F = \dfrac{\sum \left[ c' \Delta \ell + (W - u \Delta \ell \cos \alpha) \tan \phi' \right]}{\sum W \sin \alpha \left( \cos \alpha + \dfrac{\sin \alpha \tan \phi'}{F} \right)}$ | $F$ | Iterative |
| **Janbu** | $F = f_o \cdot \dfrac{\sum (c' \Delta \ell + N' \tan \phi')}{\sum W \sin \alpha}$ | $F$ | Direct + Correction |
| **Spencer** | Complex system of equations | $F$, $\theta$ | Iterative |
| **General Force Equilibrium** | Force equilibrium with explicit θ | $F$ | Iterative |
| **Corps Engineers** | Force equilibrium with constant θ | $F$ | Iterative |
| **Lowe-Karafiath** | Force equilibrium with variable θ | $F$ | Iterative |
| **Rapid Drawdown** | Multi-stage analysis | $F_{stage1}$, $F_{stage2}$, $F_{stage3}$ | Multi-stage |

## General Assessment of XSLOPE Methods

### Method Strengths and Limitations

Each method in XSLOPE offers unique advantages and limitations that engineers must consider when selecting an analysis approach. Understanding these characteristics helps ensure that the chosen method provides appropriate accuracy for the project requirements while remaining computationally efficient.

#### OMS (Ordinary Method of Slices)

The Ordinary Method of Slices provides the simplest approach to slope stability analysis, offering speed and simplicity that make it ideal for preliminary analysis and quick screening. Since it requires no iteration, results can be obtained rapidly, and the method works on any failure surface geometry. However, these advantages come with significant limitations: the method only satisfies moment equilibrium, can produce unconservative results, and completely neglects interslice forces. This makes OMS most suitable for preliminary analysis and simple geometries where quick results are more important than high accuracy.

#### Bishop's Simplified Method

Bishop's method represents a significant improvement over OMS by satisfying both moment and vertical force equilibrium. This approach provides more accurate results while maintaining reasonable computational efficiency. The method works particularly well for circular failure surfaces, where its assumptions about horizontal interslice forces are most valid. However, the requirement for circular surfaces limits its applicability, and the iterative solution process increases computational time. Bishop's method finds its niche in projects with circular failure surfaces and moderate accuracy requirements.

#### Janbu's Simplified Method

Janbu's method offers a different approach by satisfying force equilibrium rather than moment equilibrium. This makes it suitable for non-circular failure surfaces where moment equilibrium might be less critical. The method includes a correction factor to compensate for the missing moment equilibrium, though this factor is empirical and may not always provide accurate compensation. Janbu's method works well for general analysis and non-circular surfaces, providing conservative estimates that may be appropriate for preliminary design.

#### Spencer's Method

Spencer's method represents the most sophisticated approach available in XSLOPE, satisfying both force and moment equilibrium simultaneously. This comprehensive approach provides the highest accuracy and can handle both circular and non-circular failure surfaces. However, this accuracy comes at a cost: the method is the most complex to implement, requires iterative solution, and is computationally intensive. Spencer's method is best reserved for final design work, high accuracy requirements, and complex geometries where the additional computational effort is justified by the need for precise results.

#### General Force Equilibrium Method

The General Force Equilibrium Method provides a framework for explicit treatment of interslice forces using magnitude and angle parameters. This method works on any failure surface and offers good accuracy when interslice forces are important. However, it only satisfies force equilibrium, requires iterative solution, and depends on assumptions about interslice force angles. This method is most valuable when explicit interslice force treatment is needed and when complex geometries require more sophisticated analysis than simple methods can provide.

#### Corps Engineers Method

The Corps Engineers method applies the force equilibrium framework with a constant interslice angle assumption, providing good accuracy while maintaining reasonable computational efficiency. This approach works well for complex geometries and provides a good balance between accuracy and speed. The method is particularly suitable for USACE projects and general engineering practice where consistency and reliability are important. However, it still only satisfies force equilibrium and requires iterative solution.

#### Lowe-Karafiath Method

The Lowe-Karafiath method represents the most sophisticated force equilibrium approach, using variable interslice angles based on local slice geometry. This provides high accuracy for complex geometries and irregular failure surfaces where local variations significantly affect stability. However, this accuracy comes with increased complexity compared to the Corps Engineers method, and the method still only satisfies force equilibrium. Lowe-Karafiath is most valuable for projects requiring high accuracy on irregular failure surfaces.

#### Rapid Drawdown Analysis

Rapid drawdown analysis represents a specialized application that can use any of the other methods as its foundation. This approach is specifically designed for dam and levee analysis where water level changes create unique stability challenges. The method accounts for undrained conditions and provides multi-stage analysis that captures the complex behavior of soils during rapid water level changes. However, this specialization limits its applicability to specific scenarios, and the multi-stage approach increases computational complexity.

## Recommendations for Use

### When to Use Each Method

Selecting the appropriate method for a slope stability analysis requires careful consideration of project requirements, accuracy needs, and computational resources. The following guidelines help engineers make informed decisions about method selection.

For quick preliminary analysis and simple geometries, the Ordinary Method of Slices provides rapid results that can help identify potential stability issues early in the design process. When working with circular failure surfaces and moderate accuracy requirements, Bishop's Simplified Method offers a good balance of accuracy and efficiency. Janbu's Simplified Method works well for general analysis and non-circular surfaces where conservative estimates are acceptable.

When high accuracy requirements and complex geometries demand comprehensive analysis, Spencer's Method provides the most rigorous approach available. For projects requiring explicit interslice force treatment, the General Force Equilibrium Method offers a solid foundation. The Corps Engineers method is particularly well-suited for USACE projects and general engineering practice where constant interslice angle assumptions are appropriate. When variable interslice angles are needed for complex geometries, the Lowe-Karafiath method provides the necessary sophistication.

For dam and levee analysis involving water level changes, the Rapid Drawdown analysis provides specialized capabilities that address the unique challenges of these structures.

### Best Practices

Successful slope stability analysis requires more than just selecting the right method. Engineers should always verify results using multiple methods when possible, particularly for critical projects where accuracy is paramount. Soil conditions should be carefully considered when selecting analysis methods, as different soil types may respond better to different approaches.

Pore water pressure effects should always be included for accurate results, as these pressures often represent the most significant factor affecting slope stability. Safety factors should be selected based on the accuracy of the chosen method, with more sophisticated methods potentially allowing for lower safety factors due to their improved accuracy.

Finally, field observations should be used to validate analytical results whenever possible, providing a reality check that can identify issues with the analysis approach or input parameters.

## References

For detailed derivations and specific implementations, refer to the individual method documentation:
- [Ordinary Method of Slices (OMS)](oms.md)
- [Bishop's Simplified Method](bishop.md)  
- [Janbu's Simplified Method](janbu.md)
- [Spencer's Method](spencer.md)
- [General Force Equilibrium Method](force_eq.md)
- [Rapid Drawdown Analysis](rapid.md)

**Note**: The Corps Engineers and Lowe-Karafiath methods are implemented in XSLOPE but do not yet have dedicated documentation pages. These methods use force equilibrium with different approaches to interslice force angles.

Additional theoretical background can be found in the reference documents in the `docs/ref_docs_slope/` directory, particularly:
- USACE Slope Stability Manual
- Duncan (1996) State-of-the-Art Review
- Original Spencer (1967) paper

## Summary

XSLOPE provides a comprehensive suite of limit equilibrium methods for slope stability analysis, ranging from simple screening methods to sophisticated analysis techniques. The implementation covers the full spectrum of equilibrium approaches, from basic moment equilibrium methods to complete force and moment equilibrium solutions.

The simple methods (OMS, Bishop, and Janbu) serve well for quick analysis and preliminary design work, providing rapid results that help identify potential issues early in the design process. The general force equilibrium method offers explicit interslice force treatment for cases where these forces significantly affect stability. Advanced force equilibrium methods (Corps Engineers and Lowe-Karafiath) handle complex geometries with different approaches to interslice force modeling. Complete equilibrium methods like Spencer provide the highest accuracy for final design work, while specialized analysis capabilities like Rapid Drawdown address specific engineering scenarios.

The choice of method should be based on the specific requirements of the analysis, including accuracy needs, computational resources, and project specifications. For critical projects, it is recommended to verify results using multiple methods to ensure reliability and accuracy of the analysis. This comprehensive approach ensures that engineers can select the most appropriate method for their specific needs while maintaining confidence in the results.

