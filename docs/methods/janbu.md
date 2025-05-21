#Janbu Method

Janbu's Simplified Method is one of the earliest and most commonly used slope stability analysis techniques based on limit equilibrium principles. It provides an approximate factor of safety (FS) for a mass of soil sliding along a specified failure surface by assuming force equilibrium in the horizontal direction only. The method was developed by Norwegian engineer **Nils Janbu** in the 1950s to simplify slope stability analysis for circular and non-circular failure surfaces. It builds on the idea of dividing the soil mass into vertical slices and analyzing each slice individually while applying force balance in the horizontal direction. This method is especially favored for hand calculations and preliminary screening due to its simplicity. However, it is also inherently conservative because it does not enforce moment equilibrium or consider interslice side forces directly.

The factor of safety computed by Janbu's simplified method typically **underestimates** the true FS and is therefore considered conservative. A correction factor was later introduced to compensate for the oversimplifications inherent in the method. 

Key Assumptions:

- **Force equilibrium only**: Satisfies equilibrium in the **horizontal direction**; vertical and moment equilibrium are not enforced.
- **No side forces**: Interslice normal and shear forces are **ignored** (assumed to cancel out overall).
- **Base forces only**: Slice weight, pore pressure, cohesion, and base friction are considered.
- **Correction factor**: An empirical correction factor \( f_o \) is applied to account for the neglected side forces and moment imbalance.

##  Equations

The force equilibrium is applied in the horizontal direction for the entire sliding mass:

>>$\sum S_i = \sum T_i$

Where:

>>$S_i$ = Shear resistance at base of slice $i$<br>
>>$T_i$ = Driving force from weight of slice $i$

The resisting force for each slice is:

>>$S_i = c_i \ell_i + \dfrac{(N_i - u_i \ell_i) \tan\phi_i}{F}$

And:

>>$N_i = W_i \cos\alpha_i$<br>
>>$T_i = W_i \sin\alpha_i$

The factor of safety $F$ is solved iteratively from:

>>$F = \dfrac{\sum \left[ c_i \ell_i + (N_i - u_i \ell_i) \tan\phi_i / F \right]}{\sum W_i \sin\alpha_i}$

This is solved using successive substitution where the prior factor of safety is used to calculate the new factor of safety. The process is repeated until the change in factor of safety is negligible.

## Correction Factor $f_o$

To account for the neglect of interslice forces and moment equilibrium, Janbu proposed a correction factor:

The correction factor $f_o$ is based on the following relationship:

>>![janbu_correction_factor.jpg](images/janbu_correction_factor.jpg){ width=50% }

The d/L ratio is the ratio of the distance from the center of the failure surface to the point of interest (d) and the length of the failure surface (L). The correction factor is used to account for the fact that the Janbu method does not satisfy moment equilibrium. The correction factor is a function of the d/L ratio and is calculated using the following equation:

>>$f_o = 1 + b_1  * \left[\dfrac{d}{L} - 1.4 * \left(\dfrac{d}{L}\right)^2\right]$

The $b_1$ value is a function of the soils in the slope and is found as follows:


| Soil Type                           | $b_1$ |
|-------------------------------------|-------|
| c-only soil (undrained, $\phi$ = 0) | 0.67  |
| c-$\phi$ soil                       | 0.5   |
| $\phi$-only soil (no cohesion)      | 0.31  |

This correction attempts to mimic the effects of moment balance and interslice forces without modeling them directly.

## Complete Formulation with Extra Forces

For a complete implementation of Janbu's Simplified Method, we need to consider additional forces acting on the slice. The full set of forces are shown in the following figure:

![janbu_complete.png](images/janbu_complete.png){width=400px}

Where:

>$D$ = distributed load resultant force <br>
$\beta$ = inclination of the distributed load (perpendicular to slope) <br>
$kW$ = seismic force for pseudo-static seismic analysis <br>
$c.g.$ = center of gravity of the slice <br>
$P$ = reinforcement force on base of slice <br>
$T$ = tension crack water force <br>

The **distributed load** resultant force $D$ is calculated from the distributed load input which is defined as a stress along the top of the slope. It is assumed to act perpendicular to the slope, therefore the inclination of the distributed load from a vertical line is equal to the slope angle, $\beta$. The distributed load acts through point $d$ which is often the center of the slice, but it can be offset from the center, depending on how the distributed load is defined.

The **seismic force** $kW$ is calculated as a horizontal pseudo-static force acting on the slice through the center of gravity of the slice. It is assumed to act in the direction of sliding. It is equal to the seismic coefficient $k$ multiplied by the weight of the slice $W$.

The **reinforcement force** $P$ is a force on the base of the slice resisting sliding. The reinforcement force is calculated using the reinforcement lines in the input, where the user defines a longitudinal reinforcement force $F_L$ and a transverse reinforcement force $F_T$ at a series of points along each reinforcement line. With the current implementation of slope tools, only the longitudinal reinforcement force $F_L$ is used. We assume that the reinforcement is flexible and therefore bends with the sliding of the failure surface to act parallel to the bottom of the slice.

The **water force** $T$ on the side of the slice is calculated from the tension crack water input only applies if there is both a tension crack, and if the user has selected to fill the crack with water. This force only applies to the side of the uppermost slice and pushes in the direction of sliding. The force is calculated using the hydrostatic water pressure that is zero at the top of the crack (side of slice) and = $\gamma_w d_{tc}$ where $\gamma_w$ = the unit wt of water and $d_{tc}$ is the depth of the tension crack. The resultant force = $\frac{1}{2} \gamma_w d_{tc}^2$ and it acts at point $c$ which is 1/3 of the height of the slice $d_{tc}$.

### Horizontal Force Equilibrium

To revise the factor of safety equation for Janbu's method to include the $D$, $kw$, $P$, and $T$ forces, we first need to consider how these forces affect the horizontal force equilibrium. The horizontal force equilibrium equation becomes:

>$\sum F_x = 0 \Rightarrow S \cos \alpha - N \sin \alpha + kW + D \sin \beta + T = 0  \qquad (1)$

The shear force on the base of the slice remains the same as before:

>$S = \dfrac{1}{F} \left[c \Delta \ell + (N - u \Delta \ell) \tan \phi' \right]   \qquad (2)$

Substituting (2) into (1) and solving for N:

>$N = \dfrac{W \cos \alpha + D \cos \beta - P \sin \alpha - \dfrac{1}{F} \left[ c \Delta \ell - u \Delta \ell \tan \phi' \right] \sin \alpha}{\cos \alpha + \dfrac{\sin \alpha \tan \phi'}{F}}   \qquad (3)$

### Complete Factor of Safety Equation

The complete factor of safety equation for Janbu's method becomes:

>$F = \dfrac{\sum \left[ c \Delta \ell + (W \cos \alpha + D \cos \beta - P \sin \alpha - u \Delta \ell) \tan \phi' \right]}{\sum W \sin \alpha + \sum D \sin \beta + k\sum W + T}   \qquad (4)$

This equation must be solved iteratively since $F$ appears in the normal force calculation. The water force $T$ only applies to the uppermost slice, so for the summation in the numerator, the $T$ value is zero for all other slices.

Note that this equation still requires the correction factor $f_o$ to account for the neglect of interslice forces and moment equilibrium.

---

## Summary

Advantages:

- Fast and easy to compute
- Conservative (safe) results
- Effective for preliminary design and parametric studies
- Can be applied to circular and non-circular failure surfaces

Limitations:

- May significantly underestimate FS for complex geometries
- Not suitable for highly irregular slopes or layered soils with complex interactions
- Not a complete equilibrium method (no moment or interslice force balance)

