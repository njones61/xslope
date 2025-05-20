# Ordinary Method of Slices (OMS)

The Ordinary Method of Slices (also known as the Swedish Method) is one of the simplest limit equilibrium techniques used for slope stability analysis. The factor of safety to be calculated directly without iteration. The key assumption for the OMS method is that the side forces can be neglected, meaning interslice shear and normal forces cancel out and are excluded from equilibrium considerations. Thus, the forces on the slice are as follows:

>![oms_slice.png](images/oms_slice.png){ width=200px }

Further, only **moment equilibrium** about the center of the slip circle is enforced:

>$\sum M = 0$

Force equilibrium in the horizontal and vertical directions is **not** satisfied.

## Basic Formulation

To solve for the factor of safety, we need to consider the forces acting on the slice:

>$W$ = weight of the slice  <br>
$\alpha$ = inclination of the base of the slice
$\Delta \ell$ = length of the slice base  <br>
$S$ = shear force on the base of the slice = $c \Delta \ell + N \tan \phi$<br>
$c$ = cohesion  <br>
$\phi$ = friction angle <br> 
$N$ = normal force on the base of the slice <br>
$u$ = pore water pressure  <br>

The general expression for the factor of safety $FS$ is:

>$FS = \dfrac{\sum (c \Delta \ell + N \tan \phi)}{\sum W \sin \alpha}  \qquad (1)$

The normal force on the base is:

>$N = W \cos \alpha$

Therefore, we can rewrite (1) as:

>$FS = \dfrac{\sum (c \Delta \ell + W \cos \alpha \tan \phi)}{\sum W \sin \alpha}$

If $\phi = 0$, the expression simplifies to:

>$FS = \dfrac{\sum c \Delta \ell}{\sum W \sin \alpha}$

which is the same equation used in the Swedish method and the log-spiral method under the same assumptions.

To use effective stress parameters $c', \phi'$, and effective normal stress $\sigma'$, the expression becomes:

>$FS = \dfrac{\sum (c' \Delta \ell + (W \cos \alpha - u \Delta \ell) \tan \phi')}{\sum W \sin \alpha}  \qquad (2)$

## Alternate (Preferred) Formulation

This equation for FS can produce unconservative results, including negative normal stresses under high pore pressures. The following formulation is preferred. First, we define the effective weight as follows:

>$W' = W - u b, \quad \text{where} \quad b = \Delta \ell \cos \alpha$
 
>$W' = W - u  \Delta \ell \cos \alpha$

Then:

>$N' = W' \cos \alpha = (W - u  \Delta \ell \cos \alpha) \cos \alpha$

>$N' = W  \cos \alpha - u  \Delta \ell \cos^2 \alpha$

Substituting into (2):

>$F = \dfrac{\sum \left[ c' \Delta \ell + (W \cos \alpha - u \Delta \ell \cos^2 \alpha) \tan \phi' \right]}{\sum W \sin \alpha}   \qquad (3)$

This is the **preferred formulation** for OMS.

## Complete Formulation

For a complete implementation of the Ordinary Method of Slices, we need to consider some additional forces to the slice. The full set of forces acting on the slice are as follows:

>![oms_complete.png](images/oms_complete.png)

Where:

>$D$ = distributed load resultant force <br>
$\beta$ = inclination of the distributed load (perpendicular to slope) <br>
$kW$ = seismic force for pseudo-static seismic analysis <br>
$c.g.$ = center of gravity of the slice <br>
$R_f$ = reinforcement force on base of slice <br>
$water$ = tension crack water force <br>

The rest of the forces are the same as before.

### Distributed Load

The **distributed load** resultant force $D$ is calculated from the distributed load input which is defined as a stress along the top of the slope. It is assumed to act perpendicular to the slope therefore the inclination of the distributed load from a vertical line is equal to the slope angle. The distributed load acts through point $d$ which is often the center of the slice, but it can be offset from the center, depending on how the distributed load is defined. 

To recalculate the factor of safety to account for the distributed load, we start by redefining the **effective weight**. If we consider the distributed load to effectively increase the weigh of the slice, we can represent the effective weight as:

>$W' = (W + D \cos \beta) - u b$

>$W' = W + D \cos \beta - u \Delta \ell cos \alpha   \qquad (4)$

Then the effective normal force is:

>$N' = W' \cos \alpha$

Substituting (4) into this gives:

>$N' = (W + D \cos \beta - u \Delta \ell cos \alpha) \cos \alpha$

>$N' = W \cos \alpha + D \cos \beta \cos \alpha - u \Delta \ell \cos^2 \alpha  \qquad (5)$

We also need to consider the distributed load in the denominator of the factor of safety equation. Before we do this, we recall that the OMS equation is based on moments about the center of the slip circle. The moments in the original method of slices formulation included the weight of the slice and the shear force. The normal force acts through the center of the slice and therefore produces no moment. In the original equation, the limit equilibrium equation is:

>$F = \dfrac{R \sum S}{R \sum W sin \alpha}$

R is the moment arm for both $S$ and $W sin \alpha$. Before, we factored out the R value because it was in both the numerator and denominator. But now we have a distributed load $D$ that is not acting through the center of the slice. The moment arm for the distributed load is the shortest distance from the center of the circle to a line parallel to $D$ that passes through point $d$. We will call this distance $a_d$. Thus, the moment from the distributed load is:

>$M_D = a_d D  \qquad (6)$

Substituting the effective normal force defined by (5) into the numerator and (6) into the denominator of the preferred formulation (3) gives:

>$F = \dfrac{R \sum \left[ c' \Delta \ell + (W \cos \alpha + D \cos \beta \cos \alpha - u \Delta \ell \cos^2 \alpha) \tan \phi'\right]}{R \sum W \sin \alpha + \sum a_d D}   \qquad (7)$

Or if we divide everything by R, we get:

>$F = \dfrac{\sum \left[ c' \Delta \ell + (W \cos \alpha + D \cos \beta \cos \alpha - u \Delta \ell \cos^2 \alpha) \tan \phi'\right]}{\sum W \sin \alpha + \dfrac{1}{R} \sum a_d D}   \qquad (8)$

### Seismic Force

The **seismic force** $kW$ is calculated as a horizontal pseudo-static force acting on the slice through the center of gravity of the slice. It is assumed to act in the direction of sliding. It is equal to the seismic coefficient $k$ multiplied by the weight of the slice $W$. The seismic coefficient is a user-defined input, depending on the seismic conditions of the site.

The seismic force is an additional force acting in the direction of sliding and causing failure, therefore we add it to the denominator of the factor of safety equation. Since the force acts horizontally, the moment arm is equal to the vertical distance between the center of the circle and the center of gravity. If we call this distance $a_s$ then the moment arm from the seismic force is:

>$a_s = y_o - y_{cg}$

where: 

>$y_o$ = y-coordinate of the center of the circle<br>
>$y_{cg}$ = y-coordinate of the center of gravity of the slice<br>

The moment from the seismic force is then:

>$M_s = a_s kW$

Using the same method as before, we can substitute this into the denominator of the factor of safety equation. Inserting the moment from the seismic force into (8), we get:

>$F = \dfrac{\sum \left[ c' \Delta \ell + (W \cos \alpha + D \cos \beta \cos \alpha - u \Delta \ell \cos^2 \alpha) \tan \phi'\right]}{\sum W \sin \alpha + \dfrac{1}{R} \left( \sum a_d D + k\sum a_s W \right)}   \qquad (9)$


### Reinforcement Force

The **reinforcement force** $R_{f}$ is a force on the base of the slice resisting sliding. The reinforcement force is calculated using the reinforcement lines in the input, where the user defines a longitudinal reinforcement force $F_L$ and a transverse reinforcement force $F_T$ at a series of points along each reinforcement line. With the current implementation of slope tools, only the longitudinal reinforcement force $F_L$ is used. We assume that the reinforcement is flexible and therefore bends with the sliding of the failure surface to act parallel to the bottom of the slice. Thus, the reinforcement force is equal to the sum of the interpolated $F_L$ values for the reinforcement lines that intersect the base of the slice.






### Water Force

The **water force** $water$ on the side of the slice is calculated from the tension crack water input only applies if there is both a tension crack, and if the user has selected to fill the crack with water. This force only applies to the side of the upper-most slice and pushes in the direction of sliding. The force is calculated using the hydrostatic water pressure that is zero at the top of the crack (side of slice) and = $\gamma_w d_{tc}$ where $\gamma_w$ = the unit wt of water and $d_{tc}$ is the depth of the tension crack. The resultant force = $\frac{1}{2} \gamma_w d_{tc}^2$ and it acts at point $c$ which is 1/3 of the height of the slice $d_{tc}$.



## Summary

- Applicable only to **circular** slip surfaces.
- **Only moment equilibrium** is satisfied.
- **No iteration** is required.
- **Less accurate** than more complete methods (e.g., Bishop’s or Spencer’s).
- Provides the same solution as the Swedish method when $\phi = 0$.
