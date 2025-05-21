# Bishop's Simplified Method

Bishop's Simplified Method is a widely used limit equilibrium technique for analyzing slope stability, especially suitable for **circular slip surfaces**. It improves on the Ordinary Method of Slices by including interslice normal forces and satisfies both **moment** and **vertical force** equilibrium. The key assumptions are:

- Circular slip surface  
- Side forces are **horizontal** (i.e., interslice shear forces are neglected, but normal forces are included)
- Satisfies:
    - Moment equilibrium
    - Vertical force equilibrium
- Does **not** satisfy horizontal force equilibrium

With these assumptions, the forces acting on the slice are as follows:

![bishop_slice.png](images/bishop_slice.png){width=400px}

Where:

>$W$ = weight of the slice  
$\alpha$ = base inclination angle of the slice  
$\Delta \ell$ = length of the base  
$c', \phi'$ = effective cohesion and friction angle  
$u$ = pore water pressure  
$N$ = normal force on the base of the slice  
$S$ = shear force at the base

Summing forces in the vertical direction:

>$\sum F_y = 0$

>$N \cos \alpha + S \sin \alpha - W = 0  \qquad (1)$

>$S = \tau \Delta l = \dfrac{s \Delta \ell}{F}$

>$S = \dfrac{1}{F} \left[c \Delta \ell + (N - u \Delta \ell) \tan \phi' \right]   \qquad (2)$

Substituting (2) into (1):

>$N \cos \alpha + \left( \dfrac{1}{F} \left[ c \Delta \ell + (N - u \Delta \ell) \tan \phi' \right] \right) \sin \alpha - W = 0$

Solving for N:

>$N = \dfrac{W - \dfrac{1}{F} \left[ c \Delta \ell - u \Delta \ell \tan \phi' \right] \sin \alpha}{\cos \alpha + \dfrac{\sin \alpha \tan \phi'}{F}}   \qquad (3)$

From the general equation (based on moment equilibrium):

>$F = \dfrac{\sum (c + \sigma' tan \phi') \Delta \ell}{\sum W sin \alpha}$

>$\sigma' = \dfrac{N}{\Delta \ell} - u$

thus:

>$F = \dfrac{\sum (c \Delta \ell + (N - u \Delta \ell) tan \phi')}{\sum W sin \alpha}   \qquad (4)$

Substituing (3) into (4) and solving for $F$ gives:

>$F = \dfrac{\sum \left[ \dfrac{c \Delta \ell \cos \alpha + (W - u \Delta \ell \cos \alpha) \tan \phi'}{\cos \alpha + \dfrac{\sin \alpha \tan \phi'}{F}} \right]}{\sum W \sin \alpha}   \qquad (5)$

For total stress analysis:

>$F = \dfrac{\sum \left[ \dfrac{c \Delta \ell \cos \alpha + W \tan \phi'}{\cos \alpha + \dfrac{\sin \alpha \tan \phi}{F}} \right]}{\sum W \sin \alpha}$

The factor of safety $F$ appears on both sides of the equation, so it must be solved **iteratively**.

Once $F$ is determined, $N$ can be computed using equation (3) above.

## Complete Formulation

For a complete implementation of Bishop's Simplified Method, we need to consider additional forces acting on the slice. The full set of forces are shown in the following figure:

![bishop_complete.png](images/bishop_complete.png){width=400px}

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

### Vertical Force Equilibrium

To revise the factor of safety equation for Bishop's method to include the $D$, $kw$, $P$, and $T$ forces, we first need to consider how these forces affect the vertical force equilibrium. The vertical force equilibrium equation becomes:

>$N \cos \alpha + S \sin \alpha + P \sin \alpha - W - D \cos \beta = 0  \qquad (6)$

The shear force on the base of the slice remains the same as before:

>$S = \dfrac{1}{F} \left[c \Delta \ell + (N - u \Delta \ell) \tan \phi' \right]   \qquad (7)$

Substituting (7) into (6) and solving for N:

>$N \cos \alpha + \dfrac{1}{F} \left[c \Delta \ell + (N - u \Delta \ell) \tan \phi' \right] \sin \alpha + P \sin \alpha - W - D \cos \beta = 0$

>$N \cos \alpha + \dfrac{1}{F} c \Delta \ell \sin \alpha + \dfrac{1}{F} N \tan \phi' \sin \alpha - \dfrac{1}{F} u \Delta \ell \tan \phi' \sin \alpha + P \sin \alpha - W - D \cos \beta = 0$

>$N \cos \alpha + \dfrac{1}{F} N \tan \phi' \sin \alpha = W + D \cos \beta - P \sin \alpha + \dfrac{1}{F} \left[ c \Delta \ell - u \Delta \ell \tan \phi' \right] \sin \alpha$

>$N \left( \cos \alpha + \dfrac{\sin \alpha \tan \phi'}{F} \right) = W + D \cos \beta - P \sin \alpha + \dfrac{1}{F} \left[ c \Delta \ell - u \Delta \ell \tan \phi' \right] \sin \alpha$

>$N = \dfrac{W + D \cos \beta - P \sin \alpha - \dfrac{1}{F} \left[ c \Delta \ell - u \Delta \ell \tan \phi' \right] \sin \alpha}{\cos \alpha + \dfrac{\sin \alpha \tan \phi'}{F}}   \qquad (8)$

### Moment Equilibrium

The moment equilibrium equation about the center of the slip circle must also be revised to include the moments from the additional forces. The moments are:

>$F = \dfrac{R \sum (S + P)}{R \sum W \sin \alpha + \sum D \cos \beta a_{dx} - \sum D \sin \beta a_{dy} + k\sum W a_s + T a_t}   \qquad (9)$

Where:
- $a_{dx}$ = horizontal distance from center of circle to point $d$
- $a_{dy}$ = vertical distance from center of circle to point $d$
- $a_s$ = vertical distance from center of circle to center of gravity of the slice
- $a_t$ = vertical distance from center of circle to point $c$

### Complete Factor of Safety Equation

Combining (8) and (9), and dividing by R, we get:

>$F = \dfrac{\sum (S + P)}{\sum W \sin \alpha + \frac{1}{R}\sum D \cos \beta a_{dx} - \frac{1}{R}\sum D \sin \beta a_{dy} + \frac{k}{R}\sum W a_s + \frac{1}{R} T a_t}$

Substituting (7) for S:

>$F = \dfrac{\sum \left( \dfrac{1}{F} \left[c \Delta \ell + (N - u \Delta \ell) \tan \phi' \right] + P \right)}{\sum W \sin \alpha + \frac{1}{R}\sum D \cos \beta a_{dx} - \frac{1}{R}\sum D \sin \beta a_{dy} + \frac{k}{R}\sum W a_s + \frac{1}{R} T a_t}$

Substituting (8) for N:

>$F = \dfrac{\sum \left( \dfrac{1}{F} \left[c \Delta \ell + \left( \dfrac{W + D \cos \beta - P \sin \alpha - \dfrac{1}{F} \left[ c \Delta \ell - u \Delta \ell \tan \phi' \right] \sin \alpha}{\cos \alpha + \dfrac{\sin \alpha \tan \phi'}{F}} - u \Delta \ell \right) \tan \phi' \right] + P \right)}{\sum W \sin \alpha + \frac{1}{R}\sum D \cos \beta a_{dx} - \frac{1}{R}\sum D \sin \beta a_{dy} + \frac{k}{R}\sum W a_s + \frac{1}{R} T a_t}$

Simplifying the numerator:

>$F = \dfrac{\sum \left[ \dfrac{c \Delta \ell \cos \alpha + (W + D \cos \beta - P \sin \alpha - u \Delta \ell \cos \alpha) \tan \phi' + P}{\cos \alpha + \dfrac{\sin \alpha \tan \phi'}{F}} \right]}{\sum W \sin \alpha + \frac{1}{R}\sum D \cos \beta a_{dx} - \frac{1}{R}\sum D \sin \beta a_{dy} + \frac{k}{R}\sum W a_s + \frac{1}{R} T a_t}   \qquad (10)$

This is the **complete formulation** for Bishop's Simplified Method. Note that the water force $T$ only applies to the uppermost slice, so for the summation in the numerator, the $T$ value is zero for all other slices.

The factor of safety $F$ appears on both sides of the equation, so it must be solved **iteratively**, just like the basic formulation.

---

## Summary

Assumes **horizontal side forces**<br>
Satisfies **moment** and **vertical force** equilibrium<br>
Applicable to **circular slip surfaces**<br>
Requires **iteration** to solve for $F$<br>
More accurate than OMS, especially for **effective stress analysis** with high pore pressures
