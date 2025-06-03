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

Now we solve for the normal force $N$. First, we rearrange the equation:

>$N \cos \alpha + \Bigl[\dfrac{c \Delta \ell}{F} \sin \alpha + \dfrac{(N - u \Delta \ell) \tan \phi'}{F} \sin \alpha\Bigr] - W \;=\; 0$

Next we isolate all terms involving $N$:

>$N \cos \alpha + \dfrac{N \tan \phi'}{F} \sin \alpha + \left(-u \Delta \ell \tan \phi' \right)\dfrac{\sin \alpha}{F} + 
\dfrac{c\,\Delta \ell}{F} \sin \alpha - W = 0$

Now, collect the $N$ terms on the left side:

>$N \cos \alpha + \dfrac{N \tan \phi'}{F} \sin \alpha = W - \dfrac{c\,\Delta \ell}{F} \sin \alpha - \left(-u \Delta \ell \tan \phi' \right)\dfrac{\sin \alpha}{F}$

>$N\Bigl(\cos \alpha + \dfrac{\sin \alpha\,\tan \phi'}{F}\Bigr)  =
W - \dfrac{c\,\Delta \ell}{F}\,\sin \alpha + \dfrac{u\,\Delta \ell\,\tan \phi'}{F}\,\sin \alpha$

Finally, we can solve for $N$:

>$N = \dfrac{W - \dfrac{1}{F} \left[ c \Delta \ell - u \Delta \ell \tan \phi' \right] \sin \alpha}{\cos \alpha + \dfrac{\sin \alpha \tan \phi'}{F}}   \qquad (3)$

From the general equation (based on moment equilibrium):

>$F = \dfrac{\sum (c + \sigma' tan \phi') \Delta \ell}{\sum W sin \alpha}$

>$\sigma' = \dfrac{N}{\Delta \ell} - u$

thus:

>$F = \dfrac{\sum (c \Delta \ell + (N - u \Delta \ell) tan \phi')}{\sum W sin \alpha}   \qquad (4)$

Next, we substitute (3) into (4) and solve for $F$:

>$F = \dfrac{\sum \left[ \dfrac{c \Delta \ell \cos \alpha + (W - u \Delta \ell \cos \alpha) \tan \phi'}{\cos \alpha + \dfrac{\sin \alpha \tan \phi'}{F}} \right]}{\sum W \sin \alpha}   \qquad (5)$

## BEGIN INSERT


Step 1: Substitute <span>N</span> from (3) into (4):

<span>F = \frac{\sum \left[ c \Delta \ell + \left( \frac{W - \frac{1}{F} [c \Delta \ell - u \Delta \ell \tan \phi'] \sin \alpha}{\cos \alpha + \frac{\sin \alpha \tan \phi'}{F}} - u \Delta \ell \right) \tan \phi' \right]}{\sum W \sin \alpha}</span>

Step 2: Expand the numerator:

<span>= \frac{\sum \left[ c \Delta \ell + \frac{(W - \frac{1}{F} [c \Delta \ell - u \Delta \ell \tan \phi'] \sin \alpha) \tan \phi'}{\cos \alpha + \frac{\sin \alpha \tan \phi'}{F}} - u \Delta \ell \tan \phi' \right]}{\sum W \sin \alpha}</span>

Step 3: Combine like terms:

Notice <span>-u \Delta \ell \tan \phi'</span> appears both inside and outside the fraction. Combine them:

<span>= \frac{\sum \left[ \frac{c \Delta \ell \cos \alpha + (W - u \Delta \ell \cos \alpha) \tan \phi'}{\cos \alpha + \frac{\sin \alpha \tan \phi'}{F}} \right]}{\sum W \sin \alpha}</span>

Step 4: Final form (Equation 5):

<span>F = \frac{\sum \left[ \frac{c \Delta \ell \cos \alpha + (W - u \Delta \ell \cos \alpha) \tan \phi'}{\cos \alpha + \frac{\sin \alpha \tan \phi'}{F}} \right]}{\sum W \sin \alpha}</span>

## END INSERT

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

Each of these forces is described in detail in the [Ordinary Method of Slices (OMS)](oms.md) section. The forces $D$, $kW$, $P$, and $T$ are included in the Bishop's method factor of safety equation as follows:

### Vertical Force Equilibrium

First we first need to consider how these forces affect the vertical force equilibrium. The vertical force equilibrium equation becomes:

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
