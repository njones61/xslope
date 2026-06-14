#Janbu Simplified Method

Janbu's Simplified Method is one of the earliest and most commonly used slope stability analysis techniques based on limit equilibrium principles. It provides an approximate factor of safety (FS) for a mass of soil sliding along a specified failure surface by assuming partial force equilibrium . The method was developed by Norwegian engineer **Nils Janbu** in the 1950s to simplify slope stability analysis for circular and non-circular failure surfaces. The factor of safety computed by Janbu's simplified method typically **underestimates** the true FS and is therefore considered conservative. A correction factor was later introduced to compensate for the oversimplifications inherent in the method. Because it satisfies only force equilibrium and leans on this empirical correction, Janbu's method is best suited to quick or preliminary estimates — particularly on non-circular surfaces — with the result confirmed against a complete-equilibrium method such as Spencer before it is relied on for design. 

The forces on the slice considered by the Janbu method are as follows:

>>![oms_slice.png](images/oms_slice.png){width=200px}

For the Janbu method, the inter-slice **shear** forces are assumed to be zero (the inter-slice normal forces are retained). The method satisfies overall **horizontal force equilibrium** of the sliding mass but does not satisfy moment equilibrium. It can be used on both circular and non-circular failure surfaces. Because the inter-slice shear is neglected, the base normal force is found from vertical equilibrium of each slice — the same iterative relation used by Bishop's method.

**Base normal force (vertical equilibrium)** 

Summing forces vertically on a slice, neglecting the inter-slice shear, with the mobilized base shear:

>>$S = \tfrac{1}{F}\left[c\,\Delta\ell + (N - u\,\Delta\ell)\tan\phi\right]$

>>$N\cos\alpha + S\sin\alpha = W$

Solving for the (total) base normal force $N$ gives the Bishop/Janbu form, with $F$ appearing through $m_\alpha$:

>>$N = \dfrac{W - \dfrac{1}{F}\left(c\,\Delta\ell - u\,\Delta\ell\tan\phi\right)\sin\alpha}{m_\alpha}, \qquad m_\alpha = \cos\alpha + \dfrac{\sin\alpha\,\tan\phi}{F}  \qquad (1)$

The mobilized base shear capacity (resisting force along the base) uses the effective normal $N' = N - u\,\Delta\ell$:

>>$S\,F = c\,\Delta\ell + (N - u\,\Delta\ell)\tan\phi   \qquad (2)$

**Factor of safety (horizontal force equilibrium)**

Summing horizontal forces over the whole sliding mass, the inter-slice normal forces cancel and the horizontal components of the base normal ($N\sin\alpha$) and base shear ($S\cos\alpha$) must balance, $\sum S\cos\alpha = \sum N\sin\alpha$. With $S = \tfrac{1}{F}\left[c\,\Delta\ell + (N - u\,\Delta\ell)\tan\phi\right]$ this gives:

>>$F = \dfrac{\sum \left[c\,\Delta\ell + (N - u\,\Delta\ell)\tan\phi\right]\cos\alpha}{\sum N\sin\alpha}     \qquad (3)$

Because $N$ depends on $F$ through $m_\alpha$ in equation (1), equations (1) and (3) are solved by iteration (repeated substitution, starting from $F = 1$) until $F$ converges. This base value is then multiplied by Janbu's correction factor $f_o$ below. (This formulation matches the force-equilibrium factor of safety $F_f$ at zero inter-slice-shear, i.e. the Janbu Simplified value reported by codes such as GeoStudio SLOPE/W.)

## Correction Factor $f_o$

To account for the neglect of interslice forces and moment equilibrium, Janbu proposed a correction factor:

The correction factor $f_o$ is based on the following relationship:

>>![janbu_correction_factor.jpg](images/janbu_correction_factor.jpg){width=500px}

The d/L ratio is the ratio of the distance from the center of the failure surface to the point of interest (d) and the length of the failure surface (L). The correction factor is used to account for the fact that the Janbu method does not satisfy moment equilibrium. The correction factor is a function of the d/L ratio and is calculated using the following equation:

>>$f_o = 1 + b_1  * \left[\dfrac{d}{L} - 1.4 * \left(\dfrac{d}{L}\right)^2\right]    \qquad (4)$

The $b_1$ value is a function of the soils in the slope and is found as follows:


| Soil Type                           | $b_1$ |
|-------------------------------------|-------|
| c-only soil (undrained, $\phi$ = 0) | 0.67  |
| c-$\phi$ soil                       | 0.5   |
| $\phi$-only soil (no cohesion)      | 0.31  |

The polynomial in equation (4) is a fit to Janbu's design chart and is only valid up to the chart domain. It reaches a maximum at $d/L = 1/2.8 \approx 0.357$ and turns over (eventually dropping below 1) for larger $d/L$, which is unphysical. Since $d/L$ is geometrically unbounded, XSLOPE clamps it to that peak so $f_o$ saturates at its maximum rather than decreasing for very deep surfaces. Typical slip surfaces have $d/L$ well below the peak (the verification benchmarks are $\approx 0.13$–$0.20$), so this only guards pathological geometries.

This correction attempts to mimic the effects of moment balance and interslice forces without modeling them directly. Thus, the final factor of safety for Janbu's Simplified Method becomes:

>>$F_{corr} = f_o \cdot F   \qquad (5)$

Where $F_{corr}$ is the corrected factor of safety, $F$ is the base factor of safety from equation (3), and $f_o$ is the correction factor.

## Complete Formulation with Extra Forces

For a complete implementation of Janbu's Simplified Method, we need to consider additional forces acting on the slice. The full set of forces are shown in the following figure:

>>![oms_complete.png](images/oms_complete.png)

Where:

>>$D$ = distributed load resultant force <br>
$\beta$ = inclination of the distributed load (perpendicular to slope) <br>
$kW$ = seismic force for pseudo-static seismic analysis <br>
$c.g.$ = center of gravity of the slice <br>
$P$ = reinforcement force on base of slice <br>
$T$ = tension crack water force <br>
$H$ = pile/pier force at point $e$ on the failure surface <br>
$\theta_p$ = angle of pile force from horizontal (positive = counterclockwise/upward) <br>

Each of these forces is described in detail in the [Ordinary Method of Slices (OMS)](oms.md) section. The external forces enter in two places: their **vertical** components modify the base normal force (vertical equilibrium of the slice), and their **horizontal** components enter the overall horizontal force balance.

**Effective base normal (vertical equilibrium with external forces).** Adding the vertical components of the distributed load ($D\cos\beta$), reinforcement ($P\sin\alpha$), and pile force ($H\sin\theta_p$) to the vertical equilibrium of equation (1) gives the effective base normal $N'$:

>>$N'  = \dfrac{W + D\cos\beta - P\sin\alpha - H\sin\theta_p - u\,\Delta\ell\cos\alpha - \dfrac{1}{F}c\,\Delta\ell\sin\alpha}{m_\alpha}  \qquad (6)$

The total base normal is $N = N' + u\,\Delta\ell$, and the mobilized base shear capacity is $c\,\Delta\ell + N'\tan\phi$.

**Factor of safety (horizontal force equilibrium with external forces).** The horizontal force balance now includes the horizontal components of the external forces: the seismic force $kW$ (horizontal), the distributed-load horizontal component $D\sin\beta$, and the tension-crack water force $T$ are driving; the reinforcement $P$ and the pile horizontal component $H\cos\theta_p$ are resisting. The reinforcement and pile forces are known applied forces (not shear strength) and so are not divided by $F$ — they appear directly in the denominator:

>>$F = \dfrac{\sum \left[c\,\Delta\ell + N'\tan\phi\right]\cos\alpha}{\sum N\sin\alpha + \sum kW + \sum D\sin\beta + \sum T - \sum P - \sum H\cos\theta_p}   \qquad (7)$

with $N = N' + u\,\Delta\ell$. As in the basic case, equations (6) and (7) are solved together by iteration on $F$. Note that $T$ only applies to the side of the uppermost slice ($T = 0$ for all other slices).

Once again, the correction factor $f_o$ is applied to account for the neglect of inter-slice shear as shown in equation (5) above.

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

