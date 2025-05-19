# Ordinary Method of Slices (OMS)

The Ordinary Method of Slices (also known as the Swedish Method) is one of the simplest limit equilibrium techniques used for slope stability analysis. The factor of safety to be calculated directly without iteration. The key assumption for the OMS method is that the side forces can be neglected, meaning interslice shear and normal forces cancel out and are excluded from equilibrium considerations. Thus, the forces on the slice are as follows:

![oms_slice.png](images/oms_slice.png){ width=200px }

Further, only **moment equilibrium** about the center of the slip circle is enforced:

>$\sum M = 0$

Force equilibrium in the horizontal and vertical directions is **not** satisfied.

## Basic Formulation

To solve for the factor of safety, we need to consider the forces acting on the slice:

>$W$ = weight of the slice  
$\alpha$ = inclination of the base of the slice  
$\Delta \ell$ = length of the slice base  
$c$ = cohesion  
$\phi$ = friction angle  
$u$ = pore water pressure  

The normal force on the base is:

>$N = W \cos \alpha$

Shear stress is:

>$\sigma = \dfrac{N}{\Delta \ell} = \dfrac{W \cos \alpha}{\Delta \ell}$

The general expression for the factor of safety $FS$is:
>$FS = \dfrac{\sum (c \Delta \ell + W \cos \alpha \tan \phi)}{\sum W \sin \alpha}$

If $\phi = 0$, the expression simplifies to:
>$FS = \dfrac{\sum c \Delta \ell}{\sum W \sin \alpha}$

This is the same equation used in the Swedish method and the log-spiral method under the same assumptions.

To use effective stress parameters $c', \phi'$, and effective normal stress $\sigma'$, the expression becomes:

>$\sigma' = \dfrac{W \cos \alpha}{\Delta \ell} - u$

>$FS = \dfrac{\sum (c' \Delta \ell + (W \cos \alpha - u \Delta \ell) \tan \phi')}{\sum W \sin \alpha}$

## Alternate (Preferred) Formulation

This equation for FS can produce unconservative results, including negative normal stresses under high pore pressures. The following formulation is preferred. First, we define the effective weight as follows:

>$W' = W - u b, \quad \text{where} \quad b = \Delta \ell \cos \alpha$

Then:

>$N' = W' \cos \alpha = W \cos \alpha - u \Delta \ell \cos^2 \alpha$

And:

>$\sigma' = \dfrac{N'}{\Delta \ell} = \dfrac{W \cos \alpha}{\Delta \ell} - u \cos^2 \alpha$

Substituting:

>$F = \dfrac{\sum \left[ c' \Delta \ell + (W \cos \alpha - u \Delta \ell \cos^2 \alpha) \tan \phi' \right]}{\sum W \sin \alpha}$

This is the **preferred formulation** for OMS.

## Summary

- Applicable only to **circular** slip surfaces.
- **Only moment equilibrium** is satisfied.
- **No iteration** is required.
- **Less accurate** than more complete methods (e.g., Bishop’s or Spencer’s).
- Provides the same solution as the Swedish method when $\phi = 0$.
