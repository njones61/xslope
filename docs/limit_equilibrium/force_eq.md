# The Force Equilibrium Method

The force equilibrium method is a common method used to analyze the stability of slopes. It works for both circular 
and non-circular surfaces. It satisfies the force equilibrium equations for each slice, but it does not satisfy the 
moment equilibrium equations. There are several variations of the force equilibrium method, depending on the 
assumptions used for the side force inclination. In slope tools, two variations of the force equilibrium method are 
supported: the **Lowe and 
Karafaith** method and the **US Army Corps of Engineers** method. Both methods use the same equations, but they differ in 
the side force assumptions.

## Forces Acting on the Slice

To formulate the force equilibrium equations, consider the following slice diagram:

![slice_ex.png](images/slice_ex.png){width=350px }

Note that the side forces are defined by the vertical ($X_i$) and horizontal ($E_i$) components of the force acting on 
the slice. The side forces can also be defined by the magnitude ($Z_i$) and the angle of inclination ($\theta_i$) of 
the force as follows:

![slice_ztheta.png](images/slice_ztheta.png){width=350px }

The shear force ($S$) acting on the bottom of the slice is the mobilized shear strength of the soil, which is equal to:

>>$S = \tau_m  \Delta \ell$

where $\tau_m$ is the mobilized shear strength of the soil and $\Delta \ell$ is the length of the slice. The mobilized shear strength is equal to:

>>$\tau_m = \dfrac{c + (\sigma-u)tan\phi}{F}$
 
>>$\tau_m = c_m + \sigma' tan\phi_m$

where:

- $c$ = the cohesion of the soil
- $c_m$ = the mobilized cohesion of the soil = $c/F$
- $\sigma$ = the normal stress acting on the slice
- $\sigma'$ = the effective normal stress acting on the slice = $\sigma - u$
- $u$ = the pore water pressure acting on the slice
- $\phi$ = the angle of internal friction of the soil
- tan$\phi_m$ = the mobilized friction of the soil = tan$\phi/F$
- $F$ = the factor of safety

Inserting $\tau_m$ into the equation for $S$ gives:

>>$S = \left[c_m + (\sigma')  tan\phi_m\right]  \Delta \ell$

>>$S = c_m  \Delta \ell + N'  tan\phi_m$

where:

>>$N'$ = the effective normal force acting on the slice = $\sigma'  \Delta \ell$

>>$N$ = the total normal force acting on the slice = $N + u \Delta \ell$

## Solving for Unknown Forces

To satisfy force equilibrium, we need to solve for the unknown forces acting on each slice. The ultimate goal is to solve for the factor of safety, F. However, the mobilized shear strength ($c_m$ and $tan \phi_m$) is a function of the factor of safety. Therefore, we first assume a value for the factor of safety and the solve the equilbrium equations for all slices. If the forces balance, we are done. If not, we adjust the factor of safety and repeat until balance is achieved. Also, the side force inclination ($\theta$) is required for all slice boundaries. There are a number of methods for 
establishing the side force inclination. These methods will be reviewed below. For now, we will assume we have a value for the side force inclinations.

So our next step is to use the force equilibrium equations to solve for the unknown forces acting on each slice. We do this by starting on the bottom slice (slice 1) and working our way up to the top slice (slice n). For the first slice, we have:

![slice_fe_left.png](images/slice_fe_left.png){width=400px }


The unknowns in this case are the normal force ($N'$) and magnitude of the side force ($Z_{i+1}$). We do know the side 
force inclination ($\theta_{i+1}$). So we have two equations: $\sum F_x = 0$ and $\sum F_y = 0$ and two unknowns. After 
solving for the unknows, we proceed to the next slice:

![sclice_fe_next.png](images/sclice_fe_next.png){width=700px }

The side force on the left side ($Z_{i-1}$) is known from the previous slice. The side force on the right side ($Z_
{i+1}$) is unknown. The effective normal force ($N'$) is also unknown. So again, we have two equations and two unknowns. 

We can continue this process until we reach the top slice:

![slice_fe_last.png](images/slice_fe_last.png){width=400px }

At this point, there is no side force on the right side, so we have two equations and one unknown ($N'$), if both 
equilibrium equations are balanced (no residual forces), then the factor of safety is correct. If not, we need to adjust the factor of safety and repeat the process.

## Matrix Solution for Unknown Forces

As shown in the previous section, we need to solve for two unknowns at each slice. For the general case, we can set up equations to solve for the unknowns as follows. First, we sum forces in the x-direction:

>>$\sum F_x = 0 \Rightarrow \left[c_m \Delta \ell + N'  tan(\phi_m)\right]  cos(\alpha) - (N' + u \Delta \ell)  sin(\alpha) + Z_{i}
>   cos(\theta_i) - Z_{i+1}  cos(\theta_{i+1}) = 0$

>>$c_m \Delta \ell   cos(\alpha) + N'  tan(\phi_m)  cos(\alpha) - N'  sin(\alpha) - u \Delta \ell  sin(\alpha) - + Z_{i}  cos
> (\theta_{i}) - Z_
> {i+1}  cos(\theta_{i+1}) = 0$

>>$c_m \Delta \ell   cos(\alpha) + N' \left[tan(\phi_m) cos(\alpha) - sin(\alpha)\right] - u \Delta \ell sin(\alpha) + Z_{i}
>  cos
> (\theta_{i}) - Z_{i+1} cos(\theta_{i+1}) = 0$

Rearranging in terms of our two unknows ($N'$ and $Z_{i+1}$) gives:

>>$N' \left[tan(\phi_m) cos(\alpha) - sin(\alpha)\right] - Z_{i+1} cos(\theta_{i+1}) = - c_m \Delta \ell  cos(\alpha) + u 
> \Delta \ell sin(\alpha) - Z_{i} cos(\theta_i)   \qquad (1)$

Likewise, we can sum forces in the y-direction:

>>$\sum F_y = 0 \Rightarrow \left[c_m \Delta \ell + N' tan(\phi_m)\right] sin(\alpha) + (N' + u \Delta \ell) cos(\alpha) - 
> W + Z_{i} sin(\theta_{i}) - Z_{i+1} sin(\theta_{i+1}) = 0$

>>$c_m \Delta \ell  sin(\alpha) + N' tan(\phi_m) sin(\alpha) + N' cos(\alpha) + u \Delta \ell cos(\alpha) - W + Z_{i} sin
> (\theta_{i}) - 
> Z_{i+1} sin(\theta_{i+1}) = 0$

>>$c_m \Delta \ell  sin(\alpha) + N' \left[tan(\phi_m) sin(\alpha) +  cos(\alpha)\right] + u \Delta \ell cos(\alpha) - W + 
> Z_{i} sin
> (\theta_{i}) - Z_{i+1} sin(\theta_{i+1}) = 0$

Rearranging in terms of our two unknows ($N$ and $Z_{i+1}$) gives:

>>$N' \left[tan(\phi_m) sin(\alpha) + cos(\alpha)\right] - Z_{i+1} sin(\theta_{i+1}) = -c_m \Delta \ell  sin(\alpha) - u\Delta \ell cos(\alpha) + W - Z_{i} sin(\theta_{i})   \qquad (2)$ 

Now we can take equations (1) and (2) and rearrange them into a matrix form. We can write the two equations as:

>>$Ax = b$

where:

- $A$ is a 2x2 matrix of coefficients
- $x$ is a 2x1 vector of unknowns
- $c$ is a 2x1 vector of constants

The matrix $A$ is given by:

>>$A = \begin{bmatrix}tan(\phi_m) cos(\alpha) - sin(\alpha) & -cos(\theta_{i+1})\\tan(\phi_m) sin(\alpha) + cos(\alpha) & -sin(\theta_{i+1})\end{bmatrix}    \qquad (3)$

The vector $x$ is given by:

>>$x = \begin{bmatrix}N'\\Z_{i+1}\end{bmatrix}   \qquad (4)$

The vector $b$ is given by:

>>$b = \begin{bmatrix}- c_m \Delta \ell  cos(\alpha) + u \Delta \ell sin(\alpha) - Z_{i} cos(\theta_{i})\\-c_m \Delta \ell  sin(\alpha) - u\Delta \ell cos(\alpha) + W - Z_{i} sin(\theta_{i})\end{bmatrix}   \qquad (5)$

The matrix equation can then be solved for the two unknowns
($N$ and $Z_{i+1}$) using the numpy linalg method. The solution is given by:

```python
import numpy as np

N[i], Z[i + 1] = np.linalg.solve(A, b)
```

## Complete Formulation

For a complete implementation of the Force Equilibrium method, we need to consider additional forces acting on the slice. The full set of forces are shown in the following figure:

![slice_fe_complete.png](images/slice_fe_complete.png){width=400px}

where:

>$D$ = distributed load resultant force <br>
$\beta$ = inclination of the distributed load (perpendicular to slope) <br>
$kW$ = seismic force for pseudo-static seismic analysis <br>
$c.g.$ = center of gravity of the slice <br>
$P$ = reinforcement force on base of slice <br>
$T$ = tension crack water force <br>

Each of these forces is described in detail in the [Ordinary Method of Slices (OMS)](oms.md) section. The forces $D$, $kW$, $P$, and $T$ are included in the Force Equilibrium method as follows.

Once again, we begin by summing forces in the x-direction, but now we include the additional forces:

>>$\sum F_x = 0 \Rightarrow \left[c_m \Delta \ell + N' \tan (\phi_m + P)\right] \cos (\alpha) - (N' + u \Delta \ell) \sin (\alpha) + Z_{i} \cos (\theta_i) - Z_{i+1} \cos (\theta_{i+1}) + D \sin \beta - kW - T = 0$

>>$c_m \Delta \ell \cos (\alpha) + N' \tan (\phi_m) \cos (\alpha) + P \cos (\alpha) - N' \sin (\alpha) - u \Delta \ell \sin (\alpha) + Z_{i} \cos (\theta_i) - Z_{i+1} \cos (\theta_{i+1}) + D \sin \beta - kW - T= 0$

>>$c_m \Delta \ell \cos (\alpha) + N' \left[\tan (\phi_m) \cos (\alpha) - \sin (\alpha)\right] + P \cos (\alpha) - u \Delta \ell \sin (\alpha) + Z_{i} \cos (\theta_i) - Z_{i+1} \cos (\theta_{i+1}) + D \sin \beta -kW -T  = 0$

Rearranging in terms of our two unknows ($N'$ and $Z_{i+1}$) gives:

>>$N' \left[\tan (\phi_m) \cos (\alpha) - \sin (\alpha)\right] - Z_{i+1} \cos (\theta_{i+1}) = - c_m \Delta \ell \cos (\alpha) - P \cos (\alpha) + u \Delta \ell \sin (\alpha) - Z_{i} \cos (\theta_i) - D \sin \beta + kW + T   \qquad (6)$

It should be noted that the tension crack water force ($T$) only applies to right side of the top slice on a left-facing slope. For a right-facing slope, the tension crack water force is applied to the left side of the top slice and would act in the opposite direction. Therefore, the sign on $T$ would be negative in that case.

Likewise, we can sum forces in the y-direction:

>>$\sum F_y = 0 \Rightarrow \left[c_m \Delta \ell + N' \tan (\phi_m + P)\right] \sin (\alpha) + (N' + u \Delta \ell) \cos (\alpha) - W + Z_{i} \sin (\theta_{i}) - Z_{i+1} \sin (\theta_{i+1}) - D \cos \beta = 0$

>>$c_m \Delta \ell \sin (\alpha) + N' \tan (\phi_m) \sin (\alpha) + P \sin (\alpha) + N' \cos (\alpha) + u \Delta \ell \cos (\alpha) - W + Z_{i} \sin (\theta_{i}) - Z_{i+1} \sin (\theta_{i+1}) - D \cos \beta = 0$

>>$c_m \Delta \ell \sin (\alpha) + N' \left[\tan (\phi_m) \sin (\alpha) +  \cos (\alpha)\right] + P \sin (\alpha) + u \Delta \ell \cos (\alpha) - W + Z_{i} \sin (\theta_{i}) - Z_{i+1} \sin (\theta_{i+1}) - D \cos \beta = 0$

Rearranging in terms of our two unknows ($N'$ and $Z_{i+1}$) gives:

>>$N' \left[\tan (\phi_m) \sin (\alpha) +  \cos (\alpha)\right] - Z_{i+1} \sin (\theta_{i+1}) = -c_m \Delta \ell \sin (\alpha) - P \sin (\alpha) - u\Delta \ell \cos (\alpha) + W - Z_{i} \sin (\theta_{i}) + D \cos \beta   \qquad (7)$

Now we can take equations (6) and (7) and rearrange them into a matrix form. We can write the two equations as:

>>$Ax = b$

The matrix $A$ is given by:

>>$A = \begin{bmatrix}\tan (\phi_m) \cos (\alpha) - \sin (\alpha) & -\cos (\theta_{i+1})\\\tan (\phi_m) \sin (\alpha) +  \cos (\alpha) & -\sin (\theta_{i+1})\end{bmatrix}   \qquad (8)$

The vector $x$ is given by:

>>$x = \begin{bmatrix}N'\\Z_{i+1}\end{bmatrix}   \qquad (9)$

The vector $b$ is given by:

>>$b = \begin{bmatrix}- c_m \Delta \ell \cos (\alpha) - P \cos (\alpha) + u \Delta \ell \sin (\alpha) - Z_{i} \cos (\theta_i) - D \sin \beta + kW + T\\-c_m \Delta \ell \sin (\alpha) - P \sin (\alpha) - u\Delta \ell \cos (\alpha) + W - Z_{i} \sin (\theta_{i}) + D \cos \beta\end{bmatrix}   \qquad (10)$

Note that $A$ and $x$ are the same as before, but $b$ has changed to include the additional forces. The matrix equation can then be solved for the two unknowns ($N'$ and $Z_{i+1}$) using the numpy **linalg** method as described above. 

Once again, care must be taken to ensure that the tension crack water force ($T$) is applied correctly based on the slope direction.

## Side Force Inclination Assumptions

The side force inclination is a critical parameter in the force equilibrium method. Two solution methods are 
supported in slope tools:

### Lowe and Karafaith

The Lowe and Karafaith method assumes that the side force inclinations are equal to the average slope of ground 
surface and slip surface as defined by the top and bottom of the slice. 

### US Army Corps of Engineers

The US Army Corps of Engineers method assumes that the side force inclinations are parallel to the slope angle, 
using one of the following methods:

![uscoe_theta.png](images/uscoe_theta.png){width=500px }

In slope tools, first method shown above is used. That is, all interslice forces are parallel to a line connecting the 
bottom of the failure surface to the top of the failure surface. 