# Line of Thrust

For complete equilibrium procedures like Spencer's method, it is common to compute and display a "line of thrust" to 
assess the slope stability solution. The line of thrust is a line connecting the locations of the side forces acting 
on the sides of the slices. The line of thrust is generally expected to be about 1/3 of the way up from the base of 
the slice. However, if tension is present in the slice, the line of thrust may be higher or lower, typically at the 
top of the slope. 
After computing 
and plotting the line of thrust, if significant deviations from the expected location are observed, it may be 
necessary to add a tension crack to get rid of the tension. Soils cannot generally withstand tension, and the 
inclusion of tension is unconservative as it adds to resistance of the slope to failure and may artificially 
increase the factor of safety.

## Forces Acting on the Slice

The line of thrust is calculated by computing the location of the resultant force acting on each slice. To do this, 
we must first compute the forces acting on each slice, including the side forces. Consider the following slice diagram:

![slice_ex.png](images/slice_ex.png){width=350px }

Note that the side forces are defined by the vertical ($X_i$) and horizontal ($E_i$) components of the force acting on 
the slice. The side forces can also be defined by the magnitude ($Z_i$) and the angle of inclination ($\theta$) of 
the force as follows:

![slice_ztheta.png](images/slice_ztheta.png){width=350px }

The shear force ($S_i$) acting on the bottom of the slice is the mobilized shear strength of the soil, which is equal to:

>>$S = \tau_m\Delta l$

where $\tau_m$ is the mobilized shear strength of the soil and $\Delta l$ is the length of the slice. The mobilized shear strength is equal to:

>>$\tau_m = \dfrac{c + (\sigma-u) tan\phi}{F}$
 
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

>>$S = \left[c_m + (\sigma')tan\phi_m\right]\Delta l$

>>$S = c_m\Delta l + N'tan\phi_m$

where:

- $N'$ = the effective normal force acting on the slice = $\sigma'\Delta l$

## Solving for the Side Forces

For the **force equilibrium** method, we solve for the side forces using the iterative method used to solve for the factor of safety. We simply perform the force equilibrium calculations for each slice the factor of safety that satisifies equilibrium. This process is described on the [Force Equilibrium page](force_eq.md).

For **Spencer's method**, the process is simpler. The resultant force $Q_i$ that is the central component of Spencer's method represents the resulttant of the of the side forces on each slice. Specifically, 
 
>>$Q_i = Z_i - Z_{i+1}$

For the first slice, the left side force is zero. So we have:

>>$Z_{i} = 0$

>>$Z_{i+1} = Z_i - Q_i = - Q_i$

For each subsequent slice, we have:

>>$Z_{i+1} = Z_i - Q_i$

where $Z_i$ is the left side force computed from the prior slice. We continue this process until we reach the last slice and for the last slice, we have:

>>$Z_{i+1} = 0$

## Solving for the Thrust Line

Now that we have the forces acting on each slice, we can compute the locations of the resultant sides force acting on each slice. The side force locations are found by summing moments about the center of the base of each slice. W, S, and N all go through the base so the only components in the moment equilibrium equation are the side forces. Also, rather than working in terms of $Z$ and $theta$, we will work in terms of the vertical and horizontal components of the side forces. The horizontal and vertical components are given by:

>>$X_{i} = Z_i sin(\theta)$

>>$E_{i} = Z_i cos(\theta)$

where $X_i$ is the vertical component of the side force and $E_i$ is the horizontal component of the side force. 
Next, we define $\Delta y_{i}$ and $\Delta y_{i+1}$ as the distance from the center point of the base of the slice 
up to the left and right side forces, $E_{i}$ and $E_{i-1}$ respectively, and $\Delta x$ is the width of the slice. Assuming CCW rotation is positive (right-hand rule) we can write the moment equilibrium equation as:

>>$\sum M = 0 \Rightarrow - E_{i} \Delta y_{i} - X_{i} \dfrac{\Delta x}{2} + E_{i+1} \Delta y_{i+1} - X_{i+1} \dfrac
> {\Delta x}{2}= 0$

If we start with slice 1 on the left side, the left side force is zero so we have one unknown ($\Delta y_{i+1}$) and one equation. Then on the next slice, the left side moment arm is known and the right side moment arm is unknown. So again we have one equation and one unknown ($\Delta y_{i+1}$) which we can solve for as follows:

>>$\Delta y_{i+1} = \dfrac{E_{i} \Delta y_{i} + X_{i} \dfrac{\Delta x}{2} + X_{i+1} \dfrac{\Delta x}{2}}{E_{i+1}}$

We can continue this process until we reach the top slice where the right side moment arm is zero. On the last slice, the moment equation should balance using the known left side moment arm, but it may not close exactly due to accumulated rounding errors. Another alternative is to start from the left side and sweep to the the right side and then start from the right side and sweep to the left side. This will give two different moment arms for the same slice, but they should be very close. The average of the two moment arms can be used to compute the location of the resultant side force. 

## Solving for the Normal Force

It is useful to solve for the effective normal force at the same time that the line of thrust is calculated. These 
can be plotted (as stresses) along with the line of thrust. Negative values are another indicator of tension in the 
slice. If the force equilibrium method is used to solve for the side forces, the normal forces on the base of the slice 
are calculated at the same time. But if Spencer's method $Q_i$ values are used to get the side forces, it is 
necessary to solve separately for the effective normal force. We can do this as follows:

>>$\sum F_y = 0 \Rightarrow \left[c_m \Delta l + N' tan(\phi_m)\right] sin(\alpha) + (N' + u \Delta l) cos(\alpha) -  W + X_{i} - X_{i+1} = 0$

>>$c_m \Delta l  sin(\alpha) + N' tan(\phi_m) sin(\alpha) + N' cos(\alpha) + u \Delta l cos(\alpha) - W + X_{i} - X_{i+1} = 0$

>>$c_m \Delta l  sin(\alpha) + N' \left[tan(\phi_m) sin(\alpha) +  cos(\alpha)\right] + u \Delta l cos(\alpha) - W + X_{i} - X_{i+1} = 0$

Solving for $N'$ gives:

>>$N' = \dfrac{- c_m \Delta l  sin(\alpha) - u\Delta l cos(\alpha) + W - X_{i} + X_{i+1}}{tan(\phi_m) sin(\alpha) + cos(\alpha)}$