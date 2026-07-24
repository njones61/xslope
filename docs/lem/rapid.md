# Rapid Drawdown Analysis

When designing a dam or levee, it is crucial to consider the potential for rapid drawdown of the water level. This can occur due to various factors such as dam failure, sluice gate operation, or other emergency situations. When the water level drops rapidly, it can lead to significant changes in the forces acting on the structure, potentially causing instability or failure. Specifically, the rapid drawdown eliminates the hydrostatic pressure on the upstream side of the structure, removing the buoyancy effect of the water, which can lead to a sudden increase in effective stress and shearing under undrained conditions.

![pool_normal.png](rapid_images/pool_normal.png)

![pool_lowered.png](rapid_images/pool_lowered.png)

To simulate this scenario, we utilize a "multi-stage" approach in our analysis. In the first stage, we apply the initial conditions with the water level at its full height and use the consolidation stresses resulting from these conditions to determine the undrained strength of the soils. In the second stage, we use these undrained strengths and the loading conditions corresponding to the lowered water level and determine the factor of safety for rapid drawdown conditions. The **xslope** package provides a convenient way to perform this analysis using the `rapid_drawdown` function. The equations and methodology used in this function are described in the following sections.

## When Does Rapid Drawdown Apply?

Rapid drawdown occurs when the pool is lowered rapidly enough that pore water pressures in some of the soils cannot dissipate quickly enough to maintain stability. To determine when rapid drawdown may apply, Duncan, et al. (1992) suggested using the dimensionless time factor, T, from consolidation theory. T is defined as:

>>$T = \dfrac{c_v t}{D^2}   \qquad (1)$

where:

>> $c_v$ = coefficient of consolidation [L²/t]<br>
$t$ = time since the pool was lowered [t]<br>
$D$ = maximum drainage distance [L]

Then, we use the following rubric:

|  T   | Rapid Drawdown Applicable?                         |
|:----:|----------------------------------------------------|
| >3.0 | Drainage should be sufficient - not rapid drawdown |  
| <3.0 | Assume rapid drawdown                              |

Approximate values of $c_v$ are as follows:

| Soil Type    | $c_v$ [ft²/day] |
|--------------|----------------|
| Coarse Sand  | > 10,000       |
| Fine Sand    | 100 - 10,000   |
| Silty Sand   | 10 - 1000      |
| Silt         | 0.5 - 100      |
| Compacted Clay | 0.05 - 5       |
| Soft Clay    | < 0.02         |

_(Taken from Duncan et al, 1992)_

For example, consider the following levee. The maximum drainage distance for pore pressures to dissipate is 35 ft, and the coefficient of consolidation for the compacted clay is 2 ft²/day. If the pool is lowered over 50 days, we can calculate T as follows:

![example1.png](rapid_images/example1.png)

Since T = 0.09 < 3.0, we assume rapid drawdown applies.

Here is another example involving a dam with a high permeability shell and a low permeability core. In this case, we need to examine the core and the shell separately. The core consists of compacted clay with a maximum drainage distance of 14 ft and a coefficient of consolidation of 2 ft²/day. The shell consists of fine sand with a maximum drainage distance of 31 ft and a coefficient of consolidation of 200 ft²/day. If the pool is lowered over 50 days, we can calculate T for each soil as follows:

![example2.png](rapid_images/example2.png)

In this case T for the core = 0.5 < 3.0, so we assume rapid drawdown applies to the core. T for the shell = 10.4 > 3.0, so we  assume rapid drawdown does not apply to the shell.

## Three-Stage Analysis

In xslope, we use the three-stage analysis method developed by Duncan, Wright, and Wong (1990) to analyze rapid drawdown conditions. This technique is described in detail in Soil Strength and Slope Stability, 2nd Edition, Duncan, Wright, and Brandon. The three stages are:

>>1) Compute effective stresses based on pre-drawdown conditions<br>
2) Compute FS for post-drawdown conditions using total stress analysis. Undrained strength is a function of effective stress from stage 1.<br>
3) Compute FS for post-drawdown conditions using drained strengths

Stage 3 only ever *lowers* a slice's strength — it swaps in the drained strength where that is smaller than the undrained one — so its FS can never exceed the Stage 2 FS. The final FS for rapid drawdown is therefore the Stage 3 value when Stage 3 runs, and the Stage 2 value when it is not needed. Equivalently, and how xslope computes it: the lower of the two.

### Stage 1: Pre-Drawdown Conditions

Stage 1 analysis is performed for conditions prior to drawdown. The objective is to estimate the effective stresses along the slip surface prior to drawdown. This is also referred to as the consolidation stage.

![stage1.png](rapid_images/stage1.png)

To get the consolidation stresses, we perform a normal slope stability analysis using the conditions corresponding to the pool at the full level. The objective of this stage of the analysis is not the factor of safety, but the effective normal stresses on the base of each slice, which is one of the products of the slope stability functions (oms, bishop, spencer, etc.).

![slice_forces.png](rapid_images/slice_forces.png){width=300px}

Using the normal force ($N$) found by the solution, we can calculate the effective normal stress on the failure plane at the end of the consolidation stage as:

>>$\sigma_{fc} = \dfrac{N}{\Delta \ell} - u$

Since the values returned by the solvers are actually effective normal forces ($N'$), we can also calculate as follows:

>>$\sigma'_{fc} = \dfrac{N'}{\Delta \ell}   \qquad (2)$

The precise mechanics of how $N'$ is computed vary according to the solver used. The equations and methodology for computing $N'$ is described in the **Methods** section of this documentation.

Likewise, the shear stress on the failure plane at the end of the consolidation stage is:

>>$\tau_{fc} = \dfrac{S}{\Delta \ell}$

>>$\tau_{fc} = \dfrac{1}{F}(c' + \sigma'_{fc} \tan \phi')   \qquad (3)$

Note that this is the mobilized shear strength and $F$ is the factor of safety returned by the solver in the Stage 1 solution.

## Stage 2 - Compute FS for Post-Drawdown Conditions

In Stage 2, we perform a total stress analysis using undrained shear strengths for the low-K zones. The undrained strengths are estimated using the consolidation stresses computed in Stage 1, and the driving forces are calculated using post-drawdown conditions.

![stage2.png](rapid_images/stage2.png)

The process of finding the appropriate undrained strength for the Stage 2 analysis based on the consolidation stresses from Stage 1, involves a variation of the Mohr-Coulomb failure envelope. Consider the following diagram:

![mc_tauff.png](rapid_images/mc_tauff.png){width=610px}

The $\tau_{ff}$ is the shear stress on the failure plane at failure. $\sigma'_{1f}$ and $\sigma'_{3f}$ are major and minor principal effective stresses at failure. For our application, we need to correlate $\tau_{ff}$ with $\sigma'_{fc}$. To do this, we create a new plot as follows:

![tauff_vs_sigmafc.png](rapid_images/tauff_vs_sigmafc.png)

These two curves are function of the stress anisotropy during the consolidation stage. The $K_c = K_f$ curve with $c'$ and $\phi'$ corresponds to the case where the ratio of $\sigma'_1$ and $\sigma'_3$ is on the verge of failure. The other curve with $d$ and $\psi$ corresponds to the case where $\sigma'_1 = \sigma'_3$ during the consolidation stage. Note that the $c'$ and $\phi'$ values for the $K_c = K_f$ curve can be obtained from a $CU$ triaxial test or they can be obtained from a $CD$ or $\overline{CU}$ triaxial test. The $K_c = 1$ line represents $\tau_{ff}$ vs. $\sigma'_3$ from $CU$ triaxial tests. The $\tau_{ff}$ values are found using the equation in the previous figure above. Thus, both lines can be obtained from one set of $\overline{CU}$ triaxial tests.

When performing rapid drawdown analysis with the **xslope** package, a set of $d$ and $\psi$ values should be entered in the material properties table for each soil with poor drainage and assumed to apply to rapid drawdown conditions. The more freely draining soils are not included in the analysis, and the $d$ and $\psi$ values should be left blank.

### Anisotropy and K Interpolation

For the stresses on the base of each slice, the actual $K$ ratio should vary somewhere between $K_c = 1$ and $K_c = K_f$. For example, if $K_c = 2$, the middle line in this figure would apply.

![k_interp.png](rapid_images/k_interp.png){width=720px}

If we define $K_1$ as the stress ratio on the failure plane at the end of stage 1, we interpolate the two curves and extract the appropriate $\tau_{ff}$ value to use in the Stage 2 calculations. First of all, we compute $K_1$ as follows:

>>$K_1 = \dfrac{\sigma'_{fc} + \tau_{fc}[(\sin \phi' + 1) / \cos \phi']}{\sigma'_{fc} + \tau_{fc}[(\sin \phi' - 1) / \cos \phi']}   \qquad (4)$

Where $\sigma'_{fc}$ and $\tau_{fc}$ are the effective normal and shear stresses on the shear plane at the end of consolidation and are found using equations (2) and (3) above. Equation (4) assumes that the orientation of the principal stresses at the end of consolidation is the same as the orientation of the principal stresses at failure. Values of the undrained shear strength for the effective consolidation stress $\sigma'_{fc}$ and a consolidation ratio $K_c = K_1$ are obtained by interpolating the $K_c = 1$  and $K_c = K_f$  failure envelopes using the following equation:

>>$\tau_{ff} = \dfrac{(K_f - K_1) \tau_{ff(K_c=1)} + (K_1 - 1) \tau_{ff(K_c=K_f)}}{K_f - 1}  \qquad (5)$

The two $\tau_{ff}$ values are found by evaluating the two envelopes using the $\sigma'_{fc}$ found from stage 1. The $K_f$ value used in this equation can be found from:

>>$K_f = \dfrac{(\sigma'_{fc} + c' \cos \phi')(1 + \sin \phi')}{(\sigma'_{fc} - c' \cos \phi')(1 - \sin \phi')}   \qquad (6)$

### Negative Stresses

The $\sigma'_3$ values can become negative, leading to a negative (and meaningless) $K$ value.

![k_interp_neg.png](rapid_images/k_interp_neg.png)

In this case, we use the lower of the $K_c = 1$  and $K_c = K_f$  curves. Negative effective stresses can be checked using the following two equations:

>>$\sigma'_{3c} = \sigma'_{fc} + \tau_{fc} \dfrac{\sin \phi' - 1}{\cos \phi'}  \qquad (7)$ (for the $K_c = 1$ envelope)

>>$\sigma'_{3c} = (\sigma'_{fc} - c' \cos \phi') \dfrac{1 - \sin \phi'}{\cos^2 \phi'}  \qquad (8)$ (for the $K_c = K_f$ envelope)  

If either is negative **or zero**, no interpolation is required and we use the lower of the two strength values coming from the two curves. Zero must be excluded as well as negatives, because $\sigma'_{3c}$ from equation (7) is precisely the *denominator* of $K_1$ in equation (4).

These two equations go negative for different reasons, and only one of them involves cohesion:

- Equation (8) carries the $c' \cos \phi'$ term, so a **significant cohesion** drives it negative. This is the case illustrated above.
- Equation (7) has **no cohesion term at all**. It goes negative when $\tau_{fc}$ is large — and since $\tau_{fc} = \frac{1}{F}(c' + \sigma'_{fc} \tan \phi')$, that means a **low Stage 1 factor of safety**.

### The Stage 1 slope must be stable

Equation (5) *interpolates* between the $K_c = 1$ and $K_c = K_f$ envelopes, which bound the physically possible consolidation states. So $K_1$ must satisfy

>> $1 \le K_1 \le K_f$

$K_1$ increases monotonically with $\tau_{fc}$, and $\tau_{fc} \propto 1/F$. Substituting the failure value $\tau_{fc} = c' + \sigma'_{fc} \tan \phi'$ into equation (4) gives $K_1 = K_f$ exactly. Therefore

>> $K_1 > K_f \iff F < 1$ in Stage 1

which is to say: **the full-pool slope is already failing.** Its mobilized shear stress lies above the failure envelope, so there is no equilibrium consolidation stress state for Stage 2 to use, and equation (5) would extrapolate beyond $K_c = K_f$ rather than interpolate. Note this condition is not caught by the negative-stress check above — $\sigma'_{3c}$ from equation (7) remains positive well below $F = 1$.

Because $F$ is a single global value, $K_1 > K_f$ holds for every slice at once. xslope therefore checks the Stage 1 factor of safety directly, and `rapid_drawdown` returns an error rather than a factor of safety when $F < 1$. Rapid drawdown presupposes a slope that is stable before the pool is lowered.

### FS Calculations

After the undrained strengths are computed for low K zones, the factor of safety ($FS$) is calculated. Drained strengths are used for high K zones, and driving forces (based on weight of slices) are calculated using post-drawdown conditions. 

## Stage 3 - Check Drained Strengths

For Stage 3, undrained strengths used for slices in low K zones are compared with drained strengths. If drained strength is lower than the undrained strength, the drained strength is assigned to that slice and the FS is re-calculated. Note that some slices may used drained strength while others may use undrained. If drained > undrained for all slices, there is no need to re-calculate FS. For the drained strength calculations:

>>$\sigma' = \dfrac{N'}{\Delta \ell}   \qquad (9)$

>>$\tau = c' + \sigma' \tan \phi'   \qquad (10)$

Where $N'$ is the effective normal force found in Stage 2 using post-drawdown conditions.

If Stage 3 calculations are required (drained less than undrained in at least one slice), then the FS from Stage 3 = the rapid drawdown FS. If Stage 3 calculations are not required, the FS from Stage 2 = the rapid drawdown FS. Since Stage 3 only substitutes *lower* strengths, its FS never exceeds the Stage 2 FS, so this is the same as taking the lower of the two.

## Inputs and Calculations

To perform a rapid drawdown analysis, the following additional inputs must be provided in the input template Excel file:

| Category | Description |
|----------|-------------|
| Material Properties| $d$ and $\psi$ for all materials to be analyzed using multi-stage approach. |
| Pore Pressures | Two piezometric lines: one for Stage 1 and one for Stage 2 |
| Distributed Loads | Two sets of distributed loads, one for Stage 1 (full pool) and one for Stage 2 (lowered pool) |

In summary, the calculations are done in the following process:

**Stage 1**

1. Using the drained strength ($c'$ and $\phi'$) strength properties for all materials and using the piezometric line and distributed loads for the full pool condition, calculate the factor of safety using the selected solver. This will return a factor of safety (FS) and a set of effective normal forces ($N'$) on the base of the slice. If this Stage 1 FS is less than 1, the slope is already failing at full pool and the analysis is halted — see [The Stage 1 slope must be stable](#the-stage-1-slope-must-be-stable).<br>

2.  Calculate $\sigma'_{fc}$ and $\tau_{fc}$ using equations (2) and (3).<br>

**Stage 2**

3. Calculate $K_1$ via equation (4), using $\sigma'_{fc}$ and $\tau_{fc}$ from step 2.<br>

4. Using $\sigma'_{fc}$ from step 2, calculate $\tau_{ff(K_c=1)}$ from the $d - \psi$ curve and $\tau_{ff(K_c=K_f)}$ from the $c - \phi$ curve.<br>

5. Calculate $K_f$ using equation (6) and $\sigma'_{fc}$ from step 2.<br>

6. Calculate $\tau_{ff}$ from equation (5)<br>

7. Steps 3-6 are performed for each low K soil. For each of these soils, set $c$ = $\tau_{ff}$ and set $\phi = 0$. For the high K soils, use the normal $c - \phi$ values.

8. Using the strengths from step 7 and the pore pressures and distributed loads corresponding to drawdown (lowered pool) conditions, calculate $FS$ using the selected solver. This is the $FS$ for Stage 2.

**Stage 3**

9. For each slice with a low K soil at the bottom, use the $N'$ values calculated in Stage 2 and calculate the drained shear strength ($\tau$) using equations (9) and (10). Compare this value to the $\tau_{ff}$ value used in Stage 2. If the drained strength is smaller for any slice, replace the undrained strength values with the original $c - \phi$ values and recompute $FS$. This is the $FS$ for Stage 3.

10. Compare the $FS$ for Stages 2 & 3. The smaller of the two is the final $FS$ for rapid drawdown conditions.

## Rapid Drawdown from a Transient Solution {#transient-solution}

The three-stage method needs two pore-pressure fields along the slip surface: one for the pre-drawdown, full-pool consolidation state used in [Stage 1](#stage-1-pre-drawdown-conditions), and one for the drawn-down state used in [Stage 2](#stage-2-compute-fs-for-post-drawdown-conditions). The [Inputs and Calculations](#inputs-and-calculations) summary supplies these from two piezometric lines, and a finite-element workflow can instead supply them from two steady seepage solutions saved next to the workbook as `{base}_seep.csv` (full pool) and `{base}_seep2.csv` (lowered pool). A [**transient** seepage solve](../seep/transient.md) is the natural third source: a transient run *is* a time-history of pore-pressure fields as the reservoir is lowered on a schedule, so the two stage fields are simply two frames read from that one history rather than two independent steady solves assembled by hand.

### Stage times on the tseep sheet

A transient run is driven by the [**tseep** sheet](../usage/input_template.md#worksheet-tseep), which carries two optional control times, **stage_1** and **stage_2** (with `stage_1 < stage_2`). They tag the two frames the drawdown analysis will use — for instance the full-reservoir steady state at `stage_1` and the drawn-down state at `stage_2`. Both times are forced into the transient solver's saved-frame schedule, so each is a *computed* frame, never interpolated between steps. Set both or neither; setting one without the other is an error, as is a `stage_1` at or after `stage_2`.

### In-memory staging

Once the transient run is solved, `stage_transient_for_drawdown(slope_data, solution)` pulls the frames at `stage_1` and `stage_2` **in memory** and writes their pore-pressure fields into `slope_data['seep_u']` and `slope_data['seep_u2']` — exactly the structures the classic two-file path produces. No intermediate `seep.csv` / `seep2.csv` files are written; the two stage fields go straight into the structures the three-stage machinery already consumes, and `rapid_drawdown` then runs unchanged.

**Resolution order.** The classic two-file path remains fully supported. When `{base}_seep.csv` and `{base}_seep2.csv` sit next to the input workbook, `load_slope_data` reads their `u` columns into `seep_u` and `seep_u2` at load time. Calling `stage_transient_for_drawdown` afterward overwrites those two fields with the staged transient frames, so a transient solution carrying stage times takes precedence over the classic files. A model with no stage times cannot be staged this way — the call requires both `stage_1` and `stage_2` — so it falls back to the classic two-file path.

### When transient staging is preferable

Two independent steady solves treat the lowered-pool field as a new steady state at the new level: they assume the reservoir dropped and pore pressures then fully re-equilibrated. A transient solve is preferable whenever that limiting assumption is too crude:

- **Real drawdown schedules.** A reservoir is lowered over hours or days, not instantaneously. A transient run honors the actual lowering rate and the elapsed time, so the Stage 2 field reflects how far pore-pressure dissipation has *actually* progressed — the same physics the time factor $T$ at the [top of this page](#when-does-rapid-drawdown-apply) estimates, now resolved directly by the solve.
- **Partial drawdown.** The pool need not be lowered all the way. `stage_2` can be tagged at any intermediate reservoir level, and the staged field is the true partially-drawn-down state.
- **Intermediate times.** Because every saved frame is a computed state, `stage_2` can be tagged at any time along the drawdown to ask how the factor of safety evolves as the pool falls and pressures bleed off, rather than only at the fully-dissipated endpoint.

### See also

The [transient-seepage theory page](../seep/transient.md) documents the formulation, storage properties, and time-stepping behind these frames; its [Rapid drawdown staging](../seep/transient.md#rapid-drawdown-staging) section covers the coupling from the seepage side and shows the staging call in a runnable example. A worked dam-drawdown example that carries a lowering schedule through to a rapid-drawdown factor of safety is being added to the samples; until it lands, that staging section is the end-to-end reference.