# Reliability Analysis (FEM)

Reliability analysis quantifies the probability that a slope will *not* fail by
propagating the uncertainty in the soil parameters through to the factor of
safety. This page describes the **finite-element** reliability analysis, which is
identical in method to the [limit-equilibrium reliability analysis](../lem/reliability.md)
but computes each factor of safety with the finite-element Shear Strength
Reduction Method (SSRM) instead of a limit-equilibrium search.

Because the method, the parameter-uncertainty guidance (standard deviation and
COV estimation, typical COV values), and the reliability equations are the same,
they are **not repeated here** — see the [LEM Reliability Analysis](../lem/reliability.md)
page for the full theory. This page covers only what is specific to the FEM
version.

## Method

The analysis uses the same **Taylor Series Probability Method (TSPM)**:

1. **$F_{MLV}$** — the factor of safety at the most-likely parameter values,
   computed with `solve_ssrm`.
2. For each uncertain strength parameter, compute **$F_i^+$** and **$F_i^-$** with
   that parameter set to $MLV \pm \sigma_i$ (all others held at their most-likely
   value), then $\Delta F_i = |F_i^+ - F_i^-|$.
3. Combine: $\sigma_F = \sqrt{\sum_i (\Delta F_i / 2)^2}$, $COV_F = \sigma_F / F_{MLV}$,
   and finally the lognormal reliability index $\beta_{LN}$ and the reliability
   $R = \Phi(\beta_{LN})$.

The analysis runs $1 + 2N$ SSRM solves (one at the most-likely values plus a
$F^+$/$F^-$ pair for each of the $N$ uncertain parameters), all on a **single
shared mesh** — only the material-to-element mapping is rebuilt per perturbation.
The [auto-expanding SSRM bracket](overview.md) is what makes this practical:
each perturbation shifts the factor of safety, and the bracket adjusts itself so a
fixed `F_min`/`F_max` does not have to bracket every perturbed case in advance.

## Which parameters are perturbed

The FEM reliability perturbs the **same strength parameters as the LEM
reliability** — the cohesion $c$ and friction angle $\phi$ for a Mohr-Coulomb
(`mc`) material, or $c$ and the rate $c_p$ for the depth-varying undrained (`cp`)
model, plus the unit weight $\gamma$ for both — using the standard deviations
`sigma_c`, `sigma_phi` / `sigma_cp`, and `sigma_gamma` from the materials table.
Using the same standard deviations keeps the FEM and LEM reliability results
directly comparable.

### Why the elastic parameters E and ν are excluded

Unlike a limit-equilibrium analysis, a finite-element analysis also takes the
elastic modulus $E$ and Poisson's ratio $\nu$ as inputs. It is conventional
wisdom that these affect the computed **displacements** but have little effect on
the **factor of safety**, which is governed by the strength. We confirmed this
directly on Griffiths & Lane (1999) Example 1 ($c/\gamma H = 0.05$, $\phi = 20°$),
holding the strength fixed and varying $E$ and $\nu$ across their full plausible
ranges:

| Variation                       | FS     | ΔFS vs. base |
|---------------------------------|--------|--------------|
| Base ($E$ = 700,000, $\nu$ = 0.30) | 1.394  | —            |
| $E \times 0.5$ (350,000)        | 1.394  | 0.000        |
| $E \times 2.0$ (1,400,000)      | 1.394  | 0.000        |
| $\nu$ = 0.20                    | 1.406  | +0.013       |
| $\nu$ = 0.40                    | 1.381  | −0.013       |

**$E$ has no effect on the factor of safety** — halving *and* doubling it give an
identical FS, because SSRM converges to the same critical strength-reduction
factor regardless of the elastic stiffness. **$\nu$ changes the FS by only ~1%
across its entire plausible range** (0.2–0.4), which is at the resolution of the
bisection itself; at a realistic $\sigma_\nu$ (Poisson's ratio rarely varies by
more than about ±0.05) its contribution to $COV_F$ is negligible.

For these reasons $E$ and $\nu$ are treated as deterministic and are **not**
perturbed, and the materials table has no `sigma_E` / `sigma_nu` columns. The same
strength-parameter standard deviations drive both the LEM and the FEM reliability
analyses.

## Data Input

Provide the standard deviations for the uncertain strength parameters in the
Materials table exactly as for the [LEM reliability analysis](../lem/reliability.md#data-input)
— the main parameter values are the most-likely values, and `sigma_c`,
`sigma_phi` / `sigma_cp`, and `sigma_gamma` give their standard deviations. At
least one non-zero standard deviation is required. As in the LEM case, a standard
deviation may not exceed its mean (a negative $MLV - \sigma$ is non-physical and
the analysis stops with an error).

## Usage

```python
from xslope.fileio import load_slope_data
from xslope.advanced import reliability_fem

slope_data = load_slope_data("my_slope.xlsx")
success, result = reliability_fem(
    slope_data,
    element_type="quad8",   # quadratic elements (tri6/quad8/quad9); see Overview
    target_size=3.5,        # omit to auto-size from the domain width
    F_min=1.0, F_max=2.0,   # starting bracket; auto-expands if the guess is off
)
if success:
    print(f"F_MLV = {result['F_MLV']:.3f}")
    print(f"COV_F = {result['COV_F']:.3f}")
    print(f"Reliability = {result['reliability']*100:.2f}%")
    print(f"Probability of failure = {result['prob_failure']*100:.2f}%")
```

`reliability_fem` returns the same reliability keys as the LEM `reliability`
function (`F_MLV`, `sigma_F`, `COV_F`, `beta_ln`, `reliability`, `prob_failure`,
`param_info`), plus `mlv_solution` (the SSRM result at the most-likely values) and
the `mesh` used for every trial. Building the mesh is the same as for a single
SSRM run: it is built here from the geometry (or an attached `slope_data['mesh']`
/ a `mesh=` argument is reused), and all $1 + 2N$ trials share it so the only
thing that changes between trials is the perturbed material property.
