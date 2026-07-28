# Reliability Analysis (FEM)

Reliability analysis quantifies the probability that a slope will *not* fail by
propagating the uncertainty in the soil parameters through to the factor of
safety. This page describes the **finite-element** reliability analysis, which is
identical in method to the [Taylor Series Probability Method](taylor.md) used for
limit equilibrium, but computes each factor of safety with the finite-element Shear
Strength Reduction Method (SSRM) instead of a limit-equilibrium search.

Because the method, the parameter-uncertainty guidance (standard deviation and
COV estimation, typical COV values), and the reliability equations are the same,
they are **not repeated here** — see the [Reliability Analysis](index.md) overview
and the [Taylor Series Probability Method](taylor.md) page for the full theory.
This page covers only what is specific to the FEM version.

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
The [auto-expanding SSRM bracket](../fem/overview.md) is what makes this practical:
each perturbation shifts the factor of safety, and the bracket adjusts itself so a
fixed `F_min`/`F_max` does not have to bracket every perturbed case in advance.

!!! note "The `F_min`/`F_max` bracket is only used to find $F_{MLV}$"
    The bracket you supply is used **only for the initial most-likely-values
    solve**. Every subsequent perturbation solve brackets *automatically* — a
    window centred on $F_{MLV}$ (auto-expanding if a perturbation lands outside
    it) — so you do not size the bracket for the perturbed cases. This also means
    the bracket has essentially no effect on the reported reliability once the
    bisection tolerance is tight (see below).

### Numerical precision and reproducibility

A plain SSRM factor of safety is the midpoint of the final bisection band, so it
carries a $\pm\text{tolerance}/2$ imprecision *whose exact value depends on the
starting `F_min`/`F_max`* (different brackets subdivide the axis on different
grids). The reliability index amplifies this:
$\dfrac{d\beta}{dF} = \dfrac{1}{F\sqrt{\ln(1+COV_F^2)}}$, which is $\approx 9$ at
$COV_F = 0.1$ — so a small bracket-dependent change in $F_{MLV}$ produces a
proportionally larger change in $\beta$, enough to flip the last shown digit of the
reliability between two bracket choices.

To remove that entirely, `reliability_fem` runs each SSRM on a **fixed global grid**
(`grid = tolerance`, default 0.001; see [`solve_ssrm`](../fem/overview.md)). Instead of
halving your bracket, it locates the single global grid cell that straddles the
failure threshold — a fact of the slope and mesh, not of the bracket — so *every*
starting bracket lands in the same cell. The result is **identical to every decimal
regardless of `F_min`/`F_max`**, at the same cost (~log₂ solves) and precision
(the grid step) as ordinary bisection.

Two things still legitimately change the result and are **not** numerical noise:
the **mesh** (a finer or different-element mesh gives a slightly different factor
of safety, as with any FE analysis — use a converged mesh and report which one),
and genuinely larger parameter uncertainty (a higher $COV_F$ both lowers the
reliability and *reduces* the $\beta$ sensitivity above).

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
directly on [Griffiths & Lane (1999)](https://doi.org/10.1680/geot.1999.49.3.387) Example 1 ($c/\gamma H = 0.05$, $\phi = 20°$),
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
Materials table exactly as for [Taylor series reliability](taylor.md#data-input)
— the main parameter values are the most-likely values, and `sigma_c`,
`sigma_phi` / `sigma_cp`, and `sigma_gamma` give their standard deviations. At
least one non-zero standard deviation is required. As in the LEM case, a standard
deviation may not exceed its mean (a negative $MLV - \sigma$ is non-physical and
the analysis stops with an error).

## Usage

In **XSLOPE Studio**, build a mesh, then choose **Run FEM → Analysis: Reliability
(SSRM)**. It uses the `F min`/`F max` bracket from the dialog (auto-expanding) and
the material standard deviations from the mat sheet; the bisection tolerance,
however, is set internally to the tight reliability default (not the dialog's
single-run *Tolerance* field — see [Numerical precision](#numerical-precision-and-reproducibility)).
The FEM Results view shows the deformation at the most-likely values, the
reliability summary appears in a dialog and the status bar, and the full
per-parameter ΔF table is written to the Log pane. See
[Reliability analysis](../studio/analysis.md#reliability-analysis) for the dialog and
the engine comparison shared with the [Taylor series](taylor.md) and
[Monte Carlo](monte_carlo.md) pages.

From Python:

```python
from xslope.fileio import load_slope_data
from xslope.advanced import reliability_fem

slope_data = load_slope_data("my_slope.xlsx")
success, result = reliability_fem(
    slope_data,
    element_type="quad8",   # quadratic elements (tri6/quad8/quad9); see Overview
    target_size=3.5,        # omit to auto-size from the domain width
    F_min=1.0, F_max=2.0,   # starting bracket; auto-expands if the guess is off
    # tolerance defaults to a tight 0.001 here (see Numerical precision above)
)
# reliability_fem already prints the per-parameter ΔF table and the summary
# (F_MLV, COV_F, beta, reliability, Pf); all of these are also in `result`.
```

`reliability_fem` returns the same reliability keys as the LEM `reliability`
function (`F_MLV`, `sigma_F`, `COV_F`, `beta_ln`, `reliability`, `prob_failure`,
`param_info`), plus `mlv_solution` (the SSRM result at the most-likely values) and
the `mesh` used for every trial. Only the most-likely-values solve — the one
result that is plotted — captures the at-failure mechanism for the
deformation/vector/strain figures; the perturbation solves report only their
factor of safety. Building the mesh is the same as for a single
SSRM run: it is built here from the geometry (or an attached `slope_data['mesh']`
/ a `mesh=` argument is reused), and all $1 + 2N$ trials share it so the only
thing that changes between trials is the perturbed material property. Because the
factor of safety — and hence the reliability — depends on the mesh, report the
element type and size used for a given result.
