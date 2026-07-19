# Reliability Analysis

The foundation of the limit equilibrium method is the factor of safety (FoS), which is defined as the ratio of the resisting forces to the driving forces. This is a simple and useful metric, but it does not provide a complete picture of the reliability of the system. For example, if a structure has a factor of safety = 1.5, we have an intuitive sense of the stability, but ultimately, how safe is the structure? How likely is it to fail? An alternative approach to slope stability analysis is to calculate the reliability of a slope, which is defined as the probability that the slope will not fail. It is related to the probability of failure ($P_f$) as follows:

>>$R = 1 - P_f$

where:

>>$R$ = the reliability of the slope<br>
> $P_f$ = the probability of failure


!!! note "FEM-based reliability"
    The same method is also available with the finite-element SSRM solver, which
    computes each factor of safety by strength reduction instead of a limit-
    equilibrium search — see [Reliability Analysis (FEM)](../fem/reliability.md).
    The theory below (parameter uncertainty, the Taylor Series method, and the
    reliability equations) applies to both.

In reliability analysis, we consider the uncertainties in the parameters that affect the stability of the slope. These uncertainties can arise from various sources, such as variations in soil properties, loading conditions, and environmental factors. By incorporating these uncertainties into our analysis, we can obtain a more comprehensive understanding of the slope's stability. The **xslope** package provides a function to calculate the reliability of a slope using the limit equilibrium method. This function takes into account the uncertainties in the soil properties ($\gamma$, $c$, $\phi$, etc.) and provides a probability of failure (Pf) and reliability (R) for the slope. It can be used with any of the limit equilibrium methods implemented in the package, such as Bishop's method, Janbu's method, or Spencer's method. It can also be combined with a rapid drawdown analysis.

## Parameter Uncertainty

In order to perform reliability analysis, we need to consider the uncertainties in the parameters that affect the stability of the slope. These parameters can include soil properties such as unit weight ($\gamma$), cohesion ($c$), and angle of internal friction ($\phi$). The uncertainties in these parameters can be represented using probability distributions, such as normal or lognormal distributions. Uncertainty is typically represented by a standard deviation or a coefficient of variation (COV), which is the ratio of the standard deviation to the mean. The standard deviation provides a measure of the spread of the parameter values, while the COV provides a measure of the relative variability of the parameter. The standard deviation is defined as:

>>$\sigma = \sqrt{\frac{1}{n-1} \sum_{i=1}^{n} (x_i - \bar{x})^2}$

where:

>>$\sigma$ = standard deviation<br>
> $x_i$ = individual parameter value<br>
> $\bar{x}$ = mean of the parameter values<br>
> $n$ = number of parameter values


The COV is a dimensionless quantity that is defined as:

>>$COV = \dfrac{\sigma}{\bar{x}}$

where:

>>$COV$ = coefficient of variation<br>
> $\sigma$ = standard deviation<br>
> $\bar{x}$ = mean of the parameter values

To calculate a proper standard deviation or COV for a parameter, we need to have a set of representative values for that parameter. This can be done by performing laboratory tests on soil samples, conducting field investigations, or using empirical correlations based on previous studies. In practice, we often do not have a large number of values for a parameter, so we may need to use expert judgment or statistical methods to estimate the standard deviation or COV. For example, if we have a limited number of soil samples, we can use the range of the values to estimate the standard deviation. Alternatively, we can use a default value based on previous studies or guidelines. To use a range of values, we can use the $3\sigma$ rule, which states that approximately 99.7% of the values in a normal distribution lie within three standard deviations of the mean. Therefore, we can estimate the standard deviation as:

>>$\sigma = \dfrac{(x_{max} - x_{min})}{6}$

where:

>>$\sigma$ = standard deviation<br>
> $x_{max}$ = maximum value of the parameter<br>
> $x_{min}$ = minimum value of the parameter

Studies have shown that even experts often under-estimate max and min values (estimated range is typically too small). Thus, it is conservative to estimate the standard deviation as:

>>$\sigma = \dfrac{(x_{max} - x_{min})}{4}$

where:

>>$\sigma$ = standard deviation<br>
> $x_{max}$ = maximum value of the parameter<br>
> $x_{min}$ = minimum value of the parameter

We can also look at the literature to find typical values for the standard deviation or COV of a parameter. Here are some typical values for the  COV of common soil parameters:

| Parameter                                  | Typical COV (%) |
|--------------------------------------------|-----------------|
| Soil Unit Weight ($\gamma$)                | 3 - 10%         |
| Effective Stress Cohesion Intercept ($c'$) | 15 - 100%       |
| Effective Stress Friction Angle ($\phi'$)  | 3 - 20%         |
| Undrained Shear Strength ($S_u$)           | 15 - 50%        |

## Reliability Equation

Overall reliability can be determined once we know the following values:

>>$F_{MLV}$ = The factor of safety based on the most likely values of the parameters

>>$COV_{F}$ = The coefficient of variation of the factor of safety, which is a measure of the uncertainty in the factor of safety due to the uncertainties in the parameters

These values can be used to compute the lognormal reliability index ($\beta$) using the following equation:

>>$\beta_{LN} = \dfrac{\ln \left( \dfrac{F_{MLV}}{\sqrt{1 + COV_F^2}} \right)}{\sqrt{\ln (1 + COV_F^2)}}$

We can then calculate the reliability ($R$) using a normal distribution function. In Excel, we can use the `NORMSDIST` function to calculate the reliability:

```excel
= NORMSDIST(beta_LN)
```
In Python, this can be done using the `scipy.stats.norm` module:

```python
from scipy.stats import norm
def calculate_reliability(F_MLV, COV_F):
    beta_LN = (np.log(F_MLV / np.sqrt(1 + COV_F**2)) /
               np.sqrt(np.log(1 + COV_F**2)))
    R = norm.cdf(beta_LN)
    return R
```

But first, we have to calculate the factor of safety based on the most likely values of the parameters ($F_{MLV}$) and the coefficient of variation of the factor of safety ($COV_F$). Calculating the $F_{MLV}$ is straightforward, as it is simply the factor of safety calculated using the most likely values of the parameters. This can be done using any of the limit equilibrium methods implemented in the **xslope** package, such as Bishop's method, Janbu's method, or Spencer's method.

To calculate the $COV_F$, there are two common approaches to calculate these values: the Monte Carlo method and the Taylor Series Probability Method (TSPM).

### Monte Carlo Method

The Monte Carlo method is a statistical technique that uses random sampling to estimate the probability of failure. It involves generating a large number of random samples for the uncertain parameters. Then we do the following steps:

1. Generate N random samples for the uncertain parameters (e.g., soil unit weight, cohesion, angle of internal friction) based on their probability distributions.
2. Combine random values for each uncertain parameter to create N "model instances". Each of these instances represents a unique combination of parameter values and is considered to be equally likely.
3. Calculate the factor of safety for each model instance using a limit equilibrium method (e.g., Bishop's method, Janbu's method, or Spencer's method).
4. Calculate the coefficient of variation of the factor of safety values obtained from the model instances.

While the Monte Carlo method is a powerful and flexible approach, it can be computationally expensive for slope stability problems as it requires a large number of model runs to obtain accurate results. The accuracy of the Monte Carlo method depends on the number of samples generated, and typically, thousands of model instances are needed to achieve a reliable estimate of the probability of failure.

#### Monte Carlo in xslope

The `reliability_mc` function runs a Monte Carlo campaign directly from the same
inputs as the Taylor series: the material most-likely values and the standard
deviations already in the mat sheet. Nothing new is required — *"we just need a
number of runs, because we already have the standard deviations and the most-likely
values."*

**Sampling model.** Each uncertain parameter is drawn independently from a
distribution whose **mean is the material's base value** (the MLV in the mat sheet)
and whose **standard deviation is the matching s(·) column** — s(g) for $\gamma$,
s(c) for $c$, s(f) for $\phi$, s(c/p) for $c_p$. The default distribution is
**normal**, the same interpretation the Taylor series places on those columns; a
`distribution='lognormal'` option draws lognormal samples matched to the same mean
and standard deviation. Every draw is truncated at its physical floor — a strength
or unit weight cannot go negative — which is the one modelling assumption beyond the
raw mean and sigma. That truncation is exactly the $\phi \ge 0$ bound that governs
high-COV problems (see the VP34 worked example below); at ordinary COVs it never
activates and Monte Carlo and the Taylor series agree closely.

**Number of runs and convergence.** The default is **10,000 samples**. The sample
mean and standard deviation of the factor of safety converge quickly (a few thousand
samples), but the empirical probability of failure — a count of the tail — converges
more slowly, so a small $P_f$ needs more samples to resolve. Increase `n_samples`
until the reported $P_f$ is stable to the precision you need.

**Deterministic (seeded) sampling.** The random-number generator is seeded from a
fixed constant (never the clock), so a given input file reproduces the same $\beta$
and $P_f$ on every run — the results are regression-lockable. Pass a different
`rng_seed` to inspect sampling scatter.

**Outputs.** The function returns the sample **mean factor of safety**, its standard
deviation $\sigma_F$ and $COV_F$, the reliability index in **both conventions** — a
normal index $\beta_N = (\bar F - 1)/\sigma_F$ and the lognormal $\beta_{LN}$ from
the sample moments (same formula as the Taylor series) — and the **empirical
probability of failure**, the fraction of realizations with $F < 1$, alongside the
distribution-fitted normal and lognormal $P_f$ for comparison.

**Limit-equilibrium only.** Monte Carlo is available on the limit-equilibrium
solvers only: a campaign of $10^4$ factor-of-safety evaluations is affordable with a
limit-equilibrium solve but not with the finite-element SSRM, so FEM reliability
stays on the Taylor series (1 + 2N solves — see
[Reliability Analysis (FEM)](../fem/reliability.md)).

#### Surface treatment

The slip surface is a **decision variable, not a random variable** — its geometry is
therefore **never randomized**. Two treatments are meaningful:

- **Fixed-surface Monte Carlo (the default).** Every realization is evaluated on the
  *same* surface: the prescribed surface when `search=False`, or the
  most-likely-values critical surface when `search=True` (found once, then held
  fixed). This is the analogue of Slide2's **Global Minimum** probabilistic method,
  and it is what every published probabilistic benchmark did — Duncan's estimated
  LASH surface (VP29), Chowdhury & Xu's printed circles (VP28), Wolff & Harr's
  prescribed noncircular surface (VP34) — so it is the correct like-for-like basis
  for comparison. The honest caveat: a fixed-surface $P_f$ *understates* the
  slope-system probability of failure, because for some sampled parameter sets the
  true critical surface lies elsewhere.
- **Per-realization minimum (a follow-up mode).** The system-level analogue of
  Slide2's *Overall Slope* method re-evaluates, for each realization, the **ensemble
  of candidate surfaces already found by the deterministic search** and takes the
  minimum — *not* a fresh search per realization, which would be prohibitive. This
  mode is a documented follow-up and is not yet shipped; the fixed-surface campaign
  above is what the corpus comparisons use.

The Taylor-series side has the same distinction, already exercised in the corpus: on
the Cannon Dam benchmark (VP35, Hassan & Wolff 1999) the **surface of minimum
reliability index is not the surface of minimum factor of safety**, so a design
screened on FS alone examines the wrong surface. Randomizing geometry is never the
answer; searching on $\beta$ (or on FS) is.

#### When to use Monte Carlo versus the Taylor series

For ordinary parameter scatter the two methods agree, and the Taylor series is far
cheaper (1 + 2N solves versus $10^4$), so it is the default. Reach for Monte Carlo
when the first-order Taylor assumptions break down:

- **Large coefficients of variation.** When a standard deviation approaches or
  exceeds its mean, the Taylor series cannot evaluate $F(\text{MLV} - \sigma)$
  without a negative parameter and declines the analysis. Monte Carlo handles it by
  truncating the draw at zero. The worked example is **VP34** (Clarence Cannon Dam),
  whose Phase I fill has a friction-angle COV of 124%: the Taylor series declines,
  while Monte Carlo returns a probability of failure inside the range the published
  studies span.
- **Skewed or non-linear factor-of-safety response**, where the first-order
  (central-difference) slope is a poor summary of the true $F(\mathbf{x})$ surface.
- **Tail probabilities**, where the empirical $P_f$ from the samples is wanted
  directly rather than inferred from a fitted lognormal.

### Taylor Series Probability Method (TSPM)

The Taylor Series Probability Method (TSPM) is a more efficient approach for calculating the coefficient of variation of the factor of safety as it requires a much smaller set of model runs. It uses the first-order Taylor series expansion to approximate the factor of safety as a function of the uncertain parameters. The TSPM can be summarized in the following steps:

1. Determine the standard deviation ($\sigma_i$) for each uncertain parameter using the guidelines described above.
2. Find $F_i^+$ and $F_i^-$, for each parameter where:

>>$F_i^+$ = the factor of safety calculated using the parameter value = $MLV + \sigma_i$ with all other parameters held at the most likely value<br>
> $F_i^-$ = the factor of safety calculated using the parameter value = $MLV - \sigma_i$ with all other parameters held at the most likely value

3. Compute $\Delta F_i = |F_i^+ - F_i^-|$ for each parameter.
4. Compute $\sigma_F$ and $COV_F$ using the following equations:


>>$\sigma_F = \sqrt{\left(\dfrac{\Delta F_1}{2} \right)^2 + \left(\dfrac{\Delta F_2}{2} \right)^2 + \ldots + \left(\dfrac{\Delta F_n}{2} \right)^2}$

>>$COV_F = \dfrac{\sigma_F}{F_{MLV}}$

## Data Input

To perform reliability analysis using the **xslope** package, we simply need to provide standard deviations for the uncertain parameters in the input data. This is done in the Materials table of the input data file. The main values of the parameters in the table are treated as the most likely values. We can then call the `reliability` function to perform the analysis. The function will automatically calculate the factor of safety based on the most likely values ($F_{MLV}$) of the parameters using an automated search. It will then perturb each parameter by the standard deviation using the Taylor Series Method described above to calculate the coefficient of variation of the factor of safety ($COV_F$). Finally, it will compute the reliability of the slope based on the calculated values.

The strength parameters that are perturbed depend on each material's strength model: for Mohr-Coulomb (`mc`) materials the cohesion $c$ and friction angle $\phi$ are perturbed, while for the depth-varying undrained (`cp`) model the cohesion $c$ and the rate $c_p$ are perturbed. The unit weight $\gamma$ is perturbed in both cases. If a standard deviation exceeds its mean — so that mean $-\sigma$ would be negative — the Taylor-series analysis stops with an error, since a negative strength parameter is non-physical. This is the boundary of the Taylor-series method's domain; when a parameter's COV is that large, use the Monte Carlo function instead (`reliability_mc`, described in [Monte Carlo in xslope](#monte-carlo-in-xslope)), which handles the negative draw by truncating it at zero.

The Monte Carlo campaign is called the same way through `reliability_mc`, which reads the identical MLV and standard-deviation inputs and adds only the sampling controls (`n_samples`, `distribution`, `rng_seed`). It is a limit-equilibrium path only, for the compute reason noted above.

One of the arguments to the function is `method`, which specifies the limit equilibrium method to be used for the analysis. The available methods are 'bishop', 'janbu', and 'spencer'. The function will return the probability of failure ($P_f$) and reliability ($R$) of the slope. Either a circular or non-circular slope can be analyzed. If a circular slope is analyzed, care should be taken to select a set of starting circles in the circles table of the input file to ensure that the automated search finds the global minimum factor of safety for each analysis. It is good practice to include a circle that touches the bottom each material zone and perhaps a circle that passes through the toe of the slope. 