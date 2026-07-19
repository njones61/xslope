# Reliability Analysis

The foundation of the limit equilibrium method is the factor of safety (FoS), which is defined as the ratio of the resisting forces to the driving forces. This is a simple and useful metric, but it does not provide a complete picture of the reliability of the system. For example, if a structure has a factor of safety = 1.5, we have an intuitive sense of the stability, but ultimately, how safe is the structure? How likely is it to fail? An alternative approach to slope stability analysis is to calculate the reliability of a slope, which is defined as the probability that the slope will not fail. It is related to the probability of failure ($P_f$) as follows:

>>$R = 1 - P_f$

where:

>>$R$ = the reliability of the slope<br>
> $P_f$ = the probability of failure

In reliability analysis, we consider the uncertainties in the parameters that affect the stability of the slope. These uncertainties can arise from various sources, such as variations in soil properties, loading conditions, and environmental factors. By incorporating these uncertainties into our analysis, we can obtain a more comprehensive understanding of the slope's stability. The **xslope** package provides a function to calculate the reliability of a slope using the limit equilibrium method. This function takes into account the uncertainties in the soil properties ($\gamma$, $c$, $\phi$, etc.) and provides a probability of failure (Pf) and reliability (R) for the slope. It can be used with any of the limit equilibrium methods implemented in the package, such as Bishop's method, Janbu's method, or Spencer's method. It can also be combined with a rapid drawdown analysis.

## The reliability family

All of the methods described in this section are reached through a single front door:

>>`reliability(slope_data, method, engine='taylor')`

The `engine` argument selects the analysis engine. `engine='taylor'` (the default)
runs the Taylor Series Probability Method; `engine='mc'` runs a Monte Carlo campaign.
Because the default is the Taylor series, an existing call such as
`reliability(slope_data, 'bishop')` keeps its exact meaning. The two engines are also
public directly, under their own honest names, so a script can call whichever it wants:

| Function | Method | Solver | Cost | Page |
|----------|--------|--------|------|------|
| `reliability_taylor` | Taylor Series Probability Method | limit equilibrium | 1 + 2N searches | [Taylor Series (TSPM)](taylor.md) |
| `reliability_mc`     | Monte Carlo                      | limit equilibrium | N (≈ 10⁴) evaluations | [Monte Carlo](monte_carlo.md) |
| `reliability_fem`    | Taylor Series Probability Method | finite-element SSRM | 1 + 2N solves | [Reliability Analysis (FEM)](fem.md) |

`reliability_fem` is the finite-element counterpart of the Taylor series — each factor
of safety comes from a strength-reduction (SSRM) solve rather than a limit-equilibrium
search (see [Reliability Analysis (FEM)](fem.md)). Monte Carlo is a
limit-equilibrium path only, for the compute reason discussed under
[Monte Carlo in xslope](monte_carlo.md#monte-carlo-in-xslope). All three engines read the *same*
inputs — the material most-likely values and the standard-deviation columns of the mat
sheet — so choosing an engine never changes the data you enter.

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

To calculate the $COV_F$, there are two common approaches to calculate these values: the [Monte Carlo method](monte_carlo.md) and the [Taylor Series Probability Method (TSPM)](taylor.md).

## When to use Monte Carlo versus the Taylor series

For ordinary parameter scatter the two methods agree, and the [Taylor series](taylor.md) is far
cheaper (1 + 2N solves versus $10^4$), so it is the default. Reach for [Monte Carlo](monte_carlo.md)
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

## How commercial software does it

xslope's two-engine split mirrors what the established geotechnical packages do, and
the split falls along the solver, not the software. The commercial **limit-equilibrium**
codes do reliability by direct **Monte Carlo**: Rocscience **Slide2** and GeoStudio
**SLOPE/W** both sample the strength and pore-pressure inputs from user-specified
distributions, evaluate the factor of safety of each realization, and report the
probability of failure as the fraction of realizations with FoS < 1 — the same
fixed-surface campaign as xslope's `reliability_mc`. Both offer **Latin Hypercube**
sampling as a more efficient alternative to plain Monte Carlo (SLOPE/W added it in
GeoStudio 2024.1), and Slide2 adds a machine-learning **response surface** to
accelerate the sampling and a **spatial-variability** (random-field) mode
([Slide2 Probabilistic Analysis](https://www.rocscience.com/help/slide2/documentation/slide-model/project-settings/statistics/probabilistic-analysis);
[Advancing Slope Stability Analysis with Probabilistic Methods in Slide2 and Slide3](https://www.rocscience.com/learning/advancing-slope-stability-analysis-with-probabilistic-methods-in-slide2-and-slide3);
[Stability Modeling with SLOPE/W (Seequent)](https://files.seequent.com/PDFs/Geostudio-Stability%20Modeling-Oct2022.pdf)).
On the **finite-element** side the economics change — 10⁴ full FE solves is
impractical — so the FE codes do *not* rely on brute-force sampling. Rocscience **RS2**
offers the Rosenblueth **point-estimate method** (two point estimates per variable, at
±1σ, explicitly recommended for problems with only a few random variables, ≈ 2–6), a
**response-surface** method, and Monte Carlo / Latin Hypercube over that surface
([RS2 Probabilistic Analysis](https://www.rocscience.com/help/rs2/documentation/rs2-model/project-settings/statistics/probabilistic-analysis)).
**PLAXIS**'s reliability module likewise provides Monte Carlo, Latin Hypercube and an
(alternative) point-estimate method with sensitivity analysis
([New Developments in PLAXIS: Material Point Method and Reliability Analysis, TU Delft](https://research.tudelft.nl/en/publications/new-developments-in-plaxis-material-point-method-and-reliability-)).
xslope's `reliability_fem` sits squarely in this point-estimate family: its 1 + 2N
Taylor-series perturbation is a two-point-per-variable estimate of the factor-of-safety
variance — the same economical strategy the FE vendors adopt in place of mass sampling.

!!! note "Future: response-surface Monte Carlo over the FEM solves"
    The 1 + 2N SSRM solves that `reliability_fem` already runs trace out a first-order
    response surface of the factor of safety in the uncertain parameters. A natural
    extension — recorded here, not yet built — is to fit that surface and run a Monte
    Carlo campaign *on the surface* (thousands of samples, no additional finite-element
    solves), yielding an empirical FEM probability of failure and a full FS
    distribution. That is exactly the response-surface strategy RS2 exposes; it would
    complement, not replace, the Taylor-series index.
