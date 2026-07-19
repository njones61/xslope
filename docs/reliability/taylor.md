# Taylor Series Probability Method (TSPM)

The Taylor Series Probability Method (TSPM) is a more efficient approach for calculating the coefficient of variation of the factor of safety as it requires a much smaller set of model runs. It uses the first-order Taylor series expansion to approximate the factor of safety as a function of the uncertain parameters. The TSPM can be summarized in the following steps:

1. Determine the standard deviation ($\sigma_i$) for each uncertain parameter using the guidelines described in [Parameter Uncertainty](index.md#parameter-uncertainty).
2. Find $F_i^+$ and $F_i^-$, for each parameter where:

>>$F_i^+$ = the factor of safety calculated using the parameter value = $MLV + \sigma_i$ with all other parameters held at the most likely value<br>
> $F_i^-$ = the factor of safety calculated using the parameter value = $MLV - \sigma_i$ with all other parameters held at the most likely value

3. Compute $\Delta F_i = |F_i^+ - F_i^-|$ for each parameter.
4. Compute $\sigma_F$ and $COV_F$ using the following equations:


>>$\sigma_F = \sqrt{\left(\dfrac{\Delta F_1}{2} \right)^2 + \left(\dfrac{\Delta F_2}{2} \right)^2 + \ldots + \left(\dfrac{\Delta F_n}{2} \right)^2}$

>>$COV_F = \dfrac{\sigma_F}{F_{MLV}}$

$\sigma_F$ and $COV_F$ then feed the lognormal reliability index $\beta_{LN}$ and the
reliability $R$ from the [Reliability Equation](index.md#reliability-equation). Each
parameter is evaluated at exactly two points — $MLV + \sigma_i$ and $MLV - \sigma_i$ —
so TSPM is, per parameter, a two-point (Rosenblueth-style) point-estimate method; see
[How commercial software does it](index.md#how-commercial-software-does-it) for how
that compares with the finite-element vendors' own point-estimate options.

## Data Input

To perform reliability analysis using the **xslope** package, we simply need to provide standard deviations for the uncertain parameters in the input data. This is done in the Materials table of the input data file. The main values of the parameters in the table are treated as the most likely values. We can then call the `reliability` function to perform the analysis. The function will automatically calculate the factor of safety based on the most likely values ($F_{MLV}$) of the parameters using an automated search. It will then perturb each parameter by the standard deviation using the Taylor Series Method described above to calculate the coefficient of variation of the factor of safety ($COV_F$). Finally, it will compute the reliability of the slope based on the calculated values.

The strength parameters that are perturbed depend on each material's strength model: for Mohr-Coulomb (`mc`) materials the cohesion $c$ and friction angle $\phi$ are perturbed, while for the depth-varying undrained (`cp`) model the cohesion $c$ and the rate $c_p$ are perturbed. The unit weight $\gamma$ is perturbed in both cases. If a standard deviation exceeds its mean — so that mean $-\sigma$ would be negative — the Taylor-series analysis stops with an error, since a negative strength parameter is non-physical. This is the boundary of the Taylor-series method's domain; when a parameter's COV is that large, use the Monte Carlo function instead (`reliability_mc`, described in [Monte Carlo in xslope](monte_carlo.md#monte-carlo-in-xslope)), which handles the negative draw by truncating it at zero.

The Monte Carlo campaign is called the same way through `reliability_mc`, which reads the identical MLV and standard-deviation inputs and adds only the sampling controls (`n_samples`, `distribution`, `rng_seed`). It is a limit-equilibrium path only, for the compute reason discussed in [Monte Carlo in xslope](monte_carlo.md#monte-carlo-in-xslope).

One of the arguments to the function is `method`, which specifies the limit equilibrium method to be used for the analysis. The available methods are 'bishop', 'janbu', and 'spencer'. The function will return the probability of failure ($P_f$) and reliability ($R$) of the slope. Either a circular or non-circular slope can be analyzed. If a circular slope is analyzed, care should be taken to select a set of starting circles in the circles table of the input file to ensure that the automated search finds the global minimum factor of safety for each analysis. It is good practice to include a circle that touches the bottom each material zone and perhaps a circle that passes through the toe of the slope.
