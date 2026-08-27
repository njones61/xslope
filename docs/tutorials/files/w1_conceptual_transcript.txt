# w1_conceptual — recorded assistant session

Captured by `tools/assistant_sessions.py` against the live Studio assistant, offscreen.

- Project: `docs/tutorials/files/xslope_reinforced_slope.xlsx`
- Provider / model: `anthropic` / `claude-opus-5`
- Confirm-before-running: off (unattended capture)

## Turn 1 — How does a reliability analysis work in XSLOPE?

````text
You: How does a reliability analysis work in XSLOPE?

Ran code:
    for i,m in enumerate(slope_data['materials'],1):
        print(i, m['name'], {k:m.get(k) for k in ('sigma_gamma','sigma_c','sigma_phi','sigma_cp','sigma_d','sigma_psi')})

Output:
    1 shell {'sigma_gamma': 0.0, 'sigma_c': 0.0, 'sigma_phi': 0.0, 'sigma_cp': 0.0, 'sigma_d': 0.0, 'sigma_psi': 0.0}
    2 base {'sigma_gamma': 0.0, 'sigma_c': 0.0, 'sigma_phi': 0.0, 'sigma_cp': 0.0, 'sigma_d': 0.0, 'sigma_psi': 0.0}

Assistant: ## The idea

A reliability analysis treats the strength (and unit weight) inputs as **random variables** instead of single numbers, propagates that scatter through the LEM solver, and reports the distribution of FS rather than one value. The outputs are:

- **E[FS]**, **σ_FS**
- **β = (E[FS] − 1) / σ_FS** — the reliability index, how many standard deviations the mean FS sits above failure
- **P_f** — probability of failure, `P(FS < 1)`, from a normal (or lognormal) fit

Inputs are the **σ columns on each material**: `sigma_gamma`, `sigma_c`, `sigma_phi`, `sigma_cp`, `sigma_d`, `sigma_psi`. A σ of 0 means that property is deterministic.

## The three engines

| Call | How it samples | Failure surface |
|:--|:--|:--|
| `reliability_taylor(method='bishop')` | Taylor-series (FOSM): 1 base solve + 2 per random variable, at mean ± σ; σ_FS from the summed partial derivatives | **Re-searched on every one of the 1 + 2N solves** — the critical surface migrates with the parameter values |
| `reliability_mc(method='bishop', n_samples=…, rng_seed=…)` | Monte Carlo: random draws from each variable's distribution | **Held FIXED** — never re-searched |
| `reliability_rs(method='bishop')` | Response surface: a fitted polynomial surrogate, then sampled cheaply | **Held FIXED** |

Front door: `reliability(method='bishop', engine='taylor'|'mc'|'rs')`. Each renders its own plot.

That surface difference is reportable, not cosmetic: Taylor's β accounts for the critical surface moving as strengths vary; MC and RS give β and P_f **conditional on one surface**, which is cheaper and usually slightly optimistic. Cost runs the other way — Taylor is 1 + 2N solves (cheap), MC is `n_samples` solves (expensive but distribution-free), RS buys MC-like sampling at a handful of real solves.

## Micro-example

With `E[FS] = 1.45` and `σ_FS = 0.18`:

```
β  = (1.45 − 1)/0.18 = 2.50
P_f = Φ(−2.50) ≈ 0.6%
```

## Your model right now

Both `shell` and `base` carry **σ = 0 on every property**, so a reliability run has nothing to randomize — it would return σ_FS = 0 and an infinite β. You'd need to fill in at least `sigma_c` and `sigma_phi` (typical COVs: c ≈ 20–40%, φ ≈ 5–15%, γ ≈ 3–7%) before any engine gives a meaningful answer.

Say the word and I'll set the σ values you want and run one — or read the derivations at `reliability/`, `reliability/taylor/`, and `reliability/monte_carlo/`.
````

Tokens: input 28,581 · cached_input 27,308 · output 1,380 · calls 2 · 19.8s

## Turn 2 — How do I decide standard deviations for a reliability analy…

````text
You: How do I decide standard deviations for a reliability analysis if I only have a few tests?

Assistant: ## Short answer

With a handful of tests you don't *estimate* σ statistically — you **bound** it and check that the answer isn't sensitive to which bound you chose. Three ways, in the order you should reach for them.

### 1. Published COV × your mean

The most defensible default. σ = COV · mean, with COV from the literature (Duncan 2000, *Factors of Safety and Reliability in Geotechnical Engineering*, ASCE JGGE):

| Property | Typical COV |
|:--|:--|
| γ (total unit weight) | 3–7% |
| c′ (effective cohesion) | 20–40% |
| φ′ (friction angle) | 5–15% |
| S_u (undrained strength) | 20–40% |
| tan φ′ | 5–15% |

Use the **high** end when your few tests came from scattered locations or different labs, the low end for a uniform compacted fill with controlled QC.

### 2. The three-sigma rule

Ask yourself the highest and lowest **conceivable** values (not "likely" — conceivable), then

```
σ ≈ (HCV − LCV) / 6
```

because ±3σ spans ~99.7% of a normal distribution. The known failure mode: engineers pick ranges that are **too narrow**, so Duncan recommends dividing by **4** instead of 6 unless you're confident in the extremes. With 3–5 tests, use /4.

### 3. Small-sample statistics, inflated

If you do compute `s` from n tests with the (n−1) denominator, remember that for n < 10 the sample range badly understates the population range. Treat `s` as a **lower bound** on σ, and carry the larger of `s` and the COV-based estimate.

## Three things that trip people up

**Separate the sources of scatter.** Total variability = inherent spatial variability + measurement/test error + statistical uncertainty in the mean. Only the *inherent* part is real strength variation; lab scatter isn't a property of the ground. With few tests you can't decompose it, so using total scatter is conservative — say so.

**Spatial averaging reduces the effective σ.** A slip surface 100 ft long mobilizes the *average* strength along it, not a point value. Where the scale of fluctuation is short compared to the surface, the effective σ on the surface is well below the point σ. Ignoring this is conservative, which is usually fine — but it's why a P_f of a few percent from point-scale σ is not a literal failure probability.

**Correlation matters and the Taylor engine assumes none.** c′ and φ′ from the same tests are typically *negatively* correlated (a steeper envelope gives a lower intercept). Treating them as independent overstates σ_FS — conservative again.

## The practical workflow in Studio

Don't agonize over every σ. Find out which one matters first, then bracket only that one:

```python
# 1. which variable dominates?
parametric_sweep(plot='tornado')

# 2. bracket the dominant σ: low and high estimates
for cov in (0.15, 0.40):
    for m in slope_data['materials']:
        m['sigma_c'] = cov * m['c']
    print(cov, reliability_taylor(method='bishop', plot=False))
```

If β stays on the same side of your target (commonly β ≥ 3, P_f ≤ 0.1% for a permanent slope) across the whole bracket, your σ choice didn't decide the answer and you can report it with the range stated. If β crosses the target inside the bracket, **that is your justification for more testing** — and it tells you exactly which property to test.

Note that `reliability_taylor` needs only the mean and σ, not a distribution shape, which is why it's the right engine when your data are thin. Monte Carlo makes you commit to a distribution you can't support with five tests.

Derivations: `reliability/taylor/` and `reliability/` at `https://xslope.readthedocs.io/en/latest/`.

Want me to put a first-pass set of σ values on your two materials — say COV 30% on c and 10% on φ — and run the bracket above?
````

Tokens: input 15,503 · cached_input 13,654 · output 1,935 · calls 1 · 29.4s

## Session total

- Tokens: input 44,084 · cached_input 40,962 · output 3,315 · calls 3
- Wall time: 57.6s
