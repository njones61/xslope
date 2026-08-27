# w1_conceptual — recorded assistant session

Captured by `tools/assistant_sessions.py` against the live Studio assistant, offscreen.

- Project: `docs/tutorials/files/xslope_reinforced_slope.xlsx`
- Provider / model: `anthropic` / `claude-opus-5`
- Confirm-before-running: off (unattended capture)

## Turn 1 — How does a reliability analysis work in XSLOPE?

````text
You: How does a reliability analysis work in XSLOPE?

Ran code:
    for i, m in enumerate(slope_data['materials'], 1):
        print(i, m['name'], {k: m.get(k) for k in ('sigma_gamma','sigma_c','sigma_phi','sigma_cp','sigma_d','sigma_psi')})

Output:
    1 shell {'sigma_gamma': 0.0, 'sigma_c': 0.0, 'sigma_phi': 0.0, 'sigma_cp': 0.0, 'sigma_d': 0.0, 'sigma_psi': 0.0}
    2 base {'sigma_gamma': 0.0, 'sigma_c': 0.0, 'sigma_phi': 0.0, 'sigma_cp': 0.0, 'sigma_d': 0.0, 'sigma_psi': 0.0}

Assistant: ## The idea

You give each uncertain material property a **standard deviation** alongside its mean, and XSLOPE propagates that scatter through the LEM solver to get a distribution of FS rather than a single number. The outputs are:

- **μ_FS** — mean factor of safety
- **σ_FS** — standard deviation of FS
- **β = (μ_FS − 1) / σ_FS** — the reliability index, how many standard deviations the mean sits above failure
- **p_f** — probability of failure, `Φ(−β)` under a normal assumption (lognormal variants shift this)

β is the number that matters: FS = 1.5 with tight strength data and FS = 1.5 with wildly scattered data are not the same slope.

## Where the inputs live

On the materials, as the `sigma_*` fields — `sigma_gamma`, `sigma_c`, `sigma_phi`, `sigma_cp`, `sigma_d`, `sigma_psi`. A σ of 0 means "treated as deterministic," and **that is where your model stands right now** — both `shell` and `base` have every σ at 0.0, so a reliability run would report σ_FS = 0 and β = ∞. You'd need to set at least `sigma_c` and `sigma_phi` before any of the engines below say anything.

## The three engines

`reliability(method='bishop', engine=...)` is the front door; each engine also has its own entry point.

**Taylor series** (`reliability_taylor`) — the classic Duncan MSE method. For each variable it re-solves at mean ± σ, takes ΔFS/2 as that variable's contribution, and combines them in quadrature:

```
σ_FS = sqrt( Σ (ΔFS_i / 2)² )
```

Cost is 2n + 1 solves for n variables — cheap, and it hands you a per-variable breakdown showing which property dominates the uncertainty. It assumes FS is roughly linear in each variable over ±σ.

**Monte Carlo** (`reliability_mc(n_samples=..., rng_seed=...)`) — samples the full parameter set N times and solves each realization, so p_f is counted from the empirical distribution rather than assumed normal. No linearity assumption, but N solves; use `rng_seed` so the run is reproducible.

**Response surface** (`reliability_rs`) — fits a polynomial surrogate to a small designed set of real solves, then samples the surrogate millions of times. Monte Carlo accuracy at Taylor-ish cost, valid as long as the surface fits.

Each renders its own plot.

## A micro-example

Say Taylor gives μ_FS = 1.62, σ_FS = 0.21. Then β = (1.62 − 1)/0.21 = 2.95, and p_f = Φ(−2.95) ≈ 1.6×10⁻³ — about 1 in 600. Drop σ_FS to 0.12 with better site data and β jumps to 5.2, p_f ≈ 1×10⁻⁷, with the mean FS completely unchanged.

## Caveats worth stating

- The σ's are **your** estimate of parameter uncertainty; they don't include model error, and p_f is a notional probability, not an actuarial one.
- Variables are treated as independent unless you say otherwise — real c and φ from a common test set are usually negatively correlated, which makes independent sampling conservative on σ_FS.
- Each realization runs the critical-surface search for the chosen method, so the critical surface can migrate between samples. That's physical, not noise.

Details: `https://xslope.readthedocs.io/en/latest/reliability/`, with `reliability/taylor/` and `reliability/monte_carlo/` for the formulations.

Want me to put representative σ values on `shell` and `base` (say σ_c and σ_φ) and run a Taylor analysis with Spencer so you can see the breakdown on this model?
````

Tokens: input 26,319 · cached_input 25,044 · output 1,568 · calls 2 · 26.4s

## Turn 2 — How do I decide standard deviations for a reliability analy…

````text
You: How do I decide standard deviations for a reliability analysis if I only have a few tests?

Ran code:
    import math
    # d2: expected sample range in units of sigma (SQC table) -> range/d2 estimates sigma
    d2 = {3:1.693, 4:2.059, 5:2.326, 6:2.534, 8:2.847, 10:3.078, 20:3.735, 30:4.086}
    print("N   divide range by")
    for n,v in d2.items(): print(f"{n:<4}{v:.2f}")
    
    print("\nCOV-based sigmas for THIS model (imperial):")
    cov = {'gamma':(0.03,0.07), 'c':(0.20,0.40), 'phi':(0.05,0.13)}
    for m in slope_data['materials']:
        print(f"\n {m['name']}: gamma={m['gamma']}, c={m['c']}, phi={m['phi']}")
        for k,(lo,hi) in cov.items():
            v = m[k]
            print(f"   sigma_{k:<6} {v*lo:8.2f} .. {v*hi:8.2f}")
    
    # sampling error in s itself
    print("\nrelative uncertainty in s estimated from N tests  = 1/sqrt(2(N-1)):")
    for n in (3,5,10,20,30): print(f"  N={n:<3} +/- {100/math.sqrt(2*(n-1)):.0f}%")

Output:
    N   divide range by
    3   1.69
    4   2.06
    5   2.33
    6   2.53
    8   2.85
    10  3.08
    20  3.73
    30  4.09
    
    COV-based sigmas for THIS model (imperial):
    
     shell: gamma=130.0, c=300.0, phi=37.0
       sigma_gamma      3.90 ..     9.10
       sigma_c         60.00 ..   120.00
       sigma_phi        1.85 ..     4.81
    
     base: gamma=130.0, c=0.0, phi=37.0
       sigma_gamma      3.90 ..     9.10
       sigma_c          0.00 ..     0.00
       sigma_phi        1.85 ..     4.81
    
    relative uncertainty in s estimated from N tests  = 1/sqrt(2(N-1)):
      N=3   +/- 50%
      N=5   +/- 35%
      N=10  +/- 24%
      N=20  +/- 16%
      N=30  +/- 13%

Assistant: With three or five tests, the sample standard deviation you compute is itself nearly worthless — from N = 3 it carries about ±50% uncertainty (the last block above: relative error in s is 1/√(2(N−1))). So don't just trust `np.std` on four numbers. Use it as one of three estimates and reconcile them.

## 1. The three-sigma rule (Duncan's method)

Ask yourself the highest conceivable value (HCV) and lowest conceivable value (LCV) of the parameter — the range you'd bet the property will not escape. Then

```
σ ≈ (HCV − LCV) / 6
```

because ±3σ spans ~99.7% of a normal distribution. For your shell, if you'd say c is conceivably as low as 150 psf and as high as 450 psf, σ_c ≈ (450−150)/6 = **50 psf**.

The failure mode here is people picking too narrow a range — imagining plausible values rather than conceivable ones — which underestimates σ and inflates β. Duncan's own advice: widen your first guess.

## 2. The range method on the tests you actually have

With N tests, the observed range divided by the d₂ factor is an unbiased σ estimate, and it works far better on tiny samples than the sample s does:

| N | divide range by |
|---|---|
| 3 | 1.69 |
| 4 | 2.06 |
| 5 | 2.33 |
| 6 | 2.53 |
| 8 | 2.85 |
| 10 | 3.08 |

Four triaxials giving c = 260, 290, 320, 350 psf → range 90, σ ≈ 90/2.06 ≈ **44 psf**.

## 3. Published coefficients of variation as a floor

COV = σ/μ from the literature (Duncan 2000; Harr; Lacasse & Nadim), typical ranges:

| Property | COV |
|---|---|
| γ (unit weight) | 3–7% |
| c, Su | 20–40% |
| φ (drained) | 5–13% |
| tan φ | 5–15% |

Applied to your model's means:

| | shell | base |
|---|---|---|
| σ_γ | 3.9 – 9.1 pcf | 3.9 – 9.1 pcf |
| σ_c | 60 – 120 psf | 0 (c = 0 by definition) |
| σ_φ | 1.9 – 4.8° | 1.9 – 4.8° |

## How to reconcile them

Take the **largest** of the three unless you have a specific reason not to. A four-test range method giving σ_c = 44 psf against a 20% COV floor of 60 psf means your four samples happened to come from one homogeneous pocket — use 60. Small samples systematically under-sample the tails.

Two more things that matter for your model specifically:

- **`base` has c = 0 as a definition, not a measurement.** Leave `sigma_c = 0` there — randomizing a cohesionless material's cohesion invents strength that the soil model says doesn't exist. Put the base's uncertainty entirely in σ_φ.
- **Spatial averaging.** A long failure surface crosses high and low patches that partly cancel, so the σ relevant to FS is smaller than the point-property σ from lab tests — sometimes by half for a long surface in a spatially variable deposit. If you skip this correction (most people do), your β is conservative. Say so rather than silently claiming precision.

Sanity check afterward: run Taylor and look at the per-variable breakdown. If one property contributes 80% of σ_FS, that's the one worth more tests, and it's also the one where a sloppy σ guess is actually costing you an answer.

Want me to set σ_γ, σ_c, σ_φ on the two materials using the mid-range COVs above and run a Taylor analysis so you can see which variable dominates?
````

Tokens: input 30,160 · cached_input 25,044 · output 3,221 · calls 2 · 49.5s

## Session total

- Tokens: input 56,479 · cached_input 50,088 · output 4,789 · calls 4
- Wall time: 82.8s
