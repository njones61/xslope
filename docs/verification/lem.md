# Limit Equilibrium Benchmarks

### ACADS simple homogeneous slope (circular search)

Source: [GeoStudio SLOPE/W Verification Manual (Oct 2022)](https://files.seequent.com/PDFs/Geostudio-Slope%20Stability%20Verification%20Manual-Oct2022.pdf),
ACADS suite (Donald & Giam, 1989; Giam & Donald, 1992). Full problem
description, input file, and figures:
[LEM sample problem 12](../lem/samples.md#verification-acads-simple). Single soil, c′ = 3.0 kPa, φ′ = 19.6°,
γ = 20.0 kN/m³; 2:1 slope, 10 m high; automated critical-circle search with 50
slices per method. The ACADS consensus answer is FOS ≈ 1.00.

| Method | XSLOPE FOS | Reference | Diff |
|---|---|---|---|
| Ordinary (OMS) | 0.942 | 1.00 | −5.8% |
| Bishop's Simplified | 0.985 | 1.00 | −1.5% |
| Simplified Janbu | 0.986 | 1.00 | −1.4% |
| Corps of Engineers | 0.990 | 1.00 | −1.0% |
| Lowe & Karafiath | 0.987 | 1.00 | −1.3% |
| Spencer | 0.984 | 1.00 | −1.6% |
| Morgenstern-Price | 0.984 | 1.00 | −1.6% |

All rigorous methods fall within the ACADS accepted band; OMS reads low, as
expected for a legacy method (its known conservative bias on effective-stress
problems is why it is reported for completeness only).

### ACADS weak-layer slope (non-circular search)

Source: [SLOPE/W Verification Manual](https://files.seequent.com/PDFs/Geostudio-Slope%20Stability%20Verification%20Manual-Oct2022.pdf),
ACADS weak-layer case (sec. 2.7). Full details:
[LEM sample problem 13](../lem/samples.md#verification-acads-weak-layer). A 2:1
slope with a thin c′ = 0, φ′ = 10° interlayer and piezometric line; the
critical surface runs along the weak layer. ACADS accepted band ≈ 1.26.

| Method | XSLOPE FOS | Reference | Diff |
|---|---|---|---|
| Spencer | 1.258 | ~1.26 | −0.2% |
| Morgenstern-Price | 1.248 | ~1.26 | −1.0% |
| Corps of Engineers | 1.336 | ~1.26 | +6.0% |
| Lowe & Karafiath | 1.249 | ~1.26 | −0.9% |
| Simplified Janbu | 1.278 | ~1.26 | +1.4% |

The non-circular slip is seeded just above the base of the weak layer (the
standard placement for a weak interlayer — see the sample-problem write-up).
This is exactly where complete equilibrium earns its keep: XSLOPE's Spencer and
Morgenstern-Price (half-sine) both land within ~1% of SLOPE/W's named
Morgenstern-Price value (1.261) — Spencer at 1.258 (−0.2%) and M-P at 1.248
(−1.0%).
Corps of Engineers reads modestly high here, consistent with ground-parallel
side-force inclinations on this geometry (XSLOPE uses the standard "Corps #2"
convention — see [Force Equilibrium Methods](../lem/force_eq.md)).

### Arai & Tagyo homogeneous slope (circular search)

Source: [Arai & Tagyo (1985)](https://doi.org/10.3208/sandf1972.25.43), Soils
and Foundations 25(1); republished in Greco (1996), Malkawi et al. (2001), Kim
et al. (2002). Full details:
[LEM sample problem 14](../lem/samples.md#verification-arai-tagyo). Homogeneous 1.5:1 slope,
20 m high, c = 41.65 kPa, φ = 15.0°, γ = 18.82 kN/m³. Published FOS ≈ 1.451.

| Method | XSLOPE FOS | Reference | Diff |
|---|---|---|---|
| Ordinary (OMS) | 1.344 | 1.451 | −7.4% |
| Bishop's Simplified | 1.404 | 1.451 | −3.2% |
| Simplified Janbu | 1.411 | 1.451 | −2.8% |
| Corps of Engineers | 1.476 | 1.451 | +1.7% |
| Lowe & Karafiath | 1.438 | 1.451 | −0.9% |
| Spencer | 1.401 | 1.451 | −3.4% |
| Morgenstern-Price | 1.400 | 1.451 | −3.5% |

---
