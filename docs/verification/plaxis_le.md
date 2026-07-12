# PLAXIS LE (SVSLOPE) Corpus

!!! note "Status: stub — corpus not yet started"
    This page will track the SVSLOPE verification manual (SoilVision Systems, 2019 —
    the product Bentley now sells as **PLAXIS LE**) the way the
    [Slide2 corpus](rocscience.md) tracks its manual. Nothing is built yet.

The SVSLOPE 5 verification manual (207 pages, ~105 problems) is one of the largest published
limit-equilibrium verification sets, organized as:

| Group | Scope |
|---|---|
| ACADS models | The Donald & Giam (1989) ACADS benchmark set — overlaps the problems already built in the [LEM samples](../lem/samples.md) and the [Slide2 corpus](rocscience.md) |
| SVSLOPE Group 1 | Core method verification — many problems shared with the Slide2 and SLOPE/W manuals (Arai & Tagyo, Greco, Fredlund & Krahn, Prandtl, etc.), giving a **third independent published value** for cross-bearing |
| SVSLOPE Group 2 | Extended features — loading, reinforcement, water conditions |
| SVSLOPE Group 3 | Further method/feature combinations |
| Dynamic programming (SAFE) models | Non-circular searches via the SAFE dynamic-programming method |
| 3D feature examples | Out of scope for XSLOPE (2D) |

The near-term value of this corpus is the third-source cross-bearing: for problems the
Slide2 and SLOPE/W corpora already cover, SVSLOPE's published values arbitrate whenever
the first two disagree (as they do on, e.g., OMS with high pore pressures). Full problem
enumeration and a status table will follow when the corpus round starts.

**Seepage companion:** SoilVision's SVFLUX verification manual plays the same role for
FE seepage and would join the [Slide2 groundwater corpus](rocscience_groundwater.md) as a
second seepage verification source.
