# Rocscience Slide2 Groundwater Corpus

!!! note "Status: stub — corpus not yet started"
    This page tracks the [Slide2 Groundwater Verification Manual](https://www.rocscience.com/support/slide2/verification)
    (Rocscience, 2022; 21 problems) the same way the [Slide2 slope-stability corpus](rocscience.md)
    tracks its manual: every problem gets a row, built problems get an XSLOPE input file, a
    results section, and regression test tags. Nothing is built yet.

The manual verifies Slide2's finite-element groundwater engine against closed-form solutions
(Polubarinova-Kochina, Vedernikov, Terzaghi consolidation) and published numerical benchmarks.
It is the seepage analog of the slope-stability manual, and the natural verification target for
XSLOPE's own FE seepage solver — which the LEM corpus already exercises end-to-end on VP71,
VP72, VP76, VP77 and VP102, but only as a pore-pressure source, not against seepage-specific
quantities (flow rates, free-surface positions, pressure profiles).

Problems 1–13 are steady-state and analyzable with the current solver. Problems 15–21 are
**transient** (and 15–16 are consolidation), which XSLOPE's steady-state solver does not
support — those rows stay blocked until a transient seepage capability exists.

## Status

| # | Problem | Status | Notes |
|---|---|---|---|
| 1 | Shallow unconfined flow with rainfall | planned | Free surface between two rivers under uniform infiltration |
| 2 | Flow around cylinder | planned | Confined potential-flow anchor |
| 3 | Confined flow under dam foundation | planned | Classic flow-net problem |
| 4 | Steady unconfined flow through earth dam | planned | Free-surface dam seepage |
| 5 | Unsaturated flow behind an embankment | planned | |
| 6 | Steady-state seepage through saturated–unsaturated soils | planned | |
| 7 | Seepage within layered slope | planned | |
| 8 | Flow through ditch-drained soils | planned | |
| 9 | Seepage through dam | planned | |
| 10 | Steady unconfined flow, van Genuchten permeability | planned | Exercises the `vg` conductivity option |
| 11 | Earth/rock-fill dam, Gardner permeability function | planned | |
| 12 | Seepage from a trapezoidal ditch into a deep drainage layer | planned | Vedernikov analytical solution |
| 13 | Seepage from a triangular ditch into a deep drainage layer | planned | Vedernikov analytical solution |
| 14 | Unsaturated soil column | planned | |
| 15 | 1-D consolidation, uniform initial excess pore pressure | blocked | Transient/consolidation — no transient solver |
| 16 | Pore pressure dissipation of stratified soil | blocked | Transient/consolidation |
| 17 | Transient seepage, earth fill dam with toe drain | blocked | Transient |
| 18 | Transient seepage through an earth fill dam | blocked | Transient |
| 19 | Transient seepage through an earth fill dam (II) | blocked | Transient |
| 20 | Transient seepage in a layered slope | blocked | Transient |
| 21 | Transient seepage through a fully confined aquifer | blocked | Transient |

## Methodology

Same rules as the [Slide2 corpus](rocscience.md): problems are built from the manual's tabulated
data and coordinate-labeled figures; where a figure is unlabeled, geometry is extracted by
axis-calibrated pixel measurement and validated against printed solution quantities; every built
problem is locked into `run_tests.py` via test tags. Seepage problems will compare flow rates,
phreatic-surface positions, and head/pressure profiles rather than factors of safety, which will
need one or two new test-tag types (e.g. `seep_flowrate`, `seep_head_profile`) in the harness.
