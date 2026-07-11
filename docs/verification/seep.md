# Finite-Element Seepage Benchmarks

### Confined radial flow (analytical anchor — exact)

Full details: [seepage sample problem 7](../seep/samples.md#verification-confined-radial).

A quarter-annulus confined flow problem in **plan view** (one quadrant of the
Thiem radial-flow-to-a-well geometry): inner arc (r = 10) at head 30, outer
arc (r = 30) at head 10, straight radial edges as no-flow streamlines. Confined
saturated flow obeys Laplace's equation in head with no gravity term, so the
model-plane orientation is irrelevant to the mathematics. Steady
flow has the exact solution q = k·(π/2)·Δh / ln(r₂/r₁) = 28.596
(k = 1) with a logarithmic head profile.

| Quantity | XSLOPE (tri6) | Exact | Diff |
|---|---|---|---|
| Discharge q | 28.5961 | 28.5960 | <0.01% |
| Max nodal head error | 0.004 | 0 | 0.02% of total drop |

Mesh-converged (identical at 2k and 6k nodes); quad8 gives the same result.
The only error source is faceting of the curved arcs by the polygon boundary.

### Partially penetrating sheetpile (analytical anchor — exact)

Full details: [seepage sample problem 8](../seep/samples.md#verification-sheetpile).

Pavlovsky's conformal-mapping solution for a cutoff wall of depth s in a
confined stratum of thickness T (Harr, 1962; Polubarinova-Kochina, 1962).
Boundary heads are 30 upstream / 20 downstream (downstream head at the stratum
top, so pressures are non-negative throughout the vertical section):
q = k·H·K(λ′)/(2·K(λ)) with λ = sin(πs/2T). At s/T = ½ the modulus is
self-dual and q = k·H/2 **exactly**. The closed form was additionally verified
with an independent finite-difference solution of the same boundary-value
problem (agreement ~0.4–0.5% at three penetration ratios).

| Case | XSLOPE q | Exact q | Diff | Head below wall tip |
|---|---|---|---|---|
| s/T = 0.50 | 5.010 | 5.000 | +0.20% | 25.0000 (exact: 25) |
| s/T = 0.75 | 3.412 | 3.403 | +0.27% | 25.0000 (exact: 25) |

The error halves with mesh refinement (set by the r^−½ singularity at the wall
tip) and converges to the exact value from above. The head on the wall plane
below the tip equals (h₁+h₂)/2 exactly, an antisymmetry property of the exact
solution that the FE solution reproduces to four decimals.

### SEEP2D cross-check — Johnson Reservoir (established code)

Full details: [seepage sample problem 4](../seep/samples.md#johnson-reservoir).

The Johnson Reservoir zoned earth dam (permeable shell, low-permeability core,
foundation; reservoir at el. 160 ft, tailwater at el. 100 ft) was exported to a
SEEP2D input file — the **exact same tri3 mesh topology, boundary conditions,
and material parameters** — and solved with the original USACE/WES SEEP2D
Fortran program. Identical-mesh comparison over all 2,604 nodes:

| Quantity | XSLOPE | SEEP2D | Diff |
|---|---|---|---|
| Total discharge q (ft³/day per ft) | 1.9575 | 1.9603 | −0.14% |
| Nodal heads | RMS Δh = 0.105 ft | (60-ft head range) | 0.18% |

The largest local head difference (~2 ft) occurs adjacent to the free surface,
where the two codes' unsaturated relative-permeability treatments differ in
detail; the bulk flow field agrees to about 0.1 ft.

---
