# Seepage Solver Refactor: Phreatic Surface Handling

## Problem

The stream function (flow net) solution near the phreatic surface is irregular, especially on coarse meshes. Flowlines bunch up at the phreatic exit point and do not intersect equipotential lines at right angles. The SEEP2D Fortran code handles this transition more smoothly.

## Root Cause

Both the head equation and the stream function equation evaluate the relative permeability `kr` at the **element centroid** and apply it uniformly to the entire element stiffness matrix:

```python
# Head equation
p_elem = (p_nodes[i] + p_nodes[j] + p_nodes[k]) / 3.0
kr_elem = kr_frontal(p_elem, kr0, h0)
ke = kr_elem * area * grad.T @ Kmat @ grad

# Stream function equation
ke = (1.0 / kr_elem) * area * grad.T @ Kmat_flow @ grad
```

For an element straddling the phreatic surface (some nodes saturated, some not), the centroid pressure determines whether the **entire element** is treated as saturated (`kr=1`) or unsaturated (`kr=kr0`). This creates an abrupt, element-sized jump at the phreatic surface rather than a smooth transition.

The effect is worse on coarse meshes where elements are larger relative to the phreatic transition zone.

## SEEP2D Approach

The SEEP2D Fortran code (`../xslope_private/ref_docs/ref_docs_seep/seep2d_fortran/src/seep2d.f`, subroutine `qdflow`) evaluates `kr` at each **Gauss point** individually during stiffness assembly:

```fortran
c     Compute pressure heads at the integration points (lines 1548-1564)
      do 260 nt = 1, 4
      do 250 ns = 1, 4
        pp = 0.d0
        do 240 i = 1, 4
          fni = (sss * si + 1.d0) * (ttt * ti + 1.d0) * 0.25d0
          pp = (hlast(nd1) - y(nd1)) * fni + pp
        continue
        pres(nt,ns) = pp          ! pressure at each Gauss point
      continue

c     Stiffness assembly with per-Gauss-point kr (lines 1604-1618)
      do 410 nt = 1, 4
      do 410 ns = 1, 4
        fsfact = fkrel(pres(nt, ns), mat)     ! kr at THIS Gauss point
        if (itime .eq. 2) then                 ! stream function solve
          fsfact = 1.0d0 / fsfact              ! invert for stream function
        end if
        ! ... each Gauss point contributes its own kr-weighted stiffness
```

Key differences from our code:

1. **Per-Gauss-point kr**: Each of the 4x4=16 Gauss points gets its own `kr` based on the interpolated pressure at that specific point. Within a single element straddling the phreatic, Gauss points below the phreatic get `kr=1` while points above get `kr=kr0`.

2. **Same assembly for head and stream function**: SEEP2D uses the same `qdflow` subroutine for both, with `itime=1` for head and `itime=2` for stream function. The only difference is inverting `fsfact` for the stream function. This ensures consistency.

3. **`1/kr` cap**: SEEP2D caps at `1e10` (not `1e12` like our code).

## Proposed Changes

### Phase 1: Gauss-Point kr Evaluation

Modify the stiffness assembly in both `solve_unsaturated` and `solve_flow_function_unsaturated` to evaluate `kr` at each Gauss point rather than the element centroid.

For **tri3** elements (constant gradient), this means:
- Evaluate pressure at 3 Gauss points using shape function interpolation
- Compute `kr` at each Gauss point
- The effective stiffness contribution at each Gauss point is weighted by its own `kr`
- Since the gradient is constant for tri3, this simplifies to: `ke = kr_avg * base_stiffness` where `kr_avg = mean(kr at Gauss points)`. Note: `mean(kr(p_i))` is NOT the same as `kr(mean(p_i))` because `kr` is nonlinear.

For **tri6** elements (linear gradient), this is more involved:
- Evaluate pressure at each Gauss point using quadratic shape function interpolation
- Compute `kr` at each Gauss point
- Each Gauss point contributes `w * kr * gradN^T @ K @ gradN * |detJ|` to the element stiffness
- This requires moving `kr` inside the Gauss quadrature loop

For the **stream function**, the same change applies but using `1/kr` instead of `kr`.

### Phase 2: Consistency

Ensure both the head and stream function solves use identical `kr` evaluation. Factor out a common stiffness assembly function that accepts a flag for head vs stream function mode (like SEEP2D's `itime` parameter).

### Phase 3: Validation

- Compare head solution and flow net against SEEP2D for the earth dam benchmark
- Verify that fine and coarse meshes produce similar flow nets (convergence test)
- Check that the orthogonality error (grad(h) perpendicular to grad(phi)) decreases with mesh refinement

## Files to Modify

- `xslope/seep.py`: `solve_unsaturated()`, `solve_flow_function_unsaturated()`, and their quad element equivalents
- Consider refactoring element stiffness computation into shared helper functions

## References

- `../xslope_private/ref_docs/ref_docs_seep/seep2d_fortran/src/seep2d.f` — subroutine `qdflow` (line 1458), function `fkrel` (line 2562)
