"""
Independent, from-scratch implementation of Bishop's Simplified Method.

Derived from published theory (Bishop, 1955):
  Moment equilibrium of the whole sliding mass about the circle center,
  with the simplifying assumption that interslice shear forces are zero
  (so the resultant interslice force is horizontal).

For each slice, vertical force equilibrium gives the effective normal force N'
on the base. The mobilized shear from Mohr-Coulomb (divided by FS) is:

    S_m = (1/FS) * [ c*b/cos(a) + (N') * tan(phi) ]   ... along base length dl=b/cos(a)

Solving vertical equilibrium for N' and substituting into moment equilibrium
about the center (where the lever arm for shear is R and the driving term is
W*sin(a)*R, R cancels) yields the classic iterative form:

    FS = SUM[ (c*b + (W - u*b)*tan(phi)) / m_a ] / SUM[ W*sin(a) ]

with
    m_a = cos(a) * (1 + tan(a)*tan(phi)/FS)

Here a = base inclination (alpha), b = horizontal slice width (dx),
W = slice weight, u = pore pressure at base center, c, phi = base strength.

This file deliberately re-derives the FS computation and does NOT consult
xslope's solve.bishop().
"""

import numpy as np


def bishop_independent(df, fs_init=1.5, tol=1e-10, max_iter=200):
    a = np.radians(df['alpha'].to_numpy(float))   # base inclination, radians
    phi = np.radians(df['phi'].to_numpy(float))   # friction angle, radians
    c = df['c'].to_numpy(float)                    # cohesion
    W = df['w'].to_numpy(float)                    # slice weight
    u = df['u'].to_numpy(float)                    # pore pressure at base center
    b = df['dx'].to_numpy(float)                   # horizontal slice width

    tan_phi = np.tan(phi)
    sin_a = np.sin(a)
    cos_a = np.cos(a)
    tan_a = np.tan(a)

    # Driving (denominator) term is independent of FS.
    denom = np.sum(W * sin_a)

    # Cohesive + frictional numerator pieces (independent of FS).
    cohesive = c * b
    frictional = (W - u * b) * tan_phi

    fs = fs_init
    for _ in range(max_iter):
        m_a = cos_a * (1.0 + tan_a * tan_phi / fs)
        numer = np.sum((cohesive + frictional) / m_a)
        fs_new = numer / denom
        if abs(fs_new - fs) < tol:
            fs = fs_new
            break
        fs = fs_new
    return fs


if __name__ == "__main__":
    import sys
    sys.path.insert(0, "/Users/njones/python_projects/xslope")
    from xslope.fileio import load_slope_data
    from xslope.search import circular_search
    from xslope import solve

    xlsx = "/Users/njones/python_projects/xslope/docs/lem/files/xslope_acads_simple.xlsx"
    slope_data = load_slope_data(xlsx)

    fs_cache, converged, _, _ = circular_search(
        slope_data, "bishop", num_slices=50
    )
    best = fs_cache[0]
    df = best["slices"] if "slices" in best else best.get("df_slices")
    # Fallback: keys may differ; print to discover.
    print("best keys:", list(best.keys()))
    print("converged:", converged)
    print("search best FS:", best.get("FS"))

    # Find the slice DataFrame in the cache entry.
    slice_df = None
    for k, v in best.items():
        try:
            import pandas as pd
            if isinstance(v, pd.DataFrame):
                slice_df = v
                print("Found slice df under key:", k)
                break
        except Exception:
            pass

    if slice_df is None:
        raise SystemExit("Could not locate slice DataFrame in cache entry")

    ok, res = solve.bishop(slice_df)
    fs_xslope = res["FS"]
    fs_mine = bishop_independent(slice_df)

    print(f"xslope bishop FS : {fs_xslope:.8f}")
    print(f"independent  FS  : {fs_mine:.8f}")
    print(f"abs diff         : {abs(fs_xslope - fs_mine):.2e}")
    print(f"rel diff         : {abs(fs_xslope - fs_mine)/fs_xslope:.2e}")
