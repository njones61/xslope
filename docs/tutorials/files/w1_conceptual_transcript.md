# w1_conceptual — recorded assistant session

Captured by `tools/assistant_sessions.py` against the live Studio assistant, offscreen.

- Project: `docs/lem/files/xslope_simple_mult_layers.xlsx`
- Provider / model: `anthropic` / `claude-opus-5`
- Confirm-before-running: off (unattended capture)

## Turn 1 — Why do all of the limit equilibrium methods give the same f…

````text
You: Why do all of the limit equilibrium methods give the same factor of safety on this model?

Ran code:
    rows=[]
    for m in ['oms','bishop','janbu','corps','lowe','spencer','mprice']:
        r = run_lem(method=m, search=False, plot=False)
        rows.append((m, r['FS'], r.get('warnings')))
    for m,fs,w in rows:
        print(f"{m:8s} FS={fs:.6f}  warn={w if w else ''}")
    print("phi values:", [(mat['name'], mat['c'], mat['phi']) for mat in slope_data['materials']])

Output:
    OMS: FS=1.247
    oms: FS = 1.247
    Bishop: FS=1.247
    bishop: FS = 1.247
    Janbu Corrected FS=1.321, fo=1.09
    janbu: FS = 1.321
    Corps Engineers: FS=1.363, theta=18.48
      Corps admissibility warning: interslice tension (min Z = -1497.8 vs max compression 8408.4)
    corps: FS = 1.363
      admissibility: interslice tension (min Z = -1497.8 vs max compression 8408.4)
    Lowe & Karafiath: FS=1.299
      Lowe & Karafiath admissibility warning: interslice tension (min Z = -994.4 vs max compression 7816.8)
    lowe: FS = 1.299
      admissibility: interslice tension (min Z = -994.4 vs max compression 7816.8)
    Spencer: FS=1.247, theta=7.91
      Spencer admissibility warning: interslice tension (min Z = -1504.0 vs max compression 7773.3)
      Spencer admissibility warning: line of thrust outside the slice on 15% of boundaries
    spencer: FS = 1.247
      admissibility: interslice tension (min Z = -1504.0 vs max compression 7773.3)
      admissibility: line of thrust outside the slice on 15% of boundaries
    Morgenstern-Price (half_sine): FS=1.247, lambda=0.151
      Morgenstern-Price admissibility warning: interslice tension (min Z = -1724.1 vs max compression 7698.9)
      Morgenstern-Price admissibility warning: line of thrust outside the slice on 21% of boundaries
    mprice: FS = 1.247
      admissibility: interslice tension (min Z = -1724.1 vs max compression 7698.9)
      admissibility: line of thrust outside the slice on 21% of boundaries
    oms      FS=1.247119  warn=
    bishop   FS=1.247119  warn=
    janbu    FS=1.321488  warn=
    corps    FS=1.363136  warn=['interslice tension (min Z = -1497.8 vs max compression 8408.4)']
    lowe     FS=1.298788  warn=['interslice tension (min Z = -994.4 vs max compression 7816.8)']
    spencer  FS=1.247119  warn=['interslice tension (min Z = -1504.0 vs max compression 7773.3)', 'line of thrust outside the slice on 15% of boundaries']
    mprice   FS=1.247119  warn=['interslice tension (min Z = -1724.1 vs max compression 7698.9)', 'line of thrust outside the slice on 21% of boundaries']
    phi values: [('embankment', 400.0, 0.0), ('foundation', 800.0, 0.0)]

Assistant: Not all of them do — the ones that satisfy **moment** equilibrium agree exactly, and the force-equilibrium methods don't. Solved on the surface the model already defines (Xo=20, Yo=40, R=40, 40 slices, no search):

| Method | FS | Equilibrium satisfied |
|:---|---:|:---|
| OMS | 1.247 | moment only |
| Bishop | 1.247 | moment (+ vertical) |
| Spencer | 1.247 | force + moment |
| Morgenstern-Price | 1.247 | force + moment |
| Lowe & Karafiath | 1.299 | force only |
| Janbu (corrected, fo=1.09) | 1.321 | force only (+ empirical correction) |
| Corps of Engineers | 1.363 | force only |

## Why the moment methods collapse onto one number

Both materials are **φ = 0** (c = 400 psf embankment, c = 800 psf foundation). On a circular arc with φ = 0 the moment equation degenerates:

```
FS = (sum of c * L_base * R) / (sum of W * x_arm)
```

The resisting moment depends only on cohesion and arc length — the **base normal force N drops out entirely**. Since the only thing the interslice force assumption ever changes is how N is distributed among the slices, every method that takes moments about the circle center gets the identical answer no matter what it assumes about the side forces. That is why OMS = Bishop exactly, and why Spencer and Morgenstern-Price converge to that same 1.247 (their λ/θ just satisfies force equilibrium without moving the moment answer: Spencer θ = 7.91°, MP λ = 0.151).

This is expected behavior, not a bug.

## Why the force-only methods differ

Janbu, Corps and Lowe never write the moment equation, so N is set by horizontal/vertical force balance on each slice, and there the interslice-force assumption *does* change N — and with φ = 0, N changes the driving side of the ratio. Their spread (1.299–1.363) is the assumption spread, and none of them is enforcing overall moment balance, so on a circular φ = 0 problem the 1.247 from Spencer/Bishop is the one to trust.

## Warnings you should know about

Every complete- and force-equilibrium run flagged **interslice tension** (e.g. Spencer min Z = −1504 vs max compression 7773), and Spencer and Morgenstern-Price also report the **line of thrust falling outside the slice** on 15% and 21% of boundaries. These are report-only admissibility screens (Duncan & Wright) — they don't change FS, but they say the internal force distribution is not fully admissible on this surface. With a φ = 0 clay and no tension crack modeled, that usually points at tension near the crest; adding `tcrack_depth` would be the normal fix.

Background: `lem/overview/` for the method comparison, `lem/bishop/` and `lem/spencer/` for the individual formulations.

Want me to run the automated search (`search=True`) to find the actual critical circle, or add a tension crack and re-check?
````

Tokens: input 33,819 · cached_input 15,863 · output 1,672 · calls 2 · 25.5s

## Turn 2 — Both soils are undrained. What should the pore-pressure opt…

```text
You: Both soils are undrained. What should the pore-pressure option be, and what would change in the analysis if I added a water table?

Ran code:
    import copy
    base_mats = copy.deepcopy(slope_data['materials'])
    base_piezo = copy.deepcopy(slope_data.get('piezo_line', []))
    gs = np.array(slope_data['ground_surface'].coords) if hasattr(slope_data.get('ground_surface'), 'coords') else None
    print("ground surface pts:", gs[:, :2].tolist() if gs is not None else None)
    print("gamma_sat:", [(m['name'], m.get('gamma'), m.get('gamma_sat')) for m in slope_data['materials']])
    
    try:
        # water table 3 ft below ground everywhere -> no ponded water anywhere
        wt = [(float(x), float(y) - 3.0) for x, y in gs[:, :2]]
        slope_data['piezo_line'] = wt
        for m in slope_data['materials']:
            m['u'] = 'none'
        r_none = run_lem(method='bishop', search=False, plot=False)
        for m in slope_data['materials']:
            m['u'] = 'piezo'
        r_piezo = run_lem(method='bishop', search=False, plot=False)
        print(f"\nFS bishop, u='none'  + piezo line present: {r_none['FS']:.4f}")
        print(f"FS bishop, u='piezo' + same piezo line   : {r_piezo['FS']:.4f}")
    finally:
        slope_data['materials'] = base_mats
        slope_data['piezo_line'] = base_piezo
        print("\nrestored ->", [(m['name'], m['u']) for m in slope_data['materials']], "piezo pts:", len(slope_data['piezo_line']))

Output:
    ground surface pts: [[-30.0, 0.0], [0.0, 0.0], [40.0, 20.0], [90.0, 20.0]]
    gamma_sat: [('embankment', 130.0, None), ('foundation', 135.0, None)]
    Bishop: FS=1.247
    bishop: FS = 1.247
    Bishop: FS=1.247
    bishop: FS = 1.247
    
    FS bishop, u='none'  + piezo line present: 1.2472
    FS bishop, u='piezo' + same piezo line   : 1.2472
    
    restored -> [('embankment', 'none'), ('foundation', 'none')] piezo pts: 0

Assistant: `u = 'none'` is the correct option, and it should stay that way. Here is why, and what a water table would and would not do.

## The pore-pressure option for undrained soils

Both materials are total-stress undrained: c = Su (400 and 800 psf) with **φ = 0**. In a total-stress analysis the strength is already the undrained strength at the in-situ effective stress — the pore pressure is *inside* Su, not applied on top of it. So:

- `u = 'none'` — correct for φ = 0 undrained (short-term / end-of-construction).
- `u = 'piezo'` or `'seep'` — for **effective**-stress runs, where you give c' and φ' and the solver subtracts u from the base normal.

Setting `u = 'piezo'` on a φ = 0 material is not just unnecessary, it is inert. I tested it: adding a piezometric line 3 ft below the ground surface and switching both materials to `u = 'piezo'`, Bishop gave

| Case | FS (Bishop, model's circle) |
|:---|---:|
| `u = 'none'`, piezo line present | 1.2472 |
| `u = 'piezo'`, same piezo line | 1.2472 |

Identical to four decimals. The Mohr-Coulomb base strength is `c + (N/L - u) tan φ`; with φ = 0 the `tan φ` term is zero, so u multiplies out of the equation entirely. **I restored the model afterwards** — both materials are back to `u = 'none'` and the piezo line is empty, unchanged from what you started with.

## What adding a water table *would* change

Three separate effects, only two of which apply here:

1. **Base pore pressure — no effect.** As measured above, φ = 0 kills it. This is the whole reason total-stress undrained analyses are run without a water definition.

2. **Slice weight — an effect, but only if you fill in `gamma_sat`.** Both materials currently have `gamma_sat = None`, so the solver uses γ = 130 / 135 everywhere and the weights are unchanged by a water table. If you set `gamma_sat`, the slice weight is split at the water table and the driving moment changes. Note the rule: give the **total** saturated unit weight, never a buoyant one — mixing a buoyant γ with an explicit water definition double-counts buoyancy.

3. **Ponded water above the ground — a real effect.** The toe here sits at (0, 0) with flat ground running out to x = −30 at elevation 0. Any water table above elevation 0 on that side stands *on* the ground, and standing water carries a hydrostatic load. `water_loads` ships on **auto**, so xslope derives that load from the water definition itself — you add no `dloads` row (one on top of auto counts the pool twice). That load stabilizes the toe and would raise FS. This applies even in a φ = 0 total-stress run: the water load is external, not a pore pressure.

So: below the ground surface, a water table does nothing to this model as configured; above it, it does, via the surface load.

## If you want the long-term case instead

That is a different analysis, not a different `u` setting on this one: give the materials effective-stress parameters c' and φ', set `u = 'piezo'` (or run `run_seep()` and use `u = 'seep'`), and define the water table. Then the pore pressure bites, the methods will no longer agree the way they did (the φ = 0 degeneracy in the moment equation disappears), and Spencer vs Bishop becomes a meaningful comparison.

Reference: `seep/seep_slope/` for how a seepage solution becomes slice pore pressure, `lem/overview/` for where u enters the slice equations.

Want me to add a reservoir at a level you name and re-run, or set up the effective-stress version alongside this one?
```

Tokens: input 38,030 · cached_input 31,726 · output 3,393 · calls 2 · 48.1s

## Session total

- Tokens: input 71,849 · cached_input 47,589 · output 5,065 · calls 4
- Wall time: 85.6s
