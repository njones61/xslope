# w1_conceptual — recorded assistant session

Captured by `tools/assistant_sessions.py` against the live Studio assistant, offscreen.

- Project: `docs/lem/files/xslope_simple_mult_layers.xlsx`
- Provider / model: `anthropic` / `claude-opus-5`
- Confirm-before-running: off (unattended capture)

## Turn 1 — Why do all of the limit equilibrium methods give the same f…

```text
You: Why do all of the limit equilibrium methods give the same factor of safety on this model?

Ran code:
    
    rows=[]
    for m in ['oms','bishop','janbu','corps','lowe','spencer','mprice']:
        r = run_lem(method=m, search=False, plot=False)
        rows.append((m, r['FS'], r.get('warnings')))
    for m,fs,w in rows:
        print(f"{m:8s} FS={fs:.6f}  warn={w if w else ''}")
    print("phi:", [(mm['name'], mm['c'], mm['phi']) for mm in slope_data['materials']])
    print("surface:", {k:slope_data['circles'][0][k] for k in ('Xo','Yo','R','Depth')})

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
    phi: [('embankment', 400.0, 0.0), ('foundation', 800.0, 0.0)]
    surface: {'Xo': 20.0, 'Yo': 40.0, 'R': 40.0, 'Depth': 0.0}

Assistant: They don't all agree — the **moment-equilibrium** family does, and the force-equilibrium methods differ. Here is what your model actually gives on the defined circle (Xo=20, Yo=40, R=40):

| Method | FS | Type |
|:--|:--|:--|
| OMS | 1.247 | moment |
| Bishop | 1.247 | moment |
| Spencer | 1.247 | moment + force |
| Morgenstern-Price | 1.247 | moment + force |
| Lowe & Karafiath | 1.299 | force only |
| Janbu (corrected) | 1.321 | force only |
| Corps of Engineers | 1.363 | force only |

## Why the moment methods coincide

Both materials are φ = 0 (c = 400 psf embankment, c = 800 psf foundation). On a **circular** surface with φ = 0, take moments about the center:

> FS = (c · L_arc · R) / (Σ W·x) 

The resisting moment is `c · L · R` — it depends only on the arc length and cohesion, **not on the normal stress**. Interslice forces are internal to the free body and cancel in the moment equation about the center, so whatever assumption a method makes about them (none for OMS, horizontal for Bishop, constant inclination for Spencer, a half-sine function for Morgenstern-Price) never enters the answer. All four must return the identical number. This is expected behavior, not a bug — see `https://xslope.readthedocs.io/en/latest/lem/overview/`.

## Why the force methods don't

Janbu, Corps and Lowe satisfy **force** equilibrium instead, which requires the base normal N on each slice, and N *does* depend on the assumed interslice force inclination. Those assumptions differ (Janbu: horizontal + fo correction, here fo = 1.09; Corps: a single θ = 18.48° from the surface geometry; Lowe: average of base and surface slopes), so they spread over 1.25–1.36.

## Warnings you should know about

Every non-OMS/Bishop run flagged the Duncan & Wright admissibility screen:

- **Interslice tension** on all four (Spencer min Z = −1504 vs max compression 7773; Corps −1498; Lowe −994; M-P −1724).
- **Line of thrust outside the slice** on 15% of boundaries (Spencer) and 21% (Morgenstern-Price).

These don't change FS, but they say the converged solution is not fully statically admissible — typical of a steep crest-side entry with no tension crack modeled. Adding `tcrack_depth` usually clears it.

Also note: **this was a solve on the circle already in the file, not a search.** FS = 1.247 is that surface, not the critical one. Want me to run `run_lem(search=True)` to find the critical circle?
```

Tokens: input 33,813 · cached_input 31,632 · output 1,455 · calls 2 · 21.5s

## Turn 2 — Both soils are undrained. What should the pore-pressure opt…

```text
You: Both soils are undrained. What should the pore-pressure option be, and what would change in the analysis if I added a water table?

Ran code:
    
    import copy
    from shapely.geometry import LineString
    gs = list(slope_data['ground_surface'].coords)
    print("ground:", gs)
    base = run_lem(search=False, plot=False)['FS']
    
    saved = copy.deepcopy({k: slope_data[k] for k in ('materials','piezo_line')})
    # piezo line that never rises above ground -> pore pressure only, no ponded load
    pz = [(x, min(y, 8.0)) for (x, y) in gs]
    slope_data['piezo_line'] = pz
    for m in slope_data['materials']:
        m['u'] = 'piezo'
    resync_geometry()
    fs_piezo = run_lem(search=False, plot=False)['FS']
    
    # same, but give the soils a saturated unit weight (heavier below the table)
    for m in slope_data['materials']:
        m['gsat'] = m['gamma'] + 5.0
    resync_geometry()
    fs_gsat = run_lem(search=False, plot=False)['FS']
    
    for k, v in saved.items():
        slope_data[k] = v
    resync_geometry()
    print(f"\nu='none'  FS = {base:.4f}")
    print(f"u='piezo' FS = {fs_piezo:.4f}   (delta {fs_piezo-base:+.4f})")
    print(f"u='piezo' + gsat FS = {fs_gsat:.4f}   (delta {fs_gsat-base:+.4f})")
    print("restored u:", [m['u'] for m in slope_data['materials']], "piezo pts:", len(slope_data['piezo_line']))

Output:
    ground: [(-30.0, 0.0), (0.0, 0.0), (40.0, 20.0), (90.0, 20.0)]
    Spencer: FS=1.247, theta=7.91
      Spencer admissibility warning: interslice tension (min Z = -1504.0 vs max compression 7773.3)
      Spencer admissibility warning: line of thrust outside the slice on 15% of boundaries
    spencer: FS = 1.247
      admissibility: interslice tension (min Z = -1504.0 vs max compression 7773.3)
      admissibility: line of thrust outside the slice on 15% of boundaries
    Spencer: FS=1.247, theta=7.90
      Spencer admissibility warning: interslice tension (min Z = -1508.5 vs max compression 7772.2)
      Spencer admissibility warning: line of thrust outside the slice on 15% of boundaries
    spencer: FS = 1.247
      admissibility: interslice tension (min Z = -1508.5 vs max compression 7772.2)
      admissibility: line of thrust outside the slice on 15% of boundaries
    Spencer: FS=1.247, theta=7.90
      Spencer admissibility warning: interslice tension (min Z = -1508.5 vs max compression 7772.2)
      Spencer admissibility warning: line of thrust outside the slice on 15% of boundaries
    spencer: FS = 1.247
      admissibility: interslice tension (min Z = -1508.5 vs max compression 7772.2)
      admissibility: line of thrust outside the slice on 15% of boundaries
    
    u='none'  FS = 1.2471
    u='piezo' FS = 1.2472   (delta +0.0001)
    u='piezo' + gsat FS = 1.2472   (delta +0.0001)
    restored u: ['none', 'none'] piezo pts: 0

Assistant: **Keep `u = 'none'`.** That is the correct option for an undrained, total-stress (φ = 0, c = Su) analysis — and it is what the model already has.

## Why

In a total-stress analysis the shear strength is `s = Su`, full stop. The undrained strength already contains the effect of the pore water that will exist at failure; you are not computing `c' + (σ − u)·tan φ'`. So there is no `u` to subtract, and supplying one would be double-counting the water. `u = 'piezo'` or `u = 'seep'` belong to **effective-stress** runs, where the materials carry `c'` and `φ' > 0`.

## What a water table would actually do — measured

I temporarily added a piezometric line at elevation 8 (below ground everywhere, so no ponding), switched both materials to `u = 'piezo'`, and re-solved the file's circle with Spencer:

| Configuration | FS |
|:--|:--|
| `u = 'none'` (as-built) | 1.2471 |
| `u = 'piezo'`, table at el. 8 | 1.2472 |
| `u = 'piezo'` + `gsat = γ + 5` | 1.2472 |

Essentially nothing — the +0.0001 is slice-discretization noise. With φ = 0 the resisting term is `c · L · R` about the circle center, which contains no normal stress, so subtracting pore pressure from N cannot change FS. The `gsat` case also came back unchanged at four decimals; I did not chase that further, so treat it as untested rather than as proof that the weight split does nothing.

**I restored the model** — `u` is back to `'none'` on both materials and the piezo line is empty again. Nothing was left changed.

## The one case where water *would* matter

A **ponded / reservoir** water surface above the ground is a different thing entirely: it applies a hydrostatic **load** on the slope face, and that changes the driving moment even at φ = 0. In xslope that comes from the water definition with `water_loads = 'auto'` (do not also add a `dloads` row). So:

- water table **below** ground, undrained → no effect, leave `u = 'none'`
- water **above** ground (submerged toe, reservoir) → define it, and FS will change

If you plan to run a long-term/drained case as well, that is the analysis where you switch the materials to `c'`, `φ'` and `u = 'piezo'`, with a real phreatic surface. Want me to set up either the ponded case or a drained companion model?
```

Tokens: input 38,278 · cached_input 31,632 · output 2,925 · calls 2 · 41.9s

## Session total

- Tokens: input 72,091 · cached_input 63,264 · output 4,380 · calls 4
- Wall time: 70.4s
