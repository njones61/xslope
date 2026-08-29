# w1_elastic_fem — recorded assistant session

Captured by `tools/assistant_sessions.py` against the live Studio assistant, offscreen.

- Project: `docs/lem/files/xslope_simple_mult_layers.xlsx`
- Provider / model: `anthropic` / `claude-opus-5`
- Confirm-before-running: off (unattended capture)

## Turn 1 — Suggest elastic properties for both soils and enter them.

```text
You: Suggest elastic properties for both soils and enter them.

Assistant: I'll classify both materials from their strengths and write the values in.

Ran code:
    s = suggest_elastic()
    print(s)

Output:
      embankment           Soft Clay    E = 167,100   nu = 0.45
      foundation           Medium Clay  E = 668,300   nu = 0.4
    [{'material': 'embankment', 'soil_type': 'Soft Clay', 'E': 167100.0, 'nu': 0.45, 'unit_system': 'imperial', 'reason': 'embankment classified as Soft Clay from c = 400, phi = 0; E is the midpoint of the published range for that soil type and nu its typical value, in psf. Last-resort fill — a stated value outranks it.'}, {'material': 'foundation', 'soil_type': 'Medium Clay', 'E': 668300.0, 'nu': 0.4, 'unit_system': 'imperial', 'reason': 'foundation classified as Medium Clay from c = 800, phi = 0; E is the midpoint of the published range for that soil type and nu its typical value, in psf. Last-resort fill — a stated value outranks it.'}]

Ran code:
    for rec in s:
        m = next(m for m in slope_data['materials'] if m['name'] == rec['material'])
        m['E'] = rec['E']; m['nu'] = rec['nu']
    for m in slope_data['materials']:
        print(m['name'], m['E'], m['nu'])

Output:
    embankment 167100.0 0.45
    foundation 668300.0 0.4
    === MODEL CHECKS: clean ===

Assistant: Entered. Both values come from `suggest_elastic()`, which classifies each material **from its strength** — these are last-resort estimates, not stated data. If the source gives E or ν, those override.

| Material | c (psf) | Classified as | E (psf) | ν |
|:---|---:|:---|---:|---:|
| embankment | 400 | Soft Clay | 167,100 | 0.45 |
| foundation | 800 | Medium Clay | 668,300 | 0.40 |

Both are in psf, the model's own stress unit (imperial). Model checks came back clean.

Two things worth your judgment:

- These are φ = 0 **total-stress** materials, i.e. undrained. For a strictly undrained FEM run a Poisson's ratio nearer 0.49 is the conventional choice (no volume change); 0.45 is the tool's typical-value fill for the soil class. Say the word and I'll set either or both to 0.49.
- E and ν affect **displacements** and the stress path, not the LEM factor of safety — nothing about your existing circular-search setup changes.

Want me to run the FEM/SSRM analysis with these? That takes several minutes, so I'll only start it if you ask.
```

Tokens: input 49,821 · cached_input 47,448 · output 892 · calls 3 · 15.1s

## Turn 2 — Run a strength reduction analysis and compare the result wi…

```text
You: Run a strength reduction analysis and compare the result with Spencer.

Assistant: Both runs, Spencer first then SSRM — the strength-reduction solve takes several minutes, so expect a wait.

Ran code:
    sp = run_lem(method='spencer', search=True)
    print('Spencer FS =', sp['FS'])
    print('  Xo,Yo,R =', sp['Xo'], sp['Yo'], sp['R'], ' Depth =', sp['Depth'])
    print('  x_entry/x_exit =', sp.get('x_entry'), sp.get('x_exit'))
    print('  warnings:', sp.get('warnings'))

Output:
    spencer (auto search, circular): FS = 1.244 on the circle Xo=18.50, Yo=43.75, R=43.75
      admissibility: interslice tension (min Z = -1436.3 vs max compression 7382.1)
      admissibility: line of thrust outside the slice on 15% of boundaries
    Ignoring fixed y limits to fulfill fixed data aspect with adjustable data limits.
    Spencer FS = 1.2441076756911937
      Xo,Yo,R = 18.5 43.75 43.75  Depth = 0.0
      x_entry/x_exit = 55.242346141747674 4.560317233248684
      warnings: ['interslice tension (min Z = -1436.3 vs max compression 7382.1)', 'line of thrust outside the slice on 15% of boundaries']

Assistant: Spencer is done (FS = 1.244). Now the strength-reduction run — this is the long one.

Ran code:
    fem = run_fem(analysis='ssrm')
    print('SSRM FS =', fem['FS'])

Output:
    Built a tri6 mesh (target size 1.2): 8849 nodes, 4304 elements.
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    Warning: Poisson's ratio 0.450 is close to incompressible limit (0.5)
    Consider using nu <= 0.4 for better numerical stability
    === SSRM Analysis (Non-Convergence Method) ===
      Bisection range: [1.00, 2.00], tolerance: 0.01
      Verifying lower bound F=1.00 converges...
        -> Converged in 427 iters (F_min=1.00)
      Verifying upper bound F=2.00 does not converge...
        -> Did NOT converge (191 iters, F_max=2.00)
    
      SSRM step 1: F = 1.5000  [1.0000, 2.0000]
        -> Did NOT converge (FAILED, max|u| = 8.06x elastic) (1131 iters)
    
      SSRM step 2: F = 1.2500  [1.0000, 1.5000]
        -> Converged (41615 iters)
    
      SSRM step 3: F = 1.3750  [1.2500, 1.5000]
        -> Did NOT converge (FAILED, max|u| = 8.02x elastic) (2771 iters)
    
      SSRM step 4: F = 1.3125  [1.2500, 1.3750]
        -> Did NOT converge (FAILED, max|u| = 8.01x elastic) (6551 iters)
    
      SSRM step 5: F = 1.2812  [1.2500, 1.3125]
        -> Did NOT converge (FAILED, max|u| = 8.00x elastic) (17081 iters)
    
      SSRM step 6: F = 1.2656  [1.2500, 1.2812]
    
    SSRM: trial F = 1.2656 is inconclusive at the iteration ceiling (50000 iterations, out-of-balance still falling) - raise max_iterations_ceiling to decide it. It is NOT counted as a failure: the bracket's upper edge carries this trial as an uncertainty rather than a measured failure, and the factor of safety is reported as the bracket midpoint, as on any other run.
    
      SSRM step 7: F = 1.2578  [1.2500, 1.2656]
    
    SSRM: trial F = 1.2578 is inconclusive at the iteration ceiling (50000 iterations, out-of-balance still falling) - raise max_iterations_ceiling to decide it. It is NOT counted as a failure: the bracket's upper edge carries this trial as an uncertainty rather than a measured failure, and the factor of safety is reported as the bracket midpoint, as on any other run.
    
      SSRM result: FS = 1.2539
      Final interval: [1.2500, 1.2578]
      Capturing at-failure mechanism: solve at F=1.442 (FS x 1.15, failed edge 1.258; cap off, early-exit off, up to 12000 iters)…
        at-failure field: converged=False iters=12000 max_disp=4.32
        failure snapshot: trial F=1.442 (margin 0.15), 12000 iterations, max disp 4.32
      SSRM completed in 511.8 seconds
    SSRM: FS = 1.254
      SSRM: trial F = 1.2578 is inconclusive at the iteration ceiling (50000 iterations, out-of-balance still falling) - raise max_iterations_ceiling to decide it. It is NOT counted as a failure: the bracket's upper edge carries this trial as an uncertainty rather than a measured failure, and the factor of safety is reported as the bracket midpoint, as on any other run.
    SSRM FS = 1.25390625

Assistant: Both runs are done.

## Results

| Method | FS | Basis |
|:---|---:|:---|
| Spencer (LEM, auto search) | **1.244** | circle Xo=18.50, Yo=43.75, R=43.75, Depth=0.0 |
| SSRM (FEM strength reduction) | **1.254** | tri6 mesh, 8,849 nodes / 4,304 elements |

Difference: 0.010, about 0.8% — close agreement between an independent limit-equilibrium solution and a continuum strength-reduction solution.

**Warnings you should know about, both of them real:**

1. **Spencer's admissibility screen flagged two items** (report-only; they do not change FS):
   - Interslice tension, min Z = −1436.3 against max compression 7382.1
   - Line of thrust outside the slice on 15% of slice boundaries

   The number converged, but that surface is not fully admissible in the Duncan & Wright sense.

2. **The SSRM bracket did not close cleanly.** Trials at F = 1.2656 and F = 1.2578 hit the 50,000-iteration ceiling with the out-of-balance force still falling — **inconclusive, not failures.** FS = 1.254 is the midpoint of the final bracket [1.2500, 1.2578], so the true SSRM value could sit anywhere in that band or slightly above it. Raising `max_iterations_ceiling` would decide those two trials; I have not done that.

3. **Poisson's ratio 0.45 drew a numerical-stability warning on every FEM assembly** — it's near the incompressible limit, and the engine suggests ν ≤ 0.4. That is the `suggest_elastic` value we entered last turn for the embankment. It did not stop the run, but it is one plausible contributor to the slow convergence near the bracket. I have not tested that — re-running SSRM at ν = 0.40 would measure it, and costs another ~8 minutes.

The SSRM run took 511.8 s. Say the word if you want the ceiling raised, ν changed, or a comparison against Bishop/Janbu (LEM only, seconds each).
```

Tokens: input 298,517 · cached_input 31,632 · output 1,457 · calls 3 · 544.5s

## Session total

- Tokens: input 348,338 · cached_input 79,080 · output 2,349 · calls 6
- Wall time: 568.0s
- Workbook written: `docs/tutorials/files/w1_elastic_fem_after.xlsx`
