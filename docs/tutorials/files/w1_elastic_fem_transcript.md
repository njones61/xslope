# w1_elastic_fem — recorded assistant session

Captured by `tools/assistant_sessions.py` against the live Studio assistant, offscreen.

- Project: `docs/lem/files/xslope_simple_mult_layers.xlsx`
- Provider / model: `anthropic` / `claude-opus-5-5`
- Confirm-before-running: off (unattended capture)

## Turn 1 — Suggest elastic properties for both soils and enter them.

```text
You: Suggest elastic properties for both soils and enter them.

Ran code:
    out = []
    for i, m in enumerate(slope_data['materials']):
        s = suggest_elastic(m['name'])
        m['E'] = s['E']; m['nu'] = s['nu']
        out.append((m['name'], s['soil_type'], s['E'], s['nu'], s['unit_system'], s['reason']))
    for o in out: print(o)
    print([(m['name'], m['E'], m['nu']) for m in slope_data['materials']])

Output:
      embankment: Soft Clay — E = 167,100, nu = 0.45 (imperial)
      foundation: Medium Clay — E = 668,300, nu = 0.4 (imperial)
    ('embankment', 'Soft Clay', 167100.0, 0.45, 'imperial', 'embankment classified as Soft Clay from c = 400, phi = 0; E is the midpoint of the published range for that soil type and nu its typical value, in psf. Last-resort fill — a stated value outranks it.')
    ('foundation', 'Medium Clay', 668300.0, 0.4, 'imperial', 'foundation classified as Medium Clay from c = 800, phi = 0; E is the midpoint of the published range for that soil type and nu its typical value, in psf. Last-resort fill — a stated value outranks it.')
    [('embankment', 167100.0, 0.45), ('foundation', 668300.0, 0.4)]
    === MODEL CHECKS: clean ===

Assistant: I added E and ν to both materials. The model didn't have any, so I used `suggest_elastic` as a fallback. It picks a soil type from each material's strength, then uses the midpoint of the published E range for that type and its typical ν.

| Material | Strength | Classified as | E (psf) | ν |
|:--|:--|:--|--:|--:|
| embankment | c = 400 psf, φ = 0 | Soft Clay | 167,100 | 0.45 |
| foundation | c = 800 psf, φ = 0 | Medium Clay | 668,300 | 0.40 |

The model checks came back clean.

These are generic textbook values, not site data. If you have measured or specified moduli, they should replace these. In an SSRM run, E and ν mostly affect the displacements it reports, not the factor of safety.

The model has no mesh yet. If you want to run a FEM/SSRM analysis, `run_fem` will build one from the model's settings. Expect it to take a few minutes.
```

Tokens: input 39,602 · cached_input 38,578 · output 628 · calls 2 · 7.3s

## Turn 2 — Run a strength reduction analysis and compare it with the l…

```text
You: Run a strength reduction analysis and compare it with the limit equilibrium result.

Assistant: Before starting: the SSRM run will take a few minutes. I'll run the LEM critical search first, then the strength reduction analysis.

Ran code:
    lem = run_lem(search=True)
    print('LEM method:', slope_data.get('lem_method'), 'FS', lem['FS'], 'Xo',lem['Xo'],'Yo',lem['Yo'],'R',lem['R'],'Depth',lem['Depth'],'entry',lem['x_entry'],'exit',lem['x_exit'])
    print('LEM warnings:', lem.get('warnings'))
    fem = run_fem(analysis='ssrm')
    print('SSRM FS', fem['FS'])
    print({k:v for k,v in fem.items() if k not in ('solution','fem_data')})

Output:
    spencer (auto search, circular): FS = 1.244 on the circle Xo=18.50, Yo=43.75, R=43.75
      admissibility: interslice tension on 15% of interior boundaries (min Z = -1436.3 vs max compression 7382.1)
      admissibility: line of thrust outside the slice on 15% of boundaries
      admissibility: Tension on the base of 2 slices near the crest. Consider adding a tension crack.
    Ignoring fixed y limits to fulfill fixed data aspect with adjustable data limits.
    LEM method: None FS 1.2441076756911937 Xo 18.5 Yo 43.75 R 43.75 Depth 0.0 entry 55.242346141747674 exit 4.560317233248684
    LEM warnings: ['interslice tension on 15% of interior boundaries (min Z = -1436.3 vs max compression 7382.1)', 'line of thrust outside the slice on 15% of boundaries', 'Tension on the base of 2 slices near the crest. Consider adding a tension crack.']
    Built a tri6 mesh (target size 1.2): 8849 nodes, 4304 elements.
    === SSRM Analysis (Non-Convergence Method) ===
      Bisection range: [1.00, 2.00], tolerance: 0.01
      Max iterations per trial: 12,000; convergence acceleration off
      Verifying lower bound F=1.00 converges...
        -> Converged in 303 iters (F_min=1.00)
      Verifying upper bound F=2.00 does not converge...
        -> Did NOT converge (361 iters, F_max=2.00)
    
      SSRM step 1: F = 1.5000  [1.0000, 2.0000]
        -> Did NOT converge (FAILED, max|u| = 15.02x elastic) (2271 iters)
    
      SSRM step 2: F = 1.2500  [1.0000, 1.5000]
        -> Converged (321 iters)
    
      SSRM step 3: F = 1.3750  [1.2500, 1.5000]
        -> Did NOT converge (FAILED, max|u| = 15.01x elastic) (5701 iters)
    
      SSRM step 4: F = 1.3125  [1.2500, 1.3750]
        -> FAILED: over the last 6,000 iterations the movement did not slow (each block of 1,200 iterations moved the slope 99% as far as the one before), so the trial was counted as sliding, max|u| = 13.04x elastic (12000 iters)
    
      SSRM step 5: F = 1.2812  [1.2500, 1.3125]
        -> FAILED: over the last 6,000 iterations the movement did not slow (each block of 1,200 iterations moved the slope 95% as far as the one before), so the trial was counted as sliding, max|u| = 6.48x elastic (12000 iters)
    
      SSRM step 6: F = 1.2656  [1.2500, 1.2812]
        -> Converged (326 iters)
    
      SSRM step 7: F = 1.2734  [1.2656, 1.2812]
        -> FAILED: over the last 6,000 iterations the movement did not slow (each block of 1,200 iterations moved the slope 93% as far as the one before), so the trial was counted as sliding, max|u| = 5.22x elastic (12000 iters)
    
      SSRM result: FS = 1.2695
      Final interval: [1.2656, 1.2734]
      Capturing at-failure mechanism: solve at F=1.460 (FS x 1.15, failed edge 1.273; cap off, early-exit off, up to 12000 iters)…
        at-failure field: converged=False iters=12000 max_disp=4.86
        failure snapshot: trial F=1.460 (margin 0.15), 12000 iterations, max disp 4.86
        capture growth: largest single step 1.0189 at iteration 1, longest run above 1.30x = 0 ending at iteration None (stops at 5 steps, or one step past 2.00x)
      SSRM completed in 399.4 seconds
    
    The factor of safety is 1.270, the midpoint of the bracket F = 1.2656 to 1.2734. At F = 1.2656 the slope reached equilibrium in 326 iterations. At F = 1.2734 it did not: over the last 6,000 iterations the movement did not slow (each block of 1,200 iterations moved the slope 93% as far as the one before), so the trial was counted as sliding. The run took 6 min 39 s.
    SSRM: FS = 1.270
    SSRM FS 1.26953125
    {'failure_solution': {'converged': False, 'stable': False, 'max_yield_violation': 0.4129282444137395, 'n_yield_above_1pct': 2262, 'max_tension_violation': None, 'max_yield_at': (1897, 0), 'yield_flagged': False, 'gate_failed': False, 'gate_deferrals': 0, 'gate_boundary': None, 'corrector_attempts': [], 'corrector': None, 'verdict': 'FAILED', 'u_ratio': 59.19215107384598, 'u_growth': 14.451319815463565, 'u_elastic_scale': 0.08207891255893245, 'exit_reason': 'iteration_cap', 'stop_reading': None, 'creep_reading': {'window': 6000, 'block': 1200, 'increments': [5.806422055328248, 5.804332215596151, 5.802329307547649, 5.800294924923166, 5.79843460104455], 'ratios': [0.9996400813250943, 0.9996549287714579, 0.9996493851836646, 0.9996792708124854], 'ratio': 0.9996559164181639, 'moved': 29.011813104439764, 'u_elastic_scale': 0.08207891255893245, 'slip_frac': None, 'slip_ratio': None, 'trend': 'steady', 'iteration': 12000}, '_k0_state': None, 'plateau_iteration': None, 'plateau_ratio': None, 'capture_truncated': False, 'capture_truncated_at': None, 'capture_truncated_reason': None, 'capture_truncated_kind': None, 'capture_truncated_max_u': None, 'capture_growth_max_step': 1.018869401281681, 'capture_growth_max_step_at': 1, 'capture_growth_longest_run': 0, 'capture_growth_longest_run_at': None, 'diverging_iteration': None, 'diverging_signal': None, 'early_exit_suppressed': False, 'budget_extensions': 0, 'iteration_budget': 12000, 'failure_criterion': 'hybrid', 'iterations': 12000, 'displacements': array([-1.77857774e-02, -9.76134323e-05, -2.57113647e+00, ...,
           -2.05608606e-02, -2.50536284e-04, -1.61398802e-03], shape=(17698,)), 'displacements_elastic': array([-0.00959641, -0.0034747 , -0.03479714, ..., -0.02068712,
           -0.00012822, -0.0016024 ], shape=(17698,)), 'stresses': array([[-1.58878980e+02,  2.40124139e+02, -2.71225461e+01,
             3.48749217e+02],
           [-2.55139870e+02,  2.10496019e+02, -3.15303778e+01,
             4.06939876e+02],
           [-2.06064285e+02,  1.85542051e+02, -2.08647637e+01,
             3.41062618e+02],
           ...,
           [ 1.72898262e+03,  2.68401257e+03, -1.99352987e+00,
             9.37453382e+02],
           [ 8.07674992e+02,  1.10170039e+03, -2.93170133e+01,
             3.22294595e+02],
           [ 8.45939796e+02,  1.18822394e+03, -3.09160167e+01,
             3.63475509e+02]], shape=(4304, 4)), 'strains': array([[ 1.69591255e-03, -1.76641254e-03, -4.70708460e-04,
             1.74708769e-03],
           [ 2.03963315e-03, -2.00089372e-03, -5.47205839e-04,
             2.03870618e-03],
           [ 1.70797400e-03, -1.69016596e-03, -3.62105413e-04,
             1.70868923e-03],
           ...,
           [ 7.58665836e-05, -1.92479469e-03, -8.35236216e-06,
             1.00033936e-03],
           [-9.20167217e-05, -7.07960993e-04, -1.22830521e-04,
             3.14036097e-04],
           [-6.76103908e-05, -7.84650338e-04, -1.29529922e-04,
             3.64322757e-04]], shape=(4304, 4)), 'vp_shear_strain': array([0., 0., 0., ..., 0., 0., 0.], shape=(4304,)), 'plastic_elements': array([False, False, False, ..., False, False, False], shape=(4304,)), 'yield_function': array([ -72.64313843,  -39.03662077,  -77.06822973, ...,  -70.44072958,
           -377.53741006, -359.29564435], shape=(4304,)), 'max_displacement': np.float64(4.86199636513332), 'plastic_strains': {0: array([[0., 0., 0., 0.],
           [0., 0., 0., 0.],
           [0., 0., 0., 0.]]), 1: array([[0., 0., 0., 0.],
           [0., 0., 0., 0.],
           [0., 0., 0., 0.]]), 2: array([[0., 0., 0., 0.],
           [0., 0., 0., 0.],
           [0., 0., 0., 0.]]), 3: array([[ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
             0.00000000e+00],
           [ 3.52903743e-02, -3.52903743e-02,  2.59896091e-02,
            -6.04552670e-20],
           [ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
             0.00000000e+00]]), 4: array([[ 6.39488822e-02, -6.39488822e-02,  5.89979084e-02,
            -8.51243550e-20],
           [ 9.27642845e-02, -9.27642845e-02,  7.07969495e-02,
             8.84508191e-20],
           [ 8.81894954e-02, -8.81894954e-02,  7.54350559e-02,
             9.35548001e-20]]), 5: array([[0., 0., 0., 0.],
           [0., 0., 0., 0.],
           [0., 0., 0., 0.]]), 6: array([[0., 0., 0., 0.],
           [0., 0., 0., 0.],
           [0., 0., 0., 0.]]), 7: array([[0., 0., 0., 0.],
           [0., 0., 0., 0.],
           [0., 0., 0., 0.]]), 8: array([[0., 0., 0., 0.],
           [0., 0., 0., 0.],
           [0., 0., 0., 0.]]), 9: array([[0., 0., 0., 0.],
           [0., 0., 0., 0.],
           [0., 0., 0., 0.]]), 10: array([[0., 0., 0., 0.],
           [0., 0., 0., 0.],
           [0., 0., 0., 0.]]), 11: array([[0., 0., 0., 0.],
           [0., 0., 0., 0.],
           [0., 0., 0., 0.]]), 12: array([[-5.38569857e-04,  5.38569857e-04, -1.35863582e-05,
             2.36425687e-21],
           [ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
             0.00000000e+00],
           [-5.69698528e-06,  5.69698528e-06, -3.49424368e-06,
             7.57317051e-24]]), 13: array([[ 6.78229180e-02, -6.78229180e-02,  5.82316537e-02,
             3.61229166e-20],
           [ 9.05534282e-02, -9.05534282e-02,  6.42524740e-02,
            -8.74265819e-20],
           [ 6.42125919e-02, -6.42125919e-02,  4.94231343e-02,
            -5.62924760e-20]]), 14: array([[0., 0., 0., 0.],
           [0., 0., 0., 0.],
           [0., 0., 0., 0.]]), 15: array([[0., 0., 0., 0.],
           [0., 0., 0., 0.],
           [0., 0., 0., 0.]]), 16: array([[0., 0., 0., 0.],
           [0., 0., 0., 0.],
           [0., 0., 0., 0.]]), 17: array([[0., 0., 0., 0.],
           [0., 0., 0., 0.],
           [0., 0., 0., 0.]]), 18: array([[0., 0., 0., 0.],
           [0., 0., 0., 0.],
           [0., 0., 0., 0.]]), 19: array([[0., 0., 0., 0.],
           [0., 0., 0., 0.],
           [0., 0., 0., 0.]]), 20: array([[0., 0., 0., 0.],
           [0., 0., 0., 0.],
           [0., 0., 0., 0.]]), 21: array([[ 1.34056497e-01, -1.34056497e-01, -1.01767569e-02,
            -1.03406525e-20],
           [ 1.53816304e-01, -1.53816304e-01, -1.08488028e-02,
             4.28088880e-21],
           [ 1.71440315e-01, -1.71440315e-01, -5.56930704e-03,
            -1.92490487e-19]]), 22: array([[0., 0., 0., 0.],
           [0., 0., 0., 0.],
           [0., 0., 0., 0.]]), 23: array([[ 2.07812611e-01, -2.07812611e-01, -8.98120768e-03,
             2.77386861e-20],
           [ 2.25444958e-01, -2.25444958e-01, -5.43829736e-03,
            -3.07848219e-20],
           [ 2.38298457e-01, -2.38298457e-01,  1.11534979e-03,
            -2.87516937e-19]]), 24: array([[0., 0., 0., 0.],
           [0., 0., 0., 0.],
           [0., 0., 0., 0.]]), 25: array([[0., 0., 0., 0.],
           [0., 0., 0., 0.],
           [0., 0., 0., 0.]]), 26: array([[ 2.44837384e-01, -2.44837384e-01,  3.35857577e-03,
             5.26814150e-20],
           [ 2.47927105e-01, -2.47927105e-01,  1.38330578e-02,
            -5.44693555e-20],
    [... 803,915 characters (15,831 lines) of printed output cut here. The snippet printed 819,896 characters, and only the first and last parts are shown. Print less: the values you need, not whole results or arrays.]
     'nr_predictor_iterations': 0}, {'F': 1.25, 'role': 'bisect', 'converged': True, 'stable': True, 'verdict': 'CONVERGED', 'u_ratio': 3.1072130066811385, 'growth': None, 'max_displacement': 0.25503666467735875, 'exit_reason': 'converged', 'stop_reading': None, 'creep_reading': None, 'plateau_iteration': None, 'budget_extensions': 0, 'diverging_iteration': None, 'ee_suppressed': False, 'iterations': 321, 'nr_force_evals': 76, 'nr_predictor_iterations': 0, 'corrector': {'driver_of_record': 'corrector', 'checkpoint': 'vp300', 'vp_iterations': 300, 'nr_iterations': 21, 'nr_force_evals': 76, 'oob': 9.593952335772276e-10, 'force_tol': 0.001, 'yield_violation': 1.2824876804984342e-08, 'yield_tol': 1e-06, 'max_disp_deep': 0.25503666467735875, 'disp_limit': 3.0, 'disp_frac_height': 0.008501222155911958}}, {'F': 1.375, 'role': 'bisect', 'converged': False, 'stable': False, 'verdict': 'FAILED', 'u_ratio': 15.007422282849333, 'growth': 3.3707075470484518, 'max_displacement': 1.2317929012889648, 'exit_reason': 'diverging', 'stop_reading': {'rule': 'diverging', 'signal': 'runaway', 'iteration': 5701, 'u_ratio': 15.007422282849333, 'window': 2000, 'gain': 4.753895474403051}, 'creep_reading': None, 'plateau_iteration': None, 'budget_extensions': 0, 'diverging_iteration': 5700, 'ee_suppressed': False, 'iterations': 5701, 'nr_force_evals': 0, 'nr_predictor_iterations': 0}, {'F': 1.3125, 'role': 'bisect', 'converged': False, 'stable': False, 'verdict': 'FAILED', 'u_ratio': 13.040736353046427, 'growth': 2.744196651922792, 'max_displacement': 1.0710427464040613, 'exit_reason': 'not_slowing', 'stop_reading': {'window': 6000, 'block': 1200, 'increments': [1.131824866759321, 1.1206600051448092, 1.111255336284139, 1.1025506049850342, 1.096425332235977], 'ratios': [0.9901355218970588, 0.9916079196031852, 0.9921667586063414, 0.9944444520538444], 'ratio': 0.9920874537192401, 'moved': 5.56271614540928, 'u_elastic_scale': 0.08207891255893245, 'slip_frac': None, 'slip_ratio': None, 'trend': 'steady', 'iteration': 12000, 'rule': 'not_slowing'}, 'creep_reading': {'window': 6000, 'block': 1200, 'increments': [1.131824866759321, 1.1206600051448092, 1.111255336284139, 1.1025506049850342, 1.096425332235977], 'ratios': [0.9901355218970588, 0.9916079196031852, 0.9921667586063414, 0.9944444520538444], 'ratio': 0.9920874537192401, 'moved': 5.56271614540928, 'u_elastic_scale': 0.08207891255893245, 'slip_frac': None, 'slip_ratio': None, 'trend': 'steady', 'iteration': 12000}, 'plateau_iteration': None, 'budget_extensions': 0, 'diverging_iteration': None, 'ee_suppressed': False, 'iterations': 12000, 'nr_force_evals': 0, 'nr_predictor_iterations': 0}, {'F': 1.28125, 'role': 'bisect', 'converged': False, 'stable': False, 'verdict': 'FAILED', 'u_ratio': 6.484764109269235, 'growth': 1.003001456315103, 'max_displacement': 0.5324974841747216, 'exit_reason': 'not_slowing', 'stop_reading': {'window': 6000, 'block': 1200, 'increments': [0.47674873625253156, 0.4502829362674528, 0.426728163179187, 0.4063374560308641, 0.38927655316246645], 'ratios': [0.9444869005985997, 0.9476889502331151, 0.9522161673220507, 0.9580129702168], 'ratio': 0.9505876636902753, 'moved': 2.149373844892502, 'u_elastic_scale': 0.08207891255893245, 'slip_frac': None, 'slip_ratio': None, 'trend': 'steady', 'iteration': 12000, 'rule': 'not_slowing'}, 'creep_reading': {'window': 6000, 'block': 1200, 'increments': [0.47674873625253156, 0.4502829362674528, 0.426728163179187, 0.4063374560308641, 0.38927655316246645], 'ratios': [0.9444869005985997, 0.9476889502331151, 0.9522161673220507, 0.9580129702168], 'ratio': 0.9505876636902753, 'moved': 2.149373844892502, 'u_elastic_scale': 0.08207891255893245, 'slip_frac': None, 'slip_ratio': None, 'trend': 'steady', 'iteration': 12000}, 'plateau_iteration': None, 'budget_extensions': 0, 'diverging_iteration': None, 'ee_suppressed': False, 'iterations': 12000, 'nr_force_evals': 0, 'nr_predictor_iterations': 0}, {'F': 1.265625, 'role': 'bisect', 'converged': True, 'stable': True, 'verdict': 'CONVERGED', 'u_ratio': 7.727916785760038, 'growth': None, 'max_displacement': 0.6342990061211045, 'exit_reason': 'converged', 'stop_reading': None, 'creep_reading': None, 'plateau_iteration': None, 'budget_extensions': 0, 'diverging_iteration': None, 'ee_suppressed': False, 'iterations': 326, 'nr_force_evals': 107, 'nr_predictor_iterations': 0, 'corrector': {'driver_of_record': 'corrector', 'checkpoint': 'vp300', 'vp_iterations': 300, 'nr_iterations': 26, 'nr_force_evals': 107, 'oob': 2.503641547134051e-06, 'force_tol': 0.001, 'yield_violation': 1.0339746179255371e-08, 'yield_tol': 1e-06, 'max_disp_deep': 0.6342990061211045, 'disp_limit': 3.0, 'disp_frac_height': 0.021143300204036816}}, {'F': 1.2734375, 'role': 'bisect', 'converged': False, 'stable': False, 'verdict': 'FAILED', 'u_ratio': 5.219298129362977, 'growth': 0.698542953081054, 'max_displacement': 0.42855424762938016, 'exit_reason': 'not_slowing', 'stop_reading': {'window': 6000, 'block': 1200, 'increments': [0.350889238127536, 0.3249632216169573, 0.304010793319243, 0.28409201900768966, 0.26732237585187507], 'ratios': [0.9261133893734425, 0.9355236934399566, 0.9344800423232454, 0.9409710867120111], 'ratio': 0.9342569041244739, 'moved': 1.531277647923301, 'u_elastic_scale': 0.08207891255893245, 'slip_frac': None, 'slip_ratio': None, 'trend': 'steady', 'iteration': 12000, 'rule': 'not_slowing'}, 'creep_reading': {'window': 6000, 'block': 1200, 'increments': [0.350889238127536, 0.3249632216169573, 0.304010793319243, 0.28409201900768966, 0.26732237585187507], 'ratios': [0.9261133893734425, 0.9355236934399566, 0.9344800423232454, 0.9409710867120111], 'ratio': 0.9342569041244739, 'moved': 1.531277647923301, 'u_elastic_scale': 0.08207891255893245, 'slip_frac': None, 'slip_ratio': None, 'trend': 'steady', 'iteration': 12000}, 'plateau_iteration': None, 'budget_extensions': 0, 'diverging_iteration': None, 'ee_suppressed': False, 'iterations': 12000, 'nr_force_evals': 0, 'nr_predictor_iterations': 0}]}}

Assistant: I ran the limit equilibrium search first and then the strength reduction analysis. The two agree within about 2%.

| Method | FS | Notes |
|:--|--:|:--|
| LEM, Spencer, automated circular search | **1.244** | Critical circle: Xo = 18.50, Yo = 43.75, R = 43.75 |
| FEM strength reduction (SSRM) | **1.270** | Bracketed between 1.2656 (stable) and 1.2734 (failed); tri6 mesh, 8,849 nodes; took 6 min 39 s |

**How they compare**
- The SSRM result is 0.026 higher than the LEM result, about 2.0%. That is a normal level of agreement for a φ = 0 slope.
- Even the highest factor that still stood in SSRM (1.2656) is above the LEM value. So the difference is real, not just bisection tolerance.
- The LEM critical circle bottoms out at elevation 0.0, the top of the foundation clay, without dipping into it. It runs from x = 55.2 on the crest side to x = 4.6 on the toe side.
- I have not compared the SSRM failure zone with that circle, because I didn't extract where the plastic strain is concentrated. The SSRM plot shows the mechanism if you want to check it by eye.

**Warnings from the Spencer run** (they don't change FS, but they matter for whether the solution is physically acceptable):
- Tension between slices on 15% of the inner slice boundaries.
- The line of thrust falls outside the slice on 15% of boundaries.
- Two slices near the crest have tension on their base. The solver suggests adding a tension crack. That would usually lower the LEM FS somewhat, but I didn't run that case.

**Other things to know**
- The model doesn't name an LEM method, so the run used Spencer by default.
- The SSRM uses the suggested elastic values entered earlier (generic Soft Clay and Medium Clay values). For φ = 0 soils they mainly affect the reported displacements, not the factor of safety.

I can rerun the LEM with a tension crack at the crest if you'd like to see how much the gap changes.
```

Tokens: input 50,196 · cached_input 39,130 · output 1,358 · calls 2 · 420.1s

## Session total

- Tokens: input 89,798 · cached_input 77,708 · output 1,986 · calls 4
- Wall time: 434.3s
- Workbook written: `docs/tutorials/files/w1_elastic_fem_after.xlsx`
