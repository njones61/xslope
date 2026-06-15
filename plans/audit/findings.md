# xslope Pre-Release Audit — Consolidated Findings

This is the synthesis of the comprehensive pre-release audit (Stages 1–6) that
gates the JGGE Technical Note. Every confirmed defect, its resolution, and the
regression test that now guards it are listed below. Detailed per-stage notes
live in `plans/plan_comprehensive_audit.md`; this file is the citable summary.

**Bottom line:** all three solver engines (limit-equilibrium, finite-element
seepage, finite-element SSRM) were line-audited against their documented
methodology and verified correct, with a set of latent defects found and fixed.
**No fix changed a published sample/benchmark number** except where a number was
itself wrong (negative Corps FS, earth_dam_up OMS/Janbu marked n/a). Every fix is
regression-checked; the suite is 94+ public samples/benchmarks plus a private
CE 544 fixture repo (`xslope_private_tests`).

## 1. Confirmed defects found and fixed

| ID | Area | Severity | Defect | Resolution / guard |
|----|------|----------|--------|--------------------|
| TENSION-CRACK | LEM (6 methods) | Real (5–8% FS) | Water-filled tension-crack thrust wrong-signed on right-facing slopes for OMS/Bishop/Janbu/Corps/Lowe (Spencer correct); crack water *raised* FS instead of lowering it | Store T as a magnitude (like seismic kw); add `V=-V` to Spencer's facing swap. Guarded by `test/tension_crack_symmetry_check.py` (mirror symmetry + water-lowers-FS) |
| CP-FILEIO | fileio | Latent | `cp`/`s(c/p)`/`s(psi)` read by the wrong header (`cp` vs `c/p`) → cp strength, σ_cp, σ_psi always 0 (cp option unreadable from any template) | Read the actual headers. Guarded by the cp reliability fixture |
| CP-LEM | slice.py | Latent | cp undrained strength `Su=(r_elev−y)·cp` omitted the base cohesion `c` and the clamp | `Su = c + cp·max(0, r_elev−y)`. Guarded by the cp fixture FS |
| CP-FEM | fem.py | Latent | cp branch dead — read `strength_option`/`cp_ratio` (fileio stores `option`/`cp`) → cp material got constant/zero strength | Read correct keys + add base c. Verified by element-strength check |
| CP-RELIABILITY | advanced.py | Unconservative | TSPM perturbed only γ/c/φ → cp-material strength uncertainty (σ_cp) dropped from COV_F → β overestimated | Select perturbed params by material option (mc→{γ,c,φ}, cp→{γ,c,cp}). Guarded by cp reliability fixture |
| RELIABILITY-NEGSIGMA | advanced.py | Garbage output | When σ>mean the perturbed search returned the fs_fail sentinel (9999), used as a real FS → β garbage | Validate up front, abort with a clear message |
| RD-F2 | advanced.py | Robustness | τ_ff not clamped ≥0 (negative undrained cohesion possible) | `max(0, τ_ff)` |
| RD-F3 | advanced.py | Edge | Degenerate K1/Kf denominators `continue`d, silently keeping drained strength in the undrained Stage-2 solve | Route to the lower-of-two-curves fallback |
| RD-F8 | advanced.py | Latent | rapid_drawdown mutated the caller's slice_df in place | Operate on `df.copy()` |
| SEARCH-ADMISSIBILITY | search.py | Wrong minima | Search drawn to degenerate circles (near-zero driving, high base tension, non-positive FS) — produced negative shipped Corps FS | Three extent-based guards (driving<1%ΣW, >25% base tension, FS≤0). Fixed shipped negative Corps FS |
| SEEP-CONV | seep.py | Silent-wrong | Unsaturated solve could exhaust iterations and return a flowrate with no failure signal | Return `converged`/`closure_error`; seep test fails on non-convergence |
| FEM-WARNINGS | fem.py | Latent crash | `warnings.warn` used but `warnings` never imported (NameError on those paths) | Add the import; upgrade seep/FEM mesh-mismatch to a proper warning |
| SSRM-CEILING | fem.py | Silent FS=None | "FS exceeds F_max" / "FS<F_min" only printed at debug≥1 | Print unconditionally with an actionable message |

(Stages 1–2, earlier in the campaign: Janbu reformulated from OMS+f0 to a
Bishop-base horizontal-force balance; Corps/Lowe-Karafiath interslice sign made
facing-aware; F10 dload moment-arm sign in OMS/Bishop; battered-pile right-facing
moment sign in OMS/Bishop/Spencer. All facing fixes guarded by mirror-symmetry
tests (`test/pile_symmetry_check.py`, `test/tension_crack_symmetry_check.py`).)

## 2. Engines verified correct (no defect)

- **LEM force_equilibrium / Spencer**: complete-equilibrium and force-equilibrium
  formulations match the documented equations.
- **Seepage FE engine**: K-tensor anisotropy (exact rotation-invariance + directional),
  linear-front kr, per-Gauss-point kr (all element types), edge-by-edge seepage-face
  active set, 3-part convergence, inflow/outflow balance; tri3/tri6/quad4/8/9 agree
  within discretization error; `import_seep2d` validated on the SEEP2D reference files.
- **FEM viscoplastic SSRM** (Griffiths & Lane 1999): elastic-strain yield check
  (ε_el = ε − ε_vp), MOCOUF/MOCOUQ with Lode-corner freeze, Δt=4(1+ν)/(3E), two-part
  convergence, c/F + atan(tanφ/F) reduction; benchmark G&L Ex.1 = 1.377.

## 3. Method limitations documented (not bugs)

- **earth_dam_up OMS/Janbu**: on a fully-submerged face the simplified equations cannot
  balance the reservoir water load (only Janbu collapses; every other method is sensible
  on the same circle). Reported as *n/a* and documented in oms.md/janbu.md.
- **Wet-dam FEM/LEM gap**: investigated and dismissed — ±4% on a clean dam (Griffiths
  Ex.6 dry/full/seep), no systematic bias; larger historical gaps trace to a pathological
  high-k-seam earth dam (abnormal LEM) and the LEM submerged-slope limits above.
- **tri6 flow net** (prior concern): not broken — tri6 is the *best* element for flow nets;
  the φ_range/flowrate gap is a sheet-pile-tip singularity that hurts tri3 more.

## 4. Test infrastructure added

- `reliability` test type in `run_tests.py` (checks the lognormal reliability index β).
- `test/tension_crack_symmetry_check.py`, `test/pile_symmetry_check.py` — facing mirror-symmetry guards.
- `xslope_private_tests` repo (private GitHub, sibling-path discovery, skip-if-absent):
  CE 544 exam/homework fixtures — LEM (circular + non-circular), seismic, reinforcement,
  tension-crack details, seepage, FEM SSRM, reliability, and a derived cp fixture.

## 5. Documentation

The full docs set (LEM, seepage, FEM, usage) was swept against the current code: factual/
equation errors (e.g. ΣF_y, pile stiffness 12EI/L³, reliability β snippet, Janbu pile terms),
non-running code examples (imports, kwargs, the circular_search 4-tuple, undefined quickstart
variables), and ~50 typos were fixed; stale source line-number citations were stripped; the
empty `features.md` was removed. The deep methodology sections were verified accurate.
