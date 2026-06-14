# Plan: Comprehensive Pre-Release Audit of xslope Code and Documentation

## 1. Purpose and Scope

Before the JGGE manuscript is finalized and xslope is released publicly, conduct a
systematic audit of (a) the correctness of every computation that affects xslope's
results, and (b) the accuracy, concision, and coherence of the documentation.

This plan is informed by the June 2026 V&V benchmark campaign, which demonstrated
both that the methodology below finds real defects and that they exist: that work
surfaced a fragile Corps of Engineers side-force convention, a CW-winding bug that
silently zeroed tri3 seepage stiffness, a sign error in a hand-derived flow-vector
formula (caught by an independent finite-difference check), a viscoplastic
timestep that drove a limit cycle, a non-transferable SSRM failure criterion, and
documentation that contradicted the code's actual behavior. Every one of those was
found by one of four techniques, which this plan generalizes:

1. **External anchors** — compare against exact solutions or published benchmarks.
2. **Independent re-derivation** — re-derive the math from the cited source and
   check the code against the derivation, never against its own comments; where
   practical, verify numerically with a from-scratch implementation that shares
   no code (the "FD-check" pattern).
3. **Cross-method consistency** — different methods on the same input expose
   input-handling bugs and method-implementation bugs separately (e.g. Spencer
   validating a dam model that the FEM then mis-solved).
4. **Adversarial verification** — every reported finding is independently
   re-checked by a verifier that did not produce it, before any fix is written.

**Out of scope:** new features, performance work (except where a defect causes
wrong results), aesthetic code style.

---

## 2. Architecture: How Sub-Agents Are Used

### Roles

| Role | Tooling | Responsibility |
|---|---|---|
| **Coordinator** (main session) | full | Owns the plan, dispatches agents, triages findings, applies fixes, runs regression gates, commits |
| **Mapper** (Explore agents) | read-only | Stage 0 inventory: call graphs, dead code, doc-page-to-module map |
| **Auditor** (one per method/module) | read-only + scratch scripts | Re-derive theory, line-by-line review of one assigned unit, produce structured findings |
| **Verifier** (one per finding) | read-only + scratch scripts | Adversarially confirm or refute a finding it did not author; must attempt to *refute* |
| **Harness builder** | read-write (worktree) | Property-based and benchmark test scripts in `benchmarks/`/`run_tests.py` |
| **Docs auditor** (one per doc page) | read-only | Match-to-code check + prose review of one page |
| **Docs editor** | read-write | Applies approved doc edits page by page |

### Coordination rules (lessons from the pilot)

- **Finders never fix.** Auditors and verifiers are read-only; all source edits go
  through the coordinator (or a single editor agent per stage) so the diff stays
  reviewable and regression-gated. No parallel writers in `xslope/`.
- **Structured findings.** Every finding is one JSON/markdown block:
  `{id, file:line, severity (validity-critical / wrong-but-bounded / cosmetic),
  claim, evidence (the failing number or derivation step), proposed fix, verification recipe}`.
  A claim without a reproducible evidence artifact is returned to the auditor.
- **Adversarial gate.** A finding only reaches the fix queue when a separate
  verifier, prompted to refute it, fails to. (In the pilot, the FD check both
  caught a real sign error and prevented a wrong "fix" of a correct formula.)
- **Regression gate between stages.** `run_tests.py` (all modes) + the
  `benchmarks/` runners must pass before the next stage begins. Every fix adds or
  tightens a test annotation so it cannot regress silently.
- **One conversation thread per subsystem.** Use background agents for the wide
  fan-out (per-method audits, per-page doc reviews) and keep deep investigations
  (anything requiring iteration) in the main session, where judgment calls can be
  escalated to Norm. The pilot showed deep dives (e.g. the SSRM criterion) need
  human decisions; fan-out work (per-page doc review) does not.
- **Escalation triggers.** An agent must stop and report (not improvise) when: a
  fix would change a default or a published number; the cited literature is
  ambiguous; two verifiers disagree.

---

## 3. Stages

### Stage 0 — Inventory and Test Lockdown (foundation; ~1 session)

Goal: know exactly what exists, and freeze current behavior so every later change
is intentional.

- Mapper agents produce: module/function inventory with call graphs; list of dead
  or orphaned code; map of doc pages ↔ code units; inventory of all input-template
  sheets/options and which code paths consume them.
- Harness builder snapshots current outputs for every sample problem in
  `docs/*/files/` (FS per method, flowrates, SSRM FS) into a lockdown table —
  not as "correct" values but as *change detectors*.
- Exit: inventory docs in `plans/audit/`; lockdown table committed; full suite green.

### Stage 1 — LEM Engine Audit (the core; ~2-3 sessions)

The highest-stakes stage: each solution method audited independently against its
primary source.

**1a. Per-method theory audits (parallel, one auditor each):**

| Method | Primary source to re-derive from |
|---|---|
| OMS | Duncan, Wright & Brandon; Fellenius (1936) |
| Bishop's Simplified | Bishop (1955); DW&B |
| Simplified Janbu (+ fo correction) | Janbu (1968); correction-factor chart digitization |
| Spencer | Wright's UTEXAS formulation (the basis of `spencer()`) |
| Corps of Engineers | EM 1110-2-1902; force_equilibrium engine |
| Lowe & Karafiath | Lowe & Karafiath (1960); force_equilibrium engine |

Each auditor: (1) writes the method's equations from the source, (2) maps every
equation term to a line in `solve.py`, (3) checks sign conventions, effective vs
total stress handling, and every term of the slice force budget (w, u·dl, dloads
+ beta, kw seismic, tension crack t, reinforcement p, pile h/theta_p), (4) runs
the method on 3-4 canonical problems and compares to hand calculations on a
2-3-slice toy problem small enough to verify by hand. The shared
`force_equilibrium` engine gets its own auditor (it serves two methods).

> **Seed finding — Janbu base formulation does not match the standard Janbu
> Simplified (validity-critical; found June 2026 while preparing a SLOPE/W
> comparison for the paper).** `janbu()` in `solve.py` computes its base
> (pre-correction) factor of safety with the **Ordinary/Fellenius normal force**
> `N_eff = W·cosα − u·dl` and a moment-style `Σ(W·sinα)` denominator. The result
> is numerically **identical to OMS** for circular surfaces — verified on the
> benchmark critical surfaces: ACADS Simple OMS = Janbu-base = 0.9437; Arai &
> Tagyo OMS = Janbu-base = 1.3451 (exact match in both). The standard Janbu
> Simplified (and SLOPE/W's `F_f` at λ=0) instead uses the **iterative
> vertical-equilibrium normal** with `m_α = cosα + sinα·tanφ/F` and horizontal
> force equilibrium; a from-scratch implementation of that form gives 0.9376 and
> 1.3181 on the same surfaces. So xslope's "Simplified Janbu" is effectively
> **OMS × f₀**, not Janbu's method, which is why xslope's *uncorrected* base does
> not match SLOPE/W's *uncorrected* Janbu. Note SLOPE/W reports the uncorrected
> base (no f₀ applied; see the GeoStudio *Limit Equilibrium Formulation*), so the
> mismatch is purely in the base formulation, not the correction factor.
> Audit task: re-derive Janbu (1968/1973) base equation, decide whether to adopt
> the iterative `m_α` normal + horizontal-force denominator (changes published
> Janbu numbers — escalate to Norm), and either way reconcile the `janbu()` ≡
> `oms()` coincidence. Captured at `solve.py:janbu()` (the `N_eff` / numerator /
> denominator block) and `docs/lem/janbu.md`.

**1b. Cross-method consistency battery (one harness agent):**
- All methods on a frictionless (φ=0) circular problem → OMS, Bishop, Spencer
  must agree exactly (analytical property).
- c=0 dry slope → all FS must equal tanφ/tanβ for shallow surfaces.
- Symmetry: mirrored geometry must give identical FS (catches sign/direction bugs).
- Unit invariance: same problem in SI and English units → identical FS.
- Published method-ordering (Fredlund & Krahn 1977): OMS low, Bishop≈Spencer.

**1c. Independent re-implementation check (the FD-check pattern):**
For Bishop and Spencer, a from-scratch ~50-line implementation written by an
agent that has *not seen* `solve.py`, run on the ACADS benchmark; agreement to
4+ digits with `solve.py` on identical slice data, else investigate.

### Stage 2 — Geometry & Slice Pipeline Audit (~2 sessions)

Everything upstream of the solvers — where an error poisons every method equally
(the pilot's CW-winding bug was exactly this class).

- `generate_slices()`: slice boundary generation at material/piezo/dload/tcrack
  transition points; circular vs non-circular interpolation; alpha/dl computation;
  slice weight integration with multiple soil layers; right-facing vs left-facing
  slopes.
- Pore pressure: piezo-line u (vertical-distance rule, clamping, multiple piezo
  lines, piezo above ground); seepage-solution u interpolation onto slice bases.
- Distributed loads: resolution onto slice tops, normal vs vertical conventions,
  interpolation between dload points, the "ponded water as dload" path.
- Failure-surface machinery: circle generation from the three Options
  (Depth/Intercept/Radius), non-circular surface handling, Free vs Horiz vs Fixed
  endpoint semantics.
- Search: `circular_search` grid/refinement logic (does it terminate at a local
  minimum? does the starting-circle dependence matter?), `noncircular_search`
  move/shrink logic — characterize, at minimum, sensitivity to defaults.

  > **Seed item — solution admissibility filtering in the search (force methods).**
  > A free non-circular search can drive a force-equilibrium method onto a
  > *numerically degenerate* surface and report a spurious minimum. Demonstrated:
  > Corps **variant 1** (single crest-toe chord θ) on the ACADS weak layer searches
  > down to FS = 1.048 — 18% below Spencer (1.279) and below variant 1's value on
  > the genuine critical surface (1.360). The 1.048 surface is non-physical: 3/52
  > negative effective base normals, 37/52 tensile interslice forces (min −212),
  > steep base reversals (α down to −63.8°), and Spencer scores the SAME surface at
  > 1.52 — confirming it is a variant-1 artifact, not a mechanism. (This is the
  > reason variant 2 / per-slice ground-parallel is the search default; variant 1
  > and the abs() v1 are kept for fixed-surface "Corps #1" reproduction.)
  > Candidate fix: reject trial surfaces whose force-equilibrium solution is grossly
  > inadmissible during the search. **CRITICAL design constraint (Norm, from years
  > of UTEXAS use): do NOT reject any single occurrence of base tension or a
  > non-monotonic line of thrust — those are normal in valid solutions (small crest
  > slices routinely show base tension; the thrust line legitimately rises and
  > falls).** The filter must key on EXTENT/DEGREE of non-physicality (a large
  > fraction of slices in tension, multiple negative normals, gross thrust-line
  > excursion beyond the mass), not presence. Note this would HARDEN the search for
  > all force methods but would NOT make variant 1 a good default (the single-θ
  > assumption is weak on irregular surfaces; an admissibility-filtered search just
  > converges to the admissibility boundary). Also needs a "no admissible surface in
  > this neighborhood" policy so the search does not dead-end. Consider geometric
  > contortion limits (per-slice α-change caps) as a complementary, cheaper guard.
- `fileio.py`: every sheet parser vs the template; defaults when cells are blank;
  the validation rules (the pilot found seep-only models exempted but FEM-only
  models not — look for more inconsistencies).
- Same techniques: toy problems verifiable by hand, symmetry/unit invariance,
  property tests added to the harness.

### Stage 3 — Advanced Features Audit (~1-2 sessions)

- Rapid drawdown (3-stage logic vs Duncan-Wright-Brandon procedure).
- Reinforcement (tension distribution along lines, p resolution into slices).
- Piles (force + moment handling per `plans/plan_piles.md` vs implementation).
- Reliability (`advanced.py`): distribution handling, COV math.
- Seismic (kw) and tension crack (t, zw) terms end-to-end.
- Each: one auditor + verifier; cross-check LEM vs FEM where both support the
  feature (e.g. reinforced slope LEM vs `reinforce_fem`).

### Stage 4a — Seepage Engine Audit (~1-2 sessions)

The June 2026 campaign anchored the *outputs* (exact confined-radial and
sheetpile benchmarks, Kozeny free-surface bracketing) but did not line-audit the
solver internals. Full audit of `seep.py`, same techniques as Stage 1:

- **Assembly:** `solve_confined` and `solve_unsaturated` element stiffness for
  every element type (tri3/tri6/quad4/quad8/quad9); anisotropy handling (k1/k2
  rotation matrix — verify with an analytically-rotated exact case: an
  anisotropic problem solved in rotated coordinates must match the isotropic
  solution transformed back).
- **Unsaturated machinery:** the linear-frontal kr function and the per-GP kr
  refactor (tri3 centroid kr vs higher-order per-GP kr); the exit-face
  active-set iteration (cycling? mesh-dependence of the active set?);
  relaxation behavior; convergence-tolerance sensitivity (the closure-error
  lesson: document which reported quantities are tolerance-sensitive).
- **Postprocessing:** flow-potential solvers (confined + unsaturated paths) —
  including the documented open defect that tri6 phi_range disagrees with
  flowrate; `compute_velocity`; boundary-flux extraction
  (`create_flow_potential_bc_from_elements` — the inflow/outflow normalization
  logic); total flowrate integration.
- **Interfaces:** `import_seep2d` fidelity (BC type mapping, unit weight, the
  repeat-last-node triangle convention); `export_seep_solution`; the u-transfer
  path from seepage solutions onto LEM slice bases and FEM Gauss points
  (interpolation correctness at material boundaries).
- **Cross-checks:** tri3 vs tri6 vs quad8 on every seepage sample (must agree
  within discretization error — the winding bug showed how silently one element
  family can break); SEEP2D comparison when Norm's GMS reference run lands
  (SEEP-2); mesh-refinement convergence on flowrate for each sample.

### Stage 4b — FEM Spot-Audit (~1 session)

Largely exercised by the June 2026 campaign (flow rule FD-verified, dt regime
characterized, criterion validated on G&L Ex.1); remaining holes:
- `build_fem_data` material/element mapping, BC auto-assignment edges,
  dload-to-nodal-force conversion (tributary lengths, inward-normal convention
  on non-convex boundaries).
- 1D truss integration and the pile beam elements (newest, least-tested code);
  `compute_strains`.
- Known open item: the systematic ~6-11% FEM-above-LEM gap on wet dams
  (Ex.6 +10%, Johnson +11%) — root-cause or characterize for the manuscript.
- The ductile-reinforced-slope question: what should SSRM report when there is
  no displacement catastrophe (reinforce_fem grows smoothly to F=2.6)?

### Stage 5 — Documentation Audit (~2 sessions, highly parallel)

Two passes, different agent types:

**5a. Accuracy pass (one read-only auditor per doc page, fan-out):**
For each page in `docs/`: every formula vs the code; every function
signature/parameter/default vs current source; every claimed behavior actually
exhibited (run the snippet if cheap); every screenshot/figure still producible;
sample-file links resolve; test annotations match current baselines. Output:
per-page findings list. (Pilot examples of this class: the "non-convergence is
default" text, the stale "elastic-normalized" criterion description, sample
docs lagging API changes.)

**Documentation change control (applies to both passes):**

Docs are the public face and easy to damage at scale — every docs stage is
gated by a human-reviewable diff product:

- All docs edits land on a dedicated branch (`audit-docs`), one commit per page,
  never mixed with code commits.
- For each edited page the editor agent emits a **review digest** to
  `plans/audit/docs_review/<page>.md`: for every change, the old excerpt, the
  new excerpt, and a one-line rationale tagged with its finding id (or
  "editorial"). Digests are scannable in minutes, unlike raw diffs of
  long markdown lines.
- A combined `git diff --stat` + per-page word-diff file is generated for the
  whole branch as a backstop.
- Norm reviews and approves digests page-by-page before the branch merges;
  anything tone- or judgment-sensitive (engagement rewrites, section
  reordering) is explicitly held for his call rather than merged by default.
- Mechanical gate before merge: `mkdocs build` passes with no broken
  links/anchors; all `<!-- test: -->` annotations still parse and pass.

**5b. Editorial pass (after 5a fixes land):**
- One agent per chapter for: verbosity/redundancy (flag any paragraph that
  restates another, any concept explained twice in one page), narrative flow
  (does each section motivate the next?), stale artifacts ("changelog-style"
  passages, references to removed options, old TODO/plan remnants).
  **House rule (per Norm):** docs describe current behavior only — no
  "earlier versions did X" passages; history lives in commit messages.
- One cross-page agent for: terminology consistency (e.g. one name per concept
  across LEM/seep/FEM pages), consistent math notation, ordering of pages,
  duplicated content between overview pages and method pages.
- Editor agent applies approved edits page by page; Norm reviews tone-sensitive
  rewrites (engagement is a taste call).

### Stage 6 — Synthesis and Release Gate (~1 session)

- Consolidated findings report: every finding, severity, resolution, and the
  test that now guards it (goes in `plans/audit/findings.md` — the manuscript's
  V&V narrative can cite the process).
- Full regression: `run_tests.py` all modes + all `benchmarks/` runners + the
  sample-output lockdown diff (every intentional change accounted for).
- Fresh-eyes smoke test: one agent follows the public README/docs install-and-run
  path from scratch and reports friction.
- Release checklist: version bump, CHANGELOG (the audit history belongs here,
  not in docs), PyPI build per the established workflow.

---

## 4. Sequencing and Effort

```
Stage 0  ──►  Stage 1 (LEM)  ──►  Stage 2 (pipeline)  ──►  Stage 3 (advanced)
                                                                │
              Stage 4a (seepage) + 4b (FEM spot)  ◄─────────────┘
                       │
              Stage 5a (docs accuracy)  ──►  Stage 5b (docs editorial)
                       │
              Stage 6 (synthesis + release gate)
```

- Stages 1-2 are strictly ordered (no point auditing solvers against a slice
  pipeline that later changes under them — but note the converse risk too:
  pipeline fixes invalidate method-level lockdown numbers, hence the re-gate
  after each stage).
- Stage 5a can begin in parallel with Stage 3 for doc pages whose code is
  already audited (LEM pages after Stage 1).
- Rough effort: 10-14 working sessions of coordinator time; agent fan-out per
  stage is 4-8 parallel read-only auditors plus verifiers, which is cheap
  compared to coordinator attention — the bottleneck is triage and fix review,
  which stays serial deliberately.

## 5. Deliverables

1. `plans/audit/` — inventory, lockdown table, per-stage findings, final report.
2. Expanded `run_tests.py` + `benchmarks/` coverage: every audited property
   becomes a regression annotation (target: every solver term exercised by at
   least one locked test).
3. Fix commits, one logical defect per commit, each citing its finding id.
4. Audited documentation set.
5. A V&V section narrative for the manuscript that can honestly describe the
   audit process and its quantitative outcomes.
