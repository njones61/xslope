# Contributing to xslope

Thanks for your interest in contributing. This page covers the basics for submitting changes.

## Getting Started

Fork the repo on GitHub, then clone your fork and add the upstream remote:

```bash
git clone https://github.com/YOUR_USERNAME/xslope.git
cd xslope
git remote add upstream https://github.com/njones61/xslope.git
```

Set up a dev environment:

```bash
python -m venv venv
source venv/bin/activate  # Windows: venv\Scripts\activate
pip install -e .[fem]  # [fem] adds gmsh, needed for mesh/seepage code
pip install -r docs/requirements.txt  # for building the docs
```

## Making Changes

Create a feature branch from an up-to-date `main`:

```bash
git fetch upstream
git checkout -b feature/your-feature-name upstream/main
```

Branch prefixes: `feature/`, `fix/`, `docs/`, `refactor/`.

## Code Style

- Follow PEP 8, 4-space indent, ~100 char lines.
- `snake_case` for functions/variables, `PascalCase` for classes, `UPPER_CASE` for constants.
- Type hints encouraged on new code.
- Docstrings on public functions — describe parameters, returns, and raises.
- Import order: stdlib, third-party, local.

## Testing

xslope has a regression test suite in `run_tests.py`. It auto-discovers test cases tagged in the docs sample pages (`docs/{lem,fem,seep}/samples.md` and `docs/seep/seep_slope.md`) and compares computed factors of safety and flowrates against expected values. Before submitting, run it and make sure it passes:

```bash
python run_tests.py                  # all tests
python run_tests.py --lem            # only LEM tests
python run_tests.py --fem            # only FEM tests
python run_tests.py --seep           # only seepage tests
python run_tests.py --skip-benchmarks  # skip slow verification benchmarks
python run_tests.py --reference-only    # strict: reference kernel only for FEM SSRM rows
python run_tests.py --gate           # the release gate: every row, every lock re-proved
```

If your change adds a new sample, add a `<!-- test: ... -->` tag to the sample so it becomes part of the suite automatically.

### How much of the suite to run

A strength-reduction lock is **proved** once, when it is cut: a bisection on two meshes at a budget that decided every trial, run by the corpus builder or the round that published it. What the suite owes it afterwards is a **check** that it is still reproducible. Separating the two is what keeps a suite that carries two hundred strength-reduction rows runnable, because the most expensive of them are hours apiece.

**`--standard` is the default**, and it is what a change is checked against. It runs every row except the ones the gate tier holds: a row whose tag says `tier=gate`, and a strength-reduction row this machine has already timed over ten minutes in the mode a standard run would use it. A row nothing has timed runs — unmeasured is unknown, not known-heavy — so a fresh clone checks everything and learns what things cost by running them. What it measured is written to `test/row_timings.json`, which is local and never committed; the run prints its slowest rows and names every row it held back, so nothing disappears silently.

**`--gate` is the release gate**, on the owner's word. It runs every row including the held-back ones, and forces the full bisection on every strength-reduction lock, so each one is re-proved rather than checked. Naming rows with `--benchmark` also overrides the hold — asking for a row by name is asking for it.

**`--quick` is the cheapest run.** It collapses each LEM problem's method list to one check, and it rides the standard tier, so it holds back the same rows.

### Checking a lock on its bracket edges

A bisection ends on a bracket: the highest trial factor at which the model stood, and the lowest at which it failed. The locked factor of safety is that bracket's midpoint, so the pair **is** the lock, and a lock that is still reproducible is one whose model still stands at the first factor and still fails at the second. A tag that carries both — `f_stand`, `f_fail`, `check=edges` — is checked that way: two solves instead of nine.

Both trials must be **decided**. A trial that ran out of iteration budget — `STABLE_STUCK`, `AMBIGUOUS`, `INCONCLUSIVE`, or simply at the ceiling — has not shown the model standing or failing anywhere, and it fails the check rather than passing it: a bracket edge nothing ruled on is exactly how a factor of safety ends up being a statement about the budget. The row reports the midpoint of its two re-solved edges as its computed value, which is the factor of safety those two trials bound.

The kernel scheme above still applies: the two trials run on the fast kernel first and on the reference kernel after, and only then, if an edge has genuinely flipped on the oracle, does the full bisection run — once, deciding the row against `expected_fs` exactly as a bracket-mode row is decided. The row then prints which edge moved and in which direction. The `--fem` summary counts all three: how many locks were checked on their edges, how many fell back to the bisection, and how many were held for the gate.

The two factors are never written by hand. `tools/lock_edges.py` reads the trial record the figure producers persist beside each model (`ssrm_run_record` → `*_fem_meta.json`) and appends the pair to the tag, refusing any record whose bracket was not decided at both edges, does not straddle the lock, or is wider than twice the tag's tolerance. A lock with no such record stays in bracket mode until the next run that cuts it writes one; `python tools/lock_edges.py --missing` lists those. The suite's `lock_edges` row re-checks the last three properties on every `check=edges` tag, so a pair that stops matching its lock fails the suite instead of quietly checking the wrong thing.

**FEM SSRM rows use a two-tier kernel scheme.** When the optional compiled Mohr-Coulomb kernel is built (`setup_kernel.py`), each `type=fem_ssrm` row is first solved with the fast kernel. It passes the row there only by reproducing the lock **exactly** — its factor of safety within half of the lock's last printed decimal, so that it prints as the locked number at every digit the lock records (0.0005 for a lock of 1.418) — annotated *via fast kernel*. Anything else, including a factor of safety that sits comfortably inside the row's pass tolerance, falls through: the suite re-solves the same row on the pure reference kernel — the oracle that defines every locked value — and that reference verdict is final, pass or fail. Such a row is annotated *kernel read X, verified via reference*, and when the reference goes on to pass, the annotation also carries the kernel-versus-reference gap, which is the drift the lock's tolerance absorbed.

The exact gate is what makes the fallback able to catch anything. A tag's `tolerance` is used twice — it is `solve_ssrm`'s bisection stopping width as well as the lock's comparison tolerance — so the closing bracket is never wider than the tolerance, and two kernels that bisect to adjacent intervals differ by at most one tolerance. A gate one tolerance wide can therefore never see a one-bracket-step divergence. A knife-edge model does exactly that: the two kernels' rounding lands one trial on either side of a verdict, and they bisect to adjacent intervals. The row falls through, the reference decides it, and the kernel's number is printed beside it. The corpus carries no such row today; `RS2-40-d20` (`vp077b.xlsx`) was one while its file carried a modulus contrast the vendor model does not.

The run summary reports how many rows the fast kernel decided, how many the reference re-solved, and how many failed. A row whose locked value sits farther from today's reference answer than its own printed precision falls through every run and pays both solves. Those rows are recognizable by a printed kernel-versus-reference gap of zero: the two kernels agree and it is the lock that is behind, which is a lock to re-record rather than a kernel to fix. The corpus carries none — every `fem_ssrm` lock reproduces today's reference answer at the precision it is quoted with. A *rise* in the fallback count, or a fallback row whose printed gap is not zero, is the drift signal. When the compiled kernel is absent the rows are verified on the reference kernel only. Pass `--reference-only` to force the pure reference verdict for every row regardless of the fast kernel; use it for strict runs such as a pre-release check or immediately after a change to the constitutive physics. The `kernel_xcheck` gate, which compares the two kernels directly on small cases, is the companion guard that keeps this scheme sound and should not be removed while fast-first is the default.

### Documentation checks

The numbers the docs print are checked by the same suite that checks the solver. `tools/verification_checks/` holds the checkers; its README describes each in full.

**Verification pages.** The six pages under `docs/verification` run as one suite row. Every printed percentage and absolute difference is re-derived from two numbers the page prints in the same sentence or row; every value a `<!-- test: -->` tag locks must be printed in the section carrying it, every number the section attributes to XSLOPE must agree with the lock it restates, and every value presented as locked must have a tag behind it; captions are checked against their figures. The row is change-gated on a committed content hash, so an unchanged page costs one file read. A page you edit stays a failure until you run `python -m tools.verification_checks.certify --recertify <page>` and commit `certified.json` with the edit.

**Tutorials.** `python run_tests.py --tutorials` sweeps the tutorial pages for factors of safety they attribute to a run of their own, and scores each against the locks in scope — the page's own tags, plus every tag anywhere under `docs/` on a model file the page links, which is how LEM-3 inherits the seven method locks its workbook carries on `docs/lem/samples.md`. A number is *guarded* when it restates a lock verbatim or correctly rounded, *disagreeing* when its column header, row label or sentence names a method whose lock says otherwise, and *unguarded* when nothing in the docs locks it. The row reports the per-page tally and passes; it fails only if the checker raises. A tutorial legitimately re-runs a sample under settings the sample's tag does not use, so read a finding before treating it as a defect.

**Before a lock moves.** A re-lock round changes a tag and re-measures the verification section that carries it; the tutorials print the same answer in prose, in a results table and in a console transcript, and go stale silently. Run `tools/tutorial_quotes.py` first and fix them in the same round:

```bash
python tools/tutorial_quotes.py --benchmark FEM-1-ssrm   # that benchmark's locks
python tools/tutorial_quotes.py --value 1.3633           # a value on its way out
python tools/tutorial_quotes.py --since <ref>            # every lock a diff moved
```

It matches a value at the precision the tutorials write it — the lock's digits and the lock rounded down to `--dp` places, three by default — and marks a hit on a page that links the workbook the lock is recorded on.

Mutation fixtures for all of this live in `tools/verification_checks/mutations.py`: each plants one defect and requires the checks to catch it. Run `python -m tools.verification_checks.mutations` after any change to the check logic.

## Documentation

If your change affects user-facing behavior, update the matching file under `docs/` and preview locally with `mkdocs serve`.

**Sample files are linked as workbooks and published as packages.** Write a link to a sample the plain way — `[name.xlsx](files/name.xlsx)` — and the build packages that project into `name.xslz` and turns the link into the pair *Download · Open in Studio*. Nothing is hand-maintained and no package is committed; a sample that gains a sidecar gains it in its package on the next build. Where a page deliberately offers the bare workbook instead, mark the link `{: .raw-file }` and it is left exactly as written.

## Submitting a Pull Request

Push your branch and open a PR against `njones61/xslope:main`. Include:

- A short description of what changed and why
- How you tested it
- Any related issue numbers

Keep PRs focused on a single change. Respond to review comments by pushing follow-up commits to the same branch.

## Reporting Bugs / Requesting Features

Open a GitHub issue. For bugs, include a minimal reproduction, the error/traceback, and your Python version and OS.

## License

Contributions are licensed under Apache License 2.0, matching the project.
