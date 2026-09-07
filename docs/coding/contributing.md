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
```

If your change adds a new sample, add a `<!-- test: ... -->` tag to the sample so it becomes part of the suite automatically.

**FEM SSRM rows use a two-tier kernel scheme.** When the optional compiled Mohr-Coulomb kernel is built (`setup_kernel.py`), each `type=fem_ssrm` row is first solved with the fast kernel. It passes the row there only by reproducing the lock **exactly** — its factor of safety within half of the lock's last printed decimal, so that it prints as the locked number at every digit the lock records (0.0005 for a lock of 1.418) — annotated *via fast kernel*. Anything else, including a factor of safety that sits comfortably inside the row's pass tolerance, falls through: the suite re-solves the same row on the pure reference kernel — the oracle that defines every locked value — and that reference verdict is final, pass or fail. Such a row is annotated *kernel read X, verified via reference*, and when the reference goes on to pass, the annotation also carries the kernel-versus-reference gap, which is the drift the lock's tolerance absorbed.

The exact gate is what makes the fallback able to catch anything. A tag's `tolerance` is used twice — it is `solve_ssrm`'s bisection stopping width as well as the lock's comparison tolerance — so the closing bracket is never wider than the tolerance, and two kernels that bisect to adjacent intervals differ by at most one tolerance. A gate one tolerance wide can therefore never see a one-bracket-step divergence. One row in the corpus has exactly that: on `RS2-40-d20` (`vp077b.xlsx`) the reference returns 1.4180 and the fast kernel 1.4352, against a lock of 1.418 with a tolerance of 0.02. The row now falls through, the reference decides it, and the kernel's number is printed beside it.

The run summary reports how many rows the fast kernel decided, how many the reference re-solved, and how many failed. A row whose locked value sits farther from today's reference answer than its own printed precision falls through every run and pays both solves — `FEM-1-ssrm` is one, locked at 1.3633 where both kernels now return 1.3711, inside its 0.01 tolerance but nowhere near its fourth decimal. Those rows are recognizable by a printed kernel-versus-reference gap of zero: the two kernels agree and it is the lock that is behind. A *rise* in the fallback count, or a fallback row whose printed gap is not zero, is the drift signal. When the compiled kernel is absent the rows are verified on the reference kernel only. Pass `--reference-only` to force the pure reference verdict for every row regardless of the fast kernel; use it for strict runs such as a pre-release check or immediately after a change to the constitutive physics. The `kernel_xcheck` gate, which compares the two kernels directly on small cases, is the companion guard that keeps this scheme sound and should not be removed while fast-first is the default.

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
