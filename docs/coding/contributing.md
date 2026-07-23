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

**FEM SSRM rows use a two-tier kernel scheme.** When the optional compiled Mohr-Coulomb kernel is built (`setup_kernel.py`), each `type=fem_ssrm` row is first solved with the fast kernel; if its factor of safety matches the expected value within tolerance the row passes immediately, annotated *via fast kernel*. On a miss the suite automatically re-solves the same row with the pure reference kernel — the oracle that defines every locked value — and that reference verdict is final (it either passes, annotated *via reference (fast missed by d=…)*, or fails). The run summary reports how many rows were decided by the fast kernel, how many needed the reference fallback, and how many failed, so a rising fallback count flags fast-kernel drift over time. When the compiled kernel is absent the rows are verified on the reference kernel only, exactly as before. Pass `--reference-only` to force the pure reference verdict for every row regardless of the fast kernel; use it for strict runs such as a pre-release check or immediately after a change to the constitutive physics. The `kernel_xcheck` gate, which compares the two kernels directly on small cases, is the companion guard that keeps this scheme sound and should not be removed while fast-first is the default.

## Documentation

If your change affects user-facing behavior, update the matching file under `docs/` and preview locally with `mkdocs serve`.

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
