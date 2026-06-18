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
```

If your change adds a new sample, add a `<!-- test: ... -->` tag to the sample so it becomes part of the suite automatically.

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
