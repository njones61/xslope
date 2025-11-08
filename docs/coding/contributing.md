# Contributing to xslope

Thank you for your interest in contributing to xslope! This guide will help you get started with contributing to the project.

## Table of Contents

- [Code of Conduct](#code-of-conduct)
- [Getting Started](#getting-started)
- [Development Workflow](#development-workflow)
- [Code Standards](#code-standards)
- [Testing](#testing)
- [Documentation](#documentation)
- [Submitting Changes](#submitting-changes)
- [Review Process](#review-process)

## Code of Conduct

By participating in this project, you agree to maintain a respectful and collaborative environment. We expect all contributors to:

- Be respectful and considerate in communication
- Welcome newcomers and help them get started
- Focus on what is best for the community
- Show empathy towards other community members

## Getting Started

### Prerequisites

- Python 3.9 or higher
- Git
- Familiarity with limit equilibrium methods and/or finite element analysis (helpful but not required)

### Fork and Clone the Repository

1. **Fork the repository** on GitHub:
   - Go to https://github.com/njones61/xslope
   - Click the "Fork" button in the top right corner
   - This creates your own copy of the repository

2. **Clone your fork** locally:

```bash
git clone https://github.com/YOUR_USERNAME/xslope.git
cd xslope
```

3. **Add the upstream repository** as a remote:

```bash
git remote add upstream https://github.com/njones61/xslope.git
```

4. **Verify your remotes**:

```bash
git remote -v
# Should show:
# origin    https://github.com/YOUR_USERNAME/xslope.git (fetch)
# origin    https://github.com/YOUR_USERNAME/xslope.git (push)
# upstream  https://github.com/njones61/xslope.git (fetch)
# upstream  https://github.com/njones61/xslope.git (push)
```

### Set Up Development Environment

1. **Create a virtual environment**:

```bash
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate
```

2. **Install the package in editable mode** with dependencies:

```bash
pip install -e .
```

3. **Install development dependencies** (if you plan to build documentation):

```bash
pip install mkdocs mkdocstrings mkdocstrings-python
```

4. **Verify installation**:

```bash
python -c "import xslope; print(xslope.__version__)"
```

## Development Workflow

### Keep Your Fork Updated

Before starting new work, sync your fork with the upstream repository:

```bash
# Switch to main branch
git checkout main

# Fetch upstream changes
git fetch upstream

# Merge upstream changes
git merge upstream/main

# Push to your fork
git push origin main
```

### Create a Feature Branch

Always create a new branch for your changes:

```bash
# Create and switch to a new branch
git checkout -b feature/your-feature-name

# Or for bug fixes:
git checkout -b fix/bug-description
```

Branch naming conventions:
- `feature/description` - for new features
- `fix/description` - for bug fixes
- `docs/description` - for documentation updates
- `refactor/description` - for code refactoring

### Make Your Changes

1. **Write your code** following the project's code standards (see below)
2. **Test your changes** thoroughly
3. **Update documentation** if needed
4. **Add comments** where the code might be unclear

### Commit Your Changes

Write clear, concise commit messages:

```bash
git add <modified_files>
git commit -m "Brief description of changes

More detailed explanation if needed. Describe what changed
and why you made these changes.

Fixes #123"  # Reference issue numbers if applicable
```

Commit message guidelines:
- Use present tense ("Add feature" not "Added feature")
- First line should be 50 characters or less
- Add more detail in subsequent lines if needed
- Reference issues and pull requests when relevant

### Push Your Changes

```bash
git push origin feature/your-feature-name
```

## Code Standards

### Code Style

- Follow [PEP 8](https://pep8.org/) style guidelines for Python code
- Use 4 spaces for indentation (no tabs)
- Maximum line length: 100 characters (flexible for readability)
- Use descriptive variable names

### Naming Conventions

- **Functions and variables**: `snake_case`
- **Classes**: `PascalCase`
- **Constants**: `UPPER_CASE_WITH_UNDERSCORES`
- **Private methods/attributes**: prefix with single underscore `_private_method`

### Type Hints

While not required, type hints are encouraged for new code:

```python
def generate_slices(
    failure_surface: LineString,
    ground_surface: LineString,
    n_slices: int
) -> pd.DataFrame:
    """Generate analysis slices."""
    # implementation
```

### Docstrings

All public functions, classes, and methods should have docstrings:

```python
def load_slope_data(filepath: str) -> dict:
    """
    Load slope geometry and material properties from Excel file.

    Parameters:
        filepath (str): Path to the Excel input template file.

    Returns:
        dict: Dictionary containing all slope data including profile lines,
              materials, failure surfaces, and analysis parameters.

    Raises:
        ValueError: If required input data is missing or invalid.
        FileNotFoundError: If the specified file does not exist.

    Example:
        >>> data = load_slope_data('inputs/slope/example.xlsx')
        >>> print(data['gamma_water'])
        9.81
    """
    # implementation
```

### Imports

Organize imports in this order:
1. Standard library imports
2. Third-party imports
3. Local application imports

```python
# Standard library
import os
import pickle

# Third-party
import numpy as np
import pandas as pd
from shapely.geometry import LineString, Point

# Local
from .mesh import import_mesh_from_json
from .global_config import EPSILON
```

## Testing

### Running Tests

```bash
# Run specific test scripts
python test_griff_ex1.py

# Run main analysis
python main.py

# Test file I/O
python main_fileio_test.py

# Test mesh operations
python main_mesh.py
```

### Manual Testing

When making changes, verify:
1. **Basic functionality** - Can the code run without errors?
2. **Edge cases** - Does it handle unusual inputs gracefully?
3. **Backwards compatibility** - Do existing examples still work?
4. **Performance** - Does your change significantly impact speed?

### Adding Tests

If adding new features, consider adding test scripts in the root directory following the pattern of existing `test_*.py` and `main_*.py` files.

## Documentation

### Building Documentation Locally

```bash
# Install mkdocs if needed
pip install mkdocs mkdocstrings mkdocstrings-python

# Serve documentation locally
mkdocs serve

# View at http://127.0.0.1:8000/
```

### Updating Documentation

When making code changes that affect user-facing functionality:

1. **Update relevant markdown files** in the `docs/` directory
2. **Update API documentation** if function signatures change
3. **Add examples** for new features
4. **Update CLAUDE.md** if changing core architecture

### Documentation Structure

- `docs/usage/` - General usage guides
- `docs/lim_eq/` - Limit equilibrium method documentation
- `docs/seepage/` - Seepage analysis documentation
- `docs/fem/` - Finite element method documentation
- `docs/api/` - API reference documentation

## Submitting Changes

### Create a Pull Request

1. **Push your branch** to your fork (if not already done):

```bash
git push origin feature/your-feature-name
```

2. **Go to the original repository** on GitHub: https://github.com/njones61/xslope

3. **Click "Pull Requests"** then **"New Pull Request"**

4. **Click "compare across forks"**

5. **Select your fork and branch**:
   - Base repository: `njones61/xslope` base: `main`
   - Head repository: `YOUR_USERNAME/xslope` compare: `feature/your-feature-name`

6. **Click "Create Pull Request"**

### Pull Request Description

Write a clear description of your changes:

```markdown
## Description
Brief description of what this PR does and why.

## Changes Made
- Added function X to handle Y
- Modified function Z to improve performance
- Updated documentation for feature A

## Testing
Describe how you tested your changes:
- Ran test_griff_ex1.py successfully
- Tested with custom input file
- Verified backwards compatibility with existing examples

## Related Issues
Fixes #123
Related to #456

## Checklist
- [ ] Code follows project style guidelines
- [ ] Comments added for complex sections
- [ ] Documentation updated
- [ ] Tests pass
- [ ] Backwards compatible (or breaking changes documented)
```

### Pull Request Best Practices

- **One feature per PR** - Keep PRs focused on a single change
- **Small PRs are better** - Easier to review and less likely to have conflicts
- **Write descriptive titles** - "Add Spencer's method convergence check" not "Update solve.py"
- **Be responsive** - Address reviewer comments promptly
- **Keep it up to date** - Rebase or merge main if your branch falls behind

## Review Process

### What to Expect

1. **Automated checks** may run (if configured)
2. **Maintainer review** - A project maintainer will review your code
3. **Feedback** - You may receive comments or change requests
4. **Iteration** - Make requested changes and push updates
5. **Approval** - Once approved, your PR will be merged

### Addressing Review Comments

1. **Read feedback carefully** and ask questions if unclear
2. **Make requested changes** in your branch
3. **Commit and push** the updates:

```bash
git add <modified_files>
git commit -m "Address review comments: fix X and update Y"
git push origin feature/your-feature-name
```

4. **Respond to comments** - Let reviewers know you've addressed their feedback

### After Your PR is Merged

1. **Delete your feature branch** (optional but recommended):

```bash
# Delete local branch
git branch -d feature/your-feature-name

# Delete remote branch
git push origin --delete feature/your-feature-name
```

2. **Update your main branch**:

```bash
git checkout main
git pull upstream main
git push origin main
```

3. **Celebrate!** 🎉 You've contributed to xslope!

## Types of Contributions

### Bug Reports

Found a bug? Please open an issue:
- Use a clear title
- Describe the bug and how to reproduce it
- Include error messages and tracebacks
- Specify your Python version and OS
- Provide a minimal example if possible

### Feature Requests

Have an idea for a new feature?
- Open an issue to discuss it first
- Explain the use case and why it would be valuable
- Be open to feedback and alternative approaches

### Code Contributions

Areas where contributions are especially welcome:
- Bug fixes
- Performance improvements
- Additional solution methods
- Enhanced visualization options
- Additional test cases and examples
- Documentation improvements
- Code refactoring for clarity

### Documentation Contributions

Improvements to documentation are always appreciated:
- Fix typos and clarify confusing sections
- Add examples and tutorials
- Improve API documentation
- Translate documentation (if applicable in future)

## Getting Help

If you need assistance:
- **Open an issue** on GitHub for bugs or questions
- **Check existing documentation** at https://xslope.readthedocs.io/
- **Review closed PRs** for examples of similar changes
- **Read CLAUDE.md** for project architecture overview

## License

By contributing to xslope, you agree that your contributions will be licensed under the Apache License 2.0, the same license as the project.

## Questions?

Don't hesitate to ask! Open an issue if you're unsure about anything. We're here to help and appreciate your interest in improving xslope.

Thank you for contributing! 🙏
