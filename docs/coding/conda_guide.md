# Conda Publishing Guide for xslope

This guide provides instructions for publishing the xslope package to conda-forge, the community-driven conda package repository.

## Overview

Unlike PyPI, publishing to conda requires submitting a recipe to conda-forge, which is a community-maintained collection of conda packages. The conda-forge team reviews and maintains the infrastructure for building packages across multiple platforms.

## Prerequisites

- Your package must already be published on PyPI
- GitHub account
- Basic knowledge of YAML configuration files

## Option 1: Publishing to conda-forge (Recommended)

Conda-forge is the preferred method for distributing conda packages. It provides automated builds across multiple platforms and handles dependencies well.

### Step 1: Fork the Staged Recipes Repository

1. Go to https://github.com/conda-forge/staged-recipes
2. Click "Fork" to create your own copy
3. Clone your fork locally:

```bash
git clone https://github.com/YOUR_USERNAME/staged-recipes.git
cd staged-recipes
```

### Step 2: Create Your Recipe

1. Create a new directory for your package:

```bash
cd recipes
mkdir xslope
cd xslope
```

2. Create a `meta.yaml` file with the following content:

```yaml
{% set name = "xslope" %}
{% set version = "0.1.0" %}

package:
  name: {{ name|lower }}
  version: {{ version }}

source:
  url: https://pypi.io/packages/source/{{ name[0] }}/{{ name }}/xslope-{{ version }}.tar.gz
  sha256: REPLACE_WITH_SHA256_HASH

build:
  noarch: python
  script: {{ PYTHON }} -m pip install . -vv
  number: 0

requirements:
  host:
    - python >=3.9
    - pip
    - setuptools >=68
  run:
    - python >=3.9
    - numpy
    - pandas
    - matplotlib
    - scipy
    - shapely
    - openpyxl

test:
  imports:
    - xslope
  commands:
    - pip check
  requires:
    - pip

about:
  home: https://github.com/njones61/xslope
  summary: Slope stability analysis (limit equilibrium and FEM) in Python.
  description: |
    xslope is a Python package for limit equilibrium slope stability analysis.
    It provides multiple solution methods including Ordinary Method of Slices,
    Bishop's Simplified Method, Janbu Method, Spencer's Method, and more. The
    package also includes seepage analysis capabilities and integrates with
    finite element mesh analysis.
  license: Apache-2.0
  license_family: APACHE
  license_file: LICENSE
  doc_url: https://xslope.readthedocs.io/en/latest/
  dev_url: https://github.com/njones61/xslope

extra:
  recipe-maintainers:
    - njones61
```

### Step 3: Get the SHA256 Hash

You need to get the SHA256 hash of your PyPI source distribution:

```bash
# Download the source from PyPI
wget https://pypi.io/packages/source/x/xslope/xslope-0.1.0.tar.gz

# Calculate SHA256 hash
# On macOS/Linux:
shasum -a 256 xslope-0.1.0.tar.gz

# On Windows (PowerShell):
# Get-FileHash xslope-0.1.0.tar.gz -Algorithm SHA256
```

Copy the hash and replace `REPLACE_WITH_SHA256_HASH` in the `meta.yaml` file.

### Step 4: Test Your Recipe Locally (Optional)

Install conda-build and test your recipe:

```bash
conda install conda-build
conda build recipes/xslope
```

This will build the package locally and run tests.

### Step 5: Submit Pull Request

1. Commit your changes:

```bash
git add recipes/xslope/meta.yaml
git commit -m "Add recipe for xslope"
git push origin main
```

2. Go to https://github.com/conda-forge/staged-recipes
3. Click "New Pull Request"
4. Select your fork and branch
5. Submit the pull request with a clear description

### Step 6: Address Review Comments

The conda-forge team and automated bots will review your submission:

- **@conda-forge-admin**: Automated checks for recipe format
- **@conda-forge/help-python**: Human reviewers for Python packages

Common issues to watch for:
- Missing dependencies
- Incorrect license specification
- Test failures
- Platform-specific build issues

Address any comments and push updates to your branch. The PR will automatically update.

### Step 7: Merge and Activation

Once approved:
1. A conda-forge team member will merge your PR
2. The package will be automatically built for all platforms
3. Within a few hours, your package will be available via:

```bash
conda install -c conda-forge xslope
```

### Step 8: Become a Maintainer

After your recipe is merged, you'll be added as a maintainer of the xslope-feedstock repository. This allows you to:
- Publish updates
- Manage package configuration
- Review PRs for your package

## Updating Your Conda Package

Once your package is on conda-forge, updates are managed through the feedstock repository.

### Automated Updates with Regro-cf-autotick-bot

Conda-forge has bots that automatically detect new PyPI releases and create PRs:

1. Publish new version to PyPI
2. Wait a few hours for the bot to detect the new version
3. A PR will be automatically created at `https://github.com/conda-forge/xslope-feedstock`
4. Review the PR and merge if everything looks good

### Manual Updates

If you need to update manually:

1. Clone the feedstock repository:

```bash
git clone https://github.com/conda-forge/xslope-feedstock.git
cd xslope-feedstock
```

2. Create a new branch:

```bash
git checkout -b update_version
```

3. Edit `recipe/meta.yaml`:
   - Update the version number
   - Update the SHA256 hash
   - Increment the build number if making non-version changes

4. Commit and push:

```bash
git add recipe/meta.yaml
git commit -m "Update to version 0.1.1"
git push origin update_version
```

5. Create a PR on GitHub
6. Once CI passes, merge the PR

## Option 2: Publishing to Your Own Conda Channel

If you want to distribute conda packages without going through conda-forge:

### Step 1: Build the Package

```bash
conda install conda-build anaconda-client

# Build the package
conda build recipes/xslope

# Find the built package
conda build recipes/xslope --output
```

### Step 2: Upload to Anaconda Cloud

```bash
# Login to Anaconda Cloud
anaconda login

# Upload the package
anaconda upload /path/to/built/package.tar.bz2
```

### Step 3: Install from Your Channel

Users can install with:

```bash
conda install -c YOUR_USERNAME xslope
```

**Note:** This approach requires you to build packages for each platform (linux, osx, windows) separately or use `noarch: python` for pure Python packages.

## Best Practices

### Recipe Maintenance

- **Keep dependencies minimal**: Only include runtime requirements
- **Use version constraints wisely**: Be specific only when necessary
- **Test thoroughly**: Ensure the package works on all supported platforms
- **Document changes**: Use clear commit messages

### Version Updates

When updating to a new version:

1. Update version number in `meta.yaml`
2. Update SHA256 hash
3. Reset build number to 0
4. Update dependencies if needed
5. Test the build locally before submitting

### Build Numbers

Increment the build number when:
- Fixing recipe bugs without version change
- Updating dependencies
- Adding new build variants
- Fixing platform-specific issues

Reset to 0 when publishing a new version.

## Troubleshooting

### Recipe Validation Errors

If the conda-forge linter fails:
- Check YAML syntax
- Verify all required fields are present
- Ensure license is correctly specified
- Check that maintainer GitHub usernames are correct

### Build Failures

Common causes:
- Missing dependencies in the recipe
- Test imports failing (package not properly installed)
- Platform-specific issues (especially Windows)

### SHA256 Mismatch

If you get SHA256 errors:
- Ensure you're downloading from the correct URL
- Verify you're using the `.tar.gz` source distribution, not the wheel
- Recalculate the hash and update `meta.yaml`

### Import Errors During Testing

If `test: imports:` fails:
- Verify the package name is correct
- Check that all dependencies are listed in `requirements: run:`
- Test locally with `conda build`

## Resources

- Conda-forge Documentation: https://conda-forge.org/docs/
- Conda Build Documentation: https://docs.conda.io/projects/conda-build/
- Staged Recipes Repository: https://github.com/conda-forge/staged-recipes
- Anaconda Cloud: https://anaconda.org/

## Quick Reference

```bash
# Get SHA256 hash from PyPI tarball
wget https://pypi.io/packages/source/x/xslope/xslope-VERSION.tar.gz
shasum -a 256 xslope-VERSION.tar.gz

# Test recipe locally
conda build recipes/xslope

# Install from conda-forge (after package is published)
conda install -c conda-forge xslope

# Update to latest version
conda update -c conda-forge xslope
```

## Timeline Expectations

- **Initial submission review**: 1-7 days
- **After merge, package availability**: 2-6 hours
- **Automated version updates**: Usually within 24 hours of PyPI release
- **Manual update PR review**: 1-3 days

## Support

If you encounter issues:
- Check conda-forge documentation: https://conda-forge.org/docs/
- Ask in the conda-forge Gitter channel: https://gitter.im/conda-forge/conda-forge.github.io
- Open an issue in staged-recipes (before recipe is merged)
- Open an issue in your feedstock repository (after recipe is merged)
