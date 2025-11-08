# PyPI Publishing Guide for xslope

This guide provides step-by-step instructions for publishing the xslope package to PyPI, including testing on TestPyPI first.

## Prerequisites

Install the required build and upload tools:

```bash
python -m pip install --upgrade build twine
```

## Step 1: Test Upload to TestPyPI

TestPyPI is a separate instance of PyPI for testing and experimentation. Always test here first before publishing to the real PyPI.

### 1.1 Create TestPyPI Account

1. Go to https://test.pypi.org/account/register/
2. Create an account and verify your email
3. Set up two-factor authentication (recommended)

### 1.2 Generate API Token for TestPyPI

1. Go to https://test.pypi.org/manage/account/token/
2. Click "Add API token"
3. Name it (e.g., "xslope-test")
4. Set scope to "Entire account" (or specific to xslope project once created)
5. Copy the token (starts with `pypi-`) - **save it somewhere safe, you won't see it again**

### 1.3 Build the Package

Clean any previous builds and create new distribution files:

```bash
# Clean previous builds
rm -rf dist/ build/ *.egg-info

# Build the package
python -m build
```

This creates two files in the `dist/` directory:
- `xslope-0.1.0.tar.gz` (source distribution)
- `xslope-0.1.0-py3-none-any.whl` (wheel distribution)

### 1.4 Upload to TestPyPI

```bash
python -m twine upload --repository testpypi dist/*
```

When prompted:
- Username: `__token__`
- Password: paste your TestPyPI API token (including the `pypi-` prefix)

### 1.5 Test Installation from TestPyPI

Create a new virtual environment and test install:

```bash
# Create and activate a new virtual environment
python -m venv test_env
source test_env/bin/activate  # On Windows: test_env\Scripts\activate

# Install from TestPyPI (need to use --extra-index-url for dependencies)
pip install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple/ xslope

# Test the import
python -c "import xslope; from xslope._version import __version__; print(f'xslope version: {__version__}')"

# Deactivate when done
deactivate
```

**Note:** The `--extra-index-url` flag is needed because TestPyPI doesn't contain all dependencies (numpy, pandas, etc.), so pip needs to fall back to the real PyPI for those.

## Step 2: Upload to Production PyPI

Once you've verified everything works on TestPyPI, you can publish to the real PyPI.

### 2.1 Create PyPI Account

1. Go to https://pypi.org/account/register/
2. Create an account and verify your email
3. Set up two-factor authentication (required for new projects)

### 2.2 Generate API Token for PyPI

1. Go to https://pypi.org/manage/account/token/
2. Click "Add API token"
3. Name it (e.g., "xslope-prod")
4. Set scope to "Entire account" initially (you can create a project-specific token after first upload)
5. Copy the token - **save it securely**

### 2.3 Upload to PyPI

```bash
python -m twine upload dist/*
```

When prompted:
- Username: `__token__`
- Password: paste your PyPI API token

### 2.4 Verify the Upload

1. Check your package page: https://pypi.org/project/xslope/
2. Test installation:

```bash
pip install xslope
python -c "import xslope; from xslope._version import __version__; print(f'xslope version: {__version__}')"
```

## Step 3: Publishing Updates

When you make changes to the code and want to release a new version:

### 3.1 Update the Version Number

Edit `xslope/_version.py`:

```python
__all__ = ["__version__"]
__version__ = "0.1.1"  # Increment version following semantic versioning
```

**Semantic Versioning Guide:**
- `MAJOR.MINOR.PATCH` (e.g., `1.2.3`)
- **MAJOR**: Incompatible API changes (breaking changes)
- **MINOR**: Add functionality in a backwards-compatible manner (new features)
- **PATCH**: Backwards-compatible bug fixes

Examples:
- Bug fix: `0.1.0` → `0.1.1`
- New feature: `0.1.1` → `0.2.0`
- Breaking change: `0.2.0` → `1.0.0`

### 3.2 Update CHANGELOG (Optional but Recommended)

Create/update a `CHANGELOG.md` file documenting what changed in each version.

### 3.3 Commit Version Changes

```bash
git add xslope/_version.py
git commit -m "Bump version to 0.1.1"
git tag -a v0.1.1 -m "Release version 0.1.1"
git push origin main --tags
```

### 3.4 Clean and Rebuild

```bash
# Remove old builds
rm -rf dist/ build/ xslope.egg-info/

# Rebuild
python -m build
```

### 3.5 Upload the New Version

Test on TestPyPI first (optional but recommended):

```bash
python -m twine upload --repository testpypi dist/*
```

Then upload to production PyPI:

```bash
python -m twine upload dist/*
```

## Configuration Files (Optional)

### Save Credentials with .pypirc

Instead of entering credentials each time, you can create a `~/.pypirc` file:

```ini
[distutils]
index-servers =
    pypi
    testpypi

[pypi]
username = __token__
password = pypi-YOUR_ACTUAL_PYPI_TOKEN_HERE

[testpypi]
repository = https://test.pypi.org/legacy/
username = __token__
password = pypi-YOUR_ACTUAL_TESTPYPI_TOKEN_HERE
```

**Important:** Keep this file secure! Set proper permissions:

```bash
chmod 600 ~/.pypirc
```

Add `~/.pypirc` to your `.gitignore` if it's not already ignored globally.

## Troubleshooting

### Error: File already exists

You cannot upload the same version twice to PyPI. You must increment the version number in `xslope/_version.py` before uploading again.

### Error: Invalid authentication credentials

- Double-check your API token
- Make sure you're using `__token__` as the username (with double underscores)
- Ensure you copied the entire token including the `pypi-` prefix

### Error: Package name already taken

If you're uploading for the first time and the name is taken, you'll need to choose a different name. Check availability at https://pypi.org/project/your-package-name/

### Build errors

If `python -m build` fails, check:
- `pyproject.toml` syntax is valid
- All required files exist (README.md, LICENSE, etc.)
- The version in `xslope/_version.py` is accessible

## Quick Reference Commands

```bash
# Clean, build, and upload to TestPyPI
rm -rf dist/ build/ *.egg-info && python -m build && python -m twine upload --repository testpypi dist/*

# Clean, build, and upload to PyPI
rm -rf dist/ build/ *.egg-info && python -m build && python -m twine upload dist/*

# Test local installation
pip install -e .

# Uninstall
pip uninstall xslope
```

## Additional Resources

- PyPI Help: https://pypi.org/help/
- Python Packaging Guide: https://packaging.python.org/
- Semantic Versioning: https://semver.org/
- Twine Documentation: https://twine.readthedocs.io/
