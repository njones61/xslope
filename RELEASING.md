# Releasing a new version of xslope

The version number is the single label shared by PyPI, the git tag, the GitHub
Release, and the Zenodo archive. Keep them identical. Three things consume that
version independently, so a release is not done until all the relevant steps are
finished.

## When to do each step

PyPI and the GitHub Release/DOI serve different audiences, and they do **not**
have to move together:

- **PyPI is for users.** Any change you want `pip install xslope` to pick up
  needs a PyPI upload (step 2), which needs a fresh version number (step 1).
  This includes packaging, docs, and metadata fixes.
- **The GitHub Release + Zenodo DOI is for citation.** Mint a DOI only for
  versions someone might cite in a paper to reproduce results — i.e. releases
  that **change computed results or add real capability** (a new analysis
  method, a fix that changes factors of safety, a meaningful feature). Skip the
  GitHub Release (step 3) for packaging, docs, typos, and metadata-only fixes.

Rule of thumb: minor/major bumps usually get a DOI; patch bumps usually do not —
but decide by *what changed*, not the version number. The concept DOI in
`CITATION.cff` always resolves to the newest release, so anyone who just wants
to "cite xslope" has a working link even for versions where no per-version DOI
was minted.

## Steps

1. **Bump the version.** Edit `xslope/_version.py` (`__version__`), e.g.
   `0.1.51` -> `0.1.52`. Commit and push. PyPI will not accept a number that has
   already been used, so this bump is mandatory for every release.

2. **Publish to PyPI** (makes `pip install xslope` get the new code; does NOT
   touch Zenodo):
   ```
   python -m build
   twine upload dist/*
   ```

3. **Create a GitHub Release** (THIS is what mints the Zenodo DOI):
   - GitHub -> Releases -> Draft a new release.
   - New tag `vX.Y.Z` on `main` (match `_version.py`).
   - Title `xslope X.Y.Z`, add notes, Publish.
   - Zenodo (integration already enabled for `njones61/xslope`) archives the
     release zip and mints a **new version DOI** automatically within a minute.

## DOIs

- Each release gets its own **version DOI**. Anything citing a specific,
  validated release (e.g. a paper) should pin that version DOI.
- The **concept DOI** ("Cite all versions" on the Zenodo record) always resolves
  to the newest version. Put the concept DOI in `CITATION.cff` once a second
  version exists, so the citation metadata tracks the latest release.
- Update the `version`, `date-released`, and `doi` fields in `CITATION.cff` after
  each release.

## Notes

- PyPI and Zenodo are independent: a version bump alone does nothing for either;
  PyPI needs the `twine upload`; Zenodo needs the GitHub Release.
- The Zenodo integration only archives releases created while the repo toggle is
  ON in Zenodo's GitHub settings. It is currently ON.
- `CITATION.cff` edits made after a release tag are not in that release's
  archive; that is expected and harmless.
