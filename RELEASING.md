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

Every tag `vX.Y.Z` pushed to GitHub runs the installer workflow, which publishes
the installers as a GitHub Release; Zenodo archives every GitHub Release and
mints a version DOI for it. So every released version, patch or not, carries its
own DOI (accepted 2026-09-09). The concept DOI in `CITATION.cff` always resolves
to the newest release; update the `version`, `date-released` and `doi` fields
after each release once Zenodo shows the new record.

## Before a minor or major release

- [ ] Run the tutorial execution pass: every tutorial page executed as written,
      each printed number compared to what the run produces, each linked file
      and image checked. The brief that scopes it is `BRIEF.md` under
      `reports/campaign_1_0_0_readiness/tutorial_pass/` in the private repo.
      `python run_tests.py --tutorials` is the cheap standing layer under it and
      is not a substitute: it reads the numbers a tag can reach, not the steps.

Patch releases skip this.

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

3. **Tag the release** — `git tag -a vX.Y.Z -m "xslope X.Y.Z"` and push the tag.
   The push runs `.github/workflows/release-installers.yml`, which builds and signs
   the macOS and Windows installers and publishes them as the GitHub Release
   `XSLOPE vX.Y.Z` (with `latest.json` for Studio's update check). Zenodo
   (integration enabled for `njones61/xslope`) archives that release and mints
   the version DOI within a minute. Add release notes on the GitHub Release
   afterwards if wanted. If the macOS leg fails on Apple's timestamp service the
   step retries on its own; re-run the workflow from the Actions tab otherwise.

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
