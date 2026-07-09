# Purging reference PDFs from the public repo's git history

**Status:** 🟡 **Planned, not executed.** Reviewed and approved in principle 2026-07-09;
the rewrite itself is a deliberate, destructive step and has not been run. *Not legal
advice.*

## 1. Why

34 third-party PDFs — journal articles under publisher copyright (Géotechnique/ICE, ASCE,
Elsevier), vendor documents (Rocscience), and agency manuals — were committed to the public
`njones61/xslope` repo. Commit `d779408` removed them from `HEAD` and moved them to the
private sibling repo `../xslope_private/ref_docs/`.

**That is not enough.** Git history is immutable and additive: every one of those blobs is
still reachable from older commits, and GitHub serves blobs by SHA indefinitely. Anyone can
run `git log --all -- '*.pdf'` on a clone, or hit
`github.com/njones61/xslope/raw/<old-sha>/…`, and get the file. Until history is rewritten,
the repo still publishes them.

This matters now because the package is heading for a broader public release and an academic
paper. It is also simply inconsistent: `plan_import.md` goes to considerable lengths to avoid
republishing one vendor's manuals while a stack of paywalled journal PDFs sits in history.

## 2. Facts established (2026-07-09)

- **34 distinct PDF paths** across all history.
- **They lived under two different prefixes.** Originally `docs/…`, later moved to
  `ref_docs/…` by commit `cf8cdbf` ("misc reorganization of repo"). **A purge keyed only on
  `ref_docs/` would miss every `docs/…` copy.** Both must be targeted.
- `.git` is **301 MB**; GitHub reports `diskUsage` 223 MB. Nearly all of it is these PDFs.
  The rewrite should shrink the repo by roughly an order of magnitude.
- **3 forks, 3 stars.** Forks are the leak (see §5).
- **1 tag** (`v0.1.51`) with a GitHub release. The tag will be rewritten and must be re-pushed.
- ✅ **PyPI is clean.** `MANIFEST.in` includes only `LICENSE` and `README.md`, and
  `packages = ["xslope", "studio"]`. No sdist or wheel ever carried a PDF. Nothing to do there.
- ⚠️ **Read the Docs is an open question.** `mkdocs.yml` sets
  `site_url: https://xslope.readthedocs.io/en/latest/`, and the PDFs sat *inside* `docs/`
  before `cf8cdbf`. If MkDocs copied them into the built site, they were published there too.
  **Check RTD's older builds/versions before declaring this done.**
- `git-filter-repo` is **not installed** (`pip install git-filter-repo`).

## 3. Not everything here is equally exposed

Worth sorting before deciding how aggressive to be. Several of these are US federal works or
public agency releases and carry little or no restriction:

- **Low concern:** `usace_slope_stability_manual.pdf`, the UTEXAS manuals (`ADA207044`,
  `IR-GL-87-1-Volume-4`), `stabl5_documentation.pdf`, `s2dprimer.pdf`,
  `beidenharn tracy - seep2d - 1987.pdf`.
- **The actual exposure:** `Griffiths and Lane 1999.pdf`, `dawson et al - geot_1999…pdf`,
  the seven-part `smith and griffiths 1998` series, `original spencer paper.pdf`,
  `potts 1990 progressive failure.pdf`, `ito_matsui_paper.pdf`, `reinforcement_duncan_wright.pdf`,
  and the Rocscience article.

Purging all 34 is simpler than curating, costs nothing extra, and the low-concern documents are
freely re-downloadable from their agencies anyway. **Recommend: purge all PDFs.**

## 4. Procedure

Run from a scratch directory. **Do not run this in the working repo.**

```bash
pip install git-filter-repo

# 1. Back up. Non-negotiable — this is irreversible.
git clone --mirror https://github.com/njones61/xslope.git xslope-backup.git
cp -a xslope-backup.git xslope-rewrite.git   # keep the backup pristine

# 2. Enumerate exactly what will be dropped, and eyeball it.
cd xslope-rewrite.git
git log --all --pretty=format: --name-only --diff-filter=A \
  | grep -i '\.pdf$' | sort -u

# 3. Rewrite. A filename-callback is safer than --path-glob here, because the
#    files lived under BOTH docs/ and ref_docs/ at different times.
git filter-repo --force --filename-callback '
  if filename.endswith(b".pdf") and (
      filename.startswith(b"ref_docs/") or filename.startswith(b"docs/")):
      return None
  return filename
'

# 4. Verify: expect zero hits and a much smaller repo.
git log --all --pretty=format: --name-only | grep -i '\.pdf$' | sort -u   # -> empty
du -sh .

# 5. Confirm the code that must survive is still there.
git log --all --oneline -- ref_docs/ref_docs_seep/seep2d_fortran/src/seep2d.f | tail -1
```

Only once every check passes:

```bash
# 6. Push the rewritten history and tags. Destructive.
git push --force --mirror https://github.com/njones61/xslope.git
```

Then re-clone fresh locally; **do not** keep working in the pre-rewrite clone.

## 5. What the rewrite does *not* fix, and what to do about it

- **Forks (3 of them).** GitHub stores forks in a shared object network. Force-pushing your
  copy does not touch a fork's refs, and unreachable objects may still be served from the
  network. **Open a ticket with GitHub Support** asking them to garbage-collect unreachable
  objects and to confirm the fork network is clean. Forks that have their own copy of the old
  commits must be deleted by their owners — you cannot do it for them.
- **Existing clones.** Anyone who cloned already has the PDFs. Nothing recovers those.
- **Caches and archives.** GitHub's `codeload` archives, Software Heritage, and the Wayback
  Machine may hold copies. Software Heritage archives public repos wholesale and honors
  takedown requests; check `archive.softwareheritage.org` for `njones61/xslope`.
- **Read the Docs.** See §2. If the PDFs were published there, purge the affected builds and
  versions in the RTD admin.
- **Collaborators.** The tag `v0.1.51` and every commit SHA change. Anyone with a clone must
  re-clone. In practice this is just the author.

## 6. Decision gate

The rewrite is cheap to plan and expensive to undo. Before running §4:

1. Confirm the 3 forks are yours or dormant, and decide whether to ask their owners to delete.
2. Check Read the Docs for published copies.
3. Decide whether the release `v0.1.51` needs recreating after the tag is rewritten.
4. Take the mirror backup and keep it somewhere durable and **private** — it still contains
   every PDF.

A reasonable alternative, if the rewrite feels disproportionate: leave history alone and accept
the exposure. Publishers rarely pursue individual academics over a handful of PDFs in a small
repo. But the exposure is real, the fix is well-understood, and the repo is about to get much
more visible — so the recommendation is to do it, once, before the paper.
