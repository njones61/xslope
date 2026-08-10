"""Standing checks on the docs' sample links and the ``xslope://`` scheme behind them.

Two halves of one road: the docs build packages every sample and links it twice, and
Studio answers the link. Both halves have failure modes nobody would see by looking:

  A. THE VERB GATE. One verb is registered, ``open``, and the handler must refuse
     everything else BY NAME rather than guess — a future ``xslope://docs?...`` link
     clicked on an old build has to say so, not do something surprising with it. The
     gate is mutation-tested: widen ``VERBS`` and the same link is accepted, which is
     how we know the refusal comes from the gate and not from a parse failure.
  B. THE ALLOWLIST. The scheme is registered system-wide, so any page anywhere can
     put an ``xslope://`` link in front of a user. Only the docs site and the project
     repository may be fetched from, only over https, and never a local path — a
     ``file:`` URL through this door would let a web page hand Studio a file on the
     user's disk. Mutation-tested the same way: add the hostile host to
     ``ALLOWED_HOSTS`` and it is accepted, so the constant is the whole decision.
  C. REDIRECTS ARE RE-CHECKED. ``urlopen`` follows redirects itself, so a check made
     once on the URL in the link says nothing about where the bytes came from. An
     open redirect on an allowlisted host is the standard way around a bare check.
  D. THE SAVED NAME COMES FROM THE URL'S PATH, and is a plain file name — a package
     downloaded from a link may not name a directory, a parent, or an extension of
     the server's choosing.
  E. THE DOCS BUILD PACKAGES AND PAIRS. A real MkDocs build over a scratch docs tree
     has to produce a ``.xslz`` per sample workbook, holding its sidecars, and a
     Download · Open in Studio pair per link whose two hrefs name the SAME package —
     the one property the whole "emit both from one URL" design exists for. The
     escape hatch (``.raw-file``) and non-sample links must come through untouched,
     and the build must refuse to finish if a page ever links a package nobody built.
  F. STUDIO'S FLOW, IN ORDER. Refused links never reach the network; the confirmation
     comes BEFORE the download, and cancelling it downloads nothing; and what a link
     opens is the ordinary unpack-then-open path, so the document ends up on the
     extracted workbook rather than on the package.

No network: the download is stubbed in F and the redirect handler is exercised
directly in C. The MkDocs build in E writes only into a temporary directory.
"""
import contextlib
import io
import os
import re
import shutil
import sys
import tempfile

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(_HERE)
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import xslope                                                    # noqa: E402
from xslope.package import PACKAGE_EXT, package_contents         # noqa: E402
from studio import urlscheme                                     # noqa: E402

#: A real project with sidecars — the scratch docs tree is built from it, so the
#: packaged samples in E hold a mesh and a solved field rather than a lone workbook.
WITH_SIDECARS = os.path.join(_REPO, "docs/seep/files/xslope_rface_SEEP_KEY.xlsx")
#: A project that is a workbook and nothing else (single-file samples package too).
SINGLE_FILE = os.path.join(_REPO, "docs/inputs/slope/xslope_simple1.xlsx")

#: A URL on the docs site — what the docs build emits.
GOOD = "https://xslope.readthedocs.io/en/latest/lem/files/xslope_simple1.xslz"
#: The docs URL of the sidecar-carrying sample the Studio flow leg serves.
FLOW = ("https://xslope.readthedocs.io/en/latest/seep/files/"
        + os.path.splitext(os.path.basename(WITH_SIDECARS))[0] + PACKAGE_EXT)

_TMPDIRS = []


def _quiet(fn, *a, **k):
    with contextlib.redirect_stdout(io.StringIO()):
        return fn(*a, **k)


def _tmp():
    d = tempfile.mkdtemp(prefix="xslz_links_")
    _TMPDIRS.append(d)
    return d


# ------------------------------------------------------------------- A. the verb
def test_verb_gate():
    """One verb, and everything else refused by name."""
    fails = []
    verb, url = urlscheme.parse_request(urlscheme.build_url(GOOD))
    if (verb, url) != ("open", GOOD):
        fails.append(f"the docs' own link parses to {(verb, url)}, not ('open', the URL)")

    # (link, what the message must say, whether it must name the verb we DO know)
    for bad, wanted, names_verb in [
        (f"xslope://docs?url={GOOD}", "docs", True),
        (f"xslope://export?url={GOOD}", "export", True),
        (f"xslope:import?url={GOOD}", "import", True),   # the schemeless spelling
        ("xslope://", "action", True),
        ("xslope://open", "url=", False),
        (f"https://xslope.readthedocs.io/?url={GOOD}", "not an xslope: link", False),
    ]:
        try:
            urlscheme.parse_request(bad)
        except urlscheme.SchemeError as exc:
            if wanted not in str(exc).lower():
                fails.append(f"the refusal of {bad!r} does not say {wanted!r}: {exc}")
            if names_verb and "'open'" not in str(exc):
                fails.append(f"the refusal of {bad!r} does not name the verb it does "
                             f"understand: {exc}")
        else:
            fails.append(f"{bad!r} was accepted")

    # MUTATION: the gate is what refuses, not the parser. Widen it and the same
    # link parses; the URL comes through untouched.
    saved = urlscheme.VERBS
    urlscheme.VERBS = ("open", "docs")
    try:
        verb, url = urlscheme.parse_request(f"xslope://docs?url={GOOD}")
    except urlscheme.SchemeError as exc:
        fails.append(f"with 'docs' in VERBS the link is still refused ({exc}), so the "
                     f"refusal is not the verb gate")
    else:
        if (verb, url) != ("docs", GOOD):
            fails.append(f"a widened gate parsed {(verb, url)}")
    finally:
        urlscheme.VERBS = saved
    if urlscheme.VERBS != ("open",):
        fails.append(f"this build registers {urlscheme.VERBS}; the ruling is 'open' only")
    return fails


# -------------------------------------------------------------- B. the allowlist
def test_allowlist():
    """Only XSLOPE's own sites, only https, never a local file."""
    fails = []
    for host in urlscheme.ALLOWED_HOSTS:
        url = f"https://{host}/x{PACKAGE_EXT}"
        try:
            urlscheme.check_url(url)
        except urlscheme.SchemeError as exc:
            fails.append(f"an allowlisted host was refused: {exc}")

    hostile = [
        "https://evil.example/x.xslz",
        "https://xslope.readthedocs.io.evil.example/x.xslz",   # suffix, not the host
        "https://evilxslope.readthedocs.io/x.xslz",            # nor a prefix
        "https://xslope.readthedocs.io@evil.example/x.xslz",   # userinfo, not the host
        "http://xslope.readthedocs.io/x.xslz",                 # https only
        "file:///Users/me/secret.xslz",                        # never a local path
        "data:application/zip;base64,UEsDBA==",
        "/Users/me/secret.xslz",
    ]
    for url in hostile:
        try:
            urlscheme.check_url(url)
        except urlscheme.SchemeError as exc:
            if url not in str(exc):
                fails.append(f"the refusal of {url} does not name it: {exc}")
        else:
            fails.append(f"{url} was accepted for download")

    # MUTATION: the list is the whole decision. Add the hostile host and it passes.
    saved = urlscheme.ALLOWED_HOSTS
    urlscheme.ALLOWED_HOSTS = saved + ("evil.example",)
    try:
        urlscheme.check_url("https://evil.example/x.xslz")
    except urlscheme.SchemeError as exc:
        fails.append(f"with evil.example allowlisted the URL is still refused ({exc}), "
                     f"so ALLOWED_HOSTS is not what refuses it")
    finally:
        urlscheme.ALLOWED_HOSTS = saved
    if "xslope.readthedocs.io" not in urlscheme.ALLOWED_HOSTS:
        fails.append("the docs site is not on the allowlist, so no docs link works")
    return fails


# --------------------------------------------------------------- C. redirects
def test_redirect_guard():
    """A redirect off an allowlisted host is refused, in the handler itself."""
    import urllib.request

    fails = []
    handler = urlscheme._AllowlistRedirects()
    req = urllib.request.Request(GOOD)
    headers = {}

    try:
        handler.redirect_request(req, None, 302, "Found", headers,
                                 "https://evil.example/x.xslz")
    except urlscheme.SchemeError as exc:
        if "evil.example" not in str(exc):
            fails.append(f"the redirect refusal does not name the host: {exc}")
    else:
        fails.append("a redirect to evil.example was followed")

    same = "https://xslope.readthedocs.io/en/stable/lem/files/x.xslz"
    try:
        out = handler.redirect_request(req, None, 302, "Found", headers, same)
    except urlscheme.SchemeError as exc:
        fails.append(f"a redirect within the docs site was refused: {exc}")
    else:
        if out is None or out.full_url != same:
            fails.append("a redirect within the docs site was not followed")
    return fails


# ------------------------------------------------------------ D. the saved name
def test_saved_name():
    """What a link may name on the user's disk: a plain package file, nothing else."""
    fails = []
    cases = {
        GOOD: "xslope_simple1.xslz",
        "https://xslope.readthedocs.io/a/b/model.xslz?v=2": "model.xslz",
        "https://xslope.readthedocs.io/a/%2E%2E%2F%2E%2E%2Fpasswd.xslz": "passwd.xslz",
        "https://xslope.readthedocs.io/a/b/model.exe": "model.xslz",
        "https://xslope.readthedocs.io/": "project.xslz",
        "https://xslope.readthedocs.io/a/..": "project.xslz",
    }
    for url, want in cases.items():
        got = urlscheme.package_name(url)
        if got != want:
            fails.append(f"{url} would be saved as {got!r}, not {want!r}")
        if os.path.basename(got) != got or os.path.isabs(got):
            fails.append(f"{url} would be saved under a path, not a name: {got!r}")
    return fails


# ------------------------------------------------------- E. the docs build
def _scratch_docs():
    """A miniature docs tree: two samples, and a page linking them every way."""
    root = _tmp()
    docs = os.path.join(root, "docs")
    files = os.path.join(docs, "lem", "files")
    os.makedirs(files)
    for src in xslope.project_files(WITH_SIDECARS):
        shutil.copy2(src, files)
    shutil.copy2(SINGLE_FILE, files)
    with_sidecars = os.path.splitext(os.path.basename(WITH_SIDECARS))[0]
    single = os.path.splitext(os.path.basename(SINGLE_FILE))[0]

    os.makedirs(os.path.join(docs, "usage"))
    with open(os.path.join(docs, "lem", "samples.md"), "w") as fh:
        fh.write(
            f"# Samples\n\n"
            f"Excel input file: [{with_sidecars}.xlsx](files/{with_sidecars}.xlsx)\n\n"
            f"The run uses [{single}.xlsx](files/{single}.xlsx), which is a workbook "
            f"and nothing else.\n\n"
            f"* in a list: [{single}.xlsx](files/{single}.xlsx)\n\n"
            f"The bare workbook: "
            f"[{single}.xlsx](files/{single}.xlsx){{: .raw-file }}\n\n"
            f"Not a sample: [the template](../usage/input_template.xlsx)\n")
    # A page whose first pair is NOT in a paragraph: there is nowhere to put a note
    # paragraph without breaking the list, so it goes to the top of the page.
    with open(os.path.join(docs, "lem", "listed.md"), "w") as fh:
        fh.write(f"# Listed\n\n* [{single}.xlsx](files/{single}.xlsx)\n"
                 f"* nothing else\n")
    with open(os.path.join(docs, "usage", "index.md"), "w") as fh:
        fh.write("# Usage\n\nNo samples here.\n")
    open(os.path.join(docs, "usage", "input_template.xlsx"), "wb").close()

    site = os.path.join(root, "site")
    hook = os.path.join(_REPO, "hooks", "docs_packages.py").replace("\\", "/")
    with open(os.path.join(root, "mkdocs.yml"), "w") as fh:
        fh.write(f"site_name: Scratch\n"
                 f"site_url: https://xslope.readthedocs.io/en/latest/\n"
                 f"docs_dir: {docs}\n"
                 f"site_dir: {site}\n"
                 f"use_directory_urls: true\n"
                 f"markdown_extensions:\n  - attr_list\n"
                 f"hooks:\n  - {hook}\n")
    return root, docs, site, with_sidecars, single


def test_docs_build():
    """One real MkDocs build: the packages, and the pairs in the rendered HTML."""
    try:
        from mkdocs.commands.build import build as mkdocs_build
        from mkdocs.config import load_config
    except Exception as exc:                       # pragma: no cover
        print(f"  (mkdocs not installed — docs-build leg skipped: {exc})")
        return []

    fails = []
    root, docs, site, with_sidecars, single = _scratch_docs()
    _quiet(mkdocs_build, load_config(os.path.join(root, "mkdocs.yml")))

    # --- the packages ---------------------------------------------------------
    for stem in (with_sidecars, single):
        pkg = os.path.join(site, "lem", "files", stem + PACKAGE_EXT)
        if not os.path.isfile(pkg):
            fails.append(f"the build wrote no package for {stem}")
            continue
        got = sorted(package_contents(pkg))
        want = sorted(os.path.basename(p) for p in
                      xslope.project_files(os.path.join(docs, "lem", "files",
                                                        stem + ".xlsx")))
        if got != want:
            fails.append(f"{stem}{PACKAGE_EXT} holds {got}, the project is {want}")
    if len(xslope.project_files(WITH_SIDECARS)) < 2:
        fails.append("the sidecar fixture has no sidecars, so the leg proves nothing")
    if os.path.exists(os.path.join(docs, "lem", "files", single + PACKAGE_EXT)):
        fails.append("the build wrote a package into docs/ — nothing may be committed")
    if os.path.exists(os.path.join(site, "usage", "input_template" + PACKAGE_EXT)):
        fails.append("a workbook outside a files/ directory was packaged")

    # --- the pairs ------------------------------------------------------------
    html = open(os.path.join(site, "lem", "samples", "index.html")).read()
    pairs = re.findall(
        r'<a class="xslz-download" href="([^"]+)"[^>]*>Download</a> · '
        r'<a class="xslz-studio" href="xslope://open\?url=([^"]+)"[^>]*>'
        r'Open in Studio</a>', html)
    if len(pairs) != 3:
        fails.append(f"the page rendered {len(pairs)} pairs, not 3 (one per "
                     f"non-escaped sample link)")
    for href, scheme_arg in pairs:
        _verb, url = urlscheme.parse_request("xslope://open?url=" + scheme_arg)
        # THE property: both halves of a pair name the same package. Resolved
        # against the page's own URL, the relative href IS the absolute one.
        from urllib.parse import urljoin
        want = urljoin("https://xslope.readthedocs.io/en/latest/lem/samples/", href)
        if url != want:
            fails.append(f"a pair disagrees: Download is {want}, Open in Studio is {url}")
        local = os.path.normpath(os.path.join(site, "lem", "samples", href))
        if not os.path.isfile(local):
            fails.append(f"a Download link points at {href}, which the build did not "
                         f"write")
        if not href.endswith(PACKAGE_EXT):
            fails.append(f"a Download link points at {href}, not a package")

    # --- the escape hatch, and everything that is not a sample link -----------
    if f'href="../files/{single}.xlsx"' not in html:
        fails.append("the .raw-file escape hatch link did not survive the build")
    if html.count(f'href="../files/{single}.xlsx"') != 1:
        fails.append("a link that is not marked .raw-file kept its .xlsx href")
    if 'href="../../usage/input_template.xlsx"' not in html:
        fails.append("a workbook outside files/ was rewritten anyway")

    # --- the note, once, next to the first pair -------------------------------
    if html.count("xslz-note") != 1:
        fails.append(f"the page carries {html.count('xslz-note')} explanatory notes, "
                     f"not exactly one")
    note_at, first_pair = html.find("xslz-note"), html.find("xslz-download")
    if not 0 <= note_at < first_pair:
        fails.append("the note is not above the first pair")
    note = html[note_at:first_pair]
    for wanted in ("pip", "Download", "Open in Studio"):
        if wanted not in note:
            fails.append(f"the note above the first pair does not mention {wanted!r}")
    listed = open(os.path.join(site, "lem", "listed", "index.html")).read()
    if listed.count("xslz-note") != 1:
        fails.append("a page whose first pair is inside a list got no note")
    elif not 0 <= listed.find("xslz-note") < listed.find("xslz-download"):
        fails.append("the note on the list page is not above the pair")
    if not re.search(r"<li>[^<]*\(<a class=\"xslz-download\"", listed):
        fails.append("the pair broke the list item it was written in")
    usage = open(os.path.join(site, "usage", "index.html")).read()
    if "xslz-note" in usage:
        fails.append("a page with no samples carries the note anyway")

    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "docs_packages_leg", os.path.join(_REPO, "hooks", "docs_packages.py"))
    hook = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(hook)
    scratch_config = {"docs_dir": docs, "site_dir": site}

    # --- a sample that gains a sidecar gains it in its package ----------------
    # Rebuilding leaves an up-to-date package alone (mkdocs serve would otherwise
    # repack the whole corpus on every save), so the property that matters is that a
    # CHANGED project is not left alone: nothing about a package is hand-maintained.
    pkg = os.path.join(site, "lem", "files", single + PACKAGE_EXT)
    before = os.path.getmtime(pkg)
    _quiet(hook.on_post_build, scratch_config)
    if os.path.getmtime(pkg) != before:
        fails.append("an unchanged sample was repacked, so every serve rebuild "
                     "repacks the whole corpus")
    with open(os.path.join(docs, "lem", "files", single + "_mesh.json"), "w") as fh:
        fh.write("{}")
    _quiet(hook.on_post_build, scratch_config)
    if single + "_mesh.json" not in package_contents(pkg):
        fails.append("a sample that gained a sidecar kept its old package")

    # --- the build refuses to ship a link it did not honour -------------------
    hook._linked.add("lem/files/never_built" + PACKAGE_EXT)
    try:
        _quiet(hook.on_post_build, scratch_config)
    except hook.DocsPackageError as exc:
        if "never_built" not in str(exc):
            fails.append(f"the missing-package refusal does not name it: {exc}")
    else:
        fails.append("a page could link a package nobody built and the build passed")
    finally:
        hook._linked.clear()
    return fails


# ------------------------------------------------------------------ F. the flow
def test_studio_flow():
    """Refuse, then ask, then fetch, then the ordinary open — in that order."""
    from PySide6.QtWidgets import QApplication, QDialog, QMessageBox

    from studio.dialogs import UnpackPackageDialog
    from studio.main_window import MainWindow

    QApplication.instance() or QApplication([])
    fails = []

    folder = _tmp()
    for src in xslope.project_files(WITH_SIDECARS):
        shutil.copy2(src, folder)
    book = os.path.join(folder, os.path.basename(WITH_SIDECARS))
    served = xslope.pack(book, dest=os.path.join(_tmp(), "served" + PACKAGE_EXT))

    mw = MainWindow()
    mw._add_recent = lambda path: None          # never touch the user's settings
    downloads = _tmp()
    mw.download_dir = lambda: downloads

    asked, fetched, warned = [], [], []

    def _fake_download(url, dest_dir, **kw):
        fetched.append((url, dest_dir))
        out = os.path.join(dest_dir, urlscheme.package_name(url))
        shutil.copy2(served, out)
        return out

    _question, _warning, _critical = (QMessageBox.question, QMessageBox.warning,
                                      QMessageBox.critical)
    _download = urlscheme.download_package
    QMessageBox.warning = staticmethod(
        lambda *a, **k: (warned.append(a[2]), QMessageBox.Ok)[1])
    QMessageBox.critical = staticmethod(
        lambda *a, **k: (warned.append(a[2]), QMessageBox.Ok)[1])
    urlscheme.download_package = _fake_download

    dest = os.path.join(_tmp(), "from_link")
    _exec = UnpackPackageDialog.exec

    def _auto_exec(self):
        self.dest.setText(dest)
        self._mode = "extract"
        return QDialog.Accepted

    # The destination dialog is stubbed for the WHOLE leg, not only for the phase
    # that expects to reach it: a regression that let a refused link through would
    # otherwise stop at a modal dialog nobody can answer instead of failing.
    UnpackPackageDialog.exec = _auto_exec
    try:
        # 1. A refused link never reaches the network, and says why.
        QMessageBox.question = staticmethod(
            lambda *a, **k: (asked.append(a[2]), QMessageBox.Yes)[1])
        for bad in (f"xslope://open?url=https://evil.example/x{PACKAGE_EXT}",
                    f"xslope://docs?url={GOOD}",
                    "xslope://open?url=file:///etc/passwd"):
            _quiet(mw.open_scheme_url, bad)
        if fetched:
            fails.append(f"a refused link still fetched {fetched}")
        if asked:
            fails.append("a refused link still asked the user to confirm a download")
        if len(warned) != 3:
            fails.append(f"{len(warned)} of 3 refused links were shown to the user")
        if warned and not any("evil.example" in str(w) for w in warned):
            fails.append("the refusal shown to the user does not name the host")

        # 2. The confirmation comes BEFORE the download, and Cancel means nothing
        #    was fetched.
        asked.clear()
        QMessageBox.question = staticmethod(
            lambda *a, **k: (asked.append(a[2]), QMessageBox.Cancel)[1])
        out = _quiet(mw.open_scheme_url, urlscheme.build_url(FLOW))
        if not asked:
            fails.append("a link downloaded without asking")
        if fetched:
            fails.append("cancelling the confirmation still downloaded the package")
        if out is not None:
            fails.append(f"a cancelled link returned {out}")
        if asked:
            text = str(asked[0])
            for wanted in (urlscheme.package_name(FLOW), "xslope.readthedocs.io",
                           downloads):
                if wanted not in text:
                    fails.append(f"the confirmation does not say {wanted!r}: {text}")

        # 3. Confirmed: fetched once, then the ordinary unpack-then-open path.
        asked.clear()
        QMessageBox.question = staticmethod(
            lambda *a, **k: (asked.append(a[2]), QMessageBox.Yes)[1])
        pkg = _quiet(mw.open_scheme_url, urlscheme.build_url(FLOW))
        if len(fetched) != 1:
            fails.append(f"a confirmed link fetched {len(fetched)} times")
        elif fetched[0][1] != downloads:
            fails.append(f"the package was saved in {fetched[0][1]}, not the "
                         f"downloads folder")
        if pkg != os.path.join(downloads, urlscheme.package_name(FLOW)):
            fails.append(f"the link saved the package as {pkg}")
        if not mw.doc.is_open:
            fails.append("a confirmed link opened no project")
        elif os.path.normpath(mw.doc.path) != os.path.normpath(
                os.path.join(dest, os.path.basename(book))):
            fails.append(f"the link left the document on {mw.doc.path}, not the "
                         f"extracted workbook")
        if mw.doc.slope_data.get("mesh") is None:
            fails.append("the project opened from a link carries no mesh, so the "
                         "sidecars did not travel")
    finally:
        QMessageBox.question, QMessageBox.warning = _question, _warning
        QMessageBox.critical = _critical
        urlscheme.download_package = _download
        UnpackPackageDialog.exec = _exec
        mw.doc._dirty = False
        mw.close()
    return fails


CHECKS = [
    ("A. one verb, and it is named", test_verb_gate),
    ("B. only XSLOPE's own sites", test_allowlist),
    ("C. redirects are re-checked", test_redirect_guard),
    ("D. the saved name is a name", test_saved_name),
    ("E. the docs build packages and pairs", test_docs_build),
    ("F. refuse, ask, fetch, open", test_studio_flow),
]


def run():
    """Every check; returns a list of failure strings (empty = pass)."""
    try:
        import PySide6                                    # noqa: F401
        checks = CHECKS
    except Exception:
        print("docs links: PySide6 not installed — Studio check skipped.")
        checks = [c for c in CHECKS if c[1] is not test_studio_flow]
    failures = []
    try:
        for name, fn in checks:
            try:
                fs = fn()
            except Exception as exc:
                import traceback
                traceback.print_exc()
                fs = [f"{name} raised: {exc!r}"]
            print(f"  {name:44s} {'ok' if not fs else f'FAIL ({len(fs)})'}")
            failures += fs
    finally:
        for d in _TMPDIRS:
            shutil.rmtree(d, ignore_errors=True)
        _TMPDIRS.clear()
    return failures


def main():
    print("docs sample links and the xslope:// scheme:")
    failures = run()
    if failures:
        print("\nFAILURES:")
        for f in failures:
            print(f"  - {f}")
        raise SystemExit(1)
    print("\nAll docs-link checks passed.")


if __name__ == "__main__":
    import matplotlib
    matplotlib.use("Agg")
    try:
        from PySide6.QtWidgets import QApplication, QMessageBox
        _app = QApplication.instance() or QApplication([])
        QMessageBox.warning = staticmethod(lambda *a, **k: QMessageBox.Ok)
        QMessageBox.information = staticmethod(lambda *a, **k: QMessageBox.Ok)
        QMessageBox.critical = staticmethod(lambda *a, **k: QMessageBox.Ok)
    except Exception:
        pass
    main()
