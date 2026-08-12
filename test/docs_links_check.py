"""Standing checks on the docs' sample links and the ``xslope://`` scheme behind them.

Two halves of one road: the docs build packages every sample and links it twice, and
Studio answers the link. Both halves have failure modes nobody would see by looking:

  A. THE VERB GATE. One verb is registered, ``open``, and the handler must refuse
     everything else BY NAME rather than guess — a future ``xslope://docs?...`` link
     clicked on an old build has to say so, not do something surprising with it. The
     gate is mutation-tested: widen ``VERBS`` and the same link is accepted, which is
     how we know the refusal comes from the gate and not from a parse failure.
  B. THE ALLOWLIST. The scheme is registered system-wide, so any page anywhere can
     put an ``xslope://`` link in front of a user. One host may be fetched from — the
     documentation site, where every emitted link points — over https, on port 443,
     and never a local path: a ``file:`` URL through this door would let a web page
     hand Studio a file on the user's disk. GitHub is deliberately absent and checked
     to be: those hosts serve every repository anyone has pushed, so allowing them
     allows an attacker's as readily as ours. Mutation-tested: add a hostile host to
     ``ALLOWED_HOSTS`` and it is accepted, so the constant is the whole decision.
  C. REDIRECTS ARE RE-CHECKED. ``urlopen`` follows redirects itself, so a check made
     once on the URL in the link says nothing about where the bytes came from. An
     open redirect on an allowlisted host is the standard way around a bare check.
  D. THE SAVED NAME COMES FROM THE URL'S PATH, and is a plain file name — a package
     downloaded from a link may not name a directory, a parent, a drive, or an
     extension of the server's choosing. Naming also writes the confirmation dialog,
     before anything is fetched, so it may not raise on ANY url: a percent-encoded
     NUL is a name ``open()`` refuses, which is a crash arriving from a web page.
  E. THE DOCS BUILD PACKAGES AND PAIRS. A real MkDocs build over a scratch docs tree
     has to produce a ``.xslz`` per sample workbook, holding its sidecars, and a
     Download · Open in Studio pair per link whose two hrefs name the SAME package —
     the one property the whole "emit both from one URL" design exists for. The
     escape hatch (``.raw-file``) and non-sample links must come through untouched,
     and the build must refuse to finish if a page ever links a package nobody built.
     A sample link the hook does NOT pair — an href that needs decoding, one carrying
     a fragment, one whose file has moved — is warned about rather than passed over
     in silence, which is how a page keeps a bare workbook link nobody notices.
  F. STUDIO'S FLOW, IN ORDER, AND NEVER AN EXCEPTION. Refused links never reach the
     network; the confirmation comes BEFORE the download, and cancelling it downloads
     nothing; and what a link opens is the ordinary unpack-then-open path, so the
     document ends up on the extracted workbook rather than on the package. A link is
     delivered by the operating system, so anything RAISED out of the handler reaches
     Qt and takes the window down with the user's unsaved work in it: malformed URLs
     have to produce dialogs, not tracebacks.
  G. BOTH ARRIVALS END IN ONE CALL. Windows and Linux deliver a link and a
     double-clicked file as argv; macOS delivers both as an event, to a copy of
     Studio that may already be running — or, when the click is what launched it, to
     an application that has no window yet, which must hold the request rather than
     drop it.

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

    # A link malformed enough to break the parser is refused in the module's own
    # currency. Anything else — main_window catches SchemeError, and what it does not
    # catch reaches Qt from an OS-delivered link and closes the window.
    for bad in ("xslope://[evil/?url=" + GOOD,
                "xslope://open:notaport?url=" + GOOD,
                "xslope://open?url=%%%"):
        try:
            urlscheme.parse_request(bad)
        except urlscheme.SchemeError:
            pass
        except Exception as exc:
            fails.append(f"parse_request({bad!r}) raised {type(exc).__name__} rather "
                         f"than a SchemeError: {exc}")

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
        # GitHub is multi-tenant: these hosts serve every repository anyone has ever
        # pushed, so a host name there is not an identity and none of them may be
        # allowlisted. No emitted link points at them.
        "https://github.com/njones61/xslope/raw/main/x.xslz",
        "https://github.com/attacker/payload/raw/main/x.xslz",
        "https://raw.githubusercontent.com/attacker/payload/main/x.xslz",
        "https://objects.githubusercontent.com/x.xslz",
        # A port is part of the address: another port on the same host is another
        # service (a readthedocs.io user's tunnel, a proxy, anything).
        "https://xslope.readthedocs.io:8080/x.xslz",
        "https://xslope.readthedocs.io:80/x.xslz",
    ]
    hostile += [
        # Not URLs at all. urlparse raises on these, and a raise is not a refusal:
        # main_window catches SchemeError, and what it does not catch reaches Qt.
        "https://[/x.xslz",
        "https://xslope.readthedocs.io:notaport/x.xslz",
        "https://[not:an:ipv6]/x.xslz",
    ]
    for url in hostile:
        try:
            urlscheme.check_url(url)
        except urlscheme.SchemeError as exc:
            if url not in str(exc):
                fails.append(f"the refusal of {url} does not name it: {exc}")
        except Exception as exc:
            fails.append(f"check_url({url!r}) raised {type(exc).__name__} rather than "
                         f"a SchemeError: {exc}")
        else:
            fails.append(f"{url} was accepted for download")

    # The port the docs are actually served on is not refused for being written out.
    try:
        urlscheme.check_url("https://xslope.readthedocs.io:443/x.xslz")
    except urlscheme.SchemeError as exc:
        fails.append(f"port 443 written out was refused: {exc}")

    # MUTATION: the list is the whole decision. Add a hostile host and it passes —
    # once for an invented host, once for GitHub, which is what was removed.
    saved = urlscheme.ALLOWED_HOSTS
    for host, url in [("evil.example", "https://evil.example/x.xslz"),
                      ("github.com",
                       "https://github.com/attacker/payload/raw/main/x.xslz")]:
        urlscheme.ALLOWED_HOSTS = saved + (host,)
        try:
            urlscheme.check_url(url)
        except urlscheme.SchemeError as exc:
            fails.append(f"with {host} allowlisted, {url} is still refused ({exc}), "
                         f"so ALLOWED_HOSTS is not what refuses it")
        finally:
            urlscheme.ALLOWED_HOSTS = saved

    if urlscheme.ALLOWED_HOSTS != ("xslope.readthedocs.io",):
        fails.append(f"this build fetches from {urlscheme.ALLOWED_HOSTS}; only the "
                     f"docs site publishes packages, and every emitted link is on it")
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

    # Naming is also done to WRITE THE DIALOG, before anything is fetched, so it may
    # not raise on any URL at all. A percent-encoded NUL is the one that bites: open()
    # rejects it with a ValueError of its own, after the download.
    for url in ("https://xslope.readthedocs.io/a/%00bad.xslz",
                "https://xslope.readthedocs.io/a/b%2F%2E%2E%2Fc.xslz",
                "https://xslope.readthedocs.io/a/D%3Apwned.xslz",
                "https://[/x.xslz", "https://x.io:notaport/x.xslz", "", "::::"):
        try:
            got = urlscheme.package_name(url)
        except Exception as exc:
            fails.append(f"package_name({url!r}) raised {type(exc).__name__}: {exc}")
            continue
        if "\x00" in got or os.sep in got or ":" in got:
            fails.append(f"{url!r} would be saved as {got!r}, which is not a plain name")
        try:                          # the name has to be one the filesystem takes
            open(os.path.join(_tmp(), got), "wb").close()
        except Exception as exc:
            fails.append(f"{url!r} yields {got!r}, which cannot be written: {exc}")

    # Cleaning up after a failed download must not raise on its way out of one:
    # whatever went wrong is what the user has to be told about.
    try:
        urlscheme._unlink("\x00")
        urlscheme._unlink(None)
    except Exception as exc:
        fails.append(f"cleanup raised {type(exc).__name__} ({exc}), replacing the "
                     f"refusal it was cleaning up after")
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
    # A page whose first pair is NOT in a paragraph: there is nowhere above it to put
    # a note paragraph without breaking the list, so it goes below the list.
    with open(os.path.join(docs, "lem", "listed.md"), "w") as fh:
        fh.write(f"# Listed\n\n* [{single}.xlsx](files/{single}.xlsx)\n"
                 f"* nothing else\n")
    # The tutorial shape, in full: title, the problem sketch the reader chooses the
    # page by, the objectives, then the header table whose Completed-model row is the
    # page's first pair. This is the page that pins where the fallback puts the note.
    with open(os.path.join(docs, "lem", "tabled.md"), "w") as fh:
        fh.write(f"# Tabled\n\n"
                 f"![problem sketch](sketch.png)\n\n"
                 f"**Objectives** — build the thing.\n\n"
                 f"| | |\n|---|---|\n"
                 f"| **Analysis** | Limit equilibrium |\n"
                 f"| **Completed model** | [{single}.xlsx](files/{single}.xlsx) |\n\n"
                 f"Prose, after the table.\n")
    open(os.path.join(docs, "lem", "sketch.png"), "wb").close()
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


def _content(page_html):
    """The page's own content, without the theme's chrome around it.

    Placement checks are about reading order, and the navigation sidebar renders a
    heading and half a dozen closed lists before the content starts. Measured
    against the whole file, "the note comes after a </ul>" is true of every page
    ever built and proves nothing.
    """
    at = page_html.find('role="main"')
    return page_html[at:] if at != -1 else page_html


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
    listed = _content(open(os.path.join(site, "lem", "listed", "index.html")).read())
    if listed.count("xslz-note") != 1:
        fails.append("a page whose first pair is inside a list got no note")
    elif not 0 <= listed.find("</ul>") < listed.find("xslz-note"):
        fails.append("the note on the list page is not below the list holding the "
                     "pair, so it was written inside the list or above it")
    if not re.search(r"<li>[^<]*\(<a class=\"xslz-download\"", listed):
        fails.append("the pair broke the list item it was written in")

    # WHERE the fallback puts it, pinned as the whole reading order. A tutorial's
    # only pair is the Completed-model row of the header table, and a reader decides
    # whether this is the page they want by looking at the sketch — so the note goes
    # BELOW the header, next to the links it explains, and nothing this hook adds
    # comes between the title and the picture.
    tabled = _content(open(os.path.join(site, "lem", "tabled", "index.html")).read())
    if tabled.count("xslz-note") != 1:
        fails.append(f"a page whose first pair is inside a table carries "
                     f"{tabled.count('xslz-note')} notes, not exactly one")
    else:
        order = [("the page's <h1>", tabled.find("</h1>")),
                 ("the problem sketch", tabled.find("sketch.png")),
                 ("the header table", tabled.find("<table>")),
                 ("the first pair", tabled.find("xslz-download")),
                 ("the end of the table", tabled.find("</table>")),
                 ("the note", tabled.find("xslz-note"))]
        missing = [name for name, at in order if at < 0]
        if missing:
            fails.append(f"the table page rendered no {', '.join(missing)}")
        else:
            offsets = [at for _, at in order]
            if offsets != sorted(offsets):
                fails.append(
                    "the table page does not read title -> sketch -> table -> pair "
                    "-> note; it reads " +
                    " -> ".join(name for name, _ in sorted(order, key=lambda p: p[1])))
    if not re.search(r"<td>[^<]*\(<a class=\"xslz-download\"", tabled):
        fails.append("the pair broke the table cell it was written in")
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

    # --- a link that is NOT turned into a pair is said out loud ---------------
    # Silence is the expensive answer here: the page keeps a bare .xlsx link and
    # nobody notices it stopped being a pair.
    # Hrefs as MkDocs writes them for a page at lem/samples/ — one that needs
    # decoding (and whose file really is there), one with a fragment, one whose file
    # has moved, and one pointing at the package by hand.
    page = (f'<h1>x</h1><p><a href="../files/{single}%2Exlsx">encoded</a> '
            f'<a href="../files/{single}.xlsx#sheet2">fragment</a> '
            f'<a href="../files/gone.xlsx">missing</a> '
            f'<a href="../files/{single}{PACKAGE_EXT}">by hand</a></p>')
    _html, linked, warnings = hook.rewrite_links(
        page, "lem/samples/", docs, "https://xslope.readthedocs.io/en/latest/")
    for wanted in (f"{single}%2Exlsx", "#sheet2", "gone.xlsx"):
        if not any(wanted in w for w in warnings):
            fails.append(f"a link that was left alone ({wanted}) was skipped in "
                         f"silence: {warnings}")
    if f"lem/files/{single}{PACKAGE_EXT}" not in linked:
        fails.append("a page linking a package by hand does not register it, so the "
                     "build never checks that the package exists")
    if "xslz-download" in _html:
        fails.append("a link that could not be resolved was paired anyway")

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
        refused = [
            f"xslope://open?url=https://evil.example/x{PACKAGE_EXT}",
            f"xslope://docs?url={GOOD}",
            "xslope://open?url=file:///etc/passwd",
            f"xslope://open?url=https://github.com/attacker/payload/x{PACKAGE_EXT}",
            f"xslope://open?url=https://xslope.readthedocs.io:8080/x{PACKAGE_EXT}",
            # Malformed enough to raise out of urlparse rather than fail a test. A
            # link arrives from the OS: an exception here reaches Qt and takes the
            # window down, with whatever the user had open in it.
            "xslope://open?url=https://[/x.xslz",
            "xslope://[evil/?url=https://xslope.readthedocs.io/x.xslz",
            "xslope://open?url=https://xslope.readthedocs.io:notaport/x.xslz",
        ]
        for bad in refused:
            try:
                _quiet(mw.open_scheme_url, bad)
            except Exception as exc:
                fails.append(f"{bad} raised {type(exc).__name__} out of the handler "
                             f"({exc}) — from a link, that is Studio closing")
        if fetched:
            fails.append(f"a refused link still fetched {fetched}")
        if asked:
            fails.append("a refused link still asked the user to confirm a download")
        if len(warned) != len(refused):
            fails.append(f"{len(warned)} of {len(refused)} refused links were shown "
                         f"to the user")
        for wanted in ("evil.example", "github.com", "8080"):
            if not any(wanted in str(w) for w in warned):
                fails.append(f"no refusal shown to the user names {wanted}")

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


# --------------------------------------------------------------- G. arrival
def test_arrival():
    """Both ways a link reaches Studio end in the same call.

    Windows and Linux deliver a clicked ``xslope://`` link and a double-clicked file
    the same way — as ``argv[1]`` — and macOS delivers both as a ``QEvent.FileOpen``,
    to a copy of Studio that may have been running for hours. A link that arrives
    before the window exists (the click that launched the app) has to be held rather
    than dropped.
    """
    from PySide6.QtCore import QEvent
    from PySide6.QtWidgets import QApplication

    from studio import app as studio_app

    QApplication.instance() or QApplication([])
    fails = []

    class _Window:
        def __init__(self):
            self.links, self.paths = [], []

        def open_scheme_url(self, uri):
            self.links.append(uri)

        def open_path(self, path):
            self.paths.append(path)

    link = urlscheme.build_url(GOOD)
    here = os.path.abspath(__file__)

    win = _Window()
    studio_app.open_request(win, link)
    studio_app.open_request(win, here)
    studio_app.open_request(win, os.path.join(_tmp(), "no_such_file.xlsx"))
    if win.links != [link]:
        fails.append(f"an xslope:// argument reached {win.links}, not the handler")
    if win.paths != [here]:
        fails.append(f"a file argument reached {win.paths}, not the normal open")

    # The macOS event, including one that arrives before there is a window.
    class _FakeApp:                       # the real methods, off a real QApplication
        event = studio_app.StudioApplication.event
        set_window = studio_app.StudioApplication.set_window
        _deliver = studio_app.StudioApplication._deliver

        def __init__(self):
            self._window, self._pending = None, []

    class _Url:
        def __init__(self, text):
            self.text = text

        def scheme(self):
            return self.text.split(":")[0]

        def toString(self):
            return self.text

    class _FileOpen:
        def __init__(self, text, file=""):
            self.url_, self.file_ = _Url(text), file

        def type(self):
            return QEvent.FileOpen

        def url(self):
            return self.url_

        def file(self):
            return self.file_

    early = _FakeApp()
    if early.event(_FileOpen(link)) is not True:
        fails.append("a FileOpen event was not accepted")
    win2 = _Window()
    if win2.links or win2.paths:
        fails.append("the stand-in window started dirty")
    early.set_window(win2)
    if win2.links != [link]:
        fails.append(f"a link that arrived before the window was dropped: {win2.links}")

    running = _FakeApp()
    running.set_window(win2)
    running.event(_FileOpen("file://" + here, file=here))
    if win2.paths != [here]:
        fails.append(f"a double-clicked file reached {win2.paths}")
    return fails


CHECKS = [
    ("A. one verb, and it is named", test_verb_gate),
    ("B. only XSLOPE's own sites", test_allowlist),
    ("C. redirects are re-checked", test_redirect_guard),
    ("D. the saved name is a name", test_saved_name),
    ("E. the docs build packages and pairs", test_docs_build),
    ("F. refuse, ask, fetch, open", test_studio_flow),
    ("G. argv and FileOpen end in one call", test_arrival),
]


def run():
    """Every check; returns a list of failure strings (empty = pass)."""
    try:
        import PySide6                                    # noqa: F401
        checks = CHECKS
    except Exception:
        print("docs links: PySide6 not installed — Studio check skipped.")
        checks = [c for c in CHECKS
                  if c[1] not in (test_studio_flow, test_arrival)]
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
