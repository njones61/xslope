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
  H. STUDIO READS THE WAY THE TUTORIALS SAY IT DOES. A tutorial hands the reader a
     label to press, and a rewording in the app turns that sentence into an
     instruction to press something that is not there. The labels Tutorial 0 quotes
     for the file lifecycle and the package flow are pinned here, the way
     ``circles_editor_check.GENERATE_LABEL`` pins LEM-1's generator button.

No network: the download is stubbed in F and the redirect handler is exercised
directly in C. The MkDocs build in E writes only into a temporary directory.
"""
import contextlib
import glob
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
#: Tutorial LEM-2's completed model, which its label pins in H are read against —
#: the distributed-loads editor has a Direction to report only on a model that
#: carries a load.
LEM02_FILE = os.path.join(_REPO, "docs/lem/files/xslope_crest_surcharge.xlsx")
#: Tutorial LEM-1's completed model — one material, which is the state its Studio
#: path leaves the materials editor in and the state its figure shows.
LEM01_FILE = os.path.join(_REPO, "docs/lem/files/xslope_simple_embankment.xlsx")
#: Tutorial LEM-3's completed model — two materials and two profile lines, which is
#: what its label pins need: the materials table read as a pair of rows, and a
#: profile editor with a material to assign each line to.
LEM03_FILE = os.path.join(_REPO, "docs/lem/files/xslope_simple_mult_layers.xlsx")
#: Tutorial LEM-5's completed model — four material zones with a soft seam among
#: them, which is what its pins need: a non-circular surface to read the vertex
#: table against, and a weak zone for the generator to find.
LEM05_FILE = os.path.join(_REPO, "docs/lem/files/xslope_noncircular.xlsx")
#: Tutorial LEM-4's completed model — three piezo-reading materials, an eight-point
#: piezometric line, and one specified circle: what its pins need is a model whose
#: Run LEM Surface row is the fixed *Circular* label, with the γ_sat column filled
#: as the page's weight-split section reads it.
LEM04_FILE = os.path.join(_REPO,
                          "docs/lem/files/xslope_method_slices_problem.xlsx")
#: Tutorial LEM-6's completed model — two material-zone polygons on a dipping
#: base, which is what its pins need: a polygon-based project (the only kind
#: whose Polygons row opens an editor) and a domain whose floor a circle can be
#: pushed below.
LEM06_FILE = os.path.join(_REPO, "docs/lem/files/xslope_sloping_bottom.xlsx")
#: Tutorial LEM-8's completed model — six reinforcement lines, which is what its
#: pins need: an Inputs tree row that counts them and an editor with a line to
#: read the capacity fields off.
LEM08_FILE = os.path.join(_REPO, "docs/lem/files/xslope_reinforce.xlsx")
#: Tutorial LEM-9's completed model — two Anchor-type reinforcement lines and one
#: soldier pile, which is what its pins need: a reinforcement row whose Type preset
#: has filled Dir/Appl, and a Piles row in the Inputs tree with a pile to read the
#: editor's fields off.
LEM09_FILE = os.path.join(_REPO,
                          "docs/verification/files/rocscience/vp049.xlsx")
#: Tutorial LEM-10's completed model — a cohesionless embankment on soft clay,
#: the section its generator quote and its two Run LEM search options are read on.
LEM10_FILE = os.path.join(_REPO, "docs/lem/files/xslope_mult_min_KEY.xlsx")
#: Tutorial LEM-11's completed model — the submerged clay slope, the only tutorial
#: section that carries reliability standard deviations, which is what gates the
#: dialog controls and the plot type its pins are read on.
LEM11_FILE = os.path.join(_REPO,
                          "docs/lem/files/xslope_prob_submerged_KEY.xlsx")
#: Tutorial LEM-12's completed model — the pile-stabilized clay slope, the only
#: tutorial model whose pile force is left to be computed, which is the state its
#: pins are read in.
LEM12_FILE = os.path.join(_REPO, "docs/lem/files/xslope_piles.xlsx")
#: Tutorial SEEP-1's completed model — the confined sheetpile problem, whose two
#: specified-head boundaries and absent exit face are what its pins are read on.
SEEP01_FILE = os.path.join(_REPO, "docs/seep/files/xslope_clay_blanket.xlsx")
#: The editable master template, whose ``reinforce`` sheet carries the support-type
#: lookup block LEM-8 reproduces as a table.
TEMPLATE_FILE = os.path.join(_REPO, "docs/inputs/input_template.xlsx")

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


# ------------------------------------------ H. the labels the tutorials quote
#: The File-menu actions Tutorial 0 tells the reader to use. These are the SOURCE
#: strings with Qt's ``&`` accelerator markers removed — what the menu renders — not
#: the page's spelling of them, which drops the trailing ellipsis on the familiar
#: File verbs. What is guarded is that the app still calls these things what the
#: tutorial says it calls them; a rewording there turns those sentences into
#: instructions to press something that is not there, the same failure
#: ``circles_editor_check.GENERATE_LABEL`` guards for LEM-1, and nothing else in the
#: suite would notice it.
T0_FILE_ACTIONS = {
    "act_new": "New",
    "act_open": "Open…",
    "act_save": "Save",
    "act_save_as": "Save As…",
    "act_export_pkg": "Export Project Package…",
}

#: The buttons the package-open dialog offers when its destination already exists —
#: the choice Tutorial 0 walks the reader through.
T0_PACKAGE_BUTTONS = ("Change…", "Open Existing", "Extract Fresh")

#: The Inputs-tree rows Tutorial 0 names as the way into an editor.
T0_INPUT_CATEGORIES = ("Global parameters", "Materials", "Profile lines", "Circles")

#: The Inputs-tree rows LEM-2 sends the reader to. Both are added unconditionally,
#: so an empty model still lists them — which is the state the tutorial's reader is
#: in when they go looking.
LEM02_INPUT_CATEGORIES = ("Distributed loads", "Line loads")

#: The Run-menu actions LEM-2 tells the reader to press — the run itself, and the
#: sweep its design section ends on. Same source-string convention as
#: ``T0_FILE_ACTIONS``: Qt's ``&`` removed, the trailing ellipsis kept (the page
#: writes **Run LEM…** and **Parametric…** too, since neither is one of the
#: familiar File verbs that drop it).
LEM02_RUN_ACTIONS = {"act_run": "Run LEM…", "act_sensitivity": "Parametric…"}

#: The two buttons LEM-2's Studio path presses to get a load into the model, and
#: the Direction reading it leaves selected. The direction wording is the pinned
#: one because the page's face-load section turns on the difference between the two
#: options — a reworded combo would leave that section describing a choice the
#: reader cannot find.
LEM02_DLOAD_BUTTONS = ("Add load", "Add row")
LEM02_DLOAD_DIRECTION = "Normal (perpendicular to the line)"

#: The two tabs the page names in the distributed-loads editor, and the label on
#: the Direction chooser. The page tells the reader which tab the ordinary load
#: goes on and which one is the rapid-drawdown stage that stays empty, so a
#: renamed tab leaves those sentences pointing at nothing.
LEM02_DLOAD_TABS = ("Set 1", "Set 2 (rapid drawdown)")
LEM02_DLOAD_DIRECTION_LABEL = "Direction:"

#: The Global-parameters row LEM-2's seismic section sends the reader to.
LEM02_GLOBAL_ROW = "Seismic coefficient k"

#: The line-loads table's column headers, in the order the page dictates values
#: into them (**Label** `footing`, **x** `30`, **y** `20`, **P** `7500`,
#: **Angle** `-90`).
LEM02_LLOAD_HEADERS = ("Label", "x", "y", "P", "Angle")

#: The Parametric dialog in the state LEM-2 walks: the Design mode entry, every
#: form label the page names in order, and the button that starts the sweep.
LEM02_DESIGN_MODE = "Design (FS target)"
LEM02_DESIGN_ROWS = ("Mode", "Method", "Number of slices", "Material", "Property",
                     "Sweeping", "From", "To", "Steps", "Target FS")
LEM02_DESIGN_RUN = "Run"

#: The checkbox the design step tells the reader to leave ticked, and whose
#: wording the step's explanation of when re-searching matters is built on.
LEM02_DESIGN_SEARCH = "Re-search the critical surface at each step"


#: The materials-editor controls LEM-1's Studio path uses, and the view it opens on.
#: The page sends the reader from the view the editor opens on to the one its figure
#: shows, so the order is the guarded part: **Table view** is where a first open
#: lands (``editors._LAST_MATERIALS_VIEW`` starts there), **List view** is the
#: per-material form the page's numbered fields belong to, and **Add** is that
#: form's button — the table's reads "Add row", so a step that pressed the wrong one
#: would be pressing a button that is not in the view it just asked for.
#: The two Global-parameters rows LEM-1's first step names by label — the one the
#: reader sets, and the one the step tells them to leave empty. Both are read off
#: the same form the page's figure photographs, so a renamed row is a step pointing
#: at a control the picture no longer shows.
LEM01_GLOBAL_ROWS = ("Units", "Time")

LEM01_MATERIALS_OPENS_ON = "table"
LEM01_MATERIALS_VIEWS = ("Table view", "List view")
LEM01_MATERIALS_ADD = "Add"

#: The materials-editor controls LEM-3's Studio path uses. It builds two materials
#: in the view LEM-1 does not use — the table, where the row order is what fixes the
#: Mat IDs the profile lines reference — so both the view toggle and the button that
#: adds a row are labels the page tells the reader to press.
LEM03_MATERIALS_BUTTONS = ("Table view", "Add row")

#: The profile-lines editor's controls, in the order LEM-3 drives them: add a line,
#: give it a material, add its vertices — plus the bottom-boundary field, whose
#: wording carries the fact that the value is an elevation.
LEM03_PROFILE_BUTTONS = ("Add line", "Add row")
LEM03_PROFILE_LABELS = ("Material:", "Max depth (bottom boundary elevation):")

#: The material choices the page dictates, as the combo spells them. A profile line
#: names its material by ID and name together, which is the notation the page's
#: numbered steps quote.
LEM03_PROFILE_MATERIALS = ("1: embankment", "2: foundation")

#: What the starting-circle generator reports on LEM-3's two-layer section, quoted
#: verbatim on the page beside a capture of the same line. The count and the phrasing
#: are both load-bearing: the page's audit table walks the three circles it announces
#: (one through the toe, one at the base of each layer) row by row.
LEM03_GENERATE_SUMMARY = "3 on the left-facing face (toe at x = 0, height 20)"
LEM03_GENERATE_DEPTHS = (-4.72136, -10.0, 0.0)

#: The Inputs-tree rows LEM-5 sends the reader to. Both are added unconditionally,
#: so an empty model still lists them — the state its reader is in when they go
#: looking for the failure surface.
LEM05_INPUT_CATEGORIES = ("Non-circular pts", "Piezometric lines")

#: The non-circular editor's controls, in the order LEM-5 drives them: add the four
#: vertices, then the generator that builds a surface from the weak zone instead.
LEM05_NONCIRC_BUTTONS = ("Add row", "Generate from the weak zone…")

#: The vertex table's columns and the Movement settings the page dictates into
#: them. Movement is pinned option by option because the page's search reading
#: measures what each one does to the search -- `Horiz` sliding the seam points
#: along the clay, `Free` walking the ends along the ground, `Fixed` pinning them
#: -- so a renamed option leaves that reading describing choices the reader cannot
#: make.
LEM05_NONCIRC_HEADERS = ("X", "Y", "Movement")
LEM05_NONCIRC_MOVEMENTS = ("Free", "Horiz", "Fixed")

#: What the weak-zone generator reports on LEM-5's section, quoted verbatim on the
#: page beside a capture of the same line — the zone it seeds on, the two strengths
#: it compared to get there, and the ramp angles it built. The point count is
#: pinned with it because the summary opens by announcing it. Four on this section:
#: the seam is flat, so the track between the two end ramps is one straight
#: horizontal segment and a vertex partway along it would be one the surface does
#: not bend at and the search cannot move.
LEM05_GENERATE_SUMMARY = (
    "seeding on 'Soft Clay' -- mobilizable strength 200 against 570 for the next "
    "weakest ('Sand Fill')")
LEM05_GENERATE_RAMPS = ("a 28 degree ramp to the ground at the toe",
                        "a 60 degree ramp to the ground at the crest")
LEM05_GENERATE_POINTS = 4

#: The Run LEM readings LEM-5's run step quotes. BOTH analyses are pinned because
#: the page names both by label: it runs *Auto search* first, as the normal way to
#: run a model, and switches to *Single surface* for the comparisons it announces
#: as holding the surface still. With them, the fixed Surface label a model
#: carrying a non-circular surface and no circles produces (the mirror of the fixed
#: "Circular" LEM-2's model shows).
LEM05_RUN_ANALYSES = ("Auto search", "Single surface")
LEM05_RUN_SURFACE = "Non-circular"

#: What a reader sees when a non-circular search is seeded with an end ramp
#: steeper than the search's own ``max_base_angle``: the Run box's message (the
#: ``AnalysisError`` Studio's LEM runner turns into a "LEM run failed" box) and the
#: Log pane's line under it. The page quotes both verbatim in its end-ramp section,
#: demonstrated on the crest point pulled in to x = 35 -- a 72.4 degree exit ramp,
#: past the 65 degree limit -- so the seed is pinned with them. The starting
#: surface is rejected before any iteration runs, which is why this pin costs a
#: rejected slice-generation rather than a search.
LEM05_STEEP_EXIT_X = 35.0
LEM05_SEARCH_FAILED = "Search produced no valid surfaces."
LEM05_SEARCH_LOG = ("the starting surface is not viable (slice generation or the "
                    "solver failed on it)")

#: The Run LEM readings LEM-4's pivot step quotes: the Single-surface analysis it
#: chooses on a circles-only model, whose Surface row is the fixed *Circular* label
#: the page describes.
LEM04_RUN_ANALYSIS = "Single surface"
LEM04_RUN_SURFACE = "Circular"

#: The piezometric-lines editor as LEM-4's Studio path drives it: the tab the
#: water table goes on, the rapid-drawdown tab the page says stays empty, and the
#: button pressed eight times to enter the points.
LEM04_PIEZO_TABS = ("Line 1", "Line 2 (rapid drawdown)")
LEM04_PIEZO_ADD = "Add row"

#: The circles-editor columns LEM-4's pin step types the found circle into, and
#: the dock the page tells the reader to read the search's own statement of that
#: circle out of.
LEM04_CIRCLE_COLUMNS = ("Xo", "Yo", "Option", "Depth")
LEM04_LOG_DOCK = "Log"

#: The polygons editor as LEM-6's Studio step drives it: the two buttons that
#: build a zone, and the three fields the step names beside them. **Type:** is
#: pinned because the step tells the reader to leave it where it is, which is an
#: instruction about a control that has to exist.
LEM06_POLYGON_BUTTONS = ("Add polygon", "Add row")
LEM06_POLYGON_LABELS = ("Type:", "Material:", "Size:")

#: The zone list's own wording, quoted on the page as the item the reader clicks,
#: and the closed-region rule quoted from the editor's help line. Both are the
#: dialog's text rather than the page's paraphrase of it.
LEM06_POLYGON_ITEM = "Polygon 2  (mat 2: foundation)"
LEM06_POLYGON_HELP = ("Each polygon is a closed region (the ring closes "
                      "automatically, so list each vertex once)")

#: The profile-lines field that must be ABSENT from the polygons editor. The page
#: makes a teaching point of it — a polygon model's bottom boundary is drawn, not
#: typed — so a Max depth field appearing here would leave that paragraph wrong.
LEM06_POLYGON_NO_MAX_DEPTH = "Max depth (bottom boundary elevation):"

#: The Run LEM checkbox LEM-6's composite section tells the reader to tick, and
#: what a reader sees when a circle is pushed below the domain floor without it:
#: the message the ``AnalysisError`` Studio's LEM runner turns into a "LEM run
#: failed" box. The page quotes both verbatim, demonstrated on the file's deeper
#: starting circle taken 1.2 ft further down, so that circle is pinned with them.
LEM06_COMPOSITE_CHECKBOX = "Composite surfaces (truncate circles at the base)"
LEM06_DEEP_CIRCLE = {"Xo": 20.0, "Yo": 40.0, "Depth": -12.0, "R": 52.0}
LEM06_DOMAIN_REFUSAL = "Failure surface extends outside the domain polygon"

#: The reinforcement editor as LEM-8's Studio step drives it: the two buttons
#: that build the list of lines, the view the page sends the reader to for a
#: tiered wall, and the five form groups the page walks in order. The group
#: titles are pinned because the page teaches a line *as* those groups —
#: identity, geometry, capacity, anchorage, type.
LEM08_REINF_BUTTONS = ("Add", "Remove")
LEM08_REINF_VIEWS = ("Table view", "List view")
LEM08_REINF_GROUPS = ("Identity", "Geometry", "Capacity", "Anchorage", "Type")

#: The capacity fields the page names, as label *prefixes*: the list view appends
#: the units the model implies ("Tmax (per unit width, lb/ft)"), which the page
#: quotes for Tmax and which must therefore stay attached to it.
LEM08_REINF_FIELDS = ("Tmax", "Lp1", "Lp2", "Tend1", "Tend2", "Spacing")
LEM08_REINF_TMAX_LABEL = "Tmax (per unit width, lb/ft)"

#: The support-type preset table the page reproduces, as the template's own
#: lookup block (reinforce!Z8:AB11) — the values the Type drop-down offers and
#: the Dir/Appl its formula fills for each. The page's table is this table.
LEM08_TYPE_PRESETS = (("Geosynthetic", "Tangent", "Active"),
                      ("Nail", "Axial", "Passive"),
                      ("Tieback", "Axial", "Active"),
                      ("Anchor", "Axial", "Active"))

#: The Inputs tree row Tutorial LEM-9's Studio path sends the reader to for the
#: soldier pile, and the piles editor as that step drives it: the four form groups
#: the page names in order, and the three capacity fields this problem fills. The
#: group titles are pinned because the page walks the form by them.
LEM09_INPUT_CATEGORY = "Piles"
LEM09_PILE_GROUPS = ("Identity", "Geometry", "Capacity / design", "Behavior")
LEM09_PILE_FIELDS = ("H", "D", "S")
#: What LEM-9's page asserts the Anchor preset does: the Type column reads
#: ``anchor`` on both lines, and Dir/Appl come back filled without either being
#: typed. Read through the loader, which is where the preset is applied.
LEM09_ANCHOR_PRESET = ("anchor", "axial", "active")


#: The two Run LEM search options LEM-10 tells the reader to use, and the field
#: that carries the depth. The page's last section is written as instructions to
#: tick these — a rewording would leave it naming controls that are not there,
#: and the numbers beside them (1.327 seeded from the grid, 1.376 with the filter)
#: are measured with exactly these two settings.
LEM10_GRID_CHECKBOX = "Grid search (auto-seed the circular search)"
LEM10_SKIN_CHECKBOX = "Ignore surficial (skin) failures"
LEM10_MIN_SLIP_LABEL = "Min slip depth"

#: What the starting-circle generator reports on LEM-10's section, quoted verbatim
#: on the page. The skim clause is the load-bearing half: the page takes the
#: generated circle at the base of the embankment as its shallow seed, and the
#: summary is how a reader knows the generator saw a cohesionless face at all.
LEM10_GENERATE_SUMMARY = ("4 on the left-facing face (toe at x = 0, height 15), "
                          "one of them skimming its 24 degree cohesionless face")


def _lem10_run_labels():
    """The Run LEM search options Tutorial LEM-10's last section drives.

    Read on LEM-10's own model, because the dialog enables both options only for
    a circular auto-search — which is the run the page is describing — and the
    generator summary the page quotes is a property of this section's geometry.
    """
    from PySide6.QtWidgets import QCheckBox, QLabel

    from xslope.fileio import load_slope_data
    from xslope.generators import generate_starting_circles

    from studio.dialogs import RunLemDialog

    fails = []
    data = _quiet(load_slope_data, LEM10_FILE)

    run = RunLemDialog(defaults={}, slope_data=data)
    boxes = {b.text() for b in run.findChildren(QCheckBox)}
    for label in (LEM10_GRID_CHECKBOX, LEM10_SKIN_CHECKBOX):
        if label not in boxes:
            fails.append(f"Run LEM has no {label!r} checkbox; Tutorial LEM-10 "
                         f"tells the reader to tick it. Its checkboxes read "
                         f"{sorted(boxes)}")
    labels = {lab.text() for lab in run.findChildren(QLabel)}
    if LEM10_MIN_SLIP_LABEL not in labels:
        fails.append(f"Run LEM has no {LEM10_MIN_SLIP_LABEL!r} field; Tutorial "
                     f"LEM-10 names it as where the depth goes. Its labels read "
                     f"{sorted(l for l in labels if l)}")
    run.deleteLater()

    summary = _quiet(generate_starting_circles, data, report=True).get("summary")
    if summary != LEM10_GENERATE_SUMMARY:
        fails.append(f"the circle generator reports {summary!r} on LEM-10's "
                     f"section, not {LEM10_GENERATE_SUMMARY!r} — the line the "
                     f"page quotes")
    return fails


#: The Run-menu action LEM-11 sends the reader to, in the same source-string
#: convention as ``LEM02_RUN_ACTIONS`` (Qt's ``&`` removed, the ellipsis kept).
LEM11_RUN_ACTIONS = {"act_reliability": "Reliability…"}

#: The Reliability dialog as LEM-11 walks it: the two engines by the names the page
#: writes in its two section headings, the search toggle the page tells the reader
#: to leave ticked, the σ summary box it reads the run's inputs out of, and the
#: three Monte Carlo controls it names one by one. A rewording of any of them turns
#: those steps into instructions to press something that is not there.
LEM11_ENGINES = ("Taylor series (TSPM)", "Monte Carlo")
LEM11_SEARCH_CHECKBOX_TAYLOR = ("Search for the critical surface at each solve "
                                "(MLV and MLV ± σ)")
LEM11_SEARCH_CHECKBOX_MC = "Search for the critical surface at the mean values"
LEM11_SIGMA_GROUP = "Standard deviations in this file"
LEM11_MC_ROWS = ("MC samples", "MC seed", "MC distribution")

#: The materials-editor usage toggle whose default is OFF and whose ticking is the
#: page's first step — the σ boxes do not exist in the editor until it is on — and
#: the Parametric plot type the variance section selects, which the dialog offers
#: only for a model carrying standard deviations.
LEM11_RELIABILITY_TOGGLE = "Reliability"
LEM11_VARIANCE_PLOT = "Variance Pareto (σ)"

#: Tutorial LEM-12 sends the reader to the Analysis Report to read a pile force
#: that is computed rather than entered — the force appears nowhere else a reader
#: looks. Its File-menu action in the same source-string convention as
#: ``LEM02_RUN_ACTIONS``, the word the Piles table prints where the input would be,
#: and the slice-table column that carries the force with the legend sentence the
#: page quotes.
LEM12_REPORT_ACTION = ("act_report", "Generate Report…")
LEM12_COMPUTED_CELL = "computed"
LEM12_HP_HEADER = "H_p (lb/ft)"
LEM12_HP_LEGEND = ("Pile resistance mobilized at the slice base, per unit "
                   "thickness.")
#: The two capacity fields LEM-12 quotes in full, unit suffix included: the suffix
#: is what says a capacity belongs to one shaft rather than to a foot of slope,
#: which is the distinction the page's capacity section is built on.
LEM12_PILE_FIELDS = ("Vcap (per element, lb)", "Mcap (per element, lb·ft)")


#: The Inputs-tree row Tutorial SEEP-1 sends the reader to. Added for every file,
#: not only a seepage one, which is what the page's Studio path says.
SEEP01_INPUT_CATEGORIES = ("Seep BC",)

#: The two Run-menu actions SEEP-1 names, in the same source-string convention as
#: ``LEM02_RUN_ACTIONS`` (Qt's ``&`` removed, the ellipsis kept). The Run action's
#: text follows the Mode switch, so the seepage label is read with the window put
#: into seepage mode and the mode restored afterwards.
SEEP01_BUILD_MESH_ACTION = "Build Mesh…"
SEEP01_RUN_ACTION_SEEP = "Run Seep…"

#: The Build Mesh dialog as SEEP-1 walks it: the two element types the page tells
#: the reader to choose between, the auto-size pair the target size comes from, the
#: refinement pair the tip step turns on, the thin-zone guarantee the page says to
#: leave alone, and the button that builds. A rewording of any of them turns those
#: steps into instructions to press something that is not there.
SEEP01_MESH_ELEMENTS = ("Linear triangles (tri3)", "Quadratic triangles (tri6)")
SEEP01_MESH_ROWS = ("Element type", "Size divisions", "Target element size",
                    "Refinement factor")
SEEP01_MESH_CHECKBOXES = ("Auto-size from geometry", "Refine near features",
                          "Refine thin zones")
SEEP01_MESH_BUILD = "Build"

#: The Run Seepage dialog: its title, the two controls the page explains one at a
#: time, and the button. The title is pinned because the page calls the dialog by
#: it while the menu action that opens it is worded differently.
SEEP01_RUN_TITLE = "Run Seepage"
SEEP01_RUN_ROWS = ("BC set", "Convergence tol")
SEEP01_RUN_BUTTON = "Run"

#: The seepage boundary-condition editor as the page's build step drives it: the
#: button that adds a head boundary, the two head types the page contrasts, the two
#: sets, and the value field — whose label carries the declared length unit, which
#: is how the page tells the reader their Time and Units declarations took.
SEEP01_BC_BUTTONS = ("Add head", "Add row")
SEEP01_BC_TYPES = ("head", "reservoir")
SEEP01_BC_TABS = ("Set 1", "Set 2 (rapid drawdown)")
SEEP01_BC_VALUE_LABEL = "Head value (m):"

#: The materials-editor usage toggle SEEP-1 tells the reader to leave ticked while
#: unticking the rest, and the two conductivity headers the table then shows — with
#: the unit suffix the model's declared Time puts on them, which the page says to
#: check for.
SEEP01_MATERIALS_TOGGLE = "Seepage"
SEEP01_MATERIALS_COLUMNS = ("k1 (m/day)", "k2 (m/day)")

#: Global parameters' Time row, which is inert on the LEM pages and load-bearing
#: here: it is what puts the unit suffixes above on the forms and the plots.
SEEP01_GLOBAL_TIME_ROW = "Time"

#: The seep-solution Display panel's contour-level control, and the value it opens
#: at. SEEP-1 tells the reader to change it FROM that value TO 10, because the
#: channel count the renderer derives only lands on a whole number at 10 here and
#: the page's read-back arithmetic is quoted at it. A changed default silently
#: turns that step into a no-op or into the wrong number, and the flow net the page
#: describes would stop being the one on the reader's screen.
SEEP01_LEVELS_ROW = "Contour levels"
SEEP01_LEVELS_DEFAULT = 20
SEEP01_LEVELS_STEP = 10

#: The Parametric dialog in its seepage mode, as the conductivity sweep uses it:
#: the mode whose target is a discharge, the parameter picker's own label for a
#: mode that can sweep a boundary as well as a material, and the four design rows
#: the page dictates values into.
SEEP01_PARAMETRIC_MODE = "Design (q target)"
SEEP01_PARAMETRIC_PICKER = "Material / BC"
SEEP01_PARAMETRIC_ROWS = ("From", "To", "Steps", "Target q")


def _seep01_editor_labels(mw):
    """The seepage controls Tutorial SEEP-1 drives, read on its own model.

    On its own model rather than on any other because every seepage control the
    page names is gated by what the file carries: the BC editor's list has an entry
    per boundary, the materials table's conductivity headers carry the unit the
    file's Time declaration supplies, and the Parametric dialog offers a seepage
    parameter set only for a model with a conductivity in it.

    The Run action is the one label here that belongs to the window rather than to
    a dialog, and its text follows the Mode switch. The window's mode is therefore
    moved to seepage, read, and put back — the leg leaves the project it was handed
    exactly as it found it, because the legs after it read their own models through
    the same window.
    """
    from PySide6.QtWidgets import (QCheckBox, QComboBox, QDialogButtonBox,
                                   QFormLayout, QLabel, QPushButton, QTabWidget,
                                   QTableWidget)

    from xslope.fileio import load_slope_data

    from studio.dialogs import BuildMeshDialog, RunSeepDialog, SensitivityDialog
    from studio.display_panels import SeepDisplayPanel
    from studio.editors import GlobalEditor, MaterialsEditor, SeepBcEditor

    fails = []
    data = _quiet(load_slope_data, SEEP01_FILE)

    action = getattr(mw, "act_build_mesh", None)
    if action is None:
        fails.append(f"MainWindow has no act_build_mesh, which Tutorial SEEP-1 "
                     f"calls {SEEP01_BUILD_MESH_ACTION!r}")
    elif action.text().replace("&", "") != SEEP01_BUILD_MESH_ACTION:
        fails.append(f"the act_build_mesh action reads "
                     f"{action.text().replace('&', '')!r}, not "
                     f"{SEEP01_BUILD_MESH_ACTION!r} — the label Tutorial SEEP-1 "
                     f"tells the reader to click")
    was = mw._mode
    try:
        mw._mode = "seep"
        mw._update_run_actions()
        got = mw.act_run.text().replace("&", "")
    finally:
        mw._mode = was
        mw._update_run_actions()
    if got != SEEP01_RUN_ACTION_SEEP:
        fails.append(f"the Run action reads {got!r} in seepage mode, not "
                     f"{SEEP01_RUN_ACTION_SEEP!r} — the label Tutorial SEEP-1 quotes")

    mesh = BuildMeshDialog(defaults={})
    elements = {mesh.element_type.itemText(i)
                for i in range(mesh.element_type.count())}
    for label in SEEP01_MESH_ELEMENTS:
        if label not in elements:
            fails.append(f"Build Mesh offers no {label!r} element type; Tutorial "
                         f"SEEP-1 tells the reader to choose between them. It "
                         f"offers {sorted(elements)}")
    rows = {lab.text() for lab in mesh.findChildren(QLabel)}
    for label in SEEP01_MESH_ROWS:
        if label not in rows:
            fails.append(f"Build Mesh has no {label!r} row; Tutorial SEEP-1 names "
                         f"it. Its rows read {sorted(r for r in rows if r)}")
    boxes = {b.text() for b in mesh.findChildren(QCheckBox)}
    for label in SEEP01_MESH_CHECKBOXES:
        if label not in boxes:
            fails.append(f"Build Mesh has no {label!r} checkbox; Tutorial SEEP-1 "
                         f"tells the reader what to do with it. Its checkboxes read "
                         f"{sorted(boxes)}")
    build = next((b for bb in mesh.findChildren(QDialogButtonBox)
                  for b in (bb.button(QDialogButtonBox.Ok),) if b is not None), None)
    if build is None or build.text() != SEEP01_MESH_BUILD:
        fails.append(f"Build Mesh's accept button reads "
                     f"{None if build is None else build.text()!r}, not "
                     f"{SEEP01_MESH_BUILD!r} — the label Tutorial SEEP-1 tells the "
                     f"reader to press")
    mesh.deleteLater()

    run = RunSeepDialog(defaults={}, slope_data=data)
    if run.windowTitle() != SEEP01_RUN_TITLE:
        fails.append(f"the seepage run dialog is titled {run.windowTitle()!r}, not "
                     f"{SEEP01_RUN_TITLE!r} — the name Tutorial SEEP-1 calls it by")
    rows = {lab.text() for lab in run.findChildren(QLabel)}
    for label in SEEP01_RUN_ROWS:
        if label not in rows:
            fails.append(f"the seepage run dialog has no {label!r} row; Tutorial "
                         f"SEEP-1 explains it. Its rows read "
                         f"{sorted(r for r in rows if r)}")
    if run._ok.text() != SEEP01_RUN_BUTTON:
        fails.append(f"the seepage run dialog's accept button reads "
                     f"{run._ok.text()!r}, not {SEEP01_RUN_BUTTON!r}")
    run.deleteLater()

    bc = SeepBcEditor().build(data, None, select=(0, 0))
    buttons = {b.text() for b in bc.findChildren(QPushButton)}
    for label in SEEP01_BC_BUTTONS:
        if label not in buttons:
            fails.append(f"the Seep BC editor has no {label!r} button; Tutorial "
                         f"SEEP-1 tells the reader to press it. Its buttons read "
                         f"{sorted(buttons)}")
    types = {combo.itemText(i) for combo in bc.findChildren(QComboBox)
             for i in range(combo.count())}
    for label in SEEP01_BC_TYPES:
        if label not in types:
            fails.append(f"the Seep BC editor offers no {label!r} head type; "
                         f"Tutorial SEEP-1 contrasts the two. It offers "
                         f"{sorted(types)}")
    tabs = [t.tabText(i) for t in bc.findChildren(QTabWidget)
            for i in range(t.count())]
    for label in SEEP01_BC_TABS:
        if label not in tabs:
            fails.append(f"the Seep BC editor has no {label!r} tab; Tutorial "
                         f"SEEP-1 names it. Its tabs read {tabs}")
    labels = {lab.text() for lab in bc.findChildren(QLabel)}
    if SEEP01_BC_VALUE_LABEL not in labels:
        fails.append(f"the Seep BC editor's head value is labeled something other "
                     f"than {SEEP01_BC_VALUE_LABEL!r}, which Tutorial SEEP-1 quotes "
                     f"with its unit to show the Time and Units declarations took")
    bc.deleteLater()

    glob = GlobalEditor().build(data, None)
    rows = {lab.text() for lab in glob.findChildren(QLabel)}
    if SEEP01_GLOBAL_TIME_ROW not in rows:
        fails.append(f"Global parameters has no {SEEP01_GLOBAL_TIME_ROW!r} row; "
                     f"Tutorial SEEP-1 tells the reader to set it. Its rows read "
                     f"{sorted(r for r in rows if r)}")
    glob.deleteLater()

    panel = SeepDisplayPanel(data.get("materials") or [])
    rows = {lab.text() for lab in panel.findChildren(QLabel)}
    if SEEP01_LEVELS_ROW not in rows:
        fails.append(f"the seep solution's Display panel has no "
                     f"{SEEP01_LEVELS_ROW!r} row; Tutorial SEEP-1 tells the reader "
                     f"to set it before reading the flow net. Its rows read "
                     f"{sorted(r for r in rows if r)}")
    elif panel.levels.value() != SEEP01_LEVELS_DEFAULT:
        fails.append(f"the Display panel's {SEEP01_LEVELS_ROW!r} opens at "
                     f"{panel.levels.value()}, not {SEEP01_LEVELS_DEFAULT} — "
                     f"Tutorial SEEP-1 says so and tells the reader to change it to "
                     f"{SEEP01_LEVELS_STEP}, at which the flow net's channel count "
                     f"comes out whole")
    elif not (panel.levels.minimum() <= SEEP01_LEVELS_STEP
              <= panel.levels.maximum()):
        fails.append(f"the Display panel's {SEEP01_LEVELS_ROW!r} cannot be set to "
                     f"{SEEP01_LEVELS_STEP}, which Tutorial SEEP-1 tells the reader "
                     f"to type; it accepts {panel.levels.minimum()} to "
                     f"{panel.levels.maximum()}")
    panel.deleteLater()

    # Read against the app's own toggle defaults, then ticked the way the page's
    # step leaves them: Seepage on, everything else off.
    with _default_editor_toggles():
        mats = MaterialsEditor().build(data, None)
        toggles = getattr(mats, "_toggles", None) or {}
        if SEEP01_MATERIALS_TOGGLE not in {cb.text() for cb in toggles.values()}:
            fails.append(f"the materials editor has no {SEEP01_MATERIALS_TOGGLE!r} "
                         f"usage toggle; Tutorial SEEP-1 tells the reader to leave "
                         f"it ticked. Its toggles read "
                         f"{sorted(cb.text() for cb in toggles.values())}")
        for tag, cb in toggles.items():
            cb.setChecked(cb.text() == SEEP01_MATERIALS_TOGGLE)
        mats._set_mode("table")
        headers = [t.horizontalHeaderItem(i).text()
                   if t.horizontalHeaderItem(i) else ""
                   for t in mats.findChildren(QTableWidget)
                   for i in range(t.columnCount())]
        for label in SEEP01_MATERIALS_COLUMNS:
            if label not in headers:
                fails.append(f"the materials table has no {label!r} column with the "
                             f"Seepage toggle on; Tutorial SEEP-1 tells the reader "
                             f"to type into it and to check the unit suffix. Its "
                             f"columns read {[h for h in headers if h]}")
        mats.deleteLater()

    sens = SensitivityDialog(defaults={"mode": "design"}, slope_data=data,
                             app_mode="seep")
    modes = {sens.mode.itemText(i) for i in range(sens.mode.count())}
    if SEEP01_PARAMETRIC_MODE not in modes:
        fails.append(f"the Parametric dialog offers no {SEEP01_PARAMETRIC_MODE!r} "
                     f"mode in seepage mode, which Tutorial SEEP-1 selects. It "
                     f"offers {sorted(modes)}")
    rows = set()
    for form in sens.findChildren(QFormLayout):
        for r in range(form.rowCount()):
            label = form.itemAt(r, QFormLayout.LabelRole)
            if label is not None and label.widget() is not None:
                rows.add(label.widget().text().replace("&", ""))
    for label in (SEEP01_PARAMETRIC_PICKER,) + SEEP01_PARAMETRIC_ROWS:
        if label not in rows:
            fails.append(f"the Parametric dialog has no {label!r} row in seepage "
                         f"mode; Tutorial SEEP-1 dictates a value into it. Its rows "
                         f"read {sorted(rows)}")
    properties = {sens.prop.itemText(i) for i in range(sens.prop.count())}
    if "k1" not in properties:
        fails.append(f"the Parametric dialog offers no 'k1' property for this "
                     f"model's soil; Tutorial SEEP-1 sweeps it. It offers "
                     f"{sorted(properties)}")
    sens.deleteLater()
    return fails


#: SEEP-2's model — the zoned dam whose three materials, exit face and unsaturated
#: parameters every pin below is read on. Its own file rather than SEEP-1's,
#: because SEEP-1's has one material, no exit face, and no unsaturated parameters
#: to show.
SEEP02_FILE = os.path.join(_REPO, "docs/seep/files/xslope_johnson_res.xlsx")

#: The materials table with the Seepage toggle on, as SEEP-2 reproduces it: the
#: three columns SEEP-1 does not teach — the model selector and both of the
#: linear front's parameters — plus the pair the van Genuchten and Gardner models
#: share. h0 carries the declared length unit; the conductivity columns do not,
#: because this file leaves Time blank, which the page tells the reader to expect.
SEEP02_MATERIALS_COLUMNS = ("k1 (ft/day)", "k2 (ft/day)", "alpha", "unsat", "kr0", "h0 (ft)",
                            "vg_a", "vg_n")

#: The three unsaturated models, as the `unsat` cell offers them. The page runs all
#: three and names each by the value that selects it, so a renamed option turns
#: three whole sections into instructions to pick something that is not there.
SEEP02_UNSAT_MODELS = ("lf", "vg", "gard")

#: The materials editor's list view, where SEEP-2 changes the model: the group the
#: page tells the reader to scroll to, the selector inside it, and the two views
#: the page names. The list view's row labels differ from the table's headers —
#: `unsat` is a column and **Unsat model** is a form row — so both are pinned.
SEEP02_MATERIALS_GROUP = "Conductivity"
SEEP02_MATERIALS_UNSAT_ROW = "Unsat model"
SEEP02_MATERIALS_VIEWS = ("Table view", "List view")

#: The exit face as the Seep BC editor lists it. SEEP-1 leaves this entry empty and
#: says so; SEEP-2 selects it, reads its two points off it, and builds a section
#: around what it does — so the wording of the list entry is load-bearing here in a
#: way it is not there.
SEEP02_BC_EXIT_ENTRY = "Exit face"

#: Where a specified-head boundary's VALUE is shown. SEEP-2 tells the reader that
#: the head rides in the list entry and that the table beside it holds points and
#: nothing else, so a value moved out of the list entry into a field of its own
#: makes that sentence wrong about the editor on the screen.
SEEP02_BC_HEAD_ENTRIES = ("Head 1  (h = 160.0)", "Head 2  (h = 100.0)")

#: Two dialog controls SEEP-2 explains and SEEP-1 does not name. Both are visible
#: in SEEP-2's own screenshots, which is why the page covers them: Quadrilateral
#: style is the group that greys out at a triangular element type, and Model checks
#: is the preflight panel filling half the Run Seepage dialog. Its all-clear line is
#: pinned too, because the page quotes it as what this model reports.
SEEP02_MESH_QUAD_GROUP = "Quadrilateral style"
SEEP02_RUN_CHECKS_GROUP = "Model checks"
SEEP02_RUN_CHECKS_OK = "No problems found for this run."


def _seep02_editor_labels():
    """The seepage controls Tutorial SEEP-2 drives, read on the zoned dam.

    Only the controls SEEP-1's pins do not already cover: that page and this one
    share the Build Mesh and Run Seepage dialogs, and re-reading their labels here
    would say the same thing twice and fail twice for one rewording.
    """
    from PySide6.QtWidgets import (QComboBox, QGroupBox, QLabel, QListWidget,
                                   QPushButton, QTableWidget)

    from xslope.fileio import load_slope_data

    from studio.dialogs import BuildMeshDialog, RunSeepDialog
    from studio.editors import MaterialsEditor, SeepBcEditor

    fails = []
    data = _quiet(load_slope_data, SEEP02_FILE)

    with _default_editor_toggles():
        mats = MaterialsEditor().build(data, None)
        for tag, cb in (getattr(mats, "_toggles", None) or {}).items():
            cb.setChecked(cb.text() == SEEP01_MATERIALS_TOGGLE)
        mats._set_mode("table")
        headers = [t.horizontalHeaderItem(i).text()
                   if t.horizontalHeaderItem(i) else ""
                   for t in mats.findChildren(QTableWidget)
                   for i in range(t.columnCount())]
        for label in SEEP02_MATERIALS_COLUMNS:
            if label not in headers:
                fails.append(f"the materials table has no {label!r} column with the "
                             f"Seepage toggle on; Tutorial SEEP-2 tabulates this "
                             f"model's rows against those headers. Its columns read "
                             f"{[h for h in headers if h]}")
        mats._set_mode("list")
        groups = {g.title() for g in mats.findChildren(QGroupBox)}
        if SEEP02_MATERIALS_GROUP not in groups:
            fails.append(f"the materials list view has no {SEEP02_MATERIALS_GROUP!r} "
                         f"group; Tutorial SEEP-2 tells the reader to scroll to it. "
                         f"Its groups read {sorted(groups)}")
        rows = {lab.text() for lab in mats.findChildren(QLabel)}
        if SEEP02_MATERIALS_UNSAT_ROW not in rows:
            fails.append(f"the materials list view has no "
                         f"{SEEP02_MATERIALS_UNSAT_ROW!r} row; Tutorial SEEP-2 "
                         f"calls it the selector for the three unsaturated models")
        for label in SEEP02_MATERIALS_VIEWS:
            if label not in {b.text() for b in mats.findChildren(QPushButton)}:
                fails.append(f"the materials editor has no {label!r} button; "
                             f"Tutorial SEEP-2 tells the reader to press it")
        options = {combo.itemText(i) for combo in mats.findChildren(QComboBox)
                   for i in range(combo.count())}
        for label in SEEP02_UNSAT_MODELS:
            if label not in options:
                fails.append(f"the materials editor offers no {label!r} unsaturated "
                             f"model; Tutorial SEEP-2 runs all three and names each "
                             f"by the value that selects it")
        mats.deleteLater()

    bc = SeepBcEditor().build(data, None, select=(0, 2))
    entries = [lw.item(i).text() for lw in bc.findChildren(QListWidget)
               for i in range(lw.count())]
    if SEEP02_BC_EXIT_ENTRY not in entries:
        fails.append(f"the Seep BC editor's list has no {SEEP02_BC_EXIT_ENTRY!r} "
                     f"entry for this model; Tutorial SEEP-2 tells the reader to "
                     f"select it and reads its points off it. Its entries read "
                     f"{entries}")
    for label in SEEP02_BC_HEAD_ENTRIES:
        if label not in entries:
            fails.append(f"the Seep BC editor's list has no {label!r} entry; "
                         f"Tutorial SEEP-2 says a specified head's value rides in "
                         f"the list entry and the table beside it holds points "
                         f"only. Its entries read {entries}")
    bc.deleteLater()

    mesh2 = BuildMeshDialog(defaults={})
    if SEEP02_MESH_QUAD_GROUP not in {g.title()
                                      for g in mesh2.findChildren(QGroupBox)}:
        fails.append(f"Build Mesh has no {SEEP02_MESH_QUAD_GROUP!r} group; "
                     f"Tutorial SEEP-2 tells the reader why it is greyed out at a "
                     f"triangular element type")
    mesh2.deleteLater()

    run2 = RunSeepDialog(defaults={}, slope_data=data)
    if SEEP02_RUN_CHECKS_GROUP not in {g.title()
                                       for g in run2.findChildren(QGroupBox)}:
        fails.append(f"the seepage run dialog has no {SEEP02_RUN_CHECKS_GROUP!r} "
                     f"group; Tutorial SEEP-2 says the preflight report fills the "
                     f"right half of it")
    if SEEP02_RUN_CHECKS_OK not in {lab.text() for lab in run2.findChildren(QLabel)}:
        fails.append(f"the seepage run dialog's checks do not read "
                     f"{SEEP02_RUN_CHECKS_OK!r} on this model; Tutorial SEEP-2 "
                     f"quotes that as what it reports for this file")
    run2.deleteLater()
    return fails


@contextlib.contextmanager
def _default_editor_toggles():
    """Read an editor against the app's OWN toggle defaults, not this machine's.

    Every usage toggle is written to QSettings as soon as it is clicked — including
    by an earlier check in this same process — so a pin on a coded default has to
    take the remembered state out of the way first and put it back afterwards,
    whatever happens.
    """
    from PySide6.QtCore import QSettings

    settings = QSettings("XSlope", "XSlope Studio")
    settings.beginGroup("editor_toggles")
    stashed = {k: settings.value(k) for k in settings.allKeys()}
    for key in stashed:
        settings.remove(key)
    settings.endGroup()
    settings.sync()
    try:
        yield
    finally:
        settings.beginGroup("editor_toggles")
        for key, value in stashed.items():
            settings.setValue(key, value)
        settings.endGroup()
        settings.sync()


def _lem11_reliability_labels():
    """The reliability controls Tutorial LEM-11 drives, read on its own model.

    On its own model rather than on any other, because every one of them is gated
    by the file carrying a standard deviation: the Monte Carlo engine, the σ
    summary's contents, and the Parametric dialog's variance plot type all appear
    only for a model with a σ in it, and LEM-11's is the only tutorial section that
    has one. Read on a σ-less file the pins would pass for the wrong reason or fail
    for one the page is not making.
    """
    from PySide6.QtWidgets import QCheckBox, QGroupBox, QLabel

    from xslope.fileio import load_slope_data

    from studio.dialogs import ReliabilityDialog, SensitivityDialog
    from studio.editors import MaterialsEditor

    fails = []
    data = _quiet(load_slope_data, LEM11_FILE)

    dlg = ReliabilityDialog(defaults={}, slope_data=data, app_mode="lem")
    engines = [dlg.engine.itemText(i) for i in range(dlg.engine.count())]
    for label in LEM11_ENGINES:
        if label not in engines:
            fails.append(f"the Reliability dialog offers no {label!r} engine; "
                         f"Tutorial LEM-11 gives it a section. It offers {engines}")
    boxes = {b.text() for b in dlg.findChildren(QCheckBox)}
    if LEM11_SEARCH_CHECKBOX_TAYLOR not in boxes:
        fails.append(f"the Reliability dialog has no "
                     f"{LEM11_SEARCH_CHECKBOX_TAYLOR!r} checkbox; Tutorial "
                     f"LEM-11 tells the reader to leave it ticked. Its "
                     f"checkboxes read {sorted(boxes)}")
    # The label is dynamic: switching the engine to Monte Carlo must restate
    # the checkbox in the sampling engines' terms — the page teaches the
    # difference off exactly this text.
    _idx = dlg.engine.findData("mc")
    if _idx >= 0:
        dlg.engine.setCurrentIndex(_idx)
        if dlg.search.text() != LEM11_SEARCH_CHECKBOX_MC:
            fails.append(f"on the Monte Carlo engine the search checkbox reads "
                         f"{dlg.search.text()!r}, not the "
                         f"{LEM11_SEARCH_CHECKBOX_MC!r} the page teaches")
        # The sampling combo must open on Latin hypercube (the page says the
        # dialog shows the default) and both items must map to values the
        # engines accept — the display/value order was once reversed, which
        # showed the raw 'lhs' and handed the engines an invalid string.
        if dlg.mc_sampling.currentText() != "Latin hypercube":
            fails.append(f"the MC sampling combo opens on "
                         f"{dlg.mc_sampling.currentText()!r}, not the "
                         f"'Latin hypercube' default Tutorial LEM-11 shows")
        _vals = {dlg.mc_sampling.itemData(i)
                 for i in range(dlg.mc_sampling.count())}
        if _vals != {"lhs", "random"}:
            fails.append(f"the MC sampling combo's item data is {_vals}, not "
                         f"the {{'lhs', 'random'}} the engines accept")
        dlg.engine.setCurrentIndex(dlg.engine.findData("taylor"))
    groups = {g.title() for g in dlg.findChildren(QGroupBox)}
    if LEM11_SIGMA_GROUP not in groups:
        fails.append(f"the Reliability dialog has no {LEM11_SIGMA_GROUP!r} box; "
                     f"Tutorial LEM-11 reads the run's inputs out of it. Its "
                     f"boxes read {sorted(groups)}")
    labels = {lab.text() for lab in dlg.findChildren(QLabel)}
    for row in LEM11_MC_ROWS:
        if row not in labels:
            fails.append(f"the Reliability dialog has no {row!r} row; Tutorial "
                         f"LEM-11 names it in its Monte Carlo step")
    # The page says the three Monte Carlo controls "come live" on the MC engine and
    # are inert on the Taylor series. That is the sentence, so that is the pin.
    if dlg.n_samples.isEnabled():
        fails.append("the MC sample count is live on the Taylor-series engine; "
                     "Tutorial LEM-11 says it comes live when Monte Carlo is "
                     "chosen")
    dlg.engine.setCurrentIndex(dlg.engine.findData("mc"))
    if not dlg.n_samples.isEnabled():
        fails.append("the MC sample count stays disabled on the Monte Carlo "
                     "engine; Tutorial LEM-11 tells the reader to read it there")
    dlg.deleteLater()

    # The toggle states persist in QSettings the moment anyone clicks one, so the
    # coded default is only readable with this machine's remembered state out of
    # the way — the same reason tools/capture_tutorial_screenshots._app_defaults
    # clears it. Without this the pin reads whoever ran Studio last, not the app.
    with _default_editor_toggles():
        mat = MaterialsEditor().build(data, None)
        toggles = {cb.text()
                   for cb in (getattr(mat, "_toggles", None) or {}).values()}
        if LEM11_RELIABILITY_TOGGLE not in toggles:
            fails.append(f"the materials editor has no "
                         f"{LEM11_RELIABILITY_TOGGLE!r} toggle; Tutorial LEM-11's "
                         f"first step is to tick it. Its toggles read "
                         f"{sorted(toggles)}")
        else:
            rel = (getattr(mat, "_toggles", None) or {}).get("rel")
            if rel is not None and rel.isChecked():
                fails.append("the materials editor opens with Reliability already "
                             "ticked; Tutorial LEM-11 tells the reader it is the "
                             "one that is off by default")
        mat.deleteLater()

    sens = SensitivityDialog(slope_data=data, app_mode="lem", defaults={})
    plots = [sens.plot_type.itemText(i) for i in range(sens.plot_type.count())]
    if LEM11_VARIANCE_PLOT not in plots:
        fails.append(f"the Parametric dialog offers no {LEM11_VARIANCE_PLOT!r} "
                     f"plot type on a model carrying sigmas; Tutorial LEM-11 "
                     f"tells the reader to select it. It offers {plots}")
    # The sigma-based plots ignore the sweep table, and the page says the
    # dialog shows that by graying it out — so it must actually gray out, and
    # come back when a table-driven plot type is selected again.
    _v = sens.plot_type.findData("variance")
    if _v >= 0:
        sens.plot_type.setCurrentIndex(_v)
        if sens.add_btn.isEnabled() or sens.table.isEnabled():
            fails.append("on the variance Pareto the sweep table stays live; "
                         "Tutorial LEM-11 says it grays out because the plot "
                         "ignores it")
        sens.plot_type.setCurrentIndex(sens.plot_type.findData("tornado"))
        if not (sens.add_btn.isEnabled() and sens.table.isEnabled()):
            fails.append("back on the tornado the sweep table stays disabled; "
                         "the gating must follow the plot type both ways")
    sens.deleteLater()
    return fails


def _lem01_editor_labels():
    """The materials editor as Tutorial LEM-1 drives it: opened, switched, added to.

    Driven rather than read, because the page's step is a *route* — the editor
    opens on one view and the reader is sent to the other — and only the route is
    wrong when the two views' buttons drift apart. The module remembers the last
    view for the session, so the remembered value is saved and restored: this
    check must not decide which view the next editor built here opens on.
    """
    from PySide6.QtWidgets import QLabel, QPushButton

    from xslope.fileio import load_slope_data

    import studio.editors as editors_mod
    from studio.editors import MaterialsEditor

    fails = []
    remembered = editors_mod._LAST_MATERIALS_VIEW
    try:
        if remembered != LEM01_MATERIALS_OPENS_ON:
            fails.append(f"the materials editor's first open lands on "
                         f"{remembered!r}, not {LEM01_MATERIALS_OPENS_ON!r} — "
                         f"Tutorial LEM-1 tells the reader which view they arrive "
                         f"in before sending them to the other one")
        data = _quiet(load_slope_data, LEM01_FILE)
        mats = MaterialsEditor().build(data, None)
        if mats._mode != remembered:
            fails.append(f"the materials editor opened on {mats._mode!r} with "
                         f"{remembered!r} remembered; Tutorial LEM-1's step assumes "
                         f"the view it opens on is the remembered one")
        # Both panes live in the same stack, so every button exists whichever view
        # is up; what the page's route turns on is which of them the reader can
        # actually press, which is visibility within the dialog.
        def shown():
            return {b.text() for b in mats.findChildren(QPushButton)
                    if b.text() and b.isVisibleTo(mats)}

        opened = shown()
        for label in LEM01_MATERIALS_VIEWS:
            if label not in opened:
                fails.append(f"the materials editor has no {label!r} button; "
                             f"Tutorial LEM-1 names both views. Its buttons read "
                             f"{sorted(opened)}")
        if LEM01_MATERIALS_ADD in opened:
            fails.append(f"the view the materials editor opens on already offers "
                         f"{LEM01_MATERIALS_ADD!r}; Tutorial LEM-1 sends the reader "
                         f"to the list view for it. Its buttons read {sorted(opened)}")
        mats.set_view_mode("list")
        listed = shown()
        if LEM01_MATERIALS_ADD not in listed:
            fails.append(f"the materials editor's list view has no "
                         f"{LEM01_MATERIALS_ADD!r} button; Tutorial LEM-1 tells the "
                         f"reader to press it there. Its buttons read "
                         f"{sorted(listed)}")
        mats.deleteLater()

        # The global-parameters form the page's first step fills, and photographs.
        from studio.editors import GlobalEditor

        glob = GlobalEditor().build(data, None)
        rows = {lab.text() for lab in glob.findChildren(QLabel)}
        for label in LEM01_GLOBAL_ROWS:
            if label not in rows:
                fails.append(f"Global parameters has no {label!r} row; Tutorial "
                             f"LEM-1's first step names it. Its rows read "
                             f"{sorted(r for r in rows if r)}")
        glob.deleteLater()
    finally:
        editors_mod._LAST_MATERIALS_VIEW = remembered
    return fails


def _lem03_editor_labels():
    """The three editors Tutorial LEM-3 drives, read for what it tells a reader to
    press — and the generator run, since the page prints its answer.

    The circles editor is driven the way the tutorial's step drives it: the model
    arrives with its circles dropped (the state the geometry steps leave), and the
    generator is *pressed*, because the summary line the page quotes exists only as
    the answer to a press.
    """
    from PySide6.QtWidgets import QComboBox, QLabel, QPushButton

    from xslope.fileio import load_slope_data

    from studio.editors import CirclesEditor, MaterialsEditor, ProfileEditor

    fails = []
    data = _quiet(load_slope_data, LEM03_FILE)

    mats = MaterialsEditor().build(data, None)
    buttons = {b.text() for b in mats.findChildren(QPushButton)}
    for label in LEM03_MATERIALS_BUTTONS:
        if label not in buttons:
            fails.append(f"the materials editor has no {label!r} button; Tutorial "
                         f"LEM-3 tells the reader to press it. Its buttons read "
                         f"{sorted(buttons)}")
    mats.deleteLater()

    prof = ProfileEditor().build(data, None)
    buttons = {b.text() for b in prof.findChildren(QPushButton)}
    for label in LEM03_PROFILE_BUTTONS:
        if label not in buttons:
            fails.append(f"the profile-lines editor has no {label!r} button; "
                         f"Tutorial LEM-3 tells the reader to press it. Its buttons "
                         f"read {sorted(buttons)}")
    labels = {lab.text() for lab in prof.findChildren(QLabel)}
    for name in LEM03_PROFILE_LABELS:
        if name not in labels:
            fails.append(f"the profile-lines editor has no {name!r} field; Tutorial "
                         f"LEM-3 names it. Its labels read "
                         f"{sorted(l for l in labels if l)}")
    choices = {combo.itemText(i) for combo in prof.findChildren(QComboBox)
               for i in range(combo.count())}
    for name in LEM03_PROFILE_MATERIALS:
        if name not in choices:
            fails.append(f"no profile-line material reads {name!r}, which Tutorial "
                         f"LEM-3 tells the reader to select. The choices read "
                         f"{sorted(choices)}")
    prof.deleteLater()

    circles_data = dict(data)
    circles_data["circles"], circles_data["circular"] = [], False
    circles = CirclesEditor().build(circles_data, None)
    _quiet(circles._run_generate)
    summary = " ".join(lab.text() for lab in circles.findChildren(QLabel))
    if LEM03_GENERATE_SUMMARY not in summary:
        fails.append(f"the starting-circle generator no longer reports "
                     f"{LEM03_GENERATE_SUMMARY!r} on LEM-3's section, which the page "
                     f"quotes. It reports {summary.strip()!r}")
    depths = tuple(round(float(c["Depth"]), 5) for c in circles.result_rows())
    if depths != LEM03_GENERATE_DEPTHS:
        fails.append(f"the generator proposes circles at depths {depths} on LEM-3's "
                     f"section; the page's audit table walks "
                     f"{LEM03_GENERATE_DEPTHS}")
    circles.deleteLater()
    return fails


def _lem05_editor_labels():
    """The non-circular editor and the Run LEM dialog, as Tutorial LEM-5 drives them.

    The generator is *pressed* rather than its rows pre-loaded, for the reason
    LEM-1's and LEM-3's shots document: the summary line the page quotes exists
    only as the answer to a press. It is pressed on a table emptied of the file's
    own surface, which is the state a reader who has entered nothing is in — and
    the state the page's figure shows.

    Run LEM is read on the file as it ships, because the fixed **Surface** label
    the page quotes is produced by a model carrying a non-circular surface and no
    circles; a model with both would show a chooser instead.
    """
    from PySide6.QtWidgets import QComboBox, QLabel, QPushButton, QTableWidget

    from xslope.fileio import load_slope_data

    from studio.dialogs import RunLemDialog
    from studio.editors import NonCircEditor

    fails = []
    data = _quiet(load_slope_data, LEM05_FILE)

    dlg = NonCircEditor().build(data, None)
    buttons = {b.text() for b in dlg.findChildren(QPushButton)}
    for label in LEM05_NONCIRC_BUTTONS:
        if label not in buttons:
            fails.append(f"the non-circular editor has no {label!r} button; "
                         f"Tutorial LEM-5 tells the reader to press it. Its buttons "
                         f"read {sorted(buttons)}")
    headers = [t.horizontalHeaderItem(i).text() if t.horizontalHeaderItem(i) else ""
               for t in dlg.findChildren(QTableWidget)
               for i in range(t.columnCount())]
    for name in LEM05_NONCIRC_HEADERS:
        if name not in headers:
            fails.append(f"the non-circular table has no {name!r} column; Tutorial "
                         f"LEM-5 dictates a value into it. Its columns read "
                         f"{headers}")
    movements = {combo.itemText(i) for combo in dlg.findChildren(QComboBox)
                 for i in range(combo.count())}
    for name in LEM05_NONCIRC_MOVEMENTS:
        if name not in movements:
            fails.append(f"no Movement setting reads {name!r}, which Tutorial LEM-5 "
                         f"tells the reader to choose. The options read "
                         f"{sorted(movements)}")
    dlg.deleteLater()

    empty = dict(data)
    empty["non_circ"] = []
    gen = NonCircEditor().build(empty, None)
    _quiet(gen._run_generate)
    summary = " ".join(lab.text() for lab in gen.findChildren(QLabel))
    for quoted in (LEM05_GENERATE_SUMMARY,) + LEM05_GENERATE_RAMPS:
        if quoted not in summary:
            fails.append(f"the weak-zone generator no longer reports {quoted!r} on "
                         f"LEM-5's section, which the page quotes. It reports "
                         f"{summary.strip()!r}")
    built = len(gen.result_rows())
    if built != LEM05_GENERATE_POINTS:
        fails.append(f"the weak-zone generator builds {built} points on LEM-5's "
                     f"section; the page quotes a summary announcing "
                     f"{LEM05_GENERATE_POINTS}")
    gen.deleteLater()

    run = RunLemDialog(defaults={}, slope_data=data)
    analyses = {run.analysis.itemText(i) for i in range(run.analysis.count())}
    for name in LEM05_RUN_ANALYSES:
        if name not in analyses:
            fails.append(f"Run LEM offers no {name!r} analysis, which Tutorial "
                         f"LEM-5 tells the reader to choose. It offers "
                         f"{sorted(analyses)}")
    if run.surface is not None:
        fails.append("Run LEM offers a Surface chooser on LEM-5's model; the page "
                     "says the row is a fixed label, which is what a model with one "
                     "surface family produces")
    labels = {lab.text() for lab in run.findChildren(QLabel)}
    if LEM05_RUN_SURFACE not in labels:
        fails.append(f"Run LEM's fixed Surface label does not read "
                     f"{LEM05_RUN_SURFACE!r} on LEM-5's model, which the page "
                     f"quotes. Its labels read {sorted(l for l in labels if l)}")
    run.deleteLater()

    # The end-ramp refusal, in the two texts the page quotes: what the Run box
    # says and what the Log pane says under it. Run through the same
    # ``run_lem_analysis`` Studio's LEM runner calls, on the same steep seed the
    # page describes, because the Run box shows exactly the AnalysisError it
    # raises.
    from xslope.search import AnalysisError, run_lem_analysis

    steep = dict(data)
    steep["non_circ"] = [dict(p) for p in data["non_circ"]]
    steep["non_circ"][-1]["X"] = LEM05_STEEP_EXIT_X
    log = io.StringIO()
    try:
        with contextlib.redirect_stdout(log):
            run_lem_analysis(steep, "spencer", analysis="auto_search",
                             surface="noncircular", num_slices=40)
        fails.append(f"a non-circular search seeded with a "
                     f"{LEM05_STEEP_EXIT_X:g}-exit surface (a 72.4 degree end ramp, "
                     f"past the search's 65 degree limit) now runs; Tutorial LEM-5 "
                     f"shows it being refused")
    except AnalysisError as e:
        if str(e) != LEM05_SEARCH_FAILED:
            fails.append(f"a search refused for a too-steep end ramp reports "
                         f"{str(e)!r}; Tutorial LEM-5 quotes {LEM05_SEARCH_FAILED!r} "
                         f"as the message the Run box shows")
    if LEM05_SEARCH_LOG not in log.getvalue():
        fails.append(f"the Log pane no longer carries {LEM05_SEARCH_LOG!r} when a "
                     f"search is seeded with a too-steep end ramp, which Tutorial "
                     f"LEM-5 quotes. It reads {log.getvalue().strip()!r}")
    return fails


def _lem04_editor_labels(mw):
    """The piezo and circles editors, the Log dock and Run LEM, as Tutorial
    LEM-4 drives them.

    Everything is read on the file as it ships. The model carries both unit
    weights, which the page's weight-split section reads as the file's own
    gsat column; the fixed **Circular** label the page describes is what a
    model carrying circles and no non-circular surface produces; and the
    circles editor's columns are the ones the page's pin step types the
    search's critical circle into.
    """
    from PySide6.QtWidgets import QLabel, QPushButton, QTableWidget, QTabWidget

    from xslope.fileio import load_slope_data

    from studio.dialogs import RunLemDialog
    from studio.editors import CirclesEditor, PiezoEditor

    fails = []
    data = _quiet(load_slope_data, LEM04_FILE)

    if not all(m.get("gamma_sat") for m in data["materials"]):
        fails.append("LEM-4's model no longer states a saturated unit weight for "
                     "every material; the page's weight-split section reads the "
                     "gsat column as the file ships it")

    circles = CirclesEditor().build(data, None)
    headers = set()
    for table in circles.findChildren(QTableWidget):
        headers.update(table.horizontalHeaderItem(i).text()
                       for i in range(table.columnCount())
                       if table.horizontalHeaderItem(i) is not None)
    for name in LEM04_CIRCLE_COLUMNS:
        if name not in headers:
            fails.append(f"the circles editor has no {name!r} column; Tutorial "
                         f"LEM-4 tells the reader to type the search's critical "
                         f"circle into it. Its columns read {sorted(headers)}")
    circles.deleteLater()

    dock = getattr(mw, "log_dock", None)
    if dock is None or dock.windowTitle() != LEM04_LOG_DOCK:
        fails.append(f"the main window has no {LEM04_LOG_DOCK!r} dock; Tutorial "
                     f"LEM-4 tells the reader to read the search's own statement "
                     f"of its critical circle there")

    piezo = PiezoEditor().build(data, None)
    tabs = [t.tabText(i) for t in piezo.findChildren(QTabWidget)
            for i in range(t.count())]
    for name in LEM04_PIEZO_TABS:
        if name not in tabs:
            fails.append(f"the piezometric-lines editor has no {name!r} tab; "
                         f"Tutorial LEM-4 names it. Its tabs read {tabs}")
    buttons = {b.text() for b in piezo.findChildren(QPushButton)}
    if LEM04_PIEZO_ADD not in buttons:
        fails.append(f"the piezometric-lines editor has no {LEM04_PIEZO_ADD!r} "
                     f"button; Tutorial LEM-4 tells the reader to press it eight "
                     f"times. Its buttons read {sorted(buttons)}")
    piezo.deleteLater()

    run = RunLemDialog(defaults={}, slope_data=data)
    analyses = {run.analysis.itemText(i) for i in range(run.analysis.count())}
    if LEM04_RUN_ANALYSIS not in analyses:
        fails.append(f"Run LEM offers no {LEM04_RUN_ANALYSIS!r} analysis on "
                     f"LEM-4's model, which its pivot step chooses. It offers "
                     f"{sorted(analyses)}")
    if run.surface is not None:
        fails.append("Run LEM offers a Surface chooser on LEM-4's model; the page "
                     "says the row is a fixed label, which is what a circles-only "
                     "model produces")
    labels = {lab.text() for lab in run.findChildren(QLabel)}
    if LEM04_RUN_SURFACE not in labels:
        fails.append(f"Run LEM's fixed Surface label does not read "
                     f"{LEM04_RUN_SURFACE!r} on LEM-4's model, which the page "
                     f"quotes. Its labels read {sorted(l for l in labels if l)}")
    run.deleteLater()
    return fails


def _lem06_editor_labels(mw):
    """The polygons editor and the composite option, as Tutorial LEM-6 drives them.

    Read on a polygon-based model, because that is the only kind whose Polygons
    row opens an editor — which is itself what the page's Studio path says, so
    the tree row is checked in both states: editable on this model, and inert on
    a project that has no zones yet.

    The refusal is *provoked* rather than quoted from a constant: the page's
    composite section is built on a circle that will not fit, so the check pushes
    the file's deeper starting circle 1.2 ft below the base and confirms the run
    still stops on the message the page prints, and still solves with the option
    the page tells the reader to tick.
    """
    from PySide6.QtWidgets import QCheckBox, QComboBox, QLabel, QPushButton

    from xslope.fileio import load_slope_data
    from xslope.search import AnalysisError, run_lem_analysis

    from studio.dialogs import RunLemDialog
    from studio.editors import PolygonEditor

    fails = []

    def _polygons_row():
        """The Inputs tree's Polygons row: (count, editor category | None)."""
        from studio.main_window import CATEGORY_ROLE as ROLE
        tree = mw.inputs_tree
        for i in range(tree.topLevelItemCount()):
            item = tree.topLevelItem(i)
            if item.text(0) == "Polygons":
                return item.text(1), item.data(0, ROLE)
        return None, None

    mw.doc._dirty = False
    mw.new_project()
    _, category = _polygons_row()
    if category != "polygons":
        fails.append("a project started from File > New does not open the "
                     "Polygons editor; Tutorial LEM-6's Studio path builds the "
                     "model from scratch through that row")
    # The guard the empty-project liveness must NOT have removed: on a
    # profile-based file, polygons are derived and the row stays inert.
    mw.doc._dirty = False
    _quiet(mw.open_path, os.path.join(_REPO,
           "docs/lem/files/xslope_simple_mult_layers.xlsx"))
    _, category = _polygons_row()
    if category is not None:
        fails.append("the Polygons row opens an editor on a profile-based "
                     "model; derived polygons must stay read-only there")
    mw.doc._dirty = False
    _quiet(mw.open_path, LEM06_FILE)
    count, category = _polygons_row()
    if category != "polygons":
        fails.append("the Polygons row does not open an editor on LEM-6's own "
                     "model; its Studio path is written on the model being open")
    if count != "2":
        fails.append(f"the Polygons row reads {count!r} on LEM-6's model; the page "
                     f"says it reads '2' once the model is open")

    data = _quiet(load_slope_data, LEM06_FILE)
    if len(data.get("polygons") or []) != 2 or data.get("profile_lines"):
        fails.append("LEM-6's model is no longer two polygons and no profile "
                     "lines; the page's whole geometry section reads it as the "
                     "polygon-input case")

    poly = PolygonEditor().build(data, None, select=1)
    buttons = {b.text() for b in poly.findChildren(QPushButton)}
    for label in LEM06_POLYGON_BUTTONS:
        if label not in buttons:
            fails.append(f"the polygons editor has no {label!r} button; Tutorial "
                         f"LEM-6 tells the reader to press it. Its buttons read "
                         f"{sorted(buttons)}")
    labels = {lab.text() for lab in poly.findChildren(QLabel)}
    for name in LEM06_POLYGON_LABELS:
        if name not in labels:
            fails.append(f"the polygons editor has no {name!r} field; Tutorial "
                         f"LEM-6 names it. Its labels read "
                         f"{sorted(l for l in labels if l)}")
    if LEM06_POLYGON_NO_MAX_DEPTH in labels:
        fails.append(f"the polygons editor now offers "
                     f"{LEM06_POLYGON_NO_MAX_DEPTH!r}; Tutorial LEM-6 tells the "
                     f"reader it is absent, because a polygon model's bottom "
                     f"boundary is drawn rather than typed")
    help_text = " ".join(lab.text() for lab in poly.findChildren(QLabel))
    if LEM06_POLYGON_HELP not in help_text:
        fails.append(f"the polygons editor no longer states "
                     f"{LEM06_POLYGON_HELP!r}, which Tutorial LEM-6 quotes as the "
                     f"closed-region rule")
    items = [poly.list.item(i).text() for i in range(poly.list.count())]
    if LEM06_POLYGON_ITEM not in items:
        fails.append(f"the polygons editor's zone list has no "
                     f"{LEM06_POLYGON_ITEM!r} entry, which Tutorial LEM-6 tells "
                     f"the reader to click. It reads {items}")
    poly.deleteLater()

    run = RunLemDialog(defaults={}, slope_data=data)
    boxes = {b.text() for b in run.findChildren(QCheckBox)}
    if LEM06_COMPOSITE_CHECKBOX not in boxes:
        fails.append(f"Run LEM has no {LEM06_COMPOSITE_CHECKBOX!r} checkbox; "
                     f"Tutorial LEM-6 tells the reader to tick it. Its checkboxes "
                     f"read {sorted(boxes)}")
    run.deleteLater()

    # The circle the page pushes below the base: refused as an ordinary circle,
    # truncated against the base when composite surfaces are allowed.
    deep = dict(data)
    deep["circles"] = [dict(LEM06_DEEP_CIRCLE)]
    try:
        _quiet(run_lem_analysis, deep, "spencer", analysis="single_surface",
               surface="circular", num_slices=40, composite=False)
        fails.append("a circle below LEM-6's domain floor now solves without the "
                     "composite option; the page's composite section is built on "
                     "its refusal")
    except AnalysisError as exc:
        if LEM06_DOMAIN_REFUSAL not in str(exc):
            fails.append(f"a circle below LEM-6's domain floor is refused with "
                         f"{str(exc)!r}, not {LEM06_DOMAIN_REFUSAL!r} — the "
                         f"message the page quotes")
    bundle = _quiet(run_lem_analysis, deep, "spencer", analysis="single_surface",
                    surface="circular", num_slices=40, composite=True)
    if (bundle.get("results") or {}).get("FS") is None:
        fails.append("the same circle no longer solves with composite surfaces "
                     "on; Tutorial LEM-6 reports the truncated surface's factor "
                     "of safety")
    return fails


def _lem08_editor_labels(mw):
    """The reinforcement editor and the Type presets, as Tutorial LEM-8 uses them.

    Read on LEM-8's own model, because the page's Studio step is written on a
    project that already carries the six lines: the Inputs tree row counts them,
    and the editor's list has one entry per line. Runs after the LEM-6 leg, since
    both change which project the window holds.

    The preset table is read from the template's own lookup block rather than
    from a constant of ours, so the page's table and the sheet's drop-down cannot
    drift apart without this failing.
    """
    import openpyxl
    from PySide6.QtWidgets import QComboBox, QGroupBox, QLabel, QPushButton

    from xslope.fileio import load_slope_data

    from studio.editors import ReinforcementEditor

    fails = []

    mw.doc._dirty = False
    _quiet(mw.open_path, LEM08_FILE)
    count, category = None, None
    from studio.main_window import CATEGORY_ROLE as ROLE
    tree = mw.inputs_tree
    for i in range(tree.topLevelItemCount()):
        item = tree.topLevelItem(i)
        if item.text(0) == "Reinforcement lines":
            count, category = item.text(1), item.data(0, ROLE)
            break
    if category != "reinforce":
        fails.append("the Inputs tree has no 'Reinforcement lines' row opening an "
                     "editor on LEM-8's model; its Studio path tells the reader to "
                     "click it")
    if count != "6":
        fails.append(f"the Reinforcement lines row reads {count!r} on LEM-8's "
                     f"model; the page builds six geogrid layers")

    data = _quiet(load_slope_data, LEM08_FILE)
    lines = data.get("reinforcement_lines") or []
    if len(lines) != 6:
        fails.append(f"LEM-8's model carries {len(lines)} reinforcement lines, not "
                     f"the six the page enters and reads its crossings against")

    dlg = ReinforcementEditor().build(data, None)
    if dlg.windowTitle() != "Reinforcement":
        fails.append(f"the reinforcement editor is titled {dlg.windowTitle()!r}, "
                     f"not 'Reinforcement'")
    buttons = {b.text() for b in dlg.findChildren(QPushButton)}
    for label in LEM08_REINF_BUTTONS + LEM08_REINF_VIEWS:
        if label not in buttons:
            fails.append(f"the reinforcement editor has no {label!r} button; "
                         f"Tutorial LEM-8 names it. Its buttons read "
                         f"{sorted(buttons)}")
    groups = {g.title() for g in dlg.findChildren(QGroupBox)}
    for title in LEM08_REINF_GROUPS:
        if title not in groups:
            fails.append(f"the reinforcement editor's list view has no {title!r} "
                         f"group; Tutorial LEM-8 walks the form group by group. "
                         f"Its groups read {sorted(groups)}")
    labels = [lab.text() for lab in dlg.findChildren(QLabel)]
    for field in LEM08_REINF_FIELDS:
        if not any(text == field or text.startswith(field + " ")
                   for text in labels):
            fails.append(f"the reinforcement editor has no {field!r} field; "
                         f"Tutorial LEM-8 tells the reader what to put in it")
    if LEM08_REINF_TMAX_LABEL not in labels:
        fails.append(f"the reinforcement editor no longer labels Tmax "
                     f"{LEM08_REINF_TMAX_LABEL!r}; Tutorial LEM-8 quotes that "
                     f"label as how a per-element capacity shows up as the wrong "
                     f"quantity")
    combo_items = set()
    for combo in dlg.findChildren(QComboBox):
        combo_items.update(combo.itemText(i) for i in range(combo.count()))
    for choice in ("geosynthetic", "nail", "tieback", "anchor",
                   "tangent", "axial", "active", "passive"):
        if choice not in combo_items:
            fails.append(f"the reinforcement editor offers no {choice!r} choice; "
                         f"Tutorial LEM-8's preset table names it")
    dlg.deleteLater()

    # The Type preset table, read from the template's own lookup block.
    sheet = openpyxl.load_workbook(TEMPLATE_FILE)["reinforce"]
    presets = tuple(tuple(str(sheet.cell(row=r, column=c).value)
                          for c in (26, 27, 28))          # Z, AA, AB
                    for r in range(8, 12))
    if presets != LEM08_TYPE_PRESETS:
        fails.append(f"the template's support-type presets read {presets}, not "
                     f"{LEM08_TYPE_PRESETS} — the table Tutorial LEM-8 reproduces")
    return fails


def _lem09_editor_labels(mw):
    """The piles editor and the Anchor preset, as Tutorial LEM-9 uses them.

    Read on LEM-9's own model, because the page's Studio step is written on a
    project that already carries the two anchors and the soldier pile: the Inputs
    tree counts them, and each editor has a row to read its fields off. Runs with
    the other legs that change which project the window holds.

    The preset is checked through the loader rather than against a constant of
    ours, so the page's claim — pick Anchor and Dir/Appl answer on their own —
    fails here if the preset stops filling them.
    """
    from PySide6.QtWidgets import QGroupBox, QLabel

    from xslope.fileio import load_slope_data

    from studio.editors import PilesEditor

    fails = []

    mw.doc._dirty = False
    _quiet(mw.open_path, LEM09_FILE)
    from studio.main_window import CATEGORY_ROLE as ROLE
    tree = mw.inputs_tree
    count, category = None, None
    for i in range(tree.topLevelItemCount()):
        item = tree.topLevelItem(i)
        if item.text(0) == LEM09_INPUT_CATEGORY:
            count, category = item.text(1), item.data(0, ROLE)
            break
    if category != "piles":
        fails.append(f"the Inputs tree has no {LEM09_INPUT_CATEGORY!r} row opening "
                     f"an editor on LEM-9's model; its Studio path tells the reader "
                     f"to click it")
    if count != "1":
        fails.append(f"the Piles row reads {count!r} on LEM-9's model; the page "
                     f"enters one soldier pile")

    data = _quiet(load_slope_data, LEM09_FILE)
    lines = data.get("reinforcement_lines") or []
    if len(lines) != 2:
        fails.append(f"LEM-9's model carries {len(lines)} reinforcement lines, not "
                     f"the two tiebacks the page enters")
    for line in lines:
        got = (str(line.get("type")), str(line.get("dir")), str(line.get("appl")))
        if got != LEM09_ANCHOR_PRESET:
            fails.append(f"an Anchor line on LEM-9's model reads {got}, not "
                         f"{LEM09_ANCHOR_PRESET} — the page tells the reader to pick "
                         f"the Type and let Dir and Appl fill themselves")

    dlg = PilesEditor().build(data, None)
    if dlg.windowTitle() != "Piles":
        fails.append(f"the piles editor is titled {dlg.windowTitle()!r}, not 'Piles'")
    groups = {g.title() for g in dlg.findChildren(QGroupBox)}
    for title in LEM09_PILE_GROUPS:
        if title not in groups:
            fails.append(f"the piles editor's list view has no {title!r} group; "
                         f"Tutorial LEM-9 names the form's groups in order. Its "
                         f"groups read {sorted(groups)}")
    labels = [lab.text() for lab in dlg.findChildren(QLabel)]
    for field in LEM09_PILE_FIELDS:
        if not any(text == field or text.startswith(field + " ")
                   for text in labels):
            fails.append(f"the piles editor has no {field!r} field; Tutorial LEM-9 "
                         f"tells the reader what to put in it")
    dlg.deleteLater()
    return fails


def _lem12_pile_force_labels(mw):
    """Where Tutorial LEM-12 sends the reader to READ the computed pile force.

    The force the Ito & Matsui calculation produces is not on any plot and not in
    the Log pane: the page's "Where the computed force appears" section sends the
    reader to the Analysis Report for it, and quotes three strings verbatim — the
    menu action that writes the report, the word the Piles table prints in the H
    cell of a row whose force is computed, and the slice table's own column header
    and legend sentence for that force. Each is pinned at its source rather than
    by building a report, which would cost a search and a document per run.

    The two capacity fields are pinned with their unit suffixes, which the page
    quotes in full because they are what say the capacities belong to one shaft
    and not to a foot of slope.
    """
    from PySide6.QtWidgets import QLabel

    from xslope.columns import SLICE_COLUMNS, header
    from xslope.fileio import load_slope_data
    from xslope.report import _pile_fields

    from studio.editors import PilesEditor

    fails = []

    action = getattr(mw, LEM12_REPORT_ACTION[0], None)
    if action is None:
        fails.append(f"MainWindow has no {LEM12_REPORT_ACTION[0]}, which Tutorial "
                     f"LEM-12 calls {LEM12_REPORT_ACTION[1]!r}")
    elif action.text().replace("&", "") != LEM12_REPORT_ACTION[1]:
        fails.append(f"the {LEM12_REPORT_ACTION[0]} action reads "
                     f"{action.text().replace('&', '')!r}, not "
                     f"{LEM12_REPORT_ACTION[1]!r} — the label Tutorial LEM-12 quotes")

    # The Piles table's H cell on a row that states no force. Driven through the
    # report's own formatter, so a change from "computed" to anything else is the
    # failure rather than a constant of ours agreeing with itself.
    fmt = _pile_fields({"unit_system": "imperial"})["H"][2]
    got = fmt({"H": None})
    if got != LEM12_COMPUTED_CELL:
        fails.append(f"the report's Piles table prints {got!r} for a blank H, not "
                     f"{LEM12_COMPUTED_CELL!r} — the word Tutorial LEM-12 tells the "
                     f"reader to look for")

    column = next((c for c in SLICE_COLUMNS if c.key == "h_pile"), None)
    got = header(column, {"force_per_len": "lb/ft"}) if column else None
    if got != LEM12_HP_HEADER:
        fails.append(f"the slice table's pile-force column is headed {got!r}, not "
                     f"{LEM12_HP_HEADER!r} — the column Tutorial LEM-12 names")
    legend = column.description if column else None
    if legend != LEM12_HP_LEGEND:
        fails.append(f"the slice table's pile-force legend reads {legend!r}, not "
                     f"{LEM12_HP_LEGEND!r} — the sentence Tutorial LEM-12 quotes")

    data = _quiet(load_slope_data, LEM12_FILE)
    if not all(p.get("H") is None for p in (data.get("pile_lines") or [])):
        fails.append("LEM-12's model states a pile force on some row; the page's "
                     "whole subject is that both rows leave H blank")
    dlg = PilesEditor().build(data, None)
    labels = [lab.text() for lab in dlg.findChildren(QLabel)]
    for field in LEM12_PILE_FIELDS:
        if field not in labels:
            fails.append(f"the piles editor has no {field!r} field label; Tutorial "
                         f"LEM-12 quotes it to say the capacity is per shaft. Its "
                         f"labels read {sorted(l for l in labels if l)}")
    dlg.deleteLater()
    return fails


def _lem02_editor_labels(mw):
    """The two dialogs Tutorial LEM-2 drives, read for the labels it quotes.

    Both are built on the tutorial's own model rather than on an empty one: the
    distributed-loads editor's Direction combo belongs to a selected load block,
    and there is nothing to select until the model carries one.
    """
    from PySide6.QtWidgets import (QCheckBox, QComboBox, QFormLayout, QLabel,
                                   QPushButton, QTableWidget, QTabWidget)

    from xslope.fileio import load_slope_data

    from studio.dialogs import SensitivityDialog
    from studio.editors import DloadsEditor, GlobalEditor, LineLoadsEditor

    fails = []
    data = _quiet(load_slope_data, LEM02_FILE)

    dlg = DloadsEditor().build(data, None)
    buttons = {b.text() for b in dlg.findChildren(QPushButton)}
    for label in LEM02_DLOAD_BUTTONS:
        if label not in buttons:
            fails.append(f"the distributed-loads editor has no {label!r} button; "
                         f"Tutorial LEM-2 tells the reader to press it. Its buttons "
                         f"read {sorted(buttons)}")
    directions = {combo.itemText(i)
                  for combo in dlg.findChildren(QComboBox)
                  for i in range(combo.count())}
    if LEM02_DLOAD_DIRECTION not in directions:
        fails.append(f"no distributed-load Direction reads "
                     f"{LEM02_DLOAD_DIRECTION!r}, which Tutorial LEM-2 quotes. The "
                     f"options read {sorted(directions)}")
    tabs = [t.tabText(i) for t in dlg.findChildren(QTabWidget)
            for i in range(t.count())]
    for name in LEM02_DLOAD_TABS:
        if name not in tabs:
            fails.append(f"the distributed-loads editor has no {name!r} tab; "
                         f"Tutorial LEM-2 names it. Its tabs read {tabs}")
    dl_labels = {lab.text() for lab in dlg.findChildren(QLabel)}
    if LEM02_DLOAD_DIRECTION_LABEL not in dl_labels:
        fails.append(f"the distributed-loads editor labels its direction chooser "
                     f"something other than {LEM02_DLOAD_DIRECTION_LABEL!r}, which "
                     f"Tutorial LEM-2 tells the reader to leave alone")
    dlg.deleteLater()

    glob = GlobalEditor().build(data, None)
    rows = {lab.text() for lab in glob.findChildren(QLabel)}
    if LEM02_GLOBAL_ROW not in rows:
        fails.append(f"Global parameters has no {LEM02_GLOBAL_ROW!r} row; Tutorial "
                     f"LEM-2 tells the reader to set it. Its rows read "
                     f"{sorted(r for r in rows if r)}")
    glob.deleteLater()

    # The line-loads editor is built on the model with its distributed load
    # replaced by the page's line load — the state the step it photographs ends in.
    ll_data = dict(data)
    ll_data["dloads"], ll_data["dload_dirs"] = [], []
    ll_data["line_loads"] = [{"x": 30.0, "y": 20.0, "P": 7500.0, "angle": -90.0,
                              "label": "footing"}]
    lloads = LineLoadsEditor().build(ll_data, None)
    headers = [t.horizontalHeaderItem(i).text() if t.horizontalHeaderItem(i) else ""
               for t in lloads.findChildren(QTableWidget)
               for i in range(t.columnCount())]
    for name in LEM02_LLOAD_HEADERS:
        if name not in headers:
            fails.append(f"the line-loads table has no {name!r} column; Tutorial "
                         f"LEM-2 dictates a value into it. Its columns read "
                         f"{headers}")
    lloads.deleteLater()

    sens = SensitivityDialog(defaults={"mode": "design"}, slope_data=data)
    modes = {sens.mode.itemText(i) for i in range(sens.mode.count())}
    if LEM02_DESIGN_MODE not in modes:
        fails.append(f"the Parametric dialog offers no {LEM02_DESIGN_MODE!r} mode, "
                     f"which Tutorial LEM-2 selects. It offers {sorted(modes)}")
    rows = set()
    for form in sens.findChildren(QFormLayout):
        for r in range(form.rowCount()):
            label = form.itemAt(r, QFormLayout.LabelRole)
            if label is not None and label.widget() is not None:
                rows.add(label.widget().text().replace("&", ""))
    for name in LEM02_DESIGN_ROWS:
        if name not in rows:
            fails.append(f"the Parametric dialog has no {name!r} row; Tutorial LEM-2 "
                         f"names it. Its rows read {sorted(rows)}")
    boxes = {c.text() for c in sens.findChildren(QCheckBox)}
    if LEM02_DESIGN_SEARCH not in boxes:
        fails.append(f"the Parametric dialog has no {LEM02_DESIGN_SEARCH!r} "
                     f"checkbox; Tutorial LEM-2 tells the reader to leave it "
                     f"ticked. Its checkboxes read {sorted(boxes)}")
    if sens._ok.text() != LEM02_DESIGN_RUN:
        fails.append(f"the Parametric dialog's accept button reads "
                     f"{sens._ok.text()!r}, not {LEM02_DESIGN_RUN!r} — the label "
                     f"Tutorial LEM-2 tells the reader to press")
    sens.deleteLater()
    return fails


def test_tutorial_labels():
    """Studio still reads the way the tutorials say it does."""
    from PySide6.QtWidgets import QApplication

    from studio.dialogs import UnpackPackageDialog
    from studio.main_window import MainWindow

    QApplication.instance() or QApplication([])
    fails = []

    mw = MainWindow()
    mw._add_recent = lambda path: None            # never touch the user's settings
    try:
        for attr, label in T0_FILE_ACTIONS.items():
            action = getattr(mw, attr, None)
            if action is None:
                fails.append(f"MainWindow has no {attr}, which Tutorial 0 calls "
                             f"{label!r}")
                continue
            got = action.text().replace("&", "")
            if got != label:
                fails.append(f"the {attr} action reads {got!r}, not {label!r} — the "
                             f"label Tutorial 0 quotes")

        _quiet(mw.open_path, SINGLE_FILE)
        tree = mw.inputs_tree
        rows = [tree.topLevelItem(i).text(0) for i in range(tree.topLevelItemCount())]
        for name in T0_INPUT_CATEGORIES:
            if name not in rows:
                fails.append(f"the Inputs tree has no {name!r} row; Tutorial 0 tells "
                             f"the reader to click it. It reads {rows}")
        for name in LEM02_INPUT_CATEGORIES:
            if name not in rows:
                fails.append(f"the Inputs tree has no {name!r} row; Tutorial LEM-2 "
                             f"tells the reader to click it. It reads {rows}")
        for name in LEM05_INPUT_CATEGORIES:
            if name not in rows:
                fails.append(f"the Inputs tree has no {name!r} row; Tutorial LEM-5 "
                             f"tells the reader to click it. It reads {rows}")
        # Read on the same limit-equilibrium file as the rows above, which is the
        # point: SEEP-1's Studio path says the tree is the same in every mode.
        for name in SEEP01_INPUT_CATEGORIES:
            if name not in rows:
                fails.append(f"the Inputs tree has no {name!r} row; Tutorial SEEP-1 "
                             f"tells the reader to click it. It reads {rows}")

        for attr, label in LEM02_RUN_ACTIONS.items():
            action = getattr(mw, attr, None)
            if action is None:
                fails.append(f"MainWindow has no {attr}, which Tutorial LEM-2 calls "
                             f"{label!r}")
            elif action.text().replace("&", "") != label:
                fails.append(f"the {attr} action reads "
                             f"{action.text().replace('&', '')!r}, not {label!r} — "
                             f"the label Tutorial LEM-2 quotes")

        fails += _lem01_editor_labels()
        fails += _lem02_editor_labels(mw)
        fails += _lem03_editor_labels()
        fails += _lem05_editor_labels()
        fails += _lem04_editor_labels(mw)
        for attr, label in LEM11_RUN_ACTIONS.items():
            action = getattr(mw, attr, None)
            if action is None:
                fails.append(f"MainWindow has no {attr}, which Tutorial LEM-11 "
                             f"calls {label!r}")
            elif action.text().replace("&", "") != label:
                fails.append(f"the {attr} action reads "
                             f"{action.text().replace('&', '')!r}, not {label!r} — "
                             f"the label Tutorial LEM-11 quotes")

        fails += _lem10_run_labels()
        fails += _lem11_reliability_labels()
        fails += _lem12_pile_force_labels(mw)
        fails += _seep01_editor_labels(mw)
        fails += _seep02_editor_labels()
        # Last of the editor legs: these are the ones that change which project
        # the window holds, and each opens the model its own pins are read on.
        fails += _lem06_editor_labels(mw)
        fails += _lem08_editor_labels(mw)
        fails += _lem09_editor_labels(mw)

        from PySide6.QtWidgets import QPushButton

        dlg = UnpackPackageDialog(os.path.join(_tmp(), "sample" + PACKAGE_EXT), mw)
        # The dialog offers the already-exists choice only for a destination that is
        # there, which is the case the tutorial describes. (Hidden rather than
        # visible is the question: nothing here is on a screen, so isVisible() is
        # False for every button in an unshown dialog.)
        dlg.dest.setText(_tmp())
        dlg._refresh()
        buttons = {b.text(): b for b in dlg.findChildren(QPushButton)}
        for label in T0_PACKAGE_BUTTONS:
            button = buttons.get(label)
            if button is None:
                fails.append(f"the package dialog has no {label!r} button; Tutorial 0 "
                             f"names it. Its buttons read {sorted(buttons)}")
            elif button.isHidden():
                fails.append(f"{label!r} is hidden on a destination that already "
                             f"exists, which is the case the tutorial describes")
        dlg.deleteLater()
    finally:
        mw.doc._dirty = False
        mw.close()
    return fails


def test_split_code_spans():
    """No inline code span may wrap across a source line: the renderer keeps the
    break inside the backticks, splitting phrases like `Latin hypercube` mid-word
    on the published page. Fenced blocks are skipped; only single-backtick spans
    that open and close on different lines fail."""
    fails = []
    for path in sorted(glob.glob(os.path.join(_REPO, "docs", "**", "*.md"),
                                 recursive=True)):
        text = open(path, encoding="utf-8").read()
        # strip fenced code blocks so their backticks cannot pair across lines
        text = re.sub(r"```.*?```", "", text, flags=re.S)
        for m in re.finditer(r"(?<!`)`(?!`)([^`]*)`(?!`)", text):
            if "\n" in m.group(1):
                rel = os.path.relpath(path, _REPO)
                line = text[:m.start()].count("\n") + 1
                fails.append("%s:%d code span wraps across a line: `%s`"
                             % (rel, line, m.group(1).replace("\n", "\\n")))
    return fails


CHECKS = [
    ("A. one verb, and it is named", test_verb_gate),
    ("B. only XSLOPE's own sites", test_allowlist),
    ("C. redirects are re-checked", test_redirect_guard),
    ("D. the saved name is a name", test_saved_name),
    ("E. the docs build packages and pairs", test_docs_build),
    ("F. refuse, ask, fetch, open", test_studio_flow),
    ("G. argv and FileOpen end in one call", test_arrival),
    ("H. Studio reads as the tutorials say", test_tutorial_labels),
    ("I. code spans stay on one line", test_split_code_spans),
]


def run():
    """Every check; returns a list of failure strings (empty = pass)."""
    try:
        import PySide6                                    # noqa: F401
        checks = CHECKS
    except Exception:
        print("docs links: PySide6 not installed — Studio check skipped.")
        checks = [c for c in CHECKS
                  if c[1] not in (test_studio_flow, test_arrival,
                                  test_tutorial_labels)]
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
