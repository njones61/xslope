"""Sample projects on the docs site: one package per workbook, one pair of links.

A sample problem is a *set* of files — the workbook plus the mesh and result
sidecars written beside it — so the docs hand out the set, not the workbook alone.
This hook does both halves of that at build time:

* **Packaging** (``on_post_build``): every ``docs/*/files/**/*.xlsx`` is collected
  into a ``.xslz`` project package written to the matching place in the built site.
  The packing is ``xslope.package.pack`` — the same code Studio's File → Export
  Project Package runs, so a sample that gains a sidecar gains it in the package on
  the next build with nothing to edit. Nothing is written into ``docs/``; no
  package is committed to the repository.

* **The paired link** (``on_page_content``): a link to a sample workbook is
  rewritten to the pair **Download** · **Open in Studio**. Open in Studio always
  names the package — the deep link unpacks one, and a package holding a lone
  workbook is a perfectly good package. Download names the package too, *except*
  when there is nothing beside the workbook to package: a sample whose project is
  the ``.xlsx`` and nothing else is downloaded as that ``.xlsx``, because a page
  that tells the reader to open the file in Excel has to hand them a file Excel
  opens. Either way the package is built and registered, and the packaging step
  above is asked afterwards to prove that every package a page linked was actually
  built.

``Open in Studio`` is ``xslope://open?url=<package URL>``: the OS hands it to
whatever registered the ``xslope`` scheme (Studio, when installed by an installer),
which confirms, downloads and opens it. A pip install registers no scheme, so the
Download link beside it is what those readers use — the note this hook drops on
each page says so once.

**Escape hatch.** An anchor carrying the class ``raw-file`` is left exactly as
written, workbook link and all. Authored with ``attr_list``::

    [input_template.xlsx](files/input_template.xlsx){: .raw-file }

for a page that deliberately offers the bare workbook rather than the project.
"""

import html
import logging
import os
import posixpath
import re
import shutil
import sys
import time
from urllib.parse import unquote, urljoin

log = logging.getLogger("mkdocs.hooks.docs_packages")

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
# The repo root goes FIRST: an editable install of xslope elsewhere on the machine
# would otherwise answer this import, and the packages on the site would be built
# by code that is not the code in this checkout.
if sys.path[:1] != [REPO_ROOT]:
    sys.path.insert(0, REPO_ROOT)

from xslope.package import PACKAGE_EXT, pack, project_files    # noqa: E402
# The link the site emits is built by the module that parses it, so the two cannot
# drift apart. studio.urlscheme is pure stdlib + xslope.package: importing it here
# pulls in no PySide6 and nothing a docs build does not already have.
from studio.urlscheme import build_url              # noqa: E402

#: Where a sample lives. Every workbook under a ``files/`` directory in the docs
#: tree is a sample project, and the same predicate decides both what gets packaged
#: and what gets a paired link — so a link can never name a package nobody built.
FILES_DIR = "files"

#: Anchors marked with this class are left alone (see the module docstring).
RAW_CLASS = "raw-file"

#: Said once per page, next to the first pair.
NOTE_HTML = (
    '<p class="xslz-note"><em>Each input file below is linked twice. '
    '<strong>Download</strong> saves the sample\'s project package '
    '(<code>.xslz</code>) — the Excel workbook together with the mesh and results '
    'stored beside it — or the workbook itself, for a sample that has nothing '
    'stored beside it. <strong>Open in Studio</strong> opens the package in an '
    'installed copy of XSLOPE Studio, which unpacks it first. Open in Studio needs '
    'a Studio installed from an installer; a pip install registers no '
    '<code>xslope://</code> handler, so use Download.</em></p>\n'
)

_TAG = re.compile(r"<(/?)([a-zA-Z][a-zA-Z0-9]*)\b[^>]*?(/?)>")

#: Tags that never enclose anything, so they never open a block the note has to
#: clear. The HTML the hook sees is rendered Markdown, not arbitrary input.
_VOID_TAGS = frozenset("area base br col embed hr img input link meta param "
                       "source track wbr".split())

_ANCHOR = re.compile(r"<a\b[^>]*>.*?</a>", re.DOTALL | re.IGNORECASE)
_HREF = re.compile(r'href\s*=\s*"([^"]*)"', re.IGNORECASE)
_CLASS = re.compile(r'class\s*=\s*"([^"]*)"', re.IGNORECASE)
_OPEN_TAG = re.compile(r"<a\b[^>]*>", re.DOTALL | re.IGNORECASE)

#: Every package a page linked, docs-relative (``lem/files/x.xslz``). Checked
#: against what was built, in ``on_post_build``.
_linked = set()


class DocsPackageError(Exception):
    """A build-stopping problem: a link the site could not honour."""


# ---------------------------------------------------------------------------
# What is a sample workbook
# ---------------------------------------------------------------------------

def is_sample_workbook(docs_rel):
    """True if ``docs_rel`` (a docs-relative posix path) names a sample workbook.

    The rule is positional, not a list: a ``.xlsx`` somewhere under a ``files/``
    directory. ``docs/inputs/input_template.xlsx`` is therefore not one — it is a
    blank template, not a project — and neither is a workbook quoted in a code
    block, which is not a link at all.
    """
    parts = docs_rel.split("/")
    name = parts[-1]
    return (name.lower().endswith(".xlsx") and not name.startswith("~$")
            and FILES_DIR in parts[:-1])


def sample_workbooks(docs_dir):
    """Every sample workbook in the docs tree, docs-relative, in sorted order."""
    found = []
    for root, dirs, names in os.walk(docs_dir):
        dirs.sort()
        rel = os.path.relpath(root, docs_dir).replace(os.sep, "/")
        rel = "" if rel == "." else rel
        for name in sorted(names):
            docs_rel = f"{rel}/{name}" if rel else name
            if is_sample_workbook(docs_rel):
                found.append(docs_rel)
    return found


# ---------------------------------------------------------------------------
# The paired link
# ---------------------------------------------------------------------------

def _package_href(href):
    """The package's href, from the workbook's. One substitution, so the two links
    the pair emits are the same string by construction."""
    return href[: -len(".xlsx")] + PACKAGE_EXT


def _workbook_href(href):
    """The workbook's href, from the package's — what a page should have linked."""
    return href[: -len(PACKAGE_EXT)] + ".xlsx"


def _pair_html(inner, dl_href, pkg_href, pkg_url):
    """The rendered pair: the page's own words, then Download · Open in Studio.

    The anchor's text is kept verbatim and stops being a link. It is often the file
    name, but not always — pages name a workbook mid-sentence ("the SSRM runs on
    vp027_fem.xlsx, which reconstructs…") and label variants by their parameter
    ("$c_u2/c_u1 = 0.8$"). Replacing that text with the word "Download" would make
    those sentences unreadable, so the pair follows the text instead of replacing
    it.

    ``dl_href`` is what Download hands over — the package for a sample that has
    sidecars, the workbook itself for one that does not — and ``pkg_href`` is always
    the package, which is what Open in Studio carries. The two are the same string
    for every sample with something stored beside its workbook; they differ only
    where a package would have held the workbook alone, and the title says which
    of the two the reader is about to get.
    """
    name = posixpath.basename(dl_href)
    what = ("the workbook" if dl_href.lower().endswith(".xlsx")
            else "the workbook and its results in one file")
    scheme_url = build_url(pkg_url)
    return (
        f'{inner} (<a class="xslz-download" href="{html.escape(dl_href, quote=True)}" '
        f'title="Download {html.escape(name, quote=True)} — {what}">Download</a> · '
        f'<a class="xslz-studio" href="{html.escape(scheme_url, quote=True)}" '
        f'title="Open {html.escape(posixpath.basename(pkg_href), quote=True)} in '
        f'XSLOPE Studio">Open in Studio</a>)'
    )


def _enclosing_block_end(page_html, at):
    """Where the outermost element still open at ``at`` closes, or ``None``.

    A note is a paragraph, and a paragraph cannot be written inside a table row or
    a list item. The nearest offset at which one *can* be written is just past the
    close of the outermost element the pair is nested in — the ``<table>`` around
    the ``<td>``, the ``<ul>`` around the ``<li>``.

    Counted with a tag stack rather than matched with a pattern: the pair's own
    cell is several elements deep (``table > tbody > tr > td``), and the depth
    varies with what the author wrote.
    """
    stack = []
    for m in _TAG.finditer(page_html, 0, at):
        closing, name, self_closing = m.group(1), m.group(2).lower(), m.group(3)
        if name in _VOID_TAGS or self_closing:
            continue
        if not closing:
            stack.append(name)
            continue
        for i in range(len(stack) - 1, -1, -1):    # unwind to the tag being closed,
            if stack[i] == name:                   # forgiving of anything left open
                del stack[i:]
                break
    if not stack:
        return None
    outer = stack[0]
    depth = stack.count(outer)                     # nested copies close first
    for m in _TAG.finditer(page_html, at):
        if m.group(2).lower() != outer or m.group(3):
            continue
        if not m.group(1):
            depth += 1
        else:
            depth -= 1
            if depth == 0:
                return m.end()
    return None


def _insert_note(page_html, first_pair_at):
    """Put the one-time note next to the first pair.

    Next to it means before the paragraph the pair sits in — a note at the top of a
    two-thousand-line corpus page is not near anything.

    When the first pair is *not* in a paragraph (a list item, a table cell) there is
    nowhere above it to put a paragraph without breaking the enclosing structure, so
    the note goes immediately after that structure closes. On a tutorial page the
    pair is the Completed-model row of the header table, and the note lands under
    the table: title, then the problem sketch, then the objectives, then the header
    table, then the note. A reader decides whether this is the tutorial they want by
    looking at the picture, so nothing this hook adds may come between the two —
    which is what putting the note at the top of the page, or under the H1, did.
    """
    para = page_html.rfind("<p>", 0, first_pair_at)
    if para != -1 and "</p>" not in page_html[para:first_pair_at]:
        return page_html[:para] + NOTE_HTML + page_html[para:]
    end = _enclosing_block_end(page_html, first_pair_at)
    if end is not None:
        return page_html[:end] + "\n" + NOTE_HTML + page_html[end:]
    return NOTE_HTML + page_html


def rewrite_links(page_html, page_url, docs_dir, site_url):
    """Turn every sample-workbook link on one page into a pair.

    Returns ``(html, linked, warnings)``: the rewritten page, the set of docs-relative
    package paths it now points at, and anything the caller must say out loud.
    ``page_url`` is the page's URL within the site (``lem/samples/``), which is what
    the hrefs in ``page_html`` are relative to — MkDocs has already rewritten them by
    the time a hook sees the rendered content.

    Every sample paired here is registered as a linked package, whatever its Download
    link names: Open in Studio carries the package for all of them, and a package the
    build did not write is a broken link either way.

    A link that LOOKS like a sample link and is not rewritten is reported rather than
    passed over. Silence there is the expensive failure: the page keeps a bare ``.xlsx``
    link that nobody notices has stopped being a pair, and the reason is always
    something small — a name that had to be percent-encoded, a fragment on the end, a
    file that has moved. The same goes for a page that links a ``.xslz`` itself: that
    link is registered here, or the build's own does-it-exist check never sees it.
    """
    base = posixpath.dirname(page_url)
    linked = set()
    warnings = []
    out = []
    last = 0
    first_pair_at = None
    for m in _ANCHOR.finditer(page_html):
        tag = _OPEN_TAG.match(m.group(0)).group(0)
        href_m = _HREF.search(tag)
        if not href_m:
            continue
        href = html.unescape(href_m.group(1))
        cls = _CLASS.search(tag)
        raw = bool(cls and RAW_CLASS in cls.group(1).split())
        bare = href.split("#")[0].split("?")[0]
        decoded = unquote(bare)
        if decoded.lower().endswith(PACKAGE_EXT):
            # A page linking a package directly, by hand. It is still a link this
            # build has to honour, so it goes on the list the build checks.
            pkg_rel = posixpath.normpath(posixpath.join(base, decoded))
            linked.add(pkg_rel)
            warnings.append(
                f"{pkg_rel} is linked as a package by hand; the pair is emitted from "
                f"a link to the workbook, so link {_workbook_href(pkg_rel)} instead")
            continue
        if not decoded.lower().endswith(".xlsx"):
            continue
        if raw:
            continue                               # the page wants the bare workbook
        docs_rel = posixpath.normpath(posixpath.join(base, decoded))
        if not is_sample_workbook(docs_rel):
            continue
        if decoded != href:
            # The href needs decoding, or carries a query or a fragment. The pair is
            # built by substituting the extension in the href, which is exact for a
            # plain relative path and guesswork for anything else, so it is not
            # attempted — and not attempted quietly is how a page loses its pair.
            warnings.append(
                f"{href} was left as a bare workbook link: an href with a fragment, "
                f"a query or percent-encoding is not turned into a pair (rename the "
                f"file, or mark the link {{: .{RAW_CLASS} }})")
            continue
        src = os.path.join(docs_dir, *docs_rel.split("/"))
        if not os.path.isfile(src):
            warnings.append(
                f"{href} was left as a bare workbook link: there is no "
                f"{docs_rel} in the docs tree to package")
            continue
        inner = m.group(0)[len(tag):-len("</a>")]
        pkg_href = _package_href(href)
        pkg_url = urljoin(urljoin(site_url, page_url), pkg_href)
        # What Download hands over. A package holding the workbook and nothing else
        # is a package the reader has to unpack to get at the file the page just
        # told them to open in Excel, so that sample is downloaded as its workbook.
        # The classification is the packer's own, so a sample that gains a sidecar
        # goes back to a package Download with nothing to edit here.
        dl_href = href if len(project_files(src)) == 1 else pkg_href
        out.append(page_html[last:m.start()])
        if first_pair_at is None:
            first_pair_at = sum(len(s) for s in out)
        out.append(_pair_html(inner, dl_href, pkg_href, pkg_url))
        last = m.end()
        linked.add(_package_href(docs_rel))
    if first_pair_at is None:
        return page_html, linked, warnings
    out.append(page_html[last:])
    return _insert_note("".join(out), first_pair_at), linked, warnings


# ---------------------------------------------------------------------------
# MkDocs hook entry points
# ---------------------------------------------------------------------------

def on_pre_build(config, **kwargs):
    _linked.clear()


def on_page_content(page_html, page, config, files, **kwargs):
    site_url = config.get("site_url")
    if not site_url:
        raise DocsPackageError(
            "site_url is not set in mkdocs.yml, so the Open in Studio links have no "
            "absolute package URL to carry. Set site_url, or drop this hook.")
    new_html, linked, warnings = rewrite_links(
        page_html, page.url, config["docs_dir"], site_url)
    _linked.update(linked)
    for message in warnings:
        # A MkDocs warning, so `mkdocs build --strict` fails on it rather than the
        # page quietly losing its pair.
        log.warning("docs_packages: %s: %s", page.file.src_uri, message)
    return new_html


def cache_dir(site_dir):
    """Where built packages are kept between builds — beside the site directory.

    MkDocs empties the site directory at the start of every build, ``mkdocs serve``
    included, and packing four hundred projects takes some twenty seconds. Kept
    outside the site, the packages survive that and a rebuild costs a link. The
    location follows the site rather than being fixed, so a docs tree built
    somewhere else (the test's scratch tree) gets its own cache and cannot collide
    with this one.
    """
    return os.path.join(os.path.dirname(os.path.abspath(site_dir)), ".xslz_cache")


def _is_current(package, files):
    """True if ``package`` was written after every file that belongs in it.

    The one staleness rule: a sample that gains a sidecar, or whose workbook is
    edited, is older than nothing and is packed again. Nothing about a package is
    hand-maintained, and no cache entry outlives what it was made from.
    """
    try:
        stamp = os.path.getmtime(package)
    except OSError:
        return False
    return all(os.path.getmtime(f) <= stamp for f in files)


def _place(cached, out):
    """Put the cached package into the site: a hard link where the filesystem
    allows one (no bytes copied, and the site's copy is the cache's), a plain copy
    where it does not."""
    os.makedirs(os.path.dirname(out), exist_ok=True)
    if os.path.exists(out):
        os.unlink(out)
    try:
        os.link(cached, out)
    except OSError:
        shutil.copy2(cached, out)


def on_post_build(config, **kwargs):
    """Write a package for every sample workbook, into the built site."""
    docs_dir = config["docs_dir"]
    site_dir = config["site_dir"]
    cache = cache_dir(site_dir)
    started = time.time()
    built, packed = set(), 0
    for docs_rel in sample_workbooks(docs_dir):
        src = os.path.join(docs_dir, *docs_rel.split("/"))
        pkg_rel = _package_href(docs_rel)
        cached = os.path.join(cache, *pkg_rel.split("/"))
        built.add(pkg_rel)
        if not _is_current(cached, project_files(src)):
            pack(src, dest=cached, overwrite=True)
            packed += 1
        _place(cached, os.path.join(site_dir, *pkg_rel.split("/")))
    missing = sorted(_linked - built)
    if missing:
        raise DocsPackageError(
            "These pages link project packages that were not built: "
            + ", ".join(missing))
    print(f"docs_packages: {len(built)} project packages, {packed} rebuilt "
          f"({len(_linked)} linked) in {time.time() - started:.1f}s")
