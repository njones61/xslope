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
  rewritten to the pair **Download** · **Open in Studio**. Both come from one
  string — the package's URL — so they cannot point at different things, and the
  packaging step above is asked afterwards to prove that every package a page
  linked was actually built.

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
import os
import posixpath
import re
import sys
import time
from urllib.parse import urljoin

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
# The repo root goes FIRST: an editable install of xslope elsewhere on the machine
# would otherwise answer this import, and the packages on the site would be built
# by code that is not the code in this checkout.
if sys.path[:1] != [REPO_ROOT]:
    sys.path.insert(0, REPO_ROOT)

from xslope.package import PACKAGE_EXT, pack        # noqa: E402
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
    '<p class="xslz-note"><em>Each input file below is linked as a project '
    'package (<code>.xslz</code>) — the Excel workbook together with the mesh and '
    'results stored beside it. <strong>Download</strong> saves the package; '
    '<strong>Open in Studio</strong> opens it in an installed copy of XSLOPE '
    'Studio, which unpacks it first. Open in Studio needs a Studio installed from '
    'an installer; a pip install registers no <code>xslope://</code> handler, so '
    'use Download.</em></p>\n'
)

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


def _pair_html(inner, pkg_href, pkg_url):
    """The rendered pair: the page's own words, then Download · Open in Studio.

    The anchor's text is kept verbatim and stops being a link. It is often the file
    name, but not always — pages name a workbook mid-sentence ("the SSRM runs on
    vp027_fem.xlsx, which reconstructs…") and label variants by their parameter
    ("$c_u2/c_u1 = 0.8$"). Replacing that text with the word "Download" would make
    those sentences unreadable, so the pair follows the text instead of replacing
    it.
    """
    name = posixpath.basename(pkg_href)
    scheme_url = build_url(pkg_url)
    return (
        f'{inner} (<a class="xslz-download" href="{html.escape(pkg_href, quote=True)}" '
        f'title="Download {html.escape(name, quote=True)} — the workbook and its '
        f'results in one file">Download</a> · '
        f'<a class="xslz-studio" href="{html.escape(scheme_url, quote=True)}" '
        f'title="Open {html.escape(name, quote=True)} in XSLOPE Studio">'
        f'Open in Studio</a>)'
    )


def _insert_note(page_html, first_pair_at):
    """Put the one-time note next to the first pair.

    Next to it means before the paragraph the pair sits in — a note at the top of a
    two-thousand-line corpus page is not near anything. When the first pair is not
    in a paragraph (a list item, a table cell) there is nowhere to put a paragraph
    without breaking the enclosing structure, so the note goes to the top of the
    page instead.
    """
    para = page_html.rfind("<p>", 0, first_pair_at)
    if para != -1 and "</p>" not in page_html[para:first_pair_at]:
        return page_html[:para] + NOTE_HTML + page_html[para:]
    return NOTE_HTML + page_html


def rewrite_links(page_html, page_url, docs_dir, site_url):
    """Turn every sample-workbook link on one page into a pair.

    Returns ``(html, linked)`` where ``linked`` is the set of docs-relative package
    paths the page now points at. ``page_url`` is the page's URL within the site
    (``lem/samples/``), which is what the hrefs in ``page_html`` are relative to —
    MkDocs has already rewritten them by the time a hook sees the rendered content.
    """
    base = posixpath.dirname(page_url)
    linked = set()
    out = []
    last = 0
    first_pair_at = None
    for m in _ANCHOR.finditer(page_html):
        tag = _OPEN_TAG.match(m.group(0)).group(0)
        href_m = _HREF.search(tag)
        if not href_m:
            continue
        href = html.unescape(href_m.group(1))
        if not href.lower().endswith(".xlsx"):
            continue
        cls = _CLASS.search(tag)
        if cls and RAW_CLASS in cls.group(1).split():
            continue                               # the page wants the bare workbook
        docs_rel = posixpath.normpath(posixpath.join(base, href))
        if not is_sample_workbook(docs_rel):
            continue
        if not os.path.isfile(os.path.join(docs_dir, *docs_rel.split("/"))):
            continue                               # a link that was already broken
        inner = m.group(0)[len(tag):-len("</a>")]
        pkg_href = _package_href(href)
        pkg_url = urljoin(urljoin(site_url, page_url), pkg_href)
        out.append(page_html[last:m.start()])
        if first_pair_at is None:
            first_pair_at = sum(len(s) for s in out)
        out.append(_pair_html(inner, pkg_href, pkg_url))
        last = m.end()
        linked.add(_package_href(docs_rel))
    if first_pair_at is None:
        return page_html, linked
    out.append(page_html[last:])
    return _insert_note("".join(out), first_pair_at), linked


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
    new_html, linked = rewrite_links(page_html, page.url, config["docs_dir"], site_url)
    _linked.update(linked)
    return new_html


def on_post_build(config, **kwargs):
    """Write a package for every sample workbook, into the built site."""
    docs_dir = config["docs_dir"]
    site_dir = config["site_dir"]
    started = time.time()
    built = set()
    for docs_rel in sample_workbooks(docs_dir):
        src = os.path.join(docs_dir, *docs_rel.split("/"))
        pkg_rel = _package_href(docs_rel)
        out = os.path.join(site_dir, *pkg_rel.split("/"))
        pack(src, dest=out, overwrite=True)
        built.add(pkg_rel)
    missing = sorted(_linked - built)
    if missing:
        raise DocsPackageError(
            "These pages link project packages that were not built: "
            + ", ".join(missing))
    print(f"docs_packages: {len(built)} project packages "
          f"({len(_linked)} linked) in {time.time() - started:.1f}s")
