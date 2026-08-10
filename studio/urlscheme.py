"""The ``xslope://`` URL scheme — what a web link is allowed to ask Studio to do.

The docs put an **Open in Studio** link beside every sample: ``xslope://open?url=``
followed by the URL-encoded address of that sample's ``.xslz`` package. The OS hands
the link to whichever application registered the scheme, which is Studio when it was
installed by an installer.

**The registration is system-wide, so the link can come from anywhere.** Any page on
the web can contain an ``xslope://`` link, and the handler therefore trusts nothing
in it:

* **One verb.** ``open`` and nothing else. An unknown verb is refused by name rather
  than guessed at, so a future verb can be added without any old build doing
  something surprising with it.
* **An allowlisted host, over https.** :data:`ALLOWED_HOSTS` is the docs site and the
  project's own repository — the places XSLOPE's own packages are published. Any
  other host is refused with a message naming the URL. Redirects are re-checked
  against the same list (:class:`_AllowlistRedirects`), because a redirect from an
  allowlisted host to somewhere else is exactly the way a bare check would be walked
  around.
* **No local paths, ever.** ``file:``, ``data:`` and a bare filesystem path are not
  URLs this scheme accepts. Opening a local file is what File → Open is for; a web
  page must not be able to name one.
* **Nothing happens before the user says so.** The caller confirms what will be
  fetched and from where BEFORE :func:`download_package` is called, so a drive-by
  link can make Studio ask a question and nothing more.

Everything in this module is pure logic and one download; the dialogs, and the
opening itself, live in :mod:`studio.main_window`, which reuses the same
unpack-then-open path a double-clicked package takes.
"""

from __future__ import annotations

import os
import urllib.request
from urllib.parse import parse_qs, quote, unquote, urlparse

from xslope.package import PACKAGE_EXT

__all__ = [
    "SCHEME", "VERBS", "ALLOWED_HOSTS", "TIMEOUT", "MAX_BYTES",
    "SchemeError", "is_scheme_url", "parse_request", "check_url",
    "package_name", "download_package", "build_url",
]

#: The registered scheme. The installers write this string into the OS (macOS
#: ``CFBundleURLTypes``, the Windows registry), and the docs build spells it into
#: every Open in Studio link.
SCHEME = "xslope"

#: Every verb this build understands. One, ruled 2026-08-10: ``open``. Adding a verb
#: later breaks no existing link; accepting an unknown one today would.
VERBS = ("open",)

#: The only hosts a package may be fetched from — the documentation site and the
#: project's repository, which are where XSLOPE publishes packages. This list is the
#: whole of the trust decision: it is not read from the URL, from a setting, or from
#: anything the link can influence.
ALLOWED_HOSTS = (
    "xslope.readthedocs.io",        # the docs site (mkdocs.yml site_url)
    "github.com",                   # the repository (mkdocs.yml repo_url)
    "raw.githubusercontent.com",    # files served straight out of the repository
    "objects.githubusercontent.com",  # where github.com release assets redirect
)

#: Seconds before a fetch gives up. A link that hangs is a link that looks broken.
TIMEOUT = 20.0

#: A ceiling on what one link may make Studio write to disk (200 MB). The largest
#: sample package is a few tens of MB; this is the runaway guard, not a policy.
MAX_BYTES = 200 * 1024 * 1024


class SchemeError(Exception):
    """A request that will not be carried out, phrased for the person who clicked.

    Every refusal in this module raises one, and every message names the thing that
    was refused — the verb, or the URL — because the reader's next question is
    always "what was wrong with it?".
    """


def is_scheme_url(text):
    """True if ``text`` is one of our URLs.

    Used to tell a command-line argument that is a scheme invocation (Windows and
    Linux deliver them as ``argv[1]``) from one that is a file to open.
    """
    return str(text).strip().lower().startswith(SCHEME + ":")


def parse_request(uri):
    """Return ``(verb, url)`` for an ``xslope://`` request, or raise SchemeError.

    ``xslope://open?url=<percent-encoded package URL>``. The verb is the host part
    of the URI (``xslope://open?...``) or its path (``xslope:open?...``) — both
    spellings reach the OS depending on who wrote the link, and they mean the same
    thing.

    The URL is only *extracted* here. What it is allowed to point at is
    :func:`check_url`'s question, kept separate so the two refusals can be tested,
    and read, apart.
    """
    parts = urlparse(str(uri).strip())
    if parts.scheme.lower() != SCHEME:
        raise SchemeError(f"{uri} is not an {SCHEME}: link.")
    verb = (parts.netloc or parts.path).strip("/").lower()
    if not verb:
        raise SchemeError(
            f"{uri} names no action. The only one this version of XSLOPE Studio "
            f"understands is '{VERBS[0]}', as in {SCHEME}://{VERBS[0]}?url=...")
    if verb not in VERBS:
        raise SchemeError(
            f"XSLOPE Studio does not understand '{verb}'. The only action it "
            f"understands is '{VERBS[0]}', as in {SCHEME}://{VERBS[0]}?url=...")
    values = parse_qs(parts.query).get("url") or []
    url = values[0].strip() if values else ""
    if not url:
        raise SchemeError(
            f"{uri} carries no url= parameter, so there is nothing to open.")
    return verb, url


def check_url(url):
    """Return ``url`` if a package may be fetched from it; otherwise raise.

    Two questions, both answered from the URL's own parts and never from its text:
    is the transport https, and is the host one of :data:`ALLOWED_HOSTS`? A host
    that merely *ends with* an allowed name is not allowed — ``xslope.readthedocs.io
    .evil.example`` and ``evilxslope.readthedocs.io`` are other hosts.
    """
    parts = urlparse(url)
    if parts.scheme.lower() != "https":
        raise SchemeError(
            f"XSLOPE Studio will not open {url}: an Open in Studio link may only "
            f"fetch a project package over https, and never a file on this "
            f"computer. Use File → Open for a local file.")
    host = (parts.hostname or "").lower()
    if host not in ALLOWED_HOSTS:
        raise SchemeError(
            f"XSLOPE Studio will not open {url}: {host or 'that address'} is not "
            f"one of the sites XSLOPE publishes projects from "
            f"({', '.join(ALLOWED_HOSTS)}). Download the file yourself and open it "
            f"with File → Open if you trust it.")
    return url


def package_name(url):
    """The file name a package fetched from ``url`` should be saved under.

    Taken from the URL's path, never from a header the server controls, and forced
    to a plain name with the package extension — a name is a place on the user's
    disk, so it may not carry a directory, a ``..``, or a suffix of the server's
    choosing.
    """
    name = os.path.basename(unquote(urlparse(url).path))
    name = name.replace("\\", "/").split("/")[-1].strip()
    stem = os.path.splitext(name)[0].strip(". ")
    if not stem:
        stem = "project"
    return stem + PACKAGE_EXT


class _AllowlistRedirects(urllib.request.HTTPRedirectHandler):
    """Re-check the allowlist at every redirect.

    ``urlopen`` follows redirects on its own, so a check made once on the URL in the
    link would say nothing about where the bytes finally came from. This runs
    :func:`check_url` on each hop instead, which makes an open redirect on an
    allowlisted host a refusal rather than a way through it.
    """

    def redirect_request(self, req, fp, code, msg, headers, newurl):
        check_url(newurl)
        return super().redirect_request(req, fp, code, msg, headers, newurl)


def _opener():
    return urllib.request.build_opener(_AllowlistRedirects())


def download_package(url, dest_dir, progress=None, timeout=TIMEOUT,
                     max_bytes=MAX_BYTES):
    """Fetch the package at ``url`` into ``dest_dir`` and return its path.

    The URL is checked again here, so no caller can reach the network by skipping
    :func:`check_url`. A partial download is deleted rather than left behind to be
    opened as a truncated zip. ``progress(done, total)`` is called as bytes arrive
    (``total`` is 0 when the server sends no length).

    The caller must already have the user's confirmation: by the time this runs, a
    request has been made of a remote host.
    """
    check_url(url)
    os.makedirs(dest_dir, exist_ok=True)
    dest = os.path.join(dest_dir, package_name(url))
    req = urllib.request.Request(url, headers={"User-Agent": "XSLOPE Studio"})
    done = 0
    try:
        with _opener().open(req, timeout=timeout) as resp:
            check_url(resp.geturl())        # the hop the bytes actually came from
            total = int(resp.headers.get("Content-Length") or 0)
            if total and total > max_bytes:
                raise SchemeError(
                    f"{url} is {total / 1e6:.0f} MB, which is larger than a project "
                    f"package should be, so it was not downloaded.")
            with open(dest, "wb") as out:
                while True:
                    block = resp.read(1 << 18)
                    if not block:
                        break
                    done += len(block)
                    if done > max_bytes:
                        raise SchemeError(
                            f"{url} kept sending after "
                            f"{max_bytes / 1e6:.0f} MB, which is larger than a "
                            f"project package should be, so it was discarded.")
                    out.write(block)
                    if progress is not None:
                        progress(done, total)
    except SchemeError:
        _unlink(dest)
        raise
    except Exception as exc:
        _unlink(dest)
        raise SchemeError(f"Could not download {url}: {exc}") from exc
    return dest


def _unlink(path):
    try:
        os.unlink(path)
    except OSError:
        pass


def build_url(package_url, verb=VERBS[0]):
    """The ``xslope://`` link for a package URL — what the docs build emits.

    Kept here so the string has one definition on the reading side and the writing
    side both; ``hooks/docs_packages.py`` builds the same link for the site.
    """
    return f"{SCHEME}://{verb}?url={quote(package_url, safe='')}"
