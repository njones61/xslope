"""Where Studio sends people: the documentation pages, and the support address.

Every outward-facing address the application shows — Help → Documentation, the
welcome window's links, the address a user writes to — is spelled once, here, so
that re-pointing the documentation site or the support alias is one edit rather
than a search through the GUI.

The addresses are ordinary strings and nothing reads them back: a link opened from
a menu goes to the browser through ``QDesktopServices``, and the welcome window's
links are rich text a ``QLabel`` hands to the browser itself. The one place a URL
is a trust decision is the ``xslope://`` handler, which keeps its own allowlist in
:mod:`studio.urlscheme` — deliberately not read from here, because what a link may
be fetched from is a security decision and not a display constant. The two must
still agree about the documentation host, and ``test/welcome_window_check.py``
checks that they do, along with every page below still existing in ``docs/``.
"""

from __future__ import annotations

__all__ = ["DOCS_URL", "GETTING_STARTED_URL", "SAMPLES_URL", "VERIFICATION_URL",
           "SUPPORT_EMAIL", "SUPPORT_MAILTO"]

#: The documentation root — ``mkdocs.yml``'s ``site_url``. Read The Docs serves the
#: current release at ``/en/latest/``, and every URL below extends this one.
DOCS_URL = "https://xslope.readthedocs.io/en/latest/"

#: Installing Studio, first launch, and where the application's own files land.
GETTING_STARTED_URL = DOCS_URL + "getting_started/install/"

#: The sample problems: a worked example per model, each with its input file. This
#: page IS the sample browser — Studio keeps no second list of its own.
SAMPLES_URL = DOCS_URL + "lem/samples/"

#: Verification and validation: the benchmark and vendor-corpus comparisons.
VERIFICATION_URL = DOCS_URL + "verification/"

#: The address in the product. An alias rather than a person, so the destination
#: can move without a new build going out, and so nobody's personal inbox appears
#: in a screenshot.
SUPPORT_EMAIL = "xslope@jonesgeo.com"
SUPPORT_MAILTO = "mailto:" + SUPPORT_EMAIL
