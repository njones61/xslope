"""Per-page configs for the verification-page checkers.

One module per page under docs/verification.  ``PAGES`` maps the page stem to
its config; ``config_for`` accepts either a stem or a path.
"""
import os

from . import (geostudio, published, rocscience,  # noqa: F401
               rocscience_groundwater, rs2, seep, ssrm)

PAGES = {
    "rs2": rs2.CONFIG,
    "rocscience": rocscience.CONFIG,
    "geostudio": geostudio.CONFIG,
    "rocscience_groundwater": rocscience_groundwater.CONFIG,
    "ssrm": ssrm.CONFIG,
    "seep": seep.CONFIG,
    "published": published.CONFIG,
}

#: Page order used by the runner and the manifest.
ORDER = ["rs2", "rocscience", "geostudio", "rocscience_groundwater",
         "ssrm", "seep", "published"]


def config_for(page):
    """Config for a page stem or a path to the page's .md file."""
    stem = os.path.splitext(os.path.basename(page))[0]
    if stem not in PAGES:
        raise KeyError(f"no verification-page config named {stem!r}; "
                       f"known pages: {', '.join(ORDER)}")
    return PAGES[stem]
