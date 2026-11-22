"""The full changelog on the What's New page, read from CHANGELOG.md at build time.

CHANGELOG.md, at the top of the repository, is the one record of what each release
changed. The What's New page carries its own short highlights of the latest release and
a ``<!-- changelog -->`` marker below them; this hook replaces the marker with the
changelog's releases, so the page and the file can never disagree.

The changelog's own title and opening paragraph are dropped (the page has its own), and
every heading is demoted one level, so each release sits under the page's "Full
Changelog" heading. CHANGELOG.md uses only absolute links, which read the same from
the docs site.
"""

import os
import re

PAGE = "whats_new.md"
MARKER = "<!-- changelog -->"


def on_page_markdown(markdown, page, config, files):
    if page.file.src_uri != PAGE or MARKER not in markdown:
        return markdown
    path = os.path.join(os.path.dirname(config.config_file_path), "CHANGELOG.md")
    with open(path, encoding="utf-8") as fh:
        text = fh.read()
    first = re.search(r"^## ", text, flags=re.M)
    if first is None:
        raise ValueError(f"{path}: no release heading ('## ...') to show on {PAGE}")
    body = re.sub(r"^(#{2,5}) ", lambda m: "#" + m.group(1) + " ", text[first.start():], flags=re.M)
    return markdown.replace(MARKER, body)
