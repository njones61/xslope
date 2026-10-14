"""Print one release's notes from CHANGELOG.md.

The release workflow puts this text on the GitHub Release page, which is also
where Studio's "Release notes" link in the update notice goes. A version with no
section, or an empty one, is an error, so a release cannot go out with only the
list of commits GitHub generates.

    python packaging/release_notes.py 1.0.3 > notes.md
"""
import argparse
import os
import re
import sys

CHANGELOG = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                         "CHANGELOG.md")

_HEADING = re.compile(r"^## +v?(\S+)")


def section(text, version):
    """The body of the ``## <version>`` section of a changelog, without its heading.

    Raises ValueError when the version has no section or the section is empty.
    """
    version = version.lstrip("v")
    lines, inside = [], False
    for line in text.splitlines():
        m = _HEADING.match(line)
        if m:
            if inside:
                break
            inside = m.group(1) == version
            continue
        if inside:
            lines.append(line)
    body = "\n".join(lines).strip()
    if not body:
        raise ValueError(f"CHANGELOG.md has no notes for {version}: add a "
                         f"'## {version} — <date>' section before tagging.")
    return body


def unwrap(markdown):
    """Join the lines a paragraph or a list item is wrapped across.

    CHANGELOG.md is wrapped for reading in an editor; a GitHub Release page shows
    each line break as a break, so the published notes put each paragraph and
    each item on one line. Headings, blank lines, and new items start fresh.
    """
    out = []
    for line in markdown.splitlines():
        text = line.strip()
        continues = (text and out and out[-1]
                     and not text.startswith(("#", "- ", "* "))
                     and not out[-1].startswith("#"))
        if continues:
            out[-1] += " " + text
        else:
            out.append(line.rstrip())
    return "\n".join(out)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("version", help="the release, e.g. 1.0.3 or v1.0.3")
    ap.add_argument("--changelog", default=CHANGELOG)
    args = ap.parse_args(argv)
    with open(args.changelog, encoding="utf-8") as f:
        text = f.read()
    # The notes are UTF-8 (→, —, σ). A Windows runner's stdout defaults to the
    # cp1252 code page, which cannot encode "→", and the 1.0.3 tag build stopped
    # here on "File → Open".
    sys.stdout.reconfigure(encoding="utf-8")
    try:
        print(unwrap(section(text, args.version)))
    except ValueError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
