"""Standing check on CHANGELOG.md and the script that lifts one release's notes out
of it for the GitHub Release page (packaging/release_notes.py).

The release build publishes the section for the tag and stops if there is none;
this check catches the same gaps between releases: a released version with no
notes, a heading the script would not find, a section that runs into the next
one, or headings other than New, Fixed, and Changed. Pure Python, no network; the
released-versions leg reads the git tags and is skipped where there is no git.
"""
import os
import re
import subprocess
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(_HERE)
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)

CHANGELOG = os.path.join(_REPO, "CHANGELOG.md")
#: The first release the changelog covers; earlier ones are on the releases page.
FIRST = (1, 0, 0)
SUBHEADINGS = {"New", "Fixed", "Changed"}


def _released_versions():
    """Release tags from FIRST on, as 'X.Y.Z' strings, or None without git."""
    try:
        out = subprocess.run(["git", "-C", _REPO, "tag", "--list", "v*"],
                             capture_output=True, text=True, check=True).stdout
    except (OSError, subprocess.CalledProcessError):
        return None
    versions = []
    for tag in out.split():
        m = re.fullmatch(r"v(\d+)\.(\d+)\.(\d+)", tag)
        if m and tuple(int(g) for g in m.groups()) >= FIRST:
            versions.append(tag[1:])
    return versions


def _load_script():
    """packaging/release_notes.py, loaded by path: the folder is not a package, and
    its name is taken by the installed ``packaging`` library."""
    import importlib.util
    path = os.path.join(_REPO, "packaging", "release_notes.py")
    spec = importlib.util.spec_from_file_location("xslope_release_notes", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def run():
    script = _load_script()
    cli, section = script.main, script.section

    fails = []
    with open(CHANGELOG, encoding="utf-8") as f:
        text = f.read()

    # 1. Every released version has notes, and there is a place for the next one.
    released = _released_versions()
    if released is None:
        print("  no git here; released-versions leg skipped")
    else:
        for v in released:
            try:
                section(text, v)
            except ValueError:
                fails.append(f"released version {v} has no notes in CHANGELOG.md")
    if not re.search(r"^## Unreleased\s*$", text, re.M):
        fails.append("CHANGELOG.md has no '## Unreleased' section for changes "
                     "that have landed since the last release")

    # 2. Inside a section, only the three kinds of entry.
    for heading in re.findall(r"^### +(.+?)\s*$", text, re.M):
        if heading not in SUBHEADINGS:
            fails.append(f"'### {heading}' is not one of New, Fixed, or Changed")

    # 3. The script takes one section, whole, and nothing of the next; a 'v' on
    #    the tag is accepted; a version with no section, or an empty one, is an
    #    error, both in the function and as the script's exit status.
    sample = ("# Changelog\n\n## Unreleased\n\n## 2.0.0 — 2027-01-01\n\n"
              "### New\n\n- one\n\n### Fixed\n\n- two\n\n"
              "## 1.9.0 — 2026-12-01\n\n### Fixed\n\n- three\n")
    got = section(sample, "v2.0.0")
    if got != "### New\n\n- one\n\n### Fixed\n\n- two":
        fails.append(f"the 2.0.0 section reads {got!r}")
    for missing in ("1.8.0", "Unreleased"):
        try:
            section(sample, missing)
            fails.append(f"{missing!r} gave notes instead of an error")
        except ValueError:
            pass
    # 4. On the release page each item is one line: GitHub shows a line break
    #    inside an item as a break, so the wrapped lines of CHANGELOG.md are joined.
    wrapped = ("### Fixed\n\n- a fix that is\n  wrapped across lines\n- another\n\n"
               "A paragraph\nin two lines.")
    want = ("### Fixed\n\n- a fix that is wrapped across lines\n- another\n\n"
            "A paragraph in two lines.")
    if script.unwrap(wrapped) != want:
        fails.append(f"the published notes read {script.unwrap(wrapped)!r}")

    devnull = open(os.devnull, "w")
    real_err, sys.stderr = sys.stderr, devnull
    try:
        status = cli(["0.0.1", "--changelog", CHANGELOG])
    finally:
        sys.stderr = real_err
        devnull.close()
    if status != 1:
        fails.append(f"the script exits {status} for a version with no notes, not 1")
    return fails


def main():
    print("Release notes (CHANGELOG.md):")
    failures = run()
    if failures:
        print("\nFAILURES:")
        for f in failures:
            print(f"  - {f}")
        raise SystemExit(1)
    print("  every release from 1.0.0 has notes, and the script lifts one section "
          "whole")
    print("\nAll release-notes checks passed.")


if __name__ == "__main__":
    main()
