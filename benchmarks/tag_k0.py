"""The K0 initial stress a corpus file must DECLARE, read off its own test tag.

A verification model that is scored under a K0 initial stress has to say so in the
workbook. The RS2/Slide2 block is transcribed that way -- RS2 builds every model's
in-situ state from the overburden with sigma_h = K0 * sigma_v, and its slope models
use K0 = 1 -- and until this module existed that convention lived ONLY in the doc
page's ``<!-- test: ... k0=1 ... -->`` tag. The suite passed it as a solver keyword,
so the locked number was reproducible from the tag and from nothing else: the same
file opened in Studio, or loaded with ``load_slope_data`` and solved the way a user
would, initialized by gravity turn-on instead and answered a different question.
107 of the 108 tagged files were in that state.

The fix is that the FILE carries the convention (``main!D16``, loaded as
``slope_data['k0']``), and the tag keeps it only as a redundant restatement. To make
the two incapable of drifting, the builders do not hold a table of their own -- they
read the tag. This module is that reader:

    apply_tag_k0(slope_data, out_path)     # before writing out_path

sets ``slope_data['k0']`` to the value the docs tag for that file declares, and to
``None`` when no tag declares one. The clearing half is not optional. Nearly every
corpus builder starts from ``docs/lem/files/xslope_acads_simple.xlsx`` as a geometry
donor, and that file is RS2-1 and now carries K0 = 1; without the clear its K0 would
ride into every unrelated problem derived from it, exactly the way the donor's
tensile cap and elastic constants once did (see ``vendor_tcut.apply_vendor_t_cut``).

Keyed by file BASENAME, like ``VENDOR_T_CUT``, because ``benchmarks/verify_rebuild.py``
rebuilds each target into a scratch directory: a full-path key would read as
"untagged" there and the rebuild guard would report a spurious K0 diff on all 108.

``run_tests.py``'s ``tag_k0`` row is the other half of the contract: it fails if any
tagged file's ``main!D16`` disagrees with the k0 its tag names, so neither a builder
that forgets this call nor a hand-edited tag can ship the mismatch this module exists
to remove.
"""

import glob
import os
import re

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_TAG_RE = re.compile(r'<!--\s*test:\s*(.*?)-->', re.S)

_MAP = None


def tag_k0_map():
    """{file basename: K0} over every ``k0=`` test tag in ``docs/**/*.md``.

    Raises if two tags name the same file with different values -- that is a
    corpus contradiction, and silently picking one of them is how a locked
    number ends up unreproducible.
    """
    global _MAP
    if _MAP is not None:
        return _MAP
    found = {}
    for md in sorted(glob.glob(os.path.join(_ROOT, 'docs', '**', '*.md'),
                               recursive=True)):
        with open(md, encoding='utf-8') as fh:
            text = fh.read()
        for match in _TAG_RE.finditer(text):
            params = {}
            for pair in match.group(1).split(','):
                key, sep, value = pair.partition('=')
                if sep:
                    params[key.strip()] = value.strip()
            if 'k0' not in params or 'file' not in params:
                continue
            name = os.path.basename(params['file'])
            value = float(params['k0'])
            if name in found and found[name] != value:
                raise ValueError(
                    f"conflicting k0 test tags for {name}: "
                    f"{found[name]:g} and {value:g}")
            found[name] = value
    _MAP = found
    return _MAP


def tag_k0(path):
    """The K0 declared for ``path`` by the docs test tags, or None."""
    return tag_k0_map().get(os.path.basename(str(path)))


def apply_tag_k0(slope_data, path):
    """Set ``slope_data['k0']`` for the file about to be written to ``path``.

    Always assigns -- to the tag's value, or to None when nothing declares one --
    so a K0 inherited from a loaded donor file cannot survive into a problem that
    does not claim it. Returns the value written.
    """
    value = tag_k0(path)
    slope_data['k0'] = value
    return value
