"""Page config: docs/verification/published.md (worked published problems).

Figure mode: ``structural``.  Each entry on this page carries a single
whole-model render of its input file, and the caption says what the model holds
rather than how many axes the image occupies, so panel counting would test a
claim the page does not make.  The check holds each figure to what its caption
does assert: the referenced image exists and the caption describes it.

Tag keys: the pullout rows lock a table rather than a factor of safety, so
``expected_pullout`` (the resistance a layer develops beyond the assumed failure
surface) and ``expected_envelope`` (the full capacity envelope at the same
station) join the standard key list.
"""
from ..config import PageConfig

CONFIG = PageConfig(
    name="published",

    bounds=[
        # the match-dot legend, not a comparison
        ('+3', 'within 3% of the vendor and/or reference figure'),
        ('+6', '| 🟡 | 3–6% |'),
        ('+6', 'more than 6% |'),
        # the summary row states a bound over the whole eleven-layer table; the
        # per-layer pairs it summarises are re-derived row by row in the entry
        ('+0.1', 'reproduced within 0.1% of the manual'),
    ],

    tag_value_keys=['expected_fs*', 'fs_*', 'expected_pullout',
                    'expected_envelope', 'expected_beta', 'expected_kc',
                    'points', 'expected'],
    # the summary rows quote a bound over a table of resistances in lb/ft rather
    # than a fixed 3-dp factor of safety, so there is no locked-value pattern to
    # sweep back against the tags
    locked_value_re=None,

    figure_mode="structural",
)
