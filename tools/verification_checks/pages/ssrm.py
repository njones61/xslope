"""Page config: docs/verification/ssrm.md (Griffiths & Lane SSRM examples).

Figure mode: ``structural``.  Every figure here is a single-axes plot emitted by
the FEM sample scripts (inputs, mesh, results, sweep) and its alt text is the
file name, so the caption asserts nothing about layout that could be compared
against the image.  Panel classification would therefore be checking a claim the
page does not make; the check holds each figure to the claim it does make — the
referenced image exists and carries alt text.
"""
from ..config import PageConfig

CONFIG = PageConfig(
    name="ssrm",

    bounds=[
        # an agreement bound over a whole sweep, not over one pair
        ('+2.2', '15 curve within 2.2% at every one of its five stations'),
    ],

    whitelist=[
        # summary-row cell: the L/H station value (0.7) sits between the two
        # factors the delta is measured from, so the pair is named explicitly
        ('+0.8', 'minimum 1.31 vs their FE 1.30', '1.31', '1.30'),
        # Hammah et al. Hoek-Brown row: the derived constant s is printed in
        # scientific notation (2.605e-5 against the paper's 2.5 × 10⁻⁵), which
        # the pairer does not read, so the pair is named explicitly
        ('+4.2', "the paper's 2.5 × 10⁻⁵", '2.605', '2.5'),
    ],

    untagged_allow=[
        # Example 6's reservoir effect is a RATIO of the section's own two
        # locked factors of safety — 1.867 wet over 2.422 dry, both printed in
        # the table above it — against Griffiths & Lane's 1.9 / 2.4.  A ratio
        # is not a factor of safety and the deltas whitelist re-derives
        # percentages and absolute differences, not ratios, so it is named here
        # with both of its operands.
        ('0.77', 'Reservoir effect, wet/dry'),
    ],

    auth_hdr_extra=[r"Griffiths", r"Lane", r"Taylor", r"Bishop & Morgenstern",
                    r"Fig\. 15", r"FE \(Fig"],

    tag_value_keys=['expected_fs*', 'fs_*', 'expected_beta', 'expected_kc',
                    'points', 'expected'],
    # this page quotes SSRM locks at the 2-dp precision of Griffiths & Lane's
    # own figures (expected_fs=1.822 is printed 1.82), so a correctly rounded
    # restatement satisfies the tag
    tag_round_dp=2,
    # No tag on this page needs an exemption. Two coverage locks used to have
    # one — Example 1's tri6/quad8/quad9 element-type sweep and Example 2's
    # coarse tri6 confirmation, both locking 1.36, neither backing a number the
    # page publishes. Example 1's own quad8 result is now 1.36 and is printed,
    # so the checker finds that string and an exemption for it can never fire.
    # A dead exemption is a check failure here, deliberately, so they are gone
    # rather than kept as decoration. If Example 1's published value moves off
    # 1.36 the two coverage locks will report as missing from their section,
    # which is the checker asking for the exemptions back.
    tag_exempt=[],
    # the summary rows quote factors in prose ("Submerged plateau 1.86 vs ...")
    # at the 2-dp precision of the source figure, not the 3-dp lock value, so
    # there is no locked-value pattern that could be swept against the tags
    locked_value_re=None,

    figure_mode="structural",
)
