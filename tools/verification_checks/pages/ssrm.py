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
        # the match-dot legend, not a comparison
        ('+3', 'within 3% of the vendor and/or reference figure'),
        ('+6', '| 🟡 | 3–6% |'),
        ('+6', 'more than 6% |'),
        # the quad8 spread across the three refined Example 5 locks
        # (−1.6% / −1.5% / −3.6% against Griffiths & Lane's Fig. 15 points)
        ('+3.6', 'the three refined quad8 locks read'),
        ('+3.6', "reading 1.5–3.6% below the paper's FE points there"),
        ('+3.6', 'The three refined quad8 locks read'),
        # agreement bounds over a whole sweep, not over one pair
        ('+1', 'FE point for that case (1.45) to within'),
        ('+4', "the paper's wedge solution for the governing mechanism to within"),
        ('+1', 'reproduces the relative jump between them to the same'),
        ('+1', '15 curve within 1% at every one of the eight stations'),
        ('+1', "within 1% of the paper's own plotted FE value at every station"),
        ('+1', 'the location the paper states and within'),
    ],

    whitelist=[
        # summary-row cell: the L/H station value (0.7) sits between the two
        # factors the delta is measured from, so the pair is named explicitly
        ('+0.8', 'minimum 1.31 vs their FE 1.30', '1.31', '1.30'),
    ],

    auth_hdr_extra=[r"Griffiths", r"Lane", r"Taylor", r"Bishop & Morgenstern",
                    r"Fig\. 15", r"FE \(Fig"],

    tag_value_keys=['expected_fs*', 'fs_*', 'expected_beta', 'expected_kc'],
    # this page quotes SSRM locks at the 2-dp precision of Griffiths & Lane's
    # own figures (expected_fs=1.822 is printed 1.82), so a correctly rounded
    # restatement satisfies the tag
    tag_round_dp=2,
    # Coverage locks that exercise a code path rather than back a published
    # number, so the page prints no value for them: the tri6/quad8/quad9
    # element-type sweep on Example 1, and the coarse tri6 confirmation on
    # Example 2 that the foundation layer leaves the toe-failure FS unchanged.
    tag_exempt=[
        ('1.36', 'type=fem_elements'),
        ('1.36', 'xslope_griffiths2.xlsx, type=fem_ssrm, expected_fs=1.36'),
    ],
    # the summary rows quote factors in prose ("Submerged plateau 1.86 vs ...")
    # at the 2-dp precision of the source figure, not the 3-dp lock value, so
    # there is no locked-value pattern that could be swept against the tags
    locked_value_re=None,

    figure_mode="structural",
)
