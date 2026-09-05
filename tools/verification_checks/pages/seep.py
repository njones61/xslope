"""Page config: docs/verification/seep.md (steady-state seepage verification).

Figure mode: ``structural``.  This page's figures are whole-model renders
produced by different sample scripts — a head field with a flow net, an
inputs/solution pair — and their captions describe what is plotted rather than
how many axes it occupies, so panel counting would test nothing the caption
claims.  The check therefore holds each figure to what its caption does assert:
the referenced image exists and the caption is a real description.
"""
from ..config import PageConfig

CONFIG = PageConfig(
    name="seep",

    bounds=[
        # the match-dot legend, not a comparison
        ('+3', 'within 3% of the vendor and/or reference figure'),
        ('+6', '| 🟡 | 3–6% |'),
        ('+6', 'more than 6% |'),
        # tri3-vs-tri6 agreement on the radial problem: the tri3 discharge is
        # not printed, so the tolerance is a bound, not a pair
        ('+0.01', 'tri3 linear elements agree to within'),
        # the independent finite-difference confirmation of the Pavlovsky form
        # factor: the FD discharges themselves are not printed
        ('+0.5', 'agreement at three penetration ratios'),
        # the linear-front problems as a set, against SEEP2D
        ('+0.15', 'agree to better than'),
    ],

    whitelist=[
        # the 3.5–4.7% band is quoted ahead of the two case pairs it summarises;
        # 4.7% is gw009a 2.307e-5 against SEEP2D's 2.421e-5 (3.5% is gw010,
        # 6.07 against 6.29, and reconciles in scope)
        ('+4.7', "reads 3.5–4.7% below SEEP2D", '2.307', '2.421'),
    ],

    # Percentages expressed as a share of a head range or a total drop
    not_a_comparison_extra=[
        r"of (?:the )?total drop",
        r"of (?:the |a )?\d+[\s-]?ft head range",
    ],

    auth_hdr_extra=[r"SEEP2D", r"Exact", r"Thiem", r"Pavlovsky", r"Harr"],

    tag_value_keys=['expected_fs*', 'fs_*', 'expected_flowrate*',
                    'expected_head*', 'expected_beta', 'expected_kc',
                    'points', 'expected'],
    # the summary rows quote discharges and head errors in units, not a fixed
    # 3-dp factor of safety, so there is no locked-value pattern to sweep
    locked_value_re=None,

    figure_mode="structural",
)
