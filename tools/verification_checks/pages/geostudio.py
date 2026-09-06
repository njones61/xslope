"""Page config: docs/verification/geostudio.md (GeoStudio SLOPE/W + SEEP/W).

Figure mode: ``panel``.  The corpus figures are the same two-panel composite the
Slide2 page uses — inputs on the left, the representative solution on the right —
and the caption says so, which is a layout claim the classifier can read back off
the PNG.  The five SEEP/W transient figures at the end are time-series plots
whose captions name what is plotted against what and claim no layout, so each is
named in ``caption_exempt`` rather than left in a silent "other" bucket.
"""
from ..config import PageConfig

CONFIG = PageConfig(
    name="geostudio",

    bounds=[
        # the match-dot legend, not a comparison
        # agreement bounds over a whole set of surfaces or cases, where no
        # single pair exists to re-derive
        ('+1', "Taylor σ_F within ≈ 1% of SLOPE/W's Monte Carlo"),
        ('+0.19', "Morgenstern-Price within 0.19% of SLOPE/W's own M-P"),
        ('+0.5', 'Bishop within 0.5% of Slide2 Table'),
        ('+0.17', 'Bishop within 0.5% of Slide2 Table'),
        ('+0.17', 'per-slice table to'),
        ('+1.5', 'biases the factor of safety roughly'),
        ('+372', 'too much weight'),
        ('+0.01', 'to within 0.01%, and the Baker'),
        ('+0.004', "Total sliding weight matches SLOPE/W's to"),
        ('+10', "above the paper's 0.155"),
        # the reinforced-fill residuals: quoted at the precision of the
        # measurement in prose, at one decimal in the table.  Neither operand
        # is published — the base is SLOPE/W's own reinforced factor of safety
        # on its own critical circle, which §2.24's column header names — so
        # there is no pair to re-derive and the claim is a bound.
        ('−0.27', "reproduces its factor of safety to −0.27% (clay)"),
        ('−0.64', "reproduces its factor of safety to −0.27% (clay)"),
        ('−0.27', "reproduces SLOPE/W's reinforced factor of safety to within"),
        ('−0.64', "reproduces SLOPE/W's reinforced factor of safety to within"),
        ('−0.3', 'Clay fill, reinforced (imported'),
        ('−0.6', 'Sand fill, reinforced (imported'),
        # T04's worst mid-fill offset, a share of the 6.5 m pond head
        ('+2.1', '2.1% at worst'),
    ],

    whitelist=[
        # Hassan & Wolff G/H transposition: the two reversed-row residuals,
        # each against the paper's own printed value for that surface (the two
        # un-reversed residuals reconcile from the same paragraph unaided)
        ('−3.7', "moves XSLOPE's residuals against them from", '6.059', '6.293'),
        ('+9.3', "moves XSLOPE's residuals against them from", '11.561', '10.576'),
        # vp063 critical k_c: the Spencer / Bishop pair shares one cell, so the
        # delta against Slide2's Spencer needs its two operands named
        ('+10.6', '| Critical k꜀ (Spencer / Bishop) |', '0.167', '0.151'),
        # gs2_45: the Note cell and the paragraph below restate the row's −0.1%
        # to two decimals, from the same two factors the row prints
        ('−0.09', 'the two agree to −0.09%', '1.173', '1.174'),
        # the two orientation-dependent strength formulations, one per row
        ('+0.4', 'the two formulations differ by', '1.118', '1.113'),
    ],

    # Percentages expressed as a share of a reference quantity
    not_a_comparison_extra=[
        r"consolidation",
        r"of the 10 kPa initial excess",
        r"of the 8 m suction step",
        r"of the 8 m drawdown",
        r"of the 8 m column",
        r"of the 6\.5 m pond head",
        r"of the interslice forces",
    ],
    # Quantities identified by what PRECEDES them: a probability of failure and
    # a reported degree of consolidation are not comparisons.
    not_a_comparison_prefix=[
        r"\bPF (?:\d+(?:\.\d+)? → )?",
        r"degree of\s+consolidation is \d+/\d+/",
        r"Terzaghi's \d+/\d+/",
    ],
    # the Terzaghi consolidation table's U column holds degrees of
    # consolidation, not comparisons
    share_hdr_extra=[r"U \(Terzaghi\)"],

    auth_hdr_extra=[r"SLOPE/W", r"SEEP/W", r"Slide\b", r"Book \(", r"Loukidis",
                    r"Baker & Leshchinsky", r"Terzaghi"],

    tag_value_keys=['expected_fs*', 'fs_*', 'expected_beta', 'expected_kc',
                    'expected_flowrate*', 'expected_head*', 'points',
                    'expected'],
    tag_round_dp=2,

    # How much of each probe set the section publishes, where it does not
    # publish all of it.  Every other list lock on this page is printed in
    # full and machine-checked element by element — the two rapid-drawdown
    # factor-of-safety curves (11 steps each), T02's and T04's and T05's and
    # T07's head tables — so a wrong digit in a printed cell fails.
    tag_list_published=[
        # T01 tabulates the excess pore pressure γ_w(h − 100) in kPa at the
        # column center against Terzaghi Eq 17.3; the probes lock the
        # datum-offset total head, which the section never prints.
        ('benchmark=SEEPW-CONS-t', 0),
        # T03's head table publishes a selection of the four locked stations
        # per state: three of the initial condition — (35, 2) appears only on
        # the instantaneous end-state row —
        ('benchmark=SEEPW-RDD-inst-ic', 3),
        # two of the instantaneous end state, (20, 5) and (35, 2),
        ('benchmark=SEEPW-RDD-inst-end', 2),
        # two of the slow mid-drawdown frame, (20, 5) and (30, 3),
        ('benchmark=SEEPW-RDD-slow-mid', 2),
        # and none of the slow end state, which the section reports as the
        # same drained state the instantaneous case ends at.
        ('benchmark=SEEPW-RDD-slow-end', 0),
    ],

    locked_value_re=None,

    figure_mode="panel",
    caption_rules=[
        ("inputs and", "two"),
    ],
    # SEEP/W transient comparisons: the caption names the series plotted against
    # SEEP/W's (or, for the T03 stability curve, SLOPE/W's) and claims no panel
    # layout.
    caption_exempt={
        # The transient problems' inputs figures: a section with its seepage
        # boundary conditions and a legend below, which the quadrant ink test
        # reads as a two-panel gutter; the caption claims no panel layout.
        "images/gs2_cons_inputs.png",
        "images/gs2_infil_inputs.png",
        "images/gs2_rdd_inst_inputs.png",
        "images/gs2_rdd_slow_inputs.png",
        "images/gs2_pond_inputs.png",
        "images/gs2_heap_inputs.png",
        "images/gs2_mso_inputs.png",
        "images/gs2_cons.png",
        "images/gs2_infil.png",
        "images/gs2_rdd.png",
        "images/gs2_rdd_fs.png",
        "images/gs2_pond.png",
        "images/gs2_heap.png",
    },
)
