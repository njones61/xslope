"""Page config: docs/verification/rocscience.md (Rocscience Slide2 corpus).

Figure mode: ``panel``.  The corpus figures are a fixed two-panel composite —
model inputs on the left, the representative solution on the right — and every
caption says so ("inputs and ... solution", "inputs and the shallow mechanism"),
which is a layout claim the classifier reads back off the PNG.  The one
transient rapid-drawdown curve at VP102 is a single time-series axes whose
caption makes no layout claim, so it is named in ``caption_exempt``.
"""
from ..config import PageConfig

CONFIG = PageConfig(
    name="rocscience",

    bounds=[
        # the match-dot legend, not a comparison
        ('+3', 'within 3% of the vendor and/or reference figure'),
        ('+6', '| 🟡 | 3–6% |'),
        ('+6', 'more than 6% |'),
        # the conventions note's own illustration of the sign convention
        ('+2', 'Every printed difference is'),

        # ---- agreement bounds over a set of surfaces, cases or stations ----
        ('+0.3', 'Janbu/Spencer within 0.3% of Slide'),
        ('+0.3', 'Spencer within 0.3% on both'),
        ('−2.8', 'runs from −2.8% to +0.2%'),
        ('+0.2', 'runs from −2.8% to +0.2%'),
        ('−0.5', "benchmark diagnostic (−0.5% at Slide2's K"),
        ('−3.0', 'with the plain static-head `u=piezo`'),
        ('+0.3', 'Deterministic FS agrees within'),
        ('+0.7', 'and it lands within 0.7% of Slide on both circles'),
        ('+0.6', 'within 0.6% of Bishop'),
        ('+3', 'brackets the factor of safety within 1–3%'),
        ('+1', "Morgenstern–Price lands within 1% of W&H's 2.36"),
        ('+0.19', "matches SLOPE/W's own to 0.19% on all nine"),
        ('+0.5', "its Bishop matches Slide2's Table 35.2 column to 0.5%"),
        ('+372', 'which moves the sliding weight'),
        # the G/H transposition residuals: this page states the result of the
        # comparison made on the GeoStudio page, whose section prints the four
        # factors (6.059 / 6.293 and 11.561 / 10.576) and re-derives each there
        ('−42.7', "brings the residuals from"),
        ('+83.7', "brings the residuals from"),
        ('−3.7', "brings the residuals from"),
        ('+9.3', "brings the residuals from"),
        ('+0.6', "reproduces Slide's Bishop values within 0.6%"),
        ('+1.6', "Ng & Shi's own published Bishop values within 1.6%"),
        ('+2', 'halving or doubling a moves FS by'),
        ('+15', 'a lock could tolerate'),
        ('+19', 'FS with FE pore pressures is 14–19% lower'),
        ('+13', 'guts the artesian pressures and reads'),
        ('+40', 'u at the toe 40% above hydrostatic'),
        ('+65', 'and 65% at 5 ft depth'),
        ('+1', 'agrees with the Slide-figure trace within'),
        ('+0.5', 'Every method within 0.5%'),
        ('+23', 'and 23% high with no warning'),
        ('+0.6', 'The FE case lands within 0.6% of Slide'),
        ('+6', 'the critical circle is a shallow toe surface'),
        ('+2', 'on the exact vertices both methods land within'),
        ('+0.7', 'the whole family tracks Slide within 0.7%'),
        ('+1', 'and D&W within 1%'),
        ('+1', 'lands within 1% of all three Slide numbers'),
        ('+0.5', 'XSLOPE lands within 0.5% of Slide (1.527 vs 1.534)'),
        ('+3', 'worth about 3% of the factor of safety'),
        ('+3', 'every station sits within 3% of the Slide2 Spencer column'),
        ('+0.6', 'Every optimized value sits within 0.6% of Slide2'),
        ('+0.1', 'The two published columns agree to about 0.1%'),
        ('+0.8', 'Every other case agrees with Slide within 0.8%'),
        ('+2.4', 'and with the originating paper within 2.4%'),
        ('+0.5', "within 0.5% of Slide's limit-filtered search"),
        ('+0.7', "Slide's constrained block search lands within 0.7%"),

        # ---- probabilities of failure and reliability spreads ----
        ('+26.2', '26.2% against 25.3%'),
        ('+25.3', '26.2% against 25.3%'),
        ('+21.8', '21.8% against 21.0% and 21.4%'),
        ('+21.0', '21.8% against 21.0% and 21.4%'),
        ('+21.4', '21.8% against 21.0% and 21.4%'),
        ('+0.01', 'is at its resolution limit'),
        ('+0.04', 'is at its resolution limit'),
        ('+18', "reproduces Duncan's own 18% almost exactly"),
        ('+18', "sits below Duncan's 18%"),
        ('+18', '18% (Duncan 2000 TSPM)'),
        ('+33', '30–33% (D&W 2014'),
        ('+2.2', 'gives PF ≈ 2.2%'),
        ('+3.5', 'and PF ≈ 3.5%'),
        ('+0.15', 'above El-Ramly et al.'),
        # the two estimators run on one surface; neither estimate is printed
        ('+1', 'Taylor series and Monte Carlo on one surface agree to'),
        ('+3', 'an empirical probability of failure of about 2%'),
        ('+6.2', 'the three published estimates span'),
        ('+0.36', 'the three published estimates span'),
    ],

    abs_bounds=[
        # the crack cost is XSLOPE's own cracked-vs-uncracked difference; the
        # uncracked factor is not printed, only the published 1.77 / 1.69 pair
        # it is compared against
        ('−0.090', 'the crack costs'),
    ],

    whitelist=[
        # VP20: the theory value XSLOPE reproduces exactly is the table's 2.500,
        # and these two sentences measure the OTHER programs against it
        ('+1.4', 'reads ~1.4% high', '2.534', '2.500'),
        ('+3.6', 'sits ~3.6% below theory', '2.41', '2.500'),
        # Slide2's own Ordinary result against its own Spencer
        ('+44', 'sits 44%', '0.859', '0.596'),
        # theory bracket: Slide2 against the exact 1.0 the sentence states
        ('−5.9', 'against theory', '0.941', '1.0'),
        # Bishop/Spencer pair cells: the delta pair is printed as one cell, so
        # the FIRST member's two operands are named (the second reconciles from
        # the last number before it)
        ('−1.8', '| Piezometric line |', '1.070', '1.090'),
        ('−0.6', '| tangent 15 ft |', '1.389', '1.398'),
        ('−0.6', '| I: c<sub>u</sub> = 200 + 15·z |', '1.305', '1.313'),
        ('−0.5', '| II: c<sub>u</sub> = 300 |', '1.328', '1.335'),
        ('−0.4', '| II | 5 |', '0.905', '0.909'),
        ('−0.3', '| III | 10 |', '1.042', '1.045'),
        # L&H reinforcement rows: the vendor ran Ta = 10, so the delta is
        # against XSLOPE's own Ta = 10 value printed in the previous cell
        ('+0.9', '| 89 | L=4.2 m', '0.980', '0.971'),
        ('+0.2', '| 92 | Water hw=3 m', '1.039', '1.037'),
        ('+0.3', '| 93 | Surcharge q=20', '0.961', '0.958'),
        # MMO rows: the delta is against the PS+SA value, with the uni-modal
        # value printed between it and the delta
        ('+0.5', '| 1.4 | deep |', '1.221', '1.215'),
        # Newmark: the multi-modal row's own integration against Slide2's
        ('−0.5', 'and it agrees to −0.5%', '5.015', '5.042'),
        # Cai & Ugai D1/D = 3: the two published values against each other
        ('+4.4', 'Slide sits 4.4% above the paper', '1.43', '1.37'),
        # VP108: the free grid search against Slide's equivalent-cohesion Bishop
        ('+1.5', '1.5% below Slide', '1.761', '1.787'),
        # VP102 transient: the 300 h frame against Slide2's own published value,
        # with the committed field's factor printed between them
        ('−2.0', "against Slide2's 2.092", '2.051', '2.092'),
    ],

    worded_ok=[
        # a qualitative characterisation of the sensitivity endpoints, not a
        # comparison of two printed factors
        "within about a percent at every endpoint",
    ],

    # Percentages expressed as a share of a quantity, or as a parameter range
    not_a_comparison_extra=[
        r"of their means",
        r"of the pore force",
        r"of the fill strength",
        r"of the gabion fill",
        r"of the displacement",
        r"filled",
    ],
    # Quantities identified by what PRECEDES them: a coefficient of variation
    # is a value, not a difference.
    not_a_comparison_prefix=[
        r"COV(?: of)? ",
    ],
    # Reliability tables: the beta / probability-of-failure columns hold
    # probabilities, not comparisons
    share_hdr_extra=[r"\bPF\b", r"σ_F", r"RI_ln"],
    # Rows whose cells are percent CHANGES or probabilities rather than
    # comparisons; the sensitivity rows say so in their own Note cell
    share_rows=[
        '| ΔFS over the A range',
        '| ΔFS over the b range',
        '| β_ln → PF |',
    ],

    auth_hdr_extra=[r"Slide\b", r"Slide2", r"SLOPE/W", r"UTEXAS4", r"WINSTABL",
                    r"XSTABL", r"Cai & Ugai", r"Sheahan", r"El-Ramly",
                    r"Borges & Cardoso", r"C&X", r"W&H", r"Chowdhury",
                    r"Pockoski", r"theory"],

    tag_value_keys=['expected_fs*', 'fs_*', 'expected_beta', 'expected_kc',
                    'expected_flowrate*', 'expected_head*'],
    tag_round_dp=2,
    locked_value_re=None,

    figure_mode="panel",
    caption_rules=[
        ("inputs and", "two"),
        ("FS vs A and b", "two"),
    ],
    # the one time-series figure on the page: a factor-of-safety curve against
    # Slide2's own table, claiming no panel layout
    caption_exempt={
        "images/vp102t_curve.png",
    },
)
