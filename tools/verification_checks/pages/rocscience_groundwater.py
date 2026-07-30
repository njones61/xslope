"""Page config: docs/verification/rocscience_groundwater.md (RS2 groundwater).

Figure mode: ``panel``.  The corpus figures are a fixed two-panel composite —
mesh on the left, solved heads on the right — and the caption says exactly that,
so the caption makes a layout claim the classifier can read back off the PNG.
The transient/profile figures at the end of the page (isochrones, head along a
query line, an erfc profile) are plotted by different scripts at whatever panel
count the comparison needs and their captions make no layout claim, so each is
named in ``caption_exempt`` rather than left in a silent "other" bucket.
"""
from ..config import PageConfig

CONFIG = PageConfig(
    name="rocscience_groundwater",

    bounds=[
        # the match-dot legend, not a comparison
        ('+3', 'within 3% of the vendor and/or reference figure'),
        ('+6', '| 🟡 | 3–6% |'),
        ('+6', 'more than 6% |'),
        # gw4: the Q spread across the six-step target_size ladder; the
        # individual discharges are not printed, only the y₁ readings are
        ('+0.6', 'while Q moves'),
        # gw5: the same ladder's discharge spread (8.194 / 8.167 / 8.141)
        ('+0.6', 'a 0.6% spread, and'),
        # gw5: the render calibration, an agreement between two image scales
        ('+0.1', 'scales agree to'),
        # van Genuchten vs SEEP2D on the linear-front problems, as a set
        ('+0.15', 'agree on discharge to better than'),
    ],

    whitelist=[
        # gw5: "+1.8% on the closed form at the finest" — the finest-mesh
        # discharge 8.141 against the closed form 8.0, both printed in the
        # paragraph (and reachable from the summary row through its #gw5 link)
        ('+1.8', 'on the closed form at the finest', '8.141', '8.0'),
        ('+1.8', 'on the finest mesh', '8.141', '8.0'),
        # the 3.5–4.7% van Genuchten band, quoted ahead of the two case pairs
        # it summarises (gw009a 2.299 against SEEP2D's 2.412)
        ('+4.7', 'the total discharge reads 3.5–4.7% below SEEP2D', '2.299', '2.412'),
    ],

    # Percentages expressed as a share of a reference quantity
    not_a_comparison_extra=[
        r"of \$u_0\$",
        r"of u₀",
        r"of the driving head",
        r"of the profile's head range",
        r"of the 0\.25 m divide mound",
        r"of the 0\.56 m head range",
        r"of the 45 m dam",
    ],

    auth_hdr_extra=[r"Slide\b", r"SEEP2D", r"Eqs? 4\.", r"Vedernikov", r"Clement",
                    r"Kozeny", r"Terzaghi", r"Pyrah", r"Ferris"],

    tag_value_keys=['expected_fs*', 'fs_*', 'expected_flowrate*',
                    'expected_head*', 'expected_beta', 'expected_kc'],
    # Regression locks on XSLOPE's own discharge, which the page deliberately
    # does not publish as a comparison: GW2 says "flowrate is locked as a
    # regression value" (the manual publishes head points, not Q), and GW11
    # says the tag "locks XSLOPE's own flowrate rather than a published one"
    # because the manual does not print k_s for that case.
    tag_exempt=[
        ('4.534e-06', 'benchmark=GW2-q'),
        ('7.814e-07', 'benchmark=GW11-q'),
    ],
    # the summary rows quote discharges, heads and offsets in physical units
    # rather than a fixed-precision factor of safety, so there is no
    # locked-value pattern to sweep against the tags
    locked_value_re=None,

    figure_mode="panel",
    caption_rules=[
        ("mesh and solved heads", "two"),
    ],
    # Series and profile figures: the caption names what is plotted against
    # what, and claims no panel layout.
    caption_exempt={
        "images/gw015.png",
        "images/gw016.png",
        "images/gw017.png",
        "images/gw018.png",
        "images/gw019.png",
        "images/gw020.png",
        "images/gw021.png",
    },
)
