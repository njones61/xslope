"""Page config: docs/verification/rs2.md (RS2 finite-element corpus).

Figure mode: ``panel``.  Every figure on this page is produced by one figure
builder with a fixed set of layouts — a 2x2 composite (inputs, mesh, max shear
strain, displacement vectors), a side-by-side LEM pair, or an inputs-only single
axes — so the caption makes a checkable claim about how many panels the reader
will see, and the classifier reads that claim back off the PNG.
"""
from ..config import PageConfig

CONFIG = PageConfig(
    name="rs2",

    # Cross-scope comparisons: the two numbers live in different sentences, or
    # in a table the sentence refers to.  Each entry spells both out and is
    # re-derived.  Eleven are deliberate vendor-against-vendor statements the
    # page makes in its own voice rather than XSLOPE-against-source.
    #   (printed, distinctive substring of the line, value-for, value-against)
    whitelist=[
        ('+2.1', '[RS2-63](#rs2-63)', '1.409', '1.38'),
        ('+4.8', 'the other φ = 0 foundation problem', '1.477', '1.41'),
        ('−1.7', 'All six cases land within', '1.091', '1.11'),
        ('−2.3', 'is the widest of the five', '1.026', '1.05'),
        ('+9.0', 'is not a search that stopped early', '0.169', '0.155'),
        ('+2.9', 'measurement of its zone here is', '1.07', '1.04'),
        ('−2.3', 'case 4 (−2.3%) is the widest', '1.026', '1.05'),
        ('−1.7', 'case 4 (−1.7%) is the widest', '1.091', '1.11'),
        ('+2.7', 'all six within ±2.7%', '0.339', '0.33'),
        ('−1.3', 'which XSLOPE sits', '1.323', '1.34'),
        ('+10.1', 'a facing-column discretization residual', '0.944', '1.05'),
        ('+2.1', 'The same uncapped machinery is within 2.1%', '1.381', '1.41'),
        ('+3.5', 'baseline is within 3.5% at every frame', '1.987', '2.06'),
        ('+7.7', '| Ordinary (OMS) |', '1.568', '1.456'),
        ('+7.3', '| Bishop simplified |', '1.601', '1.492'),
        ('+7.3', '| Spencer | 1.600', '1.600', '1.491'),
        ('+1.1', 'The two vendor numbers landing within', '1.9', '1.88'),
        ('+2.5', 'the three answers agree within', '1.25', '1.219'),
        ('+1.9', 'that rounding is worth', '0.997', '0.978'),
        ('+15.6', "above Slide2's own Spencer on an identical slope", '1.11', '0.960'),
        ('−4.7', 'at the tagged 12.4 ft mesh', '1.384', '1.452'),
        ('+3.0', 'one count in the last place', '0.34', '0.33'),
        ('+23.2', 'GEO FEM reads 23.2% above', '1.17', '0.95'),
        ('+11.3', 'above PLAXIS on #58 case 3', '0.59', '0.53'),
        ('+11.6', 'overshoots all', '7.836', '7.02'),
        ('+12.6', 'overshoots all', '7.836', '6.96'),
        ('−2.9', 'recovers 1.156', '1.156', '1.19'),
        ('+7.9', 'an 8 m cutoff takes', '1.219', '1.13'),
        ('+4.2', "the paper's 2.5 × 10⁻⁵", '2.605', '2.5'),
    ],

    # Absolute factor-of-safety differences whose partner value the page does
    # not print (a measurement stated on its own).
    abs_bounds=[
    ],

    # Percentages spelled in words.  Each is a QUALITATIVE characterisation, not
    # a comparison of two printed factors — adjudicated once and named here, so
    # that a new worded percentage fails until someone reads it.
    worded_ok=[
        "the rest move a percent or two",
        "percent or two a strength-reduction factor drifts under refinement",
        "agreeing to a quarter of a percent",
        "costs about a percent, not several",
    ],

    # Statements of a BOUND or a SHARE rather than a comparison of two printed
    # factors: no pair exists to re-derive, so each is named here explicitly.
    bounds=[
        ('+7', 'the two published factors differ by up to 7%'),
        ('+5', 'neither mesh is within 5%'),
        ('−2.0', '`RS2 SSRM 1.33'),
        ('+4.2', 'lighter by section area'),
        ('+1.8', "the polygon's area is within"),
        ('+3', 'within 3% of the vendor and/or reference figure'),
        ('+6', '| 🟡 | 3–6% |'),
        ('+6', 'more than 6% |'),
        ('+14', 'by 14% at 12.5 kPa'),
        ('+25', 'by 25% at 5 kPa'),
        ('+43', 'by 43% at 1 kPa'),
        ('+3', 'worth ~2–3% here'),
        ('+0.25', 'to within 0.25%'),
        ('+39.2', 'larger than the Mohr-Coulomb corridor'),
        ('+3.5', 'Every case lands within 3.5% of unity'),
        ('+1', 'lands four of the five deep stations within 1%'),
        ('+1.6', 'stages land within 1.6%'),
        ('+1', 'and four of the six within 1%'),
        ('+0.8', 'land within 0.8%'),
        ('+2.2', 'that still lands within 2.2%'),
        ('+1.2', "locked, within 1.2% of RS2's own SSR"),
        ('+2.3', "land within ±2.3% of RS2's M-C"),
        ('+1.7', "land within ±1.7% of RS2's M-C"),
        ('+2.1', 'lands within 2.1% of'),
        ('+2.3', "All five within ±2.3% of RS2's M-C"),
        ('+1.7', "All six within ±1.7% of RS2's M-C"),
        ('+2.7', 'all 17 cases land within ±2.7% of RS2'),
        ('+6', 'about 6% apart on cases 1 and 2'),
        ('+20', 'runs about 20% low'),
    ],

    locked_value_re=r"(?:SSRM|Spencer|Bishop|k꜀|FS)\s+\**(\d+\.\d{3})",

    figure_mode="panel",
    # Figures that are not model/mechanism plots: their captions describe a
    # diagnostic, so no panel form is claimed and none is required.
    caption_exempt={
        "images/rs2_67_fielddiff.png",
    },
)
