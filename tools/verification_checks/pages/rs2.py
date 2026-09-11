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
    # re-derived.  Several are deliberate vendor-against-vendor statements the
    # page makes in its own voice rather than XSLOPE-against-source.
    #   (printed, distinctive substring of the line, value-for, value-against)
    whitelist=[
        ('+0.8', '[RS2-63](#rs2-63)', '1.391', '1.38'),
        ('+5.5', 'the other φ = 0 foundation problem', '1.488', '1.41'),
        ('+9.0', 'is not a search that stopped early', '0.169', '0.155'),
        ('+2.7', 'the three locks within ±2.7%', '0.339', '0.33'),
        # RS2 against itself: the section's two vendor columns, read out of the
        # table the sentence sits under.
        ('+6.1', 'vendor answers are 6.1% apart', '1.05', '0.99'),
        ('+1.8', 'The same uncapped machinery is within 1.8%', '1.669', '1.64'),
        ('+3.2', 'baseline is within 3.2% at every frame', '1.713', '1.77'),
        ('+1.1', 'The two vendor numbers landing within', '1.9', '1.88'),
        ('+2.5', 'the three answers agree within', '1.25', '1.219'),
        ('+1.9', 'that rounding is worth', '0.997', '0.978'),
        ('+15.6', "above Slide2's own Spencer on an identical slope", '1.11', '0.960'),
        ('+3.0', 'one count in the last place', '0.34', '0.33'),
        ('+23.2', 'GEO FEM reads 23.2% above', '1.17', '0.95'),
        ('+11.3', 'above PLAXIS on #58 case 3', '0.59', '0.53'),
        ('+4.2', "the paper's 2.5 × 10⁻⁵", '2.605', '2.5'),
    ],

    # Absolute factor-of-safety differences whose partner value the page does
    # not print (a measurement stated on its own).
    # Published values quoted in prose beside the row's own lock: RS2's native
    # shallow rebuild of problem 43 is a different model's number.
    untagged_allow=[
        ('1.19', 'The 1.19 published for this problem'),
        # RS2-49 is reported without a lock, so its section prints no factor of
        # its own; these two are the published spread it is reported against.
        ('1.08', "against its own SSR of 1.08"),
        ('0.99', 'the referee at 0.99 between them'),
    ],

    abs_bounds=[
    ],

    # Percentages spelled in words.  Each is a QUALITATIVE characterisation, not
    # a comparison of two printed factors — adjudicated once and named here, so
    # that a new worded percentage fails until someone reads it.
    worded_ok=[
        "agreeing to a quarter of a percent",
        "costs about a percent, not several",
    ],

    # Statements of a BOUND or a SHARE rather than a comparison of two printed
    # factors: no pair exists to re-derive, so each is named here explicitly.
    bounds=[
        ('+4.2', 'lighter by section area'),
        ('+1.8', "the polygon's area is within"),
        ('+14', 'by 14% at 12.5 kPa'),
        ('+25', 'by 25% at 5 kPa'),
        ('+43', 'by 43% at 1 kPa'),
        ('+3', 'worth ~2–3% here'),
        ('+39.2', 'larger than the Mohr-Coulomb corridor'),
        ('+3.5', 'Every case lands within 3.5% of unity'),
        ('+1', 'stages land within 1% of'),
        ('+2.3', 'within 2.3% of RS2 at all three thicknesses'),
        ('+0.8', 'land within 0.8%'),
        ('+2.2', 'that still lands within 2.2%'),
        ('+2.1', 'lands within 2.1% of'),
        ('+20', 'runs about 20% low'),
        # RS2-18: a bound over both cases against the scored vendor column.
        ('+0.8', "lands within 0.8% of RS2's own model"),
    ],

    # Flac3D is one of the four programs the RS2-62 tables carry as published
    # columns; the base authority vocabulary does not name it.
    auth_hdr_extra=["Flac3D"],

    locked_value_re=r"(?:SSRM|Spencer|Bishop|k꜀|FS)\s+\**(\d+\.\d{3})",

    # Two section headings are named by summary rows carrying different dots.
    # In both cases the heading covers more than one row of the manual, so the
    # dot it opens with is the one the section's own locked comparisons set.
    heading_dot_multi=[
        # "RS2-39/41/43" builds problems 41 and 43, both 🟢; problem 39 (and its
        # Part IV twin, row 76) is deferred with no lock of its own, which is
        # what its ⊘ says.  The section's locked comparisons are the 🟢 pair.
        ('rs2-39', '🟢'),
        # "RS2-68" carries the three seismic cases together.  Parts I-III score
        # the problem 🔴 on case 3, and Part IV's rows split the same section:
        # case 1 (row 62) 🟢 and case 3 (row 63) 🔴.  The worst locked case sets
        # the dot, here and in the summary table.
        ('rs2-68', '🔴'),
        # "RS2 Part IV VP65 / VP66" carries the two upstream-pool dams of one
        # family in a single section, because what separates them is one
        # argument about how each is watered.  VP66 (row 66) is the section's
        # locked comparison, 🟢; VP65 (row 65) is reported against a
        # zone-constrained vendor factor and carries no lock, which is its ⊘.
        ('p4-vp65', '🟢'),
    ],

    figure_mode="panel",
    # Figures that are not model/mechanism plots: their captions describe a
    # diagnostic, so no panel form is claimed and none is required.
    caption_exempt={
        "images/rs2_67_fielddiff.png",
    },
)
