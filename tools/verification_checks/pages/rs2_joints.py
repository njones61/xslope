"""Page config: docs/verification/rs2_joints.md (RS2 joint-analysis corpus).

Figure mode: ``panel``.  The page's figures come from the same builder the RS2
corpus page uses (``make_rs2_joint_figures`` delegates to ``make_rs2_figures``),
so every caption makes the same checkable claim about how many panels the reader
will see.
"""
from ..config import PageConfig

CONFIG = PageConfig(
    name="rs2_joints",

    whitelist=[
    ],

    # Problem 16 is scored in a tilt angle, and the numbers its section prints
    # are seismic COEFFICIENTS and one unit-weight scale factor -- inputs to the
    # sweep, shaped like factors of safety and not one.
    untagged_allow=[
        ('0.1798', 'The stack stands at k'),
        ('0.1814', 'The stack stands at k'),
        ('0.25', 'fails at every one above'),
        ('1.016', 'scaled by 1.016 and then by 0.5'),
        # the plowing methodology bullet: Alejano's Eq. (7) and the rigid-block
        # ceiling on problem 11 are referee values, not XSLOPE results
        ('1.76', 'on problem 11 it gives 1.76'),
        ('1.21', 'against a ceiling of 1.21'),
        # problem 2: the 0.76 is what Alejano & Alonso PRINT for their own
        # Goodman & Bray recursion; the page recomputes it as 0.7734.
        ('0.76', 'against the 0.76 Alejano'),
        ('0.76', 'their Goodman & Bray 0.76,'),
        # RJ-20 stands at 2.512 and its search does not close above it: the
        # value is the standing edge in the committed run record (rj020_fem_meta),
        # shown as a lower bound by ruling, never tagged.
        ('2.512', 'unconfirmed'),
        ('2.512', 'The slope stands at 2.512'),
        ('2.531', 'the trial at 2.531'),
        # RJ-20's block density: blocks per square meter of the vendor's network
        # (525 blocks over 4,400 m2), not a factor of safety.
        ('0.119', 'blocks per square meter'),
        # RJ-3 on the chapter's printed list (rj003_chapter, RJ-3c) stands at 1.027
        # and its search does not close above it (the trial at 1.047 is
        # undecided; 1.125 fails): trials in the committed run record
        # (rj003_chapter_fem_meta), shown as a lower bound by ruling, never tagged.
        ('1.027', 'unconfirmed*, at least 1.027'),
        ('1.125', 'and fails at 1.125'),
        ('1.047', 'the trial at 1.047'),
        # The joint-cohesion-removed readings of problems 3 and 5 (problem 7's
        # equals its tagged chapter-list value): scratch runs
        # of the vendor's files with every joint's c set to zero, kept in the
        # private reports (r44_data/zero_cohesion), shown by ruling, never tagged.
        ('1.115', "everything else as the vendor's file: 1.115"),
        ('1.799', "everything else as the vendor's file: 1.799"),
    ],

    abs_bounds=[
    ],

    worded_ok=[
    ],

    bounds=[
        # the vendor's token toe force on problem 1, as a share of the block's
        # own weight: one quantity against one, not two printed factors
    ],

    # Goodman & Bray, Alejano, Lorig & Varona and UDEC are the referees this
    # corpus scores against; the base authority vocabulary does not name them.
    auth_hdr_extra=["UDEC", "Goodman", "Alejano", "Lorig", "Barla", "Hammah"],

    locked_value_re=r"(?:SSRM|FS)\s+\**(\d+\.\d{3})",

    figure_mode="panel",

    # A jointed model's composite names its own panels: the displacement-vector
    # panel is a scaled deformed section, and the strain panel carries the joint
    # slip overlay. The base classifier reads "shear strain" + "displacement
    # vector" + "mesh" as the four-panel form and does not know either phrase, so
    # the two the jointed figures use are declared here.
    caption_rules=[
        ("joint slip at the critical srf", "four"),
        ("deformed section", "four"),
    ],
)
