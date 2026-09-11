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

    untagged_allow=[
    ],

    abs_bounds=[
    ],

    worded_ok=[
    ],

    bounds=[
        # the vendor's token toe force on problem 1, as a share of the block's
        # own weight: one quantity against one, not two printed factors
        ('+0.05', 'which is 0.05% of it'),
    ],

    # Goodman & Bray, Alejano, Lorig & Varona and UDEC are the referees this
    # corpus scores against; the base authority vocabulary does not name them.
    auth_hdr_extra=["UDEC", "Goodman", "Alejano", "Lorig", "Barla", "Hammah"],

    locked_value_re=r"(?:SSRM|FS)\s+\**(\d+\.\d{3})",

    figure_mode="panel",
)
