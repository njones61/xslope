"""Per-page configuration for the verification-page checkers.

A page config carries everything that is a property of ONE verification page:
its exemption lists (each naming both operands of a comparison the arithmetic
cannot reach), the extra vocabulary its tables use, and the figure conventions
its images follow.  The check modules themselves are page-agnostic.

An exemption is a claim about the page, not a way to silence the checker: every
entry names the two values being compared, and the whitelist entries are still
re-derived arithmetically, so a wrong exemption fails like anything else.  An
exemption that never fires is itself reported as a failure (a "dead exemption"),
which is what keeps the lists honest as pages are edited.
"""
from dataclasses import dataclass, field
from typing import List, Tuple


@dataclass
class PageConfig:
    """Everything the three checkers need to know about one page."""

    #: page identity — the file name under docs/verification without .md
    name: str

    # ---------------------------------------------------------------- deltas
    #: Cross-scope percentage comparisons: the two numbers live in different
    #: sentences, or in a table the sentence refers to.  Each entry is
    #: ``(printed, distinctive substring of the line, value-for, value-against)``
    #: and is re-derived; both operands must still be printed somewhere the
    #: claim can reach, or the entry is reported as orphaned.
    whitelist: List[Tuple[str, str, str, str]] = field(default_factory=list)

    #: Statements of a BOUND or a SHARE rather than a comparison of two printed
    #: factors ("all six land within 1.6%"): no single pair exists to re-derive,
    #: so each is named by ``(printed, distinctive substring)``.
    bounds: List[Tuple[str, str]] = field(default_factory=list)

    #: Absolute factor-of-safety differences whose partner value the page does
    #: not print ("the cutoff is worth −0.025").  ``(printed, substring)``.
    abs_bounds: List[Tuple[str, str]] = field(default_factory=list)

    #: Percentages spelled in words.  No arithmetic can be attached to them, so
    #: each is adjudicated once and named here by a distinctive substring; a new
    #: worded percentage fails until someone reads it.
    worded_ok: List[str] = field(default_factory=list)

    #: Extra regex alternatives appended to the "this percentage is not a
    #: comparison" pattern (domain shares, mesh fractions, ...).  Each must be a
    #: phrase, never a bare word.
    not_a_comparison_extra: List[str] = field(default_factory=list)

    #: Extra regex alternatives, anchored to the text IMMEDIATELY BEFORE a
    #: percentage, that mark it as not a comparison — a probability of failure,
    #: a degree of consolidation.  Some quantities are qualified by what precedes
    #: them rather than by what follows.  Each must be a phrase, never a bare
    #: word, and each is matched with `$` against the 40 characters before the
    #: token.
    not_a_comparison_prefix: List[str] = field(default_factory=list)

    #: Table rows whose percentages are values rather than comparisons — a
    #: row of percent CHANGES in a sensitivity sweep, a beta-to-probability
    #: row.  Each entry is a distinctive substring of the row (normally its
    #: label cell), so the exemption names the row it excuses.
    share_rows: List[str] = field(default_factory=list)

    #: Extra regex alternatives appended to the share-column header pattern: a
    #: table column whose header declares it holds shares (a degree of
    #: consolidation, a percentage of a total), not comparisons.
    share_hdr_extra: List[str] = field(default_factory=list)

    #: Extra regex alternatives appended to the authority-column header pattern
    #: (vendor and reference names this page's tables use as the source column).
    auth_hdr_extra: List[str] = field(default_factory=list)

    # ------------------------------------------------------------------ tags
    #: Tag keys whose values must appear verbatim in the section carrying the
    #: tag.  Prefixes end with '*'.
    tag_value_keys: List[str] = field(
        default_factory=lambda: ['expected_fs*', 'fs_*', 'expected_beta',
                                 'expected_kc', 'expected_flowrate*',
                                 'expected_head*'])

    #: Decimal places a section may restate a tag value to.  ``None`` demands
    #: the tag value verbatim (the rs2 convention: the page prints the lock).
    #: Set to 2 on a page whose prose quotes locks at the source figure's own
    #: precision, so ``expected_fs=1.822`` is satisfied by a printed ``1.82``.
    #: Only rounding DOWN in precision is accepted, and only to this many
    #: places, so a genuinely different value still fails.
    tag_round_dp: int = None

    #: Tag values the page deliberately does not print — coverage locks that
    #: exercise a code path rather than back a published number.  Each entry is
    #: ``(value, distinctive substring of the tag line)``; an entry that never
    #: fires is reported, like a dead delta exemption.
    tag_exempt: List[Tuple[str, str]] = field(default_factory=list)

    #: Regex naming the values a summary row presents as locked, for the
    #: reverse sweep (every locked value must have a tag behind it).  ``None``
    #: disables the reverse sweep on pages that carry no dotted summary table.
    locked_value_re: str = None

    #: Substring identifying the summary block the reverse sweep reads.
    summary_marker: str = 'corpus-summary'

    # --------------------------------------------------------------- figures
    #: 'panel'      — classify each PNG's panel layout and compare it against
    #:                the panel form the caption declares.
    #: 'structural' — the page's figures have no single meaningful panel form,
    #:                so check that every referenced image exists and that no
    #:                caption is empty.
    figure_mode: str = 'panel'

    #: Images whose captions describe a diagnostic and so claim no panel form.
    caption_exempt: set = field(default_factory=set)

    #: Caption phrases that force a panel form, checked before the generic
    #: classifier.  ``[(substring, form)]`` with form in {'one','two','four'}.
    caption_rules: List[Tuple[str, str]] = field(default_factory=list)
