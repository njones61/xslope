# Verification-page checks

Standing checks on the six verification pages under `docs/verification`
(`rs2`, `rocscience`, `geostudio`, `rocscience_groundwater`, `ssrm`, `seep`).
They run as one test in the default `run_tests.py` set, and they are cheap:
a page that has not changed since it was last certified costs one file read.

## What each check does

**Deltas** (`deltas.py`). Every percentage the page prints — `(+1.2%)`,
`−7%`, `**+2.1%**`, `<0.01%` — and every absolute factor-of-safety difference
("worth +0.013", "costs 0.18") is re-derived from two numbers the page itself
prints in the same sentence or the same table row:

    delta = (XSLOPE − source) / source × 100,  half-up, one decimal

Pairing is sign-aware and nearest-first, so a sign flip or a wrong last digit
cannot be certified by some other pair further down the paragraph. The check
also enforces the page convention that a signed delta is printed to one
decimal unless it is explicitly hedged ("to within −0.27%"), and reports
percentages spelled in words, which no arithmetic can reach.

Percentages that are not comparisons — a share of the mesh, a probability of
failure, a degree of consolidation — are excluded by named phrases, never by a
bare word, so a real comparison cannot hide behind an incidental noun.

**Tags** (`tags.py`). Three passes. Forward: every value a `<!-- test: -->`
tag locks must be printed in the section that carries the tag. Restatements:
every number a section attributes to XSLOPE must agree with a lock that
attribution reaches. Reverse: every value the page presents as locked — a factor
of safety in the Results column of a summary row whose Match cell carries a
colour dot — must have a tag behind it.
The tag is truth: when a tag and the text disagree, the text is what changes.

The forward pass is existential — it asks whether the lock is printed, not
whether every printing of it is right — so a section that states its locked
factor of safety twice could drift in one place with the other left standing.
The restatement pass is universal over the numbers the section claims as its
own results, which is a narrower set than "every number shaped like a factor of
safety":

* a cell in a results table's XSLOPE column, where the row's label names the
  locked method — the method a tag names in its key (`fs_bishop`), in
  `method=`, or in a `type=` that says SSRM, narrowed further by the element
  type it meshes on and by a tangency constraint on the search;
* a number a sentence introduces with "XSLOPE's" or "this model's", up to the
  next sentence boundary or the next mention of a source.

A number the page attributes to a source, and a number it attributes to nobody,
are left where they already are — to the delta check and to the untagged sweep.
A table with several XSLOPE columns (`XSLOPE composite` beside `XSLOPE
circles-only`) publishes a row's method under variants the row label does not
name, and a number a cell puts in parentheses (`1.415 *(search 1.411)*`) is a
second quantity stated there; neither is read. The locks in scope are the
section's own, a page-level bank tag whose input file the section links, any
page's tag on a file the section links, and the tags of a section this one
cross-references by anchor — the link being the page's own statement that the
two present the same problem. A restatement passes when it is the lock rounded
to a precision the page allows (`tag_round_dp`) or inside the tag's own
tolerance.

A tag value that is a semicolon **list** locks one value per element, and every
element is checked the same way — `expected=1.686;0.941;…` (a factor of safety
per march step) and `points=20:5:7.166;…` (a solved head per station, the
coordinates being inputs like `time=`). So a wrong digit in one cell of an
eleven-row table fails. Two restatements count as printing the lock, because
each is the same number read the way the comparison is made: the value rounded
to the page's own precision (`tag_round_dp`), and, for a `points` element, the
pressure head ψ = h − y at that station.

Most locks are factors of safety, and two are not. A `slip_depth` tag locks
`expected_depth`, the depth of the critical slip surface below a stated point,
in length units; a `support_force` tag locks `expected_force`, the horizontal
support force per meter that holds a slope at a target factor of safety. Both
keys are read exactly like `expected_fs`: the value must be printed in the
section carrying the tag, at a precision `tag_round_dp` allows, and the untagged
sweep treats a number printed under one of them as guarded. Neither is read by
the restatement pass, which is about factors of safety; a depth or a force
printed in an XSLOPE column belongs to the forward pass and to the delta check.
The delta check reads no tag keys at all, so the two need nothing from it — a
percentage beside a depth is re-derived from the two depths the row prints.

A page may publish only part of a probe set — a table of the stations that carry
the argument, or none of a set that guards a field's shape while the section
publishes something else. `tag_list_published` declares how many elements the
page prints, per lock and with the reason; the count is exact, so a mistyped
printed value fails just as an undeclared one does, and a declaration that
matches no tag is reported dead.

**Voice** (`voice.py`). These pages are documentation for a stranger who wants
to know whether XSLOPE reproduces a published problem and what the caveats are.
They are not addressed to the maintainer and they are not a record of how the
work went. Prose written mid-investigation reads differently — "the obvious
suspect was...", "their factors are withdrawn", "four candidate causes have been
measured" — and a reader who was not there gets nothing from it. The check
carries a curated list of phrasings that are campaign voice in essentially any
context, in four groups: first person, project process, investigation narrative,
and time measured from the project rather than from the problem. Each hit names
the page, the line, the phrase and the line's text.

It is deliberately narrow. It flags phrasings, not everything that could be
written better, and a phrase that has an ordinary descriptive reading is left
out rather than exempted case by case — "the firm base now sits at depth D" is
not campaign voice, so a bare "now sits" is not banned. Prose only: fenced
blocks, HTML comments (which is what a test tag is), inline code spans, link
targets and bare URLs are removed before matching, so a material named `us` or a
citation URL cannot fire. Where a page genuinely needs a banned phrase,
`voice_allow` names it — `(phrase, distinctive substring of the line)`, both
required, and an allowance that never fires is reported dead.

**Capability negations** (`capabilities.py`). A row that cannot be built has to
say why, and the cheapest thing to write is that XSLOPE has no input for it.
Written from memory that sentence is wrong about as often as it is right: the
`lloads` sheet has carried a concentrated line load (x, y, P, Angle) since v14
and `fem.py` applies it, and two public texts still said the loads sheet could
not carry one. A reader takes such a sentence as the capability inventory it
looks like, and a blocked row built on it never gets revisited.

So every sentence that negates a capability — *does not carry*, *cannot model*,
*has no counterpart*, *not supported*, *not implemented*, *has no input* — must
name a capability on `ABSENT`, the one list of things XSLOPE genuinely does not
have. Each entry is `(capability, page, marker, evidence)`, and the evidence
field records where the capability was looked for and not found, so the next
person to doubt an entry re-runs the grep instead of the reasoning. A sentence
that reaches no entry fails: either the absence is real and its entry is
missing, or the page is wrong about XSLOPE and the prose is what changes. An
entry no sentence cites is reported dead, like every other exemption list here.

A negation fires only when it is about XSLOPE. The sentence must name XSLOPE or
one of its artifacts — a sheet of the input template, the loader, the importer,
an engine — except for the phrasings that are about XSLOPE by construction (a
row's *not supported* verdict, "XSLOPE has no ...", "has no input"). A negation
whose subject is the source ("a quantity the example does not carry", "storage
the vendor's model does not have") and one whose object is the verification work
("cannot carry the comparison", "no published counterpart", "cannot be locked")
are both left alone, and so is the status legend that merely lists the verdict
words. Vendors run engines with XSLOPE's names, so "RS2's own SSRM does not
model it" is read as the vendor's.

**Untagged numbers** (`untagged.py`). A section prints three kinds of number
that carry an argument: the factor of safety a tag locks, the value the source
published, and the comparison between them. Anything else shaped like a factor
of safety — a mesh-sweep row, a depth-cutoff row, a with/without variant, a
reading taken off a field at the critical strength — is a companion measurement
that nothing regenerates and nothing defends when a lock moves. The check reads
each section and reports every factor-of-safety-shaped number that is neither
within the tolerance of a tag the section carries (its own, a page-level bank
tag naming a file the section links, or another page's tag on that same file)
nor printed in a column whose header names the source. Inputs, dimensions,
percentages, figure and section numbers, code, math and link targets are taken
out of the running first.

It **reports and returns zero** while `untagged.ENFORCING` is False. Each flag
is a sentence someone has to read — a companion measurement to trim, a number
worth a tag of its own, or a quantity that only looks like a factor of safety
and belongs in `untagged_allow`. Enforcing it before that reading is done would
push the pages toward blanket allowances, which is the opposite of the point.

**Figures** (`figures.py`). Two modes, chosen per page.

* `panel` — the panel layout is read directly off the PNG (an ink-profile test
  for inter-panel gutters) and compared against the layout the caption declares.
  Used where the caption makes a layout claim: rs2's four-panel composites,
  and the "inputs and ... solution" / "mesh and solved heads" two-panel
  composites on `rocscience`, `geostudio` and `rocscience_groundwater`.
  A caption the classifier cannot read is a failure, not a pass, unless the
  figure is named in the page's `caption_exempt`.
* `structural` — used on `ssrm` and `seep`, whose figures are single-axes plots
  captioned with their file name: the caption makes no layout claim, so panel
  classification would be testing a claim the page does not make. Every
  referenced image must exist and carry non-empty alt text.

**Heading dots** (`dots.py`).  The match dot a summary table gives a problem —
🟢 🟡 🔴 🟣 or the ⊘ that means no data, scored as
`docs/verification/index.md#how-the-match-dots-are-scored` describes — also opens
that problem's section heading, so a reader skimming the page sees the status
without going back to the table:

    ### 🟢 RS2-1: Simple slope stability assessment {#rs2-1}

The table stays the single source: the check reads each row's anchor and dot and
requires the heading carrying that anchor to start with it, written in a heading
as the bare ⊘ rather than the table's styled span.  Only headings with an
explicit `{#anchor}` that a row names are checked; a heading no table names is
left alone.  A row linking another page names a section that page owns, and a row
whose anchor is an inline `<a id=...>` has no heading of its own to open with a
dot — the first is skipped, the second reported as a note.  A row whose anchor
the page defines nowhere is a broken link and fails.  `--fix` rewrites the
headings, touching only the leading dot.

Several rows may name one heading — a section covering three problems of the
manual, a catalog row piggybacking on the section another row built.  Where they
agree the heading takes their dot; where they disagree `heading_dot_multi` names
which of those rows speaks for the section.

**Section order** (`order.py`).  The sections a summary table links appear down
the page in the order the table lists them, so following the table is reading the
page rather than jumping about.  The table decides the order, as it decides the
dots — the check only says the page agrees with it — and the order drifts the
moment a row is written last and appended at the end, which is how it happened.

Only anchored sections a row names are read, so a methodology section or a shared
discussion may sit anywhere; a row pointing at another page names a section that
page owns; a section several rows name is ordered by the first of them; and a row
whose Notes carry the status term *covered* is a cross-reference rather than the
owner, so it orders nothing.  The report names the pair that is out of order and
both line numbers.  There is no `--fix`: moving a section means moving its prose,
its test tag and its figure together, and a mechanical reshuffle is how a caption
ends up under the wrong image.

## What a strength-reduction tag says

A `type=fem_ssrm` tag is a locked factor of safety plus everything needed to
reproduce it: the model (`file`), the mesh (`element_type`, `target_size`, and
the `refine_*` keys), the search (`f_min`, `f_max`, `tolerance`, `max_iter`),
and the solver options the vendor model implies (`k0`, `tension_srf`,
`ssr_exclude`, `ssr_zone`, `elastic_materials`, `min_slip_depth`,
`suction_phi_b`, `suction_cap`). `benchmark` names the row. Everything in the
list reaches `solve_ssrm` through one function — `run_tests.build_fem_ssrm_case`
— which is also what the figure producers call, so a figure and its lock cannot
solve different problems.

Four keys say how the row is RUN rather than what it solves:

* **`check=edges`** checks the lock with the two trials its bisection closed on
  instead of re-running the bisection. It requires both of the next two.
* **`f_stand`** is the highest trial factor at which the model stood when the
  lock was cut, **`f_fail`** the lowest at which it failed. The locked value is
  that bracket's midpoint, so the pair is the lock; a run in which the model
  still stands at the one and still fails at the other has reproduced it.
  Written by `tools/lock_edges.py` from the trial record the figure producers
  persist (`*_fem_meta.json`), never by hand — a pair invented from the lock and
  the tolerance is a guess at where the bisection closed, and a guess one bracket
  step out passes while the lock moves underneath it.
* **`tier=gate`** holds the row out of the standard run: it is checked at the
  release gate (`run_tests.py --gate`) or when it is named with `--benchmark`.
  It is for rows that cost hours — the RS2 joint corpus — and it does not hold
  back a row that also carries an edge pair, since two trials is cheap enough to
  run routinely and holding it back would leave the lock unchecked between
  releases.

The `lock_edges` suite row re-checks every `check=edges` tag: both factors
present, `f_stand < expected_fs <= f_fail`, and the pair no wider than twice the
tag's `tolerance` — the bisection stops inside the tolerance, so a wider pair is
not a final bracket at all. `tools/ssrm_trial_audit.py` reads the same records
for the other question they answer: which locks were cut at the iteration
ceiling rather than at a mechanism.

## Tutorial restatements

`tutorials.py` sweeps the 28 tutorial pages under `docs/tutorials`, which carry
no configs and no certification manifest. A tutorial walks the reader through
runs it prints the answers to, and those answers drift the moment a solver round
moves them: the sample table and the verification section are re-measured
because a tag guards them, and the tutorial keeps printing what it printed the
day it was written.

The locks in scope for a page are its own tags plus every tag anywhere under
`docs/` on a model file the page links, so LEM-3 inherits the seven method locks
its workbook carries on `docs/lem/samples.md` and LEM-9 inherits the Rocscience
lock on the vendor model it borrows. Every factor-of-safety-shaped number the
page attributes to a run of its own then lands in one of three buckets:

- **guarded** — it restates a lock in scope, verbatim or correctly rounded to
  fewer places. Agreement is by printed form, not by the tag's tolerance: a page
  printing 1.313 where the lock reads 1.314 is restating it wrongly however
  small the gap, and that digit is the failure the sweep exists for.
- **disagreeing** — the number's own column header, row label or sentence names
  a method (or the strength-reduction run), and none of the locks it names is
  what it prints.
- **unguarded** — nothing in the docs locks the number. A sweep row, a variant
  run, a reading off a solved field: reproducible only by hand.

Two rules keep the identity honest where a tutorial's tables are read. A column
header binds its cells only where the table has ONE body row, because several
rows are variants of each other and a lock belongs to one of them; and a table
whose columns after the first are headed by numbers is a sweep, whose row label
names the method for the whole curve rather than for any column. Inside a code
fence — a verbatim console transcript — only the numbers the log itself labels
`FS`, `FOS` or `SRF` are read, and lines reporting an intermediate search
iteration are skipped.

The sweep REPORTS. `run_tests.py --tutorials` prints the per-page tally and
passes; it fails only if the checker raises. A tutorial legitimately re-runs a
sample under settings the sample's own tag does not use, so every finding is a
sentence someone reads before it is a defect.

## Running them

```bash
python tools/lock_edges.py                                   # what can be checked on its edges
python tools/lock_edges.py --missing                         # and what cannot, with reasons
python tools/lock_edges.py --write                           # write the pairs into the tags
python tools/ssrm_trial_audit.py --all                       # which locks were cut at the ceiling
python -m tools.verification_checks.certify                  # all six pages
python -m tools.verification_checks.certify rs2 seep         # named pages
python -m tools.verification_checks.certify --force          # ignore the manifest
python -m tools.verification_checks.mutations                # the mutation suite
python -m tools.verification_checks.tutorials                # the 28 tutorials
python -m tools.verification_checks.tutorials docs/tutorials/lem03_layered_slope.md
```

Each check can also be run on its own page for a detailed report:

```bash
python -m tools.verification_checks.deltas docs/verification/rs2.md
python -m tools.verification_checks.tags docs/verification/rs2.md
python -m tools.verification_checks.figures docs/verification/rs2.md
python -m tools.verification_checks.voice docs/verification/rs2.md
python -m tools.verification_checks.untagged docs/verification/rs2.md
python -m tools.verification_checks.dots docs/verification/rs2.md
python -m tools.verification_checks.dots --fix rs2      # rewrite the headings
```

## The recertify workflow

`certified.json` records the SHA-256 of each page's content as it stood when
the checks last passed **and a developer signed off**. A page whose hash still
matches is reported "unchanged, certified" and nothing else runs; a page whose
hash differs is re-checked in full, and even if every check passes it stays a
failure until the manifest is updated — certifying a page is a deliberate act,
because the developer is the one who read the flags while editing.

So: edit a page, fix whatever the checks raise, then run
`python -m tools.verification_checks.certify --recertify <page>` and commit the
updated `certified.json` in the same commit as the page change.

## Adding an exemption honestly

Each page's config lives in `pages/<page>.py`. An exemption is a claim about
the page, not a way to silence the checker, so:

* **Name both operands.** A `whitelist` entry is
  `(printed, distinctive substring, value-for, value-against)` and the checker
  still does the arithmetic — a wrong whitelist entry fails like anything else.
  It also re-checks that both operands are still printed where the claim can
  reach them, so editing one of them cannot silently orphan the comparison.
* **Say why the pairing is legitimate** in a comment above the entry: the two
  numbers sit in different sentences, or the delta is measured against a value
  printed in a table the sentence links to, or the page is deliberately
  comparing one vendor against another in its own voice.
* **Use `bounds` only where no pair exists.** "All six stages land within 1.6%"
  is a claim about a whole set; no single pair can check it, so it is
  adjudicated once, by hand, and named. The same goes for `abs_bounds` (an
  absolute FS difference whose partner value the page does not print) and
  `worded_ok` (a percentage spelled in words).
* **`not_a_comparison_extra` / `not_a_comparison_prefix` / `share_hdr_extra` /
  `share_rows` name quantities that are values, not differences** — a share of
  the domain, a probability of failure, a row of percent changes. Each must be
  a phrase or a named row, never a bare word.
* **`voice_allow` names a line where a banned phrase is legitimate** — usage
  prose addressed to the reader, a quoted source, a name that collides with a
  banned word. Say which of the four groups the phrase is in and why this line
  is not that. Adding one because the lint complained, without reading the
  sentence, is the failure the dead-allowance report exists to prevent.
* **`untagged_allow` names a number that only looks like a factor of safety** —
  a strength ratio, a stability number, a published quantity quoted in prose
  rather than tabulated. `(the number as printed, distinctive substring of its
  line)`. It is not a place to park a companion measurement: those are trimmed,
  or given a tag of their own.
* **`heading_dot_multi` names the dot a heading carries where the summary rows
  naming it disagree** — `(anchor, dot)`, and the dot must be one those rows
  actually give the anchor.  Say which of them the section's own locked
  comparisons come from.
* **`ABSENT` in `capabilities.py` is the one list of capabilities XSLOPE does
  not have**, and it is not a page config: an absence is a fact about the
  software, so the same capability carries one entry per page that cites it and
  the entries sit together. Every entry is written AFTER the grep, and the grep
  goes in the `evidence` field — which module was read, which template sheet was
  looked at, what was found instead. "I do not think we have it" is what this
  check exists to stop.
* **`tag_exempt` names a coverage lock the page deliberately does not print** —
  a tag that exercises a code path rather than backing a published number. The
  page normally says so in prose; quote that reason in the comment.
* **`tag_list_published` names how much of a list lock the page prints** — the
  count, not a licence to skip the lock. Say which elements are published and
  why the rest are not (a regression probe set, a station the table does not
  tabulate, a quantity the section publishes in another form). Raising a count
  because the checker complained, without knowing which value moved, is exactly
  the failure the exact count exists to prevent.

An exemption that never fires is reported as a **dead exemption** and fails the
check. That is deliberate: it is what stops the lists silently accumulating
entries for text that no longer exists.

## Mutation suite

`mutations.py` plants one defect at a time — a wrong last digit, a flipped
sign, an operand moved out from under a certified claim, a caption that no
longer matches its figure, a section heading whose dot disagrees with its summary
row, a tagged value dropped from its section, one element
of an eleven-value row corrupted on the page or in the tag, a slip-surface
depth or a support force moved on the page or in the tag, a planted dead
exemption, a tutorial number moved off the lock it restates, a capability
negation no ABSENT entry covers — and requires the
checks to catch every one. It also plants edits
that must **not** be flagged (a value reprinted at a different, correct
precision), because a check that fails on those would push the pages toward
printing tag values verbatim instead of at the precision each comparison is read
at. Run it after any change to the check logic; a gate that certifies a wrong
number is worse than no gate.

The tutorial block (`T1`–`T4`) moves one number the sweep reads as guarded and
names the bucket it must fall into — `disagreeing` where the page ties it to a
method, `unguarded` where nothing does — and requires it to leave the guarded
bucket either way. Its two controls run the other way: `N-T1` restates a lock
one place coarser and must add no finding, and `N-T2` moves a flagged number onto
a lock in scope and must remove exactly one.
