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

**Tags** (`tags.py`). Both directions. Forward: every value a `<!-- test: -->`
tag locks must be printed in the section that carries the tag. Reverse: every
value the page presents as locked — a factor of safety in the Results column of
a summary row whose Match cell carries a colour dot — must have a tag behind it.
The tag is truth: when a tag and the text disagree, the text is what changes.

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

## Running them

```bash
python -m tools.verification_checks.certify                  # all six pages
python -m tools.verification_checks.certify rs2 seep         # named pages
python -m tools.verification_checks.certify --force          # ignore the manifest
python -m tools.verification_checks.mutations                # the mutation suite
```

Each check can also be run on its own page for a detailed report:

```bash
python -m tools.verification_checks.deltas docs/verification/rs2.md
python -m tools.verification_checks.tags docs/verification/rs2.md
python -m tools.verification_checks.figures docs/verification/rs2.md
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
* **`tag_exempt` names a coverage lock the page deliberately does not print** —
  a tag that exercises a code path rather than backing a published number. The
  page normally says so in prose; quote that reason in the comment.

An exemption that never fires is reported as a **dead exemption** and fails the
check. That is deliberate: it is what stops the lists silently accumulating
entries for text that no longer exists.

## Mutation suite

`mutations.py` plants one defect at a time — a wrong last digit, a flipped
sign, an operand moved out from under a certified claim, a caption that no
longer matches its figure, a tagged value dropped from its section, a planted
dead exemption — and requires the checks to catch every one. Run it after any
change to the check logic; a gate that certifies a wrong number is worse than
no gate.
