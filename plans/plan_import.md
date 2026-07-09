# Importing Models from Other Slope-Stability Software — Exploration & Plan

A high-value direction beyond DXF: **importing problems from other slope-stability
software** — GeoStudio **SLOPE/W**, Rocscience **Slide2**, and similar — mapping their
geometry, materials, piezometric/water conditions, and (where sensible) failure
surfaces onto `slope_data`, the same target the DXF importer and the Studio assistant
already populate. This would let users bring existing models into XSlope without
re-drawing them.

**Status:** 🟢 **Committed — validation unblocked, importer still gated on a license.**
The author is acquiring licenses to GeoStudio and Rocscience Slide2 (and similar) to
reverse-engineer the formats, import and re-solve the problems in xslope as a form of
**validation** against the commercial codes' answers, and wire it into the Studio
interface. An initial scout (2026-06-27, §5) confirmed SLOPE/W `.gsz` is parseable and
the strong first target.

A second scout (2026-07-09, §6) found that **both vendors publish free, no-login
verification manuals with reference factor-of-safety values**. That decouples the two
halves of this plan: the **benchmarking** work (Phase B1) needs no license, no vendor
files, and no importer, and can start immediately. Only the **importer** (Phase A) still
waits on file-format access. Extracted from `plan_studio.md` so it can grow on its own.

---

## 1. Goal

Let users open a model authored in another LEM/seepage package and get a populated
`slope_data` — geometry (material zones), materials, water/pore-pressure conditions,
and analysis definitions — plus a list of caveats for anything that doesn't map 1:1.
Surfaced in Studio as **File → Import → <format>**, mirroring the DXF import path
(engine-side parser → `slope_data` + caveats → live document → user reviews → Save As).

## 2. Target formats

- **GeoStudio SLOPE/W (`.gsz`)** — the first target (zipped XML; see findings below).
  Also the de-facto interchange format (Slide2 itself imports `.gsz`).
- **Rocscience Slide2 (`.slmd`)** — harder; undocumented/proprietary, samples not freely
  downloadable.
- Others (PLAXIS, Geostudio SEEP/W for seepage, etc.) — later, demand-driven.

## 3. Unknowns to resolve first

- **File formats & docs.** What does each package's native file look like, and is the
  format documented or reverse-engineerable? Some are XML-ish or text; others are
  opaque/binary or project bundles. SLOPE/W (`.gsz` — a zipped XML project) and Slide
  are the obvious first targets to investigate.
- **Sample files.** We need a corpus of real `.gsz` / Slide files (ideally with known
  answers) to develop and regression-test an importer.
- **Semantic mapping.** Their material models, pore-pressure definitions (Ru, piezo
  lines, pressure grids), and surface conventions don't map 1:1 onto `slope_data`;
  scope what's faithfully importable vs. flagged-for-review.
- **Per-importer module.** Likely one `xslope` importer per format (engine-side, so
  notebooks benefit too), surfaced in Studio as File → Import → <format>, each
  returning a populated `slope_data` + a list of caveats — mirroring the DXF path.

## 4. Approach & phased roadmap

With licensed access to the source software and its tutorial corpus, the work proceeds
in three phases (per format, easiest format — `.gsz` — first):

- **Phase A — Reverse-engineer the format.** Acquire a license, and map the on-disk
  structure to `slope_data`. `.gsz` is zipped XML (tractable — see §5); `.slmd` is
  proprietary/binary (harder — may lean on Slide2's own `.gsz` export as a bridge, or
  reverse-engineer the binary). Prefer learning the schema from **`.gsz` files the author
  authored himself** in a licensed GeoStudio rather than from vendor sample files — see
  §6 "Author your own corpus." Output: a per-format engine-side parser returning a
  populated `slope_data` + a list of caveats (mirroring the DXF path).

- **Phase B1 — Benchmark against published verification manuals. ⭐ No license required;
  start now.** Both vendors publish verification manuals as free, no-login PDFs carrying
  **reference FS values** for a large bank of classical problems (§6). Re-create those
  problems directly as xslope `.xlsx` inputs and compare xslope's FS against the
  published number. This needs no vendor software, no `.gsz` file, and no importer, and
  the resulting corpus is **wholly ours** — the `.xlsx` inputs are our own authored
  content, and the reference numbers are publicly published values we cite rather than
  redistribute. Everything can live in the **public** repo and run in public CI. Capture
  it the way the existing `run_tests.py` benchmarks do, via a markdown
  `<!-- test: ... -->` tag carrying `expected_fs` / `expected_flowrate` with a citation to
  the manual and example number. *This is the external validation deferred elsewhere
  (e.g. [`plan_vg.md`](plan_vg.md) §7) — Slide2's Groundwater Verification Manual supplies
  seepage reference answers, and both slope-stability manuals supply LEM ones.*

- **Phase B2 — Import fidelity (license-gated).** Once the Phase A parser exists, import
  each `.gsz` and check the parsed `slope_data` round-trips faithfully, then re-solve and
  compare. Distinct from B1: B1 validates the **solvers** against published answers; B2
  validates the **parser** against files. B2 needs real `.gsz` files, so it is the only
  part that touches the licensing questions in §4/§6.
- **Phase C — Wire into Studio.** Surface as **File → Import → <format>** with a wizard
  that shows the parsed entities + caveats, populates the live document (renders on the
  canvas for review), and leaves it unsaved for the user to complete and Save As — the
  same pattern as the DXF importer.

**Corpus & licensing — revised 2026-07-09.** The earlier decision ("put vendor `.gsz`
files in the private sibling repo `xslope_private_tests`; private ≠ redistribution")
**does not survive a close read of the Seequent User Terms.** C.2.1 governs downloaded
sample files as "Materials" and requires they be used "solely for your personal and
non-commercial purposes and **not copied, modified or distributed in any way**." The
operative prohibition is *copying*, not *public* distribution — so a private repo that
collaborators or CI can read is not obviously safer than a public one. Revised policy:

- **Never commit a vendor-authored `.gsz`/`.slmd` anywhere**, public or private. Keep any
  such file in a **local, git-ignored** dev folder only.
- **Author your own corpus.** Section A.5.3 provides that content submitted to a Product
  remains owned by the Customer, and A.1's definition makes the author the Customer where
  no written agreement applies. So **models the author builds himself in a licensed
  GeoStudio are his own content** — free of the Materials clause, safe to commit to the
  *public* repo, and a much cleaner basis for learning the schema than RE'ing vendor
  samples. Combined with Phase B1's published reference answers, this makes the whole
  corpus + validation table redistributable and dissolves the private-repo problem.
- `xslope_private_tests` remains available as a fallback for anything that can't be
  authored from scratch, but should end up holding little or nothing for this plan.
- The public repo publishes the **derived validation table** (problem name, summary,
  published reference answer + citation, xslope answer, Δ/tolerance) — our own measured
  numbers and citations to public documents, never vendor files.

**Licensing & compliance — GeoStudio User Terms (reviewed 2026-06-30).** *(Not legal
advice; confirm with counsel / Seequent before public release.)* Key constraints from the
Seequent User Terms v2.1:

- **Use an academic (or paid) license for development — NOT the free trial.** The terms
  define **"trial" as a Tech Preview** (Section A definitions), which triggers Section B.3:
  Tech Preview use is limited to **"testing or evaluation" only** (B.3.1) and imposes
  **confidentiality** on anything derived from it (B.3.3). Building/validating an importer
  under a trial would strain both. An academic license (Section B.1, a normal Product
  license) carries neither restriction. **Decision: pursue an academic license.**
- **Format reverse-engineering — will likely use vendor-provided sample files.** The
  `.gsz` structure is easiest to learn from the example problems shipped with GeoStudio.
  A.2.5(a) prohibits reverse-engineering the *Product* (the software); the `.gsz` *file
  format* is a defensible interop grey area, but RE'ing from license-provided files is
  more exposed than the public-samples-only path (PyGeoStudio's open samples remain a
  fallback). This is exactly what the written-consent email below asks them to bless.
- **Never redistribute — or even copy — vendor files.** Sample `.gsz`/tutorials are
  "Materials" (Section C.2.1) — personal, non-commercial, "not copied, modified or
  distributed in any way," retaining Bentley's copyright attribution. See the revised
  corpus policy above: the private-repo answer was wrong; local git-ignored only.

- **No bulk downloading of the example library — scripted *or* manual.** A.2.5(c) bars
  using "any spider, robot, scraper, or otherwise" to "copy any Product or Online Service,
  or overwrite, monitor, collect or harvest any data, or information or content found in
  any Product or Online Service (**or manually carry out the same**)." Harvesting the
  tutorial corpus is squarely the target. Take individual files by hand, as needed.

- **Don't discuss the file format on the Seequent forums.** C.2.3 forbids copying or
  disclosing information shared through the Forums, and C.3.3 grants Bentley a
  "perpetual, royalty-free, unrestricted licence" over anything posted there. The forum is
  the natural place to ask about `.gsz` internals and the worst place to do it — mining
  forum threads for format details, or posting our own, both cut against these clauses.

- **Keep the consent email a narrow permission request.** A.9.1 (Unsolicited Idea
  Submission) states that unsolicited ideas for "technologies, processes ... or new
  products" sent to Bentley "**will become our property** and you agree that all
  intellectual property rights therein are transferred to us." A request for permission is
  not an idea submission, but a detailed importer roadmap edges toward one. Trim the draft
  below: ask whether reading `.gsz` and publishing a benchmark table is permitted; don't
  enumerate the design or the SLOPE/W → SEEP/W → SIGMA/W product sequence.

- **These User Terms are not the license grant.** B.1.2 conditions all Product use on the
  *Customer Terms*, defined in A.1 as the access agreement + these User Terms + the
  **Product Terms** (https://www.seequent.com/legal-privacy/). The decisive question —
  whether an academic license permits building and open-sourcing an importer — is answered
  in the Product Terms and the academic EULA, **neither of which we have read.** Read
  those before relying on any conclusion here.

- **A free Seequent ID buys Online Services, not a Product.** That is a genuinely better
  footing for the format work: A.2.5(a)'s reverse-engineering ban attaches to *Products*
  (the installed software and its documentation), while a sample file pulled from a portal
  is a *Material* under C.2.1 — which forbids copying/modifying/distributing but says
  nothing about analyzing. It is **not clean**, though: A.2.5(a) also names "Online
  Services" in the RE list, and C.1.3 bars copying or storing "any of the Online Services"
  without prior written consent. Still grey; still gated on the written consent below.

- **Open source ≠ automatically "academic use."** These are orthogonal: "open source"
  governs *our* code's redistribution; "academic license" restricts *who* uses GeoStudio
  and *for what purpose* (typically teaching/research by the licensee). Some academic
  EULAs permit personal coursework/research but **not producing a tool for general public
  distribution** — a scope question independent of the fact that XSlope is free. This is
  the **open item to resolve in writing** before public release.
- **When in doubt, get written consent (action item).** The terms repeatedly gate on
  "prior written consent." Ask whether an academic license permits building +
  open-sourcing a `.gsz` importer and publishing a validation table. A written yes
  converts grey to black-and-white (vendors often approve, since an importer drives format
  adoption). ⚠️ **Status: this question has not actually been asked yet.** What was
  submitted (~2026-07-02, both vendors) was an academic-license **quote request** through
  the web forms — a sales enquiry, silent after a week. Sales cannot answer a scope
  question, and the email below has never been sent. Two separate threads to run:
  **(a)** follow up the quote requests to get a named rep; **(b)** once that rep exists,
  put the permission question to them and ask them to route it to legal, cc
  legal@seequent.com. Do not treat vendor silence on (a) as any signal about (b).

  Revised draft (narrowed per A.9.1 — a permission request, not a design disclosure;
  addressed to the academic-program contact, cc legal@seequent.com):

  > **Subject:** Academic license — permitted scope for file reading and published benchmarks
  >
  > Hello,
  >
  > I'm a university faculty member applying for an academic GeoStudio license, primarily to
  > study the software's capabilities to inform my teaching in geotechnical/slope-stability
  > courses.
  >
  > I also maintain XSlope, a free, open-source Python package for slope-stability and
  > seepage analysis that I originally built for my own course. Before I rely on an academic
  > license for any of that work, I'd like written confirmation on two points of scope:
  >
  > 1. Whether an academic license permits me to **read GeoStudio's saved model files** in my
  >    own separately-written software, and to release that software under an open-source
  >    license. It would contain no GeoStudio source or other Seequent IP, and I would not
  >    redistribute any GeoStudio files or materials.
  > 2. Whether I may **publish a comparison table** of results — for each benchmark problem,
  >    the factor of safety reported in your published Slope Stability Verification Manual,
  >    the value my software computes, and the difference — with appropriate citation to that
  >    manual. I would publish only my own code and my own measured numbers.
  >
  > Could you let me know whether these are permitted and any conditions that would apply, or
  > route this to whoever is best placed to answer?
  >
  > Thank you,
  > [Name, title, institution]

## 5. Findings from an initial scout (2026-06-27)

- **SLOPE/W `.gsz` is the strong first target — confirmed parseable.** A `.gsz` is a
  ZIP holding the model as plain XML (`GSIData` root) + a `mesh_*.ply` + result CSVs.
  The XML maps closely onto `slope_data`: `<Geometries>` → **Points** `(X,Y)` / **Lines**
  / **Regions** (material zones → `polygons` + `mat_id`); `<Materials>` → `<Material>`
  (strength/hydraulic) → `materials`; `<WaterItems>` → pore-pressure/piezo;
  `<StabilityItems>` → slip surfaces; `<Analyses>` → analysis defs. (Verified by
  unzipping `Rapid drawdown.gsz` and walking the XML.)
- **Reference implementation:** [PyGeoStudio](https://github.com/MoiseRousseau/PyGeoStudio)
  — a Python `.gsz` reader/writer to study (or depend on). ⚠️ **No LICENSE file** as of
  the scout — check rights before reusing its code or redistributing its samples.
- **Sample files:** PyGeoStudio's `examples/GeoStudio_files/` has 5 real `.gsz` files,
  incl. `Rapid drawdown.gsz` and `Reinforcement with Anchors.gsz` (both relevant here).
  Official GeoStudio examples (Seequent/GeoSlope) exist but sit behind a Seequent-ID
  login. ⚠️ These are Seequent's example files — keep them in a **git-ignored** dev
  folder, don't commit, until licensing is clear.
- **Slide2 `.slmd`:** harder — samples ship with the install
  (`C:\Users\Public\Documents\Rocscience\Slide2 Examples`), not freely downloadable;
  format undocumented/proprietary. Notably Slide2 *imports* SLOPE/W `.gsz`, reinforcing
  `.gsz` as the de-facto interchange format and first target.
- **Recommendation:** start with `.gsz` → `slope_data` (Regions→polygons, Materials,
  WaterItems). ⚠️ *Superseded in part by §6:* study PyGeoStudio's **parser** for schema
  understanding, but do not vendor its code or redistribute its **samples** — the repo
  carries no license at all (all rights reserved). Prefer spiking against `.gsz` files we
  author ourselves.

## 6. Access routes & the license-free benchmark corpus (scout, 2026-07-09)

Prompted by a week of silence on the academic-license quote requests submitted to both
vendors' web forms (~2026-07-02), and by a close read of the Seequent User Terms v2.1
(findings folded into §4 above). A Seequent/Bentley account now exists; no license yet.

**⭐ The reference answers are free and public — from both vendors.** This is the finding
that unblocks Phase B1. Neither manual requires a login; both carry per-problem reference
FS values that we may **cite** (citing a published number is not redistributing a
"Material"):

- **GeoStudio Slope Stability Verification Manual** (Oct 2022) —
  `https://files.seequent.com/PDFs/Geostudio-Slope Stability Verification Manual-Oct2022.pdf`
- **Stability Modeling with SLOPE/W** engineering book (Oct 2022, 263 pp) —
  `https://files.seequent.com/PDFs/Geostudio-Stability Modeling-Oct2022.pdf`
  (the whole engineering-book series lives under `files.seequent.com/PDFs/`; a SEEP/W
  seepage volume very likely follows the same URL pattern — **unverified**)
- **SLOPE/W Getting Started tutorial** (no login) —
  `https://files.seequent.com/GeoStudio/SlopeW/SLOPE Tutorial.pdf`
- **Rocscience Slide2 verification manuals** (no login) —
  https://www.rocscience.com/help/slide2/verification-theory/verification-manuals —
  Slope Stability (Examples 1–111, incl. hand-calculation benchmarks giving Bishop
  Simplified FS), a Sarma manual (Problems 1–9), and a **Groundwater Verification Manual**
  (Examples 1–21) that supplies our seepage reference answers.

  ⚠️ The two GeoStudio PDFs are image/scanned streams — the FS values need OCR or manual
  transcription. Rocscience's are HTML + PDF and read cleanly. **Start with Rocscience.**

**Access routes (GeoStudio).** Bentley acquired Seequent, so the free academic tier is
delivered through the Bentley Education Program:

- **Seequent Academic Program** — https://www.seequent.com/products-solutions/academic-program/
  Three tiers: **Free** (PLAXIS + GeoStudio, delivered via Bentley Education),
  **Classroom WorkSuite** (paid), **Research** (paid, reduced rate, for faculty/PhD).
  "Apply now" form; applications assessed case-by-case → **this is the route that produces
  a named human to ask.**
- **Bentley Education** — https://www.bentley.com/education/ ·
  sign-up https://www.bentley.com/edu-sign-up/ — free learning licenses to students *and*
  educators; catalog explicitly lists **GeoStudio 2D/3D** and PLAXIS.
- 💡 **The free tier needs no quote and no sales reply.** Bentley Education is a self-serve
  eligibility sign-up, not a quoted sale — so it is not blocked on the pending quote
  requests and could put GeoStudio on the machine this week. Weigh that against the
  restriction below before accepting it.
- ⚠️ **The free tier's restriction may be the tightest of all:** "for educational purposes
  only, not for school operations or **commercially-funded research**." Building a
  generally-distributed importer is plausibly outside "educational purposes." The paid
  **Research** tier's terms (non-commercial, institutional email, findings must be
  published, Seequent acknowledged in publications) read as a *better* fit for this work
  than the free tier. Do not assume the free tier is the safe choice merely because it is
  the most restricted.

**Access routes (Slide2).** No free faculty or student tier exists. The route is the
**Academic Bundle** — https://www.rocscience.com/plans-pricing/academic-bundle — **$1,250/yr**
(or $5,000/5 yr), 20 Rocscience programs incl. Slide2, unlimited users per institution,
educational-use-only. Given that Slide2's verification manual is free and `.slmd` remains
proprietary and low-priority, **there is no near-term reason to buy this.**

**Author your own corpus (the unlock).** See §4. Models authored by us in a licensed
GeoStudio are our content under A.5.3, not vendor "Materials" — making them both the
cleanest basis for reverse-engineering the schema and freely committable to the public
repo. Rebuilding verification-manual problems ourselves therefore yields a corpus where
the input files are ours and the reference answers are published, with no rights encumbrance
on either.

**PyGeoStudio rights — resolved, and the answer is unfavorable.** The repo has
`license: null` and no `LICENSE` file (confirmed via the GitHub API), so its code *and* its
five bundled `.gsz` files are **all rights reserved by default**. Worse, the sample files
appear to be adapted from Seequent's own examples, so their presence in a third-party repo
launders nothing. **Study the parser for schema understanding; do not vendor its code or
redistribute its samples.**

## 7. Connections to other work

- **van Genuchten support** ([`plan_vg.md`](plan_vg.md)) — native vG makes a `.gsz`/Slide
  importer **lossless** for hydraulic functions (their unsaturated curves are vG/Fredlund),
  instead of needing a fitted conversion to the linear-front model.
- **DXF importer** (`plan_studio.md` Phase 6) — the template: an engine-side parser →
  `slope_data` + caveats → live document → user fills placeholders → Save As. A `.gsz`
  importer should mirror that structure (and the Studio File → Import wizard).
- **Soils database** ([`plan_soils_db.md`](plan_soils_db.md)) — imported materials with
  missing properties could be backfilled from the soils-DB presets.
