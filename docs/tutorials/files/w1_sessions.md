# W-1 recorded sessions — what was asked, and what came back

Eight conversations with the Studio assistant, recorded against the live provider
(`anthropic` / `claude-opus-5`) by `tools/assistant_sessions.py`. Every turn was
audited afterwards: the workbook the session left behind was reloaded from disk
and re-solved independently, and each number the assistant reported was compared
against that run. This file is the producer's record — the prompts in order, the
files each session produced, and what the audit found.

The base model is the LEM-8 reinforced slope
(`xslope_reinforced_slope.xlsx`): a 24 ft sand embankment at 1.25:1 with a 2 ft
cohesive face band, a 240 psf crest surcharge, and six geogrid layers at 800
lb/ft. Its published answer is **Spencer FS = 1.587** on a circle centered at
(−5.13, 46.98) with R = 47.26, and Bishop = 1.594.

Every session ran under both budget guards (400k tokens, 15 minutes). The
longest was 294 s.

---

## 1. `w1_build_from_image` — build from the drawing

**Start:** a new, empty project (File → New).
**Attached:** `images/lem08_problem_sketch.png`.

> Build this model. Use the dimensions and properties on the drawing. Unit
> system: US customary (ft, psf, pcf). Add a starting circle and run Spencer with
> a search.

**Files:** `images/w1_build_from_image_1.png` (560 × 5552),
`w1_build_from_image_transcript.md`, `w1_build_from_image_after.xlsx`.
**Cost:** 6 completions, 105,747 in (62,610 cached) / 17,181 out, 254 s.

**Right.** The model is complete and it solves. Both materials carry the right γ,
c and φ. The six geogrid layers are at the right elevations (0, 4, 8, 12, 16, 20
above the toe), each starts *on* the face, each is 20 ft long, and each carries
Tmax = 800 lb/ft with a 4 ft pullout length at both ends — the envelope
breakpoints in the saved file sit exactly 4 ft in from each tip. The surcharge is
240 psf across the 70 ft crest. All three starting circles bottom at or above the
base. Reloading the saved workbook and running Spencer independently reproduces
**FS = 1.6281**, exactly the number the assistant reported.

**Wrong or worth a reader's eye.**

- **The face band was read as 2 ft measured normal to the face, not 2 ft
  horizontally.** That is the whole of the difference from LEM-8. Rebuilding the
  saved model with a 2-ft-horizontal band and changing nothing else gives
  **1.5881** against the assistant's **1.6281** — LEM-8's own 1.587. The
  assistant's reading widens the cohesive band to 3.2 ft measured horizontally
  and buys 2.5% of factor of safety. Both readings are defensible from the
  drawing; only one matches the published model.
- **It widened the section past the stated dimension.** The drawing gives 130 ft
  overall; the assistant built 150 ft, extending the toe flat from 30 ft to 50
  ft, citing the extent rule. Measured: this is worth **nothing at all** (1.6281
  either way). Harmless, but it overrode a dimension the drawing states.
- **It built the geometry as polygons, not profile lines.** The saved file's
  `profile` sheet is empty and the geometry lives on the `polygon` sheet, so a
  reader who opens it in Studio finds Profile lines empty. `max_depth` reads
  blank as a consequence — that is by design (max_depth has no meaning for
  polygon input), so the assistant's `slope_data['max_depth'] = 0.0` was inert.
- **It wrote E = 0 and ν = 0 on both materials.** Harmless for LEM, a singular
  stiffness matrix for any later finite element run.
- **It spent its first two completions dumping engine source** (`inspect.getsource`
  on `xslope.solve`, then 4,000 characters of `xslope.fileio`) looking for the
  reinforcement schema that is already in its instructions. The first of those
  two produced no output at all. Roughly a third of the turn's cost bought
  nothing.
- **Good behavior worth showing:** the preflight MODEL CHECKS caught its
  `adhesion = 0 / delta = 0` (which selects the overburden pullout law with a
  zero friction angle), it read the finding, cleared both to blank, and the
  checks came back clean.

---

## 2. `w1_modify` — three edits in one conversation

**Start:** the completed LEM-8 workbook. One session, three turns, cumulative.

> (a) Change the slope face to 2:1 and rerun the search.
> (b) Add a distributed load of 500 psf on the crest from x = 60 to x = 90 and rerun.
> (c) Extend all the reinforcement lines 5 ft to the right and rerun.

**Files:** `images/w1_modify_1.png` … `_3.png`, `w1_modify_transcript.md`,
`w1_modify_after_1.xlsx`, `_2.xlsx`, `_3.xlsx` (one per turn).
**Cost:** 13 completions, 218,249 in (162,786 cached) / 11,956 out, 210 s.

**This is the session that earns the tutorial its title.**

**Turn (a) did not do what was asked, and said it had.** The assistant rebuilt
`slope_data['polygons']` for the 2:1 face. But this model is profile-line
native, and Studio's automatic resync rebuilds `polygons` *from* `profile_lines`
after every snippet — so the edit was silently discarded and the face stayed at
1.25:1. The saved `w1_modify_after_1.xlsx` still carries
`(0,0) (30,24) (32,24)`. Everything *else* in that turn did apply: the six
reinforcement lines were moved to start at x = 2y, and the 240 psf block was
moved to start at x = 48. The result is an internally inconsistent model — a
1.25:1 face with 2:1 reinforcement hanging inside it — and the MODEL CHECKS
called it clean, because after the revert the geometry was valid.

Measured cost of the failure: re-running the saved file gives **1.2276**, which
is what the assistant reported (1.228). Doing the same edit *correctly* — through
`profile_lines`, with a resync — gives **Bishop 1.948 / Spencer 1.948**. The
reported answer is 37% low.

**It also ran the wrong method.** "Rerun the search" on a model whose published
answer is Spencer got **Bishop**, because `run_lem`'s default method is bishop.
It offered Spencer as a cross-check afterwards. The same thing happens in
session 4.

**Turns (b) and (c) are the recovery, and they are excellent.** Both edits were
applied correctly (the 500 psf block lands as a second `dloads` block from
(60, 24) to (90, 24); every `x2` gains exactly 5 ft). Both reported "FS did not
move" — and rather than accept it, the assistant traced where the critical circle
meets the ground, discovered the geometry mismatch **itself** in turn (b), stated
it plainly, showed in turn (c) that not one reinforcement layer intersects the
critical surface, refused to guess which geometry the user meant, and offered two
specific repairs. Its diagnosis is correct in every particular.

Both "no change" conclusions are artifacts of the broken geometry. On a correctly
laid-back 2:1 face, the 500 psf block moves Bishop from 1.9481 to **1.9157**, and
the 5 ft tails move it to **2.0315**.

**What a reader should verify:** after any geometry edit, open the profile-line
editor (or reload the saved file) and confirm the vertices actually changed.
"Model checks clean" does not mean the edit took.

---

## 3. `w1_sweep_builtin` — a parametric sweep

**Start:** the completed workbook.

> Sweep the geogrid Tmax from 500 to 3000 lb/ft in 6 steps with a search at each
> step and plot FS against Tmax.

**Files:** `images/w1_sweep_builtin_1.png`, `w1_sweep_builtin_transcript.md`. No
workbook — the sweep restored the project, and the harness correctly wrote
nothing.
**Cost:** 2 completions, 26,868 in (25,044 cached) / 1,577 out, 104 s.

**Right, and cheap.** It used the preloaded `sensitivity(values, apply, ...)`
helper with `method='spencer', search=True`, wrote one CSV and one plot, and left
the model at Tmax = 800. Note that `list_params` carries **no reinforcement
refs** at all — Tmax is not a sweepable parameter reference — so the callback
helper is the correct tool here, not `design_sweep(param=...)`. Three of its six
steps were re-run independently and match to four figures:

| Tmax | assistant | independent re-run |
| ---: | ---: | ---: |
| 500 | 1.425 | 1.4252 |
| 1000 | 1.641 | 1.6407 |
| 3000 | 1.681 | 1.6813 |

**The explanation for the plateau is not supported.** The assistant wrote that
past ~2000 lb/ft "the mobilizable force is capped by pullout/anchorage (the
`lp1`/`lp2` embedment terms)". Removing the pullout limit entirely (Lp1 = Lp2 =
0, fully anchored) at Tmax = 2000 raises FS only from 1.681 to **1.715** — 2%,
where the claim predicts the plateau should lift. What actually happens is that
the critical circle migrates: its center moves from (−7.74, 48.06) at Tmax = 500
to (5.69, 48.64) at Tmax = 3000, as the search finds a mechanism that engages
less reinforcement. The numbers are right and the shape of the curve is right;
the mechanism offered for it is a plausible story that a measurement contradicts.

---

## 4. `w1_sweep_adhoc` — a study with no mode behind it

**Start:** the completed workbook.

> Run the analysis with 2, 3, 4, 5 and 6 geogrid layers (removing the top layers
> first), searching each time, and tabulate FS against the number of layers.

**Files:** `images/w1_sweep_adhoc_1.png`, `w1_sweep_adhoc_transcript.md`. No
workbook — the model was restored, and the harness confirms it by writing none.
**Cost:** 3 completions, 40,530 in (37,566 cached) / 1,594 out, 55 s.

**Right.** A hand-written loop with a `try/finally` restore, `plot=False` inside
it, one search per row. Two rows re-run independently, and they match the
assistant's table exactly — FS *and* critical circle:

| layers | assistant | independent re-run |
| ---: | ---: | ---: |
| 3 | 1.274, centre (2.73, 43.07), R 35.80 | 1.2739, (2.73, 43.07), R 35.80 |
| 5 | 1.519, centre (−6.00, 47.74), R 48.11 | 1.5195, (−6.00, 47.74), R 48.11 |
| 6 | 1.594 | 1.5937 — LEM-8's published Bishop |

Its reading of the table is sound: each row is its own search, the 3-layer case
legitimately finds a shallower, smaller-radius mechanism, and FS ≥ 1.5 first
occurs at 5 layers.

**Again Bishop, not Spencer** — "run the analysis" took `run_lem`'s default. The
values are Bishop values throughout, and it said so, but a reader comparing them
against the tutorial's Spencer numbers would be comparing two different methods.

---

## 5. `w1_elastic_fem` — stiffnesses, then a finite element run

**Start:** the completed workbook. Two turns.

> (a) Suggest values of Young's modulus and Poisson's ratio for these materials
> so I can run a finite element analysis, and explain your choice.
> (b) Enter them, build a quadratic mesh at 2 ft, and run the strength reduction
> analysis.

**Files:** `images/w1_elastic_fem_1.png`, `_2.png`,
`w1_elastic_fem_transcript.md`, `w1_elastic_fem_after.xlsx`.
**Cost:** 5 completions, 71,104 in (62,610 cached) / 2,654 out, 175 s.

**Right, and the SSRM number is exact.** It called `suggest_elastic()` and
reported the classifier's answers verbatim (shell → Medium Clay, E = 668,300 psf,
ν = 0.40; base → Dense Sand, E = 2,861,300 psf, ν = 0.30), and it noticed that
both materials already carried a placeholder E = 1,000,000 / ν = 0.3. It then
argued — correctly — that a c = 300, φ = 37° material is a c–φ fill rather than a
clay, and entered ν = 0.30 with E = 1.2 × 10⁶ for the shell instead of the
classifier's clay value, saying so. Rebuilding the mesh from the saved file at
the same settings gives an identical mesh (4,364 nodes / 2,101 elements, tri6 at
2 ft) and an identical **SSRM FS = 1.47265625**.

Its statement that SSRM factor of safety is insensitive to E and ν, and that ν
should stay ≤ 0.45 in a drained run, is correct.

**Worth a reader's eye.** It flagged that it could not tell whether the six
reinforcement lines were contributing to the FEM run and asked the user to
confirm — a good instinct, but one snippet would have answered it. And
`slope_data['mesh']` reports `element_type` as `None`, so its own read-back of
the mesh type came back empty.

---

## 6. `w1_conceptual` — two questions that change nothing

**Start:** the completed workbook. Two turns. No workbook written (correct).

> (a) How does a reliability analysis work in XSLOPE?
> (b) How do I decide standard deviations for a reliability analysis if I only
> have a few tests?

**Files:** `images/w1_conceptual_1.png`, `_2.png`,
`w1_conceptual_transcript.md`.
**Cost:** 4 completions, 56,479 in (50,088 cached) / 4,789 out, 83 s.

**Right.** It grounded the answer in the open model (it printed the `sigma_*`
fields first and told the user every one of them is zero). The Taylor formula it
quotes — σ_F = √Σ(ΔF_i/2)² over 1 + 2N solves — matches `reliability/taylor/`
exactly. Its description of the response-surface engine (quadratic surrogate over
a designed set of solves, then mass sampling) matches the implementation. Turn
(b) is sound practice content: Duncan's three-sigma rule, the range ÷ d₂ estimator
for tiny samples, published COV bands, the 1/√(2(N−1)) uncertainty in s itself,
and the point that `base` has c = 0 by definition rather than by measurement.

**Wrong.** It states as a caveat that "each realization runs the critical-surface
search … so the critical surface can migrate between samples." That is true of
the Taylor engine (1 + 2N searches) and **false of Monte Carlo**, where
`xslope.reliability.reliability_mc` holds the surface fixed by design: *"The
failure surface is never randomized."* Stated where it is, next to the Monte
Carlo description, it is a wrong caveat about the engine the reader is most
likely to run.

**Missed.** The links it gives — `reliability/`, `reliability/taylor/`,
`reliability/monte_carlo/` — are real pages and are correct. But it never called
`corpus_index` and named no worked example, though the corpus has plenty:
VP28, VP29 (Duncan's LASH terminal, TSPM vs Monte Carlo), VP33, VP34, RS2-25
(Syncrude tailings dyke), LEM sample 15 and FEM sample 4. It also never mentioned
`reliability/fem/` or the LEM-11 tutorial.

---

## 7. `w1_diagnose` — a model broken on purpose

**Start:** `w1_diagnose_start.xlsx` — the completed LEM-8 model with three
transcription errors written into it by the producer:

1. material 2 (`base`) φ = **3°** instead of 37 (a dropped digit);
2. the crest surcharge = **2400 psf** instead of 240 (a decimal slip);
3. `max_depth` = **−100** instead of −10 (a dropped digit).

Bishop with a search returns 0.083 on it.

> This model gives a factor of safety below 1. Can you find what is wrong?

**Files:** `images/w1_diagnose_1.png` (560 × 6086),
`w1_diagnose_transcript.md`. No workbook — it restored everything it touched and
said so, and the harness confirms it by writing none.
**Cost:** 9 completions, 151,873 in (112,698 cached) / 14,802 out, 295 s.

**Found, with the right ranking.** It isolated the causes one at a time rather
than guessing, and it led with the real one: `base` at c = 0, φ = 3° "cannot hold
any slope — it's weaker than wet mud", and it correctly read it as a dropped
digit. Its measurements are exactly reproducible: broken model **0.0830**, φ = 30
alone **0.9174** — both match its report to three figures. It also noticed on its
own that φ = 30 does not get the slope back and said the base strength probably
needs more; it does — the true value is 37, which alone gives 1.1836.

**Fault 3 found and correctly de-ranked.** It flagged `max_depth = −100` as "a
lot of invented room", and — importantly — *measured* that it is not what is
driving FS below 1 before saying so. Correct: fixing max_depth changes nothing
(1.5937 either way).

**Fault 2 flagged but not identified.** It called out the 2,400 psf crest load as
"~18.5 ft of soil equivalent" and asked whether it was real or meant as ponded
water — but did not read it as a decimal slip from 240. It is the second real
bug: with φ restored to 37, dropping the surcharge back to 240 psf returns
**1.5937**, LEM-8's published Bishop value exactly.

**One invented bug.** It reports as "the second bug" that the `shell` zone is a
2 ft sliver rather than the whole embankment, and that it "was almost certainly
meant to be the full embankment". That is wrong — the thin cohesive face band is
the model's design and the point of the published problem. The claim is stated
with the same confidence as the true findings.

**One rough edge:** a snippet crashed with `KeyError: 'Xo'` reading a circle
centre off `run_lem`'s result dict, which does not carry one. It recovered on the
next snippet.

---

## 8. `w1_report` — the analysis report

**Start:** the completed workbook. Two turns, because the report documents what
the session has solved: the first turn is the run the second writes up.
`report/finalize` was pinned off for the recording so the unattended capture
could not take over the machine's copy of Word.

> (a) Run Spencer with a search.
> (b) Generate the analysis report for this model.

**Files:** `images/w1_report_1.png`, `_2.png`, `w1_report_transcript.md`,
`w1_report_after.docx`.
**Cost:** 4 completions, 53,779 in (50,088 cached) / 868 out, 39 s.

**Turn (a) is exactly right:** Spencer with a search returns **FS = 1.587** on
the circle (−5.13, 46.98), R = 47.26 — LEM-8's published answer to four figures,
and identical to the independent baseline run.

**Turn (b) produced a report that does not contain it.** The .docx exists, opens,
and renders (checked by converting it to PDF with LibreOffice headless and
reading the pages). It is **three pages**, and its only sections are
*1 Traceability* and *2 Project Definition*. **There is no results section and no
factor of safety anywhere in it.** The assistant nevertheless wrote: "the results
section covers that run — FS = 1.587 on the critical circle at (−5.13, 46.98),
R = 47.26". That claim is false, and a reader who does not open the file will not
learn it.

The cause is in the kernel, not the model: `run_lem` never stores its bundle on
`doc.results['lem_solution']`, though `run_fem` stores `fem_solution` and
`run_seep` stores `seep_solutions`. `MainWindow.report_solutions()` reads
`doc.results['lem_solution']`, finds nothing, and the report is built with no LEM
run to document. A second, smaller gap: the kernel's `generate_report` never
passes `input_path`, so the traceability page reads **"Input file: not saved to a
file"** and carries no SHA-256, even for a project opened from disk.

The report also lands at `<model>_report.docx` beside the project — here that put
`xslope_reinforced_slope_report.docx` into `docs/tutorials/files/`; the producer
renamed it to `w1_report_after.docx`.

---

## Totals

| Session | Turns | Completions | Tokens in (cached) | Tokens out | Wall |
| --- | ---: | ---: | ---: | ---: | ---: |
| build_from_image | 1 | 6 | 105,747 (62,610) | 17,181 | 254 s |
| modify | 3 | 13 | 218,249 (162,786) | 11,956 | 210 s |
| sweep_builtin | 1 | 2 | 26,868 (25,044) | 1,577 | 104 s |
| sweep_adhoc | 1 | 3 | 40,530 (37,566) | 1,594 | 55 s |
| elastic_fem | 2 | 5 | 71,104 (62,610) | 2,654 | 175 s |
| conceptual | 2 | 4 | 56,479 (50,088) | 4,789 | 83 s |
| diagnose | 1 | 9 | 151,873 (112,698) | 14,802 | 295 s |
| report | 2 | 4 | 53,779 (50,088) | 868 | 39 s |
| **total** | **13** | **46** | **724,629 (563,490)** | **55,421** | **1,215 s** |

At `claude-opus-5` list rates ($5.00 / MTok input, $0.50 / MTok cache read,
$25.00 / MTok output), that is roughly **$2.50** for the whole capture — about
$0.81 of uncached input, $0.28 of cache reads and $1.39 of output, plus whatever
the cache writes cost (the accumulator does not break them out).

## The pattern across all eight

Every number the assistant reported was reproducible: eleven separate
FS values, the SSRM factor of safety to eight digits, and two critical-circle
centres all matched an independent re-run exactly. **Not one arithmetic or
solver-level error.**

Every failure was upstream of the arithmetic — in what was modelled, which
method was run, or what was claimed about the output:

- an edit that was silently reverted and reported as done (session 2);
- Bishop run where the model's method is Spencer, twice, from a default
  (sessions 2 and 4);
- a mechanism offered for a curve that a measurement contradicts (session 3);
- a wrong statement about how one of the engines treats the failure surface
  (session 6);
- a confidently-stated bug that is not a bug (session 7);
- a description of a deliverable that the deliverable does not contain
  (session 8).

That is the lesson the page has to teach: the assistant computes correctly and
narrates optimistically. Check the model, check the method, and open the file.
