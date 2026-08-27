# W-1 recorded sessions — what was asked, and what came back

Eight conversations with the Studio assistant, recorded against the live provider
(`anthropic` / `claude-opus-5`) by `tools/assistant_sessions.py`. Every turn was
audited afterwards: the workbook the session left behind was reloaded from disk
and re-solved independently, and each number the assistant reported was compared
against that run. This file is the producer's record — the prompts in order, the
files each session produced, and what the audit found.

Five of the eight — `modify`, `sweep_adhoc`, `conceptual`, `diagnose` and
`report` — were recorded a second time after the kernel changes of `aa63d45f`:
`run_lem` now stores its bundle on `doc.results['lem_solution']`, defaults to the
method the model declares, and returns the surface it solved on; a polygon edit
on a profile-line model is warned about; `generate_report` passes the input path;
and the chat renders markdown. Their prompts are unchanged. Each of those five
sections ends with **what changed** against the first recording. The other three
— `build_from_image`, `sweep_builtin`, `elastic_fem` — are the first recording
and stand as they were.

The base model is the LEM-8 reinforced slope
(`xslope_reinforced_slope.xlsx`): a 24 ft sand embankment at 1.25:1 with a 2 ft
cohesive face band, a 240 psf crest surcharge, and six geogrid layers at 800
lb/ft. Its published answer is **Spencer FS = 1.587** on a circle centered at
(−5.13, 46.98) with R = 47.26, and Bishop = 1.594. The model declares no LEM
method of its own, so `run_lem` called without one now runs Spencer.

Every session ran under both budget guards (400k tokens, 15 minutes). The
longest was 254 s.

---

## 1. `w1_build_from_image` — build from the drawing

**Start:** a new, empty project (File → New).
**Attached:** `images/lem08_problem_sketch.png`.

> Build this model. Use the dimensions and properties on the drawing. Unit
> system: US customary (ft, psf, pcf). Add a starting circle and run Spencer with
> a search.

**Files:** `images/w1_build_from_image_1.png` (560 × 5552),
`w1_build_from_image_transcript.txt`, `w1_build_from_image_after.xlsx`.
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

**Files:** `images/w1_modify_1.png` … `_3.png`, `w1_modify_transcript.txt`,
`w1_modify_after_1.xlsx`, `_2.xlsx`, `_3.xlsx` (one per turn).
**Cost:** 12 completions, 215,159 in (163,848 cached) / 7,243 out, 187 s.

**All three edits took, and all three numbers reproduce.** Every workbook was
reloaded from disk and re-searched independently with Spencer:

| Turn | assistant | independent re-run | critical circle |
| --- | ---: | ---: | :--- |
| (a) 2:1 face | 1.948 | 1.9479 | (7.39, 59.83), R = 61.14 |
| (b) + 500 psf | 1.916 | 1.9158 | (2.69, 84.39), R = 84.39 |
| (c) + 5 ft tails | 2.029 | 2.0292 | (8.65, 72.99), R = 74.80 |

**Turn (a) is the one to show a reader, because it corrects itself.** The
assistant's first move was to rebuild `slope_data['polygons']` for the 2:1 face —
the same move that was silently discarded in the first recording. The snippet's
result now carries the line

> `WARNING: polygons were edited on a profile-line model and have been rebuilt
> from profile_lines; edit profile_lines instead …`

It read that, printed `profile_lines` back out, rewrote them
(`(0,0) (48,24) (50,24)` for the shell band and `(−30,0) (2,0) (50,24) (100,24)`
for the base), resynced, printed the vertices again to confirm, and only then
ran. The saved `w1_modify_after_1.xlsx` carries the 2:1 face.

Everything the edit implies was carried with it: the six reinforcement layers
moved to start at x = 2y on the new face, each still 20 ft long with 4 ft of
pullout at both ends; the 240 psf block moved to start at the new crest break
x = 48; the toe starting circle was rebuilt through the toe (Xo = 24, Yo = 48,
R = 53.67) and the deep circle recentred. It stated all of it in a before/after
table.

**It ran Spencer**, without being told to and without naming a method, because
that is what the model resolves to — and it said so: "Spencer (the method the
model declares)".

**Turns (b) and (c) are ordinary competent work.** The 500 psf block lands as a
second `dloads` entry from (60, 24) to (90, 24); the assistant noticed on its own
that it overlaps the existing 240 psf block, said the crest therefore carries
740 psf over x = 60–90, and asked whether replacement was meant instead. In turn
(c) every `x2` gains exactly 5 ft and the pullout breakpoints stay 4 ft from each
tip. It reported the Spencer admissibility note it got on the last run (line of
thrust outside the slice on 13% of boundaries) rather than dropping it.

**What changed vs. the first run.** Everything that was wrong here is gone. The
first recording rebuilt the polygons on a profile-line model, had the edit
reverted with nothing said, reported **1.228** for a slope that was still at
1.25:1, and ran Bishop by default; turns (b) and (c) then reported "FS did not
move" and spent themselves diagnosing the mismatch. This run edits the right
source, reports 1.948 — the value the first audit measured as the correct
answer — and the following two turns measure real changes instead. Cost is
about the same (12 completions against 13), and the wall time is shorter.

---

## 3. `w1_sweep_builtin` — a parametric sweep

**Start:** the completed workbook.

> Sweep the geogrid Tmax from 500 to 3000 lb/ft in 6 steps with a search at each
> step and plot FS against Tmax.

**Files:** `images/w1_sweep_builtin_1.png`, `w1_sweep_builtin_transcript.txt`. No
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

**Files:** `images/w1_sweep_adhoc_1.png`, `w1_sweep_adhoc_transcript.txt`. No
workbook — the model was restored, and the harness confirms it by writing none.
**Cost:** 7 completions, 114,999 in (95,578 cached) / 5,603 out, 140 s.

**Right, and now in the model's own method.** A hand-written loop with a
`try/finally` restore, `plot=False` inside it, one search per row, and
`ensure_reinforce_pullout` / `build_reinforce_lines` called after each edit so the
pullout envelopes are rebuilt for the layers that remain. Three rows were re-run
independently — FS *and* critical circle — and match exactly:

| layers | assistant | independent re-run |
| ---: | ---: | ---: |
| 2 | 1.210, center (−5.71, 48.23), R 45.51 | 1.2097, (−5.71, 48.23), R 45.51 |
| 3 | 1.272, center (2.73, 43.07), R 35.80 | 1.2722, (2.73, 43.07), R 35.80 |
| 5 | 1.512, center (−6.00, 47.74), R 48.11 | 1.5117, (−6.00, 47.74), R 48.11 |
| 6 | 1.587 | 1.5867 — LEM-8's published Spencer |

**Its mechanism claim survives the measurement.** It reads the 3→4 step as the
point where the critical surface changes character: with 2–3 layers the minimum
is a circle bottoming at el. +2.7 and +7.3, above the reinforcement, so the grids
are bypassed. Tested by deleting the reinforcement and re-solving *on the same
circles*: FS is identical to four figures either way, and the reinforcement force
summed over the slices is exactly zero. Those circles do bypass the grids. What
the grids buy at 2–3 layers is the mechanism itself — with every layer removed
the search drops to **1.1674**, so the reinforcement is holding the deeper
surface shut and forcing the minimum up above itself.

**One unsupported sentence.** It writes that "each added layer buys progressively
less" and then lists increments that say otherwise: 2→3 buys 0.062, 3→4 buys
0.161, then 0.079 and 0.075. The largest step is the fourth layer, which is the
step its own next paragraph explains. The claim and the evidence are in adjacent
sentences.

**What changed vs. the first run.** The method is right — the first recording ran
Bishop from `run_lem`'s old default and produced a table a reader would have
compared against the tutorial's Spencer numbers. The rows are the same shape and
the readings are the same. What got worse is the cost: 7 completions and 115k
tokens against 3 and 41k, because four of the seven went on schema exploration —
`inspect.getsource`, `dir()`, `pkgutil.iter_modules`, `inspect.signature` — before
any analysis ran, including one snippet that raised `AttributeError: module
'xslope' has no attribute 'solve'`. The brief tells it not to open a turn that
way. This is the same waste as session 1.

---

## 5. `w1_elastic_fem` — stiffnesses, then a finite element run

**Start:** the completed workbook. Two turns.

> (a) Suggest values of Young's modulus and Poisson's ratio for these materials
> so I can run a finite element analysis, and explain your choice.
> (b) Enter them, build a quadratic mesh at 2 ft, and run the strength reduction
> analysis.

**Files:** `images/w1_elastic_fem_1.png`, `_2.png`,
`w1_elastic_fem_transcript.txt`, `w1_elastic_fem_after.xlsx`.
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
`w1_conceptual_transcript.txt`.
**Cost:** 3 completions, 44,084 in (40,962 cached) / 3,315 out, 58 s.

**Right on the engines.** It grounded the answer in the open model (it printed
the `sigma_*` fields first and told the user every one of them is zero). Its
three-engine table matches the implementation, including the distinction the
first recording got wrong: Taylor re-searches on every one of its 1 + 2N solves,
Monte Carlo and the response surface hold the failure surface fixed. Every helper
it names is real and correctly spelled — `reliability_taylor`, `reliability_mc`,
`reliability_rs`, the `reliability(engine=…)` front door, and
`parametric_sweep(plot='tornado')` in turn (b). Turn (b) is sound practice
content: Duncan's three-sigma rule with the /4 caution, published COV bands,
sample `s` as a lower bound with few tests, spatial averaging, and the negative
c′–φ′ correlation the Taylor engine does not model.

**Wrong: the reliability index it defines is not the one xslope reports.** It
gives β = (E[FS] − 1)/σ_FS and works a micro-example to β = 2.50, P_f ≈ 0.6% for
E[FS] = 1.45, σ_FS = 0.18. `reliability_taylor` returns only the **lognormal**
index, `beta_ln`; for those same two numbers it prints **β = 2.943 and
P_f = 0.16%**. A reader who follows the answer and then runs the tool gets a
different number for the same inputs.

**Wrong: what an all-zero-σ model does.** It says a run "would return σ_FS = 0
and an infinite β". Run on this model, `reliability_taylor` returns no result at
all — it refuses with a message naming the blank σ columns and the sheet they
live on.

**Still missed.** The links it gives — `reliability/`, `reliability/taylor/`,
`reliability/monte_carlo/` — are real pages and are correct. But it never called
`corpus_index` and named no worked example, though the corpus has plenty:
VP28, VP29 (Duncan's LASH terminal, TSPM vs Monte Carlo), VP33, VP34, RS2-25
(Syncrude tailings dyke), LEM sample 15 and FEM sample 4. It also never mentioned
`reliability/fem/` or the LEM-11 tutorial.

**What changed vs. the first run.** The Monte Carlo error is gone — the brief now
states the surface treatment of each engine, and the answer states it correctly
and calls it out as reportable. The corpus is still not consulted, unprompted. In
its place are two new factual slips, both about numbers it did not run: the β
definition and the zero-σ behavior. It is a shorter, denser answer (3 completions
against 4) and it costs less.

---

## 7. `w1_diagnose` — a model broken on purpose

**Start:** `w1_diagnose_start.xlsx` — the completed LEM-8 model with three
transcription errors written into it by the producer:

1. material 2 (`base`) φ = **3°** instead of 37 (a dropped digit);
2. the crest surcharge = **2400 psf** instead of 240 (a decimal slip);
3. `max_depth` = **−100** instead of −10 (a dropped digit).

Spencer with a search returns 0.071 on it.

> This model gives a factor of safety below 1. Can you find what is wrong?

**Files:** `images/w1_diagnose_1.png` (560 × 2950),
`w1_diagnose_transcript.txt`. No workbook — it changed nothing and said so, and
the harness confirms it by writing none.
**Cost:** 3 completions, 47,432 in (40,962 cached) / 5,764 out, 129 s.

**Its measurements are exact and its conclusion is wrong.** The search it ran
reproduces to four figures — **FS = 0.0709** on the circle (−3.09, 49.09) with
R = 45.93 — as do the zone areas it quotes (48 ft² and 14,992 ft²).

**Fault 1 was missed.** It read `base` as c = 0, φ = 3°, computed the
infinite-slope ratio tan 3° / tan 38.7° = 0.066, matched it against the FS it had
just measured, and concluded: "That is not a numerical problem; the model really
is that weak." The planted dropped digit was printed on screen and accepted as
the design. It is the fault that matters: restoring φ = 37 alone takes Spencer
from **0.0709 to 1.1805**, and restoring the 240 psf surcharge as well gives
**1.5867**, the published Spencer answer exactly.

**The invented bug of the first recording is now the headline.** Having accepted
φ = 3°, it needed a culprit and named the geometry: the `shell` zone is "a
degenerate sliver", the embankment body is therefore built out of foundation
clay, and the shell polygon "should be" the whole fill above the toe — toe (0,0)
→ crest (30,24) → (100,24) → (100,0). It offered to redraw it. The 2 ft cohesive
face band is the model's design and the point of the published problem, and the
proposed redraw would replace the problem with a different one.

**Fault 2 flagged, not identified.** It called the 2,400 psf crest load "large —
the equivalent of 38.5 ft of water" and asked "intentional surcharge, or a units
slip?" It offers the right reading without committing to it, and leaves it as a
question rather than measuring it. Measured: on the broken model the surcharge
alone is worth 0.0709 → 0.0964; with φ restored it is worth 1.1805 → 1.5867.

**Fault 3 flagged and correctly de-ranked**, on the argument that the critical
circle bottoms at el. +3.2 and so cannot be reaching a base at −100. That is
right — max_depth is inert here (1.5867 either way once the other two are fixed)
— but it is an argument, not the measurement the first recording made.

**One claim worth checking, and it holds.** It says the six geogrid layers "do
essentially nothing" in a φ = 3° soil. On its own critical circle, FS is 0.0709
with the reinforcement and 0.0605 without.

**What changed vs. the first run — this is the session that got worse.** The
first recording isolated the faults one at a time over nine completions, led with
the real one, read it as a dropped digit, and measured that φ = 30 alone still
does not recover the slope. This run spent three completions — two reads and a
narrative — never varied a single input, missed the controlling fault entirely,
and promoted the first run's one invented bug to the answer. It is also the only
session whose failure is more confident than its predecessor's: the first
recording called the sliver "the second bug"; this one calls it the cause and
titles the answer with it.

**A rendering defect this session exposes.** Now that the chat renders markdown,
the assistant reaches for LaTeX, and the dock does not render math: the display
equation comes through as literal
`$$\frac{\tan\varphi}{\tan\beta} = \frac{\tan 3°}{\tan 38.7°} = 0.066$$` in the
figure. Qt's markdown is the GitHub dialect — headings, tables, lists and fenced
code — and nothing else.

---

## 8. `w1_report` — the analysis report

**Start:** the completed workbook. Two turns, because the report documents what
the session has solved: the first turn is the run the second writes up.
`report/finalize` was pinned off for the recording so the unattended capture
could not take over the machine's copy of Word.

> (a) Run Spencer with a search.
> (b) Generate the analysis report for this model.

**Files:** `images/w1_report_1.png`, `_2.png`, `w1_report_transcript.txt`,
`w1_report_after.docx`.
**Cost:** 4 completions, 59,170 in (40,962 cached) / 1,059 out, 50 s.

**Turn (a) is exactly right:** Spencer with a search returns **FS = 1.587** on
the circle (−5.13, 46.98), R = 47.26 — LEM-8's published answer to four figures,
and identical to the independent baseline run. It read the circle straight off
the result dict (`Xo`, `Yo`, `R`, `Depth`, `x_entry`, `x_exit`) and laid the six
numbers out in a small table.

**Turn (b) produced a report that contains it.** The .docx is **12 pages** and
its contents are *1 Traceability*, *2 Project Definition* and *3 Limit
Equilibrium Analysis*, with 3.1 Analysis Inputs, 3.2 Materials, 3.3 Loads, 3.4
Reinforcement and 3.5 Spencer's Method — the last split into the search, the
results, the slice table and the calculations, six figures in all. Section 3.5.2
reads "Spencer's Method gives a factor of safety of 1.587 on the critical surface
of Figure 4." Checked two ways: by converting the file to PDF with LibreOffice
headless and reading the pages, and by extracting `word/document.xml` and
grepping it.

**Traceability now identifies the file.** The stamp reads
`xslope_reinforced_slope.xlsx` with SHA-256
`43bb96b090d647c594102c40142e4eed9acd181a6f1802318b132e0e8b412bb3`, which is the
digest of `docs/tutorials/files/xslope_reinforced_slope.xlsx`.

Everything the assistant said about the deliverable is true of the deliverable:
6 figures, the model as built, and the one analysis solved this session.

**Two rough edges, neither the assistant's.** The report still lands at
`<model>_report.docx` beside the project — here that put
`xslope_reinforced_slope_report.docx` into `docs/tutorials/files/`, and the
producer renamed it to `w1_report_after.docx`. And inside the report, the search
summary says "Trial surfaces evaluated **96**" while the same run's own console
line says it evaluated **838** trial surfaces. Both counts are real and they
count different things: 96 is the number of grid centers held in the search's
`fs_cache`, 838 the number of trial surfaces actually solved across all depths at
those centers. The report's label names the wrong one
(`xslope/report.py`, `len(fs_cache)`).

**What changed vs. the first run.** This is the failure the changes were made
for. The first recording produced a three-page report whose only sections were
Traceability and Project Definition — no results, no factor of safety anywhere in
it — while the assistant wrote that "the results section covers that run —
FS = 1.587". `run_lem` now stores its bundle where `report_solutions()` looks for
it, so the run the session made is the run the report documents; and
`generate_report` now passes the input path, so the traceability page names the
file and its digest instead of reading "not saved to a file". Same two prompts,
same cost, a report four times as long that says what the assistant says it says.

---

## Totals

| Session | Turns | Completions | Tokens in (cached) | Tokens out | Wall |
| --- | ---: | ---: | ---: | ---: | ---: |
| build_from_image | 1 | 6 | 105,747 (62,610) | 17,181 | 254 s |
| modify | 3 | 12 | 215,159 (163,848) | 7,243 | 187 s |
| sweep_builtin | 1 | 2 | 26,868 (25,044) | 1,577 | 104 s |
| sweep_adhoc | 1 | 7 | 114,999 (95,578) | 5,603 | 140 s |
| elastic_fem | 2 | 5 | 71,104 (62,610) | 2,654 | 175 s |
| conceptual | 2 | 3 | 44,084 (40,962) | 3,315 | 58 s |
| diagnose | 1 | 3 | 47,432 (40,962) | 5,764 | 129 s |
| report | 2 | 4 | 59,170 (40,962) | 1,059 | 50 s |
| **total** | **13** | **42** | **684,563 (532,576)** | **44,396** | **1,097 s** |

At `claude-opus-5` list rates ($5.00 / MTok input, $0.50 / MTok cache read,
$25.00 / MTok output), that is roughly **$2.14** for the whole capture — about
$0.76 of uncached input, $0.27 of cache reads and $1.11 of output, plus whatever
the cache writes cost (the accumulator does not break them out). The five
re-recorded sessions are about $1.26 of that.

## The pattern across all eight

Every number the assistant reported was reproducible. Across the eight sessions
that is thirteen separate factors of safety, the SSRM factor of safety to eight
digits, and eight critical-circle centers, every one of them matching an
independent re-run. **Not one arithmetic or solver-level error.**

Every failure was upstream of the arithmetic — in what was modelled, which
method was run, or what was claimed about the output. Four of the six failures
the first recording found are gone, and they were the four the tool could fix:

- the edit that was silently reverted and reported as done is now warned about,
  and the assistant reads the warning and repairs the edit itself (session 2);
- Bishop-where-the-model-means-Spencer is gone from both sessions it appeared in
  (2 and 4);
- the wrong statement about how Monte Carlo treats the failure surface is
  corrected (session 6);
- the report that did not contain the result it was said to contain now contains
  it (session 8).

What is left is the class no interface change reaches — a confident account of
something the assistant did not measure:

- a mechanism offered for a curve that a measurement contradicts (session 3);
- a summary sentence its own numbers contradict (session 4);
- a reliability index defined the textbook way rather than the way the engine
  reports it, and a prediction of what a zero-σ run does that the run does not do
  (session 6);
- a planted dropped digit read as the design, and a bug invented to explain the
  consequence (session 7).

Session 7 is the one to read closely. It is cheaper than its first recording by a
factor of three and worse by every other measure, and the difference is entirely
that it stopped varying inputs: it read the model twice, ran one search, and
reasoned. The first recording changed one thing at a time and got the answer.

That is the lesson the page has to teach: the assistant computes correctly and
narrates optimistically. Check the model, check the method, and open the file.
