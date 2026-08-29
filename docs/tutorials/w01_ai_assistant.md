---
title: "Tutorial W-1 — The AI Assistant"
description: "Studio's built-in assistant put through eight conversations on one layered slope: building the model from a drawing, editing its face, its water and its strengths, two kinds of parameter study, stiffnesses for a finite element run, two questions about the analysis, a diagnosis of a broken file, and the analysis report — with what to check after each one."
---

# Tutorial W-1 — The AI Assistant

Studio carries an assistant that does the work the rest of this series does by
hand. It can build a model from a drawing or from a description in words, edit
any part of one afterwards, run every analysis XSLOPE offers, sweep a parameter
across a range, answer a question out of the built-in documentation, and write
the analysis report — all in the open project, on the same undo stack the
editors use.

In this tutorial we first set the assistant up — choosing a provider and a
model, entering a key, and seeing what a conversation costs to run. Then we work
through a series of sample use cases on one slope: building the model from a
drawing, modifying it, running two kinds of sweep, crossing to the finite
element engine, asking two questions, diagnosing a broken file, and writing the
report. Each one closes on the checking it deserves — what to look at, and what
it should say.

<div class="tut-glance" markdown>
<div class="tgt-row">
<div class="tgt-tile"><span class="tg-label">Analysis</span><p>Any</p></div>
<div class="tgt-tile"><span class="tg-label">Read &amp; try</span><p>~30 min</p></div>
</div>
<div class="tgm-obj" markdown>
**Objectives** — Learn how to set the assistant up and what to ask it for, and
learn how to check each kind of answer it gives back.
</div>
<p><span class="tg-pill">AI assistant</span><span class="tg-pill">provider setup</span><span class="tg-pill">images</span><span class="tg-pill">model edits</span><span class="tg-pill">parameter studies</span><span class="tg-pill">strength reduction</span><span class="tg-pill">documentation</span><span class="tg-pill">analysis report</span></p>
<div class="tgm-model" markdown>**Model** — the layered slope of
[LEM-3](lem03_layered_slope.md) —
[xslope_simple_mult_layers.xlsx](../lem/files/xslope_simple_mult_layers.xlsx)</div>
</div>

---

## Setting it up

The assistant lives in a dock on the right side of the Studio window, titled
**Assistant**. It opens with the window; if it has been closed, **View →
Assistant** brings it back, and like every Studio dock it can be dragged to
another edge, floated, or widened by dragging its inner border.

![The Assistant dock as it opens, with the Files, New chat and Settings buttons, the provider and model caption, the transcript, and the input box](images/w01_chat_dock.png){width=430}

Three buttons sit across the top. **Files…** opens the folder where the
assistant saves the plots and files it generates, **New chat** starts a fresh
conversation, and **Settings…** opens the dialog below. Under them is a caption
naming the provider and model in use — *Claude (Anthropic) · claude-opus-5* in
the figure. The transcript fills the middle, the input box is below it, and
**Send** and **Stop** are at the bottom.

### The settings dialog

You bring your own model: the assistant does nothing at all until it has a
provider and credentials. To supply them, click **Settings…**.

![The Assistant settings dialog on a fresh install: Provider set to Claude (Anthropic), Model to claude-opus-5, an empty API key field, a disabled Base URL field, and Confirm before running code unchecked](images/w01_settings.png){width=460}

The dialog has four fields and a checkbox:

**Provider** — the services Studio can talk to. This release offers
**Claude (Anthropic)**, **OpenAI**, **Kimi (Moonshot AI)** and **Z.ai (GLM)**,
which are hosted and bill for what you use, and **Ollama**, which runs a model
on your own machine, free, and sends nothing anywhere; the list may change from
release to release. Where a provider's catalogue mixes models that read images
with text-only ones — OpenAI, Kimi, Z.ai and Ollama — only the models that can
read an image are listed, and the caption under the box says so; every Claude
model on the list reads images already. Handing the assistant a photograph or a
sketch of a cross section is the first thing this tutorial does.

**Model** — a list of what the chosen provider currently offers, with one marked
as recommended. The box also accepts free text, so a model id released after
your copy of Studio was built can be typed in. **Refresh** beside it asks the
provider for its list again; until a key is entered there is nothing to ask
with, and the caption under the box says the list is the one this version
shipped with.

**API key** — the key from your provider's account page, pasted in. It is masked
as you type and stored in your operating system's keychain, not in a settings
file and not in plain text. Ollama needs no key, and the field is disabled when
Ollama is selected.

**Base URL** — the address the provider is reached at, filled in for you.
Leave it alone unless your provider gave you a different endpoint — a Z.ai
coding-plan key, or an Ollama running on another port or another machine. It
is disabled for Claude and OpenAI, whose addresses are fixed.

**Confirm before running code** — unchecked by default. The assistant works by
writing small Python snippets and running them in Studio's own process; with
this box checked, every snippet is shown to you and waits for your approval
before it executes. Leaving it unchecked is fine for ordinary use: the
assistant's edits change the model in memory, and **Edit → Undo** or closing
without saving takes them back.

Below the fields, two captions report where the model list came from and what
the selected model can do — whether it supports running code, whether it can
read an image, and whether the provider caches the prompt.

### What it costs

A hosted provider bills per token, and both the tokens read and the tokens
written count. A routine question — one asked of a built model, answered with a
run — costs about **33,000 input tokens** on Claude Opus 5, nearly all of them
served from the provider's prompt cache at a tenth the price of fresh input, for
a few hundred to a few thousand tokens of reply. A request that has to measure
something costs several times that, because it is many runs: the diagnosis below
spends 108,000 across six calls.

Models differ in both price and capability, and the two do not track each other
perfectly: a cheaper model may answer a question about the documentation as well
as an expensive one and still fail at a long autonomous build. The line under
the input box reports the tokens each turn spent and the running total for the
conversation:

```
this turn: 52,241 in (47,448 cached) / 2,750 out · session: 52,241 in (47,448 cached) / 2,750 out
```

One request is usually several calls to the model, so that first number climbs
while the assistant works. **New chat** starts the count over.

The eight conversations below cost this much between them, measured as they ran:

| Session | Turns | Model calls | Tokens in (cached) | Tokens out | Wall |
| --- | :---: | :---: | :---: | :---: | :---: |
| Building a model from a drawing | 1 | 3 | 52,241 (47,448) | 2,750 | 47 s |
| Modifying the model | 3 | 13 | 245,907 (205,608) | 7,160 | 137 s |
| A sweep with the helper | 1 | 2 | 33,157 (31,632) | 1,308 | 65 s |
| A sweep written ad hoc | 1 | 4 | 68,299 (63,264) | 2,587 | 55 s |
| Stiffnesses and a strength reduction run | 2 | 6 | 104,870 (63,264) | 2,182 | 576 s |
| Two documentation questions | 2 | 4 | 71,849 (47,589) | 5,065 | 86 s |
| A broken file | 1 | 6 | 108,379 (94,896) | 6,994 | 133 s |
| Generating the report | 2 | 4 | 67,460 (63,264) | 954 | 32 s |
| **Total** | **13** | **42** | **752,162 (616,965)** | **29,000** | **1,130 s** |

Thirteen turns, 42 calls to the model, about 752,000 input tokens of which some
617,000 came from the prompt cache, and 29,000 output tokens: **$1.71** at
Anthropic's list prices on 2026-08-29 ($5.00 per million input tokens, a cache
read at a tenth of that, $25.00 per million output tokens). Wall time and price
do not track each other: the strength reduction run takes half the total time for
under a fifth of the cost, because 528 of its 576 seconds are spent inside the
solver, where no tokens are read or written.

The same eight sessions were then played, unchanged, on other models and
scored against the answers above. **Everything in this table is a snapshot: the
models as they behaved and the list prices as they stood in late August 2026,
when this tutorial was written.** Models are replaced and retrained and prices
are cut or restructured every few months, so treat the numbers as an example of
how to compare models, not as a standing verdict on any of them — the sessions
are easy to repeat on whatever is current.

Each session is scored twice, because an answer goes wrong in two different
ways and they are counted apart. **Model and numbers** asks whether the work was done: the
model it built or edited is the one the request describes, and every factor of
safety it reported is right — for the broken file, whether it found the faults
and its repair returns the published answer. **Explanation** asks whether what
it said about the mechanism, the finding or the reason is right — engineering
content only; a remark about a document's layout is noted where it occurs and
not scored. An answer counted wrong there is usually mostly right with one
sentence that is not, and the face-slope sweep below is that case in full.

| Model | Calls | Tokens in (cached) | Tokens out | Cost | Model and numbers | Explanation |
| --- | :---: | :---: | :---: | :---: | :---: | :---: |
| Claude Opus 5 | 42 | 752,162 (616,965) | 29,000 | $1.71 | 8 of 8 | 7 of 8 |
| Kimi K3 (Moonshot AI) | 54 | 727,011 (603,482) | 72,502 | $1.64 | 7 of 8 | 5 of 8 |
| OpenAI gpt-5.5 | 33 | 414,968 (350,592) | 23,237 | $1.19 | 6 of 8 | 6 of 8 |
| Claude Sonnet 5 | 51 | 1,004,052 (763,043) | 53,449 | $1.17 | 7 of 8 | 7 of 8 |
| GLM-5V-Turbo (Z.ai) | 46 | 564,440 (536,661) | 16,563 | ~$0.23 | 6 of 8 | 3 of 8 |

Opus 5 is the only one right on every model and every number, and both of its
misses are in the second column: the face-slope sweep below, where the table is
right and the sentence explaining it is not, and the report, whose summary
describes a document it did not write. Kimi K3 built and edited the model
exactly as Opus did and was the only one of the other four to test its answer
to the water-table question by running the model, but took three times as long
and loses the second column five times over — mechanisms named that its own
output contradicts, and the 100 ft base depth in the broken file examined and
pronounced legitimate. gpt-5.5 was the leanest run and exact on every number it
computed, which is what makes its two misses in the first column stand out: it
built the model and then stopped to ask permission before running the search,
and in the broken file it fixed the wrong thing — it deleted the foundation by
moving the base to elevation 0 and reported 2.396, nearly twice the answer.
Sonnet 5 matched Opus on every number it computed, including all three edits;
the 100 ft base depth is the one thing it printed and never came back to, and
its one wrong sentence puts the rise under submergence down to pore-pressure
relief, which at φ = 0 reaches nothing. GLM-5V-Turbo, a fifth of the cost of
the next cheapest, put the top of the foundation at the rock instead of the
ground and never noticed the model had lost a layer, and in the broken file it
invented a unit weight rather than reading the one the tutorial gives; it calls
the same contact circle deep in three separate sessions, which is what leaves
it at 2. Z.ai publishes no price for that model; its cost here is at the rate a
reseller lists for it.

Prices change; the current ones are on each provider's pricing page:
[Anthropic](https://www.anthropic.com/pricing),
[OpenAI](https://openai.com/api/pricing/),
[Moonshot AI](https://platform.moonshot.ai/docs/pricing) and
[Z.ai](https://docs.z.ai/guides/overview/pricing); Ollama is free. The turns,
calls, tokens and wall times behind the table above are listed session by
session in [w1_sessions.txt](files/w1_sessions.txt).

---

## What to expect

Everything below was produced with **Claude Opus 5**. Another model will do the
same work differently, and so will the same model on another day — the wording
will differ, the code it writes will differ, and it may take a different number
of steps to arrive at the same place. Treat the transcripts as an illustration
of what the requests look like and what comes back, not as a script to match
line for line.

The assistant is good at this work, but everything it produces should still be
checked. Checking is usually easy, and it is the same checking a model built by
hand deserves: read the model checks Studio runs after every edit, look at the
section on the canvas, and do one hand calculation against the number that came
back. Where the five models went wrong most often: a layer boundary read off the
drawing at the wrong elevation, a mechanism described without the surface ever
being extracted, and a broken file "repaired" by changing the wrong input — the
model checks stay clean through all three, because none of them makes a model
inconsistent, only wrong.

---

## Building a model from a drawing

The assistant can build a model straight from a drawing, which is the request
that covers the most ground at once. We will paste the problem drawing from
[LEM-3](lem03_layered_slope.md) into an empty project and ask for the model it
describes, with nothing to go on but the units and where the section should end.
One request covers reading dimensions off a picture, entering two materials and
two profile lines, setting the rigid base, generating starting circles, and
searching for the critical surface.

We start from **File → New**. Right-clicking the drawing at the top of that page
and choosing **Copy image** puts it on the clipboard; from there we click in the
chat box, press Ctrl/Cmd+V to attach it, and type:

<div class="prompt-block" markdown>
```text
Build this model. Unit system: US customary (ft, psf, pcf). Both soils are undrained. Put the toe at x = 0, run the ground 30 ft past the toe and 50 ft behind the crest break, and treat the rock as rigid. Add starting circles and run Spencer with a search.
```
</div>

![The whole build conversation in the dock: the attached drawing, the snippet that builds the model, the three generated starting circles, the clean model checks, the searched Spencer run, and the closing summary of geometry, materials and result](images/w1_build_from_image_1.png){width=560}

One snippet builds the model, and every dimension in it comes off the drawing or
out of the request. The toe sits at (0, 0), the 2:1 face runs up 20 ft to the
crest break at (40, 20), and the crest runs out to x = 90 — the 50 ft behind the
break the prompt asked for. The foundation's line runs from x = −30 to x = 90 at
elevation 0, so the ground continues 30 ft past the toe, and `max_depth` = −10
puts the rigid rock 10 ft below the toe as an elevation rather than a thickness.
The embankment reads γ = 130 pcf and c = 400 psf, the foundation γ = 135 pcf and
c = 800 psf, both with φ = 0 and the pore pressure option `none` — the undrained
reading the prompt states, entered as a total-stress strength. The materials are
listed fill first, which is the order the profile lines reference.

Three starting circles came from `generate_starting_circles` at the shared
center (20, 40): Depth −4.72 through the toe, Depth −10 tangent to the rock, and
Depth 0 tangent to the contact — the same three LEM-3's Studio walkthrough
documents. The model checks came back clean on the first pass.

Spencer with a search then returned **FS = 1.244** on the circle Xo = 18.27,
Yo = 43.84, R = 43.84, bottoming at elevation 0, entering at x = 55.07 behind the
crest break and leaving at x = 4.46 just above the toe. It reported both notes
the run printed: interslice tension (a most tensile −1,434 lb/ft against a
largest compression of 7,313) and a line of thrust outside the slice on 18% of
the interslice boundaries. It read the surface correctly as one that stays in the
embankment and skims the top of the foundation rather than cutting into the
stronger clay below it, and it offered a Bishop or Ordinary Method of Slices
cross-check.

The session saved
[w1_build_from_image_after.xlsx](files/w1_build_from_image_after.xlsx), and the
whole exchange is in
[w1_build_from_image_transcript.md](files/w1_build_from_image_transcript.md).

<!-- test: file=files/w1_build_from_image_after.xlsx, type=circular_search, method=spencer, num_slices=40, expected_fs=1.244, tolerance=0.005 -->

### Check its work

- **Reopen the saved workbook and run Spencer with a search: 1.2441** on
  (18.27, 43.84) with R = 43.84, bottoming at elevation 0 — the 1.244 the session
  reported, from a model read off disk with no session behind it.
- **Read it against the published answer.** [LEM-3](lem03_layered_slope.md)
  gives 1.244 on (18.50, 43.75) with R = 43.75, also tangent to the contact. The
  two searches start from different circle sets — three generated here against
  the two LEM-3's file carries — and end on the same mechanism and the same
  factor of safety to four figures.
- **Open the Materials editor.** γ = 130 pcf with c = 400 psf on the embankment
  and γ = 135 pcf with c = 800 psf on the foundation, φ = 0 and pore pressure
  `none` on both, fill listed first. The material order fixes the Mat IDs the
  profile lines reference, so a foundation entered first would be a section built
  upside down.
- **Open the Profile lines editor.** Line 1 reads (−30, 0), (0, 0), (40, 20),
  (90, 20) on the embankment and line 2 reads (−30, 0), (90, 0) on the
  foundation, with **Max depth** at −10. Line 1 carries one vertex more than
  LEM-3's, which starts its first line at the toe; the extra vertex gives the
  embankment zero thickness out to x = −30, so the two sections are the same
  shape.
- **Open the Circles editor.** Three circles at (20, 40) with Depth −4.72, −10
  and 0 — one through the toe and one at the base of each layer, which is what a
  layered section needs before a search can compare the two mechanisms.
- **Both materials carry E = 0 and ν = 0.** The snippet wrote the stiffness
  fields as zeros, which is harmless for a limit equilibrium run and matters the
  moment a finite element run enters, as the strength reduction section below
  does.
- **Read the units on the admissibility note.** The reply labels the interslice
  force in psf; it is a force per foot of thickness, lb/ft, as the report's own
  `Z (lb/ft)` column and [LEM-3](lem03_layered_slope.md) both write it.

---

## Modifying the model

We will take the finished slope and change three things in one conversation: the
face from 2:1 to 3:1, a piezometric line high enough to submerge the whole
section, and the foundation's cohesion down to 250 psf. Each edit is a request
in plain language against a model already built, and each is rerun so we can
watch the factor of safety move — and, on the third, watch which circle governs
move with it.

We open LEM-3's completed model,
[xslope_simple_mult_layers.xlsx](../lem/files/xslope_simple_mult_layers.xlsx),
and send three requests in one conversation. First the face:

<div class="prompt-block" markdown>
```text
Change the face to 3:1 and rerun the search.
```
</div>

Then the water:

<div class="prompt-block" markdown>
```text
Add a piezometric line at elevation 30 across the whole section so the slope is fully submerged. Use the piezo pore-pressure option on both soils and enter saturated unit weights of 135 pcf for the embankment and 140 pcf for the foundation. Rerun the search.
```
</div>

Then the strength:

<div class="prompt-block" markdown>
```text
Reduce the foundation cohesion to 250 psf and rerun the search. Which circle governs now, and why?
```
</div>

![The three-turn conversation in the dock: the geometry read back and rewritten for the 3:1 face, then the piezometric line and the saturated unit weights, then the cohesion change, each with a before-and-after table and a searched factor of safety](images/w1_modify_3.png){width=560}

The first turn reads before it writes. Two snippets print the polygons, the
profile lines and the ground surface; the third rewrites profile line 1 as
(0, 0), (60, 20), (110, 20) — the toe held, the crest elevation held, only the
break moved out to 3:1 — and calls `resync_geometry()`. Regenerating the starting
circles takes it two more snippets: the first call's return value went unassigned
and the printout came back with the old circles in it, which it noticed and
fixed, landing the pair at the new center (30, 40). Then it runs. **FS = 1.546**
on (32.98, 53.33)
with R = 53.33, bottoming at elevation 0 again, with a line of thrust outside the
slice on 10% of boundaries.

**That turn also widened the section, which nothing asked for.** It rewrote
profile line 2 from (−30, 0)–(90, 0) to (−40, 0)–(110, 0), and said so in the
same table it reported the face in: flat ground extended on both sides so that
trial circles daylight inside the model, at least twice the 20 ft height beyond
the toe and beyond the crest. The reasoning holds, and the change is disclosed
rather than hidden — but it is a change to the model, and it means this 1.546 is
not on the extents the file was built with.

The second turn adds the piezometric line at (−40, 30) to (110, 30), sets
`u = piezo` on both materials and enters γsat = 135 and 140 pcf. **FS = 2.767**
on (32.98, 53.76) with R = 53.76, no admissibility notes at all. The assistant
added no distributed load for the pool and named the reason: `water_loads` is
already `auto`, so the hydrostatic load on the submerged profile is derived from
the piezometric line itself, and a `dloads` row would count the same water twice.
It gave the φ = 0 reason for the rise as well — pore pressure cannot reduce an
undrained strength, while the water standing on the face and crest adds a large
stabilizing moment — and cautioned that a drained case would need c′ and φ′
instead.

The third turn drops the foundation to 250 psf and returns **FS = 1.418** on
(30.27, 48.85) with R = 58.85, bottoming at **elevation −10** and exiting at
x = −2.55, out past the toe. The mechanism has moved: with the foundation at
800 psf the critical surface stayed in the embankment, and at 250 psf the
foundation is the weakest material in the section, so the search drives the
surface as deep as the model allows and it stops at −10 because `max_depth` is
the floor. That is the switch [LEM-3](lem03_layered_slope.md) makes its point
out of. The turn closed by naming two things it had not tested: no second method,
and no re-seeded starting circles after the strength change.

One workbook was saved per turn —
[w1_modify_after_1.xlsx](files/w1_modify_after_1.xlsx),
[w1_modify_after_2.xlsx](files/w1_modify_after_2.xlsx) and
[w1_modify_after_3.xlsx](files/w1_modify_after_3.xlsx) — and the conversation is
in [w1_modify_transcript.md](files/w1_modify_transcript.md).

<!-- test: file=files/w1_modify_after_1.xlsx, type=circular_search, method=spencer, num_slices=40, expected_fs=1.546, tolerance=0.005 -->
<!-- test: file=files/w1_modify_after_2.xlsx, type=circular_search, method=spencer, num_slices=40, expected_fs=2.767, tolerance=0.005 -->
<!-- test: file=files/w1_modify_after_3.xlsx, type=circular_search, method=spencer, num_slices=40, expected_fs=1.418, tolerance=0.005 -->

### Check its work

- **Reload each saved workbook and run the search yourself: 1.5460, 2.7673 and
  1.4184** — the three numbers reported, to four figures, on the three circles
  reported.
- **Read the lowest point beside each factor of safety**, not only the factor of
  safety: elevation 0, then 0, then −10. The third run is a different mechanism,
  not the same surface re-solved, and the answer says which circle governs and
  why.
- **Open the first workbook's Profile lines editor.** Line 2 runs from x = −40 to
  x = 110, wider on both sides than the model we opened. The widening was
  reported but not requested, so decide whether to keep it before comparing this
  1.546 against anything computed on the original extents.
- **Open the second workbook's Materials editor.** Both materials read
  `u = piezo`, with γsat = 135 and 140 pcf beside their total unit weights of 130
  and 135 — the water is explicit and the unit weights are total, which is how
  XSLOPE wants a submerged section stated. A saturated weight that lands under
  some other name is caught rather than ignored: the model checks flag a material
  field no analysis reads and name the nearest real one — *"Did you mean
  'gamma_sat'?"*
- **Check the method.** The file names no method, so each run took XSLOPE's
  default of Spencer; the first turn's heading calls Spencer "the method the
  model declares", which is the one loose phrase in the three turns.

---

## A sweep with the helper

We will ask for the foundation's cohesion swept from 200 to 800 psf in steps of
100, with a search at every step. `sensitivity()` — a sweep helper preloaded in
the assistant's kernel, the counterpart of Studio's **Parametric…** dialog — runs
a value list through a callback and searches at each step, so the request tests
whether the assistant reaches for the machinery already there rather than
rebuilding it.

<div class="prompt-block" markdown>
```text
Use the sensitivity helper to sweep the foundation cohesion from 200 to 800 psf in steps of 100 psf, Spencer with a search at each step, and give me the factor of safety at each. Leave the model as it was.
```
</div>

![The sweep in the dock: one code block calling the sensitivity helper with an apply closure, the seven searched steps printed as they ran, and the finished table with its reading](images/w1_sweep_builtin_1.png){width=560}

It called the helper with an `apply` closure that writes the foundation's `c`,
`method='spencer'` and `search=True`, wrote `sensitivity.csv` and
`sensitivity.png`, and printed the restored cohesion — `restored c = 800.0` — at
the end. **Files…** opens the folder holding the two files.

| Foundation c (psf) | FS (Spencer, searched) |
| :---: | :---: |
| 200 | 0.640 |
| 300 | 0.792 |
| 400 | 0.964 |
| 500 | 1.135 |
| 600 | 1.244 |
| 700 | 1.244 |
| 800 | 1.244 |

Its reading of the curve holds up. The rise from 200 to 500 psf is close to
linear at about +0.165 per 100 psf; the last three rows are identical to six
decimals, 1.244108 each; the factor of safety passes 1.0 at about 420 psf by
interpolation between the 400 and 500 psf rows; and the model as it stands, at
800 psf, sits on that plateau. It also said what it had not done: it measured the
flat factor of safety but did not extract the critical surfaces, so it could not
say from this run where the governing surface sits.

Finishing that last point by hand takes one more solve per step, and it is what
turns the flat tail of the table into a statement about the slope.

### Check its work

- **Re-run any of the seven steps**: 0.639688, 0.792223, 0.963558, 1.135100 and
  1.244108 three times over — the reported values to six decimals.
- **Read two rows against the published answers.** The 300 psf row is 0.792,
  which is [LEM-3](lem03_layered_slope.md)'s weak-foundation answer, and the
  800 psf row is its 1.244; the sweep reproduces both ends of that tutorial from
  one request.
- **Pull the critical circle at each step, which the answer says it did not.**
  Below 600 psf every critical surface bottoms at elevation −10, tangent to the
  rock; at 600 psf and above it bottoms at elevation 0, tangent to the contact.
  Across the plateau the mechanism has left the foundation, which is why adding
  foundation strength buys nothing there.
- **Open the Materials editor when it finishes.** The foundation reads 800 psf,
  the helper restored the project, and no workbook was written.
- **The prose names `materials[2]['c']` where the code changed
  `materials[1]['c']`** — the same material, named once by its 1-based Mat ID and
  once by its 0-based position in the list.

---

## A sweep written ad hoc

We will ask for the factor of safety at face slopes of 2:1, 2.5:1 and 3:1, with
the toe and the crest elevation held and a search each time. What varies here is
the shape of the section rather than a value in a cell, so no dialog offers the
study: whatever runs it has to rebuild the geometry between steps.

<div class="prompt-block" markdown>
```text
Compute the factor of safety for face slopes of 2:1, 2.5:1 and 3:1, keeping the crest elevation and the toe where they are, Spencer with a search each time. Give me a table. Leave the model as it was.
```
</div>

![The face-slope sweep in the dock: the geometry read back, the closure that moves the crest break and regenerates the circles, the three searched steps, the finished table, and the restored geometry printed at the end](images/w1_sweep_adhoc_1.png){width=560}

It reads the geometry first, then hands the same `sensitivity()` helper an
`apply` that does the whole geometry edit: it moves the crest break to x = 2H,
2.5H and 3H for a 20 ft height, holds the toe at (0, 0) and the crest at
elevation 20, calls `resync_geometry()`, and regenerates the starting circles so
each face gets a set built from its own section. It also gave all three cases the
same extents, ground from x = −50 to x = 110, so that every trial circle
daylights on real ground under conditions the three runs share.

| Face slope | Crest break | FS (Spencer, searched) |
| :---: | :---: | :---: |
| 2:1 | (40, 20) | 1.244 |
| 2.5:1 | (50, 20) | 1.396 |
| 3:1 | (60, 20) | 1.546 |

The restoration is printed rather than promised: the closing snippet prints the
profile lines, the circles and the polygons back at their original values, and no
workbook was written.

Its closing account of the study is right about everything it ran, and wrong
about the one thing it did not:

> - Each case used the same wide extents (ground from x = −50 to x = 110, ~2.5 H
>   beyond toe and crest) so every trial circle daylights on real ground, and the
>   hard base stayed at the model's own `max_depth = −10`.
> - Starting circles were regenerated for each geometry with
>   `generate_starting_circles` (toe circle plus one tangent to the foundation
>   base at −10), then Spencer ran its own automated search each time.
> - FS rises about 0.15 per half-unit of flattening over this range — essentially
>   linear here, as expected for a φ = 0 two-layer section where the critical
>   circle stays deep-seated in the foundation clay.

The first two bullets are the run exactly, and the rise really is about 0.15 per
half-unit of flattening. Not one of the three critical circles is deep-seated in
the foundation: all three bottom out at elevation 0, tangent to the contact and
contained in the embankment, which is the fact LEM-3 makes its own result out
of. The sweep reported factors of safety and no surfaces, so what settles the
sentence is one re-run — the second check below extracts the three circles, and
every one comes back at elevation 0.

### Check its work

- **Re-run the three steps**: 1.244146, 1.396287 and 1.545993, the reported
  values to six decimals.
- **Read the critical circle at each step, not only the factor of safety.**
  (18.27, 43.84) with R = 43.84, (25.86, 47.73) with R = 47.73 and (32.98, 53.33)
  with R = 53.33 — every one bottoming at elevation 0, which is what settles the
  closing sentence.
- **Cross the 3:1 row against the previous section.** The modify conversation
  searched its own 3:1 geometry and returned 1.5460 on the same circle, so two
  conversations agree on that face.
- **Read the 2:1 row against the published answer.** It returns 1.244146 against
  LEM-3's 1.244108 — the same mechanism, and the difference is the starting
  circles the sweep regenerates at each step, not the wider ground: the file's
  own two circles return 1.244108 on either set of extents.
- **Open the Profile lines and Circles editors when it finishes.** The section is
  back at (0, 0)–(40, 20)–(90, 20) over (−30, 0)–(90, 0) with the original two
  circles, and no workbook was written.

---

## Stiffnesses and a strength reduction run

We will ask what elastic properties these two soils should carry, have the
assistant enter them, and then run the strength reduction analysis and set it
beside Spencer. A limit equilibrium analysis needs no stiffnesses and the finite
element engine will not start without them, so this is the gap every model
crossing from one engine to the other has to close.

<div class="prompt-block" markdown>
```text
Suggest elastic properties for both soils and enter them.
```
</div>

Then, once they are in:

<div class="prompt-block" markdown>
```text
Run a strength reduction analysis and compare the result with Spencer.
```
</div>

![Both turns in the dock: the classifier's two suggestions written onto the materials, then the Spencer search, the mesh, the seven bisection steps with their two inconclusive trials, and the comparison table](images/w1_elastic_fem_2.png){width=560}

The first turn calls `suggest_elastic()`, which classifies each material from its
strength: the embankment at c = 400 psf and φ = 0 as **Soft Clay**, at
E = 167,100 psf with ν = 0.45, and the foundation at c = 800 psf as **Medium
Clay**, at E = 668,300 psf with ν = 0.40. It wrote both onto the materials and
quoted the tool's own reason back rather than dressing the numbers up as
measurements — E the midpoint of the published range for that soil type, ν its
typical value, a last-resort fill that any stated value outranks. It added that
the two fields reach the finite element run only, since limit equilibrium ignores
them, and it warned that a strength reduction run takes several minutes before
offering to start one.

The second turn runs both engines. Spencer with a search comes back at
**FS = 1.244** on (18.50, 43.75) with R = 43.75, and `run_fem(analysis='ssrm')`
builds its own mesh first — tri6, six-node triangles, at a target size of 1.2 ft,
giving 8,849 nodes and 4,304 elements — then bisects from [1.00, 2.00] in seven
steps to a shear strength reduction factor of safety of **1.254**, 0.010 above
Spencer, about 0.8%. The run took 528 s.

It reported three things about those numbers without being asked. **The bracket
did not close:** the trials at F = 1.2656 and F = 1.2578 hit the
50,000-iteration ceiling with the out-of-balance force still falling, which
counts as inconclusive rather than as failure, so 1.254 is the midpoint of a
final interval of [1.2500, 1.2578] whose upper edge is undecided, and raising
`max_iterations_ceiling` is what would settle it. **Spencer's two admissibility
notes carry into the comparison** — interslice tension and a line of thrust
outside the slice on 15% of boundaries — which it read as a reason to want the
finite element cross-check rather than as a reason to distrust it. And **the
stiffnesses are estimates**: it named the sensitivity of the strength reduction
answer to E and ν as untested and priced the test at several minutes per point.

The session saved
[w1_elastic_fem_after.xlsx](files/w1_elastic_fem_after.xlsx), and the exchange is
in [w1_elastic_fem_transcript.md](files/w1_elastic_fem_transcript.md).

### Check its work

- **Open the Materials editor in the saved workbook.** The embankment carries
  E = 167,100 psf with ν = 0.45 and the foundation E = 668,300 psf with ν = 0.40 —
  the classifier's own values, entered unedited, with the strengths untouched.
- **Read the strength reduction answer as a bracket.** The final interval is
  [1.2500, 1.2578] with two undecided trials on its upper edge, and 1.254 is the
  midpoint; the value could sit anywhere in that interval.
- **Check the number it is compared against.** Spencer's 1.244 on
  (18.50, 43.75) is [LEM-3](lem03_layered_slope.md)'s published answer, so the
  0.8% is a gap between two engines on one model rather than between two models.
- **The mesh came from defaults, not from the file.** The summary calls it "built
  automatically from the model's own settings", but this model declares neither
  an element type nor a mesh size; tri6 and the 1.2 ft target are the mesh
  builder's own defaults, the second being the 120 ft ground-surface width over
  100. Rebuilding at another size is how to find out what the answer owes to it.
- **Ask for the stiffness test rather than assuming its result.** A second
  strength reduction solve at a different E takes minutes and settles how far the
  classifier's estimates move the answer.

---

## Two documentation questions

We will ask two questions about the analysis rather than for a change to the
model: why every limit equilibrium method returns the same factor of safety on
this slope, and what pore pressure option two undrained soils should carry. The
first question's premise is false, which makes it the more useful of the two.

<div class="prompt-block" markdown>
```text
Why do all of the limit equilibrium methods give the same factor of safety on this model?
```
</div>

Then, in the same conversation:

<div class="prompt-block" markdown>
```text
Both soils are undrained. What should the pore-pressure option be, and what would change in the analysis if I added a water table?
```
</div>

![Both answers in the dock: the seven methods solved on the circle the file defines and tabulated by the equilibrium each satisfies, the moment argument for the tie at phi = 0, then the pore-pressure answer with its two measured Bishop runs, the restored model, and the three effects a water table would have](images/w1_conceptual_2.png){width=560}

**It refused the premise.** One snippet ran all seven methods on the first of the
two circles the file carries, and the answer opens "Not all of them do — the ones
that satisfy **moment** equilibrium agree exactly, and the force-equilibrium
methods don't". It named the basis in the same breath: solved on the surface the
model already defines, Xo = 20, Yo = 40, R = 40, 40 slices, no search.

| Method | Equilibrium satisfied | FS |
| :--- | :--- | :---: |
| Ordinary Method of Slices (OMS) | moment only | 1.247 |
| Bishop | moment, and vertical force | 1.247 |
| Spencer | force and moment | 1.247 |
| Morgenstern-Price | force and moment | 1.247 |
| Lowe & Karafiath | force only | 1.299 |
| Janbu (corrected, f₀ = 1.09) | force only | 1.321 |
| Corps of Engineers | force only | 1.363 |

The argument under it is the right one. With φ = 0 on a circular arc, the
resisting moment about the center is the cohesion times the arc length times the
radius, so "the base normal force N drops out entirely"; the only thing an
interslice assumption changes is how N is distributed, and every method that
takes moments about the center therefore returns the same number whatever it
assumes. It added the two constants that make the point concrete — Spencer
settles at θ = 7.91° and Morgenstern-Price at λ = 0.151, satisfying force
equilibrium without moving the moment answer. The force-only methods never write
the moment equation, so their assumptions reach N and they spread from 1.299 to
1.363. Three documentation pages are cited for the reading behind that:
`lem/overview/`, `lem/bishop/` and `lem/spencer/`.

The second answer keeps `u = none`, which is what the model already carries and
what a total-stress undrained analysis takes: the undrained strength already
contains the pore water, so subtracting a pore pressure would count it twice.
Then it measured that rather than asserting it — a piezometric line 3 ft below
the ground surface everywhere, so nothing ponds, and Bishop on the first of the
file's two circles with the option off and then on:

| Configuration | FS (Bishop, the file's first circle) |
| :--- | :---: |
| `u = none`, piezometric line present | 1.2472 |
| `u = piezo`, same line | 1.2472 |

"Identical to four decimals," it wrote, with the reason in the strength equation:
the base strength is c + (N/L − u)·tan φ, and at φ = 0 the term the pore pressure
sits in is multiplied by zero. It printed the model restored — both materials
back at `u = none`, the piezometric line empty.

It then separated the three things a water table would do, which is the part of
the answer that goes past the question. The base pore pressure does nothing at
φ = 0, as measured. The slice weight changes only if a saturated unit weight is
entered, and it read both materials to check — `gamma_sat` is empty on this model
— with the caution to give the total saturated unit weight and never a buoyant
one, since a buoyant weight beside an explicit water definition counts buoyancy
twice. And water standing above the ground is a real effect: a hydrostatic load
on the toe-side ground, derived from the water definition itself because
`water_loads` ships on `auto`, and never also entered as a `dloads` row.

### Check its work

- **Check that nothing moved.** Neither turn edited the model, and no workbook
  was written; the pore pressure option is back at `none` on both materials and
  the piezometric line is empty again.
- **Reproduce the two rows.** Both configurations return 1.247168 on the file's
  first circle, against 1.247119 with no piezometric line at all — the line's
  presence moves the fifth decimal, and switching the option on moves nothing.
- **Measure the effect it named but did not run.** Enter γsat = γ + 5 on both
  materials, leave `u = piezo` and re-solve the same circle: **1.2127** against
  1.2472. The weight split accounts for 0.034 here; only the pore pressure is inert
  at φ = 0.
- **Read the seven methods against the published table.**
  [LEM-3](lem03_layered_slope.md) searches each method its own critical circle
  and reports OMS, Bishop, Spencer and Morgenstern-Price at 1.244 with Lowe at
  1.285, Janbu at 1.313 and Corps at 1.326. The four-way tie and the force
  methods sitting above it are the same finding; the numbers differ because these
  seven are one fixed circle.
- **Solve that circle with Analysis = Single surface**: 1.247, not the 1.244 a
  search returns. The answer states which of the two it computed, which is what
  makes the table readable at all.
- **Follow the three documentation pages.** `lem/overview/`, `lem/bishop/` and
  `lem/spencer/` are real pages, and the derivation there is what the answer
  summarizes. The transcript is
  [w1_conceptual_transcript.md](files/w1_conceptual_transcript.md).

---

## A broken file

We will open a copy of the slope with three transcription errors written into it
and say only that the factor of safety looks wrong. No request in this tutorial
carries less information, and none depends more on testing the answer that comes
back.

That copy, [w1_diagnose_start.xlsx](files/w1_diagnose_start.xlsx), carries three
faults typed into LEM-3's model. The material rows are **swapped** — foundation
first, embankment second — while the profile lines still point at material 1 for
the ground surface and material 2 for the layer beneath, so the strong clay sits
in the fill and the weak clay in the foundation. The embankment's unit weight
reads **13 pcf** instead of 130, a dropped digit. And the maximum depth reads
**−100** instead of −10, another dropped digit, putting a rigid base 100 ft below
a 20 ft section. Spencer with a search returns 1.004 on it, on a circle
bottoming 24.6 ft below the toe.

<div class="prompt-block" markdown>
```text
This model was built from the LEM-3 tutorial, but the factor of safety looks wrong. Find what is wrong, fix it, and rerun.
```
</div>

![The diagnosis in the dock: the model read once, the baseline search, the snippet that varies one input at a time and re-searches each, the repair, the rerun at 1.244, and the closing section on the input it did not change](images/w1_diagnose_1.png){width=560}

It answers by changing things rather than by reading them. After two reads of
the model and a baseline search, one snippet fixes the unit weight alone,
un-swaps the material assignment alone, then does both — searching after each
change and putting the model back at the end:

| Model state | FS (Spencer, searched) |
| :--- | :---: |
| As received | 1.004 |
| Unit weight 13 → 130 pcf only | 0.994 |
| Material assignment un-swapped only | 12.441 |
| Both fixed | 1.244 |

Both of those faults are named precisely. Profile line 1 is the embankment
surface and carried the foundation's Mat ID while the base line carried the
embankment's: "The strong clay was in the fill and the weak clay in the
foundation, exactly backwards." The 13 pcf is read as a dropped digit, and the
absurd 12.441 of the swap-only run is explained by it — "at 13 pcf the fill was
essentially weightless" — so the two faults had been pulling the answer in
opposite directions, which is why the file arrived at a plausible-looking
1.004. It then applied both fixes and reran: **FS = 1.244** on
(18.50, 43.75) with R = 43.75, bottoming at elevation 0.

The third fault gets a section of its own, headed *"One thing I did not
change"*: `max_depth = −100` puts a rigid base 100 ft below a 20 ft section,
nothing in the model describes a foundation that thick, it does not affect this
answer because the critical circle bottoms at elevation 0 — and if the tutorial
states a firm-layer elevation, say so and it will be set, because it will not
guess a base depth. Both halves of that are true, and naming a fault while
refusing to invent a value for it is the right way to leave one open.

**It repaired the swap the other way round.** Rather than reordering the material
rows, it re-pointed the profile lines' Mat IDs at the rows as they stand. The
section, the strengths along it and the answer are LEM-3's; the material table's
row order is not, and the saved
[w1_diagnose_after.xlsx](files/w1_diagnose_after.xlsx) still lists the foundation
first. A defensible fix, and one to notice before the next person reads that
table top down.

<!-- test: file=files/w1_diagnose_start.xlsx, type=circular_search, method=spencer, num_slices=40, expected_fs=1.004, tolerance=0.005 -->
<!-- test: file=files/w1_diagnose_after.xlsx, type=circular_search, method=spencer, num_slices=40, expected_fs=1.244, tolerance=0.005 -->

### Check its work

- **Reload the repaired workbook and search: 1.2441** on (18.50, 43.75) with
  R = 43.75 — the circle [LEM-3](lem03_layered_slope.md) publishes, from the file
  the session left behind.
- **Open its Materials editor before trusting the row order.** The foundation is
  still listed first and the embankment second, with profile line 1 pointing at
  material 2 and line 2 at material 1.
- **Test the input it left alone.** The repaired model still carries
  `max_depth = −100` and returns 1.2441076757; LEM-3's file carries −10 and
  returns 1.2441076757. The base depth is wrong and does not touch this answer,
  which is what the closing section says.
- **Read the baseline's search notes, which the answer quotes only in part.**
  Spencer could not solve 79 of the 533 trial surfaces on the broken model and 5
  of those rank below the reported minimum, so the 1.004 the file came in with is
  not a number to carry anywhere.
- **Read the admissibility line carefully.** It reports a line of thrust outside
  the slice on "26% → 15% of boundaries", which reads as one quantity moving; the
  26% belongs to the broken model and the 15% to the repaired one.

---

## Generating the report

We will run Spencer with a search and then ask for the analysis report in one
sentence. What comes back is the same document **File → Generate Report…**
writes, built from the same run, so this request shortcuts the dialog rather
than making a second kind of report.

The report documents what the session has solved, so the run comes first:

<div class="prompt-block" markdown>
```text
Run Spencer with a search.
```
</div>

Then the document:

<div class="prompt-block" markdown>
```text
Write the analysis report.
```
</div>

![Both turns in the dock: the search returning 1.244 with its circle laid out in a table, then the one call to generate_report, the path of the finished document, and the document offered as an attachment](images/w1_report_2.png){width=560}

The first turn returns **FS = 1.244** on the circle (18.50, 43.75) with R = 43.75
— [LEM-3](lem03_layered_slope.md)'s published answer, circle included. It read
the surface straight off the result, laid `Xo`, `Yo`, `R`, `Depth` and the entry
and exit stations out in a small table, and read the lowest point at elevation 0
as sitting well above the rigid base at −10, so the floor is not controlling the
answer. It reported both admissibility notes.

The second turn makes one call to `generate_report`. The document carries six
numbered figures and three numbered tables across *1 Traceability*,
*2 Project Definition* and
*3 Limit Equilibrium Analysis*, the last holding 3.1 Analysis Inputs,
3.2 Materials, 3.3 Loads and 3.4 Spencer's Method, itself split into the search,
the results, the slice table and the calculations. It lands in the folder
**Files…** opens, and the dock shows it as an attachment with a *show in folder*
link beside it. The session saved
[w1_report_after.docx](files/w1_report_after.docx).

### Check its work

- **Open the report at §3.4.2.** It reads "Spencer's Method gives a factor of
  safety of 1.244 on the critical surface of Figure 4" — the run made in the
  first turn, not a second analysis the report performed for itself.
- **Read the traceability page.** It names the file the session ran on and its
  digest:

```text
xslope_simple_mult_layers.xlsx
d0307449074ab1094764c075b45ec241a098fbe10731a5cf497b1846581da1bb
```

  Hashing `docs/lem/files/xslope_simple_mult_layers.xlsx` returns the same 64
  characters, so the document names the model it was written from.
- **Read what it says about water and loads.** The report states that the model
  defines no groundwater and no external water and that it carries no distributed
  loads, which is this section as built.
- **Read §3.4.4 against the run.** The calculations section reports the iteration
  converging at F = 1.244 with θ = 8.48°, the pair Spencer solves for.
- **One line of the summary is about the document, not the analysis.** The
  assistant reports the document as arriving "with the contents page finished
  and page-numbered", where the third paragraph of the document reads "Page
  numbers appear when the table is updated in Word (right-click → Update Field,
  or select all and press F9)". (Its "two-zone polygon geometry" is right: the
  two profile lines are converted to two material polygons when the model is
  built.)

---

## What it will not do

**It does not convert units.** A model is in one unit system throughout, and
nothing in XSLOPE converts one into another: the numbers entered are the numbers
the engine reads, under whatever system the project declares. Convert a drawing
in meters and kilonewtons yourself before handing it to a project in feet, and
check the entered values against the drawing afterwards.

**Asked for a value the model does not state, it named the gap rather than
filling it.** The diagnosis above found a maximum depth of −100 under a 20 ft
section, said what was wrong with it, measured that the repaired answer does not
depend on it, and then left it alone and asked for the firm-layer elevation. An
input that comes back unchanged with a question attached is usually this, not an
oversight. Treat it as behavior to check rather than as a guarantee: asked the
same thing, another of the models above invented a unit weight instead of
reading the one the tutorial gives.

**It can explain a result it never measured.** The face-slope sweep returned
three correct factors of safety and closed by attributing them to a circle
"deep-seated in the foundation clay" when all three surfaces are tangent to the
contact and never enter the foundation. Nothing in the reply reads as a guess,
the numbers around the sentence are right, and no model check has anything to say
about it. Reading the critical surface beside each factor of safety is what finds
it.

---

## A harder problem

As the slope and the prompts get more complicated, the opportunity for error
rises. The same eight tasks were also run on the reinforced slope of [LEM-8](lem08_reinforced_slope.md) — a
2 ft cohesive band along the face, six geogrid layers with pullout at both ends,
and a surcharge across the crest — scored the same two ways. Claude Opus 5 is
again right on every model and every number, at $2.23 against the $1.71 the
layered slope cost it, and again loses three of the second column. The five models
above score 5 to 8 of 8 on model and numbers there and 3 to 6 of 8 on
explanation, and a sixth and smaller one, Claude Haiku 4.5, gets 2 and 3. Late
August 2026, on the same list prices as the table above:

| Model | Cost | Model and numbers | Explanation |
| :--- | :---: | :---: | :---: |
| Claude Opus 5 | $2.23 | 8 of 8 | 5 of 8 |
| OpenAI gpt-5.5 | $1.12 | 5 of 8 | 5 of 8 |
| Kimi K3 | $1.20 | 5 of 8 | 3 of 8 |
| Claude Sonnet 5 | $1.31 | 5 of 8 | 6 of 8 |
| Claude Haiku 4.5 | $0.53 | 2 of 8 | 3 of 8 |
| GLM-5V-Turbo | ~$0.24 | 5 of 8 | 4 of 8 |

Three things did most of the damage. The 2 ft cohesive band along the face is a
thin zone read off a drawing — the reference build measures it perpendicular to
the face — and a model that missed it entirely gave the whole embankment one
strength and carried that into every number afterwards. Moving the crest break
for a flatter face left the crest surcharge behind on the old geometry in more
than one session. And in the diagnosis, several models measured the effect of a
planted fault and then filed it under something other than a fault — which is
why the checks in this tutorial re-run the measurements rather than read the
conclusions.

---

## Conclusion

This tutorial covered:

- Setting the assistant up: a provider, a model and a key in **Settings…**, and
  13 turns of work for about $1.71.
- Building a layered model from a drawing and editing it three ways — face,
  water and strength — with every reported factor of safety reproducing from the
  workbook the session saved.
- Two sweeps, a strength reduction run, two questions answered without touching
  the model, a diagnosis that found all three planted faults, and a report
  carrying the run it names.
- Checking what came back: reload the workbook and re-run the search, read the
  critical surface and not only the factor of safety, and re-measure any
  mechanism the answer explains — one sweep's closing sentence contradicts the
  surfaces its own run found.

**Where to go next:** The [AI Assistant reference](../studio/assistant.md)
documents the helpers the assistant calls, the checks that run after every edit
it makes, and what each provider can do. [LEM-3](lem03_layered_slope.md) builds
this same slope by hand three ways, and [LEM-8](lem08_reinforced_slope.md) builds
the reinforced slope the comparison above runs on.
