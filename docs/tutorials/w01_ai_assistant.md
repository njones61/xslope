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
release to release. Where a provider's catalog mixes models that read images
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
spends 177,000 across nine calls.

Models differ in both price and capability, and the two do not track each other
perfectly: a cheaper model may answer a question about the documentation as well
as an expensive one and still fail at a long autonomous build. The line under
the input box reports the tokens each turn spent and the running total for the
conversation:

```
this turn: 52,505 in (31,702 cached) / 2,715 out · session: 52,505 in (31,702 cached) / 2,715 out
```

One request is usually several calls to the model, so that first number climbs
while the assistant works. **New chat** starts the count over.

The eight conversations below cost this much between them, measured as they ran:

| Session | Turns | Model calls | Tokens in (cached) | Tokens out | Wall |
| --- | :---: | :---: | :---: | :---: | :---: |
| Building a model from a drawing | 1 | 3 | 52,505 (31,702) | 2,715 | 47 s |
| Modifying the model | 3 | 12 | 230,865 (190,212) | 6,198 | 122 s |
| A sweep with the helper | 1 | 2 | 33,236 (31,702) | 860 | 56 s |
| A sweep written ad hoc | 1 | 5 | 87,126 (79,255) | 2,759 | 54 s |
| Stiffnesses and a strength reduction run | 2 | 4 | 69,743 (47,553) | 1,485 | 556 s |
| Two documentation questions | 2 | 5 | 92,550 (79,255) | 5,250 | 88 s |
| A broken file | 1 | 9 | 176,923 (142,659) | 8,544 | 173 s |
| Generating the report | 2 | 4 | 67,664 (63,404) | 839 | 30 s |
| **Total** | **13** | **44** | **810,612 (665,742)** | **28,650** | **1,125 s** |

Thirteen turns, 44 calls to the model, about 811,000 input tokens of which some
666,000 came from the prompt cache, and 28,650 output tokens: **\$1.77** at
Anthropic's list prices on 2026-08-29 (\$5.00 per million input tokens, a cache
read at a tenth of that, \$25.00 per million output tokens). Wall time and price
do not track each other: the strength reduction run takes half the total time for
under a tenth of the cost, because 521 of its 556 seconds are spent inside the
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
sentence that is not, and the diagnosis below is that case in full.

| Model | Calls | Tokens in (cached) | Tokens out | Cost | Model and numbers | Explanation |
| --- | :---: | :---: | :---: | :---: | :---: | :---: |
| Claude Opus 5 | 44 | 810,612 (665,742) | 28,650 | \$1.77 | 8 of 8 | 7 of 8 |
| Claude Sonnet 5 | 45 | 935,563 (684,517) | 28,060 | \$0.92 | 7 of 8 | 7 of 8 |
| OpenAI gpt-5.5 | 29 | 364,339 (325,120) | 17,437 | \$0.88 | 7 of 8 | 6 of 8 |
| Kimi K3 (Moonshot AI) | 37 | 449,769 (343,143) | 44,953 | \$1.10 | 7 of 8 | 4 of 8 |
| GLM-5V-Turbo (Z.ai) | 40 | 563,429 (459,660) | 13,480 | ~\$0.29 | 6 of 8 | 4 of 8 |

All five reproduce LEM-3's published 1.244 on the contact circle — from a build
off the drawing, from the file, and inside the report document. The differences
are everywhere else. Opus 5 is the only one right on every model and every
number, and its one miss is in the second column: the diagnosis below, where the
repair is right and the sentence placing the repaired surface is not. Sonnet 5
ties it on explanation at half the price and loses only the diagnosis — it named
all three faults, then applied the material fix to the polygons Studio derives
from the profile lines, which Studio rebuilds and says it has rebuilt, so the
fix never reached the model and 1.208 went out as the repaired answer. gpt-5.5
ran the leanest of the five, 29 calls and 364,000 input tokens for the same
eight conversations, and it is exact on every number it computes; it misses the
material swap in the broken file outright and reports that same 1.208. Kimi K3
built and edited the model exactly as Opus did and ties GLM-5V-Turbo at the
bottom of the explanation column, three of its four misses being a sentence its
own output contradicts: flattening the face called a drop in a factor of safety
its table shows rising, a rise attributed to lighter soil by a snippet that made
the soil heavier, and the top of the foundation named as its
base. GLM-5V-Turbo, at a third of the cost of the next cheapest, builds the
section correctly and then inverts the slope ratio in the ad-hoc sweep — 20 ft
of rise over a 10 ft run for 2:1 — so it answers a different problem, about
0.965 at every face against 1.244 to 1.546. Z.ai publishes no price for that
model; its cost here is at the rate a reseller lists for it.

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
back. Where the five models went wrong most often: a material swap in a broken
file missed outright, or undone by editing the polygons Studio derives rather
than the profile lines they come from; a slope ratio inverted inside the code
written for a study; and a sentence contradicted by the table printed above it.
The model checks stay clean through all three, because none of them makes a
model inconsistent, only wrong.

---

## Building a model from a drawing

The assistant can build a model straight from a drawing, which is the request
that covers the most ground at once. We will paste the problem drawing from
[LEM-3](lem03_layered_slope.md) into an empty project and ask for the model it
describes, naming the two zones and saying where the section should end — the
two things the drawing does not carry. One request covers reading dimensions off
a picture, entering two materials and two profile lines, setting the rigid base,
generating starting circles, and searching for the critical surface.

We start from **File → New**. Right-clicking the drawing at the top of that page
and choosing **Copy image** puts it on the clipboard; from there we click in the
chat box, press Ctrl/Cmd+V to attach it, and type:

<div class="prompt-block" markdown>
```text
Build this model: a 20 ft embankment on a 10 ft foundation layer over rigid rock. Both soils are undrained. Put the toe at x = 0, run the ground 30 ft past the toe and 50 ft behind the crest break. Add starting circles and run Spencer with a search.
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
reading the prompt states, entered as a total-stress strength. Nothing in the
request names a unit system; the snippet declares the project imperial, which is
the system the drawing's pcf and psf are in. The materials are listed fill
first, which is the order the profile lines reference.

Three starting circles came from `generate_starting_circles` at the shared
center (20, 40): Depth −4.72 through the toe, Depth −10 tangent to the rock, and
Depth 0 tangent to the contact — the same three LEM-3's Studio walkthrough
documents. The model checks came back clean on the first pass.

Spencer with a search then returned **FS = 1.244** on the circle Xo = 18.28,
Yo = 43.84, R = 43.84, bottoming at elevation 0, entering at x = 55.07 behind the
crest break and leaving at x = 4.46 just above the toe. It reported both notes
the run printed: interslice tension (a most tensile −1,434 lb/ft against a
largest compression of 7,313) and a line of thrust outside the slice on 18% of
the interslice boundaries. It read the surface off the Depth its own run
printed — *"Depth = 0.0, i.e. tangent to the top of the foundation"* — and said
the surface stays in the embankment rather than cutting into the stronger clay
below it, then offered a tension crack or a Bishop cross-check for the two notes.

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
- **Open the Profile lines editor.** Line 1 reads (0, 0), (40, 20), (90, 20) on
  the embankment — LEM-3's line exactly — and line 2 reads (−30, 0), (0, 0),
  (90, 0) on the foundation, with **Max depth** at −10. Line 2 carries one vertex
  more than LEM-3's: the point at the toe sits on the straight run from (−30, 0)
  to (90, 0) and changes nothing about the shape.
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

![The three-turn conversation in the dock: the geometry read back and rewritten for the 3:1 face, then the piezometric line and the saturated unit weights, then the cohesion change, each with the edit reported and a searched factor of safety](images/w1_modify_3.png){width=560}

The first turn reads before it writes, and its first write goes to the wrong
place. One snippet prints the polygons, the maximum depth and the circles; the
next rewrites the polygons — and Studio derives those from the profile lines on
this model, so it rebuilt them and printed *"The polygon edit did not take."*
The session read that, printed the profile lines, and redid the edit there: line
1 becomes (0, 0), (60, 20), (110, 20), the toe held, the crest elevation held,
only the break moved out to 3:1. It staged three circles for the new face at the
center (30, 40) — one tangent to the contact, one tangent to the rock, one
through the toe at (25, 40) with R = 47.17 — and ran. **FS = 1.546** on
(32.98, 53.33) with R = 53.33, bottoming at elevation 0 again, with a line of
thrust outside the slice on 10% of boundaries.

**That turn also widened the section, which nothing asked for.** It rewrote
profile line 2 from (−30, 0)–(90, 0) to (−50, 0)–(110, 0), and said so in the
same table it reported the face in: flat ground extended 50 ft beyond the toe and
50 ft beyond the crest break, so that trial circles daylight inside the model.
Its own note calls that twice the 20 ft height; against a break now at x = 60 it
is two and a half times. The change is disclosed rather than hidden — but it is a
change to the model, and it means this 1.546 is not on the extents the file was
built with.

The second turn adds the piezometric line at (−50, 30) to (110, 30), sets
`u = piezo` on both materials and enters γsat = 135 and 140 pcf. **FS = 2.767**
on (32.98, 53.76) with R = 53.76, no admissibility notes at all. The assistant
added no distributed load for the pool and named the reason: `water_loads` is
already `auto`, so the hydrostatic load on the submerged profile is derived from
the piezometric line itself. That is why a `dloads` row here would count the same
water twice.
It gave the φ = 0 reason for the rise as well — pore pressures cannot reduce an
undrained strength, so the change comes from the 10 ft of water standing on the
slope, whose weight acts as a stabilizing surcharge — and offered to isolate the
two effects by re-solving with `water_loads` set to manual.

The third turn drops the foundation to 250 psf and returns **FS = 1.418** on
(30.27, 48.85) with R = 58.85, bottoming at **elevation −10** and exiting at
x = −2.55, out past the toe. The mechanism has moved: with the foundation at
800 psf the critical surface stayed in the embankment, and at 250 psf the
foundation is the weakest material in the section, so the search drives the
surface as deep as the model allows and it stops at −10 because `max_depth` is
the floor. That is the switch [LEM-3](lem03_layered_slope.md) makes its point
out of. The turn closed by naming the consequence it had not tested: the circle
is pinned against `max_depth`, so if the clay carries on below −10, the true
factor of safety is lower than 1.418.

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
- **Open the first workbook's Profile lines editor.** Line 2 runs from x = −50 to
  x = 110, wider on both sides than the model we opened. The widening was
  reported but not requested, so decide whether to keep it before comparing this
  1.546 against anything computed on the original extents.
- **Open the second workbook's Materials editor.** Both materials read
  `u = piezo`, with γsat = 135 and 140 pcf beside their total unit weights of 130
  and 135 — the water is explicit and the unit weights are total, which is how
  XSLOPE wants a submerged section stated.
- **Check which method ran.** The file names no method, so every
  `run_lem(search=True)` in this conversation came back headed
  `spencer (auto search, circular)` — XSLOPE's default.

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
`sensitivity.png`, and printed the restored cohesion at the end —
`foundation c now: [800.0]`. **Files…** opens the folder holding the two files.

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
linear at about +0.165 per 100 psf, and the last three rows are identical to six
decimals, 1.244108 each, which it put down to the embankment governing above
600 psf. That last step is the one it did not measure, and it said so rather
than leaving the reader to find out: *"I did not vary the embankment cohesion to
confirm that, so treat the cause as inferred from the flat segment, not
measured."*

Settling it by hand takes one more solve per step.

### Check its work

- **Re-run any of the seven steps**: 0.639688, 0.792223, 0.963558, 1.135100 and
  1.244108 three times over — the reported values to six decimals.
- **Read two rows against the published answers.** The 300 psf row is 0.792,
  which is [LEM-3](lem03_layered_slope.md)'s weak-foundation answer, and the
  800 psf row is its 1.244; the sweep reproduces both ends of that tutorial from
  one request.
- **Pull the critical circle at each step, which is the reading the answer
  flagged as inferred.** Below 600 psf every critical surface bottoms at
  elevation −10, tangent to the rock; at 600 psf and above the search settles on
  (18.50, 43.75) with R = 43.75, bottoming at elevation 0. Across the plateau the
  mechanism has left the foundation, which is why adding foundation strength buys
  nothing there.
- **Open the Materials editor when it finishes.** The foundation reads 800 psf,
  the helper restored the project, and no workbook was written.

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

![The face-slope sweep in the dock: the geometry read back, the loop that moves the crest break and regenerates the circles, the three searched steps, the finished table, and the restored geometry printed at the end](images/w1_sweep_adhoc_1.png){width=560}

It reads the geometry twice — the polygons, the circles and `lem_method None`
first, then the profile lines the polygons come from — and writes its own loop
rather than calling a helper, which is what a study with no mode behind it
needs. Each pass moves the crest break to
x = 2H, 2.5H and 3H for a 20 ft height, holds the toe at (0, 0) and the crest at
elevation 20, runs the flat ground from x = −50 out to 50 ft beyond the break,
calls `resync_geometry()`, and regenerates the starting circles so each face
gets a set built from its own section. It kept a deep copy of the profile lines
and the circles before the first pass, which is what it puts back at the end.

| Face slope | Crest break | FS (Spencer, searched) |
| :---: | :---: | :---: |
| 2:1 | (40, 20) | 1.244 |
| 2.5:1 | (50, 20) | 1.396 |
| 3:1 | (60, 20) | 1.546 |

Its closing account of the study reads the Depth column its own runs printed:

> All three critical circles bottom out at elevation 0.0 — tangent to the top of
> the foundation, not cutting into it. With the foundation (c = 800) twice as
> strong as the embankment (c = 400), the critical mechanism stays in the
> embankment in every case; none of the deeper circles tangent to
> max_depth = -10 governed.

That is the fact LEM-3 makes its own result out of, and every number in the
sentence is on the screen above it. It passed on the admissibility notes for
each face as well: interslice tension on the 2:1 and 2.5:1 runs, and a line of
thrust outside the slice on 18%, 13% and 10% of boundaries.

The restoration is printed rather than promised: the closing snippet prints the
profile lines, the circles and the ground surface back at their original values.
A workbook was still written,
[w1_sweep_adhoc_after.xlsx](files/w1_sweep_adhoc_after.xlsx) — the final
`resync_geometry()` leaves the document dirty, so Studio has an edit to save
even though the model is back where it started.

### Check its work

- **Re-run the three steps**: 1.244146, 1.396287 and 1.545993, the reported
  values to six decimals.
- **Read the critical circle at each step, not only the factor of safety.**
  (18.28, 43.84) with R = 43.84, (25.86, 47.73) with R = 47.73 and (32.98, 53.33)
  with R = 53.33 — every one bottoming at elevation 0, which is what the closing
  paragraph claims.
- **Cross the 3:1 row against the previous section.** The modify conversation
  searched its own 3:1 geometry and returned 1.5460 on the same circle, so two
  conversations agree on that face.
- **Read the 2:1 row against the published answer.** It returns 1.244146 against
  LEM-3's 1.244108 — the same mechanism on a slightly different search. The sweep
  regenerates the starting circles at each step — three on the 2:1 and 2.5:1
  sections, two on the 3:1 — and widens the ground; the saved workbook, back on
  the file's own two circles and its own extents, returns 1.244108.
- **Reload the saved workbook.** The section is back at (0, 0)–(40, 20)–(90, 20)
  over (−30, 0)–(90, 0) with the original two circles, so the restoration is in
  the file and not only in the printout.

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
gave them the classifier's own standing rather than dressing them up as
measurements: *"a last-resort fallback, since the problem states no
stiffnesses"*, classified from each material's strength. The `reason` the tool
prints beside each value carries the rest — E the midpoint of the published range
for that soil type, ν its typical value, and any stated value outranking both.
The turn closed by asking for measured or design stiffnesses if there are any.

The second turn opens by saying the solve is the long part, then runs both
engines. Spencer with a search comes back at **FS = 1.244** on (18.50, 43.75)
with R = 43.75, and `run_fem(analysis='ssrm')` builds its own mesh first — tri6,
six-node triangles, at a target size of 1.2 ft, giving 8,849 nodes and 4,304
elements — then bisects from [1.00, 2.00] in seven steps to a shear strength
reduction factor of safety of **1.254**, 0.010 above Spencer, about 0.8%. The
strength reduction solve took 521 s.

It reported two things about those numbers without being asked. **The bracket
did not close:** the trials at F = 1.2656 and F = 1.2578 hit the
50,000-iteration ceiling with the out-of-balance force still falling, which
counts as inconclusive rather than as failure, so 1.254 is the midpoint of a
final interval of [1.2500, 1.2578] whose upper edge is undecided, and raising
`max_iterations_ceiling` is what would settle it. And **Spencer's two
admissibility notes carry into the comparison** — interslice tension and a line
of thrust outside the slice on 15% of boundaries — passed on as a caveat on the
limit equilibrium number, converged but not fully admissible in the
Duncan-Wright sense. The stiffnesses behind the finite element side stay what
the first turn called them: a last-resort fill.

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
  from the model's declared settings", but this model declares neither an element
  type nor a mesh size; tri6 and the 1.2 ft target are the mesh builder's own
  defaults, the second being the 120 ft ground-surface width over 100. Rebuilding
  at another size is how to find out what the answer owes to it.
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

![Both answers in the dock: the seven methods solved on the circle the file defines and tabulated by the equilibrium each satisfies, the moment argument for the tie at phi = 0, then the pore-pressure answer with its four measured Spencer runs, the restored model, and the three effects a water table would have](images/w1_conceptual_2.png){width=560}

**It refused the premise.** One snippet ran all seven methods on the first of the
two circles the file carries, and the answer opens "They're not *all* the same —
the **moment-equilibrium** methods are identical, and the force-equilibrium ones
differ". It named the basis in the same breath: solved on the surface the model
already defines, Xo = 20, Yo = 40, R = 40, no search.

| Method | Equilibrium satisfied | FS |
| :--- | :--- | :---: |
| Ordinary Method of Slices (OMS) | moment only | 1.2471 |
| Bishop | moment, and vertical force | 1.2471 |
| Spencer | moment and force | 1.2471 |
| Morgenstern-Price | moment and force | 1.2471 |
| Lowe & Karafiath | force only | 1.2988 |
| Janbu (corrected, f₀ = 1.09) | force only | 1.3215 |
| Corps of Engineers | force only | 1.3631 |

The argument under it is the right one. With φ = 0 on a circular arc the base
shear resistance is c·ΔL, which does not depend on the base normal force N at
all; the only thing an interslice assumption changes is how N is distributed, so
moment equilibrium about the center reduces to Σ c·ΔL·R / Σ W·x for every method
that writes it. The force-only methods resolve along an assumed interslice
inclination, which enters the force balance directly rather than through N, so
their assumption survives φ = 0 and spreads them from 1.2988 to 1.3631. It
closed by offering the search it had not run, noting that the surfaces and the
spread can both shift once each method searches on its own.

The second answer keeps `u = none`, which is what the model already carries and
what a total-stress undrained analysis takes: the undrained strength already
contains the pore water, so subtracting a pore pressure would count it twice.
Then it measured that rather than asserting it — four Spencer runs on the first
of the file's two circles, and the model put back at the end:

| Configuration | FS (Spencer, the file's first circle) |
| :--- | :---: |
| `u = none`, no water | 1.2471 |
| Piezometric line at elevation 0, `u = piezo` | 1.2471 |
| Piezometric line at elevation 8, `u = piezo` | 1.3387 |
| Piezometric line at elevation 8, plus γsat = γ + 5 pcf | 1.3309 |

It separated three effects across those four rows, each with a row behind it.
The base pore pressure does nothing at φ = 0: a water surface at elevation 0
sits at or below the ground everywhere and returns 1.2471 bit for bit, because
the strength is c and the term the pore pressure sits in is multiplied by zero.
Raising the surface to elevation 8 puts it above the ground out to about x = 16,
and since `water_loads` ships on `auto` the hydrostatic load on the submerged
ground is derived from the line itself — the whole of the +0.092 is that load
pushing on the toe. And a saturated unit weight costs a little here, 1.3387 down
to 1.3309, because heavier soil adds driving weight without adding
strength at φ = 0. It closed by warning against a buoyant unit weight, since
XSLOPE already carries the water explicitly. One documentation page is cited for
the reading behind both answers: `lem/overview/`.

### Check its work

- **Check that nothing moved.** Neither turn edited the model, and no workbook
  was written; the session printed `restored: ['none', 'none'] []` and re-solved
  to 1.2471 to prove it.
- **Reproduce the four rows.** 1.247119, 1.247119, 1.338699 and 1.330862 on the
  file's first circle, which is the whole of the argument in the answer.
- **Test the attribution rather than reading it.** Re-solve the elevation 8 row
  with the pore pressure option back at `none`: **1.338699**, the same number to
  six decimals. The rise really is the ponded water load and nothing else, which
  is what the answer says and what a wrong answer here would have contradicted.
- **Read the seven methods against the published table.**
  [LEM-3](lem03_layered_slope.md) searches each method its own critical circle
  and reports OMS, Bishop, Spencer and Morgenstern-Price at 1.244 with Lowe at
  1.285, Janbu at 1.313 and Corps at 1.326. The four-way tie and the force
  methods sitting above it are the same finding; the numbers differ because these
  seven are one fixed circle.
- **Solve that circle with Analysis = Single surface**: 1.247, not the 1.244 a
  search returns. The answer states which of the two it computed, which is what
  makes the table readable at all.
- **Follow the documentation page.** `lem/overview/` is a real page, and the
  derivation there is what the answer summarizes. The transcript is
  [w1_conceptual_transcript.md](files/w1_conceptual_transcript.md).

---

## A broken file

We will open a copy of the slope with three transcription errors written into
it, give the assistant the tutorial's own inputs, and ask it to check the file
against them. Nothing says where the faults are or how many there are, and
nothing in the file refuses to solve.

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
This model was built from the LEM-3 tutorial, but the factor of safety looks wrong. The tutorial's inputs are: embankment 130 pcf, c = 400 psf; foundation 135 pcf, c = 800 psf; both undrained; rigid rock 10 ft below the top of the foundation. Check the file against them, fix anything that does not match, and rerun.
```
</div>

![The diagnosis in the dock: the model read once, the baseline search, the snippet that varies one input at a time and re-searches each, the repair, the rerun at 1.244, and the closing table setting each fault against the factor of safety](images/w1_diagnose_1.png){width=560}

It answers by changing things rather than by reading them. After the model and a
baseline search, one snippet varies each fault alone and re-searches — and the
material-assignment row comes back at 1.004, unchanged from the baseline,
because that edit went to the polygons Studio derives from the profile lines and
was rebuilt away. The kernel said so, the session read it, and it re-measured
the assignment on the profile lines themselves before reporting anything:

| Model state | FS (Spencer, searched) |
| :--- | :---: |
| As received | 1.004 |
| Unit weight 13 → 130 pcf only | 0.994 |
| `max_depth` −100 → −10 only | 1.220 |
| Assignment left the file's way, the other two fixed | 1.208 |
| All three fixed | 1.244 |

All three faults are named against the inputs the request supplied. Profile line
1 draws the embankment surface and carried the foundation's Mat ID while the base
line carried the embankment's, which it read as the section built upside down;
the 13 pcf is read as a dropped digit; and `max_depth` = −100 becomes "90 ft of
invented foundation", set back to −10 because the request states the rock sits
10 ft below the top of the foundation. It ranked them too: the false base was
the dominant error because it let the search reach circles that cannot exist
above rock, and the nearly weightless fill was pulling the other way, which is
why the file arrived at a plausible-looking 1.004 rather than at something
obviously wrong. It also restaged the trial circles for the corrected base,
though only one of the two moved: the first became a toe circle at R = 44.72,
Depth = −4.72, and the second was rewritten to the Depth = −10, R = 50 it already
carried. Its account of that step says the old deep circle bottomed 24.6 ft below
the rock line; it bottomed at −10, on the rock line, and −24.58 is where the
*searched* surface of the broken model bottomed.

The rerun returns **FS = 1.244** on (18.51, 42.98) with R = 42.98, bottoming at
elevation 0, with both admissibility notes passed on.

**It repaired the swap the other way round.** Rather than reordering the material
rows, it re-pointed the profile lines' Mat IDs at the rows as they stand. The
section, the strengths along it and the answer are LEM-3's; the material table's
row order is not, and the saved
[w1_diagnose_after.xlsx](files/w1_diagnose_after.xlsx) still lists the foundation
first. A defensible fix, and one to notice before the next person reads that
table top down.

**Two sentences in that summary put a feature at the wrong elevation.** One is
the circle description above; the other closes the answer:

> **FS = 1.244** (Spencer, critical circle Xo = 18.51, Yo = 42.98, R = 42.98,
> Depth = 0.0, entry x = 54.83, exit x = 4.61). The critical circle now bottoms
> at the top of the rock rather than cutting through it.

That `Depth = 0.0` gives the elevation of the surface's lowest point, and the
rock in this model sits at −10, so the circle bottoms at the
embankment/foundation contact, 10 ft above the rock. Every number in that
paragraph is right and the feature named beside them is not, which is the shape
most wrong explanations in these eight sessions take.

<!-- test: file=files/w1_diagnose_start.xlsx, type=circular_search, method=spencer, num_slices=40, expected_fs=1.004, tolerance=0.005 -->
<!-- test: file=files/w1_diagnose_after.xlsx, type=circular_search, method=spencer, num_slices=40, expected_fs=1.244, tolerance=0.005 -->

### Check its work

- **Reload the repaired workbook and search: 1.2441** on (18.51, 42.98) with
  R = 42.98, bottoming at elevation 0. [LEM-3](lem03_layered_slope.md) publishes
  1.2441 on (18.50, 43.75) with R = 43.75; the repaired file searches from its own
  rebuilt circles and lands on a slightly tighter one, the same mechanism to four
  figures.
- **Open its Materials editor before trusting the row order.** The foundation is
  still listed first and the embankment second, with profile line 1 pointing at
  material 2 and line 2 at material 1.
- **Check the base depth in the repaired file.** It reads −10, so all three
  faults are out; the file that came in read −100.
- **Re-read the one-at-a-time table beside the transcript.** The
  material-assignment row is 1.004 in the first pass and 1.208 in the second,
  because the first edit went to the derived polygons and was rebuilt away. The
  1.208 is the measured one, and it is only there because the session re-ran the
  test after reading the kernel's warning.
- **Read the baseline's search notes, which the answer does not quote.**
  Spencer could not solve 79 of the 533 trial surfaces on the broken model and 5
  of those rank below the reported minimum, so the 1.004 the file came in with is
  not a number to carry anywhere.

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

![Both turns in the dock: the search returning 1.244 with its circle laid out in a table, then the one call to generate_report, and the path of the finished document](images/w1_report_2.png){width=560}

The first turn returns **FS = 1.244** on the circle (18.50, 43.75) with R = 43.75
— [LEM-3](lem03_layered_slope.md)'s published answer, circle included. It read
the surface straight off the result, laid `Xo`, `Yo`, `R`, `Depth` and the entry
and exit stations out in a small table, and read the lowest point at elevation 0
as sitting well above the rigid base at −10, so the floor is not controlling the
answer. It reported both admissibility notes and offered a tension crack at the
crest as the usual way to clear them.

The second turn makes one call to `generate_report` and answers in a sentence:
the file name, six figures, and the run the document covers. The document itself
carries six numbered figures and three numbered tables across *1 Traceability*,
*2 Project Definition* and *3 Limit Equilibrium Analysis*, the last holding
3.1 Analysis Inputs, 3.2 Materials, 3.3 Loads and 3.4 Spencer's Method, itself
split into the search, the results, the slice table and the calculations. The
snippet printed the path it was written to, and the reply points at the chat and
at the dock's **Files…** button for opening it. The session saved
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
- **Count the figures against the sentence.** The reply claims six, and the
  document embeds six images; §3.4.1 also reports the search that produced them,
  529 trial surfaces over seven refinement stages.

---

## What it will not do

**It does not convert units.** A model is in one unit system throughout, and
nothing in XSLOPE converts one into another: the numbers entered are the numbers
the engine reads, under whatever system the project declares. Convert a drawing
in meters and kilonewtons yourself before handing it to a project in feet, and
check the entered values against the drawing afterwards.

**An edit can go to a copy the model rebuilds.** A section described by profile
lines carries polygons derived from them, and writing to those polygons is
undone the moment Studio resyncs the geometry: the kernel prints *"The polygon
edit did not take"* and rebuilds. The trap is that a snippet printing its own
work reads the values back before the rebuild, so the printout shows the edit and
the model does not carry it. The diagnosis above walked into this and got out by
re-measuring after reading the warning; two of the other four models reported the
repair as applied and delivered 1.208. Re-running the number after an edit is
what tells the two apart.

**It can name the wrong feature beside a right number.** The diagnosis above
closed on a critical circle that "bottoms at the top of the rock" while printing
`Depth = 0.0` in the same sentence, on a model whose rock sits at −10. Nothing in
the reply reads as a guess, the numbers around the sentence are right, and no
model check has anything to say about it. Reading the elevation the run printed,
against the elevations the model declares, is what finds it.

---

## A harder problem

As the slope and the prompts get more complicated, the opportunity for error
rises. The same eight tasks were also run on the reinforced slope of [LEM-8](lem08_reinforced_slope.md) — a
2 ft cohesive band along the face, six geogrid layers with pullout at both ends,
and a surcharge across the crest — scored the same two ways. Those sessions were
recorded earlier, under different wording of the same requests, so read the two
sets as two problems rather than as one model measured twice. Claude Opus 5 is
again right on every model and every number, at \$2.23 against the \$1.77 the
layered slope cost it, and this time loses three of the second column against the
one it lost on the layered slope. The five models above score 5 to 8 of 8 on
model and numbers there and 4 to 7 of 8 on explanation, and a sixth and smaller
one, Claude Haiku 4.5, gets 2 and 3. Late August 2026, on the same list prices as
the table above:

| Model | Cost | Model and numbers | Explanation |
| :--- | :---: | :---: | :---: |
| Claude Opus 5 | \$2.23 | 8 of 8 | 5 of 8 |
| OpenAI gpt-5.5 | \$1.12 | 5 of 8 | 7 of 8 |
| Kimi K3 | \$1.20 | 5 of 8 | 5 of 8 |
| Claude Sonnet 5 | \$1.31 | 5 of 8 | 6 of 8 |
| Claude Haiku 4.5 | \$0.53 | 2 of 8 | 3 of 8 |
| GLM-5V-Turbo | ~\$0.24 | 5 of 8 | 4 of 8 |

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
  13 turns of work for about \$1.77.
- Building a layered model from a drawing and editing it three ways — face,
  water and strength — with every reported factor of safety reproducing from the
  workbook the session saved.
- Two sweeps, a strength reduction run, two questions answered without touching
  the model, a diagnosis that found and fixed all three planted faults, and a
  report carrying the run it names.
- Checking what came back: reload the workbook and re-run the search, read the
  critical surface and not only the factor of safety, and re-measure any
  mechanism the answer explains — the diagnosis puts a right factor of safety on
  a surface it places 10 ft too low.

**Where to go next:** The [AI Assistant reference](../studio/assistant.md)
documents the helpers the assistant calls, the checks that run after every edit
it makes, and what each provider can do. [LEM-3](lem03_layered_slope.md) builds
this same slope by hand three ways, and [LEM-8](lem08_reinforced_slope.md) builds
the reinforced slope the comparison above runs on.
