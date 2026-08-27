---
title: "Tutorial W-1 — The AI Assistant"
description: "Studio's built-in assistant put through eight requests on one reinforced slope: building the model from a drawing, editing its geometry and loads, two kinds of parameter study, material properties for a finite element run, a question answered out of the documentation, a diagnosis of a broken file, and the analysis report — with what to check after each one."
---

# Tutorial W-1 — The AI Assistant

Studio carries an assistant that does the work the rest of this series does by
hand. It can build a model from a drawing or from a description in words, edit
any part of one afterwards, run every analysis XSLOPE offers, sweep a parameter
across a range, answer a question out of the built-in documentation, and write
the analysis report — all in the open project, on the same undo stack the
editors use. This tutorial takes one slope and asks for each of those in turn,
and after every request it shows what to look at to decide whether the answer is
right.

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
<div class="tgm-model" markdown>**Model** — the reinforced slope of
[LEM-8](lem08_reinforced_slope.md) —
[xslope_reinforced_slope.xlsx](files/xslope_reinforced_slope.xlsx)</div>
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

![The Assistant settings dialog on a fresh install: Provider set to Claude (Anthropic), Model to claude-opus-5, an empty API key field, a disabled Base URL field, and Confirm before running code checked](images/w01_settings.png){width=460}

The dialog has four fields and a checkbox:

**Provider** — one of five: **Claude (Anthropic)**, **OpenAI**, **Kimi
(Moonshot AI)**, **Z.ai (GLM)**, and **Ollama (local, free)**. Every one of them
can read an image, which is a requirement rather than a coincidence: handing the
assistant a photograph or a sketch of a cross section is the first thing this
tutorial does, and a model that cannot see turns that request into a
conversation about what the picture shows. Where a provider sells text-only
models alongside models that can see, only the second kind is listed. The first
four are hosted services that bill for what you use; Ollama runs a model on your
own machine, free, and sends nothing anywhere.

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

**Base URL** — the address the provider is reached at. It is enabled for Kimi,
Z.ai and Ollama, whose endpoints move — a Z.ai coding-plan key uses a different
base from a standard one, and a local Ollama may be on another port — and
disabled for Claude and OpenAI, whose addresses are fixed.

**Confirm before running code** — checked by default. The assistant works by
writing small Python snippets and running them in Studio's own process, so with
this box checked every snippet is shown to you and waits for your approval
before it executes. Unchecked, the assistant runs code without asking, which
goes faster on a long build. Read what that means before you uncheck it: the
code runs as you, with your file access.

Below the fields, two captions report where the model list came from and what
the selected model can do — whether it supports running code, whether it can
read an image, and whether the provider caches the prompt.

### What it costs

A hosted provider bills per token, and both the tokens read and the tokens
written count. A routine question — one asked of a built model, answered with a
run — costs about **26,000 input tokens** on Claude Opus 5, roughly half of them
served from the provider's prompt cache at a fraction of the price of fresh
input, for a few hundred tokens of reply. Building
a model from a drawing costs more, because it takes more turns. In writing this
documentation, $20 of API credit covered weeks of use.

Models differ in both price and capability, and the two do not track each other
perfectly: a cheaper model may answer a question about the documentation as well
as an expensive one and still fail at a long autonomous build. The line under
the input box reports the tokens each turn spent and the running total for the
conversation:

```
this turn: 26,292 in (12,522 cached) / 443 out · session: 26,292 in / 443 out
```

One request is usually several calls to the model, so that first number climbs
while the assistant works. **New chat** starts the count over.

---

## What to expect

Everything below was produced with **Claude Opus 5**. Another model will do the
same work differently, and so will the same model on another day — the wording
will differ, the code it writes will differ, and it may take a different number
of steps to arrive at the same place. Treat the transcripts as an illustration
of what the requests look like and what comes back, not as a script to match
line for line.

The assistant is good at this work, and everything it produces should still be
checked. Checking is usually easy, and it is the same checking a model built by
hand deserves: read the model checks Studio runs after every edit, look at the
section on the canvas, and do one hand calculation against the number that came
back. Each use case below ends with **Check its work** — what to look at, and
what it should say.

---

## Building a model from a drawing

We will paste the problem drawing from [LEM-8](lem08_reinforced_slope.md) into
an empty project and ask for the model it describes, with nothing else to go on
but the units. One request exercises everything at once here: reading dimensions
off a picture, entering materials, geometry, a load and six reinforcement lines,
and then searching for the critical circle.

<!-- SESSION: w1_build_from_image — the prompt and the drawing, the transcript,
     the dock screenshots, the model it built against the model LEM-8 specifies,
     the factor of safety it returned, and the "Check its work" close. -->

---

## Modifying the model

We will take the finished slope and change three things in one conversation: the
face from 1.25:1 to 2:1, a 500 psf load added across part of the crest, and
every reinforcement line extended 5 ft to the right. Each edit is a request in
plain language against a model already built, and each is rerun so we can watch
the factor of safety move.

<!-- SESSION: w1_modify — the three prompts, the transcript, the canvas after
     each edit, the undo steps each one left, the three factors of safety, and
     the "Check its work" close. -->

---

## A sweep the Parametric dialog offers

We will ask for the geogrid tensile capacity swept from 500 to 3000 lb/ft in six
steps, with a search at every step, plotted against the factor of safety. Studio
has a dialog for exactly this study, so the request tests whether the assistant
reaches for the machinery already there rather than rebuilding it.

<!-- SESSION: w1_sweep_builtin — the prompt, the transcript, the plot it
     produced, the six factors of safety against the same sweep run from the
     Parametric dialog, and the "Check its work" close. -->

---

## A sweep it does not

We will ask for the same slope with two, three, four, five and six geogrid
layers, removing from the top down and searching each time. No dialog offers
that study — the parameter being varied is the number of rows in a table, not a
number in a cell — so the assistant has to write the loop itself.

<!-- SESSION: w1_sweep_adhoc — the prompt, the transcript, the table it
     produced, the layer counts and factors of safety, a check that the right
     layers were removed, and the "Check its work" close. -->

---

## Material properties for a finite element run

We will ask what Young's modulus and Poisson's ratio to use for these two soils,
and why, then have the assistant enter them, mesh the slope and run the strength
reduction analysis. A limit equilibrium model carries no stiffnesses, and the
finite element engine refuses to start without them, so this is the gap every
model crossing from one engine to the other has to close.

<!-- SESSION: w1_elastic_fem — the two prompts, the transcript, the values it
     suggested and its reasoning, the mesh, the strength reduction factor of
     safety against the Spencer answer, and the "Check its work" close. -->

---

## Asking how something works

We will ask two questions that change nothing: how a reliability analysis works
in XSLOPE, and how to choose standard deviations from a handful of tests. The
assistant carries this documentation as context, so a question like this is
answered from the same pages the rest of the site is written on, with a link to
follow.

<!-- SESSION: w1_conceptual — the two prompts, the transcript, the pages it
     linked, the worked problems it pointed at, a check that no edit was made,
     and the "Check its work" close. -->

---

## Finding what is wrong with a model

We will open a copy of the slope with three transcription errors written into it
and say only that the factor of safety is below 1. No request here carries less
information, and none leans harder on the model checks Studio hands back to the
assistant after every edit it makes.

<!-- SESSION: w1_diagnose — the prompt, the transcript, the three errors against
     the three it found, the factor of safety before and after, and the "Check
     its work" close. -->

---

## Generating the report

We will run Spencer with a search and then ask for the analysis report in one
sentence. What comes back is the same document the Report dialog writes, built
from the same run, so this request shortcuts the dialog rather than making a
second kind of report.

<!-- SESSION: w1_report — the two prompts, the transcript, the report it wrote,
     the figures and numbers in it against the run behind it, and the "Check its
     work" close. -->

---

## What it will not do

<!-- PLACEHOLDER: the three limits below are stated from the shipped behavior;
     the measured examples come from the sessions above. -->

**It does not convert units.** A model is in one unit system throughout, and the
assistant works in whatever the project already uses. Give it feet and pounds
and it enters feet and pounds; give it a drawing in meters and kilonewtons for a
project in feet and it will say so rather than convert on its own. Convert the
numbers yourself before you hand them over.

**It asks when a request is ambiguous.** A request that could reasonably mean two
things gets a question back rather than a guess.

**It can misread a dimension.** Reading a number off a drawing leaves more room
for error than anything else in this tutorial — a dimension line that points at
two places at once, a label sitting nearer the wrong feature, a decimal point
lost to image resolution. Nothing downstream will catch it: the model checks
test whether a model is consistent, not whether it matches the picture. That is
why every use case ends with a check, and why the first one checks the built
model against the drawing line by line.

---

## Conclusion

<!-- PLACEHOLDER: the summary bullets are written once the eight sessions have
     run, from what they actually returned. -->

This tutorial covered:

- Setting the assistant up: choosing a provider and model, storing an API key,
  and what a turn costs.
- <!-- one bullet per use case, written from the recorded sessions -->

**Where to go next:** [Building Models Three Ways](building_models.md) puts the
assistant beside the other two ways into a model, and the
[AI Assistant reference](../studio/assistant.md) documents the helpers it calls,
the checks that run after every edit it makes, and what each provider can do.
The [tutorials index](index.md) lists the series.
