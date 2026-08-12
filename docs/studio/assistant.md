---
title: "AI assistant for slope stability modeling — XSLOPE Studio"
description: "The AI assistant built into XSLOPE Studio: build a model from a sketch or a photograph of a cross section, run analyses, script parameter studies, and get results explained in plain language."
---

# AI Assistant

Studio includes a dockable **chat assistant** that drives the app and the engine
with natural language (and images). It goes well beyond the standalone
[Claude Code skill](../usage/claude/index.md): rather than only building an input
file from a sketch, it can perform *any* interaction with the project and the
`xslope` Python API — build or edit inputs, run analyses, script parameter studies,
post-process results, and explain or debug what happened.

Examples of what you can ask:

- *"Add a piezo line 2 ft below the crest and re-run Spencer."*
- *"Vary the slope angle from 20° to 30° and plot FS vs. angle."*
- *"Why did the SSRM not converge?"*
- *"Build this slope from the attached sketch."* (with an image)

![The Assistant dock](images/assistant_dock.png)

---

## Ask it anything

The assistant is also a **tutor**. Not every request is a model: you can ask why a
method behaves the way it does, whether XSLOPE can handle a particular case, where
an input lives, or what a result means — and it answers in conversation, without
building or running anything unless the question needs it. Theory and how-to answers
end with a link into this documentation, and links in the transcript **open in your
browser**, so a chat reply is a starting point rather than the last word. When a
question matches a topic the [verification corpus](../verification/index.md) covers,
it points at the worked problems too.

*The exchanges below are representative, not screenshots. Wording varies with the
model you select; answers are grounded in this documentation and in the model's own
knowledge, and are shortened here to show the pattern rather than to stand in for
the answer you will get.*

**Theory** — *"Why does raising the water table lower my factor of safety?"*

> Strength on the slip surface is effective-stress based: τ_f = c' + (σ_n − u)·tan φ'.
> Pore water carries part of the normal load but contributes no shear strength, so
> raising u shrinks the frictional term. Take a slice with σ_n = 2000 psf, c' = 100 psf,
> φ' = 30°. Dry (u = 0): τ_f = 100 + 2000(0.577) = 1255 psf. Put the water table 10 ft
> above that slice base (u = 62.4 × 10 = 624 psf) and τ_f = 100 + 1376(0.577) = 894 psf
> — 29% less resistance, with the driving weight essentially unchanged. Cohesion is
> untouched, which is why a φ = 0 undrained analysis barely responds to the water table
> at all.
>
> Slice forces and the full derivation:
> [LEM overview](https://xslope.readthedocs.io/en/latest/lem/overview/)

**Capability** — *"Can I run Bishop on a non-circular surface?"*

> No, and it isn't an implementation gap. `capabilities(slope_data)["lem_method"]["bishop"]`
> comes back unavailable with the reason: Bishop sums moments about a circle center, so
> it needs a circular surface — the same is true of OMS. The other five methods (Janbu,
> Corps of Engineers, Lowe & Karafiath, Spencer, Morgenstern-Price) take either family.
> Spencer is the usual choice on a non-circular surface, since it satisfies force and
> moment equilibrium both.
>
> Theory: [Bishop's method](https://xslope.readthedocs.io/en/latest/lem/bishop/)

**How-to** — *"Where do I set a tension crack?"*

> **Global parameters** in Studio, or cells `main!D11` (crack depth) and `main!D12`
> (water depth in the crack) in the workbook; from a script,
> `slope_data['tcrack_depth']` and `['tcrack_water']`. It is a limit-equilibrium
> construction only — the finite element engine ignores it and represents tension
> constitutively, through the material tensile cutoff `t_cut`, so one file with a crack
> depth poses two different problems to the two engines. The input checks warn when a
> crack is entered but no failure surface intersects it while its water thrust still
> applies.
>
> Input checks:
> [preflight](https://xslope.readthedocs.io/en/latest/usage/preflight/)

**Worked problems** — *"How do I model a tieback wall?"*

> Each tieback is a line in the reinforcement table with `Type = tieback`, which presets
> the direction to **axial** (force along the bar, not tangent to the slip surface) and
> the application to **active** (an allowable force, not divided by F). Its capacity is
> the smallest of the bar tension `Tmax`, the pullout developed over the bond length,
> and the connection capacity at the wall face (`Tend1`); soldier piles or a facing that
> carry shear go in the piles table instead. Enter per-element capacities plus a spacing
> and the loader divides to per-unit-width.
>
> The verification corpus has worked tieback problems verified against published
> results:
> [VP49 — retaining wall, grouted tiebacks, soldier piles](https://xslope.readthedocs.io/en/latest/verification/rocscience/#vp49),
> [VP58 — tied-back wall in layered soil](https://xslope.readthedocs.io/en/latest/verification/rocscience/#vp58),
> [VP59 — tieback wall in sand with a drawdown water table](https://xslope.readthedocs.io/en/latest/verification/rocscience/#vp59).
> Theory:
> [reinforcement](https://xslope.readthedocs.io/en/latest/lem/reinforcement/)

**Interpretation** — *"Why did my FS drop when the reservoir drew down?"*

> Drawdown removes the water load that was helping hold the upstream slope up, but pore
> pressures inside a low-permeability fill do not dissipate at the same rate — the slope
> carries something close to its old pore pressures under the new, smaller external
> load. XSLOPE handles this with the three-stage Duncan-Wright-Wong procedure: stage 1
> recovers the pre-drawdown effective stresses on the slip surface, stage 2 computes FS
> after drawdown with undrained strengths derived from those stresses, and stage 3
> re-checks with drained strengths; the reported FS is the lower of stages 2 and 3. So a
> drop is the expected behaviour, not a defect — what is worth checking is whether rapid
> drawdown applies at all, via the time factor T = c_v·t/D².
>
> Procedure and the T rubric:
> [rapid drawdown](https://xslope.readthedocs.io/en/latest/lem/rapid/)

---

## Getting the assistant

The assistant needs the provider library, which ships in the `ai` extra. How you get
it depends on how you installed Studio. In the **packaged app** (`.dmg` / `.msi`) the
provider library is bundled in the installer, so the assistant works out of the box
and there is nothing to add. **From Python / pip**, add the `ai` extra:
`pip install "xslope[gui,ai]"`.

Either way you still choose a provider and enter credentials in **Settings…**
(below). Without the library, Studio runs normally but the dock reports that the
dependency is missing. See [Installation](index.md#installation).

---

## How it works

The assistant is an **agent with code execution**, not a chatbot. Its core tool is
an in-process Python kernel: the model writes small Python snippets, and Studio runs
them with `xslope` imported and the **live project in scope** — `doc` (the project),
`slope_data` (its inputs), and `results` are all preloaded, along with helpers like
`run_lem()` and `sensitivity()`. Most requests reduce to "write and run a little
Python, show me the figure or number," exactly like the notebook workflow.

Because edits go through the same path as the editors, anything the assistant changes:

- renders **immediately** on the canvas,
- lands on the **[undo stack](editing.md#undo-and-redo)** as a labeled step
  (*"Assistant: materials, dloads"*), so you can revert it, and
- invalidates stale results/mesh just like a manual edit.

The assistant **builds into the live document**, not a file. Review the result on
the canvas, then persist it with **Save As**. (An empty project is opened
automatically if none is, so the first snippet works.)

---

## Choosing a model

The assistant is **bring-your-own-model**. A **Settings…** button in the dock opens
a dialog to pick the provider and model and store credentials:

| Provider | Notes |
| --- | --- |
| **Claude** (Anthropic) | Tool use + vision; prompt caching keeps the large system prompt cheap on repeat turns. |
| **OpenAI / GPT** | Tool use + vision. |
| **Ollama** (local) | Runs models on your machine — free, no API key, fully offline. Capability depends on the chosen model. |
| **DeepSeek** | Tool use; vision per model. |
| **Z.ai (GLM)** | Tool use; vision per model. |

![Assistant Settings dialog](images/assistant_settings.png)

API keys are stored in the **OS keychain** (not in plaintext), and the Ollama base
URL is configurable for local models. The dock shows the active provider · model,
and a caption warns when the selected model can't (or may not) run code or accept
images, so the UI degrades gracefully.

### Models

The model list is not fixed when Studio is installed. Opening the dialog — and
pressing **Refresh** — asks the selected provider for its own list of models,
using the key you have stored, and shows what comes back. If the provider can't
be reached, the list falls back to the last one it returned, and then to the
models this version shipped with, so there is always something to choose offline.
A caption under the box says which of the three you are looking at.

A curated set of recommendations is published alongside XSLOPE releases and read
at most once a day. It marks one model per provider as **recommended** — the one
a new install starts on — sorts a few good choices to the top with a one-line
label, and marks models the provider has superseded. Superseded models stay in
the list and can still be chosen.

The box also accepts free text: type any model id the provider takes, including
one released after your copy of Studio, and it is remembered like any other
choice.

---

## Every edit is checked

When the assistant changes the model, the change is validated before you see the
reply. Studio rebuilds the derived geometry and runs the same
[input checks](../usage/preflight.md) the Run dialogs use, then hands the findings
back to the assistant as part of its own result — errors and warnings both. The
assistant is required to resolve what comes back, or to tell you why it stands,
before it reports the model ready.

This matters most for the changes whose consequences are somewhere else. Correct
the base elevation of a model and the circles that were tangent to the old base
are now underneath the new one; nothing about editing that one field says so, and
the model would otherwise look finished until a run failed on it. The checks run
on the edit that caused it, so the stranded circles are named in the same reply.

A read-only question — anything that reads the model without changing it — skips
the checks entirely.

---

## Autonomy: confirm vs. auto

Running model-written code is powerful and carries the same trust model as any
coding agent — the code runs as you, and can read and write files. Studio gives you
an **autonomy mode**, switchable in the dock:

- **Confirm** (default) — you review and approve each code run before it executes.
- **Auto** — the assistant runs code without prompting.

![Confirm-before-run](images/assistant_confirm.png)

Every action is visible in the transcript as a collapsible "ran code" block, and a
**Stop** button halts the agent at any time. Because the kernel runs in Studio's own
process, a bad snippet could hang the app — the confirm mode and Stop button are
your guardrails.

!!! warning "Network egress"
    Hosted models (Claude, OpenAI, DeepSeek, Z.ai) send your prompts and any
    attached images off your machine. **Ollama stays local.** Choose accordingly
    for sensitive work.

---

## Vision and chat UX

![Building from a sketch (vision)](images/assistant_vision.png)

- **Images** — paste or drop an image into the chat (e.g. a hand sketch or a screen
  capture) for models that support vision.
- **Inline figures** — plots the assistant produces appear inline in the transcript.
- **New chat** — starts a fresh conversation and resets the kernel.
- Press **Enter** to send; the transcript renders markdown with collapsible code
  blocks and surfaces actionable error messages.

---

## Relationship to the Claude Code skill

The standalone [`/xslope` Claude Code skill](../usage/claude/index.md) remains a
first-class, file-first workflow for the CLI/IDE. Studio's assistant reuses that
skill's **schema knowledge** (the sheet/category layout, geometry rules, examples)
as system context, but takes a **document-first** path: it populates the live
`slope_data` rather than writing an `.xlsx`. The two coexist — the file-first skill
and the document-first Studio assistant are two front ends to the same engine.
