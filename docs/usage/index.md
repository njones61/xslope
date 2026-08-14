---
title: "Ways to Use XSLOPE"
description: "The three ways to run XSLOPE — the Studio desktop app, Colab notebooks, and the Python package — and which one to pick."
---

# Ways to Use XSLOPE

XSLOPE can be used three ways. They share the same engine and the same
[Excel input template](input_template.md), so a model built in one moves to
the others unchanged.

**XSLOPE Studio** — the desktop application, and the recommended way to use
XSLOPE for most work. Studio wraps the full engine in a graphical interface:
editors for every input, a live section view, one-click analyses with input
checking, report generation, and an [AI assistant](../studio/assistant.md)
that can build and edit models from a description or a drawing. Installers
for Windows and macOS are on the [Install](../getting_started/install.md)
page, and the [Studio section](../studio/index.md) documents the interface.
The [tutorials](../tutorials/index.md) teach it alongside the other entry
paths.

**The Python package** — `pip install xslope` for scripting, automation, and
anything the interfaces don't expose: batch runs, custom plots, parameter
studies, integration into your own tools. The
[package install page](installation.md) covers setup and the
[API reference](../api/fileio.md) documents the modules.

**Colab notebooks** — run analyses in the browser with nothing installed.
Each notebook loads an input file, offers the analysis options as form
inputs, and plots the results — well suited to quick checks, teaching, and
trying XSLOPE against the [sample problems](../lem/samples.md). See
[Colab Notebooks](notebooks.md).

If you are new to XSLOPE: install [Studio](../getting_started/install.md),
then start with the [first tutorial](../tutorials/lem01_simple_embankment.md).
