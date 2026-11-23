---
title: "What's new in XSLOPE"
description: "What's new in XSLOPE 1.1.0, the free slope stability program: seismic analysis (pseudostatic, two-stage, yield coefficient, and Newmark displacement), a strength ratio material, and the full changelog of every release."
---

# What's New

## XSLOPE 1.1.0

- **Seismic analysis.** A new [Seismic Analysis](seismic/overview.md) section covers the pseudostatic factor of
  safety, now with a vertical seismic coefficient as well as a horizontal one, and the
  [two-stage analysis](seismic/two_stage.md) for a slope that has consolidated. It also covers the
  [yield coefficient](seismic/yield_coefficient.md), on one slip surface or by search, and the
  [Newmark displacement](seismic/newmark.md) of the sliding mass under an acceleration record. Studio runs each
  of them, and the Analysis Report and the AI assistant cover them.
- **Two new seismic tutorials.** [LEM-14](tutorials/lem14_pseudostatic.md) runs a pseudostatic analysis of a
  fill slope, first in one stage and then as a two-stage seismic analysis, and
  [LEM-15](tutorials/lem15_yield_newmark.md) finds the same slope's seismic yield coefficient and its Newmark
  displacement under an earthquake record.
- **A strength ratio material.** The [`ratio` option](lem/overview.md#undrained-strength-ratio) gives an
  undrained strength that grows with the vertical effective stress, as in a normally consolidated clay. Part C of
  [Tutorial LEM-7](tutorials/lem07_strength_envelopes.md) uses it.
- **Seepage on fine meshes.** Steady unconfined seepage converges on fine meshes, in far fewer iterations.

Installers for every release are on the [releases page](https://github.com/njones61/xslope/releases). The full
changelog below lists every change.

## Full Changelog

<!-- changelog -->
