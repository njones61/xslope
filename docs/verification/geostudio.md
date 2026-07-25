# GeoStudio (SLOPE/W) Verification Corpus

The [GeoStudio slope stability verification manual](https://files.seequent.com/GeoStudio/Manuals/Slope%20Stability%20Verification%20Manual.pdf)
(Seequent) contains 47 verification problems solved with SLOPE/W. Many are the same published benchmarks as the
Rocscience corpus above, which makes the manual doubly useful: its geometry figures are coordinate-labeled where
Slide's sometimes are not (it supplied the ACADS 1(c) interfaces and the inclined ACADS 4 seam above), and its
SLOPE/W numbers give an independent second program to compare against. Where a problem coincides with a built
Rocscience entry, the row below points at the **same XSLOPE input file** — those files are already
regression-tagged in the table above, so no tags are duplicated here. XSLOPE's Janbu carries the fo correction;
SLOPE/W's "Janbu" column is the uncorrected force solution, so those columns are compared via the correction
factor where noted.

Full bibliographic details for the author-year citations on this page are on the
shared [References](references.md) page.

## Building these problems

Most of the entries below were built by transcribing the manual's geometry figures by hand. The few rows that are
not built are held up by a strength model XSLOPE does not yet have (§2.47) or a measured pressure grid with no flow
field to regenerate (§2.10) — not by a missing source document.

That is no longer the only route. Seequent publishes the **models** behind the manual, not just the figures, on a
public CDN — no login, no license:

```
https://files.seequent.com/GeoStudio/SlopeW/<Manual Section Name>.gsz
```

Confirmed for e.g. [ACADS Simple Slope](https://files.seequent.com/GeoStudio/SlopeW/ACADS%20Simple%20Slope.gsz),
[ACADS External Load](https://files.seequent.com/GeoStudio/SlopeW/ACADS%20External%20Load.gsz),
[Prandtl Bearing Capacity](https://files.seequent.com/GeoStudio/SlopeW/Prandtl%20Bearing%20Capacity.gsz),
[Tandjiria Geosynthetic Reinforced Embankment](https://files.seequent.com/GeoStudio/SlopeW/Tandjiria%20Geosynthetic%20Reinforced%20Embankment.gsz),
[Sheahan Amhearst Soil Nails](https://files.seequent.com/GeoStudio/SlopeW/Sheahan%20Amhearst%20Soil%20Nails.gsz) and
[Cannon Dam](https://files.seequent.com/GeoStudio/SlopeW/Cannon%20Dam.gsz).

XSLOPE reads them directly — see [GeoStudio Import/Export](../usage/geostudio.md) — so a problem's exact geometry
and material properties arrive with no transcription step, which removes the ambiguity blocking the remaining
rows:

```python
from xslope.geostudio import read_gsz, list_analyses, gsz_to_slope_data

gsz = read_gsz("ACADS Simple Slope.gsz")
for a in list_analyses(gsz):
    print(a["id"], a["name"], a["method"])

slope_data, caveats = gsz_to_slope_data(gsz, analysis_id=1)
```

That also changes what "verification" can mean here. A **solved** `.gsz` records SLOPE/W's factor of safety for
*every trial surface it evaluated*, not just the critical one — several hundred per analysis — along with the
weight of each sliding mass. So XSLOPE can be run on the *identical* surfaces and the two programs compared with
no difference in search to explain away. See [Importer verification](#importer-verification) below.

!!! note "The model files are Seequent's, and stay Seequent's"
    The `.gsz` files are Seequent's copyrighted material. XSLOPE links to them and reads them; it does not
    redistribute them, and none are committed to this repository. Import them yourself from the links above if
    you want to reproduce this work. The XSLOPE `.xlsx` inputs derived from the published problem data are ours,
    and are published here as usual.

## Importer verification {#importer-verification}

Verifying the *importer* is a different job from verifying the solver, and it needs a different rig. Reading a
`.gsz` back through XSLOPE's own reader proves nothing — a reader and a writer that share the same wrong idea of
the schema agree perfectly with each other. So the importer is scored against **SLOPE/W's own answers**: import
the model, re-solve every trial circle SLOPE/W saved, and compare, with the surface and the method held fixed.
The rig is `tools/gsz_corpus.py` (it takes a path to a local folder of `.gsz` files, which are not
redistributable and so are not in the repo).

Two numbers are reported per analysis, and the split is what makes it diagnostic:

- **weight** — XSLOPE's sliding mass against the weight SLOPE/W recorded for the same surface. This tests the
  imported polygons, the material assignment and the unit weights **with the solver taken out of the way**.
- **FS** — the whole model: strengths, pore pressure, loads, cracks, reinforcement.

A geometry bug moves both. A strength, water or load bug moves only the second.

| Model | Analysis | Surfaces | Weight | FS | Within 1% |
|---|---|---:|---:|---:|---:|
| ACADS Simple Slope | a. simple slope | 355 | +0.19% | −0.10% | 91% |
| ACADS Simple Slope | b. tension crack | 152 | +0.15% | −0.12% | 100% |
| ACADS External Load | a. base case | 116 | +0.03% | +0.50% | 70% |
| Tandjiria | Clay – circular | 210 | +0.11% | −0.10% | 100% |
| Tandjiria | Sand – circular | 338 | +0.19% | −0.07% | 99% |
| Tandjiria | Clay – circular – **reinforced** | 1 | +0.17% | −0.27% | 100% |
| Tandjiria | Sand – circular – **reinforced** | 1 | +0.45% | −0.64% | 100% |
| Rapid drawdown | 2a – after rapid drawdown (SEEP/W) | 2576 | +0.16% | −0.13% | 98% |
| Rapid drawdown | 3a – during slow drawdown (SEEP/W) | 2536 | +0.15% | −0.11% | 97% |

**All 9 scorable analyses agree with SLOPE/W to within 1%** (median 0.12%) over 6285 slip surfaces. Every model's
geometry and unit weights import correctly too — the weight column never exceeds +0.5% anywhere.

The two drawdown analyses are the SEEP/W-coupled ones, and they are scored at **every one of their 11 saved time
steps**, against the pore-pressure field and the trial surfaces SLOPE/W solved at that step. They were −13.2%
until the [SEEP/W field and the reservoir it implies](../usage/geostudio.md#pore-pressure-from-a-seepw-analysis)
were both imported.

Non-circular analyses are reported as not comparable rather than scored. SLOPE/W writes a centre and a radius even
for a block or fully-specified surface, but that circle is *fitted* to the surface, not the surface it solved —
rebuilding it as a circle silently scores the wrong geometry (it reads as a −26% error, which is the rig lying,
not the importer).

??? note "What the corpus caught"
    Every one of these was silent when wrong, and none would have been found by a round-trip test:

    - **Piezometric surfaces index the analysis's *local* point list**, not the shared geometry `<Points>` table.
      Read against the geometry, one model's water table doubled back on itself and the pore pressures came out
      low: **FS +9.7%**. Cannon Dam's water table would have landed 6 m low and sloping the wrong way. The
      exporter had the mirror bug, so a `.gsz` written by XSLOPE lost its water table in GeoStudio.
    - **A tension crack is switched off by *dropping* `<TensionOption>`** — GeoStudio keeps the crack's geometry
      regardless. Treating the leftover points as a live crack put a 2 m water-filled crack into models that had
      none. Invisible in the φ=0 clay analyses (pore pressure cannot touch an undrained strength) and worth
      **−2%** in the c′=0 sand ones. No single file would have shown this; the pair did.
    - **A surcharge's `<Pressure>` is a *unit weight*.** The load is the weight of the fill between the drawn line
      and the ground, so it varies with depth — not a uniform pressure. Worth **+4%**. It looked right at first
      only because the test wedge happened to be exactly 1 m deep, where the two readings coincide.
    - **Reinforcement acts along the bar (axial), not tangent to the slip surface.** Measured, not reasoned:
      axial reproduces SLOPE/W to 0.2% on Tandjiria's reinforced clay, where tangent fails to converge at all.
    - **SLOPE/W derives the reservoir from the SEEP/W head field.** A SEEP/W-coupled analysis records **no water
      surface anywhere** — and still loads the submerged face. The rule it must be using, and is: water stands to
      `y + u/γ_w` at the ground surface. That reproduces SLOPE/W's own per-slice surcharge forces at every time
      step of the drawdown example, receding 627 → 566 → 480 → 363 → 217 → 67 → 0 kN and vanishing at exactly the
      step SLOPE/W's does. Worth **−13.2%** together with the pore pressures — and *most of it was the reservoir,
      not the pressures*. The corpus is what made the split visible: the pore-pressure-only fix left a −7% bias
      that lived only in the partially-submerged steps of the *slow* drawdown, where the water is still receding.

## Corpus status

All 47 sections of the manual. **built** = an XSLOPE input file with verified results; *covered* = the same
problem is built in the [Rocscience corpus](rocscience.md) and regression-tagged there; *partial*, *blocked*,
*planned* as in that corpus. Rows that share a problem with Rocscience link to it, and it links back.

**Completeness.** Where a section cannot be reproduced, the row records why. The one *no lock possible*
row (§2.10 Lanester) prints a measured loading-induced pressure grid with no flow field behind it, so no
seepage solution can regenerate it — a final gap it shares with the Slide2 corpus. The one remaining
*blocked* row (§2.47) is tracked against a named gap XSLOPE does not yet have: a dip-relative
anisotropic/compound strength model. §2.16's level-ground Prandtl mechanism no longer
belongs on that list: the `right_facing` override built for Rocscience VP26 unblocks the flat-arc facing
guard, so it is covered by that shared build. Everything else is built and verified, or covered by the
regression-locked Rocscience build; the corpus is complete relative to what is independently verifiable.

**Manual edition.** The manual tracked here is the **2025.2 edition**. Its 47
two-dimensional problems (chapter 2, "Verifications – 2D", §2.1–2.47) are carried
unchanged from the October 2022 edition, so the section numbering below is valid against
either. The 2025.2 edition also adds an 11-problem "Verifications – 3D" chapter
(§3.1–3.11) verifying SLOPE3D — ellipsoidal and wedge surfaces, 3D seismic and
piezometric cases, and the Kettleman Hills case history. XSLOPE is a 2D formulation, so
the 3D chapter is out of scope and not tracked here.

<div class="corpus-summary" markdown>

| § | Problem | Status | XSLOPE file / results vs SLOPE/W |
|---:|---|---|---|
| 2.1 | ACADS Simple Slope | **built** | [xslope_acads_simple.xlsx](../lem/files/xslope_acads_simple.xlsx): Bishop 0.985 vs SLOPE/W 0.963, Slide 0.987; ACADS reference 1.00. — [details](#gs-2-1) |
| 2.2 | ACADS Tension Crack | **built** | [vp002.xlsx](files/rocscience/vp002.xlsx): Bishop 1.589 / M-P 1.586 vs SLOPE/W 1.664 / 1.660 and Slide 1.596 / 1.592; ACADS reference 1.65-1.70. SLOPE/W sits closer to the ACADS band; the difference is tension-crack water handling and search. — [details](#gs-2-2) |
| 2.3 | ACADS Non-Homogeneous | **built** | [vp003.xlsx](files/rocscience/vp003.xlsx): Bishop 1.403 / M-P 1.371 vs SLOPE/W 1.414 / 1.382; ACADS 1.39. — [details](#gs-2-3) |
| 2.4 | ACADS Non-Homogeneous + Seismic | **built** | [vp004.xlsx](files/rocscience/vp004.xlsx): Bishop 1.013 / M-P 0.987 vs SLOPE/W 1.02 / 0.989; ACADS 1.00. — [details](#gs-2-4) |
| 2.5 | ACADS Talbingo Dam – Dry | **built** | [vp005.xlsx](files/rocscience/vp005.xlsx): 1.955 (all methods, infinite-slope mechanism) vs SLOPE/W 1.951. — [details](#gs-2-5) |
| 2.6 | ACADS Talbingo – Specified Surface | **built** | [vp006.xlsx](files/rocscience/vp006.xlsx): Bishop 2.206 / M-P 2.299 vs SLOPE/W 2.207 / 2.299 — exact. — [details](#gs-2-6) |
| 2.7 | ACADS Weak Layer | **built** | [xslope_acads_weak_layer.xlsx](../lem/files/xslope_acads_weak_layer.xlsx): Spencer 1.258 / M-P 1.248 vs SLOPE/W Bishop 1.269 / M-P 1.261 — [details](#acads-weak-layer). |
| 2.8 | ACADS Weak Layer – Specified Surface | **built** | [vp008.xlsx](files/rocscience/vp008.xlsx): M-P 1.260 vs SLOPE/W 1.261 — exact; Janbu(corr) 1.294 vs SLOPE/W force 1.197 (×fo ≈ 1.29). — [details](#gs-2-8) |
| 2.9 | ACADS External Loading | **built** | [vp009.xlsx](files/rocscience/vp009.xlsx): Spencer 0.724 / Janbu(corr) 0.718 vs SLOPE/W Bishop 0.699 / M-P 0.689 — search-sensitive problem, published spread 0.67-0.81. — [details](#gs-2-9) |
| 2.10 | Lanester Embankment | *no lock possible* | Same problem as [Rocscience #12](rocscience.md): the printed pore-pressure grid is measured loading-induced pressure, not a flow field, so no seepage solution can reproduce it. |
| 2.11 | Arai & Tagyo Homogeneous | **built** | [xslope_arai_tagyo.xlsx](../lem/files/xslope_arai_tagyo.xlsx) vs SLOPE/W Bishop 1.417 / M-P 1.414; A&T 1.451. — [details](#gs-2-11) |
| 2.12 | Arai & Tagyo Pore-Water Pressure | **built** | [vp016.xlsx](files/rocscience/vp016.xlsx): Bishop 1.112 vs Slide 1.118, A&T 1.138 — SLOPE/W reports 1.190, the outlier of the four sources. — [details](#gs-2-12) |
| 2.13 | Greco Layered Slope | **built** | [vp019.xlsx](files/rocscience/vp019.xlsx): circular Spencer 1.429 vs SLOPE/W M-P 1.389, Greco 1.40-1.42. — [details](#gs-2-13) |
| 2.14 | Greco Weak Layer | **built** | [vp020.xlsx](files/rocscience/vp020.xlsx): noncircular Spencer 1.082, circular 1.091 vs SLOPE/W Spencer 1.054, Greco 1.08. — [details](#gs-2-14) |
| 2.15 | Chen & Shao Frictionless Slope | covered | Same problem as [Rocscience #25](rocscience.md#vp25) — [vp025.xlsx](files/rocscience/vp025.xlsx), Prandtl mechanism on a 60° weightless slope: Spencer 1.052 vs SLOPE/W 1.036, Slide 1.051, theory 1.0. — [details](#gs-2-15) |
| 2.16 | Prandtl Bearing Capacity | covered | Same problem as [Rocscience #26](rocscience.md#vp26) — [vp026.xlsx](files/rocscience/vp026.xlsx), the level-ground Prandtl mechanism unblocked by the `right_facing` override: Spencer 1.043 vs theory 1.0, Slide2 0.941. SLOPE/W's own model solves the exact fully-specified (non-circular) mechanism at 0.960 by Morgenstern-Price, read from the file — the `.gsz` importer does not yet rebuild a fully-specified surface, so that number is quoted, not reproduced. — [details](#gs-2-16) |
| 2.17 | Chowdhury & Xu (1995), 5 examples | covered | Same problem as [Rocscience #28](rocscience.md#vp28) — Congress St. Cut + embankment on soft clay, 3 of 10 cases built and reliability-tagged. This 20-analysis `.gsz` is the **provenance source** for VP28's never-printed inputs (sand γ = 21, clay γ = 22, the below-clay-3 `bed` c′ = 200/φ′ = 35) and carries SLOPE/W's own solved FS + PF for all ten cases; XSLOPE reproduces them on the imported circles (Taylor σ_F within ≈1% of SLOPE/W's Monte Carlo). See the [VP28 section](rocscience.md#vp28). |
| 2.18 | Borges & Cardoso Geosynthetic Emb. #2 | **built** | [gs2_18.xlsx](files/geostudio/gs2_18.xlsx) — geosynthetic-reinforced embankment on depth-varying soft clay (SFnDepth → XSLOPE `cp`). M-P 1.153 vs SLOPE/W 1.171 (−1.5%) and Borges & Cardoso 1.15 (+0.3%). — [details](#gs-2-18) |
| 2.19 | Borges & Cardoso Geosynthetic Emb. #3 | covered | Same problem as [Rocscience #32](rocscience.md#vp32) / [RS2 #24](rs2.md#rs2-24) — [vp032a](files/rocscience/vp032a.xlsx) (7 m) / [vp032c](files/rocscience/vp032c.xlsx) (8.75 m): identical materials, foundation and embankment geometry (verified to <1 cm); the reinforcement-friction difference between vendor sources (39.6° vs 31.0°) is immaterial because the fully-embedded bar develops its full 200 kN/m either way. SLOPE/W's own solves (1.229 / 0.972) bracket the vp032 locks (1.218 / 0.981). |
| 2.20 | Probabilistic – Syncrude Dyke | covered | Same problem as [Rocscience #33](rocscience.md#vp33) — El-Ramly et al. (2003), built deterministically (composite surface through the presheared clay-shale). |
| 2.21 | Cannon Dam | covered | Same problem as [Rocscience #34](rocscience.md#vp34) — Clarence Cannon Dam (Wolff & Harr 1987) on the noncircular surface. The SLOPE/W model is a public download (see above). |
| 2.22 | Cannon Dam #2 | covered | Same problem as [Rocscience #35](rocscience.md#vp35) — Hassan & Wolff (1999), the min-β ≠ min-FS benchmark. |
| 2.23 | Li & Lumb – Reliability Index | **built** | [vp036.xlsx](files/rocscience/vp036.xlsx): deterministic Bishop 1.333, β_ln 2.263 vs H&W 1.334 / 2.336; GS reports minimum β 2.04 at FS 1.190 across surfaces. — [details](#gs-2-23) |
| 2.24 | Tandjiria – Geosynthetic Reinforced Emb. | **built** | Same problem as [Rocscience #39](rocscience.md#vp39) — [vp039a-d](files/rocscience/vp039a.xlsx). Also the reinforcement benchmark for the **importer**: on SLOPE/W's own circles the imported geosynthetic reproduces its FS to −0.27% (clay) and −0.64% (sand) — see [Importer verification](#importer-verification). — [details](#gs-2-24) |
| 2.25 | Baker & Leshchinsky – Earth Dam | covered | Same problem as [Rocscience #42](rocscience.md#vp42) — B&L (2001) safety-map clay-core dam. On SLOPE/W's own solved circle XSLOPE's Spencer (1.939) matches SLOPE/W (1.934); on Slide's circle 1.926 vs Slide 1.925, and on Baker's surface 1.882 vs B&L 1.91 — the published cluster reproduced. B&L Fig. 5(a) supplies the phreatic. — [details](#gs-2-25) |
| 2.26 | Baker – Planar Homogeneous | **built** | [gs2_26.xlsx](files/geostudio/gs2_26.xlsx) — Spencer/Janbu 1.352 vs SLOPE/W 1.352 and Baker ≈1.35, essentially exact; the model pins the crest offset at 2.5 m, resolving the Rocscience #43 build's geometry ambiguity. — [details](#gs-2-26) |
| 2.27 | Sheahan – Amherst Soil Nails | covered | Same problem as [Rocscience #47](rocscience.md#vp47) — [vp047.xlsx](files/rocscience/vp047.xlsx), Sheahan & Ho (2003) Amherst test wall: Janbu 0.899 vs Slide 0.890 / Sheahan 0.887. The SLOPE/W model is a public download (see above). — [details](#gs-2-27) |
| 2.28 | Sheahan – Clouterre Test Wall | covered | Same problem as [Rocscience #48](rocscience.md#vp48) — [vp048.xlsx](files/rocscience/vp048.xlsx), Clouterre full-scale test wall, 7 nail rows. |
| 2.29 | Snailz – Reinforced Slope | covered | Same problem as [Rocscience #49](rocscience.md#vp49) — [vp049.xlsx](files/rocscience/vp049.xlsx), SNAILZ soldier-pile tieback wall. |
| 2.30 | Snailz – Geotextile Layers | **built** | [vp050.xlsx](files/rocscience/vp050.xlsx) (Rocscience #50, same SNAILZ model): Janbu(corr) 1.448 vs SLOPE/W force 1.354 (×fo ≈ 1.44) and SNAILZ 1.46; M-P/Spencer 1.576 vs SLOPE/W M-P 1.606. — [details](#gs-2-30) |
| 2.31 | Zhu – Four Layer Slope | **built** | [vp051.xlsx](files/rocscience/vp051.xlsx) (Rocscience #51): Bishop 1.278 vs SLOPE/W 1.284; Spencer 1.294 vs 1.299; M-P 1.304 vs 1.310; Lowe 1.296 vs 1.283; Corps 1.404 vs 1.368. Zhu paper source in `ref_docs_lim_eq/`. — [details](#gs-2-31) |
| 2.32 | Zhu & Lee – Heterogeneous Slope | **built** | [vp052a/b.xlsx](files/rocscience/vp052a.xlsx) (Rocscience #52): wet deep-family Spencer 1.189 matches Slide exactly; see that row. — [details](#gs-2-32) |
| 2.33 | Priest – Rigid Blocks | **built** | [gs2_33.xlsx](files/geostudio/gs2_33.xlsx) — Janbu 1.049 / M-P 1.049 vs Priest (hand) 1.049 and SLOPE/W 1.049, exact. — [details](#gs-2-33) |
| 2.34 | Yamagami – Stabilizing Piles | **built** | [vp054a/b.xlsx](files/rocscience/vp054a.xlsx) (Rocscience #54): no-pile Bishop 1.100 vs SLOPE/W 1.102 — exact; with-pile 1.185 vs SLOPE/W 1.223, Slide 1.193, Yamagami 1.20 (pile-force conventions differ program-to-program). — [details](#gs-2-34) |
| 2.35 | Pockoski & Duncan – Tie-Back Wall | covered | Same problem as [Rocscience #58](rocscience.md#vp58) — P&D (2000) 8-layer tied-back wall. XSLOPE Bishop 1.142 / Spencer 1.140 vs Slide 1.147 / 1.145, UTEXAS4 1.14; the vendor .gsz was saved unsolved. — [details](#gs-2-35) |
| 2.36 | Pockoski & Duncan – Reinforcement | covered | Same problem as [Rocscience #59](rocscience.md#vp59) — P&D (2000) single-row tieback in sand, under-designed (FS < 1). XSLOPE Janbu 0.579 / Corps 0.577 vs SLOPE/W's own Janbu 0.575 / Lowe 0.587. — [details](#gs-2-36) |
| 2.37 | Pockoski & Duncan – Soil Nails | covered | Same problem as [Rocscience #60](rocscience.md#vp60) — P&D (2000) soil-nailed wall. XSLOPE Spencer 1.010 vs SLOPE/W's own 1.000 and Slide 1.009. — [details](#gs-2-37) |
| 2.38 | Loukidis – Seismic Coefficient | **built** | [vp062a/b.xlsx](files/rocscience/vp062a.xlsx) (Rocscience #62): Spencer 1.001 (both cases) vs SLOPE/W 1.00 — exact. — [details](#gs-2-38) |
| 2.39 | Loukidis – Seismic Coefficient #2 | covered | Same problem as [Rocscience #63](rocscience.md#vp63) / [RS2 #68 Case 3](rs2.md#rs2-68) — Loukidis (2003) example 2, three-layer seismic slope. A critical-k꜀ problem: XSLOPE FS 1.001 at the paper's k꜀ = 0.155; critical_kc harness returns k꜀ 0.167 (Spencer) / 0.169 (Bishop). — [details](#gs-2-39) |
| 2.40 | Rapid Drawdown – Walter Bouldin Dam | **built** | Same problem as Slide [VP98](rocscience.md#vp98): xslope DWW 3-stage 1.046 vs SLOPE/W Bishop 1.016 / Spencer 1.02, DWW 1.04. — [details](#gs-2-40) |
| 2.41 | Rapid Drawdown – USACE Benchmark | **built** | [vp096.xlsx](files/rocscience/vp096.xlsx) (Rocscience #96): 3-stage Spencer 1.434 / Bishop 1.432 vs published 1.44. — [details](#gs-2-41) |
| 2.42 | Rapid Drawdown – Pumped Storage Dam | **built** | Same problem as Slide [VP99](rocscience.md#vp99): xslope 1.527 vs SLOPE/W 1.550, DWW 1.56 (geometry re-pinned from this .gsz). — [details](#gs-2-42) |
| 2.43 | Rapid Drawdown – Pilarcitos Dam | **built** | [vp097.xlsx](files/rocscience/vp097.xlsx) (Rocscience #97): Spencer 1.044 / Bishop 1.042. — [details](#gs-2-43) |
| 2.44 | Probability – James Bay Case History | covered (deterministic); prob. is the named future case | Same problem as [Rocscience #75](rocscience.md#vp75) — the James Bay dyke (El-Ramly, Morgenstern & Cruden, after Christian–Ladd–Baecher). The deterministic base case is VP75 (SLOPE/W Bishop 1.46). The 8 Bishop analyses are a **pure spatial-averaging study**: they differ *only* in a `SamplingDistance` (autocorrelation length) — every-slice (0 m) / 30 / 40 / 50 / 80 / 100 m / none — over one set of plain σ's (marine-clay c σ = 8.14, lacustrine c σ = 8.65, fill γ & φ σ = 1). SLOPE/W's own results trace the variance-reduction curve directly: mean FS ≈ 1.46 throughout, **σ_F 0.065 → 0.215** and **PF 0 → 1.45%** as the averaging length grows. Only the "no spatial consideration" point-variance case (σ_F 0.215, PF 1.45%) is reproducible now — the rest need the correlation-length treatment [VP33](rocscience.md#vp33) names. The held `.gsz` is the ready benchmark for that future feature; no new document is required to verify it. |
| 2.45 | Eurocode 7 – Cutting in Clay | **built** | [gs2_45.xlsx](files/geostudio/gs2_45.xlsx) — DA3 partial factors baked into the material (c′* = 8.0, φ′* = 23.04°); Spencer 1.173 vs SLOPE/W ODF 1.174 (−0.07%). — [details](#gs-2-45) |
| 2.46 | Eurocode 7 – Earth Dam | **built** | [gs2_46.xlsx](files/geostudio/gs2_46.xlsx) — DA1-C2 factors + XSLOPE FE seepage; on SLOPE/W's own circle M-P 1.099 vs 1.101 (−0.19%), free-search minimum 1.073 ≈ the Smith textbook 1.07. — [details](#gs-2-46) |
| 2.47 | Compound Strength vs Anisotropic Function | blocked | Needs an orientation-dependent (dip-relative) strength model XSLOPE does not have: AnisotropicFn/CompoundStrength interpolate strength by slice-base angle against a discontinuity dip over angle ranges A/B (the .gsz importer already flags exactly this on import). SLOPE/W: Anisotropic 1.113, Compound 1.118 — solved on the *same* surface these differ by 0.4%, so the gap is the strength model, not the search. The 21-material faulted section has no printed coordinate table (a secondary blocker). |

</div>

## Problem details

### 2.1 — ACADS Simple Slope {#gs-2-1}

The headline ACADS limit-equilibrium benchmark: a simple homogeneous 2:1 slope analyzed with a circular search, with an accepted consensus FOS ≈ 1.00. It shares its XSLOPE input with the Rocscience corpus, so full geometry and material inputs are deferred to the linked VP1 detail.

**Input:** [xslope_acads_simple.xlsx](../lem/files/xslope_acads_simple.xlsx) · **Rocscience detail:** [VP1](rocscience.md#vp1)

![xslope_acads_simple: inputs and representative solution](images/xslope_acads_simple.png)

| Method | XSLOPE | SLOPE/W | Reference |
|---|---|---|---|
| Bishop | 0.985 | 0.963 | Slide 0.987; ACADS 1.00 |

XSLOPE's Bishop FOS (0.985) matches Slide (0.987) and sits close to the ACADS consensus of 1.00, and reads slightly above SLOPE/W (0.963).

**Sources:** GeoStudio SLOPE/W Verification Manual §2.1; Donald & Giam (1989), Giam & Donald (1992).

### 2.2 — ACADS Tension Crack {#gs-2-2}

ACADS problem 1(b): a homogeneous slope with a water-filled tension crack, verifying tension-crack handling in the limit-equilibrium solution. It shares the `vp002.xlsx` XSLOPE input with the Rocscience corpus, so full geometry and material inputs are deferred to the linked Rocscience detail.

**Input:** [vp002.xlsx](files/rocscience/vp002.xlsx) · **Rocscience detail:** [VP2](rocscience.md#vp2)

![vp002: inputs and representative solution](images/vp002.png)

| Method | XSLOPE | SLOPE/W | Reference |
|---|---|---|---|
| Bishop | 1.589 | 1.664 | Slide 1.596 |
| Morgenstern-Price | 1.586 | 1.660 | Slide 1.592 |

SLOPE/W sits closer to the ACADS reference band (1.65–1.70) than XSLOPE or Slide2; the difference traces to tension-crack water handling and search.

**Sources:** GeoStudio SLOPE/W Verification Manual §2.2; Giam & Donald (1989).

### 2.3 — ACADS Non-Homogeneous {#gs-2-3}

ACADS benchmark 1(c): a non-homogeneous three-layer slope analyzed for its critical circular surface. It shares the XSLOPE input (vp003.xlsx) with the Rocscience corpus, so the full geometry and material inputs are deferred to the linked Rocscience detail.

**Input:** [vp003.xlsx](files/rocscience/vp003.xlsx) · **Rocscience detail:** [VP3](rocscience.md#vp3)

![vp003: inputs and representative solution](images/vp003.png)

| Method | XSLOPE | SLOPE/W | Reference |
|---|---|---|---|
| Bishop | 1.403 | 1.414 | 1.39 (ACADS) |
| Morgenstern-Price | 1.371 | 1.382 | 1.39 (ACADS) |

XSLOPE agrees with SLOPE/W to within about 1% for both methods and brackets the ACADS consensus reference of 1.39.

**Sources:** GeoStudio SLOPE/W Verification Manual §2.3; Donald, I.B. & Giam, P. (1989), ACADS.

### 2.4 — ACADS Non-Homogeneous + Seismic {#gs-2-4}

ACADS 1(d): the three-material non-homogeneous slope of §2.3 with a horizontal pseudo-static seismic coefficient of 0.15 added, verifying seismic loading over layered strengths. It shares the XSLOPE input `vp004.xlsx` with the Rocscience corpus, so the full geometry, material properties, and slip surface are documented in the linked Rocscience detail.

**Input:** [vp004.xlsx](files/rocscience/vp004.xlsx) · **Rocscience detail:** [VP4](rocscience.md#vp4)

![vp004: inputs and representative solution](images/vp004.png)

| Method | XSLOPE | SLOPE/W | Reference |
|---|---|---|---|
| Bishop | 1.013 | 1.02 | 1.00 |
| Morgenstern-Price | 0.987 | 0.989 | 1.00 |

XSLOPE's Bishop (1.013) and Morgenstern-Price (0.987) bracket the ACADS consensus of 1.00 and agree with SLOPE/W to within ~1%.

**Sources:** GeoStudio SLOPE/W Verification Manual §2.4; Donald, I.B. & Giam, P. (1989), ACADS suite.

### 2.5 — ACADS Talbingo Dam – Dry {#gs-2-5}

ACADS benchmark 2(a): the four-zone Talbingo Dam at end of construction (dry), searched for the critical circular surface, whose minimum collapses to a shallow infinite-slope mechanism parallel to the steeper upstream face. It shares the XSLOPE input `vp005.xlsx` with the Rocscience corpus, so the full geometry, zone properties, and per-method breakdown are deferred to the linked Rocscience detail.

**Input:** [vp005.xlsx](files/rocscience/vp005.xlsx) · **Rocscience detail:** [VP5](rocscience.md#vp5)

![vp005: inputs and representative solution](images/vp005.png)

| Method | XSLOPE | SLOPE/W | Reference |
|---|---|---|---|
| All methods (infinite-slope mechanism) | 1.955 | 1.951 | — |

Because the critical mechanism is the infinite-slope limit, every method converges to the same XSLOPE factor of safety of 1.955, in close agreement with SLOPE/W's 1.951.

**Sources:** GeoStudio SLOPE/W Verification Manual §2.5; Giam & Donald (1989), ACADS benchmark 2(a).

### 2.6 — ACADS Talbingo – Specified Surface {#gs-2-6}

ACADS benchmark 2(b): the Talbingo Dam, a four-material embankment evaluated on a single specified circular slip surface, verifying that XSLOPE and SLOPE/W agree on a fixed surface. It shares the XSLOPE input vp006.xlsx with the Rocscience corpus, so the full geometry and material inputs are deferred to the linked VP6 detail.

**Input:** [vp006.xlsx](files/rocscience/vp006.xlsx) · **Rocscience detail:** [VP6](rocscience.md#vp6)

![vp006: inputs and representative solution](images/vp006.png)

| Method | XSLOPE | SLOPE/W | Reference |
|---|---|---|---|
| Bishop | 2.206 | 2.207 | — |
| Morgenstern-Price | 2.299 | 2.299 | — |

XSLOPE and SLOPE/W match essentially exactly: Bishop 2.206 vs 2.207 and Morgenstern-Price 2.299 vs 2.299.

**Sources:** GeoStudio SLOPE/W Verification Manual §2.6; ACADS benchmark (Giam & Donald).

### 2.7 — ACADS Weak Layer (non-circular) {#acads-weak-layer}

The ACADS weak-layer case
([SLOPE/W Verification Manual](https://files.seequent.com/GeoStudio/Manuals/Slope%20Stability%20Verification%20Manual.pdf)
sec. 2.7): a 2:1 slope
crossed by a thin low-strength interlayer with a piezometric line at its base.
The critical surface is non-circular, sliding along the weak layer with a back
scarp to the crest — this is the non-circular search verification test. The
ACADS accepted band is FOS ≈ 1.26.

| Property | Soil 1 | Weak layer |
|---|---|---|
| Cohesion, $c'$ (kPa) | 28.5 | 0 |
| Friction angle, $\phi'$ | 20° | 10° |
| Unit weight, $\gamma$ (kN/m³) | 18.84 | 18.84 |

Excel input file: [xslope_acads_weak_layer.xlsx](../lem/files/xslope_acads_weak_layer.xlsx)

![xslope_acads_weak_layer: inputs and representative solution](images/xslope_acads_weak_layer.png)

Results for the methods applicable to non-circular surfaces:

| Method | XSLOPE FOS | Reference | Diff |
|---|---|---|---|
| Spencer | 1.258 | ~1.26 | −0.2% |
| Morgenstern-Price | 1.248 | ~1.26 | −1.0% |
| Corps of Engineers | 1.336 | ~1.26 | +6.0% |
| Lowe & Karafiath | 1.249 | ~1.26 | −0.9% |
| Simplified Janbu | 1.278 | ~1.26 | +1.4% |

The two interior surface points are seeded just above the base of the weak
layer (base $y=26.5$, top $y=27.0$). Because the non-circular search moves
``Horiz`` points horizontally only, that seed elevation *is* the sliding plane:
placing it near the base of a weak interlayer is standard practice and matches
the reference, whereas seeding it at the layer center biases the factor of
safety roughly 1.5% high. With the base placement, XSLOPE's complete-equilibrium
methods land within ~1% of SLOPE/W's Morgenstern-Price value (1.261): Spencer at
1.258 (−0.2%) and Morgenstern-Price (half-sine) at 1.248 (−1.0%). Corps of
Engineers reads modestly
high here, consistent with ground-parallel side-force inclinations on a surface
with a steep back scarp (XSLOPE uses the standard "Corps #2" convention — see
[Force Equilibrium Methods](../lem/force_eq.md)).

**Sources:** GeoStudio [SLOPE/W Verification Manual](https://files.seequent.com/GeoStudio/Manuals/Slope%20Stability%20Verification%20Manual.pdf),
sec. 2.7; Donald, I.B. & Giam, P. (1989), ACADS.

<!-- fs-table -->
**Factor of safety by method** (each method's own critical surface):

| OMS | Bishop | Janbu | Corps | Lowe | Spencer | M-P |
|---:|---:|---:|---:|---:|---:|---:|
| — | — | 1.278 | 1.336 | 1.249 | 1.258 | 1.248 |
<!-- /fs-table -->

<!-- test: file=../lem/files/xslope_acads_weak_layer.xlsx, type=noncircular_search, num_slices=50, fs_janbu=1.278, fs_corps=1.336, fs_lowe=1.249, fs_spencer=1.258, fs_mprice=1.248, benchmark=LEM-2 -->

### 2.8 — ACADS Weak Layer – Specified Surface {#gs-2-8}

The ACADS 3(b) weak-layer slope analyzed on a fully specified non-circular slip surface, verifying XSLOPE against SLOPE/W for a two-material section with a thin weak seam. It shares its XSLOPE input (vp008.xlsx) with the Rocscience corpus, so the full geometry and material inputs are deferred to the linked Rocscience detail.

**Input:** [vp008.xlsx](files/rocscience/vp008.xlsx) · **Rocscience detail:** [VP8](rocscience.md#vp8)

![vp008: inputs and representative solution](images/vp008.png)

| Method | XSLOPE | SLOPE/W | Reference |
|---|---|---|---|
| Morgenstern-Price | 1.260 | 1.261 | — |
| Janbu (corrected) | 1.294 | 1.197 (force, ×fo ≈ 1.29) | — |

XSLOPE's Morgenstern-Price matches SLOPE/W's M-P essentially exactly (1.260 vs 1.261); XSLOPE's corrected Janbu (1.294) reproduces SLOPE/W's uncorrected force solution (1.197) once the Janbu fo correction (×fo ≈ 1.29) is applied.

**Sources:** GeoStudio SLOPE/W Verification Manual §2.8; Donald, I.B. & Giam, P. (1989), *Soil slope stability programs review*, ACADS, Melbourne.

### 2.9 — ACADS External Loading {#gs-2-9}

The ACADS benchmark (Slide #9): a two-material slope with a weak layer, a piezometric water table, and two surcharge strips, solved with a non-circular search. It shares the XSLOPE input vp009.xlsx with the Rocscience corpus, so the full geometry and inputs are deferred to the linked Rocscience detail.

**Input:** [vp009.xlsx](files/rocscience/vp009.xlsx) · **Rocscience detail:** [VP9](rocscience.md#vp9)

![vp009: inputs and representative solution](images/vp009.png)

| Method | XSLOPE | SLOPE/W | Reference |
|---|---|---|---|
| Spencer | 0.724 | — | — |
| Janbu (corrected) | 0.718 | — | — |
| Bishop | — | 0.699 | — |
| Morgenstern-Price | — | 0.689 | — |

This is a search-sensitive benchmark with a wide published spread (roughly 0.67–0.81); XSLOPE's Spencer (0.724) and corrected Janbu (0.718) sit just above the SLOPE/W Bishop (0.699) and Morgenstern-Price (0.689) values, with all four results falling inside the published band.

**Sources:** GeoStudio SLOPE/W Verification Manual §2.9; ACADS benchmark study (Giam & Donald).

### 2.11 — Arai & Tagyo Homogeneous {#gs-2-11}

A homogeneous 20 m, 1.5:1 slope in total-stress soil (Arai & Tagyo 1985, example 1), used to verify the automated critical-circle search against SLOPE/W and the published reference. It shares the XSLOPE input file with the Rocscience corpus, so the full geometry and material inputs are deferred to the linked VP14 detail.

**Input:** [xslope_arai_tagyo.xlsx](../lem/files/xslope_arai_tagyo.xlsx) · **Rocscience detail:** [VP14](rocscience.md#vp14)

![xslope_arai_tagyo: inputs and representative solution](images/xslope_arai_tagyo.png)

| Method | XSLOPE | SLOPE/W | Reference |
|---|---|---|---|
| Bishop's Simplified | see [VP14](rocscience.md#vp14) | 1.417 | 1.451 |
| Morgenstern-Price | see [VP14](rocscience.md#vp14) | 1.414 | 1.451 |

SLOPE/W's Bishop (1.417) and Morgenstern-Price (1.414) agree closely and both sit modestly below the Arai & Tagyo reference of 1.451; XSLOPE's full six-method results for this shared input are tabulated in the linked Rocscience detail.

**Sources:** GeoStudio SLOPE/W Verification Manual §2.11; Arai & Tagyo (1985).

### 2.12 — Arai & Tagyo Pore-Water Pressure {#gs-2-12}

Arai & Tagyo (1985) example 3: a homogeneous slope with a water table, verifying pore-water pressure handling on a circular search. It shares the XSLOPE input `vp016.xlsx` with the Rocscience corpus, so the full geometry and inputs are given in the linked Rocscience detail.

**Input:** [vp016.xlsx](files/rocscience/vp016.xlsx) · **Rocscience detail:** [VP16](rocscience.md#vp16)

![vp016: inputs and representative solution](images/vp016.png)

| Method | XSLOPE | SLOPE/W | Reference |
|---|---|---|---|
| Bishop | 1.112 | 1.190 | Slide 1.118; A&T 1.138 |

XSLOPE's Bishop 1.112 agrees closely with Slide 1.118 and the Arai & Tagyo reference 1.138; SLOPE/W's 1.190 is the outlier of the four sources.

**Sources:** GeoStudio SLOPE/W Verification Manual §2.12; Arai & Tagyo (1985).

### 2.13 — Greco Layered Slope {#gs-2-13}

A four-layer slope with no water (Greco 1996, example 4 / Yamagami & Ueta 1988), verifying the critical-surface search against a shallow non-circular benchmark optimum. It shares the vp019.xlsx XSLOPE input with the Rocscience corpus, so full geometry and material inputs are deferred to the linked Rocscience detail.

**Input:** [vp019.xlsx](files/rocscience/vp019.xlsx) · **Rocscience detail:** [VP19](rocscience.md#vp19)

![vp019: inputs and representative solution](images/vp019.png)

| Method | XSLOPE | SLOPE/W | Reference |
|---|---|---|---|
| Spencer / M-P | 1.429 | 1.389 | Greco 1.40–1.42 |

XSLOPE's circular Spencer (1.429) agrees with SLOPE/W's Morgenstern–Price solution (1.389) and brackets the Greco reference range (1.40–1.42).

**Sources:** GeoStudio SLOPE/W Verification Manual §2.13; Greco (1996), Yamagami & Ueta (1988).

### 2.14 — Greco Weak Layer {#gs-2-14}

A four-layer slope with a 0.5 m weak seam running along the inclined model base beneath a water table, verifying the search for a noncircular slip surface through a weak layer. It shares the XSLOPE input vp020.xlsx with the Rocscience corpus, so full geometry and inputs are deferred to the linked Rocscience detail.

**Input:** [vp020.xlsx](files/rocscience/vp020.xlsx) · **Rocscience detail:** [VP20](rocscience.md#vp20)

![vp020: inputs and representative solution](images/vp020.png)

| Method | XSLOPE | SLOPE/W | Reference |
|---|---|---|---|
| Spencer (noncircular) | 1.082 | 1.054 | Greco 1.08 |
| Spencer (circular) | 1.091 | — | — |

XSLOPE's noncircular Spencer solution (1.082) matches the published Greco value (1.08) and runs slightly above SLOPE/W's Spencer result (1.054); the circular Spencer surface gives a marginally higher 1.091.

**Sources:** GeoStudio SLOPE/W Verification Manual §2.14; Greco (1996) ex. 5; Chen & Shao (1988).

### 2.15 — Chen & Shao Frictionless Slope {#gs-2-15}

The classical Prandtl bearing mechanism on a weightless, frictionless 60° slope under a critical strip load, verifying limit-equilibrium against the theoretical FS = 1.0. This is the same problem as the Rocscience corpus, sharing XSLOPE input [vp025.xlsx](files/rocscience/vp025.xlsx); see the linked Rocscience detail for the full geometry, materials, and analytically constructed slip surface.

**Input:** [vp025.xlsx](files/rocscience/vp025.xlsx) · **Rocscience detail:** [VP25](rocscience.md#vp25)

![vp025: inputs and representative solution](images/vp025.png)

| Method | XSLOPE | SLOPE/W | Reference |
|---|---|---|---|
| Spencer | 1.052 | 1.036 | Slide 1.051; theory 1.0 |

XSLOPE's Spencer FS of 1.052 matches Slide (1.051) and sits just above SLOPE/W (1.036), all close to the theoretical value of 1.0 for this frictionless mechanism.

**Sources:** GeoStudio SLOPE/W Verification Manual §2.15; Chen & Shao (1988).

### 2.16 — Prandtl Bearing Capacity {#gs-2-16}

The classical Prandtl bearing mechanism on level ground: a weightless (γ ≈ 0), c = 20 kPa, φ = 0 soil under a strip surcharge of 102.83 kPa — exactly c·N<sub>c</sub>, so the theoretical factor of safety is 1.0 by construction. Both ground crossings of the slip surface sit at the same elevation, which used to trip the flat-arc facing guard; the `right_facing` override built for this problem (Rocscience VP26) resolves that, so the surface — an active wedge, a log-spiral/circular fan, and a passive wedge — solves cleanly. This is the same problem as the Rocscience corpus, sharing XSLOPE input [vp026.xlsx](files/rocscience/vp026.xlsx); see the linked Rocscience detail for the full geometry and slip-surface construction.

SLOPE/W's own file was also checked directly: its "Fully Specified" analysis solves the identical, non-circular Prandtl surface at FS = 0.960 (Morgenstern-Price) — read straight from the model's saved results via `read_gsz_results`. That number is quoted rather than reproduced, because the `.gsz` importer does not yet rebuild a fully-specified (non-circular) surface into an XSLOPE slip surface; that is a separate, narrower gap than the flat-arc guard, and unrelated to it.

**Input:** [vp026.xlsx](files/rocscience/vp026.xlsx) · **Rocscience detail:** [VP26](rocscience.md#vp26)

| Method | XSLOPE | SLOPE/W | Reference |
|---|---|---|---|
| Spencer | 1.043 | 0.960* | theory 1.0; Slide2 0.941 |

\* SLOPE/W's own Morgenstern-Price solve of its fully-specified surface, not an XSLOPE re-solve — see above.

XSLOPE's Spencer (1.043) and SLOPE/W's own Morgenstern-Price (0.960) bracket the theoretical FS of 1.0 from
opposite sides, about 4% apart each way — the same interslice-convention spread already documented for VP26
against Slide2 (0.941), not a new discrepancy.

**Sources:** GeoStudio SLOPE/W Verification Manual §2.16; Prandtl (1921).

### 2.18 — Borges & Cardoso – Geosynthetic Embankment #2 {#gs-2-18}

Borges & Cardoso (2002)'s Case 2: a 5 m geosynthetic-reinforced embankment (c′ = 0,
φ′ = 35°, γ = 20) over four soft-clay layers whose undrained strength increases with
depth. The two deeper clays are GeoStudio "Su-varies-with-depth" (SFnDepth) materials,
which map exactly onto XSLOPE's `cp` option (Su = c + cp·max(0, r_elev − y)): Clay 3 is
16 + 0.96·z and Clay 4 is 18.4 + 1.314·z below their layer tops (the mapped rates
reproduce the manual's printed layer-bottom strengths to the digit). The geosynthetic
(Tmax = 200 kN/m, interface friction 33.7°, unanchored) is laid at y = 1.0 across the
embankment base as an axial geosynthetic; on the critical surface it is fully embedded,
so the full 200 kN/m develops and the factor of safety is insensitive to the bond
length. This is also [Slide2 VP31](rocscience.md) (the same Borges & Cardoso Case 2,
materials matching to rounding), which is *covered* by this build; the full geometry
is in the builder (`benchmarks/geostudio/build_gs2_18.py`).

**Input:** [gs2_18.xlsx](files/geostudio/gs2_18.xlsx)

![gs2_18: inputs and representative solution](images/gs2_18.png)

| Method | XSLOPE | SLOPE/W | Reference |
|---|---|---|---|
| Morgenstern-Price | 1.153 | 1.171 | Borges & Cardoso 1.15 |
| Bishop | 1.154 | 1.170 | — |
| Janbu | 1.327 | 1.233 | — |

Run on SLOPE/W's own critical circle (which a grid search confirms is also XSLOPE's
critical), XSLOPE's Morgenstern-Price 1.153 sits 1.5% below SLOPE/W's 1.171 and 0.3%
above the Borges & Cardoso (2002) published 1.15; Bishop tracks at −1.4%. Janbu is
higher because XSLOPE reports the f₀-corrected value where the manual's Janbu column is
uncorrected. Without the geosynthetic the M-P factor falls to 1.011, so the
reinforcement carries a large share of the resistance.

**Sources:** GeoStudio SLOPE/W Verification Manual §2.18; Borges & Cardoso (2002).

<!-- test: file=files/geostudio/gs2_18.xlsx, type=single_circle, num_slices=60, fs_bishop=1.155, fs_spencer=1.148, fs_mprice=1.153, fs_janbu=1.330, benchmark=GS-2.18 -->

### 2.23 — Li & Lumb – Reliability Index {#gs-2-23}

The Li & Lumb (1987) / Hassan & Wolff (1999) reliability benchmark: a homogeneous slope with an r_u pore-pressure ratio for which both the deterministic Bishop factor of safety and the lognormal reliability index β are computed from the variability in c′, φ′, and γ. It shares the vp036.xlsx input with the Rocscience corpus, so the full geometry and material statistics are deferred to the linked Rocscience detail.

**Input:** [vp036.xlsx](files/rocscience/vp036.xlsx) · **Rocscience detail:** [VP36](rocscience.md#vp36)

![vp036: inputs and representative solution](images/vp036.png)

| Method | XSLOPE | SLOPE/W | Reference |
|---|---|---|---|
| Bishop (deterministic FS) | 1.333 | 1.190 | 1.334 |
| Reliability index β_ln | 2.263 | 2.04 | 2.336 |

XSLOPE's deterministic Bishop FS (1.333) and its β on that surface (2.263) match the Hassan & Wolff reference (1.334 / 2.336); SLOPE/W instead searches for the minimum reliability index across surfaces, reporting β 2.04 at FS 1.190, so its β and FS are not evaluated on the deterministic surface.

**Sources:** GeoStudio SLOPE/W Verification Manual §2.23; Li & Lumb (1987), Hassan & Wolff (1999).

### 2.24 — Tandjiria – Geosynthetic Reinforced Emb. {#gs-2-24}

Tandjiria (2002)'s required-reinforcement half-embankment on soft clay, evaluated as a clay fill and as a sand fill with a geosynthetic at the embankment base. It shares the XSLOPE input (vp039a) with the Rocscience corpus, so the full geometry and inputs are deferred to the linked Rocscience detail; here the section's distinct role is as the geosynthetic-reinforcement benchmark for the SLOPE/W importer.

**Input:** [vp039a.xlsx](files/rocscience/vp039a.xlsx) · **Rocscience detail:** [VP39](rocscience.md#vp39)

![vp039a: inputs and representative solution](images/vp039a.png)

| Method | XSLOPE | SLOPE/W | Reference |
|---|---|---|---|
| Clay fill, reinforced (imported, on SLOPE/W's own circles) | −0.27% | baseline | — |
| Sand fill, reinforced (imported, on SLOPE/W's own circles) | −0.64% | baseline | — |

Run on SLOPE/W's own critical circles, the imported geosynthetic reproduces SLOPE/W's reinforced factor of safety to within −0.27% (clay fill) and −0.64% (sand fill), isolating the reinforcement handling with no search difference to explain away.

**Sources:** GeoStudio SLOPE/W Verification Manual §2.24; Tandjiria (2002).

### 2.25 — Baker & Leshchinsky – Earth Dam {#gs-2-25}

Baker & Leshchinsky (2001)'s safety-map clay-core dam: granular fill (c′ = 0, φ′ = 40°,
γ = 21.5) around a diamond core (c′ = 20, φ′ = 20°, γ = 20) on a hard base (c′ = 200,
φ′ = 45°), a half-full upstream reservoir, a phreatic surface dropping through the core
to the downstream toe, and a 5-m cracked crest layer modeled as a dry tension crack. It
shares the `vp042.xlsx` XSLOPE input with the Rocscience corpus, so the full geometry and
inputs are deferred to the linked Rocscience detail — where this problem's
verification is documented in full. The XSLOPE input is tiled directly as the
material-zone polygons of this model's own `.gsz` region set.

**Input:** [vp042.xlsx](files/rocscience/vp042.xlsx) · **Rocscience detail:** [VP42](rocscience.md#vp42)

![vp042: inputs and representative solution](images/vp042.png)

| Method | XSLOPE (rigorous statics) | SLOPE/W (own solve) | Reference |
|---|---|---|---|
| Spencer, Slide's critical circle | 1.926 | — | Slide 1.925 |
| Spencer, SLOPE/W's own circle | 1.939 | 1.934 | — |
| Spencer, Baker's surface | 1.882 | — | Baker & Leshchinsky (2001) 1.91 |

XSLOPE reproduces the tightly clustered published references on all three surfaces. On
SLOPE/W's own solved critical circle (read from the `.gsz`) the two programs agree to
within 0.005 in Spencer FS on the identical surface, geometry, and water, with XSLOPE's
total sliding-mass weight ≈ 56,020 against SLOPE/W's 56,127. Both sides use total unit
weight + phreatic pore pressure (confirmed verbatim from B&L 2001 and by inspection of
this model's own materials); the reservoir-load statics are exact against a buoyant-weight
oracle, and XSLOPE's phreatic matches this model's own piezometric line to within ~0.2 m
along the failure surface. The XSLOPE values are regression-locked as VP42. The B&L (2001)
paper — held in the reference library — supplies the phreatic geometry the vendor manuals
leave unlabeled (its Fig. 5(a)) and reports the dam's global minimum Fmin = 1.91 by
Spencer's method, computed with total unit weight and pore pressure taken from the vertical
distance to the phreatic surface.

**Sources:** GeoStudio SLOPE/W Verification Manual §2.25; Baker & Leshchinsky (2001).

### 2.26 — Baker – Planar Homogeneous {#gs-2-26}

Baker (2001)'s planar-slip-surface benchmark: a homogeneous, dry c′-φ′ slope (H = 10 m,
face angle 76.0°, c′ = 30, φ′ = 30°, γ = 20) evaluated on planes through the toe, with
FS plotted against the daylight point's normalized position X = x/H on the backslope.
The critical plane sits at X = 0.85. The SLOPE/W model pins the exact geometry (crest
offset 2.5 m) — which also resolves the Rocscience corpus's VP43 analog, whose build
inferred 3 m and read FS ≈ 1.43 against the ≈ 1.35 references.

**Input:** [gs2_26.xlsx](files/geostudio/gs2_26.xlsx)

![gs2_26: inputs and representative solution](images/gs2_26.png)

| Method | XSLOPE | SLOPE/W | Reference |
|---|---|---|---|
| Spencer / Janbu | 1.352 | 1.352 | Baker & Leshchinsky (2001) ≈ 1.35 |

XSLOPE's Spencer and Janbu match SLOPE/W's own solved value (1.352, on the identical
toe-to-(8.5, 10) plane; sliding-mass weight 600 kN/m agrees exactly) to within 0.01%,
and the Baker & Leshchinsky reference to 0.15%. XSLOPE's Morgenstern-Price declines
this surface: on a single straight plane α is constant for every slice, and the
unconstrained λ-search reaches equilibrium only by driving ~70% of the interslice
forces into tension — the solver's admissibility guard rejects that as unphysical
(Corps and Lowe-Karafiath hit the same guard), leaving Spencer as the rigorous method
that converges.

**Sources:** GeoStudio SLOPE/W Verification Manual §2.26; Baker (2001); Baker &
Leshchinsky (2001).

<!-- test: file=files/geostudio/gs2_26.xlsx, type=single_noncirc, num_slices=50, fs_spencer=1.352, fs_janbu=1.352, benchmark=GS-2.26 -->

### 2.27 — Sheahan – Amherst Soil Nails {#gs-2-27}

The Amherst test wall — a 6 m vertical cut in undrained clay retained by two rows of soil nails and a shotcrete facing, which failed in the field test — evaluated on planar surfaces through the toe. It reuses the same XSLOPE input as the Rocscience corpus, so the full geometry, nail capacities, and loading are deferred to the linked Rocscience detail.

**Input:** [vp047.xlsx](files/rocscience/vp047.xlsx) · **Rocscience detail:** [VP47](rocscience.md#vp47)

![vp047: inputs and representative solution](images/vp047.png)

| Method | XSLOPE | SLOPE/W | Reference |
|---|---|---|---|
| Janbu | 0.899 | — | Slide 0.890; Sheahan 0.887 |

XSLOPE's Janbu FS of 0.899 agrees closely with Slide's 0.890 and Sheahan's trial-wedge 0.887; the SLOPE/W model itself is a public download rather than a tabulated FS here.

**Sources:** GeoStudio SLOPE/W Verification Manual §2.27; Sheahan & Ho (2003).

### 2.30 — Snailz – Geotextile Layers {#gs-2-30}

A reinforced-slope problem from the Caltrans SNAILZ reference manual — a multi-row nail/geotextile-reinforced wall evaluated on a predefined deep wedge surface — verifying XSLOPE's reinforcement handling against SLOPE/W and the published SNAILZ solution. It shares the `vp050.xlsx` input with the Rocscience corpus, so full geometry and reinforcement inputs are deferred to the linked VP50 detail.

**Input:** [vp050.xlsx](files/rocscience/vp050.xlsx) · **Rocscience detail:** [VP50](rocscience.md#vp50)

![vp050: inputs and representative solution](images/vp050.png)

| Method | XSLOPE | SLOPE/W | Reference |
|---|---|---|---|
| Janbu (corrected) | 1.448 | 1.354 (force, ×fo ≈ 1.44) | SNAILZ 1.46 |
| Morgenstern-Price / Spencer | 1.576 | 1.606 (M-P) | — |

XSLOPE's Janbu (corrected) 1.448 matches SLOPE/W's uncorrected force solution 1.354 once the fo correction (×fo ≈ 1.44) is applied and agrees with the SNAILZ 1.46 reference; the Morgenstern-Price/Spencer pair agrees at 1.576 vs SLOPE/W's 1.606.

**Sources:** GeoStudio SLOPE/W Verification Manual §2.30; Caltrans SNAILZ reference manual.

### 2.31 — Zhu – Four Layer Slope {#gs-2-31}

Zhu, Lee & Jiang's four-layer slope with a water table, a 5 m dry tension crack, and seismic loading (k=0.1), analyzed on a specified circle. It shares the vp051.xlsx input with the Rocscience corpus, so the full geometry, material properties, and phreatic-line calibration are deferred to the linked Rocscience detail.

**Input:** [vp051.xlsx](files/rocscience/vp051.xlsx) · **Rocscience detail:** [VP51](rocscience.md#vp51)

![vp051: inputs and representative solution](images/vp051.png)

| Method | XSLOPE | SLOPE/W | Reference |
|---|---|---|---|
| Bishop | 1.278 | 1.284 | — |
| Spencer | 1.294 | 1.299 | — |
| Morgenstern-Price | 1.304 | 1.310 | — |
| Lowe-Karafiath | 1.296 | 1.283 | — |
| Corps of Engineers | 1.404 | 1.368 | — |

XSLOPE tracks SLOPE/W within a few thousandths across the rigorous methods (Bishop, Spencer, Morgenstern-Price), with the wider spread confined to the Lowe-Karafiath and Corps side-force assumptions.

**Sources:** GeoStudio SLOPE/W Verification Manual §2.31; Zhu, Lee & Jiang (2003).

### 2.32 — Zhu & Lee – Heterogeneous Slope {#gs-2-32}

The Zhu & Lee heterogeneous benched slope (four materials, water table, tension crack), verifying that an unconstrained circular search lands in the governing deep failure family. This shares its XSLOPE input with the Rocscience corpus (Rocscience #52) — see the linked detail for full geometry and inputs.

**Input:** [vp052a.xlsx](files/rocscience/vp052a.xlsx) · **Rocscience detail:** [VP52](rocscience.md#vp52)

![vp052a: inputs and representative solution](images/vp052a.png)

| Method | XSLOPE | SLOPE/W | Reference |
|---|---|---|---|
| Spencer (wet, deep family) | 1.189 | — | Slide 1.189 |

XSLOPE's wet deep-family Spencer FS of 1.189 matches Slide exactly.

**Sources:** GeoStudio SLOPE/W Verification Manual §2.32; Zhu & Lee (2002).

### 2.33 — Priest – Rigid Blocks {#gs-2-33}

Priest (1993), *Discontinuity Analysis for Rock Engineering*: a classical rock-slope
planar failure — a single rigid block translating on one straight failure plane
(despite the section title, not a multi-block assembly) — with a 15 m tension crack at
the crest, 25% water-filled. One Mohr-Coulomb material: c′ = 20, φ′ = 30°, γ = 25. The
face rises at 60° from the toe (0, 0) to the crest (17.32, 30) with a flat crest
beyond. The failure plane dips at 30° from the toe; XSLOPE specifies it untruncated (to
its natural daylight point) and lets `tcrack_depth = 15` truncate it, which reproduces
SLOPE/W's own solved endpoint (25.98, 15.0) exactly.

**Input:** [gs2_33.xlsx](files/geostudio/gs2_33.xlsx)

![gs2_33: inputs and representative solution](images/gs2_33.png)

| Method | XSLOPE | SLOPE/W | Reference |
|---|---|---|---|
| Janbu | 1.049 | 1.049 | Priest (hand) 1.049 |
| Morgenstern-Price | 1.049 | 1.049 | — |

All XSLOPE methods agree to four decimals on this surface and Janbu's f₀ correction is
exactly 1.00 — expected for a single straight plane, and the same behaviour the
manual's FS-vs-λ figure shows for SLOPE/W (moment and force factors coincide regardless
of the interslice assumption). Total sliding weight matches SLOPE/W's to 0.004%.

**Sources:** GeoStudio SLOPE/W Verification Manual §2.33; Priest (1993).

<!-- test: file=files/geostudio/gs2_33.xlsx, type=single_noncirc, num_slices=50, fs_janbu=1.049, fs_mprice=1.049, fs_spencer=1.049, benchmark=GS-2.33 -->

### 2.34 — Yamagami – Stabilizing Piles {#gs-2-34}

A homogeneous slope stabilized by a row of micro-piles, verifying pile-reinforced limit equilibrium against the unreinforced baseline. It shares the XSLOPE input with the Rocscience corpus (VP54), so full geometry and material inputs are deferred to the linked Rocscience detail.

**Input:** [vp054a.xlsx](files/rocscience/vp054a.xlsx) · **Rocscience detail:** [VP54](rocscience.md#vp54)

![vp054a: inputs and representative solution](images/vp054a.png)

| Method | XSLOPE | SLOPE/W | Reference |
|---|---|---|---|
| Bishop (no pile) | 1.100 | 1.102 | — |
| Bishop (with piles) | 1.185 | 1.223 | Slide 1.193; Yamagami 1.20 |

The unreinforced case matches SLOPE/W essentially exactly (1.100 vs 1.102), while the reinforced FS spreads because pile-force conventions differ program-to-program.

**Sources:** GeoStudio SLOPE/W Verification Manual §2.34; Yamagami.

### 2.35 — Pockoski & Duncan – Tie-Back Wall {#gs-2-35}

Pockoski & Duncan (2000)'s tied-back excavation wall: a 44-ft wall in eight horizontal
soil layers (granular and cohesive fills over organic silt, an over-consolidated crust,
marine clays and glaciomarine deposits), water table at grade in front and el. 102.5
behind, retained by three tieback rows at 20° whose capacity is bond-governed at 40,000
lb/ft of wall. It shares the `vp058.xlsx` XSLOPE input with the Rocscience corpus, so the
full geometry, layer strengths and anchor layout are deferred to the linked Rocscience
detail.

**Input:** [vp058.xlsx](files/rocscience/vp058.xlsx) · **Rocscience detail:** [VP58](rocscience.md#vp58)

![vp058: inputs and representative solution](images/vp058.png)

| Method | XSLOPE | SLOPE/W | Reference |
|---|---|---|---|
| Bishop | 1.142 | — (file saved unsolved) | Slide 1.147; UTEXAS4 1.14; P&D SLOPE/W 1.14 |
| Spencer | 1.140 | — | Slide 1.145; UTEXAS4 1.14 |

The vendor `.gsz` was saved unsolved, so SLOPE/W's own critical factor of safety is not in
the file; the SLOPE/W reference is the value Pockoski & Duncan themselves reported for it
(1.14). XSLOPE's Bishop and Spencer (1.142 / 1.140) sit within the published cluster. The
P&D (2000) report — held in the reference library — supplies the eight-layer section,
strengths and anchor capacities that de-stale this row.

**Sources:** GeoStudio SLOPE/W Verification Manual §2.35; Pockoski & Duncan (2000).

### 2.36 — Pockoski & Duncan – Reinforcement {#gs-2-36}

Pockoski & Duncan (2000)'s single-row tieback wall in homogeneous sand (c′ = 0, φ′ = 30°)
with the water table drawn down to the wall face — under-designed on purpose, so every
published factor of safety is below 1. It shares the `vp059.xlsx` XSLOPE input with the
Rocscience corpus, so the full geometry and inputs are deferred to the linked Rocscience
detail.

**Input:** [vp059.xlsx](files/rocscience/vp059.xlsx) · **Rocscience detail:** [VP59](rocscience.md#vp59)

![vp059: inputs and representative solution](images/vp059.png)

| Method | XSLOPE | SLOPE/W (own solve) | Reference |
|---|---|---|---|
| Janbu simplified | 0.579 | 0.575 | Slide 0.583; P&D SLOPE/W 0.61 |
| Corps / Lowe-Karafiath | 0.577 | 0.587 (Lowe) | Slide 0.588 |

XSLOPE's force-family lock (Janbu 0.579, Corps 0.577) matches SLOPE/W's own current solve
(Janbu 0.575, Lowe-Karafiath 0.587) to within 1–2% — tighter than the historical
inter-program spread this problem was built to expose (published Bishop values alone span
0.56–0.74). Spencer and Morgenstern-Price are inadmissible on the prescribed non-circular
surface in XSLOPE (base normals near the wall go into tension), so the force-equilibrium
methods carry the lock; SLOPE/W's own search rides a circular critical surface where its
Spencer converges at 0.564 and Bishop at 0.531. The P&D (2000) report supplies the sand
strengths, drawn-down water table and anchor layout.

**Sources:** GeoStudio SLOPE/W Verification Manual §2.36; Pockoski & Duncan (2000).

### 2.37 — Pockoski & Duncan – Soil Nails {#gs-2-37}

Pockoski & Duncan (2000)'s 25-ft soil-nailed wall in undrained sandy clay (c = 800 psf,
φ = 0) under a 250-psf crest surcharge plus a 500-psf strip, with a dry 7-ft tension crack
and five passive nail rows at 15°. It shares the `vp060.xlsx` XSLOPE input with the
Rocscience corpus, so the full geometry, surcharges and nail capacities are deferred to
the linked Rocscience detail.

**Input:** [vp060.xlsx](files/rocscience/vp060.xlsx) · **Rocscience detail:** [VP60](rocscience.md#vp60)

![vp060: inputs and representative solution](images/vp060.png)

| Method | XSLOPE | SLOPE/W (own solve) | Reference |
|---|---|---|---|
| Spencer | 1.010 | 1.000 | Slide 1.009; P&D SLOPE/W 1.02 |
| Janbu simplified | 1.043 | — | Slide 1.041; P&D SLOPE/W 1.07 |
| Bishop | — | 0.995 | — |

XSLOPE's Spencer (1.010) matches SLOPE/W's own solved Spencer (1.000) to 1% and Slide
(1.009) essentially exactly; the nailed wall is designed for a factor of safety near 1.
The P&D (2000) report supplies the surcharge pattern, tension-crack depth and five-row
nail layout that de-stale this row.

**Sources:** GeoStudio SLOPE/W Verification Manual §2.37; Pockoski & Duncan (2000).

### 2.38 — Loukidis – Seismic Coefficient {#gs-2-38}

A homogeneous slope loaded pseudo-statically at its critical seismic coefficient (dry and ru = 0.5 cases), where the factor of safety is driven to ≈ 1. It shares the XSLOPE vp062a/b inputs with the Rocscience corpus, so the full geometry and inputs are deferred to the linked Rocscience detail.

**Input:** [vp062a.xlsx](files/rocscience/vp062a.xlsx) · **Rocscience detail:** [VP62](rocscience.md#vp62)

![vp062a: inputs and representative solution](images/vp062a.png)

| Method | XSLOPE | SLOPE/W | Reference |
|---|---|---|---|
| Spencer | 1.001 | 1.00 | — |

XSLOPE's Spencer factor of safety is 1.001 in both cases, matching SLOPE/W's 1.00 exactly.

**Sources:** GeoStudio SLOPE/W Verification Manual §2.38; Loukidis, Bandini & Salgado (2003).

### 2.39 — Loukidis – Seismic Coefficient #2 {#gs-2-39}

Loukidis, Bandini & Salgado (2003)'s second example: a three-layer dry slope — a light
c = 4 kPa cap (φ = 30°, γ = 17), a weak c = 25 kPa / **φ = 15°** middle band (γ = 19) that
the mechanism rides, and a strong φ = 45° base (c = 15, γ = 19). The target is not a
factor of safety but the **critical seismic coefficient** k꜀, the horizontal pseudo-static
coefficient at which the searched minimum FS = 1. It shares the Loukidis example-2 geometry
with the Rocscience and RS2 corpora, so the full geometry is deferred to the linked details.

**Input:** [vp063.xlsx](files/rocscience/vp063.xlsx) · **Rocscience detail:** [VP63](rocscience.md#vp63) · **RS2 detail (critical_kc):** [RS2-68 Case 3](rs2.md#rs2-68)

![vp063: inputs and representative solution](images/vp063.png)

| Target | XSLOPE | SLOPE/W | Reference |
|---|---|---|---|
| FS at the paper's k꜀ = 0.155 (Spencer) | 1.001 | — (file saved unsolved) | Loukidis 1.000 (by definition of k꜀) |
| Critical k꜀ (Spencer / Bishop) | 0.167 / 0.169 | — | Loukidis Spencer 0.155, FEM 0.161, UB 0.172 / LB 0.148; Slide2 0.151 / 0.155 |

Two locks cover this problem. [VP63](rocscience.md#vp63) fixes k at the paper's 0.155 and
confirms the slope is just stable there (noncircular Spencer 1.001). The
[RS2-68 Case 3](rs2.md#rs2-68) `critical_kc` bisection harness instead solves for k꜀
directly and returns 0.167 (Spencer) / 0.169 (Bishop) — ~10% above the paper's 0.155
because the governing surface rides the thin φ = 15° band, which is intrinsically
non-circular, and a circular search cannot follow it as tightly as the log-spiral. The
XSLOPE values still fall inside the paper's rigorous upper/lower-bound bracket [0.148,
0.172] and sit on the reference FEM (0.161). The vendor `.gsz` was saved unsolved, so
SLOPE/W's own k꜀ is not in the file. The Loukidis (2003) paper — held in the reference
library — supplies the three-layer section (Fig. 9(a): a 15 m upper block over a 20 m
lower face, φ = 15° band) and Table 3's k꜀ values that de-stale the "interface anchors
ambiguous" note.

**Sources:** GeoStudio SLOPE/W Verification Manual §2.39; Loukidis, Bandini & Salgado (2003).

### 2.40 — Rapid Drawdown – Walter Bouldin Dam {#gs-2-40}

The Walter Bouldin Dam is a rolled earthfill dam that failed during a rapid drawdown in 1975; this case verifies XSLOPE's Duncan-Wright-Wong three-stage rapid drawdown procedure. It shares the `vp098.xlsx` input with the Rocscience corpus, so the full geometry, zones, and undrained envelopes are deferred to the linked Rocscience detail.

**Input:** [vp098.xlsx](files/rocscience/vp098.xlsx) · **Rocscience detail:** [VP98](rocscience.md#vp98)

![vp098: inputs and representative solution](images/vp098.png)

| Method | XSLOPE | SLOPE/W | Reference |
|---|---|---|---|
| DWW 3-stage (Spencer) | 1.046 | 1.02 | 1.04 |
| DWW 3-stage (Bishop) | — | 1.016 | — |

XSLOPE's DWW three-stage FS (1.046, Spencer circular search) agrees with SLOPE/W's Spencer (1.02) and the published DWW value (1.04) to within about 1%, and brackets SLOPE/W's Bishop result (1.016).

**Sources:** GeoStudio SLOPE/W Verification Manual §2.40; Duncan, Wright & Wong (1990).

### 2.41 — Rapid Drawdown – USACE Benchmark {#gs-2-41}

A homogeneous embankment dam rapid-drawdown benchmark from USACE EM 1110-2-1902 (2003) Appendix G, evaluated with the Duncan-Wright-Wong 3-stage procedure on the specified circle. It shares the XSLOPE input vp096.xlsx with the Rocscience corpus, so full geometry and material inputs are deferred to the linked Rocscience detail.

**Input:** [vp096.xlsx](files/rocscience/vp096.xlsx) · **Rocscience detail:** [VP96](rocscience.md#vp96)

![vp096: inputs and representative solution](images/vp096.png)

| Method | XSLOPE | SLOPE/W | Reference |
|---|---|---|---|
| Bishop | 1.432 | — | 1.44 |
| Spencer | 1.434 | — | 1.44 |

XSLOPE's 3-stage Spencer (1.434) and Bishop (1.432) both agree with the published USACE benchmark of 1.44 to within about half a percent.

**Sources:** GeoStudio SLOPE/W Verification Manual §2.41; USACE EM 1110-2-1902 (2003) Appendix G (Duncan, Wright & Wong 3-stage procedure).

### 2.42 — Rapid Drawdown – Pumped Storage Dam {#gs-2-42}

Rapid drawdown (285 → 120 ft) of a hypothetical pumped-storage dam — silty-clay core and free-draining rockfill shells — analyzed with the Duncan-Wright-Wong 3-stage procedure. It shares its XSLOPE input with the Rocscience corpus (VP99), so the full geometry and material inputs are deferred to the linked Rocscience detail.

**Input:** [vp099.xlsx](files/rocscience/vp099.xlsx) · **Rocscience detail:** [VP99](rocscience.md#vp99)

![vp099: inputs and representative solution](images/vp099.png)

| Method | XSLOPE | SLOPE/W | Reference |
|---|---|---|---|
| Spencer (DWW 3-stage) | 1.527 | 1.550 | 1.56 |

The geometry is re-pinned from this very .gsz (read with `xslope.geostudio.read_gsz`); the earlier eyeball trace of Slide's unlabeled figure had made the dam ≈19 ft short and left XSLOPE ~7% low. With the vendor geometry XSLOPE reads 1.527, inside the Slide / SLOPE/W / DWW band (1.53–1.56).

**Sources:** GeoStudio SLOPE/W Verification Manual §2.42; Duncan, Wright & Wong (1990).

### 2.43 — Rapid Drawdown – Pilarcitos Dam {#gs-2-43}

A three-stage rapid-drawdown analysis of Pilarcitos Dam — a homogeneous earthfill embankment (Duncan, Wright & Wong 1990) drawn down from 72 to 37 ft, a documented drawdown-failure case. It verifies XSLOPE's staged drawdown procedure and shares its XSLOPE input with the Rocscience corpus, so the full geometry and strength envelopes are deferred to the linked Rocscience detail (VP97).

**Input:** [vp097.xlsx](files/rocscience/vp097.xlsx) · **Rocscience detail:** [VP97](rocscience.md#vp97)

![vp097: inputs and representative solution](images/vp097.png)

| Method | XSLOPE | SLOPE/W | Reference |
|---|---|---|---|
| Bishop | 1.042 | — | — |
| Spencer | 1.044 | — | — |

XSLOPE's Spencer (1.044) and Bishop (1.042) staged factors of safety are essentially identical; the geostudio entry reports the shared VP97 XSLOPE run, with the published Slide and Duncan–Wright–Wong reference values documented in the linked Rocscience detail.

**Sources:** GeoStudio SLOPE/W Verification Manual §2.43; Duncan, Wright & Wong (1990).

### 2.45 — Eurocode 7: Cutting in Clay {#gs-2-45}

A 1:2 cutting in a homogeneous clay (characteristic c′ = 10, φ′ = 28°, γ = 20) with a
water table and a 35 kPa permanent crest surcharge, checked to Eurocode 7 Design
Approach 3 (after the *Designers' Guide to Eurocode 7*). XSLOPE has no native
partial-factor feature, so the DA3 material factors (set M2) are baked into the
authored material — c′* = 10/1.25 = 8.0 kPa, φ′* = atan(tan 28°/1.25) = 23.04° — with
soil weight and the permanent surcharge at their factor of 1.0. An ordinary analysis of
the factored model then reproduces GeoStudio's Overdesign Factor (ODF).

**Input:** [gs2_45.xlsx](files/geostudio/gs2_45.xlsx)

![gs2_45: inputs and representative solution](images/gs2_45.png)

| Method | XSLOPE | SLOPE/W (ODF) | Reference |
|---|---|---|---|
| Spencer | 1.173 | 1.174 | book (Bishop) 1.193 |
| Bishop | 1.172 | 1.173 | — |

XSLOPE's Spencer ODF matches SLOPE/W's to −0.07%, and on SLOPE/W's own critical circle
the two agree to −0.09% — the factored-parameter emulation reproduces the EC7 design
check at method level. Janbu is not compared (XSLOPE reports the f₀-corrected value;
the manual's column is uncorrected).

**Sources:** GeoStudio SLOPE/W Verification Manual §2.45; *Designers' Guide to
Eurocode 7*.

<!-- test: file=files/geostudio/gs2_45.xlsx, type=circular_search, num_slices=40, fs_spencer=1.173, fs_bishop=1.172, benchmark=GS-2.45 -->

### 2.46 — Eurocode 7: Earth Dam {#gs-2-46}

A homogeneous clay dam (characteristic c′ = 12, φ′ = 20°, γ = 19.2) with an upstream
reservoir at 5.1 m total head, analysed as coupled steady-state seepage + slope
stability under Eurocode 7 Design Approach 1, Combination 2 (after *Smith's Elements of
Soil Mechanics*, ex. 5.12). Pore pressures come from an XSLOPE finite-element seepage
solve (u = 'seep' with mesh/solution sidecars) whose phreatic surface reproduces
SEEP/W's to ~0.1–0.2 m; the DA1-C2 material factors (set M2) are baked in — c′* = 9.6
kPa, φ′* = 16.23° — so an ordinary analysis yields the Overdesign Factor.

**Input:** [gs2_46.xlsx](files/geostudio/gs2_46.xlsx) (+ `gs2_46_mesh.json`,
`gs2_46_seep.csv`)

![gs2_46: inputs and representative solution](images/gs2_46.png)

| Method | XSLOPE | SLOPE/W (ODF) | Reference |
|---|---|---|---|
| Morgenstern-Price | 1.073 | 1.091 | book (Bishop) 1.07 |
| Bishop | 1.074 | 1.092 | — |

On SLOPE/W's own critical circle XSLOPE's Morgenstern-Price (1.099) matches SLOPE/W
(1.101) to −0.2%; XSLOPE's free search finds a slightly more critical downstream
surface at 1.073, landing on the Smith textbook value (1.07), where SLOPE/W reports its
base-truncated composite minimum at 1.091 — a search-coverage difference, not a model
one. Janbu is not compared (corrected vs uncorrected convention).

**Sources:** GeoStudio SLOPE/W Verification Manual §2.46; *Smith's Elements of Soil
Mechanics* (8th ed.), ex. 5.12.

<!-- test: file=files/geostudio/gs2_46.xlsx, type=circular_search, num_slices=40, fs_mprice=1.073, fs_bishop=1.074, fs_spencer=1.074, benchmark=GS-2.46 -->

## Transient seepage (SEEP/W) {#transient-seepage}

The rows above are SLOPE/W limit-equilibrium problems. GeoStudio also ships a **SEEP/W**
example library, and a handful of those examples are pure (uncoupled) transient-seepage
verifications — the same physics XSLOPE's [transient solver](../seep/transient.md)
(`run_transient_seepage`) implements: two-dimensional saturated/unsaturated flow with a
storage term, `div(kr K grad h) + Q = S dh/dt`. The ones built here range from the clean,
small, closed-form-anchored consolidation/infiltration columns (T01, T02) through the
two-dimensional reservoir-drawdown dam (T03) that exercises the falling-series-head
reservoir face, the two-dimensional clay-lined pond (T04) whose water table rises through
an exit face, and the leach column (T05) driven by a specified-flux boundary. Each is
locked with `type=tseep_head` tags (the same transient head-profile check the
[Rocscience groundwater](rocscience_groundwater.md) corpus uses): the solver marches the
transient problem and the total head is sampled at named points at a save time. The
stepped-suction multistep-outflow column (T07) exercises the time-varying **head**
(plain-Dirichlet) series. One further example (T06) is out of the current transient
envelope and is documented as blocked below.

As with the SLOPE/W corpus, the `.gsz` model files are Seequent's and are not committed;
their solved per-timestep `node.csv` results are read as an oracle and the compared
values distilled into the tables and locks below.

| Example | Physics tested | Status | XSLOPE file / result |
|---|---|---|---|
| Simulating consolidation with SEEP/W | saturated storage `Ss` (Terzaghi) | **built** | [gs2_cons.xlsx](files/geostudio/gs2_cons.xlsx): centre $u_e$ within 0.02 kPa of Terzaghi at t = 150/604/1460 s — [details](#seepw-t01) |
| Verification – infiltration into dry soil | unsaturated storage `C(ψ)` + `kr(ψ)` (van Genuchten) | **built** | [gs2_infil.xlsx](files/geostudio/gs2_infil.xlsx): wetted-zone head within 0.05 m of SEEP/W at t = 46 800 s — [details](#seepw-t02) |
| Rapid drawdown | falling-reservoir series head + potential seepage face on a dam (unsaturated) | **built** | [gs2_rdd_inst.xlsx](files/geostudio/gs2_rdd_inst.xlsx) / [gs2_rdd_slow.xlsx](files/geostudio/gs2_rdd_slow.xlsx): interior total head tracks SEEP/W within ~0.3–0.6 m through the 30-day drawdown — [details](#seepw-t03) |
| Leakage from pond with clay liner | unconfined 2-D water-table rise through an exit face (linear mesh) | **built** | [gs2_pond.xlsx](files/geostudio/gs2_pond.xlsx): interior head tracks SEEP/W within ~0.1 m mid-fill, ~0.35 m at the near-steady leaking state — [details](#seepw-t04) |
| Mineral heap leaching | specified-**flux** (Neumann) top BC + unsaturated storage (1-D column) | **built** | [gs2_heap.xlsx](files/geostudio/gs2_heap.xlsx): column head within ~0.04 m of SEEP/W at the IC/early frames, ~0.12 m at the high-rate near-steady — [details](#seepw-t05) |
| Infiltration into multi-layered system | ponded infiltration with a per-layer imposed IC + unit-gradient base (1-D, 14 layers) | *blocked* | out of the transient envelope — the imposed non-steady per-layer initial condition cannot come from the lock runner's steady IC, and the unit-gradient base BC is not in the solver's BC set — [details](#seepw-t06) |
| GeoStudio-PEST – Multistep Outflow | stepped-**suction** base head via a time-varying **head** (plain-Dirichlet) series (unsat drainage) | **built** | [gs2_mso.xlsx](files/geostudio/gs2_mso.xlsx): the drained column reaches the stepped base suction at each stage, reproducing SEEP/W's −0.07…−0.22 m pressure field — [details](#seepw-t07) |

### SEEPW-T01 — Simulating consolidation with SEEP/W {#seepw-t01}

The SEEP/W example that reproduces **Terzaghi one-dimensional consolidation** as an
uncoupled seepage run: a 0.05 m saturated clay column (K = 1×10⁻⁸ m/s,
m_v = 0.005 /kPa, so S_s = m_v γ_w = 0.04905 /m and c_v = K/S_s = 2.04×10⁻⁷ m²/s)
carries a uniform initial excess pore pressure u₀ = 10 kPa and drains at both ends
(double drainage, drainage path H = 0.025 m). It is the ideal test of the **saturated
storage term**: the problem is linear, so there is a closed-form target and no
unsaturated nonlinearity to confound it.

The excess pore pressure is carried on a datum offset of 100 m (the same device the
[GW15–21](rocscience_groundwater.md) consolidation ports use): the total head is held
uniform at t = 0 (the t = 0 steady solve sets the uniform initial condition) and both
faces step down for t > 0, so the excess head follows Terzaghi Eq. 17.3 and the excess
pressure γ_w(h − 100) is directly comparable to SEEP/W's `node.csv` in kPa.

**Input:** [gs2_cons.xlsx](files/geostudio/gs2_cons.xlsx)

![SEEPW-T01: Terzaghi consolidation isochrones](images/gs2_cons.png)

| t (s) | U (Terzaghi) | Terzaghi centre $u_e$ | XSLOPE centre | SEEP/W centre |
|---:|---:|---:|---:|---:|
| 150 | 25% | 9.97 kPa | 9.96 kPa | 9.86 kPa |
| 604 | 50% | 7.78 kPa | 7.78 kPa | 7.86 kPa |
| 1460 | 75% | 3.93 kPa | 3.95 kPa | 4.77 kPa |

XSLOPE sits on the Terzaghi closed form to within 0.02 kPa at every point and time — the
locked values are the analytical excess heads. SEEP/W is a second, independent oracle: it
agrees at early time but lags the closed form by t = 1460 s (its own reported degree of
consolidation is 24/48/69 % against Terzaghi's 25/50/75 %), because it runs only ten
exponential time steps, where XSLOPE steps finely. This is visible in the figure as the
green (t = 1460 s) SEEP/W dots standing off the analytical curve while XSLOPE's open
circles lie on it.

**Sources:** GeoStudio SEEP/W example "Simulating Consolidation with SEEP/W"; Terzaghi
(1943), *Theoretical Soil Mechanics*, Eq. 17.3.

<!-- test: file=files/geostudio/gs2_cons.xlsx, type=tseep_head, target_size=0.001, time=150, max_head_change_frac=0.02, points=0.005:0.0125:100.9073;0.005:0.025:101.0165;0.005:0.0375:100.9073, tolerance=0.02, benchmark=SEEPW-CONS-t150 -->
<!-- test: file=files/geostudio/gs2_cons.xlsx, type=tseep_head, target_size=0.001, time=604, max_head_change_frac=0.02, points=0.005:0.0125:100.5683;0.005:0.025:100.7928;0.005:0.0375:100.5683, tolerance=0.02, benchmark=SEEPW-CONS-t604 -->
<!-- test: file=files/geostudio/gs2_cons.xlsx, type=tseep_head, target_size=0.001, time=1460, max_head_change_frac=0.02, points=0.005:0.0125:100.2834;0.005:0.025:100.4008;0.005:0.0375:100.2834, tolerance=0.02, benchmark=SEEPW-CONS-t1460 -->

### SEEPW-T02 — Verification: infiltration into dry soil {#seepw-t02}

The canonical hard-infiltration benchmark, and the **unsaturated** counterpart to T01: a
1 m column, initially dry at a uniform pressure head ψ = −8 m, is ponded at the surface
(ψ = 0) and a sharp wetting front advances downward over 46 800 s. It exercises the parts
of the transient path T01 does not — the van Genuchten moisture-capacity storage C(ψ)
(θ_s = 0.363, θ_r = 0.186 so the drainable porosity S_y = 0.177; GeoStudio A = 9.81 kPa
maps to α = 1.0 /m, n = 1.53) and the van Genuchten–Mualem relative permeability kr(ψ),
on Ksat = 1×10⁻⁶ m/s. The dry initial state is set without an explicit initial-head
field: a uniform ψ with a unit downward gradient is itself a steady state, so holding the
top and bottom at ψ = −8 at t = 0 makes the steady initial solve return it, and for
t > 0 the top steps to the ponded head.

**Input:** [gs2_infil.xlsx](files/geostudio/gs2_infil.xlsx)

![SEEPW-T02: infiltration front vs SEEP/W](images/gs2_infil.png)

The primary oracle is SEEP/W's own `node.csv` pressure field at the final time (the
published external reference is the Warrick, Lomen & Yates (1985) semi-analytical
profile). XSLOPE reproduces the **wetted zone** behind the front — the physically
meaningful, water-bearing part of the profile — to within 0.05 m of head:

| y (m) | SEEP/W ψ | XSLOPE ψ | Δ head |
|---:|---:|---:|---:|
| 0.6 | −0.166 m | −0.215 m | 0.049 m |
| 0.7 | −0.064 m | −0.082 m | 0.018 m |
| 0.8 | −0.020 m | −0.025 m | 0.005 m |
| 0.9 | −0.003 m | −0.003 m | 0.000 m |

The mid-front (ψ = −4 m) crossing lands at y ≈ 0.41 in XSLOPE against y ≈ 0.37 in
SEEP/W — a ~0.04 m offset that is the expected **lumped-mass vs consistent-mass front
diffusion**: XSLOPE uses a lumped HRZ mass matrix, which damps the front oscillations
SEEP/W's own dry-soil write-up warns about at the cost of a slightly smeared front
position. The deep dry zone below the front differs by up to 0.3 m of head for a
different reason: SEEP/W holds ψ = −8 there (its explicit initial condition, essentially
frozen since kr is tiny), where XSLOPE's steady-solve initial condition relaxes to
hydrostatic; on the flat dry tail of the retention curve this barely changes the water
deficit the front must fill, which is why the front position and wetted profile still
agree. The lock is on the four wetted-zone points above, at a tolerance (0.08 m) set
deliberately wider than the saturated groundwater page's norm (0.05 m) to carry the
front-diffusion offset.

**Sources:** GeoStudio SEEP/W example "Verification – Infiltration into Dry Soil";
Warrick, Lomen & Yates (1985), *Soil Sci. Soc. Am. J.* 49.

<!-- test: file=files/geostudio/gs2_infil.xlsx, type=tseep_head, target_size=0.01, time=46800, max_head_change_frac=0.01, points=0.025:0.6:0.4344;0.025:0.7:0.6358;0.025:0.8:0.7803;0.025:0.9:0.8967, tolerance=0.08, benchmark=SEEPW-INF -->

### SEEPW-T03 — Rapid drawdown {#seepw-t03}

The flagship for XSLOPE's **reservoir feature set**. A silty-clay embankment (base
x = 3–47 m, crest 10 m tall, upstream slope carrying a reservoir to el. 8, a free-draining
toe drain below the downstream toe) sits at steady state under a full reservoir, then is
drained two ways: **instantaneously** (the reservoir is removed at t = 0, head steps
8 → 0) and **slowly** (head falls 8 → 0 m linearly over 5 days). Both run 30 days. This is
exactly the pair the **submerged-only Dirichlet series head** path was built for: a node
on the upstream face is held at the reservoir head h(t) only while submerged
(y ≤ h(t)) and flips to a potential seepage (exit) face once the falling water line drops
below it, so the emerging upstream face and the trapped-then-dissipating interior
pressures are resolved by the same active-set the steady unconfined solver uses.

The dam fill is a silty clay (Ksat = 1×10⁻⁶ m/s, θ_s = 0.4, θ_r = 0.05 so S_y = 0.35,
m_v = 1×10⁻⁴ /kPa so S_s = 9.81×10⁻⁴ /m). The vendor's SampleMatls "Silty Clay" retention
curve is mapped to van Genuchten by a least-squares fit of its 20 tabulated points
(suction → pressure head), giving α = 0.338 /m, n = 1.85 (RMS 0.007 in effective
saturation). The toe drain is ~11× more permeable (Ksat = 1.157×10⁻⁵ m/s) and its boundary
is pinned at total head 0. The initial condition is the t = 0 steady solve with the
reservoir series at 8 m (the repeated-time step-series idiom the [GW15–21](rocscience_groundwater.md)
and T01/T02 ports use).

**Input:** [gs2_rdd_inst.xlsx](files/geostudio/gs2_rdd_inst.xlsx) (instantaneous),
[gs2_rdd_slow.xlsx](files/geostudio/gs2_rdd_slow.xlsx) (slow)

![SEEPW-T03: interior total head vs time, XSLOPE vs SEEP/W](images/gs2_rdd.png)

The example's *published* answer is a factor-of-safety-vs-time curve from a downstream
SLOPE/W coupling — out of scope here — so the seepage oracle is SEEP/W's own solved
`node.csv` pore-pressure field, and the locked values are XSLOPE's own solved total heads
at four interior stations, checked against the vendor at the initial state, mid-drawdown,
and the near-drained end state:

| station (x, y) | state | XSLOPE h | SEEP/W h | Δ head |
|---|---|---:|---:|---:|
| (20, 5) | IC (full reservoir) | 7.166 m | 7.443 m | −0.28 m |
| (30, 3) | IC (full reservoir) | 4.818 m | 5.241 m | −0.42 m |
| (20, 5) | slow, t = 1.2 d | 6.427 m | 7.038 m | −0.61 m |
| (30, 3) | slow, t = 1.2 d | 4.755 m | 5.208 m | −0.45 m |
| (20, 5) | instantaneous, t = 30 d | 3.714 m | 3.971 m | −0.26 m |
| (35, 2) | instantaneous, t = 30 d | 2.236 m | 2.386 m | −0.15 m |

Across the whole 30-day drawdown the interior seepage field tracks SEEP/W to within about
0.3–0.6 m of head (the figure shows the XSLOPE dissipation curves running just below the
SEEP/W markers at every station and both drawdown rates). The residual is the recurring
**SWCC-mapping caveat** (the van Genuchten fit reproduces the retention curve closely but
shifts the unsaturated-zone drainage timing slightly against SEEP/W's tabulated spline)
plus a convention difference at the upstream face itself: XSLOPE's series-head exit face
drops to pressure head 0 as the reservoir falls, where SEEP/W's zero-flux *review* face
releases the trapped upstream pressures more gradually — so the two diverge most in the
first day near the upstream slope and converge as the dam drains. The locks are on the
robust interior stations at the IC, mid-drawdown, and end state (SWCC timing barely enters
the two near-steady end members), at a 0.03 m regression tolerance on XSLOPE's own values.

**Sources:** GeoStudio SEEP/W example "Rapid Drawdown" (Seequent); the SLOPE/W factor-of-
safety coupling is documented on the [importer verification](#importer-verification) rows
above but is not part of this seepage lock.

<!-- test: file=files/geostudio/gs2_rdd_inst.xlsx, type=tseep_head, target_size=0.7, time=0, max_head_change_frac=0.05, points=20:5:7.166;25:5:6.030;30:3:4.818;35:2:3.216, tolerance=0.03, benchmark=SEEPW-RDD-inst-ic -->
<!-- test: file=files/geostudio/gs2_rdd_inst.xlsx, type=tseep_head, target_size=0.7, time=2592000, max_head_change_frac=0.05, points=20:5:3.714;25:5:3.743;30:3:3.167;35:2:2.236, tolerance=0.03, benchmark=SEEPW-RDD-inst-end -->
<!-- test: file=files/geostudio/gs2_rdd_slow.xlsx, type=tseep_head, target_size=0.7, time=103638, max_head_change_frac=0.05, points=20:5:6.427;25:5:5.857;30:3:4.755;35:2:3.227, tolerance=0.03, benchmark=SEEPW-RDD-slow-mid -->
<!-- test: file=files/geostudio/gs2_rdd_slow.xlsx, type=tseep_head, target_size=0.7, time=2592000, max_head_change_frac=0.05, points=20:5:3.750;25:5:3.778;30:3:3.200;35:2:2.261, tolerance=0.03, benchmark=SEEPW-RDD-slow-end -->

### SEEPW-T04 — Leakage from pond with clay liner {#seepw-t04}

The two-dimensional unconfined transient with an exit face — the census flagship for a
storage-driven **water-table rise**. A clay-lined pond on a hillside (a symmetric
half-model, x = 0 the pond centre-line) is filled to a constant level; water leaks down
through the low-permeability liner into the embankment, and the phreatic surface rises
over 240 days from its initial far-field position toward a new near-steady leaking state
that drains out the downstream seepage face. The census asks for a **linear** mesh here
(the solver's quadratic-exit-face caveat), so the model is meshed with tri3 elements.

The embankment fill is Ksat = 1.157×10⁻⁶ m/s (θ_s = 0.35, θ_r = 0.032 so S_y = 0.318);
the clay liner is ~12× less permeable (Ksat = 9.259×10⁻⁸ m/s, θ_s = 0.45, θ_r = 0.131).
Both vendor retention curves are mapped to van Genuchten by a least-squares fit of their
20-point splines (fill α = 0.657 /m, n = 1.91; liner α = 0.160 /m, n = 1.97). The initial
condition is the pre-fill steady state — the pond series held at head 4 m, below the pond
floor at y = 10 m, so the floor nodes are unsubmerged (inactive exit faces) and only the
far-field water table at el. 4 sets the field (uniform total head 4). For t > 0 the pond
series steps to head 10.5 m (the floor submerges to a Dirichlet head) and the pond leaks —
the same submerged-only reservoir series idiom the rapid-drawdown port ([T03](#seepw-t03))
uses, run in the filling direction.

**Input:** [gs2_pond.xlsx](files/geostudio/gs2_pond.xlsx)

![SEEPW-T04: water-table rise vs time, XSLOPE vs SEEP/W](images/gs2_pond.png)

The example's *published* answer is a water-table-rise-vs-time graph (no closed form), so
the seepage oracle is SEEP/W's own solved `node.csv`, and the locked values are XSLOPE's
own solved total heads at interior stations at the initial state and the near-steady
leaking end state:

| station (x, y) | state | XSLOPE h | SEEP/W h | Δ head |
|---|---|---:|---:|---:|
| (5, 2) | IC (pre-fill) | 4.000 m | 3.999 m | +0.00 m |
| (5, 2) | t = 24 d (filling) | 4.228 m | 4.133 m | +0.10 m |
| (3, 5) | t = 24 d (filling) | 4.357 m | 4.252 m | +0.11 m |
| (5, 2) | t = 240 d (near-steady) | 6.815 m | 6.462 m | +0.35 m |
| (10, 3) | t = 240 d (near-steady) | 6.198 m | 5.899 m | +0.30 m |
| (20, 4) | t = 240 d (near-steady) | 4.442 m | 4.427 m | +0.02 m |

The initial condition is exact (uniform head 4). Through the filling the interior heads
track SEEP/W within about 0.1 m, and the downstream toe — where the field is pinned by the
far water table and the seepage face — matches to 0.02 m at every time. At the fully
developed leaking state the interior water table stands ~0.3–0.35 m higher in XSLOPE than
in SEEP/W: the recurring **SWCC-mapping caveat** (the van Genuchten fit of the liner and
fill retention curves shifts the storage-vs-suction relationship, so the drained-then-
refilling volume and thus the equilibrium table differ slightly) plus the linear-mesh /
lumped-mass discretization. The figure shows the XSLOPE rise curves running just above the
SEEP/W markers at the interior stations and converging on them near the toe. The locks are
on the interior stations at the two near-steady end members (IC and t = 240 d), at a 0.03 m
regression tolerance on XSLOPE's own values.

**Sources:** GeoStudio SEEP/W example "Leakage from Pond with Clay Liner" (Seequent).

<!-- test: file=files/geostudio/gs2_pond.xlsx, type=tseep_head, target_size=0.5, time=0, max_head_change_frac=0.05, points=5:2:4.0000;10:3:4.0000;15:4:4.0000, tolerance=0.03, benchmark=SEEPW-POND-ic -->
<!-- test: file=files/geostudio/gs2_pond.xlsx, type=tseep_head, target_size=0.5, time=2.0736e+07, max_head_change_frac=0.05, points=5:2:6.8152;10:3:6.1982;15:4:5.4564, tolerance=0.03, benchmark=SEEPW-POND-end -->

### SEEPW-T05 — Mineral heap leaching {#seepw-t05}

A one-dimensional column of leach ore under a surface irrigation flux — the test of
XSLOPE's **specified-flux (Neumann) boundary** path and the van Genuchten moisture-capacity
storage in a gravity-drained unsaturated column. The uniform silty-sand ore column (8 m,
Ksat = 5×10⁻⁶ m/s, θ_s = 0.5, θ_r = 0.025 so S_y = 0.475; the vendor retention spline
maps to van Genuchten α = 1.39 /m, n = 1.90) drains freely at its base (specified head 0)
and takes a downward irrigation flux at its top, bound to a `tseep` series. The initial
condition is the low-rate steady state (q = 3×10⁻⁷ m/s); at t = 0 the series steps to the
high rate (q = 3×10⁻⁶ m/s) and the extra water works its way down. (The example's *other*,
non-transient analyses layer the coarse/fine ore; the transient rate case runs on the
single uniform ore, as read from the vendor's own material assignment.)

The low-rate initial condition is a gravity-drained **unit-gradient** profile: away from
the base the steady pressure head is the constant ψ with kr(ψ)·Ksat = q, i.e. kr = 0.06 so
ψ ≈ −0.78 m, relaxing to ψ = 0 at the drained base. That nonlinear unsaturated steady
state is only found by XSLOPE's unsaturated (exit-face) initial-condition solver, not the
linear confined one, so the column's impermeable upper side is declared a potential
seepage-exit face. It never activates (ψ < 0 over the whole run, so it is a no-flow
boundary in fact) but it routes the t = 0 solve through the unsaturated solver, which
returns the unit-gradient profile; without it the linear confined IC would return a dry
hydrostatic column and the front would never advance.

**Input:** [gs2_heap.xlsx](files/geostudio/gs2_heap.xlsx)

![SEEPW-T05: pressure-head profile vs time, XSLOPE vs SEEP/W](images/gs2_heap.png)

The published answer is a graphical volumetric-water-content / flow-rate response (no
closed form), so the seepage oracle is SEEP/W's own solved `node.csv`. XSLOPE reproduces
the low-rate initial and early profiles within ~0.04 m of pressure head; at the high-rate
near-steady end state it reaches a flatter, slightly wetter unit-gradient profile than
SEEP/W (up to ~0.12 m of head at the deep stations), the SWCC-mapping timing caveat — the
van Genuchten kr(ψ) wets the column a little faster than SEEP/W's tabulated conductivity
spline. The figure shows the XSLOPE markers on the SEEP/W low-rate profile at the IC and
early frames and standing off it at the high-rate near-steady. The locks are XSLOPE's own
solved total heads at interior elevations at the IC and end state, at a 0.03 m regression
tolerance.

**Sources:** GeoStudio SEEP/W example "Mineral Heap Leaching" (Seequent).

<!-- test: file=files/geostudio/gs2_heap.xlsx, type=tseep_head, target_size=0.2, time=0, max_head_change_frac=0.05, points=0.1:2:1.2471;0.1:4:3.2540;0.1:6:5.2570, tolerance=0.03, benchmark=SEEPW-HEAP-ic -->
<!-- test: file=files/geostudio/gs2_heap.xlsx, type=tseep_head, target_size=0.2, time=345600, max_head_change_frac=0.05, points=0.1:2:1.8479;0.1:4:3.8634;0.1:6:5.8670, tolerance=0.03, benchmark=SEEPW-HEAP-end -->

### SEEPW-T06 — Infiltration into multi-layered system (blocked) {#seepw-t06}

A laboratory column ponds water on a fourteen-layer soil profile and watches a wetting
front descend (the "Infiltration" analysis), then lets it drain (a "Drainage" child). The
census flagged the drainage leg as out of envelope up front — it is hysteretic, and XSLOPE
carries a single retention curve per material — and planned to port the infiltration leg
only. Two further, harder walls appeared in the infiltration leg itself, and neither is a
builder change:

1. **A non-steady, per-layer initial condition.** The vendor imposes the initial state
   layer by layer through a per-material initial pore-pressure attribute, and it is not a
   hydrostatic or steady field — it is a measured profile with a −50 kPa suction *spike*
   wedged between −5 and −30 kPa layers near the surface. No steady solve returns that
   field. XSLOPE's transient solver *can* take an arbitrary initial field through its
   `h_init` argument, but the corpus lock runner computes the initial condition as a t = 0
   **steady** solve of the boundary configuration; it has no way to inject a per-node
   initial head from a tag, so the forward solve — expressible in code — cannot be locked.

2. **A unit-gradient (free-drainage) base boundary.** The base is a unit-gradient outlet
   (q = kr·Ksat out under gravity). XSLOPE's seepage BC set is specified head, specified
   flux and potential seepage face — there is no unit-gradient boundary, and an exit face
   clamps the base to ψ = 0 rather than letting it drain under a unit gradient.

Either wall alone blocks a faithful lock; together they put the infiltration leg out of
the current envelope. The van Genuchten storage and kr path this example would exercise is
already covered against clean oracles by [SEEPW-T02](#seepw-t02) (infiltration into dry
soil) and [SEEPW-T05](#seepw-t05) (the leach column), so nothing is lost by parking it. The
faithful model is captured in the parked builder `benchmarks/geostudio/build_gs2_mlayer.py`.

**Sources:** GeoStudio SEEP/W example "Infiltration into Multi-Layered System" (Seequent);
field / HYDRUS-1D references (Zettl 2011, Huang 2011).

### SEEPW-T07 — GeoStudio-PEST Multistep Outflow {#seepw-t07}

A multistep-outflow laboratory experiment: a coarse-sand sample (van Genuchten
a = 8.91 /m, n = 10.19, S_y = 0.319) sits on a saturated porous ceramic plate (two
orders of magnitude less permeable than the sand, so it meters the outflow), and the
base **suction** is stepped progressively more negative in five stages over ~61 hours.
It exercises the unsaturated storage term C(ψ) under a stepped specified-pressure-head
boundary.

The base head is a *suction* — a specified **pressure head that is negative at every
stage** (IC −0.073 m, stepping to −0.093 … −0.175 m). Because the base polyline sits at
y = 0, its total head equals its pressure head, and it is carried by a time-varying
**head** (plain-Dirichlet) series: the boundary is held at h(t) at every node of the
polyline at all times, so the negative-pressure Dirichlet is applied faithfully (unlike
the submerged-only *reservoir* series, which would flip an unsubmerged node to a
pressure-head-0 exit face and drop the suction). The IC is the uniform H = −0.073 m
column reached by the t = 0 steady solve, set with the repeated-time step-series idiom.

The published external answer is the lab outflow curve (the scalar the example's PEST
loop fits), not a seepage headline number, so the lock is XSLOPE's own solved total-head
field as a regression guard; the SEEP/W `node.csv` pore-water pressures are read as the
comparison oracle. The high-conductivity sample equilibrates within the sample in
seconds, so at each reporting time the column has drained to the current base suction and
the total head is uniform at that stage value — the hydrostatic profile ψ(y) = h − y then
matches SEEP/W's stepped-suction field to within the read-off precision of the published
`node.csv` (the residual is storage-discretization; the van Genuchten fit reproduces the
retention curve to RMS 1×10⁻⁴, so it is not SWCC mapping).

**Input:** [gs2_mso.xlsx](files/geostudio/gs2_mso.xlsx)

| t (s) | Stage base suction (m) | XSLOPE total head (m) | ψ at y = 0.02 … 0.10 m |
|---|---|---|---|
| 46 000 | −0.093 | −0.0932 | −0.113 … −0.193 |
| 132 000 | −0.134 | −0.1341 | −0.154 … −0.234 |
| 219 600 | −0.175 | −0.1749 | −0.195 … −0.275 |

<!-- test: file=files/geostudio/gs2_mso.xlsx, type=tseep_head, target_size=0.004, time=46000, points=0.003:0.02:-0.0932;0.003:0.06:-0.0932;0.003:0.1:-0.0932, tolerance=0.01, benchmark=SEEPW-T07-t46000 -->
<!-- test: file=files/geostudio/gs2_mso.xlsx, type=tseep_head, target_size=0.004, time=132000, points=0.003:0.02:-0.1341;0.003:0.06:-0.1341;0.003:0.1:-0.1341, tolerance=0.01, benchmark=SEEPW-T07-t132000 -->
<!-- test: file=files/geostudio/gs2_mso.xlsx, type=tseep_head, target_size=0.004, time=219600, points=0.003:0.02:-0.1749;0.003:0.06:-0.1749;0.003:0.1:-0.1749, tolerance=0.01, benchmark=SEEPW-T07-t219600 -->

**Sources:** GeoStudio SEEP/W example "GeoStudio-PEST – Multistep Outflow" (Seequent).
