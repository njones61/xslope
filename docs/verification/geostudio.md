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

## Building these problems

Most of the entries below were built by transcribing the manual's geometry figures by hand — which is why several
rows are still `planned` or `partial`, blocked on a figure that is ambiguous or a source paper we do not have.

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
seepage solution can regenerate it — a final gap it shares with the Slide2 corpus. The *blocked* / *partial*
rows are each tracked against a named gap: a strength or geometry model XSLOPE does not yet have (§2.16's
flat-ground bearing mechanism behind the flat-arc guard, §2.47's dip-relative anisotropic/compound strength),
or a source document not yet in hand (§2.25's B&L paper, the §2.35–2.37/2.39 Pockoski & Duncan and Loukidis
reports). Everything else is built and verified, or covered by the regression-locked Rocscience build; the
corpus is complete relative to what is independently verifiable.

**Manual edition.** The manual tracked here is the **2025.2 edition**. Its 47
two-dimensional problems (chapter 2, "Verifications – 2D", §2.1–2.47) are carried
unchanged from the October 2022 edition, so the section numbering below is valid against
either. The 2025.2 edition also adds an 11-problem "Verifications – 3D" chapter
(§3.1–3.11) verifying SLOPE3D — ellipsoidal and wedge surfaces, 3D seismic and
piezometric cases, and the Kettleman Hills case history. XSLOPE is a 2D formulation, so
the 3D chapter is out of scope and not tracked here.

| § | Problem | Status | XSLOPE file / results vs SLOPE/W |
|---:|---|---|---|
| 2.1 | ACADS Simple Slope | **built** | [xslope_acads_simple.xlsx](../lem/files/xslope_acads_simple.xlsx): Bishop 0.985 vs SLOPE/W 0.963, Slide 0.987; ACADS reference 1.00. — [details](#gs-2-1) |
| 2.2 | ACADS Tension Crack | **built** | [vp002.xlsx](../files/rocscience/vp002.xlsx): Bishop 1.589 / M-P 1.586 vs SLOPE/W 1.664 / 1.660 and Slide 1.596 / 1.592; ACADS reference 1.65-1.70. SLOPE/W sits closer to the ACADS band; the difference is tension-crack water handling and search. — [details](#gs-2-2) |
| 2.3 | ACADS Non-Homogeneous | **built** | [vp003.xlsx](../files/rocscience/vp003.xlsx): Bishop 1.403 / M-P 1.371 vs SLOPE/W 1.414 / 1.382; ACADS 1.39. — [details](#gs-2-3) |
| 2.4 | ACADS Non-Homogeneous + Seismic | **built** | [vp004.xlsx](../files/rocscience/vp004.xlsx): Bishop 1.013 / M-P 0.987 vs SLOPE/W 1.02 / 0.989; ACADS 1.00. — [details](#gs-2-4) |
| 2.5 | ACADS Talbingo Dam – Dry | **built** | [vp005.xlsx](../files/rocscience/vp005.xlsx): 1.955 (all methods, infinite-slope mechanism) vs SLOPE/W 1.951. — [details](#gs-2-5) |
| 2.6 | ACADS Talbingo – Specified Surface | **built** | [vp006.xlsx](../files/rocscience/vp006.xlsx): Bishop 2.206 / M-P 2.299 vs SLOPE/W 2.207 / 2.299 — exact. — [details](#gs-2-6) |
| 2.7 | ACADS Weak Layer | **built** | [xslope_acads_weak_layer.xlsx](../lem/files/xslope_acads_weak_layer.xlsx): Spencer 1.258 / M-P 1.248 vs SLOPE/W Bishop 1.269 / M-P 1.261 — [details](#acads-weak-layer). |
| 2.8 | ACADS Weak Layer – Specified Surface | **built** | [vp008.xlsx](../files/rocscience/vp008.xlsx): M-P 1.260 vs SLOPE/W 1.261 — exact; Janbu(corr) 1.294 vs SLOPE/W force 1.197 (×fo ≈ 1.29). — [details](#gs-2-8) |
| 2.9 | ACADS External Loading | **built** | [vp009.xlsx](../files/rocscience/vp009.xlsx): Spencer 0.724 / Janbu(corr) 0.718 vs SLOPE/W Bishop 0.699 / M-P 0.689 — search-sensitive problem, published spread 0.67-0.81. — [details](#gs-2-9) |
| 2.10 | Lanester Embankment | *no lock possible* | Same problem as [Rocscience #12](rocscience.md): the printed pore-pressure grid is measured loading-induced pressure, not a flow field, so no seepage solution can reproduce it. |
| 2.11 | Arai & Tagyo Homogeneous | **built** | [xslope_arai_tagyo.xlsx](../lem/files/xslope_arai_tagyo.xlsx) vs SLOPE/W Bishop 1.417 / M-P 1.414; A&T 1.451. — [details](#gs-2-11) |
| 2.12 | Arai & Tagyo Pore-Water Pressure | **built** | [vp016.xlsx](../files/rocscience/vp016.xlsx): Bishop 1.112 vs Slide 1.118, A&T 1.138 — SLOPE/W reports 1.190, the outlier of the four sources. — [details](#gs-2-12) |
| 2.13 | Greco Layered Slope | **built** | [vp019.xlsx](../files/rocscience/vp019.xlsx): circular Spencer 1.429 vs SLOPE/W M-P 1.389, Greco 1.40-1.42. — [details](#gs-2-13) |
| 2.14 | Greco Weak Layer | **built** | [vp020.xlsx](../files/rocscience/vp020.xlsx): noncircular Spencer 1.082, circular 1.091 vs SLOPE/W Spencer 1.054, Greco 1.08. — [details](#gs-2-14) |
| 2.15 | Chen & Shao Frictionless Slope | covered | Same problem as [Rocscience #25](rocscience.md#vp25) — [vp025.xlsx](../files/rocscience/vp025.xlsx), Prandtl mechanism on a 60° weightless slope: Spencer 1.052 vs SLOPE/W 1.036, Slide 1.051, theory 1.0. — [details](#gs-2-15) |
| 2.16 | Prandtl Bearing Capacity | blocked | Same as [Rocscience #26](rocscience.md) — the flat-ground bearing mechanism is rejected by the flat-arc guard (a feature gap, not a data gap). The SLOPE/W model is a public download (see above), and its critical surface is fully specified, so it will import once the guard is relaxed. |
| 2.17 | Chowdhury & Xu (1995), 5 examples | covered | Same problem as [Rocscience #28](rocscience.md#vp28) — Congress St. Cut + embankment on soft clay, 3 of 10 cases built and reliability-tagged. |
| 2.18 | Borges & Cardoso Geosynthetic Emb. #2 | **built** | [gs2_18.xlsx](../files/geostudio/gs2_18.xlsx) — geosynthetic-reinforced embankment on depth-varying soft clay (SFnDepth → XSLOPE `cp`). M-P 1.153 vs SLOPE/W 1.171 (−1.5%) and Borges & Cardoso 1.15 (+0.3%). — [details](#gs-2-18) |
| 2.19 | Borges & Cardoso Geosynthetic Emb. #3 | covered | Same problem as [Rocscience #32](rocscience.md#vp32) / [RS2 #24](rs2.md#rs2-24) — [vp032a](../files/rocscience/vp032a.xlsx) (7 m) / [vp032c](../files/rocscience/vp032c.xlsx) (8.75 m): identical materials, foundation and embankment geometry (verified to <1 cm); the reinforcement-friction difference between vendor sources (39.6° vs 31.0°) is immaterial because the fully-embedded bar develops its full 200 kN/m either way. SLOPE/W's own solves (1.229 / 0.972) bracket the vp032 locks (1.218 / 0.981). |
| 2.20 | Probabilistic – Syncrude Dyke | covered | Same problem as [Rocscience #33](rocscience.md#vp33) — El-Ramly et al. (2003), built deterministically (composite surface through the presheared clay-shale). |
| 2.21 | Cannon Dam | covered | Same problem as [Rocscience #34](rocscience.md#vp34) — Clarence Cannon Dam (Wolff & Harr 1987) on the noncircular surface. The SLOPE/W model is a public download (see above). |
| 2.22 | Cannon Dam #2 | covered | Same problem as [Rocscience #35](rocscience.md#vp35) — Hassan & Wolff (1999), the min-β ≠ min-FS benchmark. |
| 2.23 | Li & Lumb – Reliability Index | **built** | [vp036.xlsx](../files/rocscience/vp036.xlsx): deterministic Bishop 1.333, β_ln 2.263 vs H&W 1.334 / 2.336; GS reports minimum β 2.04 at FS 1.190 across surfaces. — [details](#gs-2-23) |
| 2.24 | Tandjiria – Geosynthetic Reinforced Emb. | **built** | Same problem as [Rocscience #39](rocscience.md#vp39) — [vp039a-d](../files/rocscience/vp039a.xlsx). Also the reinforcement benchmark for the **importer**: on SLOPE/W's own circles the imported geosynthetic reproduces its FS to −0.27% (clay) and −0.64% (sand) — see [Importer verification](#importer-verification). — [details](#gs-2-24) |
| 2.25 | Baker & Leshchinsky – Earth Dam | partial | Same as Rocscience #42: phreatic line through the core is unlabeled in both manuals — needs the B&L (2001) paper. |
| 2.26 | Baker – Planar Homogeneous | **built** | [gs2_26.xlsx](../files/geostudio/gs2_26.xlsx) — Spencer/Janbu 1.352 vs SLOPE/W 1.352 and Baker ≈1.35, essentially exact; the model pins the crest offset at 2.5 m, resolving the Rocscience #43 build's geometry ambiguity. — [details](#gs-2-26) |
| 2.27 | Sheahan – Amherst Soil Nails | covered | Same problem as [Rocscience #47](rocscience.md#vp47) — [vp047.xlsx](../files/rocscience/vp047.xlsx), Sheahan & Ho (2003) Amherst test wall: Janbu 0.899 vs Slide 0.890 / Sheahan 0.887. The SLOPE/W model is a public download (see above). — [details](#gs-2-27) |
| 2.28 | Sheahan – Clouterre Test Wall | covered | Same problem as [Rocscience #48](rocscience.md#vp48) — [vp048.xlsx](../files/rocscience/vp048.xlsx), Clouterre full-scale test wall, 7 nail rows. |
| 2.29 | Snailz – Reinforced Slope | covered | Same problem as [Rocscience #49](rocscience.md#vp49) — [vp049.xlsx](../files/rocscience/vp049.xlsx), SNAILZ soldier-pile tieback wall. |
| 2.30 | Snailz – Geotextile Layers | **built** | [vp050.xlsx](../files/rocscience/vp050.xlsx) (Rocscience #50, same SNAILZ model): Janbu(corr) 1.448 vs SLOPE/W force 1.354 (×fo ≈ 1.44) and SNAILZ 1.46; M-P/Spencer 1.576 vs SLOPE/W M-P 1.606. — [details](#gs-2-30) |
| 2.31 | Zhu – Four Layer Slope | **built** | [vp051.xlsx](../files/rocscience/vp051.xlsx) (Rocscience #51): Bishop 1.278 vs SLOPE/W 1.284; Spencer 1.294 vs 1.299; M-P 1.304 vs 1.310; Lowe 1.296 vs 1.283; Corps 1.404 vs 1.368. Zhu paper source in `ref_docs_lim_eq/`. — [details](#gs-2-31) |
| 2.32 | Zhu & Lee – Heterogeneous Slope | **built** | [vp052a/b.xlsx](../files/rocscience/vp052a.xlsx) (Rocscience #52): wet deep-family Spencer 1.189 matches Slide exactly; see that row. — [details](#gs-2-32) |
| 2.33 | Priest – Rigid Blocks | **built** | [gs2_33.xlsx](../files/geostudio/gs2_33.xlsx) — Janbu 1.049 / M-P 1.049 vs Priest (hand) 1.049 and SLOPE/W 1.049, exact. — [details](#gs-2-33) |
| 2.34 | Yamagami – Stabilizing Piles | **built** | [vp054a/b.xlsx](../files/rocscience/vp054a.xlsx) (Rocscience #54): no-pile Bishop 1.100 vs SLOPE/W 1.102 — exact; with-pile 1.185 vs SLOPE/W 1.223, Slide 1.193, Yamagami 1.20 (pile-force conventions differ program-to-program). — [details](#gs-2-34) |
| 2.35 | Pockoski & Duncan – Tie-Back Wall | partial | P&D series — blocked on the CGPR report (see Rocscience #55-63 note). |
| 2.36 | Pockoski & Duncan – Reinforcement | partial | Same. |
| 2.37 | Pockoski & Duncan – Soil Nails | partial | Same. |
| 2.38 | Loukidis – Seismic Coefficient | **built** | [vp062a/b.xlsx](../files/rocscience/vp062a.xlsx) (Rocscience #62): Spencer 1.001 (both cases) vs SLOPE/W 1.00 — exact. — [details](#gs-2-38) |
| 2.39 | Loukidis – Seismic Coefficient #2 | partial | See Rocscience #63 — outline pinned from the paper, interface anchors still ambiguous. |
| 2.40 | Rapid Drawdown – Walter Bouldin Dam | **built** | Same problem as Slide [VP98](rocscience.md#vp98): xslope DWW 3-stage 1.046 vs SLOPE/W Bishop 1.016 / Spencer 1.02, DWW 1.04. — [details](#gs-2-40) |
| 2.41 | Rapid Drawdown – USACE Benchmark | **built** | [vp096.xlsx](../files/rocscience/vp096.xlsx) (Rocscience #96): 3-stage Spencer 1.434 / Bishop 1.432 vs published 1.44. — [details](#gs-2-41) |
| 2.42 | Rapid Drawdown – Pumped Storage Dam | **built** | Same problem as Slide [VP99](rocscience.md#vp99): xslope 1.390 vs SLOPE/W 1.550, DWW 1.56 (~7% low; geometry to be re-pinned from the .gsz). — [details](#gs-2-42) |
| 2.43 | Rapid Drawdown – Pilarcitos Dam | **built** | [vp097.xlsx](../files/rocscience/vp097.xlsx) (Rocscience #97): Spencer 1.044 / Bishop 1.042. — [details](#gs-2-43) |
| 2.44 | Probability – James Bay Case History | covered | Same problem as [Rocscience #75](rocscience.md#vp75) — the James Bay dyke (El-Ramly et al.). |
| 2.45 | Eurocode 7 – Cutting in Clay | **built** | [gs2_45.xlsx](../files/geostudio/gs2_45.xlsx) — DA3 partial factors baked into the material (c′* = 8.0, φ′* = 23.04°); Spencer 1.173 vs SLOPE/W ODF 1.174 (−0.07%). — [details](#gs-2-45) |
| 2.46 | Eurocode 7 – Earth Dam | **built** | [gs2_46.xlsx](../files/geostudio/gs2_46.xlsx) — DA1-C2 factors + XSLOPE FE seepage; on SLOPE/W's own circle M-P 1.099 vs 1.101 (−0.19%), free-search minimum 1.073 ≈ the Smith textbook 1.07. — [details](#gs-2-46) |
| 2.47 | Compound Strength vs Anisotropic Function | blocked | Needs an orientation-dependent (dip-relative) strength model XSLOPE does not have: AnisotropicFn/CompoundStrength interpolate strength by slice-base angle against a discontinuity dip over angle ranges A/B (the .gsz importer already flags exactly this on import). SLOPE/W: Anisotropic 1.113, Compound 1.118 — solved on the *same* surface these differ by 0.4%, so the gap is the strength model, not the search. The 21-material faulted section has no printed coordinate table (a secondary blocker). |

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

**Input:** [vp002.xlsx](../files/rocscience/vp002.xlsx) · **Rocscience detail:** [VP2](rocscience.md#vp2)

![vp002: inputs and representative solution](images/vp002.png)

| Method | XSLOPE | SLOPE/W | Reference |
|---|---|---|---|
| Bishop | 1.589 | 1.664 | Slide 1.596 |
| Morgenstern-Price | 1.586 | 1.660 | Slide 1.592 |

SLOPE/W sits closer to the ACADS reference band (1.65–1.70) than XSLOPE or Slide2; the difference traces to tension-crack water handling and search.

**Sources:** GeoStudio SLOPE/W Verification Manual §2.2; Giam & Donald (1989).

### 2.3 — ACADS Non-Homogeneous {#gs-2-3}

ACADS benchmark 1(c): a non-homogeneous three-layer slope analyzed for its critical circular surface. It shares the XSLOPE input (vp003.xlsx) with the Rocscience corpus, so the full geometry and material inputs are deferred to the linked Rocscience detail.

**Input:** [vp003.xlsx](../files/rocscience/vp003.xlsx) · **Rocscience detail:** [VP3](rocscience.md#vp3)

![vp003: inputs and representative solution](images/vp003.png)

| Method | XSLOPE | SLOPE/W | Reference |
|---|---|---|---|
| Bishop | 1.403 | 1.414 | 1.39 (ACADS) |
| Morgenstern-Price | 1.371 | 1.382 | 1.39 (ACADS) |

XSLOPE agrees with SLOPE/W to within about 1% for both methods and brackets the ACADS consensus reference of 1.39.

**Sources:** GeoStudio SLOPE/W Verification Manual §2.3; Donald, I.B. & Giam, P. (1989), ACADS.

### 2.4 — ACADS Non-Homogeneous + Seismic {#gs-2-4}

ACADS 1(d): the three-material non-homogeneous slope of §2.3 with a horizontal pseudo-static seismic coefficient of 0.15 added, verifying seismic loading over layered strengths. It shares the XSLOPE input `vp004.xlsx` with the Rocscience corpus, so the full geometry, material properties, and slip surface are documented in the linked Rocscience detail.

**Input:** [vp004.xlsx](../files/rocscience/vp004.xlsx) · **Rocscience detail:** [VP4](rocscience.md#vp4)

![vp004: inputs and representative solution](images/vp004.png)

| Method | XSLOPE | SLOPE/W | Reference |
|---|---|---|---|
| Bishop | 1.013 | 1.02 | 1.00 |
| Morgenstern-Price | 0.987 | 0.989 | 1.00 |

XSLOPE's Bishop (1.013) and Morgenstern-Price (0.987) bracket the ACADS consensus of 1.00 and agree with SLOPE/W to within ~1%.

**Sources:** GeoStudio SLOPE/W Verification Manual §2.4; Donald, I.B. & Giam, P. (1989), ACADS suite.

### 2.5 — ACADS Talbingo Dam – Dry {#gs-2-5}

ACADS benchmark 2(a): the four-zone Talbingo Dam at end of construction (dry), searched for the critical circular surface, whose minimum collapses to a shallow infinite-slope mechanism parallel to the steeper upstream face. It shares the XSLOPE input `vp005.xlsx` with the Rocscience corpus, so the full geometry, zone properties, and per-method breakdown are deferred to the linked Rocscience detail.

**Input:** [vp005.xlsx](../files/rocscience/vp005.xlsx) · **Rocscience detail:** [VP5](rocscience.md#vp5)

![vp005: inputs and representative solution](images/vp005.png)

| Method | XSLOPE | SLOPE/W | Reference |
|---|---|---|---|
| All methods (infinite-slope mechanism) | 1.955 | 1.951 | — |

Because the critical mechanism is the infinite-slope limit, every method converges to the same XSLOPE factor of safety of 1.955, in close agreement with SLOPE/W's 1.951.

**Sources:** GeoStudio SLOPE/W Verification Manual §2.5; Giam & Donald (1989), ACADS benchmark 2(a).

### 2.6 — ACADS Talbingo – Specified Surface {#gs-2-6}

ACADS benchmark 2(b): the Talbingo Dam, a four-material embankment evaluated on a single specified circular slip surface, verifying that XSLOPE and SLOPE/W agree on a fixed surface. It shares the XSLOPE input vp006.xlsx with the Rocscience corpus, so the full geometry and material inputs are deferred to the linked VP6 detail.

**Input:** [vp006.xlsx](../files/rocscience/vp006.xlsx) · **Rocscience detail:** [VP6](rocscience.md#vp6)

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
[Force Equilibrium Methods](../lem/force_eq.md)). This benchmark also appears on the
[Verification](../verification/lem.md) page.

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

**Input:** [vp008.xlsx](../files/rocscience/vp008.xlsx) · **Rocscience detail:** [VP8](rocscience.md#vp8)

![vp008: inputs and representative solution](images/vp008.png)

| Method | XSLOPE | SLOPE/W | Reference |
|---|---|---|---|
| Morgenstern-Price | 1.260 | 1.261 | — |
| Janbu (corrected) | 1.294 | 1.197 (force, ×fo ≈ 1.29) | — |

XSLOPE's Morgenstern-Price matches SLOPE/W's M-P essentially exactly (1.260 vs 1.261); XSLOPE's corrected Janbu (1.294) reproduces SLOPE/W's uncorrected force solution (1.197) once the Janbu fo correction (×fo ≈ 1.29) is applied.

**Sources:** GeoStudio SLOPE/W Verification Manual §2.8; Donald, I.B. & Giam, P. (1989), *Soil slope stability programs review*, ACADS, Melbourne.

### 2.9 — ACADS External Loading {#gs-2-9}

The ACADS benchmark (Slide #9): a two-material slope with a weak layer, a piezometric water table, and two surcharge strips, solved with a non-circular search. It shares the XSLOPE input vp009.xlsx with the Rocscience corpus, so the full geometry and inputs are deferred to the linked Rocscience detail.

**Input:** [vp009.xlsx](../files/rocscience/vp009.xlsx) · **Rocscience detail:** [VP9](rocscience.md#vp9)

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

**Input:** [vp016.xlsx](../files/rocscience/vp016.xlsx) · **Rocscience detail:** [VP16](rocscience.md#vp16)

![vp016: inputs and representative solution](images/vp016.png)

| Method | XSLOPE | SLOPE/W | Reference |
|---|---|---|---|
| Bishop | 1.112 | 1.190 | Slide 1.118; A&T 1.138 |

XSLOPE's Bishop 1.112 agrees closely with Slide 1.118 and the Arai & Tagyo reference 1.138; SLOPE/W's 1.190 is the outlier of the four sources.

**Sources:** GeoStudio SLOPE/W Verification Manual §2.12; Arai & Tagyo (1985).

### 2.13 — Greco Layered Slope {#gs-2-13}

A four-layer slope with no water (Greco 1996, example 4 / Yamagami & Ueta 1988), verifying the critical-surface search against a shallow non-circular benchmark optimum. It shares the vp019.xlsx XSLOPE input with the Rocscience corpus, so full geometry and material inputs are deferred to the linked Rocscience detail.

**Input:** [vp019.xlsx](../files/rocscience/vp019.xlsx) · **Rocscience detail:** [VP19](rocscience.md#vp19)

![vp019: inputs and representative solution](images/vp019.png)

| Method | XSLOPE | SLOPE/W | Reference |
|---|---|---|---|
| Spencer / M-P | 1.429 | 1.389 | Greco 1.40–1.42 |

XSLOPE's circular Spencer (1.429) agrees with SLOPE/W's Morgenstern–Price solution (1.389) and brackets the Greco reference range (1.40–1.42).

**Sources:** GeoStudio SLOPE/W Verification Manual §2.13; Greco (1996), Yamagami & Ueta (1988).

### 2.14 — Greco Weak Layer {#gs-2-14}

A four-layer slope with a 0.5 m weak seam running along the inclined model base beneath a water table, verifying the search for a noncircular slip surface through a weak layer. It shares the XSLOPE input vp020.xlsx with the Rocscience corpus, so full geometry and inputs are deferred to the linked Rocscience detail.

**Input:** [vp020.xlsx](../files/rocscience/vp020.xlsx) · **Rocscience detail:** [VP20](rocscience.md#vp20)

![vp020: inputs and representative solution](images/vp020.png)

| Method | XSLOPE | SLOPE/W | Reference |
|---|---|---|---|
| Spencer (noncircular) | 1.082 | 1.054 | Greco 1.08 |
| Spencer (circular) | 1.091 | — | — |

XSLOPE's noncircular Spencer solution (1.082) matches the published Greco value (1.08) and runs slightly above SLOPE/W's Spencer result (1.054); the circular Spencer surface gives a marginally higher 1.091.

**Sources:** GeoStudio SLOPE/W Verification Manual §2.14; Greco (1996) ex. 5; Chen & Shao (1988).

### 2.15 — Chen & Shao Frictionless Slope {#gs-2-15}

The classical Prandtl bearing mechanism on a weightless, frictionless 60° slope under a critical strip load, verifying limit-equilibrium against the theoretical FS = 1.0. This is the same problem as the Rocscience corpus, sharing XSLOPE input [vp025.xlsx](../files/rocscience/vp025.xlsx); see the linked Rocscience detail for the full geometry, materials, and analytically constructed slip surface.

**Input:** [vp025.xlsx](../files/rocscience/vp025.xlsx) · **Rocscience detail:** [VP25](rocscience.md#vp25)

![vp025: inputs and representative solution](images/vp025.png)

| Method | XSLOPE | SLOPE/W | Reference |
|---|---|---|---|
| Spencer | 1.052 | 1.036 | Slide 1.051; theory 1.0 |

XSLOPE's Spencer FS of 1.052 matches Slide (1.051) and sits just above SLOPE/W (1.036), all close to the theoretical value of 1.0 for this frictionless mechanism.

**Sources:** GeoStudio SLOPE/W Verification Manual §2.15; Chen & Shao (1988).

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

**Input:** [gs2_18.xlsx](../files/geostudio/gs2_18.xlsx)

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

<!-- test: file=../files/geostudio/gs2_18.xlsx, type=single_circle, num_slices=60, fs_bishop=1.155, fs_spencer=1.148, fs_mprice=1.153, fs_janbu=1.330, benchmark=GS-2.18 -->

### 2.23 — Li & Lumb – Reliability Index {#gs-2-23}

The Li & Lumb (1987) / Hassan & Wolff (1999) reliability benchmark: a homogeneous slope with an r_u pore-pressure ratio for which both the deterministic Bishop factor of safety and the lognormal reliability index β are computed from the variability in c′, φ′, and γ. It shares the vp036.xlsx input with the Rocscience corpus, so the full geometry and material statistics are deferred to the linked Rocscience detail.

**Input:** [vp036.xlsx](../files/rocscience/vp036.xlsx) · **Rocscience detail:** [VP36](rocscience.md#vp36)

![vp036: inputs and representative solution](images/vp036.png)

| Method | XSLOPE | SLOPE/W | Reference |
|---|---|---|---|
| Bishop (deterministic FS) | 1.333 | 1.190 | 1.334 |
| Reliability index β_ln | 2.263 | 2.04 | 2.336 |

XSLOPE's deterministic Bishop FS (1.333) and its β on that surface (2.263) match the Hassan & Wolff reference (1.334 / 2.336); SLOPE/W instead searches for the minimum reliability index across surfaces, reporting β 2.04 at FS 1.190, so its β and FS are not evaluated on the deterministic surface.

**Sources:** GeoStudio SLOPE/W Verification Manual §2.23; Li & Lumb (1987), Hassan & Wolff (1999).

### 2.24 — Tandjiria – Geosynthetic Reinforced Emb. {#gs-2-24}

Tandjiria (2002)'s required-reinforcement half-embankment on soft clay, evaluated as a clay fill and as a sand fill with a geosynthetic at the embankment base. It shares the XSLOPE input (vp039a) with the Rocscience corpus, so the full geometry and inputs are deferred to the linked Rocscience detail; here the section's distinct role is as the geosynthetic-reinforcement benchmark for the SLOPE/W importer.

**Input:** [vp039a.xlsx](../files/rocscience/vp039a.xlsx) · **Rocscience detail:** [VP39](rocscience.md#vp39)

![vp039a: inputs and representative solution](images/vp039a.png)

| Method | XSLOPE | SLOPE/W | Reference |
|---|---|---|---|
| Clay fill, reinforced (imported, on SLOPE/W's own circles) | −0.27% | baseline | — |
| Sand fill, reinforced (imported, on SLOPE/W's own circles) | −0.64% | baseline | — |

Run on SLOPE/W's own critical circles, the imported geosynthetic reproduces SLOPE/W's reinforced factor of safety to within −0.27% (clay fill) and −0.64% (sand fill), isolating the reinforcement handling with no search difference to explain away.

**Sources:** GeoStudio SLOPE/W Verification Manual §2.24; Tandjiria (2002).

### 2.26 — Baker – Planar Homogeneous {#gs-2-26}

Baker (2001)'s planar-slip-surface benchmark: a homogeneous, dry c′-φ′ slope (H = 10 m,
face angle 76.0°, c′ = 30, φ′ = 30°, γ = 20) evaluated on planes through the toe, with
FS plotted against the daylight point's normalized position X = x/H on the backslope.
The critical plane sits at X = 0.85. The SLOPE/W model pins the exact geometry (crest
offset 2.5 m) — which also resolves the Rocscience corpus's VP43 analog, whose build
inferred 3 m and read FS ≈ 1.43 against the ≈ 1.35 references.

**Input:** [gs2_26.xlsx](../files/geostudio/gs2_26.xlsx)

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

<!-- test: file=../files/geostudio/gs2_26.xlsx, type=single_noncirc, num_slices=50, fs_spencer=1.352, fs_janbu=1.352, benchmark=GS-2.26 -->

### 2.27 — Sheahan – Amherst Soil Nails {#gs-2-27}

The Amherst test wall — a 6 m vertical cut in undrained clay retained by two rows of soil nails and a shotcrete facing, which failed in the field test — evaluated on planar surfaces through the toe. It reuses the same XSLOPE input as the Rocscience corpus, so the full geometry, nail capacities, and loading are deferred to the linked Rocscience detail.

**Input:** [vp047.xlsx](../files/rocscience/vp047.xlsx) · **Rocscience detail:** [VP47](rocscience.md#vp47)

![vp047: inputs and representative solution](images/vp047.png)

| Method | XSLOPE | SLOPE/W | Reference |
|---|---|---|---|
| Janbu | 0.899 | — | Slide 0.890; Sheahan 0.887 |

XSLOPE's Janbu FS of 0.899 agrees closely with Slide's 0.890 and Sheahan's trial-wedge 0.887; the SLOPE/W model itself is a public download rather than a tabulated FS here.

**Sources:** GeoStudio SLOPE/W Verification Manual §2.27; Sheahan & Ho (2003).

### 2.30 — Snailz – Geotextile Layers {#gs-2-30}

A reinforced-slope problem from the Caltrans SNAILZ reference manual — a multi-row nail/geotextile-reinforced wall evaluated on a predefined deep wedge surface — verifying XSLOPE's reinforcement handling against SLOPE/W and the published SNAILZ solution. It shares the `vp050.xlsx` input with the Rocscience corpus, so full geometry and reinforcement inputs are deferred to the linked VP50 detail.

**Input:** [vp050.xlsx](../files/rocscience/vp050.xlsx) · **Rocscience detail:** [VP50](rocscience.md#vp50)

![vp050: inputs and representative solution](images/vp050.png)

| Method | XSLOPE | SLOPE/W | Reference |
|---|---|---|---|
| Janbu (corrected) | 1.448 | 1.354 (force, ×fo ≈ 1.44) | SNAILZ 1.46 |
| Morgenstern-Price / Spencer | 1.576 | 1.606 (M-P) | — |

XSLOPE's Janbu (corrected) 1.448 matches SLOPE/W's uncorrected force solution 1.354 once the fo correction (×fo ≈ 1.44) is applied and agrees with the SNAILZ 1.46 reference; the Morgenstern-Price/Spencer pair agrees at 1.576 vs SLOPE/W's 1.606.

**Sources:** GeoStudio SLOPE/W Verification Manual §2.30; Caltrans SNAILZ reference manual.

### 2.31 — Zhu – Four Layer Slope {#gs-2-31}

Zhu, Lee & Jiang's four-layer slope with a water table, a 5 m dry tension crack, and seismic loading (k=0.1), analyzed on a specified circle. It shares the vp051.xlsx input with the Rocscience corpus, so the full geometry, material properties, and phreatic-line calibration are deferred to the linked Rocscience detail.

**Input:** [vp051.xlsx](../files/rocscience/vp051.xlsx) · **Rocscience detail:** [VP51](rocscience.md#vp51)

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

**Input:** [vp052a.xlsx](../files/rocscience/vp052a.xlsx) · **Rocscience detail:** [VP52](rocscience.md#vp52)

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

**Input:** [gs2_33.xlsx](../files/geostudio/gs2_33.xlsx)

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

<!-- test: file=../files/geostudio/gs2_33.xlsx, type=single_noncirc, num_slices=50, fs_janbu=1.049, fs_mprice=1.049, fs_spencer=1.049, benchmark=GS-2.33 -->

### 2.34 — Yamagami – Stabilizing Piles {#gs-2-34}

A homogeneous slope stabilized by a row of micro-piles, verifying pile-reinforced limit equilibrium against the unreinforced baseline. It shares the XSLOPE input with the Rocscience corpus (VP54), so full geometry and material inputs are deferred to the linked Rocscience detail.

**Input:** [vp054a.xlsx](../files/rocscience/vp054a.xlsx) · **Rocscience detail:** [VP54](rocscience.md#vp54)

![vp054a: inputs and representative solution](images/vp054a.png)

| Method | XSLOPE | SLOPE/W | Reference |
|---|---|---|---|
| Bishop (no pile) | 1.100 | 1.102 | — |
| Bishop (with piles) | 1.185 | 1.223 | Slide 1.193; Yamagami 1.20 |

The unreinforced case matches SLOPE/W essentially exactly (1.100 vs 1.102), while the reinforced FS spreads because pile-force conventions differ program-to-program.

**Sources:** GeoStudio SLOPE/W Verification Manual §2.34; Yamagami.

### 2.38 — Loukidis – Seismic Coefficient {#gs-2-38}

A homogeneous slope loaded pseudo-statically at its critical seismic coefficient (dry and ru = 0.5 cases), where the factor of safety is driven to ≈ 1. It shares the XSLOPE vp062a/b inputs with the Rocscience corpus, so the full geometry and inputs are deferred to the linked Rocscience detail.

**Input:** [vp062a.xlsx](../files/rocscience/vp062a.xlsx) · **Rocscience detail:** [VP62](rocscience.md#vp62)

![vp062a: inputs and representative solution](images/vp062a.png)

| Method | XSLOPE | SLOPE/W | Reference |
|---|---|---|---|
| Spencer | 1.001 | 1.00 | — |

XSLOPE's Spencer factor of safety is 1.001 in both cases, matching SLOPE/W's 1.00 exactly.

**Sources:** GeoStudio SLOPE/W Verification Manual §2.38; Loukidis, Bandini & Salgado (2003).

### 2.40 — Rapid Drawdown – Walter Bouldin Dam {#gs-2-40}

The Walter Bouldin Dam is a rolled earthfill dam that failed during a rapid drawdown in 1975; this case verifies XSLOPE's Duncan-Wright-Wong three-stage rapid drawdown procedure. It shares the `vp098.xlsx` input with the Rocscience corpus, so the full geometry, zones, and undrained envelopes are deferred to the linked Rocscience detail.

**Input:** [vp098.xlsx](../files/rocscience/vp098.xlsx) · **Rocscience detail:** [VP98](rocscience.md#vp98)

![vp098: inputs and representative solution](images/vp098.png)

| Method | XSLOPE | SLOPE/W | Reference |
|---|---|---|---|
| DWW 3-stage (Spencer) | 1.046 | 1.02 | 1.04 |
| DWW 3-stage (Bishop) | — | 1.016 | — |

XSLOPE's DWW three-stage FS (1.046, Spencer circular search) agrees with SLOPE/W's Spencer (1.02) and the published DWW value (1.04) to within about 1%, and brackets SLOPE/W's Bishop result (1.016).

**Sources:** GeoStudio SLOPE/W Verification Manual §2.40; Duncan, Wright & Wong (1990).

### 2.41 — Rapid Drawdown – USACE Benchmark {#gs-2-41}

A homogeneous embankment dam rapid-drawdown benchmark from USACE EM 1110-2-1902 (2003) Appendix G, evaluated with the Duncan-Wright-Wong 3-stage procedure on the specified circle. It shares the XSLOPE input vp096.xlsx with the Rocscience corpus, so full geometry and material inputs are deferred to the linked Rocscience detail.

**Input:** [vp096.xlsx](../files/rocscience/vp096.xlsx) · **Rocscience detail:** [VP96](rocscience.md#vp96)

![vp096: inputs and representative solution](images/vp096.png)

| Method | XSLOPE | SLOPE/W | Reference |
|---|---|---|---|
| Bishop | 1.432 | — | 1.44 |
| Spencer | 1.434 | — | 1.44 |

XSLOPE's 3-stage Spencer (1.434) and Bishop (1.432) both agree with the published USACE benchmark of 1.44 to within about half a percent.

**Sources:** GeoStudio SLOPE/W Verification Manual §2.41; USACE EM 1110-2-1902 (2003) Appendix G (Duncan, Wright & Wong 3-stage procedure).

### 2.42 — Rapid Drawdown – Pumped Storage Dam {#gs-2-42}

Rapid drawdown (285 → 120 ft) of a hypothetical pumped-storage dam — silty-clay core and free-draining rockfill shells — analyzed with the Duncan-Wright-Wong 3-stage procedure. It shares its XSLOPE input with the Rocscience corpus (VP99), so the full geometry and material inputs are deferred to the linked Rocscience detail.

**Input:** [vp099.xlsx](../files/rocscience/vp099.xlsx) · **Rocscience detail:** [VP99](rocscience.md#vp99)

![vp099: inputs and representative solution](images/vp099.png)

| Method | XSLOPE | SLOPE/W | Reference |
|---|---|---|---|
| Spencer (DWW 3-stage) | 1.390 | 1.550 | 1.56 |

XSLOPE runs ~7% low against both SLOPE/W and the DWW reference; the discrepancy is attributed to the axis-calibrated core geometry, and the model is to be re-pinned from the vendor .gsz when available.

**Sources:** GeoStudio SLOPE/W Verification Manual §2.42; Duncan, Wright & Wong (1990).

### 2.43 — Rapid Drawdown – Pilarcitos Dam {#gs-2-43}

A three-stage rapid-drawdown analysis of Pilarcitos Dam — a homogeneous earthfill embankment (Duncan, Wright & Wong 1990) drawn down from 72 to 37 ft, a documented drawdown-failure case. It verifies XSLOPE's staged drawdown procedure and shares its XSLOPE input with the Rocscience corpus, so the full geometry and strength envelopes are deferred to the linked Rocscience detail (VP97).

**Input:** [vp097.xlsx](../files/rocscience/vp097.xlsx) · **Rocscience detail:** [VP97](rocscience.md#vp97)

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

**Input:** [gs2_45.xlsx](../files/geostudio/gs2_45.xlsx)

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

<!-- test: file=../files/geostudio/gs2_45.xlsx, type=circular_search, num_slices=40, fs_spencer=1.173, fs_bishop=1.172, benchmark=GS-2.45 -->

### 2.46 — Eurocode 7: Earth Dam {#gs-2-46}

A homogeneous clay dam (characteristic c′ = 12, φ′ = 20°, γ = 19.2) with an upstream
reservoir at 5.1 m total head, analysed as coupled steady-state seepage + slope
stability under Eurocode 7 Design Approach 1, Combination 2 (after *Smith's Elements of
Soil Mechanics*, ex. 5.12). Pore pressures come from an XSLOPE finite-element seepage
solve (u = 'seep' with mesh/solution sidecars) whose phreatic surface reproduces
SEEP/W's to ~0.1–0.2 m; the DA1-C2 material factors (set M2) are baked in — c′* = 9.6
kPa, φ′* = 16.23° — so an ordinary analysis yields the Overdesign Factor.

**Input:** [gs2_46.xlsx](../files/geostudio/gs2_46.xlsx) (+ `gs2_46_mesh.json`,
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

<!-- test: file=../files/geostudio/gs2_46.xlsx, type=circular_search, num_slices=40, fs_mprice=1.073, fs_bishop=1.074, fs_spencer=1.074, benchmark=GS-2.46 -->
