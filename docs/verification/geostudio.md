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

That is no longer the only route. Seequent publishes the **models** behind the manual, not just the figures, and
XSLOPE can now read them directly — see [GeoStudio Import/Export](../usage/geostudio.md). Importing a problem's
`.gsz` gives its exact geometry and material properties with no transcription step, which removes the ambiguity
that blocks the remaining rows. Better still, a solved `.gsz` carries **SLOPE/W's own trial surfaces and their
factors of safety**, so XSLOPE can be run on the *identical circle* and the two programs compared with no
difference in search to explain away.

The XSLOPE `.xlsx` inputs linked below are our own files and are published here. The GeoStudio files they were
derived from are Seequent's, and are not redistributed — import them yourself from the manual's companion
downloads if you want to reproduce this work.

| § | Problem | Status | XSLOPE file / results vs SLOPE/W |
|---:|---|---|---|
| 2.1 | ACADS Simple Slope | **built** | [xslope_acads_simple.xlsx](../lem/files/xslope_acads_simple.xlsx): Bishop 0.985 vs SLOPE/W 0.963, Slide 0.987; ACADS reference 1.00. — [details](rocscience.md#vp1) |
| 2.2 | ACADS Tension Crack | **built** | [vp002.xlsx](../files/rocscience/vp002.xlsx): Bishop 1.589 / M-P 1.586 vs SLOPE/W 1.664 / 1.660 and Slide 1.596 / 1.592; ACADS reference 1.65-1.70. SLOPE/W sits closer to the ACADS band; the difference is tension-crack water handling and search. — [details](rocscience.md#vp2) |
| 2.3 | ACADS Non-Homogeneous | **built** | [vp003.xlsx](../files/rocscience/vp003.xlsx): Bishop 1.403 / M-P 1.371 vs SLOPE/W 1.414 / 1.382; ACADS 1.39. — [details](rocscience.md#vp3) |
| 2.4 | ACADS Non-Homogeneous + Seismic | **built** | [vp004.xlsx](../files/rocscience/vp004.xlsx): Bishop 1.013 / M-P 0.987 vs SLOPE/W 1.02 / 0.989; ACADS 1.00. — [details](rocscience.md#vp4) |
| 2.5 | ACADS Talbingo Dam – Dry | **built** | [vp005.xlsx](../files/rocscience/vp005.xlsx): 1.955 (all methods, infinite-slope mechanism) vs SLOPE/W 1.951. — [details](rocscience.md#vp5) |
| 2.6 | ACADS Talbingo – Specified Surface | **built** | [vp006.xlsx](../files/rocscience/vp006.xlsx): Bishop 2.206 / M-P 2.299 vs SLOPE/W 2.207 / 2.299 — exact. — [details](rocscience.md#vp6) |
| 2.7 | ACADS Weak Layer | **built** | [xslope_acads_weak_layer.xlsx](../lem/files/xslope_acads_weak_layer.xlsx): Spencer 1.258 / M-P 1.248 vs SLOPE/W Bishop 1.269 / M-P 1.261 — [details](#acads-weak-layer). |
| 2.8 | ACADS Weak Layer – Specified Surface | **built** | [vp008.xlsx](../files/rocscience/vp008.xlsx): M-P 1.260 vs SLOPE/W 1.261 — exact; Janbu(corr) 1.294 vs SLOPE/W force 1.197 (×fo ≈ 1.29). — [details](rocscience.md#vp8) |
| 2.9 | ACADS External Loading | **built** | [vp009.xlsx](../files/rocscience/vp009.xlsx): Spencer 0.724 / Janbu(corr) 0.718 vs SLOPE/W Bishop 0.699 / M-P 0.689 — search-sensitive problem, published spread 0.67-0.81. — [details](rocscience.md#vp9) |
| 2.10 | Lanester Embankment | planned | GS-only problem; figure to be assessed. |
| 2.11 | Arai & Tagyo Homogeneous | **built** | [xslope_arai_tagyo.xlsx](../lem/files/xslope_arai_tagyo.xlsx) vs SLOPE/W Bishop 1.417 / M-P 1.414; A&T 1.451. — [details](rocscience.md#vp14) |
| 2.12 | Arai & Tagyo Pore-Water Pressure | **built** | [vp016.xlsx](../files/rocscience/vp016.xlsx): Bishop 1.112 vs Slide 1.118, A&T 1.138 — SLOPE/W reports 1.190, the outlier of the four sources. — [details](rocscience.md#vp16) |
| 2.13 | Greco Layered Slope | **built** | [vp019.xlsx](../files/rocscience/vp019.xlsx): circular Spencer 1.429 vs SLOPE/W M-P 1.389, Greco 1.40-1.42. — [details](rocscience.md#vp19) |
| 2.14 | Greco Weak Layer | **built** | [vp020.xlsx](../files/rocscience/vp020.xlsx): noncircular Spencer 1.082, circular 1.091 vs SLOPE/W Spencer 1.054, Greco 1.08. — [details](rocscience.md#vp20) |
| 2.15 | Chen & Shao Frictionless Slope | planned | Same problem as Rocscience #25 (SLOPE/W Spencer 1.036, Slide 1.051, theory 1.0) — see that row. |
| 2.16 | Prandtl Bearing Capacity | blocked | Same as Rocscience #26: flat-ground mechanism rejected by the flat-arc guard (feature gap). |
| 2.17 | Chowdhury & Xu (1995), 5 examples | planned | Probabilistic examples (labeled figures in the manual) — candidates for the reliability pipeline. |
| 2.18 | Borges & Cardoso Geosynthetic Emb. #2 | planned | GS-only reinforced embankment. |
| 2.19 | Borges & Cardoso Geosynthetic Emb. #3 | planned | GS-only reinforced embankment. |
| 2.20 | Probabilistic – Syncrude Dyke | planned | Probabilistic case history. |
| 2.21 | Cannon Dam | planned | Probabilistic case history. |
| 2.22 | Cannon Dam #2 | planned | Probabilistic case history. |
| 2.23 | Li & Lumb – Reliability Index | **built** | [vp036.xlsx](../files/rocscience/vp036.xlsx): deterministic Bishop 1.333, β_ln 2.263 vs H&W 1.334 / 2.336; GS reports minimum β 2.04 at FS 1.190 across surfaces. — [details](rocscience.md#vp36) |
| 2.24 | Tandjiria – Geosynthetic Reinforced Emb. | planned | GS-only reinforced embankment. |
| 2.25 | Baker & Leshchinsky – Earth Dam | partial | Same as Rocscience #42: phreatic line through the core is unlabeled in both manuals — needs the B&L (2001) paper. |
| 2.26 | Baker – Planar Homogeneous | planned | Planar-surface comparison (Rocscience #43 analog). |
| 2.27 | Sheahan – Amherst Soil Nails | planned | Nail-wall case history. |
| 2.28 | Sheahan – Clouterre Test Wall | planned | Nail-wall case history. |
| 2.29 | Snailz – Reinforced Slope | planned | SNAILZ nail example (companion to 2.30). |
| 2.30 | Snailz – Geotextile Layers | **built** | [vp050.xlsx](../files/rocscience/vp050.xlsx) (Rocscience #50, same SNAILZ model): Janbu(corr) 1.448 vs SLOPE/W force 1.354 (×fo ≈ 1.44) and SNAILZ 1.46; M-P/Spencer 1.576 vs SLOPE/W M-P 1.606. — [details](rocscience.md#vp50) |
| 2.31 | Zhu – Four Layer Slope | **built** | [vp051.xlsx](../files/rocscience/vp051.xlsx) (Rocscience #51): Bishop 1.278 vs SLOPE/W 1.284; Spencer 1.294 vs 1.299; M-P 1.304 vs 1.310; Lowe 1.296 vs 1.283; Corps 1.404 vs 1.368. — [details](rocscience.md#vp51) |
| 2.32 | Zhu & Lee – Heterogeneous Slope | **built** | [vp052a/b.xlsx](../files/rocscience/vp052a.xlsx) (Rocscience #52): wet deep-family Spencer 1.189 matches Slide exactly; see that row. — [details](rocscience.md#vp52) |
| 2.33 | Priest – Rigid Blocks | planned | Block-mechanism comparison. |
| 2.34 | Yamagami – Stabilizing Piles | **built** | [vp054a/b.xlsx](../files/rocscience/vp054a.xlsx) (Rocscience #54): no-pile Bishop 1.100 vs SLOPE/W 1.102 — exact; with-pile 1.185 vs SLOPE/W 1.223, Slide 1.193, Yamagami 1.20 (pile-force conventions differ program-to-program). — [details](rocscience.md#vp54) |
| 2.35 | Pockoski & Duncan – Tie-Back Wall | partial | P&D series — blocked on the CGPR report (see Rocscience #55-63 note). |
| 2.36 | Pockoski & Duncan – Reinforcement | partial | Same. |
| 2.37 | Pockoski & Duncan – Soil Nails | partial | Same. |
| 2.38 | Loukidis – Seismic Coefficient | **built** | [vp062a/b.xlsx](../files/rocscience/vp062a.xlsx) (Rocscience #62): Spencer 1.001 (both cases) vs SLOPE/W 1.00 — exact. — [details](rocscience.md#vp62) |
| 2.39 | Loukidis – Seismic Coefficient #2 | partial | See Rocscience #63 — outline pinned from the paper, interface anchors still ambiguous. |
| 2.40 | Rapid Drawdown – Walter Bouldin Dam | **built** | Same problem as Slide VP98 ([details](rocscience.md#vp98)): xslope DWW 3-stage 1.046 vs SLOPE/W Bishop 1.016 / Spencer 1.02, DWW 1.04. |
| 2.41 | Rapid Drawdown – USACE Benchmark | **built** | [vp096.xlsx](../files/rocscience/vp096.xlsx) (Rocscience #96): 3-stage Spencer 1.434 / Bishop 1.432 vs published 1.44. — [details](rocscience.md#vp96) |
| 2.42 | Rapid Drawdown – Pumped Storage Dam | **built** | Same problem as Slide VP99 ([details](rocscience.md#vp99)): xslope 1.390 vs SLOPE/W 1.550, DWW 1.56 (~7% low; geometry to be re-pinned from the .gsz). |
| 2.43 | Rapid Drawdown – Pilarcitos Dam | **built** | [vp097.xlsx](../files/rocscience/vp097.xlsx) (Rocscience #97): Spencer 1.044 / Bishop 1.042. — [details](rocscience.md#vp97) |
| 2.44 | Probability – James Bay Case History | planned | Classic reliability case history (El-Ramly et al.) — strong candidate for the reliability pipeline. |
| 2.45 | Eurocode 7 – Cutting in Clay | planned | Partial-factor design check. |
| 2.46 | Eurocode 7 – Earth Dam | planned | Partial-factor design check. |
| 2.47 | Compound Strength vs Anisotropic Function | planned | Anisotropic-strength comparison (feature assessment needed). |

## Problem details

### ACADS weak-layer slope (non-circular) {#acads-weak-layer}

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

