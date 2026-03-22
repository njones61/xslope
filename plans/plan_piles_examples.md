# Plan: Pile Sample Problems

## Motivation

The current pile sample problem (1.5H:1V slope, D=1.0 ft, S=6.0 ft, c=200 psf, phi=10°)
has a shallow failure surface at the pile (7.3 ft depth), which produces a modest LEM pile
force and a large gap between LEM (FS=1.08) and FEM (FS=1.37). We need additional
examples that:

1. Show LEM and FEM converging when the pile is placed where the failure surface is deep
2. Demonstrate multi-layer Ito & Matsui integration
3. Provide independent verification against published results
4. Cover a range of soil types (cohesive, frictional, c-phi)

## Candidate Problems

### 1a. Ausilio — Hassiotis Slope (deep pile, low phi)

**Source**: Ausilio, Conte & Dente (2001), "Stability analysis of slopes reinforced with piles",
Computers and Geotechnics, 28: 591-611. Paper in ref_docs/ref_docs_lim_eq/.

**Slope geometry**:
- H = 13.7 m, slope angle β = 30° from horizontal (1.73H:1V)
- Slope toe at origin, crest at (H/tan30° = 23.7 m, 13.7 m)

**Soil properties** (single layer):
- c = 23.94 kPa, φ = 10°, γ = 19.63 kN/m³

**Without piles**:
- FS = 1.11 (Ausilio, upper bound limit analysis)
- FS = 1.08 (Hassiotis et al., friction circle method)
- FS = 1.12 (Hull & Poulos, Bishop's method)

**With piles** (at xF = 13.7 m from toe, roughly mid-slope):
- Target FS = 1.50
- Required stabilizing force F = 515 kN/m (from limit analysis)
- Height of pile above failure surface h = 12.7 m (DEEP — this is the key)
- Pile length ≈ 2h = 25 m
- Force distribution: linear from ground surface to failure surface (m = 1/3)
- Critical surface with piles is deeper and passes beneath toe

**Why this is ideal**:
- Deep failure surface at pile (h = 12.7 m) vs our current problem (7.3 ft = 2.2 m)
- Should produce a large Ito & Matsui force, making LEM pile contribution dominant
- φ = 10° same as our current sample — direct comparison of shallow vs deep pile
- Published FS values from 3 independent methods for verification
- Published F = 515 kN/m to compare against Ito & Matsui auto-computation

**To implement**: Need to pick pile D and S, compute Ito & Matsui H, run LEM and FEM,
compare with the published F = 515 kN/m and target FS = 1.50.

### 1b. Ausilio — Parametric Slope (mid-range phi)

**Source**: Same paper, Fig. 7 example.

**Slope geometry**:
- Shown in Fig. 7 (need to read dimensions from figure)

**Soil properties** (single layer):
- c = 4.7 kPa, φ = 25°, γ = 20 kN/m³

**Without piles**:
- FS = 1.19

**With piles**:
- Various pile positions from toe to crest studied (Figs. 8-10)
- Multiple improvement ratios (FS/FS0) and m values (0, 1/3, 1/2)
- Optimal position depends on m and target FS

**Why this is useful**:
- Higher φ (25°) exercises the c-phi formula in a different range
- Parametric data on pile position allows validating optimal placement
- Can compare Ito & Matsui force vs the paper's limit analysis force at different positions

### 2. Rocscience Slide2 Tutorial 21

**Source**: https://www.rocscience.com/help/slide2/tutorials/tutorials-overview/pile-resistance-using-rspile

- 3-layer problem: medium sand (phi=40°, gamma=18 kN/m³), dense sand (gamma=20 kN/m³),
  soft clay (c_u=82 kPa, gamma=17 kN/m³)
- Steel pipe pile: D=0.61 m, wall thickness=0.02 m, E=200,000 MPa, length=21 m, spacing=5 m
- FS without piles ≈ 1.4, target FS = 1.5
- Uses RSPile for p-y based pile resistance (not Ito & Matsui)
- **Good for**: Multi-layer exercise; comparing Ito & Matsui auto-computation against
  p-y analysis; realistic multi-material geometry

**Status**: Tutorial available online but soil property details may be incomplete.
Need to reconstruct full geometry from the tutorial.

### 3. Parametric Study from Emerald/ScienceDirect Paper

**Source**: "Factor of safety of pile-stabilised slopes: an algorithm incorporating soil-arching effect",
Geotechnical Research, 8(4), 117.

- c-phi soil: c=30 kPa, phi=20°, gamma=20 kN/m³
- 2H:1V slope
- Pile D=1.0 m, various spacings
- Optimal pile position near mid-slope
- FS improvement up to 1.4x
- **Good for**: Direct comparison with Ito & Matsui (same theory basis);
  parametric study of spacing effects; mid-slope pile placement

**Status**: Behind paywall. Need access to full paper.

### 4. Hassiotis, Chameau & Gunaratne (1997)

**Source**: "Design method for stabilization of slopes with piles", J. Geotech. Geoenviron. Eng., 123(4), 314-323.

- Already cited in our references
- Uses Ito & Matsui theory with friction circle method
- Includes worked example with optimum pile location
- **Good for**: Direct Ito & Matsui verification from the most authoritative secondary reference

**Status**: Behind paywall. Already in our references list.

### 5. Custom Deep-Pile Problem (construct ourselves)

Design a problem specifically to demonstrate LEM-FEM convergence:
- Tall slope (e.g., 15-20 m) with a deep circular failure surface
- Pile placed behind the crest or at mid-slope where failure surface depth is 15+ ft
- c-phi soil so both A1 and A2 contribute significantly
- Same pile modeled in both LEM (Ito & Matsui) and FEM (beam element)
- **Good for**: Directly testing the claim in our documentation that LEM and FEM
  converge when the failure surface is deep at the pile

## Priority

1. **Ausilio — Hassiotis slope (1a)** — highest priority. Deep pile (h=12.7 m) directly
   tests the LEM-FEM convergence claim. Published FS from 3 methods and F=515 kN/m
   for independent verification. Paper is available in ref_docs.
2. **Ausilio — Parametric slope (1b)** — exercises higher phi range and pile position effects
3. **Slide2 Tutorial 21** — useful for multi-layer and comparison with commercial software
4. **Custom deep-pile problem** — fallback if the Ausilio examples don't adequately
   demonstrate LEM-FEM convergence
