# Morgenstern–Price Method — Implementation Plan / Discussion

Status: **DISCUSSION DRAFT** — nothing implemented yet. This document is the
working space for deciding *how* to add Morgenstern–Price (M-P) to xslope before
any code is written.

### Decisions so far (2026-06-23)

- **`f(x)` library at v1:** constant (`=1`, Spencer regression), half-sine, and
  clipped/trapezoidal sine. No user-defined table for v1.
- **`f(x)` selection:** function argument only (e.g. `f_type='half_sine'`); no
  Excel-template cell for now.
- **Solver structure:** flesh out *both* the two-curve intersection and the direct
  2-D Newton (see §3) before committing — they share a building block, so we build
  that first and can switch outer solvers.
- **`F_m` (moment) formulation:** open — work it out together (see §4).

## 1. What M-P is, and how it relates to what we already have

Morgenstern–Price (1965) is a complete-equilibrium method: it satisfies force
equilibrium (horizontal + vertical) *and* moment equilibrium on every slice. In
that sense it is in the same family as Spencer's method, which xslope already
implements.

The single difference between Spencer and M-P is the **interslice force
assumption**:

| Method | Interslice force inclination | Unknowns |
|--------|------------------------------|----------|
| Spencer | constant: `X_i / E_i = tan θ` (same θ everywhere) | `F`, `θ` |
| Morgenstern–Price | varies along the surface: `X_i / E_i = λ · f(x_i)` | `F`, `λ` |

where:

- `E_i` = horizontal interslice force, `X_i` = vertical (shear) interslice force,
- `f(x)` = a prescribed, dimensionless **interslice force function** of position,
- `λ` = a scalar that scales `f(x)`.

**Spencer is exactly the special case `f(x) = 1`** (then `tan θ = λ` is constant).
This is the single most important fact for both implementation and validation:
if we implement M-P correctly and run it with `f(x) = 1`, it *must* reproduce
Spencer's FS and θ to numerical tolerance. That is our primary regression test.

Common choices for `f(x)`:

- **Constant** `f(x) = 1` → Spencer (sanity/regression case).
- **Half-sine** `f(x) = sin(π · (x − x_L)/(x_R − x_L))` → the textbook default; forces
  the interslice shear to zero at both ends of the surface, which is physically
  reasonable and is GeoStudio's default.
- **Clipped/trapezoidal sine**, **user-defined table**, etc. (later, optional).

In practice the computed FS is famously *insensitive* to the choice of `f(x)`
(typically <1% spread), which is the classic argument for M-P being "rigorous yet
robust." λ, however, does depend on `f(x)`.

## 2. Existing infrastructure we can reuse

Two functions in `xslope/solve.py` matter:

1. **`force_equilibrium(slice_df, theta_list, fs_guess, ...)`** (line ~586).
   Marches slices left→right, solving a 2×2 system per slice for the base normal
   `N_i` and the right-side interslice resultant `Z_{i+1}`, given a per-boundary
   inclination array `theta_list` (length n+1). It then root-finds `F` so the
   right-end interslice force `Z[n] = 0` (global force closure). It already stores
   `n_eff` and `z`, and has an admissibility guard (rejects surfaces with >50%
   base-tension or interslice-tension).
   **This is, essentially, the force-equilibrium half of M-P already written.**
   Corps of Engineers and Lowe–Karafiath are thin wrappers that build `theta_list`
   from geometry and delegate to it — they satisfy force equilibrium only.

2. **`spencer(slice_df, ...)`** (line ~890). Uses Steve Wright's UTEXAS lumped-`Q`
   formulation: all side forces on a slice collapse to one resultant `Q_i` acting
   at the *constant* angle θ, and the two governing equations are
   `R1 = Σ Q = 0` (force) and `R2 = Σ Q (x_b sin θ − y_Q cos θ) = 0` (moment),
   solved by a hand-coded Newton iteration with analytic Jacobian, plus
   scipy/bishop fallbacks, a tension guard, and line-of-thrust output.

## 3. The key design decision

### Why we can't just generalize Spencer's `Q` solver

Wright's lumped-`Q` derivation depends on θ being **constant**. With constant θ,
`Σ Q = 0` is equivalent to *both* horizontal (`Σ Q cos θ = 0`) and vertical
(`Σ Q sin θ = 0`) force balance at once — one scalar equation covers both
components. The instant θ varies per slice (`θ_i = arctan(λ f_i)`), that collapse
fails: horizontal and vertical balance become two *separate* equations, so we'd
have 3 global equations (2 force + 1 moment) for only 2 unknowns (`F`, λ) —
over-determined. The reason real M-P works is that it does **not** lump globally;
it carries the interslice force as a per-slice unknown through a recurrence, so
equilibrium is enforced slice-by-slice. So generalizing `spencer()`'s `Q`
machinery is the wrong path.

### The shared building block: a per-slice march evaluator

Both solver structures below sit on top of one routine we build first. Call it

```
_mp_march(slice_df, lam, f_vals, F)  ->  N[ ], Z[ ], force_res, moment_res
```

Given a trial scale `lam` and trial factor of safety `F`:

1. Form per-boundary inclinations `θ_j = arctan(lam · f(x_j))` for `j = 0..n`
   (so `tan θ = λ f`, exactly the M-P interslice assumption; `f` evaluated at each
   slice boundary, normalized over `[x_L, x_R]`).
2. March slices left→right solving the same per-slice 2×2 system that
   `force_equilibrium` already uses, producing the base normal `N_i` and the
   right-side interslice resultant `Z_{i+1}` at this fixed `F`. (We do **not**
   root-find `F` inside the march — we evaluate at the given `F`.)
3. Return two residuals:
   - **`force_res = Z[n]`** — the leftover interslice force at the right end. Zero ⇒
     global force equilibrium. (This is the quantity `force_equilibrium` already
     root-finds on.)
   - **`moment_res`** — the overall moment imbalance about a reference point,
     computed from the `N_i`, base shears `S_i`, weights and external loads (see
     §4). Zero ⇒ global moment equilibrium.

Crucially, the *only* thing that distinguishes M-P from the existing
force-equilibrium methods (Corps, L-K) is that those prescribe `theta_list` from
geometry and ignore `moment_res`; M-P chooses `(F, λ)` to drive **both** residuals
to zero. So `_mp_march` is mostly a refactor of the existing
`force_equilibrium` inner loop plus a moment accumulator.

`F_f(λ)` ("force FS at this λ") = the `F` root of `force_res(F; λ) = 0`.
`F_m(λ)` ("moment FS at this λ") = the `F` root of `moment_res(F; λ) = 0`.

### Approach A — two-curve intersection (Fredlund–Krahn / GLE style)

Outer loop over `λ`, two inner 1-D `F` root finds per `λ`:

```
g(λ) = F_f(λ) − F_m(λ)        # find λ* where g = 0
```

`F_f(λ)` reuses the existing `force_equilibrium` root find essentially verbatim;
`F_m(λ)` is the same march root-found on `moment_res` instead. Root-find `λ` on
`g`. The two FS-vs-λ curves classically cross once and nearly linearly, so a few
secant/brentq steps converge. The crossing value is the M-P factor of safety;
`λ*` fixes the interslice force distribution.

- **Pros:** transparent; reproduces the canonical `F` vs `λ` plot for free (great
  for docs/debugging); reuses the most-tested engine; degrades gracefully (if no
  crossing exists we can report the curves).
- **Cons:** nested root finds → more solver calls; slowest of the two.

### Approach B — direct 2-D Newton on `(F, λ)`

Solve the coupled system
```
force_res(F, λ) = 0
moment_res(F, λ) = 0
```
simultaneously with a 2-D Newton (numeric Jacobian to start; analytic later if
warranted), mirroring how `spencer()` solves its `(F, θ)` system.

- **Pros:** fewest evaluations once it's converging; same spirit as the existing
  Spencer solver; natural home for Bishop-FS / λ=0 seeding and a tension guard.
- **Cons:** needs good initial guesses and safeguards (the `spencer()` solver
  carries ~150 lines of line-search / fallback / tension-rejection machinery for
  exactly this reason); harder to debug a 2-D divergence than a 1-D curve crossing.

### Recommendation

Build `_mp_march` first, then **Approach A** for the reference implementation and
validation (it gives the `F`-vs-`λ` curve, which we need anyway to prove
correctness and to make the docs figure). Keep **Approach B** as a fast path that
must agree with A on every benchmark. Seed both from Bishop FS (like Spencer) and
`λ = 0` (pure force-equilibrium starting point).

## 4. Moment equilibrium `F_m` — working it out

This is the one genuinely new piece of math (force equilibrium is already done in
`force_equilibrium`). Two candidate formulations:

### Option (a) — overall moment about the origin  ← proposed

**Sum moments about the coordinate origin `(0,0)`, exactly as Spencer does** (see
Spencer eq (16): `Σ Q(x_b sin θ − y_Q cos θ) = 0`). The origin is a fixed point
that works identically for circular and non-circular surfaces — there is **no**
reliance on a circle center, so this is uniform across all surface types. (Using
the circle center would collapse the base-normal arms nicely for circles but has
no analogue for non-circular surfaces; the origin avoids that special-casing
entirely.)

**Key fact that makes this clean:** the interslice forces are *internal*. The
force on the boundary between slice `i` and `i+1` acts on both adjacent slices with
equal magnitude, the same line of action, and opposite sign. So when we sum
moments over the *whole* sliding mass about the origin, **every interslice
contribution cancels exactly** — regardless of `λ` and regardless of where on the
boundary the force acts (the thrust line).

Consequences:
- The overall moment equation contains only the base normals `N_i`, base shears
  `S_i`, weights `W_i`, and the external loads xslope already tracks (distributed
  load `D`, seismic `kW`, tension-crack water `V`, reinforcement `R`, piles
  `H, θ_p`). We do **not** need the interslice force positions to get `F_m`.
- `λ` still enters `F_m` — but only *indirectly*, through `N_i`: the per-slice
  vertical balance contains the net interslice shear `X_{i+1} − X_i`
  (`= λ(f_{i+1}E_{i+1} − f_i E_i)`), so changing `λ` redistributes the base
  normals, which moves the moment sum. This indirect coupling is exactly why
  `F_m(λ)` and `F_f(λ)` are different curves that cross at the true solution.

Sketch (moments about the origin, CCW positive, matching Spencer's sign
conventions; `S_i` is the mobilized base shear
`S_i = [c·dl_i + (N_i − u_i dl_i)tanφ_i]/F`):

```
Σ M_origin  =  Σ [ W_i, N_i, S_i and external-load moments about (0,0) ]  = 0
            →  solve for F_m
```

Because the base shear `S_i` carries the `1/F` factor (mobilized strength), the
moment balance is an equation in `F` for each trial `λ` — its root is `F_m(λ)`.
Build it from the `_mp_march` output: take the per-slice `N_i` (which already
embed the per-slice interslice forces) and `S_i`, `W_i`, plus the external loads,
and form their moment about the origin. We do **not** reuse Spencer's lumped-`Q`
collapse (that step assumes a single constant `θ` and does not survive variable
`θ_i`). Spencer's eq (3) `M_o` and its load arms are still the right *reference*
for the external-load moment terms, since those forces and arms are unchanged.

### 4a. Moment-about-origin sum, term by term

**Convention.** Counter-clockwise positive (right-hand rule), matching Spencer.
The moment about the origin of a force with global components `(F_x, F_y)` acting
at point `(x, y)` is

>  `M = x · F_y − y · F_x`

A purely vertical force then needs only its `x`; a purely horizontal force only its
`y`. The component signs below are taken **directly from Spencer eqs (1)–(2)** so
the two methods share one convention; the column names are the actual `slice_df`
columns (verified). `α` = base inclination, `β` = top/load inclination,
`ψ = α` = reinforcement inclination (flexible, parallel to base), `θ_p` = pile
inclination.

The forces, points of action, and angle conventions below are exactly those of the
Spencer slice free-body diagram — refer to it for every symbol:

![Spencer slice forces](../docs/lem/images/spencer3_forces.png)

(Source: [`docs/lem/spencer.md`](../docs/lem/spencer.md) → "Slice Geometry and
Forces"; `W, N, S, E, X, kW, P/D, T, R, V, H, θ_p, β, α, Δℓ, Δx` are all defined
there.)

| # | Force (per slice) | `(F_x, F_y)` | acts at `(x, y)` | Moment about origin `M = x F_y − y F_x` |
|---|---|---|---|---|
| 1 | Weight `W` | `(0, −W)` | `(x_c, y_cg)` | `−W · x_c` |
| 2 | Effective base normal `N′` (`n_eff`) | `(−N′ sinα, N′ cosα)` | `(x_c, y_cb)` | `N′ (x_c cosα + y_cb sinα)` |
| 2b | Pore-water uplift `U = u·dl` | `(−U sinα, U cosα)` | `(x_c, y_cb)` | `U (x_c cosα + y_cb sinα)` |
| 3 | Base shear `S` | `(S cosα, S sinα)` | `(x_c, y_cb)` | `S (x_c sinα − y_cb cosα)` |
| 4 | Surface load `D` (`dload`) | `(D sinβ, −D cosβ)` | `(d_x, d_y)` | `−D (d_x cosβ + d_y sinβ)` |
| 5 | Seismic `kW` | `(−kW, 0)` | `(x_c, y_cg)` | `kW · y_cg` |
| 6 | Tension-crack water `V` (`t`, last slice) | `(−V, 0)` | `(x_v, y_t)` | `V · y_t` |
| 7 | Reinforcement `R` (`p`) | `(R cosα, R sinα)` | `(x_c, y_cb)` | `R (x_c sinα − y_cb cosα)` |
| 8 | Pile `H` (`h_pile`) | `(H cosθ_p, H sinθ_p)` | `(x_pile, y_pile)` | `H (x_pile sinθ_p − y_pile cosθ_p)` |

Rows 2 + 2b are the total base reaction: the march solves the **effective** normal
`N′` (= `n_eff`) and carries the pore-water uplift `U = u·dl` as a separate force
in the *same* `(−sinα, cosα)` direction, so the total normal moment is
`(N′ + u·dl)(x_c cosα + y_cb sinα)` — keep them split so the term maps directly
onto the march's `n_eff` output. The shear magnitude is
`S = c_i·dl_i + (N′_i)·tanφ_i)/F` (`= [c·dl + (N_total − u·dl)tanφ]/F`).

The **interslice forces `Z_i` do not appear** — they are internal and cancel
exactly in the whole-mass sum (§4, the key fact). That is what lets us sum about
the origin without knowing the thrust line.

**Reconciliation against `force_equilibrium` (done).** The `N`/`S` directions
above are derived from `force_equilibrium`'s 2×2 matrix: the `N`-coefficient
`tan_phi_m·cosα − sinα` factors as `(−sinα) + tan_phi_m·(cosα)`, giving
`N → (−sinα, cosα)` and `S → (cosα, sinα)`; row 1 confirms the `y`-components.
Every external sign in the table matches the engine's `b`-vector
(`+D sinβ`, `−kW`, `−V`, `+R cosα`, `+H cosθ_p`). The one correction surfaced by
this check is the effective-vs-total normal (rows 2/2b): the engine's `N` is
effective, so the moment uses `N′ + u·dl`.

**Assembly.** Moment equilibrium of the whole mass:

>  `M_origin(F, λ) = Σ_i [ M(W) + M(N) + M(S) + M(D) + M(kW) + M(V) + M(R) + M(H) ]_i = 0`

The mobilized base shear is `S_i = [ c_i·dl_i + (N_i − u_i·dl_i)·tanφ_i ] / F`, so
`F` enters **explicitly** through term 3 (the only `1/F` term) and **implicitly**
through every `N_i`, which the march produces for the trial `(F, λ)`. Hence
`F_m(λ)` is the root in `F` of `M_origin(F, λ) = 0` — exactly the `moment_res`
that `_mp_march` returns (§3). Grouping the explicit `1/F` gives the familiar
"resisting strength moments / driving moments" shape, but because `N_i` comes from
the numerical march we solve the residual directly rather than as a closed form.

**Notes / things to confirm in code:**
- Weight and seismic use `x_c` for the horizontal line of action (no `x_cg`
  column exists; thin-slice assumption, same as Spencer's `x_b = x_c`).
- `V` acts only on the last slice (`t` is zero elsewhere); `x_v` doesn't enter
  (horizontal force) — only `y_t` matters.
- The `N`/`S` global-component signs (rows 2–3) are **reconciled** against
  `force_equilibrium`'s 2×2 matrix and `b`-vector (see the reconciliation note
  above): `N → (−sinα, cosα)`, `S → (cosα, sinα)`, and the engine's `N` is the
  **effective** normal — the moment term must use `N′ + u·dl`. The
  `f(x)=1 ≡ Spencer` test (§7.1) remains the backstop for any residual error.
- Right-facing slopes: apply the same reflection `force_equilibrium`/`spencer`
  already use (negate the inclinations / flip the relevant signs) so `M_origin`
  stays mirror-consistent. Validate with the `f(x)=1 ≡ Spencer` test run on a
  right-facing geometry specifically.
- A second surface load (`dload2`/`d_x2`/`d_y2`) exists for rapid drawdown;
  include it by analogy with row 4 if/when M-P supports rapid drawdown.

### Option (b) — per-slice line-of-thrust recurrence (original M-P 1965)

March the per-slice moment equation (generalize Spencer eqs (68)–(69) to
per-boundary θ) to track the thrust-line position `y_t`, and require moment
closure / an admissible thrust line. This is the historical formulation and
yields the thrust line as the primary product.

- **Cons:** more delicate numerically; the thrust line can run outside the slice
  on marginal surfaces; defining "the" FS from it is fiddlier than option (a).

### Proposed split

Use **option (a)** to *define* `F_m` (it's simpler, matches modern GLE codes, and
its crossing with `F_f` is the M-P answer). Then compute the **line of thrust as a
post-process diagnostic** via the option-(b) recurrence once `(F, λ)` are known —
giving us the same plotting/admissibility output Spencer produces, without letting
the delicate thrust march drive the FS. **This split is the main thing to confirm
before coding.**

## 5. The interslice force functions `f(x)`

**Normalized position.** `f` is evaluated at each slice boundary `x_j`
(`j = 0..n`) on a normalized coordinate spanning the slip surface:

>  `s_j = (x_j − x_L) / (x_R − x_L) ∈ [0, 1]`

where `x_L`, `x_R` are the left/right ends of the slip surface (`x_j` at `j=0`
and `j=n`). The interslice inclination at each boundary is then
`θ_j = arctan(λ · f(s_j))`.

**Scale-invariance (important).** Because `tan θ = λ·f`, multiplying `f` by any
constant just rescales `λ` — the `(F, λ·f)` *product* is what the solution fixes.
So **only the shape of `f` matters, not its amplitude**; normalizing every `f` to
a peak of 1 is pure convention. This is also why FS is famously insensitive to the
choice of `f` (only the *relative variation* of the interslice angle along the
surface moves the answer, and that effect is small).

**v1 library** (selected by `f_type`):

1. **`constant`** — `f(s) = 1`.
   Recovers Spencer (`tan θ = λ`, constant). This is the regression case (§7.1).

2. **`half_sine`** — `f(s) = sin(π s)`.
   Single arch, peak 1 at mid-surface, `f = 0` at both ends so the interslice
   *shear* vanishes where the slip surface daylights — physically reasonable and
   the textbook / GeoStudio default. This is the workhorse for design.

3. **`clipped_sine`** — quarter-sine tapers with a flat unit top, taper fraction
   `t ∈ (0, 0.5]` (default `t = 0.25`):

   ```
   f(s) = sin( π s / (2t) )          0 ≤ s < t        (rise 0 → 1)
        = 1                          t ≤ s ≤ 1 − t    (flat)
        = sin( π (1 − s) / (2t) )    1 − t < s ≤ 1    (fall 1 → 0)
   ```
   `C¹`-continuous at the joins (slope 0 there). `t → 0` approaches `constant`;
   `t = 0.5` is the half-sine. A flat-topped function that still tapers to zero at
   the ends.

   `trapezoidal` is the same shape with **linear** ramps instead of sine tapers
   (corner fractions `a`, `b`; default `a = b = 0.25`):
   ```
   f(s) = s / a            0 ≤ s < a
        = 1                a ≤ s ≤ 1 − b
        = (1 − s) / b      1 − b < s ≤ 1
   ```
   Offered as a sub-variant of `clipped_sine` (or its own `f_type` — minor).

All four taper to zero at the ends except `constant`. Since FS is insensitive to
`f` and the design recommendation is `half_sine`, `clipped_sine`/`trapezoidal`
are provided mainly for completeness and inter-code comparison; their exact corner
parameters are **not** tuned to match any one external code (see §8).

**Right-facing slopes.** `force_equilibrium` already negates `theta_list` and
  flips the tension-sign convention for right-facing surfaces; M-P's λ-driven
  angles must follow the same convention so results stay mirror-symmetric. `f(x)`
  itself is symmetric, so this is just the existing sign handling applied to
  `arctan(λ f)`.
- **Initial guesses & robustness.** Seed `F` from Bishop (as Spencer does) and
  `λ = 0` (pure force-equilibrium starting point), ramp `λ` up; reuse Spencer's
  tension-admissibility rejection idea on the resulting `N_i`.

## 6. Integration checklist (after the math is settled)

- `solve.py`: new `morgenstern_price(slice_df, f_type='half_sine', ...)` returning
  `(success, {'method','FS','lambda', 'f_type', ...})`, following the
  `(success, result)` contract.
- `solve_selected()` / `solve_all()`: register the method + a print line
  (`FS`, λ, f-type).
- Interslice-function selection: `f_type` function argument (decided — function
  arg only for v1; no input-template cell).
- `plot.py`: show λ / interslice function and line of thrust (mirror Spencer).
- `search.py`: confirm the circular/non-circular search can call it like Spencer.
- Docs: new `docs/lem/morgenstern_price.md` (derivation in the Spencer style),
  link from `docs/lem/overview.md`, update `verification.md`.
- Per the repo's "docs track solver" rule, this is one work unit:
  code + docstrings + docs page + verification + sample.

## 7. Validation plan

### 7.1 Internal consistency

1. **`f(x)=1` must equal Spencer** on every existing LEM benchmark (FS and the
   implied constant `θ` — i.e. `arctan(λ)` should match Spencer's `θ`). Tight
   tolerance. This is the make-or-break test and confirms the whole formulation.
2. **`f(x)` insensitivity** — FS should vary <~1% across constant / half-sine /
   clipped on a normal slope; flag if it doesn't (a classic M-P property).
3. **Force-only tie-out** — `F_f(λ)` from `_mp_march` must reproduce the existing
   `force_equilibrium` FS for the same `theta_list`, as a check on the refactor.

### 7.2 Published Morgenstern–Price benchmarks (already in our V&V suite)

These reuse the LEM benchmark models already built by `benchmarks/build_lem.py`
and documented in `docs/verification.md` / `docs/benchmark-run-plan.md`. The
reference M-P/GLE values come from the same published sources we already cite.

| Case | Surface | M-P / GLE reference | xslope Spencer | Source |
|---|---|---|---|---|
| ACADS simple homogeneous (LEM-1) | circular | ≈ 1.00 (ACADS consensus) | 0.984 | SLOPE/W Verification Manual (Oct 2022) |
| **ACADS weak layer (LEM-2)** | **non-circular** | **1.261** (named Morgenstern–Price, SLOPE/W) | 1.258 | SLOPE/W Verification Manual (Oct 2022) |
| Arai & Tagyo 1985 (LEM-2b) | circular | ≈ 1.451 (Spencer/GLE; records vary 1.40–1.45) | 1.401 | Arai & Tagyo (1985); Greco (1996); Malkawi et al. (2001) |
| Rocscience / Slide Ex. 1 (optional) | circular + foundation | GLE = 0.987 | 0.987 | Rocscience "FEM applied to slope stability" |

**Primary M-P target: LEM-2 (ACADS weak layer), non-circular, SLOPE/W M-P =
1.261.** This is the most valuable check — non-circular + complete equilibrium is
exactly where M-P earns its keep, and SLOPE/W reports the named M-P value (not
just Spencer). xslope Spencer already sits at 1.258 here, so M-P with half-sine
should land within ~1% of 1.261.

Acceptance: half-sine M-P within ~1% of each reference FS; FS spread across
`f(x)` choices <~1%; `f(x)=1` reproduces our own Spencer to tolerance.

**To source if possible (improves the benchmark):** a published case that also
reports **λ** for a stated `f(x)` (the SLOPE/W manual and Fredlund–Krahn 1977
both plot/tabulate λ). Matching λ — not just FS — would validate the interslice
distribution, not only the scalar result. Open item — see §8.

## 8. Still to settle before coding

1. **`F_m` formulation (§4):** confirm the proposed split — option (a)
   overall-moment-**about-the-origin** *defines* `F_m` (same reference point as
   Spencer, works for circular and non-circular alike), option (b) thrust-line
   recurrence is a post-process diagnostic. This is the main open item.
2. **Clipped/trapezoidal default parameters** (taper/corner fractions `t`, `a`,
   `b`) — defined in §5 with sensible defaults; only needs a final value choice,
   and they're deliberately *not* matched to a specific external code.
3. **A published case reporting λ** (not just FS) for a stated `f(x)`, to validate
   the interslice distribution — the SLOPE/W manual and Fredlund–Krahn (1977) are
   the likely sources (§7.2). Nice-to-have, not a blocker.

## 9. Resolved (from discussion 2026-06-23)

- f(x) library = constant + half-sine + clipped/trapezoidal sine; no user table v1.
- f(x) selected via `f_type` function argument; no Excel cell v1.
- Build the shared `_mp_march` evaluator, then implement **both** Approach A
  (two-curve intersection, reference impl) and Approach B (2-D Newton, fast path
  that must agree with A).
- `F_m` sums moments about the **origin** (not the circle center), matching
  Spencer's eq (16) so the same formulation covers circular and non-circular
  surfaces.
