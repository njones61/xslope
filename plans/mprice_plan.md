# Morgenstern–Price Method — Implementation Plan / Discussion

Status: **BUILD SPEC — design decisions locked (2026-06-24), nothing implemented
yet.** All blocking decisions are resolved; see the Decisions block below, the
resolved list (§9), and the build order (§10). Remaining items (§8) are
non-blocking. Ready to implement.

### Decisions (locked 2026-06-24) — this section is now a build spec

- **`f(x)` library at v1:** constant (`=1`, Spencer regression) and half-sine
  *only*. Clipped/trapezoidal sine and a user-defined table are both deferred to
  "later/optional" (§5a) — FS is insensitive to `f(x)`, so they buy ~no accuracy
  and only add parameters/tests; trivial to add once `_mp_march` exists.
- **`f(x)` selection:** function argument only (e.g. `f_type='half_sine'`); no
  Excel-template cell for now.
- **`F_m` (moment) formulation:** **locked** — overall moment about the **origin**
  (§4a, fully derived and sign-reconciled). Thrust line is a post-process
  diagnostic (option b). The `f(x)=1 ≡ Spencer` test is the validation backstop.
- **Shared building block:** **refactor `force_equilibrium`** to expose a fixed-`F`
  per-slice march returning `(force_res, moment_res)` (no inner `F` root-find);
  Corps/L-K/M-P all call it. Guard Corps/L-K FS bit-stability with the existing
  regression across the refactor.
- **Solver structure:** **Approach A (two-curve crossing) is the reference oracle,
  the `F`-vs-`λ` docs figure, and the stability fallback — not the shipped path.**
  **Approach B (2-D Newton on `(F, λ)`) is the production solver** (Bishop seed,
  `λ=0` start, line search, tension guard — mirroring `spencer()`); on divergence
  it falls back to A's bracketed crossing. Start B with a **numeric** Jacobian:
  unlike Spencer's closed-form lumped-`Q`, M-P's residuals come from the march, so
  the analytic Jacobian means differentiating the recurrence — defer it unless the
  numeric Newton proves slow/unstable on the benchmarks (the Spencer history says
  it might, which is exactly why A exists as the robust fallback).
- **v1 load scope:** all single-stage terms in the §4a moment table
  (`W, dload, kW, V, R, H`). **Rapid drawdown needs no solver work** — it is an
  outer 3-stage wrapper (`advanced.rapid_drawdown`) that swaps each stage's loads /
  pore pressures into the standard `dload`/`u` columns and re-calls the solver, so
  M-P inherits it for free once it solves a standard `slice_df`. (There is **no**
  `dload2` term in the M-P moment sum — the solver never sees the `2` columns.)
- **Result contract:** mirror `spencer()`'s output keys exactly, plus `lambda` and
  `f_type`, so `plot.py`, the summary, and the search reuse Spencer's rendering
  unchanged: `(success, {method, FS, lambda, f_type, theta_i (per boundary),
  line-of-thrust `y_t`, `n_eff`, `z`, ...})`.
- **No-solution / robustness policy:** mirror Spencer — Bishop `F` seed, `λ=0`
  start, reuse `force_equilibrium`'s tension-admissibility guard on `N_i`, and on
  non-convergence or no `F_f`/`F_m` crossing return `(False, reason)` so the
  automated search just skips the surface. Approach-A `λ` search bracket:
  `λ ∈ [−1.5, 1.5]` with a sign-change bracket on `g(λ)=F_f−F_m`.

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
  reasonable and is GeoStudio's/SLOPE-W's default.

These two are the v1 library. **Clipped/trapezoidal sine** and a **user-defined
table** are deferred (§5a): because FS is insensitive to `f(x)`, they add solver
parameters and test surface without moving any benchmark, and they are easy to
bolt on later once the engine exists.

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

### Recommendation (locked)

Build `_mp_march` first (the `force_equilibrium` refactor). Then:

- **Approach A is the reference oracle, the docs `F`-vs-`λ` figure, and the
  stability fallback** — bracketed and robust-but-slow. Build it first so it can
  validate everything else and back up B on divergence.
- **Approach B (2-D Newton on `(F, λ)`) is the shipped solver.** It must agree with
  A on every benchmark; on divergence it falls back to A's bracketed crossing.
  Start with a numeric Jacobian (analytic march-Jacobian deferred — see the
  Decisions block). Seed both from Bishop FS and `λ = 0`.

Lesson carried over from Spencer: a numeric/curve-crossing solver on these LEM
residuals was slow and unstable until an analytic-Jacobian Newton replaced it.
Here, A is deliberately retained as the robust fallback so B can be the fast path
without first paying for the (harder, march-based) analytic Jacobian.

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
- **No `dload2` term.** Rapid drawdown is handled entirely by the outer wrapper
  `advanced.rapid_drawdown`, which swaps each stage's `dload2`/`u2`/… into the
  standard `dload`/`u` columns and re-calls the solver. The M-P solver only ever
  sees one set of loads, so the moment sum needs the single `dload` term (row 4)
  and nothing more — M-P inherits rapid drawdown for free.

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

**v1 library** (selected by `f_type`) — **two functions**:

1. **`constant`** — `f(s) = 1`.
   Recovers Spencer (`tan θ = λ`, constant). This is the regression case (§7.1).

2. **`half_sine`** — `f(s) = sin(π s)`.
   Single arch, peak 1 at mid-surface, `f = 0` at both ends so the interslice
   *shear* vanishes where the slip surface daylights — physically reasonable and
   the textbook / GeoStudio / SLOPE-W default. This is the workhorse for design.

`half_sine` tapers to zero at the ends; `constant` does not. These two span the
practical range (Spencer + canonical M-P), which is exactly why FS-vs-`f(x)`
insensitivity makes a larger library unnecessary for v1.

**Right-facing slopes.** `force_equilibrium` already negates `theta_list` and
  flips the tension-sign convention for right-facing surfaces; M-P's λ-driven
  angles must follow the same convention so results stay mirror-symmetric. `f(x)`
  itself is symmetric, so this is just the existing sign handling applied to
  `arctan(λ f)`.
- **Initial guesses & robustness.** Seed `F` from Bishop (as Spencer does) and
  `λ = 0` (pure force-equilibrium starting point), ramp `λ` up; reuse Spencer's
  tension-admissibility rejection idea on the resulting `N_i`.

### 5a. Deferred `f(x)` options (not v1)

Kept here so the work isn't lost, but **out of the v1 library** — FS is insensitive
to `f(x)`, so these move no benchmark; they exist only for completeness /
inter-code comparison and are a few lines to add once `_mp_march` and the `f_type`
dispatch exist. Their corner parameters are deliberately **not** tuned to match any
one external code.

- **`clipped_sine`** — quarter-sine tapers with a flat unit top, taper fraction
  `t ∈ (0, 0.5]` (default `t = 0.25`):

  ```
  f(s) = sin( π s / (2t) )          0 ≤ s < t        (rise 0 → 1)
       = 1                          t ≤ s ≤ 1 − t    (flat)
       = sin( π (1 − s) / (2t) )    1 − t < s ≤ 1    (fall 1 → 0)
  ```
  `C¹`-continuous at the joins (slope 0 there). `t → 0` approaches `constant`;
  `t = 0.5` is the half-sine.

- **`trapezoidal`** — same flat-topped shape with **linear** ramps instead of sine
  tapers (corner fractions `a`, `b`; default `a = b = 0.25`):
  ```
  f(s) = s / a            0 ≤ s < a
       = 1                a ≤ s ≤ 1 − b
       = (1 − s) / b      1 − b < s ≤ 1
  ```

- **user-defined table** — a `(s, f)` lookup interpolated to the boundaries. The
  most flexible and the most input-plumbing; deferred indefinitely.

## 6. Integration checklist (after the math is settled)

- `solve.py`: new `mprice(slice_df, f_type='half_sine', ...)` returning
  `(success, {'method','FS','lambda', 'f_type', ...})`, following the
  `(success, result)` contract.
- `solve_selected()` / `solve_all()`: register the method + a print line
  (`FS`, λ, f-type).
- Interslice-function selection: `f_type` function argument (decided — function
  arg only for v1; no input-template cell).
- `plot.py`: show λ / interslice function and line of thrust (mirror Spencer).
- `search.py`: confirm the circular/non-circular search can call it like Spencer.
- Docs: new `docs/lem/mprice.md` (derivation in the Spencer style); **add it to the
  `mkdocs.yml` nav** under the LEM section, right after "Spencer's Method"; link
  from `docs/lem/overview.md`; update `verification.md`.
- Per the repo's "docs track solver" rule, this is one work unit:
  code + docstrings + docs page + verification + sample.

## 7. Validation plan

### 7.1 Internal consistency

1. **`f(x)=1` must equal Spencer** on every existing LEM benchmark (FS and the
   implied constant `θ` — i.e. `arctan(λ)` should match Spencer's `θ`). Tight
   tolerance. This is the make-or-break test and confirms the whole formulation.
2. **`f(x)` insensitivity** — FS from `constant` vs `half_sine` should differ
   <~1% on a normal slope; flag if it doesn't (a classic M-P property). (Extends
   to the deferred `f(x)` shapes if/when they're added.)
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

Acceptance: half-sine M-P within ~1% of each reference FS; `constant` vs
`half_sine` FS spread <~1%; `f(x)=1` reproduces our own Spencer to tolerance.

**To source if possible (improves the benchmark):** a published case that also
reports **λ** for a stated `f(x)` (the SLOPE/W manual and Fredlund–Krahn 1977
both plot/tabulate λ). Matching λ — not just FS — would validate the interslice
distribution, not only the scalar result. Open item — see §8.

## 8. Open items

All design decisions are locked (see the Decisions block at the top and §9). The
only remaining items are non-blocking and can be pursued during/after coding:

1. **A published case reporting λ** (not just FS) for a stated `f(x)`, to validate
   the interslice distribution — the SLOPE/W manual and Fredlund–Krahn (1977) are
   the likely sources (§7.2). Nice-to-have, not a blocker.
2. **Analytic march-Jacobian for Approach B** — only if the numeric Jacobian proves
   slow/unstable on the benchmarks (the A fallback covers stability in the
   meantime). Deferred, not blocking.

## 9. Resolved decisions (locked 2026-06-24)

- **f(x) library v1** = constant + half-sine only; clipped/trapezoidal sine and the
  user-defined table deferred (§5a) — FS is insensitive to f(x).
- **f(x) selection** via `f_type` function argument; no Excel cell v1.
- **Shared building block:** refactor `force_equilibrium` into a fixed-`F` march
  returning `(force_res, moment_res)`; Corps/L-K/M-P all call it (guard their FS
  stability with the existing regression).
- **`F_m`** sums moments about the **origin** (§4a, derived + sign-reconciled),
  covering circular and non-circular alike; thrust line is a post-process
  diagnostic. **No `dload2` term** — rapid drawdown is an outer wrapper.
- **Approach A** = reference oracle + `F`-vs-`λ` docs figure + stability fallback.
  **Approach B (2-D Newton on `(F, λ)`)** = shipped solver, numeric Jacobian to
  start, falls back to A on divergence. Both seeded from Bishop FS and `λ=0`.
- **v1 load scope** = single-stage `W, dload, kW, V, R, H`; rapid drawdown inherited
  for free via the outer wrapper (no solver work).
- **Result contract** mirrors `spencer()`'s keys + `lambda` + `f_type`.
- **No-solution policy** mirrors Spencer: tension guard, return `(False, reason)`
  on non-convergence; Approach-A `λ` bracket `[−1.5, 1.5]`.

## 10. Build order (ready to implement)

The **core** is a strictly sequential chain on `solve.py` (S0→S5). Each stage lists
its **deliverable**, **depends-on**, and **exit criterion** (the gate that must pass
before the next stage starts). Stages marked ⛔ are **hard stops** — do not proceed
past them on a failure; a defect there poisons everything downstream. Independent
work that can run in parallel with the core is in §11.

**S0 — Recon (no code).**
- *Deliverable:* confirmed map of `force_equilibrium`'s return signature + its
  Corps/L-K callers, and a verified checklist that every `slice_df` column the §4a
  moment table names exists with the assumed sign (`n_eff`, `y_cb`, `d_x/d_y`, `t`,
  `p`, `h_pile`, `theta_p`, …).
- *Depends-on:* nothing. *Exit:* checklist complete; any column/sign surprise
  reconciled against §4a before S1. (Good fit for an `Explore` subagent — §11.)
- **Status (2026-06-24): DONE — gate met, no column/sign mismatches vs §4a.**
  Findings to carry into S1/S2:
  - **Extractable march** = the `force_equilibrium` inner loop **`solve.py:648–682`**
    — the per-slice 2×2 solve (`A·[N_i, Z_{i+1}] = b`) run at a fixed `F`, *before*
    the `scipy.optimize.newton` root-find of `F` on residual `Z[n]`. It writes
    `n_eff` (effective normal) and `z` back to `slice_df`. Refactor target is clean.
  - **`spencer()` success dict is minimal:** `{'method','FS','theta'}`. Line of
    thrust (`yt_l`,`yt_r`), `n_eff`, `z`, `theta` are stored in `slice_df`, not
    returned. → M-P result dict = those keys **+ `lambda` + `f_type`** (§contract).
  - ⚠️ **Right-facing conventions differ between the two engines** (highest-risk
    reconciliation): `force_equilibrium` expects the *caller* to negate `theta_list`
    and uses the `right_facing` flag only for the tension-sign guard — it does **not**
    flip external-load signs internally; `spencer()` instead flips
    `alpha,beta,psi,R,c,kw,V,tan_p,H_cos` *internally* (solve.py:963–991). The M-P
    march is `force_equilibrium`-based, so the moment accumulator (S2) must hold **one**
    consistent convention; run the S3 `f(x)=1 ≡ Spencer` gate on a **right-facing**
    geometry specifically to catch a sign slip here.
  - **Gotcha:** `theta_p` (pile) is stored in **radians**; all other `slice_df`
    angles (`alpha,phi,beta,theta`) are **degrees**.
  - Total-normal moment term uses **`n_eff + u·dl`** (the engine's `n_eff` is the
    effective normal).

**S1 — Refactor `force_equilibrium` into the shared fixed-`F` march.**
- *Deliverable:* `_mp_march(slice_df, lam, f_vals, F) -> (N[], Z[], force_res,
  moment_res=None)`; Corps/L-K re-pointed through it.
- *Depends-on:* S0. *Exit ⛔:* every Corps/L-K LEM benchmark FS **bit-stable**
  (identical pre/post refactor). Do not start S2 until green.
- **Status (2026-06-24): DONE — gate met, bit-identical.** Design refinement vs the
  deliverable above: the shared low-level march was extracted as
  **`_equilibrium_march(alpha, phi, c, w, u, dl, D, beta, kw, T, P, H_pile,
  theta_p, theta, FS) -> (N, Z)`** (solve.py, just above `force_equilibrium`). It
  takes the *per-boundary `theta`* array (general — Corps/L-K build it from
  geometry; M-P will pass `arctan(λ·f)`), not `(lam, f_vals)`, and returns just
  `(N, Z)`. `force_res = Z[n]` and the `moment_res` accumulator belong to the M-P
  wrapper built in **S2/S3** on top of this march — keeping the low-level march
  free of any residual/`λ` concepts. `force_equilibrium.residual()` now delegates
  to it; **Corps/L-K were not touched** (they call `force_equilibrium`, so they are
  re-pointed transitively). The march carries a prominent SIGN/RIGHT-FACING comment
  block (caller-negates-theta convention, contrasted with `spencer()`'s internal
  flips). Verified bit-identical FS to full float precision over 16 Corps/L-K
  evaluations (piles, reinforcement, submerged, non-circular, multi-layer).

**S2 — Moment accumulator.**
- *Deliverable:* `moment_res` populated in the march from the §4a table.
- *Depends-on:* S1. *Exit:* each of the 8 terms unit-checked in isolation (sign +
  arm) against §4a; assembled `moment_res` finite on a sample surface.
- **Status (2026-06-24): DONE.** Added `_moment_about_origin(...)` (the §4a sum,
  convention-agnostic — uses whatever arrays it's given) and
  `_mp_march(slice_df, lam, f_vals, FS) -> (N, Z, force_res, moment_res)`
  (left-facing; builds `θ=arctan(λ·f)`, runs `_equilibrium_march`, then the moment
  sum). `right_facing=True` raises `NotImplementedError` (finalized in S3). All 8
  moment terms unit-checked to 1e-9 in isolation (sign + arm). **Bonus result that
  de-risks S3:** at Spencer's solution on a left-facing benchmark (acads_simple,
  FS=1.4462, θ=15.35°) with `f≡1`, `λ=tan θ`, BOTH residuals vanish to numerical
  precision (`force_res≈2e-8`, `moment_res≈-6e-7` vs a moment scale ~1.8e5, i.e.
  ~3e-12 relative). So the `f(x)=1 ≡ Spencer` *formulation* is already confirmed for
  left-facing — S3's left-facing gate should pass immediately; only the right-facing
  convention remains to work out and validate.

**S3 — Approach A + the Spencer gate.**
- *Deliverable:* `F_f(λ)`, `F_m(λ)` root-finds and the `λ` crossing; the `F`-vs-`λ`
  curve as a diagnostic.
- *Depends-on:* S2. *Exit ⛔ (make-or-break):* with `f(x)=1`, M-P reproduces
  `spencer()` **FS and `arctan(λ)` vs Spencer θ** to tight tolerance on every LEM
  benchmark, **and** the force-only tie-out (`F_f` vs `force_equilibrium`) matches.
  Nothing downstream is trustworthy until this passes.
- **Status (2026-06-24): DONE ✅ — gate passed, S3a + S3b.** Added public
  `mprice(slice_df, f_type='half_sine', ...)` + `_mp_f_vals`
  (constant / half_sine).
  - **S3a robustness finding:** the moment-only FS curve `F_m(λ)` is *multivalued*
    (an asymptote where its branch flips), so a naive `F_f`−`F_m` crossing jumps
    branches and fakes/misses crossings. Fix: root-find on
    **`h(λ) = moment_res(F_f(λ), λ)`** — the moment residual evaluated AT the
    force-equilibrium FS. `force_res` is monotonic in FS so `F_f(λ)` is smooth and
    single-valued, making `h(λ)` a smooth scalar in λ alone whose root is the M-P
    solution; it sidesteps `F_m`'s multivaluedness entirely.
  - **S3b right-facing:** determined empirically (both residuals vanish at Spencer's
    solution) that the right-facing convention is exactly Spencer's internal flip
    set — `_mp_march` flips `{alpha, beta, phi, c, P, kw, V}` (+ pile `θp→π−θp` by
    analogy, no right-facing+pile benchmark exists) when `right_facing`, applied
    once before both march and moment so they stay consistent.
  - **Result: 20/20 benchmarks (all left + right facing) reproduce Spencer EXACTLY**
    (FS to ±1e-6, `arctan(λ)` = Spencer θ); half_sine FS sits <~1.4% from constant
    (the §7.1 insensitivity property). The force march is unchanged from S1, so
    Corps/L-K remain bit-stable (force-only tie-out holds).

**S4 — `half_sine` + published benchmarks.**
- *Deliverable:* `half_sine` `f_type`; the §7.2 benchmark runs.
- *Depends-on:* S3 (gate passed). *Exit:* half-sine M-P within ~1% of each
  reference (primary LEM-2, SLOPE/W M-P = 1.261); `constant`-vs-`half_sine` FS
  spread <~1%.
- **Status (2026-06-24): DONE ✅** (validation only — `half_sine` landed in S3).
  Evaluated on each benchmark's critical surface:
  | Case | ref | xslope Spencer | MP half_sine | MP const |
  |---|---|---|---|---|
  | LEM-1 ACADS simple (circ) | 1.00 | 0.984 | 0.984 | 0.984 |
  | **LEM-2 weak layer (non-circ)** | **1.261** | 1.259 | **1.251 (−0.8%)** | 1.259 |
  | LEM-2b Arai & Tagyo (circ) | 1.451 | 1.401 | 1.400 | 1.401 |
  Primary target met: half_sine M-P within **−0.8%** of SLOPE/W's named non-circular
  M-P value. `constant`-vs-`half_sine` spread <1% everywhere (insensitivity holds).
  LEM-1/LEM-2b sit >1% from the *external* reference, but xslope's own Spencer is
  offset by the same amount (pre-existing V&V gap, not an M-P error) — M-P matches
  xslope Spencer exactly in every case.

**S5 — Approach B (shipped solver).**
- *Deliverable:* 2-D Newton on `(F, λ)` (numeric Jacobian, Bishop+`λ=0` seed, line
  search, tension guard), falling back to A on divergence.
- *Depends-on:* S4. *Exit:* B matches A on **every** benchmark; no-solution path
  returns `(False, reason)`; search runs M-P over many surfaces without hangs.
- **Status (2026-06-24): DONE ✅.** Approach B is the plan's **2-D Newton** (scipy
  `hybr`) on the **scaled** `(force_res, moment_res)` residuals over `(FS, λ)`,
  seeded at `(Bishop-FS, 0)`. Both residuals come from ONE `_mp_march` per eval.
  Scaling the two residuals to O(1) (they differ ~1e3×) makes `hybr`
  well-conditioned, and seeding near the physical solution keeps it off the
  multivalued moment-only branch. `mprice(..., solver=)`: `'auto'`
  (default) runs B then falls back to the robust Approach-A grid scan; `'A'`/`'B'`
  force one. **Validation: A vs B agree to max |ΔFS|=1.5e-9, |Δλ|=7e-10 over 41
  cases, 0 fallbacks.** (An earlier 1-D `h(λ)` secant also worked but was ~3× slower.)
- **Speed optimization (2026-06-24).** An early "≈ Spencer speed" claim was wrong —
  it came from one input (acads_simple) that happened to show parity. The real
  baseline on the main_lem problem (sloping_bottom, 40 slices) was an M-P
  circular_search at **14 s vs Spencer 3.3 s (~4×)**. Profiling drove four fixes:
  1. **Closed-form 2×2** in `_equilibrium_march` (replace per-slice `np.linalg.solve`).
  2. Approach B **returns its FS**, skipping a redundant final `F_f` re-solve.
  3. **Vectorized march**: the interslice recurrence `Z[i+1]=p[i]+q[i]·Z[i]` is
     LINEAR, so the whole march is `cumprod`/`cumsum` — no Python loop. Stable (the
     multiplier `q`≈1, cumprod ~1 even at 120 slices; scalar-scan fallback if
     non-finite). Speeds up Corps/L-K too (unchanged to 3.1e-14).
  4. **Extract-once** (the big one): split `_mp_march` into `_mp_extract` (pull ~25
     DataFrame columns once, with flips) + `_mp_residuals` (pure-numpy, per
     iteration). Profiling showed `DataFrame.__getitem__` (re-read every Newton
     step, ~13×/solve) was ~75% of the solve.
  **Result: M-P circular_search 3.7 s vs Spencer 3.3 s (~1.12×) — essentially
  parity**, all gates still pass.
- **Analytic Jacobian — tried, reverted (no benefit).** Implemented a gradient-checked
  (to 1e-8) tangent-linear Jacobian + hand-coded 2-D Newton. Correct and gates pass,
  but **no speedup**: the bottleneck was never the iteration count — it was the
  per-call work (DataFrame extraction, then the march), so a tangent march that does
  MORE per-slice arithmetic just trades fewer iterations for costlier ones. Reverted
  to keep the simpler scipy `hybr` Approach B. (The extract-once + vectorized-march
  wins above are what actually mattered.)

**S6 — Integration (code).**
- *Deliverable:* register `mprice` in `solve_selected`/`solve_all`
  (+ a print line: FS, λ, `f_type`); `plot.py` shows λ / the interslice function and
  the line of thrust (mirror Spencer); confirm the circular/non-circular search
  calls it exactly like Spencer.
- *Depends-on:* S5. *Exit:* `solve_all` and an automated search run M-P end-to-end
  and plot correctly.
- **Status (2026-06-24): DONE ✅.** Registered in `solve_selected` (prints
  `Morgenstern-Price (f_type): FS=..., lambda=...`) and appended to `solve_all`.
  `plot_solution` gains an M-P title (`FS`, `λ`); base-stress bars already work
  from the stored `n_eff`. Search needs no change — `circular_search` /
  `noncircular_search` call it via `getattr(solve, method)` exactly like Spencer
  (confirmed: acads_simple circular_search critical FS=0.984). Verified `solve_all`,
  `solve_selected`, and `plot_solution` end-to-end on circular + non-circular.
- **Line of thrust DONE (2026-06-24).** `_mp_line_of_thrust` generalizes Spencer's
  per-slice moment recurrence (eqs 68-69) to a per-boundary θ, using the march's
  `Z` and the same external-load `Mo` + right-facing reflection Spencer uses; for a
  constant θ it reduces exactly to Spencer's. Stored as `yt_l`/`yt_r` and drawn by
  `plot_solution`. Validated: for `f(x)=1` it reproduces Spencer's thrust line to
  ~1e-6–1e-11 across left/right-facing, piles, reinforcement, submerged, and
  non-circular.

**S7 — Documentation.**
- *Deliverable:*
  - new **`docs/lem/mprice.md`** — derivation in the Spencer style, reusing the §4a
    moment table and the `F`-vs-`λ` figure (P2 drafts the prose during S1–S5; final
    numbers/figure dropped in here);
  - **`mkdocs.yml` nav entry** under the LEM section, immediately **after "Spencer's
    Method"**: `- Morgenstern–Price Method: lem/mprice.md`;
  - link the page from `docs/lem/overview.md`;
  - add the M-P row(s) to `docs/verification.md` (P3 scaffolds);
  - add a sample + a `run_tests.py` regression tag.
- *Depends-on:* S5 (numbers/figure), and P2/P3 drafts (§11). *Exit:* `mkdocs build`
  clean **and the page is reachable from the nav**; the new regression tag passes in
  `run_tests.py`; the verification table matches the benchmark run. Per the "docs
  track solver" rule, S6 + S7 land together as one work unit.

## 11. Parallelization & subagents

The core (S1–S5) is a single-threaded chain on `solve.py` and must **not** be
fanned out — the stages depend on each other and edit the same code. But four
streams are genuinely independent and can run in parallel with, or ahead of, the
core, plus the per-gate verification can be delegated.

**Independent streams (spawn up front, run alongside the core):**

| # | Subagent task | Type | Depends only on | Runs during |
|---|---|---|---|---|
| P1 | **Recon** (S0): map `force_equilibrium` contract + callers; verify §4a columns/signs | `Explore` (read-only) | current code | before S1 |
| P2 | **Docs draft** `docs/lem/mprice.md` from this plan + `spencer.md` as style template (derivation, moment-table prose, Spencer/M-P comparison) | general/docs | the plan | all of S1–S5; figures/numbers + nav entry land at S7 |
| P3 | **Benchmark + verification scaffolding**: regression tags + `verification.md` M-P table/acceptance from the §7.2 spec | general | the **fixed result contract** (locked) | S2 onward |
| P4 | **λ-source research**: find a published case reporting λ for a stated `f(x)` (Fredlund–Krahn 1977 / SLOPE-W manual; §8.1) | research (web) | nothing | anytime |

**Delegatable verification gates (sequential, but offloaded so the main thread keeps moving):**

- **After S1:** a `verify` agent runs the LEM benchmark suite pre/post refactor and reports Corps/L-K bit-stability (the S1 exit gate).
- **After S3 and S5:** a `verify` agent runs the `f=1 ≡ Spencer` comparison (S3 gate) and the A↔B agreement across all benchmarks (S5 exit).

**Orchestration shape:** fan out P1–P4 as independent agents; build S1→S5 sequentially in the main thread; drop a verify agent at the S1, S3, and S5 gates. P2/P3's drafts land into S6 with only figures/final numbers to fill once the solver is green. This parallelizes essentially all non-core work while keeping the engine a clean sequential chain.
