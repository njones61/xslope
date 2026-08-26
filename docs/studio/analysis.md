# Running Analyses

Studio runs the same LEM, seepage, and FEM analyses as the `xslope` library, each
behind a run-options dialog. Long solves run on a background thread, so the window
stays responsive, with live output in the [Log pane](interface.md#the-log-pane)
and a Cancel button.

Pick the analysis type on the toolbar's **Mode** switch — **LEM**, **Seepage** or
**FEM**, one click or `Ctrl+1` / `Ctrl+2` / `Ctrl+3` — then click **Run** (its label
follows the mode: **Run LEM…**, **Run Seep…**, **Run FEM…**).

---

## Model checks before a run

Every run dialog checks the model against what the analysis it is about to start
actually needs, and shows the result in a **Model checks** column beside the run
controls. The point is the answer that looks fine: a blank pore-pressure ratio, a
hydraulic conductivity of zero, a material with no tensile cap, a boundary set that
drives no flow — each of those runs to completion and returns a number, and the
number is wrong. The checks say so before the solve rather than after the report.

The column is a **list**: one line per finding, marked with its severity, with the
full text of the selected line underneath it. The dialog opens on the first error,
or on the first line when there is none, so what is on screen is always the thing
most worth reading. A check that fires once per material — four materials with no
tensile cap, say — takes **one line, not four**: the line says how many and which
(*"… — 4 materials: 1 ('Shell'), 2 ('Core'), 3 ('Clay'), 4 ('Sand')"*), and the
detail states the explanation they share once, with each material's own numbers
under it. Anything the checks offer to fix is offered once, on that one entry.

Findings come in three severities, and they behave differently:

- An **error** means the run would crash, or would produce a provably wrong
  answer. **Run** is disabled while one stands, and the message names the sheet
  and the field it is about.
- A **warning** means the run will proceed and the answer may well be fine, but
  the model matches a pattern that has produced wrong answers before. A warning
  **never blocks a run**. It is on screen beside the button so it informs the
  decision instead of annotating the regret.
- A **note** records a default that was applied or an input that is inert. Notes
  are collapsed behind a *n notes* line, available without being in the way.

Here the circles sheet is empty, so there is no surface to analyse — an error, and
**Run** is greyed out until it is resolved:

![Run LEM with an error and a remedy button](images/analysis_run_lem_preflight.png)

The button under the finding's text is a **remedy**: a fix the check can offer for
what it found. A remedy is always offered and never applied by itself, because a change you
did not ask for is the same disease in a helpful disguise — the model would quietly
stop matching the file you typed. Pressing the button first shows exactly what the
change would be (*"Add 6 starting circles to the circles sheet: …"*); applying it
lands on the normal undo stack, is saved like any other edit, and re-runs the checks
so the list reflects the model that would now run. A remedy that cannot fully
succeed on this model — a piezometric line whose x values rise and then fall cannot
simply be reversed — shows its button dimmed, with the reason as its tooltip, rather
than failing when pressed.

The checks offer these fixes today: reverse a load or piezometric line entered
right-to-left, generate a starting set of circles from the slope geometry, add
standing water as a distributed load (manual water-load models), and switch a model
to [automatic water loads](../usage/input_template.md#worksheet-main), which is the better fix of the
last two — a written-in load is a snapshot that goes stale the moment the pool
moves, and a derived one cannot.

The same checks drive what the dialog lets you pick. OMS and Bishop take moments
about a circle center, so on a non-circular surface they are dimmed, with the reason
on the item itself; change the **Surface** selector and the method list re-filters on
the spot:

![Method list with OMS and Bishop dimmed](images/analysis_run_lem_methods.png)

A model that defines both a circular and a non-circular surface has no way, in the
file's geometry, of saying which one it means, and the circles simply win with no
message. The **Surface** selector is that choice, and it is written back to the model
when you run — into the **main** sheet's **Surface family** cell, so it survives saving
and reopening the file — and the next run, the plots and the saved file all read the
same answer.

---

## Limit equilibrium (LEM)

In **LEM** mode, **Run LEM…** opens a dialog with:

![Run LEM dialog](images/analysis_run_lem_dialog.png)

- **Method** — OMS, Bishop, Janbu, Corps of Engineers, Lowe & Karafiath, Spencer,
  or Morgenstern–Price.
- **Analysis** — single surface or automated search. (Probabilistic
  **reliability** analysis has its own toolbar button — see
  [Reliability analysis](#reliability-analysis) below.)
- **Surface** — circular or non-circular. The selector appears when the file
  defines both, opens on the family the file states (the **main** sheet's
  **Surface family** cell, which the last run's choice was written into), and
  re-filters the **Method** list as you change it (see
  [Model checks](#model-checks-before-a-run)).
- **Number of slices**, the **rapid drawdown** flag, and a **diagnostic** toggle.
- **Composite surfaces** — lets a circle deeper than the bottom of the model be
  truncated at it and run along the base between the two crossings (see
  [Composite Failure Surfaces](../lem/overview.md#composite-failure-surfaces)).
  Off by default, and only available for circular surfaces. Turn it on when the
  base of the model is a real impenetrable boundary — bedrock, or a weak seam
  resting on it — because the critical mechanism there follows the base and no
  ordinary circle can reach it. Leave it off when the bottom of the model is
  simply how deep you chose to look.
- **Grid search** — seeds the circular search from an automatic
  grid-and-tangent sweep instead of (only) the circles sheet (see
  [Grid Seeding](../lem/search.md#grid-seeding-global-search)). Off by default,
  available for circular auto-search. Turn it on to protect
  against the local-minimum trap of a single starting circle — a seed in the
  wrong family can converge 20% or more too high with no warning — or when you
  have no idea where the critical circle is (the circles sheet may even be
  empty). It reports the most critical surface anywhere in the model; leave it
  off to interrogate a specific mechanism with your own circles.
- **Ignore surficial (skin) failures** — rejects trial surfaces shallower than the
  **Min slip depth** below the ground surface. Off by default, and available for the
  auto-search. On a cohesionless face the critical surface is an infinitely shallow
  skin slide, and without the filter the search chases it instead of the deep-seated
  mechanism a design wants; see [surficial failures and the minimum-slip-depth
  filter](../fem/overview.md#surficial-skin-failures-and-the-minimum-slip-depth-filter)
  for how to choose the depth.
- **Search tolerances** (`fs_tol`, `tol`, `max_iter`) — enabled for the
  search-driven analyses.
- **Seepage time** and the **rapid-drawdown stage times** — shown when the model
  carries transient seepage. See [Seepage time](#seepage-time).

The result depends on the analysis type:

| Analysis | Result tabs |
| --- | --- |
| **Single surface** | **LEM · Solution** — the surface with slices, base stresses, and thrust line. |
| **Automated search** | **LEM · Search** (all trial surfaces + critical + search path) *and* **LEM · Solution** (the critical surface). |

![LEM Search view](images/analysis_lem_search.png)

![LEM Solution view](images/analysis_lem_solution.png)

When the accepted solution carries admissibility defects — base tension on a
cohesionless slice, interslice tension, or a line of thrust that leaves the
slices — an amber strip across the top of the **LEM · Solution** view lists
them, so an inadmissible FS is never mistaken for a clean success:

![LEM Solution with admissibility warnings](images/analysis_solution_warnings.png)

The warnings never change the factor of safety — they say the internal force
distribution behind it is strained, and not all of them are equally alarming —
see [Interpreting the Admissibility
Warnings](../lem/spencer.md#interpreting-the-admissibility-warnings). A clean
solution shows no strip.

Both circular and non-circular surfaces are supported for single solves and
searches. Search iteration progress streams to the Log pane, and a search can be
**cancelled** from the status bar.

---

## Parametric study

A second entry point sits beside Run — **Parametric…** — on both the **Run menu**
and the main **toolbar**, in all three modes (LEM, Seepage, FEM). It is available whenever
a model is open — in Seepage and FEM mode it additionally needs a built mesh, exactly like
Run — and is disabled while another analysis is running. Every version of it drives the
same sweep engine as the library — see
[Parametric Studies](../parametric/index.md) for the engine, and the
[`/xslope` skill](../usage/claude/index.md) for the scripted recipes.

One dialog covers four study modes, chosen by its **Mode** selector — **Sensitivity**,
**Design**, **Back-Analysis**, and **Factor of safety vs time**. In **LEM** three controls
are shared by all of them:
**Method** (any of the seven LEM methods), **Number of slices**, and a **Parameter** picker —
a **Material** dropdown (each material plus a *k_seismic (global)* entry) and a **Property**
dropdown listing that material's option-aware sweepable fields (both drawn from the engine's
`list_params`). A **Re-search the critical surface at each step** checkbox applies throughout:
on by default — the honest setting, since the critical surface moves as the parameter changes
— and off re-solves the entered surface only (much faster, but right only for that prescribed
surface).

**Sensitivity** sweeps several parameters and visualizes how FS responds. A **Plot type**
selector chooses the view — a **tornado** (the default), **scaled-sensitivity bars** (with a
**Scaling** sub-choice: elasticity, per-1%, or per-σ), a **spider** plot, and — only when the
model carries reliability standard deviations — a **variance Pareto** and a **Monte Carlo rank
correlation** (with an **MC samples** count). The tornado, scaled, and spider plots sweep the
parameters listed in the table; the variance and rank plots instead use every σ-carrying
material automatically. The plots themselves are shown, with worked examples, on the
[Sensitivity](../parametric/sensitivity.md#sensitivity-plots) engine page.

![Sensitivity dialog](images/analysis_sensitivity_dialog.png)

- A **Default ±%** and a **Points** count (points per parameter's FS-vs-value curve; the
  tornado uses only the curve's two endpoints, the spider draws the whole curve).
- **Add parameter** appends the currently picked material/property to the table. Each row
  shows the parameter reference, an editable **±%** overriding the default for that row, a
  **σ** button, and a remove (✕) button.
- The **σ** preset swaps that row's ±% range for a ±one-standard-deviation range built from
  the model's reliability `sigma_*` columns — the same standard deviations the
  [reliability analysis](../reliability/index.md) uses — so a sweep can mirror a reliability
  input band with one click. The button is disabled for a property that carries no `sigma_*`.

**Design** sweeps the one picked parameter toward a target FS:

![Design dialog](images/analysis_sensitivity_dialog_design.png)

- **From** / **To** bound the swept value (seeded to ±50% of the current value the first time
  you pick a property), **Steps** sets the number of solves, and **Target FS** is the factor
  of safety to locate.

![Back-Analysis dialog](images/analysis_sensitivity_dialog_backanalysis.png)

**Back-Analysis** is the same single-parameter sweep as Design, framed for a failure
investigation: because a slide has occurred, the target defaults to **FS = 1.0**, and the
result is read as the parameter value consistent with the observed failure (the
back-calculated strength, most commonly). The controls are identical to Design.

### Factor of safety vs time

The fourth mode sweeps the saved instants of a
[transient seepage](../seep/transient.md) march instead of a parameter. Each point solves
the same model against that instant's pore pressures — no input changes between them, so
the axis is time — and the reservoir load is re-derived from the pool as it stood at that
moment. The parameter picker steps aside and a **Saved frames** checklist takes its place,
listing every instant the march stored with all of them ticked; **All** and **None** set the
whole list, and unticking samples a long march, each frame being a full stability run. The
**Method**, **Number of slices** and **Re-search the critical surface at each step**
controls apply as they do everywhere else, and the circles sheet's search window is
applied, which is what keeps the curve on one mechanism rather than letting it jump
families between instants:

![Parametric dialog in Factor of safety vs time mode](images/analysis_sensitivity_fs_time_dialog.png)

A **Rapid drawdown at each time** checkbox below the frames list turns every ticked instant
into a three-stage [rapid drawdown](../lem/rapid.md) instead of a single-stage
analysis: stage 1 is the march's initial state (the `tseep` sheet's `stage_1`, normally
t = 0 at full pool), stage 2 is that instant's drawn-down state, and stage 3 re-checks the
same section with drained strengths. The reported value is the drawdown's own — the lower of
stages 2 and 3 — so the curve answers *how safe is this slope if the pool falls to where it
stands at t*, instant by instant. The box is available only on a model whose materials carry
the drawdown strengths `d` and ψ; without them it is greyed and says so, because all three
stages would read the same strengths. A drawdown point is always searched from the starting
circle, so **Re-search** is held on and greyed while the box is ticked.

The mode is available in **LEM** and **FEM** mode on a model that carries a transient
seepage solution and at least one material reading `u = seep`. Without one of those it is
disabled and names the reason — *Run a transient seepage analysis first*, *No material takes
its pore pressure from the seepage solution*, or, in Seepage mode, that the seepage solution
is this run's input rather than its output. The drawdown option is LEM only: the three-stage
procedure is a limit-equilibrium construction with no SSRM equivalent. The engine page
[Factor of safety versus time](../parametric/sensitivity.md#factor-of-safety-versus-time)
carries the sweep itself.

### Running and cancelling

Clicking **Run** launches the sweep on a background thread, so the window stays responsive.
The progress bar tracks each solve (`done/total`, with the current swept value echoed in the
status bar), the sweep log streams to the [Log pane](interface.md#the-log-pane), and the
status-bar **Cancel** button aborts cooperatively — the in-flight solve finishes, then the
sweep stops and the app is left consistent, with no partial result stored.

### Results

**Sensitivity** opens a **Sensitivity** tab with the selected plot. For the default
**tornado** — one horizontal bar per parameter, widest on top, with the base-case FS drawn as
a labelled vertical reference line:

![Sensitivity tornado](images/analysis_sensitivity_tornado.png)

On the tornado, **double-click a bar** to open a companion **Sensitivity · Curve** tab
showing that parameter's full FS-vs-value curve (the base case marked, and any
critical-surface jump drawn as an open circle):

![Per-parameter FS curve](images/analysis_sensitivity_curve.png)

The other plot types render into the same **Sensitivity** tab — the scaled-sensitivity bars,
the spider plot, the variance-contribution Pareto, and the Monte Carlo rank-correlation bars
(the [Sensitivity](../parametric/sensitivity.md#sensitivity-plots) page shows each with a
worked example). The double-click click-through is a tornado affordance; the other plots do
not offer it.

**Design** opens a **Design** tab with the FS-vs-value curve, the FS = 1 and target-FS guide
lines, and — when the target is bracketed — the interpolated crossing marked with a green
diamond and annotated *property = value for FS = target* (a **Back-Analysis** run renders the
same tab, with the target at FS = 1.0):

![Design curve with crossing](images/analysis_sensitivity_design_curve.png)

**Factor of safety vs time** opens an **FS vs Time** tab: the factor of safety at each
evaluated instant as a line with a marker per point, the lowest of them ringed and
annotated with its own time, and the model's `tseep` time series — the drawdown schedule
that drives the curve — drawn faintly behind it on a second axis. When the march is done the
per-instant table (time, factor of safety, the critical circle) is printed to the
[Log pane](interface.md#the-log-pane), and an instant that produced no result appears there
with its reason rather than as a gap in the line. A **rapid drawdown** run opens the same tab
under the name **Drawdown vs Time**: its table carries the three stage factors of safety and
which of stages 2 and 3 governed each instant, and the plot draws those stages as thin dashed
lines behind the reported curve:

![FS vs Time result tab](images/analysis_sensitivity_fs_time.png)

When the swept range never reaches the target, the result is honest about it: no crossing is
drawn, and an amber note reports the FS span and which way to widen the range — the GUI face
of the engine's
[never-extrapolate discipline](../parametric/design.md#honest-about-misses):

![Design honest miss](images/analysis_sensitivity_design_miss.png)

### Sweeps in FEM and Seepage mode

The Parametric dialog has a version for every mode. In **LEM** it sweeps the
limit-equilibrium analyses, and the swept output is the factor of safety. In **FEM**
each swept point is a full SSRM solve (output: factor of safety) — expect minutes per
step, so it runs in the background and is cancellable. In **Seepage** the output is
the **total discharge q** through the section. FEM and seepage sweeps run on the mesh,
so build one first. The variance Pareto and Monte Carlo rank plots reuse the LEM
Taylor-series and Monte-Carlo reliability, so they are offered only for an LEM study
that carries sigmas.

The dialog's solver rows follow the app mode: in **FEM** they become the SSRM knobs
(`F_min` / `F_max`, tolerance, failure criterion), and in **Seepage** they become the BC
set and convergence tolerance, with the design target a discharge **q** rather than an FS:

![Parametric dialog in FEM mode](images/analysis_sensitivity_dialog_fem.png)

![Parametric dialog in Seepage mode](images/analysis_sensitivity_dialog_seep.png)

---

## Reliability analysis

A **Reliability…** button sits beside **Parametric…** on the **Run menu** and the
**toolbar** — its probabilistic sibling. Where the Parametric study answers
deterministic what-ifs, Reliability turns the material standard deviations (the
`s(·)` columns of the mat sheet) into a reliability index **β** and a **probability
of failure**. It is available in **LEM** and **FEM** modes (not Seepage); the FEM run
needs a built mesh, like Run.

The dialog offers a **Method** selector with three engines:

![Reliability dialog (LEM)](images/analysis_reliability_dialog_lem.png)

- **Taylor series (TSPM)** — the mean-value factor of safety plus a ±σ perturbation
  of each uncertain parameter (`1 + 2N` solves). Available in both modes; in **FEM**
  each factor of safety comes from an SSRM solve.
- **Monte Carlo** — samples every uncertain parameter and evaluates the factor of
  safety of each realization on a fixed surface, reported as an FS histogram. Monte
  Carlo needs ~10⁴ solves, which is affordable with a limit-equilibrium solve but not
  with the finite-element SSRM, so it is **disabled in FEM mode** (a one-line note
  explains why, and FEM reliability stays on the Taylor series):
- **Response surface (RS)** — the same sampling with the factor of safety taken from
  a quadratic surrogate, fitted to a few dozen real solves and measured against
  held-out ones before it is used. Ten million realizations at no sampling cost, so
  a small probability of failure is resolved; a surrogate that fails its accuracy
  gate refuses the run rather than answering (the [Monte Carlo
  page](../reliability/monte_carlo.md#sampling-a-fitted-response-surface) gives the
  design, the gate and the measured accuracy). Limit-equilibrium only, for the same
  reason as Monte Carlo:

![Reliability dialog (FEM)](images/analysis_reliability_dialog_fem.png)

Below the engine controls, a read-only **Standard deviations in this file** summary
lists every `s(·)` column with its value, σ, and COV, so you can confirm what the run
will vary. **LEM** adds the solver method, surface, slice count, rapid-drawdown flag,
and a *search the critical surface at the mean values* toggle; **Monte Carlo** adds an
**MC samples** count, a **seed** (fixed by default, so the result is reproducible),
a normal / lognormal **distribution** choice, an **MC sampling** choice (random,
or Latin hypercube — measured at roughly 3× the information per realization on
the sample models), and an optional convergence stop —
**Stop when P_f converges** with a **P_f tolerance (±)** as a percentage *of P_f
itself* ends the campaign once the 95% confidence half-width on the empirical
probability of failure is inside that fraction, with the samples count as the cap
(the [Monte Carlo page](../reliability/monte_carlo.md) gives the rule, its
rare-event guard, and the convergence plot). The **Response surface** engine shares
the seed and distribution controls and ignores the sample count and the convergence
stop — its realization count is fixed at ten million, which costs seconds of
arithmetic rather than solves, and its accuracy is set by the fit rather than by the
count. **FEM** shows the SSRM `F_min` / `F_max` bracket and a tight reliability
tolerance. Settings are remembered for the session.

The run reports β, the probability of failure, and the mean / σ of the factor of
safety in a summary, with the per-parameter table in the Log pane. The result view
follows the engine:

| Engine | Result tab |
| --- | --- |
| **Taylor series (LEM)** | **LEM · Reliability** — the most-likely-value surface with the F⁺/F⁻ perturbation surfaces. |
| **Monte Carlo (LEM)** | **Reliability · MC** — the FS histogram with the FS = 1 line, the mean, and fitted normal / lognormal overlays (a display-panel toggle). β in both conventions and the probability of failure are in the title. |
| **Response surface (LEM)** | **Reliability · RS** — the same FS histogram, drawn from a 10,000-realization subsample of the surrogate draws. The fit and gate summary — design solves, $R^2$, RMS error, the failure-region check — is in the Log pane and the run summary. |
| **Taylor series (FEM)** | **FEM · Results** — the deformation at the most-likely values; β and the probability of failure are in the run summary. |

![LEM Reliability view](images/analysis_lem_reliability.png)

![Monte Carlo reliability histogram](images/analysis_reliability_mc_histogram.png)

The engines are the same ones the library exposes — `reliability` (the front door),
`reliability_taylor`, `reliability_mc`, `reliability_rs`, and `reliability_fem`; see
[Reliability Analysis](../reliability/index.md) for the theory and worked examples.

---

## Building a mesh

Seepage and FEM run on a finite-element mesh, which you build explicitly. In
**Seepage** or **FEM** mode, **Build Mesh** opens a dialog with:

![Build Mesh dialog](images/analysis_build_mesh_dialog.png)

- **Element type** — `tri3`, `tri6`, `quad4`, `quad8`, or `quad9`.
- **Target size** — entered directly, or auto-sized as the slope width divided by a
  number of divisions.
- **1D element size** — the element size along the reinforcement and pile lines,
  blank by default, which meshes them at the target size like everything else. A
  value refines the beam and bar elements *and* the soil elements that share their
  nodes, growing back to the target size away from the lines, so a member can be
  discretized finely without paying for a finer mesh across the whole section. It is
  a model input — **1D element size** on the main sheet — so the box opens on
  whatever the file states and an entry is written back to the model, saved with the
  file and undone like any other edit. A value at or above the target size is
  ignored, since sizes compose by taking the smaller and a coarser request could
  never bind.
- **Quadrilateral style** — **Free (recommended)**, or **Structured where possible**,
  which sweeps rows and columns of quads through the grid-like zones and free-meshes
  the rest. Available for the quad element types only, and dimmed for the triangular
  ones. It is a per-run choice and is not saved with the model; see
  [Quadrilateral Meshing Styles](../fem/mesh.md#quadrilateral-meshing-styles) for what
  each style produces.
- **Refine near features** — off by default. When checked, elements shrink near model
  features (reinforcement/pile lines, crack tips, thin material zones) and grow back to
  the target size away from them; the **Refinement factor** spinbox (default 3.0) sets
  the local size to *target size ÷ factor*. Leaving it off meshes at the target size
  everywhere. Refinement is detected automatically from the geometry — there is nothing
  to place by hand. (Selecting individual feature classes is available in the Python API
  via `refine_features`; the dialog refines near all of them.)
- **Refine thin zones** — on by default. A material zone too thin for the mesh to
  resolve does not fail: it solves, and returns a factor of safety that is too high,
  because a shear band cannot form across a single element. Checked, every thin zone
  is meshed with about four element rows across its local width, down to a limit of six
  times the target size — a very thin band in a large section would otherwise multiply
  the node count from a box nobody ticked on this run. A zone the limit leaves under
  three element rows is named in the model checks, which measure the mesh that was
  built; giving it the full size it asks for is then a **Size** on the zone, which is
  not capped, or a finer target size. Thickness is measured on the **material**: a layer
  stored as several polygons because the ground surface steps through it is measured as
  the one layer it is, so a layer the global size already resolves is left alone. A
  section with no thin zone meshes exactly as it would with the box clear, and the Log
  names every zone that was refined, the local size it got and the rows that buys. The
  **Refinement factor** above does not enter — a thin zone is sized by its own
  thickness. See [Thin material zones](../fem/mesh.md#thin-material-zones).

The mesh is built on a background thread (it includes reinforcement and pile
constraint lines, so it serves FEM too), shown in a **Mesh** tab, and written to a
`{stem}_mesh.json` sidecar. Seep/FEM **Run** stays disabled until a mesh exists; a
geometry edit that invalidates the mesh clears it and re-gates Run.

![Mesh view](images/analysis_mesh_view.png)

Mesh generation uses **gmsh**, which is installed by the `fem` extra
(`pip install "xslope[gui,fem]"`) — see [Installation](index.md#installation).

---

## Seepage

In **Seepage** mode, **Run Seep…** opens a dialog with the solve parameters —
the **convergence tolerance** and the iteration ceiling, **Max iterations** —
above the
[model checks](#model-checks-before-a-run). A steady run solves every boundary
set the file defines: set 1 alone for most files, and both sets for a
rapid-drawdown pair (a caption in the dialog says so when a second set is
present). Display choices — the
plotted variable, contour levels, flow lines, vectors, fill, the phreatic
surface — are not run options; they live on the
[Display panel](#display-options-per-view) of the solution view and re-render
the cached solution without re-solving.

![Run Seep dialog](images/analysis_run_seep_dialog.png)

The run produces two tabs — **Seep · Data** (mesh + boundary conditions) and
**Seep · Solution** (the chosen variable, contours, phreatic surface, flow lines /
vectors) — and writes the solution to `{stem}_seep.csv` (`_seep2.csv` for BC 2).
The convergence trace streams to the Log pane. Each BC set keeps its own tab pair,
so rapid-drawdown problems show BC 1 and BC 2 together.

![Seep Data view](images/analysis_seep_data.png)

![Seep Solution view](images/analysis_seep_solution.png)

---

## Transient seepage

When the open file carries a **tseep** sheet (a transient-seepage time series —
see the [Input Template](../usage/input_template.md)), the **Run Seep…** dialog
grows a **Run type** selector with a **Transient (time-dependent)** choice
alongside **Steady**. A transient solve marches the variably-saturated flow
equation through time and saves a sequence of frames, each a full seepage
solution at one instant. The theory — the storage formulation, time-stepping, the
submerged-only boundary rule, and the initial condition — is on the
[Transient Seepage](../seep/transient.md) page, whose
[Studio](../seep/transient.md#studio) section walks the same workflow from the
theory side.

A transient run drives the flow with BC set 1's time series. The run's other
parameters — the duration, save schedule,
rapid-drawdown stage times, and the time series themselves — are model inputs, edited
under **Inputs → Transient** (the [Transient editor](editing.md#transient-seepage)),
so the dialog carries no stage fields, only a caption pointing there.

![Run Seep dialog in Transient mode](images/analysis_run_seep_transient.png)

Because the march covers a known duration, the run reports **determinate progress** —
the status bar's progress bar tracks the simulated-time fraction (`t / duration`) —
and a **Cancel** button beside it stops the march cleanly; a cancelled run stores no
partial result and returns the UI to idle.

The run produces a single **Seep · Transient** tab that shows one frame at a time
through the *same* solution renderer as the steady **Seep · Solution** view — so
contours, flow lines, velocity vectors, the phreatic surface, and the colorbar all
look identical. Each frame's title reads **Seepage Solution — t = … ** with a smaller
second line carrying the frame's boundary **Inflow / Outflow** (they differ under
transient storage exchange). A pure storage-release frame — a drawdown instant with no
through-flow — has no flow net, so its flow lines are omitted and the second line
reads *no through-flow — flow lines undefined* instead. Along the bottom of the plot
sits a **play bar**:

- **Transport buttons** — first, previous, play/pause, next, and last frame.
- A **frame slider** to scrub the saved frames.
- A **t =** readout that shows the current frame time and doubles as a
  **jump-to-time** entry — type a time, press Enter, and the view snaps to the
  nearest saved frame.
- A **Speed** selector (0.5× – 4×) for playback, which advances the frames in order.

![Seep · Transient view with the play bar](images/analysis_seep_transient_playbar.png)

The frame bundle is written next to the model as a `{stem}_tseep.csv` (plus a
`{stem}_tseep_meta.json` ledger) and restored into the Seep · Transient tab on the
next Open, so a saved transient run reloads without re-solving. Its
[Display panel](#display-options-per-view) is the seep-solution panel tailored for a
frame series: **Water levels** default on (the pool visibly drops as it plays), and
the flow-net-only **Flow lines** and **Base material** controls are omitted — a
transient state has no flow net (see the note below). Changing an option re-renders
the shown frame.

An **LEM** or **FEM** run consumes ONE frame: with `u = seep`, the instant named in
the Run dialog's **Seepage time** group supplies the pore pressures. A
**rapid-drawdown** LEM run instead stages the `stage_1` and `stage_2` frames into the
drawdown analysis. See [Rapid Drawdown from a Transient
Solution](../lem/rapid.md#transient-solution) for how the staged frames enter the
three-stage calculation.

### Seepage time {#seepage-time}

A stability run against a transient seepage solution reports a factor of safety for a
single instant, so the **Run LEM** and **Run FEM** dialogs carry a **Seepage time**
group whenever a transient solution is loaded. It offers three ways to name the
instant:

- **Saved frame** — a dropdown of the times this solution actually saved. Instant:
  the pore pressures already exist.
- **Frame shown in the results viewer** — the time the Seep · Transient play bar is
  displaying. Also instant, and offered only while that tab is open.
- **Another time** — any instant within the run duration. The pore pressures for it
  do not exist yet and are never interpolated between frames, so choosing it
  **re-marches** the transient solution with that instant added to the save schedule,
  then starts the analysis. The dialog says so before you commit, the re-march reports
  progress, and it can be cancelled — cancelling it abandons the analysis too.

The group opens on the model's own `stability_time` when the tseep sheet declares one,
and otherwise on the **last saved frame** — usually the drained end state, which is
what a blank `stability_time` means. The note under the controls states which instant
will be read either way.

The choice governs that run only. Tick **Save as the model's stability time** to write
it into the tseep sheet as well, which is what makes a scripted or headless re-run read
the same frame.

On a model with a tseep sheet but no solution loaded the group is present but
disabled, carrying the reason; a steady model has one pore-pressure field and no
group at all.

The [model checks](#model-checks-before-a-run) read the group. A material with
`u = seep` needs a solved pore-pressure field, and a transient model does not carry
one in the file — the chosen frame is written into the model when the run starts. So
while a solution is loaded, the checks take the **Seepage time** choice as the answer
and note which instant it is (*Pore pressures come from the transient seepage
solution, at t = …*) rather than refusing a run for a field that is one step away.
A **rapid drawdown** names both stage times the same way, and those two frames are
the stage-1 / stage-2 pair the drawdown check asks for — one instant is not, because
it supplies stage 1 only. With no solution loaded there is nothing to stage, and the
checks refuse as they would on any other model: run the seepage analysis first.

### Rapid-drawdown stage times {#stage-times}

Ticking **Rapid drawdown** in the Run LEM dialog replaces the single-instant selector
with **Stage 1 time** and **Stage 2 time**, pre-filled from the model. These are the
two instants the drawdown stages read out of the march — pure extraction parameters,
since the drawdown schedule itself lives in the boundary conditions — so they are
edited here, at the point of use, as well as under
[Inputs → Transient](editing.md#transient-seepage). Both places write the same two
values on the tseep sheet, and an edit here lands in the model immediately.

Run refuses stage times it cannot use — one blank, stage 2 at or before stage 1, or a
stage beyond the run duration — and says which. Stage times the loaded solution never
saved trigger the same re-march as a free-entry seepage time.

---

## Finite element (FEM)

In **FEM** mode, **Run FEM…** offers a **single trial** or an **SSRM** run (the
Shear Strength Reduction Method), with `F` (or `F_min`/`F_max`), a tolerance, and
the failure criterion.

![Run FEM dialog](images/analysis_run_fem_dialog.png)

**Max iterations per trial** (12000) is the viscoplastic budget each trial *starts*
with, not a hard stop. A trial that uses it up with its out-of-balance force still
falling is given another budget's worth, repeatedly, until either it settles or it
reaches the **Iteration ceiling** (50000). That is what keeps the answer from
depending on the budget: a model needing 16,000 iterations returns the same factor
of safety whether the budget was typed as 3000 or 12000. A trial that reaches the
ceiling while still improving is *inconclusive* — neither settled nor failed — and
the run says so in the Log rather than counting it as a failure. The factor of safety
is still the final bracket's midpoint, as on any other run; what changes is that the
bracket's upper edge is an undecided trial rather than a measured failure, which the
Log states beside the answer. Raise the ceiling, or loosen the SSRM tolerance, when
that happens.

The [model checks](#model-checks-before-a-run) in the dialog's second column are the
finite-element ones: a blank Poisson's ratio (which reads as 0.0 and moved the
strength-reduction factor of safety by a third on the reference model), a modulus of
zero, a mesh that references a material the table does not define, and — above — a
material with no tensile cap, which grants it unbounded tension and raises the factor
of safety with nothing else on screen to show it.

When the seismic coefficient is nonzero the checks also carry a note about what its
**sign** means here, because it does not mean the same thing in both engines. The
finite-element engine reads `main!D13` as a vector: `+k` pushes in `+x` and `−k` in
`−x`, and since both faces of an embankment are analysed at once, choosing the
direction is a modelling decision the engine will not make for you — a pseudo-static
factor of safety can legitimately come out *above* the static one for the face the
shaking stabilises. The limit-equilibrium engine reads the same cell as a magnitude
and orients it itself, and its own dialog says so.

**Side BC** chooses what holds the left and right edges of the model. **Rollers** —
the default, and what every model that does not say otherwise means — fixes the
horizontal component and leaves the vertical free, so the truncated ground can still
settle under its own weight. **Fixed** clamps both components, which is what RS2 does
on its side boundaries. Fixed is a vendor-parity option rather than a better model: it
adds shear restraint the real ground does not have, and stiffens a domain truncated
close to the slope. The setting is part of how the model is restrained rather than
part of the strength reduction, so it applies to a single trial and an SSRM alike. Like
the other run options it is seeded from the open file (`main!D22`) and remembered for
the session, but a choice made here is not written back into the file.

When the model carries a transient seepage solution the dialog also carries the
**Seepage time** group, which names the instant the pore pressures are read from —
the same control the Run LEM dialog uses, described under
[Seepage time](#seepage-time).

For an SSRM run, the **SSR exclusions…** button opens a checkbox picker — one row per
material zone in the model, checked (included) by default:

![SSR exclusions dialog](images/analysis_ssr_exclude_dialog.png)

Unchecking a zone holds it at full strength through every trial factor instead of
dividing its c and tan(φ) like the rest of the model — RS2's per-material *Apply_SSR*
flag / SSR Exclusion Area. The mechanism is pushed up and out of an excluded zone, which
is useful for keeping a non-participating zone (a stiff foundation, say) from carrying
the failure, and for reproducing a vendor analysis that constrains the mechanism the same
way — see [SSR Exclusion Zones](../fem/overview.md#ssr-exclusion-zones) for the
engineering rationale and a worked comparison against RS2. The button and the summary
label next to it are gated to the SSRM analysis; the choice is a run option, not a model
property, so it lives with the rest of the dialog's settings (remembered for the session
to prefill the next run) rather than being saved into the input file.

An SSRM run also offers **failure-state capture**. Once the factor-of-safety
bracket resolves, **Capture failure-state mechanism** (on by default) re-solves
once just beyond the critical factor with the displacement cap off, so the
unconverged field develops the actual collapse *mechanism* the deformation and
displacement-vector figures render — rather than the diffuse settlement of the last
converged trial. **Capture margin** sets how far beyond critical that snapshot is
solved (a fraction of FS), and an optional **capture iteration budget** overrides
the automatic ceiling. Turning capture off skips the extra solve; the factor of
safety and the bracket are unaffected either way. The controls are gated to the
SSRM analysis, since a single trial has no bracket to capture beyond.

The run produces **FEM · Data** (mesh + boundary conditions + reinforcement) and
**FEM · Results** (deformation, shear strain, displacement vectors). An SSRM run
reports the factor of safety and can be **cancelled** mid-run. The solution is
exported alongside the model so it can be restored on the next Open without
re-solving — including the at-failure mechanism snapshot (a second CSV pair) and,
for a model with reinforcement or piles, the per-element structural results
(`{stem}_fem_reinf.csv` / `{stem}_fem_piles.csv`), so a reloaded solution redraws
the reinforcement-force and pile-shear plots without solving again.

![FEM Data view](images/analysis_fem_data.png)

![FEM Results view](images/analysis_fem_results.png)

The FEM · Results [Display panel](#display-options-per-view) carries the controls
that shape these plots. When SSRM captured an at-failure mechanism, a **Field
state** switch chooses which field every panel renders — **At failure** (the
default: the developed collapse mechanism) or **Last converged** (the sub-critical
converged solution) — and it applies to the deformation, shear-strain, and
displacement-vector plots alike, so they always tell the same story. The
deformation plot adds its own controls: the **Original mesh** reference (dashed
outline, full grid, or off), the **Deformed color** of the displaced grid, and
the exaggeration pair — **Scale ×**, the displacement multiplier the plot title
prints, whose **Auto** default picks whatever draws the largest displacement at
the **Auto size** percent of the mesh height; entering an explicit Scale ×
dims Auto size until the box returns to Auto. The displacement-vector plot can **color arrows by magnitude** (with
a colorbar) instead of solid black.

The FEM · Results toolbar also carries **1D Details…**, which opens a non-modal panel
listing every reinforcement line and pile in the model with a utilization badge, and
plotting the selected member's profiles along its own length: mobilized axial force
against its capacity envelope for a reinforcement line, and lateral displacement,
shear, moment, and mobilized soil reaction against depth for a pile. Under the list is
a map of the section with the selected member picked out and named, so a row of the list
is a place on the slope; it is the same drawing the report prints above its member
details. Its own **Field
state** control switches those profiles between the at-failure mechanism and the last
converged solution, exactly as the one on the results view does. **Export** in that
panel writes the current view as a PNG and its plotted series as a CSV, both named for
the field state they were taken at. The button is
dimmed, with a tooltip saying why, for a model that carries neither reinforcement lines
nor piles. See [FEM Reinforcement](../fem/reinforcement.md#inspecting-the-results) and
[FEM Piles](../fem/piles.md#inspecting-the-results) for what the profiles show.

Each reinforcement row also carries the state the line is in — *within capacity*, *near
capacity*, *pullout*, *yielded*, *softened*, *ruptured* or *inactive* — and so does the
line under the plot, with its meaning in the tooltip. Two lines both standing at 100% are
told apart by that word and not by the badge: *pullout* is an end slipping at what its
embedment can develop, *yielded* is the middle of the line at its full tensile capacity.
The states are defined in [The state of a
line](../fem/reinforcement.md#the-state-of-a-line), and they are the words
`print_reinforcement_summary()` prints and a generated report writes.

---

## Display options per view

Each result view has its own **Display** panel (in the left dock) exposing the
options the underlying plot accepts — slice numbers and seep contours on the LEM
solution; nodes / labels / padding on the mesh; variable, levels, vector scale,
flow-line toggles, and a **Water levels** overlay (each head/reservoir boundary's
level for the shown frame) on the seep solution; plot type, deformation controls, and
the converged/at-failure **Field state** switch on FEM results; legend column layout
on every view. Changing an option re-renders the
**cached** result instantly — there's no re-solve. See
[The Display dock](interface.md#the-display-dock).

One option is steady-only: **Flow lines**. A flow net requires divergence-free
through-flow, so it exists only for a **steady** seepage solution; a **transient**
frame is a storage-release state with no stream function, so its Display panel omits
the **Flow lines** and **Base material** controls and turns **Water levels** on by
default. Read a transient frame's flow direction with **velocity vectors** instead.
See [Transient outputs](../seep/transient.md#outputs).

![A Display panel](images/analysis_display_panel.png)

---

## Exporting views: images and DXF

Every canvas has a **Save…** button that exports the **current view**:

![Save view dialog](images/analysis_save_dialog.png)

- **Image** — PNG, PDF, or SVG. PNG prompts for a DPI; vector formats don't. The
  figure is saved at its true inch size, so resolution is independent of the
  on-screen zoom.
- **DXF (rendered view)** — the drawn picture as a layered DXF, good for dropping
  into a CAD document. This is *lossy* as a re-import source; for that, use the
  structured geometry export below.

---

## DXF import and export

Studio exchanges model geometry with CAD through DXF.

**Export geometry** — **File → Export Geometry (DXF)…** writes the structured,
layered DXF: material zones on per-material layers, and profile lines, circles,
reinforcement, distributed loads, and piezo lines on reserved feature layers. This
is the clean companion the importer reads (distinct from the per-view *rendered*
DXF above).

**Import** — **File → Import DXF…** opens a wizard that lists every layer in the
drawing and lets you map each one to an input feature — material zone, profile
line, piezo line, distributed load, reinforcement, failure circles, or *ignore* —
with a material column for zones and profiles.

![DXF import wizard](images/analysis_dxf_import_wizard.png)

Defaults are seeded from xslope's own export layer names and the geometry kind, but
you can override anything, so a DXF drawn in external CAD with arbitrary layer names
maps just as well. Geometry populates the features; non-geometric properties (load
magnitudes, reinforcement strengths, material properties, circle depth) come in as
editable **placeholders** to fill in afterward. The import replaces the current
project (you're prompted to discard first) and is left unsaved so you can complete
it and Save As.

The layer naming and entity conventions are shared with the library's `xslope.cad`
module — see [DXF Import/Export](../usage/dxf.md) for the full layer table and format
details. Reading and writing DXF uses the **ezdxf** package, installed with the `gui`
extra; if it's missing, the import/export actions show an actionable install message.

**Water loads follow what the drawing carries.** A DXF can hold both a `DLOADS` layer
and a `PIEZO` layer, so the import decides the model's
[Water loads](../usage/input_template.md#worksheet-main) mode from the geometry rather
than assuming one. A load block lying along the stretch of ground the piezo line covers
*is* ponded water somebody drew, so the model imports **manual** and that block carries
the reservoir once you give it its pressures — switching to `auto` without deleting it
would count the water twice. A drawing whose only statement of the water is the piezo
line imports **auto**, and the engine derives the reservoir itself. A surcharge
elsewhere on the ground is user data: it is kept, and it does not make the model manual.
Whichever way it goes, the post-import notes say so.

---

## GeoStudio (SLOPE/W) import and export

Studio also exchanges whole models — not just geometry — with other slope-stability
packages: GeoStudio SLOPE/W in both directions, and Rocscience Slide2 and RS2 as imports.

**Import** — **File → Import GeoStudio (SLOPE/W)…** reads a `.gsz`. Unlike DXF, there
is no mapping wizard to work through: a `.gsz` already knows what its geometry means, so
material zones, strengths, water conditions and the seismic coefficient all arrive
identified. The one prompt is **which analysis** to import, because a GeoStudio file
usually holds several over the same geometry — and they can differ in *materials*, not
just in slip surface, so the choice changes the model you get:

Standing water crosses as a **water definition**, not as a load. GeoStudio has no
ponded-water object — SLOPE/W carries the reservoir's weight from the water surface
itself — and so the import brings the piezometric surface across, sets
[Water loads](../usage/input_template.md#worksheet-main) to `auto`, and leaves the
dloads sheet to the model's real loads. The only exception is an analysis whose water is
a solved **SEEP/W field**: that states no water surface anywhere, so the reservoir is
recovered from the head field and written as a load, with the mode set to `manual`. Both
cases are named in the post-import notes.

![GeoStudio import — choosing an analysis](images/analysis_gsz_import_dialog.png)

The other two vendor importers follow the same shape. **File → Import Slide2…** reads a
Rocscience Slide2 model (`.sli` / `.slim` / `.slmd`) and asks the same one question with a
scenario chooser, since a `.slmd` routinely bundles several scenarios over the same
geometry. Slide2 holds its ponded water implicitly too, from the water table, so it
imports the same way GeoStudio's piezometric surface does: the water table is carried,
the mode is `auto`, and the reservoir is derived at solve time.

![Slide2 import — choosing a scenario](images/analysis_slide2_import_dialog.png)

**File → Import RS2 (.fez)…** reads a Rocscience RS2 finite-element model. A `.fez` holds
exactly one model, so the only prompt is the file picker; geometry, materials and water
conditions import directly — including RS2's distributed and ponded-water loads, which
are the one vendor water that stays a **load**. RS2 stores ponded water as an explicit
load object, and its piezometric surface is a whole-domain surface rather than a water
table drawn on the ground, so measuring ground against piezo would invent a plateau of
water the model never had. RS2 models therefore import with
[Water loads](../usage/input_template.md#worksheet-main) on `manual`, carrying RS2's own
load objects, and the notes say why. They arrive as XSLOPE distributed loads with the
[Direction](../usage/input_template.md#worksheet-dloads) RS2 gave them: a *normal* load
perpendicular to the loaded surface, a *vertical* one (a dead-weight surcharge) straight
down. A load RS2 aims straight down by *global angle* is vertical too, and crosses only
when the model's own solved edge tractions confirm which way it points — RS2 writes that
same downward load two contradictory ways, so the file's solved answer settles it, and a
model that was never solved leaves the load out rather than guessing. Materials bring
their Young's modulus and Poisson's ratio across as well, for the FEM. Whatever RS2
defines that cannot cross (its SSR settings, joints, reinforcement, line loads, and loads
at an arbitrary angle or reversed to pull away from the boundary) comes back in the
post-import notes dialog, named one by one.
Each material's pore-pressure source (piezometric line,
r<sub>u</sub> ratio, pressure grid, or a groundwater analysis) is read per material — see
[the mat worksheet's pore-pressure options](../usage/input_template.md#worksheet-mat) for exactly what
each source becomes and which ones import as zero with a warning. RS2's stability
result is an SSR field rather than a slip surface, so the import never carries a failure
surface — you define circles afterward.

**Export** — **File → Export to GeoStudio (SLOPE/W)…** writes the current model out as a
`.gsz`. It needs a polygon-based model (material zones), since a profile-line model has
no regions to map onto.

Both directions replace nothing silently: whatever cannot cross the format boundary —
SLOPE/W's search definition, reinforcement, piles, loads, non-Mohr-Coulomb strengths — is
listed in a notes dialog and in the Log pane, so you know exactly what to re-create by
hand.

Units are one thing a `.gsz` does not carry: the format has no unit-system field, so
XSLOPE infers it from the unit weight of water and refuses to guess when that value is
neither metric nor imperial. See
[GeoStudio Import/Export](../usage/geostudio.md) for the full mapping table, the
per-analysis materials wrinkle, and the limits of export.
