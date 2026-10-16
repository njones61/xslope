# Changelog

What is new, fixed, and changed in each release of XSLOPE: the `xslope` Python
package and XSLOPE Studio. The newest release is first. Installers for every
release are on the [releases page](https://github.com/njones61/xslope/releases).

## Unreleased

### Fixed

- The strength plot in the Materials list view draws τ and σ′ at the same
  scale, so the envelope rises at the friction angle entered. This applies to
  Mohr-Coulomb, power-curve, and Hoek-Brown materials.
- In Studio's tables, `Delete` (`delete` on a Mac keyboard) empties the selected
  cells. Before, a cell had to be opened with a double-click first.
- The Materials list view dims its property fields until a material exists.
  Before, an empty list left the fields open to typing with nothing to hold it.

### Changed

- The RS2 verification page and the FEM samples page give ten strength
  reduction results as the current solver computes them; each moved by 0.4 to
  1.8%, and every comparison keeps its rating.
- A profile line or polygon added in Studio starts with one vertex row to type
  or paste into, and on the material of its own number (line 2 on material 2),
  as the workbook's profile and polygon sheets are laid out.
- In Studio's tables, `Tab` from the last cell of the last row adds a row and
  moves into it, so a table can be typed row after row without reaching for
  **Add row**.

## 1.0.3 — 2026-10-14

### New

- Drop an input workbook (`.xlsx`) or a project package (`.xslz`) onto the Studio
  window to open it, just as **File → Open** would.

### Fixed

- The assistant in the installed app failed on every request with "Unknown
  encoding cl100k_base". It works again.
- The installed Studio no longer rebuilds its font list every time it starts. On
  macOS that held the window back 15 seconds or more.
- Parametric plots: a back-analysis lists FS = 1 once in its legend, and when a
  design study misses its target, the note saying so stays inside the plot and
  clear of the target line.

### Changed

- Studio started without a file opens on a blank project you can edit right
  away, instead of a window where nothing could be edited until **File → New**.
  **File → New** also clears the Log pane.
- The Log pane remembers its font size, and starts at 9 pt on Windows.
- The app icon is redrawn.

## 1.0.2 — 2026-10-12

### Fixed

- In the installed Mac app, every download over https failed with "unable to get
  local issuer certificate": **Open in Studio** links, the update check, and the
  assistant's model lists. The app now carries its own list of trusted
  certificates.

## 1.0.1 — 2026-10-12

### New

- **Help → Report a Problem…** opens the bug-report form on GitHub with your
  XSLOPE version and platform already filled in.
- Bug reports, feature requests, and questions each have their own form on
  GitHub.

### Fixed

- **Open in Studio** links in the documentation were refused with "xslope.org is
  not one of the sites XSLOPE publishes projects from". They work again.

### Changed

- The macOS disk image's window shows an arrow from the app to the Applications
  folder and says what to do.
- Studio's Help links, its reports, and the assistant point to the documentation
  of the released version.

## 1.0.0 — 2026-10-11

### New

- **Slip joints in the finite element analysis.** A reinforcement line can now be
  marked as a joint (the new Joint column on the reinforce sheet), so the soil can
  slide along a sheet such as a smooth liner or geotextile instead of being glued to
  it. Each jointed line takes its own normal and shear stiffness (kn, ks), and Jred
  chooses whether its strength is reduced in a strength reduction run.
- **A joints sheet and a Joints editor.** Slip joints that carry no reinforcement,
  such as rock joints, bedding planes, and the contacts between the courses of a
  block wall, are entered on a new joints sheet or in Studio's Joints editor (a new
  row in the Inputs tree). Each joint has a cohesion, friction angle, tension
  cutoff, and stiffnesses, plus an optional residual strength and dilation angle;
  clicking a joint on the Inputs plot opens its editor.
- **Joint networks.** The Build network button in the Joints editor generates a
  whole joint set (parallel, cross-jointed, or Voronoi blocks) inside the section, a
  material, or a named polygon of the new type `joints`, with a live preview. A
  generated set can be reopened and edited or removed with Remove set, and scripts
  can do the same with `parallel_set()`, `cross_jointed()`, and `voronoi()` in
  `xslope.joints`.
- **Joint results.** On a jointed model the deformed mesh panel draws the blocks
  with each joint face colored by how far it slipped, with Show joints, Color by
  block, Joint width, and Element edges controls. 1D Details plots the normal
  stress, the shear stress against its Mohr-Coulomb limit, and the slip along each
  joint; the report adds a table of joint states; and a packed project keeps the
  joint results.
- **Model checks for joints.** Before a run, the model checks name joint layouts the
  mesher cannot build (jointed lines that meet, a member ending on a joint, a point
  load on one), a jointed line with no strength entered, a bonded sheet that looks
  likely to act as a slip surface, and an iteration limit too short for a jointed
  model to settle. Studio opens the Run FEM dialog on the Hybrid criterion for a
  model with joints.
- **Two new finite element tutorials.** FEM-3, *A Block Wall on Slip Joints*, models
  a segmental block wall course by course, with and without geogrid, and shows when
  a sheet should be jointed rather than bonded. FEM-5, *A Rock Slope on Its Joints*,
  analyzes a jointed rock face, first with two typed joints and then with a
  generated joint network.
- **Displacement vs F plot and a closing summary.** Studio's FEM results gain a
  fourth plot, Displacement vs F, which shows each strength reduction trial and
  whether it stood or failed; the report includes it too. Every strength reduction
  run now ends with a plain-language summary of how each end of the factor-of-safety
  bracket was decided.
- **Continue an undecided search.** When a strength reduction search stops because a
  trial was still undecided at its iteration limit, "Continue with a higher limit…"
  on the FEM results toolbar resumes it from where it stopped instead of starting
  over. From Python, pass the earlier result back with `solve_ssrm(resume=...)`.
- **Lower-bound answers.** When a strength reduction search finds no failure within
  its iteration limit, the run says so and reports the factor of safety as a lower
  bound (FS >= X) in the plots, the report, and Studio's status bar.
- **Accelerate convergence.** The Run FEM dialog has an "Accelerate convergence"
  option that speeds up slow-settling trials and switches itself off when it cannot
  help. It is on by default for models with joints.
- **The AI Assistant reads the documentation.** The full documentation now ships
  with the package, and the assistant searches it before answering a "does XSLOPE
  support…" or "how do I…" question, then cites the page and section it used.
- **Clearer limit equilibrium output.** A search now says in plain words why some
  trial surfaces could not be solved. When the critical surface shows tension near
  the crest and the model has no tension crack, it suggests adding one. Warnings
  print before the results.
- **Warnings in sweeps.** The `sensitivity()` sweep and the AI Assistant's sweep
  helpers now report each step's warnings (unsolved trial surfaces, unreasonable
  solutions, and crest tension) in a `warnings` column and in the printout.
- **Reinforcement section in the documentation.** A new Reinforcement section has an
  Overview and a Reinforcement Types page with how-to entries for each kind of
  support and typical values cited to the FHWA manuals, alongside the LEM and FEM
  reinforcement pages.
- **Piles and Walls section in the documentation.** A new Piles and Walls section
  covers drilled shafts, micropiles, sheet-pile, soldier-pile, segmental block, and
  gravity walls, and piles carrying structural load near a slope, with figures,
  typical values, and which analysis (LEM or FEM) suits each.
- **New finite element reference pages.** The FEM documentation gains a Solver page,
  covering run settings, how trials are decided, and the convergence options, and a
  Joints and Interface Elements page.
- **RS2 Joint Corpus verification page.** A new verification page compares XSLOPE's
  jointed-rock results with the problems in Rocscience's RS2 joint verification
  manual, including block and flexural toppling, plane failure, step-path failure,
  and a Voronoi rock mass.
- **Tutorial additions.** Tutorial LEM-1 gains a step that fills the tension crack
  with water. Seven tutorials now end with the factor of safety by method, and every
  tutorial with a problem drawing opens on it, including new sketches for LEM-1,
  LEM-3, LEM-4, and SEEP-1.

### Fixed

- In Studio, a non-circular Run LEM search used the circular search's convergence
  tolerance and iteration limit. Each surface type now uses its own defaults, so
  Studio and a script give the same answer.
- A limit equilibrium search on a model with an elastic zone, such as a wall drawn
  as an elastic polygon, stopped with an error when a trial surface ran along the
  underside of that zone. The surface now takes the strength of the material
  beneath, so the search runs.
- LEM piles: the pile force direction is now taken with the higher end as the head,
  so a pile entered tip-first no longer pushes with the slide, and a battered pile
  and its mirror image give the same result. The report's pile table lists the
  higher end as the top and gives a battered pile's batter.
- FEM meshing: reinforcement and pile lines that meet (a tieback on a pile, or
  crossing bars) now share a node and are connected. A line end within the snap
  tolerance of a polygon corner snaps onto it, and a reinforcement line lying along
  a material boundary no longer makes the mesher fail or fall back to a coarse mesh.
- A second FEM run in the same Studio session no longer leaves Run FEM disabled for
  the rest of the session.
- In the reinforcement and line-load editors, clearing an optional cell now leaves
  it blank, as an empty cell reads in the workbook (for example, a blank Tres means
  no residual rather than a brittle bar, and a blank line-load angle means straight
  down). The editor's FEM view now shows Tmax, Lp1, Lp2, Tend1, and Tend2, which the
  FEM uses.
- Sensitivity F vs time in Studio now resets its saved-frame checkboxes when a
  different transient solution is loaded, and a frame you untick stays unticked when
  the dialog reopens.
- In the Global parameters dialog, the unit labels now follow the Units and Time
  choices while the dialog is open.
- The notes shown after importing a GeoStudio, Slide2, RS2, or DXF model now name
  the correct Water loads cell on the main sheet.
- A profile-line model with no bottom elevation was saved with a bottom of 0. The
  cell is now left blank, the workbook still opens, and the model checks ask for the
  bottom elevation before any run.
- The Auto size deformation scale now takes effect on at-failure FEM plots, which
  stayed at 1x before.
- Section plots no longer print matplotlib's "Ignoring fixed y limits" message to
  the console.
- The model checks no longer warn about a development length longer than its
  reinforcement line, which is how a pullout-governed member is entered. They warn
  only when the entered lengths keep the line far below its Tmax everywhere, and
  Tres above Tmax is now a warning only on FEM runs.
- AI Assistant: current OpenAI and Z.ai models are listed and work, very long
  printouts are trimmed before they are sent to the model, and tables and lists in
  the chat panel display correctly. A plot saved over an earlier file of the same
  name now appears in the chat, and a new, empty project is no longer described as
  having a surface type.

### Changed

- The input template is now version 27. It adds the joints sheet, the Joint, kn, ks,
  and Jred columns on the reinforce sheet, room for 30 reinforcement lines (up from
  20), and a `joints` polygon type for the region a joint network is generated in.
- A workbook now always opens, even when it is incomplete or holds values an
  analysis would reject, such as no failure surface, a Mat ID with no material, or a
  blank reinforcement length. The model checks report what to fix when you run an
  analysis, so any workbook Studio saves can be reopened, and seepage or FEM models
  no longer need a placeholder failure surface.
- FEM result plots: panel titles show the factor of safety to three decimals and say
  "at Last Converged" when they show the last trial that converged. The at-failure
  picture now stops once the slope has moved 20% of the model height, so the failure
  stays readable.
- A strength reduction trial that is still moving slowly when it reaches its
  iteration limit is now judged by the trend of its movement, and the solver checks
  whether it comes to rest before counting it as a failure. The summary and the Log
  name the reading that decided each trial and the run's iteration settings.
- The documentation is reorganized. Reinforcement and Piles and Walls are now their
  own sections, and the old page addresses redirect. Installation is the first page
  of the Usage Guide, Seepage Analysis comes before the Limit Equilibrium Method,
  and tables and quoted program output wrap within the page.
- Sample problems that repeated a tutorial were removed (ten limit equilibrium and
  three seepage samples), and their material moved into the tutorials. The remaining
  sample sketches were redrawn to scale, and the piles tutorial is now FEM-4.
- The earth dam rapid-drawdown sample workbooks were corrected to the source
  problem's values (clay psi of 18 and a base at El. 182), which slightly lowers
  their rapid-drawdown factors of safety.
- Tutorial W-1 (The AI Assistant) was re-recorded on current models, with a check of
  each session's work and a section comparing other models. It and the AI Assistant
  reference now say where to create an API key for each provider.
- Program messages and documentation use American spelling throughout (for example,
  the strength editor's "Mobilizable strength" column). The preferred citation is
  now the SoftwareX paper (Jones and Williams, 2026).

## Earlier releases

Releases before 1.0.0 are listed, with their commits, on the
[releases page](https://github.com/njones61/xslope/releases).
