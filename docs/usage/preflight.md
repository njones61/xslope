# Preflight Input Checks

Before an analysis starts, XSLOPE checks that the model actually carries what that
analysis needs. The check lives in the `xslope.preflight` module, and every solver
entry point runs it.

The reason it exists is not the crashes. A run that stops with a confusing message
is a bad afternoon; a run that *completes* on a model missing an input is a wrong
number in a report. A blank pore pressure ratio, a hydraulic conductivity of zero,
a boundary-condition set with nothing driving flow, a stored seepage field computed
on a different mesh — each of those used to run to completion and return an answer.
Preflight names them, in the vocabulary of the input template, before any of it
happens.

## What it reports

Every finding carries one of three severities.

| Severity | Meaning | Effect |
|----------|---------|--------|
| **ERROR** | The run would crash, or its answer is provably wrong | The run is refused |
| **WARNING** | The run proceeds, but the model matches a pattern that has produced wrong answers | Reported, never blocks |
| **INFO** | A default was applied, or an input is inert | Reported |

The bar for an ERROR is deliberately high: there has to be a measurement showing
the answer is wrong, not an opinion that the model is unusual. Anything short of
that is a warning, so that a legitimate but unusual model is never blocked.

Every message names the sheet and the field the way the template does — *"mat
sheet, material 3 ('Core'), column k1"*, not an internal key — because the point of
the message is to tell you which cell to open.

## Running it yourself

```python
from xslope.fileio import load_slope_data
from xslope.preflight import preflight

slope_data = load_slope_data("my_slope.xlsx")
report = preflight(slope_data, "lem", {"surface": "circular", "method": "bishop"})

print(report.format())
if not report.ok:
    for finding in report.errors:
        print(finding.message)
```

`preflight` takes the model, the analysis type, and an optional `selection`
describing what the run chose where the model alone does not say — the LEM method,
the surface family (`"circular"` or `"noncircular"`), whether the run is an
automated search. The analysis types are `lem`, `rapid`, `seep`, `tseep`, `fem`,
`ssrm`, `sensitivity` and `reliability`. Where the model *does* say — the
[**Surface family**](input_template.md#run-options) cell, for the rare file defining
both a circular and a non-circular surface — that statement stands in for a
`selection` the run omits, and a run that states a family of its own is taken instead.

Composite analyses inherit. A rapid-drawdown run must satisfy every limit-equilibrium
rule plus its own; a transient seepage run must satisfy every steady-seepage rule; a
reliability run inherits whichever base analysis it is sweeping (`{"base": "fem"}`,
defaulting to `lem`).

One more `selection` key belongs to time-dependent models. With `u = seep` against a
[transient seepage](../seep/transient.md) march the pore pressures are not in the file
at all: one frame of the march is written into the model immediately before the solver
starts. A script does that by calling `apply_transient_stability_frame` first, so the
entry point's own gate already sees a model carrying the field. An interface cannot —
it has to decide whether the run is startable *before* the frame is staged — so it
states the fact instead:

```python
preflight(slope_data, "lem", {"surface": "circular",
                              "seep_frame": {"times": [30.0]}})
```

*A frame is coming, and this is which one.* The missing-field check then reports the
instant it will read instead of refusing, and both orders end in the same place.
Without the key — which is every ordinary run — a model that needs a seepage field
and carries none is refused, as it should be. Two instants say a rapid drawdown will
stage both its frames, which is what the `{base}_seep2.csv` requirement is really
asking for; one instant supplies stage 1 only and leaves that requirement standing.
A staged frame never stands in for a missing **mesh**: there would be nothing to read
it onto.

The report is a small object: `report.ok` is `True` when nothing would refuse the
run, `report.errors` / `.warnings` / `.infos` split the findings, `report.format()`
renders them as text, and `report.raise_for_errors()` turns an error into a
`PreflightError`.

## When the check runs

**At run time only, never when a file is opened.** Loading a file keeps its own
structural checks — the template version, the sheet parsing, the option
vocabularies — which catch a *corrupt* file rather than an *incomplete* one. A
half-built model must always open, because that is how a model gets built: you draw
the geometry, then the materials, then the water, and at no point in between should
the file refuse to load.

The gate is at the solver entry points instead:

- `generate_slices` — every limit-equilibrium path, single-surface or search
- `build_seep_data` — every seepage analysis
- `sensitivity`, `design`, `back_analysis` — the parametric sweeps
- `reliability` (both engines) and `reliability_fem`

All of them take a `check_inputs` argument, `True` by default:

```python
from xslope.slice import generate_slices

success, result = generate_slices(slope_data, circle=circle)                    # checked
success, result = generate_slices(slope_data, circle=circle, check_inputs=False) # not checked
```

`check_inputs=False` is the escape hatch, for the cases where a refusal would be
wrong rather than helpful. The automated searches use it: they check once at their
own entry and then skip the check on each of the thousands of trial surfaces, since
the model does not change between them.

### Sweeps and reliability runs check twice

A parametric sweep or a reliability analysis is the same model solved many times with
one number changed, so it checks the **base model once** at the door — a defect there
is a defect in every point, and naming it once is both cheaper and clearer than
failing identically nine times.

Then each substituted value is re-checked against **only the rules that read the field
being changed**, which keeps a tornado over eight parameters at nine points cheap. A
step whose value carries an error is **skipped with its reason stated** — a
`success=False` row carrying the rule's own sentence — and the sweep continues; it
never takes the surrounding points down with it. A reliability run, which cannot skip
a perturbation and still form an index, refuses instead, and refuses *before* the
critical-surface search rather than minutes into it.

The reliability engines add their own rules on top of the base analysis: a model with
no standard deviations at all, a standard deviation set on a column the material's
strength model does not read, a deviation larger than its own mean (which the Taylor
series cannot take below zero, though Monte Carlo can, by truncating at the physical
floor), and a standard deviation on an `elastic` material — which has no strength for
it to move, so the uncertainty could never reach the factor of safety. A deterministic
`elastic` zone carrying *no* deviation is an ordinary part of a probabilistic model and
is not refused.

!!! warning "An error means the answer would be wrong, not that the model is unusual"

    Reaching for `check_inputs=False` to get past a refusal skips a check that fired
    because something in the model would make the result incorrect. The message names
    the sheet and the cell; fixing the cell is nearly always faster than working
    around the check, and it is the only route that changes the answer.

## Which options a model can run

The same rules answer a second question: which analyses and which methods this
model *could* run, and why not, for each one it cannot.

```python
from xslope.preflight import capabilities

caps = capabilities(slope_data)
caps["lem_method"]["oms"].available        # False on a non-circular-only model
caps["lem_method"]["oms"].reason           # the sentence explaining why
caps["analysis"]["seep"].available         # False without a mesh
```

Each entry carries availability **together with a reason string**, so an interface
can show an option as unavailable *and say why* rather than accepting the choice and
rejecting the run afterwards. The reason comes from the rule itself, which is what
guarantees that the explanation and the refusal are the same sentence.

The clearest case is method against surface family. The Ordinary Method of Slices
and Bishop's Simplified Method sum moments about a circle centre, so they cannot be
used on a non-circular surface — the centre does not exist, and any value supplied
for it is arbitrary. Janbu, Spencer, Morgenstern-Price, Corps of Engineers and Lowe
& Karafiath take either family. (Composite surfaces — a circle truncated at bedrock
— belong to the circular family and remain valid for all seven.)

For the same reason, `solve_all` skips a method it cannot apply and states why,
rather than refusing the whole run: on a straight-plane problem, refusing everything
because one method is inapplicable would suppress the answers from the six methods
that are.

## What is checked

Rules are grouped by family. The list below is the shape of the coverage rather than
a complete enumeration; `xslope.preflight.rules()` returns every rule with its id,
severity and one-line summary.

| Family | What it looks at |
|--------|------------------|
| **Water** | The unit weight of water; a material reading a piezometric line that does not exist or stops short of the section; a line no material reads; standing water above the ground surface with no distributed load carrying its weight; and, on a model with [automatic water loads](#automatic-water-loads), a transcribed block the engine would derive a second time, a pool the derivation could not measure, and two water definitions that disagree |
| **Materials** | The pore-pressure, unsaturated-model and strength-model vocabularies; a material inside the geometry with no strength model, no unit weight, or no strength at all; `u = ru` with no ratio; `option = cp` with no undrained strength |
| **Main sheet** | A blank seismic coefficient; a coefficient outside the plausible range, or entered with a sign the limit-equilibrium engine cannot use; crack water deeper than the crack that holds it |
| **Surfaces** | A model with no failure surface at all; a method that cannot use the selected surface family; a model carrying both families where the run did not say which; a circle whose `Depth` sits below the base of the model and that cuts no failure surface inside it |
| **Model domain** | A domain whose boundary crosses or retraces itself, or encloses no area — the shape every analysis is bounded by, derived from the zones rather than typed, so a defect in it is invisible on the sheet you are looking at. A max depth left at the elevation of the toe produces one: the base of the model runs back along the ground surface, and slicing, meshing and searching all fail on it with a geometry error naming no field |
| **Ordering** | A load or piezometric polyline entered right to left, or one whose x values rise and then fall |
| **Units** | The unit weight of water and the soil unit weights against the declared system, and against each other when nothing is declared |
| **Mesh** | An element type the seepage solver does not support; a mesh referencing a material the mat sheet does not define; a stored pore-pressure field whose node count does not match the mesh it is used with; a zone element size that is not finer than the global target; and, before a finite element run, a material zone too thin to fit three element rows across its width — which cannot develop a shear band, so the run returns a factor of safety that is too high rather than failing. The zone's own **Size** and the Build-mesh dialog's **Refine thin zones** both answer it, and where a mesh is attached the check measures that mesh rather than inferring anything |
| **Seepage** | A conductivity of zero; `k2` greater than `k1`; missing unsaturated parameters on an unconfined model; a boundary set with no boundary conditions, no specified head, or no gradient |
| **Transient seepage** | A missing time unit; a specific storage or specific yield of zero; a missing or non-positive duration; stage times that are half-set or out of order; a save schedule that reaches past the end of the run; a driving series with no value at t = 0 |
| **Finite element** | A blank or non-positive Young's modulus or Poisson's ratio; a blank tensile cap; K0 with no zone geometry to integrate the overburden through; a strength-reduction zone that contains no mesh elements |
| **Rapid drawdown** | The stage-2 water source each pore-pressure option needs; the `d`/`psi` pair; a post-drawdown pool standing higher than the full pool, or above the ground with no stage-2 load; a stage-2 load that repeats stage 1 |
| **Tension cracks** | A crack at or below the base of the slope; a crack that intersects no failure surface while its water thrust still applies; a depth far past the theoretical `2c/γ` |
| **Reinforcement and piles** | Pile spacing that is blank, zero or negative wherever the run divides by it; a pile or reinforcement line the finite element engine cannot build; a pullout length longer than its own line, or negative; an element that crosses no failure surface |
| **Plausibility** | A modulus far outside the band for its own soil type; a Poisson's ratio below any real soil; a structural modulus outside the range from geosynthetic to steel |

### The cross-analysis findings

Some findings are about a *difference between the engines* rather than a missing
input, and they are the ones easiest to miss by reading a single result:

- A **tension crack** is a limit-equilibrium construction. The finite element engine
  represents a crack constitutively, through tensile strength (`mat!t_cut`), never
  geometrically, and its answer is provably independent of `main!D11`. One file with
  a crack depth therefore poses two different problems, and comparing its own LEM and
  SSRM numbers compares a cracked surface with an uncracked continuum.
- The **seismic coefficient** means different things to the two engines. The
  limit-equilibrium engine takes the magnitude of `main!D13` and applies it in the
  failure-driving direction, ignoring the sign; the finite element engine reads the
  sign as direction, `+k` pushing `+x`. Both are correct for their formulation —
  a continuum code analyses both faces at once and cannot know which one you are
  checking — so preflight states the convention of the engine you are running rather
  than trying to choose for you.
- **Reinforcement** complete for one engine can be incomplete for the other. The
  limit-equilibrium engine applies the `Tmax`/`Lp` capacity envelope directly; the
  finite element engine models each line as a bar element and needs `E` and `Area`.
  A file with capacities but no stiffness runs in the LEM and is refused by the FEM,
  and the warning says so while the LEM answer is being read, not afterwards.

### A note on unit systems

XSLOPE never converts between unit systems — it reads the numbers as entered and
returns results in the same system — so a mislabelled system means every number in
the model is being read in the wrong one. Two signals separate the systems reliably:
the unit weight of water, which physics pins to about 9.81 or 62.4, and the soil
unit weights, whose plausible ranges sit far apart. Preflight warns when either
disagrees with the declaration on the main sheet, or when the two disagree with each
other on a model that declares nothing.

Young's modulus is deliberately *not* used as a unit-system signal. Its plausible
ranges genuinely overlap between the systems — a rock modulus of 5,000,000 kPa is
ordinary SI and also a plausible Imperial value — so no threshold separates them
without flagging correct models. A blank or non-positive modulus is still reported,
because that is a real finite-element input fault rather than a units question.

### A note on the plausibility checks

The plausibility rules ask a different question from the unit-system ones. Those ask
*"is this value in the system the file declares?"*; these ask *"given that system, is
this value in the engineering band for its own field type?"* — which is the version
of the question that can actually be answered, because a material's expected modulus
depends on what kind of material it is. Preflight classifies each material by its own
strength and compares its modulus against that soil type's range.

Two properties keep these honest. The bands are **deliberately loose**, calibrated so
that no correct, reproduced model in the verification corpus trips them: a real
model's modulus legitimately sits anywhere over about two orders of magnitude around
its soil type's midpoint. And a value that *matches* a typical default is never
evidence of anything — `E = 100,000` with `ν = 0.3` is both a common fallback pair and
a perfectly ordinary thing to specify deliberately. These rules report an implausible
magnitude. They never claim a value was left unset.

## Remedies

A rule may name a **remedy**: a fix it could offer for what it found. A remedy is
always offered and never applied on its own. A fix you did not ask for is the same
disease as a silent default — it leaves the model disagreeing with the file you
typed, and you would have no record of the change. Findings carry the remedy's name;
applying one is always an explicit act.

In a script, that act is naming the remedy:

```python
from xslope.preflight import preflight

report = preflight(slope_data, "lem", remedies=["add_ponded_water_load"])
report.applied      # ('add_ponded_water_load:stage1',)
report.model        # the remedied model -- `slope_data` itself is untouched
```

The model that comes back is a copy, and it is the copy that must be handed to the
solver. There is no blanket "fix everything" switch, and there is no mode in which
a remedy applies itself.

Five of them are built:

| Remedy | What it does |
|--------|--------------|
| `reverse_polyline` | Reverses a piezometric line or a distributed-load block entered right to left |
| `add_ponded_water_load` | Adds the standing water as a distributed load, derived from the model's own water definition |
| `switch_to_auto_water` | Sets the main sheet's **Water loads** to `auto` and removes the transcribed blocks the derivation reproduces, keeping every block that is not water |
| `generate_starting_circles` | Fills an empty circles sheet with a starting set derived from the slope geometry |
| `generate_noncircular_surface` | Fills an empty non-circ sheet with a surface tracking the model's weak zone |

A fault can have more than one sensible repair, and an empty surface sheet is the
case: `surface.none_defined` offers both generators, because which one is right
depends on what controls the mechanism rather than on anything the rule can see.
Where the slope's own geometry does, the circles are the answer; where a weak layer
does, no circle can follow the seam. The weak-zone generator is offered only where
it picks a zone on its own — a model with two comparable candidate seams needs the
question *which one* asked, and asking is a dialog rather than a remedy (see
[Studio's zone picker](#a-non-circular-surface-tracking-a-weak-layer)).

Four properties hold for all of them, and each is worth knowing because it decides
what you can rely on.

**What it will do is computed before it does it.** A proposal states the change in
the template's own vocabulary — *"Add 1 distributed-load block, 4680 peak, over
x = -150 to 225, derived from seep bc sheet, head boundary #1 at elevation 302"* —
while you are still deciding. The same computation produces the change, so the
description and the result cannot drift apart.

```python
from xslope.remedies import remedy_proposals, propose, apply_remedy

for p in remedy_proposals(slope_data):
    print(p.key, p.available, p.description or p.reason)

model, finding = apply_remedy(slope_data, "add_ponded_water_load", "stage1")
```

**Applying one produces a finding, not silence.** The error or warning it resolves
is replaced by an INFO recording what changed and which rule asked for it, and that
INFO stays in the report. A model that ran on a synthesised load says so wherever
its factor of safety is reported.

**A remedy declines rather than half-applying.** Every proposal carries
`available` together with a `reason`, so an interface dims the button and explains
why instead of failing when it is pressed — the same behaviour as `capabilities()`,
and `remedy_capabilities()` returns the same shape. The declines are as informative
as the offers: a piezometric line whose x values rise and then fall is not a
reversed line, so the reversal remedy refuses it rather than sorting it into a
third shape; a stage that already carries a load over the same reach is refused a
second one, because that is the double count.

**Nothing is guessed.** A remedy that depends on another rule already being
satisfied says so and stops. The water derivation measures along lines whose x
values increase, so on a line entered backwards it names the reversal remedy rather
than quietly deriving nothing.

### Where the water load comes from

The `add_ponded_water_load` remedy does not invent a load: it derives it from the
model's own statement of where the water stands, in a fixed order of precedence.

1. **The seepage boundary conditions**, wherever a seepage analysis is defined —
   `seep bc` for stage 1, `seep bc (2)` for stage 2. A reservoir or head boundary
   traced along the ground surface states the pool elevation directly, so no
   seepage solution has to be run first. Where the level is a `tseep` time series,
   it is evaluated through the transient march's own interpolation, at t = 0 for
   stage 1 and at the stage-2 time for stage 2 — so a derived load and the seepage
   field it accompanies cannot disagree about where the pool was.
2. **Otherwise the piezometric line** — Line 1, or Line 2 for stage 2.

Only a boundary drawn *on the ground surface* is read as a pool. A head boundary
along a deep aquifer or a model side is a groundwater condition, and reading its
level as a water surface would flood the section with a reservoir nobody drew.

The load itself is the weight of the water between that surface and the ground,
tapering to zero where the two meet, which is the same operation the vendor
importers use to recover a reservoir the source program stored implicitly. On a
model that carries a hand-entered reservoir, the derivation reproduces it: on
`xslope_dam.xlsx` the derived and transcribed loads agree to better than a
thousandth of a percent, which is the check the regression suite runs.

The load follows the ground's own shape, not a column of water measured at each
station. A **vertical face** — a stepped wall, a bench, the front of a gabion
stack — carries the hydrostatic pressure over its wetted height, and the pool ends
at the **shoreline**, the point where the ground actually crosses the water level,
rather than at whichever vertex comes next above it.

### What is not a pool

Two things a water line can do are read as geometry rather than as standing water.

A water surface that is meant to *meet* the ground rarely meets it exactly: a
phreatic surface exiting at a toe, a piezometric line drawn along a flat foreshore,
a line whose tail is carried a little past the edge of the section. What is left is
a wedge a few millimetres or centimetres deep, and it is the coordinates' own
precision rather than water. A block shallower than a **thousandth of the section's
vertical relief** is discarded, and the derivation reports what it discarded. The
fence is wide: across every water-carrying file in the verification corpus the
deepest such residual is eight ten-thousandths of its section's height, and the
shallowest real pool — 4.2 ft of tailwater on a 106 ft dam — is fifty times deeper.

A water surface is also a function of position along the section, so it can hold
only one elevation at a given station. Where a section carries **two pools at
different levels with a vertical face between them** — a dam apron with the
reservoir behind and the tailrace in front, a dewatered trench cut into a seabed —
the derivation cannot tell which pool wets the face. It loads the face from the
higher surface and says so, and a model of that shape should keep its water loads
typed in.

### Automatic water loads

The same derivation can run at every solve instead of once at a button press. The
main sheet's **Water loads** cell decides:

| Mode | Who supplies the weight of standing water |
|------|-------------------------------------------|
| `auto` | The engine, derived from the water definition each time the model is solved. The dloads sheets carry **non-water** loads only — a surcharge, a footing, traffic |
| `manual` | You, on the dloads sheets. This is what every file written before template version 22 means, whether or not it says so |

The default follows the template version, and that is a correctness requirement rather
than a preference: an older file already carries its water load typed in, and deriving a
second one under it would count the reservoir twice.

Which rules apply follows the mode, because the two modes have opposite failure
conditions. In `manual` mode the risk is water with no load, so the ponded-water
warnings and the stage-2 load rules apply. In `auto` mode those cannot arise, and three
rules take their place:

| Rule | What it catches |
|------|-----------------|
| `water.auto_dload_double_count` | A block on a dloads sheet that **is** the derived water — the same load entered twice, once by you and once by the engine. Detected by the derivation itself: the block is compared against the derived load over the same reach |
| `water.auto_derivation_empty` | Water standing above the ground surface from which the engine derived nothing, and why — a line entered right to left, a blank unit weight of water. In automatic mode this is how a reservoir goes missing quietly |
| `water.sources_disagree` | Seepage head boundaries and a piezometric line that describe **different** pools. They should say the same thing, so one of them is stale — and the boundary conditions are the one the run uses |

Switching a model over is the `switch_to_auto_water` remedy, and it is the better of the
two water remedies wherever it applies: writing blocks into a sheet records a snapshot
that goes stale the moment the pool moves, while the mode is recomputed at every solve.
It states what it will remove before it removes anything —

```
Set main!D23 (Water loads) to auto and remove 1 block from the dloads sheet (#1),
which the derivation from piezo sheet, Piezometric Line 1 reproduces to within 0%
(resultant 11886.4).
```

— it keeps every block that is *not* water verbatim, and it declines rather than removing
a block the derivation does not reproduce. That mismatch is a finding about the file:
either the transcription or the water definition is wrong, and a remedy is not the place
to settle which.

### What an imported model arrives as

Vendor formats hold standing water the way XSLOPE's automatic mode does — implicitly,
from the water surface — so an import carries the water definition rather than a load,
and the mode it arrives in is decided by what the vendor's file actually states:

| Import | Mode | Why |
|--------|------|-----|
| [GeoStudio `.gsz`](geostudio.md) with a **piezometric surface** | `auto` | SLOPE/W has no ponded-water object; it carries the reservoir from the surface, and so does XSLOPE |
| GeoStudio `.gsz` fed by a **SEEP/W field** | `manual` | The file states no water surface anywhere. The reservoir is recovered from the head field and written as a load, because nothing downstream can re-derive one from an imported field |
| **Slide2** (`.sli` / `.slim` / `.slmd`) | `auto` | Same implicit model as SLOPE/W: the water table is the water definition |
| **RS2** (`.fez`) | `manual` | The opposite case. RS2 stores ponded water as an explicit load object, and its piezometric surface is a whole-domain surface — measuring ground against it would invent a plateau of water the model never had |
| **DXF** (Studio's wizard) | per drawing | A load block tracing the pool is water somebody drew, so the model is `manual`; a piezo line and nothing else is `auto` |

An imported model's caveat list names the mode and the reason in every case. Nothing
imports `auto` while also carrying the water on its dloads sheet, which is what the
`water.auto_dload_double_count` rule above would catch.

## Generating a starting surface

A limit-equilibrium search has to start somewhere, and where it starts decides what
it finds — the adaptive search refines whatever neighbourhood its starting surface
puts it in. `xslope.generators` builds that starting surface from the geometry the
model already carries, in either family: a set of trial circles, or a non-circular
polyline tracking a weak layer.

### Circles

```python
from xslope.generators import generate_starting_circles, slope_geometry

geom = slope_geometry(slope_data)     # toe, crest, height, face segments
circles = generate_starting_circles(slope_data)
slope_data["circles"] = circles
```

The rules are standard practice, applied per slope face:

- **Centre** — `Xo` halfway between toe and crest, `Yo` at the toe elevation plus
  twice the slope height.
- **A toe circle** — one passing *through* the toe, `R = dist(centre, toe)`. This is
  not the same surface as a circle whose bottom sits at the toe *elevation*.
- **A circle tangent to the base of each distinct material layer**, with `Depth` set
  to that layer's base elevation. `Depth` is an elevation, not a depth below ground;
  the radius is `Yo - Depth`.
- **A skimming circle on a cohesionless face.** Where the material exposed on the
  face has `c = 0`, the critical surface is an arbitrarily shallow face-parallel
  slide whose factor of safety is `tan φ / tan β`, independent of depth. A set
  seeded only with toe and base circles cannot reach it and converges to a deep
  local minimum with a non-conservatively high answer. A large-radius circle
  approximates the plane; its centre lands far outside the model, which is expected.
  The steepest face *segment* governs, not a crest-to-toe chord, which on a benched
  face averages the benches away.

A dam reports two faces sharing one crest and is seeded on both, since either can be
the critical one.

Every generated circle has to be one a search could actually be handed, so a
candidate is kept only if it daylights on the ground surface **inside** the model —
never at a vertical edge — and stays inside the domain polygon. Where a section is
transcribed too narrow for the standard centre, the centre is lowered until a circle
fits, and the result says so: a section that needs it is cropped, and the real repair
is to widen the geometry by about twice the slope height beyond the toe and the
crest. Where nothing fits at all, the generator returns nothing and gives that
reason, rather than offering a surface the slicer would refuse.

The generated set is a starting set, not an answer. It exists so that the search has
a seed in every family that could win.

### A non-circular surface, tracking a weak layer

Some slopes fail along a weak layer rather than along their own geometry, and no
circle passes through that mechanism: the surface runs flat inside the seam for most
of its length and turns up sharply at each end. A circular search cannot find that
shape at all, so a model with a weak seam needs a non-circular starting surface, and
until now the only way to get one was to read it off a drawing by hand.

```python
from xslope.generators import generate_noncircular_surface, rank_weak_zones

result = generate_noncircular_surface(slope_data, report=True)
if result["surface"]:
    print(result["summary"])              # which zone, and why
    slope_data["non_circ"] = result["surface"]
else:
    for zone in result["candidates"]:     # no zone was clearly the weakest
        print(zone.name, zone.tau)
```

**Which layer is the weak one** is decided by ranking every material zone on the
shear strength it can mobilise *at the stress it actually carries*:

```
tau = c + sigma'_n · tan(phi)
```

with `sigma'_n` taken from the soil column above the zone, less pore pressure. That
is the only quantity comparable between materials. Neither cohesion nor friction
angle is: a `c = 0, φ = 35°` sand outranks a `c = 50 kPa, φ = 0` clay on cohesion and
loses on friction, and neither answer means anything. It also spans every strength
option on the `mat` sheet, because each reduces to a strength at a normal stress —
undrained `cp` is already one, a Hoek-Brown rock mass is linearised at that stress by
`hoekbrown.hb_tangent`, and a power envelope is evaluated on it. An `elastic`
material cannot fail, so it is not a candidate.

**When one zone is clearly the weakest** — its strength at or below 0.60 of the next
weakest — the generator seeds on it and states the choice and the reason:
*"seeding on 'Weak Layer' — mobilisable strength 22.3 against 67.1 for the next
weakest ('Soil 1')"*. **When two are comparable it returns the ranked candidates
instead of guessing**, and the caller asks. That is not a shortfall of the ranking;
it is a property of the sections. Guo & Griffiths' embankment-over-foundation pair
fails deep on one set of numbers and shallow on another with the *same* two zones,
and both are in the verification corpus. Passing `zone=` (a polygon index or a
material name) names the zone explicitly and skips the question, which is also how a
script overrides an automatic pick.

The 0.60 threshold is measured, not chosen. Every corpus weak-seam row was seeded on
each of its ranked zones and searched, and the converged factor of safety compared
against the row's standing value: every ranking at or under 0.56 picks a zone whose
seeded search reaches or beats it, and the first ranking that picks the *wrong* zone
is 0.63.

**The shape** is a track with a ramp at each end:

- **The track** runs at a tenth of the zone's local thickness above its base — not on
  the base, because a surface lying exactly on a material boundary is a slicing
  hazard, and not far above it, because the base is where the strength contrast the
  mechanism exploits actually is. On a thick zone the offset is capped, so "just
  above the base" means the same thing at any thickness.
- **Its ends** sit under the toe and under the crest, pulled in to where the zone
  actually exists. A zone is used over its longest *continuous* run: a seam that
  pinches out mid-section, or surfaces and re-enters, is two zones as far as a
  sliding mass is concerned, and the run that covers most of the slope is the one
  that carries a mechanism.
- **Where the zone reaches the ground**, the track simply runs out to it and that is
  the entry or exit point — the zone is its own ramp, which is the shape of Griffiths
  & Lane's dipping band and the reason their published critical surface is nothing
  but the band.
- **Where it does not**, a straight ramp at the Rankine wedge angle carries the
  surface to the ground: `45 + φ/2` for the scarp dropping in behind the crest and
  `45 − φ/2` for the wedge pushing out at the toe, using the `φ` of the *overburden*
  the ramp cuts rather than the seam's. The toe-side end daylights **at the toe**
  rather than out on the flat beyond it: a surface emerging past the toe has to shear
  a block of ground that resists it and drives nothing.
- **Every point carries an explicit Y and an explicit Movement** — `Free` at the two
  ends, `Horiz` along the track. A blank Y reaches the slicer as a `TypeError`, and a
  blank Movement silently means `Fixed`, which would freeze the search on the surface
  it was handed.

As with the circles, the generator refuses rather than offering a surface the slicer
would reject, and says why: a model with one material zone has no weak layer to
track (its mechanism follows the slope's own geometry, which is what a circular
search finds), a vertical wall has no horizontal run for a track to follow, and a
section with no room beyond the slope has nowhere for a ramp to daylight.

Like the circles, the result is a *starting* surface. Hand it to
`noncircular_search`, which moves the points to the critical position from there.
