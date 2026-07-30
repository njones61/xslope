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
`ssrm`, `sensitivity` and `reliability`.

Composite analyses inherit. A rapid-drawdown run must satisfy every limit-equilibrium
rule plus its own; a transient seepage run must satisfy every steady-seepage rule; a
reliability run inherits whichever base analysis it is sweeping (`{"base": "fem"}`,
defaulting to `lem`).

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

Both take a `check_inputs` argument, `True` by default:

```python
from xslope.slice import generate_slices

success, result = generate_slices(slope_data, circle=circle)                    # checked
success, result = generate_slices(slope_data, circle=circle, check_inputs=False) # not checked
```

`check_inputs=False` is the escape hatch, for the cases where a refusal would be
wrong rather than helpful. The automated searches use it: they check once at their
own entry and then skip the check on each of the thousands of trial surfaces, since
the model does not change between them. The parametric and reliability engines use
it because a swept or sampled value is a deliberate perturbation of the model, not
a mistake to be refused.

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
| **Water** | The unit weight of water; a material reading a piezometric line that does not exist or stops short of the section; a line no material reads; standing water above the ground surface with no distributed load carrying its weight |
| **Materials** | The pore-pressure, unsaturated-model and strength-model vocabularies; a material inside the geometry with no strength model, no unit weight, or no strength at all; `u = ru` with no ratio; `option = cp` with no undrained strength |
| **Main sheet** | A blank seismic coefficient; a coefficient outside the plausible range, or entered with a sign the limit-equilibrium engine cannot use; crack water deeper than the crack that holds it |
| **Surfaces** | A model with no failure surface at all; a method that cannot use the selected surface family; a model carrying both families where the run did not say which |
| **Ordering** | A load or piezometric polyline entered right to left, or one whose x values rise and then fall |
| **Units** | The unit weight of water and the soil unit weights against the declared system, and against each other when nothing is declared |
| **Mesh** | An element type the seepage solver does not support; a mesh referencing a material the mat sheet does not define; a stored pore-pressure field whose node count does not match the mesh it is used with; a zone element size that is not finer than the global target |
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
