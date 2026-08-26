# Copyright 2025 Norman L. Jones
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Preflight -- declarative dependency checking for an analysis run.

A model that is missing a required input either crashes with a message naming no
field, or runs and returns a confidently wrong number. The second is worse than the
first: a crash is a bad afternoon, a silent default is a wrong answer in a report.
This module is the single place that decides, *before* a solver starts, whether a
model carries what the analysis it is about to run actually needs.

Two entry points, one rule set
------------------------------

- :func:`preflight` is the **run gate**. It returns a :class:`PreflightReport` of
  findings; an ``ERROR`` finding means the run would crash or would produce a
  provably wrong number, and the solver entry points refuse on it.
- :func:`capabilities` is the **interface query**. It returns per-option
  availability *together with the reason string*, so a dialog can dim an
  unavailable option and explain why in a tooltip rather than letting the user
  choose it and then rejecting the run.

Both read the same registry, so an option that is dimmed is dimmed *because* of the
rule that would have refused the run, and the two can never disagree.

The severity contract
---------------------

``ERROR``
    The run would crash, or would produce a provably wrong number. Refused. The
    test is falsifiability: there has to be a measurement showing the answer is
    wrong, not an opinion that the model is unusual.
``WARNING``
    The run proceeds and the answer may well be fine, but the model matches a
    pattern that has produced wrong answers before. Never blocks in the library.
``INFO``
    A default was applied, or an input is inert. Worth stating in a report; never
    interrupts.

Two properties follow, and both are load-bearing:

- **Validation is loud, never silent.** A rule that fires produces a finding the
  caller can see, rather than a ``warnings.warn`` that a GUI or a scripted run
  never surfaces.
- **A rule names the PARAMETER first, in the vocabulary the interface labels it
  with**, and carries the spreadsheet coordinates at the end as a locator --
  "Material 3 ('Core') has no hydraulic conductivity: k1 is 0 ... (Materials table;
  mat sheet)", never ``materials[2]['k1']`` and never a bare "mat!k1" at the front.
  A Studio user edits a named parameter in a named editor; the cell reference is
  what they need only once they go looking for the sheet.

When the gate runs
------------------

**At run time only, never at file load.** Loading keeps its existing structural
checks -- template version, parse, option vocabulary -- which catch *corruption*,
not *incompleteness*. A half-built file must always load, because a model is built
incrementally; it is only when the user asks for an answer that the question "does
this model support that answer?" becomes meaningful.

Writing a rule
--------------

Rules are declared with the :func:`rule` decorator and carry an id, a severity,
the analysis types they apply to, and a check::

    @rule("mat.ru_zero", ERROR, ("lem", "fem"),
          "u = ru with no pore pressure ratio")
    def _check_ru(ctx):
        for i, m in ctx.strength_materials():
            if m.get("u") == "ru" and not _pos(m.get("ru")):
                yield (f"{ctx.mat_label(i)} selects u = ru, but ru is blank or "
                       f"zero, so no pore pressure would be generated.")

A check receives a context object and returns ``None``, a single message string, or
an iterable of messages. A message may be paired with a severity override --
``yield (WARNING, "...")`` -- for the rules whose severity depends on magnitude. A
rule that raises is reported as an INFO naming itself and never blocks a run: a
broken rule must not be able to stop a valid model.

A rule may also carry a **remedy** name: a named fix it can offer but never apply on
its own. The name is declared here (:data:`REMEDIES`, and the ``remedy`` field on a
finding) and the transformation lives in :mod:`xslope.remedies`, which computes what
it would do before it does it. The governing principle is one line: *a remedy is
always offered and never applied silently* -- a fix the user did not ask for is the
silent-default disease wearing a helpful face. In a script that means naming it::

    report = preflight(slope_data, "lem", remedies=["add_ponded_water_load"])
    report.model        # the remedied copy; the original is untouched
    report.applied      # what was applied, and to which target
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass, field
from typing import Callable, Optional, Sequence, Tuple

# ---------------------------------------------------------------------------
# Vocabulary
# ---------------------------------------------------------------------------

#: Severity: the run would crash or is provably wrong. Refused at the solver entry.
ERROR = "error"
#: Severity: the run proceeds, but the model matches a known bad pattern.
WARNING = "warning"
#: Severity: a default was applied or an input is inert. Never interrupts.
INFO = "info"

#: Severities in decreasing order of urgency, for sorting and filtering.
SEVERITY_ORDER = (ERROR, WARNING, INFO)
_SEVERITY_RANK = {s: i for i, s in enumerate(SEVERITY_ORDER)}

# ---------------------------------------------------------------------------
# Where an input lives: the trailing locator every message ends with.
#
# The message itself leads with the parameter's own name, in the words the
# interface labels it with -- "Tension crack depth", "Seismic coefficient k",
# "Material 3 ('Core')". These say WHERE to go and change it: the Studio editor
# that owns the input, then the sheet cell for anyone working in the spreadsheet
# instead. Both, in that order, because a Studio user never sees a cell reference
# and a spreadsheet user never sees an editor title.
# ---------------------------------------------------------------------------

#: Editor titles, verbatim from the dialogs (``studio.editors.CATEGORY_EDITORS``).
_AT_MAT = "(Materials table; mat sheet)"
_AT_CIRCLES = "(Circles; circles sheet)"
_AT_NONCIRC = "(Non-circular surface; non-circ sheet)"
_AT_PIEZO = "(Piezometric lines; piezo sheet)"
_AT_DLOADS = "(Distributed loads; dloads sheet)"
_AT_DLOADS2 = "(Distributed loads, Set 2; dloads (2) sheet)"
_AT_SEEPBC = "(Seep BC; seep bc sheet)"
_AT_PILES = "(Piles; piles sheet)"
_AT_REINF = "(Reinforcement; reinforce sheet)"
_AT_PROFILE = "(Profile lines; profile sheet)"
_AT_POLYGON = "(Polygons; polygon sheet)"
_AT_TSEEP = "(Transient seepage; tseep sheet)"
#: The mesh is built by a dialog rather than typed on a sheet.
_AT_MESH = "(Build mesh)"


def _at_global(cell):
    """The locator for one Global parameters field, e.g. ``main D11``."""
    return f"(Global parameters; main {cell})"


def _at_dloads(stage):
    """The locator for a distributed-load block of ``stage``."""
    return _AT_DLOADS if stage == 1 else _AT_DLOADS2


#: Every analysis type the registry knows. ``rapid``, ``sensitivity`` and
#: ``reliability`` are *composite*: they inherit a base analysis's requirements and
#: add their own (see :func:`expand_analysis`).
ANALYSES = ("lem", "rapid", "seep", "tseep", "fem", "ssrm",
            "sensitivity", "reliability")

#: Composition. A composite analysis satisfies every rule of the analyses it
#: inherits, plus its own. ``sensitivity`` and ``reliability`` inherit whichever
#: base analysis they are sweeping, which the caller names as
#: ``selection={"base": "fem"}`` (default ``lem``).
_INHERITS = {
    "rapid": ("lem",),
    "tseep": ("seep",),
    "ssrm": ("fem",),
}

#: The LEM methods, in the template's own vocabulary (``main!D14``). Kept in step
#: with ``fileio.LEM_METHODS`` minus the ``all`` meta-selection.
LEM_METHOD_OPTIONS = ("oms", "bishop", "janbu", "corps", "lowe", "spencer", "mprice")

#: Methods whose formulation takes moments about a circle center and which are
#: therefore valid on the circular surface family ONLY (composite surfaces included
#: -- a circle truncated at bedrock still carries ``Xo``/``Yo``/``R``).
CIRCULAR_ONLY_METHODS = ("oms", "bishop")

#: Display names for the messages, so a refusal reads as the dialog's own label.
_METHOD_NAMES = {
    "oms": "the Ordinary Method of Slices",
    "bishop": "Bishop's Simplified Method",
    "janbu": "the Janbu Method",
    "corps": "the Corps of Engineers Method",
    "lowe": "the Lowe & Karafiath Method",
    "spencer": "Spencer's Method",
    "mprice": "the Morgenstern-Price Method",
}

#: Named remedies a rule may offer. A rule declares the remedy it would offer and
#: a finding carries the name, so that the Studio button and the opt-in script
#: form (``preflight(sd, "lem", remedies=[...])``) have one vocabulary to bind to.
#: The transformations live in :mod:`xslope.remedies`; a name here with no
#: implementation is a declaration of intent, never a silent no-op, because
#: nothing applies a remedy without being asked and asking for an unbuilt one
#: raises.
REMEDIES = {
    "set_seismic_zero": "Set the Seismic coefficient k to 0 for a static run "
                        "(Global parameters; main D13).",
    "reverse_polyline": "Reverse a right-to-left load or piezometric line.",
    "add_ponded_water_load": "Add the standing water as a distributed load.",
    "switch_to_auto_water": "Set Water loads to auto (the main sheet's Water loads row) so the engine "
                            "derives them, removing the transcribed blocks they "
                            "replace.",
    "extend_piezo_line": "Extend the piezometric line across the section.",
    "generate_starting_circles": "Generate a starting set of trial circles from "
                                 "the slope geometry.",
    "generate_noncircular_surface": "Generate a non-circular starting surface "
                                    "along the model's weak zone.",
}

#: Capability groups and their options, for :func:`capabilities`. A group's options
#: are the choices an interface offers; a rule that carries ``capability=<group>``
#: is evaluated once per option with that option pinned into the selection, and a
#: rule that fires makes the option unavailable with its message as the reason.
CAPABILITY_GROUPS = {
    "analysis": ("lem", "rapid", "seep", "tseep", "fem", "ssrm",
                 "sensitivity", "reliability"),
    "lem_method": LEM_METHOD_OPTIONS,
}

#: Field names a parametric sweep or a reliability perturbation can substitute, in
#: the vocabulary :func:`xslope.sensitivity.resolve_param` addresses them by. A rule
#: declares the fields it reads (``fields=("gamma",)``) and :func:`preflight` can
#: then be filtered to *only* the rules a substituted value can invalidate --
#: which is what makes a per-step re-check inside a sweep cheap enough to run on
#: every point. A rule that declares no fields is never selected by that filter:
#: an undeclared dependency is not a silent one, it is simply not claimed.
SWEEPABLE_FIELDS = (
    "gamma", "gamma_sat", "c", "phi", "cp", "ru", "d", "psi",
    "pow_a", "pow_b", "pow_c", "pow_d", "hb_sci", "hb_gsi", "hb_mi", "hb_d",
    "k1", "k2", "alpha", "kr0", "h0", "head",
    "t_max", "t_res", "lp1", "lp2", "tend1", "tend2", "spacing",
    "adhesion", "delta",
    "H", "theta", "D", "S", "V_cap", "M_cap",
    "k_seismic", "tcrack_depth", "tcrack_water", "piezo",
)


class PreflightError(ValueError):
    """Raised by a solver entry point when preflight found an ``ERROR``.

    Subclasses :class:`ValueError` so that callers already catching a bad-input
    ``ValueError`` keep working unchanged, while a caller that wants the structured
    findings can catch this type and read :attr:`report`.
    """

    def __init__(self, message, report=None):
        super().__init__(message)
        #: The full :class:`PreflightReport`, including warnings and infos.
        self.report = report


# ---------------------------------------------------------------------------
# Findings and reports
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class Finding:
    """One thing preflight found, in the template's vocabulary.

    Attributes
    ----------
    rule_id : str
        The rule that produced it (``"mat.ru_zero"``), stable across releases so a
        caller can suppress or test one rule by name.
    severity : {"error", "warning", "info"}
    message : str
        User-facing copy, naming the sheet and field. This is the text a dialog
        shows and the text a refusal carries -- there is no second copy in the GUI.
    remedy : str or None
        The name of the PRIMARY remedy this finding could offer (see
        :data:`REMEDIES`). Never applied automatically.
    remedies : tuple of str
        Every remedy the finding's rule offers, primary first — the same tuple as
        :attr:`Rule.remedies`, carried here so an interface can offer all of them
        without going back to the registry to find out there were more. A fault with
        one repair carries a one-item tuple, and ``remedy`` is always its first
        entry.
    """

    rule_id: str
    severity: str
    message: str
    remedy: Optional[str] = None
    remedies: Tuple[str, ...] = ()

    def __str__(self):
        return f"{self.severity.upper()}: {self.message}  [{self.rule_id}]"


@dataclass
class PreflightReport:
    """The result of a :func:`preflight` run: findings, plus what was checked.

    ``model`` is the model the findings describe. It is the model that was passed
    in, unless remedies were asked for -- in which case it is the remedied copy,
    and :attr:`applied` names what was changed. The caller's model is never
    modified either way, so a run that consumes a remedy has to take
    ``report.model`` deliberately.
    """

    analysis: str
    findings: Sequence[Finding] = field(default_factory=tuple)
    selection: dict = field(default_factory=dict)
    model: Optional[dict] = None
    applied: Tuple[str, ...] = ()

    @property
    def errors(self):
        """Every ``ERROR`` finding -- the ones that refuse a run."""
        return [f for f in self.findings if f.severity == ERROR]

    @property
    def warnings(self):
        """Every ``WARNING`` finding."""
        return [f for f in self.findings if f.severity == WARNING]

    @property
    def infos(self):
        """Every ``INFO`` finding."""
        return [f for f in self.findings if f.severity == INFO]

    @property
    def ok(self):
        """True when nothing would refuse the run (warnings and infos are fine)."""
        return not self.errors

    def __bool__(self):
        return self.ok

    def __len__(self):
        return len(self.findings)

    def __iter__(self):
        return iter(self.findings)

    def format(self, min_severity=INFO):
        """Render the findings as text, most urgent first.

        Parameters
        ----------
        min_severity : {"error", "warning", "info"}
            Drop anything less urgent than this. The default shows everything.
        """
        cut = _SEVERITY_RANK[min_severity]
        rows = [f for f in self.findings if _SEVERITY_RANK[f.severity] <= cut]
        rows.sort(key=lambda f: _SEVERITY_RANK[f.severity])
        if not rows:
            return f"preflight ({self.analysis}): no findings."
        head = (f"preflight ({self.analysis}): {len(self.errors)} error(s), "
                f"{len(self.warnings)} warning(s), {len(self.infos)} info.")
        return "\n".join([head] + [f"  - {r}" for r in rows])

    def raise_for_errors(self):
        """Raise :class:`PreflightError` when any ``ERROR`` finding is present.

        The message carries every error, numbered, so a user fixing a model sees
        the whole list rather than one at a time.
        """
        errs = self.errors
        if not errs:
            return
        lines = [f"  {i + 1}. {e.message}  [{e.rule_id}]"
                 for i, e in enumerate(errs)]
        raise PreflightError(
            f"This model cannot be run as a {self.analysis} analysis "
            f"({len(errs)} problem(s) found):\n" + "\n".join(lines)
            + "\n(Fix the inputs above, or pass check_inputs=False to bypass "
              "these checks.)",
            report=self)


@dataclass(frozen=True)
class Capability:
    """Whether one option is available, and -- when it is not -- why.

    ``reason`` is user-facing copy that comes from the rule definition, so a dimmed
    control's tooltip and the refusal the run would have printed are the same
    sentence.
    """

    available: bool
    reason: str = ""

    def __bool__(self):
        return self.available


# ---------------------------------------------------------------------------
# The registry
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class Rule:
    """One declarative check.

    Attributes
    ----------
    id : str
        Dotted, stable identifier (``"seep.k1_nonpositive"``).
    severity : {"error", "warning", "info"}
        The default severity of the rule's findings; a check may override it per
        finding by yielding ``(severity, message)``.
    analyses : tuple of str
        Which analysis types the rule applies to. ``("*",)`` means every analysis.
    check : callable
        ``check(ctx)`` returning ``None``, a message, or an iterable of messages
        (optionally ``(severity, message)`` pairs).
    summary : str
        One line describing the rule, for documentation and rule listings.
    remedy : str or None
        Name of the PRIMARY remedy findings from this rule may offer -- the one a
        :class:`Finding` carries, and the one an interface reaches for first.
    remedies : tuple of str
        Every remedy this rule offers, primary first. A fault can have more than
        one sensible repair: an empty surface sheet takes either a starting set of
        circles or a surface tracked along a weak zone, and which is right depends
        on what controls the mechanism rather than on anything the rule can see.
        Defaults to ``(remedy,)``, so a rule offering one remedy declares it once.
    capability : str or None
        The :data:`CAPABILITY_GROUPS` group this rule gates, if any.
    availability : bool
        True for a rule that asks *can this analysis be started at all?* rather
        than *is this model right?* -- a missing mesh is the example. The workflow
        satisfies an availability rule before the run gate is ever reached (the
        mesh is an argument to the seepage entry point), so :func:`preflight`
        leaves these to :func:`capabilities`, which is what an interface dims from.
        Pass ``include_availability=True`` to evaluate them anyway.
    fields : tuple of str
        The sweepable fields this rule reads (see :data:`SWEEPABLE_FIELDS`). A
        parametric sweep substitutes one field per step and re-checks only the
        rules that name it, so this is the declaration that makes the per-step
        gate both cheap and honest.
    """

    id: str
    severity: str
    analyses: Tuple[str, ...]
    check: Callable
    summary: str
    remedy: Optional[str] = None
    capability: Optional[str] = None
    availability: bool = False
    fields: Tuple[str, ...] = ()
    remedies: Tuple[str, ...] = ()

    def applies_to(self, analyses):
        return "*" in self.analyses or bool(set(self.analyses) & analyses)

    def reads(self, fields):
        return bool(set(self.fields) & set(fields))


_REGISTRY = []
_REGISTRY_IDS = set()


def rule(id, severity, analyses, summary, remedy=None, capability=None,
         availability=False, fields=()):
    """Decorator registering a rule's check function.

    Parameters
    ----------
    id : str
        Stable dotted identifier, unique across the registry.
    severity : {"error", "warning", "info"}
    analyses : tuple of str
        Analysis types the rule applies to, or ``("*",)`` for all of them.
    summary : str
        A one-line description, used by :func:`rules` and the docs.
    remedy : str or tuple of str, optional
        Key into :data:`REMEDIES`, or several of them where the fault has more than
        one sensible repair. The first is the primary -- the one a :class:`Finding`
        carries; see :class:`Rule`.
    capability : str, optional
        Key into :data:`CAPABILITY_GROUPS`: this rule gates that group's options.
    availability : bool, optional
        Mark the rule as an availability question rather than a validity one; see
        :class:`Rule`.
    fields : tuple of str, optional
        Sweepable fields this rule reads; see :data:`SWEEPABLE_FIELDS`.
    """
    if severity not in _SEVERITY_RANK:
        raise ValueError(f"unknown severity {severity!r}")
    for f in fields:
        if f not in SWEEPABLE_FIELDS:
            raise ValueError(f"rule {id!r} names unknown sweepable field {f!r}")
    for a in analyses:
        if a != "*" and a not in ANALYSES:
            raise ValueError(f"rule {id!r} names unknown analysis {a!r}")
    if capability is not None and capability not in CAPABILITY_GROUPS:
        raise ValueError(f"rule {id!r} names unknown capability group {capability!r}")
    offered = ((remedy,) if isinstance(remedy, str)
               else tuple(remedy) if remedy is not None else ())
    for name in offered:
        if name not in REMEDIES:
            raise ValueError(f"rule {id!r} names unknown remedy {name!r}")
    if id in _REGISTRY_IDS:
        raise ValueError(f"duplicate preflight rule id {id!r}")

    def _register(fn):
        _REGISTRY.append(Rule(id=id, severity=severity, analyses=tuple(analyses),
                              check=fn, summary=summary,
                              remedy=offered[0] if offered else None,
                              capability=capability, availability=availability,
                              fields=tuple(fields), remedies=offered))
        _REGISTRY_IDS.add(id)
        return fn

    return _register


def rules(analysis=None, fields=None):
    """The registered rules, optionally filtered to one analysis type.

    Returns them in declaration order, which groups related rules together and is
    the order findings come back in. ``fields`` narrows further to the rules that
    read one of the named sweepable fields -- the per-step subset a parametric
    sweep re-checks after substituting a value.
    """
    out = list(_REGISTRY)
    if analysis is not None:
        wanted = expand_analysis(analysis)
        out = [r for r in out if r.applies_to(wanted)]
    if fields is not None:
        want = set(fields)
        out = [r for r in out if r.reads(want)]
    return out


#: Rules a MODEL EDIT cannot answer: their cure is running the upstream analysis
#: that produces the input they are asking for.
#:
#: A pore-pressure field is not a value anybody types -- it is the output of a
#: seepage solve read onto this mesh -- so telling a user (or an assistant) to
#: "fix" one of these invites the wrong repair: silence it by changing the
#: pore-pressure option, which does not supply the missing physics, it removes
#: the requirement and quietly changes what is being analysed. An interface that
#: reports findings while a model is being BUILT needs to tell the two apart, so
#: the distinction is declared here rather than inferred from wording.
#:
#: Derivation, and how it stays true: these are exactly the rules whose message
#: names an analysis to run ("Run the seepage analysis first", "Re-run the
#: seepage analysis on this mesh"). The guard in
#: ``test/assistant_guardrails_check.py`` reads this module's own source for that
#: phrasing and fails when a rule carries it undeclared, so a sibling added later
#: cannot slip past as a bare error.
STAGED_BY_RUN = {
    "seep_field.missing": "the seepage analysis",
    "seep_field.node_count_mismatch": "the seepage analysis, on this mesh",
}


def rules_for_field(field, analysis="lem", base=None):
    """The rule ids a substituted value of ``field`` can invalidate.

    The per-step gate of a parametric sweep or a reliability perturbation: after a
    setter writes one field, only these rules can have changed their answer, so
    only these are re-evaluated. An unknown field selects nothing, which is the
    honest result -- a field no rule claims to read has no per-step check, and the
    post-step sentinel screen is what still catches it.
    """
    wanted = expand_analysis(analysis, base)
    return [r.id for r in _REGISTRY
            if r.reads({field}) and r.applies_to(wanted) and not r.availability]


def expand_analysis(analysis, base=None):
    """Resolve an analysis name to the full set of rule tags that apply.

    Composite analyses inherit: a rapid-drawdown run must satisfy every LEM rule
    plus its own, and a reliability-on-FEM run must satisfy every FEM rule plus the
    reliability ones. ``base`` names the underlying analysis for the two composites
    that do not imply one (``sensitivity``, ``reliability``); it defaults to
    ``"lem"``.
    """
    if analysis not in ANALYSES:
        raise ValueError(
            f"unknown analysis {analysis!r}; expected one of {', '.join(ANALYSES)}")
    out = set()
    stack = [analysis]
    if analysis in ("sensitivity", "reliability"):
        stack.append(base or "lem")
    while stack:
        a = stack.pop()
        if a in out:
            continue
        out.add(a)
        stack.extend(_INHERITS.get(a, ()))
    return frozenset(out)


# ---------------------------------------------------------------------------
# Small value helpers -- blank-vs-zero is the whole game
# ---------------------------------------------------------------------------

def _num(v):
    """``v`` as a float, or ``None`` when it is blank / non-numeric / non-finite."""
    if v is None:
        return None
    if isinstance(v, str) and not v.strip():
        return None
    try:
        f = float(v)
    except (TypeError, ValueError):
        return None
    if math.isnan(f) or math.isinf(f):
        return None
    return f


def _pos(v):
    """True when ``v`` is a finite number greater than zero."""
    f = _num(v)
    return f is not None and f > 0


def _finite(v):
    """True when ``v`` is a finite number (zero and negatives included)."""
    return _num(v) is not None


def _coords(obj):
    """Coordinate list from a shapely geometry, a list of pairs, or ``None``."""
    if obj is None:
        return []
    if hasattr(obj, "coords"):
        try:
            return [(float(x), float(y)) for x, y in obj.coords]
        except Exception:
            return []
    out = []
    for p in obj:
        try:
            if isinstance(p, dict):
                out.append((float(p["X"]), float(p["Y"])))
            else:
                out.append((float(p[0]), float(p[1])))
        except (TypeError, ValueError, KeyError, IndexError):
            continue
    return out


def _extent(pts):
    """``(xmin, xmax)`` of a coordinate list, or ``None`` when it has none."""
    if not pts:
        return None
    xs = [p[0] for p in pts]
    return (min(xs), max(xs))


def _y_at(pts, x):
    """Linear interpolation of ``y`` at ``x`` along a polyline, or ``None``.

    Points are sorted by x first, so a line entered right-to-left still reads
    correctly here -- the ordering rules report that separately rather than letting
    this silently return nothing.
    """
    if len(pts) < 2:
        return None
    s = sorted(pts, key=lambda p: p[0])
    if x < s[0][0] or x > s[-1][0]:
        return None
    for a, b in zip(s, s[1:]):
        if a[0] <= x <= b[0]:
            if b[0] == a[0]:
                return max(a[1], b[1])
            t = (x - a[0]) / (b[0] - a[0])
            return a[1] + t * (b[1] - a[1])
    return s[-1][1]


#: How far a polyline may backtrack in x, as a fraction of its own span, before it
#: counts as genuinely out of order rather than as a near-vertical riser. A water
#: table following a steep face is routinely transcribed with a vertex or two that
#: step back a little -- rs2_65 backtracks 0.83 over a 225-unit span -- and calling
#: that a zig-zag would refuse a correct model. A real reversal is not subtle.
_ORDER_TOL_FRAC = 0.01


def _monotonic_x(pts):
    """``"ascending"``, ``"descending"``, or ``"mixed"`` for a polyline's x values.

    Backtracking below :data:`_ORDER_TOL_FRAC` of the polyline's own span is read as
    a near-vertical riser rather than a reversal, so a water table hugging a steep
    face is not mistaken for a line entered out of order.
    """
    xs = [p[0] for p in pts]
    if len(xs) < 2:
        return "ascending"
    span = max(xs) - min(xs)
    tol = _ORDER_TOL_FRAC * span if span > 0 else 0.0
    up = sum(b - a for a, b in zip(xs, xs[1:]) if b > a)
    down = sum(a - b for a, b in zip(xs, xs[1:]) if a > b)
    if down <= tol:
        return "ascending"
    if up <= tol:
        return "descending"
    return "mixed"


def _fmt(v):
    """Format a number the way the messages do -- short, and never ``1.0000000001``."""
    f = _num(v)
    return "blank" if f is None else f"{f:g}"


def _defect_point(poly):
    """Where a polygon's boundary meets itself, as ``(x, y)``, or ``None``.

    Shapely names the offending coordinate in its validity explanation (``"Ring
    Self-intersection[0 0]"``); it is read out here so the message can point at a
    vertex the user can find on the sheet. A format this cannot parse costs the
    location, never the finding.
    """
    try:
        from shapely.validation import explain_validity
        text = explain_validity(poly)
    except Exception:
        return None
    m = re.search(r"\[\s*([-\d.eE+]+)[\s,]+([-\d.eE+]+)", str(text))
    if not m:
        return None
    try:
        return (float(m.group(1)), float(m.group(2)))
    except ValueError:
        return None


def _ring_defect(poly):
    """``(kind, point)`` for a polygon that is not a simple region, else ``None``.

    ``kind`` is ``"zero_area"`` when the ring encloses nothing, or
    ``"self_intersection"`` when it crosses or retraces itself -- the two ways a
    polygon can look like a domain on the sheet and be unusable as one. The test is
    shapely's own ``is_valid`` plus ``is_simple`` on the ring, which is what every
    consumer downstream is implicitly relying on.
    """
    try:
        from shapely.geometry import LinearRing
        area = float(poly.area)
        valid = bool(poly.is_valid)
        simple = bool(LinearRing(poly.exterior.coords).is_simple)
    except Exception:
        return None
    if not area > 0:
        return ("zero_area", None)
    if valid and simple:
        return None
    return ("self_intersection", _defect_point(poly))


# ---------------------------------------------------------------------------
# The evaluation context
# ---------------------------------------------------------------------------

class _Ctx:
    """What a check sees: the model, the run, and the derived facts rules share.

    Derived quantities (ground-surface coordinates, material references, load
    extents) are computed once per :func:`preflight` call and cached here, so a
    registry of several dozen rules costs one pass over the model rather than one
    pass per rule.
    """

    def __init__(self, slope_data, analysis, selection=None):
        self.sd = slope_data or {}
        self.analysis_name = analysis
        self.selection = dict(selection or {})
        self.analyses = expand_analysis(analysis, self.selection.get("base"))
        self._cache = {}

    # -- run selection -----------------------------------------------------
    @property
    def method(self):
        """The LEM method for this run, from the selection or ``main!D14``."""
        m = self.selection.get("method")
        if m is None:
            m = self.sd.get("lem_method")
        return str(m).strip().lower() if m else None

    @property
    def surface_family(self):
        """``"circular"``, ``"noncircular"``, or ``None`` when nothing states one.

        A deck may carry both families, and two things can state which one is meant:
        the model itself (the main sheet's Surface family row) and this run. The file's statement is the
        standing answer -- it is what a run dialog opens on, what a scripted run that
        selects nothing gets, and what every consumer outside a dialog reads -- and a
        run that states a family of its own is a change to it, taken live here and
        written back to the cell when the run starts, so the two never disagree for
        longer than one dialog. With neither stating one this is ``None``, which is
        exactly the condition rule ``surface.family_ambiguous`` reports.
        """
        s = self.selection.get("surface") or self.sd.get("surface_family")
        if s:
            s = str(s).strip().lower()
            return "noncircular" if s.startswith("non") else "circular"
        return None

    @property
    def effective_surface_family(self):
        """The family that would actually be analysed, applying today's precedence.

        With no explicit selection the circular family wins whenever circles are
        present -- ``generate_slices`` takes the circle branch first -- which is
        the behaviour the ambiguity warning exists to make visible.
        """
        s = self.surface_family
        if s:
            return s
        if self.has_circles:
            return "circular"
        if self.has_non_circ:
            return "noncircular"
        return None

    @property
    def is_search(self):
        """True when the run is an automated search rather than a single surface."""
        return bool(self.selection.get("search"))

    @property
    def surface_supplied(self):
        """True when the caller handed the run its failure surface directly.

        ``generate_slices(slope_data, circle=...)`` and its non-circular twin take
        the surface as an **argument**; the sheet is then not the source, and a
        model whose circles sheet is empty is not missing anything. Every caller
        that builds a surface itself passes through here -- a search evaluating a
        trial circle, a sensitivity sweep re-slicing one surface per step, a plot
        drawing a result -- so the sheet-emptiness rule must ask this before it
        fires. Without it the rule refuses precisely the runs that need it least.
        """
        return bool(self.selection.get("surface_supplied"))

    @property
    def water_loads_mode(self):
        """``"auto"`` or ``"manual"``. Everything at template v21 and earlier is
        manual: the water load is whatever the user typed on the dloads sheet."""
        mode = self.sd.get("water_loads")
        if mode:
            return str(mode).strip().lower()
        return "manual"

    # -- geometry ----------------------------------------------------------
    def _c(self, key, fn):
        if key not in self._cache:
            self._cache[key] = fn()
        return self._cache[key]

    @property
    def ground(self):
        """Ground-surface coordinates as ``[(x, y), ...]`` (possibly empty)."""
        return self._c("ground", lambda: _coords(self.sd.get("ground_surface")))

    @property
    def ground_extent(self):
        return self._c("ground_ext", lambda: _extent(self.ground))

    @property
    def piezo(self):
        return self._c("piezo", lambda: _coords(self.sd.get("piezo_line")))

    @property
    def piezo2(self):
        return self._c("piezo2", lambda: _coords(self.sd.get("piezo_line2")))

    @property
    def domain(self):
        """The model domain -- the union of the material zones -- or ``None``.

        Everything downstream is bounded by it: a slice is cut between the ground
        surface and a failure surface that has to stay inside this polygon, and a
        mesh fills it. It is derived at load from the zones, never typed.
        """
        return self.sd.get("domain_polygon")

    @property
    def domain_floor(self):
        """The elevation of the lowest point of the domain, or ``None``.

        The same quantity ``search.circular_search`` calls ``depth_floor``
        (``domain_polygon.bounds[1]``), so a rule about a circle's depth and the
        search's own bound are reading one number.
        """
        d = self.domain
        try:
            return float(d.bounds[1])
        except (AttributeError, TypeError, IndexError, ValueError):
            return None

    @property
    def has_circles(self):
        return bool(self.sd.get("circles"))

    def circle_label(self, i):
        """``"Circle 2"`` -- the name the Circles editor gives the row.

        A message leads with this and carries the sheet cell at the end, because a
        Studio user edits a named circle in a named editor and only meets the
        spreadsheet coordinates if they go looking for them.
        """
        return f"Circle {i + 1}"

    def circle_slices(self, circle):
        """True when a stored circle produces slices that stay inside the domain.

        Asked of :func:`xslope.slice.generate_slices` itself rather than guessed
        at, because that function's containment test is what actually decides
        whether a trial surface is admissible -- and it is asked BOTH ways a run
        can consume a circle: as a plain arc, and as a composite surface truncated
        at the domain floor. Either one succeeding means the circle is usable, so a
        rule built on this cannot fire on a model that runs.

        ``check_inputs=False`` because the gate is what is calling: re-entering it
        here would recurse.
        """
        from .slice import generate_slices
        for composite in (False, True):
            try:
                ok, _res = generate_slices(self.sd, circle=dict(circle),
                                           composite=composite, debug=False,
                                           check_inputs=False)
            except Exception:
                ok = False
            if ok:
                return True
        return False

    @property
    def has_non_circ(self):
        return bool(self.sd.get("non_circ"))

    @property
    def materials(self):
        return list(self.sd.get("materials") or [])

    @property
    def mesh(self):
        """The mesh this run will use.

        The selection wins when it carries a ``mesh`` key, because a seepage or
        finite element run is handed its mesh as an argument -- the mesh a run
        uses is frequently built in memory and never stored on the model.
        """
        if "mesh" in self.selection:
            return self.selection["mesh"]
        return self.sd.get("mesh")

    def mat_label(self, i):
        """``"Material 3 ('Core')"`` -- the row as the Materials table names it."""
        try:
            name = self.materials[i].get("name")
        except (IndexError, AttributeError):
            name = None
        base = f"Material {i + 1}"
        return f"{base} ('{name}')" if name else base

    @property
    def referenced_mat_ids(self):
        """Indices of materials a geometry zone actually references.

        A rule about strength or unit weight has to key on this: an unreferenced
        row may legitimately be a seepage-only material or a leftover, and refusing
        a run over a row nothing analyses would be refusing the wrong thing.
        """
        def _build():
            out = set()
            for poly in self.sd.get("polygons") or []:
                mid = poly.get("mat_id") if isinstance(poly, dict) else None
                if mid is not None:
                    out.add(int(mid))
            return out
        return self._c("refmats", _build)

    def strength_materials(self):
        """``(index, material)`` for every material a failure surface could cross.

        Falls back to every material when the model carries no polygon material
        references at all, so a dict built by hand is still checked.
        """
        refs = self.referenced_mat_ids
        for i, m in enumerate(self.materials):
            if refs and i not in refs:
                continue
            yield i, m

    def fem_materials(self):
        """``(index, material)`` for every material the finite element engine reads.

        Deliberately EVERY row, not just the referenced ones: ``build_fem_data``
        loops over the whole material table and demands E, nu and gamma from each,
        so a leftover row with a blank modulus stops a finite element run whether
        or not any zone points at it. That is the behaviour the rules must predict.
        """
        return list(enumerate(self.materials))

    def seepage_materials(self):
        """``(index, material)`` for every material the seepage mesh can reach."""
        mesh = self.mesh
        ids = set()
        if mesh is not None:
            try:
                ids = {int(v) - 1 for v in mesh.get("element_materials", [])}
            except (TypeError, ValueError, AttributeError):
                ids = set()
        if not ids:
            ids = set(self.referenced_mat_ids)
        for i, m in enumerate(self.materials):
            if ids and i not in ids:
                continue
            yield i, m

    # -- water and loads ---------------------------------------------------
    @property
    def uses_piezo(self):
        """True when any material takes its pore pressure from a piezometric line."""
        return any(str(m.get("u", "")).strip().lower() == "piezo"
                   for _, m in self.strength_materials())

    @property
    def uses_seep_u(self):
        """True when any material takes its pore pressure from a seepage field."""
        return any(str(m.get("u", "")).strip().lower() == "seep"
                   for _, m in self.strength_materials())

    def dload_blocks(self, stage=1):
        """``(label, points)`` for every distributed-load block of a stage.

        The label is the one the Distributed loads editor shows -- ``"Distributed
        load #2"``, or ``"Stage-2 distributed load #2"`` for the drawdown set --
        and the sheet it lives on is carried at the end of the message instead
        (:data:`_AT_DLOADS` / :data:`_AT_DLOADS2`).
        """
        key = "dloads" if stage == 1 else "dloads2"
        prefix = "Distributed load" if stage == 1 else "Stage-2 distributed load"
        out = []
        for n, blk in enumerate(self.sd.get(key) or []):
            pts = _coords(blk)
            out.append((f"{prefix} #{n + 1}", pts))
        return out

    def dload_resultant(self, stage=1):
        """Total load a stage's distributed-load blocks carry, per unit width.

        The trapezoidal integral of the block's magnitude along its own length --
        the same quantity the slicer distributes onto slice tops, so two stages can
        be compared without caring how many points each was drawn with.
        """
        key = "dloads" if stage == 1 else "dloads2"
        total = 0.0
        for blk in self.sd.get(key) or []:
            pts = []
            for p in blk or []:
                try:
                    pts.append((float(p["X"]), float(p["Y"]), float(p["Normal"])))
                except (TypeError, ValueError, KeyError, IndexError):
                    continue
            for (xa, ya, na), (xb, yb, nb) in zip(pts, pts[1:]):
                total += 0.5 * (na + nb) * math.hypot(xb - xa, yb - ya)
        return total

    @property
    def seep_bc(self):
        return self.sd.get("seepage_bc") or {}

    def bc_counts(self, which=1):
        """``(n_heads, n_fluxes, n_exit_points)`` for a boundary-condition set."""
        bc = self.sd.get("seepage_bc" if which == 1 else "seepage_bc2") or {}
        return (len(bc.get("specified_heads") or []),
                len(bc.get("specified_fluxes") or []),
                len(bc.get("exit_face") or []))

    @property
    def unconfined(self):
        """True when either boundary set draws an exit face (desaturation possible)."""
        return bool(self.bc_counts(1)[2] or self.bc_counts(2)[2])

    # -- transient seepage -------------------------------------------------
    @property
    def tseep(self):
        """The ``tseep`` control block, or ``None`` on a steady model.

        Absent OR present-but-empty both parse to ``None`` at load, so this is the
        one test for "is this model transient?".
        """
        return self.sd.get("tseep") or None

    @property
    def seep_frame(self):
        """The transient frame this run will stage, or ``None`` when nothing will.

        A stability run against a transient march does not read a stored
        pore-pressure field: it reads ONE instant out of the march and writes it into
        ``seep_u`` immediately before the solver starts. Scripted runs do that by
        calling :func:`xslope.seep.apply_transient_stability_frame` before the entry
        point whose gate then sees a model carrying the field. An interface cannot --
        it has to decide whether to let the user press Run *before* the frame is
        staged -- so it states the fact instead, as
        ``selection={"seep_frame": {"times": [...]}}``: *a transient solution is
        attached and one of its frames will be staged.*

        Presence is the claim, and ``times`` only names which instant(s) for the
        message. An empty ``times`` still means a frame is coming: an interface whose
        time entry is momentarily unreadable refuses the run in its own words (that
        time is what is unreadable, not the seepage solution), and answering "no
        field at all" there would be the wrong sentence.

        ``None`` -- no key, or an explicit ``None`` -- is the ordinary case: nothing
        will supply a field, so a model that needs one and has none is incomplete.
        """
        f = self.selection.get("seep_frame")
        if f is None:
            return None
        return dict(f) if isinstance(f, dict) else {}

    @property
    def seep_frame_times(self):
        """The instants :attr:`seep_frame` names, as floats. Possibly empty."""
        out = []
        for t in (self.seep_frame or {}).get("times") or []:
            v = _num(t)
            if v is not None:
                out.append(v)
        return out

    def series_at_start(self, name):
        """The first value of a named tseep series, or ``None`` when it has none."""
        ts = self.tseep or {}
        vals = (ts.get("series") or {}).get(name)
        if not vals:
            return None
        return _num(vals[0])

    # -- structural elements -----------------------------------------------
    @property
    def piles(self):
        return list(self.sd.get("pile_lines") or [])

    @property
    def reinforcement(self):
        """The raw reinforcement lines -- endpoints plus per-unit-width properties."""
        return list(self.sd.get("reinforcement_lines") or [])

    def pile_label(self, i):
        """``"Pile 2 ('pile row')"`` -- the row as the Piles editor names it."""
        try:
            name = self.piles[i].get("label")
        except (IndexError, AttributeError):
            name = None
        base = f"Pile {i + 1}"
        return f"{base} ('{name}')" if name else base

    def reinf_label(self, i):
        try:
            name = self.reinforcement[i].get("label")
        except (IndexError, AttributeError):
            name = None
        base = f"Reinforcement line {i + 1}"
        return f"{base} ('{name}')" if name else base

    @property
    def surfaces(self):
        """The failure-surface geometries this deck DEFINES, for engagement tests.

        A circle contributes its full boundary and the non-circular sheet its
        polyline. This is what the deck states, not what a search would explore --
        which is why every rule reading it stands down on a search run, where the
        trial set is generated and the deck's circles are only seeds.
        """
        def _build():
            from shapely.geometry import LineString, Point
            out = []
            for c in self.sd.get("circles") or []:
                try:
                    xo, yo, r = float(c["Xo"]), float(c["Yo"]), float(c["R"])
                except (KeyError, TypeError, ValueError):
                    continue
                if r > 0:
                    out.append(Point(xo, yo).buffer(r, 128).boundary)
            pts = []
            for p in self.sd.get("non_circ") or []:
                try:
                    pts.append((float(p["X"]), float(p["Y"])))
                except (KeyError, TypeError, ValueError):
                    continue
            if len(pts) >= 2:
                out.append(LineString(pts))
            return out
        return self._c("surfaces", _build)

    def segment_engages_surface(self, x1, y1, x2, y2):
        """True when a structural line crosses any surface the deck defines.

        ``None`` when there is nothing to test against, so a caller can tell "it
        misses" from "there was no surface".
        """
        from shapely.geometry import LineString
        surf = self.surfaces
        if not surf:
            return None
        try:
            seg = LineString([(float(x1), float(y1)), (float(x2), float(y2))])
        except (TypeError, ValueError):
            return None
        return any(seg.intersects(s) for s in surf)

    @property
    def slope_height(self):
        """The ground surface's vertical relief, or ``None`` when it has none."""
        ys = [p[1] for p in self.ground]
        return (max(ys) - min(ys)) if len(ys) >= 2 else None

    # -- parametric sweep and reliability run selection ---------------------
    #
    # A sweep is described by the run, never by the file (sensitivity.py:19-22),
    # so these read the selection. They are the only place the sweep vocabulary is
    # interpreted, and every sensitivity/reliability rule reads them rather than
    # reaching into ``selection`` itself.

    @property
    def swept_ref(self):
        """``(kind, name, field)`` for the parameter a sweep varies, or ``None``.

        Accepts the same shapes :func:`xslope.sensitivity.resolve_param` does: a
        canonical ``"mat:Clay:c"`` string, a tuple, or a two-part
        ``"global:k_seismic"``.
        """
        def _build():
            ref = self.selection.get("param")
            if ref is None:
                return None
            parts = list(ref) if isinstance(ref, (tuple, list)) \
                else str(ref).split(":")
            if len(parts) == 2:
                return (str(parts[0]).strip().lower(), "", str(parts[1]).strip())
            if len(parts) == 3:
                return (str(parts[0]).strip().lower(), parts[1], str(parts[2]).strip())
            return None
        return self._c("swept_ref", _build)

    @property
    def swept_field(self):
        """The field name a sweep substitutes, or ``None`` for a ``modify=`` sweep."""
        ref = self.swept_ref
        return ref[2] if ref else None

    @property
    def swept_material(self):
        """``(index, material)`` for the material a sweep varies, or ``(None, None)``.

        Resolved by name, case-insensitively, the way ``resolve_param`` resolves it,
        or by 1-based index when the reference carries one.
        """
        ref = self.swept_ref
        if not ref or ref[0] not in ("mat", "seep"):
            return (None, None)
        name = ref[1]
        if isinstance(name, int) and not isinstance(name, bool):
            i = name - 1
            if 0 <= i < len(self.materials):
                return (i, self.materials[i])
            return (None, None)
        want = str(name).strip().lower()
        for i, m in enumerate(self.materials):
            if str(m.get("name", "")).strip().lower() == want:
                return (i, m)
        return (None, None)

    @property
    def swept_values(self):
        """The values a sweep will step through, as floats (possibly empty)."""
        out = []
        for v in self.selection.get("values") or ():
            f = _num(v)
            if f is not None:
                out.append(f)
        return out

    @property
    def sweep_mode(self):
        """``"lem"``, ``"fem"`` or ``"seep"`` -- which engine evaluates each point."""
        return str(self.selection.get("mode") or "lem").strip().lower()

    @property
    def reliability_engine(self):
        """``"taylor"`` (default), ``"mc"`` or ``"rs"`` -- which reliability engine runs.

        They have genuinely different preconditions: the Taylor series must
        decline a standard deviation larger than its mean, while the two sampling
        engines handle it by truncating the draws at the physical floor.
        """
        key = str(self.selection.get("engine") or "taylor").lower().replace("-", "_")
        if key in ("mc", "monte_carlo", "montecarlo"):
            return "mc"
        if key in ("rs", "rsm", "response_surface", "responsesurface"):
            return "rs"
        return "taylor"

    @property
    def reliability_samples_sigmas(self):
        """True when the selected reliability engine DRAWS from the sigmas rather
        than perturbing them -- Monte Carlo and the response surface both do, and
        every rule about truncation, drawdown or the fixed surface applies to both.
        """
        return self.reliability_engine in ("mc", "rs")

    # -- mesh --------------------------------------------------------------
    def element_centroids(self):
        """``[(x, y), ...]`` per mesh element, or ``[]`` when there is no usable mesh.

        The mean of an element's active nodes, which is exactly what
        ``fem._element_centroids`` computes -- and centroid membership is how the
        finite element engine decides which elements a zone overlay contains, so a
        rule about a zone has to ask the question the same way or it would be
        predicting a different run.
        """
        def _build():
            mesh = self.mesh
            if mesh is None:
                return []
            try:
                nodes = mesh["nodes"]
                elements = mesh["elements"]
                types = mesh["element_types"]
            except (KeyError, TypeError):
                return []
            out = []
            for i, elem in enumerate(elements):
                try:
                    n = int(types[i])
                    active = [nodes[int(k)] for k in list(elem)[:n]]
                    if not active:
                        continue
                    out.append((sum(float(p[0]) for p in active) / len(active),
                                sum(float(p[1]) for p in active) / len(active)))
                except (IndexError, TypeError, ValueError):
                    continue
            return out
        return self._c("centroids", _build)


# ---------------------------------------------------------------------------
# Entry points
# ---------------------------------------------------------------------------

def preflight(slope_data, analysis, selection=None, ids=None, skip=None,
              include_availability=False, remedies=None, fields=None):
    """Check a model against everything the named analysis needs.

    Parameters
    ----------
    slope_data : dict
        A model as returned by :func:`xslope.fileio.load_slope_data`, or built in
        memory. Nothing is mutated.
    analysis : str
        One of :data:`ANALYSES`. Composite types inherit their base analysis's
        rules (see :func:`expand_analysis`).
    selection : dict, optional
        What the run chose, where the model alone does not say: ``method`` (a LEM
        method name), ``surface`` (``"circular"`` / ``"noncircular"``),
        ``surface_supplied`` (bool -- the caller passed the surface in rather than
        reading it from the sheet), ``search`` (bool), ``base`` (the underlying
        analysis for a sensitivity or reliability run), ``seep_frame`` (a transient
        seepage frame this run will stage into ``seep_u`` before the solver starts;
        see :attr:`_Ctx.seep_frame`).
    ids : iterable of str, optional
        Evaluate only these rule ids -- for testing one rule in isolation.
    skip : iterable of str, optional
        Suppress these rule ids.
    include_availability : bool, optional
        Also evaluate the availability rules -- the ones asking whether the
        analysis can be started at all (is there a mesh?) rather than whether the
        model is right. They are off by default because the workflow satisfies them
        before the gate is reached: the mesh is an argument to the seepage entry
        point, not a property of the model dict. :func:`capabilities` always
        evaluates them, which is where an interface reads them from.
    remedies : iterable of str, optional
        Named remedies to apply before the model is checked -- opt-in, named one at
        a time, and never a blanket "fix everything". Each is applied to a
        **copy**, which comes back as ``report.model``; the finding it resolves is
        replaced by an INFO recording what changed, so a model that ran on a
        synthesised input says so wherever its answer is reported. A remedy that
        cannot fully succeed declines: it changes nothing, its finding stands, and
        the reason is reported beside it. See :mod:`xslope.remedies`.
    fields : iterable of str, optional
        Evaluate only the rules that read one of these sweepable fields (see
        :data:`SWEEPABLE_FIELDS`). This is the per-step gate a parametric sweep
        uses: a full check of the base model once, then only the rules the
        substituted value can have moved.

    Returns
    -------
    PreflightReport
        Findings in registry order. ``report.ok`` is True when nothing would refuse
        the run; ``report.raise_for_errors()`` turns an error into a
        :class:`PreflightError`.

    Examples
    --------
    >>> report = preflight(slope_data, "lem", {"method": "bishop"})  # doctest: +SKIP
    >>> if not report.ok:                                            # doctest: +SKIP
    ...     print(report.format())

    Ask for a fix explicitly, and run the model that comes back::

        report = preflight(sd, "lem", remedies=["add_ponded_water_load"])
        generate_slices(report.model, circle=circle)      # doctest: +SKIP
    """
    if remedies:
        return _preflight_with_remedies(slope_data, analysis, selection, ids, skip,
                                        include_availability, remedies, fields)
    ctx = _Ctx(slope_data, analysis, selection)
    only = set(ids) if ids else None
    want_fields = set(fields) if fields is not None else None
    drop = set(skip) if skip else set()
    findings = []
    for r in _REGISTRY:
        if only is not None and r.id not in only:
            continue
        if want_fields is not None and not r.reads(want_fields):
            continue
        if r.id in drop:
            continue
        if r.availability and not include_availability and only is None:
            continue
        if not r.applies_to(ctx.analyses):
            continue
        try:
            out = r.check(ctx)
        except Exception as exc:                       # pragma: no cover - defensive
            findings.append(Finding(
                r.id, INFO,
                f"the preflight check {r.id!r} could not be evaluated on this model "
                f"({type(exc).__name__}: {exc}). The run was not blocked by it."))
            continue
        findings.extend(_as_findings(r, out))
    return PreflightReport(analysis=analysis, findings=findings,
                           selection=dict(selection or {}), model=slope_data)


def _preflight_with_remedies(slope_data, analysis, selection, ids, skip,
                             include_availability, remedies, fields=None):
    """:func:`preflight` with named remedies applied to a copy first.

    The findings that come back describe the model that would actually run, not
    the file as saved, which is the only honest way to report a remedied model:
    the registry is re-evaluated *after* the change, so a remedy that resolved
    less than it promised leaves its finding visible. Each applied remedy adds the
    INFO that records it, and a remedy that declined adds an INFO saying so
    without touching anything.
    """
    from . import remedies as _rem

    model = slope_data
    infos, applied = [], []
    for name in remedies:
        target = None
        if ":" in str(name):
            name, target = str(name).split(":", 1)
        proposals = [p for p in _rem.remedy_proposals(model, [name])
                     if target is None or p.target == target]
        if not proposals:
            infos.append(Finding(
                f"remedy.{name}", INFO,
                f"The remedy '{name}' was requested, but this model has nothing for "
                f"it to act on. Nothing was changed."))
            continue
        for proposal in proposals:
            if not proposal.available:
                infos.append(Finding(
                    f"remedy.{proposal.remedy}", INFO,
                    f"The remedy '{proposal.remedy}' was requested for "
                    f"{proposal.label}, and declined: "
                    f"{proposal.reason.rstrip('. ')}. Nothing was changed.",
                    proposal.remedy))
                continue
            model, done = _rem.apply_remedy(model, proposal.remedy, proposal.target)
            infos.append(done)
            applied.append(proposal.key)

    report = preflight(model, analysis, selection, ids, skip, include_availability,
                       fields=fields)
    return PreflightReport(analysis=analysis,
                           findings=list(report.findings) + infos,
                           selection=dict(selection or {}), model=model,
                           applied=tuple(applied))


def _as_findings(r, out):
    """Normalize a check's return value into findings."""
    if out is None:
        return []
    if isinstance(out, str):
        out = [out]
    made = []
    for item in out:
        if isinstance(item, Finding):
            made.append(item)
            continue
        if isinstance(item, tuple):
            severity, message = item
        else:
            severity, message = r.severity, item
        made.append(Finding(r.id, severity, str(message), r.remedy, r.remedies))
    return made


def capabilities(slope_data, selection=None):
    """Which options this model can run, and the reason for each that it cannot.

    The query an interface renders from. For every option in every group of
    :data:`CAPABILITY_GROUPS` it returns a :class:`Capability` carrying
    availability **and** a reason string, so a control is dimmed with an
    explanation rather than being live and rejecting the user afterwards. The
    reason is the rule's own message, so the tooltip and the refusal the run would
    have printed are the same sentence.

    Parameters
    ----------
    slope_data : dict
    selection : dict, optional
        Context the query should assume (e.g. ``{"surface": "noncircular"}`` while
        the user has that family chosen in the dialog).

    Returns
    -------
    dict
        ``{group: {option: Capability(available, reason)}}``.

    Examples
    --------
    Look up one option and read both halves of the answer::

        caps = capabilities(slope_data)
        oms = caps.get("lem_method").get("oms")
        oms.available   # False on a non-circular-only model
        oms.reason      # the sentence saying why, ready to use as a tooltip
    """
    base = dict(selection or {})
    out = {}
    for group, options in CAPABILITY_GROUPS.items():
        gated = [r for r in _REGISTRY if r.capability == group]
        out[group] = {}
        for option in options:
            sel = dict(base)
            analysis = base.get("analysis") or "lem"
            if group == "lem_method":
                sel["method"] = option
            elif group == "analysis":
                analysis = option
            ctx = _Ctx(slope_data, analysis, sel)
            reason = ""
            for r in gated:
                if not r.applies_to(ctx.analyses):
                    continue
                try:
                    found = _as_findings(r, r.check(ctx))
                except Exception:                      # pragma: no cover - defensive
                    continue
                blocking = [f for f in found if f.severity == ERROR]
                if blocking:
                    reason = blocking[0].message
                    break
            out[group][option] = Capability(not reason, reason)
    return out


def check_model(slope_data, analysis, selection=None):
    """Run :func:`preflight` and raise :class:`PreflightError` on any error.

    The one-liner the solver entry points call. Returns the report so a caller can
    also read the warnings it did not refuse on.
    """
    report = preflight(slope_data, analysis, selection)
    report.raise_for_errors()
    return report


def method_surface_reason(method, surface_family):
    """Why ``method`` cannot run on ``surface_family``, or ``""`` when it can.

    Exposed separately so ``solve_all`` can state *why* it skipped a method,
    in the same words the gate would have used, without re-deriving the rule.
    """
    m = str(method or "").strip().lower()
    if m in CIRCULAR_ONLY_METHODS and surface_family == "noncircular":
        return _circular_only_message(m, search=False)
    return ""


# ---------------------------------------------------------------------------
# Reliability refusal messages -- one definition, two consumers
#
# Each of these is BOTH the text of a preflight rule and the text
# :mod:`xslope.reliability` refuses with at run time. They live here, next to the
# rules, so the two can never drift: the sentence a user reads before pressing Run
# is character-for-character the sentence they would have read afterwards.
# ---------------------------------------------------------------------------

#: The standard-deviation columns the reliability engines actually read, in the
#: template's own vocabulary. ``s(d)`` and ``s(psi)`` are deliberately absent: they
#: are loaded and written but no engine reads them (see ``reliability.dead_sigma``).
SIGMA_COLUMNS = {"sigma_gamma": "s(g)", "sigma_c": "s(c)",
                 "sigma_phi": "s(f)", "sigma_cp": "s(c/p)"}


def no_sigmas_message():
    """Why a reliability run cannot start when nothing carries a sigma."""
    return ("Reliability analysis requires standard deviations for at least one "
            "material property. None were provided: s(g), s(c), s(f) and s(c/p) are "
            "blank or zero on every material (Materials table, Standard Deviations; "
            "mat sheet).")


def elastic_sigma_message(mat_name, stated):
    """Why an ``elastic`` material's own standard deviation cannot be honoured.

    ``stated`` is the list of sigma keys the material carries. An elastic material
    with NO sigma is perfectly legal in a probabilistic model (a deterministic
    bedrock layer is the ordinary case) and this message is not produced for it.
    """
    return (f"Reliability analysis cannot vary material '{mat_name}': it is "
            f"elastic (it cannot fail, and a slip surface cannot enter it), so the "
            f"standard deviation(s) {', '.join(stated)} set on it cannot reach the "
            f"factor of safety. Clear them, or give the material a real strength "
            f"model (Materials table; mat sheet).")


def unsupported_option_message(mat_name, option):
    """Why a strength model outside mc/cp has no reliability perturbation set."""
    return (f"Reliability analysis does not support the strength option "
            f"option='{option}' on material '{mat_name}'. Supported: mc, cp "
            f"(Materials table; mat sheet).")


def sigma_exceeds_mean_message(details):
    """Why the Taylor series declines a standard deviation larger than its mean."""
    return ("Reliability: the standard deviation exceeds the mean (COV > 100%) for "
            f"{details}. mean - sigma is negative, which is non-physical. Reduce the "
            "standard deviation(s) so mean - sigma >= 0 (Materials table; mat "
            "sheet).")


def _circular_only_message(method, search):
    name = _METHOD_NAMES.get(method, method)
    valid = ", ".join(_METHOD_NAMES[k].replace("the ", "")
                      for k in LEM_METHOD_OPTIONS
                      if k not in CIRCULAR_ONLY_METHODS)
    if search:
        return (f"A non-circular search cannot be run with {name}: it takes moments "
                f"about a circle center, which a non-circular surface does not have, "
                f"so every trial surface would be rejected. Methods valid for a "
                f"non-circular surface: {valid}.")
    return (f"{name} takes moments about a circle center, so it cannot be used with "
            f"a non-circular surface. Methods valid for a non-circular surface: "
            f"{valid}.")


# ===========================================================================
# RULES
#
# Grouped by family. Every message leads with the PARAMETER, in the words the
# interface labels it with, and ends with the locator saying where to change it --
# the Studio editor, then the sheet cell. Every ERROR carries either a probe or a
# corpus instance behind it.
# ===========================================================================

# ---------------------------------------------------------------------------
# Family: water and unit weight of water
# ---------------------------------------------------------------------------

@rule("water.gamma_water_missing", ERROR, ("lem", "seep", "fem"),
      "Unit weight of water (main!D10) must be present and positive.")
def _gamma_water_missing(ctx):
    g = ctx.sd.get("gamma_water")
    if _pos(g):
        return None
    # Only a model that actually uses water needs it -- a dry LEM deck does not.
    if "lem" in ctx.analyses and not (ctx.uses_piezo or ctx.uses_seep_u
                                      or ctx.sd.get("piezo_line")
                                      or ctx.sd.get("dloads")):
        return None
    return (f"Unit weight of water is {_fmt(g)}. Enter a positive value -- there is "
            f"deliberately no default, because substituting a canonical value of the "
            f"wrong system would rescale every pore pressure in the model "
            f"{_at_global('D10')}.")


@rule("piezo.line_missing", ERROR, ("lem",),
      "A material with u = piezo needs a piezometric line to read.",
      remedy="extend_piezo_line",
      fields=("piezo",))
def _piezo_line_missing(ctx):
    if not ctx.uses_piezo:
        return None
    if len(ctx.piezo) >= 2:
        return None
    names = [ctx.mat_label(i) for i, m in ctx.strength_materials()
             if str(m.get("u", "")).strip().lower() == "piezo"]
    # The wording is the migrated guard's own (slice.py / fem.py), deliberately: the
    # phrase "no piezometric line" is the one users and the guard benchmark know.
    return (f"{names[0]} takes pore pressure from a piezometric line (u = piezo), "
            f"but this model defines no piezometric line -- Piezometric Line 1 "
            f"carries fewer than two points. Every slice base in that material would "
            f"read zero pore pressure -- an unsafe, non-conservative factor of "
            f"safety. Draw the piezometric line, or set u = none for a dry model "
            f"(or 'seep' / 'ru' for the other pore-pressure sources) {_AT_PIEZO}.")


@rule("piezo.extent_short", WARNING, ("lem",),
      "A piezometric line that stops short of the section refuses the slices past it.",
      remedy="extend_piezo_line",
      fields=("piezo",))
def _piezo_extent_short(ctx):
    if not ctx.uses_piezo or len(ctx.piezo) < 2:
        return None
    pz = _extent(ctx.piezo)
    gs = ctx.ground_extent
    if pz is None or gs is None:
        return None
    tol = max(1e-9, 1e-6 * (gs[1] - gs[0]))
    gaps = []
    if pz[0] > gs[0] + tol:
        gaps.append(f"x < {pz[0]:.6g}")
    if pz[1] < gs[1] - tol:
        gaps.append(f"x > {pz[1]:.6g}")
    if not gaps:
        return None
    return (f"Piezometric Line 1 spans x = [{pz[0]:.6g}, {pz[1]:.6g}], which does "
            f"not cover the whole section ([{gs[0]:.6g}, {gs[1]:.6g}]): it stops "
            f"short at {', '.join(gaps)}. A failure surface reaching into that reach "
            f"is refused, because a material with u = piezo would have no "
            f"piezometric surface to read there. This is fine when every trial "
            f"surface stays inside the line {_AT_PIEZO}.")


@rule("piezo.no_consumer", WARNING, ("lem",),
      "A piezometric line no material reads produces no pore pressure.",
      fields=("piezo",))
def _piezo_no_consumer(ctx):
    if len(ctx.piezo) < 2:
        return None
    if ctx.uses_piezo:
        return None
    if any(str(m.get("u", "")).strip().lower() in ("seep", "ru")
           for _, m in ctx.strength_materials()):
        return None
    if any(m.get("gamma_sat") is not None for _, m in ctx.strength_materials()):
        # The line is doing a second, legitimate job: it is the water table that
        # splits moist from saturated unit weight (the gamma_sat "sidecar" model),
        # which is independent of any material's pore-pressure option.
        return None
    if ctx.water_loads_mode == "auto":
        # A third legitimate job (v22): under Water loads = auto the line IS the
        # reservoir definition — the engine derives the ponded-water surface
        # loads from it at solve time. A submerged total-stress model carries
        # exactly this shape: u = none everywhere and the water arriving as
        # load, defined once by the line.
        return None
    return (f"Piezometric Line 1 is drawn, but no material uses u = piezo, so the "
            f"line produces no pore pressure anywhere and the model runs dry. Set "
            f"u = piezo on the materials below the water table {_AT_MAT}, or delete "
            f"the line if the analysis is deliberately total stress {_AT_PIEZO}.")


@rule("water.ponded_no_dload", WARNING, ("lem",),
      "Standing water above the ground surface needs a distributed load.",
      remedy="add_ponded_water_load",
      fields=("piezo",))
def _ponded_no_dload(ctx):
    if ctx.water_loads_mode != "manual":
        return None            # in auto mode the engine derives the load itself
    pz, gs = ctx.piezo, ctx.ground
    if len(pz) < 2 or len(gs) < 2:
        return None
    if _monotonic_x(gs) == "mixed" or _monotonic_x(pz) == "mixed":
        return None            # the ordering rules report this; measuring here would lie
    xs = sorted({p[0] for p in pz} | {p[0] for p in gs})
    wet = []
    depth = 0.0
    for x in xs:
        yp, yg = _y_at(pz, x), _y_at(gs, x)
        if yp is None or yg is None:
            continue
        if yp - yg > 1e-9:
            wet.append(x)
            depth = max(depth, yp - yg)
    if not wet:
        return None
    a, b = min(wet), max(wet)
    # Ignore a hairline overlap: a piezometric line that touches the ground within
    # drawing precision, or crosses it at a single vertex, is not standing water.
    # Measured against the model's own height so the threshold travels between
    # unit systems.
    height = max((p[1] for p in gs), default=0.0) - min((p[1] for p in gs), default=0.0)
    if b - a <= 0.0 or depth < max(1e-9, 1e-4 * height):
        return None
    covered = False
    for _, pts in ctx.dload_blocks(1):
        ext = _extent(pts)
        if ext and ext[0] <= b and ext[1] >= a:
            covered = True
            break
    if covered:
        return None
    return (f"Piezometric Line 1 rises up to {depth:.3g} above the ground surface "
            f"between x = {a:.6g} and x = {b:.6g}, but no distributed load covers "
            f"that reach. Standing water must be entered as a distributed load, or "
            f"its weight is missing from the analysis {_AT_DLOADS}.")


# ---------------------------------------------------------------------------
# Family: automatic water loads (main sheet, Water loads = auto)
#
# The mirror image of the manual-mode rules above. There, the risk is a water
# load nobody entered; here it is a water load entered twice, or a derivation
# that quietly produced nothing, or two water definitions that disagree about
# where the pool is. All three read the engine's own derivation rather than
# re-deriving anything, so a rule and the solve it gates cannot differ.
# ---------------------------------------------------------------------------

def _derived(ctx, stage):
    """The engine's derivation for one stage, computed once per context."""
    from .water import derive_water_loads
    key = f"_derived{stage}"
    if key not in ctx._cache:
        ctx._cache[key] = derive_water_loads(ctx.sd, stage=stage)
    return ctx._cache[key]


@rule("water.auto_dload_double_count", WARNING, ("lem", "rapid"),
      "In auto mode a user dload that IS the derived water counts the pool twice.")
def _auto_double_count(ctx):
    if ctx.water_loads_mode != "auto":
        return None
    from .water import match_water_blocks
    for stage, prefix in ((1, "Distributed load"), (2, "Stage-2 distributed load")):
        derived = _derived(ctx, stage)
        if not derived["blocks"]:
            continue
        blocks = ctx.sd.get("dloads" if stage == 1 else "dloads2") or []
        found = match_water_blocks(blocks, derived["blocks"])
        for i, _j, measures in found["pairs"]:
            yield (f"Water loads is set to auto, so the engine derives the standing "
                   f"water from {derived['source']} itself -- but {prefix} "
                   f"#{i + 1} is that same load, reproducing it to within "
                   f"{100 * measures['worst']:.2g}%. It will be counted twice. "
                   f"Delete the block {_at_dloads(stage)}, or set Water loads to "
                   f"manual if this file is meant to carry its water loads "
                   f"explicitly (the main sheet's Water loads row).")


@rule("water.auto_derivation_empty", WARNING, ("lem", "rapid"),
      "In auto mode, water standing above the ground must actually derive a load.")
def _auto_derivation_empty(ctx):
    if ctx.water_loads_mode != "auto":
        return None
    pz, gs = ctx.piezo, ctx.ground
    if len(pz) < 2 or len(gs) < 2:
        return None
    depth = 0.0
    for x in sorted({p[0] for p in pz} | {p[0] for p in gs}):
        yp, yg = _y_at(pz, x), _y_at(gs, x)
        if yp is not None and yg is not None:
            depth = max(depth, yp - yg)
    from .water import POOL_MIN_DEPTH_FRAC
    # The same fence the derivation uses, deliberately: a hairline the derivation
    # discards as coordinate round-off must not be reported here as a reservoir
    # gone missing. One number, both consumers.
    height = ctx.slope_height or 0.0
    if depth < max(1e-9, POOL_MIN_DEPTH_FRAC * height):
        return None            # no standing water to lose
    derived = _derived(ctx, 1)
    if derived["blocks"]:
        return None
    why = derived["reason"] or "the derivation produced nothing"
    return (f"Water loads is set to auto and Piezometric Line 1 stands up to "
            f"{depth:.3g} above the ground surface, but the engine derived no water "
            f"load at all: {why}. The reservoir's weight is missing from the "
            f"analysis (the main sheet's Water loads row).")


@rule("water.sources_disagree", WARNING, ("lem", "rapid"),
      "Seepage head boundaries and a piezometric line must describe the same pool.")
def _water_sources_disagree(ctx):
    """Both sheets say where the water stands, and they do not agree.

    The seepage boundary conditions win where a seepage analysis is defined, so a
    disagreement is not ambiguous -- it is a stale line, and the run will use the
    one the user is less likely to be looking at.
    """
    from .water import seep_bc_water_line
    pz = ctx.piezo
    if len(pz) < 2 or len(ctx.ground) < 2:
        return None
    bc = ctx.sd.get("seepage_bc") or {}
    if not bc.get("specified_heads"):
        return None
    if _monotonic_x(ctx.ground) == "mixed" or _monotonic_x(pz) == "mixed":
        return None
    try:
        line = seep_bc_water_line(ctx.sd.get("ground_surface"), bc, ctx.sd, 0.0)
    except Exception:
        return None
    if len(line) < 2:
        return None
    worst, where = 0.0, None
    for x, y in line:
        yp = _y_at(pz, x)
        yg = _y_at(ctx.ground, x)
        if yp is None or yg is None:
            continue
        if y <= yg + 1e-9 and yp <= yg + 1e-9:
            continue          # dry ground under both readings: nothing to compare
        if abs(y - yp) > worst:
            worst, where = abs(y - yp), x
    height = ctx.slope_height or 1.0
    if where is None or worst <= max(1e-6, 0.01 * height):
        return None
    return (f"The seepage head boundaries and Piezometric Line 1 describe "
            f"different pools: they differ by {worst:.3g} near x = {where:.6g}. "
            f"They should say the same thing, so one of them is stale. The boundary "
            f"conditions are what the run uses wherever a seepage analysis is "
            f"defined, including for the automatic water load (Seep BC, seep bc "
            f"sheet; Piezometric lines, piezo sheet).")


# ---------------------------------------------------------------------------
# Family: material vocabulary and strength
# ---------------------------------------------------------------------------

@rule("mat.u_vocabulary", ERROR, ("*",),
      "The pore-pressure option u must be one of none / piezo / seep / ru.")
def _mat_u_vocabulary(ctx):
    ok = ("none", "piezo", "seep", "ru")
    for i, m in enumerate(ctx.materials):
        u = m.get("u")
        if u is None or (isinstance(u, str) and not u.strip()):
            continue
        if str(u).strip().lower() not in ok:
            yield (f"{ctx.mat_label(i)} sets u = {u!r}, which is not a "
                   f"pore-pressure option. Expected one of: none, piezo, seep, ru. "
                   f"An unrecognized value used to fall through to zero pore "
                   f"pressure, silently inflating the factor of safety {_AT_MAT}.")


@rule("mat.unsat_vocabulary", ERROR, ("seep", "fem"),
      "The unsaturated model unsat must be one of lf / vg / gard.")
def _mat_unsat_vocabulary(ctx):
    ok = ("lf", "vg", "gard")
    for i, m in ctx.seepage_materials():
        v = m.get("unsat")
        if v is None or (isinstance(v, str) and not v.strip()):
            continue
        if str(v).strip().lower() not in ok:
            yield (f"{ctx.mat_label(i)} sets unsat = {v!r}, which is not an "
                   f"unsaturated conductivity model. Expected one of: lf (linear "
                   f"front), vg (van Genuchten), gard (Gardner). A typo here "
                   f"silently runs the linear-front model instead of the one you "
                   f"asked for {_AT_MAT}.")


@rule("mat.option_vocabulary", ERROR, ("lem", "fem"),
      "The strength model option must be one of mc / cp / pow / hb / elastic.")
def _mat_option_vocabulary(ctx):
    ok = ("mc", "cp", "pow", "hb", "elastic")
    for i, m in ctx.strength_materials():
        o = m.get("option")
        if o is None or (isinstance(o, str) and not o.strip()):
            continue
        if str(o).strip().lower() not in ok:
            yield (f"{ctx.mat_label(i)} sets option = {o!r}, which is not a "
                   f"strength model. Expected one of: mc, cp, pow, hb, elastic "
                   f"{_AT_MAT}.")


@rule("mat.option_missing", ERROR, ("lem",),
      "A material a failure surface can cross needs a strength model.")
def _mat_option_missing(ctx):
    for i, m in ctx.strength_materials():
        o = m.get("option")
        if o is None or (isinstance(o, str) and not str(o).strip()):
            yield (f"{ctx.mat_label(i)} has no strength model: option is blank, "
                   f"and the material is inside the geometry, so a failure surface "
                   f"can pass through it. Choose a strength model (mc, cp, pow, hb "
                   f"or elastic) {_AT_MAT}.")


@rule("mat.gamma_nonpositive", ERROR, ("lem",),
      "A material carrying strength data needs a positive unit weight.",
      fields=("gamma", "gamma_sat"))
def _mat_gamma_nonpositive(ctx):
    for i, m in enumerate(ctx.materials):
        opt = str(m.get("option") or "").strip().lower()
        if opt not in ("mc", "cp", "pow", "hb"):
            continue
        if _pos(m.get("gamma")):
            continue
        yield (f"{ctx.mat_label(i)} has no unit weight: g is "
               f"{_fmt(m.get('gamma'))}. The row carries a strength model "
               f"(option = {opt}), so it is a slope-stability material and needs a "
               f"positive unit weight -- a gamma of zero produces zero slice weights "
               f"and a meaningless factor of safety {_AT_MAT}.")


@rule("mat.no_shear_strength", WARNING, ("lem",),
      "option = mc with c = 0 and phi = 0 is a material with no shear strength.",
      fields=("c", "phi"))
def _mat_no_shear_strength(ctx):
    for i, m in ctx.strength_materials():
        if str(m.get("option") or "").strip().lower() != "mc":
            continue
        c, phi = _num(m.get("c")) or 0.0, _num(m.get("phi")) or 0.0
        if c == 0.0 and phi == 0.0:
            yield (f"{ctx.mat_label(i)} has no shear strength at all: option = mc "
                   f"with c = 0 and φ = 0. If it is cohesionless, enter its friction "
                   f"angle -- nothing can tell 'cohesionless' from 'no data' "
                   f"{_AT_MAT}.")


@rule("mat.cp_zero_strength", ERROR, ("lem",),
      "option = cp with both c and c/p zero gives zero undrained strength.",
      fields=("c", "cp"))
def _mat_cp_zero(ctx):
    for i, m in ctx.strength_materials():
        if str(m.get("option") or "").strip().lower() != "cp":
            continue
        c, cp = _num(m.get("c")) or 0.0, _num(m.get("cp")) or 0.0
        if c == 0.0 and cp == 0.0:
            yield (f"{ctx.mat_label(i)} carries no strength: option = cp but both "
                   f"c and c/p are zero, so the undrained strength Su is zero "
                   f"everywhere in this material {_AT_MAT}.")


@rule("mat.ru_zero", ERROR, ("lem",),
      "u = ru with a blank or zero ratio generates no pore pressure.",
      fields=("ru",))
def _mat_ru_zero(ctx):
    for i, m in ctx.strength_materials():
        if str(m.get("u") or "").strip().lower() != "ru":
            continue
        if not _pos(m.get("ru")):
            yield (f"{ctx.mat_label(i)} selects u = ru, but the pore pressure "
                   f"ratio ru is {_fmt(m.get('ru'))}. No pore pressure would be "
                   f"generated, which is identical to a dry model and "
                   f"non-conservative. Enter a positive pore pressure ratio, or set "
                   f"u = none {_AT_MAT}.")


# ---------------------------------------------------------------------------
# Family: main-sheet scalars
# ---------------------------------------------------------------------------

@rule("main.seismic_missing", ERROR, ("lem",),
      "Seismic coefficient k (main!D13) must be a number.",
      remedy="set_seismic_zero",
      fields=("k_seismic",))
def _seismic_missing(ctx):
    if _finite(ctx.sd.get("k_seismic")):
        return None
    return (f"Seismic coefficient k is blank or not a number -- enter 0 for a "
            f"static analysis. Left blank it reaches every slice as NaN, and the "
            f"Ordinary Method of Slices then reports success with a factor of safety "
            f"of NaN {_at_global('D13')}.")


@rule("main.seismic_magnitude", WARNING, ("lem", "fem"),
      "A Seismic coefficient k outside the usual 0.05-0.3 range is probably a "
      "percent/fraction slip.",
      fields=("k_seismic",))
def _seismic_magnitude(ctx):
    k = _num(ctx.sd.get("k_seismic"))
    if k is None or k == 0.0:
        return None
    a = abs(k)
    if a > 1.0:
        return (f"Seismic coefficient k = {k:g}. k is a fraction of gravity, "
                f"typically 0.05-0.3 -- enter 0.15 for 15%g. At this magnitude some "
                f"methods return a negative factor of safety through the success "
                f"path and others fail without naming k {_at_global('D13')}.")
    if a > 0.5:
        return (INFO, f"Seismic coefficient k = {k:g}, which is high for a "
                      f"pseudo-static analysis (typically 0.05-0.3) "
                      f"{_at_global('D13')}.")
    return None


@rule("main.seismic_negative_lem", WARNING, ("lem",),
      "The LEM takes the magnitude of k and orients it itself; the sign is ignored.",
      fields=("k_seismic",))
def _seismic_negative_lem(ctx):
    k = _num(ctx.sd.get("k_seismic"))
    if k is None or k >= 0:
        return None
    return (f"Seismic coefficient k = {k:g} is negative. The limit-equilibrium "
            f"engine applies the coefficient in the failure-driving direction "
            f"automatically and ignores the sign, so this run is identical to "
            f"k = {abs(k):g}. Enter a magnitude -- in the finite element engine the "
            f"sign does set the direction {_at_global('D13')}.")


@rule("main.crack_water_deeper_than_crack", ERROR, ("lem",),
      "Water in crack (main!D12) cannot exceed Tension crack depth (main!D11).",
      fields=("tcrack_depth", "tcrack_water"))
def _crack_water(ctx):
    w = _num(ctx.sd.get("tcrack_water")) or 0.0
    d = _num(ctx.sd.get("tcrack_depth")) or 0.0
    if w <= 0:
        return None
    if d <= 0:
        return (f"Water in crack = {w:g}, but Tension crack depth is "
                f"{_fmt(ctx.sd.get('tcrack_depth'))}. There is no crack for the "
                f"water to stand in, yet the full crack-water thrust is applied to "
                f"the surface {_at_global('D11-D12')}.")
    if w > d + 1e-12:
        return (f"Water in crack = {w:g} exceeds Tension crack depth = {d:g}. The "
                f"water cannot be deeper than the crack that holds it "
                f"{_at_global('D11-D12')}.")
    return None


@rule("main.lem_method_unknown", ERROR, ("lem",),
      "The named LEM method must be one the package implements.")
def _lem_method_unknown(ctx):
    m = ctx.method
    if not m or m == "all":
        return None
    if m in LEM_METHOD_OPTIONS:
        return None
    return (f"LEM method is {m!r}, which is not a method this package implements. "
            f"Choose one of: {', '.join(LEM_METHOD_OPTIONS)}, or 'all' "
            f"(Run LEM; main D14).")


# ---------------------------------------------------------------------------
# Family: surface family and method compatibility
# ---------------------------------------------------------------------------

@rule("surface.none_defined", ERROR, ("lem",), capability="analysis",
      summary="A limit-equilibrium run needs a circular or non-circular surface.",
      # Two repairs, because an empty surface sheet has two right answers and the
      # model picks between them: circles where the slope's own geometry controls
      # the mechanism, a surface tracked along the weak zone where a seam does --
      # a shape no circle can take. The second stands down unless the generator
      # picks a zone on its own; see remedies._propose_noncirc.
      remedy=("generate_starting_circles", "generate_noncircular_surface"))
def _surface_none(ctx):
    if ctx.has_circles or ctx.has_non_circ:
        return None
    if ctx.surface_supplied:
        # The run brought its own surface (see _Ctx.surface_supplied). An empty
        # sheet is then not an absence -- it is simply not where this surface came
        # from -- and refusing here would block the caller that needs the sheet
        # least.
        return None
    return ("This model defines no failure surface: it carries neither a circle "
            "nor a non-circular surface. A limit-equilibrium run needs at least one "
            "starting surface (Circles, circles sheet; Non-circular surface, "
            "non-circ sheet).")


@rule("surface.method_requires_circle", ERROR, ("lem",), capability="lem_method",
      summary="OMS and Bishop take moments about a circle center (circular family only).")
def _method_requires_circle(ctx):
    m = ctx.method
    if m not in CIRCULAR_ONLY_METHODS:
        return None
    family = ctx.effective_surface_family
    if family != "noncircular":
        return None
    return _circular_only_message(m, search=ctx.is_search)


@rule("surface.family_ambiguous", WARNING, ("lem",),
      "A deck carrying both surface families needs the run to say which one.")
def _family_ambiguous(ctx):
    if not (ctx.has_circles and ctx.has_non_circ):
        return None
    if ctx.surface_family is not None:
        return None
    return ("This model defines both a circular and a non-circular surface, and "
            "the run did not state which to analyse. The circular surface will be "
            "used; the non-circular surface will be ignored -- set Surface family "
            "to say which (Circles, circles sheet; Non-circular surface, non-circ "
            "sheet; the main sheet's Surface family row).")


# ---------------------------------------------------------------------------
# Family: the model domain, and the surfaces that have to fit inside it
#
# The domain is the one piece of geometry no analysis can do without, and it is
# DERIVED -- nobody types it -- so a defect in it is invisible on the sheet the
# user is looking at. xslope_reliability.xlsx shipped for weeks with a domain ring
# that retraced its own first edge, because Max Depth had been left at the
# elevation of the toe: the base of the domain ran back along the ground surface,
# leaving zero depth left of the toe. Nothing said so until generate_slices failed
# on a shapely error naming no field.
# ---------------------------------------------------------------------------

@rule("domain.degenerate_ring", ERROR, ("*",), capability="analysis",
      summary="The model domain must be a simple polygon: no retraced edges, "
              "no zero-area lobes.")
def _domain_degenerate(ctx):
    # The test is on the DOMAIN, not on the zones one at a time. A single zone
    # touching itself at a point is ordinary and correct -- a layer that wedges out
    # is drawn as a zone pinched to zero thickness where two profile lines meet,
    # which shapely calls a ring self-intersection on that zone and which dozens of
    # locked-green corpus files carry. What none of them carries is a defect that
    # SURVIVES the union: unary_union of sound zones is always a sound polygon, so
    # a domain that is not is the one condition that leaves the model unusable.
    defect = _ring_defect(ctx.domain) if ctx.domain is not None else None
    if defect is None:
        return None
    kind, pt = defect
    if kind == "zero_area":
        what = "encloses no area at all: its boundary has no inside"
    else:
        where = f" at ({_fmt(pt[0])}, {_fmt(pt[1])})" if pt else ""
        what = f"crosses or retraces itself{where}"
    if ctx.sd.get("profile_lines"):
        fix = ("This is what ground running AT the Max depth elevation produces -- "
               "the base of the model runs back along the ground surface itself, "
               "leaving no depth beneath it. Two cures, and which is right depends "
               "on the problem: where soil genuinely continues below, lower Max "
               "depth below the section (Profile lines; profile sheet B2); where "
               "the base is a stated rigid foundation at the toe elevation, keep "
               "Max depth there and END the ground surface at the toe instead -- "
               "soil of zero thickness beyond the toe is not modeled, and lowering "
               "the base would invent depth the problem does not describe.")
    else:
        fix = ("Redraw the zones' base so it stays below the ground it supports, and "
               "so no zone boundary retraces another (Polygons; polygon sheet).")
    return (f"The model domain -- the union of the material zones -- {what}, so it "
            f"is not a simple polygon. Nothing downstream can use it: slice "
            f"generation, meshing and the search all fail on it with a geometry "
            f"error naming no field. {fix}")


@rule("surface.circle_below_domain_floor", WARNING, ("lem",),
      summary="A circle whose Depth is below the domain floor must still cut a "
              "surface inside it; where that circle IS the run's surface, it is "
              "an error.")
def _circle_below_domain_floor(ctx):
    # `Depth` is the elevation of the circle's LOWEST POINT -- its nadir at x = Xo
    # -- and a nadir below the domain floor is NOT by itself an error: a skimming
    # circle has a very large radius whose arc runs nearly parallel to the face, so
    # its nadir can sit far below the model while the surface it cuts stays
    # entirely inside (search.py, since the depth clamp was removed). What decides
    # the matter is whether the circle slices, so that is what is asked -- plain,
    # and truncated at the floor as a composite surface. Only a circle that can do
    # neither is reported, and that one used to be silently clamped up to the floor
    # and solved as a different circle.
    #
    # SEVERITY is decided per finding, because what the same defect COSTS depends
    # on which run consumes the circle, and the two answers are far apart. The
    # question was settled by measurement rather than by reading:
    #
    #   A SEARCH survives it. The circles sheet is a set of SEEDS, and the
    #   refinement moves the center and the tangent depth away from them -- so a
    #   deck whose every circle is 200 below the floor still converges on the same
    #   critical surface as the sound deck (xslope_dam.xlsx, bishop: 1.8147 either
    #   way). Refusing that run would refuse a run that works, so a search keeps the
    #   WARNING: the seed is lost, the answer is not.
    #
    #   A SINGLE-SURFACE run does not survive it. That path takes circles[0] and
    #   nothing else (search.run_lem_analysis), so a dead first circle is the whole
    #   analysis. Today it arrives as a slicing failure naming no field, several
    #   layers below the sheet the user would fix -- which is precisely what this
    #   gate exists to pre-empt, so that one is an ERROR.
    #
    # WHAT AN ERROR REACHES, precisely: the gates whose surface comes from the
    # SHEET -- the Run dialogs, and Studio's post-edit input checks. A caller that
    # hands generate_slices its own circle sets surface_supplied and is not asked
    # at all (first line below), so run_lem_analysis's single-surface branch still
    # reaches slicing on its own. Gating that path is a separate change; the
    # severity here does not claim to have made it.
    #
    # This is the edit-cascade made visible: Max depth is raised or the base
    # redrawn, the circle that was tangent to the old base now sits under the new
    # one, and nothing about that edit says so on its own.
    if ctx.surface_supplied:
        return None
    if ctx.effective_surface_family == "noncircular":
        return None
    floor = ctx.domain_floor
    if floor is None:
        return None
    deeper = ("lower Max depth (Profile lines; profile sheet B2)"
              if ctx.sd.get("profile_lines") else
              "lower the base of the zones (Polygons; polygon sheet)")
    out = []
    for i, c in enumerate(ctx.sd.get("circles") or []):
        depth = _num(c.get("Depth")) if isinstance(c, dict) else None
        if depth is None or depth >= floor:
            continue
        if ctx.circle_slices(c):
            continue
        # circles[0] IS the surface of a single-surface run; any other circle, and
        # every circle of a search, is one candidate among several. Which is also
        # why a single-surface run must not be told its "other circles" are fine:
        # it reads none of them, so the honest thing to say about circle 2 there is
        # that this run never looks at it and a search would.
        fatal = (i == 0 and not ctx.is_search)
        if fatal:
            cost = (" This run analyses that circle and no other, so it has no "
                    "failure surface to work on.")
        elif ctx.is_search:
            cost = " The search's other circles are unaffected; this one is lost."
        else:
            cost = (f" This run analyses {ctx.circle_label(0)} only and does not "
                    f"read this one, so it does not affect the answer -- but a "
                    f"search seeded from this sheet would lose it.")
        out.append((
            ERROR if fatal else WARNING,
            f"{ctx.circle_label(i)} has Depth = {_fmt(depth)}, below the bottom of "
            f"the model domain at y = {_fmt(floor)}, and the circle it defines cuts "
            f"no failure surface that stays inside the domain -- no slices can be "
            f"generated from it, either as a circle or truncated along the base, so "
            f"this circle contributes nothing to the run.{cost} Raise Depth to at or "
            f"above y = {_fmt(floor)} {_AT_CIRCLES}, or deepen the model so the "
            f"circle fits inside it -- {deeper}."))
    return out


# ---------------------------------------------------------------------------
# Family: polyline ordering (the remedy path depends on it)
# ---------------------------------------------------------------------------

@rule("order.piezo_reversed", WARNING, ("lem",),
      "A piezometric line entered right-to-left is sorted before use.",
      remedy="reverse_polyline",
      fields=("piezo",))
def _piezo_reversed(ctx):
    out = []
    for label, pts in (("Piezometric Line 1", ctx.piezo),
                       ("Piezometric Line 2", ctx.piezo2)):
        if len(pts) >= 2 and _monotonic_x(pts) == "descending":
            out.append(f"{label}: the points run right to left. The slicer sorts "
                       f"them before use, so the analysis is unaffected, but check "
                       f"the line is what you intended {_AT_PIEZO}.")
    return out


@rule("order.piezo_nonmonotonic", ERROR, ("lem",),
      "A piezometric line whose x values rise then fall is not a surface.",
      fields=("piezo",))
def _piezo_nonmonotonic(ctx):
    out = []
    for label, pts in (("Piezometric Line 1", ctx.piezo),
                       ("Piezometric Line 2", ctx.piezo2)):
        if len(pts) >= 2 and _monotonic_x(pts) == "mixed":
            out.append(f"{label}: the X values are not in order (they rise, then "
                       f"fall). Sorting them would produce a different line from the "
                       f"one drawn -- enter the points along the water surface in "
                       f"one direction {_AT_PIEZO}.")
    return out


@rule("order.dload_reversed", WARNING, ("lem", "fem"),
      "A distributed-load block entered right-to-left is repaired silently at load.",
      remedy="reverse_polyline")
def _dload_reversed(ctx):
    out = []
    for stage in (1, 2):
        for label, pts in ctx.dload_blocks(stage):
            if len(pts) >= 2 and _monotonic_x(pts) == "descending":
                out.append(f"{label}: the points run right to left. They are "
                           f"sorted before use, but check the block is what you "
                           f"intended {_at_dloads(stage)}.")
    return out


@rule("order.dload_nonmonotonic", ERROR, ("lem", "fem"),
      "A distributed-load block whose x values rise then fall cannot be interpolated.")
def _dload_nonmonotonic(ctx):
    out = []
    for stage in (1, 2):
        for label, pts in ctx.dload_blocks(stage):
            if len(pts) >= 2 and _monotonic_x(pts) == "mixed":
                out.append(f"{label}: the X values are not in order (they rise, "
                           f"then fall). Enter the points along the loaded surface "
                           f"in one direction {_at_dloads(stage)}.")
    return out


# ---------------------------------------------------------------------------
# Family: unit-system plausibility (migrated from units.units_check)
# ---------------------------------------------------------------------------

def _declared_system(ctx):
    from .units import normalize_unit_system
    return normalize_unit_system(ctx.sd.get("unit_system"))


@rule("units.gamma_water_off_band", WARNING, ("*",),
      "The unit weight of water is the one quantity physics pins to a system.")
def _units_gamma_water(ctx):
    from .units import GAMMA_W, _GAMMA_W_DIVERGENCE
    system = _declared_system(ctx)
    if system is None:
        return None
    gw = _num(ctx.sd.get("gamma_water"))
    if gw is None:
        return None
    canonical = GAMMA_W[system]
    if abs(gw - canonical) / canonical <= _GAMMA_W_DIVERGENCE:
        return None
    return (f"Units declares {system!r}, whose unit weight of water is "
            f"{canonical:g}, but Unit weight of water carries {gw:g} "
            f"({abs(gw - canonical) / canonical * 100:.1f}% off). If that is a "
            f"deliberate seawater or brine value this is fine; otherwise check the "
            f"unit system -- xslope never converts, so a mislabelled system means "
            f"every number is read in the wrong one (Global parameters; main D8 and "
            f"D10).")


@rule("units.material_gamma_off_band", WARNING, ("*",),
      "Soil unit weights sit in far-apart bands in the two systems.",
      fields=("gamma", "gamma_sat"))
def _units_material_gamma(ctx):
    from .units import infer_unit_system
    system = _declared_system(ctx)
    if system is None:
        return None
    inferred = infer_unit_system(ctx.materials)
    if inferred is None or inferred == system:
        return None
    return (f"Units declares {system!r}, but the material unit weights look "
            f"{inferred!r} (SI soils are about 15-25 kN/m3, Imperial soils about "
            f"90-160 pcf). xslope does not convert, so a mislabelled system means "
            f"the numbers are being read in the wrong one (Global parameters, main "
            f"D8; Materials table, mat sheet).")


@rule("units.signals_disagree", WARNING, ("*",),
      "On an undeclared model, the water and soil unit weights should agree.")
def _units_signals_disagree(ctx):
    from .units import infer_unit_system, infer_system_from_gamma_water
    if _declared_system(ctx) is not None:
        return None            # the two rules above cover a declared model
    from_water = infer_system_from_gamma_water(ctx.sd.get("gamma_water"))
    from_soil = infer_unit_system(ctx.materials)
    if from_water is None or from_soil is None or from_water == from_soil:
        return None
    return (f"Units is blank, so this model declares no unit system. The unit "
            f"weight of water ({_fmt(ctx.sd.get('gamma_water'))}) looks "
            f"{from_water!r} while the material unit weights look {from_soil!r}. Two "
            f"independent signals disagree, so one set of numbers is in the wrong "
            f"system -- xslope never converts (Global parameters, main D8; Materials "
            f"table, mat sheet).")


# The Young's-modulus MAGNITUDE signal is deliberately absent from this family.
#
# ``units._E_BOUNDS`` cannot fire for the case it was written for: its SI and
# Imperial bands overlap on [20 000, 15 000 000], which is the whole plausible range
# of both, so the classic wrong-system E draws no warning under either declaration.
# The obvious repair -- make the bands disjoint -- does not survive contact with real
# models: a rock modulus of 5e6 kPa is legitimate SI and sits inside any Imperial
# band wide enough to hold a soft soil at 5e4 psf. The two systems' plausible moduli
# genuinely overlap, so no pair of disjoint bands exists that does not flag correct
# models. The signal is therefore dropped rather than repaired, and the unit-system
# question is carried by the two signals that DO separate: the unit weight of water
# (near-decisive) and the soil unit weights (strong).
#
# What survives is the other half of the defect: ``units_check`` skips ``E <= 0``
# entirely, so a blank modulus passes silently. That is a real fault, but it is a
# finite-element input fault rather than a unit-system one, and it gets its own rule.


@rule("mat.E_unusable", ERROR, ("fem",),
      "A blank or non-positive Young's modulus reaches the solver as a singular matrix.")
def _mat_E_unusable(ctx):
    # Every row, not just the referenced ones: build_fem_data demands E of the whole
    # material table (fem.py:772-784) and refuses the run on the first bad one.
    for i, m in ctx.fem_materials():
        raw = m.get("E")
        if raw is None or (isinstance(raw, str) and not str(raw).strip()):
            yield (f"{ctx.mat_label(i)} has no Young's modulus: E is blank. The "
                   f"finite element engine needs a positive modulus for every "
                   f"material it meshes {_AT_MAT}.")
            continue
        E = _num(raw)
        if E is None or E <= 0:
            yield (f"{ctx.mat_label(i)} has an unusable Young's modulus: E is "
                   f"{_fmt(raw)}. The finite element engine needs a positive modulus "
                   f"for every material it meshes; a zero or non-finite E reaches "
                   f"the solver as 'Factor is exactly singular' {_AT_MAT}.")


# ---------------------------------------------------------------------------
# Family: mesh and stored seepage fields
# ---------------------------------------------------------------------------

@rule("mesh.missing", ERROR, ("seep", "fem"), capability="analysis",
      availability=True,
      summary="A seepage or finite element run needs a mesh.")
def _mesh_missing(ctx):
    if ctx.mesh is not None:
        return None
    return ("This model carries no finite element mesh, and a seepage or finite "
            "element run has nothing to solve on without it. Build the mesh first "
            "-- Studio's Build mesh dialog, or a '{base}_mesh.json' file beside the "
            ".xlsx (Build mesh).")


@rule("mesh.element_type_unsupported", ERROR, ("seep",),
      "The seepage solver supports element type codes 3, 4, 6, 8 and 9.")
def _mesh_element_type(ctx):
    mesh = ctx.mesh
    if mesh is None:
        return None
    try:
        codes = sorted({int(v) for v in mesh.get("element_types", [])})
    except (TypeError, ValueError, AttributeError):
        return None
    bad = [c for c in codes if c not in (3, 4, 6, 8, 9)]
    if not bad:
        return None
    return (f"The mesh contains element type code(s) {', '.join(map(str, bad))}, "
            f"which the seepage solver does not support (supported: 3, 4, 6, 8, 9). "
            f"Rebuild the mesh with a supported element type {_AT_MESH}.")


@rule("mesh.element_type_linear", WARNING, ("fem",),
      "The stress solvers must run on quadratic elements; linear ones lock.")
def _mesh_element_type_linear(ctx):
    """Warn a finite element or strength-reduction run off a linear mesh.

    A warning rather than an error on purpose. Solving a linear mesh is how the
    locking is demonstrated -- the element-type comparison on the FEM overview page
    is exactly that run -- and a warning says so without blocking it.

    Two ways to see a linear mesh: one is attached, or none is yet but the model
    declares its element type on ``main!D18``. Both are checked, so the warning
    reaches Studio's run gate before a mesh has been built as well as after.
    """
    mesh = ctx.mesh
    if mesh is not None:
        # An attached mesh is what the run solves, so it decides -- whatever the
        # sheet declares. A mixed mesh (quad-dominant meshes carry stray triangles)
        # is named by every linear type in it.
        try:
            codes = {int(v) for v in mesh.get("element_types", [])}
        except (TypeError, ValueError, AttributeError):
            return None
        linear = sorted(codes & {3, 4})
        if not linear:
            return None
        names = " and ".join({3: "tri3", 4: "quad4"}[c] for c in linear)
        where = f"The mesh is built from linear {names} elements"
    else:
        declared = str(ctx.sd.get("element_type") or "").strip().lower()
        if declared not in ("tri3", "quad4"):
            return None
        where = (f"Element type is {declared}, so the mesh this run builds will "
                 f"be linear")
    return (f"{where}. Linear elements lock volumetrically under the nearly "
            f"incompressible plastic flow of a Mohr-Coulomb collapse: they resist "
            f"the failure mechanism more than the soil does, so the factor of "
            f"safety comes back too high, in the unconservative direction -- +21% "
            f"(tri3) and +11% (quad4) against the reference on the Griffiths & Lane "
            f"Example 1 benchmark. Use tri6, quad8 or quad9 for a finite element or "
            f"strength-reduction run; tri3 is for seepage, where the field is scalar "
            f"and nothing can lock (Build mesh; main D18).")


@rule("mesh.material_id_out_of_range", ERROR, ("seep", "fem"),
      "A mesh element cannot reference a material the Materials table does not "
      "define.")
def _mesh_material_range(ctx):
    mesh = ctx.mesh
    if mesh is None or not ctx.materials:
        return None
    try:
        ids = {int(v) for v in mesh.get("element_materials", [])}
    except (TypeError, ValueError, AttributeError):
        return None
    n = len(ctx.materials)
    bad = sorted(i for i in ids if i < 1 or i > n)
    if not bad:
        return None
    return (f"The mesh references material ID(s) {', '.join(map(str, bad))}, but "
            f"the Materials table defines {n} material(s). The mesh was built "
            f"against a different material table -- rebuild it after editing the "
            f"materials {_AT_MESH}.")


@rule("seep_field.node_count_mismatch", ERROR, ("lem", "fem"),
      "A stored pore-pressure field must have been computed on this mesh.")
def _seep_field_node_count(ctx):
    mesh = ctx.mesh
    if mesh is None:
        return None
    try:
        n_nodes = len(mesh.get("nodes"))
    except (TypeError, AttributeError):
        return None
    out = []
    for key, sheet in (("seep_u", "*_seep.csv"),
                       ("seep_u2", "*_seep2.csv")):
        u = ctx.sd.get(key)
        if u is None:
            continue
        try:
            n_u = len(u)
        except TypeError:
            continue
        if n_u != n_nodes:
            out.append(
                f"The stored pore-pressure field (the '{sheet}' companion "
                f"file) carries {n_u} node "
                f"value(s) but the mesh has {n_nodes} node(s), so it was computed on "
                f"a different mesh. Re-run the seepage analysis on this mesh.")
    return out


@rule("seep_field.no_consumer", WARNING, ("lem", "fem"),
      "A solved seepage field no material reads produces no pore pressure.")
def _seep_field_no_consumer(ctx):
    # The twin of piezo.no_consumer: the field exists -- solved, or loaded from
    # its sidecar -- and no material has u = seep, so the run computes with
    # whatever the materials do say (none, piezo, ru) and the seepage solve is
    # never read. Silent on a model with no field to read.
    if ctx.sd.get("seep_u") is None or ctx.uses_seep_u:
        return None
    return ("This model carries a solved seepage field, but no material takes its "
            "pore pressure from it: every material's u option is none, piezo or "
            "ru, so the run never reads the seepage solution. Set u = seep on the "
            "materials that should read it, or leave them as they are if the "
            "seepage run was not meant to feed this analysis.")


@rule("seep_field.missing", ERROR, ("lem",),
      "A material with u = seep needs a solved pore-pressure field.")
def _seep_field_missing(ctx):
    if not ctx.uses_seep_u:
        return None
    # A stored field satisfies this; so does a transient march with one of its frames
    # about to be staged into ``seep_u`` (see :attr:`_Ctx.seep_frame`) -- that field
    # is not missing, it is one step away, and refusing it would refuse precisely the
    # run the transient controls exist to start. ``seep_field.transient_frame`` then
    # says which instant, as an INFO. Neither substitutes for the MESH the field is
    # read onto: without one, nothing can supply pore pressures at all.
    if ctx.mesh is not None and (ctx.sd.get("seep_u") is not None
                                 or ctx.seep_frame is not None):
        return None
    names = [ctx.mat_label(i) for i, m in ctx.strength_materials()
             if str(m.get("u", "")).strip().lower() == "seep"]
    missing = "mesh" if ctx.mesh is None else "pore-pressure field"
    return (f"{names[0]} takes pore pressure from a seepage solution (u = seep), but "
            f"this model carries no {missing}. Run the seepage analysis first, or "
            f"choose a different pore-pressure option -- every slice base in that "
            f"material would otherwise read zero pore pressure.")


@rule("seep_field.transient_frame", INFO, ("lem", "fem"),
      "u = seep against a transient march reads one named instant of it.")
def _seep_field_transient_frame(ctx):
    """Name the instant the run will read out of the attached transient solution.

    The counterpart of the clause above: where a stored field says nothing about
    WHEN, a transient one carries a whole sequence, and a factor of safety belongs to
    exactly one of them. Stated as an INFO because it is a choice being reported, not
    a defect -- and stated at all because a run whose pore pressures came from an
    unnamed instant is a number nobody can reproduce.
    """
    if not ctx.uses_seep_u or ctx.seep_frame is None:
        return None
    times = ctx.seep_frame_times
    if not times:
        return None             # a frame is coming; which one is not readable yet
    unit = str(ctx.sd.get("time_unit") or "").strip()
    unit = f" {unit}" if unit else ""
    if len(times) >= 2:
        return (f"Pore pressures come from the transient seepage solution: the "
                f"rapid-drawdown stages read t = {_fmt(times[0])}{unit} and "
                f"t = {_fmt(times[1])}{unit}. Those two frames are read into the "
                f"model when the run starts, in place of any stored field.")
    return (f"Pore pressures come from the transient seepage solution, at "
            f"t = {_fmt(times[0])}{unit}. That frame is read into the model when the "
            f"run starts, in place of any stored field.")


# ---------------------------------------------------------------------------
# Family: steady seepage -- material properties
# ---------------------------------------------------------------------------

@rule("seep.k1_nonpositive", ERROR, ("seep",),
      "The major hydraulic conductivity k1 must be greater than zero.",
      fields=("k1",))
def _seep_k1(ctx):
    for i, m in ctx.seepage_materials():
        if _pos(m.get("k1")):
            continue
        yield (f"{ctx.mat_label(i)} has no hydraulic conductivity: k1 is "
               f"{_fmt(m.get('k1'))} and must be greater than zero. Model an "
               f"impermeable zone with a small conductivity (for example 1e-6 times "
               f"the largest k), never with zero -- a zero makes the conductivity "
               f"matrix singular and destroys the flow net while the head field "
               f"still looks plausible {_AT_MAT}.")


@rule("seep.k2_nonpositive", ERROR, ("seep",),
      "The minor hydraulic conductivity k2 must be greater than zero.",
      fields=("k2",))
def _seep_k2(ctx):
    for i, m in ctx.seepage_materials():
        if _pos(m.get("k2")):
            continue
        yield (f"{ctx.mat_label(i)} has no minor hydraulic conductivity: k2 is "
               f"{_fmt(m.get('k2'))} and must be greater than zero. For an isotropic "
               f"material set k2 = k1 {_AT_MAT}.")


@rule("seep.k2_greater_than_k1", WARNING, ("seep",),
      "k1 is the MAJOR principal conductivity; k2 > k1 rotates the ellipse 90 degrees.",
      fields=("k1", "k2"))
def _seep_k2_gt_k1(ctx):
    for i, m in ctx.seepage_materials():
        k1, k2 = _num(m.get("k1")), _num(m.get("k2"))
        if k1 is None or k2 is None or k1 <= 0 or k2 <= 0:
            continue
        if k2 > k1:
            yield (f"{ctx.mat_label(i)} has k2 = {k2:g}, which is greater than "
                   f"k1 = {k1:g}. k1 is the MAJOR principal conductivity, so this "
                   f"silently rotates the conductivity ellipse by 90 degrees -- swap "
                   f"the two values, or set alpha = 90 {_AT_MAT}.")


@rule("seep.anisotropy_angle_unset", INFO, ("seep",),
      "An anisotropic material with alpha blank puts k1 along +x.",
      fields=("alpha", "k1", "k2"))
def _seep_alpha(ctx):
    for i, m in ctx.seepage_materials():
        k1, k2 = _num(m.get("k1")), _num(m.get("k2"))
        if k1 is None or k2 is None or k1 <= 0 or k2 <= 0 or k1 == k2:
            continue
        if (_num(m.get("alpha")) or 0.0) == 0.0:
            yield (f"{ctx.mat_label(i)} is anisotropic (k1 = {k1:g}, k2 = {k2:g}) "
                   f"with alpha blank, so the major conductivity is taken along +x "
                   f"{_AT_MAT}.")


@rule("seep.unsat_params_missing", ERROR, ("seep",),
      "An unconfined model needs the unsaturated parameters its unsat model uses.",
      fields=("kr0", "h0"))
def _seep_unsat_params(ctx):
    _, _, n_exit = ctx.bc_counts(1)
    _, _, n_exit2 = ctx.bc_counts(2)
    if not (n_exit or n_exit2):
        return None                # confined: the kr branch is never evaluated
    out = []
    for i, m in ctx.seepage_materials():
        model = str(m.get("unsat") or "lf").strip().lower()
        if model == "lf":
            kr0, h0 = _num(m.get("kr0")), _num(m.get("h0"))
            if kr0 is None or kr0 <= 0 or h0 is None or h0 >= 0:
                out.append(
                    f"{ctx.mat_label(i)} uses unsat = lf on an unconfined model (an "
                    f"exit face is defined), so it needs kr0 > 0 and h0 < 0; it "
                    f"carries kr0 = {_fmt(m.get('kr0'))}, "
                    f"h0 = {_fmt(m.get('h0'))} {_AT_MAT}.")
        elif model in ("vg", "gard"):
            a, n = _num(m.get("vg_a")), _num(m.get("vg_n"))
            need = "n > 1" if model == "vg" else "n > 0"
            bad = (a is None or a <= 0 or n is None
                   or (n <= 1 if model == "vg" else n <= 0))
            if bad:
                out.append(
                    f"{ctx.mat_label(i)} uses unsat = {model} on an unconfined "
                    f"model, so it needs a > 0 and {need}; it carries "
                    f"a = {_fmt(m.get('vg_a'))}, n = {_fmt(m.get('vg_n'))} "
                    f"{_AT_MAT}.")
    return out


@rule("seep.confined_unsat_unused", INFO, ("seep",),
      "With no exit face the solve is confined and the unsaturated inputs are inert.")
def _seep_confined_unsat(ctx):
    _, _, n_exit = ctx.bc_counts(1)
    _, _, n_exit2 = ctx.bc_counts(2)
    if n_exit or n_exit2:
        return None
    named = [ctx.mat_label(i) for i, m in ctx.seepage_materials()
             if str(m.get("unsat") or "lf").strip().lower() != "lf"]
    if not named:
        return None
    return (f"No exit face is defined, so this solve is confined and the "
            f"unsaturated inputs (unsat, kr0, h0, a, n) are not used -- including on "
            f"{named[0]} (Seep BC, seep bc sheet; Materials table, mat sheet).")


# ---------------------------------------------------------------------------
# Family: steady seepage -- boundary conditions
# ---------------------------------------------------------------------------

@rule("seep.no_boundary_conditions", ERROR, ("seep",), capability="analysis",
      summary="A seepage run needs at least one specified head or an exit face.")
def _seep_no_bcs(ctx):
    n_head, n_flux, n_exit = ctx.bc_counts(1)
    if n_head or n_flux or n_exit:
        return None
    return ("This model defines no seepage boundary conditions at all. A seepage "
            "analysis needs at least one specified head, or an exit face. With none "
            "the solve converges to zero head everywhere and reports success, which "
            "reaches a stability run as zero pore pressure (Seep BC; seep bc "
            "sheet).")


@rule("seep.no_dirichlet", ERROR, ("seep",),
      "A flux-only model defines head only up to an additive constant.")
def _seep_no_dirichlet(ctx):
    n_head, n_flux, n_exit = ctx.bc_counts(1)
    if n_head or n_exit or not n_flux:
        return None
    return ("This model defines specified-flux boundary conditions but no "
            "specified head and no exit face. The head is then defined only up to an "
            "additive constant and the system is singular -- add a specified head, "
            "or an exit face (Seep BC; seep bc sheet).")


@rule("seep.exit_face_without_head", ERROR, ("seep",),
      "An exit face alone has nothing driving flow through it.")
def _seep_exit_only(ctx):
    n_head, n_flux, n_exit = ctx.bc_counts(1)
    if not n_exit or n_head or n_flux:
        return None
    return ("This model defines an exit face but no specified head and no flux. "
            "Once the exit face deactivates there is no boundary left to fix the "
            "head, and the solve becomes singular. Add the reservoir or tailwater "
            "head that drives the flow (Seep BC; seep bc sheet).")


@rule("seep.no_gradient", ERROR, ("seep",),
      "Two heads of the same value, or one head alone, drive no flow.")
def _seep_no_gradient(ctx):
    # An ERROR, not a warning (Norm 2026-08-26): a set of boundary conditions
    # that is nothing but specified heads at one value poses no flow problem —
    # there is no gradient and no solution, and the unconfined sweep runs to its
    # limit looking for one. An exit face does not rescue it: it is a place
    # water may leave, not a head that drives any.
    bc = ctx.seep_bc
    heads = bc.get("specified_heads") or []
    n_flux = len(bc.get("specified_fluxes") or [])
    if n_flux or not heads:
        return None
    # An exit face IS an outlet where it lies below the water the heads set:
    # there the pressure is zero and the head equals the elevation, which is a
    # gradient (the corpus's dam-with-a-seepage-face cases). An exit face that
    # sits entirely above every specified head drains nothing and rescues
    # nothing.
    exit_face = bc.get("exit_face") or []
    finite_heads = [float(b.get("head")) for b in heads if _finite(b.get("head"))]
    if exit_face and finite_heads:
        try:
            lowest_exit = min(float(pt[1]) for pt in exit_face)
        except (TypeError, ValueError, IndexError):
            lowest_exit = None
        if lowest_exit is not None and lowest_exit < max(finite_heads) - 1e-9:
            return None
    if ctx.sd.get("tseep") is not None:
        return None            # a transient run drives the heads from its series
    if any(isinstance(b.get("head"), str) for b in heads):
        return None            # bound to a tseep series, so the value varies in time
    values = {round(float(b.get("head")), 12) for b in heads
              if _finite(b.get("head"))}
    if len(values) >= 2:
        return None
    return ("The specified heads all carry the same value, so this model has no "
            "gradient and no flow: the problem has no solution to find. Give the "
            "boundary set at least two distinct head values, or a flux boundary "
            "(Seep BC; seep bc sheet).")


@rule("seep.bc_polyline_too_short", ERROR, ("seep",),
      "A boundary-condition polyline needs at least two points to bind mesh nodes.")
def _seep_bc_short(ctx):
    bc = ctx.seep_bc
    out = []
    for n, b in enumerate(bc.get("specified_heads") or []):
        if len(_coords(b.get("coords"))) < 2:
            out.append(f"Specified head #{n + 1} has fewer than two points, so "
                       f"it binds no mesh nodes and is silently skipped "
                       f"{_AT_SEEPBC}.")
    for n, b in enumerate(bc.get("specified_fluxes") or []):
        if len(_coords(b.get("coords"))) < 2:
            out.append(f"Specified flux #{n + 1} has fewer than two points, so "
                       f"it binds no mesh nodes {_AT_SEEPBC}.")
    exit_face = bc.get("exit_face") or []
    if 0 < len(_coords(exit_face)) < 2:
        out.append(f"Exit face has fewer than two points, so it binds no mesh "
                   f"nodes and the model would silently solve as fully confined "
                   f"{_AT_SEEPBC}.")
    return out


# ---------------------------------------------------------------------------
# Family: finite element elastic properties and the tensile cap
#
# The asymmetry these rules exist to close: a blank E is caught at assembly with
# a clear message, a blank nu is NOT -- it reads as 0.0, passes the [0, 0.5)
# range check, and moves the strength-reduction factor by a third. And a blank
# t_cut is silent in both engines while granting the material unbounded tension,
# which is the transcription defect the RS2 campaign traced a 2x factor to.
# ---------------------------------------------------------------------------

@rule("mat.nu_unusable", ERROR, ("fem",),
      "Poisson's ratio nu must be in (0, 0.5); a blank one reads as 0.0.")
def _mat_nu_unusable(ctx):
    for i, m in ctx.fem_materials():
        nu = _num(m.get("nu"))
        if nu is not None and 0.0 < nu < 0.5:
            continue
        if nu is None:
            what = ("is blank or not a number. The finite element engine needs a "
                    "Poisson's ratio for every material it meshes")
        elif nu == 0.0:
            what = ("is 0 -- which is what a BLANK cell reads as. A soil with no "
                    "lateral coupling is not a soil, and the value is not inert: on "
                    "the reference model it moved the strength-reduction factor of "
                    "safety by a third against nu = 0.3")
        else:
            what = (f"is {nu:g}, outside the admissible range [0, 0.5). The finite "
                    f"element engine refuses it at assembly")
        yield (f"{ctx.mat_label(i)} has an unusable Poisson's ratio: nu {what}. "
               f"Enter a value -- typically 0.2-0.45 for soil, 0.15-0.3 for rock "
               f"{_AT_MAT}.")


@rule("mat.tensile_cap_missing", WARNING, ("fem",),
      "A blank tensile strength t_cut grants the material unbounded tension.")
def _mat_tensile_cap_missing(ctx):
    for i, m in ctx.fem_materials():
        opt = str(m.get("option") or "").strip().lower()
        if opt in ("", "elastic"):
            continue                  # no strength model, or t_cut is meaningless
        raw = m.get("t_cut")
        if raw is not None and _num(raw) is not None:
            continue
        phi = _num(m.get("phi")) or 0.0
        c = _num(m.get("c")) or 0.0
        if phi <= 0:
            what = ("without limit -- with φ = 0 the Mohr-Coulomb cone has no apex, "
                    "so nothing bounds the tensile stress the material can carry")
        else:
            apex = c / math.tan(math.radians(phi)) if phi < 90 else 0.0
            what = (f"up to the Mohr-Coulomb cone apex, c/tan(φ) = {apex:.4g}, which "
                    f"may be far above the material's real tensile strength")
        yield (f"{ctx.mat_label(i)} has no tensile cutoff: a blank t_cut allows "
               f"tension {what}. Enter t_cut, or 0 for a no-tension soil. A "
               f"dropped cap is not visible in any result: it raises the "
               f"strength-reduction factor of safety silently {_AT_MAT}.")


@rule("main.tension_srf_unset", INFO, ("fem",),
      "A blank Tension SRF (main!D17) means the tensile cap IS reduced by the "
      "trial factor.")
def _tension_srf_unset(ctx):
    if ctx.sd.get("tension_srf") is not None:
        return None
    capped = [ctx.mat_label(i) for i, m in ctx.fem_materials()
              if _num(m.get("t_cut")) is not None]
    if not capped:
        return None
    return (f"Tension SRF (FEM) is blank, so the engine default applies: the "
            f"tensile cap is divided by the trial strength-reduction factor "
            f"alongside c and tan(φ). {len(capped)} material(s) carry a t_cut. Set "
            f"it to NO to hold the cap at its entered value instead "
            f"{_at_global('D17')}.")


# ---------------------------------------------------------------------------
# Family: finite element initial stress and strength-reduction zones
# ---------------------------------------------------------------------------

@rule("fem.k0_without_zone_geometry", ERROR, ("fem",),
      "K0 initialization integrates the overburden through the material zones.")
def _k0_without_zones(ctx):
    if _num(ctx.sd.get("k0")) is None:
        return None
    from .mesh import get_material_polygons
    try:
        zones = get_material_polygons(ctx.sd)
    except Exception:
        zones = []
    if zones:
        return None
    return ("K0 initial stress (FEM) is set, but this model carries no "
            "material-zone geometry for the engine to integrate the overburden "
            "through, so the initial stress state cannot be built. Supply the "
            "profile lines or polygons, or leave K0 blank for the gravity turn-on "
            "initialization (Global parameters; main D16).")


@rule("fem.k0_pre_equilibration", INFO, ("fem",),
      "A K0 run on non-level ground adds a full-strength pre-equilibration solve.")
def _k0_pre_equilibration(ctx):
    k0 = _num(ctx.sd.get("k0"))
    if k0 is None:
        return None
    h = ctx.slope_height
    if h is None or h <= 1e-9:
        return None            # level ground: the K0 field is already in equilibrium
    return (f"K0 initial stress (FEM) = {k0:g} on ground that is not level "
            f"(relief {h:.4g}). The K0 stress field does not equilibrate itself "
            f"there, so the engine runs one extra full-strength solve first and "
            f"hands the settled state to every trial. The factor of safety moves "
            f"with K0 without any further notice {_at_global('D16')}.")


@rule("ssr.zone_contains_no_elements", ERROR, ("fem",),
      "An SSR zone that catches no mesh element constrains nothing it claims to.")
def _ssr_zone_no_elements(ctx):
    zones = ctx.sd.get("ssr_zones") or []
    if not zones:
        return None
    cen = ctx.element_centroids()
    if not cen:
        return None                    # no mesh yet: mesh.missing reports that
    from shapely.geometry import Polygon, Point
    labels = {"reduce": "SSR reduce", "hold": "SSR hold",
              "hold_elastic": "SSR elastic"}
    out = []
    for n, z in enumerate(zones):
        try:
            poly = Polygon([(float(x), float(y)) for x, y in z.get("polygon") or []])
        except (TypeError, ValueError):
            continue
        if not poly.is_valid:
            poly = poly.buffer(0)
        if poly.is_empty:
            continue
        inside = sum(1 for x, y in cen if poly.contains(Point(x, y)))
        kind = str(z.get("kind") or "")
        name = z.get("label") or labels.get(kind, kind)
        where = f"{name} zone #{n + 1}"
        if inside == 0:
            if kind == "reduce":
                out.append(
                    f"{where} contains no mesh elements, so strength reduction would "
                    f"apply nowhere: every element stays at full strength, the model "
                    f"never fails, and the run blames F_max instead. Move the zone "
                    f"onto the mesh, or delete it to reduce the whole domain "
                    f"{_AT_POLYGON}.")
            else:
                out.append((WARNING,
                            f"{where} contains no mesh elements, so it excludes "
                            f"nothing. The constraint it was drawn for is not in "
                            f"the run {_AT_POLYGON}."))
        elif kind == "reduce" and inside == len(cen):
            out.append((INFO,
                        f"{where} contains every element in the mesh, so it is "
                        f"equivalent to no reduce zone at all {_AT_POLYGON}."))
    return out


@rule("mesh.zone_size_not_finer", WARNING, ("fem", "seep"),
      "A zone element size at or above the global target refines nothing.")
def _zone_size_not_finer(ctx):
    target = _num(ctx.sd.get("target_size"))
    if target is None or target <= 0:
        return None
    out = []
    for n, z in enumerate(ctx.sd.get("refine_zones") or []):
        s = _num(z.get("size"))
        if s is not None and s >= target:
            out.append(f"Refine zone #{n + 1} declares Size = {s:g}, which is "
                       f"not finer than the global target element size ({target:g}). "
                       f"The zone is inert -- a refine polygon's only effect is a "
                       f"SMALLER local element size (Polygons, polygon sheet; Build "
                       f"mesh, main D19).")
    for n, p in enumerate(ctx.sd.get("polygons") or []):
        s = _num(p.get("size") if isinstance(p, dict) else None)
        if s is not None and s >= target:
            out.append(f"Polygon #{n + 1} declares Size = {s:g}, which is not "
                       f"finer than the global target element size ({target:g}), so "
                       f"it changes nothing (Polygons, polygon sheet; Build mesh, "
                       f"main D19).")
    return out


#: Element rows a material zone needs across its local width before a shear band can
#: form through it. Three is the mesher's own definition of a resolved thin zone
#: (``mesh._REFINE_THIN_MIN_ELEMS``), so the advisory and the mesher agree on what
#: "too thin" means.
_THIN_ZONE_MIN_ROWS = 3.0
#: Rows the two repairs deliver, and so the number the suggested Size is derived
#: from -- one more than the minimum, which is what the Build-mesh dialog's toggle
#: asks for (``studio.runners.THIN_ZONE_ELEMS``).
_THIN_ZONE_WANT_ROWS = 4.0


def _element_sizes(ctx):
    """One length per mesh element -- the square root of its area, which is the same
    length scale the mesher's target size names -- aligned with
    :meth:`_Ctx.element_centroids`. ``[]`` when there is no usable mesh.

    Corner nodes only: a quadratic element's mid-side nodes lie on the same edges
    and would not change the area.
    """
    def _build():
        mesh = ctx.mesh
        if mesh is None:
            return []
        try:
            nodes = mesh["nodes"]
            elements = mesh["elements"]
            types = mesh["element_types"]
        except (KeyError, TypeError):
            return []
        out = []
        for i, elem in enumerate(elements):
            try:
                t = int(types[i])
                corners = [nodes[int(k)]
                           for k in list(elem)[:(3 if t in (3, 6) else 4)]]
            except (IndexError, TypeError, ValueError):
                out.append(0.0)
                continue
            n = len(corners)
            if n < 3:
                out.append(0.0)
                continue
            try:
                area = 0.5 * abs(sum(
                    float(corners[j][0]) * float(corners[(j + 1) % n][1])
                    - float(corners[(j + 1) % n][0]) * float(corners[j][1])
                    for j in range(n)))
            except (TypeError, ValueError):
                area = 0.0
            out.append(area ** 0.5 if area > 0 else 0.0)
        return out
    return ctx._c("element_sizes", _build)


def _median(values):
    s = sorted(v for v in values if v > 0)
    return s[len(s) // 2] if s else None


def _thin_zone_geometry(ctx):
    """``[(index, label, width, declared_size, ring), ...]`` for every material zone
    that is thin EVERYWHERE, at the element size this model is meshed at.

    That size is main!D19 where the file states one, and otherwise the median element
    size of the attached mesh -- the size the mesh was in fact built at, which is the
    honest stand-in for a blank cell and is what makes the rule live on a model whose
    mesh came from the Build-mesh dialog rather than from D19.

    Restricted to whole-thin zones -- a band, a seam, a liner -- because those are
    the ones whose local width IS the zone. A local pinch inside an otherwise thick
    polygon is a corner of the geometry rather than a layer the mechanism has to run
    through, and warning about every acute corner in a section would bury the case
    this rule exists for.

    A zone whose material is ``option = elastic`` is left out entirely. The whole
    premise is that a zone too coarse to mesh cannot develop a shear band through
    it -- and an elastic material cannot yield at ANY element size, so there is no
    band to lose and no factor of safety the mesh could inflate. Only elasticity
    earns this: a thin zone that is merely STRONG still fails somewhere in the end,
    which is exactly the case the rule exists to catch.
    """
    def _build():
        target = _num(ctx.sd.get("target_size"))
        if target is None or target <= 0:
            target = _median(_element_sizes(ctx))
        if target is None or target <= 0:
            return []
        from .mesh import detect_thin_zones, get_material_polygons
        try:
            zones = get_material_polygons(ctx.sd)
        except Exception:
            return []
        rings, declared, mat_ids = [], [], []
        for z in zones:
            rings.append(list(z.get("coords") or []) if isinstance(z, dict) else list(z))
            declared.append(_num(z.get("size")) if isinstance(z, dict) else None)
            mat_ids.append(z.get("mat_id") if isinstance(z, dict) else None)
        # Measured on the MATERIAL, not on each polygon: a layer the section geometry
        # splits into pieces (a bench, a step in the ground surface) is not thin just
        # because its pieces are, and warning that it is would send the reader looking
        # for a band that is not there. Same grouping the mesher's refinement uses, so
        # the warning and the refinement always agree about which zones are thin.
        group_ids = ([m if m is not None else ("_", i) for i, m in enumerate(mat_ids)]
                     if any(m is not None for m in mat_ids) else None)
        try:
            detected = detect_thin_zones(rings, target, group_ids=group_ids)
        except Exception:
            return []
        on_sheet = bool(ctx.sd.get("polygons"))
        out = []
        for d in detected:
            if d.get("kind") != "whole":
                continue
            i = d.get("poly_index")
            if not isinstance(i, int) or not (0 <= i < len(rings)):
                continue
            mat = None
            try:
                mat = ctx.materials[int(mat_ids[i])]
            except (IndexError, TypeError, ValueError):
                pass
            if not isinstance(mat, dict):
                mat = {}
            if str(mat.get("option") or "").strip().lower() == "elastic":
                continue
            name = mat.get("name")
            where = (f"Polygon #{i + 1}" if on_sheet
                     else f"Material zone #{i + 1}")
            label = f"{where} ('{name}')" if name else where
            # A layer stored as several polygons is reported, measured and mitigated
            # as the one zone it is: the merged ring is where the mesh's element rows
            # are counted, and a Size on ANY member is the user having sized it.
            members = d.get("poly_indices") or [i]
            ring = list(d.get("ring") or rings[i])
            decl = next((declared[j] for j in members
                         if 0 <= j < len(declared) and declared[j] is not None), None)
            out.append((i, label, float(d["width"]), decl, ring))
        return out
    return ctx._c("thin_zones", _build)


def _mesh_rows_across(ctx, ring, width):
    """Element rows the ATTACHED mesh puts across a zone -- the zone's width over the
    median size of the elements whose centroid falls inside it, or None when there is
    no mesh or it places no element there.

    Measured rather than inferred, because the two things that would otherwise have
    to be inferred -- an adequate Size, a mesh built with thin-zone refinement -- both
    show up here as elements that are small enough, and a mesh that resolves the zone
    by any route at all must not be warned about.
    """
    sizes = _element_sizes(ctx)
    centroids = ctx.element_centroids()
    if not sizes or len(sizes) != len(centroids):
        return None
    from shapely.geometry import Point, Polygon
    from shapely.prepared import prep
    try:
        poly = Polygon([(float(x), float(y)) for x, y in ring])
        if not poly.is_valid:
            poly = poly.buffer(0)
        if poly.is_empty or poly.area <= 0:
            return None
    except (TypeError, ValueError):
        return None
    xmin, ymin, xmax, ymax = poly.bounds
    inside = prep(poly)
    local = [sizes[i] for i, (x, y) in enumerate(centroids)
             if xmin <= x <= xmax and ymin <= y <= ymax
             and inside.contains(Point(x, y))]
    median = _median(local)
    if median is None:
        return None
    return float(width) / float(median)


@rule("mesh.thin_zone_unresolved", WARNING, ("fem",),
      "A material zone too thin to mesh cannot carry a shear band through it.")
def _thin_zone_unresolved(ctx):
    thin = _thin_zone_geometry(ctx)
    if not thin:
        return None
    target = _num(ctx.sd.get("target_size"))
    out = []
    for _i, label, width, declared, ring in thin:
        want = width / _THIN_ZONE_WANT_ROWS
        fix = (f"Enter Size = {want:.3g} for the zone (Polygons; polygon sheet), "
               f"or rebuild with the Build mesh dialog's 'Refine thin zones' checked "
               f"-- either one meshes it at about {_THIN_ZONE_WANT_ROWS:g} element "
               f"rows.")
        rows = _mesh_rows_across(ctx, ring, width)
        if rows is not None:
            # A mesh is attached: what it actually did settles the question, and
            # every mitigation there could be shows up in it.
            if rows >= _THIN_ZONE_MIN_ROWS:
                continue
            out.append(
                f"{label} is about {width:.3g} wide and the mesh puts {rows:.1f} "
                f"element rows across it. A zone this thin cannot develop a shear "
                f"band, so the run does not fail on it -- it returns a factor of "
                f"safety that is too high, with nothing in the result to say the "
                f"mesh was the reason. {fix}")
            continue
        # No mesh to measure: the declared inputs are all there is to go on, and the
        # width above is a proxy (twice the area over the perimeter), sharp enough
        # to find a zone nobody sized and not sharp enough to overrule a number
        # somebody entered. So a Size that actually refines the zone IS the
        # mitigation here, taken at its word; a Size at or above the global target
        # refines nothing and mesh.zone_size_not_finer already says so.
        if target is None or target <= 0:
            continue
        if declared is not None and 0 < declared < target:
            continue
        rows = width / target
        if rows >= _THIN_ZONE_MIN_ROWS:
            continue
        out.append(
            f"{label} is about {width:.3g} wide, which is {rows:.1f} element rows at "
            f"the global target element size ({target:g}, main D19), and the zone "
            f"declares no Size of its own. A zone this thin cannot develop a "
            f"shear band, so the run does not fail on it -- it returns a factor of "
            f"safety that is too high, with nothing in the result to say the mesh "
            f"was the reason. {fix}")
    return out


# ---------------------------------------------------------------------------
# Family: transient seepage
#
# The storage rules are where the silence is worst: the loader refuses a MISSING
# Ss or Sy and says so clearly, but an explicit ZERO passes both guards and runs.
# With Sy = 0 the phreatic surface tracks a falling pool instantly, which is the
# least conservative drawdown pore pressure the model can produce -- and nothing
# says a word.
# ---------------------------------------------------------------------------

@rule("tseep.time_unit_missing", ERROR, ("tseep",),
      "A transient run makes the time unit load-bearing, so it must be declared.")
def _tseep_time_unit(ctx):
    if ctx.tseep is None:
        return None
    unit = ctx.sd.get("time_unit")
    if unit is not None and str(unit).strip():
        return None
    return ("Time is blank, so this model declares no time unit, and a transient "
            "seepage run is defined. Conductivity is length/time, storage is "
            "1/length and every transient time is in that unit, so the run cannot be "
            "read without it. It is never guessed: a wrong time label is worse than "
            "none (Global parameters; main D9).")


@rule("tseep.storage_nonpositive", ERROR, ("tseep",),
      "Specific storage Ss (and specific yield Sy on an unconfined model) must "
      "be greater than zero.")
def _tseep_storage(ctx):
    if ctx.tseep is None:
        return None
    unconfined = ctx.unconfined
    for i, m in ctx.seepage_materials():
        if not _pos(m.get("Ss")):
            yield (f"{ctx.mat_label(i)} has no specific storage: Ss is "
                   f"{_fmt(m.get('Ss'))} on a transient run. A zero removes the "
                   f"entire storage term where the material is saturated, and with "
                   f"Sy zero as well the residual storage floor collapses too and "
                   f"the time march no longer terminates {_AT_MAT}.")
        if unconfined and not _pos(m.get("Sy")):
            yield (f"{ctx.mat_label(i)} has no specific yield: Sy is "
                   f"{_fmt(m.get('Sy'))}, and this model is unconfined (an exit face "
                   f"is drawn). With no drainable porosity the water table follows "
                   f"the boundary INSTANTLY, so a falling pool produces the least "
                   f"conservative pore pressures the model can return -- silently "
                   f"{_AT_MAT}.")


@rule("tseep.retention_band_substituted", WARNING, ("tseep",),
      "h0 = 0 on a linear-front or Gardner material invents a one-unit drainage band.")
def _tseep_retention_band(ctx):
    if ctx.tseep is None:
        return None
    for i, m in ctx.seepage_materials():
        model = str(m.get("unsat") or "lf").strip().lower()
        if model not in ("lf", "gard"):
            continue
        h0 = _num(m.get("h0"))
        if h0 is not None and h0 != 0.0:
            continue
        yield (f"{ctx.mat_label(i)} has h0 = {_fmt(m.get('h0'))} with "
               f"unsat = {model}. On a transient run h0 is the pressure band the "
               f"soil drains over, and a zero is silently replaced by -1 -- one "
               f"length unit of the model's own system, invented rather than "
               f"entered. Enter the negative pressure head at which the material "
               f"reaches its residual conductivity {_AT_MAT}.")


@rule("tseep.duration_invalid", ERROR, ("tseep",),
      "Duration must be present and greater than zero.")
def _tseep_duration(ctx):
    ts = ctx.tseep
    if ts is None:
        return None
    if _pos(ts.get("duration")):
        return None
    return (f"Duration is {_fmt(ts.get('duration'))}. A transient run needs a "
            f"positive duration -- it is the length of the march, and there is "
            f"deliberately no default {_AT_TSEEP}.")


@rule("tseep.save_interval", INFO, ("tseep",),
      "A blank Save interval defaults to duration/50; one above the duration "
      "saves nothing between.")
def _tseep_save_interval(ctx):
    ts = ctx.tseep
    if ts is None:
        return None
    dur = _num(ts.get("duration"))
    si = _num(ts.get("save_interval"))
    if si is None or si <= 0:
        if dur is None:
            return None
        return (f"Save interval is blank, so the run saves every {dur / 50:.6g} "
                f"(duration/50, 51 frames). This affects only which frames are "
                f"stored -- the time step is chosen adaptively and a coarse save "
                f"interval cannot coarsen the answer {_AT_TSEEP}.")
    if dur is not None and si > dur:
        return (WARNING,
                f"Save interval = {si:g} is longer than the Duration ({dur:g}), so "
                f"only the endpoint and the scheduled breakpoints are saved. A "
                f"stability run asking for an intermediate time would have to "
                f"re-march to reach it {_AT_TSEEP}.")
    return None


@rule("tseep.stage_times", ERROR, ("tseep", "rapid"),
      "Rapid-drawdown staging needs Stage 1 time and Stage 2 time, in that "
      "order.")
def _tseep_stage_times(ctx):
    ts = ctx.tseep
    if ts is None:
        # A steady rapid-drawdown model has no stage TIMES: its stage-2 state is the
        # second boundary set, read directly. Only a transient model interpolates.
        return None
    s1, s2 = _num(ts.get("stage_1")), _num(ts.get("stage_2"))
    rapid = "rapid" in ctx.analyses
    if s1 is None and s2 is None:
        if not rapid:
            return None            # a plain transient march needs no stage times
        return ("This is a rapid-drawdown run on a transient seepage model, so it "
                "needs Stage 1 time (the full pool) and Stage 2 time (the drawn-down "
                "state). Neither is set (Transient seepage; tseep sheet).")
    if s1 is None or s2 is None:
        missing = "Stage 1 time" if s1 is None else "Stage 2 time"
        return (f"{missing} is blank. Rapid-drawdown staging needs both stage times "
                f"-- one alone says nothing about which two states the drawdown runs "
                f"between {_AT_TSEEP}.")
    if s1 >= s2:
        return (f"Stage 1 time = {s1:g} is not earlier than Stage 2 time = {s2:g}. "
                f"The full-pool state must precede the drawn-down state. Today the "
                f"march runs to completion and reports success; the ordering only "
                f"surfaces if the staging is invoked afterwards {_AT_TSEEP}.")
    dur = _num(ts.get("duration"))
    if dur is not None and s2 > dur:
        return (WARNING,
                f"Stage 2 time = {s2:g} is beyond the Duration ({dur:g}), so that "
                f"state is never reached and no frame is saved there {_AT_TSEEP}.")
    return None


@rule("tseep.save_times_beyond_duration", WARNING, ("tseep",),
      "A scheduled save time past the end of the run is never reached.")
def _tseep_save_times(ctx):
    ts = ctx.tseep
    if ts is None:
        return None
    dur = _num(ts.get("duration"))
    if dur is None:
        return None
    over = [t for t in (ts.get("save_times") or []) if _num(t) is not None
            and _num(t) > dur]
    if not over:
        return None
    return (f"{len(over)} extra save time(s) lie beyond the Duration ({dur:g}) -- "
            f"{', '.join(f'{_num(t):g}' for t in over[:4])}. The march stops at the "
            f"duration, so no frame is saved at those times {_AT_TSEEP}.")


@rule("tseep.initial_condition", INFO, ("tseep",),
      "The initial condition is a steady solve at the t = 0 boundary configuration.")
def _tseep_initial_condition(ctx):
    ts = ctx.tseep
    if ts is None:
        return None
    times = [t for t in (ts.get("times") or []) if _num(t) is not None]
    t0 = min(times) if times else None
    out = []
    for b in (ctx.seep_bc.get("specified_heads") or []) + \
             (ctx.seep_bc.get("specified_fluxes") or []):
        name = b.get("head", b.get("flux"))
        if not isinstance(name, str):
            continue
        if ctx.series_at_start(name) is None:
            out.append(
                f"Time series {name!r} has no value at the first time "
                f"({'t = %g' % t0 if t0 is not None else 'the start of the run'}), "
                f"so it holds flat at its first defined breakpoint until that point. "
                f"The initial condition is a steady solve at the t = 0 boundary "
                f"configuration, so give the driving series its t = 0 value to "
                f"control the state the march starts from {_AT_TSEEP}.")
    return out


@rule("seep.steady_reads_series_at_t0", INFO, ("seep",),
      "A steady run reads each series-bound boundary at its t = 0 value.")
def _steady_reads_series_at_t0(ctx):
    # Fires on the STEADY run only. A transient run inherits every seep rule
    # (analyses expansion), but there the series drive the march itself -- guard
    # on the exact analysis, not the expanded set.
    if ctx.analysis_name != "seep":
        return None
    bc = ctx.seep_bc
    names = [b.get("head") for b in bc.get("specified_heads") or []]
    names += [b.get("flux") for b in bc.get("specified_fluxes") or []]
    bound = list(dict.fromkeys(v for v in names if isinstance(v, str) and v.strip()))
    if not bound:
        return None
    parts = []
    for s in bound:
        v = ctx.series_at_start(s)
        parts.append(f"'{s}'" if v is None else f"'{s}' = {v:g}")
    listed = ", ".join(parts)
    return (f"Boundary values bound to a time series are read at their t = 0 "
            f"values for a steady run ({listed}) -- the initial-condition "
            f"snapshot, the same field a transient march starts from. The Log "
            f"states each reading {_AT_SEEPBC}.")


@rule("tseep.reservoir_face_above_level", INFO, ("seep",),
      "A reservoir boundary is submerged-only: nodes above the level become exit faces.")
def _reservoir_above_level(ctx):
    out = []
    for n, b in enumerate(ctx.seep_bc.get("specified_heads") or []):
        if str(b.get("kind") or "").strip().lower() != "reservoir":
            continue
        pts = _coords(b.get("coords"))
        if not pts:
            continue
        top = max(p[1] for p in pts)
        head = b.get("head")
        level = _num(head)
        if level is None and isinstance(head, str):
            level = ctx.series_at_start(head)
            if level is None:
                continue
        if level is None or top <= level + 1e-9:
            continue
        out.append(f"Specified head #{n + 1} is a reservoir face rising to "
                   f"y = {top:g}, above its level ({level:g}). A reservoir boundary "
                   f"is submerged-only: nodes at or below the level are held at it, "
                   f"and the nodes above become seepage-exit faces. That is how a "
                   f"rising pool is modeled -- it is only wrong if the face was "
                   f"meant to be held at the level throughout {_AT_SEEPBC}.")
    return out


# ---------------------------------------------------------------------------
# Family: rapid drawdown
#
# Rapid drawdown inherits every limit-equilibrium rule and adds a second water
# state and a multi-stage strength envelope. Its characteristic failure is not a
# crash: a missing stage-2 water source reads as u2 = 0, which is a pore pressure
# that vanished at the moment of drawdown, and the factor of safety comes back
# HIGH. Every ERROR here is on that unconservative side.
# ---------------------------------------------------------------------------

@rule("rapid.stage2_water_missing", ERROR, ("rapid",),
      "Each pore-pressure option needs its own stage-2 source.")
def _rapid_stage2_water(ctx):
    out = []
    piezo_mats = [ctx.mat_label(i) for i, m in ctx.strength_materials()
                  if str(m.get("u") or "").strip().lower() == "piezo"]
    if piezo_mats and len(ctx.piezo2) < 2:
        out.append(
            f"{piezo_mats[0]} takes pore pressure from a piezometric line, so this "
            f"rapid-drawdown run needs Piezometric Line 2 for the drawn-down pool. "
            f"It carries fewer than two points, and the stage-2 pore pressure would "
            f"read zero everywhere -- a drawdown that removed the water AND its "
            f"pressure, and a factor of safety that is too high (Piezometric lines, "
            f"Line 2; piezo sheet).")
    seep_mats = [ctx.mat_label(i) for i, m in ctx.strength_materials()
                 if str(m.get("u") or "").strip().lower() == "seep"]
    # Two named instants of a transient march ARE the two staged fields: the run
    # writes the stage_1 and stage_2 frames into seep_u / seep_u2 before the solver
    # starts, which is the file pair's in-memory equivalent and the whole point of
    # coupling a drawdown to a march. One instant is not -- that supplies stage 1
    # only -- so the requirement stands unless both are named.
    staged_2 = len(ctx.seep_frame_times) >= 2
    if seep_mats and ctx.sd.get("seep_u2") is None and not staged_2:
        out.append(
            f"{seep_mats[0]} takes pore pressure from a seepage solution, so this "
            f"rapid-drawdown run needs the drawn-down field '{{base}}_seep2.csv' "
            f"beside the .xlsx. It is not loaded, and the stage-2 pore pressure "
            f"would read zero everywhere.")
    return out


@rule("rapid.ru_has_no_stage2", WARNING, ("rapid",),
      "u = ru has no post-drawdown variant: the pore pressure does not change.")
def _rapid_ru(ctx):
    named = [ctx.mat_label(i) for i, m in ctx.strength_materials()
             if str(m.get("u") or "").strip().lower() == "ru"]
    if not named:
        return None
    return (f"{named[0]} uses u = ru, which has no staged variant: the engine "
            f"carries the same pore pressure into stage 2, so drawdown does not "
            f"change the pore pressure in that material at all. Use u = piezo with "
            f"Piezometric Line 2, or u = seep with a drawn-down field {_AT_MAT}.")


@rule("rapid.dloads2_missing", ERROR, ("rapid",),
      "The stage-2 distributed loads carry the post-drawdown water load.")
def _rapid_dloads2(ctx):
    if ctx.sd.get("dloads2"):
        return None
    if ctx.water_loads_mode != "manual":
        return None            # in auto mode the engine derives the stage-2 load
    return ("The stage-2 distributed loads are empty, so the post-drawdown water "
            "load is zero. A rapid-drawdown run needs the stage-2 load: with the "
            "reservoir load simply absent the slope is analysed as though the water "
            "vanished without ever having pressed on it. Enter zero loads "
            "deliberately if the pool drains completely (Distributed loads, Set 2; "
            "dloads (2) sheet).")


@rule("rapid.d_psi_incomplete", ERROR, ("rapid",),
      "The Kc = 1 envelope needs BOTH d and psi.")
def _rapid_d_psi(ctx):
    for i, m in ctx.strength_materials():
        d, psi = _num(m.get("d")), _num(m.get("psi"))
        d = 0.0 if d is None else d
        psi = 0.0 if psi is None else psi
        if (d != 0.0) == (psi != 0.0):
            continue
        have, missing = ("d", "psi") if d != 0.0 else ("psi", "d")
        yield (f"{ctx.mat_label(i)} carries {have} but not {missing}. The "
               f"Duncan & Wright Kc = 1 envelope needs both, and a material with "
               f"only one of them is silently treated as free-draining: it keeps "
               f"its drained c/f through stage 2, which on the reference model was "
               f"worth +12.8% on the factor of safety. Fill both, or clear both to "
               f"declare the material free-draining {_AT_MAT}.")


@rule("rapid.no_low_k_material", ERROR, ("rapid",),
      "Rapid drawdown needs at least one material with a Kc = 1 envelope.")
def _rapid_no_low_k(ctx):
    mats = list(ctx.strength_materials())
    if not mats:
        return None
    if any((_num(m.get("d")) or 0.0) != 0.0 or (_num(m.get("psi")) or 0.0) != 0.0
           for _, m in mats):
        return None
    return ("No material carries the d / psi pair, so there is no "
            "low-permeability material for the three-stage procedure to work on and "
            "the run is refused. Enter d and psi for the materials that do not "
            "drain during drawdown (Materials table; mat sheet).")


@rule("rapid.free_draining_materials", WARNING, ("rapid",),
      "A material with no d/psi keeps its drained strength through stage 2.")
def _rapid_free_draining(ctx):
    mats = list(ctx.strength_materials())
    if not mats:
        return None
    free = [ctx.mat_label(i) for i, m in mats
            if (_num(m.get("d")) or 0.0) == 0.0 and (_num(m.get("psi")) or 0.0) == 0.0]
    if not free or len(free) == len(mats):
        return None            # none, or all -- the all case is its own ERROR
    return (f"{len(free)} of {len(mats)} material(s) a failure surface can cross "
            f"carry no d / psi and are treated as free-draining through the "
            f"drawdown, keeping their drained strength: {free[0]}"
            + (f" and {len(free) - 1} other(s)" if len(free) > 1 else "")
            + f" {_AT_MAT}.")


@rule("rapid.cp_ignores_d_psi", WARNING, ("rapid",),
      "An undrained (option = cp) material cannot also carry a Kc = 1 envelope.")
def _rapid_cp(ctx):
    for i, m in ctx.strength_materials():
        if str(m.get("option") or "").strip().lower() != "cp":
            continue
        if (_num(m.get("d")) or 0.0) == 0.0 and (_num(m.get("psi")) or 0.0) == 0.0:
            continue
        yield (f"{ctx.mat_label(i)} uses option = cp, so its d / psi are set to "
               f"zero for the drawdown: an undrained material can never be "
               f"low-permeability in the Duncan & Wright sense whatever the "
               f"materials table says. The stage-2 strength comes from c/p alone "
               f"{_AT_MAT}.")


@rule("rapid.curved_envelope_unsupported", ERROR, ("rapid",),
      "option = pow and option = hb are not supported in rapid drawdown.")
def _rapid_curved(ctx):
    for i, m in ctx.strength_materials():
        opt = str(m.get("option") or "").strip().lower()
        if opt not in ("pow", "hb"):
            continue
        name = "power-curve" if opt == "pow" else "Hoek-Brown"
        yield (f"{ctx.mat_label(i)} uses option = {opt} ({name}), which rapid "
               f"drawdown does not support: the staged procedure overwrites each "
               f"slice's c and f, and the curved envelope computes them itself. "
               f"The two would clobber each other, so the run is refused "
               f"{_AT_MAT}.")


@rule("rapid.pool_rises", ERROR, ("rapid",),
      "The drawn-down pool cannot stand higher than the full pool.")
def _rapid_pool_rises(ctx):
    p1, p2 = ctx.piezo, ctx.piezo2
    if len(p1) < 2 or len(p2) < 2:
        return None
    if _monotonic_x(p1) == "mixed" or _monotonic_x(p2) == "mixed":
        return None            # the ordering rules report that; measuring would lie
    worst, where = 0.0, None
    for x in sorted({p[0] for p in p1} | {p[0] for p in p2}):
        y1, y2 = _y_at(p1, x), _y_at(p2, x)
        if y1 is None or y2 is None:
            continue
        if y2 - y1 > worst:
            worst, where = y2 - y1, x
    height = ctx.slope_height or 1.0
    if where is None or worst <= max(1e-9, 1e-4 * height):
        return None
    return (f"Piezometric Line 2 stands {worst:.4g} ABOVE Line 1 near "
            f"x = {where:.6g}. Line 2 is the drawn-down pool and Line 1 the full "
            f"pool, so the post-drawdown water cannot be higher. Check the two lines "
            f"have not been entered the wrong way round {_AT_PIEZO}.")


@rule("rapid.stage2_water_without_load", WARNING, ("rapid",),
      "Stage-2 water standing above the ground needs a stage-2 distributed load.",
      remedy="add_ponded_water_load")
def _rapid_stage2_ponded(ctx):
    if ctx.water_loads_mode != "manual":
        return None            # auto mode derives the stage-2 load from the pool
    pz, gs = ctx.piezo2, ctx.ground
    if len(pz) < 2 or len(gs) < 2:
        return None
    if _monotonic_x(gs) == "mixed" or _monotonic_x(pz) == "mixed":
        return None
    wet, depth = [], 0.0
    for x in sorted({p[0] for p in pz} | {p[0] for p in gs}):
        yp, yg = _y_at(pz, x), _y_at(gs, x)
        if yp is None or yg is None:
            continue
        if yp - yg > 1e-9:
            wet.append(x)
            depth = max(depth, yp - yg)
    if not wet:
        return None
    a, b = min(wet), max(wet)
    height = ctx.slope_height or 0.0
    if b - a <= 0.0 or depth < max(1e-9, 1e-4 * height):
        return None
    for _, pts in ctx.dload_blocks(2):
        ext = _extent(pts)
        if ext and ext[0] <= b and ext[1] >= a:
            return None
    return (f"Piezometric Line 2 rises up to {depth:.3g} above the ground surface "
            f"between x = {a:.6g} and x = {b:.6g}, but no stage-2 distributed load "
            f"covers that reach. The drawn-down pool still presses on the slope, and "
            f"stage 2 is missing its weight {_AT_DLOADS2}.")


@rule("rapid.stage2_load_repeats_stage1", WARNING, ("rapid",),
      "A stage-2 load equal to stage 1 while the pool dropped is the copied-set "
      "error.")
def _rapid_stage2_repeat(ctx):
    if ctx.water_loads_mode != "manual":
        return None
    p1, p2 = ctx.piezo, ctx.piezo2
    if len(p1) < 2 or len(p2) < 2:
        return None
    if not ctx.sd.get("dloads") or not ctx.sd.get("dloads2"):
        return None
    # How far the pool actually fell, measured where both lines are defined.
    drops = []
    for x in sorted({p[0] for p in p1} | {p[0] for p in p2}):
        y1, y2 = _y_at(p1, x), _y_at(p2, x)
        if y1 is not None and y2 is not None:
            drops.append(y1 - y2)
    if not drops:
        return None
    fall = max(drops)
    height = ctx.slope_height or 1.0
    if fall <= max(1e-9, 0.01 * height):
        return None            # the pool did not meaningfully move
    r1, r2 = ctx.dload_resultant(1), ctx.dload_resultant(2)
    if r1 <= 0:
        return None
    if abs(r2 - r1) > 0.02 * r1:
        return None
    return (f"The stage-2 distributed loads carry the same resultant as stage 1 "
            f"({r2:.6g} against {r1:.6g}) while Piezometric Line 2 sits up to "
            f"{fall:.4g} below Line 1. The reservoir that drained is still pressing "
            f"on the slope at full height -- the usual cause is a stage-2 set copied "
            f"from stage 1 and never re-drawn {_AT_DLOADS2}.")


# ---------------------------------------------------------------------------
# Family: seismic direction -- one cell, two meanings
#
# The limit-equilibrium engine takes the MAGNITUDE of k and orients it itself
# (slice.py:1689, kw = abs(k) * W); the finite element engine reads the SIGN as
# direction (fem.py:6911, a signed horizontal body force). Both are documented,
# each on its own engine's page, and neither page mentions the other. These two
# INFOs state the convention of the engine actually being run, which is also the
# copy a Run dialog shows.
# ---------------------------------------------------------------------------

@rule("main.seismic_direction_lem", INFO, ("lem",),
      "The LEM applies k in the failure-driving direction and ignores its sign.",
      fields=("k_seismic",))
def _seismic_direction_lem(ctx):
    k = _num(ctx.sd.get("k_seismic"))
    if k is None or k == 0.0:
        return None
    return (f"Seismic coefficient k = {k:g} is applied at magnitude {abs(k):g} in "
            f"the failure-driving direction for each surface: the direction is "
            f"automatic in the limit-equilibrium engine, so the sign is ignored. "
            f"Enter a magnitude -- the finite element engine reads the same value as "
            f"a vector, where +k pushes +x and -k pushes -x {_at_global('D13')}.")


@rule("main.seismic_direction_fem", INFO, ("fem",),
      "In the FEM the sign of k sets the direction of the body force.",
      fields=("k_seismic",))
def _seismic_direction_fem(ctx):
    k = _num(ctx.sd.get("k_seismic"))
    if k is None or k == 0.0:
        return None
    pushes = "+x" if k > 0 else "-x"
    return (f"Seismic coefficient k = {k:g} pushes in {pushes}. In the finite "
            f"element engine the sign IS the direction: both faces are analysed at "
            f"once, so choose the direction that drives the face you are checking -- "
            f"the engine will not choose for you, and a pseudo-static factor of "
            f"safety can legitimately come out above the static one for the face the "
            f"shaking stabilises. The limit-equilibrium engine reads the same value "
            f"as a magnitude and orients it itself {_at_global('D13')}.")


# ---------------------------------------------------------------------------
# Family: tension cracks
#
# The crack is a limit-equilibrium construction: fem.py contains no reference to
# tcrack of any kind, and the crack surface is written to its own slope_data key
# without touching ground_surface or the polygons, so a finite element result is
# provably independent of it. One file therefore poses two different problems the
# moment a crack depth is entered, and the cross-analysis rule says so.
# ---------------------------------------------------------------------------

def _theoretical_crack_depth(ctx):
    """``(depth, label)`` for the largest 2c/gamma among the surface materials."""
    best, label = None, None
    for i, m in ctx.strength_materials():
        c, g = _num(m.get("c")), _num(m.get("gamma"))
        if c is None or g is None or g <= 0:
            continue
        d = 2.0 * c / g
        if best is None or d > best:
            best, label = d, ctx.mat_label(i)
    return best, label


@rule("crack.deeper_than_slope", ERROR, ("lem",),
      "A crack at or below the base of the slope cannot form.",
      fields=("tcrack_depth",))
def _crack_deeper_than_slope(ctx):
    d = _num(ctx.sd.get("tcrack_depth"))
    h = ctx.slope_height
    if d is None or d <= 0 or h is None or h <= 0:
        return None
    if d < h:
        return None
    return (f"Tension crack depth = {d:g} reaches at or below the base of a slope "
            f"{h:g} high. No crack can form, and the surface generator silently "
            f"falls back to the uncracked geometry -- while the thrust from Water in "
            f"crack is still applied, because that force is computed from the water "
            f"depth alone and never checks the crack {_at_global('D11')}.")


@rule("crack.exceeds_theoretical_depth", INFO, ("lem",),
      "A crack far deeper than 2c/gamma is a geometric feature, not a Rankine estimate.",
      fields=("tcrack_depth", "c", "phi", "gamma"))
def _crack_theoretical(ctx):
    d = _num(ctx.sd.get("tcrack_depth"))
    if d is None or d <= 0:
        return None
    theo, label = _theoretical_crack_depth(ctx)
    if theo is None or theo <= 0:
        return None
    if d <= 3.0 * theo:
        return None
    return (f"Tension crack depth = {d:g}, against a theoretical depth "
            f"2c/gamma = {theo:.4g} for {label} -- {d / theo:.1f} times it. That is "
            f"legitimate when the crack is a stated feature of the problem rather "
            f"than a Rankine estimate; it is worth checking that the depth is the "
            f"one intended {_at_global('D11')}.")


@rule("crack.cohesionless_materials", INFO, ("lem",),
      "The theoretical crack depth is zero in a cohesionless material.",
      fields=("tcrack_depth", "c"))
def _crack_cohesionless(ctx):
    d = _num(ctx.sd.get("tcrack_depth"))
    if d is None or d <= 0:
        return None
    mats = list(ctx.strength_materials())
    if not mats:
        return None
    if any((_num(m.get("c")) or 0.0) > 0 for _, m in mats):
        return None
    return (f"Tension crack depth = {d:g}, but every material a failure surface "
            f"can cross is cohesionless (c = 0), where the theoretical crack depth "
            f"2c/gamma is zero. The crack is a procedural device here rather than a "
            f"material property {_at_global('D11')}.")


@rule("crack.no_surface_intersection", WARNING, ("lem",),
      "A crack that misses every surface applies no crack -- but still applies its water.",
      fields=("tcrack_depth",))
def _crack_no_intersection(ctx):
    d = _num(ctx.sd.get("tcrack_depth"))
    if d is None or d <= 0 or ctx.is_search:
        return None
    gs = ctx.ground
    if len(gs) < 2 or not ctx.surfaces:
        return None
    from shapely.geometry import LineString
    crack = LineString([(x, y - d) for x, y in gs])
    if any(crack.intersects(s) for s in ctx.surfaces):
        return None
    w = _num(ctx.sd.get("tcrack_water")) or 0.0
    where = _at_global("D11")
    if w <= 0:
        return (f"Tension crack depth = {d:g} produces a crack surface that "
                f"intersects none of the failure surfaces this file defines, so no "
                f"crack is applied to the analysis {where}.")
    return (f"Tension crack depth = {d:g} produces a crack surface that intersects "
            f"none of the failure surfaces this file defines, so no crack is applied "
            f"to the analysis. The thrust from Water in crack ({w:g} deep) is "
            f"applied anyway: it is computed from the water depth alone and never "
            f"checks that a crack was formed {where}.")


@rule("crack.ignored_by_fem", WARNING, ("fem",),
      "A tension crack is a limit-equilibrium construction; the FEM does not model it.",
      fields=("tcrack_depth",))
def _crack_ignored_by_fem(ctx):
    d = _num(ctx.sd.get("tcrack_depth"))
    if d is None or d <= 0:
        return None
    return (f"Tension crack depth is {d:g}, but this run does not include the "
            f"crack: the finite element engine represents a crack CONSTITUTIVELY, "
            f"through tensile strength (t_cut on the Materials table), never "
            f"geometrically, so its result is independent of the depth. Comparing "
            f"this run against the limit-equilibrium answer for the same file "
            f"compares a cracked surface with an uncracked continuum "
            f"{_at_global('D11')}.")


# ---------------------------------------------------------------------------
# Family: piles
#
# S (pile spacing) converts between the per-pile capacities and the per-unit-width
# forces a plane-strain analysis needs, and it reaches four division sites in two
# files with four different failure modes -- including one that is silent and
# unconservative: a NEGATIVE S passes the capacity cap through unchanged
# (min(F*S, Vcap) / S == F when S < 0), so the structural capacity is disabled
# while the pile still contributes its full force.
# ---------------------------------------------------------------------------

def _pile_uses_spacing(ctx, p):
    """Why this pile row needs S, or ``None`` when nothing on this run divides by it."""
    if "fem" in ctx.analyses:
        return ("the finite element engine divides the pile's EA and EI by it for "
                "the per-unit-width beam stiffness")
    if _num(p.get("H")) is None:
        return ("H is blank, so the pile force is computed by the Ito & Matsui "
                "method, which needs the clear spacing between piles")
    if p.get("V_cap") is not None or p.get("M_cap") is not None:
        return ("the structural capacity cap (Vcap / Mcap) is per pile, and S "
                "converts it to the per-unit-width force the slices carry")
    return None


@rule("pile.spacing_invalid", ERROR, ("lem", "fem"),
      "Pile spacing S must be positive wherever the run divides by it.",
      fields=("S",))
def _pile_spacing_invalid(ctx):
    for i, p in enumerate(ctx.piles):
        s = _num(p.get("S"))
        need = _pile_uses_spacing(ctx, p)
        if s is not None and s > 0:
            continue
        if s is None:
            if need is None:
                yield (INFO, f"{ctx.pile_label(i)} has no spacing: S is blank. "
                             f"Nothing on this run divides by it, so it is inert "
                             f"here -- but a finite element run of the same file "
                             f"would need it {_AT_PILES}.")
                continue
            yield (f"{ctx.pile_label(i)} has no spacing: S is blank, and {need}. "
                   f"Enter the center-to-center spacing {_AT_PILES}.")
            continue
        if s == 0:
            yield (f"{ctx.pile_label(i)} has a spacing of 0. S is a divisor -- "
                   f"{need or 'the engines divide the per-pile quantities by it'} -- "
                   f"and the run stops with a bare division-by-zero that names "
                   f"neither the input nor the row {_AT_PILES}.")
            continue
        yield (f"{ctx.pile_label(i)} has a negative spacing: S = {s:g}. That is not "
               f"merely wrong, it is silently unconservative: the capacity cap is "
               f"applied as min(F*S, Vcap)/S, which returns F unchanged when S is "
               f"negative, so Vcap and Mcap stop limiting the pile force while the "
               f"pile keeps contributing it in full {_AT_PILES}.")


@rule("pile.spacing_not_greater_than_diameter", ERROR, ("lem",),
      "Ito & Matsui needs a clear gap between piles: S must exceed D.",
      fields=("S", "D"))
def _pile_spacing_vs_diameter(ctx):
    for i, p in enumerate(ctx.piles):
        if _num(p.get("H")) is not None:
            continue               # H given: the arching theory is never invoked
        s, d = _num(p.get("S")), _num(p.get("D_pile"))
        if s is None or d is None or s <= 0 or d <= 0 or s > d:
            continue
        yield (f"{ctx.pile_label(i)} has S = {s:g}, which is not greater than "
               f"D = {d:g}. With H blank the pile force comes from the Ito & Matsui "
               f"arching theory, which needs a clear gap between piles -- a "
               f"continuous wall (S <= D) has no soil arching and the theory does "
               f"not apply. Enter H directly for a wall, or increase the spacing "
               f"{_AT_PILES}.")


@rule("pile.sd_ratio_outside_band", WARNING, ("lem",),
      "Ito & Matsui is calibrated for S/D between about 2 and 8.",
      fields=("S", "D"))
def _pile_sd_ratio_band(ctx):
    for i, p in enumerate(ctx.piles):
        if _num(p.get("H")) is not None:
            continue               # H given: the arching theory is never invoked
        s, d = _num(p.get("S")), _num(p.get("D_pile"))
        if s is None or d is None or d <= 0 or s <= d:
            continue               # blank or wall-like: other rules carry those
        ratio = s / d
        if 2.0 <= ratio <= 8.0:
            continue
        if ratio < 2.0:
            yield (f"{ctx.pile_label(i)} has S/D = {ratio:g}, below the S/D of "
                   f"about 2 to 8 the Ito & Matsui theory is applicable for. "
                   f"This close, the row acts more like a continuous wall than "
                   f"a line of piles soil can flow between, and the computed "
                   f"force can be unrealistically large. Widen the spacing, or "
                   f"enter H directly for a wall {_AT_PILES}.")
        else:
            yield (f"{ctx.pile_label(i)} has S/D = {ratio:g}, above the S/D of "
                   f"about 2 to 8 the Ito & Matsui theory is applicable for. "
                   f"This far apart, the soil arching the theory is built on "
                   f"fades away and the computed force overestimates what the "
                   f"row can develop. Tighten the spacing, or enter a measured "
                   f"H directly {_AT_PILES}.")


@rule("pile.fem_incomplete", ERROR, ("fem",),
      "A pile the FEM cannot build is deleted from the model without a word.")
def _pile_fem_incomplete(ctx):
    for i, p in enumerate(ctx.piles):
        E = _num(p.get("E"))
        if E is None or E <= 0:
            yield (f"{ctx.pile_label(i)} has no modulus: E is blank, so the finite "
                   f"element engine drops every element of this pile from the model "
                   f"-- the pile is still drawn, and the run reports no problem. The "
                   f"FEM models a pile as a beam element and needs E with either D "
                   f"or an Area and I; the LEM does not, because it applies the pile "
                   f"force envelope directly, so this file can run in the LEM but "
                   f"not honestly in the FEM until E is filled in {_AT_PILES}.")
            continue
        if _num(p.get("D_pile")) is not None:
            continue               # Area and I are derived from the diameter
        area, I = _num(p.get("area")), _num(p.get("I"))
        if area is not None and area > 0 and I is not None and I > 0:
            continue
        if (area is None or area <= 0) and (I is None or I <= 0):
            yield (f"{ctx.pile_label(i)} has neither D nor an Area/I pair, so the "
                   f"pile's axial and bending stiffness are both zero and the finite "
                   f"element solve fails with 'Factor is exactly singular'. Enter "
                   f"the pile diameter, or the section properties {_AT_PILES}.")
        elif area is None or area <= 0:
            yield (f"{ctx.pile_label(i)} has I but neither D nor Area, so the axial "
                   f"stiffness EA is zero -- an axially free beam that solves cleanly "
                   f"and carries nothing. The result looks like a run without the "
                   f"pile at all {_AT_PILES}.")
        else:
            yield (f"{ctx.pile_label(i)} has Area but neither D nor I, so the pile "
                   f"has no bending stiffness. Enter the diameter, or the second "
                   f"moment of area {_AT_PILES}.")


@rule("pile.units_convention", INFO, ("lem", "fem"),
      "H is per unit width; D, S, Vcap, Mcap, E, I and Area are per pile.")
def _pile_units_convention(ctx):
    rows = [i for i, p in enumerate(ctx.piles)
            if _num(p.get("H")) is not None
            and (p.get("V_cap") is not None or p.get("M_cap") is not None)]
    if not rows:
        return None
    return (f"{ctx.pile_label(rows[0])} carries both H and a structural capacity. "
            f"They are in different conventions and the table cannot show it: H is a "
            f"force PER UNIT WIDTH of slope, while Vcap and Mcap are PER PILE and "
            f"are divided by S before they are compared with it {_AT_PILES}.")


@rule("pile.section_derived", INFO, ("fem",),
      "A pile with a diameter and no section gets Area and I from the circle.")
def _pile_section_derived(ctx):
    rows = []
    for i, p in enumerate(ctx.piles):
        d = _num(p.get("D_pile"))
        if d is None or d <= 0:
            continue
        if _num(p.get("area")) is None or _num(p.get("I")) is None:
            rows.append(i)
    if not rows:
        return None
    return (f"{len(rows)} pile row(s) carry a diameter with no Area or I, so the "
            f"engine derives the solid circular section: Area = pi*D^2/4 and "
            f"I = pi*D^4/64. Enter Area and I directly for any other section "
            f"{_AT_PILES}.")


# ---------------------------------------------------------------------------
# Family: reinforcement
#
# The reinforcement path already carries the best-worded input error in the
# package (fem.py:1145) -- the one that tells the user this file runs in the LEM
# but not the FEM. These rules bring that sentence forward to preflight, and add
# its mirror image: the same file, checked BEFORE a finite element run is asked
# for, so the gap is visible while the limit-equilibrium answer is being read.
# ---------------------------------------------------------------------------

def _reinf_line_length(r):
    try:
        return math.hypot(float(r["x2"]) - float(r["x1"]),
                          float(r["y2"]) - float(r["y1"]))
    except (KeyError, TypeError, ValueError):
        return None


@rule("reinforce.dir_vocabulary", ERROR, ("*",),
      "Type, Dir and Appl each speak a fixed vocabulary.")
def _reinf_vocabulary(ctx):
    types = ("", "geosynthetic", "nail", "tieback", "anchor")
    for i, r in enumerate(ctx.reinforcement):
        t = str(r.get("type") or "").strip().lower()
        if t not in types:
            yield (f"{ctx.reinf_label(i)} sets Type = {r.get('type')!r}, which is "
                   f"not a support type. Expected one of: geosynthetic, nail, "
                   f"tieback, anchor (or blank for a generic tensile line) "
                   f"{_AT_REINF}.")
        d = str(r.get("dir") or "").strip().lower()
        if d and d not in ("tangent", "axial"):
            yield (f"{ctx.reinf_label(i)} sets Dir = {r.get('dir')!r}, which is "
                   f"not a direction. Expected: tangent (along the failure surface) "
                   f"or axial (along the bar) {_AT_REINF}.")
        a = str(r.get("appl") or "").strip().lower()
        if a and a not in ("active", "passive"):
            yield (f"{ctx.reinf_label(i)} sets Appl = {r.get('appl')!r}, which is "
                   f"not an application. Expected: active (allowable force) or "
                   f"passive (ultimate capacity, divided by the factor of safety) "
                   f"{_AT_REINF}.")


@rule("reinforce.tmax_nonpositive", ERROR, ("*",),
      "Tmax is the capacity the whole pullout envelope is built from.",
      fields=("t_max",))
def _reinf_tmax(ctx):
    for i, r in enumerate(ctx.reinforcement):
        t = _num(r.get("t_max"))
        if t is not None and t > 0:
            continue
        yield (f"{ctx.reinf_label(i)} has no tensile capacity: Tmax is "
               f"{_fmt(r.get('t_max'))}. The pullout envelope is built by scaling "
               f"the anchorage ramps against Tmax, so a zero divides by zero while "
               f"the line is being built -- a bare ZeroDivisionError naming no sheet "
               f"and no row -- and a line with no capacity reinforces nothing in any "
               f"case {_AT_REINF}.")


@rule("reinforce.fem_incomplete", ERROR, ("fem",),
      "The FEM models reinforcement as a bar element and needs E and Area.")
def _reinf_fem_incomplete(ctx):
    for i, r in enumerate(ctx.reinforcement):
        E, area = _num(r.get("E")), _num(r.get("area"))
        if E is not None and E > 0 and area is not None and area > 0:
            continue
        yield (f"{ctx.reinf_label(i)} has no usable axial stiffness "
               f"(E = {_fmt(r.get('E'))}, Area = {_fmt(r.get('area'))}). The finite "
               f"element engine models reinforcement as a bar element, so it needs "
               f"both. The limit-equilibrium engine does not -- it applies the "
               f"tensile capacity envelope (Tmax/Lp) directly -- so this file runs "
               f"in the LEM but not in the FEM until E and Area are filled in "
               f"{_AT_REINF}.")


# INFO, not WARNING: on an LEM run this model is complete — a missing FEM-only
# input is not a defect of the analysis being run (owner ruling 2026-08-14).
# The FEM's own preflight refuses the file with the full explanation.
@rule("reinforce.fem_incomplete_on_lem", INFO, ("lem",),
      "Reinforcement complete for the LEM can still be incomplete for the FEM.")
def _reinf_incomplete_cross(ctx):
    bad = []
    for i, r in enumerate(ctx.reinforcement):
        E, area = _num(r.get("E")), _num(r.get("area"))
        if E is None or E <= 0 or area is None or area <= 0:
            bad.append(i)
    if not bad:
        return None
    return (f"{ctx.reinf_label(bad[0])}"
            + (f" and {len(bad) - 1} other line(s)" if len(bad) > 1 else "")
            + f" carr{'y' if len(bad) > 1 else 'ies'} a tensile capacity but no "
              f"axial stiffness (E and Area). That is complete for this "
              f"limit-equilibrium run, which applies the Tmax/Lp envelope directly, "
              f"and incomplete for a finite element run of the same file, which "
              f"models each line as a bar element and would refuse. The two "
              f"analyses are not approximating the same reinforcement "
              f"{_AT_REINF}.")


@rule("reinforce.type_blank", WARNING, ("lem", "fem"),
      "A blank Type defaults to tangent/active -- 13% off a soil nail.")
def _reinf_type_blank(ctx):
    rows = [i for i, r in enumerate(ctx.reinforcement)
            if not str(r.get("type") or "").strip()]
    if not rows:
        return None
    return (f"{len(rows)} reinforcement line(s) leave Type blank, starting at "
            f"{ctx.reinf_label(rows[0])}, so they default to Dir = tangent, "
            f"Appl = active. That is right for a generic tensile line and wrong for "
            f"a soil nail, which is axial/passive: on the reference model the "
            f"default read 13% higher than the nail preset. Set Type, or set Dir "
            f"and Appl explicitly {_AT_REINF}.")


@rule("reinforce.pullout_negative", ERROR, ("lem", "fem"),
      "A negative pullout length is read as fully anchored -- the opposite of the intent.",
      fields=("lp1", "lp2", "tend1", "tend2"))
def _reinf_pullout_negative(ctx):
    for i, r in enumerate(ctx.reinforcement):
        for key, col in (("lp1", "Lp1"), ("lp2", "Lp2")):
            v = _num(r.get(key))
            if v is not None and v < 0:
                yield (f"{ctx.reinf_label(i)} has {col} = {v:g}. A non-positive "
                       f"pullout length is read as FULLY ANCHORED, so a sign typo "
                       f"turns the weakest pullout case into the strongest and "
                       f"raises the factor of safety. Enter the anchorage length, "
                       f"or 0 for a fully anchored end {_AT_REINF}.")


@rule("reinforce.pullout_law_incomplete", ERROR, ("lem", "fem"),
      "Adhesion and Delta describe one pullout law: both, or neither.",
      fields=("adhesion", "delta"))
def _reinf_pullout_law_incomplete(ctx):
    for i, r in enumerate(ctx.reinforcement):
        a, d = _num(r.get("adhesion")), _num(r.get("delta"))
        if (a is None) == (d is None):
            continue
        have, missing = ("Adhesion", "Delta") if a is not None else ("Delta", "Adhesion")
        yield (f"{ctx.reinf_label(i)} fills {have} but leaves {missing} blank. The "
               f"two are one law -- pullout resistance 2*(Adhesion + sigma'v*tan "
               f"Delta) per unit length -- and half of it is not a law: the missing "
               f"half would have to be read as zero, which is either no interface "
               f"friction at all or no adhesion at all, neither of them what a "
               f"half-filled row means. Fill both to use the overburden-dependent "
               f"law, or clear both to use the development lengths Lp1/Lp2 "
               f"{_AT_REINF}.")


@rule("reinforce.pullout_delta_range", ERROR, ("lem", "fem"),
      "The interface friction angle is an angle between 0 and 90 degrees.",
      fields=("delta",))
def _reinf_pullout_delta_range(ctx):
    for i, r in enumerate(ctx.reinforcement):
        d = _num(r.get("delta"))
        if d is None or 0.0 < d < 90.0:
            continue
        yield (f"{ctx.reinf_label(i)} has Delta = {d:g} degrees. The interface "
               f"friction angle must lie strictly between 0 and 90: at 0 the "
               f"overburden contributes nothing and the law reduces to adhesion "
               f"alone (leave Delta blank and use Lp instead if that is the "
               f"intent), and at 90 tan Delta is infinite. A value entered as a "
               f"tangent or a percentage rather than degrees lands here "
               f"{_AT_REINF}.")


@rule("reinforce.pullout_lp_ignored", INFO, ("lem", "fem"),
      "Lp1/Lp2 play no part on a line that uses the overburden-dependent law.",
      fields=("adhesion", "delta", "lp1", "lp2"))
def _reinf_pullout_lp_ignored(ctx):
    rows = []
    for i, r in enumerate(ctx.reinforcement):
        a, d = _num(r.get("adhesion")), _num(r.get("delta"))
        if a is None or d is None:
            continue
        if _num(r.get("lp1")) or _num(r.get("lp2")):
            rows.append(i)
    if not rows:
        return []
    return [f"{len(rows)} reinforcement line(s) carry both a development length "
            f"and the overburden-dependent pullout law, starting at "
            f"{ctx.reinf_label(rows[0])}. Adhesion and Delta govern there: the "
            f"capacity is the interface resistance integrated from each end, and "
            f"Lp1/Lp2 are not read. Clear Adhesion and Delta to go back to the "
            f"development lengths {_AT_REINF}."]


@rule("reinforce.envelope_inconsistent", WARNING, ("lem", "fem"),
      "Tres above Tmax, or a pullout length longer than the line itself.",
      fields=("t_max", "t_res", "lp1", "lp2", "tend1", "tend2"))
def _reinf_envelope(ctx):
    out = []
    for i, r in enumerate(ctx.reinforcement):
        tmax, tres = _num(r.get("t_max")), _num(r.get("t_res"))
        if tmax is not None and tres is not None and tres > tmax:
            out.append(f"{ctx.reinf_label(i)} has Tres = {tres:g}, above Tmax = "
                       f"{tmax:g}, so the residual capacity is above the peak the "
                       f"bar is supposed to drop from {_AT_REINF}.")
        L = _reinf_line_length(r)
        if L is None or L <= 0:
            continue
        for key, col in (("lp1", "Lp1"), ("lp2", "Lp2")):
            v = _num(r.get(key))
            if v is not None and v > L:
                out.append(
                    f"{ctx.reinf_label(i)} has {col} = {v:g}, longer than the line "
                    f"itself ({L:.4g}), so the pullout envelope never reaches Tmax "
                    f"anywhere along it. A misplaced decimal here silently "
                    f"annihilates the reinforcement while the line stays on the "
                    f"drawing {_AT_REINF}.")
    return out


@rule("reinforce.no_surface_engagement", WARNING, ("lem",),
      "Reinforcement that crosses no failure surface contributes nothing.")
def _reinf_no_engagement(ctx):
    if ctx.is_search:
        return None            # the deck's circles are seeds, not the trial set
    lines = ctx.reinforcement
    if not lines or not ctx.surfaces:
        return None
    missed = [i for i, r in enumerate(lines)
              if ctx.segment_engages_surface(r.get("x1"), r.get("y1"),
                                             r.get("x2"), r.get("y2")) is False]
    if len(missed) < len(lines):
        return None            # at least one line is engaged: the model is doing work
    return (f"None of the {len(lines)} reinforcement line(s) crosses any failure "
            f"surface this file defines, so the reinforcement contributes nothing to "
            f"this run -- the answer is the same as with no reinforcement at all. "
            f"Check the line positions against the surface {_AT_REINF}.")


@rule("pile.no_surface_engagement", WARNING, ("lem",),
      "A pile that crosses no failure surface contributes nothing.")
def _pile_no_engagement(ctx):
    if ctx.is_search:
        return None
    piles = ctx.piles
    if not piles or not ctx.surfaces:
        return None
    missed = [i for i, p in enumerate(piles)
              if ctx.segment_engages_surface(p.get("x1"), p.get("y1"),
                                             p.get("x2"), p.get("y2")) is False]
    if len(missed) < len(piles):
        return None
    return (f"None of the {len(piles)} pile line(s) crosses any failure surface "
            f"this file defines, so no pile force enters this run. A pile that stops "
            f"above the surface, or sits outside its x-range, is silently inert "
            f"{_AT_PILES}.")


# ---------------------------------------------------------------------------
# Family: magnitude plausibility -- the sniff tests
#
# These ask a different question from the unit-system rules above. Those ask "is
# this value in the system the file declares?"; these ask "given that system, is
# this value in the engineering band for its own field type?" -- unit-CONDITIONAL
# rather than unit-DISCRIMINATING, which is the distinction that makes them
# answerable at all.
#
# The bands are calibrated against the verification corpus and are deliberately
# loose. A material's modulus legitimately ranges over two orders of magnitude
# around its soil type's expectation (the corpus spans 0.12x to 12.5x), so a
# tighter band would refuse locked, reproduced models. And a value that MATCHES a
# classifier default is never evidence of anything: (100 000, 0.3) is both this
# package's fallback pair and a perfectly ordinary thing for a vendor to specify.
# These rules flag implausible magnitudes. They never claim a value is unset.
# ---------------------------------------------------------------------------

#: How far one material's Young's modulus may sit from its soil type's expectation
#: before it is reported. The corpus's own legitimate spread is [0.12x, 12.5x]
#: (excluding Hoek-Brown rows, where the classifier has no strength to read), so
#: this band is ~4x looser than the widest correct model in it.
_E_BAND_FACTOR = 50.0

#: The same test applied to the whole material table at once. Every material off
#: in the SAME direction is the wrong-stress-unit signature, and it survives a much
#: tighter threshold: the corpus's widest whole-model divergence is 6.25x.
_E_BAND_MODEL_FACTOR = 15.0


def _e_expectations(ctx):
    """``(index, material, soil, E_expected, E)`` per material with a usable E."""
    from .units import classify_elastic
    out = []
    for i, m in ctx.fem_materials():
        E = _num(m.get("E"))
        if E is None or E <= 0:
            continue
        if str(m.get("option") or "").strip().lower() in ("hb", "elastic"):
            # Hoek-Brown carries no c/phi for the classifier to read, so its
            # expectation is the generic rock midpoint and means nothing; an
            # elastic zone is not a soil at all.
            continue
        soil, E_exp, _nu_exp = classify_elastic(
            m, declared_system=ctx.sd.get("unit_system"))
        if E_exp and E_exp > 0:
            out.append((i, m, soil, E_exp, E))
    return out


@rule("mat.E_off_soil_type_band", WARNING, ("fem",),
      "A modulus far from its own soil type's band, in the declared unit system.")
def _mat_E_band(ctx):
    rows = _e_expectations(ctx)
    if not rows:
        return None
    system = ctx.sd.get("unit_system") or "undeclared"
    out = []
    for i, m, soil, E_exp, E in rows:
        ratio = E / E_exp
        if 1.0 / _E_BAND_FACTOR <= ratio <= _E_BAND_FACTOR:
            continue
        out.append(f"{ctx.mat_label(i)} has E = {E:g}. Its strength reads as "
                   f"{soil.lower()}, whose modulus is around {E_exp:g} in this "
                   f"model's units ({system}) -- the entered value is "
                   f"{ratio:.3g} times that. Check the stress unit E was entered "
                   f"in; xslope never converts {_AT_MAT}.")
    if out:
        return out
    if len(rows) < 2:
        return None
    logs = [math.log(E / E_exp) for _i, _m, _s, E_exp, E in rows]
    g = math.exp(sum(logs) / len(logs))
    if 1.0 / _E_BAND_MODEL_FACTOR <= g <= _E_BAND_MODEL_FACTOR:
        return None
    if not all(x > 0 for x in logs) and not all(x < 0 for x in logs):
        return None            # not a coherent shift: individual values differ
    return (f"Every material's Young's modulus sits {g:.3g} times its soil type's "
            f"expectation, in the same direction. One unusual modulus is ordinary; a "
            f"whole table shifted by one factor is the signature of E entered in the "
            f"wrong stress unit -- the declared system is {system}, and xslope never "
            f"converts {_AT_MAT}.")


@rule("mat.nu_implausible", WARNING, ("fem",),
      "A Poisson's ratio below 0.1 is outside the range of any soil or rock.")
def _mat_nu_band(ctx):
    for i, m in ctx.fem_materials():
        nu = _num(m.get("nu"))
        if nu is None or nu <= 0 or nu >= 0.5:
            continue           # blank / out of range is the ERROR above
        if nu >= 0.1:
            continue
        yield (f"{ctx.mat_label(i)} has nu = {nu:g}. Soils sit at 0.2-0.45 and "
               f"rock at 0.15-0.3; below 0.1 the material has almost no lateral "
               f"coupling, which changes the stress field the strength reduction "
               f"acts on {_AT_MAT}.")


#: Plausible Young's-modulus band for a structural element -- a geosynthetic sheet
#: at the bottom, steel at the top -- expressed in kPa and converted for an
#: Imperial model. The corpus's own structural moduli span 2e4 to 2.5e7 kPa, so
#: this band clears the softest by 20x and the stiffest by 400x: it catches a
#: capacity typed into a stiffness cell or a modulus off by orders of magnitude,
#: and nothing else.
_STRUCT_E_BAND_KPA = (1.0e3, 1.0e10)


def _struct_e_band(ctx):
    from .units import KPA_TO_PSF, normalize_unit_system
    lo, hi = _STRUCT_E_BAND_KPA
    if normalize_unit_system(ctx.sd.get("unit_system")) == "imperial":
        return lo * KPA_TO_PSF, hi * KPA_TO_PSF, "psf"
    return lo, hi, "kPa"


@rule("structural.modulus_off_band", WARNING, ("lem", "fem"),
      "A reinforcement or pile modulus outside the band of any structural material.")
def _structural_modulus_band(ctx):
    lo, hi, unit = _struct_e_band(ctx)
    out = []
    for i, r in enumerate(ctx.reinforcement):
        E = _num(r.get("E"))
        if E is None or E <= 0 or lo <= E <= hi:
            continue
        out.append(f"{ctx.reinf_label(i)} has E = {E:g}. A structural modulus in "
                   f"this model's units runs from about {lo:g} (a geosynthetic) to "
                   f"{hi:g} ({unit}, steel); this value is outside that range "
                   f"altogether. Check it is a modulus and not a capacity "
                   f"{_AT_REINF}.")
    for i, p in enumerate(ctx.piles):
        E = _num(p.get("E"))
        if E is None or E <= 0 or lo <= E <= hi:
            continue
        out.append(f"{ctx.pile_label(i)} has E = {E:g}, outside the plausible "
                   f"structural range of about {lo:g} to {hi:g} {unit} for this "
                   f"model's declared units {_AT_PILES}.")
    return out


# ===========================================================================
# Family: material field ranges
#
# These are the rules a SUBSTITUTED value trips. A sweep or a Taylor perturbation
# writes one field and then runs the engine; without a range check the engine
# accepts the value and answers anyway. Every message below carries the probe that
# measured the wrong answer, because a range rule with no measurement behind it is
# an opinion about what a model should look like.
# ===========================================================================

@rule("mat.strength_negative", ERROR, ("lem", "fem"),
      "A negative cohesion or friction angle is not a strength.",
      fields=("c", "phi", "cp"))
def _mat_strength_negative(ctx):
    # c/p is deliberately NOT checked. It is a GRADIENT -- the rate at which the
    # undrained strength changes with depth below r_elev -- and a negative one is a
    # real soil profile, not a mistake: vp030's Middle Clay carries c/p = -4.71 for
    # the construction-consolidated top metre of Borges & Cardoso's case, where su
    # falls from 8.49 to 4.73 with depth. Only the strengths themselves are signed
    # quantities here.
    for i, m in ctx.strength_materials():
        opt = str(m.get("option") or "").strip().lower()
        cols = {"mc": (("c", "c"), ("phi", "f")),
                "cp": (("c", "c"),)}.get(opt, ())
        for key, col in cols:
            v = _num(m.get(key))
            if v is None or v >= 0:
                continue
            yield (f"{ctx.mat_label(i)} has {col} = {v:g}. A strength cannot be "
                   f"negative. The solvers do not refuse it -- they return a "
                   f"negative factor of safety through the success path (a sweep "
                   f"stepping c to -500 reports FS = -1.279, and phi to -30 reports "
                   f"FS = -0.493), so nothing downstream can tell the answer is "
                   f"meaningless {_AT_MAT}.")


@rule("mat.phi_out_of_range", ERROR, ("lem", "fem"),
      "A friction angle of 90 degrees or more has no finite tangent.",
      fields=("phi",))
def _mat_phi_range(ctx):
    for i, m in ctx.strength_materials():
        if str(m.get("option") or "").strip().lower() != "mc":
            continue
        v = _num(m.get("phi"))
        if v is None or v < 90.0:
            continue
        yield (f"{ctx.mat_label(i)} has f (friction angle) = {v:g} degrees. The "
               f"Mohr-Coulomb strength is c + σ′·tan(φ), which has no finite value "
               f"at 90 degrees and changes sign above it. Enter an angle below 90 "
               f"{_AT_MAT}.")


# ===========================================================================
# Family: parametric sweeps (sensitivity, design, back-analysis, tornado)
#
# A sweep's inputs live on the RUN, not in the template, so every rule here reads
# the selection: which parameter is swept, over what values, in which engine mode.
# The three ERRORs below are silent no-ops -- probed bit-identical answers at every
# swept value -- which is the worst failure mode a study can have, because the
# resulting flat line reads as "this parameter does not matter".
# ===========================================================================

@rule("sensitivity.ru_without_consumer", ERROR, ("sensitivity",), capability="analysis",
      summary="Sweeping ru on a material whose u option does not read it is a no-op.")
def _sens_ru_no_consumer(ctx):
    if ctx.swept_field != "ru":
        return None
    i, m = ctx.swept_material
    if m is None:
        return None
    u = str(m.get("u") or "").strip().lower()
    if u == "ru":
        return None
    return (f"The sweep varies ru on {ctx.mat_label(i)}, but that material's u is "
            f"'{u or 'blank'}', so ru is never read and the factor of safety does "
            f"not move: probed at ru = 0.0, 0.3, 0.9 and 5.0 the answer is "
            f"bit-identical. Set u = ru on the material, or sweep something it reads "
            f"{_AT_MAT}.")


@rule("sensitivity.rapid_only_field", ERROR, ("sensitivity",), capability="analysis",
      summary="d and psi only act in a rapid drawdown analysis, which a sweep cannot run.")
def _sens_rapid_only(ctx):
    field = ctx.swept_field
    if field not in ("d", "psi"):
        return None
    i, _m = ctx.swept_material
    where = ctx.mat_label(i) if i is not None else "the materials table"
    label = "d" if field == "d" else "psi"
    return (f"The sweep varies {label} on {where}. That parameter is read only by "
            f"the rapid drawdown analysis, and a parametric study evaluates the "
            f"single-stage model -- probed at d = 0, 10 and 100 (and psi = 0, 10 "
            f"and 40) the factor of safety is identical at every value {_AT_MAT}.")


@rule("sensitivity.range_crosses_zero", WARNING, ("sensitivity",),
      capability="analysis",
      summary="A swept range that crosses zero takes a positive-only field negative.")
def _sens_range_sign(ctx):
    field = ctx.swept_field
    positive_only = {"gamma": "unit weight", "gamma_sat": "saturated unit weight",
                     "k1": "hydraulic conductivity", "k2": "hydraulic conductivity",
                     "t_max": "tensile capacity", "S": "pile spacing",
                     "D": "pile diameter"}
    if field not in positive_only:
        return None
    vals = ctx.swept_values
    if not vals or min(vals) > 0:
        return None
    lo, hi = min(vals), max(vals)
    doomed = sum(1 for v in vals if v <= 0)
    tail = (f"{field} must stay positive: at zero the engine returns its "
            f"no-solution sentinel and records it as a success, and below zero it "
            f"returns a signed answer. Sweep a strictly positive range.")
    if hi <= 0:
        # Nothing in this range is analysable, so there is no sweep to run and no
        # partial result to preserve -- the only honest answer is to refuse it.
        return (ERROR, f"The sweep of {field} ({positive_only[field]}) runs from "
                       f"{lo:g} to {hi:g}, so every value in it is zero or "
                       f"negative. {tail}")
    return (f"The sweep of {field} ({positive_only[field]}) runs from {lo:g} to "
            f"{hi:g}, which reaches zero or below. {doomed} of {len(vals)} points "
            f"will be skipped. {tail}")


@rule("sensitivity.seismic_sign_discarded", WARNING, ("sensitivity",),
      "The LEM ignores the sign of k, so a range spanning zero is folded in half.")
def _sens_seismic_sign(ctx):
    if ctx.swept_field != "k_seismic" or ctx.sweep_mode != "lem":
        return None
    vals = ctx.swept_values
    if not vals or min(vals) >= 0:
        return None
    return (f"The sweep of k_seismic runs from {min(vals):g} to {max(vals):g}. The "
            f"limit-equilibrium engine applies the magnitude in the failure-driving "
            f"direction and discards the sign, so the negative half of this range "
            f"repeats the positive half: -0.5 and +0.5 return the same factor of "
            f"safety. Sweep magnitudes from 0 upward. (The finite element engine "
            f"does read the sign as a direction.)")


# ===========================================================================
# Family: reliability
#
# The reliability engines refuse several of these at run time already -- but they
# refuse AFTER the critical-surface search, which on a real model is minutes of
# compute thrown away, and one of them raises instead of returning the
# ``(success, result)`` pair its callers handle. Every message here is the one the
# engine itself uses (the builders above), so checking early changes when the user
# is told, never what they are told.
# ===========================================================================

@rule("reliability.no_standard_deviations", ERROR, ("reliability",),
      capability="analysis",
      summary="A reliability run needs at least one non-zero standard deviation.")
def _rel_no_sigmas(ctx):
    for _i, m in enumerate(ctx.materials):
        for key in SIGMA_COLUMNS:
            if _pos(m.get(key)):
                return None
    return no_sigmas_message()


@rule("reliability.elastic_carries_sigma", ERROR, ("reliability",),
      capability="analysis",
      summary="An elastic material has no strength for a standard deviation to move.")
def _rel_elastic_sigma(ctx):
    for i, m in enumerate(ctx.materials):
        if str(m.get("option") or "").strip().lower() != "elastic":
            continue
        stated = sorted(k for k in SIGMA_COLUMNS if _pos(m.get(k)))
        if stated:
            yield elastic_sigma_message(
                m.get("name", f"Material_{i + 1}"), stated)


@rule("reliability.option_unsupported", ERROR, ("reliability",),
      capability="analysis",
      summary="Only mc and cp carry a reliability perturbation set.")
def _rel_option_unsupported(ctx):
    for i, m in enumerate(ctx.materials):
        opt = str(m.get("option") or "").strip()
        if opt.lower() in ("", "mc", "cp", "elastic"):
            continue
        yield unsupported_option_message(m.get("name", f"Material_{i + 1}"), opt)


@rule("reliability.sigma_ignored_by_option", WARNING, ("reliability",),
      "A standard deviation the material's strength model does not read.")
def _rel_sigma_ignored(ctx):
    used = {"mc": ("sigma_gamma", "sigma_c", "sigma_phi"),
            "cp": ("sigma_gamma", "sigma_c", "sigma_cp")}
    for i, m in enumerate(ctx.materials):
        opt = str(m.get("option") or "").strip().lower()
        if opt not in used:
            continue
        for key in SIGMA_COLUMNS:
            if key in used[opt] or not _pos(m.get(key)):
                continue
            yield (f"{ctx.mat_label(i)} has {SIGMA_COLUMNS[key]} = "
                   f"{_fmt(m.get(key))}, but its strength model is option = {opt}, "
                   f"which does not read that parameter. The standard deviation is "
                   f"silently ignored and contributes nothing to the reliability "
                   f"index {_AT_MAT}.")


@rule("reliability.dead_sigma_columns", ERROR, ("reliability",),
      "s(d) and s(psi) are read by no reliability engine.")
def _rel_dead_sigmas(ctx):
    dead = {"sigma_d": "s(d)", "sigma_psi": "s(psi)"}
    live = any(_pos(m.get(k)) for m in ctx.materials for k in SIGMA_COLUMNS)
    for i, m in enumerate(ctx.materials):
        stated = [lbl for key, lbl in dead.items() if _pos(m.get(key))]
        if not stated:
            continue
        text = (f"{ctx.mat_label(i)} carries {' and '.join(stated)}, which no "
                f"reliability engine reads -- the perturbation set is built from "
                f"s(g), s(c), s(f) and s(c/p) only.")
        if live:
            yield (WARNING, text + " This uncertainty is silently dropped from the "
                                   f"reliability index {_AT_MAT}.")
        else:
            yield (text + " They are the only standard deviations in this model, so "
                          "the run is refused for carrying none at all -- a message "
                          "that blames the wrong thing. Enter the uncertainty where "
                          f"an engine reads it {_AT_MAT}.")


@rule("reliability.sigma_exceeds_mean", ERROR, ("reliability",),
      capability="analysis",
      summary="The Taylor series cannot take a parameter below zero at MLV - sigma.")
def _rel_sigma_exceeds_mean(ctx):
    if ctx.reliability_samples_sigmas:
        # The sampling engines truncate every draw at its physical floor, which IS
        # how the published high-COV treatments handle the bound, so the same model
        # that the Taylor series must decline is admissible there.
        return None
    used = {"mc": (("gamma", "sigma_gamma"), ("c", "sigma_c"), ("phi", "sigma_phi")),
            "cp": (("gamma", "sigma_gamma"), ("c", "sigma_c"), ("cp", "sigma_cp"))}
    bad = []
    for i, m in enumerate(ctx.materials):
        opt = str(m.get("option") or "").strip().lower()
        for param, key in used.get(opt, ()):
            std = _num(m.get(key))
            mlv = _num(m.get(param))
            if std is None or std <= 0 or mlv is None:
                continue
            if (mlv - std) < 0:
                bad.append(f"material {i + 1} {param} (mean={mlv:.3g}, "
                           f"sigma={std:.3g})")
    if bad:
        return sigma_exceeds_mean_message("; ".join(bad))
    return None


@rule("reliability.rapid_has_no_effect", ERROR, ("reliability",),
      "Rapid drawdown is dropped silently by two of the reliability modes.")
def _rel_rapid_noop(ctx):
    if not ctx.selection.get("rapid"):
        return None
    if ctx.reliability_samples_sigmas:
        name = ("Monte Carlo" if ctx.reliability_engine == "mc"
                else "response-surface")
        return (f"This run asks for rapid drawdown, but the {name} engine "
                "applies it only to the single critical-surface search at the "
                "most-likely values: every sampled realization is then solved "
                "WITHOUT drawdown, so the probability of failure describes the "
                "wrong analysis. Use the Taylor series engine for a rapid drawdown "
                "reliability, or clear the rapid flag.")
    if ctx.selection.get("search") is False:
        return ("This run asks for rapid drawdown with search=False, and the "
                "fixed-surface evaluation path does not perform the drawdown at "
                "all -- the run reports 'Rapid drawdown: True' and solves the "
                "single-stage model. Run with search=True, or clear the rapid flag.")
    return None


@rule("reliability.mc_samples_too_few", ERROR, ("reliability",),
      "Monte Carlo needs at least two realizations to have a standard deviation.")
def _rel_mc_samples(ctx):
    if ctx.reliability_engine != "mc":
        return None
    n = _num(ctx.selection.get("n_samples"))
    if n is None or n >= 2:
        return None
    return (f"Monte Carlo reliability needs at least 2 samples; this run asks for "
            f"{n:g}. Below two realizations there is no sample standard deviation, "
            f"and the run is refused with a message that blames the model rather "
            f"than the sample count.")


@rule("reliability.fem_mesh_generated", WARNING, ("reliability",),
      "A finite element reliability with no mesh builds one, and the mesh moves beta.")
def _rel_fem_mesh(ctx):
    if "fem" not in ctx.analyses:
        return None
    if ctx.mesh is not None:
        return None
    ext = ctx.ground_extent
    if not ext:
        return None
    size = (ext[1] - ext[0]) / 100.0
    return (f"This finite element reliability run carries no mesh, so one will be "
            f"generated automatically at a target element size of {size:g} (the "
            f"section width divided by 100). Every factor of safety in the "
            f"reliability index comes off that mesh, so the answer depends on a "
            f"discretization nobody chose. Build the mesh explicitly and pass it in "
            f"{_AT_MESH}.")
