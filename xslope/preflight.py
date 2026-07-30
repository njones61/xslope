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
- **A rule names the sheet and the field in the template's own vocabulary** --
  "mat sheet, material 3 ('Core'), column k1", never ``materials[2]['k1']``.

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
its own. Wave 2 builds the contract (:data:`REMEDIES` and the ``remedy`` field on a
finding); the remedy functions themselves arrive with the wave that needs them. The
governing principle is one line: *a remedy is always offered and never applied
silently* -- a fix the user did not ask for is the silent-default disease wearing a
helpful face.
"""

from __future__ import annotations

import math
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

#: Methods whose formulation takes moments about a circle centre and which are
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

#: Named remedies a rule may offer. Wave 2 defines the *contract* -- a rule
#: declares the remedy it would offer and a finding carries the name -- so that the
#: Studio button and the opt-in script form (``preflight(sd, "lem",
#: remedies=[...])``) have one vocabulary to bind to. The transformations
#: themselves arrive with the waves that need them; a name here with no
#: implementation is a declaration of intent, never a silent no-op, because nothing
#: applies a remedy without being asked.
REMEDIES = {
    "set_seismic_zero": "Set main!D13 (Seismic coefficient) to 0 for a static run.",
    "reverse_polyline": "Reverse a right-to-left load or piezometric line.",
    "add_ponded_water_load": "Add the standing water as a distributed load.",
    "extend_piezo_line": "Extend the piezometric line across the section.",
}

#: Capability groups and their options, for :func:`capabilities`. A group's options
#: are the choices an interface offers; a rule that carries ``capability=<group>``
#: is evaluated once per option with that option pinned into the selection, and a
#: rule that fires makes the option unavailable with its message as the reason.
CAPABILITY_GROUPS = {
    "analysis": ("lem", "rapid", "seep", "tseep", "fem", "ssrm"),
    "lem_method": LEM_METHOD_OPTIONS,
}


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
        The name of a remedy this finding could offer (see :data:`REMEDIES`).
        Never applied automatically.
    """

    rule_id: str
    severity: str
    message: str
    remedy: Optional[str] = None

    def __str__(self):
        return f"{self.severity.upper()}: {self.message}  [{self.rule_id}]"


@dataclass
class PreflightReport:
    """The result of a :func:`preflight` run: findings, plus what was checked."""

    analysis: str
    findings: Sequence[Finding] = field(default_factory=tuple)
    selection: dict = field(default_factory=dict)

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
        Name of the remedy findings from this rule may offer.
    capability : str or None
        The :data:`CAPABILITY_GROUPS` group this rule gates, if any.
    availability : bool
        True for a rule that asks *can this analysis be started at all?* rather
        than *is this model right?* -- a missing mesh is the example. The workflow
        satisfies an availability rule before the run gate is ever reached (the
        mesh is an argument to the seepage entry point), so :func:`preflight`
        leaves these to :func:`capabilities`, which is what an interface dims from.
        Pass ``include_availability=True`` to evaluate them anyway.
    """

    id: str
    severity: str
    analyses: Tuple[str, ...]
    check: Callable
    summary: str
    remedy: Optional[str] = None
    capability: Optional[str] = None
    availability: bool = False

    def applies_to(self, analyses):
        return "*" in self.analyses or bool(set(self.analyses) & analyses)


_REGISTRY = []
_REGISTRY_IDS = set()


def rule(id, severity, analyses, summary, remedy=None, capability=None,
         availability=False):
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
    remedy : str, optional
        Key into :data:`REMEDIES`.
    capability : str, optional
        Key into :data:`CAPABILITY_GROUPS`: this rule gates that group's options.
    availability : bool, optional
        Mark the rule as an availability question rather than a validity one; see
        :class:`Rule`.
    """
    if severity not in _SEVERITY_RANK:
        raise ValueError(f"unknown severity {severity!r}")
    for a in analyses:
        if a != "*" and a not in ANALYSES:
            raise ValueError(f"rule {id!r} names unknown analysis {a!r}")
    if capability is not None and capability not in CAPABILITY_GROUPS:
        raise ValueError(f"rule {id!r} names unknown capability group {capability!r}")
    if remedy is not None and remedy not in REMEDIES:
        raise ValueError(f"rule {id!r} names unknown remedy {remedy!r}")
    if id in _REGISTRY_IDS:
        raise ValueError(f"duplicate preflight rule id {id!r}")

    def _register(fn):
        _REGISTRY.append(Rule(id=id, severity=severity, analyses=tuple(analyses),
                              check=fn, summary=summary, remedy=remedy,
                              capability=capability, availability=availability))
        _REGISTRY_IDS.add(id)
        return fn

    return _register


def rules(analysis=None):
    """The registered rules, optionally filtered to one analysis type.

    Returns them in declaration order, which groups related rules together and is
    the order findings come back in.
    """
    if analysis is None:
        return list(_REGISTRY)
    wanted = expand_analysis(analysis)
    return [r for r in _REGISTRY if r.applies_to(wanted)]


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


def _monotonic_x(pts):
    """``"ascending"``, ``"descending"``, or ``"mixed"`` for a polyline's x values."""
    xs = [p[0] for p in pts]
    if len(xs) < 2:
        return "ascending"
    ups = sum(1 for a, b in zip(xs, xs[1:]) if b > a)
    downs = sum(1 for a, b in zip(xs, xs[1:]) if b < a)
    if downs == 0:
        return "ascending"
    if ups == 0:
        return "descending"
    return "mixed"


def _fmt(v):
    """Format a number the way the messages do -- short, and never ``1.0000000001``."""
    f = _num(v)
    return "blank" if f is None else f"{f:g}"


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
        """``"circular"``, ``"noncircular"``, or ``None`` when the run did not say.

        A deck may carry both families; the run is what chooses. When nothing was
        stated this is ``None``, which is exactly the condition rule
        ``surface.family_ambiguous`` reports.
        """
        s = self.selection.get("surface")
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
    def has_circles(self):
        return bool(self.sd.get("circles"))

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
        """``"mat sheet, material 3 ('Core')"`` -- the template's own vocabulary."""
        try:
            name = self.materials[i].get("name")
        except (IndexError, AttributeError):
            name = None
        base = f"mat sheet, material {i + 1}"
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
        """``(label, points)`` for every distributed-load block of a stage."""
        key = "dloads" if stage == 1 else "dloads2"
        sheet = "dloads" if stage == 1 else "dloads (2)"
        out = []
        for n, blk in enumerate(self.sd.get(key) or []):
            pts = _coords(blk)
            out.append((f"{sheet} sheet, Distributed Load #{n + 1}", pts))
        return out

    @property
    def seep_bc(self):
        return self.sd.get("seepage_bc") or {}

    def bc_counts(self, which=1):
        """``(n_heads, n_fluxes, n_exit_points)`` for a boundary-condition set."""
        bc = self.sd.get("seepage_bc" if which == 1 else "seepage_bc2") or {}
        return (len(bc.get("specified_heads") or []),
                len(bc.get("specified_fluxes") or []),
                len(bc.get("exit_face") or []))


# ---------------------------------------------------------------------------
# Entry points
# ---------------------------------------------------------------------------

def preflight(slope_data, analysis, selection=None, ids=None, skip=None,
              include_availability=False):
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
        method name), ``surface`` (``"circular"`` / ``"noncircular"``), ``search``
        (bool), ``base`` (the underlying analysis for a sensitivity or reliability
        run).
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
    """
    ctx = _Ctx(slope_data, analysis, selection)
    only = set(ids) if ids else None
    drop = set(skip) if skip else set()
    findings = []
    for r in _REGISTRY:
        if only is not None and r.id not in only:
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
                           selection=dict(selection or {}))


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
        made.append(Finding(r.id, severity, str(message), r.remedy))
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
    >>> caps = capabilities(slope_data)                     # doctest: +SKIP
    >>> caps["lem_method"]["oms"].available                 # doctest: +SKIP
    False
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


def _circular_only_message(method, search):
    name = _METHOD_NAMES.get(method, method)
    valid = ", ".join(_METHOD_NAMES[k].replace("the ", "")
                      for k in LEM_METHOD_OPTIONS
                      if k not in CIRCULAR_ONLY_METHODS)
    if search:
        return (f"A non-circular search cannot be run with {name}: it takes moments "
                f"about a circle centre, which a non-circular surface does not have, "
                f"so every trial surface would be rejected. Methods valid for a "
                f"non-circular surface: {valid}.")
    return (f"{name} takes moments about a circle centre, so it cannot be used with "
            f"a non-circular surface. Methods valid for a non-circular surface: "
            f"{valid}.")


# ===========================================================================
# RULES
#
# Grouped by family. Every message names the sheet and the field in the
# template's vocabulary; every ERROR carries either a probe or a corpus instance
# behind it.
# ===========================================================================

# ---------------------------------------------------------------------------
# Family: water and unit weight of water
# ---------------------------------------------------------------------------

@rule("water.gamma_water_missing", ERROR, ("lem", "seep", "fem"),
      "main!D10 (Unit weight of water) must be present and positive.")
def _gamma_water_missing(ctx):
    g = ctx.sd.get("gamma_water")
    if _pos(g):
        return None
    # Only a model that actually uses water needs it -- a dry LEM deck does not.
    if "lem" in ctx.analyses and not (ctx.uses_piezo or ctx.uses_seep_u
                                      or ctx.sd.get("piezo_line")
                                      or ctx.sd.get("dloads")):
        return None
    return (f"main sheet D10 (Unit weight of water) is {_fmt(g)}. Enter a positive "
            f"unit weight of water -- there is deliberately no default, because "
            f"substituting a canonical value of the wrong system would rescale "
            f"every pore pressure in the model.")


@rule("piezo.line_missing", ERROR, ("lem",),
      "A material with u = piezo needs a piezometric line to read.",
      remedy="extend_piezo_line")
def _piezo_line_missing(ctx):
    if not ctx.uses_piezo:
        return None
    if len(ctx.piezo) >= 2:
        return None
    names = [ctx.mat_label(i) for i, m in ctx.strength_materials()
             if str(m.get("u", "")).strip().lower() == "piezo"]
    return (f"{names[0]} takes pore pressure from a piezometric line (u = piezo), "
            f"but the piezo sheet carries fewer than two points, so there is no "
            f"line to read. Every slice base in that material would read zero pore "
            f"pressure -- an unsafe, non-conservative factor of safety. Draw the "
            f"piezometric line, or set u = none for a dry model.")


@rule("piezo.extent_short", WARNING, ("lem",),
      "A piezometric line that stops short of the section refuses the slices past it.",
      remedy="extend_piezo_line")
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
    return (f"piezo sheet, Piezometric Line 1 spans x = [{pz[0]:.6g}, {pz[1]:.6g}], "
            f"which does not cover the whole section "
            f"([{gs[0]:.6g}, {gs[1]:.6g}]): it stops short at {', '.join(gaps)}. A "
            f"failure surface reaching into that reach is refused, because a "
            f"material with u = piezo would have no piezometric surface to read "
            f"there. This is fine when every trial surface stays inside the line.")


@rule("piezo.no_consumer", WARNING, ("lem",),
      "A piezometric line no material reads produces no pore pressure.")
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
    return ("The piezo sheet defines a piezometric line, but no material uses "
            "u = piezo, so the line produces no pore pressure anywhere and the "
            "model runs dry. Set u = piezo on the materials below the water "
            "table, or delete the line if the analysis is deliberately total "
            "stress.")


@rule("water.ponded_no_dload", WARNING, ("lem",),
      "Standing water above the ground surface needs a distributed load.",
      remedy="add_ponded_water_load")
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
    return (f"The piezometric line rises up to {depth:.3g} above the ground surface "
            f"between x = {a:.6g} and x = {b:.6g}, but the dloads sheet has no load "
            f"there. Standing water must be entered as a distributed load, or its "
            f"weight is missing from the analysis.")


# ---------------------------------------------------------------------------
# Family: material vocabulary and strength
# ---------------------------------------------------------------------------

@rule("mat.u_vocabulary", ERROR, ("*",),
      "mat!u must be one of none / piezo / seep / ru.")
def _mat_u_vocabulary(ctx):
    ok = ("none", "piezo", "seep", "ru")
    for i, m in enumerate(ctx.materials):
        u = m.get("u")
        if u is None or (isinstance(u, str) and not u.strip()):
            continue
        if str(u).strip().lower() not in ok:
            yield (f"{ctx.mat_label(i)}, column u: {u!r} is not a pore-pressure "
                   f"option. Expected one of: none, piezo, seep, ru. An unrecognized "
                   f"value used to fall through to zero pore pressure, silently "
                   f"inflating the factor of safety.")


@rule("mat.unsat_vocabulary", ERROR, ("seep", "fem"),
      "mat!unsat must be one of lf / vg / gard.")
def _mat_unsat_vocabulary(ctx):
    ok = ("lf", "vg", "gard")
    for i, m in ctx.seepage_materials():
        v = m.get("unsat")
        if v is None or (isinstance(v, str) and not v.strip()):
            continue
        if str(v).strip().lower() not in ok:
            yield (f"{ctx.mat_label(i)}, column unsat: {v!r} is not an unsaturated "
                   f"conductivity model. Expected one of: lf (linear front), vg "
                   f"(van Genuchten), gard (Gardner). A typo here silently runs the "
                   f"linear-front model instead of the one you asked for.")


@rule("mat.option_vocabulary", ERROR, ("lem", "fem"),
      "mat!option must be one of mc / cp / pow / hb / elastic.")
def _mat_option_vocabulary(ctx):
    ok = ("mc", "cp", "pow", "hb", "elastic")
    for i, m in ctx.strength_materials():
        o = m.get("option")
        if o is None or (isinstance(o, str) and not o.strip()):
            continue
        if str(o).strip().lower() not in ok:
            yield (f"{ctx.mat_label(i)}, column option: {o!r} is not a strength "
                   f"model. Expected one of: mc, cp, pow, hb, elastic.")


@rule("mat.option_missing", ERROR, ("lem",),
      "A material a failure surface can cross needs a strength model.")
def _mat_option_missing(ctx):
    for i, m in ctx.strength_materials():
        o = m.get("option")
        if o is None or (isinstance(o, str) and not str(o).strip()):
            yield (f"{ctx.mat_label(i)}, column option is blank, but the material is "
                   f"inside the geometry, so a failure surface can pass through it. "
                   f"Choose a strength model (mc, cp, pow, hb or elastic).")


@rule("mat.gamma_nonpositive", ERROR, ("lem",),
      "A material carrying strength data needs a positive unit weight.")
def _mat_gamma_nonpositive(ctx):
    for i, m in enumerate(ctx.materials):
        opt = str(m.get("option") or "").strip().lower()
        if opt not in ("mc", "cp", "pow", "hb"):
            continue
        if _pos(m.get("gamma")):
            continue
        yield (f"{ctx.mat_label(i)}, column g (unit weight) is {_fmt(m.get('gamma'))}. "
               f"The row carries a strength model (option = {opt}), so it is a "
               f"slope-stability material and needs a positive unit weight -- a "
               f"gamma of zero produces zero slice weights and a meaningless factor "
               f"of safety.")


@rule("mat.no_shear_strength", WARNING, ("lem",),
      "option = mc with c = 0 and phi = 0 is a material with no shear strength.")
def _mat_no_shear_strength(ctx):
    for i, m in ctx.strength_materials():
        if str(m.get("option") or "").strip().lower() != "mc":
            continue
        c, phi = _num(m.get("c")) or 0.0, _num(m.get("phi")) or 0.0
        if c == 0.0 and phi == 0.0:
            yield (f"{ctx.mat_label(i)}: option = mc with c = 0 and f = 0, so the "
                   f"material has no shear strength at all. If it is cohesionless, "
                   f"enter its friction angle; the sheet cannot tell "
                   f"'cohesionless' from 'no data'.")


@rule("mat.cp_zero_strength", ERROR, ("lem",),
      "option = cp with both c and c/p zero gives zero undrained strength.")
def _mat_cp_zero(ctx):
    for i, m in ctx.strength_materials():
        if str(m.get("option") or "").strip().lower() != "cp":
            continue
        c, cp = _num(m.get("c")) or 0.0, _num(m.get("cp")) or 0.0
        if c == 0.0 and cp == 0.0:
            yield (f"{ctx.mat_label(i)}: option = cp but both c and c/p are zero, so "
                   f"the undrained strength Su is zero everywhere in this material.")


@rule("mat.ru_zero", ERROR, ("lem",),
      "u = ru with a blank or zero ratio generates no pore pressure.")
def _mat_ru_zero(ctx):
    for i, m in ctx.strength_materials():
        if str(m.get("u") or "").strip().lower() != "ru":
            continue
        if not _pos(m.get("ru")):
            yield (f"{ctx.mat_label(i)} selects u = ru, but column ru is "
                   f"{_fmt(m.get('ru'))}. No pore pressure would be generated, which "
                   f"is identical to a dry model and non-conservative. Enter a "
                   f"positive pore pressure ratio, or set u = none.")


# ---------------------------------------------------------------------------
# Family: main-sheet scalars
# ---------------------------------------------------------------------------

@rule("main.seismic_missing", ERROR, ("lem",),
      "main!D13 (Seismic coefficient) must be a number.",
      remedy="set_seismic_zero")
def _seismic_missing(ctx):
    if _finite(ctx.sd.get("k_seismic")):
        return None
    return ("main sheet D13 (Seismic coefficient) is blank or not a number -- enter "
            "0 for a static analysis. Left blank it reaches every slice as NaN, and "
            "the Ordinary Method of Slices then reports success with a factor of "
            "safety of NaN.")


@rule("main.seismic_magnitude", WARNING, ("lem", "fem"),
      "main!D13 outside the usual 0.05-0.3 range is probably a percent/fraction slip.")
def _seismic_magnitude(ctx):
    k = _num(ctx.sd.get("k_seismic"))
    if k is None or k == 0.0:
        return None
    a = abs(k)
    if a > 1.0:
        return (f"main sheet D13 (Seismic coefficient) = {k:g}. k is a fraction of "
                f"gravity, typically 0.05-0.3 -- enter 0.15 for 15%g. At this "
                f"magnitude some methods return a negative factor of safety through "
                f"the success path and others fail without naming k.")
    if a > 0.5:
        return (INFO, f"main sheet D13 (Seismic coefficient) = {k:g}, which is high "
                      f"for a pseudo-static analysis (typically 0.05-0.3).")
    return None


@rule("main.seismic_negative_lem", WARNING, ("lem",),
      "The LEM takes the magnitude of k and orients it itself; the sign is ignored.")
def _seismic_negative_lem(ctx):
    k = _num(ctx.sd.get("k_seismic"))
    if k is None or k >= 0:
        return None
    return (f"main sheet D13 (Seismic coefficient) = {k:g} is negative. The "
            f"limit-equilibrium engine applies the seismic coefficient in the "
            f"failure-driving direction automatically and ignores the sign, so this "
            f"run is identical to k = {abs(k):g}. Enter a magnitude. (In the finite "
            f"element engine the sign does set the direction.)")


@rule("main.crack_water_deeper_than_crack", ERROR, ("lem",),
      "main!D12 (Depth of water in crack) cannot exceed main!D11 (Tension crack depth).")
def _crack_water(ctx):
    w = _num(ctx.sd.get("tcrack_water")) or 0.0
    d = _num(ctx.sd.get("tcrack_depth")) or 0.0
    if w <= 0:
        return None
    if d <= 0:
        return (f"main sheet D12 (Depth of water in crack) = {w:g}, but D11 (Tension "
                f"crack depth) is {_fmt(ctx.sd.get('tcrack_depth'))}. There is no "
                f"crack for the water to stand in, yet the full crack-water thrust "
                f"is applied to the surface.")
    if w > d + 1e-12:
        return (f"main sheet D12 (Depth of water in crack) = {w:g} exceeds D11 "
                f"(Tension crack depth) = {d:g}. The water cannot be deeper than "
                f"the crack that holds it.")
    return None


@rule("main.lem_method_unknown", ERROR, ("lem",),
      "The named LEM method must be one the package implements.")
def _lem_method_unknown(ctx):
    m = ctx.method
    if not m or m == "all":
        return None
    if m in LEM_METHOD_OPTIONS:
        return None
    return (f"main sheet D14 (LEM method) is {m!r}, which is not a method this "
            f"package implements. Choose one of: "
            f"{', '.join(LEM_METHOD_OPTIONS)} (or 'all').")


# ---------------------------------------------------------------------------
# Family: surface family and method compatibility
# ---------------------------------------------------------------------------

@rule("surface.none_defined", ERROR, ("lem",), capability="analysis",
      summary="A limit-equilibrium run needs a circular or non-circular surface.")
def _surface_none(ctx):
    if ctx.has_circles or ctx.has_non_circ:
        return None
    return ("This model defines no failure surface: the circles sheet and the "
            "non-circ sheet are both empty. A limit-equilibrium run needs at least "
            "one starting surface.")


@rule("surface.method_requires_circle", ERROR, ("lem",), capability="lem_method",
      summary="OMS and Bishop take moments about a circle centre (circular family only).")
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
    return ("This model defines both a circular surface (circles sheet) and a "
            "non-circular surface (non-circ sheet), and the run did not state which "
            "to analyse. The circular surface will be used; the non-circular "
            "surface will be ignored.")


# ---------------------------------------------------------------------------
# Family: polyline ordering (the remedy path depends on it)
# ---------------------------------------------------------------------------

@rule("order.piezo_reversed", WARNING, ("lem",),
      "A piezometric line entered right-to-left is sorted before use.",
      remedy="reverse_polyline")
def _piezo_reversed(ctx):
    out = []
    for label, pts in (("Piezometric Line 1", ctx.piezo),
                       ("Piezometric Line 2", ctx.piezo2)):
        if len(pts) >= 2 and _monotonic_x(pts) == "descending":
            out.append(f"piezo sheet, {label}: the points run right to left. The "
                       f"slicer sorts them before use, so the analysis is "
                       f"unaffected, but check the line is what you intended.")
    return out


@rule("order.piezo_nonmonotonic", ERROR, ("lem",),
      "A piezometric line whose x values rise then fall is not a surface.")
def _piezo_nonmonotonic(ctx):
    out = []
    for label, pts in (("Piezometric Line 1", ctx.piezo),
                       ("Piezometric Line 2", ctx.piezo2)):
        if len(pts) >= 2 and _monotonic_x(pts) == "mixed":
            out.append(f"piezo sheet, {label}: the X values are not in order (they "
                       f"rise, then fall). Sorting them would produce a different "
                       f"line from the one drawn -- enter the points along the water "
                       f"surface in one direction.")
    return out


@rule("order.dload_reversed", WARNING, ("lem", "fem"),
      "A distributed-load block entered right-to-left is repaired silently at load.",
      remedy="reverse_polyline")
def _dload_reversed(ctx):
    out = []
    for stage in (1, 2):
        for label, pts in ctx.dload_blocks(stage):
            if len(pts) >= 2 and _monotonic_x(pts) == "descending":
                out.append(f"{label}: the points run right to left. They are sorted "
                           f"before use, but check the block is what you intended.")
    return out


@rule("order.dload_nonmonotonic", ERROR, ("lem", "fem"),
      "A distributed-load block whose x values rise then fall cannot be interpolated.")
def _dload_nonmonotonic(ctx):
    out = []
    for stage in (1, 2):
        for label, pts in ctx.dload_blocks(stage):
            if len(pts) >= 2 and _monotonic_x(pts) == "mixed":
                out.append(f"{label}: the X values are not in order (they rise, then "
                           f"fall). Enter the points along the loaded surface in one "
                           f"direction.")
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
    return (f"main sheet D8 declares {system!r} units, whose unit weight of water is "
            f"{canonical:g}, but D10 carries {gw:g} "
            f"({abs(gw - canonical) / canonical * 100:.1f}% off). If that is a "
            f"deliberate seawater or brine value this is fine; otherwise check the "
            f"unit system -- xslope never converts, so a mislabelled system means "
            f"every number is read in the wrong one.")


@rule("units.material_gamma_off_band", WARNING, ("*",),
      "Soil unit weights sit in far-apart bands in the two systems.")
def _units_material_gamma(ctx):
    from .units import infer_unit_system
    system = _declared_system(ctx)
    if system is None:
        return None
    inferred = infer_unit_system(ctx.materials)
    if inferred is None or inferred == system:
        return None
    return (f"main sheet D8 declares {system!r} units, but the material unit weights "
            f"on the mat sheet look {inferred!r} (SI soils are about 15-25 kN/m3, "
            f"Imperial soils about 90-160 pcf). xslope does not convert, so a "
            f"mislabelled system means the numbers are being read in the wrong one.")


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
    return (f"main sheet D8 declares no unit system. The unit weight of water "
            f"({_fmt(ctx.sd.get('gamma_water'))}) looks {from_water!r} while the "
            f"material unit weights look {from_soil!r}. Two independent signals "
            f"disagree, so one set of numbers is in the wrong system -- xslope "
            f"never converts.")


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


@rule("mat.E_unusable", WARNING, ("fem",),
      "A blank or non-positive Young's modulus reaches the solver as a singular matrix.")
def _mat_E_unusable(ctx):
    for i, m in ctx.strength_materials():
        raw = m.get("E")
        if raw is None or (isinstance(raw, str) and not str(raw).strip()):
            yield (f"{ctx.mat_label(i)}, column E: Young's modulus is blank. The "
                   f"finite element engine needs a positive modulus for every "
                   f"material it meshes.")
            continue
        E = _num(raw)
        if E is None or E <= 0:
            yield (f"{ctx.mat_label(i)}, column E: Young's modulus is {_fmt(raw)}. "
                   f"The finite element engine needs a positive modulus for every "
                   f"material it meshes; a zero or non-finite E reaches the solver "
                   f"as 'Factor is exactly singular'.")


# ---------------------------------------------------------------------------
# Family: mesh and stored seepage fields
# ---------------------------------------------------------------------------

@rule("mesh.missing", ERROR, ("seep", "fem"), capability="analysis",
      availability=True,
      summary="A seepage or finite element run needs a mesh.")
def _mesh_missing(ctx):
    if ctx.mesh is not None:
        return None
    return ("This model carries no finite element mesh. Build the mesh first "
            "(a '{base}_mesh.json' file beside the .xlsx, or Build Mesh in Studio) "
            "-- a seepage or finite element run has nothing to solve on without it.")


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
            f"Rebuild the mesh with a supported element type.")


@rule("mesh.material_id_out_of_range", ERROR, ("seep", "fem"),
      "A mesh element cannot reference a material the mat sheet does not define.")
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
    return (f"The mesh references material ID(s) {', '.join(map(str, bad))}, but the "
            f"mat sheet defines {n} material(s). The mesh was built against a "
            f"different material table -- rebuild it after editing the mat sheet.")


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
    for key, sheet in (("seep_u", "{base}_seep.csv"),
                       ("seep_u2", "{base}_seep2.csv")):
        u = ctx.sd.get(key)
        if u is None:
            continue
        try:
            n_u = len(u)
        except TypeError:
            continue
        if n_u != n_nodes:
            out.append(
                f"The stored pore-pressure field '{sheet}' carries {n_u} node "
                f"value(s) but the mesh has {n_nodes} node(s), so it was computed on "
                f"a different mesh. Re-run the seepage analysis on this mesh.")
    return out


@rule("seep_field.missing", ERROR, ("lem",),
      "A material with u = seep needs a solved pore-pressure field.")
def _seep_field_missing(ctx):
    if not ctx.uses_seep_u:
        return None
    if ctx.sd.get("seep_u") is not None and ctx.mesh is not None:
        return None
    names = [ctx.mat_label(i) for i, m in ctx.strength_materials()
             if str(m.get("u", "")).strip().lower() == "seep"]
    missing = "mesh" if ctx.mesh is None else "pore-pressure field"
    return (f"{names[0]} takes pore pressure from a seepage solution (u = seep), but "
            f"this model carries no {missing}. Run the seepage analysis first, or "
            f"choose a different pore-pressure option -- every slice base in that "
            f"material would otherwise read zero pore pressure.")


# ---------------------------------------------------------------------------
# Family: steady seepage -- material properties
# ---------------------------------------------------------------------------

@rule("seep.k1_nonpositive", ERROR, ("seep",),
      "mat!k1 (major conductivity) must be greater than zero.")
def _seep_k1(ctx):
    for i, m in ctx.seepage_materials():
        if _pos(m.get("k1")):
            continue
        yield (f"{ctx.mat_label(i)}, column k1: hydraulic conductivity is "
               f"{_fmt(m.get('k1'))} and must be greater than zero. Model an "
               f"impermeable zone with a small conductivity (for example 1e-6 times "
               f"the largest k), never with zero -- a zero makes the conductivity "
               f"matrix singular and destroys the flow net while the head field "
               f"still looks plausible.")


@rule("seep.k2_nonpositive", ERROR, ("seep",),
      "mat!k2 (minor conductivity) must be greater than zero.")
def _seep_k2(ctx):
    for i, m in ctx.seepage_materials():
        if _pos(m.get("k2")):
            continue
        yield (f"{ctx.mat_label(i)}, column k2: hydraulic conductivity is "
               f"{_fmt(m.get('k2'))} and must be greater than zero. For an isotropic "
               f"material set k2 = k1.")


@rule("seep.k2_greater_than_k1", WARNING, ("seep",),
      "k1 is the MAJOR principal conductivity; k2 > k1 rotates the ellipse 90 degrees.")
def _seep_k2_gt_k1(ctx):
    for i, m in ctx.seepage_materials():
        k1, k2 = _num(m.get("k1")), _num(m.get("k2"))
        if k1 is None or k2 is None or k1 <= 0 or k2 <= 0:
            continue
        if k2 > k1:
            yield (f"{ctx.mat_label(i)}: k2 = {k2:g} is greater than k1 = {k1:g}. "
                   f"k1 is the MAJOR principal conductivity, so this silently "
                   f"rotates the conductivity ellipse by 90 degrees -- swap the two "
                   f"values, or set alpha = 90.")


@rule("seep.anisotropy_angle_unset", INFO, ("seep",),
      "An anisotropic material with alpha blank puts k1 along +x.")
def _seep_alpha(ctx):
    for i, m in ctx.seepage_materials():
        k1, k2 = _num(m.get("k1")), _num(m.get("k2"))
        if k1 is None or k2 is None or k1 <= 0 or k2 <= 0 or k1 == k2:
            continue
        if (_num(m.get("alpha")) or 0.0) == 0.0:
            yield (f"{ctx.mat_label(i)} is anisotropic (k1 = {k1:g}, k2 = {k2:g}) "
                   f"with alpha blank, so the major conductivity is taken along +x.")


@rule("seep.unsat_params_missing", ERROR, ("seep",),
      "An unconfined model needs the unsaturated parameters its unsat model uses.")
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
                    f"{ctx.mat_label(i)} uses unsat = lf on an unconfined model "
                    f"(the seep bc sheet defines an exit face), so it needs kr0 > 0 "
                    f"and h0 < 0; it carries kr0 = {_fmt(m.get('kr0'))}, "
                    f"h0 = {_fmt(m.get('h0'))}.")
        elif model in ("vg", "gard"):
            a, n = _num(m.get("vg_a")), _num(m.get("vg_n"))
            need = "n > 1" if model == "vg" else "n > 0"
            bad = (a is None or a <= 0 or n is None
                   or (n <= 1 if model == "vg" else n <= 0))
            if bad:
                out.append(
                    f"{ctx.mat_label(i)} uses unsat = {model} on an unconfined "
                    f"model, so it needs a > 0 and {need}; it carries "
                    f"a = {_fmt(m.get('vg_a'))}, n = {_fmt(m.get('vg_n'))}.")
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
    return (f"The seep bc sheet defines no exit face, so this solve is confined and "
            f"the unsaturated inputs (unsat, kr0, h0, a, n) are not used -- "
            f"including on {named[0]}.")


# ---------------------------------------------------------------------------
# Family: steady seepage -- boundary conditions
# ---------------------------------------------------------------------------

@rule("seep.no_boundary_conditions", ERROR, ("seep",), capability="analysis",
      summary="A seepage run needs at least one specified head or an exit face.")
def _seep_no_bcs(ctx):
    n_head, n_flux, n_exit = ctx.bc_counts(1)
    if n_head or n_flux or n_exit:
        return None
    return ("The seep bc sheet defines no boundary conditions at all. A seepage "
            "analysis needs at least one specified head, or an exit face. With none "
            "the solve converges to zero head everywhere and reports success, which "
            "reaches a stability run as zero pore pressure.")


@rule("seep.no_dirichlet", ERROR, ("seep",),
      "A flux-only model defines head only up to an additive constant.")
def _seep_no_dirichlet(ctx):
    n_head, n_flux, n_exit = ctx.bc_counts(1)
    if n_head or n_exit or not n_flux:
        return None
    return ("The seep bc sheet defines specified-flux boundary conditions but no "
            "specified head and no exit face. The head is then defined only up to an "
            "additive constant and the system is singular -- add a specified head, "
            "or an exit face.")


@rule("seep.exit_face_without_head", ERROR, ("seep",),
      "An exit face alone has nothing driving flow through it.")
def _seep_exit_only(ctx):
    n_head, n_flux, n_exit = ctx.bc_counts(1)
    if not n_exit or n_head or n_flux:
        return None
    return ("The seep bc sheet defines an exit face but no specified head and no "
            "flux. Once the exit face deactivates there is no boundary left to fix "
            "the head, and the solve becomes singular. Add the reservoir or "
            "tailwater head that drives the flow.")


@rule("seep.no_gradient", WARNING, ("seep",),
      "Two heads of the same value, or one head alone, drive no flow.")
def _seep_no_gradient(ctx):
    bc = ctx.seep_bc
    heads = bc.get("specified_heads") or []
    n_flux = len(bc.get("specified_fluxes") or [])
    n_exit = len(bc.get("exit_face") or [])
    if n_flux or n_exit or not heads:
        return None
    if ctx.sd.get("tseep") is not None:
        return None            # a transient run drives the heads from its series
    if any(isinstance(b.get("head"), str) for b in heads):
        return None            # bound to a tseep series, so the value varies in time
    values = {round(float(b.get("head")), 12) for b in heads
              if _finite(b.get("head"))}
    if len(values) >= 2:
        return None
    return ("The seep bc sheet defines a confined model whose specified heads all "
            "carry the same value, so there is no gradient and no flow. A confined "
            "model needs at least two distinct head values, or a flux boundary.")


@rule("seep.bc_polyline_too_short", ERROR, ("seep",),
      "A boundary-condition polyline needs at least two points to bind mesh nodes.")
def _seep_bc_short(ctx):
    bc = ctx.seep_bc
    out = []
    for n, b in enumerate(bc.get("specified_heads") or []):
        if len(_coords(b.get("coords"))) < 2:
            out.append(f"seep bc sheet, specified head #{n + 1} has fewer than two "
                       f"points, so it binds no mesh nodes and is silently skipped.")
    for n, b in enumerate(bc.get("specified_fluxes") or []):
        if len(_coords(b.get("coords"))) < 2:
            out.append(f"seep bc sheet, specified flux #{n + 1} has fewer than two "
                       f"points, so it binds no mesh nodes.")
    exit_face = bc.get("exit_face") or []
    if 0 < len(_coords(exit_face)) < 2:
        out.append("seep bc sheet, Exit Face has fewer than two points, so it binds "
                   "no mesh nodes and the model would silently solve as fully "
                   "confined.")
    return out
