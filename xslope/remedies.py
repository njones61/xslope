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

"""Remedies -- the fixes preflight can offer for what it finds.

A preflight rule says what is wrong. Some of those faults have exactly one
sensible repair, and this module computes it. The governing principle is one
line, and everything else here follows from it:

> **A remedy is always offered and never applied silently.**

A fix nobody asked for is the silent-default disease wearing a helpful face. It
would be worse than the fault it repairs, because the model would then disagree
with the file the user typed and nothing would say so.

The contract
------------

**Two surfaces, both explicit.** In a script a remedy is named, one at a time::

    report = preflight(slope_data, "lem", remedies=["add_ponded_water_load"])
    report.applied            # what was applied, and to which target
    report.model              # the remedied model -- the ORIGINAL is untouched

In an interface it is a button beside the finding. There is no "fix everything"
switch in either place.

**The description is computed before the change.** :func:`propose` returns what a
remedy *would* do -- *"Add 1 distributed-load block, 262 peak, over
x = 346.7 to 594"* -- while the user is still deciding. Applying runs the same
computation, so the label on the button and the change that lands cannot drift
apart.

**Applying produces a finding, not silence.** The ERROR or WARNING is replaced by
an INFO recording what changed and which rule asked for it, and that INFO travels
with the report. A model that ran on a synthesised load says so wherever its
factor of safety is reported.

**Remedies are pure and total.** ``slope_data -> slope_data`` on a copy, no file
writes, and never a partial application. A remedy that cannot fully succeed
declines and leaves the finding standing.

**Preconditions are checked, not assumed.** A remedy that depends on another rule
already being satisfied says so and refuses. The water-load derivation is the
example: it measures along lines whose x values increase, so on a line entered
right to left it declines and names the reversal remedy rather than quietly
deriving nothing.

**A remedy is capability-gated, never clickable into an error.** Every proposal
carries ``available`` together with a ``reason``, so a button dims with a tooltip
naming what is missing rather than being live and failing when pressed -- the same
behaviour, and the same reason string, as the :func:`xslope.preflight.capabilities`
layer that dims an invalid analysis method. :func:`remedy_capabilities` returns
exactly that shape.

What is implemented
-------------------

``reverse_polyline``
    Reverse a piezometric line or a distributed-load block entered right to left.
    Declines on a line whose x values rise and then fall: that is not a reversed
    line, and reversing it would produce a different shape rather than the one
    drawn.
``add_ponded_water_load``
    Add the standing water as a distributed load, derived from the model's own
    water definition through :mod:`xslope.water` -- seepage boundary conditions
    where a seepage analysis is defined, otherwise the stage's piezometric line.
``switch_to_auto_water``
    Hand the water load to the engine: set ``main!D23`` to ``auto`` and remove the
    transcribed blocks the derivation reproduces, keeping every block that is not
    water. It is the better of the two water remedies wherever it applies, because
    writing blocks into a sheet records a snapshot that goes stale the moment the
    pool moves, while the mode is recomputed at every solve. It declines rather
    than removing a block the derivation does not reproduce: that mismatch is a
    finding about the file, and a remedy is not the place to settle it.
``generate_starting_circles``
    Fill an empty circles sheet with a starting set derived from the slope
    geometry (:mod:`xslope.generators`).
``generate_noncircular_surface``
    Fill an empty non-circ sheet with a surface that tracks the model's weak
    zone. The second answer to an empty-surface model, and the right one where a
    weak layer -- not the slope's own geometry -- controls the mechanism, since no
    circle can follow a seam. Offered only where the generator picks a zone on its
    own: a model with two comparable candidates needs a question asked, and a
    remedy states its change and applies it, it does not open a dialog. Studio's
    zone picker is that path.

:data:`xslope.preflight.REMEDIES` also names remedies whose transformation has not
been built. Those are declarations of intent, and asking for one raises rather
than doing nothing: a remedy that silently no-ops is the failure this module
exists to prevent.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Optional

from . import preflight as _pf


class RemedyDeclined(ValueError):
    """Raised when a remedy is asked for and cannot fully succeed.

    Carries the same sentence the proposal's ``reason`` does, so a script and a
    dimmed button explain a refusal identically.
    """


@dataclass(frozen=True)
class RemedyProposal:
    """What a remedy would do, computed before anything is changed.

    Attributes
    ----------
    remedy : str
        The remedy's name, a key of :data:`xslope.preflight.REMEDIES`.
    target : str
        Which object it would change (``"stage1"``, ``"piezo1"``,
        ``"dload:1:2"``). Empty for a remedy that changes the model as a whole.
    key : str
        ``"<remedy>:<target>"`` -- the stable handle a button binds to and a
        script passes back to :func:`apply_remedy`.
    label : str
        The target in the template's own vocabulary
        (``"piezo sheet, Piezometric Line 1"``).
    available : bool
        False when the remedy cannot fully succeed on this model.
    reason : str
        Why it cannot -- user-facing copy, ready to use as a tooltip on the dimmed
        button. Empty when the remedy is available.
    description : str
        Exactly what applying would do, stated in the template's vocabulary.
    rule_ids : tuple of str
        The rules that offer this remedy, so a caller can tie a button back to the
        finding it belongs beside.
    """

    remedy: str
    target: str
    key: str
    label: str
    available: bool
    reason: str = ""
    description: str = ""
    rule_ids: tuple = ()
    _apply: Optional[Callable] = field(default=None, repr=False, compare=False)

    def __str__(self):
        head = self.description if self.available else f"unavailable: {self.reason}"
        return f"{self.key}: {head}"


def rule_ids_for(remedy):
    """The ids of every rule that offers ``remedy``.

    Read from the registry rather than restated here, so a rule that starts or
    stops offering a remedy cannot leave a stale list behind. A rule offering
    several answers to one fault is matched on all of them, not only on the primary
    the finding carries.
    """
    return tuple(r.id for r in _pf.rules() if remedy in r.remedies)


# ---------------------------------------------------------------------------
# Small shared helpers
# ---------------------------------------------------------------------------

def _coords(obj):
    from .water import _as_coords
    return _as_coords(obj)


def _order(pts):
    """``"ascending"``, ``"descending"`` or ``"mixed"`` for a polyline's x values.

    The same reading the ordering rules take, from the same helper, so a rule that
    fires and the remedy beside it can never disagree about which case this is.
    """
    return _pf._monotonic_x(pts)


def _shallow(slope_data):
    """A copy safe to write the touched sheets into.

    The model's other entries -- geometry, materials, a mesh -- are shared with the
    original by reference and are never mutated: a remedy replaces whole lists
    rather than editing them in place, so the caller's model comes back unchanged
    whatever happens.
    """
    return dict(slope_data or {})


# ---------------------------------------------------------------------------
# Remedy: reverse a backwards polyline
# ---------------------------------------------------------------------------

def _reverse_line_targets(sd):
    """Every piezometric line and load block whose points do not run left to right."""
    out = []
    for key, target, label in (("piezo_line", "piezo1", "piezo sheet, Piezometric Line 1"),
                               ("piezo_line2", "piezo2", "piezo sheet, Piezometric Line 2")):
        pts = _coords(sd.get(key))
        if len(pts) >= 2 and _order(pts) != "ascending":
            out.append((target, label, pts, key, None))
    for stage, key, sheet in ((1, "dloads", "dloads"), (2, "dloads2", "dloads (2)")):
        for n, blk in enumerate(sd.get(key) or []):
            pts = _coords(blk)
            if len(pts) >= 2 and _order(pts) != "ascending":
                out.append((f"dload:{stage}:{n + 1}",
                            f"{sheet} sheet, Distributed Load #{n + 1}",
                            pts, key, n))
    return out


def _propose_reverse(sd):
    proposals = []
    rules = rule_ids_for("reverse_polyline")
    for target, label, pts, key, index in _reverse_line_targets(sd):
        if _order(pts) == "mixed":
            proposals.append(RemedyProposal(
                remedy="reverse_polyline", target=target,
                key=f"reverse_polyline:{target}", label=label, available=False,
                reason=(f"{label} does not simply run backwards: its X values rise "
                        f"and then fall. Reversing it would produce a different line "
                        f"from the one drawn -- re-enter the points along the surface "
                        f"in one direction."),
                rule_ids=rules))
            continue
        proposals.append(RemedyProposal(
            remedy="reverse_polyline", target=target,
            key=f"reverse_polyline:{target}", label=label, available=True,
            description=(f"Reverse the {len(pts)} points of {label}, so they run "
                         f"left to right from x = {pts[-1][0]:.6g} to "
                         f"x = {pts[0][0]:.6g}."),
            rule_ids=rules,
            _apply=lambda model, key=key, index=index: _apply_reverse(model, key, index)))
    return proposals


def _apply_reverse(sd, key, index):
    """Reverse one line, keeping the shape the model stores it in.

    A piezometric line is a list of coordinate pairs on one model and a shapely
    line on another; handing back the other form would be a second, unannounced
    change to the model.
    """
    out = _shallow(sd)
    if index is None:
        obj = sd.get(key)
        if hasattr(obj, "coords"):
            from shapely.geometry import LineString
            out[key] = LineString(list(reversed(list(obj.coords))))
        else:
            out[key] = list(reversed(list(obj)))
    else:
        blocks = list(sd.get(key) or [])
        blocks[index] = list(reversed(blocks[index]))
        out[key] = blocks
    return out


# ---------------------------------------------------------------------------
# Remedy: add the ponded water as a distributed load
# ---------------------------------------------------------------------------

def _is_staged(sd):
    """True when the model carries a second stage -- a rapid-drawdown deck."""
    bc2 = sd.get("seepage_bc2") or {}
    return bool(_coords(sd.get("piezo_line2")) or sd.get("dloads2")
                or bc2.get("specified_heads") or bc2.get("specified_fluxes")
                or bc2.get("exit_face"))


def _stage_sheets(stage):
    if stage == 1:
        return "dloads", "dload_dirs", "dloads"
    return "dloads2", "dload2_dirs", "dloads (2)"


def _propose_ponded(sd):
    from .water import derive_water_loads
    rules = rule_ids_for("add_ponded_water_load")
    proposals = []
    mode = str(sd.get("water_loads") or "manual").strip().lower()
    for stage in (1, 2):
        key, _dirs_key, sheet = _stage_sheets(stage)
        if stage == 2 and not _is_staged(sd):
            # A model with no second stage at all has no stage-2 sheet to be missing
            # a load from; offering the button there would be answering a question
            # nobody asked.
            continue
        target = f"stage{stage}"
        base = dict(remedy="add_ponded_water_load", target=target,
                    key=f"add_ponded_water_load:{target}",
                    label=f"{sheet} sheet", rule_ids=rules)
        if mode == "auto":
            # Nothing to write into: the engine derives the load at solve time, and
            # a block written here would be counted a second time.
            proposals.append(RemedyProposal(
                available=False,
                reason=("Water loads are set to auto on the main sheet, so the engine "
                        "derives this load itself -- adding it to the sheet would count "
                        "the water twice."), **base))
            continue
        derived = derive_water_loads(sd, stage=stage)
        if not derived["blocks"]:
            if not derived["reason"]:
                continue
            proposals.append(RemedyProposal(available=False, reason=derived["reason"],
                                            **base))
            continue
        a, b = derived["reach"]
        overlapping = _overlapping_blocks(sd.get(key) or [], a, b)
        if overlapping:
            proposals.append(RemedyProposal(
                available=False,
                reason=(f"The {sheet} sheet already carries "
                        f"{'a load' if len(overlapping) == 1 else 'loads'} over "
                        f"x = {a:.6g} to {b:.6g} (Distributed Load "
                        f"{', '.join('#%d' % i for i in overlapping)}). Adding the "
                        f"derived water there would count the reservoir twice; check "
                        f"the existing block instead."), **base))
            continue
        n = len(derived["blocks"])
        proposals.append(RemedyProposal(
            available=True,
            description=(f"Add {n} distributed-load block{'' if n == 1 else 's'} to the "
                         f"{sheet} sheet, {derived['peak']:.6g} peak, over "
                         f"x = {a:.6g} to {b:.6g}, derived from {derived['source']}."),
            _apply=lambda model, stage=stage: _apply_ponded(model, stage), **base))
    return proposals


def _overlapping_blocks(blocks, a, b):
    """1-based indices of existing blocks whose x-range overlaps ``[a, b]``."""
    hits = []
    for n, blk in enumerate(blocks):
        pts = _coords(blk)
        if not pts:
            continue
        lo, hi = min(p[0] for p in pts), max(p[0] for p in pts)
        if lo <= b and hi >= a:
            hits.append(n + 1)
    return hits


def _apply_ponded(sd, stage):
    from .water import derive_water_loads
    key, dirs_key, _sheet = _stage_sheets(stage)
    derived = derive_water_loads(sd, stage=stage)
    if not derived["blocks"]:
        raise RemedyDeclined(derived["reason"] or "no water load could be derived")
    out = _shallow(sd)
    out[key] = list(sd.get(key) or []) + [list(b) for b in derived["blocks"]]
    dirs = list(sd.get(dirs_key) or [])
    while len(dirs) < len(sd.get(key) or []):
        dirs.append("normal")
    out[dirs_key] = dirs + list(derived["dirs"])
    return out


# ---------------------------------------------------------------------------
# Remedy: switch a model to automatic water loads
# ---------------------------------------------------------------------------

def _water_change(sd):
    """What switching this model to auto would do, per stage, or why it cannot.

    Returns ``(changes, reason)``. ``changes`` is ``{stage: (keep, drop, derived,
    match)}`` -- the blocks that survive, the transcribed water blocks the
    derivation replaces, and the measured agreement between them.

    The whole judgement is one question asked per stage: does the derivation
    reproduce what was transcribed? Where it does, the transcribed block is the
    reservoir and the mode change is an exact substitution. Where it does not, the
    remedy declines rather than removing a block whose replacement would be a
    different load -- that mismatch is a finding about the file, and quietly
    resolving it in the direction of the derivation is precisely what a remedy
    must not do.
    """
    from .water import derive_water_loads, match_water_blocks
    changes = {}
    for stage in (1, 2):
        key, _dirs_key, sheet = _stage_sheets(stage)
        blocks = list(sd.get(key) or [])
        derived = derive_water_loads(sd, stage=stage)
        found = match_water_blocks(blocks, derived["blocks"])
        if derived["blocks"] and found["missing"]:
            # Water the engine would add that nothing in the sheet stands for. The
            # switch is still safe when the sheet is empty (there was nothing to
            # double-count), but not when it carries something over the same reach.
            over = [i for i in found["unmatched"]
                    if _overlapping_blocks([blocks[i]], *derived["reach"])]
            if over:
                return None, (
                    f"The {sheet} sheet carries a load over the pool's reach "
                    f"(x = {derived['reach'][0]:.6g} to {derived['reach'][1]:.6g}) "
                    f"that the derived water load does not reproduce. Switching to "
                    f"auto would either double it or replace it with a different "
                    f"load, so the difference has to be resolved first: check "
                    f"Distributed Load #{over[0] + 1} against "
                    f"{derived['source']}.")
        changes[stage] = ([b for i, b in enumerate(blocks) if i in found["unmatched"]],
                          [i for i, _j, _m in found["pairs"]], derived, found)
    if not any(c[1] for c in changes.values()):
        blocks1 = list(sd.get("dloads") or [])
        derived1 = changes[1][2]
        if blocks1 and not derived1["blocks"]:
            # The case that must not be waved through: the sheet carries loads and
            # the engine would derive nothing, so switching would leave the model
            # with whatever those blocks are and no water of its own. Name why the
            # derivation is empty -- that is the thing to fix, if anything is.
            why = derived1["reason"] or "the derivation produced nothing"
            return None, (
                f"Switching to auto would derive no water load at all for this "
                f"model: {why}. The {len(blocks1)} block(s) on the dloads sheet "
                f"would be left standing as ordinary loads, so the mode change "
                f"would not express what the file already says.")
        return None, ("This model has no transcribed water load for the derivation "
                      "to take over: setting D23 to auto would add the standing "
                      "water rather than replace anything, which is a change to the "
                      "analysis rather than a change of expression.")
    return changes, ""


def _propose_switch_to_auto(sd):
    rules = rule_ids_for("switch_to_auto_water")
    base = dict(remedy="switch_to_auto_water", target="water_loads",
                key="switch_to_auto_water:water_loads",
                label="main sheet, D23 (Water loads)", rule_ids=rules)
    if str(sd.get("water_loads") or "manual").strip().lower() == "auto":
        return []                       # already automatic: nothing to switch
    changes, reason = _water_change(sd)
    if changes is None:
        return [RemedyProposal(available=False, reason=reason, **base)]
    parts = []
    for stage, (keep, drop, derived, found) in sorted(changes.items()):
        if not drop:
            continue
        _key, _dirs, sheet = _stage_sheets(stage)
        worst = found["worst"]
        parts.append(
            f"remove {len(drop)} block{'' if len(drop) == 1 else 's'} from the "
            f"{sheet} sheet ("
            + ", ".join(f"#{i + 1}" for i in drop)
            + f"), which the derivation from {derived['source']} reproduces to "
              f"within {100 * worst:.2g}% "
              f"(resultant {derived['resultant']:.6g})"
            + (f", keeping {len(keep)} other block{'' if len(keep) == 1 else 's'}"
               if keep else ""))
    return [RemedyProposal(
        available=True,
        description=("Set main!D23 (Water loads) to auto and " + "; and ".join(parts)
                     + ". The engine then derives the water load from the model's "
                       "water definition at every solve, so it can never go stale "
                       "when the pool moves."),
        _apply=_apply_switch_to_auto, **base)]


def _apply_switch_to_auto(sd):
    changes, reason = _water_change(sd)
    if changes is None:
        raise RemedyDeclined(reason)
    out = _shallow(sd)
    out["water_loads"] = "auto"
    for stage, (keep, drop, _derived, _found) in changes.items():
        key, dirs_key, _sheet = _stage_sheets(stage)
        blocks = list(sd.get(key) or [])
        dirs = list(sd.get(dirs_key) or [])
        while len(dirs) < len(blocks):
            dirs.append("normal")
        out[key] = keep
        out[dirs_key] = [d for i, d in enumerate(dirs) if i not in set(drop)]
    return out


# ---------------------------------------------------------------------------
# Remedy: generate a starting set of circles
# ---------------------------------------------------------------------------

def _propose_circles(sd):
    from .generators import generate_starting_circles
    rules = rule_ids_for("generate_starting_circles")
    base = dict(remedy="generate_starting_circles", target="circles",
                key="generate_starting_circles:circles", label="circles sheet",
                rule_ids=rules)
    if sd.get("circles"):
        return []                       # there is already a starting surface
    result = generate_starting_circles(sd, report=True)
    if not result["circles"]:
        return [RemedyProposal(available=False, reason=result["reason"], **base)]
    n = len(result["circles"])
    return [RemedyProposal(
        available=True,
        description=(f"Add {n} starting circle{'' if n == 1 else 's'} to the circles "
                     f"sheet: {result['summary']}."),
        _apply=_apply_circles, **base)]


def _apply_circles(sd):
    from .generators import generate_starting_circles
    circles = generate_starting_circles(sd)
    if not circles:
        raise RemedyDeclined("no starting circle could be derived from the geometry")
    out = _shallow(sd)
    out["circles"] = list(sd.get("circles") or []) + circles
    out["circular"] = True
    return out


# ---------------------------------------------------------------------------
# Remedy: generate a non-circular surface along the weak zone
# ---------------------------------------------------------------------------

def _propose_noncirc(sd):
    """Offer the weak-zone surface, but only where the generator picked a zone itself.

    An empty-surface model has two repairs, not one, and which is right depends on
    what controls the mechanism. Where the slope's own geometry does, a starting set
    of circles is the answer; where a weak layer does, no circle can follow the seam
    and the tracked polyline is the only usable seed. Both are therefore offered on
    ``surface.none_defined``, and the model decides which is live.

    The gate is the generator's own confidence, not a second opinion computed here.
    :func:`~xslope.generators.generate_noncircular_surface` picks a zone without
    asking only when the weakest is at or under ``_SEPARATION`` times the next
    weakest; where two zones are comparable it returns the ranked candidates and no
    surface. That case is not offered AT ALL -- not offered-and-dimmed -- because
    nothing is wrong with the model and there is nothing for a button to say. The
    honest response to two comparable seams is to ask which one, and asking is a
    dialog: Studio's zone picker. A remedy states its exact change before applying
    it, so a remedy that had to raise a question first would be a different thing
    wearing this one's contract.

    The single-zone model is silent for the same reason: with one material the
    mechanism follows the slope's geometry, which is what the circles remedy beside
    this one produces.

    A zone that WAS picked and whose track could not be built is a different case,
    and that one dims with the generator's reason -- the model was understood and the
    geometry refused, which is worth saying.
    """
    from .generators import generate_noncircular_surface
    rules = rule_ids_for("generate_noncircular_surface")
    base = dict(remedy="generate_noncircular_surface", target="non_circ",
                key="generate_noncircular_surface:non_circ", label="non-circ sheet",
                rule_ids=rules)
    if sd.get("non_circ"):
        return []                       # there is already a non-circular surface
    result = generate_noncircular_surface(sd, report=True)
    if not result["surface"]:
        if result["zone"] is None:
            # No zone was chosen: no zones at all, one zone, or two comparable ones.
            # None of those is this remedy's business; see the note above.
            return []
        return [RemedyProposal(available=False, reason=result["reason"], **base)]
    n = len(result["surface"])
    return [RemedyProposal(
        available=True,
        description=(f"Add a {n}-point non-circular surface to the non-circ sheet, "
                     f"tracking the '{result['zone'].name}' zone: "
                     f"{result['summary']}."),
        _apply=_apply_noncirc, **base)]


def _apply_noncirc(sd):
    """The same call the proposal made, so the button's label and the change agree."""
    from .generators import generate_noncircular_surface
    result = generate_noncircular_surface(sd, report=True)
    if not result["surface"]:
        raise RemedyDeclined(result["reason"] or
                             "no weak-zone surface could be derived from the geometry")
    out = _shallow(sd)
    out["non_circ"] = list(sd.get("non_circ") or []) + list(result["surface"])
    out["circular"] = False
    return out


# ---------------------------------------------------------------------------
# The registry, and the three entry points
# ---------------------------------------------------------------------------

#: ``name -> proposal builder``. A remedy named in
#: :data:`xslope.preflight.REMEDIES` but absent here is declared, not built:
#: asking for it raises rather than quietly doing nothing.
REMEDY_BUILDERS = {
    "reverse_polyline": _propose_reverse,
    "add_ponded_water_load": _propose_ponded,
    "switch_to_auto_water": _propose_switch_to_auto,
    "generate_starting_circles": _propose_circles,
    "generate_noncircular_surface": _propose_noncirc,
}


def implemented():
    """The names of the remedies that have a transformation behind them."""
    return tuple(sorted(REMEDY_BUILDERS))


def remedy_proposals(slope_data, remedies=None):
    """Every remedy this model could be offered, available or dimmed with a reason.

    One proposal per *target*: the ponded-water remedy is evaluated per stage and
    the reversal remedy per line, because a button is offered per object and an
    availability answer that averaged over several would be wrong for each of them.

    Parameters
    ----------
    slope_data : dict
        Nothing is mutated.
    remedies : iterable of str, optional
        Restrict to these remedy names. Unknown names raise.

    Returns
    -------
    list of RemedyProposal
    """
    names = list(remedies) if remedies is not None else list(REMEDY_BUILDERS)
    out = []
    for name in names:
        if name not in _pf.REMEDIES:
            raise ValueError(f"unknown remedy {name!r}; known remedies: "
                             f"{', '.join(sorted(_pf.REMEDIES))}")
        builder = REMEDY_BUILDERS.get(name)
        if builder is None:
            raise ValueError(
                f"the remedy {name!r} is declared but not implemented, so there is "
                f"nothing it can do: {_pf.REMEDIES[name]}")
        out.extend(builder(slope_data or {}))
    return out


def propose(slope_data, remedy, target=None):
    """One remedy's proposal for one target, computed without changing anything.

    ``target`` may be omitted when the remedy has exactly one target on this model.
    Raises :class:`RemedyDeclined` when there is nothing for the remedy to act on,
    since a proposal that described no change would be indistinguishable from one
    that had already been applied.
    """
    found = [p for p in remedy_proposals(slope_data, [remedy])
             if target is None or p.target == target or p.key == target]
    if not found:
        raise RemedyDeclined(
            f"the remedy {remedy!r} has nothing to act on in this model"
            + (f" (target {target!r})" if target else ""))
    if len(found) > 1:
        raise RemedyDeclined(
            f"the remedy {remedy!r} has {len(found)} targets on this model "
            f"({', '.join(p.target for p in found)}); name the one you mean")
    return found[0]


def remedy_capabilities(slope_data, remedies=None):
    """``{key: Capability}`` for every remedy this model could be offered.

    The interface query for remedy buttons, in the same shape and with the same
    reason strings as :func:`xslope.preflight.capabilities`: a button binds to a
    proposal's ``key``, shows it dimmed when ``available`` is False, and puts
    ``reason`` in the tooltip.
    """
    return {p.key: _pf.Capability(p.available, p.reason)
            for p in remedy_proposals(slope_data, remedies)}


def apply_remedy(slope_data, remedy, target=None):
    """Apply one remedy and report what it did.

    Returns ``(model, finding)``: a copy of the model with the change applied, and
    the INFO finding that records it. The input model is not touched.

    Raises :class:`RemedyDeclined`, carrying the proposal's own reason, when the
    remedy cannot fully succeed. Nothing is half-applied: either the whole change
    lands or the model comes back unchanged and the original finding stands.
    """
    proposal = propose(slope_data, remedy, target)
    if not proposal.available:
        raise RemedyDeclined(proposal.reason)
    model = proposal._apply(slope_data)
    return model, applied_finding(proposal)


def applied_finding(proposal):
    """The INFO that records an applied remedy, for the report it travels with."""
    rules = ", ".join(proposal.rule_ids)
    return _pf.Finding(
        rule_id=f"remedy.{proposal.remedy}",
        severity=_pf.INFO,
        message=(f"Remedy applied to this model, on request: {proposal.description} "
                 f"This was offered by "
                 f"{'the rule ' + rules if rules else 'preflight'}, and the model that "
                 f"ran is not the file as saved."),
        remedy=proposal.remedy)
