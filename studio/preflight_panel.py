"""The model checks a Run dialog shows before its Run button.

One widget, shared by every Run dialog, so the interface can never disagree with
the gate: :class:`PreflightPanel` calls :func:`xslope.preflight.preflight` on the
model the run would actually use and renders what comes back, in the registry's own
words. Nothing here writes user-facing copy — an error's text, a warning's text and
a dimmed control's tooltip are all the rule's own message, taken from the finding.

The behaviour is the one the design settled on:

* **ERROR refuses.** The dialog's Run button is disabled while any error stands,
  and the error is on screen with the message the run would have failed with.
* **WARNING never blocks**, in Studio or anywhere else. Warnings are shown
  prominently *above* the button so they inform the decision rather than annotate
  the regret.
* **INFO is available, not intrusive** — collapsed behind a disclosure line, and
  counted so it is visibly there.
* **A remedy is offered, never applied silently.** A finding that carries one gets
  a button; the button is capability-gated (dimmed with the reason when the remedy
  cannot fully succeed, never live-and-failing); pressing it shows exactly what the
  change will be *before* it happens, applies it through the document's undo stack
  on confirm, then re-runs the checks so the list reflects the model that would now
  run.

The panel re-evaluates on demand (``refresh()``), which the dialogs call whenever a
control that feeds the selection changes — the method, the surface family, the
analysis type. That is what makes the findings track the run being configured
rather than the file as opened.
"""

from __future__ import annotations

from PySide6.QtCore import Qt, Signal
from PySide6.QtWidgets import (
    QFrame, QGroupBox, QHBoxLayout, QLabel, QMessageBox, QPushButton,
    QScrollArea, QSizePolicy, QToolButton, QVBoxLayout, QWidget,
)

#: The two seismic-direction infos, which the Run dialogs show as prominent notes
#: rather than behind the info disclosure: one template cell means a magnitude to
#: the limit-equilibrium engine and a vector to the finite element engine, and a
#: user carrying a habit from one to the other has no other way to learn it.
SEISMIC_NOTE_LEM = "main.seismic_direction_lem"
SEISMIC_NOTE_FEM = "main.seismic_direction_fem"

#: Prefix and color per severity. Colors are chosen to read on both the light and
#: dark palettes Studio ships; the prefix is what carries the meaning if they do not.
_SEVERITY_STYLE = {
    "error": ("Error", "#b3261e"),
    "warning": ("Warning", "#8a5a00"),
    "info": ("Note", "#5a5a5a"),
}


def apply_capabilities(combo, group, keep_current=False):
    """Dim every item of ``combo`` its model cannot run, with the gate's own reason.

    ``group`` is one group of a :func:`xslope.preflight.capabilities` result --
    ``{option: Capability}``. An unavailable option's item is disabled and carries
    the rule's message as its tooltip, which is the same sentence the run would have
    been refused with. When the current selection is one of the dimmed ones the
    combo moves to the first available option, so a dialog never opens on a choice
    it would refuse; pass ``keep_current=True`` to leave it alone.

    Returns the reason the *current* option is unavailable, or ``""``.
    """
    model = combo.model()
    for i in range(combo.count()):
        cap = group.get(combo.itemData(i))
        item = model.item(i)
        if cap is None or item is None:
            continue
        item.setEnabled(bool(cap.available))
        item.setToolTip(cap.reason or "")
    cap = group.get(combo.currentData())
    if cap is not None and not cap.available and not keep_current:
        for i in range(combo.count()):
            other = group.get(combo.itemData(i))
            if other is None or other.available:
                combo.setCurrentIndex(i)
                return ""
    return "" if cap is None or cap.available else cap.reason


class PreflightPanel(QGroupBox):
    """Findings for one run, rendered above a Run button.

    Parameters
    ----------
    analysis : str or callable
        The analysis type to check against, or a callable returning it — a
        callable because the Run LEM dialog's rapid-drawdown checkbox and the Run
        Seepage dialog's steady/transient selector change which analysis the
        registry should be asked about.
    slope_data : dict
        The model. Held by reference: a remedy replaces its contents in place, so
        the dialog and the panel never drift apart.
    document : studio.document.ProjectDocument, optional
        When present a remedy is applied through ``begin_edit`` / ``commit_edit``,
        so it lands on the normal undo stack and is saved like any other edit.
        Without one (the test harness) the change is made to ``slope_data``
        directly.
    selection_fn : callable, optional
        Returns the ``selection`` dict for :func:`~xslope.preflight.preflight` —
        the method, the surface family, whether this is a search. Called on every
        refresh, so it reads live widget state.
    notes : iterable of str, optional
        Rule ids to lift OUT of the collapsed info list and show as prominent
        notes. The seismic direction infos are the reason this exists: they are
        education rather than defects, and the design calls for them beside the
        button rather than behind a disclosure.
    """

    #: Emitted after every evaluation, including one triggered by a remedy.
    changed = Signal()

    def __init__(self, analysis, slope_data=None, document=None, selection_fn=None,
                 notes=(), parent=None, title="Model checks"):
        super().__init__(title, parent)
        self._analysis = analysis
        self._sd = slope_data if slope_data is not None else {}
        self._doc = document
        self._selection_fn = selection_fn
        self._notes = tuple(notes or ())
        self._report = None
        self._rendered_remedies = set()

        outer = QVBoxLayout(self)
        outer.setContentsMargins(8, 6, 8, 6)
        outer.setSpacing(4)

        self._body = QWidget()
        self._rows = QVBoxLayout(self._body)
        self._rows.setContentsMargins(0, 0, 0, 0)
        self._rows.setSpacing(6)

        # Long lists scroll rather than growing the dialog past the screen; short
        # ones (the common case) size to their content, so the box does not carry
        # empty space when a model has one warning. The ceiling is a share of the
        # screen rather than a pixel count, so the same panel behaves on a laptop
        # and on a large display.
        self._scroll = QScrollArea()
        self._scroll.setWidgetResizable(True)
        self._scroll.setFrameShape(QFrame.NoFrame)
        self._scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        self._scroll.setWidget(self._body)
        outer.addWidget(self._scroll)

        self.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Maximum)
        self.refresh()

    # -- state ------------------------------------------------------------
    @property
    def model(self):
        """The model the findings describe (the same dict a remedy updates)."""
        return self._sd

    @property
    def report(self):
        """The last :class:`~xslope.preflight.PreflightReport`."""
        return self._report

    @property
    def blocked(self):
        """True when an ERROR stands, so the run must be refused."""
        return bool(self._report is not None and self._report.errors)

    def block_reason(self):
        """The first error's message — what a disabled Run button explains."""
        errs = self._report.errors if self._report is not None else []
        return errs[0].message if errs else ""

    def analysis_name(self):
        a = self._analysis
        return a() if callable(a) else a

    # -- evaluation -------------------------------------------------------
    def refresh(self):
        """Re-run preflight and rebuild the list. Emits :attr:`changed`."""
        from xslope import preflight as pf

        selection = dict(self._selection_fn() if self._selection_fn else {})
        try:
            self._report = pf.preflight(self._sd, self.analysis_name(), selection)
        except Exception as exc:                    # an unknown analysis, say
            self._report = pf.PreflightReport(
                analysis=str(self.analysis_name()),
                findings=[pf.Finding("preflight.unavailable", pf.INFO,
                                     f"The model checks could not be run "
                                     f"({type(exc).__name__}: {exc}).")])
        self._rebuild()
        self.changed.emit()

    # -- rendering --------------------------------------------------------
    def _clear(self):
        while self._rows.count():
            item = self._rows.takeAt(0)
            w = item.widget()
            if w is not None:
                w.setParent(None)
                w.deleteLater()

    def _rebuild(self):
        self._clear()
        self._rendered_remedies = set()
        report = self._report
        errors = list(report.errors)
        warnings = list(report.warnings)
        notes = [f for f in report.infos if f.rule_id in self._notes]
        infos = [f for f in report.infos if f.rule_id not in self._notes]

        for finding in errors + warnings + notes:
            self._rows.addWidget(self._finding_row(finding))

        if infos:
            self._rows.addWidget(self._info_disclosure(infos))

        if not (errors or warnings or notes or infos):
            ok = QLabel("No problems found for this run.")
            ok.setWordWrap(True)
            self._rows.addWidget(ok)
        self._rows.addStretch(1)

        self.setTitle(self._title_for(errors, warnings))
        self._fit()

    def _fit(self):
        """Size the list to its content, up to a share of the screen.

        The body is inside a resizable scroll area, so its word-wrapped labels only
        know their height once they have been given the viewport's width -- hence
        the explicit activate() before the hint is read.
        """
        self._body.layout().activate()
        self._body.adjustSize()
        wanted = self._body.sizeHint().height()
        self._scroll.setMaximumHeight(min(wanted, self._ceiling()))

    @staticmethod
    def _ceiling():
        """How tall the findings list may grow: about a third of the screen."""
        from PySide6.QtGui import QGuiApplication
        screen = QGuiApplication.primaryScreen()
        available = screen.availableGeometry().height() if screen else 900
        return max(180, int(available * 0.33))

    @staticmethod
    def _title_for(errors, warnings):
        bits = []
        if errors:
            bits.append(f"{len(errors)} error{'' if len(errors) == 1 else 's'}")
        if warnings:
            bits.append(f"{len(warnings)} warning{'' if len(warnings) == 1 else 's'}")
        return "Model checks" + (f" — {', '.join(bits)}" if bits else "")

    def _finding_row(self, finding):
        row = QWidget()
        lay = QVBoxLayout(row)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.setSpacing(2)
        prefix, color = _SEVERITY_STYLE.get(finding.severity, ("", "#000000"))
        label = QLabel(f"<b style='color:{color}'>{prefix}:</b> {finding.message}")
        label.setWordWrap(True)
        label.setTextFormat(Qt.RichText)
        label.setObjectName(f"finding_{finding.rule_id}")
        lay.addWidget(label)
        for button in self._remedy_buttons(finding):
            strip = QHBoxLayout()
            strip.setContentsMargins(0, 0, 0, 0)
            strip.addWidget(button)
            strip.addStretch(1)
            holder = QWidget()
            holder.setLayout(strip)
            lay.addWidget(holder)
        return row

    def _info_disclosure(self, infos):
        """The infos, behind one line that says how many there are."""
        box = QWidget()
        lay = QVBoxLayout(box)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.setSpacing(2)
        toggle = QToolButton()
        toggle.setText(f"{len(infos)} note{'' if len(infos) == 1 else 's'}")
        toggle.setCheckable(True)
        toggle.setChecked(False)
        toggle.setToolButtonStyle(Qt.ToolButtonTextBesideIcon)
        toggle.setArrowType(Qt.RightArrow)
        toggle.setAutoRaise(True)
        lay.addWidget(toggle)
        body = QWidget()
        blay = QVBoxLayout(body)
        blay.setContentsMargins(12, 0, 0, 0)
        blay.setSpacing(4)
        for finding in infos:
            blay.addWidget(self._finding_row(finding))
        body.setVisible(False)
        lay.addWidget(body)

        def _toggled(on):
            body.setVisible(on)
            toggle.setArrowType(Qt.DownArrow if on else Qt.RightArrow)

        toggle.toggled.connect(_toggled)
        box.setObjectName("preflight_infos")
        return box

    # -- remedies ---------------------------------------------------------
    def _remedy_buttons(self, finding):
        """The buttons for one finding's remedies, capability-gated.

        A remedy is offered per *target* — one line, one stage — so a rule that
        fires several times shares one set of buttons rather than repeating them:
        a proposal already rendered above is not rendered again.

        A fault can have more than one right repair, and the rule does not choose
        between them: an empty surface sheet takes either a starting set of circles
        or a surface tracked along a weak zone, and which is right depends on what
        controls the mechanism. So every remedy the finding carries gets its own
        button, in the rule's own order with the primary first. A remedy whose
        conditions are not met simply produces no proposal and no button — a
        question the panel cannot answer is not asked as a dimmed control.
        """
        names = list(finding.remedies) or ([finding.remedy] if finding.remedy else [])
        if not names:
            return []
        from xslope import remedies as rem
        proposals = []
        for name in names:
            try:
                proposals += rem.remedy_proposals(self._sd, [name])
            except Exception:
                continue   # declared but not built: it offers nothing, the rest stand
        out = []
        for proposal in proposals:
            if proposal.key in self._rendered_remedies:
                continue
            self._rendered_remedies.add(proposal.key)
            btn = QPushButton(self._button_text(proposal))
            btn.setEnabled(bool(proposal.available))
            btn.setToolTip(proposal.reason if not proposal.available
                           else proposal.description)
            btn.setObjectName(f"remedy_{proposal.key}")
            btn.clicked.connect(lambda _=False, key=proposal.key,
                                name=proposal.remedy, target=proposal.target:
                                self._apply_remedy(name, target, key))
            out.append(btn)
        return out

    @staticmethod
    def _button_text(proposal):
        from xslope.preflight import REMEDIES
        text = {
            "reverse_polyline": "Reverse it…",
            "add_ponded_water_load": "Add the water load…",
            "switch_to_auto_water": "Switch to automatic water loads…",
            "generate_starting_circles": "Generate starting circles…",
            "generate_noncircular_surface": "Generate a surface along the weak zone…",
            "set_seismic_zero": "Set k to 0…",
        }.get(proposal.remedy)
        if text is None:                            # a remedy built after this list
            text = REMEDIES.get(proposal.remedy, proposal.remedy) + "…"
        return text

    def _apply_remedy(self, remedy, target, key):
        """Show what the remedy would do, apply on confirm, then re-check."""
        from xslope import remedies as rem

        try:
            proposal = rem.propose(self._sd, remedy, target)
        except Exception as exc:
            QMessageBox.information(self, "Nothing to change", str(exc))
            self.refresh()
            return
        if not proposal.available:
            QMessageBox.information(self, "Cannot be applied automatically",
                                    proposal.reason)
            self.refresh()
            return
        box = QMessageBox(self)
        box.setIcon(QMessageBox.Question)
        box.setWindowTitle("Apply this fix?")
        box.setText(proposal.label)
        box.setInformativeText(proposal.description)
        box.setStandardButtons(QMessageBox.Apply | QMessageBox.Cancel)
        box.setDefaultButton(QMessageBox.Apply)
        if box.exec() != QMessageBox.Apply:
            return
        self.apply_proposal(proposal)

    def apply_proposal(self, proposal):
        """Apply a proposal the user has confirmed, then re-check.

        Separate from the confirmation so the change itself can be exercised
        without a modal in the way -- the interaction is the dialog's, the
        transformation is the engine's.
        """
        from xslope import remedies as rem

        try:
            model, done = rem.apply_remedy(self._sd, proposal.remedy, proposal.target)
        except Exception as exc:
            QMessageBox.warning(self, "The fix was not applied", str(exc))
            self.refresh()
            return
        self._commit(model, f"Preflight fix: {proposal.remedy}")
        print(done.message)                    # the record travels to the Log pane
        self.refresh()

    def _commit(self, model, label):
        """Land a remedied model on the document, keeping the dict's identity.

        Every dialog and editor holds ``slope_data`` by reference, so the contents
        are replaced rather than the object: a remedy applied from the Run dialog is
        visible to the canvas and the editors without anything being re-wired.
        """
        if self._doc is not None:
            self._doc.begin_edit(label)
            target = self._doc.slope_data
        else:
            target = self._sd
        target.clear()
        target.update(model)
        if self._doc is not None:
            self._doc.commit_edit()
        self._sd = target
