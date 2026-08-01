"""The model checks a Run dialog shows beside its Run button.

One widget, shared by every Run dialog, so the interface can never disagree with
the gate: :class:`PreflightPanel` calls :func:`xslope.preflight.preflight` on the
model the run would actually use and renders what comes back, in the registry's own
words. Nothing here writes user-facing copy — an error's text, a warning's text and
a dimmed control's tooltip are all the rule's own message, taken from the finding or
from the rule's own one-line summary.

The behaviour is the one the design settled on:

* **ERROR refuses.** The dialog's Run button is disabled while any error stands,
  and the error is on screen with the message the run would have failed with.
* **WARNING never blocks**, in Studio or anywhere else. Warnings are shown
  prominently beside the button so they inform the decision rather than annotate
  the regret.
* **INFO is available, not intrusive** — collapsed behind a disclosure line, and
  counted so it is visibly there.
* **A remedy is offered, never applied silently.** A finding that carries one gets
  a button; the button is capability-gated (dimmed with the reason when the remedy
  cannot fully succeed, never live-and-failing); pressing it shows exactly what the
  change will be *before* it happens, applies it through the document's undo stack
  on confirm, then re-runs the checks so the list reflects the model that would now
  run.

**A list, not a wall.** The findings are a one-line-per-entry list with a detail
area beneath it: selecting a line shows that finding's full text, and its remedy
buttons, underneath. A dialog with five findings is then five lines tall instead of
five paragraphs tall, which is what lets the run controls and the checks sit side by
side (:func:`two_pane`) instead of stacked.

**One rule, one line.** A rule that fires once per material fired four paragraphs of
near-identical text at the user — same explanation, same remedy, four times over. So
findings that share a rule id are rendered as ONE entry here: the line names the
count and the members, and the detail states the shared explanation once with the
per-finding specifics as a sub-list. This is a presentation rule and lives here, in
the panel: :mod:`xslope.preflight` still reports one finding per material, and every
scripted consumer still sees them individually.

The panel re-evaluates on demand (``refresh()``), which the dialogs call whenever a
control that feeds the selection changes — the method, the surface family, the
analysis type. That is what makes the findings track the run being configured
rather than the file as opened.
"""

from __future__ import annotations

import html
import re

from PySide6.QtCore import QSize, Qt, Signal
from PySide6.QtWidgets import (
    QAbstractItemView, QFrame, QGridLayout, QGroupBox, QHBoxLayout, QLabel,
    QLayout, QListWidget, QListWidgetItem, QMessageBox, QPushButton, QScrollArea,
    QSizePolicy, QStackedWidget, QStyle, QToolButton, QVBoxLayout, QWidget,
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

#: Which standard icon marks each severity in the list.
_SEVERITY_ICON = {
    "error": QStyle.SP_MessageBoxCritical,
    "warning": QStyle.SP_MessageBoxWarning,
    "info": QStyle.SP_MessageBoxInformation,
}

#: At most this many members are named on an aggregate's one line; the rest are
#: counted. The full set is always in the detail.
_MAX_NAMED = 5


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


def two_pane(dialog, controls, panel, buttons, checks_width=360, fixed=False):
    """Lay a Run dialog out as run controls | model checks, buttons beneath both.

    ``controls`` is the layout the dialog filled with its own widgets, built
    UNPARENTED (``QVBoxLayout()``, not ``QVBoxLayout(self)``) so it can be moved
    into the left column here. ``panel`` is the :class:`PreflightPanel` and
    ``buttons`` the dialog's button box; both are simply placed, never re-wired.

    The controls keep their natural height at the top of the left column and the
    checks fill the dialog's height beside them, so a model with five findings
    makes the dialog WIDER rather than taller. ``fixed`` reproduces the
    ``SetFixedSize`` behaviour of a dialog whose rows appear and disappear.
    """
    left = QWidget(dialog)
    controls.setContentsMargins(0, 0, 0, 0)
    left.setLayout(controls)

    grid = QGridLayout(dialog)
    grid.setContentsMargins(10, 10, 10, 10)
    grid.setHorizontalSpacing(12)
    grid.setVerticalSpacing(8)
    grid.addWidget(left, 0, 0, Qt.AlignTop)
    grid.addWidget(panel, 0, 1)
    grid.addWidget(buttons, 1, 0, 1, 2)
    grid.setColumnStretch(0, 0)
    grid.setColumnStretch(1, 1)
    grid.setRowStretch(0, 1)

    panel.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Expanding)
    panel.setMinimumWidth(checks_width)
    if fixed:
        grid.setSizeConstraint(QLayout.SetFixedSize)
    return grid


# ---------------------------------------------------------------------------
# Grouping and wording
#
# Everything below reads the findings' own text; none of it writes any. A group's
# line is the rule's registered summary (or, for a lone finding, its own first
# sentence) plus a count and the members it names; a group's detail is every
# message, with the sentences they all share stated once instead of per member.
# ---------------------------------------------------------------------------

_SENTENCE = re.compile(r"(?<=[.!?])\s+(?=[A-Z(\"'])")

#: What identifies one member once the text its group shares has been taken off the
#: front: a number or a name, with the parenthesized name that may follow it --
#: ``1 ('Core')`` out of ``Material 1 ('Core') has no tensile strength: ...``.
_MEMBER = re.compile(r"[^\s,;:]+(?:\s*\([^)]*\))?")

_SUMMARIES = {}


def rule_summary(rule_id):
    """The registry's one-line summary of a rule, or ``""`` if it has none."""
    if not _SUMMARIES:
        try:
            from xslope.preflight import rules
            _SUMMARIES.update({r.id: r.summary for r in rules()})
        except Exception:                          # pragma: no cover - import guard
            pass
    return _SUMMARIES.get(rule_id, "")


def _sentences(text):
    return [s.strip() for s in _SENTENCE.split(str(text).strip()) if s.strip()]


def _first_sentence(text):
    got = _sentences(text)
    return got[0] if got else str(text).strip()


def _plural(word):
    if word.endswith("s"):
        return word
    if word.endswith("y") and word[-2:-1] not in "aeiou":
        return word[:-1] + "ies"
    return word + "s"


def _members(messages):
    """``(noun, [member, ...])`` for a group's messages.

    A rule that fires per material opens every message the same way and then names
    which one: ``Material 1 ('Sand') has no tensile strength: ...``. So the shared
    opening is what the members have in common, the last word of it is what they ARE
    (``Material`` -> ``4 materials``, lowercased because it is being counted mid-line
    rather than opening a sentence), and each member is the identifier that follows
    it (``1 ('Sand')``).

    Where that reading does not hold -- the messages share nothing, or what follows
    is the same word every time and so names nothing -- the count stands alone and
    the members are left to the detail.
    """
    shared = messages[0]
    for m in messages[1:]:
        n = 0
        while n < min(len(shared), len(m)) and shared[n] == m[n]:
            n += 1
        shared = shared[:n]
    cut = shared.rfind(" ")
    shared = shared[:cut + 1] if cut > 0 else ""
    words = shared.strip().rstrip(",").split()
    noun = words[-1] if words and words[-1].isalpha() and len(words[-1]) > 2 else None
    members = []
    for m in messages:
        hit = _MEMBER.match(m[len(shared):].strip())
        members.append(hit.group(0).rstrip(".,;:") if hit else "")
    if not (noun and all(members) and len(set(members)) == len(members)):
        return "findings", []
    return _plural(noun.lower()), members


def _shared_affix(messages):
    """``(head, middles, tail)`` -- the sentences every message opens and closes with.

    The middles are what is left of each message once those are lifted out: the
    per-finding specifics. When lifting them would leave any message empty (the
    messages differ only in their shared text) nothing is lifted and the messages
    stand whole.
    """
    parts = [_sentences(m) for m in messages]
    n_head = 0
    while all(len(p) > n_head + 1 for p in parts) and \
            len({p[n_head] for p in parts}) == 1:
        n_head += 1
    n_tail = 0
    while all(len(p) > n_head + n_tail + 1 for p in parts) and \
            len({p[len(p) - 1 - n_tail] for p in parts}) == 1:
        n_tail += 1
    middles = [" ".join(p[n_head:len(p) - n_tail]) for p in parts]
    if not all(middles):
        return "", list(messages), ""
    head = " ".join(parts[0][:n_head])
    tail = " ".join(parts[0][len(parts[0]) - n_tail:]) if n_tail else ""
    return head, middles, tail


class _LineLabel(QLabel):
    """A label that draws its one line elided, and still reports the whole of it.

    The list is as wide as the dialog's second column, not as wide as the longest
    finding; a line too long for it ends in an ellipsis rather than being cut mid
    word. :meth:`text` is untouched, so anything reading the panel — the tooltip,
    the tests — sees the line the panel meant.
    """

    def paintEvent(self, event):
        from PySide6.QtGui import QFontMetrics, QPainter
        painter = QPainter(self)
        elided = QFontMetrics(self.font()).elidedText(
            self.text(), Qt.ElideRight, self.width())
        self.style().drawItemText(painter, self.rect(), int(self.alignment()),
                                  self.palette(), self.isEnabled(), elided,
                                  self.foregroundRole())


class _DetailArea(QScrollArea):
    """A scroll area that asks for a few lines, not for all of its content.

    The detail of one finding is a paragraph; asking for the whole paragraph would
    make every dialog as tall as its longest message, which is the shape this
    layout exists to get away from. It asks for the few lines :meth:`~PreflightPanel
    ._fit` gives it and scrolls past the rest, and grows when the dialog does.
    """

    def sizeHint(self):
        hint = super().sizeHint()
        return QSize(hint.width(), self.minimumHeight() or hint.height())


class FindingGroup:
    """The findings of one rule, rendered as one entry.

    A group of one is an ordinary finding; a group of several is the aggregate the
    panel shows in its place. ``findings`` keeps every member, in the order the
    registry reported them, so nothing is lost to the grouping.
    """

    def __init__(self, findings):
        self.findings = list(findings)
        self.rule_id = self.findings[0].rule_id
        self.severity = self.findings[0].severity

    def __len__(self):
        return len(self.findings)

    @property
    def messages(self):
        return [f.message for f in self.findings]

    def line(self):
        """The one line this group takes in the list."""
        msgs = self.messages
        if len(msgs) == 1:
            return _first_sentence(msgs[0])
        head = (rule_summary(self.rule_id) or _first_sentence(msgs[0])).strip(" .")
        noun, members = _members(msgs)
        if not members:
            return f"{head} — {len(msgs)} {noun}"
        shown = members[:_MAX_NAMED]
        more = len(members) - len(shown)
        named = ", ".join(shown) + (f", and {more} more" if more else "")
        return f"{head} — {len(msgs)} {noun}: {named}"

    def detail(self):
        """The full text of this group, as rich text: the explanation once."""
        msgs = self.messages
        if len(msgs) == 1:
            return f"<p>{html.escape(msgs[0])}</p>"
        head, middles, tail = _shared_affix(msgs)
        summary = rule_summary(self.rule_id)
        out = []
        if summary:
            out.append(f"<p><b>{html.escape(summary)}</b></p>")
        if head:
            out.append(f"<p>{html.escape(head)}</p>")
        out.append("<ul style='margin-left:0px;-qt-list-indent:1;'>"
                   + "".join(f"<li>{html.escape(m)}</li>" for m in middles)
                   + "</ul>")
        if tail:
            out.append(f"<p>{html.escape(tail)}</p>")
        return "".join(out)


def group_findings(findings):
    """Findings grouped by rule id, in first-appearance order.

    The grouping is by rule, and severity follows the rule, so an aggregate's
    severity is its members' severity. Callers pass one severity bucket at a time,
    which keeps a rule that overrides its own severity per finding from merging
    an error with a warning.
    """
    order, by_id = [], {}
    for f in findings:
        if f.rule_id not in by_id:
            by_id[f.rule_id] = []
            order.append(f.rule_id)
        by_id[f.rule_id].append(f)
    return [FindingGroup(by_id[rid]) for rid in order]


class PreflightPanel(QGroupBox):
    """Findings for one run: a list of lines, with the selected one's detail beneath.

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
        self._infos_open = False
        self._groups = []
        self._first_info_row = 0

        outer = QVBoxLayout(self)
        outer.setContentsMargins(8, 6, 8, 6)
        outer.setSpacing(4)

        # The list: one line per (aggregated) finding, dense on every platform —
        # the row height is the font's, plus a couple of pixels, not the style's
        # idea of a comfortable list.
        self._list = QListWidget()
        self._list.setObjectName("preflight_list")
        self._list.setSelectionMode(QAbstractItemView.SingleSelection)
        self._list.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        self._list.setTextElideMode(Qt.ElideRight)
        self._list.setSpacing(0)
        self._list.setUniformItemSizes(False)
        self._list.setFrameShape(QFrame.StyledPanel)
        # The selected line is marked with the palette's own light band rather than
        # the highlight colour: the lines are coloured by severity, and an error's
        # red on a full-strength selection blue is the one pairing that cannot be
        # read. palette() keeps this right in the dark theme too.
        self._list.setStyleSheet(
            "QListWidget::item:selected { background: palette(midlight); }")
        self._list.currentRowChanged.connect(self._show_selected)
        outer.addWidget(self._list)

        # The info disclosure: one line saying how many notes there are, which
        # shows and hides them in the list above.
        self._info_toggle = QToolButton()
        self._info_toggle.setObjectName("preflight_infos")
        self._info_toggle.setCheckable(True)
        self._info_toggle.setChecked(False)
        self._info_toggle.setToolButtonStyle(Qt.ToolButtonTextBesideIcon)
        self._info_toggle.setArrowType(Qt.RightArrow)
        self._info_toggle.setAutoRaise(True)
        self._info_toggle.toggled.connect(self._toggle_infos)
        outer.addWidget(self._info_toggle)

        # The detail: the selected entry's full text and its remedy buttons. One
        # page per entry, so every remedy button a model offers exists (and can be
        # inspected) whether or not its line is the one selected.
        self._detail = QStackedWidget()
        self._detail.setObjectName("preflight_detail")
        self._detail_scroll = _DetailArea()
        self._detail_scroll.setWidgetResizable(True)
        self._detail_scroll.setFrameShape(QFrame.NoFrame)
        self._detail_scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        self._detail_scroll.setWidget(self._detail)
        outer.addWidget(self._detail_scroll, 1)

        self.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Preferred)
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

    def groups(self):
        """The entries in the list, row by row, most urgent first.

        One per list row, so ``groups()[row]`` is what that row is about; the
        all-clear line a model with nothing to report gets is a ``None``.
        """
        return list(self._groups)

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
        # Row widgets are un-parented as they go, not merely scheduled for deletion:
        # a caller that looks the findings up by object name between two refreshes
        # (the tests do, and so does anything reading the panel without an event
        # loop) must never see the previous evaluation's rows still hanging on.
        for row in range(self._list.count()):
            item = self._list.item(row)
            widget = self._list.itemWidget(item)
            if widget is not None:
                self._list.removeItemWidget(item)
                widget.setParent(None)
                widget.deleteLater()
        self._list.clear()
        while self._detail.count():
            page = self._detail.widget(0)
            self._detail.removeWidget(page)
            page.setParent(None)
            page.deleteLater()

    def _rebuild(self):
        self._list.blockSignals(True)
        self._clear()
        self._rendered_remedies = set()
        report = self._report
        errors = list(report.errors)
        warnings = list(report.warnings)
        notes = [f for f in report.infos if f.rule_id in self._notes]
        infos = [f for f in report.infos if f.rule_id not in self._notes]

        shown = (group_findings(errors) + group_findings(warnings)
                 + group_findings(notes))
        hidden = group_findings(infos)
        for group in shown:
            self._add_entry(group, hidden=False)
        rows = list(shown)
        if not shown:
            # Nothing that is a problem: say so on the one line, and leave the
            # notes (if any) where they belong, behind their disclosure.
            self._add_plain("No problems found for this run.", "preflight_ok")
            rows.append(None)
        self._first_info_row = len(rows)
        for group in hidden:
            self._add_entry(group, hidden=not self._infos_open)
        self._groups = rows + hidden

        self._info_toggle.setVisible(bool(hidden))
        if hidden:
            self._info_toggle.setText(
                f"{len(infos)} note{'' if len(infos) == 1 else 's'}")
        self._list.blockSignals(False)
        self._select_default()
        self._show_selected(self._list.currentRow())

        self.setTitle(self._title_for(shown))
        self._fit()

    def _add_entry(self, group, hidden):
        """One list line plus its detail page."""
        row = QWidget()
        lay = QHBoxLayout(row)
        lay.setContentsMargins(4, 1, 4, 1)
        lay.setSpacing(6)
        icon = QLabel()
        which = _SEVERITY_ICON.get(group.severity)
        if which is not None:
            side = max(12, self.fontMetrics().height() - 2)
            icon.setPixmap(self.style().standardIcon(which).pixmap(side, side))
        lay.addWidget(icon, 0, Qt.AlignVCenter)

        text = group.line()
        label = _LineLabel(text)
        label.setObjectName(f"finding_{group.rule_id}")
        label.setWordWrap(False)
        label.setToolTip(text)
        # Ignored, so a long line never widens the dialog: the list elides it and
        # the whole of it is in the detail below.
        label.setSizePolicy(QSizePolicy.Ignored, QSizePolicy.Preferred)
        prefix, color = _SEVERITY_STYLE.get(group.severity, ("", "#000000"))
        if group.severity != "info":
            label.setStyleSheet(f"color:{color};")
        lay.addWidget(label, 1)

        item = QListWidgetItem()
        item.setSizeHint(self._row_size(row))
        item.setToolTip(f"{prefix}: {text}")
        item.setData(Qt.UserRole, self._detail.count())
        self._list.addItem(item)
        self._list.setItemWidget(item, row)
        item.setHidden(bool(hidden))

        self._detail.addWidget(self._detail_page(group))

    def _row_size(self, row):
        """A dense row: the widget's own height, floored at one line and a little.

        Explicit, because the platform styles disagree about what a comfortable
        list row is by more than a title bar's worth over five findings.
        """
        return QSize(0, max(row.sizeHint().height(), self.fontMetrics().height() + 6))

    def _add_plain(self, text, name):
        """A line that is not a finding (the all-clear), with a matching page."""
        row = QWidget()
        lay = QHBoxLayout(row)
        lay.setContentsMargins(4, 1, 4, 1)
        label = QLabel(text)
        label.setObjectName(name)
        lay.addWidget(label, 1)
        item = QListWidgetItem()
        item.setSizeHint(self._row_size(row))
        item.setData(Qt.UserRole, self._detail.count())
        item.setFlags(Qt.ItemIsEnabled)        # nothing to open: not a selection
        self._list.addItem(item)
        self._list.setItemWidget(item, row)

        self._detail.addWidget(QWidget())      # its line says all there is to say

    def _detail_page(self, group):
        """The full text of one entry, with its remedy buttons under it."""
        page = QWidget()
        lay = QVBoxLayout(page)
        lay.setContentsMargins(2, 2, 2, 2)
        lay.setSpacing(4)
        prefix, color = _SEVERITY_STYLE.get(group.severity, ("", "#000000"))
        body = QLabel(f"<b style='color:{color}'>{prefix}:</b> {group.detail()}")
        body.setObjectName(f"detail_{group.rule_id}")
        body.setWordWrap(True)
        body.setTextFormat(Qt.RichText)
        body.setTextInteractionFlags(Qt.TextSelectableByMouse)
        body.setAlignment(Qt.AlignTop)
        lay.addWidget(body)
        for button in self._remedy_buttons(group.findings[0]):
            strip = QHBoxLayout()
            strip.setContentsMargins(0, 0, 0, 0)
            strip.addWidget(button)
            strip.addStretch(1)
            holder = QWidget()
            holder.setLayout(strip)
            lay.addWidget(holder)
        lay.addStretch(1)
        return page

    def _toggle_infos(self, on):
        """Show or hide the note lines; the arrow follows, and so does selection."""
        self._infos_open = bool(on)
        self._info_toggle.setArrowType(Qt.DownArrow if on else Qt.RightArrow)
        for row in range(self._first_info_row, self._list.count()):
            item = self._list.item(row)
            if item is not None:
                item.setHidden(not on)
        if not on and self._list.currentRow() >= self._first_info_row:
            self._select_default()
        self._fit()

    def _select_default(self):
        """Open on the first ERROR, else the first line — never on nothing."""
        for row in range(self._list.count()):
            item = self._list.item(row)
            if item is None or item.isHidden():
                continue
            group = self._groups[row] if row < len(self._groups) else None
            if group is not None and group.severity == "error":
                self._list.setCurrentRow(row)
                return
        for row in range(self._list.count()):
            if not self._list.item(row).isHidden():
                self._list.setCurrentRow(row)
                return

    def _show_selected(self, row):
        """Show the selected line's detail.

        The detail area stays in place whatever is selected -- it is the layout's
        one growing element, and taking it away would leave the list and the notes
        line drifting apart in the space neither of them wants. What an all-clear
        line shows there is nothing: the sentence is already on the line, and
        printing it twice is how an empty model ends up looking busy.
        """
        item = self._list.item(row) if row >= 0 else None
        if item is None:
            return
        page = item.data(Qt.UserRole)
        if page is not None and 0 <= int(page) < self._detail.count():
            self._detail.setCurrentIndex(int(page))

    def _fit(self):
        """Size the list to the lines it is showing, up to a share of the screen."""
        visible = [r for r in range(self._list.count())
                   if not self._list.item(r).isHidden()]
        if not visible:
            return
        row_h = max(self._list.sizeHintForRow(visible[0]),
                    self.fontMetrics().height() + 6)
        rows = min(max(len(visible), 1), self._max_rows(row_h))
        frame = 2 * self._list.frameWidth() + 2
        self._list.setMinimumHeight(rows * row_h + frame)
        self._list.setMaximumHeight(rows * row_h + frame)
        # Room for a finding's full text under the list — but a model with nothing
        # to report is not given six empty lines to say so.
        lines = 6 if any(g is not None for g in self._groups) else 2
        self._detail_scroll.setMinimumHeight(lines * (self.fontMetrics().height() + 2))

    @staticmethod
    def _max_rows(row_h):
        """How many lines the list may show before it scrolls: a share of the screen."""
        from PySide6.QtGui import QGuiApplication
        screen = QGuiApplication.primaryScreen()
        available = screen.availableGeometry().height() if screen else 900
        return max(4, int(available * 0.25 / max(row_h, 1)))

    #: What the header calls each severity, singular and plural.
    _TITLE_WORDS = (("error", "errors"), ("warning", "warnings"), ("note", "notes"))

    @staticmethod
    def _title_for(shown):
        """The header line: what the LIST is showing, plus how to open a line.

        Counted over the ENTRIES, not over the raw findings, because the entries are
        what the user can see: a rule that fired on four materials is one line here
        and says "— 4 materials" on that line itself. Counting the findings instead
        made the header claim five warnings above a list of two, and the multiplicity
        it was counting was already stated one line down.

        The collapsed infos are deliberately absent: they are not in the list, and
        their own disclosure counts them.
        """
        counts = {"error": 0, "warning": 0, "info": 0}
        for group in shown:
            if group is not None and group.severity in counts:
                counts[group.severity] += 1
        bits = []
        for severity, (one, many) in zip(("error", "warning", "info"),
                                         PreflightPanel._TITLE_WORDS):
            n = counts[severity]
            if n:
                bits.append(f"{n} {one if n == 1 else many}")
        title = "Model checks" + (f" — {' · '.join(bits)}" if bits else "")
        # The list-then-detail interaction is not self-evident: a user who does not
        # know a line opens has a one-line summary and no way to reach the rest of
        # the message. It is only offered where there is something to select -- the
        # all-clear line is not selectable and says all it has to say.
        return title + ("   (select a line for details)" if bits else "")

    # -- remedies ---------------------------------------------------------
    def _remedy_buttons(self, finding):
        """The buttons for one entry's remedies, capability-gated.

        A remedy is offered per *target* — one line, one stage — so a rule that
        fires several times shares one set of buttons rather than repeating them:
        a proposal already rendered above is not rendered again. With the findings
        of a rule aggregated into one entry, the entry asks for its rule's remedies
        once and gets every target's proposal, which is the same set the repeated
        rows used to share.

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
