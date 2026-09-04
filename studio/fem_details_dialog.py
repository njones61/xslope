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

"""Solution details for the FEM's reinforcement lines and piles.

The FEM results view colours whole members by how hard they are working, which
answers "which one should I look at" but not "what is it doing along its
length". This dialog answers the second question: pick a member on the left,
read its profiles on the right.

Under the list is a map of the section with the selected member picked out — a
profile is read along a member, and a name in a list is not a place on a slope.
It is the same drawing the report prints above its member details
(:func:`xslope.plot_fem_details.plot_member_map`), so the panel and the page
show one picture.

It carries the results view's own Field state switch, so the profiles can be
read on the at-failure mechanism or on the last converged solution, and the two
views can be set to the same instant of the analysis.

It is non-modal, so it can stay open beside the results view while the display
options change; it holds the solution bundle it was opened with and does not
re-solve. Everything it draws comes from :mod:`xslope.fem_details`, so a fresh
solve and a solution reloaded from its sidecars produce the same figures.
"""

import os

from PySide6.QtCore import QTimer, Qt, QSize
from PySide6.QtGui import QColor, QIcon, QPainter, QPixmap
from PySide6.QtWidgets import (QCheckBox, QComboBox, QDialog, QFileDialog,
                               QHBoxLayout, QLabel, QListWidget,
                               QListWidgetItem, QMessageBox, QPushButton,
                               QSplitter, QVBoxLayout, QWidget)

from .canvas import MplCanvas, BASE_DPI

# Badge colours, keyed by the band names fem_details assigns.
_BADGE_RGB = {
    "green": "#2e9e4f",
    "amber": "#e08a1e",
    "red": "#c0392b",
    "none": "#9aa0a6",
}


def _badge_icon(band, px=12):
    """A filled dot in the band's colour — the utilization badge on a list row."""
    pm = QPixmap(px, px)
    pm.fill(Qt.transparent)
    p = QPainter(pm)
    p.setRenderHint(QPainter.Antialiasing)
    p.setBrush(QColor(_BADGE_RGB.get(band, _BADGE_RGB["none"])))
    p.setPen(QColor("#00000040"))
    p.drawEllipse(1, 1, px - 2, px - 2)
    p.end()
    return QIcon(pm)


def _row_text(entry):
    """``label   utilization   verdict`` — the row a member is read off.

    The verdict is on the row because the list is where the members are
    COMPARED, and the distinction that matters between two lines both standing
    at 100% is which of them is yielding at its full tensile capacity and which
    is slipping at what its embedment can develop. A colored dot and a
    percentage cannot carry that; the word can, and it is the same word
    :func:`xslope.fem.print_reinforcement_summary` prints.
    """
    util = entry.get("utilization")
    bits = [entry["label"]]
    if util is not None and util == util:            # not None / NaN
        bits.append(f"{util:.0%}")
    # Reinforcement only. A pile's state names the capacity it was measured
    # against ("at capacity (moment vs Mcap)"), which is longer than the list is
    # wide and would be read truncated; it stays on the row's tooltip, where the
    # reinforcement meanings are too.
    if entry.get("kind") == "reinforcement" and entry.get("status"):
        bits.append(entry["status"])
    return "   ".join(bits)


class FemDetailsDialog(QDialog):
    """List of reinforcement lines and piles on the left, the selected member's
    detail profiles on the right.

    Parameters
    ----------
    fem_data, solution : the solved FEM bundle, as Studio holds it in
        ``doc.results['fem_solution']``.
    slope_data : the model. Supplies the reinforcement capacity envelope's own
        inputs and the pile diameter/spacing the limiting-resistance envelope
        needs; without it the detail views fall back to what ``fem_data`` alone
        can state.
    model_path : the project file, used only for the export dialog's default
        filename.
    failure_solution : the at-failure snapshot SSRM captured, which Studio holds
        beside the solution in the same bundle — the field the Field state
        control switches to. None (a single solve, an SSRM run with no capture,
        or a sidecar saved before the snapshot existed) dims that control.
    """

    def __init__(self, fem_data, solution, slope_data=None, model_path=None,
                 parent=None, failure_solution=None):
        super().__init__(parent)
        from xslope import fem_details

        self._fem_data = fem_data
        self._solution = solution
        self._slope_data = slope_data
        self._model_path = model_path
        self._details = fem_details
        self._failure_solution = fem_details.failure_snapshot(
            solution, failure_solution)
        self._profile = None
        self._mapped = None

        self.setWindowTitle("Reinforcement and pile details")
        self.setModal(False)

        self.list = QListWidget()
        self.list.setIconSize(QSize(12, 12))
        self.list.currentRowChanged.connect(self._on_select)

        # Where the selected member is. A profile says what a member is doing
        # along its length and nothing about where on the slope that is, and the
        # list names members a user has to place before the profile means
        # anything. Same drawing the report's locator prints, with the selection
        # picked out — one helper, so the panel and the page agree.
        self.map_canvas = MplCanvas(self)
        self.map_canvas.setToolTip(
            "Where the selected member lies in the section.")

        left = QSplitter(Qt.Vertical)
        listed = QWidget()
        lv = QVBoxLayout(listed)
        lv.setContentsMargins(0, 0, 0, 0)
        lv.addWidget(self.list, 1)
        left.addWidget(listed)
        left.addWidget(self.map_canvas)
        left.setSizes([360, 240])

        self.canvas = MplCanvas(self)

        # The results view's Field state switch, on the panel that reads the same
        # solution: same two states, same names, same default. It selects the
        # field the member forces and displacements are read from; the capacity
        # envelopes are the model's and the shear band's crossing marks are the
        # captured mechanism's, so neither moves with it. Dimmed — and held on
        # the field there is — when no at-failure snapshot was captured.
        self.field_state = QComboBox()
        for key, label in (("failure", "At failure"), ("converged", "Last converged")):
            self.field_state.addItem(label, key)
        self.field_state.setToolTip(
            "Field the profiles are read from when SSRM captured an at-failure "
            "(unconverged) mechanism:\n"
            "At failure — the member forces in the developed collapse mechanism "
            "(default).\n"
            "Last converged — the sub-critical converged solution instead.\n"
            "The capacity envelopes and the marks are the "
            "same in both.")
        if self._failure_solution is None:
            self.field_state.setCurrentIndex(self.field_state.findData("converged"))
        self.field_state.setEnabled(self._failure_solution is not None)
        self.field_state.currentIndexChanged.connect(self._on_field_state)

        self.bond_chk = QCheckBox("Bond transfer strip")
        self.bond_chk.setChecked(True)
        self.bond_chk.setToolTip(
            "Show the bond transfer rate dT/ds beneath the force profile.")
        self.bond_chk.toggled.connect(self._rerender)
        self.status = QLabel("")
        self.status.setStyleSheet("color: gray;")
        export_btn = QPushButton("Export…")
        export_btn.setToolTip(
            "Save this view as a PNG image and its plotted profiles as a CSV.")
        export_btn.clicked.connect(self.export_current)
        close_btn = QPushButton("Close")
        close_btn.clicked.connect(self.close)

        foot = QHBoxLayout()
        foot.addWidget(QLabel("Field state"))
        foot.addWidget(self.field_state)
        foot.addWidget(self.bond_chk)
        foot.addWidget(self.status, 1)
        foot.addWidget(export_btn)
        foot.addWidget(close_btn)

        right = QWidget()
        rv = QVBoxLayout(right)
        rv.setContentsMargins(0, 0, 0, 0)
        rv.addWidget(self.canvas, 1)
        rv.addLayout(foot)

        splitter = QSplitter(Qt.Horizontal)
        splitter.addWidget(left)
        splitter.addWidget(right)
        splitter.setStretchFactor(0, 0)
        splitter.setStretchFactor(1, 1)
        splitter.setSizes([300, 780])

        outer = QVBoxLayout(self)
        outer.addWidget(splitter)
        self.resize(1080, 620)

        self._entries = []
        self._populate()

    # --- field state -----------------------------------------------------
    def current_field_state(self):
        """``"failure"`` or ``"converged"`` — the field the panel is reading."""
        return self.field_state.currentData()

    def _on_field_state(self, *_):
        """Re-read the whole panel on the newly selected field: the badges and
        the status words are measured on it too, not only the drawn profile."""
        keep = (self._profile or {}).get("kind"), (self._profile or {}).get("index")
        self._populate()
        for row in range(self.list.count()):
            entry = self.list.item(row).data(Qt.UserRole)
            if entry is not None and (entry["kind"], entry["index"]) == keep:
                self.list.setCurrentRow(row)
                break

    # --- list ------------------------------------------------------------
    def _populate(self):
        """Fill the member list: a non-selectable header per kind, then that
        kind's members with their utilization badges."""
        self.list.clear()
        self._entries = self._details.list_lines(
            self._fem_data, self._solution, self._slope_data,
            field_state=self.current_field_state(),
            failure_solution=self._failure_solution)

        seen = set()
        for entry in self._entries:
            if entry["kind"] not in seen:
                seen.add(entry["kind"])
                head = QListWidgetItem(
                    "Reinforcement" if entry["kind"] == "reinforcement" else "Piles")
                head.setFlags(Qt.NoItemFlags)
                f = head.font()
                f.setBold(True)
                head.setFont(f)
                self.list.addItem(head)
            item = QListWidgetItem(_badge_icon(entry["badge"]), _row_text(entry))
            item.setData(Qt.UserRole, entry)
            item.setToolTip(f"{entry['label']} — {entry['status']}")
            self.list.addItem(item)

        for row in range(self.list.count()):
            if self.list.item(row).data(Qt.UserRole) is not None:
                self.list.setCurrentRow(row)
                break

    def entries(self):
        """The member rows this dialog is showing (checks read this)."""
        return list(self._entries)

    def current_profile(self):
        """The profile dict currently drawn, or None."""
        return self._profile

    # --- detail ----------------------------------------------------------
    def _on_select(self, row):
        entry = self.list.item(row).data(Qt.UserRole) if 0 <= row < self.list.count() else None
        if entry is None:
            return
        state = self.current_field_state()
        if entry["kind"] == "reinforcement":
            self._profile = self._details.reinforcement_profile(
                self._fem_data, self._solution, entry["index"], self._slope_data,
                field_state=state, failure_solution=self._failure_solution)
        else:
            self._profile = self._details.pile_profile(
                self._fem_data, self._solution, entry["index"], self._slope_data,
                field_state=state, failure_solution=self._failure_solution)
        self.bond_chk.setVisible(entry["kind"] == "reinforcement")
        self._rerender()
        self._render_map(entry)

    def _render_map(self, entry):
        """Redraw the inset for the selected member.

        The map is of the member's OWN kind — a reinforcement line is placed
        among the reinforcement lines and a pile among the piles — and the
        selected one is the one picked out.
        """
        from xslope.plot_fem_details import plot_member_map
        fem_data, slope_data = self._fem_data, self._slope_data
        kind, index = entry["kind"], entry["index"]

        def _draw(fig):
            # Only the selected member is named. The panel is a column beside a
            # list, and six names in it land on each other's members; the list
            # is where the others are read.
            plot_member_map(fem_data, slope_data, kind, highlight=index,
                            fig=fig, labels=False)

        self._mapped = (kind, index)
        self.map_canvas.render_figure(_draw)

    def mapped_member(self):
        """``(kind, index)`` the inset is showing picked out, or None (checks
        read this)."""
        return getattr(self, "_mapped", None)

    def _fit_window_to_profile(self):
        """Size the window to the profile figure, not the figure to the window.

        The profile keeps the printed-strip rule (its height follows its width),
        so a window taller than that strip shows blank canvas below the plot.
        After the first render, any surplus viewport height is taken off the
        window; the layout's own minimums (the member list and the map on the
        left) bound how far it can shrink, so nothing is clipped.
        """
        if not self.isVisible():
            return
        try:
            self.canvas.render_now()
            fig_h_px = float(self.canvas.figure.get_size_inches()[1]) * BASE_DPI
            vp_h = float(self.canvas.view.viewport().height())
        except Exception:
            return
        surplus = int(vp_h - fig_h_px)
        if surplus > 8:
            self.resize(self.width(), max(self.minimumSizeHint().height(),
                                          self.height() - surplus))

    def _rerender(self, *_):
        prof = self._profile
        if prof is None:
            return
        from xslope.plot_fem_details import plot_detail
        show_bond = self.bond_chk.isChecked()

        def _draw(fig):
            if prof["kind"] == "reinforcement":
                plot_detail(prof, fig=fig, show_bond=show_bond)
            else:
                plot_detail(prof, fig=fig)

        self.canvas.render_figure(_draw)
        QTimer.singleShot(0, self._fit_window_to_profile)
        # The verdict, with what it means one hover away: the panel is where a
        # reader meets these words, and a word with no gloss beside it is a
        # word they have to go and look up.
        from xslope.fem_details import reinforcement_state_meaning
        if prof.get("capture_failed"):
            # There is no at-failure field: this panel is the last converged one,
            # standing in because the capture ran away or never finished. The figure
            # title says exactly that; this line does not contradict it with a
            # verdict the failure state never produced.
            from xslope.plot_fem_details import (CAPTURE_FALLBACK_NOTE,
                                                 CAPTURE_FALLBACK_STATE)
            self.status.setText(f"{CAPTURE_FALLBACK_STATE} — {CAPTURE_FALLBACK_NOTE}")
            self.status.setToolTip(str(prof.get("capture_failed")))
        elif prof.get("capture_truncated"):
            # The at-failure field this panel is reading was stopped mid-runaway
            # (see xslope.fem_details.capture_truncated), so the forces on it are
            # not a state the model settled at and the verdict drawn from them is
            # not one the run made. The figure says the same thing in its title;
            # the two must not disagree over one field.
            from xslope.plot_fem_details import CAPTURE_TRUNCATED_NOTE
            self.status.setText(CAPTURE_TRUNCATED_NOTE)
            self.status.setToolTip(
                "The at-failure capture stopped before the mechanism finished "
                "developing, so no utilization is stated for it. The converged "
                "field is the one to read a verdict from.")
        else:
            self.status.setText(prof.get("status", ""))
            meaning = reinforcement_state_meaning(prof.get("status_key"))
            self.status.setToolTip(f"This {prof['kind']} line {meaning}."
                                   if meaning else "")

    # --- export ----------------------------------------------------------
    def default_export_stem(self):
        """``<model>_<member>_<state>`` — the export dialog's starting filename.

        The state is in the name because the two states of one member are two
        different results: exporting both must not have the second silently
        overwrite the first, and a file found later has to say which it is.
        """
        stem = (os.path.splitext(os.path.basename(self._model_path))[0]
                if self._model_path else "fem")
        label = (self._profile or {}).get("label", "detail")
        safe = "".join(c if c.isalnum() or c in "-_" else "_" for c in str(label))
        state = (self._profile or {}).get("field_state", "converged")
        return f"{stem}_{safe}_{'at_failure' if state == 'failure' else 'converged'}"

    def export_current(self, _checked=False, path=None):
        """Write the current detail as a PNG and its plotted profiles as a CSV.

        Both files share a stem, so the picture and the numbers behind it stay
        together. Returns ``(png_path, csv_path)``, or ``(None, None)`` if the
        user cancelled. ``path`` skips the file dialog (the checks pass it).
        """
        if self._profile is None:
            return None, None
        if path is None:
            start = self.default_export_stem() + ".png"
            if self._model_path:
                start = os.path.join(os.path.dirname(self._model_path), start)
            path, _sel = QFileDialog.getSaveFileName(
                self, "Export detail", start, "PNG image (*.png)")
            if not path:
                return None, None
        base = os.path.splitext(path)[0]
        png_path, csv_path = base + ".png", base + ".csv"
        try:
            self.canvas.render_now()
            self.canvas.figure.savefig(png_path, dpi=200, bbox_inches="tight")
            self._details.write_profile_csv(self._profile, csv_path)
        except Exception:
            import traceback
            traceback.print_exc()
            QMessageBox.warning(self, "Export failed",
                                "Could not write the files — see the Log pane.")
            return None, None
        win = self.parent().window() if self.parent() else None
        if win is not None and hasattr(win, "statusBar"):
            win.statusBar().showMessage(
                f"Exported {os.path.basename(png_path)} and "
                f"{os.path.basename(csv_path)}")
        return png_path, csv_path
