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

"""Generate Report… — where the user composes the Analysis Report.

The dialog is a front end to :func:`xslope.report.generate_report` and holds no
report logic of its own: it collects an output path, a format, the title-page
metadata, which methods the report documents in full, and a checkbox tree of
content, and hands them over as the options dict that function documents.

Two things persist between runs, and they are different kinds of thing. The
organization and the author are properties of the person using the program, so
they are remembered in QSettings and pre-filled every time. The content
selections are remembered the same way, because a user who wants the slice table
left out wants it left out every week. The project title and number are
properties of the PROJECT and are pre-filled from the model's own file name, never
from the last project the user reported on.

The layout follows the Run dialogs: controls in a left column, the content tree
in a right one (:func:`studio.preflight_panel.two_pane`), buttons beneath both.
"""

from __future__ import annotations

import os
import zipfile

from PySide6.QtCore import QSettings, Qt, QUrl
from PySide6.QtGui import QDesktopServices
from PySide6.QtWidgets import (
    QApplication, QCheckBox, QComboBox, QDialog, QDialogButtonBox, QFileDialog,
    QFormLayout, QGroupBox, QHBoxLayout, QLineEdit, QListWidget, QListWidgetItem,
    QMainWindow, QPushButton, QTreeWidget, QTreeWidgetItem, QVBoxLayout, QWidget,
)

from .preflight_panel import two_pane

#: QSettings prefix for everything this dialog remembers.
SETTINGS_PREFIX = "report/"

#: Turn the Word finish (:func:`finalize_document`) off. There is no checkbox
#: for it: a report that arrives complete is not a preference, and the key is
#: here for the one machine where driving Word is unwelcome. On macOS the store
#: :func:`finalization_enabled` reads is the domain Qt writes for this
#: organization and application, so the switch is::
#:
#:     defaults write "com.xslope.XSlope Studio" report.finalize -bool NO
#:
#: The domain is one argument and contains a space, so the quotes are required;
#: without them ``defaults`` reads "XSlope Studio" as the key and refuses the
#: rest of the line.
FINALIZE_KEY = SETTINGS_PREFIX + "finalize"

#: The content tree, as ``(option key, label, tooltip, [children])``. This is the
#: single declaration of what the report can contain: the tree is built from it,
#: the options dict is read back off it, and adding a section means adding a row.
CONTENT_TREE = [
    ("traceability", "Traceability",
     "The xslope version, the input file and its digest, and when the analysis "
     "was run — enough for someone else to reproduce it.", []),
    ("project_definition", "Project definition",
     "The model every analysis in the report is run on.", [
         ("pd_figure", "Model figure",
          "One plot of the shared model: geometry, materials, water surfaces, "
          "loads, reinforcement and piles."),
         ("pd_coords", "Point coordinates on the model figure",
          "Label every geometry point on the model figure with its (x, y), so "
          "the section can be read off the figure. Turn it off where the model "
          "has many closely spaced points and the labels crowd the plot."),
         ("pd_units", "Units statement",
          "Which unit system the numbers in the report are in."),
     ]),
    ("lem", "Limit equilibrium analysis",
     "The stability analysis: its inputs, its search, and its results.", [
         ("lem_inputs_figure", "Limit equilibrium model plot",
          "The section as the method of slices reads it: the materials, the "
          "water and loads on it, and the surface family the analysis is run "
          "over."),
         ("lem_materials", "Materials and pore pressures",
          "The referenced materials with the strength option each is analyzed "
          "under, and where its pore pressure comes from."),
         ("lem_loads", "Loads",
          "Distributed loads as entered, and the seismic coefficient."),
         ("lem_members", "Reinforcement and pile properties",
          "Geometry, capacities, pullout lengths and direction for every "
          "reinforcement line, and the lateral force, capacities and "
          "inclination for every pile — the properties the method of slices "
          "reads. Each gets its own section, and only where the model has one."),
         ("lem_search", "Search documentation",
          "How many trial surfaces the search evaluated, over how many "
          "refinement stages, and within what window."),
         ("lem_search_figure", "Search results plot",
          "Every trial surface the search evaluated, with the critical one "
          "highlighted."),
         ("lem_solution_figure", "Critical surface plot",
          "The critical surface with its slices and base stresses."),
         ("lem_slice_table", "Slice table",
          "Slice geometry, forces and strengths, on a landscape page, with a "
          "legend defining every column."),
         ("lem_slice_key", "Slice key figure",
          "A plot framed on the sliced mass with every slice numbered, placed "
          "immediately above the slice table it keys."),
         ("lem_calculations", "Calculations",
          "The factor of safety worked through: the method's own equation "
          "with the converged numbers in it, and the per-slice terms as "
          "slice-table columns."),
         ("lem_rapid", "Rapid drawdown detail",
          "The three-stage factors of safety, when the run was a rapid "
          "drawdown analysis."),
     ]),
    ("seep", "Seepage analysis",
     "The seepage analysis: the mesh and conductivities it was solved on, the "
     "flow net, and the computed flow. Only where a seepage solution was run.", [
         ("seep_inputs_figure", "Seepage model plot",
          "The section as the flow solver reads it: the material zones, and the "
          "water surface each specified-head boundary states."),
         ("seep_materials", "Seepage material properties",
          "Major and minor conductivity, the angle of the major axis, and the "
          "unsaturated parameters, per material."),
         ("seep_kr_figure", "Unsaturated conductivity plots",
          "Every material's relative conductivity on one set of axes, drawn "
          "twice: against matric suction, and against the pressure head the "
          "solver works in. Only where a material carries an unsaturated "
          "model."),
         ("seep_mesh_figure", "Mesh and boundary condition plot",
          "The mesh colored by material, with every specified-head and exit-face "
          "node marked, for each boundary condition set solved."),
         ("seep_flownet", "Flow net plot",
          "Head contours with the flow lines, the phreatic surface where the "
          "problem is unconfined, and the water level each specified-head "
          "boundary holds."),
         ("seep_variable_figures", "Pore pressure, velocity and gradient plots",
          "The other three fields the solve produced, one plot each: pore "
          "pressure, velocity magnitude with the velocity vectors over it, and "
          "hydraulic gradient magnitude — the same fields the seepage results "
          "view offers. Only the ones the saved solution carries."),
         ("seep_transient_figures", "Transient frame plots",
          "The section at selected states of a transient march — the first, the "
          "last, and the rest spaced evenly through the ones between — each "
          "drawn for the same fields a steady solution is, on one colour scale "
          "across every state. Only where a transient march was run."),
         ("seep_transient_history", "Transient history plot",
          "The march on one time axis: the level the reservoir boundary is held "
          "at, the phreatic elevation and the top of the seepage face above, and "
          "the boundary inflow and outflow below. Only where a transient march "
          "was run."),
     ]),
    ("fem", "Deformation and strength reduction",
     "The finite element analysis: the mesh and properties it was solved on, the "
     "deformation and shear strain it produced, and the strength reduction "
     "factor of safety. Only where an FEM run was made.", [
         ("fem_inputs_figure", "Finite element model plot",
          "The section as the finite element solver reads it: the material "
          "zones, the strength reduction zones, and the members it carries."),
         ("fem_materials", "Finite element material properties",
          "Unit weight, Mohr-Coulomb strength, Young's modulus and Poisson's "
          "ratio, per material."),
         ("fem_loads", "Loads",
          "The distributed loads the deformation analysis carries."),
         ("fem_members", "Reinforcement and pile properties",
          "The modulus, area and second moment of area the members are "
          "assembled from, with the capacities they yield and soften at — the "
          "properties the finite element analysis reads, which are not the ones "
          "the method of slices reads. Only where the model has a member."),
         ("fem_mesh_figure", "Finite element mesh plot",
          "The mesh the section was discretized onto, colored by the material "
          "each element carries, with the boundary conditions on it — the same "
          "plot as the FEM data view."),
         ("fem_figure", "Shear strain, deformed mesh and displacement plots",
          "The viscoplastic shear strain, the deformed mesh and the "
          "displacement vectors — the three plots the finite element results "
          "view offers, in that order."),
         ("fem_state_failure", "Result plots at failure",
          "Draw the result plots from the mechanism at failure — the developed "
          "state a strength reduction run captures beyond its critical factor. "
          "Where the run captured none, the last converged field is drawn "
          "instead and the report says so."),
         ("fem_state_converged", "Result plots at the last converged state",
          "Draw the result plots from the last trial that reached equilibrium. "
          "Tick it alongside the one above to get both, each variable on one "
          "scale across the pair so the two states are read against each other."),
         ("fem_convergence_figure", "Strength reduction search plot",
          "Every trial the strength reduction run solved, at the factor it was "
          "solved at and marked by whether the section stood under it, with the "
          "interval still open after each. Only where the run kept a record of "
          "its trials."),
         ("fem_reinforcement", "Reinforcement forces",
          "Declared capacity, and the axial force and utilization at the point "
          "each reinforcement line is most utilized. Only where the model "
          "carries reinforcement."),
         ("fem_reinforcement_figure", "Reinforcement detail plots",
          "Axial force over the capacity envelope, with the bond transfer "
          "beneath it, along each reinforcement line — one plot per line."),
         ("fem_piles", "Pile forces",
          "Length, peak shear and moment, head movement and utilization for "
          "every pile the run solved. Only where the model carries piles."),
         ("fem_piles_figure", "Pile detail plots",
          "Displacement, shear, bending moment and mobilized soil reaction "
          "with depth, along each pile — one plot per pile."),
     ]),
]

#: Formats offered, in order. Only the enabled ones can be picked; the rest are
#: listed so the roadmap is visible where the choice is made.
FORMATS = [
    ("docx", "Word document (.docx)", True, ""),
    ("pdf", "PDF (.pdf)", False, "PDF reports are not available yet."),
]


def open_output(path, fmt, status=None, finalize=True):
    """Show the finished report to the user, and say what was shown.

    A document is finished in Word first (:func:`finalize_document`) and then
    opened in whatever the system uses for it. Returns ``"document"``.

    ``status`` is where the progress line goes; the default finds Studio's
    status bar. ``finalize=False`` skips the Word finish for a caller that has
    already done it — Studio's report runs on a worker thread and finishes the
    document there, so by the time the document is opened Word has been and
    gone.
    """
    if finalize:
        finalize_document(path, status=status)
    QDesktopServices.openUrl(QUrl.fromLocalFile(path))
    return "document"


def carries_contents_field(path):
    """True when ``path`` is a Word document holding a table-of-contents field.

    Word is only ever handed a document that has something to update. A file
    that is not a ``.docx`` package, or is one without a contents field, is not
    opened in Word to find out — an unreadable file put in front of Word is a
    modal dialog on the user's screen, which is exactly what a cosmetic step
    must never cause.
    """
    try:
        with zipfile.ZipFile(path) as package:
            xml = package.read("word/document.xml").decode("utf-8", "replace")
    except Exception:
        return False
    return "<w:instrText" in xml and " TOC " in xml


def finalization_enabled(settings=None):
    """Whether the Word finish is switched on — it is, unless
    :data:`FINALIZE_KEY` says otherwise."""
    # The app's own store, named as studio.editors names it.
    s = settings or QSettings("XSlope", "XSlope Studio")
    return _as_bool(s.value(FINALIZE_KEY, True))


def word_finish(path, enabled=True):
    """The Word finish itself, with no Qt in it.

    Everything here is a zip read and a subprocess, so this is what runs on the
    report's worker thread (``studio.runners.ReportRunner``). ``enabled`` is
    :func:`finalization_enabled`'s answer, read by the caller where QSettings
    lives. Returns the ``(bool, message)`` of the attempt; the refusals are
    results, not errors.
    """
    if not enabled:
        return False, "Finishing reports in Word is switched off."
    if not carries_contents_field(path):
        return False, "The document has no contents field to update."

    from xslope.report_finalize import finalize_with_word

    ok, msg = finalize_with_word(path)
    if not ok:
        print(f"Report: the page numbers were left for Word to build — {msg}")
    return ok, msg


def finalize_document(path, settings=None, status=None):
    """Let Word rebuild the report's page numbers, and make no fuss if it can't.

    The report's table of contents is a real Word field whose page numbers only
    exist once a page layout engine has laid the document out. Rather than ship
    a contents page that asks the reader to press F9 — or, worse, guess the page
    numbers — the finished document is handed to Word, which updates its fields
    and saves it. What the user then opens is complete.

    Every way this can fail is ordinary: no Word, no permission to drive it, a
    Word that does not answer. Each one leaves the document exactly as it was
    written, which is a report with an unpaginated contents page, and says so in
    one line on the console. None of them is worth a dialog.

    Returns the ``(bool, message)`` of the attempt.

    This is the synchronous form, for a caller already on the GUI thread: it
    blocks for as long as Word takes, so it holds the wait cursor up and posts a
    line saying what is happening. Studio's own report does not come through
    here — it runs on a worker thread and calls :func:`word_finish`, so the
    window stays live and the phase is a labelled stretch of the progress bar.
    """
    enabled = finalization_enabled(settings)
    if not enabled or not carries_contents_field(path):
        return word_finish(path, enabled)   # the same two refusals, unannounced

    say = status or _status_line
    say("Finalizing the report in Word…")
    QApplication.setOverrideCursor(Qt.WaitCursor)
    try:
        return word_finish(path, True)
    finally:
        QApplication.restoreOverrideCursor()


def _status_line(text):
    """Post a progress line where the user is already looking.

    Word takes a few seconds, and a window that sits still for them looks stuck.
    Studio has one main window, so the line is posted to its status bar rather
    than threaded through every caller; with no window it is printed.
    """
    app = QApplication.instance()
    window = None
    if app is not None:
        window = next((w for w in app.topLevelWidgets()
                       if isinstance(w, QMainWindow)), None)
    if window is None:
        print(text)
        return
    window.statusBar().showMessage(text)
    app.processEvents()


def default_output_path(model_path, fmt="docx"):
    """``<model stem>_report.docx`` beside the model, or in the home directory
    for a project that has never been saved."""
    from xslope.report import FORMATS as REPORT_FORMATS
    suffix = REPORT_FORMATS.get(fmt, {}).get("suffix", ".docx")
    if model_path:
        stem = os.path.splitext(model_path)[0]
        return f"{stem}_report{suffix}"
    return os.path.join(os.path.expanduser("~"), f"xslope_report{suffix}")


class ReportDialog(QDialog):
    """Compose an Analysis Report.

    Parameters
    ----------
    slope_data : dict
        The model the report describes. Read only for its defaults (the title).
    solutions : dict
        What has been solved, in :func:`xslope.report.build_report`'s shape. The
        Methods list marks which of them were run.
    model_path : str
        The project file, used for the default output path and the traceability
        stamp.
    default_method : str
        The method the results view is showing — the list opens ticked on it.
    settings : QSettings
        Where the remembered fields live. None means nothing is remembered,
        which is what the offscreen checks construct with.
    """

    def __init__(self, parent=None, slope_data=None, solutions=None,
                 model_path=None, default_method=None, settings=None):
        super().__init__(parent)
        self.setWindowTitle("Generate Report")
        self._sd = slope_data or {}
        self._solutions = solutions or {}
        self._model_path = model_path
        self._settings = settings

        layout = QVBoxLayout()                   # the left column; see two_pane

        # --- output ---------------------------------------------------------
        out_box = QGroupBox("Output")
        out_form = QFormLayout(out_box)

        self.format = QComboBox()
        for key, label, enabled, tip in FORMATS:
            self.format.addItem(label, key)
            i = self.format.count() - 1
            if not enabled:
                self.format.model().item(i).setEnabled(False)
                self.format.setItemData(i, tip, Qt.ToolTipRole)
        self.format.currentIndexChanged.connect(self._on_format)
        out_form.addRow("Format", self.format)

        row = QHBoxLayout()
        self.path = QLineEdit(default_output_path(model_path))
        # Wide enough for a real path in whatever font and scaling the platform
        # is using, measured rather than guessed in pixels.
        self.path.setMinimumWidth(
            self.path.fontMetrics().horizontalAdvance("n" * 48))
        browse = QPushButton("Browse…")
        browse.clicked.connect(self._browse)
        row.addWidget(self.path, 1)
        row.addWidget(browse)
        holder = QWidget()
        holder.setLayout(row)
        row.setContentsMargins(0, 0, 0, 0)
        out_form.addRow("Save as", holder)

        layout.addWidget(out_box)

        # --- what is reported ------------------------------------------------
        run_box = QGroupBox("Analysis")
        run_form = QFormLayout(run_box)
        self.methods = QListWidget()
        self.methods.setToolTip(
            "Which methods the report documents in full — a results statement, a "
            "critical-surface plot, a slice table and the calculations, for each "
            "one ticked. Tick two or more and the report also carries a table of "
            "their factors of safety side by side.")
        self._fill_methods(default_method)
        self.methods.itemChanged.connect(self._on_method_changed)
        run_form.addRow("Methods", self.methods)
        layout.addWidget(run_box)

        # --- title page ------------------------------------------------------
        meta_box = QGroupBox("Title page")
        meta = QFormLayout(meta_box)
        self.title = QLineEdit(self._default_title())
        self.project_number = QLineEdit()
        self.organization = QLineEdit()
        self.author = QLineEdit()
        meta.addRow("Project title", self.title)
        meta.addRow("Project number", self.project_number)
        meta.addRow("Organization", self.organization)
        meta.addRow("Author", self.author)
        self.signature_lines = QCheckBox("Signature lines (prepared by / checked by)")
        self.signature_lines.setToolTip(
            "Add ruled signature lines to the title page.")
        meta.addRow("", self.signature_lines)
        layout.addWidget(meta_box)
        layout.addStretch(1)

        # --- content tree ----------------------------------------------------
        self.contents = QTreeWidget()
        self.contents.setHeaderLabel("Contents")
        self.contents.setMinimumWidth(340)
        self._items = {}
        for key, label, tip, children in CONTENT_TREE:
            parent_item = self._add_item(self.contents, key, label, tip)
            for child_key, child_label, child_tip in children:
                self._add_item(parent_item, child_key, child_label, child_tip)
        self.contents.expandAll()
        self.contents.itemChanged.connect(self._on_item_changed)

        bb = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        self._ok = bb.button(QDialogButtonBox.Ok)
        self._ok.setText("Generate")
        bb.accepted.connect(self.accept)
        bb.rejected.connect(self.reject)

        two_pane(self, layout, self.contents, bb, checks_width=340)

        self.restore()
        self._on_format()
        for key in self._items:
            self._sync_children(key)

    # --- construction helpers -------------------------------------------
    def _add_item(self, parent, key, label, tip):
        """One checkbox, opening on the option's own default.

        The default comes from :data:`xslope.report.DEFAULT_OPTIONS` rather than
        from a second list here: a box that opens checked while the builder's
        default is off would produce a report nobody asked for, and the two would
        drift the first time one of them changed.
        """
        from xslope.report import DEFAULT_OPTIONS

        item = QTreeWidgetItem(parent, [label])
        item.setFlags(item.flags() | Qt.ItemIsUserCheckable)
        item.setCheckState(0, Qt.Checked if DEFAULT_OPTIONS.get(key, True)
                           else Qt.Unchecked)
        item.setToolTip(0, tip)
        item.setData(0, Qt.UserRole, key)
        self._items[key] = item
        return item

    def _default_title(self):
        """The project's own name, not the last project's: the model file's stem,
        tidied into words."""
        if not self._model_path:
            return "Slope Stability Analysis"
        stem = os.path.splitext(os.path.basename(self._model_path))[0]
        return stem.replace("_", " ").replace("-", " ").strip().title()

    def _fill_methods(self, default_method):
        """One checkable row per method the solver offers, opening on the one the
        results view is showing.

        EVERY method is offered, not only the ones that have been run, and each
        row says which it is. Ticking a method that has not been run RUNS it,
        exactly as ticking it in the Run dialog would have: a search for its own
        critical surface where this model's analysis searched, the specified
        surface where it did not (the owner's ruling;
        :func:`~xslope.report.run_requested_methods` does the running). So every
        tick means the same thing, and the rows that say "not run yet" are saying
        how long the report will take, not that they are worth less. A method
        that cannot run on this surface family at all is listed and disabled,
        with the reason on it, rather than quietly missing.
        """
        from xslope.preflight import method_surface_reason
        from xslope.report import (DEFAULT_METHOD, analysis_options, method_label,
                                   solved_methods, supported_methods,
                                   surface_family)

        family = surface_family(self._sd, self._solutions)
        run = solved_methods(self._solutions)
        # What ticking an unrun method costs, in the terms of THIS model's
        # analysis: the report runs it the way the run was made
        # (:func:`~xslope.report.analysis_options`), which is a search on a model
        # that searched and one solve on a model that did not.
        will_run = ("Not run yet; the report will run it — a search for this "
                    "method's own critical surface, which takes as long as a run."
                    if analysis_options(self._sd, self._solutions)["analysis"]
                    == "auto_search" else
                    "Not run yet; the report will run it, on the surface the "
                    "input specifies.")

        # What opens ticked: the results view's own method, else a method that
        # was run, else the default — the first of those this SURFACE can take.
        # On a non-circular surface the results view may be showing a method that
        # is dimmed here, and the box that opens ticked has to be one the user
        # could have ticked themselves.
        usable = [n for n in supported_methods()
                  if not method_surface_reason(n, family)]
        wanted = [(default_method or "").lower()] + run + [DEFAULT_METHOD]
        opening = next((n for n in wanted if n in usable),
                       usable[0] if usable else "")
        self._opening = opening

        self.methods.blockSignals(True)
        for name in supported_methods():
            item = QListWidgetItem(method_label(name))
            item.setData(Qt.UserRole, name)
            reason = method_surface_reason(name, family)
            if reason:
                item.setFlags(Qt.NoItemFlags)
                item.setCheckState(Qt.Unchecked)
                item.setToolTip(reason)
            else:
                item.setFlags(item.flags() | Qt.ItemIsUserCheckable)
                item.setCheckState(Qt.Checked if name == opening else Qt.Unchecked)
                item.setToolTip("Solved in this run." if name in run else will_run)
            self.methods.addItem(item)
        self.methods.blockSignals(False)
        self._ensure_a_method()
        # Tall enough for every method, so the list never scrolls a choice out of
        # sight: measured from the rows themselves rather than set in pixels.
        rows = self.methods.count()
        if rows:
            self.methods.setFixedHeight(
                self.methods.sizeHintForRow(0) * rows
                + 2 * self.methods.frameWidth())

    def _method_items(self):
        return [self.methods.item(i) for i in range(self.methods.count())]

    def _on_method_changed(self, _item):
        self._ensure_a_method()

    def _ensure_a_method(self):
        """A report documents at least one method: unticking the last one ticks
        it straight back rather than producing a report with no results in it.

        What comes back is the method the list opened on — which is also what a
        remembered selection falls back to when this surface can take none of the
        methods it names.
        """
        items = [i for i in self._method_items()
                 if i.flags() & Qt.ItemIsUserCheckable]
        if any(i.checkState() == Qt.Checked for i in items) or not items:
            return
        want = next((i for i in items
                     if i.data(Qt.UserRole) == getattr(self, "_opening", None)),
                    items[0])
        self.methods.blockSignals(True)
        want.setCheckState(Qt.Checked)
        self.methods.blockSignals(False)

    def selected_methods(self):
        """The methods ticked, in the order the report will document them."""
        return [i.data(Qt.UserRole) for i in self._method_items()
                if i.checkState() == Qt.Checked]

    # --- behavior --------------------------------------------------------
    def _on_format(self, *_):
        """Keep the output path's suffix in step with the chosen format."""
        from xslope.report import FORMATS as REPORT_FORMATS
        fmt = self.format.currentData()
        suffix = REPORT_FORMATS.get(fmt, {}).get("suffix", ".docx")
        current = self.path.text().strip()
        if current:
            self.path.setText(os.path.splitext(current)[0] + suffix)

    def _browse(self):
        fmt = self.format.currentData()
        label = dict((k, l) for k, l, _e, _t in FORMATS).get(fmt, "Document")
        suffix = os.path.splitext(self.path.text())[1] or ".docx"
        path, _sel = QFileDialog.getSaveFileName(
            self, "Save report as", self.path.text(), f"{label} (*{suffix})")
        if path:
            self.path.setText(path)

    def _on_item_changed(self, item, _column):
        key = item.data(0, Qt.UserRole)
        if key in self._items and self._items[key] is item:
            self._sync_children(key)

    def _sync_children(self, key):
        """A section that is not in the report has nothing to configure, so its
        sub-items dim with it — visibly off rather than silently ignored."""
        item = self._items[key]
        on = item.checkState(0) == Qt.Checked
        for i in range(item.childCount()):
            item.child(i).setDisabled(not on)

    # --- what the caller reads -------------------------------------------
    def output_path(self):
        return self.path.text().strip()

    def output_format(self):
        return self.format.currentData()

    def options(self):
        """The options dict for :func:`xslope.report.generate_report`."""
        opts = {
            "title": self.title.text().strip(),
            "project_number": self.project_number.text().strip(),
            "organization": self.organization.text().strip(),
            "author": self.author.text().strip(),
            "signature_lines": self.signature_lines.isChecked(),
            "method": self.selected_methods(),
            "input_path": self._model_path,
        }
        for key, item in self._items.items():
            opts[key] = item.checkState(0) == Qt.Checked
        # A sub-item of a section that is off is off, whatever its own box says.
        for key, _label, _tip, children in CONTENT_TREE:
            if not opts.get(key, True):
                for child_key, _l, _t in children:
                    opts[child_key] = False
        return opts

    # --- persistence ------------------------------------------------------
    #: Fields remembered between runs: who the user is, and what they want in it.
    _REMEMBERED_TEXT = ("organization", "author")

    def restore(self):
        s = self._settings
        if s is None:
            return
        for name in self._REMEMBERED_TEXT:
            value = s.value(SETTINGS_PREFIX + name, "")
            if value:
                getattr(self, name).setText(str(value))
        fmt = s.value(SETTINGS_PREFIX + "format", "")
        i = self.format.findData(fmt) if fmt else -1
        if i >= 0 and self.format.model().item(i).isEnabled():
            self.format.setCurrentIndex(i)
        self.signature_lines.setChecked(
            _as_bool(s.value(SETTINGS_PREFIX + "signature_lines", False)))
        # The methods are remembered like the content boxes: an engineer who
        # reports Spencer and Bishop side by side reports them side by side every
        # week. The results view's own method still opens ticked where nothing was
        # stored, and a stored method this surface cannot take is skipped.
        stored = s.value(SETTINGS_PREFIX + "methods", None)
        if stored:
            want = set(stored if isinstance(stored, (list, tuple))
                       else str(stored).split(","))
            self.methods.blockSignals(True)
            for item in self._method_items():
                if item.flags() & Qt.ItemIsUserCheckable:
                    item.setCheckState(
                        Qt.Checked if item.data(Qt.UserRole) in want else Qt.Unchecked)
            self.methods.blockSignals(False)
            self._ensure_a_method()
        for key, item in self._items.items():
            stored = s.value(SETTINGS_PREFIX + "content/" + key, None)
            if stored is not None:
                item.setCheckState(0, Qt.Checked if _as_bool(stored) else Qt.Unchecked)

    def remember(self):
        """Write the remembered fields back. Called on accept."""
        s = self._settings
        if s is None:
            return
        for name in self._REMEMBERED_TEXT:
            s.setValue(SETTINGS_PREFIX + name, getattr(self, name).text().strip())
        s.setValue(SETTINGS_PREFIX + "format", self.format.currentData())
        s.setValue(SETTINGS_PREFIX + "signature_lines",
                   self.signature_lines.isChecked())
        s.setValue(SETTINGS_PREFIX + "methods", self.selected_methods())
        for key, item in self._items.items():
            s.setValue(SETTINGS_PREFIX + "content/" + key,
                       item.checkState(0) == Qt.Checked)

    def accept(self):
        if not self.output_path():
            from PySide6.QtWidgets import QMessageBox
            QMessageBox.warning(self, "Generate Report",
                                "Choose where the report should be saved.")
            return
        self.remember()
        super().accept()


def _as_bool(value):
    """QSettings hands back "true"/"false" strings on some platforms and real
    bools on others."""
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in ("true", "1", "yes")
