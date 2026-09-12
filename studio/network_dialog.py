"""The joints editor's "Build network…" dialog.

A jointed rock mass is described as SETS — bedding at a dip and a spacing, a
conjugate set crossing it, a blocky mass with no preferred orientation — and not
as a hundred separate lines. This dialog is that description: a kind, its
parameters, where the set exists, and one set of joint properties for every trace
in it. OK writes the rows onto the joints sheet.

The rows carry the description back. Every row's Label holds the set record (see
:mod:`xslope.joints`), so selecting one of a set's rows in the joints editor and
pressing the button again REOPENS this dialog on that set, filled in as it was
built, and OK regenerates it in place: changing a spacing is an edit, not a
delete-and-redo.

The preview draws what OK would write, before it writes it — the traces in the
joint color over the section, with the model's existing joint lines behind them
in gray — and the count under it says how many lines the parameters resolve to.
A set that cannot be generated at all (a block size larger than its region, two
sets at one dip) says so there instead of raising on OK.
"""

from __future__ import annotations

from PySide6.QtCore import Qt, QTimer
from PySide6.QtGui import QGuiApplication
from PySide6.QtWidgets import (QComboBox, QDialog, QDialogButtonBox, QFormLayout,
                               QGridLayout, QGroupBox, QHBoxLayout, QLabel,
                               QLineEdit, QScrollArea, QSplitter, QStackedWidget,
                               QVBoxLayout, QWidget)

from xslope.joints import JointSet, set_name

#: The kinds, in the order the combo offers them: (word, JointSet kind).
KIND_ITEMS = [("Parallel set", "parallel"),
              ("Cross-jointed (two sets)", "cross"),
              ("Voronoi (blocky mass)", "voronoi")]

#: The parameter fields of each kind: (key, label, tooltip). What the dialog
#: shows is exactly what :class:`xslope.joints.JointSet` takes.
KIND_FIELDS = {
    "parallel": [
        ("dip", "Dip (deg)",
         "The traces' inclination, counter-clockwise from horizontal: 0 is "
         "horizontal, 90 is vertical, and a set descending to the right is a "
         "NEGATIVE angle."),
        ("spacing", "Spacing",
         "The perpendicular distance between neighboring traces."),
        ("offset", "Offset",
         "Where the set sits across its own normal. One trace passes through "
         "the origin offset by this much, so 0 puts a trace through (0, 0) "
         "extended and half a spacing shifts the whole set. It moves the set; "
         "it does not change how many traces the region gets."),
        ("trace_len", "Trace length",
         "A discontinuous set: each trace exists for this length, stops, and "
         "resumes. Blank is a fully persistent set, every trace continuous "
         "across the region."),
        ("gap", "Gap",
         "The rock bridge between one trace segment and the next. Read only "
         "when a trace length is given."),
    ],
    "cross": [
        ("dip", "Dip 1 (deg)", "The first set's inclination, as above."),
        ("spacing", "Spacing 1", "The first set's perpendicular spacing."),
        ("offset", "Offset 1", "Where the first set sits across its normal."),
        ("dip2", "Dip 2 (deg)",
         "The second set's inclination. A conjugate pair is one dip and its "
         "negative, 60 and -60. The two dips must differ: two sets at one dip "
         "lie on one another, and the mesh split has no material between them."),
        ("spacing2", "Spacing 2", "The second set's perpendicular spacing."),
        ("offset2", "Offset 2", "Where the second set sits across its normal."),
        ("trace_len", "Trace length",
         "A discontinuous network: each trace of BOTH sets exists for this "
         "length, stops, and resumes. Blank is fully persistent."),
        ("gap", "Gap", "The rock bridge between one trace segment and the next."),
    ],
    "voronoi": [
        ("block_size", "Block size",
         "The characteristic cell size: a rock mass with no through-going set "
         "is described by how big its blocks are rather than by a dip and a "
         "spacing. Every cell wall becomes a joint."),
        ("seed", "Seed",
         "The random seed. Required, not optional: the same block size and "
         "seed reproduce the same network exactly, and a tessellation nobody "
         "can reproduce is not an input."),
    ],
}

#: Which parameters a kind may leave blank, and what blank means.
_OPTIONAL = {"offset", "offset2", "trace_len", "gap"}

#: The joint properties every trace of the set gets, laid out two per row in the
#: joints sheet's own pairings. The help text is the joints editor's, imported
#: rather than restated.
_PROP_ROWS = [("c", "phi"), ("c_res", "phi_res"), ("dil", "t_cut"),
              ("kn", "ks"), ("jred", None)]


def _f(text):
    """A number typed into a field, or ``None`` for a blank one."""
    text = (text or "").strip()
    if not text:
        return None
    return float(text)


def _num_or_none(text):
    try:
        return _f(text)
    except ValueError:
        return None


class BuildNetworkDialog(QDialog):
    """Build (or rebuild) one joint set.

    Parameters
    ----------
    slope_data : dict
        The model, read for its geometry, its materials and its joint regions.
    rows : list of dict
        The joint lines as the editor currently holds them — what the preview
        draws behind the new set, and what a name collision is checked against.
    jset : JointSet, optional
        An existing set to reopen. ``None`` builds a new one.
    parent : QWidget, optional
    """

    def __init__(self, slope_data, rows=None, jset=None, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Build joint network")
        self._sd = slope_data or {}
        self._rows = [dict(r) for r in (rows or [])]
        self._editing = jset.name if jset is not None else None
        self._generated = []
        self._error = ""
        self._help_by_widget = {}
        # Regeneration is debounced (see _schedule); the interval is the preview
        # pane's own, so the picture and the count arrive together.
        self._timer = QTimer(self)
        self._timer.setSingleShot(True)
        self._timer.setInterval(160)
        self._timer.timeout.connect(self._refresh)

        layout = QVBoxLayout(self)
        from .editors import _help_label
        layout.addWidget(_help_label(
            "A joint network is described as a SET — a dip and a spacing, two "
            "sets crossing, or a block size — and resolved into the individual "
            "traces the joints sheet holds, clipped to where the set exists. "
            "Every row it writes carries the set's parameters in its label, so "
            "this dialog reopens on the set and regenerates it when something "
            "changes. The preview shows what OK would write."))

        split = QSplitter(Qt.Horizontal)
        split.addWidget(self._build_form())
        split.addWidget(self._build_preview())
        split.setStretchFactor(0, 0)
        split.setStretchFactor(1, 1)
        layout.addWidget(split, 1)

        self._status = QLabel("")
        self._status.setWordWrap(True)
        layout.addWidget(self._status)

        self._buttons = QDialogButtonBox(QDialogButtonBox.Ok
                                         | QDialogButtonBox.Cancel)
        self._buttons.accepted.connect(self.accept)
        self._buttons.rejected.connect(self.reject)
        layout.addWidget(self._buttons)

        # The same context-sensitive help strip the editors carry, over the same
        # property help text — a reader meets one explanation of kn, not two.
        from .editors import attach_help
        attach_help(self, self._help_text(), self._help_key_for)

        # The preview opens at least as wide as the form beside it: a picture
        # narrower than the fields it illustrates is not a preview of anything.
        self._preview.setMinimumWidth(self._form_scroll.minimumWidth())
        split.setSizes([self._form_scroll.minimumWidth(),
                        self._form_scroll.minimumWidth()])

        if jset is not None:
            self._load(jset)
        else:
            self._name.setText(self._next_name())
        self._on_kind()
        self._refresh()
        # Opened at the size its own contents ask for rather than at a
        # remembered one: the form's measured height, and a preview wide enough
        # to show the whole section at that height — a cross section is wide and
        # short, and an equal-aspect picture of one in a narrow pane is mostly
        # empty. Capped by the screen, which is the only limit that is not a
        # number somebody chose.
        self.adjustSize()
        try:
            x0, y0, x1, y1 = self._sd["domain_polygon"].bounds
            aspect = max((x1 - x0) / max(y1 - y0, 1e-9), 1.0)
        except Exception:
            aspect = 1.0
        screen = self.screen() or QGuiApplication.primaryScreen()
        want = self._form_scroll.minimumWidth() + int(
            aspect * self._form_scroll.minimumHeight())
        if screen is not None:
            want = min(want, int(0.9 * screen.availableGeometry().width()))
        self.resize(max(want, self.width()), self.height())

    # -- construction ------------------------------------------------------
    def _edit(self, key, tooltip):
        w = QLineEdit()
        w.setToolTip(tooltip)
        w.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
        w.textChanged.connect(self._schedule)
        self._help_by_widget[w] = key
        return w

    def _build_form(self):
        from .editors import JOINTS_HELP
        holder = QWidget()
        form = QVBoxLayout(holder)
        form.setContentsMargins(0, 0, 0, 0)

        # -- the set itself
        g_set = QGroupBox("Set")
        f_set = QFormLayout(g_set)
        self._name = QLineEdit()
        self._name.setToolTip(
            "The set's name. Every row it writes is labeled with it, so the set "
            "can be found, reopened and regenerated as one thing.")
        self._name.textChanged.connect(self._schedule)
        f_set.addRow("Name", self._name)
        self._kind = QComboBox()
        for word, _kind in KIND_ITEMS:
            self._kind.addItem(word)
        self._kind.setToolTip(
            "What kind of network this is. A parallel set is one orientation at "
            "one spacing (bedding, a joint set); cross-jointed is two of them "
            "crossing; Voronoi is a blocky mass with no preferred orientation.")
        self._kind.currentIndexChanged.connect(self._on_kind)
        f_set.addRow("Kind", self._kind)
        form.addWidget(g_set)

        # -- the parameters of whichever kind is chosen
        g_par = QGroupBox("Parameters")
        v_par = QVBoxLayout(g_par)
        self._stack = QStackedWidget()
        self._param_edits = {}
        for _word, kind in KIND_ITEMS:
            page = QWidget()
            pf = QFormLayout(page)
            edits = {}
            for key, label, tip in KIND_FIELDS[kind]:
                w = self._edit(key, tip)
                if key in _OPTIONAL:
                    w.setPlaceholderText("0" if key in ("offset", "offset2", "gap")
                                         else "continuous")
                edits[key] = w
                pf.addRow(label, w)
            self._param_edits[kind] = edits
            self._stack.addWidget(page)
        v_par.addWidget(self._stack)
        form.addWidget(g_par)

        # -- where the set exists
        g_reg = QGroupBox("Region")
        f_reg = QFormLayout(g_reg)
        self._region = QComboBox()
        self._region_specs = []
        for text, spec in self._region_items():
            self._region.addItem(text)
            self._region_specs.append(spec)
        self._region.setToolTip(
            "Where the set exists. The whole section, one material (every zone "
            "carrying it), or a polygon of Type 'joints' drawn on the polygon "
            "sheet — which is how a set is confined to ground no material "
            "boundary draws.")
        self._region.currentIndexChanged.connect(self._schedule)
        f_reg.addRow("Within", self._region)
        band = QWidget()
        hb = QHBoxLayout(band)
        hb.setContentsMargins(0, 0, 0, 0)
        self._band_lo = self._edit("band", _BAND_HELP)
        self._band_hi = self._edit("band", _BAND_HELP)
        self._band_lo.setPlaceholderText("open")
        self._band_hi.setPlaceholderText("open")
        hb.addWidget(QLabel("from"))
        hb.addWidget(self._band_lo, 1)
        hb.addWidget(QLabel("to"))
        hb.addWidget(self._band_hi, 1)
        f_reg.addRow("Elevation band", band)
        form.addWidget(g_reg)

        # -- one property set for the whole network
        g_prop = QGroupBox("Joint properties (every trace in the set)")
        grid = QGridLayout(g_prop)
        self._prop_edits = {}
        for r, pair in enumerate(_PROP_ROWS):
            for half, key in enumerate(pair):
                if key is None:
                    continue
                col = 2 * half
                grid.addWidget(QLabel(_PROP_LABELS[key]), r, col)
                if key == "jred":
                    w = QComboBox()
                    w.addItems(["", "Yes", "No"])
                    w.setToolTip(JOINTS_HELP[key])
                    w.currentIndexChanged.connect(self._schedule)
                else:
                    w = self._edit(key, JOINTS_HELP[key])
                self._help_by_widget[w] = key
                self._prop_edits[key] = w
                grid.addWidget(w, r, col + 1)
        grid.setColumnStretch(1, 1)
        grid.setColumnStretch(3, 1)
        form.addWidget(g_prop)
        form.addStretch(1)

        scroll = QScrollArea()
        scroll.setWidget(holder)
        scroll.setWidgetResizable(True)
        # The form's own measured size, so no field of it is cut off and none of
        # the dialog is spent on space the form does not use. The height is
        # capped at a share of the screen rather than at a number: on a short
        # display the form scrolls, on a tall one it opens whole.
        hint = holder.sizeHint()
        scroll.setMinimumWidth(hint.width()
                               + scroll.verticalScrollBar().sizeHint().width())
        screen = self.screen()
        room = screen.availableGeometry().height() if screen is not None else 0
        scroll.setMinimumHeight(min(hint.height(), int(0.7 * room)) if room
                                else hint.height())
        self._form_scroll = scroll
        return scroll

    def _build_preview(self):
        from .canvas import PreviewPane
        self._preview = PreviewPane(
            self._draw, caption="Preview shows the traces this set resolves to "
                                "(joint color) over the model's existing joint "
                                "lines (gray). Nothing is written until OK.")
        return self._preview

    # -- the model's own vocabulary ---------------------------------------
    def _region_items(self):
        """(what the combo says, what the record stores) for every region."""
        items = [("Whole section", None)]
        mats = self._sd.get("materials") or []
        used = {int(p.get("mat_id")) for p in (self._sd.get("polygons") or [])
                if p.get("mat_id") is not None}
        for i, m in enumerate(mats):
            name = str(m.get("name") or "").strip()
            if not name or i not in used:
                continue
            items.append((f"Material: {name}", f"mat:{name}"))
        for n, z in enumerate(self._sd.get("joint_zones") or []):
            label = str(z.get("label") or "").strip()
            items.append((f"Joint region: {label}" if label
                          else f"Joint region {n + 1} (unnamed)",
                          f"poly:{label}" if label else f"poly:{n + 1}"))
        return items

    def _next_name(self):
        """A name no set in the model is using: set1, set2, …"""
        taken = {set_name(r.get("label")) for r in self._rows}
        n = 1
        while f"set{n}" in taken:
            n += 1
        return f"set{n}"

    # -- help --------------------------------------------------------------
    def _help_text(self):
        from .editors import JOINTS_HELP
        out = dict(JOINTS_HELP)
        out["band"] = _BAND_HELP
        for kind, fields in KIND_FIELDS.items():
            for key, _label, tip in fields:
                out.setdefault(key, tip)
        return out

    def _help_key_for(self, widget):
        return self._help_by_widget.get(widget)

    # -- state -------------------------------------------------------------
    def _kind_at(self):
        i = self._kind.currentIndex()
        return KIND_ITEMS[i][1] if 0 <= i < len(KIND_ITEMS) else "parallel"

    def _on_kind(self):
        self._stack.setCurrentIndex(max(self._kind.currentIndex(), 0))
        self._refresh()

    def _band(self):
        lo, hi = _num_or_none(self._band_lo.text()), _num_or_none(self._band_hi.text())
        if lo is None and hi is None:
            return None
        return (lo, hi)

    def _props(self):
        """The property set, blanks left out so the sheet's own defaults apply."""
        out = {}
        for key, w in self._prop_edits.items():
            if key == "jred":
                text = w.currentText().strip()
                if text:
                    out[key] = text
                continue
            value = _num_or_none(w.text())
            if value is not None:
                out[key] = value
        return out

    def result_set(self):
        """The set the dialog describes, or ``None`` when it does not describe one."""
        kind = self._kind_at()
        params = {}
        for key, w in self._param_edits[kind].items():
            value = _num_or_none(w.text())
            if value is None:
                if key in _OPTIONAL:
                    continue
                return None
            params[key] = value
        try:
            return JointSet(self._name.text(), kind, params,
                            region=self._region_specs[self._region.currentIndex()],
                            band=self._band(), props=self._props())
        except (ValueError, IndexError):
            return None

    def _load(self, jset):
        """Fill the form in from an existing set."""
        from .editors import _display_number
        self._name.setText(jset.name)
        for i, (_word, kind) in enumerate(KIND_ITEMS):
            if kind == jset.kind:
                self._kind.setCurrentIndex(i)
        for key, w in self._param_edits[jset.kind].items():
            w.setText(_display_number(jset.params.get(key)))
        spec = jset.region if isinstance(jset.region, str) else None
        if spec in self._region_specs:
            self._region.setCurrentIndex(self._region_specs.index(spec))
        elif spec is not None:
            # A region the model no longer carries — a renamed joint region, a
            # deleted material. Offered as itself so reopening the set does not
            # silently re-point it at the whole section.
            self._region.addItem(f"{spec} (not in this model)")
            self._region_specs.append(spec)
            self._region.setCurrentIndex(len(self._region_specs) - 1)
        if jset.band is not None:
            lo, hi = jset.band
            self._band_lo.setText(_display_number(lo))
            self._band_hi.setText(_display_number(hi))
        for key, w in self._prop_edits.items():
            value = jset.props.get(key)
            if key == "jred":
                text = str(value or "").strip().capitalize()
                w.setCurrentIndex(max(0, w.findText(text)) if text else 0)
            else:
                w.setText(_display_number(value))

    # -- the live preview --------------------------------------------------
    def _schedule(self, *_args):
        """A field changed: regenerate shortly, not on this keystroke.

        Generating is the expensive part — a Voronoi network over a section is a
        tessellation of thousands of cells — and a half-typed block size is a
        network nobody asked for. The debounce is the one the preview pane
        already uses for drawing, applied one step earlier so the typing itself
        stays responsive. Anything that READS the result flushes it first, so
        nothing ever sees a stale count.
        """
        self._timer.start()

    def _flush(self):
        """Regenerate now if a change is still pending."""
        self._timer.stop()
        self._refresh()

    def accept(self):
        """OK: never on a set the pending edits have not been generated yet."""
        self._flush()
        if not self._generated:
            return
        super().accept()

    def _refresh(self, *_args):
        """Regenerate the set for the preview and say what it came to."""
        jset = self.result_set()
        self._generated, self._error = [], ""
        if jset is None:
            self._error = ("Fill in the set's name and every parameter above.")
        else:
            model = dict(self._sd)
            model["joint_lines"] = self._rows
            try:
                self._generated = jset.generate(model)
            except Exception as exc:                # a refusal, shown not raised
                self._error = str(exc)
        clash = self._name_clash(jset)
        ok = bool(self._generated) and not clash
        self._buttons.button(QDialogButtonBox.Ok).setEnabled(ok)
        if clash:
            self._status.setText(clash)
        elif self._error:
            self._status.setText(self._error)
        else:
            note = ""
            if "phi" not in self._props():
                note = ("  The friction angle is blank; preflight refuses a run "
                        "on a joint that has none.")
            self._status.setText(
                f"{len(self._generated)} joint line"
                f"{'' if len(self._generated) == 1 else 's'}"
                f"{'' if self._editing is None else ', replacing the set'}."
                + note)
        self._preview.schedule()

    def _name_clash(self, jset):
        """The message for a name another set is already using, or ``''``."""
        if jset is None:
            return ""
        taken = {set_name(r.get("label")) for r in self._rows} - {""}
        if self._editing is not None:
            taken -= {self._editing}
        if jset.name in taken:
            return (f"This model already has a joint set named "
                    f"{jset.name!r}. Give this one a different name, or select "
                    f"one of that set's rows and press Build network to edit it.")
        return ""

    def _draw(self, ax):
        from xslope.plot import (JOINT_COLOR, draw_joint_ticks, joint_linestyle,
                                 joint_ticks_wanted, plot_base_geometry,
                                 plot_joint_zones)
        from xslope.style import resolve_style
        from .editors import _doc_style, _finish_preview_axes, _preview_span
        style = resolve_style(_doc_style(self.parent()))
        try:
            plot_base_geometry(ax, self._sd, labels=False, style=style)
            plot_joint_zones(ax, self._sd, style=style)
        except Exception:
            pass
        span = _preview_span(ax, self._sd)
        # The model's existing joint lines, behind: what is already there, so a
        # second set is placed against the first rather than in the dark.
        mine = self._editing
        for r in self._rows:
            if mine is not None and set_name(r.get("label")) == mine:
                continue                     # this set is being replaced
            try:
                xs = [float(r["x1"]), float(r["x2"])]
                ys = [float(r["y1"]), float(r["y2"])]
            except (KeyError, TypeError, ValueError):
                continue
            ax.plot(xs, ys, color="gray", linewidth=1.0, alpha=0.5, zorder=4)
        ticks = joint_ticks_wanted(len(self._generated))
        dash = joint_linestyle(len(self._generated))
        for r in self._generated:
            xs, ys = [r["x1"], r["x2"]], [r["y1"], r["y2"]]
            ax.plot(xs, ys, color=JOINT_COLOR, linewidth=1.6, linestyle=dash,
                    alpha=0.95, zorder=6)
            if ticks:
                draw_joint_ticks(ax, xs, ys, span, JOINT_COLOR, alpha=0.95,
                                 zorder=6)
        _finish_preview_axes(ax)

    # -- results -----------------------------------------------------------
    def result_rows(self):
        """The rows OK writes: the set, generated and labeled with its record."""
        self._flush()
        return [dict(r) for r in self._generated]

    def editing(self):
        """The name of the set being replaced, or ``None`` for a new one."""
        return self._editing


#: Shared by the band's two fields — one explanation, two widgets.
_BAND_HELP = ("An elevation band the region is cut to, which is how 'the "
              "sandstone above elevation 40' is asked for. Either end may be "
              "left blank for open. A trace lying along the band's own edge is "
              "kept: a band is a line drawn through material, and the mesh "
              "split has material on both sides of it.")

#: The property fields' labels — the joints sheet's own column headers.
_PROP_LABELS = {"c": "c", "phi": "phi", "c_res": "c_res", "phi_res": "phi_res",
                "dil": "dil", "t_cut": "t_cut", "kn": "kn", "ks": "ks",
                "jred": "Jred"}
