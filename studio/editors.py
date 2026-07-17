"""Input editors — modal dialogs launched from the Inputs tree (Phase 2).

Each editable ``slope_data`` category has a small ``CategoryEditor`` that builds a
dialog from the current data and writes the edited values back. Most categories
reuse the generic ``TableEditorDialog`` (a column spec + add/remove rows);
scalars use ``FormEditorDialog``. Editors only expose the commonly-edited fields
of each record and **preserve any unshown keys** (e.g. material reliability
sigmas), so a round-trip through an editor never drops data.
"""

from __future__ import annotations

import copy
import math

from PySide6.QtCore import Qt, QTimer
from PySide6.QtGui import QColor, QIcon, QPixmap
from PySide6.QtWidgets import (
    QAbstractItemView, QButtonGroup, QCheckBox, QComboBox, QDialog,
    QDialogButtonBox, QFormLayout, QGroupBox, QHBoxLayout, QLabel, QLineEdit,
    QListWidget, QListWidgetItem, QPushButton, QScrollArea, QSplitter,
    QStackedWidget, QTableWidget, QTableWidgetItem, QTabWidget, QVBoxLayout,
    QWidget,
)

# Column "usage" tags: which analysis a field applies to. Header text is colored
# to mirror the input template's header coloring (red = LEM-specific inputs,
# blue = FEM-specific inputs); seepage/reliability extend the same idea.
USAGE_COLOR = {"lem": "#c00000", "fem": "#0432ff", "seep": "#2e7d32", "rel": "#9a7d0a"}
USAGE_NAME = {"lem": "LEM only", "fem": "FEM only",
              "seep": "Seepage only", "rel": "Reliability"}
# Short labels for the show/hide toggles at the top of an editor.
USAGE_TOGGLE_LABEL = {"lem": "LEM", "fem": "FEM", "seep": "Seepage", "rel": "Reliability"}


# --------------------------------------------------------------------------- #
# Field spec + dialogs
# --------------------------------------------------------------------------- #
class Field:
    """One editable column/field. kind ∈ {'float','optfloat','int','str','choice'}.

    'optfloat' is an optional float: a blank cell reads back as None (not 0.0), so
    optional fields like a pile's H or capacities stay unset instead of becoming
    zero (which the loader would reject)."""

    _BLANK = {"float": 0.0, "optfloat": None, "int": 0, "str": "", "choice": ""}

    def __init__(self, key, header, kind="float", choices=None, default=None,
                 usage=None, applies=None):
        self.key = key
        self.header = header
        self.kind = kind
        self.choices = [str(c) for c in (choices or [])]
        # `applies`: the set of analyses this field is used in — drives the
        # show/hide column toggles; None = universal (geometry/identity, always
        # shown). `usage`: the single analysis a field is *specific* to, used for
        # the mirrored header color — only set when the field belongs to exactly
        # one analysis (shared fields like c/φ stay uncolored but still hide when
        # neither of their analyses is toggled on).
        if applies is None and usage is not None:
            applies = {usage}
        self.applies = frozenset(applies) if applies is not None else None
        self.usage = usage if usage is not None else (
            next(iter(self.applies)) if (self.applies and len(self.applies) == 1) else None)
        if default is None:
            default = self.choices[0] if kind == "choice" and self.choices else self._BLANK[kind]
        self.default = default

    def from_text(self, text):
        text = (text or "").strip()
        if self.kind == "float":
            try:
                return float(text)
            except ValueError:
                return 0.0
        if self.kind == "optfloat":
            if text == "" or text.lower() == "none":
                return None
            try:
                return float(text)
            except ValueError:
                return None
        if self.kind == "int":
            try:
                return int(float(text))
            except ValueError:
                return 0
        return text


class FormEditorDialog(QDialog):
    """Simple key/value form for scalar parameters."""

    def __init__(self, title, fields, values, parent=None):
        super().__init__(parent)
        self.setWindowTitle(title)
        self._fields = fields
        self._edits = {}
        layout = QVBoxLayout(self)
        form = QFormLayout()
        for f in fields:
            edit = QLineEdit(str(values.get(f.key, f.default)))
            self._edits[f.key] = edit
            form.addRow(f.header, edit)
        layout.addLayout(form)
        bb = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        bb.accepted.connect(self.accept)
        bb.rejected.connect(self.reject)
        layout.addWidget(bb)

    def result_values(self):
        return {f.key: f.from_text(self._edits[f.key].text()) for f in self._fields}


class _EditableTable(QWidget):
    """A table over a list of dict records with Add/Remove rows. Unshown keys are
    preserved. Reused standalone (TableEditorDialog) and per-tab (TabbedTableEditorDialog)."""

    def __init__(self, fields, rows, new_row, parent=None, swatch_state=None,
                 on_change=None, on_select=None):
        super().__init__(parent)
        self._fields = fields
        self._new_row = new_row
        # Optional live-edit hooks (used by the editor previews): on_change fires on
        # any data edit (cell text, combo, add/remove); on_select on a row-selection
        # change. Suppressed while the table populates so construction is silent.
        self._on_change = on_change
        self._on_select = on_select
        self._suppress_notify = True
        self._bases = [dict(r) for r in rows]  # keep originals to preserve extra keys
        # Optional leading display-color swatch column (Materials editor). It is a
        # *display* column, not a data field: kept as the appended LOGICAL column
        # (len(fields)) so every field/column index elsewhere — result_rows,
        # apply_usage_filter, callers indexing by field position — is unchanged, then
        # moved to VISUAL position 0 so it reads as a leading swatch. Slot-keyed via
        # the shared _MaterialColorState; committed to the style delta on OK.
        self._swatch = swatch_state
        ncols = len(fields) + (1 if swatch_state is not None else 0)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        self.table = QTableWidget(len(rows), ncols)
        self.table.setHorizontalHeaderLabels([f.header for f in fields])
        self.table.setSelectionBehavior(QAbstractItemView.SelectRows)
        for j, f in enumerate(fields):  # color usage-tagged headers (red=LEM, blue=FEM, …)
            if f.usage:
                hi = self.table.horizontalHeaderItem(j)
                if hi is not None:
                    hi.setForeground(QColor(USAGE_COLOR[f.usage]))
                    fnt = hi.font()
                    fnt.setBold(True)
                    hi.setFont(fnt)
        if swatch_state is not None:
            col = len(fields)
            hdr = QTableWidgetItem("\U0001F3A8")   # palette glyph — no data header text
            hdr.setToolTip("Display color on the Inputs plot — not a 'mat' worksheet "
                           "column. Click a swatch to change it; pick “Default” "
                           "to reset to the palette color.")
            self.table.setHorizontalHeaderItem(col, hdr)
            self.table.horizontalHeader().moveSection(col, 0)   # -> leading (visual)
        layout.addWidget(self.table)
        for i, row in enumerate(rows):
            self._set_row(i, row)
        if swatch_state is not None:
            self._rebuild_swatches()
        # Notifications are wired only AFTER the initial population, and item edits
        # go through a suppress flag, so building the table fires nothing.
        self.table.itemChanged.connect(lambda *_: self._emit_change())
        self.table.itemSelectionChanged.connect(self._emit_select)
        self._suppress_notify = False

        bar = QHBoxLayout()
        add = QPushButton("Add row")
        add.clicked.connect(self._add_row)
        rem = QPushButton("Remove selected")
        rem.clicked.connect(self._remove_rows)
        bar.addWidget(add)
        bar.addWidget(rem)
        bar.addStretch(1)
        layout.addLayout(bar)

    def _emit_change(self):
        if not self._suppress_notify and self._on_change is not None:
            self._on_change()

    def _emit_select(self):
        if not self._suppress_notify and self._on_select is not None:
            self._on_select()

    def selected_row(self):
        """Index of the currently selected row, or -1 if none."""
        return self.table.currentRow()

    def apply_usage_filter(self, enabled):
        """Show only columns whose ``applies`` set intersects ``enabled`` (a set of
        analysis tags). Universal columns (applies=None) are always shown. Hidden
        columns keep their data, so toggling never drops values on save."""
        for j, f in enumerate(self._fields):
            visible = f.applies is None or bool(f.applies & enabled)
            self.table.setColumnHidden(j, not visible)

    def _set_row(self, i, row):
        for j, f in enumerate(self._fields):
            val = row.get(f.key, f.default)
            if f.kind == "choice":
                combo = QComboBox()
                combo.addItems(f.choices)
                if str(val) in f.choices:
                    combo.setCurrentText(str(val))
                combo.currentIndexChanged.connect(lambda *_: self._emit_change())
                self.table.setCellWidget(i, j, combo)
            else:
                self.table.setItem(i, j, QTableWidgetItem("" if val is None else str(val)))

    def _rebuild_swatches(self):
        """(Re)build the leading display-color swatch button for every row, bound to
        the row's CURRENT index. Rebuilt after any add/remove so each button targets
        its live slot — colors are slot-keyed (they follow the row position, matching
        the canvas), so after a removal the remaining rows show their slot's color."""
        if self._swatch is None:
            return
        from .styles_dialog import ColorButton, MATERIAL_PALETTE
        col = len(self._fields)
        for i in range(self.table.rowCount()):
            btn = ColorButton(self._swatch.resolved_hex(i),
                              default_hex=self._swatch.default_hex(i),
                              palette=MATERIAL_PALETTE)
            btn.colorChanged.connect(lambda h, idx=i: self._swatch.set(idx, h))
            self.table.setCellWidget(i, col, btn)

    def _add_row(self):
        i = self.table.rowCount()
        self.table.insertRow(i)
        base = self._new_row()
        self._bases.append(base)
        self._set_row(i, base)
        self._rebuild_swatches()
        self._emit_change()

    def _remove_rows(self):
        for r in sorted({idx.row() for idx in self.table.selectedIndexes()}, reverse=True):
            self.table.removeRow(r)
            if r < len(self._bases):
                self._bases.pop(r)
        self._rebuild_swatches()
        self._emit_change()

    def result_rows(self):
        out = []
        for i in range(self.table.rowCount()):
            base = dict(self._bases[i]) if i < len(self._bases) else self._new_row()
            for j, f in enumerate(self._fields):
                widget = self.table.cellWidget(i, j)
                if isinstance(widget, QComboBox):
                    base[f.key] = widget.currentText()
                else:
                    item = self.table.item(i, j)
                    base[f.key] = f.from_text(item.text() if item else "")
            out.append(base)
        return out


def _help_label(text):
    lbl = QLabel(text)
    lbl.setWordWrap(True)
    return lbl


def _usage_legend(fields):
    """A colored legend for any usage-tagged columns, or None if there are none."""
    present = [u for u in ("lem", "fem", "seep", "rel")
               if any(getattr(f, "usage", None) == u for f in fields)]
    if not present:
        return None
    spans = [f'<b><span style="color:{USAGE_COLOR[u]}">&#9632; {USAGE_NAME[u]}</span></b>'
             for u in present]
    lbl = QLabel("Header colors:&nbsp;&nbsp; " + "&nbsp;&nbsp;&nbsp; ".join(spans))
    lbl.setTextFormat(Qt.RichText)
    return lbl


def _ok_cancel(dialog, layout):
    bb = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
    bb.accepted.connect(dialog.accept)
    bb.rejected.connect(dialog.reject)
    layout.addWidget(bb)


# --------------------------------------------------------------------------- #
# Editor preview drawing — the full cross-section with the object being edited
# HIGHLIGHTED, shared by the profile-line / polygon / starting-circle editors
# (Materials list-view plots generalized). One emphasis style + one dim style,
# consistent across editors: emphasis is the single deliberate color choice; "dim"
# is just reduced alpha on each feature's own (material/style-derived) color, so the
# dimmed context still reads as itself. Each drawer renders from the editor's PENDING
# rows; the debounced PreviewPane keeps the last good frame if a half-typed row makes
# a drawer raise.
# --------------------------------------------------------------------------- #
_PREVIEW_EMPH = "#ff6d00"        # emphasis color for the selected object (saturated
                                 # orange — stands out against the tab10 material
                                 # palette and the red trial surfaces alike)
_PREVIEW_EMPH_LW = 3.0           # emphasized line width
_PREVIEW_DIM_ALPHA = 0.35        # alpha applied to non-selected context geometry


def _doc_style(parent):
    """The project's resolved style delta reached via the parent window's document,
    so a preview renders with the same material colors / feature styles as the Inputs
    plot. None (defaults) in headless/round-trip use where there is no document."""
    doc = getattr(parent, "doc", None)
    return getattr(doc, "style", None) if doc is not None else None


def _finish_preview_axes(ax):
    """Match plot_inputs' final touches so the preview frames the section the same
    way: equal aspect (data-lim adjustable, so circles read round), no grid, a touch
    of headroom. Titles/legends are omitted — the preview is a focused thumbnail."""
    ax.set_aspect("equal", adjustable="datalim")
    ax.grid(False)
    y0, y1 = ax.get_ylim()
    if y1 > y0:
        pad = 0.05 * (y1 - y0)
        ax.set_ylim(y0, y1 + pad)


def _material_color(style, mat_id, fallback_idx):
    from xslope.style import material_style
    idx = mat_id if mat_id is not None else fallback_idx
    return material_style(style, idx)["color"]


def _draw_input_features(ax, slope_data, style, skip=()):
    """Draw the input *overlay* features on ``ax`` — the piezo line, distributed
    loads, tension crack, reinforcement, piles, line loads and the trial failure
    surface(s): everything ``plot_inputs`` layers OVER the base geometry. Uses the
    engine's per-feature ``plot_*`` helpers (which already take an Axes) so a preview
    shows the FULL model as context around the object being edited. Each feature is
    drawn in its own try/except so one malformed record (a half-edited preview) can't
    blank the rest of the model. ``skip`` omits named features — e.g. the circles
    editor skips ``{"circles"}`` and draws its own live-edited circles on top.

    Editor-side (not in xslope.plot) so it composes the display helpers without
    touching the engine's plotting module."""
    from xslope import plot as _p

    def _try(name, fn):
        if name in skip:
            return
        try:
            fn()
        except Exception:
            pass

    materials = slope_data.get("materials") or []
    if any(m.get("u") == "piezo" for m in materials):
        _try("piezo", lambda: _p.plot_piezo_line(ax, slope_data, style=style))
    _try("dloads", lambda: _p.plot_dloads(ax, slope_data, style=style))
    _try("tcrack", lambda: _p.plot_tcrack_surface(ax, slope_data, style=style))
    _try("reinforcement", lambda: _p.plot_reinforcement_lines(ax, slope_data, style=style))
    _try("piles", lambda: _p.plot_piles(ax, slope_data, style=style))
    _try("line_loads", lambda: _p.plot_line_loads(ax, slope_data, style=style))
    if slope_data.get("circular"):
        _try("circles", lambda: _p.plot_circles(ax, slope_data, style=style))
    elif slope_data.get("non_circ"):
        _try("non_circ", lambda: _p.plot_non_circ(ax, slope_data["non_circ"], style=style))


def _draw_profile_preview(ax, lines, selected, max_depth, slope_data, style):
    """Preview for the profile-line editor: the PENDING profile lines drawn over the
    otherwise-current model. The selected line is bold (emphasis color) with vertex
    markers; the others keep their material color, thin and dimmed. The max-depth
    base and the light surface overlays (piezo / loads / reinforcement / piles) give
    context; trial surfaces are skipped as orthogonal clutter. Zone fills are NOT
    re-derived here — rebuilding them from a half-edited line is fragile — so this is
    a line preview (see the caption)."""
    from xslope.plot import plot_max_depth
    from xslope.style import resolve_style
    rstyle = resolve_style(style)
    _draw_input_features(ax, slope_data, rstyle, skip={"circles", "non_circ"})
    if max_depth is not None:
        drawable = [{"coords": ln["coords"]} for ln in lines if ln.get("coords")]
        if drawable:
            try:
                plot_max_depth(ax, drawable, max_depth, style=rstyle)
            except Exception:
                pass
    for i, ln in enumerate(lines):
        coords = ln.get("coords") or []
        if not coords:
            continue
        xs = [c[0] for c in coords]
        ys = [c[1] for c in coords]
        color = _material_color(rstyle, ln.get("mat_id"), i)
        if i == selected:
            ax.plot(xs, ys, color=_PREVIEW_EMPH, linewidth=_PREVIEW_EMPH_LW,
                    marker="o", markersize=5, markerfacecolor=_PREVIEW_EMPH,
                    markeredgecolor="white", zorder=20)
        else:
            ax.plot(xs, ys, color=color, linewidth=1.2, alpha=_PREVIEW_DIM_ALPHA,
                    zorder=6)
    _finish_preview_axes(ax)


def _draw_polygon_preview(ax, polys, selected, slope_data, style):
    """Preview for the polygon editor: the PENDING material zones. The selected zone
    is filled (its material color, higher opacity + hatch) with a bold emphasis edge
    and vertex markers; the others are dimmed fills. The domain base and light surface
    overlays give context; trial surfaces are skipped. Rings are closed for display."""
    from xslope.plot import plot_domain_base
    from xslope.style import resolve_style
    rstyle = resolve_style(style)
    _draw_input_features(ax, slope_data, rstyle, skip={"circles", "non_circ"})
    try:
        plot_domain_base(ax, slope_data.get("domain_polygon"), style=rstyle)
    except Exception:
        pass
    for i, pg in enumerate(polys):
        coords = pg.get("coords") or []
        if len(coords) < 2:
            continue
        xs = [c[0] for c in coords] + [coords[0][0]]   # close the ring for display
        ys = [c[1] for c in coords] + [coords[0][1]]
        color = _material_color(rstyle, pg.get("mat_id"), i)
        if i == selected:
            ax.fill(xs, ys, facecolor=color, alpha=0.55, hatch="//",
                    edgecolor=_PREVIEW_EMPH, linewidth=_PREVIEW_EMPH_LW, zorder=20)
            ax.plot(xs, ys, color=_PREVIEW_EMPH, linewidth=_PREVIEW_EMPH_LW,
                    marker="o", markersize=4, markerfacecolor=_PREVIEW_EMPH,
                    markeredgecolor="white", zorder=21)
        else:
            ax.fill(xs, ys, facecolor=color, alpha=0.12, zorder=4)
            ax.plot(xs, ys, color=color, linewidth=1.0, alpha=_PREVIEW_DIM_ALPHA,
                    zorder=5)
    _finish_preview_axes(ax)


def _circle_radius_depth(c):
    """(Xo, Yo, R, Depth) from a pending circle row, resolving R/Depth from the row's
    Option EXACTLY as CirclesEditor.apply does — so the preview tracks the option live
    (a Radius/Intercept row previews from R, a Depth row from Depth)."""
    def f(key):
        try:
            return float(c.get(key, 0) or 0)
        except (TypeError, ValueError):
            return 0.0
    Xo, Yo = f("Xo"), f("Yo")
    opt = str(c.get("Option", "Depth"))
    if opt == "Radius":
        R = f("R"); depth = Yo - R
    elif opt == "Intercept":
        R = ((f("Xi") - Xo) ** 2 + (f("Yi") - Yo) ** 2) ** 0.5
        depth = Yo - R
    else:
        depth = f("Depth"); R = Yo - depth
    return Xo, Yo, R, depth


def _draw_circles_preview(ax, circles, selected, slope_data, style):
    """Preview for the starting-circles editor: the full cross-section (base geometry
    + overlays) with the PENDING circles over it. Each circle's clipped failure arc is
    drawn as the engine does; the selected one is bold (emphasis color) with a center
    marker, a radius line and a depth line, the others faint. The center/radius are
    annotation-layer artists (in_layout=False, clipped) so a center far above the
    section can't inflate the framed view — matching plot_circles."""
    import numpy as np
    from matplotlib.lines import Line2D
    from xslope.plot import plot_base_geometry
    from xslope.slice import generate_failure_surface
    from xslope.style import resolve_style
    rstyle = resolve_style(style)
    plot_base_geometry(ax, slope_data, labels=False, style=rstyle)
    _draw_input_features(ax, slope_data, rstyle, skip={"circles", "non_circ"})
    ground_surface = slope_data.get("ground_surface")
    tcrack_depth = slope_data.get("tcrack_depth", 0)

    def _annot_line(xs, ys, **kw):
        ln = Line2D(xs, ys, **kw)
        ln.set_in_layout(False)          # keep an off-section center out of the
        ln.set_clip_box(ax.bbox)         # tight-bbox layout + view autoscale
        ln.set_clip_on(True)
        ax.add_artist(ln)

    for i, c in enumerate(circles):
        Xo, Yo, R, depth = _circle_radius_depth(c)
        if R <= 0:
            continue
        emph = (i == selected)
        arc = None
        if ground_surface is not None and not getattr(ground_surface, "is_empty", False):
            try:
                ok, res = generate_failure_surface(
                    ground_surface, circular=True,
                    circle={"Xo": Xo, "Yo": Yo, "Depth": depth, "R": R},
                    tcrack_depth=tcrack_depth)
                if ok:
                    clipped = res[4]
                    from shapely.geometry import LineString
                    if not isinstance(clipped, LineString):
                        clipped = LineString(clipped)
                    arc = list(clipped.coords)
            except Exception:
                arc = None
        if arc:
            ax_ = [p[0] for p in arc]; ay_ = [p[1] for p in arc]
            if emph:
                ax.plot(ax_, ay_, color=_PREVIEW_EMPH, linewidth=_PREVIEW_EMPH_LW,
                        zorder=20)
            else:
                ax.plot(ax_, ay_, color="red", linestyle="--", linewidth=1.0,
                        alpha=_PREVIEW_DIM_ALPHA, zorder=6)
        if emph:
            # Center +, radius line (down to the circle bottom) and depth line, all
            # as annotation artists so they never blow up the framed section.
            _annot_line([Xo], [Yo], marker="+", color=_PREVIEW_EMPH, linestyle="None",
                        markersize=12, markeredgewidth=2)
            _annot_line([Xo, Xo], [Yo, depth], color=_PREVIEW_EMPH, linewidth=1.2)
            _annot_line([Xo - R, Xo + R], [depth, depth], color=_PREVIEW_EMPH,
                        linewidth=1.0, linestyle=":")
    _finish_preview_axes(ax)


class TableEditorDialog(QDialog):
    """Editable table over a list of dict records.

    ``usage_toggles`` (a list of analysis tags, e.g. ``["lem", "fem"]``) adds a row
    of checkboxes that show/hide the columns specific to each analysis, so the user
    sees only the inputs relevant to what they're doing. The toggle state persists
    per dialog. When omitted, a static color legend is shown instead.

    ``preview_draw`` (a hook ``draw(ax, rows, selected_index)``) attaches a live
    preview pane on the right behind a splitter — the highlight follows the selected
    table row. This keeps a pure-table editor a table (no list-view rewrite); only the
    editors that pass a hook (currently starting circles) grow a preview."""

    def __init__(self, title, fields, rows, new_row, parent=None, help_text=None,
                 usage_toggles=None, preview_draw=None, preview_caption=None):
        super().__init__(parent)
        self.setWindowTitle(title)
        self._title = title
        self._preview_draw = preview_draw
        self._preview = None
        self.resize(min(1200, 160 + 110 * len(fields)), 460)
        layout = QVBoxLayout(self)
        if help_text:
            layout.addWidget(_help_label(help_text))
        on_change = self._schedule_preview if preview_draw is not None else None
        on_select = self._schedule_preview if preview_draw is not None else None
        self._editable = _EditableTable(fields, rows, new_row,
                                        on_change=on_change, on_select=on_select)
        if usage_toggles:
            layout.addLayout(self._build_toggle_bar(usage_toggles))
        else:
            legend = _usage_legend(fields)
            if legend:
                layout.addWidget(legend)
        if preview_draw is not None:
            from .canvas import PreviewPane
            self._preview = PreviewPane(
                lambda ax: self._preview_draw(ax, self._editable.result_rows(),
                                              self._editable.selected_row()),
                caption=preview_caption)
            split = QSplitter(Qt.Horizontal)
            split.addWidget(self._editable)
            split.addWidget(self._preview)
            split.setStretchFactor(0, 1)
            split.setStretchFactor(1, 1)
            split.setSizes([560, 500])
            layout.addWidget(split, 1)
            self.resize(min(1280, 620 + 110 * len(fields)), 520)
        else:
            layout.addWidget(self._editable)
        _ok_cancel(self, layout)
        if usage_toggles:
            self._apply_toggles()      # set initial column visibility
        if self._preview is not None:
            self._preview.refresh_now()

    def _schedule_preview(self, *_):
        if self._preview is not None:
            self._preview.schedule()

    def _build_toggle_bar(self, tags):
        from PySide6.QtWidgets import QCheckBox
        from PySide6.QtCore import QSettings
        s = QSettings("XSlope", "XSlope Studio")
        bar = QHBoxLayout()
        bar.addWidget(QLabel("Show columns for:"))
        self._toggles = {}
        for t in tags:
            cb = QCheckBox(USAGE_TOGGLE_LABEL[t])
            default = (t != "rel")     # reliability σ columns hidden by default (niche)
            cb.setChecked(bool(s.value(f"editor_toggles/{self._title}/{t}",
                                       default, type=bool)))
            cb.setStyleSheet(f"color:{USAGE_COLOR[t]}; font-weight:bold;")
            cb.toggled.connect(self._apply_toggles)
            self._toggles[t] = cb
            bar.addWidget(cb)
        bar.addStretch(1)
        return bar

    def _apply_toggles(self):
        from PySide6.QtCore import QSettings
        enabled = {t for t, cb in self._toggles.items() if cb.isChecked()}
        self._editable.apply_usage_filter(enabled)
        s = QSettings("XSlope", "XSlope Studio")
        for t, cb in self._toggles.items():
            s.setValue(f"editor_toggles/{self._title}/{t}", cb.isChecked())

    def result_rows(self):
        return self._editable.result_rows()


class TabbedTableEditorDialog(QDialog):
    """A tabbed dialog, one editable table per tab (e.g. the two piezo lines)."""

    def __init__(self, title, tabs, parent=None, help_text=None):
        # tabs: list of (tab_title, fields, rows, new_row)
        super().__init__(parent)
        self.setWindowTitle(title)
        max_cols = max(len(fields) for _, fields, _, _ in tabs)
        self.resize(min(1200, 200 + 110 * max_cols), 480)
        layout = QVBoxLayout(self)
        if help_text:
            layout.addWidget(_help_label(help_text))
        self._tabs = QTabWidget()
        self._editables = []
        for tab_title, fields, rows, new_row in tabs:
            et = _EditableTable(fields, rows, new_row)
            self._tabs.addTab(et, tab_title)
            self._editables.append(et)
        layout.addWidget(self._tabs)
        _ok_cancel(self, layout)

    def result_rows(self, index):
        return self._editables[index].result_rows()


class _BlockListWidget(QWidget):
    """Master/detail over a list of blocks, where each block is a list of dict
    rows: a list of blocks (left) + the selected block's row table (right).
    Used for distributed loads (and reusable for other block-structured inputs)."""

    def __init__(self, blocks, fields, new_row, block_label="Load", parent=None):
        super().__init__(parent)
        self._fields = fields
        self._new_row = new_row
        self._block_label = block_label
        self._blocks = [[dict(r) for r in blk] for blk in (blocks or [])]
        self._cur = -1
        self.table = None

        body = QHBoxLayout(self)
        body.setContentsMargins(0, 0, 0, 0)
        left = QVBoxLayout()
        body.addLayout(left)
        self.list = QListWidget()
        self.list.currentRowChanged.connect(self._on_select)
        left.addWidget(self.list)
        lb = QHBoxLayout()
        b_add = QPushButton(f"Add {block_label.lower()}")
        b_add.clicked.connect(self._add_block)
        b_rem = QPushButton("Remove")
        b_rem.clicked.connect(self._remove_block)
        lb.addWidget(b_add)
        lb.addWidget(b_rem)
        left.addLayout(lb)
        self._holder = QVBoxLayout()
        body.addLayout(self._holder, 1)

        self._refresh_list()
        if self._blocks:
            self.list.setCurrentRow(0)

    def _refresh_list(self):
        self.list.blockSignals(True)
        self.list.clear()
        for i in range(len(self._blocks)):
            self.list.addItem(f"{self._block_label} {i + 1}")
        self.list.blockSignals(False)

    def _commit_current(self):
        if 0 <= self._cur < len(self._blocks) and self.table is not None:
            self._blocks[self._cur] = self.table.result_rows()

    def _load(self, idx):
        if self.table is not None:
            self.table.setParent(None)
            self.table = None
        if not (0 <= idx < len(self._blocks)):
            return
        self.table = _EditableTable(self._fields, self._blocks[idx], self._new_row)
        self._holder.addWidget(self.table)

    def _on_select(self, idx):
        self._commit_current()
        self._cur = idx
        self._load(idx)

    def _add_block(self):
        self._commit_current()
        self._blocks.append([])
        self._refresh_list()
        self.list.setCurrentRow(len(self._blocks) - 1)

    def _remove_block(self):
        idx = self.list.currentRow()
        if idx < 0:
            return
        self._blocks.pop(idx)
        self._cur = -1
        self._refresh_list()
        if self._blocks:
            self.list.setCurrentRow(min(idx, len(self._blocks) - 1))
        else:
            self._load(-1)

    def result_blocks(self):
        self._commit_current()
        return [list(b) for b in self._blocks]


# --------------------------------------------------------------------------- #
# Per-category editors
# --------------------------------------------------------------------------- #
class CategoryEditor:
    label = ""

    def build(self, slope_data, parent):
        raise NotImplementedError

    def apply(self, slope_data, dlg):
        raise NotImplementedError


class GlobalEditor(CategoryEditor):
    label = "Global parameters"
    FIELDS = [
        Field("gamma_water", "Unit weight of water"),
        Field("tcrack_depth", "Tension crack depth"),
        Field("tcrack_water", "Water in crack"),
        Field("k_seismic", "Seismic coefficient k"),
    ]

    def build(self, slope_data, parent):
        return FormEditorDialog("Global parameters", self.FIELDS, slope_data, parent)

    def apply(self, slope_data, dlg):
        slope_data.update(dlg.result_values())
        # Rebuild the derived tension-crack line (loader does the same from depth).
        from shapely.geometry import LineString
        td = slope_data.get("tcrack_depth", 0)
        gs = slope_data.get("ground_surface")
        if td and td > 0 and gs is not None and not gs.is_empty:
            slope_data["tcrack_surface"] = LineString([(x, y - td) for x, y in gs.coords])
        else:
            slope_data["tcrack_surface"] = None


def _new_material():
    # Same key set load_slope_data produces (gsat blank -> None; the rest zeroed).
    return {"name": "", "gamma": 0.0, "gamma_sat": None, "option": "mc",
            "c": 0.0, "phi": 0.0, "cp": 0.0, "r_elev": 0.0, "d": 0.0, "psi": 0.0,
            "pow_a": 0.0, "pow_b": 0.0, "pow_c": 0.0, "pow_d": 0.0,
            "hb_sci": 0.0, "hb_gsi": 0.0, "hb_mi": 0.0, "hb_d": 0.0,
            "u": "none", "ru": 0.0,
            "sigma_gamma": 0.0, "sigma_c": 0.0, "sigma_phi": 0.0, "sigma_cp": 0.0,
            "sigma_d": 0.0, "sigma_psi": 0.0, "k1": 0.0, "k2": 0.0, "alpha": 0.0,
            "unsat": "lf", "kr0": 0.0, "h0": 0.0, "vg_a": 0.0, "vg_n": 0.0,
            "E": 0.0, "nu": 0.0}


# --------------------------------------------------------------------------- #
# Materials editor — Table view (the Excel-mirroring table) + List view (a
# per-material form with strength/kr confirmation plots). Both bind to the SAME
# underlying rows, so switching mid-edit is lossless.
# --------------------------------------------------------------------------- #

# Remembered within the session (module-global); first open defaults to Table.
_LAST_MATERIALS_VIEW = "table"

# Base fields that carry a paired sigma_* (shown beside the value under Reliability).
_MAT_SIGMA = {"gamma": "sigma_gamma", "c": "sigma_c", "phi": "sigma_phi",
              "cp": "sigma_cp", "d": "sigma_d", "psi": "sigma_psi"}
# Strength option -> the option-specific fields shown in the list-view form (only
# the selected option's are visible; blank shows none, per bc0999d).
_MAT_OPTION_FIELDS = {"mc": ["c", "phi"], "cp": ["c", "cp", "r_elev"],
                      "pow": ["pow_a", "pow_b", "pow_c", "pow_d"],
                      "hb": ["hb_sci", "hb_gsi", "hb_mi", "hb_d"], "": []}
_MAT_ALL_OPTION_FIELDS = ["c", "phi", "cp", "r_elev", "pow_a", "pow_b", "pow_c",
                          "pow_d", "hb_sci", "hb_gsi", "hb_mi", "hb_d"]
# Unsaturated model -> its curve params. gard reuses vg's a/n columns (fileio.py).
_MAT_UNSAT_FIELDS = {"lf": ["kr0", "h0"], "vg": ["vg_a", "vg_n"],
                     "gard": ["vg_a", "vg_n"]}
_MAT_ALL_UNSAT_FIELDS = ["kr0", "h0", "vg_a", "vg_n"]


class _MaterialColorState:
    """Working per-material display-color overrides for the Materials editor, shared
    by both views and committed to the document's style delta on OK.

    Seeded from the incoming sparse style delta's material *color* entries; stays
    sparse — an index only appears when its color differs from the tab10 palette
    default (so "reset to default" simply drops the key). Keyed by material INDEX
    (slot), exactly as ``xslope.style`` / the canvas resolve material color
    (``str(mat_id)`` == index), so an override stays with the *slot*, not the
    material, when materials are added/removed/reordered — matching the Styles
    dialog and the Inputs plot. Hatch/alpha (owned by the Styles dialog) are not
    touched here; the OK merge carries them through untouched."""

    def __init__(self, style):
        from matplotlib.colors import to_hex
        self._over = {}          # idx -> hex string
        for k, ov in ((style or {}).get("materials") or {}).items():
            if isinstance(ov, dict) and ov.get("color"):
                try:
                    self._over[int(k)] = to_hex(ov["color"])
                except (ValueError, TypeError):
                    pass

    def default_hex(self, idx):
        from matplotlib.colors import to_hex
        from xslope.plot import get_material_color
        return to_hex(get_material_color(idx))

    def resolved_hex(self, idx):
        return self._over.get(idx) or self.default_hex(idx)

    def has_override(self, idx):
        return idx in self._over

    def set(self, idx, hex_color):
        """Record a color for ``idx``, staying sparse: a value equal to the palette
        default drops the override (that's how "Default" in the picker resets)."""
        if not hex_color or hex_color.lower() == self.default_hex(idx).lower():
            self._over.pop(idx, None)
        else:
            self._over[idx] = hex_color

    def reset(self, idx):
        self._over.pop(idx, None)


def _material_swatch(idx, color=None):
    """A color swatch QIcon for material ``idx``. When ``color`` (a hex / mpl color)
    is given it is drawn directly — the editor passes the pending override or the
    resolved default so a list entry tracks in-dialog color edits. Otherwise it falls
    back to the SAME palette the canvas zones use (``xslope.style.material_style`` ->
    tab10 fallback), so a list entry reads as the same color as its zone on the
    Inputs plot."""
    from matplotlib.colors import to_rgba
    if color is None:
        from xslope.style import material_style, resolve_style
        color = material_style(resolve_style(None), idx)["color"]
    r, g, b, a = to_rgba(color)
    pm = QPixmap(16, 16)
    pm.fill(Qt.transparent)
    from PySide6.QtGui import QPainter, QPen
    p = QPainter(pm)
    p.setRenderHint(QPainter.Antialiasing, True)
    p.setBrush(QColor.fromRgbF(r, g, b, a))
    p.setPen(QPen(QColor("#666666")))
    p.drawRoundedRect(1, 1, 13, 13, 3, 3)
    p.end()
    return QIcon(pm)


class _MaterialListView(QWidget):
    """List view for the materials editor: three splitter panes —
    [materials list | property form | confirmation plots]. Binds to a list of
    material dicts; edits write through to those dicts immediately, so
    ``result_rows`` (and a view switch) never lose an in-progress edit. All 36
    loader keys have a widget (hidden ones keep their value), so a round-trip
    preserves every field exactly as the table view does."""

    _STRENGTH_KEYS = _MAT_ALL_OPTION_FIELDS + ["d", "psi", "E", "nu"]

    def __init__(self, fields, rows, new_row, reliability_on, parent=None,
                 color_state=None):
        super().__init__(parent)
        self._field_by_key = {f.key: f for f in fields}
        self._new_row = new_row
        self._rows = [dict(r) for r in rows]
        self._reliability_on = bool(reliability_on)
        self._color = color_state if color_state is not None else _MaterialColorState({})
        self._loading = False     # suppress color write-through while populating a row
        self._cur = -1
        self._edits = {}          # key -> QLineEdit / QComboBox
        self._cell_widgets = {}   # key -> the labeled cell QWidget (for show/hide)
        self._sigma_widgets = {}  # base key -> (sigma label, sigma edit)

        self._plot_timer = QTimer(self)
        self._plot_timer.setSingleShot(True)
        self._plot_timer.setInterval(160)     # debounce so plots don't lag typing
        self._plot_timer.timeout.connect(self._refresh_plots)

        splitter = QSplitter(Qt.Horizontal)
        splitter.addWidget(self._build_list_pane())
        splitter.addWidget(self._build_form_pane())
        splitter.addWidget(self._build_plots_pane())
        splitter.setStretchFactor(0, 0)
        splitter.setStretchFactor(1, 1)
        splitter.setStretchFactor(2, 1)
        splitter.setSizes([200, 470, 530])

        lay = QVBoxLayout(self)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.addWidget(splitter)

        self._refresh_list()
        self._update_sigma_visibility()
        if self._rows:
            self.list.setCurrentRow(0)
        else:
            self._load(-1)

    # --- panes -----------------------------------------------------------
    def _build_list_pane(self):
        w = QWidget()
        v = QVBoxLayout(w)
        v.setContentsMargins(0, 0, 0, 0)
        self.list = QListWidget()
        self.list.currentRowChanged.connect(self._on_select)
        v.addWidget(self.list, 1)
        bar = QHBoxLayout()
        add = QPushButton("Add")
        add.clicked.connect(self._add)
        rem = QPushButton("Remove")
        rem.clicked.connect(self._remove)
        bar.addWidget(add)
        bar.addWidget(rem)
        v.addLayout(bar)
        return w

    def _make_edit(self, key):
        f = self._field_by_key[key]
        if f.kind == "choice":
            w = QComboBox()
            w.addItems(f.choices)
        else:
            w = QLineEdit()
            w.editingFinished.connect(self._on_edit)
        self._edits[key] = w
        return w

    def _cell(self, label, key, sigma=False, label_w=80):
        """A [label | value | (± σ | sigma value)] cell as a QWidget, registered by
        key for show/hide + load/commit. Widths are kept tight so even the widest
        row (a d|ψ pair, both carrying σ) fits the form pane without a horizontal
        scrollbar."""
        cell = QWidget()
        h = QHBoxLayout(cell)
        h.setContentsMargins(0, 2, 0, 2)
        h.setSpacing(4)
        lab = QLabel(label)
        lab.setMinimumWidth(label_w)
        h.addWidget(lab)
        edit = self._make_edit(key)
        edit.setMinimumWidth(44)
        h.addWidget(edit, 1)
        if sigma and key in _MAT_SIGMA:
            skey = _MAT_SIGMA[key]
            slab = QLabel("± σ")
            slab.setStyleSheet(f"color:{USAGE_COLOR['rel']};")
            sedit = self._make_edit(skey)
            sedit.setFixedWidth(48)
            h.addWidget(slab)
            h.addWidget(sedit)
            self._sigma_widgets[key] = (slab, sedit)
        self._cell_widgets[key] = cell
        return cell

    @staticmethod
    def _pair(cell_a, cell_b):
        w = QWidget()
        h = QHBoxLayout(w)
        h.setContentsMargins(0, 0, 0, 0)
        h.addWidget(cell_a, 1)
        h.addWidget(cell_b, 1)
        return w

    def _build_color_row(self):
        """The Identity display-color control: a swatch button (same picker as the
        Styles dialog, seeded with the resolved override-or-default) beside a small
        Reset button that drops the override so the palette default shows through."""
        from .styles_dialog import ColorButton, MATERIAL_PALETTE
        cell = QWidget()
        h = QHBoxLayout(cell)
        h.setContentsMargins(0, 2, 0, 2)
        h.setSpacing(4)
        lab = QLabel("Display color")
        lab.setMinimumWidth(80)
        lab.setToolTip("Color of this material's zone on the Inputs plot. Stored as a "
                       "style override, not a 'mat' property.")
        h.addWidget(lab)
        self._color_btn = ColorButton("#000000", palette=MATERIAL_PALETTE)
        self._color_btn.colorChanged.connect(self._on_color_changed)
        h.addWidget(self._color_btn)
        reset = QPushButton("Reset")
        reset.setToolTip("Remove the override — show the default palette color.")
        reset.clicked.connect(self._on_color_reset)
        h.addWidget(reset)
        h.addStretch(1)
        return cell

    def _on_color_changed(self, hexc):
        if self._loading or not (0 <= self._cur < len(self._rows)):
            return
        self._color.set(self._cur, hexc)
        self._refresh_list_item(self._cur)   # track the change in the list swatch

    def _on_color_reset(self):
        if not (0 <= self._cur < len(self._rows)):
            return
        self._color.reset(self._cur)
        self._loading = True
        self._color_btn.blockSignals(True)
        self._color_btn.set_hex(self._color.default_hex(self._cur))
        self._color_btn.blockSignals(False)
        self._loading = False
        self._refresh_list_item(self._cur)

    def _build_form_pane(self):
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        body = QWidget()
        v = QVBoxLayout(body)

        # Identity — name plus the display color (swatch button + reset), the same
        # override the Styles dialog edits and the Inputs plot resolves.
        g = QGroupBox("Identity")
        gv = QVBoxLayout(g)
        gv.addWidget(self._cell("Name", "name"))
        gv.addWidget(self._build_color_row())
        v.addWidget(g)

        # Unit weights (γ | γsat side by side; γ carries its σ)
        g = QGroupBox("Unit weights")
        gv = QVBoxLayout(g)
        gv.addWidget(self._pair(self._cell("γ", "gamma", sigma=True, label_w=28),
                                self._cell("γ_sat", "gamma_sat", label_w=44)))
        v.addWidget(g)

        # Strength: option combo, then only the selected option's fields, then
        # dilation d/ψ and elastic E/ν.
        g = QGroupBox("Strength")
        gv = QVBoxLayout(g)
        opt_cell = QWidget()
        oh = QHBoxLayout(opt_cell)
        oh.setContentsMargins(0, 2, 0, 2)
        olab = QLabel("Model (option)")
        olab.setMinimumWidth(96)
        oh.addWidget(olab)
        self._opt_combo = self._make_edit("option")
        self._opt_combo.currentIndexChanged.connect(self._on_option_changed)
        oh.addWidget(self._opt_combo, 1)
        gv.addWidget(opt_cell)
        for key in _MAT_ALL_OPTION_FIELDS:
            gv.addWidget(self._cell(self._label_for(key), key,
                                    sigma=key in _MAT_SIGMA))
        gv.addSpacing(4)
        gv.addWidget(self._pair(self._cell("d (dil.)", "d", sigma=True, label_w=52),
                                self._cell("ψ", "psi", sigma=True, label_w=18)))
        gv.addWidget(self._pair(self._cell("E", "E", label_w=18),
                                self._cell("ν", "nu", label_w=18)))
        v.addWidget(g)

        # Pore pressure (ru only when u = 'ru')
        g = QGroupBox("Pore pressure")
        gv = QVBoxLayout(g)
        u_cell = QWidget()
        uh = QHBoxLayout(u_cell)
        uh.setContentsMargins(0, 2, 0, 2)
        ulab = QLabel("Model (u)")
        ulab.setMinimumWidth(96)
        uh.addWidget(ulab)
        self._u_combo = self._make_edit("u")
        self._u_combo.currentIndexChanged.connect(self._on_u_changed)
        uh.addWidget(self._u_combo, 1)
        gv.addWidget(u_cell)
        gv.addWidget(self._cell("ru", "ru"))
        v.addWidget(g)

        # Conductivity: k1/k2/alpha, then unsat model + its curve params.
        g = QGroupBox("Conductivity")
        gv = QVBoxLayout(g)
        gv.addWidget(self._pair(self._cell("k1", "k1", label_w=22),
                                self._cell("k2", "k2", label_w=22)))
        gv.addWidget(self._cell("alpha", "alpha"))
        us_cell = QWidget()
        ush = QHBoxLayout(us_cell)
        ush.setContentsMargins(0, 2, 0, 2)
        uslab = QLabel("Unsat model")
        uslab.setMinimumWidth(96)
        ush.addWidget(uslab)
        self._unsat_combo = self._make_edit("unsat")
        self._unsat_combo.currentIndexChanged.connect(self._on_unsat_changed)
        ush.addWidget(self._unsat_combo, 1)
        gv.addWidget(us_cell)
        for key in _MAT_ALL_UNSAT_FIELDS:
            gv.addWidget(self._cell(self._label_for(key), key))
        v.addWidget(g)

        v.addStretch(1)
        scroll.setWidget(body)
        return scroll

    def _build_plots_pane(self):
        from .canvas import MplCanvas
        pane = QSplitter(Qt.Vertical)
        self._strength_canvas = MplCanvas()
        self._kr_canvas = MplCanvas()
        pane.addWidget(self._strength_canvas)
        pane.addWidget(self._kr_canvas)
        pane.setSizes([320, 320])
        return pane

    # List view shows friendly symbols where the table mirrors the terse 'mat'
    # sheet headers (Norm: "use φ in the list view"). Only display labels — the
    # underlying keys/headers are untouched.
    _FRIENDLY = {"f": "φ", "psi": "ψ", "s(f)": "σ(φ)", "s(g)": "σ(γ)",
                 "s(c)": "σ(c)", "s(c/p)": "σ(c/p)", "s(d)": "σ(d)",
                 "s(psi)": "σ(ψ)"}

    def _label_for(self, key):
        f = self._field_by_key.get(key)
        header = f.header if f is not None else key
        return self._FRIENDLY.get(header, header)

    # --- list ------------------------------------------------------------
    def _item_text(self, i):
        name = self._rows[i].get("name") or ""
        return f"{i + 1}.  {name}" if name else f"{i + 1}."

    def _refresh_list(self):
        self.list.blockSignals(True)
        self.list.clear()
        for i in range(len(self._rows)):
            item = QListWidgetItem(_material_swatch(i, self._color.resolved_hex(i)),
                                   self._item_text(i))
            self.list.addItem(item)
        self.list.blockSignals(False)

    def _refresh_list_item(self, i):
        if 0 <= i < self.list.count():
            it = self.list.item(i)
            it.setText(self._item_text(i))
            it.setIcon(_material_swatch(i, self._color.resolved_hex(i)))

    # --- load / commit ---------------------------------------------------
    def _load(self, idx):
        self._cur = -1                       # suppress write-through while populating
        for key, w in self._edits.items():
            val = self._rows[idx].get(key) if (0 <= idx < len(self._rows)) else None
            if isinstance(w, QComboBox):
                w.blockSignals(True)
                txt = "" if val is None else str(val)
                j = w.findText(txt)
                w.setCurrentIndex(j if j >= 0 else 0)
                w.blockSignals(False)
            else:
                w.blockSignals(True)
                w.setText("" if val is None else str(val))
                w.blockSignals(False)
        ok = 0 <= idx < len(self._rows)
        self.setEnabled(True)
        self._loading = True
        self._color_btn.blockSignals(True)
        if ok:
            self._color_btn.set_default(self._color.default_hex(idx))
            self._color_btn.set_hex(self._color.resolved_hex(idx))
        self._color_btn.setEnabled(ok)
        self._color_btn.blockSignals(False)
        self._loading = False
        if ok:
            self._cur = idx
            self._update_option_visibility()
            self._update_u_visibility()
            self._update_unsat_visibility()
            self._update_sigma_visibility()
            self._refresh_plots()
        else:
            self._clear_plots()

    def _commit(self):
        if not (0 <= self._cur < len(self._rows)):
            return
        row = self._rows[self._cur]
        for key, w in self._edits.items():
            f = self._field_by_key[key]
            if isinstance(w, QComboBox):
                row[key] = w.currentText()
            else:
                row[key] = f.from_text(w.text())

    # --- visibility ------------------------------------------------------
    def _update_option_visibility(self):
        opt = self._opt_combo.currentText().strip().lower()
        shown = set(_MAT_OPTION_FIELDS.get(opt, []))
        for key in _MAT_ALL_OPTION_FIELDS:
            self._cell_widgets[key].setVisible(key in shown)

    def _update_u_visibility(self):
        self._cell_widgets["ru"].setVisible(
            self._u_combo.currentText().strip().lower() == "ru")

    def _update_unsat_visibility(self):
        um = self._unsat_combo.currentText().strip().lower()
        shown = set(_MAT_UNSAT_FIELDS.get(um, []))
        for key in _MAT_ALL_UNSAT_FIELDS:
            self._cell_widgets[key].setVisible(key in shown)

    def _update_sigma_visibility(self):
        for _key, (slab, sedit) in self._sigma_widgets.items():
            slab.setVisible(self._reliability_on)
            sedit.setVisible(self._reliability_on)

    def set_reliability(self, on):
        self._reliability_on = bool(on)
        self._update_sigma_visibility()

    # --- events ----------------------------------------------------------
    def _on_select(self, idx):
        self._commit()
        self._load(idx)

    def _on_edit(self):
        self._commit()
        self._refresh_list_item(self._cur)   # name edits update the list label
        self._plot_timer.start()

    def _on_option_changed(self):
        self._commit()
        self._update_option_visibility()
        self._plot_timer.start()

    def _on_u_changed(self):
        self._commit()
        self._update_u_visibility()          # u affects neither plot; no refresh

    def _on_unsat_changed(self):
        self._commit()
        self._update_unsat_visibility()
        self._plot_timer.start()             # kr plot depends on the unsat model

    def _add(self):
        self._commit()
        self._rows.append(self._new_row())
        self._cur = -1
        self._refresh_list()
        self.list.setCurrentRow(len(self._rows) - 1)

    def _remove(self):
        idx = self.list.currentRow()
        if not (0 <= idx < len(self._rows)):
            return
        self._rows.pop(idx)
        self._cur = -1
        self._refresh_list()
        if self._rows:
            self.list.setCurrentRow(min(idx, len(self._rows) - 1))
        else:
            self._load(-1)

    # --- plots -----------------------------------------------------------
    def _refresh_plots(self):
        from xslope.plot import plot_material_strength, plot_material_kr
        if not (0 <= self._cur < len(self._rows)):
            self._clear_plots()
            return
        mat = dict(self._rows[self._cur])    # snapshot for the deferred draw
        self._strength_canvas.render_axes(lambda ax: plot_material_strength(ax, mat))
        self._kr_canvas.render_axes(lambda ax: plot_material_kr(ax, mat))

    def _clear_plots(self):
        self._strength_canvas.clear()
        self._kr_canvas.clear()

    # --- external API (used by MaterialsDialog on view switch) -----------
    def set_rows(self, rows, reliability_on):
        self._reliability_on = bool(reliability_on)
        self._rows = [dict(r) for r in rows]
        self._cur = -1
        self._refresh_list()
        if self._rows:
            self.list.setCurrentRow(0)
        else:
            self._load(-1)

    def result_rows(self):
        self._commit()
        return [dict(r) for r in self._rows]


class MaterialsDialog(QDialog):
    """The materials editor dialog: a segmented Table/List view toggle over a
    QStackedWidget. Both views bind to ``self._rows`` (the single source of truth):
    switching harvests the active view's edits into it and rebuilds the target, so
    a switch mid-edit is lossless. The Reliability toggle drives σ columns (table)
    and σ fields (list) in both views."""

    def __init__(self, title, fields, rows, new_row, parent=None, help_text=None,
                 usage_toggles=None, style=None, doc=None):
        super().__init__(parent)
        self.setWindowTitle(title)
        self._title = title
        self._fields = fields
        self._new_row = new_row
        self._rows = [dict(r) for r in rows]
        self._mode = None
        self._table = None
        self._list_view = None
        # Display-color overrides live in the project style delta (NOT the material
        # dicts). Seed a working state from the incoming delta, shared by both views;
        # result_style() folds it back into the delta on OK. doc is the document to
        # commit to (None in headless/round-trip use — colors just aren't committed).
        self._doc = doc
        self._orig_style = copy.deepcopy(style or {})
        self._color = _MaterialColorState(self._orig_style)
        self.resize(1180, 640)

        layout = QVBoxLayout(self)
        if help_text:
            layout.addWidget(_help_label(help_text))

        top = QHBoxLayout()
        seg_group = QButtonGroup(self)
        seg_group.setExclusive(True)
        self._seg = {}
        for mode, lbl in (("table", "Table view"), ("list", "List view")):
            b = QPushButton(lbl)
            b.setCheckable(True)
            b.clicked.connect(lambda _checked=False, m=mode: self._set_mode(m))
            seg_group.addButton(b)
            self._seg[mode] = b
            top.addWidget(b)
        top.addSpacing(24)
        top.addWidget(QLabel("Show columns for:"))
        from PySide6.QtCore import QSettings
        s = QSettings("XSlope", "XSlope Studio")
        self._toggles = {}
        for t in (usage_toggles or []):
            cb = QCheckBox(USAGE_TOGGLE_LABEL[t])
            default = (t != "rel")
            cb.setChecked(bool(s.value(f"editor_toggles/{self._title}/{t}",
                                       default, type=bool)))
            cb.setStyleSheet(f"color:{USAGE_COLOR[t]}; font-weight:bold;")
            cb.toggled.connect(self._on_toggle)
            self._toggles[t] = cb
            top.addWidget(cb)
        top.addStretch(1)
        layout.addLayout(top)

        self._stack = QStackedWidget()
        self._table_holder = QWidget()
        self._table_lay = QVBoxLayout(self._table_holder)
        self._table_lay.setContentsMargins(0, 0, 0, 0)
        self._list_holder = QWidget()
        self._list_lay = QVBoxLayout(self._list_holder)
        self._list_lay.setContentsMargins(0, 0, 0, 0)
        self._stack.addWidget(self._table_holder)     # index 0
        self._stack.addWidget(self._list_holder)      # index 1
        layout.addWidget(self._stack, 1)

        _ok_cancel(self, layout)

        initial = _LAST_MATERIALS_VIEW if _LAST_MATERIALS_VIEW in ("table", "list") else "table"
        self._set_mode(initial)

    # --- toggles ---------------------------------------------------------
    def _reliability_on(self):
        cb = self._toggles.get("rel")
        return bool(cb.isChecked()) if cb is not None else False

    def _enabled_usage(self):
        return {t for t, cb in self._toggles.items() if cb.isChecked()}

    def _on_toggle(self):
        from PySide6.QtCore import QSettings
        s = QSettings("XSlope", "XSlope Studio")
        for t, cb in self._toggles.items():
            s.setValue(f"editor_toggles/{self._title}/{t}", cb.isChecked())
        if self._mode == "table" and self._table is not None:
            self._table.apply_usage_filter(self._enabled_usage())
        if self._mode == "list" and self._list_view is not None:
            self._list_view.set_reliability(self._reliability_on())

    # --- view switching --------------------------------------------------
    def _harvest(self):
        if self._mode == "table" and self._table is not None:
            self._rows = self._table.result_rows()
        elif self._mode == "list" and self._list_view is not None:
            self._rows = self._list_view.result_rows()

    def _build_table(self):
        if self._table is not None:
            self._table.setParent(None)
            self._table.deleteLater()
        self._table = _EditableTable(self._fields, self._rows, self._new_row,
                                     swatch_state=self._color)
        self._table_lay.addWidget(self._table)
        self._table.apply_usage_filter(self._enabled_usage())

    def _ensure_list(self):
        if self._list_view is None:
            self._list_view = _MaterialListView(self._fields, self._rows,
                                                self._new_row, self._reliability_on(),
                                                color_state=self._color)
            self._list_lay.addWidget(self._list_view)
        else:
            self._list_view.set_rows(self._rows, self._reliability_on())

    def _set_mode(self, mode):
        if mode not in ("table", "list"):
            return
        if mode == self._mode:
            self._seg[mode].setChecked(True)
            return
        self._harvest()                      # pull current edits into self._rows
        if mode == "table":
            self._build_table()
            self._stack.setCurrentIndex(0)
        else:
            self._ensure_list()
            self._stack.setCurrentIndex(1)
        self._mode = mode
        self._seg[mode].setChecked(True)
        global _LAST_MATERIALS_VIEW
        _LAST_MATERIALS_VIEW = mode

    def set_view_mode(self, mode):
        """Programmatic view switch (used by the round-trip guard)."""
        self._set_mode(mode)

    def result_rows(self):
        self._harvest()
        return self._rows

    def result_style(self):
        """The project style delta with material *display colors* reconciled from the
        editor's working overrides. Sparse: an override is written only where the
        color differs from the palette default; hatch/alpha (owned by the Styles
        dialog) are carried through untouched. Slot-keyed by material index, so
        entries for slots past the current material count are dropped. Applying with
        no color change returns a delta equal to the one passed in."""
        from matplotlib.colors import to_hex
        n = len(self.result_rows())
        out = copy.deepcopy(self._orig_style)
        mats = dict(out.get("materials") or {})
        new_mats = {}
        for idx in range(n):
            key = str(idx)
            entry = dict(mats.get(key) or {})
            if self._color.has_override(idx):
                new_hex = self._color.resolved_hex(idx)
                old = entry.get("color")
                # Keep the original color string when it's the same color (avoids
                # churning e.g. "red" -> "#ff0000" on an untouched override).
                if old is None or to_hex(old).lower() != new_hex.lower():
                    entry["color"] = new_hex
            else:
                entry.pop("color", None)
            if entry:
                new_mats[key] = entry
        if new_mats:
            out["materials"] = new_mats
        else:
            out.pop("materials", None)
        return out


class MaterialsEditor(CategoryEditor):
    label = "Materials"
    # Columns mirror the 'mat' worksheet in order: name, g, gsat, option, c, f,
    # c/p, r-elev, d, psi, pow_a..pow_d, hb_sci/hb_gsi/hb_mi/hb_d, u, ru, s(g),
    # s(c), s(f), s(c/p), s(d), s(psi), k1, k2, alpha, unsat, kr0, h0, vg_a, vg_n,
    # E, n.
    # `applies` tags mirror the template's analysis usage (input_template.md):
    # the strength block (g, gsat, option, c, f, c/p, r-elev, the pow_* power-curve
    # and hb_* Hoek-Brown envelope parameters, u, ru) is shared by LEM+FEM; d/psi
    # are rapid-drawdown (LEM); s(...) are reliability; k1..vg_n seepage; E/n FEM.
    # gsat is optional (blank -> fall back to g), so it reads back as None when
    # left empty rather than 0.0. The alternate-envelope columns (pow_*/hb_*) are
    # always shown; option-driven show/hide grouping is a later UX pass.
    LF = {"lem", "fem"}
    FIELDS = [
        Field("name", "name", "str"),
        Field("gamma", "g", applies=LF),
        Field("gamma_sat", "gsat", "optfloat", applies=LF),
        # A BLANK option is valid for seep-only material rows (the loader keeps ''
        # via _choice; document._blank_material produces it for DXF imports). Offer
        # it as an empty combo entry so the editor round-trips it instead of
        # normalizing blank -> 'mc'. Kept last so the default (first choice) stays 'mc'.
        Field("option", "option", "choice", choices=["mc", "cp", "pow", "hb", ""], applies=LF),
        Field("c", "c", applies=LF), Field("phi", "f", applies=LF),
        Field("cp", "c/p", applies=LF), Field("r_elev", "r-elev", applies=LF),
        Field("d", "d", usage="lem"), Field("psi", "psi", usage="lem"),
        Field("pow_a", "pow_a", applies=LF), Field("pow_b", "pow_b", applies=LF),
        Field("pow_c", "pow_c", applies=LF), Field("pow_d", "pow_d", applies=LF),
        Field("hb_sci", "hb_sci", applies=LF), Field("hb_gsi", "hb_gsi", applies=LF),
        Field("hb_mi", "hb_mi", applies=LF), Field("hb_d", "hb_d", applies=LF),
        Field("u", "u", "choice", choices=["none", "piezo", "seep", "ru"], applies=LF),
        Field("ru", "ru", applies=LF),
        Field("sigma_gamma", "s(g)", usage="rel"), Field("sigma_c", "s(c)", usage="rel"),
        Field("sigma_phi", "s(f)", usage="rel"), Field("sigma_cp", "s(c/p)", usage="rel"),
        Field("sigma_d", "s(d)", usage="rel"), Field("sigma_psi", "s(psi)", usage="rel"),
        Field("k1", "k1", usage="seep"), Field("k2", "k2", usage="seep"),
        Field("alpha", "alpha", usage="seep"),
        # Unsaturated model: lf (linear front -> kr0/h0), vg (van Genuchten) or gard
        # (Gardner); vg/gard share the vg_a/vg_n curve pair.
        Field("unsat", "unsat", "choice", choices=["lf", "vg", "gard"], usage="seep"),
        Field("kr0", "kr0", usage="seep"), Field("h0", "h0", usage="seep"),
        Field("vg_a", "vg_a", usage="seep"), Field("vg_n", "vg_n", usage="seep"),
        Field("E", "E", usage="fem"), Field("nu", "n", usage="fem"),
    ]

    def build(self, slope_data, parent):
        # The display-color swatch edits the project's style delta, which lives on the
        # document (reached via the parent window); headless/round-trip callers pass a
        # bare parent, so both are optional and colors simply aren't committed then.
        doc = getattr(parent, "doc", None)
        style = getattr(doc, "style", None) if doc is not None else None
        return MaterialsDialog(
            "Materials", self.FIELDS, slope_data.get("materials", []), _new_material, parent,
            help_text="Table view mirrors the 'mat' worksheet (row order = Mat ID "
                      "order). List view edits one material at a time as a form with "
                      "strength- and conductivity-model plots that confirm the "
                      "selected options. Both views edit the same rows, so switching "
                      "is lossless. In list view only 'Reliability' applies — it shows "
                      "each value's σ; the other toggles hide table columns. The color "
                      "swatch sets the material's display color on the Inputs plot.",
            usage_toggles=["lem", "seep", "fem", "rel"], style=style, doc=doc)

    def apply(self, slope_data, dlg):
        slope_data["materials"] = dlg.result_rows()
        # Display colors are style deltas, not material data — commit them via the
        # document's style store (same path as the Styles dialog), only when changed,
        # so the sparse delta stays sparse and an unchanged apply is a no-op.
        doc = getattr(dlg, "_doc", None)
        if doc is not None:
            new_style = dlg.result_style()
            if new_style != doc.style:
                doc.set_style(new_style)


def _new_circle():
    return {"Xo": 0.0, "Yo": 0.0, "Option": "Depth", "Depth": 0.0,
            "Xi": 0.0, "Yi": 0.0, "R": 0.0}


class CirclesEditor(CategoryEditor):
    label = "Circles"
    # Mirror the 'circles' worksheet: Xo, Yo, Option, Depth, Xi, Yi, R.
    FIELDS = [
        Field("Xo", "Xo"), Field("Yo", "Yo"),
        Field("Option", "Option", "choice", choices=["Depth", "Radius", "Intercept"]),
        Field("Depth", "Depth"), Field("Xi", "Xi"), Field("Yi", "Yi"), Field("R", "R"),
    ]

    def build(self, slope_data, parent):
        style = _doc_style(parent)

        def preview(ax, rows, selected):
            _draw_circles_preview(ax, rows, selected, slope_data, style)

        return TableEditorDialog(
            "Circles", self.FIELDS, slope_data.get("circles", []), _new_circle, parent,
            help_text="Option sets how each circle's size is defined (only the matching "
                      "field is used):\n"
                      "  • Depth — elevation of the circle bottom at the center (R = Yo − Depth)\n"
                      "  • Radius — the circle radius R directly (Depth = Yo − R)\n"
                      "  • Intercept — a point (Xi, Yi) the circle passes through",
            preview_draw=preview,
            preview_caption="Preview shows the starting circles on the section "
                            "(selected circle bold with center, radius and depth "
                            "lines; others faint).")

    def apply(self, slope_data, dlg):
        rows = dlg.result_rows()
        for c in rows:
            opt = str(c.get("Option", "Depth"))
            if opt == "Radius":
                c["Depth"] = c["Yo"] - c["R"]
            elif opt == "Intercept":
                c["R"] = ((c["Xi"] - c["Xo"]) ** 2 + (c["Yi"] - c["Yo"]) ** 2) ** 0.5
                c["Depth"] = c["Yo"] - c["R"]
            else:  # Depth
                c["R"] = c["Yo"] - c["Depth"]
        slope_data["circles"] = rows
        slope_data["circular"] = len(rows) > 0


def _new_ncpt():
    return {"X": 0.0, "Y": 0.0, "Movement": "Free"}


class NonCircEditor(CategoryEditor):
    label = "Non-circular surface"
    FIELDS = [Field("X", "X"), Field("Y", "Y"),
              Field("Movement", "Movement", "choice", choices=["Free", "Horiz", "Fixed"])]

    def build(self, slope_data, parent):
        return TableEditorDialog(
            "Non-circular surface", self.FIELDS, slope_data.get("non_circ", []), _new_ncpt, parent,
            help_text="Points ordered left→right. Entry/exit points use Movement='Free'.")

    def apply(self, slope_data, dlg):
        slope_data["non_circ"] = dlg.result_rows()


def _new_pt():
    return {"x": 0.0, "y": 0.0}


class PiezoEditor(CategoryEditor):
    label = "Piezometric lines"
    FIELDS = [Field("x", "x"), Field("y", "y")]

    def _rows(self, slope_data, key):
        return [{"x": x, "y": y} for (x, y) in (slope_data.get(key) or [])]

    def build(self, slope_data, parent):
        return TabbedTableEditorDialog(
            "Piezometric lines",
            [("Line 1", self.FIELDS, self._rows(slope_data, "piezo_line"), _new_pt),
             ("Line 2 (rapid drawdown)", self.FIELDS, self._rows(slope_data, "piezo_line2"), _new_pt)],
            parent,
            help_text="Points ordered left→right. Line 2 is only used for rapid-drawdown "
                      "analysis (the drawdown / second water table).")

    def apply(self, slope_data, dlg):
        slope_data["piezo_line"] = [(r["x"], r["y"]) for r in dlg.result_rows(0)]
        slope_data["piezo_line2"] = [(r["x"], r["y"]) for r in dlg.result_rows(1)]


# --- distributed loads (two sets; each a list of point blocks) -------------- #
DLOAD_FIELDS = [Field("X", "X"), Field("Y", "Y"), Field("Normal", "Normal")]


def _new_dload_pt():
    return {"X": 0.0, "Y": 0.0, "Normal": 0.0}


def _apply_set_selection(widgets, tabs, select):
    """Activate the tab and list row encoded in ``select = (set_idx, row)`` — shared
    by the dload and seep-BC editors (both two-set tabbed master/detail dialogs), so
    a canvas double-click jumps straight to the picked item. No-op for other
    `select` shapes (e.g. None, or a plain index)."""
    if not (isinstance(select, tuple) and len(select) == 2):
        return
    s, row = select
    if s not in (0, 1):
        return
    tabs.setCurrentIndex(s)
    w = widgets[s]
    if row is not None and 0 <= row < w.list.count():
        w.list.setCurrentRow(row)


class DloadsEditor(CategoryEditor):
    label = "Distributed loads"

    def build(self, slope_data, parent, select=None):
        dlg = QDialog(parent)
        dlg.setWindowTitle("Distributed loads")
        dlg.resize(640, 520)
        layout = QVBoxLayout(dlg)
        layout.addWidget(_help_label(
            "Each load is a left→right series of points (≥2). Select a load to edit its "
            "points. Set 2 is the second rapid-drawdown stage."))
        tabs = QTabWidget()
        w1 = _BlockListWidget(slope_data.get("dloads"), DLOAD_FIELDS, _new_dload_pt, "Load")
        w2 = _BlockListWidget(slope_data.get("dloads2"), DLOAD_FIELDS, _new_dload_pt, "Load")
        tabs.addTab(w1, "Set 1")
        tabs.addTab(w2, "Set 2 (rapid drawdown)")
        layout.addWidget(tabs)
        _ok_cancel(dlg, layout)
        dlg._sets = (w1, w2)
        _apply_set_selection((w1, w2), tabs, select)
        return dlg

    def apply(self, slope_data, dlg):
        w1, w2 = dlg._sets
        slope_data["dloads"] = w1.result_blocks()
        slope_data["dloads2"] = w2.result_blocks()


# --- seep BC (two sets; each: a list of specified-head BCs + an exit face) --- #
XY_FIELDS = [Field("x", "x"), Field("y", "y")]


class _SeepBcSetWidget(QWidget):
    """One seepage BC set as a single master/detail: the left list holds each
    specified-head boundary, each specified-flux boundary, AND the exit face; the
    right side shows one full-height point table plus a value field — "Head value:"
    for heads, "Flux value:" for fluxes, both hidden for the exit face.

    List order is heads … fluxes … exit face; picking.py mirrors this so a canvas
    double-click jumps to the right row."""

    def __init__(self, bc, parent=None):
        super().__init__(parent)
        bc = bc or {}
        self._heads = [{"head": h.get("head", 0.0),
                        "coords": [tuple(c) for c in h.get("coords", [])]}
                       for h in (bc.get("specified_heads") or [])]
        self._fluxes = [{"flux": f.get("flux", 0.0),
                         "coords": [tuple(c) for c in f.get("coords", [])]}
                        for f in (bc.get("specified_fluxes") or [])]
        self._exit = [tuple(c) for c in (bc.get("exit_face") or [])]
        self._cur = -1
        self.table = None

        body = QHBoxLayout(self)
        body.setContentsMargins(0, 0, 0, 0)
        left = QVBoxLayout()
        body.addLayout(left)
        self.list = QListWidget()
        self.list.currentRowChanged.connect(self._on_select)
        left.addWidget(self.list)
        lb = QHBoxLayout()
        b_add = QPushButton("Add head")
        b_add.clicked.connect(self._add_head)
        b_rem = QPushButton("Remove head")
        b_rem.clicked.connect(self._remove_head)
        lb.addWidget(b_add)
        lb.addWidget(b_rem)
        left.addLayout(lb)
        lb2 = QHBoxLayout()
        b_addf = QPushButton("Add flux")
        b_addf.clicked.connect(self._add_flux)
        b_remf = QPushButton("Remove flux")
        b_remf.clicked.connect(self._remove_flux)
        lb2.addWidget(b_addf)
        lb2.addWidget(b_remf)
        left.addLayout(lb2)

        right = QVBoxLayout()
        body.addLayout(right, 1)
        hrow = QHBoxLayout()
        self.head_label = QLabel("Head value:")
        hrow.addWidget(self.head_label)
        self.head_edit = QLineEdit()
        hrow.addWidget(self.head_edit, 1)
        right.addLayout(hrow)
        frow = QHBoxLayout()
        self.flux_label = QLabel("Flux value:")
        frow.addWidget(self.flux_label)
        self.flux_edit = QLineEdit()
        frow.addWidget(self.flux_edit, 1)
        right.addLayout(frow)
        self._holder = QVBoxLayout()
        right.addLayout(self._holder, 1)

        self._refresh()
        self.list.setCurrentRow(0)

    # Index scheme: heads occupy [0, n_heads); fluxes [n_heads, n_heads+n_flux);
    # the exit face is the single trailing row.
    def _is_head(self, idx):
        return 0 <= idx < len(self._heads)

    def _is_flux(self, idx):
        return len(self._heads) <= idx < len(self._heads) + len(self._fluxes)

    def _flux_idx(self, idx):
        return idx - len(self._heads)

    def _is_exit(self, idx):
        return idx == len(self._heads) + len(self._fluxes)

    def _refresh(self):
        self.list.blockSignals(True)
        self.list.clear()
        for i, h in enumerate(self._heads):
            self.list.addItem(f"Head {i + 1}  (h = {h['head']})")
        for i, f in enumerate(self._fluxes):
            self.list.addItem(f"Flux {i + 1}  (q = {f['flux']})")
        self.list.addItem("Exit face")
        self.list.blockSignals(False)

    def _commit(self):
        if self._cur < 0:
            return
        coords = [(r["x"], r["y"]) for r in self.table.result_rows()] if self.table else []
        if self._is_exit(self._cur):
            self._exit = coords
        elif self._is_flux(self._cur):
            f = self._fluxes[self._flux_idx(self._cur)]
            try:
                f["flux"] = float(self.flux_edit.text() or 0)
            except ValueError:
                f["flux"] = 0.0
            f["coords"] = coords
        elif self._is_head(self._cur):
            try:
                self._heads[self._cur]["head"] = float(self.head_edit.text() or 0)
            except ValueError:
                self._heads[self._cur]["head"] = 0.0
            self._heads[self._cur]["coords"] = coords

    def _load(self, idx):
        if self.table is not None:
            self.table.setParent(None)
            self.table = None
        if idx < 0:
            return
        is_head = self._is_head(idx)
        is_flux = self._is_flux(idx)
        self.head_label.setVisible(is_head)
        self.head_edit.setVisible(is_head)
        self.flux_label.setVisible(is_flux)
        self.flux_edit.setVisible(is_flux)
        if is_head:
            self.head_edit.setText(str(self._heads[idx]["head"]))
            rows = [{"x": x, "y": y} for (x, y) in self._heads[idx]["coords"]]
        elif is_flux:
            f = self._fluxes[self._flux_idx(idx)]
            self.flux_edit.setText(str(f["flux"]))
            rows = [{"x": x, "y": y} for (x, y) in f["coords"]]
        else:  # exit face
            rows = [{"x": x, "y": y} for (x, y) in self._exit]
        self.table = _EditableTable(XY_FIELDS, rows, _new_pt)
        self._holder.addWidget(self.table)

    def _on_select(self, idx):
        self._commit()
        self._cur = idx
        self._load(idx)

    def _add_head(self):
        self._commit()
        self._heads.append({"head": 0.0, "coords": []})
        self._cur = -1                     # rows below shifted; avoid a stale commit
        self._refresh()
        self.list.setCurrentRow(len(self._heads) - 1)

    def _remove_head(self):
        idx = self.list.currentRow()
        if not self._is_head(idx):         # only a head row is removable here
            return
        self._heads.pop(idx)
        self._cur = -1
        self._refresh()
        self.list.setCurrentRow(min(idx, self.list.count() - 1))

    def _add_flux(self):
        self._commit()
        self._fluxes.append({"flux": 0.0, "coords": []})
        self._cur = -1
        self._refresh()
        self.list.setCurrentRow(len(self._heads) + len(self._fluxes) - 1)

    def _remove_flux(self):
        idx = self.list.currentRow()
        if not self._is_flux(idx):
            return
        self._fluxes.pop(self._flux_idx(idx))
        self._cur = -1
        self._refresh()
        self.list.setCurrentRow(min(idx, self.list.count() - 1))

    def result(self):
        self._commit()
        return {"specified_heads": [{"head": h["head"], "coords": list(h["coords"])}
                                    for h in self._heads],
                "specified_fluxes": [{"flux": f["flux"], "coords": list(f["coords"])}
                                     for f in self._fluxes],
                "exit_face": list(self._exit)}


class SeepBcEditor(CategoryEditor):
    label = "Seep BC"

    def build(self, slope_data, parent, select=None):
        dlg = QDialog(parent)
        dlg.setWindowTitle("Seep BC")
        dlg.resize(640, 520)
        layout = QVBoxLayout(dlg)
        layout.addWidget(_help_label(
            "Seepage boundary conditions: specified-head boundaries (each a head value + "
            "its points), specified-flux boundaries (each a flux value + its points), and "
            "an exit face (where water leaves the slope). Set 2 is used for rapid-drawdown "
            "(the second seepage solution)."))
        tabs = QTabWidget()
        w1 = _SeepBcSetWidget(slope_data.get("seepage_bc"))
        w2 = _SeepBcSetWidget(slope_data.get("seepage_bc2"))
        tabs.addTab(w1, "Set 1")
        tabs.addTab(w2, "Set 2 (rapid drawdown)")
        layout.addWidget(tabs)
        _ok_cancel(dlg, layout)
        dlg._sets = (w1, w2)
        _apply_set_selection((w1, w2), tabs, select)
        return dlg

    def apply(self, slope_data, dlg):
        w1, w2 = dlg._sets
        slope_data["seepage_bc"] = w1.result()
        slope_data["seepage_bc2"] = w2.result()
        bc2 = slope_data["seepage_bc2"]
        slope_data["has_seepage_bc2"] = bool(bc2.get("specified_heads") or bc2.get("exit_face"))


# --- piles ------------------------------------------------------------------ #
def _new_pile():
    return {"label": "Pile", "x1": 0.0, "y1": 0.0, "x2": 0.0, "y2": 0.0, "H": None,
            "theta_p": 0.0, "D_pile": None, "S": None, "E": None, "I": None,
            "area": None, "V_cap": None, "M_cap": None, "fixity": "free",
            "appl": "active"}


class PilesEditor(CategoryEditor):
    label = "Piles"
    FIELDS = [
        Field("label", "Label", "str"),
        Field("x1", "x1"), Field("y1", "y1"), Field("x2", "x2"), Field("y2", "y2"),
        Field("H", "H", "optfloat"),
        Field("D_pile", "D", "optfloat", usage="lem"), Field("S", "S", "optfloat", usage="lem"),
        Field("E", "E", "optfloat", usage="fem"), Field("I", "I", "optfloat", usage="fem"),
        Field("area", "Area", "optfloat", usage="fem"),
        Field("V_cap", "Vcap", "optfloat", usage="lem"), Field("M_cap", "Mcap", "optfloat", usage="lem"),
        # Force application (v12, LEM only): active = allowable force applied as-is;
        # passive = ultimate capacity divided by FS (loader default 'active').
        Field("appl", "Appl", "choice", choices=["active", "passive"], usage="lem"),
        Field("fixity", "Fixity", "choice", choices=["free", "fixed"], usage="fem"),
    ]

    def build(self, slope_data, parent):
        return TableEditorDialog(
            "Piles", self.FIELDS, slope_data.get("pile_lines", []), _new_pile, parent,
            help_text="Leave H blank for auto Ito & Matsui force. I / Area auto-compute "
                      "from D when blank. θ is auto-derived from the pile axis. Vcap/Mcap "
                      "require S (spacing). Appl: active = allowable force; passive = "
                      "ultimate capacity ÷ FS.",
            usage_toggles=["lem", "fem"])

    def apply(self, slope_data, dlg):
        rows = dlg.result_rows()
        for p in rows:  # θ is auto-derived from the axis (loader does the same)
            dx, dy = p["x2"] - p["x1"], p["y2"] - p["y1"]
            p["theta_p"] = math.degrees(math.atan2(dx, -dy))
        slope_data["pile_lines"] = rows


# --- geometry: profile lines & polygons (master/detail) --------------------- #
def _set_derived_geometry(slope_data, polys):
    """Store material-zone polygons and rebuild the geometry derived from them
    (ground surface, domain polygon, tension-crack line) — the resync the loader
    performs after parsing. ``polys`` is the [{'polygon': Polygon, 'mat_id': int}]
    list both the profile and polygon editors converge on."""
    from shapely.geometry import LineString
    from xslope.fileio import build_ground_surface_from_polygons

    slope_data["polygons"] = polys
    gs, dom = (build_ground_surface_from_polygons(polys) if polys else (None, None))
    slope_data["ground_surface"] = gs
    slope_data["domain_polygon"] = dom
    td = slope_data.get("tcrack_depth", 0)
    if td and td > 0 and gs is not None and not gs.is_empty:
        slope_data["tcrack_surface"] = LineString([(x, y - td) for x, y in gs.coords])
    else:
        slope_data["tcrack_surface"] = None


def _normalize_polygon(p):
    """Coerce a polygon entry to the canonical {'polygon': Polygon, 'mat_id': int}.
    Accepts entries that carry a shapely 'polygon', a 'coords' list, or a bare
    shapely Polygon — so a model that built polygons with only 'coords' still
    resyncs."""
    from shapely.geometry import Polygon
    if isinstance(p, dict):
        poly = p.get("polygon")
        if poly is None and p.get("coords"):
            poly = Polygon(p["coords"])
        return {"polygon": poly, "mat_id": p.get("mat_id")}
    return {"polygon": p, "mat_id": None}


def _resync_geometry(slope_data):
    """Rebuild polygons (and the geometry derived from them: ground surface,
    domain polygon, t-crack line) from the edited source — the resync the loader/
    design driver also do. Two source paths: profile_lines (build zone polygons
    from them), or polygons set directly (e.g. a zoned dam built without
    profile_lines), which are normalized and used as-is."""
    from shapely.geometry import Polygon
    from xslope.mesh import build_polygons

    profile_lines = slope_data.get("profile_lines") or []
    if profile_lines:
        polys = [{"polygon": Polygon(p["coords"]), "mat_id": p["mat_id"]}
                 for p in build_polygons(slope_data={"profile_lines": profile_lines,
                                                      "max_depth": slope_data.get("max_depth")})]
        _set_derived_geometry(slope_data, polys)
        return
    raw = slope_data.get("polygons") or []
    if raw:
        _set_derived_geometry(slope_data, [_normalize_polygon(p) for p in raw])


class MatGeometryDialog(QDialog):
    """Master/detail editor over a list of material-tagged coordinate sequences:
    the items list (left) + the selected item's material and vertex table (right).
    Shared by the profile-line editor (open polylines) and the polygon editor
    (closed rings); the two differ only in labels, help text, and how coords are
    extracted/rebuilt by the owning CategoryEditor."""

    XY = [Field("x", "x"), Field("y", "y")]

    def __init__(self, title, help_text, item_label, items, materials, parent=None,
                 select=None, max_depth=None, preview_draw=None, preview_caption=None,
                 slope_data=None, style=None):
        # items: list of {"mat_id": int|None, "coords": [(x, y), ...]}
        # select: row to pre-highlight (e.g. the double-clicked line); else first.
        # max_depth: when not None, show a "Max depth" field (profile sheet only —
        #   it has no meaning for polygon input); the bottom boundary elevation
        #   used when building zone polygons from the profile lines.
        # preview_draw: draw(ax, lines, selected_index, max_depth) — a per-editor
        #   hook that renders the full section with the selected item highlighted.
        #   When given, a live preview pane is added on the right behind a splitter.
        super().__init__(parent)
        self.setWindowTitle(title)
        self._item_label = item_label
        self._materials = materials
        self._lines = [{"mat_id": it.get("mat_id"),
                        "coords": [tuple(c) for c in it.get("coords", [])]}
                       for it in (items or [])]
        self._cur = -1
        self.table = None
        self._preview_draw = preview_draw
        self._preview = None

        main = QVBoxLayout(self)
        main.addWidget(_help_label(help_text))

        self._max_depth_edit = None
        if max_depth is not None:
            mdrow = QHBoxLayout()
            mdrow.addWidget(QLabel("Max depth (bottom boundary elevation):"))
            self._max_depth_edit = QLineEdit(str(max_depth))
            self._max_depth_edit.setToolTip(
                "Elevation of the model's bottom boundary, used to build the zone "
                "polygons from the profile lines.")
            self._max_depth_edit.setMaximumWidth(90)
            self._max_depth_edit.textChanged.connect(self._schedule_preview)
            mdrow.addWidget(self._max_depth_edit)
            mdrow.addStretch(1)
            main.addLayout(mdrow)

        # Master/detail (list | material+vertex table), plus an optional preview pane
        # on the right behind a splitter — the materials list-view layout pattern.
        left_w = QWidget()
        left = QVBoxLayout(left_w)
        left.setContentsMargins(0, 0, 0, 0)
        self.list = QListWidget()
        self.list.currentRowChanged.connect(self._on_select)
        left.addWidget(self.list)
        lbtns = QHBoxLayout()
        b_add = QPushButton(f"Add {item_label.lower()}")
        b_add.clicked.connect(self._add_line)
        b_rem = QPushButton(f"Remove {item_label.lower()}")
        b_rem.clicked.connect(self._remove_line)
        lbtns.addWidget(b_add)
        lbtns.addWidget(b_rem)
        left.addLayout(lbtns)

        right_w = QWidget()
        right = QVBoxLayout(right_w)
        right.setContentsMargins(0, 0, 0, 0)
        matrow = QHBoxLayout()
        matrow.addWidget(QLabel("Material:"))
        self.mat_combo = QComboBox()
        for i, m in enumerate(materials):
            self.mat_combo.addItem(f"{i + 1}: {m.get('name', '')}")
        self.mat_combo.currentIndexChanged.connect(self._on_mat_changed)
        matrow.addWidget(self.mat_combo, 1)
        right.addLayout(matrow)
        self._holder = QVBoxLayout()
        right.addLayout(self._holder, 1)

        splitter = QSplitter(Qt.Horizontal)
        splitter.addWidget(left_w)
        splitter.addWidget(right_w)
        if preview_draw is not None:
            from .canvas import PreviewPane
            self._preview = PreviewPane(
                lambda ax: self._preview_draw(ax, self._pending_lines(), self._cur,
                                              self.result_max_depth()),
                caption=preview_caption)
            splitter.addWidget(self._preview)
            splitter.setStretchFactor(0, 0)
            splitter.setStretchFactor(1, 0)
            splitter.setStretchFactor(2, 1)
            splitter.setSizes([190, 320, 470])
            self.resize(1060, 560)
        else:
            splitter.setStretchFactor(0, 0)
            splitter.setStretchFactor(1, 1)
            self.resize(680, 540)
        main.addWidget(splitter, 1)

        _ok_cancel(self, main)

        self._refresh_list()
        if self._lines:
            row = select if (select is not None and 0 <= select < len(self._lines)) else 0
            self.list.setCurrentRow(row)
        if self._preview is not None:
            self._preview.refresh_now()

    def _schedule_preview(self, *_):
        if self._preview is not None:
            self._preview.schedule()

    def _pending_lines(self):
        """The current items with the SELECTED one's coords/material taken live from
        the vertex table + material combo (they aren't committed to self._lines until
        a selection change / OK), so the preview reflects the in-progress edit."""
        lines = [{"mat_id": ln["mat_id"], "coords": list(ln["coords"])}
                 for ln in self._lines]
        if 0 <= self._cur < len(lines) and self.table is not None:
            coords = [(r["x"], r["y"]) for r in self.table.result_rows()]
            mid = self.mat_combo.currentIndex()
            lines[self._cur] = {"mat_id": mid if mid >= 0 else lines[self._cur]["mat_id"],
                                "coords": coords}
        return lines

    def _label(self, i):
        mid = self._lines[i]["mat_id"]
        if mid is not None and 0 <= mid < len(self._materials):
            return f"{self._item_label} {i + 1}  (mat {mid + 1}: {self._materials[mid].get('name', '')})"
        return f"{self._item_label} {i + 1}  (mat ?)"

    def _refresh_list(self):
        self.list.blockSignals(True)
        self.list.clear()
        for i in range(len(self._lines)):
            self.list.addItem(self._label(i))
        self.list.blockSignals(False)

    def _commit_current(self):
        if 0 <= self._cur < len(self._lines) and self.table is not None:
            self._lines[self._cur]["coords"] = [(r["x"], r["y"]) for r in self.table.result_rows()]
            if self.mat_combo.currentIndex() >= 0:
                self._lines[self._cur]["mat_id"] = self.mat_combo.currentIndex()

    def _load(self, idx):
        if self.table is not None:
            self.table.setParent(None)
            self.table = None
        if not (0 <= idx < len(self._lines)):
            return
        ln = self._lines[idx]
        rows = [{"x": x, "y": y} for (x, y) in ln["coords"]]
        self.table = _EditableTable(self.XY, rows, _new_pt,
                                    on_change=self._schedule_preview)
        self._holder.addWidget(self.table)
        self.mat_combo.blockSignals(True)
        mid = ln["mat_id"]
        self.mat_combo.setCurrentIndex(mid if (mid is not None and 0 <= mid < self.mat_combo.count()) else 0)
        self.mat_combo.blockSignals(False)

    def _on_select(self, new_idx):
        self._commit_current()
        self._cur = new_idx
        self._load(new_idx)
        self._schedule_preview()

    def _on_mat_changed(self, idx):
        if 0 <= self._cur < len(self._lines) and idx >= 0:
            self._lines[self._cur]["mat_id"] = idx
            item = self.list.item(self._cur)
            if item:
                item.setText(self._label(self._cur))
        self._schedule_preview()

    def _add_line(self):
        self._commit_current()
        self._lines.append({"mat_id": 0, "coords": []})
        self._refresh_list()
        self.list.setCurrentRow(len(self._lines) - 1)
        self._schedule_preview()

    def _remove_line(self):
        idx = self.list.currentRow()
        if idx < 0:
            return
        self._lines.pop(idx)
        self._cur = -1  # invalidate so the next select doesn't write to a stale index
        self._refresh_list()
        if self._lines:
            self.list.setCurrentRow(min(idx, len(self._lines) - 1))
        else:
            self._load(-1)
        self._schedule_preview()

    def result_lines(self):
        self._commit_current()
        return [{"coords": list(ln["coords"]), "mat_id": ln["mat_id"]} for ln in self._lines]

    def result_max_depth(self):
        """The edited max-depth value (float), or None if the field isn't shown."""
        if self._max_depth_edit is None:
            return None
        try:
            return float(self._max_depth_edit.text())
        except (TypeError, ValueError):
            return None


class ProfileEditor(CategoryEditor):
    label = "Profile lines"

    def build(self, slope_data, parent, select=None):
        style = _doc_style(parent)

        def preview(ax, lines, selected, max_depth):
            _draw_profile_preview(ax, lines, selected, max_depth, slope_data, style)

        return MatGeometryDialog(
            "Profile lines",
            "Each profile line is the top of a material layer, drawn left→right and "
            "ordered shallowest first. Select a line to edit its material and vertices.",
            "Line",
            slope_data.get("profile_lines") or [],
            slope_data.get("materials") or [], parent,
            select=select, max_depth=slope_data.get("max_depth") or 0.0,
            preview_draw=preview,
            preview_caption="Preview shows the pending profile lines over the current "
                            "model (selected line highlighted). Material zones aren't "
                            "re-filled until you save.",
            slope_data=slope_data, style=style)

    def apply(self, slope_data, dlg):
        slope_data["profile_lines"] = dlg.result_lines()
        md = dlg.result_max_depth()
        if md is not None:
            slope_data["max_depth"] = md
        _resync_geometry(slope_data)  # rebuild polygons / ground surface / t-crack


class PolygonEditor(CategoryEditor):
    """Geometry editor for polygon-based files (no profile sheet). Each polygon is
    a closed material zone; the loader closes the ring implicitly, so the editor
    shows/stores only the distinct exterior vertices."""

    label = "Polygons"

    def build(self, slope_data, parent, select=None):
        items = []
        for p in (slope_data.get("polygons") or []):
            coords = list(p["polygon"].exterior.coords)
            if len(coords) >= 2 and coords[0] == coords[-1]:
                coords = coords[:-1]                       # drop the closing duplicate
            items.append({"mat_id": p.get("mat_id"), "coords": coords})
        style = _doc_style(parent)

        def preview(ax, polys, selected, _max_depth):
            _draw_polygon_preview(ax, polys, selected, slope_data, style)

        return MatGeometryDialog(
            "Polygons",
            "Each polygon is a closed material zone (the ring is closed automatically, "
            "so list each vertex once). Select a polygon to edit its material and vertices.",
            "Polygon", items, slope_data.get("materials") or [], parent, select=select,
            preview_draw=preview,
            preview_caption="Preview shows the pending zones (selected zone filled and "
                            "outlined, others dimmed).",
            slope_data=slope_data, style=style)

    def apply(self, slope_data, dlg):
        from shapely.geometry import Polygon
        polys = []
        for it in dlg.result_lines():
            coords = it["coords"]
            if len(coords) < 3:
                continue                                   # not a valid ring; skip
            polys.append({"polygon": Polygon(coords), "mat_id": it["mat_id"]})
        _set_derived_geometry(slope_data, polys)


# --- reinforcement ---------------------------------------------------------- #
def _new_reinf():
    # A generic tensile line (blank type -> tangent/active via the loader presets,
    # spacing 1 -> no per-width division); mirrors the pre-v12 default row.
    return {"x1": 0.0, "y1": 0.0, "x2": 0.0, "y2": 0.0, "t_max": 0.0, "t_res": 0.0,
            "lp1": 0.0, "lp2": 0.0, "E": 0.0, "area": 0.0,
            "type": "", "dir": "tangent", "appl": "active",
            "tend1": 0.0, "tend2": 0.0, "spacing": 1.0}


class ReinforcementEditor(CategoryEditor):
    label = "Reinforcement"
    LF = {"lem", "fem"}
    # Support-type columns (v12). `type` choices are the loader's _TYPE_PRESETS keys
    # (fileio.py) plus a blank entry — a blank type is a generic line whose Dir/Appl
    # default via the presets; offering '' as an empty combo entry lets a blank type
    # round-trip unchanged (same treatment as the materials blank option). Dir/Appl
    # mirror the loader's accepted values.
    FIELDS = [
        Field("x1", "x1"), Field("y1", "y1"), Field("x2", "x2"), Field("y2", "y2"),
        Field("type", "Type", "choice",
              choices=["", "geosynthetic", "nail", "tieback", "anchor"], applies=LF),
        Field("dir", "Dir", "choice", choices=["tangent", "axial"], applies=LF),
        Field("appl", "Appl", "choice", choices=["active", "passive"], applies=LF),
        Field("t_max", "Tmax", usage="lem"), Field("t_res", "Tres", usage="fem"),
        Field("tend1", "Tend1", usage="lem"), Field("tend2", "Tend2", usage="lem"),
        Field("lp1", "Lp1", usage="lem"), Field("lp2", "Lp2", usage="lem"),
        Field("spacing", "Spacing", applies=LF),
        Field("E", "E", usage="fem"), Field("area", "Area", usage="fem"),
    ]

    def build(self, slope_data, parent):
        return TableEditorDialog(
            "Reinforcement", self.FIELDS, slope_data.get("reinforcement_lines", []),
            _new_reinf, parent,
            help_text="Lp1/Lp2 are the pullout lengths at each end (0 = fully anchored). "
                      "The LEM tension distribution shown on the plot is derived from these. "
                      "Type defaults Dir/Appl when blank; leave Type blank for a generic "
                      "tensile line. Tend1/Tend2 are the end-anchorage capacities; capacities "
                      "and E/Area are per-unit-width (Spacing divides discrete supports).",
            usage_toggles=["lem", "fem"])

    def apply(self, slope_data, dlg):
        from xslope.fileio import build_reinforce_lines
        rows = dlg.result_rows()
        slope_data["reinforcement_lines"] = rows
        # Rebuild the LEM display/analysis format so the canvas reflects the edit.
        slope_data["reinforce_lines"] = build_reinforce_lines(rows)


# --- line loads ------------------------------------------------------------- #
def _new_lload():
    # Same key set load_slope_data produces (fileio.py:1325-1330). Angle defaults
    # to -90° (straight down), mirroring the loader's blank-angle fallback.
    return {"x": 0.0, "y": 0.0, "P": 0.0, "angle": -90.0, "label": "Load"}


class LineLoadsEditor(CategoryEditor):
    label = "Line loads"
    # Loader keys (fileio.py:1325-1330): x, y, P (magnitude), angle, label. A line
    # load is a concentrated force per unit width acting at a point on the ground
    # surface; used by both LEM (slice.py) and FEM (fem.py), so no usage tags.
    FIELDS = [
        Field("label", "Label", "str"),
        Field("x", "x"), Field("y", "y"),
        Field("P", "P"),
        Field("angle", "Angle", default=-90.0),
    ]

    def build(self, slope_data, parent):
        return TableEditorDialog(
            "Line loads", self.FIELDS, slope_data.get("line_loads", []),
            _new_lload, parent,
            help_text="A concentrated line load (force per unit width) acting at a "
                      "point on the ground surface — e.g. the weight of a facing "
                      "plate. P is the magnitude (positive); Angle is the direction "
                      "in degrees (−90 = straight down). Each load is snapped onto "
                      "the ground surface on save, since the loader requires line "
                      "loads to act on the surface.")

    def apply(self, slope_data, dlg):
        rows = dlg.result_rows()
        # Mirror the loader (fileio.py:1315-1324), which requires a line load to act
        # ON the ground surface: snap each point onto the nearest surface location so
        # Studio can never author a load the loader would refuse. (The loader also
        # rejects points beyond a tolerance; snapping is the forgiving inverse — the
        # saved point always lands exactly on the surface.)
        gs = slope_data.get("ground_surface")
        if gs is not None and not gs.is_empty:
            from shapely.geometry import Point
            for ll in rows:
                try:
                    snapped = gs.interpolate(gs.project(Point(ll["x"], ll["y"])))
                    ll["x"], ll["y"] = float(snapped.x), float(snapped.y)
                except Exception:
                    pass
        slope_data["line_loads"] = rows


CATEGORY_EDITORS = {
    "global": GlobalEditor(),
    "materials": MaterialsEditor(),
    "circles": CirclesEditor(),
    "non_circ": NonCircEditor(),
    "piezo": PiezoEditor(),
    "dloads": DloadsEditor(),
    "seep_bc": SeepBcEditor(),
    "piles": PilesEditor(),
    "reinforce": ReinforcementEditor(),
    "line_loads": LineLoadsEditor(),
    "profile": ProfileEditor(),
    "polygons": PolygonEditor(),
}
