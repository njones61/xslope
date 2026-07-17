"""MplCanvas — renders the engine's Matplotlib figure and shows it in a
QGraphicsView so the WHOLE figure (axes, title, labels, margins) zooms and pans
uniformly, like an image/PDF viewer — not the axes-only data zoom of the
Matplotlib navigation toolbar.

The figure is still produced by the engine's plotting (``plot_inputs(fig=…)``).
To keep text crisp at any zoom (a plain pixmap pixelates as you zoom in), the
figure is re-rasterized **on demand**: the scene works in fixed *logical* units
(inches × ``BASE_DPI``), but the backing pixmap is redrawn at a DPI matched to the
current zoom (× ``SUPERSAMPLE``), so the bitmap always has enough pixels for the
on-screen size instead of being upscaled. Re-draws are debounced and the DPI is
capped to bound memory.
"""

from __future__ import annotations

import os

os.environ.setdefault("QT_API", "pyside6")

import math

from PySide6.QtCore import QEvent, QPoint, QPointF, QRectF, QSize, Qt, QTimer, Signal
from PySide6.QtGui import (
    QColor, QGuiApplication, QIcon, QImage, QPainter, QPen, QPixmap,
)
from PySide6.QtWidgets import (
    QFileDialog, QGraphicsScene, QGraphicsView, QHBoxLayout, QInputDialog,
    QLabel, QMessageBox, QToolButton, QVBoxLayout, QWidget,
)
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure

from xslope.plot import (
    plot_circular_search_results, plot_inputs, plot_mesh,
    plot_noncircular_search_results, plot_reliability_results, plot_solution,
)
from xslope.plot_seep import plot_seep_data, plot_seep_solution
from xslope.plot_fem import plot_fem_data, plot_fem_results

ZOOM_STEP = 1.25
BASE_DPI = 100        # logical scene units per inch (1 unit ≈ 1 screen px at 100%)
SUPERSAMPLE = 1.0     # bitmap px per *device* px. 1.0 renders the figure at exactly
                      # the device resolution and shows it 1:1, so Matplotlib's own
                      # text antialiasing stays crisp with no bilinear resample
                      # (>1 forces a downscale that grains thin text). The tight
                      # refine band re-rasterizes on zoom so lines don't alias.
MIN_DPI = 100
MAX_DPI = 900         # caps pixmap memory (raised for Retina, where the effective
                      # DPI is doubled by devicePixelRatio)
REFINE_MS = 90        # debounce before re-rasterizing after a zoom change
RENDER_DEBOUNCE_MS = 50   # debounce before rendering after draw/show/resize, so the
                          # file-open / drag-resize layout burst settles first (one
                          # render at the final size, no wrong-size-then-refit flash)
FIT_PAD = 0.04        # cushion around the content when fitting (fraction per side)


def _magnifier_icon(px=20, color="#2b2b2b"):
    """Draw a magnifying-glass icon (lens circle + handle) at ``px`` logical size.

    A drawn icon (rather than a font glyph like "⌕", which renders tiny inside its
    em-box) fills the button predictably and stays crisp on Retina via the 2×
    backing pixmap."""
    dpr = 2
    pm = QPixmap(px * dpr, px * dpr)
    pm.setDevicePixelRatio(dpr)
    pm.fill(Qt.transparent)
    p = QPainter(pm)
    p.setRenderHint(QPainter.Antialiasing, True)
    pen = QPen(QColor(color))
    pen.setWidthF(px * 0.12)
    pen.setCapStyle(Qt.RoundCap)
    p.setPen(pen)
    r = px * 0.30                      # lens radius
    cx = cy = px * 0.38                # lens centre
    p.drawEllipse(QPointF(cx, cy), r, r)
    a = math.radians(45)               # handle runs out from the lens at 45°
    p.drawLine(QPointF(cx + r * math.cos(a), cy + r * math.sin(a)),
               QPointF(px * 0.86, px * 0.86))
    p.end()
    return QIcon(pm)


class MplCanvas(QWidget):
    # Emitted on a double-click that lands inside the main axes: (data_x, data_y,
    # tol) where tol is ~a dozen screen px expressed in data units. The Inputs
    # view connects this to open the editor for the feature under the cursor.
    picked = Signal(float, float, float)
    # Emitted on a SINGLE click (not part of a pan/zoom-box drag) that lands inside
    # the main axes, when click-picking is enabled (the editor preview panes). Same
    # (data_x, data_y, tol) payload, with tol from ~8 screen px — so zooming in
    # shrinks the pick tolerance in data units and fine picks become precise.
    clicked = Signal(float, float, float)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.figure = Figure(figsize=(12, 7), dpi=BASE_DPI)
        self._agg = FigureCanvasAgg(self.figure)
        self._pixitem = None
        self._rendered_vp = None      # viewport (w, h) the current draw_fn was last
                                      # rasterized at; None = not yet rendered (or a
                                      # forced re-fit). _render_current renders only
                                      # when this differs from the live viewport, so a
                                      # settle/resize re-renders but a user zoom (same
                                      # size) is preserved, and a hidden tab (no size)
                                      # stays pending until shown.
        self._render_dpi = 0          # DPI the current pixmap was rasterized at
        self._content_rect = None     # tight bbox of inked content, in scene coords
        self._dxf_supported = False   # current view exports to DXF (set per render)
        self._zoom_box_mode = False   # drag-a-rectangle zoom (vs. drag-to-pan)
        self._band_scene = None       # scene rect of the in-progress rubber band
        self._pick_enabled = False    # double-click selects a feature (Inputs view)
        self._click_pick_enabled = False  # single-click emits `clicked` (preview panes)
        self._press_vp = None         # viewport point of the last left-press, so a
                                      # click can be told apart from a pan/box drag

        self.scene = QGraphicsScene(self)
        self.view = QGraphicsView(self.scene)
        self.view.setDragMode(QGraphicsView.ScrollHandDrag)
        self.view.setRenderHints(QPainter.Antialiasing | QPainter.SmoothPixmapTransform)
        self.view.setTransformationAnchor(QGraphicsView.AnchorUnderMouse)
        self.view.setBackgroundBrush(Qt.white)
        # Never show scrollbars. The figure is sized to the viewport and fit 1:1, so
        # a scrollbar only ever appears transiently (rounding/settle) — and when it
        # does it steals ~16px of viewport, shrinking it, which forces a re-raster at
        # the smaller size, which drops the scrollbar, which grows it back… a
        # feedback loop seen as a flash on load. Panning is by drag (ScrollHandDrag),
        # so hidden scrollbars lose nothing.
        self.view.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        self.view.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        self.view.viewport().installEventFilter(self)
        # In zoom-box mode the view uses RubberBandDrag; this fires as the band is
        # dragged (and once with a null rect on release) so we can fit to it.
        self.view.rubberBandChanged.connect(self._on_rubber_band)

        # Debounce re-rasterization so a flurry of wheel ticks redraws once.
        self._refine_timer = QTimer(self)
        self._refine_timer.setSingleShot(True)
        self._refine_timer.setInterval(REFINE_MS)
        self._refine_timer.timeout.connect(self._refine_resolution)

        # The figure is sized to the viewport, so any layout change (new draw, tab
        # shown, window/dock resize) needs a re-raster at the new size. Route them ALL
        # through one debounced timer: a burst of layout changes — a drag-resize, or
        # the flurry during file-open (result tabs added, Display dock switching
        # panels) — then coalesces into a SINGLE render at the settled size, instead
        # of rendering at a transient size and visibly re-rastering when it settles
        # (the load flash). The debounce restarts on every trigger, so it fires only
        # once the layout has stopped changing.
        self._draw_fn = None
        self._render_timer = QTimer(self)
        self._render_timer.setSingleShot(True)
        self._render_timer.setInterval(RENDER_DEBOUNCE_MS)
        self._render_timer.timeout.connect(self._render_current)

        bar = QHBoxLayout()
        for text, tip, slot in [
            # The figure is sized to the viewport, so "Fit" shows the whole plot at
            # an exact 1:1 (crisp); a separate "100%" would be identical, so it's gone.
            ("Fit", "Fit plot to window", self.fit),
            ("+", "Zoom in", self.zoom_in),
            ("−", "Zoom out", self.zoom_out),
        ]:
            btn = QToolButton()
            btn.setText(text)
            btn.setToolTip(tip)
            btn.clicked.connect(slot)
            bar.addWidget(btn)
        # Zoom-to-box: a checkable mode that swaps drag-to-pan for drag-a-rectangle.
        self._zoom_box_btn = QToolButton()
        self._zoom_box_btn.setIcon(_magnifier_icon())
        self._zoom_box_btn.setToolTip("Zoom to box — drag a rectangle to zoom into it")
        self._zoom_box_btn.setCheckable(True)
        self._zoom_box_btn.toggled.connect(self._toggle_zoom_box)
        # Match the text buttons' height; an icon button is otherwise taller. Size
        # the glyph to fill that height (minus a small margin) so it stays prominent.
        _bh = btn.sizeHint().height()
        self._zoom_box_btn.setFixedHeight(_bh)
        self._zoom_box_btn.setIconSize(QSize(_bh - 4, _bh - 4))
        bar.addWidget(self._zoom_box_btn)
        # Pick hint, just right of the zoom tools — shown only on the Inputs view
        # (pick enabled) and hidden while the zoom-box tool is active.
        bar.addSpacing(10)
        self._hint_label = QLabel("(double-click on a feature to edit)")
        self._hint_label.setStyleSheet("color: gray;")
        self._hint_label.setVisible(False)
        bar.addWidget(self._hint_label)
        bar.addStretch(1)
        # Export the current figure to an image file (per-view, so it sits next to
        # the image it saves and works for views without a Display panel).
        save_btn = QToolButton()
        save_btn.setText("Save…")
        save_btn.setToolTip("Export this view to a PNG / PDF / SVG / DXF file")
        save_btn.clicked.connect(self.save_image)
        bar.addWidget(save_btn)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addLayout(bar)
        layout.addWidget(self.view)

    # --- rendering -------------------------------------------------------
    def render_inputs(self, slope_data, mode="lem", opts=None, style=None):
        opts = opts or {}
        self._draw(lambda fig: plot_inputs(
            slope_data, fig=fig, mode=mode,
            mat_table=opts.get("mat_table", False),
            tab_loc=opts.get("tab_loc", "top"),
            label_coordinates=opts.get("label_coordinates", False),
            coord_label_size=opts.get("coord_label_size", 7),
            coord_arrows=opts.get("coord_arrows", False),
            legend_ncol=opts.get("legend_ncol", "auto"), legend_frame=opts.get("legend_frame", False), show_title=opts.get("show_title", True), show_legend=opts.get("show_legend", True), style=style))

    def render_solution(self, slope_data, slice_df, failure_surface, results, opts=None, style=None):
        opts = opts or {}
        self._draw(lambda fig: plot_solution(
            slope_data, slice_df, failure_surface, results,
            slice_numbers=opts.get("slice_numbers", False),
            seep_contours=opts.get("seep_contours", True),
            legend_ncol=opts.get("legend_ncol", "auto"), legend_frame=opts.get("legend_frame", False), show_title=opts.get("show_title", True), show_legend=opts.get("show_legend", True), fig=fig, style=style))

    def render_search(self, slope_data, search, opts=None, style=None):
        """Render auto-search results (all trial surfaces + critical + search path)."""
        opts = opts or {}
        highlight = opts.get("highlight_fs", True)
        ncol = opts.get("legend_ncol", "auto")
        if search["kind"] == "circular":
            self._draw(lambda fig: plot_circular_search_results(
                slope_data, search["fs_cache"], search["search_path"],
                circle_cache=search["circle_cache"], highlight_fs=highlight,
                legend_ncol=ncol, legend_frame=opts.get("legend_frame", False), show_title=opts.get("show_title", True), show_legend=opts.get("show_legend", True), fig=fig, style=style))
        else:
            self._draw(lambda fig: plot_noncircular_search_results(
                slope_data, search["fs_cache"], search["search_path"],
                highlight_fs=highlight, legend_ncol=ncol, legend_frame=opts.get("legend_frame", False), show_title=opts.get("show_title", True), show_legend=opts.get("show_legend", True), fig=fig, style=style))

    def render_reliability(self, slope_data, reliability_data, opts=None, style=None):
        opts = opts or {}
        self._draw(lambda fig: plot_reliability_results(
            slope_data, reliability_data,
            legend_ncol=opts.get("legend_ncol", "auto"), legend_frame=opts.get("legend_frame", False), show_title=opts.get("show_title", True), show_legend=opts.get("show_legend", True), fig=fig, style=style))

    def render_mesh(self, mesh, materials=None, opts=None, style=None):
        opts = opts or {}
        self._draw(lambda fig: plot_mesh(
            mesh, materials=materials, show_nodes=opts.get("show_nodes", True),
            label_elements=opts.get("label_elements", False),
            label_nodes=opts.get("label_nodes", False),
            pad_frac=opts.get("pad_frac", 0.05),
            legend_ncol=opts.get("legend_ncol", "auto"), legend_frame=opts.get("legend_frame", False), show_title=opts.get("show_title", True), show_legend=opts.get("show_legend", True), fig=fig, style=style))

    def render_seep_data(self, seep_data, opts=None, style=None):
        opts = opts or {}
        self._draw(lambda fig: plot_seep_data(
            seep_data, show_nodes=opts.get("show_nodes", False),
            show_bc=opts.get("show_bc", True),
            label_elements=opts.get("label_elements", False),
            label_nodes=opts.get("label_nodes", False),
            alpha=opts.get("alpha", 0.6),
            legend_ncol=opts.get("legend_ncol", "auto"), legend_frame=opts.get("legend_frame", False), show_title=opts.get("show_title", True), show_legend=opts.get("show_legend", True), fig=fig, style=style))

    def render_seep_solution(self, seep_data, solution, opts, style=None):
        opts = opts or {}
        self._draw(lambda fig: plot_seep_solution(
            seep_data, solution, variable=opts.get("variable", "head"),
            levels=opts.get("levels", 20), base_mat=opts.get("base_mat", 1),
            alpha=opts.get("alpha", 0.4),
            vector_scale=opts.get("vector_scale", 0.05),
            pad_frac=opts.get("pad_frac", 0.05),
            flowlines=opts.get("flowlines", True),
            vectors=opts.get("vectors", False),
            fill_contours=opts.get("fill_contours", False),
            phreatic=opts.get("phreatic", True),
            cmap=opts.get("cmap", "Spectral_r"),
            cbar_shrink=opts.get("cbar_shrink", 0.8),
            legend_ncol=opts.get("legend_ncol", "auto"),
            legend_frame=opts.get("legend_frame", False), show_title=opts.get("show_title", True), show_legend=opts.get("show_legend", True), mesh=False, fig=fig, style=style))

    def render_fem_data(self, fem_data, opts=None, style=None):
        opts = opts or {}
        self._draw(lambda fig: plot_fem_data(
            fem_data, show_nodes=opts.get("show_nodes", False),
            show_bc=opts.get("show_bc", True),
            label_elements=opts.get("label_elements", False),
            label_nodes=opts.get("label_nodes", False),
            alpha=opts.get("alpha", 0.6),
            bc_symbol_size=opts.get("bc_symbol_size", 0.03),
            legend_ncol=opts.get("legend_ncol", "auto"), legend_frame=opts.get("legend_frame", False), show_title=opts.get("show_title", True), show_legend=opts.get("show_legend", True), fig=fig, style=style))

    def render_fem_results(self, fem_data, solution, opts):
        opts = opts or {}
        self._draw(lambda fig: plot_fem_results(
            fem_data, solution,
            plot_type=[opts.get("plot_type", "shear_strain")],
            cmap=opts.get("cmap"),
            cbar_shrink=opts.get("cbar_shrink"),
            deform_percent=opts.get("deform_percent", 15),
            show_mesh=opts.get("show_mesh", True),
            show_reinforcement=opts.get("show_reinforcement", True),
            label_elements=opts.get("label_elements", False),
            plot_boundary=opts.get("plot_boundary", True),
            plot_nodes=opts.get("plot_nodes", False),
            plot_elements=opts.get("plot_elements", False),
            scale_vectors=opts.get("scale_vectors", True),
            displacement_tolerance=opts.get("displacement_tolerance", 0.5),
            legend_ncol=opts.get("legend_ncol", "auto"), legend_frame=opts.get("legend_frame", False), show_title=opts.get("show_title", True), show_legend=opts.get("show_legend", True), fig=fig))

    def render_axes(self, plot_fn, dxf=False):
        """Render a single-axes figure by calling ``plot_fn(ax)`` — a lightweight
        entry point for property plots (e.g. the Materials list view's strength /
        kr panels) that want one Axes rather than a full engine figure. Rides the
        same deferred, viewport-sized, 1:1 render path as the engine plots, so the
        panel stays crisp and needs no hand-tuned sizing. ``dxf`` is False by
        default (these plots aren't geometry, so the Save dialog omits DXF)."""
        def _draw(fig):
            ax = fig.add_subplot(111)
            plot_fn(ax)
            try:
                fig.tight_layout()
            except Exception:
                pass
        self._draw(_draw, dxf=dxf)

    # --- export ----------------------------------------------------------
    def save_image(self, _checked=False, suggested_name=""):
        """Export the current figure to a file. PNG prompts for a DPI; PDF/SVG/DXF
        are vector and need none. Raster/vector image formats save at the figure's
        true size in inches, so output resolution is independent of the on-screen
        zoom. DXF re-emits the rendered geometry layer-by-layer via the engine's
        ``axes_to_dxf`` (the same path as the plot functions' ``save_dxf`` kwarg),
        offered only for views whose engine plot supports it (§ not Reliability)."""
        filters = ["PNG image (*.png)", "PDF document (*.pdf)", "SVG image (*.svg)"]
        if self._dxf_supported:
            filters.append("DXF drawing (*.dxf)")
        path, sel = QFileDialog.getSaveFileName(
            self, "Save view", suggested_name, ";;".join(filters))
        if not path:
            return
        # The on-screen render is debounced; if it hasn't fired yet, populate the
        # figure now so we export the current content, not stale/empty axes.
        if self._render_timer.isActive() or (
                self._draw_fn is not None and self._rendered_vp is None):
            self._render_timer.stop()
            self._render_current()
        ext = os.path.splitext(path)[1].lower()
        if not ext:                       # user typed no extension — infer from filter
            ext = {"PDF document (*.pdf)": ".pdf",
                   "SVG image (*.svg)": ".svg",
                   "DXF drawing (*.dxf)": ".dxf"}.get(sel, ".png")
            path += ext
        try:
            if ext == ".dxf":
                ax = self._main_axes()
                if ax is None:
                    raise RuntimeError("no geometry axes to export")
                from xslope.cad import axes_to_dxf
                axes_to_dxf(ax, path)
            else:
                dpi = 300
                if ext not in (".pdf", ".svg"):   # raster formats need a resolution
                    dpi, ok = QInputDialog.getInt(self, "Image resolution",
                                                  "DPI:", 300, 50, 1200, 50)
                    if not ok:
                        return
                # Restore the on-screen render DPI afterwards so the live pixmap is
                # unaffected by the export DPI.
                self.figure.savefig(path, dpi=dpi, bbox_inches="tight")
                self.figure.set_dpi(self._render_dpi or BASE_DPI)
        except Exception:
            import traceback
            traceback.print_exc()
            QMessageBox.warning(self, "Save view failed",
                                "Could not write the file — see the Log pane.")
            return
        # Visible confirmation via the window's status bar if available.
        win = self.window()
        if hasattr(win, "statusBar"):
            win.statusBar().showMessage(f"Saved {os.path.basename(path)}")

    def _main_axes(self):
        """The largest axes on the figure — the main geometry plot, which is what
        the engine's ``axes_to_dxf`` expects (legends/colorbars/tables are smaller
        sibling axes). Returns None for an empty figure."""
        axs = self.figure.get_axes()
        if not axs:
            return None
        def area(a):
            bb = a.get_position()
            return bb.width * bb.height
        return max(axs, key=area)

    def _draw(self, draw_fn, dxf=True):
        """Store ``draw_fn(fig)`` as this canvas's content and schedule a render.

        ``dxf`` records whether this view's geometry can be exported to DXF (every
        engine plot here can except an empty canvas), gating the DXF option in the
        Save dialog.

        The raster itself is done by ``_render_current`` (debounced), at the canvas's
        real, settled size. New content resets ``_rendered_vp`` so it fits afresh; if
        the canvas is hidden or still laying out when the timer fires, the render
        stays pending until it's shown/settled — so the figure is never rasterized at
        the wrong size and then visibly re-fitted (the "loads wrong, flashes, refits"
        flash)."""
        self._dxf_supported = dxf
        self._draw_fn = draw_fn          # remembered so a resize can re-layout it
        self._rendered_vp = None         # new content: render/fit afresh
        self._schedule_render()

    def _schedule_render(self):
        """(Re)start the debounce so a burst of layout changes coalesces into one
        render at the settled size. Every trigger (draw, show, resize, tab-change)
        routes through here."""
        self._render_timer.start()

    def _render_current(self):
        """Rasterize the stored draw_fn at the current viewport size, 1:1.

        Renders only when visible and laid out (real size) AND the size differs from
        the last raster — so a settle/resize re-renders, a hidden/zero-size tab stays
        pending (``_rendered_vp`` left None; a later show/resize reschedules), and a
        re-trigger at an unchanged size is a no-op that preserves the user's zoom."""
        if self._draw_fn is None:
            return
        vp = self.view.viewport()
        w, h = vp.width(), vp.height()
        if not self.isVisible() or w <= 1 or h <= 1:
            return                       # can't render yet; show/resize will retry
        if self._rendered_vp == (w, h):
            return                       # unchanged size: keep the pixmap and any zoom
        # Reset the view to 1:1 BEFORE rasterizing so (a) the pixmap lands at the
        # fitted transform from its first frame — not briefly at the previous view's
        # zoom/pan — and (b) _target_dpi() reads scale 1.0, so the first raster is
        # already at the final DPI and the refine pass doesn't swap it. Both of those
        # single-frame swaps would read as a flash.
        self.view.resetTransform()
        self.figure.set_size_inches(w / BASE_DPI, h / BASE_DPI, forward=False)
        self.figure.clear()
        try:
            self._draw_fn(self.figure)
        except Exception:
            return
        self._rasterize(self._target_dpi())
        self._restore_pan_cursor()
        self._rendered_vp = (w, h)
        self._schedule_refine()

    def ensure_fitted(self):
        """Re-render at the current size if it changed since the last raster — e.g.
        this tab was last drawn at a different viewport size, or was drawn while
        hidden and is now shown. A no-op when the size is unchanged, so a user's zoom
        is preserved. Debounced. Called when a tab becomes the current one."""
        self._schedule_render()

    def reset_fit(self):
        """Force the next render to re-fit to the window (e.g. when a different file
        is loaded), discarding any zoom from the previous view."""
        self._rendered_vp = None

    def clear(self):
        self.figure.clear()
        self._draw_fn = None
        self._rendered_vp = None
        self._dxf_supported = False
        self._content_rect = None
        self._rasterize(BASE_DPI)

    def _rasterize(self, dpi):
        """(Re)draw the figure to a pixmap at ``dpi`` and place it at its fixed
        logical footprint (inches × BASE_DPI), so a higher DPI yields more pixels
        for the same on-screen size rather than changing the layout."""
        dpi = max(MIN_DPI, min(MAX_DPI, int(dpi)))
        self.figure.set_dpi(dpi)
        self._agg.draw()
        w, h = self._agg.get_width_height()
        img = QImage(bytes(self._agg.buffer_rgba()), w, h, QImage.Format_RGBA8888)
        pix = QPixmap.fromImage(img)  # fromImage copies, so the buffer can be freed
        if self._pixitem is None:
            self._pixitem = self.scene.addPixmap(pix)
            self._pixitem.setTransformationMode(Qt.SmoothTransformation)
        else:
            self._pixitem.setPixmap(pix)
        # Shrink the dpi-resolution pixmap to the logical footprint so scene
        # coordinates stay independent of the render DPI.
        self._pixitem.setScale(BASE_DPI / dpi)
        w_in, h_in = self.figure.get_size_inches()
        self.scene.setSceneRect(QRectF(0, 0, w_in * BASE_DPI, h_in * BASE_DPI))
        self._render_dpi = dpi
        self._content_rect = self._tight_content_rect(w_in, h_in)

    def _tight_content_rect(self, w_in, h_in):
        """The inked content's bounding box in scene coords (inches × BASE_DPI),
        so Fit frames the actual plot rather than the whole figure — which, with
        equal-aspect axes, leaves wide margins on one side. Mirrors what Save's
        bbox_inches='tight' crops to. Returns None if it can't be computed (Fit
        then falls back to the full figure)."""
        try:
            tb = self.figure.get_tightbbox(self._agg.get_renderer())  # inches, btm-left
        except Exception:
            return None
        # Matplotlib bbox origin is bottom-left; the scene's is top-left, so flip y.
        rect = QRectF(tb.x0 * BASE_DPI, (h_in - tb.y1) * BASE_DPI,
                      tb.width * BASE_DPI, tb.height * BASE_DPI)
        return rect.intersected(QRectF(0, 0, w_in * BASE_DPI, h_in * BASE_DPI))

    def _device_pixel_ratio(self):
        """The display's device-pixel ratio (2.0 on Retina). Read from the screen
        rather than the widget's devicePixelRatioF(), which returns 1.0 when the
        widget has no backing store yet — e.g. a result tab rendered while hidden,
        which is exactly when result canvases are first drawn."""
        scr = self.screen() or QGuiApplication.primaryScreen()
        return scr.devicePixelRatio() if scr is not None else 1.0

    def _target_dpi(self):
        """DPI needed so the pixmap has SUPERSAMPLE bitmap px per *device* px at
        the current zoom. Includes the display's devicePixelRatio so text stays
        crisp on Retina/HiDPI screens (where one logical px is 2+ device px)."""
        scale = self.view.transform().m11() or 1.0
        return BASE_DPI * scale * self._device_pixel_ratio() * SUPERSAMPLE

    def _schedule_refine(self):
        self._refine_timer.start()

    def _refine_resolution(self):
        """After a zoom/fit change, re-rasterize if the pixmap no longer matches
        the on-screen size: too coarse (pixelated/upscaled) OR too fine (a large
        downscale, which bilinear sampling renders grainy — and wastes memory).
        The initial render happens before Fit, so the post-fit scale almost always
        differs; keep the band tight so the bitmap tracks it."""
        if self._pixitem is None:
            return
        need = max(MIN_DPI, min(MAX_DPI, self._target_dpi()))
        ratio = need / self._render_dpi
        # Tight band so the bitmap tracks the on-screen size closely — any residual
        # fractional up/down-scale by the view transform softens text. Re-rasterize
        # (debounced) on a small zoom-in especially, where upscaling blurs most.
        if ratio > 1.02 or ratio < 0.88:
            self._rasterize(need)

    # --- zoom / pan ------------------------------------------------------
    def _can_pan(self):
        """True when the scene is larger than the viewport in either axis, i.e.
        there is somewhere to scroll to. After a Fit the content fills the window
        with no slack, so the hand cursor would be misleading — there's nothing to
        pan. The scrollbar ranges are the authoritative signal for ScrollHandDrag."""
        h = self.view.horizontalScrollBar()
        v = self.view.verticalScrollBar()
        return h.minimum() != h.maximum() or v.minimum() != v.maximum()

    def _restore_pan_cursor(self):
        """Set the viewport cursor for the active mode. fitInView / resetTransform
        / scale otherwise leave it as a plain arrow until the next mouse press.
        Box-zoom → cross; a pick-enabled view (Inputs) → a standard arrow (select)
        rather than the pan grab; otherwise pan → open hand when there's room to
        pan (a fully-fitted view shows the normal arrow instead)."""
        if self._zoom_box_mode:
            cursor = Qt.CrossCursor
        elif self._pick_enabled:
            cursor = Qt.ArrowCursor
        elif self._can_pan():
            cursor = Qt.OpenHandCursor
        else:
            cursor = Qt.ArrowCursor
        self.view.viewport().setCursor(cursor)

    def set_pick_enabled(self, on):
        """Enable double-click-to-select on this canvas (the Inputs view). Shows a
        standard arrow (select) cursor instead of the pan/grab hand, plus a hint."""
        self._pick_enabled = bool(on)
        self._restore_pan_cursor()
        self._update_hint()

    def set_click_pick_enabled(self, on):
        """Enable single-click-to-select on this canvas (the editor preview panes):
        a click that isn't part of a pan/zoom-box drag emits `clicked`. Panning
        (drag-to-scroll) and box-zoom still work — only a stationary click picks."""
        self._click_pick_enabled = bool(on)

    def _update_hint(self):
        """The double-click hint is shown only when picking is active — i.e. on
        the Inputs view and not while the zoom-box tool has taken over the drag."""
        self._hint_label.setVisible(self._pick_enabled and not self._zoom_box_mode)

    def _toggle_zoom_box(self, on):
        """Switch between drag-to-pan (ScrollHandDrag) and drag-a-rectangle zoom
        (RubberBandDrag). Stays active until toggled off, so several box-zooms can
        be done in a row; the user toggles it off to pan again."""
        self._zoom_box_mode = on
        self.view.setDragMode(QGraphicsView.RubberBandDrag if on
                              else QGraphicsView.ScrollHandDrag)
        self._restore_pan_cursor()
        self._update_hint()

    def _on_rubber_band(self, rect, from_pt, to_pt):
        """Capture the rubber-band extent while dragging and, on release, zoom to
        it. On release the rect is null and the end-points come back null too, so
        the scene rect is remembered from the last non-null update. Tiny bands
        (an accidental click) are ignored."""
        if not rect.isNull():
            self._band_scene = QRectF(from_pt, to_pt).normalized()
            return
        band, self._band_scene = self._band_scene, None
        if self._zoom_box_mode and band is not None \
                and band.width() > 2 and band.height() > 2:
            self.view.fitInView(band, Qt.KeepAspectRatio)
            self._schedule_refine()

    def fit(self):
        # The figure is sized to the viewport, so the whole plot already spans the
        # canvas at 100%. Force a fresh 1:1 render at the current size (resetTransform
        # happens inside _render_current) — not fitInView, whose fractional scale
        # (~0.99) smooth-scales the bitmap and softens it. Immediate (not debounced):
        # the Fit button should respond at once.
        if self._pixitem is None:
            return
        self._render_timer.stop()
        self._rendered_vp = None        # force a re-render even at the current size
        self._render_current()
        self._restore_pan_cursor()

    def reset_100(self):
        self.view.resetTransform()
        self._schedule_refine()
        self._restore_pan_cursor()

    def zoom_in(self):
        self.view.scale(ZOOM_STEP, ZOOM_STEP)
        self._schedule_refine()
        self._restore_pan_cursor()

    def zoom_out(self):
        self.view.scale(1 / ZOOM_STEP, 1 / ZOOM_STEP)
        self._schedule_refine()
        self._restore_pan_cursor()

    def showEvent(self, event):
        super().showEvent(event)
        # First time shown (hidden result tab) or shown at a new size: schedule a
        # (debounced) render. _render_current no-ops if the size is unchanged, so
        # switching to an already-rendered tab keeps its pixmap and any zoom.
        self._schedule_render()

    def resizeEvent(self, event):
        super().resizeEvent(event)
        # The figure is sized to the viewport, so a resize re-lays it out at the new
        # size (debounced, coalescing a drag) — see _schedule_render / _render_current.
        self._schedule_render()

    def eventFilter(self, obj, event):
        if obj is self.view.viewport() and event.type() == QEvent.NativeGesture:
            # Trackpad pinch-to-zoom (macOS): value() is the incremental
            # magnification — positive when spreading fingers (zoom in), negative
            # when pinching (zoom out), anchored under the cursor.
            if event.gestureType() == Qt.ZoomNativeGesture:
                factor = 1.0 + event.value()
                if factor > 0:
                    self.view.scale(factor, factor)
                    self._schedule_refine()
                    self._restore_pan_cursor()
                return True
            return False
        if obj is self.view.viewport() and event.type() == QEvent.Wheel:
            # A real mouse wheel (angle-only, no pixelDelta) zooms; a trackpad
            # two-finger scroll (pixelDelta set) falls through to the view's
            # normal scrolling so it pans instead of zooming. Pinch arrives as a
            # NativeGesture, handled above.
            if not event.pixelDelta().isNull():
                return super().eventFilter(obj, event)
            factor = ZOOM_STEP if event.angleDelta().y() > 0 else 1 / ZOOM_STEP
            self.view.scale(factor, factor)
            self._schedule_refine()
            self._restore_pan_cursor()
            return True
        if (obj is self.view.viewport()
                and event.type() == QEvent.MouseButtonDblClick
                and event.button() == Qt.LeftButton and not self._zoom_box_mode):
            self._emit_pick(event.position().toPoint())
            return True
        if (obj is self.view.viewport() and self._click_pick_enabled
                and event.type() == QEvent.MouseButtonPress
                and event.button() == Qt.LeftButton):
            # Remember where the press landed so the matching release can tell a
            # stationary click (a pick) from a pan/box drag. Don't consume — the
            # view still needs the press to start a pan (ScrollHandDrag).
            self._press_vp = event.position().toPoint()
        if (obj is self.view.viewport() and self._click_pick_enabled
                and event.type() == QEvent.MouseButtonRelease
                and event.button() == Qt.LeftButton and not self._zoom_box_mode):
            press, self._press_vp = self._press_vp, None
            rel = event.position().toPoint()
            # A click that barely moved is a pick; a drag (pan) is not. Threshold in
            # viewport px so it's zoom-independent.
            if press is not None and (rel - press).manhattanLength() <= 4:
                self._emit_click_pick(rel)
            # Don't consume — let the view finish its own release handling.
        if (obj is self.view.viewport() and self._pick_enabled
                and event.type() == QEvent.MouseButtonRelease):
            # ScrollHandDrag resets the cursor to the open hand on release; re-assert
            # the select cursor after Qt's own handler runs. Don't consume the event.
            QTimer.singleShot(0, self._restore_pan_cursor)
        return super().eventFilter(obj, event)

    # --- pick (double-click to edit) -------------------------------------
    def _emit_pick(self, vp_point):
        """Map a viewport double-click to axes data coords and emit `picked`,
        with a tolerance set from a ~12px offset so picking stays pixel-accurate
        at any zoom (more zoomed in → tighter tolerance)."""
        hit = self._viewport_to_data(vp_point)
        if hit is None:
            return
        x, y = hit
        off = self._viewport_to_data(vp_point + QPoint(12, 0))
        tol = math.hypot(x - off[0], y - off[1]) if off else 0.0
        self.picked.emit(x, y, tol)

    def _emit_click_pick(self, vp_point, tol_px=8):
        """Map a single-click viewport point to axes data coords and emit `clicked`,
        with a tolerance from an `tol_px`-pixel offset — so the ~8px pick radius is
        expressed in data units at the CURRENT zoom (zoomed in → tighter tolerance →
        precise fine picks). Silent when the click falls outside the main axes."""
        hit = self._viewport_to_data(vp_point)
        if hit is None:
            return
        x, y = hit
        off = self._viewport_to_data(vp_point + QPoint(tol_px, 0))
        tol = math.hypot(x - off[0], y - off[1]) if off else 0.0
        self.clicked.emit(x, y, tol)

    def _viewport_to_data(self, vp_point):
        """Map a point in viewport pixels to (x, y) in the main axes' data coords,
        or None if there's no figure or the click falls outside the axes. Goes
        viewport → scene → figure fraction → display px → data; the transFigure /
        transData pair is DPI-independent so it works at the current render DPI."""
        ax = self._main_axes()
        if ax is None or self._pixitem is None:
            return None
        scene_pt = self.view.mapToScene(vp_point)
        w_in, h_in = self.figure.get_size_inches()
        W, H = w_in * BASE_DPI, h_in * BASE_DPI
        if W <= 0 or H <= 0:
            return None
        # scene origin is top-left, y down; Matplotlib figure fraction is y up.
        disp = self.figure.transFigure.transform(
            (scene_pt.x() / W, 1.0 - scene_pt.y() / H))
        try:
            bb = ax.get_window_extent()
            pad = 6.0
            if not (bb.x0 - pad <= disp[0] <= bb.x1 + pad
                    and bb.y0 - pad <= disp[1] <= bb.y1 + pad):
                return None
        except Exception:
            pass
        x, y = ax.transData.inverted().transform(disp)
        return float(x), float(y)


class PreviewPane(QWidget):
    """A debounced live preview for the input editors: an ``MplCanvas`` plus an
    optional caption, driven by a per-editor ``draw(ax)`` hook that renders the full
    cross-section with the object being edited highlighted (the shared mechanism
    behind the profile-line / polygon / starting-circle previews — the Materials
    list-view plots generalized).

    Editors call :meth:`schedule` on any field / table / selection change; renders
    are debounced (so the preview never lags typing) and ride ``render_axes``'
    deferred, viewport-sized path. The ``draw`` hook reads the editor's LIVE pending
    state, so the picture always reflects the latest edit. If the hook raises — a
    half-typed, unparseable row mid-edit — ``MplCanvas`` keeps the last good pixmap,
    so a blank cell never blanks or errors the preview."""

    # Forwarded from the inner canvas: a single click (not a pan/box drag) inside the
    # axes, as (data_x, data_y, tol). Editors connect this to their pick resolver so a
    # click in the preview selects the object/row under it.
    clicked = Signal(float, float, float)

    def __init__(self, draw_fn, caption=None, debounce_ms=140, dxf=True, parent=None):
        super().__init__(parent)
        self._draw_fn = draw_fn
        self._dxf = dxf
        self.canvas = MplCanvas()
        self.canvas.set_click_pick_enabled(True)   # click-to-select in the preview
        self.canvas.clicked.connect(self.clicked)   # re-expose at the pane level
        v = QVBoxLayout(self)
        v.setContentsMargins(0, 0, 0, 0)
        v.addWidget(self.canvas, 1)
        if caption:
            cap = QLabel(caption)
            cap.setWordWrap(True)
            cap.setStyleSheet("color: gray; font-size: 11px;")
            v.addWidget(cap)
        self._timer = QTimer(self)
        self._timer.setSingleShot(True)
        self._timer.setInterval(debounce_ms)
        self._timer.timeout.connect(self._render)

    def schedule(self):
        """(Re)start the debounce so a burst of edits coalesces into one render."""
        self._timer.start()

    def refresh_now(self):
        """Render immediately (skip the debounce) — used for the first paint."""
        self._timer.stop()
        self._render()

    def _render(self):
        self.canvas.render_axes(self._draw_fn, dxf=self._dxf)
