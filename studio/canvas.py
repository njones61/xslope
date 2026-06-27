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

from PySide6.QtCore import QEvent, QRectF, Qt, QTimer
from PySide6.QtGui import QImage, QPainter, QPixmap
from PySide6.QtWidgets import (
    QFileDialog, QGraphicsScene, QGraphicsView, QHBoxLayout, QInputDialog,
    QMessageBox, QToolButton, QVBoxLayout, QWidget,
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
SUPERSAMPLE = 2.0     # render this many device px per on-screen px (crisp text)
MIN_DPI = 100
MAX_DPI = 600         # caps pixmap memory (12×6in @600dpi RGBA ≈ 100 MB)
REFINE_MS = 90        # debounce before re-rasterizing after a zoom change


class MplCanvas(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.figure = Figure(figsize=(12, 6), dpi=BASE_DPI)
        self._agg = FigureCanvasAgg(self.figure)
        self._pixitem = None
        self._fitted = False
        self._render_dpi = 0          # DPI the current pixmap was rasterized at
        self._content_rect = None     # tight bbox of inked content, in scene coords
        self._dxf_supported = False   # current view exports to DXF (set per render)

        self.scene = QGraphicsScene(self)
        self.view = QGraphicsView(self.scene)
        self.view.setDragMode(QGraphicsView.ScrollHandDrag)
        self.view.setRenderHints(QPainter.Antialiasing | QPainter.SmoothPixmapTransform)
        self.view.setTransformationAnchor(QGraphicsView.AnchorUnderMouse)
        self.view.setBackgroundBrush(Qt.white)
        self.view.viewport().installEventFilter(self)

        # Debounce re-rasterization so a flurry of wheel ticks redraws once.
        self._refine_timer = QTimer(self)
        self._refine_timer.setSingleShot(True)
        self._refine_timer.setInterval(REFINE_MS)
        self._refine_timer.timeout.connect(self._refine_resolution)

        bar = QHBoxLayout()
        for text, tip, slot in [
            ("Fit", "Fit figure to window", self.fit),
            ("100%", "Actual size", self.reset_100),
            ("+", "Zoom in", self.zoom_in),
            ("−", "Zoom out", self.zoom_out),
        ]:
            btn = QToolButton()
            btn.setText(text)
            btn.setToolTip(tip)
            btn.clicked.connect(slot)
            bar.addWidget(btn)
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
    def render_inputs(self, slope_data, mode="lem", opts=None):
        opts = opts or {}
        self._draw(lambda fig: plot_inputs(
            slope_data, fig=fig, mode=mode,
            mat_table=opts.get("mat_table", False),
            tab_loc=opts.get("tab_loc", "top"),
            legend_ncol=opts.get("legend_ncol", "auto")))

    def render_solution(self, slope_data, slice_df, failure_surface, results, opts=None):
        opts = opts or {}
        self._draw(lambda fig: plot_solution(
            slope_data, slice_df, failure_surface, results,
            slice_numbers=opts.get("slice_numbers", False),
            seep_contours=opts.get("seep_contours", True),
            legend_ncol=opts.get("legend_ncol", "auto"), fig=fig))

    def render_search(self, slope_data, search, opts=None):
        """Render auto-search results (all trial surfaces + critical + search path)."""
        opts = opts or {}
        highlight = opts.get("highlight_fs", True)
        ncol = opts.get("legend_ncol", "auto")
        if search["kind"] == "circular":
            self._draw(lambda fig: plot_circular_search_results(
                slope_data, search["fs_cache"], search["search_path"],
                circle_cache=search["circle_cache"], highlight_fs=highlight,
                legend_ncol=ncol, fig=fig))
        else:
            self._draw(lambda fig: plot_noncircular_search_results(
                slope_data, search["fs_cache"], search["search_path"],
                highlight_fs=highlight, legend_ncol=ncol, fig=fig))

    def render_reliability(self, slope_data, reliability_data, opts=None):
        opts = opts or {}
        self._draw(lambda fig: plot_reliability_results(
            slope_data, reliability_data,
            legend_ncol=opts.get("legend_ncol", "auto"), fig=fig))

    def render_mesh(self, mesh, materials=None, opts=None):
        opts = opts or {}
        self._draw(lambda fig: plot_mesh(
            mesh, materials=materials, show_nodes=opts.get("show_nodes", True),
            label_elements=opts.get("label_elements", False),
            label_nodes=opts.get("label_nodes", False),
            pad_frac=opts.get("pad_frac", 0.05),
            legend_ncol=opts.get("legend_ncol", "auto"), fig=fig))

    def render_seep_data(self, seep_data, opts=None):
        opts = opts or {}
        self._draw(lambda fig: plot_seep_data(
            seep_data, show_nodes=opts.get("show_nodes", False),
            show_bc=opts.get("show_bc", True),
            label_elements=opts.get("label_elements", False),
            label_nodes=opts.get("label_nodes", False),
            alpha=opts.get("alpha", 0.4),
            legend_ncol=opts.get("legend_ncol", "auto"), fig=fig))

    def render_seep_solution(self, seep_data, solution, opts):
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
            phreatic=opts.get("phreatic", True), mesh=False, fig=fig))

    def render_fem_data(self, fem_data, opts=None):
        opts = opts or {}
        self._draw(lambda fig: plot_fem_data(
            fem_data, show_nodes=opts.get("show_nodes", False),
            show_bc=opts.get("show_bc", True),
            label_elements=opts.get("label_elements", False),
            label_nodes=opts.get("label_nodes", False),
            alpha=opts.get("alpha", 0.4),
            bc_symbol_size=opts.get("bc_symbol_size", 0.03),
            legend_ncol=opts.get("legend_ncol", "auto"), fig=fig))

    def render_fem_results(self, fem_data, solution, opts):
        opts = opts or {}
        self._draw(lambda fig: plot_fem_results(
            fem_data, solution,
            plot_type=[opts.get("plot_type", "shear_strain")],
            deform_percent=opts.get("deform_percent", 15),
            show_mesh=opts.get("show_mesh", True),
            show_reinforcement=opts.get("show_reinforcement", True),
            label_elements=opts.get("label_elements", False),
            plot_boundary=opts.get("plot_boundary", True),
            plot_nodes=opts.get("plot_nodes", False),
            plot_elements=opts.get("plot_elements", False),
            scale_vectors=opts.get("scale_vectors", True),
            displacement_tolerance=opts.get("displacement_tolerance", 0.5),
            legend_ncol=opts.get("legend_ncol", "auto"), fig=fig))

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
        """Populate the embedded figure via ``draw_fn(fig)`` and rasterize it.

        ``dxf`` records whether this view's geometry can be exported to DXF (every
        engine plot here can except an empty canvas), gating the DXF option in the
        Save dialog.

        Every (re)render re-arms the fit (``_fitted = False``) so the new content
        is fitted to the window — a fresh solve/result always autofits, not just
        the first plot ever drawn. The fit itself is deferred to ``ensure_fitted``
        (and to ``showEvent`` for a canvas drawn while its tab is hidden), so it
        runs once the viewport has a real size rather than at zero size."""
        self._dxf_supported = dxf
        draw_fn(self.figure)
        self._rasterize(self._target_dpi())
        self._fitted = False
        QTimer.singleShot(0, self.ensure_fitted)
        self._schedule_refine()

    def ensure_fitted(self):
        """Fit the figure to the window once the view has a real size.

        A result canvas is often rendered while its tab is hidden, and when the
        tab is first shown the viewport may not be laid out yet (width 0). So:
        only act while visible, and if not laid out yet, retry shortly — that
        covers the "click the LEM Solution tab" case where the synchronous fit on
        tab-change runs before the page is sized. The retry stops once fitted, and
        never runs for a hidden tab (it waits for showEvent)."""
        if self._fitted or self._pixitem is None or not self.isVisible():
            return
        if self.view.viewport().width() > 1:
            self.fit()
            self._fitted = True
        else:
            QTimer.singleShot(40, self.ensure_fitted)

    def reset_fit(self):
        """Re-arm the one-shot fit so the next render fits to the window (e.g. when
        a different file is loaded), rather than keeping the previous view."""
        self._fitted = False

    def clear(self):
        self.figure.clear()
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

    def _target_dpi(self):
        """DPI needed so the pixmap has SUPERSAMPLE device px per on-screen px at
        the current zoom level."""
        scale = self.view.transform().m11() or 1.0
        return BASE_DPI * scale * SUPERSAMPLE

    def _schedule_refine(self):
        self._refine_timer.start()

    def _refine_resolution(self):
        """After a zoom change, re-rasterize if the pixmap is now too coarse for
        the on-screen size (pixelated) or wastefully fine (free memory)."""
        if self._pixitem is None:
            return
        need = max(MIN_DPI, min(MAX_DPI, self._target_dpi()))
        if need > self._render_dpi * 1.05 or need < self._render_dpi * 0.6:
            self._rasterize(need)

    # --- zoom / pan ------------------------------------------------------
    def fit(self):
        if self._content_rect is not None and not self._content_rect.isEmpty():
            self.view.fitInView(self._content_rect, Qt.KeepAspectRatio)
            self._schedule_refine()
        elif self._pixitem is not None:
            self.view.fitInView(self._pixitem, Qt.KeepAspectRatio)
            self._schedule_refine()

    def reset_100(self):
        self.view.resetTransform()
        self._schedule_refine()

    def zoom_in(self):
        self.view.scale(ZOOM_STEP, ZOOM_STEP)
        self._schedule_refine()

    def zoom_out(self):
        self.view.scale(1 / ZOOM_STEP, 1 / ZOOM_STEP)
        self._schedule_refine()

    def showEvent(self, event):
        super().showEvent(event)
        # A result canvas is usually drawn while its tab is hidden (no viewport
        # size yet), so the initial fit is skipped. Re-attempt it the moment the
        # view becomes visible — deferred one cycle so the page has been laid out.
        # Guarded by _fitted, so revisiting a tab won't clobber the user's zoom.
        if not self._fitted:
            QTimer.singleShot(0, self.ensure_fitted)

    def resizeEvent(self, event):
        super().resizeEvent(event)
        # Re-fit the whole figure to the new window size (re-rasterize debounced).
        self.fit()

    def eventFilter(self, obj, event):
        if obj is self.view.viewport() and event.type() == QEvent.Wheel:
            factor = ZOOM_STEP if event.angleDelta().y() > 0 else 1 / ZOOM_STEP
            self.view.scale(factor, factor)
            self._schedule_refine()
            return True
        return super().eventFilter(obj, event)
