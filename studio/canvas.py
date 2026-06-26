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
    QGraphicsScene, QGraphicsView, QHBoxLayout, QToolButton, QVBoxLayout, QWidget,
)
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure

from xslope.plot import (
    plot_circular_search_results, plot_inputs, plot_noncircular_search_results,
    plot_solution,
)

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

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addLayout(bar)
        layout.addWidget(self.view)

    # --- rendering -------------------------------------------------------
    def render_inputs(self, slope_data, mode="lem", mat_table=False):
        self._draw(lambda fig: plot_inputs(slope_data, fig=fig, mode=mode,
                                           mat_table=mat_table))

    def render_solution(self, slope_data, slice_df, failure_surface, results):
        self._draw(lambda fig: plot_solution(slope_data, slice_df, failure_surface,
                                             results, fig=fig))

    def render_search(self, slope_data, search):
        """Render auto-search results (all trial surfaces + critical + search path)."""
        if search["kind"] == "circular":
            self._draw(lambda fig: plot_circular_search_results(
                slope_data, search["fs_cache"], search["search_path"],
                circle_cache=search["circle_cache"], fig=fig))
        else:
            self._draw(lambda fig: plot_noncircular_search_results(
                slope_data, search["fs_cache"], search["search_path"], fig=fig))

    def _draw(self, draw_fn):
        """Populate the embedded figure via ``draw_fn(fig)`` and rasterize it.
        Fitting is deferred to ``ensure_fitted`` so a canvas drawn while its tab is
        hidden (no viewport size yet) still fits when the tab is first shown."""
        draw_fn(self.figure)
        self._rasterize(self._target_dpi())
        QTimer.singleShot(0, self.ensure_fitted)
        self._schedule_refine()

    def ensure_fitted(self):
        """Fit the figure to the window the first time the view has a real size."""
        if not self._fitted and self._pixitem is not None \
                and self.view.viewport().width() > 1:
            self.fit()
            self._fitted = True

    def reset_fit(self):
        """Re-arm the one-shot fit so the next render fits to the window (e.g. when
        a different file is loaded), rather than keeping the previous view."""
        self._fitted = False

    def clear(self):
        self.figure.clear()
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
        if self._pixitem is not None:
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
