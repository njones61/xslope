"""MplCanvas — renders the engine's Matplotlib figure and shows it in a
QGraphicsView so the WHOLE figure (axes, title, labels, margins) zooms and pans
uniformly, like an image/PDF viewer — not the axes-only data zoom of the
Matplotlib navigation toolbar.

The figure is still produced by the engine's plotting (``plot_inputs(fig=…)``);
we render it to a pixmap and let the view handle zoom (wheel / buttons, anchored
at the cursor), pan (drag), and fit-to-window.
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

from xslope.plot import plot_inputs

ZOOM_STEP = 1.25


class MplCanvas(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.figure = Figure(figsize=(12, 6), dpi=150)
        self._agg = FigureCanvasAgg(self.figure)
        self._pixitem = None
        self._fitted = False

        self.scene = QGraphicsScene(self)
        self.view = QGraphicsView(self.scene)
        self.view.setDragMode(QGraphicsView.ScrollHandDrag)
        self.view.setRenderHints(QPainter.Antialiasing | QPainter.SmoothPixmapTransform)
        self.view.setTransformationAnchor(QGraphicsView.AnchorUnderMouse)
        self.view.setBackgroundBrush(Qt.lightGray)
        self.view.viewport().installEventFilter(self)

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
        plot_inputs(slope_data, fig=self.figure, mode=mode, mat_table=mat_table)
        self._update_pixmap()

    def clear(self):
        self.figure.clear()
        self._update_pixmap()

    def _update_pixmap(self):
        self._agg.draw()
        w, h = self._agg.get_width_height()
        img = QImage(bytes(self._agg.buffer_rgba()), w, h, QImage.Format_RGBA8888)
        pix = QPixmap.fromImage(img)  # fromImage copies, so the buffer can be freed
        if self._pixitem is None:
            self._pixitem = self.scene.addPixmap(pix)
        else:
            self._pixitem.setPixmap(pix)
        self.scene.setSceneRect(QRectF(pix.rect()))
        if not self._fitted:
            # Fit once the view has a real size (after layout settles).
            QTimer.singleShot(0, self.fit)
            self._fitted = True

    # --- zoom / pan ------------------------------------------------------
    def fit(self):
        if self._pixitem is not None:
            self.view.fitInView(self._pixitem, Qt.KeepAspectRatio)

    def reset_100(self):
        self.view.resetTransform()

    def zoom_in(self):
        self.view.scale(ZOOM_STEP, ZOOM_STEP)

    def zoom_out(self):
        self.view.scale(1 / ZOOM_STEP, 1 / ZOOM_STEP)

    def eventFilter(self, obj, event):
        if obj is self.view.viewport() and event.type() == QEvent.Wheel:
            factor = ZOOM_STEP if event.angleDelta().y() > 0 else 1 / ZOOM_STEP
            self.view.scale(factor, factor)
            return True
        return super().eventFilter(obj, event)
