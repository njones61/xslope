"""MplCanvas — an embedded Matplotlib canvas with the standard nav toolbar.

Phase 1 renders the Inputs view by handing the engine's ``plot_inputs`` our own
Figure (the ``fig=`` parameter added to plot.py), so the GUI reuses the existing
plotting code unchanged. The Matplotlib NavigationToolbar provides zoom / pan /
home for free; scroll-wheel zoom is layered on top.
"""

from __future__ import annotations

import os

# Ensure Matplotlib's Qt backend binds to PySide6 (must precede backend import).
os.environ.setdefault("QT_API", "pyside6")

from PySide6.QtWidgets import QVBoxLayout, QWidget
from matplotlib.backends.backend_qtagg import (
    FigureCanvasQTAgg as FigureCanvas,
    NavigationToolbar2QT as NavigationToolbar,
)
from matplotlib.figure import Figure

from xslope.plot import plot_inputs


class MplCanvas(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.figure = Figure(figsize=(12, 6))
        self.canvas = FigureCanvas(self.figure)
        self.toolbar = NavigationToolbar(self.canvas, self)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(self.toolbar)
        layout.addWidget(self.canvas)

        self.canvas.mpl_connect("scroll_event", self._on_scroll)

    def render_inputs(self, slope_data, mode="lem"):
        """Draw the Inputs view for slope_data into the embedded figure."""
        plot_inputs(slope_data, fig=self.figure, mode=mode, mat_table=True)
        self.canvas.draw_idle()

    def clear(self):
        self.figure.clear()
        self.canvas.draw_idle()

    # --- scroll-wheel zoom about the cursor ------------------------------
    def _on_scroll(self, event):
        ax = event.inaxes
        if ax is None or event.xdata is None or event.ydata is None:
            return
        scale = 0.9 if event.button == "up" else 1.1
        x0, x1 = ax.get_xlim()
        y0, y1 = ax.get_ylim()
        ax.set_xlim(event.xdata - (event.xdata - x0) * scale,
                    event.xdata + (x1 - event.xdata) * scale)
        ax.set_ylim(event.ydata - (event.ydata - y0) * scale,
                    event.ydata + (y1 - event.ydata) * scale)
        self.canvas.draw_idle()
