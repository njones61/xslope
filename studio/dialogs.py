"""Run-options dialogs (Phase 3).

``RunLemDialog`` collects the options for an LEM solve. The thin first slice
covers a **single circular surface**; the analysis-type / surface-type selectors
are present but fixed to that case for now, and widen as the later increments
(auto-search, reliability, non-circular) land.
"""

from __future__ import annotations

from PySide6.QtWidgets import (
    QCheckBox, QComboBox, QDialog, QDialogButtonBox, QFormLayout, QLabel,
    QSpinBox, QVBoxLayout,
)

LEM_METHODS = [
    ("oms", "Ordinary Method of Slices (OMS)"),
    ("bishop", "Bishop's Simplified"),
    ("janbu", "Janbu (Corrected)"),
    ("corps", "Corps of Engineers"),
    ("lowe", "Lowe & Karafiath"),
    ("spencer", "Spencer"),
    ("mprice", "Morgenstern-Price"),
]


class RunLemDialog(QDialog):
    """Options for a single-surface circular LEM solve."""

    def __init__(self, parent=None, defaults=None):
        super().__init__(parent)
        self.setWindowTitle("Run LEM")
        defaults = defaults or {}

        layout = QVBoxLayout(self)
        form = QFormLayout()

        self.method = QComboBox()
        for key, label in LEM_METHODS:
            self.method.addItem(label, key)
        self._select_method(defaults.get("method", "bishop"))
        form.addRow("Method", self.method)

        self.num_slices = QSpinBox()
        self.num_slices.setRange(5, 500)
        self.num_slices.setValue(int(defaults.get("num_slices", 40)))
        form.addRow("Number of slices", self.num_slices)

        self.rapid = QCheckBox("Rapid drawdown")
        self.rapid.setChecked(bool(defaults.get("rapid", False)))
        form.addRow("", self.rapid)

        layout.addLayout(form)
        note = QLabel("Analyzes the first circular surface (single-surface). "
                      "Auto-search and reliability arrive in the next increment.")
        note.setWordWrap(True)
        layout.addWidget(note)

        bb = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        bb.button(QDialogButtonBox.Ok).setText("Run")
        bb.accepted.connect(self.accept)
        bb.rejected.connect(self.reject)
        layout.addWidget(bb)

    def _select_method(self, key):
        idx = self.method.findData(key)
        if idx >= 0:
            self.method.setCurrentIndex(idx)

    def options(self):
        return {
            "method": self.method.currentData(),
            "num_slices": self.num_slices.value(),
            "rapid": self.rapid.isChecked(),
        }
