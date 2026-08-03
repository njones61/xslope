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

"""Generate the reinforcement / pile details dialog screenshots used in the docs.

Offscreen, like ``tools/capture_studio_screenshots.py``: the same
``QT_QPA_PLATFORM=offscreen`` + show/settle/grab sequence, writing PNGs the
Markdown embeds directly. It lives apart from that script only because it solves
two FEM models first, which the Studio capture script never does; the capture
helpers here are deliberately identical to it so the two can be merged into one
script once the solve step is worth carrying there.

    python tools/capture_1d_details.py
"""

import os
import sys
import time

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import matplotlib                                                   # noqa: E402
matplotlib.use("Agg")

try:
    from PySide6.QtWidgets import QApplication
except Exception:                      # engine-only install — no studio layer
    print("capture_1d_details: PySide6 not installed — skipped.")
    sys.exit(0)

OUT_DIR = os.path.join(REPO_ROOT, "docs", "fem", "images")
FILES_DIR = os.path.join(REPO_ROOT, "docs", "fem", "files")

_app = QApplication.instance() or QApplication([])

# Each sample is solved at its own SSRM factor of safety — the state its page
# already documents and plots — so the numbers in the screenshot match the ones
# the rest of the sample page reports.
CASES = [
    ("xslope_reinforce_fem", 1.49, "Line 4", "reinforce_fem_details.png", (1140, 660)),
    ("xslope_piles_fem", 1.36, None, "piles_fem_details.png", (1140, 660)),
]


def _settle(cycles=14):
    """Pump the event loop so deferred layout + the canvas's debounced render fire."""
    for _ in range(cycles):
        _app.processEvents()
        time.sleep(0.02)


def _grab(dlg, name):
    """Show, settle and grab a dialog to ``docs/fem/images/<name>``."""
    dlg.show()
    _settle()
    out = os.path.join(OUT_DIR, name)
    dlg.grab().save(out)
    dlg.close()
    return out


def _solve(stem, F):
    from xslope.fem import build_fem_data, solve_fem
    from xslope.fileio import load_slope_data
    path = os.path.join(FILES_DIR, f"{stem}.xlsx")
    slope_data = load_slope_data(path)
    fem_data = build_fem_data(slope_data, slope_data.get("mesh"))
    solution = solve_fem(fem_data, F=F, debug_level=0)
    return path, slope_data, fem_data, solution


def capture_details(stem, F, select_label, name, size):
    from PySide6.QtCore import Qt
    from studio.fem_details_dialog import FemDetailsDialog

    path, slope_data, fem_data, solution = _solve(stem, F)
    dlg = FemDetailsDialog(fem_data, solution, slope_data, model_path=path)
    dlg.resize(*size)
    if select_label:
        for row in range(dlg.list.count()):
            entry = dlg.list.item(row).data(Qt.UserRole)
            if entry and entry["label"] == select_label:
                dlg.list.setCurrentRow(row)
                break
    return _grab(dlg, name)


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    for stem, F, label, name, size in CASES:
        out = capture_details(stem, F, label, name, size)
        print(f"wrote {os.path.relpath(out, REPO_ROOT)}")


if __name__ == "__main__":
    main()
