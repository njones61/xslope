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
Markdown embeds directly. It lives apart from that script only because it loads
two solved FEM models first, which the Studio capture script never does; the
capture helpers here are deliberately identical to it so the two can be merged
into one script once that step is worth carrying there.

Each sample is read back from the companions committed beside it
(``tools/make_fem_docs_sidecars.py`` writes them from the page's own tag), which
is the run the sample page reports and the mechanism it reports it at. It used to
solve each model here at a fixed strength reduction instead, which reached the
right factor and captured no mechanism at all — so the panel was grabbed on a
converged field, without the shear band's crossing or the stretches of capacity
the at-failure state is where the members hold.

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

# The member each sample's screenshot opens on, and the size the panel is grabbed
# at. The run itself is whatever ships beside the model, so the numbers in the
# screenshot are the ones the rest of the sample page reports.
CASES = [
    ("xslope_reinforce_fem", "Line 4", "reinforce_fem_details.png", (1140, 660)),
    ("xslope_piles_fem", None, "piles_fem_details.png", (1140, 660)),
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


def _shipped(stem):
    """The solved run committed beside one sample, as Studio's own bundle."""
    from xslope.fileio import load_slope_data
    from xslope.report import solutions_from_sidecars
    path = os.path.join(FILES_DIR, f"{stem}.xlsx")
    slope_data = load_slope_data(path)
    bundle = solutions_from_sidecars(path, slope_data, None).get("fem")
    if not bundle:
        raise RuntimeError(f"{stem}: no solved run beside the model — run "
                           f"tools/make_fem_docs_sidecars.py first")
    return path, slope_data, bundle


def capture_details(stem, select_label, name, size):
    from PySide6.QtCore import Qt
    from studio.fem_details_dialog import FemDetailsDialog

    path, slope_data, bundle = _shipped(stem)
    dlg = FemDetailsDialog(bundle["fem_data"], bundle["solution"], slope_data,
                           model_path=path,
                           failure_solution=bundle.get("failure_solution"))
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
    for stem, label, name, size in CASES:
        out = capture_details(stem, label, name, size)
        print(f"wrote {os.path.relpath(out, REPO_ROOT)}")


if __name__ == "__main__":
    main()
