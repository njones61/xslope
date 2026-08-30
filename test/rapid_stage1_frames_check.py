"""The Parametric dialog's FS-vs-time frame list: with Rapid drawdown ticked, the
frames at or before stage 1 are unticked and dimmed (they are the state the
others fall from, not drawdowns), All/None leave them alone, and unticking Rapid
drawdown gives them back ticked."""
import os
import sys

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from PySide6.QtCore import Qt                                 # noqa: E402
from PySide6.QtWidgets import QApplication                    # noqa: E402

FILE = os.path.join(ROOT, "docs/tutorials/files/xslope_johnson_fs_time.xlsx")
TIMES = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 80, 100, 130, 170,
         220, 300, 400, 500]


def main():
    app = QApplication.instance() or QApplication(sys.argv)   # noqa: F841
    from xslope.fileio import load_slope_data
    from studio.dialogs import SensitivityDialog
    sd = load_slope_data(FILE)
    dlg = SensitivityDialog(defaults={"method": "spencer", "mode": "fs_vs_time"},
                            slope_data=sd, app_mode="lem",
                            transient={"times": TIMES})
    first = dlg.times.item(0)

    def enabled():
        return bool(first.flags() & Qt.ItemIsEnabled)

    assert enabled() and first.checkState() == Qt.Checked
    assert len(dlg.selected_times()) == len(TIMES)
    dlg.rapid.setChecked(True)
    assert not enabled(), "t = 0 stays enabled with Rapid drawdown on"
    assert first.checkState() == Qt.Unchecked
    assert dlg.selected_times() == [float(t) for t in TIMES if t > 0]
    assert "Stage 1" in first.toolTip()
    dlg._set_all_times(True)
    assert first.checkState() == Qt.Unchecked, "All ticked the stage-1 frame"
    dlg.rapid.setChecked(False)
    assert enabled() and first.checkState() == Qt.Checked
    assert len(dlg.selected_times()) == len(TIMES)
    print("rapid_stage1_frames_check: OK")


if __name__ == "__main__":
    main()
