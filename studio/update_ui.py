"""Help → Check for Updates: the Qt half of the updater.

:mod:`studio.updater` decides everything (where the manifest is, whether a
version is newer, which artifact belongs to this machine, what command an
install runs); this module is only the wiring — a controller the main window
owns, two background runners from :mod:`studio.runners`, and the dialogs.

Two entry points, and the difference between them is the whole design:

``check_manual()``
    The user clicked the menu item, so there is ALWAYS exactly one dialog at the
    end — up to date, a new version, or "could not check". A wait cursor and a
    status-bar line stand in for progress; the check itself is a few kilobytes.

``check_at_startup()``
    Nobody asked. It runs at most once a day (``updates/last_check`` in the app's
    QSettings), only when the preference under Help is on, and it is SILENT in
    every outcome except one: a genuinely newer version, which raises a
    non-modal notification the user can ignore. An offline laptop, a proxy, a
    404 — none of them ever produce a dialog on this path.

The update action itself branches on ``sys.frozen``. A pip install is shown the
``pip install -U`` line with a copy button and nothing is executed. A frozen app
downloads its platform artifact with the same status-bar progress bar and Cancel
button a solve uses, verifies the sha256 against ``latest.json``, and only then
runs the platform's install command.
"""

from __future__ import annotations

import os
import tempfile

from PySide6.QtCore import QObject, Qt, QTimer, QUrl
from PySide6.QtGui import QDesktopServices
from PySide6.QtWidgets import (
    QApplication, QDialog, QDialogButtonBox, QHBoxLayout, QLabel, QLineEdit,
    QMessageBox, QPushButton, QVBoxLayout,
)

from . import updater
from .runners import UpdateCheckRunner, UpdateDownloadRunner


class UpdateAvailableDialog(QDialog):
    """"A newer version exists" — the one dialog both entry points share.

    Non-modal when it comes from the startup check (it appeared unbidden, so it
    must not block the window), modal when the user asked. The action row is
    whatever this install can actually do: **Download & Install** for a frozen
    app, the literal upgrade command with a **Copy** button for a pip install,
    and **Open Releases Page** when a newer version exists that this install may
    not take in place (``minimum_version``, or no artifact for this platform).
    """

    def __init__(self, info, frozen, parent=None):
        super().__init__(parent)
        self.info = info
        self.accepted_action = None      # "install" | "releases" | None
        self.setWindowTitle("Update Available")

        lay = QVBoxLayout(self)
        headline = QLabel(f"<b>XSLOPE Studio {info.latest} is available.</b><br>"
                          f"You are running {info.current}.")
        headline.setTextFormat(Qt.RichText)
        lay.addWidget(headline)

        if info.detail:
            note = QLabel(info.detail)
            note.setWordWrap(True)
            lay.addWidget(note)

        notes = QLabel(f'<a href="{info.notes_url}">Release notes</a>')
        notes.setTextFormat(Qt.RichText)
        notes.setOpenExternalLinks(True)      # QLabel hands the URL to the browser
        lay.addWidget(notes)

        pip_row = None
        if not frozen and info.status == "update-available":
            hint = QLabel("Update it the way you installed it:")
            hint.setWordWrap(True)
            lay.addWidget(hint)
            pip_row = QHBoxLayout()
            self.cmd_edit = QLineEdit(updater.PIP_UPGRADE_COMMAND)
            self.cmd_edit.setReadOnly(True)
            self.cmd_edit.setObjectName("pip_command")
            font = self.cmd_edit.font()
            font.setFamily("Menlo, Consolas, monospace")
            self.cmd_edit.setFont(font)
            self.copy_btn = QPushButton("Copy")
            self.copy_btn.setObjectName("copy_command")
            self.copy_btn.clicked.connect(self._copy_command)
            pip_row.addWidget(self.cmd_edit, 1)
            pip_row.addWidget(self.copy_btn)
            lay.addLayout(pip_row)

        buttons = QDialogButtonBox(self)
        self.install_btn = None
        self.releases_btn = None
        if frozen and info.status == "update-available":
            self.install_btn = buttons.addButton("Download && Install",
                                                 QDialogButtonBox.AcceptRole)
            self.install_btn.setObjectName("download_install")
            self.install_btn.setDefault(True)
        elif info.status == "update-manual":
            self.releases_btn = buttons.addButton("Open Releases Page",
                                                  QDialogButtonBox.AcceptRole)
            self.releases_btn.setObjectName("open_releases")
        later = buttons.addButton("Later", QDialogButtonBox.RejectRole)
        later.setObjectName("later")
        buttons.accepted.connect(self._on_accept)
        buttons.rejected.connect(self.reject)
        lay.addWidget(buttons)

    def _copy_command(self):
        cb = QApplication.clipboard()
        if cb is not None:
            cb.setText(updater.PIP_UPGRADE_COMMAND)
        self.copy_btn.setText("Copied")

    def _on_accept(self):
        self.accepted_action = "install" if self.install_btn is not None else "releases"
        self.accept()


class UpdateController(QObject):
    """Owns the check, the preference, the download and the install hand-off.

    The main window holds one of these and calls two methods on it; everything
    else — runner lifetimes, the throttle stamp, which dialog appears — lives
    here so the window's Help menu stays three lines long.
    """

    def __init__(self, window):
        super().__init__(window)
        self.win = window
        self.settings = window.settings
        self._runner = None
        self._manual = False
        self._dialog = None            # kept alive while non-modal
        self._info = None
        self._tmpdir = None

    # --- the preference --------------------------------------------------
    def startup_check_enabled(self):
        """Default ON: a classroom install that is never updated is the failure
        mode this feature exists to prevent."""
        v = self.settings.value(updater.SETTING_STARTUP_CHECK, True)
        return v if isinstance(v, bool) else str(v).lower() in ("true", "1")

    def set_startup_check_enabled(self, on):
        self.settings.setValue(updater.SETTING_STARTUP_CHECK, bool(on))

    def _stamp_check(self):
        self.settings.setValue(updater.SETTING_LAST_CHECK, updater.utc_stamp())

    # --- entry points ----------------------------------------------------
    def check_manual(self):
        """Help → Check for Updates…. Always ends in exactly one dialog."""
        if self._runner is not None and self._runner.isRunning():
            self.win.statusBar().showMessage("Already checking for updates…")
            return
        self._start(manual=True)

    def check_at_startup(self):
        """The silent check. Returns True if it actually started one."""
        if not self.startup_check_enabled():
            return False
        last = self.settings.value(updater.SETTING_LAST_CHECK, "")
        if not updater.due_for_check(last):
            return False
        self._start(manual=False)
        return True

    def _start(self, manual):
        self._manual = manual
        if manual:
            self.win.statusBar().showMessage("Checking for updates…")
            QApplication.setOverrideCursor(Qt.WaitCursor)
        self._runner = UpdateCheckRunner(parent=self)
        self._runner.checked.connect(self._on_checked)
        self._runner.finished.connect(self._on_check_finished)
        self._runner.start()

    # --- the answer ------------------------------------------------------
    def _on_checked(self, info):
        self._info = info
        # Stamp whatever the outcome was, including a failure: an offline machine
        # must not retry on every launch, and the manual item is always there.
        self._stamp_check()
        manual = self._manual

        if info.status == "error":
            if manual:
                QMessageBox.information(
                    self.win, "Check for Updates",
                    f"{info.detail}\n\nYou can download the newest version from\n"
                    f"{updater.RELEASES_URL}")
            self.win.statusBar().showMessage(
                "Update check failed." if manual else "")
            return

        if info.status == "up-to-date":
            if manual:
                QMessageBox.information(
                    self.win, "Check for Updates",
                    f"XSLOPE Studio {info.current} is up to date.")
            self.win.statusBar().showMessage(
                f"XSLOPE Studio {info.current} is up to date." if manual else "")
            return

        self.win.statusBar().showMessage(f"XSLOPE Studio {info.latest} is available.")
        self.show_update_dialog(info, modal=manual)

    def _on_check_finished(self):
        if self._manual:
            QApplication.restoreOverrideCursor()
        if self._runner is not None:
            self._runner.deleteLater()
            self._runner = None

    def show_update_dialog(self, info, modal=True):
        """The update notification. Modal from the menu, non-modal from startup."""
        frozen = updater.is_frozen()
        dlg = UpdateAvailableDialog(info, frozen, parent=self.win)
        self._dialog = dlg                       # a non-modal dialog needs an owner
        dlg.finished.connect(lambda *_: self._on_dialog_done(dlg))
        if modal:
            dlg.exec()
        else:
            dlg.setModal(False)
            dlg.setWindowModality(Qt.NonModal)
            dlg.show()
        return dlg

    def _on_dialog_done(self, dlg):
        action = dlg.accepted_action
        self._dialog = None
        if action == "releases":
            QDesktopServices.openUrl(QUrl(dlg.info.notes_url))
        elif action == "install":
            self.start_download(dlg.info)

    # --- download + install (frozen only) --------------------------------
    def start_download(self, info):
        """Fetch the platform artifact into a temp directory, with the status-bar
        progress bar and Cancel button a run uses."""
        artifact = info.artifact
        if not artifact:
            QMessageBox.warning(self.win, "Update",
                                "This release has no download for this platform.")
            return
        self._tmpdir = tempfile.mkdtemp(prefix="xslope-update-")
        runner = UpdateDownloadRunner(artifact, self._tmpdir,
                                      label=f"XSLOPE Studio {info.latest}",
                                      parent=self)
        self.win._update_dl_runner = runner
        runner.progress.connect(self.win._on_run_progress)
        runner.succeeded.connect(lambda path: self._on_downloaded(info, path))
        runner.failed.connect(self._on_download_failed)
        runner.cancelled.connect(self._on_download_cancelled)
        runner.finished.connect(self._on_download_finished)
        self.win.progress_bar.setRange(0, 100)
        self.win.progress_bar.setValue(0)
        self.win.progress_bar.setVisible(True)
        self.win.cancel_btn.setEnabled(True)
        self.win.cancel_btn.setVisible(True)
        self.win.statusBar().showMessage(f"Downloading XSLOPE Studio {info.latest}…")
        runner.start()

    def _on_download_failed(self, message):
        QMessageBox.critical(self.win, "Update failed", message)
        self.win.statusBar().showMessage("Update failed.")

    def _on_download_cancelled(self):
        self.win.statusBar().showMessage("Update cancelled.")

    def _on_download_finished(self):
        self.win.progress_bar.setVisible(False)
        self.win.progress_bar.setRange(0, 100)
        self.win.cancel_btn.setVisible(False)
        runner = getattr(self.win, "_update_dl_runner", None)
        if runner is not None:
            runner.deleteLater()
            self.win._update_dl_runner = None

    def _on_downloaded(self, info, path):
        """The artifact is on disk AND its sha256 matched. Hand off per platform."""
        try:
            plan = updater.plan_install(info, path)
        except Exception as exc:
            QMessageBox.critical(self.win, "Update failed", str(exc))
            return
        self.run_plan(plan)

    def run_plan(self, plan):
        """Execute a prepared :class:`~studio.updater.InstallPlan`.

        Windows quits first and installs after: the silent NSIS upgrade kills
        ``XSLOPE Studio.exe`` and wipes the directory we are executing from, so
        the helper is spawned only once the window has actually agreed to close
        (a cancelled save prompt cancels the update with it). macOS mounts the
        dmg and states the two steps; pip and the no-installer platforms only
        ever show text.
        """
        if plan.kind == "windows-installer":
            box = QMessageBox(self.win)
            box.setIcon(QMessageBox.Question)
            box.setWindowTitle("Install Update")
            box.setText(plan.message)
            ok = box.addButton("Install and Restart", QMessageBox.AcceptRole)
            box.addButton("Cancel", QMessageBox.RejectRole)
            box.setDefaultButton(ok)
            box.exec()
            if box.clickedButton() is not ok:
                self.win.statusBar().showMessage("Update cancelled.")
                return
            # close() runs the unsaved-changes prompt and returns False if the
            # user backed out — spawn the helper only after the app is committed
            # to exiting, or it would install underneath a still-running Studio.
            if not self.win.close():
                self.win.statusBar().showMessage("Update cancelled.")
                return
            updater.launch_detached(plan.argv)
            QTimer.singleShot(0, QApplication.quit)
            return

        if plan.kind in ("macos-dmg", "macos-replace"):
            updater.launch_detached(plan.argv)
            QMessageBox.information(self.win, "Finish Updating", plan.message)
            return

        QMessageBox.information(self.win, "Update", plan.message)
