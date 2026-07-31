"""Frozen entry point for the XSLOPE Studio desktop app.

PyInstaller freezes THIS module, not ``studio/app.py``, so that the shipped
binary gains two packaging-owned behaviours the library entry point does not
need:

* ``--version`` — print ``xslope.__version__`` and exit. The auto-updater
  (Stage B) and CI both use this as a cheap probe of what a bundle contains,
  without starting Qt.
* ``--self-test [FILE]`` — headless smoke test of the *frozen* bundle: load a
  sample input, build a tri6 mesh (which exercises gmsh's bundled shared
  library), run one limit-equilibrium solve, then construct the real Studio
  main window offscreen and open the sample through it. Exits 0 only if every
  step succeeds.

Any other argument list is handed to ``studio.app.main`` unchanged, so the
installed app behaves exactly like ``xslope-studio``.
"""

from __future__ import annotations

import os
import sys

# Matplotlib must bind to PySide6 before any backend import (studio.app does the
# same; setting it here covers the --self-test path, which imports the engine
# before the GUI).
os.environ.setdefault("QT_API", "pyside6")

SELF_TEST_SAMPLE = "xslope_simple1.xlsx"

# A windowed build has no console on Windows, so stdout is not attached to the
# shell that launched it. Every diagnostic line therefore goes to stdout AND, if
# XSLOPE_SELFTEST_LOG names a file, to that file — which is what build_app.py and
# CI read back.
_LOG_PATH = os.environ.get("XSLOPE_SELFTEST_LOG")
_log_fh = None


def emit(msg):
    global _log_fh
    print(msg)
    sys.stdout.flush()
    if _LOG_PATH:
        if _log_fh is None:
            _log_fh = open(_LOG_PATH, "w", encoding="utf-8")
        _log_fh.write(str(msg) + "\n")
        _log_fh.flush()


def bundle_dir():
    """Directory holding the bundle's data files (``sys._MEIPASS`` when frozen,
    the repository root when running from a source checkout)."""
    if getattr(sys, "frozen", False):
        return getattr(sys, "_MEIPASS", os.path.dirname(sys.executable))
    return os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def samples_dir():
    """Directory holding the bundled sample inputs and the blank template."""
    frozen = os.path.join(bundle_dir(), "samples")
    if os.path.isdir(frozen):
        return frozen
    # Source checkout: the samples come straight from docs/.
    return os.path.join(bundle_dir(), "docs", "inputs", "slope")


def _version():
    from xslope._version import __version__

    return __version__


def _self_test(path=None):
    """Exercise the parts of the bundle that freezing most often breaks.

    Returns 0 on success. Every step prints a ``[ok]`` line so the build script
    can assert on the transcript rather than on the exit code alone.
    """
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

    if path is None:
        path = os.path.join(samples_dir(), SELF_TEST_SAMPLE)
    path = os.path.abspath(path)
    if not os.path.exists(path):
        emit(f"[fail] sample not found: {path}")
        return 2

    emit(f"[ok] version {_version()}")
    emit(f"[ok] frozen={bool(getattr(sys, 'frozen', False))} bundle={bundle_dir()}")
    emit(f"[ok] sample {path}")

    # --- engine: read an .xlsx (openpyxl + lxml + the packaged template) ------
    from xslope.fileio import load_slope_data

    slope_data = load_slope_data(path)
    emit(f"[ok] loaded {len(slope_data.get('materials', []))} material(s)")

    from xslope.fileio import default_template_path

    template = default_template_path()
    if not os.path.exists(template):
        emit(f"[fail] packaged template missing: {template}")
        return 3
    emit(f"[ok] packaged template {os.path.basename(template)}")

    # --- SAVE, not just read ------------------------------------------------
    # Reading the packaged template proves it is present and parseable. It does
    # not prove a save works, and a save is the more demanding path: it copies
    # the template, rewrites cells in the workbook XML, and re-zips the archive.
    # That last step used to shell out to the `zip` command, which does not exist
    # on Windows — so every Save in the frozen Windows app failed with
    # "[WinError 2] The system cannot find the file specified" while the build's
    # own smoke test stayed green. A round trip through the writer, inside the
    # bundle, is what would have caught it.
    import tempfile

    from xslope.fileio import save_slope_data_to_xlsx

    with tempfile.TemporaryDirectory() as td:
        out = os.path.join(td, "selftest_save.xlsx")
        save_slope_data_to_xlsx(slope_data, out)
        reread = load_slope_data(out)
        if len(reread.get("materials", [])) != len(slope_data.get("materials", [])):
            emit(f"[fail] saved file lost materials: "
                 f"{len(reread.get('materials', []))} of "
                 f"{len(slope_data.get('materials', []))}")
            return 6
        if abs(float(reread["gamma_water"]) - float(slope_data["gamma_water"])) > 1e-9:
            emit("[fail] saved file did not carry gamma_water")
            return 6
        # And in place, onto a file that already exists — the ordinary Save.
        save_slope_data_to_xlsx(reread, out)
        again = load_slope_data(out)
        if len(again.get("materials", [])) != len(slope_data.get("materials", [])):
            emit("[fail] in-place re-save lost materials")
            return 6
        emit(f"[ok] save round-trip {os.path.getsize(out)} bytes, "
             f"{len(again.get('materials', []))} material(s)")

    # --- packaged skill prompt (the AI assistant reads it via importlib) ------
    from importlib import resources

    skill = (resources.files("xslope") / "resources" / "xslope_skill.md").read_text(
        encoding="utf-8")
    emit(f"[ok] packaged skill prompt, {len(skill)} chars")

    # --- gmsh: build a tri6 mesh (loads gmsh's bundled shared library) --------
    import gmsh  # noqa: F401  (import is itself the test: it dlopens libgmsh)

    from xslope.mesh import build_mesh_from_polygons, get_material_polygons

    polygons = get_material_polygons(slope_data)
    mesh = build_mesh_from_polygons(polygons, 8.0, element_type="tri6")
    emit(f"[ok] gmsh {gmsh.GMSH_API_VERSION} tri6 mesh: "
         f"{len(mesh['nodes'])} nodes, {len(mesh['elements'])} elements")

    # --- one quick LEM solve (numpy/scipy/shapely/pandas path) ---------------
    from xslope.slice import generate_slices
    from xslope.solve import bishop

    circle = slope_data["circles"][0] if slope_data.get("circles") else None
    non_circ = None if circle else slope_data.get("non_circ")
    ok, result = generate_slices(slope_data, circle=circle, non_circ=non_circ,
                                 num_slices=20)
    if not ok:
        emit(f"[fail] generate_slices: {result}")
        return 4
    slice_df, _surface = result
    ok, res = bishop(slice_df)
    if not ok or not isinstance(res, dict):
        emit(f"[fail] bishop: {res}")
        return 5
    emit(f"[ok] bishop FS = {res['FS']:.4f} over {len(slice_df)} slices")

    # --- the compiled fast kernel (installers ship it, PyPI does not) --------
    try:
        from xslope import _fem_kernel

        emit(f"[ok] compiled fast kernel "
             f"{os.path.basename(getattr(_fem_kernel, '__file__', '?'))}")
    except ImportError as exc:
        emit(f"[warn] compiled fast kernel not bundled ({exc}); "
             f"SSRM runs use the pure-NumPy path")

    # --- the GUI itself, offscreen (PySide6 plugins + matplotlib Qt backend) --
    from PySide6 import QtCore
    from PySide6.QtWidgets import QApplication

    from studio.main_window import MainWindow

    app = QApplication.instance() or QApplication(sys.argv[:1])
    emit(f"[ok] Qt {QtCore.qVersion()} platform {app.platformName()}")
    win = MainWindow()
    win.show()
    win.open_path(path)
    app.processEvents()
    emit(f"[ok] window title {win.windowTitle()!r}")
    win.close()
    app.processEvents()

    emit("[ok] self-test passed")
    return 0


def main(argv=None):
    argv = list(sys.argv if argv is None else argv)
    args = argv[1:]

    if args and args[0] in ("--version", "-V"):
        emit(_version())
        return 0
    if args and args[0] == "--self-test":
        try:
            return _self_test(args[1] if len(args) > 1 else None)
        except Exception:
            import traceback
            emit("[fail] self-test raised:\n" + traceback.format_exc())
            return 9

    from studio.app import main as studio_main

    return studio_main(argv)


if __name__ == "__main__":
    raise SystemExit(main())
