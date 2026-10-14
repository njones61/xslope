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


def use_bundled_certificates():
    """Point Python's https checks at the certificate list the bundle carries.

    Studio downloads with urllib (Open in Studio links, the update check, the AI
    providers' model lists), and urllib checks a site's certificate against
    OpenSSL's default certificate file. In the frozen app that default is a path on
    the BUILD machine (``/Library/Frameworks/Python.framework/...`` on macOS), which
    a user's Mac does not have, so every one of those downloads failed with
    "unable to get local issuer certificate" (1.0.1, on xslope.org). certifi's list
    ships inside the bundle; ``SSL_CERT_FILE`` makes it the default. On Windows the
    system certificate store is still read as well. A user's own ``SSL_CERT_FILE``
    is left alone, and a source checkout is not touched.
    """
    if not getattr(sys, "frozen", False) or os.environ.get("SSL_CERT_FILE"):
        return
    try:
        import certifi

        path = certifi.where()
    except Exception:
        return
    if os.path.isfile(path):
        os.environ["SSL_CERT_FILE"] = path


def use_persistent_matplotlib_cache():
    """Give Matplotlib a config folder that outlives the launch.

    PyInstaller's Matplotlib runtime hook points ``MPLCONFIGDIR`` at a new temporary
    folder on every start and deletes it at exit, so Matplotlib never finds its font
    list and rebuilds it each time. On macOS that rebuild runs ``system_profiler``
    over every installed font: 15 s on a warm start and over 40 s on a cold one,
    all before the window appears. A folder in the user's cache directory keeps the
    list from one launch to the next. Matplotlib reads the variable at import, so
    this runs before anything imports it. A source checkout is not touched.
    """
    if not getattr(sys, "frozen", False):
        return
    if sys.platform == "darwin":
        base = os.path.expanduser("~/Library/Caches")
    elif sys.platform.startswith("win"):
        base = os.environ.get("LOCALAPPDATA") or os.path.expanduser("~")
    else:
        base = os.environ.get("XDG_CACHE_HOME") or os.path.expanduser("~/.cache")
    path = os.path.join(base, "XSLOPE Studio", "matplotlib")
    try:
        os.makedirs(path, exist_ok=True)
    except OSError:
        return          # unwritable: keep the hook's temporary folder
    os.environ["MPLCONFIGDIR"] = path


use_bundled_certificates()
use_persistent_matplotlib_cache()

SELF_TEST_SAMPLE = "xslope_simple1.xlsx"

#: Fetched over https by the self-test through Open in Studio's own code. The Read
#: the Docs address is the one older links name: the fetch follows its forward to
#: xslope.org, so it checks the allowlist and both sites' certificate chains.
SELF_TEST_PACKAGE_URL = ("https://xslope.readthedocs.io/en/stable/lem/files/"
                         "xslope_simple_mult_layers.xslz")

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

        # The file must be RELEASED, not merely written. Reading a workbook used to
        # leave it open, which POSIX forgives -- an open file can still be replaced
        # and unlinked -- and Windows does not: the user cannot save over the file
        # they just opened, delete it, or open it in Excel while Studio is running.
        # Deleting it here states that as a check with a message on it, instead of
        # leaving it to surface as a PermissionError inside temp-directory cleanup.
        try:
            os.remove(out)
        except OSError as exc:
            emit(f"[fail] the saved file is still held open after loading it: {exc}")
            return 7
        emit("[ok] saved file released (no handle held after load)")

    # --- packaged skill prompt (the AI assistant reads it via importlib) ------
    from importlib import resources

    skill = (resources.files("xslope") / "resources" / "xslope_skill.md").read_text(
        encoding="utf-8")
    emit(f"[ok] packaged skill prompt, {len(skill)} chars")

    # --- the assistant's tokenizer (litellm counts tokens with tiktoken) -------
    # tiktoken finds cl100k_base through its tiktoken_ext plugin package, which the
    # 1.0.2 bundle left out: every assistant request failed. litellm points tiktoken
    # at the encoding files it ships, so this runs offline.
    try:
        from litellm.litellm_core_utils.default_encoding import encoding
        n = len(encoding.encode("slope stability"))
    except Exception as exc:
        emit(f"[fail] assistant tokenizer: {exc}")
        return 10
    emit(f"[ok] assistant tokenizer {encoding.name} ({n} tokens)")

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

    # --- https downloads (Open in Studio, the update check, model lists) -------
    # The frozen app must check certificates against the list it carries, not a
    # path on the build machine: a build runner that happens to have that path
    # would pass a live fetch while every user's install fails it.
    import ssl
    import tempfile

    import certifi

    cafile = ssl.get_default_verify_paths().cafile or ""
    n_ca = ssl.create_default_context().cert_store_stats().get("x509_ca", 0)
    if getattr(sys, "frozen", False) and (
            os.path.realpath(cafile) != os.path.realpath(certifi.where())):
        emit(f"[fail] https certificates come from {cafile or 'nowhere'!r}, not "
             f"from the bundled list {certifi.where()!r}")
        return 8
    if n_ca < 50:
        emit(f"[fail] only {n_ca} https root certificates are loaded")
        return 8
    emit(f"[ok] https certificates: {n_ca} roots from {os.path.basename(cafile)}")
    from studio import updater, urlscheme

    def _twice(fetch):
        try:
            return fetch()
        except Exception:
            return fetch()              # one retry for a passing network hiccup

    with tempfile.TemporaryDirectory() as tmp:
        pkg = _twice(lambda: urlscheme.download_package(SELF_TEST_PACKAGE_URL, tmp))
        emit(f"[ok] https download {os.path.basename(pkg)}, "
             f"{os.path.getsize(pkg)} bytes, via xslope.org")
    manifest = _twice(updater.fetch_manifest)
    emit(f"[ok] https update check: latest release {manifest.get('version')}")

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
