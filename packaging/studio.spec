# -*- mode: python ; coding: utf-8 -*-
"""PyInstaller spec for the XSLOPE Studio desktop app (onedir, windowed).

Build it through ``packaging/build_app.py`` (which compiles the fast kernel,
renders the platform icon and smoke-tests the result) rather than calling
PyInstaller directly. If you do call it directly, run from the repository root:

    pyinstaller --noconfirm --distpath packaging/dist --workpath packaging/build \
        packaging/studio.spec

Onedir, not onefile: onefile unpacks ~400 MB to a temp directory on every
launch, and the auto-updater wants a stable, patchable install tree.

Environment variables honoured (all optional, all set by build_app.py):
    XSLOPE_ICON        path to the platform icon (.icns / .ico)
    XSLOPE_APP_NAME    override the app name (default "XSLOPE Studio")
"""

import os
import sys
from pathlib import Path

from PyInstaller.utils.hooks import (collect_data_files, collect_submodules,
                                     copy_metadata)

REPO = Path(SPECPATH).parent  # noqa: F821  (SPECPATH is injected by PyInstaller)
APP_NAME = os.environ.get("XSLOPE_APP_NAME", "XSLOPE Studio")
ICON = os.environ.get("XSLOPE_ICON") or None

sys.path.insert(0, str(REPO))
from xslope._version import __version__ as VERSION  # noqa: E402


def _log(msg):
    print(f"[studio.spec] {msg}")


# ---------------------------------------------------------------------------
# gmsh — the pip package is a single top-level ``gmsh.py`` plus a large shared
# library installed OUTSIDE site-packages (``<prefix>/lib/libgmsh.4.15.dylib``,
# ``<prefix>/Lib/gmsh-4.15.dll`` on Windows). PyInstaller's module graph only
# sees the .py wrapper, so the library has to be collected by hand.
#
# gmsh.py locates the library at import time and keeps the resolved path on the
# ctypes handle (``gmsh.lib._name``), so we ask the installed package where its
# library actually is instead of guessing. It searches ``<module dir>/<libname>``
# first, and for a frozen app the module dir is ``sys._MEIPASS`` — hence the
# library is placed at the bundle root.
# ---------------------------------------------------------------------------
def gmsh_binaries():
    import gmsh

    lib = getattr(getattr(gmsh, "lib", None), "_name", None)
    if lib and os.path.exists(lib):
        _log(f"gmsh library: {lib}")
        return [(str(lib), ".")]

    # Fallback: replicate gmsh.py's own search over the install prefix.
    if sys.platform == "win32":
        patterns = ["gmsh-*.dll", "gmsh.dll"]
        roots = [Path(sys.prefix) / d for d in ("Lib", "lib", "bin", ".")]
    elif sys.platform == "darwin":
        patterns = ["libgmsh.*.dylib", "libgmsh.dylib"]
        roots = [Path(sys.prefix) / d for d in ("lib", "bin", ".")]
    else:
        patterns = ["libgmsh.so.*", "libgmsh.so"]
        roots = [Path(sys.prefix) / d for d in ("lib", "lib64", "bin", ".")]
    for root in roots:
        for pat in patterns:
            hits = sorted(root.glob(pat))
            if hits:
                _log(f"gmsh library (fallback): {hits[-1]}")
                return [(str(hits[-1]), ".")]
    raise SystemExit(
        "studio.spec: could not locate the gmsh shared library. Install the "
        "gmsh wheel into the build environment (pip install 'xslope[fem]').")


# ---------------------------------------------------------------------------
# Compiled fast kernel — ``xslope/_fem_kernel*.{so,pyd}``. PyPI ships pure
# Python; the installers ship the compiled kernel. build_app.py compiles it
# in place (setup_kernel.py build_ext --inplace) before invoking PyInstaller,
# so here we only have to name it as a hidden import: PyInstaller then treats
# it as a normal extension module and places it inside the frozen xslope
# package, where ``from xslope import _fem_kernel`` finds it.
# ---------------------------------------------------------------------------
def kernel_hiddenimports():
    built = sorted((REPO / "xslope").glob("_fem_kernel*.so")) + \
        sorted((REPO / "xslope").glob("_fem_kernel*.pyd"))
    if built:
        _log(f"fast kernel: {built[-1].name}")
        return ["xslope._fem_kernel"]
    _log("fast kernel: NOT BUILT — bundling the pure-Python path only")
    return []


# ---------------------------------------------------------------------------
# Data files
# ---------------------------------------------------------------------------
def app_datas():
    datas = []

    # Packaged engine resources: the skill prompt the AI assistant loads and the
    # blank input template Studio copies on "Save As". Both are found at runtime
    # relative to the xslope package directory, so keep that layout.
    for name in sorted(os.listdir(REPO / "xslope" / "resources")):
        datas.append((str(REPO / "xslope" / "resources" / name), "xslope/resources"))

    # Studio's app icon (loaded from studio/resources/icon.png at startup).
    datas.append((str(REPO / "studio" / "resources" / "icon.png"), "studio/resources"))

    # Sample projects + the editable master template, so a user who has never
    # seen a .xlsx template has something to open on first launch. Companion
    # files (mesh / seepage) are discovered by naming convention next to the
    # .xlsx, so they must travel with it.
    samples = REPO / "docs" / "inputs" / "slope"
    for f in sorted(samples.iterdir()):
        if f.suffix.lower() in (".xlsx", ".json", ".csv") and not f.name.startswith("~$"):
            datas.append((str(f), "samples"))
    template = REPO / "docs" / "inputs" / "input_template.xlsx"
    if template.exists():
        datas.append((str(template), "samples"))

    # Third-party data files that are not covered by a PyInstaller hook.
    datas += collect_data_files("litellm")          # model price/context tables

    # importlib.metadata lookups performed at runtime.
    for dist in ("xslope", "litellm", "keyring", "tiktoken"):
        try:
            datas += copy_metadata(dist)
        except Exception as exc:  # pragma: no cover - optional dependency
            _log(f"metadata for {dist} unavailable ({exc})")
    return datas


hiddenimports = []
hiddenimports += collect_submodules("xslope")
hiddenimports += collect_submodules("studio")
hiddenimports += kernel_hiddenimports()
hiddenimports += [
    "gmsh",
    # Matplotlib: Qt only. The rest are excluded below.
    "matplotlib.backends.backend_qtagg",
    "matplotlib.backends.backend_agg",
    # openpyxl/lxml are imported at module scope by fileio but keep the writer's
    # optional paths visible to the graph.
    "openpyxl", "lxml", "lxml._elementpath", "lxml.etree",
    # AI assistant.
    "litellm", "keyring",
]
# keyring finds its backend through entry points; name the platform backend so
# the module graph keeps it.
if sys.platform == "darwin":
    hiddenimports += ["keyring.backends.macOS", "keyring.backends.macOS.api"]
elif sys.platform == "win32":
    hiddenimports += ["keyring.backends.Windows"]
else:
    hiddenimports += ["keyring.backends.SecretService", "keyring.backends.chainer"]
hiddenimports += collect_submodules("litellm")

excludes = [
    # No Tk anywhere: it drags in a whole Tcl/Tk runtime for nothing.
    "tkinter", "_tkinter", "Tkinter", "tcl", "tk",
    "matplotlib.backends.backend_tkagg", "matplotlib.backends.backend_tkcairo",
    "matplotlib.backends._backend_tk",
    "matplotlib.backends.backend_webagg", "matplotlib.backends.backend_webagg_core",
    "matplotlib.backends.backend_nbagg", "matplotlib.backends.backend_gtk3",
    "matplotlib.backends.backend_gtk4", "matplotlib.backends.backend_wx",
    "matplotlib.backends.backend_wxagg", "matplotlib.backends.backend_macosx",
    "matplotlib.testing",
    # One Qt binding only.
    "PyQt5", "PyQt6", "PySide2", "qtpy", "shiboken2",
    # Qt modules Studio never touches (each is tens of MB).
    "PySide6.QtWebEngineCore", "PySide6.QtWebEngineWidgets", "PySide6.QtWebEngineQuick",
    "PySide6.QtQuick", "PySide6.QtQuick3D", "PySide6.QtQuickWidgets", "PySide6.QtQml",
    "PySide6.Qt3DCore", "PySide6.Qt3DRender", "PySide6.Qt3DInput", "PySide6.Qt3DAnimation",
    "PySide6.Qt3DExtras", "PySide6.QtCharts", "PySide6.QtDataVisualization",
    "PySide6.QtMultimedia", "PySide6.QtMultimediaWidgets", "PySide6.QtBluetooth",
    "PySide6.QtNfc", "PySide6.QtPositioning", "PySide6.QtLocation",
    "PySide6.QtSerialPort", "PySide6.QtSerialBus", "PySide6.QtWebSockets",
    "PySide6.QtWebChannel", "PySide6.QtDesigner", "PySide6.QtHelp", "PySide6.QtTest",
    "PySide6.QtSql", "PySide6.QtSensors", "PySide6.QtRemoteObjects",
    "PySide6.QtScxml", "PySide6.QtSpatialAudio", "PySide6.QtTextToSpeech",
    "PySide6.QtUiTools", "PySide6.QtNetworkAuth", "PySide6.QtHttpServer",
    # Notebook / test / docs stacks that some scientific packages pull in.
    "IPython", "ipykernel", "jupyter", "jupyter_client", "jupyter_core", "notebook",
    "nbformat", "nbconvert", "pytest", "_pytest", "nose", "sphinx", "docutils",
    "wx", "PIL.ImageQt",
    # Never freeze the repo's own scripts/tests.
    "test", "tests", "xslope.tests",
]

a = Analysis(
    [str(REPO / "packaging" / "entry_studio.py")],
    pathex=[str(REPO)],
    binaries=gmsh_binaries(),
    datas=app_datas(),
    hiddenimports=sorted(set(hiddenimports)),
    hookspath=[],
    hooksconfig={
        # Matplotlib backends: collect only the Qt (and headless Agg) backends.
        "matplotlib": {"backends": ["QtAgg", "Agg"]},
    },
    runtime_hooks=[],
    excludes=excludes,
    noarchive=False,
    optimize=0,
)

pyz = PYZ(a.pure)  # noqa: F821

exe = EXE(  # noqa: F821
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name=APP_NAME,
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=False,
    console=False,              # windowed app
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,     # signing happens after the build, in CI
    entitlements_file=None,
    icon=ICON,
    version=os.environ.get("XSLOPE_WIN_VERSION_FILE") or None,
)

coll = COLLECT(  # noqa: F821
    exe,
    a.binaries,
    a.datas,
    strip=False,
    upx=False,
    upx_exclude=[],
    name=APP_NAME,
)

if sys.platform == "darwin":
    app = BUNDLE(  # noqa: F821
        coll,
        name=f"{APP_NAME}.app",
        icon=ICON,
        bundle_identifier="org.xslope.studio",
        version=VERSION,
        info_plist={
            "CFBundleName": APP_NAME,
            "CFBundleDisplayName": APP_NAME,
            "CFBundleShortVersionString": VERSION,
            "CFBundleVersion": VERSION,
            "NSHighResolutionCapable": True,
            "NSRequiresAquaSystemAppearance": False,
            "LSMinimumSystemVersion": "11.0",
            "NSHumanReadableCopyright": "Copyright (c) Norman L. Jones. Apache-2.0.",
            "CFBundleDocumentTypes": [
                {
                    "CFBundleTypeName": "XSLOPE project",
                    "CFBundleTypeRole": "Editor",
                    # Alternate, not Owner: Excel keeps the .xlsx association.
                    "LSHandlerRank": "Alternate",
                    "LSItemContentTypes": [
                        "org.openxmlformats.spreadsheetml.sheet"],
                }
            ],
        },
    )
