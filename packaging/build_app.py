"""Build the XSLOPE Studio desktop app for the current platform.

One command, five steps:

1. **verify** the build environment (every import the bundle needs, plus the
   PyInstaller/Cython toolchain), printing the versions it found;
2. **compile the fast kernel** (``setup_kernel.py build_ext --inplace``) — PyPI
   stays pure-Python, installers carry the compiled Mohr-Coulomb SSRM kernel;
3. **freeze** with ``packaging/studio.spec`` (onedir, windowed);
4. **smoke-test the frozen build**: run the shipped binary offscreen so it opens
   a bundled sample, meshes it with the bundled gmsh, solves it, and reports its
   version — a build that cannot do that is not shipped;
5. **package**: on macOS, assemble a drag-to-Applications ``.dmg`` with hdiutil
   and verify it mounts. On Windows the installer is built by NSIS in CI
   (``packaging/windows/installer.nsi``) from the directory produced here.

Usage (from a virtualenv that has ``xslope[gui,ai,fem]`` + pyinstaller + Cython)::

    python packaging/build_app.py                 # full build + smoke test + dmg
    python packaging/build_app.py --skip-kernel   # no compiler available
    python packaging/build_app.py --no-dmg        # skip the macOS disk image
    python packaging/build_app.py --smoke-only    # re-run the smoke test only
    python packaging/build_app.py --clean         # delete build outputs and exit

Outputs land in ``packaging/dist`` (artifacts) and ``packaging/build`` (scratch),
both git-ignored.
"""

from __future__ import annotations

import argparse
import glob
import os
import plistlib
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
PACKAGING = REPO / "packaging"
DIST = PACKAGING / "dist"
WORK = PACKAGING / "build"
APP_NAME = "XSLOPE Studio"

IS_MAC = sys.platform == "darwin"
IS_WIN = sys.platform == "win32"

# Every one of these must import in the build environment, because the frozen
# app imports it. Missing entries surface here rather than as a runtime crash in
# a 400 MB bundle.
REQUIRED_IMPORTS = [
    "xslope", "studio", "numpy", "scipy", "pandas", "matplotlib", "shapely",
    "openpyxl", "lxml", "tabulate", "gmsh", "PySide6", "litellm", "keyring",
    "ezdxf", "PyInstaller",
]


# ---------------------------------------------------------------------------
# small helpers
# ---------------------------------------------------------------------------
def step(title):
    print("\n" + "=" * 72)
    print(f"== {title}")
    print("=" * 72, flush=True)


def run(cmd, **kw):
    print("$ " + " ".join(str(c) for c in cmd), flush=True)
    return subprocess.run([str(c) for c in cmd], check=True, **kw)


def du(path):
    """Human-readable size of a file or directory tree."""
    path = Path(path)
    if path.is_file():
        n = path.stat().st_size
    else:
        n = sum(f.stat().st_size for f in path.rglob("*")
                if f.is_file() and not f.is_symlink())
    for unit in ("B", "KB", "MB", "GB"):
        if n < 1024 or unit == "GB":
            return f"{n:.1f} {unit}"
        n /= 1024.0


def version():
    sys.path.insert(0, str(REPO))
    from xslope._version import __version__
    return __version__


def app_paths():
    """(bundle root to ship, executable to run) for this platform."""
    if IS_MAC:
        bundle = DIST / f"{APP_NAME}.app"
        return bundle, bundle / "Contents" / "MacOS" / APP_NAME
    exe_name = f"{APP_NAME}.exe" if IS_WIN else APP_NAME
    bundle = DIST / APP_NAME
    return bundle, bundle / exe_name


# ---------------------------------------------------------------------------
# 1. environment
# ---------------------------------------------------------------------------
def verify_env():
    step("verify build environment")
    print(f"python      {sys.version.split()[0]}  ({sys.executable})")
    print(f"platform    {sys.platform} {os.uname().machine if hasattr(os, 'uname') else ''}")
    missing = []
    for name in REQUIRED_IMPORTS:
        try:
            mod = __import__(name)
            ver = getattr(mod, "__version__", getattr(mod, "GMSH_API_VERSION", ""))
            print(f"  ok  {name:<12} {ver}")
        except Exception as exc:
            missing.append(f"{name}: {exc}")
            print(f"  MISSING  {name}: {exc}")
    if missing:
        raise SystemExit(
            "Build environment incomplete. Install the app's dependencies:\n"
            "    pip install '.[gui,ai,fem,cad]' pyinstaller Cython\n  - "
            + "\n  - ".join(missing))
    print(f"\nxslope version {version()}")


# ---------------------------------------------------------------------------
# 2. compiled fast kernel
# ---------------------------------------------------------------------------
def kernel_artifacts():
    pats = ["_fem_kernel*.so", "_fem_kernel*.pyd"]
    out = []
    for p in pats:
        out += [Path(f) for f in glob.glob(str(REPO / "xslope" / p))]
    return out


def compile_kernel():
    """Compile ``xslope/_fem_kernel.pyx`` in place.

    Installers ship the compiled kernel while PyPI stays pure-Python, so this is
    a packaging step, not a library step. Built in place (next to the .pyx) so
    that PyInstaller — which analyses the repository checkout — collects it as a
    normal extension module of the ``xslope`` package. If the interpreter's
    ``xslope`` resolves somewhere else (a non-editable install), the artifact is
    copied there too, so the build works either way.
    """
    step("compile the fast kernel")
    try:
        import Cython  # noqa: F401
    except ImportError:
        raise SystemExit(
            "Cython is required to compile the fast kernel (pip install Cython).\n"
            "Pass --skip-kernel to build a pure-Python bundle instead.")

    for stale in kernel_artifacts():
        stale.unlink()
    run([sys.executable, "setup_kernel.py", "build_ext", "--inplace"], cwd=REPO)

    built = kernel_artifacts()
    if not built:
        raise SystemExit("setup_kernel.py reported success but produced no "
                         "_fem_kernel extension module.")
    print(f"built {built[0].name} ({du(built[0])})")

    # Mirror into an out-of-tree install of xslope, if that is what will be
    # analysed (site-packages install rather than the checkout).
    import xslope
    installed = Path(xslope.__file__).resolve().parent
    if installed != (REPO / "xslope").resolve():
        shutil.copy2(built[0], installed / built[0].name)
        print(f"copied into installed package at {installed}")

    # Prove it imports before spending ten minutes freezing it.
    out = subprocess.run(
        [sys.executable, "-c",
         "from xslope import _fem_kernel; print(_fem_kernel.__file__)"],
        cwd=REPO, capture_output=True, text=True)
    if out.returncode != 0:
        raise SystemExit(f"compiled kernel does not import:\n{out.stderr}")
    print(f"imports as {out.stdout.strip()}")


# ---------------------------------------------------------------------------
# 3. icons + version resources, then freeze
# ---------------------------------------------------------------------------
def make_icon():
    """Render the platform icon from studio/resources/icon.png.

    macOS uses ``sips``/``iconutil`` (both are part of the OS, no extra deps).
    Windows needs a .ico, which Pillow writes if it is available; without Pillow
    the build simply ships without a custom icon rather than failing.
    """
    src = REPO / "studio" / "resources" / "icon.png"
    if not src.exists():
        return None
    WORK.mkdir(parents=True, exist_ok=True)
    if IS_MAC:
        icns = WORK / "xslope.icns"
        iconset = WORK / "xslope.iconset"
        if iconset.exists():
            shutil.rmtree(iconset)
        iconset.mkdir()
        for size in (16, 32, 64, 128, 256, 512):
            for scale, suffix in ((1, ""), (2, "@2x")):
                px = size * scale
                run(["sips", "-z", px, px, src, "--out",
                     iconset / f"icon_{size}x{size}{suffix}.png"],
                    stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        run(["iconutil", "-c", "icns", iconset, "-o", icns])
        return icns
    if IS_WIN:
        ico = WORK / "xslope.ico"
        try:
            from PIL import Image
        except ImportError:
            print("Pillow not installed - building without a Windows icon.")
            return None
        img = Image.open(src)
        img.save(ico, sizes=[(16, 16), (32, 32), (48, 48), (64, 64),
                             (128, 128), (256, 256)])
        return ico
    return None


def make_win_version_file():
    """Write a PyInstaller version-resource file so the .exe carries the version
    in its Windows file properties (the NSIS installer stamps the same value)."""
    if not IS_WIN:
        return None
    v = version()
    parts = [int(x) for x in v.split(".")[:3]] + [0, 0, 0, 0]
    quad = ", ".join(str(x) for x in parts[:4])
    path = WORK / "win_version.txt"
    WORK.mkdir(parents=True, exist_ok=True)
    path.write_text(f"""VSVersionInfo(
  ffi=FixedFileInfo(
    filevers=({quad}), prodvers=({quad}),
    mask=0x3f, flags=0x0, OS=0x40004, fileType=0x1, subtype=0x0,
    date=(0, 0)),
  kids=[
    StringFileInfo([StringTable('040904B0', [
      StringStruct('CompanyName', 'Norman L. Jones'),
      StringStruct('FileDescription', '{APP_NAME}'),
      StringStruct('FileVersion', '{v}'),
      StringStruct('InternalName', 'xslope-studio'),
      StringStruct('LegalCopyright', 'Apache-2.0'),
      StringStruct('OriginalFilename', '{APP_NAME}.exe'),
      StringStruct('ProductName', '{APP_NAME}'),
      StringStruct('ProductVersion', '{v}')])]),
    VarFileInfo([VarStruct('Translation', [1033, 1200])])
  ]
)
""", encoding="utf-8")
    return path


def freeze(icon, win_version_file):
    step("freeze with PyInstaller")
    env = dict(os.environ)
    env["XSLOPE_APP_NAME"] = APP_NAME
    if icon:
        env["XSLOPE_ICON"] = str(icon)
    if win_version_file:
        env["XSLOPE_WIN_VERSION_FILE"] = str(win_version_file)
    t0 = time.time()
    run([sys.executable, "-m", "PyInstaller", "--noconfirm", "--clean",
         "--log-level", "WARN",
         "--distpath", DIST, "--workpath", WORK,
         PACKAGING / "studio.spec"], cwd=REPO, env=env)
    bundle, exe = app_paths()
    if not exe.exists():
        raise SystemExit(f"PyInstaller finished but {exe} is missing.")
    print(f"\nfroze in {time.time() - t0:.0f}s -> {bundle} ({du(bundle)})")


# ---------------------------------------------------------------------------
# 4. smoke test the frozen build
# ---------------------------------------------------------------------------
EXPECTED_MARKERS = [
    "[ok] packaged template",
    # A save, not only a read: the writer copies the packaged template, edits the
    # workbook XML and re-zips the archive, and only the first of those three is
    # exercised by reading. The re-zip once shelled out to `zip`, which Windows does
    # not have, and the frozen app's every Save failed there with a green smoke test.
    "[ok] save round-trip",
    "[ok] packaged skill prompt", "[ok] gmsh",
    "[ok] bishop FS", "[ok] window title", "[ok] self-test passed",
]


def _run_frozen(args, log, timeout=600):
    """Run the frozen binary, returning (returncode, transcript).

    The transcript is read from the log file the app writes, because a windowed
    build on Windows has no console attached to stdout.
    """
    _bundle, exe = app_paths()
    env = dict(os.environ)
    env["QT_QPA_PLATFORM"] = "offscreen"
    env["XSLOPE_SELFTEST_LOG"] = str(log)
    if log.exists():
        log.unlink()
    proc = subprocess.run([str(exe), *args], capture_output=True, text=True,
                          env=env, timeout=timeout)
    transcript = log.read_text(encoding="utf-8") if log.exists() else proc.stdout
    return proc, transcript


def smoke_test(require_kernel=True):
    step("smoke-test the frozen build")
    _bundle, exe = app_paths()
    if not exe.exists():
        raise SystemExit(f"nothing to test: {exe} does not exist")
    WORK.mkdir(parents=True, exist_ok=True)
    log = WORK / "smoke.log"

    proc, out = _run_frozen(["--version"], log, timeout=180)
    reported = out.strip().splitlines()[-1] if out.strip() else ""
    print(f"--version -> {reported!r} (exit {proc.returncode})")
    if proc.returncode != 0 or reported != version():
        print(proc.stdout[-4000:])
        print(proc.stderr[-4000:], file=sys.stderr)
        raise SystemExit(f"frozen app reported version {reported!r}, "
                         f"expected {version()!r}")

    proc, out = _run_frozen(["--self-test"], log)
    print(out)
    if proc.returncode != 0:
        print(proc.stderr[-8000:], file=sys.stderr)
        raise SystemExit(f"self-test failed (exit {proc.returncode})")
    expected = list(EXPECTED_MARKERS)
    if require_kernel:
        expected.append("[ok] compiled fast kernel")
    missing = [m for m in expected if m not in out]
    if missing:
        raise SystemExit("self-test transcript is missing: " + ", ".join(missing))
    print("smoke test PASSED")


# ---------------------------------------------------------------------------
# 5. macOS disk image
# ---------------------------------------------------------------------------
def make_dmg():
    step("assemble the .dmg")
    bundle, _exe = app_paths()
    v = version()
    dmg = DIST / f"XSLOPE-Studio-{v}-macos-{os.uname().machine}.dmg"
    if dmg.exists():
        dmg.unlink()

    staging = Path(tempfile.mkdtemp(prefix="xslope-dmg-"))
    try:
        # Drag-to-Applications layout: the app next to a symlink to /Applications.
        run(["cp", "-R", bundle, staging / bundle.name])
        os.symlink("/Applications", staging / "Applications")
        run(["hdiutil", "create", "-quiet", "-volname", f"{APP_NAME} {v}",
             "-srcfolder", staging, "-fs", "HFS+", "-format", "UDZO",
             "-ov", dmg])
    finally:
        shutil.rmtree(staging, ignore_errors=True)

    # Verify it actually mounts, and that the app is inside.
    mnt = Path(tempfile.mkdtemp(prefix="xslope-mnt-"))
    try:
        run(["hdiutil", "attach", "-quiet", "-nobrowse", "-readonly",
             "-mountpoint", mnt, dmg])
        contents = sorted(p.name for p in mnt.iterdir() if not p.name.startswith("."))
        print(f"mounted {dmg.name}: {contents}")
        if bundle.name not in contents or "Applications" not in contents:
            raise SystemExit("dmg does not have the drag-to-Applications layout")
    finally:
        subprocess.run(["hdiutil", "detach", "-quiet", str(mnt)],
                       stderr=subprocess.DEVNULL)
        shutil.rmtree(mnt, ignore_errors=True)
    print(f"dmg {dmg} ({du(dmg)})")
    return dmg


# ---------------------------------------------------------------------------
def clean():
    step("clean")
    for p in (DIST, WORK):
        if p.exists():
            shutil.rmtree(p)
            print(f"removed {p}")
    for f in kernel_artifacts() + [REPO / "xslope" / "_fem_kernel.c"]:
        if f.exists():
            f.unlink()
            print(f"removed {f}")
    b = REPO / "build" / "temp.kernel"
    if b.exists():
        shutil.rmtree(b, ignore_errors=True)


def report(dmg):
    step("build report")
    bundle, exe = app_paths()
    print(f"version    {version()}")
    print(f"bundle     {bundle}  ({du(bundle)})")
    if IS_MAC:
        plist = bundle / "Contents" / "Info.plist"
        with open(plist, "rb") as fh:
            info = plistlib.load(fh)
        print(f"Info.plist CFBundleShortVersionString="
              f"{info.get('CFBundleShortVersionString')} "
              f"CFBundleVersion={info.get('CFBundleVersion')} "
              f"id={info.get('CFBundleIdentifier')}")
    if dmg:
        print(f"dmg        {dmg}  ({du(dmg)})")
    print("\nartifacts are UNSIGNED: signing/notarization is applied in CI when "
          "the credentials are present (see .github/workflows/release-installers.yml).")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--skip-kernel", action="store_true",
                    help="do not compile the fast kernel (pure-Python bundle)")
    ap.add_argument("--no-dmg", action="store_true", help="skip the macOS disk image")
    ap.add_argument("--dmg-only", action="store_true",
                    help="only assemble the .dmg from an existing build (CI signs "
                         "the .app between the freeze and this step)")
    ap.add_argument("--smoke-only", action="store_true",
                    help="only re-run the smoke test against an existing build")
    ap.add_argument("--no-kernel-check", action="store_true",
                    help="do not require the compiled kernel in the smoke test")
    ap.add_argument("--clean", action="store_true",
                    help="delete packaging/dist, packaging/build and kernel artifacts")
    args = ap.parse_args(argv)

    if args.clean:
        clean()
        return 0
    want_kernel = not (args.skip_kernel or args.no_kernel_check)
    if args.smoke_only:
        smoke_test(require_kernel=want_kernel)
        return 0
    if args.dmg_only:
        if not IS_MAC:
            raise SystemExit("--dmg-only is macOS only")
        report(make_dmg())
        return 0

    verify_env()
    if args.skip_kernel:
        step("compile the fast kernel — SKIPPED (--skip-kernel)")
    else:
        compile_kernel()
    freeze(make_icon(), make_win_version_file())
    smoke_test(require_kernel=want_kernel)
    dmg = None
    if IS_MAC and not args.no_dmg:
        dmg = make_dmg()
    report(dmg)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
