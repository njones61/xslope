"""VP104 — the Newmark permanent displacement of Slide2 verification problem #104.

Slide2 verification problem #104 is built on the vendor's own *Tutorial 28 — Seismic
Analysis with the Newmark Method*, and its fourth scenario reports a permanent
sliding displacement: the rigid block above the critical surface slides whenever the
ground acceleration exceeds the slope's yield acceleration K_y, and its relative
velocity is integrated to the displacement the slope accumulates over the record.

XSLOPE computes this with :func:`xslope.seismic.newmark_displacement`, and the
corpus tag ``VP104-newmark`` on docs/verification/rocscience.md locks it. This
script prints the comparison behind that tag: the displacement at each of Slide2's
two published yield accelerations, the yield acceleration each published
displacement implies, both directions of the record, and XSLOPE's own K_y carried
through.

The acceleration record
-----------------------
The record is not printed in the verification manual, but it ships with the
product: the tutorial model ``Tutorial 28 Seismic.slmd`` (a ZIP of one ``.slsc``
per scenario) stores it inside the "Newmark Displacement" scenario as parallel
blocks of times and accelerations (``m_newmark_user_times``,
``m_newmark_user_accelerations``), with ``m_example_earthquake_text:
"Mammoth Lakes-1 1980: CVK-090"`` alongside: 5862 samples at dt = 0.005 s, peaks
-0.41648 g and +0.341544 g.

The model is Rocscience's file and is not in this repository. Rocscience publishes
it in its deprecated-tutorials zip (``MODEL_URL``). The tags and this script read
the record straight out of the model, found by name under the local vendor archive
(``$XSLOPE_VENDOR_FILES``, default ``~/python_projects/vendor_files``), or from a
path given on the command line::

    python3 benchmarks/rocscience/newmark_vp104.py "/path/to/Tutorial 28 Seismic.slmd"

Without the model the tags skip and this script stops with a message saying where
to get it.

Published reference values (Slide2 verification manual, Table 104.1)
-------------------------------------------------------------------
=====================  ==========  ====================
Search                 K_y         Newmark displacement
=====================  ==========  ====================
uni-modal              0.140       5.081 cm
multi-modal (MMO)      0.139       5.042 cm
=====================  ==========  ====================

The pair is non-monotone in K_y (the larger yield acceleration is printed with the
larger displacement), which no integration of one record can produce, since the
displacement falls as K_y rises. The two come from two critical surfaces whose K_y
is printed to three decimals, and the run output inverts each published
displacement back to the K_y that gives it.
"""

import os
import sys
import zipfile

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..")
sys.path.insert(0, ROOT)

from xslope.seismic import newmark_displacement, _record_out  # noqa: E402

#: Root of the local archive of vendor files, which are never committed here.
VENDOR_FILES = os.environ.get(
    "XSLOPE_VENDOR_FILES", os.path.expanduser("~/python_projects/vendor_files"))
#: The Slide2 tutorial model that carries the record.
MODEL_NAME = "Tutorial 28 Seismic.slmd"
#: Rocscience's public download that holds the model.
MODEL_URL = ("https://static.rocscience.cloud/assets/model-files/Slide2/"
             "Slide2-Deprecated-Tutorials-v9027.zip")

# Slide2 verification manual Table 104.1 — {label: (K_y, displacement in cm)}.
SLIDE2_PUBLISHED = {
    "uni-modal": (0.140, 5.081),
    "multi-modal (MMO)": (0.139, 5.042),
}

# XSLOPE's own yield acceleration for this slope, the lower of its circular and
# non-circular searches on docs/verification/files/rocscience/vp104a.xlsx: the
# VP104-kc-noncirc lock (critical_kc, Spencer). The circles give 0.1445 (VP104-kc).
XSLOPE_KY = 0.1436


# --- the record, from the vendor model ---------------------------------------
def find_model(root=None):
    """The path of ``Tutorial 28 Seismic.slmd`` under the vendor archive, or None."""
    for dirpath, _dirs, files in os.walk(root or VENDOR_FILES):
        if MODEL_NAME in files:
            return os.path.join(dirpath, MODEL_NAME)
    return None


def load_record(slmd_path=None):
    """VP104's record, in the form :func:`xslope.seismic.read_record` returns,
    read from the vendor model (found with :func:`find_model` when no path is
    given). None when the model is not there."""
    path = slmd_path or find_model()
    if path is None:
        return None
    times, accel, description = extract_record(path)
    return _record_out(path, "Slide2", [description], times, accel,
                       description=description)


def write_record(path, slmd_path=None):
    """Write VP104's record as a two-column text file, for a tool that needs the
    record on disk. Returns ``path``, or None when the model is not there."""
    rec = load_record(slmd_path)
    if rec is None:
        return None
    with open(path, "w") as fh:
        fh.write("# %s\n" % rec["description"])
        for t, a in zip(rec["time"], rec["accel"]):
            fh.write("%.4f  %.8g\n" % (t, a))
    return path


def extract_record(slmd_path):
    """Pull the acceleration record out of a Slide2 ``.slmd`` tutorial model.

    Returns ``(times_s, accel_g, description)`` from the scenario that carries a
    POPULATED ``m_newmark_user_accelerations`` block — every scenario in the
    tutorial declares the block, but only "Newmark Displacement" fills it. The
    description is the scenario's ``m_example_earthquake_text``. The ``.slmd`` is
    a ZIP holding one ``.slsc`` text file per scenario; the block is a ``num:``
    count followed by ``index: value`` lines between ``<name> start:`` and
    ``<name> end:``.
    """
    with zipfile.ZipFile(slmd_path) as zf:
        for name in zf.namelist():
            if not name.lower().endswith(".slsc"):
                continue
            lines = [ln.rstrip("\r")
                     for ln in zf.read(name).decode("latin-1").split("\n")]
            if "  m_newmark_user_accelerations start:" not in lines:
                continue
            accel = _read_block(lines, "m_newmark_user_accelerations")
            if not accel:
                continue
            text = [ln.split(":", 1)[1].strip().strip('"') for ln in lines
                    if ln.strip().startswith("m_example_earthquake_text:")]
            return (_read_block(lines, "m_newmark_user_times"), accel,
                    text[0] if text else MODEL_NAME)
    raise ValueError("%s carries no Newmark acceleration record."
                     % os.path.basename(slmd_path))


def _read_block(lines, name):
    i = lines.index("  %s start:" % name)
    j = lines.index("  %s end:" % name)
    return [float(ln.split(":", 1)[1]) for ln in lines[i + 1:j]
            if not ln.strip().startswith("num")]


def yield_acceleration_for(record, displacement_cm):
    """The K_y at which the governing direction gives this displacement."""
    lo, hi = 0.001, 0.5
    for _ in range(60):
        mid = (lo + hi) / 2.0
        if newmark_displacement(record, mid)["D"] > displacement_cm:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2.0


def main(argv):
    rec = load_record(argv[1] if len(argv) > 1 else None)
    if rec is None:
        print("%s is not under %s. It is in Rocscience's deprecated-tutorials zip,\n"
              "%s;\nunzip it there or give the model's path as the argument."
              % (MODEL_NAME, VENDOR_FILES, MODEL_URL))
        return 1
    print("VP104 — Newmark permanent displacement")
    print("=" * 78)
    print("record   : %s, from %s" % (rec["description"], rec["path"]))
    print("           %d samples, dt = %.4f s, duration = %.3f s, peak %.5g g"
          % (rec["npts"], rec["dt"], rec["duration"], rec["pga"]))

    print()
    print("At the published yield accelerations (the larger direction governs)")
    print("-" * 78)
    for label, (ky, published) in sorted(SLIDE2_PUBLISHED.items()):
        r = newmark_displacement(rec, ky)
        print("%-20s K_y %.3f  %7.3f cm  vs Slide2 %.3f cm  %+5.1f%%   "
              "(other direction %.3f cm)"
              % (label, ky, r["D"], published,
                 100.0 * (r["D"] - published) / published,
                 r["directions"][-r["governing"]]["D"]))

    print()
    print("The K_y each published displacement implies")
    print("-" * 78)
    for label, (ky, published) in sorted(SLIDE2_PUBLISHED.items()):
        implied = yield_acceleration_for(rec, published)
        print("%-20s printed K_y %.3f  ->  %.3f cm implies K_y = %.5f  (%s)"
              % (label, ky, published, implied,
                 "within the printed rounding"
                 if abs(implied - ky) <= 0.0005 else
                 "outside the printed rounding"))

    print()
    r = newmark_displacement(rec, XSLOPE_KY)
    print("XSLOPE's own K_y = %.4f (VP104-kc-noncirc) -> %.3f cm" % (XSLOPE_KY, r["D"]))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
