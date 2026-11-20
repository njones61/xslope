"""Newmark records: every way a record file is written reads to the same record,
and the integration does what a rigid block does to it.

The corpus tags read records as two columns of time and acceleration. A PEER
.AT2 file never reaches them, so this check writes one record as two columns,
as the NGA-West2 file the PEER database serves, and again in the older PEER
layout, and requires the same samples and the same displacement from each. The
record is made here: a few sine waves under a rise-and-decay envelope, peaking
near 0.4 g and stronger one way than the other, so both directions slide by
different amounts. It also writes the run-together negative
values some fixed-width files carry, and the files the reader must refuse: a
record in cm/s², and times that run backward.

Two properties of the integration close it. Reversing the record's sign swaps
the two directions, which is the whole reason both are run. And the
displacement scales with g and nothing else between unit systems: the same
record and k_y give D in inches that is the D in cm times 32.174/9.81 × 12/100.
"""
import os
import sys
import tempfile

import numpy as np

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, ROOT)

from xslope.seismic import (newmark_displacement, read_record,  # noqa: E402
                            scale_for_pga)

KY = 0.139
DT = 0.005


def _synthetic_accel():
    """The test record's accelerations, in g, at DT for 20 s."""
    t = np.arange(0.0, 20.0 + DT / 2, DT)
    env = np.minimum(t / 2.0, 1.0) * np.exp(-np.maximum(t - 8.0, 0.0) / 4.0)
    a = (0.25 * np.sin(2 * np.pi * 1.3 * t)
         + 0.12 * np.sin(2 * np.pi * 2.9 * t + 1.0)
         + 0.08 * np.sin(2 * np.pi * 0.7 * t + 0.4)
         + 0.06 * np.sin(2 * np.pi * 0.35 * t) ** 2)
    return t, np.round(env * a, 8)


def _record_file(folder):
    """The test record as a two-column text file in ``folder``."""
    t, a = _synthetic_accel()
    return _write(folder, "test_record.txt", "# A synthetic test record\n" + "".join(
        f"{ti:.3f}  {ai:.8f}\n" for ti, ai in zip(t, a)))


def _record():
    with tempfile.TemporaryDirectory() as tmp:
        return read_record(_record_file(tmp))


def _write(folder, name, text):
    path = os.path.join(folder, name)
    with open(path, "w") as fh:
        fh.write(text)
    return path


def _at2_lines(accel, per_line=5, fmt="{:15.7E}"):
    rows = []
    for i in range(0, len(accel), per_line):
        rows.append("".join(fmt.format(a) for a in accel[i:i + per_line]))
    return rows


def check_at2_matches_the_text_file():
    fails = []
    text = _record()
    accel = list(text["accel"])
    with tempfile.TemporaryDirectory() as tmp:
        new = _write(tmp, "TEST090.AT2", "\n".join(
            ["PEER NGA STRONG MOTION DATABASE RECORD",
             "Test Event, 1/1/2000, Test Station, 90",
             "ACCELERATION TIME SERIES IN UNITS OF G",
             f"NPTS= {len(accel):5d}, DT=   .0050 SEC"]
            + _at2_lines(accel)) + "\n")
        old = _write(tmp, "TEST090_old.at2", "\n".join(
            ["PACIFIC EARTHQUAKE ENGINEERING RESEARCH CENTER",
             "TEST EVENT 01/01/00, TEST STATION, 090",
             "ACCELERATION TIME HISTORY IN UNITS OF G",
             f"  {len(accel)}    .00500    NPTS, DT"]
            + _at2_lines(accel, per_line=8, fmt="{:10.6f}")) + "\n")
        d_text = newmark_displacement(text, KY)["D"]
        for path, label, exact in ((new, "NGA-West2 .AT2", True),
                                   (old, "older PEER .AT2", False)):
            rec = read_record(path)
            if rec["format"] != "AT2":
                fails.append(f"{label}: read as {rec['format']}")
            if rec["npts"] != text["npts"] or rec["dt"] != 0.005:
                fails.append(f"{label}: {rec['npts']} points at dt {rec['dt']}")
            worst = float(np.max(np.abs(rec["accel"] - text["accel"])))
            # the older layout's six decimals round each value by up to 5e-7
            if worst > (1e-7 if exact else 5.01e-7):
                fails.append(f"{label}: accelerations differ by up to {worst:.2g} g")
            if "Test Station" not in rec["description"] and \
                    "TEST STATION" not in rec["description"]:
                fails.append(f"{label}: description {rec['description']!r}")
            d = newmark_displacement(rec, KY)["D"]
            if abs(d - d_text) > (1e-6 if exact else 1e-3):
                fails.append(f"{label}: D = {d:.6f} cm against the text file's "
                             f"{d_text:.6f}")
    return fails


def check_run_together_values():
    with tempfile.TemporaryDirectory() as tmp:
        path = _write(tmp, "packed.at2", "\n".join(
            ["A FIXED-WIDTH FILE", "STATION X", "UNITS OF G",
             "NPTS=     6, DT=   .0100 SEC",
             "-0.1234E-02-0.5000E-01 0.2500E+00",
             "-0.1000E+00-.2000D+00 0.3000E-00"]) + "\n")
        rec = read_record(path)
    want = [-0.001234, -0.05, 0.25, -0.1, -0.2, 0.3]
    if not np.allclose(rec["accel"], want):
        return [f"run-together values read as {list(rec['accel'])}"]
    return []


def check_refusals():
    fails = []
    with tempfile.TemporaryDirectory() as tmp:
        cms = _write(tmp, "cms.at2", "\n".join(
            ["PEER", "STATION", "ACCELERATION TIME SERIES IN UNITS OF CM/S/S",
             "NPTS=     2, DT=   .0100 SEC", "1.0 2.0"]) + "\n")
        back = _write(tmp, "back.txt", "0.00 0.1\n0.02 0.2\n0.01 0.3\n")
        short = _write(tmp, "short.at2", "\n".join(
            ["PEER", "STATION", "UNITS OF G",
             "NPTS=     5, DT=   .0100 SEC", "1.0 2.0"]) + "\n")
        for path, what in ((cms, "a record in cm/s²"),
                           (back, "times that run backward"),
                           (short, "fewer values than the header states")):
            try:
                read_record(path)
            except ValueError:
                continue
            fails.append(f"{what} was read")
    return fails


def check_uneven_steps():
    with tempfile.TemporaryDirectory() as tmp:
        path = _write(tmp, "uneven.txt", "# a pulse\ntime, accel\n"
                      "0.0, 0.3\n0.5, 0.3\n0.500000001, 0\n5.0, 0\n")
        rec = read_record(path)
    fails = []
    if rec["dt"] is not None:
        fails.append(f"uneven steps reported dt = {rec['dt']}")
    if rec["header"] != ["a pulse", "time, accel"]:
        fails.append(f"header read as {rec['header']}")
    return fails


def check_reversal_swaps_the_directions():
    rec = _record()
    fwd = newmark_displacement(rec, KY)
    rev = newmark_displacement(dict(rec, accel=-rec["accel"]), KY)
    fails = []
    for sign in (+1, -1):
        a, b = fwd["directions"][sign]["D"], rev["directions"][-sign]["D"]
        if abs(a - b) > 1e-12:
            fails.append(f"direction {sign:+d}: {a} reversed reads {b}")
    if fwd["D"] != rev["D"] or fwd["governing"] != -rev["governing"]:
        fails.append("reversing the record changed the governing D")
    return fails


def check_units_and_scale():
    rec = _record()
    si = newmark_displacement(rec, KY, units="si")
    us = newmark_displacement(rec, KY, units="imperial")
    want = si["D"] / 100.0 * 32.174 / 9.81 * 12.0
    fails = []
    if us["D_unit"] != "in" or abs(us["D"] - want) > 1e-9 * want:
        fails.append(f"imperial D {us['D']} {us['D_unit']}, expected {want} in")
    s = scale_for_pga(rec, 0.3)
    scaled = newmark_displacement(rec, KY, scale=s)
    if abs(scaled["pga"] - 0.3) > 1e-12:
        fails.append(f"scaled to 0.3 g, the peak reads {scaled['pga']}")
    return fails


def run():
    """Failures as a list, for run_tests.py."""
    fails = []
    for check in (check_at2_matches_the_text_file, check_run_together_values,
                  check_refusals, check_uneven_steps,
                  check_reversal_swaps_the_directions, check_units_and_scale):
        try:
            fails += check()
        except Exception as exc:          # a crash is a failure, not a pass
            fails.append(f"{check.__name__} raised {exc!r}")
    return fails


if __name__ == "__main__":
    failures = run()
    if failures:
        print("FAILURES:")
        for f in failures:
            print(f"  - {f}")
        raise SystemExit(1)
    print("Newmark record checks passed.")
