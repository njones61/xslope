"""Build docs/tutorials/files/anderson_dam_1989_243.txt, the acceleration record
used by the LEM-15 tutorial and shown on the Newmark docs page.

The record is the U.S. Geological Survey's corrected acceleration for the 1989
Loma Prieta earthquake at Anderson Dam, California (USGS National Strong Motion
Project station 1652), from the SMA-1 accelerograph downstream of the dam,
horizontal component 243 degrees. The USGS publishes it in its SMC format, in
cm/s², in the Loma Prieta data set:

    https://escweb.wr.usgs.gov/nsmp-data/data_sets/891018_1.html
    (corrected SMA-1 archive corr_sma-1.zip, file smdata/891018_1/sma-1/vol_2/2803c.smc)

This script converts that file to two columns of time in s and acceleration in g,
the form :func:`xslope.seismic.read_record` reads, with a header that names the
record and its source. Give the downloaded SMC file's path::

    python -I tools/build_anderson_dam_record.py /path/to/2803c.smc

The SMC layout: 11 text lines, 6 lines of integers (8 to a line, 10 characters
each), 10 lines of reals (5 to a line, 15 characters each), the comment lines the
integer header counts, then the samples, 8 to a line, 10 characters each. The
integer header gives the number of samples (line 3, first value) and the number
of comment lines (line 2, last value); the real header gives the sampling rate in
samples per second (line 1, second value).
"""

import os
import sys

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "docs",
                   "tutorials", "files", "anderson_dam_1989_243.txt")
G_CM_S2 = 980.665          # standard gravity, cm/s²

HEADER = """\
# Loma Prieta 1989, Anderson Dam downstream, 243 degrees
# The Loma Prieta, California, earthquake of October 17, 1989 (local time), recorded
# downstream of Anderson Dam by U.S. Geological Survey National Strong Motion Project
# station 1652, horizontal component 243 degrees. Corrected acceleration (Brady and
# Mork 1990, USGS Open-File Report 90-247), converted here from cm/s2 to g.
# Source: https://escweb.wr.usgs.gov/nsmp-data/data_sets/891018_1.html
# time (s), acceleration (g)
"""


def _fields(line, width):
    line = line.rstrip("\n")
    return [line[i:i + width] for i in range(0, len(line), width)
            if line[i:i + width].strip()]


def read_smc(path):
    """``(station_line, component_line, rate, samples_cm_s2)`` from an SMC file."""
    with open(path, encoding="latin-1") as fh:
        lines = fh.read().splitlines()
    text, ints, reals = lines[:11], lines[11:17], lines[17:27]
    if not text[0].lstrip().startswith("2 CORRECTED ACCELEROGRAM"):
        raise ValueError(f"{path}: not a corrected accelerogram ({text[0].strip()!r})")
    int_values = [int(v) for ln in ints for v in _fields(ln, 10)]
    real_values = [float(v) for ln in reals for v in _fields(ln, 15)]
    npts, n_comments = int_values[16], int_values[15]
    rate = real_values[1]
    data_lines = lines[27 + n_comments:]
    samples = [float(v) for ln in data_lines for v in _fields(ln, 10)]
    if len(samples) != npts:
        raise ValueError(f"{path}: {len(samples)} samples, the header states {npts}")
    station = next(ln.strip() for ln in text if ln.lstrip().startswith("station"))
    return station, rate, samples


def main(argv):
    if len(argv) != 2:
        print(__doc__)
        return 2
    station, rate, samples = read_smc(argv[1])
    if "Anderson Dam - Downstream" not in station or "243" not in station:
        raise SystemExit(f"{argv[1]} is not the downstream 243 component: {station!r}")
    dt = 1.0 / rate
    with open(OUT, "w") as fh:
        fh.write(HEADER)
        for i, a in enumerate(samples):
            fh.write(f"{i * dt:.3f}  {a / G_CM_S2:.6g}\n")
    peak = max(abs(a) for a in samples)
    print(f"wrote {os.path.relpath(OUT)}: {len(samples)} samples at {dt:g} s, "
          f"peak {peak:.1f} cm/s2 = {peak / G_CM_S2:.4f} g")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
