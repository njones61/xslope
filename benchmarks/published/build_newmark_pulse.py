"""Build the rectangular-pulse acceleration record for Newmark's (1965) closed form.

Newmark's Rankine lecture opens with a rigid block on a base that is given one
rectangular pulse of acceleration A·g lasting t0. While the pulse lasts the block
gains velocity on the base at (A − k)·g, where k·g is its yield acceleration; when
the pulse ends the base stops accelerating and the block decelerates at k·g until
it comes to rest. The displacement it accumulates is

    D = (g t0² / 2) · A (A − k) / k

(Newmark's single-pulse result V²/(2gN) · (1 − N/A), with V = A g t0 and N = k).

The record is written as two columns of time (s) and acceleration (g), sampled
every 0.01 s. The pulse's end is one sample at t0 holding A and the next a
nanosecond later holding zero, so the record, taken as linear between samples,
is the rectangle to within that nanosecond. Its other sign is zero throughout,
so the block slides one way only.

Run from the repo root:
    python benchmarks/published/build_newmark_pulse.py
"""
import os

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..")
OUT = os.path.join(ROOT, "docs", "verification", "files", "published")

A = 0.3          # pulse amplitude, g
T0 = 0.5         # pulse duration, s
T_END = 3.0      # the record runs on at zero until the block has stopped
DT = 0.01


def main():
    rows = []
    n_pulse = round(T0 / DT)
    for i in range(n_pulse + 1):
        rows.append((i * DT, A))
    rows.append((T0 + 1e-9, 0.0))
    n_end = round(T_END / DT)
    for i in range(n_pulse + 1, n_end + 1):
        rows.append((i * DT, 0.0))
    lines = ["# Rectangular acceleration pulse (Newmark 1965): "
             f"{A:g} g for {T0:g} s, then zero to {T_END:g} s",
             "# time (s), acceleration (g)"]
    lines += [f"{t:.9f}\t{a:g}" for t, a in rows]
    path = os.path.join(OUT, "newmark_pulse.txt")
    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    return path


if __name__ == "__main__":
    print(main())
