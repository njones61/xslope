#!/usr/bin/env python3
"""At-failure capture admissibility — what may be drawn as a failure state.

The verification figures render their two right-hand panels from a SEPARATE
solve a margin beyond critical (``result['failure_solution']``), exported beside
the case as ``{stem}_fem_failure_*``. That capture is a moment inside a section
coming apart, so it is never an equilibrium — but it must at least BE the section
coming apart. A capture the finite guard stopped in its first sweep is an elastic
increment of an unsettled state: magnified by a deformation scale of a couple of
thousand it draws a picture of a mechanism that never formed, and on RJ-1d it
drew the blocks moving UP the slope while the row's text describes them toppling
out of it.

Three readings decide it, each with its threshold taken from the corpus's own
records (132 committed captures, read 2026-09-13):

``iterations``
    The sound captures run 3,000 iterations (RS2-51, the lowest) to 250,000; the
    three defective ones stopped at 2. ``MIN_ITERATIONS`` sits at 100 — more than
    an order above every defect and more than an order below every sound record,
    so it separates "never left the first sweep" from "ran".

``max_displacement``
    A capture BEYOND critical must be moving more than the last standing trial
    below it. Sound captures move 2.6x to 4 orders more; the three defects moved
    0.001x to 0.02x — less than the state they are supposed to have left. The
    ratio floor is therefore 1.0, and nothing in the corpus sits near it.

``vp_shear_strain``
    A numerically zero viscoplastic strain field is a field in which nothing
    yielded. It disqualifies a capture only on a row with NO joint state: ten
    sound joint captures have exactly zero viscoplastic strain because their
    blocks are elastic and the whole mechanism is slip on the joints, so on a
    jointed row the zero is the model, not a dead solve. The displacement reading
    above is what catches a dead jointed capture.

The producers (``benchmarks/rocscience/make_rs2_figures.py`` and its joint
sibling) import :func:`refusal_reason` so the state they REFUSE to draw and the
state this module refuses to accept as committed are one rule. ``--audit`` on
either producer scans the committed sidecars through :func:`scan_stem`.

Usage: python -m tools.verification_checks.captures <stem> [<stem> ...]
"""
import csv
import json
import os
import sys

#: Iterations below which a capture never left its first sweep (see module doc).
MIN_ITERATIONS = 100

#: A capture must move at least this multiple of the last standing field's
#: max|u|. 1.0 = "moving at all more than the state it left".
MIN_DISPLACEMENT_RATIO = 1.0


def refusal_reason(failure, standing=None, strain_max=None, has_joint_state=False):
    """Why this at-failure capture must not be drawn, or None if it may be.

    ``failure`` / ``standing`` are readings of the two fields — either the meta
    dicts written beside the case or the live solution dicts, since both carry
    ``iterations`` and ``max_displacement``. ``strain_max`` is the largest
    viscoplastic shear strain in the captured field when it is known (None = not
    read, not zero). ``has_joint_state`` says the row fails on joints, which is
    what makes a zero strain field legitimate.
    """
    if not failure:
        return None                      # no capture at all is a different case
    at = failure.get('capture_truncated_at')
    it = failure.get('iterations')
    if at is not None and at <= MIN_ITERATIONS:
        kind = failure.get('capture_truncated_kind') or 'runaway'
        return (f'capture stopped at iteration {at} ({kind}), at or below the '
                f'{MIN_ITERATIONS}-iteration floor')
    if it is None or it < MIN_ITERATIONS:
        return (f'capture ran {it} iterations, below the {MIN_ITERATIONS}-iteration '
                f'floor')
    u_f = failure.get('max_displacement')
    u_s = (standing or {}).get('max_displacement')
    if u_f is not None and u_s:
        if u_f < MIN_DISPLACEMENT_RATIO * u_s:
            return (f'capture max|u| = {u_f:.3e} is below the last standing '
                    f'field\'s {u_s:.3e} — the state it is meant to have left')
    if strain_max is not None and strain_max == 0 and not has_joint_state:
        return 'viscoplastic shear-strain field is numerically zero (nothing yielded)'
    return None


def _read_rows(path):
    """CSV rows of a sidecar, past its leading ``# units:`` comment lines."""
    with open(path) as fh:
        return list(csv.DictReader([l for l in fh if not l.startswith('#')]))


def _max_abs(path, column):
    try:
        rows = _read_rows(path)
    except OSError:
        return None
    vals = [abs(float(r[column])) for r in rows
            if r.get(column) not in (None, '')]
    return max(vals) if vals else None


def _json(path):
    try:
        with open(path) as fh:
            return json.load(fh)
    except OSError:
        return None


def scan_stem(stem):
    """Refusal reason for the capture committed at ``stem``, or None.

    None also covers a stem with no at-failure capture at all — a row whose
    figure stands on its converged field is not this check's business.
    """
    failure = _json(f'{stem}_fem_failure_meta.json')
    if failure is None:
        return None
    standing = _json(f'{stem}_fem_meta.json')
    strain_max = _max_abs(f'{stem}_fem_failure_elements.csv', 'vp_shear_strain')
    joints = f'{stem}_fem_failure_joints.csv'
    has_joints = os.path.exists(joints) and bool(_read_rows(joints))
    return refusal_reason(failure, standing, strain_max, has_joints)


def scan(stems, report=print):
    """Scan many stems; returns [(stem, reason)] for every inadmissible capture."""
    bad = []
    for stem in stems:
        why = scan_stem(stem)
        if why:
            bad.append((stem, why))
            report(f'  {os.path.basename(stem):22s} {why}')
    return bad


if __name__ == '__main__':
    sys.exit(1 if scan(sys.argv[1:]) else 0)
