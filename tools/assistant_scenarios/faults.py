"""Planted faults — a sound model, broken on purpose, saved to a scratch copy.

A diagnosis scenario needs two things the real world does not hand you: a model
whose defects are known exactly, and a model that is otherwise RIGHT, so that
anything else the assistant accuses is a false positive by construction. Each
fault here is a transcription error of a kind that actually happens — a dropped
digit, a decimal slip, a deleted water definition, a starting circle left over
from another section — applied in memory to the loaded model and written to a copy
under the suite's scratch directory. The originals are never touched.
"""

from __future__ import annotations

import contextlib
import glob
import io
import os
import re
import shutil

from .core import Fault, load_model


def _material(slope_data, name):
    for mat in slope_data.get("materials") or []:
        if str(mat.get("name", "")).strip().lower() == name.strip().lower():
            return mat
    raise KeyError("no material named %r in the model" % name)


# --------------------------------------------------------------------------- #
# The faults
# --------------------------------------------------------------------------- #
def phi_typo(material, wrong, right):
    """A friction angle with a digit dropped (37 typed as 3)."""
    def apply(sd):
        _material(sd, material)["phi"] = float(wrong)
    return Fault(
        "phi_%s" % material.replace(" ", "_"),
        "material %r friction angle %s instead of %s" % (material, wrong, right),
        apply,
        # Naming the fault is naming the FIELD as wrong. Guessing the value that
        # was meant is a bonus, not the test: a session that read phi = 3 as a
        # transposed 33 found the fault, measured it, and led with it.
        re.compile(r"(?:phi|φ|friction angle)[^\n]{0,140}"
                   r"(?:\b%s\b|transpos|dropped|mis-?typ|typo|should be|"
                   r"is wrong|not credible|implausib)|"
                   r"\b%s\b[^\n]{0,140}(?:phi|φ|friction angle)"
                   % (int(right), int(right)), re.I))


def load_decimal_slip(factor, right, wrong):
    """A surface load an order of magnitude too large (240 psf typed as 2400)."""
    def apply(sd):
        for row in sd.get("dloads") or []:
            for point in row:
                if isinstance(point, dict) and point.get("Normal") is not None:
                    point["Normal"] = float(point["Normal"]) * float(factor)
    return Fault(
        "load_x%g" % factor,
        "the surface load is %s instead of %s" % (wrong, right),
        apply,
        re.compile(r"\b%s\b[^\n]{0,160}\b%s\b|\b%s\b[^\n]{0,160}\b%s\b"
                   % (wrong, right, right, wrong), re.I))


def max_depth_typo(wrong, right):
    """The bottom of the model an order of magnitude too deep."""
    def apply(sd):
        sd["max_depth"] = float(wrong)
    return Fault(
        "max_depth",
        "max_depth %s instead of %s" % (wrong, right),
        apply,
        re.compile(r"max_?depth[^\n]{0,160}(?:%s\b|dropped|typo|should be|"
                   r"order of magnitude)" % re.escape(str(right)), re.I))


def gamma_typo(material, wrong, right):
    """A unit weight with the decimal point in the wrong place."""
    def apply(sd):
        _material(sd, material)["gamma"] = float(wrong)
    return Fault(
        "gamma_%s" % material.replace(" ", "_"),
        "material %r unit weight %s instead of %s" % (material, wrong, right),
        apply,
        re.compile(r"(?:unit weight|gamma|γ)[^\n]{0,160}\b%s\b|"
                   r"\b%s\b[^\n]{0,160}(?:unit weight|gamma|γ)"
                   % (re.escape(str(wrong)), re.escape(str(wrong))), re.I))


def water_deleted(what="piezometric line"):
    """The water definition removed — every pore pressure with it."""
    def apply(sd):
        sd["piezo_line"] = []
        sd["piezo_line2"] = []
        for mat in sd.get("materials") or []:
            if str(mat.get("u", "")).lower() == "piezo":
                mat["u"] = "none"
    return Fault(
        "water_deleted",
        "the %s is gone (no pore pressure anywhere)" % what,
        apply,
        # Either half counts: the definition named as absent, or the consequence
        # named — the model is being analysed dry, with every material on u=none.
        re.compile(r"(?:piezometric|water table|phreatic|pore[- ]?(?:water )?"
                   r"pressure)[^\n]{0,160}(?:missing|none|gone|absent|blank|"
                   r"not defined|no longer|empty|removed|deleted)|"
                   r"analy[sz]ed dry|being analy[sz]ed dry|no pore[- ]?(?:water )?"
                   r"pressure|u\s*=\s*'?none'?[^\n]{0,80}(?:all|every|four|dry)",
                   re.I))


def circle_off_slope(index=0, dx=400.0):
    """A starting circle centered off the section, left over from another model."""
    def apply(sd):
        circles = sd.get("circles") or []
        if not circles:
            raise KeyError("the model has no starting circles")
        circle = circles[index]
        circle["Xo"] = float(circle.get("Xo", 0.0)) + float(dx)
    return Fault(
        "circle_off_slope",
        "starting circle %d is centered %g units off the section" % (index + 1, dx),
        apply,
        re.compile(r"(?:starting )?circle[^\n]{0,180}(?:off the|outside|"
                   r"does not|doesn'?t|beyond|away from|no longer|not on)", re.I))


# --------------------------------------------------------------------------- #
# Writing the broken copy
# --------------------------------------------------------------------------- #
def companions(path):
    """The sidecars a workbook is auto-discovered with (mesh, seepage, transient)."""
    base = os.path.splitext(path)[0]
    return sorted(set(glob.glob(base + "_*.json") + glob.glob(base + "_*.csv")))


def plant(src, faults, out_dir):
    """Write a copy of ``src`` under ``out_dir``, with every fault applied.

    Called with no faults it is simply the copy — which is how every scenario
    gets a sandbox of its own: a session that writes a file beside the project it
    opened (``generate_report`` does) writes it beside the copy.

    The copy carries the original's own sheets — it is written by copying the file
    and saving the edited model back into it — so the workbook the assistant opens
    differs from the sound one in exactly the planted fields and nothing else. The
    sidecars travel with it, or a seepage model would arrive with no field to read.
    """
    from xslope.fileio import save_slope_data_to_xlsx

    os.makedirs(out_dir, exist_ok=True)
    dest = os.path.join(out_dir, os.path.basename(src))
    shutil.copy2(src, dest)
    for extra in companions(src):
        shutil.copy2(extra, os.path.join(out_dir, os.path.basename(extra)))
    if not faults:
        # Nothing to plant: the copy is the original, byte for byte. Every
        # scenario opens a copy — so that a session can never write beside the
        # repository's own models — and a scenario with no faults must not pay a
        # load/save round trip for that, which could move something quietly.
        return dest
    slope_data = load_model(dest)
    if slope_data is None:
        raise RuntimeError("planting failed: %s does not load" % src)
    for fault in faults:
        if fault.apply is not None:
            fault.apply(slope_data)
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        save_slope_data_to_xlsx(slope_data, dest)
    return dest
