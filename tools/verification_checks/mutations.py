#!/usr/bin/env python3
"""Mutation test for the verification-page checks.

A gate that certifies a wrong number is worse than no gate.  Each mutation
below plants one defect of a class the checks are supposed to catch, runs the
relevant check against the mutated page, and requires it to FAIL.  A mutation
that slips through is reported as MISSED and fails this suite.

The first block is the rs2 set the checks were originally hardened against.
The rest cover the check forms added when the other five pages were calibrated:
upper-bound percentages, scientific notation, thousands separators, share
lists, share rows, prefix-qualified quantities, the hedged-precision escape,
rounded tag restatements, scientific-notation tag values, and the structural
figure mode.

Usage: python -m tools.verification_checks.mutations
"""
import os
import shutil
import sys
import tempfile

from . import certify, deltas, figures, tags
from .pages import PAGES

#: (page, check, name, old, new).  `check` selects which check must catch it.
MUTATIONS = [
    # ------------------------------------------------------------------ rs2
    ("rs2", "deltas", "M1 last digit, matrix table",
     "| SSRM | 1.347 | 1.36 (−1.0%)", "| SSRM | 1.347 | 1.36 (−1.1%)"),
    ("rs2", "deltas", "M2 sign flip, summary row",
     "SSRM 0.986 vs RS2 SSRM 0.99 (−0.4%)", "SSRM 0.986 vs RS2 SSRM 0.99 (+0.4%)"),
    ("rs2", "deltas", "M3 last digit, summary row",
     "SSRM 1.280 vs RS2 SSRM 1.26 (+1.6%)", "SSRM 1.280 vs RS2 SSRM 1.26 (+1.7%)"),
    ("rs2", "deltas", "M4 prose delta",
     "1.881**, −1.0% on", "1.881**, −1.2% on"),
    ("rs2", "deltas", "M5 authority-column table",
     "| 2 | 1.13 | 1.194 (+5.7%)", "| 2 | 1.13 | 1.194 (+5.8%)"),
    ("rs2", "deltas", "M6 range bound",
     "1.24–1.27 (+0.8%)", "1.24–1.27 (+0.9%)"),
    ("rs2", "deltas", "M7 whitelisted cross-ref",
     "cf. [RS2-63](#rs2-63)\n+2.1%)", "cf. [RS2-63](#rs2-63)\n+2.4%)"),
    ("rs2", "deltas", "M8 XSLOPE value moved",
     "| SSRM | 1.347 |", "| SSRM | 1.357 |"),
    ("rs2", "deltas", "M10 absolute FS difference",
     "is worth +0.012 (the 1.006 above", "is worth +0.038 (the 1.006 above"),
    ("rs2", "deltas", "M11 hedged unsigned prose %",
     "runs about 20% low", "runs about 60% low"),
    ("rs2", "deltas", "M12 bound claim",
     "All six stages land within 1.6%", "All six stages land within 0.3%"),
    ("rs2", "deltas", "M13 unsigned integer % in a table cell",
     "| SSRM | 1.409 | RS2 SSRM 1.38 (+2.1%) |",
     "| SSRM | 1.409 (44%) | RS2 SSRM 1.38 (+2.1%) |"),
    ("rs2", "deltas", "M15 range end that matches no operand",
     "is +1.9 / +1.9 /\n+2.9%", "is +1.9 / +1.9 /\n+3.4%"),
    ("rs2", "deltas", "M16 percent whose operand the page never prints",
     "gravity turn-on reads **0.944**, −3.7% on that same 0.98",
     "gravity turn-on reads −3.7% on that same 0.98"),
    ("rs2", "deltas", "N21a whitelist-operand drift (RS2-22 load direction)",
     "1.534", "1.499", "all"),
    ("rs2", "deltas", "N21d whitelist-operand drift (RS2-62 Spencer)",
     "on the same mechanism (0.784)", "on the same mechanism (0.794)"),
    ("rs2", "figures", "M14 two-panel caption in the 'other' bucket",
     "![RS2-1: FEM inputs, mesh, max shear strain and displacement vectors at the "
     "critical SRF](images/RS2-1.png)",
     "![RS2-1: FEM model and maximum shear strain contours at the critical SRF]"
     "(images/RS2-1.png)"),
    ("rs2", "figures", "M17 four-panel caption on the inputs-only figure",
     "![RS2-50: shortened 4.2 m geotextile layers (vp089, Ta = 11.4 kN/m) — FEM inputs.",
     "![RS2-50: shortened 4.2 m geotextile layers (vp089, Ta = 11.4 kN/m) — FEM inputs, "
     "mesh, max shear strain and displacement vectors at the critical SRF."),

    # ------------------------------------------- upper-bound percentages ---
    ("seep", "deltas", "U1 operand drift breaks a `<` bound",
     "| Discharge q | 28.5961 | 28.5960 (<0.01%) | |",
     "| Discharge q | 28.6200 | 28.5960 (<0.01%) | |"),
    ("seep", "deltas", "U2 `<` bound tightened below the actual difference",
     "| Discharge q | 28.5961 | 28.5960 (<0.01%) | |",
     "| Discharge q | 28.5961 | 28.5960 (<0.0001%) | |"),

    # ------------------------------------------------ scientific notation ---
    ("rocscience_groundwater", "deltas", "S1 mantissa drift in a ×10ⁿ comparison",
     "| Q (m³/s per m) | 6.070×10⁻⁵ | 6.066×10⁻⁵ (+0.1%) |",
     "| Q (m³/s per m) | 6.120×10⁻⁵ | 6.066×10⁻⁵ (+0.1%) |"),
    ("rocscience_groundwater", "deltas", "S2 delta digit in a ×10ⁿ comparison",
     "| Q (m³/s per m) | 6.070×10⁻⁵ | 6.066×10⁻⁵ (+0.1%) |",
     "| Q (m³/s per m) | 6.070×10⁻⁵ | 6.066×10⁻⁵ (+0.3%) |"),

    # ------------------------------------------------ thousands separators ---
    ("geostudio", "deltas", "T1 comma-grouped operand drift",
     "| Sliding-mass weight, SLOPE/W's own circle | ≈ 56,020 | 56,127 (−0.2%) |",
     "| Sliding-mass weight, SLOPE/W's own circle | ≈ 55,020 | 56,127 (−0.2%) |"),
    ("geostudio", "deltas", "T2 comma-grouped delta digit",
     "| Sliding-mass weight, SLOPE/W's own circle | ≈ 56,020 | 56,127 (−0.2%) |",
     "| Sliding-mass weight, SLOPE/W's own circle | ≈ 56,020 | 56,127 (−0.4%) |"),

    # ------------------------------------------------------- share lists ---
    ("rocscience_groundwater", "deltas", "L1 share list loses its 'of u₀' tail",
     "within 0.13%, 0.31% and 0.14% of $u_0$",
     "within 0.13%, 0.31% and 0.14% of the analytical series"),
    ("rocscience_groundwater", "deltas", "L2 a real delta beside a share list",
     "within 0.13%, 0.31% and 0.14% of $u_0$",
     "within 0.13%, 0.31% and 0.14% of $u_0$; Bishop 1.351 vs 1.349 (+7.7%)"),

    # -------------------------------------------------------- share rows ---
    ("rocscience", "deltas", "R1 a declared value row is renamed",
     "| ΔFS over the A range (±15%) — **sweep result** |",
     "| FS change over the A range (±15%) — **sweep result** |"),
    ("rocscience", "deltas", "R2 the value-row exemption does not blanket its table",
     "0.944 (−1.5% against XSLOPE's simplified",
     "0.944 (−1.9% against XSLOPE's simplified"),

    # ------------------------------------------- prefix-qualified values ---
    ("rocscience", "deltas", "P1 a COV loses its qualifier",
     "on a COV 124% variable", "on a 124% variable"),
    ("geostudio", "deltas", "P2 a probability of failure loses its qualifier",
     "**PF 0 → 1.45%**", "**0 → 1.45%**"),

    # ------------------------------------------------ precision + hedging ---
    ("geostudio", "deltas", "H1 unhedged two-decimal signed delta",
     "| Clay fill, reinforced (imported, on SLOPE/W's own circles) | −0.3% |",
     "| Clay fill, reinforced (imported, on SLOPE/W's own circles) | −0.27% |"),
    ("rocscience", "deltas", "H2 unhedged integer signed delta",
     "brings the residuals from −42.7% / +83.7%",
     "brings the residuals from −43% / +83.7%"),

    # ------------------------------------------------- zero at zero digits ---
    ("ssrm", "deltas", "Z1 a zero percentage that is not zero",
     "| Spencer, circles forced tangent to the foundation base (false base circle) "
     "| 1.70 | the proprietary slip-circle program's **1.7** (0%) |",
     "| Spencer, circles forced tangent to the foundation base (false base circle) "
     "| 1.80 | the proprietary slip-circle program's **1.7** (0%) |"),

    # -------------------------------------------------------------- tags ---
    ("ssrm", "tags", "G1 rounded tag restatement is a different value",
     "equilibrium criterion) | 1.33 |", "equilibrium criterion) | 1.43 |"),
    ("geostudio", "tags", "G2 tagged value dropped from its section",
     "| Janbu | 1.330 | 1.233 |", "| Janbu | — | 1.233 |"),
    ("rocscience_groundwater", "tags", "G3 scientific-notation tag mantissa drifts",
     "6.070×10⁻⁵", "6.170×10⁻⁵", "all"),
    ("rocscience_groundwater", "tags", "G4 scientific-notation tag exponent drifts",
     "4.137×10⁻⁴", "4.137×10⁻³", "all"),

    # ----------------------------------------------------------- figures ---
    ("rocscience_groundwater", "figures", "F1 two-panel caption on a one-panel figure",
     "![gw017: steady total-head field vs Fig 19-5](images/gw017.png)",
     "![gw017: mesh and solved heads](images/gw017.png)"),
    ("rocscience_groundwater", "figures", "F2 unclassifiable caption, not exempt",
     "![gw001: mesh and solved heads](images/gw001.png)",
     "![gw001](images/gw001.png)"),
    ("ssrm", "figures", "F3 structural mode: referenced image does not exist",
     "![griffiths1_mesh.png](../fem/images/griffiths1_mesh.png)",
     "![griffiths1_mesh.png](../fem/images/griffiths1_mesh_missing.png)"),
    ("ssrm", "figures", "F4 structural mode: empty alt text",
     "![griffiths1_mesh.png](../fem/images/griffiths1_mesh.png)",
     "![](../fem/images/griffiths1_mesh.png)"),
]

#: Exemptions planted to prove that an exemption which never fires is itself a
#: failure — the mechanism that keeps the lists honest as pages are edited.
DEAD_EXEMPTIONS = [
    ("rs2", "bounds", ("+9.9", "a phrase that occurs nowhere")),
    ("rs2", "whitelist", ("+9.9", "a phrase that occurs nowhere", "1.1", "1.0")),
    ("ssrm", "tag_exempt", ("9.999", "a tag line that occurs nowhere")),
]

DOCS = os.path.join(certify.REPO, "docs")


def _stage(tmp, page, text):
    """Write the mutated page where its relative image paths still resolve."""
    ver = os.path.join(tmp, "verification")
    os.makedirs(ver, exist_ok=True)
    for name in os.listdir(DOCS):
        src = os.path.join(DOCS, name)
        dst = os.path.join(tmp, name)
        if name != "verification" and not os.path.exists(dst):
            os.symlink(src, dst)
    for name in os.listdir(os.path.join(DOCS, "verification")):
        if name.endswith(".md"):
            continue
        link = os.path.join(ver, name)
        if not os.path.exists(link):
            os.symlink(os.path.join(DOCS, "verification", name), link)
    path = os.path.join(ver, page + ".md")
    with open(path, "w") as fh:
        fh.write(text)
    return path


def _quiet(*a):
    pass


def _run(check, path, cfg):
    if check == "deltas":
        return deltas.run(path, cfg, report=_quiet)[0]
    if check == "tags":
        return tags.run(path, cfg, report=_quiet)
    return figures.run(path, cfg, report=_quiet)


def main():
    fails = []
    total = 0
    for mut in MUTATIONS:
        page, check, name, old, new = mut[:5]
        how = mut[5] if len(mut) > 5 else "first"
        total += 1
        base = open(certify.page_path(page)).read()
        if base.count(old) < 1:
            print(f"  ANCHOR  {name} ({page}) — mutation anchor not found")
            fails.append(name)
            continue
        text = base.replace(old, new) if how == "all" else base.replace(old, new, 1)
        tmp = tempfile.mkdtemp()
        try:
            n = _run(check, _stage(tmp, page, text), PAGES[page])
        finally:
            shutil.rmtree(tmp, ignore_errors=True)
        print(f"  {'CAUGHT ' if n else 'MISSED '} {name} ({page}/{check})")
        if not n:
            fails.append(name)

    for page, field, entry in DEAD_EXEMPTIONS:
        total += 1
        cfg = PAGES[page]
        original = list(getattr(cfg, field))
        setattr(cfg, field, original + [entry])
        try:
            check = "tags" if field == "tag_exempt" else "deltas"
            n = _run(check, certify.page_path(page), cfg)
        finally:
            setattr(cfg, field, original)
        name = f"dead {field} exemption ({page})"
        print(f"  {'CAUGHT ' if n else 'MISSED '} {name}")
        if not n:
            fails.append(name)

    print(f"mutations: {total}  missed: {len(fails)}")
    for f in fails:
        print(f"  MISSED: {f}")
    return len(fails)


if __name__ == "__main__":
    sys.exit(main())
