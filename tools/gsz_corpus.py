"""Score the GeoStudio importer against SLOPE/W's own answers, over a whole corpus.

A solved .gsz carries SLOPE/W's factor of safety for **every trial slip surface it
scored**, not just the critical one. That makes it a rig, not just a sample: import the
model, re-solve each of SLOPE/W's circles with xslope's matching method, and compare.
The slip surface, the method and the geometry are then all held fixed, so any
disagreement is the *import* -- a dropped load, a mis-read strength, a crack that should
not be there.

That is how the tension-crack bug was caught. GeoStudio keeps a switched-off crack's
geometry, so the importer applied a 2 m water-filled crack to a model that had none. It
was invisible in the phi=0 clay analyses (pore pressure cannot touch an undrained
strength) and showed up only as a -2% bias in the c'=0 sand ones. No single file would
have told you; the corpus did.

Usage:
    python tools/gsz_corpus.py <folder-of-gsz-files> [--slices 30] [--csv out.csv]

GeoStudio's example files are Seequent's copyright and are NOT redistributable, so they
live in a local, git-ignored folder and this script takes its path. It is a development
tool, not a test: run_tests.py --gsz covers the importer with a synthetic fixture we
author ourselves.
"""

import argparse
import glob
import os
import statistics
import sys
import warnings

warnings.filterwarnings("ignore")

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from xslope.geostudio import (read_gsz, gsz_to_slope_data, read_gsz_results,  # noqa: E402
                              gsz_seep_steps)
from xslope.slice import generate_slices  # noqa: E402
from xslope import solve  # noqa: E402

# GeoStudio's method name -> the xslope solver that means the same thing. Comparing
# xslope's Spencer to SLOPE/W's Morgenstern-Price is a ~1% error all by itself.
METHODS = {
    "Spencer":            lambda df: solve.spencer(df),
    "MorgensternPrice":   lambda df: solve.mprice(df, f_type="half_sine"),
    "Morgenstern-Price":  lambda df: solve.mprice(df, f_type="half_sine"),
    "Bishop":             lambda df: solve.bishop(df),
    "BishopSimplified":   lambda df: solve.bishop(df),
    "Ordinary":           lambda df: solve.oms(df),
    "Fellenius":          lambda df: solve.oms(df),
    "Janbu":              lambda df: solve.janbu(df),
    "CorpsEngineers1":    lambda df: solve.corps_engineers(df),
    "LoweKarafiath":      lambda df: solve.lowe_karafiath(df),
}


def score_analysis(gsz, analysis, slices):
    """xslope vs SLOPE/W on every circle SLOPE/W scored. Returns a result dict.

    A transient analysis -- one whose water comes from a SEEP/W run that marched in time
    -- is scored at EVERY saved time step, with the pore-pressure field of that step and
    the surfaces SLOPE/W solved at that step. Scoring one step against another's water
    would compare two different problems.
    """
    name, method = analysis["name"], (analysis["method"] or "")
    out = {"analysis": name, "method": method, "n": 0, "skipped": 0, "noncirc": 0,
           "median": None, "worst": None, "wmedian": None, "note": "", "steps": 1}

    solver = METHODS.get(method)
    if solver is None:
        out["note"] = f"no xslope equivalent of '{method}'"
        return out

    # No seepage field -> one state, imported as-is (step stays None).
    steps = [s["step"] for s in gsz_seep_steps(gsz, name)] or [None]
    out["steps"] = len(steps)

    errs, werrs = [], []
    for step in steps:
        try:
            slope_data, caveats = gsz_to_slope_data(gsz, analysis["id"],
                                                    critical_surface=False, step=step)
        except Exception as e:
            out["note"] = f"import failed: {e}"
            return out
        out["caveats"] = caveats

        rows = read_gsz_results(gsz, name, step)
        if not rows:
            continue

        for row in rows:
            if not row["circular"]:
                out["noncirc"] += 1       # a block/specified surface: not a circle at all
                continue
            circle = {"Xo": row["xo"], "Yo": row["yo"], "R": row["r"], "Depth": None,
                      "Y": None, "Option": "Depth", "Movement": "Left to Right"}
            try:
                ok, res = generate_slices(slope_data, circle=circle, num_slices=slices,
                                          debug=False)
                if not ok:
                    out["skipped"] += 1
                    continue
                df = res[0]
                # Does the mass we built weigh what SLOPE/W says it weighs? This checks
                # the imported polygons and unit weights with the solver taken out of the
                # way -- and it tells a genuine circle from one *fitted* to a block.
                if row["weight"] > 0:
                    w = float(df["w"].sum())
                    werrs.append(100.0 * (w - row["weight"]) / row["weight"])
                success, r = solver(df)
                if not success:
                    out["skipped"] += 1
                    continue
                errs.append(100.0 * (r["FS"] - row["fs"]) / row["fs"])
            except Exception:
                out["skipped"] += 1

    if not errs:
        if out["noncirc"]:
            out["note"] = f"{out['noncirc']} non-circular surfaces (not rebuildable)"
        else:
            out["note"] = "not solved (no trial surfaces saved)"
        return out

    out["n"] = len(errs)
    out["median"] = statistics.median(errs)
    out["worst"] = max(errs, key=abs)
    out["within1"] = 100.0 * sum(1 for e in errs if abs(e) < 1.0) / len(errs)
    out["wmedian"] = statistics.median(werrs) if werrs else None
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("folder", help="folder of .gsz files (git-ignored; not redistributable)")
    ap.add_argument("--slices", type=int, default=30, help="slice count (SLOPE/W uses 30)")
    ap.add_argument("--csv", help="also write the per-analysis table here")
    args = ap.parse_args()

    files = sorted(glob.glob(os.path.join(args.folder, "*.gsz")))
    if not files:
        sys.exit(f"no .gsz files in {args.folder}")

    print(f"{len(files)} file(s), {args.slices} slices, xslope method matched to SLOPE/W's")
    print("weight = xslope's sliding mass vs SLOPE/W's (geometry + unit weights alone);"
          " FS = the whole model\n")
    hdr = (f"{'file / analysis':<50} {'method':<15} {'n':>4} {'weight':>8} "
           f"{'FS':>8} {'worst':>8} {'<1%':>5}")
    print(hdr)
    print("-" * len(hdr))

    table = []
    for path in files:
        try:
            gsz = read_gsz(path)
        except Exception as e:
            print(f"{os.path.basename(path):<50} READ FAILED: {e}")
            continue
        print(f"{os.path.basename(path)}")
        for a in gsz["analyses"]:
            r = score_analysis(gsz, a, args.slices)
            r["file"] = os.path.basename(path)
            table.append(r)
            label = f"  {r['analysis'][:46]}"
            if r["n"]:
                flag = "" if abs(r["median"]) < 1.0 else "   <-- OFF"
                wm = f"{r['wmedian']:>+7.2f}%" if r["wmedian"] is not None else f"{'-':>8}"
                print(f"{label:<50} {r['method'][:14]:<15} {r['n']:>4} {wm} "
                      f"{r['median']:>+7.2f}% {r['worst']:>+7.2f}% {r['within1']:>4.0f}%{flag}")
            else:
                print(f"{label:<50} {r['method'][:14]:<15} {'-':>4} {r['note']}")

    scored = [r for r in table if r["n"]]
    if scored:
        meds = sorted(abs(r["median"]) for r in scored)
        good = sum(1 for r in scored if abs(r["median"]) < 1.0)
        print(f"\n{len(scored)} analyses scored over "
              f"{sum(r['n'] for r in scored)} slip surfaces")
        print(f"   |median error| < 1%: {good}/{len(scored)} analyses")
        print(f"   median of |median error|: {statistics.median(meds):.2f}%")
        bad = [r for r in scored if abs(r["median"]) >= 1.0]
        if bad:
            print("\n   still off:")
            for r in bad:
                print(f"      {r['file']} / {r['analysis']}: {r['median']:+.2f}%")

    if args.csv:
        import csv
        with open(args.csv, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["file", "analysis", "method", "n", "skipped",
                        "median_pct", "worst_pct", "note"])
            for r in table:
                w.writerow([r["file"], r["analysis"], r["method"], r["n"], r["skipped"],
                            f"{r['median']:.4f}" if r["median"] is not None else "",
                            f"{r['worst']:.4f}" if r["worst"] is not None else "",
                            r["note"]])
        print(f"\nwrote {args.csv}")


if __name__ == "__main__":
    main()
