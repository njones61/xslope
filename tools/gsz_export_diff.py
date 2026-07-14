"""What does the exporter LOSE? Import a real GeoStudio file, export it, and diff.

The reference here is **GeoStudio's own XML**, not our own reader. That is the whole point.
A round-trip through our own reader proves nothing -- reader and writer share the same
mental model of the schema, agree perfectly with each other, and pass while GeoStudio
rejects the file or quietly disables half its UI. Diffing our output against the vendor's
input gives an oracle that lives on this machine.

Usage:
    python tools/gsz_export_diff.py <folder-of-gsz-files> [--verbose]

For every element GeoStudio writes, this reports the child tags and attributes it ALWAYS
writes that we never do. Some of those are genuinely not ours to write (solved results,
search definitions, window scroll state). The ones that matter are the ones GeoStudio
needs to consider the model DEFINED -- and the only way to know which is which is to look
at the list, decide, and write the decision down. Hence _EXPECTED_ABSENT below: everything
not in it is a finding.

Vendor .gsz files are Seequent's and are not redistributable, so this takes a path to a
local, git-ignored folder. run_tests.py --gsz runs the same check against the synthetic
fixture we author ourselves.
"""

import argparse
import collections
import glob
import io
import os
import sys
import xml.etree.ElementTree as ET
import zipfile

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from xslope.geostudio import read_gsz, gsz_to_slope_data, export_gsz  # noqa: E402

# Things GeoStudio writes that xslope deliberately does not, with the reason. Anything
# NOT in here is reported as a finding -- silence has to be unreachable.
_EXPECTED_ABSENT = {
    # Solved results and the search that produced them: not part of the model.
    "SlipSurface": "the search definition — no xslope equivalent",
    "InitialInputInfo": "solver state",
    "ResultInputInfo": "solver state",
    "LambdaSettings": "solver setting",
    "Convergence": "solver setting",
    "UnderRelaxationCriteria": "solver setting (GeoStudio defaults it)",
    "SampleSettings": "probabilistic sampling — not written",
    "ParentID": "analysis inheritance — we write one analysis",
    # Pure presentation, which GeoStudio regenerates or defaults.
    "View": "display preferences",
    "ResultGraphs": "canned graph definitions — GeoStudio recreates them",
    "SketchItems": "annotations",
    "Window": "the saved scroll/zoom state — uncomputable, and guessing it framed the "
              "model badly; GeoStudio fits the view when it is absent",
    "PointLabel": "per-point label visibility",
    "Precision": "display precision",
    "SourceId": "geometry provenance",
    "Locations": "3D locations",
    "Transformation": "geometry transform",
    "Settings": "geometry settings",
    "MeshDefaultEdgeLength": "meshing — SLOPE/W does not mesh",
    "MeshConstraintsByGGID": "meshing",
    "MeshId": "meshing",
    "Mesh": "meshing",
    # Physics xslope does not export.
    "BCs": "boundary conditions — SEEP/W",
    "Hydraulic": "SEEP/W material properties",
    "SeepModel": "SEEP/W material model",
    "GeometryUsesHydraulicBCs": "SEEP/W boundary conditions",
    "ExcludeInitDeformation": "SIGMA/W",
    "Description": "analysis description",
}


def paths_of(root):
    """Every tag path and tag@attr in a document, as a set."""
    out = set()

    def walk(el, prefix):
        here = prefix + "/" + el.tag
        out.add(here)
        for a in el.attrib:
            out.add(f"{here}@{a}")
        for c in el:
            walk(c, here)

    walk(root, "")
    return out


def xml_of(path):
    zf = zipfile.ZipFile(path)
    name = next(n for n in zf.namelist() if n.endswith(".xml"))
    return ET.fromstring(zf.read(name))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("folder", help="folder of vendor .gsz files (git-ignored)")
    ap.add_argument("--verbose", action="store_true", help="list every lost path")
    args = ap.parse_args()

    files = sorted(glob.glob(os.path.join(args.folder, "*.gsz")))
    lost_counts = collections.Counter()
    lost_where = collections.defaultdict(set)
    scored = 0

    for src in files:
        try:
            gsz = read_gsz(src)
            if not gsz["analyses"]:
                continue
            aid = gsz["analyses"][0]["id"]
            slope_data, _ = gsz_to_slope_data(gsz, aid)
        except Exception as e:
            print(f"{os.path.basename(src):<48} import failed: {e}")
            continue
        if not slope_data.get("polygons"):
            print(f"{os.path.basename(src):<48} skipped (no polygons to export)")
            continue

        out = os.path.join(args.folder, "_roundtrip_tmp.gsz")
        try:
            export_gsz(slope_data, out, analysis_name="roundtrip", method="Spencer")
            theirs, ours = paths_of(xml_of(src)), paths_of(xml_of(out))
        except Exception as e:
            print(f"{os.path.basename(src):<48} export failed: {e}")
            continue
        finally:
            if os.path.exists(out):
                os.remove(out)

        scored += 1
        lost = theirs - ours
        for p in lost:
            lost_counts[p] += 1
            lost_where[p].add(os.path.basename(src))
        print(f"{os.path.basename(src):<48} they write {len(theirs):>4} paths, "
              f"we write {len(ours):>3}; we lose {len(lost):>4}")

    if not scored:
        sys.exit("nothing scored")

    # A path GeoStudio writes in EVERY file is one it always needs. Those are the findings.
    always = sorted(p for p, n in lost_counts.items() if n == scored)
    print(f"\n{len(always)} path(s) GeoStudio writes in ALL {scored} files but we never do.")

    findings = []
    for p in always:
        top = p.lstrip("/").split("/")[1] if p.count("/") > 1 else p.strip("/")
        seg = [s.split("@")[0] for s in p.lstrip("/").split("/")]
        if any(s in _EXPECTED_ABSENT and _EXPECTED_ABSENT[s] for s in seg):
            continue
        findings.append(p)

    if findings:
        print(f"\n*** {len(findings)} UNEXPLAINED — GeoStudio always writes these, we never "
              f"do, and no reason is recorded: ***")
        for p in findings:
            print("   ", p)
    else:
        print("\nEverything lost is accounted for in _EXPECTED_ABSENT.")

    if args.verbose:
        print("\nfull loss list (path: files):")
        for p in sorted(lost_counts, key=lambda k: (-lost_counts[k], k)):
            print(f"   [{lost_counts[p]}/{scored}] {p}")


if __name__ == "__main__":
    main()
