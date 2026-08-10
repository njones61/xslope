# Copyright 2025 Norman L. Jones
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Regenerate the report-page images ``docs/studio/reports.md`` embeds.

The page has to show what a report LOOKS like, and the only honest way to get
that is to generate one and photograph it. So a real report is built from a
shipped sample — the ``gsat`` levee, which carries a saved seepage solution
beside it, so the document comes out with both a Seepage section and a Limit
Equilibrium one — laid out by LibreOffice, and the wanted pages rasterized.

The pages are found by the text on them (``PAGES`` below), never by number: a
sentence added anywhere above moves every page after it, and an image captioned
"the slice table" showing whatever landed on page 12 is worse than no image.

LibreOffice lays the document out here because it can be driven headlessly;
Word's own typesetting of the same file differs in small ways. What the images
are for is the shape of the page — the running heads, the numbered captions, the
cross-references, the landscape table — and that is the same in both.

Needs LibreOffice and poppler (``pdftoppm``/``pdftotext``); exits 0 with a note
where either is missing, so a checkout without them still runs the tool suite.

Run:  python tools/capture_report_pages.py
"""

from __future__ import annotations

import glob
import os
import shutil
import subprocess
import sys
import tempfile

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)

import matplotlib

matplotlib.use("Agg")

OUT_DIR = os.path.join(REPO_ROOT, "docs", "studio", "images")

#: The model the sample report is built from: a levee with a saved seepage
#: solution beside it, so one document carries both engines.
MODEL = os.path.join(REPO_ROOT, "docs/lem/files/xslope_gsat_seep.xlsx")

#: The title-page metadata, matching the Generate Report dialog's own capture so
#: the two images read as one session.
META = {"title": "Riverside Dam — Downstream Slope",
        "project_number": "2026-114",
        "organization": "Example Engineering",
        "author": "N. L. Jones",
        "date": "August 9, 2026"}

#: The method the report documents in full.
METHOD = "spencer"

#: Which pages are wanted, as ``(file name, text that is on that page and on no
#: page above it)``.
PAGES = (
    # Not a heading: every heading is on the contents page too, which stands
    # above the page that carries it.
    ("reports_page_body.png", "digest identifies"),
    ("reports_page_slice_table.png", "Slice data"),
)

#: Rasterizing resolution, in dpi. Enough that the body text is readable at the
#: width the page embeds it, without shipping a megabyte per image.
RENDER_DPI = 110

#: Where LibreOffice is looked for, in order.
SOFFICE = (os.environ.get("XSLOPE_SOFFICE"), "soffice",
           "/Applications/LibreOffice.app/Contents/MacOS/soffice")


def _tool(names):
    """The first of ``names`` that exists on this machine, or None."""
    for name in names:
        if not name:
            continue
        found = name if os.path.isabs(name) and os.path.exists(name) \
            else shutil.which(name)
        if found:
            return found
    return None


def _report(path):
    """Generate the sample report to ``path``."""
    from xslope.fileio import load_slope_data
    from xslope.report import generate_report, solutions_from_sidecars
    from xslope.slice import generate_slices
    from xslope.solve import solve_selected

    sd = load_slope_data(MODEL)
    solutions = solutions_from_sidecars(MODEL, sd)
    ok, made = generate_slices(sd, circle=sd["circles"][0], num_slices=20)
    if not ok:
        raise RuntimeError(f"the sample would not slice: {made}")
    slice_df, surface = made
    results = solve_selected(METHOD, slice_df)
    if not isinstance(results, dict):
        raise RuntimeError(f"the sample would not solve: {results}")
    solutions["lem"] = [{"slice_df": slice_df, "failure_surface": surface,
                         "results": results, "search": None, "method": METHOD}]
    options = dict(META, input_path=MODEL, method=METHOD)
    ok, out = generate_report(sd, solutions, options, path)
    if not ok:
        raise RuntimeError(out)
    return out


def _page_of(pdf, text):
    """The 1-based number of the first page of ``pdf`` carrying ``text``."""
    pages = int(subprocess.run(["pdfinfo", pdf], capture_output=True, text=True)
                .stdout.split("Pages:")[1].split()[0])
    for n in range(1, pages + 1):
        found = subprocess.run(["pdftotext", "-f", str(n), "-l", str(n), pdf, "-"],
                               capture_output=True, text=True).stdout
        if text in found:
            return n
    raise RuntimeError(f"no page of the report carries {text!r}")


def main():
    soffice = _tool(SOFFICE)
    poppler = all(_tool((name,)) for name in ("pdftoppm", "pdftotext", "pdfinfo"))
    if soffice is None or not poppler:
        print("capture_report_pages: LibreOffice or poppler is missing — skipped.")
        return 0

    print("capture_report_pages: regenerating the report page images")
    with tempfile.TemporaryDirectory(prefix="xslope_report_pages_") as tmp:
        docx = os.path.join(tmp, "sample_report.docx")
        out = _report(docx)
        print(f"  built a {len(out['figures'])}-figure report")
        subprocess.run([soffice, "--headless", "--convert-to", "pdf",
                        "--outdir", tmp, docx],
                       check=True, capture_output=True, timeout=600)
        pdf = os.path.splitext(docx)[0] + ".pdf"
        for name, marker in PAGES:
            n = _page_of(pdf, marker)
            stem = os.path.join(tmp, os.path.splitext(name)[0])
            subprocess.run(["pdftoppm", "-r", str(RENDER_DPI), "-png",
                            "-f", str(n), "-l", str(n), "-singlefile", pdf, stem],
                           check=True)
            written = glob.glob(stem + ".png")[0]
            shutil.copyfile(written, os.path.join(OUT_DIR, name))
            print(f"  wrote docs/studio/images/{name} (page {n})")
    print("done")
    return 0


if __name__ == "__main__":
    sys.exit(main())
