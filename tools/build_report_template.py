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

"""Build the Analysis Report's default Word template.

The template is a binary that ships in the wheel, so it is built by a script
rather than edited by hand: every style, size and margin in it is readable here
and the file can be reproduced exactly.

It declares the LOOK and nothing else — page size and margins, the Title,
Subtitle, Heading, Body Text and Caption styles, and empty header and footer
frames for :mod:`xslope.report_docx` to fill. It carries no content: a report is
built onto it, never inserted into it. A company template put in its place
restyles the whole report without a line of code changing.

The sizes are deliberately restrained. A submittal is read, not presented: body
text at 10.5 pt, a first-level heading at 14 pt, and one-inch margins.

    python tools/build_report_template.py
"""

import os
import sys

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)

from docx import Document                                          # noqa: E402
from docx.enum.style import WD_STYLE_TYPE                          # noqa: E402
from docx.enum.text import WD_ALIGN_PARAGRAPH                      # noqa: E402
from docx.shared import Inches, Pt, RGBColor                       # noqa: E402

OUT = os.path.join(REPO_ROOT, "xslope", "resources", "report_template.docx")

#: The body typeface. Calibri is present on every Windows and Office install and
#: substitutes cleanly elsewhere; the template names it once, here.
BODY_FONT = "Calibri"

#: Heading colour — a dark slate that prints black-and-white without turning grey.
HEADING_RGB = RGBColor(0x1F, 0x3B, 0x57)


def _style(doc, name, kind=WD_STYLE_TYPE.PARAGRAPH, base=None):
    """A style by name, added to the template when it is not already there."""
    try:
        return doc.styles[name]
    except KeyError:
        style = doc.styles.add_style(name, kind)
        if base is not None:
            style.base_style = doc.styles[base]
        return style


def build(path=OUT):
    doc = Document()

    # --- page setup ---
    for section in doc.sections:
        section.page_width, section.page_height = Inches(8.5), Inches(11)
        section.top_margin = section.bottom_margin = Inches(1.0)
        section.left_margin = section.right_margin = Inches(1.0)
        section.header_distance = section.footer_distance = Inches(0.5)
        # The renderer fills these; creating them here puts the header and footer
        # parts in the template, so a custom template is a drop-in replacement.
        section.different_first_page_header_footer = True
        section.header.is_linked_to_previous = False
        section.footer.is_linked_to_previous = False

    # --- body text ---
    normal = doc.styles["Normal"]
    normal.font.name = BODY_FONT
    normal.font.size = Pt(10.5)
    normal.font.color.rgb = RGBColor(0x00, 0x00, 0x00)
    normal.paragraph_format.space_after = Pt(8)
    normal.paragraph_format.line_spacing = 1.08

    body = _style(doc, "Body Text")
    body.font.name = BODY_FONT
    body.font.size = Pt(10.5)
    body.paragraph_format.space_after = Pt(8)
    body.paragraph_format.line_spacing = 1.08

    # --- headings: restrained, and they stay with the text beneath them ---
    for level, size, bold, italic, before in ((1, 14, True, False, 16),
                                              (2, 12, True, False, 12),
                                              (3, 11, True, True, 10)):
        style = _style(doc, f"Heading {level}")
        style.font.name = BODY_FONT
        style.font.size = Pt(size)
        style.font.bold = bold
        style.font.italic = italic
        style.font.color.rgb = HEADING_RGB
        style.paragraph_format.space_before = Pt(before)
        style.paragraph_format.space_after = Pt(4)
        style.paragraph_format.keep_with_next = True

    # --- title page ---
    title = _style(doc, "Title")
    title.font.name = BODY_FONT
    title.font.size = Pt(22)
    title.font.bold = True
    title.font.color.rgb = HEADING_RGB
    title.paragraph_format.space_before = Pt(6)
    title.paragraph_format.space_after = Pt(4)
    title.paragraph_format.alignment = WD_ALIGN_PARAGRAPH.CENTER

    subtitle = _style(doc, "Subtitle")
    subtitle.font.name = BODY_FONT
    subtitle.font.size = Pt(12.5)
    subtitle.font.bold = False
    subtitle.font.italic = False
    subtitle.font.color.rgb = RGBColor(0x44, 0x44, 0x44)
    subtitle.paragraph_format.space_after = Pt(24)
    subtitle.paragraph_format.alignment = WD_ALIGN_PARAGRAPH.CENTER

    # --- captions ---
    caption = _style(doc, "Caption")
    caption.font.name = BODY_FONT
    caption.font.size = Pt(9)
    caption.font.bold = False
    caption.font.italic = True
    caption.font.color.rgb = RGBColor(0x33, 0x33, 0x33)
    caption.paragraph_format.space_before = Pt(4)
    caption.paragraph_format.space_after = Pt(8)
    caption.paragraph_format.keep_with_next = True

    # --- header / footer text ---
    for name in ("Header", "Footer"):
        style = _style(doc, name)
        style.font.name = BODY_FONT
        style.font.size = Pt(9)
        style.font.color.rgb = RGBColor(0x44, 0x44, 0x44)

    bullet = _style(doc, "List Bullet")
    bullet.font.name = BODY_FONT
    bullet.font.size = Pt(10.5)
    bullet.paragraph_format.space_after = Pt(3)

    # Instantiate the table style so a report built on this template finds it
    # whether or not the host Word install would have supplied it latently.
    doc.add_table(rows=1, cols=1).style = doc.styles["Table Grid"]

    # The template ships EMPTY: the renderer clears the body anyway, and a
    # template with sample content in it is a template someone will ship by
    # accident.
    body_el = doc.element.body
    from docx.oxml.ns import qn
    for child in list(body_el):
        if child.tag != qn("w:sectPr"):
            body_el.remove(child)

    os.makedirs(os.path.dirname(path), exist_ok=True)
    doc.save(path)
    return path


def main():
    path = build()
    print(f"Wrote {os.path.relpath(path, REPO_ROOT)} "
          f"({os.path.getsize(path):,} bytes)")


if __name__ == "__main__":
    main()
