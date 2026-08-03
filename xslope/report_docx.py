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

"""The Word renderer: a :class:`xslope.report.Report` as a ``.docx``.

The document is built ON a template — ``xslope/resources/report_template.docx``
by default — and inherits everything the template declares: the page size and
margins, the Title / Heading / Body / Caption styles, and the header and footer
frames. Nothing here sets a font or a margin; a company template dropped in its
place restyles the whole report (S4 makes that a dialog setting).

What this module does own is the document's *structure*: the title page, the
table-of-contents field, the header and footer content, the landscape section the
slice table needs, and the mapping from content blocks to Word objects.

**Fields, not frozen text.** The title page, the header and the footer carry real
Word fields — ``DOCPROPERTY`` for the metadata, ``PAGE``/``NUMPAGES`` for the page
count, ``TOC`` for the contents. Each is written with its result already cached,
so the document reads correctly the moment it is opened and in viewers that never
update fields; Word refreshes them on open or on print, which is what makes
"page 3 of 17" and the contents list true after an edit.

The metadata maps onto the Word core properties Word's own field names reach:

    ==================  =================
    Report metadata     Core property
    ==================  =================
    Project title       Title
    Project number      Subject
    Organization        Company / Category
    Author              Author
    ==================  =================
"""

from __future__ import annotations

import os

from docx import Document
from docx.enum.section import WD_ORIENT, WD_SECTION
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.enum.text import WD_ALIGN_PARAGRAPH, WD_BREAK
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Inches, Pt

#: The styles template every report is built on unless the caller names another.
DEFAULT_TEMPLATE = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "resources", "report_template.docx")

#: Style names the renderer asks the template for. A template that does not
#: define one falls back to the document's Normal style rather than failing.
STYLE = {
    "title": "Title",
    "subtitle": "Subtitle",
    "heading": "Heading %d",
    "body": "Body Text",
    "caption": "Caption",
    "bullet": "List Bullet",
    "table": "Table Grid",
}

#: Font size for table body text, and the narrower size a wide table falls back
#: to so the slice table's twenty columns still fit a landscape page.
TABLE_PT = 8.5
WIDE_TABLE_PT = 7.0
WIDE_TABLE_COLUMNS = 12


# ---------------------------------------------------------------------------
# Word field plumbing
# ---------------------------------------------------------------------------

def add_field(paragraph, instruction, cached="", dirty=False):
    """Append a Word field to ``paragraph``, with its result already cached.

    ``dirty`` marks the field for refresh when the document is opened, which is
    what a table of contents needs and a page number does not.
    """
    begin = paragraph.add_run()._r
    fld = OxmlElement("w:fldChar")
    fld.set(qn("w:fldCharType"), "begin")
    if dirty:
        fld.set(qn("w:dirty"), "true")
    begin.append(fld)

    instr_run = paragraph.add_run()._r
    instr = OxmlElement("w:instrText")
    instr.set(qn("xml:space"), "preserve")
    instr.text = instruction
    instr_run.append(instr)

    sep_run = paragraph.add_run()._r
    sep = OxmlElement("w:fldChar")
    sep.set(qn("w:fldCharType"), "separate")
    sep_run.append(sep)

    result = paragraph.add_run(cached)

    end_run = paragraph.add_run()._r
    end = OxmlElement("w:fldChar")
    end.set(qn("w:fldCharType"), "end")
    end_run.append(end)
    return result


def _style(doc, name):
    """A style by name, or None when the template does not define it."""
    try:
        return doc.styles[name]
    except KeyError:
        return None


def _para(doc, text="", style=None, align=None, size=None, bold=None,
          space_after=None):
    p = doc.add_paragraph()
    if style and _style(doc, style) is not None:
        p.style = doc.styles[style]
    if text:
        run = p.add_run(text)
        if size is not None:
            run.font.size = Pt(size)
        if bold is not None:
            run.font.bold = bold
    if align is not None:
        p.alignment = align
    if space_after is not None:
        p.paragraph_format.space_after = Pt(space_after)
    return p


def _repeat_header_row(row):
    """Mark a table row as a header that repeats on every page it spans."""
    tr_pr = row._tr.get_or_add_trPr()
    header = OxmlElement("w:tblHeader")
    header.set(qn("w:val"), "true")
    tr_pr.append(header)


def _cell_text(cell, text, size, bold=False):
    cell.text = ""
    p = cell.paragraphs[0]
    p.paragraph_format.space_before = Pt(1)
    p.paragraph_format.space_after = Pt(1)
    run = p.add_run(str(text))
    run.font.size = Pt(size)
    run.font.bold = bold


# ---------------------------------------------------------------------------
# Page furniture
# ---------------------------------------------------------------------------

def _write_header(section, title):
    """The running head: the project title, as a live document-property field."""
    header = section.header
    header.is_linked_to_previous = False
    for p in list(header.paragraphs[1:]):
        p._element.getparent().remove(p._element)
    p = header.paragraphs[0]
    for r in list(p.runs):
        r._r.getparent().remove(r._r)
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    add_field(p, ' DOCPROPERTY "Title" \\* MERGEFORMAT ', title)
    for run in p.runs:
        run.font.size = Pt(9)


def _write_footer(section, date):
    """The running foot: page N of M, and the report date."""
    footer = section.footer
    footer.is_linked_to_previous = False
    for p in list(footer.paragraphs[1:]):
        p._element.getparent().remove(p._element)
    p = footer.paragraphs[0]
    for r in list(p.runs):
        r._r.getparent().remove(r._r)
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    p.add_run("Page ")
    add_field(p, " PAGE ", "1")
    p.add_run(" of ")
    add_field(p, " NUMPAGES ", "1")
    if date:
        p.add_run(f"     |     {date}")
    for run in p.runs:
        run.font.size = Pt(9)


def _title_page(doc, meta):
    """Title page: organization, title, document type, the metadata block, and —
    only when asked for — prepared-by / checked-by signature lines."""
    for _ in range(3):
        _para(doc, "")

    if meta.get("organization"):
        _para(doc, "", align=WD_ALIGN_PARAGRAPH.CENTER)
        p = doc.paragraphs[-1]
        add_field(p, ' DOCPROPERTY "Company" \\* MERGEFORMAT ', meta["organization"])
        for run in p.runs:
            run.font.size = Pt(13)
            run.font.bold = True
        p.alignment = WD_ALIGN_PARAGRAPH.CENTER

    p = doc.add_paragraph()
    if _style(doc, STYLE["title"]) is not None:
        p.style = doc.styles[STYLE["title"]]
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    add_field(p, ' DOCPROPERTY "Title" \\* MERGEFORMAT ', meta.get("title", ""))

    p = _para(doc, meta.get("document_type", ""), style=STYLE["subtitle"],
              align=WD_ALIGN_PARAGRAPH.CENTER)

    _para(doc, "")
    # "Author" rather than "Prepared by": the signature block below owns that
    # phrase, and a title page must not ask the same question twice.
    rows = [("Project number", meta.get("project_number", ""), "Subject"),
            ("Author", meta.get("author", ""), "Author"),
            ("Date", meta.get("date", ""), None)]
    rows = [r for r in rows if r[1]]
    if rows:
        table = doc.add_table(rows=len(rows), cols=2)
        table.alignment = WD_TABLE_ALIGNMENT.CENTER
        for i, (label, value, prop) in enumerate(rows):
            _cell_text(table.rows[i].cells[0], label + "  ", 10.5, bold=True)
            cell = table.rows[i].cells[1]
            cell.text = ""
            cp = cell.paragraphs[0]
            if prop:
                add_field(cp, f' DOCPROPERTY "{prop}" \\* MERGEFORMAT ', value)
            else:
                cp.add_run(value)
            for run in cp.runs:
                run.font.size = Pt(10.5)

    if meta.get("signature_lines"):
        _para(doc, "")
        _para(doc, "")
        sig = doc.add_table(rows=2, cols=2)
        sig.alignment = WD_TABLE_ALIGNMENT.CENTER
        for col, label in ((0, "Prepared by"), (1, "Checked by")):
            _cell_text(sig.rows[0].cells[col], "_" * 34, 10.5)
            _cell_text(sig.rows[1].cells[col], f"{label}                    Date", 9)

    doc.add_paragraph().add_run().add_break(WD_BREAK.PAGE)


def _contents_page(doc):
    """The table of contents: a real TOC field, marked for refresh on open."""
    _para(doc, "Table of Contents", size=14, bold=True, space_after=10)
    p = doc.add_paragraph()
    add_field(p, ' TOC \\o "1-3" \\h \\z \\u ',
              "The table of contents is built when this document is opened in "
              "Word. To build it now, right-click here and choose Update Field.",
              dirty=True)
    for run in p.runs:
        run.font.size = Pt(9.5)
        run.font.italic = True
    doc.add_paragraph().add_run().add_break(WD_BREAK.PAGE)


def _set_orientation(section, landscape, margin_in):
    """Rotate a section, swapping the page dimensions with it."""
    width, height = section.page_width, section.page_height
    is_landscape = width > height
    if landscape != is_landscape:
        section.page_width, section.page_height = height, width
    section.orientation = WD_ORIENT.LANDSCAPE if landscape else WD_ORIENT.PORTRAIT
    for side in ("left_margin", "right_margin"):
        setattr(section, side, Inches(margin_in))


def ensure_orientation(doc, state, landscape):
    """Put the document into the requested orientation, starting a new Word
    section when it is not already there.

    A new section inherits the previous one's properties, including its
    "different first page" flag — which would blank the header and footer on the
    first page of every section the report opens. It is cleared here, and the
    header and footer are relinked, so the running head and the page count carry
    across a landscape page rather than stopping at it.
    """
    if landscape == state["landscape"]:
        return state["section"]
    section = doc.add_section(WD_SECTION.NEW_PAGE)
    section.different_first_page_header_footer = False
    section.header.is_linked_to_previous = True
    section.footer.is_linked_to_previous = True
    _set_orientation(section, landscape, 0.75 if landscape else 1.0)
    state["section"], state["landscape"] = section, landscape
    return section


# ---------------------------------------------------------------------------
# Blocks
# ---------------------------------------------------------------------------

def _content_width_in(section):
    return (section.page_width - section.left_margin - section.right_margin) / 914400


def _render_table(doc, block, section):
    """One table, its caption above it and its legend beneath."""
    n_cols = max(1, len(block.headers))
    size = WIDE_TABLE_PT if n_cols > WIDE_TABLE_COLUMNS else TABLE_PT

    if block.caption:
        _para(doc, f"Table {block.number}. {block.caption}", style=STYLE["caption"])

    table = doc.add_table(rows=1, cols=n_cols)
    if _style(doc, STYLE["table"]) is not None:
        table.style = doc.styles[STYLE["table"]]
    table.autofit = True
    for j, head in enumerate(block.headers):
        _cell_text(table.rows[0].cells[j], head, size, bold=True)
    _repeat_header_row(table.rows[0])
    for row in block.rows:
        cells = table.add_row().cells
        for j in range(n_cols):
            _cell_text(cells[j], row[j] if j < len(row) else "", size)

    if block.legend:
        p = doc.add_paragraph()
        p.paragraph_format.space_before = Pt(4)
        for i, (term, definition) in enumerate(block.legend):
            if i:
                p.add_run("   ")
            run = p.add_run(f"{term}: ")
            run.font.bold = True
            run.font.size = Pt(7.5)
            run = p.add_run(definition)
            run.font.size = Pt(7.5)
    _para(doc, "")


def _render_keyvalues(doc, block):
    if block.title:
        _para(doc, block.title, bold=True, space_after=2)
    table = doc.add_table(rows=len(block.items), cols=2)
    table.autofit = True
    for i, (label, value) in enumerate(block.items):
        _cell_text(table.rows[i].cells[0], label, 10)
        _cell_text(table.rows[i].cells[1], value, 10)
        table.rows[i].cells[0].paragraphs[0].runs[0].font.bold = True
    _para(doc, "")


def _render_figure(doc, block, section):
    width = min(block.width_in, _content_width_in(section) - 0.05)
    if os.path.exists(block.path):
        p = doc.add_paragraph()
        p.alignment = WD_ALIGN_PARAGRAPH.CENTER
        p.add_run().add_picture(block.path, width=Inches(width))
    _para(doc, f"Figure {block.number}. {block.caption}", style=STYLE["caption"],
          align=WD_ALIGN_PARAGRAPH.CENTER)


def _render_blocks(doc, blocks, state):
    """Render a section's blocks, opening and closing a landscape page around any
    block that asks for one."""
    for block in blocks:
        ensure_orientation(doc, state, bool(getattr(block, "landscape", False)))
        if block.kind == "prose":
            _para(doc, block.text, style=STYLE["body"])
        elif block.kind == "bullets":
            for item in block.items:
                _para(doc, item, style=STYLE["bullet"])
        elif block.kind == "keyvalues":
            _render_keyvalues(doc, block)
        elif block.kind == "figure":
            _render_figure(doc, block, state["section"])
        elif block.kind == "table":
            _render_table(doc, block, state["section"])


def _render_section(doc, section_node, level, state):
    # A heading always opens on a portrait page: a landscape table is a place the
    # report visits for one table, never a place the next section starts in.
    ensure_orientation(doc, state, False)
    style = STYLE["heading"] % min(level, 3)
    _para(doc, section_node.title, style=style)
    _render_blocks(doc, section_node.blocks, state)
    for child in section_node.children:
        _render_section(doc, child, level + 1, state)


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def render_docx(report, path, template=None):
    """Write ``report`` (a :class:`xslope.report.Report`) to ``path`` as a
    ``.docx``, built on ``template`` (the shipped default when None).

    Returns the path written.
    """
    template = template or DEFAULT_TEMPLATE
    doc = Document(template) if os.path.exists(template) else Document()

    meta = report.meta or {}
    props = doc.core_properties
    props.title = meta.get("title", "")
    props.subject = meta.get("project_number", "")
    props.author = meta.get("author", "") or props.author
    props.category = meta.get("organization", "")
    props.comments = "Generated by xslope."

    # A template opened as a starting document may carry sample body content; the
    # report replaces it entirely, and only the styles and page setup are kept.
    body = doc.element.body
    for child in list(body):
        if child.tag != qn("w:sectPr"):
            body.remove(child)

    section = doc.sections[0]
    section.different_first_page_header_footer = True
    _write_header(section, meta.get("title", ""))
    _write_footer(section, meta.get("date", ""))

    _title_page(doc, meta)
    _contents_page(doc)

    state = {"section": section, "landscape": False}
    for node in report.sections:
        _render_section(doc, node, 1, state)

    # Never end on a landscape page: the report closes in the orientation it
    # opened in, so a document appended to it starts straight.
    ensure_orientation(doc, state, False)

    doc.save(path)
    return path
