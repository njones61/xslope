"""Word-level diff of a docs page between two git revisions, as HTML (and PDF).

Usage: python tools/tutorial_diff.py <base-rev> <head-rev> <page.md> [<page.md> ...]
       --out DIR   (default: xslope_private/reports/tutorial_diffs)
Deletions are struck through in red, insertions in green; unchanged text is
gray. Test tags and image lines are kept so a moved number or figure shows.
The PDF is produced by LibreOffice headless (never Word).
"""
import difflib, html, os, re, subprocess, sys

SOFFICE = "/Applications/LibreOffice.app/Contents/MacOS/soffice"
# Inline styles: LibreOffice's HTML import ignores <style> blocks but honors these.
DEL = '<del style="color:#b00020;background:#ffecec;text-decoration:line-through">'
INS = '<ins style="color:#0a6d1f;background:#e8f7ea;text-decoration:none;font-weight:bold">'
HTML_BLOCK = re.compile(r"^\s*<(div|/div|p|span|table|tr|td)", re.M)
STYLE = """<style>body{font:13px/1.5 -apple-system,Helvetica,Arial;margin:32px;color:#444}
h1{font-size:18px}h2{font-size:15px;border-top:1px solid #ddd;padding-top:12px;margin-top:24px}
del{color:#b00020;background:#ffecec;text-decoration:line-through}
ins{color:#0a6d1f;background:#e8f7ea;text-decoration:none}
p{margin:0 0 10px 0;white-space:pre-wrap}.tag{color:#888;font:11px monospace}</style>"""

def git_show(rev, path):
    r = subprocess.run(["git", "show", f"{rev}:{path}"], capture_output=True, text=True)
    return r.stdout if r.returncode == 0 else ""

def words(text):
    return re.findall(r"\s+|[^\s]+", text)

def diff_html(old, new):
    a, b = words(old), words(new)
    out = []
    for op, i1, i2, j1, j2 in difflib.SequenceMatcher(None, a, b, autojunk=False).get_opcodes():
        if op == "equal":
            out.append(html.escape("".join(a[i1:i2])))
        elif op == "delete":
            out.append(DEL + html.escape("".join(a[i1:i2])) + "</del>")
        elif op == "insert":
            out.append(INS + html.escape("".join(b[j1:j2])) + "</ins>")
        else:
            out.append(DEL + html.escape("".join(a[i1:i2])) + "</del>"
                       + INS + html.escape("".join(b[j1:j2])) + "</ins>")
    return "".join(out)

def paragraphs(text):
    # Pure HTML furniture (the glance box, pills) is noise in a prose diff.
    return [p for p in re.split(r"\n\s*\n", text) if not HTML_BLOCK.match(p.strip()[:6])]

def page_diff(base, head, path):
    old, new = git_show(base, path), git_show(head, path)
    po, pn = paragraphs(old), paragraphs(new)
    sm = difflib.SequenceMatcher(None, po, pn, autojunk=False)
    parts = [f"<h1>{html.escape(path)} — {html.escape(base)} → {html.escape(head)}</h1>"]
    changed = 0
    for op, i1, i2, j1, j2 in sm.get_opcodes():
        if op == "equal":
            continue
        changed += 1
        old_blk, new_blk = "\n\n".join(po[i1:i2]), "\n\n".join(pn[j1:j2])
        parts.append("<p>" + diff_html(old_blk, new_blk) + "</p>")
    if changed == 0:
        parts.append("<p>(no changes)</p>")
    return "<!doctype html><html><head><meta charset='utf-8'>" + STYLE + "</head><body>" + "\n".join(parts) + "</body></html>", changed

def main(argv):
    out_dir = "/Users/njones/python_projects/xslope_private/reports/tutorial_diffs"
    if "--out" in argv:
        i = argv.index("--out"); out_dir = argv[i + 1]; del argv[i:i + 2]
    base, head, pages = argv[0], argv[1], argv[2:]
    os.makedirs(out_dir, exist_ok=True)
    for page in pages:
        doc, n = page_diff(base, head, page)
        stem = os.path.splitext(os.path.basename(page))[0] + f"_{base}_to_{head}".replace("/", "-")
        hp = os.path.join(out_dir, stem + ".html")
        open(hp, "w").write(doc)
        subprocess.run([SOFFICE, "--headless", "--convert-to", "pdf", "--outdir", out_dir, hp],
                       capture_output=True, text=True)
        print(f"{page}: {n} changed block(s) -> {os.path.join(out_dir, stem + '.pdf')}")

if __name__ == "__main__":
    main(sys.argv[1:])
