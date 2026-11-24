"""What the template grays, what Studio grays, and what the engine reads, field by
field -- the three readings ``table_paste_check`` holds to one another.

A grayed cell is a promise that nothing reads it, and a live one a promise that
something does. The template keeps its promise in conditional formatting, Studio in
its editors, and the engine in what it does with the value; when the three are kept
by hand they drift, and a check that compares only two of them passes while both are
wrong. So each reading here is taken from the thing itself:

  THE SHEET. The conditional formatting of the template's mat, circles, polygon and
     reinforce sheets is read from the workbook and evaluated by a small Excel
     formula evaluator, cell by cell, for every case -- a blank switch is a case of
     its own, a row's name cell is filled where a rule guards on it, and the
     polygon sheet's block-absolute references ($B$5) are resolved per block.
  THE ENGINE. Nothing is looked up in a table kept here. Every field of every case
     is nudged and the engine's own inputs compared before and after: the slices
     ``generate_slices`` cuts on one circle of a toy model (and again with
     ``rapid=True``), the arrays ``build_fem_data`` makes on a small mesh (with the
     inputs the solver itself holds out masked as it masks them), the seepage
     arrays put through the solvers' own conductivity and storage laws,
     reliability's parameter set (rapid and not), and for circles and polygons what
     the loader hands on from a workbook copy. A nudge is what the loader makes of
     the cell nudged, learned from a calibration load (Spacing, for one, rescales
     the per-element terms); where the loader reads a blank cell as something other
     than a typed zero (None, NaN, a default), clearing the cell is a nudge of its
     own, and a blank cell is nudged by filling it. A field no engine reads but a
     check does is REFUSES when a nudge moves an ERROR, NOTES when it moves only a
     WARNING or INFO.
  STUDIO. Every view that grays: the materials table and list view (disabled or
     hidden), the circles table, the polygons editor's Material and Size widgets
     per Type, and the reinforcement table and list view, each read off the real
     dialog.

Probing stays cheap because the engine is deterministic. A probe's cases fold to
the switch values it treats alike (a value whose outputs match another's at every
setting of the other switches joins its class), fields are nudged together and
split only where an output moves, and a field one probe has found read is not
nudged again by the next. A probe refused in every case is a probe that cannot
call the engine, and says so rather than reporting every field unread.
"""
import contextlib
import copy
import io
import itertools
import math
import os
import re
import shutil
import sys
import tempfile
import warnings

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(_HERE)
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)

#: The workbook whose formatting is the sheet's reading and whose copy is the toy
#: model. A module attribute so a run can point it at another template.
TEMPLATE = os.path.join(_REPO, "docs", "inputs", "input_template.xlsx")

_CACHE = {}


@contextlib.contextmanager
def quiet():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        with contextlib.redirect_stdout(io.StringIO()), \
                contextlib.redirect_stderr(io.StringIO()):
            yield


def _scratch():
    if "dir" not in _CACHE:
        _CACHE["dir"] = tempfile.mkdtemp(prefix="graying_probe_")
    return _CACHE["dir"]


def cleanup():
    d = _CACHE.pop("dir", None)
    if d:
        shutil.rmtree(d, ignore_errors=True)
    _CACHE.clear()


# =========================================================================== #
# 1. The sheet: an evaluator for the template's conditional formatting
# =========================================================================== #
_TOKEN = re.compile(r"""\s*(?:
    (?P<str>"(?:[^"]|"")*") |
    (?P<num>\d+(?:\.\d*)?(?:[eE][-+]?\d+)?) |
    (?P<ref>\$?[A-Z]{1,3}\$?\d+) |
    (?P<name>[A-Z][A-Z0-9_.]*) |
    (?P<op><>|<=|>=|[=<>(),&])
)""", re.VERBOSE)


def _col_num(letters):
    n = 0
    for ch in letters:
        n = n * 26 + ord(ch) - 64
    return n


def _col_letters(n):
    s = ""
    while n:
        n, r = divmod(n - 1, 26)
        s = chr(65 + r) + s
    return s


def _tokenize(text):
    pos, out = 0, []
    text = text.strip()
    while pos < len(text):
        m = _TOKEN.match(text, pos)
        if not m or m.end() == pos:
            raise ValueError(f"cannot read the formula {text!r} at {text[pos:]!r}")
        pos = m.end()
        kind = m.lastgroup
        out.append((kind, m.group(kind)))
    return out


class Formula:
    """One conditional-formatting formula, parsed once, evaluated per cell.

    The subset the template uses and a little more: string and number literals,
    cell references (relative, column-absolute or fully absolute), comparisons, &,
    and AND / OR / NOT / ISNUMBER / ISBLANK / ISTEXT. Anything else is refused by
    name, so a new construct in the template fails this check instead of being
    read wrongly. Values follow Excel: text compares without regard to case, an
    empty cell equals "" and 0, and text never equals a number."""

    FUNCS = ("AND", "OR", "NOT", "ISNUMBER", "ISBLANK", "ISTEXT")

    def __init__(self, text, anchor):
        self.text = text
        self.anchor = anchor                 # (col, row) the formula is written at
        self._toks = _tokenize(text)
        self._i = 0
        self.ast = self._expr()
        if self._i != len(self._toks):
            raise ValueError(f"cannot read the formula {text!r}")
        self.refs = []
        self._collect(self.ast)

    # --- parsing ---
    def _peek(self):
        return self._toks[self._i] if self._i < len(self._toks) else (None, None)

    def _take(self, want=None):
        tok = self._peek()
        if want is not None and tok[1] != want:
            raise ValueError(f"cannot read the formula {self.text!r}: expected {want}")
        self._i += 1
        return tok

    def _expr(self):
        left = self._concat()
        kind, val = self._peek()
        if kind == "op" and val in ("=", "<>", "<", ">", "<=", ">="):
            self._take()
            return ("cmp", val, left, self._concat())
        return left

    def _concat(self):
        node = self._atom()
        while self._peek() == ("op", "&"):
            self._take()
            node = ("cat", node, self._atom())
        return node

    def _atom(self):
        kind, val = self._take()
        if kind == "str":
            return ("lit", val[1:-1].replace('""', '"'))
        if kind == "num":
            return ("lit", float(val))
        if kind == "ref":
            m = re.fullmatch(r"(\$?)([A-Z]+)(\$?)(\d+)", val)
            return ("ref", m.group(1) == "$", _col_num(m.group(2)),
                    m.group(3) == "$", int(m.group(4)))
        if kind == "name" and val in ("TRUE", "FALSE"):
            return ("lit", val == "TRUE")
        if kind == "name":
            if val not in self.FUNCS:
                raise ValueError(f"the formula {self.text!r} uses {val}, which this "
                                 f"check does not evaluate")
            self._take("(")
            args = []
            if self._peek() != ("op", ")"):
                args.append(self._expr())
                while self._peek() == ("op", ","):
                    self._take()
                    args.append(self._expr())
            self._take(")")
            return ("fn", val, args)
        if (kind, val) == ("op", "("):
            node = self._expr()
            self._take(")")
            return node
        raise ValueError(f"cannot read the formula {self.text!r}")

    def _collect(self, node):
        if node[0] == "ref":
            self.refs.append(node)
        elif node[0] == "cmp":
            self._collect(node[2])
            self._collect(node[3])
        elif node[0] == "cat":
            self._collect(node[1])
            self._collect(node[2])
        elif node[0] == "fn":
            for a in node[2]:
                self._collect(a)

    # --- evaluation ---
    def ref_cell(self, ref, at):
        """The cell a reference names when the formula is applied at ``at``."""
        _, cabs, col, rabs, row = ref
        c0, r0 = self.anchor
        return (col if cabs else at[0] + col - c0, row if rabs else at[1] + row - r0)

    def evaluate(self, cells, at):
        return bool(self._eval(self.ast, cells, at))

    def _eval(self, node, cells, at):
        tag = node[0]
        if tag == "lit":
            return node[1]
        if tag == "ref":
            return cells.get(self.ref_cell(node, at))
        if tag == "cat":
            return (_text(self._eval(node[1], cells, at))
                    + _text(self._eval(node[2], cells, at)))
        if tag == "cmp":
            a = self._eval(node[2], cells, at)
            b = self._eval(node[3], cells, at)
            op = node[1]
            if op == "=":
                return _xl_eq(a, b)
            if op == "<>":
                return not _xl_eq(a, b)
            return _xl_order(a, b, op)
        name, args = node[1], node[2]
        vals = [self._eval(a, cells, at) for a in args]
        if name == "AND":
            return all(bool(v) for v in vals)
        if name == "OR":
            return any(bool(v) for v in vals)
        if name == "NOT":
            return not bool(vals[0])
        if name == "ISNUMBER":
            return isinstance(vals[0], (int, float)) and not isinstance(vals[0], bool)
        if name == "ISBLANK":
            return vals[0] is None
        if name == "ISTEXT":
            return isinstance(vals[0], str)
        raise ValueError(name)


def _text(v):
    if v is None:
        return ""
    if isinstance(v, float) and v.is_integer():
        return str(int(v))
    return str(v)


def _xl_eq(a, b):
    if a is None and b is None:
        return True
    if a is None:
        a, b = b, a
    if b is None:                       # an empty cell is "" and 0 at once
        return a == "" or (isinstance(a, (int, float)) and a == 0)
    if isinstance(a, str) and isinstance(b, str):
        return a.lower() == b.lower()
    if isinstance(a, str) or isinstance(b, str):
        return False
    return a == b


def _xl_order(a, b, op):
    a = 0.0 if a is None else a
    b = 0.0 if b is None else b
    if isinstance(a, str) != isinstance(b, str):
        # Excel orders every number below every text.
        res = {"<": isinstance(b, str), "<=": isinstance(b, str),
               ">": isinstance(a, str), ">=": isinstance(a, str)}
        return res[op]
    if isinstance(a, str):
        a, b = a.lower(), b.lower()
    return {"<": a < b, "<=": a <= b, ">": a > b, ">=": a >= b}[op]


def _parse_ranges(sqref):
    from openpyxl.utils import range_boundaries
    out = []
    for part in str(sqref).split():
        c0, r0, c1, r1 = range_boundaries(part)
        out.append((c0, r0, c1, r1))
    return out


def _in(ranges, cell):
    return any(c0 <= cell[0] <= c1 and r0 <= cell[1] <= r1
               for c0, r0, c1, r1 in ranges)


class Sheet:
    """A template sheet's graying: its conditional-formatting rules (every rule is
    a gray one), its dropdown (list) validations, and its cell values."""

    def __init__(self, ws):
        self.name = ws.title
        self.rules = []
        for rng, cfs in ws.conditional_formatting._cf_rules.items():
            ranges = _parse_ranges(rng.sqref)
            for cf in cfs:
                if cf.type != "expression":
                    raise ValueError(f"{ws.title}!{rng.sqref}: a {cf.type} rule, which "
                                     f"this check does not read")
                anchor = (ranges[0][0], ranges[0][1])
                self.rules.append((ranges, Formula(cf.formula[0], anchor)))
        self.dropdowns = [(_parse_ranges(dv.sqref), dv.formula1)
                          for dv in ws.data_validations.dataValidation
                          if dv.type == "list"]

    def gray(self, cells, cell):
        """Does any rule covering ``cell`` gray it, given the cell values?"""
        return any(_in(ranges, cell) and f.evaluate(cells, cell)
                   for ranges, f in self.rules)



def template_sheets(path=None):
    import openpyxl
    path = path or TEMPLATE
    key = ("sheets", path)
    if key not in _CACHE:
        wb = openpyxl.load_workbook(path)
        _CACHE[key] = {n: Sheet(wb[n]) for n in
                       ("mat", "circles", "polygon", "reinforce", "piles")}
        _CACHE[key]["_header"] = {
            n: {c.column: c.value for c in wb[n][row] if c.value is not None}
            for n, row in (("mat", 10), ("circles", 2), ("reinforce", 2))}
    return _CACHE[key]


def structural_failures(path=None, allow=()):
    """Each graying rule must cover every row its switch's dropdown offers: a rule
    whose switch has a dropdown (a list validation on the referenced column) must
    reach the dropdown's last row; one guarding on cells with no dropdown must
    reach the sheet's dropdown rows. On the polygon sheet, where the switch is one
    cell per block, every block's Type cell must be referenced by a rule.
    ``allow`` holds (sheet, column letter, first uncovered row) exceptions."""
    sheets = template_sheets(path)
    out = []
    for name in ("mat", "circles", "reinforce", "piles"):
        sh = sheets[name]
        all_rows = set()
        for ranges, _ in sh.dropdowns:
            for c0, r0, c1, r1 in ranges:
                all_rows |= set(range(r0, r1 + 1))
        for ranges, f in sh.rules:
            rule_rows = {r for c0, r0, c1, r1 in ranges for r in range(r0, r1 + 1)}
            # A row switch: a reference with a relative row, once per column.
            for col, shift in sorted({(ref[2], ref[4] - f.anchor[1])
                                      for ref in f.refs if not ref[3]}):
                need = set()
                for dranges, _ in sh.dropdowns:
                    for c0, r0, c1, r1 in dranges:
                        if c0 <= col <= c1:
                            need |= set(range(r0 - shift, r1 + 1 - shift))
                where = (f"its switch {_col_letters(col)} offers its dropdown"
                         if need else "the sheet's dropdowns reach")
                if not need:                # a cell with no dropdown of its own
                    need = all_rows
                missing = sorted(need - rule_rows)
                if missing and (name, _col_letters(col), missing[0]) not in allow:
                    out.append(f"{name}: the rule {f.text!r} grays rows "
                               f"{min(rule_rows)}-{max(rule_rows)}, but {where} down "
                               f"to row {max(need)}: rows {missing[0]}-{missing[-1]} "
                               f"are never grayed")
    sh = sheets["polygon"]
    referenced = {f.ref_cell(ref, f.anchor) for _, f in sh.rules for ref in f.refs}
    for ranges, _ in sh.dropdowns:
        for c0, r0, c1, r1 in ranges:
            for c in range(c0, c1 + 1):
                for r in range(r0, r1 + 1):
                    if (c, r) not in referenced:
                        out.append(f"polygon: the Type cell {_col_letters(c)}{r} has "
                                   f"a dropdown but no graying rule reads it")
    return out


# =========================================================================== #
# 2. Comparing what an engine was handed
# =========================================================================== #
def same(a, b):
    """Deep equality over what the probes return: numpy arrays, DataFrames,
    shapely geometry, dicts, lists, sets, floats (NaN equal to NaN)."""
    import numpy as np
    import pandas as pd
    if a is b:
        return True
    if isinstance(a, pd.DataFrame) or isinstance(b, pd.DataFrame):
        return (isinstance(a, pd.DataFrame) and isinstance(b, pd.DataFrame)
                and list(a.columns) == list(b.columns) and a.equals(b))
    if isinstance(a, np.ndarray) or isinstance(b, np.ndarray):
        try:
            a2, b2 = np.asarray(a), np.asarray(b)
        except Exception:
            return False
        if a2.shape != b2.shape:
            return False
        if a2.dtype.kind in "fc" or b2.dtype.kind in "fc":
            try:
                return bool(np.array_equal(a2, b2, equal_nan=True))
            except TypeError:
                pass
        if a2.dtype == object or b2.dtype == object:
            return all(same(x, y) for x, y in zip(a2.ravel().tolist(),
                                                   b2.ravel().tolist()))
        return bool(np.array_equal(a2, b2))
    if hasattr(a, "geom_type") or hasattr(b, "geom_type"):
        return (hasattr(a, "geom_type") and hasattr(b, "geom_type")
                and a.equals_exact(b, 0.0))
    if isinstance(a, dict) and isinstance(b, dict):
        return a.keys() == b.keys() and all(same(a[k], b[k]) for k in a)
    if isinstance(a, (list, tuple)) and isinstance(b, (list, tuple)):
        return len(a) == len(b) and all(same(x, y) for x, y in zip(a, b))
    if isinstance(a, (set, frozenset)) and isinstance(b, (set, frozenset)):
        return sorted(map(repr, a)) == sorted(map(repr, b))
    if isinstance(a, float) and isinstance(b, float) and math.isnan(a) and math.isnan(b):
        return True
    if (type(a) is type(b) and type(a).__eq__ is object.__eq__
            and hasattr(a, "__dict__")):
        return same(vars(a), vars(b))           # a plain record object
    try:
        return bool(a == b)
    except Exception:
        return repr(a) == repr(b)


def _attempt(fn, *args, **kw):
    """A probe's output, or the refusal it raised (type and message)."""
    try:
        with quiet():
            return ("ok", fn(*args, **kw))
    except Exception as exc:                        # a refusal is an output too
        return ("refused", type(exc).__name__, str(exc))


def _refused(out):
    if isinstance(out, list):
        return all(_refused(o) for o in out)
    return isinstance(out, tuple) and len(out) == 3 and out[0] == "refused"


def _must_run(name, outputs):
    """A probe refused in every case is not a reading, it is a probe that cannot
    call the engine (a changed signature, a model the engine no longer takes):
    say so rather than report every field unread."""
    if outputs and all(_refused(o) for o in outputs.values()):
        first = next(iter(outputs.values()))
        raise AssertionError(f"the {name} probe was refused in every case: {first}")


# Verdicts, from strongest to weakest.
READ, REFUSES, NOTES, UNREAD = "read", "refuses", "notes", "unread"


def expect_gray(verdict):
    """The graying a verdict asks for: a field something reads, or whose value a
    check refuses a run over, must stay live (a grayed Studio cell cannot be edited
    to clear it); one nothing reads, or that only a WARNING or INFO mentions,
    grays."""
    return verdict in (UNREAD, NOTES)


# =========================================================================== #
# 3. Preflight, for the fields no engine reads
# =========================================================================== #
#: Every analysis preflight knows, each with the selection that makes its
#: composite rules apply. The union of what any of them says is the checks' reading.
_PREFLIGHT_RUNS = (("lem", {}), ("rapid", {}), ("seismic2", {}), ("seep", {}),
                   ("tseep", {}), ("fem", {}), ("ssrm", {}),
                   ("reliability", {"base": "lem"}), ("reliability", {"base": "fem"}))


#: Left out of the checks' reading: a rule that only asks the slicer whether the
#: toy circle can be cut, with pore pressure neutralized. It reads no field the
#: slicer's own probe does not, and it re-slices the circle in every analysis.
_SKIP_RULES = ("surface.circle_cannot_be_sliced",)


def findings(sd, mesh=None, runs=None):
    """analysis -> the sorted (rule, severity, message) findings preflight makes of
    ``sd`` for it: the checks' whole reading of a model, run by run."""
    from xslope.preflight import preflight
    out = {}
    for i, (an, sel) in enumerate(_PREFLIGHT_RUNS):
        if runs is not None and i not in runs:
            continue
        sel = dict(sel)
        if mesh is not None:
            sel["mesh"] = mesh
        try:
            with quiet():
                rep = preflight(sd, an, sel, skip=_SKIP_RULES)
            out[i] = sorted((f.rule_id, f.severity, f.message) for f in rep.findings)
        except Exception as exc:
            out[i] = [("raised", "error", repr(exc))]
    return out


def _changed(base, other):
    """run -> the findings that differ between two readings, for the runs that
    differ."""
    from collections import Counter
    out = {}
    for i in other:
        a, b = Counter(base[i]), Counter(other[i])
        d = list(((a - b) + (b - a)).elements())
        if d:
            out[i] = d
    return out


def split_test(keys, moved, on_unread, on_moved):
    """Adaptive group testing. ``moved(keys, ctx)`` nudges ``keys`` together and
    returns a falsy value when nothing moved, else a context passed down to the
    halves (so a split can re-ask only what moved). A group that moves nothing is
    unread whole; one that moves is halved until single keys are left."""
    def rec(ks, ctx):
        got = moved(ks, ctx)
        if not got:
            for k in ks:
                on_unread(k)
        elif len(ks) == 1:
            on_moved(ks[0], got)
        else:
            half = len(ks) // 2
            rec(ks[:half], got)
            rec(ks[half:], got)
    if keys:
        rec(list(keys), None)


def value_classes(cases, output, eq=None):
    """How a probe groups the values of each case dimension: per dimension, value
    -> the first value that moves no output, at any setting of the other
    dimensions, differently from it. A dimension the probe ignores folds to one
    class; the words a probe treats alike (three support Types that fill the same
    Dir and Appl) fold together. ``output(case)`` is the probe's output."""
    eq = eq or same
    out = []
    for i in range(len(cases[0])):
        by_val = {}
        for c in cases:
            by_val.setdefault(c[i], {})[c[:i] + c[i + 1:]] = output(c)
        reps = []
        rep_of = {}
        for v, outs in by_val.items():
            for r in reps:
                other = by_val[r]
                if other.keys() == outs.keys() and all(eq(outs[o], other[o])
                                                       for o in outs):
                    rep_of[v] = r
                    break
            else:
                reps.append(v)
                rep_of[v] = v
        out.append(rep_of)
    return out


def classes_of(cases, reps, canon=None):
    """The cases grouped by their dimensions' value classes (see value_classes),
    each case read as ``canon`` maps it (the case the loader reads it as)."""
    out = {}
    for c in cases:
        cc = canon[c] if canon else c
        out.setdefault(tuple(r[v] for r, v in zip(reps, cc)), []).append(c)
    return list(out.values())


def factored_checks(cases, canon, keys_for, model_for, change, mesh=None):
    """What the checks make of the fields no engine reads: (case, field) ->
    (REFUSES if a nudge moves an ERROR, NOTES if it moves only a WARNING or INFO,
    else UNREAD; the rules it moved).

    Keys are (field, kind): ``keys_for(case)`` names the ones to ask about,
    ``change(case, key)`` is what that nudge does to the record, and
    ``model_for(case, changes)`` builds the model. ``canon`` maps a case to the one
    the loader reads it as. Each preflight analysis is a probe of its own, reduced
    to the case dimensions its findings depend on; the keys of one kind are nudged
    together per class and split where that analysis's findings move."""
    uniq = sorted(set(canon.values()), key=cases.index)
    base = {c: findings(model_for(c, {}), mesh) for c in uniq}
    got = {}
    for run in range(len(_PREFLIGHT_RUNS)):
        dims = value_classes(uniq, lambda c: base[c][run], eq=lambda a, b: a == b)
        for members in classes_of(cases, dims, canon):
            c0 = canon[members[0]]
            todo = sorted({key for c in members for key in keys_for(c)})

            def moved(keys, _ctx, c0=c0, run=run):
                changes = {}
                for key in keys:
                    changes.update(change(c0, key))
                return _changed({run: base[c0][run]},
                                findings(model_for(c0, changes), mesh, {run}))

            def on_moved(key, diff, members=members):
                rows = [f for d in diff.values() for f in d]
                v = REFUSES if any(s == "error" for _, s, _ in rows) else NOTES
                for c in members:
                    if key not in keys_for(c):
                        continue
                    old = got.get((c, key[0]), (UNREAD, []))
                    worse = REFUSES if REFUSES in (old[0], v) else NOTES
                    got[(c, key[0])] = (worse,
                                        sorted(set(old[1]) | {r for r, _, _ in rows}))
            for kind in sorted({key[1] for key in todo}):
                split_test([key for key in todo if key[1] == kind], moved,
                           lambda k: None, on_moved)
    return {(c, key[0]): got.get((c, key[0]), (UNREAD, []))
            for c in cases for key in keys_for(c)}


def _blank(v):
    return v is None or (isinstance(v, float) and math.isnan(v)) or (
        isinstance(v, str) and not v.strip())


def loader_changes(records, fields, cleared_as_zero=True):
    """What the loader makes of each field's cell nudged and cleared, from a
    calibration load: ``records`` is the loaded base record, then one per field
    with its cell nudged, then one per field with it cleared. Returns
    field -> {"nudge": changes, "clear": changes or None} as changes to the record
    (a cell the loader combines with others, such as Spacing, changes several).
    A cleared cell the loader reads as a typed zero is no nudge of its own."""
    base = records[0]
    n = len(fields)
    out = {}
    for i, k in enumerate(fields):
        nud = {kk: v for kk, v in records[1 + i].items() if not same(v, base.get(kk))}
        clr = {kk: v for kk, v in records[1 + n + i].items()
               if not same(v, base.get(kk))}
        if cleared_as_zero and set(clr) == {k} and not _blank(clr[k]) \
                and clr[k] == 0:
            clr = None
        out[k] = {"nudge": nud, "clear": clr or None, "base": base}
    return out


def nudge_keys(changes):
    return ([(k, "nudge") for k in changes]
            + [(k, "clear") for k, ch in changes.items() if ch["clear"]])


def apply_change(record, changes, key):
    """The record's changes for one nudge key: a filled field nudged or cleared as
    the calibration did it; a blank one filled with the calibration's value (and
    not cleared again)."""
    field, kind = key
    ch = changes[field]
    if _blank(record.get(field)):
        if kind == "clear":
            return {}
        return {kk: ch["base"][kk] for kk in ch["nudge"]}
    return dict(ch[kind] or {})


# =========================================================================== #
# 4. The toy model
# =========================================================================== #
# A 10-high slope at 1:2 over level ground; the piezometric line runs up through
# the toe and stays below the crest, so the one circle cuts slice bases below it
# (pore pressure, gsat) and above it (suction); the second line sits lower (the
# drawdown); r-elev 15 falls inside the slide; a head on the left side and an exit
# face on the slope give the seepage arrays their boundary.
_GROUND = [(0.0, 10.0), (20.0, 10.0), (40.0, 20.0), (60.0, 20.0)]
_PIEZO = [(0.0, 11.0), (30.0, 13.5), (60.0, 15.0)]
_PIEZO2 = [(0.0, 10.2), (30.0, 11.0), (60.0, 13.0)]
_CIRCLE = (35.0, 32.0, 26.0)             # Xo, Yo, R


def _y_on(line, x):
    import numpy as np
    xs, ys = zip(*line)
    return float(np.interp(x, xs, ys))


def _put(ws, row, col, value):
    """Set a cell, blank included: ``ws.cell(r, c, None)`` leaves a cell as it was."""
    ws.cell(row, col).value = value


def toy_path():
    """A copy of the template with the toy geometry, piezometric lines and circle
    written in, and one placeholder material."""
    import openpyxl
    key = ("toy_path", TEMPLATE)
    if key in _CACHE:
        return _CACHE[key]
    path = os.path.join(_scratch(), "toy.xlsx")
    shutil.copy(TEMPLATE, path)
    wb = openpyxl.load_workbook(path)
    ws = wb["profile"]
    for i, (x, y) in enumerate(_GROUND):
        _put(ws, 9 + i, 1, x)
        _put(ws, 9 + i, 2, y)
    ws = wb["piezo"]
    for i, ((x1, y1), (x2, y2)) in enumerate(zip(_PIEZO, _PIEZO2)):
        _put(ws, 5 + i, 1, x1)
        _put(ws, 5 + i, 2, y1)
        _put(ws, 5 + i, 4, x2)
        _put(ws, 5 + i, 5, y2)
    ws = wb["circles"]
    for col, v in zip("BCDH", (_CIRCLE[0], _CIRCLE[1], "Radius", _CIRCLE[2])):
        ws[f"{col}3"] = v
    ws = wb["mat"]
    for col, v in zip("BCEFGMNO", ("soil", 120, "mc", 200, 30, 1e5, 0.3, "none")):
        ws[f"{col}11"] = v
    wb.save(path)
    _CACHE[key] = path
    return path


def load(path):
    from xslope.fileio import load_slope_data
    with quiet():
        return load_slope_data(path)


def _fem_mesh(sd, target=12.0):
    """The finite element mesh of ``sd``, built the way main_fem.py builds one."""
    from xslope.mesh import (build_mesh_from_polygons, extract_constraint_line_geometry,
                             extract_joint_options, extract_size_regions,
                             get_material_polygons)
    with quiet():
        lines, _, _ = extract_constraint_line_geometry(sd)
        polys = get_material_polygons(sd, reinf_lines=lines)
        return build_mesh_from_polygons(
            polys, target_size=target, element_type="tri6", lines=lines or None,
            element_size_1d=sd.get("element_size_1d"),
            size_regions=extract_size_regions(sd),
            joint_lines=extract_joint_options(sd))


def toy_model():
    """The toy slope_data with its mesh, a signed hydrostatic seepage field from
    each piezometric line on the mesh's nodes, and seepage boundaries."""
    import numpy as np
    key = ("toy", TEMPLATE)
    if key in _CACHE:
        return _CACHE[key]
    sd = load(toy_path())
    mesh = _fem_mesh(sd)
    nodes = np.asarray(mesh["nodes"], dtype=float)
    gw = float(sd["gamma_water"])
    sd["mesh"] = mesh
    sd["seep_u"] = np.array([gw * (_y_on(_PIEZO, x) - y) for x, y in nodes])
    sd["seep_u2"] = np.array([gw * (_y_on(_PIEZO2, x) - y) for x, y in nodes])
    sd["seepage_bc"] = {"specified_heads": [{"head": 11.0, "kind": "head",
                                             "coords": [(0.0, 0.0), (0.0, 10.0)]}],
                        "specified_fluxes": [],
                        "exit_face": [(20.0, 10.0), (40.0, 20.0)]}
    from xslope.water import seep_water_table_profile
    seep_water_table_profile(sd)         # read once; every copy carries it
    _CACHE[key] = (sd, mesh)
    return _CACHE[key]


def _with(sd, **items):
    """A shallow copy of ``sd`` with ``items`` set and the slicer's caches gone --
    all but the seepage water table, which reads no material and is checked against
    the identity of the field and mesh it was read from (water.py)."""
    out = {k: v for k, v in sd.items()
           if not k.startswith("_") or k == "_water_table_profile"}
    out.update(items)
    return out


# =========================================================================== #
# 5. mat
# =========================================================================== #
#: The mat sheet's switches, in the order a case lists them. Each takes the words
#: its dropdown offers, and blank.
MAT_SWITCHES = ("option", "u", "unsat")

#: The mat sheet's header -> the field key the loader and Studio use.
MAT_KEYS = {"g": "gamma", "gsat": "gamma_sat", "option": "option", "c": "c",
            "f": "phi", "c/p": "cp", "r-elev": "r_elev", "d": "d", "psi": "psi",
            "t_cut": "t_cut", "E": "E", "nu": "nu", "u": "u", "ru": "ru",
            "phi_b": "phi_b", "s_cap": "s_cap", "pow_a": "pow_a", "pow_b": "pow_b",
            "pow_c": "pow_c", "pow_d": "pow_d", "hb_sci": "hb_sci",
            "hb_gsi": "hb_gsi", "hb_mi": "hb_mi", "hb_d": "hb_d",
            "s(g)": "sigma_gamma", "s(c)": "sigma_c", "s(f)": "sigma_phi",
            "s(c/p)": "sigma_cp", "s(d)": "sigma_d", "s(psi)": "sigma_psi",
            "k1": "k1", "k2": "k2", "alpha": "alpha", "unsat": "unsat",
            "kr0": "kr0", "h0": "h0", "a": "vg_a", "n": "vg_n", "l": "vg_l",
            "Ss": "Ss", "Sy": "Sy"}

#: Every numeric cell filled, each with a value its option can use, so a nudge to
#: any of them is a valid model: gsat above g (and above g nudged), the standard
#: deviations below their means, l off the 0.5 a blank reads as.
MAT_VALUES = {"gamma": 120.0, "gamma_sat": 135.0, "c": 200.0, "phi": 30.0,
              "cp": 0.3, "r_elev": 15.0, "d": 150.0, "psi": 20.0, "t_cut": 5.0,
              "E": 1.0e5, "nu": 0.3, "ru": 0.2, "phi_b": 15.0, "s_cap": 300.0,
              "pow_a": 1.5, "pow_b": 0.7, "pow_c": 10.0, "pow_d": 5.0,
              "hb_sci": 1000.0, "hb_gsi": 50.0, "hb_mi": 10.0, "hb_d": 0.2,
              "sigma_gamma": 2.0, "sigma_c": 20.0, "sigma_phi": 2.0,
              "sigma_cp": 0.02, "sigma_d": 10.0, "sigma_psi": 2.0,
              "k1": 0.1, "k2": 0.05, "alpha": 10.0, "kr0": 0.01, "h0": -2.0,
              "vg_a": 0.5, "vg_n": 1.8, "vg_l": 0.6, "Ss": 1.0e-4, "Sy": 0.2}
NUDGE = 1.05


def mat_switch_values(path=None):
    cols = mat_columns(path)
    return {k: [""] + dropdown_values("mat", cols[k], path) for k in MAT_SWITCHES}


def mat_cases(path=None):
    v = mat_switch_values(path)
    return list(itertools.product(*(v[k] for k in MAT_SWITCHES)))


def mat_columns(path=None):
    """Field key -> the mat sheet's column number, from its header row."""
    hdr = template_sheets(path)["_header"]["mat"]
    return {MAT_KEYS[str(v).strip()]: c for c, v in hdr.items()
            if str(v).strip() in MAT_KEYS}


def mat_sheet_gray(case, path=None):
    """The fields the mat sheet grays for a named row with this case's switches."""
    sh = template_sheets(path)["mat"]
    cols = mat_columns(path)
    row = 11
    cells = {(2, row): "soil"}
    for key, col in cols.items():
        if key in MAT_VALUES:
            cells[(col, row)] = MAT_VALUES[key]
    for key, val in zip(MAT_SWITCHES, case):
        cells[(cols[key], row)] = val or None
    return {key for key, col in cols.items() if sh.gray(cells, (col, row))}


def _mat_book():
    """One workbook copy, one load: a mat row for every case (every numeric cell
    filled), then the calibration rows mat_changes reads (a base row, each cell
    nudged by NUDGE, each cell cleared). Returns (case -> material, calibration
    materials), every one named alike so a name is never a difference."""
    import openpyxl
    key = ("mat_book", TEMPLATE)
    if key in _CACHE:
        return _CACHE[key]
    cases = mat_cases()
    cols = mat_columns()
    fields = list(MAT_VALUES)
    rows = [(dict(MAT_VALUES), case) for case in cases]
    cal = ([{}] + [{k: MAT_VALUES[k] * NUDGE} for k in fields]
           + [{k: None} for k in fields])
    rows += [(dict(MAT_VALUES, **over), ("mc", "piezo", "vg")) for over in cal]
    path = os.path.join(_scratch(), "mat.xlsx")
    shutil.copy(TEMPLATE, path)
    wb = openpyxl.load_workbook(path)
    ws = wb["mat"]
    for i, (vals, switches) in enumerate(rows):
        r = 11 + i
        _put(ws, r, 2, f"m{i}")
        for k in fields:
            _put(ws, r, cols[k], vals[k])
        for k, val in zip(MAT_SWITCHES, switches):
            _put(ws, r, cols[k], val or None)
    wb.save(path)
    mats = [dict(m, name="soil") for m in load(path)["materials"]]
    if len(mats) != len(rows):
        raise AssertionError(f"the loader kept {len(mats)} of {len(rows)} materials")
    _CACHE[key] = (dict(zip(cases, mats)), mats[len(cases):])
    return _CACHE[key]


def mat_loaded():
    """case -> the material dict the loader makes of a mat row holding that case."""
    return _mat_book()[0]


# --- the mat probes: each maps a material dict to what that engine is handed ----
def _takes(fn, name):
    import inspect
    return name in inspect.signature(fn).parameters


def _probe_rel(m):
    from xslope.reliability import _reliability_param_info as info
    if not _takes(info, "rapid"):              # no rapid variant to ask about
        return [_attempt(info, [m])]
    return [_attempt(info, [m], rapid=r) for r in (False, True)]


def _probe_seep(m):
    """The seepage arrays as the solvers use them: conductivities and anisotropy as
    built, and the unsaturated parameters only through the relative-conductivity and
    storage laws, evaluated over a range of pressure heads (the solvers read them no
    other way)."""
    import numpy as np
    from xslope import seep
    sd, mesh = toy_model()

    def build():
        d = seep.build_seep_data(mesh, _with(sd, materials=[m]), check_inputs=False)
        p = np.concatenate([np.linspace(-30.0, -0.05, 60), [0.0, 0.5, 3.0]])
        arg = {k: d[k][0] for k in ("kr0_by_mat", "h0_by_mat", "vg_a_by_mat",
                                    "vg_n_by_mat", "vg_l_by_mat", "ss_by_mat",
                                    "sy_by_mat", "unsat_by_mat")}
        model = np.full_like(p, arg["unsat_by_mat"], dtype=int)
        kr = seep.kr_relative_vec(p, arg["kr0_by_mat"], arg["h0_by_mat"],
                                  arg["vg_a_by_mat"], arg["vg_n_by_mat"],
                                  arg["vg_l_by_mat"], model)
        cap = seep.storage_capacity_vec(p, arg["ss_by_mat"], arg["sy_by_mat"],
                                        arg["h0_by_mat"], arg["vg_a_by_mat"],
                                        arg["vg_n_by_mat"], model)
        pot = seep.storage_potential_vec(p, arg["ss_by_mat"], arg["sy_by_mat"],
                                         arg["h0_by_mat"], arg["vg_a_by_mat"],
                                         arg["vg_n_by_mat"], model)
        return {"k1": d["k1_by_mat"], "k2": d["k2_by_mat"],
                "angle": d["angle_by_mat"], "kr": kr, "storage": cap,
                "stored": pot, "unsat_check": (d.get("missing_unsat_params"),
                                               d.get("missing_unsat_detail"))}
    return _attempt(build)


def _probe_lem(m, rapid=False):
    from xslope.slice import generate_slices
    sd, _ = toy_model()
    xo, yo, r = _CIRCLE
    circle = {"Xo": xo, "Yo": yo, "R": r, "Depth": yo - r}

    # A slicer with no rapid switch cut one set of slices for every run, rapid
    # drawdown included.
    kw = {"rapid": True} if rapid and _takes(generate_slices, "rapid") else {}

    def run():
        res = generate_slices(_with(sd, materials=[m]), circle=circle, num_slices=6,
                              debug=False, check_inputs=False, **kw)
        ok, df = res[0], res[1]
        if isinstance(df, tuple):
            df = df[0]
        return ok, df
    return _attempt(run)


def _probe_lem_rapid(m):
    return _probe_lem(m, rapid=True)


def _fem_view(fd):
    """build_fem_data's arrays, with the inputs the solver itself holds out masked
    the way it masks them: an elastic material's tension cutoff (its elements are
    held out of plasticity), and the suction pair taken through the solver's own
    resolution (``_resolve_suction_by_elem``: no angle where a material takes no
    suction credit, and a cap only where there is an angle)."""
    import numpy as np
    from xslope import fem
    out = dict(fd)
    elastic = set(fd.get("elastic_materials") or [])
    out["tension_cutoff_by_material"] = {
        k: v for k, v in (fd.get("tension_cutoff_by_material") or {}).items()
        if k not in elastic}
    em = np.asarray(fd["element_materials"], dtype=int)
    tanphib, scap, _ = fem._resolve_suction_by_elem(fd, None, None, em)
    out.pop("phi_b_by_mat", None)
    out.pop("s_cap_by_mat", None)
    out["suction_tanphib"] = tanphib
    out["suction_cap"] = np.where(np.asarray(tanphib) > 0, scap, np.nan)
    return out


def _probe_fem(m):
    from xslope.fem import build_fem_data
    sd, mesh = toy_model()
    return _attempt(lambda: _fem_view(build_fem_data(_with(sd, materials=[m]), mesh)))


#: The probes, cheapest first: a field one finds read is not nudged by the next.
MAT_PROBES = (("reliability", _probe_rel), ("seepage", _probe_seep),
              ("LEM", _probe_lem), ("LEM rapid drawdown", _probe_lem_rapid),
              ("FEM", _probe_fem))


def mat_changes():
    """field -> its loader changes (see loader_changes), from the calibration rows
    of _mat_book."""
    return loader_changes(_mat_book()[1], list(MAT_VALUES))


def mat_engine():
    """(case, field) -> (verdict, who): what the engine and the checks make of every
    mat field in every case. ``who`` names the probe that read it, or the rules a
    check moved."""
    key = ("mat_engine", TEMPLATE)
    if key in _CACHE:
        return _CACHE[key]
    loaded = mat_loaded()
    changes = mat_changes()
    keys = nudge_keys(changes)
    cases = mat_cases()
    values = mat_switch_values()
    fields = list(MAT_VALUES) + list(MAT_SWITCHES)
    verdict = {}

    def is_read(c, k):
        return verdict.get((c, k), (UNREAD,))[0] == READ

    # Cases the loader reads identically are one case to every probe.
    rep = {}
    for c in cases:
        rep.setdefault(repr(sorted(loaded[c].items())), c)
    canon = {c: rep[repr(sorted(loaded[c].items()))] for c in cases}
    uniq = sorted(set(canon.values()), key=cases.index)

    for pname, probe in MAT_PROBES:
        base = {c: probe(loaded[c]) for c in uniq}
        _must_run(pname, base)
        for members in classes_of(cases, value_classes(uniq, base.get), canon):
            c0 = canon[members[0]]

            def mark(k, members=members):
                for c in members:
                    if not is_read(c, k):
                        verdict[(c, k)] = (READ, [pname])

            # A switch is read where another of its words moves the output.
            for i, k in enumerate(MAT_SWITCHES):
                if not all(is_read(c, k) for c in members) and any(
                        not same(base[c0], base[canon[c0[:i] + (v,) + c0[i + 1:]]])
                        for v in values[k] if v != c0[i]):
                    mark(k)

            def moved(ks, _ctx, c0=c0):
                m = dict(loaded[c0])
                for kk in ks:
                    m.update(apply_change(loaded[c0], changes, kk))
                return not same(base[c0], probe(m))
            todo = [kk for kk in keys if not all(is_read(c, kk[0]) for c in members)]
            for kind in ("nudge", "clear"):
                split_test([kk for kk in todo if kk[1] == kind], moved,
                           lambda kk: None, lambda kk, _, mark=mark: mark(kk[0]))

    # What no engine reads, the checks may.
    sd, mesh = toy_model()

    def model_for(c, ch):
        return _with(sd, materials=[dict(loaded[c], **ch)])
    verdict.update(factored_checks(
        cases, canon, lambda c: [kk for kk in keys if not is_read(c, kk[0])],
        model_for, lambda c, kk: apply_change(loaded[c], changes, kk), mesh))
    for c in cases:
        for k in fields:
            verdict.setdefault((c, k), (UNREAD, []))
    _CACHE[key] = verdict
    return verdict


# --- Studio, materials ---------------------------------------------------------
@contextlib.contextmanager
def studio_sandbox():
    """Studio's settings go to a scratch ini while a dialog is open, through the
    suite's own mechanism (tools/studio_settings.py), which also nests inside a
    caller's redirection (table_paste_check.run is isolated as a whole)."""
    from tools import studio_settings
    with studio_settings.redirected(os.path.join(_scratch(), "studio.ini")):
        yield


def _table_gray(table, i, keys):
    """The keys whose cell in row ``i`` of an _EditableTable is grayed."""
    from PySide6.QtCore import Qt
    from PySide6.QtWidgets import QComboBox
    t = table.table
    out = set()
    for j, key in enumerate(keys):
        w = t.cellWidget(i, j)
        if isinstance(w, QComboBox):
            live = w.isEnabled()
        else:
            it = t.item(i, j)
            live = it is None or bool(it.flags() & Qt.ItemIsEditable)
        if not live:
            out.add(key)
    return out


def mat_studio():
    """case -> {"table": grayed keys, "list": grayed or hidden keys, "shown":
    {u, unsat} as the views display them}, for the row the loader makes of each
    case, read off the real Materials dialog."""
    key = ("mat_studio", TEMPLATE)
    if key in _CACHE:
        return _CACHE[key]
    from studio.editors import MaterialsEditor
    loaded = mat_loaded()
    cases = mat_cases()
    sd, _ = toy_model()
    rows = []
    for i, c in enumerate(cases):
        rows.append(dict(loaded[c], name=f"m{i}"))
    keys = [f.key for f in MaterialsEditor.FIELDS]
    out = {}
    with studio_sandbox():
        dlg = MaterialsEditor().build(_with(sd, materials=rows), None)
        dlg.set_view_mode("table")
        tab = dlg._table
        for i, c in enumerate(cases):
            shown = {}
            for k in ("option", "u", "unsat"):
                w = tab.table.cellWidget(i, keys.index(k))
                shown[k] = w.currentText() if w is not None else None
            out[c] = {"table": _table_gray(tab, i, keys), "shown_table": shown}
        dlg.set_view_mode("list")
        lv = dlg._list_view
        lv._refresh_plots = lambda *a, **k: None     # the plots gray nothing
        sig_of = dict(_sigma_pairs())
        for i, c in enumerate(cases):
            lv._load(i)
            hidden = set(getattr(lv, "_model_hidden", ()))
            gray = set()
            for k in keys:
                cell = lv._cell_widgets.get(k)
                if k in sig_of.values() or cell is None:
                    continue
                if not cell.isEnabled() or k in hidden:
                    gray.add(k)
            # A σ sits beside its value, inside the value's cell: hidden with it,
            # and grayed with it or on its own.
            for base, sk in sig_of.items():
                pair = lv._sigma_widgets.get(base)
                if pair is not None and (not pair[1].isEnabled() or base in hidden):
                    gray.add(sk)
            out[c]["list"] = gray
            out[c]["shown_list"] = {"option": lv._opt_combo.currentText(),
                                    "u": lv._u_combo.currentText(),
                                    "unsat": lv._unsat_combo.currentText()}
        dlg.deleteLater()
    _CACHE[key] = out
    return out


def _sigma_pairs():
    try:
        from studio.editors import _MAT_SIGMA
        return list(_MAT_SIGMA.items())
    except ImportError:
        return [("gamma", "sigma_gamma"), ("c", "sigma_c"), ("phi", "sigma_phi"),
                ("cp", "sigma_cp"), ("d", "sigma_d"), ("psi", "sigma_psi")]


def mat_shown_blank_failures():
    """Studio's u and unsat offer no blank. For a row the sheet leaves blank, the
    loader's reading must be what both views show, and the sheet must gray the
    blank as it grays that reading."""
    loaded, studio = mat_loaded(), mat_studio()
    out = []
    for c in mat_cases():
        for i, sw in enumerate(MAT_SWITCHES):
            if sw == "option" or c[i] != "":
                continue
            read = loaded[c][sw]
            for view in ("shown_table", "shown_list"):
                shown = studio[c][view][sw]
                if str(shown).strip().lower() != str(read).strip().lower():
                    out.append(f"mat {_case_text(c)}: the loader reads the blank {sw} "
                               f"as {read!r}, but the {view[6:]} view shows "
                               f"{shown!r}")
            as_read = c[:i] + (read,) + c[i + 1:]
            a, b = mat_sheet_gray(c), mat_sheet_gray(as_read)
            if a != b:
                out.append(f"mat {_case_text(c)}: the sheet grays a blank {sw} "
                           f"differently from {read!r}, which the loader reads it "
                           f"as: {sorted(a ^ b)}")
    return out


def _case_text(case, names=("option", "u", "unsat")):
    return ", ".join(f"{n}={v or '(blank)'}" for n, v in zip(names, case))


# =========================================================================== #
# 6. circles
# =========================================================================== #
CIRCLE_FIELDS = ("Depth", "Xi", "Yi", "R")
#: A circle whose Intercept point lies on the circle its Radius and Depth state, so
#: every Option describes the same arc.
CIRCLE_VALUES = {"Xo": 35.0, "Yo": 32.0, "Depth": 7.0, "Xi": 20.0, "Yi": 12.0,
                 "R": 25.0}


def dropdown_values(sheet_name, col, path=None):
    """The entries a sheet's dropdown on column ``col`` offers, read from the
    template: an inline list, or the cells a range names."""
    import openpyxl
    sh = template_sheets(path)[sheet_name]
    for ranges, formula in sh.dropdowns:
        if any(c0 <= col <= c1 for c0, r0, c1, r1 in ranges):
            f = str(formula)
            if f.startswith('"'):
                return [w.strip() for w in f.strip('"').split(",")]
            wkey = ("wb", path or TEMPLATE)
            if wkey not in _CACHE:
                _CACHE[wkey] = openpyxl.load_workbook(path or TEMPLATE)
            ws = _CACHE[wkey][sheet_name]
            return [str(c.value).strip() for row in ws[f.replace("$", "")]
                    for c in row if c.value is not None]
    raise ValueError(f"{sheet_name}: no dropdown on column {_col_letters(col)}")


def circle_cases(path=None):
    return [""] + dropdown_values("circles", 4, path)


def circle_sheet_gray(option, path=None):
    sh = template_sheets(path)["circles"]
    hdr = {v: c for c, v in template_sheets(path)["_header"]["circles"].items()}
    cells = {(hdr[k], 3): v for k, v in CIRCLE_VALUES.items()}
    cells[(hdr["Option"], 3)] = option or None
    return {k for k in CIRCLE_FIELDS if sh.gray(cells, (hdr[k], 3))}


def _circle_rows(path, options):
    """Write a base row and one row per nudged field for each Option into a copy of
    the toy model and load it: (the circles, which (option, field) each is), or
    the loader's refusal."""
    import openpyxl
    hdr = {v: c for c, v in template_sheets()["_header"]["circles"].items()}
    shutil.copy(toy_path(), path)
    wb = openpyxl.load_workbook(path)
    ws = wb["circles"]
    rows = []
    for opt in options:
        for nudge in (None,) + CIRCLE_FIELDS:
            vals = dict(CIRCLE_VALUES)
            if nudge:
                vals[nudge] *= NUDGE
            r = 3 + len(rows)
            for k, v in vals.items():
                _put(ws, r, hdr[k], v)
            _put(ws, r, hdr["Option"], opt or None)
            _put(ws, r, 1, len(rows) + 1)
            rows.append((opt, nudge))
    wb.save(path)
    got = _attempt(load, path)
    if got[0] != "ok":
        return None, rows, got
    circles = got[1]["circles"]
    if len(circles) != len(rows):
        raise AssertionError(f"the loader kept {len(circles)} of {len(rows)} circles")
    return circles, rows, None


def circle_engine():
    """(option -> {field: verdict}, option -> the circle the loader makes of the
    base row, or its refusal): whether nudging the field moves the circle the
    loader hands the engines. Everything downstream reads the circle and nothing
    else of the row, so a circle that does not move is a field nothing reads. A
    loader that refuses a whole workbook over one Option is asked again an Option
    at a time, and an Option it refuses reads nothing of its row."""
    key = ("circle_engine", TEMPLATE)
    if key in _CACHE:
        return _CACHE[key]
    cases = circle_cases()
    out, loaded = {}, {}
    circles, rows, refusal = _circle_rows(os.path.join(_scratch(), "circles.xlsx"),
                                          cases)
    batches = [(cases, circles, rows, refusal)] if circles is not None else [
        ((opt,),) + _circle_rows(os.path.join(_scratch(), f"circles_{i}.xlsx"), [opt])
        for i, opt in enumerate(cases)]
    for opts, circles, rows, refusal in batches:
        for opt in opts:
            if circles is None:
                out[opt] = {k: UNREAD for k in CIRCLE_FIELDS}
                loaded[opt] = refusal
                continue
            base = circles[rows.index((opt, None))]
            out[opt] = {k: (UNREAD if same(base, circles[rows.index((opt, k))])
                            else READ) for k in CIRCLE_FIELDS}
            loaded[opt] = base
    _CACHE[key] = (out, loaded)
    return _CACHE[key]


def circle_studio():
    """option -> (grayed fields, the Option the table shows): the Circles table
    with one row per Option, and the row the loader makes of a blank Option."""
    key = ("circle_studio", TEMPLATE)
    if key in _CACHE:
        return _CACHE[key]
    _CACHE[key] = _circle_studio()
    return _CACHE[key]


def _circle_studio():
    from studio.editors import CirclesEditor
    keys = [f.key for f in CirclesEditor.FIELDS]
    cases = circle_cases()
    _, loaded = circle_engine()
    rows = []
    for opt in cases:
        if opt or not isinstance(loaded[opt], dict):
            rows.append(dict(CIRCLE_VALUES, Option=opt or None))
        else:
            rows.append(dict(loaded[opt]))
    sd, _ = toy_model()
    out = {}
    with studio_sandbox():
        dlg = CirclesEditor().build(_with(sd, circles=rows), None)
        table = dlg._editable
        for i, opt in enumerate(cases):
            w = table.table.cellWidget(i, keys.index("Option"))
            shown = w.currentText() if w is not None else \
                table.table.item(i, keys.index("Option")).text()
            out[opt] = (_table_gray(table, i, keys) & set(CIRCLE_FIELDS), shown)
        dlg.deleteLater()
    return out


# =========================================================================== #
# 7. polygon
# =========================================================================== #
# The toy domain split in two material zones, and a small region inside the right
# one for the overlay Types. A material Type is tried on the right zone itself.
_POLY_LEFT = [(0.0, 0.0), (30.0, 0.0), (30.0, 15.0), (20.0, 10.0), (0.0, 10.0)]
_POLY_RIGHT = [(30.0, 0.0), (60.0, 0.0), (60.0, 20.0), (40.0, 20.0), (30.0, 15.0)]
_POLY_INSIDE = [(35.0, 5.0), (50.0, 5.0), (50.0, 12.0), (35.0, 12.0)]
POLY_FIELDS = {"Mat ID": 6, "Size": 8}           # field -> its row in a block
_POLY_VALUES = {"Mat ID": 2, "Size": 3.0}
_POLY_NUDGED = {"Mat ID": 1, "Size": 3.0 * NUDGE}


def polygon_cases(path=None):
    return [""] + dropdown_values("polygon", 2, path)


def polygon_sheet_gray(ptype, path=None):
    """block number -> the fields the polygon sheet grays on that block when its
    Type is ``ptype``: every block's rules, each read at its own Type cell."""
    sh = template_sheets(path)["polygon"]
    out = {}
    for ranges, _ in sh.dropdowns:
        for c0, r0, c1, r1 in ranges:
            for col in range(c0, c1 + 1):
                cells = {(col, r0): ptype or None, (col, 6): 2, (col, 8): 3.0}
                out[_col_letters(col)] = {f for f, r in POLY_FIELDS.items()
                                          if sh.gray(cells, (col, r))}
    return out


def _polygon_model(ptype, nudge=None):
    """The model the loader makes of _polygon_book(ptype, nudge), kept."""
    key = ("poly_model", TEMPLATE, ptype, nudge)
    if key not in _CACHE:
        _CACHE[key] = load(_polygon_book(ptype, nudge))
    return _CACHE[key]


def _polygon_book(ptype, nudge=None):
    """A template copy holding two material zones and the Type under test, with
    its Mat ID and Size (one of them nudged). One template workbook is opened and
    rewritten for each: every cell a book writes, every book writes."""
    import openpyxl
    name = f"poly_{ptype.replace(' ', '_') or 'blank'}_{(nudge or 'base').replace(' ', '_')}"
    path = os.path.join(_scratch(), name + ".xlsx")
    wkey = ("poly_wb", TEMPLATE)
    if wkey not in _CACHE:
        _CACHE[wkey] = openpyxl.load_workbook(TEMPLATE)
    wb = _CACHE[wkey]
    ws = wb["mat"]
    for r, nm in ((11, "left"), (12, "right")):
        for col, v in zip("BCEFGMNO", (nm, 120, "mc", 200, 30, 1e5, 0.3, "none")):
            ws[f"{col}{r}"] = v
    ws = wb["polygon"]
    material = ptype in ("", "material")
    blocks = [(_POLY_LEFT, "material", 1, None)]
    vals = {k: (_POLY_NUDGED[k] if k == nudge else v) for k, v in _POLY_VALUES.items()}
    if material:
        blocks.append((_POLY_RIGHT, ptype, vals["Mat ID"], vals["Size"]))
    else:
        blocks.append((_POLY_RIGHT, "material", 2, None))
        blocks.append((_POLY_INSIDE, ptype, vals["Mat ID"], vals["Size"]))
    for j in range(15):
        x, y = 1 + 3 * j, 2 + 3 * j
        for r in (5, 6, 8):
            _put(ws, r, y, None)
        if j >= len(blocks):
            _put(ws, 4, x, None)
            continue
        for r in range(10, 20):
            _put(ws, r, x, None)
            _put(ws, r, y, None)
        ring, t, mid, size = blocks[j]
        _put(ws, 4, x, f"Polygon #{j + 1}")
        _put(ws, 5, y, t or None)
        _put(ws, 6, y, mid)
        _put(ws, 8, y, size)
        for i, (px, py) in enumerate(ring):
            _put(ws, 10 + i, x, px)
            _put(ws, 10 + i, y, py)
    wb.save(path)
    return path


def _polygon_view(sd):
    """What reaches the engines from the polygon sheet: the zones the mesher and
    the slicer take, the size regions the mesher takes, and the overlays as their
    readers take them (an SSR zone by kind and ring, a joints region by name and
    ring)."""
    from xslope.mesh import extract_size_regions, get_material_polygons
    with quiet():
        return {"zones": get_material_polygons(sd),
                "sizes": extract_size_regions(sd),
                "ssr": [(z.get("kind"), z.get("polygon")) for z in sd.get("ssr_zones") or []],
                "joints": [(z.get("label"), z.get("polygon"))
                           for z in sd.get("joint_zones") or []]}


def polygon_engine():
    """Type -> {field: (verdict, who)}."""
    key = ("polygon_engine", TEMPLATE)
    if key in _CACHE:
        return _CACHE[key]
    out = {}
    for ptype in polygon_cases():
        base_sd = _polygon_model(ptype)
        base = _polygon_view(base_sd)
        out[ptype] = {}
        for f in POLY_FIELDS:
            sd = _polygon_model(ptype, f)
            if not same(base, _polygon_view(sd)):
                out[ptype][f] = (READ, ["mesher/slicer"])
                continue
            diff = _changed(findings(base_sd), findings(sd))
            rows = [x for d in diff.values() for x in d]
            if not rows:
                out[ptype][f] = (UNREAD, [])
            else:
                err = any(s == "error" for _, s, _ in rows)
                out[ptype][f] = (REFUSES if err else NOTES,
                                 sorted({r for r, _, _ in rows}))
    _CACHE[key] = out
    return out


def polygon_studio():
    """Type -> ({field: grayed}, the Type the editor shows): the Polygons editor
    on each Type's workbook, at the polygon under test."""
    key = ("polygon_studio", TEMPLATE)
    if key in _CACHE:
        return _CACHE[key]
    _CACHE[key] = _polygon_studio()
    return _CACHE[key]


def _polygon_studio():
    from studio.editors import PolygonEditor
    out = {}
    for ptype in polygon_cases():
        sd = _polygon_model(ptype)
        ring = _POLY_RIGHT if ptype in ("", "material") else _POLY_INSIDE
        with studio_sandbox():
            dlg = PolygonEditor().build(sd, None)
            idx = [i for i, ln in enumerate(dlg._lines)
                   if [tuple(map(float, c)) for c in ln["coords"]] == ring]
            if len(idx) != 1:
                raise AssertionError(f"Type {ptype!r}: the editor lists the polygon "
                                     f"under test {len(idx)} times")
            dlg.list.setCurrentRow(idx[0])
            gray = {"Mat ID": not dlg.mat_combo.isEnabled(),
                    "Size": not dlg.size_edit.isEnabled()}
            shown = dlg.type_combo.currentText() if dlg.type_combo is not None else None
            dlg.deleteLater()
        out[ptype] = (gray, shown)
    return out


# =========================================================================== #
# 8. reinforce
# =========================================================================== #
REINF_KEYS = {"Type": "type", "Tmax": "t_max", "Lp1": "lp1", "Lp2": "lp2",
              "Adhesion": "adhesion", "Delta": "delta", "Tend1": "tend1",
              "Tend2": "tend2", "Spacing": "spacing", "Tres": "t_res", "E": "E",
              "Area": "area", "Joint": "joint", "kn": "kn", "ks": "ks",
              "Jred": "jred"}
#: A nail from the slope face back into the slope, crossing the toy circle.
REINF_VALUES = {"x1": 30.0, "y1": 15.0, "x2": 55.0, "y2": 12.0, "t_max": 50.0,
                "lp1": 2.0, "lp2": 3.0, "adhesion": 50.0, "delta": 20.0,
                "tend1": 10.0, "tend2": 10.0, "spacing": 2.0, "t_res": 20.0,
                "E": 2.0e7, "area": 0.01, "kn": 1.0e5, "ks": 1.0e5, "jred": "Yes"}
PULLOUT = {"neither": (), "Adhesion only": ("adhesion",),
           "Delta only": ("delta",), "both": ("adhesion", "delta")}
REINF_SWITCHES = ("joint", "pullout", "type")


def reinforce_cases(path=None):
    sh_types = dropdown_values("reinforce", 7, path)
    joints = dropdown_values("reinforce", 21, path)
    return list(itertools.product([""] + joints, list(PULLOUT), [""] + sh_types))


def _reinf_columns(path=None):
    hdr = template_sheets(path)["_header"]["reinforce"]
    out = {}
    for c, v in hdr.items():
        v = str(v).strip()
        if v in ("x1", "y1", "x2", "y2"):
            out[v] = c
        elif v in REINF_KEYS:
            out[REINF_KEYS[v]] = c
    return out


def reinforce_sheet_gray(case, path=None):
    joint, pullout, rtype = case
    sh = template_sheets(path)["reinforce"]
    cols = _reinf_columns(path)
    cells = {(cols[k], 3): v for k, v in REINF_VALUES.items() if k in cols}
    for k in ("adhesion", "delta"):
        if k not in PULLOUT[pullout]:
            cells[(cols[k], 3)] = None
    cells[(cols["joint"], 3)] = joint or None
    cells[(cols["type"], 3)] = rtype or None
    return {k for k, c in cols.items()
            if k in REINF_KEYS.values() and sh.gray(cells, (c, 3))}


def _reinf_rows_book():
    """One workbook copy of the toy model, one load: a reinforce row for every case,
    then the calibration rows reinf_changes reads (a jointed nail with both pullout
    terms; each cell nudged, then each cleared). Returns (case -> line,
    calibration lines, the model)."""
    key = ("reinf_book", TEMPLATE)
    if key in _CACHE:
        return _CACHE[key]
    cases = reinforce_cases()
    fields = [k for k in REINF_KEYS.values() if k not in ("joint", "type")]
    rows = []
    for joint, pullout, rtype in cases:
        row = {k: v for k, v in REINF_VALUES.items()
               if k not in ("adhesion", "delta") or k in PULLOUT[pullout]}
        rows.append(dict(row, joint=joint or None, type=rtype or None))
    base = dict(REINF_VALUES, joint="Yes", type="Nail")
    nudged = {k: ("No" if k == "jred" else REINF_VALUES[k] * NUDGE) for k in fields}
    rows += ([base] + [dict(base, **{k: nudged[k]}) for k in fields]
             + [dict(base, **{k: None}) for k in fields])
    lines, sd = _reinf_book(os.path.join(_scratch(), "reinforce.xlsx"), rows)
    lines = [dict(ln, label="") for ln in lines]
    _CACHE[key] = (dict(zip(cases, lines[:len(cases)])), lines[len(cases):], sd)
    return _CACHE[key]


def reinforce_loaded():
    """(case -> the reinforcement line the loader makes of a reinforce row holding
    that case, the model it was loaded in)."""
    loaded, _, sd = _reinf_rows_book()
    return loaded, sd


def _reinf_book(path, rows):
    """Write reinforce rows (field -> cell value) into a copy of the toy model and
    load it: (the lines, the model)."""
    import openpyxl
    cols = _reinf_columns()
    shutil.copy(toy_path(), path)
    wb = openpyxl.load_workbook(path)
    ws = wb["reinforce"]
    for i, row in enumerate(rows):
        r = 3 + i
        _put(ws, r, 1, i + 1)
        for k, c in cols.items():
            _put(ws, r, c, row.get(k))
    wb.save(path)
    sd = load(path)
    # The loader finishes each line with what it derives from the model (its
    # pullout profile); a line handed on is finished again by _reinf_model.
    lines = [{k: v for k, v in ln.items() if not k.startswith("_")}
             for ln in sd["reinforcement_lines"]]
    if len(lines) != len(rows):
        raise AssertionError(f"the loader kept {len(lines)} of {len(rows)} lines")
    return lines, sd


def reinf_changes():
    """field -> its loader changes, from the calibration rows of _reinf_rows_book."""
    fields = [k for k in REINF_KEYS.values() if k not in ("joint", "type")]
    return loader_changes(_reinf_rows_book()[1], fields)


def _reinf_model(line):
    """The toy model carrying one reinforcement line, finished the way the loader
    and the editor finish one: the LEM's lines and their pullout profiles. Kept by
    the line's content (nothing that reads it changes it)."""
    from xslope.fileio import attach_reinforce_pullout, build_reinforce_lines
    key = ("reinf_model", repr(sorted(line.items())))
    if key not in _CACHE:
        sd, _ = toy_model()
        out = _with(sd, reinforcement_lines=[dict(line)],
                    reinforce_lines=build_reinforce_lines([dict(line)]))
        with quiet():
            attach_reinforce_pullout(out)
        _CACHE[key] = out
    return _CACHE[key]


def _probe_reinf_lem(line):
    from xslope.slice import generate_slices
    xo, yo, r = _CIRCLE
    circle = {"Xo": xo, "Yo": yo, "R": r, "Depth": yo - r}

    def run():
        sd = _reinf_model(line)
        res = generate_slices(_with(sd), circle=circle, num_slices=6, debug=False,
                              check_inputs=False)
        df = res[1][0] if isinstance(res[1], tuple) else res[1]
        return res[0], df, sd["reinforce_lines"]
    return _attempt(run)


def _probe_reinf_fem(line):
    """build_fem_data on a mesh made for this line: a bar, or a joint split when
    its Joint says so. Meshes are kept by what the mesher is handed."""
    from xslope.fem import build_fem_data
    from xslope.mesh import extract_constraint_line_geometry, extract_joint_options

    def run():
        sd = _with(_reinf_model(line))
        # Its own mesh, so not the toy's seepage field (no material reads it: the
        # toy soil's u is none).
        for k in ("mesh", "seep_u", "seep_u2", "_water_table_profile"):
            sd.pop(k, None)
        with quiet():
            sig = repr((extract_constraint_line_geometry(sd),
                        extract_joint_options(sd)))
        meshes = _CACHE.setdefault("reinf_meshes", {})
        if sig not in meshes:
            meshes[sig] = _fem_mesh(sd)
        return _fem_view(build_fem_data(sd, meshes[sig]))
    return _attempt(run)


REINF_PROBES = (("LEM", _probe_reinf_lem), ("FEM", _probe_reinf_fem))


def reinforce_engine():
    """(case, field) -> (verdict, who), probed as the mat fields are. A blank field
    is nudged by filling it, so Adhesion and Delta move a line between the pullout
    states the cases list."""
    key = ("reinf_engine", TEMPLATE)
    if key in _CACHE:
        return _CACHE[key]
    loaded, _ = reinforce_loaded()
    changes = reinf_changes()
    keys = nudge_keys(changes)
    cases = reinforce_cases()
    fields = list(changes)
    switches = {0: "joint", 2: "type"}
    verdict = {}

    def is_read(c, k):
        return verdict.get((c, k), (UNREAD,))[0] == READ

    def changed(c, ks):
        line = dict(loaded[c])
        for kk in ks:
            line.update(apply_change(loaded[c], changes, kk))
        return line

    for pname, probe in REINF_PROBES:
        base = {c: probe(loaded[c]) for c in cases}
        _must_run(pname, base)
        for members in classes_of(cases, value_classes(cases, base.get)):
            c0 = members[0]

            def mark(k, members=members):
                for c in members:
                    if not is_read(c, k):
                        verdict[(c, k)] = (READ, [pname])
            for i, k in switches.items():
                if any(not same(base[c0], base[c0[:i] + (v,) + c0[i + 1:]])
                       for v in {c[i] for c in cases} if v != c0[i]):
                    mark(k)

            def moved(ks, _ctx, c0=c0):
                return not same(base[c0], probe(changed(c0, ks)))
            todo = [kk for kk in keys if not all(is_read(c, kk[0]) for c in members)]
            for kind in ("nudge", "clear"):
                split_test([kk for kk in todo if kk[1] == kind], moved,
                           lambda kk: None, lambda kk, _, mark=mark: mark(kk[0]))

    # What no engine reads, the checks may.
    verdict.update(factored_checks(
        cases, {c: c for c in cases},
        lambda c: [kk for kk in keys if not is_read(c, kk[0])],
        lambda c, ch: _reinf_model(dict(loaded[c], **ch)),
        lambda c, kk: apply_change(loaded[c], changes, kk)))
    for c in cases:
        for k in fields + list(switches.values()):
            verdict.setdefault((c, k), (UNREAD, []))
    _CACHE[key] = verdict
    return verdict


def reinforce_studio():
    """case -> {"table": grayed keys, "list": grayed keys}, from the real
    Reinforcement dialog holding one line per case, as the loader made them."""
    key = ("reinforce_studio", TEMPLATE)
    if key in _CACHE:
        return _CACHE[key]
    _CACHE[key] = _reinforce_studio()
    return _CACHE[key]


def _reinforce_studio():
    from studio.editors import ReinforcementEditor
    loaded, sd = reinforce_loaded()
    cases = reinforce_cases()
    keys = [f.key for f in ReinforcementEditor.FIELDS]
    out = {}
    with studio_sandbox():
        dlg = ReinforcementEditor().build(
            _with(sd, reinforcement_lines=[loaded[c] for c in cases]), None)
        dlg.set_view_mode("table")
        for i, c in enumerate(cases):
            out[c] = {"table": _table_gray(dlg._table, i, keys)}
        dlg.set_view_mode("list")
        lv = dlg._list_view
        for i, c in enumerate(cases):
            lv.list.setCurrentRow(i)
            out[c]["list"] = {k for k, cell in lv._cells.items()
                              if not cell.isEnabled()}
        dlg.deleteLater()
    return out


# =========================================================================== #
# 9. The three readings side by side
# =========================================================================== #
def _row(sheet, case, field, label, sheet_gray, studio, verdict):
    return {"sheet": sheet, "case": case, "field": field, "label": label,
            "sheet_gray": sheet_gray, "studio": studio, "verdict": verdict[0],
            "who": verdict[1]}


def mat_rows():
    """One row per (case, field) of the mat sheet: the sheet's graying, each Studio
    view's, and the engine's verdict."""
    hdr = template_sheets()["_header"]["mat"]
    label = {MAT_KEYS[str(v).strip()]: str(v).strip() for v in hdr.values()
             if str(v).strip() in MAT_KEYS}
    engine, studio = mat_engine(), mat_studio()
    out = []
    for c in mat_cases():
        sheet = mat_sheet_gray(c)
        case = dict(zip(MAT_SWITCHES, c))
        for k in list(MAT_VALUES) + list(MAT_SWITCHES):
            out.append(_row("mat", case, k, label.get(k, k), k in sheet,
                            {"table": k in studio[c]["table"],
                             "list view (grayed or hidden)": k in studio[c]["list"]},
                            engine[(c, k)]))
    return out


def circle_rows():
    engine, _ = circle_engine()
    studio = circle_studio()
    out = []
    for opt in circle_cases():
        sheet = circle_sheet_gray(opt)
        for k in CIRCLE_FIELDS:
            out.append(_row("circles", {"Option": opt}, k, k, k in sheet,
                            {"table": k in studio[opt][0]}, (engine[opt][k], [])))
    return out


def circle_shown_blank_failures():
    """Studio's Option offers no blank: for a row the sheet leaves blank, the
    table must show what the loader makes of it, and gray what the sheet grays."""
    _, loaded = circle_engine()
    studio = circle_studio()
    gray, shown = studio[""]
    rec = loaded[""]
    out = []
    if not isinstance(rec, dict):
        return [f"circles, Option=(blank): the loader refuses the workbook "
                f"({rec[-1]}), so a row the sheet's dropdown leaves blank cannot "
                f"be opened; the Circles table would show it as {shown!r}"]
    sized = not (_blank(rec.get("Depth")) or _blank(rec.get("R")))
    if not sized and str(shown or "").strip():
        out.append(f"circles, Option=(blank): the loader keeps the circle with no "
                   f"Option and no size (Depth and R undefined, refused by the "
                   f"checks), but the Circles table shows Option {shown!r}")
    if shown in circle_cases() and circle_sheet_gray("") != circle_sheet_gray(shown):
        out.append(f"circles, Option=(blank): the sheet grays "
                   f"{sorted(circle_sheet_gray('')) or 'nothing'} on a blank Option "
                   f"and {sorted(circle_sheet_gray(shown))} on {shown!r}, which the "
                   f"Circles table shows it as")
    return out


def polygon_rows():
    engine, studio = polygon_engine(), polygon_studio()
    out = []
    for t in polygon_cases():
        blocks = polygon_sheet_gray(t)
        for k in POLY_FIELDS:
            grays = {b: k in g for b, g in blocks.items()}
            # Every block's own rule, read at its own Type cell, must say the same.
            sheet_gray = grays["B"]
            for b, g in grays.items():
                if g != sheet_gray:
                    sheet_gray = f"block {b} differs"
            out.append(_row("polygon", {"Type": t}, k, k, sheet_gray,
                            {"Polygons editor": studio[t][0][k]}, engine[t][k]))
    return out


def polygon_shown_blank_failures():
    studio = polygon_studio()
    loaded = _polygon_model("")
    out = []
    shown = studio[""][1]
    zones = loaded.get("polygons") or []
    if str(shown).strip().lower() != "material" or len(zones) != 2:
        out.append(f"polygon, Type=(blank): the loader reads {len(zones)} material "
                   f"zones (2 expected), the Polygons editor shows Type {shown!r}")
    if polygon_sheet_gray("") != polygon_sheet_gray("material"):
        out.append("polygon, Type=(blank): the sheet grays a blank Type differently "
                   "from 'material', which the loader reads it as")
    return out


def reinforce_rows():
    hdr = template_sheets()["_header"]["reinforce"]
    label = {REINF_KEYS[str(v).strip()]: str(v).strip() for v in hdr.values()
             if str(v).strip() in REINF_KEYS}
    engine, studio = reinforce_engine(), reinforce_studio()
    out = []
    for c in reinforce_cases():
        sheet = reinforce_sheet_gray(c)
        case = dict(zip(REINF_SWITCHES, c))
        for k in REINF_KEYS.values():
            out.append(_row("reinforce", case, k, label.get(k, k), k in sheet,
                            {"table": k in studio[c]["table"],
                             "list view": k in studio[c]["list"]},
                            engine[(c, k)]))
    return out
