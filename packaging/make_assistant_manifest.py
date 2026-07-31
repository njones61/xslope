"""Validate ``assistant_models.json`` and stage it as a release asset.

The assistant's model recommendations ride the same road as ``latest.json``:
attached to every tagged release and read by the app from the
``releases/latest/download/`` redirect, so the advice can be corrected without
shipping a new build. Unlike ``latest.json`` there is nothing to generate — the
file is hand-curated in ``packaging/assistant_models.json`` — so this script's
only job is to refuse to publish a copy Studio could not read.

Usage::

    python packaging/make_assistant_manifest.py --out artifacts/assistant_models.json

The checks mirror ``studio.ai.models.valid_manifest`` (schema 1, a ``providers``
map) plus the per-provider shape the dialog reads, and are deliberately strict
here: a malformed manifest is ignored silently in the app, which is the right
behaviour at runtime and the wrong one in CI, where nobody would notice.
"""

from __future__ import annotations

import argparse
import json
import os
import sys

#: The only schema the shipped app understands (studio/ai/models.py).
SUPPORTED_SCHEMA = 1
DEFAULT_SOURCE = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                              "assistant_models.json")


def validate(manifest):
    """Every problem with ``manifest``, as a list of strings (empty = good)."""
    problems = []
    if not isinstance(manifest, dict):
        return ["the manifest is not a JSON object"]
    if manifest.get("schema") != SUPPORTED_SCHEMA:
        problems.append(f"schema is {manifest.get('schema')!r}, expected "
                        f"{SUPPORTED_SCHEMA}")
    if not isinstance(manifest.get("updated"), str) or not manifest["updated"]:
        problems.append("'updated' must be a date string")
    providers = manifest.get("providers")
    if not isinstance(providers, dict) or not providers:
        return problems + ["'providers' must be a non-empty object"]

    for name, entry in providers.items():
        where = f"providers.{name}"
        if not isinstance(entry, dict):
            problems.append(f"{where} is not an object")
            continue
        rec = entry.get("recommended")
        if rec is not None and (not isinstance(rec, str) or not rec.strip()):
            problems.append(f"{where}.recommended must be a model id")
        good = entry.get("good_choices", [])
        if not isinstance(good, list):
            problems.append(f"{where}.good_choices must be a list")
            good = []
        ids = []
        for i, row in enumerate(good):
            if not isinstance(row, dict) or not isinstance(row.get("id"), str) \
                    or not row["id"].strip():
                problems.append(f"{where}.good_choices[{i}] needs a string 'id'")
                continue
            ids.append(row["id"])
            for key in ("label", "note"):
                if key in row and not isinstance(row[key], str):
                    problems.append(f"{where}.good_choices[{i}].{key} must be a string")
        dep = entry.get("deprecated", [])
        if not isinstance(dep, list) or any(not isinstance(d, str) for d in dep):
            problems.append(f"{where}.deprecated must be a list of model ids")
            dep = []
        if rec and rec in dep:
            problems.append(f"{where} recommends {rec!r} and also deprecates it")
        for dup in {i for i in ids if ids.count(i) > 1}:
            problems.append(f"{where}.good_choices lists {dup!r} twice")
    return problems


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--source", default=DEFAULT_SOURCE,
                    help="the curated manifest (default: packaging/assistant_models.json)")
    ap.add_argument("--out", default="assistant_models.json",
                    help="where to write the validated copy")
    args = ap.parse_args(argv)

    with open(args.source, encoding="utf-8") as fh:
        text = fh.read()
    try:
        manifest = json.loads(text)
    except Exception as exc:
        raise SystemExit(f"{args.source} is not valid JSON: {exc}")

    problems = validate(manifest)
    if problems:
        for p in problems:
            print(f"error: {p}", file=sys.stderr)
        raise SystemExit(f"{args.source} would not be readable by Studio - "
                         f"refusing to publish it")

    with open(args.out, "w", encoding="utf-8") as fh:
        fh.write(json.dumps(manifest, indent=2) + "\n")
    providers = ", ".join(f"{k} -> {v.get('recommended', '-')}"
                          for k, v in manifest["providers"].items())
    print(f"{args.out}: schema {manifest['schema']}, updated "
          f"{manifest['updated']} ({providers})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
