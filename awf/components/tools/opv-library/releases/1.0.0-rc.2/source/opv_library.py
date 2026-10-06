#!/usr/bin/env python3
"""Deterministic command-line access to the project OPV resource library.

This helper deliberately stays project-local.  It indexes and validates the
source-preserving library, and emits a reviewable ztop command plan; it never
executes ztop or changes a scientific input by itself.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import shlex
import sys
from pathlib import Path


def discover_bundled_library() -> Path:
    relative = Path("components/artifacts/opv-fragment-library/releases/1.0.0-rc.2/source/opv_fragment_library")
    for parent in Path(__file__).resolve().parents:
        candidate = parent / relative
        if candidate.is_dir():
            return candidate
    return Path("opv_fragment_library")


DEFAULT_LIBRARY = discover_bundled_library()
DEFAULT_ZTOP = Path(os.environ.get("ZTOP", "ztop.py"))


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        raise FileNotFoundError(path)
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_rows(library: Path) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    return (
        read_tsv(library / "index" / "fragment_catalog.tsv"),
        read_tsv(library / "index" / "molecule_catalog.tsv"),
    )


def parse_fragment_line(line: str) -> dict[str, str] | None:
    stripped = line.strip()
    if not stripped or stripped.startswith("#"):
        return None
    fields = [part.strip() for part in stripped.split(";")]
    values: dict[str, str] = {"slot": fields[0]}
    for field in fields[1:]:
        if "=" in field:
            key, value = field.split("=", 1)
            values[key.strip()] = value.strip()
    return values


def fragment_line(library: Path, row: dict[str, str]) -> str:
    path = library / row["frg_file"]
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    number = int(row["line"])
    if number < 1 or number > len(lines):
        raise ValueError(f"catalog line is outside {path}: {number}")
    line = lines[number - 1].strip()
    parsed = parse_fragment_line(line)
    if parsed is None:
        raise ValueError(f"catalog points to a blank/comment line: {path}:{number}")
    return line


def find_fragment(rows: list[dict[str, str]], fragment_id: str) -> dict[str, str]:
    matches = [row for row in rows if row["fragment_id"] == fragment_id]
    if not matches:
        raise KeyError(f"unknown fragment_id: {fragment_id}")
    return matches[0]


def bool_value(value: str) -> bool:
    return value.strip().lower() in {"1", "true", "yes", "y"}


def verify_hashes(library: Path) -> tuple[int, list[str]]:
    manifest = library / "index" / "SHA256SUMS"
    if not manifest.is_file():
        return 0, [f"missing checksum manifest: {manifest}"]
    checked = 0
    errors: list[str] = []
    for line_number, line in enumerate(manifest.read_text(encoding="utf-8").splitlines(), 1):
        if not line.strip():
            continue
        try:
            expected, relative = line.split(None, 1)
        except ValueError:
            errors.append(f"malformed checksum line {line_number}")
            continue
        path = library / relative
        if not path.is_file():
            errors.append(f"checksum target missing: {relative}")
            continue
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        checked += 1
        if digest != expected:
            errors.append(f"checksum mismatch: {relative}")
    return checked, errors


def fragment_diagnostics(library: Path, row: dict[str, str]) -> list[str]:
    source_dir = library / "sources" / row["source_id"]
    frg_path = library / row["frg_file"]
    messages: list[str] = []
    if not frg_path.is_file():
        messages.append(f"missing frg file: {row['frg_file']}")
        return messages
    line = fragment_line(library, row)
    parsed = parse_fragment_line(line) or {}
    for key, label in (("p", "topology"), ("x", "coordinate")):
        reference = Path(parsed.get(key, "")).name
        if not reference:
            messages.append(f"missing {key}= reference in {row['frg_file']}:{row['line']}")
        elif not (source_dir / reference).is_file():
            messages.append(f"missing {label} reference in source set: {reference}")
    if parsed.get("slot", "") != row.get("slot", ""):
        messages.append("catalog slot disagrees with source line")
    if parsed.get("site", "") != row.get("site", ""):
        messages.append("catalog site disagrees with source line")
    return messages


def molecule_diagnostics(library: Path, row: dict[str, str]) -> list[str]:
    """Check the minimum coordinate/topology pair without judging chemistry."""
    resource_dir = library / row["resource_dir"]
    if not resource_dir.is_dir():
        return [f"missing resource directory: {row['resource_dir']}"]
    coordinates = sorted(path.name for path in resource_dir.iterdir()
                         if path.is_file() and path.suffix.lower() in {".pdb", ".gro"})
    topologies = sorted(path.name for path in resource_dir.iterdir()
                        if path.is_file() and path.suffix.lower() in {".top", ".itp"})
    errors: list[str] = []
    if not coordinates:
        errors.append("no coordinate file (.pdb/.gro)")
    if not topologies:
        errors.append("no topology file (.top/.itp)")
    return errors


def print_rows(rows: list[dict[str, str]], fields: list[str]) -> None:
    if not rows:
        print("(no matches)")
        return
    print("\t".join(fields))
    for row in rows:
        print("\t".join(row.get(field, "") for field in fields))


def command_list(args: argparse.Namespace) -> int:
    fragments, molecules = read_rows(args.library)
    if args.kind == "fragment":
        rows = fragments
        if args.category:
            rows = [row for row in rows if row["category"] == args.category]
        if args.status:
            rows = [row for row in rows if row["status"] == args.status]
        fields = ["fragment_id", "category", "slot", "resname", "site", "status"]
    else:
        rows = molecules
        if args.category:
            rows = [row for row in rows if row["category"] == args.category]
        if args.status:
            rows = [row for row in rows if row["status"] == args.status]
        fields = ["molecule_id", "category", "resource_dir", "status"]
    print_rows(rows[: args.limit] if args.limit else rows, fields)
    return 0


def command_search(args: argparse.Namespace) -> int:
    fragments, molecules = read_rows(args.library)
    pattern = re.compile(args.query, re.IGNORECASE)
    rows: list[dict[str, str]] = []
    if args.kind in {"all", "fragment"}:
        rows.extend(row for row in fragments if pattern.search(" ".join(row.values())))
    if args.kind in {"all", "molecule"}:
        rows.extend(row for row in molecules if pattern.search(" ".join(row.values())))
    if args.kind == "molecule":
        fields = ["molecule_id", "category", "resource_dir", "status"]
    else:
        fields = ["fragment_id", "category", "slot", "resname", "site", "status"]
    print_rows(rows[: args.limit] if args.limit else rows, fields)
    return 0


def command_check(args: argparse.Namespace) -> int:
    fragments, _ = read_rows(args.library)
    row = find_fragment(fragments, args.fragment_id)
    messages = fragment_diagnostics(args.library, row)
    print(json.dumps({"fragment_id": args.fragment_id, "source_id": row["source_id"],
                      "availability": row["availability"], "status": row["status"],
                      "errors": messages}, ensure_ascii=False, indent=2))
    return 1 if messages else 0


def command_check_molecule(args: argparse.Namespace) -> int:
    _, molecules = read_rows(args.library)
    matches = [row for row in molecules if row["molecule_id"] == args.molecule_id]
    if not matches:
        print(f"ERROR: unknown molecule_id: {args.molecule_id}", file=sys.stderr)
        return 2
    row = matches[0]
    messages = molecule_diagnostics(args.library, row)
    print(json.dumps({"molecule_id": args.molecule_id, "category": row["category"],
                      "resource_dir": row["resource_dir"], "status": row["status"],
                      "errors": messages}, ensure_ascii=False, indent=2))
    return 1 if messages else 0


def command_verify(args: argparse.Namespace) -> int:
    fragments, molecules = read_rows(args.library)
    errors: list[str] = []
    checked_fragments = 0
    for row in fragments:
        checked_fragments += 1
        errors.extend(f"{row['fragment_id']}: {message}" for message in fragment_diagnostics(args.library, row))
    checked_molecules = 0
    for row in molecules:
        checked_molecules += 1
        resource_dir = args.library / row["resource_dir"]
        errors.extend(f"{row['molecule_id']}: {message}"
                      for message in molecule_diagnostics(args.library, row))
    hash_count, hash_errors = verify_hashes(args.library)
    errors.extend(hash_errors)
    result = {"fragment_records_checked": checked_fragments,
              "molecules_checked": checked_molecules,
              "checksums_checked": hash_count,
              "errors": errors,
              "status": "PASS" if not errors else "FAIL"}
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return 0 if not errors else 1


def command_plan(args: argparse.Namespace) -> int:
    fragments, _ = read_rows(args.library)
    selected = [find_fragment(fragments, fragment_id) for fragment_id in args.fragment_id]
    source_ids = {row["source_id"] for row in selected}
    if len(source_ids) != 1:
        print("ERROR: selected fragments come from different source sets; choose one source set", file=sys.stderr)
        return 2
    errors: list[str] = []
    for row in selected:
        errors.extend(f"{row['fragment_id']}: {message}" for message in fragment_diagnostics(args.library, row))
    if errors:
        print("ERROR: selected fragment references are not usable:", file=sys.stderr)
        print("\n".join(f"- {error}" for error in errors), file=sys.stderr)
        return 2
    source_id = selected[0]["source_id"]
    source_dir = args.library / "sources" / source_id
    fragment_args = []
    for row in selected:
        fragment_args.extend(["-f", fragment_line(args.library, row)])
    argv = [args.python, str(args.ztop), *fragment_args, "-b", args.build, "-o", args.output]
    plan = {
        "status": "PLAN_ONLY",
        "source_id": source_id,
        "working_directory": str(source_dir),
        "selected_fragment_ids": [row["fragment_id"] for row in selected],
        "build": args.build,
        "output": args.output,
        "argv": argv,
        "shell_command": " ".join(shlex.quote(item) for item in argv),
        "execution": "not executed by this helper; review the plan, then run it from working_directory",
    }
    text = json.dumps(plan, ensure_ascii=False, indent=2)
    if args.plan_file:
        args.plan_file.write_text(text + "\n", encoding="utf-8")
        print(f"Wrote {args.plan_file}")
    else:
        print(text)
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Inspect and plan use of the bundled OPV fragment library")
    parser.add_argument("--library", type=Path, default=DEFAULT_LIBRARY)
    subparsers = parser.add_subparsers(dest="command", required=True)

    list_parser = subparsers.add_parser("list", help="list catalog entries")
    list_parser.add_argument("kind", choices=("fragment", "molecule"))
    list_parser.add_argument("--category")
    list_parser.add_argument("--status")
    list_parser.add_argument("--limit", type=int)
    list_parser.set_defaults(func=command_list)

    search_parser = subparsers.add_parser("search", help="search both catalogs")
    search_parser.add_argument("query")
    search_parser.add_argument("--kind", choices=("all", "fragment", "molecule"), default="all")
    search_parser.add_argument("--limit", type=int)
    search_parser.set_defaults(func=command_search)

    check_parser = subparsers.add_parser("check-fragment", help="check one fragment and its references")
    check_parser.add_argument("fragment_id")
    check_parser.set_defaults(func=command_check)

    molecule_check_parser = subparsers.add_parser(
        "check-molecule", help="check one complete molecule package"
    )
    molecule_check_parser.add_argument("molecule_id")
    molecule_check_parser.set_defaults(func=command_check_molecule)

    verify_parser = subparsers.add_parser("verify", help="verify all references, molecules and checksums")
    verify_parser.set_defaults(func=command_verify)

    plan_parser = subparsers.add_parser("plan-ztop", help="emit a reviewable ztop command plan")
    plan_parser.add_argument("--fragment-id", action="append", required=True)
    plan_parser.add_argument("--build", required=True, help="ztop -b build expression, e.g. D[A]4")
    plan_parser.add_argument("--output", default="assembled.top,assembled.gro")
    plan_parser.add_argument("--ztop", type=Path, default=DEFAULT_ZTOP)
    plan_parser.add_argument("--python", default=sys.executable)
    plan_parser.add_argument("--plan-file", type=Path)
    plan_parser.set_defaults(func=command_plan)
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    args.library = args.library.resolve()
    try:
        return args.func(args)
    except (FileNotFoundError, KeyError, ValueError, OSError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
