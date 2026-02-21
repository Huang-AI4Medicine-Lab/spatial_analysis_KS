#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
from dataclasses import dataclass
from pathlib import Path


DATA_REF_RE = re.compile(r"""(?P<quote>["'])(?P<path>\.\./data/[^"'\n]+)(?P=quote)""")
ALL_SRC_PATTERNS = ("*.ipynb", "*.py", "*.R")
NOTEBOOK_PATTERN = ("*.ipynb",)


@dataclass(frozen=True)
class Ref:
    source_file: Path
    line_number: int
    raw_path: str

    @property
    def is_templated(self) -> bool:
        if "{" in self.raw_path or "}" in self.raw_path:
            return True
        name = Path(self.raw_path.rstrip("/")).name
        if name.endswith("_"):
            return True
        # R paste0-style fragments often appear as prefix strings with no extension.
        if "." not in name:
            return True
        return False

    def resolved_path(self) -> Path:
        return (self.source_file.parent / self.raw_path).resolve()


def scan_refs(src_dir: Path, patterns: tuple[str, ...]) -> list[Ref]:
    refs: list[Ref] = []
    for pattern in patterns:
        for file_path in sorted(src_dir.rglob(pattern)):
            with file_path.open("r", encoding="utf-8", errors="replace") as fh:
                for line_number, line in enumerate(fh, start=1):
                    for match in DATA_REF_RE.finditer(line):
                        raw_path = match.group("path").rstrip("\\")
                        refs.append(
                            Ref(
                                source_file=file_path,
                                line_number=line_number,
                                raw_path=raw_path,
                            )
                        )
    return refs


def load_manifest(manifest_path: Path) -> dict:
    with manifest_path.open("r", encoding="utf-8") as fh:
        return json.load(fh)


def md5sum(file_path: Path, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.md5()
    with file_path.open("rb") as fh:
        while True:
            chunk = fh.read(chunk_size)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


def rel_to_repo(path: Path, repo_root: Path) -> str:
    try:
        return path.relative_to(repo_root).as_posix()
    except ValueError:
        return str(path)


def main() -> int:
    parser = argparse.ArgumentParser(description="Validate ../data/ path usage against local files and Zenodo manifest.")
    parser.add_argument("--repo-root", default=".", help="Repository root directory (default: current directory)")
    parser.add_argument(
        "--manifest",
        default="data/zenodo_17611373_manifest.json",
        help="Path to Zenodo manifest JSON, relative to repo root",
    )
    parser.add_argument(
        "--strict-manifest",
        action="store_true",
        help="Fail when a literal referenced file basename is not present in the Zenodo manifest",
    )
    parser.add_argument(
        "--require-complete-manifest",
        action="store_true",
        help="Fail when any manifest file is missing in local data/",
    )
    parser.add_argument(
        "--check-md5",
        action="store_true",
        help="Compute and verify MD5 for local files that exist and are listed in the manifest",
    )
    parser.add_argument(
        "--manifest-only",
        action="store_true",
        help="Only check references against Zenodo manifest; skip all local data/ existence/completeness checks",
    )
    parser.add_argument(
        "--notebooks-only",
        action="store_true",
        help="Scan only src/*.ipynb files (skip .py/.R)",
    )
    args = parser.parse_args()

    repo_root = Path(args.repo_root).resolve()
    src_dir = repo_root / "src"
    manifest_path = (repo_root / args.manifest).resolve()
    patterns = NOTEBOOK_PATTERN if args.notebooks_only else ALL_SRC_PATTERNS

    if args.manifest_only and (args.check_md5 or args.require_complete_manifest):
        print("ERROR: --manifest-only cannot be combined with --check-md5 or --require-complete-manifest")
        return 2

    if not src_dir.exists():
        print(f"ERROR: src directory not found: {src_dir}")
        return 2
    if not manifest_path.exists():
        print(f"ERROR: manifest not found: {manifest_path}")
        return 2

    manifest = load_manifest(manifest_path)
    manifest_files = manifest.get("files", [])
    manifest_by_name = {item["filename"]: item for item in manifest_files}
    manifest_names = set(manifest_by_name)

    refs = scan_refs(src_dir, patterns)
    literal_refs = [r for r in refs if not r.is_templated]
    templated_refs = [r for r in refs if r.is_templated]

    missing_local: list[Ref] = []
    unknown_manifest_refs: list[Ref] = []

    for ref in literal_refs:
        resolved = ref.resolved_path()
        if not args.manifest_only and not resolved.exists():
            missing_local.append(ref)
        if resolved.name not in manifest_names:
            unknown_manifest_refs.append(ref)

    data_local_files: dict[str, Path] = {}
    missing_manifest_locally: list[str] = []
    md5_mismatches: list[tuple[str, str, str]] = []

    if not args.manifest_only:
        data_dir = repo_root / "data"
        data_local_files = {
            p.name: p
            for p in data_dir.iterdir()
            if p.is_file() and p.name not in {"README.md", "data.txt", manifest_path.name}
        }
        missing_manifest_locally = sorted(name for name in manifest_names if name not in data_local_files)

    if args.check_md5 and not args.manifest_only:
        for name, file_path in sorted(data_local_files.items()):
            if name not in manifest_by_name:
                continue
            expected = manifest_by_name[name].get("checksum_md5", "").lower()
            if not expected:
                continue
            got = md5sum(file_path).lower()
            if got != expected:
                md5_mismatches.append((name, expected, got))

    print(f"Scanned refs: {len(refs)} total ({len(literal_refs)} literal, {len(templated_refs)} templated)")
    print(f"Zenodo manifest files: {len(manifest_names)}")

    if templated_refs:
        print("\nTemplated refs (not existence-checked):")
        for ref in templated_refs:
            print(f"  - {rel_to_repo(ref.source_file, repo_root)}:{ref.line_number} -> {ref.raw_path}")

    if not args.manifest_only:
        if missing_local:
            print("\nMissing local files for literal refs:")
            for ref in missing_local:
                print(f"  - {rel_to_repo(ref.source_file, repo_root)}:{ref.line_number} -> {ref.raw_path}")
        else:
            print("\nLocal existence check: OK (all literal refs exist).")

    if unknown_manifest_refs:
        print("\nLiteral refs not found in Zenodo manifest (basename match):")
        for ref in unknown_manifest_refs:
            print(f"  - {rel_to_repo(ref.source_file, repo_root)}:{ref.line_number} -> {ref.raw_path}")
    else:
        print("\nManifest cross-check: OK (all literal ref basenames are in Zenodo manifest).")

    if not args.manifest_only:
        if missing_manifest_locally:
            print("\nManifest files missing locally in data/:")
            for name in missing_manifest_locally:
                print(f"  - {name}")
        else:
            print("\nManifest completeness: OK (all manifest files present locally).")

    if args.check_md5 and not args.manifest_only:
        if md5_mismatches:
            print("\nMD5 mismatches:")
            for name, expected, got in md5_mismatches:
                print(f"  - {name}: expected {expected}, got {got}")
        else:
            print("\nMD5 check: OK (all checked files matched).")

    should_fail = bool(missing_local) if not args.manifest_only else False
    if args.strict_manifest:
        should_fail = should_fail or bool(unknown_manifest_refs)
    if args.require_complete_manifest and not args.manifest_only:
        should_fail = should_fail or bool(missing_manifest_locally)
    if args.check_md5 and not args.manifest_only:
        should_fail = should_fail or bool(md5_mismatches)

    return 1 if should_fail else 0


if __name__ == "__main__":
    raise SystemExit(main())
