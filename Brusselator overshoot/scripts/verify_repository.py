#!/usr/bin/env python3
"""Lightweight pre-release audit for the Brusselator research repository."""

from __future__ import annotations

import csv
import hashlib
import json
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "docs" / "MANUSCRIPT_DATA_MANIFEST.csv"
REQUIRED = [
    ROOT / "README.md",
    ROOT / "CITATION.cff",
    ROOT / "requirements.txt",
    ROOT / "brusselator_default_config.m",
    ROOT / "solve_brusselator_1d_forced.m",
    ROOT / "docs" / "REPOSITORY_GUIDE.md",
    ROOT / "docs" / "REPRODUCIBILITY.md",
    MANIFEST,
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def git_candidates() -> list[str]:
    result = subprocess.run(
        ["git", "ls-files", "--cached", "--others", "--exclude-standard"],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
    return [line for line in result.stdout.splitlines() if line]


def main() -> int:
    errors: list[str] = []
    warnings: list[str] = []

    for path in REQUIRED:
        if not path.is_file():
            errors.append(f"missing required file: {path.relative_to(ROOT)}")

    for path in ROOT.glob("*.json"):
        try:
            json.loads(path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as exc:
            errors.append(f"invalid JSON {path.name}: {exc}")

    for relative in git_candidates():
        path = ROOT / relative
        if path.is_file() and path.stat().st_size > 10 * 1024 * 1024:
            errors.append(f"Git candidate exceeds 10 MiB: {relative}")
        if path.suffix.lower() in {".mat", ".fig"}:
            errors.append(f"generated binary is a Git candidate: {relative}")

    if MANIFEST.is_file():
        with MANIFEST.open(newline="", encoding="utf-8") as handle:
            for row in csv.DictReader(handle):
                path = ROOT / row["path"]
                if not path.is_file():
                    warnings.append(f"archive file not present locally: {row['path']}")
                    continue
                actual_size = path.stat().st_size
                actual_hash = sha256(path)
                if actual_size != int(row["bytes"]):
                    errors.append(f"size mismatch: {row['path']}")
                if actual_hash != row["sha256"]:
                    errors.append(f"SHA-256 mismatch: {row['path']}")

    for message in warnings:
        print(f"WARNING: {message}")
    for message in errors:
        print(f"ERROR: {message}", file=sys.stderr)

    print(
        f"Audited {len(git_candidates())} Git candidates and "
        f"{sum(1 for _ in csv.DictReader(MANIFEST.open(encoding='utf-8')))} "
        "manuscript-data entries."
    )
    return 1 if errors else 0


if __name__ == "__main__":
    raise SystemExit(main())
