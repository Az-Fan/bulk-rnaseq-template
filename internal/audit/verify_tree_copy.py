#!/usr/bin/env python3
"""Prove that every file in a read-only source tree was restored byte-for-byte."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("restored", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()

    source = args.source.resolve()
    restored = args.restored.resolve()
    output = args.output.resolve()
    if not source.is_dir() or not restored.is_dir():
        raise SystemExit("Both source and restored arguments must be existing directories")

    source_files = {p.relative_to(source).as_posix(): p for p in source.rglob("*") if p.is_file()}
    restored_files = {p.relative_to(restored).as_posix(): p for p in restored.rglob("*") if p.is_file()}
    rows = []
    for relative in sorted(set(source_files) | set(restored_files)):
        old = source_files.get(relative)
        new = restored_files.get(relative)
        old_hash = sha256(old) if old else ""
        new_hash = sha256(new) if new else ""
        status = "identical" if old and new and old_hash == new_hash else (
            "missing_from_restored" if old and not new else
            "extra_in_restored" if new and not old else "hash_mismatch"
        )
        rows.append({
            "relative_path": relative,
            "source_bytes": old.stat().st_size if old else "",
            "restored_bytes": new.stat().st_size if new else "",
            "source_sha256": old_hash,
            "restored_sha256": new_hash,
            "status": status,
        })

    output.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0]) if rows else ["relative_path", "source_bytes", "restored_bytes", "source_sha256", "restored_sha256", "status"]
    with (output / "legacy_restore_manifest.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    failures = [row for row in rows if row["status"] != "identical"]
    summary = {
        "checked_at": datetime.now(timezone.utc).isoformat(),
        "source": str(source),
        "restored": str(restored),
        "source_files": len(source_files),
        "restored_files": len(restored_files),
        "identical_files": len(rows) - len(failures),
        "differences": len(failures),
        "passed": not failures,
    }
    (output / "legacy_restore_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
