#!/usr/bin/env python3
"""Atomically publish a complete staging run and preserve the previous result."""
from __future__ import annotations

import argparse
import hashlib
import json
import shutil
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import yaml


MODULES = [
    "01_QC", "02_Differential", "03_Enrichment", "04_Regulation",
    "05_Network", "06_Motif", "07_Exploratory",
]


def sha256(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def ensure_within(path: Path, parent: Path, label: str) -> None:
    try:
        path.relative_to(parent)
    except ValueError as exc:
        raise SystemExit(f"{label} escapes the expected project path: {path}") from exc


def manifest_rows(root: Path) -> list[dict[str, object]]:
    rows = []
    manifest = root / "Provenance/publication_manifest.tsv"
    for path in sorted(root.rglob("*")):
        if path.is_file() and path != manifest:
            rows.append({"published": path.relative_to(root).as_posix(), "bytes": path.stat().st_size, "sha256": sha256(path)})
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source_run", type=Path)
    parser.add_argument("project", type=Path)
    parser.add_argument("--confirm", required=True, help="Exact token PUBLISH:<project_id>:<run_id>")
    args = parser.parse_args()

    source = args.source_run.resolve()
    project = args.project.resolve()
    staging_root = (project / "work/staging").resolve()
    history_root = (project / "work/history/published").resolve()
    destination = (project / "results").resolve()
    ensure_within(source, staging_root, "Source run")
    if not (source / "run_manifest.json").is_file():
        raise SystemExit("Source run lacks run_manifest.json")
    manifest = json.loads((source / "run_manifest.json").read_text(encoding="utf-8"))
    project_id = manifest.get("project_id")
    run_id = manifest.get("run_id") or source.name
    if manifest.get("status") != "complete":
        raise SystemExit("Only a validated complete staging run can be published")
    if project_id != project.name:
        raise SystemExit("Run project_id does not match destination project")
    if args.confirm != f"PUBLISH:{project_id}:{run_id}":
        raise SystemExit("Publication confirmation token does not match the project and run")
    if any(not (source / module / "index.html").is_file() for module in MODULES):
        raise SystemExit("Staging run lacks one or more module indexes")
    project_config = yaml.safe_load((project / "project.yml").read_text(encoding="utf-8"))
    if project_config.get("migration", {}).get("requires_legacy_parity", False):
        audit = project / "work/audits" / run_id / "legacy_inventory_summary.json"
        if not audit.is_file():
            raise SystemExit("Publication blocked: migrated project lacks the required legacy output inventory")
        audit_summary = json.loads(audit.read_text(encoding="utf-8"))
        if not audit_summary.get("passed", False):
            raise SystemExit(
                "Publication blocked: legacy output parity failed "
                f"({audit_summary.get('missing_families', 'unknown')} output families missing)"
            )

    history_root.mkdir(parents=True, exist_ok=True)
    archived = None
    if destination.exists():
        current_manifest = destination / "run_manifest.json"
        if not current_manifest.is_file():
            current_manifest = destination / "Provenance/run_manifest.json"
        previous_id = "pre-platform-results"
        if current_manifest.is_file():
            previous = json.loads(current_manifest.read_text(encoding="utf-8"))
            previous_id = previous.get("run_id") or previous.get("project_id") or previous_id
        archived = history_root / previous_id
        if archived.exists():
            timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
            archived = history_root / f"{previous_id}-{timestamp}"
        destination.replace(archived)

    try:
        source.replace(destination)
    except Exception:
        if archived and archived.exists() and not destination.exists():
            archived.replace(destination)
        raise

    provenance = destination / "Provenance"
    provenance.mkdir(parents=True, exist_ok=True)
    publication = pd.DataFrame(manifest_rows(destination))
    publication.to_csv(provenance / "publication_manifest.tsv", sep="\t", index=False)
    record = {
        "published_at": datetime.now(timezone.utc).isoformat(),
        "project_id": project_id,
        "run_id": run_id,
        "archived_previous_results": str(archived) if archived else None,
        "current_results": str(destination),
    }
    (provenance / "publication_record.json").write_text(json.dumps(record, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    # Refresh manifest after adding the publication record itself.
    pd.DataFrame(manifest_rows(destination)).to_csv(provenance / "publication_manifest.tsv", sep="\t", index=False)
    print(json.dumps(record, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
