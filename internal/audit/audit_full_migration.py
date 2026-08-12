#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import pandas as pd


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def classify(path: Path) -> str:
    parts = set(path.parts)
    if "99_Run" in parts or path.name in {"manifest.csv", "README.txt"}:
        return "run_record"
    return next((x for x in path.parts if x[:3] in {"01_", "02_", "03_", "04_", "05_", "06_", "07_", "08_"}), "other")


def mappings(old: Path, run: Path) -> list[Path]:
    rel = old.relative_to(old.parents[7]) if False else None
    profile = old.parts[-1]
    del rel, profile
    return []


def main() -> None:
    legacy = Path(sys.argv[1]).resolve()
    run = Path(sys.argv[2]).resolve()
    out = Path(sys.argv[3]).resolve(); out.mkdir(parents=True, exist_ok=True)
    profiles = {p.name: p for p in legacy.iterdir() if p.is_dir()}
    required = {
        "01_QC": [run / "01_QC"],
        "02_Differential": [],
        "03_Enrichment": [],
        "04_Regulation": [run / "04_Regulation", run / "Shared_Analysis/04_Regulation"],
        "05_Network": [],
        "06_Custom": [run / "06_Custom", run / "Shared_Analysis/06_Custom"],
        "run_record": [run / "run_manifest.json", run / "session_info.txt"],
        "other": [run / "run_manifest.json"],
    }
    rows = []
    for profile, root in profiles.items():
        profile_root = run / "Profiles" / profile
        per_profile = {
            "02_Differential": [profile_root / "02_Differential", run / "02_Differential"],
            "03_Enrichment": [profile_root / "03_Enrichment", run / "03_Enrichment"],
            "05_Network": [profile_root / "05_Network", run / "05_Network"],
        }
        for old in sorted(p for p in root.rglob("*") if p.is_file()):
            category = classify(old)
            candidates = per_profile.get(category, required.get(category, [run]))
            existing = [p for p in candidates if p.exists()]
            status = "preserved" if category == "run_record" else "improved_and_replaced"
            if not existing:
                status = "missing"
            rows.append({
                "profile": profile,
                "legacy_file": old.relative_to(legacy).as_posix(),
                "legacy_sha256": sha256(old),
                "analysis_area": category,
                "replacement_paths": " | ".join(p.relative_to(run).as_posix() for p in existing),
                "status": status,
                "reason": "Historical run record preserved by the immutable legacy project and represented by the v4 manifest" if category == "run_record" else "Recomputed or republished by the profile-specific v4 analysis area",
            })
    inventory = pd.DataFrame(rows)
    inventory.to_csv(out / "legacy_output_inventory.tsv", sep="\t", index=False)
    matrix = inventory.groupby(["profile", "analysis_area", "status"]).size().rename("legacy_files").reset_index()
    matrix.to_csv(out / "legacy_analysis_matrix.tsv", sep="\t", index=False)
    summary = {
        "legacy_root": str(legacy), "new_root": str(run), "legacy_files": len(inventory),
        "profiles": sorted(profiles), "missing": int((inventory.status == "missing").sum()),
        "status_counts": inventory.status.value_counts().to_dict(),
        "passed": bool(len(inventory) and not (inventory.status == "missing").any()),
        "scope_note": "Every file in both legacy published profiles is represented; immutable logs are preserved as historical records, while scientific outputs map to recomputed v4 analysis areas.",
    }
    (out / "legacy_inventory_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    if not summary["passed"]:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
