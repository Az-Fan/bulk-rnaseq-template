#!/usr/bin/env python3
"""Combine all independently generated migration checks into one release gate."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--audit-dir", type=Path, required=True)
    ap.add_argument("--run-dir", type=Path, required=True)
    args = ap.parse_args()
    audit = args.audit_dir
    inventory = json.loads((audit / "legacy_inventory_summary.json").read_text())
    de = json.loads((audit / "de/de_regression_summary.json").read_text())
    enrichment = pd.read_csv(audit / "enrichment_regression.tsv", sep="\t")
    regulation = json.loads((audit / "regulation_regression_summary.json").read_text())
    extended = json.loads((audit / "extended_regression_summary.json").read_text())
    motif_summary = args.run_dir / "06_Motif/Promoter/Tables/Motif_Run_Summary.tsv"
    checks = [
        {"check": "legacy_output_file_inventory", "passed": bool(inventory["passed"]), "detail": f"{inventory.get('legacy_files_in_inventory', 0)} legacy files; {inventory.get('missing_files', 0)} missing"},
        {"check": "differential_expression", "passed": bool(de["passed"]), "detail": f"old/new significant genes: {de['old_significant']}/{de['new_significant']}"},
        {"check": "ora_and_gsea", "passed": bool(enrichment.passed.all()), "detail": f"{len(enrichment)} database-level comparisons"},
        {"check": "tf_progeny_gsva", "passed": all(item["passed"] for item in regulation.values()), "detail": "CollecTRI ULM, PROGENy and GSVA"},
        {"check": "tf_gsea_custom_gsea_ppi", "passed": all(item["passed"] for item in extended.values()), "detail": "ChEA, ENCODE, GTRD, custom GSEA and PPI or its explicit historical skip"},
        {"check": "motif_extension", "passed": motif_summary.is_file(), "detail": "New promoter motif output; no legacy result is replaced"},
    ]
    passed = all(item["passed"] for item in checks)
    pd.DataFrame(checks).to_csv(audit / "legacy_regression_results.tsv", sep="\t", index=False)
    summary = {"passed": passed, "run_dir": str(args.run_dir), "checks": checks,
               "tolerance_note": "Exact ID/direction/class sets are required. Declared numerical tolerances cover floating-point/package-version differences only: GSEA ES/NES 5e-7, GSEA p/FDR 1e-6, ULM 1e-6, PROGENy 1e-7, GSVA 4e-3 (legacy GSVA 2.0.4 versus locked GSVA 2.4.4; signs and top-variable pathways must also match), PPI 1e-9."}
    (audit / "legacy_regression_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    project_label = args.run_dir.parent.parent.name
    lines = [f"# {project_label} legacy regression report", "", f"Overall: **{'PASS' if passed else 'FAIL'}**", "", summary["tolerance_note"], ""]
    lines += [f"- {'PASS' if x['passed'] else 'FAIL'} — {x['check']}: {x['detail']}" for x in checks]
    (audit / "legacy_regression_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))
    return 0 if passed else 2


if __name__ == "__main__":
    raise SystemExit(main())
