#!/usr/bin/env python3
"""Regression checks for legacy TF GSEA, custom GSEA and PPI outputs."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


def read(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t" if path.suffix in {".tsv", ".gz"} else ",")


def gsea(old_path: Path, new_path: Path, atol: float = 5e-7) -> dict:
    old, new = read(old_path), read(new_path)
    merged = old.merge(new, on="ID", suffixes=("_old", "_new"))
    result = {
        "old_rows": len(old), "new_rows": len(new),
        "id_set_equal": set(old.ID.astype(str)) == set(new.ID.astype(str)),
        "direction_equal": bool(np.array_equal(np.sign(merged.NES_old), np.sign(merged.NES_new))),
        "fdr_class_equal": bool(np.array_equal(merged["p.adjust_old"].lt(0.05), merged["p.adjust_new"].lt(0.05))),
        "absolute_tolerance": atol,
    }
    for column in ("enrichmentScore", "NES", "pvalue", "p.adjust", "qvalue", "setSize"):
        a = pd.to_numeric(merged[f"{column}_old"], errors="coerce")
        b = pd.to_numeric(merged[f"{column}_new"], errors="coerce")
        result[f"{column}_max_abs_diff"] = float(np.nanmax(np.abs(a - b)))
        column_atol = atol if column in {"enrichmentScore", "NES"} else (1e-6 if column in {"pvalue", "p.adjust", "qvalue"} else 1e-10)
        result[f"{column}_tolerance"] = column_atol
        result[f"{column}_allclose"] = bool(np.allclose(a, b, atol=column_atol, rtol=1e-8, equal_nan=True))
    result["passed"] = result["id_set_equal"] and result["direction_equal"] and result["fdr_class_equal"] and all(
        value for key, value in result.items() if key.endswith("_allclose"))
    return result


def ppi(old_path: Path, new_path: Path) -> dict:
    old, new = read(old_path), read(new_path)
    merged = old.merge(new, on="Symbol", suffixes=("_old", "_new"))
    result = {
        "old_rows": len(old), "new_rows": len(new),
        "node_set_equal": set(old.Symbol.astype(str)) == set(new.Symbol.astype(str)),
    }
    for column in ("Degree", "Betweenness", "Closeness", "PageRank", "Eigenvector", "Module", "HubScore"):
        a = pd.to_numeric(merged[f"{column}_old"], errors="coerce")
        b = pd.to_numeric(merged[f"{column}_new"], errors="coerce")
        result[f"{column}_max_abs_diff"] = float(np.nanmax(np.abs(a - b)))
        result[f"{column}_allclose"] = bool(np.allclose(a, b, atol=1e-9, rtol=1e-8, equal_nan=True))
    result["passed"] = result["node_set_equal"] and all(
        value for key, value in result.items() if key.endswith("_allclose"))
    return result


def ppi_explicit_skip(old_network: Path, new_network: Path) -> dict:
    legacy_markers = list(old_network.rglob("SKIPPED_ppi_due_to_network_or_data.txt"))
    status_path = new_network / "Provenance/05_Network_status.tsv"
    status = read(status_path) if status_path.is_file() else pd.DataFrame()
    new_skipped = (
        not status.empty
        and "status" in status.columns
        and set(status["status"].astype(str)).issubset({"skipped_by_user", "not_applicable"})
    )
    return {
        "legacy_explicit_skip_marker": bool(legacy_markers),
        "new_explicit_skip_status": bool(new_skipped),
        "passed": bool(legacy_markers) and bool(new_skipped),
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--legacy", type=Path, required=True)
    ap.add_argument("--new", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--contrast", required=True)
    args = ap.parse_args()
    old_tf = args.legacy / "04_Regulation/TF_Target_GSEA"
    new_tf = args.new / "04_Regulation/TF_Target_GSEA"
    results = {
        "tf_gsea_chea": gsea(old_tf / "ChEA/Tables/GSEA_Full_Table_ChEA_Consensus.csv", new_tf / "ChEA/Tables/GSEA_Full_Table_ChEA_Consensus.tsv"),
        "tf_gsea_encode": gsea(old_tf / "ENCODE/Tables/GSEA_Full_Table_ENCODE_2015.csv", new_tf / "ENCODE/Tables/GSEA_Full_Table_ENCODE_2015.tsv"),
        "tf_gsea_gtrd": gsea(old_tf / "GTRD/Tables/GSEA_Full_Table_GTRD_TF.csv", new_tf / "GTRD/Tables/GSEA_Full_Table_GTRD_TF.tsv"),
        "custom_gsea": gsea(args.legacy / "06_Custom/Tables/gsea/GSEA_CustomGeneSets_FullTable.csv", args.new / f"07_Exploratory/Custom_Gene_Sets/Tables/{args.contrast}_GSEA_CustomGeneSets_FullTable.tsv"),
    }
    old_ppi = args.legacy / "05_Network/PPI/Tables/PPI_Hub_Ranking.csv"
    new_ppi = args.new / f"05_Network/{args.contrast}/PPI/Tables/PPI_Hub_Ranking.tsv"
    results["ppi_hubs"] = ppi(old_ppi, new_ppi) if old_ppi.is_file() else ppi_explicit_skip(
        args.legacy / "05_Network", args.new)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(results, indent=2), encoding="utf-8")
    print(json.dumps(results, indent=2))
    return 0 if all(item["passed"] for item in results.values()) else 2


if __name__ == "__main__":
    raise SystemExit(main())
