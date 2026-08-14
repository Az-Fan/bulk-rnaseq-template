#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


def read(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t" if path.suffix == ".tsv" else ",")


def compare_table(legacy: Path, new: Path, kind: str, database: str) -> dict:
    old, current = read(legacy), read(new)
    key = "ID"
    shared = sorted(set(old[key].astype(str)) & set(current[key].astype(str)))
    old_ids, new_ids = set(old[key].astype(str)), set(current[key].astype(str))
    out = {
        "kind": kind,
        "database": database,
        "legacy_rows": len(old),
        "new_rows": len(current),
        "shared_ids": len(shared),
        "id_set_equal": old_ids == new_ids,
    }
    tolerances = {"enrichmentScore": 5e-7, "NES": 5e-7}
    if kind == "ORA":
        key_cols = ["ID", "Cluster"]
        old_keys = set(map(tuple, old[key_cols].astype(str).to_numpy()))
        new_keys = set(map(tuple, current[key_cols].astype(str).to_numpy()))
        out["id_direction_set_equal"] = old_keys == new_keys
        merged = old.merge(current, on=key_cols, suffixes=("_old", "_new"))
        columns = ["pvalue", "p.adjust", "qvalue", "Count"]
    else:
        merged = old.merge(current, on="ID", suffixes=("_old", "_new"))
        columns = ["enrichmentScore", "NES", "pvalue", "p.adjust", "qvalue", "setSize"]
        out["direction_equal"] = bool(
            np.array_equal(np.sign(merged["NES_old"]), np.sign(merged["NES_new"]))
        ) if len(merged) else False
    for column in columns:
        left, right = merged[f"{column}_old"], merged[f"{column}_new"]
        numeric = pd.to_numeric(left, errors="coerce").notna() & pd.to_numeric(right, errors="coerce").notna()
        if numeric.any():
            a = pd.to_numeric(left[numeric]); b = pd.to_numeric(right[numeric])
            out[f"{column}_max_abs_diff"] = float(np.max(np.abs(a - b)))
            atol = tolerances.get(column, 1e-10)
            out[f"{column}_tolerance"] = atol
            out[f"{column}_allclose"] = bool(np.allclose(a, b, rtol=1e-8, atol=atol, equal_nan=True))
        else:
            out[f"{column}_equal"] = bool(left.astype(str).equals(right.astype(str)))
    checks = [out["id_set_equal"]]
    checks.append(out.get("id_direction_set_equal", out.get("direction_equal", False)))
    checks.extend(value for key, value in out.items() if key.endswith("_allclose"))
    out["passed"] = bool(all(checks))
    return out


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--legacy", type=Path, required=True)
    parser.add_argument("--new", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    rows = []
    for database in ["GO_BP", "GO_CC", "GO_MF"]:
        rows.append(compare_table(
            args.legacy / "ORA/GO/Tables" / f"Enrichment_Full_{database}.csv",
            args.new / "ORA/Tables" / f"Enrichment_Full_{database}_Primary.tsv",
            "ORA", database,
        ))
    for database in ["KEGG", "Hallmark"]:
        rows.append(compare_table(
            args.legacy / f"ORA/{database}/Tables" / f"Enrichment_Full_{database}.csv",
            args.new / "ORA/Tables" / f"Enrichment_Full_{database}_Primary.tsv",
            "ORA", database,
        ))
    for database in ["Reactome", "Hallmark", "KEGG", "GO_BP"]:
        rows.append(compare_table(
            args.legacy / "GSEA/Tables" / f"GSEA_Full_Table_{database}.csv",
            args.new / "GSEA/Tables" / f"GSEA_Full_Table_{database}.tsv",
            "GSEA", database,
        ))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(args.output, sep="\t", index=False)
    print(json.dumps(rows, indent=2))
    return 0 if all(row["passed"] for row in rows) else 2


if __name__ == "__main__":
    raise SystemExit(main())
