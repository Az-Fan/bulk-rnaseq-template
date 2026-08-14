#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


def table(path: Path, sep: str, index_col=None) -> pd.DataFrame:
    return pd.read_csv(path, sep=sep, index_col=index_col)


def long_compare(old: pd.DataFrame, new: pd.DataFrame, keys: list[str], values: list[str], atol: float, contrast: str) -> dict:
    if "condition" in keys and set(old["condition"].astype(str)) == {"stat"}:
        old = old.copy(); old["condition"] = contrast
    merged = old.merge(new, on=keys, suffixes=("_old", "_new"))
    old_keys = set(map(tuple, old[keys].astype(str).to_numpy()))
    new_keys = set(map(tuple, new[keys].astype(str).to_numpy()))
    result = {"old_rows": len(old), "new_rows": len(new), "shared_rows": len(merged),
              "key_set_equal": old_keys == new_keys, "absolute_tolerance": atol}
    for value in values:
        a, b = merged[f"{value}_old"], merged[f"{value}_new"]
        result[f"{value}_max_abs_diff"] = float(np.nanmax(np.abs(a - b)))
        result[f"{value}_allclose"] = bool(np.allclose(a, b, atol=atol, rtol=1e-8, equal_nan=True))
    result["significance_class_equal"] = bool(np.array_equal(
        pd.to_numeric(merged["p_value_old"], errors="coerce").lt(0.05),
        pd.to_numeric(merged["p_value_new"], errors="coerce").lt(0.05),
    )) if "p_value" in values else True
    result["passed"] = result["key_set_equal"] and result["significance_class_equal"] and all(
        value for key, value in result.items() if key.endswith("_allclose"))
    return result


def matrix_compare(old: pd.DataFrame, new: pd.DataFrame, atol: float) -> dict:
    rows = sorted(set(old.index) & set(new.index)); cols = sorted(set(old.columns) & set(new.columns))
    a, b = old.loc[rows, cols].to_numpy(float), new.loc[rows, cols].to_numpy(float)
    result = {
        "old_shape": list(old.shape), "new_shape": list(new.shape),
        "row_set_equal": set(old.index) == set(new.index), "column_set_equal": set(old.columns) == set(new.columns),
        "max_abs_diff": float(np.nanmax(np.abs(a - b))), "absolute_tolerance": atol,
        "allclose": bool(np.allclose(a, b, atol=atol, rtol=1e-8, equal_nan=True)),
        "pearson": float(np.corrcoef(a.ravel(), b.ravel())[0, 1]),
        "sign_equal": bool(np.array_equal(np.sign(a), np.sign(b))),
    }
    n_top = min(20, len(rows))
    old_top = set(old.loc[rows, cols].std(axis=1).nlargest(n_top).index)
    new_top = set(new.loc[rows, cols].std(axis=1).nlargest(n_top).index)
    result["top_variable_set_equal"] = old_top == new_top
    result["passed"] = all((result["row_set_equal"], result["column_set_equal"], result["allclose"],
                            result["sign_equal"], result["top_variable_set_equal"]))
    return result


def main() -> int:
    ap = argparse.ArgumentParser(); ap.add_argument("--legacy", type=Path, required=True); ap.add_argument("--new", type=Path, required=True); ap.add_argument("--output", type=Path, required=True); ap.add_argument("--contrast", required=True)
    args = ap.parse_args(); out = {}
    new_tables = args.new / "Tables" if (args.new / "Tables").is_dir() else args.new / "TF_Activity/Tables"
    progeny_tables = args.new / "Tables" if (args.new / "Tables").is_dir() else args.new / "PROGENy/Tables"
    gsva_tables = args.new / "Tables" if (args.new / "Tables").is_dir() else args.new / "GSVA/Tables"
    out["tf_contrast"] = long_compare(
        table(args.legacy / "TF_Activity/Tables/TF_Activity_Contrast_Full.csv", ","),
        table(new_tables / f"TF_Activity_Contrast_{args.contrast}.tsv", "\t"),
        ["source", "condition"], ["score", "p_value"], 1e-6, args.contrast)
    out["progeny_contrast"] = long_compare(
        table(args.legacy / "PROGENy/Tables/PROGENy_Contrast_Activity.csv", ","),
        table(progeny_tables / f"PROGENy_Contrast_Activity_{args.contrast}.tsv", "\t"),
        ["source", "condition"], ["score", "p_value"], 1e-7, args.contrast)
    old_gsva = table(args.legacy / "GSVA/Tables/GSVA_Score_Matrix.csv", ",", 0)
    new_gsva = table(gsva_tables / "GSVA_Score_Matrix.tsv", "\t", 0)
    # The historical Windows run used GSVA 2.0.4; the locked Linux platform uses
    # GSVA 2.4.4. Cross-version regression established a worst-case score delta
    # below 0.004 with identical signs and top-variable pathway membership.
    out["gsva"] = matrix_compare(old_gsva, new_gsva, 4e-3)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps(out, indent=2))
    return 0 if all(result["passed"] for result in out.values()) else 2


if __name__ == "__main__": raise SystemExit(main())
