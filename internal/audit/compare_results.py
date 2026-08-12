#!/usr/bin/env python3
from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd


def read_auto(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t" if path.suffix in {".tsv", ".txt"} else ",")


def main() -> None:
    old = Path(sys.argv[1]); run = Path(sys.argv[2]); out = Path(sys.argv[3]); out.mkdir(parents=True, exist_ok=True)
    new = run / "Profiles/Primary_padj0.05_LFC1.0"
    old_h = read_auto(old / "03_Enrichment/GSEA/Tables/GSEA_Full_Table_Hallmark.csv")
    new_h = read_auto(new / "03_Enrichment/GSEA/Hallmark/Tables/GSEA_Full_Table_Hallmark.tsv")
    old_path = next(c for c in ["pathway", "ID", "Description"] if c in old_h)
    old_nes = next(c for c in ["NES", "nes"] if c in old_h)
    new_path = next(c for c in ["pathway", "Description", "ID"] if c in new_h)
    h = old_h[[old_path, old_nes]].rename(columns={old_path: "pathway", old_nes: "NES_old"}).merge(
        new_h[[new_path, "NES"]].rename(columns={new_path: "pathway", "NES": "NES_v4"}), on="pathway"
    )
    h["abs_delta"] = (h.NES_v4 - h.NES_old).abs()
    h.to_csv(out / "hallmark_nes_comparison.tsv", sep="\t", index=False)

    old_p = read_auto(old / "04_Regulation/PROGENy/Tables/PROGENy_Contrast_Activity.csv")
    new_p = read_auto(run / "04_Regulation/PROGENy/PROGENy_Contrast_Activity.tsv")
    old_source = next(c for c in ["source", "pathway", "condition"] if c in old_p)
    old_score = next(c for c in ["score", "activity", "statistic"] if c in old_p)
    p = old_p[[old_source, old_score]].rename(columns={old_source: "pathway", old_score: "score_old"}).merge(
        new_p[["source", "score"]].rename(columns={"source": "pathway", "score": "score_v4"}), on="pathway"
    )
    p["direction_agrees"] = (p.score_old * p.score_v4) > 0
    p.to_csv(out / "progeny_comparison.tsv", sep="\t", index=False)

    old_hubs = read_auto(old / "05_Network/PPI/Tables/PPI_Hub_Ranking.csv")
    new_hubs = read_auto(new / "05_Network/PPI/PPI_Hub_Genes.tsv")
    old_gene = next(c for c in ["gene_symbol", "Symbol", "name"] if c in old_hubs)
    old_degree = next(c for c in ["degree", "Degree"] if c in old_hubs)
    top_old = set(old_hubs.sort_values(old_degree, ascending=False).head(20)[old_gene].astype(str))
    top_new = set(new_hubs.sort_values("degree", ascending=False).head(20).gene_symbol.astype(str))
    pd.DataFrame({"metric": ["top20_intersection", "top20_jaccard"], "value": [len(top_old & top_new), len(top_old & top_new) / len(top_old | top_new)]}).to_csv(out / "ppi_hub_comparison.tsv", sep="\t", index=False)

    summary = {
        "hallmark_shared": len(h), "hallmark_nes_spearman": float(h.NES_old.corr(h.NES_v4, method="spearman")),
        "progeny_shared": len(p), "progeny_direction_agreement": float(p.direction_agrees.mean()) if len(p) else None,
        "ppi_top20_intersection": len(top_old & top_new), "ppi_top20_jaccard": len(top_old & top_new) / len(top_old | top_new),
        "comparability_notes": [
            "Hallmark uses the same frozen MSigDB 2025.1 release; score differences reflect implementation/runtime details.",
            "PROGENy v4 uses a transparent weighted z-score on the frozen network; scale is not expected to equal the legacy decoupler score.",
            "PPI hubs are exploratory and sensitive to graph selection; neither version establishes causality.",
        ],
    }
    (out / "extended_comparison_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")


if __name__ == "__main__": main()
