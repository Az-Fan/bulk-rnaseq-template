from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_ora_and_gsea_are_first_class_enrichment_outputs():
    module = (ROOT / "workflow/03_enrichment.R").read_text(encoding="utf-8")
    figures = (ROOT / "workflow/functions.R").read_text(encoding="utf-8")
    assert "plot_ora_publication" in module
    assert "plot_gsea_overview" in module
    assert "plot_gsea_curve" in module
    assert "universe <- unique" in module
    assert 'for (direction in c("Up", "Down"))' in module
    assert "fold_enrichment" in figures
    assert "Running enrichment score" in figures
    assert "complete DESeq2 Wald-statistic list" in figures


def test_enrichment_plot_contract_preserves_direction_and_full_tables():
    figures = (ROOT / "workflow/functions.R").read_text(encoding="utf-8")
    assert "direction_palette" in figures
    assert "every tested term remains in the tables" in figures
    assert "GeneRatio" in figures or "overlap_count" in figures
