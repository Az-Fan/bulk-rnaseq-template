from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_ora_and_gsea_are_first_class_enrichment_outputs():
    module = (ROOT / "workflow/03_enrichment.R").read_text(encoding="utf-8")
    figures = (ROOT / "workflow/functions.R").read_text(encoding="utf-8")
    assert "plot_ora_publication" in module
    assert "plot_gsea_overview" in module
    assert "plot_gsea_curve" in module
    assert "universe <- unique" in module
    assert 'cluster <- list(Up = up, Down = down, All = unique(c(up, down)))' in module
    assert 'clusterProfiler::compareCluster' in module
    assert 'clusterProfiler::GSEA' in module
    assert 'pvalueCutoff = 1' in module
    assert 'BiocParallel::SerialParam()' in module
    assert "fold_enrichment" in figures
    assert "Running enrichment score" in figures
    assert "complete DESeq2 Wald-statistic list" in figures


def test_enrichment_plot_contract_preserves_direction_and_full_tables():
    figures = (ROOT / "workflow/functions.R").read_text(encoding="utf-8")
    assert "direction_palette" in figures
    assert "complete tested terms remain" in figures
    assert "GeneRatio" in figures or "overlap_count" in figures


def test_ora_databases_and_legacy_figure_families_are_explicit():
    module = (ROOT / "workflow/03_enrichment.R").read_text(encoding="utf-8")
    for database in ("GO_BP", "GO_CC", "GO_MF", "KEGG", "Reactome", "Hallmark"):
        assert database in module
    for family in ("_Dotplot_", "_Lollipop_", "_Diverging_Bar", "_Cnet", "_Emap", "_Heat"):
        assert family in module
    assert 'file.path(ora_root, "Figures", "Sankey"' in module
    assert "Each database gets" in module


def test_legacy_audit_does_not_substitute_apear_for_cnet_or_emap():
    audit = (ROOT / "internal/lib/core.py").read_text(encoding="utf-8")
    assert '"cnet", r"Cnet", ("Cnet",)' in audit
    assert '"emap", r"Emap", ("Emap",)' in audit
    assert '"Cnet", "aPEAR"' not in audit
    assert '"Emap", "aPEAR"' not in audit
    assert "_legacy_specific_candidates" in audit
