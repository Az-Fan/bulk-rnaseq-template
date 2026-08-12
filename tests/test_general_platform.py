from __future__ import annotations

import hashlib
import json
import os
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from internal.lib.core import module_confirmation_string, validate_project


ROOT = Path(__file__).resolve().parents[1]
RSCRIPT = ROOT / ".pixi/envs/default/bin/Rscript"


def test_public_orchestrator_does_not_bypass_modules():
    text = (ROOT / "workflow/run.R").read_text(encoding="utf-8")
    assert "internal/engines" not in text
    assert "system2(" not in text
    for function in (
        "run_qc", "run_differential", "run_enrichment", "run_activity",
        "run_network", "run_motif", "run_exploratory",
    ):
        assert function in text


def test_public_modules_are_free_of_project_and_two_group_hardcoding():
    text = "\n".join(path.read_text(encoding="utf-8") for path in (ROOT / "workflow").glob("*.R"))
    for forbidden in ("FKBP1A", "Treatment_vs_Control", "levels = c(\"Control\", \"Treatment\")", "hg38-motif"):
        assert forbidden not in text


def test_publication_plot_contract_includes_classic_redesigns():
    text = (ROOT / "workflow/functions.R").read_text(encoding="utf-8")
    for function in (
        "plot_pca_publication", "plot_volcano_publication", "plot_expression_boxplot",
        "plot_ora_publication", "plot_gsea_overview", "plot_gsea_curve",
    ):
        assert f"{function} <- function" in text
    assert "every tested term remains in the tables" in text
    assert "Running enrichment score" in text
    assert "Rank in complete DESeq2 Wald-statistic list" in text


def build_project(project: Path) -> Path:
    project.mkdir()
    (project / "input").mkdir()
    groups = ["Reference"] * 3 + ["DoseLow"] * 3 + ["DoseHigh"] * 3
    batches = ["B1", "B2", "B1", "B1", "B2", "B1", "B1", "B2", "B1"]
    samples = [f"S{i + 1}" for i in range(9)]
    rng = np.random.default_rng(104729)
    matrix = rng.poisson(80, size=(80, 9))
    matrix[:20, 3:6] += 40
    matrix[:20, 6:9] += 90
    counts = pd.DataFrame(matrix, columns=samples)
    counts.insert(0, "Name", [f"Gene{i + 1}" for i in range(80)])
    counts.insert(0, "Gene_ID", [f"ENSG{i + 1:011d}.1" for i in range(80)])
    counts_path = project / "input/counts.tsv"
    counts.to_csv(counts_path, sep="\t", index=False)
    digest = hashlib.sha256(counts_path.read_bytes()).hexdigest()

    pd.DataFrame({
        "sample": samples, "canonical_sample": samples, "group": groups, "batch": batches,
        "excluded": ["false"] * 9, "exclusion_reason": [""] * 9,
    }).to_csv(project / "samples.tsv", sep="\t", index=False)
    pd.DataFrame([
        {"contrast_id": "Low_vs_Reference", "factor": "group", "numerator": "DoseLow", "denominator": "Reference", "confirmed": "true"},
        {"contrast_id": "High_vs_Reference", "factor": "group", "numerator": "DoseHigh", "denominator": "Reference", "confirmed": "true"},
    ]).to_csv(project / "contrasts.tsv", sep="\t", index=False)
    (project / "input/source_manifest.yml").write_text(yaml.safe_dump({"schema_version": 1, "sha256": digest}), encoding="utf-8")
    (project / "qc_approval.yml").write_text(yaml.safe_dump({"schema_version": 1, "approved": True, "counts_sha256": digest}), encoding="utf-8")

    enabled = {"status": "enabled", "confirmed": True}
    not_applicable = lambda reason: {"status": "not_applicable", "reason": reason, "confirmed": True}
    config = {
        "schema_version": 4, "project_id": project.name, "species": "human", "assembly": "hg38",
        "counts_provenance": {
            "source_method": "featureCounts", "feature_level": "gene", "normalization": "raw_counts",
            "technical_replicates": "none", "upstream_filtering": "none", "strandedness": "unknown",
            "genome_annotation": "test-fixture", "confirmed": True,
        },
        "design": {"formula": "~ batch + group", "factor_levels": {"batch": ["B1", "B2"], "group": ["Reference", "DoseLow", "DoseHigh"]}, "display_factor": "group", "confirmed": True},
        "inputs": {"counts_file": "input/counts.tsv", "samples_file": "samples.tsv", "contrasts_file": "contrasts.tsv", "qc_approval_file": "qc_approval.yml", "gene_id_column": "Gene_ID", "gene_symbol_column": "Name", "strip_gene_version": True},
        "upstream": {
            "status": "disabled", "confirmed": True, "provider": "xuzhougeng/omics-pipelines",
            "repository": "https://github.com/xuzhougeng/omics-pipelines.git",
            "commit": "ce4e2ec88da6663a32b7099c5850e1a51ad66952", "workflow": "rnaseq",
            "execution_host": "server", "allow_wsl": False, "target": "results/04-quant/counts.tsv",
        },
        "filtering": {"min_count": 10, "min_samples": 3},
        "thresholds": {"mode": "screening", "padj": 0.05, "abs_lfc": 1.0},
        "analysis": {"random_seed": 104729, "modules": {
            "qc": enabled, "differential": enabled,
            "enrichment": not_applicable("No frozen gene-set fixture"),
            "regulation": not_applicable("No frozen activity fixture"),
            "network": not_applicable("No frozen network fixture"),
            "motif_promoter": not_applicable("No genome fixture"),
            "motif_peaks": not_applicable("No peaks fixture"),
            "custom_gene_sets": not_applicable("No custom sets fixture"),
            "pathview": not_applicable("No pathway map fixture"),
            "personalized": not_applicable("No personalized module fixture"),
            "wgcna": not_applicable("Fewer than 20 samples"),
        }},
        "report": {"detail": "standard"}, "export": {"figures": "publication", "formats": ["pdf"]},
        "resources": {"frozen_only": True}, "compute": {"max_workers": 1, "threads_per_worker": 1},
    }
    (project / "project.yml").write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
    return project


def test_generic_batch_design_multiple_contrasts_runs_to_complete(tmp_path):
    project = build_project(tmp_path / "Generic-Batch-Test")
    validate_project(project / "project.yml")
    run = project / "work/staging/test-run"
    cfg = validate_project(project / "project.yml")
    env = os.environ.copy()
    env["BULK_RNASEQ_MODULE_CONFIRMATION"] = module_confirmation_string(cfg)
    completed = subprocess.run(
        [str(RSCRIPT), str(ROOT / "workflow/run.R"), str(project), str(run)],
        cwd=ROOT, env=env, text=True, capture_output=True, timeout=180,
    )
    assert completed.returncode == 0, completed.stdout + "\n" + completed.stderr
    manifest = json.loads((run / "run_manifest.json").read_text(encoding="utf-8"))
    assert manifest["status"] == "complete"
    assert manifest["design_formula"] == "~ batch + group"
    assert manifest["contrasts"] == ["Low_vs_Reference", "High_vs_Reference"]
    for contrast in manifest["contrasts"]:
        table = run / f"Comparisons/{contrast}/02_Differential/Tables/Differential_Full.tsv"
        assert table.is_file()
        assert len(pd.read_csv(table, sep="\t")) == 80
    assert (run / "01_QC/Figures/Final_PCA.pdf").is_file()
    assert (run / "Comparisons/High_vs_Reference/02_Differential/Figures/Volcano.pdf").is_file()
    assert not list(run.rglob("*.png"))
    assert not list(run.rglob("*.svg"))
    assert all(item["status"] in {"complete", "not_applicable", "skipped_by_user"} for item in manifest["modules"].values())


def test_publisher_archives_previous_results_and_promotes_one_complete_run(tmp_path):
    project = tmp_path / "Publish-Test"
    staging = project / "work/staging/run-2"
    old_results = project / "results"
    for module in (
        "01_QC", "02_Differential", "03_Enrichment", "04_Regulation",
        "05_Network", "06_Motif", "07_Exploratory",
    ):
        (staging / module).mkdir(parents=True, exist_ok=True)
        (staging / module / "index.html").write_text(module, encoding="utf-8")
    (staging / "Provenance").mkdir()
    (staging / "index.html").write_text("new", encoding="utf-8")
    (staging / "run_manifest.json").write_text(json.dumps({
        "status": "complete", "project_id": "Publish-Test", "run_id": "run-2",
    }), encoding="utf-8")
    (project / "project.yml").write_text(yaml.safe_dump({"migration": {"requires_legacy_parity": False}}), encoding="utf-8")
    old_results.mkdir(parents=True)
    (old_results / "index.html").write_text("old", encoding="utf-8")
    completed = subprocess.run([
        str(ROOT / ".pixi/envs/default/bin/python"),
        str(ROOT / "internal/reporting/publish_project_results.py"),
        str(staging), str(project), "--confirm", "PUBLISH:Publish-Test:run-2",
    ], cwd=ROOT, text=True, capture_output=True)
    assert completed.returncode == 0, completed.stderr
    assert (project / "results/index.html").read_text(encoding="utf-8") == "new"
    archived = list((project / "work/history/published").glob("*/index.html"))
    assert len(archived) == 1
    assert archived[0].read_text(encoding="utf-8") == "old"
    assert not staging.exists()
    assert (project / "results/Provenance/publication_manifest.tsv").is_file()


def test_publisher_blocks_migrated_project_when_legacy_parity_fails(tmp_path):
    project = tmp_path / "Legacy-Project"
    staging = project / "work/staging/run-3"
    for module in (
        "01_QC", "02_Differential", "03_Enrichment", "04_Regulation",
        "05_Network", "06_Motif", "07_Exploratory",
    ):
        (staging / module).mkdir(parents=True, exist_ok=True)
        (staging / module / "index.html").write_text(module, encoding="utf-8")
    (staging / "run_manifest.json").write_text(json.dumps({
        "status": "complete", "project_id": "Legacy-Project", "run_id": "run-3",
    }), encoding="utf-8")
    (project / "project.yml").write_text(yaml.safe_dump({"migration": {"requires_legacy_parity": True}}), encoding="utf-8")
    audit = project / "work/audits/run-3"
    audit.mkdir(parents=True)
    (audit / "legacy_inventory_summary.json").write_text(json.dumps({"passed": False, "missing_families": 10}), encoding="utf-8")
    completed = subprocess.run([
        str(ROOT / ".pixi/envs/default/bin/python"),
        str(ROOT / "internal/reporting/publish_project_results.py"),
        str(staging), str(project), "--confirm", "PUBLISH:Legacy-Project:run-3",
    ], cwd=ROOT, text=True, capture_output=True)
    assert completed.returncode != 0
    assert "legacy output parity failed" in completed.stderr
    assert staging.exists()
    assert not (project / "results").exists()
