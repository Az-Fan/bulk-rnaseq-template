from pathlib import Path
import subprocess
import sys

import pytest
import yaml

from internal.lib.core import (
    analysis_decision_rows,
    clean_preview,
    import_counts,
    inventory_legacy_outputs,
    module_confirmation_string,
    plan_confirmation_token,
    validate_project,
    validate_qc_preview,
    validate_resources,
)

ROOT = Path(__file__).resolve().parents[1]
TEMPLATE = ROOT / "projects/_template/project.yml"


def confirmed_config() -> dict:
    cfg = yaml.safe_load(TEMPLATE.read_text(encoding="utf-8"))
    cfg["project_id"] = "Generic-Test"
    cfg["counts_provenance"].update({"normalization": "raw_counts", "feature_level": "gene", "confirmed": True})
    cfg["design"]["confirmed"] = True
    cfg["upstream"].update({"status": "disabled", "confirmed": True})
    for group in cfg["decision_confirmation"]:
        cfg["decision_confirmation"][group] = True
    for decision in cfg["analysis"]["modules"].values():
        decision["confirmed"] = True
    return cfg


def test_blank_template_is_deliberately_unconfirmed():
    with pytest.raises(ValueError, match="unconfirmed|Unconfirmed"):
        validate_project(TEMPLATE)


def test_confirmed_generic_project_is_valid_and_pdf_only(tmp_path):
    path = tmp_path / "project.yml"
    path.write_text(yaml.safe_dump(confirmed_config(), sort_keys=False), encoding="utf-8")
    cfg = validate_project(path)
    assert cfg["export"]["formats"] == ["pdf"]
    assert cfg["upstream"]["status"] == "disabled"


def test_normalized_input_is_rejected(tmp_path):
    cfg = confirmed_config()
    cfg["counts_provenance"]["normalization"] = "TPM"
    path = tmp_path / "project.yml"
    path.write_text(yaml.safe_dump(cfg), encoding="utf-8")
    with pytest.raises(ValueError, match="forbidden"):
        validate_project(path)


def test_unconfirmed_module_is_rejected(tmp_path):
    cfg = confirmed_config()
    cfg["analysis"]["modules"]["network"]["confirmed"] = False
    path = tmp_path / "project.yml"
    path.write_text(yaml.safe_dump(cfg), encoding="utf-8")
    with pytest.raises(ValueError, match="Unconfirmed modules"):
        validate_project(path)


def test_module_confirmation_is_complete_and_deterministic():
    cfg = confirmed_config()
    token = module_confirmation_string(cfg)
    assert token.startswith("custom_gene_sets=skipped_by_user,differential=enabled")
    assert token.count(",") == len(cfg["analysis"]["modules"]) - 1


def test_run_plan_covers_fine_grained_decisions_and_changes_with_config(tmp_path):
    project = tmp_path / "project"
    (project / "input").mkdir(parents=True)
    cfg = confirmed_config()
    (project / "project.yml").write_text(yaml.safe_dump(cfg, sort_keys=False), encoding="utf-8")
    (project / "samples.tsv").write_text("sample\tgroup\nS1\tREFERENCE_LEVEL\nS2\tCOMPARISON_LEVEL\n", encoding="utf-8")
    (project / "contrasts.tsv").write_text(
        "contrast_id\tfactor\tnumerator\tdenominator\tconfirmed\nC\tgroup\tCOMPARISON_LEVEL\tREFERENCE_LEVEL\ttrue\n",
        encoding="utf-8",
    )
    (project / "qc_approval.yml").write_text("approved: false\n", encoding="utf-8")
    (project / "input/source_manifest.yml").write_text("sha256: pending\n", encoding="utf-8")
    rows = analysis_decision_rows(project)
    assert {row["decision"] for row in rows} >= {
        "counts_provenance", "contrasts", "gene_filter", "deg_definition",
        "ora_and_gsea", "tf_gsva_progeny", "ppi", "promoter_and_peak_motif",
        "target_genes_custom_pathview", "report_and_export",
    }
    first = plan_confirmation_token(project)
    cfg["thresholds"]["padj"] = 0.01
    (project / "project.yml").write_text(yaml.safe_dump(cfg, sort_keys=False), encoding="utf-8")
    assert plan_confirmation_token(project) != first


def test_unimplemented_module_is_rejected_instead_of_silently_completed(tmp_path):
    cfg = confirmed_config()
    cfg["analysis"]["modules"]["wgcna"] = {"status": "enabled", "confirmed": True}
    path = tmp_path / "project.yml"
    path.write_text(yaml.safe_dump(cfg), encoding="utf-8")
    with pytest.raises(ValueError, match="WGCNA execution is not implemented"):
        validate_project(path)


def test_analyze_stops_without_fresh_plan_token(tmp_path):
    project = tmp_path / "project"
    project.mkdir()
    cfg = confirmed_config()
    (project / "project.yml").write_text(yaml.safe_dump(cfg, sort_keys=False), encoding="utf-8")
    completed = subprocess.run(
        [sys.executable, str(ROOT / "internal/cli/pipeline.py"), "analyze", "--project", str(project),
         "--confirm-modules", module_confirmation_string(cfg)],
        cwd=ROOT, capture_output=True, text=True,
    )
    assert completed.returncode == 2
    assert "stopped_unconfirmed_plan" in completed.stdout


def test_qc_preview_requires_only_pre_qc_decisions(tmp_path):
    cfg = confirmed_config()
    for group in cfg["decision_confirmation"]:
        cfg["decision_confirmation"][group] = group == "input_design_filtering"
    for name, decision in cfg["analysis"]["modules"].items():
        if name != "qc":
            decision["confirmed"] = False
    path = tmp_path / "project.yml"
    path.write_text(yaml.safe_dump(cfg, sort_keys=False), encoding="utf-8")
    assert validate_qc_preview(path)["analysis"]["modules"]["qc"]["status"] == "enabled"


def test_unconfirmed_upstream_decision_is_rejected(tmp_path):
    cfg = confirmed_config()
    cfg["upstream"].update({"status": "unconfirmed", "confirmed": False})
    path = tmp_path / "project.yml"
    path.write_text(yaml.safe_dump(cfg), encoding="utf-8")
    with pytest.raises(ValueError, match="Upstream-data"):
        validate_project(path)


def test_frozen_resource_registry_metadata_is_valid_in_clean_clone():
    resources = validate_resources(ROOT / "resources/registry.yml", require_local_files=False)
    assert resources
    assert all(len(resource["sha256"]) == 64 for resource in resources)


def test_runtime_resource_validation_stays_strict(tmp_path):
    registry = tmp_path / "registry.yml"
    registry.write_text(
        yaml.safe_dump({"schema_version": 1, "resources": [{
            "resource_id": "missing_resource",
            "species": "Homo sapiens",
            "assembly": "GRCh38",
            "release": "test",
            "source_url": "https://example.org/resource.tsv",
            "downloaded_at": "2026-08-13T00:00:00+00:00",
            "sha256": "0" * 64,
            "license_or_access_note": "test only",
            "local_path": "data/resource.tsv",
        }]}),
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="local file is missing"):
        validate_resources(registry)


def test_missing_output_fails_inventory(tmp_path):
    old = tmp_path / "old"
    new = tmp_path / "new"
    (old / "02_Differential/Tables").mkdir(parents=True)
    new.mkdir()
    (old / "02_Differential/Tables/DEG_results_annotated.csv").write_text("gene_id\nA\n")
    _, _, summary = inventory_legacy_outputs(old, new, tmp_path / "audit")
    assert not summary["passed"]
    assert summary["missing_families"] == 1


def test_counts_import_rejects_non_integer_values(tmp_path):
    project = tmp_path / "project"
    project.mkdir()
    cfg = confirmed_config()
    (project / "project.yml").write_text(yaml.safe_dump(cfg), encoding="utf-8")
    (project / "samples.tsv").write_text("sample\tgroup\nS1\tREFERENCE_LEVEL\nS2\tCOMPARISON_LEVEL\n", encoding="utf-8")
    source = tmp_path / "counts.tsv"
    source.write_text("gene_id\tS1\tS2\nG1\t1.5\t2\n", encoding="utf-8")
    with pytest.raises(ValueError, match="non-integer"):
        import_counts(source, project)


def test_clean_is_scoped_and_previewable(tmp_path):
    target = tmp_path / "projects" / "demo" / "work" / "staging"
    target.mkdir(parents=True)
    (target / "result.tsv").write_text("x\n1\n", encoding="utf-8")
    preview = clean_preview("demo", ["staging"], root=tmp_path)
    assert preview[0]["files"] == 1
    assert (target / "result.tsv").exists()
