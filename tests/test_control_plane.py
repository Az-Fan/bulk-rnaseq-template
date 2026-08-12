from pathlib import Path

import pytest
import yaml

from internal.lib.core import (
    clean_preview,
    import_counts,
    inventory_legacy_outputs,
    module_confirmation_string,
    validate_project,
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


def test_unconfirmed_upstream_decision_is_rejected(tmp_path):
    cfg = confirmed_config()
    cfg["upstream"].update({"status": "unconfirmed", "confirmed": False})
    path = tmp_path / "project.yml"
    path.write_text(yaml.safe_dump(cfg), encoding="utf-8")
    with pytest.raises(ValueError, match="Upstream-data"):
        validate_project(path)


def test_frozen_resource_registry_is_valid():
    resources = validate_resources(ROOT / "resources/registry.yml")
    assert resources
    assert all(len(resource["sha256"]) == 64 for resource in resources)


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
