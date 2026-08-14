from pathlib import Path
import subprocess


ROOT = Path(__file__).resolve().parents[1]


def test_top_level_is_the_locked_multi_project_template():
    expected_dirs = {"workflow", "projects", "resources", "tests", "internal"}
    actual_dirs = {path.name for path in ROOT.iterdir() if path.is_dir() and not path.name.startswith(".")}
    assert actual_dirs == expected_dirs


def test_public_workflow_has_exactly_seven_scientific_modules():
    modules = sorted(path.name for path in (ROOT / "workflow").glob("[0-9][0-9]_*.R"))
    assert modules == [
        "01_qc.R", "02_differential.R", "03_enrichment.R", "04_activity.R",
        "05_network.R", "06_motif.R", "07_exploratory.R",
    ]
    assert (ROOT / "workflow/run.R").is_file()
    assert (ROOT / "workflow/functions.R").is_file()


def test_only_blank_project_template_is_committed():
    tracked = subprocess.run(
        ["git", "ls-files", "projects"], cwd=ROOT, check=True,
        capture_output=True, text=True,
    ).stdout.splitlines()
    assert tracked
    assert all(path.startswith("projects/_template/") for path in tracked)
    template = ROOT / "projects/_template"
    assert (template / "project.yml").is_file()
    assert (template / "input/README.md").is_file()
    assert (template / "results/README.md").is_file()
    assert not list(template.rglob("*.pdf"))
    assert not list(template.rglob("*.png"))


def test_upstream_is_pinned_but_not_embedded_in_public_workflow():
    adapter = (ROOT / "internal/upstream/adapter.py").read_text(encoding="utf-8")
    assert "ce4e2ec88da6663a32b7099c5850e1a51ad66952" in adapter
    assert 'TARGET = "results/04-quant/counts.tsv"' in adapter
    workflow = (ROOT / "workflow/run.R").read_text(encoding="utf-8")
    assert "omics-pipelines" not in workflow
