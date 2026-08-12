from __future__ import annotations

import json
import os
import platform
import shutil
import subprocess
from pathlib import Path

import yaml

from internal.lib.core import ROOT, load_yaml, sha256_file


PROVIDER = "xuzhougeng/omics-pipelines"
REPOSITORY = "https://github.com/xuzhougeng/omics-pipelines.git"
PINNED_COMMIT = "ce4e2ec88da6663a32b7099c5850e1a51ad66952"
TARGET = "results/04-quant/counts.tsv"


def _inside(project: Path, value: str) -> Path:
    path = (project / value).resolve()
    if path != project and project not in path.parents:
        raise ValueError(f"Upstream path escapes the project directory: {value}")
    return path


def _is_wsl() -> bool:
    return "microsoft" in platform.release().lower() or "WSL_DISTRO_NAME" in os.environ


def _configuration(project: Path) -> tuple[dict, dict]:
    cfg = load_yaml(project / "project.yml")
    upstream = cfg.get("upstream")
    if not isinstance(upstream, dict):
        raise ValueError("project.yml must contain an explicit upstream decision")
    return cfg, upstream


def plan_upstream(project: Path) -> dict:
    project = project.resolve()
    cfg, upstream = _configuration(project)
    token = f"UPSTREAM:{cfg['project_id']}:{PINNED_COMMIT}"
    return {
        "project_id": cfg["project_id"],
        "status": upstream.get("status", "unconfirmed"),
        "confirmed": bool(upstream.get("confirmed", False)),
        "provider": PROVIDER,
        "repository": REPOSITORY,
        "pinned_commit": PINNED_COMMIT,
        "execution_host": upstream.get("execution_host", "server"),
        "wsl_execution_allowed": bool(upstream.get("allow_wsl", False)),
        "target": TARGET,
        "scope": "fastp + FastQC + STAR + merged raw gene counts only",
        "excluded": ["upstream DESeq2", "upstream enrichment", "system R", "index construction"],
        "confirmation_token": token,
    }


def _validate_enabled(cfg: dict, upstream: dict) -> None:
    expected = {
        "status": "enabled", "confirmed": True, "provider": PROVIDER,
        "repository": REPOSITORY, "commit": PINNED_COMMIT,
        "workflow": "rnaseq", "target": TARGET,
    }
    wrong = [key for key, value in expected.items() if upstream.get(key) != value]
    if wrong:
        raise ValueError("Upstream execution is not explicitly enabled or pinned: " + ", ".join(wrong))
    if not upstream.get("samples_file") or not upstream.get("config_file"):
        raise ValueError("Enabled upstream workflow requires samples_file and config_file")


def _verify_source(source: Path) -> None:
    if not (source / "rnaseq" / "Snakefile").is_file():
        raise ValueError("Pinned upstream submodule is absent; clone with --recurse-submodules")
    commit = subprocess.run(
        ["git", "-C", str(source), "rev-parse", "HEAD"],
        check=True, capture_output=True, text=True,
    ).stdout.strip()
    if commit != PINNED_COMMIT:
        raise ValueError(f"Upstream source commit mismatch: {commit}")


def run_upstream(project: Path, confirmation: str, workers: int, dry_run: bool = False) -> dict:
    """Run only the fixed raw-count target on an explicitly approved Linux server."""
    project = project.resolve()
    cfg, upstream = _configuration(project)
    _validate_enabled(cfg, upstream)
    expected_token = f"UPSTREAM:{cfg['project_id']}:{PINNED_COMMIT}"
    if confirmation != expected_token:
        raise ValueError("Upstream confirmation token does not match this project and pinned commit")
    if workers < 1:
        raise ValueError("workers must be positive")
    if _is_wsl() and not upstream.get("allow_wsl", False):
        raise ValueError("This project forbids upstream execution in WSL; use the declared Linux server")

    source = ROOT / "internal/upstream/omics-pipelines"
    _verify_source(source)
    samples = _inside(project, upstream["samples_file"])
    config = _inside(project, upstream["config_file"])
    for path in (samples, config):
        if not path.is_file():
            raise ValueError(f"Missing upstream configuration input: {path}")

    source_rnaseq = source / "rnaseq"
    work = project / "work/upstream" / PINNED_COMMIT[:12]
    work.mkdir(parents=True, exist_ok=True)
    for name in ("Snakefile", "pixi.toml", "pixi.lock"):
        shutil.copy2(source_rnaseq / name, work / name)
    scripts = work / "scripts"
    scripts.mkdir(exist_ok=True)
    # The fixed count target needs only this script. Do not copy or expose the
    # upstream DESeq2/enrichment scripts in the project runtime.
    for old in scripts.iterdir():
        if old.is_file():
            old.unlink()
    shutil.copy2(source_rnaseq / "scripts/merge_star_counts.py", scripts / "merge_star_counts.py")
    shutil.copy2(samples, work / "samples.tsv")
    shutil.copy2(config, work / "config.yaml")

    upstream_cfg = yaml.safe_load((work / "config.yaml").read_text(encoding="utf-8"))
    index = Path(upstream_cfg["ref"]["star_index"]).expanduser()
    if not index.is_absolute():
        index = (work / index).resolve()
    if not (index / "SAindex").is_file():
        raise ValueError(
            "A pre-built STAR index is required. The adapter never builds or downloads an index: " + str(index)
        )

    install = subprocess.run(["pixi", "install", "--locked", "--all"], cwd=work, check=False)
    if install.returncode:
        raise RuntimeError(f"Pinned upstream Pixi installation failed with exit code {install.returncode}")
    command = [
        "pixi", "run", "-e", "snakemake", "snakemake", "--cores", str(workers),
        "--allowed-rules", "fastqc_raw", "fastp", "fastqc_clean", "star_align",
        "index_bam", "megadepth_bigwig", "merge_counts",
    ]
    if dry_run:
        command.append("--dry-run")
    command.append(TARGET)
    completed = subprocess.run(command, cwd=work, check=False)
    if completed.returncode:
        raise RuntimeError(f"Pinned upstream count target failed with exit code {completed.returncode}")
    output = work / TARGET
    result = {
        "status": "dry_run_complete" if dry_run else "complete",
        "provider": PROVIDER,
        "commit": PINNED_COMMIT,
        "command": command,
        "work_directory": str(work),
        "count_matrix": str(output),
    }
    if not dry_run:
        if not output.is_file():
            raise RuntimeError("Upstream command completed without the declared count matrix")
        result["sha256"] = sha256_file(output)
        (work / "adapter_result.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    return result
