#!/usr/bin/env python3
from __future__ import annotations

import argparse
import datetime as dt
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from internal.lib.core import (analysis_decision_rows, compare_de, environment_resources,
                               inventory_legacy_outputs, import_counts, import_resource_pack,
                               load_yaml, module_confirmation_string, plan_confirmation_token,
                               validate_project, validate_qc_preview, validate_resources)
from internal.lib.core import clean_preview, sync_declared_resources


def main() -> int:
    parser = argparse.ArgumentParser(prog="pipeline")
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("doctor")
    init_project = sub.add_parser("init-project")
    init_project.add_argument("--project-id", required=True)
    analyze = sub.add_parser("analyze")
    analyze.add_argument("--project", type=Path, required=True)
    analyze.add_argument("--run-id")
    analyze.add_argument("--publish", action="store_true")
    analyze.add_argument("--resume", action="store_true", help="Resume the exact existing Snakemake staging run")
    analyze.add_argument("--confirm-modules", help="Exact selection printed by `pipeline modules`")
    analyze.add_argument("--confirm-plan", help="Exact one-run token printed by `pipeline plan`")
    modules = sub.add_parser("modules")
    modules.add_argument("--project", type=Path, required=True)
    plan = sub.add_parser("plan")
    plan.add_argument("--project", type=Path, required=True)
    qc_preview = sub.add_parser("qc-preview")
    qc_preview.add_argument("--project", type=Path, required=True)
    qc_preview.add_argument("--run-id")
    qc_preview.add_argument("--confirm-plan", required=True,
                            help="Exact token printed by `pipeline plan` after confirming input/design/filtering")
    check = sub.add_parser("check")
    check.add_argument("--project", type=Path, required=True)
    inv = sub.add_parser("inventory")
    inv.add_argument("--legacy-results", type=Path, required=True)
    inv.add_argument("--new-results", type=Path)
    inv.add_argument("--output", type=Path, required=True)
    comp = sub.add_parser("compare")
    comp.add_argument("--legacy-de", type=Path, required=True)
    comp.add_argument("--new-de", type=Path, required=True)
    comp.add_argument("--output", type=Path, required=True)
    comp.add_argument("--padj", type=float, default=0.05)
    comp.add_argument("--abs-lfc", type=float, default=1.0)
    res = sub.add_parser("import-resources")
    res.add_argument("--source", type=Path, required=True)
    res.add_argument("--destination", type=Path, default=ROOT / "resources/data/imported")
    res.add_argument("--registry", type=Path, default=ROOT / "resources/registry.yml")
    imp = sub.add_parser("import-counts")
    imp.add_argument("--project", type=Path, required=True)
    imp.add_argument("--source", type=Path, required=True)
    imp.add_argument("--destination-name")
    sync = sub.add_parser("resources-sync")
    sync.add_argument("--manifest", type=Path, required=True)
    sync.add_argument("--destination", type=Path, required=True)
    sync.add_argument("--registry", type=Path, default=ROOT / "resources/registry.yml")
    clean = sub.add_parser("clean")
    clean.add_argument("--project-id", required=True)
    clean.add_argument("--scope", action="append", choices=["staging", "failed", "cache"], required=True)
    clean.add_argument("--confirm", help="Exact token DELETE:<project_id>; omit for mandatory dry-run")
    upstream_plan = sub.add_parser("upstream-plan")
    upstream_plan.add_argument("--project", type=Path, required=True)
    upstream_run = sub.add_parser("upstream-run")
    upstream_run.add_argument("--project", type=Path, required=True)
    upstream_run.add_argument("--confirm-upstream", required=True)
    upstream_run.add_argument("--workers", type=int, default=8)
    upstream_run.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    if args.command == "doctor":
        resources = environment_resources()
        resources["pipeline_root"] = str(ROOT)
        runtime_dirs = {"input", "results", "work", "data"}
        source_files = [
            p for p in ROOT.rglob("*")
            if p.is_file()
            and ".pixi" not in p.parts
            and ".git" not in p.parts
            and p.name != "pixi.lock"
            and not runtime_dirs.intersection(p.parts)
            and "internal/upstream/omics-pipelines" not in p.relative_to(ROOT).as_posix()
        ]
        hardcoded = []
        for path in source_files:
            try:
                content = path.read_text(encoding="utf-8")
            except UnicodeDecodeError:
                continue
            linux_personal = "/home/" + "azhenfan"
            windows_personal = "C:" + "\\Users\\AF"
            if linux_personal in content or windows_personal in content:
                hardcoded.append(str(path.relative_to(ROOT)))
        resources["path_is_portable"] = not hardcoded
        resources["hardcoded_path_files"] = hardcoded
        print(json.dumps(resources, indent=2))
        return 0 if not hardcoded else 2
    if args.command == "init-project":
        project_id = args.project_id
        if not project_id.replace("-", "").replace("_", "").replace(".", "").isalnum():
            raise ValueError("project-id may contain only letters, numbers, dot, dash and underscore")
        source = ROOT / "projects/_template"
        destination = ROOT / "projects" / project_id
        if destination.exists():
            raise ValueError(f"Project already exists: {destination}")
        shutil.copytree(source, destination)
        config_path = destination / "project.yml"
        cfg = load_yaml(config_path)
        cfg["project_id"] = project_id
        config_path.write_text(yaml.safe_dump(cfg, sort_keys=False, allow_unicode=True), encoding="utf-8")
        print(json.dumps({"status": "created_unconfirmed", "project": str(destination),
                          "next": "Ask the user the initialization questions in AGENTS.md; do not analyze yet."}, indent=2))
        return 0
    if args.command == "analyze":
        project = args.project.resolve()
        cfg = validate_project(project / "project.yml")
        migration = cfg.get("migration", {})
        expected_modules = module_confirmation_string(cfg)
        expected_plan = plan_confirmation_token(project, cfg)
        supplied_modules = args.confirm_modules
        if not supplied_modules and sys.stdin.isatty():
            print("Confirmed module selection required for this analysis:")
            print(expected_modules)
            supplied_modules = input("Type the exact selection to continue: ").strip()
        if supplied_modules != expected_modules:
            print(json.dumps({
                "status": "stopped_unconfirmed_modules",
                "module_selection": expected_modules,
                "required_argument": f"--confirm-modules '{expected_modules}'",
            }, indent=2))
            return 2
        supplied_plan = args.confirm_plan
        if not supplied_plan and sys.stdin.isatty():
            print("Fresh approval of the complete scientific plan is required:")
            print(expected_plan)
            supplied_plan = input("Type the exact plan token to continue: ").strip()
        if supplied_plan != expected_plan:
            print(json.dumps({
                "status": "stopped_unconfirmed_plan",
                "reason": "The token binds this run to project.yml, samples, contrasts, QC approval and counts manifest.",
                "required_command": f"pixi run plan -- --project {project}",
                "required_argument": f"--confirm-plan '{expected_plan}'",
            }, indent=2))
            return 2
        run_id = args.run_id or dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ-v4")
        if not run_id.replace("-", "").replace("_", "").replace(".", "").isalnum():
            raise ValueError("run-id may contain only letters, numbers, dot, dash and underscore")
        staging = project / "work" / "staging" / run_id
        if staging.exists() and not args.resume:
            raise ValueError(f"Staging run already exists: {staging}")
        if args.resume and not staging.exists():
            raise ValueError(f"Cannot resume a staging run that does not exist: {staging}")
        resources = environment_resources()
        counts = project / cfg["inputs"]["counts_file"]
        estimate = {
            "project_id": cfg["project_id"], "run_id": run_id,
            "recommended_workers": resources["recommended_workers"],
            "threads_per_worker": 1,
            "available_memory_gib": resources["available_memory_gib"],
            "input_bytes": counts.stat().st_size,
            "estimated_staging_bytes": max(1_000_000_000, counts.stat().st_size * 8),
            "staging": str(staging),
        }
        print(json.dumps({"analysis_estimate": estimate}, indent=2), flush=True)
        env = os.environ.copy()
        env.update({"OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1",
                    "BULK_RNASEQ_MODULE_CONFIRMATION": expected_modules,
                    "BULK_RNASEQ_PLAN_CONFIRMATION": expected_plan})
        locked_bin = ROOT / ".pixi/envs/default/bin"
        env["PATH"] = str(locked_bin) + os.pathsep + env.get("PATH", "")
        locked_snakemake = ROOT / ".pixi/envs/default/bin/snakemake"
        snakemake = str(locked_snakemake) if locked_snakemake.is_file() else shutil.which("snakemake")
        if not snakemake:
            raise ValueError("Snakemake is not available inside the locked Pixi environment")
        workers = min(cfg.get("compute", {}).get("max_workers", 3), resources["recommended_workers"], 3)
        completed = subprocess.run([
            snakemake, "--snakefile", str(ROOT / "workflow/Snakefile"),
            "--config", f"project={project}", f"run_dir={staging}",
            "--profile", str(ROOT / "workflow/profiles/local"),
            "--cores", str(workers), "--jobs", str(workers),
            "--resources", f"mem_mb={max(1024, int(resources['available_memory_gib'] * 1024 * 0.75))}",
        ], cwd=ROOT, env=env, check=False)
        if completed.returncode:
            print(json.dumps({"status": "failed_explicit", "staging": str(staging)}, indent=2))
            return completed.returncode
        if migration.get("requires_legacy_parity"):
            print(json.dumps({
                "status": "complete_staging_requires_legacy_regression",
                "staging": str(staging),
                "reason": "A migrated project may be published only after the file-level and scientific regression gate passes.",
            }, indent=2))
            return 0
        if args.publish:
            token = f"PUBLISH:{cfg['project_id']}:{run_id}"
            completed = subprocess.run([
                sys.executable, str(ROOT / "internal/reporting/publish_project_results.py"),
                str(staging), str(project), "--confirm", token,
            ], cwd=ROOT, check=False)
            return completed.returncode
        print(json.dumps({"status": "complete", "staging": str(staging),
                          "publish_command": f"pixi run publish-results -- {staging} {project} --confirm PUBLISH:{cfg['project_id']}:{run_id}"}, indent=2))
        return 0
    if args.command == "qc-preview":
        project = args.project.resolve()
        cfg = validate_qc_preview(project / "project.yml")
        expected_plan = plan_confirmation_token(project, cfg)
        if args.confirm_plan != expected_plan:
            print(json.dumps({
                "status": "stopped_unconfirmed_qc_plan",
                "required_command": f"pixi run plan -- --project {project}",
                "required_argument": f"--confirm-plan '{expected_plan}'",
            }, indent=2))
            return 2
        run_id = args.run_id or dt.datetime.now(dt.timezone.utc).strftime("qc-preview-%Y%m%dT%H%M%SZ")
        if not run_id.replace("-", "").replace("_", "").replace(".", "").isalnum():
            raise ValueError("run-id may contain only letters, numbers, dot, dash and underscore")
        staging = project / "work/staging" / run_id
        if staging.exists():
            raise ValueError(f"QC preview staging already exists: {staging}")
        env = os.environ.copy()
        env.update({
            "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1",
            "BULK_RNASEQ_QC_PREVIEW": "1", "BULK_RNASEQ_PLAN_CONFIRMATION": expected_plan,
        })
        locked_bin = ROOT / ".pixi/envs/default/bin"
        env["PATH"] = str(locked_bin) + os.pathsep + env.get("PATH", "")
        rscript = locked_bin / "Rscript"
        completed = subprocess.run([
            str(rscript), str(ROOT / "workflow/snakemake/qc_preview.R"), str(project), str(staging),
        ], cwd=ROOT, env=env, check=False)
        if completed.returncode:
            print(json.dumps({"status": "failed_explicit", "qc_preview": str(staging)}, indent=2))
            return completed.returncode
        print(json.dumps({
            "status": "awaiting_user_qc_approval", "qc_preview": str(staging),
            "next": "Inspect every QC figure, then explicitly update samples.tsv and qc_approval.yml before generating a new formal plan.",
        }, indent=2))
        return 0
    if args.command == "modules":
        project_file = args.project / "project.yml" if args.project.is_dir() else args.project
        cfg = load_yaml(project_file)
        selection = module_confirmation_string(cfg)
        rows = [{"module": name, **decision} for name, decision in cfg["analysis"]["modules"].items()]
        print(json.dumps({"project_id": cfg["project_id"], "modules": rows,
                          "confirmation_token": selection}, indent=2, ensure_ascii=False))
        return 0
    if args.command == "plan":
        project = args.project.resolve()
        project_file = project / "project.yml" if project.is_dir() else project
        cfg = load_yaml(project_file)
        project_dir = project_file.parent
        print(json.dumps({
            "project_id": cfg.get("project_id"),
            "status": "awaiting_fresh_user_confirmation",
            "project_decision_confirmation": cfg.get("decision_confirmation", {}),
            "decisions": analysis_decision_rows(project_dir, cfg),
            "module_confirmation_token": module_confirmation_string(cfg),
            "plan_confirmation_token": plan_confirmation_token(project_dir, cfg),
            "instruction": "Show every row to the user. Run only after the user confirms this exact plan in the current conversation.",
        }, indent=2, ensure_ascii=False))
        return 0
    if args.command == "check":
        cfg = validate_project(args.project / "project.yml" if args.project.is_dir() else args.project)
        resources = validate_resources(ROOT / "resources/registry.yml")
        print(json.dumps({"project": cfg["project_id"], "status": "valid", "frozen_resources": len(resources)}, indent=2))
        return 0
    if args.command == "inventory":
        _, _, summary = inventory_legacy_outputs(args.legacy_results, args.new_results, args.output)
        print(json.dumps(summary, indent=2))
        return 0 if summary["passed"] else 2
    if args.command == "compare":
        _, summary = compare_de(args.legacy_de, args.new_de, args.output, args.padj, args.abs_lfc)
        print(json.dumps(summary, indent=2))
        return 0 if summary["passed"] else 2
    if args.command == "import-resources":
        imported = import_resource_pack(args.source, args.destination, args.registry)
        print(json.dumps({"imported": len(imported), "registry": str(args.registry)}, indent=2))
        return 0
    if args.command == "import-counts":
        destination, manifest = import_counts(args.source, args.project, args.destination_name)
        print(json.dumps({"imported": str(destination), "sha256": manifest["sha256"], "genes": manifest["genes"], "samples": manifest["samples"], "limitation": manifest["limitation"]}, indent=2))
        return 0
    if args.command == "resources-sync":
        synced = sync_declared_resources(args.manifest, args.destination, args.registry)
        print(json.dumps({"synced": len(synced), "resources": [{"resource_id": x["resource_id"], "sha256": x["sha256"]} for x in synced]}, indent=2))
        return 0
    if args.command == "clean":
        preview = clean_preview(args.project_id, args.scope)
        print(json.dumps({"dry_run": args.confirm is None, "targets": preview, "confirmation_token": f"DELETE:{args.project_id}"}, indent=2))
        if args.confirm is None:
            return 0
        if args.confirm != f"DELETE:{args.project_id}":
            raise ValueError("Confirmation token does not match the project ID")
        for item in preview:
            target = Path(item["target"])
            if target.exists():
                shutil.rmtree(target)
        print(json.dumps({"deleted": [item["target"] for item in preview if item["exists"]]}, indent=2))
        return 0
    if args.command in {"upstream-plan", "upstream-run"}:
        from internal.upstream.adapter import plan_upstream, run_upstream
        if args.command == "upstream-plan":
            print(json.dumps(plan_upstream(args.project.resolve()), indent=2, ensure_ascii=False))
            return 0
        result = run_upstream(args.project.resolve(), args.confirm_upstream, args.workers, args.dry_run)
        print(json.dumps(result, indent=2, ensure_ascii=False))
        return 0
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
