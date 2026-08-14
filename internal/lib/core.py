from __future__ import annotations

import csv
import hashlib
import json
import os
import re
import shutil
import subprocess
import urllib.request
from datetime import datetime, timezone
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import pandas as pd
import yaml
from jsonschema import Draft202012Validator, FormatChecker

ROOT = Path(__file__).resolve().parents[2]

NORMALIZED_INPUTS = {"TPM", "FPKM", "CPM", "batch_corrected", "other_normalized"}
ALLOWED_LEGACY_STATUS = {
    "preserved",
    "improved_and_replaced",
    "not_applicable",
    "explicitly_skipped_by_user",
    "missing",
}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_yaml(path: Path) -> dict:
    with path.open(encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def module_confirmation_string(cfg: dict) -> str:
    """Return the exact, human-readable module selection required for every run."""
    modules = cfg["analysis"]["modules"]
    ordered = sorted(modules)
    return ",".join(f"{name}={modules[name]['status']}" for name in ordered)


DECISION_GROUPS = (
    "input_design_filtering",
    "differential",
    "enrichment",
    "regulation",
    "network",
    "motif",
    "exploratory",
    "reporting_export",
)


def plan_confirmation_token(project: Path, cfg: dict | None = None) -> str:
    """Bind one-run approval to every scientific/configuration input."""
    project = project.resolve()
    if project.is_file():
        project_file = project
        project = project.parent
    else:
        project_file = project / "project.yml"
    cfg = cfg or load_yaml(project_file)
    inputs = cfg.get("inputs", {})
    paths = [
        project_file,
        project / inputs.get("samples_file", "samples.tsv"),
        project / inputs.get("contrasts_file", "contrasts.tsv"),
        project / inputs.get("qc_approval_file", "qc_approval.yml"),
        project / "input/source_manifest.yml",
    ]
    material = [cfg.get("project_id", "UNKNOWN")]
    material.extend(sha256_file(path) if path.is_file() else f"MISSING:{path.name}" for path in paths)
    digest = hashlib.sha256("|".join(material).encode("utf-8")).hexdigest()[:24]
    return f"PLAN:{cfg.get('project_id', 'UNKNOWN')}:{digest}"


def analysis_decision_rows(project: Path, cfg: dict | None = None) -> list[dict]:
    """Return the complete user-facing decision inventory for a formal run."""
    project = project.resolve()
    if project.is_file():
        project_file = project
        project = project.parent
    else:
        project_file = project / "project.yml"
    cfg = cfg or load_yaml(project_file)
    inputs = cfg.get("inputs", {})

    def row(group: str, decision: str, value, impact: str) -> dict:
        return {"group": group, "decision": decision, "value": value, "impact": impact}

    rows = [
        row("input_design_filtering", "counts_provenance", cfg.get("counts_provenance"),
            "Determines whether DESeq2 is allowed and records upstream limitations."),
        row("input_design_filtering", "species_and_assembly",
            {"species": cfg.get("species"), "assembly": cfg.get("assembly")},
            "Controls annotation and all species-specific frozen resources."),
        row("input_design_filtering", "design", cfg.get("design"),
            "Controls the fitted model, reference levels, covariates and interpretation."),
        row("input_design_filtering", "gene_filter", cfg.get("filtering"),
            "Controls which genes enter dispersion estimation and every downstream test."),
        row("input_design_filtering", "upstream", cfg.get("upstream"),
            "Controls whether the separately gated FASTQ-to-counts adapter is used."),
        row("differential", "deg_definition", cfg.get("thresholds"),
            "Controls DEG membership and threshold-dependent ORA/network/motif inputs."),
        row("differential", "differential_display", cfg.get("differential", {}),
            "Controls labels and heatmap display only; it must not change the fitted model."),
        row("enrichment", "ora_and_gsea", cfg.get("enrichment", {}),
            "Controls ORA gene lists/FDR, databases, set sizes, GSEA reporting and display windows."),
        row("regulation", "tf_gsva_progeny", cfg.get("regulation", {}),
            "Controls TFs/pathways of interest and reported activity summaries."),
        row("network", "ppi", cfg.get("network", {}),
            "Controls the DEG profile, STRING confidence, graph size, hubs and module summaries."),
        row("motif", "promoter_and_peak_motif", cfg.get("motif", {}),
            "Controls foreground definition, promoter window, matched background and motif thresholds."),
        row("exploratory", "target_genes_custom_pathview", cfg.get("exploratory", {}),
            "Controls boxplots, custom gene sets, Pathview IDs and other declared targets."),
        row("reporting_export", "modules", cfg.get("analysis", {}).get("modules", {}),
            "Controls what is computed; every module needs an explicit enabled/NA/skipped state."),
        row("reporting_export", "report_and_export",
            {"report": cfg.get("report"), "export": cfg.get("export")},
            "Controls report depth and released figures without changing numerical analysis."),
    ]
    contrasts_path = project / inputs.get("contrasts_file", "contrasts.tsv")
    if contrasts_path.is_file():
        try:
            contrast_value = pd.read_csv(contrasts_path, sep="\t").to_dict(orient="records")
        except Exception as error:
            contrast_value = f"UNREADABLE: {error}"
    else:
        contrast_value = "MISSING"
    rows.insert(3, row("input_design_filtering", "contrasts", contrast_value,
                       "Controls numerator/denominator direction for every signed result."))
    samples_path = project / inputs.get("samples_file", "samples.tsv")
    if samples_path.is_file():
        try:
            sample_value = pd.read_csv(samples_path, sep="\t").to_dict(orient="records")
        except Exception as error:
            sample_value = f"UNREADABLE: {error}"
    else:
        sample_value = "MISSING"
    rows.insert(3, row("input_design_filtering", "samples_and_exclusions", sample_value,
                       "Controls biological replication, covariates, pairing and approved exclusions."))
    approval_path = project / inputs.get("qc_approval_file", "qc_approval.yml")
    approval_value = load_yaml(approval_path) if approval_path.is_file() else "MISSING"
    rows.insert(5, row("input_design_filtering", "qc_approval", approval_value,
                       "Confirms that the current count hash and all sample exclusions were reviewed."))
    rows.append(row("reporting_export", "frozen_resources", cfg.get("resources", {}),
                    "Pins annotation, gene sets, TF networks, STRING and motif resources by registry ID."))
    rows.append(row("reporting_export", "compute", cfg.get("compute", {}),
                    "Caps workers and threads; the doctor command may lower concurrency for available memory."))
    return rows


def validate_project(path: Path) -> dict:
    cfg = load_yaml(path)
    schema = json.loads((ROOT / "internal/schemas/project.schema.json").read_text(encoding="utf-8"))
    errors = sorted(Draft202012Validator(schema).iter_errors(cfg), key=lambda e: list(e.path))
    if errors:
        detail = "\n".join(f"- {'.'.join(map(str, e.path)) or '<root>'}: {e.message}" for e in errors)
        raise ValueError(f"Invalid project configuration:\n{detail}")
    provenance = cfg["counts_provenance"]
    if not provenance.get("confirmed"):
        raise ValueError("Counts provenance is unconfirmed")
    if not cfg["design"].get("confirmed"):
        raise ValueError("Experimental design is unconfirmed")
    if provenance["normalization"] in NORMALIZED_INPUTS:
        raise ValueError(f"DESeq2 is forbidden for {provenance['normalization']} input")
    if provenance["normalization"] != "raw_counts":
        raise ValueError("Formal differential analysis requires explicitly confirmed raw_counts; unknown/other input must be resolved first")
    if provenance["feature_level"] == "transcript":
        raise ValueError("Transcript-level input must be summarized to gene-level counts with an explicit method before DESeq2")
    if (provenance["technical_replicates"] == "present" and
            provenance.get("technical_replicate_handling") != "already_summed_upstream"):
        raise ValueError("Technical replicates cannot be treated as biological replicates; confirm upstream summation first")
    confirmations = cfg.get("decision_confirmation", {})
    undecided_groups = [name for name in DECISION_GROUPS if not confirmations.get(name)]
    if undecided_groups:
        raise ValueError("Unconfirmed scientific decision groups: " + ", ".join(undecided_groups))
    undecided = [name for name, decision in cfg["analysis"]["modules"].items() if not decision.get("confirmed")]
    if undecided:
        raise ValueError("Unconfirmed modules: " + ", ".join(undecided))
    upstream = cfg["upstream"]
    if not upstream.get("confirmed") or upstream.get("status") == "unconfirmed":
        raise ValueError("Upstream-data availability is unconfirmed")
    if upstream["status"] == "enabled" and upstream.get("provider") != "xuzhougeng/omics-pipelines":
        raise ValueError("The only supported upstream provider is xuzhougeng/omics-pipelines")
    modules = cfg["analysis"]["modules"]
    unsupported = {
        "motif_peaks": "Peak-aware motif execution is not implemented yet",
        "wgcna": "WGCNA execution is not implemented yet",
        "personalized": "Personalized analysis requires a separately reviewed project-specific executor",
    }
    enabled_unsupported = [f"{name}: {reason}" for name, reason in unsupported.items()
                           if modules[name]["status"] == "enabled"]
    if enabled_unsupported:
        raise ValueError("Enabled modules cannot be silently declared complete:\n- " + "\n- ".join(enabled_unsupported))
    migration = cfg.get("migration", {})
    if (modules["pathview"]["status"] == "enabled" and
            migration.get("scientific_baseline") == "v4_workflow"):
        raise ValueError("Generic Pathview execution is not implemented; frozen legacy map copying is allowed only for migration parity")
    if cfg["species"] != "human":
        human_only = [name for name in ("enrichment", "regulation", "network")
                      if modules[name]["status"] == "enabled"]
        if human_only:
            raise ValueError("Enabled modules currently require human frozen resources: " + ", ".join(human_only))
    profiles = cfg["thresholds"].get("profiles") or {
        "Primary": {"padj": cfg["thresholds"].get("padj"), "abs_lfc": cfg["thresholds"].get("abs_lfc")}
    }
    differential = cfg.get("differential", {})
    primary_profile = differential.get("primary_profile", next(iter(profiles)))
    if primary_profile not in profiles:
        raise ValueError(f"differential.primary_profile is not a declared threshold profile: {primary_profile}")
    if modules["enrichment"]["status"] == "enabled":
        enrichment = cfg.get("enrichment", {})
        ora_profiles = enrichment.get("ora_profiles", [primary_profile])
        unknown_ora_profiles = [name for name in ora_profiles if name not in profiles]
        if unknown_ora_profiles:
            raise ValueError("enrichment.ora_profiles contains undeclared threshold profiles: " +
                             ", ".join(unknown_ora_profiles))
        if not enrichment.get("ora_databases") or not enrichment.get("gsea_databases"):
            raise ValueError("Enabled enrichment requires explicit ora_databases and gsea_databases")
        if enrichment.get("min_size", 10) > enrichment.get("max_size", 500):
            raise ValueError("enrichment.min_size cannot exceed enrichment.max_size")
    regulation = cfg.get("regulation", {})
    if regulation.get("tf_gsea_min_size", 10) > regulation.get("tf_gsea_max_size", 500):
        raise ValueError("regulation.tf_gsea_min_size cannot exceed tf_gsea_max_size")
    if modules["network"]["status"] == "enabled":
        network_profile = cfg.get("network", {}).get("deg_profile", primary_profile)
        if network_profile not in profiles:
            raise ValueError(f"network.deg_profile is not a declared threshold profile: {network_profile}")
    if modules["motif_promoter"]["status"] == "enabled":
        motif_profile = cfg.get("motif", {}).get("deg_profile", primary_profile)
        if motif_profile not in profiles:
            raise ValueError(f"motif.deg_profile is not a declared threshold profile: {motif_profile}")
        motif = cfg.get("motif", {})
        if motif.get("streme_min_width", 6) > motif.get("streme_max_width", 15):
            raise ValueError("motif.streme_min_width cannot exceed motif.streme_max_width")
    if modules["custom_gene_sets"]["status"] == "enabled":
        custom = cfg.get("exploratory", {}).get("custom_gene_sets", {})
        if custom.get("formal_min_size", 3) > custom.get("max_size", 500):
            raise ValueError("custom_gene_sets formal_min_size cannot exceed max_size")
    return cfg


def validate_qc_preview(path: Path) -> dict:
    """Validate only the decisions needed to render pre-approval sample QC."""
    cfg = load_yaml(path)
    schema = json.loads((ROOT / "internal/schemas/project.schema.json").read_text(encoding="utf-8"))
    errors = sorted(Draft202012Validator(schema).iter_errors(cfg), key=lambda e: list(e.path))
    if errors:
        detail = "\n".join(f"- {'.'.join(map(str, e.path)) or '<root>'}: {e.message}" for e in errors)
        raise ValueError(f"Invalid project configuration:\n{detail}")
    provenance = cfg["counts_provenance"]
    if not provenance.get("confirmed") or provenance.get("normalization") != "raw_counts":
        raise ValueError("QC preview requires explicitly confirmed raw_counts provenance")
    if provenance.get("feature_level") == "transcript":
        raise ValueError("QC preview requires a confirmed gene-level count matrix")
    if (provenance.get("technical_replicates") == "present" and
            provenance.get("technical_replicate_handling") != "already_summed_upstream"):
        raise ValueError("Technical replicates must be resolved before QC preview")
    if not cfg["design"].get("confirmed"):
        raise ValueError("QC preview requires a confirmed sample design and factor levels")
    if not cfg.get("decision_confirmation", {}).get("input_design_filtering"):
        raise ValueError("QC preview input/design/filtering decisions are unconfirmed")
    if not cfg["upstream"].get("confirmed") or cfg["upstream"].get("status") == "unconfirmed":
        raise ValueError("QC preview requires an explicit upstream-data decision")
    qc = cfg["analysis"]["modules"]["qc"]
    if not qc.get("confirmed") or qc.get("status") != "enabled":
        raise ValueError("QC preview requires the QC module to be explicitly enabled and confirmed")
    return cfg


def validate_resources(registry_path: Path, *, require_local_files: bool = True) -> list[dict]:
    """Validate frozen-resource metadata and, for runnable checks, local files.

    A clean Git clone intentionally excludes licensed/large resource payloads. CI can
    therefore validate the committed registry with ``require_local_files=False``;
    every operational caller keeps the strict default and verifies existence plus
    SHA256 before counts import, analysis, or resource publication.
    """
    registry = load_yaml(registry_path) or {}
    schema = json.loads((ROOT / "internal/schemas/resource.schema.json").read_text(encoding="utf-8"))
    validator = Draft202012Validator(schema, format_checker=FormatChecker())
    problems = []
    for resource in registry.get("resources", []):
        for error in validator.iter_errors(resource):
            problems.append(f"{resource.get('resource_id', '<unknown>')}: {error.message}")
        local_path = resource.get("local_path")
        if local_path and require_local_files:
            path = (registry_path.parent / local_path).resolve()
            if not path.exists():
                problems.append(f"{resource.get('resource_id')}: local file is missing: {path}")
            elif sha256_file(path).lower() != str(resource.get("sha256", "")).lower():
                problems.append(f"{resource.get('resource_id')}: SHA256 mismatch")
    if problems:
        raise ValueError("Invalid resource registry:\n- " + "\n- ".join(problems))
    return registry.get("resources", [])


def write_tsv(rows: Iterable[dict], path: Path, fieldnames: list[str] | None = None) -> Path:
    rows = list(rows)
    path.parent.mkdir(parents=True, exist_ok=True)
    if fieldnames is None:
        fieldnames = list(rows[0]) if rows else []
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    return path


def all_files(root: Path) -> list[Path]:
    return sorted(p for p in root.rglob("*") if p.is_file())


@dataclass(frozen=True)
class OutputFamily:
    analysis_id: str
    family: str
    old_pattern: str
    new_patterns: tuple[str, ...]
    replacement: bool = False


OUTPUT_FAMILIES = (
    OutputFamily("qc", "library_size", r"Library_Size_Barplot", ("Library_Size",), True),
    OutputFamily("qc", "initial_pca", r"PCA_plot_initial", ("PCA", "Initial_PCA")),
    OutputFamily("qc", "final_pca", r"PCA_plot_batch_corrected|PCA_plot_final", ("Final_PCA", "PCA_final", "PCA_batch_corrected")),
    OutputFamily("qc", "sample_clustering", r"Sample_Hierarchical_Clustering|Sample_Clustering_Order", ("Sample_Hierarchical_Clustering", "Sample_Clustering_Order")),
    OutputFamily("qc", "correlation_qc_barplot", r"Sample_Correlation_QC_Barplot", ("Sample_Correlation_QC_Barplot",)),
    OutputFamily("qc", "sample_correlation", r"Sample_Correlation_(Heatmap|Pearson)", ("Sample_Correlation",)),
    OutputFamily("qc", "sample_distance", r"Sample_Distance_(Heatmap|Euclidean)", ("Sample_Distance",)),
    OutputFamily("differential", "full_de_table", r"DEG_results_annotated", ("Differential_Full",)),
    OutputFamily("differential", "raw_ma", r"MA.*raw|MA_Plot_Raw", ("MA_Raw", "Raw_MA")),
    OutputFamily("differential", "shrunk_ma", r"MA.*shr|MA_Plot_Shr", ("MA_Shrunk", "Shrunken_MA", "MA_Raw_and_Shrunk"), True),
    OutputFamily("differential", "volcano", r"Volcano", ("Volcano",)),
    OutputFamily("differential", "deg_heatmap", r"Heatmap_DEGs|Heatmap_Top_Genes", ("Top_DEG_Heatmap", "Heatmap_Top_Genes")),
    OutputFamily("differential", "gene_boxplot", r"Boxplot", ("Boxplot",)),
    OutputFamily("differential", "de_workbook", r"Differential_Expression\.xlsx", ("Differential_Full", "Differential_.*Significant"), True),
    OutputFamily("differential", "de_significant_table", r"DEG_results_significant|DEG_model_summary", ("Differential_.*Significant", "Design_Matrix"), True),
    OutputFamily("enrichment", "ora_full_tables", r"Enrichment_Full_", (r"Enrichment_Full_",), True),
    OutputFamily("enrichment", "gsea_full_tables", r"03_Enrichment/GSEA/Tables/GSEA_Full_Table_", ("03_Enrichment/.*/GSEA_Full",)),
    OutputFamily("enrichment", "gsea_curves", r"03_Enrichment/GSEA/.*/Figures/[0-9]+_(UP|DOWN)_", ("03_Enrichment/.*/GSEA/Figures/.*/[0-9]+_",), True),
    OutputFamily("enrichment", "gsea_overviews", r"03_Enrichment/GSEA/.*/Figures/.*NES_Overview|GSEA_Summary_Dotplots", ("NES_Overview",), True),
    OutputFamily("enrichment", "gsea_provenance", r"GSEA_Gene_Set_Provenance|GSEA_ranked_gene_list", ("Gene_Set_Provenance", "GSEA_Ranked_Gene_List"), True),
    OutputFamily("enrichment", "lollipop", r"Lollipop", (r"Lollipop",), True),
    OutputFamily("enrichment", "sankey_bubble", r"sankey_bubble", ("Sankey_Bubble",), True),
    OutputFamily("enrichment", "cnet", r"Cnet", ("Cnet",), True),
    OutputFamily("enrichment", "emap", r"Emap", ("Emap",), True),
    OutputFamily("enrichment", "apear", r"aPEAR|apear", ("aPEAR", "apear")),
    OutputFamily("enrichment", "ora_classic_plots", r"ORA/.*/Figures/.*(Dotplot|Diverging_Bar|Heat)", (r"(Dotplot|Diverging_Bar|_Heat)",), True),
    OutputFamily("enrichment", "ora_reduced_tables", r"Enrichment_Reduced_", ("Enrichment_Reduced_",), True),
    OutputFamily("enrichment", "ora_provenance", r"ORA_Gene_Set_Provenance|ORA_universe_entrez", ("Gene_Set_Provenance", "ORA_Universe_Entrez"), True),
    OutputFamily("enrichment", "sankey_tables", r"ORA/Sankey/Tables/", (r"Tables/Sankey/.*Filtered_Enrichment",), True),
    OutputFamily("regulation", "chea_gsea", r"ChEA.*GSEA|GSEA_Full_Table_ChEA", ("ChEA",)),
    OutputFamily("regulation", "encode_gsea", r"ENCODE.*GSEA|GSEA_Full_Table_ENCODE", ("ENCODE",)),
    OutputFamily("regulation", "gtrd_gsea", r"GTRD.*GSEA|GSEA_Full_Table_GTRD", ("GTRD",)),
    OutputFamily("regulation", "single_tf_report", r"04_Regulation/TF_Activity/Figures/(single_tf_plots|Global_TF_Analysis_Report)", ("04_Regulation/.*/single_tf", "04_Regulation/.*/TF_Analysis_Report")),
    OutputFamily("regulation", "gsva_tf_correlation", r"GSVA|TF_Pathway_Correlation", ("GSVA", "TF_Pathway_Correlation")),
    OutputFamily("regulation", "progeny", r"PROGENy|pathway_activity|target_genes", ("PROGENy", "Pathway_Activity", "Target_Genes")),
    OutputFamily("regulation", "tf_activity_tables", r"TF_Activity_(Contrast_Full|Sample_Matrix)", ("TF_Activity_Contrast", "TF_Activity_Sample_Matrix"), True),
    OutputFamily("regulation", "tf_gsea_summaries", r"Summary_Dotplot_(ChEA|ENCODE)", ("Summary_Directional_ChEA", "Summary_Directional_ENCODE"), True),
    OutputFamily("network", "ppi_hub", r"PPI_Hub", ("PPI_Hub",)),
    OutputFamily("network", "ppi_modules", r"PPI_Module|module_networks", ("PPI_Module", "module_networks")),
    OutputFamily("network", "ppi_module_enrichment", r"PPI_Module_Enrichment", ("PPI_Module_Enrichment",)),
    OutputFamily("network", "cytoscape", r"Network_Edges|Node_Attributes", ("Network_Edges", "Node_Attributes")),
    OutputFamily("network", "ppi_explicit_skip", r"SKIPPED_ppi_due_to_network_or_data", ("05_Network_status.tsv",), True),
    OutputFamily("custom", "custom_gsea_curves", r"06_Custom.*GSEA_", (r"07_Exploratory/Custom_Gene_Sets/Figures/Curves",), True),
    OutputFamily("custom", "custom_heatmaps", r"06_Custom.*Heatmap_", ("07_Exploratory/Custom_Gene_Sets/Figures/Heatmaps/Heatmap_",)),
    OutputFamily("custom", "custom_term2gene", r"custom_gene_sets_term2gene", ("Custom_Gene_Sets_TERM2GENE",), True),
    OutputFamily("pathview", "pathview", r"pathview_maps|hsa\d+\.\.png", ("Pathview", "pathview")),
    OutputFamily("provenance", "legacy_readmes", r"(^|/)README\.txt$|Analysis_Limitations|QC_REVIEW_REQUIRED|Excluded_Samples", ("run_manifest.json", "index.html", "status.tsv"), True),
    OutputFamily("provenance", "legacy_run_logs", r"99_Run/Logs/", ("Logs/", "Provenance/.*status.tsv", "run_manifest.json"), True),
    OutputFamily("provenance", "legacy_reports", r"99_Run/Reports/|(^|/)manifest\.csv$", ("index.html", "run_manifest.json"), True),
    OutputFamily("qc", "qc_support_tables", r"01_QC/Tables/(Gene_Filtering_Summary|Sample_QC_Metrics|metadata)", ("Gene_Filtering_Summary", "Sample_QC_Metrics", "run_manifest.json"), True),
)


def _matches(paths: list[Path], root: Path, pattern: str) -> list[str]:
    regex = re.compile(pattern, re.IGNORECASE)
    return [p.relative_to(root).as_posix() for p in paths if regex.search(p.relative_to(root).as_posix())]


def _legacy_specific_candidates(spec: OutputFamily, legacy_file: str, candidates: list[str]) -> list[str]:
    """Prevent one generic enrichment figure from satisfying unrelated legacy outputs."""
    if spec.analysis_id != "enrichment" or not candidates:
        return candidates

    def normalized(value: str) -> str:
        return re.sub(r"[^A-Z0-9]+", "_", value.upper()).strip("_")

    old = normalized(legacy_file)
    filtered = list(candidates)
    database_tokens = ("GO_BP", "GO_CC", "GO_MF", "KEGG", "HALLMARK", "REACTOME")
    for token in database_tokens:
        if token in old:
            filtered = [path for path in filtered if token in normalized(path)]
            break

    output_tokens = (("NES_OVERVIEW", "GENE_SET_PROVENANCE", "RANKED_GENE_LIST")
                     if "GSEA" in old else (
        "CNET", "EMAP", "HEAT", "DOTPLOT", "LOLLIPOP", "DIVERGING_BAR",
        "SANKEY_BUBBLE", "ENRICHMENT_FULL", "ENRICHMENT_REDUCED",
        "GENE_SET_PROVENANCE", "RANKED_GENE_LIST",
    ))
    for token in output_tokens:
        if token in old:
            filtered = [path for path in filtered if token in normalized(path)]
            break

    if any(token in old for token in ("DOTPLOT", "LOLLIPOP")):
        for direction in ("ALL", "DOWN", "UP"):
            if re.search(rf"(^|_){direction}(_|$)", old):
                filtered = [path for path in filtered
                            if re.search(rf"(^|_){direction}(_|$)", normalized(path))]
                break
    return sorted(set(filtered))


def inventory_legacy_outputs(legacy_root: Path, new_root: Path | None, output_dir: Path) -> tuple[Path, Path, dict]:
    old_files = all_files(legacy_root)
    new_files = all_files(new_root) if new_root and new_root.exists() else []
    family_rows = []
    file_rows: dict[str, dict] = {}
    analysis = {}
    for spec in OUTPUT_FAMILIES:
        old_hits = _matches(old_files, legacy_root, spec.old_pattern)
        if not old_hits:
            continue
        new_hits = []
        if new_root:
            for pattern in spec.new_patterns:
                new_hits.extend(_matches(new_files, new_root, pattern))
        new_hits = sorted(set(new_hits))
        per_file_hits = {legacy_file: _legacy_specific_candidates(spec, legacy_file, new_hits)
                         for legacy_file in old_hits}
        all_legacy_mapped = all(per_file_hits.values())
        status = ("improved_and_replaced" if spec.replacement else "preserved") if all_legacy_mapped else "missing"
        family_rows.append({
            "analysis_id": spec.analysis_id,
            "output_family": spec.family,
            "legacy_file_count": len(old_hits),
            "legacy_examples": " | ".join(old_hits[:5]),
            "new_file_count": len(new_hits),
            "new_examples": " | ".join(new_hits[:5]),
            "status": status,
            "reason": ("Every legacy file matched a database-, direction- and figure-specific output" if spec.analysis_id == "enrichment"
                       else ("Mapped to a consolidated or redesigned standard output" if spec.replacement
                             else "Matched by explicit output-family rule")) if all_legacy_mapped
                      else "At least one legacy file lacks a specific mapped v4/new output",
        })
        state = analysis.setdefault(spec.analysis_id, {"families": 0, "missing": 0, "legacy_files": 0, "new_files": 0})
        state["families"] += 1
        state["missing"] += status == "missing"
        state["legacy_files"] += len(old_hits)
        state["new_files"] += len(new_hits)
        for legacy_file in old_hits:
            mapped_hits = per_file_hits[legacy_file]
            file_status = ("improved_and_replaced" if spec.replacement else "preserved") if mapped_hits else "missing"
            existing = file_rows.get(legacy_file)
            candidate = {
                "analysis_id": spec.analysis_id,
                "output_family": spec.family,
                "legacy_file": legacy_file,
                "legacy_sha256": sha256_file(legacy_root / legacy_file),
                "new_file_count": len(mapped_hits),
                "new_files": " | ".join(mapped_hits),
                "status": file_status,
                "reason": ("Mapped to a consolidated or redesigned standard output" if spec.replacement
                           else "Matched by an explicit output-family rule") if mapped_hits else "No specific mapped v4 output was found",
            }
            rank = {"preserved": 2, "improved_and_replaced": 1, "missing": 0}
            if existing is None or rank[candidate["status"]] > rank[existing["status"]]:
                file_rows[legacy_file] = candidate
    for path in old_files:
        relative = path.relative_to(legacy_root).as_posix()
        if relative not in file_rows:
            file_rows[relative] = {
                "analysis_id": "unmapped", "output_family": "unmapped",
                "legacy_file": relative, "legacy_sha256": sha256_file(path),
                "new_file_count": 0, "new_files": "", "status": "missing",
                "reason": "The legacy file is not covered by an explicit output contract rule",
            }
    analysis_rows = []
    for analysis_id, state in sorted(analysis.items()):
        status = "missing" if state["missing"] else "preserved"
        analysis_rows.append({"analysis_id": analysis_id, **state, "status": status})
    rows = [file_rows[key] for key in sorted(file_rows)]
    inv = write_tsv(rows, output_dir / "legacy_output_inventory.tsv")
    write_tsv(family_rows, output_dir / "legacy_output_family_summary.tsv")
    matrix = write_tsv(analysis_rows, output_dir / "legacy_analysis_matrix.tsv")
    summary = {
        "legacy_root": str(legacy_root),
        "new_root": str(new_root) if new_root else None,
        "legacy_files": len(old_files),
        "new_files": len(new_files),
        "families": len(family_rows),
        "legacy_files_in_inventory": len(rows),
        "missing_files": sum(row["status"] == "missing" for row in rows),
        "missing_families": sum(row["status"] == "missing" for row in family_rows),
        "passed": all(row["status"] != "missing" for row in rows),
    }
    (output_dir / "legacy_inventory_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    return matrix, inv, summary


def _safe_spearman(x: pd.Series, y: pd.Series) -> float:
    return float(x.corr(y, method="spearman"))


def compare_de(legacy_de: Path, new_de: Path, output_dir: Path, padj: float = 0.05, abs_lfc: float = 1.0) -> tuple[Path, dict]:
    old = pd.read_csv(legacy_de)
    new = pd.read_csv(new_de, sep="\t")
    required_old = {"gene_id", "AveExpr", "stat", "adj.P.Val", "logFC"}
    required_new = {"gene_id", "baseMean", "stat", "padj", "log2FoldChange_ashr"}
    if missing := required_old - set(old):
        raise ValueError(f"Legacy DE table missing: {sorted(missing)}")
    if missing := required_new - set(new):
        raise ValueError(f"New DE table missing: {sorted(missing)}")
    joined = old[list(required_old)].merge(new[list(required_new)], on="gene_id", suffixes=("_old", "_new"))
    old_sig = (joined["adj.P.Val"] < padj) & (joined["logFC"].abs() > abs_lfc)
    new_sig = (joined["padj"] < padj) & (joined["log2FoldChange_ashr"].abs() > abs_lfc)
    away = joined["adj.P.Val"].isna() | (((joined["adj.P.Val"] - padj).abs() > 1e-6) & ((joined["logFC"].abs() - abs_lfc).abs() > 1e-3))
    metrics = [
        {"metric": "genes_compared", "value": len(joined), "threshold": None},
        {"metric": "baseMean_spearman", "value": _safe_spearman(joined["AveExpr"], joined["baseMean"]), "threshold": 0.999},
        {"metric": "wald_stat_spearman", "value": _safe_spearman(joined["stat_old"], joined["stat_new"]), "threshold": 0.999},
        {"metric": "shrunken_lfc_spearman", "value": _safe_spearman(joined["logFC"], joined["log2FoldChange_ashr"]), "threshold": 0.999},
        {"metric": "deg_classification_agreement", "value": float((old_sig[away] == new_sig[away]).mean()), "threshold": 0.99},
    ]
    for metric in metrics:
        metric["passed"] = metric["threshold"] is None or metric["value"] >= metric["threshold"]
    output_dir.mkdir(parents=True, exist_ok=True)
    metrics_path = write_tsv(metrics, output_dir / "de_regression_metrics.tsv")
    summary = {
        "legacy_de": str(legacy_de), "new_de": str(new_de), "padj": padj, "abs_lfc": abs_lfc,
        "old_significant": int(old_sig.sum()), "new_significant": int(new_sig.sum()),
        "passed": all(m["passed"] for m in metrics), "metrics": metrics,
    }
    (output_dir / "de_regression_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    joined.to_csv(output_dir / "de_regression_joined.tsv.gz", sep="\t", index=False, compression="gzip")
    return metrics_path, summary


def environment_resources() -> dict:
    mem_bytes = os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_AVPHYS_PAGES") if hasattr(os, "sysconf") else 0
    mem_gib = mem_bytes / 1024**3
    workers = max(1, min(3, int(mem_gib // 4)))
    return {"available_memory_gib": round(mem_gib, 2), "recommended_workers": workers, "threads_per_worker": 1}


def import_resource_pack(source_root: Path, destination_root: Path, registry_path: Path, species: str = "human", assembly: str = "hg38") -> list[dict]:
    """Copy explicit legacy resources into a frozen, checksummed pack.

    This is an import operation, not a network sync. Existing destination files
    are reused only when their checksums match.
    """
    source_root = source_root.resolve()
    destination_root.mkdir(parents=True, exist_ok=True)
    imported = []
    for source in all_files(source_root):
        relative = source.relative_to(source_root)
        if any(part.lower() == "logs" for part in relative.parts):
            continue
        destination = destination_root / relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        digest = sha256_file(source)
        if destination.exists() and sha256_file(destination) != digest:
            raise ValueError(f"Refusing to overwrite a different frozen resource: {destination}")
        if not destination.exists():
            shutil.copy2(source, destination)
        name = relative.as_posix()
        lower = name.lower()
        if "h.all.v2025.1" in lower:
            release, source_url, note = "MSigDB 2025.1", "https://www.gsea-msigdb.org/gsea/msigdb", "Imported from authorized historical project cache; MSigDB access terms apply"
        elif "encode_and_chea" in lower:
            release, source_url, note = "legacy frozen snapshot", "historical-project://Enrichr/ChEA-ENCODE-consensus", "Imported from historical project; upstream database terms apply"
        elif "encode_tf_chip" in lower:
            release, source_url, note = "ENCODE TF ChIP-seq 2015", "historical-project://Enrichr/ENCODE-2015", "Imported from historical project; upstream database terms apply"
        elif "collectri" in lower:
            release, source_url, note = "legacy frozen snapshot", "historical-project://OmniPath/CollecTRI", "Imported from historical project cache; OmniPath terms apply"
        elif "progeny" in lower:
            release, source_url, note = "legacy frozen top500 snapshot", "historical-project://OmniPath/PROGENy", "Imported from historical project cache; PROGENy/OmniPath terms apply"
        elif "string" in lower:
            release, source_url, note = "legacy project subnetwork", "historical-project://STRING", "Derived frozen project subnetwork; STRING terms apply"
        elif "pathview" in lower:
            release, source_url, note = "legacy frozen KEGG snapshot", "historical-project://KEGG/Pathview", "Imported for reproducibility; KEGG access and redistribution terms apply"
        elif "custom_gene_sets" in lower:
            release, source_url, note = "project-authored snapshot", "historical-project://custom-gene-sets", "User-provided project resource"
        else:
            release, source_url, note = "legacy frozen snapshot", "historical-project://legacy", "Imported from historical project; verify upstream license before redistribution"
        imported.append({
            "resource_id": re.sub(r"[^A-Za-z0-9_.-]+", "_", name),
            "species": species,
            "assembly": assembly,
            "release": release,
            "source_url": source_url,
            "downloaded_at": datetime.fromtimestamp(source.stat().st_mtime, tz=timezone.utc).isoformat(),
            "sha256": digest,
            "license_or_access_note": note,
            "local_path": str(destination.relative_to(registry_path.parent).as_posix()),
        })
    existing = load_yaml(registry_path) if registry_path.exists() else {"schema_version": 1, "resources": []}
    merged = {item["resource_id"]: item for item in (existing or {}).get("resources", [])}
    for item in imported:
        merged[item["resource_id"]] = item
    registry = {"schema_version": 1, "resources": sorted(merged.values(), key=lambda x: x["resource_id"])}
    registry_path.parent.mkdir(parents=True, exist_ok=True)
    registry_path.write_text(yaml.safe_dump(registry, sort_keys=False, allow_unicode=True), encoding="utf-8")
    validate_resources(registry_path)
    return imported


def import_counts(source: Path, project_dir: Path, destination_name: str | None = None) -> tuple[Path, dict]:
    """Copy a count matrix without modifying it and write a provenance/quality manifest."""
    project_dir = project_dir.resolve()
    # Import is intentionally allowed before the project becomes runnable.
    # The initialization workflow needs the immutable count hash in order to
    # finish project.yml and qc_approval.yml; requiring validate_project() here
    # creates a circular dependency for every new project.
    cfg = load_yaml(project_dir / "project.yml")
    schema = json.loads((ROOT / "internal/schemas/project.schema.json").read_text(encoding="utf-8"))
    errors = sorted(Draft202012Validator(schema).iter_errors(cfg), key=lambda e: list(e.path))
    if errors:
        detail = "\n".join(f"- {'.'.join(map(str, e.path)) or '<root>'}: {e.message}" for e in errors)
        raise ValueError(f"Invalid project configuration structure:\n{detail}")
    source = source.resolve()
    if not source.is_file():
        raise ValueError(f"Count matrix does not exist: {source}")
    input_dir = project_dir / "input"
    input_dir.mkdir(parents=True, exist_ok=True)
    destination = input_dir / (destination_name or source.name)
    digest = sha256_file(source)
    if destination.exists() and sha256_file(destination) != digest:
        raise ValueError(f"Refusing to overwrite a different imported count matrix: {destination}")
    if not destination.exists():
        shutil.copy2(source, destination)

    samples_path = project_dir / cfg.get("inputs", {}).get("samples_file", "samples.tsv")
    samples = pd.read_csv(samples_path, sep="\t", dtype=str, keep_default_na=False)
    if "sample" not in samples:
        raise ValueError("samples.tsv must contain a sample column")
    gene_id_column = cfg["inputs"]["gene_id_column"]
    counts = pd.read_csv(
        destination,
        sep="\t",
        usecols=lambda c: c == gene_id_column or c in set(samples["sample"]),
        low_memory=False,
    )
    if gene_id_column not in counts:
        raise ValueError(f"Count matrix is missing configured gene ID column: {gene_id_column}")
    missing = sorted(set(samples["sample"]) - set(counts.columns))
    if missing:
        raise ValueError("Count matrix is missing declared samples: " + ", ".join(missing))
    matrix = counts[list(samples["sample"])].apply(pd.to_numeric, errors="coerce")
    if matrix.isna().any().any() or (matrix < 0).any().any():
        raise ValueError("Count matrix contains missing, nonnumeric, or negative sample values")
    non_integer = int((matrix.to_numpy() % 1 != 0).sum())
    if non_integer:
        raise ValueError("Count matrix contains non-integer values; DESeq2 requires raw integer counts")
    provenance = cfg["counts_provenance"]
    if provenance["normalization"] in NORMALIZED_INPUTS:
        raise ValueError(f"DESeq2 import is forbidden for declared {provenance['normalization']} input")

    manifest = {
        "schema_version": 1,
        "source_path": str(source),
        "imported_path": str(destination.relative_to(project_dir).as_posix()),
        "imported_at": datetime.now(timezone.utc).isoformat(),
        "sha256": digest,
        "bytes": source.stat().st_size,
        "genes": int(len(matrix)),
        "samples": int(matrix.shape[1]),
        "all_values_integer": True,
        "all_values_nonnegative": True,
        "zero_fraction": float((matrix == 0).to_numpy().mean()),
        "library_sizes": {column: int(matrix[column].sum()) for column in matrix},
        "provenance": provenance,
        "limitation": "Integer values are necessary but do not prove raw-count provenance; source_method and annotation remain user-declared.",
    }
    (input_dir / "source_manifest.yml").write_text(
        yaml.safe_dump(manifest, sort_keys=False, allow_unicode=True), encoding="utf-8"
    )
    return destination, manifest


def sync_declared_resources(manifest_path: Path, destination_root: Path, registry_path: Path) -> list[dict]:
    """The sole network-enabled resource synchronization operation."""
    manifest = load_yaml(manifest_path)
    destination_root.mkdir(parents=True, exist_ok=True)
    synced = []
    for item in manifest.get("resources", []):
        destination = destination_root / item["filename"]
        destination.parent.mkdir(parents=True, exist_ok=True)
        if not destination.exists():
            partial = destination.with_suffix(destination.suffix + ".partial")
            aria2 = shutil.which("aria2c")
            if aria2:
                subprocess.run([
                    aria2, "--continue=true", "--max-connection-per-server=16", "--split=16",
                    "--min-split-size=4M", "--file-allocation=none", "--auto-file-renaming=false",
                    "--allow-overwrite=false", "--summary-interval=30", "--dir", str(destination.parent),
                    "--out", partial.name, item["source_url"],
                ], check=True)
            else:
                offset = partial.stat().st_size if partial.exists() else 0
                headers = {"User-Agent": "bulk-rnaseq-v4-resource-sync/4.0"}
                if offset:
                    headers["Range"] = f"bytes={offset}-"
                request = urllib.request.Request(item["source_url"], headers=headers)
                with urllib.request.urlopen(request, timeout=120) as response, partial.open("ab" if offset and response.status == 206 else "wb") as handle:
                    shutil.copyfileobj(response, handle, length=1024 * 1024)
            partial.replace(destination)
        digest = sha256_file(destination)
        expected = item.get("sha256")
        if expected and digest.lower() != expected.lower():
            raise ValueError(f"SHA256 mismatch for synced resource: {destination}")
        synced.append({
            "resource_id": item["resource_id"], "species": item["species"], "assembly": item["assembly"],
            "release": item["release"], "source_url": item["source_url"],
            "downloaded_at": datetime.now(timezone.utc).isoformat(), "sha256": digest,
            "license_or_access_note": item["license_or_access_note"],
            "local_path": destination.relative_to(registry_path.parent).as_posix(),
        })
    existing = load_yaml(registry_path) if registry_path.exists() else {"resources": []}
    merged = {item["resource_id"]: item for item in (existing or {}).get("resources", [])}
    merged.update({item["resource_id"]: item for item in synced})
    registry_path.write_text(yaml.safe_dump({"schema_version": 1, "resources": sorted(merged.values(), key=lambda x: x["resource_id"])}, sort_keys=False, allow_unicode=True), encoding="utf-8")
    validate_resources(registry_path)
    return synced


def clean_preview(project_id: str, scopes: list[str], root: Path = ROOT) -> list[dict]:
    project = root / "projects" / project_id
    allowed = {
        "staging": project / "work" / "staging",
        "failed": project / "work" / "history" / "failed",
        "cache": project / "work" / "cache",
    }
    unknown = sorted(set(scopes) - set(allowed))
    if unknown:
        raise ValueError("Unknown clean scopes: " + ", ".join(unknown))
    rows = []
    for scope in scopes:
        target = allowed[scope].resolve()
        if root.resolve() not in target.parents:
            raise ValueError(f"Clean target escapes pipeline root: {target}")
        files = all_files(target) if target.exists() else []
        rows.append({"scope": scope, "target": str(target), "exists": target.exists(), "files": len(files), "bytes": sum(p.stat().st_size for p in files)})
    return rows
