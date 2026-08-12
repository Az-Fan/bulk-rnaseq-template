# Optional upstream adapter

The repository pins `xuzhougeng/omics-pipelines` as a Git submodule. The
adapter exposes only its RNA-seq raw-count target (`results/04-quant/counts.tsv`).
It never calls the upstream DESeq2 or enrichment rules, and it never builds a
STAR index.

Upstream execution is blocked unless `project.yml` records a confirmed
`enabled` decision and the operator supplies the exact token printed by:

```bash
pixi run pipeline upstream-plan --project projects/<project_id>
```

The default template sets `execution_host: server` and `allow_wsl: false`.
Clone on a server with submodules:

```bash
git clone --recurse-submodules <repository-url>
```

The server must already contain the STAR index declared in the project-local
upstream config. FASTQ files and reference indices remain outside Git.
