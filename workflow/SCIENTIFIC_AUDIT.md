# Scientific audit and support matrix

Audit date: 2026-08-14

## Overall conclusion

The validated human bulk RNA-seq core (raw-count QC, DESeq2, ORA/GSEA, frozen-resource activity analysis, promoter motif, and declared exploratory views) follows standard practice when its run plan is explicitly confirmed. It must not be described as universally complete: peak-aware motif, WGCNA, generic Pathview, personalized executors, and non-human resource-dependent modules remain hard-blocked until implemented and validated.

## Module assessment

| Module | Status | Scientific contract | Remaining limitation |
|---|---|---|---|
| Raw-count and sample QC | validated core | finite nonnegative integers; explicit filtering; blind VST PCA; no automatic exclusions | upstream provenance can remain unknown but is reported as a limitation |
| Differential expression | validated core | one full-rank DESeq2 model; residual degrees of freedom; ≥2 biological samples per contrasted level; raw/Wald/BH plus ashr LFC | complex interaction contrasts still require deliberate contrast review |
| ORA | validated human/frozen resources | explicit DEG profiles; tested mapped-gene universe; Up/Down separated; BH; full tables retained | current database pack is frozen human MSigDB 7.5.1 for migration parity |
| GSEA | validated human/frozen resources | complete finite ranking; explicit databases/set sizes/FDR/seed; signed NES | gene-set release is project resource, never an implicit latest release |
| GSVA/TF/PROGENy | validated human/frozen resources | VST expression or complete signed statistic; frozen CollecTRI/PROGENy; declared target TF/pathway must exist | activity is associative, not direct binding or causality |
| PPI | conditional human | explicit DEG profile/top-N/STRING score; weighted centrality and Louvain; fixed layout seed | result depends strongly on frozen STRING subnetwork and selection settings |
| Promoter motif | conditional human | explicit DEG profile/window/minimum foreground/STREME settings; matched tested non-DE background | counts-only motif is exploratory; it is not direct binding evidence |
| Peak-aware motif | blocked | none | no validated generic executor yet |
| Target-gene boxplots | validated display | declared genes; unchanged blind=FALSE VST values; all samples shown | no inferential test is added automatically |
| Custom gene sets | conditional | explicit rank/min size/FDR/seed/curve mode and heatmap clustering | legacy projects may deliberately use legacy ashr-LFC ranking |
| Pathview | legacy-only | frozen historical map may be copied and checksum-verified | generic current-result Pathview is not implemented |
| WGCNA | blocked | none | no validated executor; sample-size suitability alone is insufficient |
| Personalized analysis | blocked by default | must consume standard outputs only | requires a separately reviewed executor |

## Safety improvements from this audit

1. Added project-level confirmation for eight decision groups and a fresh per-run plan token.
2. Bound the token to project config, samples, contrasts, QC approval, and counts source manifest.
3. Made database selection, ORA/GSEA reporting thresholds, seed and redundancy settings explicit.
4. Made DE primary profile, volcano labels and heatmap gene counts explicit.
5. Made PPI input profile/size, STRING score, module FDR and label count explicit.
6. Made promoter window, motif foreground, STREME widths/P threshold and motif display count explicit.
7. Made target TF/pathway existence a failure condition instead of silently substituting another pathway.
8. Made custom gene-set rank, set sizes, seed, FDR, curve mode and heatmap clustering explicit.
9. Added residual-degrees-of-freedom and minimum biological-replication checks.
10. Added DESeq2 size-factor, dispersion, P-value and Cook's-distance diagnostics without changing results.
11. Added hard failures for modules that previously could appear enabled without a complete executor.
12. Repaired source/document encoding that could render corrupt PDF/HTML labels.

## Decisions that are locked methods

The following are not ordinary cosmetic preferences: raw counts for DESeq2, design covariates instead of corrected-expression testing, tested mapped genes as ORA universe, complete ranks for GSEA, BH adjustment, one shared DESeq2 model, and explicit contrast direction. Changing any of them requires a new method version and regression evidence.

## Decisions that the user must make

The user chooses the experimental design and exclusions; thresholds and primary profiles; enabled modules; ORA/GSEA databases and reporting windows; target genes/TFs/pathways; network and motif parameters; custom gene sets; and report/export scope. `pixi run plan` is the authoritative current inventory.
