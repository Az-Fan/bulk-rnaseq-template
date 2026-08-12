#!/usr/bin/env bash
set -euo pipefail
if [[ $# -ne 1 ]]; then echo "Usage: pixi run run-motif -- <completed_run_dir>" >&2; exit 2; fi
root="$(pwd -P)"
run_dir="$(realpath "$1")"
resource_dir="$root/resources/data/shared/hg38-motif"
fasta="$resource_dir/GRCh38.primary_assembly.genome.fa.gz"
gtf="$resource_dir/gencode.v49.primary_assembly.annotation.gtf.gz"
jaspar="$resource_dir/JASPAR2026_CORE_vertebrates_non-redundant_pfms_meme.txt"
de="$run_dir/02_Differential/Treatment_vs_Control/Tables/Differential_Full.tsv"
for file in "$fasta" "$gtf" "$jaspar" "$de"; do [[ -s "$file" ]] || { echo "Missing frozen motif input: $file" >&2; exit 3; }; done
out="$run_dir/08_Motif/Promoter"
resume_partial=0
if [[ -e "$out" ]]; then
  if [[ -s "$out/STREME/streme.txt" && -s "$out/STREME/sites.tsv" && -s "$out/AME_JASPAR/ame.tsv" && -s "$out/Tables/Motif_Run_Summary.tsv" && ! -e "$out/motif_status.tsv" ]]; then
    cp "$out/STREME/streme.txt" "$out/Tables/STREME_De_Novo_Motifs.txt"
    cp "$out/STREME/sites.tsv" "$out/Tables/STREME_De_Novo_Sites.tsv"
    cp "$out/AME_JASPAR/ame.tsv" "$out/Tables/AME_JASPAR_Known_Motifs.tsv"
    python internal/reporting/render_motif_figures.py "$run_dir"
    printf 'motif_status\tcomplete\n' > "$out/motif_status.tsv"
    python internal/motif/finalize_motif.py "$run_dir/run_manifest.json"
    exit 0
  fi
  if [[ -s "$out/STREME/streme.xml" && ! -s "$out/AME_JASPAR/ame.tsv" && ! -e "$out/motif_status.tsv" ]]; then
    echo "Resuming validated motif state: STREME complete; AME publication incomplete" >&2
    resume_partial=1
  else
    echo "Motif output exists but is not a validated resumable state: $out" >&2
    exit 4
  fi
fi
mkdir -p "$out/Tables" "$out/Figures"
tmp="$(mktemp -d)"; trap 'rm -rf -- "$tmp"' EXIT
gzip -dc "$fasta" > "$tmp/genome.fa"
samtools faidx "$tmp/genome.fa"
awk -F '\t' 'BEGIN{OFS="\t"} $3=="gene" {id=""; if(match($9,/gene_id "([^"]+)"/,a)) id=a[1]; sub(/\.[0-9]+$/, "", id); if($7=="+") {s=$4-1001; if(s<0)s=0; e=$4+100} else {s=$5-100; if(s<0)s=0; e=$5+1001} print $1,s,e,id,".",$7}' <(gzip -dc "$gtf") > "$tmp/all_promoters.bed"
awk -F '\t' 'NR==FNR {if(NR>1) {tested[$1]=1; if($10=="TRUE") sig[$1]=1}; next} ($4 in sig)' "$de" "$tmp/all_promoters.bed" > "$tmp/foreground.bed"
awk -F '\t' 'NR==FNR {if(NR>1) {tested[$1]=1; if($10=="TRUE") sig[$1]=1}; next} (($4 in tested) && !($4 in sig))' "$de" "$tmp/all_promoters.bed" | awk 'NR%10==1' > "$tmp/background.bed"
bedtools getfasta -s -name -fi "$tmp/genome.fa" -bed "$tmp/foreground.bed" -fo "$tmp/foreground.fa"
bedtools getfasta -s -name -fi "$tmp/genome.fa" -bed "$tmp/background.bed" -fo "$tmp/background.fa"
fg_n="$(grep -c '^>' "$tmp/foreground.fa" || true)"; bg_n="$(grep -c '^>' "$tmp/background.fa" || true)"
[[ "$fg_n" -ge 10 && "$bg_n" -ge 10 ]] || { echo "Too few promoters for motif analysis: foreground=$fg_n background=$bg_n" >&2; exit 5; }
if [[ "$resume_partial" -eq 0 ]]; then
  streme --p "$tmp/foreground.fa" --n "$tmp/background.fa" --oc "$out/STREME" --dna --minw 6 --maxw 15 --thresh 0.05
fi
mkdir -p "$out/AME_JASPAR"
# Keep AME's statistics in deterministic text mode; the publication layer below
# supplies portable HTML and independent high-resolution figures.
ame --text --control "$tmp/background.fa" "$tmp/foreground.fa" "$jaspar" > "$out/AME_JASPAR/ame.tsv"
printf 'metric\tvalue\nforeground_promoters\t%s\nbackground_promoters\t%s\npromoter_definition\t-1000_to_+100_bp_from_TSS\nassembly\thg38\nannotation\tGENCODE_v49\nknown_motifs\tJASPAR_2026_CORE_vertebrates\ninterpretation\tMotif_enrichment_is_not_direct_binding_or_causality\n' "$fg_n" "$bg_n" > "$out/Tables/Motif_Run_Summary.tsv"
cp "$out/STREME/streme.txt" "$out/Tables/STREME_De_Novo_Motifs.txt"
cp "$out/STREME/sites.tsv" "$out/Tables/STREME_De_Novo_Sites.tsv"
cp "$out/AME_JASPAR/ame.tsv" "$out/Tables/AME_JASPAR_Known_Motifs.tsv"
python internal/reporting/render_motif_figures.py "$run_dir"
printf 'motif_status\tcomplete\n' > "$out/motif_status.tsv"
python internal/motif/finalize_motif.py "$run_dir/run_manifest.json"
