"""RNA-seq rendering and numeric invariance, with generic synthetic inputs."""
from pathlib import Path
import json
import subprocess

import numpy as np
import pandas as pd
from jsonschema import Draft202012Validator
import yaml

ROOT = Path(__file__).resolve().parents[1]
R = ROOT / '.pixi/envs/default/bin/Rscript'


def run_r(tmp_path, body):
    script = tmp_path / 'fixture.R'
    script.write_text(f"source('{ROOT}/workflow/functions.R')\nsource('{ROOT}/workflow/enrichment_visualization.R')\n" + body)
    result = subprocess.run([str(R), str(script)], cwd=ROOT, capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr


CONTEXT = '''
ids <- paste0('ENSG', seq_len(65))
symbols <- c(paste0('Gene', seq_len(64)), 'Gene1')
samples <- paste0('S', seq_len(7))
meta <- data.frame(condition=c('Base','Base','Perturbed','Perturbed','Other','Other','Excluded'), row.names=samples)
values <- outer(seq_len(65), seq_len(7), function(g,s) 5+g/20+s/4+sin(g*s))
values[1,] <- 7
rownames(values) <- ids; colnames(values) <- samples
# A duplicate symbol with a larger effect must not replace its first input row.
de <- data.frame(gene_id=ids, gene_symbol=symbols, entrez_id=as.character(seq_len(65)),
                 stat=c(seq(5,-5,length.out=64),100), log2FoldChange_ashr=seq(3,-3,length.out=65),
                 log2FoldChange_raw=seq(3,-3,length.out=65), pvalue=.4, padj=.8,
                 direction_Primary='Not_significant')
de <- de[order(-de$stat),]
context <- new.env()
context$root <- ROOT; context$run_dir <- OUT
context$config <- list(species='human', enrichment=list(ora_fdr=.1,gsea_reporting_fdr=.2,
    overview_terms=18L,ora_overview_terms=8L,gsea_curves=2L,redundancy_jaccard=.7),
    thresholds=list(padj=.05,abs_lfc=1), analysis=list(modules=list(enrichment=list(status='enabled',confirmed=TRUE))),
    resources=list(enrichment_gene_sets='fixture'))
context$contrasts <- data.frame(contrast_id='Perturbed_vs_Base', factor='condition', numerator='Perturbed', denominator='Base')
context$state <- list(vst=values,selected_samples=samples[1:6],samples=meta,dds=values,
                     display_factor='condition',annotation=data.frame(gene_id=ids,gene_symbol=symbols,row.names=ids))
context$module_status <- list()
'''


def test_compact_figures_use_exact_gene_effects_and_contrast_vst(tmp_path):
    body = f'ROOT <- {json.dumps(str(ROOT))}\nOUT <- {json.dumps(str(tmp_path))}\n' + CONTEXT
    body += '''
ranks <- setNames(seq(5,-5,length.out=64), symbols[1:64])
paths <- c('HALLMARK_VERY_LONG_PATHWAY_NAME_WITH_TRANSCRIPTIONAL_REGULATION_AND_CELLULAR_RESPONSE_TO_EXTERNAL_STIMULUS',
           'REACTOME_VERY_LONG_PATHWAY_NAME_WITH_TRANSCRIPTIONAL_REGULATION_AND_CELLULAR_RESPONSE_TO_EXTERNAL_STIMULUS')
sets <- setNames(list(symbols[1:35],symbols[40:64]),paths)
gsea <- data.frame(pathway=paths,NES=c(1.7,-1.9),padj=c(.0003,.12),
                   leadingEdge=c(paste(symbols[1:35],collapse='/'),paste(symbols[40:64],collapse='/')))
write_compact_enrichment(context,'Perturbed_vs_Base',de,gsea,'Hallmark','GSEA',ranks=ranks,pathways=sets)
ora <- data.frame(pathway=c(paths[1],paths[1]),direction=c('Up','Down'),padj=c(.04,.08),
                  overlap_genes=c('1/2/3/4/5','64/65'),gene_ratio=c(.5,.25),overlap_count=c(5,2))
write_compact_enrichment(context,'Perturbed_vs_Base',de,ora,'Hallmark','ORA','Primary')
write_compact_enrichment(context,'Perturbed_vs_Base',de,data.frame(),'Empty','GSEA',ranks=ranks,pathways=sets)
stopifnot(enrichment_plot_modes(list())$compact,!enrichment_plot_modes(list())$legacy,
  enrichment_plot_modes(list(migration=list(requires_legacy_parity=TRUE)))$legacy,
  !enrichment_plot_modes(list(enrichment=list(visualization='legacy')))$compact)
'''
    run_r(tmp_path, body)
    root = tmp_path / '03_Enrichment/Perturbed_vs_Base'
    display = root / 'GSEA/Tables/Display/GSEA_Hallmark'
    chosen = pd.read_csv(display / 'Identifier_Representatives_Selected.tsv', sep='\t')
    assert chosen.set_index('gene_key').loc['Gene1', 'gene_id'] == 'ENSG1'
    jobs = json.loads((display / 'Visualization_Jobs.json').read_text())
    for job in jobs:
        for field in {'data', 'pdf', 'ranks', 'members'} & job.keys():
            assert not Path(job[field]).is_absolute()
            assert (tmp_path / job[field]).is_file()
    combo = pd.read_csv(display / 'GSEA_Hallmark_01_Core_Combo_Display.tsv', sep='\t')
    assert set(combo['sample']) == {'S1', 'S2', 'S3', 'S4'}
    assert set(combo['group']) == {'Base', 'Perturbed'}
    assert combo.gene_key.nunique() == 30
    assert len(pd.read_csv(display / 'GSEA_Hallmark_01_Core_Combo_All_Core_Genes.tsv', sep='\t')) == 35
    assert combo[combo.gene_key.eq('Gene1')].zscore.isna().all()
    for _, group in combo[~combo.gene_key.eq('Gene1')].groupby('gene_key'):
        expected = (group.vst - group.vst.mean()) / group.vst.std(ddof=1)
        np.testing.assert_allclose(group.zscore, expected, atol=1e-12)
        gene = int(group.gene_key.iloc[0][4:])
        np.testing.assert_allclose(group.log2FC, np.linspace(3, -3, 65)[gene - 1], atol=1e-12)
        np.testing.assert_allclose(group.gene_padj, .8, atol=1e-12)
    matrix = pd.read_csv(root / 'ORA/Tables/Display/ORA_Hallmark_Primary/Pathway_Gene_All_Memberships.tsv', sep='\t')
    assert matrix.term_id.nunique() == 2  # same pathway, distinct Up/Down foregrounds
    assert matrix[~matrix.is_member].log2FC.isna().all()
    assert matrix.groupby('term_id').is_member.sum().tolist() == [2, 5]
    assert pd.read_csv(root / 'GSEA/Tables/Display/GSEA_Empty/Display_Status.tsv', sep='\t').status.iloc[0] == 'no_reportable_pathways'
    pdfs = list(root.rglob('*.pdf'))
    assert len(pdfs) == 6 and not list(root.rglob('*.png'))
    for pdf in pdfs:
        assert pdf.read_bytes().startswith(b'%PDF')
        count = subprocess.run(['gs', '-q', '-dNODISPLAY', '-dNOSAFER', '-c',
                                f'({pdf}) (r) file runpdfbegin pdfpagecount = quit'],
                               capture_output=True, text=True, check=True)
        assert count.stdout.strip() == '1'
    # Render every figure; inspect resulting examples separately during QA.
    qa = tmp_path / 'rendered'; qa.mkdir()
    for i, pdf in enumerate(pdfs):
        subprocess.run(['gs','-q','-dBATCH','-dNOPAUSE','-sDEVICE=pngalpha','-r90',
                        f'-sOutputFile={qa / (str(i)+".png")}',str(pdf)], check=True, capture_output=True)
    (tmp_path / 'pdf_inventory.json').write_text(json.dumps([str(p) for p in pdfs], indent=2))


def test_actual_enrichment_statistics_do_not_depend_on_display_mode(tmp_path):
    body = f'ROOT <- {json.dumps(str(ROOT))}\nOUT <- {json.dumps(str(tmp_path))}\n' + CONTEXT
    body += f"source('{ROOT}/workflow/03_enrichment.R')\n"
    body += '''
gene_sets <- data.frame(gs_cat='H',gs_subcat='',gs_name=rep(c('HALLMARK_POSITIVE','HALLMARK_NEGATIVE'),each=15),
                       gene_symbol=c(symbols[1:15],symbols[50:64]),entrez_gene=as.character(c(1:15,50:64)))
resource <- file.path(OUT,'sets.rds')
saveRDS(list(species='Homo sapiens',gene_sets=gene_sets,msigdb_release='7.5.1'),resource)
resource_path <- function(context,id) resource
context$state$differential <- list(Perturbed_vs_Base=list(result=de))
context$config$enrichment$ora_databases <- 'Hallmark'
context$config$enrichment$gsea_databases <- 'Hallmark'
context$config$enrichment$min_size <- 10L; context$config$enrichment$max_size <- 50L
context$config$enrichment$gsea_curves <- 0L
# Legacy plots need not be rendered twice; the renderer test exercises all new families.
save_publication_figure <- function(...) invisible(NULL)
save_publication_grob <- function(...) invisible(NULL)
for (mode in c('compact','legacy')) {
  context$run_dir <- file.path(OUT,mode);dir.create(context$run_dir)
  context$config$enrichment$visualization <- mode
  run_enrichment(context)
}
'''
    run_r(tmp_path, body)
    def scientific_files(mode):
        root = tmp_path / mode
        return {str(p.relative_to(root)): p.read_bytes() for p in root.rglob('*.tsv*')
                if '/Display/' not in str(p) and any(key in p.name for key in ('Full','Universe','Ranked','Provenance'))}
    compact, legacy = scientific_files('compact'), scientific_files('legacy')
    assert compact and compact == legacy
    table = pd.read_csv(tmp_path / 'compact/03_Enrichment/Perturbed_vs_Base/GSEA/Tables/GSEA_Full.tsv', sep='\t')
    assert len(table) == 2 and table.NES.min() < 0 < table.NES.max()


def test_visualization_config_is_valid_and_never_alters_legacy_gates():
    schema = json.loads((ROOT / 'internal/schemas/project.schema.json').read_text())
    config = yaml.safe_load((ROOT / 'projects/_template/project.yml').read_text())
    assert not list(Draft202012Validator(schema).iter_errors(config))
    assert config['enrichment']['visualization'] == 'compact'
    config['enrichment']['visualization'] = 'unsupported'
    assert list(Draft202012Validator(schema).iter_errors(config))
