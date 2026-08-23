import hail as hl
import pandas as pd
import argparse
import hailtop.batch as hb
import hailtop.fs as hfs

# from aou_gwas import *  # for dataproc
from utils.utils import *  # for QoB
from utils.resources import *  # for QoB
from utils.results_loading import * # for QoB

ANCESTRIES = ['AFR', 'AMR', 'EAS', 'EUR', 'MID', 'SAS', 'ALL']

def get_AF_bins(field):
    bins = (hl.case()
            .when(field < 0.0001, '[0, 0.0001)')
            .when(field < 0.001, '[0.0001, 0.001)')
            .when(field < 0.01, '[0.001, 0.01)')
            .when(field < 0.1, '[0.01, 0.1)')
            .when(field >= 0.1, '[0.1, )')
            .default(hl.missing(hl.tstr)))

    return bins

def get_coverage_bins(field):
    bins = (hl.case()
            .when(field < 10, '[0, 10)')
            .when(field < 20, '[10, 20)')
            .when(field < 30, '[20, 30)')
            .when(field >= 30, '[30, )')
            .default(hl.missing(hl.tstr)))
    return bins

def add_liftover_rg37_to_rg38_ht(ht: hl.Table):
    """
    Add liftover to Hail Table from rg37 to rg38
    :param Table ht: Hail Table to add liftover on
    :return: Hail Table
    :rtype: hl.Table
    """
    rg37 = hl.get_reference("GRCh37")
    rg38 = hl.get_reference("GRCh38")
    if not rg37.has_liftover("GRCh38"):
        rg37.add_liftover(
            "gs://hail-common/references/grch37_to_grch38.over.chain.gz", rg38
        )
    ht = ht.annotate(new_locus=hl.liftover(ht.locus, "GRCh38"))
    ht = ht.filter(hl.is_defined(ht.new_locus))
    ht = ht.key_by(locus=ht.new_locus, alleles=ht.alleles)
    return ht

def annotate_clinvar_pathogenicity_ht():
    """
    Re-define pathogenicity group for clinvar variants
    :return: Hail Table with new definitions
    :rtype: Table
    """
    from gnomad.resources.grch38.reference_data import clinvar
    clinvar_ht = clinvar.ht()
    clinvar_ht = clinvar_ht.filter(
        ~clinvar_ht.info.CLNREVSTAT[0].contains("no_assertion")
    )  # remove clinvar zero star variants
    clinvar_ht = clinvar_ht.explode(clinvar_ht.info.CLNSIG)
    clinvar_ht = clinvar_ht.annotate(
        pathogenicity=hl.case()
        .when(
            hl.literal(
                {"Pathogenic", "Likely_pathogenic", "Pathogenic/Likely_pathogenic"}
            ).contains(clinvar_ht.info.CLNSIG),
            "P/LP",
        )
        .when(
            hl.literal({"Uncertain_significance"}).contains(clinvar_ht.info.CLNSIG),
            "VUS",
        )
        .when(
            hl.literal({"Benign", "Likely_benign", "Benign/Likely_benign"}).contains(
                clinvar_ht.info.CLNSIG
            ),
            "B/LB",
        )
        .or_missing()
    )
    return clinvar_ht

def main(args):
    analysis_type = "variant" if args.single_variant_only else "gene"
    result_type = "gene" if args.load_gene_results else "variant"
    name = 'ACAF' if analysis_type == 'variant' else 'exome'
    hl.init(
        tmp_dir=TMP_BUCKET,
        driver_memory="highmem",
        driver_cores=8,
        worker_memory="highmem",
        worker_cores=1,
        default_reference="GRCh38",
        log="/aou_saige.log",
        app_name=f'aou_QC_SAIGE_{analysis_type}_results'
    )
    print(f'Analysis type: {analysis_type}')
    print(f'Result type: {result_type}')
    def get_ancestry_pruned_call_stats_ht_path(ancestry, name, pruned):
        pruned = 'pruned' if pruned else 'pre_pruning'
        return f'{DATA_PATH}/utils/call_stats/{name}_{pruned}/{ancestry.upper()}_{name}_call_stats.ht'
    if args.create_call_stats_ht:
        path = get_ancestry_pruned_call_stats_ht_path(ancestry='full', name=name, pruned=True)
        if not hfs.exists(f'{path}/_SUCCESS') or args.overwrite:
            print(f'[{name}] Exporting FULL call stats to path: {path}')
            vep_ht = hl.read_table(f'{DATA_PATH}/vep/aou_v8_vep_full_raw.ht')
            global_ht = hl.read_table(get_ancestry_pruned_call_stats_ht_path(ancestry='all', name=name, pruned=True))
            ht = global_ht.annotate(**vep_ht[global_ht.key])
            ht = ht.drop('call_stats')
            ht = ht.annotate(
                freq = hl.struct(**{f'{anc.upper()}': hl.read_table(get_ancestry_pruned_call_stats_ht_path(ancestry=anc, name=name, pruned=True))[ht.key].call_stats for anc in ANCESTRIES}),
                freq_raw = hl.struct(**{f'{anc.upper()}': hl.read_table(get_ancestry_pruned_call_stats_ht_path(ancestry=anc, name=name, pruned=False))[ht.key].call_stats for anc in ANCESTRIES})
            )
            ht.describe()
            ht.naive_coalesce(1000).checkpoint(path, overwrite=args.overwrite)
        ht = hl.read_table(path)
        ht.describe()
        ht.show()
        print(ht.count())

    if args.export_variant_qc_info:
        print(f'[{name}] Exporting variant QC info to path: {path}')
        path = f'{DATA_PATH}/utils/raw_mt/ACAF/ALL_ACAF.mt' if analysis_type == 'variant' else f'{DATA_PATH}/utils/raw_mt/Exome/ALL_Exome.mt'
        if not hl.hadoop_exists(f'{DATA_PATH}/qc/aou_{name}_rows_full.ht/_SUCCESS'):
            mt = hl.read_matrix_table(path)
            mt.describe()
            ht = mt.rows().checkpoint(f'{DATA_PATH}/qc/aou_{name}_rows_full.ht', _read_if_exists = True)
        ht = hl.read_table(f'{DATA_PATH}/qc/aou_{name}_rows_full.ht')
        ht = ht.select(approx_coverage = ht.variant_qc.gq_stats.mean/3)
        call_stats_ht = hl.read_table(get_ancestry_pruned_call_stats_ht_path(ancestry='full', name=name, pruned=True))
        call_stats_ht = call_stats_ht.filter(hl.is_defined(call_stats_ht.freq))
        ht = call_stats_ht.annotate(**ht[call_stats_ht.key])
        ht = ht.naive_coalesce(1000).checkpoint(f'{DATA_PATH}/qc/aou_{name}_variant_qc.ht', overwrite = args.overwrite, _read_if_exists = not args.overwrite)
        ht.show()

    if args.export_gene_qc_info:
        print(f'[{name}] Exporting gene QC info to path: {path}')
        ht = hl.read_table(f'{DATA_PATH}/qc/aou_exome_variant_qc.ht')
        ht = ht.annotate(vep=ht.vep.annotate(transcript_consequences=ht.vep.transcript_consequences.filter(lambda x: x.gene_id.startswith('ENSG'))))
        from gnomad.utils.vep import process_consequences
        ht = process_consequences(ht)
        ht = ht.explode(ht.vep.worst_csq_by_gene_canonical)
        ht = ht.filter(ht.vep.worst_csq_by_gene_canonical.gene_id.startswith('ENSG'))
        ht = ht.annotate(annotation=annotation_case_builder(ht.vep.worst_csq_by_gene_canonical),
                         gene_id=ht.vep.worst_csq_by_gene_canonical.gene_id,
                         gene_symbol=ht.vep.worst_csq_by_gene_canonical.gene_symbol) 
        ht = ht.annotate(annotation = hl.if_else(hl.literal({'missense', 'LC'}).contains(ht.annotation), 'missenseLC', ht.annotation))
        ht = ht.checkpoint(f'{DATA_PATH}/qc/aou_exome_variant_qc_annotated.ht', _read_if_exists=not args.overwrite, overwrite=args.overwrite)
        combined_ht = ht.filter(hl.literal(['missenseLC', 'pLoF']).contains(ht.annotation))
        combined_ht = combined_ht.group_by('gene_id', 'gene_symbol').aggregate(
            mean_coverage=hl.agg.mean(combined_ht.approx_coverage),
            CAF_global_raw=hl.agg.sum(combined_ht.freq_raw.ALL.AF[1]),
            **{
                f"mean_coverage_{max_MAF}": hl.agg.filter(combined_ht.freq.ALL.AF[1] < max_MAF,
                                                          hl.agg.mean(combined_ht.approx_coverage))
                for max_MAF in [0.01, 0.001, 0.0001]
            },
            **{
                f"CAF_global_{max_MAF}": hl.agg.filter(combined_ht.freq.ALL.AF[1] < max_MAF,
                                                       hl.agg.sum(combined_ht.freq.ALL.AF[1]))
                for max_MAF in [0.01, 0.001, 0.0001]
            },
            **{
                f"CAF_{anc}_raw": hl.agg.sum(combined_ht.freq_raw[anc.upper()].AF[1])
                for anc in ANCESTRIES if anc != 'ALL'
            },
            **{
                f"CAF_{anc}_{max_MAF}": hl.agg.filter(combined_ht.freq[anc.upper()].AF[1] < max_MAF,
                                                      hl.agg.sum(combined_ht.freq[anc.upper()].AF[1]))
                for anc in ANCESTRIES if anc != 'ALL'
                for max_MAF in [0.01, 0.001, 0.0001]
            },
        )
        combined_ht = combined_ht.annotate(annotation = 'pLoF;missenseLC').key_by('gene_id', 'gene_symbol','annotation')
        fields = ['gene_id', 'gene_symbol','annotation'] + list(combined_ht.row_value)
        combined_ht = combined_ht.key_by().select(*fields).key_by('gene_id', 'gene_symbol','annotation')
        gene_ht = ht.group_by('gene_id', 'gene_symbol','annotation').aggregate(
            mean_coverage = hl.agg.mean(ht.approx_coverage),
            CAF_global_raw = hl.agg.sum(ht.freq_raw.ALL.AF[1]),
            **{
                f"mean_coverage_{max_MAF}": hl.agg.filter(ht.freq_raw.ALL.AF[1] < max_MAF,
                                                       hl.agg.mean(ht.approx_coverage))
                for max_MAF in [0.01, 0.001, 0.0001]
            },
            **{
                f"CAF_global_{max_MAF}": hl.agg.filter(ht.freq.ALL.AF[1] < max_MAF,
                                                       hl.agg.sum(ht.freq.ALL.AF[1]))
                for max_MAF in [0.01, 0.001, 0.0001]
            },
            **{
                  f"CAF_{anc}_raw": hl.agg.sum(ht.freq_raw[anc.upper()].AF[1])
                  for anc in ANCESTRIES if anc != 'ALL'
              },
            **{
                f"CAF_{anc}_{max_MAF}": hl.agg.filter(ht.freq[anc.upper()].AF[1] < max_MAF ,hl.agg.sum(ht.freq[anc.upper()].AF[1]))
                for anc in ANCESTRIES if anc != 'ALL'
                for max_MAF in [0.01, 0.001, 0.0001]
            },
        )
        gene_ht = gene_ht.select(*list(combined_ht.row_value))
        gene_ht = gene_ht.union(combined_ht)
        gene_ht = gene_ht.checkpoint(f'{DATA_PATH}/qc/aou_exome_gene_qc.ht', overwrite=args.overwrite, _read_if_exists=not args.overwrite)
        gene_ht.describe()
        gene_ht.show()
        print(gene_ht.count())
        if not hl.hadoop_exists(f'{DATA_PATH}/qc/aou_exome_gene_qc.txt.bgz') or args.overwrite:
            gene_ht.export(f'{DATA_PATH}/qc/aou_exome_gene_qc.txt.bgz')

    if args.count_hq_genes:
        phenotype_outliers = {
            'eur':['PP_928.1', 'PP_932.11'],
        }
        for ancestry in ['afr', 'amr', 'eas', 'eur', 'mid', 'sas', 'meta']:
            for maxMAF in [0.01, 0.001, 0.0001]:
                mt_path = f'gs://aou_analysis/v8/mt_results/{ancestry.upper()}_gene.mt'
                print(mt_path)
                mt = hl.read_matrix_table(mt_path)
                mt.describe()
                mt =  mt.filter_rows(mt.annotation != 'Cauchy')

                qc_ht = hl.read_table(f'{DATA_PATH}/qc/aou_exome_gene_qc.ht')
                if ancestry == 'meta':
                    qc_ht = qc_ht.annotate(CAF=qc_ht.CAF_global_raw)
                else:
                    qc_ht = qc_ht.annotate(CAF=qc_ht[f'CAF_{ancestry.upper()}_{maxMAF}'])

                mt = mt.annotate_rows(
                    mean_coverage=qc_ht[mt.annotation, mt.gene_id, mt.gene_symbol].mean_coverage,
                    CAF=qc_ht[mt.annotation, mt.gene_id, mt.gene_symbol].CAF)

                pheno_ht = hl.read_table('gs://aou_analysis/v8/data/phenotype/summary/aou_v8_phenotype_meta_info_qced.ht')
                pheno_ht = pheno_ht.filter((pheno_ht.ancestry == ancestry.upper()) & (pheno_ht.pheno_sex == 'both')).key_by('phenoname')
                mt = mt.annotate_cols(**pheno_ht[mt.phenoname])
                min_caf = 1e-5 if ancestry not in ['mid', 'sas'] else 1e-4
                print(f'MIN CAF: {min_caf}')

                mt = mt.annotate_entries(
                    hq_exp_CAC = (mt.CAF * mt.n_cases >= 5)
                )

                mt = mt.annotate_rows(
                    hq_coverage = mt.mean_coverage > 10,
                    hq_n_var_5 = mt.total_variants >= 5,
                    hq_n_var_10 = mt.total_variants > 10,
                    hq_caf = mt.CAF > min_caf,
                    hg_gene_exp_CAC = hl.agg.count_where(mt.hq_exp_CAC) > 0,
                )
                if ancestry == 'eur':
                    mt = mt.annotate_cols(hq_phenotype_lambda = ~hl.literal(phenotype_outliers[ancestry]).contains(mt.phenoname))
                mt = mt.annotate_cols(
                    hq_phenotype_exp_CAC_01  = hl.agg.filter(mt.max_MAF == 0.01, hl.agg.count_where(mt.hq_exp_CAC) > 0),
                    hq_phenotype_exp_CAC_001=hl.agg.filter(mt.max_MAF == 0.001, hl.agg.count_where(mt.hq_exp_CAC) > 0),
                    hq_phenotype_exp_CAC_0001=hl.agg.filter(mt.max_MAF == 0.0001, hl.agg.count_where(mt.hq_exp_CAC) > 0),
                )
                ht = mt.rows().checkpoint(f'{DATA_PATH}/qc/aou_{ancestry}_{maxMAF}_gene_qc_flags.ht', _read_if_exists=True)
                ht.export(f'{DATA_PATH}/qc/aou_{ancestry}_{maxMAF}_gene_qc_flags.txt.bgz')
                ht = mt.cols().checkpoint(f'{DATA_PATH}/qc/aou_{ancestry}_{maxMAF}_phenotype_qc_flags.ht', _read_if_exists=True)
                ht.export(f'{DATA_PATH}/qc/aou_{ancestry}_{maxMAF}_phenotype_qc_flags.txt.bgz')

    if args.compute_lambda_gc:
        for ancestry in ['meta', 'sas', 'mid', 'eur', 'eas', 'amr', 'afr']:
            mt_name = 'gene'
            if result_type == 'variant':
                mt_name = 'genome_variant' if name == 'ACAF' else 'exome_variant'
            mt_path = f'gs://aou_analysis/v8/mt_results/{ancestry.upper()}_{mt_name}.mt'
            print(mt_path)
            mt = hl.read_matrix_table(mt_path)
            pheno_ht = hl.read_table('gs://aou_analysis/v8/data/phenotype/summary/aou_v8_phenotype_meta_info_qced.ht')
            pheno_ht = pheno_ht.filter((pheno_ht.ancestry == ancestry.upper()) & (pheno_ht.pheno_sex == 'both')).key_by('phenoname')
            pheno_ht = pheno_ht.drop('ancestry')
            pheno_ht.show()
            pheno_ht.describe()
            if ancestry != 'meta':
                mt = mt.select_cols('heritability', 'saige_version', 'inv_normalized')
            mt = mt.annotate_cols(**pheno_ht[mt.phenoname])
            if result_type == 'variant':
                qc_ht = hl.read_table(f'{DATA_PATH}/qc/aou_{name}_variant_qc.ht')
                if args.ancestry == 'meta':
                    qc_ht = qc_ht.annotate(AF = qc_ht.freq.ALL.AF[1])
                    qc_ht.filter(hl.is_defined(qc_ht.AF)).show()
                    qc_ht.AF.summarize()
                else:
                    qc_ht = qc_ht.annotate(AF = qc_ht.freq[args.ancestry.upper()].AF[1])
                qc_ht = qc_ht.annotate(AF_bins = get_AF_bins(qc_ht.AF),
                                    coverage_bins = get_coverage_bins(qc_ht.approx_coverage))
                AF_bins = list(qc_ht.aggregate(hl.agg.collect_as_set(qc_ht.AF_bins)))
                if not ((args.ancestry == 'meta') and (mt_name == 'exome_variant')):
                    AF_bins.remove(None)
                coverage_bins = list(qc_ht.aggregate(hl.agg.collect_as_set(qc_ht.coverage_bins)))

                mt = mt.annotate_rows(
                    AF = qc_ht[mt.row_key].AF,
                    AF_bins = qc_ht[mt.row_key].AF_bins,
                    coverage_bins = qc_ht[mt.row_key].coverage_bins)
                if args.ancestry == 'meta':
                    mt = mt.annotate_entries(Pvalue = hl.exp(mt.meta_analysis[0].Pvalue),
                                            SE =mt.meta_analysis[0].SE)
                mt = mt.filter_entries(hl.is_defined(mt.Pvalue))
                mt = mt.annotate_entries(Pvalue = hl.if_else(mt.Pvalue == 0, 1e-300, mt.Pvalue))
                mt = mt.annotate_entries(Pvalue = hl.if_else(mt.Pvalue == 1, 0.999, mt.Pvalue))
                mt = mt.annotate_cols(lambda_gc_raw = hl.methods.statgen._lambda_gc_agg(mt.Pvalue))
                if args.filtered_lambda:
                    mt = mt.filter_rows(
                        (mt.AF > 0.001)
                    )
                    mt = mt.filter_entries(
                        mt.AF * mt.n_cases >= 5
                    )
                    mt = mt.annotate_cols(lambda_gc_filtered=hl.methods.statgen._lambda_gc_agg(mt.Pvalue), )

                mt = mt.annotate_cols(
                    lambda_gc_af = hl.agg.group_by(
                        mt.AF_bins, hl.methods.statgen._lambda_gc_agg(mt.Pvalue)
                    ),
                    lambda_gc_coverage=hl.agg.group_by(
                        mt.coverage_bins, hl.methods.statgen._lambda_gc_agg(mt.Pvalue)
                    ),
                )
                mt = mt.annotate_cols(
                    **{
                        f'lambda_gc_af_{bin.split(",")[0].replace("[", "")}': mt.lambda_gc_af.get(bin)
                        for bin in AF_bins
                    },
                    **{
                        f'lambda_gc_coverage_{bin.split(",")[0].replace("[", "")}': mt.lambda_gc_coverage.get(bin)
                        for bin in coverage_bins
                    }

                )
                mt = mt.drop('lambda_gc_af', 'lambda_gc_coverage')
                mt.describe()

                filter_tag = '_filtered' if args.filtered_lambda else ''
                lambda_gc_ht = mt.cols().checkpoint(f'{DATA_PATH}/qc/aou_{args.ancestry}_{mt_name}_phenotype_lambda_gc{filter_tag}.ht',
                                                    _read_if_exists=not args.overwrite, overwrite=args.overwrite)
                lambda_gc_ht.show()
                print(lambda_gc_ht.count())
                lambda_gc_ht.export(f'{DATA_PATH}/qc/aou_{args.ancestry}_{mt_name}_phenotype_lambda_gc{filter_tag}.txt.bgz')

            if analysis_type=='gene' and result_type=='gene':
                qc_ht = hl.read_table(f'{DATA_PATH}/qc/aou_{name}_gene_qc.ht')
                if ancestry == 'meta':
                    qc_ht = qc_ht.annotate(CAF = qc_ht[f'CAF_global_{args.maxMAF}'])
                else:
                    qc_ht = qc_ht.annotate(CAF = qc_ht[f'CAF_{ancestry.upper()}_{args.maxMAF}'])

                qc_ht = qc_ht.annotate(CAF_bins=get_AF_bins(qc_ht.CAF),
                                    coverage_bins=get_coverage_bins(qc_ht.mean_coverage))
                qc_ht = qc_ht.select('CAF', 'CAF_bins', 'mean_coverage', 'coverage_bins')
                CAF_bins = list(qc_ht.aggregate(hl.agg.collect_as_set(qc_ht.CAF_bins)))
                print(CAF_bins)
                coverage_bins = list(qc_ht.aggregate(hl.agg.collect_as_set(qc_ht.coverage_bins)))
                print(coverage_bins)
                mt = mt.key_rows_by(mt.gene_id, mt.gene_symbol, mt.annotation)
                mt = mt.annotate_rows(
                    mean_coverage=qc_ht[mt.row_key].mean_coverage,
                    CAF = qc_ht[mt.row_key].CAF,
                    CAF_bins = qc_ht[mt.row_key].CAF_bins,
                    coverage_bins = qc_ht[mt.row_key].coverage_bins)
                if ancestry == 'meta':
                    mt = mt.annotate_entries(Pvalue = mt[f'META_{P_VALUE_FIELDS[args.test_type]}'])
                else:
                    mt = mt.annotate_entries(Pvalue = mt[P_VALUE_FIELDS[args.test_type]])
                mt = mt.filter_entries(hl.is_defined(mt.Pvalue))
                mt = mt.filter_rows(mt.max_MAF == hl.float64(args.maxMAF))
                mt = mt.annotate_entries(Pvalue = hl.if_else(mt.Pvalue == 0, 1e-300, mt.Pvalue))
                mt = mt.annotate_cols(lambda_gc_raw = hl.methods.statgen._lambda_gc_agg(mt.Pvalue),)
                if args.filtered_lambda:
                    pheno_ht = hl.read_table('gs://aou_analysis/v8/data/phenotype/summary/aou_v8_phenotype_meta_info_qced.ht')
                    pheno_ht = pheno_ht.filter((pheno_ht.ancestry == ancestry.upper()) & (pheno_ht.pheno_sex == 'both')).key_by('phenoname')
                    mt = mt.annotate_cols(n_cases = pheno_ht[mt.phenoname].n_cases)
                    min_caf = 1e-5 if ancestry not in ['mid', 'sas'] else 1e-4
                    print(f'MIN CAF: {min_caf}')
                    if ancestry == 'meta':
                        mt = mt.annotate_entries(
                            total_variants = hl.max(mt.summary_stats.total_variants_pheno)
                        )
                    else:
                        mt = mt.annotate_entries(
                            total_variants = mt.total_variants_pheno
                        )
                    mt = mt.filter_rows(
                        (mt.mean_coverage > 10) & (mt.CAF > min_caf)
                    )
                    mt = mt.filter_entries(
                        (mt.total_variants > 10) &  (mt.CAF * mt.n_cases >= 5)
                    )
                    mt = mt.annotate_cols(lambda_gc_filtered=hl.methods.statgen._lambda_gc_agg(mt.Pvalue), )

                mt = mt.annotate_cols(
                    lambda_gc_caf = hl.agg.group_by(
                        mt.CAF_bins, hl.methods.statgen._lambda_gc_agg(mt.Pvalue)
                    ),
                    lambda_gc_coverage=hl.agg.group_by(
                        mt.coverage_bins, hl.methods.statgen._lambda_gc_agg(mt.Pvalue)
                    ),
                    lambda_gc_annotation=hl.agg.group_by(
                        mt.annotation, hl.methods.statgen._lambda_gc_agg(mt.Pvalue)
                    ),
                )
                mt.describe()
                mt = mt.annotate_cols(
                    **{
                        f'lambda_gc_caf_{bin.split(",")[0].replace("[", "")}': mt.lambda_gc_caf.get(bin)
                        for bin in CAF_bins
                    },
                    **{
                        f'lambda_gc_coverage_{bin.split(",")[0].replace("[", "")}': mt.lambda_gc_coverage.get(bin)
                        for bin in coverage_bins
                    },
                    ** {
                        f'lambda_gc_annotation_{ann}': mt.lambda_gc_annotation.get(ann)
                        for ann in ['pLoF', 'missenseLC', 'synonymous', 'pLoF;missenseLC']
                    }

                )
                mt = mt.drop('lambda_gc_caf', 'lambda_gc_coverage', 'lambda_gc_annotation')
                mt.describe()

                filter_tag = '_filtered' if args.filtered_lambda else ''
                lambda_gc_ht = mt.cols().checkpoint(f'{DATA_PATH}/qc/aou_{ancestry}_{name}_{analysis_type}_{args.maxMAF}_phenotype_lambda_gc{filter_tag}.ht', overwrite=True)
                lambda_gc_ht.show(200)
                lambda_gc_ht.export(f'{DATA_PATH}/qc/aou_{ancestry}_{name}_{analysis_type}_{args.maxMAF}_phenotype_lambda_gc{filter_tag}.txt.bgz')

    if args.compute_gene_lambda_gc:
        phenotype_outliers = ['PP_928.1', 'PP_932.11', '3011099', '3046664', '3013861', '3052648', '3045792_3016921']
        for ancestry in ['meta', 'sas', 'mid', 'eur', 'eas', 'amr', 'afr']:
            mt_path = f'gs://aou_analysis/v8/mt_results/{ancestry.upper()}_gene.mt'
            print(mt_path)
            mt = hl.read_matrix_table(mt_path)
            qc_ht = hl.read_table(f'{DATA_PATH}/qc/aou_exome_gene_qc.ht')
            pheno_ht = hl.read_table('gs://aou_analysis/v8/data/phenotype/summary/aou_v8_phenotype_meta_info_qced.ht')
            pheno_ht = pheno_ht.filter((pheno_ht.ancestry == ancestry.upper()) & (pheno_ht.pheno_sex == 'both')).key_by('phenoname')
            pheno_ht.show()
            pheno_ht.describe()
            if ancestry != 'meta':
                mt = mt.select_cols('heritability', 'saige_version', 'inv_normalized')
            mt = mt.annotate_cols(**pheno_ht[mt.phenoname])
            if ancestry != 'meta':
                qc_ht = qc_ht.select('mean_coverage', 'CAF_global_raw', 'CAF_global_0.01', 'CAF_global_0.001', 'CAF_global_0.0001',
                                    f'CAF_{ancestry.upper()}_raw', f'CAF_{ancestry.upper()}_0.01', f'CAF_{ancestry.upper()}_0.001', f'CAF_{ancestry.upper()}_0.0001')
            else:
                qc_ht = qc_ht.select('mean_coverage', 'CAF_global_raw', 'CAF_global_0.01', 'CAF_global_0.001', 'CAF_global_0.0001')
                qc_ht = qc_ht.annotate(**{f'CAF_{ancestry.upper()}_{maxmaf}': qc_ht[f'CAF_global_{maxmaf}'] for maxmaf in ['raw', '0.01', '0.001', '0.0001']})

            mt = mt.filter_rows(hl.is_defined(mt.max_MAF))
            mt = mt.filter_cols(~hl.literal(phenotype_outliers).contains(mt.phenoname))
            mt = mt.filter_cols(~mt.phenoname.startswith('random'))
            if ancestry == 'meta':
                mt = mt.annotate_entries(
                    Pvalue = mt.META_Pvalue_SKATO,
                    Pvalue_Burden=mt.META_Pvalue_Burden,
                    Pvalue_SKAT=mt.META_Pvalue_SKAT,

                )
            mt = mt.annotate_rows(lambda_gc_skato_gene=hl.methods.statgen._lambda_gc_agg(mt.Pvalue),
                                lambda_gc_burden_gene=hl.methods.statgen._lambda_gc_agg(mt.Pvalue_Burden),
                                lambda_gc_skat_gene=hl.methods.statgen._lambda_gc_agg(mt.Pvalue_SKAT),
                                **qc_ht[mt.gene_id, mt.gene_symbol, mt.annotation],
                                )
            mt = mt.annotate_rows(CAF = hl.case()
                                    .when(mt.max_MAF == 0.01, mt[f'CAF_{ancestry.upper()}_0.01'])
                                    .when((mt.max_MAF == 0.001), mt[f'CAF_{ancestry.upper()}_0.001'])
                                    .when(mt.max_MAF == 0.0001, mt[f'CAF_{ancestry.upper()}_0.0001'])
                                    .or_missing()
            )
            mt = mt.filter_entries(mt.CAF * mt.n_cases >= 5)
            mt = mt.annotate_rows(lambda_gc_skato_gene_filtered=hl.methods.statgen._lambda_gc_agg(mt.Pvalue),
                                lambda_gc_burden_gene_filtered=hl.methods.statgen._lambda_gc_agg(mt.Pvalue_Burden),
                                lambda_gc_skat_gene_filtered=hl.methods.statgen._lambda_gc_agg(mt.Pvalue_SKAT),
                                )
            ht = mt.rows()
            ht.describe()
            ht = ht.checkpoint(f'{DATA_PATH}/qc/aou_{ancestry}_gene_lambda_gc.ht', _read_if_exists = not args.overwrite, overwrite=args.overwrite)
            ht.export(f'{DATA_PATH}/qc/aou_{ancestry}_gene_lambda_gc.txt.bgz')
            ht.show()

    if args.compute_pheno_stats:
        def merge_pheno_stats():
            hl.init(
                master='local[80]',
                tmp_dir="gs://aou_tmp",
                default_reference="GRCh38",
            )
            ht= hl.import_table(f'gs://aou_wlu/v8_analysis/phenotype_info/*.txt',
                                 delimiter='\t', impute=True, types = {'n_controls': hl.tint64})
            ht.naive_coalesce(100).write('gs://aou_wlu/v8_analysis/phenotype_info.ht', overwrite=True)

        def export_pheno_stats(pop, pheno, category, category_name):
            quantitative_categories = ['physical_measurement', 'lab_measurement']

            pheno_stats = {
                'phenoname': [],
                'pheno_sex': [],
                'trait_type': [],
                'category': [],
                'n_cases': [],
                'n_controls': [],
                'pop': []
            }
            trait_type = 'continuous' if (category in quantitative_categories) or ('continuous' in pheno) else 'binary'
            ht = hl.import_table(f'gs://aou_analysis/v8/pheno_file/{pop.upper()}/phenotype_{pheno}.tsv',
                                 delimiter='\t', impute=True)
            n_female = ht.aggregate(hl.agg.count_where(ht.sex == 0))
            n_male = ht.aggregate(hl.agg.count_where(ht.sex == 1))
            pheno_sex = 'Female' if n_male == 0 else ('Male' if n_female == 0 else 'Both')
            if trait_type == 'continuous':
                n_cases = ht.aggregate(hl.agg.count_where(hl.is_defined(ht.value)))
                n_controls = hl.missing(hl.tint64)
            else:
                n_cases = ht.aggregate(hl.agg.count_where(ht.value == 1))
                n_controls = ht.aggregate(hl.agg.count_where(ht.value == 0))
            pheno_stats['phenoname'].append(pheno)
            pheno_stats['pheno_sex'].append(pheno_sex)
            pheno_stats['trait_type'].append(trait_type)
            pheno_stats['category'].append(category_name)
            pheno_stats['n_cases'].append(n_cases)
            pheno_stats['n_controls'].append(n_controls)
            pheno_stats['pop'].append(pop)
            df = pd.DataFrame(pheno_stats)
            ht = hl.Table.from_pandas(df, key = ['phenoname', 'pop'])
            ht.export(f'gs://aou_wlu/v8_analysis/phenotype_info/{pop.upper()}_{pheno}.txt')
            print(pheno_stats)

        if not hl.hadoop_exists('gs://aou_wlu/v8_analysis/phenotype_info.ht') or args.overwrite:
            backend = hb.ServiceBackend(
                billing_project="all-by-aou", remote_tmpdir='gs://aou_tmp/'
            )

            b = hb.Batch(
                name=f"aou_export_pheno_info",
                backend=backend,
                default_storage="500Mi",
                default_cpu=8,
            )
            CATEGORIES = ['physical_measurement', 'r_drug', 'pfhh_survey', 'random_pheno', 'lab_measurement', 'mcc2_phecodex', 'mhwb_survey']
            phenos_to_run_by_pop_by_group = read_pickle_dict(f'gs://aou_analysis/v8/data/phenotype/summary/pheno_by_ancestry_by_group_dict_both.dict')
            job_lst = []
            for pop in POPS:
                print(f'----------------{pop.upper()}-----------------')
                for category in CATEGORIES:
                    category_name = category
                    print(f'----------------{category_name}-----------------')
                    phenotypes = phenos_to_run_by_pop_by_group[pop][category]
                    for i in tqdm(range(len(phenotypes))):
                        pheno = phenotypes[i]
                        if pheno == 'ehr_year': continue
                        if not hl.hadoop_exists(f'gs://aou_wlu/v8_analysis/phenotype_info/{pop.upper()}_{pheno}.txt'):
                            j = b.new_python_job(name=f"aou_{pop.upper()}_{pheno}_info")
                            j.image("hailgenetics/hail:0.2.133-py3.11")
                            j.memory('lowmem')
                            j.call(export_pheno_stats,
                                   pop=pop,
                                   pheno=pheno,
                                   category=category,
                                   category_name=category_name
                                   )
                            job_lst.append(j)

                            if args.test:
                                break
                    if args.test:                        break
                if args.test:                            break
            j_merge = b.new_python_job(name=f"aou_merge_pheno_info")
            j_merge.depends_on(*job_lst)
            j_merge.cpu(16)
            j_merge.image("hailgenetics/hail:0.2.133-py3.11")
            j_merge.memory('standard')
            j_merge.call(merge_pheno_stats)

            b.run()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--single-variant-only", help="Run single variant SAIGE", action="store_true"
    )
    parser.add_argument(
        "--compute-lambda-gc", help="compute Lambda GC", action="store_true"
    )
    parser.add_argument(
        "--create-call-stats-ht", help="Create call stats HT", action="store_true"
    )
    parser.add_argument(
        "--count-hq-genes", help="Count number of genes/phenotypes after each filter", action="store_true"
    )
    parser.add_argument(
        "--filtered-lambda", help="compute lambda after filtering", action="store_true"
    )
    parser.add_argument(
        "--compute-pheno-stats", help="compute phenotype statistics", action="store_true"
    )
    parser.add_argument(
        "--export-variant-qc-info", help="Export variant qc info", action="store_true"
    )
    parser.add_argument(
        "--export-gene-qc-info", help="Export variant qc info", action="store_true"
    )
    parser.add_argument(
        "--ancestry", help="", default="meta"
    )
    parser.add_argument(
        "--maxMAF", help="", default="0.001"
    )
    parser.add_argument(
        "--test-type", help="", default="burden"
    )
    parser.add_argument(
        "--test", help="Run a test", action="store_true"
    )
    parser.add_argument(
        "--overwrite",
        help="Overwrite results",
        action="store_true",
    )
    parser.add_argument(
        "--load-gene-results", help="Load gene level test results", action="store_true"
    )
    parser.add_argument(
        "--compute-gene-lambda-gc", help="compute gene level lambda GC", action="store_true"
    )
    args = parser.parse_args()

    main(args)