import hail as hl
import copy
from collections import Counter
import pandas as pd
import pickle
import argparse
import logging
import hailtop.batch as hb
import hailtop.fs as hfs
from hailtop.batch.resource import Resource, ResourceGroup
from hailtop.batch.job import Job
from hailtop.batch.batch import Batch
from collections import Counter
from shlex import quote as shq
import time
from tqdm import tqdm

from aou_gwas import *  # for dataproc


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
    test_type = 'saige' if analysis_type == 'variant' else 'saige_gene'
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

    if args.export_variant_qc_info:
        path = ACAF_MT_PATH if analysis_type == 'variant' else EXOME_MT_PATH
        if not args.dataproc:
            mt = hl.read_matrix_table(path)
            ht = mt.rows().checkpoint(f'{DATA_PATH}/aou_{name}_rows_full.ht', _read_if_exists = True)
        ht = hl.read_table(f'{DATA_PATH}/aou_{name}_rows_full.ht')
        ht = ht.select(approx_coverage = ht.variant_qc.gq_stats.mean/3)
        call_stats_ht = hl.read_table(get_call_stats_ht_path(pop='full', pruned=True, analysis_type=analysis_type))
        call_stats_ht = call_stats_ht.filter(hl.is_defined(call_stats_ht.freq))
        ht = call_stats_ht.annotate(**ht[call_stats_ht.key])
        ht = ht.naive_coalesce(1000).checkpoint(f'{DATA_PATH}/utils/aou_{name}_variant_qc.ht', overwrite = True)
        ht.show()

    if args.export_gene_qc_info:
        ht = hl.read_table(f'{DATA_PATH}/utils/aou_exome_variant_qc.ht')
        ht = ht.annotate(annotation = hl.if_else(hl.literal({'missense', 'LC'}).contains(ht.annotation), 'missenseLC', ht.annotation))
        combined_ht = ht.filter(hl.literal(['missenseLC', 'pLoF']).contains(ht.annotation))
        combined_ht = combined_ht.group_by('gene_id', 'gene_symbol').aggregate(
            mean_coverage=hl.agg.mean(combined_ht.approx_coverage),
            CAF_global_raw=hl.agg.sum(combined_ht.freq.ALL.AF[1]),
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
                f"CAF_{pop}_raw": hl.agg.sum(combined_ht.freq[pop.upper()].AF[1])
                for pop in POPS
            },
            **{
                f"CAF_{pop}_{max_MAF}": hl.agg.filter(combined_ht.freq[pop.upper()].AF[1] < max_MAF,
                                                      hl.agg.sum(combined_ht.freq[pop.upper()].AF[1]))
                for pop in POPS
                for max_MAF in [0.01, 0.001, 0.0001]
            },
        )
        combined_ht = combined_ht.annotate(annotation = 'pLoF;missenseLC').key_by('gene_id', 'gene_symbol','annotation')
        fields = ['gene_id', 'gene_symbol','annotation'] + list(combined_ht.row_value)
        combined_ht = combined_ht.key_by().select(*fields).key_by('gene_id', 'gene_symbol','annotation')
        gene_ht = ht.group_by('gene_id', 'gene_symbol','annotation').aggregate(
            mean_coverage = hl.agg.mean(ht.approx_coverage),
            CAF_global_raw = hl.agg.sum(ht.freq.ALL.AF[1]),
            **{
                f"mean_coverage_{max_MAF}": hl.agg.filter(ht.freq.ALL.AF[1] < max_MAF,
                                                       hl.agg.mean(ht.approx_coverage))
                for max_MAF in [0.01, 0.001, 0.0001]
            },
            **{
                f"CAF_global_{max_MAF}": hl.agg.filter(ht.freq.ALL.AF[1] < max_MAF,
                                                       hl.agg.sum(ht.freq.ALL.AF[1]))
                for max_MAF in [0.01, 0.001, 0.0001]
            },
            **{
                  f"CAF_{pop}_raw": hl.agg.sum(ht.freq[pop.upper()].AF[1])
                  for pop in POPS
              },
            **{
                f"CAF_{pop}_{max_MAF}": hl.agg.filter(ht.freq[pop.upper()].AF[1] < max_MAF ,hl.agg.sum(ht.freq[pop.upper()].AF[1]))
                for pop in POPS
                for max_MAF in [0.01, 0.001, 0.0001]
            },
        )
        gene_ht = gene_ht.select(*list(combined_ht.row_value))
        gene_ht = gene_ht.union(combined_ht)
        gene_ht = gene_ht.checkpoint(f'{DATA_PATH}/utils/aou_exome_gene_qc.ht', overwrite=True)
        gene_ht.describe()
        gene_ht.show()
        print(gene_ht.count())
        gene_ht.export('gs://aou_wlu/250k_qc/aou_exome_gene_qc.txt.bgz')

    if args.count_hq_genes:
        phenotype_outliers = {
            'afr':['654'],
            'amr':['654'],
            'eur':['638', '654', '665', 'PP_927', 'PP_927.2'],
        }
        mt_path = get_aou_final_results_path(analysis_type='gene', pop=args.pop, result_type='gene',
                                             extension=f'mt')
        print(mt_path)
        mt = hl.read_matrix_table(mt_path)
        mt.describe()
        mt =  mt.filter_rows(mt.annotation != 'Cauchy')

        qc_ht = hl.read_table(f'{RESULTS_PATH}/utils/aou_3_annotation_exome_gene_info_pruned_250k.ht')
        if args.pop == 'meta':
            qc_ht = qc_ht.annotate(CAF=qc_ht.CAF_global_raw)
        else:
            qc_ht = qc_ht.annotate(CAF=qc_ht[f'CAF_{args.pop}_{args.maxMAF}'])

        mt = mt.annotate_rows(
            mean_coverage=qc_ht[mt.annotation, mt.gene_id, mt.gene_symbol].mean_coverage,
            CAF=qc_ht[mt.annotation, mt.gene_id, mt.gene_symbol].CAF)

        pheno_ht = hl.read_table('gs://aou_results/250k/v1.1/utils/aou_phenotype_meta_info.ht')
        pheno_ht = pheno_ht.filter(pheno_ht.pop == args.pop).key_by('phenoname')
        mt = mt.annotate_cols(**pheno_ht[mt.phenoname])
        min_caf = 1e-5 if args.pop not in ['mid', 'sas'] else 1e-4
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
        if args.pop in ['afr', 'amr', 'eur']:
            mt = mt.annotate_cols(hq_phenotype_lambda = ~hl.literal(phenotype_outliers[args.pop]).contains(mt.phenoname))
        mt = mt.annotate_cols(
            hq_phenotype_exp_CAC_01  = hl.agg.filter(mt.max_MAF == 0.01, hl.agg.count_where(mt.hq_exp_CAC) > 0),
            hq_phenotype_exp_CAC_001=hl.agg.filter(mt.max_MAF == 0.001, hl.agg.count_where(mt.hq_exp_CAC) > 0),
            hq_phenotype_exp_CAC_0001=hl.agg.filter(mt.max_MAF == 0.0001, hl.agg.count_where(mt.hq_exp_CAC) > 0),
        )
        ht = mt.rows().checkpoint(f'gs://aou_wlu/250k_qc/aou_{args.pop}_{args.maxMAF}_gene_qc_flags.ht', _read_if_exists=True)
        ht.export(f'gs://aou_wlu/250k_qc/aou_{args.pop}_{args.maxMAF}_gene_qc_flags.txt.bgz')
        ht = mt.cols().checkpoint(f'gs://aou_wlu/250k_qc/aou_{args.pop}_{args.maxMAF}_phenotype_qc_flags.ht', _read_if_exists=True)
        ht.export(f'gs://aou_wlu/250k_qc/aou_{args.pop}_{args.maxMAF}_phenotype_qc_flags.txt.bgz')


    if args.compute_lambda_gc:
        mt_path = get_aou_final_results_path(analysis_type=analysis_type, pop=args.pop, result_type=result_type,
                                             extension=f'mt')
        print(mt_path)
        mt = hl.read_matrix_table(mt_path)
        pheno_ht = hl.read_table('gs://aou_results/250k/v1.1/utils/aou_phenotype_meta_info.ht')
        pheno_ht = pheno_ht.drop('lambda_gc_exome', 'lambda_gc_acaf', 'lambda_gc_gene_burden_001')
        pheno_ht = pheno_ht.filter(pheno_ht.pop == args.pop).key_by('phenoname')
        if args.pop != 'meta':
            mt = mt.select_cols('heritability', 'saige_version', 'inv_normalized')
        mt = mt.annotate_cols(**pheno_ht[mt.phenoname])
        if result_type == 'variant':
            qc_ht = hl.read_table(f'{DATA_PATH}/utils/aou_{name}_variant_qc.ht')
            if args.pop == 'meta':
                qc_ht = qc_ht.annotate(AF = qc_ht.freq.ALL.AF[1])
                qc_ht.filter(hl.is_defined(qc_ht.AF)).show()
                qc_ht.AF.summarize()
            else:
                qc_ht = qc_ht.annotate(AF = qc_ht.freq[args.pop.upper()].AF[1])
            qc_ht = qc_ht.annotate(AF_bins = get_AF_bins(qc_ht.AF),
                                   coverage_bins = get_coverage_bins(qc_ht.approx_coverage))
            AF_bins = list(qc_ht.aggregate(hl.agg.collect_as_set(qc_ht.AF_bins)))
            AF_bins.remove(None)
            coverage_bins = list(qc_ht.aggregate(hl.agg.collect_as_set(qc_ht.coverage_bins)))

            mt = mt.annotate_rows(
                AF = qc_ht[mt.row_key].AF,
                AF_bins = qc_ht[mt.row_key].AF_bins,
                coverage_bins = qc_ht[mt.row_key].coverage_bins)
            if args.pop == 'meta':
                mt = mt.annotate_entries(Pvalue = mt.meta_analysis[0].Pvalue,
                                         SE =mt.meta_analysis[0].SE)
            mt = mt.filter_entries(hl.is_defined(mt.Pvalue))
            mt = mt.annotate_entries(Pvalue = hl.if_else(mt.Pvalue == 0, 1e-300, mt.Pvalue))
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
            lambda_gc_ht = mt.cols().checkpoint(f'gs://aou_wlu/250k_qc/aou_{args.pop}_{name}_{analysis_type}_phenotype_lambda_gc{filter_tag}.ht',
                                                _read_if_exists=not args.overwrite, overwrite=args.overwrite)
            lambda_gc_ht.show()
            print(lambda_gc_ht.count())
            lambda_gc_ht.export(f'gs://aou_wlu/250k_qc/aou_{args.pop}_{name}_{analysis_type}_phenotype_lambda_gc{filter_tag}.txt.bgz')


        if analysis_type=='gene' and result_type=='gene':
            qc_ht = hl.read_table(f'{DATA_PATH}/utils/aou_{name}_gene_qc.ht')
            if args.pop == 'meta':
                qc_ht = qc_ht.annotate(CAF = qc_ht.CAF_global_raw)
            else:
                qc_ht = qc_ht.annotate(CAF = qc_ht[f'CAF_{args.pop}_{args.maxMAF}'])

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
            if args.pop == 'meta':
                mt = mt.annotate_entries(Pvalue = mt[f'META_{P_VALUE_FIELDS[args.test_type]}'])
            else:
                mt = mt.annotate_entries(Pvalue = mt[P_VALUE_FIELDS[args.test_type]])
            mt = mt.filter_entries(hl.is_defined(mt.Pvalue))
            mt = mt.filter_rows(mt.max_MAF == hl.float64(args.maxMAF))
            mt = mt.annotate_entries(Pvalue = hl.if_else(mt.Pvalue == 0, 1e-300, mt.Pvalue))
            mt = mt.annotate_cols(lambda_gc_raw = hl.methods.statgen._lambda_gc_agg(mt.Pvalue),)
            if args.filtered_lambda:
                pheno_ht = hl.read_table('gs://aou_results/250k/v1.1/utils/aou_phenotype_meta_info.ht')
                pheno_ht = pheno_ht.filter(pheno_ht.pop == args.pop).key_by('phenoname')
                mt = mt.annotate_cols(n_cases = pheno_ht[mt.phenoname].n_cases)
                min_caf = 1e-5 if args.pop not in ['mid', 'sas'] else 1e-4
                print(f'MIN CAF: {min_caf}')
                mt = mt.filter_rows(
                    (mt.total_variants > 10) & (mt.mean_coverage > 10) & (mt.CAF > min_caf)
                )
                mt = mt.filter_entries(
                    mt.CAF * mt.n_cases >= 5
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
            lambda_gc_ht = mt.cols().write(f'gs://aou_wlu/250k_qc/aou_{args.pop}_{name}_{analysis_type}_{args.test_type}_{args.maxMAF}_phenotype_lambda_gc{filter_tag}_expected_ac_5.ht', overwrite=True)
            lambda_gc_ht.show(200)
            lambda_gc_ht.export(f'gs://aou_wlu/250k_qc/aou_{args.pop}_{name}_{analysis_type}_{args.test_type}_{args.maxMAF}_phenotype_lambda_gc{filter_tag}_expected_ac_5.txt.bgz')

    if args.compute_pheno_stats:
        def merge_pheno_stats():
            hl.init(
                master='local[80]',
                tmp_dir="gs://aou_tmp",
                default_reference="GRCh38",
            )
            ht= hl.import_table(f'gs://aou_wlu/250k_qc/phenotype_info/*.txt',
                                 delimiter='\t', impute=True, types = {'n_controls': hl.tint64})
            ht.naive_coalesce(100).write('gs://aou_wlu/250k_qc/phenotype_info.ht', overwrite=True)
        def export_pheno_stats(pop, pheno, category, category_name):
            quantitative_categories = ["lab_measurements", "processed_physical_measurement_table"]

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
            ht = hl.import_table(f'gs://aou_analysis/250k/pheno_file/{pop.upper()}/phenotype_{pheno}.tsv',
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
            ht.export(f'gs://aou_wlu/250k_qc/phenotype_info/{pop.upper()}_{pheno}.txt')
            print(pheno_stats)

        def combine_lambda_gc_ht(name: str, overwrite:bool):
            ht = hl.read_table(f'gs://aou_wlu/250k_qc/aou_afr_{name}.ht')
            ht = ht.annotate(pop='afr')
            ht = ht.select(lambda_gc=ht.lambda_gc_filtered, pop=ht.pop)
            for pop in ['amr', 'eas', 'eur', 'mid', 'sas', 'meta']:
                if hl.hadoop_exists(f'gs://aou_wlu/250k_qc/aou_{pop}_{name}.ht/_SUCCESS'):
                    tmp_ht = hl.read_table(f'gs://aou_wlu/250k_qc/aou_{pop}_{name}.ht')
                    tmp_ht = tmp_ht.annotate(pop=pop)
                    tmp_ht = tmp_ht.select(lambda_gc=tmp_ht.lambda_gc_filtered, pop=tmp_ht.pop)
                    ht = ht.union(tmp_ht)
            ht = ht.key_by('phenoname', 'pop')
            ht = ht.checkpoint(f'gs://aou_wlu/250k_qc/aou_full_{name.replace("_expected_ac_5", "")}.ht', overwrite=overwrite, _read_if_exists=not overwrite)
            return ht

        def combine_phenotype_descriptions(overwrite:bool):
            pheno_info = hl.import_table('gs://aou_analysis/250k/data/phenotype/phecode_definition_cleaned.tsv',
                                         impute=True)
            lab_info_raw = hl.read_table('gs://aou_analysis/250k/data/phenotype/labs_metadata_04242024_250k.ht')
            drug_info = hl.import_table('gs://aou_analysis/250k/data/phenotype/drug_definition.tsv', impute=True)
            pfhh_info = hl.import_table('gs://aou_analysis/250k/data/phenotype/pfhh_description.tsv', impute=True,
                                        types={'phenoname': hl.tstr})
            pheno_info = pheno_info.select(phenoname=pheno_info.phecode,
                                           description=pheno_info.description,
                                           phecode_category=pheno_info.category)
            lab_info = lab_info_raw.key_by()
            lab_info = lab_info.select(phenoname=lab_info.measurement_concept_id,
                                       description=lab_info.lab_name,
                                       phecode_category=hl.missing(hl.tstr))
            drug_info = drug_info.select(phenoname=drug_info.phenoname,
                                         description=drug_info.description,
                                         phecode_category=hl.missing(hl.tstr))
            pfhh_info = pfhh_info.select(phenoname=pfhh_info.phenoname,
                                         description=pfhh_info.description,
                                         phecode_category=hl.missing(hl.tstr))
            info = pheno_info.union(lab_info).union(drug_info).union(pfhh_info).key_by('phenoname')
            info = info.annotate(description_more=lab_info_raw[info.phenoname]['\ufeffstandard_concept_name'])
            info = info.checkpoint('gs://aou_analysis/250k/data/phenotype/FINAL_phenotype_description_250k.ht', overwrite=overwrite, _read_if_exists=not overwrite)
            return info



        if not hl.hadoop_exists('gs://aou_wlu/250k_qc/phenotype_info.ht') or args.overwrite:
            backend = hb.ServiceBackend(
                billing_project="all-by-aou", remote_tmpdir='gs://aou_tmp/'
            )

            b = hb.Batch(
                name=f"aou_export_pheno_info",
                backend=backend,
                default_storage="500Mi",
                default_cpu=8,
            )
            phenos_to_run_by_pop_by_group = read_pickle_dict(
                get_phenotype_info_path(version="pheno_by_pop_by_group_dict", extension="dict"))
            job_lst = []
            for pop in POPS:
                print(f'----------------{pop.upper()}-----------------')
                for category in phenotype_categories:
                    category_name = CATEGORIES[category]
                    print(f'----------------{category_name}-----------------')
                    phenotypes = phenos_to_run_by_pop_by_group[pop][category]
                    for i in tqdm(range(len(phenotypes))):
                        pheno = phenotypes[i]
                        if pheno == 'ehr_year': continue
                        if not hl.hadoop_exists(f'gs://aou_wlu/250k_qc/phenotype_info/{pop.upper()}_{pheno}.txt'):
                            j = b.new_python_job(name=f"aou_{pop.upper()}_{pheno}_info")
                            j.image("hailgenetics/hail:0.2.130-py3.9")
                            j.memory('lowmem')
                            j.call(export_pheno_stats,
                                   pop=pop,
                                   pheno=pheno,
                                   category=category,
                                   category_name=category_name
                                   )
                            job_lst.append(j)
            j_merge = b.new_python_job(name=f"aou_merge_pheno_info")
            j_merge.depends_on(*job_lst)
            j_merge.cpu(16)
            j_merge.image("hailgenetics/hail:0.2.127-py3.9")
            j_merge.memory('standard')
            j_merge.call(merge_pheno_stats)

            b.run()

        pheno_ht = hl.read_table('gs://aou_wlu/250k_qc/phenotype_info.ht')
        pheno_ht = pheno_ht.annotate(n_cases=hl.int64(pheno_ht.n_cases))
        pheno_ht = pheno_ht.filter(~hl.literal(LAB_TO_REMOVE).contains(pheno_ht.phenoname))
        pheno_ht = pheno_ht.key_by('phenoname', 'pop')
        meta_ht = pheno_ht.group_by('phenoname','pheno_sex', 'trait_type', 'category').aggregate(
            n_cases=hl.agg.sum(pheno_ht.n_cases),
            n_controls=hl.agg.sum(pheno_ht.n_controls),
        )
        meta_ht = meta_ht.annotate(pop = 'meta').key_by('phenoname', 'pop')
        pheno_ht = pheno_ht.union(meta_ht)


        ht_exome = combine_lambda_gc_ht(name='exome_gene_phenotype_lambda_gc_filtered', overwrite=args.overwrite)
        ht_acaf = combine_lambda_gc_ht(name='ACAF_variant_phenotype_lambda_gc_filtered', overwrite=args.overwrite)
        ht_gene_001 = combine_lambda_gc_ht(name='exome_gene_burden_0.001_phenotype_lambda_gc_filtered_expected_ac_5', overwrite=args.overwrite)
        info = combine_phenotype_descriptions(overwrite=args.overwrite)

        pheno_ht = pheno_ht.annotate(**info[pheno_ht.phenoname],
                                     lambda_gc_exome=ht_exome[pheno_ht.phenoname, pheno_ht.pop].lambda_gc,
                                     lambda_gc_acaf=ht_acaf[pheno_ht.phenoname, pheno_ht.pop].lambda_gc,
                                     lambda_gc_gene_burden_001=ht_gene_001[pheno_ht.phenoname, pheno_ht.pop].lambda_gc,
                                     )
        pheno_ht = pheno_ht.annotate(
            description=hl.if_else(pheno_ht.category == 'physical_measurement', pheno_ht.phenoname,
                                   pheno_ht.description))
        pheno_ht = pheno_ht.annotate(
            description=hl.if_else(pheno_ht.category == 'random', 'random_${heritability}_${prevalence}_${index}',
                                   pheno_ht.description))


        pheno_ht = pheno_ht.naive_coalesce(50).checkpoint('gs://aou_results/250k/v1.1/utils/aou_phenotype_meta_info.ht',
                                               overwrite=True)
        pheno_ht.export('gs://aou_wlu/250k_qc/aou_phenotype_meta_info_250k.txt.bgz')

        ht = hl.read_table('gs://aou_results/250k/v1.1/utils/aou_phenotype_meta_info.ht')
        ht.describe()
        ht.show()
        print(ht.count())
        print(ht.aggregate(hl.agg.counter(ht.pop)))
        print(ht.aggregate(hl.agg.counter(ht.pop + ': ' + ht.category)))




    if args.update_additional_variant_info:
        HT_MPC = 'gs://aou_wlu/utils/mpc/mpc_grch38_deduped_2024-03-13.ht'
        HT_AM = 'gs://nnfc-fdp-konrad-public/AlphaMissense/AlphaMissense_hg38.ht'
        HT_TX = 'gs://aou_wlu/250k_analysis/tx_annotation_liftover_grch38.ht'
        HT_ANNT = 'gs://aou_wlu/250k_brava/variant_annotation_cleaned_brava_annotated.ht'
        HT_VEP = 'gs://aou_wlu/250k_brava/vep_cleaned.ht'
        # HT_TX = "gs://gcp-public-data--gnomad/papers/2019-tx-annotation/pre_computed/all.possible.snvs.tx_annotated.021520.ht"
        # tx_ht = hl.read_table(HT_TX)
        # tx_ht = tx_ht.explode(tx_ht.tx_annotation)
        # tx_ht = add_liftover_rg37_to_rg38_ht(tx_ht)
        # tx_ht = tx_ht.checkpoint('gs://aou_wlu/250k_analysis/tx_annotation_liftover_grch38.ht', _read_if_exists=True)

        if not hl.hadoop_exists(HT_VEP) or args.overwrite:
            indel_vep_ht = hl.read_table(get_aou_util_path(name="vep")).select('vep')
            snp_vep_ht = hl.read_table(
                'gs://gcp-public-data--gnomad/resources/context/grch38_context_vep_annotated.v105.ht').select('vep')
            indel_fields_to_drop = ['ancestral', 'context']
            shared_fields_to_drop = ['uniparc', 'trembl', 'swissprot']
            indel_ht = indel_vep_ht.annotate(vep=indel_vep_ht.vep.drop('minimised'))
            indel_ht = indel_ht.annotate(
                vep=indel_ht.vep.annotate(
                    intergenic_consequences=indel_ht.vep.intergenic_consequences.map(
                        lambda x: x.drop(*indel_fields_to_drop)),
                    motif_feature_consequences=indel_ht.vep.motif_feature_consequences.map(
                        lambda x: x.drop(*indel_fields_to_drop)),
                    regulatory_feature_consequences=indel_ht.vep.regulatory_feature_consequences.map(
                        lambda x: x.drop(*indel_fields_to_drop)),
                    transcript_consequences=indel_ht.vep.transcript_consequences.map(
                        lambda x: x.drop(*shared_fields_to_drop, *indel_fields_to_drop)),
                ))
            snp_ht = snp_vep_ht.annotate(vep=snp_vep_ht.vep.drop('context'))
            snp_ht = snp_ht.annotate(
                vep=snp_ht.vep.annotate(
                    intergenic_consequences=snp_ht.vep.intergenic_consequences.map(lambda x: x.drop('minimised')),
                    motif_feature_consequences=snp_ht.vep.motif_feature_consequences.map(lambda x: x.drop('minimised')),
                    regulatory_feature_consequences=snp_ht.vep.regulatory_feature_consequences.map(
                        lambda x: x.drop('minimised')),
                    transcript_consequences=snp_ht.vep.transcript_consequences.map(
                        lambda x: x.drop(*shared_fields_to_drop, 'minimised')),
                ))
            vep_ht = snp_ht.union(indel_ht, unify=True)
            from gnomad.utils.vep import process_consequences
            process_vep_ht = process_consequences(vep_ht)
            vep_ht = vep_ht.annotate(
                worst_csq_by_gene_canonical=process_vep_ht[vep_ht.key].vep.worst_csq_by_gene_canonical)
            vep_ht.describe()
            vep_ht = vep_ht.checkpoint(HT_VEP, overwrite=True)

        clinvar_ht = annotate_clinvar_pathogenicity_ht()
        mpc_ht = hl.read_table(HT_MPC).key_by('locus', 'alleles')
        am_ht = hl.read_table(HT_AM)
        tx_ht = hl.read_table(HT_TX)
        qc_ht = hl.read_table(f'{DATA_PATH}/utils/aou_{name}_variant_qc.ht')
        annt_ht = hl.read_table(HT_ANNT).drop('lof', 'most_severe_csq_variant', 'most_severe_csq_gene', 'cadd', 'gnomad_revel', 'gnomad_spliceai',
                                              'brava_mask')
        vep_ht = hl.read_table(HT_VEP)

        print(get_call_stats_ht_path(pop='full', pruned=True, analysis_type=analysis_type))
        af_ht = hl.read_table(get_call_stats_ht_path(pop='full', pruned=True, analysis_type=analysis_type))
        af_ht.describe()
        print(af_ht.count())
        af_ht = af_ht.drop('lof', 'most_severe_consequence')
        af_ht = af_ht.annotate(**vep_ht[af_ht.key],
                               tx_annotation_csq=tx_ht[af_ht.key].tx_annotation.csq,
                               mean_proportion=tx_ht[af_ht.key].tx_annotation.mean_proportion,
                               pathogenicity=clinvar_ht[af_ht.key]["pathogenicity"],
                               MPC=mpc_ht[af_ht.key].mpc,
                               AM_pathogenicity=am_ht[af_ht.key].am_pathogenicity,
                               AM_class=am_ht[af_ht.key].am_class,
                               coverage=qc_ht[af_ht.key].approx_coverage,
                               **annt_ht[af_ht.key]
                               )
        # af_ht = af_ht.annotate(coverage_bins = get_coverage_bins(af_ht.coverage))
        af_ht = af_ht.filter(hl.is_defined(af_ht.freq))
        af_ht = af_ht.checkpoint(
            f'gs://aou_results/250k/v1.1/utils/aou_all_{name}_variant_info_pruned_250k_annotated_filtered.ht', overwrite=True)
        af_ht.describe()
        print(af_ht.count())

    if args.export_variant_sig_data:
        if result_type == 'variant':
            af_ht = hl.read_table(f'gs://aou_results/250k/v1.1/utils/aou_all_{name}_variant_info_pruned_250k_annotated_filtered.ht')
            for pop in POPS:
                print(pop.upper())
                mt = hl.read_matrix_table(
                    get_aou_final_results_path(analysis_type=analysis_type, pop=pop, result_type=result_type,
                                               extension=f'mt'))
                lambda_ht = hl.read_table(f'gs://aou_wlu/250k_qc/aou_{pop}_{name}_{analysis_type}_phenotype_lambda_gc.ht')
                lambda_ht = lambda_ht.filter((lambda_ht.lambda_gc < 2) & (lambda_ht.lambda_gc > 0.5))
                pheno_lst = list(lambda_ht.aggregate(hl.agg.collect_as_set(lambda_ht.phenoname)))
                mt = mt.select_rows()
                mt = mt.filter_cols(hl.literal(pheno_lst).contains(mt.phenoname))
                mt = mt.annotate_rows(n_sig=hl.agg.count_where(mt.Pvalue < 5e-8))
                ht = mt.rows()
                sub_af_ht = af_ht.annotate(AF=af_ht.freq[pop.upper()].AF[1])
                sub_af_ht = sub_af_ht.drop('freq', 'gene_symbol', 'gene_id', 'transcript_id', 'amino_acids', 'variant_id')
                if analysis_type == 'variant':
                    sub_af_ht = sub_af_ht.drop('freq_raw')
                ht = ht.join(sub_af_ht, 'left')
                ht = ht.annotate(AF_bin=get_AF_bins(ht.AF))

                ht.describe()

                ht.export(f'gs://aou_wlu/250k_analysis/af_analysis/{pop}_{name}_variant_n_sig_af.txt.bgz')

        if (result_type == 'gene') and (analysis_type == 'gene'):
            caf_ht = hl.read_table(f'{DATA_PATH}/utils/aou_exome_gene_qc.ht')
            for pop in pops:
                print(pop.upper())
                mt = hl.read_matrix_table(get_aou_final_results_path(analysis_type='gene', pop=pop, result_type='gene', extension=f'mt'))
                mt = mt.select_rows('MAC', 'Number_rare', 'Number_ultra_rare', 'total_variants')
                mt = mt.annotate_rows(n_burden_sig = hl.agg.count_where(mt.Pvalue_Burden < 2.5e-6),
                                      n_skato_sig = hl.agg.count_where(mt.Pvalue < 2.5e-6),
                                      n_skat_sig = hl.agg.count_where(mt.Pvalue_SKAT < 2.5e-6),)
                ht = mt.rows().key_by('gene_id', 'gene_symbol', 'annotation')
                ht = ht.join(caf_ht, 'left')
                ht = ht.annotate(CAF_bin = get_AF_bins(ht[f'CAF_{pop}']),
                                 CAF_001_bin = get_AF_bins(ht[f'CAF_{pop}_0.001']))

                ht.describe()
                ht.export(f'gs://aou_wlu/250k_analysis/caf_analysis/{pop}_gene_n_sig_caf.txt.bgz')





if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--single-variant-only", help="Run single variant SAIGE", action="store_true"
    )
    parser.add_argument(
        "--compute-lambda-gc", help="compute Lambda GC", action="store_true"
    )
    parser.add_argument(
        "--count-hq-genes", help="Count number of genes/phenotypes after each filter", action="store_true"
    )
    parser.add_argument(
        "--filtered-lambda", help="compute lambda after filtering", action="store_true"
    )
    parser.add_argument(
        "--compute-expected-p", help="compute expected pvalues", action="store_true"
    )
    parser.add_argument(
        "--compute-pheno-stats", help="compute phenotype statistics", action="store_true"
    )
    parser.add_argument(
        "--export-variant-qc-info", help="Export variant qc info", action="store_true"
    )
    parser.add_argument(
        "--export-variant-sig-data", help="Export variant sig info", action="store_true"
    )
    parser.add_argument(
        "--export-gene-qc-info", help="Export variant qc info", action="store_true"
    )
    parser.add_argument(
        "--update-additional-variant-info", help="Export additional variant info", action="store_true"
    )
    parser.add_argument(
        "--dataproc", help="Export variant qc info", action="store_true"
    )
    parser.add_argument(
        "--pop", help="", default="meta"
    )
    parser.add_argument(
        "--maxMAF", help="", default="0.001"
    )
    parser.add_argument(
        "--test-type", help="", default="burden"
    )
    parser.add_argument(
        "--overwrite",
        help="Overwrite results",
        action="store_true",
    )
    parser.add_argument(
        "--load-variant-results", help="Load single-variant level test results", action="store_true"
    )
    parser.add_argument(
        "--load-gene-results", help="Load gene level test results", action="store_true"
    )
    parser.add_argument(
        "--category",
        help="phenotype category to export",
    )
    args = parser.parse_args()

    main(args)