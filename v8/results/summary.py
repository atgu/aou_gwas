#!/usr/bin/env python3

__author__ = "Wenhan Lu"

import hail as hl
import hailtop.batch as hb
import argparse
import pandas as pd
import pickle
import hailtop.fs as hfs
from tqdm import tqdm

TMP_BUCKET = 'gs://aou_tmp_30_days'
RESULT_PATH = 'gs://aou_results/414k'
ANALYSIS_PATH = 'gs://aou_analysis/v8'
SUMMARY_ROOT = 'gs://aou_wlu/v8_analysis/paper_preparation/summary'
MT_ROOT = f'{RESULT_PATH}/mt_results'
P_THRESHOLDS = {'ACAF': 5e-8, 'Exome': 5e-8, 'Gene': 6.7e-7}
ANCESTRIES = ['AFR', 'AMR', 'EAS', 'EUR', 'MID', 'SAS', 'META']
ANNOTATIONS = ['pLoF', 'missenseLC', 'synonymous', 'pLoF;missenseLC']
CATEGORIES = ['lab_measurement', 'pfhh_survey', 'mhwb_survey', 'mcc2_phecodex', 'physical_measurement', 'r_drug']
P_VALUE_FIELDS = {"skato": "Pvalue", "skat": "Pvalue_SKAT", "burden": "Pvalue_Burden", "genome_variant": 'Pvalue', "exome_variant": 'Pvalue'}
META_P_FIELDS = {'Pvalue': 'META_Pvalue_SKATO', 'Pvalue_SKAT': 'META_Pvalue_SKAT', 'Pvalue_Burden': 'META_Pvalue_Burden'}
P_VALUE_FIELDS['gene'] = 'Pvalue_Burden'


UKB_GENE_MT = 'gs://ukbb-exome-public/500k/results/results.mt'
UKB_EXOME_MT = 'gs://ukbb-exome-public/500k/results/variant_results.mt'
UKB_GENOME_MT = 'gs://ukb-diverse-pops-public/sumstats_release/meta_analysis.mt'
UKB_GENE_PHENO_HT = 'gs://ukbb-exome-public/500k/results/pheno_results.ht'
AOU_GENE_MT = f'{MT_ROOT}/META_gene.mt'
AOU_EUR_GENE_MT = f'{MT_ROOT}/EUR_gene.mt'
PHENO_MAP_TSV = 'gs://aou_wlu/v8_analysis/aou_ukb_meta/aou_ukb_matched_all_phenotype_v8.csv'

def export_variant_signals(ancestry: str, test_type: str, qc: bool = True, overwrite: bool = False):
    test_label = 'genome' if test_type == 'ACAF' else 'exome'
    mt_path = f'{MT_ROOT}/{ancestry.upper()}_{test_label}_variant.mt'
    print(f'Processing {mt_path}...')
    mt = hl.read_matrix_table(mt_path)
    mt.describe()
    mt = mt.filter_cols(hl.literal(CATEGORIES).contains(mt.category))
    mt = mt.filter_cols(mt.hq_pheno)

    if ancestry.lower() == 'meta':
        mt = mt.annotate_cols(meta_pop = mt.ancestries)
        mt = mt.select_cols('description', 'category', 'meta_pop', 'n_cases', 'n_controls')
        mt = mt.annotate_entries(Pvalue = hl.exp(mt.Pvalue))
    else:
        mt = mt.select_cols('description', 'category', 'n_cases', 'n_controls')

    mt = mt.annotate_entries(sig = mt.Pvalue < P_THRESHOLDS[test_type])
    
    qc_label = '_qced' if qc else '_raw'
    if qc:
        mt = mt.filter_rows(mt.hq_variant)
        mt = mt.filter_entries(mt.hq_exp_AC & mt.sig)

    mt = mt.filter_entries(mt.sig)
    if ancestry.lower() == 'meta':
        mt = mt.select_entries('Pvalue', 'BETA', 'AF_Allele2', 'Q')
    else:
        mt = mt.select_entries('Pvalue', 'BETA', 'AF_Allele2')
    mt = mt.select_rows()

    ht = mt.entries()
    ht = ht.checkpoint(f'{SUMMARY_ROOT}/aou_v8_{ancestry}_{test_type}_signals{qc_label}.ht', overwrite=overwrite)
    ht.describe()
    print(ht.count())
    ht.export(f'{SUMMARY_ROOT}/aou_v8_{ancestry}_{test_type}_signals{qc_label}.txt.bgz')

def export_gene_signals(ancestry: str, p_field: str, qc: bool = True, overwrite: bool = False):
    mt = hl.read_matrix_table(f'{MT_ROOT}/{ancestry.upper()}_gene.mt')

    qc_label = '_qced' if qc else '_raw'
    if qc:
        mt = mt.filter_rows(mt.hq_gene)
        mt = mt.filter_cols(mt.hq_pheno)
        mt = mt.filter_entries(mt.hq_exp_CAC)
    
    pheno_ht = hl.read_table(f'{ANALYSIS_PATH}/data/phenotype/summary/aou_v8_phenotype_meta_info_qced.ht')
    pheno_ht = pheno_ht.key_by('phenoname')
    pheno_ht = pheno_ht.select('category', 'disease_category', 'description').distinct()
    mt = mt.annotate_cols(**pheno_ht[mt.col_key])
    mt = mt.filter_cols(hl.literal(CATEGORIES).contains(mt.category))

    phenotype_outliers = ['PP_928.1', 'PP_932.11']
    phenotype_to_flag = ['3011099', '3046664', '3013861', '3052648', '3045792_3016921']
    mt = mt.filter_cols(
        (~hl.literal(phenotype_outliers).contains(mt.phenoname))
    ) 
    mt = mt.annotate_cols(hq_pheno=~hl.literal(phenotype_to_flag).contains(mt.phenoname))
    mt = mt.filter_cols(mt.hq_pheno)

    if ancestry.lower() == 'meta':
        mt = mt.annotate_entries(**{f'{p_field}': mt[META_P_FIELDS[p_field]]})
    mt = mt.filter_rows((hl.literal(ANNOTATIONS).contains(mt.annotation)))
    mt = mt.annotate_entries(sig = mt[p_field] < P_THRESHOLDS['Gene'],
                             sig_1e_5 = mt[p_field] < 1e-5)
    mt.describe()
    if ancestry.lower() == 'meta':
        mt = mt.annotate_entries(BETA_Burden = mt.META_Stats_Burden)

    if ancestry.lower() == 'meta':
        mt = mt.select_cols('description', 'category',
                            N = hl.sum(mt.pheno_data.N_1),
                            n_cases = hl.sum(mt.pheno_data.n_cases),
                            n_controls = hl.sum(mt.pheno_data.n_controls),
                            ancestries = mt.pheno_data.ancestry
                            )
    else:
        mt = mt.select_cols('description', 'category')

    mt = mt.filter_entries(mt.sig_1e_5)
    mt = mt.select_entries('sig', p_field, 'BETA_Burden')
    mt = mt.select_rows()

    ht = mt.entries()
    ht.describe()
    ht = ht.select_globals()
    ht = ht.annotate(ancestry = ancestry)
    
    ht = ht.checkpoint(f'{SUMMARY_ROOT}/aou_v8_{ancestry}_{p_field}_gene_signals{qc_label}.ht', overwrite=overwrite, _read_if_exists=not overwrite)

    for maxMAF in [0.001, 0.0001, 0.01, -1]:
        sub_ht = ht.filter((ht.max_MAF == maxMAF))
        maxMAF = maxMAF if maxMAF != -1 else 'cauchy'
        if ancestry.lower() == 'meta':
            sub_ht = sub_ht.drop('ancestry')
            explode_ht = sub_ht.explode('ancestries', name = 'ancestry')
            explode_ht.export(f'{SUMMARY_ROOT}/aou_v8_meta_gene_signals_{p_field}_{maxMAF}{qc_label}_exploded.txt.bgz')
        sub_ht.describe()
        print(f'{SUMMARY_ROOT}/aou_v8_{ancestry}_gene_signals_{p_field}_{maxMAF}{qc_label}.txt.bgz')
        sub_ht.export(f'{SUMMARY_ROOT}/aou_v8_{ancestry}_gene_signals_{p_field}_{maxMAF}{qc_label}.txt.bgz')


def export_meta_component_signals(ancestry: str, p_field: str, overwrite: bool = False):
    mt = hl.read_matrix_table(f'{MT_ROOT}/{ancestry.upper()}_gene.mt')
    pheno_ht = hl.read_table(f'{ANALYSIS_PATH}/data/phenotype/summary/aou_v8_phenotype_meta_info_qced.ht')
    pheno_ht = pheno_ht.key_by('phenoname')
    pheno_ht = pheno_ht.select('category', 'disease_category', 'description').distinct()
    mt = mt.annotate_cols(**pheno_ht[mt.col_key])
    mt = mt.filter_cols(hl.literal(CATEGORIES).contains(mt.category))

    phenotype_outliers = ['PP_928.1', 'PP_932.11']
    phenotype_to_flag = ['3011099', '3046664', '3013861', '3052648', '3045792_3016921']
    mt = mt.filter_cols(
        (~hl.literal(phenotype_outliers).contains(mt.phenoname))
    ) 
    mt = mt.annotate_cols(hq_pheno=~hl.literal(phenotype_to_flag).contains(mt.phenoname))
    mt = mt.filter_cols(mt.hq_pheno)

    mt = mt.annotate_cols(
        N_eff = hl.if_else((mt.n_controls ==0) | hl.is_missing(mt.n_controls), mt.n_cases, (4 * hl.int64(mt.n_cases) * hl.int64(mt.n_controls)) / (mt.n_cases + mt.n_controls)),
        N = hl.if_else((mt.n_controls ==0) | hl.is_missing(mt.n_controls), mt.n_cases, (mt.n_cases + mt.n_controls)),
    )
    mt = mt.annotate_entries(Pvalue_Burden = hl.if_else(mt.Pvalue_Burden > 0.99, 0.99, mt.Pvalue_Burden))
    mt = mt.annotate_entries(CAF = mt.MAC/(2*mt.N))
    mt = mt.annotate_entries(two_pq = 2*mt.CAF*(1 - mt.CAF))
    mt = mt.annotate_entries(weighted_Z = hl.sqrt(mt.N_eff*mt.two_pq) * (-hl.qnorm(mt.Pvalue_Burden/2)) * hl.sign(mt.BETA_Burden))

    mt = mt.select_cols('description', 'category', 'N', 'N_eff')
    mt = mt.select_rows()
    mt = mt.select_entries('weighted_Z', p_field, 'BETA_Burden')
    ht = mt.entries()
    ht.describe()
    ht = ht.select_globals()
    ht = ht.annotate(ancestry = ancestry)
    
    ht = ht.checkpoint(f'{SUMMARY_ROOT}/aou_v8_{ancestry}_{p_field}_meta_component_raw.ht', overwrite=overwrite, _read_if_exists=not overwrite)

    meta_signal_ht = hl.read_table(f'{SUMMARY_ROOT}/aou_v8_META_{p_field}_gene_signals_qced.ht')
    component_ht = ht.filter(hl.is_defined(meta_signal_ht[ht.key]))
    print(f'{SUMMARY_ROOT}/aou_v8_{ancestry}_meta_component_{p_field}_raw.txt.bgz')
    component_ht.export(f'{SUMMARY_ROOT}/aou_v8_{ancestry}_meta_component_{p_field}_raw.txt.bgz')


def export_meta_skat_burden_comparison(ancestry: str, qc: bool = True, overwrite: bool = False):
    mt = hl.read_matrix_table(f'{MT_ROOT}/META_gene.mt')

    qc_label = '_qced' if qc else '_raw'
    if qc:
        mt = mt.filter_rows(mt.hq_gene)
        mt = mt.filter_cols(mt.hq_pheno)
        mt = mt.filter_entries(mt.hq_exp_CAC)
    
    pheno_ht = hl.read_table('gs://aou_analysis/v8/data/phenotype/summary/aou_v8_phenotype_meta_info_qced.ht')
    pheno_ht = pheno_ht.key_by('phenoname')
    pheno_ht = pheno_ht.select('category', 'disease_category', 'description').distinct()
    mt = mt.annotate_cols(**pheno_ht[mt.col_key])
    mt = mt.filter_cols(hl.literal(CATEGORIES).contains(mt.category))

    phenotype_outliers = ['PP_928.1', 'PP_932.11']
    phenotype_to_flag = ['3011099', '3046664', '3013861', '3052648', '3045792_3016921']
    mt = mt.filter_cols(
        (~hl.literal(phenotype_outliers).contains(mt.phenoname))
    ) 
    mt = mt.annotate_cols(hq_pheno=~hl.literal(phenotype_to_flag).contains(mt.phenoname))
    mt = mt.filter_cols(mt.hq_pheno)

    if ancestry.lower() == 'meta':
        mt = mt.annotate_entries(Pvalue_Burden = mt[META_P_FIELDS['Pvalue_Burden']],
                                 Pvalue_SKAT = mt[META_P_FIELDS['Pvalue_SKAT']])
    mt = mt.filter_rows((hl.literal(['pLoF', 'missenseLC']).contains(mt.annotation)))
    mt.describe()
    if ancestry.lower() == 'meta':
        mt = mt.select_cols('category',
                            n_cases = hl.sum(mt.pheno_data.n_cases),
                            n_controls = hl.sum(mt.pheno_data.n_controls),
                            )
    else:
        mt = mt.select_cols('description', 'category')
    mt = mt.select_rows()
    mt = mt.select_entries('Pvalue_Burden', 'Pvalue_SKAT')

    ht = mt.entries()
    ht = ht.select_globals()
  
    ht = ht.checkpoint(f'{SUMMARY_ROOT}/aou_v8_{ancestry}_skat_burden_gene_signals{qc_label}.ht', overwrite=overwrite, _read_if_exists=not overwrite)
    sub_ht = ht.filter((ht.Pvalue_Burden < 0.0001) | (ht.Pvalue_SKAT < 0.0001))
    gene_list = list(sub_ht.aggregate(hl.agg.collect_as_set(sub_ht.gene_id)))
    ht = ht.filter((hl.literal(gene_list).contains(ht.gene_id))) 
    for maxMAF in [0.001, 0.0001, 0.01, -1]:
        sub_ht = ht.filter((ht.max_MAF == maxMAF))
        maxMAF = maxMAF if maxMAF != -1 else 'cauchy'
        sub_ht.describe()
        print(f'{SUMMARY_ROOT}/aou_v8_{ancestry}_gene_signals_skat_burden_{maxMAF}{qc_label}.txt.bgz')
        sub_ht.export(f'{SUMMARY_ROOT}/aou_v8_{ancestry}_gene_signals_skat_burden_{maxMAF}{qc_label}.txt.bgz')

def export_prop_assoc_data(ancestry, test_type):
    ANNOTATIONS = ['pLoF', 'missenseLC', 'synonymous']
    CATEGORIES = ['lab_measurement', 'pfhh_survey', 'mhwb_survey', 'mcc2_phecodex', 'physical_measurement', 'r_drug']
    threshold = {'genome_variant': 5e-8, 'exome_variant': 5e-8, 'gene': 6.7e-7}
    P_VALUE_FIELDS['gene'] = 'Pvalue_Burden'
    mt = hl.read_matrix_table(f'gs://aou_analysis/v8/mt_results/{ancestry.upper()}_{test_type}.mt')
    mt = mt.filter_cols(~hl.literal(['PP_928.1', 'PP_932.11']).contains(mt.phenoname))
    if ancestry == 'meta' and test_type == 'gene':
        mt = mt.annotate_entries(**{P_VALUE_FIELDS[test_type]: mt.META_Pvalue_Burden})
    if ancestry == 'meta' and test_type != 'gene':
        mt = mt.annotate_entries(**{P_VALUE_FIELDS[test_type]: mt.meta_analysis[0].Pvalue})
    if test_type == 'gene':
        mt = mt.filter_rows((mt.max_MAF == 0.001)& (hl.literal(ANNOTATIONS).contains(mt.annotation)))
    pheno_ht = hl.read_table('gs://aou_analysis/v8/data/phenotype/summary/aou_v8_phenotype_meta_info_qced.ht')
    pheno_ht = pheno_ht.key_by('phenoname')
    pheno_ht = pheno_ht.select('category', 'disease_category', 'description').distinct()
    mt = mt.annotate_cols(**pheno_ht[mt.phenoname])
    mt = mt.filter_cols(hl.literal(CATEGORIES).contains(mt.category))
    
    if test_type == 'gene':
      mt = mt.select_rows()
      mt = mt.annotate_rows(n_assoc_per_gene_burden = hl.agg.count_where(mt[P_VALUE_FIELDS[test_type]] < threshold[test_type]),
                            p_assoc_per_gene_burden = hl.agg.fraction(mt[P_VALUE_FIELDS[test_type]] < threshold[test_type]),
                            n_assoc_per_gene_burden_cat = hl.agg.group_by(mt.category,
                                                                      hl.agg.count_where(mt[P_VALUE_FIELDS[test_type]] < threshold[test_type])),
                            p_assoc_per_gene_burden_cat = hl.agg.group_by(mt.category,
                                                                           hl.agg.fraction(mt[P_VALUE_FIELDS[test_type]] < threshold[test_type])),
                            )
      mt= mt.transmute_rows(
          **{f'n_assoc_per_gene_burden_{category}': mt.n_assoc_per_gene_burden_cat.get(category) for category in CATEGORIES},
          **{f'p_assoc_per_gene_burden_{category}': mt.p_assoc_per_gene_burden_cat.get(category) for category in CATEGORIES})
  
      gene_ht = mt.rows()
      gene_ht = gene_ht.annotate(ancestry = ancestry)
      gene_ht.describe()
      gene_ht.export(f'{SUMMARY_ROOT}/aou_v8_{ancestry}_{test_type}_prop_assoc.txt.bgz')

    mt = mt.select_cols('description', 'category')
    mt = mt.annotate_cols(n_assoc_per_pheno= hl.agg.count_where(mt[P_VALUE_FIELDS[test_type]] < threshold[test_type]),
                          p_assoc_per_pheno = hl.agg.fraction(mt[P_VALUE_FIELDS[test_type]] < threshold[test_type]),
                          )
    if test_type == 'gene': 
        mt= mt.annotate_cols(
          n_assoc_per_pheno_annt = hl.agg.group_by(mt.annotation,
          hl.agg.count_where(mt[P_VALUE_FIELDS[test_type]] < threshold[test_type])),
          p_assoc_per_pheno_annt = hl.agg.group_by(mt.annotation,
          hl.agg.fraction(mt[P_VALUE_FIELDS[test_type]] < threshold[test_type])),

        )
        mt= mt.transmute_cols(
          **{f'n_assoc_per_pheno_{anno}': mt.n_assoc_per_pheno_annt.get(anno) for anno in ANNOTATIONS},
          **{f'p_assoc_per_pheno_{anno}': mt.p_assoc_per_pheno_annt.get(anno) for anno in ANNOTATIONS
        })
    pheno_ht = mt.cols()
    pheno_ht = pheno_ht.annotate(ancestry = ancestry)
    pheno_ht.describe()
    pheno_ht.export(f'{SUMMARY_ROOT}/aou_v8_{ancestry}_{test_type}_pheno_prop_assoc.txt.bgz')

def main(args):
    hl.init(default_reference="GRCh38")

    if args.summarize_associations:
        for ancestry in ANCESTRIES:
            if not hl.hadoop_exists(f'{SUMMARY_ROOT}/aou_v8_{ancestry}_ACAF_signals_qced.ht/_SUCCESS') or args.overwrite:
                export_variant_signals(ancestry, 'ACAF', qc=True, overwrite=args.overwrite)
            if not hl.hadoop_exists(f'{SUMMARY_ROOT}/aou_v8_{ancestry}_Exome_signals_qced.ht/_SUCCESS') or args.overwrite:
                export_variant_signals(ancestry, 'Exome', qc=True, overwrite=args.overwrite)
        for test in ['burden', 'skato', 'skat']:
            for ancestry in ANCESTRIES:
                p_field = P_VALUE_FIELDS[test]
                OVERWRITE = not hl.hadoop_exists(f'{SUMMARY_ROOT}/aou_v8_{ancestry}_{p_field}_gene_signals_qced.ht/_SUCCESS') or args.overwrite
                export_gene_signals(ancestry, p_field=p_field, qc=True, overwrite=OVERWRITE)
    
    if args.export_meta_components:
        for ancestry in tqdm(ANCESTRIES):
            if ancestry == 'META':
                continue
            p_field = P_VALUE_FIELDS['burden']
            OVERWRITE = not hl.hadoop_exists(f'{SUMMARY_ROOT}/aou_v8_{ancestry}_{p_field}_meta_component_raw.ht/_SUCCESS') or args.overwrite
            export_meta_component_signals(ancestry, p_field=p_field, overwrite=OVERWRITE)

    if args.compare_skat_burden:
        for ancestry in ['META']:
            export_meta_skat_burden_comparison(ancestry, qc=False, overwrite=args.overwrite)

    if args.export_prop_assoc:
        for anc in ['meta']:
            for test in ['gene','genome_variant', 'exome_variant']:
                if not hl.hadoop_exists(f'{SUMMARY_ROOT}/aou_v8_{anc}_{test}_pheno_prop_assoc.txt.bgz'):
                    print(f'------{anc.upper()}-------{test}----------')
                    export_prop_assoc_data(ancestry = anc, test_type = test)

    if args.compute_phenotype_correlation: 
        def process_pheno_mt(pheno_category):
            mt = hl.read_matrix_table(f'{ANALYSIS_PATH}/data/phenotype/{pheno_category}_annotated.mt')
            mt = mt.filter_rows(mt.hard_filter)
            mt = mt.annotate_entries(value = hl.float64(mt.value))
            mt = mt.select_rows().select_cols()
            return mt
        
        def make_pairwise_ht(mt: hl.MatrixTable, pheno_field, min_cases: int = 500, correlation: bool = False):
            mt = mt.annotate_entries(_pheno=pheno_field)
            mt = mt.add_col_index()
            index_ht = mt.cols().key_by('col_idx')
            if correlation:
                bm = hl.linalg.BlockMatrix.from_entry_expr(mt._pheno, mean_impute=True, center=True, normalize=True, axis='cols', block_size=1024)
            else:
                bm = hl.linalg.BlockMatrix.from_entry_expr(mt._pheno, block_size=1024)
            bm = bm.T @ bm
            pheno_ht = bm.entries()
            pheno_ht = pheno_ht.annotate(i_data=index_ht[pheno_ht.i], j_data=index_ht[pheno_ht.j])
            if not correlation:
                pheno_ht = pheno_ht.annotate(prop_overlap=pheno_ht.entry / pheno_ht.i_data.n_cases)
            return pheno_ht
        
        def more_cases_tie_breaker(l, r):
            """
            Tie breaker function prefering phenotypes with more cases used in hl.maximal_independent_set()
            :param l: one phenotype
            :param r: the other phenotype
            :return: integer indicating which phenotype is preferred
            :rtype: int
            """
            return (
                hl.case()
                .when(l.n_cases_defined > r.n_cases_defined, -1)
                .when(l.n_cases_defined == r.n_cases_defined, 0)
                .when(l.n_cases_defined < r.n_cases_defined, 1)
                .or_missing()
            )

        ht = hl.read_table(f'{ANALYSIS_PATH}/data/phenotype/lab_measurement/lab_measurement_full_annotated.ht')
        ht.describe()
        print(ht.count())
        print(ht.aggregate(hl.agg.counter(ht.measurement_concept_id)))
        ht = ht.filter(ht.hard_filter)
        ht = ht.select(phenoname = ht.measurement_concept_id, value = ht.median_value)
        lab_mt = ht.to_matrix_table(row_key=['person_id'], col_key=['phenoname'])

        full_mt = process_pheno_mt('mcc2_phecodex').union_cols(process_pheno_mt('r_drug'), row_join_type='outer').union_cols(process_pheno_mt('pfhh_survey'), row_join_type='outer').union_cols(process_pheno_mt('mhwb_survey'), row_join_type='outer').union_cols(process_pheno_mt('physical_measurement'), row_join_type='outer').union_cols(lab_mt, row_join_type='outer')
        full_mt.describe()
        full_mt = full_mt.checkpoint(f'{TMP_BUCKET}/full_pheno_v8.mt', overwrite=True)
        full_mt = full_mt.annotate_cols(n_cases_defined = hl.agg.count_where(hl.is_defined(full_mt.value)))

        corr = make_pairwise_ht(full_mt, pheno_field=full_mt.value, correlation=True)
        corr = corr.checkpoint('gs://aou_wlu/v8_analysis/pheno_corr.ht', overwrite = True)
        corr.describe()

        corr_ht = hl.read_table('gs://aou_wlu/v8_analysis/pheno_corr.ht')
        corr_ht = corr_ht.key_by(
            i_pheno = corr_ht.i_data.phenoname,
            j_pheno = corr_ht.j_data.phenoname,
        )
        corr_ht = corr_ht.filter(corr_ht.i < corr_ht.j)
        corr_ht = corr_ht.select('entry')
        corr_ht.describe()
        corr_ht.export('gs://aou_wlu/v8_analysis/pheno_corr.txt.bgz')

        r2_threshold = 0.5
        related = corr_ht.filter((corr_ht.entry ** 2 >= r2_threshold) & (corr_ht.i != corr_ht.j))
        pheno_to_remove = hl.maximal_independent_set(
            related.i_data, related.j_data, keep=False, tie_breaker=more_cases_tie_breaker
        )

        pheno_to_remove = pheno_to_remove.annotate(
            phenoname = pheno_to_remove.node.phenoname,
            n_cases_defined = pheno_to_remove.node.n_cases_defined
        )
        pheno_to_remove = pheno_to_remove.key_by('phenoname')
        pheno_to_remove = pheno_to_remove.drop('node')
        pheno_to_remove.describe()

        pheno_to_remove = pheno_to_remove.checkpoint('gs://aou_wlu/v8_analysis/phenos_to_remove_r2_0.5.ht', _read_if_exists = True)
        print(pheno_to_remove.count())
        pheno_to_remove.export('gs://aou_wlu/v8_analysis/phenos_to_remove_r2_0.5.txt.bgz')

    if args.extract_phenotypes:
        phecodex_mt = hl.read_matrix_table('gs://aou_analysis/v8/data/phenotype/mcc2_phecodex_annotated.mt')
        drug_mt = hl.read_matrix_table('gs://aou_analysis/v8/data/phenotype/r_drug_annotated.mt')
        phecodex_mt = phecodex_mt.filter_cols(hl.literal(['EM_202.2', 'EM_236.1']).contains(phecodex_mt.phenoname))
        phecodex_mt = phecodex_mt.filter_rows(phecodex_mt.hard_filter) 
        phecodex_ht = phecodex_mt.entries()
        phecodex_ht = phecodex_ht.select('ancestry', 'sex', 'value')
        phecodex_ht.export('gs://aou_wlu/v8_analysis/t2d_obesity_phenotype_by_sex.txt.bgz')

        drug_mt = drug_mt.filter_cols(hl.literal(['A10BJ']).contains(drug_mt.phenoname))
        drug_mt = drug_mt.filter_rows(drug_mt.hard_filter) 
        drug_ht = drug_mt.entries()
        drug_ht = drug_ht.select('ancestry', 'sex', 'value')
        drug_ht.export('gs://aou_wlu/v8_analysis/a10bj_phenotype_by_sex.txt.bgz')

    if args.compare_pfhh_phecodex:
        pfhh_map_info = hl.import_table('gs://aou_wlu/v8_analysis/ehr_vs_survey/pfhh_phecodeX_mapping_simple.csv', delimiter= ',') # Note: manually mapped
        pfhh_map_info.describe()
        pfhh_map_info.show()
        print(pfhh_map_info.count())
        pfhh_map_info = pfhh_map_info.transmute(phenoname = pfhh_map_info['\ufeffphenoname'])
        pfhh_map_info = pfhh_map_info.checkpoint('gs://aou_wlu/v8_analysis/ehr_vs_survey/pfhh_phecodeX_mapping_simple.ht')

        pfhh_map_info = hl.read_table('gs://aou_wlu/v8_analysis/ehr_vs_survey/pfhh_phecodeX_mapping_simple.ht')
        pfhh_map_info = pfhh_map_info.key_by('phenoname')
        map_dict = hl.zip(pfhh_map_info.aggregate(hl.agg.collect(pfhh_map_info.phenoname)), 
                        pfhh_map_info.aggregate(hl.agg.collect(pfhh_map_info.phecodeX)))
        map_dict = hl.dict(map_dict)

        meta_gene_mt = hl.read_matrix_table('gs://aou_analysis/v8/mt_results/META_gene_CAF.mt')
        meta_gene_mt = meta_gene_mt.select_rows().select_cols().select_entries('META_Pvalue_Burden')
        pfhh_mt = meta_gene_mt.filter_cols(hl.literal(hl.eval(map_dict.keys())).contains(meta_gene_mt.phenoname))
        pfhh_ht = pfhh_mt.entries()
        pfhh_ht = pfhh_ht.annotate(phecodeX = pfhh_map_info[pfhh_ht.phenoname].phecodeX)
        pfhh_ht = pfhh_ht.key_by('gene_id', 'gene_symbol', 'annotation', 'max_MAF', 'phecodeX')
        phecodeX_mt = meta_gene_mt.filter_cols(hl.literal(hl.eval(map_dict.values())).contains(meta_gene_mt.phenoname))
        phecodeX_ht = phecodeX_mt.entries()
        phecodeX_ht = phecodeX_ht.join(pfhh_ht)

        phecodeX_ht = phecodeX_ht.filter((phecodeX_ht.META_Pvalue_Burden < 0.05) | 
                                            (phecodeX_ht.META_Pvalue_Burden_1 < 0.05))
        phecodeX_ht = phecodeX_ht.checkpoint('gs://aou_wlu/v8_analysis/phenotype/pfhh_phecodex_signals.ht', _read_if_exists = True)
        phecodeX_ht.show()
        phecodeX_ht.export('gs://aou_wlu/v8_analysis/phenotype/pfhh_phecodex_signals.txt.bgz')

        # Define all fields to collect
        fields = [
            'pfhh_code', 'phecodex',
            'pfhh_prevalence', 'pfhh_n', 'pfhh_lambda_skato', 'pfhh_lambda_skat', 'pfhh_lambda_burden',
            'phecodex_prevalence', 'phecodex_n', 'phecodex_lambda_skato', 'phecodex_lambda_skat', 'phecodex_lambda_burden',
            'p_burden_corr', 'p_skato_corr', 'p_skat_corr',
            'p_burden_corr_pLoF', 'p_skato_corr_pLoF', 'p_skat_corr_pLoF',
            'p_burden_corr_missense', 'p_skato_corr_missense', 'p_skat_corr_missense',
            'p_burden_corr_synonymous', 'p_skato_corr_synonymous', 'p_skat_corr_synonymous',
            'stats_burden_corr', 'stats_skato_corr', 'stats_skat_corr',
            'stats_burden_corr_pLoF', 'stats_skato_corr_pLoF', 'stats_skat_corr_pLoF',
            'stats_burden_corr_missense', 'stats_skato_corr_missense', 'stats_skat_corr_missense',
            'stats_burden_corr_synonymous', 'stats_skato_corr_synonymous', 'stats_skat_corr_synonymous',
        ]
        data = {field: [] for field in fields}

        base_path = 'gs://aou_analysis/v8/ht_results/META/phenotype_{}/gene_results_exact_expected_p_0.001.ht'
        map_dict_local = hl.eval(map_dict)

        for pfhh_code in hl.eval(map_dict.keys()):
            phecodex = map_dict_local[pfhh_code]
            pfhh_path = base_path.format(pfhh_code)
            phecodex_path = base_path.format(phecodex)

            if not hl.hadoop_exists(f'{pfhh_path}/_SUCCESS') or not hl.hadoop_exists(f'{phecodex_path}/_SUCCESS'):
                continue

            print(f'{pfhh_code} -> {phecodex}')
            data['pfhh_code'].append(pfhh_code)
            data['phecodex'].append(phecodex)

            ht1 = hl.read_table(pfhh_path)
            ht2 = hl.read_table(phecodex_path)
            ht = ht1.join(ht2)

            # PFHH summary
            data['pfhh_prevalence'].append(hl.eval(ht1.N_cases) / hl.eval(ht1.N))
            data['pfhh_n'].append(hl.eval(ht1.N))
            lambda_gc_1 = hl.eval(ht1['lambda_gc_maxmaf_0.001'])
            data['pfhh_lambda_skato'].append(lambda_gc_1['lambda_gc_Pvalue'])
            data['pfhh_lambda_skat'].append(lambda_gc_1['lambda_gc_Pvalue_SKAT'])
            data['pfhh_lambda_burden'].append(lambda_gc_1['lambda_gc_Pvalue_Burden'])

            # PhecodeX summary
            data['phecodex_prevalence'].append(hl.eval(ht2.N_cases) / hl.eval(ht2.N))
            data['phecodex_n'].append(hl.eval(ht2.N))
            lambda_gc_2 = hl.eval(ht2['lambda_gc_maxmaf_0.001'])
            data['phecodex_lambda_skato'].append(lambda_gc_2['lambda_gc_Pvalue'])
            data['phecodex_lambda_skat'].append(lambda_gc_2['lambda_gc_Pvalue_SKAT'])
            data['phecodex_lambda_burden'].append(lambda_gc_2['lambda_gc_Pvalue_Burden'])

            # P-value correlations (compute group_by once per metric)
            data['p_burden_corr'].append(ht.aggregate(hl.agg.corr(ht.Pvalue_Burden, ht.Pvalue_Burden_1)))
            data['p_skato_corr'].append(ht.aggregate(hl.agg.corr(ht.Pvalue, ht.Pvalue_1)))
            data['p_skat_corr'].append(ht.aggregate(hl.agg.corr(ht.Pvalue_SKAT, ht.Pvalue_SKAT_1)))

            p_burden_by_annot = ht.aggregate(hl.agg.group_by(ht.annotation, hl.agg.corr(ht.Pvalue_Burden, ht.Pvalue_Burden_1)))
            p_skato_by_annot = ht.aggregate(hl.agg.group_by(ht.annotation, hl.agg.corr(ht.Pvalue, ht.Pvalue_1)))
            p_skat_by_annot = ht.aggregate(hl.agg.group_by(ht.annotation, hl.agg.corr(ht.Pvalue_SKAT, ht.Pvalue_SKAT_1)))
            for annot, suffix in [('pLoF', 'pLoF'), ('missenseLC', 'missense'), ('synonymous', 'synonymous')]:
                data[f'p_burden_corr_{suffix}'].append(p_burden_by_annot.get(annot))
                data[f'p_skato_corr_{suffix}'].append(p_skato_by_annot.get(annot))
                data[f'p_skat_corr_{suffix}'].append(p_skat_by_annot.get(annot))

            # Meta-stat correlations
            data['stats_burden_corr'].append(ht.aggregate(hl.agg.corr(ht.META_Stats_Burden, ht.META_Stats_Burden_1)))
            data['stats_skato_corr'].append(ht.aggregate(hl.agg.corr(ht.META_Stats_SKATO, ht.META_Stats_SKATO_1)))
            data['stats_skat_corr'].append(ht.aggregate(hl.agg.corr(ht.META_Stats_SKAT, ht.META_Stats_SKAT_1)))

            stats_burden_by_annot = ht.aggregate(hl.agg.group_by(ht.annotation, hl.agg.corr(ht.META_Stats_Burden, ht.META_Stats_Burden_1)))
            stats_skato_by_annot = ht.aggregate(hl.agg.group_by(ht.annotation, hl.agg.corr(ht.META_Stats_SKATO, ht.META_Stats_SKATO_1)))
            stats_skat_by_annot = ht.aggregate(hl.agg.group_by(ht.annotation, hl.agg.corr(ht.META_Stats_SKAT, ht.META_Stats_SKAT_1)))
            for annot, suffix in [('pLoF', 'pLoF'), ('missenseLC', 'missense'), ('synonymous', 'synonymous')]:
                data[f'stats_burden_corr_{suffix}'].append(stats_burden_by_annot.get(annot))
                data[f'stats_skato_corr_{suffix}'].append(stats_skato_by_annot.get(annot))
                data[f'stats_skat_corr_{suffix}'].append(stats_skat_by_annot.get(annot))

        # Build rows from collected data
        rows = [{field: data[field][i] for field in fields} for i in range(len(data['pfhh_code']))]

        ht_out = hl.Table.parallelize(rows, key=['pfhh_code', 'phecodex'])

        # Peek and/or write
        ht_out.show()
        ht_out.write('gs://aou_wlu/v8_analysis/phenotype/pfhh_phecodex_pairwise_metrics_gene.ht', overwrite=True)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--overwrite",
        help="Whether to overwrite existing files",
        action="store_true"
    )
    parser.add_argument(
        "--summarize-associations",
        help="Whether to summarize associations",
        action="store_true"
    )
    parser.add_argument(
        "--genebass-comparison",
        help="Whether to perform genebass comparison",
        action="store_true"
    )
    parser.add_argument(
        "--export-meta-components",
        help="Whether to export meta components for each ancestry and test type",
        action="store_true"
    )
    parser.add_argument(
        "--reformat-genebass-mt",
        help="Whether to reformat genebass mt",
        action="store_true"
    )
    parser.add_argument(
        "--reformat-aou-mt",
        help="Whether to reformat aou mt",
        action="store_true"
    )
    parser.add_argument(
        "--redo-merge",
        help="Whether to redo merge of aou and ukb mts",
        action="store_true"
    )
    parser.add_argument(
        "--compare-skat-burden",
        help="Whether to compare skat and burden results",
        action="store_true"
    )
    parser.add_argument(
        "--export-prop-assoc",
        help="Whether to export proportion of associated phenotypes/genes",
        action="store_true"
    )
    parser.add_argument(
        "--compute-phenotype-correlation",
        help="Whether to compute phenotype correlation and prune correlated phenotypes",
        action="store_true"
    )
    parser.add_argument(
        "--extract-phenotypes",
        help="Whether to extract specific phenotypes",
        action="store_true"
    )
    parser.add_argument(
        "--compare-pfhh-phecodex",
        help="Whether to compare PFHH and PhecodeX phenotypes",
        action="store_true"
    )
    args = parser.parse_args()

    main(args)