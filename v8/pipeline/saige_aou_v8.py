#!/usr/bin/env python3

__author__ = "wlu"

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
from typing import *
from tqdm import tqdm
import warnings

def range_table(log_file):
    import hail as hl
    hl.init(log=log_file)
    print(hl.utils.range_table(10)._force_count())


TRANCHE = "v8"
MY_BUCKET = 'gs://aou_wlu'
ANALYSIS_BUCKET = "gs://aou_analysis/v8"
DATA_PATH = f'{ANALYSIS_BUCKET}/data'
TMP_BUCKET = 'gs://aou_tmp/v8'

N_GENE_PER_GROUP = 100
CHUNK_SIZE = {'all': int(1.25e6),'eur': int(6.25e6), "afr": int(1.25e7), "amr": int(1.25e7), "eas": int(1.25e7), "mid": int(1.25e7), "sas": int(1.25e7)}
# old_CHUNK_SIZE = {'all': int(1.25e6),'eur': int(1.25e7), "afr": int(2.5e7), "amr": int(2.5e7), "eas": int(2.5e7), "mid": int(2.5e7), "sas": int(2.5e7)}
REFERENCE = "GRCh38"
CHROMOSOMES = list(map(str, range(1, 23))) + ["X", "Y"]

PILOT_PHENOTYPES = ["height", "heart-rate-mean", "A10BJ", "A10BJ06", "random_0.5_continuous_1", "random_0.5_0.01_1", "random_0.5_0.5_1", "random_0.5_0.2_1", 
                    "random_0.5_0.1_1", "random_0.5_0.001_1"]
PILOT_PHENOTYPES = PILOT_PHENOTYPES + [f'{pheno}_male' for pheno in PILOT_PHENOTYPES] + [f'{pheno}_female' for pheno in PILOT_PHENOTYPES] + \
                    ["Birth_to_PULMHEART_amr_ALL", "Birth_to_CAD_eas_ALL", "Birth_to_CAD_sas_ALL", "Birth_to_MI_eur_ALL", "Birth_to_DEMENTIA_eur_ALL",
                    "Birth_to_ATHSCLE_eur_ALL", "HYPTENSESS_to_MI_eur_ALL", "PARKINSON_to_DEMENTIA_eur_ALL", "T2D_to_ATHSCLE_eur_ALL", 'Birth_to_LUNGCA_afr_ALL']
PHENO_CATEGORIES = ['physical_measurement', 'r_drug', 'pfhh_survey', 'random_pheno', 'lab_measurement', 'mcc2_phecode', 'mcc2_phecodex', 'onset', 'progression']
PHENO_CATEGORIES_MT = ['physical_measurement', 'r_drug', 'pfhh_survey', 'mcc2_phecode', 'mcc2_phecodex']
PHENO_CATEGORIES_HT = ['lab_measurement', 'onset', 'progression']
GATE_CATEGORIES = ['onset', 'progression']
BINARY_CATEGORIES = ['r_drug', 'pfhh_survey', 'mcc2_phecode', 'mcc2_phecodex']
QUANTITATIVE_CATEGORIES = ['physical_measurement', 'lab_measurement']
ANCESTRIES = ['AFR', 'AMR', 'EAS', 'EUR', 'MID', 'SAS']    
N_SAMPLES_PRUNED = {'afr': 77444, 'amr': 71540, 'eas': 9488, 'eur': 227273, 'mid': 1153, 'sas': 5132, 'all': 0}


################################################ Copy-paste constants.py - Remove after docker image is built ################################################
def get_adj_expr(
   gt_expr: hl.expr.CallExpression,
   gq_expr: Union[hl.expr.Int32Expression, hl.expr.Int64Expression],
   ad_expr: hl.expr.ArrayNumericExpression,
   adj_gq: int = 30,
   adj_ab: float = 0.2,
) -> hl.expr.BooleanExpression:
   """
   Get adj genotype annotation.


   Defaults correspond to gnomAD values.
   """
   return (
       (gq_expr >= adj_gq)
       & (
           hl.case()
           .when(~gt_expr.is_het(), True)
           .when(gt_expr.is_het_ref(), ad_expr[gt_expr[1]] / hl.sum(ad_expr) >= adj_ab)
           .default(
               (ad_expr[gt_expr[0]] / hl.sum(ad_expr) >= adj_ab)
               & (ad_expr[gt_expr[1]] / hl.sum(ad_expr) >= adj_ab)
           )
       )
   )

def annotate_adj(
       mt: hl.MatrixTable,
       adj_gq: int = 30,
       adj_ab: float = 0.2,
) -> hl.MatrixTable:
   """
   Annotate genotypes with adj criteria (assumes diploid).


   Defaults correspond to gnomAD values.
   """
   if "GT" not in mt.entry and "LGT" in mt.entry:
       print("No GT field found, using LGT instead.")
       gt_expr = mt.LGT
   else:
       gt_expr = mt.GT


   if "AD" not in mt.entry and "LAD" in mt.entry:
       print("No AD field found, using LAD instead.")
       ad_expr = mt.LAD
   else:
       ad_expr = mt.AD


   return mt.annotate_entries(
       adj=get_adj_expr(
           gt_expr, mt.GQ, ad_expr, adj_gq, adj_ab
       )
   )

def get_filtered_mt(mt_type: hl.tstr,
                    sample_ids: hl.Table,
                    ancestry: str='all',
                    filter_samples: bool=True, 
                    filter_variants: bool=True,
                    prune_samples:bool=True,
                    adj_filter: bool=True):
    anc_mt_path = f'{DATA_PATH}/utils/raw_mt/{mt_type}/{ancestry.upper()}_{mt_type}.mt'
    print(anc_mt_path)
    if ancestry != 'all':
        ancestry_ht = hl.read_table(f'{DATA_PATH}/utils/aou_v8_global_pca.ht')
        ancestry_ht = ancestry_ht.filter(ancestry_ht.ancestry_pred == ancestry)
        if not hfs.exists(f'{anc_mt_path}/_SUCCESS'):
            mt_path = ACAF_MT_PATH if mt_type=='ACAF' else EXOME_MT_PATH
            print(mt_path)
            mt = hl.read_matrix_table(mt_path)
            mt.describe()
            print(mt.count())
            mt = mt.filter_entries(hl.is_missing(mt.FT) | (mt.FT == 'PASS'))
            mt = mt.filter_cols(hl.is_defined(ancestry_ht[mt.col_key]))
            mt = mt.naive_coalesce(20000).checkpoint(anc_mt_path, overwrite=True)
        mt = hl.read_matrix_table(anc_mt_path)
        print(mt._force_count_cols())
    else:
        mt_path = ACAF_MT_PATH if mt_type=='ACAF' else EXOME_MT_PATH
        mt = hl.read_matrix_table(mt_path)
        mt = mt.filter_entries(hl.is_missing(mt.FT) | (mt.FT == 'PASS'))
        mt = mt.naive_coalesce(20000).checkpoint(anc_mt_path, overwrite=True)
        print(mt.count())

    if filter_variants:
        print(f'Filtering to global AC > 0...')
        mt = mt.filter_rows(
            (mt.info.AC[0] > 0)
        )
        
    if filter_samples:
        print(f'Filtering to samples with sex, ancestry defined and age <= 100...')
        mt = mt.filter_cols(hl.is_defined(sample_ids[mt.col_key]))
        
    if prune_samples:
        print(f'[{mt_type}] Filtering to samples pruned from the PCA centroid pipeline...')
        pca_pruned_sample_tsv = f'{DATA_PATH}/utils/pca/results/aou_{ancestry.lower()}_centroid_pruned.tsv'
        pca_pruned_sample_path = f'{DATA_PATH}/utils/pca/results/{ancestry.lower()}_pca_centroid_pruned.ht'
        overwrite_pruned_ht = False
        if hfs.exists(pca_pruned_sample_tsv):
            if not hfs.exists(f'{pca_pruned_sample_path}/_SUCCESS') or overwrite_pruned_ht:
                pruned_ht = hl.import_table(pca_pruned_sample_tsv, delimiter='\t', key='s', types = {'s':hl.tstr})
                pruned_ht = pruned_ht.checkpoint(pca_pruned_sample_path, overwrite=overwrite_pruned_ht)
            else: 
                pruned_ht = hl.read_table(pca_pruned_sample_path)
        print(pruned_ht.count())
        pruned_ht = hl.read_table(pca_pruned_sample_path)
        mt = mt.filter_cols(hl.is_defined(pruned_ht[mt.col_key]))
        
    if adj_filter:
        """
       Filter genotypes to adj criteria - Default:
       GQ >= 30, het_ref_altAB >= 0.2, het_non_ref_altAB >= 0.2 , het_non_ref_refAB >= 0.2
       """
        print(f'Applying adj filters...')
        mt = annotate_adj(mt)
        mt = mt.filter_entries(mt.adj)
        mt = mt.filter_rows(hl.agg.any(mt.GT.n_alt_alleles() >0))
    
    return mt

################################################ Remove after docker image is built ################################################

logging.basicConfig(
    format="%(levelname)s (%(name)s %(lineno)s): %(message)s",
    level="INFO",
    filename="saige_pipeline.log",
)
logger = logging.getLogger("ALL_x_AoU_SAIGE")
logger.setLevel(logging.INFO)

HAIL_DOCKER_IMAGE = "hailgenetics/hail:0.2.133-py3.11"
SAIGE_DOCKER_IMAGE = "wzhou88/saige:1.4.4"  # latest
QQ_DOCKER_IMAGE = "konradjk/saige_qq:0.2"

def annotate_expected_pvalue(ht: hl.Table, method:str, p_field: str='Pvalue', k:int = 5000):
    n = ht.count()
    print(n)
    ht = ht.filter(hl.is_defined(ht[p_field]))
    keys = list(ht.key)
    if method == 'exact':
        ht = ht.key_by()
        ht = ht.order_by(ht[p_field]).add_index()
        ht = ht.key_by(*keys)
        ht = ht.annotate(**{f'{p_field}_expected': ht.idx / (n + 1)})
        if n % 2 == 1:
            median = [int(n * 0.5), int(n * 0.5) + 1]
            median_p = ht.filter((ht.idx == median[0]) | (ht.idx == median[1]) )
        else:
            median_p = ht.filter(ht.idx == int(n * 0.5))

        median_p = median_p.aggregate(hl.agg.mean(median_p[p_field]))
        ht = ht.drop('idx')

    elif method == 'approx_cdf':
        ht = ht.annotate_globals(ranks=ht.aggregate(hl.agg.approx_cdf(ht[p_field], k=k)))
        ht = ht.annotate(p_rank_upper=hl.enumerate(ht.ranks['values']).find(lambda x: x[1] >= ht[p_field]))
        ht = ht.annotate(p_rank=ht.p_rank_upper[0],
                            upper_pvalue=ht.p_rank_upper[1])
        ht = ht.annotate(lower_pvalue=ht.ranks['values'][ht.p_rank - 1],
                            lower_rank=ht.ranks['ranks'][ht.p_rank - 1],
                            upper_rank=ht.ranks['ranks'][ht.p_rank],
                            )
        ht = ht.annotate(
            rank=hl.int64(hl.floor((ht.upper_rank - ht.lower_rank) / (ht.upper_pvalue - ht.lower_pvalue) * (
                    ht.Pvalue - ht.lower_pvalue) + ht.lower_rank)))
        import time
        ht = ht.checkpoint(f'gs://aou_tmp/{time.time()}.ht')
        top_ht = ht.filter((ht.rank < n * 0.01) | (ht[p_field] < 1e-3))
        rest_ht = ht.filter((ht.rank >= n * 0.01) & (ht[p_field] >= 1e-3))

        top_ht = top_ht.key_by()
        top_ht = top_ht.order_by(top_ht[p_field]).add_index('rank')
        top_ht = top_ht.key_by(*keys)
        ht = top_ht.union(rest_ht)
        ht = ht.annotate(Pvalue_expected=ht.rank / (n + 1))
        ht = ht.drop('p_rank_upper', 'p_rank', 'upper_pvalue', 'lower_pvalue', 'lower_rank', 'upper_rank')

        from bisect import bisect
        median = hl.eval(hl.floor(n * 0.5))
        upper_idx = bisect(hl.eval(ht.ranks[1]), median)
        upper_rank = hl.eval(ht.ranks[1])[upper_idx]
        lower_rank = hl.eval(ht.ranks[1])[upper_idx - 1]
        upper_pvalue = hl.eval(ht.ranks[0])[upper_idx]
        lower_pvalue = hl.eval(ht.ranks[0])[upper_idx - 1]
        median_p = (median - lower_rank) / (upper_rank - lower_rank) * (
                    upper_pvalue - lower_pvalue) + lower_pvalue
        ht = ht.drop('ranks')

    else:
        import sys
        sys.exit("Method should be one of ['exact', 'approx_cdf']")


    lambda_gc = hl.eval(hl.qchisqtail(p=median_p, df=1, ncp=0, lower_tail=False) /
                        hl.qchisqtail(p=0.5, df=1, ncp=0, lower_tail=True))
    print(f"Lambda GC: {lambda_gc}")

    ht = ht.annotate_globals(**{f'lambda_gc_{p_field}': lambda_gc})
    ht = ht.annotate(**{f'{p_field}': hl.if_else(ht[f'{p_field}'] == 0, 1e-320, ht[f'{p_field}'])}, )
    ht = ht.annotate(**{f'{p_field}_expected_log10': -hl.log10(ht[f'{p_field}_expected'])},)

    return ht


def load_gene_data(directory: str,
                   output_ht_directory: str,
                   phenoname, gene_ht_map_path: str, quantitative_trait:bool,
                   null_glmm_log: str, saige_log: str = 'NA',
                   overwrite: bool = False):

    hl.init(
        master='local[8]',
        tmp_dir=TMP_BUCKET,
        default_reference="GRCh38",
    )

    def get_cases_and_controls_from_log(log_format):
        """
        'gs://path/to/result_chr{chrom}_000000001.variant.log'
        """
        cases = controls = -1
        for chrom in range(10, 23):
            try:
                with hfs.open(log_format.format(chrom=chrom)) as f:
                    for line in f:
                        line = line.strip()
                        if line.startswith('Analyzing'):
                            fields = line.split()
                            if len(fields) == 6:
                                try:
                                    cases = int(fields[1])
                                    controls = int(fields[4])
                                    break
                                except ValueError:
                                    logger.warn(f'Could not load number of cases or controls from {line}.')
                        elif line.endswith(
                                'samples were used in fitting the NULL glmm model and are found in sample file') or \
                                line.endswith('samples have been used to fit the glmm null model'):
                            # This is ahead of the case/control count line ("Analyzing ...") above so this should be ok
                            fields = line.split()
                            try:
                                cases = int(fields[0])
                            except ValueError:
                                logger.warn(f'Could not load number of cases or controls from {line}.')
                return cases, controls
            except:
                pass
        return cases, controls

    def get_heritability_from_log(log_file, quantitative_trait: bool = False):
        import math
        heritability = -1
        with hfs.open(log_file) as f:
            for line in f:
                if line.startswith('Final'):
                    fields = line.strip().split()
                    if len(fields) == 4:
                        try:
                            tau = float(fields[2])
                            if quantitative_trait:
                                tau1 = float(fields[1])
                                heritability = tau / (tau1 + tau)
                            else:
                                heritability = tau / (tau + math.pi ** 2 / 3)
                            break
                        except:
                            logger.warn(f'Could not load heritability from {line}.')
        return heritability

    def get_saige_version_from_log(null_glmm_log):
        version = 'NA'
        with hfs.open(null_glmm_log) as f:
            for line in f:
                if line.startswith('other attached packages:'):
                    try:
                        line2 = f.readline()
                        packages = line2.strip().split()
                        version = [x for x in packages if 'SAIGE' in x][0]
                    except:
                        logger.warning(f'Could not load version number from {line2} in {null_glmm_log}.')
        return version

    def get_inverse_normalize_status(null_glmm_log):
        status = 'Unknown'
        with hfs.open(null_glmm_log) as f:
            for line in f:
                if line.startswith('$invNormalize'):
                    try:
                        status = f.readline().strip().split()[1]
                    except:
                        logger.warning(f'Could not load inv_norm status from {line} in {null_glmm_log}.')
        return status.capitalize()

    def parse_log_p_value(x):
        log_x = hl.if_else(x.contains('E'), -hl.int(x.split("E")[1]) - hl.log10(hl.float64(x.split("E")[0])), -hl.log10(hl.float64(x)))
        return log_x

    n_cases, n_controls = get_cases_and_controls_from_log(saige_log)
    heritability = get_heritability_from_log(null_glmm_log,
                                             quantitative_trait) if null_glmm_log else -1.0
    inv_normalized = get_inverse_normalize_status(null_glmm_log) if null_glmm_log else 'NA'
    saige_version = get_saige_version_from_log(null_glmm_log) if null_glmm_log else 'NA'
    output_ht_path = f'{output_ht_directory}/gene_results.ht'
    print(f'Loading: {directory}/*.gene.txt ...')
    types = {x: hl.tint32 for x in ('MAC',	'Number_rare', 'Number_ultra_rare')}
    types.update({x: hl.tfloat64 for x in ('max_MAF', 'BETA_Burden', 'SE_Burden')})
    types.update({x: hl.tstr for x in ('Pvalue', 'Pvalue_Burden', 'Pvalue_SKAT')})
    ht = hl.import_table(f'{directory}/*.gene.txt', delimiter='\t', impute=True, types=types)
    ht.describe()
    if n_cases == -1: n_cases = hl.null(hl.tint)
    if n_controls == -1: n_controls = hl.null(hl.tint)
    if heritability == -1.0: heritability = hl.null(hl.tfloat)
    if saige_version == 'NA': saige_version = hl.null(hl.tstr)
    if inv_normalized == 'NA': inv_normalized = hl.null(hl.tstr)

    ht = ht.annotate(**{f'{p_field}_log10': parse_log_p_value(ht[p_field]) for p_field in ('Pvalue', 'Pvalue_Burden', 'Pvalue_SKAT')})
    ht = ht.annotate(**{f'{p_field}': hl.float64(ht[p_field]) for p_field in ('Pvalue', 'Pvalue_Burden', 'Pvalue_SKAT')})
    ht = ht.filter(hl.len(ht.Region.split('_')) == 2)
    fields = ht.Region.split('_')
    gene_ht = hl.read_table(gene_ht_map_path).select('interval').distinct()
    ht = ht.key_by(gene_id=fields[0], gene_symbol=fields[1], annotation=ht.Group, phenoname=phenoname).drop('Region', 'Group').naive_coalesce(10).annotate_globals(
        n_cases=n_cases, n_controls=n_controls, heritability=heritability, saige_version=saige_version, inv_normalized=inv_normalized)
    ht = ht.annotate(total_variants= ht.Number_rare + ht.Number_ultra_rare,
                     interval=gene_ht.key_by('gene_id')[ht.gene_id].interval)
    ht = ht.annotate(CHR= ht.interval.start.contig, POS=ht.interval.start.position)
    # ht = ht.group_by(*list(ht.key)).aggregate(**(hl.agg.take(ht.row_value, 1, -ht.total_variants)[0])) # TODO: add for megabase spanning bug
    ht.describe()
    ht = ht.checkpoint(output_ht_path, overwrite=overwrite, _read_if_exists=not overwrite).drop('n_cases', 'n_controls')
    for max_MAF in [0.01, 0.001, 0.0001, None]:
        name = 'cauchy' if max_MAF is None else str(max_MAF)
        if name != 'cauchy':
            ht.filter(ht.max_MAF == max_MAF).export(output_ht_path.replace('.ht', f'_{name}.txt.bgz'))
        else:
            ht.filter(hl.is_missing(ht.max_MAF)).export(output_ht_path.replace('.ht', f'_{name}.txt.bgz'))

    ht.show()
    print(ht.count())

def load_variant_data(directory: str,
                      output_ht_directory: str,
                      phenoname: str,
                      quantitative_trait:bool,
                      null_glmm_log: str,
                      extension: str = 'single.txt',
                      overwrite: bool = False,
                      variant_type: str = 'genome',
                      num_partitions: int = 1000):
    hl.init(
        master='local[32]',
        tmp_dir=TMP_BUCKET,
        gcs_requester_pays_configuration='aou-neale-gwas',
        worker_memory="highmem",
        worker_cores=8,
        default_reference="GRCh38",
    )
    ## Changelog: removed functions  get_cases_and_controls_from_log(log_format) and get_heritability_from_log(log_file, quantitative_trait: bool = False)
    ## and related lines, since saige log has changed and does not contain those information
    ## Example: gsutil cat gs://aou_analysis/v8/variant_results/result/SAS/phenotype_height/result_height_chr1_000000001.variant.log

    def get_heritability_from_log(log_file, quantitative_trait: bool = False):
        import math
        heritability = -1
        with hfs.open(log_file) as f:
            for line in f:
                if line.startswith('Final'):
                    fields = line.strip().split()
                    if len(fields) == 4:
                        try:
                            tau = float(fields[2])
                            if quantitative_trait:
                                tau1 = float(fields[1])
                                heritability = tau / (tau1 + tau)
                            else:
                                heritability = tau / (tau + math.pi ** 2 / 3)
                            break
                        except:
                            logger.warn(f'Could not load heritability from {line}.')
        return heritability

    def get_saige_version_from_log(null_glmm_log):
        version = 'NA'
        with hfs.open(null_glmm_log) as f:
            for line in f:
                if line.startswith('other attached packages:'):
                    try:
                        line2 = f.readline()
                        packages = line2.strip().split()
                        version = [x for x in packages if 'SAIGE' in x][0]
                    except:
                        logger.warning(f'Could not load version number from {line2} in {null_glmm_log}.')
        return version

    def get_inverse_normalize_status(null_glmm_log):
        status = 'Unknown'
        with hfs.open(null_glmm_log) as f:
            for line in f:
                if line.startswith('$invNormalize'):
                    try:
                        status = f.readline().strip().split()[1]
                    except:
                        logger.warning(f'Could not load inv_norm status from {line} in {null_glmm_log}.')
        return status.capitalize()

    def parse_log_p_value(x):
        log_x = hl.if_else(x.contains('E'), -hl.int(x.split("E")[1]) - hl.log10(hl.float64(x.split("E")[0])), -hl.log10(hl.float64(x)))
        return log_x

    heritability = get_heritability_from_log(null_glmm_log,
                                             quantitative_trait) if null_glmm_log else -1.0
    inv_normalized = get_inverse_normalize_status(null_glmm_log) if null_glmm_log else 'NA'
    saige_version = get_saige_version_from_log(null_glmm_log) if null_glmm_log else 'NA'
    output_ht_path = f'{output_ht_directory}/{variant_type}_variant_results.ht'
    if not hfs.exists(f'{output_ht_path}/_SUCCESS') or overwrite:
        ht = hl.import_table(f'{directory}/*.{extension}', delimiter='\t', impute=True, types = {'p.value':hl.tstr})
        print(f'Loading: {directory}/*.{extension} ...')
        if heritability == -1.0: heritability = hl.null(hl.tfloat)
        if saige_version == 'NA': saige_version = hl.null(hl.tstr)
        if inv_normalized == 'NA': inv_normalized = hl.null(hl.tstr)

        if 'N' in list(ht.row_value):
            ht.describe()
            ht = ht.annotate(
                **{'p.value.NA': hl.missing(hl.tfloat64),
                'Is.SPA': hl.missing(hl.tbool),
                'AF_case': hl.missing(hl.tfloat64),
                'AF_ctrl': hl.missing(hl.tfloat64),
                'N_case': hl.missing(hl.tint32),
                'N_ctrl': hl.missing(hl.tint32),
                }
            )
        elif 'N_event' in list(ht.row_value):
            ht.describe()
            ht = ht.annotate(
                AF_case = ht.AF_event,
                AF_ctrl = ht.AF_censor,
                N_case = ht.N_event,
                N_ctrl = ht.N_censor,
                N = ht.N_event + ht.N_censor
            )
        else:
            ht = ht.annotate(N = ht.N_case + ht.N_ctrl)

        ht = ht.key_by(locus=hl.locus(ht.CHR, ht.POS, reference_genome = 'GRCh38'), alleles=[ht.Allele1, ht.Allele2],
                    phenoname=phenoname).distinct().naive_coalesce(num_partitions)
        ht = ht.transmute(Pvalue=ht['p.value']).annotate_globals(
            heritability=heritability, saige_version=saige_version,
            inv_normalized=inv_normalized)
        ht = ht.drop('Allele1', 'Allele2')
        ht = ht.annotate(Pvalue_log10 = parse_log_p_value(ht.Pvalue))
        ht = ht.annotate(Pvalue = hl.float64(ht.Pvalue))
        ht = ht.select('CHR', 'POS', 'MarkerID', 'AC_Allele2', 'AF_Allele2', 'MissingRate', 'BETA', 'SE', 'Tstat', 'var', 
        'Pvalue', 'Pvalue_log10', 'p.value.NA', 'Is.SPA', 'AF_case', 'AF_ctrl', 'N_case', 'N_ctrl', 'N')
        ht = ht.checkpoint(output_ht_path, overwrite=overwrite, _read_if_exists=not overwrite)
    ht = hl.read_table(output_ht_path)
    ht.describe()
    ht.show()
    overflow_ht = ht.filter(hl.is_nan(ht.Pvalue) | (ht.Pvalue < 1e-1000))
    if overflow_ht.count() > 0:
        warnings.warn(f'Overflowed p-value found in {overflow_ht.count()} variants.')
    print(ht.count())

    rsid_ht = hl.experimental.load_dataset(name='dbSNP', version='154',
                                       reference_genome='GRCh38')
    
    method = 'approx_cdf'
    success = output_ht_path.replace('.ht', f'_{method}_expected_p.ht/_SUCCESS')
    if not hl.hadoop_exists(output_ht_path):
        sys.exit(f'Output HT {output_ht_path} does not exist.')
    elif (hl.hadoop_exists(output_ht_path) and not hl.hadoop_exists(success)) or overwrite:
        ht = annotate_expected_pvalue(ht=ht, method=method, p_field='Pvalue')
        ht = ht.annotate(rsID=rsid_ht[ht.locus, ht.alleles].rsid).key_by()
        ht.describe()
        ht = ht.checkpoint(output_ht_path.replace('.ht', f'_{method}_expected_p.ht'), overwrite=True)
        ht = ht.filter(hl.case()
                        .when(ht.Pvalue > 0.1, hl.rand_bool(0.001))
                        .when(ht.Pvalue > 0.01, hl.rand_bool(0.01))
                        .when(ht.Pvalue > 0.001, hl.rand_bool(0.1))
                        .default(True))
        ht.export(output_ht_path.replace('.ht', f'_{method}_expected_p.txt.bgz'))
    else:
        ht = hl.read_table(output_ht_path.replace('.ht', f'_{method}_expected_p.ht'))
    return ht

def run_saige(
    p: Batch,
    phenoname: str,
    ancestry:str,
    output_root: str,
    model_file: str,
    variance_ratio_file: str,
    sparse_grm_file: str,
    bgen_file: ResourceGroup,
    samples_file: ResourceGroup,
    docker_image: str,
    group_file: str,
    groups: str,
    trait_type: str,
    chrom: str,
    max_maf_for_group: str,
    variant_type:str,
    min_mac: int = 1,
    min_maf: float = 0,
    memory: str = "",
    storage: str = "10Gi",
    add_suffix: str = "",
):
    """
    Change log:
    use_bgen and log_pvalue are defaults now
    - Removed --IsOutputlogPforSingle, --sparseSigmaFile, --IsOutputPvalueNAinGroupTestforBinary, IsOutputBETASEinBurdenTest, IsOutputAFinCaseCtrl
    - Added --sparseGRMFile, --sparseGRMSampleIDFile, --is_fastTest, --is_output_markerList_in_groupTest
    - Changed --maxMAFforGroupTest to --maxMAF_in_groupTest, --IsSingleVarinGroupTest to --is_single_in_groupTest
    :param p:
    :param output_root:
    :param model_file:
    :param variance_ratio_file:
    :param sparse_grm_file:
    :param bgen_file:
    :param samples_file:
    :param docker_image:
    :param group_file:
    :param use_bgen:
    :param trait_type:
    :param chrom:
    :param min_mac:
    :param min_maf:
    :param max_maf_for_group:
    :param memory:
    :param storage:
    :param add_suffix:
    :return:
    """
    MKL_OFF = "export MKL_NUM_THREADS=1; export MKL_DYNAMIC=false; export OMP_NUM_THREADS=1; export OMP_DYNAMIC=false; "
    analysis_type = "gene" if group_file is not None else "variant"
    run_saige_task: Job = (
        p.new_job(name=f"run_saige{'' if analysis_type=='variant' else '_gene'}_{phenoname}_{ancestry}", attributes={"analysis_type": analysis_type, "ancestry": ancestry})
        .cpu(1)
        .storage(storage)
        .image(docker_image)
        .memory(memory)
    )  # Step 2 is single-threaded only

    if analysis_type == "gene":
        run_saige_task.declare_resource_group(
            result={
                f"{add_suffix}gene.txt": "{root}",
                f"result.singleAssoc.txt": "{root}.singleAssoc.txt",
            }
        )
    else:
        run_saige_task.declare_resource_group(result={"single_variant.txt": "{root}"})

    command = (
        f"set -o pipefail; {MKL_OFF} Rscript /usr/local/bin/step2_SPAtests.R "
        f"--bgenFile={bgen_file.bgen} " 
        f'--bgenFileIndex={bgen_file["bgen.bgi"]} '
        f"--chrom={chrom} "  
        f"--minMAF={min_maf} " 
        f"--minMAC={min_mac} " 
        f"--sampleFile={samples_file} " 
        f"--GMMATmodelFile={model_file} " 
        f"--varianceRatioFile={variance_ratio_file} " 
        f"--AlleleOrder=ref-first " 
        f"--SAIGEOutputFile={run_saige_task.result} "
    )
    if variant_type == 'exome':
        command += (f"--cateVarRatioMinMACVecExclude=0.5,1.5,2.5,3.5,4.5,5.5,10.5,15.5,20.5 " 
                    f"--cateVarRatioMaxMACVecInclude=1.5,2.5,3.5,4.5,5.5,10.5,15.5,20.5 " )
    if analysis_type == "gene":
        if trait_type == "binary":
            command += (
                f"--is_output_moreDetails=TRUE "
            )
        command += (
            f"--groupFile={group_file} " 
            f"--annotation_in_groupTest={groups} "
            f"--is_output_markerList_in_groupTest=TRUE "
            f"--maxMAF_in_groupTest={max_maf_for_group} "
            f"--is_single_in_groupTest=TRUE "
            f"--sparseGRMFile={sparse_grm_file} "
            f"--sparseGRMSampleIDFile={sparse_grm_file}.sampleIDs.txt "
        )
    command += f"--LOCO=FALSE 2>&1 | tee {run_saige_task.stdout}; "
    if analysis_type == "gene":
        command += (
            f"input_length=$(wc -l {group_file} | awk '{{print $1}}'); "
            f"output_length=$(wc -l {run_saige_task.result[f'{add_suffix}gene.txt']} | awk '{{print $1}}'); "
            f"echo 'Got input:' $input_length 'output:' $output_length | tee -a {run_saige_task.stdout}; "
            f"if [[ $input_length > 0 ]]; then echo 'got input' | tee -a {run_saige_task.stdout}; "
            f"if [[ $output_length == 1 ]]; then echo 'but not enough output' | tee -a {run_saige_task.stdout}; "
            f"rm -f {run_saige_task.result[f'{add_suffix}gene.txt']} exit 1; fi; fi"
        )
    run_saige_task.command(command)
    # Check paths:
    if analysis_type == "gene":
        run_saige_task.command(f'dir={run_saige_task.result[f"{add_suffix}gene.txt"]}')
        run_saige_task.command(f'parentdir="$(dirname "$dir")"')
        run_saige_task.command(f'ls -lh $parentdir')

    p.write_output(run_saige_task.result, output_root)
    p.write_output(run_saige_task.stdout, f"{output_root}.{analysis_type}.log")
    return run_saige_task

def fit_null_glmm(
    p: Batch,
    ancestry:str,
    output_root: str,
    phenoname: str,
    pheno_file: Resource,
    trait_type: str,
    covariates: str,
    plink_file_root: str,
    docker_image: str,
    variant_type:str,
    sparse_grm: Resource = None,
    sparse_grm_extension: str = None,
    inv_normalize: bool = True,
    skip_model_fitting: bool = False,
    n_threads: int = 16,
    storage: str = "10Gi",
    memory: str = "60G",
    non_pre_emptible: bool = False,
):
    """
    Change log:
    - Removed --minCovariateCount, sparseSigma output
    - Added --useSparseGRMforVarRatio, --useSparseGRMtoFitNULL
    :param p:
    :param output_root:
    :param phenoname:
    :param pheno_file:
    :param trait_type:
    :param covariates:
    :param plink_file_root:
    :param docker_image:
    :param sparse_grm:
    :param sparse_grm_extension:
    :param inv_normalize:
    :param skip_model_fitting:
    :param n_threads:
    :param storage:
    :param memory:
    :param non_pre_emptible:
    :return:
    """

    pheno_col = "value" if trait_type != 'survival' else "secondEvent"
    user_id_col = "userId" if phenoname.startswith('random') else "person_id"
    if trait_type == 'survival': user_id_col = "s"
    in_bfile = p.read_input_group(
        **{ext: f"{plink_file_root}.{ext}" for ext in ("bed", "bim", "fam")}
    )
    fit_null_task = (
        p.new_job(
            name=f"fit_null_model_{phenoname}_{ancestry}",
            attributes={ "ancestry": ancestry},
        )
        .storage(storage)
        .image(docker_image)
    )
    if non_pre_emptible:
        fit_null_task._cpu = None
        fit_null_task._memory = None
        fit_null_task._machine_type = "n1-highmem-8"
        fit_null_task._preemptible = False
    else:
        fit_null_task = fit_null_task.cpu(n_threads).memory(memory)
    output_files = {
        ext: f'{{root}}{ext if ext.startswith("_") else "." + ext}'
        for ext in (
            "rda",
            f"varianceRatio.txt",
        )
    }

    fit_null_task.declare_resource_group(null_glmm=output_files)
    bim_fix_command = f"perl -pi -e s/^chr// {in_bfile.bim}"

    command = (
        f"set -o pipefail; Rscript /usr/local/bin/step1_fitNULLGLMM.R "
        f"--sparseGRMFile={sparse_grm[sparse_grm_extension]} "
        f'--sparseGRMSampleIDFile={sparse_grm[f"{sparse_grm_extension}.sampleIDs.txt"]} '
        f"--plinkFile={in_bfile} "
        f"--useSparseGRMtoFitNULL=TRUE "
        f"--phenoFile={pheno_file} "
        f"--covarColList={covariates} "
        f"--phenoCol={pheno_col} "
        f"--sampleIDColinphenoFile={user_id_col} "
        f"--traitType={trait_type} "
        f"--isCateVarianceRatio=TRUE "
        f"--useSparseGRMforVarRatio=TRUE "
        f"--outputPrefix={fit_null_task.null_glmm} "
        f"--outputPrefix_varRatio={fit_null_task.null_glmm} "
        f"--skipModelFitting={str(skip_model_fitting).upper()} "
    )
    if inv_normalize:
        command += "--invNormalize=TRUE "
    if variant_type == 'exome':
        command += (f"--cateVarRatioMinMACVecExclude=0.5,1.5,2.5,3.5,4.5,5.5,10.5,15.5,20.5 " 
                    f"--cateVarRatioMaxMACVecInclude=1.5,2.5,3.5,4.5,5.5,10.5,15.5,20.5 " )
    if trait_type == 'survival':
        command += f"--eventTimeCol=TTE "
        command += f"--eventTimeBinSize=1 "
        command += f"--traceCVcutoff=0.0025 "
        command += f"--ratioCVcutoff=0.001 "
    command += f"--nThreads={n_threads} --LOCO=FALSE 2>&1 | tee {fit_null_task.stdout}"
    command = "; ".join([bim_fix_command, command])
    fit_null_task.command(command)
    # Check paths:
    # fit_null_task.command(f'dir={fit_null_task.null_glmm[f"rda"]}')
    # fit_null_task.command(f'parentdir="$(dirname "$dir")"')
    # fit_null_task.command(f'ls -lh $parentdir')

    p.write_output(fit_null_task.null_glmm, output_root)
    p.write_output(fit_null_task.stdout, f"{output_root}.log")
    # Runtimes: 8 threads: ~5 minutes of 100% CPU (~3G RAM), followed by ~9 minutes of 800% (~6G RAM)
    return fit_null_task

def gt_to_gp(mt, location: str = "GP"):
    return mt.annotate_entries(
        **{
            location: hl.or_missing(
                hl.is_defined(mt.GT),
                hl.map(
                    lambda i: hl.if_else(
                        mt.GT.unphased_diploid_gt_index() == i, 1.0, 0.0
                    ),
                    hl.range(0, hl.triangle(hl.len(mt.alleles))),
                ),
            )
        }
    )


def impute_missing_gp(mt, location: str = "GP", mean_impute: bool = True):
    mt = mt.annotate_entries(_gp=mt[location])
    if mean_impute:
        mt = mt.annotate_rows(
            _mean_gp=hl.agg.array_agg(lambda x: hl.agg.mean(x), mt._gp)
        )
        gp_expr = mt._mean_gp
    else:
        gp_expr = [1.0, 0.0, 0.0]
    return mt.annotate_entries(**{location: hl.or_else(mt._gp, gp_expr)}).drop("_gp")


def export_bgen_from_mt(
    ancestry,
    analysis_type,
    interval,
    output_dir,
    log_file:str,
    mean_impute_missing: bool=True, # Whether to mean impute missing genotypes
    variant_ac_filter: int = 0,
    variant_callrate_filter: float=0,
):
    hl.init(
        master='local[24]',
        tmp_dir=TMP_BUCKET,
        driver_memory="highmem",
        driver_cores=8,
        worker_memory="highmem",
        default_reference="GRCh38",
        log=log_file
    )
    outname = f"{analysis_type}_{interval.start.contig}_{str(interval.start.position).zfill(9)}_{interval.end.position}"
    mt_type = 'ACAF' if analysis_type == 'variant' else 'Exome'
    print('Loading sample IDs')
    sample_ids = hl.read_table(f'{DATA_PATH}/utils/aou_v8_hard_filter.ht')
    print('Loading MT')
    mt = get_filtered_mt(mt_type=mt_type, sample_ids=sample_ids, filter_variants=True, filter_samples=True, adj_filter= True, ancestry=ancestry, prune_samples=True)
    print('MT loaded')
    # Filter to interval
    mt = hl.filter_intervals(mt, [interval])

    mt = mt.select_entries("GT")
    mt = mt.filter_rows(
        hl.agg.count_where(mt.GT.is_non_ref()) > 0
    )  # Filter to non-reference sites
    mt = mt.annotate_rows(
        rsid=mt.locus.contig + ":" + hl.str(mt.locus.position) + "_" + mt.alleles[0] + "/" + mt.alleles[1]
    )  # Annotate rsid

    call_stats_ht = hl.read_table(f'{DATA_PATH}/utils/call_stats/{mt_type.lower() if mt_type == "Exome" else mt_type}_pruned/{ancestry.upper()}_{mt_type.lower() if mt_type == "Exome" else mt_type}_call_stats.ht')
    call_stats_ht = call_stats_ht.filter(
        (call_stats_ht.call_stats.AC[1] > variant_ac_filter) &
        (call_stats_ht.call_stats.AN/(2*N_SAMPLES_PRUNED[ancestry]) > variant_callrate_filter)
    )
    mt = mt.filter_rows(hl.is_defined(call_stats_ht[mt.row_key]))

    mt = mt.annotate_entries(
        GT=hl.if_else(mt.GT.is_haploid(), hl.call(mt.GT[0], mt.GT[0]), mt.GT)
    )
    mt = gt_to_gp(mt)
    mt = impute_missing_gp(mt, mean_impute=mean_impute_missing)
    hl.export_bgen(mt, f"{output_dir}/{outname}", gp=mt.GP, varid=mt.rsid)
    print(hl.utils.range_table(10)._force_count())


def export_gene_group_file(interval, ancestry, output_dir):
    gene_ht_path = f'{DATA_PATH}/utils/gene_map/aou_{ancestry.upper()}_gene_map_processed_{TRANCHE}.ht'
    gene_ht = hl.read_table(gene_ht_path)
    gene_ht = hl.filter_intervals(gene_ht, [interval])
    var_ht = gene_ht.select(gene=hl.if_else(hl.is_missing(gene_ht.gene_id), '_', gene_ht.gene_id) + '_' + hl.if_else(hl.is_missing(gene_ht.gene_symbol), '_', gene_ht.gene_symbol),
                            tag='var',
                            info=hl.delimit(gene_ht.variants, " ")
                            ).add_index().key_by().drop("start")
    anno_ht = gene_ht.select(gene=hl.if_else(hl.is_missing(gene_ht.gene_id), '_', gene_ht.gene_id) + '_' + hl.if_else(hl.is_missing(gene_ht.gene_symbol), '_', gene_ht.gene_symbol),
                             tag='anno',
                             info=gene_ht.annotation[:-1]
                             ).add_index().key_by().drop("start")
    group_ht = var_ht.union(anno_ht)
    group_ht = group_ht.group_by('gene', 'tag').aggregate(info = hl.agg.collect(group_ht.info))
    group_ht = group_ht.annotate(info  = hl.str(' ').join(group_ht.info))
    group_ht = group_ht.order_by(group_ht.gene, hl.desc(group_ht.tag))
    group_ht.show()
    # group_ht = group_ht.order_by('idx').drop('idx')
    group_ht.export(output_dir, header=False, delimiter=' ')

def index_bgen(b: hb.batch.Batch, ancestry:str, analysis_type:str, bgen: str, depend_job=None):
    file = b.read_input(bgen)
    name = bgen.split('/')[-1]
    j = b.new_job(name=f'index_{name}_{ancestry}_{analysis_type}')
    if depend_job is not None:
        j.depends_on(depend_job)
    j.image('befh/bgen:latest')
    j.command(f'bgenix -index -g {file}')
    j.command(f'ls {file}.bgi')
    j.command(f'mv {file}.bgi {j.temp}')
    b.write_output(j.temp, f'{bgen}.bgi')
    return j

def export_pheno_gate(category: str, pheno_path:str, phenoname: str, fields: str):
    ht = hl.read_table(f'{DATA_PATH}/phenotype/{category}/{phenoname}.ht')
    fields = fields.split(',')
    out_ht = ht.select(**{field: ht[field] for field in fields})
    out_ht.describe()
    print(out_ht.count())
    out_ht.export(pheno_path)

def export_pheno(category: str, pheno_path:str, fields: str, ancestry: str, phenoname: str, binary_traits: list, proportion_single_sex:float):
    def irnt(he: hl.expr.Expression, output_loc: str = 'irnt'):
        ht = he._indices.source
        n_rows = ht.aggregate(hl.agg.count_where(hl.is_defined(he)))
        print(n_rows)
        ht = ht.order_by(he).add_index()
        ht = ht.annotate(**{output_loc: hl.qnorm((ht.idx + 0.5) / n_rows)})
        ht = ht.annotate(**{output_loc: hl.or_missing(~hl.is_nan(ht[output_loc]), ht[output_loc])})
        return ht

    def filter_to_single_sex_by_proportion(ht, proportion_single_sex, tag='both'):
        # If proportion of male or female cases is less than `proportion_single_sex`
        # then filter to females and males respectively
        n_cases_female = ht.aggregate(
            hl.agg.filter(
                ht.sex == 0,
                hl.agg.count_where(ht.value == 1),
            )
        )
        n_cases_male = ht.aggregate(
            hl.agg.filter(
                ht.sex == 1,
                hl.agg.count_where(ht.value == 1),
            )
        )
        prop_female = n_cases_female / (n_cases_male + n_cases_female)
        print(f"Female proportion: {prop_female}")
        if prop_female <= proportion_single_sex:
            print(
                f"Female case proportion {prop_female} less than {proportion_single_sex}. Filtering to males..."
            )
            ht = ht.filter(ht.sex == 1)
            tag = 'male_only'
        elif prop_female >= 1 - proportion_single_sex:
            print(
                f"Female case proportion {prop_female} greater than {1 - proportion_single_sex}. Filtering to females..."
            )
            ht = ht.filter(ht.sex == 0)
            tag = 'female_only'
        log_row = {
            'phenoname': phenoname,
            'single_sex': tag,
            'female_proportion': prop_female,
            'single_sex_threshold': proportion_single_sex
        }
        log_ht = hl.Table.parallelize([log_row])
        log_ht.export(pheno_path.replace('.tsv', '.log'))
        return ht

    def load_pheno_ht(category:str, ancestry:str, phenoname:str):
        if category in PHENO_CATEGORIES_MT:
            mt = hl.read_matrix_table(f'{DATA_PATH}/phenotype/{category}_annotated.mt')
            ancestry_mt = mt.filter_rows(mt.ancestry == ancestry)
            pheno_mt = ancestry_mt.filter_cols(ancestry_mt.phenoname == phenoname)
            ht = pheno_mt.entries()
            print(f'-------{category} {ancestry} {phenoname} {ht.count()}')
        elif category == 'random_pheno':
            pheno_ht = hl.read_table(f'{DATA_PATH}/phenotype/random_pheno_annotated.ht')
            ht = pheno_ht.filter(pheno_ht.ancestry == ancestry)
            ht = ht.annotate(value=ht[phenoname])
        elif category in PHENO_CATEGORIES_HT:
            if category == 'lab_measurement':
                path_name = phenoname.split('_')[0]
                stat = phenoname.split('_')[1]
            else:
                path_name = phenoname
            ht = hl.read_table(f"{DATA_PATH}/phenotype/{category}/{'lab_measurement_' if category == 'lab_measurement' else ''}{path_name}.ht")
            ht = ht.filter(ht.ancestry == ancestry)
            ht = ht.annotate(value=ht[f'{stat}_value'])
        ht = ht.checkpoint(f'gs://aou_tmp/export_pheno/{ancestry}_{phenoname}.ht', overwrite = True)
        print(f"Exporting {phenoname} for {ancestry.upper()} (n_samples: {ht.count()})...")
        return ht
    tag = 'both'
    if phenoname.endswith('_male'):
        tag = 'male_only'
        path_name = phenoname.replace('_male', '')
    elif phenoname.endswith('_female'):
        tag = 'female_only'
        path_name = phenoname.replace('_female', '')
    else:
        path_name = phenoname
    ht = load_pheno_ht(category=category, ancestry=ancestry, phenoname=path_name)

    if phenoname in binary_traits:
        ht = ht.annotate(value=hl.int(ht.value))
        if tag == 'both':
            ht = filter_to_single_sex_by_proportion(ht=ht, proportion_single_sex=proportion_single_sex)
    else:
        print('Computing irnt values...')
        ht = ht.annotate(value = hl.float64(ht['value']))
        ht = irnt(ht['value'], 'irnt_value')
        ht = ht.annotate(value=ht.irnt_value)
    user_id_col = "userId" if phenoname.startswith('random') else "person_id"
    if category in ['onset', 'progression']: user_id_col = "s"
    ht = ht.key_by(user_id_col)
    fields = fields.split(",") + ["value"]
    out_ht = ht.select(**{field: ht[field] for field in fields})
    if tag == 'male_only':
        out_ht = out_ht.filter(out_ht.sex == 1)
    elif tag == 'female_only':
        out_ht = out_ht.filter(out_ht.sex == 0)
    out_ht.export(pheno_path)

def write_pickle_dict(output: str, dict: dict):
    with hfs.open(output, "wb") as f:
        pickle.dump(dict, f)
    f.close()

def read_pickle_dict(dict_path: str):
    dict = pd.read_pickle(dict_path)
    return dict



def main(args):
    hl.init(
        tmp_dir=TMP_BUCKET,
        driver_memory="highmem",
        driver_cores=8,
        worker_memory="highmem",
        worker_cores=1,
        default_reference="GRCh38",
        app_name=f'aou_SAIGE_{args.ancestries.replace(",","_")}' if not args.export_phenos else "aou_export_phenos"
    )

    num_pcs = 20
    start_time = time.time()
    n_threads = 8
    proportion_single_sex = 0.1
    ancestries = args.ancestries.split(",")
    PC_fields = [f"PC{x}" for x in range(1, num_pcs + 1)]
    basic_covars = [
        "sex",
        "age",
        "age2",
        "age_sex",
        "age2_sex",
    ]  
    basic_covariates = basic_covars + PC_fields
    basic_covariates = ",".join(basic_covariates)

    if args.export_phenos or (not args.skip_saige) or (not args.skip_any_null_models) or (not args.skip_load_hail_results) or (not args.skip_bgen):
        all_phenos_by_group = read_pickle_dict(f'{DATA_PATH}/phenotype/summary/pheno_dict_raw.dict')
        print(f'----------Number of phenotypes per category (RAW): --------------')
        print([(category, len(all_phenos_by_group[category])) for category in all_phenos_by_group.keys()])

        quantitative_traits = []
        binary_traits = []
        gate_traits = []
        if 'random_pheno' in all_phenos_by_group.keys():
            quantitative_traits = [pheno for pheno in all_phenos_by_group['random_pheno'] if 'continuous' in pheno]
            binary_traits = [pheno for pheno in all_phenos_by_group['random_pheno'] if 'continuous' not in pheno]
        for category in QUANTITATIVE_CATEGORIES:
            quantitative_traits = quantitative_traits + list(all_phenos_by_group[category])
        for category in BINARY_CATEGORIES:
            binary_traits = binary_traits + list(all_phenos_by_group[category])
        for category in GATE_CATEGORIES:
            gate_traits = gate_traits + list(all_phenos_by_group[category])
        print(f'Number of quantitative traits: {len(quantitative_traits)}\nNumber of binary traits: {len(binary_traits)}\nNumber of gate traits: {len(gate_traits)}')

        phenos_to_run_by_ancestry = {}
        phenos_to_run_by_ancestry_by_group = {}
        if args.pilot or (args.phenos is not None):
            phenos_to_run = list(PILOT_PHENOTYPES) if args.pilot else args.phenos.split(",")
            print(phenos_to_run)
            category_lst = []
            for phenoname in phenos_to_run:
                for category in PHENO_CATEGORIES:
                    if phenoname in all_phenos_by_group[category]:
                        category_lst.append(category)
            print(category_lst)
            category_dict = dict(zip(phenos_to_run, category_lst))
            print(category_dict)

            for ancestry in ancestries:
                phenos_to_run_by_ancestry_by_group[ancestry] = {}
                for category in set(category_lst):
                    if category in GATE_CATEGORIES:
                        phenos_to_run_by_ancestry_by_group[ancestry][category] = [k for k, v in category_dict.items() if v == category and ancestry in k]
                    else:
                        phenos_to_run_by_ancestry_by_group[ancestry][category] = [k for k, v in category_dict.items() if v == category]
                phenos_to_run_by_ancestry[ancestry] = phenos_to_run
            if len(phenos_to_run) < 20:
                print(phenos_to_run_by_ancestry)
        else: 
            phenos_to_run_by_ancestry = read_pickle_dict(f'{DATA_PATH}/phenotype/summary/pheno_by_ancestry_dict.dict')
            phenos_to_run_by_ancestry_by_group = read_pickle_dict(f'{DATA_PATH}/phenotype/summary/pheno_by_ancestry_by_group_dict.dict')
            print(f'----------Number of phenotypes per category per pop (filtered to n_cases >= 200): --------------')
            print([f"{ancestry.upper()}-{category}: {len(phenos_to_run_by_ancestry_by_group[ancestry][category])}" for ancestry in ancestries for category in list(all_phenos_by_group.keys())])
            print([f"{ancestry.upper()}: {len(phenos_to_run_by_ancestry[ancestry])}" for ancestry in ancestries])
    else:
        phenos_to_run_by_ancestry = {}
        for ancestry in ancestries:
            phenos_to_run_by_ancestry[ancestry] = None
        phenos_to_run_by_ancestry_by_group = None

    backend = hb.ServiceBackend(
        billing_project="all-by-aou", remote_tmpdir=TMP_BUCKET
    )

    if args.export_phenos:
        b = hb.Batch(
            name=f"aou_export_pheno_{args.ancestries}_{'test' if args.test or args.phenos is not None else 'all'}",
            backend=backend,
            default_storage="500Mi",
            default_cpu=n_threads,
        )
        
        for ancestry in ancestries:
            print(f'Processing {ancestry.upper()}...')
            for category in PHENO_CATEGORIES:
                print(f'Processing {category}...')
                fields = basic_covariates
                if category not in list(phenos_to_run_by_ancestry_by_group[ancestry].keys()): continue
                phenos_to_export = phenos_to_run_by_ancestry_by_group[ancestry][category]
                print(f"------------{ancestry.upper()} phenotype info: {len(phenos_to_export)} {category} phenotypes------------")
                for phenoname in tqdm(phenos_to_export):
                    if (not ancestry in phenoname) and (category in GATE_CATEGORIES):
                        continue
                    if category in GATE_CATEGORIES:
                        if category == 'onset':
                            fields = PC_fields + ['birth_year', 'TTE', 'secondEvent']
                        else:
                            fields = PC_fields + ['birth_year', 'TTE', 'secondEvent', 'age_at_first_event']
                        if phenoname.endswith("_ALL"):
                            fields = fields + ["sex"]
                        fields = ",".join(fields)
                    pheno_path = f'{ANALYSIS_BUCKET}/pheno_file/{ancestry.upper()}/phenotype_{phenoname}.tsv'
                    print(pheno_path)
                    if (not hfs.exists(pheno_path)
                            or args.overwrite_pheno_data
                    ):
                        j = b.new_python_job(f'export_phenotype_{phenoname}_for_{ancestry}',
                                             attributes={"ancestry": ancestry, "phenotype": copy.deepcopy(phenoname), "category": category})
                        j.image("hailgenetics/hail:0.2.133-py3.11")
                        if category not in GATE_CATEGORIES:
                            j.memory('highmem')
                            j.cpu(8) 
                            j.env('PYSPARK_SUBMIT_ARGS', '--driver-memory 24g --executor-memory 24g pyspark-shell')    
                            j.call(export_pheno,
                            category=category,
                            fields=fields,
                            pheno_path=pheno_path,
                            phenoname=phenoname,
                            ancestry=ancestry.upper(),
                            binary_traits=binary_traits,
                            proportion_single_sex=proportion_single_sex,
                            )
                        else:
                            j.call(export_pheno_gate,
                            category=category,
                            fields=fields,
                            pheno_path=pheno_path,
                            phenoname=phenoname,
                            )

                if args.test:
                    break
        b.run()

    if args.run_pipeline:
        analysis_type = "variant" if args.single_variant_only else "gene"
        variant_type = 'genome' if analysis_type == 'variant' else 'exome'
        print(f'Analysis type: {analysis_type}')
        print(f'Variant_type: {variant_type}')
        print(f"Docker image: {SAIGE_DOCKER_IMAGE}...")
        RESULT_ROOT = f'{ANALYSIS_BUCKET}/{analysis_type}_results'

        for ancestry in ancestries:
            if not args.skip_bgen or not args.skip_saige:
                N_GENE_PER_GROUP = 50 if ancestry=='all' else 100
                size = CHUNK_SIZE[ancestry] if analysis_type == 'variant' else N_GENE_PER_GROUP
                interval_ht_path = f'{DATA_PATH}/utils/intervals/aou_{analysis_type}_interval_size_{size}.ht'
                if not hfs.exists(interval_ht_path.replace('ht', 'pickle')):
                    interval_ht = hl.read_table(interval_ht_path)
                    interval_ht = interval_ht.filter(interval_ht.interval.start.contig != "chrM")
                    intervals = interval_ht.aggregate(hl.agg.collect(interval_ht.interval))
                    write_pickle_dict(
                        interval_ht_path.replace('ht', 'pickle'),
                        intervals)
                intervals = read_pickle_dict(
                    interval_ht_path.replace('ht', 'pickle'))
                print(intervals[0:3])
                print(f'---------Number of intervals [{ancestry.upper()}]: {len(intervals)}---------')   
            if args.category is not None:
                category = args.category
                if args.category == 'small':
                    phenos_to_run = set(phenos_to_run_by_ancestry_by_group[ancestry]["lab_measurement"] + phenos_to_run_by_ancestry_by_group[ancestry]['random_pheno'] + \
                                    phenos_to_run_by_ancestry_by_group[ancestry]['pfhh_survey_table']+ phenos_to_run_by_ancestry_by_group[ancestry]["physical_measurement"])
                elif args.category == 'quantitative':
                    phenos_to_run = set(
                        phenos_to_run_by_ancestry_by_group[ancestry]["lab_measurement"] + phenos_to_run_by_ancestry_by_group[ancestry]["physical_measurement"])
                elif args.category == 'binary':
                    phenos_to_run = set(
                        phenos_to_run_by_ancestry_by_group[ancestry]["r_drug"] + phenos_to_run_by_ancestry_by_group[ancestry]["pfhh_survey_table"])
                elif args.category == 'gate':
                    phenos_to_run = set(
                        phenos_to_run_by_ancestry_by_group[ancestry]["onset"] + phenos_to_run_by_ancestry_by_group[ancestry]["progression"])
                else:
                    category = args.category
                    phenos_to_run = phenos_to_run_by_ancestry_by_group[ancestry][category]
                print(list(phenos_to_run)[0:5])
            else:
                category = 'all'
                phenos_to_run_by_ancestry['all'] = None
                phenos_to_run = phenos_to_run_by_ancestry[ancestry.lower()]

            b = hb.Batch(
                name=f"saige_{analysis_type}_aou_{ancestry}",
                backend=backend,
                default_image=SAIGE_DOCKER_IMAGE,
                default_storage="500Mi",
                default_cpu=n_threads,
            )

            # Obtain Sparse GRMs
            relatedness_cutoff = "0.125"
            num_markers = 2000
            n_threads = 8
            sparse_grm_root = f"{DATA_PATH}/utils/grm/aou_{ancestry}"
            sparse_grm_extension = f"_relatednessCutoff_{relatedness_cutoff}_{num_markers}_randomMarkersUsed.sparseGRM.mtx"
            sparse_grm = b.read_input_group(
                **{ext: f"{sparse_grm_root}.{ext}"
                   for ext in (sparse_grm_extension, f"{sparse_grm_extension}.sampleIDs.txt",)}
            )

            overwrite_null_models = args.overwrite_null_models
            null_model_dir = f'{RESULT_ROOT}/null_glmm/{ancestry.upper()}'
            null_models_already_created = {}
            null_models = {}
            pheno_exports = {}
            if phenos_to_run is not None and not args.skip_any_null_models:
                print(f"------------{ancestry.upper()} {analysis_type} analysis null models: {len(phenos_to_run)} phenotypes------------")
                print(f'Null model directory: {null_model_dir}')
                if (not overwrite_null_models) and hfs.exists(null_model_dir):
                    null_models_already_created = {
                        x["path"] for x in hl.hadoop_ls(null_model_dir)
                    }
                print(f'Found {int(len(null_models_already_created) / 3)} Null models in directory...')
                for i in tqdm(range(len(phenos_to_run))):
                    phenoname = list(phenos_to_run)[i]
                    null_glmm_root = f"{null_model_dir}/phenotype_{phenoname}"
                    model_file_path = f"{null_glmm_root}.rda"
                    variance_ratio_file_path =  f"{null_glmm_root}.varianceRatio.txt"
                    current_trait_type = None
                    if phenoname in gate_traits:
                        if not ancestry in phenoname:
                            continue
                        current_trait_type = 'survival'
                        covariates = PC_fields + ['birth_year']
                        if not phenoname.startswith('Birth'):
                            covariates = covariates + ['age_at_first_event']
                        if phenoname.endswith('_ALL'):
                            covariates = covariates + ['sex']
                        covariates = ",".join(covariates)
                    elif phenoname in quantitative_traits:
                        current_trait_type = 'quantitative'
                        covariates = basic_covariates
                    elif phenoname in binary_traits:
                        current_trait_type = 'binary'
                        covariates = basic_covariates
                    else:
                        raise ValueError(f"Unknown trait type for {phenoname}")

                    pheno_file = b.read_input(f'{ANALYSIS_BUCKET}/pheno_file/{ancestry.upper()}/phenotype_{phenoname}.tsv')
                    pheno_exports[phenoname] = pheno_file

                    if (
                        not overwrite_null_models
                        and model_file_path in null_models_already_created
                        and variance_ratio_file_path in null_models_already_created
                    ):
                        model_file = b.read_input(model_file_path)
                        variance_ratio_file = b.read_input(variance_ratio_file_path)
                    else:
                        if args.skip_any_null_models:
                            break
                        print(f'Running null model for {ancestry.upper()} {phenoname} : {current_trait_type}')
                        fit_null_task = fit_null_glmm(
                            b,
                            ancestry=ancestry,
                            output_root=null_glmm_root,
                            phenoname=phenoname,
                            pheno_file=pheno_exports[phenoname],
                            trait_type=current_trait_type,
                            covariates=covariates,
                            plink_file_root=f'{DATA_PATH}/utils/grm/{ancestry.upper()}_grm_plink',
                            docker_image=SAIGE_DOCKER_IMAGE,
                            variant_type=variant_type,
                            sparse_grm=sparse_grm,
                            sparse_grm_extension=sparse_grm_extension,
                            inv_normalize=True,
                            skip_model_fitting=False,
                            n_threads=n_threads,
                            storage="100Gi",
                            memory="60Gi",
                            non_pre_emptible=False
                        )
                        fit_null_task.attributes.update(
                            {"phenotype": copy.deepcopy(phenoname),
                            "analysis_type": copy.deepcopy(analysis_type),
                            "trait_type": copy.deepcopy(current_trait_type)}
                        )
                        model_file = fit_null_task.null_glmm.rda
                        variance_ratio_file = fit_null_task.null_glmm[f"varianceRatio.txt"]

                    null_models[phenoname] = (model_file, variance_ratio_file)
                else:
                    print('No phenotype loaded...')

            if not args.skip_bgen or not args.skip_saige:
                print(f"------------{ancestry.upper()} {analysis_type} analysis bgen files------------")
                bgen_dir = f'{RESULT_ROOT}/bgen/{ancestry.upper()}'
                print(f'bgen directory: {bgen_dir}')
                overwrite_bgens = args.overwrite_bgens
                bgens_already_created = {}
                if not overwrite_bgens and hfs.exists(bgen_dir):
                    bgens_already_created = {x["path"] for x in hl.hadoop_ls(bgen_dir) if x["path"].endswith(".bgen")}
                print(f'Found {len(bgens_already_created)} Bgens in directory...')

                bgens = {}
                if args.test:
                    intervals = [intervals[0]]

                if args.hail_image is None:
                    image = hb.build_python_image(
                        "us-central1-docker.pkg.dev/aou-neale-gwas/hail-batch/hail_gnomad_python_3.9",
                        requirements=["hail", "gnomad"],
                        python_version="3.9",
                    )
                else:
                    image = args.hail_image

            if not args.skip_bgen:
                for interval in intervals:
                    bgen_root = f"{bgen_dir}/{analysis_type}_{interval.start.contig}_{str(interval.start.position).zfill(9)}_{interval.end.position}"
                    if (f"{bgen_root}.bgen" not in bgens_already_created) or args.overwrite_bgens:
                        bgen_task = b.new_python_job(
                            name=f"{analysis_type}_analysis_export_{str(interval)}_bgen_{ancestry}"
                        )

                        bgen_task.image(image)
                        bgen_task.memory('highmem')
                        bgen_task.always_copy_output()
                        bgen_task.env('PYSPARK_SUBMIT_ARGS', '--driver-memory 24g --executor-memory 24g pyspark-shell')
                        bgen_task.call(
                            export_bgen_from_mt,
                            ancestry=ancestry,
                            analysis_type=analysis_type,
                            interval=interval,
                            output_dir=bgen_dir,
                            log_file=bgen_task.log_file,
                            mean_impute_missing=True,
                            variant_ac_filter= args.variant_ac_filter,
                            variant_callrate_filter=args.callrate_filter,
                        )
                        b.write_output(bgen_task.log_file, f'gs://aou_tmp/000_bgen_logs/{analysis_type}_analysis_export_{str(interval)}_bgen_{ancestry}.log')
                        bgen_task.attributes["ancestry"] = ancestry
                        bgen_task.attributes["analysis_type"] = analysis_type

                        bgen_index = index_bgen(b=b,
                                                ancestry=ancestry,
                                                analysis_type=analysis_type,
                                                bgen=f"{bgen_root}.bgen",
                                                depend_job=bgen_task)
                        bgen_index.attributes["ancestry"] = ancestry
                        bgen_index.attributes["analysis_type"] = analysis_type
                    if ((f"{bgen_root}.gene.txt" not in bgens_already_created)  or args.overwrite_gene_txt)and analysis_type=='gene':
                        gene_txt_task = b.new_python_job(
                            name=f"{analysis_type}_analysis_export_{str(interval)}_gene_txt_{ancestry}"
                        )
                        gene_txt_task.image(HAIL_DOCKER_IMAGE)
                        gene_txt_task.call(
                            export_gene_group_file,
                            interval=interval,
                            ancestry=ancestry,
                            output_dir=f"{bgen_root}.gene.txt",
                        )
                        gene_txt_task.attributes["ancestry"] = ancestry
                        gene_txt_task.attributes["analysis_type"] = analysis_type
                    if (not hfs.exists(f"{bgen_root}.bgen.bgi")) and  (hfs.exists(f"{bgen_root}.bgen")):
                        bgen_index=index_bgen(b=b, ancestry=ancestry, analysis_type=analysis_type, bgen=f"{bgen_root}.bgen")
                        bgen_index.attributes["ancestry"] = ancestry
                        bgen_index.attributes["analysis_type"] = analysis_type



                    bgen_file = b.read_input_group(
                        **{
                            "bgen": f"{bgen_root}.bgen",
                            "bgen.bgi": f"{bgen_root}.bgen.bgi",
                            "sample": f"{bgen_root}.sample",
                        }
                    )
                    if analysis_type == 'gene':
                        group_file = b.read_input(f"{bgen_root}.gene.txt")
                        bgens[str(interval)] = (bgen_file, group_file)
                    else:
                        bgens[str(interval)] = bgen_file

            elif not args.skip_saige:
                for interval in intervals:
                    bgen_root = f"{bgen_dir}/{analysis_type}_{interval.start.contig}_{str(interval.start.position).zfill(9)}_{interval.end.position}"
                    bgen_file = b.read_input_group(
                        **{
                            "bgen": f"{bgen_root}.bgen",
                            "bgen.bgi": f"{bgen_root}.bgen.bgi",
                            "sample": f"{bgen_root}.sample",
                        }
                    )
                    if analysis_type == 'gene':
                        group_file = b.read_input(f"{bgen_root}.gene.txt")
                        bgens[str(interval)] = (bgen_file, group_file)
                    else:
                        bgens[str(interval)] = bgen_file


            result_dir = f'{RESULT_ROOT}/result/{ancestry.upper()}'
            print(f'result directory: {result_dir}')
            overwrite_results = args.overwrite_results
            saige_tasks = {}
            saige_ur_tasks = {}
            if not args.skip_saige:
                memory = '10G' if ancestry in ['afr', 'eur', 'all'] and analysis_type == 'variant' else 'standard'
                groups=None
                if analysis_type == "gene":
                    if args.groups is None:
                        groups = ','.join(["pLoF", "missenseLC", "synonymous"])
                    else:
                        groups = args.groups
                    print(f'Groups to run: {groups}')
                print(f"------------{ancestry.upper()} {analysis_type} analysis step 2: {len(phenos_to_run)} phenotypes------------")
                for i in tqdm(range(len(phenos_to_run))):
                    phenoname = list(phenos_to_run)[i]
                    if phenoname in gate_traits:
                        current_trait_type = 'survival'
                    elif phenoname in quantitative_traits:
                        current_trait_type = 'quantitative'
                    elif phenoname in binary_traits:
                        current_trait_type = 'binary'
                    else:
                        raise ValueError(f"Unknown trait type for {phenoname}")

                    saige_tasks[phenoname] = []
                    if analysis_type == 'gene':
                        saige_ur_tasks[phenoname] = []
                    if phenoname not in null_models.keys():
                        continue

                    model_file, variance_ratio_file = null_models[phenoname]

                    if not i % 10:
                        n_jobs = dict(Counter(map(lambda x: x.name, b.select_jobs("")))).get("run_saige", 0)
                        logger.info(f"Read {i} phenotypes ({n_jobs} new to run so far)...")

                    pheno_results_dir = f"{result_dir}/phenotype_{phenoname}"
                    results_already_created = {}
                    exome_results_already_created = {}

                    if (
                        not overwrite_results
                        and hfs.exists(pheno_results_dir)
                    ):
                        results_already_created = {x["path"] for x in hl.hadoop_ls(pheno_results_dir)}
                    if (
                            not args.overwrite_exome_results
                            and hfs.exists(pheno_results_dir) and analysis_type == 'gene'
                    ):
                        exome_results_already_created = {x["path"] for x in hl.hadoop_ls(pheno_results_dir) if
                                                        x["path"].endswith('_exome.single_variant.txt')}
                        print(
                            f'-------------- Found {int(len(exome_results_already_created))} exome results in [{ancestry}: {phenoname}] directory...')
                    factor = 2 if analysis_type == 'variant' else 5
                    print(f'-------------- Found {int(len(results_already_created) / factor)} {analysis_type} results in [{ancestry}: {phenoname}] directory...')

                    for interval in intervals:
                        results_path = f"{pheno_results_dir}/result_{phenoname}_{interval.start.contig}_{str(interval.start.position).zfill(9)}"
                        group_file = None
                        max_maf_for_group_test = None
                        sparse_grm_file = None
                        suffix = ''
                        if analysis_type == "gene":
                            bgen_file, group_file = bgens[str(interval)]
                            max_maf_for_group_test = args.max_maf_group
                            sparse_grm_file = sparse_grm[sparse_grm_extension]
                        else:
                            bgen_file = bgens[str(interval)]
                        samples_file = b.read_input(f'{DATA_PATH}/utils/grm/{ancestry.upper()}_grm_plink.samples')

                        if (
                            overwrite_results
                            or (f"{results_path}.{'gene' if analysis_type == 'gene' else 'single_variant'}.txt"
                            not in results_already_created)
                        ):
                            saige_task = run_saige(
                                p=b,
                                phenoname=phenoname,
                                ancestry=ancestry,
                                output_root=results_path,
                                model_file=model_file,
                                variance_ratio_file=variance_ratio_file,
                                sparse_grm_file=sparse_grm_file,
                                bgen_file=bgen_file,
                                samples_file=samples_file,
                                docker_image=SAIGE_DOCKER_IMAGE,
                                group_file=group_file,
                                groups=groups,
                                trait_type=current_trait_type,
                                chrom=interval.start.contig,
                                min_mac=1,
                                min_maf=0,
                                max_maf_for_group=max_maf_for_group_test,
                                add_suffix=suffix,
                                variant_type=variant_type,
                                memory=memory
                            )
                            saige_task.attributes.update(
                                {"interval": str(interval), "ancestry": ancestry, 'name': f'saige{analysis_type == "gene"}'}
                            )
                            saige_task.attributes.update(
                                {"phenotype": copy.deepcopy(phenoname)}
                            )
                            saige_tasks[phenoname].append(saige_task)

                        if analysis_type == 'gene':
                            if (
                                    args.overwrite_exome_results
                                    or (f"{results_path}_exome.single_variant.txt" not in exome_results_already_created)
                            ):
                                memory = 'highmem' if ancestry in ['eur', 'all']  else 'standard'
                                saige_ur_task = run_saige(
                                    p=b,
                                    phenoname=phenoname,
                                    ancestry=ancestry,
                                    output_root=f'{results_path}_exome',
                                    model_file=model_file,
                                    variance_ratio_file=variance_ratio_file,
                                    sparse_grm_file=None,
                                    bgen_file=bgen_file,
                                    samples_file=samples_file,
                                    docker_image=SAIGE_DOCKER_IMAGE,
                                    group_file=None,
                                    groups=None,
                                    trait_type=current_trait_type,
                                    chrom=interval.start.contig,
                                    min_mac=1,
                                    min_maf=0,
                                    memory=memory,
                                    max_maf_for_group=None,
                                    add_suffix=suffix,
                                    variant_type=variant_type,
                                )
                                saige_ur_task.attributes.update(
                                    {"interval": str(interval), "ancestry": ancestry}
                                )
                                saige_ur_task.attributes.update(
                                    {"phenotype": copy.deepcopy(phenoname)}
                                )
                                saige_ur_tasks[phenoname].append(saige_ur_task)
                    if args.test:
                        break

            if not args.skip_load_hail_results:
                results_type = 'ACAF' if analysis_type == 'variant' else 'Exome'
                root = f'{ANALYSIS_BUCKET}/ht_results/{ancestry.upper()}'
                print(f'--------------Loading results from {results_type} analysis [{ancestry.upper()}]: {root}------------')
                for i in tqdm(range(len(phenos_to_run))):
                    phenoname = list(phenos_to_run)[i]
                    if phenoname in gate_traits and not ancestry in phenoname:
                        continue
                    directory = f'{ANALYSIS_BUCKET}/{analysis_type}_results/result/{ancestry.upper()}/phenotype_{phenoname}'
                    output_ht_directory = f'{root}/phenotype_{phenoname}'
                    null_glmm_log = f'{ANALYSIS_BUCKET}/{analysis_type}_results/null_glmm/{ancestry.upper()}/phenotype_{phenoname}.log'

                    saige_log = f'{directory}/result_{phenoname}_chr1_{"000065419" if analysis_type == "gene" else "000000001"}.{analysis_type}.log'

                    quantitative_trait = phenoname in quantitative_traits

                    if not args.skip_load_gene_results and analysis_type == 'gene':
                        if (not hfs.exists(f'{output_ht_directory}/gene_results.ht/_SUCCESS')) or args.overwrite_hail_results:
                            print(f'[{ancestry.upper()}: {phenoname}] Gene table')
                            j = b.new_python_job(
                                name=f'sync_saige_{analysis_type}_gene_HT_{phenoname}_{ancestry}',
                                attributes={"analysis_type": analysis_type, "ancestry": ancestry,
                                            "phenotype": copy.deepcopy(phenoname)})
                            j.image("hailgenetics/hail:0.2.133-py3.11")
                            j.memory('standard')
                            j.cpu(16)
                            j.env('PYSPARK_SUBMIT_ARGS', '--driver-memory 8g --executor-memory 8g pyspark-shell')
                            j.call(load_gene_data,
                                   directory=directory,
                                   output_ht_directory=output_ht_directory,
                                   phenoname=phenoname,
                                   gene_ht_map_path=get_aou_gene_map_ht_path(pop=ancestry, processed=False),
                                   quantitative_trait=quantitative_trait,
                                   null_glmm_log=null_glmm_log,
                                   saige_log=saige_log,
                                   overwrite=True
                                   )
                            if phenoname in list(saige_tasks.keys()) and not args.skip_saige:
                                j.depends_on(*saige_tasks[phenoname])
                            j.always_run()

                    if not args.skip_load_variant_results:
                        extension = 'single_variant.txt'
                        variant_type = 'genome' if analysis_type == 'variant' else 'exome'
                        j = b.new_python_job(
                            name=f'sync_saige_{variant_type}_variant_HT_{phenoname}_{ancestry}',
                            attributes={"analysis_type": analysis_type, "ancestry": ancestry, "phenotype": copy.deepcopy(phenoname)})
                        j.image("hailgenetics/hail:0.2.133-py3.11")
                        j.memory('highmem')
                        j.cpu(16)
                        j.env('PYSPARK_SUBMIT_ARGS', '--driver-memory 16g --executor-memory 16g pyspark-shell')
                        SUCCESS = hfs.exists(f'{output_ht_directory}/{variant_type}_variant_results.ht/_SUCCESS')

                        print(f'[{ancestry.upper()}: {phenoname}] {variant_type} variant table')
                        j.call(load_variant_data,
                                   directory=directory,
                                   output_ht_directory=output_ht_directory,
                                   phenoname=phenoname,
                                   quantitative_trait=quantitative_trait,
                                   null_glmm_log=null_glmm_log,
                                   extension=extension,
                                   overwrite=args.overwrite_hail_results or not SUCCESS,
                                   variant_type = variant_type
                                   )
                        if variant_type == 'genome':
                            if phenoname in list(saige_tasks.keys()) and not args.skip_saige:
                                j.depends_on(*saige_tasks[phenoname])
                        else:
                            if phenoname in list(saige_ur_tasks.keys()) and not args.skip_saige:
                                j.depends_on(*saige_ur_tasks[phenoname])
                        j.always_run()
                    if args.test:
                        break

            def get_tasks_from_pipeline(p):
                return dict(Counter(map(lambda x: x.name, p.select_jobs(""))))

            logger.info(
                f'Setup took: {time.strftime("%H:%M:%S", time.gmtime(time.time() - start_time))}'
            )
            logger.info(f"Submitting: {get_tasks_from_pipeline(b)}")
            logger.info(f"Total size: {sum([len(x._pretty()) for x in b.select_jobs('')])}")
            logger.info(f"Finished: {get_tasks_from_pipeline(b)}")
            b.run(delete_scratch_on_exit=False)



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--single-variant-only", help="Run single variant SAIGE", action="store_true"
    )
    parser.add_argument("--run-pipeline", help="Run SAIGE pipeline", action="store_true")
    parser.add_argument(
        "--overwrite-pheno-data", help="Overwrite phenotype data", action="store_true"
    )
    parser.add_argument(
        "--skip-any-null-models",
        help="Skip running SAIGE null models",
        action="store_true",
    )
    parser.add_argument(
        "--skip-saige", help="Skip running SAIGE tests", action="store_true"
    )
    parser.add_argument(
        "--skip-bgen", help="Skip generating bgen files", action="store_true"
    )
    parser.add_argument(
        "--skip-any-random-pheno", help="Skip phenotypes", action="store_true"
    )
    parser.add_argument(
        "--export-phenos", help="Export phenotype tsv files", action="store_true"
    )
    parser.add_argument(
        "--skip-load-hail-results", help="Skip loading results into hail", action="store_true"
    )
    parser.add_argument(
        "--overwrite-null-models",
        help="Force creation of null models",
        action="store_true",
    )
    parser.add_argument(
        "--create-bgens", help="Force creation of Bgen files", action="store_true"
    )
    parser.add_argument(
        "--overwrite-bgens", help="Force run of Bgens", action="store_true"
    )
    parser.add_argument(
        "--overwrite-gene-txt", help="Force run of gene txt files", action="store_true"
    )
    parser.add_argument(
        "--overwrite-results", help="Force run of SAIGE tests", action="store_true"
    )
    parser.add_argument(
        "--overwrite-dict", help="Force creation of phenotype dictionary", action="store_true"
    )
    parser.add_argument(
        "--overwrite-exome-results", help="Force run of SAIGE exome variant tests", action="store_true"
    )
    parser.add_argument(
        "--overwrite-hail-results",
        help="Force run of results loading",
        action="store_true",
    )
    parser.add_argument(
        "--skip-load-variant-results", help="Skip loading single-variant level test results", action="store_true"
    )
    parser.add_argument(
        "--skip-load-gene-results", help="Skip loading gene level test results", action="store_true"
    )
    parser.add_argument(
        "--test",
        help="Test run of pipeline",
        action="store_true",
    )
    parser.add_argument(
        "--non-pre-emptible", help="Local test of pipeline", action="store_true"
    )
    parser.add_argument(
        "--phenos",
        help="Comma-separated list of phenotype names"
    )
    parser.add_argument(
        "--groups",
        help="Comma-separated list of functional groups used for gene-based test ",
        default="pLoF,missenseLC,synonymous,pLoF:missenseLC",
    )
    parser.add_argument(
        "--ancestries", help="comma-separated list", default="afr,amr,eas,eur,mid,sas"
    )
    parser.add_argument(
        "--pilot", help="Run pilot phenotypes only", action="store_true"
    )
    parser.add_argument(
        "--any-random-pheno", help="if any random phenotypes are being run", action="store_true"
    )
    parser.add_argument(
        "-hail-image", help="hail docker image to use", action="store_true",
        # default='us-central1-docker.pkg.dev/broad-mpg-gnomad/images/vrs084'
        default ='hailgenetics/hail:0.2.133-py3.11'
    )
    parser.add_argument(
        "--callrate-filter",
        help="Impose filter of specified callrate (default: 0) when exporting bgen",
        default=0.0,
        type=float,
    )
    parser.add_argument(
        "--variant-ac-filter",
        help="Impose filter of specified variant AC (default: 0) when exporting bgen",
        default=0,
        type=int,
    )
    parser.add_argument(
        "--max-maf-group",
        help="Max MAF of variants to include into the group tests",
        default='0.01,0.001,0.0001',
    )
    parser.add_argument(
        "--category",
        help="phenotype category to export",
    )
    args = parser.parse_args()

    main(args)
