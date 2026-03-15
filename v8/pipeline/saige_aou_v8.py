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
import time
from typing import *
from tqdm import tqdm
import warnings
from datetime import datetime
import math
import sys

ACAF_MT_PATH = 'gs://fc-aou-datasets-controlled/v8/wgs/short_read/snpindel/acaf_threshold/splitMT/hail.mt'
EXOME_MT_PATH = 'gs://fc-aou-datasets-controlled/v8/wgs/short_read/snpindel/exome/splitMT/hail.mt'

timestamp = datetime.now().strftime("%y%m%d_%H%M%S")

def range_table(log_file):
    import hail as hl
    hl.init(log=log_file)
    print(hl.utils.range_table(10)._force_count())


TRANCHE = "v8"
MY_BUCKET = 'gs://aou_wlu'
ANALYSIS_BUCKET = "gs://aou_analysis/v8"
DATA_PATH = f'{ANALYSIS_BUCKET}/data'
TMP_BUCKET = 'gs://aou_tmp/v8'

N_GENE_PER_GROUP = {'all': 25, 'eur': 100, 'afr': 100, 'amr': 100, 'eas': 200, 'mid': 200, 'sas': 200}
# CHUNK_SIZE = {'all': int(6.25e6),'eur': int(1.25e7), "afr": int(2.5e7), "amr": int(2.5e7), "eas": int(2.5e7), "mid": int(2.5e7), "sas": int(2.5e7)} # new for pgen
CHUNK_SIZE = {'all': int(1.6e6),'eur': int(3.125e6), "afr": int(1.25e7), "amr": int(1.25e7), "eas": int(2.5e7), "mid": int(2.5e7), "sas": int(2.5e7)} # new for smaller EUR and ALL bigger others 
REFERENCE = "GRCh38"
CHROMOSOMES = list(map(str, range(1, 23))) + ["X", "Y"]


PILOT_PHENOTYPES = ['EM_202.2', 'height']

PHENO_CATEGORIES = ['physical_measurement', 'r_drug', 'pfhh_survey', 'random_pheno', 'lab_measurement', 'mcc2_phecode', 'mcc2_phecodex', 'onset', 'progression', 'mhwb_survey']
PHENO_CATEGORIES_MT = ['physical_measurement', 'r_drug', 'pfhh_survey', 'mcc2_phecode', 'mcc2_phecodex', 'mhwb_survey']
PHENO_CATEGORIES_HT = ['lab_measurement', 'onset', 'progression']
GATE_CATEGORIES = ['onset', 'progression']
BINARY_CATEGORIES = ['r_drug', 'pfhh_survey', 'mcc2_phecode', 'mcc2_phecodex', 'mhwb_survey']
QUANTITATIVE_CATEGORIES = ['physical_measurement', 'lab_measurement']
ANCESTRIES = ['AFR', 'AMR', 'EAS', 'EUR', 'MID', 'SAS']    
N_SAMPLES_PRUNED = {'afr': 77444, 'amr': 71540, 'eas': 9488, 'eur': 227273, 'mid': 1153, 'sas': 5132, 'all': 407490}


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

def get_filtered_mt(mt_type: str,
                    sample_ids: hl.Table,
                    ancestry: str='all',
                    filter_samples: bool=True, 
                    filter_variants: bool=True,
                    prune_samples:bool=True,
                    adj_filter: bool=True):
    anc_mt_path = f'{DATA_PATH}/utils/raw_mt/{mt_type}/{ancestry.upper()}_{mt_type}.mt'
    if not hfs.exists(f'{anc_mt_path}/_SUCCESS'): 
        if ancestry != 'all':
            ancestry_ht = hl.read_table(f'{DATA_PATH}/utils/aou_v8_global_pca.ht')
            ancestry_ht = ancestry_ht.filter(ancestry_ht.ancestry_pred == ancestry)
            mt_path = ACAF_MT_PATH if mt_type=='ACAF' else EXOME_MT_PATH
            print(mt_path)
            mt = hl.read_matrix_table(mt_path)
            mt.describe()
            print(mt.count())
            mt = mt.filter_entries(hl.is_missing(mt.FT) | (mt.FT == 'PASS'))
            mt = mt.filter_cols(hl.is_defined(ancestry_ht[mt.col_key]))
            mt = mt.naive_coalesce(20000).checkpoint(anc_mt_path, overwrite=True)
        else:
            mt_path = ACAF_MT_PATH if mt_type=='ACAF' else EXOME_MT_PATH
            mt = hl.read_matrix_table(mt_path)
            mt = mt.filter_entries(hl.is_missing(mt.FT) | (mt.FT == 'PASS'))
            mt = mt.naive_coalesce(20000).checkpoint(anc_mt_path, overwrite=True)
    mt = hl.read_matrix_table(anc_mt_path)
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

        print('Filtering out chrY from Female samlples....')
        demo_ht = hl.read_table('gs://aou_analysis/v8/data/utils/aou_v8_demographic.ht')
        mt = mt.annotate_cols(sex_at_birth = demo_ht[mt.col_key].sex_at_birth)
        mt = mt.annotate_entries(female_chrY = (mt.sex_at_birth == 'Female') & (mt.locus.contig == 'chrY'))
        mt = mt.filter_entries(~mt.female_chrY)
    return mt

logging.basicConfig(
    format="%(levelname)s (%(name)s %(lineno)s): %(message)s",
    level="INFO",
    filename="saige_pipeline.log",
)
logger = logging.getLogger("ALL_x_AoU_SAIGE")
logger.setLevel(logging.INFO)

HAIL_DOCKER_IMAGE = "hailgenetics/hail:0.2.133-py3.11"
SAIGE_DOCKER_IMAGE = "wzhou88/saige:1.5.0.2"  # latest
QQ_DOCKER_IMAGE = "konradjk/saige_qq:0.2"

def annotate_expected_pvalue(
        ht: hl.Table, 
        method:str, 
        p_field: str='Pvalue', 
        k:int = 5000,
        quantiles: List[float] = [0.7, 0.5, 0.1, 0.01, 0.001]
    ):
    n = ht.count()
    print(f"Total variants: {n}")
    ht = ht.filter(hl.is_defined(ht[p_field]))
    keys = list(ht.key)
    if method == 'exact':
        ht = ht.key_by()
        ht = ht.order_by(ht[p_field]).add_index()
        ht = ht.annotate(idx = ht.idx+1)
        ht = ht.key_by(*keys)
        ht = ht.annotate(**{f'{p_field}_expected': ht.idx / (n + 1)})
        expected_p = []
        for q in quantiles:
            pos = q * n
            lo = int(math.floor(pos))
            hi = int(math.ceil(pos))
            if lo == hi:
                p_at_q = ht.filter(ht.idx == pos)
            else:
                p_at_q = ht.filter((ht.idx == lo) | (ht.idx == hi))
            p_at_q = p_at_q.aggregate(hl.agg.mean(p_at_q[p_field]))
            expected_p.append(p_at_q)
        quantile_p_dict = dict(zip(quantiles, expected_p))
        quantile_lambda = [hl.eval(hl.qchisqtail(p=p, df=1, ncp=0, lower_tail=False) /
                        hl.qchisqtail(p=q, df=1, ncp=0, lower_tail=True)) for p, q in zip(expected_p, quantiles)]
        lambda_fields = {
            f'lambda_{p_field}_q{quantiles[i]}': quantile_lambda[i]
            for i in range(len(quantiles)) 
        }
        ht = ht.annotate_globals(**lambda_fields)
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
                    ht[p_field] - ht.lower_pvalue) + ht.lower_rank)))
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
        quantile_indices = [hl.eval(hl.floor(n * q)) for q in quantiles]
        quantile_upper_idx = [bisect(hl.eval(ht.ranks[1]), idx) for idx in quantile_indices]
        quantile_upper_rank = [hl.eval(ht.ranks[1])[idx] for idx in quantile_upper_idx]
        quantile_lower_rank = [hl.eval(ht.ranks[1])[idx - 1] for idx in quantile_upper_idx]
        quantile_upper_pvalue = [hl.eval(ht.ranks[0])[idx] for idx in quantile_upper_idx]
        quantile_lower_pvalue = [hl.eval(ht.ranks[0])[idx - 1] for idx in quantile_upper_idx]
        quantile_p = [((hl.eval(hl.floor(n * quantile)) - lower_rank) / 
               (upper_rank - lower_rank)) * 
               (upper_pvalue - lower_pvalue) + lower_pvalue 
              for quantile, upper_rank, lower_rank, upper_pvalue, lower_pvalue 
              in zip(quantiles, quantile_upper_rank, quantile_lower_rank, quantile_upper_pvalue, quantile_lower_pvalue)]
        quantile_p_dict = dict(zip(quantiles, quantile_p))
        quantile_lambda = [hl.eval(hl.qchisqtail(p=p, df=1, ncp=0, lower_tail=False) /
                        hl.qchisqtail(p=q, df=1, ncp=0, lower_tail=True)) for p, q in zip(quantile_p, quantiles)]
        lambda_fields = {
            f'lambda_{p_field}_q{quantiles[i]}': quantile_lambda[i]
            for i in range(len(quantiles)) 
        }
        ht = ht.annotate_globals(**lambda_fields)
        ht = ht.drop('ranks', 'rank')

    else:
        raise ValueError("`method` must be 'exact' or 'approx_cdf'")
    
    ht = ht.annotate_globals(**lambda_fields)
    lambda_gc = hl.eval(hl.qchisqtail(p=quantile_p_dict[0.5], df=1, ncp=0, lower_tail=False) /
                        hl.qchisqtail(p=0.5, df=1, ncp=0, lower_tail=True))
    print(f"Lambda GC: {lambda_gc}")

    ht = ht.annotate_globals(**{f'lambda_gc_{p_field}': lambda_gc})
    ht = ht.annotate(**{f'{p_field}': hl.if_else(ht[f'{p_field}'] == 0, 1e-320, ht[f'{p_field}'])}, )
    ht = ht.annotate(**{f'{p_field}_expected_log10': -hl.log10(ht[f'{p_field}_expected'])},)

    return ht

def create_expected_p_ht(pheno:str, ancestry:str, table_name:str, analysis_type:str, method:str, overwrite:bool, misc:bool=False):
    # if ancestry == 'meta':
        # hl.init(backend='batch', app_name=f'expected_p_{pheno}_{ancestry}_{method}_{table_name}_QoBiB')
    if '_to_' in pheno:
        software_type = 'gate_'
    else:
        software_type = ''
    sex_tag = ''
    if pheno.endswith('male') or pheno.endswith('female'):
        sex_tag = f"_{pheno.split('_')[-1]}"
    path = f'{ANALYSIS_BUCKET}/{software_type}ht_results{sex_tag}{"" if not misc else "_misc"}/{ancestry.upper()}/phenotype_{pheno}/{table_name}.ht'
    success = path.replace('.ht', f'_{method}_expected_p.ht/_SUCCESS')
    if (hl.hadoop_exists(path) and not hl.hadoop_exists(success)) or overwrite:
        overwrite=True
        ht = hl.read_table(path)
        ht.describe()
        if analysis_type == 'variant':
            ht = annotate_expected_pvalue(ht=ht, method=method, p_field='Pvalue')
            ht.describe()
            ht = ht.filter(hl.case()
                            .when(ht.Pvalue > 0.1, hl.rand_bool(0.001))
                            .when(ht.Pvalue > 0.01, hl.rand_bool(0.01))
                            .when(ht.Pvalue > 0.001, hl.rand_bool(0.1))
                            .default(True))
            ht = ht.checkpoint(path.replace('.ht', f'_{method}_expected_p.ht'), overwrite=overwrite,
                            _read_if_exists=not overwrite)
            ht.export(path.replace('.ht', f'_{method}_expected_p.txt.bgz'))
        else:
            ht = ht.key_by('gene_id', 'gene_symbol', 'annotation')
            ht = ht.annotate(max_MAF=hl.if_else(hl.is_missing(ht.max_MAF), -1, ht.max_MAF))
            for max_MAF in [0.01, 0.001, 0.0001, -1]:
                tag = 'Cauchy' if max_MAF == -1 else str(max_MAF)
                sub_ht = ht.filter(ht.max_MAF == max_MAF)
                sub_ht = annotate_expected_pvalue(ht=sub_ht, method=method, p_field='Pvalue')
                sub_ht = annotate_expected_pvalue(ht=sub_ht, method=method, p_field='Pvalue_Burden')
                sub_ht = annotate_expected_pvalue(ht=sub_ht, method=method, p_field='Pvalue_SKAT')
                sub_ht.describe()
                sub_ht = sub_ht.transmute_globals(**{f'lambda_gc_maxmaf_{tag}':
                                                    hl.struct(lambda_gc_Pvalue = sub_ht.lambda_gc_Pvalue,
                                                            lambda_gc_Pvalue_Burden=sub_ht.lambda_gc_Pvalue_Burden,
                                                            lambda_gc_Pvalue_SKAT=sub_ht.lambda_gc_Pvalue_SKAT,
                                                            )})
                sub_ht = sub_ht.checkpoint(path.replace('.ht', f'_{method}_expected_p_{tag}.ht'), overwrite=overwrite,
                                    _read_if_exists=not overwrite)
                sub_ht.export(path.replace('.ht', f'_{method}_expected_p_{tag}.txt.bgz'))
            ht = hl.read_table(path.replace('.ht', f'_{method}_expected_p_0.01.ht'))
            for max_MAF in [0.001, 0.0001, -1]:
                tag = 'Cauchy' if max_MAF == -1 else str(max_MAF)
                tmp_ht = hl.read_table(path.replace('.ht', f'_{method}_expected_p_{tag}.ht'))
                ht = ht.union(tmp_ht)
                ht = ht.annotate_globals(
                    **{f'lambda_gc_maxmaf_{tag}': tmp_ht.index_globals()[f'lambda_gc_maxmaf_{tag}']})
            ht = ht.key_by('gene_id', 'gene_symbol', 'annotation', 'max_MAF')
            ht = ht.checkpoint(path.replace('.ht', f'_{method}_expected_p.ht'), overwrite=overwrite,
                                    _read_if_exists=not overwrite)
    return ht

def load_gene_data(directory: str,
                   output_ht_directory: str,
                   phenoname, 
                   ancestry,
                   gene_ht_map_path: str, 
                   quantitative_trait:bool,
                   null_glmm_log: str,
                   overwrite: bool = False,
                   misc: bool = False):

    hl.init(
        master='local[8]',
        tmp_dir=TMP_BUCKET,
        default_reference="GRCh38",
    )

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
    output_ht_path = f'{output_ht_directory}/gene_results.ht'
    print(f'Loading: {directory}/*.gene.txt ...')
    types = {x: hl.tint32 for x in ('MAC',	'Number_rare', 'Number_ultra_rare')}
    types.update({x: hl.tfloat64 for x in ('max_MAF', 'BETA_Burden', 'SE_Burden')})
    types.update({x: hl.tstr for x in ('Pvalue', 'Pvalue_Burden', 'Pvalue_SKAT')})
    if not phenoname.endswith('male'):
        single_variant_ht = hl.import_table(f'{directory}/result_{phenoname}_chr1_000065419.result.singleAssoc.txt', delimiter='\t', impute=True, types = {'p.value':hl.tstr})
    else:
        single_variant_ht = hl.import_table(f'{directory}/result_*_chrX_000276322.result.singleAssoc.txt', delimiter='\t', impute=True, types = {'p.value':hl.tstr})
    if 'N' in list(single_variant_ht.row_value):
        n_cases = list(single_variant_ht.aggregate(hl.agg.collect_as_set(single_variant_ht.N)))[0]
        n_controls = 0
        N = n_cases
    else:
        n_cases = list(single_variant_ht.aggregate(hl.agg.collect_as_set(single_variant_ht.N_case)))[0]
        n_controls = list(single_variant_ht.aggregate(hl.agg.collect_as_set(single_variant_ht.N_ctrl)))[0]
        N = n_cases + n_controls

    ht = hl.import_table(f'{directory}/*.gene.txt', delimiter='\t', impute=True, types=types)
    ht.describe()
    if heritability == -1.0: heritability = hl.null(hl.tfloat)
    if saige_version == 'NA': saige_version = hl.null(hl.tstr)
    if inv_normalized == 'NA': inv_normalized = hl.null(hl.tstr)

    ht = ht.annotate(**{f'{p_field}_log10': parse_log_p_value(ht[p_field]) for p_field in ('Pvalue', 'Pvalue_Burden', 'Pvalue_SKAT')})
    ht = ht.annotate(**{f'{p_field}': hl.float64(ht[p_field]) for p_field in ('Pvalue', 'Pvalue_Burden', 'Pvalue_SKAT')})
    fields = ht.Region.split('_')
    gene_ht = hl.read_table(gene_ht_map_path).select('interval').distinct()
    ht = ht.key_by(gene_id=fields[0], gene_symbol=fields[1], annotation=ht.Group, phenoname=phenoname, max_MAF=ht.max_MAF).drop('Region', 'Group').naive_coalesce(10).annotate_globals(
        heritability=heritability, saige_version=saige_version, inv_normalized=inv_normalized)
    ht = ht.annotate(total_variants= ht.Number_rare + ht.Number_ultra_rare,
                     interval=gene_ht.key_by('gene_id')[ht.gene_id].interval)
    ht = ht.annotate(CHR= ht.interval.start.contig, POS=ht.interval.start.position)
    ht = ht.annotate_globals(n_cases=n_cases, n_controls=n_controls, N=N)
    ht.describe()
    ht = ht.checkpoint(output_ht_path, overwrite=overwrite, _read_if_exists=not overwrite)

    if not phenoname.endswith('male') or phenoname.startswith('EM_202.2') or phenoname.startswith('EM_236.1') or phenoname.startswith('A10BJ'):
        method = 'exact'
        success = output_ht_path.replace('.ht', f'_{method}_expected_p.ht/_SUCCESS')
        if not hl.hadoop_exists(output_ht_path):
            sys.exit(f'Output HT {output_ht_path} does not exist.')
        elif (hl.hadoop_exists(output_ht_path) and not hl.hadoop_exists(success)) or overwrite:
            ht = create_expected_p_ht(pheno=phenoname, ancestry=ancestry, table_name='gene_results', analysis_type='gene', method='exact', overwrite=True, misc=misc)
        else:
            ht = hl.read_table(output_ht_path.replace('.ht', f'_{method}_expected_p.ht'))
    ht.show()
    print(ht.count())
    # print(ht.aggregate(hl.agg.counter(ht.annotation)))
    return ht

def load_variant_data(directory: str,
                      output_ht_directory: str,
                      phenoname: str,
                      ancestry: str,
                      quantitative_trait:bool,
                      null_glmm_log: str,
                      extension: str = 'single.txt',
                      overwrite: bool = False,
                      variant_type: str = 'genome',
                      num_partitions: int = 1000,
                      misc: bool = False):
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
            n_cases = list(ht.aggregate(hl.agg.collect_as_set(ht.N)))[0]
            n_controls = 0
            N = n_cases
            ht = ht.annotate(
                **{'p.value.NA': hl.missing(hl.tfloat64),
                'Is.SPA': hl.missing(hl.tbool),
                'AF_case': hl.missing(hl.tfloat64),
                'AF_ctrl': hl.missing(hl.tfloat64),
                }
            )
        elif 'N_event' in list(ht.row_value):
            n_cases = list(ht.aggregate(hl.agg.collect_as_set(ht.N_event)))[0]
            n_controls = list(ht.aggregate(hl.agg.collect_as_set(ht.N_censor)))[0]
            N = n_cases + n_controls
            ht = ht.annotate(
                AF_case = ht.AF_event,
                AF_ctrl = ht.AF_censor,
            )
        else:
            n_cases = list(ht.aggregate(hl.agg.collect_as_set(ht.N_case)))[0]
            n_controls = list(ht.aggregate(hl.agg.collect_as_set(ht.N_ctrl)))[0]
            N = n_cases + n_controls

        ht = ht.key_by(locus=hl.locus(ht.CHR, ht.POS, reference_genome = 'GRCh38'), alleles=[ht.Allele1, ht.Allele2],
                    phenoname=phenoname).distinct().naive_coalesce(num_partitions)
        ht = ht.transmute(Pvalue=ht['p.value']).annotate_globals(
            heritability=heritability, saige_version=saige_version,
            inv_normalized=inv_normalized)
        ht = ht.drop('Allele1', 'Allele2')
        ht = ht.annotate(Pvalue_log10 = parse_log_p_value(ht.Pvalue))
        ht = ht.annotate(Pvalue = hl.float64(ht.Pvalue))
        ht = ht.select('CHR', 'POS', 'MarkerID', 'AC_Allele2', 'AF_Allele2', 'MissingRate', 'BETA', 'SE', 'Tstat', 'var', 
        'Pvalue', 'Pvalue_log10', 'p.value.NA', 'Is.SPA', 'AF_case', 'AF_ctrl')
        ht = ht.annotate_globals(n_cases=n_cases, n_controls=n_controls, N=N)
        ht = ht.checkpoint(output_ht_path, overwrite=overwrite, _read_if_exists=not overwrite)
    ht = hl.read_table(output_ht_path)
    ht.describe()
    ht.show()
    if not phenoname.endswith('male') or phenoname.startswith('EM_202.2') or phenoname.startswith('EM_236.1') or phenoname.startswith('A10BJ'):
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
            ht = create_expected_p_ht(pheno=phenoname, ancestry=ancestry, table_name=f'{variant_type}_variant_results', analysis_type='variant', method='approx_cdf', overwrite=True, misc=misc)
            ht = ht.annotate(rsID=rsid_ht[ht.locus, ht.alleles].rsid)
        else:
            ht = hl.read_table(output_ht_path.replace('.ht', f'_{method}_expected_p.ht'))
        ht.show()
    print(ht.count())
    return ht

def run_saige(
    p: Batch,
    phenoname: str,
    ancestry:str,
    output_root: str,
    model_file: str,
    variance_ratio_file: str,
    sparse_grm_file: str,
    geno_file: ResourceGroup,
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
    pgen: bool = False,
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
        f"set -o pipefail; {MKL_OFF} pixi run --manifest-path /app/pixi.toml Rscript /usr/local/bin/step2_SPAtests.R "
        # f"set -o pipefail; {MKL_OFF} Rscript /usr/local/bin/step2_SPAtests.R "
        f"--vcfField=GT "
        f"--maxMissing=0.2 "
        f"--minMAF={min_maf} " 
        f"--minMAC={min_mac} " 
        f"--sampleFile={samples_file} " 
        f"--GMMATmodelFile={model_file} " 
        f"--varianceRatioFile={variance_ratio_file} " 
        f"--AlleleOrder=ref-first " 
        f"--SAIGEOutputFile={run_saige_task.result} "
        f"--is_fastTest=FALSE "
        f"--is_noadjCov=FALSE "
    )
    if not pgen:
        command += f"--bgenFile={geno_file.bgen} " 
        command += f'--bgenFileIndex={geno_file["bgen.bgi"]} '
        command += f"--chrom={chrom} "  
    else:
        command += f"--pgenPrefix={geno_file} " 
        command += f"--chrom={chrom} "  

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
    sex_tag: str,
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
    :param sex_tag:
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
        fit_null_task._machine_type = "n1-highmem-32"
        fit_null_task._preemptible = True
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
        f"set -o pipefail; pixi run --manifest-path /app/pixi.toml Rscript /usr/local/bin/step1_fitNULLGLMM.R "
        # f"set -o pipefail; Rscript /usr/local/bin/step1_fitNULLGLMM.R "
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
        f"--usePCGwithSparseGRM=TRUE "
    )
    if sex_tag != '':
        command += f"--sexCol=sex "
        command += "--FemaleCode=0 "
        command += "--MaleCode=1 "
        if sex_tag == '_male':
            command += "--MaleOnly=TRUE "
        elif sex_tag == '_female':
            command += "--FemaleOnly=TRUE "
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
    bgen: bool=True,
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
    prune = True if ancestry.lower() != 'all' else False
    mt = get_filtered_mt(mt_type=mt_type, sample_ids=sample_ids, filter_variants=True, filter_samples=True, adj_filter= True, ancestry=ancestry, prune_samples=prune)  
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

    tag = 'pruned' if prune else 'pre_pruning'
    call_stats_ht = hl.read_table(f'{DATA_PATH}/utils/call_stats/{mt_type.lower() if mt_type == "Exome" else mt_type}_{tag}/{ancestry.upper()}_{mt_type.lower() if mt_type == "Exome" else mt_type}_call_stats.ht')
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
    if bgen:
        hl.export_bgen(mt, f"{output_dir}/{outname}", gp=mt.GP, varid=mt.rsid)
    else:
        hl.export_plink(mt, f"{output_dir}/{outname}" , ind_id = mt.s)
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

def index_vcf(b: hb.batch.Batch, vcf: str, depend_job=None):
    file = b.read_input(vcf)
    j = b.new_job(name=f'index_vcf')
    if depend_job is not None:
        j.depends_on(depend_job)
    j.image('us.gcr.io/broad-gatk/gatk:4.2.6.1')
    j.command(f'tabix -C -p vcf {file} > {j.temp}')
    b.write_output(j.temp, f'{vcf}.csi')
    return j

def from_bgen_to_pgen(
    batch: hb.batch.Batch,
    bgen_path: str,
    sample_path: str,
    allele_order: str = 'ref-first',
    depend_job=None,
    highmem: bool = False
):
    name = bgen_path.split('/')[-1]
    bgen_file = batch.read_input(bgen_path)
    sample_file = batch.read_input(sample_path)
    j = batch.new_job(name=f'bgen_to_pgen_{name}')
    j.image('us-central1-docker.pkg.dev/aou-neale-gwas/saige/saige_plink2:0.1')
    if depend_job is not None:
        j.depends_on(depend_job)
    if highmem:
        j.storage('20G')
        j.memory('65G')
    command = (f'plink2 --bgen {bgen_file} {allele_order} '
              f'--sample {sample_file} '
              f'--make-pgen '
              f'--out output ')
    if ('chrX_000276322_47659161' in name) or ('chrX_120625793_156010817' in name) or \
        ('chrX_150000001_156040895' in name) or ('chrX_000000001_25000001' in name) or \
            ('chrX_000000001_12500001' in name) or ('chrX_000276322_21654695' in name) or\
                ('chrX_139530758_156010817' in name) or ('chrX_000000001_6250001' in name) or \
                ('chrX_000276322_12977227' in name) or ('chrX_154182596_156010817' in name) or\
                    ('chrX_000000001_3125001' in name) or ('chrX_153125001_156040895' in name) or \
                        ('chrX_000276322_6535118' in name) or ('chrX_154778684_156010817' in name) or \
                            ('chrX_000000001_1600001' in name) or ('chrX_001600001_3200001' in name) or \
                                ('chrX_155200001_156040895' in name):
        command += f' --split-par b38 '
    j.command(command)
    j.command("awk 'BEGIN {OFS=\"\\t\"} NR==1 {print; next} {$1=\"chr\"$1; print}' output.pvar > output_chr.pvar")
    j.command(f'mv output.pgen {j.pgen}')
    j.command(f'mv output_chr.pvar {j.pvar}')
    j.command(f'mv output.psam {j.psam}')
    batch.write_output(j.pgen, f'{bgen_path.replace("bgen", "pgen")}')
    batch.write_output(j.pvar, f'{bgen_path.replace("bgen", "pgen").replace(".pgen", ".pvar")}')
    batch.write_output(j.psam, f'{bgen_path.replace("bgen", "pgen").replace(".pgen", ".psam")}')
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
            ht = hl.read_table(f"{DATA_PATH}/phenotype/{category}/{'lab_measurement_' if category == 'lab_measurement' else ''}{phenoname}.ht")
            ht = ht.filter(ht.ancestry == ancestry)
            if category == 'lab_measurement':
                ht = ht.annotate(value=ht['median_value'])
        ht = ht.checkpoint(f'gs://aou_tmp/export_pheno/{ancestry}_{phenoname}_{timestamp}.ht', overwrite = True)
        print(f"Exporting {phenoname} for {ancestry.upper()} (n_samples: {ht.count()})...")
        return ht

    def load_mega_pheno_ht(category:str, phenoname:str):
        user_id_col = "userId" if phenoname.startswith('random') else "person_id"
        if category in PHENO_CATEGORIES_MT:
            mt = hl.read_matrix_table(f'{DATA_PATH}/phenotype/raw/{category}_annotated.mt')
            mt = mt.filter_rows(mt.high_quality)
            pheno_mt = mt.filter_cols(mt.phenoname == phenoname)
            ht = pheno_mt.entries()
            print(f'-------{category} {ancestry} {phenoname} {ht.count()}')
        elif category == 'random_pheno':
            pheno_ht = hl.read_table(f'{DATA_PATH}/phenotype/random_pheno_mega_annotated.ht')
            ht = pheno_ht.filter(pheno_ht.high_quality)
            ht = ht.annotate(value=ht[phenoname])
        elif category in PHENO_CATEGORIES_HT:
            ht = hl.read_table(f"{DATA_PATH}/phenotype/{category}/lab_measurement_MEGA_annotated.ht")
            ht = ht.filter((ht.high_quality) & (ht.measurement_concept_id == phenoname))
            ht = ht.annotate(value=ht['median_value'])
        ht = ht.checkpoint(f'gs://aou_tmp/export_pheno/{ancestry}_{phenoname}_{timestamp}.ht', overwrite = True)
        print(f"Exporting {phenoname} for {ancestry.upper()} (n_samples: {ht.count()})...")
        return ht

    if ancestry.lower() != 'all':
        ht = load_pheno_ht(category=category, ancestry=ancestry, phenoname=phenoname)
    else:
        ht = load_mega_pheno_ht(category=category, phenoname=phenoname)

    if phenoname in binary_traits:
        ht = ht.annotate(value=hl.int(ht.value))
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
    print(f'Number of samples: {out_ht.count()}')
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
            phenos_to_run_by_ancestry = read_pickle_dict(f'{DATA_PATH}/phenotype/summary/pheno_by_ancestry_dict_{args.sex}.dict')
            phenos_to_run_by_ancestry_by_group = read_pickle_dict(f'{DATA_PATH}/phenotype/summary/pheno_by_ancestry_by_group_dict_{args.sex}.dict')
            if 'all' not in ancestries:
                print(f'----------Number of phenotypes per category per ancestry for {args.sex} (filtered to n_cases >= 200): --------------')
                print([f"{ancestry.upper()}-{category}: {len(phenos_to_run_by_ancestry_by_group[ancestry][category])}" for ancestry in ancestries for category in list(all_phenos_by_group.keys())])
                print([f"{ancestry.upper()}: {len(phenos_to_run_by_ancestry[ancestry])}" for ancestry in ancestries])
    else:
        phenos_to_run_by_ancestry = {}
        for ancestry in ancestries:
            phenos_to_run_by_ancestry[ancestry] = None
        phenos_to_run_by_ancestry_by_group = None

    if not args.dataproc:
        backend = hb.ServiceBackend(
            billing_project="all-by-aou", remote_tmpdir=TMP_BUCKET
        )

    if args.export_phenos:
        b = hb.Batch(
            name=f"aou_export_pheno_{args.ancestries}_{'test' if args.test or args.phenos or args.pilot is not None else 'all'}",
            backend=backend,
            default_storage="500Mi",
            default_cpu=n_threads,
            default_regions=['us-central1']
        )
        
        for ancestry in ancestries:
            print(f'Processing {ancestry.upper()}...')
            if args.gate_only:
                CATEGORIES = GATE_CATEGORIES
            else:
                CATEGORIES = PHENO_CATEGORIES
            for category in CATEGORIES:
                print(f'Processing {category}...')
                fields = basic_covariates
                if category not in list(phenos_to_run_by_ancestry_by_group[ancestry].keys()): continue
                if category not in GATE_CATEGORIES:
                    phenos_to_export = phenos_to_run_by_ancestry_by_group[ancestry][category]
                else:
                    phenos_to_export = []
                    for sex_group in ['both', 'male', 'female']:
                        phenos_to_run_by_ancestry_by_group = read_pickle_dict(f'{DATA_PATH}/phenotype/summary/pheno_by_ancestry_by_group_dict_{sex_group}.dict')
                        phenos_to_export = phenos_to_export + phenos_to_run_by_ancestry_by_group[ancestry][category]
                print(f"------------{ancestry.upper()} phenotype info: {len(phenos_to_export)} {category} phenotypes------------")
                for phenoname in tqdm(phenos_to_export):
                    if (not ancestry in phenoname) and (category in GATE_CATEGORIES):
                        continue
                    if category in GATE_CATEGORIES:
                        if not args.gate_only:
                            continue
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
        software_type = "gate_" if args.gate_only else ""
        variant_type = 'genome' if analysis_type == 'variant' else 'exome'
        sex_to_run = args.sex
        sex_tag = '' if sex_to_run == 'both' else f"_{sex_to_run}"
        print(f'Analysis type: {analysis_type}')
        print(f'Variant_type: {variant_type}')
        print(f"Docker image: {SAIGE_DOCKER_IMAGE}...")
        print(f"Sex to run: {sex_to_run}")
        RESULT_ROOT = f"{ANALYSIS_BUCKET}/{software_type}{analysis_type}_results{'' if not args.misc_test else '_misc'}"
        print(f"Result root: {RESULT_ROOT}")

        for ancestry in ancestries:
            if not args.skip_bgen or not args.skip_saige:
                size = CHUNK_SIZE[ancestry] if analysis_type == 'variant' else N_GENE_PER_GROUP[ancestry]
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
                    phenos_to_run = []
                    for sex_group in ['both', 'male', 'female']:
                        for category in GATE_CATEGORIES:
                            phenos_to_run_by_ancestry_by_group = read_pickle_dict(f'{DATA_PATH}/phenotype/summary/pheno_by_ancestry_by_group_dict_{sex_group}.dict')
                            phenos_to_run = phenos_to_run + phenos_to_run_by_ancestry_by_group[ancestry][category]
                else:
                    category = args.category
                    phenos_to_run = phenos_to_run_by_ancestry_by_group[ancestry][category]
                print(list(phenos_to_run)[0:5])
            else:
                category = 'full'
                phenos_to_run_by_ancestry['full'] = None
                phenos_to_run = phenos_to_run_by_ancestry[ancestry.lower()]

            b = hb.Batch(
                name=f"saige_{analysis_type}_aou_{ancestry}{sex_tag}",
                backend=backend,
                default_image=SAIGE_DOCKER_IMAGE,
                default_storage="500Mi",
                default_cpu=n_threads,
                default_regions=['us-central1']
            )

            # Obtain Sparse GRMs
            relatedness_cutoff = "0.125"
            num_markers = 2000
            n_threads = 8
            if ancestry.lower() != 'all':
                sparse_grm_root = f"{DATA_PATH}/utils/grm/aou_{ancestry}"
                sparse_grm_extension = f"_relatednessCutoff_{relatedness_cutoff}_{num_markers}_randomMarkersUsed.sparseGRM.mtx"
            else:
                sparse_grm_root = f"gs://aou_analysis/v8/fastGRM/ALL_ACAF"
                sparse_grm_extension = "sGRM.mtx"
            sparse_grm = b.read_input_group(
                **{ext: f"{sparse_grm_root}.{ext}"
                   for ext in (sparse_grm_extension, f"{sparse_grm_extension}.sampleIDs.txt",)}
            )

            overwrite_null_models = args.overwrite_null_models
            null_model_dir = f"{RESULT_ROOT}/null_glmm/{ancestry.upper()}"
            null_models_already_created = {}
            null_models = {}
            pheno_exports = {}
            if phenos_to_run is not None and not args.skip_any_null_models:
                print(phenos_to_run[0:5])
                print(f"------------{ancestry.upper()} {analysis_type} analysis null models: {len(phenos_to_run)} phenotypes------------")
                print(f'Null model directory: {null_model_dir}')
                if (not overwrite_null_models) and hfs.exists(null_model_dir):
                    null_models_already_created = {
                        x["path"] for x in hl.hadoop_ls(null_model_dir)
                    }
                print(f'Found {int(len(null_models_already_created) / 3)} Null models in directory...')
                for i in tqdm(range(len(phenos_to_run))):
                    phenoname = list(phenos_to_run)[i]
                    null_glmm_root = f"{null_model_dir}/phenotype_{phenoname}{sex_tag}"
                    model_file_path = f"{null_glmm_root}.rda"
                    variance_ratio_file_path =  f"{null_glmm_root}.varianceRatio.txt"
                    current_trait_type = None
                    if (phenoname in gate_traits) or ('_to_' in phenoname):
                        if not args.gate_only:
                            continue
                        if not ancestry in phenoname:
                            continue
                        current_trait_type = 'survival'
                        covariates = PC_fields + ['birth_year']
                        if not phenoname.startswith('Birth'):
                            covariates = covariates + ['age_at_first_event']
                        if phenoname.endswith('_ALL'):
                            covariates = covariates + ['sex']
                        covariates = ",".join(covariates)
                    elif phenoname in quantitative_traits or phenoname.endswith('_median'):
                        current_trait_type = 'quantitative'
                        covariates = basic_covariates
                    elif (phenoname in binary_traits) or (phenoname == 'C01AB'):
                        current_trait_type = 'binary'
                        covariates = basic_covariates
                    else:
                        raise ValueError(f"Unknown trait type for {phenoname}")

                    pheno_file = b.read_input(f'{ANALYSIS_BUCKET}/pheno_file{"" if not args.misc_test else "_misc"}/{ancestry.upper()}/phenotype_{phenoname}.tsv')
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
                        if category in GATE_CATEGORIES:
                            sex_arg = ''
                        else:
                            sex_arg = sex_tag
                        fit_null_task = fit_null_glmm(
                            b,
                            ancestry=ancestry,
                            output_root=null_glmm_root,
                            phenoname=f'{phenoname}{sex_tag}',
                            pheno_file=pheno_exports[phenoname],
                            sex_tag=sex_arg,
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
                            non_pre_emptible=True
                        )
                        fit_null_task.attributes.update(
                            {"phenotype": copy.deepcopy(phenoname),
                            "pheno_sex": copy.deepcopy(sex_to_run),
                            "analysis_type": copy.deepcopy(analysis_type),
                            "trait_type": copy.deepcopy(current_trait_type)}
                        )
                        model_file = fit_null_task.null_glmm.rda
                        variance_ratio_file = fit_null_task.null_glmm[f"varianceRatio.txt"]

                    null_models[f'{phenoname}{sex_tag}'] = (model_file, variance_ratio_file)
                else:
                    print('No phenotype loaded...')

            if not args.skip_bgen or not args.skip_saige:
                folder  = 'pgen' if args.pgen else 'bgen'
                print(f"------------{ancestry.upper()} {analysis_type} analysis {folder} files------------")
                geno_dir = f'{RESULT_ROOT.replace("gate_", "")}/{folder}/{ancestry.upper()}'
                bgen_dir = f'{RESULT_ROOT.replace("gate_", "")}/bgen/{ancestry.upper()}'
                pgen_dir = f'{RESULT_ROOT.replace("gate_", "")}/pgen/{ancestry.upper()}'
                print(f'{folder} directory: {geno_dir}')
                overwrite_bgens = args.overwrite_bgens
                overwrite_pgens = args.overwrite_pgens
                bgens_already_created = {}
                pgens_already_created = {}
                if not overwrite_bgens and hfs.exists(bgen_dir):
                    bgens_already_created = {x["path"] for x in hl.hadoop_ls(bgen_dir) if x["path"].endswith(".bgen")}
                if not overwrite_pgens and hfs.exists(pgen_dir) and args.pgen:
                    pgens_already_created = {x["path"] for x in hl.hadoop_ls(pgen_dir) if x["path"].endswith(".pgen")}
                print(f'Found {len(bgens_already_created)} bgens in directory...')
                print(f'Found {len(pgens_already_created)} pgens in directory...')

                genos = {}
                if args.test:
                    if ancestry.lower() == 'eur':
                        intervals = [hl.Interval(
                            hl.Locus('chr2', 169791933, reference_genome='GRCh38'),
                            hl.Locus('chr2', 190371665, reference_genome='GRCh38'), 
                            includes_start=True,
                            includes_end=True  
                            )] # TTN
                    print(intervals)

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
                    pgen_root = f"{pgen_dir}/{analysis_type}_{interval.start.contig}_{str(interval.start.position).zfill(9)}_{interval.end.position}"
                    bgen_task = None
                    if (f"{bgen_root}.bgen" not in bgens_already_created) or args.overwrite_bgens:
                    # if (f"{pgen_root}.pgen" not in pgens_already_created) or args.overwrite_bgens or (f"{bgen_root}.bgen" not in bgens_already_created):
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
                            bgen = True
                        )
                        b.write_output(bgen_task.log_file, f'gs://aou_tmp/000_bgen_logs/{analysis_type}_analysis_export_{str(interval)}_{folder}_{ancestry}.log')
                        bgen_task.attributes["ancestry"] = ancestry
                        bgen_task.attributes["analysis_type"] = analysis_type

                        bgen_index = index_bgen(b=b,
                                                ancestry=ancestry,
                                                analysis_type=analysis_type,
                                                bgen=f"{bgen_root}.bgen",
                                                depend_job=bgen_task)
                        bgen_index.attributes["ancestry"] = ancestry
                        bgen_index.attributes["analysis_type"] = analysis_type

                    if (not hl.hadoop_exists(f"{bgen_root}.gene.txt")  or args.overwrite_gene_txt)and analysis_type=='gene':
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

                    if ((f"{pgen_root}.pgen" not in pgens_already_created) or args.overwrite_pgens) and args.pgen:
                        pgen_task = from_bgen_to_pgen(
                            batch=b,
                            bgen_path=f"{bgen_root}.bgen",
                            sample_path=f"{bgen_root}.sample",
                            allele_order='ref-first',
                            depend_job=bgen_task,
                            highmem=args.highmem
                        )
                        pgen_task.attributes["ancestry"] = ancestry
                        pgen_task.attributes["analysis_type"] = analysis_type

                    if not args.pgen:
                        geno_file = b.read_input_group(
                            **{
                                "bgen": f"{bgen_root}.bgen",
                                "bgen.bgi": f"{bgen_root}.bgen.bgi",
                                "sample": f"{bgen_root}.sample",
                            }
                        )
                    else:
                        geno_file = b.read_input_group(
                            **{
                                "pgen": f"{pgen_root}.pgen",
                                "pvar": f"{pgen_root}.pvar",
                                "psam": f"{pgen_root}.psam",
                            }
                        )

                    if analysis_type == 'gene':
                        group_file = b.read_input(f"{bgen_root}.gene.txt")
                        genos[str(interval)] = (geno_file, group_file)
                    else:
                        genos[str(interval)] = geno_file

            elif not args.skip_saige:
                for interval in intervals:
                    if not args.pgen:
                        geno_root = f"{geno_dir}/{analysis_type}_{interval.start.contig}_{str(interval.start.position).zfill(9)}_{interval.end.position}"
                        geno_file = b.read_input_group(
                            **{
                                "bgen": f"{geno_root}.bgen",
                                "bgen.bgi": f"{geno_root}.bgen.bgi",
                                "sample": f"{geno_root}.sample",
                            }
                        )
                    else:
                        geno_root = f"{geno_dir}/{analysis_type}_{interval.start.contig}_{str(interval.start.position).zfill(9)}_{interval.end.position}"
                        geno_file = b.read_input_group(
                            **{
                                "pgen": f"{geno_root}.pgen",
                                "pvar": f"{geno_root}.pvar",
                                "psam": f"{geno_root}.psam",
                            }
                        )
                    if analysis_type == 'gene':
                        gene_group_root = geno_root.replace("pgen", "bgen")
                        group_file = b.read_input(f"{gene_group_root}.gene.txt")
                        genos[str(interval)] = (geno_file, group_file)
                    else:
                        genos[str(interval)] = geno_file

            result_dir = f'{RESULT_ROOT}/result/{ancestry.upper()}'
            print(f'result directory: {result_dir}')
            overwrite_results = args.overwrite_results
            saige_tasks = {}
            saige_ur_tasks = {}
            if not args.skip_saige:
                storage = None
                memory = '6.5G' if args.highmem else 'standard'
                if ancestry == 'all' and args.highmem:
                    storage = '20G'
                if args.highhighmem:
                    memory = '26G'
                    storage = '100G'
                
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
                    if (phenoname in gate_traits) or ('_to_' in phenoname):
                        if not args.gate_only:
                            continue
                        current_trait_type = 'survival'
                    elif phenoname in quantitative_traits:
                        current_trait_type = 'quantitative'
                    elif (phenoname in binary_traits) or (phenoname == 'C01AB'):
                        current_trait_type = 'binary'
                    else:
                        raise ValueError(f"Unknown trait type for {phenoname}")

                    saige_tasks[f'{phenoname}{sex_tag}'] = []
                    if analysis_type == 'gene':
                        saige_ur_tasks[f'{phenoname}{sex_tag}'] = []
                    if f'{phenoname}{sex_tag}' not in null_models.keys():
                        continue

                    model_file, variance_ratio_file = null_models[f'{phenoname}{sex_tag}']

                    if not i % 10:
                        n_jobs = dict(Counter(map(lambda x: x.name, b.select_jobs("")))).get("run_saige", 0)
                        logger.info(f"Read {i} phenotypes ({n_jobs} new to run so far)...")

                    pheno_results_dir = f"{result_dir}/phenotype_{phenoname}{sex_tag}"
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
                            f'-------------- Found {int(len(exome_results_already_created))} exome results in [{ancestry}: {phenoname}{sex_tag}] directory...')
                    factor = 2 if analysis_type == 'variant' else 5
                    print(f'-------------- Found {int(len(results_already_created) / factor)} {analysis_type} results in [{ancestry}: {phenoname}{sex_tag}] directory...')

                    for interval in intervals:
                        if (sex_tag != '') and (interval.start.contig != 'chrX') and (phenoname not in ['EM_202.2', 'EM_236.1', 'A10BJ', 'A10BJ06']):
                            continue
                        results_path = f"{pheno_results_dir}/result_{phenoname}_{interval.start.contig}_{str(interval.start.position).zfill(9)}"
                        group_file = None
                        max_maf_for_group_test = None
                        sparse_grm_file = None
                        suffix = ''
                        if analysis_type == "gene":
                            geno_file, group_file = genos[str(interval)]
                            max_maf_for_group_test = args.max_maf_group
                            sparse_grm_file = sparse_grm[sparse_grm_extension]
                        else:
                            geno_file = genos[str(interval)]
                        samples_file = b.read_input(f'{DATA_PATH}/utils/grm/{ancestry.upper()}_grm_plink.samples')

                        if (
                            overwrite_results
                            or (f"{results_path}.{'gene' if analysis_type == 'gene' else 'single_variant'}.txt"
                            not in results_already_created)
                        ):
                            saige_task = run_saige(
                                p=b,
                                phenoname=f'{phenoname}{sex_tag}',
                                ancestry=ancestry,
                                output_root=results_path,
                                model_file=model_file,
                                variance_ratio_file=variance_ratio_file,
                                sparse_grm_file=sparse_grm_file,
                                geno_file=geno_file,
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
                                memory=memory,
                                storage=storage,
                                pgen=args.pgen,
                            )
                            saige_task.attributes.update(
                                {"interval": str(interval), "ancestry": ancestry, 'name': f'saige{analysis_type == "gene"}'}
                            )
                            saige_task.attributes.update(
                                {"phenotype": copy.deepcopy(phenoname), "pheno_sex": copy.deepcopy(sex_to_run)}
                            )
                            saige_tasks[f'{phenoname}{sex_tag}'].append(saige_task)

                        if analysis_type == 'gene':
                            if (
                                    args.overwrite_exome_results
                                    or (f"{results_path}_exome.single_variant.txt" not in exome_results_already_created)
                            ):
                                memory = '6.5G' if args.highmem else 'standard'
                                if args.highhighmem:
                                    memory = '13G'
                                saige_ur_task = run_saige(
                                    p=b,
                                    phenoname=f'{phenoname}{sex_tag}',
                                    ancestry=ancestry,
                                    output_root=f'{results_path}_exome',
                                    model_file=model_file,
                                    variance_ratio_file=variance_ratio_file,
                                    sparse_grm_file=None,
                                    geno_file=geno_file,
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
                                    pgen=args.pgen,
                                )
                                saige_ur_task.attributes.update(
                                    {"interval": str(interval), "ancestry": ancestry}
                                )
                                saige_ur_task.attributes.update(
                                    {"phenotype": copy.deepcopy(phenoname), "pheno_sex": copy.deepcopy(sex_to_run)}
                                )
                                saige_ur_tasks[f'{phenoname}{sex_tag}'].append(saige_ur_task)
                    # if args.test:
                    #     break

            if not args.skip_load_hail_results:
                results_type = 'ACAF' if analysis_type == 'variant' else 'Exome'
                root = f"{ANALYSIS_BUCKET}/{software_type}ht_results{sex_tag}{'' if not args.misc_test else '_misc'}/{ancestry.upper()}"
                print(f'--------------Loading results from {results_type} analysis [{ancestry.upper()}]: {root}------------')
                for i in tqdm(range(len(phenos_to_run))):
                    phenoname = list(phenos_to_run)[i]
                    # if phenoname in gate_traits and not ancestry in phenoname:
                    #     continue
                    if phenoname in gate_traits:
                        if not args.gate_only:
                            continue
                        sex_tag = ''
                    directory = f'{ANALYSIS_BUCKET}/{software_type}{analysis_type}_results{"" if not args.misc_test else "_misc"}/result/{ancestry.upper()}/phenotype_{phenoname}{sex_tag}'
                    output_ht_directory = f'{root}/phenotype_{phenoname}{sex_tag}'
                    null_glmm_log = f'{ANALYSIS_BUCKET}/{software_type}{analysis_type}_results{"" if not args.misc_test else "_misc"}/null_glmm/{ancestry.upper()}/phenotype_{phenoname}{sex_tag}.log'

                    quantitative_trait = phenoname in quantitative_traits
   
                    if not args.skip_load_gene_results and analysis_type == 'gene':
                        if (not hfs.exists(f'{output_ht_directory}/gene_results_exact_expected_p.ht/_SUCCESS')) or args.overwrite_hail_results:
                            print(f'[{ancestry.upper()}: {phenoname}{sex_tag}] Gene table')
                            j = b.new_python_job(
                                name=f'sync_saige_{analysis_type}_gene_HT_{phenoname}{sex_tag}_{ancestry}',
                                attributes={"analysis_type": analysis_type, "ancestry": ancestry,
                                            "phenotype": copy.deepcopy(phenoname), "pheno_sex": copy.deepcopy(sex_to_run)})
                            j.image("hailgenetics/hail:0.2.133-py3.11")
                            j.memory('standard')
                            j.cpu(16)
                            j.env('PYSPARK_SUBMIT_ARGS', '--driver-memory 8g --executor-memory 8g pyspark-shell')
                            j.call(load_gene_data,
                                   directory=directory,
                                   output_ht_directory=output_ht_directory,
                                   phenoname=f'{phenoname}{sex_tag}',
                                   ancestry=ancestry,
                                   gene_ht_map_path=f'{DATA_PATH}/utils/gene_map/aou_{ancestry.upper()}_gene_map_{TRANCHE}.ht',
                                   quantitative_trait=quantitative_trait,
                                   null_glmm_log=null_glmm_log,
                                   overwrite=True,
                                   misc=args.misc_test
                                   )
                            if phenoname in list(saige_tasks.keys()) and not args.skip_saige:
                                j.depends_on(*saige_tasks[f'{phenoname}{sex_tag}'])
                            j.always_run()

                    if not args.skip_load_variant_results:
                        extension = 'single_variant.txt'
                        variant_type = 'genome' if analysis_type == 'variant' else 'exome'
                        SUCCESS = hfs.exists(f'{output_ht_directory}/{variant_type}_variant_results.ht/_SUCCESS') and hfs.exists(f'{output_ht_directory}/{variant_type}_variant_results_approx_cdf_expected_p.txt.bgz')
                        if args.overwrite_hail_results or not SUCCESS:
                            print(f'[{ancestry.upper()}: {phenoname}{sex_tag}] {variant_type} variant table')
                            j = b.new_python_job(
                                name=f'sync_saige_{variant_type}_variant_HT_{phenoname}{sex_tag}_{ancestry}',
                                attributes={"analysis_type": analysis_type, "ancestry": ancestry, "phenotype": copy.deepcopy(phenoname), "pheno_sex": copy.deepcopy(sex_to_run)})
                            j.image("hailgenetics/hail:0.2.133-py3.11")
                            j.memory('highmem')
                            j.cpu(16)
                            j.env('PYSPARK_SUBMIT_ARGS', '--driver-memory 16g --executor-memory 16g pyspark-shell')
                            j.call(load_variant_data,
                                    directory=directory,
                                    output_ht_directory=output_ht_directory,
                                    phenoname=f'{phenoname}{sex_tag}',
                                    ancestry=ancestry,
                                    quantitative_trait=quantitative_trait,
                                    null_glmm_log=null_glmm_log,
                                    extension=extension,
                                    overwrite=True,
                                    variant_type = variant_type,
                                    misc=args.misc_test
                                    )
                            j.always_run()
                        if variant_type == 'genome':
                            if phenoname in list(saige_tasks.keys()) and not args.skip_saige:
                                j.depends_on(*saige_tasks[f'{phenoname}{sex_tag}'])
                        else:
                            if phenoname in list(saige_ur_tasks.keys()) and not args.skip_saige:
                                j.depends_on(*saige_ur_tasks[f'{phenoname}{sex_tag}'])                      
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

    if args.create_meta_analysis_ht:
        # hl.stop()
        timestamp = datetime.now().strftime("%y%m%d_%H%M%S")
        if args.load_gene_results and args.single_variant_only:
            table_name = 'exome_variant_results'
        elif args.load_gene_results and not args.single_variant_only:
            table_name = 'gene_results'
        elif not args.load_gene_results and args.single_variant_only:
            table_name = 'genome_variant_results'

        # hl.init(
        #     tmp_dir=TMP_BUCKET,
        #     driver_memory="highmem",
        #     driver_cores=8,
        #     worker_memory="highmem",
        #     worker_cores=1,
        #     default_reference="GRCh38",
        #     log="/aou_saige.log",
        #     app_name=f'meta_analyze_SAIGE_{table_name}'
        # )
        
        analysis_type = "variant" if args.single_variant_only else "gene"
        result_type = "gene" if args.load_gene_results else "variant"   
        test_type = 'genome' if analysis_type == 'variant' and result_type == 'variant' else 'exome'
        sex_tag = '' if args.sex == 'both' else f'_{args.sex}'
        print(f'Analysis type: {analysis_type}')
        print(f'Test type: {test_type}')
        print(f'Sex tag: {sex_tag}')
        print(f'Result type: {result_type}')

        def run_meta_analysis_ht(pheno:str, table_name:str=table_name, analysis_type:str=analysis_type, overwrite:bool=True):
            # hl.experimental.init(backend='batch', app_name=f'meta_analysis_ht_{pheno}_{table_name}_QoBiB')

            paths = hl.hadoop_ls(
                f'{ANALYSIS_BUCKET}/ht_results/*/phenotype_{pheno}'
            )
            paths = [path['path'] for path in paths if('ALL' not in path['path']) and ('META' not in path['path']) and (path['is_dir']) and (path['path'].endswith(f'{table_name}.ht'))]
            meta_ancestry = [path.split('/')[-3].upper() for path in paths]
            ht = hl.read_table(paths[0])
            ht = ht.annotate(ancestry = meta_ancestry[0])
            if analysis_type == 'variant':
                ht = ht.annotate(BETA = hl.float64(ht.BETA),
                                 SE = hl.float64(ht.SE))
            n_cases = hl.eval(ht.n_cases)
            n_controls = hl.eval(ht.n_controls)
            N = hl.eval(ht.N)
            meta_heritability = [hl.eval(ht.heritability)]
            meta_n_cases = [hl.eval(ht.n_cases)]
            meta_n_controls = [hl.eval(ht.n_controls)]
            meta_N = [hl.eval(ht.N)]
            ht = ht.select_globals()
            ht = ht.annotate(n_cases = n_cases,
                            n_controls = n_controls,
                            N = N
            )
            fields = list(ht.row_value)
            if 'MAC_control' in fields:
                ht = ht.annotate(MAC_control = hl.int32(ht.MAC_control))
            ht = ht.select(*fields)
            for path in paths[1:]:
                print(path)
                tmp_ht = hl.read_table(path)
                meta_heritability.append(hl.eval(tmp_ht.heritability))
                meta_n_cases.append(hl.eval(tmp_ht.n_cases))
                meta_n_controls.append(hl.eval(tmp_ht.n_controls))
                meta_N.append(hl.eval(tmp_ht.N))
                if analysis_type == 'variant':
                    tmp_ht = tmp_ht.annotate(BETA=hl.float64(tmp_ht.BETA),
                                             SE=hl.float64(tmp_ht.SE))
                if 'MAC_control' in fields:
                    tmp_ht = tmp_ht.annotate(MAC_control = hl.int32(tmp_ht.MAC_control))
                n_cases = hl.eval(tmp_ht.n_cases)
                n_controls = hl.eval(tmp_ht.n_controls)
                N = hl.eval(tmp_ht.N)
                tmp_ht = tmp_ht.select_globals()
                tmp_ht = tmp_ht.annotate(
                    ancestry = path.split('/')[-3].upper(),
                )
                tmp_ht = tmp_ht.annotate(n_cases = n_cases,
                                         n_controls = n_controls,
                                         N = N
                )
                tmp_ht = tmp_ht.select(*fields)
                ht = ht.union(tmp_ht)
            print(meta_ancestry)
            print(meta_n_cases)
            print(meta_n_controls)
            print(meta_heritability)
            print(meta_N)
            ht = ht.key_by()
            ht = ht.drop('phenoname')
            ht = ht.select_globals(
                meta_ancestry=hl.literal(meta_ancestry if len(meta_ancestry) > 0 else hl.empty_array(hl.tstr)), 
                phenoname=pheno, 
                meta_n_cases=hl.literal(meta_n_cases if len(meta_n_cases) > 0 else hl.empty_array(hl.tint)), 
                meta_n_controls=hl.literal(meta_n_controls if len(meta_n_controls) > 0 else hl.empty_array(hl.tint)), 
                meta_heritability=hl.literal(meta_heritability if len(meta_heritability) > 0 else hl.empty_array(hl.tfloat64)),
                meta_N=hl.literal(meta_N if len(meta_N) > 0 else hl.empty_array(hl.tint))
            )

            if analysis_type == 'variant':
                ht = ht.annotate(weight = 1 / (ht.SE ** 2))
                meta_ht = ht.group_by('locus', 'alleles').aggregate(
                    BETA=hl.agg.sum(ht.weight * ht.BETA) / hl.agg.sum(ht.weight),
                    SE=1 / hl.sqrt(hl.agg.sum(ht.weight))
                )
                meta_ht.describe()
                ht = ht.annotate(META_BETA=meta_ht[ht.locus, ht.alleles].BETA)
                ht.describe()
                q_ht = ht.group_by('locus', 'alleles').aggregate(
                    Het_Q=hl.agg.sum(ht.weight * (ht.BETA - ht.META_BETA) ** 2),
                    AC_Allele2 = hl.agg.sum(ht.AC_Allele2),
                    AC_case = hl.agg.sum(ht.AF_case * ht.n_cases *2),
                    AC_ctrl = hl.agg.sum(ht.AF_ctrl * ht.n_controls *2)
                )
                N = hl.sum(q_ht.meta_N)
                q_ht = q_ht.annotate(
                    AF_Allele2 = q_ht.AC_Allele2 / (N*2)
                )
                q_ht.describe()
                meta_ht = meta_ht.annotate(**q_ht[meta_ht.key])

                meta_ht = meta_ht.annotate(Pvalue_log=hl.log(2) + hl.pnorm(-hl.abs(meta_ht.BETA/ meta_ht.SE), log_p=True))
                meta_ht = meta_ht.annotate(Pvalue=hl.exp(meta_ht.Pvalue_log))
                meta_ht = meta_ht.annotate(Pvalue_log10=-hl.log10(meta_ht.Pvalue))
            else:
                ht = ht.key_by('gene_id', 'gene_symbol', 'annotation')
                ht = ht.annotate(max_MAF=hl.if_else(hl.is_missing(ht.max_MAF), -1, ht.max_MAF))
                P_FIELDS = ['Pvalue', 'Pvalue_Burden', 'Pvalue_SKAT']
                P_TESTS = {'SKATO': 'Pvalue', 'Burden': 'Pvalue_Burden', 'SKAT': 'Pvalue_SKAT'}
                def _edit_pvalue(p):
                    return hl.if_else(p > 0.99, 0.99, p)

                ht = ht.filter(
                    hl.is_defined(ht[P_FIELDS[0]]) &
                    hl.is_defined(ht[P_FIELDS[1]]) &
                    hl.is_defined(ht[P_FIELDS[2]]))
                ht = ht.annotate(
                    **{p_field: _edit_pvalue(ht[p_field]) for p_field in P_FIELDS},
                    Neff = hl.if_else(ht.n_controls ==0, hl.float64(ht.n_cases), (4.0 * ht.n_cases * ht.n_controls) / (ht.n_cases + ht.n_controls))
                )
                ht = ht.annotate(CAF = ht.MAC/(2*ht.N))
                ht = ht.annotate(two_pq = 2*ht.CAF*(1 - ht.CAF))
                print(ht.aggregate(hl.agg.counter(ht.Neff)))
                print(ht.aggregate(hl.agg.counter(ht.n_cases)))
                print(ht.aggregate(hl.agg.counter(ht.n_controls)))
                ht = ht.filter(ht.annotation != 'Cauchy')

                def _stouffer_test(ht, test):
                    print(f'Meta analyzing {test} results...')
                    two_tail = test == 'Burden'
                    beta_field = f'BETA_{test}'
                    p_field = P_TESTS[test]

                    if two_tail:
                        ht = ht.annotate(**{p_field: ht[p_field] / 2},
                                                 **{beta_field: ht[beta_field]})
                    else:
                        ht = ht.annotate(
                            **{beta_field: hl.int(hl.is_defined(ht[p_field]))})
                    ht = ht.annotate(**{f'weighted_Z_numerator_{test}': (hl.sqrt(ht.two_pq)*hl.sqrt(ht.Neff) * (-hl.qnorm(ht[p_field])) * hl.sign(ht[beta_field]))})
                    meta_ht = ht.group_by('gene_id', 'gene_symbol', 'annotation', 'max_MAF').aggregate(
                        interval=hl.agg.collect(ht.interval),
                        MAC=hl.agg.sum(ht.MAC),
                        **{f'META_Stats_{test}': hl.agg.sum(ht[f'weighted_Z_numerator_{test}']) / (hl.sqrt(hl.agg.sum(ht.Neff*ht.two_pq)))},
                    )

                    if two_tail:
                        meta_ht = meta_ht.annotate(**{f'{p_field}': 2 * hl.pnorm(hl.abs(meta_ht[f'META_Stats_{test}']), lower_tail=False)})
                    else:
                        meta_ht = meta_ht.annotate(**{f'{p_field}': hl.pnorm(meta_ht[f'META_Stats_{test}'], lower_tail=False)})
                    meta_ht = meta_ht.annotate(**{f'{p_field}_log10': -hl.log10(meta_ht[p_field])})
                    return meta_ht

                ht.describe()
                meta_ht1 = _stouffer_test(ht, 'SKATO')
                meta_ht1 = meta_ht1.checkpoint(f'{TMP_BUCKET}/meta_tmp_{pheno}_SKATO_{test_type}_{result_type}_{timestamp}_12345.ht')
                meta_ht1.show(20)
                meta_ht2 = _stouffer_test(ht, 'SKAT').drop('interval', 'MAC')
                meta_ht2 = meta_ht2.checkpoint(
                    f'{TMP_BUCKET}/meta_tmp_{pheno}_SKAT_{test_type}_{result_type}_{timestamp}_12345.ht')
                meta_ht2.show(20)
                meta_ht3 = _stouffer_test(ht, 'Burden').drop('interval', 'MAC')
                meta_ht3 = meta_ht3.checkpoint(
                    f'{TMP_BUCKET}/meta_tmp_{pheno}_Burden_{test_type}_{result_type}_{timestamp}_12345.ht')
                meta_ht3.show(20)

                meta_ht = meta_ht1.annotate(
                    **meta_ht2[meta_ht1.key],
                    **meta_ht3[meta_ht1.key],
                )

            meta_ht = meta_ht.annotate_globals(
                N = hl.sum(meta_ht.meta_N),
                N_controls = hl.sum(meta_ht.meta_n_controls),
                N_cases = hl.sum(meta_ht.meta_n_cases)
            )
            if analysis_type == 'gene':
                meta_ht = meta_ht.annotate(CHR=meta_ht.interval[0].start.contig, POS=meta_ht.interval[0].start.position)
            else:
                meta_ht = meta_ht.annotate(CHR=meta_ht.locus.contig, POS=meta_ht.locus.position)
            meta_ht = meta_ht.annotate_globals(phenoname=pheno)
            meta_path = f'{ANALYSIS_BUCKET}/ht_results/META/phenotype_{pheno}/{table_name}.ht'
            meta_ht = meta_ht.checkpoint(meta_path, overwrite=overwrite)
            meta_ht.export(meta_path.replace('.ht', '.txt.bgz'))
            if analysis_type == 'gene':
                print(meta_ht.aggregate(hl.agg.counter(meta_ht.annotation)))
            meta_ht.describe()
            print(meta_ht.count())
            return meta_ht

        if not args.dataproc:
            backend = hb.ServiceBackend(
                billing_project="all-by-aou",
                remote_tmpdir='gs://aou_tmp/v8'
            )

            b = hb.Batch(
                name=f"aou_export_meta_ht_{table_name}",
                backend=backend,
                # default_storage="500Mi",
                # default_cpu=8,
                default_regions=['us-central1']
            )
        if args.pilot:
            phenotypes = PILOT_PHENOTYPES
        elif args.phenos is not None:
            phenotypes = args.phenos.split(',')
        else:
            phenos_to_run_by_pop = read_pickle_dict(f'{DATA_PATH}/phenotype/summary/pheno_by_ancestry_dict_{args.sex}.dict')
            phenotypes = phenos_to_run_by_pop['meta']
        for pheno in tqdm(phenotypes[args.meta_idx*361:min(args.meta_idx*361+361, len(phenotypes))]):
            print(pheno)
            j = None
            if not hl.hadoop_exists(f'{ANALYSIS_BUCKET}/ht_results/META/phenotype_{pheno}{sex_tag}/{table_name}.ht/_SUCCESS') or args.overwrite_meta:
                if not args.dataproc:
                    def new_query_in_batch_job(b: Batch, name: str, env: dict[str, str] | None = None) -> Job:
                        # creates a query-on-batch job in the current batch
                        backend = b._backend
                        assert isinstance(backend, ServiceBackend)

                        run_query_pipeline = textwrap.dedent(
                            f"""
                            hailctl config set batch/backend service
                            hailctl config set batch/regions {' '.join(backend.regions or [])}
                            hailctl config set batch/remote_tmpdir {backend.remote_tmpdir}
                            hailctl config set batch/billing_project {backend._billing_project}

                            cat << EOF | python3
                            from os import getenv
                            import hail as hl

                            batch_id = int(getenv('HAIL_BATCH_ID'))
                            hl.init(backend='batch', batch_id=batch_id, app_name='{name}')
                            hl.utils.range_table(2356)._force_count()
                            EOF
                            """,
                        )

                        j = b.new_bash_job(name=name)
                        j.command(run_query_pipeline)
                        for k, v in (env or dict()).items():
                            j.env(k, v)
                        return j

                    j = b.new_python_job(name=f"aou_{pheno}{sex_tag}_{table_name}")
                    j.image("hailgenetics/hail:0.2.133-py3.11")
                    if 'variant' in table_name:
                        j._machine_type = "n1-highmem-64"
                    else:
                        j.storage("30G")
                        j.cpu(16)
                        j.memory('100G')
                    j.spot(True)
                    j.env('HAIL_QUERY_BACKEND', 'batch')    
                    j.env('HAIL_BATCH_BILLING_PROJECT', 'all-by-aou')
                    j.env('HAIL_BATCH_REMOTE_TMPDIR', 'gs://aou_tmp/v8')
                    j.env('PYSPARK_SUBMIT_ARGS', '--driver-memory 100g --executor-memory 100g pyspark-shell')
                    j.call(run_meta_analysis_ht,
                        pheno=f'{pheno}{sex_tag}',
                        table_name=table_name,
                        analysis_type=analysis_type,
                        overwrite= True
                        )
                else:
                    run_meta_analysis_ht(
                        pheno=f'{pheno}{sex_tag}',
                        table_name=table_name,
                        analysis_type=analysis_type,
                        overwrite= True
                        )
            if ((analysis_type == 'variant') and (not hl.hadoop_exists(f"{ANALYSIS_BUCKET}/ht_results/META/phenotype_{pheno}/{table_name}_{'exact' if analysis_type == 'gene' else 'approx_cdf'}_expected_p.txt.bgz") or args.overwrite_meta)) or \
            ((analysis_type == 'gene') and (not hl.hadoop_exists(f"{ANALYSIS_BUCKET}/ht_results/META/phenotype_{pheno}/{table_name}_{'exact' if analysis_type == 'gene' else 'approx_cdf'}_expected_p.ht/_SUCCESS") or args.overwrite_meta)):
                if not args.dataproc:
                    j2 = b.new_python_job(name=f"aou_{pheno}_{table_name}_expected_p")
                    j2.image("hailgenetics/hail:0.2.133-py3.11")
                    j2.storage("20G")
                    j2.env('HAIL_QUERY_BACKEND', 'batch')    
                    j2.env('HAIL_BATCH_BILLING_PROJECT', 'all-by-aou')
                    j2.env('HAIL_BATCH_REMOTE_TMPDIR', 'gs://aou_tmp/v8')
                    j2.env('PYSPARK_SUBMIT_ARGS', '--driver-memory 100g --executor-memory 100g pyspark-shell')
                    if j is not None:
                        j2.depends_on(j)
                    j2.call(create_expected_p_ht,
                            pheno=f'{pheno}{sex_tag}',
                            ancestry='META',
                            table_name=table_name,
                            analysis_type=analysis_type,
                            method = 'exact' if analysis_type == 'gene' else 'approx_cdf',
                            overwrite = True
                            )
                else:
                    create_expected_p_ht(
                            pheno=f'{pheno}{sex_tag}',
                            ancestry='META',
                            table_name=table_name,
                            analysis_type=analysis_type,
                            method = 'exact' if analysis_type == 'gene' else 'approx_cdf',
                            overwrite = True
                            )
                if args.test:
                    break
        if not args.dataproc:
            b.run()

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
        "--overwrite-bgens", help="Force run of Bgens", action="store_true"
    )
    parser.add_argument(
        "--overwrite-gene-txt", help="Force run of gene txt files", action="store_true"
    )
    parser.add_argument(
        "--overwrite-results", help="Force run of SAIGE tests", action="store_true"
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
        "--load-gene-results", help="Load gene level test results", action="store_true"
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
    parser.add_argument(
        "--create-meta-analysis-ht", help="Create meta HT", action="store_true"
    )
    parser.add_argument(
        "--overwrite-meta", help="Overwrite meta HT", action="store_true"
    )
    parser.add_argument(
        "--highmem", help="Use high memory", action="store_true"
    )
    parser.add_argument(
        "--pgen", help="Export pgen instead of bgen", action="store_true"
    )
    parser.add_argument(
        "--overwrite-pgens", help="Overwrite pgen", action="store_true"
    )
    parser.add_argument(
        "--dataproc", help="Run on dataproc", action="store_true"
    )
    parser.add_argument(
        "--sex", help="Specify sex", default="both", choices=["male", "female", "both"]
    )
    parser.add_argument(
        "--misc-test", help="Run misc test", action="store_true"
    )
    parser.add_argument(
        "--highhighmem", help="Use high memory", action="store_true"
    )
    parser.add_argument(
        "--gate-only", help="Only run TTE phenotypes", action="store_true"
    )
    parser.add_argument(
        "--meta-idx", help="Index of meta analysis", type=int
    )
    args = parser.parse_args()

    main(args)
