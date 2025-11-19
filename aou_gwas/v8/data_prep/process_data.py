#!/usr/bin/env python3

__author__ = "wlu"

import sys
import os
import pandas as pd
import time
import hail as hl
import hailtop.fs as hfs
import hailtop.batch as hb
from hailtop.batch.docker import build_python_image
import pickle
import argparse
import hailtop.batch as hb
from gnomad.utils.vep import *
from gnomad.utils.filtering import *
from tqdm import tqdm

TRANCHE = "v8"
ANCESTRIES = ['AFR', 'AMR', 'EAS', 'EUR', 'MID', 'SAS']
N_GENE_PER_GROUP = 100
CHUNK_SIZE = {'all': int(1.25e6),'eur': int(6.25e6), "afr": int(1.25e7), "amr": int(1.25e7), "eas": int(1.25e7), "mid": int(1.25e7), "sas": int(1.25e7)}
REFERENCE = "GRCh38"
CHROMOSOMES = list(map(str, range(1, 23))) + ["X", "Y"]
ANALYSIS_BUCKET = "gs://aou_analysis/v8"
MY_BUCKET = 'gs://aou_wlu'
TMP_BUCKET = 'gs://aou_tmp/v8'
DATA_PATH = f'{ANALYSIS_BUCKET}/data'
ORIGINAL_GATE_PHENO_DIR = f'{DATA_PATH}/phenotype/gate'
ORIGINAL_GENO_ROOT = 'gs://fc-aou-datasets-controlled/v8' # https://batch.hail.is/batches/8251999 
ORIGINAL_GENO_DIR = f'{ORIGINAL_GENO_ROOT}/wgs/short_read/snpindel' # https://batch.hail.is/batches/8252511/jobs/1 
ORIGINAL_VDS_PATH = f'{ORIGINAL_GENO_DIR}/vds/hail.vds' # https://batch.hail.is/batches/8256092/jobs/1
ORIGINAL_AUX_DIR = f'{ORIGINAL_GENO_DIR}/aux' # https://batch.hail.is/batches/8252513/jobs/1 
ORIGINAL_QC_PATH = f'{ORIGINAL_AUX_DIR}/qc/genomic_metrics.tsv' 
ORIGINAL_ALL_SAMPLE_PATH = f'{ORIGINAL_AUX_DIR}/qc/all_samples.tsv'
ORIGINAL_FLAGGED_SAMPLE_PATH = f'{ORIGINAL_AUX_DIR}/qc/flagged_samples.tsv'
ORIGINAL_RELATED_PATH = f'{ORIGINAL_AUX_DIR}/relatedness/relatedness_flagged_samples.tsv'
ORIGINAL_KINSHIP_PATH = f'{ORIGINAL_AUX_DIR}/relatedness/relatedness.tsv'
ORIGINAL_GLOBAL_PCA_PATH = f'{ORIGINAL_AUX_DIR}/ancestry/ancestry_preds.tsv'
ORIGINAL_GLOBAL_LOADING_PATH = f'{ORIGINAL_AUX_DIR}/ancestry/loadings.ht'
ORIGINAL_EXOME_DIR = f'{ORIGINAL_GENO_DIR}/exome' # https://batch.hail.is/batches/8253890/jobs/1
ORIGINAL_ACAF_DIR = f'{ORIGINAL_GENO_DIR}/acaf_threshold' # https://batch.hail.is/batches/8253893/jobs/1

ORIGINAL_DATA_ROOT = 'gs://allxall-phenotypes-v2' # https://batch.hail.is/batches/8254254/jobs/1
ORIGINTAL_PHENO_DIR = f'{ORIGINAL_DATA_ROOT}/2025-05-14/phenotype_data' # https://batch.hail.is/batches/8257974/jobs/1 
# gs://allxall-phenotypes-v2/2025-05-14/phenotype_data/r_drug_table.parquet
# gs://allxall-phenotypes-v2/2025-05-14/phenotype_data/mcc2_phecode_table.parquet 
# gs://allxall-phenotypes-v2/2025-05-14/phenotype_data/mcc2_phecodex_table.parquet
# gs://allxall-phenotypes-v2/2025-05-14/phenotype_data/lab_measurement_table.parquet
# gs://allxall-phenotypes-v2/2025-05-14/phenotype_data/physical_measurement_table.parquet
# gs://allxall-phenotypes-v2/2025-05-14/phenotype_data/pfhh_survey_table.parquet

# gs://allxall-phenotypes-v2/2025-05-14/phenotype_data/ancestry_preds.parquet
# gs://allxall-phenotypes-v2/2025-05-14/phenotype_data/participant_level_age_covariates.parquet 
# gs://allxall-phenotypes-v2/2025-05-14/phenotype_data/demographic_table.parquet 

# gs://allxall-phenotypes-v2/2025-05-14/phenotype_data/drug_ages.parquet
# gs://allxall-phenotypes-v2/2025-05-14/phenotype_data/drug_cases_long.parquet
# gs://allxall-phenotypes-v2/2025-05-14/phenotype_data/phecode_NAs_long.parquet
# gs://allxall-phenotypes-v2/2025-05-14/phenotype_data/phecodex_NAs_long.parquet
# gs://allxall-phenotypes-v2/2025-05-14/phenotype_data/phecode_cases_long.parquet 
# gs://allxall-phenotypes-v2/2025-05-14/phenotype_data/phecodex_cases_long.parquet 
# gs://allxall-phenotypes-v2/2025-05-14/phenotype_data/drugs_controls_long1.parquet # 1-20
# gs://allxall-phenotypes-v2/2025-05-14/phenotype_data/phecode_controls_long1.parquet # 1-9
# gs://allxall-phenotypes-v2/2025-05-14/phenotype_data/phecodex_controls_long1.parquet # 1-17

ORIGINAL_GATE_PHENO_DIR = f'{ORIGINAL_DATA_ROOT}/2025-04-16/additional_data' # https://batch.hail.is/batches/8254263/jobs/1
ORIGINAL_DEMO_PATH = f'{ORIGINAL_DATA_ROOT}/2025-04-16/demographics_table.ht'
ORIGINAL_VAT_PATH = f'{ORIGINAL_DATA_ROOT}/2025-04-17/vat'
ORIGINAL_PCA_PATH = f'{ORIGINAL_DATA_ROOT}/2025-04-18/pca' # https://batch.hail.is/batches/8251124/jobs/1
ORIGINAL_CALLSTATS_PATH = f'{ORIGINAL_DATA_ROOT}/2025-04-16/call_stats'
ORIGINAL_GRM_PATH = f'{ORIGINAL_DATA_ROOT}/2025-04-16/grm' # https://batch.hail.is/batches/8254285/jobs/1 
ORIGINAL_RANDOM_PATH = f'{ORIGINAL_DATA_ROOT}/2025-04-16/random_pheno' # https://batch.hail.is/batches/8252587/jobs/1
ORIGINAL_HARD_FILTER = f'{ORIGINAL_DATA_ROOT}/2025-04-16/hard_filtered_samples.ht'

ORIGINAL_SOME_PATH ='gs://fc-secure-a960b16c-83a7-4026-b29d-72617355519a'

EXOME_MT_PATH = f'gs://fc-aou-datasets-controlled/v8/wgs/short_read/snpindel/exome/splitMT/hail.mt'
ACAF_MT_PATH =  f'gs://fc-aou-datasets-controlled/v8/wgs/short_read/snpindel/acaf_threshold/splitMT/hail.mt'
PHENO_CATEGORIES = ['physical_measurement', 'r_drug', 'pfhh_survey', 'random_pheno', 'onset', 'progression', 'lab_measurement', 'mcc2_phecode', 'mcc2_phecodex']
GATE_CATEGORIES = ['onset', 'progression']
BINARY_CATEGORIES = ['r_drug', 'pfhh_survey', 'mcc2_phecode', 'mcc2_phecodex']
QUANTITATIVE_CATEGORIES = ['physical_measurement', 'lab_measurement']

FILENAMES = PHENO_CATEGORIES + ['ancestry_preds', 'participant_level_age_covariates', 'demographic_table', 
 'drug_ages', 'drug_cases_long', 'phecode_NAs_long', 'phecodex_NAs_long', 
 'phecode_cases_long', 'phecodex_cases_long'] + \
    [f'drugs_controls_long{i}' for i in range(1, 21)] + \
    [f'phecode_controls_long{i}' for i in range(1, 10)] + \
    [f'phecodex_controls_long{i}' for i in range(1, 18)]

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
        overwrite_pruned_ht = True
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

def group_gene_interval(gtf, n, path, overwrite=False):
    from collections import Counter
    gtf = gtf.filter((gtf.feature == 'gene') & (gtf.gene_type == 'protein_coding'))
    gtf = gtf.annotate(
        chrom=gtf.interval.start.contig,
        start_pos=gtf.interval.start.position,
        end_pos=gtf.interval.end.position,
    )
    gtf = gtf.annotate(length=gtf.end_pos - gtf.start_pos + 1)
    gtf = gtf.add_index()

    n_gene = gtf.group_by("chrom").aggregate(n_gene_per_chrom=hl.agg.count())
    n_gene = n_gene.annotate(n_gene_per_chrom_round=
                             hl.if_else(n_gene.n_gene_per_chrom % n > 0,
                                        n_gene.n_gene_per_chrom + n - n_gene.n_gene_per_chrom % n,
                                        n_gene.n_gene_per_chrom
                                        ))
    n_gene = n_gene.annotate(chrom_index=n_gene.chrom[3:])
    n_gene = n_gene.annotate(
        chrom_index=hl.case()
        .when(n_gene.chrom_index == "X", "23")
        .when(n_gene.chrom_index == "Y", "24")
        .when(n_gene.chrom_index == "M", "25")
        .default(n_gene.chrom_index)
    )
    n_gene = n_gene.annotate(chrom_index=hl.int64(n_gene.chrom_index))
    n_gene = n_gene.order_by("chrom_index")

    n_gene_lst = n_gene.aggregate(
        hl.cumulative_sum(hl.agg.collect(n_gene.n_gene_per_chrom))
    )
    n_gene_lst.insert(0, 0)
    n_gene_lst.pop()
    print('------- Cumulative number of genes per chromosome: ---------')
    print(n_gene_lst)

    n_gene_round_lst = n_gene.aggregate(hl.cumulative_sum(hl.agg.collect(n_gene.n_gene_per_chrom_round)))
    n_gene_round_lst.insert(0, 0)
    n_gene_round_lst.pop()
    print('------- Cumulative number of genes per chromosome (rounded): ---------')
    print(n_gene_round_lst)


    chrom_lst = [f"chr{i + 1}" for i in range(22)] + ["chrX", "chrY", "chrM"]
    n_gene_df = pd.DataFrame(data={'chrom': chrom_lst, 'cum_n_gene': n_gene_lst, 'cum_n_gene_round':n_gene_round_lst})
    n_gene_ht = hl.Table.from_pandas(n_gene_df, key="chrom")
    print(n_gene_ht.show(26))

    gtf = gtf.annotate(n_previous_genes = n_gene_ht[gtf.chrom].cum_n_gene,
                       n_previous_genes_round = n_gene_ht[gtf.chrom].cum_n_gene_round,)
    gtf = gtf.annotate(
        new_idx=hl.if_else(
            gtf.n_previous_genes > 0, gtf.idx + gtf.n_previous_genes_round - gtf.n_previous_genes, gtf.idx
        )
    )

    gtf = gtf.annotate(group_id=gtf.new_idx // n)
    group_cnt = gtf.group_by('group_id').aggregate(cnt=hl.agg.count(),
                                                   group_length=hl.agg.sum(gtf.length))
    group_dict = gtf.aggregate(hl.agg.counter(gtf.group_id))
    print(Counter(group_dict.values()))

    gtf = gtf.annotate(n_genes_per_group=group_cnt[gtf.group_id].cnt)
    gtf = gtf.annotate(group_id=hl.if_else(gtf.n_genes_per_group < (n/2), gtf.group_id - 1, gtf.group_id))
    group_dict = gtf.aggregate(hl.agg.counter(gtf.group_id))
    print(Counter(group_dict.values()))

    gtf = gtf.select(
        "gene_name",
        "gene_id",
        "transcript_id",
        "chrom",
        "start_pos",
        "end_pos",
        "length",
        "idx",
        "new_idx",
        "group_id",
        "n_previous_genes",
    )

    sub_gtf = gtf.checkpoint(
        f"{path[:-3]}_tmp.ht",
        _read_if_exists=not overwrite,
        overwrite=overwrite,
    )

    # sub_gtf.show()

    group_ht = sub_gtf.group_by("group_id", "chrom").aggregate(
        group_start=hl.agg.min(sub_gtf.start_pos),
        group_end=hl.agg.max(sub_gtf.end_pos),
        genes=hl.agg.collect(sub_gtf.gene_name),
    )
    group_ht = group_ht.annotate(
        interval=hl.locus_interval(
            group_ht.chrom,
            group_ht.group_start,
            group_ht.group_end,
            reference_genome="GRCh38",
            includes_end=True,
            invalid_missing=True,
        )
    )

    group_ht = group_ht.checkpoint(
        path,
        _read_if_exists=not overwrite,
        overwrite=overwrite,
    )

    return group_ht

def update_variant_interval(chunk_size, overwrite=False):
    variant_interval_path = f'{DATA_PATH}/utils/intervals/aou_variant_interval_size_{chunk_size}.ht'
    if not hfs.exists(variant_interval_path) or overwrite:
        print(f'Generating interval file for SAIGE (chunk size: {chunk_size})...')
        intervals = []
        for chrom in CHROMOSOMES:
            chromosome = f"chr{chrom}"
            CHROMOSOME_LENGTHs = hl.get_reference(REFERENCE).lengths
            chrom_length = CHROMOSOME_LENGTHs[chromosome]
            for start_pos in range(1, chrom_length, chunk_size):
                end_pos = (
                    chrom_length
                if start_pos + chunk_size > chrom_length
                else (start_pos + chunk_size)
                )
                interval = hl.Interval(
                    hl.Locus(chromosome, start_pos, reference_genome='GRCh38'),
                    hl.Locus(chromosome, end_pos, reference_genome='GRCh38')
                )

                intervals.append(interval)
        import pandas as pd
        df = pd.DataFrame({'interval': intervals})
        ht = hl.Table.from_pandas(df)
        ht.write(variant_interval_path, overwrite=overwrite)
    ht = hl.read_table(variant_interval_path)
    ht.describe()
    ht.show(10)
    print(ht.count())

def update_gene_interval(n_genes_per_group, overwrite=False):
    gene_interval_path = f'{DATA_PATH}/utils/intervals/aou_gene_interval_size_{n_genes_per_group}.ht'
    if not hfs.exists(gene_interval_path) or overwrite:
        print(f'Generating interval file for SAIGE-GENE (N genes per group: {n_genes_per_group})...')
        gtf = hl.experimental.import_gtf('gs://hail-common/references/gencode/gencode.v29.annotation.gtf.bgz',
                                            reference_genome='GRCh38',
                                            skip_invalid_contigs=True)
        ht = group_gene_interval(gtf, n=n_genes_per_group, path=gene_interval_path, overwrite=overwrite)
    ht = hl.read_table(gene_interval_path)
    ht.describe()
    ht = ht.annotate(n_genes = hl.len(ht.genes))
    ht.filter(ht.n_genes > N_GENE_PER_GROUP).show(100)
    ht.filter(ht.n_genes < N_GENE_PER_GROUP/2).show(100)
    print(ht.aggregate(hl.agg.counter(ht.n_genes)))
    print(ht.count())

def write_pickle_dict(output: str, dict: dict):
    with hfs.open(output, "wb") as f:
        pickle.dump(dict, f)
    f.close()

def read_pickle_dict(dict_path: str):
    dict = pd.read_pickle(dict_path)
    return dict

def check_path(path = ORIGINAL_DATA_ROOT):
    hl.init(
        gcs_requester_pays_configuration='aou-neale-gwas'
    )
    info = hl.hadoop_ls(path)
    paths = [i['path'] for i in info]
    for path in paths:
        print(f'------{path}------')
        tmp_info = hl.hadoop_ls(path)
        tmp_paths = [i['path'] for i in tmp_info]
        for p in tmp_paths:
            print(p)

def copy_files(original_path, target_path):
    hl.hadoop_copy(original_path, target_path)
 
def write_raw_phenotype_ht(name, overwrite=False):
    # https://batch.hail.is/batches/8259477 
    back = hl.current_backend()
    spark = back._spark_session
    if name in PHENO_CATEGORIES:
        filename = f'{name}_table'
    else:
        filename = name
    parquetFile = spark.read.parquet(f"{ORIGINTAL_PHENO_DIR}/{filename}.parquet")
    ht = hl.Table.from_spark(parquetFile)
    ht.describe()
    if name not in PHENO_CATEGORIES and name not in ['ancestry_preds', 'participant_level_age_covariates', 'demographic_table']:
        outname = f"{name.replace('drugs', 'drug').split('_')[0]}_age/{name}"
    else:
        outname = name
    ht = ht.checkpoint(f'{DATA_PATH}/phenotype/raw/{outname}.ht', overwrite=overwrite, _read_if_exists= not overwrite)
    # ht.show()
    print(f"Number of samples: {ht.count()}")
    print(f"Number of phenotypes: {len(ht.row_value)}")
    return ht

def write_raw_gate_summary_ht(overwrite=False):
    ht = hl.import_table(f'{ORIGINAL_GATE_PHENO_DIR}/lists/*_summary_table.tsv', delimiter='\t', impute=True)
    ht = ht.checkpoint(f'{DATA_PATH}/phenotype/summary/gate_summary.ht', overwrite=overwrite, _read_if_exists= not overwrite)
    ht.describe()
    ht.show()
    print(ht.count())

    ht = hl.import_table(f'{ORIGINAL_GATE_PHENO_DIR}/lists/onset_phenotypes.tsv', delimiter='\t', no_header=True)
    ht = ht.rename({'f0' : 'phecode', 'f1': 'phenotype'})
    ht = ht.checkpoint(f'{DATA_PATH}/phenotype/summary/onset_overall_summary.ht', overwrite=overwrite, _read_if_exists= not overwrite)
    ht.describe()
    ht.show()
    print(ht.count())

    ht = hl.import_table(f'{ORIGINAL_GATE_PHENO_DIR}/lists/progression_phenotypes.tsv', delimiter='\t', no_header=True)
    ht = ht.rename({'f0' : 'phecode1', 'f1': 'phenotype1', 'f2': 'phecode2', 'f3': 'phenotype2'})
    ht = ht.checkpoint(f'{DATA_PATH}/phenotype/summary/progression_overall_summary.ht', overwrite=overwrite, _read_if_exists= not overwrite)
    ht.describe()
    ht.show()
    print(ht.count())

    ht = hl.import_table(f'{ORIGINAL_GATE_PHENO_DIR}/onset/*summary.txt', delimiter='\t', no_header=True)
    ht = ht.rename({'f0' : 'phecode', 'f1': 'phenotype'})
    ht = ht.checkpoint(f'{DATA_PATH}/phenotype/summary/onset_summary.ht', overwrite=overwrite, _read_if_exists= not overwrite)
    ht.describe()
    ht.show()
    print(ht.count())

    ht = hl.import_table(f'{ORIGINAL_GATE_PHENO_DIR}/progression/*summary.txt', delimiter='\t', no_header=True)
    ht.describe()
    ht = ht.rename({'f0' : 'phecode1', 'f1': 'phenotype1', 'f2': 'phecode2', 'f3': 'phenotype2'})
    ht = ht.checkpoint(f'{DATA_PATH}/phenotype/summary/progression_summary.ht', overwrite=overwrite, _read_if_exists= not overwrite)
    ht.describe()
    ht.show()
    print(ht.count())

def write_gate_pheno_ht(b, name, threshold=200, overwrite=False):
    """
    Reads all TSV files from a specified GCS directory corresponding to a GATE category (e.g., onset, progression),
    imports each into a Hail Table, and checkpoints it.

    Args:
        name (str): The GATE category name (e.g., 'onset', 'progression'). This corresponds to the subdirectory.
        overwrite (bool): Whether to overwrite existing Hail Tables.
    """
    # Files copied here: https://batch.hail.is/batches/8254727

    def process_gate_pheno(name, phenoname, overwrite=False):
        # https://batch.hail.is/batches/8254736
        raw_pheno_root = f'{DATA_PATH}/phenotype/raw/{name}'
        raw_pheno_path = f'{raw_pheno_root}/{phenoname}.ht'
        if not hl.hadoop_exists(raw_pheno_path):      
            ht = hl.import_table(
                raw_pheno_path.replace('.ht', '.tsv'),
                key='s', 
                types={'s': hl.tstr}, 
                impute=True, 
                delimiter='\t'
                )
            ht.describe()
            ht = ht.annotate(phenoname=phenoname, sex_from_data=ht.sex)
            if name == 'onset':
                ht = ht.select('phenoname', 'sex_from_data','birth_year', 'secondEvent', 'TTE')
            else:
                ht = ht.select('phenoname', 'sex_from_data','birth_year', 'secondEvent', 'TTE', 'age_at_first_event')
            ht = ht.checkpoint(raw_pheno_path, overwrite=overwrite, _read_if_exists=not overwrite)
            print(f"Count for {phenoname}: {ht.count()}")
        meta_ht = hl.read_table(f'{DATA_PATH}/utils/aou_v8_sample_meta.ht')
        print(meta_ht.aggregate(hl.agg.counter(meta_ht.ancestry)))
        ht = hl.read_table(raw_pheno_path)
        ht = ht.join(meta_ht, 'inner')
        ht = ht.checkpoint(f'{DATA_PATH}/phenotype/{name}/{phenoname}.ht', overwrite=overwrite, _read_if_exists=not overwrite)
        print(f"Count for {phenoname}: {ht.count()}")
        ht.describe()

    def summarize_processed_phenotypes(name, phenos_to_process, overwrite=False):
        n_cases = {}
        n_controls = {}
        for phenoname in phenos_to_process:
            print(phenoname)
            ht = hl.read_table(f'{DATA_PATH}/phenotype/{name}/{phenoname}.ht')
            n_cases[phenoname] = ht.aggregate(hl.agg.count_where(ht.secondEvent == 1))
            n_controls[phenoname] = ht.aggregate(hl.agg.count_where(ht.secondEvent == 0))
        rows = [{"phenoname": k, "n_cases": v, "n_controls": n_controls[k]} for k, v in n_cases.items()]
        ht = hl.Table.parallelize(rows, key="phenoname")
        ht = ht.checkpoint(f'{DATA_PATH}/phenotype/summary/{name}_summary_processed.ht', overwrite=overwrite, _read_if_exists= not overwrite)
        ht.show()
        ht.describe()
        print(ht.count())
        ht.export(f'{DATA_PATH}/phenotype/summary/{name}_summary_processed.txt.bgz')
            
    phenos_to_process = hl.import_table(f'{DATA_PATH}/phenotype/summary/INCLUDE_GWAS_{threshold}+_cases.tsv', no_header=True)
    phenos_to_process = list(phenos_to_process.aggregate(hl.agg.collect_as_set(phenos_to_process.f0)))
    if name == 'onset':
        phenos_to_process = [pheno for pheno in phenos_to_process if pheno.startswith('Birth')]
    elif name == 'progression':
        phenos_to_process = [pheno for pheno in phenos_to_process if not pheno.startswith('Birth')]
    print(f'Number of {name} phenotypes with N_cases >= {threshold}: {len(phenos_to_process)}')
    
    pheno_job = None
    for phenoname in phenos_to_process:
        if not hl.hadoop_exists(f'{DATA_PATH}/phenotype/{name}/{phenoname}.ht/_SUCCESS'):
            pheno_job = b.new_python_job(name=f"process_{name}_{phenoname}")
            pheno_job.call(process_gate_pheno, name, phenoname, overwrite=overwrite)

    if not hl.hadoop_exists(f'{DATA_PATH}/phenotype/summary/{name}_summary_processed.ht/_SUCCESS') or overwrite:
        summary_job = b.new_python_job(name=f"summarize_{name}_phenotypes")
        summary_job.memory("highmem")
        summary_job.cpu(8)
        summary_job.env('PYSPARK_SUBMIT_ARGS', '--driver-memory 32g --executor-memory 32g pyspark-shell')
        if pheno_job is not None:
            summary_job.depends_on(pheno_job)
        summary_job.call(summarize_processed_phenotypes, name, phenos_to_process, overwrite=overwrite)
    
def process_quantitative_phenotypes(overwrite):
    path = f'{DATA_PATH}/phenotype/physical_measurement.ht'
    if not hl.hadoop_exists(f'{path}/_SUCCESS') or overwrite:
        ht = hl.read_table(f'{DATA_PATH}/phenotype/raw/physical_measurement.ht')
        ## Annotate BMI and WHR
        ht = ht.annotate(BMI=ht.weight / ((ht.height / 100) ** 2),
                         WHR=ht['waist-circumference-mean'] / ht['hip-circumference-mean'])

        ## Annotate WHRadjBMI
        beta = ht.aggregate(hl.agg.linreg(ht.WHR, ht.BMI))['beta'][0]
        ht = ht.annotate(WHRadjBMI=ht.WHR - beta * ht.BMI)

        ht.write(path, overwrite=overwrite)

    ht = hl.read_table(path)
    ht.describe()
    print(ht.count())
    if not hl.hadoop_exists(path.replace(".ht", ".txt")):
        phenos = list(ht.row_value)
        with hl.hadoop_open(
                path.replace(".ht", ".txt"), "w"
        ) as f:
            for pheno in phenos:
                f.write(pheno + "\n")
    return ht

def summarize_raw_binary_phenotypes(name, overwrite):
    # https://batch.hail.is/batches/8259474
    mt = hl.read_matrix_table(f'{DATA_PATH}/phenotype/raw/{name}.mt')
    mt = mt.annotate_cols(n_cases = hl.agg.count_where(mt.value),
                      n_controls = hl.agg.count_where(~mt.value),
                      n_missing = hl.agg.count_where(hl.is_missing(mt.value)))
    if not hl.hadoop_exists(f'{DATA_PATH}/phenotype/summary/{name}_raw_summary.ht/_SUCCESS'):
        overwrite = True
    summary_ht = mt.cols().checkpoint(f'{DATA_PATH}/phenotype/summary/{name}_raw_summary.ht', overwrite=overwrite, _read_if_exists= not overwrite)
    if not hl.hadoop_exists(f'{DATA_PATH}/phenotype/summary/{name}_raw_summary.txt.bgz'):
        summary_ht.export(f'{DATA_PATH}/phenotype/summary/{name}_raw_summary.txt.bgz')
    return summary_ht

def summarize_lab_measurement(overwrite=False):
    ht = hl.read_table(f'{DATA_PATH}/phenotype/raw/lab_measurement.ht')
    if not hl.hadoop_exists(f'{DATA_PATH}/phenotype/summary/lab_measurement_raw_summary.ht/_SUCCESS'):
        n_samples = ht.aggregate(hl.agg.counter(ht.measurement_concept_id))
        max_min_value = ht.aggregate(hl.agg.group_by(ht.measurement_concept_id, hl.agg.max(ht.min_value)))
        min_min_value = ht.aggregate(hl.agg.group_by(ht.measurement_concept_id, hl.agg.min(ht.min_value)))
        max_median_value = ht.aggregate(hl.agg.group_by(ht.measurement_concept_id, hl.agg.max(ht.median_value)))
        min_median_value = ht.aggregate(hl.agg.group_by(ht.measurement_concept_id, hl.agg.min(ht.median_value)))
        max_max_value = ht.aggregate(hl.agg.group_by(ht.measurement_concept_id, hl.agg.max(ht.max_value)))
        min_max_value = ht.aggregate(hl.agg.group_by(ht.measurement_concept_id, hl.agg.min(ht.max_value)))
        max_latest_value = ht.aggregate(hl.agg.group_by(ht.measurement_concept_id, hl.agg.max(ht.latest_value)))
        min_latest_value = ht.aggregate(hl.agg.group_by(ht.measurement_concept_id, hl.agg.min(ht.latest_value)))
        max_count_value = ht.aggregate(hl.agg.group_by(ht.measurement_concept_id, hl.agg.max(ht['count'])))
        min_count_value = ht.aggregate(hl.agg.group_by(ht.measurement_concept_id, hl.agg.min(ht['count'])))

        # assume all share the same keys
        keys = sorted(n_samples.keys())

        # build a DataFrame
        df = pd.DataFrame({
            'phenoname': keys,
            'max_min_value': [max_min_value[k] for k in keys],
            'min_min_value': [min_min_value[k] for k in keys],
            'max_median_value': [max_median_value[k] for k in keys],
            'min_median_value': [min_median_value[k] for k in keys],
            'max_max_value': [max_max_value[k] for k in keys],
            'min_max_value': [min_max_value[k] for k in keys],
            'max_latest_value': [max_latest_value[k] for k in keys],
            'min_latest_value': [min_latest_value[k] for k in keys],
            'max_count_value': [max_count_value[k] for k in keys],
            'min_count_value': [min_count_value[k] for k in keys],
        })

        raw_summary_ht = hl.Table.from_pandas(df, key='phenoname')
        raw_summary_ht.describe()
        raw_summary_ht = raw_summary_ht.checkpoint(f'{DATA_PATH}/phenotype/summary/lab_measurement_raw_summary.ht', overwrite=True)
        raw_summary_ht.export(f'{DATA_PATH}/phenotype/summary/lab_measurement_raw_summary.txt.bgz')

    meta_ht = hl.read_table(f'{DATA_PATH}/utils/aou_v8_sample_meta.ht')
    ht = ht.key_by()
    ht = ht.annotate(person_id = hl.str(ht.person_id))
    ht = ht.key_by('person_id')
    ht = ht.annotate(ancestry = meta_ht[ht.key].ancestry, 
                     sex = meta_ht[ht.key].sex)
    qc_lab_ht = ht.filter(hl.is_defined(ht.ancestry))
    qc_summary_both_ht = qc_lab_ht.group_by('ancestry', 'measurement_concept_id').aggregate(n_cases = hl.agg.count())
    qc_summary_both_ht = qc_summary_both_ht.key_by()
    qc_summary_both_ht = qc_summary_both_ht.annotate(phenoname = qc_summary_both_ht.measurement_concept_id)
    qc_summary_both_ht = qc_summary_both_ht.key_by('phenoname')
    qc_summary_both_ht.describe()
    qc_summary_both_ht = qc_summary_both_ht.checkpoint(f'{DATA_PATH}/phenotype/summary/lab_measurement_both_summary.ht', overwrite=overwrite, _read_if_exists= not overwrite)
    
    qc_summary_male_ht = qc_lab_ht.filter(qc_lab_ht.sex == 1)
    qc_summary_male_ht = qc_summary_male_ht.group_by('ancestry', 'measurement_concept_id').aggregate(n_cases = hl.agg.count())
    qc_summary_male_ht = qc_summary_male_ht.key_by()
    qc_summary_male_ht = qc_summary_male_ht.annotate(phenoname = qc_summary_male_ht.measurement_concept_id + '_male')
    qc_summary_male_ht = qc_summary_male_ht.key_by('phenoname')
    qc_summary_male_ht.describe()
    qc_summary_male_ht = qc_summary_male_ht.checkpoint(f'{DATA_PATH}/phenotype/summary/lab_measurement_male_summary.ht', overwrite=overwrite, _read_if_exists= not overwrite)

    qc_summary_female_ht = qc_lab_ht.filter(qc_lab_ht.sex == 0)
    qc_summary_female_ht = qc_summary_female_ht.group_by('ancestry', 'measurement_concept_id').aggregate(n_cases = hl.agg.count())
    qc_summary_female_ht = qc_summary_female_ht.key_by()
    qc_summary_female_ht = qc_summary_female_ht.annotate(phenoname = qc_summary_female_ht.measurement_concept_id + '_female')
    qc_summary_female_ht = qc_summary_female_ht.key_by('phenoname')
    qc_summary_female_ht.describe()
    qc_summary_female_ht = qc_summary_female_ht.checkpoint(f'{DATA_PATH}/phenotype/summary/lab_measurement_female_summary.ht', overwrite=overwrite, _read_if_exists= not overwrite)
    
    qc_summary_both_ht = qc_summary_both_ht.key_by('phenoname')
    qc_summary_male_ht = qc_summary_male_ht.key_by('phenoname')
    qc_summary_female_ht = qc_summary_female_ht.key_by('phenoname')
    qc_summary_ht = qc_summary_both_ht.union(qc_summary_male_ht)
    qc_summary_ht = qc_summary_ht.union(qc_summary_female_ht)
    qc_summary_ht = qc_summary_ht.checkpoint(f'{DATA_PATH}/phenotype/summary/lab_measurement_summary.ht', overwrite=overwrite, _read_if_exists= not overwrite)
    qc_summary_ht.describe()
    qc_summary_ht.show(100)
    qc_summary_ht.export(f'{DATA_PATH}/phenotype/summary/lab_measurement_summary.txt.bgz')

def write_annotated_lab_ht(b, overwrite = False):
    if not hl.hadoop_exists(f'{DATA_PATH}/phenotype/lab_measurement/lab_measurement_full_annotated.ht/_SUCCESS'):
        raw_lab_ht = hl.read_table(f'{DATA_PATH}/phenotype/raw/lab_measurement.ht')
        meta_ht = hl.read_table(f'{DATA_PATH}/utils/aou_v8_sample_meta.ht')
        raw_lab_ht = raw_lab_ht.key_by()
        raw_lab_ht = raw_lab_ht.annotate(person_id = hl.str(raw_lab_ht.person_id))
        raw_lab_ht = raw_lab_ht.key_by('person_id')
        raw_lab_ht = raw_lab_ht.annotate(**meta_ht[raw_lab_ht.key])
        raw_lab_ht = raw_lab_ht.filter(hl.is_defined(raw_lab_ht.ancestry))
        annotated_lab_ht = raw_lab_ht.checkpoint(f'{DATA_PATH}/phenotype/lab_measurement/lab_measurement_full_annotated.ht', _read_if_exists= True)
        annotated_lab_ht.describe()
        print(annotated_lab_ht.count())
    
    annotated_lab_ht = hl.read_table(f'{DATA_PATH}/phenotype/lab_measurement/lab_measurement_full_annotated.ht')
    LAB_CODES = annotated_lab_ht.aggregate(hl.agg.collect_as_set(annotated_lab_ht.measurement_concept_id))
    print(f'Number of LAB CODES: {len(LAB_CODES)}')
    def _export_lab_ht(code, overwrite = False):
        hl.init(app_name=f'Export lab measurement {code}', 
                gcs_requester_pays_configuration='aou-neale-gwas', 
                master='local[32]', 
                tmp_dir=TMP_BUCKET, 
                # worker_memory="highmem", 
                # worker_cores=8, 
                default_reference="GRCh38", 
                log="/export_lab_measurement.log")
        annotated_lab_ht = hl.read_table(f'{DATA_PATH}/phenotype/lab_measurement/lab_measurement_full_annotated.ht')
        ht = annotated_lab_ht.filter(annotated_lab_ht.measurement_concept_id == code)
        if not hl.hadoop_exists(f'{DATA_PATH}/phenotype/lab_measurement/lab_measurement_{code}.ht/_SUCCESS'):
            overwrite=True
        print(f'Processing {code}')
        ht = ht.checkpoint(f'{DATA_PATH}/phenotype/lab_measurement/lab_measurement_{code}.ht', overwrite=overwrite, _read_if_exists= not overwrite)

    for code in LAB_CODES:
        j = b.new_python_job(name=f'Export lab measurement {code}')
        j.call(_export_lab_ht, code, overwrite)
        
def extract_mt_from_ht(name = 'r_drug', overwrite = False):
    if name != 'physical_measurement':
        ht = hl.read_table(f'{DATA_PATH}/phenotype/raw/{name}.ht')
    else: 
        ht = hl.read_table(f'{DATA_PATH}/phenotype/{name}.ht')
    if name in ['mcc2_phecode', 'mcc2_phecodex']:
        ht = ht.key_by(person_id  = hl.format("%.0f", ht.person_id) )
    else:
        ht = ht.key_by(person_id = hl.str(ht.person_id))
    columns = list(ht.row_value)
    ht2 = ht.select(values=[hl.struct(value=ht[x])for x in columns]).annotate_globals(
        columns=hl.map(lambda x: hl.struct(phenoname=x),
                       hl.literal(columns)))
    ht2 = ht2.checkpoint(f'{TMP_BUCKET}/v8/{name}.ht', overwrite=True)
    mt = ht2._unlocalize_entries('values', 'columns', [])
    mt = mt.checkpoint(f'{DATA_PATH}/phenotype/raw/{name}.mt', overwrite=overwrite, _read_if_exists= not overwrite)
    mt.describe()
    return mt

def process_vat(overwrite=False):
    # https://batch.hail.is/batches/8250503/jobs/1
    ht = hl.read_table(f'{ORIGINAL_VAT_PATH}/aou_PARSED_SORTED_COLLECTED_vat_wlu.ht')
    ht = ht.checkpoint(f'{DATA_PATH}/vat/aou_PARSED_SORTED_COLLECTED_vat_v8.ht', overwrite=overwrite, _read_if_exists= not overwrite)
    return ht

def process_vds(overwrite=False):
    hl.init(
        app_name=f'Process_vds_v8',
        gcs_requester_pays_configuration='aou-neale-gwas',
        master='local[32]',
        tmp_dir=TMP_BUCKET,
        worker_memory="highmem",
        worker_cores=8,
        default_reference="GRCh38",
        log="/process_vds.log"
    )
    vds = hl.vds.read_vds(ORIGINAL_VDS_PATH)
    variant_data = vds.variant_data
    variant_data_ht = variant_data.rows()
    variant_data_ht = variant_data_ht.checkpoint(f'{DATA_PATH}/vat/aou_vds_variant_data_row_v8.ht', overwrite=True)
    variant_data_ht.describe()
    print(f"Number of variants: {variant_data_ht.count()}")
    
def process_pca(overwrite=False):
    # https://batch.hail.is/batches/8256149/jobs/1
    hl.init(
        app_name=f'Process_pca_v8',
        gcs_requester_pays_configuration='aou-neale-gwas',
        master='local[32]',
        tmp_dir=TMP_BUCKET,
        worker_memory="highmem",
        worker_cores=8,
        default_reference="GRCh38",
        log="/process_pca.log"
    )
    if not hl.hadoop_exists(f'{DATA_PATH}/utils/pca/pca_pruned_scores_full.ht/_SUCCESS') or overwrite:
        ht = hl.import_table(f'gs://aou_analysis/v8/data/utils/pca/results/pca_{ANCESTRIES[0]}_pruned_scores_full.txt.bgz', key='s')
        ht = ht.annotate(ancestry = ANCESTRIES[0])
        for ancestry in ANCESTRIES[1:]:
            tmp_ht = hl.import_table(f'gs://aou_analysis/v8/data/utils/pca/results/pca_{ancestry}_pruned_scores_full.txt.bgz', key='s')   
            tmp_ht = tmp_ht.annotate(ancestry = ancestry)
            ht = ht.union(tmp_ht)
        ht = ht.checkpoint(f'{DATA_PATH}/utils/pca/pca_pruned_scores_full.ht', overwrite=overwrite, _read_if_exists= not overwrite)
    ht = hl.read_table(f'{DATA_PATH}/utils/pca/pca_pruned_scores_full.ht')
    ht.describe()
    print(f"Number of samples: {ht.count()}")
    print(f"Number of ancestries: {ht.aggregate(hl.agg.counter(ht.ancestry))}")

    ht = hl.import_table(f'{ORIGINAL_GLOBAL_PCA_PATH}', impute=True, key='research_id', types={"research_id":"tstr","pca_features":hl.tarray(hl.tfloat)})
    ht = ht.checkpoint(f'{DATA_PATH}/utils/aou_v8_global_pca.ht', overwrite=overwrite, _read_if_exists= not overwrite)
    ht.describe()
    print(f"Number of samples: {ht.count()}")

    ht = hl.read_table(f'{ORIGINAL_GLOBAL_LOADING_PATH}')
    ht = ht.checkpoint(f'{DATA_PATH}/utils/aou_v8_global_loadings.ht', overwrite=overwrite, _read_if_exists= not overwrite)
    ht.describe()
    print(f"Number of samples: {ht.count()}")
    return ht

def process_hard_filter_data(overwrite=False):
    # https://batch.hail.is/batches/8252542/jobs/1
    ht = hl.read_table(ORIGINAL_HARD_FILTER)
    ht = ht.checkpoint(f'{DATA_PATH}/utils/aou_v8_hard_filter.ht', overwrite=overwrite, _read_if_exists= not overwrite)
    ht.describe()
    print(f"Number of samples: {ht.count()}")
    return ht

def process_demographic_data(overwrite=False):
    # https://batch.hail.is/batches/8252550/jobs/1
    hl.init(gcs_requester_pays_configuration='aou-neale-gwas')
    ht = hl.import_table(f'{ORIGINAL_QC_PATH}', key='research_id')
    ht = ht.checkpoint(f'{DATA_PATH}/utils/aou_v8_genomic_qc_metric.ht', overwrite=overwrite, _read_if_exists= not overwrite)
    ht.describe()
    print(f"Number of samples: {ht.count()}")

    ht = hl.read_table(f'{ORIGINAL_DEMO_PATH}')
    ht = ht.checkpoint(f'{DATA_PATH}/utils/aou_v8_demographic.ht', overwrite=overwrite, _read_if_exists= not overwrite)
    ht.describe()
    print(f"Number of samples: {ht.count()}")

    ht = hl.import_table(f'{ORIGINAL_ALL_SAMPLE_PATH}', key='s')
    ht = ht.checkpoint(f'{DATA_PATH}/utils/aou_v8_all_samples.ht', overwrite=overwrite, _read_if_exists= not overwrite)
    ht.describe()
    print(f"Number of samples: {ht.count()}")

    ht = hl.import_table(f'{ORIGINAL_FLAGGED_SAMPLE_PATH}', key='s')
    ht = ht.checkpoint(f'{DATA_PATH}/utils/aou_v8_flagged_samples.ht', overwrite=overwrite, _read_if_exists= not overwrite)
    ht.describe()
    print(f"Number of samples: {ht.count()}")
    return ht

def process_related_data(overwrite=False):
    # https://batch.hail.is/batches/8252550/jobs/2
    hl.init(gcs_requester_pays_configuration='aou-neale-gwas')
    ht = hl.import_table(f'{ORIGINAL_RELATED_PATH}', key='sample_id')
    ht = ht.checkpoint(f'{DATA_PATH}/utils/aou_v8_related_samples.ht', overwrite=overwrite, _read_if_exists= not overwrite)
    ht.describe()
    print(f"Number of samples: {ht.count()}")

    ht = hl.import_table(f'{ORIGINAL_KINSHIP_PATH}', key=['i.s', 'j.s'], types = {'i.s':hl.tstr, 'j.s':hl.tstr}, impute = True)
    ht = ht.checkpoint(f'{DATA_PATH}/utils/aou_v8_related_kinship.ht', overwrite=overwrite, _read_if_exists= not overwrite)
    ht.describe()
    print(f"Number of samples: {ht.count()}")
    return ht

def update_meta_ht(overwrite: bool):
    # https://batch.hail.is/batches/8256150/jobs/1 
    hl.init(
        app_name=f'Process_meta_v8',
        gcs_requester_pays_configuration='aou-neale-gwas',
        master='local[32]',
        tmp_dir=TMP_BUCKET,
        worker_memory="highmem",
        worker_cores=8,
        default_reference="GRCh38",
        log="/process_meta.log"
    )
    pca_ht = hl.read_table(f'{DATA_PATH}/utils/pca/pca_pruned_scores_full.ht')
    pca_ht.describe()
    print(f'Number of pruned samples: {pca_ht.count()}')
    print(f'Number of pruned samples per group: {pca_ht.aggregate(hl.agg.counter(pca_ht.ancestry))}')
    if (
        not hl.hadoop_exists(f'{DATA_PATH}/utils/aou_v8_sample_meta.ht')
        or overwrite
    ):
        hard_filter_ht = hl.read_table(f'{DATA_PATH}/utils/aou_v8_hard_filter.ht')
        hard_filter_ht = hard_filter_ht.select(
            age = hard_filter_ht.age_at_cdr,
            sex = hard_filter_ht.sex_at_birth,
            race = hard_filter_ht.race,
            age_group = hard_filter_ht.age_group
        )
        related_ht = hl.read_table(f'{DATA_PATH}/utils/aou_v8_related_samples.ht')
        meta_ht = pca_ht.annotate(
            **hard_filter_ht[pca_ht.key],
            hard_filter=hl.is_defined(hard_filter_ht[pca_ht.key]),
            related=hl.is_defined(related_ht[pca_ht.key]),
        )
        meta_ht.describe()
        meta_ht.summarize()
        meta_ht = meta_ht.annotate(
            sex=hl.if_else(
                meta_ht.sex == "Male",
                1,
                hl.if_else(meta_ht.sex == "Female", 0, hl.missing(hl.tint32)),
            ),  # https://atgu.slack.com/archives/C056MJF70J1/p1695909765150529?thread_ts=1695326201.653099&cid=C056MJF70J1
            age=meta_ht.age,
            age2=meta_ht.age**2,
        )
        meta_ht = meta_ht.annotate(
            age_sex=meta_ht.age * meta_ht.sex, age2_sex=meta_ht.age2 * meta_ht.sex
        )

        #### SAVE FOR FUTURE USE ####
        # meta_ht = meta_ht.annotate(pop=pop_ht[meta_ht.key].ancestry_pred)
        # meta_ht = meta_ht.annotate(
        #     **{f"pop_PC{i + 1}": pca_ht[meta_ht.key][f'PC{i+1}'] for i in range(50)}
        # )
        #### END ####

        meta_ht = meta_ht.naive_coalesce(100).checkpoint(
            f'{DATA_PATH}/utils/aou_v8_sample_meta.ht', overwrite=overwrite, _read_if_exists= not overwrite
        )
    meta_ht = hl.read_table(f'{DATA_PATH}/utils/aou_v8_sample_meta.ht')
    meta_ht.summarize()
    print(meta_ht.aggregate(hl.agg.counter(meta_ht.ancestry)))
    print(meta_ht.count())
    return meta_ht

def annotate_phenotype_tables(name, overwrite=False):
    # https://batch.hail.is/batches/8257165
    hl.init(
        app_name=f'Process_phenotype_v8',
        gcs_requester_pays_configuration='aou-neale-gwas',
        master='local[32]',
        tmp_dir=TMP_BUCKET,
        worker_memory="highmem",
        worker_cores=8,
        default_reference="GRCh38",
        log="/process_phenotype.log"
    )
    meta_ht = hl.read_table(f'{DATA_PATH}/utils/aou_v8_sample_meta.ht')
    if name in ['random_pheno']:
        data = hl.read_table(f'{DATA_PATH}/phenotype/{name}.ht')
        data.describe()
        data = data.key_by()
        data = data.annotate(person_id = hl.str(data.person_id))
        data = data.key_by('person_id')
        data = data.annotate(**meta_ht[data.key])
        data = data.filter(hl.is_defined(data.ancestry))
        data = data.checkpoint(f'{DATA_PATH}/phenotype/{name}_annotated.ht', overwrite=overwrite, _read_if_exists= not overwrite)
    else:
        data = hl.read_matrix_table(f'{DATA_PATH}/phenotype/raw/{name}.mt')
        data.describe()
        data = data.key_rows_by()
        if name in ['mcc_phecode', 'mcc_phecodex']:
            data = data.annotate_rows(person_id = hl.format("%.0f", data.person_id))
        else:
            data = data.annotate_rows(person_id = hl.str(data.person_id))
        data = data.key_rows_by('person_id')
        data = data.key_cols_by('phenoname')
        data = data.annotate_rows(**meta_ht[data.row_key])
        data = data.filter_rows(hl.is_defined(data.ancestry))
        data = data.checkpoint(f'{DATA_PATH}/phenotype/{name}_annotated.mt', overwrite=overwrite, _read_if_exists= not overwrite)
    return data

def summarize_quantitative_phenotypes(name, overwrite):
    # https://batch.hail.is/batches/8257165
    hl.init(
        app_name=f'Summarize_quantitative_phenotypes_v8',
        gcs_requester_pays_configuration='aou-neale-gwas',
        master='local[32]',
        tmp_dir=TMP_BUCKET,
        worker_memory="highmem",
        worker_cores=8,
        default_reference="GRCh38",
        log="/summarize_quantitative_phenotypes.log"
    )
    mt = hl.read_matrix_table(f'{DATA_PATH}/phenotype/{name}_annotated.mt')
    def _compute_summary_table(mt):
        mt = mt.annotate_cols(n_cases = hl.agg.count_where(hl.is_defined(mt.value)),
                              n_cases_by_ancestry=hl.agg.group_by(mt.ancestry, hl.agg.count_where(hl.is_defined(mt.value))),
                        )
        mt = mt.annotate_cols(**{f'n_cases_{ancestry}': mt.n_cases_by_ancestry.get(ancestry) for ancestry in ANCESTRIES})
        summary_ht = mt.cols()
        return summary_ht  
    summary_both_ht = _compute_summary_table(mt)
    summary_both_ht = summary_both_ht.key_by('phenoname')
    summary_both_ht = summary_both_ht.checkpoint(f'{DATA_PATH}/phenotype/summary/{name}_both_summary.ht', overwrite=overwrite, _read_if_exists= not overwrite)
    summary_both_ht.describe()
    summary_male_ht = _compute_summary_table(mt.filter_rows(mt.sex == 1))
    summary_male_ht = summary_male_ht.key_by()
    summary_male_ht = summary_male_ht.key_by(phenoname = summary_male_ht.phenoname + '_male')
    summary_male_ht = summary_male_ht.checkpoint(f'{DATA_PATH}/phenotype/summary/{name}_male_summary.ht', overwrite=overwrite, _read_if_exists= not overwrite)
    summary_male_ht.describe()
    summary_female_ht = _compute_summary_table(mt.filter_rows(mt.sex == 0))
    summary_female_ht = summary_female_ht.key_by()
    summary_female_ht = summary_female_ht.key_by(phenoname = summary_female_ht.phenoname + '_female')
    summary_female_ht = summary_female_ht.checkpoint(f'{DATA_PATH}/phenotype/summary/{name}_female_summary.ht', overwrite=overwrite, _read_if_exists= not overwrite)
    summary_female_ht.describe()
    summary_ht = summary_both_ht.union(summary_male_ht).union(summary_female_ht)
    summary_ht = summary_ht.checkpoint(f'{DATA_PATH}/phenotype/summary/{name}_summary.ht', overwrite=overwrite, _read_if_exists= not overwrite)
    summary_ht.describe()
    summary_ht.show(100)
    summary_ht.export(f'{DATA_PATH}/phenotype/summary/{name}_summary.txt.bgz')
    return summary_ht

def summarize_binary_phenotypes(name, overwrite):
    # https://batch.hail.is/batches/8257165
    hl.init(
        app_name=f'Summarize_binary_phenotypes_v8',
        gcs_requester_pays_configuration='aou-neale-gwas',
        master='local[32]',
        tmp_dir=TMP_BUCKET,
        worker_memory="highmem",
        worker_cores=8,
        default_reference="GRCh38",
        log="/summarize_binary_phenotypes.log"
    )
    mt = hl.read_matrix_table(f'{DATA_PATH}/phenotype/{name}_annotated.mt')
    def _compute_summary_table(mt):
        mt = mt.select_cols(n_cases = hl.agg.count_where(mt.value),
                        n_controls = hl.agg.count_where(~mt.value),
                        n_missing = hl.agg.count_where(hl.is_missing(mt.value)),
                        n_cases_by_ancestry=hl.agg.group_by(mt.ancestry, hl.agg.count_where(mt.value)),
                        n_controls_by_ancestry=hl.agg.group_by(mt.ancestry, hl.agg.count_where(~mt.value)),
                        n_samples_defined_by_ancestry=hl.agg.group_by(mt.ancestry, hl.agg.count_where(hl.is_defined(mt.value))),
                        )
        mt = mt.annotate_cols(**{f'n_cases_{ancestry}': mt.n_cases_by_ancestry.get(ancestry) for ancestry in ANCESTRIES},
                              **{f'n_controls_{ancestry}': mt.n_controls_by_ancestry.get(ancestry) for ancestry in ANCESTRIES},
                              **{f'n_samples_defined_{ancestry}': mt.n_samples_defined_by_ancestry.get(ancestry) for ancestry in ANCESTRIES},
                              )

        summary_ht = mt.cols()
        summary_ht.describe()
        if name == 'r_drug':
            # https://batch.hail.is/batches/8253446/jobs/1
            for percent in [70, 80, 90]:
                ht = hl.import_table(f'{DATA_PATH}/phenotype/summary/r_drug_curation_{percent}.csv', key='atc_code', impute=True, types={"atc_code": hl.tstr}, delimiter=',')
                summary_ht = summary_ht.annotate(**{f'r_drug_curation_{percent}': ht[summary_ht.key].is_kept})
        return summary_ht 
    
    summary_both_ht = _compute_summary_table(mt)
    summary_both_ht = summary_both_ht.key_by('phenoname')
    summary_both_ht = summary_both_ht.checkpoint(f'{DATA_PATH}/phenotype/summary/{name}_both_summary.ht', overwrite=overwrite, _read_if_exists= not overwrite)
    summary_both_ht.describe()
    
    summary_male_ht = _compute_summary_table(mt.filter_rows(mt.sex == 1))
    summary_male_ht = summary_male_ht.key_by()
    summary_male_ht = summary_male_ht.key_by(phenoname = summary_male_ht.phenoname + '_male')
    summary_male_ht = summary_male_ht.checkpoint(f'{DATA_PATH}/phenotype/summary/{name}_male_summary.ht', overwrite=overwrite, _read_if_exists= not overwrite)
    summary_male_ht.describe()

    summary_female_ht = _compute_summary_table(mt.filter_rows(mt.sex == 0))
    summary_female_ht = summary_female_ht.key_by()
    summary_female_ht = summary_female_ht.key_by(phenoname = summary_female_ht.phenoname + '_female')
    summary_female_ht = summary_female_ht.checkpoint(f'{DATA_PATH}/phenotype/summary/{name}_female_summary.ht', overwrite=overwrite, _read_if_exists= not overwrite)
    summary_female_ht.describe()
    
    summary_ht = summary_both_ht.union(summary_male_ht).union(summary_female_ht)
    summary_ht = summary_ht.checkpoint(f'{DATA_PATH}/phenotype/summary/{name}_summary.ht', overwrite=overwrite, _read_if_exists= not overwrite)
    summary_ht.describe()
    summary_ht.show(100)
    summary_ht.export(f'{DATA_PATH}/phenotype/summary/{name}_summary.txt.bgz')
    return summary_ht

def process_random_phenotype(overwrite: bool):
    # https://batch.hail.is/batches/8252601/jobs/1
    hl.init(
        app_name=f'Summarize_binary_phenotypes_v8',
        gcs_requester_pays_configuration='aou-neale-gwas',
        master='local[32]',
        tmp_dir=TMP_BUCKET,
        worker_memory="highmem",
        worker_cores=8,
        default_reference="GRCh38",
        log="/summarize_binary_phenotypes.log"
    )
    meta_ht = hl.read_table(f'{DATA_PATH}/utils/aou_v8_sample_meta.ht')
    if (
        not hl.hadoop_exists(
            f'{DATA_PATH}/phenotype/random_pheno.ht/_SUCCESS'
        )
        or overwrite
    ):
        full_ht = hl.import_table(
            f'{DATA_PATH}/phenotype/raw/random_pheno/random_pheno_{ANCESTRIES[0].lower()}.tsv',
            key="userId",
            impute=True,
            types={"userId": hl.tstr},
        )
        full_ht = full_ht.checkpoint(
            f'{TMP_BUCKET}/v8/phenotype/random_pheno_{ANCESTRIES[0]}.ht', overwrite=True
        )
        for ancestry in ANCESTRIES[1:]:
            print(ancestry)
            ht = hl.import_table(
                f'{DATA_PATH}/phenotype/raw/random_pheno/random_pheno_{ancestry.lower()}.tsv',
                key="userId",
                impute=True,
                types={"userId": hl.tstr},
            )
            full_ht = full_ht.union(ht)
            full_ht = full_ht.checkpoint(
                f'{TMP_BUCKET}/v8/phenotype/random_pheno_{ancestry}.ht', overwrite=True
            )
        full_ht = hl.read_table(
            f'{TMP_BUCKET}/v8/phenotype/random_pheno_{ANCESTRIES[-1]}.ht'
        )
        
        full_ht.describe()
        print(full_ht.count())
        full_ht = full_ht.checkpoint(
            f'{DATA_PATH}/phenotype/raw/random_pheno.ht',
            overwrite=overwrite,
            _read_if_exists=not overwrite
        )
    full_ht = full_ht.join(meta_ht, how="left")
    full_ht = full_ht.checkpoint(
        f'{DATA_PATH}/phenotype/random_pheno_annotated.ht',
        overwrite=overwrite,
        _read_if_exists=not overwrite
    )
    full_ht.describe()
    print(full_ht.count())
    print(full_ht.aggregate(hl.agg.counter(full_ht.ancestry)))
    return full_ht

def export_pheno_by_anc_dict(overwrite=False):
    # https://batch.hail.is/batches/8253815/jobs/1 
    hl.init(
        app_name=f'Export_phenotype_dictionary_v8',
        gcs_requester_pays_configuration='aou-neale-gwas',
        master='local[32]',
        tmp_dir=TMP_BUCKET,
        worker_memory="highmem",
        worker_cores=8,
        default_reference="GRCh38",
        log="/export_phenotype_dictionary.log"
    )
    if not hfs.exists(f'{DATA_PATH}/phenotype/summary/pheno_dict_raw.dict') or overwrite:
        all_phenos_by_group = {}
        raw_lab_ht = hl.read_table(f'{DATA_PATH}/phenotype/lab_measurement/lab_measurement_full_annotated.ht')
        LAB_CODES = raw_lab_ht.aggregate(hl.agg.collect_as_set(raw_lab_ht.measurement_concept_id))
        LAB_CODES = [f'{x}_{stat}' for x in LAB_CODES for stat in ['median', 'min', 'max', 'latest']]
        all_phenos_by_group['lab_measurement'] = LAB_CODES
        raw_random_pheno_ht = hl.read_table(f'{DATA_PATH}/phenotype/raw/random_pheno.ht')
        random_phenos_to_run = list(raw_random_pheno_ht.row_value)
        all_phenos_by_group['random_pheno'] = random_phenos_to_run + [f'{pheno}_male' for pheno in random_phenos_to_run] + [f'{pheno}_female' for pheno in random_phenos_to_run]
        raw_gate_pheno_ht = hl.import_table(f'{DATA_PATH}/phenotype/summary/INCLUDE_GWAS_50+_cases.tsv', no_header=True)
        all_phenos_by_group['onset'] = [pheno for pheno in raw_gate_pheno_ht.aggregate(hl.agg.collect(raw_gate_pheno_ht.f0)) if pheno.startswith('Birth')]
        all_phenos_by_group['progression'] = [pheno for pheno in raw_gate_pheno_ht.aggregate(hl.agg.collect(raw_gate_pheno_ht.f0)) if not pheno.startswith('Birth')]
        
        for category in PHENO_CATEGORIES:
            print(f'Loading {category} information......')
            if category in ['lab_measurement', 'random_pheno', 'onset', 'progression']:
                continue
            ht = hl.read_table(f'{DATA_PATH}/phenotype/summary/{category}_summary.ht')
            all_phenos_by_group[category] = set(ht.phenoname.collect())
            write_pickle_dict(output=f'{DATA_PATH}/phenotype/summary/pheno_dict_raw.dict', dict =all_phenos_by_group)

    all_phenos_by_group = read_pickle_dict(f'{DATA_PATH}/phenotype/summary/pheno_dict_raw.dict')
    print(f'----------Number of phenotypes per category (RAW): --------------')
    print([(category, len(all_phenos_by_group[category])) for category in all_phenos_by_group.keys()])

    quantitative_traits = []
    binary_traits = []
    if 'random_pheno' in all_phenos_by_group.keys():
        quantitative_traits = [pheno for pheno in all_phenos_by_group['random_pheno'] if 'continuous' in pheno]
        binary_traits = [pheno for pheno in all_phenos_by_group['random_pheno'] if 'continuous' not in pheno]
    for category in QUANTITATIVE_CATEGORIES:
        quantitative_traits = quantitative_traits + list(all_phenos_by_group[category])
    for category in BINARY_CATEGORIES:
        binary_traits = binary_traits + list(all_phenos_by_group[category])
    print(f'Number of quantitative traits: {len(quantitative_traits)}\nNumber of binary traits: {len(binary_traits)}')

    if not hfs.exists(f'{DATA_PATH}/phenotype/summary/pheno_by_ancestry_dict.dict') or overwrite:
        phenos_to_run_by_ancestry = {}
        phenos_to_run_by_ancestry_by_group = {}
        lab_summary_ht = hl.read_table(f'{DATA_PATH}/phenotype/summary/lab_measurement_summary.ht')
        for ancestry in ANCESTRIES:
            phenos_to_run_by_ancestry[ancestry.lower()] = []
            phenos_to_run_by_ancestry_by_group[ancestry.lower()] = {}
            lab_phenos_by_ancestry = lab_summary_ht.filter((lab_summary_ht.ancestry == ancestry.upper()) & (lab_summary_ht.n_cases >= 200))
            lab_phenos_by_ancestry = lab_phenos_by_ancestry.phenoname.collect()
            lab_phenos_by_ancestry = [f'{x}_{stat}' for x in lab_phenos_by_ancestry for stat in ['median', 'min', 'max', 'latest']]
            phenos_to_run_by_ancestry_by_group[ancestry.lower()]['lab_measurement'] = lab_phenos_by_ancestry
            phenos_to_run_by_ancestry[ancestry.lower()] = phenos_to_run_by_ancestry[ancestry.lower()] + lab_phenos_by_ancestry
            if 'random_pheno' in all_phenos_by_group.keys():
                random_pheno_lst = [x for x in all_phenos_by_group['random_pheno'] if x.startswith('random_0.5') and x.endswith(('_1', '_2', '_3', '_4', '_5'))]
                phenos_to_run_by_ancestry_by_group[ancestry.lower()]['random_pheno'] = random_pheno_lst
                phenos_to_run_by_ancestry[ancestry.lower()] = phenos_to_run_by_ancestry[ancestry.lower()] + random_pheno_lst
        for category in GATE_CATEGORIES:
            summary_ht = hl.read_table(f'{DATA_PATH}/phenotype/summary/{category}_summary_processed.ht')
            summary_ht = summary_ht.filter(summary_ht.n_cases >= 200) 
            summary_ht = summary_ht.checkpoint(f'{TMP_BUCKET}/v8/phenotype/summary/{category}_n_cases_over_200.ht', overwrite=True)
            phenos = list(summary_ht.phenoname.collect())
            for ancestry in ANCESTRIES:
                ancestry_phenos = [pheno for pheno in phenos if ancestry.lower() in pheno]
                phenos_to_run_by_ancestry_by_group[ancestry.lower()][category] = ancestry_phenos
                phenos_to_run_by_ancestry[ancestry.lower()] = phenos_to_run_by_ancestry[ancestry.lower()] + ancestry_phenos

        for category in PHENO_CATEGORIES:
            if category in ['lab_measurement', 'random_pheno', 'onset', 'progression']: continue
            pheno_ht = hl.read_table(f'{DATA_PATH}/phenotype/summary/{category}_summary.ht')
            for ancestry in ANCESTRIES:
                print(f'--------Loading {category} phenotypes for {ancestry.upper()}--------')
                ancestry_pheno_ht = pheno_ht.filter(pheno_ht[f'n_cases_{ancestry}'] >= 200)
                if category == 'r_drug':
                    ancestry_pheno_ht = ancestry_pheno_ht.filter(ancestry_pheno_ht.r_drug_curation_80)
                ancestry_pheno_ht = ancestry_pheno_ht.checkpoint(f'{TMP_BUCKET}/v8/phenotype/summary/{ancestry}_{category}_n_cases_over_200.ht', overwrite=True)
                pheno_lst = ancestry_pheno_ht.phenoname.collect()
                phenos_to_run_by_ancestry_by_group[ancestry.lower()][category] = pheno_lst
                phenos_to_run_by_ancestry[ancestry.lower()] = phenos_to_run_by_ancestry[ancestry.lower()] + pheno_lst
        meta_pheno_lst = []
        for ancestry in ANCESTRIES:
            meta_pheno_lst = meta_pheno_lst + phenos_to_run_by_ancestry[ancestry.lower()]
        phenos_to_run_by_ancestry['meta'] = list(set(meta_pheno_lst))
        print(f"META pheno length: {len(phenos_to_run_by_ancestry['meta'])}")

        write_pickle_dict(output=f'{DATA_PATH}/phenotype/summary/pheno_by_ancestry_by_group_dict.dict',
                          dict=phenos_to_run_by_ancestry_by_group)
        write_pickle_dict(output=f'{DATA_PATH}/phenotype/summary/pheno_by_ancestry_dict.dict',
                          dict=phenos_to_run_by_ancestry)
        
        print(f'----------Number of phenotypes per category per ancestry (filtered to n_cases >= 200): --------------')
        print([f"{ancestry.upper()}-{category}: {len(phenos_to_run_by_ancestry_by_group[ancestry.lower()][category])}" for ancestry in ANCESTRIES for category in list(all_phenos_by_group.keys())])
        print([f"{ancestry.upper()}: {len(phenos_to_run_by_ancestry[ancestry.lower()])}" for ancestry in ANCESTRIES]) 

def update_call_stats_ht(ancestry, data_type, pruned, overwrite):
    if data_type not in ['ACAF', 'Exome']:
        raise ValueError(f"Invalid data type: {data_type}")
    callstats_ht_path = f'{DATA_PATH}/utils/call_stats/{data_type.lower() if data_type == "Exome" else data_type}_{"pruned" if pruned else "pre_pruning"}/{ancestry.upper()}_{data_type.lower() if data_type == "Exome" else data_type}_call_stats.ht'
    if not hfs.exists(f'{callstats_ht_path}/_SUCCESS') or overwrite:
        print(
            f"-------------Computing call stats (Ancestry: {ancestry.upper()} {data_type} ({'pruned' if pruned else 'pre_pruning'}))-----------"
        )
        defined_ht = hl.read_table(f'{DATA_PATH}/utils/aou_v8_hard_filter.ht')
        mt = get_filtered_mt(mt_type = data_type, sample_ids=defined_ht, filter_variants=True, filter_samples=True,
                             adj_filter=True, ancestry=ancestry.lower(), prune_samples=pruned)
        call_stats_ht = mt.annotate_rows(call_stats=hl.agg.call_stats(mt.GT, mt.alleles)).rows()
        call_stats_ht = call_stats_ht.select('call_stats')
        call_stats_ht = call_stats_ht.naive_coalesce(1000).checkpoint(callstats_ht_path, overwrite=True)
        call_stats_ht.describe()
    call_stats_ht = hl.read_table(callstats_ht_path)
    return call_stats_ht


def main(args):
    try:
        if args.batch:
            backend = hb.ServiceBackend(
                billing_project="all-by-aou",
                remote_tmpdir=TMP_BUCKET,
            )
            b = hb.Batch(
                name=f"phenotype_process_v8",
                requester_pays_project="aou-neale-gwas",
                default_python_image='us-central1-docker.pkg.dev/aou-neale-gwas/hail-batch/hail-python:0.2.133',
                backend=backend,
            )
        
        if args.run_vep:
            hl.init(
                app_name=f'VEP_v8',
                tmp_dir=TMP_BUCKET,
                default_reference="GRCh38",
                log=f"/VEP_v8.log"
            )
            # variant_table_tag = 'aou_indel_vat_v8'
            variant_table_tag = 'aou_vds_variant_data_row_v8'
            if not hl.hadoop_exists(f'{DATA_PATH}/vat/{variant_table_tag}_vep.ht') or args.overwrite:
                ht = hl.read_table(f'{DATA_PATH}/vat/{variant_table_tag}.ht')
                ht = hl.split_multi(ht)
                ht = ht.checkpoint(f'{DATA_PATH}/vat/{variant_table_tag}_split_multi.ht', overwrite=args.overwrite)
                indel_ht = ht.filter(~hl.is_snp(ht.alleles[0], ht.alleles[1]))
                indel_ht = indel_ht.naive_coalesce(5000).checkpoint(f'{DATA_PATH}/vat/{variant_table_tag}_indel_split_multi.ht', overwrite=args.overwrite)
                indel_ht = hl.read_table(f'{DATA_PATH}/vat/{variant_table_tag}_indel_split_multi.ht', _n_partitions=50000)
                vep_ht = vep_or_lookup_vep(indel_ht, vep_version="vep105")
                vep_ht.describe()
                vep_ht.write(
                    f'{DATA_PATH}/vep/{variant_table_tag}_vep.ht',
                    overwrite=args.overwrite,
                )
            try:
                vep_ht = hl.read_table(f'{DATA_PATH}/vep/{variant_table_tag}_vep.ht')
                vep_ht.describe()
                print(vep_ht.count())
                vep_ht.show()
            except:
                pass
            vep_ht.show()

        if args.check_path:
            j = b.new_python_job(name=f"Check aou path")
            j.call(check_path, path=args.path_to_check)

        if args.copy_files:
            j = b.new_python_job(name=f"Copy files")
            j.call(copy_files, args.original_path, args.target_path)
        
        if args.process_vat:
            j = b.new_python_job(name=f"Export VAT HT")
            j.call(process_vat, args.overwrite)
        
        if args.process_vds:
            j = b.new_python_job(name=f"Export VDS")
            j.memory('highmem')
            j.cpu(8) 
            j.env('PYSPARK_SUBMIT_ARGS', '--driver-memory 32g --executor-memory 32g pyspark-shell')
            j.call(process_vds, args.overwrite)

        if args.process_pca_data:
            j = b.new_python_job(name=f"Export PCA HT")
            j.memory('highmem')
            j.cpu(8) 
            j.env('PYSPARK_SUBMIT_ARGS', '--driver-memory 32g --executor-memory 32g pyspark-shell')
            j.call(process_pca, args.overwrite)
        
        if args.process_hard_filter_data:
            j = b.new_python_job(name=f"Export hard filter HT")
            j.call(process_hard_filter_data, args.overwrite)
        
        if args.process_demographic_data:
            j = b.new_python_job(name=f"Export demographic HT")
            j.call(process_demographic_data, args.overwrite)

        if args.process_related_data:
            j = b.new_python_job(name=f"Export related HT")
            j.call(process_related_data, args.overwrite)

        if args.update_meta_data:
            j = b.new_python_job(name=f"Update meta HT")
            j.call(update_meta_ht, args.overwrite)

        if args.update_random_phenotype:
            j = b.new_python_job(name=f"Update random phenotype HT")
            j.memory('highmem')
            j.cpu(8) 
            j.env('PYSPARK_SUBMIT_ARGS', '--driver-memory 80g --executor-memory 80g pyspark-shell')
            j.call(process_random_phenotype, args.overwrite)

        if args.export_raw_pheno_ht:
            for name in PHENO_CATEGORIES:
            # for name in FILENAMES:
                if name in ['onset', 'progression']:
                    continue
                if not hfs.exists(f'{DATA_PATH}/phenotype/raw/{name}.ht/_SUCCESS') or args.overwrite:
                    if args.batch:
                        j = b.new_python_job(name=f"Export raw pheno {name} HT")
                        j.memory('16G')
                        j.cpu(8)
                        j.env('PYSPARK_SUBMIT_ARGS', '--driver-memory 8g --executor-memory 8g pyspark-shell')
                        j.call(write_raw_phenotype_ht, name, overwrite=True)
                    else:
                        write_raw_phenotype_ht(name, args.overwrite)
        
        if args.write_raw_gate_summary_ht:
            if args.batch:
                j = b.new_python_job(name=f"Write raw gate summary HT")
                j.call(write_raw_gate_summary_ht, args.overwrite)
            else:
                write_raw_gate_summary_ht(args.overwrite)
        
        if args.write_raw_gate_pheno_ht:
            for name in GATE_CATEGORIES:
                write_gate_pheno_ht(b, name, 200, args.overwrite)
        
        if args.process_physical_measurements:
            if not hl.hadoop_exists(f'{DATA_PATH}/phenotype/physical_measurement.ht/_SUCCESS') or args.overwrite:
                if args.batch:
                    j = b.new_python_job(name=f"Process physical measurements")
                    j.call(process_quantitative_phenotypes, args.overwrite)
                else:
                    process_quantitative_phenotypes(args.overwrite)
        
        if args.extract_mt_from_ht:
            for name in PHENO_CATEGORIES:
                print(f'---------------Extracting {name} MT from HT--------------')
                if not hl.hadoop_exists(f'{DATA_PATH}/phenotype/raw/{name}.mt/_SUCCESS') or args.overwrite:
                    if args.batch:
                        j = b.new_python_job(name=f"Extract {name} MT from HT")
                        j.memory('highmem')
                        j.cpu(8) 
                        j.env('PYSPARK_SUBMIT_ARGS', '--driver-memory 32g --executor-memory 32g pyspark-shell')
                        j.call(extract_mt_from_ht, name, args.overwrite) 
                    else:
                        extract_mt_from_ht(name, args.overwrite)
                mt = hl.read_matrix_table(f'{DATA_PATH}/phenotype/raw/{name}.mt')
                mt.describe()
                mt.show()
                print(mt.count())
                print(mt.aggregate_rows(hl.agg.counter(mt.ancestry)))
        
        if args.summarize_lab_measurement:
            if args.batch:
                j = b.new_python_job(name="Summarize lab measurement")
                j.call(summarize_lab_measurement, args.overwrite)
            else:
                summarize_lab_measurement(args.overwrite)
        
        if args.write_annotated_lab_ht:
            write_annotated_lab_ht(b=b, overwrite=args.overwrite)
              
        if args.summarize_raw_binary_phenotypes:
            for name in BINARY_CATEGORIES:
                job_depend_on = None
                if not hl.hadoop_exists(f'{DATA_PATH}/phenotype/raw/{name}.mt/_SUCCESS') or args.overwrite:
                    if args.batch:
                        j_mt = b.new_python_job(name=f"Export {name} MT")
                        j_mt.memory('highmem')
                        j_mt.cpu(8) 
                        j_mt.env('PYSPARK_SUBMIT_ARGS', '--driver-memory 32g --executor-memory 32g pyspark-shell')
                        j_mt.call(extract_mt_from_ht, name, overwrite=True) 
                        job_depend_on = j_mt
                    else:
                        extract_mt_from_ht(name, args.overwrite)
                if not hl.hadoop_exists(f'{DATA_PATH}/phenotype/summary/{name}_raw_summary.txt.bgz') or args.overwrite:
                    if args.batch:
                        j = b.new_python_job(name=f"Summarize binary phenotypes {name}")
                        if job_depend_on is not None:
                            j.depends_on(job_depend_on)
                        j.call(summarize_raw_binary_phenotypes, name, overwrite=True)
                    else:
                        summarize_raw_binary_phenotypes(name, args.overwrite)

        if args.summarize_quantitative_phenotypes_per_ancestry:
            for name in QUANTITATIVE_CATEGORIES:
                print(f'---------------Summarizing quantitative phenotypes {name}--------------')
                job_depend_on = None
                ext = '.ht' if name in ['random_pheno'] else '.mt'
                if not hl.hadoop_exists(f'{DATA_PATH}/phenotype/{name}_annotated{ext}/_SUCCESS') or args.overwrite:
                    if args.batch:
                        j_annotate = b.new_python_job(name=f"Annotate {name} {ext.upper()}")
                        j_annotate.call(annotate_phenotype_tables, name, overwrite=True)
                        job_depend_on = j_annotate
                    else:
                        annotate_phenotype_tables(name, overwrite=True)
                if not hl.hadoop_exists(f'{DATA_PATH}/phenotype/summary/{name}_summary.txt.bgz') or args.overwrite:
                    if args.batch:
                        j = b.new_python_job(name=f"Summarize quantitative phenotypes {name}")
                        if job_depend_on is not None:
                            j.depends_on(job_depend_on)
                        j.memory('highmem')
                        j.cpu(8) 
                        j.env('PYSPARK_SUBMIT_ARGS', '--driver-memory 32g --executor-memory 32g pyspark-shell')
                        j.call(summarize_quantitative_phenotypes, name, args.overwrite)
                    else:
                        summarize_quantitative_phenotypes(name, args.overwrite)
        
        if args.summarize_binary_phenotypes_per_ancestry:
            BINARY_CATEGORIES = ['mcc2_phecode', 'mcc2_phecodex']
            for name in BINARY_CATEGORIES:
                print(f'---------------Summarizing binary phenotypes {name} per ancestry--------------')
                job_depend_on = None
                if not hl.hadoop_exists(f'{DATA_PATH}/phenotype/{name}_annotated.mt/_SUCCESS') or args.overwrite:
                    if args.batch:
                        j_annotate = b.new_python_job(name=f"Annotate {name} MT")
                        j_annotate.call(annotate_phenotype_tables, name, overwrite=True)
                        job_depend_on = j_annotate
                    else:
                        annotate_phenotype_tables(name, overwrite=True)
                if args.batch:
                    j = b.new_python_job(name=f"Summarize binary phenotypes {name} per ancestry")
                    if job_depend_on is not None:
                        j.depends_on(job_depend_on)
                    j.memory('highmem')
                    j.cpu(8) 
                    j.env('PYSPARK_SUBMIT_ARGS', '--driver-memory 32g --executor-memory 32g pyspark-shell')
                    j.call(summarize_binary_phenotypes, name, args.overwrite)
                else:
                    summarize_binary_phenotypes(name, args.overwrite)
                # mt = hl.read_matrix_table(f'{DATA_PATH}/phenotype/{name}_annotated.mt') 
                # mt.describe()
                # mt.show()
                # print(mt.count())
                # print(mt.aggregate_rows(hl.agg.counter(mt.ancestry)))
                # print(mt.aggregate_cols(hl.agg.collect(mt.phenoname))) 
                
        if args.export_pheno_by_anc_dict:
            # https://batch.hail.is/batches/8254763/jobs/1
            if args.batch:
                j = b.new_python_job(name=f"Export phenotype dictionary")
                j.memory('highmem')
                j.cpu(8) 
                j.env('PYSPARK_SUBMIT_ARGS', '--driver-memory 8g --executor-memory 8g pyspark-shell')
                j.call(export_pheno_by_anc_dict, args.overwrite)
            else:
                export_pheno_by_anc_dict(args.overwrite)
        
        if args.update_saige_intervals:
            if args.batch:
                j = b.new_python_job(name=f"Update SAIGE variant intervals")
                j.memory('highmem')
                j.cpu(8) 
                j.env('PYSPARK_SUBMIT_ARGS', '--driver-memory 8g --executor-memory 8g pyspark-shell')
                j.call(update_variant_interval, CHUNK_SIZE['afr'], args.overwrite)
                j.call(update_variant_interval, CHUNK_SIZE['eur'], args.overwrite)

                j = b.new_python_job(name=f"Update SAIGE gene intervals")
                j.memory('highmem')
                j.cpu(8) 
                j.env('PYSPARK_SUBMIT_ARGS', '--driver-memory 8g --executor-memory 8g pyspark-shell')
                j.call(update_gene_interval, N_GENE_PER_GROUP, args.overwrite)
            else:
                update_variant_interval(CHUNK_SIZE['afr'], args.overwrite)
                update_variant_interval(CHUNK_SIZE['eur'], args.overwrite)
                update_gene_interval(N_GENE_PER_GROUP, args.overwrite)
        
        if args.update_call_stats_ht:
            hl.stop()
            hl.init(
                app_name=f'Update_call_stats_v8',
                gcs_requester_pays_configuration='aou-neale-gwas',
                master='local[32]',
                tmp_dir=TMP_BUCKET,
                worker_memory="highmem",
                worker_cores=8,
                default_reference="GRCh38",
                log=f"/update_call_stats_v8.log"
            )       
            for data_type in ['ACAF', 'Exome']:
                print(f'Updating call stats {data_type}...')
                for ancestry in ANCESTRIES:
                    if args.batch:
                        print(f'Updating call stats {ancestry} {data_type} pre_pruning...')
                        j1 = b.new_python_job(name=f"Update call stats {ancestry} {data_type} pre_pruning")
                        j1.memory('highmem')
                        j1.cpu(8) 
                        j1.env('PYSPARK_SUBMIT_ARGS', '--driver-memory 24g --executor-memory 24g pyspark-shell')
                        j1.call(update_call_stats_ht, ancestry, data_type, pruned=False, overwrite=args.overwrite)
                        if hfs.exists(f'{DATA_PATH}/utils/pca/results/{ancestry.lower()}_pca_centroid_pruned.ht/_SUCCESS'):
                            print(f'Updating call stats {ancestry} {data_type} pruned...')
                            j2 = b.new_python_job(name=f"Update call stats {ancestry} {data_type} pruned")
                            j2.memory('highmem')
                            j2.cpu(8) 
                            j2.env('PYSPARK_SUBMIT_ARGS', '--driver-memory 24g --executor-memory 24g pyspark-shell')
                            j2.call(update_call_stats_ht, ancestry, data_type, pruned=True, overwrite=args.overwrite)
                    else:
                        print(f'Updating call stats {ancestry} {data_type} pre_pruning...')
                        update_call_stats_ht(ancestry, data_type, pruned=False, overwrite=args.overwrite)
                        if hfs.exists(f'{DATA_PATH}/utils/pca/results/{ancestry.lower()}_pca_centroid_pruned.ht/_SUCCESS') and not hfs.exists(f'{DATA_PATH}/utils/call_stats/{data_type}_pruned/{ancestry.upper()}_{data_type}_call_stats.ht/_SUCCESS'):
                            print(f'Updating call stats {ancestry} {data_type} pruned...')
                            update_call_stats_ht(ancestry, data_type, pruned=True, overwrite=args.overwrite)

        if args.batch:
            b.run()

    finally:
        if not args.batch:
            from datetime import date

            hl.copy_log(f"{MY_BUCKET}/log/preprocessing_{date.today()}.log")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--run-vep",
        help="Run VEP",
        action="store_true",
    )
    parser.add_argument(
        "--check-path",
        help="Check the path",
        action="store_true",
    )
    parser.add_argument(
        "--copy-files",
        help="Copy files",
        action="store_true",
    )
    parser.add_argument(
        "--original-path",
        help="The original path",
        type=str,
        default=f'{ORIGINAL_PCA_PATH}/results/pca_SAS_pruned_scores_full.txt.bgz'
    )
    parser.add_argument(
        "--target-path",
        help="The target path",
        type=str,
        default=f'{DATA_PATH}/utils/pca/pca_SAS_pruned_scores_full.txt.bgz'
    )
    parser.add_argument(
        "--path-to-check",
        help="The path to check",
        type=str,
        default=ORIGINAL_DATA_ROOT
    )
    parser.add_argument(
        "--export-raw-pheno-ht",
        help="Export raw phenotype HT",
        action="store_true",
    )
    parser.add_argument(
        "--write-raw-gate-summary-ht",
        help="Write raw gate summary HT",
        action="store_true",
    )
    parser.add_argument(
        "--write-raw-gate-pheno-ht",
        help="Write raw gate pheno HT",
        action="store_true",
    )
    parser.add_argument(
        "--process-vat",
        help="Process VAT table",
        action="store_true",
    )
    parser.add_argument(
        "--process-pca-data",
        help="Process PCA data",
        action="store_true",
    )
    parser.add_argument(
        "--process-demographic-data",
        help="Process demographic data",
        action="store_true",
    )
    parser.add_argument(
        "--process-related-data",
        help="Process related data",
        action="store_true",
    )
    parser.add_argument(
        "--process-hard-filter-data",
        help="Process hard filter data",
        action="store_true",
    )
    parser.add_argument(
        "--update-meta-data",
        help="Update meta data",
        action="store_true",
    )
    parser.add_argument(
        "--update-random-phenotype",
        help="Update random phenotype HT",
        action="store_true",
    )
    parser.add_argument(
        "--extract-mt-from-ht",
        help="Extract MT from HT",
        action="store_true",
    )
    parser.add_argument(
        "--process-physical-measurements",
        help="Process physical measurements",
        action="store_true",
    )
    parser.add_argument(
        "--summarize-raw-binary-phenotypes",
        help="Summarize raw binary phenotypes",
        action="store_true",
    )
    parser.add_argument(
        "--summarize-lab-measurement",
        help="Summarize lab measurement",
        action="store_true",
    )
    parser.add_argument(
        "--write-annotated-lab-ht",
        help="Write annotated lab HT",
        action="store_true",
    )
    parser.add_argument(
        "--summarize-quantitative-phenotypes-per-ancestry",
        help="Summarize quantitative phenotypes per ancestry",
        action="store_true",
    )
    parser.add_argument(
        "--summarize-binary-phenotypes-per-ancestry",
        help="Summarize binary phenotypes per ancestry",
        action="store_true",
    )
    parser.add_argument(
        "--export-pheno-by-anc-dict",
        help="Export phenotype dictionary",
        action="store_true",
    )
    parser.add_argument(
        "--batch",
        help="Run the pipeline in Batch, not QoB or dataproc",
        action="store_true",
    )
    parser.add_argument(
        "--overwrite",
        help="Whether to overwrite the existing table",
        action="store_true",
    )
    parser.add_argument(
        "--update-saige-intervals",
        help="Update SAIGE intervals",
        action="store_true",
    )
    parser.add_argument(
        "--update-call-stats-ht",
        help="Update call stats HT",
        action="store_true",
    )
    parser.add_argument(
        "--process-vds",
        help="Process VDS",
        action="store_true",
    )
    args = parser.parse_args()

    main(args)
