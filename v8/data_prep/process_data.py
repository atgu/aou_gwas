#!/usr/bin/env python3

__author__ = "wlu"

from sre_parse import CATEGORIES
import sys
print(sys.path)
# from aou_gwas import *  # for dataproc
# from utils.utils import *  # for QoB
# from utils.resources import *  # for QoB
import hailtop.batch as hb
import argparse
import hail as hl
import hailtop.fs as hfs
from gnomad.utils.vep import *
from gnomad.utils.filtering import *

TRANCHE = "v8"
ANCESTRIES = ['AFR', 'AMR', 'EAS', 'EUR', 'MID', 'SAS']
ANALYSIS_BUCKET = "gs://aou_analysis"
MY_BUCKET = 'gs://aou_wlu'
TMP_BUCKET = 'gs://aou_tmp'
DATA_PATH = f"{ANALYSIS_BUCKET}/{TRANCHE}/data"
ORIGINAL_GENO_ROOT = 'gs://fc-aou-datasets-controlled/v8' # https://batch.hail.is/batches/8251999 
ORIGINAL_GENO_DIR = f'{ORIGINAL_GENO_ROOT}/wgs/short_read' # https://batch.hail.is/batches/8252511/jobs/1 
ORIGINAL_AUX_DIR = f'{ORIGINAL_GENO_DIR}/snpindel/aux' # https://batch.hail.is/batches/8252513/jobs/1 
ORIGINAL_DEMO_PATH = f'{ORIGINAL_AUX_DIR}/qc/genomic_metrics.tsv' 
ORIGINAL_ALL_SAMPLE_PATH = f'{ORIGINAL_AUX_DIR}/qc/all_samples.tsv'
ORIGINAL_RELATED_PATH = f'{ORIGINAL_AUX_DIR}/relatedness/relatedness_flagged_samples.tsv'
ORIGINAL_KINSHIP_PATH = f'{ORIGINAL_AUX_DIR}/relatedness/relatedness.tsv'
ORIGINAL_GLOBAL_PCA_PATH = f'{ORIGINAL_AUX_DIR}/ancestry/ancestry_preds.tsv'
ORIGINAL_GLOBAL_LOADING_PATH = f'{ORIGINAL_AUX_DIR}/ancestry/loadings.ht'

ORIGINAL_DATA_ROOT = 'gs://allxall-phenotypes-v2' # https://batch.hail.is/batches/8251082/jobs/1
ORIGINAL_PHENO_DIR = f'{ORIGINAL_DATA_ROOT}/2025-04-14'
ORIGINAL_VAT_PATH = f'{ORIGINAL_DATA_ROOT}/2025-04-17/vat'
ORIGINAL_PCA_PATH = f'{ORIGINAL_DATA_ROOT}/2025-04-18/pca' # https://batch.hail.is/batches/8251124/jobs/1
ORIGINAL_CALLSTATS_PATH = f'{ORIGINAL_DATA_ROOT}/2025-04-16/call_stats'
ORIGINAL_GRM_PATH = f'{ORIGINAL_DATA_ROOT}/2025-04-16/grm'
ORIGINAL_RANDOM_PATH = f'{ORIGINAL_DATA_ROOT}/2025-04-16/random_pheno' # https://batch.hail.is/batches/8252587/jobs/1
ORIGINAL_HARD_FILTER = f'{ORIGINAL_DATA_ROOT}/2025-04-16/hard_filtered_samples.ht'

PHENO_CATEGORIES = ['physical_measurement', 'r_drug', 'pfhh_survey']
BINARY_CATEGORIES = ['r_drug', 'pfhh_survey']

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
    back = hl.current_backend()
    spark = back._spark_session
    parquetFile = spark.read.parquet(f"{ORIGINAL_PHENO_DIR}/{name}_table.parquet")
    ht = hl.Table.from_spark(parquetFile)
    ht = ht.key_by('person_id')
    ht = ht.checkpoint(f'{DATA_PATH}/phenotype/raw/{name}.ht', overwrite=overwrite, _read_if_exists= not overwrite)
    ht.describe()
    # ht.show()
    print(f"Number of samples: {ht.count()}")
    print(f"Number of phenotypes: {len(ht.row_value)}")
    return ht

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
    mt = hl.read_matrix_table(f'{DATA_PATH}/phenotype/raw/{name}.mt')
    mt = mt.annotate_cols(n_cases = hl.agg.count_where(mt.value),
                      n_controls = hl.agg.count_where(~mt.value),
                      n_missing = hl.agg.count_where(hl.is_missing(mt.value)))
    if not hl.hadoop_exists(f'{DATA_PATH}/phenotype/summary/{name}_summary.ht/_SUCCESS'):
        overwrite = True
    summary_ht = mt.cols().checkpoint(f'{DATA_PATH}/phenotype/summary/{name}_summary.ht', overwrite=overwrite, _read_if_exists= not overwrite)
    if not hl.hadoop_exists(f'{DATA_PATH}/phenotype/summary/{name}_summary.txt.bgz'):
        summary_ht.export(f'{DATA_PATH}/phenotype/summary/{name}_summary.txt.bgz')
    return summary_ht  

def extract_mt_from_ht(name = 'r_drug', overwrite = False):
    ht = hl.read_table(f'{DATA_PATH}/phenotype/raw/{name}.ht')
    columns = list(ht.row_value)
    ht2 = ht.select(values=[hl.struct(value=ht[x])for x in columns]).annotate_globals(
        columns=hl.map(lambda x: hl.struct(phenoname=x),
                       hl.literal(columns)))
    ht2 = ht2.checkpoint(f'{TMP_BUCKET}/v8/{name}.ht', overwrite=True)
    mt = ht2._unlocalize_entries('values', 'columns', [])
    mt = mt.checkpoint(f'{DATA_PATH}/phenotype/raw/{name}.mt', overwrite=overwrite, _read_if_exists= not overwrite)
    return mt

def process_vat(overwrite=False):
    # https://batch.hail.is/batches/8250503/jobs/1
    ht = hl.read_table(f'{ORIGINAL_VAT_PATH}/aou_PARSED_SORTED_COLLECTED_vat_wlu.ht')
    ht = ht.checkpoint(f'{DATA_PATH}/vat/aou_PARSED_SORTED_COLLECTED_vat_v8.ht', overwrite=overwrite, _read_if_exists= not overwrite)
    return ht

def process_pca(overwrite=False):
    # https://batch.hail.is/batches/8252515/jobs/1
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
        ht = hl.import_table(f'{DATA_PATH}/utils/pca/pca_{ANCESTRIES[0]}_pruned_scores_full.txt.bgz', key='s')
        ht = ht.annotate(ancestry = ANCESTRIES[0])
        for ancestry in ANCESTRIES[1:]:
            tmp_ht = hl.import_table(f'{DATA_PATH}/utils/pca/pca_{ancestry}_pruned_scores_full.txt.bgz', key='s')   
            tmp_ht = tmp_ht.annotate(ancestry = ancestry)
            ht = ht.union(tmp_ht)
        ht = ht.checkpoint(f'{DATA_PATH}/utils/pca/pca_pruned_scores_full.ht', overwrite=overwrite, _read_if_exists= not overwrite)
    ht = hl.read_table(f'{DATA_PATH}/utils/pca/pca_pruned_scores_full.ht')
    ht.describe()
    print(f"Number of samples: {ht.count()}")
    print(f"Number of ancestries: {ht.aggregate(hl.agg.counter(ht.ancestry))}")

    ht = hl.import_table(f'{ORIGINAL_GLOBAL_PCA_PATH}', impute=True, key='research_id', types = {'research_id':hl.tstr})
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
    ht = hl.import_table(f'{ORIGINAL_DEMO_PATH}', key='research_id')
    ht = ht.checkpoint(f'{DATA_PATH}/utils/aou_v8_demographic.ht', overwrite=overwrite, _read_if_exists= not overwrite)
    ht.describe()
    print(f"Number of samples: {ht.count()}")

    ht = hl.import_table(f'{ORIGINAL_ALL_SAMPLE_PATH}', key='s')
    ht = ht.checkpoint(f'{DATA_PATH}/utils/aou_v8_all_samples.ht', overwrite=overwrite, _read_if_exists= not overwrite)
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

    ht = hl.import_table(f'{ORIGINAL_KINSHIP_PATH}', key=['i.s', 'j.s'])
    ht = ht.checkpoint(f'{DATA_PATH}/utils/aou_v8_related_kinship.ht', overwrite=overwrite, _read_if_exists= not overwrite)
    ht.describe()
    print(f"Number of samples: {ht.count()}")
    return ht


def update_meta_ht(overwrite: bool):
    # https://batch.hail.is/batches/8252557/jobs/2
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
    return meta_ht

def annotate_phenotype_tables(name, overwrite=False):
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
    if name in ['physical_measurement']:
        data = hl.read_table(f'{DATA_PATH}/phenotype/{name}.ht')
        data.describe()
        data = data.key_by()
        data = data.annotate(person_id = hl.str(data.person_id))
        data = data.key_by('person_id')
        data = data.annotate(**meta_ht[data.key])
        data = data.checkpoint(f'{DATA_PATH}/phenotype/{name}_annotated.ht', overwrite=overwrite, _read_if_exists= not overwrite)
    else:
        data = hl.read_matrix_table(f'{DATA_PATH}/phenotype/raw/{name}.mt')
        data.describe()
        data = data.key_rows_by()
        data = data.annotate_rows(person_id = hl.str(data.person_id))
        data = data.key_rows_by('person_id')
        data = data.annotate_rows(**meta_ht[data.row_key])
        data = data.checkpoint(f'{DATA_PATH}/phenotype/{name}_annotated.mt', overwrite=overwrite, _read_if_exists= not overwrite)
    return data

def summarize_binary_phenotypes(name, overwrite):
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
    mt = mt.annotate_cols(n_cases = hl.agg.count_where(mt.value),
                      n_controls = hl.agg.count_where(~mt.value),
                      n_missing = hl.agg.count_where(hl.is_missing(mt.value)),
                      n_cases_by_ancestry=hl.agg.group_by(mt.ancestry, hl.agg.count_where(mt.value)),
                      n_controls_by_ancestry=hl.agg.group_by(mt.ancestry, hl.agg.count_where(~mt.value)),
                      n_samples_defined_by_ancestry=hl.agg.group_by(mt.ancestry, hl.agg.count_where(hl.is_defined(mt.value))),
                      )
    mt = mt.annotate_cols(**{f'n_cases_{ancestry}': mt.n_cases_by_ancestry.get(ancestry) for ancestry in ANCESTRIES})

    if not hl.hadoop_exists(f'{DATA_PATH}/phenotype/summary/{name}_summary_per_ancestry.ht/_SUCCESS'):
        overwrite = True
    summary_ht = mt.cols().checkpoint(f'{DATA_PATH}/phenotype/summary/{name}_summary_per_ancestry.ht', overwrite=overwrite, _read_if_exists= not overwrite)
    if not hl.hadoop_exists(f'{DATA_PATH}/phenotype/summary/{name}_summary_per_ancestry.txt.bgz'):
        summary_ht.export(f'{DATA_PATH}/phenotype/summary/{name}_summary_per_ancestry.txt.bgz')
    summary_ht = hl.read_table(f'{DATA_PATH}/phenotype/summary/{name}_summary_per_ancestry.ht')
    summary_ht.describe()
    summary_ht.show(100)
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
            f'{ORIGINAL_RANDOM_PATH}/random_pheno_{ANCESTRIES[0].lower()}.tsv',
            key="userId",
            impute=True,
            types={"userId": hl.tstr},
        )
        full_ht = full_ht.checkpoint(
            f'{TMP_BUCKET}/v8/phenotype/random_pheno_{ANCESTRIES[0]}.ht', _read_if_exists=True
        )
        for ancestry in ANCESTRIES[1:]:
            print(ancestry)
            ht = hl.import_table(
                f'{ORIGINAL_RANDOM_PATH}/random_pheno_{ancestry.lower()}.tsv',
                key="userId",
                impute=True,
                types={"userId": hl.tstr},
            )
            full_ht = full_ht.union(ht)
            full_ht = full_ht.checkpoint(
                f'{TMP_BUCKET}/v8/phenotype/random_pheno_{ancestry}.ht', _read_if_exists=True
            )
        full_ht = hl.read_table(
            f'{TMP_BUCKET}/v8/phenotype/random_pheno_SAS.ht'
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
    return full_ht
    

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
                default_python_image="hailgenetics/hail:0.2.134-py3.11",
                backend=backend,
            )
        
        if args.run_vep:
            indel_ht = hl.read_table(f'{DATA_PATH}/vat/aou_indel_vat_v8.ht', _n_partitions=50000)
            # indel_ht = indel_ht.filter((indel_ht.locus.contig == 'chr21'))
            # indel_ht = hl.read_table(f'{DATA_PATH}/vat/aou_indel_vat_v8.ht')
            vep_ht = vep_or_lookup_vep(indel_ht, vep_version="vep105")
            vep_ht.describe()
            vep_ht.write(
                f'{DATA_PATH}/vep/aou_indel_vat_v8_vep.ht',
                overwrite=args.overwrite,
            )

        if args.check_path:
            j = b.new_python_job(name=f"Check aou path")
            j.call(check_path, path=args.path_to_check)

        if args.copy_files:
            j = b.new_python_job(name=f"Copy files")
            j.call(copy_files, args.original_path, args.target_path)
        
        if args.process_vat:
            j = b.new_python_job(name=f"Export VAT HT")
            j.call(process_vat, args.overwrite)

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
                if not hfs.exists(f'{DATA_PATH}/phenotype/raw/{name}.ht/_SUCCESS') or args.overwrite:
                    if args.batch:
                        j = b.new_python_job(name=f"Export raw pheno {name} HT")
                        j.call(write_raw_phenotype_ht, name, args.overwrite)
                    else:
                        write_raw_phenotype_ht(name, args.overwrite)
        
        if args.process_physical_measurements:
            if not hl.hadoop_exists(f'{DATA_PATH}/phenotype/physical_measurement.ht/_SUCCESS') or args.overwrite:
                if args.batch:
                    j = b.new_python_job(name=f"Process physical measurements")
                    j.call(process_quantitative_phenotypes, args.overwrite)
                else:
                    process_quantitative_phenotypes(args.overwrite)
        
        if args.summarize_raw_binary_phenotypes:
            for name in BINARY_CATEGORIES:
                job_depend_on = None
                if not hl.hadoop_exists(f'{DATA_PATH}/phenotype/raw/{name}.mt/_SUCCESS') or args.overwrite:
                    if args.batch:
                        j_mt = b.new_python_job(name=f"Export {name} MT")
                        j_mt.memory('highmem')
                        j_mt.cpu(8) 
                        j_mt.env('PYSPARK_SUBMIT_ARGS', '--driver-memory 32g --executor-memory 32g pyspark-shell')
                        j_mt.call(extract_mt_from_ht, name, args.overwrite) 
                        job_depend_on = j_mt
                    else:
                        extract_mt_from_ht(name, args.overwrite)
                if not hl.hadoop_exists(f'{DATA_PATH}/phenotype/summary/{name}_summary.txt.bgz') or args.overwrite:
                    if args.batch:
                        j = b.new_python_job(name=f"Summarize binary phenotypes {name}")
                        if job_depend_on is not None:
                            j.depends_on(job_depend_on)
                        j.call(summarize_raw_binary_phenotypes, name, args.overwrite)
                    else:
                        summarize_raw_binary_phenotypes(name, args.overwrite)

        if args.summarize_binary_phenotypes_per_ancestry:
            for name in PHENO_CATEGORIES:
                job_depend_on = None
                ext = '.mt' if name in BINARY_CATEGORIES else '.ht'
                if not hl.hadoop_exists(f'{DATA_PATH}/phenotype/{name}_annotated{ext}/_SUCCESS'):
                    if args.batch:
                        j_annotate = b.new_python_job(name=f"Annotate {name} {ext.upper()}")
                        j_annotate.call(annotate_phenotype_tables, name, overwrite=True)
                        job_depend_on = j_annotate
                    else:
                        annotate_phenotype_tables(name, overwrite=True)
                if name not in BINARY_CATEGORIES:
                    continue
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
        default=ORIGINAL_RANDOM_PATH
    )
    parser.add_argument(
        "--export-raw-pheno-ht",
        help="Export raw phenotype HT",
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
        "--summarize-binary-phenotypes-per-ancestry",
        help="Summarize binary phenotypes per ancestry",
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
    args = parser.parse_args()

    main(args)
