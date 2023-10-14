#!/usr/bin/env python3

__author__ = "wlu"

import hail as hl
import copy
import argparse
import logging
import hailtop.batch as hb
from hailtop.batch.resource import Resource, ResourceGroup
from hailtop.batch.job import Job
from hailtop.batch.batch import Batch
from collections import Counter
from shlex import quote as shq
import time

# from aou_gwas import *  # for dataproc
from utils.utils import *  # for QoB
from utils.resources import *  # for QoB

################################################ Remove after docker image is built ################################################
# ALL PATHs in the workbech: https://support.researchallofus.org/hc/en-us/articles/4616869437204-Controlled-CDR-Directory
# Descriptions: https://support.researchallofus.org/hc/en-us/articles/4614687617556-How-the-All-of-Us-Genomic-data-are-organized
import hail as hl

##################### ROOT PATHs #####################
TRANCHE = "250k"
BUCKET = "gs://aou_analysis"
MY_BUCKET = "gs://aou_wlu"
TMP_BUCKET = "gs://aou_tmp"
CURRENT_PATH = f"{BUCKET}/{TRANCHE}"
DATA_PATH = f"{CURRENT_PATH}/data"

##################### PARAMETERs #####################
GWAS_AF_CUTOFF = 0.001  # above
RVAS_AF_CUTOFF = 0.001  # below
CALLRATE_CUTOFF = 0.8
N_GENE_PER_GROUP = 10
MIN_CALL_RATE = {
    "all": 0.95,
    "afr": 0.9,
    "amr": 0.85,
    "eas": 0.9,
    "eur": 0.95,
    "mid": 0,
    "sas": 0.7,
}
SAIGE_DOCKER_IMAGE = "wzhou88/saige:1.3.0"
PILOT_PHENOTYPES = set(["height", "250.2", "411", "495", "585.3"])

##################### POPULATION SUMMARY (Computed in the AoU workbench):
# ANCESTRY_PATH = 'gs://fc-aou-datasets-controlled/v7/wgs/short_read/snpindel/aux/ancestry/ancestry_preds.tsv'
# pop_ht = hl.import_table(ANCESTRY_PATH, key="research_id", impute=True,
#                              types={"research_id": "tstr", "pca_features": hl.tarray(hl.tfloat)})
# pop_ht.aggregate(hl.agg.counter(pop_ht.ancestry_pred))
POPS = ("afr", "amr", "eas", "eur", "mid", "sas")
N_SAMPLES = {
    "all": 245394,
    "afr": 56913,
    "amr": 45035,
    "eas": 5706,
    "eur": 133581,
    "mid": 942,
    "sas": 3217,
}

##################### ORIGINAL PATHs #####################
ORIGINAL_GENO_ROOT = "gs://prod-drc-broad"
ORIGINAL_PHENO_ROOT = "gs://allxall-phenotypes/data"

ORIGINAL_PHENO_PATHs = [
    f"{ORIGINAL_PHENO_ROOT}/top_10_labs.csv",
    f"{ORIGINAL_PHENO_ROOT}/demographics_table.csv",
    f"{ORIGINAL_PHENO_ROOT}/physical_measurement_table.csv",
    f"{ORIGINAL_PHENO_ROOT}/r_drug_table.csv",
    f"{ORIGINAL_PHENO_ROOT}/pfhh_survey_table.csv",
    f"{ORIGINAL_PHENO_ROOT}/mcc2_phecode_table.csv",
]
TRAIT_TYPEs = ["physical_measurement", "r_drug", "pfhh_survey", "mcc2_phecode"]

ORIGINAL_CALLSET_PATH = f"{ORIGINAL_GENO_ROOT}/aou-wgs-delta-small_callsets_gq0/v7.1"
ACAF_MT_PATH = f"{ORIGINAL_CALLSET_PATH}/acaf_threshold_v7.1/splitMT/delta_basis_without_ext_aian_prod_gq0_3regions.acaf_threshold.split.mt"
EXOME_MT_PATH = f"{ORIGINAL_CALLSET_PATH}/exome_v7.1/splitMT/delta_basis_without_ext_aian_prod_gq0_3regions.exome.split.mt" # (34807589, 245394)
VDS_PATH = f"{ORIGINAL_GENO_ROOT}/v7/wgs/without_ext_aian_prod/vds/aou_srwgs_short_variants_v7_without_ext_aian_prod.vds"

ORIGINAL_UTIL_PATH = f"{ORIGINAL_GENO_ROOT}/aou-wgs-delta-aux_gq0"
ORIGINAL_POP_PATH = (
    f"{ORIGINAL_UTIL_PATH}/ancestry/delta_v1_gt_no_ext_aian_gq0_prod.ancestry_preds.tsv"
)
ORIGINAL_RELATED_INFO_PATH = (
    f"{ORIGINAL_UTIL_PATH}/relatedness/delta_v1_gt_no_ext_aian_gq0_prod.relatedness.tsv"
)
ORIGINAL_RELATED_SAMPLE_PATH = f"{ORIGINAL_UTIL_PATH}/relatedness/delta_v1_gt_no_ext_aian_gq0_prod.relatedness_flagged_samples.tsv"
SAMPLE_INFO_PATHs = [
    ORIGINAL_POP_PATH,
    ORIGINAL_RELATED_INFO_PATH,
    ORIGINAL_RELATED_SAMPLE_PATH,
]
table_keys = {
    "ancestry_preds": ["research_id"],
    "relatedness": ["i.s", "j.s"],
    "relatedness_flagged_samples": ["sample_id"],
}


##################### PHENOTYPE & SAMPLE PATHs #####################
def get_raw_phenotype_path(
    name: str, parsed: bool = True, extension: str = "ht", tranche: str = TRANCHE
):
    assert name in [
        "r_drug_table",
        "pfhh_survey_table",
        "mcc2_phecode_table",
        "demographics_table",
        "physical_measurement_table",
        "top_10_labs",
    ], 'Name has to be from ["r_drug_table", "pfhh_survey_table", "mcc2_phecode_table", "demographics_table", "physical_measurement_table", "top_10_labs"]'
    return (
        f'{DATA_PATH}/phenotype/{"" if parsed else "raw/"}{name}_{tranche}.{extension}'
    )


def get_random_phenotype_path(pop: str, extension: str = "tsv", test: bool = False):
    return f'{DATA_PATH}/phenotype/random_phenos/random_pheno_{pop}{"_test" if test else ""}.{extension}'


BINARY_HT_PATHs = [
    get_raw_phenotype_path(name="r_drug_table", parsed=True),
    get_raw_phenotype_path(name="pfhh_survey_table", parsed=True),
    get_raw_phenotype_path(name="mcc2_phecode_table", parsed=True),
]
DEMOGRAPHICS_HT_PATH = get_raw_phenotype_path(name="demographics_table", parsed=True)
QUANTITATIVE_HT_PATH = get_raw_phenotype_path(
    name="physical_measurement_table", parsed=True
)


def get_full_phenotype_path(
    annotation: bool,
    random_pheno: bool = False,
    extension: str = "ht",
    tranche: str = TRANCHE,
):
    return f'{DATA_PATH}/phenotype/aou_all{"_random" if random_pheno else ""}_phenotypes{"_annotated" if annotation else ""}_{tranche}.{extension}'


def get_sample_meta_path(
    annotation: bool, extension: str = "ht", tranche: str = TRANCHE
):
    return f'{DATA_PATH}/utils/aou_sample_meta_info{"_annotated" if annotation else ""}_{tranche}.{extension}'


def get_aou_util_path(
    name: str, parsed: bool = False, extension: str = "ht", tranche: str = TRANCHE
):
    assert name in [
        "ancestry_preds",
        "relatedness",
        "relatedness_flagged_samples",
        "vat",
        "mt_sample_qc",
        "acaf_mt_sample_qc",
        "vds_sample_qc",
        "sample_qc",
        "variant_qc"
        "sample_qc_tmp",
        "vep_corrected",
        "vep_old",
    ], 'Name has to be from ["ancestry_preds", "relatedness", "relatedness_flagged_samples", "vat", "mt_sample_qc", "vds_sample_qc","sample_qc", "variant_qc", "sample_qc_tmp", "vep_corrected", "vep_old"]'
    tag1 = (name == "ancestry_preds") and parsed
    tag2 = (name == "gene_map") and parsed
    if name == "vat":
        # gsutil -m cp -r gs://aou_wlu/utils/vat/hail_0_2_107/aou_PARSED_SORTED_vat.ht gs://aou_analysis/250k/data/utils/aou_parsed_and_sorted_vat_hail_0_2_107_250k.ht
        name = "parsed_and_sorted_vat_hail_0_2_107"
    return f'{DATA_PATH}/utils/aou_{name}{"_parsed" if tag1 else ""}{"_processed" if tag2 else ""}_{tranche}.{extension}'


def get_aou_relatedness_path(extension: str = "ht", tranche: str = TRANCHE):
    return f"{DATA_PATH}/utils/aou_ibd_relatedness_{tranche}.{extension}"

def get_aou_gene_map_ht_path(pop: str, processed=False, tranche: str = TRANCHE):
    return f"{DATA_PATH}/utils/gene_map/aou_{pop}_gene_map{'_processed' if processed else ''}_{tranche}.ht"


##################### Random Phenotype GRM PATHs #####################
def get_aou_sites_for_grm_path(
    pop: str, extension: str, pruned: bool = False, tranche: str = TRANCHE
):
    return f'{DATA_PATH}/utils/grm/aou_{pop}_sites{"_pruned" if pruned else ""}_for_grm_{tranche}.{extension}'

def get_aou_sample_file_path(
    pop: str, tranche: str = TRANCHE
):
    return f'{DATA_PATH}/utils/grm/aou_{pop}_{tranche}.samples'



##################### Result PATHs #####################

GENE_INTERVAL_PATH = (
    f"gs://aou_wlu/data/group_positions_{N_GENE_PER_GROUP}_protein_coding.ht"
)

VARIANT_INTERVAL_PATH = f"gs://aou_wlu/data/variant_intervals.ht"


def get_aou_saige_results_root(analysis_type: str, pop: str, name: str, test: bool = False):
    assert analysis_type in [
        "gene",
        "variant",
    ], "Name has to be from ['gene', 'variant']"
    assert name in [
        "bgen",
        "result",
    ], "Name has to be from ['bgen', 'result']"
    return f'{CURRENT_PATH}/{analysis_type}_results/{name}/{"test" if test else pop.upper()}'


def get_aou_saige_utils_root(pop: str, name: str, random_pheno: bool = False):
    assert name in [
        "pheno_file",
        "null_glmm",
    ], "Name has to be from ['pheno_file', 'null_glmm']"
    return f"{CURRENT_PATH}/{name}/{pop.upper()}{'/random_pheno' if random_pheno else ''}"

#################
def get_filtered_mt(analysis_type: str, pop: str, filter_samples: bool=True, filter_variants: bool=True, adj_filter: bool=True):
    from gnomad.utils.annotations import annotate_adj
    mt_path = EXOME_MT_PATH if analysis_type=='gene' else ACAF_MT_PATH
    mt = hl.read_matrix_table(mt_path)
    meta_ht = hl.read_table(get_sample_meta_path(annotation=True))
    mt = mt.annotate_cols(**meta_ht[mt.col_key])
    mt = mt.filter_entries(mt.FT != 'FAIL') # TODO: if using VDS, that's a bool
    # print(mt.aggregate_entries(hl.agg.counter(mt.FT)))
    # {'FAIL': 2323118301, 'PASS': 63243269454, None: 8308706765322}
    if pop is not None:
        mt = mt.filter_cols(mt.pop==pop)
    if filter_variants:
        mt = mt.filter_rows(
            (hl.len(mt.filters) == 0) & # TODO: update this with the VDS filters
            (mt.info.AC[0] > 0)
        )
    if filter_samples:
        mt = mt.filter_cols(mt.samples_to_keep)
    if adj_filter:
        """
        Filter genotypes to adj criteria - Default: 
        GQ >= 20, haploid_DP >= 5, else_DP >= 10, 
        het_ref_altAB >= 0.2, het_non_ref_altAB >= 0.2 , het_non_ref_refAB > =0.2
        """
        mt= mt.annotate_entries(DP=hl.sum(mt.AD))
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

# HAIL_DOCKER_IMAGE = 'gcr.io/ukbb-diversepops-neale/hail_utils:6.1'
HAIL_DOCKER_IMAGE = "hailgenetics/hail:0.2.124-py3.9"
# SAIGE_DOCKER_IMAGE = 'gcr.io/ukbb-diversepops-neale/saige:0.5'
SAIGE_DOCKER_IMAGE = "wzhou88/saige:1.3.0"  # latest
QQ_DOCKER_IMAGE = "konradjk/saige_qq:0.2"

MKL_OFF = "export MKL_NUM_THREADS=1; export MKL_DYNAMIC=false; export OMP_NUM_THREADS=1; export OMP_DYNAMIC=false; "

def run_saige(
    p: Batch,
    output_root: str,
    model_file: str,
    variance_ratio_file: str,
    sparse_grm_file: str,
    bgen_file: ResourceGroup,
    samples_file: ResourceGroup,
    docker_image: str,
    group_file: str = None,
    trait_type: str = "continuous",
    chrom: str = "chr1",
    min_mac: int = 1,
    min_maf: float = 0,
    max_maf_for_group: float = 0.5,
    memory: str = "",
    storage: str = "10Gi",
    add_suffix: str = "",
):
    """
    Change log:
    use_bgen and log_pvalue are defaults now
    - Removed --IsOutputlogPforSingle, --sparseSigmaFile
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
    analysis_type = "gene" if group_file is not None else "variant"
    run_saige_task: Job = (
        p.new_job(name=f"run_saige", attributes={"analysis_type": analysis_type})
        .cpu(1)
        .storage(storage)
        .image(docker_image)
        .memory(memory)
    )  # Step 2 is single-threaded only

    if analysis_type == "gene":
        run_saige_task.declare_resource_group(
            result={
                f"{add_suffix}gene.txt": "{root}",
                f"{add_suffix}single.txt": "{root}_single",
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
        f"--maxMAF_in_groupTest={max_maf_for_group} "
        f"--sampleFile={samples_file} " 
        f"--GMMATmodelFile={model_file} " 
        f"--varianceRatioFile={variance_ratio_file} " 
        f"--sparseGRMFile={sparse_grm_file} "
        f"--sparseGRMSampleIDFile={sparse_grm_file}.sampleIDs.txt " 
        f"--LOCO=FALSE " 
        f"--AlleleOrder=ref-first " 
        f"--SAIGEOutputFile={run_saige_task.result} "
        f"--is_fastTest=TRUE"
    )

    if analysis_type == "gene":
        if trait_type == "binary":
            command += (
                f"--IsOutputPvalueNAinGroupTestforBinary=TRUE " #
                f"--is_output_moreDetails=TRUE"
            )
        command += (
            f"--groupFile={group_file} " 
            f"--is_output_markerList_in_groupTest=TRUE "
            f"--is_single_in_groupTest=TRUE "
            f"--IsOutputBETASEinBurdenTest=TRUE " #
        )
    command += f"--IsOutputAFinCaseCtrl=TRUE 2>&1 | tee {run_saige_task.stdout}; " #
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
    p.write_output(run_saige_task.result, output_root)
    p.write_output(run_saige_task.stdout, f"{output_root}.{analysis_type}.log")
    return run_saige_task


def fit_null_glmm(
    p: Batch,
    output_root: str,
    phenoname: str,
    pheno_file: Resource,
    trait_type: str,
    covariates: str,
    plink_file_root: str,
    docker_image: str,
    sparse_grm: Resource = None,
    sparse_grm_extension: str = None,
    inv_normalize: bool = False,
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

    pheno_col = "value"
    user_id_col = "person_id"
    in_bfile = p.read_input_group(
        **{ext: f"{plink_file_root}.{ext}" for ext in ("bed", "bim", "fam")}
    )
    fit_null_task = (
        p.new_job(
            name=f"fit_null_model_{phenoname}",
            attributes={"trait_type": trait_type},
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
    # if trait_type == 'icd':
    #     bim_fix_command += (f"; zcat {pheno_file.gz} | perl -p -e 's/true/1/g' | perl -p -e 's/false/0/g' "
    #                         f"| gzip -c > {pheno_file.gz}.temp.gz; mv {pheno_file.gz}.temp.gz {pheno_file.gz}")

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
    pop,
    chromosome,
    analysis_type,
    interval,
    output_dir,
    mean_impute_missing,
    no_adj: bool=False,
    variant_af_filter: float = 0.001,
    variant_ac_filter: int = 0,
    variant_callrate_filter: float=0,
):
    mt = get_filtered_mt(analysis_type=analysis_type, filter_variants=True, filter_samples=True, adj_filter= not no_adj, pop=pop)

    outname = f"{analysis_type}_{chromosome}_{str(interval.start.position).zfill(9)}_{interval.end.position}"

    # Filter to interval
    mt = hl.filter_intervals(mt, [interval])

    mt = mt.select_entries("GT")
    mt = mt.filter_rows(
        hl.agg.count_where(mt.GT.is_non_ref()) > 0
    )  # Filter to non-reference sites
    mt = mt.annotate_rows(
        rsid=mt.locus.contig + ":" + hl.str(mt.locus.position) + "_" + mt.alleles[0] + "/" + mt.alleles[1]
    )  # Annotate rsid

    if variant_ac_filter:
        mt = mt.filter_rows(mt.info.AC[0] > variant_ac_filter)

    if variant_callrate_filter:
        mt = mt.filter_rows(mt.info.AN/(2*N_SAMPLES[pop]) > variant_callrate_filter)

    mt = mt.annotate_entries(
        GT=hl.if_else(mt.GT.is_haploid(), hl.call(mt.GT[0], mt.GT[0]), mt.GT)
    )
    mt = gt_to_gp(mt)
    mt = impute_missing_gp(mt, mean_impute=mean_impute_missing)
    hl.export_bgen(mt, f"{output_dir}/{outname}", gp=mt.GP, varid=mt.rsid)


def export_gene_group_file(interval, groups, pop, output_dir):
    gene_ht = hl.read_table(get_aou_gene_map_ht_path(pop=pop, processed=True))
    gene_ht = hl.filter_intervals(gene_ht, [interval])

    gene_ht = gene_ht.filter(hl.set(groups).contains(gene_ht.annotation))
    gene_ht.select(
        group=gene_ht.gene_id + "_" + gene_ht.gene_symbol + "_" + gene_ht.annotation
        + hl.if_else(gene_ht.common_variant, "_" + gene_ht.variants[0], ""),
        variant=hl.delimit(gene_ht.variants, "\t"),
    ).key_by().drop("start").export(output_dir, header=False)


def main(args):
    hl.init(
        tmp_dir=TMP_BUCKET,
        driver_memory="highmem",
        driver_cores=8,
        worker_memory="highmem",
        worker_cores=1,
        default_reference="GRCh38",
        log="/aou_saige.log",
        app_name=f'aou_SAIGE{"_export_phenos" if args.export_phenos else ""}{"_random" if args.random_pheno else ""}'
    )

    num_pcs = 16
    start_time = time.time()
    basic_covars = [
        "sex",
        "age",
        "age2",
        "age_sex",
        "age2_sex",
    ]  # TODO: add batch covariates for future release
    covariates = ",".join(basic_covars + [f"PC{x}" for x in range(1, num_pcs + 1)])
    n_threads = 8
    chromosomes = list(map(str, range(1, 23))) + ["X", "Y"]
    reference = "GRCh38"
    chrom_lengths = hl.get_reference(reference).lengths
    pops = args.pops.split(",") if args.pops else POPS

    pheno_ht = hl.read_table(get_full_phenotype_path(annotation=True, random_pheno=args.random_pheno))
    pheno_ht = pheno_ht.filter(pheno_ht.samples_to_keep)

    if args.pilot or (args.phenos is not None):
        phenos_to_run = PILOT_PHENOTYPES if args.pilot else args.phenos.split(",")
        quantitative_phenos = hl.eval(pheno_ht.trait_type)["continuous"]
    elif args.random_pheno:
        raw_random_pheno_ht = hl.import_table(
            get_random_phenotype_path(pop=pops[0], test=True),
            key="userId",
            impute=True,
            types={"userId": hl.tstr},
        )
        phenos_to_run = list(raw_random_pheno_ht.row_value)
        quantitative_phenos = [pheno for pheno in phenos_to_run if 'continuous' in pheno ]
    else:
        raw_pheno_ht = hl.read_table(
            get_full_phenotype_path(annotation=False, random_pheno=args.random_pheno)
        )
        phenos_to_run = list(raw_pheno_ht.row_value)
        quantitative_phenos = hl.eval(pheno_ht.trait_type)["continuous"]
        # TODO: create phenotype info table and add phenotype filters (min cases, etc)

    print(f"Got {len(phenos_to_run)} phenotypes...")
    if len(phenos_to_run) <= 20:
        print(phenos_to_run)

    if args.export_phenos:
        for pop in pops:
            ht = pheno_ht
            if pop != "all":
                ht = ht.filter(ht.pop == pop)
                print(f"N samples from {pop.upper()}: {ht.count()}...")

            for phenoname in phenos_to_run:
                pheno_path = f'{get_aou_saige_utils_root(pop=pop, name="pheno_file", random_pheno=args.random_pheno)}/phenotype_{phenoname}.tsv'
                if (not hl.hadoop_exists(pheno_path)
                        or args.overwrite_pheno_data
                ):
                    print(f"Exporting {phenoname} for {pop.upper()}...")
                    if (phenoname in quantitative_phenos):
                        ht = ht.annotate(value=ht[phenoname])
                    else:
                        ht = ht.annotate(value=hl.int(ht[phenoname]))
                        print(f"N cases for {pop.upper()} samples {phenoname} : {ht.aggregate(hl.agg.sum(ht.value))}...")
                    fields = covariates.split(",") + ["related", "pop", "value"]
                    out_ht = ht.select(**{field: ht[field] for field in fields})
                    print(f"N samples defined for {phenoname}: {out_ht.aggregate(hl.agg.count_where(hl.is_defined(out_ht.value)))}...")

                    if args.proportion_single_sex > 0:
                        n_cases_female = out_ht.aggregate(
                            hl.agg.filter(
                                out_ht.sex == "Female",
                                hl.agg.count_where(out_ht.values == 1),
                            )
                        )
                        n_cases_male = out_ht.aggregate(
                            hl.agg.filter(
                                out_ht.sex == "Male",
                                hl.agg.count_where(out_ht.values == 1),
                            )
                        )
                        prop_female = n_cases_female / (n_cases_male + n_cases_female)
                        print(f"Female proportion: {prop_female}")
                        if prop_female <= args.proportion_single_sex:
                            print(
                                f"Female case proportion {prop_female} less than {args.proportion_single_sex}. Filtering to males..."
                            )
                            out_ht = out_ht.filter_rows(out_ht.sex == 1)
                        elif prop_female >= 1 - args.proportion_single_sex:
                            print(
                                f"Female case proportion {prop_female} greater than {1 - args.proportion_single_sex}. Filtering to females..."
                            )
                            out_ht = out_ht.filter_rows(out_ht.sex == 0)
                    out_ht.export(pheno_path)

    if args.run_saige:
        analysis_type = "variant" if args.single_variant_only else "gene"
        print(f'Analysis type: {analysis_type}')

        # Obtain intervals for SAIGE
        if analysis_type == "gene":
            if args.groups is None:
                groups = ["pLoF", "missense|LC", "synonymous"]
            else:
                groups = args.groups.split(",")
            interval_ht = hl.read_table(GENE_INTERVAL_PATH)
        else:
            if not hl.hadoop_exists(VARIANT_INTERVAL_PATH):
                intervals = []
                chunk_size = int(5e5)
                for chrom in chromosomes:
                    chromosome = f"chr{chrom}"
                    chrom_length = chrom_lengths[chromosome]
                    for start_pos in range(1, chrom_length, chunk_size):
                        end_pos = (
                            chrom_length
                            if start_pos + chunk_size > chrom_length
                            else (start_pos + chunk_size)
                        )
                        interval = hl.eval(
                            hl.parse_locus_interval(
                                f"[{chromosome}:{start_pos}-{end_pos}]",
                                reference_genome="GRCh38",
                            )
                        )
                        intervals.append(interval)
                import pandas as pd
                df = pd.DataFrame({'interval': intervals})
                ht = hl.Table.from_pandas(df)
                ht.write(VARIANT_INTERVAL_PATH)
            else:
                interval_ht = hl.read_table(VARIANT_INTERVAL_PATH)

        backend = hb.ServiceBackend(
            billing_project="all-by-aou", remote_tmpdir=TMP_BUCKET
        )
        print(f"Docker image: {SAIGE_DOCKER_IMAGE}...")
        for pop in pops:
            b = hb.Batch(
                name=f"saige_aou_{pop}",
                backend=backend,
                default_image=SAIGE_DOCKER_IMAGE,
                default_storage="500Mi",
                default_cpu=n_threads,
            )

            # Obtain Sparse GRMs
            relatedness_cutoff = "0.125"
            num_markers = 2000
            n_threads = 8
            sparse_grm_root = f"{DATA_PATH}/utils/grm/aou_{pop}"
            sparse_grm_extension = f"_relatednessCutoff_{relatedness_cutoff}_{num_markers}_randomMarkersUsed.sparseGRM.mtx"
            sparse_grm = b.read_input_group(
                **{ext: f"{sparse_grm_root}.{ext}"
                   for ext in (sparse_grm_extension, f"{sparse_grm_extension}.sampleIDs.txt",)}
            )

            overwrite_null_models = args.create_null_models
            null_model_dir = get_aou_saige_utils_root(pop=pop, name="null_glmm", random_pheno=args.random_pheno)
            null_models_already_created = {}
            if not overwrite_null_models and hl.hadoop_exists(null_model_dir):
                null_models_already_created = {
                    x["path"] for x in hl.hadoop_ls(null_model_dir)
                }
            print(f'Found {int(len(null_models_already_created) / 3)} Null models in directory...')
            null_models = {}
            pheno_exports = {}
            for phenoname in phenos_to_run:
                trait_type = ("quantitative"if (phenoname in quantitative_phenos) or ("coutinuous" in phenoname) else "binary")
                null_glmm_root = f"{null_model_dir}/phenotype_{phenoname}"
                model_file_path = f"{null_glmm_root}.rda"
                variance_ratio_file_path =  f"{null_glmm_root}.varianceRatio.txt"

                pheno_file = b.read_input(f'{get_aou_saige_utils_root(pop=pop, name="pheno_file", random_pheno=args.random_pheno)}/phenotype_{phenoname}.tsv')
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
                    fit_null_task = fit_null_glmm(
                        b,
                        output_root=null_glmm_root,
                        phenoname=phenoname,
                        pheno_file=pheno_exports[phenoname],
                        trait_type=trait_type,
                        covariates=covariates,
                        plink_file_root=get_aou_sites_for_grm_path(pop=pop, extension="plink", pruned=True),
                        docker_image=SAIGE_DOCKER_IMAGE,
                        sparse_grm=sparse_grm,
                        sparse_grm_extension=sparse_grm_extension,
                        inv_normalize=False,
                        skip_model_fitting=False,
                        n_threads=n_threads,
                        storage="100Gi",
                        memory="60Gi",
                        non_pre_emptible=False
                    )
                    fit_null_task.attributes.update({"pop": pop})
                    fit_null_task.attributes.update(
                        {"phenotype": copy.deepcopy(phenoname)}
                    )
                    model_file = fit_null_task.null_glmm.rda
                    variance_ratio_file = fit_null_task.null_glmm[f"varianceRatio.txt"]

                null_models[phenoname] = (model_file, variance_ratio_file)

            print(f"------------{pop.upper()} MT for {analysis_type} analysis bgen files------------")
            bgen_dir = get_aou_saige_results_root(analysis_type=analysis_type, pop=pop, name="bgen")
            overwrite_bgens = args.create_bgens
            bgens_already_created = {}
            if not overwrite_bgens and hl.hadoop_exists(bgen_dir):
                bgens_already_created = {x["path"] for x in hl.hadoop_ls(bgen_dir)}
            print(f'Found {int(len(bgens_already_created)/2)} Bgens in directory...')

            bgens = {}
            for chrom in chromosomes:
                chromosome = f"chr{chrom}"
                sub_interval_ht = interval_ht.filter(interval_ht.interval.start.contig == chromosome)
                intervals = sub_interval_ht.aggregate(
                    hl.agg.collect(sub_interval_ht.interval)
                )
                for interval in intervals:
                    bgen_root = f"{bgen_dir}/{analysis_type}_{interval.start.contig}_{str(interval.start.position).zfill(9)}_{interval.end.position}"
                    if f"{bgen_root}.bgen" not in bgens_already_created:
                        bgen_task = b.new_python_job(
                            name=f"{analysis_type}_analysis_export_{interval}_bgen"
                        )
                        if args.hail_image is None:
                            image = hb.build_python_image(
                                "us-central1-docker.pkg.dev/aou-neale-gwas/hail-batch/hail_gnomad_python_3.9",
                                requirements=["hail", "gnomad"],
                                python_version="3.9",
                            )
                        else:
                            image = args.hail_image
                        bgen_task.image(image)
                        bgen_task.call(
                            export_bgen_from_mt,
                            pop=pop,
                            chromosome=chromosome,
                            analysis_type=analysis_type,
                            interval=interval,
                            output_dir=bgen_dir,
                            mean_impute_missing=args.mean_impute_missing,
                            no_adj=args.no_adj,
                            variant_ac_filter= 0,
                            variant_callrate_filter=0
                        )
                        bgen_task.attributes["pop"] = pop
                        bgen_task.attributes["chromosome"] = chromosome
                        bgen_task.attributes["analysis_type"] = analysis_type
                    if (f"{bgen_root}.gene.txt" not in bgens_already_created) and analysis_type=='gene':
                        gene_txt_task = b.new_python_job(
                            name=f"{analysis_type}_analysis_export_{hl.eval(hl.str(interval))}_gene_txt"
                        )
                        gene_txt_task.image(HAIL_DOCKER_IMAGE)
                        gene_txt_task.call(
                            export_gene_group_file,
                            interval=interval,
                            groups=groups,
                            pop=pop,
                            output_dir=f"{bgen_root}.gene.txt",
                        )
                        gene_txt_task.attributes["pop"] = pop
                        gene_txt_task.attributes["chromosome"] = chromosome

                    bgen_file = b.read_input_group(
                        **{
                            "bgen": f"{bgen_root}.bgen",
                            "bgen.bgi": f"{bgen_root}.bgen.bgi",
                            "sample": f"{bgen_root}.sample",
                        }
                    )
                    if analysis_type == 'gene':
                        group_file = b.read_input(f"{bgen_root}.gene.txt")
                        bgens[hl.eval(hl.str(interval))] = (bgen_file, group_file)
                    else:
                        bgens[hl.eval(hl.str(interval))] = bgen_file
                    if args.test:
                        break
                if args.test:
                    break

            result_dir = get_aou_saige_results_root(
                analysis_type=analysis_type, pop=pop, name="result"
            )
            overwrite_results = args.overwrite_results
            for i in range(len(phenos_to_run)):
                if args.skip_saige:
                    break
                phenoname = list(phenos_to_run)[i]
                if phenoname not in null_models.keys():
                    continue

                model_file, variance_ratio_file = null_models[phenoname]

                if not i % 10:
                    n_jobs = dict(Counter(map(lambda x: x.name, b.select_jobs("")))).get("run_saige", 0)
                    logger.info(f"Read {i} phenotypes ({n_jobs} new to run so far)...")

                pheno_results_dir = f"{result_dir}/phenotype_{phenoname}"
                results_already_created = {}

                if (
                    not overwrite_results
                    and not args.skip_saige
                    and hl.hadoop_exists(pheno_results_dir)
                ):
                    results_already_created = {x["path"] for x in hl.hadoop_ls(pheno_results_dir)}

                saige_tasks = []
                for chrom in chromosomes:
                    chromosome = f"chr{chrom}"
                    sub_interval_ht = interval_ht.filter(interval_ht.interval.start.contig == chromosome)
                    intervals = sub_interval_ht.aggregate(
                        hl.agg.collect(sub_interval_ht.interval)
                    )
                    for interval in intervals:
                        results_path = f"{pheno_results_dir}/result_{phenoname}_{interval.start.contig}_{str(interval.start.position).zfill(9)}"
                        if analysis_type == "gene":
                            bgen_file, group_file = bgens[hl.eval(hl.str(interval))]
                            max_maf_for_group_test = RVAS_AF_CUTOFF
                        else:
                            bgen_file = bgens[hl.eval(hl.str(interval))]
                            group_file = None
                            max_maf_for_group_test = None

                        if (
                            overwrite_results
                            or f"{results_path}.single_variant.txt"
                            not in results_already_created
                        ):
                            samples_file = b.read_input(get_aou_sample_file_path(pop=pop))
                            saige_task = run_saige(
                                p=b,
                                output_root=results_path,
                                model_file=model_file,
                                variance_ratio_file=variance_ratio_file,
                                sparse_grm_file=sparse_grm[sparse_grm_extension],
                                bgen_file=bgen_file,
                                samples_file=samples_file,
                                docker_image=SAIGE_DOCKER_IMAGE,
                                group_file=group_file,
                                trait_type=trait_type,
                                chrom=interval.start.contig,
                                min_mac=1,
                                min_maf=0,
                                max_maf_for_group=max_maf_for_group_test,
                            )
                            saige_task.attributes.update(
                                {"interval": hl.eval(hl.str(interval)), "pop": pop, "chromosome":chromosome}
                            )
                            saige_task.attributes.update(
                                {"phenotype": copy.deepcopy(phenoname)}
                            )
                            saige_tasks.append(saige_task)
                        if args.test:
                            break
                    if args.test:
                        break
            b.run()
        #     res_tasks = []
        #     if (
        #         overwrite_results
        #         or args.overwrite_hail_results
        #         or f"{pheno_results_dir}/variant_results.ht"
        #         not in results_already_created
        #         or not hl.hadoop_exists(
        #             f"{pheno_results_dir}/variant_results.ht/_SUCCESS"
        #         )
        #     ):
        #         null_glmm_root = (
        #             f"{null_model_dir}/phenotype_{phenoname}.{analysis_type}.log"
        #         )
        #
        #         prefix = f"{pheno_results_dir}/result_{phenoname}_chr{{chrom}}_{str(1).zfill(9)}"
        #         saige_log = f"{prefix}.{analysis_type}.log"
        #
        #         load_task = load_results_into_hail(
        #             b,
        #             pheno_results_dir,
        #             phenoname,
        #             saige_tasks,
        #             get_aou_util_path("vep_corrected"),
        #             HAIL_DOCKER_IMAGE,
        #             saige_log=saige_log,
        #             analysis_type=analysis_type,
        #             n_threads=n_threads,
        #             null_glmm_log=null_glmm_root,
        #             reference=reference,
        #             legacy_annotations=True,
        #             log_pvalue=True,
        #         )
        #         load_task.attributes["pop"] = pop
        #         res_tasks.append(load_task)
        #         qq_export, qq_plot = qq_plot_results(
        #             b,
        #             pheno_results_dir,
        #             res_tasks,
        #             HAIL_DOCKER_IMAGE,
        #             QQ_DOCKER_IMAGE,
        #             n_threads=n_threads,
        #         )
        #         qq_export.attributes.update({"pop": pop})
        #         qq_export.attributes.update({"phenotype": copy.deepcopy(phenoname)})
        #         qq_plot.attributes.update({"pop": pop})
        #         qq_plot.attributes.update({"phenotype": copy.deepcopy(phenoname)})
        #
        # def get_tasks_from_pipeline(p):
        #     return dict(Counter(map(lambda x: x.name, p.select_jobs(""))))
        #
        # logger.info(
        #     f'Setup took: {time.strftime("%H:%M:%S", time.gmtime(time.time() - start_time))}'
        # )
        # logger.info(f"Submitting: {get_tasks_from_pipeline(b)}")
        # logger.info(f"Total size: {sum([len(x._pretty()) for x in b.select_jobs('')])}")
        # logger.info(f"Finished: {get_tasks_from_pipeline(b)}")
        # b.run()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--single_variant_only", help="Run single variant SAIGE", action="store_true"
    )
    parser.add_argument("--run_saige", help="Run SAIGE pipeline", action="store_true")
    parser.add_argument(
        "--overwrite_pheno_data", help="Overwrite phenotype data", action="store_true"
    )
    # parser.add_argument('--sex_stratified', help='Run these phenotypes in a sex-stratified fashion (experimental)', choices=(None, 'all', 'only'))
    parser.add_argument(
        "--skip_any_null_models",
        help="Skip running SAIGE null models",
        action="store_true",
    )
    parser.add_argument(
        "--skip_saige", help="Skip running SAIGE tests", action="store_true"
    )
    parser.add_argument(
        "--export_phenos", help="Export phenotype tsv files", action="store_true"
    )
    parser.add_argument(
        "--create_null_models",
        help="Force creation of null models",
        action="store_true",
    )
    parser.add_argument(
        "--create_bgens", help="Force creation of Bgen files", action="store_true"
    )
    parser.add_argument(
        "--overwrite_results", help="Force run of SAIGE tests", action="store_true"
    )
    parser.add_argument(
        "--overwrite_hail_results",
        help="Force run of results loading",
        action="store_true",
    )
    parser.add_argument(
        "--test",
        help="Test run of pipeline, using chromosome 1 only",
        action="store_true",
    )
    parser.add_argument(
        "--non_pre_emptible", help="Local test of pipeline", action="store_true"
    )
    # parser.add_argument('--skip_case_count_filter', help='Skip running SAIGE tests', action='store_true')
    parser.add_argument(
        "--phenos",
        help="Comma-separated list of trait_type-phenocode-pheno_sex-coding-modifier regexes "
        "(e.g. continuous-50-both_sexes--,icd10-E1.*,brain_mri-.* )",
    )
    parser.add_argument(
        "--pops", help="comma-separated list", default="afr,amr,eas,eur,mid,sas"
    )
    parser.add_argument(
        "--groups", help="comma-separated list of gene groups to run RVAS"
    )
    parser.add_argument(
        "--pilot", help="Run pilot phenotypes only", action="store_true"
    )
    parser.add_argument(
        "--random_pheno", help="Use random phenotypes", action="store_true"
    )
    parser.add_argument(
        "-hail_image", help="hail docker image to use", action="store_true", default='us-central1-docker.pkg.dev/broad-mpg-gnomad/images/vrs084'
    )
    parser.add_argument(
        "--proportion_single_sex",
        help="If set and proportion of male or female cases is less than "
        "this number, then filter to females and males respectively",
        type=float,
        default=0.0,
    )
    parser.add_argument(
        "--mean_impute_missing",
        help="Whether to mean impute missing genotypes (BGEN only) "
        "(default: set to hom ref)",
        action="store_true",
    )
    parser.add_argument(
        "--no_adj",
        help="Use all genotypes instead of only high-quality ones",
        action="store_true",
    )
    parser.add_argument(
        "--callrate_filter",
        help="Impose filter of specified callrate (default: none)",
        default=0.0,
        type=float,
    )
    args = parser.parse_args()

    main(args)
