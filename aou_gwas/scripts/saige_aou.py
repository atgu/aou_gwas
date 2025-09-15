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
from tqdm import tqdm

# from aou_gwas import *  # for dataproc
from aou_gwas.utils.utils import *  # for QoB
from aou_gwas.utils.resources import *  # for QoB
from aou_gwas.utils.results_loading import * # for QoB

def range_table(log_file):
    import hail as hl
    hl.init(log=log_file)
    print(hl.utils.range_table(10)._force_count())

################################################ Copy-paste resources.py - Remove after docker image is built ################################################


################################################ Remove after docker image is built ################################################

logging.basicConfig(
    format="%(levelname)s (%(name)s %(lineno)s): %(message)s",
    level="INFO",
    filename="saige_pipeline.log",
)
logger = logging.getLogger("ALL_x_AoU_SAIGE")
logger.setLevel(logging.INFO)

HAIL_DOCKER_IMAGE = "hailgenetics/hail:0.2.124-py3.9"
SAIGE_DOCKER_IMAGE = "wzhou88/saige:1.3.0"  # latest
QQ_DOCKER_IMAGE = "konradjk/saige_qq:0.2"


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
                      saige_log: str = 'NA',
                      extension: str = 'single.txt',
                      overwrite: bool = False,
                      variant_type: str = 'genome',
                      num_partitions: int = 1000):
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
    output_ht_path = f'{output_ht_directory}/variant_results.ht'
    ht = hl.import_table(f'{directory}/*.{extension}', delimiter='\t', impute=True, types = {'p.value':hl.tstr})
    print(f'Loading: {directory}/*.{extension} ...')
    marker_id_col = 'MarkerID'
    # marker_id_col = 'MarkerID' if extension == 'single_variant.txt' else 'SNPID'
    locus_alleles = ht[marker_id_col].split('_')
    if n_cases == -1: n_cases = hl.null(hl.tint)
    if n_controls == -1: n_controls = hl.null(hl.tint)
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
        ht = ht.drop('N')

    ht = ht.key_by(locus=hl.locus(ht.CHR, ht.POS, reference_genome = 'GRCh38'), alleles=[ht.Allele1, ht.Allele2],
                   phenoname=phenoname).distinct().naive_coalesce(num_partitions)
    if marker_id_col == 'SNPID':
        ht = ht.drop('CHR', 'POS', 'SNPID', 'Allele1', 'Allele2')
    ht = ht.transmute(Pvalue=ht['p.value']).annotate_globals(
        n_cases=n_cases, n_controls=n_controls, heritability=heritability, saige_version=saige_version,
        inv_normalized=inv_normalized, log_pvalue=True)
    ht = ht.annotate(Pvalue_log10 = parse_log_p_value(ht.Pvalue))
    ht = ht.annotate(Pvalue = hl.float64(ht.Pvalue))
    ht = ht.drop('Tstat', 'N_case', 'N_ctrl')
    ht = ht.checkpoint(output_ht_path, overwrite=overwrite, _read_if_exists=not overwrite)
    ht.describe()
    ht.show()
    ht.filter(hl.is_nan(ht.Pvalue) | (ht.Pvalue < 1e-100)).show()
    print(ht.count())
    if variant_type == 'genome':
        ht = ht.filter(hl.case()
                       .when(ht.Pvalue > 0.1, hl.rand_bool(0.001))
                       .when(ht.Pvalue > 0.01, hl.rand_bool(0.01))
                       .when(ht.Pvalue > 0.001, hl.rand_bool(0.1))
                       .default(True))
    else:
        ht = ht.filter(hl.case()
                       .when(ht.Pvalue > 0.1, hl.rand_bool(0.01))
                       .when(ht.Pvalue > 0.01, hl.rand_bool(0.1))
                       .default(True))
    ht = ht.checkpoint(output_ht_path.replace('.ht', '_downsampled.ht'), overwrite=overwrite, _read_if_exists=not overwrite)
    ht.export(output_ht_path.replace('.ht', '.txt.bgz'))
    print(ht.count())


def run_saige(
    p: Batch,
    phenoname: str,
    pop:str,
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
        p.new_job(name=f"run_saige{'' if analysis_type=='variant' else '_gene'}_{phenoname}_{pop}", attributes={"analysis_type": analysis_type, "pop": pop})
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
    pop:str,
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

    pheno_col = "value"
    user_id_col = "person_id"
    in_bfile = p.read_input_group(
        **{ext: f"{plink_file_root}.{ext}" for ext in ("bed", "bim", "fam")}
    )
    fit_null_task = (
        p.new_job(
            name=f"fit_null_model_{phenoname}_{pop}",
            attributes={"trait_type": trait_type, "pop":pop},
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
    mt = get_filtered_mt(analysis_type=analysis_type, filter_variants=True, filter_samples=True, adj_filter= True, pop=pop, prune_samples=True)

    # Filter to interval
    mt = hl.filter_intervals(mt, [interval])

    mt = mt.select_entries("GT")
    mt = mt.filter_rows(
        hl.agg.count_where(mt.GT.is_non_ref()) > 0
    )  # Filter to non-reference sites
    mt = mt.annotate_rows(
        rsid=mt.locus.contig + ":" + hl.str(mt.locus.position) + "_" + mt.alleles[0] + "/" + mt.alleles[1]
    )  # Annotate rsid

    call_stats_ht = hl.read_table(get_call_stats_ht_path(pop=pop, pruned=True, analysis_type=analysis_type))
    call_stats_ht = call_stats_ht.filter(
        (call_stats_ht.call_stats.AC[1] > variant_ac_filter) &
        (call_stats_ht.call_stats.AN/(2*N_SAMPLES_PRUNED[pop]) > variant_callrate_filter)
    )
    mt = mt.filter_rows(hl.is_defined(call_stats_ht[mt.row_key]))

    mt = mt.annotate_entries(
        GT=hl.if_else(mt.GT.is_haploid(), hl.call(mt.GT[0], mt.GT[0]), mt.GT)
    )
    mt = gt_to_gp(mt)
    mt = impute_missing_gp(mt, mean_impute=mean_impute_missing)
    hl.export_bgen(mt, f"{output_dir}/{outname}", gp=mt.GP, varid=mt.rsid)
    print(hl.utils.range_table(10)._force_count())


def export_gene_group_file(interval, pop, output_dir):
    gene_ht = hl.read_table(get_aou_gene_map_ht_path(pop=pop, processed=True))
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

def index_bgen(b: hb.batch.Batch, pop:str, bgen: str, depend_job=None):
    file = b.read_input(bgen)
    name = bgen.split('/')[-1]
    j = b.new_job(name=f'index_{name}_{pop}')
    if depend_job is not None:
        j.depends_on(depend_job)
    j.image('befh/bgen:latest')
    j.command(f'bgenix -index -g {file}')
    j.command(f'ls {file}.bgi')
    j.command(f'mv {file}.bgi {j.temp}')
    b.write_output(j.temp, f'{bgen}.bgi')
    return j

def export_pheno(category: str, pheno_path:str, covariates: str, pop: str, phenoname: str, binary_traits: list, phenotype_categories_in_mt: list, proportion_single_sex:float):
    def irnt(he: hl.expr.Expression, output_loc: str = 'irnt'):
        ht = he._indices.source
        n_rows = ht.aggregate(hl.agg.count_where(hl.is_defined(he)))
        print(n_rows)
        ht = ht.order_by(he).add_index()
        ht = ht.annotate(**{output_loc: hl.qnorm((ht.idx + 0.5) / n_rows)})
        ht = ht.annotate(**{output_loc: hl.or_missing(~hl.is_nan(ht[output_loc]), ht[output_loc])})
        return ht

    def filter_to_single_sex_by_proportion(ht, proportion_single_sex):
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
        elif prop_female >= 1 - proportion_single_sex:
            print(
                f"Female case proportion {prop_female} greater than {1 - proportion_single_sex}. Filtering to females..."
            )
            ht = ht.filter(ht.sex == 0)
        return ht

    def load_pheno_ht(category:str, pop:str, phenotype_categories_in_mt:list):
        def get_raw_phenotype_path(
                name: str, parsed: bool = True, extension: str = "ht", tranche: str = TRANCHE
        ):
            return (
                f'{DATA_PATH}/phenotype/{"" if parsed else "raw/"}{name}_{tranche}.{extension}'
            )

        def get_lab_measurements_path(
                name: str, parsed: bool, extension: str = "ht", tranche: str = TRANCHE
        ):
            return (
                f'{DATA_PATH}/phenotype/lab_measurements/{"" if parsed else "raw/"}measurement_{name}_{tranche}.{extension}'
            )
        if category in phenotype_categories_in_mt:
            mt = hl.read_matrix_table(get_raw_phenotype_path(name=category, parsed=True, extension=f'hard_filtered.mt'))
            pop_mt = mt.filter_rows(mt.pop == pop)
            pheno_mt = pop_mt.filter_cols(pop_mt.phenocode == phenoname)
            ht = pheno_mt.entries()
        elif category == 'random_phenotypes':
            pheno_ht = hl.read_table(get_raw_phenotype_path(name=category, parsed=True, extension=f'{pop}.ht'))
            ht = pheno_ht.annotate(value=pheno_ht[phenoname])
        elif category == 'lab_measurements':
            ht = hl.read_table(get_lab_measurements_path(name=phenoname, parsed=True))
            ht = ht.filter(ht.pop == pop)
        overwrite = False
        if not hfs.exists(f'gs://aou_tmp/export_pheno/{pop}_{phenoname}.ht/_SUCCESS'):
            overwrite = True
        ht = ht.checkpoint(f'gs://aou_tmp/export_pheno/{pop}_{phenoname}.ht', _read_if_exists= not overwrite, overwrite = overwrite)
        # print(f"Exporting {phenoname} for {pop.upper()} (n_samples: {ht.count()})...")
        return ht

    ht = load_pheno_ht(category=category, pop=pop, phenotype_categories_in_mt=phenotype_categories_in_mt)
    ht.describe()

    if phenoname in binary_traits:
        ht = ht.annotate(value=hl.int(ht.value))
        ht = filter_to_single_sex_by_proportion(ht=ht, proportion_single_sex=proportion_single_sex)
    else:
        print('Computing irnt values...')
        ht = ht.annotate(value = hl.float64(ht['value']))
        ht = irnt(ht['value'], 'irnt_value')
        ht = ht.annotate(value=ht.irnt_value)
    ht = ht.key_by('person_id')
    fields = covariates.split(",") + ["value"]
    out_ht = ht.select(**{field: ht[field] for field in fields})
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
        app_name=f'aou_SAIGE_{args.pops.replace(",","_")}' if not args.export_phenos else "aou_export_phenos"
    )

    num_pcs = 20
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
    proportion_single_sex = 0.1
    pops = args.pops.split(",")


    if args.export_phenos or (not args.skip_saige) or (not args.skip_any_null_models) or (not args.skip_load_hail_results):
        if not hfs.exists(get_phenotype_info_path(version="pheno_dict", extension="dict")):
            all_phenos_by_group = {}
            all_phenos_by_group['lab_measurements'] = LAB_CODES
            if not args.skip_any_random_pheno:
                raw_random_pheno_ht = hl.import_table(
                    get_random_phenotype_path(pop=pops[0]),
                    key="userId",
                    impute=True,
                    types={"userId": hl.tstr},
                )
                random_phenos_to_run = list(raw_random_pheno_ht.row_value)
                all_phenos_by_group['random_phenotypes'] = random_phenos_to_run

            for category in TRAIT_TYPEs:
                print(f'Loading {category} information......')
                if category in ['lab_measurements', 'random_phenotypes']:
                    continue
                ht = hl.read_table(get_phenotype_info_path(version=category, extension=f'hard_filtered.ht'))
                all_phenos_by_group[category] = set(ht.phenocode.collect())
            write_pickle_dict(output=get_phenotype_info_path(version="pheno_dict", extension="dict"), dict =all_phenos_by_group)

        all_phenos_by_group = read_pickle_dict(get_phenotype_info_path(version="pheno_dict", extension="dict"))
        print(f'----------Number of phenotypes per category (RAW): --------------')
        print([(category, len(all_phenos_by_group[category])) for category in list(all_phenos_by_group.keys())])

        quantitative_traits = []
        binary_traits = []
        if 'random_phenotypes' in list(all_phenos_by_group.keys()):
            quantitative_traits = [pheno for pheno in all_phenos_by_group['random_phenotypes'] if 'continuous' in pheno]
            binary_traits = [pheno for pheno in all_phenos_by_group['random_phenotypes'] if 'continuous' not in pheno]
        for category in quantitative_categories:
            quantitative_traits = quantitative_traits + list(all_phenos_by_group[category])
        for category in binary_categories:
            binary_traits = binary_traits + list(all_phenos_by_group[category])
        print(f'Number of quantitative traits: {len(quantitative_traits)}\nNumber of binary traits: {len(binary_traits)}')

        phenos_to_run_by_pop = {}
        phenos_to_run_by_pop_by_group = {}
        if args.pilot or (args.phenos is not None):
            phenos_to_run = list(PILOT_PHENOTYPES) if args.pilot else args.phenos.split(",")
            category_lst = []
            for phenoname in phenos_to_run:
                for category in TRAIT_TYPEs:
                    if phenoname in all_phenos_by_group[category]:
                        category_lst.append(category)
                        break
            category_dict = dict(zip(phenos_to_run, category_lst))

            for pop in pops:
                phenos_to_run_by_pop_by_group[pop] = {}
                for category in set(category_lst):
                    phenos_to_run_by_pop_by_group[pop][category] = [k for k, v in category_dict.items() if v == category]
                phenos_to_run_by_pop[pop] = phenos_to_run
            if len(phenos_to_run) < 20:
                print(phenos_to_run_by_pop)
        else:
            if not hfs.exists(get_phenotype_info_path(version="pheno_by_pop_dict", extension="dict")) or args.overwrite_dict:
                lab_phenos_by_pop = read_pickle_dict(get_phenotype_info_path(version="lab_by_pop_dict", extension="dict"))
                for pop in pops:
                    phenos_to_run_by_pop[pop] = lab_phenos_by_pop[pop]
                    phenos_to_run_by_pop_by_group[pop] = {}
                    phenos_to_run_by_pop_by_group[pop]['lab_measurements'] = lab_phenos_by_pop[pop]
                    if 'random_phenotypes' in list(all_phenos_by_group.keys()):
                        random_pheno_lst = [x for x in all_phenos_by_group['random_phenotypes'] if x.endswith(('_1', '_2', '_3', '_4', '_5'))]
                        print(random_pheno_lst)
                        phenos_to_run_by_pop_by_group[pop]['random_phenotypes'] = random_pheno_lst
                        phenos_to_run_by_pop[pop] = phenos_to_run_by_pop[pop] + random_pheno_lst

                for category in TRAIT_TYPEs:
                    if category in ['lab_measurements', 'random_phenotypes']: continue
                    ht = hl.read_table(get_phenotype_info_path(version=category, extension=f'hard_filtered.ht'))
                    for pop in pops:
                        print(f'--------Loading {category.replace("_table", "")} phenotypes for {pop.upper()}--------')
                        pop_ht = ht.filter(ht[f'n_cases_{pop}'] >= 200)
                        pop_ht = pop_ht.checkpoint(f'gs://aou_tmp/pheno_info/{pop}_{category}_n_cases_over_200.ht', _read_if_exists=True)
                        pheno_lst = pop_ht.phenocode.collect()
                        phenos_to_run_by_pop_by_group[pop][category] = pheno_lst
                        phenos_to_run_by_pop[pop] = phenos_to_run_by_pop[pop] + pheno_lst
                meta_pheno_lst = []
                for pop in pops:
                    meta_pheno_lst = meta_pheno_lst + phenos_to_run_by_pop[pop]
                phenos_to_run_by_pop['meta'] = list(set(meta_pheno_lst))
                print(f"META pheno length: {len(phenos_to_run_by_pop['meta'])}")

                write_pickle_dict(output=get_phenotype_info_path(version="pheno_by_pop_by_group_dict", extension="dict"),
                                  dict=phenos_to_run_by_pop_by_group)
                write_pickle_dict(output=get_phenotype_info_path(version="pheno_by_pop_dict", extension="dict"),
                                  dict=phenos_to_run_by_pop)

            phenos_to_run_by_pop = read_pickle_dict(get_phenotype_info_path(version="pheno_by_pop_dict", extension="dict"))
            phenos_to_run_by_pop_by_group = read_pickle_dict(get_phenotype_info_path(version="pheno_by_pop_by_group_dict", extension="dict"))
            print(f'----------Number of phenotypes per category per pop (filtered to n_cases >= 200): --------------')
            print([f"{pop.upper()}-{category}: {len(phenos_to_run_by_pop_by_group[pop][category])}" for pop in pops for category in list(all_phenos_by_group.keys())])
            print([f"{pop.upper()}: {len(phenos_to_run_by_pop[pop])}" for pop in pops])
    else:
        phenos_to_run_by_pop = {}
        for pop in POPS:
            phenos_to_run_by_pop[pop] = None
        phenos_to_run_by_pop_by_group = None

    backend = hb.ServiceBackend(
        billing_project="all-by-aou", remote_tmpdir=TMP_BUCKET
    )

    if args.export_phenos:
        b = hb.Batch(
            name=f"aou_export_pheno_{args.pops}_{'test' if args.test or args.phenos is not None else 'all'}",
            backend=backend,
            default_storage="500Mi",
            default_cpu=n_threads,
        )
        phenotype_categories = ['lab_measurements'] # TODO: remove later after this
        for category in phenotype_categories:
            for pop in pops:
                if category not in list(phenos_to_run_by_pop_by_group[pop].keys()): continue
                phenos_to_export = phenos_to_run_by_pop_by_group[pop][category]
                print(f"------------{pop.upper()} phenotype info: {len(phenos_to_export)} {category.replace('_table', '')} phenotypes------------")
                if category == 'random_phenotypes' and not hfs.exists(get_raw_phenotype_path(name = category, parsed=True, extension=f'{pop}.ht')):
                    pheno_ht = hl.import_table(
                        get_random_phenotype_path(pop=pop),
                        impute=True,
                        types={"userId": hl.tstr},
                    )
                    pheno_ht = pheno_ht.annotate(person_id = pheno_ht["userId"])
                    pheno_ht = pheno_ht.key_by('person_id')
                    pheno_ht = pheno_ht.drop("userId")
                    meta_ht = hl.read_table(get_sample_meta_path(annotation=True, extension='ht'))
                    pca_ht = hl.read_table(get_pca_ht_path(pop='all', name=f'pruned_full_scores'))
                    pheno_ht = pheno_ht.annotate(**meta_ht[pheno_ht.key])
                    pheno_ht = pheno_ht.annotate(**pca_ht[pheno_ht.key])
                    pheno_ht.checkpoint(get_raw_phenotype_path(name = category, parsed=True, extension=f'{pop}.ht'), _read_if_exists=True)
                for phenoname in phenos_to_export:
                    pheno_path = f'{get_aou_saige_utils_root(pop=pop, name="pheno_file")}/phenotype_{phenoname}.tsv'
                    if (not hfs.exists(pheno_path)
                            or args.overwrite_pheno_data
                    ):
                        j = b.new_python_job(f'export_phenotype_{phenoname}_for_{pop}',
                                             attributes={"pop": pop, "phenotype": copy.deepcopy(phenoname), "category": category})
                        j.image("hailgenetics/hail:0.2.127-py3.9")
                        j.call(export_pheno,
                               category=category,
                               pheno_path=pheno_path,
                               covariates=covariates,
                               pop=pop,
                               phenoname=phenoname,
                               binary_traits=binary_traits,
                               phenotype_categories_in_mt=phenotype_categories_in_mt,
                               proportion_single_sex=proportion_single_sex
                               )
        b.run()

    if args.run_pipeline:
        analysis_type = "variant" if args.single_variant_only else "gene"
        variant_type = 'genome' if analysis_type == 'variant' else 'exome'
        print(f'Analysis type: {analysis_type}')
        print(f'Variant_type: {variant_type}')
        print(f"Docker image: {SAIGE_DOCKER_IMAGE}...")

        for pop in pops:
            if not args.skip_bgen or not args.skip_saige:
                N_GENE_PER_GROUP = 40 if pop=='all' else N_GENE_PER_GROUP
                size = CHUNK_SIZE[pop] if analysis_type == 'variant' else N_GENE_PER_GROUP
                if not hfs.exists(
                        get_saige_interval_path(analysis_type=analysis_type, chunk_size=size, extension='pickle')):
                    interval_ht = hl.read_table(get_saige_interval_path(analysis_type=analysis_type, chunk_size=size))
                    interval_ht = interval_ht.filter(interval_ht.interval.start.contig != "chrM")
                    intervals = interval_ht.aggregate(hl.agg.collect(interval_ht.interval))
                    write_pickle_dict(
                        get_saige_interval_path(analysis_type=analysis_type, chunk_size=size, extension='pickle'),
                        intervals)
                intervals = read_pickle_dict(
                    get_saige_interval_path(analysis_type=analysis_type, chunk_size=size, extension='pickle'))
                print(intervals[0:3])
                print(f'---------Number of intervals [{pop.upper()}]: {len(intervals)}---------')
            if args.category is not None:
                category = args.category
                if args.category == 'small':
                    phenos_to_run = set(phenos_to_run_by_pop_by_group[pop]["lab_measurements"] + phenos_to_run_by_pop_by_group[pop]['random_phenotypes'] + \
                                    phenos_to_run_by_pop_by_group[pop]['pfhh_survey_table']+ phenos_to_run_by_pop_by_group[pop]["processed_physical_measurement_table"])
                elif args.category == 'quantitative':
                    phenos_to_run = set(
                        phenos_to_run_by_pop_by_group[pop]["lab_measurements"] + phenos_to_run_by_pop_by_group[pop][
                            "processed_physical_measurement_table"])
                else:
                    category = args.category
                    phenos_to_run = phenos_to_run_by_pop_by_group[pop][category]
                print(list(phenos_to_run)[0:5])
            else:
                category = 'all'
                phenos_to_run_by_pop['all'] = None
                phenos_to_run = phenos_to_run_by_pop[pop]

            b = hb.Batch(
                name=f"saige_{analysis_type}_aou_{pop}_{category.replace('_table', '')}",
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

            overwrite_null_models = args.overwrite_null_models
            null_model_dir = get_aou_saige_results_root(analysis_type=analysis_type, pop=pop, name ='null_glmm')
            null_models_already_created = {}
            null_models = {}
            pheno_exports = {}
            if phenos_to_run is not None:
                print(f"------------{pop.upper()} {analysis_type} analysis null models: {len(phenos_to_run)} phenotypes------------")
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
                    current_trait_type = 'binary' if phenoname in binary_traits else 'quantitative'

                    pheno_file = b.read_input(f'{get_aou_saige_utils_root(pop=pop, name="pheno_file")}/phenotype_{phenoname}.tsv')
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
                        print(f'Running null model for {pop.upper()} {phenoname} : {current_trait_type}')
                        fit_null_task = fit_null_glmm(
                            b,
                            pop=pop,
                            output_root=null_glmm_root,
                            phenoname=phenoname,
                            pheno_file=pheno_exports[phenoname],
                            trait_type=current_trait_type,
                            covariates=covariates,
                            plink_file_root=get_aou_sites_for_grm_path(pop=pop, extension="plink", pruned=True),
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
                            {"phenotype": copy.deepcopy(phenoname)}
                        )
                        model_file = fit_null_task.null_glmm.rda
                        variance_ratio_file = fit_null_task.null_glmm[f"varianceRatio.txt"]

                    null_models[phenoname] = (model_file, variance_ratio_file)
            else:
                print('No phenotype loaded...')

            if not args.skip_bgen or not args.skip_saige:
                print(f"------------{pop.upper()} {analysis_type} analysis bgen files------------")
                bgen_dir = get_aou_saige_results_root(analysis_type=analysis_type, pop=pop, name="bgen")
                print(f'bgen directory: {bgen_dir}')
                overwrite_bgens = args.overwrite_bgens
                bgens_already_created = {}
                if not overwrite_bgens and hfs.exists(bgen_dir):
                    bgens_already_created = {x["path"] for x in hl.hadoop_ls(bgen_dir)}
                factor = 3 if analysis_type=='variant' else 4
                print(f'Found {int(len(bgens_already_created)/factor)} Bgens in directory...')

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
                            name=f"{analysis_type}_analysis_export_{str(interval)}_bgen_{pop}"
                        )

                        bgen_task.image(image)
                        bgen_task.memory('highmem')
                        bgen_task.always_copy_output()
                        bgen_task.env('PYSPARK_SUBMIT_ARGS', '--driver-memory 24g --executor-memory 24g pyspark-shell')
                        bgen_task.call(
                            export_bgen_from_mt,
                            pop=pop,
                            analysis_type=analysis_type,
                            interval=interval,
                            output_dir=bgen_dir,
                            log_file=bgen_task.log_file,
                            mean_impute_missing=True,
                            variant_ac_filter= args.variant_ac_filter,
                            variant_callrate_filter=args.callrate_filter,
                        )
                        b.write_output(bgen_task.log_file, f'gs://aou_tmp/000_bgen_logs/{analysis_type}_analysis_export_{str(interval)}_bgen_{pop}.log')
                        bgen_task.attributes["pop"] = pop
                        bgen_task.attributes["analysis_type"] = analysis_type

                        bgen_index = index_bgen(b=b,
                                                pop=pop,
                                                bgen=f"{bgen_root}.bgen",
                                                depend_job=bgen_task)
                        bgen_index.attributes["pop"] = pop
                        bgen_index.attributes["analysis_type"] = analysis_type
                    if ((f"{bgen_root}.gene.txt" not in bgens_already_created)  or args.overwrite_gene_txt)and analysis_type=='gene':
                        gene_txt_task = b.new_python_job(
                            name=f"{analysis_type}_analysis_export_{str(interval)}_gene_txt_{pop}"
                        )
                        gene_txt_task.image(HAIL_DOCKER_IMAGE)
                        gene_txt_task.call(
                            export_gene_group_file,
                            interval=interval,
                            pop=pop,
                            output_dir=f"{bgen_root}.gene.txt",
                        )
                        gene_txt_task.attributes["pop"] = pop
                    if (not hfs.exists(f"{bgen_root}.bgen.bgi")) and  (hfs.exists(f"{bgen_root}.bgen")):
                        bgen_index=index_bgen(b=b, pop=pop, bgen=f"{bgen_root}.bgen")
                        bgen_index.attributes["pop"] = pop
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


            result_dir = get_aou_saige_results_root(
                analysis_type=analysis_type, pop=pop, name="result"
            )
            print(f'result directory: {result_dir}')
            overwrite_results = args.overwrite_results
            saige_tasks = {}
            saige_ur_tasks = {}
            if not args.skip_saige:
                memory = 'highmem' if pop in ['eur', 'all'] and analysis_type == 'variant' else 'standard'
                groups=None
                if analysis_type == "gene":
                    if args.groups is None:
                        groups = ','.join(["pLoF", "missenseLC", "synonymous"])
                    else:
                        groups = args.groups
                    print(f'Groups to run: {groups}')
                print(f"------------{pop.upper()} {analysis_type} analysis step 2: {len(phenos_to_run)} phenotypes------------")
                for i in tqdm(range(len(phenos_to_run))):
                    phenoname = list(phenos_to_run)[i]
                    current_trait_type = 'binary' if phenoname in binary_traits else 'quantitative'
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
                            f'-------------- Found {int(len(exome_results_already_created))} exome results in [{pop}: {phenoname}] directory...')
                    factor = 2 if analysis_type == 'variant' else 5
                    print(f'-------------- Found {int(len(results_already_created) / factor)} {analysis_type} results in [{pop}: {phenoname}] directory...')

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
                        samples_file = b.read_input(get_aou_sample_file_path(pop=pop))

                        if (
                            overwrite_results
                            or (f"{results_path}.{'gene' if analysis_type == 'gene' else 'single_variant'}.txt"
                            not in results_already_created)
                        ):
                            saige_task = run_saige(
                                p=b,
                                phenoname=phenoname,
                                pop=pop,
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
                                {"interval": str(interval), "pop": pop, 'name': f'saige{analysis_type == "gene"}'}
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
                                memory = 'highmem' if pop in ['eur', 'all']  else 'standard'
                                saige_ur_task = run_saige(
                                    p=b,
                                    phenoname=phenoname,
                                    pop=pop,
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
                                    {"interval": str(interval), "pop": pop}
                                )
                                saige_ur_task.attributes.update(
                                    {"phenotype": copy.deepcopy(phenoname)}
                                )
                                saige_ur_tasks[phenoname].append(saige_ur_task)
                    if args.test:
                        break

            if not args.skip_load_hail_results:
                test_type = 'saige' if analysis_type == 'variant' else 'saige_gene'
                root = f'{RESULTS_PATH}/{test_type}_results/{pop.upper()}'
                print(f'--------------Loading results from {analysis_type} analysis [{pop.upper()}]: {root}------------')
                for i in tqdm(range(len(phenos_to_run))):
                    phenoname = list(phenos_to_run)[i]
                    directory = f'{get_aou_saige_results_root(analysis_type=analysis_type, pop=pop, name="result")}/phenotype_{phenoname}'
                    output_ht_directory = f'{root}/phenotype_{phenoname}'
                    if (hfs.exists(f'{output_ht_directory}/variant_results.ht/_SUCCESS') and not args.overwrite_hail_results):
                        continue
                    null_glmm_log = f'{get_aou_saige_results_root(analysis_type=analysis_type, pop=pop, name ="null_glmm")}/phenotype_{phenoname}.log'

                    saige_log = f'{directory}/result_{phenoname}_chr1_{"000065419" if analysis_type == "gene" else "000000001"}.{analysis_type}.log'

                    quantitative_trait = phenoname in quantitative_traits

                    if not args.skip_load_gene_results and analysis_type == 'gene':
                        if (not hfs.exists(f'{output_ht_directory}/gene_results.ht/_SUCCESS')) or args.overwrite_hail_results:
                            print(f'[{pop.upper()}: {phenoname}] Gene table')
                            j = b.new_python_job(
                                name=f'sync_saige_{analysis_type}_gene_HT_{phenoname}_{pop}',
                                attributes={"analysis_type": analysis_type, "pop": pop,
                                            "phenotype": copy.deepcopy(phenoname)})
                            j.image("hailgenetics/hail:0.2.127-py3.9")
                            j.memory('standard')
                            j.cpu(16)
                            j.env('PYSPARK_SUBMIT_ARGS', '--driver-memory 8g --executor-memory 8g pyspark-shell')
                            gene_SUCCESS = hfs.exists(f'{output_ht_directory}/gene_results.ht/_SUCCESS')
                            j.call(load_gene_data,
                                   directory=directory,
                                   output_ht_directory=output_ht_directory,
                                   phenoname=phenoname,
                                   gene_ht_map_path=get_aou_gene_map_ht_path(pop=pop, processed=False),
                                   quantitative_trait=quantitative_trait,
                                   null_glmm_log=null_glmm_log,
                                   saige_log=saige_log,
                                   overwrite=args.overwrite_hail_results or not gene_SUCCESS
                                   )
                            if phenoname in list(saige_tasks.keys()) and not args.skip_saige:
                                j.depends_on(*saige_tasks[phenoname])
                            j.always_run()

                    if not args.skip_load_variant_results:
                        extension = 'single_variant.txt'
                        variant_type = 'genome' if analysis_type == 'variant' else 'exome'
                        j = b.new_python_job(
                            name=f'sync_saige_{variant_type}_variant_HT_{phenoname}_{pop}',
                            attributes={"analysis_type": analysis_type, "pop": pop, "phenotype": copy.deepcopy(phenoname)})
                        j.image("hailgenetics/hail:0.2.127-py3.9")
                        j.memory('standard')
                        j.cpu(16)
                        j.env('PYSPARK_SUBMIT_ARGS', '--driver-memory 8g --executor-memory 8g pyspark-shell')
                        SUCCESS = hfs.exists(f'{output_ht_directory}/variant_results.ht/_SUCCESS')

                        print(f'[{pop.upper()}: {phenoname}] {variant_type} variant table')
                        j.call(load_variant_data,
                                   directory=directory,
                                   output_ht_directory=output_ht_directory,
                                   phenoname=phenoname,
                                   quantitative_trait=quantitative_trait,
                                   null_glmm_log=null_glmm_log,
                                   saige_log=saige_log,
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
        "--pops", help="comma-separated list", default="afr,amr,eas,eur,mid,sas"
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
        default ='hailgenetics/hail:0.2.130-py3.9'
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
