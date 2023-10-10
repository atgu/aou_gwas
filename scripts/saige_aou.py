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

# from aou_gwas import *  # for dataproc
# from utils.utils import *  # for QoB
# from utils.resources import *  # for QoB
import time


logging.basicConfig(
    format="%(levelname)s (%(name)s %(lineno)s): %(message)s",
    level="INFO",
    filename="saige_pipeline.log",
)
logger = logging.getLogger("ALL_x_AoU_saige")
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
    vcf_file: ResourceGroup,
    samples_file: ResourceGroup,
    docker_image: str,
    group_file: str = None,
    sparse_sigma_file: str = None,
    use_bgen: bool = True,
    trait_type: str = "continuous",
    chrom: str = "chr1",
    min_mac: int = 1,
    min_maf: float = 0,
    max_maf: float = 0.5,
    memory: str = "",
    storage: str = "10Gi",
    add_suffix: str = "",
    log_pvalue: bool = False,
):
    analysis_type = "gene" if sparse_sigma_file is not None else "variant"
    run_saige_task: Job = (
        p.new_job(name=f"run_saige", attributes={"analysis_type": analysis_type})
        .cpu(1)
        .storage(storage)
        .image(docker_image)
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
        f"--minMAF={min_maf} "
        f"--minMAC={min_mac} "
        f"--maxMAFforGroupTest={max_maf} "
        f"--sampleFile={samples_file} "
        f"--GMMATmodelFile={model_file} "
        f'{"--IsOutputlogPforSingle=TRUE " if log_pvalue else ""}'
        f"--varianceRatioFile={variance_ratio_file} "
        f"--LOCO=FALSE "
        f"--SAIGEOutputFile={run_saige_task.result} "
    )

    if use_bgen:
        command += (
            f"--bgenFile={vcf_file.bgen} " f'--bgenFileIndex={vcf_file["bgen.bgi"]} '
        )
    else:
        command += (
            f'--vcfFile={vcf_file["vcf.gz"]} '
            f'--vcfFileIndex={vcf_file["vcf.gz.tbi"]} '
            f"--chrom={chrom} "
            f"--vcfField=GT "
        )
    if analysis_type == "gene":
        if trait_type == "binary":
            command += f"--IsOutputPvalueNAinGroupTestforBinary=TRUE "
        command += (
            f"--groupFile={group_file} "
            f"--sparseSigmaFile={sparse_sigma_file} "
            f"--IsSingleVarinGroupTest=TRUE "
            f"--IsOutputBETASEinBurdenTest=TRUE "
        )
    command += f"--IsOutputAFinCaseCtrl=TRUE 2>&1 | tee {run_saige_task.stdout}; "
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
    min_covariate_count: int = 10,
    n_threads: int = 16,
    storage: str = "10Gi",
    memory: str = "60G",
    non_pre_emptible: bool = False,
):
    analysis_type = "variant" if sparse_grm is None else "gene"
    pheno_col = "value"
    user_id_col = "person_id"
    in_bfile = p.read_input_group(
        **{ext: f"{plink_file_root}.{ext}" for ext in ("bed", "bim", "fam")}
    )
    fit_null_task = (
        p.new_job(
            name=f"fit_null_model_{phenoname}",
            attributes={"analysis_type": analysis_type, "trait_type": trait_type},
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
            "_30markers.SAIGE.results.txt",
            f"{analysis_type}.varianceRatio.txt",
        )
    }
    if analysis_type == "gene":
        sparse_sigma_extension = sparse_grm_extension.replace("GRM", "Sigma")
        output_files[
            f"{analysis_type}.varianceRatio.txt{sparse_sigma_extension}"
        ] = f"{{root}}.{analysis_type}.varianceRatio.txt{sparse_sigma_extension}"
    fit_null_task.declare_resource_group(null_glmm=output_files)
    bim_fix_command = f"perl -pi -e s/^chr// {in_bfile.bim}"
    # if trait_type == 'icd':
    #     bim_fix_command += (f"; zcat {pheno_file.gz} | perl -p -e 's/true/1/g' | perl -p -e 's/false/0/g' "
    #                         f"| gzip -c > {pheno_file.gz}.temp.gz; mv {pheno_file.gz}.temp.gz {pheno_file.gz}")

    command = (
        f"set -o pipefail; Rscript /usr/local/bin/step1_fitNULLGLMM.R "
        f"--plinkFile={in_bfile} "
        f"--phenoFile={pheno_file} "
        f"--covarColList={covariates} "
        f"--minCovariateCount={min_covariate_count} "
        f"--phenoCol={pheno_col} "
        f"--sampleIDColinphenoFile={user_id_col} "
        f"--traitType={trait_type} "
        f"--outputPrefix={fit_null_task.null_glmm} "
        f"--outputPrefix_varRatio={fit_null_task.null_glmm}.{analysis_type} "
        f"--skipModelFitting={str(skip_model_fitting).upper()} "
    )
    if inv_normalize:
        command += "--invNormalize=TRUE "
    if analysis_type == "gene":
        fit_null_task.declare_resource_group(
            sparse_sigma={sparse_sigma_extension: f"{{root}}.{sparse_sigma_extension}"}
        )
        command += (
            f"--IsSparseKin=TRUE "
            f"--sparseGRMFile={sparse_grm[sparse_grm_extension]} "
            f'--sparseGRMSampleIDFile={sparse_grm[f"{sparse_grm_extension}.sampleIDs.txt"]} '
            f"--isCateVarianceRatio=TRUE "
        )
    command += f"--nThreads={n_threads} --LOCO=FALSE 2>&1 | tee {fit_null_task.stdout}"
    command = "; ".join([bim_fix_command, command])
    fit_null_task.command(command)
    p.write_output(fit_null_task.null_glmm, output_root)
    p.write_output(fit_null_task.stdout, f"{output_root}.{analysis_type}.log")
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
    mt_path,
    pop,
    analysis_type,
    interval,
    output_dir,
    mean_impute_missing,
    no_adj,
    callrate_filter,
):
    mt = hl.read_matrix_table(mt_path)
    meta_ht = hl.read_table(get_sample_meta_path(annotation=True, extension="ht"))
    meta_ht = meta_ht.filter((meta_ht.pop == pop) & (meta_ht.samples_to_keep))
    mt = mt.filter_cols(hl.is_defined(meta_ht[mt.col_key]))

    outname = f"{analysis_type}_{interval.start.contig}_{interval.start.position}_{interval.end.position}"
    print(outname)

    # Filter to interval
    mt = hl.filter_intervals(mt, [hl.parse_locus_interval(interval)])

    # if not no_adj:
    #     mt = mt.filter_entries(mt.adj)

    mt = mt.select_entries("GT")
    mt = mt.filter_rows(
        hl.agg.count_where(mt.GT.is_non_ref()) > 0
    )  # Filter to non-reference sites
    mt = mt.annotate_rows(
        rsid=mt.locus.contig + ":" + hl.str(mt.locus.position) + "_" + mt.alleles[0] + "/" + mt.alleles[1]
    )  # Annotate rsid

    if callrate_filter:
        mt = mt.filter_rows(
            hl.agg.fraction(hl.is_defined(mt.GT)) >= args.callrate_filter
        )

    mt = mt.annotate_entries(
        GT=hl.if_else(mt.GT.is_haploid(), hl.call(mt.GT[0], mt.GT[0]), mt.GT)
    )
    mt = gt_to_gp(mt)
    mt = impute_missing_gp(mt, mean_impute=mean_impute_missing)
    hl.export_bgen(mt, f"{output_dir}/{outname}", gp=mt.GP, varid=mt.rsid)


def export_gene_group_file(interval, groups, output_dir):
    gene_ht = hl.read_table(get_aou_util_path(name="gene_map", parsed=True))
    gene_ht = hl.filter_intervals(gene_ht, [hl.parse_locus_interval(interval)])

    gene_ht = gene_ht.filter(hl.set(groups).contains(gene_ht.annotation))
    if args.common_variants_only:
        gene_ht = gene_ht.filter(gene_ht.common_variant)
    gene_ht.select(
        group=gene_ht.gene_id + "_" + gene_ht.gene_symbol + "_" + gene_ht.annotation
        + hl.if_else(gene_ht.common_variant, "_" + gene_ht.variants[0], ""),
        variant=hl.delimit(gene_ht.variants, "\t"),
    ).key_by().drop("start").export(output_dir, header=False)


SCRIPT_DIR = "/ukb_common/saige"


def activate_service_account(task):
    task.env("GOOGLE_APPLICATION_CREDENTIALS", "/gsa-key/key.json")


def load_results_into_hail(
    p: Batch,
    output_root: str,
    pheno_keys,
    tasks_to_hold,
    vep_path: str,
    docker_image: str,
    gene_map_path: str = None,
    null_glmm_log: str = "",
    reference: str = "GRCh38",
    saige_log: str = "",
    analysis_type: str = "gene",
    n_threads: int = 8,
    storage: str = "10Gi",
    legacy_annotations: bool = False,
    log_pvalue: bool = False,
    overwrite: bool = True,
):
    load_data_task: Job = (
        p.new_job(name=f"load_data", attributes=copy.deepcopy(pheno_keys))
        .image(docker_image)
        .cpu(n_threads)
        .storage(storage)
        .memory("standard")
    )
    load_data_task.always_run().depends_on(*tasks_to_hold)
    pheno_dict_opts = " ".join([f"--{k} {shq(v)}" for k, v in pheno_keys.items()])
    python_command = f"""python3 {SCRIPT_DIR}/load_results.py
    --input_dir {shq(output_root)}
    {"--null_glmm_log " + shq(null_glmm_log) if null_glmm_log else ''}
    --saige_run_log_format {saige_log}
    {pheno_dict_opts}
    {"--gene_map_ht_raw_path " + gene_map_path if gene_map_path else ''}
    {"--legacy_annotations" if legacy_annotations else ""}
    --ukb_vep_ht_path {vep_path}
    {"--overwrite" if overwrite else ""} --reference {reference}
    --analysis_type {analysis_type} {"--log_pvalue" if log_pvalue else ""}
    --n_threads {n_threads} | tee {load_data_task.stdout}
    ;""".replace(
        "\n", " "
    )

    python_command = python_command.replace("\n", "; ").strip()
    command = (
        f'set -o pipefail; PYTHONPATH=$PYTHONPATH:/ PYSPARK_SUBMIT_ARGS="--conf spark.driver.memory={int(3 * n_threads)}g pyspark-shell" '
        + python_command
    )
    load_data_task.command(command)
    activate_service_account(load_data_task)
    p.write_output(
        load_data_task.stdout, f'{output_root}/{pheno_keys["phenocode"]}_loading.log'
    )
    return load_data_task


def qq_plot_results(
    p: Batch,
    output_root: str,
    tasks_to_hold,
    export_docker_image: str,
    R_docker_image: str,
    n_threads: int = 8,
    storage: str = "10Gi",
):
    qq_export_task: Job = (
        p.new_job(name="qq_export")
        .image(export_docker_image)
        .cpu(n_threads)
        .storage(storage)
    )
    qq_export_task.always_run().depends_on(*tasks_to_hold)

    python_command = f"""python3 {SCRIPT_DIR}/export_results_for_qq.py
    --input_dir {shq(output_root)}
    --output_file {qq_export_task.out}
    --n_threads {n_threads}
    ; """.replace(
        "\n", " "
    )

    command = (
        f'set -o pipefail; PYTHONPATH=$PYTHONPATH:/ PYSPARK_SUBMIT_ARGS="--conf spark.driver.memory={int(3 * n_threads)}g pyspark-shell" '
        + python_command
    )
    qq_export_task.command(command)
    activate_service_account(qq_export_task)

    qq_task: Job = (
        p.new_job(name="qq_plot")
        .image(R_docker_image)
        .cpu(n_threads)
        .storage(storage)
        .always_run()
    )
    qq_task.declare_resource_group(
        result={
            ext: f"{{root}}_Pvalue_{ext}"
            for ext in (
                "qqplot.png",
                "manhattan.png",
                "manhattan_loglog.png",
                "qquantiles.txt",
            )
        }
    )
    R_command = f"/saige-pipelines/scripts/qqplot.R -f {qq_export_task.out} -o {qq_task.result} -p Pvalue; "
    qq_task.command(R_command)

    p.write_output(qq_task.result, output_root)
    return qq_export_task, qq_task


def main(args):
    hl.init(
        # driver_memory="highmem",
        # driver_cores=8,
        # worker_memory="highmem",
        # worker_cores=1,
        default_reference="GRCh38",
        log="/saige_aou.log",
    )

    analysis_type = "variant" if args.single_variant_only else "gene"

    num_pcs = 10  # TODO: check number
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
    iteration = 1
    pops = args.pops.split(",") if args.pops else POPS

    pheno_ht = hl.read_table(
        get_full_phenotype_path(annotation=True, random_pheno=args.random_pheno)
    )
    print(f"Additional sample QC...")
    pheno_ht = pheno_ht.filter(pheno_ht.samples_to_keep)
    quantitative_phenos = hl.eval(pheno_ht.trait_type)["continuous"]

    if args.pilot or (args.phenos is not None):
        phenos_to_run = PILOT_PHENOTYPES if args.pilot else args.phenos.split(",")
    elif args.random_pheno:
        raw_random_pheno_ht = hl.import_table(
            get_random_phenotype_path(pop=pops[0], test=args.test),
            key="userId",
            impute=True,
            types={"userId": hl.tstr},
        )
        phenos_to_run = list(raw_random_pheno_ht.row_value)
    else:
        raw_pheno_ht = hl.read_table(
            get_full_phenotype_path(annotation=False, random_pheno=args.random_pheno)
        )
        phenos_to_run = list(raw_pheno_ht.row_value)
        # TODO: add phenotype filters

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
                if (not hl.hadoop_exists(get_aou_pheno_file_path(pop=pop, phenocode=phenoname))
                        or args.overwrite_pheno_data
                ):
                    print(f"Exporting {phenoname} for {pop.upper()}...")
                    if (phenoname in quantitative_phenos) or ("coutinuous" in phenoname):
                        ht = ht.annotate(value=ht[phenoname])
                    else:
                        ht = ht.annotate(value=hl.int(ht[phenoname]))
                        print(
                            f"N {phenoname} cases for {pop.upper()} samples: {ht.aggregate(hl.agg.sum(ht.value))}..."
                        )
                    fields = covariates.split(",") + ["related", "pop", "value"]
                    out_ht = ht.select(**{field: ht[field] for field in fields})
                    print(
                        f"N samples defined for {phenoname}: {out_ht.aggregate(hl.agg.count_where(hl.is_defined(out_ht.value)))}..."
                    )

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
                    out_ht.export(get_aou_pheno_file_path(pop=pop, phenocode=phenoname))

    if args.run_saige:
        backend = hb.ServiceBackend(
            billing_project="all-by-aou", remote_tmpdir=TMP_BUCKET
        )
        for pop in pops:
            logger.info(f"Docker image: {SAIGE_DOCKER_IMAGE}...")
            b = hb.Batch(
                name=f"saige_aou_{pop}",
                backend=backend,
                default_image=SAIGE_DOCKER_IMAGE,
                default_storage="500Mi",
                default_cpu=n_threads,
            )
            window = "1e7" if pop == "eur" and analysis_type == "variant" else "1e6"

            sparse_grm = None
            sparse_grm_extension = ""
            if analysis_type == "gene":
                print(
                    f"--------------Loading GRM (analysis type: {analysis_type}) -------------------"
                )
                relatedness_cutoff = "0.125"
                num_markers = 2000
                n_threads = 8
                sparse_grm_root = f"{DATA_PATH}/utils/grm/aou_{pop}"
                sparse_grm_extension = f"_relatednessCutoff_{relatedness_cutoff}_{num_markers}_randomMarkersUsed.sparseGRM.mtx"
                sparse_grm = b.read_input_group(
                    **{
                        ext: f"{sparse_grm_root}.{ext}"
                        for ext in (
                            sparse_grm_extension,
                            f"{sparse_grm_extension}.sampleIDs.txt",
                        )
                    }
                )

            overwrite_null_models = args.create_null_models
            null_model_dir = get_aou_saige_root(
                analysis_type=analysis_type, pop=pop, name="null_glmm", test=args.test
            )
            null_models_already_created = {}
            if not overwrite_null_models and hl.hadoop_exists(null_model_dir):
                null_models_already_created = {
                    x["path"] for x in hl.hadoop_ls(null_model_dir)
                }
            null_models = {}
            pheno_exports = {}
            for phenoname in phenos_to_run:
                trait_type = (
                    "quantitative"
                    if (phenoname in quantitative_phenos) or ("coutinuous" in phenoname)
                    else "binary"
                )
                null_glmm_root = f"{null_model_dir}/phenotype_{phenoname}"
                model_file_path = f"{null_glmm_root}.rda"
                variance_ratio_file_path = (
                    f"{null_glmm_root}.{analysis_type}.varianceRatio.txt"
                )
                pheno_file = b.read_input(
                    get_aou_pheno_file_path(pop=pop, phenocode=phenoname)
                )
                pheno_exports[phenoname] = pheno_file
                sparse_sigma_file = None
                if analysis_type == "gene":
                    sparse_sigma_file_path = (
                        variance_ratio_file_path
                        + sparse_grm_extension.replace("GRM", "Sigma")
                    )

                if (
                    not overwrite_null_models
                    and model_file_path in null_models_already_created
                    and variance_ratio_file_path in null_models_already_created
                ):
                    model_file = b.read_input(model_file_path)
                    variance_ratio_file = b.read_input(variance_ratio_file_path)
                    if analysis_type == "gene":
                        sparse_sigma_file = b.read_input(sparse_sigma_file_path)
                else:
                    if args.skip_any_null_models:
                        continue
                    fit_null_task = fit_null_glmm(
                        b,
                        null_glmm_root,
                        phenoname,
                        pheno_exports[phenoname],
                        trait_type,
                        covariates,
                        get_aou_sites_for_grm_path(
                            pop=pop, extension="plink", pruned=True
                        ),
                        # get_aou_saige_plink_path(
                        #     analysis_type=analysis_type,
                        #     pop=pop,
                        #     extension="plink",
                        #     pruned=True,
                        #     data_iteration=iteration,
                        #     window_size=window,
                        #     test=args.test,
                        # ), # TODO: Change this once plink files are ready (only using this for test)
                        SAIGE_DOCKER_IMAGE,
                        sparse_grm=sparse_grm,
                        sparse_grm_extension=sparse_grm_extension,
                        inv_normalize=False,
                        n_threads=n_threads,
                        min_covariate_count=1,  # TODO: check this
                        storage="100Gi",
                    )
                    fit_null_task.attributes.update({"pop": pop})
                    fit_null_task.attributes.update(
                        {"phenotype": copy.deepcopy(phenoname)}
                    )
                    model_file = fit_null_task.null_glmm.rda
                    variance_ratio_file = fit_null_task.null_glmm[
                        f"{analysis_type}.varianceRatio.txt"
                    ]
                    if analysis_type == "gene":
                        sparse_sigma_file = fit_null_task.null_glmm[
                            f'{analysis_type}.varianceRatio.txt{sparse_grm_extension.replace("GRM", "Sigma")}'
                        ]
                        null_models[phenoname] = (
                            model_file,
                            variance_ratio_file,
                            sparse_sigma_file,
                        )
                    else:
                        null_models[phenoname] = (model_file, variance_ratio_file)

            bgen_dir = get_aou_saige_root(
                analysis_type=analysis_type, pop=pop, name="bgen", test=args.test
            )
            overwrite_bgens = args.create_bgens
            bgens_already_created = {}
            if not overwrite_bgens and hl.hadoop_exists(bgen_dir):
                bgens_already_created = {x["path"] for x in hl.hadoop_ls(bgen_dir)}
            # logger.info(f'Found {len(vcfs_already_created)} VCFs in directory...')
            chunk_size = int(5e6) if pop != "eur" else int(1e6)
            mt_path = EXOME_MT_PATH if analysis_type == "gene" else ACAF_MT_PATH
            print(
                f"------------{pop.upper()} MT for {analysis_type} analysis bgen files: {mt_path}------------"
            )

            bgens = {}
            for chrom in chromosomes:
                chromosome = f"chr{chrom}"
                chrom_length = chrom_lengths[chromosome]
                if analysis_type == "gene":
                    if args.groups is None:
                        groups = ["pLoF", "missense|LC", "synonymous"]
                    else:
                        groups = args.groups.split(",")
                    interval_ht = hl.read_table(GENE_INTERVAL_PATH)
                    intervals = interval_ht.aggregate(
                        hl.agg.collect(interval_ht.interval)
                    )
                    for interval in intervals:
                        bgen_root = f"{bgen_dir}/{analysis_type}_{interval.start.contig}_{str(interval.start.position).zfill(9)}_{interval.end.position}"
                        if f"{bgen_root}.bgen" not in bgens_already_created:
                            bgen_task = b.new_python_job(
                                name=f"{analysis_type}_analysis_export_{interval}_bgen"
                            )
                            bgen_task.image(HAIL_DOCKER_IMAGE)
                            bgen_task.call(
                                export_bgen_from_mt,
                                mt_path=mt_path,
                                pop=pop,
                                analysis_type=analysis_type,
                                interval=interval,
                                output_dir=bgen_dir,
                                mean_impute_missing=args.mean_impute_missing,
                                no_adj=args.no_adj,
                                callrate_filter=args.callrate_filter,
                            )
                            bgen_task.attributes["pop"] = pop
                        if f"{bgen_root}.gene.txt" not in bgens_already_created:
                            gene_txt_task = b.new_python_job(
                                name=f"{analysis_type}_analysis_export_{interval}_gene_txt"
                            )
                            gene_txt_task.image(HAIL_DOCKER_IMAGE)
                            gene_txt_task.call(
                                export_bgen_from_mt,
                                interval=interval,
                                groups=groups,
                                outdir=f"{bgen_root}.gene.txt",
                            )
                            gene_txt_task.attributes["pop"] = pop

                        bgen_file = b.read_input_group(
                            **{
                                "bgen": f"{bgen_root}.bgen",
                                "bgen.bgi": f"{bgen_root}.bgen.bgi",
                                "sample": f"{bgen_root}.sample",
                            }
                        )
                        group_file = b.read_input(f"{bgen_root}.gene.txt")
                        bgens[interval] = (bgen_file, group_file)
                        if args.test:
                            break
                else:
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
                        bgen_root = f"{bgen_dir}/{analysis_type}_{interval.start.contig}_{str(interval.start.position).zfill(9)}_{interval.end.position}"
                        if f"{bgen_root}.bgen" not in bgens_already_created:
                            bgen_task = b.new_python_job(
                                name=f"{analysis_type}_analysis_export_{interval}_bgen"
                            )
                            bgen_task.image(HAIL_DOCKER_IMAGE)
                            bgen_task.call(
                                export_bgen_from_mt,
                                mt_path=mt_path,
                                pop=pop,
                                analysis_type=analysis_type,
                                interval=interval,
                                output_dir=bgen_dir,
                                mean_impute_missing=args.mean_impute_missing,
                                no_adj=args.no_adj,
                                callrate_filter=args.callrate_filter,
                            )
                            bgen_task.attributes["pop"] = pop

                        bgen_file = b.read_input_group(
                            **{
                                "bgen": f"{bgen_root}.bgen",
                                "bgen.bgi": f"{bgen_root}.bgen.bgi",
                                "sample": f"{bgen_root}.sample",
                            }
                        )
                        bgens[interval] = bgen_file
                        if args.test:
                            break
                if args.test:
                    break
            b.run()

            result_dir = get_aou_saige_root(
                analysis_type=analysis_type, pop=pop, name="result", test=args.test
            )
            overwrite_results = args.overwrite_results
            for i in range(len(phenos_to_run)):
                phenoname = phenos_to_run[i]
                if phenoname not in null_models.keys():
                    continue
                if analysis_type == "gene":
                    model_file, variance_ratio_file, sparse_sigma_file = null_models[
                        phenoname
                    ]
                else:
                    model_file, variance_ratio_file = null_models[phenoname]

                if not i % 10:
                    n_jobs = dict(
                        Counter(map(lambda x: x.name, b.select_jobs("")))
                    ).get("run_saige", 0)
                    logger.info(f"Read {i} phenotypes ({n_jobs} new to run so far)...")

                pheno_results_dir = f"{result_dir}/phenotype_{phenoname}"
                results_already_created = {}

                if (
                    not overwrite_results
                    and not args.skip_saige
                    and hl.hadoop_exists(pheno_results_dir)
                ):
                    results_already_created = {
                        x["path"] for x in hl.hadoop_ls(pheno_results_dir)
                    }

                saige_tasks = []
                for chrom in chromosomes:
                    chromosome = f"chr{chrom}"
                    chrom_length = chrom_lengths[chromosome]
                    if args.skip_saige:
                        break
                    if analysis_type == "gene":
                        for interval in intervals:
                            bgen_file, group_file = bgens[interval]
                            results_path = f"{pheno_results_dir}/result_{phenoname}_{chromosome}_{str(interval.start.position).zfill(9)}"
                            if (
                                overwrite_results
                                or f"{results_path}.single_variant.txt"
                                not in results_already_created
                            ):
                                samples_file = b.read_input(
                                    get_aou_samples_file_path(
                                        analysis_type=analysis_type,
                                        pop=pop,
                                        iteration=iteration,
                                    )
                                )
                                saige_task = run_saige(
                                    b,
                                    results_path,
                                    model_file,
                                    variance_ratio_file,
                                    bgen_file,
                                    samples_file,
                                    SAIGE_DOCKER_IMAGE,
                                    group_file=group_file,
                                    sparse_sigma_file=sparse_sigma_file,
                                    trait_type=trait_type,
                                    use_bgen=True,
                                    chrom=chromosome,
                                    max_maf=RVAS_AF_CUTOFF,
                                    log_pvalue=True,
                                    min_maf=0,
                                    min_mac=1,
                                )
                                saige_task.attributes.update(
                                    {"interval": interval, "pop": pop}
                                )
                                saige_task.attributes.update(
                                    {"phenotype": copy.deepcopy(phenoname)}
                                )
                                saige_tasks.append(saige_task)
                    else:
                        for start_pos in range(1, chrom_length, chunk_size):
                            end_pos = (
                                chrom_length
                                if start_pos + chunk_size > chrom_length
                                else (start_pos + chunk_size)
                            )
                            interval = f"{chromosome}:{start_pos}-{end_pos}"
                            bgen_file = bgens[interval]
                            results_path = f"{pheno_results_dir}/result_{phenoname}_{chromosome}_{str(start_pos).zfill(9)}"
                            if (
                                overwrite_results
                                or f"{results_path}.single_variant.txt"
                                not in results_already_created
                            ):
                                samples_file = b.read_input(
                                    get_aou_samples_file_path(
                                        analysis_type=analysis_type,
                                        pop=pop,
                                        iteration=iteration,
                                    )
                                )
                                saige_task = run_saige(
                                    b,
                                    results_path,
                                    model_file,
                                    variance_ratio_file,
                                    bgen_file,
                                    samples_file,
                                    SAIGE_DOCKER_IMAGE,
                                    trait_type=trait_type,
                                    use_bgen=True,
                                    chrom=chromosome,
                                    log_pvalue=True,
                                )
                                saige_task.attributes.update(
                                    {"interval": interval, "pop": pop}
                                )
                                saige_task.attributes.update(
                                    {"phenotype": copy.deepcopy(phenoname)}
                                )
                                saige_tasks.append(saige_task)
                        if args.local_test:
                            break
                    if args.local_test:
                        break

            res_tasks = []
            if (
                overwrite_results
                or args.overwrite_hail_results
                or f"{pheno_results_dir}/variant_results.ht"
                not in results_already_created
                or not hl.hadoop_exists(
                    f"{pheno_results_dir}/variant_results.ht/_SUCCESS"
                )
            ):
                null_glmm_root = (
                    f"{null_model_dir}/phenotype_{phenoname}.{analysis_type}.log"
                )

                prefix = f"{pheno_results_dir}/result_{phenoname}_chr{{chrom}}_{str(1).zfill(9)}"
                saige_log = f"{prefix}.{analysis_type}.log"

                load_task = load_results_into_hail(
                    b,
                    pheno_results_dir,
                    phenoname,
                    saige_tasks,
                    get_aou_util_path("vep_corrected"),
                    HAIL_DOCKER_IMAGE,
                    saige_log=saige_log,
                    analysis_type=analysis_type,
                    n_threads=n_threads,
                    null_glmm_log=null_glmm_root,
                    reference=reference,
                    legacy_annotations=True,
                    log_pvalue=True,
                )
                load_task.attributes["pop"] = pop
                res_tasks.append(load_task)
                qq_export, qq_plot = qq_plot_results(
                    b,
                    pheno_results_dir,
                    res_tasks,
                    HAIL_DOCKER_IMAGE,
                    QQ_DOCKER_IMAGE,
                    n_threads=n_threads,
                )
                qq_export.attributes.update({"pop": pop})
                qq_export.attributes.update({"phenotype": copy.deepcopy(phenoname)})
                qq_plot.attributes.update({"pop": pop})
                qq_plot.attributes.update({"phenotype": copy.deepcopy(phenoname)})

        def get_tasks_from_pipeline(p):
            return dict(Counter(map(lambda x: x.name, p.select_jobs(""))))

        logger.info(
            f'Setup took: {time.strftime("%H:%M:%S", time.gmtime(time.time() - start_time))}'
        )
        logger.info(f"Submitting: {get_tasks_from_pipeline(b)}")
        logger.info(f"Total size: {sum([len(x._pretty()) for x in b.select_jobs('')])}")
        logger.info(f"Finished: {get_tasks_from_pipeline(b)}")
        b.run()


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
        "--pilot", help="Run pilot phenotypes only", action="store_true"
    )
    parser.add_argument(
        "--random_pheno", help="Use random phenotypes", action="store_true"
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
