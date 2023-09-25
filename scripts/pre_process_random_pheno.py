#!/usr/bin/env python3

__author__ = "wlu"

from utils.utils import *
from utils.resources import *
import hail as hl
import argparse
import hailtop.batch as hb
from hailtop.batch.job import Job
from hailtop.batch.batch import Batch

TRANCHE = "250k"
BUCKET = "gs://aou_analysis"
MY_BUCKET = "gs://aou_wlu"
TMP_BUCKET = "gs://aou_tmp"
CURRENT_PATH = f"{BUCKET}/{TRANCHE}"
DATA_PATH = f"{CURRENT_PATH}/data"


def get_aou_sites_for_grm_path(
        pop: str, extension: str, pruned: bool = False, tranche: str = TRANCHE
):
    return f'{DATA_PATH}/utils/grm/aou_{pop}_sites{"_pruned" if pruned else ""}_for_grm_{tranche}.{extension}'


####### STEP 1: Downsample a bunch of sites based on frequency for the GRM  #######
def mac_category_case_builder(call_stats_ac_expr, call_stats_af_expr):
    return (
        hl.case()
        .when(call_stats_ac_expr <= 5, call_stats_ac_expr)
        .when(call_stats_ac_expr <= 10, 10)
        .when(call_stats_ac_expr <= 20, 20)
        .when(call_stats_af_expr <= 0.001, 0.001)
        .when(call_stats_af_expr <= 0.01, 0.01)
        .when(call_stats_af_expr <= 0.1, 0.1)
        .default(0.99)
    )


def filter_ht_for_plink(
    ht: hl.Table,
    pop: str,
    min_call_rate: float = 0.95,
    variants_per_mac_category: int = 2000,
    variants_per_maf_category: int = 10000,
):
    ht = ht.filter(
        (ht.locus.in_autosome())
        & (ht.values[f"gvs_{pop}_an"][0] >= N_SAMPLES[pop] * 2 * min_call_rate)
        & (ht.values[f"gvs_{pop}_ac"][0] > 0)
    )
    ht = ht.annotate(
        mac_category=mac_category_case_builder(
            ht.values[f"gvs_{pop}_ac"][0], ht.values[f"gvs_{pop}_af"][0]
        )
    )

    # From: https://hail.zulipchat.com/#narrow/stream/123010-Hail-Query-0.2E2-support/topic/.E2.9C.94.20randomly.20sample.20table/near/388162012
    bins = ht.aggregate(hl.agg.collect_as_set(ht.mac_category))
    ac_bins = [bin for bin in bins if bin >= 1]
    print(ac_bins)
    af_bins = [bin for bin in bins if bin < 1]
    print(af_bins)

    binned_samples_ac = ht.aggregate(
        hl.agg.array_agg(
            lambda x: hl.agg.filter(
                ht.mac_category == x,
                hl.agg._reservoir_sample(ht.key, variants_per_mac_category),
            ),
            hl.literal(ac_bins),
        )
    )
    binned_samples_af = ht.aggregate(
        hl.agg.array_agg(
            lambda x: hl.agg.filter(
                ht.mac_category == x,
                hl.agg._reservoir_sample(ht.key, variants_per_maf_category),
            ),
            hl.literal(af_bins),
        ),
    )
    binned_samples = binned_samples_ac + binned_samples_af

    samples = [sample for bin in binned_samples for sample in bin]
    ht = hl.Table.parallelize(samples).key_by(*ht.key.keys())

    return ht


####### STEP 2: Create sparse GRM  #######
def create_sparse_grm(
    p: Batch,
    output_path: str,
    plink_file_root: str,
    docker_image: str,
    relatedness_cutoff: str = "0.125",
    num_markers: int = 2000,
    n_threads: int = 8,
    storage="1500Mi",
):
    in_bfile = p.read_input_group(
        **{ext: f"{plink_file_root}.{ext}" for ext in ("bed", "bim", "fam")}
    )
    create_sparse_grm_task: Job = p.new_job(name="create_sparse_grm")
    create_sparse_grm_task.cpu(n_threads).storage(storage).image(docker_image)
    create_sparse_grm_task.declare_resource_group(
        sparse_grm={
            ext: f"{{root}}{ext}"
            for ext in (
                f"_relatednessCutoff_{relatedness_cutoff}_{num_markers}_randomMarkersUsed.sparseGRM.mtx",
                f"_relatednessCutoff_{relatedness_cutoff}_{num_markers}_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt",
            )
        }
    )
    command = (
        f"Rscript /usr/local/bin/createSparseGRM.R "
        f"--plinkFile={in_bfile} "
        f"--nThreads={n_threads} "
        f"--outputPrefix={create_sparse_grm_task.sparse_grm}	"
        f"--numRandomMarkerforSparseKin={num_markers} "
        f"--relatednessCutoff={relatedness_cutoff}"
    )
    create_sparse_grm_task.command(command)
    p.write_output(create_sparse_grm_task.sparse_grm, output_path)
    # Runtime: ~40 minutes on 8 cores
    return create_sparse_grm_task.sparse_grm


def main(args):
    hl.init(
        driver_memory="highmem",
        driver_cores=8,
        default_reference="GRCh38",
        #  log="/pre_process_random_pheno.log",
    )

    pops = args.pop.split(",")
    for pop in pops:
        if args.create_plink_file:
            if (
                not hl.hadoop_exists(
                    get_aou_sites_for_grm_path(pop=pop, extension="ht")
                )
                or args.overwrite
            ):
                print(
                    f"-------------Exporting downsampled variant HT (pop: {pop})-----------"
                )
                ht = hl.read_table(get_aou_util_path(name="vat"))
                ht = ht.collect_by_key()
                filtered_ht = filter_ht_for_plink(ht, pop=pop)
                filtered_ht.naive_coalesce(1000).checkpoint(
                    get_aou_sites_for_grm_path(pop=pop, extension="ht"),
                    _read_if_exists=not args.overwrite,
                    overwrite=args.overwrite,
                )
            filtered_ht = hl.read_table(
                get_aou_sites_for_grm_path(pop=pop, extension="ht")
            )
            filtered_ht.describe()
            print(f"N variants: {filtered_ht.count()}")  # 54,000

            print(f"-------------Exporting downsampled variant MT (pop: {pop})-----------")
            mt = hl.read_matrix_table(EXOME_MT_PATH)  # (34807589, 245394)
            mt = mt.filter_rows(hl.is_defined(filtered_ht[mt.row_key]))
            mt = mt.naive_coalesce(1000).checkpoint(
                get_aou_sites_for_grm_path(pop=pop, extension="mt"),
                _read_if_exists=not args.overwrite,
                overwrite=args.overwrite,
            )

            if args.ld_prune:
                print(f"-------------Exporting the LD pruned downsampled variant HT (pop: {pop})-----------")
                mt = mt.unfilter_entries()
                ht = hl.ld_prune(mt.GT, r2=0.1)

                ht = ht.checkpoint(
                    get_aou_sites_for_grm_path(pop=pop, extension="ht", pruned=True),
                    _read_if_exists=not args.overwrite,
                    overwrite=args.overwrite,
                )
                mt = mt.filter_rows(hl.is_defined(ht[mt.row_key]))

            if args.overwrite or not hl.hadoop_exists(
                f'{get_aou_sites_for_grm_path(pop = pop, extension="bed", pruned=args.ld_prune)}'
            ):
                print(f"-------------Exporting variant downsampled plink files (pop: {pop})-----------")
                hl.export_plink(
                    mt,
                    get_aou_sites_for_grm_path(pop=pop, extension="plink", pruned=args.ld_prune),
                )

        if args.create_sparse_grm:
            sparse_grm_root = f"{DATA_PATH}/utils/grm"
            SAIGE_DOCKER_IMAGE = "wzhou88/saige:1.3.0"
            relatedness_cutoff = "0.125"
            num_markers = 2000
            n_threads = 8


            backend = hb.ServiceBackend(
                billing_project="all-by-aou",
                remote_tmpdir=TMP_BUCKET,
            )
            b = hb.Batch(
                name=f"Create_sparse_GRM",
                requester_pays_project="aou-neale-gwas",
                backend=backend,
            )



            sparse_grm = create_sparse_grm(
                b,
                sparse_grm_root,
                get_aou_sites_for_grm_path(pop=pop, extension="plink", pruned=args.ld_prune),
                SAIGE_DOCKER_IMAGE,
                relatedness_cutoff,
                num_markers,
                n_threads=n_threads,
            )

            b.run()

    if args.run_hail_ibd:
        MAF_CUTOFF = 0.05
        mt = hl.read_matrix_table(ACAF_MT_PATH)  # (99250816, 245394)
        mt.describe()
        # print(mt.info.AF.summarize())
        # # - AF (array<float64>):
        # #   Non-missing: 99250796 (100.00%)
        # #       Missing: 20 (0.00%)
        # #      Min Size: 1
        # #      Max Size: 1
        # #     Mean Size: 1.00
        # #
        # #   - AF[<elements>] (float64):
        # #     Non-missing: 99250796 (100.00%)
        # #         Missing: 0
        # #         Minimum: 0.00
        # #         Maximum: 1.00
        # #            Mean: 0.04
        # #         Std Dev: 0.13
        # print(mt.aggregate_rows(hl.agg.count_where(mt.info.AF[0] > MAF_CUTOFF))) # 10289329
        ### Filter to 1M random common variants before running IBD

        if (not hl.hadoop_exists(get_aou_relatedness_path(extension="1Mvar.ht"))) or args.overwrite:
            print(f"-------------Downsampling to 1M common variants-----------")
            ht = mt.rows()
            sampled_variants = ht.aggregate(
                hl.agg.filter(
                    ht.info.AF[0] >= MAF_CUTOFF, hl.agg._reservoir_sample(ht.key, 1000000)
                )
            )
            variants = [variant for variant in sampled_variants]
            ht = hl.Table.parallelize(variants).key_by(*ht.key.keys())
            ht = ht.checkpoint(
                get_aou_relatedness_path(extension="1Mvar.ht"), _read_if_exists=(not args.overwrite), overwrite=args.overwrite
            )
        ht = hl.read_table(get_aou_relatedness_path(extension="1Mvar.ht"))
        #print(f'N variants: {ht.count()}')
        mt = mt.filter_rows(hl.is_defined(ht[mt.row_key]))
        ### Run IBD
        print(f"-------------Running hl.identity_by_descent() -----------")
        relatedness_ht = hl.identity_by_descent(mt, maf=mt.info.AF[0])
        if args.overwrite or (not hl.hadoop_exists(
            f'{get_aou_relatedness_path(extension="ht")}'
        )):
            print(f"-------------Writing AoU IBD HT -----------")
            relatedness_ht.write(get_aou_relatedness_path(extension="ht"), args.overwrite)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--create-plink-file", help="Create plink files for GRM", action="store_true"
    )
    parser.add_argument(
        "--create-sparse-grm", help="Create the sparse grm", action="store_true"
    )
    parser.add_argument(
        "--run-hail-ibd",
        help="Run hl.identity_by_descent() to obtain relatedness information",
        action="store_true",
    )
    parser.add_argument(
        "--ld-prune",
        help="Whether to run LD pruning on the downsampled variants",
        action="store_true",
    )
    parser.add_argument("--pop", help="Comma-separated list of pops to run")
    parser.add_argument(
        "--overwrite",
        help="Whether to overwrite the existing file",
        action="store_true",
    )
    args = parser.parse_args()

    main(args)
