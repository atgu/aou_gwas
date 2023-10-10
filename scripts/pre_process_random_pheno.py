#!/usr/bin/env python3

__author__ = "wlu"

# from aou_gwas import *  # for dataproc
from utils.utils import *  # for QoB
from utils.resources import *  # for QoB
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
def filter_ht_for_grm(
    ht: hl.Table,
    pop: str,
    n_variants_to_keep: int,
    min_call_rate: float = 0.95,
    min_maf: float = 0.01,
):
    ht = ht.filter(
        (ht.locus.in_autosome())
        & (ht.values[f"gvs_{pop}_an"][0] >= (N_SAMPLES[pop] * 2 * min_call_rate))
        & (ht.values[f"gvs_{pop}_ac"][0] > 0)
    )

    sampled_variants = ht.aggregate(
        hl.agg.filter(
            ht.values[f"gvs_{pop}_af"][0] > min_maf,
            hl.agg._reservoir_sample(ht.key, n_variants_to_keep),
        )
    )

    variants = [variant for variant in sampled_variants]
    print(f"N variants sampled: {len(variants)}")
    ht = hl.Table.parallelize(variants).key_by(*ht.key.keys())

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
    memory: str = "37.5Gi",
    storage="1500Mi",
):
    in_bfile = p.read_input_group(
        **{ext: f"{plink_file_root}.{ext}" for ext in ("bed", "bim", "fam")}
    )
    create_sparse_grm_task: Job = p.new_job(name="create_sparse_grm")
    create_sparse_grm_task.cpu(n_threads).storage(storage).image(docker_image).memory(
        memory
    )
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
        worker_memory="highmem",
        worker_cores=1,
        default_reference="GRCh38",
        log="/pre_process_random_pheno.log",
    )

    pops = args.pop.split(",")
    ###############
    # Test chunk
    if args.test:
        print(pops[0])
        ht = hl.read_table(get_aou_util_path(name="vat"))
        ht = ht.collect_by_key()
        filtered_ht = filter_ht_for_grm(
            ht, pop=pops[0], min_call_rate=0.9, n_variants_to_keep=N_SAMPLES[pops[0]]
        )
        filtered_ht.describe()
    ###############
    for pop in pops:
        if args.create_plink_file:
            if (
                not hl.hadoop_exists(
                    get_aou_sites_for_grm_path(pop=pop, extension="ht")
                )
                # or args.overwrite # only needed when the final markers keeped for GRM is not enough
            ):
                print(
                    f"-------------Exporting downsampled variant HT (pop: {pop})-----------"
                )
                n_variants_to_keep = 150000
                ht = hl.read_table(get_aou_util_path(name="vat"))
                ht = ht.collect_by_key()
                filtered_ht = filter_ht_for_grm(
                    ht,
                    pop=pop,
                    min_call_rate=MIN_CALL_RATE[pop],
                    n_variants_to_keep=n_variants_to_keep,
                )
                # Note:
                # 1) It is ideal to have N_variant == N_sample when building GRM
                # 2) It is sufficient to use just variants with MAF >= 0.01
                filtered_ht.naive_coalesce(1000).checkpoint(
                    get_aou_sites_for_grm_path(pop=pop, extension="ht"),
                    _read_if_exists=not args.overwrite,
                    overwrite=args.overwrite,
                )
            if (
                not hl.hadoop_exists(
                    get_aou_sites_for_grm_path(pop=pop, extension="mt")
                )
                # or args.overwrite
            ):
                filtered_ht = hl.read_table(
                    get_aou_sites_for_grm_path(pop=pop, extension="ht")
                )

                print(
                    f"-------------Exporting downsampled variant MT (pop: {pop})-----------"
                )
                pop_ht = hl.read_table(
                    get_aou_util_path(name="ancestry_preds", parsed=True)
                )
                mt = hl.read_matrix_table(EXOME_MT_PATH)  # (34807589, 245394)
                mt = mt.filter_rows(hl.is_defined(filtered_ht[mt.row_key]))

                duplicated_samples = hl.read_table(
                    get_aou_relatedness_path(extension="1st_degrees.ht")
                )
                print(
                    f"-------------Removing {duplicated_samples.count()} potentially duplicated samples-----------"
                )
                mt = mt.filter_cols(hl.is_missing(duplicated_samples[mt.col_key]))

                if pop != "all":
                    pop_ht = pop_ht.filter(pop_ht.ancestry_pred == pop)
                    print(f"N samples kept: {pop_ht.count()}")
                    mt = mt.filter_cols(hl.is_defined(pop_ht[mt.s]))
                mt = mt.annotate_rows(call_stats=hl.agg.call_stats(mt.GT, mt.alleles))
                mt.describe()

                mt = mt.naive_coalesce(1000).checkpoint(
                    get_aou_sites_for_grm_path(pop=pop, extension="mt"),
                    _read_if_exists=not args.overwrite,
                    overwrite=args.overwrite,
                )

            if args.ld_prune:
                mt = hl.read_matrix_table(
                    get_aou_sites_for_grm_path(pop=pop, extension="mt")
                )

                mt = mt.unfilter_entries()
                if (
                    not hl.hadoop_exists(
                        get_aou_sites_for_grm_path(pop=pop, extension="ht", pruned=True)
                    )
                ) or args.overwrite:
                    print(
                        f"-------------Exporting the LD pruned downsampled variant HT (pop: {pop})-----------"
                    )
                    ht = hl.ld_prune(mt.GT, r2=0.1)
                    ht.checkpoint(
                        get_aou_sites_for_grm_path(
                            pop=pop, extension="ht", pruned=True
                        ),
                        _read_if_exists=not args.overwrite,
                        overwrite=args.overwrite,
                    )
                ht = hl.read_table(
                    get_aou_sites_for_grm_path(pop=pop, extension="ht", pruned=True)
                )
                mt = mt.filter_rows(hl.is_defined(ht[mt.row_key]))

            if args.overwrite or not hl.hadoop_exists(
                f'{get_aou_sites_for_grm_path(pop = pop, extension="bed", pruned=args.ld_prune)}'
            ):
                print(
                    f"-------------Exporting variant downsampled plink files (pop: {pop})-----------"
                )
                hl.export_plink(
                    mt,
                    get_aou_sites_for_grm_path(
                        pop=pop, extension="plink", pruned=args.ld_prune
                    ),
                )

        if args.create_sparse_grm:
            sparse_grm_root = f"{DATA_PATH}/utils/grm/aou_{pop}"
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

            create_sparse_grm(
                b,
                sparse_grm_root,
                get_aou_sites_for_grm_path(
                    pop=pop, extension="plink", pruned=args.ld_prune
                ),
                SAIGE_DOCKER_IMAGE,
                relatedness_cutoff,
                num_markers,
                n_threads=n_threads,
            )

            b.run()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--create-plink-file", help="Create plink files for GRM", action="store_true"
    )
    parser.add_argument(
        "--create-sparse-grm", help="Create the sparse grm", action="store_true"
    )
    parser.add_argument(
        "--ld-prune",
        help="Whether to run LD pruning on the downsampled variants",
        action="store_true",
    )
    parser.add_argument(
        "--pop", help="Comma-separated list of pops to run", default="all"
    )
    parser.add_argument(
        "--overwrite",
        help="Whether to overwrite the existing file",
        action="store_true",
    )
    parser.add_argument(
        "--test",
        help="Whether to run test chunk",
        action="store_true",
    )
    args = parser.parse_args()

    main(args)
