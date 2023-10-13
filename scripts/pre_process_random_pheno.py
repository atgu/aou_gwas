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
def mac_category_case_builder(call_stats_ac_expr, call_stats_af_expr, min_maf_common_variants: float = 0.01):
    return (
        hl.case()
        .when(call_stats_ac_expr <= 5, call_stats_ac_expr)
        .when(call_stats_ac_expr <= 10, 10)
        .when(call_stats_ac_expr <= 20, 20)
        .when(call_stats_af_expr <= 0.001, 0.001)
        .when(call_stats_af_expr <= min_maf_common_variants, min_maf_common_variants)
        .default(0.99)
    )

def filter_ht_for_grm(
    ht: hl.Table, # VAT format
    pop: str,
    n_common_variants_to_keep: int=150000, # NOTE: to ensure sufficient number of variants for GRM
    min_call_rate: float = 0.95,
    min_maf_common_variants: float = 0.01,
    variants_per_mac_category: int = 2000,
    variants_per_maf_category: int = 10000,
):
    ht = ht.filter(
        (ht.locus.in_autosome())
        & (ht.values[f"gvs_{pop}_an"][0] >= (N_SAMPLES[pop] * 2 * min_call_rate))
        & (ht.values[f"gvs_{pop}_ac"][0] > 0)
    )

    ht = ht.annotate(
        mac_category=mac_category_case_builder(ht.values[f"gvs_{pop}_ac"][0], ht.values[f"gvs_{pop}_af"][0], min_maf_common_variants)
    )

    # From: https://hail.zulipchat.com/#narrow/stream/123010-Hail-Query-0.2E2-support/topic/.E2.9C.94.20randomly.20sample.20table/near/388162012
    bins = ht.aggregate(hl.agg.collect_as_set(ht.mac_category))
    ac_bins = [bin for bin in bins if (bin >= 1) and (bin<=5)]
    af_bins = [bin for bin in bins if (bin < 0.99) or (bin > 5)]

    sampled_common_variants = ht.aggregate(
        hl.agg.filter(
                ht.values[f"gvs_{pop}_af"][0] > min_maf_common_variants,
                hl.agg._reservoir_sample(ht.key, n_common_variants_to_keep),
            ),
    )
    print('Finished sampling common variants...')
    common_variants = [variant for variant in sampled_common_variants]

    binned_variants_af = ht.aggregate(
        hl.agg.array_agg(
            lambda x: hl.agg.filter(
                ht.mac_category == x,
                hl.agg._reservoir_sample(ht.key, variants_per_maf_category),
            ),
            hl.literal(af_bins),
        ),
    )
    print('Finished sampling rare variants...')

    binned_variants_ac = ht.aggregate(
        hl.agg.array_agg(
            lambda x: hl.agg.filter(
                ht.mac_category == x,
                hl.agg._reservoir_sample(ht.key, variants_per_mac_category),
            ),
            hl.literal(ac_bins),
        )
    )
    print('Finished sampling ultra-rare variants...')

    binned_rare_variants = binned_variants_ac + binned_variants_af
    rare_variants = [variant for bin in binned_rare_variants for variant in bin]

    print(f"N rare variants sampled: {len(rare_variants)}")
    print(f"N common variants sampled: {len(common_variants)}")
    rare_ht = hl.Table.parallelize(rare_variants).key_by(*ht.key.keys())
    common_ht = hl.Table.parallelize(common_variants).key_by(*ht.key.keys())
    ht = rare_ht.union(common_ht)
    ht.describe()

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
    app_name = None
    if args.create_plink_file:
        app_name = f"create_plink_file_{args.pop.replace(',', '_')}"
    if args.test:
        app_name = f"run_test"
    print(app_name)

    hl.init(
        app_name=f'Creating_plink_{args.pop.replace(",", "_")}',
        # app_name='Sample_QC_on_MT',
        tmp_dir=TMP_BUCKET,
        driver_memory="highmem",
        driver_cores=8,
        worker_memory="highmem",
        worker_cores=1,
        default_reference="GRCh38",
        log="/pre_process_saige_data.log"
    )

    pops = args.pop.split(",")
    ###############
    # Test chunk
    if args.test:
        mt = hl.read_matrix_table(EXOME_MT_PATH)
        mt.filters.show()
        # print(hl.len(mt.filters).summarize())
        # - <expr> (int32):
        #   Non-missing: 34807589 (100.00%)
        #       Missing: 0
        #       Minimum: 0
        #       Maximum: 0
        #          Mean: 0.00
        #       Std Dev: 0.00
    ###############
    for pop in pops:
        if args.create_plink_file:
            if (not hl.hadoop_exists(get_aou_sites_for_grm_path(pop=pop, extension="ht"))
                or args.overwrite_variant_ht
            ):
                print(f"-------------Exporting downsampled variant HT (pop: {pop})-----------")
                print(f'Min call rate for {pop.upper()}: {MIN_CALL_RATE[pop]}')
                ht = hl.read_table(get_aou_util_path(name="vat"))
                ht = ht.collect_by_key()
                filtered_ht = filter_ht_for_grm(
                    ht,
                    pop=pop,
                    n_common_variants_to_keep= 150000,
                    min_call_rate= MIN_CALL_RATE[pop],
                )

                filtered_ht.naive_coalesce(1000).checkpoint(
                    get_aou_sites_for_grm_path(pop=pop, extension="ht"),
                    _read_if_exists=not args.overwrite_variant_ht,
                    overwrite=args.overwrite_variant_ht,
                )
            filtered_ht = hl.read_table(get_aou_sites_for_grm_path(pop=pop, extension="ht"))
            print(f'Number of variants sampled for {pop}: {filtered_ht.count()}')
            if (not hl.hadoop_exists(get_aou_sites_for_grm_path(pop=pop, extension="mt"))
                or args.overwrite_variant_mt
            ):
                filtered_ht = hl.read_table(get_aou_sites_for_grm_path(pop=pop, extension="ht"))

                print(
                    f"-------------Exporting downsampled variant MT (pop: {pop})-----------"
                )
                mt = get_filtered_mt(analysis_type='gene', filter_samples=False, filter_variants=True, adj_filter=True, pop=pop)
                meta_ht = hl.read_table(get_sample_meta_path(annotation=True))
                meta_ht = meta_ht.filter(~meta_ht.related_0th_degree)
                mt = mt.filter_cols(hl.is_defined(meta_ht[mt.col_key]))
                mt = mt.filter_rows(hl.is_defined(filtered_ht[mt.row_key]))

                print("Removing HLA...")
                # Common inversion taken from Table S4 of https://www.ncbi.nlm.nih.gov/pubmed/27472961
                # (converted to GRCh38 by: https://liftover.broadinstitute.org/#input=chr8%3A8055789-11980649&hg=hg19-to-hg38 )
                # Also removing HLA, from https://www.ncbi.nlm.nih.gov/grc/human/regions/MHC?asm=GRCh38
                mt = mt.filter_rows(
                    ~hl.parse_locus_interval(
                        "chr8:8198267-12123140", reference_genome="GRCh38"
                    ).contains(mt.locus)
                    & ~hl.parse_locus_interval(
                        "chr6:28510120-33480577", reference_genome="GRCh38"
                    ).contains(mt.locus)
                )

                mt.naive_coalesce(1000).checkpoint(
                    get_aou_sites_for_grm_path(pop=pop, extension="mt"),
                    _read_if_exists=not args.overwrite_variant_mt,
                    overwrite=args.overwrite_variant_mt,
                )
                mt = hl.read_matrix_table(get_aou_sites_for_grm_path(pop=pop, extension="mt"))
                mt.describe()
                print(mt.count())

            if args.ld_prune:
                mt = hl.read_matrix_table(get_aou_sites_for_grm_path(pop=pop, extension="mt"))
                mt = mt.unfilter_entries()
                if (not hl.hadoop_exists(get_aou_sites_for_grm_path(pop=pop, extension="ht", pruned=args.ld_prune))
                    or args.overwrite_ld_ht
                ):
                    print(f"-------------Exporting the LD pruned downsampled variant HT (pop: {pop})-----------")
                    ht = hl.ld_prune(mt.GT,
                                     r2=0.1,
                                     bp_window_size=1e7,
                                     block_size=1024,
                                     )
                    ht.checkpoint(
                        get_aou_sites_for_grm_path(pop=pop, extension="ht", pruned=args.ld_prune),
                        _read_if_exists=not args.overwrite_ld_ht,
                        overwrite=args.overwrite_ld_ht,
                    )
                ht = hl.read_table(get_aou_sites_for_grm_path(pop=pop, extension="ht", pruned=args.ld_prune))
                mt = mt.filter_rows(hl.is_defined(ht[mt.row_key]))

            if args.overwrite_plink or not hl.hadoop_exists(
                f'{get_aou_sites_for_grm_path(pop = pop, extension="plink.bed", pruned=True)}'
            ):
                print(f"-------------Exporting variant downsampled plink files (pop: {pop})-----------")
                hl.export_plink(
                    mt,
                    get_aou_sites_for_grm_path(pop=pop, extension="plink", pruned=args.ld_prune),
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
        "--overwrite-variant-ht", help="Resample variants for GRM and plink files", action="store_true"
    )
    parser.add_argument(
        "--overwrite-variant-mt", help="Overwrite the MT filtered to selected variants and samples for GRM and plink files", action="store_true"
    )
    parser.add_argument(
        "--overwrite-ld-ht", help="Overwrite LD-pruned variant HT", action="store_true"
    )
    parser.add_argument(
        "--overwrite-plink", help="Overwrite plink files", action="store_true"
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
        "--pop", help="Comma-separated list of pops to run", default='afr,amr,eas,eur,mid,sas'
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
