import hail as hl
import argparse
import hailtop.batch as hb
from hailtop.batch.job import Job
from hailtop.batch.batch import Batch

tmp_dir = 'gs://aou_tmp'
root = 'gs://aou_wlu'
vat_root = 'gs://aou_wlu/utils/vat'

ACAF_MT_PATH = 'gs://prod-drc-broad/aou-wgs-delta-small_callsets_gq0/v7.1/acaf_threshold_v7.1/splitMT/delta_basis_without_ext_aian_prod_gq0_3regions.acaf_threshold.split.mt'
EXOME_MT_PATH = 'gs://prod-drc-broad/aou-wgs-delta-small_callsets_gq0/v7.1/exome_v7.1/splitMT/delta_basis_without_ext_aian_prod_gq0_3regions.exome.split.mt'
VDS_PATH = 'gs://prod-drc-broad/v7/wgs/without_ext_aian_prod/vds/aou_srwgs_short_variants_v7_without_ext_aian_prod.vds'
VAT_HT_PATH = f'{vat_root}/aou_PARSED_SORTED_vat.ht'
##################### Computed in the AoU workbench by running:
# ANCESTRY_PATH = 'gs://fc-aou-datasets-controlled/v7/wgs/short_read/snpindel/aux/ancestry/ancestry_preds.tsv'
# pop_ht = hl.import_table(ANCESTRY_PATH, key="research_id", impute=True,
#                              types={"research_id": "tstr", "pca_features": hl.tarray(hl.tfloat)})
# pop_ht.aggregate(hl.agg.counter(pop_ht.ancestry_pred))
N_SAMPLES = {'all':245394, 'afr':56913, 'amr':45035, 'eas':5706, 'eur':133581, 'mid':942, 'sas':3217}
###############################################################

def get_aou_sites_for_grm_path(extension: str, pruned: bool=False):
    return f'{root}/utils/aou_sites{"_pruned" if pruned else ""}_for_grm.{extension}'
def get_aou_relatedness_path(extension: str = 'ht'):
    return f'{root}/utils/aou_ibd_relatedness.{extension}'


####### STEP 1: Downsample a bunch of sites based on frequency for the GRM  #######
def mac_category_case_builder(call_stats_ac_expr, call_stats_af_expr):
    return (hl.case()
            .when(call_stats_ac_expr <= 5, call_stats_ac_expr)
            .when(call_stats_ac_expr <= 10, 10)
            .when(call_stats_ac_expr <= 20, 20)
            .when(call_stats_af_expr <= 0.001, 0.001)
            .when(call_stats_af_expr<= 0.01, 0.01)
            .when(call_stats_af_expr<= 0.1, 0.1)
            .default(0.99))

def filter_ht_for_plink(
        ht: hl.Table,
        pop:str,
        min_call_rate: float = 0.95,
        variants_per_mac_category: int = 2000,
        variants_per_maf_category: int = 10000
):
    ht = ht.filter((ht.locus.in_autosome()) &
                   (ht[f'gvs_{pop}_an'][0] >= N_SAMPLES[pop] * 2 * min_call_rate) &
                   (ht[f'gvs_{pop}_ac'][0] > 0))
    ht = ht.annotate(mac_category=mac_category_case_builder(ht[f'gvs_{pop}_ac'][0], ht[f'gvs_{pop}_af'][0]))

    # From: https://hail.zulipchat.com/#narrow/stream/123010-Hail-Query-0.2E2-support/topic/.E2.9C.94.20randomly.20sample.20table/near/388162012
    bins = ht.aggregate(hl.agg.collect_as_set(ht.mac_category))
    ac_bins = [bin for bin in bins if bin >= 1]
    af_bins = [bin for bin in bins if bin < 1]
    binned_samples = ht.aggregate(
        hl.if_else(ht.mac_category >= 1,
                   hl.agg.array_agg(lambda x: hl.agg.filter(ht.mac_category == x,
                                                            hl.agg._reservoir_sample(ht.key,
                                                                                     variants_per_mac_category)), ac_bins),
                   hl.agg.array_agg(lambda x: hl.agg.filter(ht.mac_category == x,
                                                            hl.agg._reservoir_sample(ht.key,
                                                                                     variants_per_maf_category)),
                                    af_bins)
                   )
    )
    samples = [sample for bin in binned_samples for sample in bin]
    ht = hl.Table.parallelize(samples).key_by(*ht.key.keys())

    category_counter = ht.aggregate(hl.agg.counter(ht.mac_category))
    print(category_counter)
    ht = ht.annotate_globals(category_counter=category_counter)
    return ht

####### STEP 2: Create sparse GRM  #######
def create_sparse_grm(p: Batch, output_path: str, plink_file_root: str, docker_image: str,
                      relatedness_cutoff: str = '0.125', num_markers: int = 2000,
                      n_threads: int = 8, storage = '1500Mi'):
    in_bfile = p.read_input_group(
        **{ext: f'{plink_file_root}.{ext}' for ext in ('bed', 'bim', 'fam')})
    create_sparse_grm_task: Job = p.new_job(name='create_sparse_grm')
    create_sparse_grm_task.cpu(n_threads).storage(storage).image(docker_image)
    create_sparse_grm_task.declare_resource_group(
        sparse_grm={ext: f'{{root}}{ext}' for ext in
                   (f'_relatednessCutoff_{relatedness_cutoff}_{num_markers}_randomMarkersUsed.sparseGRM.mtx',
                    f'_relatednessCutoff_{relatedness_cutoff}_{num_markers}_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt')})
    command = (f'Rscript /usr/local/bin/createSparseGRM.R '
               f'--plinkFile={in_bfile} '
               f'--nThreads={n_threads} '
               f'--outputPrefix={create_sparse_grm_task.sparse_grm}	'
               f'--numRandomMarkerforSparseKin={num_markers} '
               f'--relatednessCutoff={relatedness_cutoff}')
    create_sparse_grm_task.command(command)
    p.write_output(create_sparse_grm_task.sparse_grm, output_path)
    # Runtime: ~40 minutes on 8 cores
    return create_sparse_grm_task.sparse_grm


def main(args):
    hl.init(driver_memory='highmem', driver_cores=8, default_reference='GRCh38', log='/pre_process_random_pheno.log')

    if args.create_plink_file:
        mt = hl.read_matrix_table(EXOME_MT_PATH)  # (34807589, 245394)
        ht = hl.read_table(VAT_HT_PATH)
        ht = ht.collect_by_key()
        filtered_ht = filter_ht_for_plink(ht, pop=args.pop)
        filtered_ht = filtered_ht.naive_coalesce(1000).checkpoint(get_aou_sites_for_grm_path(extension='ht'),
                                                                  _read_if_exists=not args.overwrite,
                                                                  overwrite=args.overwrite)

        mt = mt.filter_rows(hl.is_defined(filtered_ht[mt.row_key]))
        mt = mt.naive_coalesce(1000).checkpoint(get_aou_sites_for_grm_path(extension='mt'),
                                                _read_if_exists=not args.overwrite,
                                                overwrite=args.overwrite)

        if args.ld_prune: # wei: it should be fine to skip LD pruning here
            mt = mt.unfilter_entries()
            ht = hl.ld_prune(mt.GT, r2=0.1)

            ht = ht.checkpoint(get_aou_sites_for_grm_path(extension='ht', pruned=True),
                               _read_if_exists=not args.overwrite,
                               overwrite=args.overwrite)
            mt = mt.filter_rows(hl.is_defined(ht[mt.row_key]))

        if args.overwrite or not hl.hadoop_exists(f'{get_aou_sites_for_grm_path(extension="bed", pruned=True)}'):
            hl.export_plink(mt, get_aou_sites_for_grm_path(extension="plink", pruned=True))

    if args.create_sparse_grm:
        backend = hb.ServiceBackend(
            billing_project="all-by-aou",
            remote_tmpdir=tmp_dir,
        )
        b = hb.Batch(
            name=f"Create_sparse_GRM",
            requester_pays_project='aou-neale-gwas',
            default_image='wzhou88/saige:1.2.0',
            backend=backend,
        )
        sparse_grm_root = f'{root}/utils/grm/sparse'
        SAIGE_DOCKER_IMAGE = 'wzhou88/saige:1.2.0'
        relatedness_cutoff = '0.125'
        num_markers = 2000
        n_threads = 8

        sparse_grm = create_sparse_grm(b,
                                       sparse_grm_root,
                                       get_aou_sites_for_grm_path(extension="plink", pruned=False),
                                       SAIGE_DOCKER_IMAGE,
                                       relatedness_cutoff,
                                       num_markers,
                                       n_threads=n_threads
                                       )
    if args.run_hail_ibd:
        MAF_CUTOFF = 0.05
        mt = hl.read_matrix_table(ACAF_MT_PATH) # (99250816, 245394)
        mt.describe()
        print(mt.count())
        print(mt.info.AF.summarize())
        print(mt.aggregate_rows(hl.agg.count_where(mt.info.AF > 0.05)))
        ### Filter to 1M random common variants before running IBD
        ht = mt.rows()
        ht = ht.aggregate(hl.agg.filter(ht.info.AF >= MAF_CUTOFF, hl.agg._reservoir_sample(ht.key, 1000000)))
        ht = ht.checkpoint(get_aou_relatedness_path(extension="1Mvar.ht"))
        mt = mt.filter_rows(hl.is_defined(ht[mt.row_key]))
        ### Run IBD
        relatedness_ht = hl.identity_by_descent(mt, maf=mt.info.AF)
        if args.overwrite or not hl.hadoop_exists(f'{get_aou_relatedness_path(extension="ht")}'):
            relatedness_ht.write(get_aou_relatedness_path(extension="ht"))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--create-plink-file",
        help="Create plink files for GRM",
        action="store_true"
    )
    parser.add_argument(
        "--create-sparse-grm",
        help="Create the sparse grm",
        action="store_true"
    )
    parser.add_argument(
        "--run-hail-ibd",
        help="Run hl.identity_by_descent() to obtain relatedness information",
        action="store_true"
    )
    parser.add_argument(
        "--ld-prune",
        help="Whether to run LD pruning on the downsampled variants",
        action="store_true"
    )
    parser.add_argument('--pop', help='Comma-separated list of pops to run',
                        choices=['afr', 'amr', 'eas', 'eur', 'mid', 'sas', 'all'])
    parser.add_argument('--overwrite', help='Whether to overwrite the existing file', action='store_true')
    args = parser.parse_args()

    main(args)
