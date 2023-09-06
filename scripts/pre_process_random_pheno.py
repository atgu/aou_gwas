import hail as hl
import argparse
import hailtop.batch as hb

root = 'gs://aou_wlu'
# Computed in the AoU workbench by running:
# ANCESTRY_PATH = 'gs://fc-aou-datasets-controlled/v7/wgs/short_read/snpindel/aux/ancestry/ancestry_preds.tsv'
# pop_ht = hl.import_table(ANCESTRY_PATH, key="research_id", impute=True,
#                              types={"research_id": "tstr", "pca_features": hl.tarray(hl.tfloat)})
# pop_ht.aggregate(hl.agg.counter(pop_ht.ancestry_pred))
N_SAMPLES = {'all':245394, 'afr':56913, 'amr':45035, 'eas':5706, 'eur':133581, 'mid':942, 'sas':3217}

def get_aou_sites_for_grm_path(extension: str, pruned: bool=False):
    return f'{root}/utils/aou_sites{"_pruned" if pruned else ""}_for_grm.{extension}'


mt_path = 'gs://prod-drc-broad/aou-wgs-delta-small_callsets_gq0/v7.1/exome_v7.1/splitMT/delta_basis_without_ext_aian_prod_gq0_3regions.exome.split.mt'

####### STEP 1: Downsample a bunch of sites based on frequency for the GRM  #######
# Step 1-1: build a docker image with `gnomad` module

def mac_category_case_builder(call_stats_expr):
    return (hl.case()
            .when(call_stats_expr.AC[0] <= 5, call_stats_expr.AC[0])
            .when(call_stats_expr.AC[0] <= 10, 10)
            .when(call_stats_expr.AC[0] <= 20, 20)
            .when(call_stats_expr.AF[0] <= 0.001, 0.001)
            .when(call_stats_expr.AF[0] <= 0.01, 0.01)
            .when(call_stats_expr.AF[0] <= 0.1, 0.1)
            .default(0.99))


def filter_ht_for_plink(
        ht: hl.Table,
        pop:str,
        min_call_rate: float = 0.95,
        variants_per_mac_category: int = 2000,
        variants_per_maf_category: int = 10000
):

    ht = ht.filter(ht.locus.in_autosome() &
                   (ht.info.AN >= N_SAMPLES[pop] * 2 * min_call_rate) &
                   (ht.info.AC[0] > 0))
    ht = ht.annotate(mac_category=mac_category_case_builder(ht.info))

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


###################### VAT per population random sampling ##################################
# def mac_category_case_builder(call_stats_ac_expr, call_stats_af_expr):
#     return (hl.case()
#             .when(call_stats_ac_expr <= 5, call_stats_ac_expr)
#             .when(call_stats_ac_expr <= 10, 10)
#             .when(call_stats_ac_expr <= 20, 20)
#             .when(call_stats_af_expr <= 0.001, 0.001)
#             .when(call_stats_af_expr<= 0.01, 0.01)
#             .when(call_stats_af_expr<= 0.1, 0.1)
#             .default(0.99))

# def filter_ht_for_plink(
#         ht: hl.Table,
#         pop:str,
#         min_call_rate: float = 0.95,
#         variants_per_mac_category: int = 2000,
#         variants_per_maf_category: int = 10000
# ):
#     from gnomad.utils.filtering import filter_to_autosomes
#     ht = filter_to_autosomes(ht)
#     ht = ht.filter((ht[f'gvs_{pop}_an'] >= N_SAMPLES[pop] * 2 * min_call_rate) &
#                    (ht[f'gvs_{pop}_ac']> 0))
#     ht = ht.annotate(mac_category=mac_category_case_builder(ht[f'gvs_{pop}_ac'], ht[f'gvs_{pop}_af']))
#
#     # From: https://hail.zulipchat.com/#narrow/stream/123010-Hail-Query-0.2E2-support/topic/.E2.9C.94.20randomly.20sample.20table/near/388162012
#     bins = ht.aggregate(hl.agg.collect_as_set(ht.mac_category))
#     ac_bins = [bin for bin in bins if bin >= 1]
#     af_bins = [bin for bin in bins if bin < 1]
#     binned_samples = ht.aggregate(
#         hl.if_else(ht.mac_category >= 1,
#                    hl.agg.array_agg(lambda x: hl.agg.filter(ht.mac_category == x,
#                                                             hl.agg._reservoir_sample(ht.key,
#                                                                                      variants_per_mac_category)), ac_bins),
#                    hl.agg.array_agg(lambda x: hl.agg.filter(ht.mac_category == x,
#                                                             hl.agg._reservoir_sample(ht.key,
#                                                                                      variants_per_maf_category)),
#                                     af_bins)
#                    )
#     )
#     samples = [sample for bin in binned_samples for sample in bin]
#     ht = hl.Table.parallelize(samples).key_by(*ht.key.keys())
#
#     category_counter = ht.aggregate(hl.agg.counter(ht.mac_category))
#     print(category_counter)
#     ht = ht.annotate_globals(category_counter=category_counter)
#     return ht

def main(args):
    hl.init(driver_memory='highmem', driver_cores=8, default_reference='GRCh38', log='/pre_process.log')

    if args.create_plink_file:
        mt = hl.read_matrix_table(mt_path)  # (34807589, 245394)
        ht = mt.rows()
        filtered_ht = filter_ht_for_plink(ht, pop='all')
        filtered_ht = filtered_ht.naive_coalesce(1000).checkpoint(get_aou_sites_for_grm_path(extension='ht'),
                                                                  _read_if_exists=not args.overwrite,
                                                                  overwrite=args.overwrite)

        mt = mt.filter_rows(hl.is_defined(filtered_ht[mt.row_key]))
        mt = mt.naive_coalesce(1000).checkpoint(get_aou_sites_for_grm_path(extension='mt'),
                                                _read_if_exists=not args.overwrite,
                                                overwrite=args.overwrite)

        mt = mt.unfilter_entries()
        ht = hl.ld_prune(mt.GT, r2=0.1)

        ht = ht.checkpoint(get_aou_sites_for_grm_path(extension='ht', pruned=True),
                           _read_if_exists=not args.overwrite,
                           overwrite=args.overwrite)
        mt = mt.filter_rows(hl.is_defined(ht[mt.row_key]))

        if args.overwrite or not hl.hadoop_exists(f'{get_aou_sites_for_grm_path(extension="bed", pruned=True)}'):
            hl.export_plink(mt, get_aou_sites_for_grm_path(extension="plink", pruned=True))



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--create-plink-file",
        help="Create plink files for GRM",
        action="store_true"
    )
    parser.add_argument('--overwrite', help='Whether to overwrite the existing file', action='store_true')
    args = parser.parse_args()

    main(args)

# python3 check_path.py
# hailctl config set batch/billing_project all-by-aou
# hailctl config set batch/remote_tmpdir gs://aou_wlu/tmp
# hailctl config set batch/tmp_dir gs://aou_wlu/tmp
# hailctl config list