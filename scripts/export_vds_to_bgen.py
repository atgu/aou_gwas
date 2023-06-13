#!/usr/bin/env python3

__author__ = "wlu"

import hailtop.batch as hb
import argparse
import hail as hl
import logging
import hailtop

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("AoU_VDS_to_Bgen")
logger.setLevel(logging.INFO)

root = 'gs://aou_wlu'
data_path = f'{root}/data'
test_path = f'{root}/test'

def gt_to_gp(mt, location: str = 'GP'):
    return mt.annotate_entries(**{location: hl.or_missing(
        hl.is_defined(mt.LGT),
        hl.map(lambda i: hl.if_else(mt.LGT.unphased_diploid_gt_index() == i, 1.0, 0.0),
               hl.range(0, hl.triangle(hl.len(mt.alleles)))))})


def impute_missing_gp(mt, location: str = 'GP', mean_impute: bool = True):
    mt = mt.annotate_entries(_gp = mt[location])
    if mean_impute:
        mt = mt.annotate_rows(_mean_gp=hl.agg.array_agg(lambda x: hl.agg.mean(x), mt._gp))
        gp_expr = mt._mean_gp
    else:
        gp_expr = [1.0, 0.0, 0.0]
    return mt.annotate_entries(**{location: hl.or_else(mt._gp, gp_expr)}).drop('_gp')

def export_vds_to_bgen(vds_path, interval, output_dir, mean_impute_missing):
    def read_vds(vds_path):
        vds = hl.vds.read_vds(vds_path)
        return vds
    vds = read_vds(vds_path)
    outname = f"gene_{interval.start.contig}_{interval.start.position}_{interval.end.position}"
    print(outname)
    sub_vds = hl.vds.filter_intervals(vds, hl.literal([interval]))
    mt = hl.vds.to_dense_mt(sub_vds)
    mt = mt.annotate_rows(
        rsid=mt.locus.contig + ':' + hl.str(mt.locus.position) + '_' + mt.alleles[0] + '/' + mt.alleles[1])
    mt = mt.annotate_entries(GT=hl.if_else(mt.LGT.is_haploid(), hl.call(mt.LGT[0], mt.LGT[0]), mt.LGT))
    mt = gt_to_gp(mt)
    mt = impute_missing_gp(mt, mean_impute=mean_impute_missing)
    hl.export_bgen(mt, f'{output_dir}/{outname}', gp=mt.GP, varid=mt.rsid)
    try:
        print(hl.vds.read_vds(
            'gs://fc-aou-datasets-controlled/v7/wgs/short_read/snpindel/vds/hail.vds').variant_data.count())
    except:
        print('NO ACCESS YET')



def main(args):
    try:
        hl.init(default_reference="GRCh38")

        # image = hailtop.batch.docker.build_python_image('us-central1-docker.pkg.dev/aou-neale-gwas/hail-batch/hail_python_3.9',
        #                                                 requirements=['hail'], python_version='3.9')

        backend = hb.ServiceBackend(
            billing_project="all-by-aou",
            remote_tmpdir=f"{root}/tmp",
        )
        b = hb.Batch(
            name=f"AoU_VDS_to_Bgen",
            # default_python_image='hailgenetics/hail:0.2.117',
            requester_pays_project='aou-neale-gwas',
            default_python_image='us-central1-docker.pkg.dev/aou-neale-gwas/hail-batch/wenhan-hail:0.2.117',
            backend=backend,
        )

        if args.update_vds:
            vds = hl.vds.read_vds("gs://aou_wlu/hgdp_decontaminated_3_samples.vds")
            vds2 = hl.vds.truncate_reference_blocks(vds, max_ref_block_base_pairs=5000)
            vds2.write(args.input_vds_path)
        int_ht = hl.read_table(args.input_interval_path)
        intervals = int_ht.aggregate(hl.agg.collect(int_ht.interval))
        output_dir = test_path if args.test else f"{root}/bgen/{args.data_type}"
        for interval in intervals:
            j = b.new_python_job(name=f"export_{interval}")
            j.call(export_vds_to_bgen, args.input_vds_path, interval, output_dir, args.mean_impute_missing)
            if args.test:
                break
        b.run()


    finally:
        from datetime import date

        logger.info("Copying log to logging bucket...")
        hl.copy_log(f"{root}/log/bgen_{date.today()}.log")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--data-type",
        nargs="?",
        help="Dataset to export bgens from",
        default="gene",
        choices=["gene", "variant"],
    )
    parser.add_argument(
        "--input-vds-path",
        help="Path to a VDS file",
        nargs="?",
        default="gs://aou_wlu/hgdp_decontaminated_3_samples_ref_truncated.vds"
    )
    parser.add_argument(
        "--input-interval-path",
        help="Path to a hail table with interval info stored",
        nargs="?",
        default="gs://aou_wlu/data/group_positions_10.ht"
    )
    parser.add_argument('--mean_impute_missing', help='Whether to mean impute missing genotypes (BGEN only) '
                                                      '(default: set to hom ref)', action='store_true')
    parser.add_argument(
        "--test",
        help="Whether to test with one bgen or not",
        action="store_true"
    )
    parser.add_argument(
        "--update-vds",
        help="Update VDS with truncated ref blocks",
        action="store_true"
    )
    args = parser.parse_args()

    main(args)

# python3 export_vds_to_bgen.py --test --mean_impute_missing --update-vds
# hailctl config set batch/billing_project all-by-aou
# hailctl config set batch/remote_tmpdir gs://aou_wlu/tmp