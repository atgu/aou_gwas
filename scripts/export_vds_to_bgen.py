#!/usr/bin/env python3

__author__ = "wlu"

import hailtop.batch as hb
import argparse
import hail as hl
import logging
import hailtop
import pandas as pd
from collections import Counter

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("AoU_VDS_to_Bgen")
logger.setLevel(logging.INFO)

root = 'gs://aou_wlu'
data_path = f'{root}/data'
test_path = f'{root}/test'
GTF_PATH = 'gs://hail-common/references/gencode/gencode.v29.annotation.gtf.bgz'
ANCESTRY_PATH = 'gs://fc-aou-datasets-controlled/v7/wgs/short_read/snpindel/aux/ancestry/ancestry_preds.tsv'
N_GENE_PER_GROUP = 10

def summarize_genes_per_chrom(gtf):
    # Count number of genes per chromosome
    n_gene = gtf.group_by('chrom').aggregate(n_gene_per_chrom=hl.agg.count())

    # Round the numbers of genes per chromosome to the next multiplies of N_GENE_PER_GROUP
    n_gene = n_gene.annotate(n_gene_per_chrom_round=
                             hl.if_else(n_gene.n_gene_per_chrom % N_GENE_PER_GROUP > 0,
                                        n_gene.n_gene_per_chrom + N_GENE_PER_GROUP - n_gene.n_gene_per_chrom % N_GENE_PER_GROUP,
                                        n_gene.n_gene_per_chrom
                                        ))

    # Annotate numeric chromosome index and order by chromosome
    n_gene = n_gene.annotate(chrom_index=n_gene.chrom[3:])
    n_gene = n_gene.annotate(chrom_index=hl.case()
                             .when(n_gene.chrom_index == 'X', '23')
                             .when(n_gene.chrom_index == 'Y', '24')
                             .when(n_gene.chrom_index == 'M', '25')
                             .default(n_gene.chrom_index)
                             )
    n_gene = n_gene.annotate(chrom_index=hl.int64(n_gene.chrom_index))
    n_gene = n_gene.order_by('chrom_index')

    # Compute cumulative number of genes before each chromosome
    n_gene_lst = n_gene.aggregate(hl.cumulative_sum(hl.agg.collect(n_gene.n_gene_per_chrom)))
    n_gene_lst.insert(0, 0)
    n_gene_lst.pop()

    # Compute cumulative rounded number of genes before each chromosome
    n_gene_round_lst = n_gene.aggregate(hl.cumulative_sum(hl.agg.collect(n_gene.n_gene_per_chrom_round)))
    n_gene_round_lst.insert(0, 0)
    n_gene_round_lst.pop()

    # Build Hail table with the information above
    chrom_lst = [f'chr{i + 1}' for i in range(22)] + ['chrX', 'chrY', 'chrM']
    n_gene_df = pd.DataFrame(data={'chrom': chrom_lst, 'cum_n_gene': n_gene_lst, 'cum_n_gene_round': n_gene_round_lst})
    n_gene_ht = hl.Table.from_pandas(n_gene_df, key='chrom')

    return n_gene_ht

def create_gene_group_interval():
    gtf = hl.experimental.import_gtf(GTF_PATH,
                                     reference_genome='GRCh38',
                                     skip_invalid_contigs=True)
    # Filter to protein coding genes
    gtf = gtf.filter((gtf.feature == 'gene') & (gtf.gene_type == 'protein_coding'))

    # Annotate basic interval information
    gtf = gtf.annotate(
        chrom=gtf.interval.start.contig,
        start_pos=gtf.interval.start.position,
        end_pos=gtf.interval.end.position
    )
    gtf = gtf.annotate(length=gtf.end_pos - gtf.start_pos + 1)
    gtf = gtf.add_index()

    # Compute and annotate cumulative number of genes per chromosome
    n_gene_ht = summarize_genes_per_chrom(gtf)
    gtf = gtf.annotate(n_previous_genes=n_gene_ht[gtf.chrom].cum_n_gene,
                       n_previous_genes_round=n_gene_ht[gtf.chrom].cum_n_gene_round)

    # Edit global index
    gtf = gtf.annotate(new_idx=gtf.idx + gtf.n_previous_genes_round - gtf.n_previous_genes)
    # Annotate group index
    gtf = gtf.annotate(group_id=gtf.new_idx // N_GENE_PER_GROUP)

    # If the number of genes in a group is < 5, merge this group with the previous group
    group_cnt = gtf.group_by('group_id').aggregate(cnt=hl.agg.count())
    gtf = gtf.annotate(n_genes_per_group=group_cnt[gtf.group_id].cnt)
    gtf = gtf.annotate(group_id=hl.if_else(gtf.n_genes_per_group < 5, gtf.group_id - 1, gtf.group_id))

    # Summarize and print number of genes per group
    group_dict = gtf.aggregate(hl.agg.counter(gtf.group_id))
    print(f'Number of genes per group: {Counter(group_dict.values())}')

    sub_gtf = gtf.select('gene_name', 'gene_id', 'transcript_id',
                         'chrom', 'start_pos', 'end_pos', 'length', 'idx', 'new_idx', 'group_id', 'n_previous_genes')

    return sub_gtf


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

def export_vds_to_bgen(vds_path, pop_ht, interval, output_dir, mean_impute_missing, no_adj, callrate_filter):
    def read_vds(vds_path):
        vds = hl.vds.read_vds(vds_path)
        return vds
    vds = read_vds(vds_path)
    outname = f"gene_{interval.start.contig}_{interval.start.position}_{interval.end.position}"
    print(outname)

    # Filter to interval
    sub_vds = hl.vds.filter_intervals(vds, hl.literal([interval]))
    sub_vds = hl.vds.filter_samples(sub_vds, pop_ht, keep=True, remove_dead_alleles=True)

    # Densify subset VDS
    mt = hl.vds.to_dense_mt(sub_vds)
    print(mt.describe())

    if not no_adj:
        mt = mt.filter_entries(mt.adj)

    mt = mt.select_entries('LGT') # Select LGT field
    mt = mt.filter_rows(hl.agg.count_where(mt.LGT.is_non_ref()) > 0) # Filter to non-reference sites
    mt = mt.annotate_rows(rsid=mt.locus.contig + ':' + hl.str(mt.locus.position) + '_' + mt.alleles[0] + '/' + mt.alleles[1]) # Annotate rsid

    if callrate_filter:
        mt = mt.filter_rows(hl.agg.fraction(hl.is_defined(mt.LGT)) >= args.callrate_filter)

    mt = mt.annotate_entries(GT=hl.if_else(mt.LGT.is_haploid(), hl.call(mt.LGT[0], mt.LGT[0]), mt.LGT))
    mt = gt_to_gp(mt)
    mt = impute_missing_gp(mt, mean_impute=mean_impute_missing)
    hl.export_bgen(mt, f'{output_dir}/{outname}', gp=mt.GP, varid=mt.rsid)

    try:
        hl.hadoop_ls('gs://fc-aou-datasets-controlled')
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

        if args.update_gene_intervals:
            group_ht = create_gene_group_interval()
            group_ht.write(f"gs://aou_wlu/data/group_positions_{N_GENE_PER_GROUP}_protein_coding.ht")

        int_ht = hl.read_table(args.input_interval_path)
        intervals = int_ht.aggregate(hl.agg.collect(int_ht.interval))

        pop_ht = hl.read_table(ANCESTRY_PATH)
        pop_ht = pop_ht.filter(pop_ht.ancestry_pred == args.pop)

        output_dir = test_path if args.test else f"{root}/bgen/{args.data_type}/{args.pop}"

        for interval in intervals:
            j = b.new_python_job(name=f"export_{interval}")
            j.call(export_vds_to_bgen, args.input_vds_path, pop_ht, interval, output_dir, args.mean_impute_missing, args.no_adj, args.callrate_filter)
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
        default="gs://aou_wlu/data/group_positions_10_protein_coding.ht"
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
    parser.add_argument(
        "--update-gene-intervals",
        help="Update grouped gene intervals",
        action="store_true"
    )
    parser.add_argument('--pop', help='Comma-separated list of pops to run', choices=['afr', 'amr', 'eas', 'eur', 'mid', 'sas'])
    parser.add_argument('--no_adj', help='Use all genotypes instead of only high-quality ones', action='store_true')
    parser.add_argument('--callrate_filter', help='Impose filter of specified callrate (default: none)', default=0.0, type=float)
    args = parser.parse_args()

    main(args)

# python3 export_vds_to_bgen.py --test --mean_impute_missing --update-vds
# hailctl config set batch/billing_project all-by-aou
# hailctl config set batch/remote_tmpdir gs://aou_wlu/tmp