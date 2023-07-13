#!/usr/bin/env python3

__author__ = 'konradk, wlu'

import argparse
import logging
from gnomad.utils.vep import *
from gnomad.utils.filtering import *
from ukbb_common import *
from ukbb_pan_ancestry import *

logger = logging.getLogger("ALL_x_AoU_preprocessing_SAIGE_data")

BUCKET = 'gs://aou-analysis/'
TRANCHE = '?'
GWAS_AF_CUTOFF = 0.001
RVAS_AF_CUTOFF = 0.01
CALLRATE_CUTOFF = 0.8
N_GENE_PER_GROUP = 10
CURRENT_BUCKET = f'{BUCKET}/{TRANCHE}'
POPS = ('afr', 'amr', 'eas', 'eur', 'mid', 'sas')

# ALL PATHs: https://support.researchallofus.org/hc/en-us/articles/4616869437204-Controlled-CDR-Directory
# Descriptions: https://support.researchallofus.org/hc/en-us/articles/4614687617556-How-the-All-of-Us-Genomic-data-are-organized
WGS_VDS_PATH = 'gs://fc-aou-datasets-controlled/v7/wgs/short_read/snpindel/vds/hail.vds'
WES_MT_PATH = "gs://fc-aou-datasets-controlled/v7/wgs/short_read/snpindel/exome/splitMT/hail.mt"
VAT_TSV_PATH = 'gs://fc-aou-datasets-controlled/v7/wgs/short_read/snpindel/aux/vat/vat_complete.bgz.tsv.gz'
ANCESTRY_PATH = 'gs://fc-aou-datasets-controlled/v7/wgs/short_read/snpindel/aux/ancestry/ancestry_preds.tsv'
SAMPLE_QC_PATH = 'gs://fc-aou-datasets-controlled/v7/wgs/short_read/snpindel/aux/qc/flagged_samples.tsv'
RELATEDNESS_PATH = 'gs://fc-aou-datasets-controlled/v7/wgs/short_read/snpindel/aux/relatedness/relatedness.tsv'
GENE_INTERVAL_PATH = f"gs://aou_wlu/data/group_positions_{N_GENE_PER_GROUP}_protein_coding.ht"

def do_annotate(mt, meta_ht, by_col=True):
    if by_col:
        mt = mt.annotate_cols(**meta_ht[mt.col_key])
    else:
        mt = mt.annotate_rows(**meta_ht[mt.row_key])
    return mt

def get_aou_vep_path(tranche: str = TRANCHE):
    return f'{CURRENT_BUCKET}/util/aou_vep_{tranche}.ht'

def get_aou_vat_path(tranche: str = TRANCHE):
    return f'{CURRENT_BUCKET}/util/aou_vat_{tranche}.ht'

def get_aou_grm_mt_path(pop: str, data_iteration:int, tranche: str = TRANCHE):
    suffix = f'_{data_iteration}' if data_iteration else ""
    return f'{CURRENT_BUCKET}/grm/aou_{pop}_for_grm{suffix}_{tranche}.mt'

def get_aou_grm_pruned_ht_path(pop: str, window_size: str = '1e6', tranche: str = TRANCHE):
    cut = '' if window_size == '1e6' else f'_{window_size}'
    return f'{CURRENT_BUCKET}/grm/aou_{pop}_for_grm_pruned{cut}_{tranche}.ht'

def get_aou_grm_plink_path(pop: str, data_iteration: int = 0, window_size: str = '1e6', tranche: str = TRANCHE):
    suffix = f'.{data_iteration}' if data_iteration else ""
    cut = '' if window_size == '1e6' else f'.{window_size}'
    return f'{CURRENT_BUCKET}/grm/aou_{pop}_for_grm{suffix}_pruned{cut}_{tranche}.plink'

def get_vat_field_types():
    tags = ('afr', 'amr', 'asj', 'eas', 'eur', 'fin', 'mid', 'nfr', 'sas', 'oth', 'max', 'all')
    fields = ('ac', 'an', 'sc')
    types = {}
    for tag in tags:
        types[f'gvs_{tag}_af'] = hl.tfloat64
        types[f'gnomad_{tag}_af'] = hl.tfloat64
        for field in fields:
            types[f'gvs_{tag}_{field}'] = hl.tint32
            types[f'gnomad_{tag}_{field}'] = hl.tint32
    types['is_canonical_transcript'] = hl.tbool
    types['omim_phenotypes_id'] = hl.tarray(hl.tint32)
    for x in ['position', 'gene_omim_id', 'splice_ai_acceptor_gain_distance', 'splice_ai_acceptor_loss_distance',
              'splice_ai_donor_gain_distance', 'splice_ai_donor_loss_distance']:
        types[x] = hl.tint32
    for x in ['revel', 'splice_ai_acceptor_gain_score', 'splice_ai_acceptor_loss_score', 'splice_ai_donor_gain_score',
              'splice_ai_donor_loss_score']:
        types[x] = hl.tfloat64
    for x in ['omim_phenotypes_name', 'clinvar_classification', 'clinvar_phenotype', 'consequence', 'dbsnp_rsid']:
        types[x] = hl.tarray(hl.tstr)
    return types


def main(args):
    hl.init(default_reference='GRCh38')

    # Load and VDS and split-multi
    vds = hl.vds.read_vds(WGS_VDS_PATH)
    # Pre splitting multi: vmt: (1031611675, 245394), rmt: (2853640141, 245394)
    vds = hl.vds.split_multi(vds, filter_changed_loci=False)

    # Load variant info table and re-format
    vat_ht = hl.import_table(VAT_TSV_PATH, force=True, quote='"', delimiter="\t", force_bgz=True, types = get_vat_field_types())
    vat_ht = vat_ht.annotate(locus=hl.locus(vat_ht.contig, vat_ht.position, reference_genome='GRCh38'),
                             alleles=[vat_ht.ref_allele, vat_ht.alt_allele])
    vat_ht = vat_ht.key_by('locus', 'alleles')
    vat_ht = vat_ht.checkpoint(get_aou_vat_path(), _read_if_exists=not args.overwrite, overwrite=args.overwrite)
    vat_ht = vat_ht.filter(vat_ht.gvs_all_ac > 0)
    pop_ht = hl.import_table(ANCESTRY_PATH, key="research_id", impute=True,
                             types={"research_id": "tstr", "pca_features": hl.tarray(hl.tfloat)})
    # {'afr': 56913, 'amr': 45035, 'eas': 5706, 'eur': 133581, 'mid': 942, 'sas': 3217} - Total: 245394

    if args.test:
        # filter to the first 1% partition (Total: 84648 partitions)
        n_partitions = round(vds.variant_data.n_partitions()*0.01)
        vds = hl.vds.VariantDataset(vds.reference_data._filter_partitions(range(n_partitions)),
                                    vds.variant_data._filter_partitions(range(n_partitions)))


    pops = args.pops.split(',') if args.pops else POPS

    if args.vep:
        # Vep resource from gnomAD:
        # https://github.com/broadinstitute/gnomad_methods/blob/e2f2b055a2d40991851e7852dea435b4eaef45ea/gnomad/resources/grch38/reference_data.py#L100
        vep_ht = vep_or_lookup_vep(vat_ht, vep_version='105')
        vep_ht.write(get_aou_vep_path(), args.overwrite)

    if args.create_plink_file:
        window = args.window
        # Note: 1e7 LD pruning for EUR was run with r2=0.05, and chr8 inversion and HLA were removed
        r2 = 0.1 if args.window != '1e7' else 0.05
        iteration = 1
        for pop in pops:
            # Filter samples
            sub_pop_ht = pop_ht.filter(pop_ht.ancestry_pred == pop)
            sub_vds = hl.vds.filter_samples(vds, sub_pop_ht, keep=True, remove_dead_alleles=True)
            n_samples = sub_vds.variant_data.count_cols()
            print(f'Got {n_samples} samples for {pop}...')

            # Filter variants
            call_stats_ht = vat_ht.select(
                AF=hl.float64(vat_ht[f'gvs_{pop}_af']),
                AC=vat_ht[f'gvs_{pop}_ac'],
                AN=vat_ht[f'gvs_{pop}_an'],
                N=vat_ht[f'gvs_{pop}_sc'],
            ) # with_x = False?
            if args.data_type == 'gwas':
                call_stats_ht = call_stats_ht.filter(call_stats_ht.AF > GWAS_AF_CUTOFF)
            elif args.data_type == 'rvas':
                int_ht = hl.read_table(GENE_INTERVAL_PATH)
                sub_vds = hl.vds.filter_intervals(sub_vds, int_ht, keep=True)
                call_stats_ht = call_stats_ht.filter(call_stats_ht.AF < RVAS_AF_CUTOFF)
            else:
                logger.warning('INVALID TEST TYPE...')

            call_stats_ht = call_stats_ht.filter(call_stats_ht.AN/(2*n_samples) > CALLRATE_CUTOFF)
            sub_vds = hl.vds.filter_chromosomes(sub_vds, keep_autosomes=True)
            sub_vds = hl.vds.filter_variants(sub_vds, call_stats_ht, keep=True)

            if pop == 'eur':
                print('Removing HLA...')
                # Common inversion taken from Table S4 of https://www.ncbi.nlm.nih.gov/pubmed/27472961
                # (converted to GRCh38 by: https://liftover.broadinstitute.org/#input=chr8%3A8055789-11980649&hg=hg19-to-hg38 )
                # Also removing HLA, from https://www.ncbi.nlm.nih.gov/grc/human/regions/MHC?asm=GRCh38
                sub_vds = hl.vds.filter_intervals(sub_vds,
                                                  hl.literal([hl.parse_locus_interval('chr8:8198267-12123140', reference_genome='GRCh38'),
                                                              hl.parse_locus_interval('chr6:28510120-33480577', reference_genome='GRCh38')]),
                                                  keep=False)
            mt = hl.vds.to_dense_mt(sub_vds)
            mt = mt.checkpoint(get_aou_grm_mt_path(pop, iteration), _read_if_exists=not args.overwrite, overwrite=args.overwrite)
            mt = mt.unfilter_entries()

            if not args.omit_ld_prune:
                print(f'LD-pruning to {get_aou_grm_pruned_ht_path(pop, window)}...')
                ht = hl.ld_prune(mt.LGT, r2=float(r2), bp_window_size=int(float(window)), block_size=1024)
                ht.write(get_aou_grm_pruned_ht_path(pop, window), overwrite=args.overwrite)
            ht = hl.read_table(get_aou_grm_pruned_ht_path(pop, window))
            mt = mt.filter_rows(hl.is_defined(ht[mt.row_key]))
            if pop == 'eur':
                mt = mt.filter_rows(hl.rand_bool(0.55))

            if args.overwrite or not hl.hadoop_exists(f'{get_aou_grm_plink_path(pop, iteration, window)}.bed'):
                print(f'Exporting plink to {get_aou_grm_plink_path(pop, iteration, window)}...')
                hl.export_plink(mt, get_aou_grm_plink_path(pop, iteration, window))

    if args.genotype_summary:
        variants = hl.read_table(ukb_imputed_info_ht_path)
        print(variants.count())
        variants = variants.filter(variants.info > 0.8)
        print(variants.count())
        meta_ht = hl.import_table(get_ukb_meta_pop_tsv_path(), impute=True, types={'s': hl.tstr}, key='s')

        mt = get_ukb_imputed_data('all', variant_list=variants, entry_fields=('dosage', ))
        mt = mt.annotate_cols(**meta_ht[mt.col_key])
        ht = mt.annotate_rows(af=hl.agg.group_by(mt.pop, hl.agg.mean(mt.dosage)),
                              an=hl.agg.group_by(mt.pop, hl.agg.count_where(hl.is_defined(mt.dosage)))).rows()
        ht = ht.checkpoint(get_ukb_af_ht_path(False), args.overwrite, _read_if_exists=not args.overwrite)

        mt = get_ukb_imputed_data('X', variant_list=variants, entry_fields=('dosage', ))
        mt = mt.annotate_cols(**meta_ht[mt.col_key])
        ht_x = mt.annotate_rows(af=hl.agg.group_by(mt.pop, hl.agg.mean(mt.dosage)),
                              an=hl.agg.group_by(mt.pop, hl.agg.count_where(hl.is_defined(mt.dosage)))).rows()
        ht = ht.union(ht_x)
        ht = ht.checkpoint(get_ukb_af_ht_path(), args.overwrite, _read_if_exists=not args.overwrite)
        ht = ht.naive_coalesce(1000).checkpoint(get_ukb_af_ht_path(repart=True), args.overwrite, _read_if_exists=not args.overwrite)

        print(ht.aggregate(hl.struct(
            # hist=hl.agg.hist(hl.sum(ht.an.values()), 0, total_samples, 10),  # No missing data
            # fraction_missingness=hl.agg.fraction(hl.sum(ht.an.values()) < total_samples),
            # number_sites_above_001=hl.agg.count_where(hl.any(lambda x: x / 2 > 0.001, ht.af.values())),
            # number_sites_above_005=hl.agg.count_where(hl.any(lambda x: x / 2 > 0.005, ht.af.values())),
            # number_sites_above_01=hl.agg.count_where(hl.any(lambda x: x / 2 > 0.01, ht.af.values())),
            # number_sites_above_05=hl.agg.count_where(hl.any(lambda x: x / 2 > 0.05, ht.af.values())),
            # number_sites_above_10=hl.agg.count_where(hl.any(lambda x: x / 2 > 0.1, ht.af.values()))
            number_sites_above_mac_20=hl.agg.count_where(hl.any(lambda x: ht.af[x] * ht.an[x] >= 20, hl.literal(POPS))),
            number_run_sites_above_mac_20=hl.agg.sum(hl.sum(hl.map(lambda x: hl.int(ht.af[x] * ht.an[x] >= 20), hl.literal(POPS))))
        )))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-type", nargs="?", help="Type of analysis to conduct", default="gwas", choices=["gwas", "rvas"])
    parser.add_argument("--test", help="Whether to run test on a subset of the data", action="store_true")

    parser.add_argument('--pheno_summary', action='store_true')
    parser.add_argument('--prepare_genotype_data', action='store_true')
    parser.add_argument('--genotype_summary', action='store_true')
    parser.add_argument('--overwrite', help='Overwrite everything', action='store_true')
    parser.add_argument('--omit_ld_prune', help='Overwrite everything', action='store_true')
    parser.add_argument('--window', help='Overwrite everything', default='1e6')
    parser.add_argument('--pops', help='Comma-separated list of pops to run', default=['afr', 'amr', 'eas', 'eur', 'mid', 'sas'])
    parser.add_argument('--create_plink_file', help='Overwrite everything', action='store_true')
    parser.add_argument('--vep', help='Overwrite everything', action='store_true')
    parser.add_argument('--slack_channel', help='Send message to Slack channel/user', default='@konradjk')
    args = parser.parse_args()

    # if args.slack_channel:
    #     try_slack(args.slack_channel, main, args)
    # else:
    #     main(args)
    main(args)

# hailctl dataproc submit wlu2 ~/PycharmProjects/ukbb_pan_ancestry/pre_process_saige_data.py --create_plink_file --window 1e7 --pops all