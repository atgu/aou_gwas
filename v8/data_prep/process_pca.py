from typing import Union
import hailtop.fs as hfs
import argparse
import hail as hl
import sys
import warnings

OVERWRITE = False
OVERWRITE_2 = True
CALLRATE_CUTOFF = 0.9
MIN_AF = 0.01
TMP_BUCKET = 'gs://aou_tmp/v8'
DATA_PATH = 'gs://aou_analysis/v8/data'
EXOME_MT_PATH = f'gs://fc-aou-datasets-controlled/v8/wgs/short_read/snpindel/exome/splitMT/hail.mt'
ACAF_MT_PATH =  f'gs://fc-aou-datasets-controlled/v8/wgs/short_read/snpindel/acaf_threshold/splitMT/hail.mt'
# Discard: N_SAMPLES_HARD_FILTERED = {'afr': 82578, 'amr': 78354, 'eas': 10032, 'eur': 232359, 'mid': 1521, 'sas': 5544} # computed from defined_ht below, didn't remove duplicated samples
# Discard: N_SAMPLES_HARD_FILTERED = {'afr': 81924, 'amr': 78031, 'eas': 10010, 'eur': 231440, 'mid': 1517, 'sas': 5535}  # computed from defined_ht below, removed duplicated samples but not qced samples
N_SAMPLES_HARD_FILTERED = {'afr': 81595, 'amr': 77772, 'eas': 9941, 'eur': 231199, 'mid': 1501, 'sas': 5482, 'all': 407490}
N_SAMPLES_PRUNED = {'afr': 77444, 'amr': 71540, 'eas': 9488, 'eur': 0, 'mid': 1153, 'sas': 5132, 'all': 0}

def get_adj_expr(
   gt_expr: hl.expr.CallExpression,
   gq_expr: Union[hl.expr.Int32Expression, hl.expr.Int64Expression],
   ad_expr: hl.expr.ArrayNumericExpression,
   adj_gq: int = 30,
   adj_ab: float = 0.2,
) -> hl.expr.BooleanExpression:
   """
   Get adj genotype annotation.


   Defaults correspond to gnomAD values.
   """
   return (
       (gq_expr >= adj_gq)
       & (
           hl.case()
           .when(~gt_expr.is_het(), True)
           .when(gt_expr.is_het_ref(), ad_expr[gt_expr[1]] / hl.sum(ad_expr) >= adj_ab)
           .default(
               (ad_expr[gt_expr[0]] / hl.sum(ad_expr) >= adj_ab)
               & (ad_expr[gt_expr[1]] / hl.sum(ad_expr) >= adj_ab)
           )
       )
   )




def annotate_adj(
       mt: hl.MatrixTable,
       adj_gq: int = 30,
       adj_ab: float = 0.2,
) -> hl.MatrixTable:
   """
   Annotate genotypes with adj criteria (assumes diploid).


   Defaults correspond to gnomAD values.
   """
   if "GT" not in mt.entry and "LGT" in mt.entry:
       print("No GT field found, using LGT instead.")
       gt_expr = mt.LGT
   else:
       gt_expr = mt.GT


   if "AD" not in mt.entry and "LAD" in mt.entry:
       print("No AD field found, using LAD instead.")
       ad_expr = mt.LAD
   else:
       ad_expr = mt.AD


   return mt.annotate_entries(
       adj=get_adj_expr(
           gt_expr, mt.GQ, ad_expr, adj_gq, adj_ab
       )
   )


def get_filtered_mt(mt_type: str,
                    sample_ids: hl.Table,
                    ancestry: str='all',
                    filter_samples: bool=True, 
                    filter_variants: bool=True,
                    prune_samples:bool=True,
                    adj_filter: bool=True):
    anc_mt_path = f'{DATA_PATH}/utils/raw_mt/{mt_type}/{ancestry.upper()}_{mt_type}.mt'
    if ancestry != 'all':
        ancestry_ht = hl.read_table(f'{DATA_PATH}/utils/aou_v8_global_pca.ht')
        ancestry_ht = ancestry_ht.filter(ancestry_ht.ancestry_pred == ancestry)
        if not hfs.exists(f'{anc_mt_path}/_SUCCESS'):
            mt_path = ACAF_MT_PATH if mt_type=='ACAF' else EXOME_MT_PATH
            print(mt_path)
            mt = hl.read_matrix_table(mt_path)
            mt.describe()
            print(mt.count())
            mt = mt.filter_entries(hl.is_missing(mt.FT) | (mt.FT == 'PASS'))
            mt = mt.filter_cols(hl.is_defined(ancestry_ht[mt.col_key]))
            mt = mt.naive_coalesce(20000).checkpoint(anc_mt_path, overwrite=True)
        mt = hl.read_matrix_table(anc_mt_path)
        print(mt._force_count_cols())
    else:
        mt_path = ACAF_MT_PATH if mt_type=='ACAF' else EXOME_MT_PATH
        if not hfs.exists(f'{anc_mt_path}/_SUCCESS'):
            mt = hl.read_matrix_table(mt_path)
            mt = mt.filter_entries(hl.is_missing(mt.FT) | (mt.FT == 'PASS'))
            mt = mt.naive_coalesce(20000).checkpoint(anc_mt_path, overwrite=True)
        mt = hl.read_matrix_table(anc_mt_path)
        print(mt.count())

    if filter_variants:
        print(f'Filtering to global AC > 0...')
        mt = mt.filter_rows(
            (mt.info.AC[0] > 0)
        )
        
    if filter_samples:
        print(f'Filtering to samples with sex, ancestry defined and age <= 100...')
        mt = mt.filter_cols(hl.is_defined(sample_ids[mt.col_key]))
        
    if prune_samples:
        print(f'[{mt_type}] Filtering to samples pruned from the PCA centroid pipeline...')
        pca_pruned_sample_tsv = f'{DATA_PATH}/utils/pca/results/aou_{ancestry.lower()}_centroid_pruned.tsv'
        pca_pruned_sample_path = f'{DATA_PATH}/utils/pca/results/{ancestry.lower()}_pca_centroid_pruned.ht'
        overwrite_pruned_ht = True
        if hfs.exists(pca_pruned_sample_tsv):
            if not hfs.exists(f'{pca_pruned_sample_path}/_SUCCESS') or overwrite_pruned_ht:
                pruned_ht = hl.import_table(pca_pruned_sample_tsv, delimiter='\t', key='s', types = {'s':hl.tstr})
                pruned_ht = pruned_ht.checkpoint(pca_pruned_sample_path, overwrite=overwrite_pruned_ht)
            else: 
                pruned_ht = hl.read_table(pca_pruned_sample_path)
        print(pruned_ht.count())
        pruned_ht = hl.read_table(pca_pruned_sample_path)
        mt = mt.filter_cols(hl.is_defined(pruned_ht[mt.col_key]))
        
    if adj_filter:
        """
       Filter genotypes to adj criteria - Default:
       GQ >= 30, het_ref_altAB >= 0.2, het_non_ref_altAB >= 0.2 , het_non_ref_refAB >= 0.2
       """
        print(f'Applying adj filters...')
        mt = annotate_adj(mt)
        mt = mt.filter_entries(mt.adj)
        mt = mt.filter_rows(hl.agg.any(mt.GT.n_alt_alleles() >0))
    
    return mt


def main(args):
    ANCESTRY = args.ancestry if args.ancestry is not None else 'all'
    hl.init(
        app_name=f'Process_pca_{ANCESTRY}',
        gcs_requester_pays_configuration='aou-neale-gwas',
        tmp_dir=f"{TMP_BUCKET}/pca",
        default_reference="GRCh38",
        log="/process_pca.log"
    )

    relatedness_flag_sample = hl.read_table(f'{DATA_PATH}/utils/aou_v8_related_samples.ht')
    print(f'Number of samples in relatedness flag: {relatedness_flag_sample.count()}')
    relatedness_flag_sample.describe()

    kinship_ht = hl.read_table(f'{DATA_PATH}/utils/aou_v8_related_kinship.ht')
    print(f'Number of samples in kinship: {kinship_ht.count()}')
    kinship_ht.describe()

    demographic_ht = hl.read_table(f'{DATA_PATH}/utils/aou_v8_demographic.ht')
    print(f'Number of samples in demographics: {demographic_ht.count()}')
    print(demographic_ht.aggregate(hl.agg.counter(demographic_ht.sex_at_birth)))
    demographic_ht.describe()

    ancestry_ht = hl.read_table(f'{DATA_PATH}/utils/aou_v8_global_pca.ht')
    print(f'Number of samples in global pca: {ancestry_ht.count()}')
    ancestry_ht.describe()

    qc_ht = hl.read_table(f'{DATA_PATH}/utils/aou_v8_flagged_samples.ht')
    print(f'Number of samples in qc: {qc_ht.count()}')
    qc_ht.describe()
    

    if not hfs.exists(f'{DATA_PATH}/utils/aou_v8_hard_filter.ht') or args.overwrite_hard_filters:
        duplicated_ht = kinship_ht.filter(kinship_ht.kin > 0.375)
        duplicated_ht = duplicated_ht.key_by('j.s').select().distinct()
        print(duplicated_ht.count())
        defined_ht = ancestry_ht.annotate(**relatedness_flag_sample[ancestry_ht.key],
                                          **demographic_ht[ancestry_ht.key])
        defined_ht = defined_ht.filter(hl.is_defined(defined_ht.ancestry_pred) & 
                                   (defined_ht.age_at_cdr <= 100) & 
                                   (hl.literal(['Female', 'Male']).contains(defined_ht.sex_at_birth))&
                                   (hl.is_missing(duplicated_ht[defined_ht.key])) &
                                   (hl.is_missing(qc_ht[defined_ht.key]))
                                   )
        defined_ht = defined_ht.checkpoint(f'{DATA_PATH}/utils/aou_v8_hard_filter.ht', overwrite = True)
        defined_ht.describe()
        print(defined_ht.count())
    defined_ht = hl.read_table(f'{DATA_PATH}/utils/aou_v8_hard_filter.ht')
    print(f'Number of samples post hard filter: {defined_ht.count()}')
    # defined_ht.describe()
    # N_SAMPLES_HARD_FILTERED = defined_ht.aggregate(hl.agg.counter(defined_ht.ancestry_pred))
    # print(N_SAMPLES_HARD_FILTERED)    

    callstats_ht_path = f'{DATA_PATH}/utils/call_stats/ACAF_pre_pruning/{ANCESTRY.upper()}_ACAF_call_stats.ht'
    if not hfs.exists(f'{callstats_ht_path}/_SUCCESS') or args.overwrite_call_stats:
        if not args.batch and ANCESTRY.lower() != 'eur':
            warnings.warn("Please run in QoB mode for small ancestries")
        print(
            f"-------------Computing 1st round call stats (Ancestry: {ANCESTRY.upper()})-----------"
        )
        mt = get_filtered_mt(mt_type = 'ACAF', sample_ids=defined_ht, filter_variants=True,
                            adj_filter=True, ancestry=ANCESTRY, prune_samples=False)
        call_stats_ht = mt.annotate_rows(call_stats=hl.agg.call_stats(mt.GT, mt.alleles)).rows()
        call_stats_ht = call_stats_ht.select('call_stats')
        call_stats_ht = call_stats_ht.naive_coalesce(1000).checkpoint(callstats_ht_path, overwrite=True)
        call_stats_ht.describe()
    call_stats_ht = hl.read_table(callstats_ht_path) 

    if args.batch:
        sys.exit()   

    pca_mt_path = f'{DATA_PATH}/utils/pca/pca_{ANCESTRY.upper()}.mt'
    related_mt_path = f'{DATA_PATH}/utils/pca/pca_{ANCESTRY.upper()}_related.mt'
    if not hfs.exists(f'{pca_mt_path}/_SUCCESS') or args.overwrite_pca_mt:
        print(
        f"-------------Exporting MT for PCA (Ancestry: {ANCESTRY.upper()})-----------"
    )
        mt = get_filtered_mt(mt_type = 'ACAF', sample_ids=defined_ht, filter_variants=True,
                            adj_filter=True, ancestry=ANCESTRY, prune_samples=False)
        print("------Removing HLA & inversion regions------")
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
        MIN_CR = CALLRATE_CUTOFF
        print(f'------Call rate filter: {MIN_CR}------')
        print(f'------AF filter: {MIN_AF}------')
        print(f'------Number of samples: {N_SAMPLES_HARD_FILTERED[ANCESTRY]}------')
        variants_to_keep = call_stats_ht.filter(
            (call_stats_ht.locus.in_autosome()) &
            (hl.is_snp(call_stats_ht.alleles[0], call_stats_ht.alleles[1])) &
            (call_stats_ht.call_stats.AF[1] >= MIN_AF) &
            ((call_stats_ht.call_stats.AN >= (N_SAMPLES_HARD_FILTERED[ANCESTRY] * 2 * MIN_CR)))
        )
        print('Filtering Variants...')
        mt = mt.filter_rows(hl.is_defined(variants_to_keep[mt.row_key])) # filter to high quality variants

        print(
            f"-------------Exporting related MT for PCA (Ancestry: {ANCESTRY.upper()})-----------"
        )
        print(related_mt_path)
        related_mt = mt.filter_cols(hl.is_defined(relatedness_flag_sample[mt.col_key]))
        related_mt = related_mt.naive_coalesce(5000).checkpoint(related_mt_path, overwrite=True)
        print(related_mt.count())
        
        print(f"-------------Exporting unrelated MT for PCA (Ancestry: {ANCESTRY.upper()})-----------")
        mt = mt.filter_cols(~hl.is_defined(relatedness_flag_sample[mt.col_key])) # filter to unrelated samples -> later to project
        mt = mt.unfilter_entries()
        n_partitions = 250 if ANCESTRY.lower() != 'eur' else 1000
        mt = mt.naive_coalesce(n_partitions).checkpoint(pca_mt_path, overwrite=True)
    mt = hl.read_matrix_table(pca_mt_path)
    print(mt.count())

    small_pca_mt_path = f'{DATA_PATH}/utils/pca/pca_{ANCESTRY.upper()}_downsampled.mt'
    if not hfs.exists(f'{small_pca_mt_path}/_SUCCESS') or args.overwrite_downsampled_mt:
        print(
        f"-------------Exporting Downsampled MT for PCA (Ancestry: {ANCESTRY.upper()})-----------"
        )
        mt = hl.read_matrix_table(pca_mt_path)
        print(mt.count())

        prop = round(1000000 / mt.count_rows(), 5)
        print(f'Proportion of variants: {prop}')
        print(f'Number of variants: {prop * mt.count_rows()}')
        mt = mt.sample_rows(prop)
        mt = mt.checkpoint(small_pca_mt_path, overwrite=True)
    mt = hl.read_matrix_table(small_pca_mt_path)
    print(mt.count())

    ld_prune_ht_path = f'{DATA_PATH}/utils/pca/pca_{ANCESTRY.upper()}_ld_pruned.ht'
    if not hfs.exists(f'{ld_prune_ht_path}/_SUCCESS') or args.overwrite_ld_prune_ht:
        print(f"-------------Exporting the LD pruned downsampled variant HT (ANCESTRY: {ANCESTRY.upper()})-----------")
        mt = hl.read_matrix_table(small_pca_mt_path)
        print(mt.count())

        ht = hl.ld_prune(mt.GT,
                        r2=0.1,
                        bp_window_size=int(1e7),
                        block_size=1024,
                        )
        ht = ht.checkpoint(ld_prune_ht_path, overwrite=True)
    ht = hl.read_table(ld_prune_ht_path)
    print(ht.count())

    ld_prune_mt_path = f'{DATA_PATH}/utils/pca/pca_{ANCESTRY.upper()}_ld_pruned.mt'
    repartition_mt_path = ld_prune_mt_path.replace('.mt', '_50partitions.mt')
    if not hfs.exists(f'{repartition_mt_path}/_SUCCESS') or args.overwrite_ld_prune_mt:
        print(f"-------------Exporting the LD pruned downsampled variant MT (ANCESTRY: {ANCESTRY.upper()})-----------")
        overwrite = not hfs.exists(f'{ld_prune_mt_path}/_SUCCESS')

        mt = mt.filter_rows(hl.is_defined(ht[mt.row_key]))
        mt = mt.checkpoint(ld_prune_mt_path, overwrite=overwrite, _read_if_exists=not overwrite)
        print(f'-----------Exporting LD-pruned & repartitioned MT (ANCESTRY: {ANCESTRY.upper()})-----------')
        partition = 50 if ANCESTRY.lower() != 'eur' else 100
        mt = mt.naive_coalesce(partition).checkpoint(repartition_mt_path, overwrite=True)
    mt = hl.read_matrix_table(repartition_mt_path)
    print(mt.count())


    pca_pruned_sample_tsv = f'{DATA_PATH}/utils/pca/results/aou_{ANCESTRY.lower()}_centroid_pruned.tsv'
    pca_pruned_sample_path = f'{DATA_PATH}/utils/pca/results/{ANCESTRY.lower()}_pca_centroid_pruned.ht'
    if hfs.exists(pca_pruned_sample_tsv):
        if not hfs.exists(f'{pca_pruned_sample_path}/_SUCCESS') or args.overwrite_pca_pruned_ht:
            pruned_ht = hl.import_table(pca_pruned_sample_tsv, delimiter='\t', key='s', types = {'s':hl.tstr})
            pruned_ht = pruned_ht.checkpoint(pca_pruned_sample_path, overwrite=True)
        else: 
            pruned_ht = hl.read_table(pca_pruned_sample_path)
        print(pruned_ht.count())


    pca_pruned_mt_path = f'{DATA_PATH}/utils/pca/pca_{ANCESTRY.upper()}_sample_pruned.mt'
    related_pca_pruned_mt_path = f'{DATA_PATH}/utils/pca/pca_{ANCESTRY.upper()}_related_sample_pruned.mt'
    PRUNE_SAMPLE = hfs.exists(f'{pca_pruned_sample_path}/_SUCCESS')
    TAG = ''
    pca_mt = hl.read_matrix_table(repartition_mt_path)
    print(f'Number of unrelated samples in PCA MT: {pca_mt.count()}') 
    related_mt = hl.read_matrix_table(related_mt_path)
    print(f'Number of related samples in PCA MT: {related_mt.count()}')
    if PRUNE_SAMPLE:
        TAG = '_pruned'
        if not hfs.exists(f'{pca_pruned_mt_path}/_SUCCESS') or args.overwrite_pca_pruned_mt:
            print(f'-----------Filtering to PCA pruned samples (ANCESTRY: {ANCESTRY.upper()})-----------')
            TAG = '_pruned'
            pruned_ht = hl.read_table(pca_pruned_sample_path)
            pca_mt = pca_mt.filter_cols(hl.is_defined(pruned_ht[pca_mt.col_key]))
            pca_mt = pca_mt.checkpoint(pca_pruned_mt_path, overwrite=True) 
            print(pca_mt.count())
            related_mt = related_mt.filter_cols(hl.is_defined(pruned_ht[related_mt.col_key]))
            related_mt = related_mt.checkpoint(related_pca_pruned_mt_path, overwrite=True)
            print(related_mt.count())
        pca_mt = hl.read_matrix_table(pca_pruned_mt_path)
        print(f'Number of unrelated samples in pruned PCA MT: {pca_mt.count()}')
        related_mt = hl.read_matrix_table(related_pca_pruned_mt_path)
        print(f'Number of related samples in pruned PCA MT: {related_mt.count()}')

    pca_results_ht_root = f'{DATA_PATH}/utils/pca/results/pca_{ANCESTRY.upper()}{TAG}'
    print(pca_results_ht_root)

    pca_unrelated_scores_ht_path = f'{pca_results_ht_root}_scores_unrelated.ht'
    overwrite_1=not hfs.exists(f'{pca_unrelated_scores_ht_path}/_SUCCESS')

    pca_unrelated_loadings_ht_path = f'{pca_results_ht_root}_loadings_unrelated.ht'
    overwrite_2=not hfs.exists(f'{pca_unrelated_loadings_ht_path}/_SUCCESS')

    if overwrite_1 or overwrite_2 or args.overwrite_unrelated_pca_ht:
        print(pca_mt.count()) 
        print(f'-----------Running PCA  (ANCESTRY: {ANCESTRY.upper()})-----------')
        eigenvalues, scores, loadings = hl.hwe_normalized_pca(
            pca_mt.GT,
            compute_loadings=True,
            k=50,
        )
        
        print(f'-----------Exporting PCA Hail Tables (ANCESTRY: {ANCESTRY.upper()})-----------')
        scores_ht = scores.checkpoint(pca_unrelated_scores_ht_path, overwrite=True)
        scores_ht.describe()

        loadings_ht=loadings.checkpoint(pca_unrelated_loadings_ht_path, overwrite=True)
        loadings_ht.describe()

    pca_related_scores_ht_path = f'{pca_results_ht_root}_scores_related.ht'
    if not hfs.exists(f'{pca_related_scores_ht_path}/_SUCCESS') or args.overwrite_related_pca_ht: 
        print(f'-----------Projecting PCA for related samples (ANCESTRY: {ANCESTRY.upper()})-----------')
        print(pca_related_scores_ht_path)
        pca_mt = pca_mt.annotate_rows(pca_af=hl.agg.mean(pca_mt.GT.n_alt_alleles()) / 2)

        pca_loadings = hl.read_table(pca_unrelated_loadings_ht_path)
        pca_loadings = pca_loadings.annotate(
            pca_af=pca_mt.rows()[pca_loadings.key].pca_af
        )
        related_scores_ht = hl.experimental.pc_project(
            related_mt.GT,
            pca_loadings.loadings,
            pca_loadings.pca_af,
        )
        related_scores_ht = related_scores_ht.checkpoint(pca_related_scores_ht_path, overwrite=True)
    related_scores_ht = hl.read_table(pca_related_scores_ht_path)
    print(related_scores_ht.count())


    full_pca_score_ht_path = f'{pca_results_ht_root}_scores_full.ht'
    if not hfs.exists(f'{full_pca_score_ht_path}/_SUCCESS') or args.overwrite_full_pca_ht:
        print(f"-------------Exporting Full PCA scores (ANCESTRY: {ANCESTRY.upper()})-----------")
        print(full_pca_score_ht_path)
        related_scores_ht = hl.read_table(pca_related_scores_ht_path)
        print(f"Number of related samples : {related_scores_ht.count()}")
        unrelated_scores_ht = hl.read_table(pca_unrelated_scores_ht_path)
        full_scores_ht = related_scores_ht.union(unrelated_scores_ht)
        full_scores_ht = full_scores_ht.annotate(**{f"PC{i + 1}": full_scores_ht.scores[i] for i in range(50)})

        full_scores_ht = full_scores_ht.checkpoint(full_pca_score_ht_path, overwrite=True)
        full_scores_ht.describe()
        
    full_scores_ht = hl.read_table(full_pca_score_ht_path)
    print(full_scores_ht.count())
    full_scores_ht.export(full_pca_score_ht_path.replace('.ht','.txt.bgz'))

    
    

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--ancestry', type=str)
    parser.add_argument('--overwrite-hard-filters', action='store_true')
    parser.add_argument('--overwrite-call-stats', action='store_true')
    parser.add_argument('--batch', action='store_true')
    parser.add_argument('--overwrite-pca-mt', action='store_true')
    parser.add_argument('--overwrite-downsampled-mt', action='store_true')
    parser.add_argument('--overwrite-ld-prune-ht', action='store_true')
    parser.add_argument('--overwrite-ld-prune-mt', action='store_true')
    parser.add_argument('--overwrite-pca-pruned-ht', action='store_true')
    parser.add_argument('--overwrite-pca-pruned-mt', action='store_true')
    parser.add_argument('--overwrite-unrelated-pca-ht', action='store_true')
    parser.add_argument('--overwrite-related-pca-ht', action='store_true')
    parser.add_argument('--overwrite-full-pca-ht', action='store_true')
    args = parser.parse_args()
    main(args)