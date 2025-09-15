import hail as hl
import hailtop.fs as hfs
import hailtop.batch as hb
from hailtop.batch.job import Job
from hailtop.batch.batch import Batch
from typing import Union
import argparse
import sys

SAIGE_DOCKER_IMAGE = "wzhou88/saige:1.4.4"
TMP_BUCKET = 'gs://aou_tmp/v8'
DATA_PATH = 'gs://aou_analysis/v8/data'
EXOME_MT_PATH = f'gs://fc-aou-datasets-controlled/v8/wgs/short_read/snpindel/exome/splitMT/hail.mt'
ACAF_MT_PATH =  f'gs://fc-aou-datasets-controlled/v8/wgs/short_read/snpindel/acaf_threshold/splitMT/hail.mt'


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


def get_filtered_mt(mt_type: hl.tstr,
                    sample_ids: hl.Table,
                    ancestry: str='all',
                    filter_samples: bool=True, 
                    filter_variants: bool=True,
                    prune_samples:bool=True,
                    adj_filter: bool=True):
    if ancestry != 'all':
        anc_mt_path = f'{DATA_PATH}/utils/raw_mt/{mt_type}/{ancestry.upper()}_{mt_type}.mt'
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
        print(mt.count())
    else:
        mt_path = ACAF_MT_PATH if mt_type=='ACAF' else EXOME_MT_PATH
        mt = hl.read_matrix_table(mt_path)
        mt = mt.filter_entries(hl.is_missing(mt.FT) | (mt.FT == 'PASS'))
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
    ht: hl.Table, # sample pruned call stats HT fro Exome MT
    ancestry: str,
    n_samples: int, 
    n_common_variants_to_keep: int=50000, # 100000 for per pop
    min_call_rate: float = 0.9,
    min_maf_common_variants: float = 0.01,
    variants_per_mac_category: int = 2000,
    variants_per_maf_category: int = 10000,
):
    print(f'Number of common variants to sample: {n_common_variants_to_keep}')
    ht = ht.filter(
        (ht.locus.in_autosome())
        & (ht.call_stats.AN >= (n_samples * 2 * min_call_rate))
        & (ht.call_stats.AC[1] > 0)
    )

    ht = ht.annotate(
        mac_category=mac_category_case_builder(ht.call_stats.AC[1], ht.call_stats.AF[1], min_maf_common_variants)
    )

    # From: https://hail.zulipchat.com/#narrow/stream/123010-Hail-Query-0.2E2-support/topic/.E2.9C.94.20randomly.20sample.20table/near/388162012
    bins = ht.aggregate(hl.agg.collect_as_set(ht.mac_category))
    ac_bins = [bin for bin in bins if (bin >= 1) and (bin<=5)]
    af_bins = [bin for bin in bins if (bin < 0.99) or (bin > 5)]

    sampled_common_variants = ht.aggregate(
        hl.agg.filter(
                ht.call_stats.AF[1] > min_maf_common_variants,
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

def create_sparse_grm(
    p: Batch,
    output_path: str,
    plink_file_root: str,
    docker_image: str,
    relatedness_cutoff: str = "0.125",
    num_markers: int = 2000,
    n_threads: int = 16,
    memory: str = "highmem",
    storage="1500Mi",
):
    in_bfile = p.read_input_group(
        **{ext: f"{plink_file_root}.{ext}" for ext in ("bed", "bim", "fam")}
    )
    create_sparse_grm_task: Job = p.new_job(name="create_sparse_grm")
    create_sparse_grm_task.cpu(n_threads).storage(storage).image(docker_image).memory(
        memory
    )
    create_sparse_grm_task._preemptible = False
    create_sparse_grm_task.declare_resource_group(
        sparse_grm={
            ext: f"{{root}}{ext}"
            for ext in (
                f"_relatednessCutoff_{relatedness_cutoff}_{num_markers}_randomMarkersUsed.sparseGRM.mtx",
                f"_relatednessCutoff_{relatedness_cutoff}_{num_markers}_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt",
            )
        }
    )
    # create_sparse_grm_task._machine_type = 'n1-highmem-32'
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
    ANCESTRY = args.ancestry if args.ancestry is not None else 'all'
    if args.prepare_grm_data:
        hl.init(
            app_name=f'Process_grm_{ANCESTRY}',
            gcs_requester_pays_configuration='aou-neale-gwas',
            tmp_dir=f"{TMP_BUCKET}/grm",
            default_reference="GRCh38",
            log="/process_grm.log"
        )

        sample_ids = hl.read_table(f'{DATA_PATH}/utils/aou_v8_hard_filter.ht')
        sub_mt = get_filtered_mt(
            mt_type='Exome',
            sample_ids=sample_ids,
            ancestry=ANCESTRY
        ) 


        call_stat_path = f'{DATA_PATH}/utils/call_stats/exome_pruned/{ANCESTRY.upper()}_exome_call_stats.ht'
        if not hfs.exists(f'{call_stat_path}/_SUCCESS') or args.overwrite_call_stats:
            call_stats_ht = sub_mt.annotate_rows(call_stats=hl.agg.call_stats(sub_mt.GT, sub_mt.alleles)).rows()
            call_stats_ht = call_stats_ht.naive_coalesce(1000).checkpoint(call_stat_path, overwrite=True)
        
        if args.batch:
            sys.exit()

        N_SAMPLES = sub_mt.count_cols()
        print(f"-------------Exporting downsampled variant HT (Ancestry: {ANCESTRY}| N: {N_SAMPLES})-----------")
        grm_variant_path = f'{DATA_PATH}/utils/grm/{ANCESTRY.upper()}_downsampled_variants.ht'
        if not hfs.exists(f'{grm_variant_path}/_SUCCESS') or args.overwrite_downsampled_variants:
            ht = hl.read_table(call_stat_path)
            filtered_ht = filter_ht_for_grm(
                ht=ht,
                ancestry=ANCESTRY,
                n_samples = N_SAMPLES,

            )

            filtered_ht.naive_coalesce(1000).checkpoint(
                grm_variant_path,
                overwrite=True,
            )
        filtered_ht = hl.read_table(grm_variant_path)
        print(f'Downsampled variant HT (pop: {ANCESTRY.upper()}) written to {grm_variant_path}')
        print(f'Number of variants sampled for {ANCESTRY.upper()}: {filtered_ht.count()}')


        grm_mt_path = f'{DATA_PATH}/utils/grm/{ANCESTRY.upper()}_downsampled_variants.mt'
        if not hfs.exists(f'{grm_mt_path}/_SUCCESS') or args.overwrite_downsampled_mt:
            print(
                f"-------------Exporting downsampled variant MT (Ancestry: {ANCESTRY.upper()})-----------"
            )
            print(f'Number of variants sampled for {ANCESTRY.upper()}: {filtered_ht.count()}')
            sub_mt = sub_mt.filter_rows(hl.is_defined(filtered_ht[sub_mt.row_key]))

            print("Removing HLA...")
            # Common inversion taken from Table S4 of https://www.ncbi.nlm.nih.gov/pubmed/27472961
            # (converted to GRCh38 by: https://liftover.broadinstitute.org/#input=chr8%3A8055789-11980649&hg=hg19-to-hg38 )
            # Also removing HLA, from https://www.ncbi.nlm.nih.gov/grc/human/regions/MHC?asm=GRCh38
            sub_mt = sub_mt.filter_rows(
                ~hl.parse_locus_interval(
                    "chr8:8198267-12123140", reference_genome="GRCh38"
                ).contains(sub_mt.locus)
                & ~hl.parse_locus_interval(
                    "chr6:28510120-33480577", reference_genome="GRCh38"
                ).contains(sub_mt.locus)
            )

            sub_mt.naive_coalesce(5000).checkpoint(
                grm_mt_path,
                overwrite=True,
            )
        mt = hl.read_matrix_table(grm_mt_path)
        mt.describe()
        print(mt.count())

        mt = hl.read_matrix_table(grm_mt_path)
        print(f'-------------{ANCESTRY.upper()} MT entering LD-pruning: {mt.count()}-------------')
        mt = mt.unfilter_entries()
        grm_ld_pruned_variant_path = f'{DATA_PATH}/utils/grm/{ANCESTRY.upper()}_ld_pruned_variants.ht'
        print(f"-------------Exporting the LD pruned downsampled variant HT (Ancestry: {ANCESTRY.upper()})-----------")
        if not hfs.exists(f'{grm_ld_pruned_variant_path}/_SUCCESS') or args.overwrite_ld_pruned_ht:
            ld_prune_ht = hl.ld_prune(mt.GT,
                            r2=0.1,
                            bp_window_size=int(1e7),
                            block_size=1024,
                            )
            ld_prune_ht = ld_prune_ht.checkpoint(
                grm_ld_pruned_variant_path,
                overwrite=True,
            )


        ht = hl.read_table(grm_ld_pruned_variant_path)
        mt = mt.filter_rows(hl.is_defined(ht[mt.row_key]))
        grm_plink_path = f'{DATA_PATH}/utils/grm/{ANCESTRY.upper()}_grm_plink'

        print(f"-------------Exporting variant downsampled plink files (Ancestry: {ANCESTRY.upper()})-----------")
        hl.export_plink(
            mt,
            grm_plink_path
        )

        print(f"-------------Exporting sample ID file (Ancestry: {ANCESTRY.upper()})-----------")
        mt = hl.read_matrix_table(grm_mt_path)
        grm_samples_path = f'{DATA_PATH}/utils/grm/{ANCESTRY.upper()}_grm_plink.samples'
        with hl.hadoop_open(grm_samples_path, 'w') as f:
            f.write('\n'.join(mt.s.collect()) + '\n')

    if args.export_sparse_grm:
        print(f"-------------Exporting GRM (Ancestry: {ANCESTRY.upper()})-----------")
        sparse_grm_root = f"{DATA_PATH}/utils/grm/aou_{ANCESTRY.lower()}"
        grm_plink_path = f'{DATA_PATH}/utils/grm/{ANCESTRY.upper()}_grm_plink'
        relatedness_cutoff = "0.125"
        num_markers = 2000
        n_threads = 16

        backend = hb.ServiceBackend(
                billing_project="all-by-aou",
                remote_tmpdir=TMP_BUCKET,
            )
        b = hb.Batch(
                name=f"Create_sparse_GRM_{ANCESTRY.lower()}",
                requester_pays_project="aou-neale-gwas",
                backend=backend,
            )

        create_sparse_grm(
                b,
                sparse_grm_root,
                grm_plink_path,
                SAIGE_DOCKER_IMAGE,
                relatedness_cutoff,
                num_markers,
                n_threads=n_threads,
            )

        b.run()
        

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process GRM data.')
    parser.add_argument(
        '--ancestry',
        type=str,
        required=True,
        help='Ancestry to process.'
    )
    parser.add_argument(
        '--batch',
        action='store_true',
        help='Run in batch mode.'
    )
    parser.add_argument(
        '--overwrite-call-stats',
        action='store_true',
        help='Overwrite call stats HT.'
    )
    parser.add_argument(
        '--overwrite-downsampled-variants',
        action='store_true',
        help='Overwrite downsampled variants HT.'
    )
    parser.add_argument(
        '--overwrite-downsampled-mt',
        action='store_true',
        help='Overwrite downsampled MT.'
    )
    parser.add_argument(
        '--overwrite-ld-pruned-ht',
        action='store_true',
        help='Overwrite LD pruned HT.'
    )
    parser.add_argument(
        '--export-sparse-grm',
        action='store_true',
        help='Export sparse GRM.'
    )
    parser.add_argument(
        '--prepare-grm-data',
        action='store_true',
        help='Prepare GRM data.'
    )
    
    args = parser.parse_args()
    main(args)