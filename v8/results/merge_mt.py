import hail as hl
from collections import Counter
import pandas as pd
import hailtop.batch as hb
import hailtop.fs as hfs
from hailtop.batch.resource import Resource, ResourceGroup
from hailtop.batch.job import Job
from hailtop.batch.batch import Batch
from collections import Counter
from shlex import quote as shq
from datetime import date
from tqdm import tqdm
import uuid
import argparse

from aou_gwas import *  # for dataproc
# from utils.utils import *  # for QoB
# from utils.resources import *  # for QoB
# from utils.results_loading import * # for QoB

ANCESTRIES = ['AFR', 'AMR', 'EAS', 'EUR', 'MID', 'SAS']
TMP_BUCKET = 'gs://aou_tmp_30_days'

def read_pickle_dict(dict_path: str):
    dict = pd.read_pickle(dict_path)
    return dict

def unify_saige_burden_ht_schema(ht, pop='other', sex_specific: bool = False):
    # ht.describe()
    ht = ht.key_by()
    ht = ht.annotate(max_MAF = hl.if_else(hl.is_missing(ht.max_MAF), -1, ht.max_MAF))
    ht = ht.key_by('gene_id', 'gene_symbol', 'annotation', 'max_MAF', 'phenoname')
    if pop != 'meta':
        shared = ( 'Pvalue', 'Pvalue_Burden', 'Pvalue_SKAT', 'Pvalue_log10', 'Pvalue_Burden_log10', 'Pvalue_SKAT_log10', 'CHR', 'POS',
                  'BETA_Burden', 'SE_Burden','MAC',  'Number_rare', 'Number_ultra_rare', 'total_variants', 'interval')
        new_ints = ('MAC_case', 'MAC_control')
        if 'MAC_case' not in list(ht.row):
            ht = ht.select(*shared, **{field: hl.missing(hl.tint32) for field in new_ints})
        else:
            ht = ht.select(*shared, *new_ints)
    else:
        shared = ('Pvalue', 'Pvalue_Burden', 'Pvalue_SKAT', 'Pvalue_log10', 'Pvalue_Burden_log10', 'Pvalue_SKAT_log10',
                  'META_Stats_SKATO', 'META_Stats_SKAT', 'META_Stats_Burden')
        new_ints = ('CHR', 'POS')
        if 'interval' not in list(ht.row_value):
            ht = ht.select(*shared, CHR=hl.missing(hl.tstr), POS = hl.missing(hl.tint32))
        else:
            ht = ht.select(*shared, *new_ints)
    ht = ht.annotate(MAC_control = hl.int32(ht.MAC_control))
    if sex_specific:
        print('chrXXXXXXXXXX')
        ht = ht.filter(ht.interval.start.contig == 'chrX')
    # print(ht.aggregate(hl.agg.counter(ht.annotation)))
    return ht

def main(args):
    today = date.today().strftime("%y%m%d")
    analysis_type = "variant" if args.single_variant_only else "gene"
    test_type = 'saige' if analysis_type == 'variant' else 'saige_gene'
    result_type = "gene" if args.load_gene_results else "variant"
    inner_mode = 'overwrite' if args.overwrite else '_read_if_exists'
    variant_type = 'exome' if result_type == 'gene' and analysis_type == 'variant' else 'genome'
    print(f'INNER MODE: {inner_mode}')
    ancestries = args.ancestries.split(',') if args.ancestries is not None else ANCESTRIES
    hl.init(
        tmp_dir=TMP_BUCKET,
        driver_memory="highmem",
        driver_cores=8,
        worker_memory="highmem",
        worker_cores=1,
        default_reference="GRCh38",
        log="/aou_saige.log",
        app_name=f'aou_merge_SAIGE_{analysis_type}_results' if args.merge_hail_tables else f'aou_{analysis_type}_meta_analysis'
    )
    # TODO: look at the pan ancestry repo, add merge populations
    # TODO: format the table for meta-analysis on SAIGE-GENE results
    print(f'Analysis type: {analysis_type}')
    print(f'Result type: {result_type}')
    RESULT_ROOT = 'gs://aou_analysis/v8/ht_results'
    if args.merge_hail_tables:
        for ancestry in ancestries:
            print(ancestry)
            suffix = '' if not args.sex_specific else '_sex_specific'
            all_phenos_dir = hl.hadoop_ls(f'{RESULT_ROOT}/{ancestry.upper()}')
            if args.sex_specific:
                all_female_phenos_dir = hl.hadoop_ls(f'{RESULT_ROOT}_female/{ancestry.upper()}')
                all_male_phenos_dir = hl.hadoop_ls(f'{RESULT_ROOT}_male/{ancestry.upper()}')
                all_phenos_dir = all_female_phenos_dir + all_male_phenos_dir + all_phenos_dir
            print(all_phenos_dir[0:5])
            if args.category is not None:
                suffix = f'{suffix}.{args.category}'
                phenos_to_run_by_pop_by_group = read_pickle_dict('gs://aou_analysis/v8/data/phenotype/summary/pheno_by_ancestry_by_group_dict_both.dict')
                if args.category == 'quantitative_traits':
                    phenotypes = set(
                        phenos_to_run_by_pop_by_group[ancestry]["lab_measurement"] + phenos_to_run_by_pop_by_group[ancestry][
                            "physical_measurement"])
                else:
                    phenotypes = phenos_to_run_by_pop_by_group[ancestry][args.category]

            if analysis_type == 'gene' and result_type == 'gene':
                if args.category is not None:
                    all_gene_outputs = [f'{RESULT_ROOT}/{ancestry.upper()}/phenotype_{phenoname}/gene_results.ht' for phenoname in phenotypes
                                        if hfs.exists(f'{RESULT_ROOT}/{ancestry.upper()}/phenotype_{phenoname}/gene_results.ht/_SUCCESS')]
                else:
                    all_gene_outputs = get_files_in_parent_directory(all_phenos_dir, 'gene_results.ht')
                
                print(f'Got {len(all_gene_outputs)} HT paths...')
                print(all_gene_outputs[0:5])
                if ancestry == 'meta':
                    ht_keys = ['gene_id', 'gene_symbol', 'annotation', 'max_MAF']
                else:
                    ht_keys = ['gene_id', 'gene_symbol', 'annotation', 'max_MAF', 'phenoname']

                # all_hts = list(map(lambda x: unify_saige_burden_ht_schema(hl.read_table(x).key_by(*ht_keys), pop=pop), all_gene_outputs))
                all_hts = [unify_saige_burden_ht_schema(hl.read_table(x).key_by(*ht_keys), pop=ancestry, sex_specific=args.sex_specific) for x in tqdm(all_gene_outputs)]

                row_keys = ['gene_id', 'gene_symbol', 'annotation', 'max_MAF']
                col_keys = ['phenoname']

                mt = join_pheno_hts_to_mt(all_hts=all_hts, row_keys=row_keys, col_keys=col_keys, temp_dir=f'{TMP_BUCKET}/gene_{ancestry}{suffix}' , inner_mode='overwrite')
                mt = mt.drop('CHR', 'POS', 'Pvalue_log10', 'Pvalue_Burden_log10', 'Pvalue_SKAT_log10')

                mt.describe()
                mt.max_MAF.show()

                if ancestry == 'meta':
                    shared_fields = ['CHR', 'POS']
                else:
                    shared_fields = ['interval']
                    mt=mt.transmute_entries(total_variants_pheno=mt.total_variants)

                mt = pull_out_fields_from_entries( mt=mt,
                                                   shared_fields=shared_fields,
                                                   index='rows',
                                                   agg_funcs=[lambda x: hl.agg.take(x, 1)[0]] * len(shared_fields))
                mt.describe()
                mt_path = f'{RESULT_ROOT.replace("ht_results", "mt_results")}/{ancestry.upper()}_gene{suffix}.mt'
                print(f'SAVING {ancestry}{suffix} gene MT to {mt_path}')
                mt = mt.naive_coalesce(1000).checkpoint(mt_path,
                                                        overwrite=args.overwrite, _read_if_exists=not args.overwrite)
                mt.describe()
                print(mt.count())
                print(mt.aggregate_rows(hl.agg.counter(mt.max_MAF)))
                print(mt.aggregate_rows(hl.agg.counter(mt.annotation)))


            else:
                if args.category is not None:
                    all_variant_outputs = [f'{RESULT_ROOT}/{ancestry.upper()}/phenotype_{phenoname}/{variant_type}_variant_results.ht' for phenoname in phenotypes
                                          if hfs.exists(f'{RESULT_ROOT}/{ancestry.upper()}/phenotype_{phenoname}/{variant_type}_variant_results.ht/_SUCCESS')]
                else:
                    all_variant_outputs = get_files_in_parent_directory(all_phenos_dir, f'{variant_type}_variant_results.ht')

                print(f'Got {len(all_variant_outputs)} HT paths...')
                print(all_variant_outputs[0:5])
                all_hts = list(map(lambda x: unify_saige_ht_variant_schema(hl.read_table(x), sex_specific=args.sex_specific), all_variant_outputs)) # Add back

                # all_hts = list(map(lambda x: hl.read_table(x), all_variant_outputs))

                row_keys = ['locus', 'alleles']
                col_keys = ['phenoname']
                # inner_mode = 'overwrite' if args.overwrite else '_read_if_exists'
                mt = join_pheno_hts_to_mt(all_hts=all_hts, row_keys=row_keys, col_keys=col_keys, temp_dir=f"{TMP_BUCKET}/{variant_type}_variant_{ancestry}{suffix}", repartition_final=20000, inner_mode='overwrite')

                mt.describe()
                if ancestry != 'meta':
                    mt = pull_out_fields_from_entries(mt, ['MarkerID'], 'rows',
                                                      agg_funcs=[lambda x: hl.agg.take(x, 1)[0]] * 1)

                mt.describe()
                mt_path = f'{RESULT_ROOT.replace("ht_results", "mt_results")}/{ancestry.upper()}_{variant_type}_variant{suffix}.mt'
                print(f'SAVING {ancestry}{suffix} {variant_type} {result_type} MT to {mt_path}')

                mt = mt.checkpoint(mt_path,
                    overwrite=args.overwrite, _read_if_exists=not args.overwrite)
                mt.describe()
                print(mt.count())

    if args.merge_hail_matrix_tables:
        def write_full_mt(analysis_type, sex_specific, overwrite):
            suffix = '_sex_specific' if sex_specific else ''
            mts = []
            print(ANCESTRIES)
            for ancestry in ANCESTRIES:
                if analysis_type == 'variant':
                    mt_path = f'{RESULT_ROOT.replace("ht_results", "mt_results")}/{ancestry.upper()}_{variant_type}_variant{suffix}.mt'
                else:
                    mt_path = f'{RESULT_ROOT.replace("ht_results", "mt_results")}/{ancestry.upper()}_gene{suffix}.mt'
                print(mt_path)
                mt = hl.read_matrix_table(mt_path)
                mt = mt.annotate_cols(ancestry=ancestry)
                if analysis_type == 'gene':
                    mt = mt.annotate_entries(CAF = mt.MAC/(2*mt.N))
                    mt = mt.annotate_entries(two_pq = 2*mt.CAF*(1 - mt.CAF))
                    mt = mt.annotate_entries(Pvalue_defined = hl.int64(hl.is_defined(mt.Pvalue)),
                                            Pvalue_Burden_defined = hl.int64(hl.is_defined(mt.Pvalue_Burden)),
                                            Pvalue_SKAT_defined = hl.int64(hl.is_defined(mt.Pvalue_SKAT)))
                    mt = mt.annotate_entries(two_pq_Pvalue = mt.Pvalue_defined * mt.two_pq,
                                            two_pq_Pvalue_Burden = mt.Pvalue_Burden_defined * mt.two_pq,
                                            two_pq_Pvalue_SKAT = mt.Pvalue_SKAT_defined * mt.two_pq)
                print(mt.count())
                mt.describe()

                mt = mt.select_cols(pheno_data=mt.col_value)
                mt = mt.select_entries(summary_stats=mt.entry)
                mt = mt.collect_cols_by_key()

                def get_index(mt):
                    return hl.sorted(hl.enumerate(mt.pheno_data), key=lambda x: x[1].saige_version.split('_')[1],
                                     reverse=True)[0][0]

                mt = mt.select_entries(summary_stats=mt.summary_stats[get_index(mt)])
                mt = mt.select_cols(pheno_data=mt.pheno_data[get_index(mt)])
                mts.append(mt)

            full_mt = mts[0]
            for mt in mts[1:]:
                full_mt = full_mt.union_cols(mt, row_join_type='outer')

            full_mt.describe()

            full_mt = full_mt.collect_cols_by_key()
            full_mt.pheno_data.ancestry.show()
            if analysis_type == 'variant':
                full_path = f'{RESULT_ROOT.replace("ht_results", "mt_results")}/FULL_{variant_type}_variant{suffix}.mt'
            else:
                full_path = f'{RESULT_ROOT.replace("ht_results", "mt_results")}/FULL_gene{suffix}.mt'
            full_mt = full_mt.filter_rows(full_mt.annotation != 'Cauchy')
            print(f'SAVING {full_path}')
            full_mt = full_mt.checkpoint(full_path, overwrite=overwrite, _read_if_exists=not overwrite)
            full_mt.describe()
            print(full_mt.count())
            print('Ancestries per pheno:')
            print(dict(Counter(full_mt.aggregate_cols(hl.agg.counter(hl.len(full_mt.pheno_data))))))

        write_full_mt(analysis_type, args.sex_specific, overwrite=args.overwrite)

    if args.run_meta_analysis:
        suffix = '_sex_specific' if args.sex_specific else ''
        category = f'{args.category}.' if args.category is not None else ''
        result_type = 'gene' if args.load_gene_results and (analysis_type != 'variant') else 'variant'
        print(f'Running {test_type} {result_type} meta-analysis')
        if analysis_type == 'variant':
            mt_path = f'{RESULT_ROOT.replace("ht_results", "mt_results")}/FULL_{variant_type}_variant{suffix}.mt'
        else:
            mt_path = f'{RESULT_ROOT.replace("ht_results", "mt_results")}/FULL_gene{suffix}.mt'   
        print(mt_path)
        mt = hl.read_matrix_table(mt_path)
        mt.describe()
        # suffix += f'_CAF'
        if analysis_type == 'variant':
            meta_path = f'{RESULT_ROOT.replace("ht_results", "mt_results")}/META_{variant_type}_variant{suffix}.mt'
        else:
            meta_path = f'{RESULT_ROOT.replace("ht_results", "mt_results")}/META_gene{suffix}.mt'
        print(meta_path)
        if result_type == 'variant':
            mt = run_meta_analysis(mt, beta_field='BETA', se_field = 'SE', af_allele2_field = 'AF_Allele2', af_case_field = 'AF_case', af_ctrl_field = 'AF_ctrl')
            mt.describe()
            mt.checkpoint(meta_path, _read_if_exists=not args.overwrite, overwrite = args.overwrite)
        else:
            mt = mt.filter_rows(hl.is_defined(mt.max_MAF))
            mt = run_stouffer(mt)
            mt.describe()

            mt = mt.checkpoint(meta_path, overwrite=args.overwrite, _read_if_exists=not args.overwrite)

            mt.META_Pvalue_Burden.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--single-variant-only", help="Run single variant SAIGE", action="store_true"
    )
    parser.add_argument(
        "--merge-hail-tables", help="Merge hail tables", action="store_true"
    )
    parser.add_argument(
        "--merge-hail-matrix-tables", help="Merge hail matrix tables", action="store_true"
    )
    parser.add_argument(
        "--run-meta-analysis", help="run-meta-analysis", action="store_true"
    )
    parser.add_argument(
        "--create-meta-analysis-ht", help="Generate meta analysis results in HT", action="store_true"
    )
    parser.add_argument(
        "--annotate-expected-pvalue", help="Annotate expected pvalue to HT", action="store_true"
    )
    parser.add_argument(
        "--merge-call-stats-ht", help="Merge ACAF and Exome HTs", action="store_true"
    )
    parser.add_argument(
        "--run-additional-load", help="Merge ACAF and Exome HTs", action="store_true"
    )
    parser.add_argument(
        "--force-reload", help="Force reload phenotypes", action="store_true"
    )
    parser.add_argument(
        "--load-only", help="Load only", action="store_true"
    )
    parser.add_argument(
        "--dry-run", help="Dry run", action="store_true"
    )
    parser.add_argument(
        "--test", help="Test on just one job", action="store_true"
    )
    parser.add_argument(
        "--ancestries", help="comma-separated list", default="afr,amr,eas,eur,mid,sas,meta"
    )
    parser.add_argument(
        "--phenos", help="comma-separated list", default=""
    )
    parser.add_argument(
        "--overwrite",
        help="Overwrite results",
        action="store_true",
    )
    parser.add_argument(
        "--load-variant-results", help="Load single-variant level test results", action="store_true"
    )
    parser.add_argument(
        "--load-gene-results", help="Load gene level test results", action="store_true"
    )
    parser.add_argument(
        "--category",
        help="phenotype category to export",
    )
    parser.add_argument(
        "--sex-specific", help="Sex-specific analysis", action="store_true"
    )
    args = parser.parse_args()

    main(args)
