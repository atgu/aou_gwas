import hail as hl
import copy
from collections import Counter
import pandas as pd
import pickle
import argparse
import logging
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

from aou_gwas import *  # for dataproc


def get_files_in_parent_directory(parent_dir, fname: str = 'variant_results.ht'):
    all_outputs = []
    for directory in parent_dir:
        if not directory['is_dir']:
            continue
        file_path = f'{directory["path"]}/{fname}'
        if hfs.exists(f'{file_path}/_SUCCESS'):
            all_outputs.append(file_path)
    return all_outputs

def read_pickle_dict(dict_path: str):
    dict = pd.read_pickle(dict_path)
    return dict

def get_aou_saige_results_root(analysis_type: str, pop: str, name: str, test: bool = False):
    assert analysis_type in [
        "gene",
        "variant",
    ], "Name has to be from ['gene', 'variant']"
    assert name in [
        "bgen",
        "result",
        "pilot_result",
        "plot"
    ], "Name has to be from ['bgen', 'result', 'plot']"
    return f"{CURRENT_PATH}/{analysis_type}_results/{name}/{'test/' if test else ''}{pop.upper()}"



def main(args):
    today = date.today().strftime("%y%m%d")
    analysis_type = "variant" if args.single_variant_only else "gene"
    test_type = 'saige' if analysis_type == 'variant' else 'saige_gene'
    result_type = "gene" if args.load_gene_results else "variant"
    inner_mode = 'overwrite' if args.overwrite else '_read_if_exists'
    print(f'INNER MODE: {inner_mode}')
    pops = args.pops.split(',')
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
    if args.merge_hail_tables:
        for pop in pops:
            print(pop)
            all_phenos_dir = hl.hadoop_ls(f'{RESULTS_PATH}/{test_type}_results/{pop.upper()}')
            print(all_phenos_dir[0:5])
            suffix = ''
            if args.category is not None:
                suffix = f'{args.category}.'
                phenos_to_run_by_pop_by_group = read_pickle_dict(
                    get_phenotype_info_path(version="pheno_by_pop_by_group_dict", extension="dict"))
                if args.category == 'quantitative_traits':
                    phenotypes = set(
                        phenos_to_run_by_pop_by_group[pop]["lab_measurements"] + phenos_to_run_by_pop_by_group[pop][
                            "processed_physical_measurement_table"])
                else:
                    phenotypes = phenos_to_run_by_pop_by_group[pop][args.category]

            if analysis_type == 'gene' and result_type == 'gene':
                if args.category is not None:
                    all_gene_outputs = [f'{RESULTS_PATH}/{test_type}_results/{pop.upper()}/phenotype_{phenoname}/gene_results.ht' for phenoname in phenotypes
                                        if hfs.exists(f'{RESULTS_PATH}/{test_type}_results/{pop.upper()}/phenotype_{phenoname}/gene_results.ht/_SUCCESS')]
                else:
                    all_gene_outputs = get_files_in_parent_directory(all_phenos_dir, 'gene_results.ht')

                print(f'Got {len(all_gene_outputs)} HT paths...')
                print(all_gene_outputs[0:5])
                if pop == 'meta':
                    ht_keys = ['gene_id', 'gene_symbol', 'annotation', 'max_MAF']
                else:
                    ht_keys = ['gene_id', 'gene_symbol', 'annotation', 'max_MAF', 'phenoname']
                # all_hts = list(map(lambda x: unify_saige_burden_ht_schema(hl.read_table(x).key_by(*ht_keys), pop=pop), all_gene_outputs))
                all_hts = [unify_saige_burden_ht_schema(hl.read_table(x).key_by(*ht_keys), pop=pop) for x in tqdm(all_gene_outputs)]

                row_keys = ['gene_id', 'gene_symbol', 'annotation', 'max_MAF']
                col_keys = ['phenoname']

                mt = join_pheno_hts_to_mt(all_hts=all_hts, row_keys=row_keys, col_keys=col_keys, temp_dir=f'{TMP_BUCKET}/new_gene_{pop}' , inner_mode='overwrite')

                mt.describe()
                mt.max_MAF.show()

                if pop == 'meta':
                    shared_fields = ['interval', 'CHR', 'POS']
                else:
                    shared_fields = ['MAC',  'Number_rare', 'Number_ultra_rare', 'total_variants', 'interval', 'CHR', 'POS']
                    mt=mt.annotate_entries(total_variants_pheno=mt.total_variants)

                mt = pull_out_fields_from_entries( mt=mt,
                                                   shared_fields=shared_fields,
                                                   index='rows',
                                                   agg_funcs=[lambda x: hl.agg.take(x, 1)[0]] * len(shared_fields))
                mt = mt.filter_cols(~hl.literal(LAB_TO_REMOVE).contains(mt.phenoname))
                mt.describe()
                mt_path = get_aou_final_results_path(analysis_type=analysis_type, pop=pop, result_type='gene', extension=f'{suffix}mt')
                print(f'SAVING {pop} {test_type} {result_type} MT to {mt_path}')
                mt = mt.naive_coalesce(1000).checkpoint(get_aou_final_results_path(analysis_type=analysis_type, pop=pop, result_type='gene', extension=f'{suffix}mt'),
                                                        overwrite=args.overwrite, _read_if_exists=not args.overwrite)
                mt.describe()
                print(mt.count())
                print(mt.aggregate_rows(hl.agg.counter(mt.max_MAF)))


            else:
                if args.category is not None:
                    all_variant_outputs = [f'{RESULTS_PATH}/{test_type}_results/{pop.upper()}/phenotype_{phenoname}/variant_results.ht' for phenoname in phenotypes
                                          if hfs.exists(f'{RESULTS_PATH}/{test_type}_results/{pop.upper()}/phenotype_{phenoname}/variant_results.ht/_SUCCESS')]
                else:
                    all_variant_outputs = get_files_in_parent_directory(all_phenos_dir, 'variant_results.ht')
                print(f'Got {len(all_variant_outputs)} HT paths...')
                all_hts = list(map(lambda x: unify_saige_ht_variant_schema(hl.read_table(x)), all_variant_outputs))

                row_keys = ['locus', 'alleles']
                col_keys = ['phenoname']
                # inner_mode = '_read_if_exists'
                mt = join_pheno_hts_to_mt(all_hts=all_hts, row_keys=row_keys, col_keys=col_keys, temp_dir=f"{TMP_BUCKET}/{suffix}variant_{pop}", repartition_final=None, inner_mode=inner_mode)

                mt.describe()
                mt = pull_out_fields_from_entries(mt, ['MarkerID'], 'rows',
                                                  agg_funcs=[lambda x: hl.agg.take(x, 1)[0]] * 1)


                mt_path = get_aou_final_results_path(analysis_type=analysis_type, pop=pop, result_type='variant',
                                               extension=f'{suffix}mt')
                print(f'SAVING {pop} {test_type} {result_type} MT to {mt_path}')

                mt = mt.naive_coalesce(1000).checkpoint(mt_path,
                    overwrite=args.overwrite, _read_if_exists=not args.overwrite)
                mt.describe()
                print(mt.count())


    if args.merge_hail_matrix_tables:
        category = f'{args.category}.' if args.category is not None else ''
        def write_full_mt(category, result_type, overwrite):
            col_info_ht = hl.read_table('gs://aou_results/250k/v1.1/utils/aou_phenotype_meta_info.ht')
            mts = []
            for pop in POPS:
                print(get_aou_final_results_path(analysis_type=analysis_type, pop=pop, result_type=result_type, extension=f'{category}mt'))
                col_info_ht = col_info_ht.filter(col_info_ht.pop == pop).key_by('phenoname')
                mt = hl.read_matrix_table(get_aou_final_results_path(analysis_type=analysis_type, pop=pop, result_type=result_type, extension=f'{category}mt'))
                # if pop in ['eas', 'mid']: mt = mt.drop('log_pvalue')
                mt = mt.annotate_cols(**col_info_ht[mt.col_key]).drop('lambda_gc_exome', 'lambda_gc_acaf', 'lambda_gc_gene_burden_001')
                mt = mt.annotate_cols(pop=pop)
                mt = mt.select_cols('n_cases', 'n_controls', 'heritability', 'saige_version', 'inv_normalized', 'pheno_sex', 'trait_type',
                                    'category', 'description', 'phecode_category', 'description_more',  'pop')
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

            full_mt = full_mt.collect_cols_by_key()
            full_mt.pheno_data.pop.show()
            print(get_aou_final_results_path(analysis_type=analysis_type, pop='full', result_type=result_type, extension=f'{category}mt'))
            full_mt = full_mt.checkpoint(get_aou_final_results_path(analysis_type=analysis_type, pop='full', result_type=result_type, extension=f'{category}mt'), overwrite)
            full_mt.describe()
            print(full_mt.count())
            print('Pops per pheno:')
            print(dict(Counter(full_mt.aggregate_cols(hl.agg.counter(hl.len(full_mt.pheno_data))))))

        write_full_mt(category=category, result_type=result_type, overwrite=args.overwrite)

    if args.run_meta_analysis:
        category = f'{args.category}.' if args.category is not None else ''
        result_type = 'gene' if args.load_gene_results and (analysis_type != 'variant') else 'variant'
        print(f'Running {test_type} {result_type} meta-analysis')
        mt_path = get_aou_final_results_path(analysis_type=analysis_type, pop='full', result_type=result_type, extension=f'{category}mt')
        print(mt_path)
        mt = hl.read_matrix_table(mt_path)
        mt.describe()
        meta_path = get_aou_final_results_path(analysis_type=analysis_type, pop='meta', result_type=result_type, extension=f'{category}mt')
        print(meta_path)
        if result_type == 'variant':
            mt = run_meta_analysis(mt, beta_field='BETA', se_field = 'SE', af_allele2_field = 'AF_Allele2', af_case_field = 'AF_case', af_ctrl_field = 'AF_ctrl')
            mt.describe()
            mt.checkpoint(meta_path, _read_if_exists=not args.overwrite, overwrite = args.overwrite)
        else:
            mt = mt.filter_rows(hl.is_defined(mt.max_MAF))
            mt = run_stouffer(mt)
            mt.describe()

            mt = mt.checkpoint(meta_path, overwrite=args.overwrite)

            mt.META_Pvalue_Burden.show()

    if args.create_meta_analysis_ht:
        def run_meta_analysis_ht(pheno:str, test_type:str=test_type, result_type:str=result_type, overwrite:bool=True):
            hl.init(
                master='local[16]',
                tmp_dir="gs://aou_tmp",
                default_reference="GRCh38",
            )
            paths = hl.hadoop_ls(
                f'gs://aou_results/250k/v1.1/{test_type}_results/*/phenotype_{pheno}'
            )
            paths = [path['path'] for path in paths if ('META' not in path['path']) and (path['is_dir']) and (path['path'].endswith(f'{result_type}_results.ht'))]
            meta_pop = [path.split('/')[-3].lower() for path in paths]
            ht = hl.read_table(paths[0])
            ht = ht.annotate(pop = meta_pop[0])
            if result_type == 'variant':
                ht = ht.annotate(BETA = hl.float64(ht.BETA),
                                 SE = hl.float64(ht.SE))
            fields = list(ht.row_value)
            ht = ht.select(*fields)
            for path in paths[1:]:
                print(path)
                tmp_ht = hl.read_table(path)
                tmp_ht = tmp_ht.annotate(pop = path.split('/')[-3].lower())
                if result_type == 'variant':
                    tmp_ht = tmp_ht.annotate(BETA=hl.float64(tmp_ht.BETA),
                                             SE=hl.float64(tmp_ht.SE))
                tmp_ht = tmp_ht.select(*fields)
                ht = ht.union(tmp_ht)

            if result_type == 'variant':
                ht = ht.annotate(weight = 1 / (ht.SE ** 2))
                meta_ht = ht.group_by('locus', 'alleles').aggregate(
                    BETA=hl.agg.sum(ht.weight * ht.BETA) / hl.agg.sum(ht.weight),
                    SE=1 / hl.sqrt(hl.agg.sum(ht.weight))
                )
                meta_ht.describe()
                ht = ht.annotate(META_BETA=meta_ht[ht.locus, ht.alleles].BETA)
                ht.describe()
                q_ht = ht.group_by('locus', 'alleles').aggregate(
                    Q=hl.agg.sum(ht.weight * (ht.BETA - ht.META_BETA) ** 2)
                )
                q_ht.describe()
                meta_ht = meta_ht.annotate(Het_Q=q_ht[meta_ht.key].Q)

                meta_ht = meta_ht.annotate(Pvalue=2 * hl.pnorm(-hl.abs(meta_ht.BETA) / meta_ht.SE))
                meta_ht = meta_ht.annotate(Pvalue_log10=-hl.log10(meta_ht.Pvalue))
            else:
                ht = ht.annotate(max_MAF=hl.if_else(hl.is_missing(ht.max_MAF), -1, ht.max_MAF))
                pheno_info_path = 'gs://aou_results/250k/v1.1/utils/aou_phenotype_meta_info.ht'
                pheno_ht = hl.read_table(pheno_info_path)
                pheno_ht = pheno_ht.annotate( n_controls=hl.if_else(hl.is_missing(pheno_ht.n_controls), 0, pheno_ht.n_controls))
                sub_pheno_ht = pheno_ht.filter(pheno_ht.phenoname == pheno)
                ht = ht.drop('n_cases', 'n_controls')
                ht = ht.annotate(
                    **sub_pheno_ht[ht.phenoname, ht.pop]
                )
                P_FIELDS = ['Pvalue', 'Pvalue_Burden', 'Pvalue_SKAT']
                P_TESTS = {'SKATO': 'Pvalue', 'Burden': 'Pvalue_Burden', 'SKAT': 'Pvalue_SKAT'}
                def _edit_pvalue(p):
                    return hl.if_else(p > 0.99, 0.99, p)

                ht = ht.filter(
                    hl.is_defined(ht[P_FIELDS[0]]) &
                    hl.is_defined(ht[P_FIELDS[1]]) &
                    hl.is_defined(ht[P_FIELDS[2]]))
                ht = ht.annotate(
                    **{p_field: _edit_pvalue(ht[p_field]) for p_field in P_FIELDS},
                    N = hl.if_else(ht.n_controls ==0, ht.n_cases, (4 * ht.n_cases * ht.n_controls) / (ht.n_cases + ht.n_controls))
                )

                def _stouffer_test(ht, test):
                    print(f'Meta analyzing {test} results...')
                    two_tail = test == 'Burden'
                    beta_field = f'BETA_{test}'
                    p_field = P_TESTS[test]

                    if two_tail:
                        ht = ht.annotate(**{p_field: ht[p_field] / 2},
                                                 **{beta_field: ht[beta_field]})
                    else:
                        ht = ht.annotate(
                            **{beta_field: hl.int(hl.is_defined(ht[p_field]))})
                    ht = ht.annotate(**{f'weighted_Z_numerator_{test}': (hl.sqrt(ht.N) * (-hl.qnorm(ht[p_field])) * hl.sign(ht[beta_field]))})
                    meta_ht = ht.group_by('gene_id', 'gene_symbol', 'annotation', 'max_MAF').aggregate(
                        interval=hl.agg.collect(ht.interval),
                        **{f'META_Stats_{test}': hl.agg.sum(ht[f'weighted_Z_numerator_{test}']) / (hl.sqrt(hl.agg.sum(ht.N)))},
                    )

                    if two_tail:
                        meta_ht = meta_ht.annotate(**{f'{p_field}': 2 * hl.pnorm(hl.abs(meta_ht[f'META_Stats_{test}']), lower_tail=False)})
                    else:
                        meta_ht = meta_ht.annotate(**{f'{p_field}': hl.pnorm(meta_ht[f'META_Stats_{test}'], lower_tail=False)})
                    meta_ht = meta_ht.annotate(**{f'{p_field}_log10': -hl.log10(meta_ht[p_field])})
                    return meta_ht

                meta_ht1 = _stouffer_test(ht, 'SKATO')
                meta_ht1 = meta_ht1.checkpoint(f'{TMP_BUCKET}/meta_tmp_{pheno}_SKATO_{test_type}_{result_type}_{today}_12345.ht')
                meta_ht1.show(20)
                meta_ht2 = _stouffer_test(ht, 'SKAT').drop('interval')
                meta_ht2 = meta_ht2.checkpoint(
                    f'{TMP_BUCKET}/meta_tmp_{pheno}_SKAT_{test_type}_{result_type}_{today}_12345.ht')
                meta_ht2.show(20)
                meta_ht3 = _stouffer_test(ht, 'Burden').drop('interval')
                meta_ht3 = meta_ht3.checkpoint(
                    f'{TMP_BUCKET}/meta_tmp_{pheno}_Burden_{test_type}_{result_type}_{today}_12345.ht')
                meta_ht3.show(20)

                meta_ht = meta_ht1.annotate(
                    **meta_ht2[meta_ht1.key],
                    **meta_ht3[meta_ht1.key],
                )


            meta_ht = meta_ht.annotate_globals(meta_pop=meta_pop, phenoname=pheno)
            meta_ht.describe()
            meta_ht.show()
            if result_type == 'gene':
                meta_ht = meta_ht.annotate(CHR=meta_ht.interval[0].start.contig, POS=meta_ht.interval[0].start.position)
            else:
                meta_ht = meta_ht.annotate(CHR=meta_ht.locus.contig, POS=meta_ht.locus.position)
            meta_path = f'gs://aou_results/250k/v1.1/{test_type}_results/META/phenotype_{pheno}/{result_type}_results.ht'
            meta_ht = meta_ht.checkpoint(meta_path, overwrite=overwrite)
            meta_ht.describe()
            print(meta_ht.count())


        backend = hb.ServiceBackend(
            billing_project="all-by-aou",
            remote_tmpdir='gs://aou_tmp/'
        )

        b = hb.Batch(
            name=f"aou_export_meta_ht_{test_type}_{result_type}_result",
            backend=backend,
            default_storage="500Mi",
            default_cpu=8,
        )

        phenos_to_run_by_pop = read_pickle_dict(get_phenotype_info_path(version="pheno_by_pop_dict", extension="dict"))
        phenotypes = phenos_to_run_by_pop['meta']
        phenos_to_run_by_pop_by_group = read_pickle_dict(get_phenotype_info_path(version="pheno_by_pop_by_group_dict", extension="dict"))
        phenotypes = []
        for pop in POPS:
            phenotypes = phenotypes + phenos_to_run_by_pop_by_group[pop]['lab_measurements'] + phenos_to_run_by_pop_by_group[pop]['processed_physical_measurement_table'] + \
                         [pheno for pheno in phenos_to_run_by_pop_by_group[pop]['random_phenotypes'] if 'continuous' in pheno]
        phenotypes = list(set(phenotypes))
        print(phenotypes)
        for i in tqdm(range(len(phenotypes))):
            pheno = phenotypes[i]
            if not hl.hadoop_exists(f'gs://aou_results/250k/v1.1/{test_type}_results/META/phenotype_{pheno}/{result_type}_results.ht/_SUCCESS') or args.overwrite:
                j = b.new_python_job(name=f"aou_{pheno}_{test_type}_{result_type}_result")
                j.image("hailgenetics/hail:0.2.130-py3.9")
                j.cpu(16)
                j.env('PYSPARK_SUBMIT_ARGS', '--driver-memory 100g --executor-memory 100g pyspark-shell')
                j.memory('highmem')
                j.call(run_meta_analysis_ht,
                       pheno=pheno,
                       test_type=test_type,
                       result_type=result_type,
                       overwrite= True
                       )
                if args.test:
                    break

        b.run()

    if args.annotate_expected_pvalue:
        def create_expected_p_ht(pheno:str, pop:str, result_type:str, test_type:str, method:str, overwrite:bool):
            hl.init(
                master='local[32]',
                tmp_dir="gs://aou_tmp",
                default_reference="GRCh38",
            )
            def annotate_expected_pvalue(ht: hl.Table, method:str, p_field: str='Pvalue', k:int = 5000):
                n = ht.count()
                print(n)
                ht = ht.filter(hl.is_defined(ht[p_field]))
                keys = list(ht.key)
                if method == 'exact':
                    ht = ht.key_by()
                    ht = ht.order_by(ht[p_field]).add_index()
                    ht = ht.key_by(*keys)
                    ht = ht.annotate(**{f'{p_field}_expected': ht.idx / (n + 1)})
                    if n % 2 == 1:
                        median = [int(n * 0.5), int(n * 0.5) + 1]
                        median_p = ht.filter((ht.idx == median[0]) | (ht.idx == median[1]) )
                    else:
                        median_p = ht.filter(ht.idx == int(n * 0.5))

                    median_p = median_p.aggregate(hl.agg.mean(median_p[p_field]))
                    ht = ht.drop('idx')

                elif method == 'approx_cdf':
                    ht = ht.annotate_globals(ranks=ht.aggregate(hl.agg.approx_cdf(ht[p_field], k=k)))
                    ht = ht.annotate(p_rank_upper=hl.enumerate(ht.ranks['values']).find(lambda x: x[1] >= ht[p_field]))
                    ht = ht.annotate(p_rank=ht.p_rank_upper[0],
                                     upper_pvalue=ht.p_rank_upper[1])
                    ht = ht.annotate(lower_pvalue=ht.ranks['values'][ht.p_rank - 1],
                                     lower_rank=ht.ranks['ranks'][ht.p_rank - 1],
                                     upper_rank=ht.ranks['ranks'][ht.p_rank],
                                     )
                    ht = ht.annotate(
                        rank=hl.int64(hl.floor((ht.upper_rank - ht.lower_rank) / (ht.upper_pvalue - ht.lower_pvalue) * (
                                ht.Pvalue - ht.lower_pvalue) + ht.lower_rank)))
                    import time
                    ht = ht.checkpoint(f'gs://aou_tmp/{time.time()}.ht')
                    top_ht = ht.filter((ht.rank < n * 0.01) | (ht[p_field] < 1e-3))
                    rest_ht = ht.filter((ht.rank >= n * 0.01) & (ht[p_field] >= 1e-3))

                    top_ht = top_ht.key_by()
                    top_ht = top_ht.order_by(top_ht[p_field]).add_index('rank')
                    top_ht = top_ht.key_by(*keys)
                    ht = top_ht.union(rest_ht)
                    ht = ht.annotate(Pvalue_expected=ht.rank / (n + 1))
                    ht = ht.drop('p_rank_upper', 'p_rank', 'upper_pvalue', 'lower_pvalue', 'lower_rank', 'upper_rank')

                    from bisect import bisect
                    median = hl.eval(hl.floor(n * 0.5))
                    upper_idx = bisect(hl.eval(ht.ranks[1]), median)
                    upper_rank = hl.eval(ht.ranks[1])[upper_idx]
                    lower_rank = hl.eval(ht.ranks[1])[upper_idx - 1]
                    upper_pvalue = hl.eval(ht.ranks[0])[upper_idx]
                    lower_pvalue = hl.eval(ht.ranks[0])[upper_idx - 1]
                    median_p = (median - lower_rank) / (upper_rank - lower_rank) * (
                                upper_pvalue - lower_pvalue) + lower_pvalue
                    ht = ht.drop('ranks')

                else:
                    import sys
                    sys.exit("Method should be one of ['exact', 'approx_cdf']")


                lambda_gc = hl.eval(hl.qchisqtail(p=median_p, df=1, ncp=0, lower_tail=False) /
                                    hl.qchisqtail(p=0.5, df=1, ncp=0, lower_tail=True))
                print(f"Lambda GC: {lambda_gc}")

                ht = ht.annotate_globals(**{f'lambda_gc_{p_field}': lambda_gc})
                ht = ht.annotate(**{f'{p_field}': hl.if_else(ht[f'{p_field}'] == 0, 1e-320, ht[f'{p_field}'])}, )
                ht = ht.annotate(**{f'{p_field}_expected_log10': -hl.log10(ht[f'{p_field}_expected'])},)

                return ht

            path = f'gs://aou_results/250k/v1.1/{test_type}_results/{pop.upper()}/phenotype_{pheno}/{result_type}_results.ht'
            success = path.replace('.ht', f'_{method}_expected_p.ht/_SUCCESS')
            if (hl.hadoop_exists(path) and not hl.hadoop_exists(success)) or overwrite:
                overwrite=True
                ht = hl.read_table(path)
                ht.describe()
                if result_type == 'variant':
                    ht = annotate_expected_pvalue(ht=ht, method=method, p_field='Pvalue')
                    ht.describe()
                    ht = ht.checkpoint(path.replace('.ht', f'_{method}_expected_p.ht'), overwrite=overwrite,
                                  _read_if_exists=not overwrite)
                    ht = ht.filter(hl.case()
                                   .when(ht.Pvalue > 0.1, hl.rand_bool(0.001))
                                   .when(ht.Pvalue > 0.01, hl.rand_bool(0.01))
                                   .when(ht.Pvalue > 0.001, hl.rand_bool(0.1))
                                   .default(True))

                    ht.export(path.replace('.ht', f'_{method}_expected_p.txt.bgz'))
                else:
                    if pop == 'meta':
                        ht = ht.key_by('gene_id', 'gene_symbol', 'annotation')
                    ht = ht.annotate(max_MAF=hl.if_else(hl.is_missing(ht.max_MAF), -1, ht.max_MAF))
                    for max_MAF in [0.01, 0.001, 0.0001, -1]:
                        tag = 'Cauchy' if max_MAF == -1 else str(max_MAF)
                        sub_ht = ht.filter(ht.max_MAF == max_MAF)
                        sub_ht = annotate_expected_pvalue(ht=sub_ht, method=method, p_field='Pvalue')
                        sub_ht = annotate_expected_pvalue(ht=sub_ht, method=method, p_field='Pvalue_Burden')
                        sub_ht = annotate_expected_pvalue(ht=sub_ht, method=method, p_field='Pvalue_SKAT')
                        sub_ht.describe()
                        sub_ht = sub_ht.transmute_globals(**{f'lambda_gc_maxmaf_{tag}':
                                                          hl.struct(lambda_gc_Pvalue = sub_ht.lambda_gc_Pvalue,
                                                                   lambda_gc_Pvalue_Burden=sub_ht.lambda_gc_Pvalue_Burden,
                                                                   lambda_gc_Pvalue_SKAT=sub_ht.lambda_gc_Pvalue_SKAT,
                                                                   )})
                        sub_ht = sub_ht.checkpoint(path.replace('.ht', f'_{method}_expected_p_{tag}.ht'), overwrite=overwrite,
                                           _read_if_exists=not overwrite)
                        sub_ht.export(path.replace('.ht', f'_{method}_expected_p_{tag}.txt.bgz'))
                    ht = hl.read_table(path.replace('.ht', f'_{method}_expected_p_0.01.ht'))
                    for max_MAF in [0.001, 0.0001, -1]:
                        tag = 'Cauchy' if max_MAF == -1 else str(max_MAF)
                        tmp_ht = hl.read_table(path.replace('.ht', f'_{method}_expected_p_{tag}.ht'))
                        ht = ht.union(tmp_ht)
                        ht = ht.annotate_globals(
                            **{f'lambda_gc_maxmaf_{tag}': tmp_ht.index_globals()[f'lambda_gc_maxmaf_{tag}']})
                    ht = ht.key_by('gene_id', 'gene_symbol', 'annotation', 'max_MAF')
                    ht = ht.checkpoint(path.replace('.ht', f'_{method}_expected_p.ht'), overwrite=overwrite,
                                           _read_if_exists=not overwrite)




        backend = hb.ServiceBackend(
            billing_project="all-by-aou", remote_tmpdir='gs://aou_tmp/'
        )

        b = hb.Batch(
            name=f"aou_annotate_expected_p_{test_type}_{result_type}_result",
            backend=backend,
            default_storage="500Mi",
            default_cpu=8,
        )
        # phenos_to_run_by_pop  = {'afr': ['D03BA', 'CA_106.1', '628', 'PP_900', 'GU_626.11'],
        #                          'amr': ['695.4', 'BI_179.9', 'L01C', 'SS_810', '960'],
        #                          'eur': ['GI_552.2', 'CM_769.7', '536.8', 'MB_283.3', 'GE_970']}
        # pops = ['afr', 'amr', 'eur'] # Check some inflated phenotypes
        phenos_to_run_by_pop = read_pickle_dict(get_phenotype_info_path(version="pheno_by_pop_dict", extension="dict"))
        # phenos_to_run_by_pop_by_group = read_pickle_dict(get_phenotype_info_path(version="pheno_by_pop_by_group_dict", extension="dict"))
        # phenotypes = []
        # for pop in POPS:
        #     phenotypes = phenotypes + phenos_to_run_by_pop_by_group[pop]['lab_measurements'] + \
        #                  phenos_to_run_by_pop_by_group[pop]['processed_physical_measurement_table'] + \
        #                  [pheno for pheno in phenos_to_run_by_pop_by_group[pop]['random_phenotypes'] if
        #                   'continuous' in pheno]
        # phenos_to_run_by_pop['meta'] =  list(set(phenotypes))
        # phenos_to_run_by_pop_by_group['meta'] = {}
        # phenos_to_run_by_pop_by_group['meta']['lab_measurements'] =  list(set(phenotypes))
        for pop in pops:
            print(pop.upper())
            phenotypes = phenos_to_run_by_pop[pop]
            print(phenotypes)
            # phenotypes = phenos_to_run_by_pop_by_group[pop]['lab_measurements']
            for i in tqdm(range(len(phenotypes))):
                pheno = phenotypes[i]
                method = 'exact' if result_type == 'gene' else 'approx_cdf'
                # test_type = 'saige_variant' # TODO: change
                path = f'gs://aou_results/250k/v1.1/{test_type}_results/{pop.upper()}/phenotype_{pheno}/{result_type}_results.ht'
                success = path.replace('.ht', f'_{method}_expected_p.ht/_SUCCESS')
                if not hl.hadoop_exists(success) or args.overwrite:

                    j = b.new_python_job(name=f"{pop}_{pheno}_{test_type}_{result_type}_result")
                    j.image("hailgenetics/hail:0.2.128-py3.9")
                    j.memory('standard')
                    if result_type == 'variant' and test_type == 'saige':
                        j.storage('10GiB')
                    j.cpu(8)
                    j.env('PYSPARK_SUBMIT_ARGS', '--driver-memory 8g --executor-memory 8g pyspark-shell')
                    j.call(create_expected_p_ht,
                           pheno=pheno,
                           pop=pop,
                           test_type=test_type,
                           result_type=result_type,
                           method = 'exact' if result_type == 'gene' else 'approx_cdf',
                           overwrite=args.overwrite
                           )
                    if args.test:
                        break
                if args.test:
                    break



        b.run()

    if args.merge_call_stats_ht:
        ht1 = hl.read_table('gs://aou_results/250k/v1.1/utils/aou_all_exome_variant_info_pruned_250k.ht')
        ht2 = hl.read_table('gs://aou_results/250k/v1.1/utils/aou_all_genome_variant_info_pruned_250k.ht')

        ht = ht2.annotate(exome_freq = ht1[ht2.key].freq,
                          genome_freq = ht2.freq)
        ht = ht.filter(hl.is_defined(ht.exome_freq) | hl.is_defined(ht.genome_freq))

        ht = ht.select(**{f'genome_MAF_{pop}': ht.genome_freq[pop.upper()].AF[1] for pop in POPS},
                       **{f'exome_MAF_{pop}': ht.exome_freq[pop.upper()].AF[1] for pop in POPS},
                       **{f'genome_MAC_{pop}': ht.genome_freq[pop.upper()].AC[1] for pop in POPS},
                       **{f'exome_MAC_{pop}': ht.exome_freq[pop.upper()].AC[1] for pop in POPS},
                       )
        ht = ht.annotate(**{f'MAF_{pop}': hl.if_else(hl.is_missing(ht[f'genome_MAF_{pop}']), [ht[f'exome_MAF_{pop}']],
                                                      hl.if_else(hl.is_missing(ht[f'exome_MAF_{pop}']), [ht[f'genome_MAF_{pop}']],
                                                                 [ht[f'genome_MAF_{pop}'], ht[f'exome_MAF_{pop}']])) for pop in POPS},
                         **{f'MAC_{pop}': hl.if_else(hl.is_missing(ht[f'genome_MAC_{pop}']),
                                                      [ht[f'exome_MAC_{pop}']],
                                                      hl.if_else(hl.is_missing(ht[f'exome_MAC_{pop}']),
                                                                 [ht[f'genome_MAC_{pop}']],
                                                                 [ht[f'genome_MAC_{pop}'], ht[f'exome_MAC_{pop}']]))
                            for pop in POPS}
                         )

        ht.write('gs://aou_wlu/250k_ld/call_stats_merged.ht')

    if args.run_additional_load:
        def get_all_valid_variant_results_ht_paths(pop, test_type, result_type):
            results_dir = f'{RESULTS_PATH}/{test_type}_results/{pop.upper()}'
            all_phenos_dir = hl.hadoop_ls(results_dir)
            all_variant_outputs = get_files_in_parent_directory(all_phenos_dir, f'{result_type}_results.ht')
            return [x for x in all_variant_outputs]

        def generate_sumstats_mt(outputs_to_load,
                                 # heritability_dict, # TODO: add once ready
                                 temp_dir,
                                 pop,
                                 result_type,
                                 inner_mode='overwrite',
                                 checkpoint: bool = False):
            if result_type == 'variant':
                row_keys = ['locus', 'alleles']
            else:
                row_keys = ['gene_id', 'gene_symbol', 'annotation', 'max_MAF']
            col_keys = ['phenoname']

            if result_type == 'variant':
                all_hts = [unify_saige_ht_schema(hl.read_table(x), result_type=result_type, patch_case_control_count=x) for x in tqdm(outputs_to_load)]
            else:
                all_hts = [unify_saige_burden_ht_schema(hl.read_table(x)) for x in tqdm(outputs_to_load)]

            mt = join_pheno_hts_to_mt(all_hts=all_hts, row_keys=row_keys, col_keys=col_keys,
                                      temp_dir=temp_dir,
                                      inner_mode=inner_mode, repartition_final=20000)
            if checkpoint:
                # mt = mt.checkpoint(f'{temp_dir}_staging.mt', **{inner_mode: True})
                mt.write(f'{temp_dir}_staging.mt', **{inner_mode: True})
                mt = hl.read_matrix_table(f'{temp_dir}_staging.mt')
            if mt.inv_normalized.dtype == hl.tstr:
                mt = mt.annotate_cols(inv_normalized=hl.bool(mt.inv_normalized))
            pheno_ht = hl.read_table('gs://aou_results/250k/v1.1/utils/aou_phenotype_meta_info.ht')
            pheno_ht = pheno_ht.filter(pheno_ht.pop == pop).key_by('phenoname')
            mt = mt.annotate_cols(**pheno_ht[mt.col_key],
                                  inv_normalized = hl.str(mt.inv_normalized))

            # mt = check_and_annotate_with_dict(mt, heritability_dict, key)
            mt = mt.key_rows_by(*row_keys)
            return mt


        for pop in pops:
            POP_MT_PATH = get_aou_final_results_path(analysis_type=analysis_type, pop=pop, result_type=result_type, extension=f'mt')
            all_outputs = get_all_valid_variant_results_ht_paths(pop=pop, test_type=test_type, result_type=result_type)
            # heritability_dict = get_heritability_dict(pop) # TODO: add once heritability is ready

            pop_mt =hl.read_matrix_table(POP_MT_PATH)
            all_keys = pop_mt.phenoname.collect()
            loaded_phenos = [f'phenotype_{pheno}' for pheno in all_keys]
            lab_phenos = [f'phenotype_{pheno}' for pheno in LAB_CODES]

            # Some phenotypes were manually modified later but had already made it in the release
            # loaded_phenos |= {
            #     'categorical-41229-5N8', 'categorical-20004-1503', 'categorical-20004-1563',  # CSA
            #     'continuous-5097-irnt', 'continuous-104550-104550',  # EAS
            # }

            def _matches_any_pheno(pheno_path, phenos_to_match, result_type):
                return any(x for x in phenos_to_match if f'/{x}/{result_type}_results.ht' in pheno_path)



            if args.force_reload:
                pheno_matches = set(args.force_reload.split(','))
                outputs_to_load = [x for x in all_outputs if _matches_any_pheno(x, pheno_matches, result_type)]
            else:
                additional_outputs = [x for x in all_outputs if not _matches_any_pheno(x, loaded_phenos, result_type)]
                lab_outputs = [x for x in all_outputs if _matches_any_pheno(x, lab_phenos, result_type)]
                outputs_to_load = list(set(additional_outputs + lab_outputs))


                if args.load_only:
                    pheno_matches = set(args.load_only.split(','))
                    if '' in pheno_matches:
                        print('WARNING: Empty string in pheno_matches. Might reload more than expected')
                    outputs_to_load = [x for x in all_outputs if _matches_any_pheno(x, pheno_matches, result_type)]

            print(f'Loading {len(outputs_to_load)} additional HTs...')
            if len(outputs_to_load) < 20:
                print(outputs_to_load)
            if args.dry_run:
                continue
            if not len(outputs_to_load):
                continue

            mt = generate_sumstats_mt(outputs_to_load,
                                      #heritability_dict,
                                      f'{MY_BUCKET}/{pop}/{test_type}_{result_type}_{today}',
                                      pop=pop,
                                      result_type=result_type,
                                      inner_mode=inner_mode,
                                      checkpoint=True
                                      )
            mt.describe()
            if result_type == 'variant':
                mt = pull_out_fields_from_entries(mt, ['MarkerID'], 'rows',
                                                  agg_funcs=[lambda x: hl.agg.take(x, 1)[0]] * 1)
            else:
                mt = pull_out_fields_from_entries(mt=mt.annotate_entries(total_variants_pheno=mt.total_variants),
                                                  shared_fields=['MAC', 'Number_rare', 'Number_ultra_rare',
                                                                 'total_variants', 'interval', 'CHR', 'POS'],
                                                  index='rows',
                                                  agg_funcs=[lambda x: hl.agg.take(x, 1)[0]] * 7)
            if result_type == 'variant':
                mt = mt.select_entries('AC_Allele2', 'AF_Allele2', 'MissingRate', 'BETA', 'SE', 'var', 'p.value.NA', 'Is.SPA', 'AF_case', 'AF_ctrl', 'Pvalue')
            col_fields = ['n_cases', 'n_controls', 'heritability', 'saige_version', 'inv_normalized',  'pheno_sex', 'trait_type',
                                'category', 'pop', 'description', 'phecode_category', 'lambda_gc_exome', 'lambda_gc_acaf', 'lambda_gc_gene_burden_001']
            mt = mt.select_cols(*col_fields)

            original_mt = hl.read_matrix_table(POP_MT_PATH)
            pheno_ht = hl.read_table('gs://aou_results/250k/v1.1/utils/aou_phenotype_meta_info.ht')
            pheno_ht = pheno_ht.filter(pheno_ht.pop == pop).key_by('phenoname')
            original_mt = original_mt.annotate_cols(**pheno_ht[original_mt.col_key])
            original_mt = original_mt.select_cols(*col_fields)
            original_mt = original_mt.checkpoint(f'{TMP_BUCKET}/{pop}/{test_type}_{result_type}_before_{today}.mt', overwrite=True)
            if args.force_reload:
                original = original_mt.count_cols()
                original_mt = original_mt.filter_cols(
                    hl.literal(pheno_matches).contains(original_mt.phenoname), keep=False)
                print(f'\n\nGoing from {original} to {original_mt.count_cols()}...\n\n')


            mt = original_mt.union_cols(mt, row_join_type='outer')
            mt = mt.checkpoint(get_aou_final_results_path(analysis_type=analysis_type, pop=pop, result_type=result_type, extension=f'mt'),
                               overwrite=args.overwrite,
                               _read_if_exists=not args.overwrite)
            print(mt.count())





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
        "--pops", help="comma-separated list", default="afr,amr,eas,eur,mid,sas,meta"
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
    args = parser.parse_args()

    main(args)
