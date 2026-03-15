import hail as hl
import argparse
from datetime import date

from aou_gwas import *  # for dataproc
# from utils.utils import *  # for QoB
# from utils.resources import *  # for QoB
# from utils.results_loading import * # for QoB

FINAL_BUCKET = 'gs://allxall-finaldata-v2'

def reannotate_mt_dataproc(path, ancestry, result_type, analysis_type, sex=False):

    def result_name(result_type, analysis_type):
        if analysis_type =='gene':
            name = 'gene'
        elif result_type =='gene':
            name = 'exome'
        else:
            name = 'ACAF'
        return name

    name = result_name(result_type, analysis_type)
    lambda_name = {'gene': 'exome_gene_0.001', 'exome': 'exome_variant', 'ACAF': 'genome_variant'}
    row_ht_name = {'gene': 'exome_gene', 'exome': 'exome_variant', 'ACAF': 'ACAF_variant'}
    tag = '_sex_specific' if sex else ''

    mt = hl.read_matrix_table(path)
    print(mt.count())
    today = date.today().strftime("%y%m%d")
    print(name)
    print(f'gs://aou_tmp/mt_results/{name}_{analysis_type}{tag}_{ancestry}_before_{today}.mt')
    mt = mt.checkpoint(f'gs://aou_tmp/mt_results/{name}_{analysis_type}{tag}_{ancestry}_before_{today}.mt', _read_if_exists=True)
    pheno_ht = hl.read_table('gs://aou_analysis/v8/data/phenotype/summary/aou_v8_phenotype_meta_info_qced.ht')
    if not sex:
        pheno_ht = pheno_ht.filter((pheno_ht.ancestry == ancestry.upper()) & (pheno_ht.pheno_sex == 'both')).key_by('phenoname')
        lambda_ht = hl.read_table(f'gs://aou_analysis/v8/data/qc/aou_{ancestry}_{lambda_name[name]}_phenotype_lambda_gc_filtered.ht')
        lambda_ht = lambda_ht.select(lambda_gc_raw = lambda_ht.lambda_gc_raw, lambda_gc_hq = lambda_ht.lambda_gc_filtered)
        pheno_ht = pheno_ht.drop('ancestry', 'pheno_sex')
    else:
        pheno_ht = pheno_ht.filter((pheno_ht.ancestry == ancestry.upper())).key_by('phenoname')
        pheno_ht = pheno_ht.drop('ancestry')

    # Prepare row information
    if name != 'exome':
        row_info_path = f'gs://aou_analysis/v8/data/qc/aou_{row_ht_name[name]}_qc.ht'
    else:
        row_info_path = f'gs://aou_analysis/v8/data/qc/aou_{row_ht_name[name]}_qc_annotated.ht'
    print(row_info_path)
    row_info_ht = hl.read_table(row_info_path)
    if analysis_type == 'variant':
        tag = ancestry if ancestry != 'meta' else 'all'
        row_info_ht = row_info_ht.annotate(
            AF=row_info_ht.freq[f'{tag.upper()}'].AF,
            AC=row_info_ht.freq[f'{tag.upper()}'].AC,
            AN=row_info_ht.freq[f'{tag.upper()}'].AN,
            homozygote_count=row_info_ht.freq[f'{tag.upper()}'].homozygote_count,
        ).drop('freq', 'freq_raw', 'old_locus', 'old_alleles', 'was_split', 'a_index', 'as_vets', 'filters')
    else:
        tag = ancestry.upper() if ancestry != 'meta' else 'global'
        ht = row_info_ht.select(mean_coverage_raw = row_info_ht.mean_coverage,
                                mean_coverage_max_MAF = row_info_ht['mean_coverage_0.01'],
                                CAF_raw = row_info_ht[f'CAF_{tag}_raw'],
                                CAF_max_MAF = row_info_ht[f'CAF_{tag}_0.01'])
        ht = ht.annotate(max_MAF = 0.01)
        for maxmaf in [0.001, 0.0001]:
            tmp_ht = row_info_ht.select(mean_coverage_raw=row_info_ht.mean_coverage,
                                        mean_coverage_max_MAF=row_info_ht[f'mean_coverage_{maxmaf}'],
                                        CAF_raw=row_info_ht[f'CAF_{tag}_raw'],
                                        CAF_max_MAF=row_info_ht[f'CAF_{tag}_{maxmaf}'])
            tmp_ht = tmp_ht.annotate(max_MAF=maxmaf)
            ht = ht.union(tmp_ht)
        row_info_ht = ht.key_by('gene_id', 'gene_symbol', 'annotation', 'max_MAF')

    # Annotate variant
    mt = mt.annotate_rows(**row_info_ht[mt.row_key])
    # Annotate phenotype
    mt = mt.annotate_cols(**pheno_ht[mt.col_key])
    mt = mt.annotate_cols(**lambda_ht[mt.col_key])
    phenotype_outliers = ['PP_928.1', 'PP_932.11']
    phenotype_to_flag = ['3011099', '3046664', '3013861', '3052648', '3045792_3016921']
    mt = mt.filter_cols(
        (~hl.literal(phenotype_outliers).contains(mt.phenoname))
    ) 
    mt = mt.annotate_cols(hq_pheno=~hl.literal(phenotype_to_flag).contains(mt.phenoname))

    # QC flags
    if analysis_type == 'variant':
        mt = mt.annotate_entries(hq_exp_AC = mt.n_cases * mt.AF[1] >= 5)
        mt = mt.annotate_rows(hq_AF_variant = mt.AF[1] > 0.0001,
                              hq_exp_AC_variant = hl.agg.count_where(mt.hq_exp_AC) > 0,)
        mt = mt.transmute_rows(quality_flags=hl.struct(hq_exp_AC_variant=mt.hq_exp_AC_variant),
                               quality_flags_lambda=hl.struct(hq_AF_variant=mt.hq_AF_variant),
                               hq_variant = mt.hq_exp_AC_variant,
                               hq_variant_lambda = mt.hq_exp_AC_variant & mt.hq_AF_variant)
        mt = mt.drop('MarkerID')
    else:
        min_caf = 1e-5 if ancestry not in ['mid', 'sas'] else 1e-4
        print(f'MIN CAF: {min_caf}')
        if ancestry == 'meta':
            mt = mt.annotate_entries(
                total_variants = hl.max(mt.summary_stats.total_variants_pheno)
            )
        else:
            mt = mt.annotate_entries(
                total_variants = mt.total_variants_pheno
            )
        mt = mt.annotate_entries(hq_exp_CAC = mt.n_cases * mt.CAF_max_MAF >= 5,
                                 hq_n_var_5 = mt.total_variants >= 5,
                                 hq_n_var_10 = mt.total_variants > 10)
        mt = mt.annotate_rows(hq_coverage = mt.mean_coverage_max_MAF > 10,
                              hq_CAF = mt.CAF_max_MAF > min_caf,
                              hq_exp_CAC_gene = hl.agg.count_where(mt.hq_exp_CAC) > 0,
                              hq_n_var_5_gene = hl.agg.count_where(mt.hq_n_var_5) > 0,
                              hq_n_var_10_gene = hl.agg.count_where(mt.hq_n_var_10) > 0,)
        mt = mt.transmute_rows(
            quality_flags=hl.struct(hq_coverage=mt.hq_coverage, hq_n_var_5 = mt.hq_n_var_5_gene, hq_exp_CAC_gene=mt.hq_exp_CAC_gene),
            quality_flags_lambda=hl.struct(hq_n_var_10=mt.hq_n_var_10_gene, hq_CAF = mt.hq_CAF),
            hq_gene=mt.hq_coverage & mt.hq_n_var_5_gene & mt.hq_exp_CAC_gene,
            hq_gene_lambda=mt.hq_coverage & mt.hq_n_var_10_gene & mt.hq_exp_CAC_gene & mt.hq_CAF
        )
        mt = mt.annotate_rows(
            hq_gene = hl.if_else(mt.annotation == 'Cauchy', hl.missing(hl.tbool), mt.hq_gene),
            hq_gene_lambda = hl.if_else(mt.annotation == 'Cauchy', hl.missing(hl.tbool), mt.hq_gene_lambda)
        )
        if analysis_type == 'gene' and ancestry != 'meta':
            mt = mt.drop('heritability_1', 'saige_version_1', 'inv_normalized_1', 'n_cases_1', 'n_controls_1', 'N_1')

    if ancestry == 'meta':
        inv_normalized = True
        mt = mt.annotate_cols(pheno_data = mt.pheno_data.map(lambda x: x.drop('inv_normalized', 'N')))
        if result_type == 'gene' and analysis_type == 'gene':
            mt = mt.drop('Pvalue', 'Pvalue_Burden', 'Pvalue_SKAT', 'BETA_SKATO', 'BETA_Burden', 'BETA_SKAT', 'weighted_Z_numerator_SKATO', 'weighted_Z_numerator_Burden', 'weighted_Z_numerator_SKAT', 'summary_stats')
            mt = mt.annotate_cols(ancestries = mt.pheno_data.ancestry)
        else:
            mt = mt.transmute_entries(**mt.meta_analysis[0])
            mt = mt.transmute_cols(ancestries=mt.meta_analysis_data[0].ancestry)
    else:
        inv_normalized = list(mt.aggregate_cols(hl.agg.collect_as_set(mt.inv_normalized)))[0]
        mt = mt.drop('inv_normalized', 'N')
    mt = mt.annotate_globals(ancestry=ancestry, inv_normalized=inv_normalized)
    
    mt.describe()
    mt = mt.checkpoint(path, overwrite=True)
    print(mt.count())

    mt.cols().show()

def main(args):
    hl.init(
        tmp_dir='gs://aou_tmp/',
        driver_memory="highmem",
        driver_cores=8,
        worker_memory="highmem",
        worker_cores=1,
        default_reference="GRCh38",
        app_name='Final production'
    )
    if args.reannotate_mt:
        if args.load_gene_results and args.single_variant_only:
            table_name = 'exome_variant'
        elif args.load_gene_results and not args.single_variant_only:
            table_name = 'gene'
        elif not args.load_gene_results and args.single_variant_only:
            table_name = 'genome_variant'
        for ancestry in ['meta', 'afr', 'amr', 'eas', 'eur', 'sas', 'mid']:
            path = f'gs://aou_analysis/v8/mt_results/{ancestry.upper()}_{table_name}.mt'
            print(path)
            analysis_type = "variant" if args.single_variant_only else "gene"
            result_type = "gene" if args.load_gene_results else "variant"
            reannotate_mt_dataproc(path=path, ancestry=ancestry, result_type=result_type, analysis_type=analysis_type)


    if args.minor_edits:
        for ancestry in ['afr', 'amr', 'eas', 'eur', 'sas', 'mid', 'meta']:
            path = f'gs://aou_analysis/v8/mt_results/{ancestry.upper()}_gene.mt'
            mt = hl.read_matrix_table(path)
            print(path)
            print(mt.count())
            today = date.today().strftime("%y%m%d")
            tmp_path = f'gs://aou_tmp_30_days/mt_results/gene_gene_{ancestry}_before_{today}.mt'
            mt = mt.checkpoint(tmp_path, overwrite=True)
            mt = mt.annotate_rows(
                hq_gene = hl.if_else(mt.annotation == 'Cauchy', hl.missing(hl.tbool), mt.hq_gene),
                hq_gene_lambda = hl.if_else(mt.annotation == 'Cauchy', hl.missing(hl.tbool), mt.hq_gene_lambda)
            )
            mt = mt.checkpoint(path, overwrite=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--single-variant-only", help="Run single variant SAIGE", action="store_true"
    )
    parser.add_argument(
        "--reannotate-mt", help="Finalize MTs", action="store_true"
    )
    parser.add_argument(
        "--load-gene-results", help="Load gene level test results", action="store_true"
    )
    parser.add_argument(
        "--minor-edits", help="Run minor edits", action="store_true"
    )
    args = parser.parse_args()

    main(args)
