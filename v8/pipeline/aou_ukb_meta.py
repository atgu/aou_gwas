import argparse
import hail as hl
import hailtop.batch as hb
import pickle, os

ANALYSIS_BUCKET = 'gs://aou_analysis/v8'
TMP_BUCKET = 'gs://aou_tmp/biobank_meta'
UKB_GENE_MT = 'gs://ukbb-exome-public/500k/qc/gene_qc_metrics_ukb_exomes_500k.mt'
UKB_EXOME_MT = 'gs://ukbb-exome-public/500k/results/variant_results.mt'
UKB_QC_GENOME_MT = 'gs://ukb-diverse-pops-public/sumstats_release/meta_analysis.h2_qc.mt'
UKB_RAW_GENOME_MT = 'gs://ukb-diverse-pops-public/sumstats_release/meta_analysis.mt'
UKB_GENOME_MT = 'gs://ukb-diverse-pops-public/sumstats_release/meta_analysis.mt'
UKB_GENE_PHENO_HT = 'gs://ukbb-exome-public/500k/results/pheno_results.ht'
AOU_EXOME_MT = 'gs://aou_results/414k/mt_results/META_exome_variant.mt'
AOU_GENOME_MT = 'gs://aou_results/414k/mt_results/META_genome_variant.mt'
AOU_GENE_MT = 'gs://aou_results/414k/mt_results/META_gene.mt'
AOU_EUR_GENE_MT = 'gs://aou_results/414k/mt_results/EUR_gene.mt'
GENEBASS_PHENO_MAP_TSV = 'gs://aou_wlu/v8_analysis/aou_ukb_meta/aou_ukb_matched_all_phenotype_v8.csv'
PANUKB_PHENO_MAP_TSV = 'gs://aou_wlu/v8_analysis/aou_ukb_meta/aou_pan_ukb_matched_all_phenotype_v8.csv'

def format_aou_genome_variant_mt(mt):
    v8_pheno_map = hl.import_table(PANUKB_PHENO_MAP_TSV, delimiter='\t')
    ht = v8_pheno_map.key_by('phenoname')
    mt = mt.drop('N')
    mt = mt.select_globals()
    mt = mt.annotate_cols(**ht[mt.col_key],
                          N = mt.n_cases + mt.n_controls)
    mt = mt.filter_cols(hl.is_defined(mt.ukb_phenocode))
    mt = mt.select_cols('trait_type', 'category', 'description', 'N', 'ukb_phenocode')
    mt = mt.select_rows()
    mt = mt.select_entries(
        BETA = mt.BETA,
        SE = mt.SE,
        Pvalue = hl.exp(mt.Pvalue),
        AF_cases = mt.AF_Cases,
        AF_controls = mt.AF_Controls,
        AF_Allele2 = mt.AF_Allele2
    )
    return mt

def add_liftover_rg37_to_rg38_mt(mt: hl.MatrixTable):
    """
    Add liftover to Hail Table from rg37 to rg38
    :param Table ht: Hail Table to add liftover on
    :return: Hail MatrixTable
    :rtype: hl.MatrixTable
    """
    rg37 = hl.get_reference("GRCh37")
    rg38 = hl.get_reference("GRCh38")
    if not rg37.has_liftover("GRCh38"):
        rg37.add_liftover(
            "gs://hail-common/references/grch37_to_grch38.over.chain.gz", rg38
        )
    mt = mt.annotate_rows(new_locus=hl.liftover(mt.locus, "GRCh38"))
    mt = mt.filter_rows(hl.is_defined(mt.new_locus))
    mt = mt.checkpoint('gs://aou_tmp/v8/pan_ukb_grch38_lifted_no_h2_qc.mt', overwrite=True)
    mt = mt.key_rows_by(locus=mt.new_locus, alleles=mt.alleles)
    mt = mt.drop('new_locus')
    return mt

def format_ukb_genome_variant_mt(mt):
    v8_pheno_map = hl.import_table(PANUKB_PHENO_MAP_TSV, delimiter='\t')
    ht = v8_pheno_map.key_by('ukb_phenocode')
    mt = add_liftover_rg37_to_rg38_mt(mt)
    mt = mt.checkpoint('gs://aou_tmp/v8/pan_ukb_grch38_no_h2_qc.mt', overwrite=True)
    mt = mt.key_cols_by(ukb_phenocode = mt.phenocode)
    mt = mt.annotate_cols(**ht[mt.col_key],
                          N = mt.meta_analysis_data[0].n_cases + mt.meta_analysis_data[0].n_controls) 
    mt = mt.filter_cols(hl.is_defined(mt.phenoname))
    mt = mt.key_cols_by('phenoname').select_cols('trait_type', 'category', 'description', 'N', 'ukb_phenocode')
    mt = mt.select_rows()
    mt = mt.select_entries(
        BETA = mt.meta_analysis[0].BETA,
        SE = mt.meta_analysis[0].SE,
        Pvalue = 10**(-mt.meta_analysis[0].Pvalue),
        AF_cases = mt.meta_analysis[0]['AF_Cases'],
        AF_controls = mt.meta_analysis[0]['AF_Controls'],
        AF_Allele2 = mt.meta_analysis[0].AF_Allele2
    )
    return mt

def format_aou_exome_variant_mt(mt):
    v8_pheno_map = hl.import_table(GENEBASS_PHENO_MAP_TSV, delimiter='\t')
    ht = v8_pheno_map.key_by('phenoname')
    mt = mt.drop('N')
    mt = mt.select_globals()
    mt = mt.annotate_cols(**ht[mt.col_key],
                          N = mt.n_cases + mt.n_controls)
    mt = mt.filter_cols(hl.is_defined(mt.ukb_phenocode))
    mt = mt.annotate_rows(
        call_stats = hl.struct(
            AC = mt.AC[1],
            AF = mt.AF[1],
            AN = mt.AN,
            homozygote_count = mt.homozygote_count[1]
        )
    )
    mt = mt.select_cols('trait_type', 'category', 'description', 'N', 'ukb_phenocode')
    mt = mt.select_rows('call_stats', 'gene_id', gene = mt.gene_symbol, annotation = mt.annotation)
    mt = mt.select_entries(
        BETA = mt.BETA,
        SE = mt.SE,
        Pvalue = hl.exp(mt.Pvalue),
        AF_cases = mt.AF_Cases,
        AF_controls = mt.AF_Controls,
        AF_Allele2 = mt.AF_Allele2
    )
    return mt

def format_ukb_exome_variant_mt(mt):
    v8_pheno_map = hl.import_table(GENEBASS_PHENO_MAP_TSV, delimiter='\t')
    ht = v8_pheno_map.key_by('ukb_phenocode')
    mt = mt.key_cols_by(ukb_phenocode = mt.phenocode)
    mt = mt.annotate_cols(**ht[mt.col_key],
                          N = mt.n_cases + mt.n_controls) 
    mt = mt.filter_cols(hl.is_defined(mt.phenoname))
    mt = mt.key_cols_by('phenoname').select_cols('trait_type', 'category', 'description', 'N', 'ukb_phenocode')
    mt = mt.select_rows('call_stats', gene_id = '' , gene = mt.gene, annotation = mt.annotation)
    mt = mt.select_entries(
        BETA = mt.BETA,
        SE = mt.SE,
        Pvalue = mt.Pvalue,
        AF_cases = mt['AF.Cases'],
        AF_controls = mt['AF.Controls'],
        AF_Allele2 = mt.AF
    )
    return mt

def run_meta_analysis(mt: hl.MatrixTable,
                      beta_field: str = 'BETA',
                      se_field: str = 'SE',
                      af_allele2_field: str = 'AF_Allele2',
                      af_case_field: str='AF_cases',
                      af_ctrl_field: str='AF_controls'):
    """
    Run inverse-variance fixed-effect meta-analysis for a given MatrixTable.

    :param mt: Input MatrixTable, formatted similar to `load_final_sumstats_mt()`
    ...
    :return: Result MatrixTable with `meta_analysis` entries and `meta_analysis_data` columns.
    :rtype: MatrixTable
    """

    # Run fixed-effect meta-analysis (all + leave-one-out)
    mt = mt.annotate_entries(
        unnorm_beta=mt[beta_field] / (mt[se_field] ** 2), inv_se2=1 / (mt[se_field] ** 2)
    )
    mt = mt.annotate_entries(
        sum_unnorm_beta=hl.sum(mt.unnorm_beta),
        sum_inv_se2=hl.sum(mt.inv_se2),
    )
    mt = mt.transmute_entries(
        META_BETA=mt.sum_unnorm_beta / mt.sum_inv_se2,
        META_SE = hl.sqrt(1 /mt.sum_inv_se2)
    )
    mt = mt.annotate_entries(
        META_Pvalue= hl.log(2) + hl.pnorm(-hl.abs(mt.META_BETA / mt.META_SE), log_p=True)
    )

    # Run heterogeneity test (Cochran's Q)
    mt = mt.annotate_entries(
        META_Q=hl.sum((mt[beta_field] - mt.META_BETA) ** 2 * mt.inv_se2),
        variant_exists=hl.is_defined(mt[beta_field]),
    )

    mt = mt.annotate_entries(
        META_Pvalue_het=hl.pchisqtail(mt.META_Q, 1, log_p=True)
    )

    # Add other annotations
    if af_case_field is not None or af_ctrl_field is not None:
        mt = mt.annotate_entries(
            ac_cases=mt[af_case_field] * mt.N,
            ac_controls=mt[af_ctrl_field] * mt.N,
            META_AC_Allele2=hl.sum(mt[af_allele2_field] * mt.N),
            META_N=hl.sum(mt.N),
        )


    mt = mt.annotate_entries(
        META_AF_Allele2=mt.META_AC_Allele2 / mt.META_N,
        META_AF_Cases=hl.sum(mt.ac_cases) / mt.META_N,
        META_AF_Controls=hl.sum(mt.ac_controls) / mt.META_N,
    )
    mt = mt.drop(
        "unnorm_beta", "inv_se2", "variant_exists", "ac_cases", "ac_controls", "META_AC_Allele2"
    )

    return mt

def format_aou_gene_meta_mt(mt, maxmaf=0.01):
    v8_pheno_map = hl.import_table(GENEBASS_PHENO_MAP_TSV, delimiter='\t')
    # Pheno-level (column)
    ht = v8_pheno_map.key_by('phenoname')
    mt = mt.annotate_cols(**ht[mt.col_key])
    mt = mt.filter_cols(hl.is_defined(mt.ukb_phenocode))
    mt = mt.annotate_cols(n_cases = hl.sum(mt.pheno_data.n_cases),
                          n_controls = hl.sum(mt.pheno_data.n_controls))
    mt = mt.annotate_cols(
        Neff = hl.if_else(mt.n_controls ==0, hl.float64(mt.n_cases), (4.0 * mt.n_cases * mt.n_controls) / (mt.n_cases + mt.n_controls)),
        N = (mt.n_cases + mt.n_controls)
    )
    mt = mt.select_cols('ukb_phenocode', 'Neff', 'description', 'N', 'n_cases', 'n_controls')
    # Gene-level (row)
    mt = mt.filter_rows(mt.max_MAF == maxmaf)
    mt = mt.key_rows_by('gene_id', 'gene_symbol', 'annotation')
    mt = mt.select_rows()
    # Entry level 
    mt = mt.select_entries(
        MAC = mt.META_MAC,
        Pvalue = mt.META_Pvalue_SKATO,
        Pvalue_Burden = mt.META_Pvalue_Burden,
        Pvalue_SKAT = mt.META_Pvalue_SKAT,
        BETA_SKATO = mt.META_Stats_SKATO,
        BETA_Burden = mt.META_Stats_Burden,
        BETA_SKAT = mt.META_Stats_SKAT,
    )
    mt = mt.annotate_entries(CAF = mt.MAC / (2*mt.N))
    mt = mt.annotate_entries(two_pq = 2*mt.CAF*(1-mt.CAF))
    mt = mt.annotate_entries(Pvalue_defined = hl.int64(hl.is_defined(mt.Pvalue)),
                             Pvalue_Burden_defined = hl.int64(hl.is_defined(mt.Pvalue_Burden)),
                             Pvalue_SKAT_defined = hl.int64(hl.is_defined(mt.Pvalue_SKAT)))
    mt = mt.annotate_entries(two_pq_Pvalue = mt.Pvalue_defined * mt.two_pq,
                             two_pq_Pvalue_Burden = mt.Pvalue_Burden_defined * mt.two_pq,
                             two_pq_Pvalue_SKAT = mt.Pvalue_SKAT_defined * mt.two_pq)
    return mt

def format_aou_gene_eur_mt(mt, maxmaf=0.01):
    # Pheno-level (column)
    v8_pheno_map = hl.import_table(GENEBASS_PHENO_MAP_TSV, delimiter='\t')
    ht = v8_pheno_map.key_by('phenoname')
    mt = mt.annotate_cols(**ht[mt.col_key])
    mt = mt.filter_cols(hl.is_defined(mt.ukb_phenocode))
    mt = mt.annotate_cols(
        Neff = hl.if_else(mt.n_controls ==0, hl.float64(mt.n_cases), (4.0 * mt.n_cases * mt.n_controls) / (mt.n_cases + mt.n_controls)),
        N = (mt.n_cases + mt.n_controls)
    )
    mt = mt.select_cols('ukb_phenocode', 'Neff', 'description', 'N', 'n_cases', 'n_controls')
    # Gene-level (row)
    mt = mt.filter_rows(mt.max_MAF == maxmaf)
    mt = mt.key_rows_by('gene_id', 'gene_symbol', 'annotation')
    mt = mt.select_rows()
    # Entry level 
    mt = mt.select_entries(
        Pvalue_Burden = mt.Pvalue_Burden,
        BETA_Burden = mt.BETA_Burden,
        SE_Burden = mt.SE_Burden
    )
    return mt

def format_ukb_gene_mt(ukb_mt):
    v8_pheno_map = hl.import_table(GENEBASS_PHENO_MAP_TSV, delimiter='\t')
    ht = v8_pheno_map.key_by('ukb_phenocode')
    ukb_mt = ukb_mt.key_cols_by(ukb_phenocode = ukb_mt.phenocode)
    ukb_mt = ukb_mt.annotate_cols(n_controls = hl.if_else(hl.is_missing(ukb_mt.n_controls), 0, ukb_mt.n_controls))
    ukb_mt = ukb_mt.annotate_cols(
        phenoname = ht[ukb_mt.col_key].phenoname, 
        Neff = hl.if_else(ukb_mt.n_controls ==0, hl.float64(ukb_mt.n_cases), (4.0 * ukb_mt.n_cases * ukb_mt.n_controls) / (ukb_mt.n_cases + ukb_mt.n_controls)),
        N = ukb_mt.n_cases + ukb_mt.n_controls
    )
    ukb_mt = ukb_mt.filter_cols(hl.is_defined(ukb_mt.phenoname))                            
    ukb_mt = ukb_mt.key_cols_by('phenoname').select_cols('ukb_phenocode', 'Neff', 'description', 'N', 'n_cases', 'n_controls')
    ukb_mt = ukb_mt.annotate_entries(
        Pvalue_Burden = ukb_mt.Pvalue_Burden,
        BETA_SKATO = hl.float64(hl.is_defined(ukb_mt.Pvalue)),
        BETA_SKAT = hl.float64(hl.is_defined(ukb_mt.Pvalue_SKAT)),
        MAC = hl.int64(ukb_mt.total_variants_pheno)
    )
    ukb_mt = ukb_mt.select_entries('MAC', 'Pvalue', 'Pvalue_Burden', 'Pvalue_SKAT', 'BETA_SKATO', 'BETA_Burden', 'BETA_SKAT', 'SE_Burden')
    ukb_mt = ukb_mt.select_rows().key_rows_by()
    ukb_mt = ukb_mt.annotate_rows(
        annotation = hl.if_else(ukb_mt.annotation == 'missense|LC', 'missenseLC', 
                                hl.if_else(ukb_mt.annotation == 'pLoF|missense|LC', 'pLoF;missenseLC', ukb_mt.annotation))
    ).key_rows_by('gene_id', 'gene_symbol', 'annotation')
    ukb_mt = ukb_mt.annotate_entries(CAF = ukb_mt.MAC / (2*ukb_mt.N))
    ukb_mt = ukb_mt.annotate_entries(two_pq = 2*ukb_mt.CAF*(1-ukb_mt.CAF))
    ukb_mt = ukb_mt.annotate_entries(Pvalue_defined = hl.int64(hl.is_defined(ukb_mt.Pvalue)),
                             Pvalue_Burden_defined = hl.int64(hl.is_defined(ukb_mt.Pvalue_Burden)),
                             Pvalue_SKAT_defined = hl.int64(hl.is_defined(ukb_mt.Pvalue_SKAT)))
    ukb_mt = ukb_mt.annotate_entries(two_pq_Pvalue = ukb_mt.Pvalue_defined * ukb_mt.two_pq,
                             two_pq_Pvalue_Burden = ukb_mt.Pvalue_Burden_defined * ukb_mt.two_pq,
                             two_pq_Pvalue_SKAT = ukb_mt.Pvalue_SKAT_defined * ukb_mt.two_pq)
    return ukb_mt

def run_stouffer_mt(mt):
    P_FIELDS = ['Pvalue', 'Pvalue_Burden', 'Pvalue_SKAT']
    P_TESTS = {'SKATO': 'Pvalue', 'Burden': 'Pvalue_Burden', 'SKAT':'Pvalue_SKAT'}

    for test in list(P_TESTS.keys()):
        print(f'Meta analyzing {test} results...')
        two_tail = (test == 'Burden')
        beta_field = f'BETA_{test}'
        p_field = P_TESTS[test]
        mt.describe()
        if two_tail:
            mt = mt.annotate_entries(**{p_field: mt[p_field]/2})

        mt = mt.annotate_entries(**{f'weighted_Z_numerator_{test}' :
                                    hl.map(lambda w, x, y, z: hl.sqrt(w)*hl.sqrt(x)*(-hl.qnorm(y))*  hl.sign(z),
                                           mt[f'two_pq_{p_field}'], mt.Neff, mt[p_field], mt[beta_field])})
        mt = mt.annotate_entries(**{f'META_Stats_{test}': 
                                    hl.sum(mt[f'weighted_Z_numerator_{test}']) / (hl.sqrt(hl.sum(mt[f'two_pq_{p_field}']*mt.Neff)))}, )

        if two_tail:
            mt = mt.annotate_entries(**{f'META_Pvalue_{test}': 2*hl.pnorm(hl.abs(mt[f'META_Stats_{test}']), lower_tail=False)})
        else:
            mt = mt.annotate_entries(**{f'META_Pvalue_{test}': hl.pnorm(mt[f'META_Stats_{test}'], lower_tail=False)})

    return mt

def format_aou_gene_meta_ht(ht, ukb_code, maxmaf=0.01):
    ht = ht.filter(ht.max_MAF == maxmaf)
    ht = ht.key_by('gene_id', 'gene_symbol', 'annotation')
    phenoname = hl.eval(ht.phenoname)
    n_cases = hl.int64(hl.eval(ht.N_cases))
    n_controls = hl.int64(hl.eval(ht.N_controls))
    ht = ht.select_globals()
    ht = ht.select(
        phenoname = phenoname,
        ukb_phenocode = ukb_code, 
        n_cases = n_cases,
        n_controls = n_controls,
        Pvalue = ht.Pvalue,
        Pvalue_Burden = ht.Pvalue_Burden,
        Pvalue_SKAT = ht.Pvalue_SKAT,
        BETA_SKATO = ht.META_Stats_SKATO,
        BETA_Burden = ht.META_Stats_Burden,
        BETA_SKAT = ht.META_Stats_SKAT,
        interval = ht.interval[0]
    )
    return ht

def format_ukb_gene_ht(ht, aou_code):
    ht = ht.key_by()
    if not aou_code.startswith('random'):
        n_cases = hl.int64(hl.eval(ht.n_cases))
        n_controls = hl.int64(hl.eval(ht.n_controls))
        ht = ht.select_globals()
        ht = ht.annotate(n_cases = n_cases, n_controls = n_controls)                             
    ht = ht.annotate(
        phenoname = aou_code,
        ukb_phenocode = ht.phenocode, 
        n_cases = hl.int64(ht.n_cases),
        n_controls = hl.int64(ht.n_controls),
        # Pvalue_Burden = ht.Pvalue_Burden/2, # double check
        Pvalue_Burden = ht.Pvalue_Burden,
        BETA_SKATO = hl.float64(hl.is_defined(ht.Pvalue)),
        BETA_SKAT = hl.float64(hl.is_defined(ht.Pvalue_SKAT)),
                annotation = hl.if_else(ht.annotation == 'missense|LC', 'missenseLC', 
                                hl.if_else(ht.annotation == 'pLoF|missense|LC', 'pLoF;missenseLC', ht.annotation))
    )
    ht = ht.key_by('gene_id', 'gene_symbol', 'annotation')
    ht = ht.select('phenoname', 'ukb_phenocode', 'n_cases', 'n_controls', 'Pvalue', 'Pvalue_Burden', 'Pvalue_SKAT', 'BETA_SKATO', 'BETA_Burden', 'BETA_SKAT', 'interval')
    return ht

def run_single_meta_analysis_ht(aou_ht_path:str, ukb_ht_path:str, aou_pheno:str, ukb_pheno:str, table_name:str, analysis_type:str, overwrite:bool):
    # hl.init(
    #     gcs_requester_pays_configuration='ukbb-exome-pharma'
    # )
    aou_ht = hl.read_table(aou_ht_path).key_by('locus', 'alleles')
    ukb_ht = hl.read_table(ukb_ht_path).key_by('locus', 'alleles')
    if analysis_type == 'gene':
        aou_ht = format_aou_gene_meta_ht(aou_ht, ukb_pheno)
        ukb_ht = format_ukb_gene_ht(ukb_ht, aou_pheno)
    else:
        N = hl.eval(aou_ht.N)
        N_case = hl.eval(aou_ht.N_cases)
        N_ctrl = hl.eval(aou_ht.N_controls)
        aou_ht = aou_ht.select_globals()
        aou_ht = aou_ht.annotate(N = N, N_case = N_case, N_ctrl = N_ctrl)
        aou_ht = aou_ht.select('BETA', 'SE', 'Pvalue', 'AC_Allele2', 'N', 'N_case', 'N_ctrl', 'AC_case', 'AC_ctrl')
        ukb_ht = ukb_ht.select_globals()
        ukb_ht = ukb_ht.annotate(N = ukb_ht.n_cases + ukb_ht.n_controls, 
        AC_Allele2 = hl.int64(ukb_ht.AC),
        N_case = ukb_ht.n_cases, N_ctrl = ukb_ht.n_controls)
        ukb_ht = ukb_ht.annotate(
        AC_case = ukb_ht['AF.Cases'] * ukb_ht.N_case * 2,
        AC_ctrl = ukb_ht['AF.Controls'] * ukb_ht.N_ctrl * 2)
        ukb_ht = ukb_ht.select('BETA', 'SE', 'Pvalue', 'AC_Allele2', 'N', 'N_case', 'N_ctrl', 'AC_case', 'AC_ctrl')
    ht = aou_ht.union(ukb_ht)
    path = f'{ANALYSIS_BUCKET}/ht_results/aou_ukb_META/phenotype_{aou_pheno}_{ukb_pheno}_{table_name}_raw.ht'
    overwrite = not hl.hadoop_exists(f'{path}/_SUCCESS')
    ht = ht.checkpoint(path, overwrite=overwrite, _read_if_exists=not overwrite)
    if analysis_type == 'gene':
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
            N = hl.if_else(ht.n_controls ==0, hl.float64(ht.n_cases), (4.0 * ht.n_cases * ht.n_controls) / (ht.n_cases + ht.n_controls))
        ) 
        ht.drop('interval').export(path.replace('.ht', '.txt.bgz'))  
    if analysis_type == 'variant':
        ht = ht.annotate(BETA = hl.float64(ht.BETA),
                         SE = hl.float64(ht.SE))
        ht = ht.annotate(weight = 1 / (ht.SE ** 2))
        meta_ht = ht.group_by('locus', 'alleles').aggregate(
            BETA=hl.agg.sum(ht.weight * ht.BETA) / hl.agg.sum(ht.weight),
            SE=1 / hl.sqrt(hl.agg.sum(ht.weight))
        )
        meta_ht.describe()
        ht = ht.annotate(META_BETA=meta_ht[ht.locus, ht.alleles].BETA)
        ht.describe()
        q_ht = ht.group_by('locus', 'alleles').aggregate(
            Het_Q=hl.agg.sum(ht.weight * (ht.BETA - ht.META_BETA) ** 2),
            AC_Allele2 = hl.agg.sum(ht.AC_Allele2),
            N = hl.agg.sum(ht.N),
            N_case = hl.agg.sum(ht.N_case),
            N_ctrl = hl.agg.sum(ht.N_ctrl),
            AC_case = hl.agg.sum(ht.AC_case),
            AC_ctrl = hl.agg.sum(ht.AC_ctrl)
        )
        q_ht = q_ht.annotate(
            AF_Allele2 = q_ht.AC_Allele2 / (q_ht.N*2)
        )
        q_ht.describe()
        meta_ht = meta_ht.annotate(**q_ht[meta_ht.key])

        meta_ht = meta_ht.annotate(Pvalue=2 * hl.pnorm(-hl.abs(meta_ht.BETA) / meta_ht.SE))
        meta_ht = meta_ht.annotate(Pvalue_log10=-hl.log10(meta_ht.Pvalue))
    # else:
    #     P_FIELDS = ['Pvalue', 'Pvalue_Burden', 'Pvalue_SKAT']
    #     P_TESTS = {'SKATO': 'Pvalue', 'Burden': 'Pvalue_Burden', 'SKAT': 'Pvalue_SKAT'}
    #     def _edit_pvalue(p):
    #         return hl.if_else(p > 0.99, 0.99, p)

    #     ht = ht.filter(
    #         hl.is_defined(ht[P_FIELDS[0]]) &
    #         hl.is_defined(ht[P_FIELDS[1]]) &
    #         hl.is_defined(ht[P_FIELDS[2]]))
    #     ht = ht.annotate(
    #         **{p_field: _edit_pvalue(ht[p_field]) for p_field in P_FIELDS},
    #         N = hl.if_else(ht.n_controls ==0, hl.float64(ht.n_cases), (4.0 * ht.n_cases * ht.n_controls) / (ht.n_cases + ht.n_controls))
    #     )

    #     def _stouffer_test(ht, test):
    #         print(f'Meta analyzing {test} results...')
    #         two_tail = test == 'Burden'
    #         beta_field = f'BETA_{test}'
    #         p_field = P_TESTS[test]

    #         if two_tail:
    #             ht = ht.annotate(**{p_field: ht[p_field] / 2},
    #                                      **{beta_field: ht[beta_field]})
    #         else:
    #             ht = ht.annotate(
    #                 **{beta_field: hl.int(hl.is_defined(ht[p_field]))})
    #         ht = ht.annotate(**{f'weighted_Z_numerator_{test}': (hl.sqrt(ht.N) * (-hl.qnorm(ht[p_field])) * hl.sign(ht[beta_field]))})
    #         meta_ht = ht.group_by('gene_id', 'gene_symbol', 'annotation').aggregate(
    #             interval=hl.agg.collect(ht.interval),
    #             n_cases=hl.agg.sum(ht.n_cases),
    #             n_controls=hl.agg.sum(ht.n_controls),
    #             **{f'META_Stats_{test}': hl.agg.sum(ht[f'weighted_Z_numerator_{test}']) / (hl.sqrt(hl.agg.sum(ht.N)))},
    #         )

    #         if two_tail:
    #             meta_ht = meta_ht.annotate(**{f'{p_field}': 2 * hl.pnorm(hl.abs(meta_ht[f'META_Stats_{test}']), lower_tail=False)})
    #         else:
    #             meta_ht = meta_ht.annotate(**{f'{p_field}': hl.pnorm(meta_ht[f'META_Stats_{test}'], lower_tail=False)})
    #         meta_ht = meta_ht.annotate(**{f'{p_field}_log10': -hl.log10(meta_ht[p_field])})
    #         return meta_ht

    #     from datetime import datetime

    #     timestamp = datetime.now().strftime("%y%m%d_%H%M%S")
    #     meta_ht3.write(f'{TMP_BUCKET}/ukb_meta_tmp_{aou_pheno}_Burden_{timestamp}_1234566.ht')
    #     meta_ht3 = hl.read_table(f'{TMP_BUCKET}/ukb_meta_tmp_{aou_pheno}_Burden_{timestamp}_1234566.ht')
    #     meta_ht1 = _stouffer_test(ht, 'SKATO')
    #     meta_ht1.write(f'{TMP_BUCKET}/ukb_meta_tmp_{aou_pheno}_SKATO_{timestamp}_1234566.ht')
    #     meta_ht1 = hl.read_table(f'{TMP_BUCKET}/ukb_meta_tmp_{aou_pheno}_SKATO_{timestamp}_1234566.ht')
    #     meta_ht2 = _stouffer_test(ht, 'SKAT').drop('interval',  'n_cases', 'n_controls')
    #     meta_ht2.write(f'{TMP_BUCKET}/ukb_meta_tmp_{aou_pheno}_SKAT_{timestamp}_1234566.ht')    
    #     meta_ht2 = hl.read_table(f'{TMP_BUCKET}/ukb_meta_tmp_{aou_pheno}_SKAT_{timestamp}_1234566.ht')
    #     meta_ht3 = _stouffer_test(ht, 'Burden').drop('interval',  'n_cases', 'n_controls')

    #     meta_ht = meta_ht1.annotate(
    #         **meta_ht2[meta_ht1.key],
    #         **meta_ht3[meta_ht1.key],
    #     )

    # meta_ht.describe()
    # meta_ht.show()
    if analysis_type == 'gene':
        meta_ht = meta_ht.annotate(CHR=meta_ht.interval[0].start.contig, POS=meta_ht.interval[0].start.position)
    else:
        meta_ht = meta_ht.annotate(CHR=meta_ht.locus.contig, POS=meta_ht.locus.position)
    meta_path = f'{ANALYSIS_BUCKET}/v8/ht_results/aou_ukb_META/phenotype_{aou_pheno}/{table_name}.ht'
    meta_ht.write(meta_path, overwrite=overwrite)
    meta_ht = hl.read_table(meta_path)
    meta_ht.export(meta_path.replace('.ht', '.txt.bgz'))
    meta_ht.describe()
    print(meta_ht.count())
    # meta_ht.describe()
    # print(meta_ht.count())
    # print(meta_ht.aggregate(hl.agg.counter(meta_ht.annotation)))
    # print(meta_ht.summarize())
    return meta_ht

def main(args):
    print('='*60)
    print(f'AoU x UKB cross-biobank meta-analysis')
    print(f'Analysis type: {args.analysis_type}')
    print('='*60)

    hl.init(
        tmp_dir=TMP_BUCKET,
        default_reference='GRCh38',
    )
    table_name = args.table_name
    qc_tag = '' if not args.qc_results else '_QC'
    pan_ukb_qc_tag = '' if args.pan_ukb_h2_qc else '_non_h2_qc'

    if args.map_phenotypes:
        print('\n' + '='*60)
        print('STAGE: Phenotype mapping (PhecodeX -> Phecode -> ICD-10)')
        print('='*60)
        print('Step 1: Map ICD-10 to Phecode using Pan-UKB pairwise correlations')
        ## Step 1: Map ICD-10 to Phecode using Pan-UKB pairwise correlations
        ht = hl.read_table('gs://aou_wlu/pairwise_correlations_regressed.ht')
        ht = ht.filter((ht.i_data.trait_type == 'phecode') & (ht.j_data.trait_type.startswith('icd')))
        ht = ht.filter(ht.j_data.trait_type != 'icd9')
        # Require positive correlation and r² >= 0.8
        ht = ht.filter(ht.entry > 0)
        ht = ht.annotate(r2=ht.entry ** 2)
        ht = ht.filter(ht.r2 >= 0.8)
        fields = ['trait_type', 'phenocode', 'pheno_sex', 'coding', 'modifier', 'n_cases_both_sexes', 'description']
        ht = ht.annotate(**{f'i_{field}': ht.i_data[field] for field in fields},
                        **{f'j_{field}': ht.j_data[field] for field in fields})
        ht = ht.drop('i_data', 'j_data')
        # Best ICD per phecode (keep all phecodes; final 1-to-1 dedup in R)
        map_ht = ht.select(ICD10=ht.j_phenocode, phecode=ht.i_phenocode, r2=ht.r2, r=ht.entry)
        map_ht = map_ht.group_by('phecode').aggregate(
            ICD10=hl.agg.take(map_ht.ICD10, 1, ordering=-map_ht.r2)[0],
            r=hl.agg.take(map_ht.r, 1, ordering=-map_ht.r2)[0],
            r2=hl.agg.take(map_ht.r2, 1, ordering=-map_ht.r2)[0],
        )  # best ICD per phecode
        map_ht.write('gs://aou_wlu/v8_analysis/icd_phecode_map_0.8.ht', overwrite=True)
        map_ht.export('gs://aou_wlu/v8_analysis/icd_phecode_map_0.8.txt.bgz')

        print('Step 2: Map PhecodeX to Phecode using AoU phenotype correlations')
        ## Step 2: Map PhecodeX to Phecode using AoU phenotype correlations
        phecode_mt = hl.read_matrix_table('gs://aou_analysis/v8/data/phenotype/raw/mcc2_phecode.mt')
        phecode_mt = phecode_mt.annotate_cols(category='phecode')
        phecodex_mt = hl.read_matrix_table('gs://aou_analysis/v8/data/phenotype/raw/mcc2_phecodex.mt')
        phecodex_mt = phecodex_mt.annotate_cols(category='phecodeX')
        full_disease_mt = phecode_mt.union_cols(phecodex_mt)
        print(full_disease_mt.count())
        full_disease_mt = full_disease_mt.key_cols_by('phenoname', 'category')
        full_disease_mt.describe()

        def make_pairwise_ht(mt: hl.MatrixTable, pheno_field, correlation: bool = False):
            mt = mt.annotate_entries(_pheno=pheno_field)
            mt = mt.add_col_index()
            index_ht = mt.cols().key_by('col_idx')
            if correlation:
                bm = hl.linalg.BlockMatrix.from_entry_expr(mt._pheno, mean_impute=True, center=True, normalize=True, axis='cols', block_size=1024)
            else:
                bm = hl.linalg.BlockMatrix.from_entry_expr(mt._pheno, block_size=1024)
            bm = bm.T @ bm
            pheno_ht = bm.entries()
            pheno_ht = pheno_ht.annotate(i_data=index_ht[pheno_ht.i], j_data=index_ht[pheno_ht.j])
            return pheno_ht

        corr = make_pairwise_ht(full_disease_mt, pheno_field=full_disease_mt.value, correlation=True)
        corr = corr.checkpoint('gs://aou_wlu/v8_analysis/phecode_phecodex_corr_full.ht')

        # Require positive correlation, r² >= 0.8, cross-category pairs only
        related = corr.filter((corr.entry > 0) & (corr.entry ** 2 >= 0.8) & (corr.i != corr.j))
        related = related.filter((related.i_data.category == 'phecode') & (related.j_data.category == 'phecodeX'))
        related = related.select(Phecode=related.i_data.phenoname, PhecodeX=related.j_data.phenoname, r=related.entry, r2=related.entry ** 2)
        # Best phecode per phecodeX (keep all phecodeX; final 1-to-1 dedup in R)
        related = related.group_by('PhecodeX').aggregate(
            Phecode=hl.agg.take(related.Phecode, 1, ordering=-related.r2)[0],
            r=hl.agg.take(related.r, 1, ordering=-related.r2)[0],
            r2=hl.agg.take(related.r2, 1, ordering=-related.r2)[0],
        )  # best phecode per phecodeX
        related = related.checkpoint('gs://aou_wlu/v8_analysis/phecode_phecodex_map_0.8.ht', overwrite=True)
        related.export('gs://aou_wlu/v8_analysis/phecode_phecodex_map_0.8.txt.bgz')

        print('Step 3: Combine mappings (PhecodeX -> Phecode -> ICD-10) done in supp_numbers.R')
        ## Step 3: Combine mappings to get PhecodeX -> Phecode -> ICD-10 in supp_numbers.R

    if args.analysis_type == 'gene':
        print('\n' + '='*60)
        print('STAGE: Gene-level meta-analysis')
        print('='*60)
        maxmaf = args.maxmaf

        path = f'gs://aou_wlu/v8_analysis/aou_ukb_meta/ukb_gene_formatted_eur{qc_tag}.mt'
        if args.overwrite_ukb_gene or not hl.hadoop_exists(f'{path}/_SUCCESS'):
            print('\n  Reformatting UKB Genebass gene-level MT...')
            ukb_mt = hl.read_matrix_table(UKB_GENE_MT)
            if args.qc_results:
                # Null out test-specific results where keep_gene_{test} is false
                ukb_mt = ukb_mt.annotate_entries(
                    Pvalue = hl.or_missing(ukb_mt.keep_gene_skato, ukb_mt.Pvalue),
                    Pvalue_Burden = hl.or_missing(ukb_mt.keep_gene_burden, ukb_mt.Pvalue_Burden),
                    Pvalue_SKAT = hl.or_missing(ukb_mt.keep_gene_skat, ukb_mt.Pvalue_SKAT),
                    BETA_Burden = hl.or_missing(ukb_mt.keep_gene_burden, ukb_mt.BETA_Burden),
                    SE_Burden = hl.or_missing(ukb_mt.keep_gene_burden, ukb_mt.SE_Burden),
                )
                # Apply test-agnostic gene-level QC filters
                ukb_mt = ukb_mt.filter_rows(
                    ukb_mt.keep_gene_coverage
                    # & ukb_mt.keep_gene_expected_ac
                    # & ukb_mt.keep_gene_n_var
                )
                # ukb_mt = ukb_mt.filter_entries(ukb_mt.keep_entry_expected_ac)
            ukb_mt = format_ukb_gene_mt(ukb_mt)
            overwrite = True
            ukb_mt = ukb_mt.checkpoint(path, _read_if_exists=not overwrite, overwrite=overwrite)
            ukb_mt.describe()
            print(f'  UKB Genebass gene MT count: {ukb_mt.count()}')

        path = f'gs://aou_wlu/v8_analysis/aou_ukb_meta/aou_eur_gene_formatted_0.01{qc_tag}.mt'
        if args.overwrite_aou_eur_gene or not hl.hadoop_exists(f'{path}/_SUCCESS'):
            print('\n  Reformatting AoU EUR gene-level MT...')
            aou_eur_mt = hl.read_matrix_table(AOU_EUR_GENE_MT)
            if args.qc_results:
                aou_eur_mt = aou_eur_mt.filter_rows(aou_eur_mt.hq_coverage)
                aou_eur_mt = aou_eur_mt.filter_cols(aou_eur_mt.hq_pheno)
                # aou_eur_mt = aou_eur_mt.filter_entries(aou_eur_mt.hq_exp_CAC)
            aou_eur_mt = format_aou_gene_eur_mt(aou_eur_mt, maxmaf=0.01)
            overwrite = True
            aou_eur_mt = aou_eur_mt.checkpoint(path, _read_if_exists=not overwrite, overwrite=overwrite)
            aou_eur_mt.describe()
            print(f'  AoU EUR gene MT count: {aou_eur_mt.count()}')

        path = f'gs://aou_wlu/v8_analysis/aou_ukb_meta/aou_gene_formatted_{maxmaf}{qc_tag}.mt'
        if args.overwrite_aou_gene or not hl.hadoop_exists(f'{path}/_SUCCESS'):
            print(f'\n  Reformatting AoU meta gene-level MT (maxMAF={maxmaf})...')
            aou_mt = hl.read_matrix_table(AOU_GENE_MT)
            if args.qc_results:
                aou_mt = aou_mt.filter_rows(aou_mt.quality_flags.hq_coverage)
                aou_mt = aou_mt.filter_cols(aou_mt.hq_pheno)
                # aou_mt = aou_mt.filter_entries(aou_mt.hq_exp_CAC)
            aou_mt = format_aou_gene_meta_mt(aou_mt, maxmaf=maxmaf)
            overwrite = True
            aou_mt = aou_mt.checkpoint(path, _read_if_exists=not overwrite, overwrite=overwrite)
            aou_mt.describe()
            print(f'  AoU meta gene MT (maxMAF={maxmaf}) count: {aou_mt.count()}')

        path = f'gs://aou_wlu/v8_analysis/aou_ukb_meta/full_eur_gene_formatted_0.01{qc_tag}.mt'
        if args.overwrite_merge_eur_gene or not hl.hadoop_exists(f'{path}/_SUCCESS'):
            print('\n  Merging AoU EUR + UKB gene-level MTs...')
            aou_eur_mt = hl.read_matrix_table(f'gs://aou_wlu/v8_analysis/aou_ukb_meta/aou_eur_gene_formatted_0.01{qc_tag}.mt')
            ukb_mt = hl.read_matrix_table(f'gs://aou_wlu/v8_analysis/aou_ukb_meta/ukb_gene_formatted_eur{qc_tag}.mt')
            ukb_mt = ukb_mt.select_entries('Pvalue_Burden', 'BETA_Burden', 'SE_Burden')
            mt = aou_eur_mt.union_cols(ukb_mt, row_join_type='outer')
            mt = mt.collect_cols_by_key()
            mt = mt.filter_cols(hl.len(mt.Neff) == 2)
            mt.describe()
            overwrite = True
            mt = mt.checkpoint(path, _read_if_exists=not overwrite, overwrite=overwrite)
            mt.describe()
            print(f'  Merged EUR gene MT count: {mt.count()}')

        path = f'gs://aou_wlu/v8_analysis/aou_ukb_meta/full_gene_formatted_{maxmaf}{qc_tag}.mt'
        if args.overwrite_merge_gene or not hl.hadoop_exists(f'{path}/_SUCCESS'):
            print(f'\n  Merging AoU meta + UKB gene-level MTs (maxMAF={maxmaf})...')
            aou_mt = hl.read_matrix_table(f'gs://aou_wlu/v8_analysis/aou_ukb_meta/aou_gene_formatted_{maxmaf}{qc_tag}.mt')
            ukb_mt = hl.read_matrix_table(f'gs://aou_wlu/v8_analysis/aou_ukb_meta/ukb_gene_formatted_eur{qc_tag}.mt')
            mt = aou_mt.union_cols(ukb_mt.drop('SE_Burden'), row_join_type='outer')
            mt = mt.collect_cols_by_key()
            mt = mt.filter_cols(hl.len(mt.Neff) == 2)
            mt.describe()
            overwrite = True
            mt = mt.checkpoint(path, overwrite=overwrite)
            mt.describe()
            print(f'  Merged meta gene MT (maxMAF={maxmaf}) count: {mt.count()}')

        path = f'gs://aou_wlu/v8_analysis/aou_ukb_meta/META_gene_aou_ukb_{maxmaf}{qc_tag}.mt'
        if args.overwrite_meta_gene or not hl.hadoop_exists(f'{path}/_SUCCESS'):
            print(f'\n  Running Stouffer meta-analysis (maxMAF={maxmaf})...')
            mt = hl.read_matrix_table(f'gs://aou_wlu/v8_analysis/aou_ukb_meta/full_gene_formatted_{maxmaf}{qc_tag}.mt')
            mt = run_stouffer_mt(mt)
            overwrite = True
            mt = mt.checkpoint(path, _read_if_exists=not overwrite, overwrite=overwrite)
            mt.describe()
            print(f'  Gene meta-analysis (maxMAF={maxmaf}) count: {mt.count()}')

        meta_result_qc_tag = '' if not args.qc_meta_results else '_META_QCed'
        path = f'gs://aou_wlu/v8_analysis/aou_ukb_meta/significant_in_meta_only_{maxmaf}{qc_tag}{meta_result_qc_tag}.txt.bgz'
        if args.overwrite_processed_meta_gene or not hl.hadoop_exists(f'{path}'):
            print(f'\n  Identifying significant-in-meta-only results (maxMAF={maxmaf})...')
            meta_mt = hl.read_matrix_table(f'gs://aou_wlu/v8_analysis/aou_ukb_meta/META_gene_aou_ukb_{maxmaf}{qc_tag}.mt')          
            if args.qc_meta_results:
                meta_mt = meta_mt.annotate_entries(meta_CAC = hl.sum(hl.map(lambda caf, n_cases: caf * n_cases, meta_mt.CAF, meta_mt.n_cases)))
                meta_mt = meta_mt.filter_entries((hl.sum(meta_mt.MAC) >=5)& (meta_mt.meta_CAC >=5))
            meta_mt = meta_mt.select_entries('META_Stats_SKATO', 'META_Pvalue_SKATO', 'META_Stats_Burden', 'META_Pvalue_Burden', 'META_Stats_SKAT', 'META_Pvalue_SKAT')
            print(meta_mt.count())
            full_mt = hl.read_matrix_table(f'gs://aou_wlu/v8_analysis/aou_ukb_meta/full_gene_formatted_{maxmaf}{qc_tag}.mt')
            print(full_mt.count())
            fields = ['MAC','Pvalue','Pvalue_Burden', 'Pvalue_SKAT', 'BETA_SKATO', 'BETA_Burden','BETA_SKAT']
            final_mt = meta_mt.annotate_entries(
                **{f'aou_{field}': full_mt[meta_mt.row_key, meta_mt.col_key][field][0] for field in fields},
                **{f'ukb_{field}': full_mt[meta_mt.row_key, meta_mt.col_key][field][1] for field in fields},
            )
            final_mt.describe()
            overwrite = True
            final_mt = final_mt.checkpoint(f'gs://aou_wlu/v8_analysis/aou_ukb_meta/ALL_3_for_analysis_gene_aou_ukb_meta_{maxmaf}{qc_tag}{meta_result_qc_tag}.mt', 
                                    _read_if_exists=not overwrite, overwrite=overwrite)
            
            threshold = 6.7e-7
            final_Mt = hl.read_matrix_table(f'gs://aou_wlu/v8_analysis/aou_ukb_meta/ALL_3_for_analysis_gene_aou_ukb_meta_{maxmaf}{qc_tag}{meta_result_qc_tag}.mt')
            final_Mt1 = final_Mt.annotate_entries(
                sig_in_meta_only_skato = (final_Mt.META_Pvalue_SKATO <= threshold) & (final_Mt.aou_Pvalue > threshold) & (final_Mt.ukb_Pvalue > threshold),
                sig_in_meta_only_skat = (final_Mt.META_Pvalue_SKAT <= threshold) & (final_Mt.aou_Pvalue_SKAT > threshold) & (final_Mt.ukb_Pvalue_SKAT > threshold),
                sig_in_meta_only_burden = (final_Mt.META_Pvalue_Burden <= threshold) & (final_Mt.aou_Pvalue_Burden > threshold) & (final_Mt.ukb_Pvalue_Burden > threshold),
            )
            final_Mt1 = final_Mt1.filter_entries(final_Mt1.sig_in_meta_only_skato | final_Mt1.sig_in_meta_only_skat | final_Mt1.sig_in_meta_only_burden)
            final_Ht1 = final_Mt1.entries()
            final_Ht1.describe()
            overwrite = True
            final_Ht1 = final_Ht1.checkpoint(f'gs://aou_wlu/v8_analysis/aou_ukb_meta/significant_in_meta_only_{maxmaf}{qc_tag}{meta_result_qc_tag}.ht', _read_if_exists=not overwrite, overwrite=overwrite)
            final_Ht1.export(f'gs://aou_wlu/v8_analysis/aou_ukb_meta/significant_in_meta_only_{maxmaf}{qc_tag}{meta_result_qc_tag}.txt.bgz')

    if args.analysis_type == 'exome':
        print('\n' + '='*60)
        print('STAGE: Exome variant-level meta-analysis')
        print('='*60)
        path = f'gs://aou_wlu/v8_analysis/aou_ukb_meta/aou_meta_exome_formatted{qc_tag}.mt'
        if args.overwrite_aou_exome or not hl.hadoop_exists(f'{path}/_SUCCESS'):
            print('\n  Reformatting AoU exome variant MT...')
            aou_exome_mt = hl.read_matrix_table(AOU_EXOME_MT)
            if args.qc_results:
                # aou_exome_mt = aou_exome_mt.filter_rows(aou_exome_mt.hq_variant)
                # aou_exome_mt = aou_exome_mt.filter_entries(aou_exome_mt.hq_exp_AC)
                aou_exome_mt = aou_exome_mt.filter_cols(aou_exome_mt.hq_pheno)
            aou_exome_mt = format_aou_exome_variant_mt(aou_exome_mt)
            overwrite = True
            aou_exome_mt = aou_exome_mt.checkpoint(path, _read_if_exists=not overwrite, overwrite=overwrite)
            aou_exome_mt.describe()
            print(f'  AoU exome MT count: {aou_exome_mt.count()}')

        path = f'gs://aou_wlu/v8_analysis/aou_ukb_meta/ukb_exome_formatted_eur{qc_tag}.mt'
        if args.overwrite_ukb_exome or not hl.hadoop_exists(f'{path}/_SUCCESS'):
            print('\n  Reformatting UKB exome variant MT...')
            ukb_exome_mt = hl.read_matrix_table(UKB_EXOME_MT)
            if args.qc_results:
                # ukb_exome_mt = ukb_exome_mt.annotate_entries(expected_AC=ukb_exome_mt.AF * ukb_exome_mt.n_cases_defined)
                # ukb_exome_mt = ukb_exome_mt.annotate_rows(
                #     expected_ac_row_filter=hl.agg.count_where(
                #         ukb_exome_mt.expected_AC >= 5
                #     )
                # )
                ukb_exome_mt = ukb_exome_mt.annotate_rows(
                    annotation=hl.if_else(
                        hl.literal({"missense", "LC"}).contains(ukb_exome_mt.annotation),
                        "missense|LC",
                        ukb_exome_mt.annotation,
                    ),
                    # keep_var_expected_ac=ukb_exome_mt.expected_ac_row_filter > 0,
                    keep_var_annt=hl.is_defined(ukb_exome_mt.annotation),
                )
                # ukb_exome_mt = ukb_exome_mt.annotate_entries(keep_entry_expected_ac = ukb_exome_mt.expected_AC >= 5)
                # ukb_exome_mt = ukb_exome_mt.filter_rows(ukb_exome_mt.keep_var_expected_ac & ukb_exome_mt.keep_var_annt)
                # ukb_exome_mt = ukb_exome_mt.filter_entries(ukb_exome_mt.keep_entry_expected_ac)
                # ukb_exome_mt = ukb_exome_mt.filter_rows(ukb_exome_mt.keep_var_annt)
            ukb_exome_mt = format_ukb_exome_variant_mt(ukb_exome_mt)
            overwrite = True
            ukb_exome_mt = ukb_exome_mt.checkpoint(path, _read_if_exists=not overwrite, overwrite=overwrite)
            ukb_exome_mt.describe()
            print(f'  UKB exome MT count: {ukb_exome_mt.count()}')

        path = f'gs://aou_wlu/v8_analysis/aou_ukb_meta/full_exome_formatted{qc_tag}.mt'
        if args.overwrite_merge_exome or not hl.hadoop_exists(f'{path}/_SUCCESS'):
            print('\n  Merging AoU + UKB exome MTs...')
            aou_exome_mt = hl.read_matrix_table(f'gs://aou_wlu/v8_analysis/aou_ukb_meta/aou_meta_exome_formatted{qc_tag}.mt')
            ukb_exome_mt = hl.read_matrix_table(f'gs://aou_wlu/v8_analysis/aou_ukb_meta/ukb_exome_formatted_eur{qc_tag}.mt')
            mt = aou_exome_mt.union_cols(ukb_exome_mt, row_join_type='outer')
            mt = mt.collect_cols_by_key()
            mt = mt.filter_cols(hl.len(mt.N) == 2)
            mt.describe()
            overwrite = True
            mt = mt.checkpoint(path, _read_if_exists=not overwrite, overwrite=overwrite)
            mt.describe()
            print(f'  Merged exome MT count: {mt.count()}')

        path = f'gs://aou_wlu/v8_analysis/aou_ukb_meta/META_exome_aou_ukb{qc_tag}.mt'
        if args.overwrite_meta_exome or not hl.hadoop_exists(f'{path}/_SUCCESS'):
            print('\n  Running exome IVW meta-analysis...')
            mt = hl.read_matrix_table(f'gs://aou_wlu/v8_analysis/aou_ukb_meta/full_exome_formatted{qc_tag}.mt')
            mt = run_meta_analysis(mt)
            overwrite = True
            mt = mt.checkpoint(path, _read_if_exists=not overwrite, overwrite=overwrite)
            mt.describe()
            print(f'  Exome meta-analysis count: {mt.count()}')

    if args.analysis_type == 'genome':
        print('\n' + '='*60)
        print('STAGE: Genome variant-level meta-analysis')
        print('='*60)
        path = f'gs://aou_wlu/v8_analysis/aou_ukb_meta/aou_meta_genome_formatted{qc_tag}.mt'
        if args.overwrite_aou_genome or not hl.hadoop_exists(f'{path}/_SUCCESS'):
            print('\n  Reformatting AoU genome variant MT...')
            aou_genome_mt = hl.read_matrix_table(AOU_GENOME_MT)
            if args.qc_results:
                # aou_genome_mt = aou_genome_mt.filter_rows(aou_genome_mt.hq_variant)
                # aou_genome_mt = aou_genome_mt.filter_entries(aou_genome_mt.hq_exp_AC)
                aou_genome_mt = aou_genome_mt.filter_cols(aou_genome_mt.hq_pheno)
            aou_genome_mt = format_aou_genome_variant_mt(aou_genome_mt)
            print(f'  AoU genome MT count (pre-checkpoint): {aou_genome_mt.count()}')
            overwrite = True
            aou_genome_mt = aou_genome_mt.checkpoint(path, _read_if_exists=not overwrite, overwrite=overwrite)
            aou_genome_mt.describe()
            print(f'  AoU genome MT count: {aou_genome_mt.count()}')
        
        path = f'gs://aou_wlu/v8_analysis/aou_ukb_meta/ukb_genome_formatted_eur{qc_tag}{pan_ukb_qc_tag}.mt'
        if args.overwrite_ukb_genome or not hl.hadoop_exists(f'{path}/_SUCCESS'):
            hl._set_flags(shuffle_max_branch_factor=32)
            print('\n  Reformatting UKB genome variant MT (includes liftover GRCh37->GRCh38)...')
            if args.pan_ukb_h2_qc:
                print('  Applying Pan-UKB H2 QC filters to UKB genome MT...')
                UKB_GENOME_MT = UKB_QC_GENOME_MT
            else:
                UKB_GENOME_MT = UKB_RAW_GENOME_MT
            ukb_genome_mt = hl.read_matrix_table(UKB_GENOME_MT)
            if args.qc_results:
                ukb_genome_mt = ukb_genome_mt.filter_rows(ukb_genome_mt.high_quality)
            ukb_genome_mt = format_ukb_genome_variant_mt(ukb_genome_mt)
            overwrite = True
            ukb_genome_mt = ukb_genome_mt.checkpoint(path, _read_if_exists=not overwrite, overwrite=overwrite)
            ukb_genome_mt.describe()
            print(f'  UKB genome MT count: {ukb_genome_mt.count()}')

        path = f'gs://aou_wlu/v8_analysis/aou_ukb_meta/full_genome_formatted{qc_tag}{pan_ukb_qc_tag}.mt'
        if args.overwrite_merge_genome or not hl.hadoop_exists(f'{path}/_SUCCESS'):
            print('\n  Merging AoU + UKB genome MTs...')
            aou_genome_mt = hl.read_matrix_table(f'gs://aou_wlu/v8_analysis/aou_ukb_meta/aou_meta_genome_formatted{qc_tag}.mt')
            ukb_genome_mt = hl.read_matrix_table(f'gs://aou_wlu/v8_analysis/aou_ukb_meta/ukb_genome_formatted_eur{qc_tag}{pan_ukb_qc_tag}.mt')
            mt = aou_genome_mt.union_cols(ukb_genome_mt, row_join_type='outer')
            mt = mt.collect_cols_by_key()
            mt = mt.filter_cols(hl.len(mt.N) == 2)
            mt.describe()
            overwrite = True
            mt = mt.checkpoint(path, _read_if_exists=not overwrite, overwrite=overwrite)
            mt.describe()
            print(f'  Merged genome MT count: {mt.count()}')

        path = f'gs://aou_wlu/v8_analysis/aou_ukb_meta/META_genome_aou_ukb{qc_tag}{pan_ukb_qc_tag}.mt'
        if args.overwrite_meta_genome or not hl.hadoop_exists(f'{path}/_SUCCESS'):
            print('\n  Running genome IVW meta-analysis...')
            mt = hl.read_matrix_table(f'gs://aou_wlu/v8_analysis/aou_ukb_meta/full_genome_formatted{qc_tag}{pan_ukb_qc_tag}.mt')
            mt = run_meta_analysis(mt)
            overwrite = True
            mt = mt.checkpoint(path, _read_if_exists=not overwrite, overwrite=overwrite)
            mt.describe()
            print(f'  Genome meta-analysis count: {mt.count()}')

    if args.individual_pheno_meta:
        print('\n' + '='*60)
        print('STAGE: Per-phenotype meta-analysis')
        print('='*60)
        # Load phenotype mapping as dict {aou_phenoname: ukb_phenocode}
        pheno_map_pkl = 'v8_pheno_map.pkl'
        if os.path.exists(pheno_map_pkl):
            with open(pheno_map_pkl, 'rb') as f:
                v8_pheno_map = pickle.load(f)
            print(f'Loaded {len(v8_pheno_map)} phenotypes from {pheno_map_pkl}')
        else:
            _pheno_ht = hl.import_table(PHENO_MAP_TSV, delimiter='\t')
            _pheno_ht.describe()
            print(_pheno_ht.count())
            v8_pheno_map = dict(zip(
                _pheno_ht.phenoname.collect(),
                _pheno_ht.ukb_phenocode.collect(),
            ))
            with open(pheno_map_pkl, 'wb') as f:
                pickle.dump(v8_pheno_map, f)
            print(f'Saved {len(v8_pheno_map)} phenotypes to {pheno_map_pkl}')
        if args.batch:
            print(f'Submitting jobs via Hail Batch (billing project: {args.billing_project})')
            backend = hb.ServiceBackend(
                billing_project=args.billing_project, remote_tmpdir=TMP_BUCKET
            )
            b = hb.Batch(
                name=f"aou_ukb_meta_{table_name}",
                backend=backend,
                default_image="hailgenetics/hail:0.2.133-py3.11",
                default_storage="500Mi",
                default_cpu=16,
                default_regions=['us-central1'],
            )

        for aou_code, ukb_code in v8_pheno_map.items():
            aou_ht_path = f'{ANALYSIS_BUCKET}/ht_results/META/phenotype_{aou_code}/{table_name}.ht'
            if aou_code.startswith('random'):
                ukb_ht_path = f'{ANALYSIS_BUCKET}/ht_results/aou_ukb_META/ukb_random/{aou_code}_{args.analysis_type}_results.ht'
            else:
                ukb_ht_path = f'gs://ukbb-exome-public/500k/results/per_phenotype_hail_tables/{ukb_code}/{args.analysis_type}_results.ht'

            output_check = f'{ANALYSIS_BUCKET}/ht_results/aou_ukb_META/phenotype_{aou_code}_{ukb_code}_{table_name}_raw.txt.bgz'
            if hl.hadoop_exists(output_check) and not args.overwrite:
                print(f'Skipping {aou_code} (output exists)')
                continue

            print(f'Running {aou_code} <-> {ukb_code}')
            if args.batch:
                j = b.new_python_job(name=f"aou_ukb_meta_{aou_code}_{ukb_code}")
                j.image("hailgenetics/hail:0.2.133-py3.11")
                j.cpu(16)
                j.memory('100G')
                j.storage("20G")
                j.env('PYSPARK_SUBMIT_ARGS', '--driver-memory 100g --executor-memory 100g pyspark-shell')
                j.call(run_single_meta_analysis_ht,
                    aou_ht_path=aou_ht_path,
                    ukb_ht_path=ukb_ht_path,
                    aou_pheno=aou_code,
                    ukb_pheno=ukb_code,
                    table_name=table_name,
                    analysis_type=args.analysis_type,
                    overwrite=args.overwrite,
                    )
            else:
                run_single_meta_analysis_ht(
                    aou_ht_path=aou_ht_path,
                    ukb_ht_path=ukb_ht_path,
                    aou_pheno=aou_code,
                    ukb_pheno=ukb_code,
                    table_name=table_name,
                    analysis_type=args.analysis_type,
                    overwrite=args.overwrite,
                )

        if args.batch:
            b.run()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='AoU x UKB cross-biobank meta-analysis')
    parser.add_argument('--analysis-type', choices=['variant', 'gene', 'exome', 'genome'], default='variant',
                        help='Type of analysis: variant, gene, exome, or genome (default: variant)')
    parser.add_argument('--table-name', default='variant_results',
                        help='Table name for per-phenotype HT paths (default: variant_results)')
    parser.add_argument('--phenos',
                        help='Comma-separated phenotypes. Use aou_code:ukb_code for mapping, or just name if same in both')
    parser.add_argument('--pheno-map-csv',
                        help='CSV with phenoname and ukb_phenocode columns for phenotype mapping')
    parser.add_argument('--batch', action='store_true',
                        help='Submit jobs via Hail Batch instead of running locally')
    parser.add_argument('--billing-project', default='all-by-aou',
                        help='Hail Batch billing project (default: all-by-aou)')
    parser.add_argument('--overwrite', action='store_true',
                        help='Overwrite existing results')
    parser.add_argument('--qc-results', action='store_true',
                        help='Apply QC filters to genome-wide variant MTs and only output QC-passing results (overwrites --qc)')

    # Phenotype mapping
    parser.add_argument('--map-phenotypes', action='store_true',
                        help='Run phenotype mapping (PhecodeX->Phecode->ICD-10)')

    # Gene-level overwrite flags
    parser.add_argument('--maxmaf', type=float, default=0.001,
                        help='Maximum MAF threshold for gene-level analysis (default: 0.001)')
    parser.add_argument('--overwrite-ukb-gene', action='store_true',
                        help='Force reformat UKB Genebass gene-level MT')
    parser.add_argument('--overwrite-aou-eur-gene', action='store_true',
                        help='Force reformat AoU EUR gene-level MT')
    parser.add_argument('--overwrite-aou-gene', action='store_true',
                        help='Force reformat AoU meta gene-level MT')
    parser.add_argument('--overwrite-merge-eur-gene', action='store_true',
                        help='Force redo AoU EUR + UKB gene-level merge')
    parser.add_argument('--overwrite-merge-gene', action='store_true',
                        help='Force redo AoU meta + UKB gene-level merge')
    parser.add_argument('--overwrite-meta-gene', action='store_true',
                        help='Force redo gene-level Stouffer meta-analysis')
    parser.add_argument('--qc-meta-results', action='store_true',
                        help='Apply QC filters to meta-analysis results (e.g. MAC >= 5)')
    parser.add_argument('--overwrite-processed-meta-gene', action='store_true',
                        help='Force overwrite of processed meta-analysis results with QC filters applied')

    # Exome-level overwrite flags
    parser.add_argument('--overwrite-aou-exome', action='store_true',
                        help='Force reformat AoU exome variant MT')
    parser.add_argument('--overwrite-ukb-exome', action='store_true',
                        help='Force reformat UKB exome variant MT')
    parser.add_argument('--overwrite-merge-exome', action='store_true',
                        help='Force redo AoU + UKB exome merge')
    parser.add_argument('--overwrite-meta-exome', action='store_true',
                        help='Force redo exome meta-analysis')

    # Genome-level overwrite flags
    parser.add_argument('--overwrite-aou-genome', action='store_true',
                        help='Force reformat AoU genome variant MT')
    parser.add_argument('--overwrite-ukb-genome', action='store_true',
                        help='Force reformat UKB genome variant MT')
    parser.add_argument('--overwrite-merge-genome', action='store_true',
                        help='Force redo AoU + UKB genome merge')
    parser.add_argument('--overwrite-meta-genome', action='store_true',
                        help='Force redo genome meta-analysis')
    parser.add_argument('--pan-ukb-h2-qc', action='store_true',
                        help='Apply Pan-UKB H2 QC filters to UKB genome MT (overwrites --qc_results)')

    # Per-phenotype meta-analysis
    parser.add_argument('--individual-pheno-meta', action='store_true',
                        help='Run per-phenotype meta-analysis using HailTables')

    args = parser.parse_args()
    main(args)