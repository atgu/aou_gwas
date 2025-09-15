# ALL PATHs in the workbech: https://support.researchallofus.org/hc/en-us/articles/4616869437204-Controlled-CDR-Directory
# Descriptions: https://support.researchallofus.org/hc/en-us/articles/4614687617556-How-the-All-of-Us-Genomic-data-are-organized
import hail as hl

##################### ROOT PATHs #####################
TRANCHE = "250k"
RESULTS_BUCKET = "gs://aou_results"
ANALYSIS_BUCKET = "gs://aou_analysis"
MY_BUCKET = "gs://aou_wlu"
TMP_BUCKET = "gs://aou_tmp"
PREVIOUS_PATH = f"{ANALYSIS_BUCKET}/{TRANCHE}"
CURRENT_PATH = f"{ANALYSIS_BUCKET}/{TRANCHE}/v1.1"
DATA_PATH = f"{ANALYSIS_BUCKET}/{TRANCHE}/data"
RESULTS_PATH = f'{RESULTS_BUCKET}/{TRANCHE}/v1.1'

##################### PARAMETERs #####################
CALLRATE_CUTOFF = 0.9
N_GENE_PER_GROUP = 200
CHUNK_SIZE = {'all': int(5e6),'eur': int(2.5e7), "afr": int(5e7), "amr": int(5e7), "eas": int(5e7), "mid": int(5e7), "sas": int(5e7)}
REFERENCE = "GRCh38"
CHROMOSOMES = list(map(str, range(1, 23))) + ["X", "Y"]
SAIGE_DOCKER_IMAGE = "wzhou88/saige:1.3.0"
R_DOCKER_IMAGE = "us-central1-docker.pkg.dev/aou-neale-gwas/plot/aou_man_qq:1.2"
PILOT_PHENOTYPES = set(["height", "250.2", "411", "495", "585.3", "3028288", "random_0.5_continuous_1", "random_0.5_0.01_1"])
P_VALUE_FIELDS = {"skato": "Pvalue", "skat": "Pvalue_SKAT", "burden": "Pvalue_Burden"}

phenotype_categories_in_mt = ['mcc2_phecode_table', 'new_phecode_table_Dec23', 'r_drug_table', 'pfhh_survey_table',
                              "processed_physical_measurement_table"]
phenotype_categories_in_ht = ["lab_measurements", 'random_phenotypes']
phenotype_categories = phenotype_categories_in_mt + phenotype_categories_in_ht
quantitative_categories = ["lab_measurements", "processed_physical_measurement_table"]
binary_categories = ['mcc2_phecode_table', 'new_phecode_table_Dec23', 'r_drug_table', 'pfhh_survey_table']
CATEGORIES = {'new_phecode_table_Dec23': 'phecodeX',
               'mcc2_phecode_table': 'phecode',
               'processed_physical_measurement_table': 'physical_measurement',
               'r_drug_table': 'r_drug',
               'pfhh_survey_table': 'pfhh_survey',
               'lab_measurements': 'lab_measurement',
               'random_phenotypes': 'random'}

LAB_TO_REMOVE = ['3000963', '3006906', '3023103', '3016723']

PHENOTYPE_TO_REMOVE = {
    'afr':['654'] + LAB_TO_REMOVE,
    'amr':['654'] + LAB_TO_REMOVE,
    'eas':LAB_TO_REMOVE,
    'eur':['638', '654', '665', 'PP_927', 'PP_927.2'] + LAB_TO_REMOVE,
    'mid': LAB_TO_REMOVE,
    'sas':LAB_TO_REMOVE,
}

##################### POPULATION SUMMARY (Computed in the AoU workbench):
# ANCESTRY_PATH = 'gs://fc-aou-datasets-controlled/v7/wgs/short_read/snpindel/aux/ancestry/ancestry_preds.tsv'
# pop_ht = hl.import_table(ANCESTRY_PATH, key="research_id", impute=True,
#                              types={"research_id": "tstr", "pca_features": hl.tarray(hl.tfloat)})
# pop_ht.aggregate(hl.agg.counter(pop_ht.ancestry_pred))
POPS = ("afr", "amr", "eas", "eur", "mid", "sas")
N_SAMPLES_RAW = {
    "all": 245394,
    "afr": 56913,
    "amr": 45035,
    "eas": 5706,
    "eur": 133581,
    "mid": 942,
    "sas": 3217,
}

N_SAMPLES_HARD_FILTERED = {
    "all": 239734,
    "afr": 55166,
    "amr": 44143,
    "eas": 5638,
    "eur": 130697,
    "mid": 914,
    "sas": 3176,
}

N_SAMPLES_PRUNED = {
    "all": 214216,
    "afr": 50457,
    "amr": 40437,
    "eas": 5344,
    "eur": 114343,
    "mid": 691,
    "sas": 2944,
}

##################### ORIGINAL PATHs #####################
ORIGINAL_GENO_ROOT = "gs://prod-drc-broad"
ORIGINAL_PHENO_ROOT = "gs://allxall-phenotypes/data"

ORIGINAL_PHENO_PATHs = [
    f"{ORIGINAL_PHENO_ROOT}/top_10_labs.csv",
    f"{ORIGINAL_PHENO_ROOT}/top100_labs_concept_ids.csv",
    f"{ORIGINAL_PHENO_ROOT}/demographics_table.csv",
    f"{ORIGINAL_PHENO_ROOT}/physical_measurement_table.csv",
    f"{ORIGINAL_PHENO_ROOT}/new_drug_table.csv",
    f"{ORIGINAL_PHENO_ROOT}/r_drug_table.csv",
    f"{ORIGINAL_PHENO_ROOT}/pfhh_survey_table.csv",
    f"{ORIGINAL_PHENO_ROOT}/mcc2_phecode_table.csv",
    f"{ORIGINAL_PHENO_ROOT}/new_phecode_table_Dec23.csv"
]

LAB_MEASUREMENT_PATHs = [f'{ORIGINAL_PHENO_ROOT}/measurement_3000330.csv',
 f'{ORIGINAL_PHENO_ROOT}/measurement_3001008.csv',
 f'{ORIGINAL_PHENO_ROOT}/measurement_3001122.csv',
 f'{ORIGINAL_PHENO_ROOT}/measurement_3001420.csv',
 f'{ORIGINAL_PHENO_ROOT}/measurement_3002400.csv',
 f'{ORIGINAL_PHENO_ROOT}/measurement_3004501.csv',
 f'{ORIGINAL_PHENO_ROOT}/measurement_3006923.csv',
 f'{ORIGINAL_PHENO_ROOT}/measurement_3007070.csv',
 f'{ORIGINAL_PHENO_ROOT}/measurement_3007124.csv',
 f'{ORIGINAL_PHENO_ROOT}/measurement_3007238.csv',
 f'{ORIGINAL_PHENO_ROOT}/measurement_3007359.csv',
 f'{ORIGINAL_PHENO_ROOT}/measurement_3008342.csv',
 f'{ORIGINAL_PHENO_ROOT}/measurement_3008364.csv',
 f'{ORIGINAL_PHENO_ROOT}/measurement_3008598.csv',
 f'{ORIGINAL_PHENO_ROOT}/measurement_3009201.csv',
 f'{ORIGINAL_PHENO_ROOT}/measurement_3009744.csv',
 f'{ORIGINAL_PHENO_ROOT}/measurement_3010156.csv',
 f'{ORIGINAL_PHENO_ROOT}/measurement_3011099.csv',
 f'{ORIGINAL_PHENO_ROOT}/measurement_3011904.csv',
 f'{ORIGINAL_PHENO_ROOT}/measurement_3011948.csv',
 f'{ORIGINAL_PHENO_ROOT}/measurement_3013682.csv',
 f'{ORIGINAL_PHENO_ROOT}/measurement_3013721.csv',
 f'{ORIGINAL_PHENO_ROOT}/measurement_3013861.csv',
 f'{ORIGINAL_PHENO_ROOT}/measurement_3013869.csv',
 f'{ORIGINAL_PHENO_ROOT}/measurement_3014576.csv',
 f'{ORIGINAL_PHENO_ROOT}/measurement_3015183.csv',
 f'{ORIGINAL_PHENO_ROOT}/measurement_3015632.csv',
 f'{ORIGINAL_PHENO_ROOT}/measurement_3016293.csv',
 f'{ORIGINAL_PHENO_ROOT}/measurement_3016407.csv',
 f'{ORIGINAL_PHENO_ROOT}/measurement_3016436.csv',
 f'{ORIGINAL_PHENO_ROOT}/measurement_3018311.csv',
 f'{ORIGINAL_PHENO_ROOT}/measurement_3019550.csv',
 f'{ORIGINAL_PHENO_ROOT}/measurement_3019897.csv',
 f'{ORIGINAL_PHENO_ROOT}/measurement_3020149.csv',
 f'{ORIGINAL_PHENO_ROOT}/measurement_3020460.csv',
 f'{ORIGINAL_PHENO_ROOT}/measurement_3020630.csv',
 f'{ORIGINAL_PHENO_ROOT}/measurement_3022192.csv',
 f'{ORIGINAL_PHENO_ROOT}/measurement_3022621.csv',
 f'{ORIGINAL_PHENO_ROOT}/measurement_3023314.csv',
 f'{ORIGINAL_PHENO_ROOT}/measurement_3023599.csv',
 f'{ORIGINAL_PHENO_ROOT}/measurement_3024929.csv',
 f'{ORIGINAL_PHENO_ROOT}/measurement_3026910.csv',
 f'{ORIGINAL_PHENO_ROOT}/measurement_3027114.csv',
 f'{ORIGINAL_PHENO_ROOT}/measurement_3027597.csv',
 f'{ORIGINAL_PHENO_ROOT}/measurement_3028288.csv',
 f'{ORIGINAL_PHENO_ROOT}/measurement_3034426.csv',
 f'{ORIGINAL_PHENO_ROOT}/measurement_3034884.csv',
 f'{ORIGINAL_PHENO_ROOT}/measurement_3035124.csv',
 f'{ORIGINAL_PHENO_ROOT}/measurement_3035583.csv',
 f'{ORIGINAL_PHENO_ROOT}/measurement_3035774.csv',
 f'{ORIGINAL_PHENO_ROOT}/measurement_3035995.csv',
 f'{ORIGINAL_PHENO_ROOT}/measurement_3037556.csv',
 f'{ORIGINAL_PHENO_ROOT}/measurement_3043111.csv',
 f'{ORIGINAL_PHENO_ROOT}/measurement_3044491.csv',
 f'{ORIGINAL_PHENO_ROOT}/measurement_3045716.csv',
 f'{ORIGINAL_PHENO_ROOT}/measurement_3013466_3018677.csv',
 f'{ORIGINAL_PHENO_ROOT}/measurement_3043359_3027450.csv',
 f'{ORIGINAL_PHENO_ROOT}/measurement_3045440_3023939.csv',
 f'{ORIGINAL_PHENO_ROOT}/measurement_3045792_3016921.csv',
 f'{ORIGINAL_PHENO_ROOT}/measurement_3021614_3024763_3015688.csv']


LAB_CODES = [path.replace('.csv', '').replace('gs://allxall-phenotypes/data/measurement_', '') for path in LAB_MEASUREMENT_PATHs]

NEW_AGE_PHENO_PATHs = [
    f"{ORIGINAL_PHENO_ROOT}/all_first_last_billing_ages.csv",
    f"{ORIGINAL_PHENO_ROOT}/new_mcc2_phecode_table1.csv",
    f"{ORIGINAL_PHENO_ROOT}/new_mcc2_phecode_table2.csv",
    f"{ORIGINAL_PHENO_ROOT}/new_mcc2_phecode_table3.csv",
    f"{ORIGINAL_PHENO_ROOT}/new_mcc2_phecode_table4.csv",
    f"{ORIGINAL_PHENO_ROOT}/new_mcc2_phecode_table5.csv",
    f"{ORIGINAL_PHENO_ROOT}/new_mcc2_phecode_table6.csv"
]

TRAIT_TYPEs = ["physical_measurement", "r_drug", "pfhh_survey", "phecode", "phecodeX", 'lab_measurement', 'random_phenotypes']

ORIGINAL_CALLSET_PATH = f"{ORIGINAL_GENO_ROOT}/aou-wgs-delta-small_callsets_gq0/v7.1"
ACAF_MT_PATH = f"{ORIGINAL_CALLSET_PATH}/acaf_threshold_v7.1/splitMT/delta_basis_without_ext_aian_prod_gq0_3regions.acaf_threshold.split.mt"
EXOME_MT_PATH = f"{ORIGINAL_CALLSET_PATH}/exome_v7.1/splitMT/delta_basis_without_ext_aian_prod_gq0_3regions.exome.split.mt" # (34807589, 245394)
VDS_PATH = f"{ORIGINAL_GENO_ROOT}/v7/wgs/without_ext_aian_prod/vds/aou_srwgs_short_variants_v7_without_ext_aian_prod.vds"

ORIGINAL_UTIL_PATH = f"{ORIGINAL_GENO_ROOT}/aou-wgs-delta-aux_gq0"
ORIGINAL_POP_PATH = (
    f"{ORIGINAL_UTIL_PATH}/ancestry/delta_v1_gt_no_ext_aian_gq0_prod.ancestry_preds.tsv"
)
ORIGINAL_RELATED_INFO_PATH = (
    f"{ORIGINAL_UTIL_PATH}/relatedness/delta_v1_gt_no_ext_aian_gq0_prod.relatedness.tsv"
)
ORIGINAL_0th_RELATED_PATH = (
    'gs://aou_analysis/250k/data/utils/aou_0_degree_to_remove_manual_maximal_independent_jan2024.txt'
)
ORIGINAL_RELATED_SAMPLE_PATH = f"{ORIGINAL_UTIL_PATH}/relatedness/delta_v1_gt_no_ext_aian_gq0_prod.relatedness_flagged_samples.tsv"
SAMPLE_INFO_PATHs = [
    ORIGINAL_POP_PATH,
    ORIGINAL_RELATED_INFO_PATH,
    ORIGINAL_RELATED_SAMPLE_PATH,
]
table_keys = {
    "ancestry_preds": ["research_id"],
    "relatedness": ["i.s", "j.s"],
    "relatedness_flagged_samples": ["sample_id"],
}

CENTROID_PRUNED_SAMPLES = f"{DATA_PATH}/utils/pca/aou_pops_centroid_pruned.tsv"
GTF_PATH = "gs://hail-common/references/gencode/gencode.v29.annotation.gtf.bgz"
PHECODE_PATH = f'{DATA_PATH}/phenotype/phecode_definition_cleaned.ht'

##################### PHENOTYPE & SAMPLE PATHs #####################
def get_raw_phenotype_path(
    name: str, parsed: bool = True, extension: str = "ht", tranche: str = TRANCHE
):
    assert name in [
        "age",
        "r_drug_table",
        "pfhh_survey_table",
        "mcc2_phecode_table",
        "new_phecode_table_Dec23",
        "demographics_table",
        "physical_measurement_table",
        "processed_physical_measurement_table",
        "top_10_labs",
        "top100_labs_concept_ids",
        "labs_metadata_04242024",
        "labs_metadata_n_cases",
        "all_first_last_billing_ages",
        "new_drug_table",
        "new_mcc2_phecode_table1",
        "new_mcc2_phecode_table2",
        "new_mcc2_phecode_table3",
        "new_mcc2_phecode_table4",
        "new_mcc2_phecode_table5",
        "new_mcc2_phecode_table6",
        'random_phenotypes'
    ], 'Name has to be from ["age", "r_drug_table", "pfhh_survey_table", "mcc2_phecode_table", "demographics_table", "physical_measurement_table", "top_10_labs", "top100_labs_concept_ids", "labs_metadata_04242024", "labs_metadata_n_cases", "all_first_last_billing_ages", ' \
       '"new_drug_table", "new_mcc2_phecode_table1", "new_mcc2_phecode_table2", "new_mcc2_phecode_table3", "new_mcc2_phecode_table4", "new_mcc2_phecode_table5", "new_mcc2_phecode_table6", "new_phecode_table_Dec23", "random_phenotypes"]'
    return (
        f'{DATA_PATH}/phenotype/{"" if parsed else "raw/"}{name}_{tranche}.{extension}'
    )


def get_lab_measurements_path(
    name: str, parsed: bool, extension: str = "ht", tranche: str = TRANCHE
):
    return (
        f'{DATA_PATH}/phenotype/lab_measurements/{"" if parsed else "raw/"}measurement_{name}_{tranche}.{extension}'
    )


def get_random_phenotype_path(pop: str, extension: str = "tsv"):
    return f'{DATA_PATH}/phenotype/random_phenos/random_pheno_{pop}.{extension}'


BINARY_HT_PATHs = [
    get_raw_phenotype_path(name="r_drug_table", parsed=True),
    get_raw_phenotype_path(name="pfhh_survey_table", parsed=True),
    get_raw_phenotype_path(name="mcc2_phecode_table", parsed=True),
]
DEMOGRAPHICS_HT_PATH = get_raw_phenotype_path(name="demographics_table", parsed=True)
QUANTITATIVE_HT_PATH = get_raw_phenotype_path(
    name="physical_measurement_table", parsed=True
)
NEW_QUANTITATIVE_HT_PATH = get_raw_phenotype_path(
    name="processed_physical_measurement_table", parsed=True
)


def get_full_phenotype_path(
    annotation: bool,
    random_pheno: bool = False,
    extension: str = "ht",
    tranche: str = TRANCHE,
):
    return f'{DATA_PATH}/phenotype/aou_all{"_random" if random_pheno else ""}_phenotypes{"_annotated" if annotation else ""}_{tranche}.{extension}'


def get_phenotype_info_path(
    version: str,
    extension: str = "ht",
    tranche: str = TRANCHE,
):
    assert version in [
        "R", # manully computed raw phenotype summary, which includes the original version of phenotype_definitions_v2.1
        "new_phecode_table_Dec23",
        "pfhh_survey_table",
        "mcc2_phecode_table",
        "r_drug_table",
        "processed_physical_measurement_table",
        "pheno_dict",
        "trait_type_dict",
        "category_dict",
        "pheno_by_pop_by_group_dict",
        "pheno_by_pop_dict",
        "lab_by_pop_dict"
    ], 'Name has to be from ["R", "new_phecode_table_Dec23", "new_drug_table", "pfhh_survey_table", "mcc2_phecode_table", "r_drug_table", "processed_physical_measurement_table", "trait_type_dict", "category_dict", "pheno_by_pop_by_group_dict", "pheno_by_pop_dict", "lab_by_pop_dict"]'

    return f'{DATA_PATH}/phenotype/aou_pheno_info_{version}_{tranche}.{extension}'


def get_sample_meta_path(
    annotation: bool, extension: str = "ht", tranche: str = TRANCHE
):
    return f'{DATA_PATH}/utils/aou_sample_meta_info{"_annotated" if annotation else ""}_{tranche}.{extension}'


def get_aou_util_path(
    name: str, parsed: bool = False, extension: str = "ht", tranche: str = TRANCHE
):
    assert name in [
        "ancestry_preds",
        "relatedness",
        "relatedness_flagged_samples",
        "vat",
        "vat_collected_by_key",
        "vat_collected_by_key_indel",
        "mt_sample_qc",
        "acaf_mt_sample_qc",
        "vds_sample_qc",
        "sample_qc",
        "variant_qc"
        "sample_qc_tmp",
        "vep_corrected",
        "vep_old",
        "vep", # indel vep table
        "vep_full",
        "miRNA_info",
        "miRNA_raw",
        "miRNA_split",
        "miRNA_split_repartitioned",
        "miRNA_split_filtered"
    ], 'Name has to be from ["ancestry_preds", "relatedness", "relatedness_flagged_samples", "vat", "vat_collected_by_key", "vat_collected_by_key_indel", "mt_sample_qc", "vds_sample_qc","sample_qc", "variant_qc", "sample_qc_tmp", "vep_corrected", "vep_old", "vep", "vep_full"]'
    tag1 = (name == "ancestry_preds") and parsed
    tag2 = (name == "gene_map") and parsed
    if name == "vat":
        # gsutil -m cp -r gs://aou_wlu/utils/vat/hail_0_2_107/aou_PARSED_SORTED_vat.ht gs://aou_analysis/250k/data/utils/aou_parsed_and_sorted_vat_hail_0_2_107_250k.ht
        name = "parsed_and_sorted_vat_hail_0_2_107"
    return f'{DATA_PATH}/utils/aou_{name}{"_parsed" if tag1 else ""}{"_processed" if tag2 else ""}_{tranche}.{extension}'


def get_aou_relatedness_path(extension: str = "ht", tranche: str = TRANCHE):
    return f"{DATA_PATH}/utils/aou_ibd_relatedness_{tranche}.{extension}"

def get_aou_gene_map_ht_path(pop: str, processed=False, tranche: str = TRANCHE):
    return f"{DATA_PATH}/utils/gene_map/aou_{pop}_gene_map{'_processed' if processed else ''}_{tranche}.ht"

def get_pca_ht_path(pop: str, name:str, extension: str = 'ht', tranche: str = TRANCHE):
    assert name in [
        "evals",
        "scores",
        "loadings",
        "related_scores",
        "full_scores",
        "pruned_evals",
        "pruned_scores",
        "pruned_loadings",
        "pruned_related_scores",
        "pruned_full_scores"
    ], "Name has to be from ['evals', 'scores', 'loadings', 'related_scores', 'full_scores']"
    return f"{DATA_PATH}/utils/pca/aou_{pop}_pca_{name}_{tranche}.{extension}"

def get_call_stats_ht_path(pop: str, pruned: bool, analysis_type: str, extension: str = 'ht', tranche: str = TRANCHE):
    if pop =='full':
        analysis_type = 'exome_variant' if analysis_type=='gene' else 'genome_variant'
        pruned = False
    return f"{DATA_PATH}/utils/call_stats/aou_{pop}_{analysis_type}_info{'_pruned' if pruned else ''}_{tranche}.{extension}"


def get_saige_interval_path(analysis_type, chunk_size, extension:str = 'ht', tranche: str = TRANCHE):
    return f"{DATA_PATH}/utils/intervals/aou_{analysis_type}_interval_size_{chunk_size}_{tranche}.{extension}"


##################### Random Phenotype GRM PATHs #####################
def get_aou_sites_for_grm_path(
    pop: str, extension: str, pruned: bool = False, tranche: str = TRANCHE
):
    return f'{DATA_PATH}/utils/grm/aou_{pop}_sites{"_pruned" if pruned else ""}_for_grm_{tranche}.{extension}'

def get_aou_sample_file_path(
    pop: str, tranche: str = TRANCHE
):
    return f'{DATA_PATH}/utils/grm/aou_{pop}_{tranche}.samples'



##################### Result PATHs #####################
def get_aou_saige_results_root(analysis_type: str, pop: str, name: str):
    assert analysis_type in [
        "gene",
        "variant",
    ], "Name has to be from ['gene', 'variant']"
    assert name in [
        "bgen",
        "result",
        "plot",
        "null_glmm"
    ], "Name has to be from ['bgen', 'result', 'plot', 'null_glmm']"
    return f"{CURRENT_PATH}/{analysis_type}_results/{name}/{pop.upper()}"


def get_aou_saige_utils_root(pop: str, name: str):
    assert name in [
        "pheno_file",
        "null_glmm",
    ], "Name has to be from ['pheno_file', 'null_glmm']"
    root = CURRENT_PATH if name == 'null_glmm' else PREVIOUS_PATH
    return f"{root}/{name}/{pop.upper()}"

def get_aou_final_results_path(analysis_type:str, pop:str, result_type:str='variant', extension:str='mt'):
    assert analysis_type in [
        "gene",
        "variant",
    ], "Name has to be from ['gene', 'variant']"
    assert result_type in [
        "gene",
        "variant",
    ], "Name has to be from ['gene', 'variant']"
    variant_name = 'exome_variant' if analysis_type == 'gene' else 'genome_variant'
    analysis_name = 'saige_gene' if analysis_type == 'gene' else 'saige'
    if analysis_type == 'gene' and result_type == 'gene':
        variant_name = 'gene'
    return f"{RESULTS_PATH}/{analysis_name}_results/final_mt/{pop.upper()}_{variant_name}_results.{extension}"

def get_aou_analysis_results_path(analysis_type:str, pop:str, output_tag:str, result_type:str='variant', extension:str='mt'):
    assert analysis_type in [
        "gene",
        "variant",
    ], "Name has to be from ['gene', 'variant']"
    assert result_type in [
        "gene",
        "variant",
    ], "Name has to be from ['gene', 'variant']"
    name = 'exome_variant' if analysis_type == 'gene' else 'genome_variant'
    if analysis_type == 'gene' and result_type == 'gene':
        name = 'gene'
    return f"{CURRENT_PATH}/analysis_tables/{analysis_type}_results/{pop.upper()}_{name}_{output_tag}_results.{extension}"
#################
from typing import Union
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

def get_filtered_mt(analysis_type: str, pop: str, filter_samples: bool=True, filter_variants: bool=True, adj_filter: bool=True, prune_samples: bool=False):
    name = 'EXOME MT' if analysis_type=='gene' else 'ACAF MT'
    mt_path = EXOME_MT_PATH if analysis_type=='gene' else ACAF_MT_PATH
    mt = hl.read_matrix_table(mt_path)
    meta_ht = hl.read_table(get_sample_meta_path(annotation=True))
    mt = mt.annotate_cols(**meta_ht[mt.col_key])
    mt = mt.filter_entries(hl.is_missing(mt.FT) | (mt.FT == 'PASS')) # TODO: if using VDS, that's a bool
    # print(mt.aggregate_entries(hl.agg.counter(mt.FT)))
    # {'FAIL': 2323118301, 'PASS': 63243269454, None: 8308706765322}
    if pop is not None and pop != 'all':
        print(f'[{name}] Filtering to {pop.upper()}...')
        mt = mt.filter_cols(mt.pop==pop)
    if filter_variants:
        print(f'[{name}] Filtering to global AC > 0...')
        mt = mt.filter_rows(
            (mt.info.AC[0] > 0)
        )
    if filter_samples:
        print(f'[{name}] Filtering to non-0th degree samples with sex, pop defined, age <= 100...')
        mt = mt.filter_cols(mt.samples_to_keep)

    if prune_samples:
        print(f'[{name}] Filtering to samples pruned from the PCA centroid pipeline...')
        pop_pruned_ht = hl.import_table(CENTROID_PRUNED_SAMPLES, types={'s': hl.tstr}).key_by('s')
        mt = mt.filter_cols(hl.is_defined(pop_pruned_ht[mt.col_key]))
    if adj_filter:
        """
        Filter genotypes to adj criteria - Default: 
        GQ >= 30, het_ref_altAB >= 0.2, het_non_ref_altAB >= 0.2 , het_non_ref_refAB >= 0.2
        """
        print(f'[{name}] Applying adj filters...')
        mt = annotate_adj(mt)
        mt = mt.filter_entries(mt.adj)
        mt = mt.filter_rows(hl.agg.any(mt.GT.n_alt_alleles() >0))
    return mt
