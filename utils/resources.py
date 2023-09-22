# ALL PATHs in the workbech: https://support.researchallofus.org/hc/en-us/articles/4616869437204-Controlled-CDR-Directory
# Descriptions: https://support.researchallofus.org/hc/en-us/articles/4614687617556-How-the-All-of-Us-Genomic-data-are-organized

##################### ROOT PATHs #####################
TRANCHE = "250k"
BUCKET = "gs://aou_analysis"
MY_BUCKET = "gs://aou_wlu"
TMP_BUCKET = "gs://aou_tmp"
CURRENT_PATH = f"{BUCKET}/{TRANCHE}"
DATA_PATH = f"{CURRENT_PATH}/data"

##################### PARAMETERs #####################
GWAS_AF_CUTOFF = 0.001  # above
RVAS_AF_CUTOFF = 0.0001  # below
CALLRATE_CUTOFF = 0.8
N_GENE_PER_GROUP = 10

##################### POPULATION SUMMARY (Computed in the AoU workbench):
# ANCESTRY_PATH = 'gs://fc-aou-datasets-controlled/v7/wgs/short_read/snpindel/aux/ancestry/ancestry_preds.tsv'
# pop_ht = hl.import_table(ANCESTRY_PATH, key="research_id", impute=True,
#                              types={"research_id": "tstr", "pca_features": hl.tarray(hl.tfloat)})
# pop_ht.aggregate(hl.agg.counter(pop_ht.ancestry_pred))
POPS = ("afr", "amr", "eas", "eur", "mid", "sas")
N_SAMPLES = {
    "all": 245394,
    "afr": 56913,
    "amr": 45035,
    "eas": 5706,
    "eur": 133581,
    "mid": 942,
    "sas": 3217,
}

##################### ORIGINAL PATHs #####################
ORIGINAL_GENO_ROOT = "gs://prod-drc-broad"
ORIGINAL_PHENO_ROOT = "gs://allxall-phenotypes/data"

ORIGINAL_PHENO_PATHs = [
    f"{ORIGINAL_PHENO_ROOT}/top_10_labs.csv",
    f"{ORIGINAL_PHENO_ROOT}/demographics_table.csv",
    f"{ORIGINAL_PHENO_ROOT}/physical_measurement_table.csv",
    f"{ORIGINAL_PHENO_ROOT}/r_drug_table.csv",
    f"{ORIGINAL_PHENO_ROOT}/pfhh_survey_table.csv",
    f"{ORIGINAL_PHENO_ROOT}/mcc2_phecode_table.csv",
]
TRAIT_TYPEs = ['physical_measurement', 'r_drug', 'pfhh_survey', 'mcc2_phecode']

ORIGINAL_CALLSET_PATH = f"{ORIGINAL_GENO_ROOT}/aou-wgs-delta-small_callsets_gq0/v7.1"
ACAF_MT_PATH = f"{ORIGINAL_CALLSET_PATH}/acaf_threshold_v7.1/splitMT/delta_basis_without_ext_aian_prod_gq0_3regions.acaf_threshold.split.mt"
EXOME_MT_PATH = f"{ORIGINAL_CALLSET_PATH}/exome_v7.1/splitMT/delta_basis_without_ext_aian_prod_gq0_3regions.exome.split.mt"
VDS_PATH = f"{ORIGINAL_GENO_ROOT}/v7/wgs/without_ext_aian_prod/vds/aou_srwgs_short_variants_v7_without_ext_aian_prod.vds"

ORIGINAL_UTIL_PATH = f"{ORIGINAL_GENO_ROOT}/aou-wgs-delta-aux_gq0"
ORIGINAL_POP_PATH = (
    f"{ORIGINAL_UTIL_PATH}/ancestry/delta_v1_gt_no_ext_aian_gq0_prod.ancestry_preds.tsv"
)
ORIGINAL_RELATED_INFO_PATH = (
    f"{ORIGINAL_UTIL_PATH}/relatedness/delta_v1_gt_no_ext_aian_gq0_prod.relatedness.tsv"
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


##################### PHENOTYPE & SAMPLE PATHs #####################
def get_raw_phenotype_path(name: str, parsed: bool = True, tranche: str = TRANCHE):
    assert name in [
        "r_drug_table",
        "pfhh_survey_table",
        "mcc2_phecode_table",
        "demographics_table",
        "physical_measurement_table",
        "top_10_labs",
    ], 'Name has to be from ["r_drug_table", "pfhh_survey_table", "mcc2_phecode_table", "demographics_table", "physical_measurement_table", "top_10_labs"]'
    return f'{DATA_PATH}/phenotype/{"" if parsed else "raw/"}{name}_{tranche}.ht'


BINARY_HT_PATHs = [
    get_raw_phenotype_path(name="r_drug_table", parsed=True),
    get_raw_phenotype_path(name="pfhh_survey_table", parsed=True),
    get_raw_phenotype_path(name="mcc2_phecode_table", parsed=True),
]
DEMOGRAPHICS_HT_PATH = get_raw_phenotype_path(name="demographics_table", parsed=True)
QUANTITATIVE_HT_PATH = get_raw_phenotype_path(
    name="physical_measurement_table", parsed=True
)


def get_full_phenotype_path(
    annotation: bool, extension: str = "ht", tranche: str = TRANCHE
):
    return f'{DATA_PATH}/phenotype/aou_all_phenotypes{"_annotated" if annotation else ""}_{tranche}.{extension}'


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
    ], 'Name has to be from ["ancestry_preds", "relatedness", "relatedness_flagged_samples", "vat"]'
    tag = (name == "ancestry_preds") and parsed
    if name == "vat":
        # gsutil -m cp -r gs://aou_wlu/utils/vat/hail_0_2_107/aou_PARSED_SORTED_vat.ht gs://aou_analysis/250k/data/utils/aou_parsed_and_sorted_vat_hail_0_2_107_250k.ht
        name = "parsed_and_sorted_vat_hail_0_2_107"
    return (
        f'{DATA_PATH}/utils/aou_{name}{"_parsed" if tag else ""}_{tranche}.{extension}'
    )


##### RANDOM PHENOTYPE PATHs
def get_aou_sites_for_grm_path(
    pop: str, extension: str, pruned: bool = False, tranche: str = TRANCHE
):
    return f'{DATA_PATH}/utils/aou_{pop}_sites{"_pruned" if pruned else ""}_for_grm_{tranche}.{extension}'


def get_aou_relatedness_path(extension: str = "ht", tranche: str = TRANCHE):
    return f"{DATA_PATH}/utils/aou_ibd_relatedness_{tranche}.{extension}"


#
# def get_aou_grm_mt_path(pop: str, data_iteration:int, tranche: str = TRANCHE):
#     suffix = f'_{data_iteration}' if data_iteration else ""
#     return f'{CURRENT_BUCKET}/grm/aou_{pop}_for_grm{suffix}_{tranche}.mt'
#
# def get_aou_grm_pruned_ht_path(pop: str, window_size: str = '1e6', tranche: str = TRANCHE):
#     cut = '' if window_size == '1e6' else f'_{window_size}'
#     return f'{CURRENT_BUCKET}/grm/aou_{pop}_for_grm_pruned{cut}_{tranche}.ht'
#
# def get_aou_grm_plink_path(pop: str, data_iteration: int = 0, window_size: str = '1e6', tranche: str = TRANCHE):
#     suffix = f'.{data_iteration}' if data_iteration else ""
#     cut = '' if window_size == '1e6' else f'.{window_size}'
#     return f'{CURRENT_BUCKET}/grm/aou_{pop}_for_grm{suffix}_pruned{cut}_{tranche}.plink'
