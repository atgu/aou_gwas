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
RVAS_AF_CUTOFF = 0.001  # below
CALLRATE_CUTOFF = 0.8
N_GENE_PER_GROUP = 10
MIN_CALL_RATE = {
    "all": 0.95,
    "afr": 0.9,
    "amr": 0.85,
    "eas": 0.9,
    "eur": 0.95,
    "mid": 0.55,
    "sas": 0.7,
}
SAIGE_DOCKER_IMAGE = "wzhou88/saige:1.3.0"
PILOT_PHENOTYPES = set(["height", "250.2", "411", "495", "585.3"])

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
TRAIT_TYPEs = ["physical_measurement", "r_drug", "pfhh_survey", "mcc2_phecode"]

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
def get_raw_phenotype_path(
    name: str, parsed: bool = True, extension: str = "ht", tranche: str = TRANCHE
):
    assert name in [
        "r_drug_table",
        "pfhh_survey_table",
        "mcc2_phecode_table",
        "demographics_table",
        "physical_measurement_table",
        "top_10_labs",
    ], 'Name has to be from ["r_drug_table", "pfhh_survey_table", "mcc2_phecode_table", "demographics_table", "physical_measurement_table", "top_10_labs"]'
    return (
        f'{DATA_PATH}/phenotype/{"" if parsed else "raw/"}{name}_{tranche}.{extension}'
    )


def get_random_phenotype_path(pop: str, extension: str = "tsv", test: bool = False):
    return f'{DATA_PATH}/phenotype/random_phenos/random_pheno_{pop}{"_test" if test else ""}.{extension}'


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
    annotation: bool,
    random_pheno: bool = False,
    extension: str = "ht",
    tranche: str = TRANCHE,
):
    return f'{DATA_PATH}/phenotype/aou_all{"_random" if random_pheno else ""}_phenotypes{"_annotated" if annotation else ""}_{tranche}.{extension}'


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
        "sample_qc",
        "sample_qc_tmp",
        "vep_corrected",
        "gene_map",
    ], 'Name has to be from ["ancestry_preds", "relatedness", "relatedness_flagged_samples", "vat", "sample_qc", "sample_qc_tmp", "vep_corrected", "gene_map"]'
    tag1 = (name == "ancestry_preds") and parsed
    tag2 = (name == "gene_map") and parsed
    if name == "vat":
        # gsutil -m cp -r gs://aou_wlu/utils/vat/hail_0_2_107/aou_PARSED_SORTED_vat.ht gs://aou_analysis/250k/data/utils/aou_parsed_and_sorted_vat_hail_0_2_107_250k.ht
        name = "parsed_and_sorted_vat_hail_0_2_107"
    return f'{DATA_PATH}/utils/aou_{name}{"_parsed" if tag1 else ""}{"" if tag1 else "_raw"}_{tranche}.{extension}'


def get_aou_relatedness_path(extension: str = "ht", tranche: str = TRANCHE):
    return f"{DATA_PATH}/utils/aou_ibd_relatedness_{tranche}.{extension}"


##################### Random Phenotype GRM PATHs #####################
def get_aou_sites_for_grm_path(
    pop: str, extension: str, pruned: bool = False, tranche: str = TRANCHE
):
    return f'{DATA_PATH}/utils/grm/aou_{pop}_sites{"_pruned" if pruned else ""}_for_grm_{tranche}.{extension}'


##################### Result PATHs #####################

GENE_INTERVAL_PATH = (
    f"gs://aou_wlu/data/group_positions_{N_GENE_PER_GROUP}_protein_coding.ht"
)


def get_aou_saige_root(analysis_type: str, pop: str, name: str, test: bool = False):
    assert analysis_type in [
        "gene",
        "variant",
    ], "Name has to be from ['gene', 'variant']"
    assert name in [
        "bgen",
        "null_glmm",
        "result",
        "plink",
    ], "Name has to be from ['bgen', 'null_glmm', 'result', 'plink']"
    return f'{CURRENT_PATH}/{analysis_type}_results/{name}/{"test" if test else pop.upper()}'


def get_aou_pheno_file_path(pop: str, phenocode: str, extension: str = "tsv"):
    return f"{CURRENT_PATH}/pheno_file/{pop.upper()}/phenotype_{phenocode}.{extension}"


def get_aou_saige_plink_path(
    analysis_type: str,
    pop: str,
    extension: str,
    pruned: bool,
    data_iteration: int = 0,
    window_size: str = "1e6",
    test: bool = False,
):
    suffix = f"_{data_iteration}" if data_iteration else ""
    cut = "" if window_size == "1e6" else f"_{window_size}"
    pruned = "_pruned" if pruned else ""
    return f'{get_aou_saige_root(analysis_type=analysis_type, pop=pop, name="plink", test=test)}/aou_{pop}_saige_grm{suffix}{pruned}{cut}.{extension}'


def get_aou_samples_file_path(analysis_type: str, pop: str, data_iteration: int = 0):
    suffix = f"_{data_iteration}" if data_iteration else ""
    return f'{get_aou_saige_root(analysis_type=analysis_type, pop=pop, name="plink")}/aou_{pop}{suffix}.samples'
