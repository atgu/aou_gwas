import hail as hl

BUCKET = 'gs://aou-analysis/'
TRANCHE = '250k'
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
