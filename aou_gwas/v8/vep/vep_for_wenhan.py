from typing import List, Optional

import hail as hl
from hail.expr.types import tarray, tfloat, tint32, tstr, tstruct
from hail.genetics import reference_genome
from hail.methods.qc import VEPConfig

TRANCHE = "v8"
ANALYSIS_BUCKET = "gs://aou_analysis"
MY_BUCKET = 'gs://aou_wlu'
TMP_BUCKET = 'gs://aou_tmp'
DATA_PATH = f"{ANALYSIS_BUCKET}/{TRANCHE}/data"
ORIGINAL_DATA_ROOT = 'gs://allxall-phenotypes-v2'
ORIGINAL_PHENO_DIR = f'{ORIGINAL_DATA_ROOT}/2025-04-14'
ORIGINAL_VAT_PATH = f'{ORIGINAL_DATA_ROOT}/2025-04-17/vat'

hl.init(app_name="vep_test", gcs_requester_pays_configuration="aou-neale-gwas")
# ht = hl.read_table('gs://aou_analysis/250k/data/utils/aou_vat_collected_by_key_indel_250k.ht')
# ht = ht.filter((ht.locus.contig == 'chr1') & (ht.locus.position < 10330))
# # ht = ht.filter((ht.locus.contig == 'chr21'))

# ht = hl.Table.parallelize([
#     {'locus': hl.locus('chr11', 113409605, reference_genome='GRCh38'), 'alleles': ['T', 'A']},
#     {'locus': hl.locus('chr11', 113409605, reference_genome='GRCh38'), 'alleles': ['T', 'C']},
#     {'locus': hl.locus('chr11', 113409605, reference_genome='GRCh38'), 'alleles': ['T', 'G']},
#     {'locus': hl.locus('chr11', 113409606, reference_genome='GRCh38'), 'alleles': ['C', 'A']},
#     {'locus': hl.locus('chr11', 113409606, reference_genome='GRCh38'), 'alleles': ['C', 'G']},
#     {'locus': hl.locus('chr11', 113409606, reference_genome='GRCh38'), 'alleles': ['C', 'T']},
#     {'locus': hl.locus('chr11', 113409607, reference_genome='GRCh38'), 'alleles': ['C', 'T']},
#     {'locus': hl.locus('chr11', 113409607, reference_genome='GRCh38'), 'alleles': ['C', 'A']},
#     {'locus': hl.locus('chr11', 113409607, reference_genome='GRCh38'), 'alleles': ['C', 'G']},
#     {'locus': hl.locus('chr11', 113409608, reference_genome='GRCh38'), 'alleles': ['T', 'A']}],
#     hl.tstruct(locus=hl.tlocus('GRCh38'), alleles=hl.tarray(hl.tstr)),
#     key=['locus', 'alleles'])
# ht.show()

##### 2025.4.18 export indels in VAT
# ht = hl.read_table(f'{ORIGINAL_VAT_PATH}/aou_PARSED_SORTED_COLLECTED_vat_wlu.ht')
# indel_ht = ht.filter(~hl.is_snp(ht.alleles[0], ht.alleles[1]))
# indel_ht = indel_ht.naive_coalesce(500).checkpoint(f'{DATA_PATH}/vat/aou_indel_vat_v8.ht', overwrite=True)
# indel_ht.describe()
# indel_ht.show()
# print(indel_ht.count())
# indel_ht = indel_ht.filter((indel_ht.locus.contig == 'chr1') & (indel_ht.locus.position < 10330))

indel_ht = hl.read_table(f'{DATA_PATH}/vat/aou_indel_vat_v8.ht')
# indel_ht = indel_ht.filter((indel_ht.locus.contig == 'chr1') & (indel_ht.locus.position < 10330))
indel_ht = indel_ht.filter((indel_ht.locus.contig == 'chr21'))

MY_IMAGE_WITH_VEP_105_AND_HAIL_PYTHON_VEP_SCRIPT = 'us-central1-docker.pkg.dev/aou-neale-gwas/vep/vep_105_wenhan:v2.7'

vep_json_typ = tstruct(
    assembly_name=tstr,
    allele_string=tstr,
    ancestral=tstr,
    colocated_variants=tarray(
        tstruct(
            aa_allele=tstr,
            aa_maf=tfloat,
            afr_allele=tstr,
            afr_maf=tfloat,
            allele_string=tstr,
            amr_allele=tstr,
            amr_maf=tfloat,
            clin_sig=tarray(tstr),
            end=tint32,
            eas_allele=tstr,
            eas_maf=tfloat,
            ea_allele=tstr,
            ea_maf=tfloat,
            eur_allele=tstr,
            eur_maf=tfloat,
            exac_adj_allele=tstr,
            exac_adj_maf=tfloat,
            exac_allele=tstr,
            exac_afr_allele=tstr,
            exac_afr_maf=tfloat,
            exac_amr_allele=tstr,
            exac_amr_maf=tfloat,
            exac_eas_allele=tstr,
            exac_eas_maf=tfloat,
            exac_fin_allele=tstr,
            exac_fin_maf=tfloat,
            exac_maf=tfloat,
            exac_nfe_allele=tstr,
            exac_nfe_maf=tfloat,
            exac_oth_allele=tstr,
            exac_oth_maf=tfloat,
            exac_sas_allele=tstr,
            exac_sas_maf=tfloat,
            id=tstr,
            minor_allele=tstr,
            minor_allele_freq=tfloat,
            phenotype_or_disease=tint32,
            pubmed=tarray(tint32),
            sas_allele=tstr,
            sas_maf=tfloat,
            somatic=tint32,
            start=tint32,
            strand=tint32,
        )
    ),
    context=tstr,
    end=tint32,
    id=tstr,
    input=tstr,
    intergenic_consequences=tarray(
        tstruct(
            allele_num=tint32,
            consequence_terms=tarray(tstr),
            impact=tstr,
            minimised=tint32,
            variant_allele=tstr,
        )
    ),
    most_severe_consequence=tstr,
    motif_feature_consequences=tarray(
        tstruct(
            allele_num=tint32,
            consequence_terms=tarray(tstr),
            high_inf_pos=tstr,
            impact=tstr,
            minimised=tint32,
            motif_feature_id=tstr,
            motif_name=tstr,
            motif_pos=tint32,
            motif_score_change=tfloat,
            strand=tint32,
            variant_allele=tstr,
        )
    ),
    regulatory_feature_consequences=tarray(
        tstruct(
            allele_num=tint32,
            biotype=tstr,
            consequence_terms=tarray(tstr),
            impact=tstr,
            minimised=tint32,
            regulatory_feature_id=tstr,
            variant_allele=tstr,
        )
    ),
    seq_region_name=tstr,
    start=tint32,
    strand=tint32,
    transcript_consequences=tarray(
        tstruct(
            allele_num=tint32,
            amino_acids=tstr,
            biotype=tstr,
            canonical=tint32,
            ccds=tstr,
            cdna_start=tint32,
            cdna_end=tint32,
            cds_end=tint32,
            cds_start=tint32,
            codons=tstr,
            consequence_terms=tarray(tstr),
            distance=tint32,
            domains=tarray(tstruct(db=tstr, name=tstr)),
            exon=tstr,
            gene_id=tstr,
            gene_pheno=tint32,
            gene_symbol=tstr,
            gene_symbol_source=tstr,
            hgnc_id=tstr,
            hgvsc=tstr,
            hgvsp=tstr,
            hgvs_offset=tint32,
            impact=tstr,
            intron=tstr,
            lof=tstr,
            lof_flags=tstr,
            lof_filter=tstr,
            lof_info=tstr,
            minimised=tint32,
            polyphen_prediction=tstr,
            polyphen_score=tfloat,
            protein_end=tint32,
            protein_start=tint32,
            protein_id=tstr,
            sift_prediction=tstr,
            sift_score=tfloat,
            strand=tint32,
            swissprot=tstr,
            transcript_id=tstr,
            trembl=tstr,
            uniparc=tstr,
            variant_allele=tstr,
        )
    ),
    variant_class=tstr,
)


class VEPConfigGRCh38Version105(VEPConfig):
    """
    The Hail-maintained VEP configuration for GRCh38 for VEP version 95.

    This class takes the following constructor arguments:

     - `data_bucket` (:obj:`.str`) -- The location where the VEP data is stored.
     - `data_mount` (:obj:`.str`) -- The location in the container where the data should be mounted.
     - `image` (:obj:`.str`) -- The docker image to run VEP.
     - `cloud` (:obj:`.str`) -- The cloud where the Batch Service is located.
     - `data_bucket_is_requester_pays` (:obj:`.bool`) -- True if the data bucket is set to requester pays.
     - `regions` (:obj:`.list` of :obj:`.str`) -- A list of regions the VEP jobs can run in.

    """

    def __init__(
        self,
        *,
        data_bucket: str,
        data_mount: str,
        image: str,
        regions: List[str],
        cloud: str,
        data_bucket_is_requester_pays: bool,
    ):
        self.data_bucket = data_bucket
        self.data_mount = data_mount
        self.image = image
        self.regions = regions
        self.env = {}
        self.data_bucket_is_requester_pays = data_bucket_is_requester_pays
        self.cloud = cloud
        self.batch_run_command = ['python3', '/hail-vep/vep.py', 'vep']
        self.batch_run_csq_header_command = ['python3', '/hail-vep/vep.py', 'csq_header']

        # You might need to change the type. See the schema I sent you previously to know what you might need to change
        # there's also "_drop_fields('not-needed-field-name')"
        self.json_typ = vep_json_typ._insert_field(
            'transcript_consequences',
            tarray(
                vep_json_typ['transcript_consequences'].element_type._insert_fields(
                    appris=tstr,
                    tsl=tint32,
                )
            ),
        )
    def command(
        self,
        *,
        consequence: bool,
        tolerate_parse_error: bool,
        part_id: int,
        input_file: Optional[str],
        output_file: str,
    ) -> str:
        vcf_or_json = '--vcf' if consequence else '--json'
        input_file = f'--input_file {input_file}' if input_file else ''

        return f"""bash -lc ' \
cp -r {self.data_mount} /tmp/vep && \
find /tmp/vep -type f -name "*.csi" -exec touch {{}} + && \
/vep {input_file} \
--format vcf \
{vcf_or_json} \
--everything \
--allele_number \
--no_stats \
--cache \
--dir_cache {self.data_mount} \
--offline \
--minimal \
--assembly GRCh38 \
--fasta {self.data_mount}homo_sapiens/105_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa.gz \
--plugin "LoF,loftee_path:/opt/vep/Plugins/,gerp_bigwig:{self.data_mount}/gerp_conservation_scores.homo_sapiens.GRCh38.bw,\
human_ancestor_fa:{self.data_mount}/human_ancestor.fa.gz,conservation_file:{self.data_mount}/loftee.sql" \
--dir_plugins /opt/vep/Plugins/ \
-o STDOUT'"""



#     def command(
#         self,
#         *,
#         consequence: bool,
#         tolerate_parse_error: bool,
#         part_id: int,
#         input_file: Optional[str],
#         output_file: str,
#     ) -> str:
#         vcf_or_json = '--vcf' if consequence else '--json'
#         input_file = f'--input_file {input_file}' if input_file else ''

#         return f"""/vep {input_file} \
# --format vcf \
# {vcf_or_json} \
# --everything \
# --allele_number \
# --no_stats \
# --cache \
# --dir_cache {self.data_mount} \
# --offline \
# --minimal \
# --assembly GRCh38 \
# --fasta {self.data_mount}homo_sapiens/105_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa.gz \
# --plugin "LoF,loftee_path:/opt/vep/Plugins/,gerp_bigwig:{self.data_mount}/gerp_conservation_scores.homo_sapiens.GRCh38.bw,human_ancestor_fa:{self.data_mount}/human_ancestor.fa.gz,conservation_file:{self.data_mount}/loftee.sql" \
# --dir_plugins /opt/vep/Plugins/ \
# -o STDOUT
# """


my_config_105 = VEPConfigGRCh38Version105(
        data_bucket='vep_105',
        data_mount='/opt/vep/.vep/',
        image=MY_IMAGE_WITH_VEP_105_AND_HAIL_PYTHON_VEP_SCRIPT,
        regions=['us-central1'],
        cloud='gcp',
        data_bucket_is_requester_pays=True,
    )
vep_ht = hl.vep(indel_ht, config=my_config_105)
vep_ht.describe()
# vep_ht = vep_ht.checkpoint(f'{DATA_PATH}/vep/aou_vep_chr1_indel_test.ht', overwrite=True)
# vep_ht.show()


# def test_vep_grch38_against_dataproc():
#     dataproc_result = hl.import_table(
#         'gs://vep_105/dataproc_vep_grch38_annotations.tsv.gz',
#         key=['locus', 'alleles'],
#         types={'locus': hl.tlocus('GRCh38'), 'alleles': hl.tarray(hl.tstr), 'vep': hl.tstr},
#         force=True,
#     )
#     loftee_variants = dataproc_result.select()

#     hail_vep_result = hl.vep(loftee_variants, my_config_105)
#     hail_vep_result = hail_vep_result.annotate(
#         vep=hail_vep_result.vep.annotate(
#             input=hl.str('\t').join([
#                 hail_vep_result.locus.contig,
#                 hl.str(hail_vep_result.locus.position),
#                 ".",
#                 hail_vep_result.alleles[0],
#                 hail_vep_result.alleles[1],
#                 ".",
#                 ".",
#                 "GT",
#             ])
#         )
#     )
#     hail_vep_result = hail_vep_result.select('vep')

#     def parse_lof_info_into_dict(ht):
#         def tuple2(arr):
#             return hl.tuple([arr[0], arr[1]])

#         return ht.annotate(
#             vep=ht.vep.annotate(
#                 transcript_consequences=ht.vep.transcript_consequences.map(
#                     lambda csq: csq.annotate(
#                         lof_info=hl.or_missing(
#                             csq.lof_info != 'null',
#                             hl.dict(csq.lof_info.split(',').map(lambda kv: tuple2(kv.split(':')))),
#                         )
#                     )
#                 )
#             )
#         )

#     dataproc_result = dataproc_result.annotate(vep=hl.parse_json(dataproc_result.vep, hail_vep_result.vep.dtype))

#     hail_vep_result = parse_lof_info_into_dict(hail_vep_result)
#     dataproc_result = parse_lof_info_into_dict(dataproc_result)

#     assert hail_vep_result._same(dataproc_result)


# test_vep_grch38_against_dataproc()