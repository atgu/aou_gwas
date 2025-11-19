import hail as hl
from hail.methods.qc import VEPConfigGRCh38Version95, HAIL_GENETICS_VEP_GRCH38_95_IMAGE

dataproc_result = hl.import_table(
    'gs://vep_105/dataproc_vep_grch38_annotations.tsv.gz',
    key=['locus', 'alleles'],
    types={'locus': hl.tlocus('GRCh38'), 'alleles': hl.tarray(hl.tstr), 'vep': hl.tstr},
    force=True,
)
loftee_variants = dataproc_result.select()

hail_vep_result = hl.vep(loftee_variants,  VEPConfigGRCh38Version95(
        data_bucket='hail-qob-vep-grch38-us-central1',
        data_mount='/vep_data/',
        image=HAIL_GENETICS_VEP_GRCH38_95_IMAGE,
        regions=['us-central1'],
        cloud='gcp',
        data_bucket_is_requester_pays=True,
    ))
hail_vep_result = hail_vep_result.annotate(
    vep=hail_vep_result.vep.annotate(
        input=hl.str('\t').join([
            hail_vep_result.locus.contig,
            hl.str(hail_vep_result.locus.position),
            ".",
            hail_vep_result.alleles[0],
            hail_vep_result.alleles[1],
            ".",
            ".",
            "GT",
        ])
    )
)
hail_vep_result = hail_vep_result.select('vep')

def parse_lof_info_into_dict(ht):
    def tuple2(arr):
        return hl.tuple([arr[0], arr[1]])

    return ht.annotate(
        vep=ht.vep.annotate(
            transcript_consequences=ht.vep.transcript_consequences.map(
                lambda csq: csq.annotate(
                    lof_info=hl.or_missing(
                        csq.lof_info != 'null',
                        hl.dict(csq.lof_info.split(',').map(lambda kv: tuple2(kv.split(':')))),
                    )
                )
            )
        )
    )

dataproc_result = dataproc_result.annotate(vep=hl.parse_json(dataproc_result.vep, hail_vep_result.vep.dtype))

hail_vep_result = parse_lof_info_into_dict(hail_vep_result)
dataproc_result = parse_lof_info_into_dict(dataproc_result)

assert hail_vep_result._same(dataproc_result)