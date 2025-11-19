import hail as hl
import argparse
import hailtop.batch as hb

REDO_RAW = False
REDO_PARSED = False
REDO_PARSED_SORTED = False

SORT_STRATEGY = "new-shuffle"
assert SORT_STRATEGY in {"key-by-agg", "new-shuffle", "old-shuffle"}

TAGs = (
    "afr",
    "amr",
    "asj",
    "eas",
    "eur",
    "fin",
    "mid",
    "nfr",
    "sas",
    "oth",
    "max",
    "all",
)
TAGs_gvs = ("afr", "amr", "eas", "eur", "mid", "sas", "oth", "max", "all")
TAGs_gnomad = ("afr", "amr", "asj", "eas", "fin", "nfr", "sas", "oth", "max", "all")
FIELDs = ["ac", "an", "sc"]
tmp_dir = "gs://aou_tmp"
root = "gs://aou_wlu"
vat_root = f"{root}/utils/vat/hail_0_2_107"
vat_path = "gs://prod-drc-broad/v7/wgs/without_ext_aian_prod/aou_srwgs_short_variants_v7_without_ext_aian_prod_20230829_AN_Fixed.vat_complete.bgz.tsv.gz"


def get_vat_int_fields():
    return (
        [f"gvs_{tag}_{field}" for tag in TAGs_gvs for field in FIELDs]
        + [f"gnomad_{tag}_{field}" for tag in TAGs_gnomad for field in FIELDs[:-1]]
        + [
            "position",
            "gene_omim_id",
            "splice_ai_acceptor_gain_distance",
            "splice_ai_acceptor_loss_distance",
            "splice_ai_donor_gain_distance",
            "splice_ai_donor_loss_distance",
        ]
    )


def get_vat_float_fields():
    return (
        [f"gvs_{tag}_af" for tag in TAGs_gvs]
        + [f"gnomad_{tag}_af" for tag in TAGs_gnomad]
        + [
            "revel",
            "splice_ai_acceptor_gain_score",
            "splice_ai_acceptor_loss_score",
            "splice_ai_donor_gain_score",
            "splice_ai_donor_loss_score",
        ]
    )


def get_vat_array_fields():
    return [
        "omim_phenotypes_name",
        "clinvar_classification",
        "clinvar_phenotype",
        "consequence",
        "dbsnp_rsid",
    ]


def parse_empty_missing(v, f):
    new_v = f(v)
    return hl.if_else(hl.len(v) == 0, hl.missing(new_v.dtype), new_v)


def reformat_table():
    if not hl.hadoop_exists(f"{vat_root}/aou_RAW_vat.ht") or REDO_RAW:
        print("---------------Writing RAW VAT----------------------")
        vat_table = hl.import_table(vat_path, quote='"', delimiter="\t", force_bgz=True)
        vat_table.describe()
        vat_table.write(f"{vat_root}/aou_RAW_vat.ht")

    if not hl.hadoop_exists(f"{vat_root}/aou_PARSED_vat.ht") or REDO_PARSED:
        print("---------------Writing PARSED VAT----------------------")
        ht = hl.read_table(f"{vat_root}/aou_RAW_vat.ht")
        ht = ht.annotate(
            **{
                field: parse_empty_missing(ht[field], hl.int32)
                for field in get_vat_int_fields()
            },
            **{
                field: parse_empty_missing(ht[field], hl.float64)
                for field in get_vat_float_fields()
            },
            **{
                field: parse_empty_missing(ht[field], lambda x: x.split(","))
                for field in get_vat_array_fields()
            },
            omim_phenotypes_id=parse_empty_missing(
                ht.omim_phenotypes_id, lambda x: x.split(",")
            ),
            is_canonical_transcript=parse_empty_missing(
                ht.is_canonical_transcript, hl.bool
            ),
        )
        ht.write(f"{vat_root}/aou_PARSED_vat.ht", overwrite=True)

    if (
        not hl.hadoop_exists(f"{vat_root}/aou_PARSED_SORTED_vat.ht")
        or REDO_PARSED_SORTED
    ):
        print("---------------Writing PARSED & SORTED VAT----------------------")
        ht = hl.read_table(f"{vat_root}/aou_PARSED_vat.ht")
        if SORT_STRATEGY == "new-shuffle":
            hl._set_flags(use_new_shuffle="1")
            ht = ht.key_by(
                locus=hl.locus(
                    ht.contig, hl.int32(ht.position), reference_genome="GRCh38"
                ),
                alleles=[ht.ref_allele, ht.alt_allele],
            )
        elif SORT_STRATEGY == "old-shuffle":
            ht = ht.key_by(
                locus=hl.locus(
                    ht.contig, hl.int32(ht.position), reference_genome="GRCh38"
                ),
                alleles=[ht.ref_allele, ht.alt_allele],
            )
        else:
            assert SORT_STRATEGY == "key-by-agg"
            ht = ht.key_by(
                locus=hl.locus(
                    ht.contig, hl.int32(ht.position), reference_genome="GRCh38"
                )
            )
            ht = ht.annotate(alleles=[ht.ref_allele, ht.alt_allele])
            ht = ht.collect_by_key("the_rows")
            ht = ht.annotate(the_rows=hl.sorted(ht.the_rows, lambda x: x.alleles))
            ht = ht.explode("the_rows")
            ht = ht.annotate(**ht.the_rows)
            ht = ht._key_by_assert_sorted("locus", "alleles")

        ht.write(f"{vat_root}/aou_PARSED_SORTED_vat.ht", overwrite=True)


def main():
    hl.init(tmp_dir=tmp_dir)

    reformat_table()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    args = parser.parse_args()
    main(args)
