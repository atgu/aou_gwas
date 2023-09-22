#!/usr/bin/env python3

__author__ = "wlu"

from utils.utils import *
from utils.resources import *
import hailtop.batch as hb
import argparse
import hail as hl


def update_sample_util_ht(original_path, overwrite):
    name = original_path.split(".")[-2]
    raw_path = get_aou_util_path(name=name, parsed=False, extension="ht")
    print(raw_path)

    if not hl.hadoop_exists(raw_path):
        ht = hl.import_table(original_path, quote='"', delimiter="\t", impute=True)
        ht = ht.key_by(*{ht[i] for i in table_keys[name]})
        ht.naive_coalesce(1000).checkpoint(
            raw_path, _read_if_exists=not overwrite, overwrite=overwrite
        )
        ht = hl.read_table(raw_path)
        ht.describe()

    if name == "ancestry_preds":
        parsed_path = get_aou_util_path(name=name, parsed=True, extension="ht")
        if not hl.hadoop_exists(parsed_path) or overwrite:
            ht = hl.read_table(raw_path)
            ht = ht.annotate(
                **{
                    field: hl.map(
                        lambda x: hl.float64(x),
                        parse_empty_missing(
                            ht[field],
                            lambda x: x.replace("\\[", "")
                            .replace("\\]", "")
                            .split(","),
                        ),
                    )
                    for field in ["probabilities", "pca_features"]
                }
            )
            ht.naive_coalesce(1000).checkpoint(
                parsed_path, _read_if_exists=not overwrite, overwrite=overwrite
            )
            ht = hl.read_table(parsed_path)
            ht.describe()
    return ht


def write_phenotype_ht(original_path=ORIGINAL_PHENO_PATHs, overwrite: bool = False):
    def parse_empty_missing(v, f):
        new_v = f(v)
        return hl.if_else((hl.len(v) == 0), hl.missing(new_v.dtype), new_v)

    # paths = list(path['path'] for path in hl.hadoop_ls(original_path))
    for path in original_path:
        name = path.split("/")[-1].split(".")[0]
        if not hl.hadoop_exists(get_raw_phenotype_path(name=name, parsed=False)):
            ht = hl.import_table(path, delimiter=",", missing=".")
            print(f"--------{path}--------------")
            if "person_id" in ht.row_value:
                ht = ht.annotate(person_id=hl.int32(ht.person_id))
                ht = ht.key_by("person_id")
            else:
                ht = ht.annotate(
                    measurement_concept_id=hl.int32(ht.measurement_concept_id)
                )
                ht = ht.key_by("measurement_concept_id")
            ht.checkpoint(
                get_raw_phenotype_path(name=name, parsed=False), _read_if_exists=True
            )

        if (
            not hl.hadoop_exists(get_raw_phenotype_path(name=name, parsed=True))
        ) or overwrite:
            impute = False
            if name in ["r_drug_table", "pfhh_survey_table", "mcc2_phecode_table"]:
                parse = hl.bool
            if name == "physical_measurement_table":
                parse = hl.float64
            if name in ["top_10_labs", "demographics_table"]:
                impute = True
                parse = None
            if impute:
                ht = hl.import_table(path, delimiter=",", missing=".", impute=impute)
                if "person_id" in ht.row_value:
                    ht = ht.key_by("person_id")
                else:
                    ht = ht.key_by("measurement_concept_id")
            if parse is not None:
                ht = hl.read_table(get_raw_phenotype_path(name=name, parsed=False))
                ht = ht.annotate(
                    **{
                        field: hl.if_else(ht[field] == "NA", "", ht[field])
                        for field in list(ht.row_value)
                    },
                )
                ht = ht.annotate(
                    **{
                        field: parse_empty_missing(ht[field], parse)
                        for field in list(ht.row_value)
                    },
                )
            ht = ht.repartition(1000)
            ht.checkpoint(
                get_raw_phenotype_path(name=name, parsed=True), overwrite=overwrite
            )
        else:
            ht = hl.read_table(get_raw_phenotype_path(name=name, parsed=True))
            ht.describe()
        return ht


def merge_phenotype_hts(overwrite: bool):
    if not hl.hadoop_exists(get_full_phenotype_path(annotation=False)) or overwrite:
        ht = hl.read_table(QUANTITATIVE_HT_PATH)
        for path in BINARY_HT_PATHs:
            print(path)
            tmp_ht = hl.read_table(path)
            ht = ht.join(tmp_ht, how="outer")
        ht.write(get_full_phenotype_path(annotation=False), overwrite)
    else:
        # Table info: https://batch.hail.is/batches/8041785/jobs/1
        ht = hl.read_table(get_full_phenotype_path(annotation=False))
        print(get_full_phenotype_path(annotation=False))
        print(f"N Phenotypes: {len(list(ht.row_value))}")
    return ht


def update_meta_ht(overwrite: bool):
    if (
        not hl.hadoop_exists(get_sample_meta_path(annotation=True, extension="ht"))
        or overwrite
    ):
        meta_ht = hl.read_table(DEMOGRAPHICS_HT_PATH)  # Demographic table
        meta_ht.checkpoint(
            get_sample_meta_path(annotation=False, extension="ht"), _read_if_exists=True
        )
        pop_ht = hl.read_table(get_aou_util_path(name="ancestry_preds"))
        related_ht = hl.read_table(
            get_aou_util_path(name="relatedness_flagged_samples")
        )
        meta_ht = meta_ht.annotate(
            related_sample=hl.is_defined(related_ht[meta_ht.key]),
            sex=hl.if_else(
                meta_ht.sex_at_birth == "Male",
                0,
                hl.if_else(meta_ht.sex_at_birth == "Female", 1, hl.missing(hl.tint32)),
            ),  # TODO: confirm this
            age=meta_ht.age_at_cdr,  # TODO: 8 people with age over 115
            age2=meta_ht.age_at_cdr**2,
        )
        meta_ht = meta_ht.annotate(
            age_sex=meta_ht.age * meta_ht.sex, age2_sex=meta_ht.age2 * meta_ht.sex
        )
        meta_ht = meta_ht.annotate(**pop_ht[meta_ht.key])
        meta_ht = meta_ht.annotate(
            **{f"PC{x + 1}": meta_ht.pca_features[x] for x in range(16)}
        )

        meta_ht.naive_coalesce(1000).write(
            get_sample_meta_path(annotation=True, extension="ht"), overwrite
        )
    else:
        meta_ht = hl.read_table(get_sample_meta_path(annotation=True, extension="ht"))
        meta_ht.describe()
    return meta_ht


def annotate_phenotype_ht(overwrite: bool):
    if not hl.hadoop_exists(get_full_phenotype_path(annotation=True)) or overwrite:
        global_info = {}
        for trait_type in TRAIT_TYPEs:
            tmp_ht = hl.read_table(get_raw_phenotype_path(name=f"{trait_type}_table"))
            global_info[trait_type] = list(tmp_ht.row_value)
        print(global_info)
        ht = hl.read_table(get_full_phenotype_path(annotation=False))
        meta_ht = hl.read_table(get_sample_meta_path(annotation=True, extension="ht"))
        ht = ht.annotate(**meta_ht[ht.key])
        ht = ht.annotate_globals(trait_type=global_info)
        ht.write(get_full_phenotype_path(annotation=True), overwrite)
    else:
        ht = hl.read_table(get_full_phenotype_path(annotation=True))
        ht.describe()
    return ht


def update_phenotype_mt(overwrite):
    ht = hl.read_table(get_full_phenotype_path(annotation=True))
    quant_ht = hl.read_table(QUANTITATIVE_HT_PATH)
    mt = ht.to_matrix_table_row_major(
        columns=list(ht.row_value), entry_field_name="value", col_field_name="phenocode"
    )
    mt = mt.annotate_cols(
        trait_type=hl.if_else(
            hl.literal(list(quant_ht.row_value)).contains(mt.phenocode),
            "countinuous",
            "categorical",
        )
    )
    mt.describe()
    mt.write(get_full_phenotype_path(annotation=True, extension="mt"), overwrite)


def main(args):
    try:
        hl.init(default_reference="GRCh38", tmp_dir=TMP_BUCKET)

        if args.batch:
            backend = hb.ServiceBackend(
                billing_project="all-by-aou",
                remote_tmpdir=TMP_BUCKET,
            )
            b = hb.Batch(
                name=f"AoU_sample_info_preprocess",
                requester_pays_project="aou-neale-gwas",
                default_python_image="hailgenetics/hail:0.2.123-py3.9",
                backend=backend,
            )

        if args.update_sample_util_ht:
            for path in SAMPLE_INFO_PATHs:
                update_sample_util_ht(path, args.overwrite)

        if args.update_meta_ht:
            if args.batch:
                j = b.new_python_job(name=f"Update_Sample_Meta_HT")
                j.call(update_meta_ht, args.overwrite)
            else:
                update_meta_ht(overwrite=args.overwrite)

        if args.update_raw_phenotypes:
            print("-------------Writing Raw Phenotype HTs-----------")
            write_phenotype_ht(
                original_path=ORIGINAL_PHENO_PATHs, overwrite=args.overwrite
            )

        if args.update_phenotype_ht:
            print("-------------Merging Phenotype HTs to a Main HT-----------")
            merge_phenotype_hts(overwrite=args.overwrite)

        if args.annotate_phenotype_ht:
            print("-------------Annotating Main Phenotype HT-----------")
            if args.batch:
                j = b.new_python_job(name=f"Annotate_phenotype_HT")
                j.call(annotate_phenotype_ht, args.overwrite)
            else:
                annotate_phenotype_ht(overwrite=args.overwrite)

        if args.update_phenotype_mt:  # Not Necessary
            print("-------------Converting Main Phenotype HT to an MT-----------")
            if args.batch:
                j = b.new_python_job(name=f"Create_phenotype_MT")
                j.call(update_phenotype_mt, args.overwrite)
            else:
                update_phenotype_mt(overwrite=args.overwrite)

        if args.batch:
            b.run()

    finally:
        if not args.batch:
            from datetime import date

            hl.copy_log(f"{MY_BUCKET}/log/preprocessing_{date.today()}.log")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--update-sample-util-ht",
        help="Update the sample information tables",
        action="store_true",
    )
    parser.add_argument(
        "--update-raw-phenotypes",
        help="Update the raw separate phenotype tables",
        action="store_true",
    )
    parser.add_argument(
        "--update-meta-ht",
        help="Update the sample meta info Hail table",
        action="store_true",
    )
    parser.add_argument(
        "--update-phenotype-ht",
        help="Update the full phenotype Hail table",
        action="store_true",
    )
    parser.add_argument(
        "--annotate-phenotype-ht",
        help="Annotate the full phenotype Hail table",
        action="store_true",
    )
    parser.add_argument(
        "--update-phenotype-mt",
        help="Update the processed phenotype Hail Matrixtable",
        action="store_true",
    )
    parser.add_argument(
        "--batch",
        help="Run the pipeline in Batch, not QoB or dataproc",
        action="store_true",
    )
    parser.add_argument(
        "--overwrite",
        help="Whether to overwrite the existing table",
        action="store_true",
    )
    args = parser.parse_args()

    main(args)
