#!/usr/bin/env python3

__author__ = "wlu"


# from aou_gwas import *  # for dataproc
from utils.utils import *  # for QoB
from utils.resources import *  # for QoB
import hailtop.batch as hb
import argparse
import hail as hl


def update_sample_util_ht(original_path, overwrite):
    name = original_path.split(".")[-2]
    raw_path = get_aou_util_path(name=name, parsed=False, extension="ht")
    print(raw_path)

    if not hl.hadoop_exists(raw_path) or overwrite:
        print(f"-----------Loading Raw {name} Table-----------")
        types = {k: hl.tstr for k in table_keys[name]}
        if name == "relatedness_flagged_samples":
            ht = hl.import_table(original_path, quote='"', delimiter="\t")
        else:
            ht = hl.import_table(
                original_path, quote='"', delimiter="\t", impute=True, types=types
            )
        ht = ht.key_by(*{ht[i] for i in table_keys[name]})
        ht.describe()
        ht.naive_coalesce(1000).checkpoint(
            raw_path, _read_if_exists=not overwrite, overwrite=overwrite
        )
        ht = hl.read_table(raw_path)
        ht.describe()

    if name == "ancestry_preds":
        print(f"-----------Parsing Raw {name} Table-----------")
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
    # paths = list(path['path'] for path in hl.hadoop_ls(original_path))
    # NOTE: Cannot run hl.hadoop_ls on QoB
    for path in original_path:
        name = path.split("/")[-1].split(".")[0]
        if not hl.hadoop_exists(get_raw_phenotype_path(name=name, parsed=False)):
            ht = hl.import_table(path, delimiter=",", missing=".")
            print(f"--------{path}--------------")
            if "person_id" in ht.row_value:
                ht = ht.key_by("person_id")
            else:
                ht = ht.key_by("measurement_concept_id")
            ht.checkpoint(
                get_raw_phenotype_path(name=name, parsed=False),
                _read_if_exists=not overwrite,
                overwrite=overwrite,
            )
            ht.describe()

        if (
            not hl.hadoop_exists(get_raw_phenotype_path(name=name, parsed=True))
        ) or overwrite:
            print(
                f"-------- Parsing {name} to {get_raw_phenotype_path(name=name, parsed=True)}--------------"
            )
            if name in ["r_drug_table", "pfhh_survey_table", "mcc2_phecode_table"]:
                ht = hl.read_table(get_raw_phenotype_path(name=name, parsed=False))
                # def _parse_bool(v):
                #     return hl.if_else((v == "NA") | (hl.len(v) == 0), hl.missing(hl.tbool), hl.bool(v))
                #
                # parse_bool = hl.experimental.define_function(_parse_bool, hl.tstr)
                #
                # ht = ht.annotate(
                #     **{
                #         field: parse_bool(ht[field])
                #         for field in list(ht.row_value)
                #     },
                # )

                ht = ht.annotate(
                    **{
                        field: hl.if_else(ht[field] == "NA", "", ht[field])
                        for field in list(ht.row_value)
                    },
                )
                ht = ht.annotate(
                    **{
                        field: parse_empty_missing(ht[field], hl.bool)
                        for field in list(ht.row_value)
                    },
                )  # These works with dataproc

            if name == "physical_measurement_table":
                ht = hl.read_table(get_raw_phenotype_path(name=name, parsed=False))

                def _parse_float(v):
                    return hl.if_else(
                        (v == "NA") | (hl.len(v) == 0),
                        hl.missing(hl.tfloat64),
                        hl.float64(v),
                    )

                parse_float = hl.experimental.define_function(_parse_float, hl.tstr)

                ht = ht.annotate(
                    **{field: parse_float(ht[field]) for field in list(ht.row_value)},
                )

            if name in ["top_10_labs", "demographics_table"]:
                if name == "demographics_table":
                    ht = hl.import_table(
                        path,
                        delimiter=",",
                        missing=".",
                        impute=True,
                        types={"person_id": hl.tstr},
                    )
                    ht = ht.key_by("person_id")
                else:
                    ht = hl.import_table(
                        path,
                        delimiter=",",
                        missing=".",
                        impute=True,
                        types={"measurement_concept_id": hl.tstr},
                    )
                    ht = ht.key_by("measurement_concept_id")

            ht = ht.repartition(1000)
            ht.checkpoint(
                get_raw_phenotype_path(name=name, parsed=True), overwrite=overwrite
            )
        print(name)
        ht = hl.read_table(get_raw_phenotype_path(name=name, parsed=True))
        ht.describe()
        phenos = list(ht.row_value)

        with hl.hadoop_open(
            get_raw_phenotype_path(name=name, parsed=True, extension="txt"), "w"
        ) as f:
            for pheno in phenos:
                f.write(pheno + "\n")
    return ht


def merge_phenotype_hts(overwrite: bool):
    if not hl.hadoop_exists(get_full_phenotype_path(annotation=False)) or overwrite:
        # TODO: add WHRadjBMI_custom (how to regress out BMI in hail) & sex-specific phenotypes & upcoming phenotypes
        # TODO: update new quantitative HT
        ht = hl.read_table(QUANTITATIVE_HT_PATH)
        ht = ht.annotate(
            BMI_custom=ht.weight / ((ht.height / 100) ** 2),
            WHR_custom=ht["waist-circumference-mean"] / ht["hip-circumference-mean"],
        )
        ht.summarize()
        for path in BINARY_HT_PATHs:
            print(path)
            tmp_ht = hl.read_table(path)
            ht = ht.join(tmp_ht, how="outer")
        ht.write(get_full_phenotype_path(annotation=False), overwrite)
    # Table info: https://batch.hail.is/batches/8041785/jobs/1
    ht = hl.read_table(get_full_phenotype_path(annotation=False))
    print(list(ht.row_value))
    print(f"N Phenotypes: {len(list(ht.row_value))}")
    return ht


def update_meta_ht(overwrite: bool):
    if (
        not hl.hadoop_exists(get_sample_meta_path(annotation=True, extension="ht"))
        or overwrite
    ):
        meta_ht = hl.read_table(DEMOGRAPHICS_HT_PATH)  # Demographic table
        meta_ht.checkpoint(
            get_sample_meta_path(annotation=False, extension="ht"),
            _read_if_exists=True,
        )
        pop_ht = hl.read_table(get_aou_util_path(name="ancestry_preds", parsed=True))
        related_0th_ht = hl.read_table(
            get_aou_relatedness_path(extension="1st_degrees.ht")
        )
        related_ht = hl.read_table(
            get_aou_util_path(name="relatedness_flagged_samples")
        )
        # TODO: ADD REAL Duplicated samples + 1st degree relatives
        meta_ht = meta_ht.annotate(
            related_0th_degree=hl.is_defined(related_0th_ht[meta_ht.key]),
            related=hl.is_defined(related_ht[meta_ht.key]),
            sex=hl.if_else(
                meta_ht.sex_at_birth == "Male",
                1,
                hl.if_else(meta_ht.sex_at_birth == "Female", 0, hl.missing(hl.tint32)),
            ),  # https://atgu.slack.com/archives/C056MJF70J1/p1695909765150529?thread_ts=1695326201.653099&cid=C056MJF70J1
            age=meta_ht.age_at_cdr,
            age2=meta_ht.age_at_cdr**2,
        )
        meta_ht = meta_ht.annotate(
            age_sex=meta_ht.age * meta_ht.sex, age2_sex=meta_ht.age2 * meta_ht.sex
        )
        meta_ht = meta_ht.annotate(**pop_ht[meta_ht.key])
        meta_ht = meta_ht.annotate(
            **{f"PC{x + 1}": meta_ht.pca_features[x] for x in range(16)},
            pop=meta_ht.ancestry_pred,
        )
        # TODO: add additional sample QC flags
        meta_ht = meta_ht.annotate(
            samples_to_keep=(meta_ht.age <= 100)
            & (hl.is_defined(meta_ht.sex))
            & (~meta_ht.related_0th_degree)
            & (~meta_ht.related)
            & (hl.is_defined(meta_ht.pop))
        )

        meta_ht.naive_coalesce(1000).write(
            get_sample_meta_path(annotation=True, extension="ht"), overwrite
        )
    meta_ht = hl.read_table(get_sample_meta_path(annotation=True, extension="ht"))
    # meta_ht.describe()
    meta_ht.summarize()
    return meta_ht


def annotate_phenotype_ht(overwrite: bool):
    if not hl.hadoop_exists(get_full_phenotype_path(annotation=True)) or overwrite:
        global_info = {}
        for trait_type in TRAIT_TYPEs:
            tag = "continuous" if trait_type == "physical_measurement" else trait_type
            tmp_ht = hl.read_table(get_raw_phenotype_path(name=f"{trait_type}_table"))
            global_info[tag] = list(tmp_ht.row_value)
        print(global_info)
        ht = hl.read_table(get_full_phenotype_path(annotation=False))
        meta_ht = hl.read_table(get_sample_meta_path(annotation=True, extension="ht"))
        ht = ht.annotate(**meta_ht[ht.key])
        ht = ht.annotate_globals(trait_type=global_info)
        ht.write(get_full_phenotype_path(annotation=True), overwrite)
    ht = hl.read_table(get_full_phenotype_path(annotation=True))
    # ht.export(get_full_phenotype_path(annotation=True, extension='txt.bgz'))
    return ht


def update_random_phenotype_ht(pops: str, test: bool, overwrite: bool):
    if (
        not hl.hadoop_exists(
            get_full_phenotype_path(annotation=True, random_pheno=True)
        )
        or overwrite
    ):
        meta_ht = hl.read_table(get_sample_meta_path(annotation=True, extension="ht"))
        pops = pops.split(",")
        full_ht = hl.import_table(
            get_random_phenotype_path(pop=pops[0], test=test),
            key="userId",
            impute=True,
            types={"userId": hl.tstr},
        )
        for pop in pops[1:]:
            print(pop)
            ht = hl.import_table(
                get_random_phenotype_path(pop=pop, test=test),
                key="userId",
                impute=True,
                types={"userId": hl.tstr},
            )
            full_ht = full_ht.union(ht)
        meta_ht = meta_ht.join(full_ht, how="left")
        meta_ht.summarize()
        meta_ht.describe()
        meta_ht.write(
            get_full_phenotype_path(annotation=True, random_pheno=True),
            overwrite=overwrite,
        )
    ht = hl.read_table(get_full_phenotype_path(annotation=True, random_pheno=True))
    print(ht.count())
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

        if args.update_random_phenotype_ht:
            print("-------------Merging random phenotype TSVs to a Main HT-----------")
            update_random_phenotype_ht(
                pops=args.pops, test=True, overwrite=args.overwrite
            )

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
        "--pops", help="comma-separated list", default="afr,amr,eas,eur,mid,sas"
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
        "--update-random-phenotype-ht",
        help="Update the full random phenotype Hail table",
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
