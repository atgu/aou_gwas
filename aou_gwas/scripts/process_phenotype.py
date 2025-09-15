#!/usr/bin/env python3

__author__ = "wlu"

import sys
print(sys.path)
from aou_gwas.utils.utils import *  # for QoB
from aou_gwas.utils.resources import *  # for QoB
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
        print(path)
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
            if name in ["r_drug_table", "pfhh_survey_table", "mcc2_phecode_table", "new_phecode_table_Dec23"]:
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
                ht = ht.annotate(
                    BMI_custom=ht.weight / ((ht.height / 100) ** 2),
                    WHR_custom=ht["waist-circumference-mean"] / ht["hip-circumference-mean"],
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

        if not hl.hadoop_exists(get_raw_phenotype_path(name=name, parsed=True, extension="txt")):
            phenos = list(ht.row_value)
            with hl.hadoop_open(
                get_raw_phenotype_path(name=name, parsed=True, extension="txt"), "w"
            ) as f:
                for pheno in phenos:
                    f.write(pheno + "\n")

def process_quantitative_phenotypes(quantitative_path, overwrite):
    path = get_raw_phenotype_path(name='processed_physical_measurement_table', parsed=True)
    if not hl.hadoop_exists(path) or overwrite:
        print(f'Processing {quantitative_path}')
        ht = hl.read_table(quantitative_path).key_by('person_id')
        ## Annotate BMI and WHR
        ht = ht.annotate(BMI=ht.weight / ((ht.height / 100) ** 2),
                         WHR=ht['waist-circumference-mean'] / ht['hip-circumference-mean'])

        ## Annotate WHRadjBMI
        beta = ht.aggregate(hl.agg.linreg(ht.WHR, ht.BMI))['beta'][0]
        ht = ht.annotate(WHRadjBMI=ht.WHR - beta * ht.BMI)

        ht.write(path, overwrite=overwrite)

    ht = hl.read_table(path)
    ht.describe()
    print(ht.count())
    if not hl.hadoop_exists(get_raw_phenotype_path(name='processed_physical_measurement_table', parsed=True, extension="txt")):
        phenos = list(ht.row_value)
        with hl.hadoop_open(
                get_raw_phenotype_path(name='processed_physical_measurement_table', parsed=True, extension="txt"), "w"
        ) as f:
            for pheno in phenos:
                f.write(pheno + "\n")
    return ht


def write_age_phenotype_ht(original_path=NEW_AGE_PHENO_PATHs, overwrite: bool = False):
    # paths = list(path['path'] for path in hl.hadoop_ls(original_path))
    # NOTE: Cannot run hl.hadoop_ls on QoB
    for path in original_path:
        name = path.split("/")[-1].split(".")[0]
        print(name)
        if not hl.hadoop_exists(get_raw_phenotype_path(name=name, parsed=False)) and name.startswith('new_mcc2_phecode'):
            ht = hl.import_table(path, delimiter=",", missing=".")
            print(f"--------{path}--------------")
            ht = ht.key_by("person_id")
            ht.describe()
            ht.checkpoint(
                get_raw_phenotype_path(name=name, parsed=False),
                _read_if_exists=not overwrite,
                overwrite=overwrite,
            )
        if not hl.hadoop_exists(get_raw_phenotype_path(name=name, parsed=True)) or overwrite:
            if name.startswith('new_mcc2_phecode'):
                ht = hl.read_table(get_raw_phenotype_path(name=name, parsed=False))
                ht = ht.annotate(
                    **{
                        field: parse_empty_missing(ht[field], hl.bool)
                        for field in list(ht.row_value) if not field.endswith('_age')
                    },
                )
                ht = ht.annotate(
                    **{
                        field: parse_empty_missing(ht[field], hl.float64)
                        for field in list(ht.row_value) if field.endswith('_age')
                    },
                )
            else:
                print(
                    f"-------- Parsing {name} to {get_raw_phenotype_path(name=name, parsed=True)}--------------"
                )
                ht = hl.import_table(path, delimiter=",", missing=".", impute=True, types={"person_id": hl.tstr})
                ht = ht.key_by("person_id")
                if name == 'all_first_last_billing_ages':
                    ht = ht.drop(ht[''])
            ht.describe()
            ht = ht.repartition(1000)
            ht.checkpoint(
                get_raw_phenotype_path(name=name, parsed=True), overwrite=overwrite
            )
        print(name)
        ht = hl.read_table(get_raw_phenotype_path(name=name, parsed=True))
        ht.describe()

        if not hl.hadoop_exists(get_raw_phenotype_path(name=name, parsed=True, extension="txt")):
            phenos = list(ht.row_value)

            with hl.hadoop_open(
                get_raw_phenotype_path(name=name, parsed=True, extension="txt"), "w"
            ) as f:
                for pheno in phenos:
                    if pheno.split('_')[-1] == 'age':
                        continue
                    f.write(pheno + "\n")


def merge_phenotype_hts(overwrite: bool):
    if not hl.hadoop_exists(get_full_phenotype_path(annotation=False)) or overwrite:
        ht = hl.read_table(NEW_QUANTITATIVE_HT_PATH)
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
    if not hl.hadoop_exists(get_pca_ht_path(pop='all', name=f'pruned_full_scores')):
        ht = hl.read_table(get_pca_ht_path(pop='afr', name=f'pruned_full_scores'))
        ht = ht.annotate(pop = 'afr')
        for pop in POPS[1:6]:
            print(f'Loading {pop} PCs (R2)...')
            pca_r2_ht = hl.read_table(get_pca_ht_path(pop=pop, name=f'pruned_full_scores'))
            pca_r2_ht = pca_r2_ht.annotate(pop = pop)
            ht = ht.union(pca_r2_ht)
            ht.drop('scores')
        ht.checkpoint(get_pca_ht_path(pop='all', name=f'pruned_full_scores'), overwrite=True)
    pca_ht = hl.read_table(get_pca_ht_path(pop='all', name=f'pruned_full_scores'))
    pca_ht.describe()
    print(f'Number of pruned samples: {pca_ht.count()}')
    print(f'Number of pruned samples per group: {pca_ht.aggregate(hl.agg.counter(pca_ht.pop))}')
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
        print(f'Loading 0th related sample table: {ORIGINAL_0th_RELATED_PATH}')
        related_0th_ht = hl.import_table(
            ORIGINAL_0th_RELATED_PATH,
            # Produced by ~/Dropbox (Partners HealthCare)/github_repo/aou_gwas/R/relatedness_cluster_check.R
            no_header=True)  # n=570
        print(f'Number of 0th related samples: {related_0th_ht.count()}')
        related_0th_ht = related_0th_ht.key_by('f0')
        related_ht = hl.read_table(
            get_aou_util_path(name="relatedness_flagged_samples")
        )
        pop_pruned_ht = hl.import_table(CENTROID_PRUNED_SAMPLES, types={'s': hl.tstr}).key_by('s')
        print(f'Number of PC pruned samples: {pop_pruned_ht.count()}')
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
            pruned_samples = hl.is_defined(pop_pruned_ht[meta_ht.key]),
        )
        meta_ht = meta_ht.annotate(
            age_sex=meta_ht.age * meta_ht.sex, age2_sex=meta_ht.age2 * meta_ht.sex
        )
        meta_ht = meta_ht.annotate(pop=pop_ht[meta_ht.key].ancestry_pred)
        # meta_ht = meta_ht.annotate(
        #     **{f"pop_PC{i + 1}": pca_ht[meta_ht.key][f'PC{i+1}'] for i in range(50)}
        # )

        meta_ht = meta_ht.annotate(
            samples_to_keep=(meta_ht.age <= 100)
            & (hl.is_defined(meta_ht.sex))
            & (~meta_ht.related_0th_degree)
            & (hl.is_defined(meta_ht.pop))
        )
        meta_ht.describe()

        meta_ht.naive_coalesce(1000).write(
            get_sample_meta_path(annotation=True, extension="ht"), overwrite
        )
    meta_ht = hl.read_table(get_sample_meta_path(annotation=True, extension="ht"))
    meta_ht.summarize()
    return meta_ht


def annotate_phenotype_ht(overwrite: bool):
    if not hl.hadoop_exists(get_full_phenotype_path(annotation=True)) or overwrite:
        global_info = {}
        for trait_type in TRAIT_TYPEs:
            tag = "continuous" if trait_type == "processed_physical_measurement" else trait_type
            tmp_ht = hl.read_table(get_raw_phenotype_path(name=f"{trait_type}_table"))
            global_info[tag] = list(tmp_ht.row_value)
        print(global_info)
        ht = hl.read_table(get_full_phenotype_path(annotation=False))
        meta_ht = hl.read_table(get_sample_meta_path(annotation=True, extension="ht"))
        ht = ht.annotate(**meta_ht[ht.key])
        ht = ht.annotate_globals(trait_type=global_info)
        ht.write(get_full_phenotype_path(annotation=True), overwrite)
    ht = hl.read_table(get_full_phenotype_path(annotation=True))
    ht.describe()
    ht.export(get_full_phenotype_path(annotation=True, extension='txt.bgz'))
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
    print(ht.summarize())
    ht.describe()
    return ht

def process_lab_measurements(paths, overwrite):
    codes = [path.replace('.csv', '').replace('gs://allxall-phenotypes/data/measurement_', '') for path in paths]
    if not hl.hadoop_exists(get_raw_phenotype_path(name='labs_metadata_04242024', parsed=True)):
        # path = f"{ORIGINAL_PHENO_ROOT}/top100_labs_concept_ids.csv",
        # print(path)
        # # Detail check: https://batch.hail.is/batches/8107870/jobs/1
        # ht = hl.import_table(path, delimiter='(?<!\s),(?!\s)', no_header=True, impute=True)
        # ht = ht.filter(~hl.literal(['Number of Participants', "Total"]).contains(ht.f2))
        # ht = ht.rename({
        #     'f0': 'concept_name',
        #     'f1': 'measurement_concept_id',
        #     'f2': 'n_participants_total', # Number of Participants Total
        #     'f3': 'n_participants_categorical', # Number of Participants with Categorical Lab Values
        #     'f4': 'n_participants_numeric', # Number of Participants with Numeric Lab Values
        #     'f5': 'n_participants_both', # Number of Participants with Both Types of Lab Values
        #     'f6': 'p_participants_total', # Percentage of Participants Total
        #     'f7': 'p_participants_categorical', # Percentage of Participants with Categorical Lab Values
        #     'f8': 'p_participants_numeric', # Percentage of Participants with Numeric Lab Values
        #     'f9': 'p_participants_both', # Percentage of Participants with Both Types of Lab Values
        #      })
        # ht = ht.key_by('measurement_concept_id')
        # ht = ht.annotate(
        #     **{
        #         field: parse_empty_missing(ht[field], hl.int64)
        #         for field in list(ht.row_value) if (field not in ['concept_name'])
        #     },
        # )  # These works with dataproc
        path = f"{ORIGINAL_PHENO_ROOT}/labs_metadata_04242024.csv"
        ht = hl.import_table(path, delimiter='(?<!\s),(?!\s)', impute=True, types = {'measurement_concept_id':hl.tstr})
        combined_codes = ['3013466_3018677', '3043359_3027450', '3045440_3023939', '3045792_3016921', '3021614_3024763_3015688']
        for i in range(len(combined_codes)):
            tmp_codes = combined_codes[i].split('_')
            ht = ht.annotate(measurement_concept_id = hl.if_else(hl.literal(tmp_codes).contains(ht.measurement_concept_id), combined_codes[i], ht.measurement_concept_id))
        ht = ht.key_by('measurement_concept_id').distinct()
        ht.describe()
        ht = ht.checkpoint(get_raw_phenotype_path(name='labs_metadata_04242024', parsed=True), _read_if_exists= not overwrite, overwrite=overwrite)
        print(ht.count())
        ht.show(1000)
    ht = hl.read_table(get_raw_phenotype_path(name='labs_metadata_04242024', parsed=True))
    ht = ht.filter(hl.literal(codes).contains(ht.measurement_concept_id))
    ht.show(len(codes))
    meta_ht = hl.read_table(get_sample_meta_path(annotation=True, extension="ht"))
    pca_ht = hl.read_table(get_pca_ht_path(pop='all', name=f'pruned_full_scores'))
    from tqdm import tqdm
    for i in tqdm(range(len(codes))):
        code = codes[i]
        if not hl.hadoop_exists(get_lab_measurements_path(name=code, parsed=True)) or overwrite:
            path = f'{ORIGINAL_PHENO_ROOT}/measurement_{code}.csv'
            print(path)
            ht = hl.import_table(path, delimiter='(?<!\s),(?!\s)',
                                    missing=".",
                                    impute=True,
                                    types={"person_id": hl.tstr})
            ht = ht.key_by('person_id')
            ht = ht.annotate(**meta_ht[ht.key])
            ht = ht.checkpoint(get_lab_measurements_path(name=code, parsed=False),  _read_if_exists= True)
            print(ht['count'].summarize())
            ht = ht.filter(ht.samples_to_keep & ht.pruned_samples)
            ht = ht.annotate(**pca_ht[ht.key])
            ht = ht.annotate(value = ht.median)
            ht.checkpoint(get_lab_measurements_path(name=code, parsed=True), overwrite=overwrite,
                          _read_if_exists=not overwrite)
        ht = hl.read_table(get_lab_measurements_path(name=code, parsed=True))
        ht.describe()
        print(ht.count())


def generate_lab_n_cases(overwrite: bool):
    lab_case_cnt = {}
    if not hl.hadoop_exists(get_raw_phenotype_path(name='labs_metadata_n_cases', parsed=True)) or overwrite:
        for i in tqdm(range(len(LAB_CODES))):
            code = LAB_CODES[i]
            path = get_lab_measurements_path(name=code, parsed=True)
            ht = hl.read_table(path)
            print(path)
            sub_ht = ht.filter(ht.samples_to_keep & ht.pruned_samples)
            tmp_cnt = sub_ht.aggregate(hl.agg.group_by(sub_ht.pop, hl.agg.count_where(hl.is_defined(sub_ht.value))))
            lab_case_cnt[code] = tmp_cnt
            print(tmp_cnt)
        df = pd.DataFrame(lab_case_cnt)
        print(df)
        df = df.melt(ignore_index=False).reset_index().rename(
            columns={'index': 'pop', 'variable': 'phenoname', 'value': 'n_cases'})
        print(df)
        ht = hl.Table.from_pandas(df, key='phenoname')
        info_ht = hl.read_table(get_raw_phenotype_path(name='labs_metadata_04242024', parsed=True))
        ht = ht.join(info_ht)
        ht.checkpoint(get_raw_phenotype_path(name='labs_metadata_n_cases', parsed=True), _read_if_exists=not overwrite,
                      overwrite=overwrite)
    ht = hl.read_table(get_raw_phenotype_path(name='labs_metadata_n_cases', parsed=True))
    ht = ht.select('pop', 'n_cases')
    print(ht.count())

    df = ht.to_pandas()
    df = df[df.n_cases >= 200]
    dict = df.filter(['pop', 'phenoname']).drop_duplicates().groupby('pop')['phenoname'].apply(list).to_dict()
    write_pickle_dict(output=get_phenotype_info_path(version="lab_by_pop_dict", extension="dict"), dict=dict)
    dict = read_pickle_dict(get_phenotype_info_path(version="lab_by_pop_dict", extension="dict"))
    print([f"{pop.upper()}: {len(dict[pop])}" for pop in pops])



def write_phenotype_mt(path, raw, overwrite: bool = False):
    meta_ht = hl.read_table(get_sample_meta_path(annotation=True, extension="ht"))
    pca_ht = hl.read_table(get_pca_ht_path(pop='all', name=f'pruned_full_scores'))
    # info_ht = hl.read_table(get_phenotype_info_path(version='R', extension='ht'))
    # info_ht = info_ht.select(
    #     **{field: info_ht[field] for field in list(info_ht.row_value) if not field.startswith('n_')})
    name = path.split("/")[-1].split(".")[0].replace('_250k', '')
    print(name)
    if not hl.hadoop_exists(get_raw_phenotype_path(name=name, parsed=False, extension='mt')):
        print(f"--------Exporting {path} to MT--------------")
        if not path.endswith('.ht'):
            mt = hl.import_matrix_table(path,
                                        row_fields={'person_id': hl.tstr},
                                        row_key='person_id', delimiter=",", entry_type=hl.tstr)
            mt = mt.rename({'col_id': 'phenocode', 'x': 'value'})
        else:
            ht = hl.read_table(path)
            ht.describe()
            mt = ht.to_matrix_table_row_major(
                columns=list(ht.row_value), entry_field_name="value", col_field_name="phenocode"
            ).key_cols_by()
        mt.describe()
        mt.checkpoint(
            get_raw_phenotype_path(name=name, parsed=False, extension='mt'), _read_if_exists=not overwrite, overwrite=overwrite
        )
    tag = "hard_filtered." if not raw else ""
    out_path = get_raw_phenotype_path(name=name, parsed=True, extension=f'{tag}mt')
    if not hl.hadoop_exists(out_path) or overwrite:
        print(
            f"-------- Parsing {name} to {out_path}--------------"
        )
        mt = hl.read_matrix_table(get_raw_phenotype_path(name=name, parsed=False, extension='mt'))
        mt = mt.key_cols_by()
        TRAIT_TYPEs = {'new_phecode_table_Dec23' : 'phecodeX',
                       'mcc2_phecode_table': 'phecode',
                       'processed_physical_measurement_table': 'physical_measurement',
                       'r_drug_table': 'r_drug',
                       'pfhh_survey_table': 'pfhh_survey'}
        mt = mt.annotate_cols(trait_type = TRAIT_TYPEs[name])
        # mt = mt.annotate_cols(phenocode=hl.if_else(mt.phenocode.endswith('_irnt'), mt.phenocode.split('_')[0],
        #                                            mt.phenocode))
        mt = mt.key_cols_by('trait_type', 'phenocode')
        mt = mt.annotate_rows(**meta_ht[mt.row_key])
        if not raw:
            print('Filter to pruned samples....')
            mt = mt.filter_rows(mt.samples_to_keep & mt.pruned_samples)
        if name != "processed_physical_measurement_table":
            mt = mt.annotate_entries(value = hl.if_else(mt.value == "NA", "", mt.value))
            mt = mt.annotate_entries(value=parse_empty_missing(mt.value, hl.bool))
            mt = mt.annotate_cols(
               #  **info_ht[mt.col_key],
                n_cases=hl.agg.filter(hl.is_defined(mt.pop), hl.agg.count_where(mt.value)),
                n_controls=hl.agg.filter(hl.is_defined(mt.pop), hl.agg.count_where(~mt.value)),
                n_samples_defined=hl.agg.filter(hl.is_defined(mt.pop), hl.agg.count_where(hl.is_defined(mt.value))),
                n_cases_by_sex=hl.agg.group_by(mt.sex, hl.agg.filter(hl.is_defined(mt.pop), hl.agg.count_where(mt.value))),
                n_controls_by_sex=hl.agg.group_by(mt.sex,
                                                  hl.agg.filter(hl.is_defined(mt.pop), hl.agg.count_where(~mt.value))),
                n_samples_defined_by_sex=hl.agg.group_by(mt.sex, hl.agg.filter(hl.is_defined(mt.pop),
                                                                               hl.agg.count_where(hl.is_defined(mt.value)))),
                n_cases_by_pop=hl.agg.group_by(mt.pop, hl.agg.count_where(mt.value)),
                n_controls_by_pop=hl.agg.group_by(mt.pop, hl.agg.count_where(~mt.value)),
                n_samples_defined_by_pop=hl.agg.group_by(mt.pop, hl.agg.count_where(hl.is_defined(mt.value))),
                n_cases_raw=hl.agg.count_where(mt.value),
                n_controls_raw=hl.agg.count_where(~mt.value),
                n_samples_defined_raw=hl.agg.count_where(hl.is_defined(mt.value)),
                n_cases_by_sex_raw=hl.agg.group_by(mt.sex, hl.agg.count_where(mt.value)),
                n_controls_by_sex_raw=hl.agg.group_by(mt.sex, hl.agg.count_where(~mt.value)),
                n_samples_defined_by_sex_raw=hl.agg.group_by(mt.sex, hl.agg.count_where(hl.is_defined(mt.value))),
            )
        else:
            mt = mt.annotate_cols(
               #  **info_ht[mt.col_key],
                n_cases=hl.agg.filter(hl.is_defined(mt.pop), hl.agg.count_where(hl.is_defined(mt.value))),
                n_controls=hl.missing(hl.tint64),
                n_samples_defined=hl.agg.filter(hl.is_defined(mt.pop), hl.agg.count_where(hl.is_defined(mt.value))),
                n_cases_by_sex=hl.agg.group_by(mt.sex, hl.agg.filter(hl.is_defined(mt.pop), hl.agg.count_where(hl.is_defined(mt.value)))),
                n_controls_by_sex=hl.missing(hl.tint64),
                n_samples_defined_by_sex=hl.agg.group_by(mt.sex, hl.agg.filter(hl.is_defined(mt.pop),
                                                                               hl.agg.count_where(hl.is_defined(mt.value)))),
                n_cases_by_pop=hl.agg.group_by(mt.pop, hl.agg.count_where(hl.is_defined(mt.value))),
                n_controls_by_pop=hl.missing(hl.tint64),
                n_samples_defined_by_pop=hl.agg.group_by(mt.pop, hl.agg.count_where(hl.is_defined(mt.value))),
                n_cases_raw=hl.agg.count_where(hl.is_defined(mt.value)),
                n_controls_raw=hl.missing(hl.tint64),
                n_samples_defined_raw=hl.agg.count_where(hl.is_defined(mt.value)),
                n_cases_by_sex_raw=hl.agg.group_by(mt.sex, hl.agg.count_where(hl.is_defined(mt.value))),
                n_controls_by_sex_raw=hl.missing(hl.tint64),
                n_samples_defined_by_sex_raw=hl.agg.group_by(mt.sex, hl.agg.count_where(hl.is_defined(mt.value))),
            )
        mt = mt.annotate_rows(**pca_ht[mt.row_key])
        mt = mt.annotate_cols(**{f'n_cases_{pop}': mt.n_cases_by_pop.get(pop) for pop in POPS})
        mt = mt.repartition(1000)
        mt = mt.checkpoint(
            out_path, overwrite=overwrite, _read_if_exists=not overwrite
        )
        mt.describe()
        mt.cols().write(get_phenotype_info_path(version=name, extension=f'{tag}ht'), overwrite)
        ht = hl.read_table(get_phenotype_info_path(version=name, extension=f'{tag}ht'))
        ht.describe()
        ht.export(get_phenotype_info_path(version=name, extension=f'{tag}tsv'))
        if name == 'processed_physical_measurement_table':
            mt = mt.annotate_rows(n_cases_per_sample=hl.agg.count_where(hl.is_defined(mt.value)),
                                  n_phenotypes_defined_per_sample=hl.agg.count_where(hl.is_defined(mt.value)))
        else:
            mt = mt.annotate_rows(n_cases_per_sample=hl.agg.count_where(mt.value),
                                  n_phenotypes_defined_per_sample=hl.agg.count_where(hl.is_defined(mt.value)))
        if raw:
            mt.rows().export(get_phenotype_info_path(version=name, extension='persample.tsv'))
    mt = hl.read_matrix_table(get_raw_phenotype_path(name=name, parsed=True, extension=f'{tag}mt'))
    print(mt.count())
    mt.cols().show()
    if not hl.hadoop_exists(get_raw_phenotype_path(name=name, parsed=True, extension="txt")):
        phenos = mt.phenocode.collect()
        print(phenos)

        with hl.hadoop_open(
                get_raw_phenotype_path(name=name, parsed=True, extension="txt"), "w"
        ) as f:
            for pheno in phenos:
                if pheno.split('_')[-1] == 'age':
                    continue
                f.write(pheno + "\n")

def write_age_mt(raw, overwrite):
    if not hl.hadoop_exists(get_raw_phenotype_path(name='age', parsed=False, extension='mt')):
        print('-----------Exporting raw Age MT------------')
        AGE_NAMES = ["new_mcc2_phecode_table1", "new_mcc2_phecode_table2", "new_mcc2_phecode_table3", "new_mcc2_phecode_table4", "new_mcc2_phecode_table5", "new_mcc2_phecode_table6"]
        ht = hl.read_table(get_raw_phenotype_path(name="new_drug_table", parsed=True))
        ht = ht.select(
            **{f'{field.replace("_age", "")}': hl.float64(ht[field]) for field in list(ht.row_value) if field.endswith('_age')}
        )
        r_drug_name = list(ht.row_value)
        mt = ht.to_matrix_table_row_major(
            columns=list(ht.row_value), entry_field_name="age", col_field_name="phenocode"
        )
        for name in AGE_NAMES:
            print(name)
            tmp_ht = hl.read_table(get_raw_phenotype_path(name=name, parsed=True))
            tmp_ht = tmp_ht.select(
            **{f'{field}': hl.float64(tmp_ht[field]) for field in list(tmp_ht.row_value) if field.endswith('_age')})
            tmp_ht = tmp_ht.checkpoint(f'gs://aou_tmp/phenotype_{name}.ht', _read_if_exists=True)
            tmp_mt = tmp_ht.to_matrix_table_row_major(
                columns=list(tmp_ht.row_value), entry_field_name="age", col_field_name="phenocode"
            )

            mt = mt.union_cols(tmp_mt, row_join_type='outer')
            mt = mt.checkpoint(f'gs://aou_tmp/phenotype_age_raw_{name}.mt', _read_if_exists=True)

        mt = mt.annotate_cols(
            trait_type = hl.if_else(hl.literal(r_drug_name).contains(mt.phenocode), 'r_drug', 'phecode'),
            modifier=hl.if_else(mt.phenocode.endswith('_irnt'), 'irnt', ''))
        mt = mt.key_cols_by('trait_type', 'phenocode', 'modifier')

        mt.describe()
        mt.checkpoint(
            get_raw_phenotype_path(name='age', parsed=False, extension='mt'), _read_if_exists=True
        )

    if not hl.hadoop_exists(get_raw_phenotype_path(name='age', parsed=True, extension=f'{"hard_filtered." if not raw else ""}mt')) or overwrite:
        print('-----------Annotating Age MT------------')
        mt = hl.read_matrix_table(get_raw_phenotype_path(name='age', parsed=False, extension='mt'))
        mt.show()
        print(mt.count())
        meta_ht = hl.read_table(get_sample_meta_path(annotation=True, extension="ht"))
        meta_ht = meta_ht.rename({'age': 'sample_age'})
        info_ht = hl.read_table(get_phenotype_info_path(version='R', extension='ht'))
        age_ht = hl.read_table(get_raw_phenotype_path(name='all_first_last_billing_ages', parsed = True))
        info_ht = info_ht.select(
            **{field: info_ht[field] for field in list(info_ht.row_value) if not field.startswith('n_')})
        mt = mt.annotate_rows(**meta_ht[mt.row_key], **age_ht[mt.row_key])
        if not raw:
            mt = mt.filter_rows(mt.samples_to_keep)
        mt = mt.annotate_cols(
            **info_ht[mt.col_key],
            n_samples_age_defined=hl.agg.filter(hl.is_defined(mt.pop), hl.agg.count_where(hl.is_defined(mt.age))),
            n_samples_age_defined_by_sex=hl.agg.group_by(mt.sex, hl.agg.filter(hl.is_defined(mt.pop),
                                                                           hl.agg.count_where(
                                                                               hl.is_defined(mt.age)))),
            n_samples_age_defined_by_pop=hl.agg.group_by(mt.pop, hl.agg.count_where(hl.is_defined(mt.age))),
            n_samples_age_defined_raw=hl.agg.count_where(hl.is_defined(mt.age)),
            n_samples_age_defined_by_sex_raw=hl.agg.group_by(mt.sex, hl.agg.count_where(hl.is_defined(mt.age))),
        )

        mt = mt.repartition(1000)
        mt.checkpoint(
            get_raw_phenotype_path(name='age', parsed=True, extension=f'{"hard_filtered." if not raw else ""}mt'), overwrite=overwrite,
            _read_if_exists=not overwrite
        )
    mt = hl.read_matrix_table(get_raw_phenotype_path(name='age', parsed=True, extension=f'{"hard_filtered." if not raw else ""}mt'))
    mt.describe()
    print(mt.count())


def main(args):
    try:
        hl.init(
            app_name=f'Process_phenotypes',
            tmp_dir=TMP_BUCKET,
            driver_memory="highmem",
            driver_cores=8,
            worker_memory="highmem",
            worker_cores=1,
            default_reference="GRCh38",
            log="/process_phenotypes.log"
        )

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
                # original_path=ORIGINAL_PHENO_PATHs,
                original_path= [f"{ORIGINAL_PHENO_ROOT}/new_phecode_table_Dec23.csv"],
                overwrite=args.overwrite
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

        if args.write_phenotype_mt:
            print("-------------Write phenotype MT-----------")
            PATHs = [
                    # NEW_QUANTITATIVE_HT_PATH,
                    f"{ORIGINAL_PHENO_ROOT}/new_phecode_table_Dec23.csv",
                    f"{ORIGINAL_PHENO_ROOT}/pfhh_survey_table.csv",
                    f"{ORIGINAL_PHENO_ROOT}/mcc2_phecode_table.csv",
                    f"{ORIGINAL_PHENO_ROOT}/r_drug_table.csv"
            ]
            for path in PATHs:
                if args.batch:
                    j = b.new_python_job(name=f"Create_phenotype_MT")
                    j.call(update_phenotype_mt, args.overwrite)
                else:
                    # write_phenotype_mt(path=path, raw=True, overwrite=args.overwrite)
                    write_phenotype_mt(path=path, raw=False, overwrite=args.overwrite)

        if args.write_age_mt:
            print("-------------Write age MT-----------")
            write_age_mt(overwrite=args.overwrite)

        if args.update_raw_age_phenotypes:
            print("-------------Writing Raw AGE Phenotype HTs-----------")
            if args.batch:
                j = b.new_python_job(name=f"Write_raw_age_phenotype_HTs")
                j.call(write_age_phenotype_ht, original_path=NEW_AGE_PHENO_PATHs, overwrite=args.overwrite)
            else:
                write_age_phenotype_ht(
                    # original_path=NEW_AGE_PHENO_PATHs,
                    original_path= [f"{ORIGINAL_PHENO_ROOT}/new_mcc2_phecode_table4.csv",
                                    f"{ORIGINAL_PHENO_ROOT}/new_mcc2_phecode_table5.csv",
                                    f"{ORIGINAL_PHENO_ROOT}/new_mcc2_phecode_table6.csv",
                                    f"{ORIGINAL_PHENO_ROOT}/new_mcc2_phecode_table1.csv",
                                    f"{ORIGINAL_PHENO_ROOT}/new_mcc2_phecode_table2.csv",
                                    f"{ORIGINAL_PHENO_ROOT}/new_mcc2_phecode_table3.csv"],
                    overwrite=args.overwrite
                )

        if args.update_quantitative_phenotype_ht:
            print("-------------Writing Processed Quantitative HT-----------")
            process_quantitative_phenotypes(
                quantitative_path=get_raw_phenotype_path(name='physical_measurement_table', parsed=True), overwrite=args.overwrite
            )

        if args.update_lab_measurements_ht:
            print("-------------Writing Processed lab measurement HTs-----------")
            process_lab_measurements(paths=LAB_MEASUREMENT_PATHs, overwrite=args.overwrite)
            generate_lab_n_cases(overwrite=args.overwrite)

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
        "--update-quantitative-phenotype-ht",
        help="Update the quantitative Hail table",
        action="store_true",
    )
    parser.add_argument(
        "--update-lab-measurements-ht",
        help="Update lab measurement Hail tables",
        action="store_true",
    )
    parser.add_argument(
        "--update-raw-age-phenotypes",
        help="Update the raw age phenotype Hail tables",
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
        "--write-phenotype-mt",
        help="Update the processed phenotype Hail Matrixtable",
        action="store_true",
    )

    parser.add_argument(
        "--write-age-mt",
        help="Write age information to Hail Matrixtable",
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
