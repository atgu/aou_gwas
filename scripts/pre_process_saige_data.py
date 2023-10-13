#!/usr/bin/env python3

__author__ = "wlu"

import argparse

# from aou_gwas import *  # for dataproc
from utils.utils import *  # for QoB
from utils.annotations import * # for QoB
from utils.resources import *  # for QoB
from gnomad.utils.vep import *
from gnomad.utils.filtering import *

# from ukbb_pan_ancestry import *


# ALL PATHs: https://support.researchallofus.org/hc/en-us/articles/4616869437204-Controlled-CDR-Directory
# Descriptions: https://support.researchallofus.org/hc/en-us/articles/4614687617556-How-the-All-of-Us-Genomic-data-are-organized
GENE_INTERVAL_PATH = (
    f"gs://aou_wlu/data/group_positions_{N_GENE_PER_GROUP}_protein_coding.ht"
)


def main(args):
    app_name = None
    if args.run_vep:
        app_name = 'run_vep'
    if args.run_relatedness:
        app_name = f'run_relatedness_{args.relatedness_type}'
    hl.init(
        tmp_dir='gs://aou_tmp/',
        driver_memory="highmem",
        driver_cores=8,
        worker_memory="highmem",
        worker_cores=1,
        default_reference="GRCh38",
        log="/pre_process_saige_data.log",
        app_name='gene_map'
    )

    # Load info tables
    vat_ht = hl.read_table(get_aou_util_path(name="vat"))
    vat_ht = vat_ht.collect_by_key()
    # vat_ht = vat_ht.filter((vat_ht.locus.contig == 'chr1') & (vat_ht.locus.position < 10330))
    vat_ht = vat_ht.filter((vat_ht.locus.contig == 'chr21'))

    if args.run_vep:
        # Vep resource from gnomAD:
        # https://github.com/broadinstitute/gnomad_methods/blob/e2f2b055a2d40991851e7852dea435b4eaef45ea/gnomad/resources/grch38/reference_data.py#L100
        if (not hl.hadoop_exists(get_aou_util_path(name="vep_corrected"))) or args.overwrite:
            # vep_ht = hl.vep(vat_ht, csq=False)
            vep_ht = vep_or_lookup_vep(vat_ht, vep_version="vep105")
            vep_ht.write(
                get_aou_util_path(name="vep_corrected"),
                overwrite=args.overwrite,
            )
        vep_ht = hl.read_table(get_aou_util_path(name="vep_corrected"))
        vep_ht.describe()

    if args.create_gene_mapping_file:  # TODO: run this once VEP table is ready + add a freq filter
        # ht = hl.read_table(get_aou_util_path(name="vep_corrected"))  # TODO: Write VEP table
        ht = hl.read_table('gs://gcp-public-data--gnomad/resources/context/grch38_context_vep_annotated.v105.ht') # use the gnomAD context table for now
        ht.describe()
        # qual_ht = hl.read_table(get_aou_util_path(name='variant_qc'))  # TODO: Write table from VDS
        # qual_ht = qual_ht.filter(hl.len(qual_ht.filters) == 0)
        call_stats_ht = hl.read_table(get_aou_util_path(name="vat"))
        call_stats_ht = call_stats_ht.collect_by_key()
        pops = args.pops.split(",") if (args.pops is not None) else POPS
        for pop in pops:
            print(pop)
            call_stats_ht.values[f"gvs_{pop}_af"].show()
            max_an = call_stats_ht.aggregate(
                hl.struct(
                    autosomes=hl.agg.max(call_stats_ht.values[f"gvs_{pop}_an"][0]),
                    x=hl.agg.filter(
                        call_stats_ht.locus.in_x_nonpar(),
                        hl.agg.max(call_stats_ht.values[f"gvs_{pop}_an"][0])
                    ),
                    y=hl.agg.filter(
                        call_stats_ht.locus.in_y_nonpar(),
                        hl.agg.max(call_stats_ht.values[f"gvs_{pop}_an"][0]),
                    ),
                ),
            )

            an = call_stats_ht.values[f"gvs_{pop}_an"][0]
            call_stats_ht = call_stats_ht.filter(
                hl.case()
                .when(call_stats_ht.locus.in_x_nonpar(), an > 0.8 * max_an.x)
                .when(call_stats_ht.locus.in_y_nonpar(), an > 0.8 * max_an.y)
                .default(an > 0.8 * max_an.autosomes)
            )

            ht = ht.annotate(freq= call_stats_ht[ht.key]['values'][f"gvs_{pop}_af"][0])
            # ht = ht.filter(hl.is_defined(qual_ht[ht.key]) & hl.is_defined(ht.freq))
            print('Yes')
            gene_map_ht = create_gene_map_ht(ht, freq_field=ht.freq)
            gene_map_ht.describe()
            gene_map_ht.write(
                get_aou_util_path(name="gene_map", parsed=False), args.overwrite
            )

            gene_map_ht = hl.read_table(
                get_aou_util_path(name="gene_map", parsed=False)
            )
            gene_map_ht = post_process_gene_map_ht(gene_map_ht, freq_cutoff=0.001)
            gene_map_ht.write(
                get_aou_util_path(name="gene_map", parsed=True), args.overwrite
            )
            break

    if args.get_duplicated_samples:
        # TODO: follow up on zulip https://hail.zulipchat.com/#narrow/stream/123010-Hail-Query-0.2E2-support/topic/ClassCastException/near/394696118
        if (not hl.hadoop_exists(get_aou_relatedness_path(extension="duplicates.ht"))) or args.overwrite:
            print('Generating duplicated samples...')
            degree_1st_cutoff = 0.354
            ht = hl.read_table(get_aou_util_path(name="relatedness"))
            ht = ht.filter(ht.kin > degree_1st_cutoff)
            print('--------------------Running hl.maximal_independent_set()----------------')
            duplicated_samples_to_remove = hl.maximal_independent_set(ht['i.s'], ht['j.s'], keep=False)
            duplicated_samples_to_remove.describe()
            duplicated_samples_to_remove.checkpoint(get_aou_relatedness_path(extension="duplicates.ht"),
                                                    overwrite=args.overwrite,
                                                    _read_if_exists= not args.overwrite)
        duplicated_samples_to_remove = hl.read_table(get_aou_relatedness_path(extension="duplicates.ht"))
        print(f'N duplicated samples: {duplicated_samples_to_remove.count()}')

        if (not hl.hadoop_exists(get_aou_relatedness_path(extension="1st_degrees.ht"))) or args.overwrite:
            print('Generating 0th degree samples...')
            ht = hl.read_table(get_aou_util_path(name="relatedness", parsed=True))
            ht = ht.filter(ht.kin > 0.354)
            ht = ht.key_by()
            ht_i = ht.key_by(s=ht['i.s']).select()
            ht_j = ht.key_by(s=ht['j.s']).select()
            samples_1st_degree = ht_i.union(ht_j).distinct()
            samples_1st_degree.describe()
            samples_1st_degree.checkpoint(get_aou_relatedness_path(extension="1st_degrees.ht"),
                                                    overwrite=args.overwrite,
                                                    _read_if_exists= not args.overwrite)
        samples_1st_degree = hl.read_table(get_aou_relatedness_path(extension="1st_degrees.ht"))
        print(f'N duplicated samples: {samples_1st_degree.count()}')

    if args.run_relatedness:
        MAF_CUTOFF = 0.05
        mt = hl.read_matrix_table(ACAF_MT_PATH)  # (99250816, 245394)
        mt.describe()
        # print(mt.info.AF.summarize())
        # # - AF (array<float64>):
        # #   Non-missing: 99250796 (100.00%)
        # #       Missing: 20 (0.00%)
        # #      Min Size: 1
        # #      Max Size: 1
        # #     Mean Size: 1.00
        # #
        # #   - AF[<elements>] (float64):
        # #     Non-missing: 99250796 (100.00%)
        # #         Missing: 0
        # #         Minimum: 0.00
        # #         Maximum: 1.00
        # #            Mean: 0.04
        # #         Std Dev: 0.13
        # print(mt.aggregate_rows(hl.agg.count_where(mt.info.AF[0] > MAF_CUTOFF))) # 10289329
        ### Filter to 1M random common variants before running IBD

        if (
                not hl.hadoop_exists(get_aou_relatedness_path(extension="1Mvar.ht"))
        ) or args.overwrite:
            print(f"-------------Downsampling to 1M common variants-----------")
            ht = mt.rows()
            sampled_variants = ht.aggregate(
                hl.agg.filter(
                    ht.info.AF[0] >= MAF_CUTOFF,
                    hl.agg._reservoir_sample(ht.key, 1000000),
                )
            )
            variants = [variant for variant in sampled_variants]
            ht = hl.Table.parallelize(variants).key_by(*ht.key.keys())
            ht = ht.checkpoint(
                get_aou_relatedness_path(extension="1Mvar.ht"),
                _read_if_exists=(not args.overwrite),
                overwrite=args.overwrite,
            )
        ht = hl.read_table(get_aou_relatedness_path(extension="1Mvar.ht"))
        # print(f'N variants: {ht.count()}')
        mt = mt.filter_rows(hl.is_defined(ht[mt.row_key]))
        if args.relatedness_type == "ibd":
            print(f"-------------Running hl.identity_by_descent() -----------")
            relatedness_ht = hl.identity_by_descent(mt, maf=mt.info.AF[0])
            if args.overwrite or (
                    not hl.hadoop_exists(f'{get_aou_relatedness_path(extension="ht")}')
            ):
                print(f"-------------Writing AoU IBD HT -----------")
                relatedness_ht.write(
                    get_aou_relatedness_path(extension="ibd.ht"), args.overwrite
                )
        if args.relatedness_type == "king":
            print(f"-------------Running hl.king() -----------")
            relatedness_mt = hl.king(mt.GT)
            if args.overwrite or (
                    not hl.hadoop_exists(f'{get_aou_relatedness_path(extension="king.ht")}')
            ):
                print(f"-------------Writing AoU King relatedness MT -----------")
                relatedness_mt.write(
                    get_aou_relatedness_path(extension="king.mt"), args.overwrite
                )
    if args.run_sample_qc:
        if not hl.hadoop_exists(get_aou_util_path('mt_sample_qc')):
            print('Run sample qc MT.....')
            mt = hl.read_matrix_table(ACAF_MT_PATH)
            mt = mt.filter_rows(mt.locus.in_autosome())
            # mt = mt.filter_rows(mt.locus.contig == 'chr1')
            ht = hl.sample_qc(mt, name='mt_sample_qc')
            ht.write(get_aou_util_path('mt_sample_qc'), overwrite=args.overwrite)
        mt = hl.read_matrix_table(get_aou_util_path('mt_sample_qc'))
        ht = mt.cols()
        print('Write sample QC HT.....')
        ht.write(get_aou_util_path('acaf_mt_sample_qc'))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--single_variant_only", help="Run single variant SAIGE", action="store_true"
    )
    parser.add_argument(
        "--test",
        help="Whether to run test on a subset of the data (e.g. chr1)",
        action="store_true",
    )
    parser.add_argument(
        "--dataproc",
        help="Whether to run the pipeline on dataproc",
        action="store_true",
    )
    parser.add_argument("--overwrite", help="Overwrite everything", action="store_true")
    parser.add_argument(
        "--omit_ld_prune", help="Overwrite everything", action="store_true"
    )
    parser.add_argument(
        "--pops",
        help="Comma-separated list of pops to run",
        default="afr,amr,eas,eur,mid,sas",
    )
    parser.add_argument(
        "--create_plink_file", help="Overwrite everything", action="store_true"
    )
    parser.add_argument(
        "--create_gene_mapping_file", help="Overwrite everything", action="store_true"
    )
    parser.add_argument(
        "--run-relatedness",
        help="Compute relatedness information",
        action="store_true",
    )
    parser.add_argument(
        "--get-duplicated-samples",
        help="Get duplicated samples from the original relatedness kinship table",
        action="store_true",
    )
    parser.add_argument(
        "--relatedness-type",
        help="What function to use for running relatedness",
        choices=["ibd", "king"],
    )
    parser.add_argument(
        "--run-sample-qc",
        help="Whether to run sample QC on the ACAF MT",
        action="store_true",
    )
    parser.add_argument("--run_vep", help="Run Vep on the VAT", action="store_true")
    args = parser.parse_args()
    main(args)
