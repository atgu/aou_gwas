#!/usr/bin/env python3

__author__ = "wlu"

import argparse

# from aou_gwas import *  # for dataproc
from utils.utils import *  # for QoB
from utils.resources import *  # for QoB
from gnomad.utils.vep import *
from gnomad.utils.filtering import *
from ukbb_common import create_gene_map_ht, post_process_gene_map_ht

# from ukbb_pan_ancestry import *


# ALL PATHs: https://support.researchallofus.org/hc/en-us/articles/4616869437204-Controlled-CDR-Directory
# Descriptions: https://support.researchallofus.org/hc/en-us/articles/4614687617556-How-the-All-of-Us-Genomic-data-are-organized
GENE_INTERVAL_PATH = (
    f"gs://aou_wlu/data/group_positions_{N_GENE_PER_GROUP}_protein_coding.ht"
)


def mac_category_case_builder(call_stats_ac_expr, call_stats_af_expr):
    return (
        hl.case()
        .when(call_stats_ac_expr <= 5, call_stats_ac_expr)
        .when(call_stats_ac_expr <= 10, 10)
        .when(call_stats_ac_expr <= 20, 20)
        .when(call_stats_af_expr <= 0.001, 0.001)
        .when(call_stats_af_expr <= 0.01, 0.01)
        .when(call_stats_af_expr <= 0.1, 0.1)
        .default(0.99)
    )


def filter_ht_for_plink(
    ht: hl.Table,
    pop: str,
    min_call_rate: float = 0.95,
    variants_per_mac_category: int = 2000,
    variants_per_maf_category: int = 10000,
):
    ht = ht.filter(
        (ht.locus.in_autosome())
        & (ht.info.AN >= N_SAMPLES[pop] * 2 * min_call_rate)
        & (ht.info.AC[0] > 0)
    )
    ht = ht.annotate(
        mac_category=mac_category_case_builder(ht.info.AC[0], ht.info.AF[0])
    )

    # From: https://hail.zulipchat.com/#narrow/stream/123010-Hail-Query-0.2E2-support/topic/.E2.9C.94.20randomly.20sample.20table/near/388162012
    bins = ht.aggregate(hl.agg.collect_as_set(ht.mac_category))
    ac_bins = [bin for bin in bins if bin >= 1]
    print(ac_bins)
    af_bins = [bin for bin in bins if bin < 1]
    print(af_bins)

    binned_samples_ac = ht.aggregate(
        hl.agg.array_agg(
            lambda x: hl.agg.filter(
                ht.mac_category == x,
                hl.agg._reservoir_sample(ht.key, variants_per_mac_category),
            ),
            hl.literal(ac_bins),
        )
    )
    binned_samples_af = ht.aggregate(
        hl.agg.array_agg(
            lambda x: hl.agg.filter(
                ht.mac_category == x,
                hl.agg._reservoir_sample(ht.key, variants_per_maf_category),
            ),
            hl.literal(af_bins),
        ),
    )
    binned_samples = binned_samples_ac + binned_samples_af

    samples = [sample for bin in binned_samples for sample in bin]
    ht = hl.Table.parallelize(samples).key_by(*ht.key.keys())

    return ht


def main(args):
    hl.init(default_reference="GRCh38")

    # Load info tables
    vat_ht = hl.read_table(get_aou_util_path(name="vat"))
    vat_ht = vat_ht.collect_by_key()

    if args.run_vep:
        # Vep resource from gnomAD:
        # https://github.com/broadinstitute/gnomad_methods/blob/e2f2b055a2d40991851e7852dea435b4eaef45ea/gnomad/resources/grch38/reference_data.py#L100
        if (not hl.hadoop_exists(get_aou_util_path(name="vep_corrected"))) or args.overwrite:
            # vep_ht = vep_or_lookup_vep(vat_ht, vep_version="vep105") # with init file: gs://gcp-public-data--gnomad/resources/vep/v105/vep105-init.sh
            vep_ht = hl.vep(vat_ht, csq=True)  # with init file: gs://aou_analysis/vep105-init_mod.sh
            vep_ht.write(
                get_aou_util_path(name="vep_corrected"),
                overwrite=args.overwrite,
            )
        vep_ht = hl.read_table(get_aou_util_path(name="vep_corrected"))
        vep_ht.describe()

    if args.create_gene_mapping_files:  # TODO: run this once VEP table is ready
        ht = hl.read_table(get_aou_util_path(name="vep_corrected"))
        qual_ht = hl.read_table(
            get_ukb_exomes_qual_ht_path(CURRENT_TRANCHE)
        )  # TODO: find the correct quality filter
        qual_ht = qual_ht.filter(hl.len(qual_ht.filters) == 0)
        call_stats_ht = hl.read_table(get_aou_util_path(name="vat"))
        call_stats_ht = call_stats_ht.collect_by_key()
        pops = args.pops.split(",") if (args.pops is not None) else POPS
        for pop in pops:
            max_an = call_stats_ht.aggregate(
                hl.struct(
                    autosomes=hl.agg.max(call_stats_ht.values[f"gvs_{pop}_an"]),
                    x=hl.agg.filter(
                        call_stats_ht.locus.in_x_nonpar(),
                        hl.agg.max(call_stats_ht.values[f"gvs_{pop}_an"]),
                        y=hl.agg.filter(
                            call_stats_ht.locus.in_y_nonpar(),
                            hl.agg.max(call_stats_ht.values[f"gvs_{pop}_an"]),
                        ),
                    ),
                )
            )
            an = call_stats_ht.values[f"gvs_{pop}_an"]
            call_stats_ht = call_stats_ht.filter(
                hl.case()
                .when(call_stats_ht.locus.in_x_nonpar(), an > 0.8 * max_an.x)
                .when(call_stats_ht.locus.in_y_nonpar(), an > 0.8 * max_an.y)
                .default(an > 0.8 * max_an.autosomes)
            )
            ht = ht.annotate(freq=call_stats_ht[ht.key].values[f"gvs_{pop}_af"][0])
            ht = ht.filter(hl.is_defined(qual_ht[ht.key]) & hl.is_defined(ht.freq))
            gene_map_ht = create_gene_map_ht(ht, freq_field=ht.freq)
            gene_map_ht.write(
                get_aou_util_path(name="gene_map", parsed=False), args.overwrite
            )

            gene_map_ht = hl.read_table(
                get_aou_util_path(name="gene_map", parsed=False)
            )
            gene_map_ht = post_process_gene_map_ht(gene_map_ht, freq_cutoff=0.01)
            gene_map_ht.write(
                get_aou_util_path(name="gene_map", parsed=True), args.overwrite
            )

    if args.create_plink_file:
        pop_ht = hl.read_table(get_aou_util_path(name="ancestry_preds", parsed=True))
        # {'afr': 56913, 'amr': 45035, 'eas': 5706, 'eur': 133581, 'mid': 942, 'sas': 3217} - Total: 245394
        analysis_type = "variant" if args.single_variant_only else "gene"
        if not args.dataproc:
            if analysis_type == "variant":
                mt = hl.read_matrix_table(ACAF_MT_PATH)
            else:
                mt = hl.read_matrix_table(EXOME_MT_PATH)

            if args.test:
                mt = mt.filter_rows(mt.locus.contig == "chr1")

        pops = args.pops.split(",") if (args.pops is not None) else POPS
        for pop in pops:
            # Note: 1e7 LD pruning for EUR was run with r2=0.05, and chr8 inversion and HLA were removed
            window = "1e7" if pop == "eur" else "1e6"
            r2 = 0.05 if (window == "1e7") and (analysis_type == "variant") else 0.1
            iteration = 1
            pop_mt_path = get_aou_saige_plink_path(
                analysis_type=analysis_type,
                pop=pop,
                extension="mt",
                pruned=False,
                data_iteration=iteration,
                test=args.test,
            )
            pop_ld_ht_path = get_aou_saige_plink_path(
                analysis_type=analysis_type,
                pop=pop,
                extension="ht",
                pruned=True,
                window_size=window,
                test=args.test,
            )
            pop_plink_path = get_aou_saige_plink_path(
                analysis_type=analysis_type,
                pop=pop,
                extension="plink",
                pruned=True,
                data_iteration=iteration,
                window_size=window,
                test=args.test,
            )
            if not hl.hadoop_exists(f"{pop_mt_path}/_SUCCESS"):
                # Filter samples
                sub_pop_ht = pop_ht.filter(pop_ht.ancestry_pred == pop)
                sub_mt = mt.filter_cols(hl.is_defined(sub_pop_ht[mt.col_key]))
                n_samples = sub_mt.count_cols()
                print(f"Got {n_samples} samples for {pop}...")

                # Filter variants
                # sub_mt = sub_mt.filter_rows(sub_mt.locus.in_autosome()) # TODO: double check that we want to keep sex chromosomes
                if analysis_type == "variant":
                    sub_mt = sub_mt.filter_rows(sub_mt.info.AF[0] >= 0.01)
                else:
                    sampled_variants_ht = filter_ht_for_plink(ht=sub_mt.rows(), pop=pop)
                    sampled_variants_ht = sampled_variants_ht.naive_coalesce(1000).checkpoint(
                        get_aou_saige_plink_path(
                            analysis_type=analysis_type,
                            pop=pop,
                            extension="ht",
                            pruned=False,
                            window_size="1e6",
                            test=args.test,
                        ),
                        _read_if_exists=not args.overwrite,
                        overwrite=args.overwrite,
                    )
                    sub_mt = sub_mt.filter_rows(
                        hl.is_defined(sampled_variants_ht[sub_mt.row_key])
                    )

                if pop == "eur" and analysis_type == "variant":
                    print("Removing HLA...")
                    # Common inversion taken from Table S4 of https://www.ncbi.nlm.nih.gov/pubmed/27472961
                    # (converted to GRCh38 by: https://liftover.broadinstitute.org/#input=chr8%3A8055789-11980649&hg=hg19-to-hg38 )
                    # Also removing HLA, from https://www.ncbi.nlm.nih.gov/grc/human/regions/MHC?asm=GRCh38
                    sub_mt = sub_mt.filter_rows(
                        ~hl.parse_locus_interval(
                            "chr8:8198267-12123140", reference_genome="GRCh38"
                        ).contains(sub_mt.locus)
                        & ~hl.parse_locus_interval(
                            "chr6:28510120-33480577", reference_genome="GRCh38"
                        ).contains(sub_mt.locus)
                    )
                print(f"Exporting {pop.upper()} MT to {pop_mt_path}")
                sub_mt = sub_mt.naive_coalesce(1000).checkpoint(
                    pop_mt_path,
                    _read_if_exists=not args.overwrite,
                    overwrite=args.overwrite,
                )
            sub_mt = hl.read_matrix_table(pop_mt_path)

            if (not args.omit_ld_prune) and (args.dataproc):
                # This has to be run on dataproc, QoB hits transient error:
                # Error summary: GoogleJsonResponseException: 404 Not Found @ hail==0.2.124
                print(f"LD-pruning to {pop_ld_ht_path}...")
                ht = hl.ld_prune(
                    sub_mt.GT,
                    r2=float(r2),
                    bp_window_size=int(float(window)),
                    block_size=1024,
                )
                ht.write(pop_ld_ht_path, overwrite=not args.omit_ld_prune)

            if hl.hadoop_exists(f"{pop_ld_ht_path}/_SUCCESS"):
                sub_mt = sub_mt.unfilter_entries()
                ht = hl.read_table(pop_ld_ht_path)
                sub_mt = sub_mt.filter_rows(hl.is_defined(ht[sub_mt.row_key]))

                ## TODO: check numbers of variants per population
                ## if pop == 'EUR':
                ## sub_mt = sub_mt.filter_rows(hl.rand_bool(0.55))

                if args.overwrite or not hl.hadoop_exists(f"{pop_plink_path}.bed"):
                    print(f"Exporting plink to {pop_plink_path}...")
                    hl.export_plink(sub_mt, pop_plink_path)

                # TODO: residualize GRM by first 20 PCs
                if args.overwrite or not hl.hadoop_exists(
                    get_aou_samples_file_path(
                        analysis_type=analysis_type, pop=pop, data_iteration=iteration
                    )
                ):
                    with hl.hadoop_open(
                        get_aou_samples_file_path(
                            analysis_type=analysis_type,
                            pop=pop,
                            data_iteration=iteration,
                        ), "w",) as f:
                        f.write("\n".join(sub_mt.s.collect()) + "\n")


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
    parser.add_argument("--run_vep", help="Run Vep on the VAT", action="store_true")
    args = parser.parse_args()
    main(args)
