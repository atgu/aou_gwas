#!/usr/bin/env python3

import sys
import os
import ast
import argparse
import json
import hail as hl
import hailtop.fs as hfs
import gcsfs
import pickle
import pandas as pd
import joblib

from pyspark.sql import SparkSession

# Global constants
MY_BUCKET = 'gs://aou_amc'
TMP_BUCKET = 'gs://aou_tmp'
TRANCHE = "v8"
ANALYSIS_BUCKET = "gs://aou_amc_analyses"
EXTERNAL_ANALYSIS_BUCKET = "gs://aou_analysis/"
DATA_PATH = f"{ANALYSIS_BUCKET}/data"
EXTERNAL_DATA_PATH = f"gs://aou_analysis/{TRANCHE}/data"
SNPINDEL_OUT_PATH = f'{EXTERNAL_DATA_PATH}/vep/aou_{TRANCHE}_vep_full.ht'
SNPINDEL_AMC_PATH = f'{ANALYSIS_BUCKET}/data/vep/aou_{TRANCHE}_vep_full.ht'

BRAVA_PATH = 'gs://aou_amc_analyses/data/utils/brava_annot/brava.ht'

ALL_MIS_VSM_SCALLION = f"{ANALYSIS_BUCKET}/data/utils/missense_predictions/all_missense_w_predictions_w_pct.parquet"

def load_predictions_missense(path: str) -> hl.Table:
    """Read a predictions parquet file into an annotated, keyed Hail Table,
    with all '_pct' columns grouped into a 'preds_missense' struct.
    Parameters
    ----------
    path : str
        GCS (or other) path to the predictions parquet file, e.g.
    """
    spark = SparkSession.builder.getOrCreate()
    spark_df = spark.read.parquet(path)
    pred_ht = hl.Table.from_spark(spark_df)

    pred_ht = pred_ht.annotate(
        locus=hl.locus(pred_ht.chrom, hl.int32(pred_ht.pos), reference_genome='GRCh38'),
        alleles=hl.array([pred_ht.ref, pred_ht.alt]),
    )

    pred_ht = pred_ht.key_by('locus', 'alleles', 'ensg')

    pct_cols = [f for f in pred_ht.row.dtype.fields if f.endswith('_pct')]
    pred_ht = pred_ht.annotate(preds_missense=hl.struct(**{c: pred_ht[c] for c in pct_cols}))
    pred_ht = pred_ht.drop(*pct_cols)

    ht_path = path.replace('.parquet', '.ht')
    print(f"Checkpointing missense predictions table to {ht_path}...")
    pred_ht = pred_ht.checkpoint(ht_path)

    return pred_ht

def create_brava_ht(overwrite: bool = False) -> hl.Table:
    """
    Build the BRaVa annotation Hail Table: import the per-chromosome
    variant files, union them, key by locus/alleles, and write the
    result to BRAVA_PATH.
    """
    chromosomes = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10',
                   '11', '12', '13', '14', '15', '16', '17', '18', '19', '20',
                   '21', '22', 'X']

    tables = []
    for chrom in chromosomes:
        file_path = f'{DATA_PATH}/utils/brava_annot/brava_split_chr/aou.v8.chr{chrom}.variants_only.spliceai=0.20_cadd=28.1_revel=0.773.canonical.txt.gz'
        ht_chr = hl.import_table(
            file_path,
            delimiter='\t',
            impute=True,
            comment='#',
            force=True  # Required for regular gzip files
        )
        # Optionally add chromosome info if not in the data
        ht_chr = ht_chr.annotate(chromosome=chrom)
        tables.append(ht_chr)

    ht_combined = tables[0]
    for table in tables[1:]:
        ht_combined = ht_combined.union(table)

    ht_combined = ht_combined.annotate(
        id_parts=ht_combined.ID.split(':')
    )
    ht_combined = ht_combined.annotate(
        locus=hl.locus(
            ht_combined.id_parts[0],
            hl.int32(ht_combined.id_parts[1])
        ),
        alleles=hl.array([
            ht_combined.id_parts[2],
            ht_combined.id_parts[3]
        ])
    )
    ht_combined = ht_combined.drop('id_parts')
    ht_combined = ht_combined.key_by('locus', 'alleles')

    print("Table schema after setting keys:")
    print(ht_combined.describe())

    ht_combined = ht_combined.naive_coalesce(500)

    print(f"Writing BRaVa table to {BRAVA_PATH}...")
    ht_combined = ht_combined.checkpoint(BRAVA_PATH, overwrite=overwrite)

    return ht_combined

def snp_indel_vep_concat(output_path, overwrite=False):
    from gnomad.utils.vep import process_consequences
    snp_vep_path = 'gs://gcp-public-data--gnomad/resources/context/grch38_context_vep_annotated.v105.ht'
    snp_vep_ht = hl.read_table(snp_vep_path).naive_coalesce(3000)
    shared_fields_to_drop = ['uniparc', 'trembl', 'swissprot']

    print("Processing SNP VEP annotations...")
    snp_ht = snp_vep_ht.annotate(vep=snp_vep_ht.vep.drop('context'))
    snp_ht = snp_ht.annotate(
        vep=snp_ht.vep.annotate(
            intergenic_consequences=snp_ht.vep.intergenic_consequences.map(lambda x: x.drop('minimised')),
            motif_feature_consequences=snp_ht.vep.motif_feature_consequences.map(lambda x: x.drop('minimised')),
            regulatory_feature_consequences=snp_ht.vep.regulatory_feature_consequences.map(
                lambda x: x.drop('minimised')),
            transcript_consequences=snp_ht.vep.transcript_consequences.map(
                lambda x: x.drop(*shared_fields_to_drop, 'minimised')),
        ))

    indel_vep_path = f'{EXTERNAL_DATA_PATH}/vep/aou_vds_variant_data_row_{TRANCHE}_vep.ht'
    indel_vep_ht = hl.read_table(indel_vep_path)
    indel_fields_to_drop = ['ancestral', 'context']
    print("Processing INDEL VEP annotations...")
    indel_ht = indel_vep_ht.annotate(vep=indel_vep_ht.vep.drop('minimised'))
    indel_ht = indel_ht.annotate(
        vep=indel_ht.vep.annotate(
            intergenic_consequences=indel_ht.vep.intergenic_consequences.map(
                lambda x: x.drop(*indel_fields_to_drop)),
            motif_feature_consequences=indel_ht.vep.motif_feature_consequences.map(
                lambda x: x.drop(*indel_fields_to_drop)),
            regulatory_feature_consequences=indel_ht.vep.regulatory_feature_consequences.map(
                lambda x: x.drop(*indel_fields_to_drop)),
            transcript_consequences=indel_ht.vep.transcript_consequences.map(
                lambda x: x.drop(*shared_fields_to_drop, *indel_fields_to_drop)),
        ))
    print("Merging SNP and INDEL VEP annotations...")
    
    merged_ht = snp_ht.union(indel_ht, unify=True)
    
    process_vep_ht = process_consequences(merged_ht)
    merged_ht = merged_ht.annotate(
        worst_csq_by_gene_canonical=process_vep_ht[merged_ht.key].vep.worst_csq_by_gene_canonical
    )
    
    merged_ht = merged_ht.naive_coalesce(5000)
    print(f"Writing merged VEP table to {SNPINDEL_OUT_PATH}...")
    merged_ht = merged_ht.checkpoint(SNPINDEL_OUT_PATH, overwrite=overwrite)
    
    if overwrite:
        merged_ht = merged_ht.checkpoint(output_path, overwrite=overwrite)
    
    return merged_ht

def create_raw_gene_map(pop: str, annot_type: str, overwrite: bool = False, overwrite_report: bool = False, overwrite_context: bool = False):
    """
    Create raw gene mapping file for a specific ancestry population.

    Args:
        pop: Ancestry population code (e.g., 'AFR', 'EUR')
        annot_type: Annotation type ('snp_indel' or 'brava')
        overwrite: Whether to overwrite existing files
        overwrite_report: Whether to overwrite the snpindel/missingness report files
        overwrite_context: Whether to regenerate the merged snp+indel VEP context table
    """
    from annotations import create_gene_map_ht

    if not hl.hadoop_exists(SNPINDEL_OUT_PATH) or overwrite_context:
        print(f"SNP+indel VEP table not found at {SNPINDEL_OUT_PATH} (or overwrite_context=True), generating it...")
        snp_indel_vep_ht = snp_indel_vep_concat(SNPINDEL_AMC_PATH, overwrite=overwrite_context)
    else:
        snp_indel_vep_ht = hl.read_table(SNPINDEL_OUT_PATH)
    snp_indel_vep_ht = snp_indel_vep_ht.key_by('locus', 'alleles')

    gene_map_subdir = 'brava' if annot_type == 'brava' else 'gnomad_context'
    gene_map_ht_path = f"{DATA_PATH}/utils/gene_map/{gene_map_subdir}/aou_{pop.upper()}_gene_map_{TRANCHE}.ht"
    if not overwrite and hl.hadoop_exists(gene_map_ht_path):
        print(f"Raw gene map file already exists for {pop.upper()} and overwrite=False. Skipping creation.")
        return

    report_dir = f"{DATA_PATH}/utils/gene_map/{gene_map_subdir}/report"
    report_prefix = f"{report_dir}/aou_{pop.upper()}_{TRANCHE}"
    
    call_stats_ht_path = f"{EXTERNAL_DATA_PATH}/utils/call_stats/exome_pruned/{pop.upper()}_exome_call_stats.ht"
    print(f"Loading call stats from {call_stats_ht_path}...")
    call_stats_ht = hl.read_table(call_stats_ht_path)
    call_stats_ht = call_stats_ht.filter(call_stats_ht.call_stats.AC[1] > 0)
    
    print(f'---------Generating raw gene mapping HT ({pop.upper()})-----------------')
    max_an = call_stats_ht.aggregate(
        hl.struct(
            autosomes=hl.agg.max(call_stats_ht.call_stats.AN),
            x=hl.agg.filter(
                call_stats_ht.locus.in_x_nonpar(),
                hl.agg.max(call_stats_ht.call_stats.AN)
            ),
            y=hl.agg.filter(
                call_stats_ht.locus.in_y_nonpar(),
                hl.agg.max(call_stats_ht.call_stats.AN),
            ),
        ),
    )
    
    an = call_stats_ht.call_stats.AN
    call_stats_ht = call_stats_ht.filter(
        hl.case()
        .when(call_stats_ht.locus.in_x_nonpar(), an > 0.8 * max_an.x)
        .when(call_stats_ht.locus.in_y_nonpar(), an > 0.8 * max_an.y)
        .default(an > 0.8 * max_an.autosomes)
    )
    
    snp_indel_vep_ht = snp_indel_vep_ht.annotate(
        freq = call_stats_ht[snp_indel_vep_ht.key].call_stats.AF[1]
    )

    print(f'count before handling call stats for {pop.upper()} {snp_indel_vep_ht.count()}')
    snp_indel_vep_ht = snp_indel_vep_ht.filter(
        hl.is_defined(snp_indel_vep_ht.freq)
    )

    snp_indel_vep_ht = snp_indel_vep_ht.explode(snp_indel_vep_ht.worst_csq_by_gene_canonical)
    snp_indel_vep_ht = snp_indel_vep_ht.annotate(
        ensg=snp_indel_vep_ht.worst_csq_by_gene_canonical.gene_id
    )
    snp_indel_vep_ht = snp_indel_vep_ht.filter(
        snp_indel_vep_ht.worst_csq_by_gene_canonical.gene_id.startswith('ENSG')
    )
    snp_indel_vep_ht = snp_indel_vep_ht.key_by('locus', 'alleles', 'ensg')

    tmp_filtered_vep_path = f"{DATA_PATH}/utils/gene_map/tmp/{gene_map_subdir}/aou_{pop.upper()}_snp_indel_vep_filtered_{TRANCHE}.ht"
    print(f"Checkpointing call-stats-filtered VEP table to {tmp_filtered_vep_path}...")
    snp_indel_vep_ht = snp_indel_vep_ht.checkpoint(
        tmp_filtered_vep_path, overwrite=overwrite, _read_if_exists=not overwrite
    )

    # Get missense variant weights
    all_mis_preds_ht_path = ALL_MIS_VSM_SCALLION.replace('.parquet', '.ht')
    if hl.hadoop_exists(all_mis_preds_ht_path):
        print(f"Missense predictions table already exists at {all_mis_preds_ht_path}, loading it directly.")
        all_mis_preds = hl.read_table(all_mis_preds_ht_path)
    else:
        all_mis_preds = load_predictions_missense(ALL_MIS_VSM_SCALLION)
    snp_indel_vep_ht = snp_indel_vep_ht.annotate(
        preds_missense=all_mis_preds[snp_indel_vep_ht.key].preds_missense
    )
    
    print(f'count after handling call stats for {pop.upper()} {snp_indel_vep_ht.count()}')
    if annot_type == 'brava':
        if not hl.hadoop_exists(BRAVA_PATH) or overwrite:
            print(f"BRaVa table not found at {BRAVA_PATH} (or overwrite=True), generating it...")
            brava_ht = create_brava_ht(overwrite=overwrite)
        else:
            brava_ht = hl.read_table(BRAVA_PATH)

        brava_ht = brava_ht.rename({'GENE': 'ensg'})
        brava_ht = brava_ht.key_by('locus', 'alleles', 'ensg')

        snp_indel_vep_ht = snp_indel_vep_ht.annotate(
            brava=brava_ht[snp_indel_vep_ht.key]
        )
    
    gene_map_ht = create_gene_map_ht(
        snp_indel_vep_ht, annot_type, freq_field='freq',
        report_prefix=report_prefix,
        overwrite_report=overwrite_report,
    )
    print(f'---------Exporting raw gene mapping HT ({pop.upper()})-----------------')
    gene_map_ht.checkpoint(gene_map_ht_path, overwrite=overwrite)
    print(f"Raw gene map for {pop.upper()} saved to {gene_map_ht_path}")
    
    return gene_map_ht_path

def process_gene_map(pop: str, annot_type: str, overwrite: bool = False):
    """
    Process gene mapping file for a specific ancestry population.

    Args:
        pop: Ancestry population code (e.g., 'AFR', 'EUR')
        annot_type: Annotation type ('snp_indel' or 'brava')
        overwrite: Whether to overwrite existing files
    """
    from annotations import post_process_gene_map_ht

    gene_map_subdir = 'brava' if annot_type == 'brava' else 'gnomad_context'
    gene_map_ht_path = f"{DATA_PATH}/utils/gene_map/{gene_map_subdir}/aou_{pop.upper()}_gene_map_{TRANCHE}.ht"
    processed_gene_map_ht_path = f"{DATA_PATH}/utils/gene_map/{gene_map_subdir}/aou_{pop.upper()}_gene_map_processed_{TRANCHE}.ht"

    if not overwrite and hl.hadoop_exists(processed_gene_map_ht_path):
        print(f"Processed gene map file already exists for {pop.upper()} and overwrite=False. Skipping processing.")
        return
    
    if not hl.hadoop_exists(gene_map_ht_path):
        raise FileNotFoundError(f"Raw gene map file does not exist for {pop.upper()}. Run create_raw_gene_map first.")
    
    print(f"Loading raw gene map from {gene_map_ht_path}...")
    gene_map_ht = hl.read_table(gene_map_ht_path)


    gene_map_ht = post_process_gene_map_ht(gene_map_ht, freq_cutoff=0.01, annot_type=annot_type)
    print(f'---------Adding VSM weights and SCALLION predictions({pop.upper()})-----------------')
    # gene_map_ht = add_top_decile_annotations(gene_map_ht)
    print(f'---------Exporting processed gene mapping HT ({pop.upper()})-----------------')
    gene_map_ht = gene_map_ht.checkpoint(processed_gene_map_ht_path, overwrite=overwrite)
    gene_map_ht.describe()
    gene_map_ht.show()
    print(f'Completed processing for ancestry: {pop.upper()}')
    return processed_gene_map_ht_path


def main(args):
    hl.init(
        tmp_dir=TMP_BUCKET,
        gcs_requester_pays_configuration="aou-neale-gwas",
        default_reference="GRCh38",
        log=f"/gene_map_generation_{TRANCHE}.log",
    )

    try:
        for pop in args.ancestries:
            print(f"Processing ancestry {pop}")
            if not args.process_only:
                create_raw_gene_map(pop, args.annotation_type, args.overwrite, args.overwrite_report, args.overwrite_context)
            process_gene_map(pop, args.annotation_type, args.overwrite)
    finally:
        from datetime import date
        hl.copy_log(f"{MY_BUCKET}/pipeline_{TRANCHE}_{date.today()}.log")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="AoU VEP and Gene Map Processing Pipeline")
    parser.add_argument("--overwrite",         help="Overwrite existing raw/processed gene map files",                                     action="store_true")
    parser.add_argument("--overwrite-report",  help="Overwrite existing report files (snpindel formatted HT + missingness report)",         action="store_true")
    parser.add_argument("--overwrite-context", help="Regenerate the merged snp+indel VEP context table even if it already exists",          action="store_true")
    parser.add_argument("--process-only",      help="Skip raw gene map creation and go straight to processing (raw gene map must already exist)", action="store_true")
    parser.add_argument("--annotation-type",   help="Type of annotation", choices=["brava", "snp_indel"], default="snp_indel", type=str)
    parser.add_argument("--ancestries",        help="Comma-separated ancestries to process (e.g., 'EUR,AFR,AMR')", type=lambda s: s.split(","))

    args = parser.parse_args()

    if not args.ancestries:
        args.ancestries = ["EUR", "AFR", "AMR", "EAS", "SAS", "MID"]
        print(f"No ancestries specified, using defaults: {', '.join(args.ancestries)}")

    main(args)