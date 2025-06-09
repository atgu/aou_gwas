#!/usr/bin/env python3

import sys
import os
import argparse
import hail as hl
import hailtop.batch as hb
import hailtop.fs as hfs

# Global constants
MY_BUCKET = 'gs://aou_amc'
TMP_BUCKET = 'gs://aou_tmp'

TRANCHE = "v8"
ANALYSIS_BUCKET = "gs://aou_amc_analyses/"
EXTERNAL_ANALYSIS_BUCKET = "gs://aou_analysis/"
DATA_PATH = f"{ANALYSIS_BUCKET}/brava_annot/data"
EXTERNAL_DATA_PATH = f"gs://aou_analysis/{TRANCHE}/data"

SNPINDEL_OUT_PATH = f'{DATA_PATH}/aou_{TRANCHE}_vep_full.ht'
SNP_BRAVA_VAT_OUT = f'{DATA_PATH}/aou_{TRANCHE}_snp_vep_brava.ht'


def initialize_hail(batch_mode=False, log_file="/hail_operation.log", app_name=None):
    init_params = {
        "tmp_dir": TMP_BUCKET,
        "gcs_requester_pays_configuration": "aou-neale-gwas",
        "default_reference": "GRCh38",
        "log": log_file
    }
    
    if app_name:
        init_params["app_name"] = app_name
        
    if batch_mode:
        init_params.update({
            "master": "local[32]",
            "worker_cores": 8,
            "worker_memory": "highmem"
        })
    
    hl.init(**init_params)


def merge_data(merge_type='snp_indel', overwrite=False, batch_mode=False):
    """
    Unified function to merge data based on merge_type
    
    Args:
        merge_type: 'snp_indel' for SNP-INDEL union or 'brava' for BRAVA-SNP-VAT intersection
        overwrite: Whether to overwrite existing files
        batch_mode: Whether running in batch mode
    """
    
    if batch_mode:
        log_file = f"/merge_{merge_type}_ht.log"
        initialize_hail(batch_mode=batch_mode, log_file=log_file)

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
    
    if merge_type == 'snp_indel':
        indel_vep_path = f'{EXTERNAL_ANALYSIS_BUCKET}/{TRANCHE}/data/vep/aou_vds_variant_data_row_{TRANCHE}_vep.ht'
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
        output_path = SNPINDEL_OUT_PATH
        
    elif merge_type == 'brava':
        brava_path = 'gs://aou_amc_analyses/brava_annot/data/brava_annot_processed_v2.ht'
        brava_ht = hl.read_table(brava_path).naive_coalesce(2000)
        vat_aou_path = 'gs://aou_analysis/v8/data/vat/aou_PARSED_SORTED_COLLECTED_vat_v8.ht/'
        vat_aou = hl.read_table(vat_aou_path).naive_coalesce(3000)
        
        print("Performing three-way intersection: BRAVA, SNP, and VAT AoU...")
        brava_snp_ht = brava_ht.join(snp_ht, how='inner')
        print(f"STEP1 - BRAVA-SNP intersection complete. Variants: {brava_snp_ht.count()}")
        
        print("Processing VAT data...")
        fields_to_keep = ['revel'] + [
            f for f in vat_aou.values.dtype.element_type.fields if f.startswith('splice_ai')
            ]
        sub_vat = vat_aou.annotate(
            aou_vat_annot=hl.struct(**{f: vat_aou.values[0][f] for f in fields_to_keep})
        )
        sub_vat = sub_vat.checkpoint(f'{DATA_PATH}/variant_annot_processed.ht', _read_if_exists=True)
        print("VAT data processing complete...")

        print("STEP2 - Adding VAT AoU data to the intersection...")
        merged_ht = brava_snp_ht.join(sub_vat, how='inner')
        print(f"STEP2 - Final three-way intersection complete. Variants: {merged_ht.count()}")
        output_path = SNP_BRAVA_VAT_OUT
    else:
        raise ValueError(f"Invalid merge_type: {merge_type}. Must be 'brava' or 'snp_indel'")
    
    print("Adding process_consequences annotations...")
    process_vep_ht = process_consequences(merged_ht)
    merged_ht = merged_ht.annotate(
        worst_csq_by_gene_canonical=process_vep_ht[merged_ht.key].vep.worst_csq_by_gene_canonical
    )

    if merge_type == 'brava':
        merged_ht = merged_ht.annotate(
            values=merged_ht.values.map(
                lambda val: val.annotate(
                    worst_csq_by_gene_canonical=hl.find(
                        lambda x: x.transcript_id == val.TRANSCRIPT,
                        merged_ht.worst_csq_by_gene_canonical))))
        merged_ht = merged_ht.drop('worst_csq_by_gene_canonical')

    print(f"Writing merged table to {output_path}...")
    merged_ht = merged_ht.naive_coalesce(2500)
    merged_ht = merged_ht.checkpoint(output_path, overwrite=overwrite)
    
    return merged_ht


def create_raw_gene_map(pop: str, annot_type: str, overwrite: bool = False, batch_mode=False):
    """
    Create raw gene mapping file for a specific ancestry population.
    
    Args:
        pop: Ancestry population code (e.g., 'AFR', 'EUR')
        annot_type: Annotation type ('snp_indel' or 'brava')
        overwrite: Whether to overwrite existing files
        batch_mode: Whether the function is being run in a batch job
    """
    from annotations import create_gene_map_ht

    if batch_mode:
        initialize_hail(batch_mode=batch_mode, log_file=f"/create_raw_gene_map_{pop}.log")
    
    # Load VEP table
    snp_indel_vep_path = SNPINDEL_OUT_PATH if annot_type == 'snp_indel' else SNP_BRAVA_VAT_OUT

    try:
        snp_indel_vep_ht = hl.read_table(snp_indel_vep_path)
        print(snp_indel_vep_ht.count())
        snp_indel_vep_ht.show()
    except:
        pass
        
    gene_map_ht_path = f"{DATA_PATH}/utils/gene_map/aou_{pop}_gene_map_{TRANCHE}.ht"
    if not overwrite and hl.hadoop_exists(gene_map_ht_path):
        print(f"Raw gene map file already exists for {pop} and overwrite=False. Skipping creation.")
        return
    
    call_stats_ht_path = f"{EXTERNAL_DATA_PATH}/utils/call_stats/exome_pruned/{pop}_exome_call_stats.ht"
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
    snp_indel_vep_ht = snp_indel_vep_ht.filter(
        hl.is_defined(snp_indel_vep_ht.freq)
    )
    
    gene_map_ht = create_gene_map_ht(snp_indel_vep_ht, annot_type, freq_field='freq')
    print(f'---------Exporting raw gene mapping HT ({pop.upper()})-----------------')
    gene_map_ht.checkpoint(gene_map_ht_path, overwrite=overwrite)
    print(f"Raw gene map for {pop} saved to {gene_map_ht_path}")
    return gene_map_ht_path


def process_gene_map(pop: str, overwrite: bool = False, batch_mode=False):
    """
    Process gene mapping file for a specific ancestry population.
    
    Args:
        pop: Ancestry population code (e.g., 'AFR', 'EUR')
        overwrite: Whether to overwrite existing files
        batch_mode: Whether the function is being run in a batch job
    """
    from annotations import post_process_gene_map_ht

    if batch_mode:
        initialize_hail(batch_mode=batch_mode, log_file=f"/process_gene_map_{pop}.log")
    
    gene_map_ht_path = f"{DATA_PATH}/utils/gene_map/aou_{pop}_gene_map_{TRANCHE}.ht"
    processed_gene_map_ht_path = f"{DATA_PATH}/utils/gene_map/aou_{pop}_gene_map_processed_{TRANCHE}.ht"
    
    if not overwrite and hl.hadoop_exists(processed_gene_map_ht_path):
        print(f"Processed gene map file already exists for {pop} and overwrite=False. Skipping processing.")
        return
    
    if not hl.hadoop_exists(gene_map_ht_path):
        raise FileNotFoundError(f"Raw gene map file does not exist for {pop}. Run create_raw_gene_map first.")
    
    print(f"Loading raw gene map from {gene_map_ht_path}...")
    gene_map_ht = hl.read_table(gene_map_ht_path)
    gene_map_ht = post_process_gene_map_ht(gene_map_ht, freq_cutoff=0.01)
    print(f'---------Exporting processed gene mapping HT ({pop.upper()})-----------------')
    gene_map_ht = gene_map_ht.checkpoint(processed_gene_map_ht_path, overwrite=overwrite)
    gene_map_ht.describe()
    gene_map_ht.show()
    print(f'Completed processing for ancestry: {pop}')
    return processed_gene_map_ht_path


def main(args):
    """
    Main function that runs the merged pipeline:
    1. Merge SNP and INDEL VEP annotations OR BRAVA intersection based on annotation type
    2. Create and process gene mapping files for specified ancestries
    
    Args:
        args: Command-line arguments
    """
    try:
        if args.batch:
            backend = hb.ServiceBackend(
                billing_project="all-by-aou",
                remote_tmpdir=TMP_BUCKET,
            )
            b = hb.Batch(
                name=f"aou_{TRANCHE}_pipeline",
                requester_pays_project="aou-neale-gwas",
                default_python_image="amartinezcarrasco/hailgnomad:latest",
                backend=backend,
            )
            
            if args.merge:
                merge_job = b.new_python_job(name=f"Merge {args.annotation_type}")
                merge_job.memory('highmem')
                merge_job.cpu(8)
                merge_job.call(merge_data, args.annotation_type, args.overwrite, True)
                
                if args.gene_map and args.ancestries:
                    for pop in args.ancestries:
                        job = b.new_python_job(name=f"gene_map_{pop}")
                        job.depends_on(merge_job)
                        job.attributes['ancestry'] = pop
                        
                        if not args.skip_raw_gene_map_file:
                            job.call(create_raw_gene_map, pop, args.annotation_type, args.overwrite, True)
                        job.call(process_gene_map, pop, args.overwrite, True)
            
            elif args.gene_map and args.ancestries:
                for pop in args.ancestries:
                    job = b.new_python_job(name=f"gene_map_{pop}")
                    job.attributes['ancestry'] = pop
                    
                    if not args.skip_raw_gene_map_file:
                        job.call(create_raw_gene_map, pop, args.annotation_type, args.overwrite, True)
                    job.call(process_gene_map, pop, args.overwrite, True)
            
            b.run()
            
        else:
            # Run on dataproc cluster or locally with QoB
            initialize_hail(log_file=f"/gene_map_generation_{TRANCHE}.log")
            if args.merge:
                print(f"Running {args.annotation_type} merge...")
                vep_merged_ht = merge_data(args.annotation_type, args.overwrite)
            
            if args.gene_map and args.ancestries:
                for pop in args.ancestries:
                    print(f"Processing ancestry {pop}")
                    if not args.skip_raw_gene_map_file:
                        create_raw_gene_map(pop, args.annotation_type, args.overwrite)
                    process_gene_map(pop, args.overwrite)
    finally:
        if not args.batch:
            from datetime import date
            hl.copy_log(f"{MY_BUCKET}/pipeline_{TRANCHE}_{date.today()}.log")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="AoU VEP and Gene Map Processing Pipeline")
    
    parser.add_argument(
        "--batch",
        help="Run with Hail Batch",
        action="store_true",
    )
    
    parser.add_argument(
        "--merge",
        help="Run merge operation",
        action="store_true",
    )
    
    parser.add_argument(
        "--annotation-type",
        help="Type of annotation to use: 'brava' for BRAVA intersection or 'snp_indel' for SNP-INDEL merge",
        choices=['brava', 'snp_indel'],
        default='snp_indel',
        type=str
    )
    
    parser.add_argument(
        "--gene-map",
        help="Create and process gene mapping files",
        action="store_true",
    )
    
    parser.add_argument(
        "--overwrite",
        help="Overwrite existing files",
        action="store_true"
    )
    
    parser.add_argument(
        "--ancestries",
        help="Comma-separated list of ancestries to process (e.g., 'EUR,AFR,AMR')",
        type=lambda s: s.split(',')
    )
    
    parser.add_argument(
        "--skip-raw-gene-map-file",
        help="Skip creating raw gene map file (only process existing files)",
        action="store_true"
    )
    
    args = parser.parse_args()
    
    if args.gene_map and not args.ancestries:
        args.ancestries = ['EUR', 'AFR', 'AMR', 'EAS', 'SAS', 'MID']
        print(f"No ancestries specified, using defaults: {', '.join(args.ancestries)}")
    
    main(args)