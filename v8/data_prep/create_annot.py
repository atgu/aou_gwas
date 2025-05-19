#!/usr/bin/env python3

import os
import argparse
import hail as hl
import hailtop.batch as hb

TRANCHE = "v8"
ANALYSIS_BUCKET = "gs://aou_analysis"
MY_BUCKET = 'gs://aou_amc'
TMP_BUCKET = 'gs://aou_tmp'
DATA_PATH = f"{ANALYSIS_BUCKET}/{TRANCHE}/data"


def initialize_hail(batch_mode=False, log_file="/hail_operation.log"):
    init_params = {
        "tmp_dir": TMP_BUCKET,
        "gcs_requester_pays_configuration": "aou-neale-gwas",
        "default_reference": "GRCh38",
        "log": log_file
    }
    if batch_mode:
        init_params.update({
            "master": "local[32]",
            "worker_cores": 8,
            "worker_memory": "highmem"
        })
    
    hl.init(**init_params)


def create_raw_gene_map(pop: str, overwrite: bool = False, batch_mode=False):
    """
    Create raw gene mapping file for a specific ancestry population.
    
    Args:
        pop: Ancestry population code (e.g., 'AFR', 'EUR')
        overwrite: Whether to overwrite existing files
    """
    from annotations import create_gene_map_ht

    if batch_mode:
        initialize_hail(batch_mode=batch_mode)
    
    # Load VEP table
    snp_indel_vep_path = 'gs://aou_analysis/v8/data/vep/aou_v8_vep_full.ht'
    try:
        snp_indel_vep_ht = hl.read_table(snp_indel_vep_path)
    except:
        print("Error reading VEP table, retrying with lower parallelism")
        hl._set_flags(use_new_shuffle='1')
        snp_indel_vep_ht = hl.read_table(snp_indel_vep_path)
    
    # Output path
    gene_map_ht_path = f"{DATA_PATH}/utils/gene_map/aou_{pop}_gene_map_{TRANCHE}.ht"
    
    # Check if file exists and skip if not overwriting
    if not overwrite and hl.hadoop_exists(gene_map_ht_path):
        print(f"Raw gene map file already exists for {pop} and overwrite=False. Skipping creation.")
        return
    
    # Load and process call stats
    call_stats_ht_path = f"gs://aou_analysis/v8/data/utils/call_stats/exome_pruned/{pop}_exome_call_stats.ht"
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
    
    gene_map_ht = create_gene_map_ht(snp_indel_vep_ht, freq_field='freq')
    print(f'---------Exporting raw gene mapping HT ({pop.upper()})-----------------')
    gene_map_ht.checkpoint(gene_map_ht_path, overwrite=overwrite)
    print(f"Raw gene map for {pop} saved to {gene_map_ht_path}")


def process_gene_map(pop: str, overwrite: bool = False, batch_mode=False):
    """
    Process gene mapping file for a specific ancestry population.
    
    Args:
        pop: Ancestry population code (e.g., 'AFR', 'EUR')
        overwrite: Whether to overwrite existing files
    """
    from annotations import post_process_gene_map_ht

    if batch_mode:
        initialize_hail(batch_mode=batch_mode)
    
    # Input and output paths
    gene_map_ht_path = f"{DATA_PATH}/utils/gene_map/aou_{pop}_gene_map_{TRANCHE}.ht"
    processed_gene_map_ht_path = f"{DATA_PATH}/utils/gene_map/aou_{pop}_gene_map_processed_{TRANCHE}.ht"
    
    # Check if processed file exists and skip if not overwriting
    if not overwrite and hl.hadoop_exists(processed_gene_map_ht_path):
        print(f"Processed gene map file already exists for {pop} and overwrite=False. Skipping processing.")
        return
    
    # Check if input file exists
    if not hl.hadoop_exists(gene_map_ht_path):
        raise FileNotFoundError(f"Raw gene map file does not exist for {pop}. Run create_raw_gene_map first.")
    
    # Load and process gene map
    gene_map_ht = hl.read_table(gene_map_ht_path)
    gene_map_ht = post_process_gene_map_ht(gene_map_ht, freq_cutoff=0.01)
    print(f'---------Exporting processed gene mapping HT ({pop.upper()})-----------------')
    gene_map_ht = gene_map_ht.checkpoint(processed_gene_map_ht_path, overwrite=overwrite)
    gene_map_ht.describe()
    gene_map_ht.show()
    print(f'Completed processing for ancestry: {pop}')


def main(args):
    """
    Main function that creates gene mapping files based on command-line arguments.
    Can be run locally or submitted as a Hail Batch job.
    
    Args:
        args: Command-line arguments
    """
    
    if args.batch:
        if args.batch:
            backend = hb.ServiceBackend(
                billing_project="all-by-aou",
                remote_tmpdir=TMP_BUCKET,
            )
            b = hb.Batch(
                name=f"snpindel_merge_v8",
                requester_pays_project="aou-neale-gwas",
                default_python_image="amartinezcarrasco/hailgnomad:latest",
                backend=backend,
            )
        
        # Create jobs for each requested ancestry
        for pop in args.ancestries:
            print(f"Creating batch job for ancestry: {pop}")
            j = b.new_python_job(name=f"gene_map_{pop}")
            j.attributes['ancestry'] = pop
            
            if not args.skip_raw_gene_map_file:
                j.call(create_raw_gene_map, pop, args.overwrite, True)
            j.call(process_gene_map, pop, args.overwrite, True)

        b.run()
        
    else:
        initialize_hail(batch_mode=False, log_file="/create_annotations.log")
        for pop in args.ancestries:
            print(f"Processing ancestry {pop}")
            if not args.skip_raw_gene_map_file:
                create_raw_gene_map(pop, args.overwrite)
            process_gene_map(pop, args.overwrite)



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--batch",
        help="Run using Hail Batch instead of locally",
        action="store_true"
    )
    parser.add_argument(
        "--overwrite",
        help="Overwrite existing files",
        action="store_true"
    )
    parser.add_argument(
        "--ancestries",
        help="Comma-separated list of ancestries to process (default: all)",
        type=lambda s: s.split(',')
    )
    parser.add_argument(
        "--skip-raw-gene-map-file",
        help="Skip creating raw gene map file (only process existing files)",
        action="store_true"
    )

    args = parser.parse_args()
    main(args)

