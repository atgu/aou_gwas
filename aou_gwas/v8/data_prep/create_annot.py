#!/usr/bin/env python3

import sys
import os
import argparse
import hail as hl
import hailtop.batch as hb
import hailtop.fs as hfs

# Global constants
TRANCHE = "v8"
ANALYSIS_BUCKET = "gs://aou_analysis"
MY_BUCKET = 'gs://aou_amc'
TMP_BUCKET = 'gs://aou_tmp'
DATA_PATH = f"{ANALYSIS_BUCKET}/{TRANCHE}/data"
SNPINDEL_OUT_PATH = f'{DATA_PATH}/vep/aou_{TRANCHE}_vep_full.ht'


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


def merge_snpindel(overwrite = False, batch_mode=False):
    """Merge SNP and INDEL VEP annotations"""

    if batch_mode:
        initialize_hail(batch_mode=batch_mode, log_file=f"/merge_snpindel_ht.log")

    from gnomad.utils.vep import process_consequences
       
    indel_vep_path = f'{ANALYSIS_BUCKET}/{TRANCHE}/data/vep/aou_vds_variant_data_row_{TRANCHE}_vep.ht'
    snp_vep_path = 'gs://gcp-public-data--gnomad/resources/context/grch38_context_vep_annotated.v105.ht'
    indel_vep_ht = hl.read_table(indel_vep_path)
    snp_vep_ht = hl.read_table(snp_vep_path)
    indel_fields_to_drop = ['ancestral', 'context']
    shared_fields_to_drop = ['uniparc', 'trembl', 'swissprot']
    
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
    
    print("Merging SNP and INDEL VEP annotations...")
    vep_ht = snp_ht.union(indel_ht, unify=True)
    
    process_vep_ht = process_consequences(vep_ht)
    vep_ht = vep_ht.annotate(
        worst_csq_by_gene_canonical=process_vep_ht[vep_ht.key].vep.worst_csq_by_gene_canonical
    )
    
    vep_ht = vep_ht.naive_coalesce(3000)
    print(f"Writing merged VEP table to {SNPINDEL_OUT_PATH}...")
    vep_ht = vep_ht.checkpoint(SNPINDEL_OUT_PATH, overwrite=overwrite)
    
    return vep_ht


def create_raw_gene_map(pop: str, overwrite: bool = False, batch_mode=False):
    """
    Create raw gene mapping file for a specific ancestry population.
    
    Args:
        pop: Ancestry population code (e.g., 'AFR', 'EUR')
        overwrite: Whether to overwrite existing files
        batch_mode: Whether the function is being run in a batch job
    """
    from aou_gwas.utils.annotations import create_gene_map_ht

    if batch_mode:
        initialize_hail(batch_mode=batch_mode, log_file=f"/create_raw_gene_map_{pop}.log")
    
    # Load VEP table
    snp_indel_vep_path = SNPINDEL_OUT_PATH
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
    
    call_stats_ht_path = f"{ANALYSIS_BUCKET}/{TRANCHE}/data/utils/call_stats/exome_pruned/{pop}_exome_call_stats.ht"
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
    
    gene_map_ht = create_gene_map_ht(snp_indel_vep_ht, freq_field='freq')
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
    from aou_gwas.utils.annotations import post_process_gene_map_ht

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
    1. Merge SNP and INDEL VEP annotations
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
                merge_job = b.new_python_job(name=f"Merge SNP Indels")
                merge_job.memory('highmem')
                merge_job.cpu(8)
                merge_job.call(merge_snpindel, args.overwrite, True)
                
                if args.gene_map and args.ancestries:
                    for pop in args.ancestries:
                        job = b.new_python_job(name=f"gene_map_{pop}")
                        job.depends_on(merge_job)
                        job.attributes['ancestry'] = pop
                        
                        if not args.skip_raw_gene_map_file:
                            job.call(create_raw_gene_map, pop, args.overwrite, True)
                        job.call(process_gene_map, pop, args.overwrite, True)
            
            elif args.gene_map and args.ancestries:
                for pop in args.ancestries:
                    job = b.new_python_job(name=f"gene_map_{pop}")
                    job.attributes['ancestry'] = pop
                    
                    if not args.skip_raw_gene_map_file:
                        job.call(create_raw_gene_map, pop, args.overwrite, True)
                    job.call(process_gene_map, pop, args.overwrite, True)
            
            b.run()
            
        else:
            # Run on dataproc cluster or locally with QoB
            initialize_hail(log_file=f"/gene_map_generation_{TRANCHE}.log")
            if args.merge:
                vep_merged_ht = merge_snpindel(args.overwrite)
            
            if args.gene_map and args.ancestries:
                for pop in args.ancestries:
                    print(f"Processing ancestry {pop}")
                    if not args.skip_raw_gene_map_file:
                        create_raw_gene_map(pop, args.overwrite)
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
        help="Run SNP-INDEL merge",
        action="store_true",
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