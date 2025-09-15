#!/usr/bin/env python3

import argparse
import subprocess
import pandas as pd
import numpy as np
from pathlib import Path
import logging
import sys

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    stream=sys.stdout
)

def run_saige_command(cmd):
    """Execute SAIGE command and check for errors"""
    try:
        subprocess.run(cmd, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        logging.error(f"Error running SAIGE command: {e}")
        raise

def find_top_variant(assoc_file, chrom, start_pos, end_pos, p_value_col=6):
    """Find the top variant in a region based on p-value"""
    try:
        df = pd.read_csv(assoc_file, sep='\t')
        region_df = df[
            (df['CHR'] == chrom) & 
            (df['POS'] >= start_pos) & 
            (df['POS'] <= end_pos)
        ]
        if region_df.empty:
            return None, None
        
        top_variant = region_df.sort_values(by=f'p.value', ascending=True).iloc[0]
        variant_id = f"{top_variant['CHR']}:{top_variant['POS']}_{top_variant['REF']}/{top_variant['ALT']}"
        return variant_id, top_variant[f'p.value']
    except Exception as e:
        logging.error(f"Error reading association file: {e}")
        return None, None

def get_input_cmd_args(args, input_type='common'):
    """Generate input-related command arguments based on file format"""
    if input_type == 'common':
        file_path = args.vcf_file if args.input_format == 'vcf' else args.bgen_file
        index_path = f"{file_path}.tbi" if args.input_format == 'vcf' else f"{file_path}.bgi"
    else:  # rare
        file_path = args.rare_vcf_file if args.input_format == 'vcf' else args.rare_bgen_file
        index_path = f"{file_path}.tbi" if args.input_format == 'vcf' else f"{file_path}.bgi"
    
    if args.input_format == 'vcf':
        return f"--vcfFile={file_path} --vcfFileIndex={index_path}"
    else:
        return f"--bgenFile={file_path} --bgenFileIndex={index_path}"

def conditional_analysis_workflow(args):
    # Initial parameters
    p_value_threshold = args.p_value_threshold
    independent_variants = []
    current_p_value = 0
    iteration = 0
    
    # Step 1: Run initial association test
    logging.info("Running initial association test...")
    input_args = get_input_cmd_args(args, 'common')
    
    initial_cmd = f"""
    step2_SPAtests.R \\
        {input_args} \\
        --chrom={args.chrom} \\
        --start={args.start_pos} \\
        --end={args.end_pos} \\
        --GMMATmodelFile={args.model_file} \\
        --varianceRatioFile={args.var_ratio_file} \\
        --sampleFile={args.sample_file} \\
        --SAIGEOutputFile={args.output_prefix}_initial.txt \\
        --minInfo={args.min_info} \\
        --minMAC={args.min_mac} \\
        --minMAF={args.min_maf} \\
        --LOCO=FALSE
    """
    run_saige_command(initial_cmd)
    
    current_assoc_file = f"{args.output_prefix}_initial.txt"
    
    # Iterative conditional analysis
    while True:
        # Find top variant in the region
        top_variant, current_p_value = find_top_variant(
            current_assoc_file, 
            args.chrom, 
            args.start_pos, 
            args.end_pos
        )
        
        if not top_variant or current_p_value >= p_value_threshold:
            logging.info("No more significant variants found.")
            break
            
        logging.info(f"Found significant variant: {top_variant} (p={current_p_value})")
        independent_variants.append(top_variant)
        
        # Prepare conditioning string
        condition_string = ",".join(independent_variants)
        
        # Run conditional analysis
        cond_output = f"{args.output_prefix}_cond_{iteration}.txt"
        cond_cmd = f"""
        step2_SPAtests.R \\
            {input_args} \\
            --chrom={args.chrom} \\
            --start={args.start_pos} \\
            --end={args.end_pos} \\
            --GMMATmodelFile={args.model_file} \\
            --varianceRatioFile={args.var_ratio_file} \\
            --sampleFile={args.sample_file} \\
            --SAIGEOutputFile={cond_output} \\
            --minInfo={args.min_info} \\
            --minMAC={args.min_mac} \\
            --minMAF={args.min_maf} \\
            --condition={condition_string} \\
            --LOCO=FALSE
        """
        run_saige_command(cond_cmd)
        
        current_assoc_file = cond_output
        iteration += 1
        
        if iteration >= args.max_iterations:
            logging.warning(f"Reached maximum iterations ({args.max_iterations})")
            break
    
    # Save independent variants
    if independent_variants:
        with open(f"{args.output_prefix}.independent_variants.txt", 'w') as f:
            f.write("\n".join(independent_variants))
        logging.info(f"Found {len(independent_variants)} independent variants")
    
    # Run rare variant analysis if requested
    if args.run_rare_variant:
        logging.info("Running rare variant analysis...")
        rare_input_args = get_input_cmd_args(args, 'rare')
        
        # Unconditional rare variant test
        rare_cmd = f"""
        step2_SPAtests.R \\
            {rare_input_args} \\
            --chrom={args.chrom} \\
            --start={args.start_pos} \\
            --end={args.end_pos} \\
            --GMMATmodelFile={args.model_file} \\
            --varianceRatioFile={args.var_ratio_file} \\
            --sampleFile={args.sample_file} \\
            --SAIGEOutputFile={args.output_prefix}_rare_uncond.txt \\
            --minMAF=0 \\
            --maxMAF={args.rare_maf_threshold} \\
            --LOCO=FALSE \\
            --groupFile={args.group_file}
        """
        run_saige_command(rare_cmd)
        
        # Conditional rare variant test
        if independent_variants:
            rare_cond_cmd = f"""
            step2_SPAtests.R \\
                {rare_input_args} \\
                --chrom={args.chrom} \\
                --start={args.start_pos} \\
                --end={args.end_pos} \\
                --GMMATmodelFile={args.model_file} \\
                --varianceRatioFile={args.var_ratio_file} \\
                --sampleFile={args.sample_file} \\
                --SAIGEOutputFile={args.output_prefix}_rare_cond.txt \\
                --minMAF=0 \\
                --maxMAF={args.rare_maf_threshold} \\
                --LOCO=FALSE \\
                --groupFile={args.group_file} \\
                --condition={condition_string}
            """
            run_saige_command(rare_cond_cmd)

def main():
    parser = argparse.ArgumentParser(description='SAIGE Conditional Analysis Pipeline')
    
    # Input format
    parser.add_argument('--input-format', choices=['vcf', 'bgen'], default='vcf',
                      help='Input format: VCF or BGEN (default: vcf)')
    
    # Required arguments - Common to both formats
    parser.add_argument('--chrom', required=True, help='Chromosome')
    parser.add_argument('--start-pos', type=int, required=True, help='Start position')
    parser.add_argument('--end-pos', type=int, required=True, help='End position')
    parser.add_argument('--model-file', required=True, help='GMMAT model file')
    parser.add_argument('--var-ratio-file', required=True, help='Variance ratio file')
    parser.add_argument('--sample-file', required=True, help='Sample file')
    parser.add_argument('--output-prefix', required=True, help='Output prefix for results')
    
    # Format-specific arguments
    vcf_group = parser.add_argument_group('VCF input options')
    vcf_group.add_argument('--vcf-file', help='VCF file for common variants')
    vcf_group.add_argument('--rare-vcf-file', help='VCF file for rare variants')
    
    bgen_group = parser.add_argument_group('BGEN input options')
    bgen_group.add_argument('--bgen-file', help='BGEN file for common variants')
    bgen_group.add_argument('--rare-bgen-file', help='BGEN file for rare variants')
    
    # Optional arguments
    parser.add_argument('--min-info', type=float, default=0.3, help='Minimum info score')
    parser.add_argument('--min-mac', type=int, default=10, help='Minimum minor allele count')
    parser.add_argument('--min-maf', type=float, default=0.01, help='Minimum MAF for common variants')
    parser.add_argument('--p-value-threshold', type=float, default=5e-8, help='P-value significance threshold')
    parser.add_argument('--max-iterations', type=int, default=10, help='Maximum number of conditional iterations')
    
    # Rare variant analysis arguments
    parser.add_argument('--run-rare-variant', action='store_true', help='Run rare variant analysis')
    parser.add_argument('--rare-maf-threshold', type=float, default=0.01, help='MAF threshold for rare variants')
    parser.add_argument('--group-file', help='Group file for rare variant analysis')
    
    args = parser.parse_args()
    
    # Validate input format and files
    if args.input_format == 'vcf':
        if not args.vcf_file:
            parser.error("--vcf-file is required when using VCF format")
        if args.run_rare_variant and not args.rare_vcf_file:
            parser.error("--rare-vcf-file is required for rare variant analysis with VCF format")
    else:  # bgen
        if not args.bgen_file:
            parser.error("--bgen-file is required when using BGEN format")
        if args.run_rare_variant and not args.rare_bgen_file:
            parser.error("--rare-bgen-file is required for rare variant analysis with BGEN format")
    
    # Validate rare variant arguments
    if args.run_rare_variant and not args.group_file:
        parser.error("--group-file is required when --run-rare-variant is specified")
    
    logging.info("Starting SAIGE conditional analysis pipeline...")
    conditional_analysis_workflow(args)
    logging.info("Pipeline completed successfully!")

if __name__ == "__main__":
    main()