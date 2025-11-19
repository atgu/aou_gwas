#!/usr/bin/env python3

__author__ = "Wenhan Lu"

import os
import hail as hl
import hailtop.batch as hb
import hailtop.fs as hfs
import argparse
from typing import List, Optional, Union


FastSparseGRM_root = 'gs://aou_analysis/v8/fastGRM'
fast_sparse_grm_image = "us-central1-docker.pkg.dev/aou-neale-gwas/saige/saige-king:0.0.6" # need to build on linux/amd64
hail_image="hailgenetics/hail:0.2.134-py3.11"
hard_filter_path = 'gs://aou_analysis/v8/data/utils/aou_v8_hard_filter.ht'
bed_prefix = f'{FastSparseGRM_root}/aou_v8_chrall'

def convert_mt_to_bed(mt_path: str, 
                     out_prefix: str = f'{bed_prefix}', 
                     min_maf: float = 0.05, 
                     max_missingness: float = 0.05, 
                     test: bool = True,
                     overwrite: bool = False) -> str:
    """
    Convert a Hail MatrixTable to PLINK BED format required by FastSparseGRM.
    
    Parameters:
    -----------
    mt_path : str
        Path to the Hail MatrixTable.
    out_prefix : str
        Prefix for output BED files.
    min_maf : float, optional
        Minimum minor allele frequency for variant filtering. Default is 0.05.
    max_missingness : float, optional
        Maximum missingness rate allowed for variants. Default is 0.05.
    overwrite : bool, optional
        If True, overwrite existing output files. Default is False.
    Returns:
    --------
    str
        Path to the output BED file.
    """

    hl.init(default_reference="GRCh38",
            master='local[32]',
            worker_memory="highmem",
            tmp_dir='gs://aou_tmp/v8/')
    # Load the MatrixTable
    if test:
        mt_path = 'gs://aou_wlu/250k_ld/aou_full.mt'
        mt = hl.read_matrix_table(mt_path)
        mt = mt.filter_rows(mt.locus.contig == 'chr21')
    else:
        mt = hl.read_matrix_table(mt_path)
        ht = hl.read_table(hard_filter_path)
        mt = mt.filter_cols(hl.is_defined(ht[mt.col_key]))

    mt.describe()
    print(mt.count())
    mt = mt.filter_rows(mt.locus.in_autosome())

    
    # Calculate variant QC metrics
    mt = hl.variant_qc(mt)
    
    # Filter variants based on MAF and missingness
    mt = mt.filter_rows((mt.variant_qc.AF[1] >= min_maf) & 
                       (mt.variant_qc.AF[1] <= (1 - min_maf)) &
                       (mt.variant_qc.call_rate >= (1 - max_missingness)))

    mt = mt.checkpoint(f'{out_prefix}.mt', overwrite=overwrite, _read_if_exists=not overwrite)
    print(mt.count())

    pruned_ht = hl.ld_prune(
        mt.GT,       # genotype calls
        r2=0.1,      # r² threshold
        bp_window_size=50000  # 50 kb window
    )
    pruned_mt = mt.filter_rows(hl.is_defined(pruned_ht[mt.row_key]))
    pruned_mt = pruned_mt.checkpoint(f'{out_prefix}_pruned.mt', overwrite=overwrite, _read_if_exists=not overwrite)
    print(pruned_mt.count())
    
    # Export to PLINK BED format
    hl.export_plink(pruned_mt, out_prefix)
    
    # Return paths to the BED, BIM, and FAM files
    return f"{out_prefix}.bed", f"{out_prefix}.bim", f"{out_prefix}.fam"


def main(args):
    backend = hb.ServiceBackend(
        billing_project=args.billing_project,
        remote_tmpdir=args.tmp_bucket
    )
    
    b = hb.Batch(
        name="aou_v8_fast_sparse_grm",
        backend=backend,
    )
    
    # Step 1: Create a job to convert the MatrixTable to BED format
    mt_to_bed_job = None
    if not hfs.exists(f'{bed_prefix}.bed') or args.overwrite:
        mt_to_bed_job = b.new_python_job(name=f"convert_mt_to_bed")
        mt_to_bed_job.image(hail_image)
        mt_to_bed_job.memory('highmem')
        mt_to_bed_job.cpu(32)
        mt_to_bed_job.env('PYSPARK_SUBMIT_ARGS', '--driver-memory 24g --executor-memory 24g pyspark-shell')
        
        mt_path = args.input_mt_path
        min_maf = args.min_maf
        max_missingness = args.max_missingness  
        test = args.test

        mt_to_bed_job.call(convert_mt_to_bed, mt_path, bed_prefix, min_maf, max_missingness, test, args.overwrite)

    plink_files = b.read_input_group(
        **{
            "bim": f'{bed_prefix}.bim',
            "fam": f'{bed_prefix}.fam',
            "bed": f'{bed_prefix}.bed',
        }
    )
    
    # Step 2: Run KING for IBD segments
    king_job = None
    if not hfs.exists(f'{bed_prefix}.king') or args.overwrite:
        king_job = b.new_job(name="run_king")
        king_job.image(fast_sparse_grm_image)
        if mt_to_bed_job is not None:
            king_job.depends_on(mt_to_bed_job)
    
        # Input files from previous job
        if king_job is not None:
            king_job.command(f"""
            echo "Removing 'chr' prefix via sed…"
            sed 's/^chr//' \
            {plink_files['bim']} \
            > {plink_files['bim']}.tmp && \
            mv {plink_files['bim']}.tmp \
            {plink_files['bim']}
            echo "Done."
            king -b {plink_files['bed'].replace('.bed', '')} --ibdseg --degree 4 --cpus {args.n_cpus} --prefix {king_job.ofile}
            """)
        b.write_output(king_job.ofile, f"{bed_prefix}.king")

    king_file = b.read_input(f"{bed_prefix}.king")
    
    # Step 3: Get ancestry divergence estimates
    divergence_job = None
    if not hfs.exists(f'{bed_prefix}.divergence') or args.overwrite:
        divergence_job = b.new_job(name="get_divergence")
        divergence_job.image(fast_sparse_grm_image)
        if king_job is not None:
            divergence_job.depends_on(king_job)
    
        divergence_job.command(f"""        
        cp /FastSparseGRM/inst/extdata/getDivergence_wrapper.R .
        cp /FastSparseGRM/inst/extdata/king .
        cp /FastSparseGRM/inst/extdata/plink .
        
        R CMD BATCH --vanilla '--args --prefix.in {plink_files['bed']} --file.seg {king_file} --num_threads {args.n_cpus} --degree {args.degree} --divThresh {args.div_threshold} --nRandomSNPs {args.n_random_snps} --prefix.out {divergence_job.ofile}' getDivergence_wrapper.R getDivergence.Rout
        """)
        b.write_output(divergence_job.ofile, f"{bed_prefix}.divergence")

    divergence_file = b.read_input(f"{bed_prefix}.divergence")
    
    # Step 4: Extract unrelated samples
    unrelated_job = None
    if not hfs.exists(f'{bed_prefix}.unrelated') or args.overwrite:
        unrelated_job = b.new_job(name="extract_unrelated")
        unrelated_job.image(fast_sparse_grm_image)
        if divergence_job is not None:
            unrelated_job.depends_on(divergence_job)
    
        unrelated_job.command(f"""
        # Copy necessary wrapper scripts
        cp /FastSparseGRM/inst/extdata/extractUnrelated_wrapper.R .
        
        # Extract unrelated samples
        R CMD BATCH --vanilla '--args --prefix.in {plink_files['bed']} --file.seg {king_file} --degree {args.degree} --file.div {divergence_file} --file.include {args.include_file} --prefix.out {unrelated_job.ofile}' extractUnrelated_wrapper.R extractUnrelated.Rout
        """)
        b.write_output(unrelated_job.ofile, f"{bed_prefix}.unrelated")

    unrelated_file = b.read_input(f"{bed_prefix}.unrelated")
    
    # Step 5: Run PCA
    pca_job = None
    if not hfs.exists(f'{bed_prefix}.pca') or args.overwrite:
        pca_job = b.new_job(name="run_pca")
        pca_job.image(fast_sparse_grm_image)
        if unrelated_job is not None:
            pca_job.depends_on(unrelated_job)
    
        pca_job.command(f"""
        # Copy necessary wrapper scripts
        cp /FastSparseGRM/inst/extdata/runPCA_wrapper.R .
        
        # Run PCA
        R CMD BATCH --vanilla '--args --prefix.in {plink_files['bed']} --file.unrels {unrelated_file} --prefix.out {pca_job.ofile} --no_pcs {args.n_pcs} --num_threads {args.n_cpus} --no_iter {args.n_iterations}' runPCA_wrapper.R runPCA.Rout
    """)
        b.write_output(pca_job.ofile, f"{bed_prefix}.pca")

    pca_file = b.read_input(f"{bed_prefix}.pca")
    
    # Step 6: Calculate Sparse GRM
    if not hfs.exists(f'{bed_prefix}.sparseGRM') or args.overwrite:
        sparse_grm_job = b.new_job(name="run_sparse_grm")
        sparse_grm_job.image(fast_sparse_grm_image)
        if pca_job is not None:
            sparse_grm_job.depends_on(pca_job)
    
        sparse_grm_job.command(f"""
        # Copy necessary wrapper scripts
        cp /FastSparseGRM/inst/extdata/calcSparseGRM_wrapper.R .

        # Calculate Sparse GRM
        R CMD BATCH --vanilla '--args --prefix.in {bed_prefix} --prefix.out {sparse_grm_job.ofile} --file.train {bed_prefix}.unrelated --file.score {bed_prefix}.pca --file.seg {bed_prefix}.king --num_threads {args.n_cpus} --no_pcs {args.n_pcs} --block.size {args.block_size} --max.related.block {args.max_related_block} --KINGformat.out {str(args.king_format_out).upper()} --degree {args.degree}' calcSparseGRM_wrapper.R calcSparseGRM.Rout
        """)
        b.write_output(sparse_grm_job.ofile, f"{bed_prefix}.sparseGRM")
    
    # Run all jobs
    b.run()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    # Core parameters
    parser.add_argument(
        "--input-mt-path",
        help="Path to the input Hail MatrixTable",
        default="gs://aou_wlu/v8/data/aou_v8.mt"
    )
    # FastSparseGRM parameters
    parser.add_argument(
        "--n-cpus",
        help="Number of CPU threads to use",
        type=int,
        default=24
    )
    parser.add_argument(
        "--degree",
        help="Degree of relationship to consider",
        type=int,
        default=4
    )
    parser.add_argument(
        "--min-maf",
        help="Minimum minor allele frequency",
        type=float,
        default=0.05
    )
    parser.add_argument(
        "--max-missingness",
        help="Maximum missingness rate",
        type=float,
        default=0.05
    )
    parser.add_argument(
        "--div-threshold",
        help="Divergence threshold",
        type=float,
        default=-0.02209709
    )
    parser.add_argument(
        "--n-random-snps",
        help="Number of random SNPs to use",
        type=int,
        default=0
    )
    parser.add_argument(
        "--include-file",
        help="File with sample IDs to include",
        default=""
    )
    parser.add_argument(
        "--n-pcs",
        help="Number of PCs to calculate",
        type=int,
        default=20
    )
    parser.add_argument(
        "--n-iterations",
        help="Number of iterations for PCA",
        type=int,
        default=10
    )
    parser.add_argument(
        "--block-size",
        help="Block size for sparse GRM calculation",
        type=int,
        default=5000
    )
    parser.add_argument(
        "--max-related-block",
        help="Maximum related block size",
        type=int,
        default=5000
    )
    parser.add_argument(
        "--king-format-out",
        help="Output in KING format",
        action="store_true"
    )
    parser.add_argument(
        "--test",
        help="Run with test data",
        action="store_true"
    )
    parser.add_argument(
        "--overwrite",
        help="Overwrite existing files",
        action="store_true"
    )
    # Batch parameters
    parser.add_argument(
        "--billing-project",
        help="Name of the billing project",
        default="all-by-aou"
    )
    parser.add_argument(
        "--tmp-bucket",
        help="Path to the temporary bucket",
        default='gs://aou_wlu/tmp'
    )
    
    args = parser.parse_args()
    main(args)


# Example usage:
# python3 fast_sparse_grm.py \
#     --test \
#     --n-cpus 16