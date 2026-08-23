#!/usr/bin/env python3

__author__ = "Wenhan Lu"

import hail as hl
import hailtop.batch as hb
import argparse
import pandas as pd
import pickle
import hailtop.fs as hfs
import logging
from typing import List
from tqdm import tqdm
import sys

TRANCHE = "v8"
MY_BUCKET = 'gs://aou_wlu'
ANALYSIS_BUCKET = "gs://aou_analysis/v8"
RESULT_BUCKET = "gs://aou_results/414k"
DATA_PATH = f'{ANALYSIS_BUCKET}/data'
TMP_BUCKET = 'gs://aou_tmp/v8'
BINARY_CATEGORIES = ['r_drug', 'pfhh_survey', 'mcc2_phecode', 'mcc2_phecodex', 'mhwb_survey']
QUANTITATIVE_CATEGORIES = ['physical_measurement', 'lab_measurement']


def read_pickle_dict(dict_path: str):
    dict = pd.read_pickle(dict_path)
    return dict


def write_pickle_dict(output: str, lst: list):
    with hfs.open(output, "wb") as f:
        pickle.dump(lst, f)
    f.close()


def copy_files(original_path, target_path):
    """Copy files from requester-pays bucket to non-requester-pays bucket"""
    hl.init(
        gcs_requester_pays_configuration='aou-neale-gwas'
    )
    hl.hadoop_copy(original_path, target_path)


def add_liftover_rg37_to_rg38_ht(ht: hl.Table):
    """
    Add liftover to Hail Table from rg37 to rg38
    :param Table ht: Hail Table to add liftover on
    :return: Hail Table
    :rtype: hl.Table
    """
    rg37 = hl.get_reference("GRCh37")
    rg38 = hl.get_reference("GRCh38")
    if not rg37.has_liftover("GRCh38"):
        rg37.add_liftover(
            "gs://hail-common/references/grch37_to_grch38.over.chain.gz", rg38
        )
    ht = ht.annotate(new_locus=hl.liftover(ht.locus, "GRCh38"))
    ht = ht.filter(hl.is_defined(ht.new_locus))
    ht = ht.key_by(locus=ht.new_locus, alleles=ht.alleles)
    ht = ht.drop('new_locus')
    return ht


def produce_assoc_file_plink_ld_clump(gwas_path: str, pheno: str, ancestry: str, output_path: str):
    """
    Produce association file for PLINK LD clumping from Hail Table
    """
    # hl.init(
    #     master='local[32]',
    #     tmp_dir=TMP_BUCKET,
    #     default_reference="GRCh38",
    # )
    if gwas_path.endswith('.ht'):
        acaf = hl.read_table(gwas_path)
        acaf = acaf.annotate(
            BETA = hl.float64(acaf.BETA),
            SE = hl.float64(acaf.SE),
        )
        exome = hl.read_table(gwas_path.replace('genome', 'exome'))
        exome = exome.annotate(
            BETA = hl.float64(exome.BETA),
            SE = hl.float64(exome.SE),
        )
        gwas = acaf.union(exome).distinct()
        # Filter for defined p-values
        gwas = gwas.filter(hl.is_defined(gwas.Pvalue))
        gwas = gwas.annotate(SNP=hl.variant_str(gwas.locus, gwas.alleles), P=gwas.Pvalue)
    if ancestry.lower() != 'meta':
        gwas = gwas.filter((gwas.AF_Allele2 > 0.001) & (gwas.AF_Allele2 * gwas.n_cases >= 5))
    else:
        gwas = gwas.filter((gwas.AF_Allele2 > 0.001) & (gwas.AF_Allele2 * gwas.N_cases >= 5))

    for chrom in tqdm(range(1, 25)):
        chrom_str = 'X' if chrom == 23 else ('Y' if chrom == 24 else str(chrom))
        gwas_sub = gwas.filter(gwas.locus.contig == f'chr{chrom_str}')
        n_sig = gwas_sub.aggregate(hl.agg.count_where(gwas_sub.P < 5e-8))
        if n_sig >0: 
            gwas_sub = gwas_sub.checkpoint(f'{TMP_BUCKET}/clumping/{ancestry}_{pheno}_assoc_chr{chrom_str}.ht', overwrite=True)
            gwas_sub.describe()
            gwas_sub = gwas_sub.key_by(gwas_sub.SNP)
            gwas_sub = gwas_sub.select(gwas_sub.P)
            output_file = f'{output_path}/{ancestry}_{pheno}_chr{chrom_str}.assoc'
            gwas_sub.export(output_file, header=True)
            print(f"Exported association file to: {output_file}")
        else: 
            print('No significant association')
            # sys.exit(233)


def main(args):
    logging.basicConfig(
        format="%(levelname)s (%(name)s %(lineno)s): %(message)s",
        level="INFO",
        filename="gwas_clumping.log",
    )
    logger = logging.getLogger("AoU_GWAS_Clumping")
    logger.setLevel(logging.INFO)

    backend = hb.ServiceBackend(
        billing_project=args.billing_project,
        remote_tmpdir=args.tmp_bucket
    )
    b = hb.Batch(
        name='AoU GWAS Clumping V8',
        backend=backend,
    )

    # Read phenotype list
    phenos_to_run_by_ancestry = read_pickle_dict(
        f'{DATA_PATH}/phenotype/summary/pheno_by_ancestry_dict_both.dict'
    )

    # Loop through ancestries
    ancestries = args.ancestries.split(',') if args.ancestries else ['AFR', 'AMR', 'EAS', 'EUR', 'MID', 'SAS', 'META']
    
    for anc in ancestries:
        anc = anc.upper()
        print(f'--- Processing ancestry: {anc}---')
        
        # Use 'ALL' for META ancestry, otherwise use the ancestry code
        plink_ancestry = 'ALL' if anc == 'META' else anc
        output_prefix = f'gs://aou_wlu/v8_analysis/clumped_results/{anc}'
        
        copy_plink_jobs = {}
        filter_plink_jobs = {}
        if not args.skip_data_prep:
            # Step 1: Create sample keep file for this ancestry
            fid_iid_job = None
            if not hl.hadoop_exists(f"gs://aou_wlu/v8_analysis/ACAF_plink/{plink_ancestry}.samples"):
                fid_iid_job = b.new_job(name=f"Make keep file for {anc}")
                fid_iid_job.declare_resource_group(ofile={"keep": "samples_to_keep.txt"})
                
                if anc == 'META':
                    # For META, concatenate sample IDs from all six ancestry groups
                    ancestry_groups = ['AFR', 'AMR', 'EAS', 'EUR', 'MID', 'SAS']
                    sample_files = [b.read_input(f"gs://aou_analysis/v8/data/utils/grm/{anc_group}_grm_plink.samples") 
                                    for anc_group in ancestry_groups]
                    
                    # Concatenate all sample files
                    cat_command = f"cat {' '.join([str(f) for f in sample_files])}"
                    fid_iid_job.command(f"""
                        {cat_command} | awk '{{print 0, $1}}' > {fid_iid_job.ofile['keep']}
                    """)
                else:
                    # For other ancestries, use the ancestry-specific sample file
                    iid_list_gcs = f"gs://aou_analysis/v8/data/utils/grm/{plink_ancestry}_grm_plink.samples"
                    iid_list = b.read_input(iid_list_gcs)
                    fid_iid_job.command(f"""
                        awk '{{print 0, $1}}' {iid_list} > {fid_iid_job.ofile['keep']}
                    """)
                
                b.write_output(fid_iid_job.ofile['keep'], f"gs://aou_wlu/v8_analysis/ACAF_plink/{plink_ancestry}.samples")
            
            # Step 2: Copy and filter chromosome-specific PLINK files from requester-pays bucket
            for chrom in range(1, 25):
                chrom = 'X' if chrom == 23 else chrom
                chrom = 'Y' if chrom == 24 else chrom
                if args.test:
                    chrom = 22
                if hl.hadoop_exists(f'gs://aou_wlu/v8_analysis/combined_plink/{plink_ancestry}_merged.chr{chrom}.bed'):
                    continue
                # Copy chromosome PLINK files
                copy_plink_job = None
                if not hl.hadoop_exists(f"gs://aou_wlu/v8_analysis/ACAF_plink/ACAF.chr{chrom}.bed"):
                    copy_plink_job = b.new_python_job(name=f"Copy ACAF PLINK chr{chrom}")
                    copy_plink_job.call(copy_files, 
                        f"gs://fc-aou-datasets-controlled/v8/wgs/short_read/snpindel/acaf_threshold/plink_bed/chr{chrom}.bed", 
                        f"gs://aou_wlu/v8_analysis/ACAF_plink/ACAF.chr{chrom}.bed")
                    copy_plink_job.call(copy_files, 
                        f"gs://fc-aou-datasets-controlled/v8/wgs/short_read/snpindel/acaf_threshold/plink_bed/chr{chrom}.bim", 
                        f"gs://aou_wlu/v8_analysis/ACAF_plink/ACAF.chr{chrom}.bim")
                    copy_plink_job.call(copy_files, 
                        f"gs://fc-aou-datasets-controlled/v8/wgs/short_read/snpindel/acaf_threshold/plink_bed/chr{chrom}.fam", 
                        f"gs://aou_wlu/v8_analysis/ACAF_plink/ACAF.chr{chrom}.fam")
                    copy_plink_jobs[chrom] = copy_plink_job
                
                # Filter chromosome PLINK files to this ancestry
                filter_plink_job = None
                if not hl.hadoop_exists(f"gs://aou_wlu/v8_analysis/ACAF_plink/{plink_ancestry}_filtered.chr{chrom}.bed"):
                    plink_prefix = f"gs://aou_wlu/v8_analysis/ACAF_plink/ACAF.chr{chrom}"
                    
                    filter_plink_job = b.new_job(name=f"Filter chr{chrom} PLINK to {anc}")
                    if fid_iid_job is not None:
                        filter_plink_job.depends_on(fid_iid_job)
                    if chrom in copy_plink_jobs:
                        filter_plink_job.depends_on(copy_plink_jobs[chrom])
                    fid_iid_file = b.read_input(f"gs://aou_wlu/v8_analysis/ACAF_plink/{plink_ancestry}.samples")
                    storage = '2000G' if chrom == 1 or chrom == 2 else '1000G'
                    if anc=='EUR' and chrom in [3, 4, 5, 6, 7, 8]:
                        storage = '2000G'
                    if anc=='META' and chrom in [3, 4, 5, 6, 7, 8, 9, 10, 11, 12]:
                        storage = '2000G'
                    filter_plink_job.storage(storage)
                    filter_plink_job._machine_type = "n1-highmem-16"
                    filter_plink_job.declare_resource_group(ofile={"bed": "filtered.bed", "bim": "filtered.bim", "fam": "filtered.fam"})
                    filter_plink_job.image('us-central1-docker.pkg.dev/aou-neale-gwas/saige/saige_google_plink:0.4')
                    filter_plink_job.command(f"gcloud storage cp {plink_prefix}.bed /io/original.bed")
                    filter_plink_job.command(f"gcloud storage cp {plink_prefix}.bim /io/original.bim")
                    filter_plink_job.command(f"gcloud storage cp {plink_prefix}.fam /io/original.fam")
                    filter_plink_job.command(f"""
                        cd /io
                        plink2 --bfile original \\
                            --keep {fid_iid_file} \\
                            --make-bed \\
                            --out /io/filtered
                    """)
                    filter_plink_job.command('ls')
                    filter_plink_job.command(f'mv /io/filtered.bed {filter_plink_job.ofile["bed"]}')
                    filter_plink_job.command(f'mv /io/filtered.bim {filter_plink_job.ofile["bim"]}')
                    filter_plink_job.command(f'mv /io/filtered.fam {filter_plink_job.ofile["fam"]}')
                    b.write_output(filter_plink_job.ofile['bed'], f"gs://aou_wlu/v8_analysis/ACAF_plink/{plink_ancestry}_filtered.chr{chrom}.bed")
                    b.write_output(filter_plink_job.ofile['bim'], f"gs://aou_wlu/v8_analysis/ACAF_plink/{plink_ancestry}_filtered.chr{chrom}.bim")
                    b.write_output(filter_plink_job.ofile['fam'], f"gs://aou_wlu/v8_analysis/ACAF_plink/{plink_ancestry}_filtered.chr{chrom}.fam")
                    filter_plink_jobs[chrom] = filter_plink_job
                
                if args.test:
                    break
        
        # Step 3: Merge ACAF and exome PLINK files per chromosome for this ancestry
        merge_acaf_exome_jobs = {}
        if not args.skip_merge:
            for chrom in range(1, 25):
                chrom_str = 'X' if chrom == 23 else ('Y' if chrom == 24 else str(chrom))
                if args.test:
                    chrom_str = '22'
                
                merge_job = None
                output_path = f"gs://aou_wlu/v8_analysis/combined_plink/{plink_ancestry}_merged.chr{chrom_str}.bed"
                if not hl.hadoop_exists(output_path):
                    acaf_exists = hl.hadoop_exists(f"gs://aou_wlu/v8_analysis/ACAF_plink/{plink_ancestry}_filtered.chr{chrom_str}.bed")
                    exome_exists = hl.hadoop_exists(f"gs://aou_wlu/v8_within_gene_LD/{plink_ancestry}/{plink_ancestry}_filtered.chrom{chrom_str}.bed")
                    
                    if acaf_exists and exome_exists:
                        merge_job = b.new_job(name=f"Merge ACAF+Exome chr{chrom_str} for {anc}")
                        storage = '1000G'
                        if anc == 'META':
                            if chrom in range(1,13):
                                storage = '4000G'
                            elif chrom in range(13,20):
                                storage = '2000G'
                            elif chrom == 23:
                                storage = '2000G'
                        if anc == 'EUR':
                            if chrom in range(1,3):
                                storage = '3000G'
                            elif chrom in range(3,13):
                                storage = '2000G'
                        merge_job.storage(storage)
                        merge_job._machine_type = "n1-highmem-16"
                        merge_job.image('us-central1-docker.pkg.dev/aou-neale-gwas/saige/saige_google_plink:0.5')
                        merge_job.declare_resource_group(ofile={"bed": "merged.bed", "bim": "merged.bim", "fam": "merged.fam"})
                        
                        # Depend on filter jobs if they exist
                        if chrom_str in filter_plink_jobs and filter_plink_jobs[chrom_str] is not None:
                            merge_job.depends_on(filter_plink_jobs[chrom_str])

                        if chrom == 1 and anc =='META':
                            merge_job.spot(True)
                        
                        # Copy ACAF PLINK files
                        merge_job.command(f'gsutil cp gs://aou_wlu/v8_analysis/ACAF_plink/{plink_ancestry}_filtered.chr{chrom_str}.bed acaf.bed')
                        merge_job.command(f'gsutil cp gs://aou_wlu/v8_analysis/ACAF_plink/{plink_ancestry}_filtered.chr{chrom_str}.bim acaf.bim')
                        merge_job.command(f'gsutil cp gs://aou_wlu/v8_analysis/ACAF_plink/{plink_ancestry}_filtered.chr{chrom_str}.fam acaf.fam')
                        
                        # Copy exome PLINK files
                        merge_job.command(f'gsutil cp gs://aou_wlu/v8_within_gene_LD/{plink_ancestry}/{plink_ancestry}_filtered.chrom{chrom_str}.bed exome.bed')
                        merge_job.command(f'gsutil cp gs://aou_wlu/v8_within_gene_LD/{plink_ancestry}/{plink_ancestry}_filtered.chrom{chrom_str}.bim exome.bim')
                        merge_job.command(f'gsutil cp gs://aou_wlu/v8_within_gene_LD/{plink_ancestry}/{plink_ancestry}_filtered.chrom{chrom_str}.fam exome.fam')
                        
                        # Merge with PLINK 1.9 (--bmerge is stable in PLINK 1.9)
                        merge_job.command(f"""
                            plink --bfile acaf \\
                                --bmerge exome \\
                                --make-bed \\
                                --out /io/merged_temp
                        """)
                        
                        # Remove duplicate variants (keep first occurrence) using PLINK2
                        merge_job.command(f"""
                            plink2 --bfile /io/merged_temp \\
                                --rm-dup force-first \\
                                --make-bed \\
                                --out /io/merged
                        """)
                        
                        merge_job.command(f'mv /io/merged.bed {merge_job.ofile["bed"]}')
                        merge_job.command(f'mv /io/merged.bim {merge_job.ofile["bim"]}')
                        merge_job.command(f'mv /io/merged.fam {merge_job.ofile["fam"]}')
                        
                        b.write_output(merge_job.ofile['bed'], f"gs://aou_wlu/v8_analysis/combined_plink/{plink_ancestry}_merged.chr{chrom_str}.bed")
                        b.write_output(merge_job.ofile['bim'], f"gs://aou_wlu/v8_analysis/combined_plink/{plink_ancestry}_merged.chr{chrom_str}.bim")
                        b.write_output(merge_job.ofile['fam'], f"gs://aou_wlu/v8_analysis/combined_plink/{plink_ancestry}_merged.chr{chrom_str}.fam")
                        
                        merge_acaf_exome_jobs[chrom_str] = merge_job
                
                if args.test:
                    break
        
        # Step 4: Process phenotypes for this ancestry
        phenos_to_run = phenos_to_run_by_ancestry[anc.lower()]
        print(f'Number of phenotypes: {len(phenos_to_run)}')
        print(f'First 5 phenotypes: {phenos_to_run[0:5]}')
        pheno_prep_depend_on = None
        if not args.skip_pheno_prep:
            if args.pheno_idx is not None:
                phenos_to_run = phenos_to_run[args.pheno_idx*300:min(args.pheno_idx*300+300, len(phenos_to_run))]
            for pheno in tqdm(phenos_to_run):
                if args.test:
                    pheno = 'height'
                if '_to_' in pheno:
                    continue
                # Construct paths
                input_ht_path = f'{RESULT_BUCKET}/ht_results/{anc.upper()}/phenotype_{pheno}/genome_variant_results.ht'
                
                output_assoc_path = f'{output_prefix}/assoc_files/{anc}_{pheno}.assoc'
                
                # Check if input HT exists
                if not hl.hadoop_exists(f'{input_ht_path}/_SUCCESS'):
                    print(f'Input HT does not exist for {anc} - {pheno}, skipping...')
                    continue
                
                if not hl.hadoop_exists(output_assoc_path) or args.overwrite:
                    if not args.dataproc:
                        print(f'Creating association file for {anc} - {pheno}...')
                        # Create association file job
                        j_assoc = b.new_python_job(name=f"assoc_{anc}_{pheno}")
                        j_assoc.image("hailgenetics/hail:0.2.133-py3.11")
                        j_assoc.memory('standard')
                        j_assoc.storage('1000Gi')
                        j_assoc.env('PYSPARK_SUBMIT_ARGS', '--driver-memory 32g --executor-memory 32g pyspark-shell')
                        j_assoc.call(
                            produce_assoc_file_plink_ld_clump,
                            input_ht_path,
                            pheno,
                            anc,
                            f'{output_prefix}/assoc_files'
                        )
                        pheno_prep_depend_on = j_assoc
                    else:
                        produce_assoc_file_plink_ld_clump(
                            input_ht_path,
                            pheno,
                            anc,
                            f'{output_prefix}/assoc_files'
                        )
                if args.test:
                    break

        if not args.skip_clumping:
            for pheno in tqdm(phenos_to_run):
                if args.test:
                    pheno = 'CV_401'
                if '_to_' in pheno:
                    continue
                print(f'Clumping for {anc} - {pheno}...')
                # input_assoc_file = b.read_input(f'{output_prefix}/assoc_files/{anc}_{pheno}.assoc')
                # Create clumping jobs by chromosome
                for chrom in range(1, 25):
                    chrom_str = 'X' if chrom == 23 else ('Y' if chrom == 24 else str(chrom))
                    if not hl.hadoop_exists(f'{output_prefix}/assoc_files/{anc}_{pheno}_chr{chrom_str}.assoc'):
                        print(f'{pheno}-chr{chrom_str} has no high-quality signals')
                        continue
                    input_assoc_file = b.read_input(f'{output_prefix}/assoc_files/{anc}_{pheno}_chr{chrom_str}.assoc')

                    output_clumped_path = f'{output_prefix}/results/phenotype_{pheno}/{anc}_{pheno}_chr{chrom_str}.clumped'
                    if hl.hadoop_exists(output_clumped_path) and not args.overwrite:
                        print(f'Clumped results already exist for {anc} - {pheno} - chr{chrom_str}, skipping...')
                        continue 
                    
                    j_clump = b.new_job(name=f'clump_{anc}_{pheno}_chr{chrom_str}')

                    # Use regular machines by default (cheapest), only highmem for largest files
                    # j_clump.memory('39GiB')
                    if anc in ['EUR']:
                        if chrom in [1, 2, 3, 4]:
                            j_clump.storage('800Gi')
                        else:
                            j_clump.storage('500Gi')
                    if anc in ['META']:
                        if chrom in [1, 2]:
                            j_clump.storage('1200Gi')
                            j_clump._machine_type = 'n1-highmem-8'  # Only large META files need highmem
                        elif chrom in [3, 4, 5, 6]:
                            j_clump.storage('1000Gi')
                            j_clump._machine_type = 'n1-highmem-8'  # Only large META files need highmem
                        elif chrom in [7, 8, 9, 10, 11, 12, 16, 17]:
                            j_clump.storage('800Gi')
                        else:
                            j_clump.storage('500Gi')
                    if anc in ['AFR', 'AMR']:
                        j_clump.storage('300Gi')
                    if anc in ['EAS', 'MID', 'SAS']:
                        j_clump.storage('50Gi')
                    # j_clump.image('us-docker.pkg.dev/ukbb-diversepops-neale/pan-ukbb-docker-repo/nbaya_plink:0.1')
                    j_clump.image('us-central1-docker.pkg.dev/aou-neale-gwas/saige/saige_google_plink:0.5')
                    if pheno_prep_depend_on is not None:
                        j_clump.depends_on(pheno_prep_depend_on)
                    
                    # Depend on merge job if it exists, otherwise on filter job
                    if chrom_str in merge_acaf_exome_jobs and merge_acaf_exome_jobs[chrom_str] is not None:
                        j_clump.depends_on(merge_acaf_exome_jobs[chrom_str])
                    elif chrom_str in filter_plink_jobs and filter_plink_jobs[chrom_str] is not None:
                        j_clump.depends_on(filter_plink_jobs[chrom_str])

                    j_clump.declare_resource_group(ofile={'clumped': '{root}.clumped'})
                    work_dir = f'/io'
                    
                    j_clump.command(f'mkdir -p {work_dir}')
                    # Upload settings
                    j_clump.command('gcloud config set storage/parallel_composite_upload_enabled True')
                    j_clump.command('gcloud config set storage/parallel_composite_upload_threshold 150MiB')
                    j_clump.command('gcloud config set storage/parallel_composite_upload_component_size 64MiB')
                    # Aggressive download settings for large files (500GB-1TB)
                    j_clump.command('gcloud config set storage/sliced_object_download_threshold 150MiB')
                    j_clump.command('gcloud config set storage/sliced_object_download_max_components 128')
                    j_clump.command('gcloud config set storage/sliced_object_download_component_size 64MiB')
                    # Aggressive parallelism settings
                    j_clump.command('gcloud config set storage/process_count 16')
                    j_clump.command('gcloud config set storage/thread_count 4')

                    j_clump.command(f'gcloud storage cp '
                                    f'gs://aou_wlu/v8_analysis/combined_plink/{plink_ancestry}_merged.chr{chrom_str}.bed '
                                    f'gs://aou_wlu/v8_analysis/combined_plink/{plink_ancestry}_merged.chr{chrom_str}.bim '
                                    f'gs://aou_wlu/v8_analysis/combined_plink/{plink_ancestry}_merged.chr{chrom_str}.fam '
                                    f'{work_dir}/')

                    # Run PLINK clumping with the filtered files
                    j_clump.command(
                        f'''
                        cd {work_dir}
                        plink --bfile {plink_ancestry}_merged.chr{chrom_str} \\
                        --clump {input_assoc_file} \\
                        --clump-p1 {args.clump_p1} \\
                        --clump-p2 {args.clump_p2} \\
                        --clump-r2 {args.clump_r2} \\
                        --clump-kb {args.clump_kb} \\
                        --out {j_clump.ofile}
                        ''')

                    b.write_output(j_clump.ofile, f'{output_prefix}/results/phenotype_{pheno}/{anc}_{pheno}_chr{chrom_str}') 

                    if args.test:
                        break

                if args.test:
                    break
       
            if args.test:
                break
    
        if args.test:
            break
    
    if not args.dataproc:
        b.run()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--billing-project",
        help="Name of the billing project",
        nargs="?",
        default='all-by-aou'
    )
    parser.add_argument(
        "--tmp-bucket",
        help="Path to the temporary bucket",
        nargs="?",
        default='gs://aou_tmp/v8'
    )
    parser.add_argument(
        "--batch-name",
        help="Name of the batch",
        nargs="?",
        default='aou_gwas_clumping_v8'
    )
    parser.add_argument(
        "--docker-image",
        nargs="?",
        help="Docker image to use",
        default='hailgenetics/hail:0.2.136-py3.11',
    )
    parser.add_argument(
        "--bfile-path",
        help="Path to ancestry-specific sample lists (for --keep files)",
        nargs="?",
        default='gs://aou_analysis/v8/data/utils/grm'
    )

    parser.add_argument(
        "--output-path",
        help="Path to the output directory",
        nargs="?",
        default='gs://aou_wlu/v8_analysis/clumped_results'
    )
    parser.add_argument(
        "--test",
        help="Whether to run test on a subset of the data",
        action="store_true"
    )
    parser.add_argument(
        "--overwrite",
        help="Whether to overwrite existing files",
        action="store_true"
    )
    parser.add_argument(
        "--ancestries",
        help="Comma-separated list of ancestries to process (e.g., 'EUR,AFR')",
        nargs="?",
        default=None
    )
    parser.add_argument(
        "--clump-p1",
        help="P-value threshold for index SNPs",
        type=float,
        default=5e-8
    )
    parser.add_argument(
        "--clump-p2",
        help="P-value threshold for clumped SNPs",
        type=float,
        default=5e-8
    )
    parser.add_argument(
        "--clump-r2",
        help="LD r^2 threshold for clumping",
        type=float,
        default=0.1
    )
    parser.add_argument(
        "--clump-kb",
        help="Physical distance threshold for clumping (kb)",
        type=int,
        default=500
    )
    parser.add_argument(
        "--skip-data-prep",
        help="Whether to skip data preparation steps",
        action="store_true"
    )
    parser.add_argument(
        "--skip-clumping",
        help="Whether to skip the clumping step",
        action="store_true"
    )
    parser.add_argument(
        "--skip-pheno-prep",
        help="Whether to skip phenotype preparation steps",
        action="store_true"
    )
    parser.add_argument(
        "--skip-merge",
        help="Whether to skip merging ACAF and exome PLINK files",
        action="store_true"
    )
    parser.add_argument(
        "--dataproc",
        help="Whether to run on Dataproc instead of Hail Batch",
        action="store_true"
    )
    parser.add_argument(
        '--pheno-idx',
        help='Index of the phenotype to process (for Dataproc mode)',
        type=int,
        default=None
    )
    args = parser.parse_args()

    main(args)

# Example usage:
# python gwas_clumping.py --test --ancestries EUR
# python gwas_clumping.py --ancestries EUR,AFR --clump-p1 5e-8 --overwrite
