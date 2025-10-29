# aou_gwas
Code for the All-x-All AoU project.


<div style="background-color: #ffebeb; border: 1px solid #ff0000; padding: 10px; border-radius: 5px;">
  <p style="color: #ff0000; margin: 0;">Note: The scripts in this GitHub repository are not compatible with the All of Us Researcher Workbench. We are in the process of adapting the code and developing a public workspace to enable users to reproduce the analyses within the Workbench environment. Updates will be provided once the workspace becomes available.!</p>
</div>

## Project Overview
### Description
GWAS and RVAS on the AoU data.

### Contributors
- Wenhan Lu ([@wlu04](https://github.com/wlu04))

## Resources


## Analyses
### Notes on running jobs with Hail Query-on-Batch
Before running any QoB job, make sure the configurations are set up correctly by doing:
```
hailctl config set batch/billing_project all-by-aou
hailctl config set batch/remote_tmpdir gs://aou_tmp
# hailctl config set batch/tmp_dir gs://aou_tmp
hailctl config set query/backend batch
hailctl config list
```

### R Scripts
#### random_phenos.R
- Description:
  - Generate random phenotypes using GRM from SAIGE step 0
- Usage: 
  ```
   Rscript random_phenos.R \
   -g ~/Downloads/250k_data_utils_grm_aou_afr._relatednessCutoff_0.125_2000_randomMarkersUsed.sparseGRM.mtx \
   -s ~/Downloads/250k_data_utils_grm_aou_afr._relatednessCutoff_0.125_2000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt \ 
   -p 0.01 \
   -o data/random_pheno_afr
   ```
  - for more information: Run `Rscript random_phenos.R -h`
#### manhattan_and_qq_plot.R
- Description:
  - Read locus and pvalue information from input `-f`
  - Generate QQ-plot(s)
  - (Optional) Gerate Manhattan plot 
  - Required arguments: input file `-f`, with pvalue column `-p`, chromosome column `-c` + base-pare column `-bp` OR comma-seperated locus identifier column `-i`
- Usage: 
  - single phenotype: `Rscript manhattan_and_qq_plot.R -f ~/Downloads/amr_3446_both_sexes.txt.bgz -p Pvalue -c chr -b pos`
  - multiple phenotypes: `Rscript manhattan_and_qq_plot.R -f ~/Downloads/afr_variant_results_pilot_af.txt.bgz -p Pvalue -m phenoname`
- Options: 
  - `-q TRUE`: generate qq-plot(s) only, ignore manhattan plot(s)
  - for more information: Run `Rscript manhattan_and_qq_plot.R -h`

### Python Scripts
#### saige_aou.py
- Usage: `python3 saige_aou.py --run_pipeline --phenos height,p_0.5_continuous_1 --pops eur --irnt --single_variant_only --skip_saige --skip_bgen --test`
#### reformat_vat.py
- Description: 
  - Load the original .tsv.gz VAT from the original bucket
  - Parse the non-string fields to the correct datatype
  - Sort the table and key by `locus` and  `alleles`
- Usage: `python3 reformat_vat.py`

#### process_phenotype.py
- Description: 
  - Load original .csv phenotype files to HTs and parse the non-string fields to the correct datatype (`--update-raw-phenotypes`) 
  - Load the .tsv sample information files to HTs and parse the non-string fields to the correct datatype (`--update-sample-util-ht`)
  - Merge the parsed phenotype HTs (`--update-phenotype-ht`)  and annotate sample information (`--annotate-phenotype-ht`) 
  - Generate sample meta information table with standard covariates to be used for association tests (`--update-meta-ht`)
- Usage: `python3 process_phenotype.py --update-sample-util-ht --update-raw-phenotypes --update-phenotype-ht --annotate-phenotype-ht --update-meta-ht`
- Options: `--batch` run jobs with Hail Batch in case nothing else works
#### pre_process_random_pheno.py
- Description: pre-processing data to generate the GRM for producing population-specific random phenotypes
- Usage: 
  - QoB: `python3 pre_process_random_pheno.py --create-plink-file --create-sparse-grm --pop amr --overwrite-variant-ht --overwrite-variant-mt --ld-prune --overwrite-ld-ht --overwrite-plink --overwrite-sample-file`
  - Dataproc: `hailctl dataproc submit clustername pre_process_random_pheno.py --pop eur --create-plink-file --ld-prune --overwrite-ld-ht --overwrite-plink --overwrite-sample-file --pyfiles ~/Dropbox\ \(Partners\ HealthCare\)/github_repo/aou_gwas/`
  - **Note**: 
    1. LD pruning will run out of memory on QoB for EUR, which should be run on dataproc.
    2. To run the pipeline on dataproc, use `from aou_gwas import * `
    3. To run the pipeline on QoB, use `from utils.utils import * // from utils.resources import * `
    4. hail 0.2.124 causes transient error 
- Options: 
  - `--pop` can be a comma-separated list of any combinations from `['afr', 'amr', 'eas', 'eur', 'mid', 'sas', 'all']`
  - `--ld-prune` whether to run ld pruning before exporting to plink files (required for GRM)
#### export_vds_to_bgen.py
- Description: chunking VDS into Bgen files each covers an interval of approximately `N_GENE_PER_GROUP` genes
- Usage: `python3 export_vds_to_bgen.py --test --mean_impute_missing --update-vds`