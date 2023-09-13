# aou_gwas
Code for the All-x-All AoU project.

## Project Overview
### Description
GWAS and RVAS on the AoU data.

### Contributors
- Wenhan Lu ([@wlu04](https://github.com/wlu04))


## Analyses
### Notes on running jobs with Hail Batch
Before running any Batch job, make sure the configurations are set up correctly by doing:
```
# hailctl config set batch/billing_project all-by-aou
# hailctl config set batch/remote_tmpdir gs://aou_tmp
# hailctl config set batch/tmp_dir gs://aou_tmp
hailctl config list
```
### Scripts
#### pre_process_random_phecode.py
- Description: pre-processing data to generate the GRM for producing population-specific random phenotypes
- Usage: `python3 pre_process_random_pheno.py --create-plink-file --create-sparse-grm --pop all --overwrite`
#### export_vds_to_bgen.py
- Description: chunking VDS into Bgen files each covers an interval of approximately `N_GENE_PER_GROUP` genes
- Usage: `python3 export_vds_to_bgen.py --test --mean_impute_missing --update-vds`