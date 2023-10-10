# aou_gwas
Code for the All-x-All AoU project.

## Project Overview
### Description
GWAS and RVAS on the AoU data.

### Contributors
- Wenhan Lu ([@wlu04](https://github.com/wlu04))


## Analyses
### Notes on running jobs with Hail Query-on-Batch
Before running any QoB job, make sure the configurations are set up correctly by doing:
```
# hailctl config set batch/billing_project all-by-aou
# hailctl config set batch/remote_tmpdir gs://aou_tmp
# hailctl config set batch/tmp_dir gs://aou_tmp
# hailctl config set query/backend batch
hailctl config list
```


### Scripts
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
- Usage: `python3 --update-sample-util-ht --update-raw-phenotypes --update-phenotype-ht --annotate-phenotype-ht --update-meta-ht`
- Options: `--batch` run jobs with Hail Batch in case nothing else works
#### pre_process_random_pheno.py
- Description: pre-processing data to generate the GRM for producing population-specific random phenotypes
- Usage: `python3 pre_process_random_pheno.py --create-plink-file --create-sparse-grm --run-hail-ibd --ld-prune --pop all --overwrite`
- Options: 
  - `--pop` can be a comma-separated list of any combinations from `['afr', 'amr', 'eas', 'eur', 'mid', 'sas', 'all']`
  - `--ld-prune` whether to run ld pruning before exporting to plink files (required for GRM)
#### export_vds_to_bgen.py
- Description: chunking VDS into Bgen files each covers an interval of approximately `N_GENE_PER_GROUP` genes
- Usage: `python3 export_vds_to_bgen.py --test --mean_impute_missing --update-vds`