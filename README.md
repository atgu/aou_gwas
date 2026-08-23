# aou_gwas
Code for the *All-by-All AoU* project: GWAS and rare-variant association
(RVAS) analyses on the All of Us Research Program data.

> [!CAUTION]
> The scripts in this GitHub repository are **NOT** directly executable on the
> All of Us Researcher Workbench. We are in the process of adapting the code
> and developing a public workspace to enable users to reproduce the
> analyses within the Workbench environment. Updates will be posted here
> once the workspace becomes available.

## Project Overview
GWAS and RVAS on the AoU data, with companion meta-analyses against UKB and
FinnGen.

## Repository layout
- [R/](R/) — R utilities for plotting and downstream analysis (QQ / Manhattan,
  random phenotypes, PCA, etc.).
- [scripts/](scripts/) — top-level Python entry points for the v7 / pilot
  pipeline (SAIGE driver, phenotype processing, VAT reformatting, etc.).
- [utils/](utils/) — shared Python utilities, resources, Docker assets, and
  helpers for generating test inputs.
- [v8/](v8/) — production v8 pipeline organised by stage:
  [data_prep/](v8/data_prep/), [vep/](v8/vep/), [pipeline/](v8/pipeline/),
  [results/](v8/results/), [analysis/](v8/analysis/), [figures/](v8/figures/),
  [resources/](v8/resources/).
- [data/](data/) — small reference inputs and look-up tables checked into the
  repo.

## System requirements

### Software dependencies
- **Python** 3.11 with the packages pinned in
  [requirements.txt](requirements.txt) (notably `hail==0.2.135`,
  `gnomad==0.8.2`, `pysnptools`, `google-cloud-storage`, `scipy`, `numpy`).
- **R** ≥ 4.2 with `tidyverse`, `data.table`, `ggplot2`, `ggrepel`,
  `optparse`, `gridExtra`, `R.utils`, `readr`, `ggpubr`, `Matrix`,
  `MatrixModels`, `quantreg`, `scattermore`, `nloptr`, `lme4`, `pbkrtest`,
  `car`, `rstatix` (full list installed by [utils/Dockerfile](utils/Dockerfile)).
- **External tools (optional, stage-dependent):** SAIGE, PLINK 1.9, KING,
  bcftools, VEP 95 (GRCh38) — all wired up through the per-stage Dockerfiles
  under [v8/](v8/).
- **Cloud runtime:** [Hail Batch / Query-on-Batch](https://hail.is/) with
  Google Cloud Storage; some steps additionally use Hail Dataproc.

### Operating systems tested
- Linux (Ubuntu 22.04 inside the production Docker images).
- macOS (Apple Silicon / Intel) for local development and plotting.

### Non-standard hardware
None required for the published-code path. Several pipeline stages assume
access to the All of Us Researcher Workbench and a Google Cloud project with
billing enabled (Hail Batch, Dataproc, Cloud Storage). LD-pruning for the
EUR cohort is run on a Dataproc cluster because it exceeds Query-on-Batch
memory limits.

## Installation guide

### Local (plotting / lightweight scripts)
```
git clone https://github.com/atgu/aou_gwas.git
cd aou_gwas

# Python deps
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt

# R deps (interactive R session)
# install.packages(c("tidyverse","data.table","ggplot2","ggrepel","optparse",
#                    "gridExtra","R.utils","readr","ggpubr","Matrix",
#                    "MatrixModels","quantreg","scattermore","nloptr","lme4",
#                    "pbkrtest","car","rstatix"))
```
Typical install time on a modern laptop: **5–15 minutes**, dominated by R
package compilation.

### Docker (reproducible pipeline)
Each pipeline stage ships its own Dockerfile and is the recommended way to
reproduce the analyses:
- [utils/Dockerfile](utils/Dockerfile) — base Hail + R image
- [v8/pipeline/Dockerfile](v8/pipeline/Dockerfile) — SAIGE / KING / PLINK /
  FastSparseGRM
- [v8/pipeline/Dockerfile_bcftools](v8/pipeline/Dockerfile_bcftools),
  [v8/pipeline/Dockerfile_plink2](v8/pipeline/Dockerfile_plink2),
  [v8/pipeline/Dockerfile_gene_LD](v8/pipeline/Dockerfile_gene_LD)
- [v8/vep/Dockerfile](v8/vep/Dockerfile) — VEP 95 GRCh38
- [v8/results/Dockerfile](v8/results/Dockerfile),
  [v8/results/Dockerfile_gnomad](v8/results/Dockerfile_gnomad)

Build a stage image with, e.g.:
```
docker build -t aou_gwas-utils -f utils/Dockerfile .
```
Typical Docker image build time: **20–40 minutes** for the heaviest stages
(SAIGE / FastSparseGRM), 5–10 minutes for the plotting image.

## Demo

The end-to-end demo runs only inside the All of Us Researcher Workbench
because the underlying genomic data is access-controlled and cannot be
redistributed. A public Workbench that reproduces the demo end-to-end is in
preparation — see the note at the top of this README.

For users with Workbench access, individual scripts can be sanity-checked
without full input data:
- [utils/generate_random_vcf.py](utils/generate_random_vcf.py) creates a tiny
  simulated VCF (≈ seconds on a laptop) that exercises the variant-loading
  paths.
- [utils/generate_saige_test_files.py](utils/generate_saige_test_files.py)
  emits minimal SAIGE inputs.
- Most pipeline entry points accept `--test`, which restricts the workflow to
  a small interval and produces output in **a few minutes** instead of hours.

Expected outputs in `--test` mode mirror the full-scale outputs (Hail
Tables / MatrixTables, BGEN chunks, SAIGE summary statistics) but at
chromosome-fragment scale.

## Instructions for use

### Hail Query-on-Batch configuration
Before running any QoB job, configure `hailctl`:
```
hailctl config set batch/billing_project all-by-aou
hailctl config set batch/remote_tmpdir gs://aou_tmp
hailctl config set query/backend batch
hailctl config list
```

### R scripts
#### [R/random_phenos.R](R/random_phenos.R)
Generate random phenotypes from the SAIGE step-0 GRM.
```
Rscript R/random_phenos.R \
  -g ~/Downloads/250k_data_utils_grm_aou_afr._relatednessCutoff_0.125_2000_randomMarkersUsed.sparseGRM.mtx \
  -s ~/Downloads/250k_data_utils_grm_aou_afr._relatednessCutoff_0.125_2000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt \
  -p 0.01 \
  -o data/random_pheno_afr
```
See `Rscript R/random_phenos.R -h` for all flags.

#### [R/manhattan_and_qq_plot.R](R/manhattan_and_qq_plot.R)
Read locus and p-value information and generate QQ (and optional Manhattan)
plots. Required args: input file `-f`, p-value column `-p`, and either
chromosome `-c` + position `-bp` or a locus identifier column `-i`.
```
# single phenotype
Rscript R/manhattan_and_qq_plot.R -f ~/Downloads/amr_3446_both_sexes.txt.bgz -p Pvalue -c chr -b pos
# multiple phenotypes
Rscript R/manhattan_and_qq_plot.R -f ~/Downloads/afr_variant_results_pilot_af.txt.bgz -p Pvalue -m phenoname
```
Options: `-q TRUE` skips Manhattan plots; `-h` lists all options.

### Python scripts
#### [scripts/saige_aou.py](scripts/saige_aou.py)
```
python3 scripts/saige_aou.py --run_pipeline --phenos height,p_0.5_continuous_1 --pops eur --irnt --single_variant_only --skip_saige --skip_bgen --test
```

#### [scripts/reformat_vat.py](scripts/reformat_vat.py)
Load the original `.tsv.gz` VAT, parse non-string fields, sort, and key by
`locus`/`alleles`.
```
python3 scripts/reformat_vat.py
```

#### [scripts/process_phenotype.py](scripts/process_phenotype.py)
Load CSV phenotype files and sample-info TSVs into Hail Tables, merge, and
build the meta table used by association tests.
```
python3 scripts/process_phenotype.py --update-sample-util-ht --update-raw-phenotypes \
    --update-phenotype-ht --annotate-phenotype-ht --update-meta-ht
```
Add `--batch` to run via Hail Batch.

#### [scripts/pre_process_random_pheno.py](scripts/pre_process_random_pheno.py)
Pre-processing to build the GRM used for population-specific random
phenotypes.
```
# Query-on-Batch
python3 scripts/pre_process_random_pheno.py --create-plink-file --create-sparse-grm \
    --pop amr --overwrite-variant-ht --overwrite-variant-mt --ld-prune \
    --overwrite-ld-ht --overwrite-plink --overwrite-sample-file

# Dataproc
hailctl dataproc submit clustername scripts/pre_process_random_pheno.py --pop eur \
    --create-plink-file --ld-prune --overwrite-ld-ht --overwrite-plink --overwrite-sample-file \
    --pyfiles ~/Dropbox\ \(Partners\ HealthCare\)/github_repo/aou_gwas/
```
Notes:
1. LD pruning runs out of memory on QoB for EUR — use Dataproc.
2. On Dataproc, import via `from aou_gwas import *`.
3. On QoB, import via `from utils.utils import *` / `from utils.resources import *`.
4. Hail 0.2.124 caused transient errors; use the pinned version above.

`--pop` accepts any comma-separated subset of
`['afr','amr','eas','eur','mid','sas','all']`.

#### [scripts/export_vds_to_bgen.py](scripts/export_vds_to_bgen.py)
Chunk a VDS into BGEN files, each covering an interval of approximately
`N_GENE_PER_GROUP` genes.
```
python3 scripts/export_vds_to_bgen.py --test --mean_impute_missing --update-vds
```

### Reproducing manuscript results
The v8 production pipeline lives under [v8/](v8/). Each subdirectory
corresponds to a stage in the manuscript (data preparation, VEP annotation,
association testing, meta-analysis, results QC, figures); the per-stage
Dockerfile pins the exact environment used for the published runs. Detailed
pseudocode of the pipeline is provided in the Methods section of the
manuscript.

## License
This project is released under the MIT License — see [LICENSE](LICENSE).
