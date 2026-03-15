#!/usr/bin/env Rscript

# Load required packages
packages = c('dplyr', 'optparse', 'tidyverse', 'Matrix', 'sparseMVN', 'spdep', 'qlcMatrix', 'ggplot2', 'stringr')
for(p in packages){
  if(!require(p, character.only = T)){
    install.packages(p, repos = "http://cran.us.r-project.org")
  }
}

# Define input arguments
option_list = list(
  make_option(c("-g", "--grm"), type="character", default=NULL,
              help="path to the sparse GRM (.mtx file)", metavar="character"),
  make_option(c("--heritabilities"), type="character", default='0.01,0.05,0.1,0.2,0.5',
              help="heritability of the phenotypes (a list of comma split values)", metavar="character"),
  make_option(c("-n", "--n"), type="integer", default=1000,
              help="number of raw phenotypes needed", metavar="integer"),
  make_option(c("-p", "--prevalence"), type="character", default='0.0001,0.001,0.01,0.02,0.05,0.1,0.2,0.5',
              help="number of cases (a list of comma split integers) or prevalences (a list of comma split value between 0 and 1)", metavar="character"),
  make_option(c("-r", "--n_reps"), type="integer", default= 20,
              help="number of draws to run", metavar="integer"),
  make_option(c("-s", "--samples"), type="character", default=NULL,
              help="path to the sample ID file", metavar="character"),
  make_option(c("-i", "--seed"), type="double", default=7993,
              help="set seed for random sampling", metavar="double"),
  make_option(c("-t", "--trait_type"), type="character", default='both',
              help="the type of phenotypes to simulate, (choose from `continuous`, `binary` or `both`)", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="random_phenos",
              help="output file name [default= %default]", metavar="character")
)
opt_parser = OptionParser(add_help_option=TRUE, option_list=option_list)
opt = parse_args(opt_parser);

# Confirm that a valid GRM is provided
if (is.null(opt$grm)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).\n", call.=FALSE)
}

# Define local temporary paths
local_grm_path <- tempfile(fileext = ".mtx")
local_samples_path <- tempfile(fileext = ".txt")

# Copy files from GCS using gsutil
system(paste0("gsutil -m cp ", opt$grm, " ", local_grm_path))
system(paste0("gsutil -m cp ", opt$samples, " ", local_samples_path))

# Read the local files
data = readMM(local_grm_path) # GRM matrix
samples = read.table(local_samples_path) # a list of sample ID

# Optional: Clean up temporary files if needed
file.remove(local_grm_path)
file.remove(local_samples_path)


n_samples = dim(data)[1] # Number of samples
print(paste0('N samples:', n_samples))

n_phenos_to_simulate = opt$n # Number of raw phenotypes to simulate
heritabilities = as.numeric(unlist(strsplit(opt$heritabilities, ','))) # List of heritabilities to use for the random phenotypes
n_heritabilities = length(heritabilities) # Number of heritabilities to use for the random phenotypes
binary_info = as.numeric(unlist(strsplit(opt$prevalence, ',')))
if(mean(binary_info) > 1){
  n_cases = binary_info # List of n_cases to use for the random phenotypes
  prevalences = n_cases/n_samples # Compute the prevalences from the n_cases
}else{
  prevalences = binary_info
  n_cases = prevalences*n_samples
}
n_prevalences = length(prevalences) # Number of prevalences to use for the random phenotypes
qnorms = qnorm(prevalences) # Compute the corresponding quantile of each prevalence in a standard normal distribution
n_draws = opt$n_reps # Number of draws to run per heritbility per n_case
n_phenos = ((n_prevalences + 1)*n_draws) # Number of phenotypes to draw per heritability

# Compute the Cholesky decomposition of the GRM
decomp = Cholesky(data)


# Build `n_phenos_to_simulate` raw random phenotypes from multivariate normal distributions
# with mean 0 and variance information from decomposition computed in the previous step
set.seed(opt$seed)
phenos = rmvn.sparse(n_phenos_to_simulate, rep(0, n_samples), decomp, prec=FALSE)

# Pivot the values in GRM to a data.frame
orig = data.frame(x=data@i + 1, y=data@j + 1, orig=data@x) %>% filter(x != y)

# Randomly sample a same number of uncorrelated pairs of samples (genetic correlation = 0)
# from a uniform distribution with min 0 and max = `n_samples`
# and merge it with the original data.frame
set.seed(opt$seed)
orig_uncorr = data.frame(x=round(runif(nrow(orig), max=n_samples)),
                         y=round(runif(nrow(orig), max=n_samples)),
                         orig=0)
orig = orig %>%
  union_all(orig_uncorr %>%
              anti_join(orig, by=c('x', 'y'))) %>%
  filter(x > 0 & y > 0)

# Annotate the corresponding phenotypic correlation from the previously simulated phenotypes
output = orig %>%
  rowwise() %>%
  mutate(pheno_corr=cor(phenos[,x], phenos[,y]))

# Plot and save a figure of the relationship between genetic and phenotypic correlation
figure <- output %>%
  ggplot + aes(x = orig, y = pheno_corr) +
  geom_point() + theme_classic() + geom_abline(intercept = 0, slope = 1, color = 'coral') + 
  labs(x = paste0('GRM (', toupper(str_split(opt$out, '_') %>% map_chr(., last)), ')'), y = 'Random Phenotype Correlation')
in2mm = 25.4
png(paste0(opt$out, '.png'), width = 100/in2mm, height = 70/in2mm, units = 'in', res = 300)
print(figure)
dev.off()
# ggsave(paste0(opt$out, '.png'))

# Create a data.frame with only userIds
output_df = samples %>%
  transmute(userId=V1)

set.seed(opt$seed)
new_h2 <- c()
new_slope <- c()
new_slope_se <- c()
new_corr <- c()
for (j in 1:n_heritabilities) {
  heritability = heritabilities[j]
  # Subset to the `n_phenos` unused raw phenotypes for this round
  start = (j - 1) * n_phenos + 1
  end = j * n_phenos
  use_phenos = t(phenos[start:end,])

  # Create `n_phenos` noise phenotypes sampled from the standard normal distribution
  noise_pheno = matrix(rnorm(length(use_phenos)), nrow=nrow(use_phenos))

  # Create `n_phenos` phenotypes with the current heritability from the selected raw and noise phenotypes
  use_phenos = use_phenos*sqrt(heritability) + noise_pheno*sqrt(1 - heritability)

  # Pivot the phenotypes to a long table and compute the phenotypic correlation and print the info
  test = orig %>%
    rowwise() %>%
    mutate(pheno_corr=cor(use_phenos[x,], use_phenos[y,]))
  print(summary(lm(pheno_corr ~ orig, test)))
  print(paste('Heritability:', heritability,
              '; Slope:', round(lm(pheno_corr ~ orig, test)$coefficients[[2]], 4),
              '; Correlation:', round(cor(test$pheno_corr, test$orig), 3)))
  new_h2 <- c(new_h2, heritability)
  new_slope <- c(new_slope, lm(pheno_corr ~ orig, test)$coefficients[[2]])
  new_slope_se <- c(new_slope_se, summary(lm(pheno_corr ~ orig, test))[[4]][2,2])
  new_corr <- c(new_corr, round(cor(test$pheno_corr, test$orig), 3))

  # Create `n_prevalence` binary phenotypes with the current heritability
  if(opt$trait_type %in% c('binary', 'both') | !is.null(n_cases)){
      for (i in 1:n_prevalences) {
        start = (i - 1) * n_draws + 1
        end = i * n_draws
        out = data.frame(use_phenos[,start:end] < qnorms[i])
        output_df = output_df %>%
          bind_cols(out %>%
                      rename_with(function(x) {if(n_draws==1){paste0('random_', heritability, '_', binary_info[i], '_1')}else{paste0('random_', heritability, '_', binary_info[i], '_', gsub('X', '', x))}}) %>%
                      mutate_all(as.numeric))
      }
  }

  # Create `n_heritabilities` continuous phenotypes
  if(opt$trait_type %in% c('continuous', 'both')){
    continuous_df = data.frame(use_phenos[,(n_prevalences*n_draws + 1):n_phenos])
    output_df = output_df %>%
      bind_cols(continuous_df %>%
                  rename_with(function(x) {if(n_draws==1){paste0('random_', heritability, '_continuous_1')}else{paste0('random_', heritability, '_continuous_', gsub('X', '', x))}}))
  }

}

plot_data <- data.frame(
  Heritability = new_h2,
  Slope = new_slope,
  Slope_SE = new_slope_se, 
  Correlation = new_corr
)

figure <- plot_data %>%
  ggplot + aes(x=Slope, y = Correlation, color = factor(Heritability)) +
  # geom_abline(intercept = 0, slope = 1, lty = 2, color = 'gray90') +
  geom_point(size = 2) + 
  geom_errorbarh(aes(xmax = Slope + Slope_SE, xmin = Slope-Slope_SE)) + 
  theme_classic() +
  scale_color_brewer(name = 'Heritability', palette = 'RdYlBu') + 
  labs(x = expression(Slope(r[pheno] * " ~ " * r[GRM])), y = paste0('Correlation  (', toupper(str_split(opt$out, '_') %>% map_chr(., last)), ')')) + 
  theme(legend.position = c(0.13, 0.7))
png(paste0(opt$out, '_heritability.png'), width = 100/in2mm, height = 70/in2mm, units = 'in', res = 300)
print(figure)
dev.off()


# Print a summary of the output random phenotype table
output_df %>% summary

# library(scales)
# binary_summary <- output_df %>%
#   pivot_longer(cols = starts_with('random'), names_to = 'label', values_to = 'value') %>%
#   mutate(heritability = as.numeric(str_split(label, '_') %>% map_chr(., 2)), 
#          prevalence = str_split(label, '_')%>% map_chr(., 3),
#          index = as.numeric(str_split(label, '_')%>% map_chr(., 4)))%>%
#   filter(prevalence != 'continuous') %>%
#   group_by(prevalence) %>%
#   dplyr::summarize(Proportion = mean(value), SE = sd(value) / sqrt(n()), n = n(), .groups = "drop")
# # Plot binary summary
# figure <- binary_summary %>%
#   mutate(prevalence = as.numeric(prevalence)) %>%
#   ggplot + aes(x = prevalence, y = Proportion, color = factor(prevalence)) +
#   geom_abline(intercept = 0, slope = 1, lty = 2, color = 'gray90') +
#   geom_point() + 
#   theme_classic() +
#   scale_x_continuous(label = percent) + 
#   scale_y_continuous(label = percent) + 
#   scale_color_brewer(name = NULL, palette = 'PuOr') + 
#   labs(x = paste0("True Prevalence  (", toupper(str_split(opt$out, '_') %>% map_chr(., last)), ')'),
#        y = "Simulated Prevalence") + 
#   theme(legend.position = c(0.2, 0.8))
# 
# png(paste0(opt$out, '_prevalence.png'), width = 100/in2mm, height = 70/in2mm, units = 'in', res = 300)
# print(figure)
# dev.off()

# Write a tsv file of the random phenotypes
# output_df %>%
#   write_tsv(paste0(opt$out, '.tsv'))
# print(paste0(opt$out, '.tsv'))
# write_delim(output_df, paste0(opt$out, '.tsv'), delim = "\t")

# Define local and GCS output paths
local_output_path <- tempfile(fileext = ".tsv")
# Define your GCS base path here - adjust if needed
gcs_base_path <- "gs://aou_analysis/v8/data/phenotype/raw/random_pheno/"
gcs_output_path <- paste0(gcs_base_path, opt$out, '.tsv')

print(paste0("Writing local file: ", local_output_path))
write_delim(output_df, local_output_path, delim = "\t")

# Copy the local file to GCS
print(paste0("Copying to GCS: ", gcs_output_path))
system(paste0("gsutil cp ", local_output_path, " ", gcs_output_path))

# Optional: Remove the local file after copying
print(paste0("Removing local file: ", local_output_path))
file.remove(local_output_path)

print("File successfully written to GCS.")


# Rscript random_pheno.R -g gs://aou_analysis/v8/data/utils/grm/aou_afr._relatednessCutoff_0.125_2000_randomMarkersUsed.sparseGRM.mtx -s gs://aou_analysis/v8/data/utils/grm/aou_afr._relatednessCutoff_0.125_2000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt  -o random_pheno_afr
# Rscript random_pheno.R -g gs://aou_analysis/v8/data/utils/grm/aou_amr._relatednessCutoff_0.125_2000_randomMarkersUsed.sparseGRM.mtx -s gs://aou_analysis/v8/data/utils/grm/aou_amr._relatednessCutoff_0.125_2000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt  -o random_pheno_amr
# Rscript random_pheno.R -g gs://aou_analysis/v8/data/utils/grm/aou_eas._relatednessCutoff_0.125_2000_randomMarkersUsed.sparseGRM.mtx -s gs://aou_analysis/v8/data/utils/grm/aou_eas._relatednessCutoff_0.125_2000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt  -o random_pheno_eas
# Rscript random_pheno.R -g gs://aou_analysis/v8/data/utils/grm/aou_eur._relatednessCutoff_0.125_2000_randomMarkersUsed.sparseGRM.mtx -s gs://aou_analysis/v8/data/utils/grm/aou_eur._relatednessCutoff_0.125_2000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt  -o random_pheno_eur
# Rscript random_pheno.R -g gs://aou_analysis/v8/data/utils/grm/aou_mid._relatednessCutoff_0.125_2000_randomMarkersUsed.sparseGRM.mtx -s gs://aou_analysis/v8/data/utils/grm/aou_mid._relatednessCutoff_0.125_2000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt  -o random_pheno_mid
# Rscript random_pheno.R -g gs://aou_analysis/v8/data/utils/grm/aou_sas._relatednessCutoff_0.125_2000_randomMarkersUsed.sparseGRM.mtx -s gs://aou_analysis/v8/data/utils/grm/aou_sas._relatednessCutoff_0.125_2000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt  -o random_pheno_sas
# Rscript random_pheno.R -g ALL_ACAF.sGRM.mtx -s ALL_grm_plink.samples  -o random_pheno_all
