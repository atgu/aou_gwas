source('~/CHARR/R/constants.R')
source('~/Partners HealthCare Dropbox/Wenhan Lu/github_repo/ukbb_exomes_pleiotropy/R/constants.R')
# library(gg.layers)
# install.packages("ggrastr", dependencies = TRUE)
# install.packages("R.utils", dependencies = TRUE)
# install.packages("tidygraph", dependencies = TRUE)
# install.packages("pROC", dependencies = TRUE)
library(emojifont)
library(gganimate)
library(fuzzyjoin)
ANCESTRIES = c('afr', 'amr', 'eas', 'eur', 'mid', 'sas', 'meta')
pops = c('AFR-like', 'AMR-like', 'EAS-like', 'EUR-like', 'MID-like', 'SAS-like', 'META')
names(pops) <- c('afr', 'amr', 'eas', 'eur', 'mid', 'sas', 'meta')
color_mde =  '#EEA9B8'

axa_themes <- theme(
  axis.title = element_text(face = 'plain', size = 15, color = 'black'),
  axis.text =  element_text(size = 14),
  legend.text = element_text(size = 13),
  legend.title = element_text(size = 13, face = 'plain')
)

pop_colors = c('afr' = color_afr,
               'amr' = color_amr,
               'eas' = color_eas,
               'eur' = color_nfe,
               'sas' = color_sas,
               'mid' = color_mde,
               'meta' = color_oth)
version = 'v8'
P_VALUE_FIELDS = c("skato"="Pvalue", "skat"="Pvalue_SKAT", "burden"="Pvalue_Burden", "genome_variant"='Pvalue', "exome_variant"='Pvalue')

## PARAMETERS
category_colors <- c('#d56f3e', '#880d1e', '#43aa8b', '#4f345a', '#4f145a', '#334195')
names(category_colors) <- c('physical_measurement', 'r_drug', 'pfhh_survey', 'phecode', 'phecodeX', 'lab_measurement')
category_labels <- c('Physical measurement', 'Drug', 'PFHH', 'Phecode', 'PhecodeX', 'Lab measurement')
category_names <- c('physical_measurement', 'r_drug', 'pfhh_survey', 'phecode', 'phecodeX', 'lab_measurement')
names(category_labels) <- category_names
categories <- c('Physical \n measurement', 'Drug', 'PFHH', 'Phecode', 'PhecodeX', 'Lab measurement')
names(categories) <- c('physical_measurement', 'r_drug', 'pfhh_survey', 'phecode', 'phecodeX', 'lab_measurement')
file_names <- c('gene'= 'exome_gene_burden_0.001', 'exome'='exome_gene', 'genome'='ACAF_variant')

label_type = labeller(ancestry = pops, annotation=annotation_names,trait_type=categories)

category_colors <- c('#334195', '#d56f3e', '#7e57c2', '#4f345a','#880d1e', '#43aa8b')
category_labels <- c('Lab measurements', 'Physical measurements', 'Phecode', 'PhecodeX', 'Prescriptions', 'Self-reported')
categories <- category_labels
names(category_colors) <- category_labels
category_color_scale<- scale_color_manual(name = NULL, values = category_colors, breaks = names(category_colors), labels = category_labels) 
category_fill_scale <- scale_fill_manual(name = NULL, values = category_colors, breaks = names(category_colors), labels = category_labels)


## DIRECTORIES
if(version == 'v8'){
  aou_qc_data_path <- '~/Dropbox (Partners HealthCare)/aou_v8/data/'
  aou_qc_fig_path <- '~/Dropbox (Partners HealthCare)/aou_v8/figures/'
  aou_util_data_path <- '~/Partners HealthCare Dropbox/Wenhan Lu/github_repo/aou_gwas/data/'
}else{
  aou_qc_data_path <- '~/Partners HealthCare Dropbox/Wenhan Lu/aou/aou_250k_lambda_gc/data/'
  aou_qc_fig_path <- '~/Partners HealthCare Dropbox/Wenhan Lu/aou/aou_250k_lambda_gc/'
  aou_util_data_path <- '~/Partners HealthCare Dropbox/Wenhan Lu/github_repo/aou_gwas/data/'
}


## UTIL DATA
pheno_info_full <- read_tsv('~/Dropbox (Partners HealthCare)/aou_v8/data/phenotype/aou_v8_phenotype_info_qced_with_meta.tsv') 

matched_pheno <- read_csv(paste0(aou_util_data_path, 'aou_ukb_matched_phenotype.csv')) %>% filter(keep)
if(version == 'v8'){
  full_pheno_sum <- read_delim('~/Dropbox (Partners HealthCare)/aou_v8/data/phenotype/aou_v8_phenotype_info_qced_with_meta_both.tsv')
}else{
  full_pheno_sum <- read_delim(paste0(aou_util_data_path, '250k_qc_aou_phenotype_meta_info_250k.txt.bgz'), delim = '\t') %>%
    dplyr::select(-lambda_gc_exome, -lambda_gc_acaf, -lambda_gc_gene_burden_001)
}

PHENO_COUNT <- data.frame(
  category = rep(c('physical_measurement', 'lab_measurement', 'self-reported', 'r_drug', 'phecode', 'phecodeX'), 2),
  name = rep(c('Physical measurements', 'Lab measurements', 'Self-reported', 'Prescriptions', 'Phecode', 'PhecodeX'), 2),
  n_pheno = c(10, 60, 135, 1054, 1857, 3524, 10, 86, 145, 1067, 1843, 3426),
  version = rep(c('v7', 'v8'), each=6)
) 

aou_data <- data.frame(
  ancestry = rep(c('afr', 'amr', 'eas', 'eur', 'mid', 'sas', 'meta'), 2),
  version = rep(c('v7', 'v8'), each=7),
  n_samples = c(50457, 40437, 5344, 114343, 691, 2944, 214216, 77444, 71540, 9488, 227273, 1153, 5132, 392030),
  n_samples_raw = c(56913, 45035, 5706, 133581, 942, 3217, 245394, 84148, 79106, 10099, 234353, 1545, 5579, 414830),
  n_snps = c(332493469, 291730812, 104529343, 557334532, 33345022, 69753882, 0,
             393945194, 383787195, 142092571, 723283093, 45614446, 92589058, 0),
  n_indels = c(51208798, 42660159, 18199781, 71601047, 8497672, 13830435, 0, 
               62115876, 55431773, 23279499, 95472977, 10477170, 17002700, 0),
  n_singleton_snps = c(130753560, 127367894, 56742383, 274772030, 13992364, 36249581, 0,
                       147714769, 163441987, 76383332, 334961537, 20312396, 48998827, 0), 
  n_singleton_indels = c(17185945, 15272638, 7034206, 30944153, 2548155, 4955494, 0,
                         21212079, 20781828, 9523458, 41183260, 3357107, 6544974, 0),
  # n_lab = c(55,53,45,58,23,45, 60),
  # p_lab = c(55,53,45,58,23,45, 60)/60, 
  # n_phecode = c(573, 477, 45, 1000, 1, 20, 1857),
  # p_phecode = c(573, 477, 45, 1000, 1, 20, 1857)/1857,
  # n_phecodeX = c(869, 776, 80, 1361, 0, 42, 3524),
  # p_phecodeX = c(869, 776, 80, 1361, 0, 42, 3524)/3524,
  # n_drug = c(758, 715, 357, 857, 42, 288, 1054),
  # p_drug = c(758, 715, 357, 857, 42, 288, 1054)/1054,
  # n_pfhh = c(20, 17, 0, 78, 0, 0, 135),
  # p_pfhh = c(20, 17, 0, 78, 0, 0, 135)/135,
  # n_physical = rep(10, 7),
  # p_physical = rep(1, 7),
  n_pheno = c(2272, 2037, 534, 3348, 95, 402, 3417, 2882, 2750, 1138, 3589, 327, 881, 3602)
) %>%
  mutate(n_variants = n_snps + n_indels,
         n_singletons = n_singleton_snps + n_singleton_indels)

# FUNCTIONS
sig_cnt_summary = function(data, sig_col = 'sig_cnt'){
  summary = data %>%
    dplyr::summarize(mean = mean(get(sig_col), na.rm = TRUE),
                     prop = sum(get(sig_col) > 0, na.rm = T) / n(),
                     sem = sd(get(sig_col), na.rm = T)/sqrt(n()),
                     sig_cnt = sum(get(sig_col) > 0, na.rm = T),
                     cnt = n()) %>%
    mutate(sd = 1.96* sqrt(prop* (1-prop)/cnt))
  return(summary)
}

get_freq_interval <- function(freq){
  interval = case_when(
    freq <= 1e-5  ~ '[0, 0.00001]',
    freq <= 1e-4  ~ '(0.00001, 0.0001]',
    freq <= 1e-3  ~ '(0.0001, 0.001]',
    freq <= 5e-3  ~ '(0.001, 0.005]',
    freq <= 1e-2  ~ '(0.005, 0.01]',
    freq <= 1e-1  ~ '(0.01, 0.1]',
    freq > 1e-1  ~ '(0.1, )',
  )
  return(interval)
}

scientific_10 <- function(x) {
  # Create scientific notation, remove "e", replace with "10^", and clean up the "1 x"
  formatted <- gsub("e\\+?", " %*% 10^", scales::scientific_format()(x))  # Replace "e" with "10^" and remove "+"
  formatted <- gsub("^1 \\\\%\\*% ", "", formatted)  # Remove the leading "1 %*% " for numbers like 1 x 10^
  parse(text = formatted)
}

output_figure <- function(figure, folder, name, height, width){
  png(paste0('~/Dropbox (Partners HealthCare)/aou_v8/figures/', folder, '/', name, '.png'), width=width, height=height, units = 'in', res = 300)
  print(figure)
  dev.off()
}

add_color_scale <- function(figure, name){
  if(name == 'category'){
    figure <- figure +
      scale_color_manual(name = 'Category', values = category_colors, breaks = names(category_colors), labels = category_labels) +
      scale_fill_manual(name = 'Category', values = category_colors, breaks = names(category_colors), labels = category_labels)
  }else{
    figure <- figure +
      scale_color_jama(name = 'Trait type', breaks = c('continuous', 'binary'), labels = c('Continuous', 'Binary')) +
      scale_fill_jama(name = 'Trait type', breaks = c('continuous', 'binary'), labels = c('Continuous', 'Binary'))
  }
  return(figure)
}
plot_cdf <- function(data, field, group=NULL, save=F, output_path = '~/Dropbox (Partners HealthCare)/aou/phenotype/figures/phenotype_cdf/'){
  if(is.null(group)){
    cdf_data <- data %>% 
      dplyr::arrange(get(field)) %>%
      dplyr::mutate(id = row_number(),
                    n_total = dplyr::n()) %>%
      dplyr::mutate(cdf = id/n_total)
    
    figure <- cdf_data %>%
      ggplot + 
      aes(x=get(field), y=cdf)
  }else{
    cdf_data <- data %>% 
      dplyr::arrange(get(field)) %>%
      dplyr::group_by(get(group)) %>%
      dplyr::mutate(id = row_number(),
                    n_total = dplyr::n()) %>%
      dplyr::mutate(cdf = id/n_total)
    figure <- cdf_data %>%
      ggplot + 
      aes(x=get(field), y=cdf, color = get(group))
    figure <- add_color_scale(figure, group)
  }
  figure <- figure +
    labs(x = str_to_title(str_replace_all(field, '_', ' ')), y = 'Proportion') +
    # scale_x_log10(label=comma, 
    #               # breaks = c(10, 100, 200, 1000, 10000, 100000)
    #               ) + 
    # geom_vline(xintercept = 239734*0.5, lty=2) +
    geom_line() + themes 
  if(save){
    png(paste0(output_path, 'cdf_', field, if_else(is.null(group), '', paste0('_by_', group)), ".png"), width=5, height=3, units = 'in', res = 300)
    print(figure)
    dev.off()
  }
  return(figure)
  
}

filter_7digit_atcs <- function(df, prevalence_cut = 0.02, overlap_cut = 0.8) {
  # Only operate on 7-digit codes
  df <- df %>%
    mutate(
      parent_code = if_else(nchar(phenoname) == 7, str_sub(phenoname, 1, 5), NA)
    )
  
  # Merge with parent prevalence
  df_parents <- df %>%
    filter(nchar(phenoname) == 5) %>%
    select(parent_code = phenoname, parent_cases = n_cases_v8)
  
  df_joined <- df %>%
    left_join(df_parents, by = "parent_code")
  
  df_filtered <- df_joined %>%
    mutate(child_parent_overlap =  (n_cases_v8 / parent_cases)) %>%
    filter(
      # keep if not 7-digit
      nchar(phenoname) != 7 |
        (
          prevalence_v8 >= prevalence_cut &                # at least 5% prevalence
            (is.na(parent_cases) |        # OR no parent prevalence to compare
               child_parent_overlap < overlap_cut)  # not ≥80% of parent
        )
    ) 
  
  return(df_filtered)
}

load_sample_filtered_pheno_summary <- function(){
  v8_pheno_info <- rbind(
    read_delim('~/Dropbox (Partners HealthCare)/aou_v8/data/phenotype/v8_data_phenotype_summary_mhwb_survey_summary.txt.bgz') %>% 
      select(-n_cases_ALL, -n_samples_defined_ALL, -n_controls_ALL)%>%
      select(1:3, 8:19) %>%
      mutate(category = 'mhwb_survey'), 
    read_delim('~/Dropbox (Partners HealthCare)/aou_v8/data/phenotype/v8_data_phenotype_summary_r_drug_summary.txt.bgz') %>% 
      select(-n_cases_ALL, -n_samples_defined_ALL, -n_controls_ALL)%>%
      select(1:3, 8:19) %>%
      mutate(category = 'r_drug') ,
    read_delim('~/Dropbox (Partners HealthCare)/aou_v8/data/phenotype/pfhh_survey_summary.txt.bgz')%>%
      select(1:3, 8:19) %>% 
      mutate(category = 'pfhh_survey'),
    read_delim('~/Dropbox (Partners HealthCare)/aou_v8/data/phenotype/physical_measurement_summary.txt.bgz') %>%
      mutate(n_controls=0) %>%
      select(1:2, n_controls, 4:9) %>%
      mutate(n_controls_AFR = 0, n_controls_AMR = 0, n_controls_EAS = 0, n_controls_EUR = 0, n_controls_MID = 0, n_controls_SAS = 0, category = 'physical_measurement'),
    read_delim('~/Dropbox (Partners HealthCare)/aou_v8/data/phenotype/mcc2_phecode_summary.txt.bgz')%>%
      select(1:3, 8:19) %>% 
      mutate(category = 'mcc2_phecode'),
    read_delim('~/Dropbox (Partners HealthCare)/aou_v8/data/phenotype/mcc2_phecodex_summary.txt.bgz')%>%
      select(1:3, 8:19) %>% 
      mutate(category = 'mcc2_phecodex'),
    read_delim('~/Dropbox (Partners HealthCare)/aou_v8/data/phenotype/v8_data_phenotype_summary_random_pheno_binary_summary.txt.bgz') %>%
      select(1:3, 8:19) %>% 
      filter(startsWith(phenoname, 'random_0.5') & grepl('(_1|_2|_3|_4|_5|_1_female|_2_female|_3_female|_4_female|_5_female|_1_male|_2_male|_3_male|_4_male|_5_male)$', phenoname))%>%
      mutate(category = 'random_phenotype'),
    read_delim('~/Dropbox (Partners HealthCare)/aou_v8/data/phenotype/v8_data_phenotype_summary_random_pheno_quantitative_summary.txt.bgz') %>%
      mutate(n_controls=0) %>%
      select(1:2, n_controls, 4:9) %>%
      mutate(n_controls_AFR = 0, n_controls_AMR = 0, n_controls_EAS = 0, n_controls_EUR = 0, n_controls_MID = 0, n_controls_SAS = 0,  category = 'random_phenotype') %>%
      filter(startsWith(phenoname, 'random_0.5') & grepl('(_1|_2|_3|_4|_5|_1_female|_2_female|_3_female|_4_female|_5_female|_1_male|_2_male|_3_male|_4_male|_5_male)$', phenoname))
  ) %>%
    select(-n_cases, -n_controls) 
  
  v8_pheno_info_long <- merge( 
    v8_pheno_info %>%
      dplyr::select(-(8:13)) %>%
      pivot_longer(cols = starts_with('n_cases'), names_to = 'ancestry', values_to = 'n_cases') %>%
      dplyr::mutate(ancestry = str_replace(ancestry, 'n_cases_', '')),
    v8_pheno_info %>%
      dplyr::select(-(2:7)) %>%
      pivot_longer(cols = starts_with('n_controls'), names_to = 'ancestry', values_to = 'n_controls')%>%
      dplyr::mutate(ancestry = str_replace(ancestry, 'n_controls_', '')),
    by = c('phenoname', 'ancestry', 'category')
  ) 
  v8_pheno_info_more <- rbind(
    v8_pheno_info_long,
    read_delim('~/Dropbox (Partners HealthCare)/aou_v8/data/phenotype/onset_summary_processed.txt.bgz') %>%
      mutate(ancestry = toupper(str_split(phenoname, '_') %>% map_chr(., 4)),
             category = 'onset') %>%
      select(phenoname, ancestry, category, n_cases, n_controls),
    read_delim('~/Dropbox (Partners HealthCare)/aou_v8/data/phenotype/progression_summary_processed.txt.bgz') %>%
      mutate(ancestry = toupper(str_split(phenoname, '_') %>% map_chr(., 4)),
             category = 'progression') %>%
      select(phenoname, ancestry, category, n_cases, n_controls),
    read_delim('~/Dropbox (Partners HealthCare)/aou_v8/data/phenotype/lab_measurement_summary.txt.bgz')  %>%
      mutate(n_controls = 0,
             category = 'lab_measurement') %>%
      select(phenoname, ancestry,category, n_cases, n_controls) 
  ) %>%
    mutate(pheno_sex = if_else(grepl('_male|_MALE', phenoname), 'male', if_else(grepl('_female|_FEMALE', phenoname), 'female', 'both'))) 
  return(v8_pheno_info_more)
}


# Fast GCS reader: gsutil -> tempfile -> read -> cleanup
# Supports: csv/tsv/txt, rds, RData, parquet, feather, fst (auto-detect by extension)
read_gcs_fast <- function(gcs_uri, format = NULL, ...) {
  gsutil <- "/Users/wlu/google-cloud-sdk/bin/gsutil"
  if (!file.exists(gsutil)) {
    stop("Real gsutil not found at: ", gsutil,
         "\nInstall Google Cloud SDK or update the path in read_gcs_fast().")
  }
  
  # Figure out extension & tempfile
  uri_path <- sub("^gs://[^/]+/", "", gcs_uri)
  ext <- tools::file_ext(uri_path)
  base_ext <- if (tolower(ext) == "gz") tools::file_ext(tools::file_path_sans_ext(uri_path)) else ext
  fmt <- tolower(if (!is.null(format)) format else if (tolower(ext) == "gz") base_ext else ext)
  tmp <- tempfile(fileext = if (nzchar(ext)) paste0(".", ext) else "")
  
  on.exit(unlink(tmp), add = TRUE)
  
  # Copy down with *real* gsutil (bypasses PATH / wrapper)
  res <- tryCatch(
    system2(gsutil, c("cp", gcs_uri, tmp), stdout = TRUE, stderr = TRUE),
    error = function(e) e
  )
  status <- attr(res, "status")
  if (inherits(res, "error") || (!is.null(status) && status != 0) || !file.exists(tmp)) {
    msg <- paste(c("gsutil copy failed.", if (length(res)) res else NULL), collapse = "\n")
    stop(msg)
  }
  
  # Reader dispatch
  switch(fmt,
         "csv" = { if (!requireNamespace("data.table", quietly = TRUE)) stop("Need data.table"); data.table::fread(tmp, ...) },
         "tsv" = { if (!requireNamespace("data.table", quietly = TRUE)) stop("Need data.table"); data.table::fread(tmp, sep = "\t", ...) },
         "txt" = { if (!requireNamespace("data.table", quietly = TRUE)) stop("Need data.table"); data.table::fread(tmp, ...) },
         "bgz" = { if (!requireNamespace("data.table", quietly = TRUE)) stop("Need data.table"); data.table::fread(tmp, ...) },
         
         "rds" = readRDS(tmp),
         "rdata" = { e <- new.env(); load(tmp, envir = e); objs <- ls(e); if (length(objs) == 1) e[[objs]] else mget(objs, envir = e) },
         
         "parquet" = { if (!requireNamespace("arrow", quietly = TRUE)) stop("Need arrow"); arrow::read_parquet(tmp, ...) },
         "feather" = { if (!requireNamespace("arrow", quietly = TRUE)) stop("Need arrow"); arrow::read_feather(tmp, ...) },
         "fst" = { if (!requireNamespace("fst", quietly = TRUE)) stop("Need fst"); fst::read_fst(tmp, ...) },
         
         stop(sprintf("Unknown or unsupported format '%s'. Pass format= ('csv','tsv','rds','parquet','feather','fst','rdata').", fmt))
  )
}



### figure function
plot_aou_sample_distribution <- function(ver, raw=T, meta = F, type = 'samples', save=T){
  title <- 'Number'
  tag <- 'qced'
  cohort <- ''
  title <- NULL
  limits <- c(0, 285000)
  if(raw){
    aou_data <- aou_data %>%
      mutate(n_samples = n_samples_raw)
    # title <- 'Raw number'
    tag <- 'raw'
  }
  if(ver == 'v7') limits <- c(0, 165000)
  if(!meta) aou_data <- aou_data %>% filter(ancestry  != 'meta')
  if(meta){
    cohort = '_all' 
    limits <- limits*1.75
    title <- paste0('Number of ', str_to_sentence(type),' (', ver, ')')
  }
  if(type == 'samples'){
    aou_data$N <- aou_data$n_samples
  }else{
    aou_data$N <- aou_data$n_pheno
    limits <- c(0, 3900)
  }
  p <- aou_data %>%
    filter(version == ver)%>%
    mutate(ancestry = factor(ancestry, levels = rev(names(pops)))) %>%
    ggplot + aes(x = ancestry, y = N, color = ancestry, fill = ancestry) +
    geom_col(width = 0.5)+
    scale_y_continuous(labels = comma, limits = limits ) +
    scale_alpha_discrete(name = NULL, range = c(0.2, 1)) + 
    scale_x_discrete(breaks = names(pops), labels = pops) + 
    scale_color_manual(name = 'Ancestry', values = pop_colors, breaks = names(pops), labels = pops)+
    scale_fill_manual(name = 'Ancestry', values = pop_colors, breaks = names(pops), labels = pops) + 
    themes + 
    labs(title = NULL, x= NULL, y=NULL) + 
    guides(colour = 'none', fill= 'none') + coord_flip() + 
    geom_text(aes(label =comma(N)),  size = 2, hjust = -0.1) + 
    theme(legend.position = c(0.9, 0.9), 
          axis.text.y = element_text(size = 9, family = 'mono'),
          axis.text.x = element_text(family = 'mono'),
          legend.text = element_text(size=6),
          plot.margin = unit(c(0.1, 0.5, 0, 0), 'cm')
    ) +
    guides(
      alpha = guide_legend(nrow = 2, byrow = TRUE, reverse = F,  keywidth = 1, keyheight = 0.3)
    )
  if(meta) p <- p + theme(axis.title.x = element_blank())
  if(save){
    output_figure(p, 'summary', paste0('aou_', substr(type, 1, nchar(type) - 1),'_', tag,'_summary_', ver, cohort), height=1, width=8/3)
  }
  return(p)
}

obtain_novelty_data_for_plot <- function() {
  join_cols <- c('gene_id', 'gene_symbol', 'annotation', 'phenoname')
  sig_threshold <- 6.7e-7
  
  # Load per-ancestry hits
  ancestries <- c('afr', 'amr', 'eas', 'eur', 'mid', 'sas')
  new_hits <- map_dfr(ancestries, ~read_gcs_fast(
    paste0('gs://aou_wlu/v8_analysis/aou_ukb_meta/aou_v8_', .x, '_interesting_hits_001.txt.bgz')
  )) %>% select(-max_MAF)
  
  # Load phenotype info and meta-analysis signal IDs
  pheno_info <- read_tsv('~/Dropbox (Partners HealthCare)/aou_v8/data/phenotype/aou_v8_phenotype_info_qced_with_meta.tsv')
  pheno_lookup <- pheno_info %>%
    filter(ancestry == 'META', pheno_sex == 'both') %>%
    distinct(phenoname, category, description)
  
  hits_id <- read_gcs_fast('gs://aou_wlu/v8_analysis/aou_ukb_meta/interesting_meta_lof_signals_v8_001.tsv') %>%
    distinct()
  
  # Identify signals significant only in meta-analysis
  hits_info <- read_gcs_fast('gs://aou_wlu/v8_analysis/aou_ukb_meta/significant_in_meta_only_0.001.txt.bgz') %>%
    distinct() %>%
    right_join(hits_id, by = join_cols) %>%
    filter(META_Pvalue_Burden < sig_threshold & ukb_Pvalue_Burden > sig_threshold & aou_Pvalue_Burden > sig_threshold) %>%
    select(-description) %>%
    left_join(pheno_lookup, by = 'phenoname')
  
  # Reshape meta-level summary stats into long format
  cohort_map <- list(
    aou_meta = c(Pvalue_Burden = 'aou_Pvalue_Burden', BETA_Burden = 'aou_BETA_Burden'),
    ukb_meta = c(Pvalue_Burden = 'ukb_Pvalue_Burden', BETA_Burden = 'ukb_BETA_Burden'),
    all_meta = c(Pvalue_Burden = 'META_Pvalue_Burden', BETA_Burden = 'META_Stats_Burden')
  )
  
  formatted_hits_info <- rbind(
    hits_info %>% 
      select(gene_id, gene_symbol, annotation, phenoname, Pvalue_Burden = aou_Pvalue_Burden, BETA_Burden = aou_BETA_Burden) %>%
      mutate(SE_Burden = 0, ancestry = 'aou_meta'),
    hits_info %>% 
      select(gene_id, gene_symbol, annotation, phenoname, Pvalue_Burden = ukb_Pvalue_Burden, BETA_Burden = ukb_BETA_Burden) %>% 
      mutate(SE_Burden = 0, ancestry = 'ukb_meta'),
    hits_info %>% 
      select(gene_id, gene_symbol, annotation, phenoname, Pvalue_Burden = META_Pvalue_Burden, BETA_Burden = META_Stats_Burden) %>%
      mutate(SE_Burden = 0, ancestry = 'all_meta')
  )
  
  # Combine per-ancestry and meta-level results
  full_new_hits <- new_hits %>%
    right_join(hits_info %>% distinct(!!!syms(join_cols)), by = join_cols) %>%
    bind_rows(formatted_hits_info) %>%
    mutate(pvalue = if_else(Pvalue_Burden < sig_threshold, 'Significant', 'Non-significant')) %>%
    left_join(pheno_info %>% distinct(phenoname, description, category), by = 'phenoname') %>%
    mutate(logP = -log10(Pvalue_Burden),
           assoc = paste0(str_to_sentence(description), ' ~ ', gene_symbol))
  
  # Load AI novelty verdicts
  ai_results <- read_tsv('~/Dropbox (Partners HealthCare)/aou_v8/data/novelty_check/novelty_analysis_results_meta_burden_001_lof.tsv') %>%
    transmute(assoc = paste0(str_to_sentence(phenotype), ' ~ ', gene), verdict = verdict_max)
  
  # Write supplementary table
  hits_info %>%
    mutate(assoc = paste0(str_to_sentence(description), ' ~ ', gene_symbol)) %>%
    left_join(ai_results, by = 'assoc') %>%
    write_tsv('~/Dropbox (Partners HealthCare)/aou_v8/tables/supplementary_data7.tsv')
  print(nrow(hits_info))
  
  full_new_hits <- full_new_hits %>%
    left_join(ai_results, by = 'assoc')
  return(full_new_hits)
}

plot_novelty <- function(data, category_name, pop_pattern, indexes = 1:9){
  pops = c('AFR-like', 'AMR-like', 'EAS-like', 'EUR-like', 'MID-like', 'SAS-like', 'AoU META', 'UKB', 'AoU UKB META')
  names(pops) <- c('afr', 'amr', 'eas', 'eur', 'mid','sas', 'aou_meta', 'ukb_meta', 'all_meta')
  label_type = labeller(ancestry = pops, annotation=annotation_names)
  color_mde =  '#EEA9B8'
  
  pop_colors = c('afr' = color_afr,
                 'amr' = color_amr,
                 'eas' = color_eas,
                 'eur' = color_nfe,
                 'mid' = color_mde,
                 'sas' = color_sas,
                 'aou_meta'= '#1E275C',
                 'ukb_meta'= '#007FA3',
                 'all_meta' = color_oth)
  figure_index <- c(lab_measurement='29', physical_measurement='30')
  full_new_hits2 <- data %>%
    filter(ancestry %in% names(pops)[indexes])%>%
    mutate(assoc = if_else(grepl('Vitamin d', assoc), str_replace(assoc, ' d ', ' D '), assoc)) %>%
    mutate(row_color = pop_colors[ancestry],
           row_pattern = 'circle',
           row_background = pop_pattern[ancestry], 
           ancestry = factor(ancestry, levels = names(pops[indexes]), labels = pops[indexes]))  %>%
    dplyr::filter(category == category_name) %>%
    dplyr::arrange(assoc)
  
  # Split assoc into two panels
  unique_assocs <- sort(unique(full_new_hits2$assoc))
  midpoint <- ceiling(length(unique_assocs) / 2)
  panel1_assocs <- unique_assocs[1:midpoint]
  panel2_assocs <- unique_assocs[(midpoint + 1):length(unique_assocs)]
  
  full_new_hits2 <- full_new_hits2 %>%
    mutate(panel = if_else(assoc %in% panel1_assocs, "1", "2"))
  
  
  figure <- full_new_hits2  %>%
    mutate(verdict = factor(verdict, levels=c( 'Novel','Hypothesized', 'Existing', 'Established'))) %>%
    mutate(BETA_Burden = as.numeric(BETA_Burden),
           log_p_scale = case_when(
             Pvalue_Burden > 0.05 ~ '> 0.05',
             Pvalue_Burden > 0.0001 ~ '> 0.0001',
             Pvalue_Burden > 6.7e-7 ~ '> 6.7e-7',
             Pvalue_Burden <= 6.7e-7 ~ 'Significant',
           )) %>%
    mutate(log_p_scale = factor(log_p_scale, levels=rev(c('> 0.05', '> 0.0001', '> 6.7e-7', 'Significant')))) %>%
    ggplot + aes(x = assoc, y = ancestry, fill)  +
    # Add row-wise background rectangles
    geom_rect(
      aes(
        xmin = -Inf, xmax = Inf,
        ymin = as.numeric(as.factor(ancestry)) - 0.5,
        ymax = as.numeric(as.factor(ancestry)) + 0.5
      ),
      fill =  full_new_hits2$row_background,
      # pattern = 'circle', 
      # pattern_fill = full_new_hits$row_color,       # Pattern color
      # pattern_color = full_new_hits$row_color, 
      # pattern_density = 0.1,        # Density of the pattern
      # pattern_size = 0.5,           # Thickness of the pattern
      inherit.aes = FALSE, alpha = 0.1
    ) +
    geom_point(aes(color = log_p_scale), size = 3) +   # Circle size = p-value, color = beta
    geom_point(data = full_new_hits2,
               aes(x = assoc, y = 0.5, fill = verdict),
               shape = 21, size = 3, color = "white") +
    scale_fill_manual(values = c("#D73027", "#F46D43", "#FDAE61", "#FEF6E8"), breaks = c( 'Not found','Hypothesized', 'Existing', 'Established'), name = 'Verdict') +
    scale_color_manual(name = 'Burden Pvalue', values = c("#165B69", "#28A4BD", "#96DBE9", "#EAF8FB"), breaks = rev(c('> 0.05', '> 0.0001', '> 6.7e-7', 'Significant'))) + 
    scale_size_manual(name = NULL, values = c("Significant" = 5, "Non-significant"=3)) + 
    geom_text(aes(label = if_else(BETA_Burden > 0, '+', '-')), color = "gray40", size = 5) +
    labs(x = NULL, y = NULL, color = expression(-log[10]~P)) +
    themes + 
    coord_flip() +
    guides(
      fill = guide_legend(keywidth = 0.5, keyheight = 0.5, order = 1),
      color = guide_legend(keywidth = 0.5, keyheight = 0.5, order = 2),
      size = 'none'
    )+
    facet_wrap(~ panel, scales = "free_y", ncol = 2)  +
    theme(axis.text.x = element_text(colour = pop_colors[indexes], size = 8, face = 'bold', angle = 90, hjust = 1),
          
          axis.ticks.x = element_blank(),
          axis.line.x = element_blank(),
          axis.ticks.y = element_line(linewidth = 0.1),
          axis.line.y = element_line(linewidth = 0.1),
          legend.position = "top",
          legend.box = "horizontal",
          legend.direction = "horizontal",
          legend.justification = "left",
          
          legend.box.margin = margin(0, 0, 0, 0),
          strip.background = element_blank(),
          strip.text = element_blank()
      
    )
  if(category_name == 'mcc2_phecodex'){
    figure <- figure + 
      theme(axis.text.y = element_text(size = 10),
            legend.title = element_text(size = 12),
            legend.text = element_text(size = 12),
            legend.margin = margin(0, 0, 0, -200),
            legend.box.spacing = unit(0.3, "cm"),
            legend.spacing.x = unit(8, "cm"), 
            
      )
    output_figure(figure, 'main', paste0('figure4b_', category_name), width=12, height=10)
  }else{
    figure <- figure + 
      theme(axis.text.y = element_text(size = 10),
            legend.title = element_text(size = 10),
            legend.text = element_text(size = 10),
            legend.margin = margin(0, 0, 0, -150),
            legend.box.spacing = unit(1, "cm"),
            legend.spacing.x = unit(8, "cm"), 
      )
    output_figure(figure, 'supplement', paste0('supp_figure', figure_index[category_name],'_', category_name), width=12, height=10)
  }
  return(figure)
}

