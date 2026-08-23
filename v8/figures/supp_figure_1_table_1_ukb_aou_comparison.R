source('~/Dropbox (Partners HealthCare)/github_repo/aou_gwas/R/constants.R')

pan_ukb_pheno <- read_delim('~/Dropbox (Partners HealthCare)/aou/preliminary/aou_pan_ancestry_phecode_corr.csv') %>%
  select(phenoname = phecode, ancestry, prevalence_pan, n_cases_pan, n_controls_pan)
pheno_full <- read_tsv('~/Dropbox (Partners HealthCare)/aou_v8/data/phenotype/aou_v8_phenotype_info_full.tsv')
pheno_full <- pheno_full %>%
  filter(pheno_sex == 'both') %>%
  select(phenoname, ancestry, category, n_cases, n_controls, in_GWAS) 
pheno_full <- pheno_full %>%
  mutate(prevalence = n_cases /(n_cases + n_controls)) %>%
  group_by(phenoname, category) %>%
  mutate(prevalence_meta = sum(n_cases)/(sum(n_controls) + sum(n_cases)))

p <- pheno_full %>%
  merge(., pan_ukb_pheno %>% mutate(ancestry = if_else(ancestry == 'CSA', 'SAS', ancestry)), by = c('phenoname', 'ancestry'),  all.y=T) %>%
  mutate(ancestry = tolower(ancestry)) %>%
  mutate(ancestry_label = pops[ancestry]) %>%
  ggplot + aes(x = prevalence_pan, y = prevalence, color = ancestry) +
  geom_point() + geom_abline(intercept = 0, slope = 1, lty=2) + 
  scale_color_manual(name = 'Ancestry', values = pop_colors, breaks = names(pops), labels = pops, guide = 'none')+
  labs(x = 'Prevalence (UK Biobank)', y = 'Prevalence (All of Us QCed)') +
  themes + facet_wrap(~ancestry_label, scale = 'free')
p
output_figure(p, 'supplement', 'supp_figure1_prevalence_ukb_vs_aou_phecode', 4, 7.5)

corr <- pheno_full %>%
  filter(category == 'mcc2_phecode') %>%
  merge(., pan_ukb_pheno %>% mutate(ancestry = if_else(ancestry == 'CSA', 'SAS', ancestry)), by = c('phenoname', 'ancestry'))
table(corr$ancestry)
# AFR  AMR  EAS  EUR  MID  SAS 
# 335   31   87 1326   81  410 
corr %>%
  group_by(ancestry) %>%
  dplyr::summarize(tidy(cor.test(n_cases, n_cases_pan))) %>%
  write_csv(., '~/Dropbox (Partners HealthCare)/aou_v8/data/phenotype/aou_ukb_n_cases_corr.csv')
test <- cor.test(corr$n_cases, corr$n_cases_pan)
test
test$p.value

corr %>%
  group_by(ancestry) %>%
  dplyr::summarize(tidy(cor.test(n_controls, n_controls_pan))) %>%
  write_csv(., '~/Dropbox (Partners HealthCare)/aou_v8/data/phenotype/aou_ukb_n_controls_corr.csv')
test <- cor.test(corr$n_controls, corr$n_controls_pan)
test
test$p.value

corr %>%
  group_by(ancestry) %>%
  dplyr::summarize(tidy(cor.test(prevalence, prevalence_pan))) %>%
  write_csv(., '~/Dropbox (Partners HealthCare)/aou_v8/data/phenotype/aou_ukb_prevalence_corr.csv')
test <- cor.test(corr$prevalence, corr$prevalence_pan)
test$p.value



# p <- pheno_full %>%
#   filter(prevalence < 1) %>%
#   mutate(ancestry = tolower(ancestry)) %>%
#   ggplot + aes(x = prevalence_meta, y = prevalence, color = ancestry) +
#   geom_point() + geom_abline(intercept = 0, slope = 1, lty=2) + 
#   scale_color_manual(name = 'Ancestry', values = pop_colors, breaks = names(pops), labels = pops, guide = 'none')+
#   labs(x = 'Prevalence (Meta)', y = 'Prevalence (Ancestry-specific)') +
#   themes + facet_wrap(~ancestry, labeller=label_type)
# p
# output_figure(p, 'summary', 'pheno_prevalence_comparison_meta_vs_ancestry_qced', 4, 7)

# p <- pheno_full %>%
#   filter(prevalence < 1 & in_GWAS) %>%
#   mutate(ancestry = tolower(ancestry)) %>%
#   ggplot + aes(x = prevalence_meta, y = prevalence, color = ancestry) +
#   geom_point() + geom_abline(intercept = 0, slope = 1, lty=2) + 
#   scale_color_manual(name = 'Ancestry', values = pop_colors, breaks = names(pops), labels = pops, guide = 'none')+
#   labs(x = 'Prevalence (Meta)', y = 'Prevalence (Ancestry-specific)') +
#   themes + facet_wrap(~ancestry, labeller=label_type)
# p
# output_figure(p, 'summary', 'pheno_prevalence_comparison_meta_vs_ancestry_qced_filtered', 4, 7)

# p <- pheno_full %>%
#   merge(., pan_ukb_pheno %>% mutate(ancestry = if_else(ancestry == 'CSA', 'SAS', ancestry)), by = c('phenoname', 'ancestry'),  all.y=T) %>%
#   mutate(ancestry = tolower(ancestry)) %>%
#   ggplot + aes(x = n_cases_pan + n_controls_pan, y = n_cases + n_controls, color = ancestry) +
#   geom_point() + geom_abline(intercept = 0, slope = 1, lty=2) + 
#   scale_color_manual(name = 'Ancestry', values = pop_colors, breaks = names(pops), labels = pops, guide = 'none')+
#   labs(x = 'N (UKB)', y = 'N (AoU QCed)') +
#   themes + facet_wrap(~ancestry, labeller=label_type, scale = 'free')
# output_figure(p, 'summary', 'pheno_sample_size_comparison_pan_ukb_vs_aou_qced', 4, 7.5)