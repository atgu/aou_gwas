source('~/Dropbox (Partners HealthCare)/github_repo/aou_gwas/R/constants.R')
figure_root <- paste0(aou_qc_fig_path, 'cost/')
data_root <- paste0(aou_qc_data_path, 'cost/')
# install.packages('gg.layers')
library(gg.layers)
# install.packages("remotes")
# remotes::install_github("dill/emoGG")
# library(emoGG)

# Search for a emoji
# https://unicode.org/emoji/charts/full-emoji-list.html

ancestries <- c('afr', 'amr', 'eas', 'eur', 'mid', 'sas')
large_ancestries <- c('afr', 'amr', 'eas', 'eur')
saige_cost <- data.frame()
for(test in c('saige', 'saige_gene')){
  variant_type <- if_else(test=='saige', 'ACAF', 'exome')
  for(pop in ancestries){
    for(sex in c('both', 'male', 'female')){
      cost = read_delim(paste0(data_root, pop,'_', test,'_cost_results_final_', sex,'.csv'), delim = '\t') %>% 
        mutate(variant_type = variant_type, pheno_sex = sex) %>%
        select("ancestry", "phenotype", "analysis_type", 'variant_type', 'pheno_sex', "states", "cost", "total_cost", "duration", "cost_Success", "count_Success", 'desired_cost')
      saige_cost <- rbind(saige_cost, cost)
      if((pop %in% large_ancestries) & (test == 'saige') & (sex == 'both')){
        cost = read_delim(paste0(data_root, pop,'_', test,'_mhwb_cost_results_final_', sex,'.csv'), delim = '\t') %>% 
          mutate(variant_type = variant_type, pheno_sex = sex) %>%
          select("ancestry", "phenotype", "analysis_type", 'variant_type', 'pheno_sex', "states", "states", "cost", "total_cost", "duration", "cost_Success", "count_Success", 'desired_cost')
        saige_cost <- rbind(saige_cost, cost)
      }
    }
  }
}

saige_cost %>%
  filter(!grepl(pattern = '_to_', x = phenotype))%>%
  group_by(ancestry) %>%
  dplyr::summarize(avg_cpu_time_per_job = mean(duration/60/count_Success, na.rm=T),
                   total_cost = sum(total_cost),
  )

cost_full <- saige_cost %>%
  filter(!grepl(pattern = '_to_', x = phenotype)) %>%
  mutate(variant_type = if_else(analysis_type == 'gene', 'gene', variant_type)) %>%
  mutate(variant_type = paste('SAIGE', variant_type, if_else(variant_type == 'gene', '', 'variant'))) 
write_csv(cost_full, paste0(data_root, 'aou_v8_batch_cost_full.csv'))


cost_full <- read_csv(paste0(data_root, 'aou_v8_batch_cost_full.csv'))

plot_cost_per_phenotype <- function(sex, data){
  if(sex == 'both'){
    data <- data %>% filter(pheno_sex == sex)
  }else{
    data <- data %>% filter(pheno_sex != 'both')
  }
  total <- data %>%
    group_by(variant_type) %>%
    dplyr::summarize(max_cost = max(total_cost),
                     total_cost = sum(total_cost),
                     total_cpu_hrs = sum(duration/3600),
                     ) %>%
    mutate(variant_type = factor(variant_type, levels = c('SAIGE ACAF variant', 'SAIGE exome variant', 'SAIGE gene'), labels = c('SAIGE ACAF variant test', 'SAIGE exome variant test', 'SAIGE-GENE+ gene-based test')))
  p <- data %>%
    mutate(analysis_type = if_else(analysis_type == 'gene', 'SAIGE-GENE', 'SAIGE') ) %>%
    mutate(variant_type = factor(variant_type, levels = c('SAIGE ACAF variant', 'SAIGE exome variant', 'SAIGE gene'), labels = c('SAIGE ACAF variant test', 'SAIGE exome variant test', 'SAIGE-GENE+ gene-based test'))) %>%
    # mutate(pheno_sex = factor(pheno_sex, levels = c('both', 'female', 'male'), labels = c('Both sexes', 'Female only', 'Male only'))) %>%
    ggplot + aes(x = ancestry, y = total_cost, color = ancestry) +
    # add_emoji(emoji = "1f4b8")+
    labs(x = NULL, y= 'Cost\n(USD per phenotype)', color='Ancestry') +
    scale_y_log10(label = comma) +
    scale_x_discrete(breaks = c('afr', 'amr', 'eas', 'eur', 'mid', 'sas','oth'), labels = c('AFR-like', 'AMR-like', 'EAS-like', 'EUR-like', 'MID-like', 'SAS-like', 'OTH')) +
    scale_color_manual(name = 'Ancestry',values = pop_colors, breaks = c('afr', 'amr', 'eas', 'eur', 'mid', 'sas','oth'), labels = c('AFR-like', 'AMR-like', 'EAS-like', 'EUR-like', 'MID-like', 'SAS-like', 'OTH'))+
    scale_fill_manual(name = 'Ancestry',values = pop_colors, breaks = c('afr', 'amr', 'eas',  'eur','mid', 'sas','oth'), labels = c('AFR-like', 'AMR-like', 'EAS-like', 'EUR-like', 'MID-like', 'SAS-like', 'OTH')) +
    geom_boxplot2() +
    themes +
    facet_wrap(~variant_type, labeller=label_type, scale = 'free') +
    theme(legend.position = 'none', strip.text.x = element_text(face='bold'), axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    geom_text_repel(
      data = total %>% filter(!is.na(total_cost)),
      aes(x = 0.5, y = max_cost*0.6, label = paste0(comma(total_cost, accuracy = 1), ' USD')),
      vjust = 0.5, fontface = 'bold', color = 'black',
      inherit.aes = FALSE
    )+
    guides(color = guide_legend(nrow = 1)) 
  
  output_figure(p, 'cost', paste0('aou_v8_cost_per_phenotype_', sex), height = 2.5, width = 10)
  return(p)
}

cost_p1 <- plot_cost_per_phenotype(sex='both', data=cost_full)
cost_p2 <- plot_cost_per_phenotype(sex='sex_specific', data=cost_full)


plot_time_per_phenotype <- function(sex, data){
  if(sex == 'both'){
    data <- data %>% filter(pheno_sex == sex)
  }else{
    data <- data %>% filter(pheno_sex != 'both')
  }
  total <- data %>%
    group_by(variant_type) %>%
    dplyr::summarize(max_cpu = max(duration/60),
                     total_cost = sum(total_cost),
                     total_cpu_hrs = sum(duration/3600))%>%
    mutate(variant_type = factor(variant_type, levels = c('SAIGE ACAF variant', 'SAIGE exome variant', 'SAIGE gene'), labels = c('SAIGE ACAF variant test', 'SAIGE exome variant test', 'SAIGE-GENE+ gene-based test')))
  p <- data %>%
    mutate(analysis_type = if_else(analysis_type == 'gene', 'SAIGE-GENE', 'SAIGE') ) %>%
    mutate(variant_type = factor(variant_type, levels = c('SAIGE ACAF variant', 'SAIGE exome variant', 'SAIGE gene'), labels = c('SAIGE ACAF variant test', 'SAIGE exome variant test', 'SAIGE-GENE+ gene-based test'))) %>%
    # mutate(pheno_sex = factor(pheno_sex, levels = c('both', 'female', 'male'), labels = c('Both sexes', 'Female only', 'Male only'))) %>%
    ggplot + aes(x = ancestry, y = duration/60, color = ancestry) +
    # add_emoji(emoji = "23f0")+
    labs(x = NULL, y= 'CPU time\n(minute per phenotype)', color='Ancestry') +
    scale_y_log10(label = comma) +
    scale_x_discrete(breaks = c('afr', 'amr', 'eas', 'eur', 'mid', 'sas','oth'), labels = c('AFR-like', 'AMR-like', 'EAS-like', 'EUR-like', 'MID-like', 'SAS-like', 'OTH')) +
    scale_color_manual(name = 'Ancestry',values = pop_colors, breaks = c('afr', 'amr', 'eas', 'eur', 'mid', 'sas','oth'), labels = c('AFR-like', 'AMR-like', 'EAS-like', 'EUR-like', 'MID-like', 'SAS-like', 'OTH'))+
    scale_fill_manual(name = 'Ancestry',values = pop_colors, breaks = c('afr', 'amr', 'eas',  'eur','mid', 'sas','oth'), labels = c('AFR-like', 'AMR-like', 'EAS-like', 'EUR-like', 'MID-like', 'SAS-like', 'OTH')) +
    geom_boxplot2() +
    themes +
    facet_wrap(~variant_type, labeller=label_type, scale = 'free') +
    theme(legend.position = 'none', strip.background =  element_blank(), strip.text = element_blank()) +
    geom_text_repel(
      data = total %>% filter(!is.na(total_cpu_hrs)),
      aes(x = 0.5, y = max_cpu*0.6, label = paste0(comma(total_cpu_hrs, accuracy = 1), ' CPU hrs')),
      vjust = 0.5, fontface = 'bold', color = 'black',
      inherit.aes = FALSE
    )+
    guides(color = guide_legend(nrow = 1)) 
  
  output_figure(p, 'cost', paste0('aou_v8_cpu_mins_per_phenotype_', sex), height = 2.5, width = 10)
  return(p)
}


time_p1 <- plot_time_per_phenotype(sex='both', data=cost_full)
time_p2 <- plot_time_per_phenotype(sex='sex_specific', data=cost_full)

p = ggpubr::ggarrange(cost_p1, time_p1, nrow=2, hjust = 0, align = 'v',
                      font.label = list(size = 10, color = "black", face = "bold", family = 'Arial')
)
output_figure(p, folder = 'supplement', 'supp_figure11_aou_v8_both_cost_summary', 4, 10)

# p = ggpubr::ggarrange(cost_p2, time_p2, nrow=2, hjust = 0, align = 'v',
#                       font.label = list(size = 10, color = "black", face = "bold", family = 'Arial')
# )
# output_figure(p, folder = 'cost', 'aou_v8_sex_specific_cost_summary', 4, 10)
# 
# 
# cost_full <- read_csv(paste0(data_root, 'aou_v8_batch_cost_full.csv')) %>%
#   mutate(phenoname = phenotype) %>%
#   merge(., read_tsv('~/Dropbox (Partners HealthCare)/0_very_often_used/aou_v8_unique_phenotype_qced.tsv'), by = 'phenoname') %>%
#   filter(pheno_sex == 'both') %>%
#   mutate(trait_type = if_else(category %in% c('lab_measurement', 'physical_measurement') | grepl('continuous', phenoname), 'quantitative', 'binary'))
# 
# total <- cost_full %>%
#   group_by(variant_type, trait_type) %>%
#   dplyr::summarize(total_cost = sum(total_cost),
#                    total_cpu_hrs = sum(duration/3600))
# p <- cost_full %>%
#   mutate(analysis_type = if_else(analysis_type == 'gene', 'SAIGE-GENE', 'SAIGE') ) %>%
#   mutate(variant_type = factor(variant_type, levels = c('SAIGE ACAF variant', 'SAIGE exome variant', 'SAIGE gene'))) %>%
#   # mutate(pheno_sex = factor(pheno_sex, levels = c('both', 'female', 'male'), labels = c('Both sexes', 'Female only', 'Male only'))) %>%
#   ggplot + aes(x = ancestry, y = total_cost, color = ancestry, alpha = factor(trait_type), fill=ancestry) +
#   # add_emoji(emoji = "1f4b8")+
#   labs(x = NULL, y= 'Cost\n(USD per phenotype)', color='Ancestry') +
#   scale_y_log10(label = comma) +
#   scale_x_discrete(breaks = c('afr', 'amr', 'eas', 'eur', 'mid', 'sas','oth'), labels = c('AFR', 'AMR', 'EAS', 'EUR', 'MID', 'SAS', 'OTH')) +
#   scale_color_manual(name = 'Ancestry',values = pop_colors, breaks = c('afr', 'amr', 'eas', 'eur', 'mid', 'sas','oth'), labels = c('AFR', 'AMR', 'EAS', 'EUR', 'MID', 'SAS', 'OTH'))+scale_fill_manual(name = 'Ancestry',values = pop_colors, breaks = c('afr', 'amr', 'eas',  'eur','mid', 'sas','oth'), labels = c('AFR', 'AMR', 'EAS', 'EUR', 'MID', 'SAS', 'OTH')) +
#   geom_boxplot2(position='dodge') +
#   scale_alpha_discrete(name = 'Trait type', range = c(0, 0.8)) +
#   themes +
#   facet_wrap(~variant_type, labeller=label_type, scale = 'free') +
#   theme(legend.position = 'top', strip.text.x = element_text(face='bold'), axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
#   # geom_text_repel(data = total, aes(x = 0.5, y = Inf, label = paste0(comma(total_cost, accuracy = 1), ' USD'), vjust = 0.5, fontface = 'bold'), color = 'black') + 
#   guides(color = guide_legend(nrow = 1),    alpha = guide_legend(
#     override.aes = list(fill = "grey40", color = "black") # ONLY legend gets fill, so alpha is visible
#   )) 
# 
# output_figure(p, 'cost', paste0('aou_v8_cost_per_phenotype_per_trait_type_both'), height = 4, width = 10)

######### [Archived] cost analysis during v8 pilot #########
# cost_full %>%
#   group_by(ancestry, analysis_type, pheno_sex) %>%
#   dplyr::summarize(cost_per_pheno = mean(total_cost)) %>%
#   group_by(ancestry) %>%
#   dplyr::summarize(mean_cost_pheno = mean(cost_per_pheno))
# 
# cost_v4 <- read_csv('~/Desktop/cost_pilot_v4.csv') %>%
#   mutate(mins_per_job = duration/60/count_Success,
#          cost_per_job = total_cost/count_Success,
#          cost_per_pheno = sum(total_cost/count_Success))
# cost_v2 <- read_csv('~/Dropbox (Partners HealthCare)/aou_v8/pipeline_analysis/pilot_cost_saige_1.4.4/cost_pilot_v2.csv') %>%
#   mutate(mins_per_job = duration/60/count_Success,
#          cost_per_job = total_cost/count_Success,
#          cost_per_pheno = sum(total_cost/count_Success))
# cost_full <- cost_v2 %>%
#   merge(., cost_v4, by = c('ancestry', 'phenotype', 'analysis_type', 'variant_type'), all.x=T, suffixes = c('.v2', '.v3'))
# 
# p <- cost_full %>%
#   # merge(., v8_pheno_info_final %>% select(phenoname, category) %>% distinct(), by.x = 'phenotype', by.y = 'phenoname') %>%
#   # filter(ancestry != 'eur') %>%
#   mutate(analysis_type = if_else(analysis_type == 'gene', 'SAIGE-GENE', 'SAIGE') ) %>%
#   mutate(variant_type = factor(variant_type, levels = c('SAIGE genome variant', 'SAIGE exome variant', 'SAIGE gene'))) %>%
#   ggplot + aes(x = desired_cost.v2, y = desired_cost.v3, color = ancestry) +
#   labs(x = 'Last step cost (1.4.4 bgen)', y= 'Last step cost (1.4.8 pgen)', color='Ancestry') +
#   # scale_y_log10(label = comma) +
#   # scale_x_discrete(breaks = c('afr', 'amr', 'eas', 'eur', 'mid', 'sas','oth'), labels = c('AFR', 'AMR', 'EAS', 'EUR', 'MID', 'SAS', 'OTH')) +
#   scale_color_manual(name = 'Ancestry',values = pop_colors, breaks = c('afr', 'amr', 'eas', 'eur', 'mid', 'sas','oth'), labels = c('AFR', 'AMR', 'EAS', 'EUR', 'MID', 'SAS', 'OTH'))+scale_fill_manual(name = 'Ancestry',values = pop_colors, breaks = c('afr', 'amr', 'eas',  'eur','mid', 'sas','oth'), labels = c('AFR', 'AMR', 'EAS', 'EUR', 'MID', 'SAS', 'OTH')) +
#   # geom_boxplot(alpha=0.5) +
#   geom_abline(slope = 1, intercept = 0, lty = 2) +
#   geom_point() +
#   themes +
#   facet_wrap(~variant_type, labeller=label_type) +
#   theme(legend.position = 'top') +
#   guides(color = guide_legend(nrow = 1))
# p
# png(paste0("~/Desktop/cost_per_phenotype_comparison.png"), width=8, height=3, units = 'in', res = 300)
# print(p)
# dev.off()
# 
# 
# library(gg.layers)
# p <- cost_full %>%
#   # merge(., v8_pheno_info_final %>% select(phenoname, category) %>% distinct(), by.x = 'phenotype', by.y = 'phenoname') %>%
#   # filter(ancestry != 'eur') %>%
#   mutate(analysis_type = if_else(analysis_type == 'gene', 'SAIGE-GENE', 'SAIGE') ) %>%
#   mutate(variant_type = factor(variant_type, levels = c('SAIGE genome variant', 'SAIGE exome variant', 'SAIGE gene '))) %>%
#   ggplot + aes(x = ancestry, y = mins_per_job, color = ancestry) +
#   labs(x = NULL, y= 'Approximate minutes per job', color='Population') +
#   scale_y_log10(label = comma) +
#   scale_x_discrete(breaks = c('afr', 'amr', 'eas', 'eur', 'mid', 'sas','oth'), labels = c('AFR', 'AMR', 'EAS', 'EUR', 'MID', 'SAS', 'OTH')) +
#   scale_color_manual(name = 'Ancestry',values = pop_colors, breaks = c('afr', 'amr', 'eas', 'eur', 'mid', 'sas','oth'), labels = c('AFR', 'AMR', 'EAS', 'EUR', 'MID', 'SAS', 'OTH'))+scale_fill_manual(name = 'Ancestry',values = pop_colors, breaks = c('afr', 'amr', 'eas',  'eur','mid', 'sas','oth'), labels = c('AFR', 'AMR', 'EAS', 'EUR', 'MID', 'SAS', 'OTH')) +
#   # geom_boxplot(alpha=0.5) +
#   geom_boxplot(alpha=0.5, ) +
#   themes +
#   facet_wrap(~variant_type, labeller=label_type, scale='free') +
#   theme(legend.position = 'top', axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
#   guides(color = guide_legend(nrow = 1))
# p
# png(paste0("~/Desktop/minutes_per_job_3_test.png"), width=8, height=3, units = 'in', res = 300)
# print(p)
# dev.off()
# 
# 
# p <- cost_full %>%
#   # merge(., v8_pheno_info_final %>% select(phenoname, category) %>% distinct(), by.x = 'phenotype', by.y = 'phenoname') %>%
#   # filter(ancestry != 'eur') %>%
#   mutate(analysis_type = if_else(analysis_type == 'gene', 'SAIGE-GENE', 'SAIGE')) %>%
#   mutate(variant_type = factor(variant_type, levels = c('SAIGE genome variant', 'SAIGE exome variant', 'SAIGE gene '))) %>%
#   ggplot + aes(x = ancestry, y = cost_per_job, color = ancestry) +
#   labs(x = NULL, y= 'Approximate cost per job', color='Ancestry') +
#   scale_y_log10(label = comma) +
#   scale_x_discrete(breaks = c('afr', 'amr', 'eas', 'eur', 'mid', 'sas','oth'), labels = c('AFR', 'AMR', 'EAS', 'EUR', 'MID', 'SAS', 'OTH')) +
#   scale_color_manual(name = 'Ancestry',values = pop_colors, breaks = c('afr', 'amr', 'eas', 'eur', 'mid', 'sas','oth'), labels = c('AFR', 'AMR', 'EAS', 'EUR', 'MID', 'SAS', 'OTH'))+scale_fill_manual(name = 'Ancestry',values = pop_colors, breaks = c('afr', 'amr', 'eas',  'eur','mid', 'sas','oth'), labels = c('AFR', 'AMR', 'EAS', 'EUR', 'MID', 'SAS', 'OTH')) +
#   # geom_boxplot(alpha=0.5) +
#   geom_boxplot(alpha=0.5, ) +
#   themes +
#   facet_wrap(~variant_type, labeller=label_type, scale='free') +
#   theme(legend.position = 'top', axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
#   guides(color = guide_legend(nrow = 1))
# p
# png(paste0("~/Desktop/cost_per_job_3_test.png"), width=8, height=3, units = 'in', res = 300)
# print(p)
# dev.off()
# 
# 
# p <- cost_full %>%
#   # merge(., v8_pheno_info_final %>% select(phenoname, category) %>% distinct(), by.x = 'phenotype', by.y = 'phenoname') %>%
#   # filter(ancestry != 'eur') %>%
#   mutate(analysis_type = if_else(analysis_type == 'gene', 'SAIGE-GENE', 'SAIGE')) %>%
#   mutate(variant_type = factor(variant_type, levels = c('SAIGE genome variant', 'SAIGE exome variant', 'SAIGE gene '))) %>%
#   ggplot + aes(x = ancestry, y = total_cost, color = ancestry) +
#   labs(x = NULL, y= 'Approximate cost per phenotype', color='Ancestry') +
#   scale_y_log10(label = comma) +
#   scale_x_discrete(breaks = c('afr', 'amr', 'eas', 'eur', 'mid', 'sas','oth'), labels = c('AFR', 'AMR', 'EAS', 'EUR', 'MID', 'SAS', 'OTH')) +
#   scale_color_manual(name = 'Ancestry',values = pop_colors, breaks = c('afr', 'amr', 'eas', 'eur', 'mid', 'sas','oth'), labels = c('AFR', 'AMR', 'EAS', 'EUR', 'MID', 'SAS', 'OTH'))+scale_fill_manual(name = 'Ancestry',values = pop_colors, breaks = c('afr', 'amr', 'eas',  'eur','mid', 'sas','oth'), labels = c('AFR', 'AMR', 'EAS', 'EUR', 'MID', 'SAS', 'OTH')) +
#   # geom_boxplot(alpha=0.5) +
#   geom_boxplot(alpha=0.5, ) +
#   themes +
#   facet_wrap(~variant_type, labeller=label_type, scale='free') +
#   theme(legend.position = 'top', axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
#   guides(color = guide_legend(nrow = 1))
# p
# png(paste0("~/Desktop/cost_per_pheno_3_test.png"), width=8, height=3, units = 'in', res = 300)
# print(p)
# dev.off()
# 
# 
# p <- saige_cost %>%
#   merge(., v8_pheno_info_final %>% select(phenoname, category) %>% distinct(), by.x = 'phenotype', by.y = 'phenoname') %>%
#   filter(ancestry != 'eur') %>%
#   mutate(analysis_type = if_else(analysis_type == 'gene', 'SAIGE-GENE', 'SAIGE'),
#          mins_per_job = duration/60/count_Success) %>%
#   mutate(variant_type = factor(variant_type, levels = c('SAIGE genome variant', 'SAIGE exome variant', 'SAIGE gene '))) %>%
#   ggplot + aes(x = ancestry, y = mins_per_job, color = ancestry) +
#   labs(x = NULL, y= 'Approximate minutes per job', color='Population') +
#   scale_y_log10(label = comma) +
#   scale_x_discrete(breaks = c('afr', 'amr', 'eas', 'eur', 'mid', 'sas','oth'), labels = c('AFR', 'AMR', 'EAS', 'EUR', 'MID', 'SAS', 'OTH')) +
#   scale_color_manual(name = 'Ancestry',values = pop_colors, breaks = c('afr', 'amr', 'eas', 'eur', 'mid', 'sas','oth'), labels = c('AFR', 'AMR', 'EAS', 'EUR', 'MID', 'SAS', 'OTH'))+scale_fill_manual(name = 'Ancestry',values = pop_colors, breaks = c('afr', 'amr', 'eas',  'eur','mid', 'sas','oth'), labels = c('AFR', 'AMR', 'EAS', 'EUR', 'MID', 'SAS', 'OTH')) +
#   # geom_boxplot(alpha=0.5) +
#   geom_boxplot(alpha=0.5, ) +
#   themes +
#   facet_wrap(~category, labeller=label_type, scale='free', nrow = 1) +
#   theme(legend.position = 'top', axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
#   guides(color = guide_legend(nrow = 1))
# p
# png(paste0("~/Desktop/minutes_per_job_by_cat.png"), width=8, height=3, units = 'in', res = 300)
# print(p)
# dev.off()
# 
# 
# p <- saige_cost %>%
#   merge(., v8_pheno_info_final %>% select(phenoname, category) %>% distinct(), by.x = 'phenotype', by.y = 'phenoname') %>%
#   filter(ancestry != 'eur') %>%
#   mutate(analysis_type = if_else(analysis_type == 'gene', 'SAIGE-GENE', 'SAIGE'),
#          cost_per_job = total_cost/count_Success) %>%
#   mutate(variant_type = factor(variant_type, levels = c('SAIGE genome variant', 'SAIGE exome variant', 'SAIGE gene '))) %>%
#   ggplot + aes(x = ancestry, y = cost_per_job, color = ancestry) +
#   labs(x = NULL, y= 'Approximate cost per job', color='Ancestry') +
#   scale_y_log10(label = comma) +
#   scale_x_discrete(breaks = c('afr', 'amr', 'eas', 'eur', 'mid', 'sas','oth'), labels = c('AFR', 'AMR', 'EAS', 'EUR', 'MID', 'SAS', 'OTH')) +
#   scale_color_manual(name = 'Ancestry',values = pop_colors, breaks = c('afr', 'amr', 'eas', 'eur', 'mid', 'sas','oth'), labels = c('AFR', 'AMR', 'EAS', 'EUR', 'MID', 'SAS', 'OTH'))+scale_fill_manual(name = 'Ancestry',values = pop_colors, breaks = c('afr', 'amr', 'eas',  'eur','mid', 'sas','oth'), labels = c('AFR', 'AMR', 'EAS', 'EUR', 'MID', 'SAS', 'OTH')) +
#   # geom_boxplot(alpha=0.5) +
#   geom_boxplot(alpha=0.5, ) +
#   themes +
#   # facet_wrap(~category, labeller=label_type, scale='free') +
#   theme(legend.position = 'top', axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
#   guides(color = guide_legend(nrow = 1))
# p
# png(paste0("~/Desktop/cost_per_job.png"), width=5, height=3, units = 'in', res = 300)
# print(p)
# dev.off()
# 
# 
# p <- saige_cost %>%
#   merge(., v8_pheno_info_final %>% select(phenoname, category) %>% distinct(), by.x = 'phenotype', by.y = 'phenoname') %>%
#   filter(ancestry != 'eur') %>%
#   mutate(analysis_type = if_else(analysis_type == 'gene', 'SAIGE-GENE', 'SAIGE'),
#          cost_per_pheno = sum(total_cost/count_Success),
#          category = factor(category )) %>%
#   mutate(variant_type = factor(variant_type, levels = c('SAIGE genome variant', 'SAIGE exome variant', 'SAIGE gene '))) %>%
#   ggplot + aes(x = ancestry, y = total_cost, color = ancestry) +
#   labs(x = NULL, y= 'Approximate cost per phenotype', color='Ancestry') +
#   scale_y_log10(label = comma) +
#   scale_x_discrete(breaks = c('afr', 'amr', 'eas', 'eur', 'mid', 'sas','oth'), labels = c('AFR', 'AMR', 'EAS', 'EUR', 'MID', 'SAS', 'OTH')) +
#   scale_color_manual(name = 'Ancestry',values = pop_colors, breaks = c('afr', 'amr', 'eas', 'eur', 'mid', 'sas','oth'), labels = c('AFR', 'AMR', 'EAS', 'EUR', 'MID', 'SAS', 'OTH'))+scale_fill_manual(name = 'Ancestry',values = pop_colors, breaks = c('afr', 'amr', 'eas',  'eur','mid', 'sas','oth'), labels = c('AFR', 'AMR', 'EAS', 'EUR', 'MID', 'SAS', 'OTH')) +
#   # geom_boxplot(alpha=0.5) +
#   geom_boxplot(alpha=0.5, ) +
#   themes +
#   facet_wrap(~category, labeller=label_type, scale='free', nrow = 1) +
#   theme(legend.position = 'top', axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
#   guides(color = guide_legend(nrow = 1))
# p
# png(paste0("~/Desktop/cost_per_pheno_by_cat.png"), width=8, height=3, units = 'in', res = 300)
# print(p)
# dev.off()
# 
# View(saige_cost %>%
#        merge(., v8_pheno_info_final %>% select(phenoname, category) %>% distinct(), by.x = 'phenotype', by.y = 'phenoname')%>%
#        group_by(category, ancestry) %>%
#        dplyr::summarize(mean_min_per_job = mean(duration/60/count_Success)))