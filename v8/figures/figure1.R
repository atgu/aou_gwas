source('~/Dropbox (Partners HealthCare)/github_repo/aou_gwas/R/constants.R')
figure_root <- paste0(aou_qc_fig_path, 'main/')
threshold = 6.7e-7
library(magick)

df_both <- tribble(
  ~Category,               ~`6`,  ~`5`,   ~`4`,   ~`3`,   ~`2`,   ~`1`,
  "Lab measurements",       57,    12,     4,      5,      4,      3,
  "PhecodeX",          9,    76,     75,     815,    169,    579,
  "Prescriptions",                 205,    446,    163,     676,    31,     79,
  "Self-reported",            1,    20,     15,     69,     9,     25,
  "Physical measurements",  10,    NA,     NA,     NA,     NA,     NA,
) # both

plot_phenotype_by_category <- function(sex, df, save = F){
  category_colors <- c('#334195', '#d56f3e', '#7e57c2', '#4f345a','#880d1e', '#43aa8b', '#D470A2', '#FFE584', 'gray')
  names(category_colors) <- c('Lab measurements', 'Physical measurements', 'Phecode', 'PhecodeX', 'Prescriptions', 'Self-reported', 'Onset', 'Progression', 'Random phenotype')
  category_labels <- c('Lab measurements', 'Physical measurements', 'Phecode', 'PhecodeX', 'Prescriptions', 'Self-reported', 'Onset', 'Progression', 'Random phenotype')
  categories <- c('Lab measurements', 'Physical measurements', 'Phecode', 'PhecodeX', 'Prescriptions', 'Self-reported', 'Onset', 'Progression', 'Random phenotype')
  category_color_scale<- scale_color_manual(name = NULL, values = category_colors, breaks = names(category_colors), labels = category_labels) 
  category_fill_scale <- scale_fill_manual(name = NULL, values = category_colors, breaks = names(category_colors), labels = category_labels)
  
  n_qced_phenotypes <- read_tsv('~/Dropbox (Partners HealthCare)/aou_v8/data/phenotype/aou_v8_phenotype_info_qced_with_meta.tsv') %>%
    filter(!category %in% c('onset', 'progression') & pheno_sex == 'both')%>%
    mutate(category = if_else(category %in% c('pfhh_survey', 'mhwb_survey'), 'self-reported', category)) %>%
    select(phenoname, category) %>%
    distinct() %>%
    group_by(category) %>%
    dplyr::summarize(n = n()) %>%
    ungroup() %>%
    merge(., data.frame(category = c('lab_measurement', 'physical_measurement', 'mcc2_phecodex', 'r_drug', 'self-reported'),
                        n_raw = c(86, 10, 3426, 4521, 145)) , by = 'category'
          ) %>% 
    mutate(category = factor(category, levels = rev(c('physical_measurement', 'lab_measurement', 'self-reported', 
                                                  'r_drug', 'mcc2_phecodex')),
                             labels = rev(c('Physical measurements', 'Lab measurements', 'Self-reported',
                                        'Prescriptions', 'PhecodeX'))))
  
  n_phenos_by_group <-  read_tsv('~/Dropbox (Partners HealthCare)/aou_v8/data/phenotype/aou_v8_phenotype_info_qced_with_meta.tsv') %>%
    filter(!category %in% c('onset', 'progression')) %>%
    mutate(category = if_else(category %in% c('pfhh_survey', 'mhwb_survey'), 'self-reported', category))%>%
    dplyr::group_by(category, ancestry, pheno_sex) %>%
    dplyr::summarize(
      n_phenos = n()
    )  %>%
    mutate(category = factor(category, levels = c('lab_measurement', 'physical_measurement', 'mcc2_phecode', 
                                                  'mcc2_phecodex', 'r_drug', 'self-reported', 'onset', 'progression', 'random_phenotype'),
                             labels = c('Lab measurements', 'Physical measurements', 'Phecode', 
                                        'PhecodeX', 'Prescriptions', 'Self-reported',
                                        'Onset', 'Progression', 'Random phenotype'))) %>%
    group_by(ancestry, pheno_sex) %>%
    mutate(n_total = sum(n_phenos))
  
  p_left<- n_phenos_by_group %>%
    filter(pheno_sex == sex) %>%
    mutate(ancestry = if_else(toupper(ancestry) == 'META', 'META', paste0(toupper(ancestry), '-like'))) %>%
    mutate(ancestry = factor(ancestry, levels = rev(c('AFR-like', 'AMR-like', 'EAS-like', 'EUR-like', 'MID-like', 'SAS-like', 'META')))) %>%
    ggplot + aes(x = ancestry, y = n_phenos, color = category, fill = category) + 
    geom_col(width = 0.5) +
    scale_y_continuous(label = comma, limits=c(0, 4100)) + 
    category_fill_scale +category_color_scale + 
    labs(y = 'Number of phenotypes', x = NULL, color = NULL, fill = NULL) + themes +
    geom_text(aes(label=comma(n_total), y = n_total), color = '#000000',  size = 5, hjust = -0.2) +
    coord_flip() +
    guides(
      fill = guide_legend(nrow = 2, byrow = F, reverse = F,  keywidth = 0.5, keyheight = 0.2), 
      color = guide_legend(nrow =2, byrow = F, reverse = F,  keywidth = 0.5, keyheight = 0.2)
    )+ 
    theme(
      legend.key = element_rect(fill = NA, colour = NA),
      legend.key.spacing.x = unit(0.1, "cm"),
      legend.key.spacing.y = unit(0.1, "cm"),
      legend.text = element_text(margin = margin(r = 0)),
      legend.margin = margin(0, 0, 0, 0),
      legend.box.margin = margin(0, 0, 0, 0),
      axis.text.y = element_text(family = "mono"),
      plot.margin = unit(c(0.5, 0,0,0), 'cm')) + axa_themes
  
  category_colors <- c('#334195', '#d56f3e', '#7e57c2', '#4f345a','#880d1e', '#43aa8b', '#D470A2', '#FFE584', 'gray')
  names(category_colors) <- c('Lab measurements', 'Physical measurements', 'Phecode', 'Disease (PhecodeX)', 'Prescriptions', 'Disease (self-reported)', 'Onset', 'Progression', 'Random phenotype')
  category_labels <- c('Lab measurements', 'Physical measurements', 'Phecode', 'Disease (PhecodeX)', 'Prescriptions', 'Disease (self-reported)', 'Onset', 'Progression', 'Random phenotype')
  category_color_scale<- scale_color_manual(name = NULL, values = category_colors, breaks = names(category_colors), labels = category_labels) 
  category_fill_scale <- scale_fill_manual(name = NULL, values = category_colors, breaks = names(category_colors), labels = category_labels)
  
  df_long <- df %>%
    pivot_longer(
      cols = `6`:`1`,
      names_to = "Num_pop_defined",
      values_to = "value"
    ) %>%
    mutate(Category = factor(Category, levels = rev(c('Physical measurements', 'Lab measurements', 
                                                      'Self-reported', 'PhecodeX', 'Prescriptions', 
                                                      'Phecode', 'Onset', 'Progression', 
                                                      'Random phenotype')),
                             labels = rev(c('Physical measurements', 'Lab measurements', 
                                            'Disease (self-reported)', 'Disease (PhecodeX)', 'Prescriptions', 
                                            'Phecode', 'Onset', 'Progression', 
                                            'Random phenotype')
                               
                             )
                             ))
  

  
  p_middle <- ggplot(df_long, aes(x = Num_pop_defined, y = Category)) +
    geom_tile(color = 'white', aes(fill = Category), alpha = 0.5) +
    # Use any color scale you prefer; below is a simple white-to-blue gradient
    geom_text(aes(label = comma(value)), color = "black", na.rm = TRUE, family = "mono", size = 5) +
    labs(
      # title = paste0('Number of Ancestry Groups Defined'),
      x = 'Number of similarity groups',
      y = NULL,
    ) +
    theme_minimal() + category_color_scale + category_fill_scale +
    # Optionally, ensure the order of columns is 6 -> 5 -> 4 -> 3 -> 2 -> 1
    scale_x_discrete(limits = c("6", "5", "4", "3", "2", "1")) + themes+
    theme(
      panel.grid = element_blank(),
      legend.position = 'none',
      axis.text.x = element_text(face = 'bold'),
      axis.text = element_text(color = 'black'),
      plot.margin = unit(c(0.5, 0,0,0), 'cm')
    )  + axa_themes
  
  p_right <- n_qced_phenotypes %>%
    mutate(category = factor(category, levels = rev(c('Physical measurements', 'Lab measurements', 
                                                      'Self-reported', 'PhecodeX', 'Prescriptions', 
                                                      'Phecode', 'Onset', 'Progression', 
                                                      'Random phenotype')),
                             labels = rev(c('Physical measurements', 'Lab measurements', 
                                            'Disease (self-reported)', 'Disease (PhecodeX)', 'Prescriptions', 
                                            'Phecode', 'Onset', 'Progression', 
                                            'Random phenotype')
                                          
                             )
    ))%>%
    ggplot + aes(x = category, y = n, color = category, fill = category) +
    geom_col(width = 0.5, size = 0)+
    scale_y_continuous(labels = comma, limits=c(0, 3000)) +
    category_color_scale + category_fill_scale +
    themes + labs(x= NULL, y=NULL) + guides(colour = 'none', fill= 'none') + 
    coord_flip() + theme_classic() + axa_themes + 
    theme(panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(), 
          axis.line.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(), 
          axis.line.x = element_blank(),
          legend.position = c(0.9, 0.9)) +
    guides(
      alpha = guide_legend(nrow = 2, byrow = TRUE, reverse = F,  keywidth = 1, keyheight = 0.3)
    ) + geom_text(aes(y = n, label = comma(n)),  size = 4, hjust = -0.1, show_guide = F)
  
  pp = ggpubr::ggarrange(p_left, 
                         p_middle,  
                         p_right + theme(plot.margin = unit(c(0,0.3, 0,0), 'cm')), 
                         nrow=1, hjust = 0, 
                          align = 'h',widths = c(1.8, 1.9, 0.7), 
                         labels = c('(A) Number of phenotypes analyzed by similarity groups', 
                                    '(B) Number of phenotypes overlapped across number of groups',
                                    ''),
                         font.label = list(size = 20, color = "black", face = "bold", family = 'sans'), 
                         common.legend = T, legend = 'bottom'
  )
  if(save){
    pp_long_bottom = ggpubr::ggarrange(p_middle,
                                       p_right + theme(legend.position = 'none'), 
                                       nrow=1, hjust = 0, 
                                       align = 'h',widths = c(2.5, 0.7),
                                       font.label = list(size = 10, color = "black", face = "bold", family = 'sans')
    )
    pp_long = ggpubr::ggarrange(p_left, 
                                pp_long_bottom + theme(plot.margin = unit(c(0, 0.3, 0,0), 'cm')), 
                                labels = c('a', 'b'),
                           ncol=1, hjust = 0, 
                           heights = c(1, 1), 
                           font.label = list(size = 20, color = "black", face = "bold", family = 'sans'), 
                           common.legend = T, legend = 'top'
    )
    output_figure(pp_long, 'summary', paste0('aou_qced_pheno_summary_v8_total_', sex), height=6, width=7.5)
  }
  return(pp_long)
}

p_pheno <- plot_phenotype_by_category('both', df_both, save = T)

pheno_corr <- read_gcs_fast('gs://aou_wlu/v8_analysis/pheno_corr.txt.bgz')
# pheno_corr <- pheno_corr %>%
#   merge(., full_pheno_sum %>% select(phenoname, category, description) %>% distinct(), by.x = 'i_pheno', by.y = 'phenoname') %>%
#   merge(., full_pheno_sum %>% select(phenoname, category, description) %>% distinct(), by.x = 'j_pheno', by.y = 'phenoname', suffixes = c('.i', '.j')) %>%
#   filter(`category.i` %in% c('pfhh_survey', 'mcc2_phecodex') & `category.j` %in% c('pfhh_survey', 'mcc2_phecodex') )


pfhh_phecodex <- read_gcs_fast('gs://aou_wlu/v8_analysis/phenotype/pfhh_phecodex_pairwise_metrics_gene.txt.bgz') 
sub_pheno_corr <- pheno_corr %>%
  merge(., pfhh_phecodex %>% mutate(pfhh_code = as.character(pfhh_code)), by.x = c('i_pheno', 'j_pheno'), by.y = c('phecodex', 'pfhh_code'))

pfhh_phecodex_annotated <- pfhh_phecodex %>%
  mutate(pfhh_code = as.character(pfhh_code)) %>%
  merge(., pheno_info_full %>% select(phenoname, description) %>% distinct(), by.x = 'pfhh_code', by.y = 'phenoname', all.x = T) %>%
  merge(., pheno_info_full %>% select(phenoname, description, disease_category) %>% distinct(), by.x = 'phecodex', by.y = 'phenoname', all.x = T, suffix = c('.pfhh', '.phecodex')) %>%
  merge(., pheno_corr, by.y = c('i_pheno', 'j_pheno'), by.x = c('phecodex', 'pfhh_code'))

disease_names <- unique(pfhh_phecodex_annotated$disease_category)
icd_colors <- icd_colors[-1]
names(icd_colors) <- disease_names
icd_colors <- icd_colors[which(!is.na(names(icd_colors)))]
disease_color_scale<- scale_color_manual(name = NULL, values = icd_colors, breaks = names(icd_colors), labels = names(icd_colors)) 
disease_fill_scale <- scale_fill_manual(name = NULL, values = icd_colors, breaks = names(icd_colors), labels = names(icd_colors))

p_pfhh <- pfhh_phecodex_annotated %>%
  mutate(corr_bin = cut(entry, breaks = c(0, 0.15, 0.3, Inf),
                        include.lowest = TRUE, right = FALSE,
                        labels = c('<0.15', '0.15-0.3', '>0.3'))) %>%
  ggplot + aes(x = pfhh_prevalence, y = phecodex_prevalence, color = corr_bin) +
  geom_point(size = 3) +
  geom_abline(slope = 1, intercept = 0, lty = 2, lwd = 0.1) + themes +
  geom_text_repel(data = . %>% filter((pfhh_prevalence/ phecodex_prevalence > 1.2 | pfhh_prevalence/ phecodex_prevalence < 0.8) & (phecodex_prevalence > 0.2 | pfhh_prevalence > 0.2)),
                  aes(label = description.phecodex), size = 5, nudge_x = 0.01, nudge_y = 0,
                  show.legend = FALSE, color = 'grey30', max.overlaps = 20) +
  scale_x_continuous(label = percent) +
  scale_y_continuous(label = percent) +
  scale_color_brewer(name = 'Phenotypic Correlation', palette = 'Reds', direction = 1) +
  labs(x = 'Self-reported Prevalence', y = 'PhecodeX Prevalence') +
  theme(legend.position = 'top') + axa_themes
p_pfhh
output_figure(p_pfhh, folder = 'phenotype', name = 'pfhh_phecode_prevalence', 3, 4)
manhattan_path <- '~/Dropbox (Partners HealthCare)/aou/axaou/v8/'
p_obesity <- ggdraw() +
  draw_image(image_read(paste0(manhattan_path, 'ACAF/META/phenotype_EM_236.1.META.Pvalue_log10.variant.manhattan.png')))+
  theme(plot.margin = unit(c(0.5, 0, 0, 0), "cm"))
p_hypertension <- ggdraw() +
  draw_image(image_read(paste0(manhattan_path, 'ACAF/META/phenotype_CV_401.1.META.Pvalue_log10.variant.manhattan.png')))+
  theme(plot.margin = unit(c(0.5, 0, 0, 0), "cm"))

pp = ggpubr::ggarrange(p_pheno, 
                       p_obesity,
                       p_hypertension, 
                       p_pfhh, 
                       ncol=1, hjust = 0, heights =c(1.5, 0.6, 0.6, 1), labels = c('', 'c', 'd', 'e'),
                       font.label = list(size = 20, color = "black", face = "bold", family = 'Arial')
)
output_figure(pp, 'main', 'figure1', height=400/25.4, width=89*2/25.4)



################### Archived ##############################
# p1 <- plot_aou_sample_distribution(ver='v8', raw=F)
p1 <- aou_data %>%
  filter(version == 'v8') %>%
  filter(ancestry != 'meta') %>%
  select(ancestry, version, n_samples, n_samples_raw) %>%
  mutate(ancestry = factor(ancestry, levels = rev(names(pops)))) %>%
  mutate(n_samples_diff = n_samples_raw - n_samples,
         label = paste0(comma(n_samples_raw), if_else(ancestry =='eur', ' → ', ' → '),  comma(n_samples))) %>% 
  select(-n_samples_raw) %>%
  pivot_longer(cols = starts_with('n_samples'), names_to = 'qc_status', values_to = 'number') %>%
  group_by(ancestry) %>%
  mutate(N = sum(number)) %>%
  mutate(qc_status = factor(if_else(qc_status == 'n_samples', 'QCed', 'Raw'), levels= c('Raw', 'QCed'))) %>%
  ggplot + aes(x = ancestry, y = number, color = ancestry, fill = ancestry, alpha = qc_status) +
  geom_col(width = 0.5, size = 0)+
  scale_y_continuous(label = comma, limits=c(0, 330000), breaks = c(0, 50000, 150000, 250000)) + 
  # scale_y_continuous(labels = scientific_10, limits=c(0, 285000)) +
  scale_alpha_discrete(name = NULL, range = c(0.2, 1)) + 
  scale_x_discrete(breaks = names(pops), labels = str_pad(toupper(names(pops)),side = 'left', width = 4)) + 
  scale_color_manual(name = NULL, values = pop_colors, breaks = names(pops), labels = toupper(names(pops)))+
  scale_fill_manual(name = NULL, values = pop_colors, breaks = names(pops), labels = toupper(names(pops))) + 
  themes + labs(x= NULL, y=NULL) + coord_flip() + theme_classic() + 
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.position = c(0.9, 0.9), 
        legend.text = element_text(size=6),
        axis.text.y = element_text(size = 10, family = "mono", color='black'),
        axis.text.x = element_text(size = 8, hjust = 0.25)) +
  guides(
    colour = 'none', fill= 'none',
    alpha = guide_legend(nrow = 2, byrow = TRUE, reverse = F,  keywidth = 1, keyheight = 0.3)
  ) + geom_text(aes(y = N, label = label),  size = 3, hjust = -0.1, show_guide = F)

p2 <- aou_data %>%
  filter(version == 'v8') %>%
  filter(ancestry != 'meta') %>%
  mutate(ancestry = factor(ancestry, levels = rev(names(pops)))) %>%
  select(ancestry, version, n_snps, n_indels, number_variants = n_variants) %>%
  pivot_longer(cols = starts_with('n_'), names_to = 'variant_type', values_to = 'number') %>%
  mutate(variant_type = if_else(variant_type == 'n_snps', 'SNPs', 'Indels')) %>%
  ggplot + aes(x = ancestry, y = number, color = ancestry, fill = ancestry, alpha = variant_type) +
  geom_col(width = 0.5, size = 0)+
  scale_y_continuous(labels = scientific_10, limits=c(0, 980000000)) +
  scale_alpha_discrete(name = NULL, range = c(0.2, 1)) + 
  scale_x_discrete(breaks = names(pops), labels = pops) + 
  scale_color_manual(name = 'Ancestry', values = pop_colors, breaks = names(pops), labels = toupper(names(pops)))+
  scale_fill_manual(name = 'Ancestry', values = pop_colors, breaks = names(pops), labels = toupper(names(pops))) + 
  themes + labs(x= NULL, y=NULL) + guides(colour = 'none', fill= 'none') + coord_flip() + theme_classic() + 
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.position = c(0.9, 0.9), 
        legend.text = element_text(size=6),
        axis.text.x = element_text(size = 8, hjust = 0.25)) +
  guides(
    alpha = guide_legend(nrow = 2, byrow = TRUE, reverse = F,  keywidth = 1, keyheight = 0.3)
  ) + geom_text(aes(y = number_variants, label =comma(number_variants)),  size = 3, hjust = -0.1, show_guide = F)

# p3 <- plot_aou_sample_distribution(ver='v8', raw=F, type = 'phenotypes')

p_top = ggpubr::ggarrange(p1 + theme(plot.margin = unit(c(0.5, 0, 0, 0), 'cm')), 
                          p2+ theme(axis.ticks.y = element_blank(), axis.text.y = element_blank(), plot.margin = unit(c(0.5, 0.5, 0, 0), 'cm')), 
                          # p3 + theme(axis.ticks.y = element_blank(), axis.text.y = element_blank(), plot.margin = unit(c(0.5, 0.5, 0, 0), 'cm')), 
                          nrow=1, hjust = 0, widths = c(1,1),align = "h", 
                          labels = c('(A) Number of participants', '(B) Number of variants'), 
                          font.label = list(size = 10, color = "black", face = "bold", family = 'Arial')
)
output_figure(p_top, 'main', 'figure1_top', height=2, width=10)

pp = ggpubr::ggarrange(p_top,
                       p_middle, 
                       ncol=1, heights =c(0.9, 1),
                       font.label = list(size = 110, color = "black", face = "bold", family = 'Arial')
)

output_figure(pp, 'main', 'figure1_middle', height=4, width=10)

figure1_bottom_left = ggpubr::ggarrange(p_bottom1 + 
                                          theme(plot.margin = unit(c(0.5, 0, 0, 0), 'cm'),
                                                axis.title = element_text(face = 'plain', size = 8)), 
                                        p_bottom2 + 
                                          theme(plot.margin = unit(c(0.5, 0, 0, 0.2), 'cm'),
                                                axis.title = element_text(face = 'plain', size = 8)), 
                                        labels = c('(E) Self-reported & PhecodeX prevalence comparison',
                                                   '(F) Self-reported p-value comparison'),
                                        nrow =1, hjust=0, vjust = 1, common.legend = T, align = 'h', legend = 'bottom', 
                                        font.label = list(size = 10, color = "black", face = "bold", family = 'Arial'))
output_figure(figure1_bottom_left, 'main', 'figure1_bottom_left', height =3, width = 6)

figure1_bottom = ggpubr::ggarrange(figure1_bottom_left, p_bottom3  + theme(plot.margin = unit(c(0.5, 0, 0, 0.2), 'cm'),
                                                                    axis.title = element_text(face = 'plain', size = 8)), 
                                   labels = c('',
                                              '(G) Self-reported & PhecodeX association comparison'),
                                   nrow =1, hjust=0, vjust =1, widths = c(2, 1),
                                   font.label = list(size = 10, color = "black", face = "bold", family = 'Arial'))
output_figure(figure1_bottom, 'main', 'figure1_bottom', height =3, width = 8)


pp = ggpubr::ggarrange(p_top, p_middle, 
                       figure1_bottom + theme(plot.margin = unit(c(0.3, 0.5, 0, 0.3), 'cm')), 
                       ncol=1, hjust = 0, heights =c(1, 1.2, 1.2), labels = c('', '',''),
                       font.label = list(size = 10, color = "black", face = "bold", family = 'Arial')
)
output_figure(pp, 'main', 'figure1', height=8, width=10)
