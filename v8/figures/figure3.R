source('~/Dropbox (Partners HealthCare)/github_repo/ukbb_exomes_pleiotropy/R/constants.R')
library(magick)
manhattan_path <- '~/Dropbox (Partners HealthCare)/aou/axaou/v8/'
threshold = 6.7e-7

axa_themes <- theme(
  axis.title = element_text(face = 'plain', size = 15),
  axis.text =  element_text(size = 14),
  legend.text = element_text(size = 13),
)

p_GIGYF1 <- ggdraw() +
  draw_image(image_read(paste0(manhattan_path, 'gene/META/phenotype_3004501.META.Pvalue_Burden_log10_0.001.gene.manhattan_pLoF.png')))+
  theme(plot.margin = unit(c(0.5, 0, 0, 0), "cm"))
p_TIMD4 <- ggdraw() +
  draw_image(image_read(paste0(manhattan_path, 'gene/META/phenotype_3022192.META.Pvalue_Burden_log10_0.001.gene.manhattan_pLoF.png')))+
  theme(plot.margin = unit(c(0.5, 0, 0, 0), "cm"))
p_TNN <- ggdraw() +
  draw_image(image_read(paste0(manhattan_path, 'gene/META/phenotype_heart-rate-mean.META.Pvalue_Burden_log10_0.001.gene.qq_and_manhattan_pLoF.png')))+
  theme(plot.margin = unit(c(0.5, 0, 0, 0), "cm"))

# source('~/Dropbox (Partners HealthCare)/github_repo/aou_gwas/R/phecode_mega.R')
# phecode <- read_gcs_fast(paste0('gs://aou_wlu/v8_analysis/phecode_meta_gene_burden.txt.bgz')) %>% 
#   filter(max_MAF == 0.001) 
# p_mega <- phecode_mega_figure(data=phecode, save_plot = F, output_path=paste0('~/Desktop/',pop, '_',  test, '_', maxmaf,'_mega_phecode.png')) + 
#   theme(axis.title = element_text(size = 8, face = 'bold'),
#         axis.text = element_text(size = 6),
#         plot.margin = unit(c(0.5,0,0,0), 'cm'))

meta_detail <- read_gcs_fast(paste0('gs://aou_wlu/v8_analysis/paper_preparation/summary/aou_v8_META_gene_signals_Pvalue_Burden_0.001_qced.txt.bgz')) %>%
  filter(!startsWith(phenoname, 'random') & !(phenoname %in% c('PP_928.1', 'PP_932.11'))) %>%
  filter(Pvalue_Burden < 6.7e-7) %>%
  dplyr::mutate(annotation = factor(if_else(annotation == 'missenseLC', 'missense|LC', if_else(annotation == 'pLoF;missenseLC', 'pLoF|missense|LC', annotation)), levels = annotation_types2))

pleiotropic_genes <- meta_detail %>%
  select(gene_id, gene_symbol, annotation, category, BETA_Burden) %>%
  group_by(gene_id, gene_symbol, annotation) %>%
  dplyr::summarize(n_domain = dplyr::n_distinct(category),
                   n_pos = sum(BETA_Burden>0),
                   n_neg = sum(BETA_Burden < 0) )%>%
  filter(n_domain > 1 & annotation != 'pLoF|missense|LC')

meta_detail_filtered <- meta_detail %>%
  group_by(gene_id, gene_symbol, annotation) %>%
  dplyr::mutate(n_domain = dplyr::n_distinct(category)) %>%
  filter(n_domain > 1 & annotation != 'pLoF|missense|LC')

pfhh_phecodex <- read_gcs_fast('gs://aou_wlu/v8_analysis/phenotype/pfhh_phecodex_pairwise_metrics_gene.txt.bgz') 

pfhh_phecodex_annotated <- pfhh_phecodex %>%
  mutate(pfhh_code = as.character(pfhh_code)) %>%
  merge(., pheno_info_full %>% select(phenoname, description) %>% distinct(), by.x = 'pfhh_code', by.y = 'phenoname', all.x = T) %>%
  merge(., pheno_info_full %>% select(phenoname, description, disease_category) %>% distinct(), by.x = 'phecodex', by.y = 'phenoname', all.x = T, suffix = c('.pfhh', '.phecodex'))

disease_names <- unique(pfhh_phecodex_annotated$disease_category)
icd_colors <- icd_colors[-1]
names(icd_colors) <- disease_names
icd_colors <- icd_colors[which(!is.na(names(icd_colors)))]
disease_color_scale<- scale_color_manual(name = NULL, values = icd_colors, breaks = names(icd_colors), labels = names(icd_colors)) 
disease_fill_scale <- scale_fill_manual(name = NULL, values = icd_colors, breaks = names(icd_colors), labels = names(icd_colors))

pfhh_phecodex_p <- read_gcs_fast('gs://aou_wlu/v8_analysis/phenotype/pfhh_phecodex_signals.txt.bgz') 
pfhh_phecodex_p <- pfhh_phecodex_p %>%
  mutate(phenoname_1 = as.character(phenoname_1)) %>%
  merge(., pheno_info_full %>% select(phenoname, description) %>% distinct(), by.x = 'phenoname_1', by.y = 'phenoname', all.x = T) %>%
  merge(., pheno_info_full %>% select(phenoname, description, disease_category) %>% distinct(), by.x = 'phenoname', by.y = 'phenoname', all.x = T, suffix = c('.pfhh', '.phecodex'))
# pfhh_phecodex_p_for_josh <- pfhh_phecodex_p %>%
#   mutate(
#     PhecodeX = phenoname,
#     PFHH = phenoname_1,
#     META_Pvalue_Burden.phecodex = META_Pvalue_Burden,
#     META_Pvalue_Burden.pfhh = META_Pvalue_Burden_1
#   ) %>%
#   select(gene_id, gene_symbol, annotation, max_MAF, PhecodeX, description.phecodex, META_Pvalue_Burden.phecodex, PFHH, description.pfhh, META_Pvalue_Burden.pfhh, disease_category)
# write_tsv(pfhh_phecodex_p_for_josh, '~/Desktop/pfhh_phecodex_disease_category_comparison.tsv')

pfhh1 <- pfhh_phecodex_p %>%
  dplyr::filter(max_MAF == 0.001) %>%
  filter(META_Pvalue_Burden_1 < 1e-4 | META_Pvalue_Burden < 1e-4) %>%
  ggplot + aes(x = -log10(META_Pvalue_Burden_1), y = -log10(META_Pvalue_Burden), color = disease_category) +
  geom_point(size = 2) + 
  geom_hline(yintercept = -log10(threshold), lty = 2, lwd = 0.1) + 
  geom_vline(xintercept = -log10(threshold), lty = 2, lwd = 0.1) + 
  geom_abline(slope = 1, intercept = 0, lty = 2, lwd = 0.1) + themes +
  scale_x_continuous(trans = locusviz::trans_loglog_p(), labels=comma, breaks = c(0, 5, 10, 30, 100, 300, 500)) +
  scale_y_continuous(trans = locusviz::trans_loglog_p(), labels=comma, breaks = c(0, 5, 10, 30, 100, 300, 500)) +
  geom_text_repel(data = pfhh_phecodex_p %>%
                    filter(max_MAF == 0.001)%>%
                    filter(META_Pvalue_Burden_1 < 1e-10 | META_Pvalue_Burden < 1e-10),
                  aes(label = paste0(description.phecodex, '-', gene_symbol, '-', annotation)), 
                  size = 3,nudge_x  = 0.01, nudge_y = 0, show.legend = FALSE) + 
  # scale_color_brewer(palette = 'Spectral') + 
  disease_color_scale + disease_fill_scale +
  labs(x = expression(Self-reported~-log[10](p)), y = expression(PhecodeX~-log[10](p))) +
  guides(color = guide_legend(ncol = 4, byrow = TRUE, keywidth = unit(0.2, "cm"), override.aes = list(size = 2)), size = 'none')+
  theme(plot.margin = unit(c(0.5, 0, 0, 0), 'cm'),
        legend.key = element_rect(fill = NA, colour = NA),
        legend.key.spacing.x = unit(0.5, "cm"),
        legend.key.spacing.y = unit(-0.1, "cm"),
        legend.text = element_text(margin = margin(r = 2)),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(0, 0, 0, 0)) + axa_themes
output_figure(pfhh1, folder = 'phenotype', name = 'pfhh_phecode_pvalue', 3, 4)


pfhh_phecodex_p <- pfhh_phecodex_p %>%
  mutate(pfhh_yes_phecodex_no = META_Pvalue_Burden_1 < threshold & META_Pvalue_Burden > threshold,
         pfhh_no_phecodex_yes = META_Pvalue_Burden_1 > threshold & META_Pvalue_Burden < threshold) %>%
  mutate(disease_category = factor(disease_category, levels = names(icd_colors)))
sum(pfhh_phecodex_p$META_Pvalue_Burden_1 < threshold & pfhh_phecodex_p$META_Pvalue_Burden > threshold, na.rm = T)
sum(pfhh_phecodex_p$META_Pvalue_Burden_1 > threshold & pfhh_phecodex_p$META_Pvalue_Burden < threshold, na.rm = T)
table(pfhh_phecodex_p %>% filter(pfhh_yes_phecodex_no) %>% select(description.phecodex, disease_category) %>% distinct() %$% disease_category)
table(pfhh_phecodex_p %>% filter(pfhh_no_phecodex_yes)  %>% select(description.phecodex, disease_category) %>% distinct()  %$% disease_category)
table(pfhh_phecodex_p %>% filter(pfhh_yes_phecodex_no) %$% disease_category)
table(pfhh_phecodex_p %>% filter(pfhh_no_phecodex_yes)  %$% disease_category)

pfhh_phecodex_p_count <- pfhh_phecodex_p %>% 
  filter(pfhh_no_phecodex_yes)%>%
  dplyr::group_by(disease_category) %>%
  dplyr::summarize(n_assoc = n()) %>%
  mutate(label = 'PhecodeX only (N = 91)') %>%
  rbind(., pfhh_phecodex_p %>% 
          filter(pfhh_yes_phecodex_no)%>%
          group_by(disease_category) %>%
          dplyr::summarize(n_assoc = n())%>%
          mutate(label = 'PFHH only (N = 93)')) %>%
  rbind(., data.frame(
    disease_category = c('Dermatological', 'Genitourinary', 'Respiratory'),
    n_assoc = c(0,0,0), label = c('PhecodeX only (N = 91)', 'PhecodeX only (N = 91)', 'PFHH only (N = 93)')
  ))
pfhh_phecodex_p_count %>%
  group_by(label) %>%
  dplyr::summarize(n = sum(n_assoc))

pfhh2 <- pfhh_phecodex_p_count %>%
  mutate(disease_category = factor(disease_category, levels = rev(names(icd_colors))))%>%
  # dplyr::filter(max_MAF == 0.001) %>%
  ggplot + aes(x = disease_category, y = n_assoc, color = disease_category, 
               fill = disease_category, alpha = label, group = interaction(label, disease_category)) + 
  geom_bar(stat = "identity", position = position_dodge(), size = 0.2, width = 0.5)  +
  scale_alpha_discrete(name = NULL, range = c(0.4, 1)) + 
  disease_color_scale + disease_fill_scale +
  labs(x = NULL, y = 'Number of associations')  + themes +
  guides(
    alpha = guide_legend(nrow = 2, byrow = TRUE, reverse = FALSE,  keywidth = 0.5, keyheight = 0.2),
    color = 'none', fill = 'none'
  )+
  theme(legend.position = c(0.7, 0.7),
        legend.box.spacing = unit(0, "cm"),
        legend.spacing.x = unit(0, "cm"),
        legend.spacing.y = unit(0, "cm")
  ) + 
  coord_flip() + axa_themes

output_figure(pfhh2, folder = 'phenotype', name = 'pfhh_phecode_association_comparison', 3, 4)



############################################
# Opposite phenotype correlation analysis: #
############################################
pheno_list <- c(meta_detail_filtered  %>%
                  select(phenoname) %>%
                  distinct() %$%
                  phenoname)

pheno_corr <- read_gcs_fast('gs://aou_wlu/v8_analysis/pheno_corr.txt.bgz')
sub_pheno_corr <- pheno_corr %>%
  filter(i_pheno %in% pheno_list & j_pheno %in% pheno_list)

result_list <- list()
for(gene in c(pleiotropic_genes %>% filter(n_pos > 0 & n_neg > 0) %$% gene_symbol)){
  print(paste0('--------------------', gene, '---------------'))
  sub_meta_detail <- meta_detail_filtered %>% filter(gene_symbol == gene) %>% select(gene_symbol, phenoname, description, annotation, BETA_Burden)
  pheno_list <- c(meta_detail %>% filter(gene_symbol == gene) %$% phenoname)
  tmp <- sub_pheno_corr %>%
    filter(i_pheno %in% pheno_list & j_pheno %in% pheno_list)%>%
    merge(., sub_meta_detail, by.x = 'i_pheno', by.y = 'phenoname', all.x = T)  %>%
    merge(., sub_meta_detail, by.x = c('gene_symbol',  'annotation', 'j_pheno'), by.y = c('gene_symbol',  'annotation', 'phenoname'), all.x = T, suffixes=c('.i', '.j')) %>%
    mutate(corr_sign = sign(entry) > 0, beta_sign = sign(BETA_Burden.i) == sign(BETA_Burden.j))
  filtered <- tmp %>% ungroup() %>% filter(corr_sign != beta_sign & entry^2 > 0.005) %>% 
    select(gene_symbol, annotation, description.i, description.j, corr=entry, BETA_Burden.i, BETA_Burden.j)
  result_list[[gene]] <- filtered
  print(filtered)
}
result_df <- bind_rows(result_list)
write_tsv(result_df, '~/Dropbox (Partners HealthCare)/aou_v8/tables/supplementary_data6.tsv')

##########################
# Pleiotropy upset plot: #
##########################
format_entries <- function(x) {
  x <- as.numeric(x)
  x[is.na(x)] <- 0
  as.integer(x > 0)
}

categories <- c('Lab measurements', 'Physical measurements', 'Phecode', 'PhecodeX', 'Prescriptions', 'Self-reported')
names(categories) <- c('lab_measurement', 'physical_measurement', 'mcc2_phecode', 'mcc2_phecodex', 'r_drug', 'self-reported')

gene_upset <- meta_detail_filtered %>%
  mutate(category = if_else(category %in% c('pfhh_survey', 'mhwb_survey'), 'self-reported', category))%>%
  mutate(category = factor(category, levels = names(categories), labels = categories)) %>%
  group_by(gene_id, gene_symbol, annotation, max_MAF, category) %>%
  dplyr::summarize(
    n_sig_per_domain = n() 
  ) %>%
  mutate(presence = if_else(n_sig_per_domain > 0, 1, 0)) %>%
  pivot_wider(., id_cols = 1:4, names_from = 'category', values_from = 'presence') %>%
  ungroup(.) %>%
  mutate(across(5:last_col(), format_entries))  %>%
  dplyr::mutate(n_domains = `Lab measurements` +
                  `Physical measurements` + PhecodeX + Prescriptions+
                  `Self-reported`) %>%
  dplyr::filter(n_domains > 1) %>%
  dplyr::select(-n_domains)


gene_upset_detail <- gene_upset %>%
  dplyr::select(1:4) %>%
  merge(., meta_detail_filtered, by = c('gene_id', 'gene_symbol', 'annotation', 'max_MAF'), all.x = T)

write_tsv(gene_upset_detail, '~/Dropbox (Partners HealthCare)/aou_v8/tables/supplementary_data5.tsv')
length(table(gene_upset_detail %>% mutate(label = paste0(gene_id, '-', annotation)) %$% label))

p_pleiotropy <- upset(
  gene_upset , 
  intersect = categories[c(1:2,4:6)],
  name = '',
  # queries = list(
  #   upset_query(set = categories['lab_measurement'], fill = category_colors[categories['lab_measurement']]),
  #   upset_query(set = categories['physical_measurement'], fill = category_colors[categories['physical_measurement']]),
  #   # upset_query(set = categories['mhwb_survey'], fill = category_colors[categories['mhwb_survey']]),
  #   upset_query(set = categories['pfhh_survey'], fill = category_colors[categories['pfhh_survey']]),
  #   upset_query(set = categories['mcc2_phecodex'], fill = category_colors[categories['mcc2_phecodex']]),
  #   upset_query(set = categories['r_drug'], fill = category_colors[categories['r_drug']])
  # ),
  queries = list(
    upset_query(set = categories['lab_measurement'], fill = 'white'),
    upset_query(set = categories['physical_measurement'], fill = 'white'),
    upset_query(set = categories['mhwb_survey'], fill = 'white'),
    upset_query(set = categories['pfhh_survey'], fill = 'white'),
    upset_query(set = categories['mcc2_phecodex'], fill = 'white'),
    upset_query(set = categories['r_drug'], fill = 'white')
  ),
  base_annotations=list(
    ' '=(
      intersection_size(
        text=list(
          vjust=-0.1,
          size = 5
        ),
        mapping=aes(fill=annotation),
        bar_number_threshold=1,  # show all numbers on top of bars
        width=0.5,   # reduce width of the bars
      )
      # add some space on the top of the bars
      + scale_y_continuous(expand=expansion(mult=c(0, 0.05)), limits = c(0, 16))
      + annotation_fill_scale
      +labs(x=NULL, y=NULL, color = NULL, fill=NULL)
      + theme(
        # hide grid lines
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        # show axis lines
        axis.line.x=element_line(colour='black', linewidth = 0.1),
        axis.line.y=element_blank(),
        axis.text.y=element_blank(),
        legend.text = element_text(size=11, color = 'black'),
        legend.position = c(0.9, 0.72),
        legend.direction = 'vertical',
        legend.title = element_blank()
      )+
        guides(
          fill = guide_legend(keywidth = 0.8, keyheight = 0.5)
        )
    )
  ),
  stripes=upset_stripes(
    geom=geom_segment(size=5),  # make the stripes larger
    colors=c('grey95', 'white')
  ),
  matrix=intersection_matrix(
    geom=geom_point(
      shape='circle filled',
      size=1.5,
      stroke=0.45
    ),
  )+ labs(x=NULL, y=NULL),
  set_sizes = FALSE,
  themes=upset_modify_themes(
    list(
      'intersections_matrix'=theme(text=element_text(size=15, color='black')),
      'overall_sizes'=theme(axis.text.x=element_text(size=15, vjust = 1))
    )
  )
) + labs(x='Phenotypic domain combinations', y = NULL) + axa_themes

output_figure(p_pleiotropy, 'main', 'figure3_pleiotropy', height =4, width =7.5)

p_pfhh_phecodex = ggpubr::ggarrange(
  pfhh1 , 
  pfhh2, 
  labels = c('c', 'd'),
  ncol =1, hjust=0, vjust = 1, common.legend = T, legend = 'top', 
  font.label = list(size = 20, color = "black", face = "bold", family = 'Arial'))

figure3 = ggpubr::ggarrange(p_GIGYF1, 
                            p_TIMD4, 
                            p_pfhh_phecodex + theme(plot.margin = unit(c(0.2,0,0,0), 'cm')) , 
                            p_pleiotropy, 
                            labels = c('a', 'b', '', 'e'), hjust=0, vjust = 1.2, heights = c(1.1, 1.1, 3, 1.3),
                            font.label = list(size = 20, color = "black", face = "bold", family = 'Arial'), 
                            ncol=1)
output_figure(figure3, 'main', 'figure3', height =370/25.4, width = 89*2/25.4)
