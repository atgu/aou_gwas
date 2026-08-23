source('~/Dropbox (Partners HealthCare)/github_repo/aou_gwas/R/constants.R')
# install.packages("BiocManager")
# BiocManager::install("biovizBase")
# BiocManager::install("plyranges")
# BiocManager::install("ensembldb")
# BiocManager::install("AnnotationDbi")
# BiocManager::install("ComplexHeatmap")
# 
# devtools::install_github("mkanai/locusviz", dependencies = c('biovizBase', 'plyranges', 'ensembldb', 'AnnotationDbi'))
library(locusviz)
library(magick)

any_burden <- read_gcs_fast('gs://aou_wlu/v8_analysis/aou_ukb_meta/significant_in_any_burden_0.001.txt.bgz') %>%
  filter(annotation != 'pLoF;missenseLC')

aou_caf <- read_gcs_fast('gs://aou_analysis/v8/data/utils/CAF_pruned/ALL_CAF.txt.bgz')
ukb_caf <- read_gcs_fast('gs://aou_analysis/v8/data/utils/CAF_pruned/genebass_gene_caf_500k.txt.bgz')  %>%
  select(gene_symbol=gene_id, annotation, CAF_genebass_01 = CAF) %>%
  mutate(annotation = if_else(annotation == 'missense|LC', 'missenseLC', annotation))

full_caf <- aou_caf %>%
  merge(., ukb_caf, by = c('gene_symbol', 'annotation')) %>%
  mutate(max_CAF = pmax(CAF_0.01, CAF_genebass_01))

fmt_dcimals <- function(decimals=0){
  # return a function responpsible for formatting the 
  # axis labels with a given number of decimals 
  function(x) as.character(round(x,decimals))
}
high_col = '#008080'
low_col = 'black'

threshold = 1e-10
log_threshold = -log10(threshold)
p1 <- any_burden %>%
  filter(ukb_Pvalue_Burden < threshold | META_Pvalue_Burden < threshold) %>%
  merge(., full_caf, by = c('gene_symbol', 'annotation', 'gene_id') ) %>%
  mutate(annotation = if_else(annotation == 'missenseLC', 'missense|LC', annotation)) %>%
  mutate(annotation = factor(annotation, levels = annotation_types, labels = annotation_names)) %>%
  mutate(col_field = log10(CAF_0.01/CAF_genebass_01)) %>%
  ggplot + aes(x = -log10(ukb_Pvalue_Burden), y = -log10(META_Pvalue_Burden), color = col_field) +
  geom_point() + geom_abline(intercept = 0, slope=1, lty=2, lwd = 0.3) + 
  labs(x = expression(-log[10](Genebass-P[Burden])), y = expression(-log[10](META-P[Burden]))) + 
  scale_x_continuous(trans = locusviz::trans_loglog_p(), labels=comma, breaks = c(0, 5, 10, 30, 100)) +
  scale_y_continuous(trans = locusviz::trans_loglog_p(), labels=comma, breaks = c(0, 5, 10, 30, 100)) +
  scale_color_gradient2(high = high_col, low = low_col, mid = 'gray90',
                        name=expression(paste(log[10], '(', frac(AxA~CAF, Genebass~CAF), ')'))) +
  geom_segment(x = 0, xend = -log10(threshold), y = -log10(threshold), yend = -log10(threshold), color = "black", lty=2, lwd = 0.3) +
  geom_segment(x = -log10(threshold), xend = -log10(threshold), y = 0, yend = -log10(threshold), color = "black", lty=2, lwd = 0.3) +
  themes +
  facet_grid(~annotation, labeller = label_type)  + 
  guides(fill = guide_colorbar(
    barwidth = 0.5,
    barheight = 3,
    title.theme = element_text(size = 3),
    label.theme = element_text(size = 3)
  )) + theme(
    legend.key.size = unit(0.3, "cm"),
    legend.box.spacing = unit(0, "cm"),
    legend.box.margin = margin(0, 0, 0, 0),
  ) + axa_themes

pop_pattern = c('afr' = 'white',
                'amr' = 'grey90',
                'eur' = 'white',
                'aou_meta'= 'grey90',
                'ukb_meta' = 'white',
                'all_meta'= 'grey90')

full_new_hits <- obtain_novelty_data_for_plot()
p2 <- plot_novelty(data= full_new_hits, category_name='mcc2_phecodex', pop_pattern = pop_pattern, indexes = c(1,2,4,7,8,9))

manhattan_path <- '~/Dropbox (Partners HealthCare)/aou/axaou/v8/'
genebass_info <- read_delim('~/Dropbox (Partners HealthCare)/0_very_often_used/genebass_pheno_results.txt.bgz')
p3 <- ggdraw() +
  draw_image(image_read(paste0(manhattan_path, 'gene/aou_ukb_META/META_CV_404_Burden.manhattan.png')))+
  theme(plot.margin = unit(c(0.7, 0, 0, 0), "cm"))
p4 <- ggdraw() +
  draw_image(image_read(paste0(manhattan_path, 'gene/aou_ukb_META/META_DE_660.12_Burden.manhattan.png')))+
  theme(plot.margin = unit(c(0.7, 0, 0, 0), "cm"))


# p_bottom = ggpubr::ggarrange(p3,
#                              p4, 
#                              labels = c('c',
#                                         'd'), 
#                              ncol=2, vjust=1, widths = c(1, 1), hjust = 0,
#                              font.label = list(size = 9, color = "black", face = "bold", family = 'Arial'))
# 
# figure = ggpubr::ggarrange( p1,
#                             p2 + theme(plot.margin = unit(c(1,0,0.3,0), 'cm')) , 
#                             p_bottom, 
#                             labels = c(
#                               'a',
#                               'b', ''), 
#                             ncol=1, vjust=1,hjust = 0, heights = c(1, 1.5, 0.6),
#                             font.label = list(size = 9, color = "black", face = "bold", family = 'Arial'))
# output_figure(figure, 'main','figure4', height =7.5, width = 8)


figure = ggpubr::ggarrange( p1,
                            p2 + theme(plot.margin = unit(c(1,0,0.3,0), 'cm')) , 
                            p3,
                            p4, 
                            labels = c('a', 'b', 'c', 'd'), 
                            ncol=1, vjust=1,hjust = 0, heights = c(0.8, 1, 0.8, 0.8),
                            font.label = list(size = 20, color = "black", face = "bold", family = 'Arial'))
output_figure(figure, 'main','figure4', height =274*1.5/25.4, width = 183*1.5/25.4)




############################ Archived TIMD4 examples ##############################################################
# timd4_new_hits <- rbind(read_gcs_fast('gs://aou_wlu/v8_analysis/aou_ukb_meta/aou_v8_afr_timd4_hits_0.01.txt.bgz'),
#                   read_gcs_fast('gs://aou_wlu/v8_analysis/aou_ukb_meta/aou_v8_amr_timd4_hits_0.01.txt.bgz'),
#                   read_gcs_fast('gs://aou_wlu/v8_analysis/aou_ukb_meta/aou_v8_eas_timd4_hits_0.01.txt.bgz'),
#                   read_gcs_fast('gs://aou_wlu/v8_analysis/aou_ukb_meta/aou_v8_eur_timd4_hits_0.01.txt.bgz'),
#                   read_gcs_fast('gs://aou_wlu/v8_analysis/aou_ukb_meta/aou_v8_mid_timd4_hits_0.01.txt.bgz'),
#                   read_gcs_fast('gs://aou_wlu/v8_analysis/aou_ukb_meta/aou_v8_sas_timd4_hits_0.01.txt.bgz')
# ) %>%
#   select(-max_MAF)
# 
# timd4_hits_info <- read_gcs_fast('gs://aou_wlu/v8_analysis/aou_ukb_meta/timd4_significant_in_meta_only_0.01.txt.bgz')  %>%
#   distinct()
# timd4_formatted_hits_info <- rbind(
#   timd4_hits_info %>% 
#     select(gene_id, gene_symbol, annotation, phenoname, Pvalue_Burden = aou_Pvalue_Burden, BETA_Burden = aou_BETA_Burden) %>%
#     mutate(SE_Burden = 0, ancestry = 'aou_meta'),
#   timd4_hits_info %>% 
#     select(gene_id, gene_symbol, annotation, phenoname, Pvalue_Burden = ukb_Pvalue_Burden, BETA_Burden = ukb_BETA_Burden) %>% 
#     mutate(SE_Burden = 0, ancestry = 'ukb_meta'),
#   timd4_hits_info %>% 
#     select(gene_id, gene_symbol, annotation, phenoname, Pvalue_Burden = META_Pvalue_Burden, BETA_Burden = META_Stats_Burden) %>%
#     mutate(SE_Burden = 0, ancestry = 'all_meta')
# )
# timd4_new_hits <- timd4_new_hits  %>%
#   merge(., timd4_hits_info %>% select(gene_id, gene_symbol, annotation, phenoname) %>% distinct(), by = c('gene_id', 'gene_symbol', 'annotation', 'phenoname'), all.y=T) %>%
#   rbind(., timd4_formatted_hits_info) %>%
#   mutate(pvalue = if_else(Pvalue_Burden < 6.7e-7, 'Significant', 'Non-significant'))
# # Add transformed p-values for heatmap visualization
# timd4_new_hits <- timd4_new_hits %>%
#   merge(., pheno_info %>% select(phenoname, description, category) %>% distinct, by = 'phenoname', all.x = T) %>%
#   mutate(logP = -log10(Pvalue_Burden),
#          assoc = paste0(str_to_sentence(description)))
# 
# pops = c('African/African American (AFR)', 'American Admixed/Latino (AMR)', 'East Asian (EAS)', 'European (EUR)', 'Middle Eastern (MID)', ' South Asian (SAS)', 'All of Us (AoU) Meta-analysis', 'UK Biobank (UKB)', 'AoU UKB Meta-analysis')
# names(pops) <- c('afr', 'amr', 'eas', 'eur', 'mid','sas', 'aou_meta', 'ukb_meta', 'all_meta')
# label_type = labeller(ancestry = pops, annotation=annotation_names)
# color_mde =  '#EEA9B8'
# 
# pop_colors = c('afr' = color_afr,
#                'amr' = color_amr,
#                'eas' = color_eas,
#                'eur' = color_nfe,
#                'mid' = color_mde,
#                'sas' = color_sas,
#                'aou_meta'= '#1E275C',
#                'ukb_meta'= '#007FA3',
#                'all_meta' = color_oth)
# pop_pattern = c('afr' = 'white',
#                 'amr' = 'grey90',
#                 'eas' = 'white',
#                 'eur' = 'grey90',
#                 'mid' = 'white',
#                 'sas'= 'grey90',
#                 'aou_meta'= 'white',
#                 'ukb_meta' = 'grey90',
#                 'all_meta'= 'white')
# 
# if (!require("ggpattern")) install.packages("ggpattern")
# timd4_new_hits2 <- timd4_new_hits %>%
#   mutate(row_color = pop_colors[ancestry],
#          row_pattern = 'circle',
#          row_background = pop_pattern[ancestry], 
#          ancestry = factor(ancestry, levels = names(pops), labels = pops))%>%
#   mutate(assoc = paste0(str_to_sentence(description), ' ~ ', gene_symbol)) %>%
#   merge(., ai_results, by = 'assoc', all.x = T) %>%
#   mutate(verdict = if_else(assoc == 'Triglycerides ~ TIMD4', 'Established', verdict)) %>%
#   mutate(verdict = factor(verdict, levels=c( 'Novel','Hypothesized', 'Existing', 'Established')))
# 
# p4_upper2 <- timd4_new_hits2  %>%
#   mutate(BETA_Burden = as.numeric(BETA_Burden),
#          log_p_scale = case_when(
#            Pvalue_Burden > 0.05 ~ '> 0.05',
#            Pvalue_Burden > 0.0001 ~ '> 0.0001',
#            Pvalue_Burden > 6.7e-7 ~ '> 6.7e-7',
#            Pvalue_Burden <= 6.7e-7 ~ 'Significant',
#          )) %>%
#   mutate(log_p_scale = factor(log_p_scale, levels=rev(c('> 0.05', '> 0.0001', '> 6.7e-7', 'Significant')))) %>%
#   ggplot + aes(x = assoc, y = ancestry, fill)  +
#   # Add row-wise background rectangles
#   geom_rect(
#     aes(
#       xmin = -Inf, xmax = Inf,
#       ymin = as.numeric(as.factor(ancestry)) - 0.5,
#       ymax = as.numeric(as.factor(ancestry)) + 0.5
#     ),
#     fill =  timd4_new_hits2$row_background,
#     # pattern = 'circle', 
#     # pattern_fill = full_new_hits$row_color,       # Pattern color
#     # pattern_color = full_new_hits$row_color, 
#     # pattern_density = 0.1,        # Density of the pattern
#     # pattern_size = 0.5,           # Thickness of the pattern
#     inherit.aes = FALSE, alpha = 0.1
#   ) +
#   geom_point(aes(color = log_p_scale), size = 3) +   # Circle size = p-value, color = beta
#   geom_point(data = timd4_new_hits2,
#              aes(x = assoc, y = 0.5, fill = verdict),
#              shape = 21, size = 3, color = "white") +
#   scale_color_manual(name = 'Burden Pvalue', values = c("#165B69", "#28A4BD", "#96DBE9", "#EAF8FB"), breaks = rev(c('> 0.05', '> 0.0001', '> 6.7e-7', 'Significant'))) + 
#   scale_fill_manual(values = c("#D73027", "#F46D43", "#FDAE61", "#FEF6E8"), breaks = c( 'Not found','Hypothesized', 'Existing', 'Established'), name = NULL) +
#   scale_size_manual(name = NULL, values = c("Significant" = 5, "Non-significant"=3)) + 
#   geom_text(aes(label = if_else(BETA_Burden > 0, '+', '-')), color = "black", size = 4) +
#   labs(x = NULL, y = NULL, color = expression(-log[10]~P)) +
#   themes + 
#   theme(axis.text.y = element_text(colour = pop_colors, size = 9, face = 'bold'),
#         axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
#         axis.ticks.y = element_blank(),
#         axis.line.y = element_blank(),
#         axis.ticks.x = element_line(linewidth = 0.1),
#         axis.line.x = element_line(linewidth = 0.1),
#         legend.position = 'top',
#         legend.title = element_text(size = 8),
#         legend.text = element_text(size = 8)
#   ) + 
#   guides(
#     color = 'none',
#     size = 'none',
#     fill = 'none'
#   )
# p4_upper2
# output_figure(p4_upper2, 'main', paste0('figure4_timd4'), height=3, width=5)
# 
# 
# 
# p4_upper = ggpubr::ggarrange(p4_upper1 + theme(plot.margin = unit(c(1,0,0.3,0), 'cm')), 
#                              p4_upper2  + theme(plot.margin = unit(c(1,1.1,0.3,0), 'cm')), 
#                            labels = c('(A) Rare variant burden significance between\n      UK Biobank and meta-analysis of two biobanks', 
#                                       '(C) Meta-analyses across cohorts and biobanks reveal\n      pleiotropic effets of pLoF variants in TIMD4'), 
#                            ncol=2, vjust=1, widths = c(1.5, 1), hjust = 0,
#                            font.label = list(size = 10, color = "black", face = "bold", family = 'Arial'))
# 
# output_figure(p4_upper, 'main', 'figure4_upper', height =4, width = 10)