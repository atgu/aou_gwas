source('~/Dropbox (Partners HealthCare)/github_repo/aou_gwas/R/constants.R')
threshold <- 6.7e-7

meta_comparison <- read_gcs_fast(paste0('gs://aou_wlu/v8_analysis/paper_preparation/summary/aou_ukb_meta_comparison/aou_ukb_gene_signals_Pvalue_Burden.txt.bgz'))
eur_comparison <- read_gcs_fast(paste0('gs://aou_wlu/v8_analysis/paper_preparation/summary/aou_ukb_meta_comparison/aou_ukb_eur_gene_signals_Pvalue_Burden.txt.bgz'))

meta_comparison <- meta_comparison %>%
  dplyr::filter(!is.na(ukb_Pvalue_Burden) & !is.na(aou_Pvalue_Burden)) %>%
  mutate(annotation = if_else(annotation == 'missenseLC', 'missense|LC', annotation)) %>%
  mutate(annotation = if_else(annotation == 'pLoF;missenseLC', 'pLoF|missense|LC', annotation)) %>%
  filter((aou_Pvalue_Burden < threshold | ukb_Pvalue_Burden < threshold) ) %>%
  mutate(annotation = factor(annotation, levels = annotation_types2),
         significance = if_else(aou_Pvalue_Burden < threshold & ukb_Pvalue_Burden < threshold, 'Both', 'Either') )


eur_comparison <- eur_comparison %>%
  dplyr::filter(!is.na(ukb_Pvalue_Burden) & !is.na(aou_Pvalue_Burden)) %>%
  mutate(annotation = if_else(annotation == 'missenseLC', 'missense|LC', annotation)) %>%
  mutate(annotation = if_else(annotation == 'pLoF;missenseLC', 'pLoF|missense|LC', annotation)) %>%
  filter((aou_Pvalue_Burden < threshold | ukb_Pvalue_Burden < threshold) ) %>%
  mutate(annotation = factor(annotation, levels = annotation_types2),
         ukb_BETA_Burden = if_else(aou_BETA_Burden < 0, -ukb_BETA_Burden, ukb_BETA_Burden),
         aou_BETA_Burden = if_else(aou_BETA_Burden < 0, -aou_BETA_Burden, aou_BETA_Burden),
         significance = if_else(aou_Pvalue_Burden < threshold & ukb_Pvalue_Burden < threshold, 'Both', 'Either') ) 

p1 <- eur_comparison %>%
  ggplot + aes(y = ukb_BETA_Burden, x=aou_BETA_Burden, color=annotation, group = annotation, alpha=significance) +
  labs(x = expression(beta[All~by~All~v8~(EUR)]), y = expression(beta[Genebass]), alpha=expression(Significance ~ (6.7 %*% 10^{-7})))+
  geom_point() +
  geom_errorbar(aes(ymin = ukb_BETA_Burden- 1.96*as.numeric(ukb_SE_Burden),ymax = ukb_BETA_Burden + 1.96*as.numeric(ukb_SE_Burden)), lwd=0.1, alpha=0.3) +
  geom_errorbarh(aes(xmin = aou_BETA_Burden- 1.96*as.numeric(aou_SE_Burden),xmax = aou_BETA_Burden + 1.96*as.numeric(aou_SE_Burden)), lwd=0.1, alpha=0.3) +
  xlim(c(0, NA))+
  geom_abline(slope = 1, intercept = 0, lty=2) + 
  geom_abline(slope = -1, intercept = 0, lty=2) + 
  scale_alpha_discrete(range = c(1, 0.3)) +
  annotation_color_scale2 + 
  annotation_fill_scale2 + themes +
  theme(legend.position = 'top',
        legend.direction = "horizontal", 
        # legend.box = "horizontal",
        legend.spacing = unit(0, "pt"),      # spacing between keys
        legend.spacing.y = unit(-2, "pt"),   # vertical spacing between legends
        legend.box.spacing = unit(0, "pt"),   # spacing between groups of legends
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 9),
        axis.title = element_text(size = 8, face = 'bold'),
        axis.text = element_text(size = 6)
  )  + 
  ylim(-0.25, 0.5) +
  guides(
    fill = guide_legend(nrow = 1, byrow = F, reverse = F,  keywidth = 0.5, keyheight = 0.2), 
    color = guide_legend(nrow = 1, byrow = F, reverse = F,  keywidth = 0.5, keyheight = 0.2),
    alpha = guide_legend(nrow = 1, byrow = F, reverse = F,  keywidth = 0.5, keyheight = 0.2)
  )+
  geom_segment(aes(yend = -0.25, xend = aou_BETA_Burden),
               arrow = arrow(length = unit(0.2, "cm")),
               data = eur_comparison %>% filter(ukb_BETA_Burden- 1.96*as.numeric(ukb_SE_Burden) < -0.25), show.legend = F) +
  geom_segment(aes(yend =0.5, xend = aou_BETA_Burden),
               arrow = arrow(length = unit(0.2, "cm")),
               data = eur_comparison %>% filter(ukb_BETA_Burden+ 1.96*as.numeric(ukb_SE_Burden) > 0.5), show.legend = F) 
p1

p2 <- meta_comparison %>%
  ggplot + aes(y = -log10(ukb_Pvalue_Burden), x=-log10(aou_Pvalue_Burden), color=annotation, group = annotation, alpha=significance) +
  labs(x = expression(P[All~by~All~v8~(META)]), y = expression(P[Genebass]))+
  geom_point() + xlim(c(0, NA))+
  geom_abline(slope = 1, intercept = 0, lty=2) + 
  # geom_abline(slope = -1, intercept = 0, lty=2) + 
  scale_alpha_discrete(name = 'Significance (1e-6)', range = c(1, 0.3)) +
  scale_x_continuous(trans = locusviz::trans_loglog_p(), labels=comma, breaks = c(0, 5, 10, 30, 100, 300, 500)) +
  scale_y_continuous(trans = locusviz::trans_loglog_p(), labels=comma, breaks = c(0, 5, 10, 30, 100, 300, 500)) +
  annotation_color_scale2 + 
  annotation_fill_scale2 + themes +
  theme(legend.position = 'top',
        legend.direction = "horizontal", 
        # legend.box = "horizontal",
        legend.spacing.y = unit(-2, 'cm'),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6),
        axis.title = element_text(size = 8, face = 'bold'),
        axis.text = element_text(size = 6)
  )  + 
  guides(
    fill = guide_legend(nrow = 1, byrow = F, reverse = F,  keywidth = 0.5, keyheight = 0.2), 
    color = guide_legend(nrow = 1, byrow = F, reverse = F,  keywidth = 0.5, keyheight = 0.2)
  )


p = ggpubr::ggarrange(p1 + theme(plot.margin = unit(c(0.5,0,0,0), 'cm')), 
                                    p2 + theme(plot.margin = unit(c(0.5,0,0,0), 'cm')), 
                                    labels = c('(A) Burden effect size comparison (Genebass vs. AoU-EUR-like)',
                                               '(B) Burden p-value comparison (Genebass vs. AoU-META)'
                                     ),
                                     hjust=0, vjust = 1, nrow=1, widths = c(0.5, 0.5), common.legend = T, legend = 'top',
                                     font.label = list(size = 8, color = "black", face = "bold", family = 'Arial'))
p
output_figure(p, 'supplement', 'supp_figure18_genebass', height =3, width =8)

