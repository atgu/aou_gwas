source('~/Dropbox (Partners HealthCare)/github_repo/aou_gwas/R/constants.R')

### GIANT
giant_colors <- c(
  "Significant in All 3" = "#5C4B8A",
  "Significant in AoU and UKB" = "#4F6BD8",
  "Significant in GIANT and UKB" = "#00A1C9",
  "Significant in AoU and GIANT" = "#D55E00",
  "Significant in AoU only" = "#CC79A7",
  "Significant in GIANT only" = "#E69F00",
  "Significant in UKB only" = "#0072B2"
)

giant_comparison_figure <- function(ancestry, ylim=0.2, xlim = 0.1){
  data <- read_gcs_fast(paste0('gs://aou_wlu/v8_analysis/paper_preparation/giant/aou_ukb_', ancestry,'_height_beta_comparison_downsample.txt.bgz'))
  plot_data <- data %>%
    mutate(significant_ukb = Pvalue < 5e-8,
           significant_aou = Pvalue_aou < 5e-8,
           significant_giant = Pvalue_giant < 5e-8) %>%
    mutate(label = case_when(
      significant_ukb & significant_aou & significant_giant ~ 'Significant in All 3',
      significant_ukb & significant_aou ~ 'Significant in AoU and UKB', 
      significant_ukb & significant_giant ~ 'Significant in GIANT and UKB',
      significant_giant & significant_aou ~ 'Significant in AoU and GIANT', 
      significant_aou ~ 'Significant in AoU only',
      significant_giant ~ 'Significant in GIANT only',
      significant_ukb ~ 'Significant in UKB only',
      TRUE ~ 'Not significant'
    )) %>%
    mutate(BETA_giant = if_else(BETA_aou<0, -BETA_giant, BETA_giant),
           BETA_aou = if_else(BETA_aou<0, -BETA_aou, BETA_aou),
    ) 
  p <- plot_data %>%
    ggplot + aes(x=BETA_aou, y=BETA_giant, color= label) +
    geom_errorbar(aes(ymin = BETA_giant- 1.96*as.numeric(SE_giant),ymax = BETA_giant + 1.96*as.numeric(SE_giant)), lwd=0.1, alpha=0.1) +
    geom_errorbarh(aes(xmin = BETA_aou- 1.96*as.numeric(SE_aou),xmax = BETA_aou + 1.96*as.numeric(SE_aou)), lwd=0.1, alpha=0.1) +
    geom_point(size = 1)+
    geom_abline(slope=1, intercept = 0, lty=2) + 
    geom_abline(slope=-1, intercept = 0, lty=2) + 
    labs(x=expression(beta[All~of~Us~v8]), y=expression(beta[GIANT]), color = NULL, pch = NULL) +
    ylim(-ylim, ylim) +
    xlim(0, xlim) +
    geom_segment(aes(yend = -ylim, xend = BETA_aou),
                 arrow = arrow(length = unit(0.2, "cm")),
                 data = plot_data %>% filter(BETA_giant- 1.96*as.numeric(SE_giant) < -ylim), show.legend = F, lwd=0.1) +
    geom_segment(aes(yend = ylim, xend = BETA_aou),
                 arrow = arrow(length = unit(0.2, "cm")),
                 data = plot_data %>% filter(BETA_giant+ 1.96*as.numeric(SE_giant) > ylim), show.legend = F, lwd=0.1) +
    geom_segment(aes(yend = BETA_giant, xend = xlim),
                 arrow = arrow(length = unit(0.2, "cm")),
                 data = plot_data %>% filter(BETA_aou+ 1.96*as.numeric(SE_aou) > xlim), show.legend = F, lwd=0.1) +
    scale_color_manual(name = NULL, values = giant_colors) + 
    scale_fill_manual(name = NULL, values =giant_colors ) +
    themes +
    guides(colour = guide_legend(nrow = 2, byrow = T))
  output_figure(p, 'supplement', paste0('giant_comparison_', ancestry), height =4, width = 6)
  return(p)
}
p1 <- giant_comparison_figure(ancestry='afr', ylim=0.2, xlim = 0.1)
p2 <- giant_comparison_figure(ancestry='eur', ylim=0.5, xlim = 0.2)

p = ggpubr::ggarrange(p2 +
                      theme(plot.margin = unit(c(0.7, 0, 0, 0), "cm")), 
                      p1 +
                      theme(plot.margin = unit(c(0.7, 0, 0, 0), "cm")), 
                      labels = c('(A) Comparison of variant effect sizes on height (EUR-like)',
                                 '(B) Comparison of variant effect sizes on height (AFR-like)'), 
                      nrow=1, vjust=1, widths = c(1, 1), hjust = 0, common.legend = T, 
                      font.label = list(size = 10, color = "black", face = "bold", family = 'Arial'))


output_figure(p, 'supplement', 'supp_figure17_giant', height =4, width = 10)


for(ancestry in c('afr', 'eas', 'eur', 'sas')){
  print(paste0('---------------------', ancestry, '--------------------'))
  data <- read_gcs_fast(paste0('gs://aou_wlu/v8_analysis/paper_preparation/giant/aou_ukb_', ancestry,'_height_beta_comparison_downsample.txt.bgz')) %>%
    filter(AF_Allele2_aou > 0.001 & AF_Allele2_giant > 0.001)
  print('Significant in both: ')
  print(sum(data$Pvalue_giant <= 5e-8 & data$Pvalue_aou <= 5e-8, na.rm=T))
  print('Significant in AoU: ')
  print(sum((data$Pvalue_giant > 5e-8 | is.na(data$Pvalue_giant)) & data$Pvalue_aou <= 5e-8, na.rm=T))
  print('Significant in GIANT: ')
  print(sum(data$Pvalue_giant <= 5e-8 & (data$Pvalue_aou > 5e-8 | is.na(data$Pvalue_aou )), na.rm=T))
}
