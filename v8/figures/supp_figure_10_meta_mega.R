source('~/Dropbox (Partners HealthCare)/github_repo/aou_gwas/R/constants.R')
library(locusviz)
library(magick)

### Meta vs. Mega
manhattan_path = "~/Partners HealthCare Dropbox/Wenhan Lu/aou/axaou/v8/pilot/ACAF_variant/"

p1 <- ggdraw() +
  draw_image(image_read(paste0(manhattan_path, 'ALL/phenotype_3022192.ALL.Pvalue_log10.variant.qq_and_manhattan.png')))+
  theme(plot.margin = unit(c(0.7, 0, 0, 0), "cm"))
p2 <- ggdraw() +
  draw_image(image_read(paste0(manhattan_path, "META/phenotype_3022192.META.Pvalue_log10.variant.qq_and_manhattan.png")))+
  theme(plot.margin = unit(c(0.7, 0, 0, 0), "cm"))

p3 <- ggdraw() +
  draw_image(image_read(paste0(manhattan_path, 'ALL/phenotype_EM_202.2.ALL.Pvalue_log10.variant.qq_and_manhattan.png')))+
  theme(plot.margin = unit(c(0.7, 0, 0, 0), "cm"))
p4 <- ggdraw() +
  draw_image(image_read(paste0(manhattan_path, "META/phenotype_EM_202.2.META.Pvalue_log10.variant.qq_and_manhattan.png")))+
  theme(plot.margin = unit(c(0.7, 0, 0, 0), "cm"))

manhattan_path = "~/Partners HealthCare Dropbox/Wenhan Lu/aou/axaou/v8/ACAF/"

p5 <- ggdraw() +
  draw_image(image_read(paste0(manhattan_path, 'ALL/phenotype_random_0.5_continuous_1.ALL.Pvalue_log10.variant.qq_and_manhattan.png')))+
  theme(plot.margin = unit(c(0.7, 0, 0, 0), "cm"))
p6 <- ggdraw() +
  draw_image(image_read(paste0(manhattan_path, "META/phenotype_random_0.5_continuous_1.META.Pvalue_log10.variant.qq_and_manhattan.png")))+
  theme(plot.margin = unit(c(0.7, 0, 0, 0), "cm"))


p_top = ggpubr::ggarrange(p5, p6, 
                          labels = c('(A) Mega-analysis of simulated quantitative phenotype',
                                     '(B) Meta-analysis of simulated quantitative phenotype'), 
                          nrow=1, vjust=1, widths = c(1, 1), hjust = 0,
                          font.label = list(size = 10, color = "black", face = "bold", family = 'Arial'))


p_middle = ggpubr::ggarrange(p3, p4, 
                             labels = c('(C) Mega-analysis of Type II diabetes',
                                        '(D) Meta-analysis of Type II diabetes'), 
                             nrow=1, vjust=1, widths = c(1, 1), hjust = 0,
                             font.label = list(size = 10, color = "black", face = "bold", family = 'Arial'))

p_bottom = ggpubr::ggarrange(p1, p2, 
                             labels = c('(E) Mega-analysis of Triglycerides',
                                        '(F) Meta-analysis of Triglycerides'), 
                             nrow=1, vjust=1, widths = c(1, 1), hjust = 0,
                             font.label = list(size = 10, color = "black", face = "bold", family = 'Arial'))



p = ggpubr::ggarrange(p_top, p_middle, p_bottom,
                      ncol=1, vjust=1, heights = c(1, 1, 1), hjust = 0,
                      font.label = list(size = 10, color = "black", face = "bold", family = 'Arial'))


output_figure(p, 'supplement', 'supp_figure10_meta_mega', height =5.5, width = 10)
