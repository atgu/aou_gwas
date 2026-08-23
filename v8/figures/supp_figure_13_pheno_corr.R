source('~/Dropbox (Partners HealthCare)/github_repo/aou_gwas/R/constants.R')

pheno_corr <- read_gcs_fast('gs://aou_wlu/v8_analysis/pheno_corr.txt.bgz')

p1 <-pheno_corr%>%
  ggplot() + 
  aes(x = entry^2, fill = after_stat(x)) +
  geom_histogram(binwidth = 0.01, alpha = 0.85, color = "white", linewidth = 0.1) +
  scale_fill_distiller(palette = "Blues", direction = 1, guide = "none") +
  scale_y_log10(labels = comma, expand = expansion(mult = c(0, 0.05))) +
  scale_x_continuous(breaks = seq(0, 1, 0.1), expand = expansion(mult = c(0.01, 0.01))) +
  labs(
    y = "Number of phenotype pairs",
    x = expression(bold(paste("Phenotypic correlation (", r^2, ")")))
  ) +
  theme_classic(base_size = 11) +
  themes +
  theme(
    axis.title = element_text(size = 9, face = "bold"),
    axis.text = element_text(size = 8, color = "grey30"),
    axis.ticks = element_line(color = "grey70", linewidth = 0.3),
    axis.line = element_line(color = "grey40", linewidth = 0.4),
    panel.grid.major.y = element_line(color = "grey92", linewidth = 0.3),
    plot.margin = margin(10, 15, 10, 10)
  )
p1

corr_count = data.frame(
  Count = c(3848, 3374, 3059, 2790, 2502, 2219, 1921, 1610, 1288), 
  Corr = seq(0.1, 0.9, 0.1))


p2 = corr_count %>%
  ggplot() + 
  aes(x = Corr, y = Count, fill = Corr) +
  geom_bar(stat = "identity", alpha = 0.85, color = "white", linewidth = 0.1) +
  scale_fill_distiller(palette = "Blues", direction = 1, guide = "none") +
  geom_text(aes(label = comma(Count)), vjust = -0.3, size = 2.5, color = "grey30") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)), label = comma) +
  scale_x_continuous(breaks = seq(0.1, 0.9, 0.1)) +
  labs(y = "Phenotypes removed", x = "Correlation threshold") +
  theme_classic(base_size = 11) +
  themes +
  theme(
    axis.title = element_text(size = 9, face = "bold"),
    axis.text = element_text(size = 8, color = "grey30"),
    axis.ticks = element_line(color = "grey70", linewidth = 0.3),
    axis.line = element_line(color = "grey40", linewidth = 0.4),
    panel.grid.major.y = element_line(color = "grey92", linewidth = 0.3),
    plot.margin = margin(10, 15, 10, 10)
  )
p2

p = ggpubr::ggarrange(p1, p2,  labels = c('(A) Number of phenotype pairs by correlation', '(B) Number of phenotypes to remove by r2 threshold'), 
                      nrow = 1, align = 'h', hjust = 0, vjust = 1, 
                      font.label = list(size = 10, color = "black", face = "bold", family = 'Arial'))
output_figure(p, 'supplement', 'supp_figure13_pheno_corr', height=3, width=8)
