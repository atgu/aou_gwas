source('~/Dropbox (Partners HealthCare)/github_repo/aou_gwas/R/constants.R')

empirical_p <- read_gcs_fast("gs://aou_wlu/v8_analysis/empirical_min_p_values_per_ancestry_per_test_max_MAF_1e_3.tsv", format = "tsv")
empirical_p %>%
  group_by(Ancestry, Test_type) %>%
  dplyr::summarize(empirical_p = 0.05*median(Min_P)) 
empirical_p_group <- empirical_p %>%
  group_by(Ancestry, Test_type) %>%
  dplyr::summarize(empirical_p = 0.05*median(Min_P)) 


p <- empirical_p_group %>%
  mutate(test = factor(Test_type, levels = c('genome_variant', 'exome_variant', 'burden', 'skato', 'skat'), 
                       labels = c('ACAF', 'Exome', 'Burden', 'SKATO', 'SKAT'))) %>%
  ggplot + aes(x = test, y = -log10(empirical_p), color = Ancestry) +
  geom_point() + 
  geom_hline(yintercept = -log10(6.7e-7), lty = 2) + 
  geom_hline(yintercept = -log10(2.5e-7), lty = 2) + 
  geom_hline(yintercept = -log10(8e-9), lty = 2) + 
  coord_cartesian(clip = "off") +
  annotate(geom = 'text', x = 6, y = -log10(6.7e-7), label = 'Genebass: Burden', vjust = -0.8, size = 3) + 
  annotate(geom = 'text', x = 6, y = -log10(2.5e-7), label = 'Genebass: SKATO', vjust = -0.8, size = 3) + 
  annotate(geom = 'text', x = 6.3, y = -log10(8e-9), label = 'Genebass: Exome variants', vjust = -0.8, size = 3) + 
  labs(x = NULL, y = expression(Empirical~P~threshold~(-log[10]))) + 
  scale_color_manual(name = NULL,values = pop_colors, breaks = c('afr', 'amr', 'eas', 'eur', 'mid', 'sas','meta'), labels =pops)+
  scale_fill_manual(name = NULL,values = pop_colors, breaks = c('afr', 'amr', 'eas',  'eur','mid', 'sas','meta'), labels = pops) +
  themes +
  theme(legend.position = 'top',legend.margin = margin(0, 0, 0, 0),
        legend.text = element_text(size = 8), legend.title = element_text(size = 8),
        plot.margin = margin(t = 0, r = 60, b = 0, l = 0, unit = "pt")) +
  guides(
    color = guide_legend(nrow = 2, byrow = TRUE, reverse = F,  keywidth = 1, keyheight = 0.3)
  ) 
p
output_figure(p, 'supplement','supp_figure16_empirical_p', height=3, width=5)
