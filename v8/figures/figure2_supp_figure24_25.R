source('~/Dropbox (Partners HealthCare)/github_repo/aou_gwas/R/constants.R')
# Install and load necessary libraries
if (!require("ComplexUpset")) install.packages("ComplexUpset")
if (!require("tidyr")) install.packages("tidyr")
library(ComplexUpset)
library(tidyr)

category_colors <- c('#334195', '#d56f3e', '#7e57c2', '#4f345a','#880d1e', '#43aa8b')
category_labels <- c('Lab\nmeasurements', 'Physical\nmeasurements', 'Phecode', 'PhecodeX', 'Prescriptions', 'Self-reported')
categories <- category_labels
names(category_colors) <- category_labels
category_color_scale<- scale_color_manual(name = NULL, values = category_colors, breaks = names(category_colors), labels = category_labels) 
category_fill_scale <- scale_fill_manual(name = NULL, values = category_colors, breaks = names(category_colors), labels = category_labels)

gene_burden_001 <- rbind(read_gcs_fast('gs://aou_wlu/v8_analysis/paper_preparation/summary/aou_v8_AFR_gene_signals_Pvalue_Burden_0.001_qced.txt.bgz'),
                         read_gcs_fast('gs://aou_wlu/v8_analysis/paper_preparation/summary/aou_v8_AMR_gene_signals_Pvalue_Burden_0.001_qced.txt.bgz'),
                         read_gcs_fast('gs://aou_wlu/v8_analysis/paper_preparation/summary/aou_v8_EAS_gene_signals_Pvalue_Burden_0.001_qced.txt.bgz'),
                         read_gcs_fast('gs://aou_wlu/v8_analysis/paper_preparation/summary/aou_v8_EUR_gene_signals_Pvalue_Burden_0.001_qced.txt.bgz'),
                         read_gcs_fast('gs://aou_wlu/v8_analysis/paper_preparation/summary/aou_v8_MID_gene_signals_Pvalue_Burden_0.001_qced.txt.bgz'),
                         read_gcs_fast('gs://aou_wlu/v8_analysis/paper_preparation/summary/aou_v8_SAS_gene_signals_Pvalue_Burden_0.001_qced.txt.bgz'),
                         read_gcs_fast('gs://aou_wlu/v8_analysis/paper_preparation/summary/aou_v8_META_gene_signals_Pvalue_Burden_0.001_qced.txt.bgz') %>% 
                           select(-N, -n_cases, -n_controls, -ancestries) %>% mutate(ancestry = 'META')
) %>%
  filter(!startsWith(phenoname, 'random') & !(phenoname %in% c('PP_928.1', 'PP_932.11'))) 

lof_upset <- gene_burden_001 %>% filter(annotation == 'pLoF') %>%
  filter(sig) %>%
  distinct()


aou_info_wide <-  read_tsv('~/Dropbox (Partners HealthCare)/aou_v8/data/phenotype/aou_v8_phenotype_info_qced_with_meta.tsv')  %>% 
  filter(pheno_sex == 'both') %>% 
  mutate(presence = 1, ancestry = tolower(ancestry)) %>%
  select(-n_cases, -n_controls) %>%
  pivot_wider(., 
              names_from = ancestry,                # Create columns from `group`
              values_from = presence,            # Use the `presence` column
              values_fill = 0  # Ensure missing values are filled with 0
  ) %>%
  mutate(all_6 = (afr+amr+eas+eur+mid+sas == 6),
         big_3_only = (afr+amr+eur == 3) & (afr+amr+eas+eur+mid+sas == 3),
         big_3_at_least = (afr+amr+eur == 3)) %>%
  select(phenoname, all_6, big_3_only, big_3_at_least)


data <- lof_upset %>% 
  merge(., pheno_info_full %>% select(phenoname, disease_category) %>% distinct(), by = 'phenoname') %>%
  dplyr::filter(sig & !(phenoname %in% c('PP_928.1', 'PP_932.11')) & annotation == 'pLoF') %>%
  dplyr::mutate(ancestry = factor(ancestry, levels=toupper(names(pops)), labels=pops), presence = 1)

# Reshape the data to wide format
data_wide <- pivot_wider(
  data %>% select(-gene_symbol, -Pvalue_Burden, -BETA_Burden) %>% distinct(),
  names_from = ancestry,                # Create columns from `group`
  values_from = presence,            # Use the `presence` column
  values_fill = 0  # Ensure missing values are filled with 0
) %>%
  mutate(category = if_else(category %in% c('pfhh_survey', 'mhwb_survey'), 'self-reported', category))%>%
  mutate(category = factor(category, levels = c('lab_measurement', 'physical_measurement', 'mcc2_phecodex', 'r_drug', 'self-reported'), 
                           labels=c("Lab\nmeasurements", "Physical\nmeasurements", "PhecodeX", "Prescriptions", 'Self-reported'))) %>%
  merge(., aou_info_wide, by = 'phenoname')

p_top <- upset(
  data_wide %>% dplyr::filter(big_3_at_least),
  intersect = c('AFR-like', 'AMR-like', 'EAS-like', 'EUR-like', 'META'),
  name = '',
  queries = list(
    upset_query(set = 'AFR-like', fill = pop_colors["afr"]),
    upset_query(set ='AMR-like', fill = pop_colors["amr"]),
    upset_query(set = 'EAS-like', fill = pop_colors["eas"]),
    upset_query(set = 'EUR-like', fill = pop_colors["eur"]),
    # upset_query(set = pops['mid'], fill = pop_colors["mid"]),
    # upset_query(set = pops['sas'], fill = pop_colors["sas"]),
    upset_query(set = 'META', fill = pop_colors["meta"])
  ),
  base_annotations=list(
    ' '=(
      intersection_size(
        text=list(
          vjust=-0.1,
          size = 6
        ),
        mapping=aes(fill=category),
        bar_number_threshold=1,  # show all numbers on top of bars
        width=0.5,   # reduce width of the bars
      )
      # add some space on the top of the bars
      + scale_y_continuous(expand=expansion(mult=c(0, 0.05)), limits = c(0, 280))
      + scale_fill_manual(name = NULL, values = category_colors, breaks = names(category_colors), labels = category_labels)
      +labs(x=NULL, y=NULL)
      + theme(
        # hide grid lines
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        # show axis lines
        axis.line.x=element_line(colour='black', linewidth = 0.1),
        axis.line.y=element_blank(),
        axis.text.y=element_blank(),
        legend.text = element_text(size=13, color = 'black'),
        legend.position = c(0.9, 0.72),
        legend.direction = 'vertical'
      )+
        guides(
          fill = guide_legend(keywidth = 0.8, keyheight = 1.5)
        )
    )
  ),
  stripes=upset_stripes(
    geom=geom_segment(size=5),  # make the stripes larger
    colors=c('grey95', 'white')
  ),
  # to prevent connectors from getting the colorured
  # use `fill` instead of `color`, together with `shape='circle filled'`
  matrix=intersection_matrix(
    geom=geom_point(
      shape='circle filled',
      size=3.5,
      stroke=0.45
    ),
  )+ labs(x=NULL, y=NULL),
  # set_sizes=(
  #   upset_set_size(geom=geom_bar(width = 0.8), position='right') + 
  #     geom_text(aes(label=..count..), hjust=1, stat='count', size = 3, color = 'white') +
  #     theme(
  #     # axis.line.x=element_line(colour='black'),
  #     # axis.ticks.x=element_line()
  #       panel.grid.major = element_blank(),
  #       panel.grid.minor = element_blank(),
  #       # strip.background = element_blank(),
  #     axis.text.x = element_blank()
  #   ) + labs(x = NULL, y= NULL) + scale_y_log10()
  # ),
  # guides='over',
  
  set_sizes= FALSE,
  themes=upset_modify_themes(
    list(
      'intersections_matrix'=theme(text=element_text(size=15, color='black')),
      'overall_sizes'=theme(axis.text.x=element_text(size=15, vjust = 1))
    )
  ),
  # annotations =list(
  #       'Phenotype category'=list(
  #           aes=aes(x=intersection, fill=category),
  #           geom=list(
  #               geom_bar(stat='count', position='fill', na.rm=TRUE),
  #               scale_fill_manual(name = NULL, values = category_colors, breaks = names(category_colors), labels = category_labels)
  #           )
  #       )
  #   )
  # sort_sets='ascending',
  # sort_intersections='ascending'
) + labs(x=NULL, y = NULL) + axa_themes
p_top
output_figure(p_top, 'main', 'figure2_upper_upset_plot_CAF_weight_v2', 4, 10)

multi_ancestry <- rbind(read_gcs_fast('gs://aou_wlu/v8_analysis/paper_preparation/summary/aou_v8_AFR_meta_component_Pvalue_Burden_raw.txt.bgz') %>% select(-N_eff),
                        read_gcs_fast('gs://aou_wlu/v8_analysis/paper_preparation/summary/aou_v8_AMR_meta_component_Pvalue_Burden_raw.txt.bgz') %>% select(-N_eff),
                        read_gcs_fast('gs://aou_wlu/v8_analysis/paper_preparation/summary/aou_v8_EAS_meta_component_Pvalue_Burden_raw.txt.bgz') %>% select(-N_eff),
                        read_gcs_fast('gs://aou_wlu/v8_analysis/paper_preparation/summary/aou_v8_EUR_meta_component_Pvalue_Burden_raw.txt.bgz') %>% select(-N_eff),
                        read_gcs_fast('gs://aou_wlu/v8_analysis/paper_preparation/summary/aou_v8_MID_meta_component_Pvalue_Burden_raw.txt.bgz') %>% select(-N_eff),
                        read_gcs_fast('gs://aou_wlu/v8_analysis/paper_preparation/summary/aou_v8_SAS_meta_component_Pvalue_Burden_raw.txt.bgz') %>% select(-N_eff)
) %>%
  distinct()
meta_explode <- read_gcs_fast('gs://aou_wlu/v8_analysis/paper_preparation/summary/aou_v8_meta_gene_signals_Pvalue_Burden_0.001_qced_exploded.txt.bgz') %>%
  mutate(meta_p = Pvalue_Burden, meta_N = N, ancestry = tolower(ancestry)) %>%
  select(-Pvalue_Burden, -N) %>%
  dplyr::filter(sig ) %>%
  distinct()
library(stringr)

truncate_sentence <- function(x, n = 4) {
  vapply(x, function(s) {
    if (is.na(s) || !nzchar(s)) return(s)
    # normalize any Unicode spaces to regular spaces, then squish
    s <- str_replace_all(s, "\\p{Z}+", " ")
    s <- str_squish(s)
    words <- str_split(s, " ", simplify = FALSE)[[1]]
    if (length(words) > n) {
      paste(paste(words[seq_len(n)], collapse = " "), "...")
    } else {
      s
    }
  }, character(1))
}


full_data <- meta_explode %>%
  dplyr::filter(annotation == 'pLoF') %>%
  mutate(ancestry = toupper(ancestry)) %>%
  merge(., multi_ancestry %>% filter(max_MAF == 0.001), all.x = T, by = c('gene_id', 'gene_symbol', 'annotation', 'phenoname', 'description', 'category', 'max_MAF', 'ancestry')) %>%
  dplyr::mutate(description = str_split(description, '\\[') %>% map_chr(., 1)) %>%
  mutate(description = trimws(description, which = "right")) %>%
  # dplyr::mutate(association = paste0(truncate_sentence(str_to_sentence(description)), if_else(!(category %in% c('physical_measurement', 'lab_measurement')), paste0('|', str_pad(phenoname, width = 9, side = "left")), '') , '|', str_pad(gene_symbol, width = 8, side = "left"))) %>%
  dplyr::mutate(association = paste0(str_to_sentence(description), if_else(!(category %in% c('physical_measurement', 'lab_measurement')), paste0('|', str_pad(phenoname, width = 9, side = "left")), '') , '|', str_pad(gene_symbol, width = 8, side = "left"))) %>%
  mutate(association = if_else(grepl('itamin d', association), str_replace(association, 'itamin d', 'itamin D'), association)) %>%
  mutate(association = if_else(grepl('itamin b', association), str_replace(association, 'itamin b', 'itamin B'), association)) %>%
  mutate(association = if_else(grepl('itamin c', association), str_replace(association, 'itamin c', 'itamin C'), association)) %>%
  mutate(association = if_else(grepl('itamin a', association), str_replace(association, 'itamin a', 'itamin A'), association)) %>%
  mutate(association = if_else(grepl('Mch', association), str_replace(association, 'Mch', 'MCH'), association)) %>%
  mutate(association = if_else(grepl('hdl', association), str_replace(association, 'hdl', 'HDL'), association)) %>%
  dplyr::filter(complete.cases(.)) %>%
  dplyr::group_by(association) %>%
  mutate(n_big_3 = sum(as.integer(ancestry %in% c('AFR', 'AMR', 'EUR')))) %>%
  filter(n_big_3 == 3) %>%
  dplyr::mutate(unweighted_Z = weighted_Z/sqrt(N),
         total_unweighted_Z = sum(unweighted_Z),
         total_Z = sum(weighted_Z),
         oppo_sign = any(weighted_Z > 0) + any(weighted_Z < 0)) %>%
  # filter(oppo_sign == 1) %>%
  dplyr::mutate(prop_Z = weighted_Z/total_Z,
         prop_unweighted_Z = unweighted_Z/total_unweighted_Z) %>%
  dplyr::mutate(total_prop_Z = sum(prop_Z),
         total_prop_unweighted_Z = sum(prop_unweighted_Z)) %>%
  merge(., aou_info_wide, by = 'phenoname') %>%
  distinct() %>%
  dplyr::arrange(desc(total_unweighted_Z))


full_data_filtered_table <- full_data %>%
  dplyr::group_by(association) %>%
  dplyr::mutate(n_significant = sum(Pvalue_Burden < 6.7e-7)) %>%
  dplyr::filter(n_significant == 0) %>%
  # dplyr::filter(meta_N > 150000) %>%
  # dplyr::filter(meta_N > 280000 | category == 'lab_measurement') %>%
  dplyr::mutate(prop_N = N/sum(N),
         prop_Z = unweighted_Z/sum(unweighted_Z))
print(nrow(full_data_filtered_table %>% select(gene_id, gene_symbol, annotation, phenoname) %>% distinct()))
print(nrow(full_data %>% select(gene_id, gene_symbol, annotation, phenoname) %>% distinct()))
write_tsv(full_data_filtered_table, '~/Dropbox (Partners HealthCare)/aou_v8/tables/supplementary_data3.tsv')

full_data_filtered <- full_data %>%
  dplyr::group_by(association) %>%
  dplyr::mutate(n_significant = sum(Pvalue_Burden < 6.7e-7)) %>%
  dplyr::filter(n_significant == 0) %>%
  # dplyr::filter(meta_N > 150000) %>%
  # dplyr::filter(meta_N > 280000 | category == 'lab_measurement') %>%
  dplyr::filter(category != 'physical_measurement') %>%
  dplyr::mutate(prop_N = N/sum(N),
                prop_Z = unweighted_Z/sum(unweighted_Z))


p_bottom <- full_data_filtered %>%
  mutate(ancestry = tolower(ancestry))%>%
  dplyr::mutate(association = factor(association, levels = full_data %>% select(association) %>% distinct() %$% association)) %>%
  dplyr::mutate(category = if_else(endsWith(category, 'measurement'), 'Lab measurements', category)) %>%
  dplyr::mutate(category = factor(category, levels = c('Lab measurements', 'mcc2_phecodex', 'r_drug', 'pfhh_survey', 'mhwb_survey'), labels=c("Lab measurements", "PhecodeX", "Prescriptions", 'Personal family\nhealth history (PFHH)', 'Mental Health Wellbeing (MHWB)'))) %>%
  dplyr::filter(category %in% c('PhecodeX', "Lab measurements")) %>%
  # dplyr::filter(meta_N > 280000 | category != 'PhecodeX')  %>%
  ggplot + aes(x = association, y = unweighted_Z, fill = ancestry) +
  geom_bar(stat = "identity") +
  labs(
    x = 'Association',
    y = "Unweighted Z score",
    fill = NULL
  )+
  geom_hline(yintercept = 0,  color = "black", lwd = 0.2) +
  theme_classic() +
  themes + 
  scale_color_manual(name = NULL, values = pop_colors, breaks = names(pops), labels = pops)+
  scale_fill_manual(name = NULL, values = pop_colors, breaks = names(pops), labels = pops) + 
  theme(
    legend.position = 'bottom',
    axis.text.y = element_text(size = 11, family = "mono"),
    axis.title = element_text(family='Arial', face = 'bold'),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.text = element_text(size=10, face = 'plain'),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12),
    plot.margin = unit(c(0.8, 0, 0,0.5), 'cm')
  ) + 
  coord_flip() + 
  facet_wrap(category~., ncol=3, scales = 'free', strip.position ="top")+
  guides(
    fill = guide_legend(nrow = 1, byrow = TRUE, reverse = FALSE,  keywidth = 1, keyheight = 0.3), 
    color = guide_legend(nrow = 1, byrow = TRUE, reverse = FALSE,  keywidth = 1, keyheight = 0.3)
  )
p_bottom

output_figure(p_bottom, 'main', 'figure2_lower_component_plot', 10, 15)

pp = ggpubr::ggarrange(p_top+ theme(plot.margin = unit(c(1.2, 0, 0, 0), 'cm')), 
                       p_bottom + theme(plot.margin = unit(c(0.8, 0, 0.3, 0), 'cm')), ncol=1, hjust = 0, heights =c(1, 2), vjust=1,
                       labels = c('a', 
                                  'b'),
                       font.label = list(size = 20, color = "black", face = "bold", family = 'Arial'), common.legend = T, legend = 'bottom'
)
output_figure(pp, 'main', 'figure2', height=350/25.4, width=183*2/25.4)





full_data %>%
  group_by(association) %>%
  dplyr::mutate(prop_N = N/sum(N),
                prop_Z = unweighted_Z/sum(unweighted_Z))%>%
  dplyr::group_by(ancestry) %>%
  dplyr::summarize(avg_Z = mean(prop_Z),
                   avg_N = mean(prop_N)) %>% 
  dplyr::summarize(tidy(cor.test(avg_Z, avg_N)))

p_avg <- full_data %>%
  group_by(association) %>%
  dplyr::mutate(
    prop_N = N / sum(N),
    prop_Z = unweighted_Z / sum(unweighted_Z)
  ) %>%
  group_by(ancestry) %>%
  dplyr::summarize(
    avg_Z = mean(prop_Z),
    avg_N = mean(prop_N)
  ) %>%
  mutate(ancestry = tolower(ancestry)) %>%
  ggplot() +
  aes(x = avg_N, y = avg_Z, color = ancestry) +
  geom_abline(
    intercept = 0, slope = 1,
    linetype = "dashed", color = "grey60", linewidth = 0.5
  ) +
  geom_point(size = 3.5, alpha = 0.85) +
  scale_color_manual(
    name   = NULL,
    values = pop_colors,
    breaks = names(pops),
    labels = pops
  ) +
  scale_fill_manual(
    name   = NULL,
    values = pop_colors,
    breaks = names(pops),
    labels = pops
  ) +
  labs(
    x = "Average proportion of sample size across phenotypes",
    y = "Average proportion of unweighted\nZ-scores across phenotypes"
  ) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, NA)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, NA)) +
  themes +
  theme(
    axis.title         = element_text(size = 11, color = "grey20"),
    axis.text          = element_text(color = "grey30"),
    plot.margin        = unit(c(0.5, 0.5, 0.3, 0.3), "cm"),
    legend.position    = c(0.72, 0.28),
    legend.background  = element_rect(fill = alpha("white", 0.8), color = NA),
    legend.text        = element_text(size = 7),
    legend.key.size    = unit(0.4, "cm"),
    legend.spacing.y   = unit(0.1, "cm")
  )

p_avg
output_figure(p_avg, 'supplement', 'supp_figure25_avg_proportion_z_scores_vs_avg_n', height=3, width=5)


p <- full_data_filtered  %>%
  mutate(ancestry = tolower(ancestry))%>%
  dplyr::mutate(association = factor(association, levels = full_data %>% select(association) %>% distinct() %$% association)) %>%
  dplyr::mutate(category = if_else(endsWith(category, 'measurement'), 'Lab measurements', category)) %>%
  dplyr::mutate(category = factor(category, levels = c('Lab measurements', 'mcc2_phecodex', 'r_drug', 'pfhh_survey', 'mhwb_survey'), labels=c("Lab measurements", "PhecodeX", "Prescriptions", 'Personal family\nhealth history (PFHH)', 'Mental Health Wellbeing (MHWB)'))) %>%
  dplyr::filter(category %in% c('PhecodeX', "Lab measurements")) %>%
  ggplot + aes(x = association, y = prop_Z, fill = ancestry) +
  geom_bar(stat = "identity") +
  labs(
    x = 'Association',
    y = "Unweighted Z",
    fill = NULL
  )+
  geom_hline(yintercept = 0,  color = "black", lwd = 0.2) +
  theme_classic() +
  themes + 
  scale_color_manual(name = NULL, values = pop_colors, breaks = names(pops), labels = pops)+
  scale_fill_manual(name = NULL, values = pop_colors, breaks = names(pops), labels = pops) + 
  theme(
    # axis.ticks.y = element_blank(),
    # axis.line.y = element_blank(),
    legend.position = 'bottom',
    axis.text.y = element_text(size = 8, family = "mono"),
    axis.title = element_text(family='Arial', face = 'bold'),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    # axis.line.y = element_line(linewidth = 0.2),
    strip.text = element_text(size=10, face = 'bold'),
    # axis.text.x = element_text(angle = 45, hjust = 1),
    # axis.ticks.x = element_line(linewidth = 0.1),
    # axis.line.x = element_line(linewidth = 0.1),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8),
    plot.margin = unit(c(0.8, 0, 0,0.5), 'cm')
  ) + 
  coord_flip() + 
  # facet_grid(~category, scales = "free_x", space = "free_x") +
  # theme_classic(base_family = "mono") +
  facet_wrap(category~., ncol=3, scales = 'free', strip.position ="top")+
  guides(
    fill = guide_legend(nrow = 1, byrow = TRUE, reverse = FALSE,  keywidth = 1, keyheight = 0.3), 
    color = guide_legend(nrow = 1, byrow = TRUE, reverse = FALSE,  keywidth = 1, keyheight = 0.3)
  ) + 
  theme(legend.position = 'bottom')
p
output_figure(p, 'supplement', 'supp_figure24_proportion_z_scores', height=12, width=12)

