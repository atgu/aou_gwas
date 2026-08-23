source('~/Dropbox (Partners HealthCare)/github_repo/aou_gwas/R/constants.R')
library(ggsci)

create_long_burden_data <- function(p, maxMAF){
  print(paste0('----------', p, '--------', maxMAF, '----------'))
  unique <- read_gcs_fast(paste0('gs://aou_wlu/v8_analysis/paper_preparation/summary/aou_v8_META_gene_signals_', p,'_', maxMAF,'_qced.txt.bgz'))  
  unique$P <- unique[[p]]
  unique <- unique %>% 
    filter(P < 6.7e-7)
  unique_long <- unique %>%
    select(gene_id, gene_symbol, annotation, max_MAF, phenoname, description, P) %>%
    mutate(ancestry = 'META')
  for(ancestry in ANCESTRIES[1:6]){
    ancestry = toupper(ancestry)
    print(paste0('#########', ancestry, ':'))
    
    tmp <- read_gcs_fast(paste0('gs://aou_wlu/v8_analysis/paper_preparation/summary/aou_v8_', ancestry,'_gene_signals_', p,'_', maxMAF,'_qced.txt.bgz')) 
    if(nrow(tmp) == 0){
      unique[[paste0('ancestry.', tolower(ancestry))]] = NA
      next
    } 
    tmp$P <- tmp[[p]]
    tmp <- tmp %>% 
      filter(P < 6.7e-7) %>%
      select(gene_id, gene_symbol, annotation, max_MAF, phenoname, description, P) %>%
      mutate(ancestry = ancestry)
    unique_long = rbind(unique_long, tmp)
  }
  unique_long <- unique_long %>%
    filter(!startsWith(phenoname, 'random') & !(phenoname %in% c('PP_928.1', 'PP_932.11')))
  return(unique_long)
}

add_color_scale <- function(figure, name) {
  if (name == "ancestry") {
    figure <- figure +
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
      guides(
        fill  = guide_legend(nrow = 1, byrow = T, keywidth = 1.2, keyheight = 0.4),
        color = guide_legend(nrow = 1, byrow = T, keywidth = 1.2, keyheight = 0.4)
      )
  } else {
    figure <- figure +
      scale_color_jama(
        name   = "Trait type",
        breaks = c("continuous", "binary"),
        labels = c("Continuous", "Binary")
      ) +
      scale_fill_jama(
        name   = "Trait type",
        breaks = c("continuous", "binary"),
        labels = c("Continuous", "Binary")
      ) +
      guides(
        fill  = guide_legend(nrow = 2, byrow = TRUE, keywidth = 1.2, keyheight = 0.4),
        color = guide_legend(nrow = 2, byrow = TRUE, keywidth = 1.2, keyheight = 0.4)
      )
  }
  
  return(figure)
}

plot_cdf <- function(data,
                     field,
                     group       = NULL,
                     save        = FALSE,
                     output_path = "~/Dropbox (Partners HealthCare)/aou/phenotype/figures/phenotype_cdf/") {
  
  if (is.null(group)) {
    cdf_data <- data %>%
      dplyr::arrange(get(field)) %>%
      dplyr::mutate(
        id      = row_number(),
        n_total = dplyr::n(),
        cdf     = id / n_total
      )
    
    figure <- cdf_data %>%
      ggplot() +
      aes(x = get(field), y = cdf)
    
  } else {
    cdf_data <- data %>%
      dplyr::arrange(get(field)) %>%
      dplyr::group_by(get(group)) %>%
      dplyr::mutate(
        id      = row_number(),
        n_total = dplyr::n(),
        cdf     = id / n_total
      )
    View(cdf_data)
    
    figure <- cdf_data %>%
      ggplot() +
      aes(x = get(field), y = cdf * 100, color = get(group))
    
    figure <- add_color_scale(figure, group)
  }
  
  figure <- figure +
    labs(
      x = str_to_title(str_replace_all(field, "_", " ")),
      y = "Proportion"
    ) +
    geom_line() +
    geom_hline(yintercept = 100, lty = 2) + 
    themes +
    labs(
      x        = "Number of Associations",
      y        = "Phenotype Percentile",
    ) +
    scale_x_log10(labels = scales::comma) +
    geom_line(linewidth = 0.9, alpha = 0.85)  +
    theme_classic()+
    theme(
      plot.title       = element_text(face = "bold", size = 14, margin = margin(b = 4)),
      plot.subtitle    = element_text(color = "grey40", size = 10, margin = margin(b = 10)),
      plot.margin      = unit(c(0.5, 0.5, 0.3, 0.3), "cm"),
      axis.title       = element_text(size = 11),
      axis.text        = element_text(color = "grey30"),
      legend.position  = "top",
      legend.title     = element_text(face = "bold", size = 10),
      legend.text      = element_text(size = 9)
    )
  
  if (save) {
    suffix <- if_else(is.null(group), "", paste0("_by_", group))
    png(
      paste0(output_path, "cdf_", field, suffix, ".png"),
      width  = 5,
      height = 3,
      units  = "in",
      res    = 300
    )
    print(figure)
    dev.off()
  }
  
  return(figure)
}


gwas_ld_pruned <- read_tsv('~/Dropbox (Partners HealthCare)/aou_v8/data/all_clumped_annotated.tsv') %>%
  merge(., full_pheno_sum %>% select(phenoname, description, category, ancestry) %>% filter(! (category %in% c('onset', 'progression'))) %>% distinct(), 
        all.y = T, by.x = c('Phenotype', 'ancestry'), by.y = c('phenoname', 'ancestry'))

gwas_ld_pruned_pheno_assoc <- gwas_ld_pruned %>%
  group_by(ancestry, Phenotype, description, category) %>%
  dplyr::summarize(
    n_assoc = sum(P < 5e-8, na.rm = T)
  ) %>%
  mutate(ancestry  = tolower(ancestry)) %>% 
  filter(!startsWith(Phenotype, 'random') & !(Phenotype %in% c('PP_928.1', 'PP_932.11')))

burden_001_qced <- create_long_burden_data('Pvalue_Burden', '0.001')

burden_001_qced_pheno_assoc <- burden_001_qced %>%
  merge(., full_pheno_sum %>% select(phenoname, description, category, ancestry) %>% filter(! (category %in% c('onset', 'progression'))) %>% distinct(), 
        all.y = T, by.x = c('phenoname', 'ancestry', 'description'), by.y = c('phenoname', 'ancestry', 'description'))%>%
  group_by(ancestry, phenoname, description, category) %>%
  dplyr::summarize(
    n_assoc = sum(P < 6.7e-7, na.rm = T)
  ) %>%
  mutate(ancestry  = tolower(ancestry)) %>% 
  filter(!startsWith(phenoname, 'random') & !(phenoname %in% c('PP_928.1', 'PP_932.11')))


p1 <- plot_cdf(
  gwas_ld_pruned_pheno_assoc,
  field       = "n_assoc",
  group       = "ancestry",
  save        = FALSE,
  output_path = ""
)  + ylim(25, 100) 

p2 <- plot_cdf(
  burden_001_qced_pheno_assoc,
  field       = "n_assoc",
  group       = "ancestry",
  save        = FALSE,
  output_path = ""
) + ylim(85, 100) 


p = ggpubr::ggarrange(p1, p2, 
                      nrow=1, hjust = 0, widths = c(1,1),align = "h", 
                      labels = c('(A) Single-variant associations', '(B) Gene burden associations'), 
                      common.legend = T, font.label = list(size = 10, color = "black", face = "bold", family = 'Arial')
)
output_figure(p, 'supplement', 'supp_figure23_pheno_assoc_cdf', height=3.5, width=8)
