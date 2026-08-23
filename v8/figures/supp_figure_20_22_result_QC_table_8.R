source('~/Dropbox (Partners HealthCare)/github_repo/aou_gwas/R/constants.R')
figure_root <- paste0(aou_qc_fig_path, 'qc/')

# for(anc in c('afr', 'amr', 'eur', 'meta')){
#   for(t in c('exome_gene_0.001', 'genome_variant', 'exome_variant')){
#     data <- read_gcs_fast(paste0('gs://aou_analysis/v8/data/qc/aou_', anc,'_',t, '_phenotype_lambda_gc_filtered.txt.bgz')) %>%
#       filter(phenoname %in% c('PP_928.1', 'PP_932.11')) %>%
#       select(phenoname, lambda_gc_raw, lambda_gc_filtered)
#     print(anc)
#     print(t)
#     print(data)
#   }
# }

plot_lambda_gc_boxplot <- function(test_type){
  plt_lst <- list()
  for(pop in c('afr', 'amr', 'eas', 'eur', 'mid', 'sas')){
    plt_name <- paste0(pop, '_', test_type, '_filtered_QCed')
    print(plt_name)
    pheno_info_full <- read_tsv('~/Dropbox (Partners HealthCare)/aou_v8/data/phenotype/aou_v8_phenotype_info_qced_with_meta.tsv') %>%
      filter(ancestry == toupper(pop) & pheno_sex == 'both') %>%
      select(phenoname, category, description, n_cases, n_controls) %>%
      distinct() 
    pheno_info_full <- pheno_info_full %>% select(-n_cases, -n_controls)
    if(!grepl('gene', test_type)){
      pheno_info_full <- pheno_info_full %>% select(-category, -description)
    }
    data <- read_gcs_fast(paste0('gs://aou_analysis/v8/data/qc/aou_', pop,'_', type, '_phenotype_lambda_gc_filtered.txt.bgz')) %>%
      merge(., pheno_info_full, by = 'phenoname')%>%
      mutate(category = if_else(category %in% c('pfhh_survey', 'mhwb_survey'), 'self-reported', category))%>%
      filter(category != 'random_phenotype') %>%
      mutate(category = factor(category, levels = c('lab_measurement', 'physical_measurement', 'mcc2_phecode', 
                                                    'mcc2_phecodex', 'r_drug', 'self-reported', 'onset', 'progression', 'random_phenotype'),
                               labels = c('Lab measurements', 'Physical measurements', 'Phecode', 
                                          'PhecodeX', 'Prescriptions', 'Self-reported',
                                          'Onset', 'Progression', 'Random phenotype')),
             lambda_gc = lambda_gc_filtered) %>%
      filter(!phenoname %in% c('PP_932.11', 'PP_928.1', '3011099', '3046664', '3013861', '3052648', '3045792_3016921'))
    p_qc <- data %>%
      filter(lambda_gc_raw < 2) %>%
      select(lambda_gc_raw, lambda_gc_filtered, category, n_cases) %>%
      pivot_longer(cols = starts_with('lambda_gc_'), names_to = 'QC', values_to = 'Lambda') %>%
      mutate(QC = factor(if_else(QC == 'lambda_gc_filtered', 'QCed', 'Raw'), levels = c('Raw', 'QCed'))) %>%
      ggplot() + 
      aes(x = QC, y = Lambda, color = category, fill = category, group = interaction(QC, category)) +
      geom_hline(yintercept = 1, linetype = 2, color = "grey50", linewidth = 0.4) +
      geom_boxplot(width = 0.7, linewidth = 0.6, alpha = 0.15, outlier.size = 0.8) +
      scale_y_continuous(limits = c(NA, 1.25), expand = expansion(mult = c(0.02, 0.02))) +
      scale_color_manual(name = NULL, values = category_colors, breaks = names(category_colors), labels = category_labels) +
      scale_fill_manual(name = NULL, values = category_colors, breaks = names(category_colors), labels = category_labels) +
      labs(y = expression(lambda[GC]), x = NULL) +
      guides(colour = guide_legend(nrow = 1, keywidth = 0.8, keyheight = 1),
             fill = guide_legend(nrow = 1, keywidth = 0.8, keyheight = 1)) +
      theme_classic(base_size = 11) +
      themes +
      theme(
        axis.title = element_text(size = 9, face = "bold"),
        axis.text = element_text(size = 10, color = "grey30"),
        axis.ticks = element_line(color = "grey70", linewidth = 0.3),
        axis.line = element_line(color = "grey40", linewidth = 0.4),
        panel.grid.major.y = element_line(color = "grey92", linewidth = 0.3),
        legend.text = element_text(size = 8, color = "grey30"),
        legend.position = "bottom",
        legend.spacing.x = unit(2, "pt"),
        legend.spacing.y = unit(0, "pt"),
        legend.box.margin = margin(0, 0, 0, 0),
        legend.key.size = unit(8, "pt"),
        plot.margin = unit(c(0.8, 0.3, 0, 0), "cm")
      )
    plt_lst <- c(plt_lst, list(plt_name = p_qc))
  }
  return(plt_lst)
}

TYPES <- c('genome_variant', 'exome_variant', 'exome_gene_0.001')
# TYPES <- c('genome_variant')
INDEX <- c('genome_variant'='21', 'exome_variant'='22', 'exome_gene_0.001'='20')
for(type in TYPES){
  plt_lst <- plot_lambda_gc_boxplot(test_type = type)
  
  p = ggpubr::ggarrange(plotlist = plt_lst, vjust = 1.2, hjust = 0, 
                                  ncol = 2, nrow = 3, common.legend = TRUE, 
                                  labels = c('AFR-like', 'AMR-like', 'EAS-like', 'EUR-like', 'MID-like', 'SAS-like'), font.label = list(size = 12 ))
  output_figure(p, 'supplement', paste0('supp_figure', INDEX[type] , '_aou_', type,'_lambda_gc'), 6,7.5)

}


##########################
### Archived: lambda GC ##
##########################
gene_data <- rbind(
  read_gcs_fast(paste0('gs://aou_analysis/v8/data/qc/aou_afr_exome_gene_0.001_phenotype_lambda_gc.txt.bgz')) %>% select(1:13),
  read_gcs_fast(paste0('gs://aou_analysis/v8/data/qc/aou_amr_exome_gene_0.001_phenotype_lambda_gc.txt.bgz')) %>% select(1:13),
  read_gcs_fast(paste0('gs://aou_analysis/v8/data/qc/aou_eas_exome_gene_0.001_phenotype_lambda_gc.txt.bgz')) %>% select(1:13),
  read_gcs_fast(paste0('gs://aou_analysis/v8/data/qc/aou_eur_exome_gene_0.001_phenotype_lambda_gc.txt.bgz')) %>% select(1:13),
  read_gcs_fast(paste0('gs://aou_analysis/v8/data/qc/aou_mid_exome_gene_0.001_phenotype_lambda_gc.txt.bgz')) %>% select(1:13),
  read_gcs_fast(paste0('gs://aou_analysis/v8/data/qc/aou_sas_exome_gene_0.001_phenotype_lambda_gc.txt.bgz')) %>% select(1:13)
)
meta_gene_data <- read_gcs_fast(paste0('gs://aou_analysis/v8/data/qc/aou_meta_exome_gene_0.001_phenotype_lambda_gc.txt.bgz'))

genome_data <- rbind(
  read_gcs_fast(paste0('gs://aou_analysis/v8/data/qc/aou_afr_genome_variant_phenotype_lambda_gc.txt.bgz')),
  read_gcs_fast(paste0('gs://aou_analysis/v8/data/qc/aou_amr_genome_variant_phenotype_lambda_gc.txt.bgz')),
  read_gcs_fast(paste0('gs://aou_analysis/v8/data/qc/aou_eas_genome_variant_phenotype_lambda_gc.txt.bgz')),
  read_gcs_fast(paste0('gs://aou_analysis/v8/data/qc/aou_eur_genome_variant_phenotype_lambda_gc.txt.bgz')),
  read_gcs_fast(paste0('gs://aou_analysis/v8/data/qc/aou_mid_genome_variant_phenotype_lambda_gc.txt.bgz')),
  read_gcs_fast(paste0('gs://aou_analysis/v8/data/qc/aou_sas_genome_variant_phenotype_lambda_gc.txt.bgz'))
)
meta_genome_data <- read_gcs_fast(paste0('gs://aou_analysis/v8/data/qc/aou_meta_genome_variant_phenotype_lambda_gc.txt.bgz'))
meta_genome_data <- read_gcs_fast(paste0('gs://aou_analysis/v8/data/qc/aou_meta_genome_variant_phenotype_lambda_gc_filtered.txt.bgz'))
meta_exome_data <- read_gcs_fast(paste0('gs://aou_analysis/v8/data/qc/aou_meta_exome_variant_phenotype_lambda_gc.txt.bgz'))
meta_exome_data <- read_gcs_fast(paste0('gs://aou_analysis/v8/data/qc/aou_meta_exome_variant_phenotype_lambda_gc_filtered.txt.bgz'))
pheno_info <- read_tsv('~/Dropbox (Partners HealthCare)/aou_v8/data/phenotype/aou_v8_phenotype_info_qced_with_meta.tsv')  %>% 
  filter(ancestry != 'META' & pheno_sex == 'both') %>%
  select(ancestry, phenoname) %>%
  merge(., meta_genome_data, by = 'phenoname') %>%
  group_by(phenoname, category, n_cases, n_controls, lambda_gc_raw, lambda_gc_af_0.001, lambda_gc_af_0.0001) %>%
  dplyr::summarize(ancestries = paste(ancestry.x, collapse = ','))
sub_pheno_info_upper <- pheno_info %>%
  filter(lambda_gc_af_0.001 < 0.8 & !startsWith(phenoname, 'random'))
sub_pheno_info_lower <- pheno_info %>%
  filter(lambda_gc_af_0.001 >= 0.8 & !startsWith(phenoname, 'random'))
table(sub_pheno_info_upper$ancestries)
# table(sub_pheno_info_lower$ancestries)
sum(sub_pheno_info_lower$ancestries == 'EUR')

sub_pheno_info_upper <- pheno_info %>%
  filter(lambda_gc_af_0.0001 < 0.6 & !startsWith(phenoname, 'random'))
sub_pheno_info_lower <- pheno_info %>%
  filter(lambda_gc_af_0.0001 >= 0.6 & !startsWith(phenoname, 'random'))
table(sub_pheno_info_upper$ancestries)
table(sub_pheno_info_lower$ancestries)
sum(sub_pheno_info_lower$ancestries == 'EUR')

exome_data <- rbind(
  read_gcs_fast(paste0('gs://aou_analysis/v8/data/qc/aou_afr_exome_variant_phenotype_lambda_gc.txt.bgz')),
  read_gcs_fast(paste0('gs://aou_analysis/v8/data/qc/aou_amr_exome_variant_phenotype_lambda_gc.txt.bgz')),
  read_gcs_fast(paste0('gs://aou_analysis/v8/data/qc/aou_eas_exome_variant_phenotype_lambda_gc.txt.bgz')),
  read_gcs_fast(paste0('gs://aou_analysis/v8/data/qc/aou_eur_exome_variant_phenotype_lambda_gc.txt.bgz')),
  #read_gcs_fast(paste0('gs://aou_analysis/v8/data/qc/aou_mid_exome_variant_phenotype_lambda_gc.txt.bgz')),
  read_gcs_fast(paste0('gs://aou_analysis/v8/data/qc/aou_sas_exome_variant_phenotype_lambda_gc.txt.bgz'))
)


for(pop in c('meta', 'afr', 'amr', 'eas', 'eur', 'mid', 'sas')){
  # plot_lambda(pop_group=pop, type='genome_variant', save=T)
  # plot_lambda(pop_group=pop, type='exome_variant', save=T)
  plot_lambda(pop_group=pop, type='exome_gene_0.001', save=T)
}


for(pop in c('meta', 'afr', 'amr', 'eas', 'eur', 'mid', 'sas')){
  # plot_lambda(pop_group=pop, type='genome_variant', save=T, qc=T, filter=T)
  # plot_lambda(pop_group=pop, type='exome_variant', save=T, qc=T, filter=T)
  plot_lambda(pop_group=pop, type='exome_gene_0.001', save=T, qc=T, filter=T)
}

plot_lambda <- function(pop_group, type, save=T, filter = F, qc = F, output_root = '~/Desktop/'){
  plt_name <- paste0(pop_group, '_lambda_gc_', type, if_else(filter, '_filtered', ''), if_else(qc, '_QCed', ''))
  pheno_info_full <- read_tsv('~/Dropbox (Partners HealthCare)/aou_v8/data/phenotype/aou_v8_phenotype_info_qced_with_meta.tsv') %>%
    filter(ancestry == toupper(pop_group) & pheno_sex == 'both') %>%
    select(phenoname, category, description, n_cases, n_controls) %>%
    distinct()
  
  if(pop_group != 'meta'| filter){
    pheno_info_full <- pheno_info_full %>% select(-n_cases, -n_controls)
  }
  if(!grepl('gene', type)){
    pheno_info_full <- pheno_info_full %>% select(-category, -description)
  }
  
  print(plt_name)
  data <- read_gcs_fast(paste0('gs://aou_analysis/v8/data/qc/aou_', pop_group,'_', type, '_phenotype_lambda_gc', if_else(filter, '_filtered', ''),'.txt.bgz')) %>%
    merge(., pheno_info_full, by = 'phenoname')%>%
    filter(category != 'random_phenotype') %>%
    mutate(category = factor(category, levels = c('lab_measurement', 'physical_measurement', 'mcc2_phecode', 
                                                  'mcc2_phecodex', 'r_drug', 'pfhh_survey', 'onset', 'progression', 'random_phenotype', 'mhwb_survey'),
                             labels = c('Lab measurements', 'Physical measurements', 'Phecode', 
                                        'PhecodeX', 'Prescriptions', 'Personal family health history (PFHH)',
                                        'Onset', 'Progression', 'Random phenotype', 'Mental health and Well-being (MHWB)'))) %>%
    filter(!phenoname %in% c('PP_932.11', 'PP_928.1', '3011099', '3046664', '3013861', '3052648', '3045792_3016921'))
  
  if(filter){
    data <- data %>% mutate(lambda_gc = lambda_gc_filtered)
  }else{
    data <- data %>% mutate(lambda_gc = lambda_gc_raw)
  }
  
  
  print('Plotting overall lambda...')
  data0 <- data
  if(filter){
    data0 <- data %>%
      filter(lambda_gc < 2 & lambda_gc > 0 & n_cases > 200)
  }
  p0 <- data0 %>%
    ggplot + aes(x = n_cases, y = lambda_gc, color = category) +
    labs(x = 'N cases', y = expression(lambda[GC])) +
    geom_point() +
    themes + scale_x_log10(label=comma) + geom_hline(yintercept = 1, lty =2) +
    scale_color_manual(name = 'Category', values = category_colors, breaks = names(category_colors), labels = category_labels) +
    scale_fill_manual(name = 'Category', values = category_colors, breaks = names(category_colors), labels = category_labels)  +
    guides(colour = guide_legend(nrow = 3)) + theme(plot.margin = unit(c(0.5,0,0,0), "cm")) 
  if(save){
    filter_tag = if_else(filter, '_filtered', '')
    output_figure(p0, 'qc', paste0(plt_name, '_overall'), 4, 6)
  }
  
  print('Plotting overall lambda with geom_smooth...')
  data1 <- data
  if(filter){
    data1 <- data %>%
      filter(lambda_gc < 2 & lambda_gc > 0 & n_cases > 200)
  }
  if(grepl( 'gene', type)){
    data1 <- data1 %>%
      mutate(lambda_gc = lambda_gc_annotation_synonymous)
  }
  p1 <- data1 %>%
    filter(category %in% c('PhecodeX', 'Prescriptions', 'Personal family health history (PFHH)', 'Mental health and Well-being (MHWB)')) %>%
    ggplot + aes(x = n_cases, y = lambda_gc, color = category) +
    labs(x = 'N cases', y = expression(lambda[GC])) +
    geom_point(alpha = 0.25) +
    geom_smooth() +
    themes + scale_x_log10(label=comma) + geom_hline(yintercept = 1, lty =2) +
    scale_color_manual(name = 'Category', values = category_colors, breaks = names(category_colors), labels = category_labels) +
    scale_fill_manual(name = 'Category', values = category_colors, breaks = names(category_colors), labels = category_labels)  +
    guides(colour = guide_legend(nrow = 2)) + theme(plot.margin = unit(c(0.5,0,0,0), "cm"))
  if(save){
    filter_tag = if_else(filter, '_filtered', '')
    output_figure(p1, 'qc', paste0(plt_name, '_geom_smooth', if_else(grepl('gene', type), '_synonymous', '')), 4, 6)
  }
  
  print('Plotting AF binned lambda...')
  data2 <- data %>%
    pivot_longer(cols = starts_with(if_else(grepl('gene', type), 'lambda_gc_caf_', 'lambda_gc_af_')), names_to = 'AF_interval', values_to = 'Lambda')
  if(filter){
    data2 <- data2 %>%
      filter(Lambda < 2 & Lambda > 0 & n_cases > 200)
  }
  p2 <- data2 %>%
    mutate(AF_interval = factor(as.numeric(str_replace(AF_interval, if_else(grepl('gene', type), 'lambda_gc_caf_', 'lambda_gc_af_'), '')),
                                levels = c(0, 0.0001, 0.001, 0.01, 0.1),
                                labels = c('[0%, 0.01%]', '(0.01%, 0.1%]', '(0.1%, 1%]', '(1%, 10%]', paste0('(10%, ', bquote("\U221E"), ' )'))))%>%
    ggplot + aes(x = n_cases, y = Lambda, color = category) +
    labs(x = 'N cases', y = expression(lambda[GC])) +
    geom_point() + themes + scale_x_log10(label=comma) + geom_hline(yintercept = 1, lty =2) +
    scale_color_manual(name = 'Category', values = category_colors, breaks = names(category_colors), labels = category_labels) +
    scale_fill_manual(name = 'Category', values = category_colors, breaks = names(category_colors), labels = category_labels) + 
    facet_grid(~AF_interval) + guides(color = guide_legend(nrow = 2))
  if(save){
    filter_tag = if_else(filter, '_filtered', '')
    output_figure(p2, 'qc', paste0(plt_name, '_af_binned'), 3, 12)
  }
  
  print('Plotting coverage binned lambda...')
  data3 <- data %>%
    pivot_longer(cols = starts_with('lambda_gc_coverage_'), names_to = 'coverage_interval', values_to = 'Lambda')
  if(filter){
    data3 <- data3 %>%
      filter(Lambda < 2 & Lambda > 0 & n_cases > 200)
  }
  p3 <- data3 %>%
    mutate(coverage_interval = factor(as.numeric(str_replace(coverage_interval, 'lambda_gc_coverage_', '')),
                                      levels = c(0, 10, 20, 30),
                                      labels = c('[0, 10]', '(10, 20]', '(20, 30]', paste0('(30, ', bquote("\U221E"), ' )'))))%>%
    ggplot + aes(x = n_cases, y = Lambda, color = category) +
    labs(x = 'N cases', y = expression(lambda[GC])) +
    geom_point() + themes + scale_x_log10(label=comma) + geom_hline(yintercept = 1, lty =2) +
    scale_color_manual(name = 'Category', values = category_colors, breaks = names(category_colors), labels = category_labels) +
    scale_fill_manual(name = 'Category', values = category_colors, breaks = names(category_colors), labels = category_labels) +
    facet_grid(~coverage_interval) + guides(color = guide_legend(nrow = 2))
  
  if(save){
    filter_tag = if_else(filter, '_filtered', '')
    output_figure(p3, 'qc', paste0(plt_name, '_coverage_binned'), 3, 12)
  }
  if(grepl( 'gene', type)){
    print('Plotting lambda by annotation...')
    data4 <- data %>%
      pivot_longer(cols = starts_with('lambda_gc_annotation_'), names_to = 'annotation', values_to = 'Lambda')
    if(filter){
      data4 <- data4 %>%
        filter(Lambda < 2 & Lambda > 0 & n_cases > 200)
    }
    p4 <- data4 %>%
      mutate(annotation = str_replace(annotation, 'lambda_gc_annotation_', ''))%>%
      mutate(annotation = if_else(annotation == 'missenseLC', 'missense|LC', annotation)) %>%
      mutate(annotation = factor(if_else(annotation == 'pLoF;missenseLC', 'pLoF|missense|LC', annotation), levels=annotation_types2)) %>%
      ggplot + aes(x = n_cases, y = Lambda, color = category) +
      labs(x = 'N cases', y = expression(lambda[GC])) +
      geom_point() + themes + scale_x_log10(label=comma) + geom_hline(yintercept = 1, lty =2) +
      scale_color_manual(name = 'Category', values = category_colors, breaks = names(category_colors), labels = category_labels) +
      scale_fill_manual(name = 'Category', values = category_colors, breaks = names(category_colors), labels = category_labels) +
      facet_grid(~annotation, labeller = labeller(annotation=annotation_names2)) + 
      guides(color = guide_legend(nrow = 2))
    
    if(save){
      filter_tag = if_else(filter, '_filtered', '')
      qc_tag = if_else(qc, '_QCed', '')
      output_figure(p4, 'qc', paste0(plt_name, '_annotation'), 4, 8)
    }
    
  }
  
  if(qc){
    p_qc <- data %>%
      filter(lambda_gc_raw < 2) %>%
      select(lambda_gc_raw, lambda_gc_filtered, category, n_cases) %>%
      pivot_longer(cols = starts_with('lambda_gc_'), names_to = 'QC', values_to = 'Lambda') %>% 
      mutate(QC = factor(if_else(QC == 'lambda_gc_filtered', 'QCed', 'Raw'), levels = c('Raw', 'QCed'))) %>%
      ggplot + aes(x = QC, y = Lambda, color = category, group = interaction(QC, category)) +
      labs(y = expression(lambda[GC]), x = NULL) +
      ylim(NA, 1.25) + 
      geom_boxplot() + themes + geom_hline(yintercept = 1, lty =2) +
      scale_color_manual(name = 'Category', values = category_colors, breaks = names(category_colors), labels = category_labels) +
      scale_fill_manual(name = 'Category', values = category_colors, breaks = names(category_colors), labels = category_labels) +
      guides(colour = guide_legend(nrow = 2))+ theme(plot.margin = unit(c(0.5,0,0,0), "cm"))
    if(save){
      filter_tag = if_else(filter, '_filtered', '')
      qc_tag = if_else(qc, '_QCed', '')
      output_figure(p_qc, 'qc', paste0(plt_name, '_QCed_overall', if_else(grepl('gene', type), '_synonymous', '')), 4, 6)
    }
    
  }
  
  print('Plotting lambda geom_hex ...')
  data5 <- data
  if(filter){
    data5 <- data %>%
      filter(lambda_gc < 2 & lambda_gc > 0 & n_cases > 200)
  }
  if(grepl( 'gene', type)){
    data5 <- data5 %>%
      mutate(lambda_gc = lambda_gc_annotation_synonymous)
  }
  p5 <- data5 %>%
    ggplot + aes(x = n_cases, y = lambda_gc) +
    labs(x = 'N cases', y = expression(lambda[GC])) +
    geom_point(alpha = 0.25) +
    geom_hex() +
    themes + scale_x_log10(label=comma) + geom_hline(yintercept = 1, lty =2) +
    guides(colour = guide_legend(nrow = 2)) + theme(plot.margin = unit(c(0.5,0,0,0), "cm"))
  if(save){
    filter_tag = if_else(filter, '_filtered', '')
    output_figure(p5, 'qc', paste0(plt_name, '_geom_hex', if_else(grepl('gene', type), '_synonymous', '')), 4, 6)
  }
  
  if(grepl( 'gene', type)){
    lst = list(p0, p1, p2, p3, p5, p4)
  }else{
    lst = list(p0, p1, p2, p3, p5)
  }
  if(qc){
    lst[['qc']] = p_qc
  }
  
  return(lst)
}

arrange_pop_plots <- function(index, output_root, qc=F, filter=F){
  plt_names <- c('1' = 'main', '2' = 'geom_smooth', '3' = 'af_bin', '4'= 'coverage_bin', '5' = 'geom_hex', 
                 '6' = 'annotations', 'qc' = 'qc'
  )
  TYPES <- c('genome_variant', 'exome_variant', 'exome_gene_0.001')
  # TYPES <- c('exome_gene_0.001')
  output_name <- plt_names[as.character(index)]
  plt_lst <- list()
  for(type in TYPES){
    for(pop in c('afr', 'amr', 'eas', 'eur', 'mid', 'sas')){
      for(filter in c(T)){
        plt_name <- paste0(pop, '_', type, if_else(filter, '_filtered', ''),  if_else(qc, '_QCed', ''))
        print(plt_name)
        if(index == 6 & grepl('variant', type)) next 
        plts <- plot_lambda(pop_group=pop, type=type, save=F, filter = filter, qc=qc)[[index]]
        if(filter){
          plt_lst <- c(plt_lst, list(plt_name = plts))
        }
      }
    }
  }
  if(index %in% c(1,2,5)){
    genome_lambda = ggpubr::ggarrange(plotlist = plt_lst[1:6], vjust = 1.2, hjust = 0,
                                      ncol = 3, nrow=2, common.legend = TRUE, labels = toupper(c('afr', 'amr', 'eas', 'eur', 'mid', 'sas')), font.label = list(size = 12 ))
    output_figure(genome_lambda, 'qc', paste0('aou_genome_lambda_all_pop_', output_name, if_else(filter, '_filtered', ''),  if_else(qc, '_QCed', '')), 4, 9)

    exome_lambda = ggpubr::ggarrange(plotlist = plt_lst[7:12], vjust = 1.2, hjust = 0,
                                     ncol = 3, nrow=2, common.legend = TRUE, labels = toupper(c('afr', 'amr', 'eas', 'eur', 'mid', 'sas')), font.label = list(size = 12 ))
    output_figure(exome_lambda, 'qc', paste0('aou_exome_lambda_all_pop_', output_name, if_else(filter, '_filtered', ''),  if_else(qc, '_QCed', '')), 4,9)
    
    annt_lambda = ggpubr::ggarrange(plotlist = plt_lst[13:18], vjust = 1.2, hjust = 0, 
                                    ncol = 3, nrow=2, common.legend = TRUE, labels = toupper(c('afr', 'amr', 'eas', 'eur', 'mid', 'sas')), font.label = list(size = 12 ))
    output_figure(annt_lambda, 'qc', paste0('aou_gene_lambda_all_pop_', output_name, if_else(filter, '_filtered', ''),  if_else(qc, '_QCed', '')), 4, 9)
  }else if(index %in% c(3,4)){
    genome_lambda = ggpubr::ggarrange(plotlist = plt_lst[1:6], vjust = 1.2, hjust = 0,
                                      ncol = 2, nrow = 3, common.legend = TRUE, labels = toupper(c('afr', 'amr', 'eas', 'eur', 'mid', 'sas')), font.label = list(size = 12 ))
    output_figure(genome_lambda, 'qc', paste0('aou_genome_lambda_all_pop_', output_name, if_else(filter, '_filtered', ''),  if_else(qc, '_QCed', '')), 8, 15)

    exome_lambda = ggpubr::ggarrange(plotlist = plt_lst[7:12], vjust = 1.2, hjust = 0,
                                     ncol = 2, nrow = 3, common.legend = TRUE, labels = toupper(c('afr', 'amr', 'eas', 'eur', 'mid', 'sas')), font.label = list(size = 12 ))
    output_figure(exome_lambda, 'qc', paste0('aou_exome_lambda_all_pop_', output_name, if_else(filter, '_filtered', ''),  if_else(qc, '_QCed', '')), 8,15)
    
    annt_lambda = ggpubr::ggarrange(plotlist = plt_lst[13:18], vjust = 1.2, hjust = 0, 
                                    ncol = 2, nrow = 3, common.legend = TRUE, labels = toupper(c('afr', 'amr', 'eas', 'eur', 'mid', 'sas')), font.label = list(size = 12 ))
    output_figure(annt_lambda, 'qc', paste0('aou_gene_lambda_all_pop_', output_name, if_else(filter, '_filtered', ''),  if_else(qc, '_QCed', '')), 8,15)
    
  }else if(index == 6){
    
    annt_lambda = ggpubr::ggarrange(plotlist = plt_lst, vjust = 1.2, hjust = 0, 
                                    ncol = 2, nrow = 3, common.legend = TRUE, 
                                    labels = toupper(c('afr', 'amr', 'eas', 'eur', 'mid', 'sas')), font.label = list(size = 12 ))
    output_figure(annt_lambda, 'qc', paste0('aou_gene_lambda_all_pop_', output_name, if_else(filter, '_filtered', ''),  if_else(qc, '_QCed', '')), 8,15)
    
  }else{
    genome_lambda = ggpubr::ggarrange(plotlist = plt_lst[1:6], vjust = 1.2, hjust = 0,
                                      ncol = 2, nrow = 3, common.legend = TRUE, labels = toupper(c('afr', 'amr', 'eas', 'eur', 'mid', 'sas')), font.label = list(size = 12 ))
    output_figure(genome_lambda, 'qc', paste0('aou_genome_lambda_all_pop_', output_name, if_else(filter, '_filtered', ''),  if_else(qc, '_QCed', '')), 8, 15)
    
    exome_lambda = ggpubr::ggarrange(plotlist = plt_lst[7:12], vjust = 1.2, hjust = 0,
                                     ncol = 2, nrow = 3, common.legend = TRUE, labels = toupper(c('afr', 'amr', 'eas', 'eur', 'mid', 'sas')), font.label = list(size = 12 ))
    output_figure(exome_lambda, 'qc', paste0('aou_exome_lambda_all_pop_', output_name, if_else(filter, '_filtered', ''),  if_else(qc, '_QCed', '')), 8,15)
    
    annt_lambda = ggpubr::ggarrange(plotlist = plt_lst[13:18], vjust = 1.2, hjust = 0, 
                                    ncol = 2, nrow = 3, common.legend = TRUE, 
                                    labels = toupper(c('afr', 'amr', 'eas', 'eur', 'mid', 'sas')), font.label = list(size = 12 ))
    output_figure(annt_lambda, 'qc', paste0('aou_gene_lambda_all_pop_', output_name, if_else(filter, '_filtered', ''),  if_else(qc, '_QCed', '')), 8,15)
  }
}


arrange_pop_plots(1)
arrange_pop_plots(2)
arrange_pop_plots(3)
arrange_pop_plots(4)
arrange_pop_plots(5)
arrange_pop_plots(6)


arrange_pop_plots(1, qc=T, filter = T)
arrange_pop_plots(2, qc=T, filter = T)
arrange_pop_plots(3, qc=T, filter = T)
arrange_pop_plots(4, qc=T, filter = T)
arrange_pop_plots(5, qc=T, filter = T)
arrange_pop_plots(6, qc=T, filter = T)
arrange_pop_plots('qc', qc=T, filter = T)


genewise_lambda_gc <- rbind(
  read_gcs_fast('gs://aou_analysis/v8/data/qc/aou_afr_gene_lambda_gc.txt.bgz') %>% mutate(ancestry = 'afr'),
  read_gcs_fast('gs://aou_analysis/v8/data/qc/aou_amr_gene_lambda_gc.txt.bgz')  %>% mutate(ancestry = 'amr'),
  read_gcs_fast('gs://aou_analysis/v8/data/qc/aou_eas_gene_lambda_gc.txt.bgz')  %>% mutate(ancestry = 'eas'),
  read_gcs_fast('gs://aou_analysis/v8/data/qc/aou_eur_gene_lambda_gc.txt.bgz')  %>% mutate(ancestry = 'eur'),
  read_gcs_fast('gs://aou_analysis/v8/data/qc/aou_mid_gene_lambda_gc.txt.bgz')  %>% mutate(ancestry = 'mid'),
  read_gcs_fast('gs://aou_analysis/v8/data/qc/aou_sas_gene_lambda_gc.txt.bgz') %>% mutate(ancestry = 'sas'),
  read_gcs_fast('gs://aou_analysis/v8/data/qc/aou_meta_gene_lambda_gc.txt.bgz') %>% mutate(ancestry = 'meta')
)

genewise_lambda_gc %>%
  mutate(annotation = factor(if_else(annotation== 'missenseLC', 'missense|LC', annotation), levels = annotation_types),
         ancestry = factor(ancestry, levels = names(pops), labels = pops)) %>%
  filter(annotation != 'Cauchy') %>%
  ggplot + aes(x = ancestry, y = lambda_gc_burden_gene, color = annotation, group = interaction(ancestry, annotation)) +
  geom_boxplot() + themes + labs(x = NULL, y = expression(lambda[GC][burden])) + annotation_color_scale + 
  ylim(0, 20) + theme(axis.text.x = element_text(angle=90, hjust =1))


genewise_lambda_gc_filtered <- rbind(
  read_gcs_fast('gs://aou_analysis/v8/data/qc/aou_afr_gene_lambda_gc.txt.bgz') %>% filter(CAF > 1e-4)  %>% mutate(ancestry = 'afr'),
  read_gcs_fast('gs://aou_analysis/v8/data/qc/aou_amr_gene_lambda_gc.txt.bgz') %>% filter(`CAF` > 1e-4)%>% mutate(ancestry = 'amr'),
  read_gcs_fast('gs://aou_analysis/v8/data/qc/aou_eas_gene_lambda_gc.txt.bgz') %>% filter(`CAF` > 1e-4)  %>% mutate(ancestry = 'eas'),
  read_gcs_fast('gs://aou_analysis/v8/data/qc/aou_eur_gene_lambda_gc.txt.bgz') %>% filter(`CAF` > 1e-4) %>% mutate(ancestry = 'eur'),
  read_gcs_fast('gs://aou_analysis/v8/data/qc/aou_mid_gene_lambda_gc.txt.bgz') %>% filter(`CAF` > 1e-5)  %>% mutate(ancestry = 'mid'),
  read_gcs_fast('gs://aou_analysis/v8/data/qc/aou_sas_gene_lambda_gc.txt.bgz') %>% filter(`CAF` > 1e-5) %>% mutate(ancestry = 'sas'),
  read_gcs_fast('gs://aou_analysis/v8/data/qc/aou_meta_gene_lambda_gc.txt.bgz') %>% filter(`CAF` > 1e-4) %>% mutate(ancestry = 'meta')
)

genewise_lambda_gc_filtered %>%
  mutate(annotation = factor(if_else(annotation== 'missenseLC', 'missense|LC', annotation), levels = annotation_types),
         ancestry = factor(ancestry, levels = names(pops), labels = pops)) %>%
  filter(annotation != 'Cauchy') %>%
  filter(total_variants > 10 & mean_coverage > 10) %>%
  ggplot + aes(x = ancestry, y = lambda_gc_burden_gene, color = annotation, group = interaction(ancestry, annotation)) +
  geom_boxplot() + themes + labs(x = NULL, y = expression(lambda[GC][burden])) + annotation_color_scale + 
  ylim(0, 20) + theme(axis.text.x = element_text(angle=90, hjust =1))


x <- read_gcs_fast('gs://aou_analysis/v8/data/qc/aou_meta_gene_lambda_gc.txt.bgz')
y <- read_gcs_fast('gs://aou_analysis/v8/data/qc/aou_sas_gene_lambda_gc.txt.bgz')
table(y$hq_gene)
table(y$hq_gene_lambda)
sum(!is.na(y$lambda_gc_burden_gene_filtered))


x %>%
  mutate(annotation = factor(if_else(annotation== 'missenseLC', 'missense|LC', annotation), levels = annotation_types)) %>%
  filter(annotation != 'Cauchy') %>%
  ggplot + aes(x =annotation, y = lambda_gc_burden_gene, color = annotation, group = annotation) +
  geom_boxplot() + themes + labs(x = NULL, y = expression(lambda[GC][burden])) + annotation_color_scale + 
  theme(axis.text.x = element_text(angle=90, hjust =1)) + scale_y_log10(label = comma)


x %>%
  mutate(annotation = factor(if_else(annotation== 'missenseLC', 'missense|LC', annotation), levels = annotation_types)) %>%
  filter(annotation != 'Cauchy') %>%
  filter(total_variants > 10 & mean_coverage > 10 & CAF > 0.0001) %>%
  ggplot + aes(x =annotation, y = lambda_gc_burden_gene, color = annotation, group = annotation) +
  geom_boxplot() + themes + labs(x = NULL, y = expression(lambda[GC][burden])) + annotation_color_scale + 
  theme(axis.text.x = element_text(angle=90, hjust =1)) + scale_y_log10(label = comma)

x %>%
  mutate(annotation = factor(if_else(annotation== 'missenseLC', 'missense|LC', annotation), levels = annotation_types)) %>%
  filter(annotation != 'Cauchy') %>%
  filter(total_variants > 10 & mean_coverage > 10 & CAF > 0.0001) %>%
  ggplot + aes(x =annotation, y = lambda_gc_burden_gene_filtered, color = annotation, group = annotation) +
  geom_boxplot() + themes + labs(x = NULL, y = expression(lambda[GC][burden])) + annotation_color_scale + 
  theme(axis.text.x = element_text(angle=90, hjust =1)) + scale_y_log10(label = comma)



plot_lambda_comparison <- function(pop_group, type1, type2, save, filter, output_root = '~/Desktop/'){
  file_names <- c('gene'= 'exome_gene_0.001', 'exome'='exome_variant', 'genome'='genome_variant')
  
  lambda1 <- read_gcs_fast(paste0('gs://aou_analysis/v8/data/qc/aou_', pop_group, '_', file_names[type1],'_phenotype_lambda_gc_filtered.txt.bgz'))%>%
    select(-n_cases, -n_controls) %>%
    mutate(lambda_gc = lambda_gc_filtered)
  lambda2 <-  read_gcs_fast(paste0('gs://aou_analysis/v8/data/qc/aou_', pop_group, '_', file_names[type2],'_phenotype_lambda_gc_filtered.txt.bgz')) %>%
    select(-n_cases, -n_controls)%>%
    mutate(lambda_gc = lambda_gc_filtered)
  if(type2 == 'gene'){
    lambda2 <- lambda2 %>%
      mutate(lambda_gc = lambda_gc_annotation_synonymous)
  }
  lambda <- lambda1 %>%
    merge(., lambda2, by = 'phenoname') %>%
    mutate(ancestry = toupper(pop_group)) %>%
    merge(., pheno_info_full %>% filter(pheno_sex == 'both'), by = c('ancestry', 'phenoname'), all.x = T) 
  if(filter){
    lambda <- lambda %>%
      filter(lambda_gc.x < 10 & lambda_gc.y < 10)
  }
  
  p <- lambda %>%
    ggplot + aes(x = lambda_gc.x, y = lambda_gc.y, color = category) +
    labs(x = paste(str_to_sentence(type1), 'Lambda GC'), y = paste(str_to_sentence(type2), 'Lambda GC', if_else(type2 == 'gene', '\n(synonymous)', ''))) +
    geom_vline(xintercept = 1, lty =2) +
    ylim(0.9, 1.15) +
    geom_point() + themes+ geom_hline(yintercept = 1, lty =2) +
    scale_color_manual(name = 'Category', values = category_colors, breaks = names(category_colors), labels = category_labels) +
    scale_fill_manual(name = 'Category', values = category_colors, breaks = names(category_colors), labels = category_labels) +
    guides(colour = guide_legend(nrow = 1)) + theme(plot.margin = unit(c(0.5,0.5,0,0.5), "cm"))
  if(save){
    filter_tag = if_else(filter, '_filtered', '')
    png(paste0(output_root,'aou_', pop_group,'_', file_names[type1],'_', file_names[type2],'_lambda_comparison', filter_tag,'.png'), height = 4, width = 6, units = 'in', res = 300)
    print(p)
    dev.off()
  }
  return(p)
}

plt_lst <- list()
for(pop in c('afr', 'amr', 'eas', 'eur', 'mid', 'sas')){
  for(filter in c(T, F)){
    plt_name <- paste0(pop, if_else(filter, '_filtered', ''))
    print(plt_name)
    plt <- plot_lambda_comparison(pop_group=pop, type1 = 'genome', type2='exome', save=F, filter = filter, output_root = aou_qc_fig_path)
    if(filter){
      plt_lst <- c(plt_lst, list(plt))
    }
  }
}


lambda_comparison = ggpubr::ggarrange(plotlist = plt_lst,
                                      ncol = 3, nrow=2, common.legend = TRUE, labels = toupper(c('afr', 'amr', 'eas', 'eur', 'mid', 'sas')), font.label = list(size = 8 ))
output_figure(lambda_comparison, 'qc', paste0('aou_lambda_comparison_genome_exome'), 4,9)

plt_lst <- list()
for(pop in c('afr', 'amr', 'eas', 'eur', 'mid', 'sas')){
  for(filter in c(T, F)){
    plt_name <- paste0(pop, if_else(filter, '_filtered', ''))
    print(plt_name)
    plt <- plot_lambda_comparison(pop_group=pop, type1 = 'genome', type2='exome', save=F, filter = filter, output_root = aou_qc_fig_path)
    if(filter){
      plt_lst <- c(plt_lst, list(plt))
    }
  }
}


lambda_comparison = ggpubr::ggarrange(plotlist = plt_lst,
                                      ncol = 3, nrow=2, common.legend = TRUE, labels = toupper(c('afr', 'amr', 'eas', 'eur', 'mid', 'sas')), font.label = list(size = 8 ))
png(paste0('~/Desktop/aou_lambda_comparison_acaf_gene_synonymous_all_pop_filtered_expected_ac_5.png'), width=9, height=4, units = 'in', res = 300)
print(lambda_comparison)
dev.off()

