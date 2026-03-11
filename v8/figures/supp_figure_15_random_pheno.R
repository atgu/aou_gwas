source('~/Dropbox (Partners HealthCare)/github_repo/aou_gwas/R/constants.R')

raw_qq_data <- function(max_MAF = 0.001, type = 'burden'){
  p_field = P_VALUE_FIELDS[type]
  print(p_field)
  qq_data <- data.frame()
  for(anc in names(pops)){
    if(anc == 'meta') next
    for(prevalence in c('1e-04', '0.001', '0.01', '0.02', '0.05', '0.1', '0.2', '0.5', 'continuous')){
      print(paste0('Loading ', anc, ' random ', prevalence, ' ...'))
      tmp_qq <- read_gcs_fast(paste0('gs://aou_analysis/v8/ht_results/', toupper(anc),'/phenotype_random_0.5_', prevalence,'_3/gene_results_exact_expected_p_', max_MAF,'.txt.bgz')) %>%
        mutate(ancestry = anc, 
               prevalence = if_else(prevalence == '1e-04', '0.0001', prevalence))
      tmp_qq$observed = tmp_qq[[paste0(p_field, '_log10')]]
      tmp_qq$expected = tmp_qq[[paste0(p_field, '_expected_log10')]]
      tmp_qq <- tmp_qq %>% select(ancestry, prevalence, observed, expected)
      qq_data <- rbind(qq_data, tmp_qq)
    }
  }
  qq_data <- qq_data %>%
    mutate(
      ancestry = factor(ancestry, levels = c( 'afr', 'amr', 'eas', 'eur', 'mid','sas')), 
      prevalence = factor(prevalence, levels = c('0.0001', '0.001', '0.01', '0.02', '0.05', '0.1', '0.2', '0.5', 'continuous'))
    ) 
  return(qq_data)
}


prepare_qq_plot_data_v8 <- function(max_MAF = '0.001', type='burden'){
  ## Generated in quick_analysis.py
  lambda_data <- read_gcs_fast(paste0('gs://aou_wlu/v8_analysis/random_phenotype_3_lambda_per_ancestry_max_MAF_1e_3_', type,'.tsv')) 
  qq_path <- paste0(aou_qc_data_path, 'preliminary/random_phenotype_3_max_MAF_', max_MAF, '_', type, '_qq_data.tsv')
  if(!file.exists(qq_path)){
    qq_data <- raw_qq_data(max_MAF = max_MAF, type = type)
    write_tsv(qq_data, qq_path)
  }else{
    qq_data <- read_tsv(qq_path) %>%
      mutate(prevalence = if_else(is.na(prevalence), '0.0001', prevalence)) %>%
      mutate(prevalence = factor(prevalence, levels = rev(c('0.0001', '0.001', '0.01', '0.02', '0.05', '0.1', '0.2', '0.5', 'continuous'))))
  }
  
  max_data <- qq_data  %>%
    group_by(prevalence) %>%
    dplyr::summarize(max = max(observed, na.rm=T)) %>%
    dplyr::mutate(pheno_max = if_else(is.infinite(max), 300, max)) 
  lambda_data <- lambda_data %>%
    mutate(prevalence = if_else(prevalence == '1e-04', '0.0001', prevalence))%>%
    merge(., max_data, by = 'prevalence', all.x = T) %>%
    mutate(prevalence = factor(prevalence, levels = rev(c('0.0001', '0.001', '0.01', '0.02', '0.05', '0.1', '0.2', '0.5', 'continuous'))))
  return(list(lambda_data, qq_data))
}

make_qq_plot <- function(plot_data, type, fast_plot = TRUE){
  if(fast_plot){
    figure <- ggplot() + scattermore::geom_scattermore(
      aes(x =expected, y=observed, color = ancestry),
      pointsize = 15, pixels = c(1024, 1024), alpha = 0.5,
      data = plot_data[[2]]
    )+
      geom_abline(intercept = 0, slope = 1, lty=2) +
      # geom_hline(yintercept = -log10(0.05), lty=2) +
      labs(x="Expected -log10(p)", y="Observed -log10(p)", color = 'Genetic similarity group') +
      scale_color_manual(values = pop_colors, breaks = c('afr', 'amr', 'eas', 'eur', 'mid', 'sas','oth'), labels = pops) +
      scale_fill_manual(values = pop_colors, breaks = c('afr', 'amr', 'eas', 'mid', 'eur', 'sas','oth'), labels = pops)+
      themes +
      theme(plot.title = element_text(size = 13)) +
      geom_text(data = plot_data[[1]], aes(x = 0, y = -as.numeric(factor(ancestry))*pheno_max/15 + 16*pheno_max/15, label = paste0(as.character(round(lambda_gc, 3))), color = ancestry), hjust = -0.1, show.legend = F)
    
  }else{
    figure <- plot_data[[2]] %>%
      ggplot+ aes(y=observed,x=expected, color = ancestry) +
      geom_point(alpha = 0.5) +
      geom_abline(intercept = 0, slope = 1, lty=2) +
      geom_hline(yintercept = -log10(0.05), lty=2) +
      labs(x="Expected -log10(p)", y="Observed -log10(p)", color = 'Genetic similarity group') +
      scale_color_manual(values = pop_colors, breaks = c('afr', 'amr', 'eas', 'eur', 'mid', 'sas','oth'), labels =pops) +
      scale_fill_manual(values = pop_colors, breaks = c('afr', 'amr', 'eas', 'mid', 'eur', 'sas','oth'), labels = pops) +
      themes +
      theme(plot.title = element_text(size = 13)) +
      geom_text(data = plot_data[[1]], aes(x = 0, y = as.numeric(factor(ancestry))*(pheno_max-3.5)/4+3.5, label = paste0(as.character(round(lambda_gc, 3))), color = ancestry), hjust = -0.1, show.legend = F)
  }
  figure <- figure +
    theme(legend.position = 'top') +
    facet_wrap(~prevalence, scales = 'free', ncol = 3)  +
    guides(colour = guide_legend(nrow = 1, byrow = TRUE), label='none')
  output_figure(figure, 'supplement', paste0('supp_figure15_aou_v8_random_3_maxmaf_0.001_', type), height = 8, width = 10)
  return(figure)
}

plot_data <- prepare_qq_plot_data_v8(type='burden')
figure <- make_qq_plot(plot_data, 'burden', fast_plot = T)
# plot_data <- prepare_qq_plot_data_v8(type='skato')
# figure <- make_qq_plot(plot_data, 'skato', fast_plot = T)
# plot_data <- prepare_qq_plot_data_v8(type='skat')
# figure <- make_qq_plot(plot_data, 'skat', fast_plot = T)