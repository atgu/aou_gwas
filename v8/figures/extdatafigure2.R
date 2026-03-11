source('~/Dropbox (Partners HealthCare)/github_repo/ukbb_exomes_pleiotropy/R/constants.R')
source('~/Dropbox (Partners HealthCare)/github_repo/aou_gwas/R/constants.R')
library(magick)
manhattan_path <- '~/Dropbox (Partners HealthCare)/aou/axaou/v8/'


sub_pheno_info <- read_tsv('~/Dropbox (Partners HealthCare)/aou_v8/data/phenotype/aou_v8_phenotype_info_qced_with_meta.tsv') %>%
  filter(ancestry == 'META' & phenoname %in% c('EM_202.2', 'EM_236.1', 'A10BJ')) %>%
  select(phenoname, category, description, n_cases, n_controls, pheno_sex) %>%
  distinct() %>%
  mutate(prevalence = n_cases/(n_cases + n_controls)) %>%
  mutate(label = paste0(phenoname, ': ', str_to_sentence(description)),
         Ntotal = n_cases + n_controls,
         pheno_sex = factor(str_to_sentence(pheno_sex), levels = c('Both', 'Female', 'Male'))) %>%
  mutate(label = str_replace_all(label, 'glp', 'GLP')) 

# sub_pheno_info %>%
#   ggplot + aes(x = pheno_sex, y = prevalence, color = pheno_sex) + 
#   geom_point(position = position_jitterdodge()) +
#   labs(x = NULL) + scale_y_continuous(label = percent, breaks = c(0, 0.05, 0.1, 0.15 ,0.2, 0.25)) + 
#   facet_wrap(~label, strip.position = "bottom") +
#   themes + 
#   theme(
#     strip.placement = "inside", # Places labels outside the plot area
#     axis.ticks.x = element_blank(),
#     axis.text.x = element_blank(),
#     legend.direction = 'vertical',
#     legend.position = 'right'
#   ) + 
#   scale_color_manual(name = NULL, values = c('#757575', '#CC3333', '#0088DE'), breaks = c('Both', 'Female', 'Male'))

extfigure2_top <- sub_pheno_info %>%
  pivot_longer(., names_to = 'Case/Control', values_to = 'N', cols = c('n_cases', 'n_controls')) %>%
  mutate(`Case/Control` = factor(str_to_sentence(str_replace(`Case/Control`, 'n_', '')), levels = c('Controls', 'Cases'))) %>%
  ggplot + aes(x = pheno_sex, y = N, color = pheno_sex, fill = pheno_sex) + 
  geom_col(aes(alpha = `Case/Control`)) +
  labs(x = NULL) + scale_y_continuous(label = comma, limits = c(0, 350000)) + 
  facet_wrap(~label, strip.position = "bottom", labeller = label_wrap_gen(width = 25)) +
  themes + 
  theme(
    strip.placement = "inside", # Places labels outside the plot area
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    legend.direction = 'vertical',
    legend.position = 'right'
  ) + 
  scale_alpha_discrete(name = NULL, range = c(0.2, 1)) + 
  scale_fill_manual(name = NULL, values = c('#757575', '#CC3333', '#0088DE'), breaks = c('Both', 'Female', 'Male')) + 
  scale_color_manual(name = NULL, values = c('#757575', '#CC3333', '#0088DE'), breaks = c('Both', 'Female', 'Male')) + 
  geom_text(data = sub_pheno_info, aes(label = percent(prevalence), y = Ntotal), size = 3.5, vjust = -1, show.legend = F)
extfigure2_top


p_drug <- ggdraw() +
  draw_image(image_read(paste0(manhattan_path, 'ACAF/META/phenotype_A10BJ.META.Pvalue_log10.variant.manhattan.png')))+
  theme(plot.margin = unit(c(0.5, 0, 0, 0), "cm"))

p_drug_male <- ggdraw() +
  draw_image(image_read(paste0(manhattan_path, 'ACAF/META/phenotype_A10BJ_male.META.Pvalue_log10.variant.manhattan.png')))+
  theme(plot.margin = unit(c(0.5, 0, 0, 0), "cm"))

p_drug_female <- ggdraw() +
  draw_image(image_read(paste0(manhattan_path, 'ACAF/META/phenotype_A10BJ_female.META.Pvalue_log10.variant.manhattan.png')))+
  theme(plot.margin = unit(c(0.5, 0, 0, 0), "cm"))

library(ggVennDiagram)

load_venn_data <- function(sex){
  if(sex == 'male'){
    T2D_only   <- 8593
    Obesity_only <- 9605
    A10BJ_only <- 224
    
    T2D_Obesity   <- 7309   
    T2D_A10BJ     <- 944   
    Obesity_A10BJ <- 470  
    
    T2D_Obesity_A10BJ <- 2882   
  }else if(sex == 'female'){
    T2D_only   <- 7640
    Obesity_only <- 25877
    A10BJ_only <- 740
    
    T2D_Obesity   <- 11284   
    T2D_A10BJ     <- 945   
    Obesity_A10BJ <- 2326  
    
    T2D_Obesity_A10BJ <- 4968   
  }else{
    T2D_only   <- 16233
    Obesity_only <- 35482
    A10BJ_only <- 964
    
    T2D_Obesity   <- 18593   
    T2D_A10BJ     <- 1889  
    Obesity_A10BJ <- 2796  
    
    T2D_Obesity_A10BJ <- 7850   
    
  }
  # Triple-overlap items
  items_ABC <- paste0("T2D_Obesity_A10BJ_", seq_len(T2D_Obesity_A10BJ))
  
  # Pairwise-only regions
  items_AB  <- paste0("T2D_Obesity_", seq_len(T2D_Obesity))
  items_AC  <- paste0("T2D_A10BJ_",   seq_len(T2D_A10BJ))
  items_BC  <- paste0("Obesity_A10BJ_", seq_len(Obesity_A10BJ))
  
  # Singleton regions
  items_A <- paste0("T2D_",    seq_len(T2D_only))
  items_B <- paste0("Obesity_", seq_len(Obesity_only))
  items_C <- paste0("A10BJ_",   seq_len(A10BJ_only))
  
  setT2D    <- c(items_A, items_AB, items_AC, items_ABC)
  setObesity <- c(items_B, items_AB, items_BC, items_ABC)
  setA10BJ  <- c(items_C, items_AC, items_BC, items_ABC)
  
  venn_list <- list(
    T2D     = setT2D,
    Obesity = setObesity,
    A10BJ   = setA10BJ
  )
  
  return(venn_list)
}



myPalette <- colorRampPalette(rev(brewer.pal(11, "YlGnBu")))

plot_venn_diagram <- function(data){
  venn <- ggVennDiagram(data, set_color = c("#92c5de","#f4a582","#b2abd2"),
                             label = "both",          # "count", "percent", "both", "none"
                             label_alpha = 0.5,          # 0 = no white boxes behind numbers
                             label_font  = "Arial",# font family for region labels
                             label_color = "black",
                             label_size  = 2, hjust = 0,
                             set_size  = 3,    # size of set names
                             
                             edge_lty   = "solid",     # circle line type
                             edge_size  = 1) +         # circle line width) +
    # scale_fill_gradientn(colors = c('#887a7e', '#8f89bb', '#bea0b6')) +
    scale_fill_gradientn(
      colours = c('white')
    )+
    coord_sf(clip = "off") + 
    theme(legend.position = "none")
  venn$data$region$label <- scales::comma(venn$data$region$count)
  
  return(venn)
}




venn_male <- plot_venn_diagram(load_venn_data('male'))
venn_female <- plot_venn_diagram(load_venn_data('female'))
venn_both <- plot_venn_diagram(load_venn_data('both'))


extfigure2_both = ggpubr::ggarrange(p_drug, venn_both + theme(plot.margin = unit(c(0.5,0.5,0.2,0), 'cm')), 
                            labels = c('b',
                                       'c'
                                       
                            ), hjust=0, vjust = 1, widths = c(2,1.2),
                            font.label = list(size = 10, color = "black", face = "bold", family = 'Arial'), 
                            nrow=1)
extfigure2_both

output_figure(extfigure2_both, 'main', 'extfigure2_both', height =2, width = 8)

extfigure2_female = ggpubr::ggarrange(p_drug_female, venn_female + theme(plot.margin = unit(c(0.5,0.5,0.2,0), 'cm')), 
                                   labels = c('d',
                                              'e'
                                   ), hjust=0, vjust = 1, widths = c(2,1.2),
                                   font.label = list(size = 10, color = "black", face = "bold", family = 'Arial'), 
                                   nrow=1)
extfigure2_female

output_figure(extfigure2_female, 'main', 'extfigure2_female', height =2, width = 8)

extfigure2_male = ggpubr::ggarrange(p_drug_male, venn_male + theme(plot.margin = unit(c(0.5,0.5,0.2,0), 'cm')), 
                                      labels = c('f',
                                                 'g'
                                      ), hjust=0, vjust = 1, widths = c(2,1.2),
                                      font.label = list(size = 10, color = "black", face = "bold", family = 'Arial'), 
                                      nrow=1)
extfigure2_male

output_figure(extfigure2_male, 'main', 'extfigure2_male', height =2, width = 8)

extfigure2 = ggpubr::ggarrange(extfigure2_top+ theme(plot.margin = unit(c(0.6,0,0,0), 'cm')),
                                extfigure2_both + theme(plot.margin = unit(c(0.2,0,0,0), 'cm')), 
                               extfigure2_female + theme(plot.margin = unit(c(0.2,0,0,0), 'cm')), 
                               extfigure2_male + theme(plot.margin = unit(c(0.2,0,0,0), 'cm')), 
                               labels = c('a',
                                          '', '', ''
                               ), hjust=0, vjust = 1.2, heights = c(1.3,1.1,1.1, 1.1),
                            font.label = list(size = 10, color = "black", face = "bold", family = 'Arial'), 
                            ncol=1)
output_figure(extfigure2, 'main', 'extfigure2', height =7, width = 7)

